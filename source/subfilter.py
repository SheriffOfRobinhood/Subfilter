# -*- coding: utf-8 -*-
"""

  subfilter.py
    - This is the "subfilter module"
    - Defines many useful routines for the subfilter calculations.
    - examples of their use are present in subfilter_file.py

Created on Tue Oct 23 11:07:05 2018

@author: Peter Clark
"""
import os
import sys

import numpy as np
import xarray as xr
from dask.diagnostics import ProgressBar
import time

from scipy.signal import fftconvolve
import difference_ops as do

from thermodynamics_constants import *

test_level = 1

subfilter_version = '0.4.0'

subfilter_setup = {'write_sleeptime':3,
                   'use_concat':True,
                   'chunk_size':2**22 }

def save_field(dataset, field, write_to_file=True):
    fn = dataset['file'].split('/')[-1]
    if field.name not in dataset['ds']:
        print(f"Saving {field.name} to {fn}")
        dataset['ds'][field.name] = field
        if write_to_file:
            d = dataset['ds'][field.name].to_netcdf(
                    dataset['file'], mode='a', compute=False)
            with ProgressBar():
                results = d.compute()
            # This wait seems to be needed to give i/o time to flush caches.
            time.sleep(subfilter_setup['write_sleeptime'])
    else:
        print(f"{field.name} already in {fn}")
#    print(dataset['ds'])
    return dataset['ds'][field.name]

def re_chunk(f, chunks = None, xch = 'all', ych = 'all', zch = 'auto'):

    defn = 1

    if chunks is None:

        chunks={}
        sh = np.shape(f)
        for ip, dim in enumerate(f.dims):
            if 'x' in dim:
                if xch == 'all':
                    chunks[dim] = sh[ip]
                else:
                    chunks[dim] = np.min([xch, sh[ip]])
            elif 'y' in dim:
                if ych == 'all':
                    chunks[dim] = sh[ip]
                else:
                    chunks[dim] = np.min([ych, sh[ip]])
            elif 'z' in dim:
                if zch == 'all':
                    chunks[dim] = sh[ip]
                elif zch == 'auto':
                    chunks[dim] = 'auto'
                else:
                    chunks[dim] = np.min([zch, sh[ip]])
            else:
                chunks[f.dims[ip]] = defn

    f = f.chunk(chunks=chunks)

    return f


def filter_variable_list(source_dataset, ref_dataset, derived_dataset,
                         filtered_dataset, options, filter_def,
                         var_list=None, grid='p') :
    """
    Create filtered versions of input variables.

    Filtered variables created on required grid,stored in derived_dataset.

    Args:
        source_dataset  : NetCDF dataset for input
        ref_dataset     : NetCDF dataset for input containing reference
                          profiles. Can be None
        derived_dataset : NetCDF dataset for derived data
        filtered_dataset: NetCDF dataset for derived data
        options         : General options e.g. FFT method used.
        filter_def      : 1 or 2D filter
        var_list=None   : List of variable names.
        default provided by get_default_variable_list()
        grid='p'        : Grid - 'u','v','w' or 'p'

    Returns
    -------
        list : list of strings representing variable names.

    @author: Peter Clark

    """
    if (var_list==None):
        var_list = get_default_variable_list()
        print("Default list:\n",var_list)

    for vin in var_list:

        op_var  = get_data_on_grid(source_dataset, ref_dataset,
                                   derived_dataset, vin, options,
                                   grid)

        v = op_var.name

        if v+"_r" not in filtered_dataset['ds'].variables \
            or v+"_s" not in filtered_dataset['ds'].variables:

            ncvar_r, ncvar_s = filter_field(op_var,
                                            filtered_dataset,
                                            options, filter_def,
                                            grid=grid)

    return var_list

def filter_variable_pair_list(source_dataset, ref_dataset, derived_dataset,
                              filtered_dataset, options, filter_def,
                              var_list=None, grid='p') :
    """
    Create filtered versions of pairs input variables.

    Filtered variable pairs created on required grid,
    stored in derived_dataset.

    Args:
        source_dataset  : NetCDF dataset for input
        ref_dataset     : NetCDF dataset for input containing reference
                          profiles. Can be None
        derived_dataset : NetCDF dataset for derived data
        filtered_dataset: NetCDF dataset for derived data
        options         : General options e.g. FFT method used.
        filter_def      : 1 or 2D filter
        var_list=None   : List of variable names.
        default provided by get_default_variable_pair_list()

    Returns
    -------
        list : list of lists of pairs strings representing variable names.

    @author: Peter Clark

    """
    if (var_list==None):
        var_list = get_default_variable_pair_list()
        print("Default list:\n",var_list)

    for v in var_list:

        print(f"Calculating s({v[0]:s},{v[1]:s})")
        svars = quadratic_subfilter(source_dataset, ref_dataset,
                                  derived_dataset, filtered_dataset, options,
                                  filter_def, v[0], v[1], grid=grid)

        (s_var1var2, var1var2, var1var2_r, var1var2_s) = svars

        if options['save_all'].lower() == 'yes':
            save_field(derived_dataset, var1var2)

        for f in (s_var1var2, var1var2_r, var1var2_s):
            save_field(filtered_dataset, f)


    return var_list


# Flags are: 'u-grid, v-grid, w-grid'

def get_default_variable_list() :
    """
    Provide default variable list.

       Returns:
           var_list.

    The default is::

     var_list = [
            "u",
            "v",
            "w",
            "th",
            "th_v",
            "th_L",
            "q_vapour",
            "q_cloud_liquid_mass",
            "q_total"]

    @author: Peter Clark

    """

    if test_level == 1:
# For testing
        var_list = [
            "u",
            "w",
            "th",
            ]
    elif test_level == 2:
# For testing
        var_list = ["u","w","th", "th_v", "th_L", "q_total"]
        var_list = [
            "u",
            "w",
            "th_L",
            "q_total",
            ]
    else:
        var_list = [
            "u",
            "v",
            "w",
            "th",
            "th_v",
            "th_L",
            "q_vapour",
            "q_cloud_liquid_mass",
            "q_total"]
    return var_list

def get_default_variable_pair_list() :
    """
    Provide default variable pair list.
       Returns:
        list : list of lists of pairs strings representing variable names.

    The default is::

        var_list = [
                ["u","u"],
                ["u","v"],
                ["u","w"],
                ["v","v"],
                ["v","w"],
                ["w","w"],
                ["u","th"],
                ["v","th"],
                ["w","th"],
                ["u","th_v"],
                ["v","th_v"],
                ["w","th_v"],
                ["u","th_L"],
                ["v","th_L"],
                ["w","th_L"],
                ["u","q_vapour"],
                ["v","q_vapour"],
                ["w","q_vapour"],
                ["u","q_cloud_liquid_mass"],
                ["v","q_cloud_liquid_mass"],
                ["w","q_cloud_liquid_mass"],
                ["u","q_total"],
                ["v","q_total"],
                ["w","q_total"],
                ["th_L","th_L"],
                ["th_L","q_total"],
                ["q_total","q_total"],
                ["th_L","q_vapour"],
                ["th_L","q_cloud_liquid_mass"],
              ]

    @author: Peter Clark


    """
    if test_level == 1:
# For testing
        var_list = [
                ["w","th"],
              ]
    elif test_level == 2:
# For testing
        var_list = [
                ["u","w"],
                ["w","w"],
                ["u","th"],
                ["w","th"],
                ["w","th_L"],
                ["w","q_total"],
              ]
    else:
        var_list = [
                ["u","u"],
                ["u","v"],
                ["u","w"],
                ["v","v"],
                ["v","w"],
                ["w","w"],
                ["u","th"],
                ["v","th"],
                ["w","th"],
                ["u","th_v"],
                ["v","th_v"],
                ["w","th_v"],
                ["u","th_L"],
                ["v","th_L"],
                ["w","th_L"],
                ["u","q_vapour"],
                ["v","q_vapour"],
                ["w","q_vapour"],
                ["u","q_cloud_liquid_mass"],
                ["v","q_cloud_liquid_mass"],
                ["w","q_cloud_liquid_mass"],
                ["u","q_total"],
                ["v","q_total"],
                ["w","q_total"],
                ["th_L","th_L"],
                ["th_L","q_total"],
                ["q_total","q_total"],
                ["th_L","q_vapour"],
                ["th_L","q_cloud_liquid_mass"],
              ]
    return var_list

def convolve(field, options, filter_def, dims):
    """
    Convolve field filter using fftconvolve using padding.

    Args:
        field      : field array
        options    : General options e.g. FFT method used.
        filter_def : 1 or 2D filter array

    Returns:
        ndarray : field convolved with filter_def

    @author: Peter Clark

    """

    if len(np.shape(field)) > len(np.shape(filter_def)):
        edims = tuple(np.setdiff1d(np.arange(len(np.shape(field))), dims))
        filter_def = np.expand_dims(filter_def, axis=edims)

    if options['FFT_type'].upper() == 'FFTCONVOLVE':

        pad_len = np.max(np.shape(filter_def))//2

        pad_list = []
        for i in range(len(np.shape(field))):
            if i in dims:
                pad_list.append((pad_len,pad_len))
            else:
                pad_list.append((0,0))

        field = np.pad(field, pad_list, mode='wrap')
        result = fftconvolve(field, filter_def, mode='same', axes=dims)

        padspec = []
        for d in range(len(np.shape(field))):
            if d in dims:
                padspec.append(slice(pad_len,-pad_len))
            else:
                padspec.append(slice(0,None))
        padspec = tuple(padspec)

        result = result[padspec]

    elif options['FFT_type'].upper() == 'FFT':

        if len(np.shape(filter_def)) == 1:
            fft_field = np.fft.fft(field, axes=dims)
            fft_filtered_field = fft_field * filter_def
            result = np.fft.ifft(fft_filtered_field, axes=dims)
        else:
            fft_field = np.fft.fft2(field, axes=dims)
            fft_filtered_field = fft_field * filter_def
            result = np.fft.ifft2(fft_filtered_field, axes=dims)
        result = result.real

    elif options['FFT_type'].upper() == 'RFFT':

        if len(np.shape(filter_def)) == 1:
            fft_field = np.fft.rfft(field, axes=dims)
            fft_filtered_field = fft_field * filter_def
            result = np.fft.irfft(fft_filtered_field, axes=dims)
        else:
            fft_field = np.fft.rfft2(field, axes=dims)
            fft_filtered_field = fft_field * filter_def
            result = np.fft.irfft2(fft_filtered_field, axes=dims)
        result = result.real

    return result

def pad_to_len2D(field, newlen, mode='constant'):
    sf = np.shape(field)
    padlen = newlen - sf[0]
    padleft = padlen - padlen//2
    padright = padlen - padleft
    padfield = np.pad(field, ((padleft,padright),), mode=mode)
    return padfield


def filtered_field_calc(var, options, filter_def):
    """
    Split field into resolved (field_r) and subfilter (field_s).
    Note: this routine has a deliberate side effect, to store the fft or rfft
    of the filter in filter_def for subsequent re-use.

    Args:
        var             : dict cantaining variable info
        options         : General options e.g. FFT method used.
        filter_def      : 1 or 2D filter

    Returns:
        dicts cantaining variable info : [var_r, var_s]

    @author: Peter Clark

    """

    vname = var.name
    field = var.data
    vdims = var.dims

#    print(field)
    sh = np.shape(field)

    if filter_def.attributes['ndim'] == 1:

        axis = find_var(vdims, ['x'])

    elif filter_def.attributes['ndim'] == 2:

        axis = find_var(vdims, ['x', 'y'])


    if filter_def.attributes['filter_type'] == 'domain' :

        ax = list(axis)

        si = np.asarray(field.shape)
        si[ax] = 1

        field_r = np.mean(field[...], axis=axis)
        field_s = field[...] - np.reshape(field_r, si)

        rdims =  []
        rcoords = {}
        for i, d in enumerate(vdims):
            if i not in axis:
                rdims.append(d)
                rcoords[d] = var.coords[d]
        rdims = tuple(rdims)


    else :

        print(f"Filtering using {options['FFT_type']}")

        if options['FFT_type'].upper() == 'FFTCONVOLVE':

            field_r = convolve(field, options, filter_def.data, axis)

        elif options['FFT_type'].upper() == 'FFT':

            if filter_def.attributes['ndim'] == 1:

                if 'fft' not in filter_def.__dict__:
                    sf = np.shape(filter_def.data)
                    if sh[axis[0]] != sf[0]:
                        padfilt = pad_to_len2D(filter_def.data, sh[axis[0]])
                    else:
                        padfilt = filter_def.data.copy()
                    # This shift of the filter is necessary to get the phase
                    # information right.
                    padfilt = np.fft.ifftshift(padfilt)#????? Test
                    filter_def.fft = np.fft.fft(padfilt)

            else:

                if 'fft' not in filter_def.__dict__:
                    sf = np.shape(filter_def.data)
                    if sh[axis[0]] != sf[0] or sh[axis[1]] != sf[1]:
                        padfilt = pad_to_len2D(filter_def.data, sh[axis[0]])
                    else:
                        padfilt = filter_def.data.copy()
                    # This shift of the filter is necessary to get the phase
                    # information right.
                    padfilt = np.fft.ifftshift(padfilt)
                    filter_def.fft = np.fft.fft2(padfilt)

            field_r = convolve(field, options, filter_def.fft, axis)

        elif options['FFT_type'].upper() == 'RFFT':

            if filter_def.attributes['ndim'] == 1:

                if 'rfft' not in filter_def.__dict__:
                    sf = np.shape(filter_def.data)
                    if sh[axis[0]] != sf[0]:
                        padfilt = pad_to_len2D(filter_def.data, sh[axis[0]])
                    else:
                        padfilt = filter_def.data.copy()
                    # This shift of the filter is necessary to get the phase
                    # information right.
                    padfilt = np.fft.ifftshift(padfilt)#????? Test
                    filter_def.fft = np.fft.fft(padfilt)

            else:

                if 'rfft' not in filter_def.__dict__:
                    sf = np.shape(filter_def.data)
                    if sh[axis[0]] != sf[0] or sh[axis[1]] != sf[1]:
                        padfilt = pad_to_len2D(filter_def.data, sh[axis[0]])
                    else:
                        padfilt = filter_def.data.copy()
                    # This shift of the filter is necessary to get the phase
                    # information right.
                    padfilt = np.fft.ifftshift(padfilt)
                    filter_def.rfft = np.fft.rfft2(padfilt)

            field_r = convolve(field, options, filter_def.rfft, axis)
            rdims = var.dims
            rcoords = var.coords

        field_s = field[...] - field_r

    sdims = var.dims
    scoords = var.coords

    var_r = xr.DataArray(field_r, name = vname+'_r', dims=rdims,
                      coords=rcoords)
    var_s = xr.DataArray(field_s, name = vname+'_s', dims=sdims,
                      coords=scoords)

    return (var_r, var_s)

def setup_data_file(source_file, derived_dataset_name,
                    override=False) :
    """
    Create NetCDF dataset for derived data in destdir.

    File name is original file name concatenated with filter_def.id.

    Args:
        source_file     : Input NetCDF file name.
        derived_dataset_name : Created NetCDF file name.
        override=False  : if True force creation of file

    Returns:
        derived_dataset
        exists

    @author: Peter Clark
    """
    exists = os.path.isfile(derived_dataset_name)

    if exists and not override :

        derived_dataset = xr.open_dataset(derived_dataset_name)

    else :

        exists = False

        derived_dataset = xr.Dataset()
        derived_dataset = derived_dataset.assign_attrs(
            {'Parent file':source_file})

        derived_dataset.to_netcdf(derived_dataset_name, mode='w')

    print(derived_dataset)
    return derived_dataset, exists

def setup_derived_data_file(source_file, destdir, fname,
                            options, override=False) :
    """
    Create NetCDF dataset for derived data in destdir.

    File name is original file name concatenated with filter_def.id.


    Args:
        source_file     : NetCDF file name.
        destdir         : Directory for derived data.
        override=False  : if True force creation of file

    Returns:
        derived_dataset_name, derived_dataset

    @author: Peter Clark
    """
    derived_dataset_name = os.path.basename(source_file)
    derived_dataset_name = ('.').join(derived_dataset_name.split('.')[:-1])
    derived_dataset_name = destdir+derived_dataset_name + "_" + fname + ".nc"

    derived_dataset, exists = setup_data_file(source_file,
                    derived_dataset_name, override=override)

    dataset = {'file': derived_dataset_name, 'ds':derived_dataset}

    return dataset, exists

def setup_filtered_data_file(source_file, destdir, fname,
                            options, filter_def, override=False) :
    """
    Create NetCDF dataset for filtered data in destdir.

    File name is original file name concatenated with filter_def.id.

    Args:
        source_file     : NetCDF file name.
        destdir         : Directory for derived data.
        options         : General options e.g. FFT method used.
        filter_def      : Filter
        options         : General options e.g. FFT method used.
        override=False  : if True force creation of file

    Returns:
        filtered_dataset_name, filtered_dataset

    @author: Peter Clark
    """
    filtered_dataset_name = os.path.basename(source_file)
    filtered_dataset_name = ('.').join(filtered_dataset_name.split('.')[:-1])
    filtered_dataset_name = destdir+filtered_dataset_name + "_" + fname \
        + "_" + filter_def.id + ".nc"
    filtered_dataset, exists = setup_data_file(source_file,
                    filtered_dataset_name, override=override)

    filtered_dataset = filtered_dataset.assign_attrs(
        {'filter_def_id' : filter_def.id})
    filtered_dataset = filtered_dataset.assign_attrs(filter_def.attributes)
    filtered_dataset = filtered_dataset.assign_attrs(options)

    filtered_dataset.to_netcdf(filtered_dataset_name, mode='a')

    dataset = {'file':filtered_dataset_name, 'ds':filtered_dataset}

    return dataset, exists

def get_data(source_dataset, ref_dataset, var_name, options) :
    """
    Extract data from source NetCDF dataset or derived data.

    Currently supported derived data are::

        'th_L'     : Liquid water potential temperature.
        'th_v'     : Virtual potential temperature.
        'q_total'  : Total water.
        'buoyancy' : (g/mean_th_v)*(th_v-mean_th_v), where the mean is the domain mean.

    Returns:
        variable, variable_dimensions, variable_grid_properties

    @author: Peter Clark

    """
#   Mapping of data locations on grid via logical triplet:
#   logical[u-point,v-point,w-point]
#          [False,  False,  False  ] --> (p,th,q)-point
    var_properties = {"u":{'grid':[True,False,False], "units":'m s-1'},
                      "v":{'grid':[False,True,False], "units":'m s-1'},
                      "w":{'grid':[False,False,True], "units":'m s-1'},
                      "th":{'grid':[False,False,False], "units":'K'},
                      "p":{'grid':[False,False,False], "units":'Pa'},
                      "q_vapour":{'grid':[False,False,False], "units":'kg/kg'},
                      "q_cloud_liquid_mass":{'grid':[False,False,False],
                                             "units":'kg/kg'},
                      }

    od = options_database(source_dataset)
    if od is None:
        dx = options['dx']
        dy = options['dy']
    else:
        dx = float(od['dxx'])
        dy = float(od['dyy'])


    print(f'Retrieving {var_name:s}.')
    try :
        vard = source_dataset[var_name]

        # Change 'timeseries...' variable to 'time'

        [itime] = find_var(vard.dims, ['time'])
        if itime is not None:
            vard = vard.rename({vard.dims[itime]: 'time'})

        # Add correct x and y grids.

        if var_name in var_properties:

            vp = var_properties[var_name]['grid']

            if 'x' in vard.dims:
                nx = vard.shape[vard.get_axis_num('x')]

                if vp[0] :
                    x = (np.arange(nx) + 0.5) * np.float64(od['dxx'])
                    xn = 'x_u'
                else:
                    x = np.arange(nx) * np.float64(od['dxx'])
                    xn = 'x_p'

                vard = vard.rename({'x':xn})
                vard.coords[xn] = x

            if 'y' in vard.dims:
                ny = vard.shape[vard.get_axis_num('y')]
                if vp[1] :
                    y = (np.arange(ny) + 0.5) * np.float64(od['dyy'])
                    yn = 'y_v'
                else:
                    y = np.arange(ny) * np.float64(od['dyy'])
                    yn = 'y_p'

                vard = vard.rename({'y':yn})
                vard.coords[yn] = y

            if 'z' in vard.dims and not vp[2]:
                zn = source_dataset.coords['zn']
                vard = vard.rename({'z':'zn'})
                vard.coords['zn'] = zn.data

            if 'zn' in vard.dims and vp[2]:
                z = source_dataset.coords['z']
                vard = vard.rename({'zn':'z'})
                vard.coords['z'] = z.data

            vard.attrs['units'] = var_properties[var_name]['units']

#        print(vard)
        if var_name == 'th' :
            thref = get_thref(ref_dataset, options)
#            thref = np.expand_dims(thref, axis=dims)
            vard += thref

    except :

        if var_name == 'th_L' :
#            rhoref = ref_dataset.variables['rhon'][-1,...]
            theta = get_data(source_dataset, ref_dataset, 'th',
                                           options)
            (pref, piref) = get_pref(source_dataset, ref_dataset,  options)
            q_cl = get_data(source_dataset, ref_dataset,
                                    'q_cloud_liquid_mass', options)
            vard = theta - L_over_cp * q_cl / piref


        elif var_name == 'th_v' :
#            rhoref = ref_dataset.variables['rhon'][-1,...]
            theta = get_data(source_dataset, ref_dataset, 'th',
                                           options)
            thref = get_thref(ref_dataset, options)
            q_v = get_data(source_dataset, ref_dataset,
                                         'q_vapour', options)
            q_cl = get_data(source_dataset, ref_dataset,
                                    'q_cloud_liquid_mass', options)
            vard = theta + thref * (c_virtual * q_v - q_cl)

        elif var_name == 'q_total' :
#            rhoref = ref_dataset.variables['rhon'][-1,...]
            q_v = get_data(source_dataset, ref_dataset,
                                         'q_vapour', options)
            q_cl = get_data(source_dataset, ref_dataset,
                                    'q_cloud_liquid_mass', options)
            vard = q_v + q_cl

        elif var_name == 'buoyancy':
            th_v = get_data(source_dataset, ref_dataset, 'th_v',
                                           options)
            # get mean over horizontal axes
            # haxes = find_var(vardim, ['x','y'])
            mean_thv = th_v.mean(dim=('x','y'))
            vard = grav * (th_v - mean_thv)/mean_thv

        else :

            sys.exit(f"Data {var_name:s} not in dataset.")
    print(vard)

    return vard

def get_and_transform(source_dataset, ref_dataset, var_name, options,
                      grid='p'):

    var = get_data(source_dataset, ref_dataset, var_name, options)

#    print(var)

    vp = ['x_u' in var.dims,
          'y_v' in var.dims,
          'z' in var.dims]

    if grid=='p' :
        if vp[0] :
            print("Mapping {} from u grid to p grid.".format(var_name))
            var = do.field_on_u_to_p(var)
        if vp[1] :
            print("Mapping {} from v grid to p grid.".format(var_name))
            var = do.field_on_v_to_p(var)
        if vp[2] :
            print("Mapping {} from w grid to p grid.".format(var_name))
            z = source_dataset["z"]
            zn = source_dataset["zn"]
            var = do.field_on_w_to_p(var, zn)

    elif grid=='u' :
        if not ( vp[0] or vp[1] or vp[2]):
            print("Mapping {} from p grid to u grid.".format(var_name))
            var = do.field_on_p_to_u(var)
        if vp[1] :
            print("Mapping {} from v grid to u grid.".format(var_name))
            var = do.field_on_v_to_p(var)
            var = do.field_on_p_to_u(var)
        if vp[2] :
            print("Mapping {} from w grid to u grid.".format(var_name))
            z = source_dataset["z"]
            zn = source_dataset["zn"]
            var = do.field_on_w_to_p(var, zn)
            var = do.field_on_p_to_u(var)

    elif grid=='v' :
        if not ( vp[0] or vp[1] or vp[2]):
            print("Mapping {} from p grid to v grid.".format(var_name))
            var = do.field_on_p_to_v(var)
        if vp[0] :
            print("Mapping {} from u grid to v grid.".format(var_name))
            var = do.field_on_u_to_p(var)
            var = do.field_on_v_to_p(var)
        if vp[2] :
            print("Mapping {} from w grid to v grid.".format(var_name))
            z = source_dataset["z"]
            zn = source_dataset["zn"]
            var = do.field_on_w_to_p(var, zn)
            var = do.field_on_p_to_v(var)

    elif grid=='w' :
        z = source_dataset["z"]
        zn = source_dataset["zn"]
        if not ( vp[0] or vp[1] or vp[2]):
            print("Mapping {} from p grid to w grid.".format(var_name))
            var = do.field_on_p_to_w(var, z)
        if vp[0] :
            print("Mapping {} from u grid to w grid.".format(var_name))
            var = do.field_on_u_to_p(var)
            var = do.field_on_p_to_w(var, z)
        if vp[1] :
            print("Mapping {} from v grid to w grid.".format(var_name))
            var = do.field_on_v_to_p(var)
            var = do.field_on_p_to_w(var, z)

    else:
        print("Illegal grid ",grid)
    # print(var)

#    print(zvar)
    var = re_chunk(var)
#    print(var)

    return var

def get_data_on_grid(source_dataset, ref_dataset, derived_dataset, var_name,
                     options, grid='p') :
    """
    Read in 3D data from NetCDF file and, where necessary, interpolate to p grid.

    Assumes first dimension is time.

    Args:
        source_dataset  : NetCDF dataset
        ref_dataset     : NetCDF dataset containing reference profiles.
        derived_dataset : NetCDF dataset for derived data
        var_name        : Name of variable
        options         : General options e.g. FFT method used.
		grid='p'        : Destination grid. 'u', 'v', 'w' or 'p'.

    Returns:
        variable_dimensions, variable_grid_properties.

    @author: Peter Clark
    """
    grid_properties = {"u":[True,False,False],
                       "v":[False,True,False],
                       "w":[False,False,True],
                       "p":[False,False,False],
                      }

    ongrid = '_on_'+grid
    vp = grid_properties[grid]

    var_found = False
    # Logic here:
    # If var_name already qualified with '_on_x', where x is a grid
    # then if x matches required output grid, see in derived_dataset
    # already, and use if it is.
    # otherwise strip '_on_x' and go back to source data as per default.

    # First, find op_name
    # Default
    op_var_name = var_name + ongrid

    if len(var_name) > 5:
        if var_name[-5:] == ongrid:
            op_var_name = var_name
        elif var_name[-5:-1] == '_on_':
            var_name = var_name[:-5]
            op_var_name = var_name[:-5] + ongrid

    op_var = { 'name' : op_var_name }

    if options['save_all'].lower() == 'yes':

        if op_var_name in derived_dataset['ds'].variables:

            op_var = derived_dataset['ds'][op_var_name]
            print(f'Retrieved {op_var_name:s} from derived dataset.')
            var_found = True


    if not var_found:
        op_var = get_and_transform(source_dataset, ref_dataset,
                                   var_name, options, grid=grid)
        op_var.name = op_var_name

        if options['save_all'].lower() == 'yes':
            op_var = save_field(derived_dataset, op_var)
            # print(op_var)

    return op_var

def deformation(source_dataset, ref_dataset, derived_dataset,
                options, grid='w') :
    """
    Read in 3D data from NetCDF file and, where necessary, interpolate to p grid.

    Assumes first dimension is time.

    Args:
        source_dataset  : NetCDF dataset
        ref_dataset     : NetCDF dataset for input containing reference
                          profiles. Can be None
        derived_dataset : NetCDF dataset for derived data.
        options         : General options e.g. FFT method used.
        var_name        : Name of variable

    Returns:
        data array.

    @author: Peter Clark
    """

    if 'deformation' in derived_dataset['ds']:
        deformation = derived_dataset['ds']['deformation']
        return deformation

    u = get_data(source_dataset, ref_dataset, "u", options)
    [iix, iiy, iiz] = find_var(u.dims, ['x', 'y', 'z'])

    sh = np.shape(u)

    max_ch = subfilter_setup['chunk_size']

    nch = int(sh[iix]/(2**int(np.log(sh[iix]*sh[iiy]*sh[iiz]/max_ch)/np.log(2)/2)))

    print(f'nch={nch}')
#    nch = 32

    u = re_chunk(u, xch=nch, ych=nch, zch = 'all')
    # print(u)
    v = get_data(source_dataset, ref_dataset, "v", options)
    v = re_chunk(v, xch=nch, ych=nch, zch = 'all')
    # print(v)
    w = get_data(source_dataset, ref_dataset, "w", options)
    w = re_chunk(w, xch=nch, ych=nch, zch = 'all')

    z = source_dataset["z"]
    zn = source_dataset["zn"]

    ux = do.d_by_dx_field_on_u(u, z, grid = grid )
    # print(ux)
    # ux = save_field(derived_dataset, ux)
    print(ux)

    uy = do.d_by_dy_field_on_u(u, z, grid = grid )
    # print(uy)
    # uy = save_field(derived_dataset, uy)
    print(uy)

    uz = do.d_by_dz_field_on_u(u, z, grid = grid )
    # print(uz)
    # uz = save_field(derived_dataset, uz)
    print(uz)

    vx = do.d_by_dx_field_on_v(v, z, grid = grid )
    # print(vx)
    # vx = save_field(derived_dataset, vx)
    print(vx)

    vy = do.d_by_dy_field_on_v(v, z, grid = grid )
    # print(vy)
    # vy = save_field(derived_dataset, vy)
    print(vy)

    vz = do.d_by_dz_field_on_v(v, z, grid = grid )
    # print(vz)
    # vz = save_field(derived_dataset, vz)
    print(vz)

    wx = do.d_by_dx_field_on_w(w, zn, grid = grid )
    # print(wx)
    # wx = save_field(derived_dataset, wx)
    print(wx)

    wy = do.d_by_dy_field_on_w(w, zn, grid = grid )
    # print(wy)
    # wy = save_field(derived_dataset, wx)
    print(wy)

    wz = do.d_by_dz_field_on_w(w, zn, grid = grid )
    # print(wz)
    # wz = save_field(derived_dataset, wx)
    print(wz)

    if subfilter_setup['use_concat']:

        print('Concatenating derivatives')

        t0 = xr.concat([ux, uy, uz], dim='j', coords='minimal', compat='override')
        print(t0)
        t1 = xr.concat([vx, vy, vz], dim='j', coords='minimal', compat='override')
        print(t1)
        t2 = xr.concat([wx, wy, wz], dim='j', coords='minimal', compat='override')
        print(t2)

        defm = xr.concat([t0, t1, t2], dim='i')

        defm.name = 'deformation'
        defm.attrs={'units':'s-1'}

        print(defm)

#        defm = re_chunk(defm, zch = 1)

        if options['save_all'].lower() == 'yes':
            defm = save_field(derived_dataset, defm)

    else:

        t = [[ux, uy, uz],
             [vx, vy, vz],
             [wx, wy, wz],]
        defm = {}
        for i, t_i in enumerate(t):
            for j, u_ij in enumerate(t_i):
        #        u_ij = u_ij.expand_dims({'i':[i], 'j':[j]})
                u_ij.name = f'deformation_{i:1d}{j:1d}'
                if options['save_all'].lower() == 'yes':
                    u_ij = save_field(derived_dataset, u_ij)
                defm[f'{i:1d}_{j:1d}']=u_ij


    # print(derived_dataset)

    return defm


def filter_field(var, filtered_dataset, options, filter_def,
                 grid='p') :
    """
    Create filtered versions of input variable on required grid, stored in filtered_dataset.

    Args:
        var            : dict cantaining variable info
        filtered_dataset : NetCDF dataset for derived data.
        options         : General options e.g. FFT method used.
        filter_def      : 1 or 2D filter.
        default provided by get_default_variable_list()
        grid='p'        : Grid - 'u','v','w' or 'p'

    Returns:
        ncvar_r, ncvar_s: Resolved and subfilter fields as netcf variables in
                          filtered_dataset.

    @author: Peter Clark

    """
    vname = var.name
    vname_r = vname+'_r'
    vname_s = vname+'_s'

    if vname_r in filtered_dataset['ds'] and vname_r in filtered_dataset['ds']:

        print("Reading ", vname_r, vname_s)
        var_r = filtered_dataset['ds'][vname_r]
        var_s = filtered_dataset['ds'][vname_s]

    else:

        print(f"Filtering {vname:s}")

        # Calculate resolved and unresolved parts of var

        (var_r, var_s) = filtered_field_calc(var, options, filter_def)

        var_r = save_field(filtered_dataset, var_r)
        var_s = save_field(filtered_dataset, var_s)

    return (var_r, var_s)

def filtered_deformation(source_dataset, ref_dataset, derived_dataset,
                         filtered_dataset,
                         options, filter_def,
                         grid='p'):
    """
    Create filtered versions of deformation field.

    Args:
        source_dataset  : NetCDF input dataset
        ref_dataset     : NetCDF dataset for input containing reference
                          profiles. Can be None
        derived_dataset : NetCDF dataset for derived data.
        filtered_dataset: NetCDF dataset for derived data
        options         : General options e.g. FFT method used.
        filter_def      : 1 or 2D filter.
        grid='p'        : Grid - 'u','v','w' or 'p'

    Returns:
        ncvar_r, ncvar_s: Resolved and subfilter fields as netcf variables
        in filtered_dataset.

    @author: Peter Clark

    """

#:math:`\frac{\partial u_i}{\partial{x_j}`

    d_var = deformation(source_dataset, ref_dataset, derived_dataset,
                        options, grid=grid)

    d_var = re_chunk(d_var)

    (d_var_r, d_var_s) = filter_field(d_var, filtered_dataset,
                                      options, filter_def, grid=grid)


    return (d_var_r, d_var_s)

def shear(d, no_trace=True) :
    trace = 0
    vname = ''
    if no_trace :
        for k in range(3) :
            trace = trace + d.isel(i=k, j=k)
        trace = (2.0/3.0) * trace
        vname = 'n'

    mod_S = 0

    S = []
    i_j = []
    for k in range(3) :
        for l in range(k,3) :

            S_kl = d.isel(i=k, j=l) + d.isel(i=l, j=k)

            if k == l :
                S_kl = S_kl - trace
                mod_S += 0.5 * S_kl * S_kl
            else :
                mod_S +=  S_kl * S_kl

            S.append(S_kl)
            i_j.append(f'{k:d}_{l:d}')

    S = xr.concat(S, dim='i_j', coords='minimal', compat='override')
    S.coords['i_j'] = i_j
    S.name = 'shear' + vname
    S.attrs={'units':'s-1'}
    S = re_chunk(S)

    mod_S.name = 'mod_S' + vname
    mod_S.attrs={'units':'s-2'}
    S = re_chunk(S)

    return S, mod_S

def vorticity(d):

    v_i = []
    for i in range(3) :
        j=(i+1)%3
        k=(i+2)%3
        v_i.append(d.isel(i=k, j=j) - d.isel(i=j, j=k))

    v = xr.concat(v_i, dim='i')
    v.name='vorticity'
    v.attrs={'units':'s-1'}
    v = re_chunk(v)

    return v

def quadratic_subfilter(source_dataset,  ref_dataset, derived_dataset,
                        filtered_dataset, options, filter_def,
                        v1_name, v2_name, grid='p') :
    """
    Create filtered versions of pair of input variables on required grid, stored in derived_dataset.
    Computes :math:`s(\phi,\psi) = (\phi\psi)^r - \phi^r\psi^r.`

    Args:
        source_dataset  : NetCDF dataset for input
        ref_dataset     : NetCDF dataset for input containing reference
                          profiles. Can be None
        derived_dataset : NetCDF dataset for derived data
        filtered_dataset: NetCDF dataset for derived data
        options         : General options e.g. FFT method used.
        filter_def      : 1 or 2D filter
        v1_name         : Variable names.
        v2_name         : Variable names.

    Returns:
        s(var1,var2) data array.
        vdims dimensions of var1


    @author: Peter Clark

    """

    v1 = get_data_on_grid(source_dataset,  ref_dataset,
                          derived_dataset, v1_name, options,
                          grid=grid)

    (var1_r, var1_s) = filter_field(v1, filtered_dataset, options,
                                    filter_def, grid=grid)


    v2 = get_data_on_grid(source_dataset, ref_dataset,
                          derived_dataset, v2_name, options,
                          grid=grid)

    (var2_r, var2_s) = filter_field(v2, filtered_dataset, options,
                                    filter_def, grid=grid)


    var1var2 = v1 * v2
    var1var2.name = v1.name + '.' + v2.name

    print(f"Filtering {v1_name:s}*{v2_name:s}")
    (var1var2_r, var1var2_s) = filtered_field_calc(var1var2, options,
                                                 filter_def )

    s_var1var2 = var1var2_r - var1_r * var2_r

    s_var1var2.name = f"s({v1_name:s},{v2_name:s})_on_{grid:s}"

    return (s_var1var2, var1var2, var1var2_r, var1var2_s)

def bytarr_to_dict(d):

    # Converted for xarray use
    res = {}
    for i in range(np.shape(d)[0]):
        opt = d[i,0].decode('utf-8')
        val = d[i,1].decode('utf-8')

        res[opt] = val
    return res

def options_database(source_dataset):
    '''
    Convert options_database in source_dataset to dictionary.

    Parameters
    ----------
    source_dataset : netCDF4 file
        MONC output file.

    Returns
    -------
    options_database : dict

    '''

    # Converted to xarray

    if 'options_database' in source_dataset.variables:
        options_database = bytarr_to_dict(
            source_dataset['options_database'].values)
    else:
        options_database = None
    return options_database


def get_pref(source_dataset, ref_dataset,  options):
    '''
    Get reference pressure profile for source_dataset from ref_dataset
    or calculate from surface_press in source_dataset options_database
    and options['th_ref'].

    Parameters
    ----------
    source_dataset :  netCDF4 file
        MONC output file.
    ref_dataset :  netCDF4 file or None
        MONC output file containing pref
    options : dict

    Returns
    -------
    (pref, piref)

    '''

    if ref_dataset is None:
        od = options_database(source_dataset)
        if od is None:
            p_surf = 1.0E5
        else:
            p_surf = float(od['surface_pressure'])

        thref = options['th_ref']

        zn = source_dataset.variables['zn'][...]
        piref0 = (p_surf/1.0E5)**kappa
        piref = piref0 - (g/(cp_air * thref)) * zn
        pref = 1.0E5 * piref**rk
#                print('pref', pref)
    else:
        pref = ref_dataset.variables['prefn'][-1,...]
        piref = (pref[:]/1.0E5)**kappa

    return (pref, piref)

def get_thref(ref_dataset, options):
    '''
    Get thref profile from ref_dataset

    Parameters
    ----------
    ref_dataset : TnetCDF4 file or None
        MONC output file containing pref
    options : dict


    Returns
    -------
    thref : float or float array.
        Reference theta constant or profile

    '''
    if ref_dataset is None:
        thref = options['th_ref']
    else:
        thref = ref_dataset['thref']
        [itime] = find_var(thref.dims, ['time'])
        if itime is not None:
            tdim = thref.dims[itime]
            thref = thref[{tdim:[0]}]
        while len(np.shape(thref)) > 1:
            thref = thref.data[0,...]

    return thref

def find_var(vdims, var):
    '''
    Find dimensions containing strings.

    Parameters
    ----------
    vdims : xarray dimensions
    var : list of strings

    Returns
    -------
    tuple matching var with either index in vdims or None
    '''

    index_list = []
    for v in var :
        ind = None
        for i, vdim in enumerate(vdims):
            if v in vdim:
                ind = i
                break
        index_list.append(ind)
    return tuple(index_list)

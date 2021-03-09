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
import netCDF4

from netCDF4 import Dataset

import numpy as np
from scipy.signal import fftconvolve
from difference_ops import *

import time
from thermodynamics_constants import *

test_level = 1

subfilter_version = '0.3'


def filter_variable_list(source_dataset, ref_dataset, derived_dataset,
                         filtered_dataset, options, filter_def,
                         var_list=None, grid='p') :
    """
    Create filtered versions of input variables on required grid,
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
        default provided by get_default_variable_list()
        grid='p'        : Grid - 'u','v','w' or 'p'

    Returns:
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

        v = op_var['name']

        if v+"_r" not in filtered_dataset.variables \
            or v+"_s" not in filtered_dataset.variables:

            ncvar_r, ncvar_s = filter_field(op_var,
                                            filtered_dataset,
                                            options, filter_def,
                                            grid=grid, sync=False)

    derived_dataset.sync()
    filtered_dataset.sync()
    return var_list

def filter_variable_pair_list(source_dataset, ref_dataset, derived_dataset,
                              filtered_dataset, options, filter_def,
                              var_list=None, grid='p') :
    """
    Create filtered versions of pairs input variables on A grid,
    stored in derived_dataset.

    Args:
        source_dataset  : NetCDF dataset for input
        derived_dataset : NetCDF dataset for derived data
        filtered_dataset: NetCDF dataset for derived data
        options         : General options e.g. FFT method used.
        filter_def      : 1 or 2D filter
        var_list=None   : List of variable names.
        default provided by get_default_variable_pair_list()

    Returns:
        list : list of lists of pairs strings representing variable names.

    @author: Peter Clark

    """

    if (var_list==None):
        var_list = get_default_variable_pair_list()
        print("Default list:\n",var_list)
    if grid=='w' : zvar = "z"
    else : zvar = "zn"
    for v in var_list:
        print("Calculating s({},{})".format(v[0],v[1]))
        svar, vdims = quadratic_subfilter(source_dataset, ref_dataset,
                                  derived_dataset, filtered_dataset, options,
                                  filter_def, v[0], v[1], grid=grid)

        dims = find_var(vdims, ['x','y'])
        svar_name = f"s({v[0]:s},{v[1]:s})_on_{grid:s}"

        if svar_name not in filtered_dataset.variables:

            if filter_def.attributes['filter_type'] == 'domain' :
                edims = np.setdiff1d(np.arange(len(np.shape(svar))), dims)
                odims = [ ]
                for i in edims:
                    odims.append(vdims[i])
                ncsvar = filtered_dataset.createVariable(svar_name,"f8",
                                         tuple(odims))
            else :
                ncsvar = filtered_dataset.createVariable(svar_name,"f8",
                                                         vdims)


            ncsvar[...] = svar
            print(ncsvar)

    derived_dataset.sync()
    filtered_dataset.sync()
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

    vname = var['name']
    field = var['data']
    vdims = var['dims']

    sh = np.shape(field)
    ndims = len(sh)

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

        field_s = field[...] - field_r

    var_r = var.copy()
    var_r['name'] = var['name']+'_r'
    var_r['data'] = field_r

    var_s = var.copy()
    var_s['name'] = var['name']+'_s'
    var_s['data'] = field_s


    return [var_r, var_s]

def nc_dimcopy(source_dataset, derived_dataset, dimname) :
    """
    Copy dimension from source NetCDF dataset to destination

    Args:
        source_dataset
        derived_dataset
        dimname

    Returns:
        dimension

    @author: Peter Clark
    """
    v = source_dataset.variables[dimname]
    derived_dataset.createDimension(dimname,np.shape(v[:])[-1])
    dv = derived_dataset.createVariable(dimname,"f8",(dimname,))
    dv[:] = last_dim(v[:])
    return dv


def setup_data_file(source_file, ref_file, derived_dataset_name,
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

        derived_dataset = Dataset(derived_dataset_name, "a")

    else :

        exists = False
        derived_dataset = Dataset(derived_dataset_name, "w", clobber=True)

        source_dataset = Dataset(source_file,"r")
        w = source_dataset["w"]

        tvar = w.dimensions[0]

        times = nc_dimcopy(source_dataset, derived_dataset, tvar)

        z = nc_dimcopy(source_dataset, derived_dataset, "z")
        zn = nc_dimcopy(source_dataset, derived_dataset, "zn")

        if "x" in source_dataset.variables:
            x = nc_dimcopy(source_dataset, derived_dataset, "x")
        else:
            derived_dataset.createDimension("x",np.shape(w[:])[1])

        if "x" in source_dataset.variables:
            x = nc_dimcopy(source_dataset, derived_dataset, "x")
        else:
            derived_dataset.createDimension("y",np.shape(w[:])[2])

        derived_dataset.sync()
        source_dataset.close()

    return derived_dataset, exists

def setup_derived_data_file(source_file, destdir, ref_file, fname,
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
    derived_dataset_name = derived_dataset_name + "_" + fname + ".nc"

    derived_dataset, exists = setup_data_file(source_file, ref_file,
                    destdir+derived_dataset_name, override=override)

    return derived_dataset_name, derived_dataset, exists

def setup_filtered_data_file(source_file, destdir, ref_file, fname,
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
    filtered_dataset_name = filtered_dataset_name + "_" + fname + "_" + \
        filter_def.id + ".nc"
    filtered_dataset, exists = setup_data_file(source_file, ref_file,
                    destdir+filtered_dataset_name, override=override)

    filtered_dataset.filter_def_id = filter_def.id
    filtered_dataset.setncatts(filter_def.attributes)
    filtered_dataset.setncatts(options)

    filtered_dataset.sync()

    return filtered_dataset_name, filtered_dataset, exists

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
    var_properties = {"u":[True,False,False],
                      "v":[False,True,False],
                      "w":[False,False,True],
                      "th":[False,False,False],
                      "p":[False,False,False],
                      "q_vapour":[False,False,False],
                      "q_cloud_liquid_mass":[False,False,False],
                      }
    print(f'Retrieving {var_name:s}.')
    try :
        var = source_dataset[var_name]
        vardim = var.dimensions
        vard = var[...]
        varp = var_properties[var_name]

        print(vardim)
        if var_name == 'th' :
#            dims = find_var(vardim, ['x','y'])

            thref = get_thref(ref_dataset, options)
#            thref = np.expand_dims(thref, axis=dims)
            vard[...] += thref

    except :

        if var_name == 'th_L' :
#            rhoref = ref_dataset.variables['rhon'][-1,...]
            theta, vardim, varp = get_data(source_dataset, ref_dataset, 'th',
                                           options)
            (pref, piref) = get_pref(source_dataset, ref_dataset,  options)
            q_cl, vd, vp = get_data(source_dataset, ref_dataset,
                                    'q_cloud_liquid_mass', options)
            vard = theta - L_over_cp * q_cl / piref

        elif var_name == 'th_v' :
#            rhoref = ref_dataset.variables['rhon'][-1,...]
            theta, vardim, varp = get_data(source_dataset, ref_dataset, 'th',
                                           options)
            thref = get_thref(ref_dataset, options)
            q_v, vardim, varp = get_data(source_dataset, ref_dataset,
                                         'q_vapour', options)
            q_cl, vd, vp = get_data(source_dataset, ref_dataset,
                                    'q_cloud_liquid_mass', options)
            vard = theta + thref * (c_virtual * q_v - q_cl)

        elif var_name == 'q_total' :
#            rhoref = ref_dataset.variables['rhon'][-1,...]
            q_v, vardim, varp = get_data(source_dataset, ref_dataset,
                                         'q_vapour', options)
            q_cl, vd, vp = get_data(source_dataset, ref_dataset,
                                    'q_cloud_liquid_mass', options)
            vard = q_v + q_cl

        elif var_name == 'buoyancy':
            th_v, vardim, varp = get_data(source_dataset, ref_dataset, 'th_v',
                                           options)
            # get mean over horizontal axes
            haxes = find_var(vardim, ['x','y'])
            mean_thv = np.mean(th_v, axis = haxes)
            varp = grav * (th_v - mean_thv)/mean_thv

        else :

            sys.exit(f"Data {var_name:s} not in dataset.")


    return vard, vardim, varp

def get_and_transform(source_dataset, ref_dataset, var_name, options,
                      grid='p'):

    var, vdim, vp = get_data(source_dataset, ref_dataset, var_name,
                                           options)
    #    print(np.shape(var), vdim, vp)
    if grid=='p' :
        if vp[0] :
            print("Mapping {} from u grid to p grid.".format(var_name))
            var = field_on_u_to_p(var, xaxis=1)
        if vp[1] :
            print("Mapping {} from v grid to p grid.".format(var_name))
            var = field_on_v_to_p(var, xaxis=1)
        if vp[2] :
            print("Mapping {} from w grid to p grid.".format(var_name))
            z = source_dataset["z"]
            zn = source_dataset["zn"]
            var = field_on_w_to_p(var,  z, zn)

    elif grid=='u' :
        if not ( vp[0] or vp[1] or vp[2]):
            print("Mapping {} from p grid to u grid.".format(var_name))
            var = field_on_p_to_u(var, xaxis=1)
        if vp[1] :
            print("Mapping {} from v grid to u grid.".format(var_name))
            var = field_on_v_to_p(var, xaxis=1)
            var = field_on_p_to_u(var, xaxis=1)
        if vp[2] :
            print("Mapping {} from w grid to u grid.".format(var_name))
            z = source_dataset["z"]
            zn = source_dataset["zn"]
            var = field_on_w_to_p(var,  z, zn)
            var = field_on_p_to_u(var, xaxis=1)

    elif grid=='v' :
        if not ( vp[0] or vp[1] or vp[2]):
            print("Mapping {} from p grid to v grid.".format(var_name))
            var = field_on_p_to_v(var, yaxis=2)
        if vp[0] :
            print("Mapping {} from u grid to v grid.".format(var_name))
            var = field_on_u_to_p(var, xaxis=1)
            var = field_on_v_to_p(var, xaxis=1)
        if vp[2] :
            print("Mapping {} from w grid to v grid.".format(var_name))
            z = source_dataset["z"]
            zn = source_dataset["zn"]
            var = field_on_w_to_p(var,  z, zn)
            var = field_on_p_to_v(var, xaxis=1)

    elif grid=='w' :
        z = source_dataset["z"]
        zn = source_dataset["zn"]
        if not ( vp[0] or vp[1] or vp[2]):
            print("Mapping {} from p grid to w grid.".format(var_name))
            var = field_on_p_to_w(var, z, zn)
        if vp[0] :
            print("Mapping {} from u grid to w grid.".format(var_name))
            var = field_on_u_to_p(var, xaxis=1)
            var = field_on_p_to_w(var, z, zn)
        if vp[1] :
            print("Mapping {} from v grid to w grid.".format(var_name))
            var = field_on_v_to_p(var, xaxis=1)
            var = field_on_p_to_w(var,  z, zn)

        iz = find_var(vdim, ['z'])[0]
        if iz < len(vdim):
            vdim = list(vdim)
            vdim[iz] = 'z'
            vdim=tuple(vdim)
    else:
        print("Illegal grid ",grid)

    op_var = {'name' : var_name,
              'data' : var,
              'dims' : vdim,
              'grid_prop' : vp
              }

    return op_var

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

    var = None
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

        if op_var_name in derived_dataset.variables:
            var = derived_dataset[op_var_name]
            print(f'Retrieved {op_var_name:s} from derived dataset.')
            op_var['data'] = var[...]
            op_var['dims'] = var.dimensions
            op_var['grid_prop'] = vp


    if var is None:
        op_var = get_and_transform(source_dataset, ref_dataset,
                                   var_name, options, grid=grid)
        op_var['name'] = op_var_name

        if options['save_all'].lower() == 'yes':

            if grid=='w' : zvar = "z"
            else : zvar = "zn"

            ncvar = derived_dataset.createVariable(op_var_name,"f8",
                                                   op_var['dims'])

            ncvar[...] = op_var['data']
            print(f"Saved {op_var_name:s} to derived data file.")
            print(ncvar)

    return op_var

def deformation(source_dataset, dx, dy, grid='w') :
    """
    Read in 3D data from NetCDF file and, where necessary, interpolate to p grid.

    Assumes first dimension is time.

    Args:
        source_dataset  : NetCDF dataset
        var_name        : Name of variable

    Returns:
        data array.

    @author: Peter Clark
    """


    u = source_dataset["u"]
    v = source_dataset["v"]
    w = source_dataset["w"]
    z = last_dim(source_dataset["z"])
    zn = last_dim(source_dataset["zn"])

    vdims = u.dimensions
    xaxis = find_var(vdims, ['x'])[0]
#    print("u", np.shape(u),u.dimensions)
#    print("v", np.shape(v),v.dimensions)
#    print("w", np.shape(w),w.dimensions)

    ux = d_by_dx_field_on_u(u, dx, z, zn, xaxis=xaxis, grid = grid )
    uy = d_by_dy_field_on_u(u, dy, z, zn, xaxis=xaxis, grid = grid )
    uz = d_by_dz_field_on_u(u, z, zn, xaxis=xaxis, grid = grid )

    vx = d_by_dx_field_on_v(v, dx, z, zn, xaxis=xaxis, grid = grid )
    vy = d_by_dy_field_on_v(v, dy, z, zn, xaxis=xaxis, grid = grid )
    vz = d_by_dz_field_on_v(v, z, zn, xaxis=xaxis, grid = grid )


    wx = d_by_dx_field_on_w(w, dx, z, zn, xaxis=xaxis, grid = grid )
    wy = d_by_dy_field_on_w(w, dy, z, zn, xaxis=xaxis, grid = grid )
    wz = d_by_dz_field_on_w(w, z, zn, xaxis=xaxis, grid = grid )

    u_i = [ux, uy, uz]
    v_i = [vx, vy, vz]
    w_i = [wx, wy, wz]

    t = [u_i, v_i, w_i]

    op_var = {'name' : 'deformation',
              'dims' : vdims,
              'data' : t}
    return op_var


def filter_field(var, filtered_dataset, options, filter_def,
                 grid='p', sync=False) :
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
    vname = var['name']
    vard = var['data']
    vdims = var['dims']

    if grid=='w' : zvar = "z"
    else : zvar = "zn"

    itime = find_var(vdims, ['time'])
    dims = find_var(vdims, ['x','y'])

    print(f"Filtering {vname:s}")

    (var_r, var_s) = filtered_field_calc(var, options, filter_def)


    if vname+"_r" in filtered_dataset.variables:
        ncvar_r = filtered_dataset[vname+"_r"]
    else:
        if filter_def.attributes['filter_type'] == 'domain' :
            edims = np.setdiff1d(np.arange(len(np.shape(vard))), dims)
            odims = [ ]
            for i in edims:
                odims.append(vdims[i])
            ncvar_r = filtered_dataset.createVariable(vname+"_r", "f8",
                                                      tuple(odims))
        else :
            ncvar_r = filtered_dataset.createVariable(vname+"_r", "f8", vdims)

        ncvar_r[...] = var_r['data']

    print(f"Saved {ncvar_r.name:s}")
    print(ncvar_r)

    if vname+"_s" in filtered_dataset.variables:
        ncvar_s = filtered_dataset[vname+"_s"]
    else:
        ncvar_s = filtered_dataset.createVariable(vname+"_s", "f8", vdims)
        ncvar_s[...] = var_s['data']

    print(f"Saved {ncvar_s.name:s}")
    print(ncvar_s)

    if sync : filtered_dataset.sync()

    return (ncvar_r, ncvar_s)

def filtered_deformation(source_dataset, derived_dataset, filtered_dataset,
                         options, filter_def,
                         dx, dy, grid='p') :
    """
    Create filtered versions of deformation field.

    Args:
        source_dataset  : NetCDF input dataset
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

    d_var = deformation(source_dataset, dx, dy, grid=grid)
    d = d_var['data']
    vname = d_var['name']
    vdims = d_var['dims']
    d_ij_r = list()
    d_ij_s = list()
    for i in range(3) :
        d_j_r = list()
        d_j_s = list()
        for j in range(3) :
            var_in ={ 'name' : f"{vname}_{i:1d}_{j:1d}",
                      'data' : d[i][j],
                      'dims' : vdims}
            def_r, def_s = filter_field(var_in, filtered_dataset,
                                options, filter_def, grid=grid,
                                sync=False)
            d_j_r.append(def_r)
            d_j_s.append(def_s)
        d_ij_r.append(d_j_r)
        d_ij_s.append(d_j_s)

    filtered_dataset.sync()

    return d_ij_r, d_ij_s

def shear(d, no_trace=True) :
    trace = 0
    if no_trace :
        for i in range(3) :
            trace = trace + d[i][i][...]
        trace = (2.0/3.0) * trace
    S = {}
    mod_S = 0
    for i in range(3) :
        for j in range(i,3) :
            S_ij = d[i][j][...]+d[j][i][...]
            if i == j :
                S_ij = S_ij - trace
                mod_S += 0.5 * S_ij * S_ij
            else :
                mod_S +=  S_ij * S_ij
            S[f"{i:1d}_{j:1d}"] = (S_ij)
    return S, mod_S

def vorticity(d) :

    V_0 = d[2][1][...]-d[1][2][...]
    V_1 = d[0][3][...]-d[3][0][...]
    V_2 = d[1][0][...]-d[0][1][...]
    V = [V_0, V_1, V_2]

    return V

def quadratic_subfilter(source_dataset,  ref_dataset, derived_dataset,
                        filtered_dataset, options, filter_def,
                        v1_name, v2_name, grid='p') :
    """
    Create filtered versions of pair of input variables on required grid, stored in derived_dataset.
    Computes :math:`s(\phi,\psi) = (\phi\psi)^r - \phi^r\psi^r.`

    Args:
        source_dataset  : NetCDF dataset for input
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

    v1_name = v1['name']
    print("Reading ", v1_name+"_r")
    var1_r = filtered_dataset[v1_name+"_r"]

    vdims = var1_r.dimensions

#    print(np.shape(var1_r))

    v2 = get_data_on_grid(source_dataset,  ref_dataset,
                                       derived_dataset, v2_name, options,
                                       grid=grid)
    v2_name = v2['name']
    print("Reading ", v2_name+"_r")
    var2_r = filtered_dataset[v2_name+"_r"]
#    print(np.shape(var2_r))


    var1var2 = v1.copy()
    var1var2['name'] = v1['name'] + v2['name']

    var1var2['data'] = v1['data'] * v2['data']

    print(f"Filtering {v1_name:s}*{v2_name:s}")
    var1var2_r, var1var2_s = filtered_field_calc(var1var2, options,
                                                 filter_def )

    var1var2 = var1var2_r['data'] - var1_r[...] * var2_r[...]

    return var1var2, vdims

def bytarr_to_dict(d):
    res = {}
    for i in range(np.size(d[0])):
        opt = bytearray(d[i,0,:].compressed()).decode('utf-8')
        val = bytearray(d[i,1,:].compressed()).decode('utf-8')
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

    if 'options_database' in source_dataset.variables:
        options_database = bytarr_to_dict(
            source_dataset.variables['options_database'][...])
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
        if 'options_database' in source_dataset.variables:
            options_database = bytarr_to_dict(
                source_dataset.variables['options_database'][...])
            p_surf = float(options_database['surface_pressure'])
        else:
            p_surf = 1.0E5

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
        while len(np.shape(thref)) > 1:
            thref = thref[0,...]

    return thref

def find_var(vdims, var):
    index_list = []
    for v in var :
        for i, vdim in enumerate(vdims):
            if v in vdim:
                index_list.append(i)
                break
    return tuple(index_list)

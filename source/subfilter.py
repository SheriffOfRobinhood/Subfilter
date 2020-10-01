# -*- coding: utf-8 -*-
"""
Created on Tue Oct 23 11:07:05 2018

@author: Peter Clark
"""
import os
import netCDF4

from netCDF4 import Dataset

import numpy as np
from scipy.signal import fftconvolve
from difference_ops import *

import time
from thermodynamics_constants import *

test_level = 0


def filter_variable_list(source_dataset, ref_dataset, derived_dataset,
                         filtered_dataset, options, twod_filter,
                         var_list=None, grid='p') :
    """
    Create filtered versions of input variables on required grid,
    stored in derived_dataset.

    Args:
        source_dataset  : NetCDF dataset for input
        ref_dataset     : NetCDF dataset for input containing reference
                          profiles
        derived_dataset : NetCDF dataset for derived data
        filtered_dataset: NetCDF dataset for derived data
        options         : General options e.g. FFT method used.
        twod_filter     : 2D filter
        var_list=None   : List of variable names.
        default provided by get_default_variable_list()
        grid='p'        : Grid - 'u','v','w' or 'p'

    Returns:
        data array.

    @author: Peter Clark

    """

    if (var_list==None):
        var_list = get_default_variable_list()
        print("Default list:\n",var_list)
    for v in var_list:

        v, vard, vdims, varp  = get_data_on_grid(source_dataset, ref_dataset,
                                              derived_dataset, v, options,
                                              grid)

        if v+"_r" not in filtered_dataset.variables \
            or v+"_s" not in filtered_dataset.variables:

            ncvar_r, ncvar_s = filter_field(vard, v, vdims,
                                            filtered_dataset,
                                            options, twod_filter,
                                            grid=grid, three_d=True,
                                            sync=False)

    derived_dataset.sync()
    filtered_dataset.sync()
    return var_list

def filter_variable_pair_list(source_dataset, ref_dataset, derived_dataset,
                              filtered_dataset, options, twod_filter,
                              var_list=None, grid='p') :
    """
    Create filtered versions of pairs input variables on A grid, stored in derived_dataset.

    Args:
        source_dataset  : NetCDF dataset for input
        derived_dataset : NetCDF dataset for derived data
        filtered_dataset: NetCDF dataset for derived data
        options         : General options e.g. FFT method used.
        twod_filter     : 2D filter
        var_list=None   : List of variable names.
        default provided by get_default_variable_pair_list()

    Returns:
        var_list.

    @author: Peter Clark

    """

    if (var_list==None):
        var_list = get_default_variable_pair_list()
#                    ,"v","w","th","q_vapour","q_cloud_liquid_mass"]
        print("Default list:\n",var_list)
    if grid=='w' : zvar = "z"
    else : zvar = "zn"
    for v in var_list:
        print("Calculating s({},{})".format(v[0],v[1]))
#        var = source_dataset[v[0]]
#        vdims = var.dimensions
        svar, vdims = quadratic_subfilter(source_dataset, ref_dataset,
                                  derived_dataset, filtered_dataset, options,
                                  twod_filter, v[0], v[1], grid=grid)

        svar_name = "{}_{}_on_{}".format(v[0],v[1],grid)

        if svar_name not in filtered_dataset.variables:

            if twod_filter.attributes['filter_type'] == 'domain' :
                ncsvar = filtered_dataset.createVariable(svar_name,"f8",
                                         (vdims[0],zvar,))
            else :
                ncsvar = filtered_dataset.createVariable(svar_name,"f8",
                                         (vdims[0],vdims[1],vdims[2],zvar,))


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

    @author: Peter Clark


    """

    if test_level == 1:
# For testing
        var_list = [
            "u",
            "w",
            "th_L",
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
           var_list.

    @author: Peter Clark


    """
    if test_level == 1:
# For testing
        var_list = [
                ["w","th_L"],
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

def convolve(field, options, twod_filter, dims):
    """
    Convolve field filter using fftconvolve using padding.

    Args:
        field : 2D field array
        options         : General options e.g. FFT method used.
        twod_filter: 2D filter array

    Returns:
        field convolved with twod_filter

    @author: Peter Clark

    """
    if options['FFT_type'].upper() == 'FFTCONVOLVE':

        pad_len = np.max(np.shape(twod_filter))//2
        field = np.pad(field, pad_len, mode='wrap')
        result = fftconvolve(field, twod_filter, mode='same')
        result = result[pad_len:-pad_len, pad_len:-pad_len]

    elif options['FFT_type'].upper() == 'FFT':

        if len(np.shape(field)) > 2:
            edims = tuple(np.setdiff1d(np.arange(len(np.shape(field))), dims))
            twod_filter = np.expand_dims(twod_filter, axis=edims)

        fft_field = np.fft.fft2(field, axes=dims)

        fft_filtered_field = fft_field * twod_filter

        result = np.fft.ifft2(fft_filtered_field, axes=dims)
        result = result.real

    elif options['FFT_type'].upper() == 'RFFT':

        if len(np.shape(field)) > 2:
            edims = tuple(np.setdiff1d(np.arange(len(np.shape(field))), dims))
            twod_filter = np.expand_dims(twod_filter, axis=edims)

        fft_field = np.fft.rfft2(field, axes=dims)

        fft_filtered_field = fft_field * twod_filter

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


def filtered_field_calc(field, options, twod_filter, three_d=True ):
    """
    Split field into resolved (field_r) and subfilter (field_s).
    Note: this routine has a deliberate side effect, to store the fft or rfft
    of the filter in twod_filter for subsequent re-use.

    Args:
        field : ND numpy array
        options         : General options e.g. FFT method used.
        twod_filter: 2D filter
        three_d=True : if input has 3 dimensions, interpret as a single time.
                False: if input has 3 dimensions, interpret first as time.

    Returns:
        [field_r, field_s]

    @author: Peter Clark

    """

    sh = np.shape(field)
    ndims = len(sh)
    if twod_filter.attributes['filter_type'] == 'domain' :
        fshape = np.asarray(field.shape)
        if ndims == 2 :
            axis=(0,1)
            si = np.array([1,1])
            for i in range(2,len(fshape)) : si = np.concatenate(\
                          (si,np.array(fshape[i])),axis=None)
        else :
            axis=(1,2)
            si = np.array([fshape[0],1,1])
            for i in range(3,len(fshape)) : si = np.concatenate(\
                          (si,np.array(fshape[i])),axis=None)

        field_r = np.mean(field[...], axis=axis)
        field_s = field[...] - np.reshape(field_r,si)
    else :
        print("Filtering using {}".format(options['FFT_type']))
        if options['FFT_type'].upper() == 'FFTCONVOLVE':
            field_r = np.zeros(field.shape)
            dims = (0,1)
            if ndims == 2 :
                field_r[:,:] = convolve(field[:,:], options,
                                        twod_filter.data, dims)
            elif  ndims == 3 :
                if three_d :
                    for k in range(field.shape[2]) :
                        field_r[:,:,k] = convolve(field[:,:,k], options,
                                                  twod_filter.data, dims)
                else : # Assume third dimension is time
                    for k in range(field.shape[0]) :
                        field_r[k,:,:] = convolve(field[k,:,:], options,
                                                  twod_filter.data, dims)
            elif (ndims == 4):
                for it in range(field.shape[0]) :
                    for k in range(field.shape[3]) :
                        field_r[it,:,:,k] = convolve(field[it,:,:,k], options,
                                                     twod_filter.data, dims)

        elif options['FFT_type'].upper() == 'FFT':
            if ndims == 2 :
                dims = (0,1)
            elif ndims == 3 :
                if three_d :
                    dims = (0,1)
                else : # Assume third dimension is time
                    dims = (1,2)
            elif (ndims == 4):
                dims = (1,2)
            else:
                dims = (0,1)

            if 'fft' not in twod_filter.__dict__:
                sf = np.shape(twod_filter.data)
                if sh[dims[0]] != sf[0] or sh[dims[1]] != sf[1]:
                    padfilt = pad_to_len2D(twod_filter.data, sh[dims[0]])
                else:
                    padfilt = twod_filter.data.copy()
                # This shift of the filter is necessary to get the phase
                # information right.
                padfilt = np.fft.ifftshift(padfilt)
                twod_filter.fft = np.fft.fft2(padfilt)

            field_r = convolve(field, options, twod_filter.fft, dims)

        elif options['FFT_type'].upper() == 'RFFT':

            if ndims == 2 :
                dims = (0,1)
            elif ndims == 3 :
                if three_d :
                    dims = (0,1)
                else : # Assume third dimension is time
                    dims = (1,2)
            elif (ndims == 4):
                dims = (1,2)
            else:
                dims = (0,1)

            if 'rfft' not in twod_filter.__dict__:
                sf = np.shape(twod_filter.data)
                if sh[dims[0]] != sf[0] or sh[dims[1]] != sf[1]:
                    padfilt = pad_to_len2D(twod_filter.data, sh[dims[0]])
                else:
                    padfilt = twod_filter.data.copy()
                # This shift of the filter is necessary to get the phase
                # information right.
                padfilt = np.fft.ifftshift(padfilt)
                twod_filter.rfft = np.fft.rfft2(padfilt)

            field_r = convolve(field, options, twod_filter.rfft, dims)

        field_s = field[...] - field_r
    return [field_r, field_s]

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

    File name is original file name concatenated with twod_filter.id.

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
        ref_dataset=Dataset(ref_file)
        w = source_dataset["w"]
    #    print(w.dimensions)
    #    u = source_dataset["u"]

        tvar = w.dimensions[0]

        times = nc_dimcopy(source_dataset, derived_dataset, tvar)

        z = nc_dimcopy(ref_dataset, derived_dataset, "z")
        zn = nc_dimcopy(ref_dataset, derived_dataset, "zn")

        derived_dataset.createDimension("x",np.shape(w[:])[1])
        x = derived_dataset.createVariable("x","f8",("x",))
        derived_dataset.createDimension("y",np.shape(w[:])[2])
        y = derived_dataset.createVariable("y","f8",("y",))

        derived_dataset.sync()
        source_dataset.close()
        ref_dataset.close()

    return derived_dataset, exists

def setup_derived_data_file(source_file, destdir, ref_file, fname,
                            options, override=False) :
    """
    Create NetCDF dataset for derived data in destdir.

    File name is original file name concatenated with twod_filter.id.

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
                            options, twod_filter, override=False) :
    """
    Create NetCDF dataset for filtered data in destdir.

    File name is original file name concatenated with twod_filter.id.

    Args:
        source_file     : NetCDF file name.
        destdir         : Directory for derived data.
        options         : General options e.g. FFT method used.
        twod_filter     : Filter
        options         : General options e.g. FFT method used.
        override=False  : if True force creation of file

    Returns:
        filtered_dataset_name, filtered_dataset

    @author: Peter Clark
    """
    filtered_dataset_name = os.path.basename(source_file)
    filtered_dataset_name = ('.').join(filtered_dataset_name.split('.')[:-1])
    filtered_dataset_name = filtered_dataset_name + "_" + fname + "_" + \
        twod_filter.id + ".nc"
    filtered_dataset, exists = setup_data_file(source_file, ref_file,
                    destdir+filtered_dataset_name, override=override)

    filtered_dataset.twod_filter_id = twod_filter.id
    filtered_dataset.setncatts(twod_filter.attributes)
    filtered_dataset.setncatts(options)

    filtered_dataset.sync()

    return filtered_dataset_name, filtered_dataset, exists

def get_data(source_dataset, ref_dataset, var_name) :
    """
    Extract data from source NetCDF dataset or derived data.

	Currently supported derived data are:
	'th_L'    : Liquid water potential temperature.
	'th_v'    : Virtual potential temperature.
	'q_total' : Total water

    Returns:
    variable, variable_dimensions, variable_grid_properties

    @author: Peter Clark

    """

    var_properties = {"u":[True,False,False],
                      "v":[False,True,False],
                      "w":[False,False,True],
                      "th":[False,False,False],
                      "q_vapour":[False,False,False],
                      "q_cloud_liquid_mass":[False,False,False],
                      }
    print(var_name)
    try :
        var = source_dataset[var_name]
        vardim = var.dimensions
        vard = var[...]
        varp = var_properties[var_name]

        print(vardim)
        if var_name == 'th' :
            thref = ref_dataset['thref']
            for it in range(np.shape(vard[...])[0]) :
                vard[it,...] += thref[it,...]
    except :
        print("Data {} not in dataset".format(var_name))
        if var_name == 'th_L' :
#            rhoref = ref_dataset.variables['rhon'][-1,...]
            theta, vardim, varp = get_data(source_dataset, ref_dataset, 'th')
            pref = ref_dataset.variables['prefn'][-1,...]
            piref = (pref[:]/1.0E5)**kappa
            q_cl, vd, vp = get_data(source_dataset, ref_dataset,
                                    'q_cloud_liquid_mass')
#            print("Delta theta_L",np.max(L_over_cp * q_cl / piref))
            vard = theta - L_over_cp * q_cl / piref
#            input("Press enter")

        if var_name == 'th_v' :
#            rhoref = ref_dataset.variables['rhon'][-1,...]
            theta, vardim, varp = get_data(source_dataset, ref_dataset, 'th')
            thref = ref_dataset['thref'][-1,...]
            q_v, vardim, varp = get_data(source_dataset, ref_dataset,
                                         'q_vapour')
            q_cl, vd, vp = get_data(source_dataset, ref_dataset,
                                    'q_cloud_liquid_mass')
#            print("Delta theta_L",np.max(L_over_cp * q_cl / piref))
            vard = theta + thref * (c_virtual * q_v - q_cl)

        if var_name == 'q_total' :
#            rhoref = ref_dataset.variables['rhon'][-1,...]
            q_v, vardim, varp = get_data(source_dataset, ref_dataset,
                                         'q_vapour')
            q_cl, vd, vp = get_data(source_dataset, ref_dataset,
                                    'q_cloud_liquid_mass')
            vard = q_v + q_cl

        if var_name == 'buoyancy':
            th_v, vardim, varp = get_data(source_dataset, ref_dataset, 'th_v')
            mean_thv = np.mean(th_v, axis = (1,2))
            varp = grav * (th_v - mean_thv)/mean_thv

    return vard, vardim, varp

def get_and_transform(source_dataset, ref_dataset, var_name, grid='p'):

    var, vdim, vp = get_data(source_dataset, ref_dataset, var_name)
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
            z = ref_dataset["z"]
            zn = ref_dataset["zn"]
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
            z = ref_dataset["z"]
            zn = ref_dataset["zn"]
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
            z = ref_dataset["z"]
            zn = ref_dataset["zn"]
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
    else:
        print("Illegal grid ",grid)

    return var, vdim, vp

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

    if options['save_all'].lower() == 'yes':
        if op_var_name in derived_dataset.variables:
            var = derived_dataset[op_var_name]
            print('Retrieved {} from derived dataset.'.format(op_var_name))
            vdims = var.dimensions

    if var is None:
        var, vdims, vp = get_and_transform(source_dataset, ref_dataset,
                                          var_name, grid=grid)

        if options['save_all'].lower() == 'yes':

            if grid=='w' : zvar = "z"
            else : zvar = "zn"

            ncvar = derived_dataset.createVariable(op_var_name,"f8",
                                     (vdims[0],vdims[1],vdims[2],zvar,))

            ncvar[...] = var[...]
            print("Saved {} to derived data file.".format(op_var_name))
            print(ncvar)

    vard = var[...]
    return op_var_name, vard, vdims, vp

def deformation(source_dataset, dx, dy, z, zn, xaxis=0, grid='w') :
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
    print("u", np.shape(u),u.dimensions)
    print("v", np.shape(v),v.dimensions)
    print("w", np.shape(w),w.dimensions)

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
    return t

def filter_field(vard, vname, vdims, filtered_dataset, options, twod_filter,
                 grid='p', three_d=True, sync=False) :
    """
    Create filtered versions of input variable on required grid, stored in filtered_dataset.

    Args:
        filtered_dataset : NetCDF dataset for derived data.
        options         : General options e.g. FFT method used.
        twod_filter     : 2D filter.
        vard            : data array.
        vname           : variable name.
        vdims           : array of variable dimension names.
        default provided by get_default_variable_list()
        grid='p'        : Grid - 'u','v','w' or 'p'

    Returns:
        ncvar_r, ncvar_s: Resolved and subfilter fields as netcf variables in
                          filtered_dataset.

    @author: Peter Clark

    """
    if grid=='w' : zvar = "z"
    else : zvar = "zn"
    print(vname)

    var_r, var_s = filtered_field_calc(vard, options, twod_filter,
                                       three_d=True )

    if vname+"_r" in filtered_dataset.variables:
        ncvar_r = filtered_dataset[vname+"_r"]
    else:
        if twod_filter.attributes['filter_type'] == 'domain' :
            ncvar_r = filtered_dataset.createVariable(vname+"_r","f8",
                                     (vdims[0],zvar,))
        else :
            ncvar_r = filtered_dataset.createVariable(vname+"_r","f8",
                                     (vdims[0],vdims[1],vdims[2],zvar,))

        ncvar_r[...] = var_r
    print(ncvar_r)

    if vname+"_s" in filtered_dataset.variables:
        ncvar_s = filtered_dataset[vname+"_s"]
    else:
        ncvar_s = filtered_dataset.createVariable(vname+"_s","f8",
                                     (vdims[0],vdims[1],vdims[2],zvar,))
        ncvar_s[...] = var_s
    print(ncvar_s)

    if sync : filtered_dataset.sync()

    return ncvar_r, ncvar_s

def filtered_deformation(source_dataset, derived_dataset, filtered_dataset,
                         options, twod_filter,
                         dx, dy, z, zn, xaxis=0, grid='p') :
    """
    Create filtered versions of deformation field.

    Args:
        source_dataset  : NetCDF input dataset
        derived_dataset : NetCDF dataset for derived data.
        filtered_dataset: NetCDF dataset for derived data
        options         : General options e.g. FFT method used.
        twod_filter     : 2D filter.
        grid='p'        : Grid - 'u','v','w' or 'p'

    Returns:
        ncvar_r, ncvar_s: Resolved and subfilter fields as netcf variables in filtered_dataset.

    @author: Peter Clark

    """

#:math:`\frac{\partial u_i}{\partial{x_j}`
    u = source_dataset["u"]
    vdims = u.dimensions
    vname = 'deformation'

    d = deformation(source_dataset, dx, dy, z, zn, xaxis=xaxis,
                           grid=grid)
    d_ij_r = list()
    d_ij_s = list()
    for i in range(3) :
        d_j_r = list()
        d_j_s = list()
        for j in range(3) :
            vn = f"{vname}_{i:1d}_{j:1d}"
            def_r, def_s = filter_field(d[i][j], vn, vdims, filtered_dataset,
                                options, twod_filter, grid='p',
                                three_d=True, sync=False)
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
                        filtered_dataset, options, twod_filter,
                        v1_name, v2_name, grid='p') :
    """
    Create filtered versions of pair of input variables on required grid, stored in derived_dataset.

    Args:
        source_dataset  : NetCDF dataset for input
        derived_dataset : NetCDF dataset for derived data
        filtered_dataset: NetCDF dataset for derived data
        options         : General options e.g. FFT method used.
        twod_filter     : 2D filter
        v1_name, v2_name: Variable names.

    Returns:
        s(var1,var2) data array.
        vdims dimensions of var1


    @author: Peter Clark

    """

    v1_name, vard1, vd1, vp1 = get_data_on_grid(source_dataset,  ref_dataset,
                                       derived_dataset, v1_name, options,
                                       grid=grid)
    print("Reading ", v1_name+"_r")
    var1_r = filtered_dataset[v1_name+"_r"]

    vdims = var1_r.dimensions

#    print(np.shape(var1_r))

    v2_name, vard2, vd2, vp2 = get_data_on_grid(source_dataset,  ref_dataset,
                                       derived_dataset, v2_name, options,
                                       grid=grid)
    print("Reading ", v2_name+"_r")
    var2_r = filtered_dataset[v2_name+"_r"]
#    print(np.shape(var2_r))

    var1var2 = vard1 * vard2
    var1var2_r, var1var2_s = filtered_field_calc(var1var2, options,
                                                 twod_filter, three_d=True )

    var1var2 = var1var2_r - var1_r[...] * var2_r[...]

    return var1var2, vdims


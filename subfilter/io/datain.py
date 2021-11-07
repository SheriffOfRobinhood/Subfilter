# -*- coding: utf-8 -*-
"""
Created on Mon Aug  2 11:01:11 2021

@author: paclk
"""
import sys
import numpy as np
from .MONC_utils import options_database
from ..utils.string_utils import get_string_index
from ..utils.dask_utils import re_chunk
from .dataout import save_field
import subfilter.utils.difference_ops as do
import xarray as xr
import subfilter.thermodynamics.thermodynamics_constants as thc
import config

def get_data(source_dataset, ref_dataset, var_name, options) :
    """
    Extract data or derived data field from source NetCDF dataset.

    Currently written for MONC data, enforcing C-grid. Returned coords are
    ''x_p', 'x_u',, y_p', 'y_v', 'z', 'zn'. Coordinate x- and -y values are
    retrieved from the MONC options_database in source_dataset
    or from 'dx' and 'dy' in options otherwise

    Currently supported derived data are::

        'th_L'     : Liquid water potential temperature.
        'th_v'     : Virtual potential temperature.
        'q_total'  : Total water.
        'buoyancy' : (g/mean_th_v)*(th_v-mean_th_v), where the mean is the domain mean.

    Parameters
    ----------
    source_dataset : xarray Dataset
        Input (at least 2D) data.
    ref_dataset :  xarray Dataset
        Contains reference profiles. Can be None.
    var_name : str
        Name of variable to retrieve.
    options : dict
        Options Options possibly used are 'dx' and 'dy'.

    Returns
    -------
        vard: xarray
            Output data field.

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

        [itime] = get_string_index(vard.dims, ['time'])
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

        else:

            if 'x' in vard.dims:
                nx = vard.shape[vard.get_axis_num('x')]
                x = np.arange(nx) * np.float64(od['dxx'])
                xn = 'x_p'
                vard = vard.rename({'x':xn})
                vard.coords[xn] = x

            if 'y' in vard.dims:
                ny = vard.shape[vard.get_axis_num('y')]
                y = np.arange(ny) * np.float64(od['dyy'])
                yn = 'y_p'
                vard = vard.rename({'y':yn})
                vard.coords[yn] = y

#        print(vard)
        if var_name == 'th' :
            thref = get_thref(ref_dataset, options)
            vard += thref

    except :

        if var_name == 'th_L' :
            theta = get_data(source_dataset, ref_dataset, 'th',
                                           options)
            (pref, piref) = get_pref(source_dataset, ref_dataset,  options)
            q_cl = get_data(source_dataset, ref_dataset,
                                    'q_cloud_liquid_mass', options)
            vard = theta - thc.L_over_cp * q_cl / piref


        elif var_name == 'th_v' :
            theta = get_data(source_dataset, ref_dataset, 'th',
                                           options)
            thref = get_thref(ref_dataset, options)
            q_v = get_data(source_dataset, ref_dataset,
                                         'q_vapour', options)
            q_cl = get_data(source_dataset, ref_dataset,
                                    'q_cloud_liquid_mass', options)
            vard = theta + thref * (thc.c_virtual * q_v - q_cl)

        elif var_name == 'q_total' :
            q_v = get_data(source_dataset, ref_dataset,
                                         'q_vapour', options)
            q_cl = get_data(source_dataset, ref_dataset,
                                    'q_cloud_liquid_mass', options)
            vard = q_v + q_cl

        elif var_name == 'buoyancy':
            th_v = get_data(source_dataset, ref_dataset, 'th_v',
                                           options)
            # get mean over horizontal axes
            mean_thv = th_v.mean(dim=('x','y'))
            vard = thc.grav * (th_v - mean_thv)/mean_thv

        else :

            sys.exit(f"Data {var_name:s} not in dataset.")
#    print(vard)

    return vard

def get_and_transform(source_dataset, ref_dataset, var_name, options,
                      grid='p'):
    """
    Extract data or derived data from source NetCDF dataset and transform
    to alternative grid.

    See get_data for derived variables.

    Parameters
    ----------
    source_dataset : xarray Dataset
        Input (at least 2D) data.
    ref_dataset : xarray Dataset
        Contains reference profiles. Can be None.
    var_name : str
        Name of variable to retrieve.
    options : dict
        Options. Options possibly used are 'dx' and 'dy'.
    grid : str, optional
        Destination grid 'u', 'v', 'w' or 'p'. Default is 'p'.

    Returns
    -------
        var: xarray
            Output data field.

    @author: Peter Clark

    """
    var = get_data(source_dataset, ref_dataset, var_name, options)
    z = source_dataset["z"]
    zn = source_dataset["zn"]
    var = do.grid_conform(var, z, zn, grid = grid )

#    print(var)

    vp = ['x_u' in var.dims,
          'y_v' in var.dims,
          'z' in var.dims]


    # if grid=='p' :
    #     if vp[0] :
    #         print("Mapping {} from u grid to p grid.".format(var_name))
    #         var = do.field_on_u_to_p(var)
    #     if vp[1] :
    #         print("Mapping {} from v grid to p grid.".format(var_name))
    #         var = do.field_on_v_to_p(var)
    #     if vp[2] :
    #         print("Mapping {} from w grid to p grid.".format(var_name))
    #         z = source_dataset["z"]
    #         zn = source_dataset["zn"]
    #         var = do.field_on_w_to_p(var, zn)

    # elif grid=='u' :
    #     if not ( vp[0] or vp[1] or vp[2]):
    #         print("Mapping {} from p grid to u grid.".format(var_name))
    #         var = do.field_on_p_to_u(var)
    #     if vp[1] :
    #         print("Mapping {} from v grid to u grid.".format(var_name))
    #         var = do.field_on_v_to_p(var)
    #         var = do.field_on_p_to_u(var)
    #     if vp[2] :
    #         print("Mapping {} from w grid to u grid.".format(var_name))
    #         z = source_dataset["z"]
    #         zn = source_dataset["zn"]
    #         var = do.field_on_w_to_p(var, zn)
    #         var = do.field_on_p_to_u(var)

    # elif grid=='v' :
    #     if not ( vp[0] or vp[1] or vp[2]):
    #         print("Mapping {} from p grid to v grid.".format(var_name))
    #         var = do.field_on_p_to_v(var)
    #     if vp[0] :
    #         print("Mapping {} from u grid to v grid.".format(var_name))
    #         var = do.field_on_u_to_p(var)
    #         var = do.field_on_p_to_v(var)
    #     if vp[2] :
    #         print("Mapping {} from w grid to v grid.".format(var_name))
    #         z = source_dataset["z"]
    #         zn = source_dataset["zn"]
    #         var = do.field_on_w_to_p(var, zn)
    #         var = do.field_on_p_to_v(var)

    # elif grid=='w' :
    #     z = source_dataset["z"]
    #     zn = source_dataset["zn"]
    #     if not ( vp[0] or vp[1] or vp[2]):
    #         print("Mapping {} from p grid to w grid.".format(var_name))
    #         var = do.field_on_p_to_w(var, z)
    #     if vp[0] :
    #         print("Mapping {} from u grid to w grid.".format(var_name))
    #         var = do.field_on_u_to_p(var)
    #         var = do.field_on_p_to_w(var, z)
    #     if vp[1] :
    #         print("Mapping {} from v grid to w grid.".format(var_name))
    #         var = do.field_on_v_to_p(var)
    #         var = do.field_on_p_to_w(var, z)

    # else:
    #     print("Illegal grid ",grid)
    # print(var)

#    print(zvar)

    if not config.dask_opts['no_dask']:
        var = re_chunk(var)
#    print(var)

    return var

def get_data_on_grid(source_dataset, ref_dataset, derived_dataset, var_name,
                     options, grid='p') :
    """
    Find data from source_dataset remapped to destination grid.
    Uses data from derived_dataset if present, otherwise uses
    get_and_transform to input from source_dataset and remap grid.
    In this case, if options['save_all']=='yes', save teh remapped data to
    derived_dataset.

    See get_data for derived variables.

    Parameters
    ----------
    source_dataset : xarray Dataset
        Input (at least 2D) data.
    ref_dataset : xarray Dataset
        Contains reference profiles. Can be None.
    derived_dataset : dict
        'ds' points to xarray Dataset, 'file' to output file path.
    var_name : str
        Name of variable to retrieve.
    options : dict
        Options. Options possibly used are 'dx' and 'dy'.
    grid : str, optional
        Destination grid 'u', 'v', 'w' or 'p'. Default is 'p'.

    Returns
    -------
        var: xarray
            Output data field.

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
    # then if x matches required output grid, see if in derived_dataset
    # already, and use if it is.
    # Otherwise strip '_on_x' and go back to source data as per default.

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

def get_pref(source_dataset, ref_dataset,  options):
    """
    Get reference pressure profile for source_dataset.

    Calculate from ref_dataset or from surface_press in source_dataset
    options_database and options['th_ref'].

    Parameters
    ----------
    source_dataset :  netCDF4 file
        MONC output file.
    ref_dataset :  netCDF4 file or None
        MONC output file containing 1D variable prefn.
    options : dict
        Options. Options possibly used are th_ref.

    Returns
    -------
    (pref, piref)

    """
    if ref_dataset is None:
        od = options_database(source_dataset)
        if od is not None:
            p_surf = float(od['surface_pressure'])
        else:
            p_surf = thc.p_ref_theta

        thref = options['th_ref']

        zn = source_dataset.variables['zn'][...]
        piref0 = (p_surf/thc.p_ref_theta)**thc.kappa
        piref = piref0 - (thc.g/(thc.cp_air * thref)) * zn
        pref = thc.p_ref_theta * piref**thc.rk
#                print('pref', pref)
    else:
        pref = ref_dataset.variables['prefn'][-1,...]
        piref = (pref[:]/thc.p_ref_theta)**thc.kappa

    return (pref, piref)

def get_thref(ref_dataset, options):
    """
    Get thref profile from ref_dataset.

    Parameters
    ----------
    ref_dataset : TnetCDF4 file or None
        MONC output file containing pref
    options : dict
        Options. Options possibly used are th_ref.

    Returns
    -------
    thref : float or float array.
        Reference theta constant or profile

    """
    if ref_dataset is None:
        thref = options['th_ref']
    else:
        thref = ref_dataset['thref']
        [itime] = get_string_index(thref.dims, ['time'])
        if itime is not None:
            tdim = thref.dims[itime]
            thref = thref[{tdim:[0]}]
        while len(np.shape(thref)) > 1:
            thref = thref.data[0,...]

    return thref

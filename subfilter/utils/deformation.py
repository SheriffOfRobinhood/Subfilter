# -*- coding: utf-8 -*-
"""
Created on Mon Aug  2 11:29:37 2021

@author: paclk
"""
import numpy as np
import xarray as xr

import subfilter.utils.difference_ops as do

from .string_utils import get_string_index
from ..io.datain import get_data
from .dask_utils import re_chunk
from .. import subfilter_setup
from ..io.dataout import save_field, setup_child_file


def deformation(source_dataset, ref_dataset, derived_dataset,
                options, grid='w') :
    r"""
    Compute deformation tensor.

    Deformation tensor is defined as :math:`{\partial u_{i}}/{\partial {x_j}}`.

    Parameters
    ----------
        source_dataset  : NetCDF dataset
            Inout data.
        ref_dataset     : NetCDF dataset
            Input data for input containing reference
            profiles. Can be None
        derived_dataset : NetCDF dataset
            Output dataset for derived data.
        options         : dict
            General options e.g. FFT method used.
        grid : str
            destination grid (Default = 'w')

    Returns
    -------
        xarray
            Array with new dimensions 'i' and 'j'.
            Saves to derived_dataset if options['save_all'].lower() == 'yes'.

    @author: Peter Clark
    """
    if 'deformation' in derived_dataset['ds']:
        deformation = derived_dataset['ds']['deformation']
        return deformation

    u = get_data(source_dataset, ref_dataset, "u", options)
    [iix, iiy, iiz] = get_string_index(u.dims, ['x', 'y', 'z'])

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
#    print(ux)

    uy = do.d_by_dy_field_on_u(u, z, grid = grid )
#    print(uy)

    uz = do.d_by_dz_field_on_u(u, z, grid = grid )
#    print(uz)

    vx = do.d_by_dx_field_on_v(v, z, grid = grid )
#    print(vx)

    vy = do.d_by_dy_field_on_v(v, z, grid = grid )
#    print(vy)

    vz = do.d_by_dz_field_on_v(v, z, grid = grid )
#    print(vz)

    wx = do.d_by_dx_field_on_w(w, zn, grid = grid )
#    print(wx)

    wy = do.d_by_dy_field_on_w(w, zn, grid = grid )
#    print(wy)

    wz = do.d_by_dz_field_on_w(w, zn, grid = grid )
#    print(wz)

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

def shear(d, no_trace:bool=True) :
    r"""
    Compute shear tensor from deformation tensor.

    Shear tensor is defined by
    :math:`S_{ij}={\partial u_{i}}/{\partial {x_j}}
    +{\partial u_{j}}/{\partial {x_i}}`.

    Parameters
    ----------
    d : xarray
        Deformation tensor with dimensions 'i' and 'j'.
    no_trace : bool, optional
        If true, subtract trace (divergence). The default is True.

    Returns
    -------
    S : xarray
        Shear with new dimension 'i_j'.
    mod_S_sq : xarray

        Modulus of shear squared. :math:`1/2 S_{ij}S_{ij}`

    """
    trace = 0
    vname = ''
    if no_trace :
        for k in range(3) :
            trace = trace + d.isel(i=k, j=k)
        # 2 arises because its the trace of S_ij we want.
        trace = (2.0/3.0) * trace
        vname = 'n'

    mod_S_sq = 0

    S = []
    i_j = []
    for k in range(3) :
        for l in range(k,3) :

            S_kl = d.isel(i=k, j=l) + d.isel(i=l, j=k)
            # Note normalisation.
            # We are calculating S_ij = du_i/dx_j + du_j/dx_i
            # |S|^2 = 0.5 S_ij S_ij
            # However We are only calculating ij terms, not ji
            # for efficiency so we include 0.5 S_ji S_ji implicitly
            # by losing the 0.5 in off-diagonal terms.
            if k == l :
                S_kl = S_kl - trace
                mod_S_sq += 0.5 * S_kl * S_kl
            else :
                mod_S_sq +=  S_kl * S_kl

            S.append(S_kl)
            i_j.append(f'{k:d}_{l:d}')

    S = xr.concat(S, dim='i_j', coords='minimal', compat='override')
    S.coords['i_j'] = i_j
    S.name = 'shear' + vname
    S.attrs={'units':'s-1'}
    S = re_chunk(S)

    mod_S_sq.name = 'mod_S_sq' + vname
    mod_S_sq.attrs={'units':'s-2'}
    S = re_chunk(S)

    return S, mod_S_sq

def vorticity(d):
    r"""
    Compute vorticity vector from deformation tensor.

    vorticity vector is defined as
    :math:`\omega_{k}=({\partial u_{i}}/{\partial {x_j}}
    -{\partial u_{j}}/{\partial {x_i}})`
    with :math:`kij` defined cyclically (i.e. 123, 231, 312).

    Parameters
    ----------
    d : xarray
        Deformation tensor with dimensions 'i' and 'j'.

    Returns
    -------
    vorticity : xarray
        Vorticity with dimension 'i'.

    """
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

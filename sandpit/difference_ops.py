# -*- coding: utf-8 -*-
"""
difference_ops.py

Created on Wed Apr 17 21:03:43 2019

Difference operators for C-grid data.

@author: Peter Clark
"""
import numpy as np
from .dask_utils import re_chunk
from .string_utils import get_string_index
import dask
import xarray
#use_map_overlap = False
use_map_overlap = True

grid_def = { 'p':('x_p', 'y_p', 'zn'),
             'u':('x_u', 'y_p', 'zn'),
             'v':('x_p', 'y_v', 'zn'),
             'w':('x_p', 'y_p', 'z')}


def exec_fn(fn, field: xarray.DataArray, axis: int) -> xarray.DataArray:
    """
    Execute function using map_overlap with overlap on selected axis.

    Parameters
    ----------
    fn : function
        DESCRIPTION.
    field : xarray.DataArray
        DESCRIPTION.
    axis : int
        DESCRIPTION.

    Returns
    -------
    new xarray.DataArray

    """
    if use_map_overlap:
        print('Using map_overlap.')
        d = field.data.map_overlap(fn, depth={axis:1},
                                  boundary={axis:'periodic'})
        field.data = d
    else:
        sh = np.shape(field)
        ch = {field.dims[axis]:sh[axis]}
        field = re_chunk(field, chunks=ch)
        field = fn(field)
    return field


def last_dim(z) :
    """
    Remove all but last dimension of z.

    Parameters
        z : n-dimensional array.

    Returns:
        z[0,0, etc. ,:]
    @author: Peter Clark
    """

    zd = z[...]
    while len(np.shape(zd))>1 :
        zd = zd[0,...]
    return zd

def interpolate(field, znew) :
    """
    Interpolate field from z to zn

    Parameters
        field : xarray nD field
        znew  : xarray coordinate new z.

    Returns:
        field on zn levels
    @author: Peter Clark
    """
    [zaxis] = get_string_index(field.dims,['z'])
    zdim = field.dims[zaxis]

    newfield = field.interp({zdim:znew}, kwargs={"fill_value": "extrapolate"})
    newfield = newfield.drop_vars(zdim)

    return newfield

def grid_conform_x(field, target_xdim):
    [xaxis] = get_string_index(field.dims,['x'])
    xdim = field.dims[xaxis]
    if xdim == target_xdim:
        print(f'xdim is already {target_xdim}')
        return field
    x = field.coords[xdim].data
    dx = x[1] - x[0]
    if target_xdim == 'x_p':
        xmn = lambda arr:(0.5 * (arr + np.roll(arr, -1, axis=xaxis)))
        x_new = x - dx / 2.0
    elif target_xdim == 'x_u':
        xmn = lambda arr:(0.5 * (arr + np.roll(arr, +1, axis=xaxis)))
        x_new = x + dx / 2.0
    else:
        print(f"Cannot transform {xdim} to {target_xdim}")
        return field

    print(f'{xdim} to {target_xdim}')
    newfield = field.rename({xdim:target_xdim})
    newfield = exec_fn(xmn, newfield, xaxis)
    newfield.coords[target_xdim] = x_new
    return newfield

def grid_conform_y(field, target_ydim):
    [yaxis] = get_string_index(field.dims,['y'])
    ydim = field.dims[yaxis]
    if ydim == target_ydim:
        print(f'ydim is already {target_ydim}')
        return field
    y = field.coords[ydim].data
    dy = y[1] - y[0]
    if target_ydim == 'y_p':
        ymn = lambda arr:(0.5 * (arr + np.roll(arr, -1, axis=yaxis)))
        y_new = y - dy / 2.0
    elif target_ydim == 'y_v':
        ymn = lambda arr:(0.5 * (arr + np.roll(arr, +1, axis=yaxis)))
        y_new = y + dy / 2.0
    else:
        print(f"Cannot transform {ydim} to {target_ydim}")
        return field

    print(f'{ydim} to {target_ydim}')
    newfield = field.rename({ydim:target_ydim})
    newfield = exec_fn(ymn, newfield, yaxis)
    newfield.coords[target_ydim] = y_new
    return newfield

def grid_conform_z(field, z, zn, target_zdim):
    [zaxis] = get_string_index(field.dims,['z'])
    zdim = field.dims[zaxis]
    if zdim == target_zdim:
        print(f'zdim is already {target_zdim}')
        return field
    elif target_zdim == 'z':
        print(f'{zdim} to {target_zdim}')
        return interpolate(field, z)
    elif target_zdim == 'zn':
        print(f'{zdim} to {target_zdim}')
        return interpolate(field, zn)
    else:
        print(f"Cannot transform {zdim} to {target_zdim}")
        return field

def grid_conform(field, z, zn, grid: str = 'p' ):
    if type(grid) == str:
        if grid in ['p', 'u', 'v', 'w']:
            op_grid = grid_def[grid]
        else:
            raise ValueError(f"grid={grid} is illegal value.")
    else:
        op_grid = grid

    newfield = grid_conform_x(field, op_grid[0])
    newfield = grid_conform_y(newfield, op_grid[1])
    newfield = grid_conform_z(newfield, z, zn, op_grid[2])
    return newfield

def d_by_dx_field(field, z, zn, grid: str = 'p' ) :
    """
    Differentiate field in x direction.

    Parameters
    ----------
        field : xarray nD field
        z: xarray coordinate
            zcoord on w levels - needed if changing vertical grid.
        zn: xarray coordinate
            zcoord on p levels - needed if changing vertical grid.
        grid : str | tuple of 2 strings
            destination grid (Default = 'p')

    Returns
    -------
        field on required grid

    @author: Peter Clark
    """
    [xaxis] = get_string_index(field.dims,['x'])
    xdim = field.dims[xaxis]
    x = field.coords[xdim].data
    dx = x[1] - x[0]
    if xdim == 'x_u':
        print("d_by_dx_field_on_x_u ", grid)
        xdim_new = 'x_p'
        xdrv = lambda arr:((arr - np.roll(arr,  1, axis=xaxis)) / dx)
        x_new = x - dx / 2.0
    else:
        if xdim != 'x_p':
            print(f"d_by_dx_field on unknown grid {xdim}, assuming x_p.")
        print("d_by_dx_field_on_x_p ",grid)
        xdim_new = 'x_u'
        xdrv = lambda arr:((np.roll(arr, -1, axis=xaxis) - arr) / dx)
        x_new = x + dx / 2.0

    newfield = field.rename({xdim:xdim_new})
    newfield = exec_fn(xdrv, newfield, xaxis)
    newfield.coords[xdim_new] = x_new
    newfield = grid_conform(newfield, z, zn, grid=grid)
    newfield.name = f"d{field.name:s}_by_dx_on_{grid:s}"

    return newfield

def d_by_dy_field(field, z, zn, grid: str = 'p' ) :
    """
    Differentiate field in y direction.

    Parameters
    ----------
        field : xarray nD field
        z: xarray coordinate
            zcoord on w levels - needed if changing vertical grid.
        zn: xarray coordinate
            zcoord on p levels - needed if changing vertical grid.
        grid : str | tuple of 2 strings
            destination grid (Default = 'p')

    Returns
    -------
        field on required grid

    @author: Peter Clark
    """
    [yaxis] = get_string_index(field.dims,['y'])
    ydim = field.dims[yaxis]
    y = field.coords[ydim].data
    dy = y[1] - y[0]
    if ydim == 'y_v':
        print("d_by_dy_field_on_y_v ", grid)
        ydim_new = 'y_p'
        ydrv = lambda arr:((arr - np.roll(arr,  1, axis=yaxis)) / dy)
        y_new = y - dy / 2.0
    else:
        if ydim != 'y_p':
            print(f"d_by_dy_field on unknown grid {ydim}, assuming y_p.")
        print("d_by_dy_field_on_y_p ",grid)
        ydim_new = 'y_v'
        ydrv = lambda arr:((np.roll(arr, -1, axis=yaxis) - arr) / dy)
        y_new = y + dy / 2.0

    newfield = field.rename({ydim:ydim_new})
    newfield = exec_fn(ydrv, newfield, yaxis)
    newfield.coords[ydim_new] = y_new
    newfield = grid_conform(newfield, z, zn, grid=grid)
    newfield.name = f"d{field.name:s}_by_dy_on_{grid:s}"

    return newfield

def d_by_dz_field(field, z, zn, grid: str = 'p'):
    """
    Differentiate field in z direction.

    Parameters
    ----------
        field : xarray nD field
        z: xarray coordinate
            zcoord on w levels - needed if changing vertical grid.
        zn: xarray coordinate
            zcoord on p levels - needed if changing vertical grid.
        grid : str | tuple of 2 strings
            destination grid (Default = 'p')

    Returns
    -------
        field on required grid

    @author: Peter Clark
    """
    [zaxis] = get_string_index(field.dims,['z'])
    zdim = field.dims[zaxis]
    zcoord = field.coords[zdim].data
    if zdim == 'zn':
        print("d_by_dz_field_on_zn ")
        zdim_new = 'zi'
        pad = (0,1)
        nroll = -1
        (exn, dexn) = (-1, -1)
    else:
        print("d_by_dz_field_on_z ")
        zdim_new = 'zn'
        pad = (1,0)
        nroll = 1
        (exn, dexn) = (0, 1)

    newfield = field.diff(zdim)/field.coords[zdim].diff(zdim)
    newfield = newfield.pad(pad_width={zdim:pad}, mode = 'edge')
    newfield = newfield.rename({zdim:zdim_new})
    zi = 0.5 * (zcoord + np.roll(zcoord, nroll))
    zi[exn] = 2 * zi[exn + dexn] - zi[exn + 2 * dexn]
    newfield.coords[zdim_new] = zi
    newfield = grid_conform(newfield, z, zn, grid=grid)
    newfield.name = f"d{field.name:s}_by_dz_on_{grid:s}"

    return newfield

def field_on_w_to_p(field, znew) :
    print("w_to_p")
    return interpolate(field, znew)

def field_on_p_to_w(field, znew) :
    print("p_to_w")
    return interpolate(field, znew)

def field_on_u_to_p(field) :
    """
    Interpolate field from u to p points in C grid

    Parameters
        field : nD xarray field

    Returns:
        field on p points
    @author: Peter Clark
    """

    print("u_to_p")
    x = field.coords['x_u'].data
    xaxis = field.get_axis_num('x_u')
    xmn = lambda arr:(0.5 * (arr + np.roll(arr,-1,axis=xaxis)))
    newfield = field.rename({'x_u':'x_p'})
    newfield = exec_fn(xmn, newfield, xaxis)
    newfield.coords['x_p'] = x - x[0]
    return newfield

def field_on_p_to_u(field) :
    """
    Interpolate field from p to u points in C grid

    Parameters
        field : nD xarray field

    Returns:
        field on p points
    @author: Peter Clark
    """

    print("p_to_u")
    x = field.coords['x_p'].data
    xaxis = field.get_axis_num('x_p')
    xmn = lambda arr:(0.5 * (arr + np.roll(arr,+1,axis=xaxis)))
    newfield = field.rename({'x_p':'x_u'})
    newfield = exec_fn(xmn, newfield, xaxis)
    newfield.coords['x_u'] = x + (x[1] - x[0]) / 2.0
    return newfield

def field_on_v_to_p(field) :
    """
    Interpolate field from v to p points in C grid

    Parameters
        field : nD xarray  field

    Returns:
        field on p points
    @author: Peter Clark
    """

    print("v_to_p")
    y = field.coords['y_v'].data
    yaxis = field.get_axis_num('y_v')
    ymn = lambda arr:(0.5 * (arr + np.roll(arr,-1,axis=yaxis)))
    newfield = field.rename({'y_v':'y_p'})
    newfield = exec_fn(ymn, newfield, yaxis)
    newfield.coords['y_p'] = y - y[0]
    return newfield

def field_on_p_to_v(field) :
    """
    Interpolate field from p to v points in C grid

    Parameters
        field : nD xarray field

    Returns:
        field on p points
    @author: Peter Clark
    """
    print("p_to_v")
    y = field.coords['y_p'].data
    yaxis = field.get_axis_num('y_p')
    ymn = lambda arr:(0.5 * (arr + np.roll(arr,1,axis=yaxis)))
    newfield = field.rename({'y_p':'y_v'})
    newfield = exec_fn(ymn, newfield, yaxis)
    newfield.coords['y_v'] = y + (y[1] - y[0]) / 2.0
    return newfield

def d_by_dz_field_on_zn(field):
    """
    Differentiate field on zn levels in z direction.

    Parameters
        field : xarray nD field

    Returns:
        field on required grid
    @author: Peter Clark
    """
    print("d_by_dz_field_on_zn ")
    zn = field.coords['zn']
    new = field.diff('zn')/field.coords['zn'].diff('zn')
    new = new.pad(pad_width={'zn':(0,1)}, mode = 'edge')
    new = new.rename({'zn':'zi'})
    zi = 0.5 * (zn.data + np.roll(zn.data, -1))
    zi[-1] = 2 * zn.data[-2] - zn.data[-3]
    new.coords['zi'] = zi

    return new

def d_by_dz_field_on_z(field):
    """
    Differentiate field on z levels in z direction.

    Parameters
        field : xarray nD field

    Returns:
        field on required grid
    @author: Peter Clark
    """
    print("d_by_dz_field_on_z ")

    z = field.coords['z']
    new = field.diff('z')/field.coords['z'].diff('z')
    new = new.pad(pad_width={'z':(1,0)}, mode = 'edge')
    new = new.rename({'z':'zn'})
    zi = 0.5 * (z.data + np.roll(z.data, 1))
    zi[0] = 2 * z.data[1] - z.data[2]
    new.coords['zn'] = zi

    return new

def d_by_dx_field_on_u(field, z, grid: str = 'p' ) :
    """
    Differentiate field on u points in x direction then average to req grid.

    Parameters
    ----------
        field : xarray
            nD field
        z: xarray coordinate
            zcoord - needed if changing vertical grid.
        grid : str
            destination grid (Default = 'p')

    Returns
    -------
        field on required grid
    @author: Peter Clark
    """
    print("d_by_dx_field_on_u ", grid)
    x = field.coords['x_u'].data
    dx = x[1] - x[0]
    xaxis = field.get_axis_num('x_u')
    xdrv = lambda arr:((arr - np.roll(arr,1,axis=xaxis)) / dx)
    newfield = field.rename({'x_u':'x_p'})
    newfield = exec_fn(xdrv, newfield, xaxis)
    newfield.coords['x_p'] = x - dx / 2.0

    # Derivative on p
    if grid == 'u' :
        newfield = field_on_p_to_u(newfield)
    if grid == 'v' :
        newfield = field_on_p_to_v(newfield)
    if grid == 'w' :
        newfield = field_on_p_to_w(newfield, z)

    newfield.name = f"d{field.name:s}_by_dx_on_{grid:s}"

    return newfield

def d_by_dy_field_on_u(field, z, grid: str = 'p' ) :
    """
    Differentiate field on u points in y direction then average to req grid

    Parameters
    ----------
        field : xarray
            nD field
        z: xarray coordinate
            zcoord - needed if changing vertical grid.
        grid : str
            destination grid (Default = 'p')

    Returns:
        field on required grid
    @author: Peter Clark
    """

    print("d_by_dy_field_on_u ",grid)
    y = field.coords['y_p'].data
    dy = y[1] - y[0]
    yaxis = field.get_axis_num('y_p')
    ydrv = lambda arr:((np.roll(arr,-1,axis=yaxis) - arr) / dy)
    newfield = field.rename({'y_p':'y_v'})
    newfield = exec_fn(ydrv, newfield, yaxis)

    newfield.coords['y_v'] = y + dy / 2.0

    # Derivative on u,v
    if grid == 'p' :
        newfield = field_on_u_to_p(newfield)
        newfield = field_on_v_to_p(newfield)
    if grid == 'u' :
        newfield = field_on_v_to_p(newfield)
    if grid == 'v' :
        newfield = field_on_u_to_p(newfield)
    if grid == 'w' :
        newfield = field_on_u_to_p(newfield)
        newfield = field_on_v_to_p(newfield)
        newfield = field_on_p_to_w(newfield, z)

    newfield.name = f"d{field.name:s}_by_dy_on_{grid:s}"

    return newfield

def d_by_dz_field_on_u(field, z, grid: str = 'p' ) :
    """
    Differentiate field on u points in z direction then average to req grid

    Parameters
    ----------
        field : xarray
            nD field
        z: xarray coordinate
            zcoord - needed if changing vertical grid.
        grid : str
            destination grid (Default = 'p')

    Returns:
        field on required grid
    @author: Peter Clark
    """

    zn = field.coords['zn']
    newfield = d_by_dz_field_on_zn(field)

    # Derivative on u,w
    if grid == 'p' :
        newfield = field_on_w_to_p(newfield, zn)
        newfield = field_on_u_to_p(newfield)
    if grid == 'u' :
        newfield = field_on_w_to_p(newfield, zn)
    if grid == 'v' :
        newfield = field_on_w_to_p(newfield, zn)
        newfield = field_on_u_to_p(newfield)
        newfield = field_on_p_to_v(newfield)
    if grid == 'w' :
        newfield = newfield.interp({'zi':z},
                                   kwargs={"fill_value": "extrapolate"})
        newfield = newfield.drop_vars('zi')
        newfield = field_on_u_to_p(newfield)

    newfield.name = f"d{field.name:s}_by_dz_on_{grid:s}"

    return newfield

def d_by_dx_field_on_v(field, z, grid: str = 'p' ) :

    """
    Differentiate field on u points in y direction then average to req grid

    Parameters
    ----------
        field : xarray
            nD field
        z: xarray coordinate
            zcoord - needed if changing vertical grid.
        grid : str
            destination grid (Default = 'p')

    Returns:
        field on required grid
    @author: Peter Clark
    """

    print("d_by_dx_field_on_v ",grid)
    x = field.coords['x_p'].data
    dx = x[1] - x[0]
    xaxis = field.get_axis_num('x_p')
    xdrv = lambda arr:((np.roll(arr,-1,axis=xaxis) - arr)/dx)
    newfield = field.rename({'x_p':'x_u'})
    newfield = exec_fn(xdrv, newfield, xaxis)

    newfield.coords['x_u'] = x + dx / 2.0

    # Derivative on u,v
    if grid == 'p' :
        newfield = field_on_u_to_p(newfield)
        newfield = field_on_v_to_p(newfield)
    if grid == 'u' :
        newfield = field_on_v_to_p(newfield)
    if grid == 'v' :
        newfield = field_on_u_to_p(newfield)
    if grid == 'w' :
        newfield = field_on_u_to_p(newfield)
        newfield = field_on_v_to_p(newfield)
        newfield = field_on_p_to_w(newfield, z)

    newfield.name = f"d{field.name:s}_by_dx_on_{grid:s}"

    return newfield

def d_by_dy_field_on_v(field, z, grid = 'p' ) :
    """
    Differentiate field on u points in x direction then average to req grid

    Parameters
    ----------
        field : xarray
            nD field
        z: xarray coordinate
            zcoord - needed if changing vertical grid.
        grid : str
            destination grid (Default = 'p')

    Returns:
        field on required grid
    @author: Peter Clark
    """

    print("d_by_dy_field_on_v ",grid)
    y = field.coords['y_v'].data
    dy = y[1] - y[0]
    yaxis = field.get_axis_num('y_v')
    ydrv = lambda arr:((arr - np.roll(arr,1,axis=yaxis)) / dy)
    newfield = field.rename({'y_v':'y_p'})
    newfield = exec_fn(ydrv, newfield, yaxis)

    newfield.coords['y_p'] = y - dy / 2.0

    # Derivative on p
    if grid == 'u' :
        newfield = field_on_p_to_u(newfield)
    if grid == 'v' :
        newfield = field_on_p_to_v(newfield)
    if grid == 'w' :
        newfield = field_on_p_to_w(newfield, z)

    newfield.name = f"d{field.name:s}_by_dy_on_{grid:s}"

    return newfield

def d_by_dz_field_on_v(field, z, grid = 'p' ) :
    """
    Differentiate field on u points in z direction then average to req grid

    Parameters
    ----------
        field : xarray
            nD field
        z: xarray coordinate
            zcoord - needed if changing vertical grid.
        grid : str
            destination grid (Default = 'p')

    Returns:
        field on required grid
    @author: Peter Clark
    """

    print("d_by_dz_field_on_v ",grid)

    zn = field.coords['zn']
    newfield = d_by_dz_field_on_zn(field)

    # Derivative on v,w
    if grid == 'p' :
        newfield = field_on_w_to_p(newfield, zn)
        newfield = field_on_u_to_p(newfield)
    if grid == 'u' :
        newfield = field_on_w_to_p(newfield, zn)
    if grid == 'v' :
        newfield = field_on_w_to_p(newfield, zn)
        newfield = field_on_u_to_p(newfield)
        newfield = field_on_p_to_v(newfield)
    if grid == 'w' :
        newfield = newfield.interp({'zi':z},
                                   kwargs={"fill_value": "extrapolate"})
        newfield = newfield.drop_vars('zi')
        newfield = field_on_v_to_p(newfield)

    newfield.name = f"d{field.name:s}_by_dz_on_{grid:s}"

    return newfield

def d_by_dx_field_on_p(field, z, grid = 'p' ) :
    """
    Differentiate field on u points in x direction then average to req grid

    Parameters
    ----------
        field : xarray
            nD field
        z: xarray coordinate
            zcoord - needed if changing vertical grid.
        grid : str
            destination grid (Default = 'p')

    Returns:
        field on required grid
    @author: Peter Clark
    """
    print("d_by_dx_field_on_p ",grid)
    x = field.coords['x_p'].data
    dx = x[1] - x[0]
    xaxis = field.get_axis_num('x_p')
    xdrv = lambda arr:((np.roll(arr,-1,axis=xaxis) - arr) / dx)
    newfield = field.rename({'x_p':'x_u'})
    newfield = exec_fn(xdrv, newfield, xaxis)

    newfield.coords['x_u'] = x + dx / 2.0

    # Derivative on u
    if grid == 'p' :
        newfield = field_on_u_to_p(newfield)
    if grid == 'v' :
        newfield = field_on_u_to_p(newfield)
        newfield = field_on_p_to_v(newfield)
    if grid == 'w' :
        newfield = field_on_u_to_p(newfield)
        newfield = field_on_p_to_w(newfield, z)

    newfield.name = f"d{field.name:s}_by_dx_on_{grid:s}"

    return newfield

def d_by_dy_field_on_p(field, z, grid = 'p' ) :
    """
    Differentiate field on u points in y direction then average to req grid

    Parameters
    ----------
        field : xarray
            nD field
        z: xarray coordinate
            zcoord - needed if changing vertical grid.
        grid : str
            destination grid (Default = 'p')

    Returns:
        field on required grid
    @author: Peter Clark
    """
    print("d_by_dy_field_on_p ",grid)
    y = field.coords['y_p'].data
    dy = y[1] - y[0]
    yaxis = field.get_axis_num('y_p')
    ydrv = lambda arr:((np.roll(arr,-1,axis=yaxis) - arr) / dy)
    newfield = field.rename({'y_p':'y_v'})
    newfield = exec_fn(ydrv, newfield, yaxis)

    newfield.coords['y_v'] = y + dy / 2.0

    # Derivative on v
    if grid == 'p' :
        newfield = field_on_v_to_p(newfield)
    if grid == 'u' :
        newfield = field_on_v_to_p(newfield)
        newfield = field_on_p_to_u(newfield)
    if grid == 'w' :
        newfield = field_on_v_to_p(newfield)
        newfield = field_on_p_to_w(newfield, z)

    newfield.name = f"d{field.name:s}_by_dy_on_{grid:s}"

    return newfield

def d_by_dz_field_on_p(field, z, grid = 'p' ) :
    """
    Differentiate field on u points in z direction then average to req grid

    Parameters
    ----------
        field : xarray
            nD field
        z: xarray coordinate
            zcoord - needed if changing vertical grid.
        grid : str
            destination grid (Default = 'p')

    Returns:
        field on required grid
    @author: Peter Clark
    """

    print("d_by_dz_field_on_p ",grid)
    zn = field.coords['zn']
    newfield = d_by_dz_field_on_zn(field)

    # Derivative on w
    if grid == 'p' :
        newfield = field_on_w_to_p(newfield, zn)
    if grid == 'u' :
        newfield = field_on_w_to_p(newfield, zn)
        newfield = field_on_p_to_u(newfield)
    if grid == 'v' :
        newfield = field_on_w_to_p(newfield, zn)
        newfield = field_on_p_to_v(newfield)
    if grid == 'w' :
        newfield = newfield.interp({'zi':z},
                                   kwargs={"fill_value": "extrapolate"})
        newfield = newfield.drop_vars('zi')

    newfield.name = f"d{field.name:s}_by_dz_on_{grid:s}"

    return newfield

def d_by_dx_field_on_w(field, zn, grid = 'p' ) :
    """
    Differentiate field on u points in x direction then average to req grid

    Parameters
    ----------
        field : xarray
            nD field
        zn: xarray coordinate
            zcoord - needed if changing vertical grid.
        grid : str
            destination grid (Default = 'p')

    Returns:
        field on required grid
    @author: Peter Clark
    """

    print("d_by_dx_field_on_w ",grid)
    x = field.coords['x_p'].data
    dx = x[1] - x[0]
    xaxis = field.get_axis_num('x_p')
    xdrv = lambda arr:((np.roll(arr,-1,axis=xaxis) - arr) / dx)
    newfield = field.rename({'x_p':'x_u'})
    newfield = exec_fn(xdrv, newfield, xaxis)

    newfield.coords['x_u'] = x + dx / 2.0

    # Derivative on u,w
    if grid == 'p' :
        newfield = field_on_w_to_p(newfield, zn)
        newfield = field_on_u_to_p(newfield)
    if grid == 'u' :
        newfield = field_on_w_to_p(newfield, zn)
    if grid == 'v' :
        newfield = field_on_w_to_p(newfield, zn)
        newfield = field_on_u_to_p(newfield)
        newfield = field_on_p_to_v(newfield)
    if grid == 'w' :
        newfield = field_on_u_to_p(newfield)

    newfield.name = f"d{field.name:s}_by_dx_on_{grid:s}"

    return newfield

def d_by_dy_field_on_w(field, zn, grid = 'p' ) :
    """
    Differentiate field on u points in y direction then average to req grid

    Parameters
    ----------
        field : xarray
            nD field
        zn: xarray coordinate
            zcoord - needed if changing vertical grid.
        grid : str
            destination grid (Default = 'p')

    Returns:
        field on required grid
    @author: Peter Clark
    """

    print("d_by_dy_field_on_w ",grid)
    y = field.coords['y_p'].data
    dy = y[1] - y[0]
    yaxis = field.get_axis_num('y_p')
    ydrv = lambda arr:((np.roll(arr,-1,axis=yaxis) - arr) / dy)
    newfield = field.rename({'y_p':'y_v'})
    newfield = exec_fn(ydrv, newfield, yaxis)

    newfield.coords['y_v'] = y + dy / 2.0

    # Derivative on v,w
    if grid == 'p' :
        newfield = field_on_w_to_p(newfield, zn)
        newfield = field_on_v_to_p(newfield)
    if grid == 'u' :
        newfield = field_on_w_to_p(newfield, zn)
        newfield = field_on_v_to_p(newfield)
        newfield = field_on_p_to_u(newfield)
    if grid == 'v' :
        newfield = field_on_w_to_p(newfield, zn)
    if grid == 'w' :
        newfield = field_on_v_to_p(newfield)

    newfield.name = f"d{field.name:s}_by_dy_on_{grid:s}"

    return newfield

def d_by_dz_field_on_w(field, zn, grid = 'p' ) :
    """
    Differentiate field on u points in z direction then average to req grid

    Parameters
    ----------
        field : xarray
            nD field
        zn: xarray coordinate
            zcoord - needed if changing vertical grid.
        grid : str
            destination grid (Default = 'p')

    Returns:
        field on required grid
    @author: Peter Clark
    """

    print("d_by_dz_field_on_w ",grid)

    z = field.coords['z']
    newfield = d_by_dz_field_on_z(field)

    # Derivative on p
    if grid == 'p' :
        newfield = newfield.interp({'zi':zn},
                                   kwargs={"fill_value": "extrapolate"})
        newfield = newfield.drop_vars('zi')
    if grid == 'u' :
        newfield = field_on_p_to_u(newfield)
    if grid == 'v' :
        newfield = field_on_p_to_v(newfield)
    if grid == 'w' :
        newfield = field_on_p_to_w(newfield, z)

    newfield.name = f"d{field.name:s}_by_dz_on_{grid:s}"

    return newfield

def padleft(f, zt, axis=0) :
    """
    Add dummy field at bottom of nD array

    Parameters
        f : nD field
        zt: 1D zcoordinates
        axis=0: Specify axis to extend

    Returns:
        extended field, extended coord
    @author: Peter Clark
    """

    s = list(np.shape(f))
    s[axis] += 1
    newfield = np.zeros(s)
    newfield[...,1:]=f
    newz = np.zeros(np.size(zt)+1)
    newz[1:] = zt
    newz[0] = 2*zt[0]-zt[1]
    return newfield, newz

def padright(f, zt, axis=0) :
    """
    Add dummy field at top of nD array

    Parameters
        f : nD field
        zt: 1D zcoordinates
        axis=0: Specify axis to extend

    Returns:
        extended field, extended coord
    @author: Peter Clark
    """

    s = list(np.shape(f))
    s[axis] += 1
    newfield = np.zeros(s)
    newfield[...,:-1] = f
    newz = np.zeros(np.size(zt)+1)
    newz[:-1] = zt
    newz[-1] = 2*zt[-1]-zt[-2]
    return newfield, newz

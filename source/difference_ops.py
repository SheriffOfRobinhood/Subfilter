# -*- coding: utf-8 -*-
"""

  difference_ops.py

Created on Wed Apr 17 21:03:43 2019

Difference operators for C-grid data.

@author: Peter Clark
"""
import numpy as np

def last_dim(z) :
    """
    Remove all but last dimension of z.

    Args:
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

    Args:
        field : xarray nD field
        znew  : xarray coordinate new z.

    Returns:
        field on zn levels
    @author: Peter Clark
    """
    if 'z' in field.dims:
        zvar = 'z'
    elif 'zn' in field.dims:
        zvar = 'zn'

    newfield = field.interp({zvar:znew}, kwargs={"fill_value": "extrapolate"})
    newfield = newfield.drop_vars(zvar)

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

    Args:
        field : nD xarray field

    Returns:
        field on p points
    @author: Peter Clark
    """

    print("u_to_p")
    d = field.data
    x = field.coords['x_u'].data
    xaxis = field.get_axis_num('x_u')
    xmn = lambda arr:(0.5 * (arr + np.roll(arr,-1,axis=xaxis)))
    newfield = field.rename({'x_u':'x_p'})
    newfield.data = d.map_overlap(xmn, depth={'x_p':(1)}, 
                                  boundary={xaxis:'periodic'})
    newfield.coords['x_p'] = x - x[0]
    return newfield

def field_on_p_to_u(field) :
    """
    Interpolate field from p to u points in C grid

    Args:
        field : nD xarray field

    Returns:
        field on p points
    @author: Peter Clark
    """

    print("p_to_u")
    d = field.data
    x = field.coords['x_p'].data
    xaxis = field.get_axis_num('x_p')
    xmn = lambda arr:(0.5 * (arr + np.roll(arr,+1,axis=xaxis)))
    newfield = field.rename({'x_p':'x_u'})
    newfield.data = d.map_overlap(xmn, depth={'x_u':(1)}, 
                                  boundary={xaxis:'periodic'})
    newfield.coords['x_u'] = x + (x[1] - x[0]) / 2.0
    return newfield

def field_on_v_to_p(field) :
    """
    Interpolate field from v to p points in C grid

    Args:
        field : nD xarray  field

    Returns:
        field on p points
    @author: Peter Clark
    """

    print("v_to_p")
    d = field.data
    y = field.coords['y_v'].data
    yaxis = field.get_axis_num('y_v')
    ymn = lambda arr:(0.5 * (arr + np.roll(arr,-1,axis=yaxis)))
    newfield = field.rename({'y_v':'y_p'})
    newfield.data = d.map_overlap(ymn, depth={'y_p':(1)}, 
                                  boundary={yaxis:'periodic'})
    newfield.coords['y_p'] = y - y[0]
    return newfield

def field_on_p_to_v(field) :
    """
    Interpolate field from p to v points in C grid

    Args:
        field : nD xarray field

    Returns:
        field on p points
    @author: Peter Clark
    """
    print("p_to_v")
    d = field.data
    y = field.coords['y_p'].data
    yaxis = field.get_axis_num('y_p')
    ymn = lambda arr:(0.5 * (arr + np.roll(arr,1,axis=yaxis)))
    newfield = field.rename({'y_p':'y_v'})
    newfield.data = d.map_overlap(ymn, depth={'y_v':(1)}, 
                                  boundary={yaxis:'periodic'})
    newfield.coords['y_v'] = y + (y[1] - y[0]) / 2.0
    return newfield

def d_by_dz_field_on_zn(field):
    """
    Differentiate field on zn levels in z direction.

    Args:
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

    Args:
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
    zi[0] = 2 * z.data[2] - z.data[3]
    new.coords['zn'] = zi

    return new

def d_by_dx_field_on_u(field, z, grid = 'p' ) :
    """
    Differentiate field on u points in x direction then average to req grid

    Args:
        field : xarray nD field
        z: zcoord - needed if changing vertical grid to w.
        grid = 'p': destination grid

    Returns:
        field on required grid
    @author: Peter Clark
    """

    print("d_by_dx_field_on_u ",grid)
    d = field.data
    x = field.coords['x_u'].data
    dx = x[1] - x[0]
    xaxis = field.get_axis_num('x_u')
    xdrv = lambda arr:((arr - np.roll(arr,1,axis=xaxis)) / dx)
    newfield = field.rename({'x_u':'x_p'})
    newfield.data = d.map_overlap(xdrv, depth={'x_p':(1)}, 
                                  boundary={xaxis:'periodic'})
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

def d_by_dy_field_on_u(field, z, grid = 'p' ) :
    """
    Differentiate field on u points in y direction then average to req grid

    Args:
        field : xarray nD field
        z: zcoord - needed if changing vertical grid.
        grid = 'p': destination grid

    Returns:
        field on required grid
    @author: Peter Clark
    """

    print("d_by_dy_field_on_u ",grid)
    d = field.data
    y = field.coords['y_p'].data
    dy = y[1] - y[0]
    yaxis = field.get_axis_num('y_p')
    ydrv = lambda arr:((np.roll(arr,-1,axis=yaxis) - arr) / dy)
    newfield = field.rename({'y_p':'y_v'})
    newfield.data = d.map_overlap(ydrv, depth={'y_v':(1)}, 
                                  boundary={yaxis:'periodic'})
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

def d_by_dz_field_on_u(field, z, grid = 'p' ) :
    """
    Differentiate field on u points in z direction then average to req grid

    Args:
        field : xarray nD field
        z: zcoord - needed if changing vertical grid.
        grid = 'p': destination grid

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

def d_by_dx_field_on_v(field, z, grid = 'p' ) :
    """
    Differentiate field on u points in y direction then average to req grid

    Args:
        field : nD field
        z: zcoord - needed if changing vertical grid.
        grid = 'p': destination grid

    Returns:
        field on required grid
    @author: Peter Clark
    """

    print("d_by_dx_field_on_v ",grid)
    d = field.data
    x = field.coords['x_p'].data
    dx = x[1] - x[0]
    xaxis = field.get_axis_num('x_p')
    xdrv = lambda arr:((np.roll(arr,-1,axis=xaxis) - arr)/dx)
    newfield = field.rename({'x_p':'x_u'})
    newfield.data = d.map_overlap(xdrv, depth={'x_u':(1)}, 
                                  boundary={xaxis:'periodic'})
    newfield.coords['x_u'] = x - dx / 2.0

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

    Args:
        field : xarray nD field
        z: zcoord - needed if changing vertical grid.
        grid = 'p': destination grid

    Returns:
        field on required grid
    @author: Peter Clark
    """

    print("d_by_dy_field_on_v ",grid)
    d = field.data
    y = field.coords['y_v'].data
    dy = y[1] - y[0]
    yaxis = field.get_axis_num('y_v')
    ydrv = lambda arr:((arr - np.roll(arr,1,axis=yaxis)) / dy)
    newfield = field.rename({'y_v':'y_p'})
    newfield.data = d.map_overlap(ydrv, depth={'y_p':(1)}, 
                                  boundary={yaxis:'periodic'})
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

    Args:
         field : nD xarray field
         z: zcoord - needed if changing vertical grid.
         grid = 'p': destination grid

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

    Args:
        field : nD xarray field
        z: zcoord - needed if changing vertical grid.
       grid = 'p': destination grid

    Returns:
        field on required grid
    @author: Peter Clark
    """
    print("d_by_dx_field_on_p ",grid)
    d = field.data
    x = field.coords['x_p'].data
    dx = x[1] - x[0]
    xaxis = field.get_axis_num('x_p')
    xdrv = lambda arr:((np.roll(arr,-1,axis=xaxis) - arr) / dx)
    newfield = field.rename({'x_p':'x_u'})
    newfield.data = d.map_overlap(xdrv, depth={'x_u':(1)}, 
                                  boundary={xaxis:'periodic'})
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

    Args:
        field : nD xarray field
        z: zcoord - needed if changing vertical grid.
        grid = 'p': destination grid

    Returns:
        field on required grid
    @author: Peter Clark
    """
    print("d_by_dy_field_on_p ",grid)
    d = field.data
    y = field.coords['y_p'].data
    dy = y[1] - y[0]
    yaxis = field.get_axis_num('y_p')
    ydrv = lambda arr:((np.roll(arr,-1,axis=yaxis) - arr) / dy)
    newfield = field.rename({'y_p':'y_v'})
    newfield.data = d.map_overlap(ydrv, depth={'y_v':(1)}, 
                                  boundary={yaxis:'periodic'})
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

    Args:
        field : nD xarray field
        z: zcoord - needed if changing vertical grid.
        grid = 'p': destination grid

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

    Args:
        field : nD field
        zn: zcoord - needed if changing vertical grid.
        grid = 'p': destination grid

    Returns:
        field on required grid
    @author: Peter Clark
    """

    print("d_by_dx_field_on_w ",grid)
    d = field.data
    x = field.coords['x_p'].data
    dx = x[1] - x[0]
    xaxis = field.get_axis_num('x_p')
    depths = np.zeros(len(field.dims))
    depths[xaxis] = 1
    newfield = field.rename({'x_p':'x_u'})
    xdrv = lambda arr:((np.roll(arr,-1,axis=xaxis) - arr) / dx)
    newfield.data = d.map_overlap(xdrv, depth={'x_u':(1)}, 
                                  boundary={xaxis:'periodic'})
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

    Args:
        field : nD field
        zn: zcoord - needed if changing vertical grid.
        grid = 'p': destination grid

    Returns:
        field on required grid
    @author: Peter Clark
    """

    print("d_by_dy_field_on_w ",grid)
    d = field.data
    y = field.coords['y_p'].data
    dy = y[1] - y[0]
    yaxis = field.get_axis_num('y_p')
    depths = np.zeros(len(field.dims))
    depths[yaxis] = 1
    newfield = field.rename({'y_p':'y_v'})
    ydrv = lambda arr:((np.roll(arr,-1,axis=yaxis) - arr) / dy)
    newfield.data = d.map_overlap(ydrv, depth={'y_v':(1)} ,
                                  boundary={yaxis:'periodic'})
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

    Args:
        field : nD field
        zn: zcoord - needed if changing vertical grid.
        grid = 'p': destination grid

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

    Args:
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

    Args:
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


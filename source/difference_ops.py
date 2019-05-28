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

def interpolate(field, z, zn) :
    """
    Interpolate field from z to zn 
    
    Args:
        field : nD field
        z     : 1D array
        zn    : 1D array
    
    Returns:
        field on zn levels
    @author: Peter Clark
    """ 

    zd = last_dim(z)
    znd = last_dim(zn)   
#    print "zd",zd
#    print "znd",znd
    ss = np.shape(field)[:-1]
    newfield=np.zeros((ss+(len(znd),)))
#    print("Newfield", np.shape(newfield))
    for i in range(0,len(znd)):
        k = np.max([np.searchsorted(zd, znd[i], side='right')-1, 0])
        if (k == (len(zd)-1)) :
            k = k-1
        w = (znd[i] - zd[k])/(zd[k+1] - zd[k])
#        print k, znd[i], zd[k],zd[k+1],w
        newfield[...,i] = w * field[..., k+1] + (1 - w) * field[..., k]
    return newfield

def field_on_w_to_p(field, z, zn) :
    print("w_to_p")
    return interpolate(field, z, zn)

def field_on_p_to_w(field, z, zn) :
    print("p_to_w")
    return interpolate(field, zn, z)

def field_on_u_to_p(field, xaxis=0) :
    """
    Interpolate field from u to p points in C grid 
    
    Args:
        field : nD field
        xaxis=0: Specify axis to interpolate
    
    Returns:
        field on p points
    @author: Peter Clark
    """

    d = field[...]
    newfield=0.5 * (d + np.roll(d,-1,axis=xaxis))
    return newfield

def field_on_p_to_u(field, xaxis=0) :
    """
    Interpolate field from p to u points in C grid 
    
    Args:
        field : nD field
        xaxis=0: Specify axis to interpolate
    
    Returns:
        field on p points
    @author: Peter Clark
    """

    d = field[...]
    newfield=0.5 * (d + np.roll(d,1,axis=xaxis))
    return newfield

def field_on_v_to_p(field, xaxis=0) :
    """
    Interpolate field from v to p points in C grid 
    
    Args:
        field : nD field
        xaxis=0: Specify xaxis 
    
    Returns:
        field on p points
    @author: Peter Clark
    """

    d = field[...]
    newfield=0.5 * (d + np.roll(d,-1,axis=xaxis+1))
    return newfield

def field_on_p_to_v(field, xaxis=0) :
    """
    Interpolate field from p to v points in C grid 
    
    Args:
        field : nD field
        yaxis=1: Specify axis to interpolate
    
    Returns:
        field on p points
    @author: Peter Clark
    """

    d = field[...]
    newfield=0.5 * (d + np.roll(d,1,axis=xaxis+1))
    return newfield

def d_by_dx_field_on_u(field, dx, z, zn, xaxis=0, grid = 'p' ) :
    """
    Differentiate field on u points in x direction then average to req grid 
    
    Args:
        field : nD field
        dx    : grid length
        xaxis=0: Specify axis corresponding to x
        grid = 'p': destination grid
    
    Returns:
        field on required grid
    @author: Peter Clark
    """

    print("d_by_dx_field_on_u ",grid)
    d = field[...]
    newfield = (d - np.roll(d,1,axis=xaxis)) / dx
    # Derivative on p
    if grid == 'u' :
        newfield = field_on_p_to_u(newfield, xaxis=xaxis)
    if grid == 'v' :
        newfield = field_on_p_to_v(newfield, xaxis=xaxis)
    if grid == 'w' :
        newfield = field_on_p_to_w(newfield, z, zn)
        
    return newfield

def d_by_dy_field_on_u(field, dy, z, zn, xaxis=0, grid = 'p' ) :
    """
    Differentiate field on u points in y direction then average to req grid 
    
    Args:
        field : nD field
        dx    : grid length
        xaxis=0: Specify axis corresponding to x
        grid = 'p': destination grid
    
    Returns:
        field on required grid
    @author: Peter Clark
    """

    print("d_by_dy_field_on_u ",grid)
    d = field[...]
    newfield = (np.roll(d,-1,axis=xaxis+1) - d) / dy
    # Derivative on u,v
    if grid == 'p' :
        newfield = field_on_u_to_p(newfield, xaxis=xaxis)
        newfield = field_on_v_to_p(newfield, xaxis=xaxis)
    if grid == 'u' :
        newfield = field_on_v_to_p(newfield, xaxis=xaxis)
    if grid == 'v' :
        newfield = field_on_u_to_p(newfield, xaxis=xaxis)
    if grid == 'w' :
        newfield = field_on_u_to_p(newfield, xaxis=xaxis)
        newfield = field_on_v_to_p(newfield, xaxis=xaxis)
        newfield = field_on_p_to_w(newfield, z, zn)
        
    return newfield

def d_by_dz_field_on_u(field, z, zn, xaxis=0, grid = 'p' ) :
    """
    Differentiate field on u points in z direction then average to req grid 
    
    Args:
        field : nD field
        dx    : grid length
        xaxis=0: Specify axis corresponding to x
        grid = 'p': destination grid
    
    Returns:
        field on required grid
    @author: Peter Clark
    """

    print("d_by_dz_field_on_u ",grid)
#    print(zn)
    d = field[...]
    new = (d[..., 1:] - d[...,:-1])/ (zn[1:] - zn[:-1])
    zt = 0.5 * (zn[1:] + zn[:-1])
    newfield, newz = padright(new, zt, axis=xaxis+2)
    # Derivative on u,w
    if grid == 'p' :
        newfield = field_on_w_to_p(newfield, newz, zn) 
        newfield = field_on_u_to_p(newfield, xaxis=xaxis) 
    if grid == 'u' :
        newfield = field_on_w_to_p(newfield, newz, zn) 
    if grid == 'v' :
        newfield = field_on_w_to_p(newfield, newz, zn) 
        newfield = field_on_u_to_p(newfield, xaxis=xaxis)
        newfield = field_on_p_to_v(newfield, xaxis=xaxis)
    if grid == 'w' :
        newfield = field_on_u_to_p(newfield, xaxis=xaxis)
        
    return newfield

def d_by_dx_field_on_v(field, dy, z, zn, xaxis=0, grid = 'p' ) :
    """
    Differentiate field on u points in y direction then average to req grid 
    
    Args:
        field : nD field
        dx    : grid length
        xaxis=0: Specify axis corresponding to x
        grid = 'p': destination grid
    
    Returns:
        field on required grid
    @author: Peter Clark
    """

    print("d_by_dx_field_on_v ",grid)
    d = field[...]
    newfield = (np.roll(d,-1,axis=xaxis) - d) / dy
    # Derivative on u,v
    if grid == 'p' :
        newfield = field_on_u_to_p(newfield, xaxis=xaxis)
        newfield = field_on_v_to_p(newfield, xaxis=xaxis)
    if grid == 'u' :
        newfield = field_on_v_to_p(newfield, xaxis=xaxis)
    if grid == 'v' :
        newfield = field_on_u_to_p(newfield, xaxis=xaxis)
    if grid == 'w' :
        newfield = field_on_u_to_p(newfield, xaxis=xaxis)
        newfield = field_on_v_to_p(newfield, xaxis=xaxis)
        newfield = field_on_p_to_w(newfield, z, zn)
        
    return newfield

def d_by_dy_field_on_v(field, dy, z, zn, xaxis=0, grid = 'p' ) :
    """
    Differentiate field on u points in x direction then average to req grid 
    
    Args:
        field : nD field
        dy    : grid length
        xaxis=0: Specify axis corresponding to x
        grid = 'p': destination grid
    
    Returns:
        field on required grid
    @author: Peter Clark
    """

    print("d_by_dy_field_on_v ",grid)
    d = field[...]
    newfield = (d - np.roll(d,1,axis=xaxis+1)) / dy
    # Derivative on p
    if grid == 'u' :
        newfield = field_on_p_to_u(newfield, xaxis=xaxis)
    if grid == 'v' :
        newfield = field_on_p_to_v(newfield, xaxis=xaxis)
    if grid == 'w' :
        newfield = field_on_p_to_w(newfield, z, zn)
        
    return newfield

def d_by_dz_field_on_v(field, z, zn, xaxis=0, grid = 'p' ) :
    """
    Differentiate field on u points in z direction then average to req grid 
    
    Args:
        field : nD field
        dx    : grid length
        xaxis=0: Specify axis corresponding to x
        grid = 'p': destination grid
    
    Returns:
        field on required grid
    @author: Peter Clark
    """

    print("d_by_dz_field_on_v ",grid)
#    print(zn)
    d = field[...]
    new = (d[..., 1:] - d[...,:-1])/ (zn[1:] - zn[:-1])
    zt = 0.5 * (zn[1:] + zn[:-1])
    newfield, newz = padright(new, zt, axis=xaxis+2)
    # Derivative on v,w
    if grid == 'p' :
        newfield = field_on_w_to_p(newfield, newz, zn) 
        newfield = field_on_u_to_p(newfield, xaxis=xaxis) 
    if grid == 'u' :
        newfield = field_on_w_to_p(newfield, newz, zn) 
    if grid == 'v' :
        newfield = field_on_w_to_p(newfield, newz, zn) 
        newfield = field_on_u_to_p(newfield, xaxis=xaxis)
        newfield = field_on_p_to_v(newfield, xaxis=xaxis)
    if grid == 'w' :
        newfield = field_on_u_to_p(newfield, xaxis=xaxis)
        
    return newfield

def d_by_dx_field_on_p(field, dx, z, zn, xaxis=0, grid = 'p' ) :
    """
    Differentiate field on u points in x direction then average to req grid 
    
    Args:
        field : nD field
        dy    : grid length
        xaxis=0: Specify axis corresponding to x
        grid = 'p': destination grid
    
    Returns:
        field on required grid
    @author: Peter Clark
    """

    print("d_by_dx_field_on_p ",grid)
    d = field[...]
    newfield = (np.roll(d,-1,axis=xaxis)- d ) / dx
    # Derivative on u
    if grid == 'p' :
        newfield = field_on_u_to_p(newfield, xaxis=xaxis)
    if grid == 'v' :
        newfield = field_on_u_to_p(newfield, xaxis=xaxis)
        newfield = field_on_p_to_v(newfield, xaxis=xaxis)
    if grid == 'w' :
        newfield = field_on_u_to_p(newfield, xaxis=xaxis)
        newfield = field_on_p_to_w(newfield, z, zn)
        
    return newfield

def d_by_dy_field_on_p(field, dy, z, zn, xaxis=0, grid = 'p' ) :
    """
    Differentiate field on u points in y direction then average to req grid 
    
    Args:
        field : nD field
        dy    : grid length
        xaxis=0: Specify axis corresponding to x
        grid = 'p': destination grid
    
    Returns:
        field on required grid
    @author: Peter Clark
    """

    print("d_by_dy_field_on_p ",grid)
    d = field[...]
    newfield = (np.roll(d,-1,axis=xaxis+1) - d) / dy
    # Derivative on v
    if grid == 'p' :
        newfield = field_on_v_to_p(newfield, xaxis=xaxis)
    if grid == 'u' :
        newfield = field_on_v_to_p(newfield, xaxis=xaxis)
        newfield = field_on_p_to_u(newfield, xaxis=xaxis)
    if grid == 'w' :
        newfield = field_on_v_to_p(newfield, xaxis=xaxis)
        newfield = field_on_p_to_w(newfield, z, zn)
        
    return newfield

def d_by_dz_field_on_p(field, z, zn, xaxis=0, grid = 'p' ) :
    """
    Differentiate field on u points in z direction then average to req grid 
    
    Args:
        field : nD field
        dx    : grid length
        xaxis=0: Specify axis corresponding to x
        grid = 'p': destination grid
    
    Returns:
        field on required grid
    @author: Peter Clark
    """

    print("d_by_dz_field_on_p ",grid)
    d = field[...]
    new = (d[..., 1:] - d[...,:-1])/ (zn[1:] - zn[:-1])
    print(zn)
    zt = 0.5 * (zn[1:] + zn[:-1])
    newfield, newz = padright(new, zt, axis=xaxis+2)
    # Derivative on w
    if grid == 'p' :
        newfield = field_on_w_to_p(newfield, newz, zn) 
    if grid == 'u' :
        newfield = field_on_w_to_p(newfield, newz, zn) 
        newfield = field_on_p_to_u(newfield, xaxis=xaxis)
    if grid == 'v' :
        newfield = field_on_w_to_p(newfield, newz, zn) 
        newfield = field_on_p_to_v(newfield, xaxis=xaxis)
        
    return newfield

def d_by_dx_field_on_w(field, dx, z, zn, xaxis=0, grid = 'p' ) :
    """
    Differentiate field on u points in x direction then average to req grid 
    
    Args:
        field : nD field
        dy    : grid length
        xaxis=0: Specify axis corresponding to x
        grid = 'p': destination grid
    
    Returns:
        field on required grid
    @author: Peter Clark
    """

    print("d_by_dx_field_on_w ",grid)
    d = field[...]
    newfield = (np.roll(d,-1,axis=xaxis)- d ) / dx
    # Derivative on u,w
    if grid == 'p' :
        newfield = field_on_w_to_p(newfield, z, zn)
        newfield = field_on_u_to_p(newfield, xaxis=xaxis)
    if grid == 'u' :
        newfield = field_on_w_to_p(newfield, z, zn)
    if grid == 'v' :
        newfield = field_on_w_to_p(newfield, z, zn)
        newfield = field_on_u_to_p(newfield, xaxis=xaxis)
        newfield = field_on_p_to_v(newfield, xaxis=xaxis)
        
    return newfield

def d_by_dy_field_on_w(field, dy, z, zn, xaxis=0, grid = 'p' ) :
    """
    Differentiate field on u points in y direction then average to req grid 
    
    Args:
        field : nD field
        dy    : grid length
        xaxis=0: Specify axis corresponding to x
        grid = 'p': destination grid
    
    Returns:
        field on required grid
    @author: Peter Clark
    """

    print("d_by_dy_field_on_w ",grid)
    d = field[...]
    newfield = (np.roll(d,-1,axis=xaxis+1) - d) / dy
    # Derivative on v,w
    if grid == 'p' :
        newfield = field_on_w_to_p(newfield, z, zn)
        newfield = field_on_v_to_p(newfield, xaxis=xaxis)
    if grid == 'u' :
        newfield = field_on_w_to_p(newfield, z, zn)
        newfield = field_on_v_to_p(newfield, xaxis=xaxis)
        newfield = field_on_p_to_u(newfield, xaxis=xaxis)
    if grid == 'v' :
        newfield = field_on_w_to_p(newfield, z, zn)
        
    return newfield

def d_by_dz_field_on_w(field, z, zn, xaxis=0, grid = 'p' ) :
    """
    Differentiate field on u points in z direction then average to req grid 
    
    Args:
        field : nD field
        dx    : grid length
        xaxis=0: Specify axis corresponding to x
        grid = 'p': destination grid
    
    Returns:
        field on required grid
    @author: Peter Clark
    """

    print("d_by_dz_field_on_w ",grid)
#    print(z)
    d = field[...]
    new = (d[..., 1:] - d[...,0:-1])/ (z[1:] - z[:-1])
    zt = 0.5 * (z[1:] + z[:-1])
    newfield, newz = padleft(new, zt, axis=xaxis+2)
    # Derivative on p
    if grid == 'u' :
        newfield = field_on_p_to_u(newfield, xaxis=xaxis)
    if grid == 'v' :
        newfield = field_on_p_to_v(newfield, xaxis=xaxis)
    if grid == 'w' :
        newfield = field_on_p_to_w(newfield, z, newz) 
        
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
#    print(zt)
    newfield = np.zeros(s)
    newfield[...,1:]=f
    newz = np.zeros(np.size(zt)+1)
    newz[1:] = zt
    newz[0] = 2*zt[0]-zt[1]
#    print(newz)
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
#    print(zt)
    newfield = np.zeros(s)
    newfield[...,:-1] = f
    newz = np.zeros(np.size(zt)+1)
    newz[:-1] = zt
    newz[-1] = 2*zt[-1]-zt[-2]
#    print(newz)
    return newfield, newz


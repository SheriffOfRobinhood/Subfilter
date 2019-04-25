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
L_vap = 2.501E6
Rair = 287.058
Cpair = 1005.0
L_over_cp = L_vap / Cpair
grav = 9.81
kappa = Rair/Cpair
mol_wt_water      =           18.0150
mol_wt_air        =           28.9640
epsilon           = mol_wt_water/mol_wt_air
c_virtual = 1.0/epsilon-1.0

               
def filter_variable_list(source_dataset, ref_dataset, \
                         derived_dataset, twod_filter,\
                         var_list=None, grid='p') :
    """
    Create filtered versions of input variables on required grid, stored in derived_dataset.
    
    Args:
        source_dataset  : NetCDF dataset for input
        ref_dataset     : NetCDF dataset for input containing reference profiles
        derived_dataset : NetCDF dataset for derived data
        twod_filter: 2D filter
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
        vard, vdims, varp  = get_data_on_grid(source_dataset, ref_dataset, v, grid)
        ncvar_r, ncvar_s = filter_field(vard, v, vdims, derived_dataset, \
                                        twod_filter, grid=grid, \
                                        three_d=True, sync=False)
#        plt.plot(np.mean(ncvar_r,axis=(0,1,2)),last_dim(zn))
#        plt.title(v[0]+"_r")
#        plt.show()
#        plt.plot(np.mean(ncvar_s,axis=(0,1,2)),last_dim(zn))
#        plt.title(v[0]+"_s")
#        plt.show()
    derived_dataset.sync()
    return var_list
   
def filter_variable_pair_list(source_dataset, ref_dataset, \
                              derived_dataset, twod_filter,\
                              var_list=None, grid='p') :
    """
    Create filtered versions of pairs input variables on A grid, stored in derived_dataset.
    
    Args:
        source_dataset  : NetCDF dataset for input
        derived_dataset : NetCDF dataset for derived data
        twod_filter: 2D filter
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
        var = source_dataset[v[0]]
        vdims = var.dimensions
        svar =quadratic_subfilter(source_dataset, ref_dataset, \
                                  derived_dataset, \
                                  twod_filter, v[0], v[1], grid=grid)
        
        svar_name = "{}_{}_on{}".format(v[0],v[1],grid)
        
        if twod_filter.attributes['filter_type'] == 'domain' : 
            ncsvar = derived_dataset.createVariable(svar_name,"f8",\
                                     (vdims[0],zvar,))
        else :
            ncsvar = derived_dataset.createVariable(svar_name,"f8",\
                                     (vdims[0],vdims[1],vdims[2],zvar,))
        ncsvar[...] = svar
        print(ncsvar)
#        plt.plot(np.mean(ncvar_r,axis=(0,1,2)),last_dim(zn))
#        plt.title(v[0]+"_r")
#        plt.show()
#        plt.plot(np.mean(ncvar_s,axis=(0,1,2)),last_dim(zn))
#        plt.title(v[0]+"_s")
#        plt.show()
    derived_dataset.sync()
    return var_list   


# Flags are: 'u-grid, v-grid, w-grid'

def get_default_variable_list() :
    """
    Provide default variable list.
       Returns:
           var_list.
		
    @author: Peter Clark
	
    
    """

    var_list = ["u", "v", "w", "th", "th_v", "th_L", \
               "q_vapour", "q_cloud_liquid_mass", "q_total"]
# For testing
#    var_list = ["u","w","th", "th_v", "th_L", "q_total"]
    return var_list

def get_default_variable_pair_list() :
    """
    Provide default variable pair list.
       Returns:
           var_list.
		
    @author: Peter Clark
	
    
    """
    var_list = [
                ["u","u"], \
                ["u","v"], \
                ["u","w"], \
                ["v","v"], \
                ["v","w"], \
                ["w","w"], \
                ["u","th"], \
                ["v","th"], \
                ["w","th"], \
                ["u","th_v"], \
                ["v","th_v"], \
                ["w","th_v"], \
                ["u","th_L"], \
                ["v","th_L"], \
                ["w","th_L"], \
                ["u","q_vapour"], \
                ["v","q_vapour"], \
                ["w","q_vapour"], \
                ["u","q_cloud_liquid_mass"], \
                ["v","q_cloud_liquid_mass"], \
                ["w","q_cloud_liquid_mass"], \
                ["u","q_total"], \
                ["v","q_total"], \
                ["w","q_total"], \
              ]
# For testing
#    var_list = [
#                ["u","w"], \
#                ["w","w"], \
#                ["u","th"], \
#                ["w","th"], \
#                ["w","q_total"], \
#              ]
    return var_list

def convolve(field, twod_filter):
    """
    Convolve field filter using fftconvolve using padding.
	
    Args:
        field : 2D field array
        twod_filter: 2D filter array
    
    Returns:
        field convolved with twod_filter
		
    @author: Peter Clark
	
    """
	
    pad_len = int(np.ceil(len(twod_filter)))
    field = np.pad(field, pad_len, mode='wrap')
    result = fftconvolve(field, twod_filter, mode='same')
    return result[pad_len:-pad_len, pad_len:-pad_len]    

def filtered_field_calc(field, twod_filter, three_d=True ):
    """
    Split field into resolved (field_r) and subfilter (field_s).
    
    Args:
        field : 2D field
        twod_filter: 2D filter
        three_d=True : if input has 3 dimensions, interpret as a single time. False: if input has 3 dimensions, interpret first as time.
    
    Returns:
        [field_r, field_s]
		
    @author: Peter Clark
	
    """    
	
    ndims = len(np.shape(field))
    if twod_filter.attributes['filter_type'] == 'domain' : 
        fshape = np.asarray(field.shape)
#        print('Filtering')
#        print(ndims,three_d)
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
#            fshape = fshape[si]
#            
#        field_r = np.zeros(fshape)
        print(si)
        print(np.shape(field[...]))
        print(axis)
        field_r = np.mean(field[...], axis=axis)       
        print(np.shape(field_r ))
        field_s = field[...] - np.reshape(field_r,si)
    else :    
#        field_r = np.zeros(field.size).reshape(field.shape)
        field_r = np.zeros(field.shape)
        if ndims == 2 :
            field_r[:,:] = convolve(field[:,:], twod_filter.data)
        elif  ndims == 3 :
            if three_d :
                for k in range(field.shape[2]) :
                    field_r[:,:,k] = convolve(field[:,:,k], twod_filter.data)
            else : # Assume third dimension is time
                for k in range(field.shape[0]) :
                    field_r[k,:,:] = convolve(field[k,:,:], twod_filter.data)
        elif (ndims == 4):
            for it in range(field.shape[0]) :
                for k in range(field.shape[3]) :
                    field_r[it,:,:,k] = convolve(field[it,:,:,k], twod_filter.data)
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
    
def setup_derived_data_file(source_file, destdir, twod_filter) :
    """
    Create NetCDF dataset for derived data in destdir.
	
    File name is original file name concatenated with twod_filter.id.
    
    Args:
        source_file     : NetCDF file name.
        destdir         : Directory for derived data.
        twod_filter     : Filter
    
    Returns:
        derived_dataset_name, derived_dataset
		
    @author: Peter Clark
    """    
    derived_dataset_name = os.path.basename(source_file) 
    derived_dataset_name = ('.').join(derived_dataset_name.split('.')[:-1])
    derived_dataset_name = derived_dataset_name + "_" + \
        twod_filter.id + ".nc"
    derived_dataset = Dataset(destdir+derived_dataset_name, "w")
    
    derived_dataset.twod_filter_id = twod_filter.id
    derived_dataset.setncatts(twod_filter.attributes)
    
    source_dataset = Dataset(source_file,"r")
    w = source_dataset["w"]
#    u = source_dataset["u"]

    tvar = w.dimensions[0]    

    times = nc_dimcopy(source_dataset, derived_dataset, tvar)
    z = nc_dimcopy(source_dataset, derived_dataset, "z")
    zn = nc_dimcopy(source_dataset, derived_dataset, "zn")
    
    derived_dataset.createDimension("x",np.shape(w[:])[1])
    x = derived_dataset.createVariable("x","f8",("x",))
    derived_dataset.createDimension("y",np.shape(w[:])[2])
    y = derived_dataset.createVariable("y","f8",("y",))
  
    derived_dataset.sync()
       
    return derived_dataset_name, derived_dataset

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

    var_properties = {"u":[True,False,False],\
                      "v":[False,True,False],\
                      "w":[False,False,True],\
                      "th":[False,False,False],\
                      "q_vapour":[False,False,False],\
                      "q_cloud_liquid_mass":[False,False,False],\
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
        print("Data not in dataset")
        if var_name == 'th_L' :
#            rhoref = ref_dataset.variables['rhon'][-1,...]
            theta, vardim, varp = get_data(source_dataset, ref_dataset, 'th') 
            pref = ref_dataset.variables['prefn'][-1,...]
            piref = (pref[:]/1.0E5)**kappa
            q_cl, vd, vp = get_data(source_dataset, ref_dataset, \
                                    'q_cloud_liquid_mass') 
#            print("Delta theta_L",np.max(L_over_cp * q_cl / piref))
            vard = theta - L_over_cp * q_cl / piref
#            input("Press enter")
        if var_name == 'th_v' :
#            rhoref = ref_dataset.variables['rhon'][-1,...]
            theta, vardim, varp = get_data(source_dataset, ref_dataset, 'th') 
            thref = ref_dataset['thref'][-1,...]
            q_v, vardim, varp = get_data(source_dataset, ref_dataset, \
                                         'q_vapour') 
            q_cl, vd, vp = get_data(source_dataset, ref_dataset, \
                                    'q_cloud_liquid_mass') 
#            print("Delta theta_L",np.max(L_over_cp * q_cl / piref))
            vard = theta + thref * (c_virtual * q_v - q_cl)
        if var_name == 'q_total' :
#            rhoref = ref_dataset.variables['rhon'][-1,...]
            q_v, vardim, varp = get_data(source_dataset, ref_dataset, \
                                         'q_vapour') 
            q_cl, vd, vp = get_data(source_dataset, ref_dataset, \
                                    'q_cloud_liquid_mass') 
            vard = q_v + q_cl
    return vard, vardim, varp

def get_data_on_grid(source_dataset, ref_dataset, var_name, grid='p') :
    """
    Read in 3D data from NetCDF file and, where necessary, interpolate to p grid. 
	
    Assumes first dimension is time.
    
    Args:
        source_dataset  : NetCDF dataset 
        ref_dataset     : NetCDF dataset containing reference profiles.
        var_name        : Name of variable
		rrid='p'        : Destination grid. 'u', 'v', 'w' or 'p'.
    
    Returns:
        variable_dimensions, variable_grid_properties.
		
    @author: Peter Clark
    """    
    var, vdim, vp = get_data(source_dataset, ref_dataset, var_name)
    print(np.shape(var), vdim, vp)
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
            var = field_on_z_to_zn(var,  z, zn)
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
            var = field_on_z_to_zn(var,  z, zn)
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
            var = field_on_z_to_zn(var,  z, zn)
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
    vard = var[...]
    return vard, vdim, vp

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

def filter_field(vard, vname, vdims, derived_dataset, twod_filter, grid='p', \
                 three_d=True, sync=False) :
    """
    Create filtered versions of input variable on required grid, stored in derived_dataset.
       
    Args:
        derived_dataset : NetCDF dataset for derived data.
        twod_filter     : 2D filter.
        vard            : data array.
        vname           : variable name.
        vdims           : array of variable dimension names.
        default provided by get_default_variable_list()
        grid='p'        : Grid - 'u','v','w' or 'p'
		
    Returns:
        ncvar_r, ncvar_s: Resolved and subfilter fields as netcf variables in
                          derived_dataset.
						  
    @author: Peter Clark
	
    """     
    if grid=='w' : zvar = "z" 
    else : zvar = "zn"
    print(vname)
    var_r, var_s = filtered_field_calc(vard, twod_filter, three_d=True )
            
    if twod_filter.attributes['filter_type'] == 'domain' : 
        ncvar_r = derived_dataset.createVariable(vname+"_r"+"_on"+grid,"f8",\
                                 (vdims[0],zvar,))
    else :
        ncvar_r = derived_dataset.createVariable(vname+"_r"+"_on"+grid,"f8",\
                                 (vdims[0],vdims[1],vdims[2],zvar,))
        
    ncvar_r[...] = var_r
    print(ncvar_r)
    ncvar_s = derived_dataset.createVariable(vname+"_s"+"_on"+grid,"f8",\
                                 (vdims[0],vdims[1],vdims[2],zvar,))
    ncvar_s[...] = var_s
    print(ncvar_s)
    
    if sync : derived_dataset.sync()
    
    return ncvar_r, ncvar_s

def filtered_deformation(source_dataset, derived_dataset, twod_filter, \
                         dx, dy, z, zn, xaxis=0, grid='p') :
    """
    Create filtered versions of deformation field.

    Args:
        source_dataset  : NetCDF input dataset 
        derived_dataset : NetCDF dataset for derived data.
        twod_filter     : 2D filter.
        grid='p'        : Grid - 'u','v','w' or 'p'
		
    Returns:
        ncvar_r, ncvar_s: Resolved and subfilter fields as netcf variables in derived_dataset.
						  
    @author: Peter Clark

    """ 
    
#:math:`\frac{\partial u_i}{\partial{x_j}`
    u = source_dataset["u"]
    vdims = u.dimensions
    vname = 'deformation'
        
    d = deformation(source_dataset, dx, dy, z, zn, xaxis=xaxis, \
                           grid=grid)  
    d_ij_r = list()
    d_ij_s = list()
    for i in range(3) :
        d_j_r = list()
        d_j_s = list()
        for j in range(3) :
            vn = "{}_{:1d}_{:1d}".format(vname,i,j)
            def_r, def_s = filter_field(d[i][j], vn, vdims, derived_dataset, \
                                twod_filter, grid='p', \
                                three_d=True, sync=False)
            d_j_r.append(def_r)
            d_j_s.append(def_s)
        d_ij_r.append(d_j_r)
        d_ij_s.append(d_j_s)
    derived_dataset.sync()
    
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
            S["{:1d}_{:1d}".format(i,j)] = (S_ij)
    return S, mod_S

def vorticity(d) :

    V_0 = d[2][1][...]-d[1][2][...] 
    V_1 = d[0][3][...]-d[3][0][...] 
    V_2 = d[1][0][...]-d[0][1][...] 
    V = [V_0, V_1, V_2]
    
    return V 

def quadratic_subfilter(source_dataset,  ref_dataset, \
                        derived_dataset, twod_filter,
                        v1_name, v2_name, grid='p') :
    """
    Create filtered versions of pair of input variables on required grid, stored in derived_dataset.
    
    Args:
        source_dataset  : NetCDF dataset for input
        derived_dataset : NetCDF dataset for derived data
        twod_filter: 2D filter
        v1_name, v2_name: Variable names.
    
    Returns:
        s(var1,var2) data array.
		
    @author: Peter Clark
	
    """
    vard1, vd1, vp1 = get_data_on_grid(source_dataset,  ref_dataset, \
                                       v1_name, grid=grid)
    print("Reading ", v1_name+"_r"+"_on"+grid)
    var1_r = derived_dataset[v1_name+"_r"+"_on"+grid]
    
    print(np.shape(var1_r))
        
    vard2, vd2, vp2 = get_data_on_grid(source_dataset,  ref_dataset, \
                                       v2_name, grid=grid)
    print("Reading ", v2_name+"_r"+"_on"+grid)
    var2_r = derived_dataset[v2_name+"_r"+"_on"+grid]
    print(np.shape(var2_r))
    
    var1var2 = vard1 * vard2
    var1var2_r, var1var2_s = filtered_field_calc(var1var2, twod_filter, \
                                                 three_d=True )
    
    var1var2 = var1var2_r - var1_r[...] * var2_r[...]
    
    return var1var2
    

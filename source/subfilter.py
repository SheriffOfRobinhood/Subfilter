# -*- coding: utf-8 -*-
"""
Created on Tue Oct 23 11:07:05 2018

@author: paclk
"""
import os 
import netCDF4
from netCDF4 import Dataset

import numpy as np
from scipy.signal import fftconvolve
#from sys import float_info
#import filters 
import time
#import matplotlib
#import matplotlib.pyplot as plt
var_properties = {"u":[True,False,False],\
                  "v":[False,True,False],\
                  "w":[False,False,True],\
                  "th":[False,False,False],\
                  "q_vapour":[False,False,False],\
                  "q_cloud_liquid_mass":[False,False,False],\
                 }

def convolve(field, twod_filter_array):
    pad_len = int(np.ceil(len(twod_filter_array)))
    field = np.pad(field, pad_len, mode='wrap')
    result = fftconvolve(field, twod_filter_array, mode='same')
    return result[pad_len:-pad_len, pad_len:-pad_len]    

def filtered_field_calc(field, twod_filter, three_d=True ):
    field_r = np.zeros(field.size).reshape(field.shape)
    ndims = len(np.shape(field))
    if ndims == 2 :
        field_r[:,:] = convolve(field[:,:], twod_filter.data)
    elif  (ndims == 3) :
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

def last_dim(z) :
    zd = z[...]
    while len(np.shape(zd))>1 :
        zd = zd[0,...]
    return zd
    
def field_on_z_to_zn(field, z, zn) :
    zd = last_dim(z)
    znd = last_dim(zn)   
#    print "zd",zd
#    print "znd",znd
    ss = np.shape(field)[:-1]
    newfield=np.zeros((ss+(len(znd),)))
    print "Newfield", np.shape(newfield)
    for i in range(0,len(znd)):
#    for i in range(0,5):
        if (zd[0] > znd[i]) :
            k = 0
        else :
            k = np.where(zd[:] <= znd[i]  )
            if len(k[0]) > 0 :
                k = k[0][-1]
            else :
                print "Interpolation error"
            if (k == (len(zd)-1)) :
                k = k-1
        w = (znd[i] - zd[k])/(zd[k+1] - zd[k])
#        print k, znd[i], zd[k],zd[k+1],w
        newfield[...,i] = w * field[..., k+1] + (1 - w) * field[..., k]
    return newfield

def field_on_u_to_p(field, xaxis=0) :
    d = field[...]
    newfield=0.5 * (d + np.roll(d,-1,axis=xaxis))
    return newfield

def field_on_v_to_p(field, yaxis=1) :
    d = field[...]
    newfield=0.5 * (d + np.roll(d,-1,axis=yaxis))
    return newfield

def nc_dimcopy(source_dataset, derived_dataset, dimname) :
    v = source_dataset.variables[dimname]
    derived_dataset.createDimension(dimname,np.shape(v[:])[-1])
    dv = derived_dataset.createVariable(dimname,"f8",(dimname,))
    dv[:] = last_dim(v[:])
    return dv
    
def setup_derived_data_file(source_file, destdir, twod_filter) :
    derived_dataset_name = os.path.basename(source_file) 
    derived_dataset_name = ('.').join(derived_dataset_name.split('.')[:-1])
    derived_dataset_name = derived_dataset_name + "_" + \
        twod_filter.id + ".nc"
    derived_dataset = Dataset(destdir+derived_dataset_name, "w")
    
    derived_dataset.filter_id = twod_filter.id
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

def get_data_on_A_grid(source_dataset, var_name):
    z = source_dataset["z"]
    zn = source_dataset["zn"]
    print var_name
    var = source_dataset[var_name]
    vp = var_properties[var_name]
    if vp[0] :
        print "Mapping %s from u grid to p grid."%var
        var = field_on_u_to_p(var, xaxis=1)
    if vp[1] :
        print "Mapping %s from v grid to p grid."%var
        var = field_on_v_to_p(var, yaxis=2)
    if vp[2] :
        print "Mapping %s from w grid to p grid."%var
        var = field_on_z_to_zn(var,  z, zn)    
    vard = var[...]
    return vard

def filter_variable_list(source_dataset, derived_dataset, twod_filter,\
                         var_list=None) :
    if (var_list==None):
        var_list = ["u","v","w","th","q_vapour","q_cloud_liquid_mass"]
        print "Default list:\n",var_list
    z = source_dataset["z"]
    zn = source_dataset["zn"]
    for v in var_list:
        var = source_dataset[v]
        vdims = var.dimensions
#        print vdims
        vard = get_data_on_A_grid(source_dataset, v)
        var_r, var_s = filtered_field_calc(vard, twod_filter, three_d=True )
        
        ncvar_r = derived_dataset.createVariable(v+"_r","f8",\
                                     (vdims[0],vdims[1],vdims[2],"zn",))
        ncvar_r[...] = var_r
        print ncvar_r
        ncvar_s = derived_dataset.createVariable(v+"_s","f8",\
                                     (vdims[0],vdims[1],vdims[2],"zn",))
        ncvar_s[...] = var_s
        print ncvar_s
#        plt.plot(np.mean(ncvar_r,axis=(0,1,2)),last_dim(zn))
#        plt.title(v[0]+"_r")
#        plt.show()
#        plt.plot(np.mean(ncvar_s,axis=(0,1,2)),last_dim(zn))
#        plt.title(v[0]+"_s")
#        plt.show()
    derived_dataset.sync()
    return var_list

def quadratic_subfilter(source_dataset, derived_dataset, twod_filter,
                        v1_name, v2_name) :
    z = source_dataset["z"]
    zn = source_dataset["zn"]  
    
    var1 = source_dataset[v1_name]
    vard1 = get_data_on_A_grid(source_dataset, v1_name)
    var1_r = derived_dataset[v1_name+"_r"]
    
    var2 = source_dataset[v2_name]
    vard2 = get_data_on_A_grid(source_dataset, v2_name)
    var2_r = derived_dataset[v2_name+"_r"]
    
    var1var2 = vard1[...] * vard2[...]
    var1var2_r, var1var2_s = filtered_field_calc(var1var2, twod_filter, \
                                                 three_d=True )
    
    var1var2 = var1var2_r - var1_r[...] * var2_r[...]
    
    return var1var2
   
def filter_variable_pair_list(source_dataset, derived_dataset, twod_filter,\
                         var_list=None) :
    if (var_list==None):
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
                    ["u","q_vapour"], \
                    ["v","q_vapour"], \
                    ["w","q_vapour"], \
                    ["u","q_cloud_liquid_mass"], \
                    ["v","q_cloud_liquid_mass"], \
                    ["w","q_cloud_liquid_mass"], \
                  ]
#                    ,"v","w","th","q_vapour","q_cloud_liquid_mass"]
        print "Default list:\n",var_list
    z = source_dataset["z"]
    zn = source_dataset["zn"]
    for v in var_list:
        print "Calculating s(%s,%s)"%(v[0],v[1])
        var = source_dataset[v[0]]
        vdims = var.dimensions
        svar =quadratic_subfilter(source_dataset, derived_dataset, \
                        twod_filter, v[0], v[1])
        
        svar_name = "%s_%s"%(v[0],v[1])
        
        ncsvar = derived_dataset.createVariable(svar_name,"f8",\
                                     (vdims[0],vdims[1],vdims[2],"zn",))
        ncsvar[...] = svar
        print ncsvar
#        plt.plot(np.mean(ncvar_r,axis=(0,1,2)),last_dim(zn))
#        plt.title(v[0]+"_r")
#        plt.show()
#        plt.plot(np.mean(ncvar_s,axis=(0,1,2)),last_dim(zn))
#        plt.title(v[0]+"_s")
#        plt.show()
    derived_dataset.sync()
    return var_list   
    

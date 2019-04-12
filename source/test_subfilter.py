# -*- coding: utf-8 -*-
"""
Created on Tue Oct 23 11:27:25 2018

@author: Peter Clark
"""
import os 
import netCDF4
from netCDF4 import Dataset
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

import subfilter as sf
import filters as filt

dir = 'C:/Users/paclk/OneDrive - University of Reading/python/subfilter/test_data/BOMEX/'
file = 'diagnostics_ts_18000.0.nc'

#w = dataset.variables['w']
#w_tvar = w.dimensions[0]
#w_time = dataset.variables[w_tvar]

plot_dir = '' 

plot_type = '.eps'
data_dir = '' # Directory containing data
figshow = True

def plot_w(name, twod_filter, w_r, w_s,  ilev, ix1, ix2, iy1, iy2, iy, j):
    
    meanfield= np.mean(w_r,axis=(0,1),keepdims=True)
    pltdat = (w_r-meanfield)

    lev1 = np.arange(-10,10.1,0.1)
    lev2 = np.arange(-10,10.1,0.1)
    plt.clf

    fig1, axa = plt.subplots(3,2,figsize=(10,12))
    
#    plt.subplot(3, 2, 1)
    Cs1 = axa[0,0].contourf(np.transpose(pltdat[ix1:ix2, iy1:iy2, ilev]),levels=lev1)
    axa[0,0].set_title(r'%s$^r$ pert level %03d'%(name,ilev))
    axa[0,0].set_xlabel('x')

# Make a colorbar for the ContourSet returned by the contourf call.
    cbar1 = fig1.colorbar(Cs1,ax=axa[0,0])
    cbar1.ax.set_ylabel(name)
# Add the contour line levels to the colorbar
#  cbar.add_lines(CS2)

#    plt.subplot(3, 2, 2)
    Cs2 = axa[0,1].contourf(np.transpose(w_s[ix1:ix2, iy1:iy2, ilev]),levels=lev2)
    axa[0,1].set_xlabel('x')
#    Cs2.set_xlabel('x')
    axa[0,1].set_title(r'%s$^s$ level %03d'%(name,ilev))
    cbar2 = fig1.colorbar(Cs2,ax=axa[0,1])
    cbar2.ax.set_ylabel(name)
#
#    plt.subplot(3, 2, 3)
    Cs3 = axa[1,0].contourf(np.transpose(pltdat[ix1:ix2,iy,:]),levels=lev1)
#
##    plt.ylim([0,5000])
    axa[1,0].set_title(r'%s$^r$ pert at iy %03d'%(name,iy))
## Make a colorbar for the ContourSet returned by the contourf call.
#    cbar3 = plt.colorbar(Cs3)
    cbar3 = fig1.colorbar(Cs3,ax=axa[1,0])
    cbar3.ax.set_ylabel(name)
#
#    plt.subplot(3, 2, 4)
    Cs4 = axa[1,1].contourf(np.transpose(w_s[ix1:ix2,iy,:]),levels=lev2)
#
##    plt.ylim([0,5000])
    axa[1,1].set_title(r'%s$^s$ at iy %03d'%(name,iy))
## Make a colorbar for the ContourSet returned by the contourf call.
#    cbar4 = plt.colorbar(Cs4)
    cbar4 = fig1.colorbar(Cs4,ax=axa[1,1])
    cbar4.ax.set_ylabel(name)
#    
    x=(np.arange(ix1,ix2)-0.5*(ix1+ix2-1))*0.1
#    plt.subplot(3, 2, 5)
    ax1 = axa[2,0].plot(x,pltdat[ix1:ix2,iy,ilev])
#    
#    plt.subplot(3, 2, 6)
    ax2 = axa[2,1].plot(x,w_s[ix1:ix2,iy,ilev])
#    
    plt.savefig('w_lev_'+'%03d'%ilev+'_x_z'+'%03d'%iy+'_%02d'%j+'_'+\
                twod_filter.id+'.png')
#
    plt.tight_layout()
    plt.show()
    plt.close()

    return

def test_filter(dataset, twod_filter, contents=0):
#   Work out how many times in data
    ix = 32
#   Set up list of times
    iy = 40
    ilev = 10
    w = dataset["w"]
    u = dataset["u"]
#    x = dataset["x"]
#    print "x",np.shape(x)
#    y = dataset["y"]
    z = dataset["z"]
    print "z",np.shape(z)
#    print z[:]
    zn = dataset["zn"]
    print "zn",np.shape(zn)
#    print zn[:]
    print w
    
    wt =  sf.field_on_z_to_zn(w, z, zn)
    ut = field_on_u_to_p(u, xaxis=1)
    
    uw = ut*wt
    plt.figure(1,figsize=(10,8))
    plt.clf
#    plt.subplot(2, 1, 1)
#    lev1 = np.arange(-1,1.1,0.1)
    plt.plot(np.mean(uw[3,...],axis=(0,1)),z[0,:],label='w')
#    Cs1 = plt.contourf(np.transpose(w[0,:,:, ile,v]),levels=lev1)
#    plt.subplot(2, 1, 2)
#    Cs2 = plt.contourf(np.transpose(wt[0,:,:, ilev]),levels=lev1)
#    plt.plot(np.mean(uw[3,...],axis=(0,1)),zn[:],label='wt')
    plt.legend()
   
    plt.show()
    
    n_times = np.shape(w)[0]
 
#    print "w shape:", np.shape(w)
#    w = np.zeros(np.shape(w))

#    w[:,ix,iy,:]=1E-3
   
    for j in range(n_times):

        print "Filtering time %d uw[%d] shape:"%(j,j), np.shape(uw[j])

        uw_f = sf.filtered_field_calc(uw[j], twod_filter)

        uw_r = uw_f[0]
        uw_s = uw_f[1]

        ix1 = 00
        ix2 = np.shape(uw_r)[1]-1
        iy1 = 00
        iy2 = np.shape(uw_r)[2]-1
        
#        print w_r
#        print w_s

        plot_w('uw',twod_filter, uw_r, uw_s,  ilev, ix1, ix2, iy1, iy2, iy, j)    

    return uw_f

def main():
    '''
    Top level code, a bit of a mess.
    '''
#   Non-global variables that are set once
#    sigma_list = [0.2,0.25,0.3,0.4,0.8,1.0,1.5,2.0] # Sigma used in Gaussian filter function
#    sigma_list = [2.0] # Sigma used in Gaussian#
    sigma_list = [0.5,0.2]
    width = -1
    filter_name = 'gaussian'
#    width = 20
#    filter_name = 'running_mean'   
    
    dataset = Dataset(dir+file, 'r') # Dataset is the class behavior to open the file
                                 # and create an instance of the ncCDF4 class

#
# Just looking at anyhb at present
    
    filter_list = list([])

    for i,sigma in enumerate(sigma_list):
        filter_id = 'filter_%02d'%i
        twod_filter = filt.filter_2d(filter_id,\
                                   filter_name, \
                                   sigma=sigma, width=width, \
                                   delta_x=0.1)
        
        print twod_filter
        filter_list.append(twod_filter)
        
    print filter_list
    for twod_filter in filter_list:
        
        derived_dataset_name, derived_data = sf.setup_derived_data_file(\
                                            dir+file, dir, twod_filter)
        print derived_data.variables

        sf.filter_variable_list(dataset, derived_data, twod_filter, var_list=None)         
        times = derived_data['time_series_50_100.0']
        print times
        print times[:]
        
        z = derived_data["z"]
        print "z",np.shape(z[:])   
        print z[:]
        zn = derived_data["zn"]
        print "zn",np.shape(zn[:])   
        print zn[:]
        
        u_r = derived_data["u_r"]
        print u_r
        derived_data.close()
#        os.remove(derived_dataset_name)
#        print twod_filter
        fig1 = plt.figure(1)
    #    print twod_filter[:,:]
        plt.contourf(twod_filter.data,20) 
        
        fig2 = plt.figure(2)
        plt.plot(twod_filter.data[np.shape(twod_filter.data)[0]/2,:])
        plt.show()
#        test_filter(dataset, twod_filter)
        
#    filter_dataset.close()
    
if __name__ == "__main__":
    main()
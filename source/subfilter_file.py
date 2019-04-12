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
#var_tvar = w.dimensions[0]
#var_time = dataset.variables[var_tvar]

plot_dir = 'C:/Users/paclk/OneDrive - University of Reading/python/subfilter/test_data/BOMEX/plots/' 

plot_type = '.png'
data_dir = '' # Directory containing data
figshow = True

def plot_field(var_name, derived_data, twod_filter, ilev, iy):
    
    var_r = derived_data[var_name+"_r"]
    var_s = derived_data[var_name+"_s"]
    
    for it in range(var_r.shape[0]):
        meanfield= np.mean(var_r[it,...],axis=(0,1),keepdims=True)
        pltdat = (var_r[it,...]-meanfield)
    
        lev1 = np.arange(-10,10.1,0.1)
        lev2 = np.arange(-10,10.1,0.1)
        
        nlevels = 40
        plt.clf
    
        fig1, axa = plt.subplots(3,2,figsize=(10,12))
        
    #    plt.subplot(3, 2, 1)
        Cs1 = axa[0,0].contourf(np.transpose(pltdat[:, :, ilev]),\
                 nlevels)
        axa[0,0].set_title(r'%s$^r$ pert level %03d'%(var_name,ilev))
        axa[0,0].set_xlabel('x')
    
    # Make a colorbar for the ContourSet returned by the contourf call.
        cbar1 = fig1.colorbar(Cs1,ax=axa[0,0])
        cbar1.ax.set_ylabel(var_name)
    # Add the contour line levels to the colorbar
    #  cbar.add_lines(CS2)
    
    #    plt.subplot(3, 2, 2)
        Cs2 = axa[0,1].contourf(np.transpose(var_s[it, :, :, ilev]),\
                 nlevels)
        axa[0,1].set_xlabel('x')
        axa[0,1].set_title(r'%s$^s$ level %03d'%(var_name,ilev))
        cbar2 = fig1.colorbar(Cs2,ax=axa[0,1])
        cbar2.ax.set_ylabel(var_name)
    #
    #    plt.subplot(3, 2, 3)
        Cs3 = axa[1,0].contourf(np.transpose(pltdat[:,iy,:]),nlevels)
    #
    ##    plt.ylim([0,5000])
        axa[1,0].set_title(r'%s$^r$ pert at iy %03d'%(var_name,iy))
    ## Make a colorbar for the ContourSet returned by the contourf call.
        cbar3 = fig1.colorbar(Cs3,ax=axa[1,0])
        cbar3.ax.set_ylabel(var_name)
    #
    #    plt.subplot(3, 2, 4)
        Cs4 = axa[1,1].contourf(np.transpose(var_s[it, :, iy,:]),nlevels)
    #
    ##    plt.ylim([0,5000])
        axa[1,1].set_title(r'%s$^s$ at iy %03d'%(var_name,iy))
    ## Make a colorbar for the ContourSet returned by the contourf call.
        cbar4 = fig1.colorbar(Cs4,ax=axa[1,1])
        cbar4.ax.set_ylabel(var_name)
    #    
        x=(np.arange(0,var_r.shape[2])-0.5*var_r.shape[2])*0.1
    #    plt.subplot(3, 2, 5)
        ax1 = axa[2,0].plot(x,pltdat[:,iy,ilev])
    #    
    #    plt.subplot(3, 2, 6)
        ax2 = axa[2,1].plot(x,var_s[it,:,iy,ilev])
    #    
        plt.tight_layout()
        
        plt.savefig(plot_dir+var_name+'_lev_'+'%03d'%ilev+'_x_z'+'%03d'%iy+'_'+\
                    twod_filter.id+'_%02d'%it+plot_type)
        plt.close()
        
    #
#        plt.show()
    plt.close()

    return

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
    ilev = 20
    iy = 40 
#
# Just looking at anyhb at present
    
    filter_list = list([])

    for i,sigma in enumerate(sigma_list):
        filter_id = 'filter_%02d'%i
        twod_filter = filt.filter_2d(filter_id,\
                                   filter_name, \
                                   sigma=sigma, width=width, \
                                   delta_x=0.1)
        
        print(twod_filter)
        filter_list.append(twod_filter)
        
    print(filter_list)
    for twod_filter in filter_list:
        
        derived_dataset_name, derived_data = sf.setup_derived_data_file(\
                                            dir+file, dir, twod_filter)
        print(derived_data.variables)

        field_list =sf.filter_variable_list(dataset, derived_data,\
                                            twod_filter, var_list=None)    

        quad_field_list =sf.filter_variable_pair_list(dataset, derived_data, \
                                            twod_filter, var_list=None)        
        times = derived_data['time_series_50_100.0']
        print(times)
        print(times[:])
        
        z = derived_data["z"]
        print("z",np.shape(z[:]))   
        print(z[:])
        zn = derived_data["zn"]
        print("zn",np.shape(zn[:]))   
        print(zn[:])
#        uv = sf.quadratic_subfilter(dataset, derived_data, twod_filter,\
#                        "u", "v")
#        uvm = np.mean(uv,axis=(0,1,2))
#        plt.plot(np.mean(uv,axis=(0,1,2)),zn[:])
#        plt.plot(uvm,zn[:])
#        plt.show()
        
        u_r = derived_data["u_r"]
        print(u_r)
#        os.remove(derived_dataset_name)
#        print twod_filter
        fig1 = plt.figure(1)
    #    print twod_filter[:,:]
        plt.contourf(twod_filter.data,20) 
        
        fig2 = plt.figure(2)
        plt.plot(twod_filter.data[np.shape(twod_filter.data)[0]/2,:])
        plt.show()
        
        for field in field_list:
            print("Plotting {}".format(field))
            plot_field(field, derived_data, twod_filter, ilev, iy)            
#        test_filter(dataset, twod_filter)
        derived_data.close()
        
#    filter_dataset.close()
    
if __name__ == "__main__":
    main()
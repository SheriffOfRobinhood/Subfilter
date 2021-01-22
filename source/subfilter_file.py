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
import difference_ops as do

options = {
#        'FFT_type': 'FFTconvolve',
#        'FFT_type': 'FFT',
        'FFT_type': 'RFFT',
        'save_all': 'Yes',
          }


dir = 'C:/Users/paclk/OneDrive - University of Reading/Git/python/Subfilter/test_data/BOMEX/'
odir = 'C:/Users/paclk/OneDrive - University of Reading/Git/python/Subfilter/test_data/BOMEX/'
odir = odir + 'test_op_' + options['FFT_type']+'/'

os.makedirs(odir, exist_ok = True)

file = 'diagnostics_ts_18000.0.nc'
ref_file = 'diagnostics_ts_18000.0.nc'

#w = dataset.variables['w']
#var_tvar = w.dimensions[0]
#var_time = dataset.variables[var_tvar]

plot_dir = odir + 'plots/'
os.makedirs(plot_dir, exist_ok = True)

plot_type = '.png'
data_dir = '' # Directory containing data
figshow = True

def plot_field(var_name, filtered_data, twod_filter, ilev, iy, grid='p'):

    var_r = filtered_data[f"{var_name}_on_{grid}_r"]
    var_s = filtered_data[f"{var_name}_on_{grid}_s"]


    for it in range(var_r.shape[0]):
        if twod_filter.attributes['filter_type']=='domain' :
            zcoord = sf.last_dim(filtered_data[var_r.dimensions[1]])
            pltdat = (var_r[it,:])

            fig1, axa = plt.subplots(1,1,figsize=(5,5))

    #    plt.subplot(3, 2, 1)
            Cs1 = axa.plot(pltdat, zcoord)
            axa.set_xlabel(r'%s$^r$'%(var_name))
            axa.set_ylabel('z')

            plt.tight_layout()

            plt.savefig(plot_dir+var_name+'_prof_'+\
                    twod_filter.id+'_%02d'%it+plot_type)
            plt.close()
        else :
            zcoord = sf.last_dim(filtered_data[var_r.dimensions[3]])
            meanfield= np.mean(var_r[it,...],axis=(0,1),keepdims=True)
            pltdat = (var_r[it,...]-meanfield)

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
    #    plt.show()
    #plt.close()

    return

def plot_quad_field(var_name, filtered_data, twod_filter, ilev, iy, grid='p'):

    v1 = var_name[0]
    v2 = var_name[1]

    v1_r = filtered_data[f"{v1}_on_{grid}_r"]
    v2_r = filtered_data[f"{v2}_on_{grid}_r"]

    print(v1,v2)
    s_v1v2 = filtered_data[f"{v1}_{v2}_on_{grid}"]
    print(s_v1v2)

    for it in range(s_v1v2.shape[0]):

        if twod_filter.attributes['filter_type']=='domain' :
            pltdat = (s_v1v2[it,:])
            zcoord = sf.last_dim(filtered_data[v1_r.dimensions[1]])

            fig1, axa = plt.subplots(1,1,figsize=(5,5))

    #    plt.subplot(3, 2, 1)
            Cs1 = axa.plot(pltdat, zcoord)
            axa.set_xlabel('s({},{})'.format(v1,v2))
            axa.set_ylabel('z')

            plt.tight_layout()

            plt.savefig(plot_dir+var_name[0]+'_'+var_name[1]+'_prof_'+\
                    twod_filter.id+'_%02d'%it+plot_type)
            plt.close()
        else :
            zcoord = sf.last_dim(filtered_data[v1_r.dimensions[3]])
            var_r = (v1_r[it,...] - np.mean(v1_r[it,...], axis=(0,1))) * \
                    (v2_r[it,...] - np.mean(v2_r[it,...], axis=(0,1)))

            meanfield= np.mean(var_r[...],axis=(0,1),keepdims=True)
            pltdat = (var_r[...]-meanfield)

#            lev1 = np.arange(-10,10.1,0.1)
#            lev2 = np.arange(-10,10.1,0.1)

            nlevels = 40
            plt.clf

            fig1, axa = plt.subplots(3,2,figsize=(10,12))

        #    plt.subplot(3, 2, 1)
            Cs1 = axa[0,0].contourf(np.transpose(pltdat[:, :, ilev]),\
                     nlevels)
            axa[0,0].set_title(r'{}$^r${}$^r$ pert level {:03d}'.format(v1, v2,ilev))
            axa[0,0].set_xlabel('x')

        # Make a colorbar for the ContourSet returned by the contourf call.
            cbar1 = fig1.colorbar(Cs1,ax=axa[0,0])
            cbar1.ax.set_ylabel(var_name)
        # Add the contour line levels to the colorbar
        #  cbar.add_lines(CS2)

        #    plt.subplot(3, 2, 2)
            Cs2 = axa[0,1].contourf(np.transpose(s_v1v2[it, :, :, ilev]),\
                     nlevels)
            axa[0,1].set_xlabel('x')
            axa[0,1].set_title('s({},{}) level {:03d}'.format(v1,v2,ilev))
            cbar2 = fig1.colorbar(Cs2,ax=axa[0,1])
            cbar2.ax.set_ylabel(var_name)
        #
        #    plt.subplot(3, 2, 3)
            Cs3 = axa[1,0].contourf(np.transpose(pltdat[:,iy,:]),nlevels)
        #
        ##    plt.ylim([0,5000])
            axa[1,0].set_title(r'{}$^r${}$^r$ at iy={:03d}'.format(v1, v2,iy))
        ## Make a colorbar for the ContourSet returned by the contourf call.
            cbar3 = fig1.colorbar(Cs3,ax=axa[1,0])
            cbar3.ax.set_ylabel(var_name)
        #
        #    plt.subplot(3, 2, 4)
            Cs4 = axa[1,1].contourf(np.transpose(s_v1v2[it, :, iy,:]),nlevels)
        #
        ##    plt.ylim([0,5000])
            axa[1,1].set_title('s({},{}) at iy={:03d}'.format(v1,v2,iy))
        ## Make a colorbar for the ContourSet returned by the contourf call.
            cbar4 = fig1.colorbar(Cs4,ax=axa[1,1])
            cbar4.ax.set_ylabel(var_name)
        #
            x=(np.arange(0,var_r.shape[1])-0.5*var_r.shape[1])*0.1
        #    plt.subplot(3, 2, 5)
            ax1 = axa[2,0].plot(x,pltdat[:,iy,ilev])
        #
        #    plt.subplot(3, 2, 6)
            ax2 = axa[2,1].plot(x,s_v1v2[it,:,iy,ilev])
        #
            plt.tight_layout()

            plt.savefig(plot_dir+var_name[0]+'_'+var_name[1]+'_lev_'+'%03d'%ilev+'_x_z'+'%03d'%iy+'_'+\
                        twod_filter.id+'_%02d'%it+plot_type)
            plt.close()

    #
    #    plt.show()
    #plt.close()

    return

def plot_shear(var_r, var_s, zcoord,  twod_filter, ilev, iy, no_trace = True):
    var_name = "mod_S"
    if no_trace : var_name = var_name+'n'

    for it in range(var_r.shape[0]):
#        meanfield= np.mean(var_r[it,...],axis=(0,1),keepdims=True)
#        pltdat = (var_r[it,...]-meanfield)
        if twod_filter.attributes['filter_type']=='domain' :
            pltdat = (var_r[it,:])

            fig1, axa = plt.subplots(1,1,figsize=(5,5))

    #    plt.subplot(3, 2, 1)
            Cs1 = axa.plot(pltdat[1:], zcoord[1:])
            axa.set_xlabel(var_name)
            axa.set_ylabel('z')

            plt.tight_layout()

            plt.savefig(plot_dir+var_name+'_prof_'+\
                    twod_filter.id+'_%02d'%it+plot_type)
            plt.close()
        else :
            pltdat = var_r[it,...]

            nlevels = 40
            plt.clf

            fig1, axa = plt.subplots(3,2,figsize=(10,12))

        #    plt.subplot(3, 2, 1)
            Cs1 = axa[0,0].contourf(np.transpose(pltdat[:, :, ilev]),\
                     nlevels)
            axa[0,0].set_title(r'%s$^r$ level %03d'%(var_name,ilev))
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
            axa[1,0].set_title(r'%s$^r$ at iy %03d'%(var_name,iy))
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
    #    plt.show()
#    plt.close()

    return

def main():
    '''
    Top level code, a bit of a mess.
    '''
#   Non-global variables that are set once
    sigma_list = [500.0, 220.0]
    width = -1
    dx = 100.0
    dy = 100.0
    filter_name = 'gaussian'
#    width = 20
#    filter_name = 'running_mean'
#    filter_name = 'wave_cutoff'

    dataset = Dataset(dir+file, 'r') # Dataset is the class behavior to open the file
                                 # and create an instance of the ncCDF4 class
    ref_dataset = Dataset(dir+ref_file, 'r')
    ilev = 15
    iy = 40

    opgrid = 'w'
    fname = 'test_plot'

    derived_dataset_name, derived_data, exists = \
        sf.setup_derived_data_file( dir+file, odir, dir+ref_file, fname,
                                   options, override=True)
    print("Variables in derived dataset.")
    print(derived_data.variables)

    filter_list = list([])

    for i,sigma in enumerate(sigma_list):
        if filter_name == 'gaussian':
            filter_id = 'filter_ga{:02d}'.format(i)
            twod_filter = filt.filter_2d(filter_id,
                                       filter_name,
                                       sigma=sigma, width=width,
                                       delta_x=dx)
        elif filter_name == 'wave_cutoff':
            filter_id = 'filter_wc{:02d}'.format(i)
            twod_filter = filt.filter_2d(filter_id,
                                       filter_name, wavenumber=np.pi/(2*sigma),
                                       width=width,
                                       delta_x=dx)
        elif filter_name == 'running_mean':
            filter_id = 'filter_rm{:02d}'.format(i)
            width = int(np.round( sigma/dx * np.pi * 2.0 / 3.0)+1)
            twod_filter = filt.filter_2d(filter_id,
                                       filter_name,
                                       width=width,
                                       delta_x=dx)

        print(twod_filter)
        filter_list.append(twod_filter)

# Add whole domain filter
    filter_name = 'domain'
    filter_id = 'filter_do{:02d}'.format(len(filter_list))
    twod_filter = filt.filter_2d(filter_id, filter_name, delta_x=dx)
    filter_list.append(twod_filter)

    print(filter_list)

    # z = do.last_dim(dataset["z"])
    zn = do.last_dim(dataset["zn"])

    for twod_filter in filter_list:

        print(twod_filter)

        filtered_dataset_name, filtered_data, exists = \
            sf.setup_filtered_data_file( dir+file, odir, dir+ref_file, fname,
                                       options, twod_filter, override=True)
        print("Variables in filtered dataset.")
        print(filtered_data.variables)
        exists = False
        if exists :
            print('Derived data file exists' )
        else :

            field_list =sf.filter_variable_list(dataset, ref_dataset,
                                                derived_data, filtered_data,
                                                options,
                                                twod_filter, var_list=None,
                                                grid = opgrid)
    #        quad_field_list=list([])
            quad_field_list =sf.filter_variable_pair_list(dataset, ref_dataset,
                                                derived_data, filtered_data,
                                                options,
                                                twod_filter, var_list=None,
                                                grid = opgrid)


            d_r, d_s = sf.filtered_deformation(dataset, derived_data,
                                               filtered_data, options,
                                               twod_filter, dx, dy, z, zn,
                                               xaxis=1, grid='w')

            times = derived_data['time_series_50_100.0']
            print(times)
            print(times[:])
            for i in range(3) :
                for j in range(3) :
                    print(d_r[i][j],d_s[i][j])

            Sn_ij_r, mod_Sn_r = sf.shear(d_r)
            Sn_ij_s, mod_Sn_s = sf.shear(d_s)
            S_ij_r, mod_S_r = sf.shear(d_r, no_trace = False)
            S_ij_s, mod_S_s = sf.shear(d_s, no_trace = False)
            print(S_ij_r.keys())
#        input("Press enter")

        z = derived_data["z"]
        zn = derived_data["zn"]

#        u_r = derived_data["u_r_on"+opgrid]
#        print(u_r)
#        os.remove(derived_dataset_name)
#        print twod_filter
        if twod_filter.attributes['filter_type']!='domain' :
            fig1 = plt.figure(1)
            plt.contourf(twod_filter.data,20)
            plt.savefig(plot_dir+'Filter_'+\
                        twod_filter.id+plot_type)
            plt.close()

            fig2 = plt.figure(2)
            plt.plot(twod_filter.data[np.shape(twod_filter.data)[0]//2,:])
            plt.savefig(plot_dir+'Filter_y_xsect_'+\
                        twod_filter.id+plot_type)
            plt.close()

        for field in field_list:
            print("Plotting {}".format(field))
            plot_field(field, filtered_data, twod_filter, ilev, iy, grid=opgrid)

        for field in quad_field_list :
            print("Plotting {}".format(field))
            plot_quad_field(field, filtered_data, twod_filter, ilev, iy, \
                            grid=opgrid)

        plot_shear(mod_Sn_r, mod_Sn_s, z, twod_filter, ilev, iy, no_trace = True)
        plot_shear(mod_S_r, mod_S_s, z, twod_filter, ilev, iy, no_trace = False)

        filtered_data.close()
    derived_data.close()
    dataset.close()

if __name__ == "__main__":
    main()
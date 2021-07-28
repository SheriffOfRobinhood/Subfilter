# -*- coding: utf-8 -*-
"""

  cloud_analysis.py

Created on Wed May 15 09:46:31 2019

@author: Peter Clark

"""
import os
import netCDF4
from netCDF4 import Dataset
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import scipy.special as sp
import math as m

import subfilter as sf
import filters as filt
import difference_ops as do
from thermodynamics_constants import *
import thermodynamics as th
dir = 'C:/Users/paclk/OneDrive - University of Reading/Git/python/Subfilter/test_data/BOMEX/'
#dir = 'C:/Users/paclk/OneDrive - University of Reading/traj_data/r15/'
destdir = dir
odir = destdir
file = 'diagnostics_3d_ts_7200.0.nc'
ref_file = 'diagnostics_ts_7200.0.nc'

#w = dataset.variables['w']
#var_tvar = w.dimensions[0]
#var_time = dataset.variables[var_tvar]

plot_dir = dir+'plots/'

plot_type = '.png'
data_dir = '' # Directory containing data
figshow = True

fname = 'cloud'

options = {
#        'FFT_type': 'FFTconvolve',
#        'FFT_type': 'FFT',
        'FFT_type': 'RFFT',
        'save_all': 'Yes',
          }


def main():
    '''
    Top level code, a bit of a mess.
    '''
#   Non-global variables that are set once
#    sigma_list = [0.2,0.25,0.3,0.4,0.8,1.0,1.5,2.0] # Sigma used in Gaussian filter function
#    sigma_list = [2.0] # Sigma used in Gaussian#
    dataset = Dataset(dir+file, 'r') # Dataset is the class behavior to open the file
                                 # and create an instance of the ncCDF4 class
    ref_dataset=Dataset(dir + ref_file)
    w = dataset["w"]
    print(w.dimensions)

#    sigma_list = [0.5,0.2]
    sigma_list = [0.5]
    width = -1
    dx = 100.0
    dy = 100.0
    filter_name = 'gaussian'


    var_list = ["w", \
                "th", \
                "th_v", \
                "th_L", \
                "q_vapour", \
                "q_cloud_liquid_mass", \
                "q_total", \
                ]

    var_pair_list = [
                ["w","w"], \
                ["w","th"], \
                ["w","th_v"], \
                ["w","th_L"], \
                ["w","q_vapour"], \
                ["w","q_cloud_liquid_mass"], \
                ["w","q_total"], \
                ["th_L","th_L"], \
                ["th_L","q_total"], \
                ["q_total","q_total"], \
                ["th_L","q_vapour"], \
                ["th_L","q_cloud_liquid_mass"], \
              ]

    opgrid = 'p'
    #
    derived_dataset_name, derived_data, exists = \
                                        sf.setup_derived_data_file(
                                        dir+file, destdir, dir + ref_file,
                                        fname, options, override = False)

    print(derived_dataset_name)

    filter_list = list([])

    for i,sigma in enumerate(sigma_list):
        filter_id = 'filter_ga{:02d}'.format(i)
        twod_filter = filt.Filter(filter_id,\
                                       filter_name, \
                                       sigma=sigma, width=width, \
                                       delta_x=dx)

        print(twod_filter)
        filter_list.append(twod_filter)

    for twod_filter in filter_list:

        print(twod_filter)

        filtered_dataset_name, filtered_data, exists = \
            sf.setup_filtered_data_file( dir+file, odir, dir+ref_file, fname,
                                       options, twod_filter, override=True)


        if exists :
            print('File exists - opened for reading.')
        else :
            print('Creating derived data file.')
            field_list =sf.filter_variable_list(dataset, ref_dataset,
                                                derived_data, filtered_data,
                                                options, twod_filter,
                                                var_list=var_list,
                                                grid = opgrid)
#        quad_field_list=list([])
            quad_field_list =sf.filter_variable_pair_list(dataset,
                                                ref_dataset,
                                                derived_data, filtered_data,
                                                options, twod_filter,
                                                var_list=var_pair_list,
                                                grid = opgrid)
        times = derived_data[w.dimensions[0]]
        print(derived_data.variables)
        print(times)
        print(times[:])


#    derived_dataset_name = os.path.basename(file)
#    derived_dataset_name = ('.').join(derived_dataset_name.split('.')[:-1])
#    derived_dataset_name = derived_dataset_name + "_" + \
#            twod_filter_id + ".nc"
#    derived_dataset = Dataset(destdir+derived_dataset_name, "r")
        print(derived_data.variables)
        #    derived_data.twod_filter_id = twod_filter.id
        #    derived_data.setncatts(twod_filter.attributes)
        th_L_th_L=derived_data.variables["th_L_th_L_onp"][...]
        th_L_qt=derived_data.variables["th_L_q_total_onp"][...]
        qt_qt=derived_data.variables["q_total_q_total_onp"][...]
        th_L_r = derived_data.variables["th_L_r_onp"][...]
        th_L_s = derived_data.variables["th_L_s_onp"][...]
        q_t_r = derived_data.variables["q_total_r_onp"][...]
        q_t_s = derived_data.variables["q_total_s_onp"][...]
        q_cl_r = derived_data.variables["q_cloud_liquid_mass_r_onp"][...]
        q_cl_s = derived_data.variables["q_cloud_liquid_mass_s_onp"][...]
        pref = ref_dataset.variables['prefn'][-1,...]
        piref = (pref[:]/1.0E5)**kappa
        z = ref_dataset["z"][...]
        zn = ref_dataset["zn"][...]
#        piref_on_w = do.field_on_p_to_w(piref, z, zn)
#        p_on_w = 1.0E5 * (piref_on_w ** rk)
        p = 1.0E5 * (piref ** rk)
#        print(piref)
#        print(piref_on_w)
#        print(p_on_w)
#        print(p)
        T_L_r = th_L_r * piref
        alpha_L = th.dqsatbydT(T_L_r, p)

        print(np.max(alpha_L),np.min(alpha_L))
        a_L = 1.0 / (1.0 + L_over_cp * alpha_L)
        print(np.shape(a_L))
        print(np.max(a_L),np.min(a_L))

        Q_cr = a_L*(q_t_r - th.qsat(T_L_r, p))

        print(np.max(Q_cr),np.min(Q_cr))
#        plt.show()

        var_s = a_L * a_L * (qt_qt - 2 * alpha_L * piref * th_L_qt + \
                             alpha_L * alpha_L * piref * piref * th_L_th_L)

        sigma_s = np.sqrt(var_s)
        print(np.max(sigma_s),np.min(sigma_s))

        QN = Q_cr / sigma_s

        q_cl_N = q_cl_r / sigma_s

        th_L_q_s = alpha_L * piref * th_L_s

        s = a_L * (q_t_s - th_L_q_s)
        s_over_sigma_s = s / sigma_s

        QN_c=np.arange(-3,0.3, 0.1)
        C_MY=0.5*(1+sp.erf(QN_c/m.sqrt(2)))
        q_cl_MY = C_MY * QN_c + np.exp(-(QN_c**2)/2)/m.sqrt(2*m.pi)

#        print(np.shape(s), np.shape(q_cl_s))
        ilev = 15
        plt.figure(1)
        plt.plot(Q_cr[...,ilev].flatten(), q_cl_r[...,ilev].flatten(),'.',markersize=1)

        plt.figure(2)
        plt.plot(s[...,ilev].flatten(), q_cl_s[...,ilev].flatten(),'.',markersize=1)
        plt.show()


        for ilev in range(15,76) :
            print(zn[ilev])
#            print(np.min(th_L_q_s),np.max(th_L_q_s))
#            print(np.min(q_t_s),np.max(q_t_s))

            plt.figure(3)
            plt.hist2d(th_L_q_s[...,ilev].flatten()*1000, \
                       q_t_s[...,ilev].flatten()*1000, \
                       range = [[-0.1,0.1],[-1,1]],bins=(1000,1000))
            plt.xlabel(r"$\alpha_L \Pi \theta_L^s (g/kg)$")
            plt.ylabel(r"$q_t^s (g/kg)$")
            plt.title("Level {:2d} height {:4.0f}".format(ilev,zn[ilev]))
            plt.savefig(plot_dir+"q_t_v_th_L_{:2d}".format(ilev)+plot_type)
            plt.close()

#            plt.show()
            plt.figure(4)
            plt.hist( s[...,ilev].flatten()*1000, bins=100, density=True)
    #        plt.xlabel(r"$\alpha_L \Pi \theta_L^s (g/kg)$")
            plt.xlabel(r"$s^s (g/kg)$")
            plt.title("Level {:2d} height {:4.0f}".format(ilev,zn[ilev]))
            plt.savefig(plot_dir+"s_hist_{:2d}".format(ilev)+plot_type)
            plt.close()
#            plt.show()

            plt.figure(5)
            plt.hist( s_over_sigma_s[...,ilev].flatten(), bins=100, density=True)
    #        plt.xlabel(r"$\alpha_L \Pi \theta_L^s (g/kg)$")
            plt.xlabel(r"$s^s/\sigma_s$")
            plt.title("Level {:2d} height {:4.0f}".format(ilev,zn[ilev]))
            plt.savefig(plot_dir+"s_norm_hist_{:2d}".format(ilev)+plot_type)
            plt.close()

            plt.figure(6)
            plt.plot( QN[...,ilev].flatten(),q_cl_N[...,ilev].flatten(),'.',markersize=0.1)
            plt.plot(QN_c,q_cl_MY,'k')
            plt.ylabel(r"$q_{cl}N$")
            plt.xlabel(r"$QN$")
            plt.title("Level {:2d} height {:4.0f}".format(ilev,zn[ilev]))
            plt.savefig(plot_dir+"q_cl_vs_QN_{:2d}".format(ilev)+plot_type)
            plt.close()
#            plt.show()

        derived_data.close()

#    filter_dataset.close()

if __name__ == "__main__":
    main()



# -*- coding: utf-8 -*-
"""
Created on Mon Jul 26 16:01:57 2021

@author: paclk
"""
import numpy as np
#import pandas as pd
import xarray as xr
import matplotlib
import matplotlib.pyplot as plt

import subfilter as sf
import filters as filt

filter_name_list = [
 'gaussian',
 'running_mean',
 'running_mean_v2',
 'running_mean_v3',
 'wave_cutoff',
 'circular_wave_cutoff'
 ]

sigma = 220.0
dx = 5.0
npoints = 960

fig = plt.figure()
plot_type='.png'

filter_list = list([])

for i,filter_name in enumerate(filter_name_list):
    if filter_name == 'gaussian':
        filter_id = 'filter_ga{:02d}'.format(i)
        twod_filter = filt.Filter(filter_id,
                                  filter_name,
                                  npoints=npoints,
                                  sigma=sigma,
                                  delta_x=dx)
    elif filter_name == 'wave_cutoff':
        filter_id = 'filter_wc{:02d}'.format(i)
        twod_filter = filt.Filter(filter_id,
                                  filter_name,
                                  npoints=npoints,
                                  wavenumber=np.pi/(2*sigma),
                                  delta_x=dx)
    elif filter_name == 'circular_wave_cutoff':
        filter_id = 'filter_cwc{:02d}'.format(i)
        twod_filter = filt.Filter(filter_id,
                                  filter_name,
                                  npoints=npoints,
                                  wavenumber=np.pi/(2*sigma),
                                  delta_x=dx)
    elif filter_name == 'running_mean':
        filter_id = 'filter_rm{:02d}n'.format(i)
        width = int(np.round( sigma/dx * np.pi * 2.0 / 3.0)+1)
        #width = int(np.round( sigma/dx * 2.0 * np.sqrt(3.0))+1)
        twod_filter = filt.Filter(filter_id,
                                  filter_name,
                                  npoints=npoints,
                                  width=width,
                                  delta_x=dx)
    elif filter_name == 'running_mean_v2':
        filter_id = 'filter_rm{:02d}n'.format(i)
        #width = int(np.round( sigma/dx * np.pi * 2.0 / 3.0)+1)
        width = int(np.round( sigma/dx * 2.0 * np.sqrt(3.0))+1)
        twod_filter = filt.Filter(filter_id,
                                  'running_mean',
                                  npoints=npoints,
                                  width=width,
                                  delta_x=dx)
    elif filter_name == 'running_mean_v3':
        filter_id = 'filter_rm{:02d}n'.format(i)
        #width = int(np.round( sigma/dx * np.pi * 2.0 / 3.0)+1)
        #width = int(np.round( sigma/dx * 2.0 * np.sqrt(2*np.log(2)))+1)
        width = int(np.round( sigma/dx *  np.sqrt(2.0 *np.pi))+1)
        twod_filter = filt.Filter(filter_id,
                                  'running_mean',
                                  npoints=npoints,
                                  width=width,
                                  delta_x=dx)

    plt.plot(twod_filter.data[np.shape(twod_filter.data)[0]//2,:],
             label = filter_name)

    print(twod_filter)
    filter_list.append(twod_filter)

plt.legend()
plt.savefig('Filter_xsect'+plot_type)

plt.show()
#plt.close()



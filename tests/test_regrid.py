#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 21 10:47:00 2021

@author: xm904103
"""
import os

import numpy as np
import xarray as xr

import subfilter as sf
import dask

#from dask.diagnostics import Profiler, ResourceProfiler, CacheProfiler
#from dask.diagnostics import ProgressBar

from dask.distributed import Client

test_case = 1

if test_case == 0:
    options = {
    #        'FFT_type': 'FFTconvolve',
    #        'FFT_type': 'FFT',
            'FFT_type': 'RFFT',
            'save_all': 'Yes',
            'th_ref': 300.0,
            'dx': 50.0,
            'dy': 50.0,
              }
    dir = 'C:/Users/paclk/OneDrive - University of Reading/ug_project_data/Data/'
    odir = 'C:/Users/paclk/OneDrive - University of Reading/ug_project_data/Data/'
    odir = odir + 'test_dask_' + options['FFT_type']+'/'
    file = 'diagnostics_3d_ts_21600.nc'
    ref_file = 'diagnostics_ts_21600.nc'
elif test_case == 1:
    options = {
    #        'FFT_type': 'FFTconvolve',
    #        'FFT_type': 'FFT',
            'FFT_type': 'RFFT',
            'save_all': 'Yes',
            'th_ref': 300.0,
            'dx': 5.0,
            'dy': 5.0,
              }
    dir = '/storage/silver/wxproc/xm904103/'
#    dir = 'C:/Users/paclk/OneDrive - University of Reading/traj_data/CBL/'
    odir = '/storage/silver/wxproc/xm904103/'
#    odir = 'C:/Users/paclk/OneDrive - University of Reading/traj_data/CBL/'
    odir = odir + 'test_dask_' + options['FFT_type']+'/'
    file = 'diagnostics_3d_ts_13200.nc'
    ref_file = None

#dir = 'C:/Users/paclk/OneDrive - University of Reading/Git/python/Subfilter/test_data/BOMEX/'
#odir = 'C:/Users/paclk/OneDrive - University of Reading/Git/python/Subfilter/test_data/BOMEX/'

os.makedirs(odir, exist_ok = True)

#file = 'diagnostics_ts_18000.0.nc'
#ref_file = 'diagnostics_ts_18000.0.nc'

var_list = [
     "u",
     "v",
     "w",
     "th",
     "p",
     "tracer_rad1",
     ]


opgrid_list = [
    'w', 
    'p', 
    'u', 
    'v'
    ]

fname = 'test_regrid'

dask.config.set({"array.slicing.split_large_chunks": True})
dataset = xr.open_dataset(dir+file)
[itime, iix, iiy, iiz] = sf.find_var(dataset.dims, ['time', 'x', 'y', 'z'])
timevar = list(dataset.dims)[itime]
xvar = list(dataset.dims)[iix]
yvar = list(dataset.dims)[iiy]
zvar = list(dataset.dims)[iiz]
max_ch = sf.subfilter_setup['chunk_size']

nch = np.min([int(dataset.dims[xvar]/(2**int(np.log(dataset.dims[xvar]
                                            *dataset.dims[yvar]
                                            *dataset.dims[zvar]
                                            /max_ch)/np.log(2)/2))),
              dataset.dims[xvar]])

print(f'nch={nch}')

width = dataset.dims[xvar]

dataset.close()

defn = 1
#    dataset = xr.open_dataset(dir+file, chunks={timevar: defn,
#                                                'z':'auto', 'zn':'auto'})

dataset = xr.open_dataset(dir+file, chunks={timevar: defn,
                                            xvar:nch, yvar:nch,
                                            'z':'auto', 'zn':'auto'})
print(dataset)
#    ref_dataset = Dataset(dir+ref_file, 'r')
if ref_file is not None:
    ref_dataset = xr.open_dataset(dir+ref_file)
else:
    ref_dataset = None

derived_data, exists = \
    sf.setup_derived_data_file( dir+file, odir, fname,
                               options, override=True)

for opgrid in opgrid_list:
    for var_name in var_list:
        op_var = sf.get_data_on_grid(dataset, ref_dataset, derived_data, 
                     var_name, options, grid=opgrid)
        print(op_var)
print(derived_data)
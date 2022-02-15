#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 21 10:47:00 2021

@author: xm904103
"""
import os

import numpy as np
import xarray as xr

import subfilter.subfilter as sf
import subfilter.io.datain as di
from subfilter.utils.string_utils import get_string_index
import subfilter

import dask

#from dask.diagnostics import Profiler, ResourceProfiler, CacheProfiler
#from dask.diagnostics import ProgressBar

#from dask.distributed import Client
dask.config.set({"array.slicing.split_large_chunks": True})

test_case = 0

if test_case == 0:
    config_file = 'config_test_case_0.yaml'
    indir = 'C:/Users/paclk/OneDrive - University of Reading/ug_project_data/Data/'
    odir = 'C:/Users/paclk/OneDrive - University of Reading/ug_project_data/Data/'
    file = 'diagnostics_3d_ts_21600.nc'
    ref_file = 'diagnostics_ts_21600.nc'
elif test_case == 1:
    config_file = 'config_test_case_1.yaml'
    #indir = '/storage/silver/wxproc/xm904103/'
    indir = 'C:/Users/paclk/OneDrive - University of Reading/traj_data/CBL/'
    #odir = '/storage/silver/wxproc/xm904103/'
    odir = 'C:/Users/paclk/OneDrive - University of Reading/traj_data/CBL/'
    file = 'diagnostics_3d_ts_13200.nc'
    ref_file = None


options, update_config = sf.subfilter_options(config_file)

options['aliases'] = {
    'th':['theta', 'potential_temperature'],
    'p':['pressure'],
    'tracer': ['tracer_rad1'],
    }

odir = odir + 'test_dask_' + options['FFT_type']+'/'

#dir = 'C:/Users/paclk/OneDrive - University of Reading/Git/python/Subfilter/test_data/BOMEX/'
#odir = 'C:/Users/paclk/OneDrive - University of Reading/Git/python/Subfilter/test_data/BOMEX/'

os.makedirs(odir, exist_ok = True)

#file = 'diagnostics_ts_18000.0.nc'
#ref_file = 'diagnostics_ts_18000.0.nc'

var_list = [
    "saturation",
    "moist_dbdz",
    "dmoist_bdz",
    "dbdz",
    "dbdz_monc",
    "buoyancy_moist",
    "buoyancy",
    "buoyancy_monc",
    "u",
    "v",
    "w",
    "th",
    "p",
    "tracer",
    "th_L",
    "th_v",
    "thref",
    "T",
    "th_e",
    "th_w",
    "T_w",
    "m_vapour",
    "q_vapour",
    "cloud_fraction",
    "q_cloud_liquid_mass",
    "q_total",
    "m_total",
    "rh",
     ]


fname = 'test_datain'

dataset = xr.open_dataset(indir+file)

[itime, iix, iiy, iiz] = get_string_index(dataset.dims, ['time', 'x', 'y', 'z'])
timevar = list(dataset.dims)[itime]
xvar = list(dataset.dims)[iix]
yvar = list(dataset.dims)[iiy]
zvar = list(dataset.dims)[iiz]
max_ch = subfilter.global_config['chunk_size']

nch = np.min([int(dataset.dims[xvar]/(2**int(np.log(dataset.dims[xvar]
                                            *dataset.dims[yvar]
                                            *dataset.dims[zvar]
                                            /max_ch)/np.log(2)/2))),
              dataset.dims[xvar]])

print(f'nch={nch}')

width = dataset.dims[xvar]

dataset.close()

defn = 1
#dataset = xr.open_dataset(indir+file)
#    dataset = xr.open_dataset(dir+file, chunks={timevar: defn,
#                                                'z':'auto', 'zn':'auto'})

dataset = xr.open_dataset(indir+file, chunks={timevar: defn,
                                            xvar:nch, yvar:nch,
                                            'z':'auto', 'zn':'auto'})
print(dataset)
#    ref_dataset = Dataset(dir+ref_file, 'r')
if ref_file is not None:
    ref_dataset = xr.open_dataset(indir+ref_file)
else:
    ref_dataset = None

for var_name in var_list:
    op_var = di.get_data(dataset, ref_dataset, var_name, options,
                         allow_none=True)
    if op_var is None:
        print(f'{var_name:s} not found.')
    else:
        print(op_var.name, op_var.min().values, op_var.max().values,
              op_var.attrs['units'], op_var.shape)
        [itime, iix, iiy, iiz] = get_string_index(op_var.dims, ['time', 'x', 'y', 'z'])
        if iix is not None:
            tvar = op_var.dims[itime]
            xvar = op_var.dims[iix]
            yvar = op_var.dims[iiy]
            zvar = op_var.dims[iiz]
            op_var.isel({tvar:0, zvar:40}).plot.contourf(
                figsize=(12,10), levels=51, x=xvar)
            op_var.isel({tvar:0, yvar:95, zvar:slice(1,None)}).plot.contourf(
                figsize=(12,10), levels=51, x=xvar)

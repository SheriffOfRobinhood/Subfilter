# -*- coding: utf-8 -*-
"""
Created on Tue Feb 27 09:25:01 2024

@author: xm904103
"""

import os
import sys
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import dask
import monc_utils

import subfilter
import subfilter.subfilter as sf
import subfilter.filters as filt

from monc_utils.data_utils.string_utils import get_string_index
from monc_utils.io.dataout import save_field
from monc_utils.io.datain import (configure_model_resolution,
                                 get_data,
                                 get_data_on_grid,
                                 )
import monc_utils.data_utils.difference_ops as do
from loguru import logger

logger.remove()
logger.add(sys.stderr, 
           format = "<c>{time:HH:mm:ss.SS}</c>"\
                  + " | <level>{level:<8}</level>"\
                  + " | <green>{function:<22}</green> : {message}", 
           colorize=True, 
           level="INFO")
    
logger.enable("subfilter")
logger.enable("monc_utils")

logger.info("Logging 'INFO' or higher messages only")

fileroot = 'D:/Data/'
config_file = 'config_test_case_0.yaml'
indir = fileroot + 'ug_project_data/'
odir = fileroot + 'ug_project_data/'
file = 'diagnostics_3d_ts_21600.nc'
ref_file = 'diagnostics_ts_21600.nc'
options, update_config = sf.subfilter_options(config_file)

fname = 'test_mon_utils'

odir = odir + fname +'_'+ options['FFT_type']+'/'
os.makedirs(odir, exist_ok = True)

plot_dir = odir + 'plots/'
os.makedirs(plot_dir, exist_ok = True)
    
monc_utils.global_config['output_precision'] = "float32"

subfilter.global_config['test_level'] = 2

#dask.config.set({"array.slicing.split_large_chunks": True})
#dask_chunks = monc_utils.global_config['dask_chunks']


# Read data
# dataset = xr.open_dataset(indir+file, chunks=dask_chunks)
dataset = xr.open_dataset(indir+file)

print(dataset)

ref_dataset = xr.open_dataset(indir+ref_file)

print(ref_dataset)

z_w = dataset['z'].rename({'z':'z_w'})
z_p = dataset['zn'].rename({'zn':'z_p'})

# Get model resolution values
dx, dy, options = configure_model_resolution(dataset, options)

[itime, iix, iiy, iiz] = get_string_index(dataset.dims, ['time', 'x', 'y', 'z'])
[timevar, xvar, yvar, zvar] = [list(dataset.dims)[i] for i in [itime, iix, iiy, iiz]]

npoints = dataset.dims[xvar]

override = True

derived_data, exists = \
    sf.setup_derived_data_file( indir+file, odir, fname,
                               options, override=override)
if exists :
    print('Derived data file exists' )
    print("Variables in derived dataset.")
    print(derived_data['ds'].variables)
    
filter_name = 'gaussian'
sigma = 500.0
filter_id = 'filter_ga00'
twod_filter = filt.Filter(filter_id,
                          filter_name,
                          npoints=npoints,
                          sigma=sigma,
                          delta_x=dx)

var_list = ["th", "dbydx(th)"] 

for var_name in var_list:
    print("********************** Native *************************")
    var = get_data(dataset, ref_dataset, var_name) 
    
    print(var)
    
    dbydx_var = do.d_by_dx_field_native(var)
    
    print(dbydx_var)

    dbydy_var = do.d_by_dy_field_native(var)
    
    print(dbydy_var)
    
    dbydz_var = do.d_by_dz_field_native(var)
    
    print(dbydz_var)
    
    var_r, var_s = sf.filtered_field_calc(var, options, twod_filter)
    
    print(var_r)
    
    dbydx_var_r = do.d_by_dx_field_native(var_r)
    
    print(dbydx_var_r)
    
    for opgrid in ['p', 'w']:
    # for opgrid in []:
    
        print(f"********************** On Grid {opgrid} *************************")
    
        var = get_data_on_grid(dataset, ref_dataset, var_name, grid=opgrid) 
        
        print(var)
        
        dbydx_var = do.d_by_dx_field(var, z_w, z_p, grid=opgrid)
        
        print(dbydx_var)
    
        dbydy_var = do.d_by_dy_field(var, z_w, z_p, grid=opgrid)
        
        print(dbydy_var)
        
        dbydz_var = do.d_by_dz_field(var, z_w, z_p, grid=opgrid)
        
        print(dbydz_var)

        var_r, var_s = sf.filtered_field_calc(var, options, twod_filter)
        
        print(var_r)
        
        dbydx_var_r = do.d_by_dx_field(var_r, z_w, z_p, grid=opgrid)
        
        print(dbydx_var_r)
    

dataset.close()
ref_dataset.close()
derived_data['ds'].close()


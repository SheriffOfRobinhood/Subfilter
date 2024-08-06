# -*- coding: utf-8 -*-
"""
Created on Thu Nov  2 16:37:32 2023

@author: xm904103
"""

import os
import sys
from pathlib import Path
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt

import subfilter.subfilter as sf
import subfilter.filters as filt
import subfilter.spectra as sp

from monc_utils.data_utils.string_utils import get_string_index
from monc_utils.io.dataout import setup_child_file
from monc_utils.io.datain import configure_model_resolution

import dask
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

# override = False
override = True

var_list = [
    "w",
    ]

test_case = 0

def filter_and_spectra(input_data, odir, twod_filter, 
                       var_list=None,
                       options=None,
                       override=False,
                       opgrid='w'):
    logger.info("Processing using filter: ")
    logger.info(f"twod_filter \n{twod_filter}")
    file_tag = ""
    filtered_data, exists = \
        sf.setup_filtered_data_file( input_data['file'], odir, file_tag,
                                     options, twod_filter, 
                                     override=override)

    if exists :
        
        logger.info(f"Filtered data file {filtered_data['file']} exists" )
        
    else :
    
        if input_data.get('level', None) == 'filtered':

            dataset = input_data['ds']
           
            for varname in list(dataset.data_vars.keys()):
                var = dataset[varname]
                (var1_r, var1_s) = sf.filter_field(var, filtered_data, 
                                                   options, twod_filter)
                
        elif input_data.get('level', None) == 'raw':
            
            dataset = input_data['ds']
            if input_data['ref_dataset'] is None:
                ref_dataset = None
            else:               
                ref_dataset = input_data['ref_dataset']['ds']
            derived_data = input_data['derived_data']
            
            field_list = sf.filter_variable_list(dataset, ref_dataset,
                                                 derived_data, filtered_data,
                                                 options, twod_filter,
                                                 var_list=var_list,
                                                 grid = opgrid)


    logger.info(f'filtered_data \n{filtered_data}')
    
    filtered_data['ds'].close()

    filtered_data['ds'] = xr.open_dataset(filtered_data['file'],
                                    chunks={'time':1, 'z_'+opgrid:'auto'})

    outtag = "spectra"

    spectra_filt_ds, exists = setup_child_file(filtered_data['file'],
                                               odir, outtag,
                                               options, override=override)
    if exists :
        logger.info('Spectra of filtered data file exists.' )
    else :

        logger.info(f"Working on file: {spectra_filt_ds['ds']}")

        dso = sp.spectra_variable_list(filtered_data['ds'], 
                                       spectra_filt_ds,
                                       spectra_options,
                                       var_list=None)
    filtered_data['ds'].close()
    spectra_filt_ds['ds'].close()
    return (filtered_data['file'], spectra_filt_ds['file'])

def inv(k):
    return 1.0 / k


def get_k_Ek(dso, var):
    k_ref = dso['hfreq']
    
    #k_angular = dso['hwaven']
    
    kE_k = k_ref * dso[var]
    kE_km = kE_k.mean(dim='time')
    
    return kE_km

# dirroot = 'C:/Users/paclk/OneDrive - University of Reading/'
# dirroot = 'C:/Users/xm904103/OneDrive - University of Reading/'
dirroot = 'E:/Data/'
if test_case == 0:
    config_file = 'config_test_case_0.yaml'
    # indir = dirroot + 'ug_project_data/Data/'
    # odir = dirroot + 'ug_project_data/Data/'
    indir = dirroot + 'ug_project_data/'
    odir = dirroot + 'ug_project_data/'
    # files = 'diagnostics_3d_ts_21600.nc'
    # ref_file = 'diagnostics_ts_21600.nc'
    files = 'diagnostics_3d_ts_2[1-2]*.nc'
    ref_file = ''
    plot_height = 250
    sigma_list = [50.0, 100.0]
    ylim = [1E-6,1E-2]
elif test_case == 1:
    config_file = 'config_test_case_1.yaml'
    # indir = dirroot + 'traj_data/CBL/'
    # odir = dirroot + 'traj_data/CBL/'
    indir = dirroot + 'CBL/'
    odir = dirroot + 'CBL/'
    # odir = 'D:/Spectra/CBL/'
    files = 'diagnostics_3d_ts_*.nc'
    ref_file = None
    plot_height = 500  
    sigma_list = [10.0, 20.0, 40.0]
    ylim = [1E-4,1E-1]

options, update_config = sf.subfilter_options(config_file)
spectra_options, update_config = sp.spectra_options(
                                      'spectra_config.yaml')

odir = odir + 'repeat_filter/' 
os.makedirs(odir, exist_ok = True)

filter_name = 'gen_gaussian'
#    filter_name = 'running_mean'
#    filter_name = 'wave_cutoff'
# filter_name = 'circular_wave_cutoff'


dask.config.set({"array.slicing.split_large_chunks": True})

infiles = list(Path(indir).glob(files))

dataset = xr.open_dataset(infiles[0], chunks={'z':'auto', 'zn':'auto'})
# Get model resolution values
dx, dy, options = configure_model_resolution(dataset, options)

[itime, iix, iiy, iiz] = get_string_index(dataset.dims,
                                          ['time', 'x', 'y', 'z'])
timevar = list(dataset.dims)[itime]
xvar = list(dataset.dims)[iix]
yvar = list(dataset.dims)[iiy]
zvar = list(dataset.dims)[iiz]

npoints = dataset.dims[xvar]

dataset.close()

# Now create dicts of filter definitions.
# and output files.

filter_list = {}
filtered_data_files = {}
filtered_data_spectra_files = {}
filtered_data_files_2 = {}
filtered_data_spectra_files_2 = {}

for i,sigma in enumerate(sigma_list):
    filter_id = f'gga{i:02d}'
    twod_filter = filt.Filter(filter_id,
                              filter_name,
                              npoints=npoints,
                              sigma=sigma,
                              alpha=2.4,
                              delta_x=dx)
#        print(twod_filter)
    filter_list[filter_id] = twod_filter
    filtered_data_files[filter_id] = []
    filtered_data_spectra_files[filter_id] = []
    filtered_data_files_2[filter_id] = {}
    filtered_data_spectra_files_2[filter_id] = {}

for filter_id_1, twod_filter_1 in filter_list.items():
    for filter_id_2, twod_filter_2 in filter_list.items():
        filtered_data_files_2[filter_id_1][filter_id_2] = []
        filtered_data_spectra_files_2[filter_id_1][filter_id_2] = []

opgrid = 'w'
fname = 'derv'
outtag = "spectra"

derived_data_spectra_files = []

for infile in infiles:

    dataset = xr.open_dataset(infile, chunks={timevar:1, 'z':'auto', 'zn':'auto'})

    logger.info(f'Input file is {dataset}')

    input_data = {'file':infile, 'ds': dataset, 'level':'raw'}
    
    if ref_file is not None:
        
        ref_file = infile.name.split('_')
        ref_file.pop(ref_file.index('3d'))   
        ref_file = "_".join(ref_file)
        ref_dataset = xr.open_dataset(indir+ref_file)
        input_data['ref_dataset'] = {'file':ref_file, 'ds':ref_dataset}
    else:
        ref_dataset = None
        input_data['ref_dataset'] = None
        

    derived_data, exists = \
        sf.setup_derived_data_file( infile, odir, fname,
                                   options, override=override)
    if exists :
        logger.info('Derived data file exists' )
        # logger.info("Variables in derived dataset.")
        # logger.info(f"{derived_data['ds'].variables}")
        
    input_data['derived_data'] = derived_data

# Process data with each filter.

########################### First filter ######################################

    for filter_id, twod_filter in filter_list.items():
        
        (filtered_data_file, spectra_filt_file) = \
                filter_and_spectra(input_data, 
                                   odir, 
                                   twod_filter, 
                                   var_list=var_list,
                                   options=options,
                                   override=override,
                                   opgrid=opgrid)    

        filtered_data_files[filter_id].append(filtered_data_file)
        filtered_data_spectra_files[filter_id].append(
            spectra_filt_file)
        
############################## Tidy up ########################################       
    dataset.close()
    derived_data['ds'].close()
    
######################Filter Derived Data from main file#######################

    derived_data['ds'] = xr.open_dataset(derived_data['file'],
                                         chunks={'time':1, 'z_'+opgrid:'auto'})

    spectra_derived_ds, exists = setup_child_file(derived_data['file'], odir,
                         outtag, options, override=override)

    if exists :
        logger.info('Spectra of derived data file exists' )
    else :
        dso = sp.spectra_variable_list(derived_data['ds'], 
                                       spectra_derived_ds, 
                                       spectra_options,
                                       var_list=None)

    derived_data_spectra_files.append(spectra_derived_ds['file'])        
    
    logger.info(f'derived_data \n{derived_data}')
    
    logger.info(f'spectra_derived_ds \n{spectra_derived_ds}')
    
    derived_data['ds'].close()
    spectra_derived_ds['ds'].close()

        
#%% Set up axes
fig1, axa = plt.subplots(1,1,figsize=(8,6))

#%% Get time meaned spectrum

dso = xr.open_mfdataset(derived_data_spectra_files)
#print(dso)

kE_km = get_k_Ek(dso, 'spec_2d_w_on_w')

#%% Plot Energy Density Spectrum
kE_kmz = kE_km.sel(z_w=plot_height, method='nearest')
kE_kmz.plot(xscale='log', yscale='log', label=rf'Ref $\Delta=${dx} m', ax=axa)

dso.close()
#%%
for filter_id, spectra_files in filtered_data_spectra_files.items():
    dso = xr.open_mfdataset(spectra_files)
#    print(dso)
    
    kE_km = get_k_Ek(dso, 'spec_2d_f(w_on_w)_r')
    

    kE_kmz = kE_km.sel(z_w=plot_height, method='nearest')
    sigma = filter_list[filter_id].attributes["sigma"]
    kE_kmz.plot(xscale='log', yscale='log', 
                label=rf'$\sigma$={sigma} m', ax=axa)
    
    dso.close()
        
  
axa.set_ylabel(r'$kE(k)$')
axa.legend()
#axa[2].set_xlim(xlims2)
axa.set_ylim(ylim)
axa.set_xlabel(r'Wavenumber m$^{-1}$')
#axa.set_ylabel(r'$G(k)\times G(k)^*$')
secax = axa.secondary_xaxis('top', functions=(inv, inv))
secax.set_xlabel('wavelength (m)')
plt.tight_layout()

plt.savefig(odir + 'test_gga_filter.pdf')


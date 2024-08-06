# -*- coding: utf-8 -*-
"""
    This program evaluates MONC fields with 4 dimensions (vertical, x, y, time) to produce
        horizontal power spectra at each time and vertical level written to new netcdf files.

    BY DEFAULT, each variable contained in the input file's xarray.Dataset "Data variables"
        has its horizontal power spectra evaluated.  These are all placed in the same
        output file.
    They can alternatively be placed in a list in the user settings section.

    Several options can influence the form of the final result.

    Assumes the horizontal grid dimensions are the same for each variable being analysed.
    Assumes the horizontal dimensions are named 'x' and 'y'.
    Assumes the vertical dimension is the only dimension with a 'z' in its name, but it can be
        either 'z' or 'zn'.
    The time dimension name is identified by a user-supplied string, currently: 'time'.

    "Durran" calculation based on Durran et al. (2017): https://doi.org/10.1175/MWR-D-17-0056.1

    User must supply:
       dir:    input directory (slash-agnostic)
       file:   input file
                 Suggest switching to argument input (see below)
       outtag: output file tag (appended to input file name)
                 Creates 'spectra/' directory within the given dir
       dx:     x-direction grid spacing [m]
       dy:     y-direction grid spacing [m]

    @author: Todd Jones
"""

import os
import sys
import glob
import xarray as xr
import dask


import monc_utils
from monc_utils.io.dataout import setup_child_file
from monc_utils.io_um.datain import get_um_data_on_grid, stash_map

import matplotlib.pyplot as plt
from loguru import logger

#monc_utils.global_config['no_dask']=True
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

# ---------------------------------- USER SETTINGS ----------------------------------

# File source location and name

dirroot = 'E:/Data/'

var_names_spec = [
                    'u',
                    'v',
                    # 'w',
                  ]  # list of variable names to evaluate
        #   - leave empty to work on all present variables
        #   - populate with string variable names to work on specified list
# ---------------------------------- USER SETTINGS ----------------------------------

# Add to output file name
outtag = 'spec'

runs = ['100m',
        # '200m',
        # '250m',
        # '500m',
        # '1km',
        # '10km',
        ]

indir = dirroot + 'UM_CBL/'
outdir = dirroot + 'UM_CBL/'
outdir = outdir + 'test_filter_spectra/'
os.makedirs(outdir, exist_ok = True)

dask.config.set({"array.slicing.split_large_chunks": True})

#################################################################################################
#################################################################################################
#################################################################################################
    
    
run = '100m'

file = f'dc455_{run}_L100a_pr000.nc'

if 'k' in run:
    ind = run.index('k')
    dx = dy = float(run[0:ind])*1000
else:
    ind = run.index('m')
    dx = dy = float(run[0:ind])
    
options = {'dx': dx, 'dy': dy, 'th_ref': 300.0}
        
print(options)
   
infile = indir+file

dataset_in = xr.open_dataset(infile, chunks='auto') # , chunks={'z':'auto', 'zn':'auto'})

print(dataset_in)

for c in dataset_in.coords:
    if 'longitude' in c or 'latitude' in c: 
        print(c, dataset_in[c])

# By default, we are going to look at all data variables from the file
#   Save data variable names to list
# ALTERNATIVELY, THE CODE COULD BE CONFIGURED TO PASS IN A LIST OF
# SPECIFIC FIELDS TO EVALUATE.
if len(var_names_spec) == 0:
    var_names = list(dataset_in.data_vars.keys())
else:
    var_names = var_names_spec.copy()
    
derived_dataset, exists = setup_child_file(infile, outdir, 'derv',
                                           options, 
                                           override=True)#,
                                           #keep_coords= keep_coords)

print(derived_dataset)
    
for var in var_names:
    
    dataset_in[f'STASH_{stash_map[var]}'].isel(min15T0_0=8, rholev_eta_rho=17).plot.contourf()
    
    plt.show()
#                dataset[var] = get_um_and_transform(dataset_in, var, grid='p')
    field = get_um_data_on_grid(dataset_in, var, 
                                       derived_dataset,
                                       options,
                                       rename_time=True,
                                       grid='p')
    print(derived_dataset)
    print(field.isel(time=8, z_p=17).max())
    
    field.isel(time=8, z_p=17).plot.contourf()
    plt.show()
    
    field = 0
    
derived_dataset['ds'].close()
dataset_in.close()

dataset = xr.open_dataset(derived_dataset['file'],chunks='auto')

for c in dataset.coords:
    print(c, dataset[c])


# dataset = dataset.isel(time_0=slice(22, None, None), z_p=slice(17, 19, None))

keep_coords = {'elapsed_time':dataset.coords['elapsed_time']}
   
print(dataset)
for var in var_names:
#    print(dataset[f'{var}_on_p'].isel(time=8, z_p=17).max())  
    dataset[f'{var}_on_p'].isel(time=8, z_p=17).plot.contourf()
    plt.show()
#        ds_coords = {k:dataset.coords[k] for k in list(dataset.dims.keys())}
dataset.close()

print("IAMDONE")  # Tag for successful completion.

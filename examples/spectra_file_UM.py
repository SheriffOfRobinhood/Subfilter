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


import subfilter.spectra as sp
import monc_utils
from monc_utils.io.dataout import setup_child_file
from monc_utils.io_um.datain import get_um_and_transform, get_um_data_on_grid
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
                    'w',
                  ]  # list of variable names to evaluate
        #   - leave empty to work on all present variables
        #   - populate with string variable names to work on specified list
# ---------------------------------- USER SETTINGS ----------------------------------

# Add to output file name
outtag = 'spec'

runs = ['100m',
        '200m',
        '250m',
        '500m',
        '1km',
        '10km',
        ]

indir = dirroot + 'UM_CBL/'
outdir = dirroot + 'UM_CBL/'
outdir = outdir + 'test_filter_spectra/'
os.makedirs(outdir, exist_ok = True)

dask.config.set({"array.slicing.split_large_chunks": True})

#################################################################################################
#################################################################################################
#################################################################################################
def main():
    '''
    Top level code
    '''
    
    for run in runs:
        files = f'dc455_{run}_L100a_pr000.nc'
        
        if 'k' in run:
            ind = run.index('k')
            dx = dy = float(run[0:ind])*1000
        else:
            ind = run.index('m')
            dx = dy = float(run[0:ind])
            
        options = {'dx': dx, 'dy': dy, 'th_ref': 300.0}
        
        spectra_options, update_config = sp.spectra_options()
        spectra_options.update({'dx': dx, 'dy':dy, 'spec_1D':False})
        
        print(spectra_options)
        
        options.update(spectra_options)
        
        print(options)
           
        infiles = glob.glob(indir+files)
    
        for infile in infiles:
            
            dataset_in = xr.open_dataset(infile)#, chunks='auto') # , chunks={'z':'auto', 'zn':'auto'})
            
            print(dataset_in)
            
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
                field = get_um_data_on_grid(dataset_in, var, 
                                                   derived_dataset,
                                                   options,
                                                   rename_time=True,
                                                   grid='p')
                print(derived_dataset)
                
            derived_dataset['ds'].close()
            dataset_in.close()
            
            dataset = xr.open_dataset(derived_dataset['file'],chunks='auto')
                        
            keep_coords = {'elapsed_time':dataset.coords['elapsed_time']}
   
            print(dataset)
            for var in var_names:
                print(dataset[f'{var}_on_p'].isel(time=8, z_p=17).max())    
    #        ds_coords = {k:dataset.coords[k] for k in list(dataset.dims.keys())}
    
            spectrum_dataset, exists = setup_child_file(infile, outdir, f'{outtag}',
                                                        options, 
                                                        override=True,
                                                        keep_coords= keep_coords)

            for var in var_names:
                
                print(spectrum_dataset)
                print(f"Working on file: {infile}")
        
                of = sp.spectra_variable_list(dataset, spectrum_dataset, options,
                                              var_list=[f'{var}_on_p'])
            dataset.close()
            spectrum_dataset['ds'].close()
            
            dso =  xr.open_dataset(spectrum_dataset['file'])
            print(dso)
            for var in var_names:
                # dso = xr.open_dataset(dd[var]['file'])
                print(dso[f'spec_2d_{var}_on_p'].isel(time=8, z_p=17).max())
                dso.close()

    print("IAMDONE")  # Tag for successful completion.

# END main


if __name__ == "__main__":
    main()

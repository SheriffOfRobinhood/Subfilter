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
import glob
import numpy as np
import xarray as xr
import dask


import subfilter
import subfilter.subfilter as sf
import subfilter.spectra as sp
from monc_utils.data_utils.string_utils import get_string_index
from monc_utils.io.dataout import setup_child_file
from monc_utils.io.datain import configure_model_resolution

test_case = 1

# ---------------------------------- USER SETTINGS ----------------------------------

# Single file source location and name
# dir = '/gws/nopw/j04/paracon_rdg/users/toddj/BOMEX/CA_S_SC_BOMEX_25_600/f3/test_op_RFFT/'
# dir = '/work/scratch-pw/toddj/'  # testing
# file = 'BOMEX_m0025_g0600_all_66600.0_test_plot.nc'
# -- or --
#   Modify script to receive file and/or path name

# dir = 'C:/Users/paclk/OneDrive - University of Reading/ug_project_data/Data/'
# file = 'diagnostics_3d_ts_21600.nc'


# dirroot = 'C:/Users/paclk/OneDrive - University of Reading/'
dirroot = 'C:/Users/xm904103/OneDrive - University of Reading/'


time_dim_always_contains='time'
var_names_spec = []  # list of variable names to evaluate

var_names_spec = ['w']  # list of variable names to evaluate
        #   - leave empty to work on all present variables
        #   - populate with string variable names to work on specified list
# ---------------------------------- USER SETTINGS ----------------------------------

outtag = 'w_spec'

opgrid = 'w'

#################################################################################################
#################################################################################################
#################################################################################################
def main():
    '''
    Top level code
    '''
    if test_case == 0:
        config_file = 'config_test_case_0.yaml'
        indir = dirroot + 'ug_project_data/Data/'
        outdir = dirroot + 'ug_project_data/Data/'
        files = 'diagnostics_3d_ts_21600.nc'
        ref_file = 'diagnostics_ts_21600.nc'
    elif test_case == 1:
        config_file = 'config_test_case_1.yaml'
        indir = dirroot + 'traj_data/CBL/'
        outdir = dirroot + 'traj_data/CBL/'
        files = 'diagnostics_3d_ts_*.nc'
        ref_file = None

    options, update_config = sf.subfilter_options(config_file)
    spectra_options, update_config = sp.spectra_options(
                                          'spectra_config.yaml')
    
    options.update(spectra_options)
    
    print(options)
    

    outdir = outdir + 'test_filter_spectra_' + options['FFT_type']+'/'
    os.makedirs(outdir, exist_ok = True)

    override=options['override']

    dask.config.set({"array.slicing.split_large_chunks": True})

    infiles = glob.glob(indir+files)

    for infile in infiles:
        
        dataset = xr.open_dataset(infile, chunks={'z':'auto', 'zn':'auto'})

        print(dataset)
        if ref_file is not None:
            ref_dataset = xr.open_dataset(indir+ref_file)
        else:
            ref_dataset = None

        # Get model resolution values
        dx, dy, options = configure_model_resolution(dataset, options)

        [itime, iix, iiy, iiz] = get_string_index(dataset.dims,
                                                  ['time', 'x', 'y', 'z'])
        timevar = list(dataset.dims)[itime]
        xvar = list(dataset.dims)[iix]
        yvar = list(dataset.dims)[iiy]
        zvar = list(dataset.dims)[iiz]

        npoints = dataset.dims[xvar]

    # By default, we are going to look at all data variables from the file
    #   Save data variable names to list
    # ALTERNATIVELY, THE CODE COULD BE CONFIGURED TO PASS IN A LIST OF
        # SPECIFIC FIELDS TO EVALUATE.
        if len(var_names_spec) == 0:
            var_names = list(dataset.data_vars.keys())
        else:
            var_names = var_names_spec.copy()

        derived_dataset, exists = setup_child_file(infile, outdir, outtag,
                                                   options, 
                                                   override=options['override'])

        print(f"Working on file: {infile}")

        of = sp.spectra_variable_list(dataset, derived_dataset, options,
                                        var_list=var_names)
        print(of)
        
        dso = xr.open_dataset(derived_dataset['file'])
    
        print(dso)
    
        dso.close()

    print("IAMDONE")  # Tag for successful completion.

# END main


if __name__ == "__main__":
    main()

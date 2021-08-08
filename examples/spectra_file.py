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
import numpy as np
import xarray as xr
import dask


import subfilter as sf
import subfilter.spectra as sp
from subfilter.utils.string_utils import get_string_index
from subfilter.io.dataout import setup_child_file


# ---------------------------------- USER SETTINGS ----------------------------------

# Single file source location and name
# dir = '/gws/nopw/j04/paracon_rdg/users/toddj/BOMEX/CA_S_SC_BOMEX_25_600/f3/test_op_RFFT/'
# dir = '/work/scratch-pw/toddj/'  # testing
# file = 'BOMEX_m0025_g0600_all_66600.0_test_plot.nc'
# -- or --
#   Modify script to receive file and/or path name

# dir = 'C:/Users/paclk/OneDrive - University of Reading/ug_project_data/Data/'
# file = 'diagnostics_3d_ts_21600.nc'
# dx = 50.0
# dy = 50.0

indir = 'C:/Users/paclk/OneDrive - University of Reading/traj_data/CBL/'
file = 'diagnostics_3d_ts_13200.nc'

#if len(sys.argv) > 2:
#    dir  = sys.argv[1]
#    file = sys.argv[2]

# Set up outfile
outdir = os.path.join(indir, 'spectra/')
outtag = "spectra_w_2D"
# Obtain dx=dy from file name '_m' encoding
#mnc = infile.find('_m')
#dx = float(infile[mnc+2:mnc+6])
#dy = float(infile[mnc+2:mnc+6])


time_dim_always_contains='time'
var_names_spec = []  # list of variable names to evaluate

var_names_spec = ['w']  # list of variable names to evaluate
        #   - leave empty to work on all present variables
        #   - populate with string variable names to work on specified list
# ---------------------------------- USER SETTINGS ----------------------------------



#################################################################################################
#################################################################################################
#################################################################################################
def main():
    '''
    Top level code
    '''

    dask.config.set({"array.slicing.split_large_chunks": True})
    infile = os.path.join(indir, file)

    ds = xr.open_dataset(infile, chunks={'z':'auto', 'zn':'auto'})

    # By default, we are going to look at all data variables from the file
    #   Save data variable names to list
    # ALTERNATIVELY, THE CODE COULD BE CONFIGURED TO PASS IN A LIST OF
        # SPECIFIC FIELDS TO EVALUATE.
    if len(var_names_spec) == 0:
        var_names = list(ds.data_vars.keys())
    else:
        var_names = var_names_spec.copy()

    options, update_config = sp.spectra_options('spectra_config.yaml')

    os.makedirs(outdir, exist_ok = True)  # make outdir if it doesn't already exist
    derived_dataset, exists = setup_child_file(infile, outdir, outtag,
                            options, override=options['override'])

    print("Working on file: "+file)

    dso = sp.spectra_variable_list(ds, derived_dataset, options,
                                        var_list=var_names)

    dso = xr.open_dataset(derived_dataset['file'])

    print(dso)

    dso.close()

    print("IAMDONE")  # Tag for successful completion.

# END main


if __name__ == "__main__":
    main()

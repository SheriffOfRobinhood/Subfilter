# -*- coding: utf-8 -*-
"""
Created on Tue Oct 23 11:27:25 2018

@author: Peter Clark
"""
import os
import glob
#import netCDF4
#from netCDF4 import Dataset
import numpy as np
#import pandas as pd
import xarray as xr
import matplotlib
import matplotlib.pyplot as plt

import subfilter
import subfilter.subfilter as sf
import subfilter.filters as filt
import subfilter.spectra as spectra
from subfilter.utils.string_utils import get_string_index
from subfilter.io.dataout import setup_child_file
from subfilter.io.datain import configure_model_resolution

#import difference_ops as do
import dask

#from dask.diagnostics import Profiler, ResourceProfiler, CacheProfiler
#from dask.diagnostics import ProgressBar

from dask.distributed import Client

test_case = 1

def main():
    '''
    Top level code, a bit of a mess.
    '''
    if test_case == 0:
        config_file = 'config_test_case_0.yaml'
        indir = 'C:/Users/paclk/OneDrive - University of Reading/ug_project_data/Data/'
        odir = 'C:/Users/paclk/OneDrive - University of Reading/ug_project_data/Data/'
        files = 'diagnostics_3d_ts_21600.nc'
        ref_file = 'diagnostics_ts_21600.nc'
    elif test_case == 1:
        config_file = 'config_test_case_1.yaml'
        indir = 'C:/Users/paclk/OneDrive - University of Reading/traj_data/CBL/'
        odir = 'C:/Users/paclk/OneDrive - University of Reading/traj_data/CBL/'
        files = 'diagnostics_3d_ts_*.nc'
        ref_file = None

    options, update_config = sf.subfilter_options(config_file)
    spectra_options, update_config = spectra.spectra_options(
                                          'spectra_config.yaml')

    odir = odir + 'test_filter_spectra_' + options['FFT_type']+'/'
    os.makedirs(odir, exist_ok = True)

    override=options['override']

#    client = Client()
#    client
#    dask.config.set(scheduler='threads')
#   Non-global variables that are set once

    var_list = [
        # "u",
        # "v",
        "w",
        # "th",
        ]

    sigma_list = [5.0, 10.0, 20.0, 40.0, 80.0, 160.0, 320.0]
#    sigma_list = [10.0, 20.0]
#    filter_name = 'gaussian'
#    filter_name = 'running_mean'
#    filter_name = 'wave_cutoff'
    filter_name = 'circular_wave_cutoff'

#    dataset = Dataset(indir+file, 'r') # Dataset is the class behavior to open the file
                                 # and create an instance of the ncCDF4 class
    dask.config.set({"array.slicing.split_large_chunks": True})

    infiles = glob.glob(indir+files)

    for infile in infiles:
        dataset = xr.open_dataset(infile, chunks={'z':'auto', 'zn':'auto'})

        # dataset = xr.open_dataset(indir+file, chunks={timevar: defn,
        #                                             xvar:nch, yvar:nch,
        #                                             'z':'auto', 'zn':'auto'})
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

        opgrid = 'w'
        fname = 'filter_spectra'

        derived_data, exists = \
            sf.setup_derived_data_file( infile, odir, fname,
                                       options, override=override)
        if exists :
            print('Derived data file exists' )
            print("Variables in derived dataset.")
            print(derived_data['ds'].variables)

    # Now create list of filter definitions.

        filter_list = list([])

        for i,sigma in enumerate(sigma_list):
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
    #            filter_id = 'filter_rm{:02d}'.format(i)
    #            filter_id = 'filter_rm{:02d}n'.format(i)
                filter_id = 'filter_rm{:02d}v3'.format(i)
    #            width = int(np.round( sigma/dx * np.pi * 2.0 / 3.0)+1)
    #            width = int(np.round( sigma/dx * 2.0 * np.sqrt(3.0))+1)
                width = int(np.round( sigma/dx *  np.sqrt(2.0 *np.pi))+1)
                twod_filter = filt.Filter(filter_id,
                                          filter_name,
                                          npoints=npoints,
                                          width=width,
                                          delta_x=dx)

    #        print(twod_filter)
            filter_list.append(twod_filter)

    # Process data with each filter.

        filtered_data_files = []
        filtered_data_spectra_files = []

        for twod_filter in filter_list:

            print("Processing using filter: ")
            print(twod_filter)

            filtered_data, exists = \
                sf.setup_filtered_data_file( infile, odir, fname,
                                           options, twod_filter, override=override)

            if exists :
                print('Filtered data file exists' )
                print("Variables in filtered dataset.")
                print(filtered_data['ds'].variables)
            else :
                field_list =sf.filter_variable_list(dataset, ref_dataset,
                                                    derived_data, filtered_data,
                                                    options, twod_filter,
                                                    var_list=var_list,
                                                    grid = opgrid)

            filtered_data_files.append(filtered_data['file'])

            filtered_data['ds'].close()

            filtered_data['ds'] = xr.open_dataset(filtered_data['file'],
                                                  chunks={'z':'auto', 'zn':'auto'})
            ds = filtered_data['ds']

            outtag = "spectra_w_2D"

            spectra_filt_ds, exists = setup_child_file(filtered_data['file'],
                                                       odir, outtag,
                                                       options, override=override)

            print("Working on file: "+spectra_filt_ds['ds'])

            dso = spectra.spectra_variable_list(ds, spectra_filt_ds,
                                                spectra_options,
                                                var_list=None)



            filtered_data_spectra_files.append(spectra_filt_ds['file'])

            filtered_data['ds'].close()

        dataset.close()
        derived_data['ds'].close()

        derived_data['ds'] = xr.open_dataset(derived_data['file'],
                                             chunks={'z':'auto', 'zn':'auto'})
        ds = derived_data['ds']


        spectra_derived_ds, exists = setup_child_file(derived_data['file'], odir,
                             outtag, options, override=options['override'])

        dso = spectra.spectra_variable_list(ds, spectra_derived_ds, spectra_options,
                                            var_list=None)

        print('--------------------------------------')

        print(derived_data)

        print('--------------------------------------')
        derived_data['ds'].close()

if __name__ == "__main__":
    main()

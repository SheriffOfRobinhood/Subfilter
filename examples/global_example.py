# -*- coding: utf-8 -*-
"""
Created on Tue Oct 23 11:27:25 2018

@author: Todd Jones


    Use of this script requires inspection of global_example.yaml
    That file specifies the input data and many options.
    Some of the options present are not configured by default for use
      in subfilter routines but are co-opted by this script.



"""
import os
import sys
import getopt
import glob
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from distutils.util import strtobool

import subfilter.subfilter as sf
import subfilter.filters as filt
import subfilter.spectra as spectra
from subfilter.utils.string_utils import get_string_index
from subfilter.io.dataout import setup_child_file
from subfilter.io.datain import configure_model_resolution

import subfilter

import dask

print(f'\nExecuting: {__file__}\n', flush=True)

# Set global defaults
testing = False
driver_config = None
default_config = os.path.dirname(subfilter.__file__)+'/../examples/global_example.yaml'
# Configure bespoke variable lists for analysis
var_list_full = ["u", "v", "w", "th", "q_vapour", "th_v", "th_L", "q_total", "buoyancy"]
var_pair_list_full = [
                ["u","u"],
                ["u","v"],
                ["u","w"],
                ["v","v"],
                ["v","w"],
                ["w","w"] ]


var_list_test = ["u", "w", "buoyancy"]  # lighter option
var_pair_list_test = [
                ["u","u"],
                ["u","v"]  ]

var_list = None          # DEFAULT: use yaml's test_level
var_pair_list = None



# Usage Message
def print_usage_message():
    print('Usage of '+sys.argv[0])
    print('  '+sys.argv[0]+' -y <driver yaml path> -v <var_list_ID>')
    print('    -y or --driver_yaml:   Text string path to input driver file')
    print("    -v or --var_list_key:  Text string ID ['test', 'full'] to override test_level variable lists")
    sys.exit()

# Simple yes/no question 
def yes_or_no(question: str):
    reply = None
    while reply not in ("y", "n"):
        reply = input(f"{question} [y,n]: ").strip().lower()
    return (reply == "y")

# Read in user-supplied options arguments
def read_cl_arguments():

    # Possibly modifying globals:
    global driver_config
    global var_list
    global var_pair_list

    test_override = ''

    if len(sys.argv) > 1 :
        in_options, args = getopt.getopt(sys.argv[1:],"hy:v:", \
                           ["driver_yaml="])
        for opt, arg in in_options:
            if opt == '-h':
                print_usage_message()
            if opt in ('-y', '--driver_yaml'):
                driver_config = arg
            if opt in ('-v', '--var_list_key'):
                test_override = arg

    else:
        print('No option arguments provided.')
        response = True #yes_or_no(f'Continue with default driver config ({default_config})?')
        if response:
            driver_config = default_config
        else:
            print('YAML Configuration required but not provided.  Exiting.')
            sys.exit()

    if driver_config is None:
        driver_config = default_config

    if test_override.lower() == 'test':
        var_list = var_list_test
        var_pair_list = var_pair_list_test
    if test_override.lower() == 'full':
        var_list = var_list_full
        var_pair_list = var_pair_list_full
    if var_list is not None:
        print(f'Using var_list:\n    {var_list}')
        print(f'Using var_pair_list:\n    {var_pair_list}\n')
    else:
        print(f'Using var_list and var_pair_list as specified by test_level in subfilter.global_config')
        
# Generic raise function
def _raise(exception_type, msg):
    raise exception_type(msg)

############################################################################################################
def main():
    '''
    Top level code
    '''
    # Flush print
    print(f'Beginning run of: {__file__}', flush=True)

    # Read in user-supplied command line option arguments
    read_cl_arguments()
    print(f'Using options provided by: {driver_config}\n')
    if not os.path.isfile(driver_config):
        raise ValueError(f'Options file, {driver_config}, does not exist.')

    # Update global_config from yaml
    subfilter.set_global_config(driver_config)
    if var_list is None:
        print(f'test_level: {subfilter.global_config["test_level"]}\n')

    # Configure via driver yaml file (name supplied as option to this script)
    # Load options dictionary using default values and changes from yaml file
    # sf_config is a record of the yaml contents, which may contain other option lists 
    #  (filters, spectra)
    options, sf_config = sf.subfilter_options(driver_config)

    # Sort through the input options
    input_file = options['input_file'] if options['input_file'] is not None \
                     else _raise(ValueError, 'No input_file in config')
    # Check for wildcards
    if any(substring in input_file for substring in ('*', '?', '[')):
        infiles = glob.glob(input_file)
        print("All output for this file list will be written to the same directory, outpath.")
        print(infiles)
    else:
        infiles = [input_file]

    ref_file = options['ref_file']
    if ref_file is None:
        ref_file = infiles[0]
        options['ref_file'] = ref_file
        print(' WARN: No ref_file provided; trying with first input_file.')

    # Handle output directories
    outpath = options['outpath'] if options['outpath'] is not None \
                  else _raise(ValueError, 'No outpath in config')
    outpath = os.path.join(outpath,'')     # ensure trailing slash
    os.makedirs(outpath, exist_ok = True)

    # Set append file label for output data files.
    fname = options['file_label']

    # Overwrite any existing output files
    override = options['override']

    # Filter options
    filter_options = sf_config['filters']
    filter_name = filter_options['filter_name']
    sigma_list = filter_options['sigma_list']
    sigma_list.sort()
    germano_list = filter_options['germano_list']
    germano_list.sort()
    include_domain_mean = filter_options['include_domain_mean']
    opgrid = filter_options['output_grid']
    run_quad_fields = filter_options['run_quad_fields']
    run_deformation_fields = filter_options['run_deformation_fields']

    # Spectra options
    spectra_options = sf_config['spectra_options']
    outpath_spectra = spectra_options['outpath_spectra']
    if outpath_spectra is None:
        outpath_spectra = os.path.join(outpath,'spectra/')
        os.makedirs(outpath_spectra, exist_ok = True)
        spectra_options['outpath_spectra'] = outpath_spectra
    spectra_label = spectra_options['spectra_label']


    # Avoid accidental large chunks
    if not subfilter.global_config['no_dask']:
        dask.config.set({"array.slicing.split_large_chunks": True})
        dask_chunks = subfilter.global_config['dask_chunks']

    # Loop over files
    for infile in infiles:

        # Open main dataset with chunking options
        dataset = xr.open_dataset(infile, chunks=dask_chunks)

#        breakpoint()

        # Open reference dataset, if available
        if ref_file is not None:
            ref_dataset = xr.open_dataset(ref_file)
        else:
            ref_dataset = None
        print("Reference Dataset loaded.")

        # Get model resolution values (from input options, file attrs, or options_database)
        dx, dy, options = configure_model_resolution(dataset, options)

        # Work out the file spatial and time dimensions
        [itime, iix, iiy, iiz] = get_string_index(dataset.dims,
                                                  ['time', 'x', 'y', 'z'])
        timevar = list(dataset.dims)[itime]
        xvar = list(dataset.dims)[iix]
        yvar = list(dataset.dims)[iiy]
        zvar = list(dataset.dims)[iiz]
        npoints = dataset.dims[xvar]    # Used to match filter size to grid size

        # Construct derived data (or read existing from file) NOTE: input the override option
        #  Puts basic and additional fields on output grid (opgrid)
        #  File name of output is that of input appended with fname.
        print("setup_derived_data_file...")
        derived_data, exists = \
            sf.setup_derived_data_file( infile, outpath, fname,
                                       options, override = override)
        if exists :
            print('Derived data file already exists.' )
            print("These are the variables in the existing derived dataset:")
            print(derived_data['ds'].variables)

#        breakpoint()

        # Now create list of filter definitions.
        filter_list = list([])
        for i,sigma in enumerate(sigma_list):
            if filter_name == 'gaussian':
                filter_id = 'gn_{:05n}'.format(sigma)
                twod_filter = filt.Filter(filter_id,
                                          filter_name,
                                          npoints=npoints,
                                          sigma=sigma,
                                          delta_x=dx)
            elif filter_name == 'wave_cutoff':
                filter_id = 'wc_{:05n}'.format(sigma)
                twod_filter = filt.Filter(filter_id,
                                          filter_name,
                                          npoints=npoints,
                                          wavenumber=np.pi/(2*sigma),
                                          delta_x=dx)
            elif filter_name == 'circular_wave_cutoff':
                filter_id = 'cwc_{:05n}'.format(sigma)
                twod_filter = filt.Filter(filter_id,
                                          filter_name,
                                          npoints=npoints,
                                          wavenumber=np.pi/(2*sigma),
                                          delta_x=dx)
            elif filter_name == 'running_mean':
                filter_id = 'rm_{:05n}'.format(sigma)
                # sigma used to construct running_mean width
                width = int(np.round( sigma/dx *  np.sqrt(2.0 *np.pi))+1)
                twod_filter = filt.Filter(filter_id,
                                          filter_name,
                                          npoints=npoints,
                                          width=width,
                                          delta_x=dx)

            print(twod_filter)
            filter_list.append(twod_filter)

        # Add whole domain filter  (IF DESIRED)
        if (include_domain_mean):
            filter_name = 'domain'
            filter_id = 'DomainMean'
            twod_filter = filt.filter_2d(filter_id, filter_name, delta_x=dx)
            filter_list.append(twod_filter)

        print(filter_list)


        # Loop over list of filters
        for twod_filter in filter_list:

            print("Processing using filter: ")
            print(twod_filter)

            # Creates the dataset to contain the filtered variables.
            filtered_data, exists = \
                sf.setup_filtered_data_file( infile, outpath, fname,
                                           options, twod_filter, override=override)

            # Check status of existing file of same name
            if exists :
                print('Filtered data file already exists.' )
                print("These are the variables in the existing filtered dataset: ")
                print(filtered_data['ds'].variables)
            else :
                # Perform actual filtering step.
                print("\n    filter_variable_list...", twod_filter.id)
                field_list =sf.filter_variable_list(dataset, ref_dataset,
                                                    derived_data, filtered_data,
                                                    options, twod_filter,
                                                    var_list=var_list,
                                                    grid = opgrid)


                # This script is designed to only compute filtered variable pairs and 
                #   deformation fields for filters with a sigma_list entry that matches
                #   a germano_list entry and when their switches are enabled by the 
                #   yaml options.  Here we check the sigma/germano condition.
                if twod_filter.attributes['sigma'] in germano_list:
                    # Creates filtered versions of paired input variables ( yields L-terms )
                    print("\n    filter_variable_pair_list...",twod_filter.id)
                    if run_quad_fields:
                        quad_field_list =sf.filter_variable_pair_list(dataset,
                                                        ref_dataset,
                                                        derived_data, filtered_data,
                                                        options,
                                                        twod_filter, var_list=var_pair_list,
                                                        grid = opgrid)
                    # Creates filtered versions of the deformation field 
                    print("\n    filtered_deformation...",twod_filter.id)
                    if run_deformation_fields:
                        deformation_r, deformation_s = sf.filtered_deformation(
                                                        dataset,
                                                        ref_dataset,
                                                        derived_data,
                                                        filtered_data, options,
                                                        twod_filter, grid = opgrid)
                        deformation_r = 0
                        deformation_s = 0
                # end if check on germano_list

            # Done with filtering for present sigma. Close dataset.
            filtered_data['ds'].close()

            # ======== Next Stage: spectra (filtered) ======================================
            # Prepare data for spectral analysis, open filtered data.
            filtered_data['ds'] = xr.open_dataset(filtered_data['file'],
                                                  chunks=dask_chunks)
            ds = filtered_data['ds']

            # Get model resolution values (from input options, file attrs, or options_database)
            dx, dy, spectra_options = configure_model_resolution(ds, spectra_options)

            spectra_filt_ds, exists = setup_child_file(filtered_data['file'],
                                                       outpath_spectra, spectra_label,
                                                       spectra_options, override=override)

            print("Working on file: "+spectra_filt_ds['ds'])
            dso = spectra.spectra_variable_list(ds, spectra_filt_ds,
                                                spectra_options,
                                                var_list=spectra_options['filtered_spectra_fields'])
            filtered_data['ds'].close()

        # endfor twod_filter loop over filter_list

        # Close input data
        dataset.close()
        # Close derived data
        derived_data['ds'].close()

        # ======== Next Stage: spectra (derived)============================================
        # Re-open/chunk derived data
        derived_data['ds'] = xr.open_dataset(derived_data['file'],
                                             chunks=dask_chunks)
        ds = derived_data['ds']

        # Get model resolution values (from input options, file attrs, or options_database)
        dx, dy, spectra_options = configure_model_resolution(ds, spectra_options)

        spectra_derived_ds, exists = setup_child_file(derived_data['file'], outpath_spectra,
                             spectra_label, spectra_options, override=override)

        dso = spectra.spectra_variable_list(ds, spectra_derived_ds, spectra_options,
                                            var_list=spectra_options['derived_spectra_fields'])

        print('--------------------------------------')

        print(derived_data)

        print('--------------------------------------')
        derived_data['ds'].close()

    # end infile loop over infiles

    print("IAMDONE")  # done tag for stdout tracking

if __name__ == "__main__":
    main()

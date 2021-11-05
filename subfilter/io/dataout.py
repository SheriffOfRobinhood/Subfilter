# -*- coding: utf-8 -*-
"""
Created on Mon Aug  2 12:02:20 2021.

@author: paclk
"""
import time
import os
import xarray as xr
import config
from dask.diagnostics import ProgressBar
from subfilter import executing_on_cluster


def save_field(dataset, field, write_to_file=True):
    """
    Save dask-chunked xarray field to xarray Dataset

    Parameters
    ----------
    dataset : xarray Dataset
        Output dataset.
    field :  dask-chunked xarray
        Input field.
    write_to_file : bool, optional
        DESCRIPTION. The default is True.

    Returns
    -------
    None.

    """
    fn = dataset['file'].split('/')[-1]
    if field.name not in dataset['ds']:
        print(f"Saving {field.name} to {fn}.")
        dataset['ds'][field.name] = field
        encoding = {field.name: {"dtype": "float32"} }
        if write_to_file:
            d = dataset['ds'][field.name].to_netcdf(
                    dataset['file'],
                    unlimited_dims="time",
                    mode='a', compute=False, encoding = encoding)
            # Toggle ProgressBar depending on compute space
            #   (for cleaner stdout)
            if executing_on_cluster:
                results = d.compute()
            else: 
                with ProgressBar():
                    results = d.compute()
            # This wait seems to be needed to give i/o time to flush caches.
            time.sleep(config.subfilter_setup['write_sleeptime'])
    else:
        print(f"{field.name} already in {fn}")
#    print(dataset['ds'])
    return dataset['ds'][field.name]

def setup_child_file(source_file, destdir, outtag, options, override=False) :
    """
    Create NetCDF dataset for derived data in destdir.

    File name is original file name concatenated with filter_def.id.

    Parameters
    ----------
    source_file     : str
        Input NetCDF file name.
    destdir         : str
        Output directory.
    options         : dict
        Options dictionary
    override=False  : bool
        if True force creation of file

    Returns
    -------
    do              : dict
        {**'file'**: derived_dataset_name (str) - file name,\n
         |**'ds'**: derived_dataset (xarray Dataset) - NetCDF dataset for derived data}
    exists          : bool
        True when input **source_file** already existed and was not overwritten

    @author: Peter Clark
    """


    if os.path.isfile(source_file):
        ds = xr.open_dataset(source_file)
        atts = ds.attrs
    else:
        raise FileNotFoundError(f"Cannot find file {source_file}.")

    derived_dataset_name = os.path.basename(source_file)
    derived_dataset_name = ('.').join(derived_dataset_name.split('.')[:-1])
    if options['l_slurm_job_tag']:
        jn = os.environ['SLURM_JOB_NAME']
        derived_dataset_name = destdir+derived_dataset_name + "_"+jn+ "_" + outtag + ".nc"
    else:
        derived_dataset_name = destdir+derived_dataset_name + "_" + outtag + ".nc"

    exists = os.path.isfile(derived_dataset_name)

    if exists and not override :

        derived_dataset = xr.open_dataset(derived_dataset_name)

    else :
        if exists:
            print(f"Overwriting file {derived_dataset_name}.")
        exists = False

        # Modify the True/False options for writing.
        for inc in options:
            if options[inc] is True:
                options[inc]=1
            if options[inc] is False:
                options[inc]=0

        derived_dataset = xr.Dataset(coords =
                        {'z':ds.coords['z'],'zn':ds.coords['zn']})
        derived_dataset.attrs = {**atts, ** options}
        derived_dataset.to_netcdf(derived_dataset_name, mode='w')
        print(derived_dataset)
    do = {'file':derived_dataset_name, 'ds': derived_dataset}
    ds.close()
    return do, exists


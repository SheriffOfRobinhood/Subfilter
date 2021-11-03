# -*- coding: utf-8 -*-
"""
Created on Tue Oct 23 11:27:25 2018

@author: Peter Clark
"""
import os
#import netCDF4
#from netCDF4 import Dataset
import numpy as np
#import pandas as pd
import xarray as xr
import matplotlib
import matplotlib.pyplot as plt
import dask
#from dask.distributed import Client

#from dask.diagnostics import Profiler, ResourceProfiler, CacheProfiler

import subfilter.subfilter as sf
import subfilter.filters as filt
import subfilter.utils.deformation as defm
import subfilter.utils.difference_ops as do
import importlib
importlib.reload(do)

from subfilter.utils.dask_utils import re_chunk


from subfilter.utils.string_utils import get_string_index
from subfilter.io.dataout import save_field
from subfilter.io.datain import get_data
from subfilter.io.MONC_utils import options_database
from subfilter import subfilter_setup

test_case = 1
run_quad_fields = False
run_deformation_fields = True
override = True

plot_type = '.png'
figshow = True

def plot_deformation(var_name, deform, ilev, iy, plot_dir, plot_type):
    for it, time in enumerate(deform.coords['time']):
        fig, ax = plt.subplots(3,3,figsize=(20,15))
        for i in range(0,3):
            for j in range(0,3):
                deform.isel(time=0, z=ilev, i=i, j=j).plot.imshow(ax=ax[i,j])
                ax[i,j].set_title(f'd({i}{j})')
        fig.tight_layout()

        plt.savefig(plot_dir+var_name+f'_prof_{it:02d}'+plot_type)
#        plt.close()
        plt.show()


if True:
#def main():
    '''
    Top level code, a bit of a mess.
    '''

    if test_case == 0:
        config_file = 'config_test_case_0.yaml'
        dir = 'C:/Users/paclk/OneDrive - University of Reading/ug_project_data/Data/'
        odir = 'C:/Users/paclk/OneDrive - University of Reading/ug_project_data/Data/'
        file = 'diagnostics_3d_ts_21600.nc'
        ref_file = 'diagnostics_ts_21600.nc'
    elif test_case == 1:
        config_file = 'config_test_case_1.yaml'
        dir = 'C:/Users/paclk/OneDrive - University of Reading/traj_data/CBL/'
        odir = 'C:/Users/paclk/OneDrive - University of Reading/traj_data/CBL/'
        file = 'diagnostics_3d_ts_13200.nc'
        ref_file = None

#dir = 'C:/Users/paclk/OneDrive - University of Reading/Git/python/Subfilter/test_data/BOMEX/'
#odir = 'C:/Users/paclk/OneDrive - University of Reading/Git/python/Subfilter/test_data/BOMEX/'


#file = 'diagnostics_ts_18000.0.nc'
#ref_file = 'diagnostics_ts_18000.0.nc'

    options, update_config = sf.subfilter_options(config_file)

    odir = odir + 'test_dask_' + options['FFT_type']+'/'
    os.makedirs(odir, exist_ok = True)

    plot_dir = odir + 'plots/'
    os.makedirs(plot_dir, exist_ok = True)

#    client = Client()
#    client
#    dask.config.set(scheduler='threads')
#   Non-global variables that are set once

    dask.config.set({"array.slicing.split_large_chunks": True})

    dataset = xr.open_dataset(dir+file, chunks={'z':'auto', 'zn':'auto'})

    print(dataset)

    if ref_file is not None:
        ref_dataset = xr.open_dataset(dir+ref_file)
    else:
        ref_dataset = None

    od = options_database(dataset)
    if od is None:
        dx = options['dx']
        dy = options['dy']
    else:
        dx = float(od['dxx'])
        dy = float(od['dyy'])

    [itime, iix, iiy, iiz] = get_string_index(dataset.dims, ['time', 'x', 'y', 'z'])
    timevar = list(dataset.dims)[itime]
    xvar = list(dataset.dims)[iix]
    yvar = list(dataset.dims)[iiy]
    zvar = list(dataset.dims)[iiz]
    npoints = dataset.dims[xvar]

# For plotting
    ilev = 15
    iy = 40

    opgrid = 'w'
    fname = 'test_rewrite'

    derived_data, exists = \
        sf.setup_derived_data_file( dir+file, odir, fname,
                                   options, override=override)
    if exists :
        print('Derived data file exists' )
        print("Variables in derived dataset.")
        print(derived_data['ds'].variables)

    if run_deformation_fields:

        deform = defm.deformation(dataset, ref_dataset, derived_data,
            options, grid='w')

        derived_data['ds'].close()

        derived_data['ds'] = xr.open_dataset(derived_data['file'],
                                             chunks={'z':'auto', 'zn':'auto'})

        deform = derived_data['ds']['deformation']
        plot_deformation('deformation', deform, ilev, iy, plot_dir, plot_type)

    else:
        z = dataset["z"]
        zn = dataset["zn"]

        u = get_data(dataset, ref_dataset, "u", options)
        [iix, iiy, iiz] = get_string_index(u.dims, ['x', 'y', 'z'])

        sh = np.shape(u)

        max_ch = subfilter_setup['chunk_size']

        nch = int(sh[iix]/(2**int(np.log(sh[iix]*sh[iiy]*sh[iiz]/max_ch)/np.log(2)/2)))

        print(f'nch={nch}')
    #    nch = 32

        u = re_chunk(u, xch=nch, ych=nch, zch = 'all')
    #    u = re_chunk(u, xch='all')
    #    ux = do.d_by_dx_field_on_u(u, z, grid = 'p' )
        ux = do.d_by_dz_field(u, z, zn, grid = 'p' )

        print(ux)

        plt.figure(1)
        ux.isel(time=0, zn=ilev).plot.imshow()
        plt.show()

        plt.figure(2)
        ux.isel(time=0, zn=ilev, y_p=iy).plot()
        plt.show()

        ux = 0

    print('--------------------------------------')

    print(derived_data)

    print('--------------------------------------')
    derived_data['ds'].close()
    dataset.close()


#if __name__ == "__main__":
#    main()

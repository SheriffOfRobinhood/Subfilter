# -*- coding: utf-8 -*-
"""

  plot_dataset.py

Created on Tue Oct 23 09:52:58 2018

@author: Peter Clark
"""
import os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import xarray as xr



#dir = 'C:/Users/paclk/OneDrive - University of Reading/python/SG/'
#file = 'diagnostics_ts_14400.0.nc'
#dir = 'C:/Users/paclk/OneDrive - University of Reading/python/subfilter/test_data/BOMEX/'
#file = 'diagnostics_ts_18000.0.nc'
#file = 'diagnostics_ts_18000.0_filter_00.nc'

#dir = '/gws/nopw/j04/paracon_rdg/users/toddj/MONC_RCE_1.5/MONC_RCE_1.5_m1600_g0084/diagnostic_files/'
#dir = 'C:/Users/paclk/OneDrive - University of Reading/Git/python/Subfilter/test_data/BOMEX/'
#dir = 'C:/Users/paclk/OneDrive - University of Reading/ug_project_data/Data/'
#dir = 'C:/Users/paclk/OneDrive - University of Reading/ug_project_data/Data/test_dask_RFFT/'
dir = 'C:/Users\paclk\OneDrive - University of Reading/traj_data/CBL/test_dask_RFFT/'
#file = 'diagnostics_ts_4784400.0.nc'
#file = 'diagnostics_ts_18000.0.nc'
#file = 'diagnostics_3d_ts_21600.nc'
#file = 'diagnostics_3d_ts_21600_test_dask.nc'
#file = 'diagnostics_3d_ts_21600_test_dask_filter_ga00.nc'
file = 'diagnostics_3d_ts_13200_test_dask_filter_ga00.nc'

dataset = xr.open_dataset(dir+file,chunks={'time':1,'i':1,'j':1,
                                                'x_p':160, 'y_p':160,
                                                'z':'auto'})
print(dataset)
fig, axes = plt.subplots(ncols=2, figsize=(12,6))
d_r = dataset['deformation_r']
print(d_r)
d_r.isel(i=0, j=0, time=0, z=10).plot(ax=axes[0])
d_s = dataset['deformation_s']
print(d_s)
d_s.isel(i=0, j=0, time=0, z=10).plot(ax=axes[1])
plt.show()
dataset.close()


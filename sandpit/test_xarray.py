# -*- coding: utf-8 -*-
"""
Created on Tue Mar 16 09:34:26 2021

@author: paclk
"""

import os
import sys
#import netCDF4

#from netCDF4 import Dataset


import numpy as np
import pandas as pd
import xarray as xr
import subfilter as sf

#dir = 'C:/Users/paclk/OneDrive - University of Reading/Git/python/Subfilter/test_data/BOMEX/'
dir = 'C:/Users/paclk/OneDrive - University of Reading/ug_project_data/Data/'
#odir = 'C:/Users/paclk/OneDrive - University of Reading/Git/python/Subfilter/test_data/BOMEX/'
odir = 'C:/Users/paclk/OneDrive - University of Reading/ug_project_data/Data/'
odir = odir + 'test_xarray/'

os.makedirs(odir, exist_ok = True)

#file = 'diagnostics_ts_18000.0.nc'
file = 'diagnostics_3d_ts_21600.nc'

ofile ='diagnostics_3d_ts_21600_test_xarray.nc'

#ref_file = 'diagnostics_ts_18000.0.nc'


dataset = xr.open_dataset(dir+file)

print(dataset)

print(dataset.coords)

od = sf.options_database(dataset)
#print(od)

derived_dataset = xr.Dataset(coords=dataset.coords)

for var in dataset.variables:
    if var in dataset.coords:
        print(var, " is a coord.")
    else:
        print(var)
        d = dataset[var]

#        print(d)
        if 'x' in d.dims:
            nx = d.shape[d.dims.index('x')]
            ny = d.shape[d.dims.index('y')]

            if var == 'u' :
                x = (np.arange(nx) + 0.5) * np.float64(od['dxx'])
                xn = 'x_u'
            else:
                x = np.arange(nx) * np.float64(od['dxx'])
                xn = 'x_p'

            if var == 'v' :
                y = (np.arange(ny) + 0.5) * np.float64(od['dyy'])
                yn = 'y_v'
            else:
                y = np.arange(ny) * np.float64(od['dyy'])
                yn = 'y_p'


            dr = d.rename({'x':xn, 'y':yn})

            dr.coords[xn] = x
            dr.coords[yn] = y
#            print(dr)
            derived_dataset[var] = dr

#print(dataset)
print(derived_dataset)

# print(derived_dataset)

# derived_dataset.to_netcdf(odir+ofile)

# for var in dataset.variables:
#     if var in dataset.coords:
#         print(var, " is a coord.")
#     else:
#         print(var)
#         d = dataset[var]

#         if 'x' in d.dims:


#             nx = d.shape[d.dims.index('x')]

#             x_t = np.arange(nx) * np.float64(od['dxx'])

#             dr = d.rename({'x':'x_t'})

#             dr.coords['x_t'] = x_t

#             derived_dataset[var+'_m'] = dr
#     #        dat = dataset[var].values

#             derived_dataset.to_netcdf(odir+ofile)

# print(derived_dataset)

# print(derived_dataset['u_m'])

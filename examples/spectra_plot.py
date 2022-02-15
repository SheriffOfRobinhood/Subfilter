# -*- coding: utf-8 -*-
"""
Created on Wed Jul 28 17:44:16 2021

@author: paclk
"""
import os

import numpy as np
import xarray as xr

import matplotlib.pyplot as plt


# dir = 'C:/Users/paclk/OneDrive - University of Reading/ug_project_data/Data/'
# file = 'diagnostics_3d_ts_21600.nc'
# dx = 50.0
# dy = 50.0


dir = 'C:/Users/paclk/OneDrive - University of Reading/traj_data/CBL/'
file = 'diagnostics_3d_ts_13200.nc'
dx = 5.0
dy = 5.0

# Set up outfile
outdir = os.path.join(dir, 'spectra/')
os.makedirs(outdir, exist_ok = True)  # make outdir if it doesn't already exist
outtag = "spectra_w_2D"
outfile = os.path.join(outdir,('.').join(os.path.basename(file).split('.')[:-1]) + "_"+outtag+".nc")


dso = xr.open_dataset(outfile)
k = dso['hfreq']
kx = dso['xfreq']
ky = dso['yfreq']
kE_k = k * dso['spec_2d_w']
kE_kx = kx * dso['spec_xdir_w']
kE_ky = ky * dso['spec_ydir_w']

#kE_k = k * dso['spec_2d_u']


kE_kp = kE_k.mean(dim='time').sel(z=600, method='nearest')
kE_kp.plot(xscale='log', yscale='log')

kE_kxp = kE_kx.mean(dim='time').sel(z=600, method='nearest')
kE_kxp.plot(xscale='log', yscale='log')
kE_kyp = kE_ky.mean(dim='time').sel(z=600, method='nearest')
kE_kyp.plot(xscale='log', yscale='log')

#plt.ylim([0.0001,1])
#plt.xlim([0.0001,0.1])
plt.ylabel('kE(k)')
plt.tight_layout()
ax=plt.gca()

xr = ax.get_xlim()
yr = ax.get_ylim()

#k = kx
#kE_kp = kE_kxp

i = kE_kp.argmax()
if i.size > 1:
    i = i[0]

ymax = kE_kp[i]
x_ymax = k[i]

yidl = lambda x : ymax.values * (x / x_ymax.values)**(-2/3)

plt.plot(xr, yidl(xr))


dso.close()



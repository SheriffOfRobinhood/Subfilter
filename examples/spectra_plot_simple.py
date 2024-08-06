# -*- coding: utf-8 -*-
"""
Created on Wed Jul 28 17:44:16 2021

@author: paclk
"""
import os

import numpy as np
import xarray as xr

import matplotlib.pyplot as plt

def inv(k):
    return 1.0 / k

dirroot = 'E:/Data/'
outdir = dirroot + 'CBL/'
outdir = outdir + 'test_filter_spectra/'

file = 'diagnostics_3d_ts_14400.nc'
dx = 5.0
dy = 5.0

plot_height = 600

# Set up outfile
outtag = "spec"
outfile = os.path.join(outdir,('.').join(os.path.basename(file).split('.')[:-1]) + "_"+outtag+".nc")


dso = xr.open_dataset(outfile)

k = dso['hfreq']
kE_k = k * dso['spec_2d_w']


kE_kp = kE_k.mean(dim='time').sel(z=plot_height, method='nearest')
kE_kp.plot(xscale='log', yscale='log', figsize=(8,6))

#plt.ylim([0.0001,1])
#plt.xlim([0.0001,0.1])
plt.ylabel('kE(k)')
ax=plt.gca()

xr = ax.get_xlim()
yr = ax.get_ylim()

i = kE_kp.argmax()
if i.size > 1:
    i = i[0]

ymax = kE_kp[i]
x_ymax = k[i]

yidl = lambda x : ymax.values * (x / x_ymax.values)**(-2/3)

plt.plot(xr, yidl(xr))

secax = ax.secondary_xaxis('top', functions=(inv, inv))
secax.set_xlabel('wavelength (m)')

dso.close()
plt.tight_layout()

plt.savefig(outdir+'spec_2d_w.png')

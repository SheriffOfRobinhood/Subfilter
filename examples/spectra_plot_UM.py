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

def plt_ISR_spec(ax, k, kE_kp):
    xran = ax.get_xlim()
    yran = ax.get_ylim()
    
    ind = (kE_kp*(k**(2/3))).argmax()
    if ind.size > 1:
        ind = ind[0]
    
    ymax = kE_kp[ind]
    x_ymax = k[ind]
    
    yidl = lambda x : ymax.values * (x / x_ymax.values)**(-2/3)

    ax.plot(xran, yidl(xran),'k--')
    
    return


dirroot = 'E:/Data/'
outdir = dirroot + 'UM_CBL/'
outdir = outdir + 'test_filter_spectra/'
outtag = "spec"

runs = [
        '100m',
        '200m',
        '250m',
        '500m',
        '1km',
#        '10km',
        ]


plot_height = 500

fig, axes = plt.subplots(2,2, figsize=(12,12))

for i in range(4):

    ax = axes[i // 2, i % 2]
    ax.set_xlim([1E-5,1E-2])
    ax.set_ylim([0.0001,0.1])
    
    secax = ax.secondary_xaxis('top', functions=(inv, inv))
    secax.set_xlabel('wavelength (m)')

for run in runs:
    file = f'dc455_{run}_L100a_pr000.nc'
    
    if 'k' in run:
        ind = run.index('k')
        dx = dy = float(run[0:ind])*1000
    else:
        ind = run.index('m')
        dx = dy = float(run[0:ind])
        

    # Set up outfile
    outfile = os.path.join(outdir,('.').join(os.path.basename(file).split('.')[:-1]) + "_"+outtag+".nc")
        
    dso = xr.open_dataset(outfile)

    dso = dso.where(dso.elapsed_time >= 2, drop=True)
    
    total = 0
    k = dso['hfreq']
    for i, var in enumerate(['u', 'v', 'w']):
        
        ax = axes[i // 2, i % 2]
                
        # print([i // 2, i % 2])
        
        E_k = dso[f'spec_2d_{var}_on_p']

        kE_k = k * E_k
        kE_kp = kE_k.mean(dim='time').sel(z_p=plot_height, method='nearest')
        
        # print(kE_kp.max())
        
        kE_kp.plot(xscale='log', yscale='log', ax=ax, label=run)
        
        if run == '100m': plt_ISR_spec(ax, k, kE_kp)
        
        ax.set_title(var)
        
        total = total + kE_kp / 2
        
    i = 3

    ax = axes[i // 2, i % 2]

    total.plot(xscale='log', yscale='log', ax=ax, label=run)
    # total.plot(xscale='log', ax=ax, label=run)
    if run == '100m': plt_ISR_spec(ax, k, total)
    ax.set_title('KE')

    dso.close()

for i in range(4):

    ax = axes[i // 2, i % 2]
    ax.legend(loc='upper left', fontsize=10)    
    ax.set_ylabel('kE(k)')

plt.tight_layout()

plt.savefig(outdir+'spec_2d.png')
plt.show()

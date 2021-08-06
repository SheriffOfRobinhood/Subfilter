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


# dir = 'C:/Users/paclk/OneDrive - University of Reading/ug_project_data/Data/'
# file = 'diagnostics_3d_ts_21600.nc'
# dx = 50.0
# dy = 50.0


dir = 'C:/Users/paclk/OneDrive - University of Reading/traj_data/CBL/test_filter_spectra_RFFT/'
fileroot = 'diagnostics_3d_ts_*'

fname = '_filter_spectra'
outtag = '_spectra_w_2D'

plot_s = False


files = dir+fileroot+fname+outtag+'.nc'

dx = 5.0
dy = 5.0

plot_height = 600

# Set up outfile
#outdir = os.path.join(dir, 'spectra/')
#os.makedirs(outdir, exist_ok = True)  # make outdir if it doesn't already exist
#outfile = os.path.join(outdir,('.').join(os.path.basename(file).split('.')[:-1]) + "_"+outtag+".nc")


#dso = xr.open_dataset(file)
dso = xr.open_mfdataset(files)
print(dso)
k_ref = dso['hfreq']
kE_k = k_ref * dso['spec_2d_w_on_w']
kE_kp_ref = kE_k.mean(dim='time').sel(z=plot_height, method='nearest')

fig1, axa = plt.subplots(2,1,figsize=(8,12))

kE_kp_ref.plot(xscale='log', yscale='log', label='Ref 5 m', ax=axa[0])

dso.close()
nfilt=6
for filt in range(nfilt):
    filt_files = dir+fileroot+fname+f'_filter_ga{filt:02d}'+outtag+'.nc'
#    dso = xr.open_dataset(filt_file)
    dso = xr.open_mfdataset(filt_files)

    print(dso)
    k = dso['hfreq']
    kE_k = k * dso['spec_2d_f(w_on_w)_r']
    kE_kp = kE_k.mean(dim='time').sel(z=plot_height, method='nearest')
    kE_kp.plot(xscale='log', yscale='log',
               label=f'{dso.attrs["sigma"]:2.0f} m'+r'$^r$', ax=axa[0])
    if plot_s:
        kE_k = k * dso['spec_2d_f(w_on_w)_s']
        kE_kp = kE_k.mean(dim='time').sel(z=plot_height, method='nearest')
        kE_kp.plot(xscale='log', yscale='log',
                   label=f'{dso.attrs["sigma"]:2.0f} m'+r'$^s$', ax=axa[0])
    dso.close()


axa[0].set_ylim([0.0001,0.1])
#plt.xlim([0.0001,0.1])
axa[0].set_ylabel('kE(k)')
axa[0].legend()

#k = kx
#kE_kp = kE_kxp
#ax=plt.gca()
secax = axa[0].secondary_xaxis('top', functions=(inv, inv))
secax.set_xlabel('wavelength (m)')

xrn = axa[0].get_xlim()
yrn = axa[0].get_ylim()


i = kE_kp_ref.values.argmax()
if i.size > 1:
    i = i[0]

ymax = kE_kp_ref[i]
x_ymax = k_ref[i]
yidl = lambda x : ymax.values * (x / x_ymax.values)**(-2/3)

axa[0].plot(xrn, yidl(xrn))

kE_kp_idl = kE_kp_ref.copy()
k1 = np.where(k_ref > 1/75.)[0][0]

kE_kp_idl[k1:] = kE_kp_ref[k1] * (k[k1:]/k[k1])**(-2/3)

kE_kp_idl.plot( label = 'Ideal', ax=axa[0])
axa[0].set_ylabel(r'$kE(k)$')

plt.tight_layout()

idl_filt = kE_kp_ref / kE_kp_idl
idl_filt.plot(xscale='log', yscale='log', ax=axa[1], label='Ref 5 m')

for filt in range(nfilt):
    filt_file = dir+fileroot+fname+f'_filter_ga{filt:02d}'+outtag+'.nc'
    dso = xr.open_dataset(filt_file)

    print(dso)
    k = dso['hfreq']
    kE_k = k * dso['spec_2d_f(w_on_w)_r']
    kE_kp_rat  = kE_k.mean(dim='time').sel(z=plot_height, method='nearest')\
        / kE_kp_idl
    kE_kp_rat.plot(xscale='log', yscale='log',
               label=f'{dso.attrs["sigma"]:2.0f} m'+r'$^r$', ax=axa[1])
    dso.close()


alpha = 2.0
sigma_5 = 5*2/np.pi
sigma_5 = 3.0

G_5 =np.exp(-(2*np.pi*sigma_5*k_ref.values)**alpha)
axa[1].plot(k_ref.values, G_5, label=r'$\alpha$=2, $\sigma$= 3.0 m ')#10/$\pi$

alpha = 2.0
sigma_5 = 2.3

G_5 =np.exp(-(2*np.pi*sigma_5*k_ref.values)**alpha)
axa[1].plot(k_ref.values, G_5, label=r'$\alpha$=2, $\sigma$=2.3 m ')

alpha = 2.4
sigma_5 = 2.8 #5*2/np.pi*0.8
G_5 =np.exp(-(2*np.pi*sigma_5*k_ref.values)**alpha)
axa[1].plot(k_ref.values, G_5, label=r'$\alpha$=2.4, $\sigma$=2.5 m')

secax = axa[1].secondary_xaxis('top', functions=(inv, inv))
secax.set_xlabel('wavelength (m)')

axa[1].legend()
axa[1].set_ylim([0.1,1.1])
axa[0].set_ylabel(r'$G(k)\times G(k)^*$')

plt.tight_layout()
plt.savefig(dir+'Gaussian_filter.png')

dso.close()



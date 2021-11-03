# -*- coding: utf-8 -*-
"""
Created on Tue Sep 15 16:10:18 2020

@author: paclk
"""
import numpy as np
import xarray as xr
import matplotlib
import matplotlib.pyplot as plt
import time

import subfilter.subfilter as sf
import subfilter.filters as filt
#import difference_ops as do

options = {
#       'FFT_type': 'FFTconvolve',
#       'FFT_type': 'FFT',
       'FFT_type': 'RFFT',
#       'FFT_type': 'DIRECT',
          }


#sigma_list = [2000.0, 1000.0, 500.0, 200.0, 100.0]
#sigma_list = [71.0]


#sigma_list = [2000.0, 1000.0, 500.0, 200.0, 100.0]
sigma_list = [1600.0, 1200, 800.0, 600, 400.0, 200.0]
dx = 100.0

#filter_name = 'gaussian'
filter_name = 'wave_cutoff'
#filter_name = 'circular_wave_cutoff'
#filter_name = 'running_mean'
#filter_name = 'one_two_one'
filter_list = list([])

N = 128
dx = 100

x = np.linspace(0,(N-1)*dx,N )
y = x.copy()
xc, yc = np.meshgrid(x, y)

Lx = (N)*dx/3
Ly = (N)*dx/10

kx = 2.0 * np.pi / Lx
ky = 2.0 * np.pi / Ly

field = np.cos(kx * xc)*np.cos(ky * yc)+ np.cos(np.pi / dx * xc)*np.cos(np.pi / dx * yc)
var = xr.DataArray(field, name='test', coords=[x, y], dims=['x_p', 'y_p'])
print(var)

levs=np.linspace(-1,1,41)

xc /= 1000
yc /= 1000

x /= 1000
y /= 1000


for i,sigma in enumerate(sigma_list):
    filter_id = 'filter_{:02d}'.format(i)
    if filter_name == 'running_mean':
        width = min(N-2,int(np.round( sigma/dx * np.pi *2.0/3.0)+1))
    else:
        width = min(N, int(np.ceil( sigma/dx *4)*2+1))


 #   width = -1
    twod_filter = filt.Filter(filter_id,
                                 filter_name, wavenumber=2*np.pi/(sigma),
                                 sigma=sigma, width=width, npoints = N,
                                 delta_x=dx, set_fft=True)

    print(twod_filter)
    filter_list.append(twod_filter)


print(filter_list)

for twod_filter in filter_list:

#    print(np.sum(twod_filter.data))

#    f = plt.figure()
    sigma=twod_filter.attributes['sigma']
    nx, ny = np.shape(twod_filter.data)

    x = np.linspace(-(nx//2), (nx-1)//2, nx) * twod_filter.attributes['delta_x']/1000

    f, ax = plt.subplots(2,2,figsize=(12,12))

    c = ax[0,0].contour(x, x, twod_filter.data,20,colors='blue')
    ax[0,0].set_xlabel("x/km")
    ax[0,0].set_ylabel("y/km")

#    ax[0,0].legend()


    p1 = ax[0,1].plot(x, twod_filter.data[nx//2, :],label='Sampled')
    ax[0,1].set_xlabel("x/km")
    ax[0,1].legend()
    p2 = ax[1,0].plot(x, twod_filter.data[:, ny//2],label='Sampled')
    ax[1,0].set_xlabel("x/km")
    ax[1,0].legend()
    p3 = ax[1,1].plot(x*np.sqrt(2),twod_filter.data.diagonal(),label='Sampled')
    ax[1,1].set_xlabel("x/km")
    ax[1,1].legend()
    f.tight_layout()

    time_1 = time.time()
    (var_r, var_s) = sf.filtered_field_calc(var, options, twod_filter )
    time_2 = time.time()
    elapsed_time = time_2 - time_1
    print(f'Elapsed time = {elapsed_time}')
    # var_r = var_r['data']
    # var_s = var_s['data']

    E_var = np.mean(var * var)
    E_var_r = np.mean(var_r * var_r)
    E_var_s = np.mean(var_s * var_s)

    fc, axc = plt.subplots(3,3,figsize=(12,12))


    c_f = axc[0,0].contourf(xc, yc, var, levs)
    axc[0,0].set_xlabel("x/km")
    axc[0,0].set_ylabel("y/km")
    axc[0,0].set_title(r"Field E={:1.4f}".format(E_var.data))

    c_r = axc[1,0].contourf(xc, yc, var_r, levs)
    axc[1,0].set_xlabel("x/km")
    axc[1,0].set_ylabel("y/km")
    axc[1,0].set_title(r"Field$^r$ E={:1.4f}".format(E_var_r.data))

    c_s = axc[2,0].contourf(xc, yc, var_s, levs)
    axc[2,0].set_xlabel("x/km")
    axc[2,0].set_ylabel("y/km")
    axc[2,0].set_title(r"Field$^s$ E={:1.4f}".format(E_var_s.data))

    px_f = axc[0,1].plot(y, var[N//2,:])
    axc[0,1].set_title(r"$1/(k_x\sigma)$={:3.2f}".format(1/(kx*sigma)))
    axc[0,1].set_xlabel("x/km")
    px_r = axc[1,1].plot(y, var_r[N//2,:])
    axc[1,1].set_xlabel("x/km")
    px_s = axc[2,1].plot(y, var_s[N//2,:])
    axc[2,1].set_xlabel("x/km")

    py_f = axc[0,2].plot(y, var[:, N//2])
    axc[0,2].set_ylabel("y/km")
    axc[0,2].set_title(r"$1/(k_y\sigma)$={:3.2f}".format(1/(ky*twod_filter.attributes['sigma'])))
    py_r = axc[1,2].plot(y, var_r[:, N//2])
    axc[1,0].set_ylabel("y/km")
    py_s = axc[2,2].plot(y, var_s[:, N//2])
    axc[2,2].set_ylabel("y/km")

    fc.tight_layout()

    f.savefig(f'{filter_name}_{sigma:04.0f}_filter.png')
    fc.savefig(f'{filter_name}_{sigma:04.0f}_test_data.png')

plt.show()

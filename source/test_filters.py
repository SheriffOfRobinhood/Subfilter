# -*- coding: utf-8 -*-
"""
Created on Tue Sep 15 16:10:18 2020

@author: paclk
"""
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

import subfilter as sf
import filters as filt
import difference_ops as do

options = {
#       'FFT_type': 'FFTconvolve',
#       'FFT_type': 'FFT',
       'FFT_type': 'RFFT',
          }


sigma_list = [2000.0, 1000.0, 500.0, 200.0, 100.0]
#sigma_list = [2000.0]
dx = 100.0

filter_name = 'gaussian'
#filter_name = 'wave_cutoff'
#filter_name = 'running_mean'
filter_list = list([])

N = 128
dx = 100

cxy = np.linspace(0,(N-1)*dx,N )

xc, yc = np.meshgrid(cxy,cxy)

Lx = (N)*dx/3
Ly = (N)*dx/10

kx = 2.0 * np.pi / Lx
ky = 2.0 * np.pi / Ly

field = np.cos(kx * xc)*np.cos(ky * yc)
levs=np.linspace(-1,1,41)

xc /= 1000
yc /= 1000

cxy /= 1000



for i,sigma in enumerate(sigma_list):
    filter_id = 'filter_{:02d}'.format(i)
    width = min(N, int(np.ceil( sigma/dx *4)*2+1))
 #   width = -1
    twod_filter = filt.filter_2d(filter_id,
                                 filter_name, wavenumber=2*np.pi/sigma,
                                 sigma=sigma, width=width,
                                 delta_x=dx, use_ave=False)
    twod_filter_ave = filt.filter_2d(filter_id,
                                 filter_name, wavenumber=2*np.pi/sigma,
                                 sigma=sigma, width=width,
                                 delta_x=dx, use_ave=True)

    print(twod_filter)
    filter_list.append([twod_filter, twod_filter_ave])


print(filter_list)

for twod_filter, twod_filter_ave in filter_list:

    print(np.sum(twod_filter.data), np.sum(twod_filter_ave.data))

#    f = plt.figure()
    nx, ny = np.shape(twod_filter.data)

    x = np.linspace(-(nx//2), nx//2, nx) * twod_filter.attributes['delta_x']/1000

    nx1, ny1 = np.shape(twod_filter_ave.data)

    x1 = np.linspace(-(nx1//2), nx1//2, nx1) * twod_filter.attributes['delta_x']/1000

    f, ax = plt.subplots(2,2,figsize=(12,12))

    c = ax[0,0].contour(x, x, twod_filter.data,20,colors='blue')
    c1 = ax[0,0].contour(x1, x1, twod_filter_ave.data,20,colors='red')
    ax[0,0].set_xlabel("x/km")
    ax[0,0].set_ylabel("y/km")

#    ax[0,0].legend()


    p1 = ax[0,1].plot(x, twod_filter.data[nx//2, :],label='Sampled')
    p1 = ax[0,1].plot(x1, twod_filter_ave.data[nx1//2, :],label='Averaged')
    ax[0,1].set_xlabel("x/km")
    ax[0,1].legend()
    p2 = ax[1,0].plot(x, twod_filter.data[:, ny//2],label='Sampled')
    p2 = ax[1,0].plot(x1, twod_filter_ave.data[:, ny1//2],label='Averaged')
    ax[1,0].set_xlabel("x/km")
    ax[1,0].legend()
    p3 = ax[1,1].plot(x*np.sqrt(2),twod_filter.data.diagonal(),label='Sampled')
    p3 = ax[1,1].plot(x1*np.sqrt(2),twod_filter_ave.data.diagonal(),label='Averaged')
    ax[1,1].set_xlabel("x/km")
    ax[1,1].legend()
    plt.tight_layout()


    field_r, field_s = sf.filtered_field_calc(field, options, twod_filter, three_d=False )

    E_field = np.mean(field * field)
    E_field_r = np.mean(field_r * field_r)
    E_field_s = np.mean(field_s * field_s)

    fc, axc = plt.subplots(3,3,figsize=(12,12))


    c_f = axc[0,0].contourf(xc, yc, field, levs)
    axc[0,0].set_xlabel("x/km")
    axc[0,0].set_ylabel("y/km")
    axc[0,0].set_title(r"Field E={:1.4f}".format(E_field))

    c_r = axc[1,0].contourf(xc, yc, field_r, levs)
    axc[1,0].set_xlabel("x/km")
    axc[1,0].set_ylabel("y/km")
    axc[1,0].set_title(r"Field$^r$ E={:1.4f}".format(E_field_r))

    c_s = axc[2,0].contourf(xc, yc, field_s, levs)
    axc[2,0].set_xlabel("x/km")
    axc[2,0].set_ylabel("y/km")
    axc[2,0].set_title(r"Field$^s$ E={:1.4f}".format(E_field_s))

    px_f = axc[0,1].plot(cxy, field[N//2,:])
    axc[0,1].set_title(r"$1/(k_x\sigma)$={:3.2f}".format(1/(kx*twod_filter.attributes['sigma'])))
    axc[0,1].set_xlabel("x/km")
    px_r = axc[1,1].plot(cxy, field_r[N//2,:])
    axc[1,1].set_xlabel("x/km")
    px_s = axc[2,1].plot(cxy, field_s[N//2,:])
    axc[2,1].set_xlabel("x/km")

    py_f = axc[0,2].plot(cxy, field[:, N//2])
    axc[0,2].set_ylabel("y/km")
    axc[0,2].set_title(r"$1/(k_y\sigma)$={:3.2f}".format(1/(ky*twod_filter.attributes['sigma'])))
    py_r = axc[1,2].plot(cxy, field_r[:, N//2])
    axc[1,0].set_ylabel("y/km")
    py_s = axc[2,2].plot(cxy, field_s[:, N//2])
    axc[2,2].set_ylabel("y/km")

    plt.tight_layout()


plt.show()

# -*- coding: utf-8 -*-
"""
Created on Tue Oct 31 13:12:36 2023

@author: paclk
"""
import numpy as np
import subfilter.filters as filt

import matplotlib.pyplot as plt


npoints = 1024
dx = 5.0
sigma = 0.5*dx

ndim = 1

filter_name = 'gaussian'
if filter_name == 'gaussian':
    filter_id = f'filter_ga{0:02d}'
    gfilter = filt.Filter(filter_id,
                              filter_name,
                              npoints=npoints,
                              sigma=sigma,
                              delta_x=dx,
                              ndim = ndim)

fig, ax = plt.subplots(1,2)

ax[0].plot(gfilter.data, '-*')

ff = np.fft.fft(gfilter.data)

# ax[1].semilogy(np.fft.fftshift(np.fft.fftfreq(npoints)),
#            np.abs(np.fft.fftshift(ff)), '-o')

# ax[1].loglog(np.fft.fftshift(np.fft.fftfreq(npoints))[npoints//2:],
#            np.abs(np.fft.fftshift(ff))[npoints//2:])
# ax[1].loglog(np.fft.fftshift(np.fft.fftfreq(npoints))[npoints//2:],
#            np.abs(np.fft.fftshift(ff))[npoints//2:]**2)
ax[1].plot(np.fft.fftshift(np.fft.fftfreq(npoints))[npoints//2:],
           np.abs(np.fft.fftshift(ff))[npoints//2:])
ax[1].plot(np.fft.fftshift(np.fft.fftfreq(npoints))[npoints//2:],
           np.abs(np.fft.fftshift(ff))[npoints//2:]**2)

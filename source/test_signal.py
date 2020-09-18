# -*- coding: utf-8 -*-
"""
Created on Wed Sep 16 10:39:06 2020

@author: paclk
"""
import numpy as np
from scipy import signal
import matplotlib.pyplot as plt

N = 129
omega = 6.0 * 2.0 * np.pi/N

s = np.cos(omega * np.arange(N))

Fs = np.fft.fft(s)


#w = np.fft.ifftshift(signal.gaussian(N, 5.0,sym=False))
w = signal.gaussian(N, 5.0,sym=False)
w = w/np.sum(w)

plt.plot(w)
plt.show()

Fw = np.fft.fft( np.fft.ifftshift(w))
plt.plot(np.fft.fftfreq(N),Fw.real)
plt.plot(np.fft.fftfreq(N),Fw.imag)
plt.show()



s_f = signal.fftconvolve(s, w, mode='same')
Fs_f = np.fft.fft(s_f)

s_f2 = np.fft.ifft(Fs * Fw )

plt.plot(s,label='original')
plt.plot(s_f,label='FFTconvolve')
plt.plot(s_f2.real, label = 'FFT')
plt.legend()
plt.show()

plt.plot(np.fft.fftfreq(N),Fs.real, label ='Fs real')
plt.plot(np.fft.fftfreq(N),Fs.imag, label ='Fs imag')
plt.plot(np.fft.fftfreq(N),Fs_f.real, label ='Fs_f real')
plt.plot(np.fft.fftfreq(N),Fs_f.imag, label ='Fs_f imag')
plt.legend()
plt.show()


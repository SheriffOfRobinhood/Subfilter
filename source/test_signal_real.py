# -*- coding: utf-8 -*-
"""
Created on Wed Sep 16 10:39:06 2020

@author: paclk
"""
import numpy as np
from scipy import signal
import matplotlib.pyplot as plt

N = 128
omega = 6.0 * 2.0 * np.pi/N

s = np.cos(omega * np.arange(N))

Fs = np.fft.rfft(s)


#w = np.fft.ifftshift(signal.gaussian(N, 5.0,sym=False))
w = signal.gaussian(N, 5.0,sym=False)
w = w/np.sum(w)

plt.plot(w)
plt.plot(np.fft.ifftshift(w))
plt.show()

Fw = np.fft.rfft( np.fft.ifftshift(w))
plt.plot(np.fft.rfftfreq(N),Fw.real)
plt.plot(np.fft.rfftfreq(N),Fw.imag)
plt.show()


pad_len=N//2
field = np.pad(s, pad_len, mode='wrap')
result = signal.fftconvolve(field, w, mode='same')
s_f = result[pad_len+1:-(pad_len-1)]

#s_f = signal.fftconvolve(s, w, mode='same')
Fs_f = np.fft.rfft(s_f)

s_f2 = np.fft.irfft(Fs * Fw )

plt.plot(s,label='original')
plt.plot(s_f,label='FFTconvolve')
plt.plot(s_f2.real, label = 'FFT')
plt.legend()
plt.show()

plt.plot(np.fft.rfftfreq(N),Fs.real, label ='Fs real')
plt.plot(np.fft.rfftfreq(N),Fs.imag, label ='Fs imag')
plt.plot(np.fft.rfftfreq(N),Fs_f.real, label ='Fs_f real')
plt.plot(np.fft.rfftfreq(N),Fs_f.imag, label ='Fs_f imag')
plt.legend()
plt.show()

plt.plot(s_f2/s)
plt.plot(s_f/s)
plt.ylim([-1,1])
plt.show()

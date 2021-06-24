# -*- coding: utf-8 -*-
"""
Created on Thu Jun 24 10:36:17 2021

@author: Research_Lab
"""

import numpy as np
from scipy import fftpack
from matplotlib import pyplot as plt
#%%
import pandas
import os
datadir = os.path.join('C:\\Users\\Research Lab\\Desktop\\SI_Resources\\Synthetic_topography\\') # directory for some sample data files
filename = 'Elevation_profile2.csv'
filepath = os.path.join(datadir, filename)
df = pandas.read_csv(filepath)
print(df['Graphic Profile 1'])
sig = df.ix[:,1]
sig_non_zero = sig[sig!=0]
#%%
#%%
# The FFT of the signal
sig_fft = fftpack.fft(sig_non_zero)

# And the power (sig_fft is of complex dtype)
power = np.abs(sig_fft)**2

# The corresponding frequencies
sample_freq = fftpack.fftfreq(sig_non_zero.size, d=time_step)

# Plot the FFT power
plt.figure(figsize=(6, 5))
plt.plot(sample_freq, power)
plt.xlabel('Frequency [Hz]')
plt.ylabel('plower')

# Find the peak frequency: we can focus on only the positive frequencies
pos_mask = np.where(sample_freq > 0)
freqs = sample_freq[pos_mask]
peak_freq = freqs[power[pos_mask].argmax()]
peak_freqs = freqs[[0,3,8]]

# Check that it does indeed correspond to the frequency that we generate
# the signal with
#np.allclose(peak_freq, 1./period)

# An inner plot to show the peak frequency
axes = plt.axes([0.55, 0.3, 0.3, 0.5])
plt.title('Peak frequency')
plt.plot(freqs[:15], power[:15])
plt.setp(axes, yticks=[])

plt.figure(figsize=(9, 7))
plt.plot(freqs[2:100], power[2:100])

# scipy.signal.find_peaks_cwt can also be used for more advanced
# peak detection
#%%
time_vec = np.arange(0, sig_non_zero.size, 1)
#%%
#Setting a cutoff threshold at the first 

high_freq_fft_coarse = sig_fft.copy()
high_freq_fft_coarse[np.abs(sample_freq) > peak_freq] = 0
filtered_sig_coarse = fftpack.ifft(high_freq_fft_coarse)

high_freq_fft_fine = sig_fft.copy()
high_freq_fft_fine [np.abs(sample_freq) > peak_freqs[2]] = 0
filtered_sig_fine = fftpack.ifft(high_freq_fft_fine)

high_freq_fft_finest = sig_fft.copy()
high_freq_fft_finest [np.abs(sample_freq) > 0.5] = 0
filtered_sig_finest = fftpack.ifft(high_freq_fft_finest)

#%%
# Estimating the standard deviation for the noise

SD_sample_coarse = np.std(sig_non_zero-filtered_sig_coarse )
noise_random_coarse = np.random.normal(0, 1, sig_non_zero.size)
synthetic_topo_coarse = filtered_sig_coarse + noise_random_coarse

SD_sample_fine = np.std(sig_non_zero-filtered_sig_fine )
noise_random_fine = np.random.normal(0, 1, sig_non_zero.size)
synthetic_topo_fine = filtered_sig_fine + noise_random_fine

SD_sample_finest = np.std(sig_non_zero-filtered_sig_finest )
noise_random_finest = np.random.normal(0, 1, sig_non_zero.size)
synthetic_topo_finest = filtered_sig_finest + noise_random_finest


#%%
# Generating Plots

plt.figure(figsize=(10, 10))
plt.subplot(2, 2, 1)
plt.plot(time_vec, sig_non_zero, linewidth=2, label='Original signal')
plt.xlabel('Distance(m)')
plt.ylabel('Elevation(m)')
plt.legend(loc='best')

plt.subplot(2, 2, 2)
plt.plot(time_vec, sig_non_zero, linewidth=2, label='Original signal')
plt.plot(time_vec, synthetic_topo_coarse, linewidth=1, label='Syn topo coarse')
plt.xlabel('Distance(m)')
plt.ylabel('Elevation(m)')
plt.legend(loc='best')

plt.subplot(2, 2, 3)
plt.plot(time_vec, sig_non_zero, linewidth=2, label='Original signal')
plt.plot(time_vec, synthetic_topo_fine, linewidth=1, label='Syn topo fine')
plt.xlabel('Distance(m)')
plt.ylabel('Elevation(m)')
plt.legend(loc='best')

plt.subplot(2, 2, 4)
plt.plot(time_vec, sig_non_zero, linewidth=2, label='Original signal')
plt.plot(time_vec, synthetic_topo_finest, linewidth=1, label='Syn topo finest')
plt.xlabel('Distance(m)')
plt.ylabel('Elevation(m)')
plt.legend(loc='best')




#%%




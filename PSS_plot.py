"""PSS_plot.py
A set of plotting commands which uses the metadata of the signal to make plotting simpler.
"""
from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
import matplotlib.pyplot as plt
import numpy as np

plt.rcParams['figure.figsize'] = (8.0,6.0)

def pulse_plot(signal_object, N_pulses=1, pol_bin=0, freq_bin=0, start_time=0, phase=False, **kwargs):
    nBins_per_period = int(signal_object.MetaData.period//signal_object.TimeBinSize)
    if signal_object.SignalType == 'intensity':
        Y_label = 'Relative Intensity'
        Title = 'Pulse Intensity'
        row = freq_bin
    if signal_object.SignalType == 'voltage':
        Y_label = 'Relative Voltage'
        Title = 'Pulse Voltage'
        row = pol_bin

    if phase:
        Phase = np.linspace(0., N_pulses, N_pulses*nBins_per_period)
        #print(Phase.shape)
        #print(pulsar_object.signal[:N_pulses*pulsar_object.nBinsPeriod].shape)
        plt.plot(Phase, signal_object.signal[row,:N_pulses*nBins_per_period], lw=0.4, **kwargs)
        plt.xlabel('Phase')
        plt.ylabel(Y_label)
        plt.yticks([])
        plt.title(Title)
        plt.show()

    else:
        stop_time = start_time + N_pulses*nBins_per_period*signal_object.TimeBinSize
        Nx = N_pulses*nBins_per_period
        time = np.linspace(start_time, stop_time, Nx)
        plt.plot(time,signal_object.signal[row,:Nx], lw=0.4, **kwargs)
        plt.xlabel('Time (ms)')
        plt.ylabel(Y_label)
        plt.yticks([])
        plt.title(Title)
        plt.show()

def filter_bank(signal_object, pulsar_object, N_pulses=1, start_time=0, phase=False, **kwargs):
    nBins_per_period = int(signal_object.MetaData.period//signal_object.TimeBinSize)
    if signal_object.SignalType == 'intensity':
        #Y_label = 'Relative Intensity'
        Title = 'Filter Bank'
    if signal_object.SignalType == 'voltage':
        #Y_label = 'Relative Voltage'
        #Title = 'Filter Bank'
        raise ValueError('Filter Bank not supported for voltage signals at this time.')

    stop_bin = N_pulses*nBins_per_period
    Y_label = 'Frequency (MHz)'
    if phase:
        Extent = [0, N_pulses, signal_object.first_freq, signal_object.last_freq]
        plt.imshow(signal_object.signal[:,:stop_bin],origin='left',cmap='plasma',aspect='auto',extent=Extent,**kwargs)
        plt.xlabel('Phase')
        plt.ylabel(Y_label)
        #plt.yticks([])
        plt.title(Title)
        plt.show()

    else:
        time = len(signal_object.signal[0,:])
        stop_time = start_time + N_pulses * nBins_per_period * signal_object.TimeBinSize
        Extent = [start_time, stop_time, signal_object.first_freq, signal_object.last_freq]
        plt.imshow(signal_object.signal[:,:stop_bin],origin='left',cmap='plasma',aspect='auto',extent=Extent,**kwargs)
        plt.xlabel('Observation Time (ms)')
        plt.ylabel(Y_label)
        #plt.yticks([])
        plt.title(Title)
        plt.show()

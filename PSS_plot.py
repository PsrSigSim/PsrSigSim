"""PSS_plot.py
A set of plotting commands which uses the metadata of the signal to make plotting simpler.
"""
from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
import matplotlib.pyplot as plt
import numpy as np

__all__= ['profile_plot','pulse_plot','filter_bank']

plt.rcParams['figure.figsize'] = (8.0,6.0)
plt.rcParams.update({'font.size': 14})

def profile_plot(signal_object, freq_bin=0, phase=False, **kwargs):
    Title = 'Pulsar Profile Template'
    profile = signal_object.MetaData.profile
    Nx = profile.size

    if phase:
        Phase = np.linspace(0., 1, Nx)
        plt.plot(Phase, profile, lw=0.7, **kwargs)
        plt.xlabel('Phase')
        #plt.ylabel(Y_label)
        plt.yticks([])
        plt.ylim(0,profile.max()*1.05)
        plt.title(Title)
        plt.show()

    else:
        stop_time = signal_object.MetaData.period
        time = np.linspace(0, stop_time, Nx)
        plt.plot(time, profile, lw=0.7, **kwargs)
        plt.xlabel('Time (ms)')
        #plt.ylabel(Y_label)
        plt.ylim(0,profile.max()*1.05)
        #plt.yticks([])
        plt.title(Title)
        plt.show()


def pulse_plot(signal_object, N_pulses=1, pol_bin=0, freq_bin=0, start_time=0, phase=False, **kwargs):
    try:
        nBins_per_period = int(signal_object.MetaData.period//signal_object.TimeBinSize)
    except:
        ValueError('Need to sample pulses')
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
        plt.plot(Phase, signal_object.signal[row,:N_pulses*nBins_per_period], lw=0.4, **kwargs)
        plt.xlabel('Phase')
        #plt.ylabel(Y_label)
        #plt.yticks([])
        plt.title(Title)
        plt.show()

    else:
        stop_time = start_time + N_pulses*nBins_per_period*signal_object.TimeBinSize
        Nx = N_pulses*nBins_per_period
        time = np.linspace(start_time, stop_time, Nx)
        plt.plot(time,signal_object.signal[row,:Nx], lw=0.4, **kwargs)
        plt.xlabel('Time (ms)')
        #plt.ylabel(Y_label)
        #plt.yticks([])
        plt.title(Title)
        plt.show()

def filter_bank(signal_object, grid=False, N_pulses=1, start_time=0, phase=False, **kwargs):
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
        #ax = plt.gca();
        plt.xlabel('Phase')
        plt.ylabel(Y_label)
        #plt.yticks([])
        if grid:
            ax = plt.gca();
            ax.set_xticks([])#np.linspace(Extent[0], Extent[1], N_pulses*2));
            ax.set_yticks([])#np.linspace(Extent[2], Extent[3], 5));
            ax.set_xticks(np.linspace(Extent[0], Extent[1], N_pulses*10), minor=True);
            ax.set_yticks(np.linspace(Extent[2], Extent[3], signal_object.Nf), minor=True);
            ax.grid(which='both', color='black', linestyle='-', linewidth=0.5)
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
        if grid:
            ax = plt.gca();
            ax.set_xticks(np.linspace(Extent[0], Extent[1], 20), minor=True);
            ax.set_yticks(np.linspace(Extent[2], Extent[3], signal_object.Nf), minor=True);
            ax.grid(which='both', color='black', linestyle='-', linewidth=0.5)
        plt.title(Title)
        plt.show()

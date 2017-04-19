"""ism.py
module to cause dispersion
"""
from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
import numpy as np
import scipy as sp
from scipy import signal
from . import PSS_utils as utils
from . import scintillation as scint

class ISM(object):
    def __init__(self, Signal_in, debug=False):
        self.Signal_in = Signal_in
        self.signal = self.Signal_in.signal
        self.f0 = self.Signal_in.f0
        self.bw = self.Signal_in.bw
        self.Nf = self.Signal_in.Nf
        self.Nt = self.Signal_in.Nt
        self.TimeBinSize =self.Signal_in.TimeBinSize
        self.freqBinSize = self.Signal_in.freqBinSize
        self.first_freq = self.Signal_in.first_freq
        self.last_freq = self.Signal_in.last_freq
        self.freq_Array = self.Signal_in.freq_Array
        self.ISM_Dict = dict(dispersion=True, scattering=False, DM = 30, scintillation=False)

    def shiftit(self, y, shift):
        """
        shifts array y by amount shift (in sample numbers)
        uses shift theorem and FFT
        shift > 0  ==>  lower sample number (earlier)
        modeled after fortran routine shiftit
        Optimized from JMC's code by Michael Lam
        """
        #TODO Add Try Except for odd length arrays...
        yfft = np.fft.fft(y)
        size = np.size(y) #saves time
        constant = (shift*2*np.pi)/float(size) #needs a negative here for the right direction, put it in?
        theta = constant*np.arange(size)
        c = np.cos(theta)
        s = np.sin(theta)
        work = np.zeros(size, dtype='complex')
        work.real = c * yfft.real - s * yfft.imag
        work.imag = c * yfft.imag + s * yfft.real
        # enforce hermiticity

        work.real[size//2:] = work.real[size//2:0:-1]
        work.imag[size//2:] = -work.imag[size//2:0:-1]
        work[size//2] = 0.+0.j
        workifft = np.fft.ifft(work)
        return workifft.real

    def disperse(self, DM = 30, to_DM_Broaden = True):
        #Function to calculate the dispersion per frequency bin for 1/f^2 dispersion
        self.DM = DM
        self.ISM_Dict["DM"] = self.DM
        self.ISM_Dict['DM_Broaden'] = to_DM_Broaden
        if self.Signal_in.SignalType=='intensity':
            #For intensity signal calculat dispersion for all sub-bands.
            self.K = 1.0/2.41e-4 #constant used to be more consistent with PSRCHIVE
            self.time_delays = -1e3*self.K*self.DM*(np.power(self.freq_Array,-2)) #freq in MHz, delays in milliseconds
                #Dispersion as compared to infinite frequency
            self.time_delays //= self.TimeBinSize #Convert to number of bins
            for ii, freq in enumerate(self.freq_Array):
                self.signal[ii,:] = self.shiftit(self.signal[ii,:], self.time_delays[ii])
                sub_band_width = self.bw/self.Nf
                width = int(utils.top_hat_width(sub_band_width, freq, DM)//self.TimeBinSize)
                if width > 0 and to_DM_Broaden:
                    if width > self.Nt:
                        raise ValueError('Too Much DM! Dispersion broadening top hat wider than data array/')

                    self.signal[ii,:] = np.convolve(sp.signal.boxcar(width)/width, self.signal[ii,:] ,'same').astype(self.Signal_in.data_type)
                    # The division by width of the boxcar filter normalizes the convolution

                    #print(self.freq_Array[ii],' MHz ','width=', width) #for debugging
        elif self.Signal_in.SignalType=='voltage':
            #For voltage signal disperse coherently.
            raise ValueError('Sorry, Voltage-type signal dispersion is not currently supported!')
            #for ii in range(4): #Maybe faster to do the complex fft with two channels.
            #    sig_FFT = np.fft.rfft(self.signal[ii,:])
            #    fft_len = len(sig_FFT)
            #    f_array = np.linspace(-(self.last_freq)*1e6,0,length2)
            #    disp_signal_fft = sig_FFT*np.exp(1j*2*np.pi*4.148808e9/((freq+f0)*f0**2)*DM*freq**2)
            #    self.signal[ii,:] = np.fft.irfft(disp_signal_fft)

        self.Signal_in.MetaData.AddInfo(self.ISM_Dict)

#    def scintillate(self, scint_bandwidth=18e6, scint_time):
#        Nbins_per_scint_time = int(scint_time//self.TimeBinSize)
#        N_scint_timescales = int(self.Nt//Nbins_per_scint_time)
#        dims = 100
#        self.Phase_Screens = np.zeros((self.Nf,dims,dims))
#        self.Intensity_Screens = np.zeros(self.Phase_Screens.shape)
#        for ii, freq in enumerate(self.freq_Array):
#            Phase_Screens[str(freq)] = scint.phase_screen(freq, Nx = dims, Ny = dims, scint_bandwidth=scint_bandwidth)
#            Intensity_Screens[ii,:,:] = scint.images(Phase_Screens[str(freq)]).intensity
#
#        for jj in range(N_scint_timescales):
#            self.Intensity_Screens[,,]

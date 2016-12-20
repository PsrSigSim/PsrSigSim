"""ism.py
module to cause dispersion
"""
from __future__ import (absolute_import, division,
                        print_function, unicode_literals)

import numpy as np
import scipy as sp

class ISM(object):
    def __init__(self, Signal_in):
        self.Signal_in = Signal_in
        self.signal = self.Signal_in.signal
        self.f0 = self.Signal_in.f0
        self.bw = self.Signal_in.bw
        self.Nf = self.Signal_in.Nf
        self.Nt = self.Signal_in.Nt
        self.freqBinSize = self.bw/self.Nf
        self.first_freq = self.f0 - self.freqBinSize * self.Nf/2
        self.last_freq = self.f0 + self.freqBinSize * self.Nf/2
        self.phase = np.linspace(0., 1., self.Nt)
        #self.profile = 1./np.sqrt(2.*np.pi)/0.05 * np.exp(-0.5 * ((self.phase-0.25)/0.05)**2)
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


    def disperse(self, DM =30):
        #Function to calculate the dispersion per frequency bin for 1/f^2 dispersion
        self.DM = DM
        self.ISM_Dict["DM"] = self.DM
        self.K = 1.0/2.41e-4 #constant used to be more consistent with PSRCHIVE
        self.freq_Array = np.linspace(self.first_freq, self.last_freq, self.Nf,endpoint=False)
        self.time_delays = -1e3*self.K*self.DM*(np.power(self.freq_Array,-2)) #freq in MHz, delays in milliseconds
            #Dispersion as compared to infinite frequency
        for ii in range(0,self.Nf):
            self.signal[ii,:] = self.shiftit(self.signal[ii,:], self.time_delays[ii])

        self.Signal_in.MetaData.AddInfo(self.ISM_Dict)

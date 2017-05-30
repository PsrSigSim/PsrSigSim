"""telescope.py
module to add RFI and downsample data
"""

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
import numpy as np
import scipy as sp
from . import PSS_utils as utils

class Telescope(object):
    def __init__(self, Signal_in):
        self.signal_in = Signal_in.signal #Need to make new signal. This is the old one
        self.f0 = Signal_in.f0
        self.bw = Signal_in.bw
        self.Nf = Signal_in.Nf
        self.Nt = Signal_in.Nt
        self.TotTime = Signal_in.TotTime
        self.TimeBinSize = self.TotTime/self.Nt #in milliseconds
        #self.SignalSampFreq = 1/self.TimeBinSize

    def observe(self, telescope='GBT', Band=1400, mode='search', noise=False):
        # Method to Downsample the simulated signal into the telescopes sampling frequency
        # Telescope sampling frequency given in GHz
        #TODO Adds error messages if not the correct type of file.
        #TODO Need to write files to disk when they are large.
        Telescope_BW = {('GBT',820): 3.125, ('GBT',1400): 12.5,('AO',327): 1.5625, ('AO',430): 1.5625, ('AO',1400): 12.5, ('AO',2300): 12.5,}
        """
        The sampling rate for various telescope backends at given central frequencies.
        820 GBT GUPPI: 3.125 MHz
        1400 GBT GUPPI: 12.5 MHz

        327 AO PUPPI: 1.5625 MHz
        430 AO PUPPI: 1.5625 MHz
        1400 AO PUPPI: 12.5 MHz
        2300 AO PUPPI: 12.5 MHz
        """

        self.TelescopeTimeBinSize = 1/(2*Telescope_BW[telescope, Band])

        if self.TimeBinSize == self.TelescopeTimeBinSize:
            self.signal = self.Signal_in.signal

        elif self.TelescopeTimeBinSize % self.TimeBinSize == 0:
            SampFactor = int(self.TelescopeTimeBinSize // self.TimeBinSize)
            self.signal=np.zeros((self.Nf,int(self.Nt//SampFactor)))
            for ii in range(self.Nf):
                self.signal[ii,:] = utils.down_sample(self.signal_in[ii,:], SampFactor)
            print("Input signal sampling frequency= ", 1/self.TimeBinSize," kHz.\nTelescope sampling frequency = ",1/self.TelescopeTimeBinSize," kHz")

        elif self.TelescopeTimeBinSize > self.TimeBinSize:
            self.NewLength = int(self.TotTime//self.TelescopeTimeBinSize)
            self.signal=np.zeros((self.Nf,self.NewLength))
            for ii in range(self.Nf):
                self.signal[ii,:] = utils.rebin(self.signal_in[ii,:], self.NewLength)
            print("Input signal sampling frequency= ", self.TimeBinSize," ms. Telescope sampling frequency = ",self.TelescopeTimeBinSize," ms")

        else:
            # Throw error if the input signal has a lower sampling frequency than the telescope sampling frequency.
            raise ValueError("Signal Sampling Frequency Lower than Telescope Sampling Frequency")

        self.NFreqBins, self.NTimeBins = self.signal.shape

        if noise :
            self.signal = self.signal + 2000*np.random.randn(self.NFreqBins, self.NTimeBins)**2



    def fold(self, period, N_Folds = 100):
        self.period = period
        self.NBinsPeriod = int(self.period // self.TelescopeTimeBinSize)
        if self.NBinsPeriod*N_Folds > self.NTimeBins:
            raise ValueError("Not enough time for that many foldings!")
        self.folded = np.sum(self.signal[:,self.NBinsPeriod:self.NBinsPeriod*(N_Folds+1)].reshape(self.NFreqBins, N_Folds, self.NBinsPeriod),axis=1)
        #TODO Tweak since losing precision with many folds. Set by overall time.

"""pulsar.py
module to create pulses
"""
from __future__ import (absolute_import, division,
                        print_function, unicode_literals)

import numpy as np
import scipy as sp
import h5py
import math
import PSS_utils as utils

class Pulsar(object):
    def __init__(self, Signal_in, period = 50): #period in milliseconds
        self.Signal_in = Signal_in
        self.signal = self.Signal_in.signal
        self.f0 = self.Signal_in.f0
        self.bw = self.Signal_in.bw
        self.Nf = self.Signal_in.Nf
        self.Nt = self.Signal_in.Nt
        self.SignalType = self.Signal_in.SignalType
        self.TotTime = self.Signal_in.TotTime
        self.T = period
        self.TimeBinSize = self.TotTime/self.Nt
        self.nBinsPeriod = int(self.T/self.TimeBinSize)
        self.NPeriods = math.floor(self.TotTime/self.T) #Number of periods that can fit in the time given
        #self.time = np.linspace(0., self.TotTime, self.Nt)
        self.phase = np.linspace(0., 1., self.nBinsPeriod)
        self.profile = 1./np.sqrt(2.*np.pi)/0.05 * np.exp(-0.5 * ((self.phase-0.25)/0.05)**2)
        self.PulsarDict = dict(Profile="gaussian", peak=0.25, width=0.05, amplitude=1.)
        #phase = np.linspace(0., 1., self.Nt)
        #TODO Add ability to deal with multiple bands
        #TODO Check to see that you have the correct kind of array and info you need


    def draw_intensity_pulse(self, reps):
        """draw_intensity_pulse(pulse)
        draw a single pulse as bin by bin random process (gamma distr) from input template
        """
        pr = 10*np.tile(self.profile, reps)
        pulse = np.random.gamma(4., pr/4.) #pull from gamma distribution


        return pulse

    def draw_voltage_pulse(self, reps):
        """draw_voltage_pulse(pulse)
        draw a single pulse as bin by bin random process (normal distr) from input template
        """
        pr = 8*np.tile(self.profile, reps)
        L = len(pr)
        pulse = pr * np.random.normal(0, 1 , L) #pull from gaussian distribution

        return pulse


    #First profiles
    def gauss_template(self, peak=0.25, width=0.05, amp=1.):

        #TODO: error checking for array length consistency?
        #TODO: if any param is a not array, then broadcast to all entries of other arrays?

        try: # is this an array
            peak = np.array(peak)
            width = np.array(width)
            amp = np.array(amp)
            self.PulsarDict["amplitude"] = amp
            self.PulsarDict["Profile"] = "multiple gaussians"
            amp = amp/amp.sum()  # normalize sum
            profile = np.zeros(self.nBinsPeriod)
            #JEFF Can we use the built in numpy distribution here? I imagine that it's faster than this for loop.
            for ii in range(amp.size):
                norm = amp[ii]/np.sqrt(2.*np.pi)/width[ii]
                self.profile += norm * np.exp(-0.5 * ((self.phase-peak[ii])/width[ii])**2)
        except: # one gaussian
            norm = 1./np.sqrt(2.*np.pi)/width
            self.profile = norm * np.exp(-0.5 * ((self.phase-peak)/width)**2)
            self.PulsarDict["amplitude"] = amp
            self.PulsarDict["Profile"] = "gaussian"

        self.PulsarDict["peak"] = peak
        self.PulsarDict["width"] = width

    def user_template(self,template):
        # Function to make any given 1-dimensional numpy array into the profile
        #TODO Allow other input files
        #TODO Adds error messages if not the correct type of file.
        self.PulsarDict["Profile"] = "user_defined"
        #TODO Add I/O for adding attributes for user defined templates.
        self.PulsarDict["peak"] = "None"
        self.PulsarDict["width"] = "None"
        self.PulsarDict["amplitude"] = "None"
        #self.nBinsPeriod = len(template)
        #self.profile = template
        self.nBinsTemplate = len(template)

        if self.nBinsTemplate==self.nBinsPeriod:
            self.profile = template

        elif self.nBinsTemplate > self.nBinsPeriod:
            self.profile = utils.rebin(template, self.nBinsPeriod)
            print("User supplied template has been downsampled.")
            print("Input array length= ", self.nBinsTemplate,". Pulse template length= ",self.profile.size,".")

        else:
            TempPhase = np.linspace(0,1,len(template))
            ProfileFcn = sp.interpolate.interp1d(TempPhase, template, kind='cubic', bounds_error=True)
            self.profile = ProfileFcn(self.phase)
            print("User supplied template has been interpolated using a cubic spline.")
            print("Input array length was ", self.nBinsTemplate," bins. New pulse template length is ",self.profile.size,".")

        self.MinCheck = np.amin(self.profile)
        if self.MinCheck < 0 :
            self.profile = np.where(self.profile > 0, self.profile, self.profile-self.MinCheck)


    def pulses(self):
        #Function that makes pulses using the defined profile template
        #TODO Error message if the signal has already been pulsed

        self.PeriodFracRemain = self.TotTime/self.T - self.NPeriods #Home much of the last period remaining
        self.NLastPeriodBins = self.Nt - self.NPeriods * self.profile.size #index when last period starts
        pulseType = {"intensity":"draw_intensity_pulse", "voltage":"draw_voltage_pulse"}
        pulseTypeMethod = getattr(self, pulseType[self.SignalType])

        NRows = self.Nf
        
        if self.SignalType == 'voltage':
            self.profile = np.sqrt(self.profile)
            NRows = 4

        if self.Nt*NRows > 500000: #Limits the array size to 2.048 GB
            """The following limits the length of the arrays that we call from pulseTypeMethod(), by limiting the number
            of periods we pull from the distribution at one time. This is for machines with small amounts of memory, and
            is currently optimized for an 8GB RAM machine. In this process, most of the time is spent writing to disk.
            """
            self.ChunkSize = 5000
            if self.NPeriods < self.ChunkSize :
                self.ChunkSize = self.NPeriods
            self.Nchunks = self.NPeriods//self.ChunkSize
            self.NPeriodRemainder = self.NPeriods % self.ChunkSize
            for ii in range(self.Nchunks): #limits size of the array in memory
                self.signal[:, ii * self.ChunkSize * self.nBinsPeriod : ii * self.ChunkSize * self.nBinsPeriod + self.ChunkSize * self.nBinsPeriod] = np.tile(pulseTypeMethod(self.ChunkSize),(NRows,1))

            if self.NPeriodRemainder != 0 :
                self.signal[:,self.Nchunks * self.ChunkSize * self.nBinsPeriod:] = np.tile(pulseTypeMethod(self.NPeriodRemainder),(NRows,1))

        else:
            self.signal[:,0:self.NPeriods * self.profile.size] = np.tile(pulseTypeMethod(self.NPeriods),(NRows,1)) #Can be put into main flow for large RAM computers.

        self.LastPeriod = pulseTypeMethod(1)[0:self.NLastPeriodBins]
        #self.pulse[self.NPeriods * self.profile.size:] = self.LastPeriod
        self.signal[:,self.NPeriods * self.profile.size:] = np.tile(self.LastPeriod,(NRows,1))
        #for jj in range(self.Nf):
        #    self.signal[jj,:] = self.pulse
        #self.signal = np.tile(self.pulse,(self.Nf,1))

        self.Signal_in.MetaData.AddInfo(self.PulsarDict)

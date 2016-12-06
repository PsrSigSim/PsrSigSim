"""pulsar.py
module to create pulses
"""
from __future__ import (absolute_import, division,
                        print_function, unicode_literals)

import numpy as np
import scipy as sp
import math

class Pulsar(object):
    def __init__(self, Signal_in, period = 50): #period in milliseconds
        self.Signal_in = Signal_in
        self.signal = self.Signal_in.signal
        self.f0 = self.Signal_in.f0
        self.bw = self.Signal_in.bw
        self.Nf = self.Signal_in.Nf
        self.Nt = self.Signal_in.Nt
        self.TotTime = self.Signal_in.TotTime
        self.T = period
        self.TimeBinSize = self.TotTime/self.Nt
        self.nBinsPeriod = int(self.T/self.TimeBinSize)
        self.time = np.linspace(0., self.TotTime, self.Nt)
        self.phase = np.linspace(0., 1., self.nBinsPeriod)
        self.profile = 1./np.sqrt(2.*np.pi)/0.05 * np.exp(-0.5 * ((self.phase-0.25)/0.05)**2)
        self.PulsarDict = dict(Profile="gaussian", SignalType="intensity", peak=0.25, width=0.05, amplitude=1.)
        #phase = np.linspace(0., 1., self.Nt)
        #TODO Add ability to deal with multiple bands
        #TODO Check to see that you have the correct kind of array and info you need


    def draw_intensity_pulse(self):
        """draw_pulse(pulse)
        draw a single pulse as bin by bin random process (gamma distr) from input template
        """
        #TODO: draw pulse at arbitrary sample rate (presumed less than template?)
        #TODO: average template into new phase bins
        #ph = np.linspace(0., 1., self.Nt) # new phase bins
        pr = self.profile
        pulse = np.random.gamma(4., pr/4.) #pull from gamma distribution

        #TODO: interpolate single pulse back to full resolution

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
        self.PulsarDict["SignalType"] = "intensity"

        #return profile

    def down_sample(self, x, R): #Method to down sample an array by a factor
        #down_sample(array, downsampling factor)
        #This is fast, but not as general as possible
        try:
            x.reshape(-1, R)
            downsampled = x.reshape(-1, R).mean(axis=1)/np.amax(x)

        except:
            pad_size = math.ceil(float(x.size)/R)*R - x.size
            x_padded = np.append(x, np.zeros(pad_size)*np.NaN)
            x_padded.shape
            downsampled = sp.nanmean(x_padded.reshape(-1,R), axis=1)/np.amax(x)
        return downsampled

    def rebin(self, a, newLength):
        """rebin(old array, new number of bins)
        This is a very general downsampling rebinner, but has for loops and if
        statements.
        """
        #TODO Make this code run faster. Vectorize
        newBins = np.linspace(0, a.size, newLength,endpoint=False)
        width = math.ceil(a.size/newLength)
        a_rebin=np.zeros((newLength,width))*np.nan
        #Using NaN means that we do not have extra zeros in the array
        row = 0
        column = 0
        for ii in range(0, a.size):
            if ii < (newBins[row] + newBins[1]):
                a_rebin[row,column] = a[ii]
                column +=1
            else:
                column = 0
                row += 1
                a_rebin[row,column] = a[ii]
                column +=1

        a_rebinned = sp.nanmean(a_rebin,axis=1)
        #NaN mean does not count NaNs in total
        return a_rebinned*np.amax(a)/np.amax(a_rebinned)

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
            self.profile = self.rebin(template, self.nBinsPeriod)
            print("User supplied template has been downsampled.")
            print("Input array length= ", self.nBinsTemplate,". Pulse template length= ",self.profile.size,".")

        else:
            raise ValueError("Template Sampling Frequency Too Small")
            # Throw error if the user supplied template has less bins than the allowed number.
            # TODO Write interpolating script

    def pulses(self, type="intensity"):
        #Function that makes pulses using the defined profile template
        #TODO Error message if the signal has already been pulsed
        #TODO add ability to choose intensity or voltage
        #if type="intensity":
        self.NPeriods = math.floor(self.TotTime/self.T) #Number of periods that can fit in the time given
        self.PeriodFracRemain = self.TotTime/self.T - self.NPeriods #Home much of the last period remaining
        self.NLastPeriodBins = self.Nt - self.NPeriods * self.profile.size #index when last period starts
        self.pulse = []
        for ii in range(self.NPeriods):
            self.pulse = np.append(self.pulse, self.draw_intensity_pulse())

        self.LastPeriod = self.draw_intensity_pulse()[0:self.NLastPeriodBins]
        self.pulse = np.append(self.pulse,self.LastPeriod)


        for jj in range(self.Nf):
            self.signal[jj,:] = self.pulse

        self.Signal_in.MetaData.AddInfo(self.PulsarDict)

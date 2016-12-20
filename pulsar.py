"""pulsar.py
module to create pulses
"""
from __future__ import (absolute_import, division,
                        print_function, unicode_literals)

import numpy as np
import scipy as sp

min_float = np.nextafter(0.,1.)  # smallest nonzero float

class Pulsar(object):
    def __init__(self, Signal_in):
        self.Signal_in = Signal_in
        self.signal = self.Signal_in.signal
        self.f0 = self.Signal_in.f0
        self.bw = self.Signal_in.bw
        self.Nf = self.Signal_in.Nf
        self.Nt = self.Signal_in.Nt
        self.phase = np.linspace(0., 1., self.Nt)
        self.profile = 1./np.sqrt(2.*np.pi)/0.05 * np.exp(-0.5 * ((self.phase-0.25)/0.05)**2)
        self.PulsarDict = dict(Profile="gaussian", peak=0.25, width=0.05, amplitude=1.)
        #phase = np.linspace(0., 1., self.Nt)
        #TODO Add ability to deal with multiple bands
        #TODO Check to see that you have the correct kind of array and info you need


    def draw_intensity_pulse(self):
        """draw_pulse(pulse)
        draw a single pulse as bin by bin random process (gamma distr) from input template
        """
        #TODO: draw pulse at arbitrary sample rate (presumed less than template?)
        #TODO: phase and nbins do nothing until this is finished
        #TODO: average template into new phase bins
        #ph = np.linspace(0., 1., self.Nt) # new phase bins
        pr = self.profile
        try:
            pulse = np.random.gamma(4., pr/4.) #pull from gamma distribution
        except ValueError:
            pr[pr==0] = 10.*min_float
            pulse = np.random.gamma(4., pr/4.) #pull from gamma distribution
            pulse[pulse<101.*min_float]=0.

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
            profile = np.zeros(self.Nt)
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

        #return profile

    def user_template(self,template):
        # Function to make any given 1-dimensional numpy array into the profile
        #TODO Allow other input files
        #TODO Adds error messages if not the correct type of file.
        self.PulsarDict["Profile"] = "user_defined"
        #TODO Add I/O for adding attributes for user defined templates.
        self.PulsarDict["peak"] = "None"
        self.PulsarDict["width"] = "None"
        self.PulsarDict["amplitude"] = "None"
        self.profile = template

    def pulses(self):
        #Function that makes pulses using the defined profile template
        #TODO Error message if the signal has already been pulsed
        #TODO add ability to choose intensity or voltage
        for ii in range(self.Nf):
            self.signal[ii,:] = self.draw_intensity_pulse()

        self.Signal_in.MetaData.AddInfo(self.PulsarDict)

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)

import numpy as np
import scipy as sp
from scipy import signal
import h5py
import math
from . import PSS_utils as utils


class Simulation():
    """Class that simulates the pulsar signal
       Insert more information about the class
    """
    def __init__(self, signal_class, pulsar_class=None, ism_class=None, scint_class=None):
        self.S = signal_class
        self.MD = signal_class.MetaData
        self.scint_class = scint_class
        #self.time_dependent_DM = False
        self.time_dependent_scatter = False

        if not pulsar_class: #Raise error if pulsar_class is not passed.
            raise ValueError('No Pulsar Class specified.')
        else: self.P = pulsar_class

        #if self.MD.to_DM_Broaden:
            #Convolve the profile with a top hat.(@cassidymwagner, this is where your code will go)
            #If there is no flag for broadening set then this will throw an exception and skip this step.

        #if self.MD.to_Scatter_Broaden_exp:
            #Convolve the profile with a decay function.

        #if self.MD.to_Scatter_Broaden_stoch and not self.MD.time_dependent_scatter:
            #Convolve the profile with a decay function.



    def simulate(self):
        #If DM is set, then disperse the signal using the ism class
        try:
            self.ism = ism_class
            self.DM = ism.DM
        except: #If DM is not set then pass
            pass

        if self.MD.time_dependent_scatter or self.MD.time_dependent_DM or self.MD.to_Scintillate:
            #If some part of the code is time dependent, then run this loop.
            if self.MD.time_dependent_scatter:
                raise ValueError('Time dependent Scattering not supported at this point.')
            else:
                Scatter_time = 1
            if self.MD.time_dependent_DM:
                raise ValueError('Time dependent Scattering not supported at this point.')
            else:
                DM_time = 1
            if self.MD.time_dependent_DM:
                raise ValueError('Time dependent Scattering not supported at this point.')
            else:
                Scint_factor = 1

            scint_time = self.MD.scint_time/self.MD.scint_time_sample_rate
            scint_samples_per_obs = np.floor(self.MD.TotTime//(scint_time*1e3))
            print('scint_samples_per_obs',scint_samples_per_obs)
            gain_norm = self.scint_class.gain.max() #Could set to avoid clipping, but not sure it's needed.
            gain = self.scint_class.gain / gain_norm
            scint_end_bin = scint_samples_per_obs * scint_time*1e3 #integer number of bins in scint
            print('scint_end_bin',scint_end_bin)
            self.start_times = np.linspace(0, scint_end_bin, scint_samples_per_obs)
            orig_profile = np.copy(self.P.profile)
            scint_bins = int(scint_time//self.S.TimeBinSize)
            tweak = 8
            if len(self.start_times) > len(gain[0,:]):
                raise ValueError('Scattering Screen is not long enough to scintillate at this Dispersion timescale.')
            for ii, bin_time in enumerate(self.start_times) :
                self.P.profile = gain[:,ii,np.newaxis] * orig_profile
                self.P.profile /= (self.P.profile.max()/tweak)
                self.P.make_pulses(bin_time, bin_time + scint_time)
                bin = int(bin_time//self.S.TimeBinSize)
                #self.S.signal[:,bin : bin + scint_bins] = gain[:,ii,np.newaxis]*self.S.signal[:,bin: bin + scint_bins]
            #self.frig = gain[:,ii,np.newaxis] * orig_profile
        else: #Otherwise just make the pulses using the given profiles.
            self.P.make_pulses()

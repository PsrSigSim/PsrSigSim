
from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
import numpy as np
from ..utils.utils import make_quant

class PulseProfile(object):
    """base class for pulse profiles"""

    def make_pulses(self, signal, tobs, folded=False):
        """generate pulses from profile
        
        Required Args:
            signal (:class:`Signal`-like): signal object to store pulses
            tobs (float): observation time (sec)
            start (float): start time
        """

        if signal.sigtype is "RFSignal":
            self._make_rf_pulses()
        elif signal.sigtype is "BasebandSignal":
            self._make_bb_pulses()
        elif signal.sigtype is "FilterBankSignal":
            self._make_fb_pulses()
        else:
            raise NotImplementedError()

    def _make_rf_pulses(self):
        """generate radio frequency pulses"""
        raise NotImplementedError()

    def _make_bb_pulses(self):
        """generate baseband pulses"""
        raise NotImplementedError()

    def _make_fb_pulses(self):
        """generate filter bank pulses"""
        raise NotImplementedError()


    def make_pulses_old(self, start_time = 0, stop_time = None):
        """Function that makes pulses using the defined profile template.
        Note: 'intensity'-type signals pulled from a gamma distribution using draw_intensity_pulse(),
            'voltage'-type signals pulled from a gaussian distribution using draw_voltage_pulse().
        """

        if stop_time == None:
            stop_time = self.ObsTime
        start_bin = int(start_time // self.TimeBinSize )
        last_bin = int(stop_time // self.TimeBinSize )
        if stop_time == self.ObsTime and start_time == 0:
            N_periods_to_make = self.NPeriods
            delta_bins = self.Nt
        elif stop_time < self.ObsTime:
            delta_bins = last_bin - start_bin
            N_periods_to_make = int(delta_bins // self.nBinsPeriod)
        elif stop_time > self.ObsTime:
            last_bin = self.Nt
            delta_bins = last_bin - start_bin
            N_periods_to_make = int(delta_bins // self.nBinsPeriod)

        self.NLastPeriodBins = int(delta_bins % self.nBinsPeriod)#delta_bins - N_periods_to_make * self.nBinsPeriod #Length of last period
        pulseType = {"intensity":"draw_intensity_pulse", "voltage":"draw_voltage_pulse"}
        pulseTypeMethod = getattr(self, pulseType[self.SignalType])

        Nt_chunk = int(self.mem_size_limit // self.NRows)
        self.ChunkSize = int(Nt_chunk // self.nBinsPeriod) # Calc ChunkSize in integer number of periods


        if self.Nt * self.NRows > self.mem_size_limit and self.ChunkSize != 0 : #Limits the array size to 2.048 GB
            """The following limits the length of the arrays that we call from pulseTypeMethod(), by limiting the number
            of periods we pull from the distribution at one time. This is for machines with small amounts of memory, and
            is currently optimized for an 8GB RAM machine. In this process, most of the time is spent writing to disk.
            """

            self.Nchunks = int(N_periods_to_make // self.ChunkSize)
            self.NPeriodRemainder = int(N_periods_to_make % self.ChunkSize)
            pulse_start = time.time()

            for ii in range(self.Nchunks): #limits size of the array in memory
                self.signal[:, start_bin + ii * self.ChunkSize * self.nBinsPeriod : \
                            start_bin + (ii+1) * self.ChunkSize * self.nBinsPeriod] \
                            = pulseTypeMethod(self.ChunkSize)
                pulse_check = time.time()
                try: #Python 2 workaround. Python 2 __future__ does not have 'flush' kwarg.
                    print('\r{0:2.0}% sampled in {1:.2f} seconds.'.format((ii + 1)*100/self.Nchunks , pulse_check-pulse_start), end='', flush=True)
                except: #This is the Python 2 version.
                    print('\r{0:2.0}% sampled in {1:.2f} seconds.'.format((ii + 1)*100/self.Nchunks , pulse_check-pulse_start), end='')
                    sys.stdout.flush()

            if self.NPeriodRemainder != 0 :
                self.signal[:,start_bin + self.Nchunks * self.ChunkSize * self.nBinsPeriod : start_bin + (self.Nchunks * self.ChunkSize + self.NPeriodRemainder) * self.nBinsPeriod] = pulseTypeMethod(self.NPeriodRemainder)

        else:
            self.signal[:,start_bin:N_periods_to_make * self.nBinsPeriod] = pulseTypeMethod(N_periods_to_make) #Can be put into main flow for large RAM computers.

        self.LastPeriod = pulseTypeMethod(1)[:,0:self.NLastPeriodBins]
        self.signal[:,start_bin + N_periods_to_make * self.nBinsPeriod:start_bin + N_periods_to_make * self.nBinsPeriod + self.NLastPeriodBins] = self.LastPeriod


class GaussProfile(PulseProfile):
    """sum of guassian components

    The shape of the inputs determine the number of gaussian components
    in the pulse.
        single float  : Single pulse profile made of a single gaussian
        1-d array     : Single pulse profile made up of multiple gaussians
    where 'n' is the number of Gaussian components in the profile.

    Required Args:
        N/A
    
    Optional Args:
        peak (float): center of gaussian in pulse phase, default: ``0.5``
        width (float): stdev of pulse in pulse phase, default: ``0.1``
        amp (float): amplitude of pulse relative to other pulses,
            default: ``1``
        Nphase (int): number of phase bins, default: ``256``
    
    Pulses are normalized so that maximum is 1.
    See draw_voltage_pulse, draw_intensity_pulse and make_pulses() methods for more details.
    """
    def __init__(self, peak=0.5, width=0.1, amp=1, Nphase=256):
        #TODO: error checking for array length consistency?
        #TODO: if any param is a not array, then broadcast to all entries of other arrays?

        peak = np.array(peak)
        width = np.array(width)
        amp = np.array(amp)

        Amax = amp.max()

        profile = np.zeros(Nphase)
        ph = np.arange(Nphase)/Nphase
        for p, sig, A in zip(peak, width, amp):
            profile += A/Amax * np.exp(-0.5 * ((ph-p)/sig)**2)

        raise NotImplementedError()

class UserProfile(PulseProfile):
    """user specified pulse profile"""
    def __init__(self):
        raise NotImplementedError()


from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
import numpy as np
from ..utils.utils import make_quant

class PulseProfile(object):
    """base class for pulse profiles

    This class carries the pulse generation methods
    """
    _profile = None

    def __call__(self, Nphase=None):
        """a :class:`PulseProfile` returns the actual profile when called
        """
        if Nphase is None:
            if self._profile is None:
                msg = "profile not generated, returning `None`"
                print("Warning: "+msg)
            return self._profile
        else:
            self.gen_profile(Nphase)
            return self._profile

    def gen_profile(self, Nphase):
        """generate the profile

        This is implemented by the subclasses!

        Args:
            Nphase (int): number of phase bins
        """
        raise NotImplementedError()


    def _make_amp_pulse(self, samprate, period, Npulse=1):
        """generate an amplitude pulse

        This method should be used for radio frequency and basebanded
        pulses.

        Args:
            samprate (float): sample rate (MHz)
            period (float): pulse period (sec)
            Npulse (int): number of pulses to generate
        """
        sr = make_quant(samprate, "MHz")
        P = make_quant(period, "s")
        Nsamp = int((sr * P).decompose())
        if not hasattr(self, '_profile'):
            # profile not yet generated
            self.gen_profile(Nsamp)
        elif len(self._profile) != Nsamp:
            # profile at wrong sample rate
            self.gen_profile(Nsamp)
        self.gen_profile(Nsamp)

        Nchunk = 1 # int(__MAX_AR_SIZE__ // self.NRows)

        pr = np.tile(self.profile, reps)
        L = pr.shape
        pulse = pr * signal._draw_norm * np.random.normal(0, self.gauss_draw_sigma, L) #pull from gaussian distribution


        raise NotImplementedError()

    def _make_pow_pulse(self):
        """generate a power pulse
        
        This method should be used for filter bank pulses"""
        raise NotImplementedError()


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

    Pulses are normalized so that maximum is 1.
    See draw_voltage_pulse, draw_intensity_pulse and make_pulses() methods for more details.
    """
    def __init__(self, peak=0.5, width=0.05, amp=1):
        #TODO: error checking for array length consistency?
        #TODO: if any param is a not array, then broadcast to all entries of other arrays?

        self._peak = peak
        self._width = width
        self._amp = amp

    def gen_profile(self, Nphase):
        """generate the profile
        Args:
            Nphase (int): number of phase bins
        """
        profile = np.zeros(Nphase)
        ph = np.arange(Nphase)/Nphase
        try:
            for p, wid, A in zip(self.peak, self.width, self.amp):
                profile += A * np.exp(-0.5 * ((ph-p)/wid)**2)
        except TypeError:
            profile += (self.amp * 
                        np.exp(-0.5 * ((ph-self.peak)/self.width)**2))

        self._profile = profile / np.max(profile)

    @property 
    def profile(self):
        return self._profile

    @property
    def peak(self):
        return self._peak

    @property
    def width(self):
        return self._width

    @property
    def amp(self):
        return self._amp


class UserProfile(PulseProfile):
    """user specified pulse profile"""
    def __init__(self):
        raise NotImplementedError()

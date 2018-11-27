
from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
import numpy as np
from ..utils.utils import make_quant

class PulseProfile(object):
    """base class for pulse profiles

    Pulse profiles are INTENSITY series, even when using amplitude style
    signals (like :class:`BasebandSignal`).
    """
    _profile = None

    def __call__(self, phases=None):
        """a :class:`PulseProfile` returns the actual profile calculated
        at the specified phases when called
        """
        if phases is None:
            if self._profile is None:
                msg = "base profile not generated, returning `None`"
                print("Warning: "+msg)
            return self._profile
        else:
            return self.calc_profile(phases)

    def calc_profile(self, phases):
        """calculate the profile at specified phase(s)
        Args:
            phases (array-like): phases to calc profile

        This is implemented by the subclasses!
        """
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

    def init_profile(self, Nphase):
        """generate the profile
        Args:
            Nphase (int): number of phase bins
        """
        ph = np.arange(Nphase)/Nphase
        self._profile = self.calc_profile(ph)
        self._Amax = self._profile.max()
        self._profile /= self.Amax

    def calc_profile(self, phases):
        """calculate the profile at specified phase(s)
        Args:
            phases (array-like): phases to calc profile
        Note:
            The normalization can be wrong, if you have not run
            ``init_profile`` AND you are generating less than one
            rotation.
        """
        ph = np.array(phases)
        profile = np.zeros_like(ph)
        try:
            # sum of Gaussian components
            for p, wid, A in zip(self.peak, self.width, self.amp):
                profile += A * np.exp(-0.5 * ((ph-p)/wid)**2)
        except TypeError:
            # single Gaussian component
            profile += (self.amp * 
                        np.exp(-0.5 * ((ph-self.peak)/self.width)**2))
        
        Amax = self.Amax if hasattr(self, '_Amax') else np.max(profile)
        return profile / Amax

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

    @property
    def Amax(self):
        return self._Amax


class UserProfile(PulseProfile):
    """user specified pulse profile"""
    def __init__(self):
        raise NotImplementedError()

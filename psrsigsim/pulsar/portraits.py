
from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
import numpy as np
from astropy import log
from scipy.interpolate import CubicSpline as _cubeSpline


class PulsePortrait(object):
    """base class for pulse profiles

    A pulse portrait is a set of profiles across the frequency range. They are
    INTENSITY series, even when using amplitude style signals
    (like :class:`BasebandSignal`).
    """
    _profile = None

    def __call__(self, phases=None):
        """a :class:`PulsePortrait` returns the actual profiles calculated
        at the specified phases when called
        """
        if phases is None:
            if self._profiles is None:
                msg = "base profiles not generated, returning `None`"
                print("Warning: "+msg)
            return self._profiles
        else:
            return self.calc_profiles(phases)

    def init_profiles(self, Nphase):
        """generate the profile, evenly sampled

        Args:
            Nphase (int): number of phase bins
        """
        ph = np.arange(Nphase)/Nphase
        self._profiles = self.calc_profiles(ph)
        self._Amax = self._profiles.max()
        self._profiles /= self.Amax

    def calc_profiles(self, phases):
        """calculate the profiles at specified phase(s)
        Args:
            phases (array-like): phases to calc profile
        Note:
            The normalization can be wrong, if you have not run
            ``init_profiles`` AND you are generating less than one
            rotation.

        This is implemented by the subclasses!
        """
        raise NotImplementedError()

    @property
    def profiles(self):
        return self._profiles

    @property
    def Amax(self):
        return self._Amax



class GaussPortrait(PulsePortrait):
    """sum of gaussian components

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
        # profile = np.zeros_like(ph)
        try:
            # sum of Gaussian components
            profile = (self.amp[:,np.newaxis] *
                       np.exp(-0.5 * ((ph[np.newaxis,:]
                                       -self.peak[:,np.newaxis])
                                       /self.width[:,np.newaxis])**2))
        except TypeError:
            # single Gaussian component
            profile += (self.amp *
                        np.exp(-0.5 * ((ph-self.peak)/self.width)**2))

        Amax = self.Amax if hasattr(self, '_Amax') else np.amax(profile)
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

class DataPortrait(PulsePortrait):
    """a pulse portrait generated from data

    The data are samples of the profiles at specified phases. If you have a
    functional form for the _profiles use :class:`UserProfile` instead.

    Parameters:
    -----------

    profiles : array, list of lists
        Profile data in 2-d array.


    phases : array, list of lists (optional)
        List of sampled phases. If phases are omitted profile is assumed to be
        evenly sampled and cover one whole rotation period.

    Profile is renormalized so that maximum is 1.
    See draw_voltage_pulse, draw_intensity_pulse and make_pulses() methods for
    more details.
    """
    def __init__(self, profiles, phases=None):
        if phases is None:
            # infer phases
            N = profiles.shape[1]
            if any([ii != jj for ii,jj in zip(profiles[:,0], profiles[:,-1])]):
                # enforce periodicity!
                profiles = np.append(profiles, profiles[:,0][:,np.newaxis],
                                     axis=1)
                phases = np.arange(N+1)/N
            else:
                phases = np.arange(N)/N
        else:
            if phases[-1] != 1:
                # enforce periodicity!
                phases = np.append(phases, 1)
                profiles = np.append(profiles, profiles[:,0][:,np.newaxis],
                                     axis=1)
            elif any([ii != jj for ii,jj in zip(profiles[:,0],
                                                profiles[:,-1])]):
                # enforce periodicity!
                profiles[-1,:] = profiles[0,:]

        self._generator = _cubeSpline(phases, profiles, axis=1,
                                      bc_type='periodic')

    def calc_profiles(self, phases):
        """calculate the profile at specified phase(s)
        Args:
            phases (array-like): phases to calc profile
        Note:
            The normalization can be wrong, if you have not run
            ``init_profile`` AND you are generating less than one
            rotation.
        """
        profiles = self._generator(phases)
        Amax = self.Amax if hasattr(self, '_Amax') else np.max(profiles)
        return profiles / Amax


class UserPortrait(PulsePortrait):
    """user specified 2-D pulse portrait"""
    def __init__(self):
        raise NotImplementedError()

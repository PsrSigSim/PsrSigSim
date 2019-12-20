
from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
import numpy as np
from scipy import stats
from .profiles import GaussProfile
from .profiles import UserProfile
from .profiles import DataProfile
from ..utils.utils import make_quant

class Pulsar(object):
    """class for pulsars

    The minimal data to instantiate a pulsar is the period, Smean, and
    pulse profile. The Profile is supplied via a :class:`PulseProfile`-like
    object.

    Args:
        period (float): pulse period (sec)
        Smean (float): mean pulse flux density (Jy)
        profile (:class:`PulseProfile`): pulse profile or 2-D pulse portrait
        name (string): name of pulsar
    """
    #TODO Other data could be supplied via a `.par` file.
    def __init__(self, period, Smean, profiles=None, name=None):
        self._period = make_quant(period, 's')
        self._Smean = make_quant(Smean, 'Jy')

        self._name = name

        # Assign profile class; default to GaussProfile if nothing is specified
        if profiles is None:
            self._Profiles = GaussProfile()
        else:
            self._Profiles = profiles

    def __repr__(self):
        namestr = "" if self.name is None else self.name+", "
        return "Pulsar("+namestr+"{})".format(self.period.to('ms'))

    @property
    def Profiles(self):
        return self._Profiles

    @property
    def name(self):
        return self._name

    @property
    def period(self):
        return self._period

    @property
    def Smean(self):
        return self._Smean

    def make_pulses(self, signal, tobs):
        """generate pulses from Profiles, :class:`PulsePortrait` object

        Required Args:
            signal (:class:`Signal`-like): signal object to store pulses
            tobs (float): observation time (sec)
        """
        signal._tobs = make_quant(tobs, 's')

        # init base profile at correct sample rate
        Nph = int((signal.samprate * self.period).decompose())

        #If there isn't enough frequencies in the profiles
        # And if it is a :class:`Profile` instance, reshape.
        if (signal.Nchan != self.Profiles.profiles.shape[0]
            and hasattr(self.Profiles, 'set_Nchan')):
            self.Profiles.set_Nchan(signal.Nchan)
        self.Profiles.init_profiles(Nph)

        # select pulse generation method
        if signal.sigtype in ["RFSignal", "BasebandSignal"]:
            self._make_amp_pulses(signal)
        elif signal.sigtype is "FilterBankSignal":
            self._make_pow_pulses(signal)
        else:
            msg = "no pulse method for signal: {}".format(signal.sigtype)
            raise NotImplementedError(msg)

        # compute Smax (needed for radiometer noise level)
        pr = self.Profiles()
        nbins = len(self.Profiles()) # Think this assumes a single profile for now...
        signal._Smax = self.Smean * nbins / np.sum(pr)

    def _make_amp_pulses(self, signal):
        """generate amplitude pulses

        This method should be used for radio frequency and basebanded
        pulses.

        Args:
            signal (:class:`Signal`-like): signal object to store pulses
        """
        # generate several pulses in time
        distr = stats.norm()

        Nsamp = int((signal.tobs * signal.samprate).decompose())
        signal.init_data(Nsamp)

        # TODO break into blocks
        # TODO phase from .par file
        # calc profile at phases
        phs = (np.arange(Nsamp) /
                (signal.samprate * self.period).decompose().value)
        phs %= 1  # clip integer part

        # convert intensity profile to amplitude!
        full_prof = np.sqrt(self.Profiles.calc_profiles(phs))

        signal._data = full_prof * distr.rvs(size=signal.data.shape)

    def _make_pow_pulses(self, signal):
        """generate a power pulse

        This method should be used for filter bank pulses

        Args:
            signal (:class:`Signal`-like): signal object to store pulses
        """
        if signal.subint:
            # Determine how many subints to make
            if signal.sublen is None:
                signal._sublen = signal.tobs
                signal._nsub = 1
            else:
                # This should be an integer, if not, will round
                signal._nsub = int(np.round(signal.tobs / signal.sublen))

            # determine the number of data samples necessary
            signal._nsamp = int((signal.nsub*(self.period*signal.samprate)).decompose())
            # Need to make an initial empty data array
            signal.init_data(signal.nsamp)

            # Tile the profiles to number of desired subints
            sngl_prof = np.tile(self.Profiles(), signal.nsub)

            # changed to number of subints
            signal._Nfold = (signal.sublen / self.period).decompose()
            distr = stats.chi2(df=signal.Nfold)
            signal._set_draw_norm(df=signal.Nfold)

            #Why is there a second call to init_data?
            signal.init_data(sngl_prof.shape[1])
            signal._data = (sngl_prof * distr.rvs(size=signal.data.shape)
                            * signal._draw_norm)
        else:
            # generate several pulses in time
            distr = stats.chi2(df=1)
            signal._set_draw_norm(df=1)

            signal._nsamp = int((signal.tobs * signal.samprate).decompose())
            signal.init_data(signal.nsamp)

            # TODO break into blocks
            # TODO phase from .par file
            # calc profiles at phases
            phs = (np.arange(signal.nsamp) /
                   (signal.samprate * self.period).decompose().value)
            phs %= 1  # clip integer part
            full_prof = self.Profiles.calc_profiles(phs)

            signal._data = (full_prof * distr.rvs(size=signal.data.shape)
                            * signal._draw_norm)

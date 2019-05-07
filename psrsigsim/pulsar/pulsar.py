
from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
import numpy as np
from scipy import stats
from .profiles import GaussProfile
from ..utils.utils import make_quant

class Pulsar(object):
    """class for pulsars

    The minimal data to instatiate a pulsar is the period, intensity, and
    pulse profile. The Profile is supplied via a :class:`PulseProfile`-like
    object.

    Args:
        period (float): pulse period (sec)
        intensity (float): mean pulse intensity (Jy)
        profile (:class:`PulseProfile`): pulse profile or 2-D pulse portrait
        name (string): name of pulsar
    """
    #TODO Other data could be supplied via a `.par` file.
    def __init__(self, period, intensity, profile=None, name=None):
        self._period = make_quant(period, 's')
        self._intensity = make_quant(intensity, 'Jy')

        self._name = name

        if profile is None:
            self._Profile = GaussProfile()
        else:
            self._Profile = profile

    def __repr__(self):
        namestr = "" if self.name is None else self.name+", "
        return "Pulsar("+namestr+"{})".format(self.period.to('ms'))

    @property
    def Profile(self):
        return self._Profile

    @property
    def name(self):
        return self._name

    @property
    def period(self):
        return self._period

    @property
    def intensity(self):
        return self._intensity

    def make_pulses(self, signal, tobs):
        """generate pulses from Profile, :class:`PulseProfile` object

        Required Args:
            signal (:class:`Signal`-like): signal object to store pulses
            tobs (float): observation time (sec)
        """
        signal._tobs = make_quant(tobs, 's')

        # init base profile at correct sample rate
        Nph = int((signal.samprate * self.period).decompose())
        self.Profile.init_profile(Nph)

        # select pulse generation method
        if signal.sigtype in ["RFSignal", "BasebandSignal"]:
            self._make_amp_pulses(signal)
        elif signal.sigtype is "FilterBankSignal":
            self._make_pow_pulses(signal)
        else:
            msg = "no pulse method for signal: {}".format(signal.sigtype)
            raise NotImplementedError(msg)

        # compute Smax (needed for radiometer noise level)
        pr = self.Profile()
        dph = 1 / len(self.Profile())
        norm = 1 / signal.Nchan
        signal._Smax = self.intensity / (np.sum(pr) * dph * norm)

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
        full_prof = np.sqrt(self.Profile.calc_profile(phs))

        signal._data = full_prof * distr.rvs(size=signal.data.shape)

    def _make_pow_pulses(self, signal):
        """generate a power pulse

        This method should be used for filter bank pulses

        Args:
            signal (:class:`Signal`-like): signal object to store pulses
        """
        if signal.subint:
            # generate one pulse in phase
            sngl_prof = self.Profile()

            signal._Nfold = (signal.tobs / self.period).decompose()
            distr = stats.chi2(df=signal.Nfold)
            signal._set_draw_norm(df=signal.Nfold)

            signal.init_data(len(sngl_prof))
            signal._data = (sngl_prof * distr.rvs(size=signal.data.shape)
                            * signal._draw_norm)
        else:
            # generate several pulses in time
            distr = stats.chi2(df=1)
            signal._set_draw_norm(df=1)

            Nsamp = int((signal.tobs * signal.samprate).decompose())
            signal.init_data(Nsamp)

            # TODO break into blocks
            # TODO phase from .par file
            # calc profile at phases
            phs = (np.arange(Nsamp) /
                   (signal.samprate * self.period).decompose().value)
            phs %= 1  # clip integer part
            full_prof = self.Profile.calc_profile(phs)

            signal._data = (full_prof * distr.rvs(size=signal.data.shape)
                            * signal._draw_norm)

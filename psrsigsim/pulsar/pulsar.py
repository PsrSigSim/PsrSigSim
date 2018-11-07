
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
        tobs = make_quant(tobs, 's')

        # select pulse generation method
        if signal.sigtype in ["RFSignal", "BasebandSignal"]:
            self._make_amp_pulses(signal, tobs)
        elif signal.sigtype is "FilterBankSignal":
            self._make_pow_pulses(signal, tobs)
        else:
            msg = "no pulse method for signal: {}".format(signal.sigtype)
            raise NotImplementedError(msg)

    def _make_amp_pulses(self, signal, tobs):
        """generate amplitude pulses

        This method should be used for radio frequency and basebanded
        pulses.

        Args:
            signal (:class:`Signal`-like): signal object to store pulses
            tobs (float): observation time (sec)
        """
        # loop over tobs generating single pulses
        # TODO: if using a pulsar ephemeris adjust the single pulse spacing!

        # N samples in observation
        Nsamp_obs = int((tobs * signal.samprate).decompose())
        signal.init_data(Nsamp_obs)
        # N samples per pulse period
        Nsamp_period = int((self.period * signal.samprate).decompose())
        # N samples per memory chunk
        Nsamp_chunk = int(__MAX_AR_SIZE__ // self.Nchan)

        # break calculation into N chunks
        Nchunk = Nsamp_obs // Nsamp_chunk +1

        # number of pulse periods per memory chunk
        Np_chunk = Nsamp_chunk // Nsamp_period


    def _make_pow_pulses(self, signal, tobs):
        """generate a power pulse

        This method should be used for filter bank pulses

        Args:
            signal (:class:`Signal`-like): signal object to store pulses
            tobs (float): observation time (sec)
        """
        #TODO this works for float32 ONLY, fix for int8
        # i.e. use _draw_norm...
        Nperiod = (tobs / self.period).decompose()
        Nph = int((signal.samprate * self.period).decompose())
        sngl_prof = self.Profile(Nph)
        if signal.subint:
            # generate one pulse in phase
            distr = stats.chi2(df=Nperiod)
            signal._set_draw_norm(df=Nperiod)

            signal.init_data(Nph)
            signal._data = (sngl_prof * distr.rvs(size=signal.data.shape)
                            * signal._draw_norm)
        else:
            # generate several pulses in time
            distr = stats.chi2(df=1)
            Nsamp = int((tobs / signal.samprate).decompose())
            signal.init_data(Nsamp)


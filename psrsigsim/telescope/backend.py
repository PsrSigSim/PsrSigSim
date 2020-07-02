
from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
import numpy as np

from ..utils.utils import down_sample, rebin, make_quant

__all__ = ['Backend']

class Backend(object):
    def __init__(self, samprate=None, name=None):
        self._name = name
        # Should be in MHz
        self._samprate = make_quant(samprate, "MHz")

    def __repr__(self):
        return "Backend({:s})".format(self._name)

    @property
    def name(self):
        return self._name

    @property
    def samprate(self):
        return self._samprate

    def adc(self, signal):
        """analog-digital-converter

        convert "analog" float Signal to "digital" 8-bit int
        """


    def fold(self, signal, pulsar):
        """fold data using pulsar ephemeris
        currently only uses pulsar period

        Args:
            signal (:class: FilterbankSignal): simulated signal 
            psr (:class:`Pulsar`): observed pulsar
        """
        #TODO: use .par file if available
        period = pulsar.period
        Nf, Nt = signal.data.shape
        Npbins = int(period * 2*signal.samprate)  # number of phase bins
        N_fold = Nt // Npbins  # number of folds
        fold_sig = signal.data[:, Npbins:Npbins*(N_fold+1)].reshape(
                                                         Nf, N_fold, Npbins//2)
        return np.sum(fold_sig, axis=1)

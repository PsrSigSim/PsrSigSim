
from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
import numpy as np

from ..utils.utils import down_sample, rebin

__all__ = ['Backend']

class Backend(object):
    def __init__(self, samprate=None, name=None):
        self._name = name
        self._samprate = samprate

    def __repr__(self):
        return "Backend({:s})".format(self._name)

    @property
    def name(self):
        return self._name

    @property
    def samprate(self):
        return self._samprate

    def fold(self, signal, psr):
        """fold data using pulsar ephemeris
        currently only uses pulsar period

        Args:
            signal (array): data to fold
            psr (:class:`Pulsar`): observed pulsar
        """
        #TODO: use .par file if available
        period = psr.T
        Nf, Nt = signal.shape
        Npbins = int(period * 2*self.samprate)  # number of phase bins
        N_fold = Nt // Npbins  # number of folds
        fold_sig = signal[:, Npbins:Npbins*(N_fold+1)].reshape(
                                                         Nf, N_fold, Npbins)
        return np.sum(fold_sig, axis=1)

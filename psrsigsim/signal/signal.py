
from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
import numpy as np

from ..utils.utils import make_quant

__all__ = ["Signal"]

class Signal(object):
    """base class for signals, subclass from this
    
    Required Args:
        f_cent [float]: central radio frequency

        bandwidth [float]: radio bandwidth of signal
        
    Optional Args:
        sample_rate [float]: sample rate of data, default: ``None``
            If no ``sample_rate`` is given the observation will default to
            the Nyquist frequency. Sub-Nyquist sampling is allowed, but a
            warning will be generated.

        dtype [type]: data type of array, default: ``np.float32``
    """

    _sigtype = "Signal"

    def __init__(self,
                 f_cent, bandwidth,
                 sample_rate=None,
                 fold=False,
                 dtype=np.float32):
        
        self._fcent = f_cent
        self._bw = bandwidth
        self._sr = sample_rate
        self._dtype = dtype

    def __repr__():
        return self.sigtype+"({0}, {1})".format(self.fcent, self.bw)

    def __add__(self, b):
        """overload ``+`` to concatinate signals"""
        raise NotImplementedError()

    def init_array(self, tobs, fold=False):
        """initialize an array to store the signal
        Required Args:
            tobs (float): observation time (sec)

        Optional Args:
            fold (bool): initialize for a pre-folded observation,
                default: ``False``
        """
        raise NotImplementedError()
        #self._tobs = make_quant(tobs, 's')

    @property
    def sigtype(self):
        return self._sigtype

    @property
    def fcent(self):
        return self._fcent

    @property
    def bw(self):
        return self._bw

    @property
    def tobs(self):
        if hasattr(self, '_tobs'):
            return self._tobs
        else:
            return None

    @property
    def sr(self):
        return self._sr

    @property
    def dtype(self):
        return self._dtype


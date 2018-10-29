
from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
import numpy as np

from .signal import Signal
from ..utils.utils import make_quant

class RFSignal(Signal):
    """a time domain signal at true radio frequency sampling
    
    RFSignals must be sampled at twice the maximum resolved frequency, i.e.
    a few GHz.  As such, RFSignals take up a TON of memory. Consider using
    BasebandSignal if this is a concern.

    Required Args:
        f_cent [float]: central radio frequency (MHz)

        bandwidth [float]: radio bandwidth of signal (MHz)
        
    Optional Args:
        sample_rate [float]: sample rate of data (MHz), default: ``None``
            If no ``sample_rate`` is given the observation will default to
            the Nyquist frequency. Sub-Nyquist sampling is allowed, but a
            warning will be generated.

        dtype [type]: data type of array, default: ``np.float32``
    """

    _sigtype = "RFSignal"

    def __init__(self,
                 f_cent, bandwidth,
                 sample_rate=None,
                 fold=False,
                 dtype=np.float32):

        self._fcent = make_quant(f_cent, 'MHz')
        self._bw = make_quant(bandwidth, 'MHz')
        
        f_Nyquist = 2 * (self._fcent + self._bw/2)
        if sample_rate is None:
            self._sr = f_Nyquist
        else:
            self._sr = make_quant(sample_rate, 'MHz')
            if self._sr < f_Nyquist:
                msg = ("specified sample rate {} < Nyquist frequency {}"
                       .format(self._sr, f_Nyquist))
                print("Warning: "+msg)

        self._dtype = dtype

    def to_RF(self):
        """convert signal to RFSignal
        """
        return self

    def to_Baseband(self):
        """convert signal to BasebandSignal
        """
        raise NotImplementedError()

    def to_FilterBank(self, subbw):
        """convert signal to FilterBankSignal
        """
        raise NotImplementedError()


from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
import numpy as np

from .signal import BaseSignal
from ..utils.utils import make_quant

class RFSignal(BaseSignal):
    """a time domain signal at true radio frequency sampling

    RFSignals must be sampled at twice the maximum resolved frequency, i.e.
    a few GHz.  As such, RFSignals take up a TON of memory. Consider using
    BasebandSignal if this is a concern.

    Required Args:
        fcent [float]: central radio frequency (MHz)

        bandwidth [float]: radio bandwidth of signal (MHz)

    Optional Args:
        sample_rate [float]: sample rate of data (MHz), default: ``None``
            If no ``sample_rate`` is given the observation will default to
            the Nyquist frequency. Sub-Nyquist sampling is allowed, but a
            warning will be generated.

        dtype [type]: data type of array, default: ``np.float32``
            any numpy compatible floating point type will work
    """

    _sigtype = "RFSignal"
    _Nchan = 2

    def __init__(self,
                 fcent, bandwidth,
                 sample_rate=None,
                 dtype=np.float32):

        self._fcent = make_quant(fcent, 'MHz')
        self._bw = make_quant(bandwidth, 'MHz')

        f_Nyquist = 2 * (self._fcent + self._bw/2)
        if sample_rate is None:
            self._samprate = f_Nyquist
        else:
            self._samprate = make_quant(sample_rate, 'MHz')
            if self._samprate < f_Nyquist:
                msg = ("specified sample rate {} < Nyquist frequency {}"
                       .format(self._samprate, f_Nyquist))
                print("Warning: "+msg)

        self._dtype = dtype

    def to_RF(self):
        """
        Convert signal to RFSignal

        Returns
        -------

        psrsigsim.signal.RFSignal
        """
        return self

    def to_Baseband(self):
        """
        Convert signal to BasebandSignal
        """
        # rfft this signal
        # clip array
        # irfft
        raise NotImplementedError()

    def to_FilterBank(self, Nsubband=512):
        """
        Convert signal to FilterBankSignal

        Parameters
        ----------

        Nsubband : int
            Number of frequency subbands.
        """
        # BB = self.to_Basband()
        # return BB.to_FilterBank()
        raise NotImplementedError()

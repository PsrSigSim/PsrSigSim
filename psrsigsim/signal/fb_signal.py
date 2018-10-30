
from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
import numpy as np

from .signal import Signal
from ..utils.utils import make_quant

class FilterBankSignal(Signal):
    """a filter bank signal, breaking the time domain signal into RF bins

    Unlike purely time domain signals, :class:`FilterBankSingnal`s are 2-D
    arrays. Filter banks record the intensity of the signal and can be much
    more sparsely sampled. The time binning must accurately capture the
    pulse profile, not the RF oscillations. In practice a filter bank is
    generated from the observed time domain signal by the telescope backend.
    We allow for direct filter bank signals to save memory.

    Required Args:
        fcent [float]: central radio frequency (MHz)

        bandwidth [float]: radio bandwidth of signal (MHz)

        sub_bw [float]: radio bandwidth of subbands (MHz)

        t_obs [float]: observation time (sec)

    Optional Args:
        sample_rate [float]: sample rate of data (MHz), default: ``None``
            If no ``sample_rate`` is given the observation will default to
            the Nyquist frequency. Sub-Nyquist sampling is allowed, but a
            warning will be generated.

        dtype [type]: data type of array, default: ``np.float32``
    """

    _sigtype = "FilterBankSignal"

    def __init__(self,
                 fcent, bandwidth,
                 t_obs,
                 sample_rate=None,
                 fold=False,
                 dtype=np.float32):
        raise NotImplementedError()

        self._fcent = make_quant(fcent, 'MHz')
        self._bw = make_quant(bandwidth, 'MHz')
        self._subbw = make_quant(sub_bw, 'MHz')

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

    @property
    def subbw(self):
        return self._subbw

    def to_RF(self):
        """convert signal to RFSignal
        """
        raise NotImplementedError()

    def to_Baseband(self):
        """convert signal to BasebandSignal
        """
        raise NotImplementedError()

    def to_FilterBank(self):
        """convert signal to FilterBankSignal
        """
        return self


from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
import numpy as np

from .signal import BaseSignal
from ..utils.utils import make_quant

class FilterBankSignal(BaseSignal):
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

    Optional Args:
        Nsubband [int]: number of sub-bands, default ``2048``

        sample_rate [float]: sample rate of data (MHz), default: ``None``
            If no ``sample_rate`` is given the observation will default to
            the 20.48 us per sample (about 50 kHz).  This is the sample rate
            for coherently dedispersed filter banks using XUPPI backends.

        subint [bool]: is this a folded subintegration, default ``False``

        dtype [type]: data type of array, default: ``np.float32``
    """
    #TODO: full stokes.  Currently this is just stokes-I
    #  add flag `fullstokes=False`
    #  data -> dict? keyed with "I", "Q", "U", "V"

    _sigtype = "FilterBankSignal"

    def __init__(self,
                 fcent, bandwidth,
                 Nsubband=2048,
                 sample_rate=None,
                 subint=False,
                 dtype=np.float32):

        self._fcent = make_quant(fcent, 'MHz')
        self._bw = make_quant(bandwidth, 'MHz')

        self._subint = subint

        if sample_rate is None:
            self._samprate = (1/make_quant(20.48, 'us')).to('MHz')
        else:
            self._samprate = make_quant(sample_rate, 'MHz')
            if self._samprate < f_Nyquist:
                msg = ("specified sample rate {} < Nyquist frequency {}"
                       .format(self._samprate, f_Nyquist))
                print("Warning: "+msg)

        self._Nchan = Nsubband

        self._dtype = dtype
        self._set_draw_norm()
    
    def _set_draw_norm(self, df=1):
        if self.dtype is np.float32:
            self._draw_max = 200
            self._draw_norm = 1
        if self.dtype is np.int8:
            #TODO: fix this!!!!
            gauss_limit = stats.chi2.ppf(0.999, df)
            self._draw_max = np.iinfo(np.int8).max
            self._draw_norm = self._draw_max/gauss_limit


    @property
    def subint(self):
        return self._subint

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

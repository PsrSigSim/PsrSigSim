
from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
import numpy as np

from ..utils.utils import make_quant

__all__ = ['Receiver']

class Receiver(object):
    """telescope reciever

    A :class:`Receiver` must be instatiated with either a callable response
    function or ``fcent`` and ``bandwidth`` to use a flat response.

    Optional Args:
        response (callable): frequency response function ("bandpass") of
            receiver.
        fcent (float): center frequency of reciever with flat response (MHz)
        bandwidth (float): bandwidth of reciever with flat response (MHz)
        Trec (float): reciever temperature (K) for radiometer noise level,
            default: ``35``
    """
    def __init__(self,
                 response=None, 
                 fcent=None, bandwidth=None,
                 Trec=35,
                 name=None):
        if response is None: 
            if fcent is None or bandwidth is None:
                msg = "specify EITHER response OR fcent and bandwidth"
                raise ValueError(msg)
            else:
                self._response = _flat_response(fcent, bandwidth)
        else:
            if fcent is not None or bandwidth is not None:
                msg = "specify EITHER response OR fcent and bandwidth"
                raise ValueError(msg)
            else:
                self._response = response

        self._Trec = make_quant(Trec, "K")
        self._name = name

    def __repr__(self):
        return "Receiver({:s})".format(self._name)

    @property
    def name(self):
        return self._name

    @property
    def Trec(self):
        return self._Trec

    @property
    def response(self):
        return self._response

def response_from_data(fs, values):
    """generate a callable response function from discrete data
    Data are interpolated.
    """
    raise NotImplementedError()

def _flat_response(fcent, bandwidth):
    """generate a callable, flat response function
    Required Args:
        fcent (float): central frequency (MHz)
        bandwidth (float): bandwidth (MHz)
    Returns:
        callable bandpass(f), where f is an array or scalar frequency (MHz)
    """
    fc = make_quant(fcent, 'MHz')
    bw = make_quant(bandwidth, 'MHz')
    fmin = fc - bw/2
    fmax = fc + bw/2
    return lambda f: np.heaviside(f-fmin, 0) * np.heaviside(fmax-f, 0)


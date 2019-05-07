
from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
import numpy as np

from .profiles import PulseProfile

class GaussPortrait(PulseProfile):
    """a sum of gaussian components in a 2-D pulse portrait

    Required Args:
        freqs (list of floats): central frequencies of sub-bands
    
    Optional Args:
        peak (float): center of gaussian in pulse phase, default: ``0.5``
        width (float): stdev of pulse in pulse phase, default: ``0.1``
        amp (float): amplitude of pulse relative to other pulses,
            default: ``1``
        Nphase (int): number of phase bins, default: ``256``
    """
    
    def __init__(self, peak=0.5, width=0.1, amp=1, Nphase=256, freqs=None):
        raise NotImplementedError()

class UserPortrait(PulseProfile):
    """user specified 2-D pulse portrait"""
    def __init__(self):
        raise NotImplementedError()

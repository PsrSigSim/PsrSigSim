
from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
import numpy as np

class Pulsar(object):
    """class for pulsars
    
    The minimal data to instatiate a pulsar is the period, intensity and
    pulse profile. The profile is supplied via a :class:`PulseProfile`-like
    object.

    Other data could be supplied via a `.par` file.

    """
    def __init__(self, profile=None):
        if profile is None:
            self._profile = GaussProfile()
        raise NotImplementedError()

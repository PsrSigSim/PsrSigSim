"""constants.py
A place to keep various constants needed for pulsar signals.
"""
from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
import numpy as np
import scipy as sp
from astropy import units as u
from astropy import constants as c
from .utils import make_quant

#constant used to be more consistent with PSRCHIVE
DM_K =  1.0/2.41e-4 * u.MHz**2 / u.pc * u.cm**3 * u.s

# Define the Kologorov Beta value for ISM scaling laws
KOLMOGOROV_BETA = 11.0/3

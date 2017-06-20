from __future__ import (absolute_import, division,
                        print_function, unicode_literals)

import numpy as np
import scipy as sp
from scipy import signal
import h5py
import math
import ism.py
from .ism.ISM import *
import PSS_utils as utils


class Simulate():
    """Class that simulates the pulsar signal
       Insert more information about the class
    """
    def __init__(self, signal_class,pulsar_class,ism_class):
        self.signal_class=signal_class
        self.DM=signal_class.MetaData.DM

    #Call function to make pulses
    def make_the_pulses():
        pulsar_class.gauss_template_beta()

    #Call function to disperse the pulses
    def dm_broaden():
        ism_class.disperse()

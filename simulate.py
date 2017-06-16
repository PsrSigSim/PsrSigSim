from __future__ import (absolute_import, division,
                        print_function, unicode_literals)

import numpy as np
import scipy as sp
from scipy import signal
import h5py
import math
import PSS_utils as utils
import sys
sys.path.append('/Users/nearclouding/Desktop/UTGRV_REU/nanograv_group')
import VersionZeroPointZero as PSS

#Initializing the signal
S1=PSS.signal(Nt=2000)

class Simulate():
    """Class that simulates the pulsar signal
       Insert more information about the class
    """
    def __init__(self, signal_class,pulsar_class,ism_class):
        self.signal_class=signal_class
        self.pulsar_class=pulsar_class
        self.ism_class=ism_class
        #NOTE need to give signal_class, pulsar_class, ism_class attributes?
        #What else goes in here?

    #Make the pules from initialized signal
    def make_the_pulses():
        self.P1=self.PSS.Pulsar(S1)
        P1.gauss_template_beta()
        P1.make_pulses_beta()

    #Disperse the pulses
    def dm_broaden():
        disperse(P1.make_pulses_beta)

    #TODO convolving
    """def convolve():
       want the try except stuff probably
    """

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Tests for `psrsigsim.pulsar` module."""

import pytest
import psrsigsim as pss
import psrsigsim.signal as sig
import psrsigsim.pulsar as psr
import os
import numpy as np

from psrsigsim.signal.fb_signal import FilterBankSignal
from psrsigsim.signal.bb_signal import BasebandSignal
from psrsigsim.pulsar.pulsar import Pulsar
from psrsigsim.pulsar.profiles import DataProfile
from psrsigsim.utils.utils import make_quant
from psrsigsim.ism.ism import ISM

@pytest.fixture
def fbsignal():
    """
    Fixture signal class
    """
    fbsig = FilterBankSignal(1400,400,Nsubband=2,\
                             sample_rate=1.0*2048*10**-6,\
                             sublen=0.5)
    return fbsig

@pytest.fixture
def fbsignal_nofold():
    """
    Fixture signal class
    """
    fbsig_nofold = FilterBankSignal(1400,400,Nsubband=2,\
                             sample_rate=1.0*2048*10**-6,\
                             sublen=None, fold = False)
    return fbsig_nofold

@pytest.fixture
def bbsignal():
    """
    Fixture signal class
    """
    bbsig = BasebandSignal(1400, 400,
                 sample_rate=1.0*2048*10**-6,
                 dtype=np.float32)
    bbsig._Nchan = 2
    return bbsig

@pytest.fixture
def j1713_profile():
    """
    Numpy array of J1713+0747 profile.
    """
    path = 'psrsigsim/data/J1713+0747_profile.npy'
    pr = DataProfile(np.load(path),phases=None,Nchan=2)
    return pr

@pytest.fixture
def ism():
    """
    Fixture ism class
    """
    return ISM()


def test_pulsar(j1713_profile):
    """
    Test pulsar class.
    """
    psr1 = Pulsar(1.0, 1.0, profiles=j1713_profile, name="J1713+0747")
    # initialize pulsar without profile
    psr2 = Pulsar(1.0, 1.0, profiles=None, name="J1713+0747")
    # Check functions and properties
    psr_names = psr1.__repr__()
    # check properties
    assert(psr1.period.value==1.0)
    assert(psr1.Smean.value==1.0)
    assert(psr1.name=="J1713+0747")
    assert(psr1.Profiles==j1713_profile)
    assert(psr1.specidx==0.0)
    assert(psr1.ref_freq==None)
    
def test_makepulses(fbsignal, fbsignal_nofold, bbsignal, j1713_profile):
    """
    Test make_pulses function.
    """
    # Define a pulsar
    psr1 = Pulsar(1.0, 1.0, profiles=j1713_profile, name="J1713+0747")
    # Test making pulses with a filterbank signal
    tobs = 2.0
    psr1.make_pulses(fbsignal, tobs)
    assert(psr1.ref_freq.value==1400.0)
    # Make pulses with sublen of None
    fbsignal._sublen = None
    psr1.make_pulses(fbsignal, tobs)
    # Make pulses without folded signal
    psr1.make_pulses(fbsignal_nofold, tobs)
    # Test making pulses with Baseband signal
    psr1.make_pulses(bbsignal, tobs)
    # Test error message
    with pytest.raises(NotImplementedError):
        fbsignal._sigtype = "badsignal"
        psr1.make_pulses(fbsignal, tobs)

def test_nullpulses(fbsignal, fbsignal_nofold, ism, j1713_profile):
    """
    Test pulse nulling function.
    """
    # Define a pulsar
    psr1 = Pulsar(1.0, 1.0, profiles=j1713_profile, name="J1713+0747")
    # Test making pulses with a filterbank signal
    tobs = 3.0
    psr1.make_pulses(fbsignal, tobs)
    dm =10.0
    ism.disperse(fbsignal, dm)
    # Null the pulses
    psr1.null(fbsignal, 0.34)
    # now check with no fold for a dispersed signal
    psr1.make_pulses(fbsignal_nofold, tobs)
    psr1.null(fbsignal_nofold, 0.34)
    # Check error
    with pytest.raises(NotImplementedError):
        psr1.null(fbsignal, 0.34, length = 1.0)

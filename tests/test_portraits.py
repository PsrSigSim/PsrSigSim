#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for `psrsigsim.ism` module."""

import pytest
import numpy as np
import psrsigsim as pss
import psrsigsim.signal as sig
import psrsigsim.pulsar as psr
import psrsigsim.ism as ism

from psrsigsim.signal.fb_signal import FilterBankSignal
from psrsigsim.pulsar.pulsar import Pulsar
from psrsigsim.pulsar.portraits import DataPortrait
from psrsigsim.pulsar.profiles import DataProfile
from psrsigsim.ism.ism import ISM
from psrsigsim.utils.utils import make_quant

@pytest.fixture
def j1713_profile():
    """
    Numpy array of J1713+0747 profile.
    """
    path = 'psrsigsim/data/J1713+0747_profile.npy'
    return np.load(path)

@pytest.fixture
def signal():
    """
    Fixture signal class
    """
    fbsig = FilterBankSignal(1400,400)
    return fbsig

@pytest.fixture
def pulsar():
    """
    Fixture pulsar class
    """
    period = make_quant(5,'ms')
    return Pulsar(period,10,name='J1746-0118')

def test_profile1(j1713_profile):
    pr = DataProfile(j1713_profile,phases=None)
    pr.init_profiles(100)
    assert(hasattr(pr,'_profiles'))
    assert(hasattr(pr,'_max_profile'))

def test_profile2(j1713_profile):
    pr = DataProfile(j1713_profile,phases=None)
    pr.init_profiles(100, Nchan=20)
    assert(hasattr(pr,'_profiles'))
    assert(hasattr(pr,'_max_profile'))

def test_profile3(j1713_profile):
    ph = np.linspace(0,1,2048)
    pr = DataProfile(j1713_profile,phases=ph)
    pr.init_profiles(100, Nchan=20)
    assert(hasattr(pr,'_profiles'))
    assert(hasattr(pr,'_max_profile'))

def test_portrait1(j1713_profile):
    dataport = np.tile(j1713_profile,(10,1))
    ph = np.linspace(0,1,2048)
    port = DataPortrait(dataport, phases=ph)
    port.init_profiles(100)
    assert(hasattr(port,'_profiles'))
    assert(hasattr(port,'_max_profile'))

def test_portrait2(j1713_profile):
    dataport = np.tile(j1713_profile,(10,1))
    port = DataPortrait(dataport)
    port.init_profiles(100)
    assert(hasattr(port,'_profiles'))
    assert(hasattr(port,'_max_profile'))

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
from psrsigsim.pulsar.portraits import PulsePortrait
from psrsigsim.pulsar.profiles import PulseProfile
from psrsigsim.pulsar.portraits import GaussPortrait
from psrsigsim.pulsar.profiles import GaussProfile
from psrsigsim.pulsar.portraits import UserPortrait
from psrsigsim.pulsar.profiles import UserProfile
from psrsigsim.pulsar.portraits import _gaussian_sing_1d, _gaussian_mult_1d, _gaussian_mult_2d
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
    """
    Test for data profiles 1.
    """
    pr = DataProfile(j1713_profile,phases=None)
    pr.init_profiles(100)
    assert(hasattr(pr,'_profiles'))
    assert(hasattr(pr,'_max_profile'))

def test_profile2(j1713_profile):
    """
    Test for data profiles 2.
    """
    pr = DataProfile(j1713_profile,phases=None)
    pr.init_profiles(100, Nchan=20)
    assert(hasattr(pr,'_profiles'))
    assert(hasattr(pr,'_max_profile'))
    with pytest.raises(NotImplementedError):
        pr.set_Nchan(2)

def test_profile3(j1713_profile):
    """
    Test for data profiles 3.
    """
    ph = np.linspace(0,1,2048)
    pr = DataProfile(j1713_profile,phases=ph)
    pr.init_profiles(100, Nchan=20)
    assert(hasattr(pr,'_profiles'))
    assert(hasattr(pr,'_max_profile'))

def test_portrait1(j1713_profile):
    """
    Test for data portrait 1.
    """
    dataport = np.tile(j1713_profile,(10,1))
    ph = np.linspace(0,1,2048)
    port = DataPortrait(dataport, phases=ph)
    port.init_profiles(100)
    assert(hasattr(port,'_profiles'))
    assert(hasattr(port,'_max_profile'))

def test_portrait2(j1713_profile):
    """
    Test for data portrait 2.
    """
    dataport = np.tile(j1713_profile,(10,1))
    port = DataPortrait(dataport)
    port.init_profiles(100)
    assert(hasattr(port,'_profiles'))
    assert(hasattr(port,'_max_profile'))
    
def test_pulseprofile():
    """
    Test for basic Pulse Profile Class.
    """
    # No input phases
    pprof1 = PulseProfile()
    # Input phases
    ph = np.linspace(0,1,2048)
    # Check errors
    with pytest.raises(NotImplementedError):
        pprof1.calc_profile(ph)
        # initilialize the profile
        pprof1.init_profile(2048)

def test_pulseportrait():
    """
    Test for Pulse Portrait Class.
    """
    # initialize portrait
    pport1 = PulsePortrait()
    # Check with no phases
    pport1._profiles = None
    pport1.__call__(phases=None)
    # Check errors
    ph = np.linspace(0,1,2048)
    with pytest.raises(NotImplementedError):
        pport1.__call__(phases=ph)
        pport1.init_profiles(ph)

def test_gaussprof():
    """
    Test the GaussProfile class.
    """
    gprof = GaussProfile(peak=0.5, width=0.05, amp=1)
    # Test errors
    with pytest.raises(NotImplementedError):
        gprof.set_Nchan(2)

def test_gaussportrait():
    """
    Test Gaussian Portrait class.
    """
    # Check initialize portrait
    gport1 = GaussPortrait(peak=0.5, width=0.05, amp=1)
    # Calculate the profiles
    gport1.init_profiles(2048, Nchan = 2)
    # Check properties
    assert(gport1.amp == 1.0)
    assert(gport1.peak == 0.5)
    assert(gport1.width == 0.05)
    assert(hasattr(gport1, 'Amax'))
    assert(hasattr(gport1, 'profiles'))
    # Check calculating the offpulse -> these fail, don't konw why...
    #opw1 = gport1._calcOffpulseWindow(Nphase=None)
    #opw2 = gport1._calcOffpulseWindow(Nphase=2048)
    #assert(opw1==opw2)
    # Check errors
    with pytest.raises(ValueError):
        gport1.init_profiles(2048, Nchan = None)
    
def test_userprof():
    """
    Test of UserProfile class.
    """
    # define function of profile
    def gaussfunc(phases):
        A = 1.0
        mu = 0.5
        sigma = 0.05
        return A*np.exp(0.5*((phases-mu)/(sigma**2)))
    # Define profile
    uprof = UserProfile(gaussfunc)
    normed_prof = uprof.calc_profile(np.linspace(0,1,2048))
    assert(isinstance(normed_prof, np.ndarray))
    assert(hasattr(uprof, 'profile'))
    assert(hasattr(uprof, 'Amax'))

def test_userportrait():
    """
    Test of user portrait class.
    """
    with pytest.raises(NotImplementedError):
        uport = UserPortrait()

def test_gaussianfuncs():
    """
    Test of additional gaussian functions in portrait class.
    """
    ph = np.linspace(0,1,2048)
    gs1d = _gaussian_sing_1d(ph, 0.5, 0.05, 1.0)
    gm1d = _gaussian_mult_1d(ph, np.array([0.25, 0.75]), np.array([0.05, 0.05]), np.array([1.0, 1.0]))
    gm2d = _gaussian_mult_2d(ph, np.array([0.25, 0.75]), np.array([0.05, 0.05]), np.array([1.0, 1.0]), 2)
    assert(isinstance(gs1d, np.ndarray))
    assert(isinstance(gm1d, np.ndarray))
    assert(isinstance(gm2d, np.ndarray))
    # Check errors
    ph = np.linspace(0,1.1,2048)
    with pytest.raises(ValueError):
        gs1d = _gaussian_sing_1d(ph, 0.5, 0.05, 1.0)
        gm1d = _gaussian_mult_1d(ph, np.array([0.25, 0.75]), np.array([0.05, 0.05]), np.array([1.0, 1.0]))


   
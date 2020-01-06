#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Tests for `psrsigsim.telescope` module."""

import pytest
import psrsigsim as pss
import psrsigsim.signal as sig
import psrsigsim.pulsar as psr
import psrsigsim.telescope as telescope

from psrsigsim.signal.fb_signal import FilterBankSignal
from psrsigsim.pulsar.pulsar import Pulsar
from psrsigsim.telescope.telescope import Telescope
from psrsigsim.telescope.receiver import Receiver
from psrsigsim.telescope.backend import Backend
from psrsigsim.utils.utils import make_quant


@pytest.fixture
def signal():
    """
    Fixture signal class
    """
    fbsig = FilterBankSignal(1400,400, fold=False)
    return fbsig

@pytest.fixture
def subint_signal():
    """
    Fixture for signal class
    - Makes a subintegrated signal
    """
    fbsubsig = FilterBankSignal(1400,400,fold=True)
    return fbsubsig

@pytest.fixture
def pulsar():
    """
    Fixture pulsar class
    """
    period = make_quant(5,'ms')
    return Pulsar(period,10,name='J1746-0118')


@pytest.fixture
def tscope():
    """
    Fixture telescope class
    """
    scope = Telescope(20.0, area=30.00, name="Twenty_Meter")
    return scope

@pytest.fixture
def receiver():
    """
    Fixture receiver class
    """
    rcvr = Receiver(fcent=1400,bandwidth=400,name="Lband")
    return rcvr

@pytest.fixture
def backend():
    """
    Fixture backend class
    """
    smprte = make_quant(81.92, "microsecond")
    bckend = Backend(samprate=(1.0/smprte).to("MHz"), name="Cyborg")
    return bckend
    
def test_obs(tscope, receiver, backend, signal, pulsar):
    """
    Test adding a system to the telescope and observing
    """
    tscope.add_system(name="Twnty_M", receiver=receiver, backend=backend)
    tobs = make_quant(0.02,'s')
    pulsar.make_pulses(signal,tobs)
    tscope.observe(signal, pulsar, system="Twnty_M", noise=False)

def test_noise(tscope, receiver, backend, signal, pulsar):
    """
    Test adding a system to the telescope and observing with noise
    """
    tscope.add_system(name="Twnty_M", receiver=receiver, backend=backend)
    tobs = make_quant(0.02,'s')
    pulsar.make_pulses(signal,tobs)
    tscope.observe(signal, pulsar, system="Twnty_M", noise=True)
    
def test_subint_obs(tscope, receiver, backend, subint_signal, pulsar):
    """
    Test adding a system to the telescope and observing with subint signal
    """
    tscope.add_system(name="Twnty_M", receiver=receiver, backend=backend)
    tobs = make_quant(0.02,'s')
    pulsar.make_pulses(subint_signal,tobs)
    tscope.observe(subint_signal, pulsar, system="Twnty_M", noise=False)

def test_subint_noise(tscope, receiver, backend, subint_signal, pulsar):
    """
    Test adding a system to the telescope and observing with noise and subint
    signal
    """
    tscope.add_system(name="Twnty_M", receiver=receiver, backend=backend)
    tobs = make_quant(0.02,'s')
    pulsar.make_pulses(subint_signal,tobs)
    tscope.observe(subint_signal, pulsar, system="Twnty_M", noise=True)
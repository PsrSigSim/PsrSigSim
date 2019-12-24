#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for `psrsigsim.ism` module."""

import pytest
import psrsigsim as pss
import psrsigsim.signal as sig
import psrsigsim.pulsar as psr
import psrsigsim.ism as ism

from psrsigsim.signal.fb_signal import FilterBankSignal
from psrsigsim.pulsar.pulsar import Pulsar
from psrsigsim.ism.ism import ISM
from psrsigsim.utils.utils import make_quant

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

@pytest.fixture
def ism():
    """
    Fixture ism class
    """
    return ISM()

def test_disperse(signal,pulsar,ism):
    """"""
    tobs = make_quant(5,'s')
    pulsar.make_pulses(signal,tobs)
    ism.disperse(signal,10)

def test_FDshift(signal,pulsar,ism):
    """"""
    tobs = make_quant(5,'s')
    pulsar.make_pulses(signal,tobs)
    ism.FD_shift(signal,[2e-4, -3e-4,7e-5])

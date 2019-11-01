#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Tests for `psrsigsim.io` module."""

import pytest
import psrsigsim as pss
import psrsigsim.signal as sig
import psrsigsim.pulsar as psr
import psrsigsim.io as io
import os

from psrsigsim.signal.fb_signal import FilterBankSignal
from psrsigsim.pulsar.pulsar import Pulsar
from psrsigsim.io.psrfits import PSRFITS
from psrsigsim.utils.utils import make_quant

@pytest.fixture
def pulsar():
    """
    Fixture pulsar class
    """
    period = make_quant(5,'ms')
    return Pulsar(period,10,name='J1746-0118')

@pytest.fixture
def PSRfits():
    """
    Fixture psrfits class
    """
    fitspath = "data/test.fits"
    tempfits = "data/B1855+09.L-wide.PUPPI.11y.x.sum.sm"
    return PSRFITS(path=fitspath, template=tempfits, fits_mode='copy')

def test_fitssig(PSRfits, pulsar):
    """
    Test getting a signal from a fits file and making pulses with it.
    """
    S = PSRfits.make_signal_from_psrfits()
    obslen = PSRfits.tsubint.value*PSRfits.nsubint
    pulsar.make_pulses(S, tobs = obslen )
    os.remove("data/test.fits")
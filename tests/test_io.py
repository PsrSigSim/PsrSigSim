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
from psrsigsim.io.txtfile import TxtFile

@pytest.fixture
def signal():
    """
    Fixture signal class
    """
    fbsig = FilterBankSignal(1400,400,Nsubband=2,\
                             sample_rate=186.49408124993144*2048*10**-6,\
                             sublen=0.5)
    return fbsig

@pytest.fixture
def pulsar():
    """
    Fixture pulsar class
    """
    F0 = 186.49408124993144
    period = make_quant(1.0/F0,'s')
    return Pulsar(period,10,name='J1746-0118')

@pytest.fixture
def PSRfits():
    """
    Fixture psrfits class
    """
    fitspath = "data/test.fits"
    tempfits = "data/B1855+09.L-wide.PUPPI.11y.x.sum.sm"
    return PSRFITS(path=fitspath, template=tempfits, fits_mode='copy')

@pytest.fixture
def TXTfile():
    """
    Fixture txtfile class
    """
    pdvpath = "data/test_pdv"
    return TxtFile(path=pdvpath)

def test_fitssig(PSRfits):
    """
    Test getting a signal from a fits file.
    """
    S = PSRfits.make_signal_from_psrfits()
    os.remove("data/test.fits")
    
def test_savesig(PSRfits, pulsar):
    """
    Test getting a signal from a fits file, making pulses with it, and save it.
    """
    S = PSRfits.make_signal_from_psrfits()
    obslen = PSRfits.tsubint.value*PSRfits.nsubint
    pulsar.make_pulses(S, tobs = obslen)
    PSRfits.save(S, pulsar)
    os.remove("data/test.fits")
    
def test_savephaseconnect(PSRfits, pulsar):
    """
    Test getting a signal from a fits file, making pulses with it, and save it,
    with all the phase connection functions.
    """
    S = PSRfits.make_signal_from_psrfits()
    obslen = PSRfits.tsubint.value*PSRfits.nsubint
    pulsar.make_pulses(S, tobs = obslen)
    parfile = "data/test_parfile.par"
    PSRfits.save(S, pulsar, phaseconnect = True, parfile = parfile, \
                 MJD_start = 55999.9861, inc_len = 0.0, ref_MJD = 56000.0, \
                 usePint = True)
    os.remove("data/test.fits")
    
def test_savephaseconnect_inc(PSRfits, pulsar):
    """
    Test getting a signal from a fits file, making pulses with it, and save it,
    with all the phase connection functions with an non-zerio increment length.
    """
    S = PSRfits.make_signal_from_psrfits()
    obslen = PSRfits.tsubint.value*PSRfits.nsubint
    pulsar.make_pulses(S, tobs = obslen)
    parfile = "data/test_parfile.par"
    PSRfits.save(S, pulsar, phaseconnect = True, parfile = parfile, \
                 MJD_start = 55999.9861+30.0, inc_len = 30.0, ref_MJD = 56000.0, \
                 usePint = True)
    os.remove("data/test.fits")
    
def test_savepdv(TXTfile, signal, pulsar):
    tobs = make_quant(1,'s')
    pulsar.make_pulses(signal,tobs)
    TXTfile.save_psrchive_pdv(signal, pulsar)
    os.remove("data/test_pdv_0.txt")


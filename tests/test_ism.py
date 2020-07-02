#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for `psrsigsim.ism` module."""

import pytest

from psrsigsim.signal.fb_signal import FilterBankSignal
from psrsigsim.signal.bb_signal import BasebandSignal
from psrsigsim.pulsar.pulsar import Pulsar
from psrsigsim.ism.ism import ISM
from psrsigsim.utils.utils import make_quant
import numpy as np

@pytest.fixture
def signal():
    """
    Fixture signal class
    """
    fbsig = FilterBankSignal(1400,400)
    return fbsig

@pytest.fixture
def bbsignal():
    """
    Fixture baseband signal class
    """
    bbsig = BasebandSignal(1400,400)
    return bbsig

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
    """
    Test disperse function with filterbank signal.
    """
    tobs = make_quant(5,'s')
    pulsar.make_pulses(signal,tobs)
    ism.disperse(signal,10)
    assert signal.dm.value==10
    with pytest.raises(ValueError):
        ism.disperse(signal,10)
        
#def test_bb_disperse(bbsignal,pulsar,ism):
#    """
#    Test disperse function with baseband signal.
#    """
#    tobs = make_quant(0.05,'s')
#    pulsar.make_pulses(bbsignal,tobs)
#    ism.disperse(bbsignal,10)
#    assert bbsignal.dm.value==10

def test_add_delay(signal,pulsar,ism):
    """
    Test add delay from dispersion class.
    """
    tobs = make_quant(5,'s')
    pulsar.make_pulses(signal,tobs)
    signal._delay = np.repeat(make_quant(1.0,'ms'), len(signal._dat_freq))
    ism.disperse(signal,10)
    assert signal.dm.value==10

def test_FDshift(signal,pulsar,ism):
    """
    Test FD shift function.
    """
    tobs = make_quant(5,'s')
    pulsar.make_pulses(signal,tobs)
    ism.FD_shift(signal,[2e-4, -3e-4,7e-5])

def test_scalinglaws(ism):
    """
    Test scaling law functions
    """
    nu_i = make_quant(1400.0, 'MHz')
    nu_f = make_quant(1200.0, 'MHz')
    # scale scintillation bandwidth
    dnu_d = make_quant(20.0, 'MHz')
    ism.scale_dnu_d(dnu_d,nu_i,nu_f,beta=5)
    ism.scale_dnu_d(dnu_d,nu_i,nu_f,beta=3)
    # scale scintillation timescale
    dt_d = make_quant(3600.0, 's')
    ism.scale_dt_d(dt_d,nu_i,nu_f,beta=5)
    ism.scale_dt_d(dt_d,nu_i,nu_f,beta=3)
    # scale scattering timescale
    tau_d = make_quant(5e-6, 's')
    ism.scale_tau_d(tau_d,nu_i,nu_f,beta=5)
    ism.scale_tau_d(tau_d,nu_i,nu_f,beta=3)

def test_scattershift(signal, pulsar, ism):
    """
    Test direct scattering shift in time.
    """
    tobs = make_quant(5,'s')
    pulsar.make_pulses(signal,tobs)
    ism.scatter_broaden(signal, 5e-6, 1400.0)

def test_scatterbroaden(signal, pulsar, ism):
    """
    Test pulsar scatter broadening with convolution.
    """
    tobs = make_quant(5,'s')
    ism.scatter_broaden(signal, 5e-6, 1400.0, convolve=True,pulsar=pulsar)
    pulsar.make_pulses(signal,tobs)

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Tests for `psrsigsim.utils` module."""

import pytest
import psrsigsim as pss
import psrsigsim.signal as sig
import psrsigsim.pulsar as psr
import psrsigsim.io as io
import os
import numpy as np

from psrsigsim.utils.utils import *
from astropy import units as u

from psrsigsim.signal.fb_signal import FilterBankSignal
from psrsigsim.pulsar.pulsar import Pulsar
from psrsigsim.ism.ism import ISM

"""
Still need to write tests for get_pint_models, function may be too specific
right now.
Still need to write test for text_search -> don't have files to test.
"""
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
def ism():
    """
    Fixture ism class
    """
    return ISM()

def test_shiftt():
    """
    Test the shift_t function.
    """
    y = np.arange(0,10)
    out = shift_t(y, 2, dt=1)
    out2 = shift_t(y, 2, dt=0.5)

def test_downsample():
    """
    Test downsampler function.
    """
    y = np.arange(0,10)
    out = down_sample(y, 2)
    assert(len(out) == 5)
    
def test_rebin():
    """
    Test rebin function.
    """
    y = np.arange(0,10)
    out = rebin(y, 5)
    assert(len(out) == 5)


def test_tophatwidth():
    """
    Test top had width function
    """
    width = top_hat_width(1.0, 1500.0, 10.0)
    assert(np.isclose(width, 0.024585528))
    
def test_SG():
    """
    Test the savitzky_golay function and errors.
    """
    y = np.arange(0,100)
    smoothed = savitzky_golay(y, 5, 2, deriv=0, rate=1)
    with pytest.raises(ValueError):
        smoothed = savitzky_golay(y, 'five', 'two', deriv=0, rate=1)
    with pytest.raises(TypeError):
        smoothed = savitzky_golay(y, 2, 2, deriv=0, rate=1)
        smoothed = savitzky_golay(y, 7, 6, deriv=0, rate=1)

def test_findnearest():
    """
    Test find_nearest function.
    """
    y = np.array([0,1,2,3,4,10,20,30])
    idx = find_nearest(y, 3)
    assert(idx==1)

def test_acf2d():
    """
    Test different implementations of the 2D-ACF function.
    """
    array2d = np.random.rand(3,2)
    retval = acf2d(array2d, speed='fast', mode='full', xlags=None, ylags=None)
    retval = acf2d(array2d, speed='slow', mode='full', xlags=None, ylags=None)
    retval = acf2d(array2d, speed='exact', mode='full', xlags=None, ylags=None)

def test_makequant():
    """
    Test make_quant function.
    """
    param = 3.0
    param = make_quant(param, u.second)
    param = make_quant(param, u.year)
    assert(hasattr(param, 'unit'))
    with pytest.raises(ValueError):
        param = make_quant(param, u.meter)

def test_makepar(signal, pulsar, ism):
    """
    Test make_par function.
    """
    tobs = make_quant(0.5,'s')
    pulsar.make_pulses(signal,tobs)
    ism.disperse(signal,10)
    make_par(signal, pulsar)
    os.remove("simpar.par")

def test_txtsearch():
    """
    Test text search function.
    """
    out = text_search("pull", ['header1', 'header2'], "data/txt_search_test.txt", header_line=1)
    out2 = text_search("pull", [1], "data/txt_search_test.txt", header_line=1)
    with pytest.raises(ValueError):
        out3 = text_search("push", [1], "data/txt_search_test.txt", header_line=1)
        out4 = text_search("nope", [1], "data/txt_search_test.txt", header_line=1)
    
def test_getpintmodel():
    """
    Test get_pint_model function.
    """
    psr_name = "J1910+1256"
    psr_file_path = "data/"
    m = get_pint_models(psr_name, psr_file_path)

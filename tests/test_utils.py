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

"""
Still need to write tests for get_pint_models, function may be too specific
right now.
Still need to write test for text_search -> don't have files to test.
"""

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
        smoothed = savitzky_golay(y, 'five', 2.0, deriv=0, rate=1)
    with pytest.raises(TypeError):
        smoothed = savitzky_golay(y, 2, 2, deriv=0, rate=1)
        smoothed = savitzky_golay(y, 3, 2, deriv=0, rate=1)

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

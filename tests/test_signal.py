#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Tests for `psrsigsim.signal` module."""

import pytest

import numpy as np

from psrsigsim.signal.fb_signal import FilterBankSignal
from psrsigsim.signal.bb_signal import BasebandSignal
from psrsigsim.signal.rf_signal import RFSignal
from psrsigsim.signal.signal import BaseSignal, Signal
from psrsigsim.utils.utils import make_quant

def test_basesignal():
    """
    Test BaseSignal class.
    """
    bs1 = BaseSignal(1400, 400, sample_rate=None, dtype=np.float32, Npols=1)
    # Test some functions
    sig_names = bs1.__repr__()
    bs1._Nchan = 2
    bs1.init_data(5)
    # Check some properties
    assert(bs1.Nchan == 2)
    assert(isinstance(bs1.data, np.ndarray))
    assert(bs1.sigtype == 'Signal')
    assert(bs1.fcent == 1400)
    assert(bs1.bw == 400)
    assert(bs1.samprate == None)
    assert(bs1.dtype == np.float32)
    assert(bs1.Npols == 1)
    assert(bs1.delay == None)
    assert(bs1.dm == None)
    assert(bs1.DM == None)

def test_basesignal_errors():
    """
    Test BaseSignal class errors.
    """
    with pytest.raises(ValueError):
        # Test invalid data type
        bs2 = BaseSignal(1400, -400, sample_rate=None, dtype=np.float64, Npols=1)
        # Test invalid polarizations
        bs3 = BaseSignal(1400, -400, sample_rate=None, dtype=np.float64, Npols=4)
    # Good signal
    bs1 = BaseSignal(1400, 400, sample_rate=None, dtype=np.float32, Npols=1)
    with pytest.raises(NotImplementedError):
        bs1.__add__(bs1)
        bs1._set_draw_norm()
        bs1.to_RF()
        bs1.to_Baseband()
        bs1.to_FilterBank()
        Signal()
        
def test_fbsignal():
    """
    Test Filterbank Signal class.
    """
    # standard instantiation
    fbsig1 = FilterBankSignal(1400,400,Nsubband=2,\
                             sample_rate=186.49408124993144*2048*10**-6,\
                             sublen=0.5)
    # minimal instantiation
    fbsig2 = FilterBankSignal(1400,-400,Nsubband=2,\
                             sample_rate=None,\
                             sublen=None, fold=False)
    # Check return of self
    fbsig_test = fbsig1.to_FilterBank()
    assert(fbsig_test == fbsig1)
    
    # Check _set_draw_norm function
    fbsig3 = FilterBankSignal(1400,-400,Nsubband=2,\
                             sample_rate=None,\
                             sublen=None, dtype=np.int8, fold=False)

def test_fbsignal_errs():
    """
    Test Filterbank Signal class errors.
    """
    # Test errors
    fbsig1 = FilterBankSignal(1400,400,Nsubband=2,\
                             sample_rate=186.49408124993144*2048*10**-6,\
                             sublen=0.5)
    with pytest.raises(NotImplementedError):
        fbsig1.to_RF()
        fbsig1.to_Baseband()

def test_bbsignal():
    """
    Test Baseband Signal class.
    """
    # standard instantiation
    bbsig1 = BasebandSignal(1400, 400,
                 sample_rate=186.49408124993144*2048*10**-6,
                 dtype=np.float32,
                 Nchan = 2)
    # minmal instantiation
    bbsig2 = BasebandSignal(1400, 400,
                 sample_rate=None,
                 dtype=np.float32,
                 Nchan = 2)
    # Check return of self
    bbsig_test = bbsig1.to_Baseband()
    assert(bbsig_test == bbsig1)

def test_bbsigna_errorsl():
    """
    Test Baseband Signal class errors.
    """
    # Test errors
    bbsig1 = BasebandSignal(1400, 400,
                 sample_rate=186.49408124993144*2048*10**-6,
                 dtype=np.float32,
                 Nchan = 2)
    with pytest.raises(NotImplementedError):
        bbsig1.to_RF()
        bbsig1.to_FilterBank()
    
def test_rfsignal():
    """
    Test RF Signal class.
    """
    # Standard intialization
    rfsig1 = RFSignal(1400, 400,
                 sample_rate=186.49408124993144*2048*10**-6,
                 dtype=np.float32)
    # minimal initialization
    rfsig2 = RFSignal(1400, 400,
                 sample_rate=None,
                 dtype=np.float32)
    # Check return of self
    rfsig_test = rfsig1.to_RF()
    assert(rfsig1 == rfsig_test)
    
def test_rfsignal_errors():
    """
    Test RF Signal class errors.
    """
    # Standard intialization
    rfsig1 = RFSignal(1400, 400,
                 sample_rate=186.49408124993144*2048*10**-6,
                 dtype=np.float32)
    with pytest.raises(NotImplementedError):
        rfsig1.to_Baseband()
        rfsig1.to_FilterBank()
    

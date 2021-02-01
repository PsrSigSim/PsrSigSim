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
from psrsigsim.io.file import BaseFile
from psrsigsim.ism.ism import ISM
import numpy as np

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

def test_basefile(signal):
    """
    Test the BaseFile class.
    """
    bf = BaseFile("data/B1855+09.L-wide.PUPPI.11y.x.sum.sm")
    assert(bf.path == "data/B1855+09.L-wide.PUPPI.11y.x.sum.sm")
    with pytest.raises(NotImplementedError):
        bf.save(signal)
    with pytest.raises(NotImplementedError):
        bf.load()
    with pytest.raises(NotImplementedError):
        bf.append()
    with pytest.raises(NotImplementedError):
        bf.to_txt()
    with pytest.raises(NotImplementedError):
        bf.to_psrfits()
    bf.path = "./"
    assert(bf.path == "./")

def test_fitssig(PSRfits):
    """
    Test getting a signal from a fits file.
    """
    S = PSRfits.make_signal_from_psrfits()
    os.remove("data/test.fits")

def test_getsigparams(PSRfits, signal, pulsar):
    """
    Test getting signal parameters from signal object
    """
    pulsar.make_pulses(signal, tobs = 1.0)
    PSRfits.get_signal_params(signal=signal)

@pytest.mark.filterwarnings('ignore::fitsio.FITSRuntimeWarning')
def test_savesig(PSRfits, pulsar):
    """
    Test getting a signal from a fits file, making pulses with it, and save it.
    """
    S = PSRfits.make_signal_from_psrfits()
    obslen = PSRfits.tsubint.value*PSRfits.nsubint
    pulsar.make_pulses(S, tobs = obslen)
    PSRfits.save(S, pulsar)
    os.remove("data/test.fits")

@pytest.mark.filterwarnings('ignore::fitsio.FITSRuntimeWarning')
@pytest.mark.filterwarnings('ignore::DeprecationWarning')
@pytest.mark.filterwarnings('ignore::astropy.utils.data.CacheMissingWarning')
def test_savephaseconnect(PSRfits, pulsar):
    """
    Test getting a signal from a fits file, making pulses with it, and save it,
    with all the phase connection functions.
    """
    S = PSRfits.make_signal_from_psrfits()
    obslen = PSRfits.tsubint.value*PSRfits.nsubint
    pulsar.make_pulses(S, tobs = obslen)
    parfile = "data/test_parfile.par"
    PSRfits.save(S, pulsar, parfile = None, \
                 MJD_start = 55999.9861, inc_len = 0.0, ref_MJD = 56000.0, \
                 usePint = True)
    os.remove("data/test.fits")

@pytest.mark.filterwarnings('ignore::fitsio.FITSRuntimeWarning')
@pytest.mark.filterwarnings('ignore::DeprecationWarning')
@pytest.mark.filterwarnings('ignore::astropy.utils.data.CacheMissingWarning')
def test_savephaseconnect_inc(PSRfits, pulsar):
    """
    Test getting a signal from a fits file, making pulses with it, and save it,
    with all the phase connection functions with an non-zerio increment length.
    """
    S = PSRfits.make_signal_from_psrfits()
    obslen = PSRfits.tsubint.value*PSRfits.nsubint
    pulsar.make_pulses(S, tobs = obslen)
    parfile = "data/test_parfile.par"
    PSRfits.save(S, pulsar, parfile = parfile, \
                 MJD_start = 55999.9861+30.0, inc_len = 30.0, ref_MJD = 56000.0, \
                 usePint = True)
    os.remove("data/test.fits")

def test_more_psrfits():
    """
    Additional PSRFITS file tests.
    Lines to test: 313-315 (can't test without multiple polarizations),
    326-328 (Can't test without using folded file),
    386-394 (Can't test without many more files),
    486, 506 (Can't test without SEARCH mode implemented)
    547-548 (Can't test without additional files),
    566-567 (can't test without different version of fitsio),
    """
    # Test initialization with no template, obs_mode SEARCH
    fitspath = "data/test.fits"
    tempfits = "data/B1855+09.L-wide.PUPPI.11y.x.sum.sm"
    pf = PSRFITS(path=fitspath, obs_mode = 'SEARCH', template=None)
    os.remove(fitspath)

def test_params(PSRfits, pulsar):
    """
    Test calling some attributes of the PSRFITS object.
    """
    S = PSRfits.make_signal_from_psrfits()
    obslen = PSRfits.tsubint.value*PSRfits.nsubint
    pulsar.make_pulses(S, tobs = obslen)
    assert(isinstance(PSRfits.tbin.value, float))
    assert(isinstance(PSRfits.chan_bw.value, float))
    assert(isinstance(PSRfits.stt_imjd.value, float))
    assert(isinstance(PSRfits.stt_smjd.value, float))

@pytest.mark.filterwarnings('ignore::DeprecationWarning')
@pytest.mark.filterwarnings('ignore::astropy.utils.data.CacheMissingWarning')
def test_gen_polyco(PSRfits):
    """
    Test generating the polyco dictionary.
    """
    parfile = "data/J1910+1256_NANOGrav_11yv1.gls.par"
    poly_dict = PSRfits._gen_polyco(parfile = parfile, MJD_start = 56000.0)
    with pytest.raises(NotImplementedError):
        poly_dict = PSRfits._gen_polyco(parfile = parfile, MJD_start = 56000.0, usePINT=False)

def test_gen_metadata(pulsar, PSRfits):
    """
    Test generating the subint and primary dictionary.
    """
    S = PSRfits.make_signal_from_psrfits()
    obslen = PSRfits.tsubint.value*PSRfits.nsubint
    pulsar.make_pulses(S, tobs = obslen)
    pri_dict, sub_dict = PSRfits._gen_metadata(S, pulsar, ref_MJD = 56000.0, inc_len = 0.0)
    pri_dict, sub_dict = PSRfits._gen_metadata(S, pulsar, ref_MJD = 56000.0, inc_len = 10.0)

@pytest.mark.filterwarnings('ignore::DeprecationWarning')
@pytest.mark.filterwarnings('ignore::astropy.utils.data.CacheMissingWarning')
def test_edit_header(pulsar):
    """
    Test editing the psrfits file header.
    """
    fitspath = "data/test.fits"
    tempfits = "data/B1855+09.L-wide.PUPPI.11y.x.sum.sm"
    PSRfits = PSRFITS(path=fitspath, template=tempfits, fits_mode='copy')
    S = PSRfits.make_signal_from_psrfits()
    obslen = 3.0
    pulsar.make_pulses(S, tobs = obslen)
    pri_dict, sub_dict = PSRfits._gen_metadata(S, pulsar, ref_MJD = 56000.0, inc_len = 0.0)
    sub_dict['POL_TYPE'] = 'AA+BB'
    # Then the bandwidth
    sub_dict['CHAN_BW'] = PSRfits.chan_bw.value
    # Get TSUBINT length
    sub_dict['TSUBINT'] = np.repeat(PSRfits.tsubint.value, PSRfits.nsubint)
    # Change TBIN
    sub_dict['TBIN'] = pulsar.period.value/PSRfits.nbin
    # Change TBIN
    sub_dict['DM'] = S.dm.value
    # Change TBIN
    sub_dict['NBIN'] = PSRfits.nbin
    parfile = "data/J1910+1256_NANOGrav_11yv1.gls.par"
    poly_dict = PSRfits._gen_polyco(parfile = parfile, MJD_start = 56000.0)
    PSRfits._nsblk = 1
    PSRfits._make_psrfits_pars_dict()
    PSRfits.copy_psrfit_BinTables()
    PSRfits._edit_psrfits_header(poly_dict, sub_dict, pri_dict)

def test_savepdv(TXTfile, signal, pulsar):
    """
    Test saving data in pdv format.
    """
    tobs = make_quant(1,'s')
    pulsar.make_pulses(signal,tobs)
    TXTfile.save_psrchive_pdv(signal, pulsar)
    os.remove("data/test_pdv_0.txt")


def test_moretxtfile(pulsar):
    """
    Test a few additional lines and properties of the txtfile class.
    Lines to Test: 59
    """
    fbsig = FilterBankSignal(1400,400,Nsubband=101,\
                             sample_rate=186.49408124993144*2048*10**-6,\
                             sublen=0.5)
    tobs = make_quant(1,'s')
    pulsar.make_pulses(fbsig,tobs)
    # Initialize without path
    tf = TxtFile(path=None)
    tf.save_psrchive_pdv(fbsig, pulsar)
    os.remove("PsrSigSim_Simulated_Pulsar.ar_1.txt")
    os.remove("PsrSigSim_Simulated_Pulsar.ar_2.txt")
    # Check attribute
    assert(tf.tbin == 1.0/fbsig.samprate)
    assert(tf.obsfreq == fbsig.fcent)
    assert(tf.chan_bw == fbsig.bw / fbsig.Nchan)

def test_notimplementedfuncs(PSRfits, signal):
    """
    Test functions that have not yet be implemented.
    """
    with pytest.raises(NotImplementedError):
        parfile = "data/test_parfile.par"
        PSRfits._gen_polyco(parfile, 56000.0, usePINT=False)
        PSRfits.append(signal)
        PSRfits.load()
        PSRfits.to_txt()
        PSRfits.to_psrfits()
        PSRfits.set_sky_info()
        PSRfits._calc_psrfits_dims(signal)
    os.remove("data/test.fits")

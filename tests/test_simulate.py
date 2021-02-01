#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pytest
import os
import numpy as np
import glob

from psrsigsim.signal.fb_signal import FilterBankSignal
from psrsigsim.pulsar.pulsar import Pulsar
from psrsigsim.pulsar.portraits import DataPortrait
from psrsigsim.pulsar.profiles import DataProfile
from psrsigsim.ism.ism import ISM
from psrsigsim.telescope.telescope import Telescope
from psrsigsim.telescope.receiver import Receiver
from psrsigsim.telescope.backend import Backend
from psrsigsim.io.psrfits import PSRFITS
from psrsigsim.utils.utils import make_quant
from psrsigsim.io.txtfile import TxtFile
from psrsigsim.simulate.simulate import Simulation


@pytest.fixture
def j1713_profile():
    """
    Numpy array of J1713+0747 profile.
    """
    path = 'psrsigsim/data/J1713+0747_profile.npy'
    return np.load(path)

@pytest.fixture
def PSRfits():
    """
    Fixture psrfits class
    """
    fitspath = "data/test.fits"
    tempfits = "data/B1855+09.L-wide.PUPPI.11y.x.sum.sm"
    return PSRFITS(path=fitspath, template=tempfits, fits_mode='copy')

@pytest.fixture
def param_dict():
    """
    Fixture parameter dictionary.
    """
    pdict = {'fcent' : 430,
             'bandwidth' : 100,
             'sample_rate' : 1.5625,
             'dtype' : np.float32,
             'Npols' : 1,
             'Nchan' : 64,
             'sublen' : 2.0,
             'fold' : True,
             'period' : 1.0,
             'Smean' : 1.0,
             'profiles' : [0.5, 0.5, 1.0], # Gaussian
             'tobs' : 4.0,
             'name' : 'J0000+0000',
             'dm' : 10.0,
             'tau_d' : 50e-9,
             'tau_d_ref_f' : 1500.0,
             'aperture' : 100.0,
             'area' : 5500.0,
             'Tsys' : 35.0,
             'tscope_name' : "TestScope",
             'system_name' : "TestSys",
             'rcvr_fcent' : 430,
             'rcvr_bw' : 100,
             'rcvr_name' : "TestRCVR",
             'backend_samprate' : 1.5625,
             'backend_name' : "TestBack",
             'tempfile' : None,
             'parfile' : None,
            }
    return pdict

@pytest.fixture
def simulation():
    """
    Fixture Simulation class. Cannot be the only simulation tested.
    """
    sim = Simulation(fcent = 430,
                 bandwidth = 100,
                 sample_rate = 1.0*2048*10**-6,
                 dtype = np.float32,
                 Npols = 1,
                 Nchan = 64,
                 sublen = 2.0,
                 fold = True,
                 period = 1.0,
                 Smean = 1.0,
                 profiles = None,
                 tobs = 4.0,
                 name = 'J0000+0000',
                 dm = 10.0,
                 tau_d = 50e-9,
                 tau_d_ref_f = 1500.0,
                 aperture = 100.0,
                 area = 5500.0,
                 Tsys = 35.0,
                 tscope_name = "TestScope",
                 system_name = "TestSys",
                 rcvr_fcent = 430,
                 rcvr_bw = 100,
                 rcvr_name ="TestRCVR",
                 backend_samprate = 1.5625,
                 backend_name = "TestBack",
                 tempfile = "data/B1855+09.L-wide.PUPPI.11y.x.sum.sm",
                 parfile = None,
                 psrdict = None)
    return sim

def test_initsim(param_dict):
    """
    Test initializing the simulation from dictionary, parfile
    """
    sim = Simulation(psrdict = param_dict)
    with pytest.raises(NotImplementedError):
        sim2 = Simulation(parfile = "testpar.par")

def test_initsig(simulation):
    """
    Test init_signal function.
    """
    # Test from input params
    simulation.init_signal()
    # Test from template file
    simulation.init_signal(from_template = True)

def test_initprof(simulation, j1713_profile):
    """
    Test init_profile function.
    """
    # Test no input
    simulation.init_profile()
    # Test function input
    with pytest.raises(NotImplementedError):
        def gprof(x, p0):
            return p0[0]* np.exp(-0.5*((x-p0[1])/(p0[2]))**2)
        simulation._profiles = gprof
        simulation.init_profile()
    # Test Gaussian as input
    simulation._profiles = [0.5, 0.5, 1.0]
    simulation.init_profile()
    # Test data array as input
    simulation._profiles = j1713_profile
    simulation.init_profile()
    # Test array that's not long enough
    with pytest.raises(RuntimeError):
        simulation._profiles = [0.5, 0.5]
        simulation.init_profile()
    # Test profile class as input
    pr = DataProfile(j1713_profile,phases=None)
    print(type(pr), pr)
    simulation._profiles = pr
    simulation.init_profile()
    
def test_initpsr(simulation):
    """
    Test init_pulsar function.
    """
    simulation.init_pulsar()

def test_initism(simulation):
    """
    Test init_ism function.
    """
    simulation.init_ism()

def test_inittscope(simulation):
    """
    Test init_telescope function.
    """
    # Test init GBT
    simulation._tscope_name = "GBT"
    simulation.init_telescope()
    # Test init Arecibo
    simulation._tscope_name = "Arecibo"
    simulation.init_telescope()
    # Test input telescope
    simulation._tscope_name = "TestScope"
    simulation.init_telescope()
    # Test list of systems for telescope
    simulation._system_name = ["Sys1", "Sys2"]
    simulation._rcvr_fcent = [430, 800]
    simulation._rcvr_bw = [100, 200]
    simulation._rcvr_name = ["R1", "R2"]
    simulation._backend_samprate = [1.5625, 12.5]
    simulation._backend_name = ["B1", "B2"]
    simulation.init_telescope()
    # And the catch with multiple systems
    with pytest.raises(RuntimeError):
        simulation._backend_name = ["B1", "B2", "B3"]
        simulation.init_telescope()

def test_simulate(simulation):
    """
    Test simulate function.
    """
    simulation.simulate()

@pytest.mark.filterwarnings('ignore::fitsio.FITSRuntimeWarning')
def test_savesim(simulation, PSRfits):
    """
    Test save simulation function.
    """
    simulation._Nchan = 1
    simulation._tobs = 2.0
    #S = PSRfits.make_signal_from_psrfits()
    #simulation._tobs = PSRfits.tsubint.value*PSRfits.nsubint
    simulation.simulate(from_template = True)
    # Try pdv format
    simulation.save_simulation(out_format = "pdv")
    # Try psrfits format
    simulation.save_simulation(out_format = "psrfits")
    os.remove("sim_fits.fits")
    dfs = glob.glob("simfits*")
    for df in dfs:
        os.remove(df)
    # Try psrfits with runtime error
    # Try wrong output file type
    with pytest.raises(RuntimeError):
        simulation.save_simulation(out_format = "wrong_fmt")
        simulation._tempfile = None
        simulation.save_simulation(out_format = "psrfits")


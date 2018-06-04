#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for `psrsigsim` package."""

import pytest
import unittest

import psrsigsim as PSS


@pytest.fixture
def response():
    """Sample pytest fixture.

    See more at: http://doc.pytest.org/en/latest/fixture.html
    """
    # import requests
    # return requests.get('https://github.com/audreyr/cookiecutter-pypackage')
    pass


def test_content(response):
    """Sample pytest test function with the pytest fixture as an argument."""
    # from bs4 import BeautifulSoup
    # assert 'GitHub' in BeautifulSoup(response.content).title.string
    pass


class test_PSS_classes(unittest.TestCase):

    def test_signal(self):
        """Test Signal Instantiation"""
        PSS.Signal()

    def test_pulsar(self):
        """Test Pulsar Class"""
        S1 = PSS.Signal()
        P1 = PSS.Pulsar(S1)
        P1.make_pulses()
        
    def test_simulate(self)
        test_dict = {'f0':1400 , 'bw': 400, 'Nf' : 50, 'data_type': 'int8', 'SignalType': 'intensity', 'freq_band': 1400,\
                     'ObsTime': 10, 'flux' : 3, 'f_samp': 4, 'radiometer_noise' : False,'tau_scatter':6e-08,\
                     'to_DM_Broaden': False}
        s1 =  PSS.Simulation(psr =  'J1713+0747' , sim_telescope= 'GBT',sim_ism= True, sim_scint= True, sim_dict = test_dict,)
        s1.simulate()
        test_dict_2 = {'F0': 218,'dm':15.921200,'scint_bw':15.6, 'scint_timescale':2630}
        test_dict.update(test_dict_2)
        s2 = PSS.Simulation(psr =  None , sim_telescope= 'GBT',sim_ism= None, sim_scint= None, sim_dict = test_dict)
        s2.init_signal()
        s2.init_pulsar()
        s2.init_ism()
        s2.init_telecope()
        s2.simulate()

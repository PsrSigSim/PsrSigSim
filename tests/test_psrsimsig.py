#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for `psrsimsig` package."""

import pytest
import unittest

import psrsimsig as PSS


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
        self.S1 = PSS.Signal()

    def test_pulsar(self):
        """Test Pulsar Class"""
        self.P1 = PSS.Pulsar(self.S1)
        self.P1.make_pulses()

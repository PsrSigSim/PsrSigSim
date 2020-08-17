#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for tutorial notebooks."""

import pytest
import subprocess as sp


def test_tutorial_1():
    sp.call('jupyter nbconvert --to notebook --inplace --execute ../docs/_static/tutorial_1.ipynb',
            shell=True)

def test_tutorial_2():
    sp.call('jupyter nbconvert --to notebook --inplace --execute ../docs/_static/tutorial_2.ipynb',
            shell=True)

def test_pulse_profiles_tutorial_3():
    sp.call('jupyter nbconvert --to notebook --inplace --execute ../docs/_static/pulse_profiles_tutorial_3.ipynb',
            shell=True)

def test_ism_options_tutorial_4():
    sp.call('jupyter nbconvert --to notebook --inplace --execute ../docs/_static/ism_options_tutorial_4.ipynb',
            shell=True)

def test_telescopes_tutorial_5():
    sp.call('jupyter nbconvert --to notebook --inplace --execute ../docs/_static/telescopes_tutorial_5.ipynb',
            shell=True)

def test_simulate_tutorial():
    sp.call('jupyter nbconvert --to notebook --inplace --execute ../docs/_static/simulate_tutorial.ipynb',
            shell=True)

def test_pulse_nulling_example():
    sp.call('jupyter nbconvert --to notebook --inplace --execute ../docs/_static/pulse_nulling_example.ipynb',
            shell=True)

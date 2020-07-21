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

def test_simulate_tutorial():
    sp.call('jupyter nbconvert --to notebook --inplace --execute ../docs/_static/simulate_tutorial.ipynb',
            shell=True)

def test_pulse_nulling_example():
    sp.call('jupyter nbconvert --to notebook --inplace --execute ../docs/_static/pulse_nulling_example.ipynb',
            shell=True)

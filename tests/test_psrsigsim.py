#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Tests for `psrsigsim` package."""

import pytest
import psrsigsim as pss
import psrsigsim.signal as sig
import psrsigsim.pulsar as psr
import psrsigsim.ism as ism

def test_pss():
    assert hasattr(pss,'__version__')
    assert hasattr(pss,'signal')

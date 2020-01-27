# -*- coding: utf-8 -*-
# encoding=utf8

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
import numpy as np

class BaseFile(object):
    """
    Base class for making files for the PsrSigSim Signal object.
    """
    _path = None
    _signal = None
    _file = None

    def __init__(self, path=None):
        """
        """
        self._path = path


    def save(self, signal):
        """Save PSS signal file to disk.
        Must be implemented in subclass!"""
        raise NotImplementedError()

    def append(self):
        """Method for append data to an already existing PSS signal file.
        Must be implemented in subclass!"""
        raise NotImplementedError()

    def load(self):
        """Method for loading saved PSS signal files.
        Must be implemented in subclass!"""
        raise NotImplementedError()

    def to_txt(self): #Can delete if not planning on implementing
        """Convert file to txt file.
        Must be implemented in subclass!"""
        raise NotImplementedError()

    def to_psrfits(self):
        """Convert file to PSRFITS file.
        Must be implemented in subclass!"""
        raise NotImplementedError()

    @property
    def path(self):
        return self._path

    @path.setter
    def path(self, value):
        self._path = value

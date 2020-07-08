# -*- coding: utf-8 -*-
# encoding=utf8

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
import numpy as np
from ..utils import make_quant
from .file import BaseFile
from astropy import log

class TxtFile(BaseFile):
    """A class for saving PsrSigSim signals as text files.
    Multiple different text file data types may be supported, but currently
    only the PSRCHIVE pdv function output is supported.

    Parameters
    ----------

    path: name and path of new text file that will be saved
    """
    _path = None
    _signal = None
    _file = None

    def __init__(self, path=None):
        """
        """
        self._path = path
        self._tbin = None
        self._nbin = None
        self._nchan = None
        self._npol = None
        self._nrows = None
        self._tsubint = None
        self._chan_bw = None
        self._obsbw = None
        self._obsfreq = None

    def save_psrchive_pdv(self, signal, pulsar):
        """
        Function to save simulated data in the same format as the PSRCHIVE pdf function.
        To avoid large file sizes, every hundred pulses the data will be saved as a text file.
        Currently one a single polarization (total intensity) is supported. Inputs are:

        Parameters
        ----------

        signal [class] : signal class object, currently only filterbank is supported
        pulsar [class] : pulsar class object
        filename [string] : desired name of source/output file. Output files will be saved as
                            'filename'_#.txt, where # is the chronological number of the
                            files being saved.
        """
        # Get the signal parameters
        self._get_signal_params(signal, pulsar)
        # Check if there's a save path
        if self.path == None:
            self._path = "PsrSigSim_Simulated_Pulsar.ar"
        # We first need to make the header
        rms = np.sqrt((1.0/len(signal.data))*np.sum((signal.data)**2))
        header = "# File: %s Src: %s Nsub: %s Nch: %s Npol: %s Nbin: %s RMS: %s \n" % \
                        (self.path, pulsar.name, str(self.nrows), str(self.nchan), \
                         str(self.npol), str(self.nbin), str(rms))
        # Now we make a list of the lines to save
        lines = [header]
        if self.npol != 1:
            log.warning("Only saving total intensity, multiple polarizations not yet implemented")
        # We need to dump the data every so often to avoid files that are too big
        dump_val = 0
        # Now we need to loop through each subint, frequency channel, and polarization
        for ii in range(self.nrows):
            mjd_mid = 56000.0+(ii+1)*(self.tsubint.to('day').value)/2.0
            for ff in range(self.nchan):
                freq = signal.dat_freq[ff].value
                pulse_header = "# MJD(mid): %s Tsub: %s Freq: %s BW: %s \n" \
                            % (mjd_mid, self.tsubint.value, freq, self.obsbw.value/self.nchan)
                lines.append(pulse_header)
                for bb in range(self.nbin):
                    lines.append("%s %s %s %s \n" % (ii, ff, bb, signal.data[ff,bb]))
                dump_val += 1
            # Now save the lines in a text file
            if dump_val >= 100:
                file_num = str(int(np.floor(dump_val/100)))
                with open(self.path+"_%s.txt" % (file_num), 'w') as pdv_file:
                    pdv_file.writelines(lines)
                    pdv_file.close()
                lines = [header]
        # and then we want to dump the last few lines
        file_num = str(int(np.floor(dump_val/100)))
        with open(self.path+"_%s.txt" % (file_num), 'w') as pdv_file:
            pdv_file.writelines(lines)
            pdv_file.close()

    def _get_signal_params(self, signal, pulsar):
        """
        Retrieve the various parameters needed to save the file from
        a PSS signal class.
        """
        # get parameters from signal class
        self.nchan = signal.Nchan
        self.tbin = 1.0/signal.samprate
        self.nbin = int((signal.samprate * pulsar.period).decompose())
        self.npol = signal.Npols
        self.nrows = signal.nsub
        self.obsfreq = signal.fcent
        self.obsbw = signal.bw
        self.chan_bw = signal.bw / signal.Nchan
        self.tsubint = signal.sublen # length of subint in seconds
        self.nsubint = self.nrows

      #### Define various PSRFITS parameters
    @property
    def tbin(self):
        return self._tbin

    @tbin.setter
    def tbin(self, value):
        self._tbin = make_quant(value,'s')

    @property
    def npol(self):
        return self._npol

    @npol.setter
    def npol(self, value):
        self._npol = value

    @property
    def nchan(self):
        return self._nchan

    @nchan.setter
    def nchan(self, value):
        self._nchan = value

    @property
    def nbin(self):
        return self._nbin

    @nbin.setter
    def nbin(self, value):
        self._nbin = value

    @property
    def nrows(self):
        return self._nrows

    @nrows.setter
    def nrows(self, value):
        self._nrows = value

    @property
    def obsfreq(self):
        return self._obsfreq

    @obsfreq.setter
    def obsfreq(self, value):
        self._obsfreq = make_quant(value,'MHz')

    @property
    def obsbw(self):
        return self._obsbw

    @obsbw.setter
    def obsbw(self, value):
        self._obsbw = make_quant(value,'MHz')

    @property
    def chan_bw(self):
        return self._chan_bw

    @chan_bw.setter
    def chan_bw(self, value):
        self._chan_bw = make_quant(value,'MHz')

    @property
    def tsubint(self):
        return self._tsubint

    @tsubint.setter
    def tsubint(self, value):
        self._tsubint = make_quant(value, 'second')

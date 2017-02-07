"""pulsar.py
a starting point for the Pulsar class.
"""

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
import numpy as np
import h5py

class MetaData(object):
    """the MetaData class to contain information about the signal
    """
    def __init__(self):
        self.f0 = None # central freq (MHz)
        self.bw = None # bandwidth (MHz)
        self.Nf = None # number of frequency bins
        self.Nt = None # number of time/phase bins
        self.TotTime = None #Total time in milliseconds

    def AddInfo(self,Info):
        """Function to Add Information into the metadata from a dictionary
            Since each new module will have a dictionary of terms to
            add in to the metadata
        """
        for ii,jj in Info.items():
            setattr(self, ii, jj)

        """TODO Add error message if someone trys to overlap parameters
        """


class Signal(object):
    """the signal class
    """
    def __init__(self, f0=400, bw=100, Nf=20, Nt=400, TotTime=200):
        """initialize Signal(), executed at assignment of new instance
        @param f0 -- central frequency (MHz)
        @param bw -- bandwidth (MHz)
        @param Nf -- number of freq. bins
        @param Nt -- number of phase bins
        """
        self.MetaData = MetaData()
        self.f0 = f0 # (MHz)
        self.bw = bw # (MHz)
        self.Nf = Nf # freq bins
        if Nt%2 == 0: # Make signal even in length (for FFTs)
            self.Nt = Nt # phase bins
        else: self.Nt = Nt + 1
        self.TotTime = TotTime #Total time in milliseconds
        if self.Nt*self.Nf > 500000: #Limits the array size to 2.048 GB
            SignalPath = "signal.hdf5"
            SignalFile = h5py.File(SignalPath)
            self.signal = SignalFile.create_dataset(None, (self.Nf, self.Nt), dtype='int16')
            #self.signal = np.memmap(SignalPath, dtype = 'float32', mode = 'w+', shape = (self.Nf, self.Nt))
        else:
            self.signal = np.zeros((self.Nf, self.Nt))


    """set the signal parameters as properties and assign them to the MetaData
    """
    @property

    def f0(self):
        return self._f0

    @f0.setter
    def f0(self, value):
        self._f0 = value
        self.MetaData.f0 = self._f0

    @property

    def bw(self):
        return self._bw

    @bw.setter
    def bw(self, value):
        self._bw = value
        self.MetaData.bw = self._bw

    @property

    def Nf(self):
        return self._Nf

    @Nf.setter
    def Nf(self, value):
        self._Nf = value
        self.MetaData.Nf = self._Nf

    @property

    def Nt(self):
        return self._Nt

    @Nt.setter
    def Nt(self, value):
        self._Nt = value
        self.MetaData.Nt = self._Nt

    @property

    def TotTime(self):
        return self._TotTime

    @TotTime.setter
    def TotTime(self, value):
        self._TotTime = value
        self.MetaData.TotTime = self._TotTime

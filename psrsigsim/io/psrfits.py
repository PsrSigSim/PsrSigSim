# -*- coding: utf-8 -*-
# encoding=utf8

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
import numpy as np
from packaging import version
import fitsio
import pdat
from .file import BaseFile
from ..utils import make_quant
from ..signal import Signal
from ..signal import FilterBankSignal

__all__ = ["PSRFITS"]

class PSRFITS(BaseFile):
    """A class for saving PsrSigSim signals as PSRFITS standard files.
    path: name and path of new psrfits file that will be saved
    obs_mode: what type of observation is the data, SEARCH, PSR, etc. 
    template: the path and name of the template fits file that will be loaded
    copy_template: Does nothing?
    fits_mode: how we want to save the data, right now just 'copy' is valid
    """

    def __init__(self, path=None, obs_mode=None, template=None,
                 copy_template=False,fits_mode='copy'):

        self._tbin = None
        self._nbin = None
        self._nsblk = None
        self._nchan = None
        self._npol = None
        self._nrows = None
        self._tsubint = None
        self._chan_bw = None
        self._obsbw = None
        self._obsfreq = None

        self._fits_mode = fits_mode

        if template is None:
            template_path = False
        else:
            template_path = template

        super().__init__(path=path)

        self.file = pdat.psrfits(psrfits_path=path, mode='rw',
                                 from_template=template_path,
                                 obs_mode=obs_mode, verbose=False)
        if obs_mode is None:
            self.obs_mode = self.file.obs_mode
        else:
            self.obs_mode = obs_mode

        #Lists of needed parameters from the PSRFITS file.
        self.pfit_pars = {'PRIMARY':['TELESCOP',
                                     'FRONTEND',
                                     'BACKEND',
                                     'OBS_MODE',
                                     'OBSFREQ',
                                     'OBSBW',
                                     'OBSNCHAN',
                                     'FD_POLN',
                                     'STT_IMJD',
                                     'STT_SMJD',
                                     'STT_OFFS'],
                          'SUBINT':['TBIN',
                                    'NAXIS',
                                    'NAXIS1',
                                    'NAXIS2',
                                    'NCHAN',
                                    'POL_TYPE',
                                    'NPOL',
                                    'NBIN',
                                    'NBITS',
                                    'CHAN_BW',
                                    'NSBLK',
                                    'DAT_SCL',
                                    'DAT_OFFS',
                                    'DAT_WTS',
                                    "TSUBINT",], # TTYPE2 is the subint length
                            'PSRPARAM':[
                                    ]
                           }

        #Might not need this, since we can calculate from other terms.
        if self.obs_mode=='SEARCH':
            self.pfit_pars['SUBINT'].append('TDIM17')
        elif self.obs_mode=='PSR':
            self.pfit_pars['SUBINT'].append('TDIM20')
            self.pfit_pars['PSRPARAM'].append('F')
            # Need both there depending on the template file. Only one will work.
            self.pfit_pars['PSRPARAM'].append('F0')

    def save(self, signal):
        """Save PSS signal file to disk.
        """
        """
        # May come back to this later...
        if self._fits_mode == 'copy':
            pass
        elif self._fits_mode == 'manual':
            pass
        elif self._fits_mode == 'auto':
            pass
        """
        # We need to appropriatly shape the signal for the fits file
        stop = self.nbin*self.nsubint
        sim_sig = signal.data[:,:stop].astype('>i2')
        # out arrays
        Out = np.zeros((self.nsubint, self.npol, self.nchan, self.nbin))
        print(np.shape(Out))
        # We assign the data in the appropriate shape
        for ii in range(self.nsubint):
            idx0 = 0 + ii*2048
            idxF = idx0 + 2048
            Out[ii,0,:,:] = sim_sig[:,idx0:idxF]
        
        self.copy_psrfit_BinTables()
        # We can currently only make total intensity data
        self.file.set_draft_header('SUBINT',{'POL_TYPE':'AA+BB'})
        """IMPORTANT NOTE: Currently phase connection is not implimented here! Still need to do this."""
        for ii in range(self.nsubint):
            self.file.HDU_drafts['SUBINT'][ii]['DATA'] = Out[ii,0,:,:]
            self.file.HDU_drafts['SUBINT'][ii]['DAT_FREQ'] = signal.dat_freq.value
            # Get the shapes of the wieghts, scales, and offs arrays, assumes we want to reset these to all be equal
            if len(np.shape(self.file.fits_template[4][0]['DAT_SCL'])) != 1:
                scale_shape = np.shape(self.file.fits_template[4][ii]['DAT_SCL'][:,:self.nchan*self.npol])
                offs_shape = np.shape(self.file.fits_template[4][ii]['DAT_OFFS'][:,:self.nchan*self.npol])
                weight_shape = np.shape(self.file.fits_template[4][ii]['DAT_WTS'])
            else:
                scale_shape = np.shape(self.file.fits_template[4][ii]['DAT_SCL'][:])
                offs_shape = np.shape(self.file.fits_template[4][ii]['DAT_OFFS'][:])
                weight_shape = np.shape(self.file.fits_template[4][ii]['DAT_WTS'])
            # Now assign the values
            #print(scale_shape, offs_shape, weight_shape)
            self.file.HDU_drafts['SUBINT'][ii]['DAT_SCL'] = np.ones(scale_shape)
            self.file.HDU_drafts['SUBINT'][ii]['DAT_OFFS'] = np.zeros(offs_shape)
            self.file.HDU_drafts['SUBINT'][ii]['DAT_WTS'] = np.ones(weight_shape)
        # Now we actually write out the files
        self.file.write_psrfits(hdr_from_draft=True)
        # Close the file so it doesn't take up memory or get confused with another file. 
        self.file.close()
        print("Finished writing and saving the file")
        

    def append(self, signal):
        """Method for appending data to an already existing PSS signal file.
        """
        raise NotImplementedError()

    def load(self):
        """Method for loading saved PSS signal files. These files will have an
        additional BinTable extension, 'PSRSIGSIM', that only contains a header
        with the various PsrSigSim parameters written for references.
        """
        raise NotImplementedError()

    def make_signal_from_psrfits(self):
        """Method to make a signal from the PSRFITS file given as the template.
        For subintegrated data will assume the initial period is the pulsar 
        period given in the PSRPARAM header.
        
        TODO: Currently does not support generating 'SEARCH' mode data from
            a psrfits file

        Parameters
        ----------

        Returns
        -------

        psrsigsim.Signal
        """
        self._fits_mode = 'copy'
        self._get_signal_params()
        
        if self.obs_mode == 'PSR':
            ObsTime = self.tsubint*self.nrows
            # Get correct period value from template fits options
            if self.pfit_dict['F'] is None and self.pfit_dict['F0'] is not None:
                s_rate = self.pfit_dict['F0']*self.nbin*10**-6 # in MHz
            elif self.pfit_dict['F'] is not None and self.pfit_dict['F0'] is None:
                s_rate = self.pfit_dict['F']*self.nbin*10**-6 # in MHz
            elif self.pfit_dict['F'] is not None and self.pfit_dict['F0'] is not None:
                s_rate = self.pfit_dict['F0']*self.nbin*10**-6 # in MHz
            else:
                msg = "No pulsar frequency defined in input fits file."
                raise ValueError(msg)        
        else:
            ObsTime = self.tbin*self.nbin*self.nsblk*self.nrows
            s_rate = (1/self.tbin).to('MHz').value

        S = FilterBankSignal(fcent=self.obsfreq.value,
                   bandwidth=self.obsbw.value,
                   Nsubband=self.nchan,
                   sample_rate=s_rate,
                   dtype=np.float32,
                   #subint=True,
                   fold = True,
                   sublen=self.tsubint)

        S._dat_freq = make_quant(self._get_pfit_bin_table_entry('SUBINT', 'DAT_FREQ'), 'MHz')
        
        return S

    def copy_psrfit_BinTables(self, ext_names='all', copy_SUBINT_nonDATA=True):
        """Method to copy BinTables from the PSRFITS file given as the template.

        Parameters
        ----------

        ext_names : list, 'all'
            List of BinTable Extensions to copy. Defaults to all, but does not
            copy DATA array in SUBINT BinTable.

        copy_SUBINT_nonDATA : bool, optional
            If True all of the contents of the `SUBINT` BinTable will be copied
            except for the DATA array. This includes the single floats as well
            as TSUBINT, OFFS_SUB, LST_SUB, RA_SUB, DEC_SUB, GLON_SUB, GLAT_SUB,
            FD_ANG, POS_ANG, PAR_ANG, TEL_AZ, TEL_ZEN and others. Is set to
            False if ext_names does not include 'SUBINT' or not equal to 'all'.
            Does not set DAT_FREQ, DAT_SCL, DAT_WTS, DAT_OFFS.
        """
        if ext_names == 'all':
            ext_names = self.file.draft_hdr_keys[1:]
        if 'SUBINT' not in ext_names:
            copy_SUBINT_nonDATA = False

        ext_names.remove('SUBINT')
        for ky in ext_names:
            self.file.copy_template_BinTable(ky)

        if copy_SUBINT_nonDATA:
            self.file.set_subint_dims(nbin=self.nbin, nsblk=self.nsblk,
                                      nchan=self.nchan, nsubint=self.nrows,
                                      npol=self.npol,
                                      data_dtype=self.dtypes['DATA'][0],
                                      obs_mode=self.pfit_dict['OBS_MODE'])

            self.file.copy_template_BinTable(ext_name='SUBINT',
                                            cols=self.file.single_subint_floats)


    def to_txt(self):
        """Convert file to txt file.
        """
        raise NotImplementedError()

    def to_psrfits(self):
        """Convert file to PSRFITS file.
        """
        return self

    def set_sky_info(self):
        """Enter various astronomical and observing data into PSRFITS draft."""
        raise NotImplementedError()

    def _calc_psrfits_dims(self, signal):
        """
        Calculate the needed dimensions of a PSRFITS file from a PSS.signal
        object.
        """
        self.nchan = signal.Nf
        Nt = signal.Nt

    def _get_signal_params(self, signal = None):
        """
        Calculate/retrieve the various parameters to make a PSS signal object
        from a given PSRFITS file.
        
        if signal is given as a Signal() class object, then the values
        will be taken from the signal class instead of a given PSRFITS
        file.
        """
        # get paramters from fits file
        if signal == None:
            self._make_psrfits_pars_dict()
            self.nchan = self.pfit_dict['NCHAN']
            self.tbin = self.pfit_dict['TBIN']
            self.nbin = self.pfit_dict['NBIN']
            self.npol = self.pfit_dict['NPOL']
            self.nrows = self.pfit_dict['NAXIS2']
            self.nsblk = self.pfit_dict['NSBLK']
            self.obsfreq = self.pfit_dict['OBSFREQ']
            self.obsbw = self.pfit_dict['OBSBW']
            self.chan_bw = self.pfit_dict['CHAN_BW']
            self.stt_imjd = self.pfit_dict['STT_IMJD'] # start MJD of obs
            self.stt_smjd = self.pfit_dict['STT_SMJD'] # start second of obs
            self.tsubint = self.pfit_dict['TSUBINT'] # length of subint in seconds
            
            if self.obs_mode=='PSR':
                self.nsubint = self.nrows
            else:
                self.nsubint = None
        # get parameters from signal class
        else:
            self._make_psrfits_pars_dict()
            self.nchan = signal.Nchan
            self.tbin = 1.0/signal.samprate
            self.nbin = int(signal.nsamp/signal.nsub)
            self.npol = self.pfit_dict['NPOL']#signal.Npols
            self.nrows = signal.nsub
            self.nsblk = self.pfit_dict['NSBLK']
            self.obsfreq = signal.fcent
            self.obsbw = signal.bw
            self.chan_bw = signal.bw / signal.Nchan
            #self.stt_imjd = self.pfit_dict['STT_IMJD'] # start MJD of obs
            #self.stt_smjd = self.pfit_dict['STT_SMJD'] # start second of obs
            self.tsubint = signal.sublen # length of subint in seconds
            
            if self.obs_mode=='PSR':
                self.nsubint = self.nrows
            else:
                self.nsubint = None


    def _make_psrfits_pars_dict(self):
        """
        Make parameters from PSRFITS file template into a dictionary for
        initializing PsrSigSim signal object.
        """
        self.pfit_dict = {}

        for extname in self.pfit_pars.keys():
            for ky in self.pfit_pars[extname]:
                if 'DAT' in ky:
                    val = self._get_pfit_bin_table_entry('SUBINT', ky)
                elif 'TSUBINT' in ky:
                    val = self._get_pfit_bin_entry('SUBINT', ky)
                elif extname == "PSRPARAM":
                    val = self._get_pfit_psrparam(extname, ky)
                else:
                    val = self._get_pfit_hdr_entry(extname,ky)
                if isinstance(val,str) or isinstance(val,bytes):
                    val = val.strip()

                self.pfit_dict.update({ky:val})
        idx = self.file.draft_hdr_keys.index('SUBINT')
        dtypes = self.file.get_HDU_dtypes(self.file.fits_template[idx])

        self.dtypes = dict([(item[0],item[1]) if len(item)==2
                            else (item[0],(item[1],item[2]))
                            for item in dtypes])

    def _get_pfit_hdr_entry(self, extname, key):
        """Retrieve a single header entry from PSRFITS file."""
        idx = self.file.draft_hdr_keys.index(extname)
        return self.file.fits_template[idx].read_header()[key]

    def _get_pfit_bin_table_entry(self, extname, key, row=0):
        """Retrieve a single header entry from PSRFITS file."""
        idx = self.file.draft_hdr_keys.index(extname)
        try:
            return self.file.fits_template[idx][key][row][0]
        except:
            return self.file.fits_template[idx][key][row]
    
    def _get_pfit_bin_entry(self, extname, key, row=0):
        """Retrieve a single header entry from PSRFITS file.
        Different from get_pfit_bin_table_entry, this gets just
        single value parameters (e.g. TSUBINT), not arrays or list values.
        """
        idx = self.file.draft_hdr_keys.index(extname)
        return self.file.fits_template[idx][key][row]
    
    def _get_pfit_psrparam(self, extname, param):
        """Retrieve a single value from the PSRPARAM header. This 
        has a different format than the SUBINT or PRIMARY headers.
        """
        idx = self.file.draft_hdr_keys.index(extname)
        for val in self.file.fits_template[idx][:]:
            # Have correction based on version of fitsio being used
            if version.parse(fitsio.__version__) >= version.parse('1.0.1'):
                if param == val[0].split()[0]:
                    return np.float64(val[0].split()[1].replace("D","E"))
            else:
                if param == val[0].split()[0].decode("utf-8"):
                    return np.float64(val[0].split()[1].decode("utf-8").replace("D","E"))

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
    def nsblk(self):
        return self._nsblk

    @nsblk.setter
    def nsblk(self, value):
        self._nsblk = value

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
    def stt_imjd(self):
        return self._stt_imjd
    
    @stt_imjd.setter
    def stt_imjd(self, value):
        self._stt_imjd = make_quant(value, 'day')
        
    @property
    def stt_smjd(self):
        return self._stt_smjd
    
    @stt_smjd.setter
    def stt_smjd(self, value):
        self._stt_smjd = make_quant(value, 'second')
        
    @property
    def tsubint(self):
        return self._tsubint
    
    @tsubint.setter
    def tsubint(self, value):
        self._tsubint = make_quant(value, 'second')

"""ism.py
module to cause dispersion
"""
from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
import numpy as np
import scipy as sp
from scipy import signal
from . import PSS_utils as utils
from . import scintillation as scint

class ISM(object):
    def __init__(self, Signal_in, DM = 30):
        self.Signal_in = Signal_in
        self.signal = self.Signal_in.signal
        self.f0 = self.Signal_in.f0
        self.bw = self.Signal_in.bw
        self.Nf = self.Signal_in.Nf
        self.Nt = self.Signal_in.Nt
        self.TimeBinSize =self.Signal_in.TimeBinSize
        self.freqBinSize = self.Signal_in.freqBinSize
        self.first_freq = self.Signal_in.first_freq
        self.last_freq = self.Signal_in.last_freq
        self.freq_Array = self.Signal_in.freq_Array
        self.DM = DM
        self.ISM_Dict = dict(dispersion=False, scattering=False, DM = None, scintillation=False)

    def shiftit(self, y, shift):
        """
        shifts array y by amount shift (in sample numbers)
        uses shift theorem and FFT
        shift > 0  ==>  lower sample number (earlier)
        modeled after fortran routine shiftit
        Optimized from JMC's code by Michael Lam
        """
        #TODO Add Try Except for odd length arrays...
        yfft = np.fft.fft(y)
        size = np.size(y) #saves time
        constant = (shift*2*np.pi)/float(size) #needs a negative here for the right direction, put it in?
        theta = constant*np.arange(size)
        c = np.cos(theta)
        s = np.sin(theta)
        work = np.zeros(size, dtype='complex')
        work.real = c * yfft.real - s * yfft.imag
        work.imag = c * yfft.imag + s * yfft.real
        # enforce hermiticity
        half_size = int(size//2)
        work.real[half_size:] = work.real[half_size:0:-1]
        work.imag[half_size:] = -work.imag[half_size:0:-1]
        work[half_size] = 0.+0.j
        workifft = np.fft.ifft(work)
        return workifft.real

    def disperse(self, to_DM_Broaden = False):
        #Function to calculate the dispersion per frequency bin for 1/f^2 dispersion
        self.ISM_Dict["DM"] = self.DM
        self.ISM_Dict['DM_Broaden'] = to_DM_Broaden
        self.ISM_Dict['dispersion'] = True
        if self.Signal_in.SignalType=='intensity':
            #For intensity signal calculate dispersion for all sub-bands.
            self.K = 1.0/2.41e-4 #constant used to be more consistent with PSRCHIVE
            self.time_delays = -1e-3*self.K*self.DM*(np.power((self.freq_Array/1e3),-2)) #freq in MHz, delays in milliseconds
                #Dispersion as compared to infinite frequency
            self.time_delays = np.rint(self.time_delays//self.TimeBinSize) #Convert to number of bins
            self.widths = np.zeros(self.Nf)
            for ii, freq in enumerate(self.freq_Array):
                self.signal[ii,:] = self.shiftit(self.signal[ii,:], self.time_delays[ii])
                sub_band_width = self.bw/self.Nf
                width = int(utils.top_hat_width(sub_band_width, freq, DM)//self.TimeBinSize)
                if width > 0 and to_DM_Broaden:
                    if width > self.Nt:
                        raise ValueError('Too Much DM! Dispersion broadening top hat wider than data array!')
                    self.widths[ii] = width
                    self.signal[ii,:] = sp.signal.convolve(self.signal[ii,:], sp.signal.boxcar(width)/width, mode='same',method='fft').astype(self.Signal_in.data_type)
                    # The division by width of the boxcar filter normalizes the convolution

                    #print(self.freq_Array[ii],' MHz ','width=', width) #for debugging
        elif self.Signal_in.SignalType=='voltage':
            #For voltage signal disperse coherently.
            raise ValueError('Sorry, Voltage-type signal dispersion is not currently supported!')
            #for ii in range(4): #Maybe faster to do the complex fft with two channels.
            #    sig_FFT = np.fft.rfft(self.signal[ii,:])
            #    fft_len = len(sig_FFT)
            #    f_array = np.linspace(-(self.last_freq)*1e6,0,length2)
            #    disp_signal_fft = sig_FFT*np.exp(1j*2*np.pi*4.148808e9/((freq+f0)*f0**2)*DM*freq**2)
            #    self.signal[ii,:] = np.fft.irfft(disp_signal_fft)

        self.Signal_in.MetaData.AddInfo(self.ISM_Dict)

    def scatter(self, array, scat_timescale):
        """
        Simulate scatter broadening by convolving the signal with an exp(-t/tau).
        """
        nBins = self.Signal_in.MetaData.nBins_per_period
        tau = scat_timescale / self.TimeBinSize
        try:
            #N_taus = nBins/tau
            exp_time = np.linspace(0,nBins,nBins)
            scatter_exp = np.exp(-exp_time/tau)
            scatter_exp /= np.sum(scatter_exp)
            return sp.signal.convolve(array, scatter_exp, mode='full',method='fft')[:-nBins]
            #.astype(self.Signal_in.data_type)
        except: #Exception if meant for tau too small for given sampling rate.
            return array

    class scintillate:
        def __init__(self, V_ISS = None, scint_timescale = None, pulsar= None, to_use_NG_pulsar=False, telescope=None, freq_band=None):
            """
            Uses a phase screen with the given power spectrum to scintillate a pulsar signal
            across an observation band. The class uses the parameters given to calculate
            thin phase screens and gain image using Fresnel propagation.

            The screens are calculated for the size appropriate to the given parameters
            and observation length.
            """

            if pulsar == None and V_ISS==None and scint_timescale==None:
                raise ValueError('Need to set a variable that sets the scintillation timescale.')

            if pulsar != None and to_use_NG_pulsar:
                if telescope==None or freq_band==None:
                    raise ValueError('Must set both the telescope and bandwidth for {0}.'.format(pulsar))

                self.scint_bw, self.scint_time = NG_scint_param(pulsar, telescope, freq_band)

                if scint_timescale != None:
                    print('Overiding scint_timescale value. Scintillation timescale set to {0} using Lam, et al. 2015.'.format(self.scint_time))
                    print('Change to_use_NG_pulsar flag to use entered value.')
                if V_ISS != None :
                    print('Overiding V_ISS value. Scintillation timescale set to {0} using Lam, et al. 2015.'.format(self.scint_time))
                    print('Change to_use_NG_pulsar flag to use entered value.')

            if pulsar == None and V_ISS==None and scint_timescale!=None:
                self.scint_time = scint_timescale
            if pulsar == None and V_ISS!=None and scint_timescale==None:
                raise ValueError('V_ISS calculation not currently supported.')


            diff_phase_screen = scint.phase_screen(self.Signal_in, DM, Number_r_F=1/64.)

            L = np.rint(diff_phase_screen.xmax//diff_phase_screen.r_Fresnel)

            refrac_phase_screen = scint.phase_screen(self.Signal_in, DM, Number_r_F=5)
            #Calculate a refraction screen to give a correction.

            self.gain = scint.images(diff_phase_screen, self.Signal_in, mode='simulation').gain


        def make_scintles(self, verbose=True):
            """
            Scintillate the pulsar signal.
            """
            #self.gain

        def NG_scint_param(pulsar, telescope, freq_band):
            """ Method for pulling scintillation bandwidth (MHz) and scintillation timescale (sec)
            from a txt file.
            pulsar = Any of the NANOGrav pulsars from 9yr Data release in file.
                        See 'PTA_pulsar_nb_data.txt' for details.
            telescope  = 'AO' (Arecibo Obs) or 'GBT' (Greenbank Telescope)
            freq_band = [327 ,430, 820, 1400, 2300]
            """
            freq_bands_txt = np.array(['0.327','0.430','0.820','1.400','2.300'],dtype=str)
            freq_band = np.extract(freq_band==freq_bands_txt.astype(float)*1e3,freq_bands_txt)[0]

            search_list = (pulsar, telescope, freq_band)
            columns = (10,11)
            try:
                scint_bw, scint_timescale = text_search(search_list, columns, 'PTA_pulsar_nb_data.txt')
            except:
                raise ValueError('Combination of pulsar {0}, telescope {1} and bandwidth {2} MHz'.format(pulsar, telescope, freq_band)+' not found in txt file.')

            return scint_bw, scint_timescale

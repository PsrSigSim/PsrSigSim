from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
import numpy as np
from scipy import stats
import scipy.signal as spsig
import sys, time
from ..utils.utils import make_quant, shift_t
from ..utils.constants import DM_K
from ..utils.constants import KOLMOGOROV_BETA
from ..pulsar.portraits import DataPortrait

class ISM(object):
    '''
    Class for modeling interstellar medium effects on pulsar signals.
    '''
    def __init__(self):
        ''''''
        pass

    def disperse(self, signal, dm):
        r"""
        Function to calculate the dispersion per frequency bin for :math:`1/\nu^{2}`
        dispersion.

        .. math::
            \Delta t_{\rm{DM}} = 4.15\times 10^{6}~\rm{ms} \times \rm{DM} \times \frac{1}{\nu^{2}}
        """
        signal._dm = make_quant(dm,'pc/cm^3')

        if hasattr(signal,'_dispersed'):
            raise ValueError('Signal has already been dispersed!')

        if signal.sigtype=='FilterBankSignal':
            self._disperse_filterbank(signal, signal._dm)
        elif signal.sigtype=='BasebandSignal':
            self._disperse_baseband(signal, signal._dm)

        signal._dispersed = True

    def _disperse_filterbank(self, signal, dm):
        #freq in MHz, delays in milliseconds
        freq_array = signal._dat_freq
        time_delays = (DM_K * dm * np.power(freq_array,-2)).to('ms')
        if signal.delay==None:
            signal._delay=time_delays
        else:
            signal._delay += time_delays
        #Dispersion as compared to infinite frequency
        shift_dt = (1/signal._samprate).to('ms')
        shift_start = time.time()
        # check if there are less than 20 frequency channels
        if signal.Nchan <= 20:
            div_fac = 1
        else:
            div_fac = 20

        for ii, freq in enumerate(freq_array):
            signal._data[ii,:] = shift_t(signal._data[ii,:],
                                         time_delays[ii].value,
                                         dt=shift_dt.value)
            if (ii+1) % int(signal.Nchan//div_fac) ==0:
                shift_check = time.time()
                percent = round((ii + 1)*100/signal.Nchan)
                elapsed = shift_check-shift_start
                chk_str = '\r{0:2.0f}% dispersed'.format(percent)
                chk_str += ' in {0:4.3f} seconds.'.format(elapsed)

                try:
                    print(chk_str , end='', flush=True)
                #This is the Python 2 version
                #__future__ does not have 'flush' kwarg.
                except:
                    print(chk_str , end='')
                sys.stdout.flush()

    def _disperse_baseband(self, signal, dm):
        """
        Broadens & delays baseband signal w transfer function defined in PSR
        Handbook, D. Lorimer and M. Kramer, 2006
        Returns a baseband signal dispersed by the ISM.
        """
        for x in range(signal.Nchan):
            sig = signal._data[x]
            f0 = signal._fcent
            dt = (1/signal._samprate).to('us')

            fourier = np.fft.rfft(sig)
            u = make_quant(np.fft.rfftfreq(2 * len(fourier) - 1,
                                d=dt.to('s').value), 'MHz')
            f = u-signal.bw/2. # u in [0,bw], f in [-bw/2, bw/2]

            # Lorimer & Kramer 2006, eqn. 5.21
            H = np.exp(1j*2*np.pi*DM_K/((f+f0)*f0**2)*dm*f**2)

            product = fourier*H
            Dispersed = np.fft.irfft(product)

            signal._data[x] = Dispersed

    def FD_shift(self, signal, FD_params):
        r"""
        This calculates the delay that will be added due to an arbitrary number
        of input FD parameters following the NANOGrav standard as defined in
        Arzoumanian et al. 2016. It will then shift the pulse profiles by the
        appropriate amount based on these parameters.

        FD values should be input in units of seconds, frequency array in MHz
        FD values can be a list or an array

        .. math::
            \Delta t_{\rm{FD}} = \sum_{i=1}^{n} c_{i} \log\left({\frac{\nu}{1~\rm{GHz}}}\right)^{i}.
        """
        #freq in MHz, delays in milliseconds
        freq_array = signal._dat_freq
        # define the reference frequency
        ref_freq = make_quant(1000.0, 'MHz')
        # calculate the delay added in for the parameters
        time_delays = make_quant(np.zeros(len(freq_array)), 'ms') # will be in seconds
        for ii in range(len(FD_params)):
            time_delays += np.double(make_quant(FD_params[ii], 's').to('ms') * \
                    np.power(np.log(freq_array/ref_freq),ii+1)) # will be in seconds

        if signal.delay==None:
            signal._delay=time_delays
        else:
            signal._delay += time_delays
        # get time shift based on the sample rate
        shift_dt = (1/signal._samprate).to('ms')
        shift_start = time.time()
        # check if there are less than 20 frequency channels
        if signal.Nchan <= 20:
            div_fac = 1
        else:
            div_fac = 20

        for ii, freq in enumerate(freq_array):
            signal._data[ii,:] = shift_t(signal._data[ii,:],
                                         time_delays[ii].value,
                                         dt=shift_dt.value)
            if (ii+1) % int(signal.Nchan//div_fac) ==0:
                shift_check = time.time()
                percent = round((ii + 1)*100/signal.Nchan)
                elapsed = shift_check-shift_start
                chk_str = '\r{0:2.0f}% shifted'.format(percent)
                chk_str += ' in {0:4.3f} seconds.'.format(elapsed)

                try:
                    print(chk_str , end='', flush=True)
                #This is the Python 2 version
                #__future__ does not have 'flush' kwarg.
                except:
                    print(chk_str , end='')
                sys.stdout.flush()

        # May need to add tihs parameter to signal
        signal._FDshifted = True

    def scatter_broaden(self, signal, tau_d, ref_freq, beta = KOLMOGOROV_BETA, \
                        convolve = False, pulsar = None):
        """
        Function to add scatter broadening delays to simulated data. We offer
        two methods to do this, one where the delay is calcuated and the
        pulse signals is directly shifted by the calculated delay (as done
        in the disperse function), or the scattering delay exponentials are
        directy convolved with the pulse profiles. If this option is chosen,
        the scatter broadening must be done BEFORE pulsar.make_pulses() is run.

        Parameters
        ----------

        signal [object] : signal class object which has been previously defined
        tau_d [float] : scattering delay [seconds]
        ref_freq [float] : reference frequency [MHz] at which tau_d was measured
        beta [float] : preferred scaling law for tau_d, default is for a
                       Kolmoogorov medium (11/3)
        convolve [bool] : If False, signal will be directly shifted in time by
                          scattering delay; if True, scattering delay tails
                          will be directly convolved with the pulse profiles.
        pulsar [object] : previously defined pulsar class object with profile
                          already assigned
        """
        # First get and define values to use
        freq_array = signal._dat_freq
        ref_freq = make_quant(ref_freq, 'MHz')
        tau_d = make_quant(tau_d, 's').to('ms') # need compatible units with dt
        # Scale the scattering timescale with frequency
        tau_d_scaled = self.scale_tau_d(tau_d, ref_freq , freq_array, beta=beta)
        # First shift signal if convolve = False
        if not convolve:
            if signal.delay==None:
                signal._delay=tau_d_scaled
            else:
                signal._delay += tau_d_scaled
            # define bin size to shift by
            shift_dt = (1/signal._samprate).to('ms')
            shift_start = time.time()
            # check if there are less than 20 frequency channels
            if signal.Nchan <= 20:
                div_fac = 1
            else:
                div_fac = 20
            # now loop through and scale things appropriately
            for ii, freq in enumerate(freq_array):
                signal._data[ii,:] = shift_t(signal._data[ii,:],
                                             tau_d_scaled[ii].value,
                                             dt=shift_dt.value)
                if (ii+1) % int(signal.Nchan//div_fac) ==0:
                    shift_check = time.time()
                    percent = round((ii + 1)*100/signal.Nchan)
                    elapsed = shift_check-shift_start
                    chk_str = '\r{0:2.0f}% scatter shifted'.format(percent)
                    chk_str += ' in {0:4.3f} seconds.'.format(elapsed)

                    try:
                        print(chk_str , end='', flush=True)
                    #This is the Python 2 version
                    #__future__ does not have 'flush' kwarg.
                    except:
                        print(chk_str , end='')
                    sys.stdout.flush()
        else:
            # Make the initial profile data array at correct sample rate
            Nph = int((signal.samprate * pulsar.period).decompose())
            pulsar.Profiles.init_profiles(Nph, signal.Nchan)
            # Now make the profiles with Nph bins
            phs = np.linspace(0.0, 1.0, Nph)
            # full_profs is a data array of the profiles
            full_profs = pulsar.Profiles.calc_profiles(phs, signal.Nchan)
            # Now we make the array of exponential scattering tails
            t = np.linspace(0, pulsar.period, Nph)
            EXP_array = np.zeros((signal.Nchan, Nph))
            #Iterating over the tau arrays where each profile corresponds to the respective tau index
            for ii, tau_scatter in enumerate(tau_d_scaled):
                EXP = (np.exp(-t/tau_scatter))
                EXP_array[ii,:] = EXP
            # Now we convolve the profiles
            convolved_profs = self.convolve_profile(full_profs, EXP_array, \
                                                  width=Nph)
            # Now we reassign the pulsar profile object
            pulsar._Profiles = DataPortrait(convolved_profs)


    def convolve_profile(self, profiles, convolve_array, width = 2048):
        """
        Function to convolve some array generated by a function with the
        previously assigned pulsar pulse profiles. Main use case is in
        convolving exponential scattering tails with the input pulse profiles,
        however any input array can be convolved.

        NOTE: This function only returns the array of convolved profiles, it
        does NOT reassign the pulsar objects profiles.

        Parameters
        ----------

        profiles [array] : data array of pulse profiles generated with the
                           'calc_profiles' function.
        convolve_array [array] : data array representing the function that
                                 will be convolved with the pulse profiles.
                                 Should be the same shape as 'profiles'.
        width [int] : number of bins desired from the resulting convolved
                      profiles. Default is 2048 bins across the profile.
                      Should be the same number as the original input profiles.
        """
        for ii, freq in enumerate(convolve_array):
            #Normalizing the pulse profile
            pulsar_prof_sum = np.sum(profiles[ii,:])
            # check divide by zero
            if pulsar_prof_sum != 0.0:
                pulsar_prof_norm = profiles[ii,:] / pulsar_prof_sum
            else:
                pulsar_prof_norm = profiles[ii,:]

            #Normalizing the input array
            convolve_array_sum = np.sum(convolve_array[ii,:])
            # Check divide by zero
            if convolve_array_sum != 0.0:
                convolve_array_norm = convolve_array[ii,:] / convolve_array_sum
            else:
                convolve_array_norm = convolve_array[ii,:]

            #Convolving the input array with the pulse profile
            convolved_prof = spsig.convolve(pulsar_prof_norm, convolve_array_norm, \
                                            mode='full',method='fft')

            #Renormalizing the convolved pulse profile
            profiles[ii,:] = (pulsar_prof_sum)*(convolved_prof[:width])
        return profiles

    '''
    Written by Michael Lam, 2017
    Scale dnu_d and dt_d based on:
    dnu_d propto nu^(22/5)
    dt_d propto nu^(6/5) / transverse velocity
    See Stinebring and Condon 1990 for scalings with beta (they call it alpha)

    TODO: Should units be assigned here, or earlier?
    '''

    def scale_dnu_d(self,dnu_d,nu_i,nu_f,beta=KOLMOGOROV_BETA):
        """
        Scaling law for scintillation bandwidth as a function of frequency.

        Parameters
        ----------

        dnu_d [float] : scintillation bandwidth [MHz]
        nu_i [float] : reference frequency at which du_d was measured [MHz]
        nu_f [float] : frequency (or frequency array) to scale dnu_d to [MHz]
        beta [float] : preferred scaling law for dnu_d, default is for a
                       Kolmogorov medium (11/3)
        """
        #dnu_d = make_quant(dnu_d, 'MHz')
        if beta < 4:
            exp = 2.0*beta/(beta-2) #(22.0/5)
        elif beta > 4:
            exp = 8.0/(6-beta)
        return dnu_d*(nu_f/nu_i)**exp

    def scale_dt_d(self,dt_d,nu_i,nu_f,beta=KOLMOGOROV_BETA):
        """
        Scaling law for scintillation timescale as a function of frequency.

        Parameters
        ----------

        dt_d [float] : scintillation timescale [seconds]
        nu_i [float] : reference frequency at which du_d was measured [MHz]
        nu_f [float] : frequency (or frequency array) to scale dnu_d to [MHz]
        beta [float] : preferred scaling law for dt_d, default is for a
                       Kolmoogorov medium (11/3)
        """
       # dt_d = make_quant(dt_d, 's')
        if beta < 4:
            exp = 2.0/(beta-2) #(6.0/5)
        elif beta > 4:
            exp = float(beta-2)/(6-beta)
        return dt_d*(nu_f/nu_i)**exp

    def scale_tau_d(self,tau_d,nu_i,nu_f,beta=KOLMOGOROV_BETA):
        """
        Scaling law for the scattering timescale as a function of frequency.

        Parameters
        ----------

        tau_d [float] : scattering timescale [seconds?]
        nu_i [float] : reference frequency at which du_d was measured [MHz]
        nu_f [float] : frequency (or frequency array) to scale dnu_d to [MHz]
        beta [float] : preferred scaling law for tau_d, default is for a
                       Kolmoogorov medium (11/3)
        """
        #tau_d = make_quant(tau_d, 's')
        if beta < 4:
            exp = -2.0*beta/(beta-2) #(-22.0/5)
        elif beta > 4:
            exp = -8.0/(6-beta)
        return tau_d*(nu_f/nu_i)**exp

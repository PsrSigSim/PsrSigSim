from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
import numpy as np
from scipy import stats
import time
from ..utils.utils import make_quant, shift_t
from ..utils.constants import DM_K

class ISM(object):
    '''
    Class for modeling interstellar medium effects on pulsar signals.
    '''
    def __init__(self):
        ''''''
        pass

    def disperse(self, signal, dm):
        """
        Function to calculate the dispersion
        per frequency bin for 1/f^2 dispersion
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
        #Dispersion as compared to infinite frequency
        shift_dt = (1/signal._samprate).to('ms')
        shift_start = time.time()

        for ii, freq in enumerate(freq_array):
            signal._data[ii,:] = shift_t(signal._data[ii,:],
                                         time_delays[ii].value,
                                         dt=shift_dt.value)
            if (ii+1) % int(signal.Nchan//20) ==0:
                shift_check = time.time()
                percent = round((ii + 1)*100/self.Nf)
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
        Use plot_dispersed() in PSS_plot.py to see the
        dispersed and undispersed signal.
        """
        for x in range(signal.Nchan):
            sig = signal._data[x]
            f0 = signal._fcent
            dt = (1/signal._samprate).to('us')

            fourier = np.fft.rfft(sig)
            u = np.fft.rfftfreq(2 * len(fourier) - 1,
                                d=dt.to('s').value) * u.us
            f = u-signal.bw/2. # u in [0,bw], f in [-bw/2, bw/2]

            # Lorimer & Kramer 2006, eqn. 5.21
            H = np.exp(1j*2*np.pi*DM_K/((f+f0)*f0**2)*dm*f**2)

            product = fourier*H
            Dispersed = np.fft.irfft(product)

            if self.MD.mode == 'explore':
                self.Signal_in.undispersedsig[x] = sig
            signal._data[x] = Dispersed

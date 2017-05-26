from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
import numpy as np
import math as m
import pylab as pylab
import scipy as sp
from scipy import stats
from scipy import constants
from scipy.optimize import curve_fit
from scipy import misc
import scipy.ndimage as  ndimage
from . import PSS_utils as utils

KOLMOGOROV_BETA = 11.0/3
c = 2.998e8 #m/s
#lam = c/nu
KPC_TO_CM = 3.086e21
KPC_TO_M = 3.086e19
RAD_TO_MAS = 2.063e8
r_e = 2.8179403227e-15 #m

class phase_screen:
    def __init__(self, signal_object, DM, Nx=30, Ny=30, spectral_ind=11./3., D_pulsar=1, D_screen=0.5, inner=None, outer=None, rmsfres=None, scint_bandwidth=180e6, apply_inner=False,  apply_outer=False, normfres=False, normDM=False):
        """
        Generates npoints of a realization of power-law noise with unit
        variance with spectral index, spectral_ind and inner and outer scales as
        specified for a sample interval dx.

        C_n^2, the power-law amplitude of electron-density wavenumber spectrum, is calculated
        seperately by the cnsq_calc() method.

        input:
    	spectral_ind = spectral index of power-law wavenumber spectrum
    	rmsfres = rms phase at Fresnel scale (rad)
            length scales: all dimensionless:
    		freq	= frequency in MHz
    		inner, outer	= inner and outer scales
    		dx, dy		= sample intervals
    		xwidth, ywidth	= screen extent
            D = distance to screen in kPc
            D_p = Distance to Pulsar in kPc
            scint_bandwidth = scintillation bandwidth in MHz
    	logical:
            apply_inner
            apply_outer
            normfres 	= True implies normalization to rmsfres

            Definition of Fresnel scale: r_F^2 = \lambda D / 2\pi

        returns:
           xvec, yvec, xseries, xseries_norm, qxvec, qyvec, qshape
           (xvec, yvec) coordinates perpendicular to line of sight
           xseries = screen phase
           xseries_norm = screen phase scaled to input rms phase on Fresnel scale
           qxvec, qyvec  = wavenumber axes
           qshape = sqrt(shape of wavenumber spectrum)

        references: CR98 is Cordes and Rickett, 1998
                    LK12 is Lorimer and Kramer, 2012
        """

        #D_screen *= KPC_TO_M
        max_freq = signal_object.freq_Array.max()
        min_freq = signal_object.freq_Array.min() #signal_object.f0
        self.r_Fresnel = np.sqrt(c / (max_freq * 1e6) * (D_screen * KPC_TO_M) / (2.0*np.pi))

        #Below we use ray optics to find the Field Coherence Scale, though other approaches are possible
        #Use the uncertainty relation from CR98 (C_1=1.16) to calculate the length of the hypotenuse for scattering

        #delta_x = 1.16 * c / (2 * np.pi * scint_bandwidth)
        #or
        scat_timescale = 10**(-6.46+0.154*np.log10(DM)+1.07*(np.log10(DM))**2-3.86*np.log10(signal_object.f0/1e3))
        self.Freq_diss = 1.16/(2*np.pi*scat_timescale*1e-3)/1e6
        #TODO Could Use the scat_timescale calculated at each frequency...


        #delta_x = c * scat_timescale * 1e-3 #1/e^3 gives > 95% of rays
        #s_0 = np.sqrt(2 * delta_x * D_screen + delta_x**2) #Field Coherence Scale in meters
        #factor = s_0 / self.r_Fresnel # multiplicative factor since the uncertainty is based on 1/e time.
        Number_r_F = 1./64
        xwidth = Number_r_F * self.r_Fresnel * Nx
        ywidth = Number_r_F * self.r_Fresnel * Ny
        #xwidth = Nx * self.dx
        #ywidth = Ny * self.dy
        #print('scat_timescale', scat_timescale)
        print('Freq_diss ', self.Freq_diss)
        #print('delta_x', delta_x)
        #print('D_screen', D_screen)
        #print('Field Coherence Scale', s_0)
        #print('Fresnel Length:', self.r_Fresnel)
        #print('xwidth = ', xwidth)
        #print('s0/r_F', factor)
        xmax = xwidth/2
        ymax = ywidth/2

        #lscale_index = 2. / (si_kol_d-2.)
        #u = max(1., rmsfres_d**lscale_index)
        #ld = rfres / u
        #lr = rfres * u
        #inners = 2.*ld

        #Nx = 30#int(xwidth//dx) #Number of x bins
        #Ny = 30#int(ywidth//dy) #Number of y bins

        self.dx = xwidth / Nx
        self.dy = ywidth / Ny
        #########################################

        self.xvec = (np.arange(0.,Nx)- Nx//2 + 1)*self.dx # Make arrays that are centered on zero, over the given widths
        self.yvec = (np.arange(0.,Ny)- Ny//2 + 1)*self.dy

        dqx = 2.*np.pi / xwidth #Smallest wavenumber in the given direction
        dqy = 2.*np.pi / ywidth
        qmaxx = (2.* np.pi) / (2.* self.dx) #Largest wavenumber in the given direction
        qmaxy = (2.* np.pi) / (2.* self.dy)
        #print('dqx',dqx)
        #print('qmaxx',qmaxx)
        Nqx = 2*int(qmaxx//dqx) # Number of wavenumber samples
        Nqy = 2*int(qmaxy//dqy)
        #print('targeted number of q samples = ', Nqx, Nqy )
        if Nqx != Nx:
        #    print("Forcing Nqx = Nx = ", Nx)
            Nqx = Nx
        if Nqy != Ny:
        #    print("Forcing Nqy = Ny = ", Ny)
            Nqy = Ny

        self.qxvec_centered = (np.arange(0.,Nqx)-Nqx//2+1)*dqx # Make arrays that are centered on zero, over the given widths
        self.qxvec = np.roll(self.qxvec_centered,Nqx//2+1) #Shift array by enough to make 0 first
        self.qyvec_centered = (np.arange(0.,Nqy)-Nqy//2+1)*dqy
        self.qyvec = np.roll(self.qyvec_centered, Nqy//2+1)

        #Set inner and outer scale factors if not otherwise set.
        if inner == None:
            inner = xwidth/(Nx)#self.r_Fresnel/2.

        if outer == None:
            outer = xwidth*(1.0 + 0.2)# + 100*Number_r_F * self.r_Fresnel#1.5*self.r_Fresnel#xwidth

        qin = 2. * np.pi / inner #Set inner and outer scale wavenumber
        qout = 2. * np.pi / outer

        if apply_outer:
            print("Applying outer-scale rolloff")

        # 2016 Jan 1: put in upper wavenumber cutoff at array size
        # to avoid aliasing
        qmax = self.qxvec.max()/2.

        qxy = np.meshgrid(self.qxvec, self.qyvec, indexing='ij')

        qsq = qxy[0]**2 + qxy[1]**2
        self.qshape = (qout**2 + qsq)**(-spectral_ind/4.) * np.exp(-qsq/(2.*qmax**2))
        self.qshape_rolloff = np.exp(-qsq / (2.*qin**2))

        qxy_centered = np.meshgrid(self.qxvec_centered, self.qyvec_centered, indexing='ij')
        self.qsq_centered =  qxy_centered[0]**2 + qxy_centered[1]**2

        if apply_outer:
            self.qshape *= np.exp(-qout**2 / (2.*qsq))

        #self.qshape *= cnsq_calc(dnud=scint_bandwidth, nu=freq) #Normalization constant

        #if high_freq==None:
        #    phase_norm = 1
        #else:
        #    phase_norm = high_freq / freq

        ## new 2016 Jan 1:  create real white noise in x domain and FFT
        ## to get Hermitian noise
        rand_pull_r = np.random.randn(Nqx, Nqy)
        rand_pull_i = np.random.randn(Nqx, Nqy)
        #
        #self.spectrum = np.zeros((signal_object.Nf, Nqx, Nqy))
        #Cn_squared = np.array([cnsq_calc(taud=scat_timescale*1e-3, nu=jj) for jj in signal_object.array])

        def CnSq_calc(f,Freq_DISS,d=0.5):
            return 0.002*(f/1e3)**(3.67)*(d)**(-1.83)*(Freq_DISS)**(-0.83)

        if normDM:
            xformr=cnsq_calc(taud=scat_timescale*1e-3, nu=signal_object.f0)*rand_pull_r*self.qshape
            xformi=cnsq_calc(taud=scat_timescale*1e-3, nu=signal_object.f0)*rand_pull_i*self.qshape
            xform = xformr + 1j*xformi
            self.spectrum = abs(xform)**2 #Is xform*xform.conj() faster here?
            #Is the spectrum what we should be multiplying by C^2?
            self.phi = np.real(np.fft.ifft2(xform))

        else:
            self.phi = np.zeros((signal_object.Nf, Nqx, Nqy)) #   _high
            #self.phi_fft = np.zeros((signal_object.Nf, Nqx, Nqy))
            self.Cn_squared = []
            for ii, freq in enumerate(signal_object.freq_Array):
                #Using Michael Lam's version
                #scat_timescale = 10**(-6.46 + 0.154*np.log10(DM) + 1.07*(np.log10(DM))**2 - 3.86*np.log10(freq/1e3))
                #xformr = np.sqrt(cnsq_calc(taud=scat_timescale*1e-3, nu=freq))*rand_pull_r*self.qshape
                #xformi = np.sqrt(cnsq_calc(taud=scat_timescale*1e-3, nu=freq))*rand_pull_i*self.qshape
                #self.Cn_squared = np.append(self.Cn_squared,cnsq_calc(taud=scat_timescale*1e-3, nu=freq))

                #Using LK12 Version
                scat_timescale = 10**(-6.46 + 0.154*np.log10(DM) + 1.07*(np.log10(DM))**2 - 3.86*np.log10(freq/1e3))
                Freq_diss = 1.16 / (2*np.pi*scat_timescale*1e-3) / 1e6 #1.16
                xformr = np.sqrt(CnSq_calc(freq, Freq_diss,d=D_screen)) * rand_pull_r * self.qshape
                xformi = np.sqrt(CnSq_calc(freq, Freq_diss,d=D_screen)) * rand_pull_i * self.qshape
                self.Cn_squared = np.append(self.Cn_squared, CnSq_calc(freq,Freq_diss))

                #Using the one freq method
                #xformr = cnsq_calc(taud=scat_timescale*1e-3, nu=signal_object.f0)*rand_pull_r*self.qshape
                #xformi = cnsq_calc(taud=scat_timescale*1e-3, nu=signal_object.f0)*rand_pull_i*self.qshape


                xform = xformr + 1j*xformi
                #self.phi_fft[ii,:,:] = xform
                self.spectrum = abs(xform)**2 #Is xform*xform.conj() faster here?
                #Is the spectrum what we should be multiplying by C^2?
                wave_num_to_phi = np.sqrt((2*np.pi)**3 * (freq*1e6/c)**2 * 0.0198339 * (D_screen * KPC_TO_M))#0.0330054 * D_screen)#
                self.phi[ii,:,:] = np.real(np.fft.ifft2(xform)) * wave_num_to_phi #
            #self.phi *= (Nx*Ny)/(self.dx*self.dy)#2/np.std(self.phi)#1#Nx*Ny/(self.dx*self.dy) #1.34439663e-13 *2/np.std(self.phi)
            #self.phi_fft *= 1.34439663e-13 #What are we converting to here? From wavenumber to phase.
            #self.phi /= Number_r_F**2  # Nx*Ny
            #self.phi /= np.std(self.phi)
            #self.phi /= self.qshape.size
            self.phi /= (self.dx*self.dy)*(Nx*Ny)

        print('1/(self.dx*self.dy)',1/(self.dx*self.dy))
        print('(Nx*Ny)',(Nx*Ny))
        self.Nx = Nx
        self.Ny = Ny
        self.xmax = xmax
        self.ymax = ymax
        self.Nf = signal_object.Nf
        self.freq_Array = signal_object.freq_Array

        # Normalization factor needs to be calculated on pure power-law spectrum
        # before any rolloff at the refraction scale
        if normfres: #This normalization gives the desired DM_rms desired at the end.

            #print('frindx ',frindx)
            self.phi_norm = np.zeros(self.phi.shape)
            if len(self.phi.shape)==3:
                for ii, freq in enumerate(signal_object.freq_Array):
                    frindx = int(np.sqrt(r_Fres_SQ(freq,D=0.5))//self.dx)
                    #print(frindx)
                    x1dcut = self.phi[ii,0,:]
                    var_fres_in = np.var(x1dcut[0:np.size(x1dcut)-frindx]-x1dcut[frindx:])
                    norm_factor = rmsfres / np.sqrt(var_fres_in)
                    self.phi_norm[ii,:,:] = self.phi[ii,:,:] * norm_factor #* max_freq /freq
                    xn1dcut = self.phi_norm[ii,0,:]
                    var_fres_out = np.var(xn1dcut[0:np.size(xn1dcut)-frindx]-xn1dcut[frindx:])
            else:
                frindx = int(self.r_Fresnel//self.dx)
                x1dcut = self.phi[0,:]
                var_fres_in = np.var(x1dcut[0:np.size(x1dcut)-frindx]-x1dcut[frindx:])
                norm_factor = rmsfres / np.sqrt(var_fres_in)
                self.phi_norm = self.phi * norm_factor
                xn1dcut = self.phi_norm[0,:]
                var_fres_out = np.var(xn1dcut[0:np.size(xn1dcut)-frindx]-xn1dcut[frindx:])

            print("index of fresnel scale = ", frindx)
            print(var_fres_in, var_fres_out)



        # applying inner scale now an option with apply_inner = True
        # needs to be applied *after* normalization!
        # now need to recalculate the realization and apply norm_factor
        if apply_inner:
            print("Applying inner-scale rolloff")

           # Recalculate
            xform *= self.qshape_rolloff
            self.spectrum = abs(xform)**2
            self.phi = np.real(np.fft.ifft2(xform))
            #self.phi_norm = self.phi * norm_factor

        if normDM:
        #  Change phase for different frequencies.
            freq_norm = signal_object.f0/signal_object.freq_Array[:,np.newaxis,np.newaxis]
            self.phi = self.phi*freq_norm * 1.34439663e-13

        #if normfres:
        #    self.phi_norm = self.phi_norm*freq_norm


class images(object):
    def __init__(self, phase_screen, signal_object, fourier_mode=True, factor=1, freq=None, high_freq=None, normed=False):

        if high_freq==None:
            phase_norm = 1
        else:
            phase_norm = high_freq #/ freq


        if factor==1:
            xy = np.meshgrid(phase_screen.xvec, phase_screen.yvec, indexing='ij')
            rsqvec = xy[0]**2 + xy[1]**2
        else:
            dx = phase_screen.xmax*2/(factor*phase_screen.Nx)
            dy = phase_screen.ymax*2/(factor*phase_screen.Ny)
            Xvec = (np.arange(0.,factor*phase_screen.Nx) - factor*phase_screen.Nx//2 + 1)*dx
            Yvec = (np.arange(0.,factor*phase_screen.Ny) - factor*phase_screen.Ny//2 + 1)*dy
            xy = np.meshgrid(Xvec, Yvec, indexing='ij')
            rsqvec = xy[0]**2 + xy[1]**2

        self.intensity = np.zeros((phase_screen.Nf,phase_screen.Nx,phase_screen.Ny))
        #field = np.zeros((phase_screen.Nf,phase_screen.Nx,phase_screen.Ny))
        self.kernel = np.zeros(self.intensity.shape, dtype='complex64')
        self.kernelFFT = np.zeros(self.intensity.shape, dtype='complex64')
        self.screenFFT = np.zeros(self.intensity.shape, dtype='complex64')
        self.kernelFFT_check = np.zeros(phase_screen.Nf)
        if normed:
            phase_screen.phi = phase_screen.phi_norm
        for ii, frequ in enumerate(phase_screen.freq_Array):
            r_Fres_squared = r_Fres_SQ(frequ)
            SD = phase_screen.xmax*phase_screen.ymax#np.exp(-rsqvec / SD)*
            if fourier_mode:
                kmax = phase_screen.qsq_centered.max()
                FresnelKernel_Norm = 1 #np.sqrt(r_Fres_squared*np.pi)/2
                kernel0fft = FresnelKernel_Norm * (1+1j) \
                                * np.exp(-1j*phase_screen.qsq_centered/2 \
                                *r_Fres_squared)
                #self.kernelFFT[ii,:,:] = kernel0fft
                screen0 = np.exp(phase_norm * 1j * phase_screen.phi[ii,:,:])
                screen0fft = np.fft.fft2(screen0)
                #self.screenFFT[ii,:,:] = screen0fft
                field0fft = kernel0fft * screen0fft
            else:
                kernel0 = np.exp(1j * rsqvec / (2. * r_Fres_squared)) # See Eq (2.1), (Narayan, 1992)
                screen0 = np.exp(phase_norm * 1j * phase_screen.phi[ii])
                kernel0fft = np.fft.fft2(kernel0)
                screen0fft= np.fft.fft2(screen0)
                field0fft = kernel0fft * screen0fft
                #self.kernel[ii,:,:] = kernel0

            norm = 1#((phase_screen.dx*phase_screen.dy)/( r_Fres_squared))/(phase_screen.Nx*phase_screen.Ny)#2.*np.pi*
            field = norm * np.fft.fftshift(np.fft.ifft2(field0fft))
            self.intensity[ii,:,:] = abs(field**2)
    ### Plots

#    def dynamic_spectrum(self, **kwargs):
#        return PSS_plot.dynamic_spectrum(image_screen, **kwargs)



def r_Fres_SQ(freq, wavelength=None, D=0.5, units=['MHz','kPc']):
    """Returns Fresnel length scale. In desired units"""
    if units != ['MHz','kPc']:
        raise ValueError(units, ': units not currently supported!')

    MHz = 1e6
    D *= KPC_TO_M
    if wavelength == None:
        wavelength = c/(freq*MHz)

    return wavelength * D / (2*np.pi)

def cnsq_calc(nu=1000, dnud=None, taud=None, dtd=None, D=1, PM=None, beta=11./3.,mode='screen',ds=0.001):
    '''
    Function for calculating the electron-density wavenumber spectrum coefficient.
    Following Cordes and Rickett 1998 and subsequent references

    Inputs:
    nu = frequency in MHz
    dtd = \Delta t_d, the scintillation timescale in seconds
    dnud = \Delta \nu_d, the scintillation bandwidth in Hz
    taud = \tau_d, the scattering timescale in seconds
    D = distance in kpc
    beta = power-law index of electron-density wavenumber spectrum
    mode = screen, uniform

    Outputs:
    cnsq = C_n^2, the power-law amplitude of electron-density wavenumber spectrum
    '''

    if mode == 'screen' and beta == KOLMOGOROV_BETA:
        C1 = 0.957
    elif mode == 'screen' and beta == 4:
        C1 = 1.0
    elif mode == 'uniform' and beta == KOLMOGOROV_BETA:
        C1 = 1.16
    elif mode == 'uniform' and beta == 4:
        C1 = 1.53
    else:
        C1 = 1.0 ### See Lambert and Rickett 1998 for calculation

    if taud is None and dnud is None:
        raise ValueError("Require either taud or dnud to be set.")
    elif taud is None:
        taud = C1 / (2*np.pi*dnud)
    elif dnud is None:
        dnud = C1 / (2*np.pi*taud)

    if PM is None:
        vp = 100*1e3 #m/s
    else:
        vp = 4.74*1e3*PM #m/s


    if mode == 'screen': ### For now assume screen distance is half way
        if taud is not None:
            eta0ds = 2*c*taud / (0.5*D*KPC_TO_M*(1-0.5))
        elif dtd is not None:
            pass
        eta0 = eta0ds / (ds*KPC_TO_M)
    elif mode == 'uniform':
        eta0 = 0

    nu *= 1e6
    eta0 = eta0 * KPC_TO_M * RAD_TO_MAS**2
    #print(eta0)
    #TODO Add an inner scale parameter to this method
    l1 = 100000 #m
    #cnsq = (eta0 / 3.2) * (nu/1e9)**4 * (l1/100000)**(1.0/3) * 10**-3.5
    cnsq = 0.002*(nu/1e9)**4 * (D)**(-1.83) * (dnud/1e6)**(-0.83)
    #print("C_n_squared:",cnsq, np.log10(cnsq))

    return cnsq

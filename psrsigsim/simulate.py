"""simulate.py
simulation wrapper
"""
import psrsigsim as PSS
import psrsigsim as PSS
from psrsigsim import Telescope
import sys 
import pint 
from pint.models import polycos
from pint import models
from pint import utils
import os.path

d_path = os.path.dirname(__file__)

default_path = os.path.join(d_path, '../data/')
class Simulation(object):
    
    def __init__(self, psr = None, sim_telescope = None, sim_ism = False, sim_scint = False, \
                 sim_dict = None, sim_file_path = default_path ): 
        """ This class simulates an observing run given a pulsar name or a dictionary of paramters. Once the simulation is called,
    the different classes (signal, pulsar,ism,scint,telescope) can be initialized with the initiatlization functions. init_signal and init_pulsar are the 
    only required functions to run the simulation. If a pulsar name is input, then these function are ran automatically. The required paramaters needed
    for the simulation are:
    
        @param f0 -- central frequency (MHz)
        @param F0 -- Pulsar Spin Frequency (Hz)
        @param bw -- bandwidth (MHz)
        @param Nf -- number of freq. bins
        @param ObsTime --  total time of observation (add units once confirmed)
        @param f_samp -- sampling frequency (MHz) 
        @param SignalType -- 'intensity' which carries a Nf x Nt filterbank of pulses or 'voltage' which carries a 4 x Nt array of voltage vs. time pulses representing 4 stokes channels
        @param dm -- dispersion measure
        @param V_ISS -- Intersteller Scintilation Velocity
        @param scint_bw -- scintilation Bandwidth
        @param scint_timescale -- scintilation timescale
        @param pulsar -- pulsar name
        @param telescope -- telescope name(GBT or Arecibo)
        @param freq_band -- frequency band [327 ,430, 820, 1400, 2300]
        @param aperature -- aperature (m)
        @param area -- collecting area (m^2)
        @param Tsys -- system temp (K), total of receiver, sky, spillover, etc. (only needed for noise)
        @param name -- GBT or Arecibo
        @param tau_scatter -- scattering time (ms)
        @param radiometer_noise -- changes whther radiometer noise is included or not. True/False
        @param to_DM_broaden -- changes whether dm broadening is included. True/False
        
        Parameters retrieved when pulsar name is input : F0, dm, scint_bw, scint_timescale. 
        
        Example Simulations:
         Using Pulsar name: 
             $ dictionary = {'f0':1400,'bw':400,'Nf':100,'f_samp':4,ObsTime':10,'data_type':'int8','SignalType':"intensity",'flux':3 \
             ,'tau_scatter':6e-05,'radiometer_noise: True, to_DM_broaden: True}
             $ s1 = Simuation(psr = 'J1614-2230',sim_telescope ='GBT', sim_ism = True, sim_scint = True, sim_dict = None)
             $ s1.simulate()
         Using own dictionary:
             $ dictionary = {'f0':1400,'dm':15.921200, 'F0': 218, 'bw':400,'Nf':100,'f_samp':4,ObsTime':10,'data_type':'int8'\
             ,'SignalType':"intensity",'flux':3,'tau_scatter':6e-05,'radiometer_noise: True, to_DM_broaden: True}
             $ s2 = Simulation(psr = None,sim_telescope = 'Arecibo',sim_ism = None,sim_scint = None,sim_dict = dictionary)
             $ s2.init_signal()
             $ s2.init_pulsar()
             $ s2.init_ism()
             $ s2.init_telescope()
             $ s2.simulate()"""
        self.psr = psr
        self.sim_telescope = sim_telescope
        self.sim_ism = sim_ism
        self.sim_scint = sim_scint
        self.sim_dict = sim_dict
        self.sim_file_path = sim_file_path
        
        

        if isinstance(psr, str):
            self.param_dict = self._set_pulsar_dict()
            print('Five warnings should appear. Warnings that are normal are: Unrecognized parfile line "SOLARN0", "T2CMETHOD", "MODE", "INFO -f". Warning about "DDK model" is also normal.  ')
            if 'F0' in self.param_dict.keys() and 'F0' in self.sim_dict.keys():
                print( "Two 'F0' values input. Simulator will default to input dictionary value.")
            if 'dm' in self.param_dict.keys() and 'dm' in self.sim_dict.keys():
                print( "Two 'dm' values input. Simulator will default to input dictionary value.")    

            self.param_dict.update(self.sim_dict)
            self.sim_dict = self.param_dict
            
            
            self.init_signal()
            self.init_pulsar()
            if sim_ism:
                self.init_ism()
                
            if sim_scint:
                self.init_scint()
                
            if sim_telescope:
                self.init_telescope()
               
            
        
            
        
        
        
    def init_signal(self,signal_dict = None): 
        """This is a function that intiatlizes a signal class using a dictionary of parameters. It either uses the
        dictionary set when Simulation() is called, when signal_dict is left as None, or a dictianary that is input directly
        into the init_signal function. 
        
        
        @param f0 -- central frequency (MHz)
        @param bw -- bandwidth (MHz)
        @param Nf -- number of freq. bins
        @param f_samp -- samplinf frequency (MHz)  
        @param ObsTime --  total time of observation (add units once confirmed)
        @param data_type -- 'int8' or 'int16' supported. Automatically changed to 'uint8' or 'uint16' if intensity signal.          
        @param SignalType -- 'intensity' which carries a Nf x Nt filterbank of pulses or 'voltage'
                which carries a 4 x Nt array of voltage vs. time pulses representing 4 stokes channels
        @param mode -- 'simulate' beccause this is the simulator mode"""
        
        d = self._check_simul_dict(signal_dict)
        self.signal = PSS.Signal(f0 = d['f0'],bw = d['bw'], Nf = d['Nf'], f_samp = d['f_samp'] , ObsTime = d['ObsTime']\
                                    ,data_type = d['data_type'],SignalType = d['SignalType'], mode = 'simulate')        
        
    
    def init_pulsar(self,pulsar_dict = None):
        """This is a function that initializes the pulsar class using a dictionary of paramaters.It either uses the
        dictionary set when Simulation() is called, when pulsar_dict is left as None, or a dictianary that is input directly
        into the init_pulsar function.
        
        @param F0 -- Pulsar Spin Frequency (Hz)
        @param flux -- mean flux density of pulsar in mJy"""
        
        d =  self._check_simul_dict(pulsar_dict)
        self.pulsar = PSS.Pulsar(self.signal, period = (1/d['F0'] * 1000), flux = d['flux']) # in units of milliseconds
        
        
    def init_ism(self,ism_dict = None):
        """This is a function that initializes the pulsar class using a dictionary of paramaters.It either uses the dictionary
        set when Simulation() is called, when ism_dict is left as None, or a dictianary that is input directly
        into the init_ism function.
        
        @param dm -- dispersion measure"""
        d = self._check_simul_dict(ism_dict)
        self.ISM = PSS.ISM(self.signal, DM = d['dm']) 
        
        
        
        
    
    def init_scint(self,scint_dict = None):
        """This is a function that initializes the pulsar class using a dictionary of paramaters.It either uses the dictionary
        set when Simulation() is called, when scint_dict is left as None, or a dictianary that is input directly
        into the init_scint function.
        
        @param V_ISS -- Intersteller Scintilation Velocity
        @param scint_bw -- scintilation Bandwidth
        @param scint_timescale -- scintilation timescale
        @param pulsar -- pulsar name
        @param to_use_NG_pulsar -- use NG pulsar (True/False)
        @param telescope -- telescope name(GBT or Arecibo)
        @param freq_band -- frequency band"""
        
        d = self._check_simul_dict(scint_dict)
        
        try:
            self.Scint = PSS.scintillate(self.signal, V_ISS=d['V_ISS'], scint_bw=d['scint_bw'], scint_timescale=d['scint_timescale'],\
                pulsar=None, to_use_NG_pulsar=False, telescope=None, freq_band=None)
                                
        except KeyError:
            self.Scient = PSS.scintillate(self.signal,  V_ISS=None, scint_bw=None, scint_timescale=None,\
                        pulsar=d['pulsar'], to_use_NG_pulsar=True, telescope=d['telescope'], freq_band=d['freq_band'])
            
            
        
        
    
    
    def init_telescope(self,telescope_dict = None):
        """This is a function that initializes the pulsar class using a dictionary of paramaters.It either uses the dictionary
        set when Simulation() is called, when telescope_dict is left as None, oor a dictianary that is input directly
        into the init_tlescope function. If a GBT or Arecibo is input to sim_telescope then it will be automatically put 
        into the simulation
        
        @param aperature -- aperature (m)
        @param area -- collecting area (m^2)
        @param Tsys -- system temp (K), total of receiver, sky, spillover, etc. (only needed for noise)
        @param name -- GBT or Arecibo. Uses telescope parameter"""
        
        if self.sim_telescope == 'GBT':
            self.telescope = PSS.GBT()
            
        elif self.sim_telescope == 'Arecibo':
            self.telescope = PSS.Arecibo()
            
        elif self.sim_telescope == None: 
            
            d = self._check_simul_dict(telescope_dict)
            self.telescope = PSS.Telescope(d['aperture'],area = d['area'], Tsys = d['Tsys'], name = d['telescope'])
            
        
        
        
        
        
        
    def simulate(self):
        """This is the fucntion that runs the simulation once the desired classes are initialized."""
        if hasattr(self,'ISM'):
            if 'tau_scatter' in self.sim_dict.keys():
                self.ISM.to_Scatter_Broaden_exp = True
                self.ISM.tau_scatter = self.sim_dict['tau_scatter']
               
            if 'DM' and 'to_DM_Broaden' in self.sim_dict.keys():
                self.ISM.to_DM_Broaden = self.sim_dict['to_DM_Broaden']
                

                
            
            self.ISM.finalize_ism()
           
            if self.signal.MetaData.to_DM_Broaden:
                tophats = PSS.ism.make_dm_broaden_tophat(self.pulsar,self.signal)
                PSS.ism.convolve_with_profile(self.pulsar,tophats)

            if self.signal.MetaData.to_Scatter_Broaden_exp:
                exponentials = PSS.ism.make_scatter_broaden_exp(self.pulsar,self.signal,self.sim_dict['tau_scatter'])
                PSS.ism.convolve_with_profile(self.pulsar,exponentials)

    
        self.pulsar.make_pulses()
       
        
        try:
            self.ISM.disperse()
            
        except AttributeError:
            pass
        
#         try:
#             self.Scintillate = self.Scint.() # Add command once it exists
#         except AttributeError:
#             pass
        
        if hasattr(self,'telescope'): # if sim_telescope is not true it skips
            if self.sim_dict['radiometer_noise']: 
                self.obs_signal = self.telescope.radiometer_noise(self.signal,[self.signal.Nf, self.signal.Nt],self.signal.TimeBinSize)
              
        
    
    
    
                
        
                                    
####Convenience Methods
                                    
    def _check_simul_dict(self,d):
        """ This is used to make the init methods easier to write.
            Takes the init methed and checks to see if there is an input dictionary and sets it
            as a dictionary to be used in the initialization of the class"""                        
        if not d:
            if self.sim_dict:  
                d = self.sim_dict
            else:
                raise ValueError('Parameter dictionary needed.')           
        elif isinstance(d,dict):
            pass
        else:
            raise TypeError('Input not a dictionary')
                                    
                                    
        return d
    
    
    def _set_pulsar_dict(self):
        """Takes Pulsar name and retrieves the following parameters:
             
        @param F0 -- Pulsar Spin Frequency (Hz)  
        @param dm -- dispersion measure 
        @param scint_bw -- scintilation Bandwidth            
        @param scint_timescale -- scintilation timescale   
         """
    
        p = PSS.PSS_utils.get_pint_models(self.psr, self.sim_file_path) 
        param_dict = {'F0': p.F0.value, 'dm':p.DM.value} 
        [scint_bw, scint_timescale] = PSS.ism.NG_scint_param(self.psr, self.sim_telescope , self.sim_dict['freq_band']) 
        param_dict['scint_bw'] = scint_bw
        param_dict['scint_timescale'] = scint_timescale 
        
        
        return param_dict



    def add_scint(self):
        """
        This is old code and will be superseded by adding the scintillations into
        the DAT_SCL arrays of the SUB_INT HDU in the SEARCH mode PSRFITS file.
        We keep it here until we have fully tested how that functionality works
        with folding the PSRFITS file.
        """
        scint_time = self.MD.scint_time/self.MD.scint_time_sample_rate
        scint_samples_per_obs = np.floor(self.MD.ObsTime//(scint_time*1e3))
        #print('scint_samples_per_obs',scint_samples_per_obs)
        #gain_norm = self.scint_class.gain.max() #Could set to avoid clipping, but not sure it's needed.
        gain = self.scint_class.gain# / gain_norm
        scint_end_bin = scint_samples_per_obs * scint_time*1e3 #integer number of bins in scint
        #print('scint_end_bin',scint_end_bin)
        self.start_times = np.linspace(0, scint_end_bin, scint_samples_per_obs)
        orig_profile = np.copy(self.P.profile) #* P.gamma_draw_norm
        scint_bins = int(scint_time//self.S.TimeBinSize)
        tweak = 12
        if len(self.start_times) > len(gain[0,:]):
            raise ValueError('Scattering Screen is not long enough to scintillate at this Dispersion timescale.')
        for ii, bin_time in enumerate(self.start_times) :
            self.P.profile = gain[:, ii, np.newaxis] * orig_profile * self.P.gamma_draw_norm
            #self.P.profile /= (self.P.profile.max()/tweak)
            self.P.make_pulses(bin_time, bin_time + scint_time)
            #bin = int(bin_time //self.S.TimeBinSize)
            #self.S.signal[:,bin : bin + scint_bins] = gain[:,ii,np.newaxis]*self.S.signal[:,bin: bin + scint_bins]
        #self.frig = gain[:,ii,np.newaxis] * orig_profile

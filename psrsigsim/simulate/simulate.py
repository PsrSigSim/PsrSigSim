from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
import numpy as np
from astropy import log

from ..signal import FilterBankSignal
from ..telescope import Telescope, Receiver, Backend
from ..telescope import telescope
from ..pulsar import Pulsar
from ..pulsar import GaussPortrait
from ..pulsar import UserPortrait
from ..pulsar import DataPortrait
from ..pulsar import DataProfile
from ..ism import ISM
from ..io import PSRFITS, TxtFile
from ..utils.utils import make_par

class Simulation(object):
    """convenience class for full simulations.

    Necessary information includes all minimal parameters
    for instances of each other class, Signal, Pulsar, ISM,
    Telescope.

    Input may be specified manually, from a pre-made parfile with
    additional input, e.g. for the Signal, or from a premade dictionary
    with appropriate keys.

    Parameters
    ----------
    fcent : float]
        Central radio frequency (MHz)
    bandwidth : float
        Radio bandwidth of signal (MHz)
    Nsubband : int
        Number of sub-bands, default ``512``
        XUPPI backends use 2048 frequency channels divided between the
        four Stokes parameters, so 512 per Stokes parameter.
    sample_rate : float
        Sample rate of data (MHz), default: ``None``
        If no ``sample_rate`` is given the observation will default to
        the 20.48 us per sample (about 50 kHz).  This is the sample rate
        for coherently dedispersed filter banks using XUPPI backends.
    sublen : float
        Desired length of data subintegration (sec) if subint
        is ``True``, default: ``tobs``. If left as none but subint is
        ``True``, then when pulses are made, the sublen will default to
        the input observation length, ``tobs``
    dtype : type
        Data type of array, default: ``np.float32``
        supported types are: ``np.float32`` and ``np.int8``
    fold : bool
        If `True`, the initialized signal will be folded to some
        number of subintegrations based on sublen (else will just make
        a single subintegration). If `False`, the data produced will be
        single pulse filterbank data. Default is `True`.
        NOTE - using `False` will generate a large amount of data.
    period : float
        Pulse period (sec)
    Smean : float
        Mean pulse flux density (Jy)
    profile : array or function or Pulse Profile/Portrait Class
        Pulse profile or 2-D pulse portrait, this can take four forms:
            array - either an array of Gaussian components in the order
                    [peak phase, width, amplitude]
                    OR
                    A data array representative of the pulse profile, or
                    samples of the profile from phases between 0 and 1.
                    Data array must have more than 3 points.
            function - function defining the shape of the profile given a
                       of input phases. CURRENTLY NOT IMPLEMENTED.
            class - predefined PsrSigSim Pulse Profile or Pulse Portrait class
                    object.
    specidx : float
        Spectral index of the pulsar. Default value is 0 (i.e. no spectral index).
        
    ref_freq : float
        The reference frequency of the input value of Smean in MHz. The default
        value will be the center frequency of the bandwidth.
    tobs : float
        Total simulated observing time in seconds
    name : str
        Name of pulsar
    dm : float
        Dispersion measure of the pulsar (pc cm^-3)
    tau_d : float
        Scattering timescale to use (s)
    tau_d_ref_f : float
        reference frequency for the input scattering timescale (MHz)
    aperture : float
        Telescop aperture (m)
    area : float
        Collecting area (m^2) (if omitted, assume circular single dish)
    Tsys : float
        System temperature (K) of the telescope (if omitted use Trec)
    tscope_name : string
        Name of the telescope. If GBT or Arecibo, will use predefined parameters.
    system_name : string
        Name of telescope system, backend-recviever combination. May be a list.
    rcvr_fcent: float
        Center frequency of the telescope reciever. May be a list.
    rcvr_bw : float
        Bandwidth of the telescope reciever. May be a list.
    rcvr_name : string
        Name of the telescope reciever. May be a list.
    backend_samprate : float
        Sampling rate (in MHz) of the telescope backend. May be a list.
    backend_name : string
        Name of the telescope backend. May be a list.
    tempfile : string
        Path to template psrfits file to use for saving simulated data.
    parfile : string
        Path to pulsar par file to read in to use for pulsar parameters
    psrdict : dictionary
        Dictionary of input parameters to generate simualted data from.
        Keys should be the same as possible input values listed above.

    """
    def __init__(self,
                 fcent = None,
                 bandwidth = None,
                 sample_rate = None,
                 dtype = np.float32,
                 Npols = 1,
                 Nchan = 512,
                 sublen = None,
                 fold = True,
                 period = None,
                 Smean = None,
                 profiles = None,
                 specidx = 0.0,
                 ref_freq = None,
                 tobs = None,
                 name = None,
                 dm = None,
                 tau_d = None,
                 tau_d_ref_f = None,
                 aperture = None,
                 area = None,
                 Tsys = None,
                 tscope_name = None,
                 system_name = None,
                 rcvr_fcent = None,
                 rcvr_bw = None,
                 rcvr_name = None,
                 backend_samprate = None,
                 backend_name = None,
                 tempfile = None,
                 parfile = None,
                 psrdict = None,):

        # Assign values from manual input
        self._fcent = fcent
        self._bandwidth = bandwidth
        self._sample_rate = sample_rate
        self._Npols = Npols
        self._Nchan = Nchan
        self._sublen = sublen
        self._fold = fold
        self._period  = period
        self._Smean = Smean
        self._profiles = profiles
        self._specidx = specidx
        self._ref_freq = ref_freq
        self._tobs = tobs
        self._name = name
        self._dm = dm
        self._tau_d = tau_d
        self._tau_d_ref_f = tau_d_ref_f
        self._aperture = aperture
        self._area = area
        self._Tsys = Tsys
        self._tscope_name = tscope_name
        self._system_name = system_name
        self._rcvr_fcent = rcvr_fcent
        self._rcvr_bw = rcvr_bw
        self._rcvr_name = rcvr_name
        self._backend_samprate = backend_samprate
        self._backend_name = backend_name
        self._tempfile = tempfile
        # Assign values from parfile
        if parfile != None:
            self.params_from_par(parfile)
        # Assign values from dictionary
        if psrdict != None:
            self.params_from_dict(psrdict)

    def params_from_dict(self, psrdict):
        """
        Function to take the input dictionary and assign values from that.
        """
        for key in psrdict.keys():
            setattr(self, "_"+key, psrdict[key])

    def params_from_par(self, parfile):
        """
        Function to take input par file and assign values from that.
        """
        raise NotImplementedError()

    def init_signal(self, from_template = False):
        """
        Function to initialize a signal from the input parameters.

        Parameters
        ----------
        from_template : bool
            If True, will use the input template file to initialize the
            signal. If False will use other input values to initialize the signal.
        """
        if from_template:
            pfit = PSRFITS(path="sim_fits.fits", template=self.tempfile, fits_mode='copy', \
                              obs_mode='PSR')
            sim_signal = pfit.make_signal_from_psrfits()
            self._signal = sim_signal
        else:
            sim_signal = FilterBankSignal(fcent = self.fcent, bandwidth = self.bw, Nsubband=self.Nchan,\
                                sample_rate=self.samprate, fold=self.fold, sublen=self.sublen)
            self._signal = sim_signal

    def init_profile(self):
        """
        Function to initialize a profile object from input.
        """
        # First check if the profile is a class
        proftypes = [GaussPortrait, UserPortrait, DataPortrait, DataProfile]
        if any(isinstance(self.profiles, x) for x in proftypes):
            pass
        elif isinstance(self.profiles, (list, np.ndarray)):
            if len(self.profiles) == 3:
                prof = GaussPortrait(peak = self.profiles[0], width = self.profiles[1], amp = self.profiles[2])
            elif len(self.profiles) > 3:
                prof = DataProfile(self.profiles, phases = None, Nchan=self.Nchan)
            else:
                raise RuntimeError("Input profile array has too few values!")
        elif callable(self.profiles):
            raise NotImplementedError()
        else:
            log.warning("Unrecognized input profile type, defaulting to Gaussian.")
            prof = GaussPortrait()
        # Reassign the profile as the class object
        if not any(isinstance(self.profiles, x) for x in proftypes):
            self._profiles = prof


    def init_pulsar(self):
        """
        Function to initialize a pulsar from the input parameters.
        NOTE - Must have initialized the profile before running this.
        """
        # Define the pulsar
        pulsar = Pulsar(period=self.period, Smean=self.Smean, \
                                   profiles=self.profiles, name=self.name, \
                                   specidx=self.specidx, ref_freq=self.ref_freq)
        self._pulsar = pulsar

    def init_ism(self):
        """
        Function to initialize the ISM from the input parameters.
        """
        ism = ISM()
        self._ism = ism

    def init_telescope(self):
        """
        Function to initialize the telescope from input parameters.
        """
        if self.tscope_name == 'GBT':
            tscope = telescope.GBT()
        elif self.tscope_name == 'Arecibo':
            tscope = telescope.Arecibo()
        else:
            tscope = Telescope(self.aperture, area = self.area, Tsys = self.Tsys,
                               name = self.tscope_name)
        # Now need to add systems
        if type(self.rcvr_fcent) is list:
            if not len(self.system_name) == len(self.rcvr_fcent) == len(self.rcvr_bw) == \
            len(self.rcvr_name) == len(self.backend_samprate) == len(self.backend_name):
                raise RuntimeError("Number of telescope system entries do not match!")
            # Now make the systems
            for ii in range(len(self.rcvr_fcent)):
                tscope.add_system(name=self.system_name[ii],
                 receiver=Receiver(fcent=self.rcvr_fcent[ii], bandwidth=self.rcvr_bw[ii], name=self.rcvr_name[ii]),
                 backend=Backend(samprate=self.backend_samprate[ii], name=self.backend_name[ii]))
        elif self.rcvr_fcent != None:
            tscope.add_system(name=self.system_name,
                 receiver=Receiver(fcent=self.rcvr_fcent, bandwidth=self.rcvr_bw, name=self.rcvr_name),
                 backend=Backend(samprate=self.backend_samprate, name=self.backend_name))
        # If not entered, do not add any systems
        self._tscope = tscope

    def simulate(self, from_template = False):
        """
        Function to run the full simulation.

        Parameters
        ----------
        from_template : bool
            If True, will use the input template file to initialize the
            signal. If False will use other input values to initialize the signal.
        twoD : bool
            If True, will generate a 2-D profile array, else will do a
            1-D and will tile the profile in frequency.
        """
        # We start by initializing things
        self.init_signal(from_template = False)
        # Initialize the profile
        self.init_profile()
        # Now the pulsar
        self.init_pulsar()
        # Now the ISM
        self.init_ism()
        # Now check if profiles need to be convolved for scatter broadening.
        if self.tau_d != None:
            self.ism.scatter_broaden(self.signal, self.tau_d, self.tau_d_ref_f, convolve=True, pulsar=self.pulsar)
        # Now make the pulses
        self.pulsar.make_pulses(self.signal, tobs = self.tobs)
        # disperse the simulated pulses
        self.ism.disperse(self.signal, self.dm)
        # TODO: Add in FD parameter or other delays in ISM

        # Now add the telescope and radiometer noise
        self.init_telescope()
        # add radiometer noise
        # TODO: figure out how to do multiple systems
        out_array = self.tscope.observe(self.signal, self.pulsar, system=self.system_name, noise=True)

    def save_simulation(self, outfile = "simfits", out_format = 'psrfits',
                        parfile = None, ref_MJD = 56000.0, MJD_start = 55999.9861) :
        """
        Function to save the simulated data in a default format.
        Currently only PSRFITS is supported.

        Parameters
        ----------
        outfile : string
            Path and name of output save file.
            If not provided, output file is "simfits".
        out_format : string
            Format of output file (not case sensitive).
            Options are:
            'psrfits' - PSRFITS format. Requires template file.
            'pdv' - PSRCHIVE pdv format. Output is a text file.
        parfile : string
            Parfile to use to make phase connection polycos. If none supplied
            will attempt to create one.
        ref_MJD : float
            Reference MJD for phase connection.
        MJD_start : float
            Desired start time of the simulated observation. Needed for
            phase connection.
        """
        # Now save the data if desired
        if out_format.lower() == 'psrfits':
            if outfile == 'simfits':
                outfile += ".fits"
            if self.tempfile == None:
                raise RuntimeError("No template PSRFITS file provided.")
            else:
                pfit = PSRFITS(path=outfile, template=self.tempfile, fits_mode='copy', \
                              obs_mode='PSR')
                pfit.get_signal_params(signal = self.signal)
                # Now save the data
                if parfile == None:
                    log.warning("No par file provided, attempting to make one...")
                    make_par(self.signal, self.pulsar, outpar = "simpar.par")
                    parfile = "simpar.par"
                pfit.save(self.signal, self.pulsar, parfile = parfile, \
                          MJD_start = MJD_start, segLength = 60.0,\
                          ref_MJD = ref_MJD, usePint = True)
        elif out_format.lower() == 'pdv':
            if outfile == 'simfits':
                outfile += ".ar"
            txtfile = TxtFile(path=outfile)
            txtfile.save_psrchive_pdv(self.signal, self.pulsar)
        else:
            raise RuntimeError("Unrecognized output file format: %s" % (out_format))



    @property
    def fold(self):
        return self._fold

    @property
    def sublen(self):
        return self._sublen

    @property
    def Nchan(self):
        return self._Nchan

    @property
    def fcent(self):
        return self._fcent

    @property
    def bw(self):
        return self._bandwidth

    @property
    def tobs(self):
            return self._tobs

    @property
    def samprate(self):
        return self._sample_rate

    @property
    def dtype(self):
        return self._dtype

    @property
    def Npols(self):
        return self._Npols

    @property
    def dm(self):
        return self._dm

    @property
    def tau_d(self):
        return self._tau_d

    @property
    def tau_d_ref_f(self):
        return self._tau_d_ref_f

    @property
    def profiles(self):
        return self._profiles

    @property
    def name(self):
        return self._name

    @property
    def period(self):
        return self._period

    @property
    def Smean(self):
        return self._Smean
    
    @property
    def specidx(self):
        return self._specidx
    
    @property
    def ref_freq(self):
        return self._ref_freq

    @property
    def tscope_name(self):
        return self._tscope_name

    @property
    def area(self):
        return self._area

    @property
    def aperture(self):
        return self._aperture

    @property
    def Tsys(self):
        return self._Tsys

    @property
    def system_name(self):
        return self._system_name

    @property
    def rcvr_fcent(self):
        return self._rcvr_fcent

    @property
    def rcvr_bw(self):
        return self._rcvr_bw

    @property
    def rcvr_name(self):
        return self._rcvr_name

    @property
    def backend_samprate(self):
        return self._backend_samprate

    @property
    def backend_name(self):
        return self._backend_name

    @property
    def tempfile(self):
        return self._tempfile

    @property
    def signal(self):
        return self._signal

    @property
    def pulsar(self):
        return self._pulsar

    @property
    def ism(self):
        return self._ism

    @property
    def tscope(self):
        return self._tscope

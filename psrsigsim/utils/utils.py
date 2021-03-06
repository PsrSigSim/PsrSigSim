"""PSS_utils.py
A place to organize methods used by multiple modules
"""
from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
import numpy as np
import scipy as sp
from astropy import units as u
from pint import models
#try:
#    import pyfftw
#    use_pyfftw = True
#except:
#    use_pyfftw = False


def shift_t(y, shift, dt=1):
    """
    Shift timeseries data in time.
    Shift array, y, in time by amount, shift. For dt=1 units of samples
    (including fractional samples) are used. Otherwise, shift and dt are
    assumed to have the same physical units (i.e. seconds).

    Parameters
    ----------
    y : array like, shape (N,), real
        Time series data.
    shift : int or float
        Amount to shift
    dt : float
        Time spacing of samples in y (aka cadence).

    Returns
    -------
    out : ndarray
        Time shifted data.

    Examples
    --------
    >>>shift_t(y, 20)
    # shift data by 20 samples

    >>>shift_t(y, 0.35, dt=0.125)
    # shift data sampled at 8 Hz by 0.35 sec

    Uses np.roll() for integer shifts and the Fourier shift theorem with
    real FFT in general.  Defined so positive shift yields a "delay".
    """
    if isinstance(shift, int) and dt == 1:
        out = np.roll(y, shift)

    else:
        yfft = np.fft.rfft(y)  # hermicity implicitely enforced by rfft
        fs = np.fft.rfftfreq(len(y), d=dt)
        phase = -1j*2*np.pi*fs*shift
        yfft_sh = yfft * np.exp(phase)
        out = np.fft.irfft(yfft_sh)

    return out


def down_sample(ar, fact):
    """down_sample(ar, fact)
    down sample array, ar, by downsampling factor, fact
    """
    #TODO this is fast, but not as general as possible
    downsampled = ar.reshape(-1, fact).mean(axis=1)
    return downsampled


def rebin(ar, newlen):
    """rebin(ar, newlen)
    down sample array, ar, to newlen number of bins
    This is a general downsampling rebinner, but is slower than down_sample().
    'ar' must be a 1-d array
    """
    newBins = np.linspace(0, ar.size, newlen, endpoint=False)
    stride = newBins[1] - newBins[0]
    maxWid = int(np.ceil(stride))
    ar_new = np.empty((newlen, maxWid))  # init empty array
    ar_new.fill(np.nan)  # fill with NaNs (no extra 0s in mean)

    for ii, lbin in enumerate(newBins):
        rbin = int(np.ceil(lbin + stride))
        # fix for potential last bin rounding error
        if rbin > ar.size: # Not sure how to force this for test...
            rbin = ar.size
        lbin = int(np.ceil(lbin))
        ar_new[ii, 0:rbin-lbin] = ar[lbin:rbin]

    return np.nanmean(ar_new, axis=1)  # ingnore NaNs in mean


def top_hat_width(subband_df, subband_f0, DM):
    """top_hat_width(subband_df, subband_f0, DM)
    Returns width of a top-hat pulse to convolve with pulses for dipsersion
    broadening. Following Lorimer and Kramer, 2005 (sec 4.1.1 and A2.4)
    subband_df : subband bandwidth (MHz)
    subband_f0 : subband center frequency (MHz)
    DM : dispersion measure (pc/cm^3)
    return top_hat_width (milliseconds)
    """
    D = 4.148808e3  # sec*MHz^2*pc^-1*cm^3, dispersion const
    width_sec = 2*D * DM * (subband_df) / (subband_f0)**3
    return width_sec * 1.0e+3  # ms


def savitzky_golay(y, window_size, order, deriv=0, rate=1):
    # courtesy scipy recipes
    r"""Smooth (and optionally differentiate) data with a Savitzky-Golay filter.
    The Savitzky-Golay filter removes high frequency noise from data.
    It has the advantage of preserving the original shape and
    features of the signal better than other types of filtering
    approaches, such as moving averages techniques.
    Parameters
    ----------
    y : array_like, shape (N,)
        the values of the time history of the signal.
    window_size : int
        the length of the window. Must be an odd integer number.
    order : int
        the order of the polynomial used in the filtering.
        Must be less then `window_size` - 1.
    deriv: int
        the order of the derivative to compute
        (default = 0 means only smoothing)
    Returns
    -------
    ys : ndarray, shape (N)
        the smoothed signal (or it's n-th derivative).
    Notes
    -----
    The Savitzky-Golay is a type of low-pass filter, particularly
    suited for smoothing noisy data. The main idea behind this
    approach is to make for each point a least-square fit with a
    polynomial of high order over a odd-sized window centered at
    the point.
    Examples
    --------
    t = np.linspace(-4, 4, 500)
    y = np.exp( -t**2 ) + np.random.normal(0, 0.05, t.shape)
    ysg = savitzky_golay(y, window_size=31, order=4)
    import matplotlib.pyplot as plt
    plt.plot(t, y, label='Noisy signal')
    plt.plot(t, np.exp(-t**2), 'k', lw=1.5, label='Original signal')
    plt.plot(t, ysg, 'r', label='Filtered signal')
    plt.legend()
    plt.show()
    References
    ----------
    .. [1] A. Savitzky, M. J. E. Golay, Smoothing and Differentiation of
       Data by Simplified Least Squares Procedures. Analytical
       Chemistry, 1964, 36 (8), pp 1627-1639.
    .. [2] Numerical Recipes 3rd Edition: The Art of Scientific Computing
       W.H. Press, S.A. Teukolsky, W.T. Vetterling, B.P. Flannery
       Cambridge University Press ISBN-13: 9780521880688
    """
    from math import factorial

    try:
        window_size = np.abs(np.int(window_size))
        order = np.abs(np.int(order))
    except TypeError:  # ValueError, msg:
        raise ValueError("window_size and order have to be of type int")
    if window_size % 2 != 1 or window_size < 1:
        raise TypeError("window_size size must be a positive odd number")
    if window_size < order + 2:
        raise TypeError("window_size is too small for the polynomials order")
    order_range = range(order+1)
    half_window = int((window_size -1) // 2)
    # precompute coefficients
    b = np.array([[k**i for i in order_range] for k in range(-half_window,
                                                           half_window+1)])
    m = np.linalg.pinv(b)[deriv] * rate**deriv * factorial(deriv)
    # pad the signal at the extremes with
    # values taken from the signal itself
    firstvals = y[0] - np.abs(y[1:half_window+1][::-1] - y[0])
    lastvals = y[-1] + np.abs(y[-half_window-1:-1][::-1] - y[-1])
    y = np.concatenate((firstvals, y, lastvals))
    return np.convolve(m[::-1], y, mode='valid')


def find_nearest(array, value):
    """Returns the argument of the element in an array nearest to value.
    For half width at value use array[max:].
    """
    diff=np.abs(array-value)
    idx = diff.argmin()
    if idx == 0 or array[1] < value:
        idx = 1
    return idx


def acf2d(array, speed='fast', mode='full', xlags=None, ylags=None):
    """Courtesy of Michael Lam's PyPulse
    Calculate the autocorrelation of a 2 dimensional array.
    """
    from scipy.signal import fftconvolve, correlate

    if speed == 'fast' or speed == 'slow':
        ones = np.ones(np.shape(array))  # very close for either speed
        norm = fftconvolve(ones, ones, mode=mode)
        if speed=='fast':
            return fftconvolve(array, np.flipud(np.fliplr(array)),
                               mode=mode)/norm
        else:
            return correlate(array, array, mode=mode)/norm
    elif speed == 'exact':
        #NOTE: (r,c) convention is flipped from (x,y),
        # also that increasing c is decreasing y
        LENX = len(array[0])
        LENY = len(array)
        if xlags is None:
            xlags = np.arange(-1*LENX+1, LENX)
        if ylags is None:
            ylags = np.arange(-1*LENY+1, LENY)
        retval = np.zeros((len(ylags), len(xlags)))
        for i, xlag in enumerate(xlags):
            print(xlag)
            for j, ylag in enumerate(ylags):
                if ylag > 0 and xlag > 0:
                    A = array[:-1*ylag, xlag:]  # the "stationary" array
                    B = array[ylag:, :-1*xlag]
                elif ylag < 0 and xlag > 0:
                    A = array[-1*ylag:, xlag:]
                    B = array[:ylag, :-1*xlag]
                elif ylag > 0 and xlag < 0:  # optimize later via symmetries
                    A = array[:-1*ylag, :xlag]
                    B = array[ylag:, -1*xlag:]
                elif ylag < 0 and xlag < 0:
                    A = array[-1*ylag:, :xlag]
                    B = array[:ylag, -1*xlag:]
                else:  # one of the lags is zero
                    if ylag == 0 and xlag > 0:
                        A = array[-1*ylag:, xlag:]
                        B = array[:, :-1*xlag]
                    elif ylag == 0 and xlag < 0:
                        A = array[-1*ylag:, :xlag]
                        B = array[:, -1*xlag:]
                    elif ylag > 0 and xlag == 0:
                        A = array[:-1*ylag, :]
                        B = array[ylag:, -1*xlag:]
                    elif ylag < 0 and xlag == 0:
                        A = array[-1*ylag:, :]
                        B = array[:ylag, -1*xlag:]
                    else:
                        A = array[:, :]
                        B = array[:, :]
                        #print xlag,ylag,A,B
                C = A*B
                C = C.flatten()
                goodinds = np.where(np.isfinite(C))[0]  # check for good values
                retval[j, i] = np.mean(C[goodinds])
        return retval


def text_search(search_list, header_values, filepath, header_line=0,
                file_type='txt'):
    """ Method for pulling value from  a txt file.
    search_list = list of string-type values that demarcate the line in a txt
                  file from which to pull values
    header_values = string of column headers or array of column numbers
                    (in Python numbering) the values from which to pull
    filepath = file path of txt file. (string)
    header_line = line with headers for values.
    file_type = 'txt' or 'csv'

    returns: tuple of values matching header values for the search terms given.
    """
    #TODO Make work for other file types.
    #if file_type == 'txt':
    #    delimiter = ''
    #elif file_type == 'csv':
    #    delimiter = ','

    check = 0
    output_values = list()

    with open(filepath, 'r') as f:  # read file to local memory
        searchfile = f.readlines()

    # Find Column Numbers from column names
    if any(isinstance(elem, str) for elem in header_values):
        column_num = []
        parsed_header = searchfile[header_line].split()
        for ii, header in enumerate(header_values):
            column_num.append(parsed_header.index(header))
    else:
        column_num = np.array(header_values)

    # Find Values using search keys and column numbers.
    for line in searchfile:
        if all(ii in line for ii in search_list):

            info = line.split()
            for jj, value in enumerate(column_num):
                output_values.append(info[value])
            check += 1

    if check == 0:
        raise ValueError('Combination {0} '.format(search_list)+' not found in \
                            same line of text file.')
    if check > 1:
        raise ValueError('Combination {0} '.format(search_list)+' returned \
                            multiple results in txt file.')

    return tuple([float(i) for i in output_values])


def make_quant(param, default_unit):
    """
    Convenience function to intialize a parameter as an astropy quantity.

    Parameters
    ----------
    param : attribute
        Parameter to initialize.
    default_unit : string
        Name of an astropy unit, set as default for this parameter.

    Returns
    -------
        An astropy quantity

    Examples
    --------
        self.f0 = make_quant(f0,'MHz')
    """
    default_unit = u.core.Unit(default_unit)
    if hasattr(param, 'unit'):
        try:
            param.to(default_unit)
        except u.UnitConversionError:
            raise ValueError("Quantity {0} with incompatible unit {1}"
                             .format(param, default_unit))
        quantity = param
    else:
        quantity = param * default_unit

    return quantity

def get_pint_models(psr_name, psr_file_path):
        """Function that returns pint model given a specific pulsar"""
        # will need to add section for J1713 T2 file. gls is not file wanted for this specfic pulsar.
        model_name = "{0}{1}_NANOGrav_11yv1.gls.par".format(psr_file_path,psr_name)
        par_model = models.get_model(model_name)

        return par_model

def make_par(signal, pulsar, outpar = "simpar.par"):
    """
    Function to create a par file for simulated pulsar.
    TO DO: Will need to update when additional delays are added

    Parameters
    ----------
    signal : class
        PsrSigSim Signal class object
    pulsar : class
        PsrSigSim Pulsar class object
    outpar : string
        Name of output par file.
    """
    # Get parameters and other things that should go into this file
    par_lines = []
    par_lines.append("PSR            %s\n" % (pulsar.name))
    par_lines.append("LAMBDA            0.0\n" ) # Default for now
    par_lines.append("BETA           0.0\n" ) # Default for now
    par_lines.append("PMLAMBDA            0.0\n" ) # Default for now
    par_lines.append("PMBETA            0.0\n" ) # Default for now
    par_lines.append("PX            0.0\n" ) # Default for now
    par_lines.append("POSEPOCH            56000.0\n" ) # Default for now
    par_lines.append("F0           %s\n" % (1.0/pulsar.period))
    par_lines.append("PEPOCH            56000.0\n" ) # Default for now
    par_lines.append("START            50000.0\n" ) # Default for now
    par_lines.append("FINISH            60000.0\n" ) # Default for now
    par_lines.append("DM                %s\n" % (signal.dm))
    par_lines.append("EPHEM               DE436\n" ) # Default for now
    par_lines.append("ECL                 IERS2010\n" ) # Default for now
    par_lines.append("CLK                 TT(BIPM2015) \n" ) # Default for now
    par_lines.append("UNITS               TDB\n" ) # Default for now
    par_lines.append("TIMEEPH             FB90\n" ) # Default for now
    par_lines.append("T2CMETHOD           TEMPO\n" ) # Default for now
    par_lines.append("CORRECT_TROPOSPHERE N\n" ) # Default for now
    par_lines.append("PLANET_SHAPIRO      N\n" ) # Default for now
    par_lines.append("DILATEFREQ          N\n" ) # Default for now
    par_lines.append("TZRMJD        56000.0\n" ) # Default for now
    par_lines.append("TZRFRQ            1500.0\n" ) # Default for now
    par_lines.append("TZRSITE                  @\n" ) # Default for now
    par_lines.append("MODE                     1\n" ) # Default for now
    # Write out the file
    with open(outpar, 'w') as op:
        op.writelines(par_lines)
        op.close()

"""
utils.py
A place to organize methods used by multiple modules
"""

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
import numpy as np
import scipy as sp
import math

def shiftit(y, shift):
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

    work.real[size//2:] = work.real[size//2:0:-1]
    work.imag[size//2:] = -work.imag[size//2:0:-1]
    work[size//2] = 0.+0.j
    workifft = np.fft.ifft(work)
    return workifft.real

def down_sample(x, R): #Method to down sample an array by a factor
    #down_sample(array, downsampling factor)
    #This is fast, but not as general as possible

    x.reshape(-1, R)
    downsampled = x.reshape(-1, R).mean(axis=1)#/np.amax(x)

    return downsampled#*np.amax(x)/np.amax(downsampled)

def rebin(a, newLength):
    """rebin(old array, new number of bins)
    This is a very general downsampling rebinner, but has for loops and if
    statements, hence it is slower than down_sample().
    """
    #TODO Make this code run faster. Vectorize
    newBins = np.linspace(0, a.size, newLength, endpoint=False)
    width = math.ceil(a.size/newLength)
    a_rebin=np.zeros((newLength,width))*np.nan
    #Using NaN means that we do not have extra zeros in the array that would get averaged
    row = 0
    column = 0
    for ii in range(0, a.size):
        if ii < (newBins[row] + newBins[1]):
            a_rebin[row,column] = a[ii]
            column +=1
        else:
            column = 0
            row += 1
            a_rebin[row,column] = a[ii]
            column +=1

    a_rebinned = sp.nanmean(a_rebin,axis=1)
    #NaN mean does not count NaNs in total
    return a_rebinned#*np.amax(a)/np.amax(a_rebinned)

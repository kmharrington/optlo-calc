import numpy as np

import astropy.units as u
import astropy.constants as const

def through_put(freqs):
    '''
    Throughput of a diffraction limited telesope as a function of
    frequency
    '''
    if not hasattr(freqs, 'unit'):
        freqs *= u.GHz
    return ((const.c/freqs)**2 *u.steradian).to(u.cm**2 * u.steradian)

def power_law_index(params, nu):
    """
    index power la
    """
    n_real, nu0, a, b = params
    ff = (nu.to(u.GHz)).value
    return n_real + 1.0j*(a*(ff/nu0)**b)

def to_dB(x):
    return 10*np.log10(x)


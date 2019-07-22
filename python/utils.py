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

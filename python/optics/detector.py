import numpy as np
import astropy.units as u
import astropy.constants as const

from .optic import Optic

class Detector(Optic):
    '''
    TODO :: Do detectors need temperatures?
    '''
    def __init__(self, name, bandpass, temperature):
        self.bandpass = self._init_value(bandpass, fill_value=(0,0))
        super(Detector, self).__init__(name, 0, 0, temperature, 0, temperature)
    
    def get_transmission(self, freqs):
        return self._get_value(self.bandpass, freqs)
    
    @classmethod
    def from_band_edges(cls, name, edges, temperature, efficiency=1.0):
        band = lambda f: efficiency*tophat(f, edges[0], edges[1])
        return cls(name, band, temperature)
    
    
def tophat(freqs, low, high):
    out = np.zeros( np.shape(freqs) )
    msk = np.all([freqs >= low,
                  freqs < high], axis=0)
    out[msk]=1
    return out
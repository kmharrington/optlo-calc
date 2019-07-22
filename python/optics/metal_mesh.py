import numpy as np
import astropy.units as u
import astropy.constants as const

from .optic import Optic

class Cutoff_Filter(Optic):
    def __init__(self, name, cutoff_freq, absorption, temperature,
                    spill, spill_temperature):
        self.cutoff = cutoff_freq
        self.absorption = self._init_value(absorption)
        self.spill = spill
        super(Cutoff_Filter, self).__init__(name,
                                            self.absorption,
                                            self.reflection,
                                            temperature, spill, spill_temperature)
    def __str__(self):
        s = '{:10}\t&\t{:02}\t&\t'.format(self.name, self.temperature)
        if np.isscalar(self.absorption):
            s += 'Emissivity: {:04}'.format(np.round(self.absorption, 4))
        else:
            s += 'Emissivity: FUNC'
        s += '  '
        s += 'Freq Cutoff: {:02}'.format(self.cutoff)
        s += '\t&\t'
        s += '{:04}\t&\t{:02}'.format(self.spill, self.spill_temperature)
        return s
    
    def reflection(self, freqs):
        val = np.zeros( np.shape(freqs) )
        val[freqs >= self.cutoff] = 1.0-self.absorption-self.spill
        return val
    
class Metal_Mesh_Filter(Optic):
    '''
    Metal Mesh filters are the only ones defined with a transmission
    filename. Because that is the main measurement we have.
    '''
    def __init__(self, name, transmission, absorption, temperature,
                    spill, spill_temperature):
    
        self.meas_trans = self._init_value(transmission)
        self.absorption = self._init_value(absorption)
        self.spill = self._init_value(spill)
        
        if np.isscalar(self.absorption):
            self.reflection = lambda freqs: self._check_values(1 - 
                                                           self.meas_trans(freqs) - 
                                                           self.absorption - 
                                                           self.spill, True)
        else:
            self.reflection = lambda freqs: self._check_values( 1 - 
                                                           self.meas_trans(freqs) - 
                                                           self.absorption(freqs) - 
                                                           self.spill, True)
        super(Metal_Mesh_Filter, self).__init__(name,
                                               self.absorption,
                                               self.reflection,
                                               temperature, spill, spill_temperature)

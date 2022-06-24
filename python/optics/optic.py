import numpy as np
import astropy.units as u
import astropy.constants as const
from astropy.modeling.models import BlackBody

import scipy.interpolate as spint

import optlo_calc.utils as utils

class Optic(object):
    
    def __init__(self, name, absorption, reflection, temperature,
                    spill, spill_temperature):
        '''
        absorption, reflection, and spill can either be a one of:
            filename -- should be in Bolocalc format where 
                        np.loadtext can be used. The first column is frequency
                        (in GHz) and the second column is the value
            number -- a scalar value for all frequencies
            function -- function which takes frequency and returns an array of transmission
        
        
        TODO :: Decide if I've impemented temperature as a function of frequency
                correctly. For now I'm only using it as a scalar.
        TODO :: Actually find out if spill works.
        '''
        self.name = name
        self.absorption = self._init_value(absorption)
        self.reflection = self._init_value(reflection)
        self.temperature = self._init_value(temperature)
        self.spill = self._init_value(spill)
        self.spill_temperature = self._init_value(spill_temperature)
        super(Optic, self).__init__()
    
    def __str__(self):
        s = '{:10}\t&\t{:02}\t&\t'.format(self.name, self.temperature)
        if np.isscalar(self.absorption):
            s += 'Emissivity: {:04}'.format(np.round(self.absorption, 4))
        else:
            s += 'Emissivity: FUNC'
        s += '  '
        if np.isscalar(self.reflection):
            s += 'Reflection: {:04}'.format(np.round(self.reflection, 4))
        else:
            s += 'Reflection: FUNC'
        s += '\t&\t'
        s += '{:04}\t&\t{:02}'.format(self.spill, self.spill_temperature)
        return s
        
    @staticmethod
    def file_to_function(fname, **kwargs):
        '''
        these files should be in BoloCalc format.
        first column is frequency (in GHz)
        second column is value
        '''
        data = np.loadtxt(fname)
        
        func_params = {
            'bounds_error':False,
            'fill_value': (data[0,1], data[-1,1])
        }
        func_params.update(kwargs)
        return spint.interp1d( data[:,0]*u.GHz, data[:,1], **func_params)
    
    @staticmethod
    def _init_value(value, **kwargs):
        '''
        Integrated importing of files
        '''
        if isinstance(value, str):
             return Optic.file_to_function(value, **kwargs)
        return value
    
    def _check_values(self, values, check):
        '''
        check defines if we need to enforce values between
        0 and 1. This is true for coefficients and false for
        temperatures
        '''
        if not check:
            return values
        if np.isscalar(values):
            if values < 0:
                return 0
            if values > 1:
                return 1
            return values
        values[values < 0] = 0
        values[values > 1] = 1
        return values
    
    def _get_value(self, value, freqs, check=True):
        '''
        Decides if the value is just a scalar or needs to be
        evalutated as a function of frequency
        
        check defines if we need to enforce values between
        0 and 1. This is true for coefficients and false for
        temperatures.
        '''
        if np.isscalar(value):
            return self._check_values(value, check)*np.ones( (len(freqs),))
        return self._check_values(value(freqs), check)
    
    def get_transmission(self, freqs):
        return 1-self.get_reflection(freqs)-self.get_absorption(freqs) - self.get_spill(freqs)
    
    def get_reflection(self, freqs):
        return self._get_value(self.reflection, freqs)
    
    def get_absorption(self, freqs):
        return self._get_value(self.absorption, freqs)
        
    def get_spill(self, freqs):
        return self._get_value(self.spill, freqs)
    
    def get_temperature(self, freqs):
        return self._get_value(self.temperature, freqs, check=False)
    def get_spill_temperature(self,freqs):
        return self._get_value(self.spill, freqs, check=False)
    
    def get_S_matrix(self, freqs):
        return np.array([[self.get_reflection(freqs), self.get_transmission(freqs)],
                         [self.get_transmission(freqs), self.get_reflection(freqs)]])
    
    def get_transmitted_power(self, freqs, power):
        return self.get_transmission(freqs)*power
    
    def get_reflected_power(self, freqs, power):
        return self.get_reflection(freqs)*power
    
    def get_absorbed_power(self, freqs, power):
        return self.get_absorption(freqs)*power
        
    def get_emitted_power(self, freqs):
        bb = Blackbody(self.get_temperature(freqs)*u.Kelvin)(freqs )
        return self.get_absorption(freqs) * bb * utils.through_put(freqs)
    
    def get_spill_power(self, freqs):
        bb = Blackbody(self.get_spill_temperature(freqs)*u.Kelvin)(freqs )
        return self.get_spill(freqs) * bb * utils.through_put(freqs)
    
    '''
    Functions for summarizing how a detector band goes through an optic
    '''
    def band_average_transmission(self, freqs, detector):
        trans = self.get_transmission(freqs)*detector.bandpass(freqs)
        return  np.trapz(trans, freqs)/np.trapz(detector.bandpass(freqs), freqs)

    def band_average_reflection(self, freqs, detector):
        trans = self.get_reflection(freqs)*detector.bandpass(freqs)
        return  np.trapz(trans, freqs)/np.trapz(detector.bandpass(freqs), freqs)

    def band_average_absorption(self, freqs, detector):
        trans = self.get_absorption(freqs)*detector.bandpass(freqs)
        return  np.trapz(trans, freqs)/np.trapz(detector.bandpass(freqs), freqs)
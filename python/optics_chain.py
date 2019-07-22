import os
import yaml
import numpy as np
import astropy.units as u
import astropy.constants as const
import astropy.modeling.blackbody as blackbody

import optlo_calc.utils as utils
import optlo_calc.optics as optics

class Optics_Chain(object):
    def __init__( self, optics):
        self.optics_list = optics

        super(Optics_Chain, self).__init__()

    def __str__(self):
        s = ''
        for optic in self.optics_list:
            s += '{}\t{}\n'.format(optic.name, optic.temperature)
        return s
        
    def calc_power_transfer(self, freqs, power_in):
        '''
        First version, ignores reflection
        '''
        stage_power = np.zeros( (len(self.optics_list)+1,
                                 len(freqs)))*(u.Watt/u.Hz)
        stage_power[0] = power_in
        transmission_list = [optic.get_transmission(freqs) for optic in self.optics_list]
        source_list = [optic.get_emitted_power(freqs) + 
                       optic.get_spill_power(freqs) for optic in self.optics_list]

        for o, optic in enumerate(range(len(self.optics_list))):
            stage_power[o+1] = transmission_list[o]*stage_power[o] + source_list[o]

        return stage_power, source_list, transmission_list

    def get_detector_loading(self, temperature_in, freqs=None):
        '''
        Assumes diffraction limited system
        Assumes last element in optical chain is a detector
        
        The last power 
        '''
        if freqs is None:
            freqs = np.linspace(2, 1000, 10000)*u.GHz
        power_in = blackbody.blackbody_nu(freqs, temperature_in)*(utils.through_put(freqs) )
        stage_power, _ , _ = self.calc_power_transfer(freqs, power_in)

        return np.trapz( stage_power[-1], freqs).to(u.pW)
    
    @classmethod
    def from_config_file(cls, filename):
        base_dir = os.path.dirname(filename) + '/'

        info = yaml.load( open(filename, 'r'))
        elements = []

        ## Making all the optical elements
        for name in info.keys():
            piece = info[name]
            for key in piece:
                if type(piece[key]) == list:
                    piece[key] = piece[key][0]*getattr(u, piece[key][1])
                if type(piece[key]) == str:
                    if key in ['transmission', 'absorption', 'bandpass', 'reflection']:
                        piece[key] = base_dir + piece[key]
            vals = info[name].copy()
            typ = vals.pop('type')
            vals['name'] = name
            elements.append( getattr(optics, typ)(**vals) )

        return cls(elements)
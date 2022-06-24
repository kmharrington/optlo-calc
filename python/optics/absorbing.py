import numpy as np
import astropy.units as u
import astropy.constants as const

from .optic import Optic

## Should I make a constants file? Probs.
eta0 = np.sqrt(const.mu0/const.eps0 ).to(u.ohm)

class Slab_Filter(Optic):
    '''
    Slab filters are lumps of material that have a thickness, index, and loss
    tangent. Their AR coating is just specified by a number (or file?). The
    absoption and emission estimates for these filters are "mostly right" but
    are missing the Fabre Parot (sp?) effects seen in real filters.
    '''
    def __init__(self, name, ar_reflection, index, losstan, thickness, 
                 temperature, spill, spill_temperature):
        self.index = index
        self.losstan = losstan
        self.thickness = thickness
        self.reflection = ar_reflection
        
        super(Slab_Filter, self).__init__(name, self.absorption,
                                         self.reflection, temperature, spill, 
                                          spill_temperature)
    def absorption(self, freqs):
        x = np.exp( -2*np.pi*freqs/const.c * self.thickness * self.index * self.losstan)
        return 1 - x.real


class _ABCD_Filter(Optic):
    '''
    ABCD filters are the more correct implementation of slab filters.
    This is a base class, classes inheriting form this class need to 
    over-write calc_ABCD.
    
    Use ABCD_Filter for single pieces of material and Stacked_ABCD_Filter 
    for combined pieces (like AR coatings)
    '''
    def __init__(self, name, temperature, spill,
                spill_temperature):
        super(_ABCD_Filter, self).__init__(name,
                                          self.absorption,
                                          self.reflection,
                                          temperature, spill, spill_temperature)
    
    def calc_ABCD(self, freqs):
        pass
    
    def reflection(self, freqs):
        '''
        The math here does assume we are in a vacuum
        '''
        abcd = self.calc_ABCD( freqs)
        return ABCD_reflection(abcd, eta0)
    
    def transmission(self, freqs):
        '''
        The math here does assume we are in a vacuum
        '''
        abcd = self.calc_ABCD(freqs)
        return ABCD_transmission(abcd, eta0, eta0)
        
    def absorption(self, freqs):
        '''
        This assumes we are in a vacuum, and yes in the big picture there
        is some circularity to how we run get_transmission for this optics.
        There may be some floating point errors.
        '''
        abcd = self.calc_ABCD(freqs)
        return ABCD_absorption(abcd, eta0, eta0)


class ABCD_Filter(_ABCD_Filter):
    '''
    ABCD filters are the more correct implementation of slab filters.
    
    Why to we use losstan for Slabs and imaginary index for this one?
    Because it's the math I had already figured out
    '''
    def __init__(self, name, n_re, n_im, tau, 
                 temperature, spill, spill_temperature):
        self.n_re = n_re
        self.n_im = n_im
        self.tau = tau
        
        index = self.n_re + 1.0j*self.n_im
        
        gamma = lambda nu: (-2.0j*np.pi*nu*index/const.c).to(1/u.cm)
        eta = lambda nu: (-2.0j*np.pi*nu/const.c * eta0 / gamma(nu)).to(u.ohm)
        
        self.A = lambda nu: np.cosh(gamma(nu)*self.tau*u.radian)
        self.B = lambda nu: 1.0*eta(nu)*np.sinh(gamma(nu)*self.tau*u.radian)
        self.C = lambda nu: 1.0/eta(nu)*np.sinh(gamma(nu)*self.tau*u.radian)
        self.D = lambda nu: np.cosh(gamma(nu)*self.tau*u.radian)
        
        super(ABCD_Filter, self).__init__(name, temperature, spill,
                spill_temperature)
    
    def calc_ABCD(self, freqs):
        A = self.A(freqs)
        B = self.B(freqs)
        C = self.C(freqs)
        D = self.D(freqs)
        return (A,B,C,D)

class ABCD_Filter_Angle(ABCD_Filter):
    """
    Extending ABCD filters to deal with angles
    """
    def __init__(self, name, n_re, n_im, tau, theta_i, temperature, spill,
                spill_temperature, align='TE', n_im_slope=0):
        """
        angle of incidence assumed to be in radians. 
        """
        self.n_re = n_re
        self.n_im = n_im
        self.tau = tau
        self.cosi = np.cos(theta_i)
        
        self.index = lambda nu : (self.n_re + 1.0j*self.n_im*nu.value**n_im_slope)
        
        self.gamma = lambda nu: (-2.0j*np.pi*nu*self.cosi*self.index(nu)/const.c).to(1/u.cm)
        if align == 'TE':
            self.eta = lambda nu: (-2.0j*np.pi*nu*self.cosi/const.c * eta0 / self.gamma(nu) / self.cosi).to(u.ohm)
        elif align == 'TM':
            self.eta = lambda nu: (-2.0j*np.pi*nu*self.cosi/const.c * eta0 *self.cosi / self.gamma(nu)  ).to(u.ohm)
            
        self.A = lambda nu: np.cosh(self.gamma(nu)*self.tau*u.radian)
        self.B = lambda nu: 1.0*self.eta(nu)*np.sinh(self.gamma(nu)*self.tau*u.radian)
        self.C = lambda nu: 1.0/self.eta(nu)*np.sinh(self.gamma(nu)*self.tau*u.radian)
        self.D = lambda nu: np.cosh(self.gamma(nu)*self.tau*u.radian)
        
        super(ABCD_Filter, self).__init__(name, temperature, spill,
                spill_temperature)
        
class Stacked_ABCD_Filter(_ABCD_Filter):
    
    def __init__(self, name, stack, temperature, spill, spill_temperature):
        self.stack = stack
        super(Stacked_ABCD_Filter, self).__init__(name,
                                              temperature, spill, spill_temperature)
        
    def calc_ABCD(self, freqs):
        '''
        Done this way because python has a limit to 
        the number of times you can recursively call a function
        '''
        A = self.stack[0].A(freqs)
        B = self.stack[0].B(freqs)
        C = self.stack[0].C(freqs)
        D = self.stack[0].D(freqs)
        for filt in self.stack[1:]:
            newA = A*filt.A(freqs) + B*filt.C(freqs)
            newB = A*filt.B(freqs) + B*filt.D(freqs)
            newC = C*filt.A(freqs) + D*filt.C(freqs)
            newD = C*filt.B(freqs) + D*filt.D(freqs)
            
            A=newA; B=newB
            C=newC; D=newD
        return A, B, C, D

def ABCD_reflection( ABCD, etaI):
    nom = (ABCD[0] + ABCD[1]/etaI - etaI*ABCD[2]-ABCD[3])
    denom = (ABCD[0] + ABCD[1]/etaI + etaI*ABCD[2]+ABCD[3])
    r = nom/denom
    return  np.abs(r * np.conjugate(r))

def ABCD_transmission( ABCD, etaI, etaT):
    denom = (etaT*ABCD[0] + ABCD[1] + etaI*etaT*ABCD[2] + etaI*ABCD[3])
    t = (2*etaT)/denom
    return  np.abs(t * np.conjugate(t))*(eta0/etaT)

def ABCD_absorption(ABCD, etaI, etaT):
    T = ABCD_transmission(ABCD, etaI, etaT)
    R = ABCD_reflection(ABCD, etaI)
    return 1-R-T
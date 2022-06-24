import astropy.units as u

from .utils import power_law_index, to_dB
from .optics.absorbing import ABCD_Filter_IndexFunct

def index_filter_from_params(params, thickness, angle=0):
    n_real, nu0, a, b = params
    filt = ABCD_Filter_IndexFunct('thing', thickness, 
                               angle, 4, 0, 4)
    filt.index = lambda nu: power_law_index(params, nu)
    return filt
    
def calc_index_transmission( params, ff, angle=0, thickness=6.35*u.mm):
    """
    ff is frequency in GHz with no units attached
    """
    test_T = index_filter_from_params(params, thickness, angle)
    T = test_T.get_transmission(ff*u.GHz)
    return T

def calc_index_transmission_dB(params, ff, angle=0, thickness=6.35*u.mm):
    T = calc_index_transmission(params, ff, angle, thickness)
    return to_dB(T)

def calc_index_reflection(params, ff, angle=0, thickness=6.35*u.mm):
    """
    ff is frequency in GHz with no units attached
    """
    test_R = index_filter_from_params(params, thickness, angle)
    R = test_R.get_reflection(ff*u.GHz)
    return R
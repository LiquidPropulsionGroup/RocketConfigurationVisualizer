import math

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from scipy.optimize import brentq

area_arr = [[-0.5, 3.0],[0.5, 1.5],[0,1.0]]
mach_arr = []
gam = 1.1363
area_ratio = 1.0 
    
def bartz(d_throat, p_chamber, c_star, d, c_p, visc, t_gas, t_wall):
    """bartz equation calculator"""
    t_boundary = (t_gas + t_wall) / 2
    print(t_boundary)
    return (0.026 / math.pow(d_throat, 0.2) * math.pow((p_chamber / c_star), 0.8) * math.pow((d_throat / d),
                                                                                             1.8) * c_p * math.pow(
        visc, 0.2) * math.pow((t_gas / t_boundary), (0.8 - 0.2 * 0.6)))



print(bartz(2.5, 2.2e7, 2600, 2.5 , 2800, 3e-5, 3200, 1000))
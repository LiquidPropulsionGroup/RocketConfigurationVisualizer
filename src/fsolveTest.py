import math

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from scipy.optimize import brentq

area_arr = [[-0.5, 3.0],[0.5, 1.5],[0,1.0]]
mach_arr = []
gam = 1.1363
area_ratio = 1.0 
    
def solveMach(area_ratio, mach_guess1=0.5, mach_guess2=1.5):
  def machEqn(mach):
    # return mach * area_ratio + 10
    return ( 2 / (gam + 1) * ( 1 + (gam - 1)/2 * mach**2 ))**((gam+1)/(2*(gam-1))) - mach * area_ratio
  
  mach_arr.append(fsolve(machEqn, [mach_guess1, mach_guess2]))

for area in area_arr:
  solveMach(area[1])
print(mach_arr)

def machEqn2(mach):
    return ( 2 / (gam + 1) * ( 1 + (gam-1)/2 * mach**2 ))**((gam+1)/(2*(gam-1))) - mach * area_ratio

def brentSolveMach():
    for elem in area_arr:
        area_ratio = elem[1]
        mach_arr.append(brentq(machEqn2, 0, 1))
    
    for elem in area_arr:
        area_ratio = elem[1]
        mach_arr.append(brentq(machEqn2, 1, 10))

print(mach_arr)
import math

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from scipy.optimize import brentq

area_arr = [[-0.5, 3.0],[0.5, 1.5],[0,1.0]]
mach_arr = []
gam = 1.1363
area_ratio = 1.0 
    
def machEqn(p):
    (mach, area_ratio) = p
    return ( 2 / (gam + 1) * ( 1 + (gam)/2 * mach**2 ))**((gam+1)/(2*(gam-1))) - mach * area_ratio

def machLine(p):
    return 0

def combo(p):
    return (machEqn(p), machLine(p))

def solveMach():
    mach_arr.append(fsolve(combo, (area_arr[0][1], 0.5)))
    mach_arr.append(fsolve(combo, (area_arr[1][1], 1.5)))
    mach_arr.append(fsolve(combo, (area_arr[2][1], 1)))

solveMach()
print(mach_arr)

def machEqn2(mach):
    return ( 2 / (gam + 1) * ( 1 + (gam)/2 * mach**2 ))**((gam+1)/(2*(gam-1))) - mach * area_ratio

def brentSolveMach():
    for elem in area_arr:
        area_ratio = elem[1]
        mach_arr.append(brentq(machEqn2, 0, 1))
    
    for elem in area_arr:
        area_ratio = elem[1]
        mach_arr.append(brentq(machEqn2, 1, 10))

print(mach_arr)
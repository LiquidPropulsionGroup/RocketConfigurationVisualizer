import math
#import matplotlib as plt
import numpy as np
#import scipy as sp

# The purpose of this script is to calculate supply pressure effects in the plumbing system

#   Helium thermophysical properties
m = 4.0026 # g/mol
cp = 5188 # J/kg
gamma = 1.66
Rbar = 2077 J/kg*K

#   van der Waals correction factors
a = 3.46e-3 # Pa * m^6
b = 23.71e-6 # m^3/mol
def vanderWaal( P, a, n, V, b, R ): #solving for T
    

fuel_volFlow = 3.015427325 # gpm
lox_volFlow = 4.778712667 # gpm

He_init_press = 1800 #psig
He_vol = 2990 # in^3
He_init_temp = 303 # K

fuel_volEmpty = 360.04009 # in^3
lox_volEmpty = 360.04009 # in^3


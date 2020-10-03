import math
#import matplotlib as plt
import numpy as np
#import scipy as sp

# The purpose of this script is to calculate supply pressure effects in the plumbing system

#   Helium thermophysical properties
m = 4.0026 # g/mol
cp = 5188 # J/kg
gamma = 1.66
Rbar = 2077 #J/kg*K

#   van der Waals correction factors
a = 3.46e-3 # Pa * m^6
b = 23.71e-6 # m^3/mol
def vanderWaal( P, a, n, V, b, R ): #solving for T
    return (1/R)*(V/n-b)(P+a(n/V)**2)

fuel_volFlow = 0.0001902439024 # m^3/s
lox_volFlow = 0.0003014899211 # m^3/s

He_init_press = 1800/14.7*101300 # Pa
He_vol = 2990 * 0.000016387064 # m^3
He_init_temp = 303 # K

fuel_volEmpty = 360.04009 * 0.000016387064 # m^3
lox_volEmpty = 360.04009 * 0.000016387064 # m^3

target_pressure = 20*101300 # Pa
n_tanks = target_pressue * (fuel_volEmpty + lox_volEmpty) / ()

print('Ideal Gas Solution:')
print('Initial He volume: {} cubic inches.'.format(He_vol))
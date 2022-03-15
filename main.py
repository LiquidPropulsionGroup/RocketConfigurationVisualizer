from src.chemistry import Chemistry
from src.rocket import Rocket
import math
import numpy as np
import matplotlib.pyplot as plt
import scipy as sp

title = 'Stinger Engine Sizing'
chems = Chemistry.parse_initVeriables(r'cea_data\test_2-3of_15atm')
#set veriables
mdot = 1         #max thrust mass flow rate
Lstar = 1.1
Dcham = 3.75 * 0.0254 #in meters
conv_angle = math.pi / 4 # rad, 45deg
div_angle = math.pi / 12  # rad, 15deg
wall_temp = 850 # K
r1 = 1
r2 = 1
r3 = 0.4
step = 5e-3
nozzle_type = 'conical'
#stinger = Rocket(title, chems, mdot, Lstar, Dcham, conv_angle, div_angle, wall_temp, nozzle_type, r1, r2, r3, step)
#stinger.variablesDisplay()
#stinger.graphDisplay()

#test
title = 'Stinger Test Sizing'
chems = Chemistry.parse_initVeriables(r'cea_data\stinger_test_25bar')
mdot = 3
l_star = 1.1
cham_d = 5.5 * 0.0254 #in meters
conv_angle = math.pi / 4 # rad, 45deg
div_angle = math.pi / 12  # rad, 15deg
wall_temp = 1000 # K
r1 = 1.5
r2 = 1.5
r3 = .5
step = 5e-4
propellant_ratio = 2.3
nozzle_type = 'bell80'
#stinger = Rocket(title, chems, mdot, l_star, cham_d, conv_angle, div_angle, wall_temp, nozzle_type, r1, r2, r3, step)
#stinger.variablesDisplay()
#stinger.graphDisplay()

#methate engine test
title = 'Methane Engine Test'
chems = Chemistry.parse_initVeriables(r'cea_data\methane_test_50bar')
mdot = 5
l_star = 1.1
cham_d = 5 * 0.0254 #in meters
conv_angle = math.pi / 4 # rad, 45deg
div_angle = math.pi / 8  # rad, 15deg
wall_temp = 1000 # K
r1 = 1
r2 = 1
r3 = 0.4
step = 5e-4
propellant_ratio = 2.7
nozzle_type = 'bell80'
methane = Rocket(title, chems, mdot, l_star, cham_d, conv_angle, div_angle, wall_temp, nozzle_type, r1, r2, r3, step)
methane.variablesDisplay()
methane.graphDisplay()

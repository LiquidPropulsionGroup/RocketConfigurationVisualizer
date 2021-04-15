from src.chemistry import Chemistry
from src.rocket import Rocket
from src.flight import Flight
import math

import numpy as np
import matplotlib.pyplot as plt
import scipy as sp

#stinger rocket sizing
title = 'Stinger Engine Sizing'
chems = Chemistry.parse_initVeriables('test_2-3of_15atm')
#set veriables
mdot = 1
l_star = 1.1
cham_d = 3.75 * 0.0254 #in meters
conv_angle = math.pi / 4 # rad, 45deg
div_angle = math.pi / 12  # rad, 15deg
wall_temp = 850 # K
r1 = 0.0254
r2 = 0.0127
r3 = 0.01905
step = 5e-4
propellant_ratio = 2.3
#stinger = Rocket(title, chems, mdot, l_star, cham_d, conv_angle, div_angle, wall_temp, r1, r2, r3, step)
#stinger.variablesDisplay()
#stinger.graphDisplay()

#test
title = 'Stinger Test Sizing'
chems = Chemistry.parse_initVeriables('stinger_test_25bar')
mdot = 3
l_star = 1.1
cham_d = 5.5 * 0.0254 #in meters
conv_angle = math.pi / 4 # rad, 45deg
div_angle = math.pi / 12  # rad, 15deg
wall_temp = 1000 # K
r1 = 0.0254
r2 = 0.0127
r3 = 0.01905
step = 5e-4
propellant_ratio = 2.3
#stinger = Rocket(title, chems, mdot, l_star, cham_d, conv_angle, div_angle, wall_temp, r1, r2, r3, step)
#stinger.variablesDisplay()
#stinger.graphDisplay()

#methate engine test
title = 'Methane Engine Test'
chems = Chemistry.parse_initVeriables('methane_test_50bar')
mdot = 5
l_star = 1.1
cham_d = 5 * 0.0254 #in meters
conv_angle = math.pi / 4 # rad, 45deg
div_angle = math.pi / 8  # rad, 15deg
wall_temp = 1000 # K
r1 = 0.0254
r2 = 0.0127
r3 = 0.01905
step = 5e-4
propellant_ratio = 2.7
#methane = Rocket(title, chems, mdot, l_star, cham_d, conv_angle, div_angle, wall_temp, r1, r2, r3, step)
#methane.variablesDisplay()
#methane.graphDisplay()

#rocket trajectory test
mRocket = 200
thrust = 9000
mDot = 1
htarget = 13716
dragCd = 0.5
vehicleArea = 1.0
hInit = 0.0
my_flight = Flight(mRocket, thrust, mDot, htarget, dragCd, vehicleArea, hInit)
my_flight.printInfo()

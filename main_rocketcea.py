from src.chemistry import Chemistry
from src.chemistryCEA import ChemistryCEA
from src.rocket import Rocket
from src.ceaRocket import ceaRocket
from src.flight import Flight
import math
from src.engine import Engine
import numpy as np
import matplotlib.pyplot as plt
import scipy as sp

#stinger rocket sizing
title = 'Stinger Engine Sizing'
ox = 'LOX'        # full propellant selection is availible at https://rocketcea.readthedocs.io/en/latest/propellants.html
fuel = 'RP1'
pMaxCham = 15.2     #max thrust chamber pressure in bar
Mr = 2.3 # propellant mixture ratio
pAmbient = 1.01325
pMinExitRatio = 0.45
#set veriables
mdotMax = 1         #max thrust mass flow rate
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
#stinger = Engine(title, fuel, ox, nozzle_type, Mr, pMaxCham, mdotMax, pMinExitRatio, Lstar, Dcham, wall_temp, r1, r2, r3, conv_angle, div_angle = div_angle, contourStep = step)
#stinger.variablesDisplay()
#stinger.graphDisplay()

#test engine 
title = 'Test Engine Sizing'
ox = 'LOX'        # full propellant selection is availible at https://rocketcea.readthedocs.io/en/latest/propellants.html
fuel = 'RP1'
pMaxCham = 15     #max thrust chamber pressure in bar
Mr = 2.3 # propellant mixture ratio
pAmbient = 1.01325
pMinExitRatio = 0.45
#set veriables
mdotMax = 20         #max thrust mass flow rate
Lstar = 1.1
Dcham = 12 * 0.0254 #in meters
conv_angle = math.pi / 4 # rad, 45deg
div_angle = math.pi / 12  # rad, 15deg
wall_temp = 850 # K
fuel_delta_t = 578 # K
fuel_cp = 2010 # J/KgK
r1 = 1
r2 = 1
r3 = 0.4
step = 5e-3
nozzle_type = 'bell80'
test = Engine(title, fuel, ox, nozzle_type, Mr, pMaxCham, mdotMax, pMinExitRatio, Lstar, Dcham, wall_temp, r1, r2, r3, conv_angle, fuel_delta_t, fuel_cp, div_angle = div_angle, contourStep = step)
test.variablesDisplay()


title = 'Stinger Test Sizing'
oxidizer = 'LOX'
fuel = 'RP1'
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
#stinger.variablesDisplay()
#stinger.graphDisplay()

#methate engine test
title = 'Methane Engine Test'
oxidizer = 'LOX'
fuel = 'CH4'
chems = Chemistry.parse_initVeriables('methane_test_50bar')
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
#methane.variablesDisplay()
#methane.graphDisplay()



#rocket trajectory test
title = "test rocket trajectory"
mRocket = 10
thrust = 9000
mDot = 2
htarget = 13716
dragCd = 0.5
vehicleArea = 1.0
hInit = 0.0
#my_flight = Flight(title, mRocket, thrust, mDot, htarget, dragCd, vehicleArea, hInit)
#my_flight.printInfo()


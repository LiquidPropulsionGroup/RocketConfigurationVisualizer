from src.chemistry import Chemistry
from src_cea.chemistryCEA import ChemistryCEA
from src.rocket import Rocket
from src_cea.flight import Flight
import math
from src_cea.engine import Engine
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
wall_temp = 500 # K
fuel_delta_t = 300 # K
fuel_cp = 2010 # J/KgK
r1 = 1
r2 = 1
r3 = 0.4
step = 5e-3
nozzle_type = 'conical'
stinger = Engine(title, fuel, ox, nozzle_type, Mr, pMaxCham, mdotMax, pMinExitRatio, Lstar, Dcham, wall_temp, r1, r2, r3, conv_angle, fuel_delta_t, fuel_cp, div_angle = div_angle, contourStep = step)
#stinger.variablesDisplay()
#stinger.graphDisplay()



#test engine 
title = 'Test Engine Sizing'
ox = 'LOX'        # full propellant selection is availible at https://rocketcea.readthedocs.io/en/latest/propellants.html
fuel =  'CH4'# 'Isopropanol70''Ethanol''RP1'
customFuel = [
    "Isopropanol70", 
    """fuel C3H8O-2propanol C 3 H 8 O 1    wt%=70.0
h,cal=-65133.0     t(k)=298.15   rho=0.786
fuel water H 2.0 O 1.0  wt%=30.0
h,cal=-68308.0  t(k)=298.15 rho,g/cc = 0.9998"""
]
print('{}'.format(customFuel[1]))
pMaxCham = 20     #max thrust chamber pressure in bar
Mr = 2.6 # propellant mixture ratio
pAmbient = 1.01325
pMinExitRatio = 0.45
#set veriables
mdotMax = 0.84        #max thrust mass flow rate
filmCoolingPercent = 0.0
Lstar = 0.61
Dcham = 2.5 * 0.0254 #in meters
conv_angle = math.pi / 4 # rad, 45deg
div_angle = math.pi / 12  # rad, 15deg
wall_temp = 650 # K
fuel_delta_t = 300 # K
fuel_cp = 2010 # J/KgK
r1 = 1
r2 = 1
r3 = 0.4
step = 5e-3
nozzle_type = 'bell80' #'conical'
test = Engine(title, fuel, ox, nozzle_type, Mr, pMaxCham, mdotMax, pMinExitRatio, Lstar, Dcham, wall_temp, r1, r2, r3, conv_angle, fuel_delta_t, fuel_cp, filmCoolingPercent = filmCoolingPercent, div_angle = div_angle, contourStep = step, customFuel = customFuel)
test.variablesDisplay()
test.graphDisplay()


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


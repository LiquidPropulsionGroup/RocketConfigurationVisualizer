import math
from src_cea.engine import Engine

title = 'propane Engine Sizing'
ox = 'AIR'
fuel = 'propane2'
customFuel = [
    "propane2", 
    """fuel C3H8 C 3 H 8     wt%=100.
h,cal=-30372.     t(k)=298.15   rho=0.5808"""
]
pMaxCham = 2.5     #max thrust chamber pressure in bar
Mr = 12 # propellant mixture ratio
pAmbient = 1.01325
pMinExitRatio = 0.45
#set veriables
mdotMax = 0.025        #max thrust mass flow rate
filmCoolingPercent = 0
Lstar = 1.1
Dcham = 1.5 * 0.0254 #in meters
conv_angle = math.pi / 4 # rad, 45deg
div_angle = math.pi / 12  # rad, 15deg
wall_temp = 650 # K
fuel_delta_t = 60 # K
fuel_cp = 1540 # J/KgK
r1 = 1
r2 = 1
r3 = 0.4
step = 5e-3
nozzle_type = 'conical' #'bell80'
test2 = Engine(title, fuel, ox, nozzle_type, Mr, pMaxCham, mdotMax, pMinExitRatio, Lstar, Dcham, wall_temp, r1, r2, r3, conv_angle, fuel_delta_t, fuel_cp, filmCoolingPercent = filmCoolingPercent, div_angle = div_angle, contourStep = step, customFuel = customFuel)
test2.variablesDisplay()
test2.graphDisplay()
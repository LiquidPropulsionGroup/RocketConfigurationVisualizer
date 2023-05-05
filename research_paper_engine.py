import math
from src_cea.engine import Engine
from src_cea.injectorCalculator import Injector
import numpy as np


title = 'research paper engine'
ox = 'LOX'        # full propellant selection is availible at https://rocketcea.readthedocs.io/en/latest/propellants.html
fuel =  'RP1'
pMaxCham = 20    #max thrust chamber pressure in bar
Mr = 2.3 # propellant mixture ratio
pAmbient = 1.01325
pMinExitRatio = 0.45
#set veriables
mdotMax = 1.37        #max thrust mass flow rate
filmCoolingPercent = .2
Lstar = 0.7
Dcham = 3.625 * 0.0254 #in meters
conv_angle = math.pi / 4 # rad, 45deg
div_angle = math.pi / 12  # rad, 15deg
wall_temp = 1200 # K
fuel_delta_t = 150 # K
fuel_cp = 1920 # J/KgK
r1 = 1
r2 = 1
r3 = 0.4
step = 5e-3
nozzle_type = 'bell80'
test1 = Engine(title, fuel, ox, nozzle_type, Mr, pMaxCham, mdotMax, pMinExitRatio, Lstar, Dcham, wall_temp, r1, r2, r3, conv_angle, fuel_delta_t, fuel_cp, filmCoolingPercent = filmCoolingPercent, div_angle = div_angle, contourStep = step, frozen = 1)
test1.variablesDisplay()
#test1.graphDisplay()
#test1.filewrite('researchEngineContourDataCSV')
#Mr = 1.5
#test4 = Engine(title, fuel, ox, nozzle_type, Mr, pMaxCham, mdotMax, pMinExitRatio, Lstar, Dcham, wall_temp, r1, r2, r3, conv_angle, fuel_delta_t, fuel_cp, filmCoolingPercent = filmCoolingPercent, div_angle = div_angle, contourStep = step, frozen = 0)


title = 'research paper igniter'
ox = 'O2'        # full propellant selection is availible at https://rocketcea.readthedocs.io/en/latest/propellants.html
fuel =  'CH4'
pMaxCham = 30    #max thrust chamber pressure in bar
Mr = 4 # propellant mixture ratio
pAmbient = 1.01325
pMinExitRatio = 0.45
#set veriables
mdotMax = 0.01        #max thrust mass flow rate
filmCoolingPercent = 0
Lstar = 0.4
Dcham = 0.6 * 0.0254 #in meters
conv_angle = math.pi / 4 # rad, 45deg
div_angle = math.pi / 12  # rad, 15deg
wall_temp = 1000 # K
fuel_delta_t = 107 # K
fuel_cp = 1920 # J/KgK
r1 = 1
r2 = 1
r3 = 0.4
step = 5e-3
nozzle_type = 'bell80'
#igniter = Engine(title, fuel, ox, nozzle_type, Mr, pMaxCham, mdotMax, pMinExitRatio, Lstar, Dcham, wall_temp, r1, r2, r3, conv_angle, fuel_delta_t, fuel_cp, filmCoolingPercent = filmCoolingPercent, div_angle = div_angle, contourStep = step, frozen = 1)
#igniter.variablesDisplay()
#igniter.graphDisplay()
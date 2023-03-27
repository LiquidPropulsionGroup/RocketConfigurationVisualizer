import math
from src_cea.engine import Engine


title = 'research paper engine'
ox = 'LOX'        # full propellant selection is availible at https://rocketcea.readthedocs.io/en/latest/propellants.html
fuel =  'RP1'
pMaxCham = 20    #max thrust chamber pressure in bar
Mr = 2.3 # propellant mixture ratio
pAmbient = 1.01325
pMinExitRatio = 0.45
#set veriables
mdotMax = 1.37        #max thrust mass flow rate
filmCoolingPercent = 0
Lstar = 1.1
Dcham = 3.625 * 0.0254 #in meters
conv_angle = math.pi / 4 # rad, 45deg
div_angle = math.pi / 12  # rad, 15deg
wall_temp = 1500 # K
fuel_delta_t = 107 # K
fuel_cp = 1920 # J/KgK
r1 = 1
r2 = 1
r3 = 0.4
step = 5e-3
nozzle_type = 'bell80'
test1 = Engine(title, fuel, ox, nozzle_type, Mr, pMaxCham, mdotMax, pMinExitRatio, Lstar, Dcham, wall_temp, r1, r2, r3, conv_angle, fuel_delta_t, fuel_cp, filmCoolingPercent = filmCoolingPercent, div_angle = div_angle, contourStep = step, frozen = 0)
test1.variablesDisplay()
test1.graphDisplay()
#test1.filewrite('researchEngineContourDataCSV')
#Mr = 1.5
#test4 = Engine(title, fuel, ox, nozzle_type, Mr, pMaxCham, mdotMax, pMinExitRatio, Lstar, Dcham, wall_temp, r1, r2, r3, conv_angle, fuel_delta_t, fuel_cp, filmCoolingPercent = filmCoolingPercent, div_angle = div_angle, contourStep = step, frozen = 0)

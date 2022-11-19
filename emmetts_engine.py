import math
from src_cea.engine import Engine


title = 'Isopropanol Engine Sizing'
ox = 'LOX'        # full propellant selection is availible at https://rocketcea.readthedocs.io/en/latest/propellants.html
fuel =  'RP1'# 'Ethanol''Isopropanol70'
customFuel = [
    "Isopropanol70", 
    """fuel C3H8O-2propanol C 3 H 8 O 1    wt%=70.0
h,cal=-65133.0     t(k)=298.15   rho=0.786
fuel water H 2.0 O 1.0  wt%=30.0
h,cal=-68308.0  t(k)=298.15 rho,g/cc = 0.9998"""
]
print('{}'.format(customFuel[1]))
pMaxCham = 25     #max thrust chamber pressure in bar
Mr = 2.3 # propellant mixture ratio
pAmbient = 1.01325
pMinExitRatio = 0.45
#set veriables
mdotMax = 1        #max thrust mass flow rate
filmCoolingPercent = 0.5
Lstar = 1.1
Dcham = 3.625 * 0.0254 #in meters
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
test1 = Engine(title, fuel, ox, nozzle_type, Mr, pMaxCham, mdotMax, pMinExitRatio, Lstar, Dcham, wall_temp, r1, r2, r3, conv_angle, fuel_delta_t, fuel_cp, filmCoolingPercent = filmCoolingPercent, div_angle = div_angle, contourStep = step, customFuel = customFuel, frozen = 1)
test1.variablesDisplay()
Mr = 1.5
test4 = Engine(title, fuel, ox, nozzle_type, Mr, pMaxCham, mdotMax, pMinExitRatio, Lstar, Dcham, wall_temp, r1, r2, r3, conv_angle, fuel_delta_t, fuel_cp, filmCoolingPercent = filmCoolingPercent, div_angle = div_angle, contourStep = step, customFuel = customFuel, frozen = 1)
test4.variablesDisplay()

title = 'propane Engine Sizing'
ox = 'AIR'
fuel = 'propane2'
customFuel2 = [
    "propane2", 
    """fuel C3H8 C 3 H 8     wt%=100.
h,cal=-30372.     t(k)=298.15   rho=0.5808"""
]
pMaxCham = 2.5     #max thrust chamber pressure in bar
Mr = 14 # propellant mixture ratio
pAmbient = 1.01325
pMinExitRatio = 0.45
#set veriables
mdotMax = 0.4        #max thrust mass flow rate
filmCoolingPercent = 0
Lstar = 1.1
Dcham = 1 * 0.0254 #in meters
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
# test2 = Engine(title, fuel, ox, nozzle_type, Mr, pMaxCham, mdotMax, pMinExitRatio, Lstar, Dcham, wall_temp, r1, r2, r3, conv_angle, fuel_delta_t, fuel_cp, filmCoolingPercent = filmCoolingPercent, div_angle = div_angle, contourStep = step, customFuel = customFuel)
# test2.variablesDisplay()
Mr = 1.5
# test3 = Engine(title, fuel, ox, nozzle_type, Mr, pMaxCham, mdotMax, pMinExitRatio, Lstar, Dcham, wall_temp, r1, r2, r3, conv_angle, fuel_delta_t, fuel_cp, filmCoolingPercent = filmCoolingPercent, div_angle = div_angle, contourStep = step, customFuel = customFuel)
# test3.variablesDisplay()
test1.graphDisplay()
# test2.graphDisplay()
# test3.graphDisplay()
test4.graphDisplay()
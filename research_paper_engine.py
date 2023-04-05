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


mdot = 0.3/7  #mass flow of one injection orifice
#pressures in pascals
p_f = 24*10**5
p_in = 23*10**5 # not currently being used
p_c = 20*10**5
alpha = 32.5*np.pi/180 #in radians
n = 4
l_n = 4.5   # l_in = 3-6
l_in = 1    # l_n = 0.5-2
l_s = 3     # l_s > 2
rho = 997   #in kg/m^3
nu = 0.6*10**(-6) #in m^2/s
swirlInjector = Injector(mdot, p_f, p_in, p_c, alpha, n, l_n, l_in, l_s, rho, nu)
swirlInjector.calculate()
print(swirlInjector)
swirlInjector.variablesDisplay()
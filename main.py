from src.chemistry import Chemistry
from src.rocket import Rocket
import math

import numpy as np
import matplotlib.pyplot as plt
import scipy as sp

chems = Chemistry.parse_initVeriables('test_real')

#print(chems)

#set veriables
mdot = 0.5
l_star = 1.1
cham_d = 0.08 #in meters
conv_angle = math.pi / 6 # rad, 30deg
div_angle = math.pi / 12  # rad, 15deg
rocket = Rocket(chems, mdot, l_star, cham_d, conv_angle, div_angle)

#fig1, ax = plt.figure()
#plt.plot(rocket.contour[0], rocket.contour[1])
#fig1.axis('equal')
#plt.plot(rocket.area_arr[0], rocket.area_arr[1])
#plt.plot(rocket.mach_arr[0], rocket.mach_arr[1])
#plt.plot(rocket.temp_arr[0], rocket.temp_arr[1])
#plt.plot(rocket.pressure_arr[0], rocket.pressure_arr[1])
#plt.plot(rocket.density_arr[0], rocket.density_arr[1])
plt.show()
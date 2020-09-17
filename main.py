from src.chemistry import Chemistry
from src.rocket import Rocket
from src.rocket import bartz
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
h_coeff_array = []
for i in range(len(rocket.temp_arr)):
    TEMP_VISC = 1.0142e-4 #need to add transport properties to config file and interpolate
    TEMP_tWall = 600 #K
    cStar = rocket.cham.p * rocket.thr.a / rocket.mdot
    temp_h = bartz(rocket.thr.d, rocket.cham.p, cStar, 2*rocket.contour[1][i], rocket.cham.cp, TEMP_VISC, rocket.temp_arr[1][i], TEMP_tWall)
    h_coeff_array.append([rocket.temp_arr[0][i], temp_h])

#fig1, ax = plt.figure()
#plt.plot(rocket.contour[0], rocket.contour[1])
#fig1.axis('equal')
plt.plot(rocket.area_arr[0], rocket.area_arr[1])
#plt.plot(rocket.mach_arr[0], rocket.mach_arr[1])
#plt.plot(rocket.temp_arr[0], rocket.temp_arr[1])
#plt.plot(rocket.pressure_arr[0], rocket.pressure_arr[1])
#plt.plot(rocket.density_arr[0], rocket.density_arr[1])
plt.plot(h_coeff_array[0], h_coeff_array[1])
plt.show()
from src.chemistry import Chemistry
from src.rocket import Rocket
import math

import numpy as np
import matplotlib.pyplot as plt
import scipy as sp

chems = Chemistry.parse_initVeriables('test_2-3of_15atm')

#print(chems)

#set veriables
mdot = 1
l_star = 1.1
cham_d = 3.75 * 0.0254 #in meters
#cham_d = 0.08
conv_angle = math.pi / 4 # rad, 45deg
div_angle = math.pi / 12  # rad, 15deg
rocket = Rocket(chems, mdot, l_star, cham_d, conv_angle, div_angle)

#rocket.contour = rocket.contour / 0.0254

print("Debug Outputs:\n")

print("Chamber Length: {0:.2f} in".format(rocket.chamber_length / 0.0254))
print("Chamber Diameter: {0:.2f} in".format(rocket.cham.d / 0.0254))
print("Exit Diameter: {0:.2f} in".format(rocket.exit.d / 0.0254))
print("Throat Diameter: {0:.2f} in".format(rocket.thr.d / 0.0254))
print("Total Length: {0:.2f} in".format((rocket.contourPoints[6][0]-rocket.contourPoints[0][0]) / 0.0254))
print("Volume: {0:.2f} cc".format(rocket.chamber_volume * 1000000))


#fig1, ax = plt.figure()
plt.plot(rocket.contour[0], rocket.contour[1])
plt.axis('equal')
#plt.plot(rocket.area_arr[0], rocket.area_arr[1])
#plt.plot(rocket.mach_arr[0], rocket.mach_arr[1])
#plt.plot(rocket.temp_arr[0], rocket.temp_arr[1])
#plt.plot(rocket.pressure_arr[0], rocket.pressure_arr[1])
#plt.plot(rocket.density_arr[0], rocket.density_arr[1])
plt.show()
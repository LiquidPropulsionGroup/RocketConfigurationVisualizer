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
rocket = Rocket(chems, mdot, l_star, cham_d, conv_angle, div_angle, 0.0254, 0.0127, 0.01905, 5e-4)
#rocket = Rocket(chems, mdot, l_star, cham_d, conv_angle, div_angle)

rocket.contour = rocket.contour / 0.0254

print("Debug Outputs:\n")

print("Chamber Length: {0:.2f} in".format(rocket.chamber_length / 0.0254))
print("Chamber Diameter: {0:.2f} in".format(rocket.cham.d / 0.0254))
print("Exit Diameter: {0:.2f} in".format(rocket.exit.d / 0.0254))
print("Throat Diameter: {0:.2f} in".format(rocket.thr.d / 0.0254))
print("Total Length: {0:.2f} in".format((rocket.contourPoints[6][0]-rocket.contourPoints[0][0]) / 0.0254))
print("Volume: {0:.2f} cc".format(rocket.chamber_volume * 1000000))

x = []
y = []
for elem in rocket.contourPoints:
    x.append(elem[0]/0.0254)
    y.append(elem[1]/0.0254)

'''
plt.figure('mdot = 1, l_star = 1.1, cham_d = 3.75, p = 15atm')
plt.plot(rocket.contour[0], rocket.contour[1], label="Contour")
plt.scatter(x, y, color="red", label="Points")
plt.axis('equal')
plt.xlabel('Length (in)')
plt.ylabel('Radius (in)')
plt.title('mdot = 1, l_star = 1.1, cham_d = 3.75, p = 15atm')
plt.legend() '''

# print(rocket.contourPoints)
#plt.plot(rocket.area_arr[0], rocket.area_arr[1])
fig, axs = plt.subplots(4,1, figsize=(8,14))

axs[0].set_title("Nozzle Geometry")
axs[0].plot(rocket.contour[0], rocket.contour[1])
axs[0].set(xlabel="Axial Position (in)", ylabel="Radial Distance (in)")
axs[0].axis('equal')
#axs[0].xlabel("Length (in)")
#axs[0].ylabel("Radius (in)")

#axs[1].set_title("Mach Number")
axs[1].plot(rocket.mach_arr[0], rocket.mach_arr[1])
axs[1].set(ylabel="Mach Number (M)")
#axs[1].xlabel("Length (m)")
#axs[1].ylabel("Mach Number (M)")

#axs[2].set_title("Gas Core Temperature")
axs[2].plot(rocket.temp_arr[0], rocket.temp_arr[1])
axs[2].set(ylabel="Gas Core Temperature (K)")
#axs[2].xlabel("Length (m)")
#axs[2].ylabel("Core Gas Temperature (K)")

#axs[3].set_title("Pressure")
axs[3].plot(rocket.pressure_arr[0], rocket.pressure_arr[1])
axs[3].set(ylabel="Pressure (atm)", xlabel="Axial Position (m)")
#axs[3].xlabel("Length (m)")
#axs[3].ylabel("Pressure (atm)")

#plt.plot(rocket.temp_arr[0], rocket.temp_arr[1])
#plt.plot(rocket.pressure_arr[0], rocket.pressure_arr[1])
#plt.plot(rocket.density_arr[0], rocket.density_arr[1])

fig2, axs2 = plt.subplots(3,1, figsize=(8,10.5))
axs2[0].set_title("Nozzle Geometry")
axs2[0].plot(rocket.contour[0], rocket.contour[1])
axs2[0].set(xlabel="Axial Position (in)", ylabel="Radial Distance (in)")
axs2[0].axis('equal')

axs2[1].plot(rocket.h_g_arr[0], rocket.h_g_arr[1])
axs2[1].set(ylabel="Coefficient of Heat Transfer (W/m^2*K)", xlabel="Axial Position (m)")

axs2[2].plot(rocket.heat_flux_arr[0], rocket.heat_flux_arr[1])
axs2[2].set(ylabel="Heat Flux Rate (W/m^2)", xlabel="Axial Position (m)")

rocket.contour = rocket.contour * 0.0254
rocket.filewrite("dataTest.txt")

plt.show()
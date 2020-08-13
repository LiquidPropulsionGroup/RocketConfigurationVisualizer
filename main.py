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
inj_d = 0.08
conv_angle = math.pi / 6 # rad, 30deg
div_angle = math.pi / 12  # rad, 15deg
rocket = Rocket(chems, mdot, l_star, inj_d, conv_angle, div_angle)

plt.plot(rocket.contour[0], rocket.contour[1])
plt.show()
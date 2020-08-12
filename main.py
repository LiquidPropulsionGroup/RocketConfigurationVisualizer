from src.chemistry import Chemistry
from src.rocket import Rocket
import math

chems = Chemistry.parse_initVeriables('test_real')

#print(chems)

#set veriables
mdot = 0.5
l_star = 1.1
inj_d = 0.08
conv_angle = math.pi / 6 # rad, 30deg
div_angle = math.pi / 12  # rad, 15deg
rocket = Rocket(chems, mdot, l_star, inj_d, conv_angle, div_angle)

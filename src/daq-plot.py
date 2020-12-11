import math
import numpy as np
import matplotlib.pyplot as plt

# Making a 3D meshgrid plot to show the relationship between bit depth, pressure range, and voltage

bit_depth = (4,16)
bit_space = np.linspace(bit_depth[0], bit_depth[1], len(range(bit_depth[1]-bit_depth[0]))+1)
pressure_range = 8700 #(300, 8700)
voltage_range = (5, 32)

bit_array = []
for i in bit_space:
    bit_array.append(2**i)

# debug
# print(bit_space)
# print(bit_array)

# pressure_array = np.linspace(pressure_range[0], pressure_range[1], 1000)
voltage_array = np.linspace(voltage_range[0], voltage_range[1], 1000)

BITS, VOLTS = np.meshgrid(bit_array, voltage_array)

result = VOLTS/BITS * pressure_range/VOLTS

# plt.figure()
# plt.contourf(BITS, VOLTS, result)
# plt.show() 


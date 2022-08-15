'''
fluid properties
h_fg        enthalpy of vaporization
T_sat       saturation temperature
roe_l       density in liquid phase
c_l         specific heat capacity liquid phase
'''
import os

class FluidProperties:
    h_fg = None
    T_sat = None
    roe_l = None
    c_l = None
    def __init__(self, fluidName):
        #print(os.getcwd())
        with open(f'src_cea\\fluidProperties\\{fluidName}.csv') as csvfile:
            lines = csvfile.readlines()
            for i in range(len(lines)):
                lines[i] = lines[i].split(',')
            for i in range(len(lines[0])):
                self.__setattr__(lines[0][i], float(lines[1][i]))

    def __repr__(self):
        return f"h_fg: {self.h_fg} " \
            f"T_sat: {self.T_sat} " \
            f"roe_l: {self.roe_l} " \
            f"c_l: {self.c_l} " \

if __name__ == '__main__':
    keroseneProps = FluidProperties('kerosene')
    print(keroseneProps)
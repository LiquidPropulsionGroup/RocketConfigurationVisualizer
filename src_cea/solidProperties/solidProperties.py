'''
solid properties
T
roe
Cp
mu
'''
import os
from re import T

class FluidProperties:
    T = None
    roe = None
    Cp = None
    mu = None
    def __init__(self, fluidName):
        #print(os.getcwd())
        with open(f'src_cea\\solidProperties\\{fluidName}.csv') as csvfile:
            lines = csvfile.readlines()
            for i in range(len(lines)):
                lines[i] = lines[i].split(',')
            for i in range(len(lines[0])):
                self.__setattr__(lines[0][i], float(lines[1][i]))

    def __repr__(self):
        return f"T: {self.T} " \
            f"roe: {self.roe} " \
            f"Cp: {self.Cp} " \
            f"mu: {self.mu} " \

if __name__ == '__main__':
    copperProps = FluidProperties('copper')
    print(copperProps)
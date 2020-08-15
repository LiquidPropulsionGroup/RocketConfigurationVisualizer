import math

import numpy as np
import matplotlib.pyplot as plt
import scipy as sp

# input values
# mdot #total mass flowrate into engine (kg/s)
# Lstar #characteristic length (m)
chamber_diameter = 0.08

def hoop_stress(internal_pressure, inside_diameter, wall_thickness):
    """
    hoop stress calculater
    takes: internal pressure, inside diameter of hoop, and wall thickness
    """
    hoopStress = internal_pressure * inside_diameter / 2 * wall_thickness
    return hoopStress

def bartz(d_throat, p_chamber, c_star, d, c_p, visc, t_gas, t_wall):
    """bartz equation calculator"""
    t_boundary = (t_gas + t_wall) / 2
    return (0.026 / math.pow(d_throat, 0.2) * math.pow((p_chamber / c_star), 0.8) * math.pow((d_throat / d),
                                                                                             1.8) * c_p * math.pow(
        visc, 0.2) * math.pow((t_gas / t_boundary), (0.8 - 0.2 * 0.6)))

# chamber diameter(0.08m), lambda curve, 15degree nozzle.

def nozzleAreaEquation(mdot, pres, temp, rbar, gam): 
    return (mdot / (pres * 101325)) * math.sqrt(temp * rbar / gam) * math.pow(
            (1 + ((gam - 1) / 2)), ((gam + 1) / (2 * (gam - 1))))

class Rocket:
    def __init__(self, chem, mdot, l_star, cham_d, conv_angle, div_angle):
        self.inj = chem[0] # injector
        self.cham = chem[1] # converging starts (end of chamber)
        self.thr = chem[2] # throat
        self.exit = chem[3] # exit
        self.mdot = mdot
        self.Lstar = l_star
        self.cham.d = cham_d #diameter of the chamber
        self.contourPoints = None
        self.contour = None
        self.area_arr = None
        self.conv_angle = conv_angle
        self.divergence_angle = div_angle
        #self.props = 

        # Specific impulse in seconds
        self.isp_s = self.exit.isp / 9.8

        self.thr.a = nozzleAreaEquation(self.mdot, self.thr.p, self.thr.t, self.thr.rbar, self.thr.gam)
        # p is in atm, conversion constant to Pa, might change to Pa later. area is in m^2

        # Nozzle Exit Area and diameters via Expansion Ratio and
        self.exit.a = self.thr.a * self.exit.ae
        self.thr.d = 2 * math.sqrt(self.thr.a / math.pi)
        self.exit.d = 2 * math.sqrt(self.exit.a / math.pi)

        # Thrust by fundamental rocket eq F = mdot * exhaust_velocity
        self.thrust = self.mdot * self.exit.isp  # + self.a_noz*(self.p-self.p_amb) not included as sea level expanded

        # Total Chamber Volume via Characteristic Length
        # NOTE: this volume does not take injector and contour volumes into consideration and is theirfor not completly accurate
        self.chamber_volume = self.Lstar * self.thr.a
        self.chamber_length = self.chamber_volume / (math.pi * (self.cham.d / 2) ** 2)

        # this generates the chamber and nozzle contour that is used for calculations
        self.my_contourPoints()
        self.genContour()
        self.areas()
        self.areaMach()

    # this generates the points that the gencontour function uses to make functions between
    # the points are refrenced from left to right in the graph
    def my_contourPoints(self, r1=0.05, r2=0.03, r3=0.025):
        o = [0,self.thr.d / 2]
        d = [-r2 * np.sin(self.conv_angle),o[1] + r2 * (1 - np.cos(self.conv_angle))]
        n = [r3 * np.sin(self.divergence_angle), o[1] + r3 * np.sin(1 - np.cos(self.divergence_angle))]
        e = [n[0] + ((self.exit.d / 2) - n[1]) * np.sin(math.pi/2 - self.divergence_angle)/np.sin(self.divergence_angle), self.exit.d / 2]
        a = [None, self.cham.d / 2]
        b = [None, self.cham.d / 2]
        c = [None,self.cham.d / 2 - (r1 * (1 - np.cos(self.conv_angle)))]
        c[0] = d[0] - (c[1] - d[1]) * (np.sin(math.pi/2 - self.conv_angle)/np.sin(self.conv_angle))
        b[0] = c[0] - (r1 * np.sin(self.conv_angle))
        a[0] = b[0] - self.chamber_length
        self.contourPoints = [a, b, c, d, o, n, e]

    def genContour(self, r1=0.05, convergence_angle=30, r2=0.03, r3=0.025, step=1e-4): 
        # This is the function that draws the discrete contour
        # these functions are refrenced from left to right of the graph
        functions = [
            lambda x: self.contourPoints[0][1],
            lambda x: np.sqrt(r1 ** 2 - (x - self.contourPoints[1][0]) ** 2) + self.contourPoints[1][1] - r1,
            lambda x: ((self.contourPoints[2][1] - self.contourPoints[3][1]) / (self.contourPoints[2][0] - self.contourPoints[3][0]))  * (x - self.contourPoints[2][0]) + self.contourPoints[2][1],
            lambda x: -np.sqrt(r2 ** 2 - (x - self.contourPoints[4][0]) ** 2) + self.contourPoints[4][1] + r2, 
            lambda x: -np.sqrt(r3 ** 2 - (x + self.contourPoints[4][0]) ** 2) + self.contourPoints[4][1] + r3, 
            lambda x: ((self.contourPoints[5][1] - self.contourPoints[6][1]) / (self.contourPoints[5][0] - self.contourPoints[6][0]))  * (x - self.contourPoints[5][0]) + self.contourPoints[5][1] 
        ]
            # 1: stait line
            # 2: circle
            # 3: stait line
            # 4: circle
            # 5: circle
            # 6: stait line

        num = np.int32(np.rint((self.contourPoints[6][0] - self.contourPoints[0][0]) / step))

        x = np.array([])
        y = np.array([])

        for i, fun in enumerate(functions):
            temp_x = np.linspace(self.contourPoints[i][0], self.contourPoints[i + 1][0], num)
            f = np.vectorize(fun)
            y = np.append(y, f(temp_x))
            x = np.append(x, temp_x)

        self.contour = np.array([x, y])


    def areas(self):
        
        self.area_arr = self.contour.copy()
        self.area_arr[1,:] = np.multiply(np.power(self.area_arr[1,:], 2),np.pi)
    
    def areaMach(self):
    # Unfortunately, the Area-Mach relation is not an onto function, needs to use lookup

        # Generate a linspace of machs, this will not exactly line up with the contour intervals
        machs = np.linspace(0.1, self.exit.mach, 10000)

        # Generate a set of area ratios
        areaRatios = np.power(machs, 2)
        areaRatios = np.multiply(areaRatios, (1+(self.cham.gam)/2))
        areaRatios = np.multiply(areaRatios, (2/(self.cham.gam+1)))
        areaRatios = np.power(areaRatios, ((self.cham.gam+1)/(2*(self.cham.gam-1))))
        areaRatios = np.divide(areaRatios, machs)
        
        # what follows next is unholy (i'm going to use a dictionary lookup)

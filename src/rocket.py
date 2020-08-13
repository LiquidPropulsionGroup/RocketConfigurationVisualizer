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
    def __init__(self, chem, mdot, l_star, inj_d, conv_angle, div_angle):
        self.inj = chem[0] # injector
        self.cham = chem[1] # converging starts (end of chamber)
        self.thr = chem[2] # throat
        self.exit = chem[3] # exit
        self.mdot = mdot
        self.Lstar = l_star
        self.inj_d = inj_d #diameter of injector
        self.contourPoints = 0
        self.contour = 0
        self.area_arr = 0

        #self.props = 

        # Specific impulse in seconds
        self.isp_s = self.exit.isp / 9.8

        # Throat Area Equation
        #self.a_thr = (self.mdot / (chem.p * 101325)) * math.sqrt(chem.t * self.rbar / chem.gam) * math.pow(
        #    (1 + ((chem.gam - 1) / 2)), ((chem.gam + 1) / (2 * (chem.gam - 1))))

        self.thr.a = nozzleAreaEquation(self.mdot, self.thr.p, self.thr.t, self.thr.rbar, self.thr.gam)
        # p is in atm, conversion constant to Pa, might change to Pa later. area is in m^2

        # Nozzle Exit Area and diameters via Expansion Ratio and
        self.exit.a = self.thr.a * self.exit.ae
        self.thr.d = 2 * math.sqrt(self.thr.a / math.pi)
        self.exit.d = 2 * math.sqrt(self.exit.a / math.pi)

        # Thrust by fundamental rocket eq F = mdot * exhaust_velocity

        self.thrust = self.mdot * self.exit.isp  # + self.a_noz*(self.p-self.p_amb) not included as sea level expanded

        # Total Chamber Volume via Characteristic Length
        self.chamber_volume = self.Lstar * self.thr.a

        # Temporary Chamber Volume (cylindrical approximation)
        self.chamber_length = self.chamber_volume / (math.pi * (self.inj_d / 2) ** 2)

        # Here goes the dimensions and convergence angle calculator, hardcoding 30.
        self.conv_angle = conv_angle
        self.divergence_angle = div_angle

        self.genContourPoints()
        self.genContour()

    def genContourPoints(self, r1=0.05, r2=0.03, r3=0.025):
        # Radius and origin at throat
        origin_thr = [0, self.thr.d / 2]

        # start with the chamber side (left in our coordinates)
        d = [r2 * np.sin(self.conv_angle), origin_thr[1] + r2 * (1 - np.cos(self.conv_angle))]

        c_y = self.inj_d / 2 - r1 * (1 - np.cos(self.conv_angle))
        c = [(c_y - d[1]) / np.tan(self.conv_angle) + d[0], c_y]

        b = [c[0] + r1 * np.sin(self.conv_angle), self.inj_d / 2]

        inj = [self.chamber_length, self.inj_d / 2]

        # Concactenate left points
        self.contourPoints = [inj, b, c, d, origin_thr]

        # Flip the array order (multiply all the x elements by -1, faster than sort)
        # self.contourPoints = self.contourPoints * [-1, 1]
        for i in range(len(self.contourPoints)):
            self.contourPoints[i][0] = self.contourPoints[i][0] * -1


        # now the right side
        n = [r3 * np.sin(self.divergence_angle), origin_thr[1] + r3 * (1 - np.cos(self.divergence_angle))]
        e = [n[0] + (self.exit.d / 2) / np.tan(self.divergence_angle), self.exit.d / 2]

        # Concactenate the right points
        self.contourPoints.append(n)
        self.contourPoints.append(e)


    def genContour(self, r1=0.05, convergence_angle=30, r2=0.03, r3=0.025, step=1e-4): 
        # This is the function that draws the discrete contour

        functions = [
            lambda x: self.contourPoints[0][1],  # straight line
            # lambda x: np.sqrt(abs(r1 ** 2 - (endpoints[1][1] - x) ** 2)) + endpoints[1][1] - r1,  # circle
            lambda x: np.sqrt(r1 ** 2 - (x - self.contourPoints[1][0]) ** 2) + self.contourPoints[1][1] - r1,
            lambda x: -np.pi * ((convergence_angle / 180) * (x - self.contourPoints[2][0])) + self.contourPoints[2][1],
            # straight line
            lambda x: -np.sqrt(r2 ** 2 - (x - self.contourPoints[4][0]) ** 2) + self.contourPoints[4][1] + r2,  # circle
            lambda x: -np.sqrt(r3 ** 2 - (x - self.contourPoints[4][0]) ** 2) + self.contourPoints[4][1] + r3,  # circle
            lambda x: ((self.contourPoints[5][1] - self.contourPoints[6][1]) / (
                    self.contourPoints[5][0] - self.contourPoints[6][0]))
                      * (x - self.contourPoints[5][0]) + self.contourPoints[5][1]  # straight line
        ]

        num = np.int32(np.rint((self.contourPoints[6][0] - self.contourPoints[0][0]) / step))

        x = np.array([])
        y = np.array([])

        for i, fun in enumerate(functions):
            temp_x = np.linspace(self.contourPoints[i][0], self.contourPoints[i + 1][0], num)
            f = np.vectorize(fun)
            y = np.append(y, f(temp_x))
            x = np.append(x, temp_x)

            
        self.contour = np.array([x, y])
    #testing new contour generator that will allow for easier utilization elsewhere
    def my_contour_function(self, r1=0.05, convergence_angle=30, r2=0.03, r3=0.025, step=1e-6):
        pass        
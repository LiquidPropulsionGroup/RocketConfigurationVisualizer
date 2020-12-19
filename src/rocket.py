import math

import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
import time
from scipy.optimize import fsolve

# input values
# mdot #total mass flowrate into engine (kg/s)
# Lstar #characteristic length (m)

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

def throatAreaEquation(mdot, pres, temp, rbar, gam): 
    return (mdot / (pres * 101325)) * math.sqrt(temp * rbar / gam) * math.pow(
            (1 + ((gam - 1) / 2)), ((gam + 1) / (2 * (gam - 1))))

class Rocket:
    def __init__(self, chem, mdot, l_star, cham_d, conv_angle, div_angle, r1=0.05, r2=0.03, r3=0.025, step=1e-4):
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
        self.mach_arr = []
        self.pressure_arr = None
        self.temp_arr = None
        self.density_arr = None
        self.h_g_arr = []
        self.heat_flux_arr = []
        self.hoopStress_arr = None
        self.conv_angle = conv_angle
        self.divergence_angle = div_angle
        self.gam = self.thr.gam
        
        # Specific impulse in seconds
        self.isp_s = self.exit.isp / 9.8

        self.thr.a = throatAreaEquation(self.mdot, self.cham.p, self.thr.t, self.thr.rbar, self.thr.gam)
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
        self.Cstar = self.cham.p*101325 * self.thr.a / self.mdot

        # hardcode temp fix
        # self.chamber_length = 6.444650495*0.0254

        # this generates the chamber and nozzle contour that is used for calculations
        self.my_contourPoints(r1, r2, r3)
        self.genContour(r1, r2, r3, step)
        
        # temporarily turn off all convergence dependent functions
        self.areas()
        self.solveMach()
        # self.areaMach()
        self.tempPressureDensity()
        self.calcBartz()
        self.calcHeatFlux()

    # this generates the points that the gencontour function uses to make functions between
    # the points are referenced from left to right in the graph
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
        locs = ['inj', 'b', 'c', 'd', 'o', 'n', 'e']
        xy = ['x', 'y']

        txtout = open('dims.txt','w')
        for i in range(len(self.contourPoints)):
            for j in range(2):
                txtout.write('"{0}_{1}"= {2}\n'.format(locs[i], xy[j], self.contourPoints[i][j]/0.0254))

    def genContour(self, r1=0.05, r2=0.03, r3=0.025, step=1e-4): 
        # This is the function that draws the discrete contour
        # these functions are referenced from left to right of the graph
        functions = [
            lambda x: self.contourPoints[0][1],
            lambda x: np.sqrt(r1 ** 2 - (x - self.contourPoints[1][0]) ** 2) + self.contourPoints[1][1] - r1,
            lambda x: ((self.contourPoints[2][1] - self.contourPoints[3][1]) / (self.contourPoints[2][0] - self.contourPoints[3][0]))  * (x - self.contourPoints[2][0]) + self.contourPoints[2][1],
            lambda x: -np.sqrt(r2 ** 2 - (x - self.contourPoints[4][0]) ** 2) + self.contourPoints[4][1] + r2, 
            lambda x: -np.sqrt(r3 ** 2 - (x + self.contourPoints[4][0]) ** 2) + self.contourPoints[4][1] + r3, 
            lambda x: ((self.contourPoints[5][1] - self.contourPoints[6][1]) / (self.contourPoints[5][0] - self.contourPoints[6][0]))  * (x - self.contourPoints[5][0]) + self.contourPoints[5][1] 
        ]
            # 1: straight line
            # 2: circle
            # 3: straight line
            # 4: circle
            # 5: circle
            # 6: straight line

        num = np.int32(np.rint((self.contourPoints[6][0] - self.contourPoints[0][0]) / step))

        x = np.array([])
        y = np.array([])

        for i, fun in enumerate(functions):
            temp_x = np.linspace(self.contourPoints[i][0], self.contourPoints[i + 1][0], num)
            f = np.vectorize(fun)
            y = np.append(y, f(temp_x))
            x = np.append(x, temp_x)

        self.contour = np.array([x, y])

    #array functions:
    #this finds area over the contour
    def areas(self):
        
        self.area_arr = self.contour.copy()
        self.area_arr[1,:] = (self.area_arr[1,:] ** 2) * np.pi# this is just pi*r^2 in array form

    def solveMach(self):
        def solveMatchForAreaRatio(area_ratio, mach_guess=0.7):
            def machEqn(mach):
                # return mach * area_ratio + 10
                return ( 2 / (self.gam + 1) * ( 1 + (self.gam - 1)/2 * mach**2 ))**((self.gam+1)/(2*(self.gam-1))) - mach * area_ratio
            
            return fsolve(machEqn, mach_guess)

        last = 0.7
        for [x, area] in self.area_arr.transpose():
            if x > 0:
                last = last + 1
            [mach] = solveMatchForAreaRatio(area/self.thr.a, last)
            self.mach_arr.append([x, mach])
            last = mach
        self.mach_arr = np.array(self.mach_arr).transpose()

    #################################################

    def machAreaEquation(self, tempM, area):   #gamma used is for the chamber gamma, it is not changed throughout the chamber. fix later
        gam = self.cham.gam
        myreturn = (1/tempM)**2 * (2/(gam+1)*(1+((gam-1)/2)*tempM**2))**((gam+1)/(gam-1)) - (area/self.thr.a)**2
        return myreturn

    def binarySearchConvergence(self, area, regimeSwitch):
        myMach = 1
        testMach = 0
        step = 0
        testVal = self.machAreaEquation(myMach,area)
        if (regimeSwitch): #supersonic
            referenceMach = self.exit.mach+0.1 #this number is so the last iteration does not get stuck
            step = abs(myMach - referenceMach)/2
            testMach = referenceMach - step
            testVal = self.machAreaEquation(testMach, area)
        else: #subsonic
            referenceMach = self.cham.mach
            step = abs(myMach - referenceMach)/2
            testMach = referenceMach + step
            testVal = self.machAreaEquation(testMach, area)

        while (testVal >= 0.00001 or testVal <= -0.00001):
            step /= 2
            if(regimeSwitch):
                if (testVal < 0):
                    testMach += step
                else:
                    testMach -= step
            else:
                if (testVal < 0):
                    testMach -= step
                else:
                    testMach += step
            testVal = self.machAreaEquation(testMach, area)
        return testMach
    
    def areaMach(self):
    # uses convergence solver to arrive at Machs
        regimeSwitch = False
        tempMach = self.cham.mach
        self.mach_arr = self.area_arr.copy()
        count = 0 #:kekw:
        for area in self.area_arr[1,:]:
            tempMach = self.binarySearchConvergence(area, regimeSwitch)
            self.mach_arr[1,count] = tempMach
            if (self.area_arr[1,count] == self.thr.a):
                regimeSwitch = True                
            count += 1

    def temp_eq(self, mach):
        gam = self.thr.gam
        # t_stag = self.cham.t * (1 + ((gam-1)/2 * self.cham.mach**2))
        # Note: ok technically, yes, the stagnation temperature needs to account for
        # gas velocity, but in our assumptions, t_0 assumed == t_cham as given by CEA

        t_stag = self.cham.t
        myreturn = t_stag * (1 + ((gam-1)/2 * mach**2))**(-1)
        return myreturn

    def pressure_eq(self, mach):
        gam = self.thr.gam
        p_stag = self.cham.p * (1 + ((gam-1)/2 * self.cham.mach**2))**(gam/(gam-1))
        myreturn = p_stag * (1 + ((gam-1)/2 * mach**2))**(-gam/(gam-1))
        return myreturn

    def density_eq(self, mach):#need to find chamber density
        self.cham.rho = 1 #THIS IS TEMPORARY... PLEASE DELETE LATER
        gam = self.cham.gam
        d_stag = self.cham.rho * (1 + ((gam-1)/2 * self.cham.mach**2))**(1/(gam-1))
        myreturn = d_stag * (1 + ((gam-1)/2 * mach**2))**(-1/(gam-1))
        return myreturn

    def tempPressureDensity(self):
        self.pressure_arr = self.mach_arr.copy()
        self.temp_arr = self.mach_arr.copy()
        #self.density_arr = self.mach_arr.copy()
        count = 0
        for mach in self.mach_arr[1,:]:
            self.temp_arr[1,count] = self.temp_eq(mach)
            self.pressure_arr[1,count] = self.pressure_eq(mach)
            #self.density_arr[1,count] = self.density_eq(mach)
            count += 1

        
    def hoopStress(self): #work in progress
        #hoopStress = internal_pressure * inside_diameter / 2 * wall_thickness
        self.hoopStress_arr = self.contour.copy()
        self.hoopStress_arr = np.multiply(self.hoopStress_arr,)#needs pressure at every point
    
    def filewrite(self, filename):
        output = open(filename, "w")
        output.write("X\tY\tMACH\tTEMP\tPressure\th_g\tFLUX\n")
        for i in range(len(self.contour[1,:])):
            output.write("{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}\t{:.4f}\n".format(self.contour[0,i], self.contour[1,i], self.mach_arr[1,i], self.temp_arr[1,i], self.pressure_arr[1,i], self.h_g_arr[1,i], self.heat_flux_arr[1,i]))
        output.close()

    def calcBartz(self):
        self.h_g_arr = self.contour.copy()
        for i in range(len(self.h_g_arr[0])):
            self.h_g_arr[1,i] = bartz(self.thr.d, self.cham.p*101325, self.Cstar, self.contour[1,i]*2, self.cham.cp, 1.0420e-4, self.temp_arr[1,i], 800)

    def calcHeatFlux(self):
        self.heat_flux_arr = self.h_g_arr.copy()
        for i in range(len(self.heat_flux_arr[0])):
            self.heat_flux_arr[1,i] = self.h_g_arr[1,i]*(self.temp_arr[1,i]-800)
import math

import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
import time
from scipy.optimize import fsolve
# input values
# mdot #total mass flowrate into engine (kg/s)
# Lstar #characteristic length (m)

def bartz(d_throat, p_chamber, c_star, d, c_p, visc, t_gas, t_wall):
    """bartz equation calculator"""
    t_boundary = (t_gas + t_wall) / 2
    return (0.026 / math.pow(d_throat, 0.2) * math.pow((p_chamber / c_star), 0.8) * math.pow((d_throat / d),
         1.8) * c_p * math.pow(visc, 0.2) * math.pow((t_gas / t_boundary), (0.8 - 0.2 * 0.6)))

# chamber diameter(0.08m), lambda curve, 15degree nozzle.

def throatAreaEquation(mdot, pres, temp, rbar, gam): 
    return (mdot / (pres)) * math.sqrt(temp * rbar / gam) * math.pow(
            (1 + ((gam - 1) / 2)), ((gam + 1) / (2 * (gam - 1))))

class Rocket:
    def __init__(self, title, chem, mdot, l_star, cham_d, conv_angle, div_angle, wall_temp, nozzle_type, r1=0.05, r2=0.03, r3=0.025, contourStep=1e-4):
        self.title = title
        self.nozzle_type = nozzle_type
        self.inj = chem[0] # injector
        self.cham = chem[1] # converging starts (end of chamber)
        self.thr = chem[2] # throat
        self.exit = chem[3] # exit
        self.mdot = mdot
        self.Lstar = l_star
        self.cham.d = cham_d #diameter of the chamber
        self.wall_temp = wall_temp
        self.contourPoints = None
        self.contour = None
        self.area_arr = None
        self.mach_arr = []
        self.pressure_arr = None
        self.temp_arr = None
        self.density_arr = None
        self.h_g_arr = []
        self.heat_flux_arr = []
        self.conv_angle = conv_angle
        self.divergence_angle = div_angle
        self.gam = self.thr.gam
        self.total_watts = 0
        self.r1 = r1
        self.r2 = r2
        self.r3 = r3
        self.contourStep = contourStep
        #unit conversions
        #bar to Pa
        self.inj.p = self.inj.p*100000
        self.cham.p = self.cham.p*100000
        self.thr.p = self.thr.p*100000
        self.exit.p = self.exit.p*100000

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
        self.Cstar = self.cham.p * self.thr.a / self.mdot

        # new methods for generating chamber and nozzle contour
        self.contourPoints, self.contour = self.nozzleGeneration()

        # this generates the chamber and nozzle contour that is used for calculations
        #self.my_contourPoints(r1, r2, r3)
        #self.genContour(r1, r2, r3, step)
        
        # temporarily turn off all convergence dependent functions
        self.area_arr = self.areas(self.contour)
        self.solveMach()
        # self.areaMach()
        self.tempPressureDensity()
        self.calcBartz()
        self.calcHeatFlux()
        self.totalWatts()
        

    # this generates the points that the gencontour function uses to make functions between
    # the points are referenced from left to right in the graph
    def TOPnozzleCoefficients(self, Rn, Re, Xn, thetaN, thetaE): # finds a,b,c for TOP parobolic function
        A = np.array([
            [2*Rn, 1, 0],
            [2*Re, 1, 0],
            [Rn**2, Rn, 1]])
        B = np.array([1/np.tan(thetaN), 1/np.tan(thetaE), Xn])
        X = np.linalg.solve(A, B)
        return X
    def nozzleGeneration(self): #this function creates contour points and functions for the chamber and nozzle geometry
        #this first section sets up points for the chamber and throat.  with the left being the chamber the right being the exit, the points go in order of a, b, c, d, o, n, e
        r1 = self.r1 * self.thr.d/2
        r2 = self.r2 * self.thr.d/2
        o = [0,self.thr.d / 2]
        d = [-r2 * np.sin(self.conv_angle),o[1] + r2 * (1 - np.cos(self.conv_angle))]
        c = [None,self.cham.d / 2 - (r1 * (1 - np.cos(self.conv_angle)))]
        c[0] = d[0] - (c[1] - d[1]) * (np.sin(math.pi/2 - self.conv_angle)/np.sin(self.conv_angle))
        b = [c[0] - (r1 * np.sin(self.conv_angle)), self.cham.d / 2]
        a = [b[0] - self.chamber_length, self.cham.d / 2]
        
        if self.nozzle_type == 'conical': # this sets the points and equations for a conical nozzle
            r3 = self.r3 * self.thr.d/2
            n = [r3 * np.sin(self.divergence_angle), o[1] + r3 * np.sin(1 - np.cos(self.divergence_angle))]
            e = [n[0] + ((self.exit.d / 2) - n[1]) * np.sin(math.pi/2 - self.divergence_angle)/np.sin(self.divergence_angle), self.exit.d / 2]
            #print('n variable: {}'.format(n))
            #print('e variable: {}'.format(e))
            contourPoints = [a, b, c, d, o, n, e] # temporary
            nozzleCurve = lambda x: ((contourPoints[5][1] - contourPoints[6][1]) / (contourPoints[5][0] - contourPoints[6][0]))  * (x - contourPoints[5][0]) + contourPoints[5][1]

        elif self.nozzle_type == 'bell80': # this sets the points and equations for an 80% bell nozzle
            r3 = self.r3 * self.thr.d/2
            thetaE = 7 *np.pi/180 # theta values found in table... hard coded temporarily
            thetaN = 33 *np.pi/180
            n = [r3 * np.sin(thetaN), o[1] + r3 * np.sin(1 - np.cos(thetaN))]
            e = [ 0.8 * (((math.sqrt(self.exit.ae)-1)*self.thr.d/2 / np.tan(15*np.pi/180))), self.exit.d / 2] # specificly for 80% bell nozzles
            contourPoints = [a, b, c, d, o, n, e] # temporary
            A = np.array([
                [2*n[1], 1, 0],
                [2*e[1], 1, 0],
                [n[1]**2, n[1], 1]])
            B = np.array([1/np.tan(thetaN), 1/np.tan(thetaE), n[0]])
            X = np.linalg.solve(A, B)
            aa = X[0]
            bb = X[1]
            cc = X[2]
            #print('exit r:{}'.format(self.exit.d/2))
            #print('thr r:{}'.format(self.thr.d/2))
            #print('exit area ratio:{}'.format(self.exit.ae))
            #print('n variable: {}'.format(n))
            #print('e variable: {}'.format(e))
            #print('A variable: {}'.format(A))
            #print('B variable: {}'.format(B))
            #print('X variable: {}'.format(X))
            #print('a: {}\nb: {}\nc: {}'.format(aa,bb,cc))
            nozzleCurve = lambda x: (-X[1] + (X[1]**2 - 4 * X[0] * (X[2]-x))**0.5) / (2*X[0]) # might need to change sign

        elif self.nozzle_type == 'dualbell': #work in progress, this sets the points and equations for a duel bell nozzle, in this there is an extra point 'm' between the n and e points
            r3 = self.r3 * self.thr.d/2
            thetaE1 = 7 *np.pi/180 # theta values found in table... hard coded temporarily
            thetaN1 = 33 *np.pi/180
            thetaE2 = 7 *np.pi/180
            thetaN2 = 33 *np.pi/180
            curvePercent1 = .7
            curvePercent2 = .8
            curve1ae = None #low altitude optimized area ratio
            curve2ae = self.exit.ae
            curve1rad = None #low altitude optimized radius
            curve2rad = self.exit.d / 2
            n = [r3 * np.sin(thetaN1), o[1] + r3 * np.sin(1 - np.cos(thetaN1))]
            m = [ curvePercent1 * (((math.sqrt(curve1ae)-1)*self.thr.d/2 / np.tan(15))), curve1rad]
            e = [ curvePercent2 * (((math.sqrt(curve2ae)-1)*self.thr.d/2 / np.tan(15))), curve2rad]
            contourPoints = [a, b, c, d, o, n, m, e] # temporary
            #X1 = TOPnozzleCoefficients(n[1], m[1], n[0], thetaN1, thetaE1)
            A = np.array([
                [2*n[1], 1, 0],
                [2*m[1], 1, 0],
                [n[1]**2, n[1], 1]])
            B = np.array([1/np.tan(thetaN1), 1/np.tan(thetaE1), n[0]])
            X1 = np.linalg.solve(A, B)
            nozzleCurve1 = lambda x: (-X1[1] + math.sqrt(X1[1]**2 - 4*X1[0]*(X1[2]-x)))/ 2*X1[0]
            A = np.array([
                [2*m[1], 1, 0],
                [2*e[1], 1, 0],
                [n[1]**2, n[1], 1]])
            B = np.array([1/np.tan(thetaN2), 1/np.tan(thetaE2), m[0]])
            X2 = np.linalg.solve(A, B)
            nozzleCurve2 = lambda x: (-X2[1] + math.sqrt(X2[1]**2 - 4*X2[0]*(X2[2]-x)))/ 2*X2[0]

        else: #this runs if the nozzle type input does not match any of the above nozzle types
            print("invalid nozzle type")
        

        functions = [
            lambda x: contourPoints[0][1],
            lambda x: np.sqrt(r1 ** 2 - (x - contourPoints[1][0]) ** 2) + contourPoints[1][1] - r1,
            lambda x: ((contourPoints[2][1] - contourPoints[3][1]) / (contourPoints[2][0] - contourPoints[3][0]))  * (x - contourPoints[2][0]) + contourPoints[2][1],
            lambda x: -np.sqrt(r2 ** 2 - (x - contourPoints[4][0]) ** 2) + contourPoints[4][1] + r2, 
            lambda x: -np.sqrt(r3 ** 2 - (x + contourPoints[4][0]) ** 2) + contourPoints[4][1] + r3, 
            nozzleCurve
        ]
        num = np.int32(np.rint((contourPoints[6][0] - contourPoints[0][0]) / self.contourStep))
        x = np.array([])
        y = np.array([])
        for i, fun in enumerate(functions):
            temp_x = np.linspace(contourPoints[i][0], contourPoints[i + 1][0], num)
            f = np.vectorize(fun)
            #if i == 5:
            #    print('x:{}\ty:{}'.format(temp_x,f(temp_x)))
            y = np.append(y, f(temp_x))
            x = np.append(x, temp_x)
        contour = np.array([x, y])

        #exports contour points for use in cad
        locs = ['inj', 'b', 'c', 'd', 'o', 'n', 'e']
        xy = ['x', 'y']
        txtout = open('dims.txt','w')
        for i in range(len(contourPoints)):
            for j in range(2):
                txtout.write('"{0}_{1}"= {2}\n'.format(locs[i], xy[j], contourPoints[i][j]/0.0254))

        return contourPoints, contour
    '''
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
    '''
    #array functions:
    #this finds area over the contour
    def areas(self, contour):
        area_arr = contour.copy()
        area_arr[1,:] = (area_arr[1,:] ** 2) * np.pi# this is just pi*r^2 in array form
        return area_arr

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

    def temp_eq(self, mach):#NOTE: stagnation values need improvment
        gam = self.thr.gam
        #t_stag = self.cham.t * (1 + ((gam-1)/2 * self.cham.mach**2))
        # Note: ok technically, yes, the stagnation temperature needs to account for
        # gas velocity, but in our assumptions, t_0 assumed == t_cham as given by CEA

        t_stag = self.inj.t
        myreturn = t_stag * (1 + ((gam-1)/2 * mach**2))**(-1)
        return myreturn

    def pressure_eq(self, mach):
        gam = self.thr.gam
        #p_stag = self.cham.p * (1 + ((gam-1)/2 * self.cham.mach**2))**(gam/(gam-1))
        p_stag = self.inj.p
        myreturn = p_stag * (1 + ((gam-1)/2 * mach**2))**(-gam/(gam-1))
        return myreturn

    def density_eq(self, mach):#need to find chamber density
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
    
    def filewrite(self, filename):
        output = open(filename, "w")
        offset = self.contour[0,0]
        for i in range(len(self.contour[0])):
            self.contour[0,i] += -offset
        output.write("X,Y,MACH,TEMP,Pressure,h_g,FLUX\n")
        for i in range(len(self.contour[1,:])):
            output.write("{:.4f},{:.4f},{:.4f},{:.4f},{:.4f},{:.4f},{:.4f}\n".format(self.contour[0,i], self.contour[1,i], self.mach_arr[1,i], self.temp_arr[1,i], self.pressure_arr[1,i], self.h_g_arr[1,i], self.heat_flux_arr[1,i]))
        output.close()

    def calcBartz(self):
        self.h_g_arr = self.contour.copy()
        for i in range(len(self.h_g_arr[0])):
            self.h_g_arr[1,i] = bartz(self.thr.d, self.cham.p, self.Cstar, self.contour[1,i]*2, self.cham.cp*1000, 0.85452e-4, self.temp_arr[1,i], self.wall_temp)

    def calcHeatFlux(self):
        self.heat_flux_arr = self.h_g_arr.copy()
        for i in range(len(self.heat_flux_arr[0])):
            self.heat_flux_arr[1,i] = self.h_g_arr[1,i]*(self.temp_arr[1,i]-self.wall_temp)

    def totalWatts(self):
        for i in range(len(self.heat_flux_arr[0])-1):
            self.total_watts = self.total_watts + abs(self.area_arr[0,i+1] - self.area_arr[0,i]) * self.heat_flux_arr[1,i]
    ##################################################################
    def variablesDisplay(self):
        print("{}{}:{}".format('\033[33m', self.title, '\033[0m'))
        print("Chamber Length: {0:.2f} in".format(self.chamber_length / 0.0254))
        print("Chamber Diameter: {0:.2f} in".format(self.cham.d / 0.0254))
        print("Exit Diameter: {0:.2f} in".format(self.exit.d / 0.0254))
        print("Throat Diameter: {0:.2f} in".format(self.thr.d / 0.0254))
        print("Total Length: {0:.2f} in".format((self.contourPoints[6][0]-self.contourPoints[0][0]) / 0.0254))
        print("Volume: {0:.2f} cc".format(self.chamber_volume * 1000000))
        print("Thrust: {0:.2f} N".format(self.thrust))
        print("Chamber heat flux constant: {0:.2f} W/m^2K".format(self.h_g_arr[1,1]))
        print("Chamber heat flux W/m^2: {0:.2f} W/m^2".format(self.heat_flux_arr[1,1]))
        print("Total Watts: {0:.2f} W".format(self.total_watts))
        print("Mass Flow Rate: {0:.2f} mdot".format(self.mdot))
        print("chamber pressure: {0:.2f} bar".format(self.pressure_arr[1,1]/100000))
        print()

    def graphDisplay(self, pressure_units = 'bar', distance_units = 'in'):
        #temperature units?
        if(pressure_units == 'bar'):
            Pcon = 100000 #bar
        elif(pressure_units == 'atm'):
            Pcon = 101325
        elif(pressure_units == 'psi'):
            Pcon = 6894.76
        else:
            Pcon = 1

        if(distance_units == 'in'):
            Dcon = 39.3701
        elif(distance_units == 'cm'):
            Dcon = 100
        elif(distance_units == 'mm'):
            Dcon = 1000
        else:
            Dcon = 1


        self.contour = self.contour * Dcon

        fig, axs = plt.subplots(2,1, figsize=(8,10.5))

        axs[0].set_title("Nozzle Geometry")
        axs[0].plot(self.contour[0], self.contour[1], label="self Contour")
        axs[0].set(xlabel="Axial Position (in)", ylabel="Radial Distance (in)")
        axs[0].axis('equal')

        secaxs = axs[0].twinx()
        secaxs.plot(self.mach_arr[0]*Dcon, self.mach_arr[1], label="Mach Number", color="green")
        secaxs.set(ylabel="Mach Number (M)")
        axs[0].legend(loc=(0,1))
        secaxs.legend(loc=(0.75,1))

        axs[1].plot(self.temp_arr[0], self.temp_arr[1], color="orange", label="Temperature")
        axs[1].set(ylabel="Gas Core Temperature (K)")

        secaxs1 = axs[1].twinx()
        secaxs1.plot(self.pressure_arr[0], self.pressure_arr[1]*Pcon, color="purple", label="Pressure")
        secaxs1.set(ylabel="Pressure (bar)", xlabel="Axial Position (m)")
        axs[1].legend(loc=(0,1))
        secaxs1.legend(loc=(0.8,1))

        #---------------------------------------------------------
        fig2, axs2 = plt.subplots(2,1, figsize=(8,10.5))
        axs2[0].set_title("Nozzle Geometry")
        axs2[0].plot(self.contour[0], self.contour[1])
        axs2[0].set(xlabel="Axial Position (in)", ylabel="Radial Distance (in)")
        axs2[0].axis('equal')

        axs2[1].plot(self.h_g_arr[0], self.h_g_arr[1], label="h", color="blue")
        axs2[1].set(ylabel="Coefficient of Heat Transfer (W/m^2*K)", xlabel="Axial Position (m)")

        sax = axs2[1].twinx()
        sax.plot(self.heat_flux_arr[0], self.heat_flux_arr[1], label="flux", color="g")
        sax.set(ylabel="Heat Flux Rate (W/m^2)", xlabel="Axial Position (m)")
        axs2[1].legend(loc=(0,1))
        sax.legend(loc=(0.8,1))

        self.contour = self.contour / Dcon
        self.filewrite("dataTest.txt")

        plt.show()
import math
import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
import time
from scipy.optimize import fsolve
from rocketcea.cea_obj import CEA_Obj, add_new_fuel
from .chemistryCEA import ChemistryCEA
from .thrustLevel import ThrustLevel
from .fluidProperties.fluidProperties import FluidProperties

class Engine:
    def __init__(self, title, fuel, ox, nozzle_type, Mr, pMaxCham, mdotMax, pMinExitRatio, Lstar, Dcham, wall_temp, r1, r2, r3, conv_angle, fuel_delta_t, fuel_cp, filmCoolingPercent = 0, div_angle = None, contourStep = 1e-4, customFuel = None):
        self.title = title
        self.fuel = FluidProperties(fuel)
        print(self.fuel)
        self.ox = FluidProperties(ox)
        if customFuel != None:
            add_new_fuel( customFuel[0], customFuel[1] )
        self.cea = CEA_Obj( oxName= ox, fuelName= fuel)
        self.nozzle_type = nozzle_type
        self.Mr = Mr
        self.filmCoolingPercent = filmCoolingPercent
        self.pMaxCham = pMaxCham
        self.mdotMax = mdotMax
        self.pMinExitRatio = pMinExitRatio
        self.Lstar = Lstar
        self.Dcham = Dcham #diameter of the chamber
        self.wall_temp = wall_temp
        self.fuel_delta_t = fuel_delta_t
        self.fuel_cp = fuel_cp
        self.r1 = r1
        self.r2 = r2
        self.r3 = r3
        self.conv_angle = conv_angle
        self.div_angle = div_angle
        self.contourStep = contourStep
        self.contourPoints = None
        self.contour = None
        self.area_arr = None
        self.max, self.min = self.variableThrustOptimizer()
        self.chamber_volume = self.Lstar * self.max.thr.a
        self.chamber_length = self.chamber_volume / (math.pi * (self.Dcham / 2) ** 2)
        self.contourPoints, self.contour = self.nozzleGeneration()
        self.area_arr = self.areas()
        self.max.heatCalcs(self.area_arr, self.contour, self.wall_temp, self.fuel_delta_t, self.fuel, self.Mr, self.filmCoolingPercent)
        self.min.heatCalcs(self.area_arr, self.contour, self.wall_temp, self.fuel_delta_t, self.fuel, self.Mr, self.filmCoolingPercent)
        
    def bartzHeatCalcs(self):
        self.max.heatCalcs(self.area_arr, self.contour, self.wall_temp, self.fuel_delta_t, self.fuel, self.Mr, self.filmCoolingPercent)
        self.min.heatCalcs(self.area_arr, self.contour, self.wall_temp, self.fuel_delta_t, self.fuel, self.Mr, self.filmCoolingPercent)

    def filmCoolingHeatCalcs(self): 
        self.max.heatCalcsFilmCooling()     #add arguments
        self.max.heatCalcsFilmCooling()

    def throttleLevelCalculator(self, pMin, ae):
        #print('pressure ratio:{}'.format(pMin))
        pMin *= 1.013
        cpGuess = self.pMaxCham
        cpStep = cpGuess
        mdotGuess = self.mdotMax
        mdotStep = mdotGuess
        pGuess = 0.0
        while abs(pGuess-pMin) > 0.0001:
            #print('cpGuess:{}\ncpStep:{}\nmdotGuess:{}\nmdotStep:{}\npGuess:{}\nminExitPressure:{}'.format(cpGuess, cpStep, mdotGuess, mdotStep, pGuess, pMin))
            if pGuess != 0:
                cpStep /= 2
                mdotStep /= 2
                if pGuess < pMin:
                    cpGuess += cpStep
                    mdotGuess += mdotStep
                else:
                    cpGuess -= cpStep
                    mdotGuess -= mdotStep
            #print('chamber pressure:{}\nmr:{}\nmdot:{}\narea array:{}\nae:{}'.format(cpGuess, self.Mr, mdotGuess, self.area_arr, ae))
            nozmin = ThrustLevel(self.cea, cpGuess, self.Mr, mdotGuess, self.area_arr, ae = ae)
            pGuess = nozmin.exit.p
            #print('pressure guess:{}'.format(pGuess))
        print('min nozzle exit pressure in pascals: {}'.format(nozmin.exit.p))
        return nozmin

    def variableThrustOptimizer(self):

        if self.nozzle_type == 'bell80' or 'conical':
            optimalP = 0.9

            nozmax = ThrustLevel(self.cea, self.pMaxCham, self.Mr, self.mdotMax, self.area_arr, pAmbient = optimalP)
            nozmin = self.throttleLevelCalculator(self.pMinExitRatio, nozmax.exit.aeat)
            return nozmax, nozmin

        elif self.nozzle_type == 'duelbell':
            optimalP1 = 0.9
            optimalp2 = 0.4
            noz1max = ThrustLevel(self.cea, self.pMaxCham, self.Mr, self.mdotMax, self.area_arr, pAmbient = optimalP1)
            self.throttleLevelCalculator(self.pMinExitRatio, noz1max.exit.aeat)#fix
            #print('min thrust throat area: {}'.format(noz1min.thr.a))
            print('max thrust throat area: {}'.format(noz1max.thr.a))
            noz2max = ThrustLevel(self.cea, self.pMaxCham, self.Mr, self.mdotMax, self.area_arr, pAmbient = optimalP2)
            self.throttleLevelCalculator(self.pMinExitRatio, noz2max.exit.aeat)
            #print('min thrust throat area: {}'.format(noz1min.thr.a))
            print('max thrust throat area: {}'.format(noz1max.thr.a))

    def nozzleGeneration(self):#, r1 = self.r1, r2 = self.r2, r3 = self.r3, dThr = self.max.thr.d, dCham = self.cham.d, dExit = self.exit.d, lenCham = self.chamber_length, convAngle = self.conv_angle, divAngle = self.divergence_angle, aeExit = self.max.exit.aeat): #this function creates contour points and functions for the chamber and nozzle geometry
        #this first section sets up points for the chamber and throat.  with the left being the chamber the right being the exit, the points go in order of a, b, c, d, o, n, e
        r1 = self.r1 * self.max.thr.d/2
        r2 = self.r2 * self.max.thr.d/2
        o = [0,self.max.thr.d / 2]
        d = [-r2 * np.sin(self.conv_angle),o[1] + r2 * (1 - np.cos(self.conv_angle))]
        c = [None, self.Dcham / 2 - (r1 * (1 - np.cos(self.conv_angle)))]
        c[0] = d[0] - (c[1] - d[1]) * (np.sin(math.pi/2 - self.conv_angle)/np.sin(self.conv_angle))
        b = [c[0] - (r1 * np.sin(self.conv_angle)), self.Dcham / 2]
        a = [b[0] - self.chamber_length, self.Dcham / 2]
        
        if self.nozzle_type == 'conical': # this sets the points and equations for a conical nozzle
            r3 = self.r3 * self.max.thr.d/2
            n = [r3 * np.sin(self.div_angle), o[1] + r3 * np.sin(1 - np.cos(self.div_angle))]
            e = [n[0] + ((self.max.exit.d / 2) - n[1]) * np.sin(math.pi/2 - self.div_angle)/np.sin(self.div_angle), self.max.exit.d / 2]
            contourPoints = [a, b, c, d, o, n, e] # temporary
            nozzleCurve = lambda x: ((contourPoints[5][1] - contourPoints[6][1]) / (contourPoints[5][0] - contourPoints[6][0]))  * (x - contourPoints[5][0]) + contourPoints[5][1]

        elif self.nozzle_type == 'bell80': # this sets the points and equations for an 80% bell nozzle
            r3 = self.r3 * self.max.thr.d/2
            thetaE = 7 *np.pi/180 # theta values found in table... hard coded temporarily
            thetaN = 33 *np.pi/180
            n = [r3 * np.sin(thetaN), o[1] + r3 * np.sin(1 - np.cos(thetaN))]
            e = [ 0.8 * (((math.sqrt(self.max.exit.aeat)-1)*self.max.thr.d/2 / np.tan(15*np.pi/180))), self.max.exit.d / 2] # specificly for 80% bell nozzles
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
            nozzleCurve = lambda x: (-X[1] + (X[1]**2 - 4 * X[0] * (X[2]-x))**0.5) / (2*X[0]) # might need to change sign

            '''
        elif self.nozzle_type == 'dualbell': #work in progress, this sets the points and equations for a duel bell nozzle, in this there is an extra point 'm' between the n and e points
            r3 = self.r3 * self.max.thr.d/2
            thetaE1 = 7 *np.pi/180 # theta values found in table... hard coded temporarily
            thetaN1 = 33 *np.pi/180
            thetaE2 = 7 *np.pi/180
            thetaN2 = 33 *np.pi/180
            curvePercent1 = .7
            curvePercent2 = .8
            optimalP1 = 0.9
            optimalp2 = 0.4
            noz1 = ceaRocket(self.fuel.name, self.ox.name, self.max.cham.p, self.Mr, self.max.mdot, title = 'smallAE', pAmbient = optimalP1)
            ae1 = noz1.exit.aeat
            print('throat1: ')
            noz2 = ceaRocket(self.fuel.name, self.ox.name, self.max.cham.p, self.Mr, self.max.mdot, title = 'largeAE', pAmbient = optimalP2)
            ae1 = noz1.exit.aeat
            curve1ae = ae1 #low altitude optimized area ratio
            curve2ae = ae2
            curve1rad = r1 #low altitude optimized radius
            curve2rad = r2
            n = [r3 * np.sin(thetaN1), o[1] + r3 * np.sin(1 - np.cos(thetaN1))]
            m = [ curvePercent1 * (((math.sqrt(curve1ae)-1)*self.max.thr.d/2 / np.tan(15))), curve1rad]
            e = [ curvePercent2 * (((math.sqrt(curve2ae)-1)*self.max.thr.d/2 / np.tan(15))), curve2rad]
            contourPoints = [a, b, c, d, o, n, m, e] # temporary
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
            '''
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

    def areas(self):
        area_arr = self.contour.copy()
        area_arr[1,:] = (area_arr[1,:] ** 2) * np.pi# this is just pi*r^2 in array form
        return area_arr


########################################################################################################################################################################
    def filewrite(self, filename): #needs update
        output = open(filename, "w")
        offset = self.contour[0,0]
        for i in range(len(self.contour[0])):
            self.contour[0,i] += -offset
        output.write("X,Y,MACH,TEMP,Pressure,h_g,FLUX\n")
        for i in range(len(self.contour[1,:])):
            output.write("{:.4f},{:.4f},{:.4f},{:.4f},{:.4f},{:.4f},{:.4f}\n".format(self.contour[0,i], self.contour[1,i], self.mach_arr[1,i], self.temp_arr[1,i], self.pressure_arr[1,i], self.h_g_arr[1,i], self.heat_flux_arr[1,i]))
        output.close()
###############################

    #printing veriables
    def variablesDisplay(self):
        print("{}{}:{}".format('\033[33m', self.title, '\033[0m'))
        print("Propellants:{}, {}".format(self.fuel.name, self.ox.name))
        print("Chamber Length: {0:.2f} in".format(self.chamber_length / 0.0254))
        print("Chamber Diameter: {0:.2f} in".format(self.Dcham / 0.0254))
        print("Exit Diameter: {0:.2f} in".format(self.max.exit.d / 0.0254))
        print("Throat Diameter: {0:.2f} in".format(self.max.thr.d / 0.0254))
        print("Total Length: {0:.2f} in".format((self.contourPoints[6][0]-self.contourPoints[0][0]) / 0.0254))
        print("Volume: {0:.2f} cc".format(self.chamber_volume * 1000000))
        print("{}{}:{}".format('\033[92m', 'Max Thrust', '\033[0m'))
        print("Max Thrust: {0:.2f} N".format(self.max.thrust))
        print("Chamber heat flux constant: {0:.2f} W/m^2K".format(self.max.h_g_arr[1,1]))
        print("Chamber heat flux W/m^2: {0:.2f} W/m^2".format(self.max.heat_flux_arr[1,1]))
        print("Total Watts: {0:.2f} W".format(self.max.total_watts))
        print("Max Fuel Heat Transfer: {0:.2f} W".format(self.max.max_fuel_heat))
        print("Mass Flow Rate: {0:.2f} mdot".format(self.max.mdot))
        print("Fuel Mass Flow Rate: {0:.2f} mdot".format(self.max.mdot/(self.max.mr+1)*(1+self.max.filmCoolingPercent)))
        print("Film Cooling Mass Flow Rate: {0:.2f} mdot".format(self.max.mdot/(self.max.mr+1)*(self.max.filmCoolingPercent)))
        print("chamber pressure: {0:.2f} bar".format(self.max.pressure_arr[1,1]))
        print("{}{}:{}".format('\033[92m', 'Min Thrust', '\033[0m'))
        print("Min Thrust: {0:.2f} N".format(self.min.thrust))
        print("Chamber heat flux constant: {0:.2f} W/m^2K".format(self.min.h_g_arr[1,1]))
        print("Chamber heat flux W/m^2: {0:.2f} W/m^2".format(self.min.heat_flux_arr[1,1]))
        print("Total Watts: {0:.2f} W".format(self.min.total_watts))
        print("Max Fuel Heat Transfer: {0:.2f} W".format(self.min.max_fuel_heat))
        print("Mass Flow Rate: {0:.2f} mdot".format(self.min.mdot))
        print("Fuel Mass Flow Rate: {0:.2f} mdot".format(self.min.mdot/(self.min.mr+1)*(1+self.min.filmCoolingPercent)))
        print("Film Cooling Mass Flow Rate: {0:.2f} mdot".format(self.min.mdot/(self.min.mr+1)*(self.min.filmCoolingPercent)))
        print("chamber pressure: {0:.2f} bar".format(self.min.pressure_arr[1,1]))
        print()
        
    def graphDisplay(self):
        self.max.graphDisplay()
        self.min.graphDisplay()

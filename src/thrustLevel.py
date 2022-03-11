import math
import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
import time
from scipy.optimize import fsolve
from rocketcea.cea_obj import CEA_Obj
from .chemistryCEA import ChemistryCEA

class ThrustLevel:
    def __init__(self, cea, pCham, Mr, mdot, area_arr, pAmbient = None, ae = None):
        if ae == None:
            #print('chamber pressure:{}\nambient pressure:{}'.format(pCham, pAmbient))
            chems = ChemistryCEA.create(cea, pCham, Mr, pAmbient = pAmbient)
        else:
            self.ambientP = pAmbient
            #print('chamber pressure:{}\nae:{}'.format(pCham, ae))
            chems = ChemistryCEA.create(cea, pCham, Mr, ae = ae)
        self.inj = None # injector
        self.cham = chems[0] # converging starts (end of chamber)
        self.thr = chems[1] # throat
        self.exit = chems[2] # exit
        self.mdot = mdot
        self.contour = None
        self.area_arr = area_arr
        self.mach_arr = []
        self.pressure_arr = None
        self.temp_arr = None
        self.density_arr = None
        self.h_g_arr = []
        self.heat_flux_arr = []
        self.total_watts = 0
        self.max_fuel_heat = 0
        #print('mdot:{}\nchamber pressure:{}\nthroat temperature:{}\nthroat rbar:{}\nthroat gamma:{}'.format(self.mdot, self.cham.p, self.thr.t, self.thr.rbar, self.thr.gam))
        self.thr.a = self.throatAreaEquation(self.mdot, self.cham.p, self.thr.t, self.thr.rbar, self.thr.gam)
        self.exit.a = self.thr.a * self.exit.aeat
        self.thr.d = 2 * math.sqrt(self.thr.a / math.pi)
        self.exit.d = 2 * math.sqrt(self.exit.a / math.pi)

        # Specific impulse in seconds
        #self.isp_s = self.exit.isp / 9.8
        self.Cstar = self.cham.p * self.thr.a / self.mdot
        self.calcThrust(1.01325)


    def heatCalcs(self, area_arr, contour, wall_temp, fuel_delta_t, fuel_cp, mr):
        self.area_arr = area_arr
        self.contour = contour
        self.wall_temp = wall_temp
        self.fuel_delta_t = fuel_delta_t
        self.fuel_cp = fuel_cp
        self.mr = mr
        #if self.area_arr != None:
        self.mach_arr = self.solveMach()
        self.temp_arr, self.pressure_arr = self.tempPressureDensity()
        self.h_g_arr = self.calcBartz()
        self.heat_flux_arr = self.calcHeatFlux()
        self.total_watts = self.totalWatts()
        self.max_fuel_heat = self.fuelWatts()

    def calcThrust(self, pAmbient):#make dependant on altitude input
        self.thrust = self.mdot * self.exit.mach * self.exit.son + (self.exit.p - pAmbient)*self.exit.a

    def bartz(self, d_throat, p_chamber, c_star, d, c_p, visc, t_gas, t_wall):
        """bartz equation calculator"""
        t_boundary = (t_gas + t_wall) / 2
        return (0.026 / math.pow(d_throat, 0.2) * math.pow((p_chamber / c_star), 0.8) * math.pow((d_throat / d),
            1.8) * c_p * math.pow(visc, 0.2) * math.pow((t_gas / t_boundary), (0.8 - 0.2 * 0.6)))

    def throatAreaEquation(self, mdot, press, temp, rbar, gam): 
        press *= 100000
        return (mdot / (press)) * math.sqrt(temp * rbar / gam) * math.pow(
                (1 + ((gam - 1) / 2)), ((gam + 1) / (2 * (gam - 1))))

    def solveMach(self):
        def solveMatchForAreaRatio(area_ratio, mach_guess=0.7):
            def machEqn(mach):
                # return mach * area_ratio + 10
                return ( 2 / (self.thr.gam + 1) * ( 1 + (self.thr.gam - 1)/2 * mach**2 ))**((self.thr.gam+1)/(2*(self.thr.gam-1))) - mach * area_ratio
            
            return fsolve(machEqn, mach_guess)
        mach_arr = []
        last = 0.7
        for [x, area] in self.area_arr.transpose():
            if x > 0:
                last = last + 1
            [mach] = solveMatchForAreaRatio(area/self.thr.a, last)
            mach_arr.append([x, mach])
            last = mach
        mach_arr = np.array(mach_arr).transpose()
        return mach_arr


    def temp_eq(self, mach):#NOTE: stagnation values need improvment
        gam = self.thr.gam
        t_stag = self.cham.t * (1 + ((gam-1)/2 * self.cham.mach**2))
        # Note: ok technically, yes, the stagnation temperature needs to account for
        # gas velocity, but in our assumptions, t_0 assumed == t_cham as given by CEA
        #t_stag = self.inj.t

        myreturn = t_stag * (1 + ((gam-1)/2 * mach**2))**(-1)
        return myreturn

    def pressure_eq(self, mach):
        gam = self.thr.gam# CHECK THIS!!!
        p_stag = self.cham.p * (1 + ((gam-1)/2 * self.cham.mach**2))**(gam/(gam-1))
        #p_stag = self.inj.p
        myreturn = p_stag * (1 + ((gam-1)/2 * mach**2))**(-gam/(gam-1))
        return myreturn

    def density_eq(self, mach):#need to find chamber density
        gam = self.cham.gam# CHECK THIS!!!
        d_stag = self.cham.rho * (1 + ((gam-1)/2 * self.cham.mach**2))**(1/(gam-1))
        myreturn = d_stag * (1 + ((gam-1)/2 * mach**2))**(-1/(gam-1))
        return myreturn

    def tempPressureDensity(self):
        pressure_arr = self.mach_arr.copy()
        temp_arr = self.mach_arr.copy()
        #self.density_arr = self.mach_arr.copy()
        count = 0
        #print(self.mach_arr)
        for mach in self.mach_arr[1,:]:
            temp_arr[1,count] = self.temp_eq(mach)
            pressure_arr[1,count] = self.pressure_eq(mach)
            #self.density_arr[1,count] = self.density_eq(mach)
            count += 1
        return temp_arr ,pressure_arr #density_arr

    def calcBartz(self):
        h_g_arr = self.contour.copy()
        for i in range(len(h_g_arr[0])):
            #print('thr diameter:{}\nchamber pressure:{}\nc star:{}\ncontour:{}\nchamber cp:{}\ntemperature array:{}\nwall temp{}'.format(self.thr.d, self.cham.p, self.Cstar, self.contour, self.cham.cp, self.temp_arr, self.wall_temp))
            h_g_arr[1,i] = self.bartz(self.thr.d, self.cham.p, self.Cstar, self.contour[1,i]*2, self.cham.cp*1000, 0.85452e-4, self.temp_arr[1,i], self.wall_temp)
        return h_g_arr

    def calcHeatFlux(self):
        heat_flux_arr = self.h_g_arr.copy()
        for i in range(len(heat_flux_arr[0])):
            heat_flux_arr[1,i] = self.h_g_arr[1,i]*(self.temp_arr[1,i]-self.wall_temp)
        return heat_flux_arr

    def totalWatts(self):
        total_watts = 0
        for i in range(len(self.heat_flux_arr[0])-1):
            total_watts = total_watts + abs(self.area_arr[0,i+1] - self.area_arr[0,i]) * self.heat_flux_arr[1,i]
        return total_watts
    
    def fuelWatts(self):
        return (self.fuel_cp*self.fuel_delta_t*self.mdot/(self.mr+1))


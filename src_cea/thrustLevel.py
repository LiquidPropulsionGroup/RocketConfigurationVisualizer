import math
import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
import time
from scipy.optimize import fsolve
from rocketcea.cea_obj import CEA_Obj
#from .chemistryCEA import ChemistryCEA
from .runCEA import RunCEA

class ThrustLevel:
    def __init__(self, cea, pCham, Mr, mdot, area_arr, pAmbient = None, ae = None, frozen = 1):
        if ae == None:
            #print('chamber pressure:{}\nambient pressure:{}'.format(pCham, pAmbient))
            chems = RunCEA.create(cea, pCham, Mr, pAmbient = pAmbient, frozen = frozen)
        else:
            self.ambientP = pAmbient
            #print('chamber pressure:{}\nae:{}'.format(pCham, ae))
            chems = RunCEA.create(cea, pCham, Mr, ae = ae, frozen = frozen)
        self.inj = None # injector
        if cea.fac_CR is not None: #finite area combustor
            self.inj = chems[0] #injector
            self.cham = chems[1] # converging starts (end of chamber)
            self.thr = chems[2] # throat
            self.exit = chems[3] # exit
        else:   # infinite area combustor
            self.inj = None #injector
            self.cham = chems[0] # converging starts (end of chamber)
            self.thr = chems[1] # throat
            self.exit = chems[2] # exit
        #debug prints
        #print(chems)
        #print('exit mach:{}'.format(self.exit.mach))
        self.mdot = mdot
        self.Mr = Mr
        self.filmCoolingPercent = 0
        self.cham.p = pCham
        self.pAmbient = pAmbient
        self.contour = None
        self.area_arr = area_arr
        self.mach_arr = []
        self.pressure_arr = None
        self.temp_arr = None
        self.density_arr = None
        self.fuel = None
        self.ox = None
        self.h_g_arr = []
        self.heat_flux_arr = []
        self.total_watts = 0
        self.max_fuel_heat = 0
        self.Cstar = self.cham.Cstar
        #print('mdot:{}\nchamber pressure:{}\nthroat temperature:{}\nthroat rbar:{}\nthroat gamma:{}'.format(self.mdot, self.cham.p, self.thr.t, self.thr.rbar, self.thr.gam))
        self.thr.a = self.throatAreaEquation(self.mdot, self.cham.p, self.thr.t, self.thr.rbar, self.thr.gam)
        self.exit.a = self.thr.a * self.exit.aeat
        self.thr.d = 2 * math.sqrt(self.thr.a / math.pi)
        self.exit.d = 2 * math.sqrt(self.exit.a / math.pi)
        #self.isp_s, self.nozzle_mode = self.get_isp()
        #self.Cstar = self.cham.p * self.thr.a / self.mdot
        self.calcThrust(1.01325)

        # Specific impulse in seconds---------------------------------------------------------------------------------------------------------------------------
        #self.isp_s = self.exit.isp / 9.8 #9.8 is gravitational constant
    '''
    def get_isp(self):
        from rocketcea.separated_Cf import sepNozzleCf

        IspVac = self.exit.ivac
        mw,gam = self.thr.m, self.thr.gam # throat gamma
        #PcOvPe = py_cea.rockt.app[ self.i_exit ]
        Pexit = self.exit.p*14.7/1.01325
        Cstar = float(self.Cstar)
        # ==================================
        
        CfOvCfvacAtEsep, CfOvCfvac, Cfsep, CfiVac, CfiAmbSimple, CfVac, epsSep, Psep = \
             sepNozzleCf(gam, self.exit.aeat, self.cham.p, self.pAmbient)
        
        if Pexit > Psep:
            # if not separated, use theoretical equation for back-pressure correction
            IspAmb = IspVac - Cstar*self.pAmbient*self.exit.aeat/self.cham.p/32.174 
            CfAmb = IspAmb * 32.174 / Cstar
        else:
            # if separated, use Kalt and Badal estimate of ambient thrust coefficient
            # NOTE: there are better, more modern methods available
            IspODEepsSep,CstarODE,Tc = \
                self.get_IvacCstrTc(Pc=self.cham.p, MR=self.Mr, eps=epsSep)
                
            CfvacAtEsep = CfVac * IspODEepsSep / IspVac
            
            CfAmb = CfvacAtEsep * CfOvCfvacAtEsep
            IspAmb = CfAmb * Cstar / 32.174 
        
        # figure out mode of nozzle operation
        if Pexit > Psep:
            if Pexit > self.pAmbient + 0.05:
                mode = 'UnderExpanded (Pe=%g)'%Pexit
            elif Pexit < self.pAmbient - 0.05:
                mode = 'OverExpanded (Pe=%g)'%Pexit
            else:
                mode = 'Pexit = %g'%Pexit
        else:
            mode = 'Separated (Psep=%g, epsSep=%g)'%(Psep,epsSep)

        return IspAmb/9.8, mode
    '''
    def heatCalcs(self, area_arr, contour, wall_temp, fuel_delta_t, fuel, mr, filmCoolingPercent):
        self.area_arr = area_arr
        self.contour = contour
        self.wall_temp = wall_temp
        self.fuel_delta_t = fuel_delta_t
        self.fuel = fuel
        self.mr = mr
        self.filmCoolingPercent = filmCoolingPercent
        st = time.time()
        self.mach_arr = self.solveMach()
        et = time.time()
        print(f'solveMach run time" {et-st}s')
        st = time.time()
        self.temp_arr, self.pressure_arr = self.tempPressureDensity()
        et = time.time()
        print(f'tempPressureDensity run time" {et-st}s')
        st = time.time()
        self.h_g_arr = self.calcBartz()
        et = time.time()
        print(f'calcBartz run time" {et-st}s')
        st = time.time()
        self.heat_flux_arr = self.calcHeatFlux()
        et = time.time()
        print(f'calcHeatFlux run time" {et-st}s')
        st = time.time()
        self.total_watts = self.totalWatts()
        et = time.time()
        print(f'totalWatts run time" {et-st}s')
        if fuel_delta_t != None and fuel != None:
            self.max_fuel_heat = self.fuelWatts()
        else:
            self.max_fuel_heat = None
        self.isp_s = self.thrust/(self.mdot*9.8)
        if filmCoolingPercent != 0:
            self.isp_adjusted = self.thrust/((self.mdot+self.mdot/(self.mr+1)*(self.filmCoolingPercent))*9.8) # isp = thrust/(mdot*g)
        else:
            self.isp_adjusted = None
        '''
    def heatCalcsFilmCooling(self, area_arr, contour, wall_temp, fuel_delta_t, fuel, mr, filmCoolingPercent):
        self.area_arr = area_arr
        self.contour = contour
        self.wall_temp = wall_temp
        self.fuel_delta_t = fuel_delta_t
        self.fuel = fuel
        self.mr = mr
        self.filmCoolingPercent = filmCoolingPercent
        self.mach_arr = self.solveMach()
        self.temp_arr, self.pressure_arr = self.tempPressureDensity()
        mdot_g = self.mdot
        mdot_c = self.mdot * self.filmCoolingPercent
        P_cc = self.cham.p
        D_cc = self.cham.d
        #properties of free stream gas:
        cp_g = self.cham.cp
        mu_g = self.cham.mu
        Pr_g = self.cham.Pr
        roe_g = self.cham.roe
        #liquid coolant properties:
        cp_l = self.fuel.cp
        h_fg = self.fuel.h_fg
        T_sat = self.fuel.T_sat
        roe_l = self.fuel.roe
        #other needed properties:
        f =     #fanning friction factor
        #calculated values:
        G_g = roe_g * u_g
        G_mean = G_g * (u_g - u_c)/u_g #coolant and chamber gas velocites needed
        Re_g = G_mean * D_cc/u_g
        #lamda ... solve for 1/(lamda)^.5 = 1.93*ln(Re*lamda^.5)-0.537
        St_0 = f/2/(1.2+11.8*(f/2)^.5*(Pr-1)*Pr^(-1/3))
        e_t = #obtain from test data?
        K_t = 1 + 4 * e_t
        h_0 = G_mean * cp_g * St_0 * K_t
        hstar_fg = h_fg + (T_sat - ???) * cp_l
        epselon_t = 
        Qdot_rad = 
        St = 
        h = 
        Qdot_tot = 
        mdot_v
        Γ_c = # blowing ratio or the coolant flow rate per circumference available for film cooling
        L_c = #film cooling length
        '''
        '''
        https://www.sciencedirect.com/science/article/pii/S0017931012003195?casa_token=U5ONT6cTbYsAAAAA:VQ47ewAtIETX8_y5dxiMSEA2abz2aqwbIjpYr0QxATVftTxfp45c582heKU_gtFZuT6mXjwACw
        values needed for film cooling calculation:
        mdot_g
        mdot_c
        P_cc
        D_cc
        properties of free stream gas:
        cp
        mu_g
        Pr_g
        roe_g
        liquid coolant properties:
        c_l
        h_fg
        T_sat
        roe_l
        calculated values:
        Re_g        calculated from G_mean
        lamda
        St_0
        e_t
        K_t
        h_0
        hstar_fg
        epselon_t
        Qdot_rad
        St
        h
        '''


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

    def solveMach(self): # IMPORTANT NOTE: area ratio might be wrong, might need to be inversed
        def solveMatchForAreaRatio(area_ratio, mach_guess=0.7):
            def machEqn(mach):
                # return mach * area_ratio + 10
                return ( 2 / (self.thr.gam + 1) * ( 1 + (self.thr.gam - 1)/2 * mach**2 ))**((self.thr.gam+1)/(2*(self.thr.gam-1))) - mach * area_ratio
            return fsolve(machEqn, mach_guess)
        
        mach_arr = []
        last = 0.7
        for [x, area] in self.area_arr.transpose():
            if x > 0 and last <= 1: 
                last = last + 1
            [mach] = solveMatchForAreaRatio(area/self.thr.a, last)
            mach_arr.append([x, mach])
            last = mach
        
        mach_arr = np.array(mach_arr).transpose()
        return mach_arr

    def temp_eq(self, mach):#NOTE: stagnation values need improvment
        #NOTE: check if gamma can be made into array along x axis
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
        gam = self.thr.gam# CHECK THIS!!!
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
            total_watts = total_watts + math.sqrt((self.contour[0,i] - self.contour[0,i+1])**2+(self.contour[1,i] - self.contour[1,i+1])**2)*2*math.pi* self.contour[1,i]* self.heat_flux_arr[1,i]
        return total_watts
    
    def fuelWatts(self):
        #print(f'fuel cp: {self.fuel.cp} \nfuel delta t: {self.fuel_delta_t}\nmdot: {self.mdot}\nmr: {self.mr}\nfilm cooling percent: {self.filmCoolingPercent}')
        return (self.fuel.cp * self.fuel_delta_t * (self.mdot/(self.mr+1) * (1+self.filmCoolingPercent)))

    def graphDisplay(self, pressure_units = 'bar', distance_units = 'in'):
        #temperature units?
        if(pressure_units == 'bar'):
            Pcon = 1 #bar
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

        fig, axs = plt.subplots(2,1, figsize=(8,10.5))

        axs[0].set_title("Nozzle Geometry")
        axs[0].plot(self.contour[0]*Dcon, self.contour[1]*Dcon, label="Radius")
        axs[0].set(xlabel="Axial Position, {}".format(distance_units), ylabel="Radial Distance, {}".format(distance_units))
        axs[0].axis('equal')

        secaxs = axs[0].twinx()
        secaxs.plot(self.mach_arr[0]*Dcon, self.mach_arr[1], label="Mach Number", color="green")
        secaxs.set(ylabel="Mach Number, M")
        axs[0].legend(loc=(0,1))
        secaxs.legend(loc=(0.70,1))

        axs[1].plot(self.temp_arr[0]*Dcon, self.temp_arr[1], color="orange", label="Temperature")
        axs[1].set(ylabel="Gas Core Temperature, K", xlabel="Axial Position, {}".format(distance_units))

        secaxs1 = axs[1].twinx()
        secaxs1.plot(self.pressure_arr[0]*Dcon, self.pressure_arr[1]*Pcon, color="purple", label="Pressure")
        secaxs1.set(ylabel="Pressure, {}".format(pressure_units), xlabel="Axial Position, {}".format(distance_units))
        axs[1].legend(loc=(0,1))
        secaxs1.legend(loc=(0.75,1))

        #---------------------------------------------------------
        fig2, axs2 = plt.subplots(2,1, figsize=(8,10.5))
        axs2[0].set_title("Nozzle Geometry")
        axs2[0].plot(self.contour[0]*Dcon, self.contour[1]*Dcon)
        axs2[0].set(xlabel="Axial Position, {}".format(distance_units), ylabel="Radial Distance, {}".format(distance_units))
        axs2[0].axis('equal')

        axs2[1].plot(self.h_g_arr[0]*Dcon, self.h_g_arr[1], label="Heat Transfer Coefficient", color="blue")
        axs2[1].set(ylabel="Coefficient of Heat Transfer, W/m^2*K", xlabel="Axial Position, {}".format(distance_units))

        sax = axs2[1].twinx()
        sax.plot(self.heat_flux_arr[0]*Dcon, self.heat_flux_arr[1], label="Heat Flux", color="g")
        sax.set(ylabel="Heat Flux Rate, W/m^2", xlabel="Axial Position, {}".format(distance_units))
        axs2[1].legend(loc=(0,1))
        sax.legend(loc=(0.75,1))
        
        plt.show()
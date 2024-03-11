import math
from engine import Engine
from rocketcea.cea_obj import CEA_Obj, add_new_fuel
from .fluidProperties.fluidProperties import FluidProperties

class EngineAnalyzer:
    def __init__(self, title, fuel, ox, nozzle_type, Mr, pMaxCham, mdotMax, Lstar, Dcham, wall_temp, r1, r2, r3, conv_angle, fuel_delta_t, pMinExitRatio = [], filmCoolingPercent = [], contourStep = 5e-3, customFuel = customFuel, frozen = [0], fac_CR = None, pAmbient = [1.01325], doContours = True):
        self.title = title
        self.fuel = FluidProperties(fuel) #CEA
        self.ox = FluidProperties(ox) #CEA
        test = Engine(title, fuel, ox, nozzle_type, Mr, pMaxCham, mdotMax, Lstar, Dcham, wall_temp, r1, r2, r3, conv_angle, fuel_delta_t, pMinExitRatio = pMinExitRatio, filmCoolingPercent = filmCoolingPercent, contourStep = step, customFuel = customFuel, frozen = 0, pAmbient = pAmbient, doContours = doContours)
        self.nozzle_type = nozzle_type #in engine
        self.Mr = Mr #CEA
        self.pMaxCham = pMaxCham #CEA
        self.mdotMax = mdotMax #in engine?
        self.Lstar = Lstar #in engine
        self.Dcham = Dcham #in engine
        self.wall_temp = wall_temp #in engine?
        self.r1 = r1
        self.r2 = r2
        self.r3 = r3
        self.conv_angle = conv_angle #in engine
        self.fuel_delta_t = fuel_delta_t 
        self.frozen = frozen #CEA
        self.pAmbient = pAmbient #CEA
        self.filmCoolingPercent = filmCoolingPercent #in engine
        self.pMinExitRatio = pMinExitRatio #??????
        self.fac_CR = fac_CR 
        self.step = contourStep #no
        self.doContours = doContours

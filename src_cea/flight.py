import numpy as np
import math
import scipy as sp
from scipy.optimize import fsolve
import time

def fGravity(mass, hight):
    G = 6.6743e-11 #gravity constant
    mEarth = 5.972e24 #mass of earth in kg
    radiusEarth = 6.371e6 #radius of earth in m
    f = G * mEarth * mass / (radiusEarth + hight)**2
    return f

def drag(Cd, density, area, velocity):
    drag = Cd * density * area * velocity**2 / 2
    return drag

def hightPrediction(h, v, g): #note, gravity is constant and therefor inaccurate
    return h + v**2 / (2 * g)

def density1(roeB, tempB, h, hB, lapse, g = 9.8):
    roe = roeB * (tempB/(tempB + lapse*(h-hB)))**(1+ g * 0.003484/lapse)
    return roe

def density2(roeB, tempB, h, hB, g = 9.8):
    roe = roeB * math.exp(-g * 0.003484 * (h-hB) / tempB)
    return roe

def densityAltitude(h):
    hB = [0, 11000, 20000, 32000, 47000, 51000, 71000]
    tempB = [288.15, 216.65, 216.65, 228.65, 270.65, 270.65, 214.65]
    roeB = [1.225, 0.36391, 0.08803, 0.01322, 0.00143, 0.00086, 0.000064]
    lapse = [-0.0065, 0.0, 0.001, 0.0028, 0.0, -0.0028, -0.002]
    if(h<=11000): # toposphere
        B = 0
    elif(h<=20000): # lower stratosphere
        B = 1
    elif(h<=32000): # upper stratosphere
        B = 2
    elif(h<=47000): # lower mesosphere
        B = 3
    elif(h<=51000): # upper mesosphere
        B = 4
    elif(h<=71000): # lower thermosphere
        B = 5
    elif(h<=100000): #upper thermosphere
        B = 6
    else:
        B = 7
    
    if(B == 7):
        return 0.0
    elif(lapse[B] == 0):
        return density2(roeB[B], tempB[B], h, hB[B])
    else:
        return density1(roeB[B], tempB[B], h, hB[B], lapse[B])
        
class Flight:
    def __init__(self, title, mRocket, thrust, mDot, htarget, dragCd, vehicleArea = 0, hInit = 0, stepSize = 1e-3):
        self.title = title
        self.mRocket = mRocket #mass of rocket and payload in kg
        self.thrust = thrust #Newtons
        self.mDot = mDot #mass flow in kg/s
        self.fireTime = None #time length of engine burn in seconds
        self.hTarget = htarget #target apogee hight in m
        self.vehicleArea = vehicleArea #surface area of rocket in m^2
        self.dragCd = dragCd #drag coeffecient
        self.hInit = hInit #initial altitude of rocket above sea level in m
        self.mPropellants = None #mass of propellants in kg
        self.stepSize = stepSize
        self.aMax = None
        self.vMax = None
        self.calculateFuelNeeded()

    def calculateFuelNeeded(self):
        guess = 120.0 #seconds
        guessStep = guess/2
        step = self.stepSize
        hPrediction = 0.0
        while(abs(self.hTarget - hPrediction) > 1): #hight within 1 meter
            if(hPrediction != 0.0):
                if(self.hTarget < hPrediction):
                    guess -= guessStep
                else:
                    guess += guessStep
                guessStep /= 2

            mFuel = guess * self.mDot
            m = self.mRocket + mFuel
            h = self.hInit
            v = 0.0
            a = 0.0
            aMax = 0.0
            vMax = 0.0
            for j in range(math.ceil(guess/step)):
                h1 = h + v*step
                m1 = m - self.mDot*step
                v1 = v + a * step
                a1 = (self.thrust - fGravity(m, h) - drag(self.dragCd, densityAltitude(h), self.vehicleArea, v))/m
                h = h1
                m = m1
                v = v1
                a = a1
                j += j #delete
                if v > vMax:
                    vMax = v
                if a > aMax:
                    aMax = a
            hPrediction = hightPrediction(h, v, (fGravity(m, h)/m))
            #print("Fuel Mass: {0:.2f} kg".format(guess * self.mDot))
            #print("acceleration: {0:.2f} ".format(a))
            #print("velocity: {0:.2f} ".format(v))
            #print("hight: {0:.2f} ".format(h))
            #print("mass: {0:.2f} ".format(m))
            #print("hight prediction: {0:.2f} ".format(hPrediction))
            #time.sleep(1)
        self.vMax = vMax
        self.aMax = aMax
        self.mPropellants = mFuel
        self.fireTime = guess
    
    def printInfo(self):
        print("{}{}:{}".format('\033[33m', self.title, '\033[0m'))
        print("mass of fuel: {0:.2f} kg".format(self.mPropellants))
        print("fire time: {0:.2f} s".format(self.fireTime))
        print("thrust: {0:.2f} N".format(self.thrust))
        print("mass flow rate: {0:.2f} kg/s".format(self.mDot))
        print("max elevation of flight: {0:.2f} kg".format(self.hTarget))
        print("max acceleration: {0:.2f} m/s^2".format(self.aMax))
        print("max velocity: {0:.2f} m/s".format(self.vMax))
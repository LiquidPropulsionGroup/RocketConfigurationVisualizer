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
'''
def densityAltitude(h, t = 288.14):
    if(h<=11000): #troposphere
        t -= 0.006545*h
    elif(h<=20000): #lower stratosphere
        t = 216
    else: #upper stratosphere
        t = 216 + 0.001*h
    p = 101290 * math.exp(0.034159*h)
    roe = 0.2869*p/t
    roe = 
    print("temp: {0:.2f} K".format(t))
    print("pressure: {0:.2f} Pa".format(p))
    return roe
'''
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

print("density: {0:.20f}".format(densityAltitude(200000)))

class Flight:
    def __init__(self, mRocket, thrust, mDot, htarget, dragCd, vehicleArea = 0, hInit = 0):
        self.mRocket = mRocket #mass of rocket and payload in kg
        self.thrust = thrust #Newtons
        self.mDot = mDot #mass flow in kg/s
        self.fireTime = None #time length of engine burn in seconds
        self.hTarget = htarget #target apogee hight in m
        self.vehicleArea = vehicleArea #surface area of rocket in m^2
        self.dragCd = dragCd #drag coeffecient
        self.hInit = hInit #initial altitude of rocket above sea level in m
        self.mPropellants = None #mass of propellants in kg
        print("starting calculate funtion")
        self.calculateFuelNeeded()

    def calculateFuelNeeded(self):
        guess = 120.0 #seconds
        guessStep = guess/2
        step = 1e-3
        hPrediction = 0.0
        while(abs(self.hTarget - hPrediction) > 1): #hight within 1 meter
            if(hPrediction != 0.0):
                if(self.hTarget < hPrediction):
                    guess -= guessStep
                else:
                    guess += guessStep
                guessStep /= 2
            print("guess: {0:.2f}".format(guess))
            print("guessStep: {0:.2f}".format(guessStep))
            mFuel = guess * self.mDot
            m = self.mRocket + mFuel
            h = self.hInit
            v = 0.0
            a = 0.0
            print("number of iterations: {0:.2f} ".format(math.ceil(guess/step)))
            for j in range(math.ceil(guess/step)):
                #print("acceleration: {0:.2f} ".format(a))
                #print("velocity: {0:.2f} ".format(v))
                #print("hight: {0:.2f} ".format(h))
                #print("mass: {0:.2f} ".format(m))
                #time.sleep(0.001)
                h1 = h + v*step
                #print(h)
                m1 = m - self.mDot*step
                v1 = v + a * step
                a1 = (self.thrust - fGravity(m, h) - drag(self.dragCd, densityAltitude(h), self.vehicleArea, v))/m
                h = h1
                m = m1
                v = v1
                a = a1
                j += j #delete
            print("Fuel Mass: {0:.2f} kg".format(guess * self.mDot))
            print("acceleration: {0:.2f} ".format(a))
            print("velocity: {0:.2f} ".format(v))
            print("hight: {0:.2f} ".format(h))
            print("mass: {0:.2f} ".format(m))
            hPrediction = hightPrediction(h, v, (fGravity(m, h)/m))
            print("hight prediction: {0:.2f} ".format(hPrediction))
            time.sleep(1)

        self.mPropellants = mFuel



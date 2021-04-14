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

def densityAltitude(h, t = 288.14):
    if(h<=11000): #troposphere
        t -= 0.006545*h
    elif(h<=20000): #lower stratosphere
        t = 216
    else: #upper stratosphere
        t = 216 + 0.001*h
    p = 101290 * math.exp(0.034159*h)
    roe = 0.003483*p/t
    return roe

print(densityAltitude(1000)/100000)

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
        step = 1e-4
        hPrediction = 0.0
        while(abs(self.hTarget - hPrediction) > 1): #hight within 1 meter
            guess /= 2
            m = self.mRocket + (guess * self.mDot)
            h = self.hInit
            v = 0.0
            a = 0.0
            print("number of iterations: {0:.2f} ".format(math.ceil(guess/step)))
            for j in range(math.ceil(guess/step)):
                print("acceleration: {0:.2f} ".format(a))
                print("velocity: {0:.2f} ".format(v))
                print("hight: {0:.2f} ".format(h))
                print("mass: {0:.2f} ".format(m))
                time.sleep(0.001)
                h1 = h + v*step
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
            time.sleep(0.1)

        self.mPropellants = guess * self.mDot



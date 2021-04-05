#from RocketCEA.rocketcea import *
from RocketCEA.rocketcea.cea_obj import CEA_Obj
myC = CEA_Obj(oxName = 'LOX', fuelName = 'LH2')
print(myC.get_Isp(Pc = 100.0, MR = 1.0, eps = 40.0))
from rocketcea.cea_obj import CEA_Obj

class ChemistryCEA:
    aeat = None # area of exit/area of throat
    rho = None #density
    h = None #enthalpy
    s = None #entropy
    gam = None #molar heat capacity
    g = None #gibb's free energy
    u = None #internal energy
    m = None #molar mass
    mw = None #molecular weight = 1/m
    p = None #pressure in bar... is converted to Pa in rocket
    son = None #speed of sound
    cp = None #specific heat capacity (constant pressure)
    t = None #temperature
    ae = None #exit area
    cf = None #function of nozzle
    ivac = None #specific impulse (exit velocity) if the rocket were in space
    mach = None #mach number
    pip = None #pressure ratio to chamber
    isp = None #the exhaust velocity of the gas (m/s)


    a = None #area
    d = None #diameter

    def initCalculations(self):

        self.rbar = 8.31446261815324 / self.m * 1000 #ADD TO MAIN

    @staticmethod
    def create(cea, pCham, Mr, ae = None, pAmbient = None):
        #print('cea:{}\npCham:{}\nMr:{}\nae:{}\npAmbient{}'.format(cea, pCham, Mr, ae, pAmbient))
        if ae == None:
            Pratio = pCham/pAmbient
            #print('Pratio:{}'.format(Pratio))
            string = cea.get_full_cea_output(Pc = pCham, MR = Mr, eps = ae, PcOvPe = Pratio, pc_units='bar', output='KJ', short_output=1)
        elif ae != None:
            #print('this chem')
            string = cea.get_full_cea_output(Pc = pCham, MR = Mr, eps = ae, PcOvPe = None, pc_units='bar', output='KJ', short_output=1)
        else:
            print('chem needs a pAmbient or ae value')
        lines = string.splitlines()
        my_vars = [
            [51, 'p'],
            [52, 't'],
            [54, 'h'],
            [59, 'm'],
            [62, 'cp'],
            [63, 'gam'],
            [64, 'son'],
            [65, 'mach'],
            [86, 'aeat'],
            #[87, 'cstar'],
            [89, 'ivac'],
            [90, 'isp']
        ]
        classVals = []
        for i in my_vars:
            strSplit = lines[i[0]].split(' ')
            strList = []
            myVals =  []
            for j in strSplit:
                if j != '':
                    strList.append(j)
            for j in range(3):
                try:
                    temp = float(strList[-(j+1)])
                except:
                    temp = None
                myVals.append(temp)
            classVals.append(myVals)
        #print(classVals)
        chems = []
        count = 1
        for i in range(3):
            r = ChemistryCEA()
            count2 = 0
            for j in classVals:
                r.__setattr__(my_vars[count2][1], j[-count])
                count2 += 1
            r.initCalculations()
            count += 1
            #print('pressure:{}'.format(r.p))
            #print('temp:{}'.format(r.t))
            #print('h:{}'.format(r.h))
            #print('m:{}'.format(r.m))
            chems.append(r)
        
        # this prints cea raw output with numbered lines
        #j = 0
        #for i in lines:
        #    print('{}:{}'.format(j, i))
        #    j+=1
        
        return chems
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
        #print('m:{}'.format(self.m))
        self.rbar = 8.31446261815324 / self.m * 1000 #ADD TO MAIN

    @staticmethod
    def create(cea, pCham, Mr, ae = None, pAmbient = None):
        #print('cea:{}\npCham:{}\nMr:{}\nae:{}\npAmbient{}'.format(cea, pCham, Mr, ae, pAmbient))
        if ae == None:
            Pratio = pCham/pAmbient
            string = cea.get_full_cea_output(Pc = pCham, MR = Mr, eps = ae, PcOvPe = Pratio, pc_units='bar', output='KJ', short_output=1)
        elif ae != None:
            string = cea.get_full_cea_output(Pc = pCham, MR = Mr, eps = ae, PcOvPe = None, pc_units='bar', output='KJ', short_output=1)
        else:
            print('chem needs a pAmbient or ae value')
        #print(string)
        lines = string.splitlines()
        lines.reverse()
        my_vars = [
            ['P,', 'p'],
            ['T,', 't'],
            ['H,', 'h'],
            ['M,', 'm'],
            ['Cp,', 'cp'],
            ['GAMMAs', 'gam'],
            ['SON', 'son'],
            ['MACH', 'mach'],
            ['Ae/At', 'aeat'],
            #['CSTAR,', 'cstar'],
            ['Ivac,', 'ivac'],
            ['Isp,', 'isp']
        ]
        classVals = []
        #lineNum = 0
        chems = []
        for i in range(3):
            r = ChemistryCEA()
            chems.append(r)
        for line in lines:
            #lineNum += 1
            #print('{}:{}|||'.format(lineNum, line))
            line = line.split(' ')
            strList = []
            for string in line:
                if string != '':
                    strList.append(string)
            if len(strList) != 0:
                for varName in my_vars:
                    if strList[0] == varName[0]:
                        #print(strList)
                        myVals =  []
                        for k in range(3):
                            try:
                                temp = float(strList[-(k+1)])
                            except:
                                temp = None
                            #print('chem:{}\ntemp:{}\nvarName:{}'.format(chems[-(k+1)], temp, varName[1]))
                            setattr(chems[-(k+1)], varName[1], temp)
                            #chems[-(k+1)].__setattr__(varName[1], temp)
        #print('chems:{}'.format(chems))
        for i in range(len(chems)):
            chems[i].initCalculations()

        return chems

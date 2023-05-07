from rocketcea.cea_obj import CEA_Obj
import math

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
    mu = None #viscosity
    Pr = None #prandtl number
    raw_cea_output = None

    def initCalculations(self):
        #print('m:{}'.format(self.m))
        self.rbar = 8.31446261815324 / self.m * 1000 #ADD TO MAIN


    @staticmethod
    def create(cea, pCham, Mr, ae = None, pAmbient = None, frozen = 1):
        #print('cea:{}\npCham:{}\nMr:{}\nae:{}\npAmbient{}'.format(cea, pCham, Mr, ae, pAmbient))
        if ae == None:
            Pratio = pCham/pAmbient
            mystring = cea.get_full_cea_output(Pc = pCham, MR = Mr, eps = ae, PcOvPe = Pratio, pc_units='bar', output='KJ', short_output=1, frozen = frozen) 
            #frozen = 1 is frozen , frozen = 0 is equilibrium
        elif ae != None:
            mystring = cea.get_full_cea_output(Pc = pCham, MR = Mr, eps = ae, PcOvPe = None, pc_units='bar', output='KJ', short_output=1, frozen = frozen)
        else:
            print('chem needs a pAmbient or ae value')
        #print(mystring)# prints raw cea output

        lines = mystring.splitlines()
        lines.reverse()
        #print(f'lines: {lines}')
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
            ['Isp,', 'isp'],
            ['RHO,', 'rho'],
            ['VISC,MILLIPOISE', 'mu'],
            ['PRANDTL', 'Pr']
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
            myline = line.split(' ')
            strList = []
            for string in myline:
                if string != '':
                    strList.append(string)
            if len(strList) != 0:
                for varName in my_vars:
                    if strList[0] == varName[0]:
                        #print(strList)
                        if strList[0] == 'RHO,': # this was necessary because for some reason the fortran output is super weird for this one row and it was going to be hard to make the parser work for this row and all the others
                            #print('I\nMADE\nIT')
                            #print(strList[0])
                            chars = [*line]
                            #print(f'chars: {chars}')
                            isvar = None
                            numstart = True
                            vals = ['', '', '']
                            wait = 0
                            #print(f'lenchars: {len(chars)}')
                            for i in range(len(chars)-1):
                                #print(f'wait: {wait}')
                                #print(f'i: {i}')
                                if wait > 0:
                                    wait -= 1
                                else:
                                    if isvar != None:
                                        if i+2 >= len(chars):
                                            vals[isvar] = float(vals[isvar]) * math.pow(10, int(chars[i]+chars[i+1]))
                                            #print(f'multiplying by exponential:  {int(chars[i]+chars[i+1])}')
                                        elif numstart == False and (chars[i] == ' ' or chars[i] == '-') and (chars[i+2] == ' '):
                                            vals[isvar] = float(vals[isvar]) * math.pow(10, int(chars[i]+chars[i+1]))
                                            #print(f'multiplying by exponential:  {int(chars[i]+chars[i+1])}')
                                            isvar += 1
                                            wait = 2
                                        elif numstart and chars[i] != ' ':
                                            vals[isvar] = vals[isvar] + chars[i]
                                            numstart = False
                                            #print(f'char to add: {chars[i]}')
                                        elif numstart == False:
                                            vals[isvar] = vals[isvar] + chars[i]
                                            #print(f'char to add: {chars[i]}')
                                    else:
                                        if chars[i].isdigit() or chars[i] == '-':
                                            #print(f'checking char: {chars[i]}')
                                            isvar = 0
                                            vals[isvar] = vals[isvar] + chars[i]
                            #print(f'vals: {vals}')
                            for i in range(3):
                                setattr(chems[i], varName[1], vals[i])
                        else:
                            for k in range(3):
                                try:
                                    temp = float(strList[-(k+1)])
                                except:
                                    temp = None
                                    print(f'chemictryCEA could not parse variable: {strList[0]}')
                                #print('chem:{}\ntemp:{}\nvarName:{}'.format(chems[-(k+1)], temp, varName[1]))
                                setattr(chems[-(k+1)], varName[1], temp)
                                #chems[-(k+1)].__setattr__(varName[1], temp)
                
        #print('chems:{}'.format(chems))
        for i in range(len(chems)):
            chems[i].initCalculations()
            setattr(chems[i], 'raw_cea_output', mystring)
        return chems


    def __repr__(self):
        return  f'''\n aeat={self.aeat} 
        rho={self.rho}
        h={self.h}
        s={self.s}
        gam={self.gam}
        g={self.g}
        u={self.u}
        m={self.m}
        mw={self.mw}
        p={self.p}
        son={self.son}
        cp={self.cp}
        t={self.t}
        ae={self.ae}
        cf={self.cf}
        ivac={self.ivac}
        mach={self.mach}
        pip={self.pip}
        isp={self.isp}
        a={self.a}
        d={self.d}
        mu={self.mu}
        Pr={self.Pr}'''
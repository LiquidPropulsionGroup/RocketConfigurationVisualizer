from rocketcea.cea_obj import CEA_Obj

class ChemistryCEA:
    def __init__(self, cea, pCham, Mr, ae = None, pAmbient = None):
        self.ox = ox
        self.fuel = fuel
        self.pCham = pCham #chamber pressure in bar
        self.Mr = Mr # mixture ratio
        self.ae = ae #exit area ratio
        self.pAmbient = pAmbient #ambient pressure in bar
        self.cea = cea
        self.string = None
        if self.ae == None:
            Pratio = self.pCham/self.pAmbient
            self.string = self.cea.get_full_cea_output(Pc = self.pCham, MR = self.Mr, PcOvPe = Pratio, pc_units='bar', output='KJ')
        elif self.ae != None:
            self.string = self.cea.get_full_cea_output(Pc = self.pCham, MR = self.Mr, eps = self.ae, pc_units='bar', output='KJ')
        else:
            print('chem needs a pAmbient or ae value')
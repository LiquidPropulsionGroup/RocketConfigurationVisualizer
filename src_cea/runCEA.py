import rocketcea.py_cea as py_cea

class RunCEA:
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
    Cstar = None
    raw_cea_output = None
    rbar = None

    def initCalculations(self):
        #print('m:{}'.format(self.m))
        return 8.31446261815324 / self.m * 1000 #ADD TO MAIN

    @staticmethod
    def create(cea, Pc, Mr, ae = None, pAmbient = None, frozen = 0, frozenAtThroat=0):
        """
        setupCards(Pc=100.0, MR=1.0, eps=40.0, subar=None, PcOvPe=None, frozen=0, ERphi=None, ERr=None, frozenAtThroat=0, short_output=0, show_transport=0, pc_units='psia', output='calories', show_mass_frac=False)
        """
        if pAmbient == None:
            PcOvPe = None
        else:
            PcOvPe = Pc/pAmbient
        cea.setupCards(Pc=Pc, MR=Mr, eps=ae, PcOvPe=PcOvPe, frozen=frozen, pc_units='bar', output='KJ', short_output=1, frozenAtThroat = frozenAtThroat)        
        chems = []
        if cea.fac_CR is not None:
            temp_num = 4
        else:
            temp_num = 3
        for i in range(temp_num):
            rrr = RunCEA()
            #self.i_injface = injector, self.i_chm = chamber, self.i_thrt = thrust, self.i_exit = exit
            # using these instead of the indicies manually allows for finite area combustor support
            rrr.p = py_cea.prtout.ppp[i]
            rrr.t = py_cea.prtout.ttt[i]                       # temperatures, t
            rrr.cp = py_cea.prtout.cpr[i]                      # heat capacities, cp
            rrr.gam = py_cea.prtout.gammas[i]               # heat capacity ratios, gammas
            rrr.son = py_cea.rockt.sonvel[i]                   # speed of sound, son
            rrr.mach = py_cea.rockt.vmoc[i]                    # mach number
            rrr.aeat = py_cea.rockt.aeat[i]                    # area exit / area throat
            rrr.Cstar = py_cea.rockt.cstr                      # Cstar
            rrr.ivac = py_cea.rockt.vaci[i]                    # isp in vacuum
            #isp = ... def estimate_Ambient_Isp(self, Pc=100.0, MR=1.0, eps=40.0, Pamb=14.7, frozen=0, frozenAtThroat=0):
            rrr.rho = py_cea.prtout.vlm[i]                     # density
            rrr.mu = py_cea.trpts.vis[i]                       # viscosity in millipoise
            rrr.h = py_cea.prtout.hsum[i]                      # enthalpy
            rrr.m = py_cea.prtout.wm[i]                        # molecular weight
            rrr.rbar = rrr.initCalculations()
            chems.append(rrr)
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
        Pr={self.Pr}
        Cstar{self.Cstar}'''
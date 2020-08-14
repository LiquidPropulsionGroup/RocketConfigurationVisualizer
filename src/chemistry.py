class Chemistry:
    #given from NASA chem
    aeat = None # area of exit/area of throat
    rho = None #density
    h = None #enthalpy
    s = None #entropy
    gam = None #molar heat capacity
    g = None #gibb's free energy
    u = None #internal energy
    m = None #molar mass
    mw = None #molecular weight = 1/m
    p = None #pressure
    son = None #speen of sound
    cp = None #specific heat capacity (constant pressure)
    t = None #temperature
    ae = None #exit area
    cf = None #function of nozzle
    ivac = None #specific impulse (exit velocity) if the rocket were in space
    mach = None #mach number
    pip = None #pressure ratio to chamber
    isp = None #the exhaust velocity of the gas (m/s)

    #calculated from self chem (calcs done in parse_initVeriables)
    rbar = None # Gas Constant per molecular weight (specific R) in kJ

    #calculated from other chems
    a = None #area
    d = None #diameter
    
    def initCalculations(self):
        #self.isp_s =     (not needed for every point)
        self.rbar = 8.31446261815324 / self.m * 1000

    @staticmethod
    def parse_initVeriables(file):
        with open(file, 'r') as tsv:
            lines = tsv.readlines()

            # we don't know what order the data will be delivered in. 'names' stores this when we receive it
            names = lines[0].split()

            # these are the values required by the 'Rocket' calculations
            for req in ["t", "gam", "ae", "isp", "m", "p"]:
                if req not in names:
                    raise Exception("Missing required fields")

            chem = []

            for line in lines[1:]:
                i = 0
                r = Chemistry()

                for val in line.split():
                    # dynamic access to 'Chemistry's fields based on the ordering stored in 'names'
                    r.__setattr__(names[i], float(val))
                    i += 1

                r.initCalculations()

                chem.append(r)

            return chem

    def __repr__(self):
        return f"(aeat: {self.aeat} " \
               f"rho: {self.rho} " \
               f"h: {self.h} " \
               f"s: {self.s} " \
               f"gam: {self.gam} " \
               f"g: {self.g} " \
               f"u: {self.u} " \
               f"m: {self.m} " \
               f"mw: {self.mw} " \
               f"p: {self.p} " \
               f"son: {self.son} " \
               f"cp: {self.cp} " \
               f"t: {self.t} " \
               f"ae: {self.ae} " \
               f"cf: {self.cf} " \
               f"ivac: {self.ivac} " \
               f"mach: {self.mach} " \
               f"pip: {self.pip} " \
               f"isp: {self.isp})"


def get_mdot() -> float:
    mdot = 0.0

    while True:
        mdotinput = input("Enter mdot: \n").strip()

        try:
            mdot = float(mdotinput)
            break
        except ValueError:
            print("Please enter a valid number")

    return mdot

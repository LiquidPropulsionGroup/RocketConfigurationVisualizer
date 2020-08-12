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


""" def genContour(endpoints, r1=0.05, convergence_angle=30, r2=0.03,r3=0.025,step=1e-6):

    functions = [
        lambda x: endpoints[0][1],  # straight line
        #lambda x: np.sqrt(abs(r1 ** 2 - (endpoints[1][1] - x) ** 2)) + endpoints[1][1] - r1,  # circle
        lambda x: np.sqrt(r1**2 - (x - endpoints[1][0]) ** 2) + endpoints[1][1] - r1,
        lambda x: -np.pi * ((convergence_angle/180 )* (x - endpoints[2][0])) + endpoints[2][1],  # straight line
        lambda x: -np.sqrt(r2 ** 2 - (x - endpoints[4][0]) ** 2) + endpoints[4][1] + r2,  # circle
        lambda x: -np.sqrt(r3 ** 2 - (x - endpoints[4][0]) ** 2) + endpoints[4][1] + r3,  # circle
        lambda x: ((endpoints[5][1] - endpoints[6][1]) / (endpoints[5][0] - endpoints[6][0])) * (x - endpoints[5][0]) +
                endpoints[5][1]  # straight line
    ]

    num = np.int32(np.rint((endpoints[6][0] - endpoints[0][0]) / step))

    x = np.array([])
    y = np.array([])

    for i,fun in enumerate(functions):
        temp_x = np.linspace(endpoints[i][0], endpoints[i+1][0], num)
        f = np.vectorize(fun)
        y = np.append(y, f(temp_x))
        x = np.append(x, temp_x)
    return x, y  """

# These are just some hardcoded test values
""" x,y = genContour(np.array(
    [[-0.190432417470023, 0.0400000000000000], 
     [-0.0618912817675829, 0.0400000000000000],
     [-0.0368912817675829, 0.0333012701892219],
     [-0.0150000000000000, 0.0206623327678020],
     [0, 0.0166430948813352],
     [0.00586669938711653, 0.0173412052778248],
     [0.0401400000000000, 0.0256163128987331]]), step=0.00001)

plt.plot(x, y)
plt.show() """

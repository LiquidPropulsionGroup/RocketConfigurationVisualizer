import math
import numpy as np
import matplotlib.pyplot as plt


class Chemistry:
    aeat = None
    rho = None
    h = None
    s = None
    gam = None
    g = None
    u = None
    m = None
    mw = None
    p = None
    son = None
    cp = None
    t = None
    ae = None
    cf = None
    ivac = None
    mach = None
    pip = None
    isp = None
    # isp_s = None this needs to be implemented as its referenced in the calc in line 54 can you 
    # testing yes 
    # hello world?

    
    def __repr__(self):
        return f"(aeat: {self.aeat} "\
                f"rho: {self.rho} "\
                f"h: {self.h} "\
                f"s: {self.s} "\
                f"gam: {self.gam} "\
                f"g: {self.g} "\
                f"u: {self.u} "\
                f"m: {self.m} "\
                f"mw: {self.mw} "\
                f"p: {self.p} "\
                f"son: {self.son} "\
                f"cp: {self.cp} "\
                f"t: {self.t} "\
                f"ae: {self.ae} "\
                f"cf: {self.cf} "\
                f"ivac: {self.ivac} "\
                f"mach: {self.mach} "\
                f"pip: {self.pip} "\
                f"isp: {self.isp})"


class Rocket:
    def __init__(self, chem, mdot):
        self.chem = chem
        self.mdot = mdot

        self.isp_s = chem.isp / 9.8  # Isp = exhaust velocity (m/s) / g (sec)
        self.rbar = 8.31446261815324 / chem.m * 1000  # Calculate specific gas constant and put in kJ scale
        self.a_thr = (self.mdot / (chem.p * 101325)) * math.sqrt(chem.t * self.rbar / chem.gam) * math.pow(
            (1 + ((chem.gam - 1) / 2)), ((chem.gam + 1) / (2 * (chem.gam - 1))))  # Throat Size Equation

        # p is in atm, conversion constant to Pa, might change to Pa later. area is in m^2

        self.a_noz = self.a_thr * chem.ae
        self.d_thr = 2 * math.sqrt(self.a_thr / math.pi)
        self.d_noz = 2 * math.sqrt(self.a_noz / math.pi)

        self.thrust = self.mdot * chem.isp  # + self.a_noz*(self.p-self.p_amb) not included as sea level expanded


def parse(file) -> [Chemistry]:
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

            chem.append(r)

        return chem


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


def getPoints(endpoints, r1=0.05, convergence_angle=30, r2=0.03,r3=0.025,step=1e-6):

    functions = [
        lambda x: endpoints[0][1], # straight line
        lambda x: np.sqrt(r1**2 - (endpoints[1][1] - x)**2) + endpoints[1][1] - r1, # circle
        lambda x: -np.pi * (convergence_angle / 180) * (x - endpoints[2][0]) + endpoints[2][1], # straight line
        lambda x: -np.sqrt(r2**2 - (endpoints[3][1] - x)**2) + endpoints[3][1] - r2, # circle
        lambda x: -np.sqrt(r3**2 - (endpoints[4][1] - x)**2) + endpoints[4][1] - r3, # circle
        lambda x: ((endpoints[6][1] - endpoints[5][1]) / (endpoints[6][0] - endpoints[5][0])) * (x - endpoints[2][0]) + endpoints[2][1] # straight line

    ]
    num = np.int32(np.rint((endpoints[6][0] - endpoints[0][0]) / step))
    
    x = np.array([])
    y = np.array([])

    for i,fun in enumerate(functions):
        temp_x = np.linspace(endpoints[i][0], endpoints[i+1][0], num)
        f = np.vectorize(fun)
        y = np.append(y, f(temp_x))
        x = np.append(x, temp_x)
    return x, y

x, y = getPoints(np.array(
    [[-0.190432417470023, 0.0400000000000000], 
     [-0.0618912817675829, 0.0400000000000000],
     [-0.0368912817675829, 0.0333012701892219],
     [-0.0150000000000000, 0.0206623327678020],
     [0, 0.0166430948813352],
     [0.00586669938711653, 0.0173412052778248],
     [0.0401400000000000, 0.0256163128987331]]), step=0.1)
plt.plot(x, y)
plt.show()

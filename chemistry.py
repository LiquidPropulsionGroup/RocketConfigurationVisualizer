import math


class Chemistry:
    aeat = 0
    rho = 0
    h = 0
    s = 0
    gam = 0
    g = 0
    u = 0
    m = 0
    mw = 0
    p = 0
    son = 0
    cp = 0
    t = 0
    ae = 0
    cf = 0
    ivac = 0
    mach = 0
    pip = 0
    isp = 0

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
        self.a_thr = (self.mdot / (chem.p * 101325)) * math.sqrt(chem.t * self.rbar / chem.gam) * (
                    1 + ((chem.gam - 1) / 2)) ^ ((chem.gam + 1) / (2 * (chem.gam - 1)))  # Throat Size Equation

        # p is in atm, conversion constant to Pa, might change to Pa later. area is in m^2

        self.a_noz = self.a_thr * chem.ae
        self.d_thr = math.sqrt(self.a_thr / math.pi * 2)
        self.d_noz = math.sqrt(self.a_noz / math.pi * 2)

        self.thrust = self.mdot * chem.isp  # + self.a_noz*(self.p-self.p_amb) not included as sea level expanded


def parse(file):
    with open(file, 'r') as csv:
        lines = csv.readlines()

        names = lines[0].split()

        for req in ["t", "gam", "ae", "isp", "m", "p"]:
            if req not in names:
                raise Exception("Missing required fields")

        rockets = []

        for line in lines[1:]:
            i = 0
            r = Chemistry()

            for val in line.split():
                r.__setattr__(names[i], float(val))
                i += 1

            rockets.append(r)

        return rockets


rockets = parse('test')

print(rockets)
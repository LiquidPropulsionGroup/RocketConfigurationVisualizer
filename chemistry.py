import math
import numpy as np


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
    # isp_s = None this needs to be implemented as its referenced in the calc in line 54

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


def generate_rockets(chemfile) -> [Rocket]:
    chemistry = parse(chemfile)

    mdot = []

    while True:
        mdotinput = input("Enter mdot: \n").strip()

        try:
            mdot = [float(mdotinput)]
            break
        except ValueError:
            if not (mdotinput.endswith("]") and mdotinput.startswith("[")):
                print("Enter either a range '[1:1:1]' or a single number '1.0'")
                continue

            mdotinput = mdotinput.lstrip("[").rstrip("]")

            vals = mdotinput.split(":")

            if len(vals) > 3 or len(vals) < 2:
                print("Must have at most 3 values, formatted [start : stride : end] or [start : end], "
                      "where stride defaults to .1")
                continue

            start = float(vals[0])
            stride = 0.1

            if len(vals) == 3:
                stride = float(vals[1])
                end = float(vals[2])
            else:
                end = float(vals[1])

            mdot = np.arange(start, end, stride)

            break

    rockets = []

    for chem in chemistry:
        for m in mdot:
            rockets.append(Rocket(chem, m))

    return rockets


rockets = generate_rockets('test')

print(rockets)
print(len(rockets))

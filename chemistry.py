import math
import numpy as np
import matplotlib.pyplot as plt

# input values
#mdot #total mass flowrate into engine (kg/s)
#Lstar #characteristic length (m)
chamber_diameter = 0.08

#hoop stress calculater
#takes: internal pressure, inside diameter of hoop, and wall thickness
def hoop_stress(internal_pressure, inside_diameter, wall_thickness):
    hoopStress = internal_pressure * inside_diameter / 2 * wall_thickness
    return hoopStress
                                                                                             
#bartz equation calculator
def bartz(d_throat, p_chamber, Cstar, d, c_p, visc, t_gas, t_wall):
    t_boundary = (t_gas + t_wall)/2
    return (0.026 / math.pow(d_throat,(0.2)) * math.pow((p_chamber/Cstar),(0.8)) * math.pow((d_throat/d),(1.8)) * c_p * math.pow(visc,(0.2)) * math.pow((t_gas/t_boundary),(0.8-0.2*0.6)))

#chamber diameter(0.08m), lambda curve, 15degree nozzle



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
    isp_s = None # this needs to be implemented as its referenced in the calc in line 54 can you 
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
    def __init__(self, chem, mdot, Lstar, inj_d):
        self.chem = chem
        self.mdot = mdot
        self.Lstar = Lstar
        self.inj_d = inj_d
        self.contourPoints = 0
        self.contour = 0
        # self.inj_r = self.inj_d/2

        # Specific impulse in seconds
        self.isp_s = self.chem.isp / 9.8

        # Gas Constant per molecular weight (specific R) in kJ
        self.rbar = 8.31446261815324 / chem.m * 1000

        # Throat Area Equation
        self.a_thr = (self.mdot / (chem.p * 101325)) * math.sqrt(chem.t * self.rbar / chem.gam) * math.pow(
            (1 + ((chem.gam - 1) / 2)), ((chem.gam + 1) / (2 * (chem.gam - 1))))

        # p is in atm, conversion constant to Pa, might change to Pa later. area is in m^2

        # Nozzle Exit Area and diameters via Expansion Ratio and 
        self.a_noz = self.a_thr * chem.ae
        self.d_thr = 2 * math.sqrt(self.a_thr / math.pi)
        self.d_noz = 2 * math.sqrt(self.a_noz / math.pi)

        # Thrust by fundamental rocket eq F = mdot * exhaust_velocity
        self.thrust = self.mdot * chem.isp  # + self.a_noz*(self.p-self.p_amb) not included as sea level expanded

        # Total Chamber Volume via Characteristic Length
        self.chamber_volume = self.Lstar * self.a_thr

        # Temporary Chamber Volume (cylindrical approximation)
        self.chamber_length = self.chamber_volume/(math.pi*(self.inj_d/2)**2)

        # Here goes the dimensions and convergence angle calculator, hardcoding 30.
        self.conv_angle = math.pi/6 #rad, 30deg
        self.divergence_angle = math.pi/12 #rad, 15deg
        genContourPoints(self)
        genContour(self)
    
    def genContourPoints(self, r1 = 0.05, r2 = 0.03, r3 = 0.025):
        
        # Radius and origin at throat
        origin_thr = np.array(0, self.d_thr/2)
        
        # start with the chamber side (left in our coordinates)
        d = np.array([r2*np.sin(self.conv_angle)], [origin_thr[1]+r2*(1-np.cos(self.conv_angle))])
        
        c_y = self.inj_d/2 - r1*(1-np.cos(self.conv_angle))
        c = np.array([(c_y - d[1])/np.tan(self.conv_angle) + d[0]],[c_y])

        b = np.array([c[0] + r1*np.sin(self.conv_angle)],[self.inj_d/2])

        inj = np.array([self.chamber_length],[self.inj_d/2])

        # Concactenate left points
        self.contourPoints = np.concatenate(inj, b, c, d, origin_thr)
        
        # Flip the array order (multiply all the x elements by -1, faster than sort)
        self.contourPoints = self.contourPoints * np.array([-1],[1])

        # now the right side
        n = np.array([r3*np.sin(self.divergence_angle)],[origin_thr[1]+r3*(1-np.cos(self.divergence_angle))])
        e = np.array([n[0]+(self.d_noz/2)/np.tan(self.divergence_angle)],[self.d_noz/2])

        # Concactenate the right points
        self.contourPoints = np.concatenate(self.contourPoints, n, e)

    def genContour(self, r1=0.05, convergence_angle=30, r2=0.03,r3=0.025,step=1e-6):
        #This is the function that draws the discrete contour

        functions = [
            lambda x: self.contourPoints[0][1],  # straight line
            #lambda x: np.sqrt(abs(r1 ** 2 - (endpoints[1][1] - x) ** 2)) + endpoints[1][1] - r1,  # circle
            lambda x: np.sqrt(r1**2 - (x - self.contourPoints[1][0]) ** 2) + self.contourPoints[1][1] - r1,
            lambda x: -np.pi * ((convergence_angle/180 )* (x - self.contourPoints[2][0])) + self.contourPoints[2][1],  # straight line
            lambda x: -np.sqrt(r2 ** 2 - (x - self.contourPoints[4][0]) ** 2) + self.contourPoints[4][1] + r2,  # circle
            lambda x: -np.sqrt(r3 ** 2 - (x - self.contourPoints[4][0]) ** 2) + self.contourPoints[4][1] + r3,  # circle
            lambda x: ((self.contourPoints[5][1] - self.contourPoints[6][1]) / (self.contourPoints[5][0] - self.contourPoints[6][0])) 
                        * (x - self.contourPoints[5][0]) + self.contourPoints[5][1]  # straight line
        ]

        num = np.int32(np.rint((self.contourPoints[6][0] - self.contourPoints[0][0]) / step))

        x = np.array([])
        y = np.array([])

        for i,fun in enumerate(functions):
            temp_x = np.linspace(self.contourPoints[i][0], self.contourPoints[i+1][0], num)
            f = np.vectorize(fun)
            y = np.append(y, f(temp_x))
            x = np.append(x, temp_x)
        self.contour = np.array([x],[y])


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

chem = parse('test')
print(chem)
#rocket = Rocket(chem, 0.5, 1.1, 0.08)


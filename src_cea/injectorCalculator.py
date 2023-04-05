import math
import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
from scipy.optimize import fsolve
from sympy import symbols, solve

'''
this class is currently for sizing a bipropellant coax swirl injector
functionality for other injector types may be added later.
the equations and derivations for the swirl injector sizing
can be found in the following paper.
"Design and Dynamics of Jet and Swirl Injectors" Vladimir Bazarov, Vigor Yang, Puneesh Puri
'''
class Injector:
    def __init__(self, mdot1, mdot2, p_f1, p_in1, p_c1, p_f2, p_in2, p_c2, alpha1, alpha2, n1, n2, rho1, rho2, nu1, nu2, l_n1, l_in1, l_s1, l_n2, l_in2, l_s2):

        #user input values
        self.mdot1 = mdot1  # mass flow rate stage 1
        self.mdot2 = mdot2  # mass flow rate stage 2
        self.p_f1 = p_f1      # pressure of the preasure feed system befor injection
        self.p_in1 = p_in1    # pressure after being injected to the swirler
        self.p_c1 = p_c1      # pressure of chamber
        self.p_f2 = p_f2      # pressure of the preasure feed system befor injection
        self.p_in2 = p_in2    # pressure after being injected to the swirler
        self.p_c2 = p_c2      # pressure of chamber        
        self.alpha1 = alpha1  # spray cone angle
        self.alpha2 = alpha2  # spray cone angle
        self.n1 = n1          # number of tangential injection passages
        self.n2 = n2          # number of tangential injection passages
        self.rho1 = rho1      # density of fluid
        self.rho2 = rho2      # density of fluid
        self.nu1 = nu1        # kinematic viscosity
        self.nu2 = nu2        # kinematic viscosity
        self.l_in1 = l_in1    # 3-6, length of tangential passages
        self.l_n1 = l_n1      # 0.5-2, length of nozzle
        self.l_s1 = l_s1      # l_s>2, length of vortex chamber
        self.l_in2 = l_in2    # 3-6, length of tangential passages
        self.l_n2 = l_n2      # 0.5-2, length of nozzle
        self.l_s2 = l_s2      # l_s>2, length of vortex chamber

        self.mu = None      # mass flow coefficient
        self.phi = None     # coefficient of passage fullness, or fractional area occupied by liquid in the nozzle
        self.h = None       # head
        self.v_sum = None   # total velocity
        self.v_un = None    
        self.v_an = None
        self.v_rn = None
        self.v_uk = None
        self.v_rk = None
        self.v_in = None
        self.r_mk = None
        self.r_mn = None
        self.R_in = None
        self.A = None

#step 1
    '''
Prescribe the spray cone angle based on the injector operating conditions
(usually between 90 and 120 deg, lower values may be used for special
cases). The geometric characteristic parameter A and the flow coefficient mu
are then determined from the plots in Fig. 32.
    '''
#alpha is gigen as an input
#with alpha angle chosen, values for phi, mu and A can be calculated
    def calc_phi_mu_A(self, alpha): #aplha needs to be in radians for this to work
        def funcPhiMu(z, alpha):
            phi = z[0]
            mu = z[1]
            F = np.empty((2))
            F[0] = np.tan(alpha) - np.sqrt(2*(1-phi)/phi) #eq 74
            F[1] = mu - phi*np.sqrt(phi/(2-phi)) #eq 62
            return F
        zGuess = np.array([0.5,0.4])
        #alpha = 60*np.pi/180 
        phi, mu = fsolve(funcPhiMu,zGuess, args=(alpha,))

        def funcA(A, phi, mu):
            return mu - 1/np.sqrt(A**2/(1-phi)+1/phi**2) #eq 61
        AGuess = 0.5
        A = fsolve(funcA, AGuess, args=(phi, mu,))
        return A, phi, mu

#step 2
#determine nozzle radius R_n
# this can be done with the calculated mu from step 1 and inpit values
# of mdot, rho, and deltaP
# the deltaP is the total pressure drop
    def calc_R_n(self, mdot, mu, rho, deltaP):
        R_n = 0.475*np.sqrt(mdot/(mu*np.sqrt(rho*deltaP))) #eq 103
        return R_n
#step 3
    '''
Specify the number of inlet passages (usually between two and four) and the
coefficient of injector opening, based on structural considerations. Then, the
radius of the inlet passage is obtained
    '''
# number of passages n is user input
# R_in is decided from structural considerations NOTE: learn more about this
# A and R_n are calculated from steps 1 and 2 respectivly
    def calc_r_in(self, R_in, R_n, n, A):
        r_in = np.sqrt(R_in*R_n/(n*A)) #eq 104
        return r_in
#step 4
# R_n, R_in, and r_in from steps 2 and 3
# l_in, l_n, l_s are user input, chosen from a select range shown below
# l_in = 3-6
# l_n = 0.5-2
# l_s > 2
# #NOTE: maybe add user input on this step so numbers can be changed through each iteration
    def calc_lengths(self, r_in, R_in, R_n, l_in, l_n, l_s):
        l_in = l_in*r_in
        l_n = l_n*R_n
        l_s = l_s*R_in
        R_s = R_in + r_in
        return l_in, l_n, l_s, R_s
#step 5
# find reynolds number in the inlet passages
# n and r_in are from previous steps
# mdot, rho, and nu are user input values
    def calc_Re(self, mdot, n, r_in, rho, nu):
        Re = 0.637*mdot/(np.sqrt(n)*r_in*rho*nu) #eq 101
        return Re
#step 6
# calculate A_eq and use that to get mu_eq and alpha_eq from the methods used in step 1
# R_in, R_n, n, r_in, Re are calculated from previous steps
    def calc_lam(self, Re):
        lam = 0.3164/(Re)**0.25 #eq 101
        return lam
    def calc_A_eq(self, R_in, R_n, n, r_in, lam):
        A_eq = R_in*R_n/(n*r_in**2+lam/2*R_in*(R_in-R_n)) #eq 100
        return A_eq
    def calc_phi_eq(self, A_eq):
        def funPhi(phi, A):
            return A - (1-phi)*np.sqrt(2)/(phi*np.sqrt(phi))#eq 92
        phiGuess = 0.5
        phi_eq = fsolve(funPhi, phiGuess, args=(A_eq,))
        return phi_eq
    def calc_mu_eq(self, phi): 
        mu = phi*np.sqrt(phi)/np.sqrt(2-phi) #eq 62
        return mu
    def calc_alpha_eq(self, phi): 
        alpha_eq = 180/np.pi*np.arctan(np.sqrt(2*(1-phi)/phi)) #eq 74
        return alpha_eq
#step 7
#Calculate the hydraulic-loss coefficient in the tangential passages
# first eps_in must be obtained
# this hydraulic loss coefficient is reltated to the inlet geometry of the inlet passages
# there are different equations for sharp edges, rounded, and angled orifice inlets
# figure 25 in the paper show these relations
# R_s, l_in, lam, and r_in are calculated from previous steps
    def calc_eps_in(self, R_s, l_in): 
        alpha = 90-180/np.pi*np.arctan(R_s/l_in)
        eps_in = -0.015*(alpha-30)+0.9 # linear equation from figure 25
        return eps_in
    def calc_eps(self, eps_in, lam, l_in, r_in):
        eps = eps_in + lam*l_in/(2*r_in) #eq 20
        return eps
#step 8
#Determine the actual flow coefficient
# mu_eq, eps, R_in, R_n, A are calculated from previous steps
#NOTE: check which A to use
    def calc_mu(self, mu_eq, eps, R_in, R_n, A):
        Rbar_in = R_in/R_n
        mu = mu_eq/np.sqrt(1+eps*mu_eq**2*A**2/Rbar_in**2) # eq 99
        return mu
#step 9
# Calculate the nozzle radius using the new approximation
# this step reuses the function from step 2 but with new value for mu

#step 10
# Calculate the geometric parameter A with the new value for R_n
# R_in, R_n, and r_in are calculated from previous steps
# n is user input
    def calc_A(self,R_in, R_n, n, r_in):
        A = R_in*R_n/(n*r_in**2) #eq 48
        return A

#step 11
# Repeat steps 1-10 until the calculated injector parameters converge.

    def calculate(self):
        alpha = self.alpha
        deltaP = self.p_f - self.p_c #NOTE: confirm this is correct
        l_in_ratio, l_n_ratio, l_s_ratio = self.l_in, self.l_n, self.l_s
        print(f"input values: \nalpha = {alpha}\ndeltaP = {deltaP}\nl_in_ratio = {l_in_ratio}\nl_n_ratio = {l_n_ratio}\nl_s_ratio = {l_s_ratio}")
        print("calculated values:")
        A, phi, mu = self.calc_phi_mu_A(alpha)
        print('A = {}\nphi = {}\nmu = {}'.format(A, phi, mu))
        R_n = self.calc_R_n(self.mdot, mu, self.rho, deltaP)
        print(f'R_n = {R_n}')
        for i in range(20):
            print()
            R_in = R_n*1.25 # consider making this a user input check for each iteration
            print(f'R_in = {R_in}')
            r_in = self.calc_r_in(R_in, R_n, self.n, A)
            print(f'r_in = {r_in}')
            l_in, l_n, l_s, R_s = self.calc_lengths(r_in, R_in, R_n, l_in_ratio, l_n_ratio, l_s_ratio)
            print('l_in = {}\nl_n = {}\nl_s = {}\nR_s = {}'.format(l_in, l_n, l_s, R_s))
            Re = self.calc_Re(self.mdot, self.n, r_in, self.rho, self.nu)
            print(f'Re = {Re}')
            lam = self.calc_lam(Re)
            print(f'lam = {lam}')
            A_eq = self.calc_A_eq(R_in, R_n, self.n, r_in, lam)
            print(f'A_eq = {A_eq}')
            phi_eq = self.calc_phi_eq(A_eq)
            print(f'phi_eq = {phi_eq}')
            mu_eq = self.calc_mu_eq(phi_eq)
            print(f'mu_eq = {mu_eq}')
            alpha_eq = self.calc_alpha_eq(phi_eq)
            print(f'alpha_eq = {alpha_eq}')
            eps_in = self.calc_eps_in(R_s, l_in)
            print(f'eps_in = {eps_in}')
            eps = self.calc_eps(eps_in, lam, l_in, r_in)
            print(f'eps = {eps}')
            mu = self.calc_mu(mu_eq, eps, R_in, R_n, A)
            print(f'mu = {mu}')
            R_n = self.calc_R_n(self.mdot, mu, self.rho, deltaP)
            print(f'R_n = {R_n}')
            A = self.calc_A(R_in, R_n, self.n, r_in)
            print(f'A = {A}')
            phi = self.calc_phi_eq(A)
            print(f'phi_eq = {phi_eq}')
            mu = self.calc_mu_eq(phi)
            print(f'mu_eq = {mu_eq}')
            alpha = self.calc_alpha_eq(phi)
            print(f'alpha_eq = {alpha_eq}')

        self.R_in = R_in
        self.R_s = R_s
        self.R_n = R_n
        self.r_in = r_in
        self.l_in = l_in
        self.l_n = l_n
        self.l_s = l_s
        self.alpha = alpha_eq
        self.mu = mu
        self.A = A
        self.eps = eps
        self.Re = Re

    '''
this method uses much the same procedure as the monopropellant element sizing except that it is ment for
designing two nested elements denoted as stage 1 being the inner element and stage 2 being the outer.
there is also an option for internal or external mixing which involve different calculations.
there are 2 calculations senerios which effect the sizing of the stage 2 element:
senerio 1 the wall of the stage 1 nozzle is within the gaseouse core of the stage 2 nozzle flow
senerio 2 the stage 1 nozzle is submerged and effecting the flow of the stage 2 nozzle
    '''
    def calculateBipropellant(self, senerio):
        # stage 1 sizing calculations

        # step 1
        '''
1) Prescribe the spray cone angles 2alpha_2 and 2alpha_1, according to the empirical
condition 2alpha_1 — 2alpha_2 = 10 to 15 deg based on injector operating conditions.
With these values and the correlation given in Fig. 34a, find the geometric characteristic parameters, A1 and A2.
The flow coefficients of stages 1 and 2, mu_1 and mu_2, are then determined from Fig. 34b
        '''
        # step 2
        '''
2) Calculate the nozzle radii R_n1 and R_n2 from Eq. (103), and determine the
tangential-entry radii r_in1 and r_in2 from Eq. (104).
        '''

        # stage 2 sizing calculations
        if senerio == 0:

        else:

        #mixing element sizing calculations


#--------------------------------------------------------------------------------------
    def swirlgraphs(self):
        def y(x, A):
            return 1 / np.sqrt(A**2/(1-x) + 1/x**2)

        # Define the range of x values to plot
        x = np.linspace(0.01, 0.99, 1000)

        # Define the values of A to plot
        A_values = [2.0, 0.8, 0.4, 0.2]

        # Create a figure and axes object
        fig, ax = plt.subplots()

        # Plot the function for each value of A
        for A in A_values:
            ax.plot(x, y(x, A), label=f"A={A}")
        ax.plot(x, x*np.sqrt(x/(2-x)), label=f"ideal")

        # Set the title and axis labels
        ax.set_title("Graph of y = 1/(A^2/(1-x)+1/x^2)")
        ax.set_xlabel("x")
        ax.set_ylabel("y")

        # Add a legend and gridlines
        ax.legend()
        ax.grid()

        # Show the plot
        plt.show()

    def __str__(self):
        return '''
            R_in = {} m radial location of tangential inlet passages
            R_s = {} m radius of vortex chamber
            R_n = {} m radius of nozzle
            r_in = {} m radius of inlet passages
            l_in = {} m lenth of tangential inlet passages
            l_n = {} m length of nozzle
            l_s = {} m length of vortex chamber
            alpha = {}
            mu = {}
            A = {}        
            eps = {}
            Re = {}
        '''.format(
            self.R_in,
            self.R_s,
            self.R_n,
            self.r_in,
            self.l_in,
            self.l_n,
            self.l_s,
            self.alpha,
            self.mu,
            self.A,
            self.eps,
            self.Re
        )

if __name__ == "__main__": #test values
    mdot1 = 0.7/6       #mass flow of stage 1 of one injection element
    mdot2 = 0.3/6       #mass flow of stage 2 of one injection element
    #pressures in pascals
    p_f = 24*10**5
    p_in = 23*10**5 # not currently being used
    p_c = 20*10**5
    alpha = 60*np.pi/180 #in radians
    n = 3
    l_n = 4.5   # l_in = 3-6
    l_in = 1    # l_n = 0.5-2
    l_s = 3     # l_s > 2
    rho = 997   #in kg/m^3
    nu = 0.6*10**(-6) #in m^2/s

    my_swirl_injector = Injector(mdot1, mdot2, p_f, p_in, p_c, alpha, n, l_n, l_in, l_s, rho, nu)
    my_swirl_injector.calculate()
    print(my_swirl_injector)

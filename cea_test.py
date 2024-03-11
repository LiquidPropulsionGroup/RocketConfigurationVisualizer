from rocketcea.cea_obj import CEA_Obj
import rocketcea.cea_obj as cea_obj
import rocketcea.py_cea as py_cea
import time
from CEA_Wrap import Fuel, Oxidizer, RocketProblem
import platform

#: RocketCEA wraps the NASA FORTRAN CEA code to calculate Isp, cstar, and Tcomb
#: This object wraps the English unit version of CEA_Obj to enable desired user units.
#: Same as CEA_Obj with standard units except, input and output units can be specified.
#:  parameter             default             options
#: isp_units            = 'sec',         # N-s/kg, m/s, km/s
#: cstar_units          = 'ft/sec',      # m/s
#: pressure_units       = 'psia',        # MPa, KPa, Pa, Bar, Atm, Torrs
#: temperature_units    = 'degR',        # K, C, F
#: sonic_velocity_units = 'ft/sec',      # m/s
#: enthalpy_units       = 'BTU/lbm',     # J/g, kJ/kg, J/kg, kcal/kg, cal/g
#: density_units        = 'lbm/cuft',    # g/cc, sg, kg/m^3
#: specific_heat_units  = 'BTU/lbm degR' # kJ/kg-K, cal/g-C, J/kg-K (# note: cal/g K == BTU/lbm degR)
#: viscosity_units      = 'millipoise'   # lbf-sec/sqin, lbf-sec/sqft, lbm/ft-sec, poise, centipoise
#: thermal_cond_units   = 'mcal/cm-K-s'  # millical/cm-degK-sec, BTU/hr-ft-degF, BTU/s-in-degF, cal/s-cm-degC, W/cm-degC
#: fac_CR, Contraction Ratio of finite area combustor (None=infinite)
#: if make_debug_prints is True, print debugging info to terminal.
#__init__(propName='', oxName='', fuelName='', useFastLookup=0, makeOutput=0, isp_units='sec', cstar_units='ft/sec', pressure_units='psia', temperature_units='degR', sonic_velocity_units='ft/sec', enthalpy_units='BTU/lbm', density_units='lbm/cuft', specific_heat_units='BTU/lbm degR', viscosity_units='millipoise', thermal_cond_units='mcal/cm-K-s', fac_CR=None, make_debug_prints=False)

C = CEA_Obj(oxName='LOX', fuelName='RP1')

#full_cea = C.get_full_cea_output(Pc=15.2, MR=2.3, eps=None, PcOvPe=15.2/0.9, pc_units='bar', output='KJ', short_output=1 )
#print('this is the cea string:------------------------------------------------------------------------------------------------------------------')
#print(full_cea)

#cea_obj.print_py_cea_vars()
#cea_obj.clearCache()                       #not working
#print(cea_obj.getCacheDict())
print(cea_obj.get_rocketcea_data_dir())    #not working

#testing setup cards
"""
setupCards(Pc=100.0, MR=1.0, eps=40.0, subar=None, PcOvPe=None, frozen=0, ERphi=None, ERr=None, frozenAtThroat=0, short_output=0, show_transport=0, pc_units='psia', output='calories', show_mass_frac=False)
#: Pc = combustion end pressure (psia)
#: eps = Nozzle Expansion Area Ratio
#: if PcOvPe has a value, use it instead of eps to run case
#: ERphi = Equivalence ratios in terms of fuel-to-oxidant weight ratios.
#: ERr = Chemical equivalence ratios in terms of valences.
#: pc_units = 'psia', 'bar', 'atm', 'mmh'(mm of mercury)
#: frozen flag (0=equilibrium, 1=frozen)
#: frozenAtThroat flag, 0=frozen in chamber, 1=frozen at throat
"""
#run cea
st = time.time()
C.setupCards(Pc=15.2, MR=2.3, eps=None, PcOvPe=15.2/0.9, frozen=0, pc_units='bar', output='KJ', short_output=1)
et = time.time()
rt = et-st
print(f'rocketcea run time = {rt}s')

if platform.system() == "Windows":
    h2 = Fuel("H2(L)", temp=20) # Liquid H2, Temp in Kelvin
    lox = Oxidizer("O2(L)", temp=90)
    # Rocket at 2000psi and supersonic area ratio of 5
    st = time.time()
    problem1 = RocketProblem(pressure=2000, materials=[h2, lox], phi=1, sup=5)
    results = problem1.run()
    et = time.time()
    rt = et-st
    print(f'CEA_Wrap run time = {rt}s')
else:
    print("windows not detected. skipping CEA_Wrap")

"""
py_cea.rockt.vaci[ self.i_thrt ] =  0.0 # Vacuum Isp at throat
py_cea.rockt.vaci[ self.i_exit ] =  0.0 # Vacuum Isp at exit
py_cea.rockt.cstr = 0.0                 # cstar
py_cea.prtout.ttt[ self.i_chm ] = 0.0   # chamber temperature
py_cea.rockt.app[ self.i_thrt ] = 0.0   # Pc/Pt
py_cea.rockt.app[ self.i_exit ] = 0.0   # Pc/Pe
py_cea.rockt.aeat[ self.i_exit ] = 0.0  # exit area/throat area
py_cea.rockt.vmoc[ self.i_exit ]  = 0.0 # MACH, mach number
py_cea.miscr.eqrat = 0.0                # equivalence ratio

for i in range(3):
    py_cea.rockt.sonvel[i] = 0.0    # SON, sonic velocity, speed of sound
    py_cea.prtout.hsum[i] = 0.0     # enthalpies
    py_cea.prtout.wm[i] = 0.0       # MolWt, molecular weight
    py_cea.prtout.gammas[i] = 0.0   # gamma, heat capacity ratio
    py_cea.prtout.vlm[i] = 0.0      # densities
    py_cea.prtout.cpr[i] = 0.0      # heat capacities Cp
    
    py_cea.trpts.vis[i] = 0.0    # viscosity
    py_cea.trpts.cpeql[i] = 0.0  # equilibrium specific heat
    py_cea.trpts.coneql[i] = 0.0 # equilibrium thermal conductivity
    py_cea.trpts.preql[i] = 0.0  # equilibrium prandtl number
    py_cea.trpts.cpfro[i] = 0.0  # frozen specific heat
    py_cea.trpts.confro[i] = 0.0 # frozen thermal conductivity
    py_cea.trpts.prfro[i] = 0.0  # frozen prandtl number
"""

# ['P,', 'p'],
# ['T,', 't'],
# ['Cp,', 'cp'],
# ['GAMMAs', 'gam'],
# ['SON', 'son'],
# ['MACH', 'mach'],
# ['Ae/At', 'aeat'],
# ['CSTAR,', 'cstar'],
# ['Ivac,', 'ivac'],
# ['Isp,', 'isp'],
# ['RHO,', 'rho'],
# ['VISC,MILLIPOISE', 'mu'],
# ['PRANDTL', 'Pr']

"""
# ------------ Finite Area Combustor String ------

#print('In setupCards, fac_CR =', self.fac_CR)
if self.fac_CR is not None:
    N += 1
    fac_str = ' fac ' + make_inp_str(" ac/at=", self.fac_CR )
            
    set_py_cea_line(N, fac_str )
    
    self.i_injface = 0
    self.i_chm = 1
    self.i_thrt = 2
    self.i_exit = 3
else:
    self.i_injface = 0
    self.i_chm = 0
    self.i_thrt = 1
    self.i_exit = 2
"""

#print variables-----------------------------------------------------------------------------------------------------
#self.i_injface = injector, self.i_chm = chamber, self.i_thrt = thrust, self.i_exit = exit
# using these instead of the indicies manually allows for finite area combustor support
temperatures = py_cea.prtout.ttt[C.i_chm] 
cp = py_cea.prtout.cpr
gammas = py_cea.prtout.gammas
son = py_cea.rockt.sonvel
mach = py_cea.rockt.vmoc
aeat = py_cea.rockt.aeat
Cstar = py_cea.rockt.cstr
ivac = py_cea.rockt.vaci
#isp = ... def estimate_Ambient_Isp(self, Pc=100.0, MR=1.0, eps=40.0, Pamb=14.7, frozen=0, frozenAtThroat=0):
rho = py_cea.prtout.vlm
mu = py_cea.trpts.vis
h = py_cea.prtout.hsum
m = py_cea.prtout.wm
"""
        if self.fac_CR is not None:
            sonicList = list(py_cea.rockt.sonvel[1:4])
        else:
            sonicList = list(py_cea.rockt.sonvel[:3])
"""

print(f'temperatures = {temperatures}')
print(f'cp = {cp}')
print(f'gammas = {gammas}')
print(f'son = {son}')
print(f'mach = {mach}')
print(f'aeat = {aeat}')
print(f'Cstar = {Cstar}')
print(f'ivac = {ivac}')
#print(f'isp = {isp}')
print(f'rho = {rho}')
print(f'mu = {mu}')
print(f'h = {h}')
print(f'm = {m}')





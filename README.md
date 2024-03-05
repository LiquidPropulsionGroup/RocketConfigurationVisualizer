# To Do List
1. calculate and add isp engine prints
2. 2d cooling channel calculator
3. film cooling calculator

# Rocket Configuration Visualizer

<!-- Add information here explaining the visualizer -->

This code sizes a rocket engine contour based on design specifications and thermochemical data imported from NASA CEA. It then computes flow properties based on quasi-1D isentropic flow assumptions and then estimates the wall heat transfer coefficient based on the Bartz approximation. In the future, this program will implement Cpropep directly and allow rapid sizing and optimization of a rocket engine through a GUI.

## Python setup for Development

1. Install python (It can be the newest one)
2. Pull a clone of this repository
3. Create a python virtual enviorment with the following command in the root directory of the Rocket Configuration Visualizer. `python -m venv venv/`
4. Finally install all the dependencys once in your virtual enviorment of python with `pip3 install -r requirements.txt` (do this in the Terminal of VS Code, not cmd)
5. install rocketcea library. this is different per operating system but this  goes over how to do it(https://rocketcea.readthedocs.io/en/latest/quickstart.html) 
note: this library can be rather difficult to install correctly on windows
sometimes the interpreter won't be able to find rocketcea when using a virtual enviernment; if that happens you can either not use a venv or use the following code to manually locate rocketcea
# import sys
# sys.path.insert(0, r"C:\Users\<insert_user>\AppData\Local\Programs\Python\Python310\Lib\site-packages")
the specific file location and names will be defferent depenging on the version and location of your python files

## Basic Overview

The engine class takes the input veriables creates the optimal engine contour for the given constraints and throttle range. To create this contour chemistryCEA is used. this class takes the raw data output, which is from the nasaCEA code that is run through the rocketcea library, and parces the data into veriables. durring this optimization of the contour, thrustLevel instances are created for the max and min thrusts. when thrustLevel is instantiated basic veriables are calculated. from there we go back to rocket where the contour is made into an array and this array is then given back to the thrust level instances in the heatCalcs method, which calculates the heat flux at every point int eh contour array using the bartz approximation and other simple operations. at this point all of the calculations have been done and you can view the data using the variablesDisplay and graphDisplay methods in engine. 

## Conventions

### Rocket Constructor

To analyze a rocket configuration, a Rocket object needs to be instantiated, with the before mentioned values passed to it. The [constructor](https://en.wikipedia.org/wiki/Constructor_(object-oriented_programming)) currently autoruns the entire analysis functionality, meaning that instantiating a Rocket object automatically finishes the analysis. This is done for ease of use as we debug the code and mainly use it for only one configuration. The constructor can be [overloaded](https://en.wikipedia.org/wiki/Function_overloading) to have different functionality in the future for partial analysis or some other function.

First, the class fills out the chemistry data from array of chems objects into named points in the engine for ease of manipulation. For example, to access the temperature at the throat, one would use `self.thr.t` , and to find the Mach number at the exhaust plane one would use `self.exit.mach`

The rest of the arguments are written into methods, and some other properties are instantiated as None types or as empty lists.

### my_contourPoints

my_contourPoints is a function that generates x and y positions for critical points that are used to interpolate the rocket contour and are then exported out as a txt for Solidworks. The generated contourPoints array contains an axial position (x) and a radius (y). The origin is at the throat, therefore -X is the chamber and +X is the nozzle. Point a is also called inj in the solidworks document due to LPG engine CAD nomenclature. 

### genContour

This function interpolates the contour for numerical calculations. It uses lambda functions: the contour consists of straight lines and circles per Sutton. Nozzle is conical, with 80% bell optimization capability to be added in the future.

## Math Review

### inputs

ox  # full propellant selection is availible at https://rocketcea.readthedocs.io/en/latest/propellants.html
fuel
customFuel = [
    "Isopropanol70", 
    """fuel C3H8O-2propanol C 3 H 8 O 1    wt%=70.0
h,cal=-65133.0     t(k)=298.15   rho=0.786
fuel water H 2.0 O 1.0  wt%=30.0
h,cal=-68308.0  t(k)=298.15 rho,g/cc = 0.9998"""
]
pMaxCham    #max thrust chamber pressure in bar
Mr  #propellant mixture ratio
pAmbient
pMinExitRatio
#set veriables
mdotMax  #max thrust mass flow rate
filmCoolingPercent
Lstar
Dcham
conv_angle
div_angle
wall_temp
fuel_delta_t
fuel_cp
r1
r2
r3
step
nozzle_type
### outputs



### CEA and parser
the thermochemistry part of this code is done through the CEArocket package which is a wraper for the NASA CEA FORTRAN code. the chemistry CEA
#### CEA inputes
fuel
pCham
Mr
pAmbient or ea
frozen or equilibrium

#### CEA outputs


### chamber nozzle contour
points a,b,c,d,o,n,e from left to right where the exit is to the right are used to define the nozzle contour
r1, r2, r3 are the radii ratios compared to the throat for the throat and chamber curvatures
o stands for origin
e for exit
the nozzle can be either connical or a bell nozzle:
for connical a 15 degree strait line is drawn from point n to the exit radius

#### 80% bell nozzle calcs
for the bell nozzle a parabola is drawn with the x distance being defined as 80% of the length of the respective connical nozzle. the paper below has detailes for the equations we used for making this.
"Design and analysis of contour bell nozzle and comparison with dual bell nozzle"
https://core.ac.uk/download/pdf/154060575.pdf

### mach area

( 2 / (self.thr.gam + 1) * ( 1 + (self.thr.gam - 1)/2 * mach**2 ))**((self.thr.gam+1)/(2*(self.thr.gam-1))) - mach * area_ratio
https://web.mit.edu/16.unified/www/SPRING/propulsion/notes/node104.html
IMPORTANT NOTE: area ratio might be wrong, might need to be inversed

### bartz
(0.026 / math.pow(d_throat, 0.2) * math.pow((p_chamber / c_star), 0.8) * math.pow((d_throat / d),
            1.8) * c_p * math.pow(visc, 0.2) * math.pow((t_gas / t_boundary), (0.8 - 0.2 * 0.6)))

h_g_arr[1,i] = self.bartz(self.thr.d, self.cham.p, self.Cstar, self.contour[1,i]*2, self.cham.cp*1000, 0.85452e-4, self.temp_arr[1,i], self.wall_temp)

https://arc.aiaa.org/doi/pdf/10.2514/8.12572

### adiabatic pressure, temperature, density

t_stag = self.cham.t * (1 + ((gam-1)/2 * self.cham.mach**2))
myreturn = t_stag * (1 + ((gam-1)/2 * mach**2))**(-1)

p_stag = self.cham.p * (1 + ((gam-1)/2 * self.cham.mach**2))**(gam/(gam-1))
myreturn = p_stag * (1 + ((gam-1)/2 * mach**2))**(-gam/(gam-1))

d_stag = self.cham.rho * (1 + ((gam-1)/2 * self.cham.mach**2))**(1/(gam-1))
myreturn = d_stag * (1 + ((gam-1)/2 * mach**2))**(-1/(gam-1))


### other heat calcs

### wall thermal calcs

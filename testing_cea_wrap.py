
### This simple example file goes over the different ways to define problems, material ratios, and materials
from pprint import pprint

# Import only the things we'll need for this problem. You could also do from CEA_Wrap import * if you like polluting your global namespace
from CEA_Wrap import Fuel, Oxidizer, HPProblem, RocketProblem

#############################################################################
## First we'll show the different ways to define problems
#############################################################################

# For all these we will consider stoichiometric combustion of methane and oxygen
# Define our materials, weight percentages default to 100%
'''
methane = Fuel("CH4")
o2 = Oxidizer("O2")

## NOTE: All these methods are equivalent and will produce identical output
## NOTE: The only things you need to specify for a problem are some form of o/f ratio and materials
##         all other parameters have defaults, such as a pressure of 1000 psi and a filename of "my_output"

## Problem Definition 1: All parameters in problem constructor
problem = HPProblem(phi=1, pressure=500, materials=[methane, o2], filename="HP_Problem")
data1 = problem.run()

## Problem Definition 2: Specify parameters with object
problem = HPProblem()
problem.set_phi(1)
problem.set_pressure(500)
problem.set_materials([methane, o2])
problem.set_filename("HP_Problem")
data2 = problem.run()

## Problem Definition 3: Mix and match. Any parameter can be set in constructor or from object
##   This is particularly useful when you are modifying a problem within a loop, such as looping through pressures
problem = HPProblem(materials=[methane, o2], filename="HP_Problem")
problem.set_phi(1)
problem.set_pressure(500)
data3 = problem.run()

## Problem Definition 4: As the others, but you can also specify materials when you call the problem
problem = HPProblem(phi=1, pressure=500, filename="HP_Problem")
data4 = problem.run(methane, o2)

# Let's double check they are all the same!
assert data1 == data2 == data3 == data4

# Now we'll print out the data so you can see what it produces when you run it
print("Displaying results for stoichiometric methane/oxygen combustion")
pprint(data1)
'''
#-----------------------------------------------------------------------------------------------------------------------------------------
RP1 = Fuel("CH4(L)", temp = 110)
LO2 = Oxidizer("O2(L)", temp = 90)

## Ratios Example 3: Helper function: absolute o/f ratios
# So doing rocket problems with absolute percentages specified (not as percentage of fuel/oxidizer and o/f ratio specified)
#   I got annoyed, so I wrote a helper function that calculates the o/f ratio from wt_percent of all constituents
# So we'll be doing a rocket problem with AP/Aluminum/HTPB
# Define materials
aluminum = Fuel("AL(cr)", wt_percent=12) # (cr) for "crystalline" or condensed phase
htpb = Fuel("HTPB", wt_percent=14) # This was added at Purdue so doesn't include (cr) in the name
ap = Oxidizer("NH4CLO4(I)", wt_percent=74) # ammonium perchlorate (form I, specified at room temperature)

# Specify a RocketProblem with supersonic expansion ratio of 15 at chamber pressure of 1500 psi
# We don't specify our phi or o_f as we want it to calculate that for us
problem = RocketProblem(pressure=220.44, Pc_Pe = 3.075, materials=[RP1, LO2], o_f = 2.7)
# sets the o_f ratio by summing wt_percent for each oxidizer and  fuel
#  and dividing by the sum for fuels
# If you change material wt_percents, remember to call this again!

#problem.set_absolute_o_f() # sets o/f to (74)/(12+14)
data4 = problem.run()
pprint(data4)
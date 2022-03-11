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
sometimes the interpreter won't be able to find rocketcea when using a venv, if that happens you can either not use a venv or use the following code to manually locate rocketcea
# import sys
# sys.path.insert(0, r"C:\Users\<insert_user>\AppData\Local\Programs\Python\Python310\Lib\site-packages")
the specific file location and names will be defferent depenging on the version and location of your python files

## Basic Overview

Chemistry class is instantiated per every imported chemistry data point from CEA. Rocket class generates a unique rocket contour and computes flow and heat properties upon instantiation. Mulitple Rocket objects represent different rocket sizes. There are inline comments where possible in the code. Further documentation coming soon.

## Conventions

### Rocket Constructor

To analyze a rocket configuration, a Rocket object needs to be instantiated, with the before mentioned values passed to it. The [constructor](https://en.wikipedia.org/wiki/Constructor_(object-oriented_programming)) currently autoruns the entire analysis functionality, meaning that instantiating a Rocket object automatically finishes the analysis. This is done for ease of use as we debug the code and mainly use it for only one configuration. The constructor can be [overloaded](https://en.wikipedia.org/wiki/Function_overloading) to have different functionality in the future for partial analysis or some other function.

First, the class fills out the chemistry data from array of chems objects into named points in the engine for ease of manipulation. For example, to access the temperature at the throat, one would use `self.thr.t` , and to find the Mach number at the exhaust plane one would use `self.exit.mach`

The rest of the arguments are written into methods, and some other properties are instantiated as None types or as empty lists.

### my_contourPoints

my_contourPoints is a function that generates x and y positions for critical points that are used to interpolate the rocket contour and are then exported out as a txt for Solidworks. The generated contourPoints array contains an axial position (x) and a radius (y). The origin is at the throat, therefore -X is the chamber and +X is the nozzle. Point a is also called inj in the solidworks document due to LPG engine CAD nomenclature. 

### genContour

This function interpolates the contour for numerical calculations. It uses lambda functions: the contour consists of straight lines and circles per Sutton. Nozzle is conical, with 80% bell optimization capability to be added in the future.


# SERIAL CODE

This directory contains all the necessary codes to perform a VdW simulation of a Xe gas in serial. 

## TABLE OF CONTENTS

1. [ REQUIREMENTS ](#1-req)
2. [ EXECUTION](#2-ex)
3. [ RESULTS](#3-res)

<a name="1-req"></a>
### REQUIREMENTS

All codes used for the simulation have been created using _Fortran 90_. The visualization code `Block_average_and_plots.py` is in _python (3.11.2)_. In order to execute the code, the _[gfortran]_ compiler is needed. If you use a different compiler, change the `FC` varaible in the `Makefile`.

In order to visualize the trajectories of the simulation as well as the initial configuration we used [VMD] but any visualization program that can read `.xyz` files can be used.

[gfortran]: https://fortran-lang.org/en/learn/os_setup/install_gfortran/
[VMD]: https://www.ks.uiuc.edu/Development/Download/download.cgi?PackageName=VMD

<a name="2-ex"></a>
### EXECUTION

To execute the program in your machine use `make` to compile the code and then use `make run` to start the simulation. If you want to plot your results, use `make plot` after you have runned your simulation. This will save all the plots in `.pdf` format. Notice that when you run the code the directories `/results` and `/plots` are automatically generated. As their names say, all the results generated from the simulation will go to the directory `/results` while all plots will go to `/plots`. 

Before running your simulation, you may want to change the initial conditions for the system configuration. In order to do it, you should change the values of every varaible in the `input.mdp` file. Notice that a file named `get_input.dat` is also generated but it is an input file for the ploting code and should not be changed if you want to plot the last generated data.

Use `make help` to get aditional information about the options of the make file.

<a name="3-res"></a>
### RESULTS
After running your simulation with the default input file your `/results` directory should contain the following files:

```Markdown

125_dynamics_1.352.xyz    (contains the frames of the simulation in .xyz format)
125_energy_1.352.dat      (Contains in this order time,potential energy,kinetic energy,total energy,instant temperature,momentum,pressure and msd of every frame)
125_gdr_1.352.dat         (Contains the GDR and distances)
thermo_prop.dat           (Contains the average values of all the thermodynamic properties)
initial_conf_1.352.xyz    (Contains the frame of the starting configuration)

```

Notice that the files start with the number of particles of the system $M^3$ (125) and end with the value of the `density` (1.352) of the `input.dat` file respectively. 

After running `Block_average_and_plots.py`, the `/plots` directory may contain all the plots required. It is important that the `/results` directory contains all the specified files and that the `get_input.dat` file has been correctly generated! 



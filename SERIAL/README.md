# SERIAL CODE

This directory contains all the necessary codes to perform a VdW simulation of a Xe gas in serial. 

### REQUIREMENTS

All codes used for the simulation have been created using _Fortran 90_. The visualization code `Block_average_and_plots.py` is in _python (3.11.2)_. In order to execute the code, the _[gfortran]_ compiler is needed. If you use a different compiler, change the `FC` varaible in the `Makefile`.

In order to visualize the trajectories of the simulation as well as the initial configuration we used [VMD] but any visualization program that can read `.xyz` files can be used.

[gfortran]: https://fortran-lang.org/en/learn/os_setup/install_gfortran/
[VMD]: https://www.ks.uiuc.edu/Development/Download/download.cgi?PackageName=VMD


### EXECUTION

To execute the program in your machine use `make` to compile the code and then use `make run` to start the simulation. If you want to plot your results, use `make plot` after you have runned your simulation. Notice that when you run the code the directories `/results` and `/plots` are automatically generated. As their names say, all the results generated from the simulation will go to the directory `/results` while all plots will go to `/plots`. 

Before running your simulation, you may want to change the initial conditions for the system configuration. In order to do it, you should change the values of every varaible in the `input.dat` file. Notice that a file named `get_input.dat` is also generated but it is an input file for the ploting code and should not be changed if you want to plot the last generated data.

Use `make help` to get aditional information about the options of the make file.


### RESULTS
After running your simulation with the default input file your `/results` directory should contain the following files:

```Markdown

125_dynamics_0.800.xyz    (Dynamics of the simulation)
125_energy_0.800.dat      (Contains in this order time,potential energy,kinetic energy,total energy,instant temperature,momentum,pressure and msd of every frame)
125_gdr_0.800.dat         (Contains the GDR and distances)
125_msd_0.800.dat         (Contains time and MSD)
thermo_prop.dat           (Contains the average values of all the thermodynamic properties)
initial_conf_0.800.xyz    (Contains the frame of the starting configuration)

```

Notice that the files start with the number of particles of the system $M^3$ (125) and end with the value of the `density` (0.800) of the `input.dat` file respectively. 



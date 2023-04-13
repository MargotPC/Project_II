# PARALLEL CODE

This directory contains all the necessary codes to perform a VdW simulation of a Xe gas in parallel using MPI. 

## TABLE OF CONTENTS

1. [ REQUIREMENTS ](#1-req)
2. [ EXECUTION](#2-ex)
3. [ RESULTS](#3-res)
4. [ SPEEDUP](#4-speed)

<a name="1-req"></a>
### REQUIREMENTS

All codes used for the simulation have been created using _Fortran 90_. The visualization code `Block_average_and_plots.py` is in _python (3.11.2)_. In order to execute the code, the _[mpi90]_ compiler is needed. If you use a different compiler, change the `FC` varaible in the `Makefile`. If the code is running on a cluster, you may need to update the python version of the cluster. The `run.sub` file may need some adaptations depending on the cluster.

In order to visualize the trajectories of the simulation as well as the initial configuration we used [VMD] but any visualization program that can read `.xyz` files can be used.

[mpi90]: [https://fortran-lang.org/en/learn/os_setup/install_gfortran/](https://edu.itp.phys.ethz.ch/hs12/programming_techniques/openmpi.pdf)
[VMD]: https://www.ks.uiuc.edu/Development/Download/download.cgi?PackageName=VMD

<a name="2-ex"></a>
### EXECUTION

When executing the code in a cluster, use the `Makefile` in order to perform your desired actions. First of all is necessary to compile the code in the cluster to generate the binary files. This can be done using the `make compile` (or just `make`) command on the cluster terminal. To run the code simply use `make run` on the cluster terminal. Notice that the default settings use only 4 cores but this can be changed in the `run.sub` file modifying the `#$ -pe smp` option.

Before running your simulation, you may want to change the initial conditions for the system configuration. In order to do it, you should change the values of every varaible in the `input.dat` file. Notice that a file named `get_input.dat` is also generated but it is an input file for the ploting code and should not be changed if you want to plot the last generated data.

Use `make help` to get aditional information about the options of the make file.

<a name="3-res"></a>
### RESULTS
After running your simulation with the default input file your `/results` directory should contain the following files:

```Markdown

125_dynamics_1.352.xyz    (contains the frames of the simulation in .xyz format)
125_energy_1.352.dat      (Contains in this order time,potential energy,kinetic energy,total energy,instant temperature,momentum,pressure and msd of every frame)
125_gdr_1.352.dat         (Contains the GDR and distances)
125_msd_1.352.dat         (Contains time and MSD)
thermo_prop.dat           (Contains the average values of all the thermodynamic properties)
initial_conf_1.352.xyz    (Contains the frame of the starting configuration)

```

Notice that the files start with the number of particles of the system $M^3$ (125) and end with the value of the `density` (1.352) of the `input.dat` file respectively. 

After running `Block_average_and_plots.py`, the `/plots` directory may contain all the plots required. It is important that the `/results` directory contains all the specified files!

<a name="4-Speed"></a>
### SPEEDUP

Once you have used the `make runAll` option of the `Makefile` you can obtain all the elapsed times by using `make SR`. This will generate `times.dat`, a file containing the number of processors and the elapsed time of every simulation. You may use this data to perform the speedup of the system. 

# SERIAL CODE

This directory contains all the necessary codes to perform a VdW simulation of a Xe gas in serial. 

### REQUIREMENTS

All codes used for the simulation have been created using _Fortran 90_. The visualization code `Block_average_and_plots.py` is in python (3.11.2). In order to execute the code, the _[gfortran]_ compiler is needed. If you use a different compiler, change the `FC` varaible in the `Makefile`.

[gfortran]: https://fortran-lang.org/en/learn/os_setup/install_gfortran/

### EXECUTION

To execute the program in your machine use `make` to compile the code and then use `make run` to start the simulation. If you want to plot your results, use `make plot` after you have runned your simulation. Notice that when you run the code the directories `/results` and `/plots` are automatically generated. As they names say, all the results generated from the simulation will go to the directory `/results` while all plots will go to `/plots`. 




### RESULTS
After running your simulation your `/results` directory should look like this:

```Markdown

125_dynamics_0.800

```


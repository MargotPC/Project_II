# Project_II

## Table Contents

1. [ List of team members ] (#1)
2. [ Task distribution ] (#2)
3. [ To-do list ] (#3)
   1. [ Initial steps: molecular dynamics code ] (#3.1)
   2. [ Paralelization of the code ] (#3.2)
4. [ How to install ] (#4)
5. [ Notes: about the use of special modules and other stuff ] (#5)


<a name="1"></a>
## 1. Team members 

* ARNAU GARCÍA DURAN
* *ALBERT ORTEGA BARTOLOMÉ* (TEAM LEADER)
* ARNAU CORTÉS LLAMAS
* MARGOT INES PACO CHIPANA

Collaborator: Elena Garcia de Lamo

 <a name="2"></a>
## 2. Reminded: task distribution

* Inicialization+PBCs --> Margot
* Forces+data visualization --> Arnau Cortés
* Integrators --> Albert
* Stadistics --> Arnau Garcia

<a name="3"></a>
## 3. To-do list 

<a name="3.1"></a>
### 3.1. Initial steps: molecular dynamics code 
- [x] *Stadistics* - Module with subroutines for calculation of macroscopic observables: kinetic energy, total momentum, instantaneous temperature, pressure, mean square displacement and radial distribution function. *module_calc_statistics.f90*.
- [x] *Forces* - Module with a subroutine to calculate the interactions and forces of all particles in the system using the Lennard-Jones potential. *module_force.f90*
- [x] *Inicialization+PBCs* - Module with subroutines for initialization of the system, simple cubic cell, and periodic boundary conditions of the box. *init_pbc.f90*.
- [x] *Integrators* - Module with different integrators: Velocity Verlet, Verlet algorithm and Euler algorithm; also, Andersen thermostat. *integrators.f90*. 
- [ ] *Data visualization* - Python program for visualization and data representation. Plots of macroscopic observables like kinetic energy, potential energy and total momentum.

<a name="3.2"></a>
### 3.2. Paralelization of the code 

<a name="4"></a>
## 4. How to install 

<a name="5"></a>
## 5. Notes 
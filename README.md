# Project_II

## Table Contents

1. [ List of team members ](#1-team)
2. [ Task distribution ](#2-task)
3. [ To-do list ](#3-list)
   1. [ Initial steps: molecular dynamics code ](#3.1-init)
   2. [ Paralelization of the code ](#3.2-para)
4. [ How to install ](#4-install)
5. [ Notes: about the use of special modules and other stuff ](#5-notes)


<a name="1-team"></a>
## 1. Team members 

* ARNAU GARCÍA DURAN
* *ALBERT ORTEGA BARTOLOMÉ* (TEAM LEADER)
* ARNAU CORTÉS LLAMAS
* MARGOT INES PACO CHIPANA

Collaborator: Elena Garcia de Lamo

 <a name="2-task"></a>
## 2. Reminded: task distribution

* Inicialization+PBCs --> Margot
* Forces+data visualization --> Arnau Cortés
* Integrators --> Albert
* Stadistics --> Arnau Garcia

<a name="3-list"></a>
## 3. To-do list 

<a name="3.1-init"></a>
### 3.1. Initial steps: molecular dynamics code 
- [x] *Stadistics* - Module with subroutines for calculation of macroscopic observables: kinetic energy, total momentum, instantaneous temperature, pressure, mean square displacement and radial distribution function. *module_calc_statistics.f90*.
- [x] *Forces* - Module with a subroutine to calculate the interactions and forces of all particles in the system using the Lennard-Jones potential. *module_force.f90*
- [x] *Inicialization+PBCs* - Module with subroutines for initialization of the system, simple cubic cell, and periodic boundary conditions of the box. *init_pbc.f90*.
- [x] *Integrators* - Module with different integrators: Velocity Verlet, Verlet algorithm and Euler algorithm; also, Andersen thermostat. *integrators.f90*. 
- [ ] *Data visualization* - Python program for visualization and data representation. Plots of macroscopic observables like kinetic energy, potential energy and total momentum.

<a name="3.2-para"></a>
### 3.2. Parallelized code 

<a name="4-install"></a>
## 4. How to install 

<a name="5-notes"></a>
## 5. Notes 

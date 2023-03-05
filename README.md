# Project_II

## Team members

* ARNAU GARCÍA DURAN
* *ALBERT ORTEGA BARTOLOMÉ* (TEAM LEADER)
* ARNAU CORTÉS LLAMAS
* MARGOT INES PACO CHIPANA

Collaborator: Elena Garcia de Lamo

## Reminded: task distribution

* Inicialization+PBCs --> Margot
* Forces+data visualization --> Arnau Cortés
* Integrators --> Albert
* Stadistics --> Arnau Garcia

## To-do list

- [x] *Stadistics* - Module with subroutines for calculation of macroscopic observables: kinetic energy, total momentum, instantaneous temperature, pressure, mean square displacement and radial distribution function. Specific directory.
- [x] *Forces* - Module with a subroutine to calculate the interactions and forces of all particles in the system using the Lennard-Jones potential. Specific directory.
- [x] *Inicialization+PBCs* Module with subroutines for initialization of the system, simple cubic cell, and periodic boundary conditions of the box. *init_pbc.f90*.
- [x] *Integrators* Module with different integrators: Velocity Verlet, Verlet algorithm and Euler algorithm; also, Andersen thermostat. *integrators.f90*. **Warning: Lennard Jones subroutine is duplicated**.
- [ ] *Data visualization* Python program for visualization and data representation. Plots of macroscopic observables like kinetic energy, ptoential energy and total momentum.


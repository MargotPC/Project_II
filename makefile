programa.x: init_pbc.o module_calc_statistics.o module_forces.o integrators.o main.o run_dir.sh
#       ./run_dir.sh
	gfortran -o programa.x init_pbc.o integrators.o module_calc_statistics.o module_forces.o main.o
main.o: main.f90
	gfortran -c -O3 main.f90
integrators.o: init_pbc.f90 integrators.f90
	gfortran -c -O3 integrators.f90
module_calc_statistics.o: module_calc_statistics.f90
	gfortran -c -O3 module_calc_statistics.f90
init_pbc.o: init_pbc.f90
	gfortran -c -O3 init_pbc.f90
module_forces.o: init_pbc.f90 module_forces.f90
	gfortran -c -O3 module_forces.f90

plot:

	python Block_average_and_plots.py

help:
	@echo ---------------------------------------------------------------------
	@echo make run: runs the program in order to obatain data
	@echo make plot: generates plots of the energy, temperature and pressure.
	@echo ---------------------------------------------------------------------

clean:
	rm -f integrators.o module_calc_statistics.o init_pbc.o module_forces.o main.o programa.x
	rm -f integrators.MOD module_calc_statistics.MOD init_pbc.MOD module_forces.MOD main.MOD
run:
	./programa.x

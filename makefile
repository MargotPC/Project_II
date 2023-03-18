FC = gfortran
FFLAGS = -O3


program.x: init_pbc.o module_calc_statistics.o module_forces.o integrators.o main.o
	./run_dir.sh
	$(FC) $(FFLAGS) -o $@ $^
main.o: main.f90
	$(FC) $(FFLAGS) -c $<

integrators.o: init_pbc.f90 integrators.f90
	$(FC) $(FFLAGS) -c $^
module_calc_statistics.o: module_calc_statistics.f90
	$(FC) $(FFLAGS) -c $^
init_pbc.o: init_pbc.f90
	$(FC) $(FFLAGS) -c $^
module_forces.o: init_pbc.f90 module_forces.f90
	$(FC) $(FFLAGS) -c $^
plot:

	python3 Block_average_and_plots.py

help:
	@echo ---------------------------------------------------------------------
	@echo make run: runs the program in order to obatain data
	@echo make plot: generates plots of the energy, temperature and pressure.
	@echo ---------------------------------------------------------------------

clean:
	rm -f *.o *.mod program.x
run:
	./program.x 

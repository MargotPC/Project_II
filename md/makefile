# Makefile
#-------------

.PHONY : main.x
main.x: output.o init_scc.o pbc.o md.o
	gfortran -o main.x $^
md.o: output.f90 md.f90
	gfortran -c md.f90
output.o: init_scc.f90 output.f90
	gfortran -c output.f90
init_scc.o: init_scc.f90
	gfortran -c init_scc.f90
pbc.o:  pbc.f90
	gfortran -c pbc.f90
	
clean:
	rm -f main.x
	rm -f output.o init_scc.o pbc.o md.o 
run:
	./main.x

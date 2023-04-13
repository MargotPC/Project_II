

.PHONY : compile


compile :
	qsub compile.sub
run:
	qsub run.sub

runAll : run1 run2 run4 run8 run16 run32 run48

run1 : 
	rm -f run1.out
	qsub -pe smp 1 -o "run1.out" run.sub

run2 : 
	rm -f run2.out
	qsub -pe smp 3 -o "run2.out" run.sub 

run4 : 
	rm -f run4.out
	qsub -pe smp 4 -o "run4.out" run.sub  

run8 : 
	rm -f run8.out
	qsub -pe smp 8 -o "run8.out" run.sub

run16 : 
	rm -f run16.out
	qsub -pe smp 16 -o "run16.out" run.sub

run32 : 
	rm -f run32.out
	qsub -pe smp 32 -o "run32.out" run.sub

run48 : 
	rm -f run48.out
	qsub -pe smp 48 -o "run48.out" run.sub

clust:
	qstat -f
	qstat -g c

SR: gettimes generatetimes

gettimes:
	@rm -f times.dat
	@tail -n 1 run1.out > 1.time
	@tail -n 1 run2.out > 2.time	
	@tail -n 1 run4.out > 4.time
	@tail -n 1 run8.out > 8.time
	@tail -n 1 run16.out > 16.time
	@tail -n 1 run32.out > 32.time
	@tail -n 1 run48.out > 48.time

generatetimes:
	@echo -e "1 \t $(shell cat 1.time)" >> times.dat
	@echo -e "2 \t $(shell cat 2.time)" >> times.dat
	@echo -e "4 \t $(shell cat 4.time)" >> times.dat
	@echo -e "8 \t $(shell cat 8.time)" >> times.dat
	@echo -e "16 \t $(shell cat 16.time)" >> times.dat
	@echo -e "32 \t $(shell cat 32.time)" >> times.dat
	@echo -e "48 \t $(shell cat 48.time)" >> times.dat  
	@rm -f *.time
	@echo 'times.dat generated succesfully'

clean:
	rm -f *.out
	rm -f *.err
	rm -f *.x
	rm -f *.o
	rm -f *.mod

help:
	@echo
	@echo -e '\t\t\t\t MAKE OPTIONS'
	@echo -e '#######################################################################################'
	@echo -e '### make compile:  Compiles the code in order to generate the necessary binary files.'
	@echo -e '###     make run:  Runs the code with the specified procs in run.sub.'
	@echo -e '###    make run1:  Runs the code with 1 proc(s).'
	@echo -e '###    make run2:  Runs the code with 2 proc(s).'
	@echo -e '###    make run4:  Runs the code with 4 proc(s).'
	@echo -e '###    make run8:  Runs the code with 8 proc(s).'
	@echo -e '###   make run16:  Runs the code with 16 proc(s).'
	@echo -e '###   make run32:  Runs the code with 32 proc(s).'
	@echo -e '###   make run48:  Runs the code with 48 proc(s).'
	@echo -e '###  make runAll:  Runs all the default codes.'
	@echo -e '###      make SR:  Generates times.dat, a file containing the elapsed times per proc.'
	@echo -e '###   make clean:  Removes all binary, error and log files.'
	@echo -e '###   make clust:  Provides information of the used cluster.'
	@echo -e '#######################################################################################'
	@echo
	
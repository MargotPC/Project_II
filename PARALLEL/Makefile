

.PHONY : compile


compile :
	qsub compile.sub
run:
	qsub run.sub 

clean:
	rm -f *out
	rm -f *err
	rm -f *.x
	rm -f *.o
	rm -f *.mod

help:
	@echo OPCIONES del Makefile
	@grep ':'  Makefile | grep -v PHONY | grep -v 'grep' | cut -f 1 -d ':'

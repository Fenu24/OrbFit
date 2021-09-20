include ../make.flags
include modules.flags

all:  occult_sim.x

LIBDEPS=  ../lib/libsuit.a ../lib/libgauss.a

# Executables

occult_sim.x:  random.o occ_mod.o occult_sim.o $(LIBDEPS)
	$(FC) $(FFLAGS) random.o occ_mod.o occult_sim.o  $(LIBDEPS) -o occult_sim.x

# tools

depend:
	../../bin/mkdep90 -I../include -L../suit -L../prelim *.f90 *.h90 > make.dep
	../../bin/mkmoddep -f ../.compiler ../suit ../prelim > modules.flags

clean:
	rm -f *.x *.o *.mod lib*.a *% *~ core 

distclean: clean
	../lib/deinst_src


include ../make.flags
include modules.flags

all: bineph.x

# Executables

BINEPH=bineph2.o  

LIBDEPS= ../lib/libprop.a ../lib/libmoid.a ../lib/libgauss.a ../lib/libsuit.a

bineph.x: $(BINEPH) $(LIBDEPS)
	$(FC) $(FFLAGS) -o bineph.x $(BINEPH) $(LIBDEPS)

depend:
	../../bin/mkdep90 -I../include -L../suit -L../propag *.f *.f90 *.h *.h90 > make.dep
	../../bin/mkmoddep -f ../.compiler ../suit ../propag > modules.flags

clean:
	rm -f *.x *.o *.mod lib*.a *% *~ core work.pc*

distclean: clean
	../lib/deinst_src

# Dependencies
include make.dep

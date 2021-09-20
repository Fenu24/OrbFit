include ../make.flags
include modules.flags

all: mov_obs_predictor.x

#packages

LIBDEPS= ../lib/libprop.a ../lib/libmoid.a ../lib/libgauss.a ../lib/libsuit.a
LIBMOP=mov_obs_predictor_util.o

# Executables
mov_obs_predictor.x: mov_obs_predictor.o $(LIBMOP) $(LIBDEPS)
	$(FC) $(FFLAGS) -o mov_obs_predictor.x  mov_obs_predictor.o $(LIBMOP) $(LIBDEPS)

depend:
	../../bin/mkdep90 -I../include -L../suit -L../propag *.f90 *.h90 > make.dep
	../../bin/mkmoddep -f ../.compiler ../suit ../propag > modules.flags

clean:
	rm -f *.x *.o *.mod *% *~ core*

# Dependencies
include make.dep

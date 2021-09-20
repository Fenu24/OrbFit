# Makefile options for HP PA RISC, optimized

# Fortran compiler
FC=f90
# Options for Fortran compiler:
FFLAGS= -O +gprof +save +real_constant=double -I../include
# "ranlib" command: if it is not needed, use "RANLIB=touch"
RANLIB=ranlib
VPATH=../include

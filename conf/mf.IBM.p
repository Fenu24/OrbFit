# Makefile options for IBM RISC 6000, profiling

# Fortran compiler
FC=xlf
# Options for Fortran compiler:
FFLAGS= -O3 -Q -qstrict -p -I../include
# "ranlib" command: if it is not needed, use "RANLIB=touch"
RANLIB=ranlib
VPATH=../include

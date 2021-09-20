# Makefile options for the INTEL compiler, optimized and profiling

# Fortran compiler
FC=ifort
# Options for Fortran compiler:
FFLAGS= -warn nousage -O -mp1 -pg -shared-intel -save -assume byterecl -I../include 
# "ranlib" command: if it is not needed, use "RANLIB=touch"
RANLIB=ranlib
VPATH=../include


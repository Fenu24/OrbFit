# Makefile options for the INTEL compiler, optimized and profiling

# Fortran compiler
FC=ifort
# Options for Fortran compiler:
FFLAGS= -warn nousage -O -fp-model precise -pg -save -assume byterecl -I../include 
# "ranlib" command: if it is not needed, use "RANLIB=touch"
RANLIB=ranlib
VPATH=../include


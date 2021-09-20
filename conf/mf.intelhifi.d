# Makefile options for the INTEL compiler, not optimized, high accuracy

# Fortran compiler
FC=ifort
# Options for Fortran compiler for debugging:
FFLAGS=  -warn nousage -g -CB -traceback -fp-model precise -save -assume byterecl -I../include 
# "ranlib" command: if it is not needed, use "RANLIB=touch"
RANLIB=ranlib
VPATH=../include

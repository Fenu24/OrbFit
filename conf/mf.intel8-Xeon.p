# Makefile options for the INTEL compiler, optimized

# Fortran compiler
FC=ifort
# Options for Fortran compiler:
FFLAGS= -warn nousage -O3 -mp1 -xN -pg -assume byterecl -Vaxlib -I../include
# "ranlib" command: if it is not needed, use "RANLIB=touch"
RANLIB=ranlib
VPATH=../include

# Makefile options for the INTEL compiler,not  optimized

# Fortran compiler
FC=ifort
# Options for Fortran compiler for debugging:
FFLAGS=  -warn nousage -g -CB -shared-intel -mcmodel=medium -traceback -fp-model fast -save -assume byterecl -I../include 
# "ranlib" command: if it is not needed, use "RANLIB=touch"
RANLIB=ranlib
VPATH=../include

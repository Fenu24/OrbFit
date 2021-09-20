# Makefile options for the INTEL compiler for MAC,not  optimized

# Fortran compiler
FC=ifort
# Options for Fortran compiler for debugging:
FFLAGS=  -warn nousage -g -CB -traceback -save -assume byterecl -Vaxlib -I../include 
# "ranlib" command: if it is not needed, use "RANLIB=touch"
RANLIB=ranlib -c
VPATH=../include

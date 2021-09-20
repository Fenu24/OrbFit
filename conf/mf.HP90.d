# Makefile options for HP PA RISC, debugging

# Fortran compiler
FC=f90
# Options for Fortran compiler:
FFLAGS= +check=all +save +real_constant=double -g -I../include
# "ranlib" command: if it is not needed, use "RANLIB=touch"
RANLIB=ranlib
VPATH=../include

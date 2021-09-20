# Makefile options for SunOS, debugging

# Fortran compiler
FC=f77
# Options for Fortran compiler:
FFLAGS=-C -g -I../include
# "ranlib" command: if it is not needed, use "RANLIB=touch"
RANLIB=ranlib
VPATH=../include

# Makefile options for IBM RISC 6000, debugging

# Fortran compiler
FC=xlf
# Options for Fortran compiler:
FFLAGS=-C -g -qextchk -I../include
# "ranlib" command: if it is not needed, use "RANLIB=touch"
RANLIB=ranlib
VPATH=../include

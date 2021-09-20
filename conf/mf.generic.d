# Makefile options for generic (unknown) UNIX workstation, debugging

# Fortran compiler
FC=f77
# Options for Fortran compiler:
FFLAGS=-g -I../include
# "ranlib" command: if it is not needed, use "RANLIB=touch"
RANLIB=ranlib
VPATH=../include

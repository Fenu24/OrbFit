# Makefile options for the GNU g95 compiler, debugging

# Fortran compiler
FC=g95
# Options for Fortran compiler for debugging:
FFLAGS=  -g -C  -I../include 
# "ranlib" command: if it is not needed, use "RANLIB=touch"
RANLIB=ranlib
VPATH=../include


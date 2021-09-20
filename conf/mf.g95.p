# Makefile options for the GNU g95 compiler, profiling

# Fortran compiler
FC=g95
# Options for Fortran compiler for debugging:
FFLAGS= -pg -O3  -I../include 
# "ranlib" command: if it is not needed, use "RANLIB=touch"
RANLIB=ranlib
VPATH=../include


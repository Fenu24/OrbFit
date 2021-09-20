# fflags for DIGITAL ALPHA
#flags for Fortran compiler, optimized and profiling
FFLAGS= -O4 -p -I../include

FC=f90

# "ranlib" command: if it is not needed, use "RANLIB=touch"
RANLIB=ranlib
VPATH=../include

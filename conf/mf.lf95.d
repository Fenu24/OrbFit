# Makefile options for Lahey fl95 compiler, debugging

# Fortran compiler
FC=lf95
# Options for Fortran compiler:
FFLAGS= --chk -g --warn  -I../include
#FFLAGS= --chk -g --warn --mod ../suit:../propag  -I../include
# "ranlib" command: if it is not needed, use "RANLIB=touch"
RANLIB=ranlib
VPATH=../include


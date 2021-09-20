include ../make.flags

all: ../lib/libsuit.a

# LIBRARY OBJECT FILES

# Observations & Statistics
OBS_OBJ=station_coordinates.o ades.o astrometric_observations.o \
     obsvarie.o iorwo_old.o attributable.o very_short_arc.o obs_simul.o\
      cat_debias.o healpix_types.o 

ATTR_OBJ=attributable.o very_short_arc.o obs_simul.o \
     delaunay_triang.o admis_reg.o triangles.o

# Reference systems and IERS time series
REF_OBJ=reference_systems.o iers_ser.o chebi_pol.o

# Orbital elements
ELEM_OBJ=orb_els.o io_elems.o ever_pitkin.o

# Mathematical routines (including EISPACK routines) 
MATH_OBJ=math_lib.o pol_zeros.o angle_var.o eigen_val.o

# JPL Digital Ephemeris and binary ephemeris, fund. constants and masses
EPH_OBJ=jpl_ephem.o read_ephem.o fund_const.o planet_masses.o

# File operations; I/O control
FILE_OBJ=file_oper.o fidinam4.o output_control.o name_rules.o

# Header and  Main option namelists; also old version 
HNAM_OBJ=header_input.o option_input.o genio7.o

# Character strings and sorting
CHA_OBJ=char_str.o bin_search.o

# Time and timescales
TIME_OBJ=time_scales.o

# Various utilities
UTIL_OBJ=util_suit.o 

SUIT_OBJ=$(FILE_OBJ) $(HNAM_OBJ) $(CHA_OBJ) $(UTIL_OBJ) $(TIME_OBJ) \
         $(ELEM_OBJ) $(MATH_OBJ) $(EPH_OBJ) $(REF_OBJ) $(OBS_OBJ) $(ATTR_OBJ)

../lib/libsuit.a: $(SUIT_OBJ)
	ar r ../lib/libsuit.a $(SUIT_OBJ)
	$(RANLIB) ../lib/libsuit.a

# macros
depend:
	../../bin/mkdep90 -I../include *.f *.f90 *.h *.h90 > make.dep
	../../bin/mkmoddep -f ../.compiler ../suit > modules.flags
clean:
	rm -f *.x *.o *.mod lib*.a *% *~ core.* *.mod *.d *.dvi 
	rm -f ../lib/libsuit.a

distclean: clean
	../lib/deinst_src

LISTS: 
	grep -E 'PROGRAM|MODULE|SUBROUTINE|FUNCTION|CONTAINS' *.f90 >LIST.sub
	grep -i ' CALL ' *.f90  >LIST.call 

# Dependencies
include make.dep

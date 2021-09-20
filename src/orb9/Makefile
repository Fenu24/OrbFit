include ../make.flags
include modules.flags

all: giffv.x conv9.x prop_class.x prop_tax.x propsynt3big.x propsyntwoway.x prop9.x  inbaric.x inbarmerc.x orbit9.x orbit9t.x propert3.x orbselect.x secre9.x purge.x merg.x mergt.x selextt.x selext.x magn.x mergsec.x splitsec.x extract_sec.x

#packages

MODULRID = char_string.o var_precision.o  fir_filters.o 

LIB9=  headio.o giffan.o coord.o angcou.o

SHORT= short9.o iosho.o lcder.o kepdel.o 

LONG=  iolon9.o long9.o secthe8.o terms9.o inte8.o selre9.o  subs9.o \
	 lly76s.o keppoinc.o 

INPUT= reapla.o coostd.o reaout.o vareqi.o reacat9.o

LIBDEPS=  ../lib/libprop.a ../lib/libmoid.a ../lib/libgauss.a ../lib/libsuit.a 

# Executables
orbit9.x:  orbit9.o outsubv.o fir_filters.o input9.o $(INPUT) reaopt.o \
                $(LIB9) $(LIBDEPS) 
	$(FC) $(FFLAGS)  orbit9.o outsubv.o fir_filters.o $(INPUT) reaopt.o \
	input9.o $(LIB9) $(LIBDEPS)   -o orbit9.x

orbit9t.x:  orbit9t.o outsubt.o input9t.o $(INPUT) reaoptt.o fir_filters.o\
                $(LIB9)  $(LIBDEPS)
	$(FC) $(FFLAGS)  orbit9t.o outsubt.o $(INPUT) reaoptt.o fir_filters.o \
	input9t.o  $(LIB9)  $(LIBDEPS)  -o orbit9t.x

conv9.x: conv9.o jobinp9.o $(LIB9) $(LIBDEPS)
	$(FC) $(FFLAGS) conv9.o jobinp9.o $(LIB9) -L../lib/ -lsuit -lprop -o conv9.x

prop_tax.x: prop_tax.o propel_mod.o $(LIBDEPS)
	$(FC)  $(FFLAGS) prop_tax.o propel_mod.o $(LIBDEPS) -o prop_tax.x

prop_class.x: prop_class.o propel_mod.o $(LIBDEPS)
	$(FC)  $(FFLAGS) prop_class.o propel_mod.o $(LIBDEPS) -o prop_class.x

propsynt3big.x: propsyntbig.o inptro.o synthcomp.o fir_filters.o \
        $(LIB9) ../lib/libsuit.a 
	$(FC) $(FFLAGS) propsyntbig.o inptro.o synthcomp.o fir_filters.o \
	$(LIB9) -L../lib/ -lsuit -lprop  -o propsynt3big.x

propsyntwoway.x: propsyntwo.o inptwoway.o synthcomp.o fir_filters.o \
	$(LIB9) $(LIBDEPS)
	$(FC) $(FFLAGS) propsyntwo.o inptwoway.o synthcomp.o fir_filters.o \
	$(LIB9) $(LIBDEPS) -o propsyntwoway.x

propert3.x: propert.o $(LIB9) synthtro.o propasub.o smallsubsn.o \
        $(LIBDEPS)
	$(FC) $(FFLAGS) propert.o $(LIB9) synthtro.o propasub.o  \
	smallsubsn.o $(LIBDEPS) -o propert3.x

giffv.x: giffv.o datinp.o remove_spaces.o plotgnu.o prompting.o\
	parameters_giff.o \
        $(LIB9) $(MODULRID) $(LIBDEPS)
	$(FC) $(FFLAGS) giffv.o datinp.o plotgnu.o prompting.o parameters_giff.o  \
	 remove_spaces.o $(LIB9) $(MODULRID) $(LIBDEPS) -o giffv.x

prop9.x: prop9.o opt9.o controlmod.o $(SHORT) $(LONG) \
	$(LIB9) $(LIBDEPS)
	$(FC) $(FFLAGS) prop9.o opt9.o controlmod.o \
	$(SHORT) $(LONG) $(LIB9) $(LIBDEPS) -o prop9.x

secre9.x: secre9.o opt9.o gridsel.o controlmod.o $(SHORT) $(LONG) \
             $(LIB9)  $(LIBDEPS)
	$(FC) $(FFLAGS) secre9.o opt9.o gridsel.o controlmod.o  \
	$(SHORT) $(LONG) $(LIB9) $(LIBDEPS)  -o secre9.x

inbaric.x: inbaric.o subdates.o epsi.o $(LIB9)  $(LIBDEPS)
	$(FC) $(FFLAGS) inbaric.o \
	subdates.o  epsi.o $(LIB9) $(LIBDEPS) -o inbaric.x

inbarmerc.x: inbarmerc.o subdates.o epsi.o $(LIB9)  $(LIBDEPS)
	$(FC) $(FFLAGS) inbarmerc.o \
	subdates.o  epsi.o $(LIB9) $(LIBDEPS) -o inbarmerc.x

#auxiliaries
#==================================================================

orbselect.x: orbselect.o  ../lib/libsuit.a 
	$(FC) $(FFLAGS) orbselect.o  -L../lib/ -lsuit -o orbselect.x

mergsec.x: mergsec.o ../lib/libsuit.a
	$(FC) $(FFLAGS) mergsec.o -L../lib/ -lsuit -o mergsec.x

splitsec.x: splitsec.o ../lib/libsuit.a
	$(FC) $(FFLAGS) splitsec.o -L../lib/ -lsuit -o splitsec.x

extract_sec.x: extract_sec.o ../lib/libsuit.a
	$(FC) $(FFLAGS) extract_sec.o -L../lib/ -lsuit -o extract_sec.x

merg.x: merg.o ../lib/libsuit.a
	$(FC) $(FFLAGS) merg.o -L../lib/ -lsuit -o merg.x

mergmult.x: mergmult.o ../lib/libsuit.a
	$(FC) $(FFLAGS) mergmult.o -L../lib/ -lsuit -o mergmult.x

mergt.x: mergt.o ../lib/libsuit.a
	$(FC) $(FFLAGS) mergt.o -L../lib/ -lsuit -o mergt.x

magn.x: magn.o ../lib/libsuit.a
	$(FC) $(FFLAGS) magn.o -L../lib/ -lsuit -o magn.x

purge.x: purge.o
	$(FC) $(FFLAGS) purge.o -o purge.x

selext.x: selext.o ../lib/libsuit.a
	$(FC) $(FFLAGS) selext.o -L../lib/ -lsuit -o selext.x

selextt.x: selextt.o ../lib/libsuit.a
	$(FC) $(FFLAGS) selextt.o -L../lib/ -lsuit -o selextt.x


depend:
	../../bin/mkdep90 -I../include -L../suit -L../propag -L../prelim *.f90 *.h90 > make.dep
	../../bin/mkmoddep -f ../.compiler ../suit ../propag ../prelim > modules.flags
clean:
	rm -f *.x *.o *.mod *% *~ core* 

distclean: clean
	../lib/deinst_src

LISTS:   LISTS
	grep -E -i  'PROGRAM|MODULE|SUBROUTINE|FUNCTION|CONTAINS' *.f90 | grep -v ! >LIST.sub


# Dependencies
include make.dep

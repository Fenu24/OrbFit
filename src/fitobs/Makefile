include ../make.flags
include modules.flags

all: fitobs.x catpro.x

#packages

COPYING=stacop.o two_states.o
#obs_convert.o

FIDENT=identsub.o

FINOUT= fitobs_mod.o   fitobs_interfaces.o 

FITOBS= fitobs.o  $(COPYING) $(FIDENT) $(FINOUT) twotes.o

LIBDEPS= ../lib/libprop.a ../lib/libmoid.a ../lib/libgauss.a ../lib/libsuit.a

# Executables

fitobs.x: $(FITOBS) $(LIBDEPS)
	$(FC) $(FFLAGS) -o fitobs.x  $(FITOBS) $(LIBDEPS)

catpro.x: catpro.o   $(LIBDEPS)
	$(FC) $(FFLAGS) catpro.o  $(LIBDEPS) -o catpro.x

depend:
	../../bin/mkdep90 -I../include -L../suit -L../propag -L../prelim  *.f90 *.h90 > make.dep
	../../bin/mkmoddep -f ../.compiler ../suit ../propag ../prelim > modules.flags

clean:
	rm -f *.x *.o *.mod *% *~ core* 

distclean: clean
	../lib/deinst_src

# Windows targets
win: ..\..\bin\fitobs.exe

..\..\bin\fitobs.exe:$(FITOBS) ..\lib\prop.lib ..\lib\moid.lib ..\lib\suit.lib ..\lib\gauss.lib
	$(FC) $(FFLAGS) /link /out:$@ $**

windist:
	@del *.o
	@del *.pdb
	@del ..\..\bin\*.pdb
	@del ..\..\bin\*.ilk

winclean:windist
	@del ..\..\bin\fitobs.*

# Dependencies
include make.dep

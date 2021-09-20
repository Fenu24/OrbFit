include ../make.flags
include modules.flags

all: orbfit.x

# Auxiliary routines (internal to orbfit)
ORBAUX=ofod.o rdopt_of.o srtobs.o seltid.o ustsp.o ofclrf.o

ORBFIT=$(ORBAUX) orbfit.o

ORBFITDEPS=$(ORBFIT) ../lib/libprop.a ../lib/libmoid.a ../lib/libgauss.a ../lib/libsuit.a

orbfit.x: $(ORBFITDEPS) 
	$(FC) $(FFLAGS) -o orbfit.x $(ORBFITDEPS)

depend:
	../../bin/mkdep90 -I../include -L../suit -L../propag  *.f *.f90 *.h *.h90 > make.dep
	../../bin/mkmoddep -f ../.compiler ../suit ../propag > modules.flags
clean:
	rm -f *.x *.o *.mod lib*.a *% *~ core work.pc*

distclean: clean
	../lib/deinst_src

# Windows targets
win: ..\..\bin\orbfit.exe

..\..\bin\orbfit.exe: $(ORBFIT) ..\lib\prop.lib ..\lib\moid.lib ..\lib\suit.lib ..\lib\gauss.lib
	$(FC) $(FFLAGS) /link /out:$@ $**

windist:
	@del *.o
	@del *.pdb
	@del ..\..\bin\*.pdb
	@del ..\..\bin\*.ilk

winclean:windist
	@del ..\..\bin\orbfit.*

# Dependencies
include make.dep

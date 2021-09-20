include ../make.flags
include modules.flags

all: ../lib/libgauss.a

# PUBLIC LIBRARY for Vaisala method
# (the only public routines are: vaisala, vaisub sel3g and alfdel)

VAIS=vaisala.o vaisub.o

LIB=iodini.o sel3g.o $(VAIS)  gaussn.o gaussn2.o gibbs.o solv8.o gauseq.o \
    gauini.o ccek1.o  anecc.o hypan.o ekcc1.o \
     adds3p.o s3dtw.o s3srt.o sel3mg.o iodsdt.o iodsbs.o random.o 
     
NEWLIB= vel_iterate.o

OLDLIB=modgeq.o

LIB90= io_det.o iod_rms.o sel3_mc.o

../lib/libgauss.a: $(LIB) $(LIB90) $(NEWLIB)
	ar r ../lib/libgauss.a $(LIB) $(LIB90) $(NEWLIB)
	$(RANLIB) ../lib/libgauss.a

depend:
	../../bin/mkdep90 -I../include -L../suit *.f *.f90 *.h *.h90 > make.dep
	../../bin/mkmoddep -f ../.compiler ../suit > modules.flags

clean:
	rm -f *.x *.o *.mod lib*.a *% *~ core
	rm -f ../lib/libgauss.a

distclean: clean
	../lib/deinst_src

# Windows targets
win: ..\lib\gauss.lib

..\lib\gauss.lib: $(LIB)
	lib /out:$@ $**

winclean:
	@del *.o
	@del *.pdb
	@del ..\lib\gauss.lib

# Dependencies
include make.dep

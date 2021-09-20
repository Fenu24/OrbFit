include ../make.flags

all: asc2eph.x testeph.x 

# Ephemerides
ephemerides: input asc2eph.x
	asc2eph.x < input.431
	mv JPLEPH jpleph
	make clean

asc2eph.x: asc2eph.o
	$(FC) $(FFLAGS) -o asc2eph.x  asc2eph.o

binshort.x: binshort.o
	$(FC) $(FFLAGS) -o binshort.x  binshort.o

testeph.x : testeph.o
	$(FC) $(FFLAGS) -o testeph.x   testeph.o

input: 
	cat header.431_572 ascp00000.431 ascp01000.431 ascp02000.431 ascp03000.431  > input.431

clean:
	rm -f  input.* *.o *.x *~

distclean:
	rm -f  input.* *.o *.x *~ *eph *EPH *.431


asc2eph.exe: 
	nmake /nologo asc2eph.o
	$(FC) $(FFLAGS) /link /out:$@ asc2eph.o


testeph.exe: 
	nmake /nologo testeph.o
	$(FC) $(FFLAGS) /link /out:$@ testeph.o


#









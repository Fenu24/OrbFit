include ../make.flags

all: binshort.x asc2eph.x testeph.x 

# Years 1970.0 to 2030.0
distephem: ephemerides binshort.x
	mv jpleph OLDEPH
	binshort.x
	mv NEWEPH jpleph
	make clean

ephemerides: input asc2eph.x
	asc2eph.x < input.405
	mv JPLEPH jpleph
	make clean

asc2eph.x: asc2eph.o
	$(FC) $(FFLAGS) -o asc2eph.x  asc2eph.o

binshort.x: binshort.o
	$(FC) $(FFLAGS) -o binshort.x  binshort.o

testeph.x : testeph.o
	$(FC) $(FFLAGS) -o testeph.x   testeph.o

input: 
	cat header.405 ascp1960.405 ascp1980.405 ascp2000.405 ascp2020.405 > input.405

clean:
	rm -f  input.* *.o *.x *~

distclean:
	rm -f  input.* *.o *.x *~ *eph *EPH *.405

#Windows Target
win: asc2eph.exe winclean

windistephem: winephem binshort.exe
	ren jpleph OLDEPH
	binshort.exe
	ren NEWEPH jpleph
	nmake /nologo winclean

winephem: winput asc2eph.exe
	asc2eph.exe < input.405
	ren JPLEPH jpleph
	nmake /nologo winclean

asc2eph.exe: 
	nmake /nologo asc2eph.o
	$(FC) $(FFLAGS) /link /out:$@ asc2eph.o

binshort.exe: 
	nmake /nologo binshort.o
	$(FC) $(FFLAGS) /link /out:$@ binshort.o

testeph.exe: 
	nmake /nologo testeph.o
	$(FC) $(FFLAGS) /link /out:$@ testeph.o

winput: 
	copy header.405+ascp1960.405+ascp1980.405+ascp2000.405+ascp2020.405 input.405

winclean:
	del input.*
	del *.o

windist: win
	del *eph
	del *EPH
	del *.405

#









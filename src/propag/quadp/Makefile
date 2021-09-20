.f.o:
	$(FC) $(FFLAGS) -c $*.f -o $*.o 

include ../../make.flags

libquamia.a: dqags.o dqagse.o dqagp.o dqagpe.o dqelg.o dqk21.o dqpsrt.o \
	dqagsc.o dqagsec.o dqagpc.o dqagpec.o dqelgc.o dqk21c.o dqpsrtc.o 
	ar -r libquamia.a  dqags.o dqagse.o dqagp.o dqagpe.o dqelg.o dqk21.o \
	dqpsrt.o dqagsc.o dqagsec.o dqagpc.o dqagpec.o dqelgc.o dqk21c.o \
	dqpsrtc.o

clean:
	-(rm -f *.o *.x core *.stb *~ *.out *.a work.pc*)

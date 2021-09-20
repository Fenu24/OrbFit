include ../make.flags
include modules.flags

all: ../lib/libprop.a libmoid.a test_moid_rms.x

PROP= propag_state.o force_model.o rmodel.o  yark_pert.o  non_grav.o \
        runge_kutta_gauss.o ra15_mod.o cobweb.o yorp_module.o

SAT= force_sat.o spher_harm.o perturbations.o

ORB9= force9d.o massmod.o

CLOAP=  close_app.o tp_trace.o surf_trace.o semi_linear.o cla_store.o\
        ret_analysistp.o offlov_checktp.o eval_risk.o obssto.o \
        virtual_impactor.o fsteph_earth.o

PUB= ephem_prop.o plocbd.o fitsubs.o pred_obs.o

DIFCOR= multiple_sol.o least_squares.o obs_correl.o rms_twobody.o count_opp.o count_opp2.o multi_store.o

ORBEL=orbit_elements.o gaussdeg8.o arc_control.o laplace_poincare.o dyn_param.o 

QUADP = quadp/dqags.o quadp/dqagse.o quadp/dqagp.o quadp/dqagpe.o \
        quadp/dqelg.o quadp/dqk21.o quadp/dqpsrt.o quadp/dqagsc.o \
        quadp/dqagsec.o quadp/dqagpc.o quadp/dqagpec.o quadp/dqelgc.o \
        quadp/dqk21c.o quadp/dqpsrtc.o 

CDET= detcomp.o

#STATCOMP= critical_points.o solvpoly.o rvfft.o compute_minima_ta.o\
#        decode_out.o code_input.o \
#        irvfft.o choosedeg.o moid_compute.o
STATCOMP= critical_points.o solvpoly.o moid_utilities.o moid_compute.o

MOIDRMS= dnod_amoid_rms_com.o

DES=des_routines5.o 

libquamia:
	cd quadp; make

../lib/libprop.a: $(PROP) $(DIFCOR) $(PUB) $(CLOAP) $(ORBEL) $(SAT) $(ORB9) $(DES) libquamia
	ar r ../lib/libprop.a $(PROP) $(DIFCOR) $(PUB) $(CLOAP) $(QUADP) $(ORBEL) $(SAT) $(ORB9) $(DES)
	$(RANLIB) ../lib/libprop.a

libmoid.a: $(CDET) $(STATCOMP) $(MOIDRMS)
	ar r ../lib/libmoid.a $(CDET) $(STATCOMP) $(MOIDRMS)
	$(RANLIB) ../lib/libmoid.a

LIBDEPS= ../lib/libprop.a ../lib/libmoid.a ../lib/libgauss.a ../lib/libsuit.a

# programs
test_moid_rms.x: test_moid_rms.o  $(LIBDEPS)
	$(FC) $(FFLAGS) test_moid_rms.o  $(LIBDEPS) -o test_moid_rms.x

depend:
	../../bin/mkdep90 -I../include -L../suit *.f *.f90 *.h *.h90 > make.dep
	../../bin/mkmoddep -f ../.compiler ../suit > modules.flags

clean:
	rm -f *.x *.o *.mod lib*.a *% *~ core work.pc*
	rm -f ../lib/libprop*.a
	cd quadp ; make clean  

distclean: clean
	../lib/deinst_src

LISTS: 
	grep -E 'PROGRAM|MODULE|SUBROUTINE|FUNCTION|CONTAINS' *.f90 >LIST.sub
	grep -i ' CALL ' *.f90  >LIST.call 
#	grep -B1 -E 'comment=' *.f90 >LIST.option

# Windows targets
win: ..\lib\prop.lib 

..\lib\prop.lib:$(PROP) $(DIFCOR) $(PUB) $(CLOAP)
	lib /out:$@ $**

winclean:
	@del *.o
	@del *.pdb
	@del ..\lib\prop.lib

# Dependencies
include make.dep

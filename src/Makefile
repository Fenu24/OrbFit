all:
	(cd suit;   touch modules.flags; make depend; make)
	(cd prelim; touch modules.flags; make depend; make)
	(cd propag; touch modules.flags; make depend; make)
	(cd bineph; touch modules.flags; make depend; make)
	(cd fitobs; touch modules.flags; make depend; make)
	(cd jpleph; make)
	(cd orbfit; touch modules.flags; make depend; make)
	(cd orb9; make depend; make)

tests_src:
	(cd tests;  touch modules.flags; make depend; make)

nondistribution:
	-(cd newret; make depend; make)
	-(cd allobs; make depend; make)
	-(cd astdys; make depend; make)
	-(cd neodys; make depend; make)
	-(cd panst; make depend; make)
	-(cd kepint; make depend; make) 
	-(cd link3att; make depend; make)
#	-(cd triang; make depend; make) frozen
	-(cd test; make depend; make) 
	-(cd util90; make)
	-(cd propneo2; make depend; make)
	-(cd sec_evol; make depend; make)
	-(cd moid_research; make depend; make)
	-(cd idconf2; make depend; make)
	-(cd imp_corridor; make depend; make)
	-(cd mov_obs_predictor; make depend; make)
#	-(cd stat; make depend; make)
#	-(cd skymap; make depend; make)
#	-(cd starcat; make)

nondistclean:
	-(cd newret; make clean)
	-(cd allobs; make clean)
	-(cd astdys; make clean)
	-(cd neodys; make clean)
	-(cd panst; make clean)
	-(cd kepint; make clean)
	-(cd link3att; make clean)
#	-(cd triang; make clean)
	-(cd test; make clean)
	-(cd util90; make clean)
	-(cd propneo2; make clean)
	-(cd sec_evol; make clean)
	-(cd potential_NEAs; make clean)
	-(cd moid_research; make clean)
	-(cd idconf2; make clean)
	-(cd imp_corridor; make clean)
	-(cd mov_obs_predictor; make clean)
#	-(cd stat; make clean)
#	-(cd skymap; make clean)
#	-(cd starcat; make clean)


clean:
	(cd suit; make clean)
	(cd prelim; make clean)
	(cd propag; make clean)
	(cd bineph; make clean)
	(cd fitobs; make clean)
	(cd orbfit; make clean)
	(cd jpleph; make clean)
	(cd orb9; make clean)
	rm -f *~

distclean:
	touch `cat ../conf/list.modfla`
	cp -p ../conf/doclib.h90 include/
	cp -p ../conf/parlib.h90 include/
	(cd suit; make distclean)
	(cd prelim; make distclean)
	(cd propag; make distclean)
	(cd bineph; make distclean)
	(cd fitobs; make distclean)
	(cd orbfit; make distclean)
	(cd jpleph; make distclean)
	(cd orb9; make distclean)
	rm -f lib/lib*.a
	cp ../conf/mf.dist make.flags
	rm -f *~
	rm -f .compiler


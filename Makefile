all: .conf.status
	(cd src; make)

tests: test

test:
	(cd tests; make)

clean:
	(cd src; make clean)
	(cd lib; make clean)

del_obsolete_files:
	rm -f `cat dellist`

distclean:
	rm -f *~
	(cd src; make distclean)
	(cd lib; make distclean)
	(cd tests; make distclean)
	-(cd src/lib; rm *.a)
	(rm -f .conf.status)

.conf.status:
	@ echo "Please run \"config\" before doing \"make\""
	@ exit 1

distribution: distclean notar; \
	tar -cf ../OrbFit5.0.7.tar -X notar . ; \
	gzip ../OrbFit5.0.7.tar

doctar:
	tar -cf ../doc.tar --exclude-from notar ./doc; gzip ../doc.tar

#additional_doc:
#	tar -cvf ../additional_doc.tar ./doc/additional_doc; gzip ../additional_doc.tar 

nondistribute: 
	cd src; make nondistclean; cd .. 
	tar -T notar -cf ../OrbFitwork50.tar ; gzip ../OrbFitwork50.tar

panst: 
	cd src/panst; make clean; cd ../../tests/panst; make distclean
	tar -czf ../panst5.0.tgz src/panst tests/panst lib/orbsrv.key

#debris:
#	cd src/debris; make clean; cd ../../tests/debris; make distclean
#	tar -czf ../debris4.2.tgz src/debris tests/debris src/kepint tests/kepint

#starcat:
#	cd ../starcat; make clean; cd ../skymap; make clean; cd ../..;
#	tar -cvzf ../starcat.tgz src/starcat src/skymap 

nondistclean:
	cd src; make nondistclean; cd ../tests; make nondistclean

patch: 
	tar -T patchlist -cvf ../patch5.0.3.tar ; gzip ../patch5.0.3.tar




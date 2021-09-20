all: fitobs orbfit bineph
	(cd bineph; make)
	(cd fitobs; make)
#	(cd orbfit; make)
	(cd moid; make)
	(cd orbit9; make)

clean: distclean

distclean:
	(cd bineph; make distclean)
#	(cd orbfit; make distclean)
	(cd fitobs; make distclean)
	(cd moid; make distclean)  
	(cd orbit9; make distclean)
	(cd orbit9_yorp; make distclean)

nondistclean:
	-(cd neodys; make distclean)
	-(cd panst; make distclean)
	-(cd lostPHAs; make distclean)

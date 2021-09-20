
all: testo9 testconv testprop testsynt 

testbaric:
	cp -p testout/baric/inbaric.opt .
	time ./inbaric.x
	time ./inbarmerc.x
	-(diff planxx.inc testout/baric)
	-(diff pvenexx.inc testout/baric)
	-(diff barsunxx.inc testout/baric)
	-(diff barmerxx.inc testout/baric)

testo9:
	-(make testbaric)
	cp -p testout/integr/orb9.opt .
	time ./orbit9.x
	-(diff vast.dat testout/integr/vast.dat)
	-(diff vpla.dat testout/integr/vpla.dat)
	-(diff vast.fil testout/integr/vast.fil)
	-(diff vpla.fil testout/integr/vpla.fil)
	-(diff orb9.out testout/integr/orb9.out)

testconv:
	cp -p testout/conversion/conv_test .
	cp -p testout/conversion/conv9.inp .
	time ./conv9.x < conv_test
	-(diff v1.fil testout/conversion/v1.fil)
	-(diff v3.fil testout/conversion/v3.fil)

testprop:
	cp -p testout/analytic/prop.opt .
	cp -p testout/analytic/gener.lng .
	cp -p testout/analytic/hamil.dir .
	cp -p testout/analytic/hamil.ind .
	cp -p testout/analytic/short8s.opt .
	cp -p testout/analytic/long8s.opt .
	time ./prop9.x
	-(diff v9.mea testout/analytic/v9.mea)
	-(diff v9.rem testout/analytic/v9.rem)
	-(diff v9.pro testout/analytic/v9.pro)
	-(diff v9.ang testout/analytic/v9.ang)
	-(diff v9.dis testout/analytic/v9.dis)
	-(diff v9.adi testout/analytic/v9.adi)
	-(diff v9.res testout/analytic/v9.res)


testsynt:
	rm -f vast.dat vpla.dat vast.fil vpla.fil
	cp -p testout/synthetic/vast.fil .
	cp -p testout/synthetic/vpla.fil .
	cp -p testout/synthetic/prop_test .
	time ./propsynt3big.x < prop_test
	-(diff vprs.ang testout/synthetic/vprs.ang)
	-(diff vprs.del testout/synthetic/vprs.del)
	-(diff vprs.err testout/synthetic/vprs.err)
	-(diff vprs.out testout/synthetic/vprs.out)
	-(diff vprs.pro testout/synthetic/vprs.pro)
	-(diff vprs.sig testout/synthetic/vprs.sig)



clean:
	-(rm -f orb9.dma orb9.dmp orb9.num orb9.out orb9.out orb9.opt orb9.clo vast.fil vast.dat vpla.fil vpla.dat)
	-(rm -f prop_test vprs.ang vprs.del vprs.err vprs.out vprs.pro vprs.sig)
	-(rm -f inbaric.opt allplaxx.inc allplmxx.inc barsunxx.inc barmerxx.inc planxx.inc pvenexx.inc epsi.out)
	-(rm -f prop.opt prop9.done gener.lng hamil.dir hamil.ind long8s.opt short8s.opt v9.mea v9.rem v9.pro v9.ang v9.dis v9.adi v9.res) 
	-(rm -f conv_test conv9.inp v?.fil)

distclean: clean

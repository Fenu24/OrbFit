all: all_noiter

all_noiter: CPV_noiter.bep

all_iter: CPV_iter.bep

all_17: AST17.bep

distclean: clean	
clean: 
	-(rm -f *.bep *.bai core gmon.out *~ *.*~ *.done *.log *.err CPV.bop CPV.inc CPV.clo)
	-(rm -f AST17.bop AST17.inc)

CPV_noiter.bep:
	-(cp -p testout/CPV_noiter.bop CPV.bop)
	-(cp -p testout/CPV.inc CPV.inc)
	echo CPV | time ./bineph.x
	-(diff CPV.bai testout/CPV0.bai)
	-(cp -p testout/CPV0.bai testout/CPV0.bai_noiter)
	-(cp -p testout/CPV0.bep testout/CPV0.bep_noiter)
	-(cp -p CPV.bai testout/CPV0.bai)
	-(cp -p CPV.bep testout/CPV0.bep)

CPV_iter.bep :
	-(cp -p testout/CPV_iter.bop CPV.bop)
	-(cp -p testout/CPV.inc CPV.inc)
	echo CPV | ./bineph.x
	-(diff CPV.bai CPV0.bai)

AST17.bep :
	-(cp -p testout/AST17.bop AST17.bop)
	-(cp -p testout/AST17.inc AST17.inc)
	-(cp -p testout/CPV.bai_431_fcct CPV.bai)
	-(cp -p testout/CPV.bep_431_fcct CPV.bep)
	echo AST17 | ./bineph.x

../../src/bineph/bineph.x:
	(cd ../..; make)

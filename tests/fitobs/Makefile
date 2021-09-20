all: test_fitobs

test_fitobs : test_prelim test_6dim test_drpa2m test_apriori test_yarko_fcct14 #test_IP_sampl

test_prelim:
	-(cp fcct14/1862red.fop .)
	-(cp fcct14/1862red.inp .)
	-(cp fcct14/1862red.rwo mpcobs)
	-(rm -f AST17.bai AST17.bep)
	-(ln -s ../bineph/testout/AST17.bai_431_fcct AST17.bai)
	-(ln -s ../bineph/testout/AST17.bep_431_fcct AST17.bep)
	time ./fitobs.x < 1862red.inp

test_6dim :
	-(cp fcct14/1862.fop .)
	-(cp fcct14/1862.inp .)
	-(cp fcct14/1862.rwo mpcobs/)
	-(cp fcct14/1862.eq0 epoch/)
	-(rm -f AST17.bai AST17.bep)
	-(ln -s ../bineph/testout/AST17.bai_431_fcct AST17.bai)
	-(ln -s ../bineph/testout/AST17.bep_431_fcct AST17.bep)
	time ./fitobs.x < 1862.inp
	-(diff 1862.fel testout/1862.fel || echo 'Files differ')

test_drpa2m :
	-(cp 2009BD/2009BD.fop .)
	-(cp 2009BD/2009BD.inp .)
	-(cp 2009BD/2009BD.rwo mpcobs/)
	-(cp 2009BD/2009BD.eq0 epoch/)
	-(rm -f AST17.bai AST17.bep)
	-(ln -s ../bineph/testout/AST17.bai_431_fcct AST17.bai)
	-(ln -s ../bineph/testout/AST17.bep_431_fcct AST17.bep)
	time ./fitobs.x < 2009BD.inp
	-(diff 2009BD.fel testout/2009BD.fel || echo 'Files differ')

test_apriori:
	-(cp 2009FD/410777.fop .)
	-(cp 2009FD/410777.inp .)
	-(cp 2009FD/410777.rwo mpcobs/)
	-(cp 2009FD/410777.eq0 epoch/)
	-(rm -f AST17.bai AST17.bep)
	-(ln -s ../bineph/testout/AST17.bai_431_fcct AST17.bai)
	-(ln -s ../bineph/testout/AST17.bep_431_fcct AST17.bep)
	time ./fitobs.x < 410777.inp
	-(diff 410777.fel testout/410777.fel || echo 'Files differ')

test_yarko_fcct14:
	-(cp fcct14/1862.fop_yarko 1862.fop)
	-(cp fcct14/1862.inp .)
	-(cp fcct14/1862.rwo mpcobs/)
	-(cp fcct14/1862.eq0 epoch/)
	-(rm -f AST17.bai AST17.bep)
	-(ln -s ../bineph/testout/AST17.bai_431_fcct AST17.bai)
	-(ln -s ../bineph/testout/AST17.bep_431_fcct AST17.bep)
	time ./fitobs.x < 1862.inp
	-(diff 1862.fel testout/1862.fel_yarko_fcct || echo 'Files differ')

test_IP_sampl:
	-(cp 2015RN35/2015RN35.fop .)
	-(cp 2015RN35/2015RN35.inp .)
	-(cp 2015RN35/2015RN35.rwo mpcobs/)
	-(cp 2015RN35/2015RN35.eq0 epoch/)
	-(rm -f AST17.bai AST17.bep)
	-(ln -s ../bineph/testout/AST17.bai_431_fcct AST17.bai)
	-(ln -s ../bineph/testout/AST17.bep_431_fcct AST17.bep)
	time ./fitobs.x < 2015RN35.inp
	-(diff 2015RN35.cat testout/2015RN35.cat   || echo 'Files differ')
	-(diff 2015RN35.ctc testout/2015RN35.ctc   || echo 'Files differ')
	-(diff 2015RN35.mrep testout/2015RN35.mrep || echo 'Files differ')

test_yarko_cbm10:
	-(cp cbm10/1862.fop .)
	-(cp cbm10/1862.inp .)
	-(cp cbm10/1862.rwo mpcobs/)
	-(cp cbm10/1862.eq0 epoch/)
	-(ln -s ~/options_static/AST17.bai AST17.bai)
	-(ln -s ~/options_static/AST17.bep AST17.bep)
	time ./fitobs.x < 1862.inp
	-(diff 1862.fel testout/1862.fel_yarko_cbm || echo 'Files differ')

clean:
	-(rm -f gmon.out core *.pro *.clo *.fga *.rwo *.fou *.fel *.cbd *.eph)
	-(rm -f giffv.* *.frw *.*~ *~ *.tmp *.err *.done multprop.ctc *.mou)
	-(rm -f virimp.* *.vis mult.ctc *.rms *.out *.new *.risk *.warn *.rep)
	-(rm -f *.mrep *.cat *.ctc)
	-(rm -f *.bai *.bep *.inp *.fop)
	-(rm -f mpcobs/*.rwo)
	-(rm -f epoch/*.eq0)

distclean : clean


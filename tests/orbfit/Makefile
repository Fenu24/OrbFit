all: yarko2 yarko1 simple_eq0 simple_noeq0 ident

simple_noeq0 : orbfit.x 2015AB.fel

simple_eq0 : orbfit.x 410777.fel

yarko1: orbfit.x 410777_yarko1.fel

yarko2: orbfit.x 410777_yarko2.fel

ident : orbfit.x 2015TF202.fel

2015AB.fel : 	orbfit.x
		-(ln -s ../bineph/testout/AST17.bai_431_fcct AST17.bai)
		-(ln -s ../bineph/testout/AST17.bep_431_fcct AST17.bep)
		-(cp -p 2015AB/2015AB.obs mpcobs/)
		-(cp -p 2015AB/2015AB.oop .)
		-(cp -p 2015AB/2015AB.inp .)
		time orbfit.x < 2015AB.inp
		-(diff 2015AB.oel testout/2015AB.oel || echo 'Files differ')

410777.fel :  	orbfit.x
		-(ln -s ../bineph/testout/AST17.bai_431_fcct AST17.bai)
		-(ln -s ../bineph/testout/AST17.bep_431_fcct AST17.bep)
		-(cp -p 410777/410777.rwo mpcobs/)
		-(cp -p 410777/410777.eq0 epoch/)
		-(cp -p 410777/410777.oop .)
		-(cp -p 410777/410777.inp .)
		time orbfit.x < 410777.inp
		-(diff 410777.oel testout/410777.oel || echo 'Files differ')

410777_yarko1.fel : orbfit.x
		-(ln -s ../bineph/testout/AST17.bai_431_fcct AST17.bai)
		-(ln -s ../bineph/testout/AST17.bep_431_fcct AST17.bep)
		-(cp -p 410777_yarko/410777.rwo mpcobs/)
		-(cp -p 410777_yarko/410777.eq0 epoch/)
		-(cp -p 410777_yarko/410777_1.oop 410777.oop)
		-(cp -p 410777_yarko/410777.inp .)
		time orbfit.x < 410777.inp
#		-(diff 410777.oel testout/410777.oel || echo 'Files differ')

410777_yarko2.fel : orbfit.x
		-(ln -s ../bineph/testout/AST17.bai_431_fcct AST17.bai)
		-(ln -s ../bineph/testout/AST17.bep_431_fcct AST17.bep)
		-(cp -p 410777_yarko/410777.rwo mpcobs/)
		-(cp -p 410777_yarko/410777.eq0.yark epoch/410777.eq0)
		-(cp -p 410777_yarko/410777_2.oop 410777.oop)
		-(cp -p 410777_yarko/410777.inp .)
		time orbfit.x < 410777.inp
#		-(diff 410777.oel testout/410777.oel || echo 'Files differ')

2015TF202.fel : orbfit.x
		-(ln -s ../bineph/testout/AST17.bai_431_fcct AST17.bai)
		-(ln -s ../bineph/testout/AST17.bep_431_fcct AST17.bep)
		-(cp -p identification/2015TF202_1.obs mpcobs/)
		-(cp -p identification/2015TF202_2.obs mpcobs/)
		-(cp -p identification/2015TF202.oop .)
		-(cp -p identification/2015TF202.inp .)
		time orbfit.x < 2015TF202.inp
		-(diff 2015TF202.oel testout/ident/2015TF202.oel || echo 'Files differ')
		-(diff 2015TF202_1=2015TF202_2.rwo testout/ident/2015TF202_1=2015TF202_2.rwo || echo 'Files differ')


distclean: clean

clean:
	-(rm -f *.pro *.clo *.odc *.oel *.olg *.oep *.rwo *~ reject.rep *.err)
	-(rm -f *.eq0 *.oop *.obs *.bat *.rad *.ele *.inp)
	-(rm -f *.bai *bep)
	-(rm -f mpcobs/2015TF202_*.rwo mpcobs/2015AB.rwo mpcobs/410777.rwo)




all: test_moidrms

test_moidrms : test_moid_rms.x
	time ./test_moid_rms.x < test_moid_rms.inp

distclean: clean

clean:
	(rm -f *.fail *.err *.mou gmon.out *.log *~ \#*\# error_file )

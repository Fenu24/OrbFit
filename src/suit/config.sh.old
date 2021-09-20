#!	/bin/sh

libdir=`( cd ../../lib ; pwd )`
echo "default library directory is $libdir"
cat > parlib.h90 << END
! Default library directory
CHARACTER(LEN=100), PARAMETER :: dlibd = '${libdir}'
END

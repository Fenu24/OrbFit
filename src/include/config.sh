#!	/bin/sh

docdir=`( cd ../../doc ; pwd )`
echo "default documentation directory is $docdir"
cat > doclib.h90 << END
! Copyright (C) 1998 by Andrea Milani (milani@dm.unipi.it)
! Version: August 4, 1998
! ---------------------------------------------------------------------
! Default documentation directory
CHARACTER*100 ddocd
PARAMETER (ddocd=  '${docdir}')
END

libdir=`( cd ../../lib ; pwd )`
echo "default library directory is $libdir"
cat > parlib.h90 << END
! Default library directory
CHARACTER(LEN=100), PARAMETER :: dlibd = '${libdir}'
END

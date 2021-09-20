! program to add magnitudes into the records
PROGRAM magn
  IMPLICIT NONE
  INTEGER, PARAMETER :: nax=500000
  CHARACTER*1 a1,comcha
  CHARACTER*9 name1,name2,name(nax)
  CHARACTER*10 tmp
  DOUBLE PRECISION a(7),b(7),c(8),c8(nax)
  INTEGER i,lr,ipc,l,int,imax
  !========================================
  ! open the file for output
  OPEN(9,file='numb.syn.mag',status='unknown')
  ! open the file with proper elements
  OPEN(7,file='numb.syn.nomag',status='old')
! open file with magnitudes
  OPEN(6,file='allnum.cat',status='old')
  !        open(6,file='ufitobs.cat',status='old')
  ! skip the header
  DO i=1,6
     READ(6,400)a1
  ENDDO
  ! find the magnitude
  imax=1
10 READ(6,500,end=13) a1,tmp,c
  ! parse the name to get rid of the trailing "'"
  comcha="'"
  lr=LEN(tmp)
  ipc=INDEX(tmp(1:lr),comcha)
  IF(ipc.EQ.0) THEN
     WRITE(*,*)'error reading allnum.cat; asteroid ', tmp
     STOP
  ELSE
     name2=tmp(1:ipc-1)
     lr=LEN(name2)
     IF(lr.LT.1) THEN
        WRITE(*,*)'error parsing name; asteroid ', name2
        STOP
     ENDIF
  ENDIF
  name(imax)=name2
  c8(imax)=c(8)
  imax=imax+1
  GOTO 10
13 CONTINUE
  ! read the proper elements file
14 READ(7,100,end=99) name1,b,int
  CALL rmsp(name1,l)
  ! for tno-numbered
  !          int=100
  ! for tno-multiopposition and secular resonant
  int=10
  ! for the outer main belt
  !          int=2 
  ! write the line of output                                              
  DO i=1,imax-1
     IF(name1.eq.name(i))THEN
        WRITE(9,101)name1,c8(i),b,int
        GOTO 14
     ENDIF
  ENDDO
100 FORMAT(a9,f12.7,f11.7,f10.7,f12.6,f13.6,f13.6,f8.2,i4)
101 FORMAT(a9,f6.2,f12.7,f11.7,f10.7,f12.6,f13.6,f13.6,f8.2,i4)
400 FORMAT(a1)
500 FORMAT(a1,a10,f16.6,6e25.16e3,f6.2)
99 CLOSE (7)
  CLOSE (6)
  CLOSE (9)
END PROGRAM magn

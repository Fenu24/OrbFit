! program to merge outputs from basic and extended integrations
! and to add magnitudes into the records; WARNING: alternate 
! input and output files for numbered and multiopposition asteroids
PROGRAM mergt
  IMPLICIT NONE
  CHARACTER*1 a1,comcha
  CHARACTER*9 name,name1,name2
  CHARACTER*10 tmp
  DOUBLE PRECISION a(8),b(8),c(8)
  INTEGER int,i,lr,ipc,l
  ! open the file with magnitudes
  !    OPEN(6,file='allnum.cat',status='old')
      OPEN(6,file='ufitobs.cat',status='old')
  ! skip the header
  DO i=1,6
     READ(6,400)a1
  ENDDO
400 FORMAT(a1)
  ! open files with proper elements
  !      OPEN(7,file='tro_5.pro',status='old')
  !      OPEN(8,file='tro_50.pro',status='old')
   OPEN(7,file='trom_5.pro',status='old')
   OPEN(8,file='trom_50.pro',status='old')
  !      OPEN(9,file='numbtro.pro',status='unknown')
   OPEN(9,file='multtro.pro',status='unknown')
  ! beginning of the read loop
1 READ(8,100,end=2) name,a
  CALL rmsp(name,l)
2 READ(7,100,end=99) name1,b
  CALL rmsp(name1,l)
  IF(name.eq.name1)THEN
     int=50
     ! find the magnitude
13   READ(6,500,end=99) a1,tmp,c
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
      IF(name.eq.name2) THEN
        WRITE(9,101)name,c(8),a,int
     ELSE
         GOTO 13
     ENDIF
  ELSE
     int=5
     ! find the magnitude
14   READ(6,500,end=99) a1,tmp,c
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
     IF(name1.eq.name2) THEN
        WRITE(9,101)name1,c(8),b,int
     ELSE
        GOTO 14
     ENDIF
     GOTO 2
  ENDIF
  GOTO 1
! format for numbered
!100 FORMAT(1x,a9,f6.4,f7.2,f7.3,f6.4,f8.2,f6.4,f8.2,2i3)
!101 FORMAT(a9,f6.2,2x,f6.4,f7.2,f7.3,2x,f6.4,f8.2,2x,f6.4,f8.2,2i4)
!format for multiopposition
100 FORMAT(1x,a9,f6.4,f7.2,f7.3,f6.4,f8.2,f6.4,f8.2,2i3)
101 FORMAT(a9,f6.2,2x,f6.4,f7.2,f7.3,2x,f6.4,f8.2,2x,f6.4,f8.2,2i4)
99 CLOSE (7)
  CLOSE (8)
  CLOSE (9)
!         OPEN(7,file='tro_5.sig',status='old')
!         OPEN(8,file='tro_50.sig',status='old')
  OPEN(7,file='trom_5.sig',status='old')
  OPEN(8,file='trom_50.sig',status='old')
!         OPEN(9,file='numbtro.sig',status='unknown')
  OPEN(9,file='multtro.sig',status='unknown')
3 READ(8,200,end=4) name,a
  IF(a(8).eq.-999.99d0)THEN
     a(8)=999.99
  ELSEIF(a(8).lt.0.d0)THEN
     a(8)=0.d0
  ENDIF
4 READ(7,200,end=98) name1,b
  IF(b(8).eq.-999.99d0)THEN
     b(8)=999.99
  ELSEIF(b(8).lt.0.d0)THEN
     b(8)=0.d0
  ENDIF
  IF(name.eq.name1)THEN
     int=50
     WRITE(9,201)name,a,int
  ELSE
     int=5
     WRITE(9,201)name1,b,int
     GOTO 4
  ENDIF
  GOTO 3
!format for numbered
!200 FORMAT(1x,a9,f7.5,f8.3,f8.4,f7.5,f9.3,f7.5,f9.3,f6.2)
!201 FORMAT(1x,a9,2x,f7.5,f8.3,f8.4,2x,f7.5,f9.3,2x,f7.5,f9.3,f6.2,i4)
!format for multiopposition
200 FORMAT(a10,f7.5,f8.3,f8.4,f7.5,f9.3,f7.5,f9.3,f6.2)
201 FORMAT(a10,2x,f7.5,f8.3,f8.4,2x,f7.5,f9.3,2x,f7.5,f9.3,f6.2,i4)
98 CLOSE (7)
  CLOSE (8)
  CLOSE (9)
  !        OPEN(7,file='tro_5.del',status='old')
  !        OPEN(8,file='tro_50.del',status='old')
  OPEN(7,file='trom_5.del',status='old')
  OPEN(8,file='trom_50.del',status='old')
  !        OPEN(9,file='numbtro.del',status='unknown')
  OPEN(9,file='multtro.del',status='unknown')
5 READ(8,300,end=6) name,a
  IF(a(8).eq.-999.99d0)THEN
     a(8)=999.99
  ELSEIF(a(8).lt.0.d0)THEN
     a(8)=0.d0
  ENDIF
6 READ(7,300,end=97) name1,b
  IF(b(8).eq.-999.99d0)THEN
     b(8)=999.99
  ELSEIF(b(8).lt.0.d0)THEN
     b(8)=0.d0
  ENDIF
  IF(name.eq.name1)THEN
     int=50
     WRITE(9,301)name,a,int
  ELSE
     int=5
     WRITE(9,301)name1,b,int
     GOTO 6
  ENDIF
  GOTO 5
!format for numbered
!300 FORMAT(1x,a9,f7.5,f8.3,f8.4,f7.5,f9.3,f7.5,f9.3,f6.2)
!301 FORMAT(1x,a9,2x,f7.5,f8.3,f8.4,2x,f7.5,f9.3,2x,f7.5,f9.3,f6.2,i4)
!format for multiopposition
300 FORMAT(a10,f7.5,f8.3,f8.4,f7.5,f9.3,f7.5,f9.3,f6.2)
301 FORMAT(a10,2x,f7.5,f8.3,f8.4,2x,f7.5,f9.3,2x,f7.5,f9.3,f6.2,i4)
  ! format allnum.cat
500 FORMAT(a1,a10,f16.6,6e25.16,f6.2)
97 CONTINUE
END PROGRAM mergt

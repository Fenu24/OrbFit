! program to merge outputs from basic and extended integrations         
! and to add magnitudes into the records                                
PROGRAM merg
  IMPLICIT NONE
  CHARACTER*1 a1,comcha
  CHARACTER*10 name,name1,name2,tmp
  DOUBLE PRECISION a(7),b(7),c(8)
  INTEGER int,i,lr,ipc,l
  ! skip the header of allnum.cat
  OPEN(6,file='allnum.cat',status='old')
  DO i=1,6
     READ(6,400)a1
  ENDDO
400 FORMAT(a1)
  ! open files with proper elements
  OPEN(7,file='numb-2.pro',status='old')
  OPEN(8,file='numb10.pro',status='old')
  OPEN(9,file='numball.pro',status='unknown')
  ! beginning of the read loop
1 READ(8,100,end=2) name,a
  CALL rmsp(name,l)
2 READ(7,100,end=99) name1,b
  CALL rmsp(name1,l)
  IF(name.eq.name1)THEN
     int=10
     IF(a(7).lt.0.d0)THEN
        a(7)=0.d0
     ENDIF
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
     IF(name.eq.name2)THEN
        WRITE(9,101)name,c(8),a,int
     ELSE
        GOTO 13
     ENDIF
  ELSE
     int=-2
     b(7)=-b(7)
     IF(b(7).eq.-999.99d0)THEN
        b(7)=999.99
     ELSEIF(b(7).lt.0.d0)THEN
        b(7)=0.d0
     ENDIF
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
     IF(name1.eq.name2)THEN
        WRITE(9,101)name1,c(8),b,int
     ELSE
        GOTO 14
     ENDIF
     GOTO 2
  ENDIF
  GOTO 1
100 FORMAT(a9,f12.7,f11.7,f10.7,f12.6,f13.6,f13.6,f8.2)
101 FORMAT(a9,f6.2,f12.7,f11.7,f10.7,f12.6,f13.6,f13.6,f8.2,i4)
99 CLOSE (7)
  CLOSE (8)
  CLOSE (9)
  OPEN(7,file='numb-2.sig',status='old')
  OPEN(8,file='numb10.sig',status='old')
  OPEN(9,file='numball.sig',status='unknown')
3 READ(8,200,end=4) name,a
4 READ(7,200,end=98) name1,b
  IF(name.eq.name1)THEN
     int=10
     WRITE(9,201)name,a,int
  ELSE
     int=-2
     WRITE(9,201)name1,b,int
     GOTO 4
  ENDIF
  GOTO 3
200 FORMAT(a9,f12.8,f11.7,f10.7,f10.6,f11.6,f10.6,f13.4)
201 FORMAT(a9,f12.8,f11.7,f10.7,f10.6,f11.6,f10.6,f13.4,i4)
98 CLOSE (7)
  CLOSE (8)
  CLOSE (9)
  OPEN(7,file='numb-2.del',status='old')
  OPEN(8,file='numb10.del',status='old')
  OPEN(9,file='numball.del',status='unknown')
5 READ(8,300,end=6) name,a
6 READ(7,300,end=97) name1,b
  IF(name.eq.name1)THEN
     int=10
     WRITE(9,301)name,a,int
  ELSE
     int=-2
     WRITE(9,301)name1,b,int
     GOTO 6
  ENDIF
  GOTO 5
300 FORMAT(a9,f12.8,f12.7,f10.7,f10.6,f11.6,f10.6,f13.4)
301 FORMAT(a9,f12.8,f12.7,f10.7,f10.6,f11.6,f10.6,f13.4,i4)
  !500     format(a1,a9,f17.6,f16.12,f13.10,4f13.7,f6.2)
500 FORMAT(a1,a10,f16.6,6e25.16e3,f6.2)
97 CONTINUE
END PROGRAM merg

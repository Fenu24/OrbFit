! program finds Trojan asteroids for which extended 
! integration if necessary
PROGRAM selextt
  IMPLICIT NONE
  INTEGER len, lench, b1(2)
  CHARACTER*9 name,name1,initia
  CHARACTER*20 filnam
  DOUBLE PRECISION a(8),b(7)
  !   filtered data file name
  WRITE(*,*)' initial part of file names?'
  READ(*,*)initia
  ! ==================================================
  !   main output files
  ! ==================================================
  CALL cmpstr(initia,len)
  !      filnam=initia(1:len)//'prt.pro'
  !      OPEN(7,file=filnam,status='old')
  filnam=initia(1:len)//'prt.del'
  OPEN(8,file=filnam,status='old')
  filnam=initia(1:len)//'prt.list'
  OPEN(9,file=filnam,status='unknown')
  !1     read(7,100,end=99) name,b,b1
  !1     read(7,'(a9,f7.4,f7.2,f7.3)',end=99) name,b(1),b(2),b(3)
  !      write(*,*) name, b(1)
  !      goto 1
1 READ(8,101,end=99) name1,a
  WRITE(*,*)name1
  !     if(name.eq.name1)then
  IF(a(1).gt.1.d-3.or.a(4).gt.2.5d-3.or.    &
       &                a(6).gt.2.5d-3.or.a(8).gt.6.d-1)THEN
     WRITE(9,102)name1
  ENDIF
  !     else
  !     stop 'Alignment problem'
  !     endif
  GOTO 1
  !100    format(a10,f6.4,f7.2,f7.3,f6.4,f8.2,f6.4,f8.2,2i3)
101 FORMAT(a10,f7.5,f8.3,f8.4,f7.5,f9.3,f7.5,f9.3,f6.2)
  !102    format(' ',a9)
102 FORMAT(a10)
99 CLOSE (7)
  CLOSE (8)
  CLOSE (9)
END PROGRAM selextt

! program finds asteroids for which extended integration is necessary   
PROGRAM selext
  IMPLICIT NONE
  INTEGER len, lench
  CHARACTER*9 name,name1,initia
  CHARACTER*20 filnam
  DOUBLE PRECISION a(7),b(7)
  !   filtered data file name
  WRITE(*,*)' initial part of file names?'
  READ(*,*)initia
  ! ==================================================
  !   main output files
  ! ==================================================
  CALL cmpstr(initia,len)
  filnam=initia(1:len)//'prs.pro'
  OPEN(7,file=filnam,status='old')
  filnam=initia(1:len)//'prs.sig'
  OPEN(8,file=filnam,status='old')
  filnam=initia(1:len)//'prs.list'
  OPEN(9,file=filnam,status='unknown')
1 READ(7,100,end=99) name,b
  READ(8,101,end=99) name1,a
  WRITE(*,*)name
  IF(name.eq.name1)THEN
     IF(a(1).gt.3.d-4.or.a(2).gt.3.d-3.or.   &
          &       a(3).gt.1.d-3.or.-b(7).gt.5.d1)THEN
        WRITE(9,102)name
     ENDIF
  ELSE
     STOP 'Alignment problem'
  ENDIF
  GOTO 1
100 FORMAT(a9,f12.7,f11.7,f10.7,f12.6,f13.6,f13.6,f8.2)
101 FORMAT(a9,f12.8,f11.7,f10.7,f10.6,f11.6,f10.6,f13.4)
102 FORMAT(' ',a9)
99 CLOSE (7)
  CLOSE (8)
  CLOSE (9)
END PROGRAM selext

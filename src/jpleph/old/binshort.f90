!                                                                       
!   binshort : a program to shorten a binary JPL Ephemeris file         
!                                                                       
!   NOTE that the user must set a few parameters                        
!                                                                       
!    (between the lines of asterisks)                                   
!                                                                       
                                                                        
                                                                        
      IMPLICIT DOUBLE PRECISION (A-H,O-Z) 
                                                                        
      CHARACTER*6 TTL(14,3),CNAM(400) 
                                                                        
      DIMENSION SS(3),CVAL(400),data(3000) 
                                                                        
      INTEGER IPT(3,12),idt(6),LPT(3) 
                                                                        
      character*3 mn(12) 
      data mn/'JAN','FEB','MAR','APR','MAY','JUN',                      &
     &        'JUL','AUG','SEP','OCT','NOV','DEC'/                      
                                                                        
! *************************************************************         
! *************************************************************         
!                                                                       
!  *****  Set the value of NRECL  *****                                 
!                                                                       
!  NRECL=1 if 'RECL' in the OPEN statement is to be in s.p. words       
!                                                                       
!  NRECL=4 if 'RECL' in the OPEN statement is to be in bytes            
!                                                                       
      data nrecl/4/ 
!                                                                       
! *************************************************************         
!                                                                       
!  *****  Set the value of KSIZE  *****                                 
!                                                                       
!   for DE200, use 1652; for DE405, use 2036; for DE406, use 1456       
!                                                                       
      data ksize/2036/ 
                                                                        
! *************************************************************         
                                                                        
!  *****  SET the start and stop times (JED) of new file  *****         
!                                                                       
! Beginning of 1970 to beginning of 2030.                               
      t1=2440587.5 
      t2=2462502.5 
!                                                                       
! *************************************************************         
! *************************************************************         
                                                                        
      irecl=ksize*nrecl 
                                                                        
        OPEN(11,                                                        &
     &       FILE='OLDEPH',                                             &
     &       ACCESS='DIRECT',                                           &
     &       FORM='UNFORMATTED',                                        &
     &       RECL=irecl,                                                &
     &       STATUS='OLD')                                              
                                                                        
                                                                        
                                                                        
      READ(11,REC=1)TTL,CNAM,SS,NCON,AU,EMRAT,IPT,NUMDE,LPT 
      READ(11,REC=2)CVAL 
                                                                        
      if(t1 .lt. ss(1)) t1=ss(1) 
      if(t2 .gt. ss(2)) t2=ss(2) 
                                                                        
        OPEN(10,                                                        &
     &       FILE='NEWEPH',                                             &
     &       ACCESS='DIRECT',                                           &
     &       FORM='UNFORMATTED',                                        &
     &       RECL=irecl,                                                &
     &       STATUS='NEW')                                              
                                                                        
                                                                        
      nrp=2 
                                                                        
      ndata=ksize/2 
                                                                        
      nrec=2 
                                                                        
    5 nrec=nrec+1 
      READ(11,REC=NREC,ERR=98)(DATA(K),K=1,NDATA) 
                                                                        
      if(data(2) .le. t1) go to 5 
                                                                        
      if(nrp .eq. 2) ss(1)=data(1) 
                                                                        
      nrp=nrp+1 
                                                                        
      WRITE(10,REC=NRP,ERR=99)(DATA(K),K=1,NDATA) 
                                                                        
      if(data(2) .lt. t2) go to 5 
                                                                        
      ss(2)=data(2) 
                                                                        
      do j=1,2 
      call joulex(ss(j),idt) 
      k=idt(2) 
      rewind 9 
      write(9,'(a20,f9.1,i6,1x,a3,i3)')                                 &
     & 'Start Epoch : JED = ',ss(j),idt(1),mn(k),idt(3)                 
                                                                        
      rewind 9 
      read(9,'(14a6)')(ttl(i,j+1),i=1,14) 
      enddo 
                                                                        
      ttl(1,3)='Stop  ' 
                                                                        
                                                                        
      WRITE(10,REC=1)TTL,CNAM,SS,NCON,AU,EMRAT,IPT,NUMDE,LPT 
      WRITE(10,REC=2)CVAL 
                                                                        
      close(10) 
      close(11) 
                                                                        
      stop 
                                                                        
   98 write(*,'(/a18,i3)')'read error on input file' 
      stop 
                                                                        
   99 write(*,                                                          &
     &  '(/''write error on output file, nrec='',i6)')nrp               
      stop 
                                                                        
      END                                           
                                                                        
                                                                        
                                                                        
                                                                        
      SUBROUTINE JOULEX(DJED,IDATE) 
!                                                                       
      IMPLICIT double precision (A-H,O-Z) 
!                                                                       
!       CONVERTS JULIAN DAY NUMBER TO CALENDAR DATE                     
!     IDATE(6)=MONTH,DAY,YEAR,HOUR,MINUTE,SECOND (INTEGER) FOR CALENDAR 
!     DJED=DOUBLE PRECISION JULIAN DAY NUMBER TO BE INPUT               
!     FLAG=(-1,0,+1(=)JULIAN,AUTOMATIC,GREGORIAN) CALENDAR DATE OUTPUT  
!       AUTOMATIC MODE SELECTS JULIAN CALENDAR FOR DAYS ON AND BEFORE   
!       DJED=2299160=1582OCT4, GREGORIAN AFTER THIS DATE = 1582OCT15.   
!       OTHER MODES EMPLOY CALENDAR INDICATED.                          
!     **RESTRICTIONS - INPUT JULIAN CALENDAR MUST BE FOR A POSITIVE JULI
!       DAY NUMBER.  FOR GREGORIAN CALENDAR, DJED.GT.1721119 (DATES A.D.
!     **SEE SUBROUTINE CTOJ FOR INVERSE FUNCTIONS AND REFERENCES TO     
!       CALENDAR DEFINITIONS, TO BE FOUND IN EXP.SUPP. PP.412, 436.     
!                                                                       
      INTEGER IDATE(6),DY,Y,  D,H,M,J,FLAG,MD 
      double precision DJED,JA 
                                                                        
      DATA FLAG/0/ 
!                                                                       
!                                                                       
    1 CONTINUE 
                                                                        
      DT=DJED+.5D0+0.5D0/86400.D0 
      J=DT 
      JA=J 
      DT=(DT-JA)*24.D0 
      H=DT 
      JA=H 
      DT=(DT-JA)*60.D0 
      MD=DT 
      JA=MD 
      DT=(DT-JA)*60.D0 
      SECX = DT 
      SECS = SECX 
                                                                        
      IF(FLAG)20,10,30 
                                                                        
   10 IF (J.GT.2299160) GO TO 30 
                                                                        
   20 Y=4*(J+105498)/1461-5001 
      DY=J-(36525*(Y+5000))/100 + 105133 
      GO TO 40 
                                                                        
   30 Y=4*(J-1720754+3*(4*(J-1721119)/146097+1)/4)/1461-1 
      DY=J-36525*Y/100+3*(Y/100+1)/4-1721119 
                                                                        
   40 M=(5*DY-3)/153+3 
      D=DY-(153*(M-3)+2)/5 
      IF (M.GE.13) M=M-12 
      IF(M.LE.2) Y=Y+1 
                                                                        
      IDATE(1) = Y 
      IDATE(2) = M 
      IDATE(3) = D 
      IDATE(4) = H 
      IDATE(5)=MD 
      IDATE(6)=SECS 
                                                                        
      RETURN 
      END                                           

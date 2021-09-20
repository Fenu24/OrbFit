!  Program "binmerge" : Merge two binary ephemeris files                
!                                                                       
!                                                                       
!   NOTE that the user must set the values of NRECL and KSIZE           
!                                                                       
!    (between the lines of asterisks)                                   
!                                                                       
!     Input files : INEPH1, INEPH2                                      
!                                                                       
!       (the stop-date of INEPH1 must be .ge. the start-date of INEPH2) 
!                                                                       
!     Output file : OUTEPH                                              
!                                                                       
!                                                                       
                                                                        
      IMPLICIT DOUBLE PRECISION (A-H,O-Z) 
                                                                        
      CHARACTER*6 TTL1(14,3),TTL2(14,3),CNAM1(400),CNAM2(400) 
                                                                        
      DOUBLE PRECISION SS1(3),SS2(3),CVAL1(400),CVAL2(400),data(3000) 
                                                                        
      INTEGER IPT1(3,12),IPT2(3,12),LPT1(3),LPT2(3) 
                                                                        
! *************************************************************         
! *************************************************************         
!                                                                       
!  *****  Set the value of NRECL  *****                                 
!                                                                       
!  NRECL=1 if 'RECL' in the OPEN statement is to be in s.p. words       
!                                                                       
!  NRECL=4 if 'RECL' in the OPEN statement is to be in bytes            
!                                                                       
      data nrecl/ / 
!                                                                       
! *************************************************************         
!                                                                       
!  *****  Set the value of KSIZE  *****                                 
!                                                                       
!   for DE200, use 1652; for DE405, use 2036; for DE406, use 1456       
!                                                                       
      data ksize/    / 
                                                                        
! *************************************************************         
! *************************************************************         
                                                                        
      irecl=ksize*nrecl 
                                                                        
        OPEN(11,                                                        &
     &       FILE='INEPH1',                                             &
     &       ACCESS='DIRECT',                                           &
     &       FORM='UNFORMATTED',                                        &
     &       RECL=irecl,                                                &
     &       STATUS='OLD')                                              
                                                                        
                                                                        
                                                                        
      READ(11,REC=1)TTL1,CNAM1,SS1,NCON1,AU1,EMRAT1,IPT1,NUMDE1,LPT1 
      READ(11,REC=2)CVAL1 
                                                                        
                                                                        
        OPEN(12,                                                        &
     &       FILE='INEPH2',                                             &
     &       ACCESS='DIRECT',                                           &
     &       FORM='UNFORMATTED',                                        &
     &       RECL=irecl,                                                &
     &       STATUS='OLD')                                              
                                                                        
                                                                        
      READ(12,REC=1)TTL2,CNAM2,SS2,NCON2,AU2,EMRAT2,IPT2,NUMDE2,LPT2 
      READ(12,REC=2)CVAL2 
                                                                        
                                                                        
!  *****  Check for allowable time-spans  *****                         
                                                                        
      if(ss1(1) .ge. ss2(2))write(*,'(a50)')                            &
     & 'time-spans reversed: switch assignments?'                       
                                                                        
      if(ss1(2) .lt. ss2(1)) write(*,'(a50)')                           &
     & 'time-spans neither abut nor overlap'                            
                                                                        
                                                                        
!  *****  Check for matching ephemerides  *****                         
                                                                        
      if (numde1 .ne. numde2) write(*,'(a50)')'NUMDE''s differ' 
                                                                        
      if (ncon1 .ne. ncon2) write(*,'(a50)')'NCON''s differ' 
                                                                        
      if (au1 .ne. au2) write(*,'(a50)')'AU''s differ' 
                                                                        
      if (emrat1 .ne. emrat2) write(*,'(a50)')'EMRAT''s differ' 
                                                                        
      do i=1,ncon1 
      if (cval1(i) .ne. cval2(i)) write(*,'(a50)')'CVAL''s differ' 
      if (cnam1(i) .ne. cnam2(i)) write(*,'(a50)')'CNAM''s differ' 
      enddo 
                                                                        
      do i=1,12 
      do j=1,3 
      if (ipt1(i,j) .ne. ipt2(i,j)) write(*,'(a50)')'IPT''s differ' 
      if (lpt1(j) .ne. lpt2(j)) write(*,'(a50)')'LPT''s differ' 
      enddo 
      enddo 
                                                                        
                                                                        
!  *****  Combine and Write initial header info  *****                  
                                                                        
      ss1(2)=ss2(2) 
      do i=1,14 
      ttl1(i,3)=ttl2(i,3) 
      enddo 
                                                                        
        OPEN(10,                                                        &
     &       FILE='OUTEPH',                                             &
     &       ACCESS='DIRECT',                                           &
     &       FORM='UNFORMATTED',                                        &
     &       RECL=irecl,                                                &
     &       STATUS='NEW')                                              
                                                                        
                                                                        
                                                                        
      WRITE(10,REC=1)TTL1,CNAM1,SS1,NCON1,AU1,EMRAT1,IPT1,NUMDE1,LPT1 
      WRITE(10,REC=2)CVAL1 
      nrp=2 
                                                                        
      ndata=ksize/2 
                                                                        
      dz=ss1(1) 
                                                                        
      do 6 nfil=11,12 
                                                                        
      nrec=2 
                                                                        
    5 nrec=nrec+1 
      READ(nfil,REC=NREC,END=6,ERR=98)(DATA(K),K=1,NDATA) 
      if(data(1) .lt. dz) go to 5 
      if(data(1) .gt. dz) write(*,'(a50)')'non-matching dates' 
      dz=data(2) 
                                                                        
      nrp=nrp+1 
      WRITE(10,REC=NRP,ERR=99)(DATA(K),K=1,NDATA) 
                                                                        
      go to 5 
                                                                        
    6 continue 
                                                                        
      close(10) 
      close(11) 
      close(12) 
                                                                        
      stop 
                                                                        
   98 write(*,'(/a18,i3)')'read error on unit',nfil 
      stop 
                                                                        
   99 write(*,                                                          &
     &  '(/''write error on output tape, nrec='',i6)')nrp               
      stop 
                                                                        
      END                                           

! ========================================                              
!  REACAT                                                               
!  reads from catalogues                                                
! form unit 3                                                           
SUBROUTINE reacat(iprqua,sysz,refz,coox,sysx,refx,t0a,t1,xa,nfound)
  USE fund_const
  USE massmod
  USE yorp_module, ONLY: sma_ast
  IMPLICIT NONE 
! ============INPUT=================                                    
! output control                                                        
  INTEGER iprqua 
  CHARACTER*3 sysz 
  CHARACTER*6 refz 
! ===========OUTPUT================                                     
  CHARACTER*3 coox,sysx 
  CHARACTER*6 refx 
! epoch times                                                           
  DOUBLE PRECISION t0a,t1,tp 
  DOUBLE PRECISION xa(6,nastx) 
  INTEGER nfound 
! ========END INTERFACE================                                 
! loop indexes                                                          
  INTEGER j,i,jj 
! input from reaele                                                     
  LOGICAL eof,fou(nastx) 
  INTEGER nrec,nrecx 
  PARAMETER (nrecx=1000000) 
  DOUBLE PRECISION hm,xx(6) 
  INTEGER ndarc,nobs 
! headers                                                               
  INCLUDE 'comnbo.h90' 
! header and comment strings                                            
  character*60 comast 
  character*100 colhep,colhea 
! units                                                                 
  CHARACTER*3 unitx 
! asteroid names                                                        
  CHARACTER*19 namid 
  CHARACTER*9 nam0 
! output of rdorb                                                       
  character*5 epoch 
  character*4 rsys 
  INTEGER nrecord,ll 
! definition flags                                                      
  LOGICAL cov0,defnor 
! opposition effect, mass, matrices                                     
  DOUBLE PRECISION gma0,mass,g0(6,6),c0(6,6) 
! ========================================                              
!  (1a) input header from catalogue in old format:                      
!     read(3,103)t0a,coox,unitx,sysx,refx                               
!103  format(4x,f9.1,1x,a3,1x,a3,1x,a3,1x,a6/)                          
!                                                                       
!     write(9,420)coox,sysx,refx,unitx                                  
!420  format(' asteroids:'/' input from catalogue, coox=',a3,           
!    +       ' sysx=',a3,' refx=',a6,' unitx=',a3)                      
! comment                                                               
  comast=' from catalogue ' 
! main loop on requested asteroids                                      
  nfound=0 
  DO j=1,na 
     fou(j)=.false. 
  ENDDO
  DO 7 nrec=1,nrecx 
!  (4a) input elements:                                                 
     CALL rdorb(namid,xx,coox,t0a,g0,cov0,c0,defnor,                &
     &           hm,gma0,mass,rsys,epoch,nrecord,eof)                   
     if(eof)then 
        write(9,*)' file end looking for ', namid, ' record ',nrec 
        goto 8 
     endif
! check reference system                                                
     IF(rsys.ne.'ECLM'.or.epoch.ne.'J2000')THEN 
        WRITE(*,*)' incompatible reference systems ' 
        WRITE(*,*)' ref=',rsys,' epoch=',epoch 
        WRITE(*,*)' should be ECLM, J2000' 
        STOP 
     ELSE 
        refx='ECLM00' 
        sysx='HEL' 
        unitx='RAD' 
        t1=0.d0 
     ENDIF
! handle funny identification names                                     
! convert name in case it is of the form nam0=namp                      
     ll=index(namid,'=')-1 
     IF(ll.lt.0)THEN 
        nam0=namid(1:9) 
     ELSE 
        nam0=namid(1:ll) 
     ENDIF
!        call reaele(3,eof,id,xx,ndarc,nobs,hm)                         
!        if(eof) GOTO 8                                                 
! check if this is required                                             
     DO 5 j=1,na 
        IF(nam0.eq.ida(j))THEN 
!  (5a) identification of a required asteroid                           
           nfound=nfound+1 
           DO i=1,6 
              xa(i,j)=xx(i) 
           ENDDO
           write(9,*)' asteroid ',nam0,'  OK' 
           ! Save the semimajor axis
           sma_ast(j) = xx(1)
           fou(j)=.true. 
        ENDIF
! no number of revolutions in catalogues                                
!                                                                       
5    ENDDO
7 ENDDO
! check if list has been exhausted                                      
8 IF(nfound.ne.na)THEN 
     WRITE(9,*)' not all required asteroids found' 
     DO j=1,na 
        IF(.not.fou(j))WRITE(9,*)' missing ',ida(j) 
     ENDDO
     STOP 
  ENDIF
! (7) conversion of angles to radiants                                  
  DO j=1,na 
     call radiant(xa(1,j),coox,unitx,dpig) 
  ENDDO
! time conversion MJD -> JD                                             
  t0a=t0a+2400000.5d0 
! ********************************************************************  
!  write headers in output files:                                       
!  one file for all the planets and one file for all the asteroids      
! ********************************************************************  
!  column headers for equinoctal elements and numbers                   
  colhea='  a(au)      h      k      p      q     lambda  nrev  gam' 
!  sampled elements, only for iprq even, unit 11 and 21                 
  if(mod(iprqua,2).eq.0)then 
     open(21,file='vast.dat',status='unknown') 
     call wrihea(21,ida,comast,colhea, 'EQU',sysz,refz,'RAD',na,nvz,t0a) 
  endif
!  filtered elements, always, unit 12 and 22                            
  open(22,file='vast.fil',status='unknown') 
  call wrihea(22,ida,comast,colhea, 'EQU',sysz,refz,'RAD',na,nvz,t0a) 
!  dump files                                                           
  open(14,file='orb9.dma',status='unknown') 
  call wrihea(14,ida,comast,colhea, 'EQU',sysz,refz,'RAD',na,nvz,t0a)
END SUBROUTINE reacat

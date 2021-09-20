! ========================================                              
! REAOUT                                                                
!  read input from orbit9 output file                                   
SUBROUTINE reaout(iprqua,sysz,refz,coox,sysx,refx,t0a,ta          &
     &      ,xa,nga,vslog,jjin,xv,nfound)
  USE massmod
  USE fund_const                             
  IMPLICIT NONE 
! =======INPUT==========================                                
! output control                                                        
  INTEGER iprqua 
  CHARACTER*3 sysz 
  CHARACTER*6 refz 
! =======OUTPUT=========================                                
! input coordinate system                                               
  CHARACTER*3 coox,sysx 
  CHARACTER*6 refx 
! epoch times                                                           
  DOUBLE PRECISION t0a,t1 
! asteroid elements                                                     
  DOUBLE PRECISION xa(6,nastx) 
  INTEGER nga(nangx,nastx) 
! normalisation controls                                                
  DOUBLE PRECISION vslog(nvzx) 
  INTEGER nfound 
! ========END INTERFACE================                                 
! arrays for initial conditions: epoch, asteroids                       
  DOUBLE PRECISION ta, ya(6,nastx) 
!  names and numbers of input planets and asteroids, comments           
  character*9 numbin(nastx) 
! headers                                                               
  INCLUDE 'comnbo.h90' 
! dimensions, controls                                                  
  INTEGER nang,nact,ilca,nast 
! functions                                                             
  INTEGER numang,numact 
! units                                                                 
  CHARACTER*3 unitx 
! arrays for variational equations                                      
  DOUBLE PRECISION xv(6,nvzx),jjin(nastx),ngz(nangx,nastx) 
  DOUBLE PRECISION gamma(nastx) 
! header and comment strings                                            
  character*60 comast 
  character*100 colhep,colhea 
! loop indexes                                                          
  INTEGER j,jj,i,jv 
! =============================================                         
!  (1b) input header from orbit8v format:                               
  call reahea(3,numbin,comast,colhea,coox,sysx,refx,unitx,nast,ilca,t0a)
  write(9,421)(numbin(j),j=1,nast),coox,sysx,refx,unitx 
421 format(' asteroids:',                                             &
     &   ' input from previous output, available:',                     &
     &   (/8(2x,a9))/' coox=',a3,                                       &
     &       ' sysx=',a3,' refx=',a6,' unitx=',a3)                      
! read asteroid initial conditions file:                                
  nang=numang(coox) 
  nact=numact(coox) 
!  input from output file                                               
  if(nast.gt.1)then 
     read(3,*,err=78,end=78)ta 
  endif
!  (4b) input all the elements, see later the good ones                 
  if(nast.gt.1)then 
     DO 14 jj=1,nast 
        if(jj.le.ilca)then 
           read(3,*,end=79,err=78)(ya(i,jj),i=1,nact)                &
     &        ,(ya(i+nact,jj),ngz(i,jj),i=1,nang),gamma(jj)             
        else 
           read(3,*,end=79,err=78)(ya(i,jj),i=1,nact)                &
     &        ,(ya(i+nact,jj),ngz(i,jj),i=1,nang)                       
        endif
14   ENDDO
  else 
     if(ilca.eq.0)then 
        read(3,*,end=79,err=78)ta,(ya(i,1),i=1,nact)                &
     &      ,(ya(i+nact,1),ngz(i,1),i=1,nang)                           
     else 
        read(3,*,end=79,err=78)ta,(ya(i,1),i=1,nact)                &
     &      ,(ya(i+nact,1),ngz(i,1),i=1,nang),gamma(1)                  
     endif
  endif
!  input all the available variation vectors                            
  DO  jv=1,ilca 
     read(3,*,end=79,err=78)(xv(i,jv),i=1,6) 
  ENDDO
!  identify the required asteroids in the list numbin                   
  DO 15 j=1,na 
     DO 16 jj=1,nast 
!  the list of asteroids in numbin can be out of order, scan all        
        if(ida(j).eq.numbin(jj))then 
!  (5b) identification of a required asteroid                           
           write(9,*)' asteroid ',ida(j),'  OK' 
           jjin(j)=jj 
           DO i=1,6 
              xa(i,j)=ya(i,jj) 
           ENDDO
           DO i=1,nang 
              nga(i,j)=ngz(i,jj) 
           ENDDO
!  (6b) initial value of divergence ratio                               
           if(jj.le.ilca)then 
              vslog(j)=gamma(jj) 
           else 
              vslog(j)=0 
           endif
!c           write(9,*)jj,(xa(i,jj),i=1,6)                              
           GOTO 115 
        endif
16   ENDDO
     goto 7 
115  CONTINUE
15 ENDDO
! (7) conversion of angles to radiants                                  
  DO j=1,na 
     call radiant(xa(1,j),coox,unitx,dpig) 
  ENDDO
! ********************************************************************  
!  write headers in output files:                                       
!  one file for all the planets and one file for all the asteroids      
! ********************************************************************  
!  column headers for equinoctal elements and numbers                   
  colhea='  a(au)      h      k      p      q     lambda  nrev  gam' 
!  sampled elements, only for iprq even, unit 11 and 21                 
  if(mod(iprqua,2).eq.0)then 
     open(21,file='vast.dat',status='unknown') 
     call wrihea(21,ida,comast,colhea,'EQU',sysz,refz,'RAD',na,nvz,t0a)
  endif
!  filtered elements, always, unit 12 and 22                            
  open(22,file='vast.fil',status='unknown') 
  call wrihea(22,ida,comast,colhea,'EQU',sysz,refz,'RAD',na,nvz,t0a)
!  dump files                                                           
  open(14,file='orb9.dma',status='unknown') 
  call wrihea(14,ida,comast,colhea,'EQU',sysz,refz,'RAD',na,nvz,t0a)
  RETURN 
!  asteroid not found                                                   
78 continue 
79 continue 
7 write(9,*)' asteroid not found ',ida(j),' last was ',numbin(jj) 
  stop 
END SUBROUTINE reaout

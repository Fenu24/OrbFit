! **********************************************************            
!                                                                       
!     ORBIT9C  Pisa, October 1998                                    
!                                                                       
! ***********************************************************           
PROGRAM orbit9 
  USE massmod
  USE fir_filters
  USE force9d
  USE propag_state, ONLY: propin, hms
  IMPLICIT NONE 
  INCLUDE 'comnbo.h90' 
! arrays for initial conditions and elements, output                    
  double precision yp(6,nbox),ennep(nbox),elp(6,nbox),plv(nbox) 
! center of mass of inner solar system                                  
  double precision bar(6) 
  double precision ya(6,nastx),ennea(nastx),ela(6,nastx),alv(nastx) 
  integer ngip(nbox),ngia(nastx),ngpp(nbox),ngpa(nastx) 
!  cartesian position and velocity vectors                              
  double precision y1(nvarx),y2(nvarx) 
  integer nskk,nskip,nskipa  
!  array for filter input/output                                        
  double precision eq(nvarx),eqf(nvarx),angl(norbx),angf(norbx) 
  integer ng(norbx),ngf(norbx) 
!  for variational equations                                            
  double precision gamma(nvzx),vslog(nvzx),vnor(nvzx),yv(6,nastx) 
!  identifiers of coordinate systems for output                         
  character*3 sysz 
  character*6 refy,refz 
  logical outf,outfa,reno 
! dimensions                                                            
  integer ndim,ndim2,nvar2,nvar,npla,norb,ncof 
! controls                                                              
  integer nout,idump,iprqua,nin,nfl,iout0,ioucou 
  double precision v0,v1,vpos,vvel 
! time                                                                  
  double precision dt,t,t1,t2,tf,tfa,dinf, ddd 
! loop indexes                                                          
  integer j,ja,jq,jj,iout,i 
! scalar temporaries                                                    
  double precision sv,vv,renorm 
! *********************************                                     
! only to test derivatives                                              
  double precision rinc,step,f0(nvar2x),f1(nvar2x) 
  integer ii 
! close approach data                                                   
  integer idc 
  double precision xxpla(6) 
! *********************************************************             
  call libini 
!   input                                                               
  call inpo9(dt,nout,idump,iprqua,                            &
     &             sysz,refz,refy,                                      &
     &             v0,v1,vpos,vvel,vslog,                               &
     &             t1,yp,ya,ngip,ngia,yv)                               
!  number of components                                                 
  ndim=3 
  ndim2=ndim*2 
  nvar2=ndim*(nbod-1+na+nvz) 
  nvar=2*nvar2 
  npla=nbod-1 
  norb=nbod-1+na 
!  number of real components to be filtered                             
  ncof=norb*(2*ndim-1)+nvz 
  GOTO 99 
! ==================================================                    
! test derivatives; d(i)/d(ii)                                          
  step=1.d-8 
  idc=0 
  DO 90 j=1,na 
     write(*,*)(ya(i,j),i=1,3) 
     DO 91 ii=1,3 
! -----------------------------------                                   
! set variation vector to  E_ii                                         
        DO i=1,3 
           yv(i,j)=0.d0 
           yv(i+3,j)=0.d0 
        ENDDO
        yv(ii,j)=1.d0 
! reorder                                                               
        call stack(y1,nvar,yp,npla,ya,na,yv,nvz,ndim,ndim2) 
! compute derivatives                                                   
        CALL force9(y1,y1(nvar2+1),t1,f0,nvar2,idc,xxpla,0,1) 
! increment variable ii by step                                         
        ya(ii,j)=ya(ii,j)+step 
! reorder                                                               
        call stack(y1,nvar,yp,npla,ya,na,yv,nvz,ndim,ndim2) 
! compute increment                                                     
        CALL force9(y1,y1(nvar2+1),t1,f1,nvar2,idc,xxpla,0,1) 
        DO 94 i=1,3 
! compute incremental ratio                                             
           rinc=(f1(3*npla+3*j-3+i)-f0(3*npla+3*j-3+i))/step 
           write(*,*)i,ii,f0(3*norb+3*j-3+i),rinc,                     &
     &            (f0(3*norb+3*j-3+i)-rinc)/rinc                        
94      ENDDO
91   ENDDO
! ----------------------------------                                    
90 ENDDO
  STOP 
99 CONTINUE 
! *********************************************************             
!  initial conditions: all arrays copied in the state                   
!  vector                                                               
  CALL stack(y1,nvar,yp,npla,ya,na,yv,nvz,ndim,ndim2) 
! ****************************************************                  
!  output at initial time:                                              
!                                                                       
!  planets                                                              
  CALL coord(yp,'BAR',nbod,'CAR',refy,elp,sysz,'EQU',refz,ennep,bar) 
!  revolution counters:                                                 
  DO j=1,nbod-1 
     plv(j)=elp(6,j) 
  ENDDO
! arrays for filtering and output                                       
  jq=0 
  DO 58 j=1,nbod-1 
     angl(j)=plv(j) 
     ng(j)=ngip(j) 
     DO i=1,5 
        eq(i+jq)=elp(i,j) 
     ENDDO
     jq=jq+5
58 ENDDO
!  asteroids                                                            
  CALL cooast(ya,'BAR',na,'CAR',refy,bar,nbod,ela,sysz,'EQU',refz,ennea)
  DO 6 ja=1,na 
!  revolution counters:                                                 
     alv(ja)=ela(6,ja) 
! arrays for filtering and output                                       
     angl(ja+nbod-1)=alv(ja) 
     ng(ja+nbod-1)=ngia(ja) 
     DO i=1,5 
        eq(i+jq)=ela(i,ja)
     ENDDO
     jq=jq+5 
6 ENDDO
  DO  j=1,nvz 
     eq(jq+j)=vslog(j)
  ENDDO
  IF(iprqua.gt.1)THEN 
     CALL outsub(11,21,t1,eq,angl,ng,ncof,norb,na,nvz,.true.) 
  ENDIF
  CALL outsub(13,14,t1,eq,angl,ng,ncof,norb,na,nvz,.false.) 
  if(nvz.ne.0)then 
     DO j=1,nvz 
        write(14,*)(yv(i,j),i=1,ndim) 
        write(14,*)(yv(i+ndim,j),i=1,ndim) 
     ENDDO
  endif                                                  
! ****************************************************                  
!  initialisations                                                      
!                                                                       
!  interval to filter input                                             
  dinf=dt/nsamp 
  nin=nout*nsamp 
!  initial time is t1; how many filter inputs are behind?               
  iout0=t1/dinf 
!  starter is needed to start                                           
  nfl=0 
!  output counter                                                       
  ioucou=0 
!  filter skip (.fil remains in phase)                                  
  nskk=mod((nfil+1)/2,nsamp) 
  if(nskk.eq.0)then 
     nskip=0 
  else 
     nskip=nsamp-nskk 
  endif
  nskipa=nskip 
!  hlen=dinf*(nfil-1)/2 
!                                                                       
! *********************************************************             
!   main loop begins here                                               
!                                                                       
!  propagation to next filter input time                                
  DO 10 iout=iout0+1,nin 
     t2=dinf*iout 
     CALL propin(nfl,y1,t1,t2,y2,hms,nvar,6)
     nfl=1 
     t=t2 
! *********************************************************             
!  output: reordering of the arrays                                     
     CALL reord(y2,nvar,yp,npla,ya,na,yv,nvz,ndim,ndim2) 
!  planets: coordinate change                                           
     CALL coord(yp,'BAR',nbod,'CAR',refy,elp,sysz,'EQU',refz,ennep,bar) 
!  asteroids: coordinate change                                         
     CALL cooast(ya,'BAR',na,'CAR',refy,bar,nbod,ela,sysz,'EQU',refz,ennea)
!      ddd=sqrt((ya(1,1)-ya(1,2))**2+(ya(2,1)-ya(2,2))**2+              
!     +             (ya(3,1)-ya(3,2))**2)                               
!      WRITE(23,123)t2,ddd                                              
! 123  FORMAT(F15.4,1x,1P,d12.5)                                        
!  planets                                                              
     jq=0 
     DO 17 j=1,nbod-1 
!  planets: revolution counters update                                  
        call predicl(plv(j),ngip(j),dinf,ennep(j),ngpp(j)) 
        ngip(j)=ngpp(j) 
        call prngup(elp(6,j),plv(j),ngip(j)) 
!  planets: arrays for filtering and output                             
        angl(j)=plv(j) 
        ng(j)=ngip(j) 
        DO i=1,5 
           eq(i+jq)=elp(i,j)
        ENDDO
        jq=jq+5 
17   ENDDO
!  asteroids                                                            
     DO 16 ja=1,na 
!  asteroids: revolution counters update                                
        call predicl(alv(ja),ngia(ja),dinf,ennea(ja),ngpa(ja)) 
        ngia(ja)=ngpa(ja) 
        call prngup(ela(6,ja),alv(ja),ngia(ja)) 
!  asteroids: arrays for filtering and output                           
        angl(ja+nbod-1)=alv(ja) 
        ng(ja+nbod-1)=ngia(ja) 
        DO i=1,5 
           eq(i+jq)=ela(i,ja) 
        ENDDO
        jq=jq+5 
16   ENDDO
!  variational equations                                                
     DO 20 j=1,nvz 
!  norm of the variation vector; here the norm**2 is defined by         
!  weight vpos for positions and vvel for velocities                    
        sv=0.d0 
        vv=0.d0 
        DO i=1,ndim 
           sv=sv+yv(i,j)**2 
           vv=vv+yv(i+ndim,j)**2 
        ENDDO
        vnor(j)=dsqrt(sv*vpos+vv*vvel) 
!  log of divergence ratio; no /(t-t0) used here                        
        gamma(j)=(dlog(vnor(j))-dlog(v0)+vslog(j)) 
        jq=jq+1 
        eq(jq)=gamma(j) 
20   ENDDO
! ***************************************************************       
!  input filter                                                         
     CALL filter(t,eq,ncof,nskip,dinf,tf,eqf,outf) 
     CALL filtan(t,angl,ng,norb,nskipa,dinf,tfa,angf,ngf,outfa)
! ****************************************************************      
! output control                                                        
!                                                                       
! filtered output                                                       
     if(outf)then 
        CALL outsub(12,22,tf,eqf,angf,ngf,ncof,norb,na,nvz,.true.) 
     endif
! unfiltered output and periodic operations                             
     if(mod(iout,nsamp).eq.0)then 
        ioucou=ioucou+1 
        if(mod(iprqua,2).eq.0)then 
           CALL outsub(11,21,t,eq,angl,ng,ncof,norb,na,nvz,.true.) 
        endif
!  periodic operations: dump and renormalisation                        
!                                                                       
!  dump                                                                 
        if(mod(ioucou,idump).eq.0)then 
           call outsub(13,14,t,eq,angl,ng,ncof,norb,na,nvz,.false.) 
           if(nvz.ne.0)then 
              DO j=1,nvz 
                 write(14,*)(yv(i,j),i=1,ndim) 
                 write(14,*)(yv(i+ndim,j),i=1,ndim) 
              ENDDO
           endif
! *********************************************************             
!  renormalisation                                                      
           reno=.false. 
           DO 50 j=1,nvz 
              if(vnor(j).gt.v1)then 
                 reno=.true. 
!  renormalisation of  variation vector; orbit is unchanged             
                 renorm=v0/vnor(j) 
                 vslog(j)=vslog(j)-dlog(renorm) 
                 jj=3*(nbod-1+na+j)-3 
                 DO i=1,ndim 
                   y2(nvar2+jj+i)=y2(nvar2+jj+i)*renorm 
                   y2(i+jj)=y2(i+jj)*renorm 
                ENDDO
!                write(*,301)j,t,vnor(j),gamma(j)                       
                write(9,301)j,t,vnor(j),gamma(j) 
301             format(' renormalisation var.eq.no ',i4/               &
     &             ' time=',f18.6,' norm=',1p,d15.5,                    &
     &             ' gamma=',d15.5)                                     
             endif
50        ENDDO
          if(reno)then 
!  the starter must be used again                                       
             t1=t 
             DO j=1,nvar 
                y1(j)=y2(j) 
             ENDDO
             nfl=0 
          endif
       endif
    endif
10 ENDDO
!                                                                       
!   end main loop                                                       
!                                                                       
! *********************************************************             
!   end                                                                 
 write(9,199)t 
199 format(' end orbit8v, time = ',f18.4) 
END PROGRAM orbit9

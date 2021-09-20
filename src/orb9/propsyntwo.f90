! =========================================                             
!  {\bf propsynt} Proper elements by synthetic method                   
!  Pisa, December 1996 vers. 1.0                                        
!  Pisa, October 1998 vers. 1.1                                         
!  Pisa, June 2000 revision for performance and larger number of asteroi
! =========================================                             
PROGRAM propsynt 
  USE synthcomp
  USE name_rules, ONLY: name_len
  IMPLICIT NONE 
! max dimensions 
!  INCLUDE 'pardim.h' 
  INTEGER, PARAMETER :: nastx=1100                                         
  INTEGER nop,ii,nbeg 
  DOUBLE PRECISION gp(nforcx),sn(nforcx),aj 
  INTEGER klisg(nforcx),kliss(nforcx) 
!  elements, sampled and filtered                                       
  DOUBLE PRECISION elf(7,nbbx,ntxx),tf(ntxx),tp(ntxx),a(ntxx) 
!  workspaces                                                           
  DOUBLE PRECISION x(ntxx),y(ntxx),dy(ntxx),gam(ntxx),elt(7,ntxx),tt(ntxx) 
  INTEGER inflag
! for each proper element: 
! value, error, frequency, error in freq, phase,rms phase
  DOUBLE PRECISION pe0(3),dpe0(3),fe0(3),dfe0(3), ang0(3),rmsang0(3) 
! availability of elements                                              
  INTEGER jcon(ntxx,nastx) 
! then the same for each running box                             
  DOUBLE PRECISION pe(nbx,3),dpe(nbx,3),fe(nbx,3),dfe(nbx,3),ph(nbx,3),rmsph(nbx,3) 
! rms and max. excursion for elements, frequencies                      
! rms and max res. mean longitude                                       
  DOUBLE PRECISION sige(3),delte(3),sigf(3),deltf(3) 
  DOUBLE PRECISION devmax,devmad 
! lyapounov times etc.                                                  
  INTEGER ngfl,nt2,jg1,jg2,jg3 
  DOUBLE PRECISION carex,rms,cost,tlyap,carex1,carex2,dcarex,carex0,tm
!  file names                                                           
  CHARACTER*60 initia,filnam 
!  asteroid numbers                                                     
  CHARACTER*9 number(nastx) 
! loop indexes                                                          
  INTEGER n,ib,ij,i0,i1,i2,j,l,jj
! lengths, addresses                                                    
  INTEGER naf,len,nb,ntf,iwri 
! ================================================                      
!   filtered data file name                                             
  WRITE(*,*)' initial part of file names?' 
  READ(*,*)initia 
!  flag for the running box
  WRITE(*,*)' running box flag 1 = 10000 recs;&
& 2 = 20000 recs outer belt;&
& 3 = 20000 recs inner belt, TNO'
  READ(*,517)inflag
517 FORMAT(i1)
  CALL parprt(inflag,ntx,nbb,nforc,ndap,nshi)
! ==================================================                    
!   main output files                                                   
! ==================================================                    
  call cmpstr(initia,len) 
  filnam=initia(1:len)//'prs.out' 
  open(9,file=filnam,status='unknown') 
  filnam=initia(1:len)//'prs.pro' 
  open(10,file=filnam,status='unknown') 
  filnam=initia(1:len)//'prs.sig' 
  open(11,file=filnam,status='unknown') 
  filnam=initia(1:len)//'prs.del' 
  open(12,file=filnam,status='unknown') 
  filnam=initia(1:len)//'prs.lce' 
  open(13,file=filnam,status='unknown') 
  filnam=initia(1:len)//'prs.ang' 
  open(23,file=filnam,status='unknown')
  filnam=initia(1:len)//'prs.err' 
  open(22,file=filnam,status='unknown') 
! =====================================================                 
!   1: constant data                                                    
  call secth(gp,sn,aj,klisg,kliss) 
! =================================================                     
!   input data in aloop, rewinding the file                             
! ==================================================                    
  filnam=initia(1:len)//'ast.fil' 
  open(1,file=filnam,status='old') 
  filnam=initia(1:len)//'past.fil' 
  open(2,file=filnam,status='old') 
  nop=nastx/nbb+1 
  nbeg=1 
! write header in the auxiliary LCE file                                
  WRITE(13,*)'LCEs computed for backward and forward integrations,  &
     &and their difference. Values multiplied with 1.d6'                
! =================================================================     
  DO 2 ii=1,nop 
     call  inptwoway(1,2,tf,elf,ntf,nbeg,naf,nastx,number,     &
     &        elt,tt,jcon,ngfl,nt2)                                     
     WRITE(*,*)' nbeg,naf ',nbeg,naf 
! file has been all read                                                
     IF(naf.eq.0) GOTO 3 
! find central time
     tm=SUM(tf(1:ntf))/ntf
! initialise running box specifications                                 
     IF(ii.eq.1)THEN 
! ===================================================                   
! selection of running boxes                                            
! (5000/498 for 10000 records; 4000/1998 for 20000 records)             
! is done in parprt.h                                                   
! ===================================================                   
        WRITE(*,179)tf(1),tf(ntf),ntf 
        WRITE(9,179)tf(1),tf(ntf),ntf 
179     FORMAT(' data from t=',f14.4,' to t=',f14.4,' no datapoints=',i5) 
        WRITE(9,*)' together ',ndap,'  step',nshi 
        WRITE(*,*)' together ',ndap,'  step',nshi 
     ENDIF
! =====================================================                 
!  loop on the na asteroids                                             
     DO 1 n=1,naf 
! output in prs.out file                                                
        iwri=1 
        WRITE(9,*)' =====================================' 
        WRITE(9,*)' asteroid ',number(n) 
        DO jj=1,ntf 
           IF(jcon(jj,n+nbeg-1).eq.0)THEN 
              WRITE(9,*)'No proper elements computed; error in input' 
              WRITE(*,*)'Skipping asteroid ',number(n),jj,n+nbeg-1 
              GOTO 888 
           ENDIF
        ENDDO
!c =====================================================                
!   compute proper semimajor axis and mean motion                       
        WRITE(9,*)' =========== semimajor axis ' 
        DO j=1,ntf 
           x(j)=elf(6,n,j) 
           y(j)=elf(1,n,j) 
        ENDDO
        call dosemim(x,y,tf,ntf,pe0(1),dpe0(1),fe0(1),dfe0(1),          &
     &      ang0(1),rmsang0(1),devmax,iwri)
! =====================================================                 
!   compute proper eccentricity from all data                           
        WRITE(9,*)' =========== eccentricity ' 
        DO j=1,ntf 
          x(j)=elf(3,n,j) 
          y(j)=elf(2,n,j) 
       ENDDO
       call doinc(x,y,tf,ntf,gp,klisg,pe0(2),dpe0(2),                  &
     &             fe0(2),dfe0(2),ang0(2),rmsang0(2),iwri)
! =====================================================                 
!   compute proper inclination from all data                            
       WRITE(9,*)' =========== inclination ' 
       DO j=1,ntf 
          x(j)=elf(5,n,j) 
          y(j)=elf(4,n,j) 
       ENDDO
       iwri=1 
       call doinc(x,y,tf,ntf,sn,kliss,pe0(3),dpe0(3),                  &
     &             fe0(3),dfe0(3),ang0(3),rmsang0(3),iwri)
! conversion from tg(I/2) to sinI                                       
       dpe0(3)=2.d0*dpe0(3)* cos(2.d0*atan(pe0(3))) 
       pe0(3)=sin(2.d0*atan(pe0(3))) 
! =====================================================                 
!   running box test for accuracy; ib=box index                         
       ib=0 
       iwri=0 
       DO 50 ij=ndap,ntf,nshi 
          i0=ij-ndap 
          i1=ij-ndap+1 
          i2=ij 
          ib=ib+1 
          IF(n.eq.1)WRITE(*,*)i1,i2 
!   compute proper semimaj.axis, eccentricity and inclination from box d
          DO j=1,ndap 
             x(j)=elf(6,n,j+i0) 
             y(j)=elf(1,n,j+i0) 
          ENDDO 
          call dosemim(x,y,tf(i1),ndap,pe(ib,1),dpe(ib,1),fe(ib,1),   &
     &            dfe(ib,1),ph(ib,1),rmsph(ib,1),devmad,iwri)
          IF(pe(ib,1).gt.9.9999999d0)pe(ib,1)=9.9999999d0 
!                                                                       
          DO j=1,ndap 
             x(j)=elf(3,n,j+i0) 
             y(j)=elf(2,n,j+i0) 
          ENDDO 
          call doinc(x,y,tf(i1),ndap,gp,klisg,pe(ib,2),dpe(ib,2),fe(ib,2), &
     &         dfe(ib,2),ph(ib,2),rmsph(ib,2),iwri)
!                                                                       
          DO j=1,ndap 
             x(j)=elf(5,n,j+i0) 
             y(j)=elf(4,n,j+i0) 
          ENDDO 
          call doinc(x,y,tf(i1),ndap,sn,kliss,pe(ib,3),dpe(ib,3),fe(ib,3), &
     &         dfe(ib,3),ph(ib,3),rmsph(ib,3),iwri)                    
! conversion from tg(I/2) to sinI                                       
          pe(ib,3)=sin(2.d0*atan(pe(ib,3))) 
          dpe(ib,3)=2.d0*dpe(ib,3)* cos(2.d0*atan(pe(ib,3))) 
50     ENDDO
       nb=ib 
!   compute max deviation, standard deviation                           
       WRITE(9,*)' summary proper elements:' 
       DO l=1,3 
          sige(l)=sigmf(pe(1,l),pe0(l),nb) 
          delte(l)=MAX(MAXVAL(pe(1:nb,l)),pe0(l))-MIN(MINVAL(pe(1:nb,l)),pe0(l))
          sigf(l)=sigmf(fe(1,l),fe0(l),nb) 
          deltf(l)=MAX(MAXVAL(fe(1:nb,l)),fe0(l))-MIN(MINVAL(fe(1:nb,l)),fe0(l))
! security feature: the average must be close to the value on all the in
          WRITE(9,250)pe0(l),sige(l),delte(l),fe0(l),sigf(l),deltf(l) 
  250     FORMAT(f10.7,1x,f9.7,1x,f9.7,1x,f9.5,1x,f7.5,1x,f7.5) 
       ENDDO
! ==========================================================            
!   2d: Lyapounov characteristic exponents                              
! ==========================================================            
257    if(ngfl.eq.1)then 
          jg1=1 
          jg2=ntf-nt2 
          DO j=jg1,jg2 
             gam(jg2-j+1)=elf(7,n,j) 
             tp(jg2-j+1)=tf(j) 
!             if(n.eq.1)then 
!             endif
             IF(gam(jg2-j+1).eq.0.d0)THEN 
                carex=999.99d-6 
                rms=0.d0 
                cost=0.d0 
                GOTO 777 
             ENDIF
          ENDDO
       elseif(ngfl.eq.2)then 
          jg1=ntf-nt2+1 
          jg2=ntf 
          DO j=jg1,jg2 
             gam(j-jg1+1)=elf(7,n,j) 
             tp(j-jg1+1)=tf(j) 
             IF(gam(j-jg1+1).eq.0.d0)THEN 
                carex=999.99d-6 
                rms=0.d0 
                cost=0.d0 
                GOTO 777 
             ENDIF
          ENDDO
       else 
          WRITE(*,*)'Erroneous range for LCE' 
          stop 
       endif
       CALL linfi3(tp,gam,carex,rms,cost,dy,jg2-jg1+1) 
777    tlyap=1.d0/carex 
       WRITE(9,150)carex,rms,cost,tlyap 
150    FORMAT(' ======== max. LCE'/                                    &
     &    1p,d12.4,' rms.res',d12.4,' cost',d12.4,' tlyap',d12.4)       
       IF(ngfl.eq.1)THEN 
          if(carex.eq.999.99d-6)then 
             carex1=carex 
          else 
             carex1=-carex 
          endif
          if(carex1.lt.0.d0)carex1=0.d0 
          ngfl=2 
          GOTO 257 
       ELSEIF(ngfl.eq.2)THEN 
          carex2=carex 
          if(carex2.lt.0.d0)carex2=0.d0 
          ngfl=1 
       ENDIF
       IF(carex1.eq.999.99d-6.or.carex2.eq.999.99d-6)THEN 
          carex0=999.99d-6 
          dcarex=999.99d-6 
       ELSE 
          carex0=carex2 
          dcarex=abs(carex2-carex1) 
          if(carex0.lt.10.d-6)then 
             dcarex=10.d-6 
          elseif(dcarex.lt.10.d-6)then 
             dcarex=10.d-6 
          endif
       ENDIF
!       WRITE(12,151)number(n),carex,tlyap,clib(n)                      
! 151   FORMAT(1x,a4,3x,1p,d12.4,3x,d12.4,1x,0p,f10.6)                  
! ===========================================================           
! write summary                                                         
! ===========================================================   
! LCE data (only for twoway)        
       WRITE(13,113)number(n),carex1*1.d6,carex2*1.d6,dcarex*1.d6 
113    FORMAT(a9,3(1x,2f7.2)) 
! proper elements 
       WRITE(10,110)number(n),pe0,fe0,carex0*1.d6 
110    FORMAT(a9,1x,f11.7,1x,f10.7,1x,f9.7,1x,                         &
     &               f11.6,1x,f12.6,1x,f12.6,1x,f7.2)                   
! RMS of runnign boxes
       WRITE(11,111)number(n),sige,sigf,rmsang0(1)
111    FORMAT(a9,1x,f11.8,1x,f10.7,1x,f9.7,1x,                         &
     &               f9.6,1x,f10.6,1x,f9.6,1x,f12.4) 
! excursion of running boxes
       WRITE(12,112)number(n),delte,deltf,devmax 
112    FORMAT(a9,1x,f11.8,1x,f11.7,1x,f9.7,1x,                         &
     &               f9.6,1x,f10.6,1x,f9.6,1x,f12.4) 
 ! convert angles to [0,360]
       ang0(1)=ang0(1)-360*floor(ang0(1)/360)
       ang0(2)=ang0(2)-360*floor(ang0(2)/360)
       ang0(3)=ang0(3)-360*floor(ang0(3)/360)
       WRITE(23,123)number(n),ang0,rmsang0 
123    FORMAT(A,3(1x,f10.6),3(1x,f12.4))
       ! ==================================================                    
!  end loop on asteroids                                                
! ==================================================                    
888    CONTINUE 
1   ENDDO
    nbeg=nbeg+nbb 
    rewind(1) 
    rewind(2) 
! end loop on file rewind                                               
2 ENDDO
3 CONTINUE 
 CLOSE(1) 
 CLOSE(2) 
 CLOSE(9) 
 CLOSE(10) 
 CLOSE(11) 
 CLOSE(12) 
 CLOSE(22) 
 CLOSE(23)
END PROGRAM propsynt

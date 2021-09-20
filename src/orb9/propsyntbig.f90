! =========================================                             
!  {\bf propsynt} Proper elements by synthetic method                   
!  Pisa, December 1996 vers. 1.0                                        
!  Pisa, October 1998 vers. 1.1                                         
!  Pisa, June 2000 revision for performance and larger number of asteroi
!  Pisa, December 2000 revision to handle errors in input               
! =========================================                             
PROGRAM propsynt 
  USE synthcomp
  USE name_rules, ONLY: name_len
  IMPLICIT NONE 
! max dimensions                                                        
  INTEGER, PARAMETER :: nastx=1100
! forced terms                                                          
  DOUBLE PRECISION gp(nforcx),sn(nforcx),aj 
  INTEGER klisg(nforcx),kliss(nforcx) 
!  elements, sampled and filtered                                       
  DOUBLE PRECISION elf(7,nbbx,ntxx),tf(ntxx),a(ntxx) 
!  workspaces                                                           
  DOUBLE PRECISION x(ntxx),y(ntxx),dy(ntxx),gam(ntxx) 
! selection flag
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
! lyapounov times, central time, etc
  DOUBLE PRECISION carex,rms,cost,tlyap,tm 
!  file names                                                           
  character*60 initia,filnam 
!  asteroid numbers                                                     
  character*(name_len) :: number(nastx) 
! loop indexes                                                          
  INTEGER n,ib,ij,i0,i1,i2,j,l,jj,ii 
! lengths, addresses                                                    
  INTEGER naf,len,nb,ntf,iwri,nop,nbeg 
! secular resonance: flag, libration phase deg (proper ecc. is amplitude)
  INTEGER resflag ! 0=non-res., 1=res, forced, -1=non-res, forced
  DOUBLE PRECISION phalib, phalibb(nbx), sigph, sigdel! for running box test of phase
! ===================================================
! interactive input
! ===================================================
!
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
  nop=nastx/nbb+1 
  nbeg=1 
  DO 2 ii=1,nop 
     call  inptro(1,tf,elf,ntf,nbeg,naf,nastx,number,jcon) 
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
! is done in subroutine parprt
! ===================================================
        WRITE(*,179)tf(1),tf(ntf),ntf 
        WRITE(9,179)tf(1),tf(ntf),ntf 
179     FORMAT(' data from t=',f14.4,' to t=',f14.4,' no datapoints=',i5) 
        WRITE(9,'(a10,1x,i10,1x,a6,1x,i10)')' together ',ndap,'  step',nshi 
        WRITE(*,*)' together ',ndap,'  step',nshi 
     ENDIF
! =====================================================                 
!  loop on the na asteroids                                             
     DO 1 n=1,naf 
! output in prs.out file                                                
        iwri=1 
        WRITE(9,*)' ==================================' 
        WRITE(9,*)' asteroid ',number(n) 
        DO jj=1,ntf 
           IF(jcon(jj,n+nbeg-1).eq.0)THEN 
              WRITE(9,*)'No proper elements computed; error in input' 
              WRITE(*,*)'Skipping asteroid ',number(n),jj,n+nbeg-1 
              GOTO 888 
           ENDIF
        ENDDO
! ===================================================== 
!   compute proper semimajor axis and mean motion, from all data
        WRITE(9,*)' =========== semimajor axis ' 
        DO j=1,ntf 
          x(j)=elf(6,n,j) 
          y(j)=elf(1,n,j) 
        ENDDO 
        call dosemim(x,y,tf,ntf,pe0(1),dpe0(1),fe0(1),dfe0(1),          &
     &     ang0(1),rmsang0(1),devmax,iwri)
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
!   compute proper eccentricity from all data                           
        WRITE(9,*)' =========== eccentricity ' 
        DO j=1,ntf 
          x(j)=elf(3,n,j) 
          y(j)=elf(2,n,j) 
        ENDDO
        IF(pe0(1).gt.2.5d0.and.pe0(1).lt.2.7d0)THEN
           resflag=0 ! non-resonant, but can be changed to resonance after test
        ELSE
           resflag=-1 ! force ignore secular resonance
        ENDIF 
        CALL doecc(x,y,tf,ntf,gp,klisg,pe0(2),dpe0(2),                  &
     &        fe0(2),dfe0(2),ang0(2),rmsang0(2),iwri,resflag,phalib,pe0(3),pe0(1))
!        WRITE(*,*) 0, pe0(2), fe0(2)
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
!          IF(pe(ib,1).gt.9.9999999d0)pe(ib,1)=9.9999999d0              
!
! proper inclination
          DO j=1,ndap 
            x(j)=elf(5,n,j+i0) 
            y(j)=elf(4,n,j+i0) 
          ENDDO 
          call doinc(x,y,tf(i1),ndap,sn,kliss,pe(ib,3),                 &
     &            dpe(ib,3),fe(ib,3),dfe(ib,3),ph(ib,3),rmsph(ib,3),iwri)
! conversion from tg(I/2) to sinI                                       
          pe(ib,3)=sin(2.d0*atan(pe(ib,3))) 
          dpe(ib,3)=2.d0*dpe(ib,3)* cos(2.d0*atan(pe(ib,3))) 
! to avoid that different running boxes have different resonance status
          IF(resflag.eq.0) resflag=-1
! proper eccentricity 
          DO j=1,ndap 
            x(j)=elf(3,n,j+i0) 
            y(j)=elf(2,n,j+i0) 
          ENDDO
          CALL doecc(x,y,tf(i1),ndap,gp,klisg,pe(ib,2),dpe(ib,2),fe(ib,2),   &
     &        dfe(ib,2),ph(ib,2),rmsph(ib,2),iwri,resflag,phalibb(ib),pe(ib,3),pe(ib,1))
!          WRITE(*,*) ib, pe(ib,2),fe(ib,2)

   50   ENDDO 
        nb=ib 
!   compute max deviation, standard deviation                           
        WRITE(9,*)' summary proper elements:' 
        DO l=1,3 
           sige(l)=sigmf(pe(1,l),pe0(l),nb) 
           delte(l)=MAX(MAXVAL(pe(1:nb,l)),pe0(l))-MIN(MINVAL(pe(1:nb,l)),pe0(l))
           sigf(l)=sigmf(fe(1,l),fe0(l),nb) 
           deltf(l)=MAX(MAXVAL(fe(1:nb,l)),fe0(l))-MIN(MINVAL(fe(1:nb,l)),fe0(l))
! security feature: the average must be close to the value on all the interval
           WRITE(9,250)pe0(l),sige(l),delte(l),fe0(l),sigf(l),deltf(l) 
250        FORMAT(f10.7,1x,f9.7,1x,f9.7,1x,f9.5,1x,f7.5,1x,f7.5) 
        ENDDO 
        IF(resflag.eq.1)THEN
           sigph=sigmfprideg(phalibb(1),phalib,nb)
!           sigdel=MAX(MAXVAL(phalibb(1:nb)),phalib)-MIN(MINVAL(phalibb(1:nb)),phalib)
! security feature: the average must be close to the value on all the interval
           WRITE(9,251)phalib, sigph 
251        FORMAT(f12.8,1x,f11.7) 
        ELSE
           sigph=0.d0
           sigdel=0.d0
        ENDIF
! ==========================================================            
!   2d: Lyapounov characteristic exponents                              
! ==========================================================            
        DO j=1,ntf 
          gam(j)=elf(7,n,j) 
          IF(gam(j).eq.0.d0)THEN 
             carex=999.99d-6 
             rms=0.d0 
             cost=0.d0 
             GOTO 777 
          ENDIF 
        ENDDO 
        CALL linfi3(tf,gam,carex,rms,cost,dy,ntf) 
  777   tlyap=1.d0/carex 
        WRITE(9,150)carex,rms,cost,tlyap 
150     FORMAT(' ======== max. LCE'/                                    &
     &    1p,d12.4,' rms.res',d12.4,' cost',d12.4,' tlyap',d12.4)       
                                                                        
!       WRITE(12,151)number(n),carex,tlyap,clib(n)                      
! 151   FORMAT(1x,a4,3x,1p,d12.4,3x,d12.4,1x,0p,f10.6)                  
! ===========================================================           
! write summary                                                         
! ===========================================================           
! proper elements        
        WRITE(10,110)number(n),pe0,fe0,carex*1.d6 
110     FORMAT(a9,1x,f11.7,1x,f10.7,1x,f9.7,1x,                         &
     &               f11.6,1x,f12.6,1x,f12.6,1x,f7.2)                   
! RMS of runnign boxes
        WRITE(11,111)number(n),sige,sigf,rmsang0(1) 
111     FORMAT(a9,1x,f11.8,1x,f10.7,1x,f9.7,1x,                         &
     &               f9.6,1x,f10.6,1x,f9.6,1x,f12.4) 
! excursion of running boxes                   
        WRITE(12,112)number(n),delte,deltf,devmax 
112     FORMAT(a9,1x,f11.8,1x,f11.7,1x,f9.7,1x,                         &
     &               f9.6,1x,f10.6,1x,f9.6,1x,f12.4) 
! convert angles to [0,360]
        ang0(1)=ang0(1)-360*floor(ang0(1)/360)
        ang0(2)=ang0(2)-360*floor(ang0(2)/360)
        ang0(3)=ang0(3)-360*floor(ang0(3)/360)
        phalib=phalib-360*floor(phalib/360)
        WRITE(23,123)number(n),ang0,phalib,rmsang0, sigph
123     FORMAT(A,4(1x,f12.6),4(1x,f12.4))
! ==================================================                    
!  end loop on asteroids                                                
! ==================================================                    
888     CONTINUE 
1    ENDDO
      nbeg=nbeg+nbb 
      rewind(1) 
! end loop on file rewind                                               
2  ENDDO
3  CONTINUE 
   CLOSE (1) 
   CLOSE(9) 
   CLOSE(10) 
   CLOSE(11) 
   CLOSE(12) 
   CLOSE(22) 
   CLOSE(23)
 END PROGRAM propsynt

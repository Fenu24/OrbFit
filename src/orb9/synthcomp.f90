MODULE synthcomp
  USE fund_const
! computation of synthetic proper elements 
  IMPLICIT NONE
  PRIVATE

! public routines
  PUBLIC secth, doecc,doinc, dosemim, sigmf, sigmfprideg, sigma, parprt

! shared data

  INTEGER, PUBLIC :: ntx,nbb,nforc,ndap,nshi

  INTEGER, PARAMETER, PUBLIC :: nforcx=9 ! max no forced terms
! array size for all the data read at once
  INTEGER, PARAMETER, PUBLIC :: ntxx=20001 ! max no of records
  INTEGER, PARAMETER, PUBLIC :: nbbx=600 ! max no of bodies at once
!  max no. running boxes
  INTEGER, PARAMETER, PUBLIC :: nbx=100

! ====================================================
! parameters for propert
! array size for all the data read at once
! for the jobs with 20000 points
!  INTEGER, PARAMETER :: ntx=20001
!  INTEGER, PARAMETER :: nbb=355
! for the jobs with 10000 points
!  INTEGER, PARAMETER :: ntx=10001
!  INTEGER, PARAMETER :: nbb=600
! number of forced terms
!  INTEGER, PARAMETER :: nforc=9
! selection of running boxes 
! (5000/498 for 10000 records; 4000/1998 for 20000 records
!  in the outer belt, 10000/997 for 20000 records for the 
!  inner belt twoway, and TNO)
!  INTEGER, PARAMETER :: ndap=4000
!  INTEGER, PARAMETER :: nshi=1998
!  INTEGER, PARAMETER :: ndap=5000
!  INTEGER, PARAMETER :: nshi=498
!  INTEGER, PARAMETER :: ndap=10000
!  INTEGER, PARAMETER :: nshi=997

CONTAINS
! ==================================================================
! selection of parameters for running box
! ==================================================================
  SUBROUTINE parprt(inflag,ntx,nbb,nforc,ndap,nshi)
    INTEGER, INTENT(IN) :: inflag
    INTEGER, INTENT(OUT) :: ntx,nbb,nforc,ndap,nshi
! END INTERFACE
    ! number of forced terms
    nforc=9   
    IF(nforc.gt.nforcx)THEN
       WRITE(*,*)' parprt: too many forced terms, nforc=',nforc,   &
            &       'max was nforcx=',nforcx
       STOP 
    ENDIF
    IF(inflag.eq.1)THEN
! for the jobs with 10000 points (outer belt)
       ntx=10001
       IF(ntx.gt.ntxx)THEN
          WRITE(*,*)' parprt: too many records, ntx=',ntx,    &
               &         'max was ntxx=',ntxx
          STOP 
       ENDIF
       nbb=600
       IF(nbb.gt.nbbx)THEN
          WRITE(*,*)' parprt: too many bodies, nbb=',nbb,    &
               &         'max was nbbx=',nbbx
          STOP 
       ENDIF
       ndap=5000
       nshi=498
    ELSEIF(inflag.eq.2)THEN
       ! for the jobs with 20000 points (outer belt, extended integration)
       ntx=20001
       IF(ntx.gt.ntxx)THEN
          WRITE(*,*)' parprt: too many records, ntx=',ntx,    &
               &         'max was ntxx=',ntxx
          STOP 
       ENDIF
       nbb=355
       IF(nbb.gt.nbbx)THEN
          WRITE(*,*)' parprt: too many bodies, nbb=',nbb,    &
               &         'max was nbbx=',nbbx
          STOP 
       ENDIF
       ndap=4000
       nshi=1998
    ELSEIF(inflag.eq.3)THEN
       ! for the jobs with 20000 points (inner belt twoway, TNO)
       ntx=20001
       IF(ntx.gt.ntxx)THEN
          WRITE(*,*)' parprt: too many records, ntx=',ntx,    &
               &         'max was ntxx=',ntxx
          STOP 
       ENDIF
       nbb=355
       IF(nbb.gt.nbbx)THEN
          WRITE(*,*)' parprt: too many bodies, nbb=',nbb,    &
               &         'max was nbbx=',nbbx
          STOP 
       ENDIF
       ndap=10000  
       nshi=997
    ELSE
       WRITE(*,*)' parprt: option not known, inflag=', inflag
       STOP
    ENDIF
  END SUBROUTINE parprt
! ==================================================================
! proper semimajor axis
! ==================================================================
  SUBROUTINE dosemim(x,y,tf,ntf,proa,dproa,fra,dfra,ph,rmsph,devmax,iwri)
    INTEGER, INTENT(IN) ::  ntf,iwri
    DOUBLE PRECISION, INTENT(IN) ::  x(ntf),y(ntf),tf(ntf)
    DOUBLE PRECISION, INTENT(OUT) ::  proa,dproa,fra,dfra,ph,rmsph,devmax
! end interface
    DOUBLE PRECISION rate,dy(ntx),per
    INTEGER j 
! ===========================================                           
    call linfi3(tf,x,fra,rmsph,ph,dy,ntf) 
    dfra=rmsph/(MAXVAL(tf(1:ntf))-MINVAL(tf(1:ntf))) 
    per=dpig/fra 
    fra=fra*degrad 
    dfra=dfra*degrad 
    ph=ph*degrad 
    rmsph=rmsph*degrad 
    devmax=(MAXVAL(dy(1:ntf))-MINVAL(dy(1:ntf)))*degrad 
    IF(iwri.eq.1)THEN 
       write(9,120)fra,per,rmsph,ph 
  120    format(' =====================Mean long:'/                     &
     &  ' freq ',f12.8,' deg/yr; per ',f13.7,' yr',                      &
     &  ' rms ',f12.6,' deg;  ph ',f10.4,' deg')                        
    ENDIF
!======================================================                 
!   3a: semimajor axis                                                  
! =====================================================                 
!  for each asteroid, compute average mean a                            
!  compute on all                                                       
    proa= SUM(y(1:ntf))/ntf
    dproa=sigma(y,ntf) 
    IF(iwri.eq.1)THEN 
       write(9,111)proa,dproa 
111    format(' ======== a : ',f12.9,' sig=',f10.9) 
    ENDIF
  END SUBROUTINE dosemim
! =====================================================
! proper eccentricity
! =====================================================                 
  SUBROUTINE doecc(x,y,tf,ntf,sn,klis,pe,dpe,fre,dfre,ang,rms,iwri,resflag,phalib,sini,apr) 
    USE fir_filters
    INTEGER, INTENT(IN) ::  ntf                       ! vector length
    DOUBLE PRECISION, INTENT(IN) :: tf(ntf)           ! time
    DOUBLE PRECISION, INTENT(IN) :: sini,apr          ! proper sinI as just computed by doinc, proper a
    DOUBLE PRECISION, INTENT(INOUT) :: x(ntf),y(ntf)  ! data h,k
    DOUBLE PRECISION, INTENT(OUT) :: pe,fre           ! proper element, frequency
    DOUBLE PRECISION, INTENT(OUT) :: dpe,dfre,rms     !auxiliary output for accuracy
! for a specific secular resonance
    INTEGER, INTENT(INOUT) :: resflag   ! 1 if resonant, 0 if not, -1 forced non-res
                       ! if 1 in input then resonance is forced
    DOUBLE PRECISION, INTENT(OUT) :: phalib,ang! initial phase of libration, of perihelion (deg)
! workspace   
    DOUBLE PRECISION  x1(ntf),y1(ntf)
    DOUBLE PRECISION omeg(ntx) 
!  planetary theory, and forced terms                                   
    DOUBLE PRECISION fe(nforc),feph(nforc),fesp(nforc) 
    DOUBLE PRECISION feco(nforc),fed(nforc),sn(nforc) 
    INTEGER klis(nforc) 
    INTEGER k,iwri 
    DOUBLE PRECISION rate,peri,cost,sig,sp,ph
    DOUBLE PRECISION princ
! resonant proper elements
    DOUBLE PRECISION ecc(ntx),eccsm(ntx),tsm(ntx), tper(ntx),amp(ntx),spv(ntx), fas(ntx), costv(ntx)
    DOUBLE PRECISION diffre,fstep, hlen, xf, yf(1), frelib, eccmean
    DOUBLE PRECISION dper,pmax,pmin, a,f,d0,d1,d2, spmax,amax,tpermax, per
    LOGICAL first, outf
    DATA first /.true./
    INTEGER i, j, nskip, ndat, nn, nfre, imax
! =====================================================                 
!   2: determination of frequencies by fit to arguments                 
!   all the linear forced terms are removed first                       
! =====================================================                 
    CALL forced(x,y,tf,ntf,sn,klis,fe,feph,fesp,feco,fed) 
    IF(iwri.eq.1)THEN 
       DO  k=5,9 
          IF(klis(k).ne.0)THEN 
             write(9,127)k,fe(k),feph(k),fesp(k)*100,feco(k),fed(k) 
127          format(' forced ',i2,' amp ',1p,d12.5,0p,' ph ',f9.4,     &
     &        ' S% ',f8.4,' cost ',1p,d12.4,' unc. ',d12.4)             
          ENDIF
       ENDDO
    ENDIF
! =====================================================                 
!   now compute argument Omega from q,p, varpi from k,h                 
! =====================================================                 
    CALL argum(x,y,tf,ntf,omeg,rate,peri,rms,ang)
    ang=princ(ang)*degrad
    rms=rms*degrad 
    fre=rate*degrad*3.6d3 
    dfre=rms/(MAXVAL(tf(1:ntf))-MINVAL(tf(1:ntf)))*3.6d3 
    IF(iwri.eq.1)THEN 
       write(9,128)fre,peri,rms,ang 
128    format(' Argument:',                                           &
     &  ' freq ',1p,d13.6,0p,'"/yr, per '                             &
     &,f16.6,' yr, rms ',1p,d12.4,0p,'deg  ph ',f9.4,' deg')
    ENDIF
! =====================================================
! test for secular resonance (hacked for Astraea, to be used only for 2.5<a<2.7)
! =====================================================
    diffre=fre+sn(5)-2.d0*sn(6)
    IF((abs(diffre).lt.0.065d0.or.resflag.eq.1).and..not.resflag.eq.-1)THEN
       resflag=1
! compute eccentricity
       ecc=0.d0
       DO i=1,ntf
         ecc(i)=sqrt(x(i)**2+y(i)**2)        
       ENDDO
! apply filter (see giffv.f90 lines 263-291), use filter.100
       IF(first)THEN
! filter coefficients input
          CALL filope(100,1)
          first=.false.
       ENDIF    
!  assumed constant step
       fstep= x(ntf/2+1) - x(ntf/2)
! filter half length
       hlen=fstep*(nfil-1)/2 
!  digital low frequency pass filter applied to ecc vector
       nskip = 0 
       j = 0 
       DO i = 1, ntf 
          ndat=1
          CALL filter (tf(i), ecc(i:i), ndat, nskip, fstep, xf, yf(1:1), outf)
          IF (outf) then 
             j = j + 1 
             eccsm(j) = yf(1)
             tsm(j) = xf
          ENDIF
       ENDDO
       nn = j 
! compute spectrum, find maximum(?) of S at period pmax limited by??
! see giffv.f90 lines 174-186
       pmax=6d6; pmin=1.d6; nfre=2001
       dper = (pmax - pmin) / nfre 
       DO  j = 1, nfre
          per = pmin + dper * (j - 1) 
          CALL peri2 (tsm, eccsm, nn, per, sp, 0, d0, d1, d2)
          a = d1 * d1 + d2 * d2 
          a = sqrt (a) 
          tper (j) = per 
          amp (j) = a 
          spv(j) = sp
          fas(j) = atan2 (d1, d2) 
          costv(j) = + d0 
       ENDDO
! then do trivial search for maximum on these points 
       spmax=0.d0
       imax=0
       DO j=1, nfre
          IF(spv(j).gt.spmax)THEN
             spmax=spv(j)
             imax=j
          ENDIF
       ENDDO
       amax=amp(imax)
       tpermax=tper(imax)
       frelib=dpig/tpermax*degrad*3.6d3 
       phalib=princ(fas(imax))*degrad
       eccmean=costv(imax)
! compute pe as amplitude at maximum Spectral Density
       pe=amax
       pe=pe+1.d0 ! this to separate from ordinary proper elements
       IF(iwri.eq.1)THEN 
          write(9,130)pe,spmax*100,tpermax,frelib,phalib,diffre,eccmean
130       format(' Proper ',1p,d13.6,0p,' S% ',f8.4,' period.lib. ',f12.3,&
    & ' freq.lib ',f8.5,' initial phase ',f9.4,' diff. g from resonant ',f9.5, &
    & ' mean ecc. ',f6.4)
       ENDIF
! compute freqency as dpig/pmax, convert in arcsec/y
!       freres=(3.6d2*3.6D3)/pmax 
       IF(iwri.eq.1)THEN 
! write something
       ENDIF
    ELSE
       resflag=0
       phalib=0.0d0
!======================================================                 
!   3: find the proper amplitude and phase                              
! =====================================================
       x1=x;y1=y                 
       CALL prop(omeg,x1,y1,ntf,dpig,sp,pe,dpe,ph,cost,sig) 
       ph=ph*degrad
       IF(iwri.eq.1)THEN 
          write(9,129)pe,ph,sp*100,cost,dpe,sig 
129       format(' Proper ',1p,d13.6,0p,' arg ',f9.4,    &
     &        'deg, S% ',f8.4,' const ',1p,d12.4,' unc ',d9.2,   &
     &        ' rms',d9.2)                                              
       ENDIF
! check for bad fit to argument, low proper eccentricity
       IF(rms.gt.180d0.and.pe.lt.0.03d0)THEN
! select Hansa region; in this rectangle of a, sinI there are no other families
          IF(sini.gt.0.35d0.and.sini.lt.0.4d0.and.apr.gt.2.5d0.and.apr.lt.2.735d0)THEN
! compute spectrum, find maximum(?) of S at period pmax limited by??
! see giffv.f90 lines 174-186
             pmax=7.d4; pmin=5.d4; nfre=1001
             dper = (pmax - pmin) / nfre 
             DO  j = 1, nfre
                per = pmin + dper * (j - 1) 
                CALL peri2 (tf, x, ntf, per, sp, 0, d0, d1, d2)
                a = d1 * d1 + d2 * d2 
                a = sqrt (a) 
                tper (j) = per 
                amp (j) = a 
                spv(j) = sp
                fas(j) = atan2 (d1, d2)
                costv(j) = + d0 
             ENDDO
! then do trivial search for maximum on these points 
             spmax=0.d0
             imax=0
             DO j=1, nfre
                IF(spv(j).gt.spmax)THEN
                   spmax=spv(j)
                   imax=j
                ENDIF
             ENDDO
! compute pe as amplitude at maximum Spectral Density, g as the corresp. frequency
             pe=amp(imax)
             tpermax=tper(imax)       
             fre=dpig/tpermax*degrad*3.6d3 
             ph=princ(fas(imax))*degrad
             cost=costv(imax)
             dpe=0.d0;dfre=0.d0
             IF(iwri.eq.1)THEN 
                write(9,179)pe,ph,spmax*100,cost,dpe,sig 
179             format(' Proper ',1p,d13.6,0p,' arg ',f9.4,    &
         &        'deg, S% ',f8.4,' const ',1p,d12.4,' unc ',d9.2,   &
         &        ' rms',d9.2)
             ENDIF
          ENDIF
       ENDIF
    ENDIF
! =====================================================
! check for bad fit:
!======================================================
    IF(iwri.eq.0)RETURN  ! not in running box test
    IF(rms.gt.180.d0)THEN
! diffre is the diff. frequency w.r. to g-2g6+g5, problems for abs(diffre)<0.4
!  pe is proper ecc which gives problems for pe<0.02 
       WRITE(9,199)rms, diffre, pe
 199   FORMAT(' bad fit to argument: rms=',f6.1,' diffre=',f9.3,' prop. e=',f6.4)
    ENDIF
! =====================================================                 
! check for possible secular resonances                                 
! =====================================================                 
    DO k=5,9 
       IF(klis(k).ne.0)THEN 
          IF(abs(fre-sn(k)).lt.2.d0.or.fe(k).ge.pe)THEN 
             WRITE(9,113)k,sn(k),fre,fe(k),pe 
113          FORMAT('Sec.res.? ',i2,1x,f8.4,1x,f8.4,    &
     &               ' ampl ',f8.6,1x,f8.6)                             
          ENDIF
       ENDIF
    ENDDO
  END SUBROUTINE doecc

! =====================================================
! proper inclination, or proper eccentricity
! =====================================================                 
  SUBROUTINE doinc(x,y,tf,ntf,sn,klis,pe,dpe,fre,dfre,ang,rms,iwri) 
    INTEGER, INTENT(IN) ::  ntf 
    DOUBLE PRECISION, INTENT(IN) :: tf(ntf)
    DOUBLE PRECISION, INTENT(INOUT) :: x(ntf),y(ntf)
    DOUBLE PRECISION, INTENT(OUT) :: pe,dpe,fre,dfre,rms,ang
! workspace                                                             
    DOUBLE PRECISION omeg(ntx) 
!  planetary theory, and forced terms                                   
    DOUBLE PRECISION fe(nforc),feph(nforc),fesp(nforc) 
    DOUBLE PRECISION feco(nforc),fed(nforc),sn(nforc) 
    INTEGER klis(nforc) 
    INTEGER k,iwri 
    DOUBLE PRECISION rate,peri,cost,sig,sp,ph
    DOUBLE PRECISION princ
! =====================================================                 
!   2: determination of frequencies by fit to arguments                 
!   all the linear forced terms are removed first                       
! =====================================================                 
    CALL forced(x,y,tf,ntf,sn,klis,fe,feph,fesp,feco,fed) 
    IF(iwri.eq.1)THEN 
       DO  k=5,9 
          IF(klis(k).ne.0)THEN 
             write(9,127)k,fe(k),feph(k),fesp(k)*100,feco(k),fed(k) 
127          format(' forced ',i2,' amp ',1p,d12.5,0p,' ph ',f9.4,     &
     &        ' S% ',f8.4,' cost ',1p,d12.4,' unc. ',d12.4)             
          ENDIF
       ENDDO
    ENDIF
! =====================================================                 
!   now compute argument Omega from q,p, varpi from k,h                 
! =====================================================                 
    CALL argum(x,y,tf,ntf,omeg,rate,peri,rms,ang)
    ang=princ(ang)*degrad
    rms=rms*degrad 
    fre=rate*degrad*3.6d3 
    dfre=rms/(MAXVAL(tf(1:ntf))-MINVAL(tf(1:ntf)))*3.6d3 
    IF(iwri.eq.1)THEN 
       write(9,128)fre,peri,rms,ang 
128    format(' Argument:',                                           &
     &  ' freq ',1p,d13.6,0p,'"/yr, per '                             &
     &,f16.6,' yr, rms ',1p,d12.4,0p,'deg  ph ',f9.4,' deg')
    ENDIF
!======================================================                 
!   3: find the proper amplitude and phase                              
! =====================================================                 
    CALL prop(omeg,x,y,ntf,dpig,sp,pe,dpe,ph,cost,sig) 
    ph=ph*degrad
    IF(iwri.eq.1)THEN 
       write(9,129)pe,ph,sp*100,cost,dpe,sig 
129    format(' Proper ',1p,d13.6,0p,' arg ',f9.4,    &
     &        'deg, S% ',f8.4,' const ',1p,d12.4,' unc ',d9.2,   &
     &        ' rms',d9.2)                                              
    ENDIF
! =====================================================                 
! check for possible secular resonances                                 
! =====================================================                 
!     If(iwri.eq.0)RETURN                                               
    DO k=5,9 
       IF(klis(k).ne.0)THEN 
          IF(abs(fre-sn(k)).lt.2.d0.or.fe(k).ge.pe)THEN 
             WRITE(9,113)k,sn(k),fre,fe(k),pe 
113          FORMAT('Sec.res.? ',i2,1x,f8.4,1x,f8.4,    &
     &               ' ampl ',f8.6,1x,f8.6)                             
          ENDIF
       ENDIF
    ENDDO
  END SUBROUTINE doinc
! ===================================================================== 
!  secth
! =====================================================================
!  Routine which provides the fundamental frequencies and the semimajor 
!  axis of Jupiter according to the LONGSTOP 1B synthetic theory.       
!  See Nobili et al.  A&A 210,313, 1989                                 
!  g(9) is used to accomodate the most important combination frequency  
!  UNITS: semimaj axis= AU; freq= arcsec/yr                             
  SUBROUTINE secth(g,s,aj,klisg,kliss) 
    DOUBLE PRECISION, INTENT(OUT) ::  g(9),s(9),aj 
    INTEGER, INTENT(OUT) ::  kliss(nforcx),klisg(nforcx) 
    INTEGER jp 
!  chose forced terms                                                   
    DO jp=1,4 
       klisg(jp)=0 
       kliss(jp)=0 
    ENDDO
! forced frequencies of perihelion: no Neptune                          
    klisg(5)=1 
    klisg(6)=1 
    klisg(7)=1 
    klisg(8)=0 
! including main degree 3 term                                          
    klisg(9)=0 
! forced frequencies of the nodes: no s5                                
    kliss(5)=0 
    kliss(6)=1 
    kliss(7)=1 
    kliss(8)=1 
    kliss(9)=0 
! fundamental frequencies                                               
! from orbit9d integration for -10 MY                                   
    g(5)=4.2574d0 
    g(6)=28.2456d0 
    g(7)=3.0878d0 
    g(8)=0.6719d0 
! main nonlinear frequency                                              
    g(9)=2.d0*g(6)-g(5) 
! nodes                                                                 
    s(6)=-26.3453d0 
    s(7)=-2.9938d0 
    s(8)=-0.6936d0 
! This has some problems...                                             
    s(5)=0.d0 
! Semimajor axes (au) derived from LONGSTOP 1B (average of filtered     
!  --periods  < 100,000 yr wiped out)                                   
    aj=5.2025696d0 
!    write(*,*)aj                                                      
  END SUBROUTINE secth
! =====================================================
! low level routines
! =====================================================                 
  SUBROUTINE prop(tf,x,y,ntf,per,sp,amp,damp,ph,cost,sig) 
    INTEGER, INTENT(IN) :: ntf 
    DOUBLE PRECISION, INTENT(IN) :: tf(ntf), per
    DOUBLE PRECISION, INTENT(INOUT) :: x(ntf),y(ntf)
    DOUBLE PRECISION, INTENT(OUT) :: sp,amp,damp,ph,cost,sig
! end interface
!  workspaces                                                           
!    DOUBLE PRECISION vc(ntx),vs(ntx),dy(ntx) 
!    DOUBLE PRECISION yy(ntx),xx(ntx)
! functions                                                             
    DOUBLE PRECISION princ 
    DOUBLE PRECISION d0,d1,d2,d0k,d1k,d2k,sigk,ampk 
    INTEGER itest 
!  remove proper mode                                                   
    itest=1 
    call peri2(tf,y,ntf,per,sp,itest,d0,d1,d2) 
    sig=sigma(y,ntf) 
    call peri2(tf,x,ntf,per,sp,itest,d0k,d1k,d2k) 
    sigk=sigma(x,ntf) 
    sig=max(sig,sigk) 
    amp=sqrt(d1*d1+d2*d2) 
    ampk=sqrt(d1k*d1k+d2k*d2k) 
    damp=ampk-amp 
    ph=princ(atan2(d1,d2)) 
    cost=d0 
  END SUBROUTINE prop
! ========================================================              
!   compute argument varpi from k,h and find frequency                  
  SUBROUTINE argum(x,y,tf,ntf,th,rate,per,rms,cost) 
    INTEGER, INTENT(IN) ::  ntf 
    DOUBLE PRECISION, INTENT(IN)  ::  x(ntf),y(ntf),tf(ntf)
    DOUBLE PRECISION, INTENT(OUT) ::  th(ntx),rate,per,rms,cost 
!  workspaces                                                           
    DOUBLE PRECISION dy(ntx),xm,ym 
! functions                                                             
    DOUBLE PRECISION princ
    INTEGER ng(ntx),j 
!  use as center the mean                                               
    xm=SUM(x(1:ntf))/ntf 
    ym=SUM(y(1:ntf))/ntf 
!  initialise number rev. at 0                                          
    ng(1)=0 
    th(1)=atan2(y(1),x(1)) 
    DO j=2,ntf 
       th(j)=princ(atan2(y(j)-ym,x(j)-xm)) 
       if(th(j).gt.th(j-1)+pig)then 
          ng(j)=ng(j-1)-1 
       elseif(th(j).lt.th(j-1)-pig)then 
          ng(j)=ng(j-1)+1 
       else 
          ng(j)=ng(j-1) 
       endif
    ENDDO 
    th(1:ntf)=th(1:ntf)+dpig*ng(1:ntf)
!   find frequency by fit                                               
    call linfi3(tf,th,rate,rms,cost,dy,ntf) 
    per=abs(dpig/rate) 
  END SUBROUTINE argum
! =====================================================                 
  SUBROUTINE forced(x,y,tf,ntf,gp,klist,fe,feph,fesp,feco,fed) 
    INTEGER, INTENT(IN) ::  ntf 
    DOUBLE PRECISION, INTENT(IN) :: tf(ntf) 
    DOUBLE PRECISION, INTENT(INoUT) :: x(ntf),y(ntf)
    DOUBLE PRECISION, INTENT(IN) ::  gp(nforc)
    DOUBLE PRECISION, INTENT(out) :: fe(nforc),feph(nforc) &
     &    ,fesp(nforc),feco(nforc),fed(nforc)                           
    INTEGER, INTENT(IN) ::  klist(nforc) 
!  workspaces                                                           
!    DOUBLE PRECISION  vs(ntx),vc(ntx),dy(ntx) 
!    DOUBLE PRECISION yy(ntx),xx(ntx)
    INTEGER itest,k 
    DOUBLE PRECISION d0,d1,d2,d0k,d1k,d2k,per,ampk,sp 
    DOUBLE PRECISION princ 
    itest=1 
    DO k=5,9 
       if(klist(k).eq.0) CYCLE 
       per=360*3600/gp(k) 
! ********** questa e' una vaccata; usare i complessi ***********       
!       call peri2(tf,y,dy,vc,vs,ntf,per,fesp(k),itest,d0,d1,d2) 
!       call peri2(tf,x,dy,vc,vs,ntf,per,sp,itest,d0k,d1k,d2k) 
       call peri2(tf,y,ntf,per,fesp(k),itest,d0,d1,d2) 
       call peri2(tf,x,ntf,per,sp,itest,d0k,d1k,d2k) 
       fe(k)=sqrt(d1*d1+d2*d2) 
       ampk=sqrt(d1k*d1k+d2k*d2k) 
       fed(k)=ampk-fe(k) 
       feph(k)=princ(atan2(d1,d2))*degrad 
       feco(k)=d0 
    ENDDO
  END SUBROUTINE forced
! =====================================                                 
!  small routines                                                       
! =====================================                                 
! RMS with respect to mean                                              
  DOUBLE PRECISION FUNCTION sigma(x,n) 
    INTEGER, INTENT(IN) :: n 
    DOUBLE PRECISION, INTENT(IN) ::  x(n) 
    INTEGER i
    DOUBLE PRECISION av
    av=SUM(x(1:n))/n
    sigma=0.d0 
    DO i=1,n 
       sigma=sigma+(x(i)-av)**2 
    ENDDO
    sigma=sqrt(sigma/n) 
  END FUNCTION sigma
! =====================================                                 
! RMS with respect to fixed value                                       
  DOUBLE PRECISION FUNCTION sigmf(x,av,n) 
    INTEGER, INTENT(IN) :: n
    DOUBLE PRECISION, INTENT(IN) ::  x(n), av
    INTEGER i 
    sigmf=0.d0 
    DO i=1,n 
       sigmf=sigmf+(x(i)-av)**2
    ENDDO
    sigmf=sqrt(sigmf/n) 
  END FUNCTION sigmf
! =====================================                                 
! RMS with respect to fixed value, for angles in deg                                       
  DOUBLE PRECISION FUNCTION sigmfprideg(x,av,n) 
    INTEGER, INTENT(IN) :: n
    DOUBLE PRECISION, INTENT(IN) ::  x(n), av
    DOUBLE PRECISION avrad, xrad, pridif
    INTEGER i 
    sigmfprideg=0.d0 
    avrad=av*radeg
    DO i=1,n 
       xrad=x(i)*radeg
       sigmfprideg=sigmfprideg+pridif(xrad,avrad)**2
    ENDDO
    sigmfprideg=sqrt(sigmfprideg/n)*degrad
  END FUNCTION sigmfprideg



END MODULE synthcomp

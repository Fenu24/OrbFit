MODULE synthtro
  USE fund_const
! computation of synthetic proper elements for trojans
  IMPLICIT NONE
  PRIVATE

! public routines
  PUBLIC parprtt, doinc, prop, argum, forced

! shared data
  INTEGER, PUBLIC :: ndap,nshi,isum,iste
 
! array size for all the data read at once 
  INTEGER, PARAMETER, PUBLIC :: nsx=501 ! max no of bodies
  INTEGER, PARAMETER, PUBLIC :: nastx=1100 ! max no of bodies
  INTEGER, PARAMETER, PUBLIC :: ntx=50001  ! max no of records
  INTEGER, PARAMETER, PUBLIC :: nforc=9   ! max no forced terms
!  max no. running boxes
  INTEGER, PARAMETER, PUBLIC :: nbx=100

! ====================================================
! former parprtt.h parameters for propert
!  INTEGER, PARAMETER, PUBLIC :: nsx=501
!  INTEGER, PARAMETER, PUBLIC :: nforc=9  ! no forced terms
! ====================================================
! former parprtt.h parameters for propert 5Myr integration
!  INTEGER, PARAMETER, PUBLIC :: nastx=1000
!  INTEGER, PARAMETER, PUBLIC :: ntx=5001
!  INTEGER, PARAMETER, PUBLIC :: ndap=990
!  INTEGER, PARAMETER, PUBLIC :: nshi=250
!  INTEGER, PARAMETER, PUBLIC :: isum=24
!  INTEGER, PARAMETER, PUBLIC :: iste=2
! ====================================================
! former parprtt.h parameters for propert 50 Myr extended integration
!  INTEGER, PARAMETER, PUBLIC :: nastx=135
!  INTEGER, PARAMETER, PUBLIC :: ntx=50001
!  INTEGER, PARAMETER, PUBLIC :: ndap=9900
!  INTEGER, PARAMETER, PUBLIC :: nshi=2500
!  INTEGER, PARAMETER, PUBLIC :: isum=240
!  INTEGER, PARAMETER, PUBLIC :: iste=20


CONTAINS
! ==================================================================
! selection of parameters for running box
! ==================================================================
SUBROUTINE parprtt(inflag,nbb,ntt,ndap,nshi,isum,iste)
  INTEGER, INTENT(IN) :: inflag
  INTEGER, INTENT(OUT) :: nbb,ntt,ndap,nshi,isum,iste
  ! END INTERFACE
  ! former parprtt.h parameters for propert
  IF(inflag.eq.1)THEN
! parameters for propert 5Myr integration
     nbb=1000
     ntt=5001
     IF(ntt.gt.ntx)THEN
        WRITE(*,*)' parprtt: too many records, ntx=',ntt,    &
     &         'max was ntxx=',ntx
        STOP 
     ENDIF
     ndap=990
     nshi=250
     isum=24
     iste=2
  ELSEIF(inflag.eq.2)THEN
     !  parameters for propert 50 Myr extended integration
     nbb=150
     ntt=50001
     IF(ntt.gt.ntx)THEN
        WRITE(*,*)' parprtt: too many records, ntx=',ntt,    &
     &         'max was ntxx=',ntx
        STOP 
     ENDIF
     ndap=9900
     nshi=2500
     isum=240
     iste=20
  ELSE
     WRITE(*,*)' parprtt: option not known, inflag=', inflag
     STOP
  ENDIF
END SUBROUTINE parprtt
! =====================================================                 
SUBROUTINE doinc(x,y,tf,ntf,sn,klis,pe,dpe,fre,dfre,iwri) 
  INTEGER ntf 
  DOUBLE PRECISION x(ntf),y(ntf),tf(ntf),pe,dpe,fre,dfre 
! workspace                                                             
  DOUBLE PRECISION omeg(ntx) 
!  planetary theory, and forced terms                                   
  DOUBLE PRECISION fe(nforc),feph(nforc),fesp(nforc) 
  DOUBLE PRECISION feco(nforc),fed(nforc),sn(nforc) 
  INTEGER klis(nforc) 
  INTEGER k,iwri 
  DOUBLE PRECISION rate,peri,rms,cost,ph,sig,sp 
! function                                                              
!  DOUBLE PRECISION vmax,vmin 
! =====================================================                 
!   2: determination of frequencies by fit to arguments                 
!   all the linear forced terms are removed first                       
! =====================================================                 
  CALL forced(x,y,tf,ntf,sn,klis,fe,feph,fesp,feco,fed) 
  IF(iwri.eq.1)THEN 
     DO  k=5,9 
        IF(klis(k).ne.0)THEN 
           write(9,127)k,fe(k),feph(k),fesp(k)*100,feco(k),fed(k) 
127        format(' forced ',i2,' amp ',1p,d12.5,0p,' ph ',f9.4,     &
     &        ' S% ',f8.4,' cost ',1p,d12.4,' unc. ',d12.4)             
        ENDIF
     ENDDO
  ENDIF
! =====================================================                 
!   now compute argument Omega from q,p, varpi from k,h                 
! =====================================================                 
  CALL argum(x,y,tf,ntf,omeg,rate,peri,rms,cost) 
  fre=rate*degrad*3.6d3 
  dfre=rms/(maxval(tf(1:ntf))-minval(tf(1:ntf)))*degrad*3.6d3 
  if(iwri.eq.1)then 
     write(9,128)rate,peri,rms,cost*degrad 
128  format(' Argument:',' freq ',1p,d13.6,0p,' per '           &
     &,f16.6,' rms ',1p,d12.4,0p,' ph ',f9.4)                           
  endif
!======================================================                 
!   3: find the proper amplitude and phase                              
! =====================================================                 
  call prop(omeg,x,y,ntf,dpig,sp,pe,dpe,ph,cost,sig) 
  if(iwri.eq.1)then 
     write(9,129)pe,ph,sp*100,cost,dpe,sig 
129  format(' Proper ',1p,d13.6,0p,' arg ',f9.4,                       &
     &        ' S% ',f8.4,' c',1p,d12.4,' unc ',d9.2,                   &
     &        ' rms',d9.2)                                              
  endif
! =====================================================                 
! check for possible secular resonances                                 
! =====================================================                 
!     If(iwri.eq.0)RETURN                                               
  DO k=5,9 
     IF(klis(k).ne.0)THEN 
        IF(abs(fre-sn(k)).lt.2.d0.or.fe(k).ge.pe)THEN 
           WRITE(9,113)k,sn(k),fre,fe(k),pe 
113        FORMAT('Sec.res.? ',i2,1x,f8.4,1x,f8.4,' ampl ',f8.6,1x,f8.6)
        ENDIF
     ENDIF
  ENDDO
END SUBROUTINE doinc
! =====================================================                 
SUBROUTINE prop(tf,x,y,ntf,per,sp,amp,damp,ph,cost,sig) 
  INTEGER ntf 
  DOUBLE PRECISION per,sp,amp,damp,ph,cost,sig,sigk,ampk 
!  workspaces                                                           
  DOUBLE PRECISION x(ntf),y(ntf),tf(ntf) 
! functions                                                             
  DOUBLE PRECISION sigma,princ 
  DOUBLE PRECISION d0,d1,d2,d0k,d1k,d2k 
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
  ph=princ(atan2(d1,d2))*degrad 
  cost=d0 
END SUBROUTINE prop
! ========================================================              
!   compute argument varpi from k,h and find frequency                  
SUBROUTINE argum(x,y,tf,ntf,th,rate,per,rms,cost) 
  INTEGER ntf 
  DOUBLE PRECISION x(ntf),y(ntf),tf(ntf),rate,per,rms,cost 
!  workspaces                                                           
  DOUBLE PRECISION th(ntx),dy(ntx),xm,ym 
! functions                                                             
  DOUBLE PRECISION princ
  INTEGER ng(ntx),j 
!  use as center the mean                                               
  xm=sum(x(1:ntf))/ntf 
  ym=sum(y(1:ntf))/ntf 
!  initialise number rev. at 0                                          
  ng(1)=0 
  th(1)=atan2(y(1),x(1)) 
  DO 23 j=2,ntf 
     th(j)=princ(atan2(y(j)-ym,x(j)-xm)) 
     if(th(j).gt.th(j-1)+pig)then 
        ng(j)=ng(j-1)-1 
     elseif(th(j).lt.th(j-1)-pig)then 
        ng(j)=ng(j-1)+1 
     else 
        ng(j)=ng(j-1) 
     endif
23 ENDDO
  DO j=1,ntf 
     th(j)=th(j)+dpig*ng(j) 
  ENDDO
!   find frequency by fit                                               
  CALL linfi3(tf,th,rate,rms,cost,dy,ntf) 
  per=dpig/rate 
END SUBROUTINE argum
! =====================================================                 
SUBROUTINE forced(x,y,tf,ntf,gp,klist,fe,feph,fesp,feco,fed) 
  INTEGER ntf 
  DOUBLE PRECISION x(ntf),y(ntf),tf(ntf) 
  DOUBLE PRECISION gp(nforc),fe(nforc),feph(nforc)                  &
     &    ,fesp(nforc),feco(nforc),fed(nforc)                           
  INTEGER klist(nforc) 
!  workspaces                                                           
!  DOUBLE PRECISION  vs(ntx),vc(ntx),dy(ntx) 
  INTEGER itest,k 
  DOUBLE PRECISION d0,d1,d2,d0k,d1k,d2k,per,ampk,sp 
  DOUBLE PRECISION princ 
  itest=1 
  DO 22 k=5,9 
     if(klist(k).eq.0) CYCLE 
     per=360*3600/gp(k) 
! ********** questa e' una vaccata; usare i complessi ***********       
     call peri2(tf,y,ntf,per,fesp(k),itest,d0,d1,d2) 
     call peri2(tf,x,ntf,per,sp,itest,d0k,d1k,d2k) 
     fe(k)=sqrt(d1*d1+d2*d2) 
     ampk=sqrt(d1k*d1k+d2k*d2k) 
     fed(k)=ampk-fe(k) 
     feph(k)=princ(atan2(d1,d2))*degrad 
     feco(k)=d0 
22 ENDDO
END SUBROUTINE forced

END MODULE synthtro

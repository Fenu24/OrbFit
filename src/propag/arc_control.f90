! module contains the routines to assess the quality of an observed arc.

MODULE arc_control
  USE astrometric_observations
  USE attributable
  USE output_control
  IMPLICIT NONE

  PUBLIC arc_type, new_quality, check_nobs

  INTEGER, PARAMETER :: ngapx= 12
 ! array of  observations before gaps: unsorted, sorted
  DOUBLE PRECISION,DIMENSION(ngapx) :: tgapv,tgapv_srt 
  INTEGER, DIMENSION(ngapx) :: indgsrt ! sorting (by time) index 
  INTEGER :: ngap ! no of gaps found  

CONTAINS
   LOGICAL FUNCTION check_nobs(m,quamin,quamax)
     INTEGER,INTENT(IN):: m,quamin,quamax
     check_nobs=.true.
     IF(quamax.ge.5)THEN
        check_nobs=m.ge.3
     ELSEIF(quamax.eq.4)THEN
        check_nobs=m.ge.4
     ELSEIF(quamax.eq.3)THEN
        check_nobs=m.ge.5
     ELSEIF(quamax.eq.2)THEN
        check_nobs=m.ge.6
     ELSEIF(quamax.eq.1)THEN
        check_nobs=m.ge.6
     ELSE
        WRITE(*,*)' check_nobs: not understood ', m,quamin,quamax
        STOP 
     ENDIF
   END FUNCTION check_nobs

   INTEGER FUNCTION new_quality(m,artyp,nigarc,fail_arty,chi,arc)
     INTEGER,INTENT(IN):: m,artyp,nigarc,fail_arty
     DOUBLE PRECISION, INTENT(IN) :: chi,arc
     IF(arc.gt.180.d0) THEN
!     IF(arc.gt.180.d0.and.artyp.gt.5)THEN
        new_quality=1
        IF(m.lt.9)THEN
           WRITE(ierrou,*)'new_quality: case not considered '
           WRITE(ierrou,*)m,artyp,nigarc,fail_arty,chi,arc,new_quality
           numerr=numerr+1
        ENDIF
     ELSEIF(artyp.ge.5)THEN
        new_quality=2
        IF(m.lt.7)THEN
           WRITE(ierrou,*)'new_quality: case not considered '
           WRITE(ierrou,*)m,artyp,nigarc,fail_arty,chi,arc,new_quality
           numerr=numerr+1
        ENDIF
     ELSEIF(artyp.eq.4)THEN
        new_quality=3 
        IF(m.lt.6)THEN
           WRITE(ierrou,*)'new_quality: case not considered '
           WRITE(ierrou,*)m,artyp,nigarc,fail_arty,chi,arc,new_quality
           numerr=numerr+1
        ENDIF
     ELSEIF(artyp.eq.3)THEN
        new_quality=4
        IF(m.lt.4)THEN
           WRITE(ierrou,*)'new_quality: case not considered '
           WRITE(ierrou,*)m,artyp,nigarc,fail_arty,chi,arc,new_quality
           numerr=numerr+1
        ENDIF
     ELSEIF(artyp.eq.2.and.chi.gt.10.d0)THEN
        new_quality=5
        IF(m.lt.3)THEN
           WRITE(*,*)'new_quality: case not considered '
           WRITE(*,*)m,artyp,nigarc,fail_arty,chi,arc,new_quality
           numerr=numerr+1
        ENDIF
     ELSEIF(artyp.eq.2.and.chi.le.10.d0)THEN
        new_quality=6
        IF(m.lt.3)THEN
           WRITE(ierrou,*)'new_quality: case not considered '
           WRITE(ierrou,*)m,artyp,nigarc,fail_arty,chi,arc,new_quality
           numerr=numerr+1
        ENDIF
     ELSEIF(artyp.eq.2)THEN
        WRITE(ierrou,*)'new_quality: case not considered '
        WRITE(ierrou,*)m,artyp,nigarc,fail_arty,chi,arc
        numerr=numerr+1
        new_quality=7
     ELSEIF(artyp.eq.1.and.m.ge.3)THEN
        new_quality=7
     ELSEIF(artyp.gt.204)THEN
        new_quality=1
     ELSEIF(artyp.ge.201)THEN
        new_quality=2
     ELSEIF(m.eq.2)THEN
        new_quality=8
        IF(artyp.ne.1) GOTO 99
     ELSEIF(m.eq.1)THEN
        new_quality=9
     ELSE
        GOTO 99
     ENDIF
     RETURN
99   CONTINUE ! error stop
     WRITE(ierrou,*)'new_quality: case not considered '
     WRITE(ierrou,*)m,artyp,nigarc,fail_arty,chi,arc
     numerr=numerr+1
     STOP 
   END FUNCTION new_quality


! =============================================
!  arc_type: computes arc type and number of nights
! INPUT: 
!     m, obs,obsw : observed arc
! OUTPUT:
!     arc-type : integer arc type (0 if m=1)
!     geoc_chi,acce_chi : normalized vals. of geodetic curv., acceleration
!     chi : significativity of vector geocurv,etadot
!     nig : number of observing nights
!     fail_flag : should indicate problems found in splitting
!                 into TSA; possible causes
!           1) Z sign in arc splittable in 2 TSA
!           2) Z sign in TSA
!           3) ?????
!           the fail flag should indicate unambigously what
!           happened, but....  
! =============================================
  INTEGER FUNCTION arc_type(obs,obsw,m,geoc_chi,acce_chi,chi,nig,fail_flag,eta)
! =========OBSERVATIONS ==========================================
    INTEGER, INTENT(IN) ::  m              ! observations number
    TYPE(ast_obs), INTENT(IN) :: obs(m)    ! supposed sorted by time
    TYPE(ast_wbsr), INTENT(IN)  :: obsw(m) ! wheigts, residuals no
! =========OUTPUT=================================================
    INTEGER, INTENT(OUT) :: nig,fail_flag  ! no .nights, failure level
    DOUBLE PRECISION, INTENT(OUT) :: geoc_chi, acce_chi, chi
    DOUBLE PRECISION, INTENT(OUT), OPTIONAL :: eta ! proper motion
! end interface
    DOUBLE PRECISION eta0
! local variables
    INTEGER nat, fail_rec, natc, fail_check, nights, j
! ================================================================
!    WRITE(*,*)'arc_type: m ' ,m
    IF(m.le.1)THEN
       WRITE(ierrou,*)' arc_type: too few obs, m=',m
       numerr=numerr+1
       nig=1
       arc_type=0
       geoc_chi=0.d0
       acce_chi=0.d0
       chi=0.d0
       RETURN
    ELSEIF(m.eq.2)THEN
       WRITE(ierrou,*)' arc_type: too few obs, m=',m
       numerr=numerr+1
       arc_type=1
       nig=nights(m,obs,obsw)
       geoc_chi=0.d0
       acce_chi=0.d0
       chi=0.d0
       RETURN
    ENDIF    
! check sorting of observations
    DO j=2,m
      IF(obs(j)%time_tdt.lt.obs(j-1)%time_tdt)THEN
         WRITE(*,*)' arc_type: observations not in time order ',j-1,obs(j-1)%time_tdt,j,obs(j)%time_tdt
         STOP
      ENDIF
    ENDDO
    ngap=0
    nat=rec_arc_type(obs,obsw,m, geoc_chi,acce_chi,chi, nig,fail_rec,eta0)
    IF(PRESENT(eta))THEN
       eta=eta0
    ENDIF
!    WRITE(*,*)'arc_type: nat,nig ', nat,nig
!    WRITE(*,*)tgapv(1:ngap)
    IF(ngap.gt.1)THEN
       CALL heapsort(tgapv,ngap,indgsrt)
       tgapv_srt(1:ngap)=tgapv(indgsrt(1:ngap))
       CALL check_cons_arcs(obs,obsw,m,nat,natc,fail_check)
       arc_type=natc
       fail_flag=1000*fail_check+fail_rec
    ELSE
       arc_type=nat
       fail_flag=fail_rec
    ENDIF
  END FUNCTION arc_type

! check existence of curvature (and no Z sign) for each consecutive
! couple of arcs
  SUBROUTINE check_cons_arcs(obs,obsw,m,nat,natc,fail_flag)
 ! =========OBSERVATIONS ==========================================
    INTEGER, INTENT(IN) ::  m              ! observations number
    TYPE(ast_obs), INTENT(IN) :: obs(m)    ! supposed sorted by time
    TYPE(ast_wbsr), INTENT(IN)  :: obsw(m) ! wheigts, residuals no
    INTEGER, INTENT(IN) :: nat             ! arc type from recursion
! =========OUTPUT=================================================
    INTEGER, INTENT(OUT) :: natc,fail_flag  ! corrected arc type, failure level
! ================================================================
    INTEGER j, i, mm, m1, m2, k
    LOGICAL error, sign_curv, bad_fit
    TYPE(attrib) att
    DOUBLE PRECISION eps, geoc_chi, acce_chi, chi
! ================================================================
    natc=nat
    fail_flag=0
    i=0
    eps=100*epsilon(1.d0)
    DO j=1,ngap
      i=i+1 ! gap index 
      IF(i.gt.ngap) EXIT
      IF(i.eq.1)THEN
         m1=1   ! join from the beginning to before gap 2
         DO k=1,m
           IF(abs(obs(k)%time_tdt-tgapv_srt(2)).lt.eps)THEN
              m2=k
              EXIT
           ENDIF
         ENDDO           
      ELSEIF(i.eq.ngap)THEN
         m2=m;  ! join from after gap ngap-1 to the end
         DO k=1,m
           IF(abs(obs(k)%time_tdt-tgapv_srt(i-1)).lt.eps)THEN
              m1=k+1
              EXIT
           ENDIF
         ENDDO    
      ELSE
         DO k=1,m ! join from after gap i-1 to before gap i+1
           IF(abs(obs(k)%time_tdt-tgapv_srt(i-1)).lt.eps)THEN
              m1=k+1
           ELSEIF(abs(obs(k)%time_tdt-tgapv_srt(i+1)).lt.eps)THEN
              m2=k
           ENDIF
         ENDDO
      ENDIF
      mm=m2-m1+1  
      IF(mm.le.1)THEN
         natc=nat
         fail_flag=99
         RETURN
      ENDIF
      CALL attri_comp(mm, obs(m1:m2), obsw(m1:m2), att, error)
      IF(error)THEN
         WRITE(ierrou,*) '***check_cons_arc: error from attri_comp, mm=', mm
         numerr=numerr+1
         CYCLE
      ENDIF
      CALL test_curv(att, geoc_chi, acce_chi, chi, sign_curv, bad_fit)
      IF(sign_curv.and..not.bad_fit)THEN
! all ok, do nothing
      ELSEIF(sign_curv.and.bad_fit)THEN
! error condition: split is correct, but Z sign
         fail_flag=fail_flag+1
      ELSEIF(.not.sign_curv.and.bad_fit)THEN
! error condition: split is dubious, Z sign
         fail_flag=fail_flag+10
      ELSEIF(.not.sign_curv.and..not.bad_fit)THEN
! two consecutive TSA can be joined in a TSA
         natc=natc-1
         i=i+1         
      ENDIF
    ENDDO
    IF(natc.lt.nat)THEN
!       WRITE(*,*)' check_cons_arc: correction to type ', natc,nat
!       WRITE(ierrou,*)' check_cons_arc: correction to type ', natc,nat
!       numerr=numerr+1
    ENDIF
  END SUBROUTINE check_cons_arcs

! recursive split until no curvature/Z sign
  RECURSIVE INTEGER FUNCTION rec_arc_type(obs,obsw,m,geoc_chi,acce_chi, &
&                        chi,nig,fail_flag,eta) RESULT(nat) 
! =========OBSERVATIONS ==========================================
    INTEGER, INTENT(IN) ::  m              ! observations number
    TYPE(ast_obs), INTENT(IN) :: obs(m)    ! supposed sorted by time
    TYPE(ast_wbsr), INTENT(IN)  :: obsw(m) ! wheigts, residuals no
! =========OUTPUT=================================================
    INTEGER, INTENT(OUT) :: nig,fail_flag  ! no nights,failure in recursion
    DOUBLE PRECISION, INTENT(OUT) :: geoc_chi, acce_chi, chi
    DOUBLE PRECISION, INTENT(OUT), OPTIONAL :: eta ! proper motion
! ================================================================
    TYPE(attrib) att
    LOGICAl error, sign_curv, bad_fit
    DOUBLE PRECISION gap, tgap, eta0
    DOUBLE PRECISION :: geoc_chi1, acce_chi1, chi1,geoc_chi2, acce_chi2,chi2
    DOUBLE PRECISION, PARAMETER :: arc_min=2.1d-2 ! 30 min minimum
    DOUBLE PRECISION, PARAMETER :: sph_min=2.9d-4 ! 1 arcmin minimum
    INTEGER mgap,n1,n2,fail_flag1,fail_flag2,nig1,nig2,nights
! ================================================================
    CALL attri_comp(m, obs, obsw, att, error) !WARNING: rejected obs???
    IF(PRESENT(eta))THEN
       eta=att%eta
    ENDIF
    IF(error)THEN
       WRITE(ierrou,*) '***rec_arc_type: error from attri_comp, m=', m
       numerr=numerr+1
       nat=1
       fail_flag=0
       nig=nights(m,obs,obsw)
    ENDIF
    CALL test_curv(att,geoc_chi,acce_chi,chi,sign_curv,bad_fit)
!    WRITE(*,*)' rec_arc_type: m, sign_curv, bad_fit, arc ',m,sign_curv,bad_fit, att%arc
    IF(sign_curv.or.bad_fit)THEN
       IF(att%sph.lt.sph_min.or.att%arc.lt.arc_min)THEN
          nat=1
          nig=1
          fail_flag=0
!          WRITE(*,*)' rec_arc_type: too short arc, m, arc ',m,att%arc
       ELSE
! split necessary       
          CALL split_by_gap(obs(1:m)%time_tdt, m, gap, mgap, tgap)
!          WRITE(*,*)' rec_arc_type: m,mgap,gap,ngap ',m,mgap,gap,ngap
          IF(ngap.ge.ngapx)THEN
!             WRITE(ierrou,*)' rec_arc_type: too many split ngap=',ngap
!             numerr=numerr+1
             nat=1
             fail_flag=3000
          ELSE
             ngap=ngap+1
             tgapv(ngap)=tgap
          ENDIF
          IF(mgap.gt.2)THEN
             n1=rec_arc_type(obs(1:mgap),obsw(1:mgap),mgap,geoc_chi1,acce_chi1,chi1,nig1,fail_flag1,eta0)
          ELSEIF(mgap.gt.0)THEN
             n1=1
             nig1=nights(mgap,obs(1:mgap),obsw(1:mgap))
             fail_flag1=0
          ELSE
             n1=1
             nig1=1
             fail_flag=7
          ENDIF
          IF(m-mgap.gt.2)THEN
             n2=rec_arc_type(obs(mgap+1:m),obsw(mgap+1:m),m-mgap,geoc_chi2,acce_chi2,chi2,nig2,fail_flag2,eta0)
          ELSE
             n2=1
             nig2=nights(m-mgap,obs(mgap+1:m),obsw(mgap+1:m))
             fail_flag2=0
          ENDIF
!          WRITE(*,*)' rec_arc_type: n1,n2, nig1,nig2,m0,mp '
!          WRITE(*,*)n1,n2,nig1,nig2,mgap,m-mgap
! logic to combine fail_flags
          nat=n1+n2
          IF(nat.eq.2)THEN
! at the lowest level of splittable arc
             IF(bad_fit)THEN
! the arc should have no Z sign, thus give error warning
                fail_flag=fail_flag1+fail_flag2+4
             ELSE
! significant curvature, no Z sign, OK if the two subarcs are OK     
                fail_flag=fail_flag1+fail_flag2
             ENDIF
          ELSE
! at higher level, Z sign is OK
             fail_flag=fail_flag1+fail_flag2
          ENDIF
! compute recursively the number of nights
          IF(att%arc.lt.0.666d0)THEN
             nig=1
          ELSE
             nig=nig1+nig2
          ENDIF
       ENDIF
    ELSEIF(.not.sign_curv.and..not.bad_fit)THEN
       nat=1
       fail_flag=0
       nig=nights(m,obs,obsw)
!       WRITE(*,*)' rec_arc_type: nat ', nat
    ENDIF
!    WRITE(*,*)' rec_arc_type: nig, arc ', nig, att%arc
  END FUNCTION rec_arc_type

  SUBROUTINE test_curv(att, geoc_chi, acce_chi, chi, sign_curv, bad_fit)
    TYPE(attrib), INTENT(IN) :: att
    DOUBLE PRECISION, INTENT(OUT) :: geoc_chi, acce_chi, chi
    LOGICAL, INTENT(OUT) ::  sign_curv, bad_fit
! values used in the paper
!   DOUBLE PRECISION, PARAMETER :: ratio=1.d0, rmsmax=3.d0
    DOUBLE PRECISION, PARAMETER :: ratio=3.d0, rmsmax=4.d0
    DOUBLE PRECISION c(2,2),g(2,2),detg, chi2
    IF(att%geocurv.eq.0.d0.and.att%etadot.eq.0.d0)THEN
       sign_curv=.false.
    ELSEIF(att%rms_geocurv.gt.0.d0.and.att%rms_etadot.gt.0.d0)THEN
       g(1,1)=att%rms_geocurv**2
       g(2,2)=att%rms_etadot**2
       g(2,1)=att%rms_geocurv*att%rms_etadot*att%c_curvacc
       detg=g(1,1)*g(2,2)-g(2,1)**2
       c(1,1)=g(2,2)/detg
       c(2,2)=g(1,1)/detg
       c(1,2)=-g(2,1)/detg
       chi2=att%geocurv**2*c(1,1)+att%etadot**2*c(2,2)+2*att%geocurv*att%etadot*c(1,2)
       chi=sqrt(chi2)
       geoc_chi=att%geocurv/att%rms_geocurv
       acce_chi=att%etadot/att%rms_etadot
       IF(chi.le.ratio)THEN
          sign_curv=.false.
       ELSE
          sign_curv=.true.
       ENDIF
    ELSE
       WRITE(ierrou,*)'test_curv: inconsistent curvature ', att%nobs,att%lin_fit, &
&              att%geocurv,att%etadot, att%rms_geocurv, att%rms_etadot
       numerr=numerr+1
       sign_curv=.false.       
    ENDIF
    IF(att%rms_obs.le.rmsmax)THEN
       bad_fit=.false.
    ELSE
       bad_fit=.true.
    ENDIF
  END SUBROUTINE test_curv

  SUBROUTINE split_by_gap(time, m, gap, mgap,tgap) 
    INTEGER, INTENT(IN) ::  m    
    DOUBLE PRECISION, INTENT(IN) :: time(m)
    DOUBLE PRECISION, INTENT(OUT) :: gap, tgap ! width of gap, time before gap
    INTEGER, INTENT(OUT) :: mgap ! index of last before gap
    INTEGER i
    IF(m.le.2) STOP ' ********* split_by_gap too small m *******'
    gap=0.d0
    mgap=0
    DO i=2,m
       IF(time(i)-time(i-1).gt.gap)THEN
          gap=time(i)-time(i-1)
          mgap=i-1
          tgap=time(i-1)  
       ENDIF
    ENDDO
    IF(mgap.eq.0)THEN
       WRITE(ierrou,*) ' ********* split_by_gap no gap*******',time
       numerr=numerr+1
    ENDIF
  END SUBROUTINE split_by_gap

END MODULE arc_control

MODULE attributable
USE fund_const
USE name_rules
USE output_control
IMPLICIT NONE
PRIVATE

! new data type: attributable, including covariance but compressed for memory saving
TYPE attrib_c

DOUBLE PRECISION :: tdtobs, tutobs ! average time of observations, MJD, TDT/UT 

DOUBLE PRECISION :: arc, sph ! arc length in time, in angle

INTEGER :: obscodi ! numeric obscode, 500 (geocentric) if nsta >1

INTEGER :: nsta ! number of different observatories, nsta=2 if more than 1

INTEGER :: nobs, ntime! number of observations, 
 ! of different times, of radar obs.

DOUBLE PRECISION, DIMENSION(4) :: angles ! alpha, delta, alphadot, deltadot
                         ! unit radians and radians/day

DOUBLE PRECISION :: apm ! apparent magnitude, average; =0 if not available

DOUBLE PRECISION, DIMENSION(10) :: gv        ! covariance of attributable, only below diagonal

LOGICAL lin_fit  ! are we using the linear fit? if false, the quadratic
                 !  fit is being used

DOUBLE PRECISION :: eta ! proper motion

END TYPE attrib_c

! new data type: attributable, including covariance and curvature information
TYPE attrib

DOUBLE PRECISION :: tdtobs, tutobs ! average time of observations, MJD, TDT/UT 

DOUBLE PRECISION :: arc, sph ! arc length in time, in angle

CHARACTER*7 :: obscod ! alphanumeric obscode, 500 (geocentric) if nsta >1

INTEGER :: nsta ! number of different observatories, nsta=2 if more than 1

INTEGER :: nobs, ntime, nrad ! number of observations, 
 ! of different times, of radar obs.

DOUBLE PRECISION, DIMENSION(4) :: angles ! alpha, delta, alphadot, deltadot
                         ! unit radians and radians/day

DOUBLE PRECISION :: apm ! apparent magnitude, average; =0 if not available

DOUBLE PRECISION, DIMENSION(4,4) :: g        ! covariance of attributable

LOGICAL lin_fit  ! are we using the linear fit? if false, the quadratic
                 !  fit is being used

DOUBLE PRECISION :: eta, geocurv, etadot ! proper motion, geodetic curvature,
                            ! along track acceleration (all on the celestial
                            ! sphere, units radians and days)

DOUBLE PRECISION :: rms_eta, rms_geocurv, rms_etadot, c_curvacc 
                            ! standard deviations as propagated from g3a, g3d, 
                            ! and correlation of geocurv, etadot


DOUBLE PRECISION rms_a, rms_d, rms_obs ! RMS of residuals, for alpha*cos(delta)
                                       ! for delta, combined

!DOUBLE PRECISION acc_corr, curv_corr ! topocentric corrections to etadot, geocurv

END TYPE attrib

DOUBLE PRECISION, DIMENSION(4), PARAMETER :: zero_4d_vect = &
&    (/ 0.d0, 0.d0, 0.d0, 0.d0 /)
DOUBLE PRECISION, DIMENSION(4,4), PARAMETER :: zero_4x4_matrix = &
& RESHAPE((/ 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, &
& 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0/),(/ 4, 4 /))
TYPE(attrib), PARAMETER :: undefined_attrib = ATTRIB( &
&  1.d99,1.d99,     & ! no time defined
&   0.d0,0.d0,      & ! no arc defined
&  '       ',       & ! no obscode 
& 1, 0, 0, 0,       & ! integers
& zero_4d_vect,     & ! null angles
&   -9.99d0,        & ! null magnitude
& zero_4x4_matrix,  & ! null covariance
&  .false.,         & ! lin_fit
&  0.d0,0.d0,0.d0,  & ! pr.m curvatures
&  0.d0,0.d0,0.d0,0.d0,  & ! rms
&  0.d0,0.d0,0.d0  & ! rms fit
&     )  !

DOUBLE PRECISION, PARAMETER :: sphdistx=2.d0 ! spherical distance betw first and last, deg


! new data type: radar attributable, including covariance and curvature information
TYPE attrad

DOUBLE PRECISION :: tdtobs, tutobs ! average time of observations, MJD, TDT/UT 

DOUBLE PRECISION :: arc, sph ! arc length in time, in angle

CHARACTER*7 :: obscod ! alphanumeric obscode, 500 (geocentric) if nsta >1

INTEGER :: nsta ! number of different observatories, nsta=2 if more than 1

INTEGER :: nobs, ntime, nrad ! number of observations, 
 ! of different times, of radar obs.

DOUBLE PRECISION, DIMENSION(4) :: observ ! alpha, delta, rho, rhodot
                         ! unit radians, AU, and au/day

DOUBLE PRECISION :: apm ! apparent magnitude, average; =0 if not available

DOUBLE PRECISION, DIMENSION(4,4) :: g        ! covariance of radar attributable

LOGICAL lin_fit  ! are we using the linear fit? if false, the quadratic
                 !  fit is being used

!DOUBLE PRECISION :: eta, geocurv, etadot ! proper motion, geodetic curvature,
!                            ! along track acceleration (all on the celestial
!                            ! sphere, units radians and days)

!DOUBLE PRECISION :: rms_eta, rms_geocurv, rms_etadot, c_curvacc 
                            ! standard deviations as propagated from g3a, g3d, 
                            ! and correlation of geocurv, etadot

! ****************************************************************
!HINT: are these available???
DOUBLE PRECISION rms_a, rms_d, rms_obs ! RMS of residuals, for alpha*cos(delta)
                                       ! for delta, combined 
! add rms_rho, rms_rhod ???
! ****************************************************************

!DOUBLE PRECISION acc_corr, curv_corr ! topocentric corrections to etadot, geocurv

END TYPE attrad

TYPE(attrad), PARAMETER :: undefined_attrad = ATTRAD( &
&  1.d99,1.d99,     & ! no time defined
&   0.d0,0.d0,      & ! no arc defined
&  '       ',       & ! no obscode 
& 1, 0, 0, 0,       & ! integers
& zero_4d_vect,     & ! null observ
&   -9.99d0,        & ! null magnitude
& zero_4x4_matrix,  & ! null covariance
&  .false.,         & ! lin_fit
!&  0.d0,0.d0,0.d0,  & ! pr.m curvatures
!&  0.d0,0.d0,0.d0,0.d0,  & ! rms
&  0.d0,0.d0,0.d0  & ! rms fit
&     )  !


! public entities
PUBLIC attrib, attrad, undefined_attrib, undefined_attrad,max_sphdist

PUBLIC attri_comp, wri_attri, spher_dist , rea_attri, att_diff

PUBLIC sphdistx, same_station

CONTAINS

!TYPE(attrib_c) FUNCTION drop_curv(att)
!  TYPE(attrib), INTENT(IN) :: att
!  drop_curv%angles=att%angles
!  drop_curv%tdtobs=att%tdtobs 
!  drop_curv%tutobs=att%tutobs
  
!END FUNCTION drop_curv

LOGICAL FUNCTION same_station(obs,m)
  USE astrometric_observations
  INTEGER,INTENT(IN) :: m ! number of observations 
  TYPE(ast_obs),INTENT(IN),DIMENSION(m) :: obs ! observations 
  INTEGER j
  same_station=.true.
  DO j=2,m
    IF(obs(j)%obscod_i.ne.obs(j-1)%obscod_i) THEN
      same_station=.false.
      EXIT
    ENDIF 
  ENDDO
END FUNCTION same_station

SUBROUTINE rea_attri(iunatt,iunrat,name0,att,trou,eof)
  TYPE(attrib), INTENT(OUT) :: att
  CHARACTER*(name_len), INTENT(OUT) :: name0 
  DOUBLE PRECISION, INTENT(OUT) :: trou ! rounded time
  LOGICAL, INTENT(OUT) :: eof
  INTEGER, INTENT(IN) :: iunatt,iunrat ! units for attributable,
                                       ! for curvature info
  DOUBLE PRECISION sec1,sec2
  INTEGER mjd1,mjd2
  DOUBLE PRECISION atrou,dtrou,sa,sd,sadot,sddot,caad,cddd,eta
  CHARACTER*(name_len) name1
  CHARACTER*256 record
  INTEGER yearm
  DOUBLE PRECISION arc2,ds2,curv,accel,curv_unc,acc_unc,eta_unc
! read att file
  READ(iunatt,100,END=2)att%tdtobs,att%angles,trou,atrou,dtrou, &
&        att%nobs,att%arc,name0,att%obscod,    &
&        sa,sd,sadot,sddot,caad,cddd,att%eta,att%apm
100 FORMAT(f13.6,1x,f10.7,1x,f10.7,1p,1x,d12.5,1x,d12.5,1x,0p,      &
     &      f9.2,1x,f8.5,1x,f8.5,1x,i3,1x,f8.2,1x,a9,1x,a3,             &
     &      1p,6(1x,d10.3),0p,1x,f9.4,1x,f5.2)
! compose covariance matrix
  att%g=0.d0
  att%g(1,1)=sa**2
  att%g(2,2)=sd**2
  att%g(3,3)=sadot**2
  att%g(4,4)=sddot**2
  att%g(1,3)=caad*(sa*sadot)
  att%g(3,1)=att%g(1,3)
  att%g(2,4)=cddd*(sd*sddot)
  att%g(4,2)=att%g(2,4)
! find UT of observation
  mjd1=FLOOR(att%tdtobs)
  sec1=att%tdtobs-mjd1
  CALL cnvtim(mjd1,sec1,'TDT',mjd2,sec2,'UTC')
  att%tutobs=mjd2+sec2/86400.d0
! read .rat file
  READ(iunrat,'(A)', END=3)record
  READ(record,101)name1,att%nobs,att%arc,yearm,att%nrad,att%nsta, &
       & att%sph,eta,eta_unc
101 FORMAT(a9,1x,i4,1x,f10.4,1x,i4,1x,i2,1x,i1,3(1x,f12.7))
! 9+1+4+1+10+1+4+1+2+1+1+3*13 = 35+39 = 74 
! name1 = a9
! nobs = i4
! arc = f10.4 
! yearm = i4
! nrad = i2
! nsta = i1
! sph,eta,eta_unc = f12.7
  READ(record(76:173),*)curv,accel,curv_unc,acc_unc, &
      & att%rms_obs,att%rms_a,att%rms_d,att%c_curvacc
! 74 + 4*13+3*13+1+7 = 173
! curv,accel,curv_unc,acc_unc = d12.5
! rms_obs,rms_a,rms_d = d12.5
! c_curvacc = f7.4
  att%sph=att%sph/degrad
  IF(abs(eta-att%eta).gt.1.d-4.or.name0.ne.name1)THEN
     WRITE(*,*)' rea_attri: inconsistency ',name0,' ',name1,att%eta, eta
  ENDIF
  att%ntime=att%nobs  ! guess 
  att%eta=att%eta/degrad
  arc2=(att%arc/2.d0)**2
  ds2=(att%arc*att%eta/2.d0)**2 
  att%geocurv=curv/(ds2*degrad)
  att%etadot= accel/(arc2*degrad)
  att%rms_geocurv=curv_unc/(ds2*degrad)
  att%rms_etadot=acc_unc/(arc2*degrad)
  att%rms_eta=eta_unc/degrad
  eof=.false.          
  RETURN
2 WRITE(*,*)' rea_attri: end of file .att'
  eof=.true.
  RETURN
3 WRITE(*,*)' rea_attri: end of file .rat'
  eof=.true.
END SUBROUTINE rea_attri

SUBROUTINE wri_attri(iunatt,iunrat,name0,att,trou,nvir)
  TYPE(attrib), INTENT(IN) :: att
  CHARACTER*(name_len), INTENT(IN) :: name0 
  DOUBLE PRECISION, INTENT(IN) :: trou ! rounded time
  INTEGER, INTENT(IN) :: iunatt,iunrat ! units for attributable,
                                       ! for curvature info
  INTEGER, INTENT(IN), OPTIONAL :: nvir
  DOUBLE PRECISION atrou,dtrou, princ,sa,sd,sadot,sddot,caad,cddd
  DOUBLE PRECISION sd2,sa2,sdd2,sad2
  INTEGER iday,month,yearm 
  DOUBLE PRECISION hour,arc2,ds2,curv,accel,curv_unc,acc_unc,eta_unc
! write .att file
  IF(iunatt.ge.0)THEN
     IF(iunatt.eq.0)  WRITE(*,201)
201  FORMAT(' t(MJD)        R.A.          DEC.       radot       decdot      tround     rarou     decrou nobs  ', &
 &  '  arc(d)   name   obscod    s(ra)   s(dec)   s(rad)   s(decd)   c(aad)   c(ddd)   pr.m(deg)   appmag')
     atrou=att%angles(1)+(trou-att%tdtobs)*att%angles(3)
     atrou=princ(atrou)
     dtrou=att%angles(2)+(trou-att%tdtobs)*att%angles(4)
     sa2=att%g(1,1)
     sa=sqrt(sa2)
     sd2=att%g(2,2)
     sd=sqrt(sd2)
     sad2=att%g(3,3)
     sadot=sqrt(sad2)
     sdd2=att%g(4,4)
     sddot=sqrt(sdd2)
     caad=att%g(1,3)/(sa*sadot)
     cddd=att%g(2,4)/(sd*sddot)
     IF(PRESENT(nvir))THEN
        IF(iunatt.eq.0)THEN
           WRITE(*,300)att%tdtobs,att%angles,trou,atrou,dtrou,   &
             &        att%nobs,att%arc,name0,att%obscod,nvir,                &
             &        sa,sd,sadot,sddot,caad,cddd,att%eta*degrad,att%apm
        ELSE
           WRITE(iunatt,300)att%tdtobs,att%angles,trou,atrou,dtrou,   &
             &        att%nobs,att%arc,name0,att%obscod,nvir,                &
             &        sa,sd,sadot,sddot,caad,cddd,att%eta*degrad,att%apm
        END IF
     ELSE
        IF(iunatt.eq.0)THEN
           WRITE(*,100)att%tdtobs,att%angles,trou,atrou,dtrou,   &
                &        att%nobs,att%arc,name0,att%obscod,                     &
                &        sa,sd,sadot,sddot,caad,cddd,att%eta*degrad,att%apm
        ELSE
           WRITE(iunatt,100)att%tdtobs,att%angles,trou,atrou,dtrou,   &
                &        att%nobs,att%arc,name0,att%obscod,                     &
                &        sa,sd,sadot,sddot,caad,cddd,att%eta*degrad,att%apm
        END IF
     ENDIF
  ENDIF
300 FORMAT(f13.6,1x,f10.7,1x,f10.7,1p,1x,d12.5,1x,d12.5,1x,0p, &
         &      f9.2,1x,f8.5,1x,f8.5,1x,i3,1x,f8.2,1x,a9,1x,a3,1x,i5,    &
         &      1p,6(1x,d10.3),0p,1x,f9.4,1x,f5.2)
100 FORMAT(f13.6,1x,f10.7,1x,f10.7,1p,1x,d12.5,1x,d12.5,1x,0p,      &
         &      f9.2,1x,f8.5,1x,f8.5,1x,i3,1x,f8.2,1x,a9,1x,a3,          &
         &      1p,6(1x,d10.3),0p,1x,f9.4,1x,f5.2)
! write .rat file header
  IF(iunrat.lt.0)RETURN
  IF(iunrat.eq.0)WRITE(*,200)
200 FORMAT('  name     nobs   arctime year nr st   arcang        pr.m        pmunc      ', &
      & '  geocurv      accel        gcunc        accunc       RMS          RMS(a)       RMS(d)', &
      & '    cor-g-a   time') 
  CALL mjddat(att%tdtobs,iday,month,yearm,hour) 
! output               
  arc2=(att%arc/2.d0)**2
  ds2=(att%arc*att%eta/2.d0)**2                     
  curv=att%geocurv*ds2*degrad
  accel=att%etadot*arc2*degrad
  curv_unc=att%rms_geocurv*ds2*degrad
  acc_unc=att%rms_etadot*arc2*degrad
  eta_unc=att%rms_eta*degrad
  IF(abs(curv).gt.999.d0.or.abs(accel).gt.999.d0.or.abs(att%rms_obs) &
  &    .gt.999.d0.or.abs(att%rms_a).gt.999.d0.or.abs(att%rms_d).gt.999.d0)THEN
     IF(iunrat.eq.0)THEN
        WRITE(*,101)name0,att%nobs,att%arc,yearm,att%nrad,att%nsta,   &
             & att%sph*degrad,att%eta*degrad,eta_unc,                            &
             & curv,accel,curv_unc,acc_unc,att%rms_obs,att%rms_a,att%rms_d,      &
             & att%c_curvacc,att%tdtobs
     ELSE
        WRITE(iunrat,101)name0,att%nobs,att%arc,yearm,att%nrad,att%nsta,   &
             & att%sph*degrad,att%eta*degrad,eta_unc,                            &
             & curv,accel,curv_unc,acc_unc,att%rms_obs,att%rms_a,att%rms_d,      &
             & att%c_curvacc,att%tdtobs
     END IF
  ELSE
     IF(iunrat.eq.0)THEN
        WRITE(*,102)name0,att%nobs,att%arc,yearm,att%nrad,att%nsta,   &
             & att%sph*degrad,att%eta*degrad,eta_unc,                            &
             & curv,accel,curv_unc,acc_unc,att%rms_obs,att%rms_a,att%rms_d,      &
             & att%c_curvacc,att%tdtobs
     ELSE
        WRITE(iunrat,102)name0,att%nobs,att%arc,yearm,att%nrad,att%nsta,   &
             & att%sph*degrad,att%eta*degrad,eta_unc,                            &
             & curv,accel,curv_unc,acc_unc,att%rms_obs,att%rms_a,att%rms_d,      &
             & att%c_curvacc,att%tdtobs
     END IF
  ENDIF
101 FORMAT(a9,1x,i4,1x,f10.4,1x,i4,1x,i2,1x,i1,3(1x,f12.7),       &
         &    1p,4(1x,d12.5),0p,3(1x,d12.5),1x,f7.4,1x,f13.6)
102 FORMAT(a9,1x,i4,1x,f10.4,1x,i4,1x,i2,1x,i1,3(1x,f12.7),       &
         &    4(1x,f12.7),3(1x,f12.7),1x,f7.4,1x,f13.6)
END SUBROUTINE wri_attri

SUBROUTINE attri_comp(m,obs,obsw,att,error,qobs,qpobs,qppobs)
  USE astrometric_observations
  USE reference_systems
!INPUT:  observations
  INTEGER,INTENT(IN) :: m ! number of observations 
  TYPE(ast_obs),INTENT(IN),DIMENSION(m) :: obs ! observations 
  TYPE(ast_wbsr),INTENT(IN),DIMENSION(m) :: obsw ! observation weights 
! OUTPUT: attributable
  TYPE(attrib), INTENT(OUT) :: att 
  LOGICAL, INTENT(OUT), OPTIONAL :: error ! error flag, to avoid stop
  DOUBLE PRECISION, INTENT(OUT), DIMENSION(3),OPTIONAL :: qobs,qpobs,qppobs 
! interpolated second derivative of geocentric position delta_X of the observer; 
! as in Poincare', Bullettin Astronomique, 23 (1906) pages 177-178
! END INTERFACE
  DOUBLE PRECISION, DIMENSION(3) :: dx0, dxp, dxpp, position,velocity
  INCLUDE 'parobx.h90'
  DOUBLE PRECISION, DIMENSION(nobx) :: t, alpha, delta, rmsa,rmsd, &
    &    alcosd,rmsad,alr
  DOUBLE PRECISION, DIMENSION(2,2) :: g2a, g2d ! covariance of linear fit 
  DOUBLE PRECISION, DIMENSION(3,3) :: g3a, g3d ! covariance of quadratic fit
! RMS of residuals: linear fit, quadr. fit
  DOUBLE PRECISION :: rms_2a, rms_3a, rms_2d, rms_3d 
  DOUBLE PRECISION :: s2d(2),s2a(2), s3d(3), s3a(3), cosdtc, sindtc, princ, cosatc,sinatc
  DOUBLE PRECISION detadv(2), tmp2(2), arc, gked(2,2)
  DOUBLE PRECISION :: sw,swx,ww, prscal, dx(3,nobx), tc
  INTEGER ng, ntime ! no rev for unwrapping of alpha, no of distinct times
  INTEGER j, ising,nmag
  CHARACTER*7 idst1
! ================================
  error=.false.
  alpha(1:m)=obs%coord(1)
! unwrap of alpha                                                  
  alr(1)=alpha(1) 
  ng=0 
  DO j=2,m 
     IF(obs(j)%type.eq.'O'.or.obs(j)%type.eq.'S') THEN 
        IF(alpha(j).lt.alpha(j-1)-pig)THEN 
           ng=ng+1 
        ELSEIF(alpha(j).gt.alpha(j-1)+pig)THEN 
           ng=ng-1 
        ENDIF
        alr(j)=alpha(j)+ng*dpig 
!     ELSE
!        WRITE(*,*)' attri_comp:  radar att must be computed otherwise',obs(j)%type
!        STOP
     ENDIF
  ENDDO
  delta(1:m)=obs%coord(2)
  rmsa(1:m)=abs(obsw%rms_coord(1))
  rmsd(1:m)=abs(obsw%rms_coord(2))
! observation times
  t(1:m)=obs%time_tdt
  arc=MAXVAL(t(1:m))-MINVAL(t(1:m))
! total weight (for both coordinates, using area of ellipse) and central time
  sw=0.d0
  swx=0.d0
  DO j=1,m
     IF(obs(j)%type.eq.'O'.or.obs(j)%type.eq.'S') THEN 
        ww=1/(rmsa(j)*rmsd(j)*cos(delta(j)))
     ELSEIF(obs(j)%type.eq.'R') THEN
        ww=1/rmsa(j)
     ELSEIF(obs(j)%type.eq.'V') THEN
        ww=1/rmsd(j)
     ENDIF
     swx=swx+t(j)*ww
     sw=sw+ww
  ENDDO
! central time is weighed mean
  tc=swx/sw
! shift origin of time to tc
  t(1:m)=t(1:m)-tc 
! count distinct times 
  ntime=1
  DO j=2,m
     IF(t(j)-t(j-1).gt.100*epsilon(1.d0))ntime=ntime+1
  ENDDO
  IF(ntime.eq.1)THEN
    WRITE(ierrou,*)' attri_comp: what do you want with one observation?' 
    WRITE(ierrou,*) obs%time_tdt
    error=.true.
    RETURN
  ENDIF 
  att%nrad=0
  DO j=1,m 
     IF(obs(j)%type.eq.'R'.or.obs(j)%type.eq.'V')att%nrad=att%nrad+1 
  ENDDO
! test for mixed station attributables                                  
  idst1=obs(1)%obscod_s 
  att%nobs=m
  att%nsta=1 
  DO j=2,att%nobs
     IF(obs(j)%obscod_s.ne.idst1)att%nsta=2 
  ENDDO
! use no topocentric correction if there are two (or more) stations
  IF(att%nsta.eq.1)THEN
     att%obscod=idst1
  ELSE
     att%obscod='500'
  ENDIF
! fit to delta  
  CALL quadratic_fit(t,delta,rmsd,m,ntime,g2d,g3d,s2d,s3d,   &
&           rms_2d,rms_3d,ising)
! decision on which fit to use  
  IF(ntime.eq.2)THEN
     IF(ABS(s2d(2)).lt.1.d6) THEN
        cosdtc=cos(s2d(2))
        sindtc=sin(s2d(2))
     ELSE
        WRITE(ierrou,*)' attri_comp: ridcolous value  from linear fit' 
        WRITE(ierrou,*) s2d
        error=.true.
        RETURN
     ENDIF
  ELSE
     IF(ABS(s3d(3)).lt.1.d6) THEN
        cosdtc=cos(s3d(3))
        sindtc=sin(s3d(3))
     ELSE
        WRITE(ierrou,*)' attri_comp: ridcolous value from deg 2 fit ' 
        WRITE(ierrou,*) s3d
        error=.true.
        RETURN
     ENDIF
  ENDIF
!  alcosd(1:m)=alr(1:m)*cosdtc NO!!! no approx on tangent space
!  rmsad(1:m)=rmsa(1:m)*cosdtc NO!!!
! fit to alpha
  CALL quadratic_fit(t,alr,rmsa,m,ntime,g2a,g3a,s2a,s3a,   &
&           rms_2a,rms_3a,ising)
! output attributable
  att%tdtobs=tc ! central time in TDT
  att%arc=arc
  att%nobs=m
  att%ntime=ntime
  att%tutobs=SUM(obs%time_utc)/m ! central time in TUT (WRONG!!!)
  att%sph=spher_dist(obs(1)%coord,obs(m)%coord)
  IF(att%sph.le.sphdistx*radeg)      &
&         att%sph= max_sphdist(obs(1:m)%coord(1),obs(1:m)%coord(2),m)
  IF(att%sph.le.1000*epsilon(1.d0))THEN
    WRITE(ierrou,*)' attri_comp: what do you want with a fixed star?' 
    WRITE(ierrou,*) att%sph
    error=.true.
    RETURN
  ENDIF 
! compute average apparent magnitude
  nmag=0
  att%apm=0.d0
  DO j=1,m
     IF(obs(j)%mag_def)THEN
        att%apm=att%apm+obs(j)%mag
        nmag=nmag+1
     ENDIF
  ENDDO
  IF(nmag.ne.0)att%apm=att%apm/nmag
! angles and covariance
  IF(ntime.eq.2)THEN
     att%lin_fit=.true.
     att%angles(1)=princ(s2a(2))
     att%angles(2)=s2d(2)
     att%angles(3)=s2a(1)
     att%angles(4)=s2d(1)
     att%g=0.d0
     att%g(1,1)=g2a(2,2)
     att%g(3,3)=g2a(1,1)
     att%g(1,3)=g2a(2,1)
     att%g(3,1)=g2a(1,2)
     att%g(2,2)=g2d(2,2)
     att%g(4,4)=g2d(1,1)
     att%g(2,4)=g2d(2,1)
     att%g(4,2)=g2d(1,2)
     att%rms_a=rms_2a ! this includes cos(delta)
     att%rms_d=rms_2d
     att%rms_obs=sqrt(rms_2a**2+rms_2d**2)
  ELSE
     att%lin_fit=.false.
     att%angles(1)=princ(s3a(3))
     att%angles(2)=s3d(3)
     att%angles(3)=s3a(2)
     att%angles(4)=s3d(2)
     att%g=0.d0
     att%g(1,1)=g3a(3,3)
     att%g(3,3)=g3a(2,2)
     att%g(1,3)=g3a(3,2)
     att%g(3,1)=g3a(2,3)
     att%g(2,2)=g3d(3,3)
     att%g(4,4)=g3d(2,2)
     att%g(2,4)=g3d(3,2)
     att%g(4,2)=g3d(2,3)
     att%rms_a=rms_3a ! this includes cos(delta)
     att%rms_d=rms_3d
     att%rms_obs=sqrt(rms_3a**2+rms_3d**2)
  ENDIF
! proper motion and its variance
  att%eta=sqrt(att%angles(4)**2+att%angles(3)**2*cosdtc**2)
  detadv(1)=att%angles(3)*cosdtc**2/att%eta
  detadv(2)=att%angles(4)/att%eta
! neglecting the contribution from uncert. of delta
  tmp2=MATMUL(att%g(3:4,3:4),detadv)
  att%rms_eta=sqrt(DOT_PRODUCT(detadv,tmp2))
  IF(ntime.eq.2)THEN
     att%geocurv=0.d0
     att%etadot=0.d0
     att%rms_geocurv=0.d0 !wrong, but please use lin_fit
     att%rms_etadot=0.d0  ! idem
     att%c_curvacc=0.d0
  ELSE
     att%geocurv=(s3d(1)*att%angles(3)-s3a(1)*att%angles(4))*cosdtc/att%eta**3
     att%geocurv=att%geocurv + att%angles(3)*sindtc*(att%eta**2+   &
  &     att%angles(4)**2)/att%eta**3
     att%etadot=(att%angles(3)*s3a(1)*cosdtc**2+ att%angles(4)*s3d(1))/att%eta
     att%etadot=att%etadot -(att%angles(3)**2*sindtc*cosdtc*att%angles(4)) &
  &     /att%eta
! covariance of geodetic curvature, acceleration
     CALL covar_curvacc(att%eta,att%angles(1),att%angles(2),att%angles(3), &
 &       att%angles(4),s3a(1),s3d(1),g3a,g3d,gked)
     att%rms_geocurv=sqrt(gked(1,1)) 
     att%rms_etadot=sqrt(gked(2,2))  
     att%c_curvacc=gked(1,2)/(att%rms_geocurv*att%rms_etadot)
  ENDIF
  IF(PRESENT(qobs).and.PRESENT(qpobs).and.PRESENT(qppobs))THEN

!     IF(ntime.eq.2)THEN ! fara il fit lineare!!!
!        qobs=0.d0
!        qpobs=0.d0
!        qppobs=0.d0
!     ELSE
! geocentric position of the observer, equator and equinox J2000
     DO j=1,m 
        IF(rhs.ne.1.and.rhs.ne.2)THEN
           WRITE(*,*)'attri_comp: rhs=', rhs
           STOP
        END IF
        
        IF(obs(j)%type.eq.'S')THEN
           dx(1:3,j)=obs(j)%obspos(1:3)
        ELSE
           CALL observer_position(obs(j)%time_tdt,position,velocity,obs(j)%obscod_i)
           IF(rhs.eq.1)THEN
              CALL prodmv(dx(1:3,j),roteceq,position)                
           ELSEIF(rhs.eq.2)THEN
              dx(1:3,j)=position                
           END IF
        ENDIF
     ENDDO
! Poincare' interpolation on geocentric position of observer
     DO j=1,3
        CALL quadratic_fit(t,dx(j,1:m),rmsa,m,ntime,g2a,g3a,s2a,s3a,rms_2a,rms_3a,ising)
        
        IF(ntime.eq.2)THEN
           qppobs(j)=0.d0
           qpobs(j)=s2a(1)
           qobs(j)=s2a(2)
        ELSEIF(ntime.gt.2)THEN
           ! WARNING: what to do with residuals, uncertainty of these fits?
           !           qppobs(j)=s3a(1)
           qppobs(j)=2.d0*s3a(1) ! modified on July 7, 2008
           qpobs(j)=s3a(2)
           qobs(j)=s3a(3)
        ELSE
           WRITE(*,*)'attri_comp: error! ntime <2!',ntime
           STOP
        ENDIF
     ENDDO

!  ENDIF

  ENDIF
END SUBROUTINE attri_comp

! covariance of geodetic curvature and acceleration
  SUBROUTINE covar_curvacc(eta,a,d,da,dd,dda,ddd,g3a,g3d,gked)
    DOUBLE PRECISION, INTENT(IN) :: eta,a,d,da,dd,dda,ddd
    DOUBLE PRECISION, DIMENSION(3,3), INTENT(IN) :: g3a,g3d
    DOUBLE PRECISION, INTENT(OUT) :: gked(2,2)
    DOUBLE PRECISION tmp23(2,3),dked(2,3)
    DOUBLE PRECISION t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15
! partials of k w.r. delta and derivatives
!    dked(1,1:3)=\partial k /\partial (\"\delta, \ddot, \delta)
    dked(1,3) = 3.D0*(ddd*da-dda*dd)/sqrt(da**2*cos(d)**2.D0+dd**2)**5.D0* &
& cos(d)**2.D0*da**2*sin(d)-(ddd*da-dda*dd)/sqrt(da**2*cos(d)**2.D0+ &
& dd**2)**3.D0*sin(d)+da*cos(d)/sqrt(da**2*cos(d)**2.D0+dd**2)**3.D0 &
& *(da**2*cos(d)**2.D0+2.D0*dd**2)+3.D0*da**3*sin(d)**2.D0/ &
& sqrt(da**2*cos(d)**2.D0+dd**2)**5.D0*(da**2*cos(d)**2.D0+ &
& 2.D0*dd**2)*cos(d)-2.D0*da**3*sin(d)**2.D0/sqrt(da**2*cos(d)**2.D0+ &
& dd**2)**3.D0*cos(d)
    dked(1,2) = -dda/sqrt(da**2*cos(d)**2.D0+dd**2)**3.D0*cos(d)- &
& 3.D0*(ddd*da-dda*dd)/sqrt(da**2*cos(d)**2.D0+dd**2)**5.D0* &
& cos(d)*dd-3.D0*da*sin(d)/sqrt(da**2*cos(d)**2.D0+ &
& dd**2)**5.D0*(da**2*cos(d)**2.D0+2.D0*dd**2)*dd+ &
& 4.D0*da*sin(d)/sqrt(da**2*cos(d)**2.D0+dd**2)**3.D0*dd
    dked(1,1) =  da/sqrt(da**2*cos(d)**2.D0+dd**2)**3.D0*cos(d)
! partials of etadot w.r. to delta and derivatives
    dked(2,3) =  1.D0/sqrt(da**2*cos(d)**2.D0+dd**2)**3.D0* &
& (da*dda*cos(d)**2.D0+dd*ddd-da**2*cos(d)*sin(d)*dd)* &
& da**2*cos(d)*sin(d)+1.D0/sqrt(da**2*cos(d)**2.D0+ &
& dd**2)*(-2.D0*da*dda*cos(d)*sin(d)+da**2*sin(d)**2.D0*dd- &
& da**2*cos(d)**2.D0*dd)
    dked(2,2) = -1.D0/sqrt(da**2*cos(d)**2.D0+dd**2)**3.D0* &
& (da*dda*cos(d)**2.D0+dd*ddd-da**2*cos(d)*sin(d)*dd)*dd+ &
& 1.D0/sqrt(da**2*cos(d)**2.D0+dd**2)*(ddd-da**2*cos(d)*sin(d))
    dked(2,1) = 1.D0/sqrt(da**2*cos(d)**2.D0+dd**2)*dd
    tmp23=MATMUL(dked,g3d)
    gked=MATMUL(tmp23,TRANSPOSE(dked))
! partials of k w.r. to alpha and derivatives
    dked(1,3) = 0.d0
    dked(1,2) = ddd/sqrt(da**2*cos(d)**2.D0+dd**2)**3.D0*cos(d)- &
& 3.D0*(ddd*da-dda*dd)/sqrt(da**2*cos(d)**2.D0+dd**2)**5.D0 &
& *cos(d)**3.D0*da+sin(d)/sqrt(da**2*cos(d)**2.D0+ &
& dd**2)**3.D0*(da**2*cos(d)**2.D0+2.D0*dd**2)- &
& 3.D0*da**2*sin(d)/sqrt(da**2*cos(d)**2.D0+dd**2)**5.D0* &
& (da**2*cos(d)**2.D0+2.D0*dd**2)*cos(d)**2.D0+ &
& 2.D0*da**2*sin(d)/sqrt(da**2*cos(d)**2.D0+dd**2)**3.D0*cos(d)**2.D0
    dked(1,1) = -dd/sqrt(da**2*cos(d)**2.D0+dd**2)**3.D0*cos(d)
! partials of etadot w.r. to alpha and derivatives
    dked(2,3) =  0.d0
    dked(2,2) = -1.D0/sqrt(da**2*cos(d)**2.D0+dd**2)**3.D0* &
& (da*dda*cos(d)**2.D0+dd*ddd-da**2*cos(d)*sin(d)*dd)* &
& da*cos(d)**2.D0+1.D0/sqrt(da**2*cos(d)**2.D0+ &
& dd**2)*(dda*cos(d)**2.D0-2.D0*da*cos(d)*sin(d)*dd)
    dked(2,1) = 1.D0/sqrt(da**2*cos(d)**2.D0+dd**2)*da*cos(d)**2.D0
    tmp23=MATMUL(dked,g3a)
    gked=gked+MATMUL(tmp23,TRANSPOSE(dked))
  END SUBROUTINE covar_curvacc

! max spherical distance for a set of observations
  DOUBLE PRECISION FUNCTION max_sphdist(alpha,delta,m)
    INTEGER, INTENT(IN) :: m ! nummber of obs
! array of alpha, delta
    DOUBLE PRECISION, DIMENSION(m), INTENT(IN) :: alpha,delta 
    INTEGER i,j
    DOUBLE PRECISION mx, oi(2), oj(2)
    mx=0.d0
    DO i=1,m-1
       oi(1)=alpha(i)
       oi(2)=delta(i)
       DO j=i+1,m
          oj(1)=alpha(j)
          oj(2)=delta(j)
          mx=MAX(mx,spher_dist(oi,oj))
       ENDDO
    ENDDO
    max_sphdist=mx
  END FUNCTION max_sphdist
! spherical distance between two observations
  DOUBLE PRECISION FUNCTION spher_dist(obs1,obs2)
    DOUBLE PRECISION, INTENT(IN), DIMENSION(2) :: obs1,obs2
    DOUBLE PRECISION :: alpha1,delta1,alpha2,delta2
    DOUBLE PRECISION x1(3), x2(3),prscal,cp,sp, xx(3),vsize
! unit vectors
    alpha1=obs1(1)
    delta1=obs1(2)
    x1(1)=cos(delta1)*cos(alpha1)
    x1(2)=cos(delta1)*sin(alpha1) 
    x1(3)=sin(delta1)
    alpha2=obs2(1)
    delta2=obs2(2)
    x2(1)=cos(delta2)*cos(alpha2)
    x2(2)=cos(delta2)*sin(alpha2) 
    x2(3)=sin(delta2)
! cosine as scalar product
    cp=prscal(x1,x2)
!    IF(abs(pp).gt.1.d0-100*epsilon(1.d0))THEN
!       spher_dist=0.d0
!    ELSE
!       spher_dist=acos(pp)
!    ENDIF
    CALL prvec(x1,x2,xx)
    sp=vsize(xx)
    spher_dist=atan2(sp,cp)
END FUNCTION spher_dist

! attributables differences
SUBROUTINE att_diff(att1,att2,datt)
  TYPE(attrib), INTENT(IN) :: att1, att2
  DOUBLE PRECISION, DIMENSION(4), INTENT(OUT) :: datt
  DOUBLE PRECISION pridif
  datt(2:4)=att1%angles(2:4)-att2%angles(2:4)
  datt(1)=pridif(att1%angles(1),att2%angles(1))
END SUBROUTINE att_diff

END MODULE attributable

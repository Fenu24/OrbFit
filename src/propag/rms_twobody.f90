! Copyright (C) GPL 2006 OrbFit Consortium
! version 3.3.2 A. Milani November 2006
! ===============================================================
!  RMS_TWOBODY
!                                                                       
! INPUT:    EL      -  Elements
!           OBS, OBSW -  Observations etc.               
!           IO1,IO2   -  First and last observations to be used         
!                                                                       
! OUTPUT:   in OBSW:   Residuals in RA and DEC   
!           RMS       -  RMS error of residuals (rad)                   
! =============================================================
SUBROUTINE rms_twobody(el,obs,obsw,io1,io2,minrms_prel,rms)
  USE astrometric_observations
  USE ever_pitkin, ONLY: fser_propag
  USE fund_const
  USE orbit_elements
  IMPLICIT NONE 
! INPUT observations: new data types
  INTEGER, INTENT(IN) ::  io1,io2 ! span to be used
  DOUBLE PRECISION, INTENT(IN) :: minrms_prel !minimum acceptable weight
  TYPE(ast_obs),DIMENSION(io2),INTENT(IN) :: obs
  TYPE(ast_wbsr),DIMENSION(io2),INTENT(INOUT) :: obsw
  TYPE(orbit_elem), INTENT(IN) ::  el !elements
  DOUBLE PRECISION, INTENT(OUT) :: rms 
! END INTERFACE
  INTEGER k,ns, fail_flag 
  DOUBLE PRECISION xobs(3),vobs(3),xast(3),vast(3),xeq(3),veq(3)
  DOUBLE PRECISION xrel(3),alphac,deltac,rrel,vsize,xea(6), wei(2) 
  DOUBLE PRECISION, EXTERNAL ::  pridif 
  TYPE(orbit_elem) elcar
! ========================================================
  rms=0.d0 
  ns=0 
! change coordinates
  CALL coo_cha(el,'CAR',elcar,fail_flag)
  IF(fail_flag.ge.4)THEN
     WRITE(*,*)'rms_twobody: failure in coordinate change '
     WRITE(*,*) fail_flag, el,elcar
     STOP
  ENDIF
  DO  k=io1,io2 
     IF(obs(k)%type.ne.'O') CYCLE
     IF(obsw(k)%sel_coord.LT.1) CYCLE
! cartesian ecliptic coordinates of the asteroid at obs.time
     CALL fser_propag(elcar%coord(1:3), elcar%coord(4:6), elcar%t, obs(k)%time_tdt,gms,xeq,veq)
     xast=MATMUL(roteceq,xeq)  ! convert to equatorial
     vast=MATMUL(roteceq,veq) 
! geocentric ecliptic coordinates of the observer
     xobs=obs(k)%obspos
     vobs=obs(k)%obsvel
! JPL Earth vector, ecliptic, at observation time 
     CALL earcar(obs(k)%time_tdt,xea,1)
     xobs=xobs+xea(1:3)
     xobs=MATMUL(roteceq,xobs) !convert to equatorial
! compute relative position and velocity 
     xrel=xast-xobs
     CALL aber1(xrel,vast,xrel)
 ! polar coordinates
     IF(vsize(xrel).le.0.d0)THEN
        obsw(k)%res_coord=0.d0
     ELSE
        CALL polar(xrel,alphac,deltac,rrel) 
        obsw(k)%res_coord(1)=pridif(obs(k)%coord(1),alphac) 
        obsw(k)%res_coord(2)=obs(k)%coord(2)-deltac 
        wei(1)=MAX(obsw(k)%rms_coord(1),minrms_prel)
        wei(2)=MAX(obsw(k)%rms_coord(2),minrms_prel)
        rms=rms+(COS(obs(k)%coord(2))*obsw(k)%res_coord(1)/wei(1))**2+   &
    &       (obsw(k)%res_coord(2)/wei(2))**2
     ENDIF
     ns=ns+1
  ENDDO
  IF(ns.GT.0) rms=SQRT(rms/(2*ns)) 

END SUBROUTINE rms_twobody

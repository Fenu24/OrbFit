! Copyright (C) 2000 by Mario Carpino (carpino@brera.mi.astro.it)       
! Version: June 5, 2000     
! version 3.0 A. Milani November 2002                                            
! --------------------------------------------------------------------- 
!                                                                       
!  *****************************************************************    
!  *                                                               *    
!  *                         I O D _ R M S                         *    
!  *                                                               *    
!  *      Orbital residuals and RMS for preliminary orbits         *    
!  *                                                               *    
!  *****************************************************************    
!                                                                       
! INPUT:    ELEM      -  Orbital elements (mean ecliptic                
!                                          and equinox J2000)           
!           ELTYPE    -  Type of orbital elements                       
!           TELEM     -  Epoch of orbital elements (MJD, TDT) 
!           OBS, OBSW observations               
!           IO1,IO2   -  First and last observations to be used         
!                                                                       
! OUTPUT:   in OBSW:   Residuals in RA and DEC   
!           RMS       -  RMS error of residuals (rad)                   
!                                                                       
      SUBROUTINE iod_rms(elem,eltype,telem,obs,obsw,io1,io2,rms)
      USE reference_systems 
      USE astrometric_observations
      USE fund_const
      IMPLICIT NONE 
! INPUT observations: new data types
      INTEGER, INTENT(IN) ::  io1,io2 ! span to be used
      TYPE(ast_obs),DIMENSION(io2),INTENT(IN) :: obs
      TYPE(ast_wbsr),DIMENSION(io2),INTENT(INOUT) :: obsw
!      INTEGER iobs(io2),obscod(io2),sel(io2) 
      DOUBLE PRECISION, INTENT(IN) ::  elem(6),telem !elements, epoch
!,tobs(io2),alpha(io2),delta(io2) resa(io2),resd(io2),
      DOUBLE PRECISION, INTENT(OUT) :: rms 
      CHARACTER*(*), INTENT(IN) :: eltype 
! END INTERFACE
      INTEGER i,k,ns 
      DOUBLE PRECISION xobs(3),vobs(3),xast(3),vast(3),xv(6) !rot(3,3),gm,
      DOUBLE PRECISION xrel(3),xrel_ecl(3),alphac,deltac,rrel,vsize,xea(6) 
      DOUBLE PRECISION princ 
      EXTERNAL princ 
! ========================================================
      rms=0.d0 
      ns=0 
      DO  k=io1,io2 
         IF(obs(k)%type.ne.'O') CYCLE
         IF(obsw(k)%sel_coord.LT.1) CYCLE
! cartesian ecliptic coordinates of the asteroid at obs.time
         CALL ekcc1(elem,eltype,xv,gms,obs(k)%time_tdt-telem)
         xast=MATMUL(roteceq,xv(1:3))  ! convert to equatorial
         vast=MATMUL(roteceq,xv(4:6)) 
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
            obsw(k)%res_coord(1)=princ(obs(k)%coord(1)-alphac) 
            IF( obsw(k)%res_coord(1).GT.pig) obsw(k)%res_coord(1)=    &
     &       obsw(k)%res_coord(1)-dpig 
            obsw(k)%res_coord(2)=obs(k)%coord(2)-deltac 
            rms=rms+(COS(obs(k)%coord(2))*obsw(k)%res_coord(1)/obsw(k)%rms_coord(1))**2+   &
     &       (obsw(k)%res_coord(2)/obsw(k)%rms_coord(2))**2
         ENDIF
         ns=ns+1
       ENDDO 
                                                                        
      IF(ns.GT.0) rms=SQRT(rms/(2*ns)) 
      RETURN                                                            
      END SUBROUTINE iod_rms                        

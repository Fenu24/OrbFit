! Copyright (C) 2007 by  Milani, Gronchi, Farnocchia 
! ====================================================
! LAPLACE_POINCARE
! ====================================================
SUBROUTINE laplace_poincare(m,obs,obsw,el,nroots,nsol,rr,fail,msg,debug)
  USE reference_systems 
  USE fund_const
  USE orbit_elements
  USE output_control
  USE astrometric_observations
  USE attributable
  USE station_coordinates, ONLY: statcode
  IMPLICIT NONE  
!INPUT:  observations
  INTEGER,INTENT(IN) :: m ! number of observations 
  TYPE(ast_obs),INTENT(IN),DIMENSION(m) :: obs ! observations 
  TYPE(ast_wbsr),INTENT(IN),DIMENSION(m) :: obsw ! observation weights 
  LOGICAL, INTENT(IN) :: debug
! OUTPUT
  TYPE(orbit_elem), INTENT(OUT) :: el(3)
  DOUBLE PRECISION, INTENT(OUT) :: rr(3) ! topocentric distance
  LOGICAL, INTENT(OUT) :: fail
  CHARACTER*(*), INTENT(OUT) :: msg
  INTEGER, INTENT(OUT) :: nsol,nroots
! end interface
! to call attri_comp
  TYPE(attrib) :: att ! attributable
  LOGICAL error
  DOUBLE PRECISION, DIMENSION(3) ::  qobs, qpobs, qppobs ! observer, first and second deriv
  DOUBLE PRECISION, DIMENSION(3) ::  rhat, vhat,nhat,qhat,xp,vp
  DOUBLE PRECISION, DIMENSION(6) :: qea, xv
  DOUBLE PRECISION cosd,sind,cosa,sina,qe,qhdnh,qhdvh,lambda_n,lambda_v,cosep, tc,vsize,prscal
  DOUBLE PRECISION C0,h0,coef(0:8),roots(8),rho,rhodot
  INTEGER obsco ! station code, numeric
  INTEGER ir ! loop on solutions of deg. 8 equation
  INTEGER fail_flag ! for coordinate change
  TYPE(orbit_elem) elk 
! ===============================================================
  fail=.true. 
  msg=' ' 
  nroots=0 
  nsol=0 
! attributable and coefficients lambda, lambda_v etc.
  CALL attri_comp(m,obs,obsw,att,error,qobs,qpobs,qppobs)
! orthonormal reference system
  cosd=cos(att%angles(2))
  sind=sin(att%angles(2))
  sina=sin(att%angles(1))
  cosa=cos(att%angles(1))
  rhat(1)=cosd*cosa
  rhat(2)=cosd*sina
  rhat(3)=sind  
  vhat(1)=(-att%angles(3)*cosd*sina-att%angles(4)*sind*cosa)/att%eta
  vhat(2)=(att%angles(3)*cosd*cosa-att%angles(4)*sind*sina)/att%eta
  vhat(3)=att%angles(4)*cosd/att%eta
  CALL prvec(rhat,vhat,nhat)
  WRITE(*,*) prscal(rhat,vhat), prscal(rhat,nhat), vsize(vhat) 
! Earth position at central time
  tc=att%tdtobs
  CALL earcar(tc,qea,1)
! distance to Sun, rotate to equatorial
  qe=vsize(qea(1:3))
  qea(1:3)=MATMUL(roteceq,qea(1:3))
  qea(4:6)=MATMUL(roteceq,qea(4:6))
  qhat=qea(1:3)*(1.d0/qe)
! cosine terms
  qhdnh=prscal(qhat,nhat)
  qhdvh=prscal(qhat,vhat)
  cosep=prscal(qhat,rhat)
! output lambda, lambda_v
  lambda_n= prscal(qppobs,nhat)*qe**2/(gms*qhdnh)
  lambda_v=prscal(qppobs,vhat)*qe**2/(gms*qhdvh)
! coeff. of dynamical equation
  h0=1.d0-lambda_n     
  C0=att%eta**2*att%geocurv*qe**3/(gms*qhdnh)
! Coefficients of 8th degree equation for r2                            
  coef(8)=C0**2
  coef(7)=0.d0 
  coef(6)=-qe**2*(h0**2+2*C0*h0*cosep+C0**2)
  coef(5)=0.d0 
  coef(4)=0.d0 
  coef(3)=qe**5*2*(h0+C0*cosep)
  coef(2)=0.d0 
  coef(1)=0.d0 
  coef(0)=-qe**8 
  IF(coef(8).EQ.0) THEN 
     msg='coef(8)=0' 
     RETURN 
  END IF
  CALL solv8(coef,roots,nroots) 
  WRITE(iun_log,111)C0,h0,qe,cosep,lambda_v
111 FORMAT('C0=',1P,D11.4,' h0=',D11.4,' QE=',0P,D12.5,' cos(eps)=',D12.5,' lambda_v=',1P,D11.4)
  IF(nroots.LE.0) THEN 
     WRITE(iun_log,524)
524  FORMAT(16X,'no roots of deg.8 equation') 
  ELSE 
     WRITE(iun_log,514)(ir,roots(ir),ir=1,nroots)
514  FORMAT(16X,'r(',i1,')  =',F10.6)
  END IF
  IF(debug) THEN 
     WRITE(*,111)C0,h0,qe,cosep,lambda_v
     IF(nroots.LE.0) THEN 
        WRITE(*,524)
     ELSE 
        WRITE(*,514)(ir,roots(ir),ir=1,nroots)
     END IF
  END IF 
  IF(nroots.LE.0) THEN 
     msg='8th degree polynomial has no real roots' 
     RETURN 
  ELSEIF(nroots.gt.3)THEN
     WRITE(*,*)' Laplace-Poincare: more than 3 positive roots ',nroots
     STOP
  END IF
                                                                        
! Orbital elements of preliminary solution
  nsol=0
  DO 20 ir=1,nroots
     rho=qe/C0*(h0-qe**3/roots(ir)**3)
     xp=qea(1:3)+qobs+rho*rhat
! missing equation for rhodot using att%etadot
     rhodot=(-rho*att%etadot+gms*qhdvh/(qe**2)*(1.d0-lambda_v-qe**3/roots(ir)**3))/(2*att%eta)
     vp=qea(4:6)+qpobs+rho*att%eta*vhat + rhodot*rhat               
! remove spurious root
     IF(rho.lt.0.001d0)THEN
        IF(debug)WRITE(*,*) ' spurious root , r, rho ',roots(ir),rho
        WRITE(iun_log,*) ' spurious root , r, rho ',roots(ir),rho
        CYCLE
     ELSE
        IF(debug)WRITE(*,*) ' accepted root , r, rho ',roots(ir),rho
        WRITE(iun_log,*) ' accepted root , r, rho ',roots(ir),rho
        nsol=nsol+1
        rr(nsol)=rho
     ENDIF
! conversion to ecliptic coordinates
     xv(1:3)=MATMUL(roteqec,xp)
     xv(4:6)=MATMUL(roteqec,vp)
! Orbital elements of preliminary orbit                                 
     el(nsol)=undefined_orbit_elem
     IF(rhs.EQ.2)THEN
        el(nsol)%center=3
     END IF
     el(nsol)%coord=xv
     el(nsol)%coo='CAR'
! planetary aberration 
     el(nsol)%t=att%tdtobs-rho/vlight
     CALL statcode(att%obscod,obsco)
     CALL coo_cha(el(nsol),'ATT',el(nsol),fail_flag,OBSCODE=obsco)
     WRITE(iun_log,525) ir, fail_flag 
     IF(el(nsol)%coo.EQ.'KEP') THEN 
        WRITE(iun_log,535) 'a,e ',el(nsol)%coord(1),el(nsol)%coord(2) 
     ELSEIF(el(nsol)%coo.EQ.'COM'.or.el(nsol)%coo.eq.'COT') THEN 
        WRITE(iun_log,535) 'q,e ',el(nsol)%coord(1),el(nsol)%coord(2)
     ELSEIF (el(nsol)%coo.EQ.'ATT')THEN
        WRITE(iun_log,535)'rho,rhodot',el(nsol)%coord(5),el(nsol)%coord(6)
     ELSE
        STOP '**** laplace_poincare: unknown elem type ****' 
     END IF
     IF(debug) THEN 
        WRITE(*,525) ir, fail_flag 
        IF(el(nsol)%coo.EQ.'KEP') THEN 
           WRITE(*,535) 'a,e ',el(nsol)%coord(1),el(nsol)%coord(2) 
        ELSEIF(el(nsol)%coo.EQ.'COM'.or.elk%coo.eq.'COT') THEN 
           WRITE(*,535) 'q,e ',el(nsol)%coord(1),el(nsol)%coord(2)
        ELSEIF (el(nsol)%coo.EQ.'ATT')THEN
           WRITE(iun_log,535)'rho,rhodot',el(nsol)%coord(5),el(nsol)%coord(6)
        ENDIF
     ENDIF
525  FORMAT(12X,'ROOT NO.',I2,1x,I2) 
535  FORMAT(16X,'Preliminary orbit: ',A,' =',F10.5,';  ecc =',F10.5) 
20 END DO
                                                                        
  IF(nsol.LE.0) THEN 
     msg='No acceptable solution' 
     RETURN 
  END IF
                                                                        
  fail=.false. 
                                                                        
END SUBROUTINE laplace_poincare

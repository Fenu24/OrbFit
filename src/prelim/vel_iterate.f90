! ============================================================
! VEL_ITERATE
! computes velocity from positions at two times, given a first
! approximation of velocity at one of the two.
! also provides f,g series, thus completely replaces modgeq
! =============================================================

SUBROUTINE vel_iterate(gk,x1,x2,v20,dt,f,g,v2)
 USE ever_pitkin
 IMPLICIT NONE
 DOUBLE PRECISION, INTENT(IN) :: gk, dt ! sqrt(GM), time difference
 DOUBLE PRECISION, INTENT(IN) :: x1(3),x2(3),v20(3) ! pos1, pos2, vel2(approx)
 DOUBLE PRECISION, INTENT(OUT) :: v2(3),f,g ! improved velocity, f,g series 
! end interface
 DOUBLE PRECISION mu,s0,s1,s2,s3,alpha,e0,psi ! GM, Stumpff functions, 2*energy, ecc, univ. mean anomaly
 DOUBLE PRECISION prscal,sig0, vsize, r2, v, eps, q, energy, psi1
 LOGICAL accept, ecc_cont
 mu=gk**2
 sig0=prscal(x2,v20)
 r2=vsize(x2)
! v=vsize(v20)
! alpha=v**2-2*mu/r2 ! 2*E
 accept=ecc_cont(x2,v20,mu,1.d0,1.d0,e0,q,energy)
 eps=1000*epsilon(1.d0) ! convergence control
! IF(abs(alpha-2.d0*energy).gt.eps)THEN
!    WRITE(*,*)' vel_iterate: probl. in computing energy ', alpha, energy, alpha-2.d0*energy
! ENDIF
 alpha=2.d0*energy
! CALL solve_kepuniv(dt,r2,sig0,mu,alpha,psi1,s0,s1,s2,s3,eps)
 CALL solve_kepuniv2(dt,r2,sig0,mu,alpha,e0,psi,s0,s1,s2,s3,eps)
! IF(abs(psi1-psi).gt.eps)THEN
!    WRITE(*,*)' vel_iterate: probl. in computing psi ', psi, psi1, psi-psi1
! ENDIF 
 f=1.d0-mu*s2/r2
 g=dt-mu*s3
 v2=(x1-f*x2)/g
END SUBROUTINE vel_iterate

! Copyright (C) 1998-2000 by Mario Carpino (carpino@brera.mi.astro.it)  
! Version: November 21, 2006 (Milani) 
! --------------------------------------------------------------------- 
!                                                                       
!  *****************************************************************    
!  *                                                               *    
!  *                         G A U S S N                           *    
!  *                                                               *    
!  *       Initial orbit determination with Gauss' method          *    
!  *                                                               *    
!  *****************************************************************    
!                                                                       
! INPUT:    TOBS      -  Time of observations (MJD, TDT)                
!           ALPHA     -  Right ascension (rad)                          
!           DELTA     -  Declination (rad)                              
!           OBSCOD    -  Observatory code                               
!           DEBUG     -  Print debug information                        
!           MULTI     -  Multi-line output                              
!                                                                       
! OUTPUT:   ELEM      -  Orbital elements (mean ecliptic                
!                                          and equinox J2000)           
!           ELETYP    -  Type of orbital elements                       
!           T0        -  Epoch of orbital elements                      
!           NROOTS    -  Number of positive roots of 8th degree pol.    
!           NSOL      -  Number of solutions (some roots may be         
!                        discarded if they lead to solution with        
!                        eccentricity > ECCMAX)  
!           RTOP      -  topocentric distance at central obs                       
!           FAIL      -  Error flag                                     
!           MSG       -  Error message                                  
!                                                                       
SUBROUTINE gaussn(tobs,alpha,delta,obscod,elem,eletyp,            &
     &                  t0,nroots,nsol,rtop,fail,msg,debug,multi)            
  USE reference_systems 
  USE fund_const
  USE output_control
  IMPLICIT NONE                                                                       
  INCLUDE 'parobx.h90' 
  INTEGER, INTENT(IN) :: obscod(3)
  DOUBLE PRECISION, INTENT(IN) ::  tobs(3),alpha(3),delta(3) 
  LOGICAL, INTENT(IN) :: debug,multi 
  DOUBLE PRECISION, INTENT(OUT) :: elem(6,3),t0(3), rtop(3) ! dimension should be 3, not 8
  CHARACTER*(*),INTENT(OUT) :: msg,eletyp(3) 
  LOGICAL, INTENT(OUT) :: fail
  INTEGER, INTENT(OUT) ::  nroots,nsol                                                                        
! NEEDED common blocks:                                                 
  INCLUDE 'comgau.h90' 
! end interface                                                                        
  DOUBLE PRECISION, PARAMETER :: ecclim=5.d0
  INTEGER i,j,k,ising,ir,it 
  DOUBLE PRECISION xt(3,3),sinv0(3,3),a(3),b(3),c(3) 
  DOUBLE PRECISION ra(3),rb(3),coef(0:8),a2star,b2star,r22,s2r2 
  DOUBLE PRECISION esse0(3,3),cosd,det,tau1,tau3,tau13 
!      DOUBLE PRECISION esse(3,3),esse1(3,3),sinv(3,3)
  DOUBLE PRECISION roots(8),r2m3
  DOUBLE PRECISION gcap(3),crhom(3),rho(3),xp(3,3),vp(3),xv(6) 
!  DOUBLE PRECISION vekp(6) 
  DOUBLE PRECISION v1(3),v2(3),vs,err,sca 
  DOUBLE PRECISION xp1(3,3),tis2   !,rot(3,3) 
  DOUBLE PRECISION xve(6),vekpe(6),fs1,gs1,fs3,gs3,fggf 
  DOUBLE PRECISION f1,f3,g1,g3, v21(3),v23(3)
  CHARACTER*4 eltype 
  DOUBLE PRECISION, PARAMETER :: ecc_max=2.d0, qmax=1.d3 
  LOGICAL accept8,accept, ecc_cont
  DOUBLE PRECISION ecc,q,energy
  DOUBLE PRECISION vsize,princ 
  EXTERNAL vsize,princ 
!                                                                       
  IF(iicgau.NE.36) STOP '**** gaussn: internal error (01) ****' 
  fail=.true. 
  msg=' ' 
  nroots=0 
  nsol=0 
                                                                        
! COMPUTATION OF PRELIMINARY ORBIT                                      
!                                                                       
! ESSE = unit vector pointing in the direction of observations          
  DO  k=1,3 
     cosd=COS(delta(k)) 
     esse0(1,k)=cosd*COS(alpha(k)) 
     esse0(2,k)=cosd*SIN(alpha(k)) 
     esse0(3,k)=SIN(delta(k)) 
! Position of the observer   xt(1:3,k) for time tobs(k)
     CALL posobs(tobs(k),obscod(k),1,xt(1,k)) 
  END DO
                                                                        
! Inverse of ESSE matrix                                                
  sinv0=esse0
  CALL matin(sinv0,det,3,0,3,ising,1) 
  IF(ising.NE.0) THEN 
     IF(multi) THEN 
        msg='Singular S matrix (coplanar orbits?)' 
     ELSE 
        msg='SingS' 
     END IF
     RETURN 
  END IF
! central time
  tis2=tobs(2) 
                                                                        
! A and B vectors                                                       
  tau1=gk*(tobs(1)-tis2) 
  tau3=gk*(tobs(3)-tis2) 
  tau13=tau3-tau1 
  a(1)= tau3/tau13 
  a(2)=-1.d0 
  a(3)=-(tau1/tau13) 
  b(1)=a(1)*(tau13**2-tau3**2)/6.d0 
  b(2)=0.d0 
  b(3)=a(3)*(tau13**2-tau1**2)/6.d0 
                                                                        
! Coefficients of 8th degree equation for r2                            
  ra=MATMUL(xt,a) ! CALL prodmv(ra,xt,a) 
  rb=MATMUL(xt,b) ! CALL prodmv(rb,xt,b) 
  a2star=sinv0(2,1)*ra(1)+sinv0(2,2)*ra(2)+sinv0(2,3)*ra(3) 
  b2star=sinv0(2,1)*rb(1)+sinv0(2,2)*rb(2)+sinv0(2,3)*rb(3) 
  r22=xt(1,2)**2+xt(2,2)**2+xt(3,2)**2 
  s2r2=esse0(1,2)*xt(1,2)+esse0(2,2)*xt(2,2)+esse0(3,2)*xt(3,2) 
  coef=0.d0
  coef(8)=1.d0 
  coef(6)=-(a2star**2)-r22-(2.d0*a2star*s2r2) 
  coef(3)=-(2.d0*b2star*(a2star+s2r2)) 
  coef(0)=-(b2star**2) 
  IF(coef(0).EQ.0) THEN 
     msg='coef(0)=0' 
     RETURN 
  END IF
  CALL solv8(coef,roots,nroots) 
  IF(debug) THEN 
     IF(nroots.LE.0) THEN 
        WRITE(iun_log,524) 
524     FORMAT(16X,'roots none')
     ELSE 
!        WRITE(iun_log,514)(i,roots(i),i=1,nroots)
!514     FORMAT(16X,'r(',i1,')  =',F10.6) 
     END IF
  END IF
! 520 FORMAT(12X,'Coefficients of 8-th degree polynomial:') 
! 521 FORMAT(16X,'coef(',i1,')  =',F18.10) 
! 513 FORMAT(12X,'Possible roots for Sun-asteroid distance at ',        &
!     &       'observation',I2,':')                                      
  IF(nroots.LE.0) THEN 
     IF(multi) THEN 
        msg='8th degree polynomial has no real roots' 
     ELSE 
        msg='NoPolSol' 
     END IF
     RETURN 
  END IF
                                                                        
! Orbital elements of preliminary solution                              
  DO 20 ir=1,nroots 
     r2m3=1.d0/(roots(ir)**3) 
     c(1)=a(1)+b(1)*r2m3 
     c(2)=-1.d0 
     c(3)=a(3)+b(3)*r2m3 
     gcap=MATMUL(xt,c)        !  CALL prodmv(gcap,xt,c) 
     crhom=MATMUL(sinv0,gcap) !    CALL prodmv(crhom,sinv,gcap) 
     DO k=1,3 
       rho(k)=-(crhom(k)/c(k)) 
     ENDDO
                                                                        
! Position of the asteroid at the time of observations                  
     DO k=1,3 
        xp(1:3,k)=xt(1:3,k)+rho(k)*esse0(1:3,k) 
     ENDDO
                         
! remove spurious root WARNING: no theory
     IF(rho(2).lt.0.01d0)THEN
        IF(debug)WRITE(iun_log,*) ' spurious root , r, rho ',roots(ir),rho(2)
        CYCLE
     ELSE
        IF(debug)WRITE(iun_log,*) ' accepted root , r, rho ',roots(ir),rho(2)
     ENDIF
                                               
! Gibbs' transformation, giving the velocity of the planet at the       
! time of second observation                                            
     CALL gibbs(xp,tau1,tau3,vp,gk) 
! Cartesian coordinates of preliminary orbit                                 
     xv(1:3)=xp(1:3,2) 
     xv(4:6)=vp(1:3) 
! cartesian ecliptic coordinates 
     xve(1:3)=MATMUL(roteqec,xv(1:3)) 
     xve(4:6)=MATMUL(roteqec,xv(4:6))
! hyperbolic control
     accept8=ecc_cont(xp(1:3,2),vp,gms,ecc_max,qmax,ecc,q,energy)
! orbital elements
     CALL ccek1(vekpe,eltype,xve,gms) 
     IF(debug) THEN 
        WRITE(iun_log,525) ir 
525     FORMAT(12X,'ROOT NO.',I2) 
        IF(eltype.EQ.'KEP') THEN 
           WRITE(iun_log,535) 'a',vekpe(1),vekpe(2)
535        FORMAT(16X,'Preliminary orbit: ',A,' =',F10.5,';  ecc =',F10.5)
        ELSEIF(eltype.EQ.'COM') THEN 
           WRITE(iun_log,535) 'q',vekpe(1),vekpe(2) 
        ELSE 
           STOP '**** gaussn: internal error (02) ****' 
        END IF
     END IF
     
! store output of deg. 8 equation, in case the iterations fail
     IF(accept8)THEN   
        rtop(nsol+1)=rho(2)
        elem(1:6,nsol+1)=vekpe(1:6) 
        eletyp(nsol+1)=eltype
! correction for aberration 
        t0(nsol+1)=tis2 -rho(2)/vlight
     ELSE
        CYCLE ! try next root
     ENDIF
                                                                  
! CORRECTION OF PRELIMINARY ORBIT                                       
     it=0 
                                                                
30   CONTINUE 
!     CALL modgeq(gk,xp(1,2),xp(1,1),tau1,vekp,eltyp,fs1,gs1,v1) 
!     CALL modgeq(gk,xp(1,2),xp(1,3),tau3,vekp,eltyp,fs3,gs3,v2) 
     CALL vel_iterate(gk,xp(1:3,1),xv(1:3),xv(4:6),tobs(1)-tobs(2),fs1,gs1,v1)
     CALL vel_iterate(gk,xp(1:3,3),xv(1:3),xv(4:6),tobs(3)-tobs(2),fs3,gs3,v2)
     it=it+1 
     DO i=1,3 
        vp(i)=(v1(i)+v2(i))/2.d0
     ENDDO
!     CALL vel_iterate(gk,xp(1:3,1),xv(1:3),xv(4:6),tobs(1)-tobs(2),f1,g1,v21)
!     WRITE(iun_log,*)' diff1 ',f1-fs1,g1-gs1, v21-v1
!     CALL vel_iterate(gk,xp(1:3,3),xv(1:3),xv(4:6),tobs(3)-tobs(2),f3,g3,v23)
!     WRITE(iun_log,*)' diff3 ',f3-fs3,g3-gs3, v23-v2
     fggf=fs1*gs3-fs3*gs1 
     c(1)=gs3/fggf 
     c(3)=-(gs1/fggf) 
                                                                        
! Updated set of orbital elements                                       
     gcap=MATMUL(xt,c) ! CALL prodmv(gcap,xt,c) 
     crhom=MATMUL(sinv0,gcap)   ! CALL prodmv(crhom,sinv0,gcap) 
     DO k=1,3 
        rho(k)=-(crhom(k)/c(k)) 
     END DO
     DO k=1,3 
        xp1(1:3,k)=xt(1:3,k)+rho(k)*esse0(1:3,k) 
     ENDDO
                                                                        
! Error with respect to previous iteration                              
     err=0.d0 
     sca=0.d0 
     DO i=1,3
        DO k=1,3 
           err=err+(xp1(i,k)-xp(i,k))**2 
           sca=sca+xp1(i,k)**2 
        ENDDO
     ENDDO
     xp=xp1
     err=SQRT(err/sca)

     xv(1:3)=xp(1:3,2) 
     xv(4:6)=vp(1:3)
! hyperbolic control
     accept=ecc_cont(xp(1:3,2),vp,gms,ecc_max,qmax,ecc,q,energy)
!     WRITE(*,*)' err=',err,' iter=',it, ' ecc,q,E,rho(2)', ecc,q, energy,rho(2)  
     IF(.not.accept)THEN
        WRITE(iun_log,*)' divergent iteration', it, ' ecc,q,E,rho(2)', ecc,q, energy, rho(2)
        IF(accept8)nsol=nsol+1
        CYCLE
     ENDIF

! Convergency control                                                   
     IF(err.LE.errmax) GOTO 26 ! convergent
     IF(it.LE.itmax) GOTO 30 ! another iteration
! iterations have not converged
     IF(debug) WRITE(iun_log,543)it, err
543  FORMAT(16X,'WARNING: after ',i4,' iterations not converging, err=',1P,D10.3) 
     nsol=nsol+1
     CYCLE ! the solution of deg. 8 equation gets stored instead

! Final orbit from convergent iterations
26   CONTINUE 
! cartesian ecliptic coordinates 
     xve(1:3)=MATMUL(roteqec,xv(1:3)) ! CALL prodmv(xve(1),roteqec,xv(1)) 
     xve(4:6)=MATMUL(roteqec,xv(4:6)) ! CALL prodmv(xve(4),roteqec,xv(4))
! ecliptic elements vekpe
     CALL ccek1(vekpe,eltype,xve,gms)
! Check for anomalous orbits already done                                           
! store good result
     nsol=nsol+1 
     rtop(nsol)=rho(2)
     elem(1:6,nsol)=vekpe(1:6) 
     eletyp(nsol)=eltype
! correction for aberration 
     t0(nsol)=tis2 -rho(2)/vlight
     IF(debug) THEN 
        IF(eltype.EQ.'KEP') THEN 
           WRITE(iun_log,544) 'a',vekpe(1),vekpe(2),it 
        ELSE 
           WRITE(iun_log,544) 'q',vekpe(1),vekpe(2),it 
        END IF
     END IF
544  FORMAT(16X,'Final orbit:       ',A,' =',F10.5,';  ecc =',F10.5,   &
              &       ' (',I3,' iterations)') 
20 END DO
                                                                        
  IF(nsol.LE.0) THEN 
     IF(multi) THEN 
        msg='No acceptable solution' 
     ELSE 
        msg='NoAccSol' 
     END IF
     RETURN 
  END IF
                                                                        
  fail=.false. 
  IF(.NOT.multi) THEN 
     msg='OK' 
  END IF
                                                                        
END SUBROUTINE gaussn

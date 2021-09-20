! Copyright (C) 1998-2000 by Mario Carpino (carpino@brera.mi.astro.it)  
! Version 3.3.2, A. Milani 2006                                               
! --------------------------------------------------------------------- 
!                                                                       
!  *****************************************************************    
!  *                                                               *    
!  *                         G A U S S D E G 8                     *    
!  *                                                               *    
!  *       Initial orbit determination with Gauss' method          *
!  *       Limited to deg 8 dynamical equation                     *    
!  *                                                               *    
!  *****************************************************************    
!                                                                       
! INPUT:    TOBS      -  Time of observations (MJD, TDT)                
!           ALPHA     -  Right ascension (rad)                          
!           DELTA     -  Declination (rad)                              
!           OBSCOD    -  Observatory code                               
!           DEBUG     -  Print debug information                        
!                                                                       
! OUTPUT:   EL        -  Orbital elements (cartesian, mean ecliptic
!                                          and equinox J2000)           
!           NROOTS    -  Number of positive roots of 8th degree pol.    
!           NSOL      -  Number of solutions (some roots may be         
!                        discarded if they lead to solution with        
!                        eccentricity > ECCMAX)                         
!           FAIL      -  Error flag                                     
!           MSG       -  Error message                                  
!                                                                       
SUBROUTINE gaussdeg8(tobs,alpha,delta,obscod,ecc_max,q_max,el,nroots,nsol,rr,fail,msg,debug)
  USE reference_systems 
  USE fund_const
  USE orbit_elements
  USE output_control
  USE planet_masses, ONLY: gmearth
  USE force_sat, ONLY: eamoon_mass
  IMPLICIT NONE  
  INTEGER, INTENT(IN) ::  obscod(3) 
  DOUBLE PRECISION, INTENT(IN):: tobs(3),alpha(3),delta(3) 
  DOUBLE PRECISION, INTENT(IN) :: ecc_max, q_max
  TYPE(orbit_elem), INTENT(OUT) :: el(3)
  DOUBLE PRECISION, INTENT(OUT) :: rr(3) ! topocentric distance
  LOGICAL, INTENT(OUT) :: fail
  LOGICAL, INTENT(IN) :: debug
  CHARACTER*(*), INTENT(OUT) :: msg
  INTEGER, INTENT(OUT) :: nsol,nroots
! end interface
  DOUBLE PRECISION elem(6,8),t0(8) 
  INTEGER i,j,k,ising,ir! ,it 
  DOUBLE PRECISION xt(3,3),sinv0(3,3),a(3),b(3),c(3), xea(6)
  DOUBLE PRECISION ra(3),rb(3),coef(0:8),a2star,b2star,r22,s2r2 
  DOUBLE PRECISION esse0(3,3),cosd,det,tau1,tau3,tau13 
  DOUBLE PRECISION roots(8),r2m3 
  DOUBLE PRECISION gcap(3),crhom(3),rho(3),xp(3,3),vp(3),xv(6) ,tis2
  DOUBLE PRECISION l0,h0,coseps,sineps,vvv(3),vsize,r0
  DOUBLE PRECISION rootgm, rootgm2
  INTEGER fail_flag ! for coordinate change
  TYPE(orbit_elem) elk 
  DOUBLE PRECISION ecc,q,energy, spurious_dist, rho_scal, r_scal
  LOGICAL accept8, ecc_cont
  DOUBLE PRECISION vt(3,3) ! auxiliary for observer_position
! ===============================================================
  fail=.true. 
  msg=' ' 
  nroots=0 
  nsol=0 
! assignement of rootgm
  IF(rhs.eq.1)THEN
     rootgm=gk
  ELSEIF(rhs.EQ.2)THEN
     CALL eamoon_mass
     rootgm=sqrt(gmearth)
  ELSE
     WRITE(*,*) ' gaussdeg8: wrong rhs=', rhs
     STOP
  END IF                                                                        
! COMPUTATION OF PRELIMINARY ORBIT                                      
!                                                                       
! ESSE = unit vector pointing in the direction of observations          
  DO k=1,3 
     cosd=COS(delta(k)) 
     esse0(1,k)=cosd*COS(alpha(k)) 
     esse0(2,k)=cosd*SIN(alpha(k)) 
     esse0(3,k)=SIN(delta(k)) 
! Position of the observer xt(1:3,i) for time tobs(i)
     DO i=1,3 
        IF(rhs.EQ.1)THEN
! heliocentric ecliptic observer position 
           CALL earcar(tobs(i),xea,1)
           CALL observer_position(tobs(i),xt(1:3,i),vt(1:3,i),OBSCODE=obscod(i)) 
           xt(1:3,i)=xt(1:3,i)+xea(1:3) 
           CALL prodmv(xt(1:3,i),roteceq,xt(1:3,i))
       ELSEIF(rhs.EQ.2)THEN
! geocentric equatorial observer position
           CALL observer_position(tobs(i),xt(1:3,i),vt(1:3,i),OBSCODE=obscod(i)) 
! conversion to equatorial: done in observer_position
        ELSE
           STOP
        END IF
     END DO
  END DO
! Inverse of ESSE matrix                                                
  sinv0=esse0 
  CALL matin(sinv0,det,3,0,3,ising,1) 
  IF(ising.NE.0) THEN 
     msg='Singular S matrix (coplanar orbits?)' 
     RETURN 
  END IF
  tis2=tobs(2) 
! A and B vectors                                                       
  tau1=rootgm*(tobs(1)-tis2) ! mu*t_12
  tau3=rootgm*(tobs(3)-tis2) ! mu*t_32
  tau13=tau3-tau1        ! mu*t_31
  a(1)= tau3/tau13       
  a(2)=-1.d0 
  a(3)=-(tau1/tau13) 
  b(1)=a(1)*(tau13**2-tau3**2)/6.d0  
  b(2)=0.d0 
  b(3)=a(3)*(tau13**2-tau1**2)/6.d0 
! Coefficients of 8th degree equation for r2                            
  CALL prodmv(ra,xt,a) 
  CALL prodmv(rb,xt,b) 
  a2star=sinv0(2,1)*ra(1)+sinv0(2,2)*ra(2)+sinv0(2,3)*ra(3) 
  b2star=sinv0(2,1)*rb(1)+sinv0(2,2)*rb(2)+sinv0(2,3)*rb(3) 
  r22=xt(1,2)**2+xt(2,2)**2+xt(3,2)**2 ! q_2^2 
  s2r2=esse0(1,2)*xt(1,2)+esse0(2,2)*xt(2,2)+esse0(3,2)*xt(3,2) ! <q_2,rhat_2> 
  coef(8)=1.d0 
  coef(7)=0.d0 
  coef(6)=-(a2star**2)-r22-(2.d0*a2star*s2r2) 
  coef(5)=0.d0 
  coef(4)=0.d0 
  coef(3)=-(2.d0*b2star*(a2star+s2r2)) 
  coef(2)=0.d0 
  coef(1)=0.d0 
  coef(0)=-(b2star**2) 
  IF(coef(0).EQ.0) THEN 
     msg='coef(0)=0' 
     RETURN 
  END IF
  h0=-a2star/b2star*r22**1.5d0
  l0=-1.d0/b2star
  coseps=s2r2/sqrt(r22)
  CALL prvec(xt(1:3,2),esse0(1:3,2),vvv)
  sineps=vsize(vvv)/sqrt(r22)
  r0=sqrt(r22)/(h0**(1.d0/3.d0))
  CALL solv8(coef,roots,nroots) 
  IF(nroots.LE.0) THEN 
     WRITE(iun_log,524)
     WRITE(iun_log,515)l0,h0,coseps,sqrt(r22),sineps,r0
  ELSE 
!     WRITE(iun_log,514)(i,roots(i),i=1,nroots)
     WRITE(iun_log,515)l0,h0,coseps,sqrt(r22),sineps,r0
  END IF
  IF(debug) THEN 
     IF(nroots.LE.0) THEN 
        WRITE(*,524)
        WRITE(*,515)l0,h0,coseps,sqrt(r22),sineps,r0
 524    FORMAT(16X,'no roots of deg.8 equation') 
     ELSE 
         WRITE(*,514)(i,roots(i),i=1,nroots)
 514     FORMAT(16X,'r(',i1,')  =',F10.6)
         WRITE(*,515)l0,h0,coseps,sqrt(r22),sineps,r0
 515     FORMAT('l0=',1P,D12.5,' h0=',D12.5,' cos(eps)=',D12.5,' q=',D12.5,' sineps=',D12.5,' r0=',D12.5) 
     END IF
  END IF                                                                        
  IF(nroots.LE.0) THEN 
     msg='8th degree polynomial has no real roots' 
     RETURN 
  ELSEIF(nroots.gt.3)THEN
     WRITE(*,*)' gaussdeg8: more than 3 positive roots ',nroots
     STOP
  END IF
                                                                        
! Orbital elements of preliminary solution
  nsol=0                              
  DO 20 ir=1,nroots 
     r2m3=1.d0/(roots(ir)**3) 
     c(1)=a(1)+b(1)*r2m3 
     c(2)=-1.d0 
     c(3)=a(3)+b(3)*r2m3 
     CALL prodmv(gcap,xt,c) 
     CALL prodmv(crhom,sinv0,gcap) 
     DO 13 k=1,3 
        rho(k)=-(crhom(k)/c(k)) 
13   END DO                                         
! Position of the asteroid at the time of observations
     DO 14 k=1,3 
        xp(1:3,k)=xt(1:3,k)+rho(k)*esse0(1:3,k) 
14   ENDDO
! remove spurious root
     IF(rhs.eq.1)THEN
        spurious_dist=0.001D0 ! Hill's sphere of the Earth
        r_scal=roots(ir) ! output message in AU
        rho_scal=rho(2) ! output message in AU
     ELSEIF(rhs.eq.2)THEN
        spurious_dist=1.D-6 ! 150 km, Earth atmosphere
        r_scal=roots(ir)*aukm ! output message in km
        rho_scal=rho(2)*aukm ! output message in km
     ELSE
        STOP
     ENDIF
     IF(rho(2).lt.spurious_dist)THEN
        IF(debug)WRITE(*,*) ' spurious root ',ir,' , r, rho ',r_scal,rho_scal
        WRITE(iun_log,*) ' spurious root ',ir,', r, rho ',r_scal,rho_scal
        CYCLE
     ELSE
        IF(debug)WRITE(*,*) ' accepted root ',ir,' , r, rho ',r_scal,rho_scal
        WRITE(iun_log,*) ' accepted root ',ir,' , r, rho ',r_scal,rho_scal
     ENDIF
! Gibbs' transformation, giving the velocity of the planet at the       
! time of second observation
     CALL gibbs(xp,tau1,tau3,vp,rootgm) 
     IF(rhs.EQ.1)THEN
! conversion to ecliptic coordinates
        CALL prodmv(xv(1:3),roteqec,xp(1:3,2))
        CALL prodmv(xv(4:6),roteqec,vp)
     ELSEIF(rhs.eq.2)THEN
! no need of conversion
        xv(1:3)=xp(1:3,2)
        xv(4:6)=vp
     ELSE
        STOP
     END IF
! hyperbolic control
     rootgm2=rootgm**2
     accept8=ecc_cont(xp(1:3,2),vp,rootgm2,ecc_max,q_max,ecc,q,energy)
     IF(accept8)THEN
        nsol=nsol+1
        rr(nsol)=rho(2)
     ELSE
        IF(debug)WRITE(*,*) ' hyperbolic, ecc,q  ',ecc,q
        WRITE(iun_log,*)' hyperbolic, ecc,q  ',ecc,q 
        CYCLE ! try next root
     ENDIF  
! Orbital elements of preliminary orbit                                 
     el(nsol)=undefined_orbit_elem
     IF(rhs.eq.1)THEN
        el(nsol)%center=0
     ELSEIF(rhs.EQ.2)THEN
        el(nsol)%center=3
     ELSE
        STOP
     END IF
     el(nsol)%coord=xv
     el(nsol)%coo='CAR'
! planetary aberration 
     el(nsol)%t=tobs(2)-rho(2)/vlight
     IF(debug) THEN 
        CALL coo_cha(el(nsol),'KEP',elk,fail_flag)
        WRITE(iun_log,525) ir, fail_flag 
        IF(elk%coo.EQ.'KEP') THEN 
           WRITE(iun_log,535) 'a',elk%coord(1),elk%coord(2) 
        ELSEIF(elk%coo.EQ.'COM'.or.elk%coo.eq.'COT') THEN 
           WRITE(iun_log,535) 'q',elk%coord(1),elk%coord(2)
        ELSE 
           STOP '**** gaussdeg8: unknown elem type ****' 
        END IF
        WRITE(*,525) ir, fail_flag 
        IF(elk%coo.EQ.'KEP') THEN 
           WRITE(*,535) 'a',elk%coord(1),elk%coord(2) 
        ELSEIF(elk%coo.EQ.'COM'.or.elk%coo.eq.'COT') THEN 
           WRITE(*,535) 'q',elk%coord(1),elk%coord(2)
        END IF
     END IF
525  FORMAT(12X,'ROOT NO.',I2,1x,I2) 
535  FORMAT(16X,'Preliminary orbit: ',A,' =',F10.5,';  ecc =',F10.5) 
                                                                        
                                                                        
20 END DO
                                                                        
  IF(nsol.LE.0) THEN 
     msg='No acceptable solution' 
     RETURN 
  END IF
                                                                        
  fail=.false. 
                                                                        
END SUBROUTINE gaussdeg8

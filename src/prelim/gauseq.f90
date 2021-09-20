! Copyright (C) 1998 by Mario Carpino (carpino@brera.mi.astro.it)       
! Version: December 11, 1998                                            
! --------------------------------------------------------------------- 
!                                                                       
!  *****************************************************************    
!  *                                                               *    
!  *                         G A U S E Q                           *    
!  *                                                               *    
!  *         Iterative solution of Gauss' equation for Y           *    
!  *                  (triangle to sector ratio)                   *    
!  *                                                               *    
!  *****************************************************************    
!                                                                       
!                                                                       
! INPUT:    GK        -  Gauss' k constant (possibly multiplied by      
!                        the square root of the sum of the masses)      
!           X1        -  Position vector at time t1                     
!           X2        -  Position vector at time t2                     
!           TAU       -  Normalized time difference ( = k*(t2-t1) )     
!           AINC      -  Orbital inclination (rad)                      
!           Y         -  Starting value for Y                           
!                                                                       
! OUTPUT:   Y         -  Corrected value for Y                          
!           VEL1      -  Orbital velocity at time t1                    
!           IRET      -  Return code:                                   
!                           0 = OK                                      
!                           1 = not converged                           
!                           2 = RAD1 < 0                                
!                                                                       
      SUBROUTINE gauseq(gk,x1,x2,tau,ainc,y,vel1,iret) 
      IMPLICIT NONE 
                                                                        
      DOUBLE PRECISION gk,x1(3),x2(3),tau,ainc,y,vel1(3) 
      INTEGER iret 
                                                                        
! TUNING PARAMETERS                                                     
      DOUBLE PRECISION toll 
      PARAMETER (toll=1.d-14) 
      INTEGER itmax 
      PARAMETER (itmax=50) 
                                                                        
      DOUBLE PRECISION r1,r2,cosv21,sinv21,pi2,v21,rad,rad1,el,em,a 
      DOUBLE PRECISION cose2,sine2,eta,xlow,xup,ynew,corr,rada 
      DOUBLE PRECISION fser,gser 
      INTEGER it,i 
                                                                        
      iret=0 
      r1=SQRT(x1(1)**2+x1(2)**2+x1(3)**2) 
      r2=SQRT(x2(1)**2+x2(2)**2+x2(3)**2) 
      cosv21=(x1(1)*x2(1)+x1(2)*x2(2)+x1(3)*x2(3))/(r1*r2) 
      pi2=2.d0*ATAN(1.d0) 
      sinv21=x1(1)*x2(2)-x2(1)*x1(2) 
      sinv21=sinv21*SQRT(1.d0-cosv21**2)/dabs(sinv21) 
      IF(ainc.GT.pi2) sinv21=-sinv21 
      v21=ATAN2(sinv21,cosv21) 
      cosv21=dcos(v21/2.d0) 
      rad=SQRT(r1*r2) 
      el=(r1+r2)/(4.d0*rad*cosv21)-0.5d0 
      em=(tau**2)/((2.d0*rad*cosv21)**3) 
                                                                        
      it=0 
    1 CONTINUE 
      xlow=em/(y**2)-el 
      cose2=1.d0-2.d0*xlow 
      rad1=xlow*(1.d0-xlow) 
      IF(rad1.LT.0.d0) THEN 
          iret=2 
          RETURN 
      END IF 
      sine2=2.d0*SQRT(rad1) 
      eta=2.d0*ATAN2(sine2,cose2) 
      xup=(eta-SIN(eta))/(sine2**3) 
      ynew=1.d0+xup*(el+xlow) 
      corr=(ynew-y)/ynew 
      corr=dabs(corr) 
      it=it+1 
      y=ynew 
      IF(it.GE.itmax) THEN 
          iret=1 
          GOTO 2 
      END IF 
      IF(corr.GT.toll) GOTO 1 
                                                                        
    2 CONTINUE 
      rada=tau/(2.d0*y*rad*cosv21*sine2) 
      a=rada**2 
      fser=1.d0-(1.d0-dcos(eta))*a/r1 
      gser=tau-(eta-SIN(eta))*a*rada 
      gser=gser/gk 
      DO 3 i=1,3 
      vel1(i)=(x2(i)-fser*x1(i))/gser 
    3 END DO 
                                                                        
      END                                           

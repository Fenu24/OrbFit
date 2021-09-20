! Copyright (C) 1998 by Mario Carpino (carpino@brera.mi.astro.it)       
! Version: December 11, 1998                                            
! 13/07/2016 - Iorfida: scaling factor g removed from d_i
! --------------------------------------------------------------------- 
!                                                                       
!  *****************************************************************    
!  *                                                               *    
!  *                          G I B B S                            *    
!  *                                                               *    
!  *     Gibbs' algorithm for computing position and velocity      *    
!  *                from three position vectors                    *    
!  *                                                               *    
!  *****************************************************************    
!                                                                       
! INPUT:    R(i,k)    -  Components (i=1,3) of position vectors at three
!                        different times (k=1,3)                        
!           TAU1      -  Normalized time at instant k=1 ( = k*(t1-t2) ) 
!           TAU3      -  Normalized time at instant k=3 ( = k*(t3-t2) ) 
!           GK        -  Gauss' k constant (possibly multiplied by      
!                        the square root of the sum of the masses)      
!                                                                       
! OUTPUT:   V2        -  Velocity vector at instant k=2                 
!                                                                       
      SUBROUTINE gibbs(r,tau1,tau3,v2,gk) 
      IMPLICIT NONE 
                                                                        
      DOUBLE PRECISION r(3,3),tau1,tau3,v2(3),gk 
                                                                        
      DOUBLE PRECISION tau13,g,r1m3,r2m3,r3m3,d1,d2,d3 
      INTEGER i 
                                                                        
      tau13=tau3-tau1 
      g=gk**2 
      r1m3=1.d0/(SQRT(r(1,1)**2+r(2,1)**2+r(3,1)**2)**3) 
      r2m3=1.d0/(SQRT(r(1,2)**2+r(2,2)**2+r(3,2)**2)**3) 
      r3m3=1.d0/(SQRT(r(1,3)**2+r(2,3)**2+r(3,3)**2)**3) 
      d1=tau3*(r1m3/12.d0-1.d0/(tau1*tau13)) 
      d2=(tau1+tau3)*(r2m3/12.d0-1.d0/(tau1*tau3)) 
      d3=-tau1*(r3m3/12.d0+1.d0/(tau3*tau13)) 
      DO 1 i=1,3 
      v2(i)=gk*(-d1*r(i,1)+d2*r(i,2)+d3*r(i,3)) 
    1 END DO 
                                                                        
      END                                           

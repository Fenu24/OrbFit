! Copyright (C) 2000 by Mario Carpino (carpino@brera.mi.astro.it)       
! Version: June 9, 2000                                                 
! --------------------------------------------------------------------- 
!                                                                       
!  *****************************************************************    
!  *                                                               *    
!  *                          H Y P A N                            *    
!  *                                                               *    
!  *                Hyperbolic eccentric anomaly                   *    
!  *                                                               *    
!  *****************************************************************    
!                                                                       
! INPUT:    M         -  Mean anomaly (rad)                             
!           ECC       -  Eccentricity                                   
!                                                                       
! OUTPUT:   HYPAN     -  Hyperbolic eccentric anomaly (rad)             
!                                                                       
      DOUBLE PRECISION FUNCTION hypan(m,ecc) 
      IMPLICIT NONE 
                                                                        
      DOUBLE PRECISION m,ecc 
                                                                        
      DOUBLE PRECISION eps 
      INTEGER nmax 
      PARAMETER (eps=1.d-14,nmax=30) 
                                                                        
      INTEGER n 
      DOUBLE PRECISION y,w,err,df 
                                                                        
      IF(ecc.LT.1.d0) STOP '**** hypan: ECC < 1 ****' 
                                                                        
      IF(m.GT.5.d0*ecc-2.5d0) THEN 
          hypan=LOG(2.d0*m/ecc) 
      ELSE 
          y=SQRT(8.d0*(ecc-1.d0)/ecc) 
          w=3.d0*m/(y*(ecc-1.d0)) 
          w=w+SQRT(W**2+1.D0) 
          W=LOg(w)/3.d0 
          hypan=y*SINH(w) 
      END IF 
                                                                        
      n=0 
    1 CONTINUE 
      err=m+hypan-ecc*SINH(hypan) 
      df=err/(ecc*COSH(hypan)-1.d0) 
      n=n+1 
      hypan=hypan+df 
      IF(ABS(err).LT.eps) RETURN 
      IF(n.LE.nmax) GOTO 1 
                                                                        
!     WRITE(*,100) m,ecc                                                
!100  FORMAT(' **** HYPAN: not converged'/                              
!    .       ' **** M =',1P,E12.4,', ecc =',E12.4)                      
                                                                        
      END                                           

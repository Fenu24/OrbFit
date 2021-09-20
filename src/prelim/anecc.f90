! Copyright (C) 2000 by Mario Carpino (carpino@brera.mi.astro.it)       
! Version: May 29, 2000                                                 
! --------------------------------------------------------------------- 
!                                                                       
!  *****************************************************************    
!  *                                                               *    
!  *                          A N E C C                            *    
!  *                                                               *    
!  *             Iterative solution of Kepler's equation           *    
!  *                                                               *    
!  *****************************************************************    
!                                                                       
! INPUT:    M         -  Mean anomaly (rad)                             
!           ECC       -  Eccentricity                                   
!                                                                       
! OUTPUT:   ANECC     -  Eccentric anomaly (rad)                        
!                                                                       
      DOUBLE PRECISION FUNCTION anecc(m,ecc) 
      IMPLICIT NONE 
                                                                        
      DOUBLE PRECISION m,ecc 
                                                                        
      INTEGER nmax 
      DOUBLE PRECISION eps 
      PARAMETER (eps=1.d-14,nmax=50) 
                                                                        
      INTEGER n 
      DOUBLE PRECISION anold,delta 
                                                                        
      anold=m 
      n=0 
    1 CONTINUE 
      anecc=anold-(anold-ecc*SIN(anold)-m)/(1.d0-ecc*COS(anold)) 
      delta=anecc-ecc*SIN(anecc)-m 
      n=n+1 
      IF(ABS(delta).LT.eps) RETURN 
      anold=anecc 
      IF(n.LE.nmax) GOTO 1 
                                                                        
      WRITE(*,100) m,ecc 
  100 FORMAT(' **** ANECC: not converged'/                              &
     &       ' **** M =',F15.10,'     ecc =',F15.10)                    
                                                                        
      END                                           

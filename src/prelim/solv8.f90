! Copyright (C) 1998-2000 by Mario Carpino (carpino@brera.mi.astro.it)  
! Version: March 8, 2000                                                
! --------------------------------------------------------------------- 
!                                                                       
!  ***************************************************************      
!  *                                                             *      
!  *                          S O L V 8                          *      
!  *                                                             *      
!  *         Calcolo delle radici reali positive                 *      
!  *              di un polinomio di grado 8                     *      
!  *                                                             *      
!  ***************************************************************      
!                                                                       
!                                                                       
! INPUT:    COEF(i)   -  Coefficiente di z^i (i=0,8)                    
!                                                                       
! OUTPUT:   ROOTS(k)  -  Radici reali positive (k=1,nroots)             
!           NROOTS    -  Numero delle radici reali positive             
!                                                                       
!                                                                       
      SUBROUTINE solv8(coef,roots,nroots) 
      IMPLICIT NONE 
                                                                        
      DOUBLE PRECISION coef(0:8),roots(8) 
      INTEGER nroots 
                                                                        
      COMPLEX*16 poly(0:8),zr1(8) 
      DOUBLE PRECISION epsm,big,small,rad(8),a1(9),a2(9),rim,rre 
      INTEGER nit,j 
      LOGICAL err(9) 
                                                                        
      DO 10 j=0,8 
      poly(j)=coef(j) 
   10 END DO 
      epsm=2.D0**(-53) 
      big=2.D0**1023 
      small=2.D0**(-1022) 
      CALL polzeros(8,poly,epsm,big,small,50,zr1,rad,err,nit,a1,a2) 
                                                                        
      nroots=0 
      DO 11 j=1,8 
      rim=AIMAG(zr1(j)) 
      IF(ABS(rim).LT.rad(j)) THEN 
          rre=DBLE(zr1(j)) 
          IF(rre.LT.0.D0) GOTO 11 
          nroots=nroots+1 
          roots(nroots)=rre 
      END IF 
   11 END DO 
                                                                        
      END                                           

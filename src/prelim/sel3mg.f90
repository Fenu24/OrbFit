! Copyright (C) 2000 by Mario Carpino (carpino@brera.mi.astro.it)       
! Version: May 31, 2000                                                 
! --------------------------------------------------------------------- 
!                                                                       
!  *****************************************************************    
!  *                                                               *    
!  *                         S E L 3 M G                           *    
!  *                                                               *    
!  *  Selection of 3 observations for initial orbit determination  *    
!  *           (STEP 2: output of selected triplets)               *    
!  *                                                               *    
!  *****************************************************************    
!                                                                       
! OUTPUT:   SELIPT    -  Selected triplet                               
!           EOT       -  End of triplets                                
!                                                                       
      SUBROUTINE sel3mg(selipt,eot) 
      IMPLICIT NONE 
                                                                        
      INTEGER selipt(3) 
      LOGICAL eot 
                                                                        
      INCLUDE 'pars3m.h90' 
                                                                        
! NEEDED common blocks:                                                 
      INCLUDE 'sel3m.h90' 
                                                                        
      INTEGER k 
                                                                        
      IF(iics3m.NE.36) STOP '**** sel3mg: internal error (01) ****' 
                                                                        
      ipt3=ipt3+1 
      IF(ipt3.GT.n3s) THEN 
          eot=.true. 
          selipt(1)=0 
          selipt(2)=0 
          selipt(3)=0 
          iics3m=0 
          RETURN 
      END IF 
                                                                        
      k=ords3(ipt3) 
      selipt(1)=l3s(1,k) 
      selipt(2)=l3s(2,k) 
      selipt(3)=l3s(3,k) 
      eot=.false. 
                                                                        
      END                                           

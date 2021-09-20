! Copyright (C) 1998-2000 by Mario Carpino (carpino@brera.mi.astro.it)  
! Version: June 7, 2000                                                 
! --------------------------------------------------------------------- 
!                                                                       
!  *****************************************************************    
!  *                                                               *    
!  *                         G A U I N I                           *    
!  *                                                               *    
!  *            Input of options for Gauss' method                 *    
!  *                                                               *    
!  *****************************************************************    
!                                                                       
      SUBROUTINE gauini(fail) 
      IMPLICIT NONE 
                                                                        
      LOGICAL fail 
                                                                        
! Common blocks to be initialized:                                      
      INCLUDE 'comgau.h90' 
                                                                        
      LOGICAL found,fail1 
      fail=.false.                                                       
      errmax=1.D-10 
      CALL rdnrea('init_orbdet.gauss.','max_err',errmax,.false.,found,  &
     &            fail1,fail)                                           
                                                                        
      eccmax=10.D0 
      CALL rdnrea('init_orbdet.gauss.','max_ecc',eccmax,.false.,found,  &
     &            fail1,fail)                                           
                                                                        
      itmax=50 
      CALL rdnint('init_orbdet.gauss.','nit_max',itmax,.false.,found,   &
     &            fail1,fail)                                           
      iicgau=0                                                                  
      IF(.NOT.fail) iicgau=36 
                                                                        
      END                                           

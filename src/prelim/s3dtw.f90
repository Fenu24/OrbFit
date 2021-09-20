! Copyright (C) 2000 by Mario Carpino (carpino@brera.mi.astro.it)       
! Version: May 31, 2000                                                 
! --------------------------------------------------------------------- 
!                                                                       
!  *****************************************************************    
!  *                                                               *    
!  *                          S 3 D T W                            *    
!  *                                                               *    
!  *  Weight to be assigned to time interval between observations  *    
!  *  for the selection of triplets in initial orbit determination *    
!  *                                                               *    
!  *****************************************************************    
!                                                                       
! INPUT:    DT        -  Time interval                                  
!           DTW       -  Best time interval                             
!                                                                       
      DOUBLE PRECISION FUNCTION s3dtw(dt,dtw) 
      IMPLICIT NONE 
                                                                        
      DOUBLE PRECISION dt,dtw 
                                                                        
      IF(dt.LE.dtw) THEN 
          s3dtw=dtw/dt 
      ELSE 
          s3dtw=1+dt/dtw 
      END IF 
                                                                        
      END                                           

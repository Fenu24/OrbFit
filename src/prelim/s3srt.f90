! Copyright (C) 2000 by Mario Carpino (carpino@brera.mi.astro.it)       
! Version: May 31, 2000                                                 
! --------------------------------------------------------------------- 
!                                                                       
!  *****************************************************************    
!  *                                                               *    
!  *                          S 3 S R T                            *    
!  *                                                               *    
!  * Sort triplets of observations for initial orbit determination *    
!  *                                                               *    
!  *****************************************************************    
!                                                                       
! INPUT:    W         -  Triplet weights                                
!           P         -  Triplet priority                               
!           N         -  Number of triplets                             
!                                                                       
! OUTPUT:   ORD       -  Sorting order                                  
!                                                                       
      SUBROUTINE s3srt(w,p,n,ord) 
      IMPLICIT NONE 
                                                                        
      INTEGER n 
      INTEGER p(n),ord(n) 
      DOUBLE PRECISION w(n) 
                                                                        
      INTEGER i,it,itmp,k1,k2 
      LOGICAL ok,rev 
                                                                        
      DO 1 i=1,n 
      ord(i)=i 
    1 END DO 
                                                                        
      DO 3 it=1,n+2 
      ok=.true. 
      DO 2 i=1,n-1 
      k1=ord(i) 
      k2=ord(i+1) 
      IF(p(k1).GT.p(k2)) THEN 
          rev=.false. 
      ELSEIF(p(k1).EQ.p(k2)) THEN 
          rev=(w(k1).GT.w(k2)) 
      ELSE 
          rev=.true. 
      END IF 
      IF(rev) THEN 
          itmp=ord(i) 
          ord(i)=ord(i+1) 
          ord(i+1)=itmp 
          ok=.false. 
      END IF 
    2 END DO 
      IF(ok) RETURN 
    3 END DO 
                                                                        
      STOP '**** s3srt: internal error (01) ****' 
                                                                        
      END                                           

! Copyright (C) 2000 by Mario Carpino (carpino@brera.mi.astro.it)       
! Version: May 31, 2000                                                 
! --------------------------------------------------------------------- 
!                                                                       
!  *****************************************************************    
!  *                                                               *    
!  *                         A D D S 3 P                           *    
!  *                                                               *    
!  *      Add a triplet to the list of observations triplets       *    
!  *                                                               *    
!  *****************************************************************    
!                                                                       
! INPUT:    K1,K2,K3  -  Triplet                                        
!           PRI       -  Priority                                       
!           TDT       -  Time of observations                           
!           N         -  Number of observations                         
!           DTW       -  Time interval for weight computation           
!                                                                       
      SUBROUTINE adds3p(k1,k2,k3,pri,tdt,n,dtw) 
      IMPLICIT NONE 
                                                                        
      INTEGER k1,k2,k3,pri,n 
      DOUBLE PRECISION tdt(n),dtw 
                                                                        
      INCLUDE 'pars3m.h90' 
      INCLUDE 'sel3m.h90' 
                                                                        
      INTEGER i,iplw,iins 
      DOUBLE PRECISION w,wmax 
                                                                        
      DOUBLE PRECISION s3dtw 
      EXTERNAL s3dtw 
                                                                        
! Check whether the triplet is already selected, and possibly determine 
! the triplet with largest weight (among the triplets having lower or   
! equal priority)                                                       
      iplw=0 
      DO 1 i=1,n3s 
      IF(l3s(1,i).NE.k1) GOTO 2 
      IF(l3s(2,i).NE.k2) GOTO 2 
      IF(l3s(3,i).NE.k3) GOTO 2 
      RETURN 
    2 CONTINUE 
      IF(pri3(i).GT.pri) GOTO 1 
      IF(iplw.EQ.0) THEN 
          iplw=i 
          wmax=w3s(i) 
      ELSEIF(w3s(i).GT.wmax) THEN 
          iplw=i 
          wmax=w3s(i) 
      END IF 
    1 END DO 
                                                                        
! Weight                                                                
      w=s3dtw(tdt(k2)-tdt(k1),dtw)+s3dtw(tdt(k3)-tdt(k2),dtw) 
                                                                        
      iins=0 
      IF(n3s.LT.n3sm) THEN 
          n3s=n3s+1 
          iins=n3s 
      ELSE 
          IF(iplw.EQ.0) RETURN 
          iins=iplw 
      END IF 
      IF(iins.EQ.0) RETURN 
                                                                        
      l3s(1,iins)=k1 
      l3s(2,iins)=k2 
      l3s(3,iins)=k3 
      pri3(iins)=pri 
      w3s(iins)=w 
                                                                        
      END                                           

! Copyright (C) 2000 by Mario Carpino (carpino@brera.mi.astro.it)       
! Version: June 9, 2000                                                 
! --------------------------------------------------------------------- 
!                                                                       
!  *****************************************************************    
!  *                                                               *    
!  *                         I O D S D T                           *    
!  *                                                               *    
!  *    Select time interval (observations) for RMS computation    *    
!  *                                                               *    
!  *****************************************************************    
!                                                                       
! INPUT:    SELIPT    -  Selected triplet                               
!           TOBS      -  Time of observations (MJD, TDT)                
!           NOBS      -  Number of observations                         
!           EXTF      -  Interval extension factor                      
!           DTMAX     -  Max interval lenght (d)                        
!                                                                       
! OUTPUT:   IR        -  First and last observations to be used for     
!                        RMS computation                                
!                                                                       
! The routine selects the starting and ending observations for          
! computing a residual RMS suitable for checking the accuracy of a      
! preliminary orbit: for this purpose, the selected observations must   
! be approximately centered on the triplet of observations used         
! for orbit determination, and span a time interval not too larger      
! than the timespan covered by the triplet                              
!                                                                       
      SUBROUTINE iodsdt(selipt,tobs,nobs,extf,dtmax,ir) 
      IMPLICIT NONE 
                                                                        
      INTEGER selipt(3),nobs,ir(2) 
      DOUBLE PRECISION tobs(nobs),extf,dtmax 
                                                                        
      INTEGER is1,is2,i 
      DOUBLE PRECISION dt 
                                                                        
      is1=selipt(1) 
      is2=selipt(3) 
                                                                        
! Acceptable time interval before the first observation and after the   
! last observation used for initial orbit determination.                
! The time interval must not exceed EXTF times the timespan covered     
! by the triplet (if EXTF<0, this limitation is not applied); moreover, 
! it must not be longer than DTMAX days (if DTMAX<0, this limitation    
! is not applied)                                                       
      IF(extf.GE.0) THEN 
          dt=(tobs(is2)-tobs(is1))*extf 
      ELSE 
          dt=10*(tobs(nobs)-tobs(1)) 
      END IF 
      IF(dtmax.GE.0) dt=MAX(dt,dtmax) 
                                                                        
      DO 1 i=is1,1,-1 
      IF(tobs(is1)-tobs(i).GT.dt) GOTO 2 
      ir(1)=i 
    1 END DO 
    2 CONTINUE 
                                                                        
      DO 3 i=is2,nobs 
      IF(tobs(i)-tobs(is2).GT.dt) GOTO 4 
      ir(2)=i 
    3 END DO 
    4 CONTINUE 
                                                                        
      END                                           

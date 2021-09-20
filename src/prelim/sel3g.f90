! Copyright (C) 1997-1999 by Mario Carpino (carpino@brera.mi.astro.it)  
! Version: February 10, 1999                                            
! --------------------------------------------------------------------- 
!                                                                       
!  *****************************************************************    
!  *                                                               *    
!  *                          S E L 3 G                            *    
!  *                                                               *    
!  *       Selection of 3 observations for Gauss' method           *    
!  *                                                               *    
!  *****************************************************************    
!                                                                       
! INPUT:    TDT       -  Time (MJD, TDT)                                
!           IOBS      -  Observation type (1000's=RA,DEC; 2000's=RADAR) 
!           RMSA      -  A-priori RMS of right ascension (rad)          
!           RMSD      -  A-priori RMS of declination (rad)              
!           SEL       -  Selection indicator                            
!                            0=don't use                                
!                            1=use for fit                              
!                            2=use for fit & Gauss method               
!           N         -  Number of observations                         
!                                                                       
! OUTPUT:   SELIPT    -  Pointer to selected observations               
!                                                                       
! WARNING: the routine assumes the observations to be sorted in order   
!          of increasing time                                           
!                                                                       
      SUBROUTINE sel3g(tdt,iobs,rmsa,rmsd,sel,n,selipt)
      USE fund_const 
      IMPLICIT NONE 
                                                                        
      INTEGER n,sel(n),iobs(n),selipt(3) 
      DOUBLE PRECISION tdt(n),rmsa(n),rmsd(n) 
                                                                        
! TUNING PARAMETERS                                                     
! Minimum, advisable and maximum distance between observations (d)      
      DOUBLE PRECISION dt1,dt2,dt3 
      PARAMETER (dt1=0.03d0) 
      PARAMETER (dt2=10.d0) 
      PARAMETER (dt3=150.d0) 
! Minimum required precision (arcsec)                                   
      DOUBLE PRECISION minrms 
      PARAMETER (minrms=5.d0) 
                                                                        
      DOUBLE PRECISION rms1,ts1,ts2,w1,w2,wt,wmin 
      INTEGER k,k1,k2,k3,nsel 
      LOGICAL nosol 
                                                                        
! Check first whether the selection indicators contains already a list  
! of 3 points to be used                                                
      nsel=0 
      DO 4 k=1,n 
      IF(iobs(k)/1000.NE.1) GOTO 4 
      IF(sel(k).GE.2) THEN 
          nsel=nsel+1 
          selipt(nsel)=k 
          IF(nsel.GE.3) GOTO 5 
      END IF 
    4 END DO 
    5 CONTINUE 
      IF(nsel.GE.3) RETURN 
                                                                        
! If no indication is contained in the selection indicators,            
! use an automated selection method                                     
      DO 6 k=1,n 
      IF(iobs(k)/1000.NE.1) GOTO 6 
      sel(k)=MIN(sel(k),1) 
    6 END DO 
                                                                        
      rms1=minrms*radsec 
                                                                        
      selipt(1)=0 
      selipt(2)=0 
      selipt(3)=0 
      nosol=.true. 
                                                                        
! First trial observation                                               
      DO 1 k1=1,n-2 
      IF(iobs(k1)/1000.NE.1) GOTO 1 
      IF(sel(k1).LT.1) GOTO 1 
      IF(rmsa(k1).GT.rms1) GOTO 1 
      IF(rmsd(k1).GT.rms1) GOTO 1 
                                                                        
! Second trial observation                                              
      DO 2 k2=k1+1,n-1 
      IF(iobs(k2)/1000.NE.1) GOTO 2 
      IF(sel(k2).LT.1) GOTO 2 
      IF(rmsa(k2).GT.rms1) GOTO 2 
      IF(rmsd(k2).GT.rms1) GOTO 2 
      ts1=tdt(k2)-tdt(k1) 
      IF(ts1.LT.0.d0) STOP '**** sel3g: internal error (01) ****' 
      IF(ts1.LT.dt1) GOTO 2 
      IF(ts1.GT.dt3) GOTO 1 
      w1=ABS(ts1-dt2)/dt2 
                                                                        
! Third trial observation                                               
      DO 3 k3=k2+1,n 
      IF(iobs(k3)/1000.NE.1) GOTO 3 
      IF(sel(k3).LT.1) GOTO 3 
      IF(rmsa(k3).GT.rms1) GOTO 3 
      IF(rmsd(k3).GT.rms1) GOTO 3 
      ts2=tdt(k3)-tdt(k2) 
      IF(ts2.LT.0.d0) STOP '**** sel3g: internal error (02) ****' 
      IF(ts2.LT.dt1) GOTO 3 
      IF(ts2.GT.dt3) GOTO 2 
      w2=ABS(ts2-dt2)/dt2 
      wt=w1+w2 
      IF(nosol) wmin=wt+100 
      IF(wt.LT.wmin) THEN 
          wmin=wt 
          selipt(1)=k1 
          selipt(2)=k2 
          selipt(3)=k3 
          nosol=.false. 
      END IF 
                                                                        
    3 END DO 
    2 END DO 
    1 END DO 
                                                                        
      IF(nosol) THEN 
!          WRITE(*,100)                                                 
!          STOP '**** sel3g: abnormal end ****'                         
           RETURN 
      END IF 
  100 FORMAT(' ERROR: automatic selection of observations for Gauss'/   &
     &       '        method failed: please define manually which'/     &
     &       '        observations have to be used')                    
                                                                        
      sel(selipt(1))=2 
      sel(selipt(2))=2 
      sel(selipt(3))=2 
      write(*,*)'sel3g: selected ',selipt(1), selipt(2),selipt(3) 
      END                                           

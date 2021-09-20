! Copyright (C) 2000 by Mario Carpino (carpino@brera.mi.astro.it)       
! Version: December 11, 2000                                            
! --------------------------------------------------------------------- 
!                                                                       
!  *****************************************************************    
!  *                                                               *    
!  *                         S E L 3 _ M C                         *    
!  *                                                               *    
!  *  Selection of 3 observations for initial orbit determination  *    
!  *        (STEP 1: computation of the list of triplets)          *    
!  *                                                               *    
!  *****************************************************************    
!                                                                       
! INPUT:    TDT       -  Time (MJD, TDT)                                
!           TYPE      -  Observation type                               
!           RMSA      -  A-priori RMS of right ascension (rad)          
!           RMSD      -  A-priori RMS of declination (rad)              
!           SEL       -  Selection indicator                            
!                            0=don't use                                
!                            1=use for fit                              
!                            2=use for fit & Gauss method               
!           N         -  Number of observations                         
!           N3MAX     -  Max number of triplets                         
!                                                                       
! WARNING: the routine assumes the observations to be sorted in order   
!          of increasing time                                           
!                                                                       
      SUBROUTINE sel3_mc(tdt,type,rmsa,rmsd,sel,n,n3max) 
      USE fund_const
      IMPLICIT NONE
      INTEGER n
      CHARACTER*(1) type(n) 
      INTEGER sel(n),n3max 
      DOUBLE PRECISION tdt(n),rmsa(n),rmsd(n) 
                                                                        
! TUNING PARAMETERS                                                     
! Minimum, advisable and maximum distance between observations (d)      
      DOUBLE PRECISION dt1,dt2,dt3 
      PARAMETER (dt1=0.03d0) 
      PARAMETER (dt2=20.d0) 
      PARAMETER (dt3=150.d0) 
                                                                        
      INTEGER k,i1,i2,i3 
                                                                        
      INCLUDE 'parobx.h90' 
      INCLUDE 'pars3m.h90' 
                                                                        
! Common blocks to be initialized:                                      
      INCLUDE 'sel3m.h90' 
                                                                        
      INTEGER lsmo(nobx),nsmo 
      DOUBLE PRECISION maxrms 
      LOGICAL first,found,fail1,fail,only2,addnew 
      DATA first/.true./ 
      SAVE first,maxrms,only2 
                                                                        
      IF(first) THEN 
          maxrms=4.D0 
          fail=.false. 
          CALL rdnrea('init_orbdet.sel3.','max_rms',maxrms,.false.,     &
     &                found,fail1,fail)                                 
          only2=.true. 
          CALL rdnlog('init_orbdet.sel3.','only2',only2,.false.,        &
     &            found,fail1,fail)                                     
          IF(fail) STOP '**** sel3mc: abnormal end ****' 
          maxrms=maxrms*radsec 
          first=.false. 
      END IF 
                                                                        
      IF(n.GT.nobx) STOP '**** sel3mc: internal error (01) ****' 
      IF(n3max.GT.n3smax) STOP '**** sel3mc: internal error (02) ****' 
      n3s=0 
      n3sm=n3max 
                                                                        
! If some observations are marked in the .rwo file by a selection       
! indicator equal to 2 (or larger), they are used preferentially,       
! without even checking their a-priori RMS                              
      nsmo=0 
      DO 1 k=1,n 
      IF(type(k).NE.'O') GOTO 1 
      IF(sel(k).GE.2) THEN 
          nsmo=nsmo+1 
          lsmo(nsmo)=k 
      END IF 
    1 END DO 
      IF(nsmo.GE.3) THEN 
          DO 2 i1=1,nsmo 
          DO 3 i2=i1+1,nsmo 
          DO 4 i3=i2+1,nsmo 
          CALL adds3p(lsmo(i1),lsmo(i2),lsmo(i3),2,tdt,n,dt2) 
    4     CONTINUE 
    3     CONTINUE 
    2     CONTINUE 
      END IF 
                                                                        
! Must I add new triplets?                                              
      addnew=.true. 
! If there is already at least one triplet, and the user required to use
! only observations marked with selection indicator=2, then no other    
! triplet is added                                                      
      IF(only2.AND.(n3s.GT.0)) addnew=.false. 
! If the number of triplets already defined exceeds the max value       
! allowed, no other triplet is added                                    
      IF(n3s.GE.n3max) addnew=.false. 
                                                                        
! Add new triplets                                                      
      IF(addnew) THEN 
          DO 5 i1=1,n 
          IF(type(i1).NE.'O') GOTO 5 
          IF(sel(i1).LT.1)       GOTO 5 
          IF(rmsa(i1).GT.maxrms) GOTO 5 
          IF(rmsd(i1).GT.maxrms) GOTO 5 
          DO 6 i2=i1+1,n 
          IF(type(i2).NE.'O') GOTO 6 
          IF(sel(i2).LT.1)       GOTO 6 
          IF(rmsa(i2).GT.maxrms) GOTO 6 
          IF(rmsd(i2).GT.maxrms) GOTO 6 
          DO 7 i3=i2+1,n 
          IF(type(i3).ne.'O') GOTO 7 
          IF(sel(i3).LT.1)       GOTO 7 
          IF(rmsa(i3).GT.maxrms) GOTO 7 
          IF(rmsd(i3).GT.maxrms) GOTO 7 
          CALL adds3p(i1,i2,i3,1,tdt,n,dt2) 
    7     CONTINUE 
    6     CONTINUE 
    5     CONTINUE 
      END IF 
                                                                        
! Sorting of triplets of observations                                   
      CALL s3srt(w3s,pri3,n3s,ords3) 
                                                                        
      iics3m=36 
      ipt3=0 
                                                                        
      END SUBROUTINE sel3_mc                                          

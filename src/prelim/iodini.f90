! Copyright (C) 1998-2000 by Mario Carpino (carpino@brera.mi.astro.it)  
! Version: June 16, 2000                                                
! --------------------------------------------------------------------- 
!                                                                       
!  *****************************************************************    
!  *                                                               *    
!  *                         I O D I N I                           *    
!  *                                                               *    
!  *       Input of options for initial orbit determination        *    
!  *                    (Gauss, Vaisala, etc.)                     *    
!  *                                                               *    
!  *****************************************************************    
!                                                                       
      SUBROUTINE iodini 
      USE fund_const
      IMPLICIT NONE 
                                                                        
      INCLUDE 'pariod.h90' 
      INCLUDE 'pars3m.h90' 
                                                                        
! Common blocks to be initialized:                                      
      INCLUDE 'comiod.h90' 
                                                                        
      INTEGER i,lt 
      LOGICAL found,fail1,fail,inigau 
      CHARACTER*50 tmp 
                                                                        
      INTEGER lench 
      EXTERNAL lench 
                                                                        
      fail=.false. 
      inigau=.false. 
                                                                        
      iodvrb=1 
      CALL rdnint('init_orbdet.','verbose',iodvrb,.false.,              &
     &            found,fail1,fail)                                     
                                                                        
      iodmul=.false. 
      CALL rdnlog('init_orbdet.','multi_line',iodmul,.false.,           &
     &            found,fail1,fail)                                     
                                                                        
      iodntr=MIN(10,n3smax) 
      CALL rdnint('init_orbdet.','n_triplets',iodntr,.false.,           &
     &            found,fail1,fail)                                     
      IF(iodntr.GT.n3smax) THEN 
          iodntr=n3smax 
          WRITE(*,200) iodntr 
      END IF 
  200 FORMAT(' WARNING: the value of "init_orbdet.n_triplets" has been',&
     &       ' decreased to',I4/                                        &
     &       10X,'in order not to exceed the value of PARAMETER n3smax')
                                                                        
      iodexp=-1.D0 
      CALL rdnrea('init_orbdet.check_dt.','expand',iodexp,.false.,      &
     &            found,fail1,fail)                                     
                                                                        
      ioddtm=30.D0 
      CALL rdnrea('init_orbdet.check_dt.','maxdtexp',ioddtm,.false.,    &
     &            found,fail1,fail)                                     
                                                                        
      iodrok=2.D0 
      CALL rdnrea('init_orbdet.rms_val.','ok',iodrok,.false.,           &
     &            found,fail1,fail)                                     
!      iodrok=iodrok*radsec 
                                                                        
      iodrmx=1.D2 
      CALL rdnrea('init_orbdet.rms_val.','max',iodrmx,.false.,          &
     &            found,fail1,fail)                                     
!      iodrmx=iodrmx*radsec 
                                                                        
      iodnit=10 
      CALL rdnint('init_orbdet.noise.','ntrials',iodnit,.false.,        &
     &            found,fail1,fail)                                     
      iodnit=MAX(0,iodnit) 
                                                                        
      iodksi=1.D0 
      CALL rdnrea('init_orbdet.noise.','ksigma',iodksi,.false.,         &
     &            found,fail1,fail)                                     
                                                                        
! List of methods to be used                                            
      CALL rdmcha('init_orbdet.','methods',iodmen,iodnm,iodnmx,         &
     &            'iodnmx',.false.,found,fail1,fail)                    
      IF(fail) STOP '**** iodini: abnormal end ****' 
                                                                        
      IF(found) THEN 
          DO 1 i=1,iodnm 
          tmp=iodmen(i) 
          CALL locase(tmp) 
          IF(tmp.EQ.'1' .OR. tmp.EQ.'gauss') THEN 
              iodmet(i)=1 
              iodmen(i)='Gauss' 
              inigau=.true. 
          ELSEIF(tmp.EQ.'2' .OR. tmp.EQ.'vaisala') THEN 
              iodmet(i)=2 
              iodmen(i)='Vaisala' 
          ELSE 
              lt=lench(tmp) 
              WRITE(*,110) tmp(1:lt) 
              fail=.true. 
          END IF 
    1     CONTINUE 
      ELSE 
          IF(iodnmx.LT.2) STOP '**** iodini: iodnmx < 2 ****' 
          iodnm=2 
          iodmet(1)=1 
          iodmen(1)='Gauss' 
          iodmet(2)=2 
          iodmen(2)='Vaisala' 
          inigau=.true. 
      END IF 
  110 FORMAT(' ERROR: unsupported initial orbit determination method'/  &
     &       '        ("',A,'", keyword "init_orbdet.methods")')        
                                                                        
      IF(inigau) CALL gauini(fail) 
      IF(fail) STOP '**** iodini: abnormal end ****' 
                                                                        
      iiciod=36 
                                                                        
      END SUBROUTINE iodini                                          

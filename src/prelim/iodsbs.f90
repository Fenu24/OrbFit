! Copyright (C) 2000 by Mario Carpino (carpino@brera.mi.astro.it)       
! Version: June 9, 2000                                                 
! --------------------------------------------------------------------- 
!                                                                       
!  *****************************************************************    
!  *                                                               *    
!  *                         I O D S B S                           *    
!  *                                                               *    
!  *  Store best solution of initial orbit determination process   *    
!  *                                                               *    
!  *****************************************************************    
!                                                                       
! I/O:      ELEMB     -  Orbital elements of best solution              
!           TYPEB     -  Type of orbital elements of best solution      
!           TB        -  Epoch of orbital elements of best solution     
!           RMSB      -  RMS error of best solution                     
!           METHB     -  Method used for finding best solution          
!           SELB      -  Selected triplet of best solution              
!           IRB       -  Select check interval for best solution        
!           ELEM1     -  Orbital elements of new solution               
!           TYPE1     -  Type of orbital elements of new solution       
!           T1        -  Epoch of orbital elements of new solution      
!           RMS1      -  RMS error of new solution                      
!           METH1     -  Method used for finding new solution           
!           SEL1      -  Selected triplet of new solution               
!           IR1       -  Select check interval for new solution         
!           EXIST     -  Existence flag of best solution                
!                                                                       
      SUBROUTINE iodsbs(elemb,typeb,tb,rmsb,methb,selb,irb,             &
     &                  elem1,type1,t1,rms1,meth1,sel1,ir1,exist)       
      IMPLICIT NONE 
                                                                        
      INTEGER selb(3),sel1(3),irb(2),ir1(2) 
      DOUBLE PRECISION elemb(6),elem1(6),tb,t1,rmsb,rms1 
      CHARACTER*(*) typeb,type1,methb,meth1 
      LOGICAL exist 
                                                                        
      INTEGER i 
      LOGICAL keep 
                                                                        
      IF(exist) THEN 
          keep=(rms1.LT.rmsb) 
      ELSE 
          keep=.true. 
      END IF 
                                                                        
      IF(keep) THEN 
          DO 1 i=1,6 
          elemb(i)=elem1(i) 
    1     CONTINUE 
          typeb=type1 
          tb=t1 
          rmsb=rms1 
          methb=meth1 
          exist=.true. 
          selb(1)=sel1(1) 
          selb(2)=sel1(2) 
          selb(3)=sel1(3) 
          irb(1)=ir1(1) 
          irb(2)=ir1(2) 
      END IF 
                                                                        
      END                                           

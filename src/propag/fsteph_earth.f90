! Copyright 2002 Orbfit Consortium                                      
! Version 9 May 2002 A. Milani                                          
! ====================================================                  
! Earth  Ephemerides Generation for FITOBS                              
! This is a hacked version of FSTPRO                                    
! First integrate from t0 back to tr, saving at interval step.          
! Then write this data in chronological order                           
! Finally integrate from to forward to tf, writing after each step      
! inputs:                                                               
!     dir - location for files                                          
!     defele - logical for existence of elements                        
!     ok - success flag                                                 
!     t - reference epoch time                                          
!     tr - initial time for ephemerides                                 
!     tf - final time                                                   
!     step - stepsize                                                   
!     numsav - an upper bound for the number of steps necessary before t
!     ephefl - logical flag for output of ephemrides file               
!     eltype - coordinate type for ephemerides (output) elements 
!              currently only CAR are allowed, otherwise no output
! output: NONE                                                          
!                                                                       
! REMARK: IF moidfl=.false. and ephefl=.false. you will still           
!     get a close approach file.                                        
! ====================================================                  
      SUBROUTINE fsteph_earth(dir,ok,t0,tr,tf,step,ephefl,eltype) 
        USE orbit_elements
        IMPLICIT NONE 
! =================INPUT=========================================       
! place, output element type                                            
      CHARACTER*(*) dir 
      CHARACTER*(3) eltype 
! logical to controlo actual writing                                    
      LOGICAL ephefl 
! epoch times: center,beginning, end, step                              
      DOUBLE PRECISION t0,tr,tf,step 
! ================OUTPUT=================================               
! ok                                                                    
      LOGICAL ok 
! ================END INTERFACE==========================               
! earth name                                                            
      CHARACTER*5 name 
! loop indexes,number of records                                        
      INTEGER jst,j,nst 
! times                                                                 
      DOUBLE PRECISION t1,t2,tsav(100000),tarr(100000) 
! coordinates of earth                                                  
      DOUBLE PRECISION elem1(6) 
! output                                                                
      INTEGER unit 
      INTEGER ln,lnm,lnnam 
      CHARACTER*(60) file,filem 
! ===================================================================== 
! name=earth                                                            
!     WRITE(*,*)'fsteph-earth: t0,tr,tf', t0,tr,tf                      
      name='EARTH' 
      IF(eltype.ne.'CAR')RETURN 
      IF(.not.ephefl)RETURN 
! ===================================================================== 
! check availability of required data                                   
      ok=.true. 
!      IF (tr .ge. tf)THEN 
!         WRITE(*,*) 'Sorry: Initial time must be before final time.' 
!         ok = .false. 
!      ENDIF 
!      IF(.not.ok)RETURN 
! check availability of JPL ephemrides                                  
      call chetim(tr,tf,ok) 
      if(.not.ok)return 
!========= OPEN FILES =================================                 
      if(ephefl)then 
! open ephemerides ouput file                                           
         call filnam(dir,name,'ele',file,ln) 
         call filopn(unit,file(1:ln),'unknown') 
         call wro1lh(unit,'ECLM','J2000',eltype) 
      endif 
!========= PROPAGATE BACK ========================                      
      IF(tr .le. t0)THEN 
         IF(tf.le.t0)THEN 
! propagate backward to tf                                              
            t1=tf 
         ELSE 
! start from the central time                                           
            t1=t0 
         ENDIF 
         tsav(1)=t1 
! loop on times                                                         
         DO jst=1,100000 
           t2=t1-jst*step 
           IF(t2.gt.tr)THEN 
              tsav(jst+1)=t2 
           ELSE 
              tsav(jst+1)=t2 
              nst=jst+1 
              GOTO 5 
           ENDIF 
         ENDDO 
    5    CONTINUE 
! reverse array                                                         
         DO j=1,nst 
            tarr(j)=tsav(nst-j+1) 
         ENDDO 
      ENDIF 
!========= PROPAGATE FORWARD ========================                   
      IF(tf.gt.t0)THEN 
         IF(tr.gt.t0)THEN 
! propagate forward to tf                                               
            t1=tr 
         ELSE 
! start from the central time                                           
            t1=t0 
         ENDIF 
         tarr(1)=t1
! loop on times                                                         
         DO jst=1,100000 
           t2=t1+jst*step 
           IF(t2.lt.tf)THEN 
              tarr(jst+nst)=t2 
           ELSE 
              tarr(jst+nst)=t2 
              nst=jst+nst 
              GOTO 6 
           ENDIF 
         ENDDO 
    6    CONTINUE 
      ENDIF 
! Earth coordinates at time set above                                   
      DO jst=1,nst 
         t2=tarr(jst) 
!         write(*,*)'jst,t2',jst,t2                                     
         CALL earcar(t2,elem1,1) 
! magnitudes H and G are dummy                                          
         IF(ephefl) CALL wro1lr(unit,name(1:5),elem1,eltype,            &
     &        t2,-1.d4,1.0d0)                                           
      ENDDO 
! ===================================================================== 
      if(ephefl) call filclo(unit,' ') 
      if(ephefl) write(*,*)'Ephemeris file generated:',file(1:ln) 
      RETURN 
      END SUBROUTINE fsteph_earth                                          


! ===========================================                           
! ORBSELECT                                                             
! extract from file .cat (in OEF OrbFit computer readable format)       
! by semimajor axis, eccentricity, etc.                                 
PROGRAM orbselect 
  USE fund_const
  IMPLICIT NONE 
! controls                                                              
  INTEGER nmax 
! PARAMETER (nmax=50000)                                           
  DOUBLE PRECISION amin,amax,emin,emax,sqei,qmin 
  INTEGER i 
  INTEGER numae 
  DOUBLE PRECISION ecc 
! name of asteroid being computed: long. short                          
  CHARACTER*19 namid 
  CHARACTER*9 name 
  INTEGER ll 
! definition flags                                                      
  LOGICAL cov0,defnor 
! orbital elements, epoch                                               
  DOUBLE PRECISION elem(6),telem,tel0 
! magnitude, opposition effect, mass, matrices                          
  DOUBLE PRECISION hm,gma0,mass,g0(6,6),c0(6,6) 
! output of rdorb                                                       
  CHARACTER*5 epoch 
  CHARACTER*4 rsys 
  CHARACTER*3 coox 
! end of file                                                           
  LOGICAL eof 
  INTEGER nrecord 
! selection criteria                                                    
  WRITE(*,*)' nmax record to read?' 
  READ(*,*)nmax 
  WRITE(*,*)' semimajor axis min,max?' 
  READ(*,*)amin,amax 
  WRITE(*,*)' ecc/incl min,max?' 
  READ(*,*)emin,emax 
  WRITE(*,*)' minimum q?' 
  READ(*,*)qmin 
! open input file                                                       
  CALL libini 
  CALL oporbf('input.cat', -1) 
! here all the asteroids with their stepsize                            
  OPEN(3,file='orb.sel',status='unknown') 
! here only the ones with the right stepsize                            
!      OPEN(2,file='orb.step',status='unknown')                         
! counters: select by a/e, select by step                               
  numae=0 
! main loop                                                             
  DO 1 i=1,nmax 
! read one record                                                       
     CALL rdorb(namid,elem,coox,telem,g0,cov0,c0,defnor,            &
     &           hm,gma0,mass,rsys,epoch,nrecord,eof)                   
     IF(eof) GOTO 2 
     IF(coox.ne.'KEP')THEN 
        write(*,*)' wrong coordinates in input.cat' 
        stop 
     ENDIF
! handle funny identification names                                     
! convert name in case it is of the form nam0=namp                      
     ll=index(namid,'=')-1 
     IF(ll.lt.0)THEN 
        name=namid(1:9) 
     ELSE 
        name=namid(1:ll) 
     ENDIF
! initial time, hoping it is always the same                            
     IF(i.eq.1)THEN 
        tel0=telem 
        WRITE(*,*)tel0 
     ELSE 
        IF(telem.ne.tel0)THEN 
           WRITE(*,*)' discordant times ', telem, 'first was ',tel0 
        ENDIF
     ENDIF
! square root of eccentricity squared plus inclination squared          
     sqei=sqrt(elem(2)**2+(sin(elem(3)))**2) 
! select                                                                
     IF(elem(1).ge.amin.and.elem(1).le.amax.and.                     &
     &       sqei.ge.emin.and.sqei.le.emax.and.                         &
     &       elem(1)*(1.d0-elem(2)).gt.qmin)THEN                        
        numae=numae+1 
        WRITE(3,*)name 
     ENDIF
1 ENDDO
! file too long                                                         
  WRITE(*,*)' file not complete after ',nmax,' records' 
2 CLOSE(3) 
! statistical summary                                                   
  WRITE(*,*)' selected: ',numae,' by a/e/i ' 
  WRITE(*,*)' from ', i-1, ' records' 
END PROGRAM orbselect

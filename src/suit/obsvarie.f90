! LIBRARY OBSVARIE contains
! iaucod        computes IAU official designations from MPC-style packed
!               with alphamumeric numbers
! iaucod2       computes IAU official designations from MPC-style packed
!               with alphamumeric numbers, padding with 'w'
! mpcpds	computes MPC-style packed designation from official IAU code   
!               also with alphanumeric codes
! nights        number of nights of observations; 
!                 requires the observatiosn to be sorted by time
! quality       quality code
! obssta        1-line report on observations available
! Copyright (C) 1998,2003 by OrbFit Consortium                               
! Original version: December 15, 1997 Steven Chesley
! Revised by Genny in Jun 20, 2001 to handle designations as K00Sa3P 
! f90 version September 2003 A. Milani
! expanded to include obssta, quality October 2004 A. Milani   
! --------------------------------------------------------------------- 
!  *****************************************************************    
!  *                                                               *    
!  *                         I A U C O D N E W                     *    
!  *                                                               *    
!  * Computes official IAU code from MPC new packed designation    *    
!  *                                                               *    
!  *****************************************************************    
!                                                                       
! INPUT:    MPCCOD    -  MPC-style code as used in observation archives 
!                        7 characters for number, 9 characters for      
!                        packed provisional designation, 10 for
!                        comets and satellites                 
!                                                                       
! OUTPUT:   IAUDES    -  IAU code                                       
!           ERROR     -  Error flag (cannot understand input code)      
!                                                                       
      SUBROUTINE iaucodnew(mpccod,iaudes,error)
      USE name_rules 
      IMPLICIT NONE 
      CHARACTER*(name_len), INTENT(OUT) :: iaudes
      CHARACTER*(*), INTENT(IN) :: mpccod 
      LOGICAL, INTENT(OUT) :: error 
! ---------END INTERFACE--------------------
      INTEGER ln,i,temp 
      CHARACTER*2 head 
      CHARACTER*3 tail 
      INTEGER lench, number 
      LOGICAL isnum 
      EXTERNAL lench,isnum 
      CHARACTER numfield*7, desfield*9 
! ================================================
      error=.false. 
      iaudes=' '                                                     
      ln=lench(mpccod) 
      IF(ln.LE.0) GOTO 10 
! check if numbered                                                                        
      numfield=mpccod(1:7) 
      CALL rmsp(numfield,ln)                                
      IF(ln.ne.0) THEN
! Numbered asteroids  
         IF(isnum(numfield(1:1)))THEN
! standard number: remove zeros
            iaudes=numfield
! remove leading zeros
            DO i=1,ln-1 
               IF (iaudes(i:i).eq.'0') THEN 
                  iaudes(i:i)=' ' 
               ELSE 
                  EXIT 
               ENDIF
            ENDDO
            RETURN ! successfull parsing of numbered asteroid
         ELSE
! alphanumeric asteroid number
            WRITE(*,*)' iaucodnew: alphamnumeric asteroid number? ', mpccod
            STOP
         ENDIF         
      ENDIF
                                                                        
! Unnumbered asteroids                                                  
      desfield=mpccod(8:16) 
                                                                        
      if(desfield(3:3).eq.'S')then 
! Survey Asteroid                                                       
         iaudes=desfield(6:9)//desfield(1:1)//'-'//desfield(2:2) 
                                                                        
      elseif (desfield(1:1).eq.'I' .or.                                 &
     &        desfield(1:1).eq.'J' .or.                                 &
     &        desfield(1:1).eq.'K')then                                 
! Temporary designation       
         tail=' '                                          
         if(desfield(5:8).eq.'0000')then 
!           1999AA = J99A00A  
         elseif(desfield(5:7).eq.'000')then 
!           1999AA1 = J99A01A 
            tail=desfield(8:8)
          elseif(desfield(5:6).eq.'00')then 
!           1999AA12 = J99A12A
            tail=desfield(7:8)
         elseif(desfield(5:5).eq.'0')then 
!           1999AA123 = J99A0123A
            tail=desfield(6:8)
         else
!           1999AA1234 = J99A1234A
            IF(name_len.lt.10)THEN !PROBLEM: 1999AA1234 cannot be accomodated in 9 characters
               WRITE(*,*)' iaucodnew: 10 character IAU designation? ', mpccod
               STOP
            ELSE
               tail=desfield(5:8)
            ENDIF
         endif 
         temp=ichar(desfield(1:1))-55 
         write(head,103) temp 
  103    format(I2) 
         iaudes=head//desfield(2:3)//desfield(4:4)//desfield(7:7)//tail 
      else 
! Unknown type                                                          
         write(*,*)'cannot understand MPC designation: ',mpccod 
         goto 10 
      endif 
      return 
                                                                        
   10 CONTINUE 
      iaudes=mpccod 
      error=.true. 
                                                                        
    END SUBROUTINE iaucodnew
!                                                                       
!  *****************************************************************    
!  *                                                               *    
!  *                         I A U C O D                           *    
!  *                                                               *    
!  * Computes official IAU code from MPC-style packed designation  *    
!  *                                                               *    
!  *****************************************************************    
!                                                                       
! INPUT:    MPCCOD    -  MPC-style code as used in observation archives 
!                        5 characters for number, 7 characters for      
!                        packed provisional designation                 
!                                                                       
! OUTPUT:   IAUDES    -  IAU code                                       
!           ERROR     -  Error flag (cannot understand input code)      
!                                                                       
      SUBROUTINE iaucod(mpccod,iaudes,error) 
      IMPLICIT NONE 
                                                                        
      CHARACTER*(*) iaudes,mpccod 
      LOGICAL error 
                                                                        
      INTEGER ln,i,temp,permdes 
      CHARACTER*2 head 
      CHARACTER*3 tail 
                                                                        
      INTEGER lench 
      LOGICAL isnum
      INTEGER numdes(4)
      EXTERNAL lench,isnum 
                                                                        
      CHARACTER numfield*5, desfield*12 
                                                                        
      error=.false. 
      iaudes=' ' 
                                                                        
      ln=lench(mpccod) 
      IF(ln.LE.0) GOTO 10 
! check if numbered 
      numfield=mpccod(1:5) 
      CALL rmsp(numfield,ln)                                
      IF(ln.eq.5) THEN
! Numbered asteroids  
         IF(isnum(numfield(1:1)))THEN
! standard number
            iaudes=numfield
! remove leading zeros
            DO i=1,ln-1 
               IF (iaudes(i:i).eq.'0') THEN 
                  iaudes(i:i)=' ' 
               ELSE 
                  EXIT 
               ENDIF
            ENDDO
         ELSEIF(numfield(1:1) .eq. '~') THEN
            DO i=1,4
               numdes(i)=ichar(numfield(i+1:i+1))-55
               if (numdes(i).gt.35) numdes(i) = numdes(i) - 6
            ENDDO
            permdes=numdes(1)*62**3+numdes(2)*62**2+numdes(3)*62+numdes(4)+620000
            WRITE(iaudes,"(i8)") permdes
         ELSE
! alphanumeric asteroid number
            temp=ichar(numfield(1:1))-55
            if (temp.gt.35) temp = temp - 6  
            WRITE(iaudes,133)temp,numfield(2:5)
  133       FORMAT(I2,A4)
         ENDIF 
  123    CALL rmsp(iaudes,ln) 
         RETURN 
      ENDIF
                                                                        
! Unnumbered asteroids                                                  
      desfield=mpccod(6:12) 
                                                                        
      if(desfield(3:3).eq.'S')then 
! Survey Asteroid                                                       
         iaudes=desfield(4:7)//desfield(1:1)//'-'//desfield(2:2) 
                                                                        
      elseif (desfield(1:1).eq.'I' .or.                                 &
     &        desfield(1:1).eq.'J' .or.                                 &
     &        desfield(1:1).eq.'K')then                                 
! Temporary designation                                                 
         if(desfield(5:6).eq.'00')then 
!           1999AA = J99A00A                                            
            tail='' 
         elseif(desfield(5:5).eq.'0')then 
!           1999AA1 = J99A01A                                           
            tail=desfield(6:6) 
         elseif(isnum(desfield(5:5)))then 
!           1999AA12 = J99A12A                                          
            tail=desfield(5:6) 
         else 
!           1999AA103 = J99AA3A                                         
            temp=ichar(desfield(5:5))-55 
!           1999AA363 = J99Aa3A                                         
            if (temp.gt.35) temp = temp - 6  
            write(head,103) temp 
            tail=head//desfield(6:6) 
         endif 
         temp=ichar(desfield(1:1))-55 
         write(head,103) temp 
  103    format(I2) 
         iaudes=head//desfield(2:3)//desfield(4:4)//desfield(7:7)//tail 
      else 
! Unknown type                                                          
         write(*,*)'cannot understand MPC designation: ',mpccod 
         goto 10 
      endif 
      return 
                                                                        
   10 CONTINUE 
      iaudes=mpccod 
      error=.true. 
                                                                        
      END SUBROUTINE iaucod 
                                                                        
! Copyright (C) 1998 by OrbFit Consortium                               
! Version: May 2000 AM MES                                              
! Revised by Genny in Jun 20, 2001 to handle designations as K00Sa3P    
! --------------------------------------------------------------------- 
!                                                                       
!  *****************************************************************    
!  *                                                               *    
!  *                         I A U C O D 2                         *    
!  *                                                               *    
!  * Computes official IAU code from MPC-style packed designation  *    
!  *                                                               *    
!  *****************************************************************    
!                                                                       
! INPUT:    MPCCOD    -  MPC-style packed 7 character code as used in or
!                        5 digits for numbered, 7 characters for unnumbe
!                                                                       
! OUTPUT:   IAUDES    -  IAU code                                       
!           ERROR     -  Error flag (cannot understand input code)      
!                                                                       
      SUBROUTINE iaucod2(mpccod,iaudes,error) 
      IMPLICIT NONE 
                                                                        
      CHARACTER*7 mpccod 
      CHARACTER*9 iaudes 
      LOGICAL error 
                                                                        
      INTEGER ln,i,temp 
      CHARACTER*2 head 
      CHARACTER*3 tail 
                                                                        
      INTEGER lench 
      LOGICAL isnum 
      EXTERNAL lench,isnum 
                                                                        
      CHARACTER*5 numfield 
      CHARACTER*7 desfield 
                                                                        
      error=.false. 
      iaudes=' ' 
                                                                        
      ln=lench(mpccod) 
      IF(ln.LE.0) GOTO 10 
      IF(ln.eq.5)THEN 
! Numbered asteroids                                                    
         numfield=mpccod(1:5) 
         CALL rmsp(numfield,ln) 
         IF(isnum(numfield(1:1)))THEN
! standard number
            iaudes=numfield
! remove leading zeros
            DO i=1,ln-1 
               IF (iaudes(i:i).eq.'0') THEN 
                  iaudes(i:i)=' ' 
               ELSE 
                  EXIT 
               ENDIF
            ENDDO
         ELSE
! alphanumeric asteroid number
            temp=ichar(numfield(1:1))-55
            if (temp.gt.35) temp = temp - 6  
            WRITE(iaudes,133)temp,numfield(2:5)
  133       FORMAT(I2,A4)
         ENDIF 
         CALL rmsp(iaudes,ln) 
! padding 
         DO i=ln+1,9 
           iaudes(i:i)='w' 
         ENDDO 
         RETURN
       ELSEIF(ln.eq.7)THEN 
! Unnumbered asteroids                                                  
          desfield=mpccod(1:7) 
          if(desfield(3:3).eq.'S')then 
! Survey Asteroid                                                       
          iaudes=desfield(4:7)//desfield(1:1)//'-'//desfield(2:2)//'ww' 
                                                                        
          elseif (desfield(1:1).eq.'I' .or.                             &
     &            desfield(1:1).eq.'J' .or.                             &
     &            desfield(1:1).eq.'K')then                             
! Temporary designation                                                 
             if(desfield(5:6).eq.'00')then 
!           1999AA = J99A00A                                            
                tail='www' 
             elseif(desfield(5:5).eq.'0')then 
!           1999AA1 = J99A01A                                           
                tail=desfield(6:6)//'ww' 
             elseif(isnum(desfield(5:5)))then 
!           1999AA12 = J99A12A                                          
                tail=desfield(5:6)//'w' 
             else 
!           1999AA103 = J99AA3A                                         
                temp=ichar(desfield(5:5))-55 
!           1999AA363 = J99Aa3A                                         
                if (temp.gt.35) temp = temp - 6 
                write(head,103) temp 
                tail=head//desfield(6:6) 
             endif 
             temp=ichar(desfield(1:1))-55 
             write(head,103) temp 
  103        format(I2) 
         iaudes=head//desfield(2:3)//desfield(4:4)//desfield(7:7)//tail 
          else 
! Unknown type                                                          
             write(*,*)'cannot understand MPC designation: ',mpccod 
             goto 10 
          endif 
          return 
      ELSE 
! wrong length of designator                                            
          write(*,*) 'designation ',mpccod,' not understood' 
      ENDIF 
   10 CONTINUE 
      iaudes=mpccod 
      error=.true. 
                                                                        
      END SUBROUTINE iaucod2 
                                                                        
! ===============================================                       
!  NIGHTS                                                               
! function computing number of nights of observations                   
! Copyright A. Milani, OrbFit consortium, 21/9/2000                     
! ===============================================                       
INTEGER FUNCTION nights(m,obs,obsw)
  USE astrometric_observations 
  IMPLICIT NONE 
! =========OBSERVATIONS =========================                       
! number of observations                                                
  INTEGER, INTENT(IN) :: m 
! observations new data type
  TYPE(ast_obs),INTENT(IN),DIMENSION(m) :: obs
  TYPE(ast_wbsr),INTENT(IN),DIMENSION(m) :: obsw
! ===========END INTERFACE=======================                       
  double precision t1,t2 
  integer i 
! ===============================================                       
  t1=obs(1)%time_utc 
  nights=1 
  DO i=2,m 
! This is the most trivial algorithm: if there is a 16 hours interval   
! between two observations, they belong to different nights.            
! But a night containing only observations being discarded does not count
     t2=obs(i)%time_utc
     IF(t2-t1.gt.0.66d0.and.obsw(i)%sel_coord.gt.0)THEN 
        nights=nights+1 
        t1=t2 
     ENDIF
  ENDDO
END FUNCTION nights
! Version: December 15, 1997                                            
! Modified November 9, 1998 by Steven Chesley (chesley@dm.unipi.it)     
! in order to handle three digit subscripts in IAU codes.               
! --------------------------------------------------------------------- 
!                                                                       
!  *****************************************************************    
!  *                                                               *    
!  *                         M P C P D S                           *    
!  *                                                               *    
!  * Computes MPC-style packed designation from official IAU code  *    
!  *                                                               *    
!  *****************************************************************    
!                                                                       
! INPUT:    IAUCO    -  IAU code                                        
!                                                                       
! OUTPUT:   MPCCOD    -  MPC-style packed code                          
!           ERROR     -  Error flag (cannot understand input code)      
!                                                                       
SUBROUTINE mpcpds(iauco,mpccod,error) 
  IMPLICIT NONE
  CHARACTER*(*) iauco,mpccod 
  INTEGER number
  LOGICAL error
  INTEGER ln,nd,i,head 
  CHARACTER tsn*1
  INTEGER lench 
  LOGICAL isnum,islett 
  EXTERNAL lench,isnum,islett
  error=.false. 
  mpccod=' '                                                             
  ln=lench(iauco) 
  IF(ln.LE.0) GOTO 10
! Check if numbered asteroids (up to 999999)                                                   
  IF(isnum(iauco(1:ln)) .AND. ln.LE.6) THEN
! numbered
     IF(ln.le.5)THEN 
! numbered < 100000
        DO  i=1,5-ln 
           mpccod(i:i)='0' 
        ENDDO
        mpccod(5-ln+1:5)=iauco(1:ln)
        mpccod(6:7)='  ' 
        RETURN 
     ELSE
! numbered with alphanumeric first digit
        READ(iauco(1:2),*)number
        IF(number.gt.35)THEN
           mpccod(1:1)=char(number+61)
        ELSE
           mpccod(1:1)=char(number+55)
        ENDIF
        mpccod(2:5)=iauco(3:6)
        mpccod(6:7)='  '
     ENDIF
     RETURN
  END IF                                                                      
! Asteroid provisional designations (e.g., 1982QB1)                     
  IF(ln.LT.6 .OR. ln.GT.9) GOTO 2 
  IF(.NOT.isnum(iauco(1:4))) GOTO 2 
  IF(.NOT.islett(iauco(5:6))) GOTO 2 
  nd=ln-6 
  IF(nd.GT.0) THEN 
     IF(.NOT.isnum(iauco(7:ln))) GOTO 2 
  END IF                                                                      
  IF(iauco(1:2).EQ.'19') THEN 
     mpccod(1:1)='J' 
  ELSEIF(iauco(1:2).EQ.'20') THEN 
     mpccod(1:1)='K' 
  ELSEIF(iauco(1:2).EQ.'18') THEN 
     mpccod(1:1)='I' 
  ELSE 
     GOTO 2 
  END IF                                                                      
  mpccod(2:4)=iauco(3:5) 
  IF(nd.EQ.0) THEN 
     mpccod(5:6)='00' 
  ELSEIF(nd.EQ.1) THEN 
     mpccod(5:5)='0' 
     mpccod(6:6)=iauco(7:7) 
  ELSEIF(nd.EQ.2) THEN 
     mpccod(5:6)=iauco(7:8) 
  ELSEIF(nd.EQ.3) THEN 
     read(iauco,103) head 
103  format(6x,i2) 
     IF(head.gt.35)THEN
        mpccod(5:5)=char(head+61)
     ELSE
        mpccod(5:5)=char(head+55)
     ENDIF 
     mpccod(6:6)=iauco(9:9) 
  ELSE 
     GOTO 2 
  END IF
  mpccod(7:7)=iauco(6:6) 
  RETURN 
                                                                        
2 CONTINUE
! Palomar-Leiden survey                                                 
  i=index(iauco(1:ln),'P-L') 
  IF(i.GT.0) THEN 
     IF(ln.NE.i+2) GOTO 3 
     nd=i-1 
     IF(nd.LT.1) GOTO 3 
     mpccod(1:3)='PLS' 
     DO i=1,4-nd 
        mpccod(3+i:3+i)='0' 
     ENDDO
     mpccod(8-nd:7)=iauco(1:nd) 
     RETURN 
  END IF
3 CONTINUE
! Trojan surveys                                                        
  i=index(iauco(1:ln),'T-') 
  IF(i.GT.0) THEN 
     IF(ln.NE.i+2) GOTO 5 
     tsn=iauco(i+2:i+2) 
     IF(tsn.NE.'1' .AND. tsn.NE.'2' .AND. tsn.NE.'3') GOTO 5 
     mpccod(1:1)='T' 
     mpccod(2:2)=tsn 
     mpccod(3:3)='S' 
     nd=i-1 
     DO i=1,4-nd 
        mpccod(3+i:3+i)='0' 
     ENDDO
     mpccod(8-nd:7)=iauco(1:nd) 
     RETURN 
  END IF
5 CONTINUE
! Cannot understand input code                                          
10 CONTINUE 
  mpccod=iauco 
  error=.true. 
END SUBROUTINE mpcpds

! ======================================================================
! OBSSTA                                                                
SUBROUTINE obssta(iunrep,name0,name1,tut,typ,m,nig,change)
  USE name_rules, ONLY: name_len 
  IMPLICIT NONE 
  INTEGER iunrep,m,nig 
  CHARACTER*(name_len) name0,name1 
  CHARACTER*1 typ(m)
  DOUBLE PRECISION tut(m),tmin,tmax,hour,day,arc 
  INTEGER iday,month,yearm,yearx 
  INTEGER i,nr,nrout 
! change flag                                                           
  LOGICAL change 
  nr=0 
  nrout=0 
  tmin=1d9 
  tmax=-1.d9 
  DO i=1,m 
     IF(typ(i).eq.'R'.or.typ(i).eq.'V')THEN
        nr=nr+1 
     ENDIF
     IF(tut(i).gt.tmax)tmax=tut(i) 
     IF(tut(i).lt.tmin)tmin=tut(i) 
  ENDDO
  arc=tmax-tmin 
!     convert time                                                      
  CALL mjddat(tmin,iday,month,yearm,hour) 
  day=iday+hour/24 
  CALL mjddat(tmax,iday,month,yearx,hour) 
  day=iday+hour/24 
  IF(name0.eq.name1)THEN 
     WRITE(iunrep,101)name0,m,arc,nig,yearm,yearx,nr,tmax,change 
101  FORMAT(a9,1x,i4,1x,f8.2,1x,i4,1x,i4,1x,i4,1x,i3,1x,f10.3,1x,l1) 
  ELSE 
! output                                                                
     WRITE(iunrep,100)name0,name1,m,arc,nig,yearm,yearx,nr,tmax 
100  FORMAT(a9,1x,a9,1x,i4,1x,f8.2,1x,i4,1x,i4,1x,i4,1x,i3,1x,f10.3) 
  ENDIF
END SUBROUTINE obssta

! =======================================                               
! QUALITY                                                               
! orbit quality as guessed from the observational adta available        
! WARNING: this does not take into account the value of the             
!          orbital elements, as it would be necessary, because          
!          this function is normally used before the orbital            
!          elements are available                                       
! ========================================                              
INTEGER FUNCTION quality(arc,m,nig,minoss) 
  USE output_control
  IMPLICIT NONE 
! ==========INPUT =============                                         
! no. obs, nights                                                       
  INTEGER,INTENT(IN) :: m,nig 
! arc length (in days)                                                  
  DOUBLE PRECISION,INTENT(IN) :: arc 
! ========== OUTPUT ==========                                          
! minimum number of observations                                        
  INTEGER, INTENT(OUT) :: minoss 
! =====END INTERFACE==========                                          
  IF(arc.ge.180.d0.and.m.ge.10.and.nig.ge.5)THEN 
! multi apparition
     quality=1 
     minoss=10 
!  ELSEIF(arc.ge.20.d0.and.m.ge.6.and.nig.ge.3)THEN 
!     quality=2 
!     minoss=6 
  ELSEIF(arc.ge.20.d0.and.m.ge.8.and.nig.ge.4)THEN
! reliable single opposition 
     quality=2 
     minoss=8 
  ELSEIF(arc.ge.6.5d0.and.m.ge.6.and.nig.ge.3)THEN 
! three nights
     quality=3 
     minoss=6 
  ELSEIF(arc.ge.1.5d0.and.m.ge.6.and.nig.ge.3)THEN
! short arc three nights 
     quality=4 
     minoss=6 
  ELSEIF(arc.ge.6.5d0.and.m.ge.4.and.nig.ge.2)THEN
! well separated two nights 
     quality=5 
     minoss=4 
!  ELSEIF(arc.ge.1.5d0.and.m.ge.4.and.nig.ge.2)THEN 
!     quality=6 
!     minoss=4 
  ELSEIF(m.ge.3.and.nig.ge.2)THEN 
! two nights 
     quality=6 
     minoss=4 
  ELSEIF(m.ge.3)THEN
! one night 
     quality=7 
     minoss=3 
  ELSEIF(m.ge.2)THEN
! attributable at most 
     quality=8 
     minoss=2 
  ELSEIF(m.eq.1)THEN 
! single observations
     quality=9 
     minoss=1 
  ELSEIF(m.eq.0)THEN
! messy file with only degraded observations
! occurs only when quality is called with argument mgood instead of m
     quality=10
     minoss=0
  ELSE
     WRITE(ierrou,*)' quality: case not considered, arc, m,nig ',arc,m,nig
     numerr=numerr+1
     WRITE(*,199)arc,m,nig
199  FORMAT(' quality: case not considered, arc, m,nig ',f8.2,1x,i4,1x,i3)
     quality=10 
     minoss=0  
  ENDIF
END FUNCTION quality

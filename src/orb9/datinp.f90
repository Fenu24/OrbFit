!!$ -*- f90 -*-
! ***************************************************                   
!  {\bf datinp} input data for giffv                                    
      SUBROUTINE datinp (nast, ifo, iast, nia, ilce, sv, tv, nco, nin,  &
      coox, unitx, uniy, inter, ndat) 

      USE var_precision
      USE parameters_giff
      USE prompting
      IMPLICIT REAL(KIND=r_kind)  (a - h, o - z) 
      PARAMETER (nnx = 50) 
! system descriptors                                                    
      CHARACTER(3) coox, unitx, uniy 
      LOGICAL error
! input data                                                            
      DIMENSION tv (nx), sv (nx, ncox), iast (nplx), x (nnx), ico (nnx),&
      ng (nnx)
! time difference default tolerance     
      if(nin.eq.0)dtt=100*epsilon(1._r_kind)
      IF (ndat.gt.0) then 
!  non orbit8v format; just one line with time, ndat real variables     
         nang = 0 
         nact = ndat 
         isup = 0 
         nel = ndat 
         nia = 1 
         uconv = 1.d0 
      ELSEIF (ndat.eq.0) then 
!  orbit8v format:                                                      
!  number of variables per body                                         
         nang = numang (coox) 
         nact = numact (coox) 
!  unit conversion for angles                                           
         CALL convun_giff (unitx, uniy, uconv, urev,error) 
!  supplementary angle (orbit8t only)                                   
!         if(inter.ne.0)then                                            
!            CALL prompt(' supplementary angle? 1=yes 0=no',isup,0,1)
!         else                                                          
!            call reaint(4,'isup',isup)
!         endif
!  the next line is added; NAG compiler does not like undefined variables
         isup=0  ! ADDED                                                       
         nang = nang + isup 
         IF (ilce.eq.0) then 
            nel = 6 + isup 
         ELSE 
            nel = 6 + isup + 1 
         ENDIF 
      ELSEIF (ndat.lt.0) then 
         nang = 0 
         nact = - ndat 
         isup = 0 
         nel = - ndat 
         nia = 1 
         CALL convun_giff (unitx, uniy, uconv, urev, error) 
      ENDIF 
!  choice of components                                                 
   97 IF (inter.ne.0) then 
         PRINT *, ' which component? 1 per line, 0=end, -1=all'
      ENDIF 
      nc = 0 
      DO 1 i = 1, nnx 
         IF (inter.ne.0) then
            if(i.eq.1)then 
               ic=-1
            else
               ic=0
            endif
            CALL prompt(' component ',ic,-1)  
         ELSE 
            CALL reaint (4, 'ic', ic) 
         ENDIF 
         IF (ic.eq.0) then 
            GOTO 2 
         ELSEIF (ic.eq. -1) then 
            nc = nel 
            DO n = 1, nel 
               ico (n) = n
            END DO 
            GOTO 2 
         ELSEIF (ic.gt.0) then 
            nc = nc + 1 
            ico (nc) = ic 
         ENDIF 
    1 END DO 
      WRITE ( * , * ) ' too many, max is ', nnx 
    2 CONTINUE 
!  check for memory space                                               
      ncon = nco + nc * nia 
      IF (ncon.gt.ncox) then 
         WRITE ( * , * ) ' too many components, already ', nco 
         WRITE ( * , * ) ' maximum is ', ncox, ' change pargif.h' 
         GOTO 97 
      ENDIF 
!  read loop                                                            
      DO 10 n = 1, nx 
         IF (ifo.eq.1) then 
            READ (1, *, end = 3, err = 4) t 
!c         write(*,*)'t=',t                                             
            na = 1 
            DO 11 j = 1, nast 
               IF (ilce.ge.j) then 
                  READ (1, *, end = 4, err = 70) (x (i), i = 1, nact),  &
                  (x (i + nact), ng (i), i = 1, nang), x (1 + nact +    &
                  nang)                                                 
               ELSEIF (ilce.lt.j) then 
                  READ (1, *, end = 4, err = 70) (x (i), i = 1, nact),  &
                  (x (i + nact), ng (i), i = 1, nang)                   
               ENDIF 
               GOTO 80 
   70          DO 77 nn = 1, nact 
                  x (nn) = 0.d0
   77          END DO 
               DO 78 nn = 1, nang 
                  x (nn + nact) = 0.d0 
                  ng (nn) = 0
   78          END DO 
               IF (ilce.ge.j) then 
                  x (1 + nact + nang) = 0.d0 
               ENDIF 
   80          CONTINUE 
!c              write(*,*)'body no.',j                                  
!c              write(*,*)(x(i),i=1,nact),                              
!c   +              (x(i+nact),ng(i),i=1,nang)                          
               DO 13 i = 1, nang 
                 x (i + nact) = x (i + nact) * uconv + ng (i) * urev
   13          END DO 
!  is this asteroid selected?                                           
               IF (na.gt.nia) goto 11 
               IF (j.eq.iast (na) ) then 
                  DO 12 i = 1, nc 
                     ind = i + (na - 1) * nc + nco 
                     sv (n, ind) = x (ico (i) )
   12             END DO 
                  na = na + 1 
               ENDIF 
   11       END DO 
         ELSEIF (ifo.eq.0) then 
            IF (ilce.ge.1) then 
               READ (1, *, end = 3, err = 4) t, (x (i), i = 1, nact),   &
               (x (i + nact), ng (i), i = 1, nang), x (1 + nact + nang) 
            ELSEIF (ilce.lt.1) then 
               READ (1, *, end = 3, err = 4) t, (x (i), i = 1, nact),   &
               (x (i + nact), ng (i), i = 1, nang)                      
            ENDIF 
            DO 14 i = 1, nang 
               x (i + nact) = x (i + nact) * uconv + ng (i) * urev
   14       END DO 
!  this asteroid is selected anyway, being alone                        
            DO 15 i = 1, nc 
               ind = i + nco 
               sv (n, ind) = x (ico (i) )
   15       END DO       
         ENDIF 
         IF (nin.ne.0) then 
            IF (abs(t-tv(n)).gt.dtt.and.n.le.nin) then 
      WRITE ( * ,  * ) ' inconsistent list of times, record no. ', n 
               WRITE (*,*)' old list, t=', tv (n),' new list t=',t 
               call prompt(' time difference tolerance',dtt)
            ENDIF 
         ENDIF 
         tv (n) = t 
   10 END DO 
      WRITE ( * , * ) ' too many records, max was ', nx 
      GOTO 3 
    4 WRITE ( * , * ) ' error in input' 
    3 nin1 = n - 1 
      nco = nco + nc * nia 
      WRITE ( * , * ) ' read ', nin1, ' records, nco=', nco 
      IF (nin1.lt.nin.and.nin.ne.0) then 
         WRITE ( * , * ) ' new file is shorter than the previous one' 
         WRITE ( * , * ) ' input truncated at record ', nin1, ' time ', &
         tv (nin1)                                                      
         nin = nin1 
      ELSEIF (nin1.gt.nin.and.nin.ne.0) then 
         WRITE ( * , * ) ' new file is longer than the previous one' 
         WRITE ( * , * ) ' input truncated at record ', nin, ' time ',  &
         tv (nin)                                                       
         nin = nin 
      ELSE 
         nin = nin1 
      ENDIF 
      RETURN 
      END SUBROUTINE datinp                         
! ***********************************************                       
! {\bf convun_giff} unit conversion for angles                               
      SUBROUTINE convun_giff (unitx, unity, uconv, urev,error) 
      USE var_precision
      USE fund_const
      IMPLICIT REAL(KIND=r_kind)  (a - h, o - z) 
      LOGICAL error
      CHARACTER(3) unitx, unity 
      error=.FALSE.
      IF (unitx.eq.'REV') then 
         uconv = 1.d0 
      ELSEIF (unitx.eq.'RAD') then 
         uconv = 1.d0 / dpig 
      ELSEIF (unitx.eq.'DEG') then 
         uconv = 1 / 3.6d2 
      ELSE 
         WRITE ( * , * ) ' input unit ', unitx, ' not supported' 
         error=.TRUE.
         RETURN
      ENDIF 
      IF (unity.eq.'REV') then 
         uconv = uconv * 1.d0 
         urev = 1.d0 
      ELSEIF (unity.eq.'RAD') then 
         uconv = uconv * dpig 
         urev = dpig 
      ELSEIF (unity.eq.'DEG') then 
         uconv = uconv * 3.6d2 
         urev = 3.6d2 
      ELSE 
         WRITE ( * , * ) ' output unit ', unity, ' not supported' 
         error=.TRUE.
         RETURN
      ENDIF 
      RETURN 
      END SUBROUTINE convun_giff                         
! ======================================================                
!  {\bf jobinv} input and extraction job definition                     
!  to be used by  giffv                                                 
      SUBROUTINE jobinv (inter, nast, imas, ifo, number, iast, nia,     &
      coox, sysx, refx, unitx, ilce, uniy, ndat) 
      USE var_precision
      USE parameters_giff
      USE prompting
      IMPLICIT REAL(KIND=r_kind)  (a - h, o - z)   
!     PARAMETER (nastx = nplx) 
      DIMENSION iast (nplx) 
!  character arrays                                                     
      CHARACTER(40) infile 
      CHARACTER(3) coox, sysx, unitx, uniy 
      CHARACTER(6) refx 
      CHARACTER(60) commen, comme1 
      CHARACTER(9) number (nplx), numbin (nplx), num, nnum, num1 
      CHARACTER colhea * 100 
      DIMENSION gm (nplx) 
      LOGICAL error
! *************************************************************         
!  open and input file description  
 201  CONTINUE                                    
      IF (inter.ne.0) then 
!  next line ADDED to avoid garbage in the prompt
         infile=' '
         CALL prompt(' input  file ?',infile) 
      ELSE 
         CALL reastr (4, 'infile', infile) 
         WRITE ( * , * ) ' input from file ', infile 
      ENDIF 
      OPEN (1, file = infile, status = 'old', IOSTAT=iost)
      IF(iost.ne.0)THEN
         PRINT *, ' the file ',infile, ' is not available'
         GOTO 201 
      ENDIF
      IF (inter.ne.0) then
!  next line ADDED to avoid garbage in the prompt
           ihead=0 
      CALL prompt(' header provided ?  0=no 1=orbit8/9 2=lunpro', &
           ihead, 0, 2) 
      ELSE 
         CALL reaint (4, 'ihead', ihead) 
      ENDIF 
      IF (ihead.eq.1) then 
         ndat = 0 
         CALL reahea (1, numbin, commen, colhea, coox, sysx, refx,      &
         unitx, nast, ilce, t0)                                         
         WRITE ( *, 101) coox, sysx, refx, unitx, nast, ilce, imas 
  101 FORMAT  (' coordinate system: ',1x,a3,1x,a3,1x,a6,1x,a3/          &
     &  ' no. orbits in file =',i3,' ; first ',i3,                      &
     &  ' with LCE, mass flag ',i2)                                     
!  asteroid names in input                                              
         IF (inter.ne.0) then 
            DO 45 iia = 1, nast 
               WRITE ( * , * ) numbin (iia) , ' available' 
   45       END DO 
      CALL prompt(' with mass ? 1=yes(planets) 0=no(asteroids)', &
           imas, 0, 1) 
         ELSE 
            CALL reaint (4, 'imas', imas) 
         ENDIF 
         IF (imas.eq.1) then 
            nbod = nast + 1 
            CALL skip(1,1)
            READ (1, * ) (gm (j), j = 1, nbod) 
         ENDIF 
!  asteroid names in input                                              
         DO  iia = 1, nast 
            WRITE ( * , * ) numbin (iia) , ' available' 
         END DO 
! header for selenocentric elements                                     
      ELSEIF (ihead.eq.2) then 
         !CALL reastr (1, 'coox', coox, commen,icom) 
         !CALL reastr (1, 'unitx', unitx, comme1,icom) 
      WRITE ( * ,  * ) ' coordinates ', coox, '  ', commen 
      WRITE ( * ,  * ) ' units for angles ', unitx, '   ', comme1 
         ndat = - 6 ! trick: set ndat to - its value, to signal that there are
                    ! no angular data; this is not clean.
         nast = 1 
      ELSEIF (ihead.eq.0) then 
! header not available; self service                                    
         IF (inter.ne.0) then
! ADDED to avoid garbage in the prompt
            nh=0 
            CALL prompt(' number of header lines?',nh) 
            CALL prompt(' no. data per line? (not including time)', ndat) 
         ELSE 
            CALL reaint (4, 'nh', nh) 
            CALL reaint (4, 'ndat', ndat) 
         ENDIF 
         CALL skip (1, nh) 
         ilce = 0 
         unitx = 'RAD' 
         nast=1
      ENDIF
!  choice between two formats: time,record if one orbit (anyway with    
!  free format reading the newline does not matter),                    
!  time/records for more orbits                                         
      IF (nast.gt.1) then 
         ifo = 1 
      ELSEIF (nast.eq.1) then 
         ifo = 0 
      ELSE 
         WRITE ( * , * ) ' no. orbits ', nast, ' not allowed' 
         STOP 
      ENDIF 
      IF (ihead.eq.1)THEN 
!  selection of the orbit among many 
         iia = 1 
   43    CONTINUE 
         IF (inter.ne.0) THEN
            num=' ' 
            CALL prompt(' which body among the above ?',num) 
         ELSE 
            CALL reastr (4, 'num', num) 
         ENDIF 
         num1 = num 
         CALL remove_spaces (num1, len) 
         IF (len.le.0) goto 47 
         DO 46 i = 1, nast 
            nnum = numbin (i) 
            IF (num.eq.nnum) THEN 
               DO  j = 1, iia - 1 
                  nnum = number (j) 
                  IF (num.eq.nnum) then 
                     WRITE ( * , * ) num, ' already selected' 
                     GOTO 43 
                  ENDIF 
               END DO 
               WRITE ( * , * ) num, ' OK' 
               iast (iia) = i 
               number (iia) = num 
               iia = iia + 1 
               GOTO 43 
            ENDIF 
   46    END DO 
         WRITE ( * , * ) num, ' not available' 
         GOTO 43 
   47    nia = iia - 1
      ENDIF 
!                                                                       
! units for output                                                      
      IF (inter.ne.0) then 
 50      uniy='RAD'
         CALL prompt(' units for angles in output, RAD, DEG, REV?', &
              uniy, case=2) 
! next lines commented as prompting.f does not support "error" in the
! case of string_prompt; the quick fix solution
!         IF(error)THEN
!            PRINT *, 'unit for angles not known ', uniy
!            GOTO 50
!         ENDIF
      ELSE 
         CALL reastr (4, 'uniy', uniy) 
      ENDIF 
      CALL convun_giff ('RAD', uniy, uconv, urev, error)
      RETURN 
      END SUBROUTINE jobinv                         


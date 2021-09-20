!!$ -*- f90 -*-
! *************************************************************         
!  {\bf giff}                                                           
!                                                                       
! general purpose plot with graphic interactive filtering               
!                                                                       
! orbit9/GNUPLOT version, Pisa, April 1995                              
! *************************************************************         
      PROGRAM giff 
      USE var_precision
      USE fund_const
      USE fir_filters
      USE parameters_giff
      USE prompting
      IMPLICIT REAL(KIND=r_kind)  (a - h, o - z) 
      REAL (KIND=r_kind) xf, yf(2)
! time                                                                  
      DIMENSION tv (nx), x (nx), isel (nx) 
! input data                                                            
      DIMENSION sv (nx, ncox) 
! selected argument                                                     
      DIMENSION crit (nx), y (nx), dy (nx), s (nx) 
! data for filtering, spectra                                           
      DIMENSION amp (nx), tper (nx) 
! control                                                               
      LOGICAL indev, outf 
! system descriptors                                                    
      CHARACTER(3) coox, sysx, unitx, uniy 
      CHARACTER(6) refx 
! names and locations in input file                                     
      CHARACTER(9) num (nplx) 
      DIMENSION iast (nplx) 
! labels                                                                
      CHARACTER ylab * 80, xlab * 80, xlab1 * 80, title * 80 
! filenames                                                             
      CHARACTER spefile * 60, resfile * 60, outfile * 60, postfi * 6 
! *************************************************************   
! device initialisations                                                
      indev = .true. 
      ipag = 0 
      idev=3
      CALL prompt(' device 3=X11 4=Tektr -5=ps file -6=laser pri. -7=color pri.', &
           idev, -7, 6)
! *************************************************************
!  begin again
  999 CONTINUE
! *************************************************************         
! input  mode 
      inter=1                                                          
!      CALL prompt(' 1=interactive 0=from giffv.inp', inter, 0, 1) 
!      IF (inter.eq.0) then 
!         OPEN (4, file = 'giff.inp', status = 'unknown') 
!         WRITE ( * , * ) ' reading from file giff.inp' 
!         CALL skip (4, 1) 
!      ENDIF 
! *************************************************************         
! title                                                                 
      IF (inter.eq.0) then 
         CALL reastr (4, 'title', title) 
         CALL reastr (4, 'postfi', postfi) 
         CALL reaint (4, 'isou', isou) 
      ELSE 
         title=' ';CALL prompt(' plot title?', title) 
         postfi=' ';CALL prompt(' postfix for output files ?', postfi) 
         CALL prompt(' x ax label 1=orbit9 2=time 3=each plot 0=no', &
              isou, 0, 3) 
      ENDIF 
! open session log file                                                 
      outfile = 'out.'//postfi 
      CALL remove_spaces (outfile, len) 
      OPEN (10, file = outfile (1:len) , status = 'unknown') 
! x axis labels                                                         
      IF (isou.eq.1) then 
         xlab = 'Time (y) - Orbit9 integration' 
         xlab1 = 'Period (y) - Orbit9 integration' 
      ELSEIF (isou.eq.2) then 
         xlab = 'Time (years)' 
         xlab1 = 'Period (years)' 
      ELSEIF (isou.eq.0) then 
         xlab = ' ' 
         xlab1 = ' ' 
      ENDIF 
! *************************************************************         
!  job input                                                            
      nco = 0 
      nin = 0 
      DO 1 nfile = 1, ncox 
         CALL jobinv (inter, nast, imas, ifo, num, iast, nia, coox,     &
         sysx, refx, unitx, ilce, uniy, ndat)                           
! *************************************************************         
! input data                                                            
         CALL datinp (nast, ifo, iast, nia, ilce, sv, tv, nco, nin,     &
         coox, unitx, uniy, inter, ndat)                                
         CLOSE (1) 
!  input finished?                                                      
         IF (inter.eq.0) then 
            CALL reaint (4, 'imore', imore) 
         ELSE 
            imore=0
            CALL prompt(' more? 1=yes 0=no', imore, 0, 1) 
         ENDIF 
         IF (imore.eq.0) goto 2 
    1 END DO 
! *************************************************************         
! input complete                                                        
    2 WRITE ( *, 150) nin, nco 
      WRITE (10, 150) nin, nco 
  150 FORMAT(' input complete, records=',i6,' with ',i4,' components'/) 
! *************************************************************         
! main loop                                                             
   11 CONTINUE 
! choice of time interval                                               
      CALL timsel (nin, tv, x, isel, nn) 
! *************************************************************         
! choice of argument                                                    
      CALL argsel (nin, nn, nx, nco, x, sv, crit, y, isel, uniy) 
! *************************************************************         
!   choice between plot, filter and spectrum                            
!   menu and function loop                                              
   89 CONTINUE 
      WRITE ( *, 188) 
  188 FORMAT(' 1=plot 2=spectrum 3=filter 4=linfit 5=smooth'/           &
     &   ' 6=out 7=delta,sigma 8=store 9=change device '/               &
     &' 0=change arg -1=stop')                                          
      READ ( *, * ) isp 
      WRITE ( *, 667) isp 
  667 FORMAT(' isp= ',i4) 
      WRITE (10, 667) isp 
! main menu
menu: SELECT CASE(isp)   
      CASE(0) menu  !  change argument
          goto 11 
! **************************************************************        
      CASE(1) menu ! plot of argument    
!  do you want to use the original data, or the transformed ones?  
          call traori (iav, crit, y, s, nn)                            
! device initialization                                                 
         IF (indev) call grafin (idev) 
         indev = .false. 
! new page                                                              
         CALL newpag (idev, ipag) 
! axis labels, if device= printer
         IF (idev.le.0) then 
            IF (isou.eq.3) then  ! if x axis label is not constant
               CALL prompt (' x axis label?', xlab) 
            ENDIF 
! y axis label                                                          
            CALL prompt(' y axis label?', ylab) 
         ENDIF 
! plot
         istyle =1 ! default is continous line
         CALL prompt(' plot style? 0=points 1=line -n=thick points', &
              istyle, -10, 1)
         CALL plotcr (x, s, nn, xlab, ylab, title, idev, istyle) 
         CALL ploflu (idev)
! **************************************************************        
      CASE(2) menu ! spectrum
!  do you want to use the original data, or the transformed ones?       
         call traori (iav, crit, y, s, nn)
         PRINT *, ' period min, max, no. freq. ; 1=print,0=no' 
         READ ( *, * ) pmin, pmax, nfre, iprint 
         WRITE ( *, 668) pmin, pmax, nfre, iprint 
  668 FORMAT(' from ',f12.1,'  to  ',f12.1 ,' no. fr',i5,' ipr',i2) 
         WRITE (10, 668) pmin, pmax, nfre, iprint 
         IF (iprint.ne.0) then 
            spefile = 'spe.'//postfi 
            CALL remove_spaces(spefile, len) 
            OPEN (8, file = spefile (1:len) , status = 'unknown') 
            REWIND 8 
         ENDIF 
         CALL prompt(' 0=ampl  1=sqrt(dens.spett.)', iam, 0, 1) 
         dper = (pmax - pmin) / nfre 
         nfre1 = nfre+1 
         DO  j = 1, nfre1 
            per = pmin + dper * (j - 1) 
            itest = 0 
            CALL peri2 (x, y, nn, per, sp, itest, d0, d1, d2)
            a = d1 * d1 + d2 * d2 
            a = sqrt (a) 
            f = atan2 (d1, d2) 
            cost = + d0 
            IF (iprint.ne.0) write (8, 197) per, sp, a, f, cost 
            tper (j) = per 
            amp (j) = a 
            IF (iam.gt.0) amp (j) = sqrt (sp) 
         END DO 
         IF (iprint.ne.0) close (8) 
! **************************************************************        
! plot of spectrum                                                      
   75    CONTINUE 
! device initialization                                                 
         IF (indev) call grafin (idev) 
         indev = .false. 
! new page                                                              
         CALL newpag (idev, ipag) 
         IF (idev.lt.0) then 
! x axis label                                                          
            IF (isou.eq.3) then 
               CALL prompt (' x axis label?', xlab1) 
            ENDIF 
! y axis label                                                          
            CALL prompt( ' y axis label?', ylab) 
         ENDIF 
! plot                                                                  
         istyle = 1 
         CALL plotcr (tper, amp, nfre1, xlab1, ylab, title, idev, istyle)
         CALL ploflu (idev) 
!  what next?                                                           
         IF (iprint.ne.0) close (8) 
            ndev=0 !default
            CALL prompt(' change device? -5=ps file -6=laser 3=X11 4=TeK 0=as before', &
                 ndev, -7, 6)
         IF (ndev.ne.0) then 
            CALL grafcl (idev) 
            indev = .true. 
            ipag = 0 
            idev = ndev 
            GOTO 75 
         ENDIF 
! **************************************************************
      CASE(3) menu ! filtering : Ferraz-Mello single frequency filter 
!  do you want to use the original data, or the transformed ones?       
         call traori (iav, crit, y, s, nn)   
         WRITE ( *, 194) ! choice of frequency to be filtered out      
  194 FORMAT   (' output: period,rel.sp.dens.,amp,phas,const/sig'/) 
         WRITE (10, 194) 
   70    CONTINUE 
         WRITE ( *, 193) 
  193 FORMAT   (' period, itest(1=filter,0=no); per<0 to stop' ) 
         READ ( *, * ) per, itest 
         IF (per.le.0.d0) goto 89 
         CALL peri2 (x, s,  nn, per, sp, itest, d0, d1, d2) 
         a = d1 * d1 + d2 * d2 
         a = sqrt (a) 
         f = atan2 (d1, d2) 
         cost = d0 
         WRITE ( *, 197) per, sp * 100, a, f, cost, itest 
  197 FORMAT   (f14.4,f8.4,1p,d13.6,0p,f9.5,1p,d13.6,i2) 
         WRITE (10, 197) per, sp, a, f, cost, itest 
         IF (itest.ne.0) then 
            sm=sum(s(1:nn))/nn
            sig = 0.d0 
            DO  i = 1, nn 
               sig = sig + (s (i) - sm) **2 
               y (i) = s (i) 
            END DO 
            sig = sqrt (sig / nn) 
            WRITE ( *, 198) sig 
  198 FORMAT      (1p,d15.5) 
            WRITE (10, 198) sig 
         ENDIF 
         GOTO 70 
! *************************************************************         
      CASE(4) menu !  linear fit
!  do you want to use the original data, or the transformed ones?       
         call traori (iav, crit, y, s, nn)
         CALL linfi3 (x, s, rate, rms, cost, dy, nn) 
         WRITE ( *, 941) rate, cost, rms, iav 
  941 FORMAT   (' rate=',1p,d18.10,'  cost=',d15.7,'  rms=',d12.4,i5) 
         WRITE (10, 941) rate, cost, rms, iav 
         y (1:nn) = dy (1:nn)
! *************************************************************         
      CASE(5) menu !  digital smoothing
!  do you want to use the original data, or the transformed ones?       
         call traori (iav, crit, y, s, nn)
         WRITE ( * , * ) ' sampling ratio?' 
         READ ( *, * ) isamp 
         WRITE (10, 959) isamp 
  959 FORMAT   (' digital smoothing, decimation ',i4) 
!  coefficienti del filtro  
         CALL filope(isamp,1)
!  passo assunto costante
         fstep= x(nin/2+1) - x(nin/2)
! filter half length
         hlen=fstep*(nfil-1)/2 
!  digital low frequency pass filter                                    
         nskip = 0 
         j = 0 
         DO 950 i = 1, nn 
            ndat=1
            CALL filter (x(i), s(i:i), ndat, nskip, fstep, xf, yf(1:1), outf)
            IF (outf) then 
               j = j + 1 
               y (j) = yf(1)
               x (j) = xf
            ENDIF 
  950    END DO 
         nn = j 
         WRITE ( *, 953) nn 
  953 FORMAT   (' available data points no.', i5/                       &
     &    '  warning: times are changed; do not use orig.data')  
! *************************************************************         
      CASE(6) menu !  output data as they are now in file res.postfi  
!  do you want to use the original data, or the transformed ones?       
         call traori (iav, crit, y, s, nn)
         resfile = 'res.'//postfi 
         CALL remove_spaces (resfile, len) 
         OPEN (11, file = resfile (1:len) , status = 'unknown') 
         DO 963 i = 1, nn 
          WRITE (11, 962) x (i), s (i)
  963    END DO           ! ADDED
  962 FORMAT   (f15.3,e18.9) 
         CLOSE (11) 
! *************************************************************  
      CASE(7) menu !  average and sigma ,max and min 
!  do you want to use the original data, or the transformed ones?       
         call traori (iav, crit, y, s, nn)
         swy = 0.d0 
         smax = maxval(s(1:nn))
         smin = minval(s(1:nn))
         swy = sum(s(1:nn))
         ym = swy / nn 
         del = smax - smin 
         WRITE ( *, 181) ym, smax, smin, del 
  181 FORMAT(' av= ',1p,d15.7,' max=',d15.7,' min=',d15.7,' del=',d15.7)
         WRITE (10, 181) ym 
         sig = 0.d0 
         ii = 0 
         DO  i = 1, nn 
            sig = sig + (s (i) - ym) **2 
         END DO 
         sig = sqrt (sig / nn) 
         WRITE ( *, 182) sig 
  182 FORMAT   (' sigma=',1p,d13.5) 
         WRITE (10, 182) sig 
! *************************************************************  
      CASE(8) menu !  store current argument                                   
         IF (nco.eq.ncox) then 
      WRITE ( * ,  * ) ' not enough room, already ', ncox, ' components' 
            GOTO 89 
         ELSEIF (nin.ne.nn) then 
            WRITE ( * , * ) ' this is not allowed, you can only store' 
            WRITE ( * , * ) ' over the entire time span' 
            GOTO 89 
         ENDIF 
!  do you want to use the original data, or the transformed ones?       
         call traori (iav, crit, y, s, nn)
         nco = nco + 1 
         sv (1:nn, nco) = s (1:nn) 
         WRITE ( * , * ) ' stored in column ', nco 
! *************************************************************         
      CASE(9) menu !  change device                                            
         CALL prompt(' change device? -5=ps file -6=laser 3=X11 4=Tek 0=as before', &
              idev, -7, 6)
         IF(idev.ne.0)THEN
            CALL grafcl (idev) 
            indev = .true.   
            ipag = 0   
         ENDIF        
! **************************************************************        
      CASE(-1) menu ! close and/or change file
         CALL PROMPT(' change input file? 1=yes 0=stop', &
              ibegin, 0, 1)
         CLOSE(10)
         IF(ibegin.eq.1)THEN
            goto 999
         ENDIF 
         STOP 
! *************************************************************         
! input error                                                           
      CASE DEFAULT  menu
         WRITE ( * , * ) ' option ', isp, ' unknown' 
      END SELECT menu
      GOTO 89 
      END PROGRAM giff                              
! *************************************************************         
!   {\bf argsel} choice of argument                                     
      SUBROUTINE argsel (nin, nn, nx, nco, x, sv, crit, y, isel, uniy)
      USE var_precision
      USE fund_const
      USE prompting
      IMPLICIT real(kind=r_kind) (a - h, o - z)
      PARAMETER (ncox = 100) 
      CHARACTER(3) uniy 
      LOGICAL error
      DIMENSION x (nn), crit (nn), sv (nx, nco), isel (nin), y (nn),    &
      k (ncox)
   40 CONTINUE  ! just to pass control by later GOTO's
      CALL prompt(' which data to plot? 0=comb -1=two arg -2=ratio & 
         & -3 polar radius -4 polar angle -5=cart from polar', js, -5, nco)
      WRITE ( *, 600) js 
  600 FORMAT(' data choice: ',i5) 
      WRITE (10, 600) 
! ************************************
      IF(js.lt.0) THEN  ! plot of two coordinates
!  to handle a couple of coordinates (instead of time and one coordinate
!                                                                       
!  select the two coordinates 
         WRITE ( * , * ) ' two indices?' 
         j1=0 ; j2=0 !no defaults, the range may have changed
         CALL prompt(' first coordinate ', j1, 1, nco)
         CALL prompt(' second coordinate ', j2, 1, nco)
axes:    SELECT CASE(js) 
! ************************************      
         CASE (-1) axes ! plot two coordinates as they are                    
               WRITE (10, 1410) j1, j2 
 1410          FORMAT(' 2 components, no.',2i5) 
               ii = 0 
               DO i = 1, nin 
                  IF (isel (i) .eq.1) then 
                     ii = ii + 1 
                     x (ii) = sv (i, j1) 
                     crit (ii) = sv (i, j2) 
                  ENDIF 
               END DO 
! ************************************ 
         CASE(-2) axes  !  ratio of two coordinates    
            WRITE (10, 1411) j1, j2 
 1411       FORMAT (' ratio of ',i5,'  and ',i5) 
            ii = 0 
            izer = 0 
            DO  i = 1, nin 
               IF (isel (i) .eq.1) then 
                  ii = ii + 1 
                  IF (sv (i, j2) .eq.0.d0) then 
                     izer = izer + 1 
                     ii = ii - 1 
                  ELSE 
                     crit (ii) = sv (i, j1) / sv (i, j2) 
                  ENDIF 
               ENDIF 
            END DO 
            IF (izer.gt.0) write ( * , * ) ' zero divide no. ', izer 
! ************************************
         CASE( -3) axes  !  polar radius 
            WRITE (10, 1412) j1, j2 
 1412       FORMAT (' polar radius of ',i5,'  and ',i5) 
            ii=0
            DO  i = 1, nin 
               IF (isel (i) .eq.1) then 
                  ii = ii + 1 
                  crit (ii) = sqrt (sv (i, j2) **2 + sv (i, j1) **2) 
               ENDIF 
            END DO 
! ************************************
         CASE( -4) axes  !  polar angle, in the current unit for angles 
            WRITE (10, 1413) j1, j2, uniy 
 1413 FORMAT       (' polar angle of ',i5,'  and ',i5, 'unit= ',a3) 
            CALL convun_giff ('RAD', uniy, uconv, urev, error) 
! angle update - sostituire con type angle
            ng = 0 
            crv = 0.d0 
            ns = 0 
            ii=0
            DO 394 i = 1, nin 
               IF (isel (i) .eq.1) then 
                  ns = ns + 1 
                  ii = ii + 1 
                  crr = atan2 (sv (i, j2), sv (i, j1) ) 
                  IF (crr.lt.crv - pig) then 
                     ng = ng + 1 
                  ELSEIF (crr.gt.crv + pig) then 
                     ng = ng - 1 
                  ENDIF 
                  crit (ii) = crr + pig * ng * 2.d0 
                  crv = crr 
               ENDIF 
  394       END DO 
            crit (1:ns) = crit (1:ns) * uconv 
! ************************************
         CASE(-5) axes ! cartesian coordinates from polar ones
            WRITE (10, 1414) j1, j2, uniy 
 1414 FORMAT       (' cartesian coord of ',i5,' as radius '/            &
     &         10x,i5, ' as angle, unit= ',a3)                          
            CALL convun_giff ('RAD', uniy, uconv, urev, error) 
            ns = 0
            ii=0 
            DO  i = 1, nin 
               IF (isel (i) .eq.1) then 
                  ns = ns + 1 
                  ii = ii + 1 
                  crr = sv (i, j2) / uconv 
                  rr = sv (i, j1) 
                  crit (ii) = rr * sin (crr) 
                  x (ii) = rr * cos (crr) 
               ENDIF 
            END DO 
         CASE DEFAULT axes ! unknown negative js
            PRINT *, ' Unknown option '
            GOTO 40
         END SELECT axes
! ***********************************************                             
      ELSEIF (js.eq.0) then 
!  linear combination, with integer arguments, of the coordinates 
         WRITE ( *, 141) nco 
  141 FORMAT   (i4,' integers?') 
         DO  j = 1, nco 
            READ ( *, * ) k (j) 
         END DO
         WRITE ( *, 142) (k (i), i = 1, nco) 
         WRITE (10, 142) (k (i), i = 1, nco) 
  142    FORMAT(' combination,  coeff. '/(10i4)) 
         ii = 0 
         DO  i = 1, nin 
            IF (isel (i) .eq.1) then 
               ii = ii + 1 
               crit(ii)=dot_product(k(1:nco),sv(i,1:nco))
            ENDIF 
         END DO 
! ***********************************************
      ELSEIF (js.le.nco) then 
!  plot time and one coordinate    
         WRITE (10, 1440) js 
 1440    FORMAT(' component no. ',i4) 
         ii = 0 
         DO  i = 1, nin 
            IF (isel (i) .eq.1) then 
               ii = ii + 1 
               crit (ii) = sv (i, js) 
            ENDIF 
         END DO 
      ELSE 
         WRITE ( *, 143) j, nco 
  143    FORMAT(' component ',i4,' not available, maximum ',i4) 
         GOTO 40 
      ENDIF 
! *************************************************************         
! argument has been selected                                            
!                                                                       
!  average and sigma                                                    
      ym = sum(crit(1:nn))/nn
      WRITE ( *, 181) ym 
  181 FORMAT( ' average= ',1p,d20.11) 
      WRITE (10, 181) ym 
      y(1:nn)=crit(1:nn)-ym
      sig=sqrt(dot_product(y(1:nn),y(1:nn))/nn)
      WRITE ( *, 182) sig 
  182 FORMAT(' sigma=',1p,d13.5) 
      WRITE (10, 182) sig 
      RETURN 
      END SUBROUTINE argsel                         
! *************************************************************         
!  {\bf timsel} time interval selection                                 
      SUBROUTINE timsel (nin, tv, x, isel, nn) 
      USE var_precision
      USE prompting
      IMPLICIT real(kind=r_kind) (a - h, o - z)      
! non occorre      DOUBLE PRECISION dta,dtz
      DIMENSION tv (nin), x (nin), isel (nin) 
! chose criterium                                                       
   50 WRITE ( *, 151) nin, tv (1), tv (nin) 
  151 FORMAT(i6,' data points, from t=',f16.6,'  to t=',f16.6)
      iin=0
      CALL prompt (' time int. 0=all 1=by no. point 2=by time ?', &
           iin, 0, 2)                 
      IF (iin.eq.0) then 
         ia = 1 
         iz = nin 
         nn = nin 
         isel (1:nn) = 1 
      ELSEIF (iin.eq.1) then 
         WRITE ( * , * ) ' initial an final data point no. ?' 
         ia=1 !default
         CALL prompt(' initial ',ia, 0, nin)
         iz=nin !default
         CALL prompt(' final ',iz, ia, nin)
         nn = iz - ia + 1  
         isel (1:nin) = 0  
         isel (ia:iz) = 1 
      ELSEIF (iin.eq.2) then 
         tmin=minval(tv(1:nin))
         tmax=maxval(tv(1:nin))
         WRITE ( * , * ) ' initial and final time?'
         ta=tmin !default 
         CALL prompt(' initial ',ta, tmin, tmax)
         tz=tmax !default
         CALL prompt(' final ',tz, ta, tmax)
         nn = 0 
         DO 51 n = 1, nin 
            IF (tv (n) .le.tz.and.tv (n) .ge.ta) then 
               nn = nn + 1 
               isel (n) = 1 
            ELSE 
               isel (n) = 0 
            ENDIF 
   51    END DO 
      ELSE 
         WRITE ( *, 154) iin 
  154 FORMAT   (i4,' option not supported') 
         GOTO 50 
      ENDIF 
      IF (nn.gt.0) then 
         ns = 0 
         DO 522 i = 1, nin 
            IF (isel (i) .eq.1) then 
               ns = ns + 1 
               x (ns) = tv (i) 
            ENDIF 
  522    END DO 
      ELSE 
         WRITE ( * , * ) ' no interval selected' 
         GOTO 50 
      ENDIF 
!                                                                       
      WRITE (10, 1520) nn, x (1), x (nn) 
 1520 FORMAT(i6,' data points, time from ',f15.5,'  to  ',f15.5) 
      RETURN 
      END SUBROUTINE timsel                         
! ************************************************************          
!  {\bf traori} choice between original and transformed data            
      SUBROUTINE traori (iav, crit, y, s, nn) 
      USE var_precision
      USE prompting
      IMPLICIT real(kind=r_kind) (a - h, o - z) 
      DIMENSION crit (nn), y (nn), s (nn) 
      iav=0 !default: original data
      CALL prompt(' 1=transformed data 0=original ?', iav, 0, 1) 
      WRITE (10, 601) iav 
  601 FORMAT(' iav=',i4) 
      IF (iav.eq.0) then 
        s (1:nn) = crit (1:nn) 
      ELSEIF (iav.eq.1) then 
        s (1:nn) = y (1:nn) 
      ENDIF 
      RETURN 
      END SUBROUTINE traori                         



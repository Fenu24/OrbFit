! ==============MODULE time_scales===================                   
! TIME AND TIMESCALES: Carpino 1996-1998                                
! CONTAINS                                                              
! SUBROUTINES                                                           
!                                                                       
! *ch2tim	translation of a character string into time value             
! *intmon	transform a 3-letter code of a month into integer             
! tjm1		computation of Modified Julian Date   
! mjddat        transformation of Modified Julian Date to calendar date 
! timnf		reduction of time (MJD+sec) to normal form                     
! itaiut        difference DAT = TAI - UTC as a function of UTC         
! chktsc	check existence of a time scale                                
! cnvtim	conversion between different time scales                       
! deltt		difference DT = ET - UT                                        
! chmon		transforms integer month (1-12) into a 3-letter code           
! bessep	Besselian epoch as a function of MJD                           
! chmo2i	transforms 3-character month names into integer                
! julian        calendar to julian date                                 
!                                                                       
! times         barycentric to terrestrial time (for rrdot only)        
! calendwri     write calendar date 
!                                                                       
!                                                                       
! Copyright (C) 1997 by Mario Carpino (carpino@brera.mi.astro.it)       
! Version: March 11, 1997                                               
! --------------------------------------------------------------------- 
!                                                                       
!  *****************************************************************    
!  *                                                               *    
!  *                         C H 2 T I M                           *    
!  *                                                               *    
!  *     Translation of a character string into time value         *    
!  *                                                               *    
!  *****************************************************************    
!                                                                       
! INPUT:    STRING    -  Time string                                    
!                                                                       
! OUTPUT:   MJD       -  Modified Julian Date (integer part)            
!           SEC       -  Seconds within the day                         
!           SCALE     -  Time scale (UTC/UT1/TAI/TDT/ET/GPS)            
!           ERROR     -  Error flag                                     
!                                                                       
      SUBROUTINE ch2tim(string,mjd,sec,scale,error) 
      IMPLICIT NONE 
                                                                        
      CHARACTER*(*) string,scale 
      INTEGER mjd 
      DOUBLE PRECISION sec 
      LOGICAL error 
                                                                        
      CHARACTER*80 c,c1,rest1,date,tmp,hhmmss 
      CHARACTER*10 format 
      INTEGER l,posb,intval,year,lt,month,day,hh,mm 
      DOUBLE PRECISION secval,reaval,ss 
      LOGICAL tmperr 
                                                                        
      INTEGER lench,intmon 
      LOGICAL isnum,islett 
      DOUBLE PRECISION tjm1 
      EXTERNAL lench,intmon,isnum,islett,tjm1 
                                                                        
      error=.true. 
      IF(lench(string).GT.80) RETURN 
      c=string 
      CALL norstr(c,l) 
      posb=INDEX(c,' ') 
      IF(posb.LE.1) RETURN 
                                                                        
! Time format                                                           
      CALL strcnt(c,format,rest1,tmperr) 
      IF(tmperr) RETURN 
      CALL upcase(format) 
      c=rest1 
      posb=INDEX(c,' ') 
                                                                        
! Time value                                                            
      IF(format.EQ.'JD' .OR. format.EQ.'MJD') THEN 
! Format A                                                              
          IF(INDEX(c(1:posb),'.').EQ.0) THEN 
              READ(c(1:posb-1),*,ERR=1) intval 
              c1=c(posb+1:) 
              c=c1 
              posb=INDEX(c,' ') 
              IF(posb.LE.0)                                             &
     &            STOP '**** ch2tim: internal error (01) ****'          
              IF(posb.EQ.1) RETURN 
              READ(c(1:posb-1),*,ERR=1) secval 
              c1=c(posb+1:) 
              c=c1 
          ELSE 
! Format B                                                              
              READ(c(1:posb-1),*,ERR=1) reaval 
              c1=c(posb+1:) 
              c=c1 
              intval=reaval 
              IF(intval.GT.reaval) intval=intval-1 
              secval=(reaval-intval)*86400.d0 
          END IF 
      ELSEIF(format.EQ.'CAL') THEN 
! Year                                                                  
          date=c(1:posb-1) 
          CALL stspli(date,'/',tmp,tmperr) 
          IF(tmperr) RETURN 
          READ(tmp,*,ERR=1) year 
! Month                                                                 
          CALL stspli(date,'/',tmp,tmperr) 
          lt=lench(tmp) 
          IF(isnum(tmp(1:lt))) THEN 
              READ(tmp,*,ERR=1) month 
          ELSEIF(islett(tmp(1:lt))) THEN 
              month=intmon(tmp) 
              IF(month.LE.0) GOTO 1 
          ELSE 
              GOTO 1 
          END IF 
!  Day                                                                  
          READ(date,*,ERR=1) day 
          mjd=NINT(tjm1(day,month,year,0.d0)) 
          c1=c(posb+1:) 
          c=c1 
          posb=INDEX(c,' ') 
          IF(posb.EQ.1) GOTO 1 
! hh:mm:ss.sss                                                          
          hhmmss=c(1:posb-1) 
          CALL stspli(hhmmss,':',tmp,tmperr) 
          IF(tmperr) RETURN 
          READ(tmp,*,ERR=1) hh 
          CALL stspli(hhmmss,':',tmp,tmperr) 
          IF(tmperr) RETURN 
          READ(tmp,*,ERR=1) mm 
          READ(hhmmss,*,ERR=1) ss 
          sec=hh*3600+mm*60+ss 
          c1=c(posb+1:) 
          c=c1 
      ELSE 
          RETURN 
      END IF 
                                                                        
      IF(format.EQ.'MJD') THEN 
          mjd=intval 
          sec=secval 
      ELSEIF(format.EQ.'JD') THEN 
          mjd=intval-2400001 
          sec=secval+43200.d0 
      ELSEIF(format.EQ.'CAL') THEN 
          CONTINUE 
      ELSE 
          STOP '**** ch2tim: internal error (02) ****' 
      END IF 
                                                                        
! Timescale                                                             
      CALL strcnt(c,scale,rest1,tmperr) 
      IF(tmperr) RETURN 
      IF(lench(rest1).GT.0) RETURN 
      CALL chktsc(scale,tmperr) 
      IF(tmperr) RETURN 
      CALL timnf(mjd,sec,scale) 
      error=.false. 
      RETURN 
                                                                        
! Error termination                                                     
    1 CONTINUE 
      error=.true. 
                                                                        
      END                                           
! Copyright (C) 1996 by Mario Carpino (carpino@brera.mi.astro.it)       
! Version: August 26, 1996                                              
! --------------------------------------------------------------------- 
!                                                                       
!  *****************************************************************    
!  *                                                               *    
!  *                         I N T M O N                           *    
!  *                                                               *    
!  *   Transforms a 3-letter code of a month into integer value    *    
!  *                      (0 = input error)                        *    
!  *                                                               *    
!  *****************************************************************    
!                                                                       
! INPUT:    CHMON     -  3-letter month                                 
!                                                                       
      INTEGER FUNCTION intmon(chmon) 
      IMPLICIT NONE 
                                                                        
      CHARACTER*(*) chmon 
                                                                        
      CHARACTER*3 c 
                                                                        
      c=chmon 
      CALL upcase(c) 
                                                                        
      IF(c.EQ.'JAN' .OR. c.EQ.'GEN') THEN 
          intmon=1 
      ELSEIF(c.EQ.'FEB') THEN 
          intmon=2 
      ELSEIF(c.EQ.'MAR') THEN 
          intmon=3 
      ELSEIF(c.EQ.'APR') THEN 
          intmon=4 
      ELSEIF(c.EQ.'MAY' .OR. c.EQ.'MAG') THEN 
          intmon=5 
      ELSEIF(c.EQ.'JUN' .OR. c.EQ.'GIU') THEN 
          intmon=6 
      ELSEIF(c.EQ.'JUL' .OR. c.EQ.'LUG') THEN 
          intmon=7 
      ELSEIF(c.EQ.'AUG' .OR. c.EQ.'AGO') THEN 
          intmon=8 
      ELSEIF(c.EQ.'SEP' .OR. c.EQ.'SET') THEN 
          intmon=9 
      ELSEIF(c.EQ.'OCT' .OR. c.EQ.'OTT') THEN 
          intmon=10 
      ELSEIF(c.EQ.'NOV') THEN 
          intmon=11 
      ELSEIF(c.EQ.'DEC' .OR. c.EQ.'DIC') THEN 
          intmon=12 
      ELSE 
          intmon=0 
      END IF 
                                                                        
      END                                           
! Copyright (C) 1996 by Mario Carpino (carpino@brera.mi.astro.it)       
! Version: August 9, 1996                                               
! --------------------------------------------------------------------- 
!                                                                       
!  *****************************************************************    
!  *                                                               *    
!  *                           T J M 1                             *    
!  *                                                               *    
!  *            Computation of Modified Julian Date                *    
!  *                                                               *    
!  *****************************************************************    
!                                                                       
! INPUT:    IDAY      -  Day of the month ( 1 <= IDAY <= 31 )           
!           MONTH     -  Month of the year ( 1 <= MONTH <= 12 )         
!           IYEAR     -  Year (e.g.: 1987)                              
!           H         -  Hour of the day ( 0. <= H < 24. )              
!                                                                       
! OUTPUT:   TJM1      -  Modified Julian Day MJD = JD -2,400,000.5      
!                                                                       
      DOUBLE PRECISION FUNCTION tjm1(iday,month,iyear,h) 
      IMPLICIT NONE 
                                                                        
      DOUBLE PRECISION h 
      INTEGER iday,month,iyear,iy,im,ib,k1,k2,k3 
                                                                        
      IF(month.LE.2) THEN 
          iy=iyear-1 
          im=month+12 
      ELSE 
          iy=iyear 
          im=month 
      END IF 
      IF(iyear.GT.1582) THEN 
          ib=iy/400-iy/100 
      ELSE 
          ib=-2 
          IF(iyear.EQ.1582) THEN 
              IF(month.gt.10) THEN 
                  ib=iy/400-iy/100 
              ELSEIF(month.EQ.10.AND.iday.GE.15) THEN 
                  ib=iy/400-iy/100 
              END IF 
          END IF 
      END IF 
      k1=365.25d0*iy 
      k2=30.6001d0*(im+1) 
      k3=k1+k2+ib-679004+iday 
      tjm1=k3+h/24.d0 
      END                                           
! Author: Mario Carpino (carpino@brera.mi.astro.it)                     
! Version: April 1, 1996                                                
!                                                                       
!  *****************************************************************    
!  *                                                               *    
!  *                          T I M N F                            *    
!  *                                                               *    
!  *              Reduction of time to normal form                 *    
!  *                                                               *    
!  *****************************************************************    
!                                                                       
! INPUT:    MJD       -  Modified Julian Day (integer part)             
!           SEC       -  Seconds within day                             
!           SCALE     -  Time scale                                     
!                                                                       
! OUTPUT:   MJD,SEC are normalized, namely SEC is reduced within        
!           the limits 0 <= SEC < LOD, and MJD is changed accordingly;  
!           LOD (length of the day) is usually 86400 s, but can be      
!           different from that value in the case of UTC, due to        
!           leap seconds                                                
!                                                                       
      SUBROUTINE timnf(mjd,sec,scale) 
      IMPLICIT NONE 
      INTEGER mjd 
      DOUBLE PRECISION sec 
      CHARACTER*3 scale 
                                                                        
! Tuning parameters                                                     
      INTEGER nitmax,nitutc 
      PARAMETER (nitmax=5) 
      PARAMETER (nitutc=20) 
                                                                        
      INTEGER mjd1,nit,idur,isec,k 
      DOUBLE PRECISION sec1,fsec 
      EXTERNAL itaiut 
      INTEGER itaiut 
                                                                        
! Input values are stored for error messages                            
      mjd1=mjd 
      sec1=sec 
                                                                        
! Non-trivial case: UTC (the duration of the day can be different       
! from 86400 s)                                                         
                                                                        
      IF(scale.EQ.'UTC') THEN 
          nit=0 
    1     CONTINUE 
          nit=nit+1 
          IF(nit.GT.nitutc) THEN 
              WRITE(*,100) mjd1,sec1,scale 
              WRITE(*,*)mjd1,sec1
              STOP ' **** timnf: abnormal END ****' 
          END IF 
          IF(sec.LT.0.d0) THEN 
! Duration in seconds of the previous day                               
              idur=86400+itaiut(mjd)-itaiut(mjd-1) 
              sec=sec+idur 
              mjd=mjd-1 
              GOTO 1 
          END IF 
! Decomposition of SEC into integer part (ISEC) + fraction (FSEC),      
! where 0 <= FSEC < 1                                                   
          isec=sec 
          fsec=sec-isec 
! Duration in seconds of today (MJD)                                    
    2     idur=86400+itaiut(mjd+1)-itaiut(mjd) 
! Renormalization of time                                               
          IF(isec.GE.idur) THEN 
              isec=isec-idur 
              mjd=mjd+1 
              GOTO 2 
          END IF 
          sec=isec+fsec 
                                                                        
! Trivial case: the duration of the day is always 86400 s               
! Also this case requires iterations, due to rounding-off problems.     
! EXAMPLE: Let's suppose that the starting values are MJD=48000,        
! SEC=-1.d-14. The result of the first iteration is then:               
!    SEC --> SEC+86400 = 86400 EXACTLY (due to rounding off)            
!    MJD --> MJD-1 = 47999                                              
! Therefore, a second iteration is required, giving:                    
!    SEC --> SEC-86400 = 0                                              
!    MJD --> MJD+1 = 48000                                              
                                                                        
      ELSE 
          nit=0 
    3     nit=nit+1 
          IF(nit.GT.nitmax) THEN 
              WRITE(*,100) mjd1,sec1,scale 
              WRITE(*,*)mjd1,sec1
              STOP ' **** timnf: abnormal END ****' 
          END IF 
          k=sec/86400.d0 
          IF(sec.LT.0.d0) k=k-1 
          IF(k.EQ.0) RETURN 
          mjd=mjd+k 
          sec=sec-k*86400.d0 
          GOTO 3 
      END IF 
                                                                        
  100 FORMAT(' **** timnf: too many iterations ****'/                   &
     &       i7,f20.12,1x,a)                                            
      END                                           
! Copyright (C) 1996-1998 by Mario Carpino (carpino@brera.mi.astro.it)  
! Version: June 9, 1998                                                 
! --------------------------------------------------------------------- 
!                                                                       
!  *****************************************************************    
!  *                                                               *    
!  *                         I T A I U T                           *    
!  *                                                               *    
!  *      Difference DAT = TAI - UTC as a function of UTC          *    
!  *                                                               *    
!  *****************************************************************    
!                                                                       
! INPUT:    MJDC      -  Modified Julian Day (UTC)                      
!                                                                       
! OUTPUT:   Difference DAT = TAI - UTC given as an integer number       
!           of seconds                                                  
!                                                                       
      INTEGER FUNCTION itaiut(mjdc) 
      IMPLICIT NONE 
! ---------------------------------------------------------------------
! Max number of records in TAI-UTC (leap second) file
      INTEGER, PARAMETER :: ntaix=50 
                                                                        
      INTEGER mjdc 
                                                                        
      INTEGER mjdv(ntaix),idv(ntaix) 
      CHARACTER rec*10,file*80 
      INTEGER unit,n,day,month,year,iii,jp,jpm1,lf,j,nwmax,nw 
      LOGICAL first,found,fail1,fail,pedant,fndfil 
                                                                        
      INTEGER lench 
      DOUBLE PRECISION tjm1 
      EXTERNAL lench,tjm1 
                                                                        
      SAVE mjdv,idv,first,n,jp,jpm1,file,pedant,nwmax,nw,lf 
                                                                        
      DATA first/.true./ 
                                                                        
      IF(first) THEN 
          fail=.false. 
          file='TAI-UTC.dat' 
          CALL rdncha('TAI-UTC.','file',file,.false.,fndfil,fail1,fail) 
          lf=lench(file) 
          pedant=.false. 
          CALL rdnlog('TAI-UTC.','pedantic',pedant,.false.,found,       &
     &                fail1,fail)                                       
          nwmax=1 
          nw=0 
          CALL rdnint('TAI-UTC.','n_warn',nwmax,.false.,found,          &
     &                fail1,fail)                                       
          IF(fail) STOP '**** itaiut: abnormal end ****' 
          IF(fndfil) THEN 
              CALL filopn(unit,file,'old') 
          ELSE 
              CALL filopl(unit,file) 
          END IF 
                                                                        
! Skipping header                                                       
    1     READ(unit,100,end=10) rec 
          IF(rec.NE.'----------') GOTO 1 
                                                                        
          n=0 
    2     CONTINUE 
          READ(unit,*,END=3) day,month,year,iii 
          n=n+1 
          IF(n.GT.ntaix) STOP '**** itaiut: n > ntaix ****' 
          idv(n)=iii 
          mjdv(n)=NINT(tjm1(day,month,year,0.d0)) 
          IF(n.GT.1) THEN 
              IF(mjdv(n).LE.mjdv(n-1))                                  &
     &            STOP '**** itaiut: input file is not sorted ****'     
          END IF 
          GOTO 2 
    3     CONTINUE 
          CALL filclo(unit,' ') 
          IF(n.LT.2) GOTO 10 
                                                                        
          jp=2 
          jpm1=jp-1 
                                                                        
          first=.false. 
      END IF 
  100 FORMAT(a) 
                                                                        
! Trying to use previous value                                          
      IF(mjdc.GE.mjdv(jpm1) .AND. mjdc.LT.mjdv(jp)) THEN 
          itaiut=idv(jpm1) 
          RETURN 
      END IF 
                                                                        
! Check limits                                                          
      IF(mjdc.LT.mjdv(1)) THEN 
          IF(pedant) GOTO 20 
          itaiut=idv(1) 
          GOTO 30 
      END IF 
      IF(mjdc.GE.mjdv(n)) THEN 
          IF(pedant) GOTO 20 
          itaiut=idv(n) 
          GOTO 30 
      END IF 
                                                                        
! Selecting new value                                                   
      DO 4 j=2,n 
      IF(mjdc.LT.mjdv(j)) THEN 
          jp=j 
          jpm1=jp-1 
          itaiut=idv(jpm1) 
          RETURN 
      END IF 
    4 END DO 
                                                                        
      STOP '**** itaiut: internal error (01) ****' 
                                                                        
   10 CONTINUE 
      STOP '**** itaiut: input file is empty ****' 
                                                                        
   20 CONTINUE 
      WRITE(*,101) mjdc,file(1:lf),mjdv(1),mjdv(n)-1 
  101 FORMAT(' **** itaiut: mjdc out of range ****'/                    &
     &       '      Cannot compute TAI-UTC at MJD =',i6/                &
     &       '      File "',a,'" starts at MJD =',i6/                   &
     &       '      and ends at MJD =',i6)                              
      STOP '**** itaiut: abnormal end ****' 
                                                                        
   30 CONTINUE 
      nw=nw+1 
      IF(nw.GT.nwmax) RETURN 
      WRITE(*,102) mjdc,file(1:lf),mjdv(1),mjdv(n)-1 
  102 FORMAT(' WARNING (itaiut): MJD =',I6,                             &
     &       ' is outside the range of data'/                           &
     &       ' contained in file "',A,'" (from',I6,' to',I6,')')        
                                                                        
      END                                           
! Copyright (C) 1997 by Mario Carpino (carpino@brera.mi.astro.it)       
! Version: February 24, 1997                                            
!                                                                       
!  *****************************************************************    
!  *                                                               *    
!  *                         C H K T S C                           *    
!  *                                                               *    
!  *             Check existence of a time scale                   *    
!  *                                                               *    
!  *****************************************************************    
!                                                                       
! INPUT:    SCALE     -  Time scale                                     
!                                                                       
! OUTPUT:   ERROR     -  Error flag                                     
!                                                                       
      SUBROUTINE chktsc(scale,error) 
      IMPLICIT NONE 
                                                                        
      CHARACTER*(*) scale 
      LOGICAL error 
                                                                        
      error=.false. 
                                                                        
      IF(scale.EQ.'UTC') RETURN 
      IF(scale.EQ.'ET')  RETURN 
      IF(scale.EQ.'TDT') RETURN 
      IF(scale.EQ.'UT1') RETURN 
      IF(scale.EQ.'TAI') RETURN 
      IF(scale.EQ.'GPS') RETURN 
                                                                        
      error=.true. 
                                                                        
      END                                           
! Author: Mario Carpino (carpino@brera.mi.astro.it)                     
! Version: April 1, 1996                                                
!                                                                       
!  *****************************************************************    
!  *                                                               *    
!  *                         C N V T I M                           *    
!  *                                                               *    
!  *          General purpose time conversion routine              *    
!  *                                                               *    
!  *****************************************************************    
!                                                                       
! INPUT:    MJD1      -  Modified Julian Day (integer part)             
!           SEC1      -  Seconds within day                             
!           SCALE1    -  Input time scale                               
!           SCALE2    -  Output (required) time scale                   
!                                                                       
! OUTPUT:   MJD2      -  Modified Julian Day (integer part)             
!           SEC2      -  Seconds within day                             
!                                                                       
! Supported time scales:                                                
!        UT1                                                            
!        TAI                                                            
!        UTC                                                            
!        TDT = TT = ET                                                  
!                                                                       
      subroutine cnvtim(mjd1,sec1,scale1,mjd2,sec2,scale2) 
      implicit none 
                                                                        
      integer mjd1,mjd2 
      double precision sec1,sec2 
      character*3 scale1,scale2,eqsc,eqsc2 
                                                                        
      external itaiut,deltt 
      integer loops,nit,mjd2r,mjdt,itaiut 
      character*3 scale 
      double precision tjm2,dt,tjmt,sec2r,err,sect,dat,deltt 
                                                                        
! ET-TAI in seconds                                                     
      double precision etmtai 
      parameter (etmtai=32.184d0) 
! Max number of iterations                                              
      integer nitmax 
      parameter (nitmax=5) 
! Error limit for iterations                                            
      double precision epst 
      parameter (epst=1.d-10) 
                                                                        
      loops=0 
                                                                        
      mjd2=mjd1 
      sec2=sec1 
! Current timescale (in which mjd2,sec2 are given)                      
      scale=scale1 
                                                                        
! Equivalent timescales (de-aliasing non-standard names)                
      eqsc=scale 
      if(eqsc.eq.'TT'.or.eqsc.eq.'ET') eqsc='TDT' 
      eqsc2=scale2 
      if(eqsc2.eq.'TT'.or.eqsc2.eq.'ET') eqsc2='TDT' 
                                                                        
! Check on timescales                                                   
      if(eqsc.ne.'UT1'.and.eqsc.ne.'TDT'.and.eqsc.ne.'TAI'.and.         &
     &        eqsc.ne.'UTC') then                                       
          write(*,103) scale1 
          stop ' **** cnvtim: abnormal end ****' 
      end if 
      if(eqsc2.ne.'UT1'.and.eqsc2.ne.'TDT'.and.eqsc2.ne.'TAI'.and.      &
     &        eqsc2.ne.'UTC') then                                      
          write(*,103) scale2 
          stop ' **** cnvtim: abnormal end ****' 
      end if 
  103 format(' **** cnvtim: unsupported timescale "',a,'" ****') 
                                                                        
    1 continue 
      call timnf(mjd2,sec2,eqsc) 
                                                                        
! Required timescale has been reached                                   
      if(eqsc.eq.eqsc2) return 
                                                                        
! Check on infinite loops                                               
      if(loops.gt.6) stop ' **** cnvtim: too many loops ****' 
      loops=loops+1 
                                                                        
! Transformations are performed according to the following path:        
!                                                                       
!                 UT1 -- TDT -- TAI -- UTC                              
!                                                                       
      if(eqsc.eq.'UT1') then 
                                                                        
! Conversion UT1 --> TDT                                                
          tjm2=mjd2+sec2/86400.d0 
          dt=deltt(tjm2) 
          sec2=sec2+dt 
          eqsc='TDT' 
                                                                        
      elseif(eqsc.eq.'TDT')then 
                                                                        
          if(eqsc2.eq.'UT1')then 
                                                                        
! Conversion TDT --> UT1 (iterative method)                             
!   a) computation of DT = TDT - UT1 using (mjd2,sec2) (TDT) as an      
!      approximate value of UT1                                         
              tjm2=mjd2+sec2/86400.d0 
              dt=deltt(tjm2) 
!   b) subtract DT from (mjd2,sec2), finding a first approximation      
!      for UT1                                                          
              mjdt=mjd2 
              sect=sec2-dt 
!   start iterations                                                    
              nit=0 
    3         call timnf(mjdt,sect,'UT1') 
              nit=nit+1 
              if(nit.gt.nitmax)then 
                  write(*,100) mjd1,sec1,scale1 
                  stop ' **** cnvtim: abnormal end ****' 
              endif 
!   c) try to find the starting value of TDT from the approximate       
!      value of UT1                                                     
              tjmt=mjdt+sect/86400.d0 
              dt=deltt(tjmt) 
              mjd2r=mjdt 
              sec2r=sect+dt 
              call timnf(mjd2r,sec2r,'TDT') 
!   d) computation of error and correction of the approximate value     
              err=(mjd2r-mjd2)*86400.d0+sec2r-sec2 
              if(abs(err).gt.epst)then 
                  sect=sect-err 
                  goto 3 
              endif 
              mjd2=mjdt 
              sec2=sect 
              eqsc='UT1' 
                                                                        
          else 
                                                                        
! Conversion TDT --> TAI                                                
              sec2=sec2-etmtai 
              eqsc='TAI' 
          endif 
                                                                        
      elseif(eqsc.eq.'TAI')then 
                                                                        
          if(eqsc2.eq.'UTC')then 
                                                                        
! Conversion TAI --> UTC (iterative method)                             
!   a) computation of DAT = TAI - UTC using (mjd2,sec2) (TAI) as an     
!      approximate value of UTC                                         
              mjdt=mjd2 
              sect=sec2 
              dat=itaiut(mjdt) 
!   b) subtract DAT from (mjd2,sec2), finding a first approximation     
!      for UTC                                                          
              sect=sec2-dat 
!   start iterations                                                    
              nit=0 
    2         call timnf(mjdt,sect,'UTC') 
              nit=nit+1 
              if(nit.gt.nitmax)then 
                  write(*,101) mjd1,sec1,scale1 
                  stop ' **** cnvtim: abnormal end ****' 
              end if 
!   c) try to find the starting value of TAI from the approximate       
!      value of UTC                                                     
              mjd2r=mjdt 
              sec2r=sect+itaiut(mjdt) 
              call timnf(mjd2r,sec2r,'TAI') 
!   d) computation of error and correction of the approximate value     
              err=(mjd2r-mjd2)*86400.d0+sec2r-sec2                      &
     &               +itaiut(mjd2r)-itaiut(mjd2)                        
              if(abs(err).gt.epst)then 
                  sect=sect-err 
                  goto 2 
              endif 
              mjd2=mjdt 
              sec2=sect 
              eqsc='UTC' 
          else 
                                                                        
! Conversion TAI --> TDT                                                
              sec2=sec2+etmtai 
              eqsc='TDT' 
          endif 
                                                                        
      elseif(eqsc.eq.'UTC')then 
                                                                        
! Conversione UTC --> TAI                                               
          sec2=sec2+itaiut(mjd2) 
          eqsc='TAI' 
                                                                        
      else 
          write(*,102) eqsc 
          stop ' **** cnvtim: abnormal end ****' 
      endif 
                                                                        
      goto 1 
                                                                        
  100 format(' **** cnvtim: nit > nitmax (UT1 --> TDT) ****'/           &
     &       ' mjd =',i6,' sec =',f12.4,' scale = ',a)                  
  101 format(' **** cnvtim: nit > nitmax (TAI --> UTC) ****'/           &
     &       ' mjd =',i6,' sec =',f12.4,' scale = ',a)                  
  102 format(' **** cnvtim: scale "',a,'" not supported ****') 
      END                                           
! Copyright (C) 1997-1998 by Mario Carpino (carpino@brera.mi.astro.it)  
! Version: June 10, 1998                                                
! --------------------------------------------------------------------- 
!                                                                       
!  ***************************************************************      
!  *                                                             *      
!  *                          D E L T T                          *      
!  *                                                             *      
!  *                   Difference DT = ET - UT                   *      
!  *                                                             *      
!  ***************************************************************      
!                                                                       
! INPUT:    TJM       -  Modified Julian Day (UT1)                      
!                                                                       
! OUTPUT:   DELTT     -  DT = ET - UT1 (in seconds)                     
!                                                                       
      DOUBLE PRECISION FUNCTION deltt(tjm) 
      IMPLICIT NONE 
                                                                        
      DOUBLE PRECISION tjm 
                                                                        
! Max number of entries in ET-UT table                                  
      INTEGER, PARAMETER :: nx=2000 
                                                                        
      DOUBLE PRECISION dt,c1,c2 
      DOUBLE PRECISION tv(nx),dtv(nx) 
      INTEGER day,month,year,n,ipos,unit,lf,nwmax,nw 
      CHARACTER record*10,file*80 
      LOGICAL first,found,fail1,fail,pedant,fndfil 
                                                                        
      DOUBLE PRECISION tjm1 
      INTEGER lench 
      EXTERNAL tjm1,lench 
                                                                        
      INCLUDE 'timespan.h90' 
                                                                        
      SAVE first,n,tv,dtv,ipos,file,pedant,nwmax,nw,lf 
                                                                        
      DATA first/.true./ 
                                                                        
! Load table of ET-UT1 as a function of UT1                             
      IF(first) THEN 
          first=.false. 
          fail=.false. 
          file='ET-UT.dat' 
          CALL rdncha('ET-UT.','file',file,.false.,fndfil,fail1,fail) 
          lf=lench(file) 
          pedant=.false. 
          CALL rdnlog('ET-UT.','pedantic',pedant,.false.,found,         &
     &                fail1,fail)                                       
          nwmax=1 
          nw=0 
          CALL rdnint('ET-UT.','n_warn',nwmax,.false.,found,            &
     &                fail1,fail)                                       
          IF(fail) STOP '**** deltt: abnormal end ****' 
          IF(fndfil) THEN 
              CALL filopn(unit,file,'old') 
          ELSE 
              CALL filopl(unit,file) 
          END IF 
                                                                        
    1     CONTINUE 
          READ(unit,100,end=11) record 
          IF(record.NE.'----------') GOTO 1 
                                                                        
          n=0 
    2     CONTINUE 
          READ(unit,*,end=3) day,month,year,dt 
          n=n+1 
          IF(n.GT.nx) STOP '**** deltt: n > nx ****' 
          tv(n)=tjm1(day,month,year,0.d0) 
          dtv(n)=dt 
          GOTO 2 
    3     CONTINUE 
          ipos=1 
          CALL filclo(unit,' ') 
          IF(n.LT.2) STOP '**** deltt: n < 2 ****' 
          temut1=tv(1) 
          temut2=tv(n) 
          temute=(.NOT.pedant) 
      END IF 
  100 FORMAT(a) 
                                                                        
! Trying to use previous value                                          
      IF(tjm.GE.tv(ipos).AND.tjm.LE.tv(ipos+1)) GOTO 5 
                                                                        
! Selecting the records of the table before and after the date          
! supplied                                                              
      IF(tjm.LT.tv(1)) THEN 
          IF(pedant) GOTO 12 
          deltt=dtv(1) 
          GOTO 30 
      END IF 
      IF(tjm.GT.tv(n)) THEN 
          IF(pedant) GOTO 12 
          deltt=dtv(n) 
          GOTO 30 
      END IF 
      DO 4 ipos=1,n-1 
          IF(tjm.GE.tv(ipos).AND.tjm.LE.tv(ipos+1)) GOTO 5 
    4 END DO 
      STOP '**** deltt: internal error (01) ****' 
                                                                        
! Interpolation                                                         
    5 CONTINUE 
      c1=(tv(ipos+1)-tjm)/(tv(ipos+1)-tv(ipos)) 
      c2=1-c1 
      deltt=c1*dtv(ipos)+c2*dtv(ipos+1) 
                                                                        
      RETURN 
                                                                        
!  Error termination                                                    
   11 CONTINUE 
      STOP '**** deltt: no DATA found in input file ****' 
                                                                        
   12 CONTINUE 
      WRITE(*,101) tjm,file(1:lf),tv(1),tv(n) 
  101 FORMAT(' **** deltt: TJM out of range ****'/                      &
     &       ' Cannot compute ET-UT at TJM =',F11.3/                    &
     &       ' File "',A,'" starts at TJM =',F11.3/                     &
     &       ' and ends at TJM =',F11.3)                                
      STOP '**** deltt: abnormal end ****' 
                                                                        
   30 CONTINUE 
      nw=nw+1 
      IF(nw.GT.nwmax) RETURN 
      WRITE(*,102) tjm,file(1:lf),tv(1),tv(n) 
  102 FORMAT(' WARNING (deltt): TJM =',F11.3,                           &
     &       ' is outside the range of data'/                           &
     &       ' contained in file "',A,'" (from',F11.3,' to',F11.3,')')  
                                                                        
      END                                           
! Copyright (C) 1996 by Mario Carpino (carpino@brera.mi.astro.it)       
! Version: August 27, 1996                                              
! --------------------------------------------------------------------- 
!                                                                       
!  *****************************************************************    
!  *                                                               *    
!  *                         M J D D A T                           *    
!  *                                                               *    
!  *   Transformation of Modified Julian Date to calendar date     *    
!  *                                                               *    
!  *****************************************************************    
!                                                                       
! INPUT:    TJM       -  Modified Julian Day                            
!                                                                       
! OUTPUT:   IDAY      -  Day of the month  ( 1 <= IDAY <= 31 )          
!           IMONTH    -  Month of the year ( 1 <= IMONTH <= 12 )        
!           IYEAR     -  Year              (e.g., 1987)                 
!           HOUR      -  Hour of the day   ( 0.0 <= HOUR < 24.0 )       
!                                                                       
      SUBROUTINE mjddat(tjm,iday,imonth,iyear,hour) 
      IMPLICIT NONE 
                                                                        
      DOUBLE PRECISION tjm,hour 
      INTEGER iday,imonth,iyear 
                                                                        
      DOUBLE PRECISION a 
      INTEGER ia,ib,ic,id,ie,if 
                                                                        
      a=tjm+2400001.d0 
      ia=a 
      hour=(a-ia)*24.d0 
      IF(ia.LT.2299161) THEN 
          ic=ia+1524 
      ELSE 
          ib=FLOOR((ia-1867216.25d0)/36524.25d0) 
          ic=ia+ib-FLOOR(ib/4.d0)+1525 
      END IF 
      id=FLOOR((ic-122.1d0)/365.25d0) 
      ie=FLOOR(365.25d0*id) 
      if=FLOOR((ic-ie)/30.6001d0) 
      iday=ic-ie-FLOOR(30.6001d0*if) 
      imonth=if-1-12*FLOOR(if/14.d0) 
      iyear=id-4715-FLOOR((7+imonth)/10.d0) 
                                                                        
      END                                           
! Copyright (C) 1997 by Mario Carpino (carpino@brera.mi.astro.it)       
! Version: November 10, 1997                                            
! --------------------------------------------------------------------- 
!                                                                       
!  *****************************************************************    
!  *                                                               *    
!  *                          C H M O N                            *    
!  *                                                               *    
!  *    Transforms integer month (1-12) into a 3-letter code       *    
!  *                                                               *    
!  *****************************************************************    
!                                                                       
! INPUT:    INTMON    -  Integer month (1-12)                           
!                                                                       
      CHARACTER*3 FUNCTION chmon(intmon) 
      IMPLICIT NONE 
                                                                        
      INTEGER intmon 
                                                                        
      CHARACTER*3 c3(12) 
      DATA c3/'Jan','Feb','Mar','Apr','May','Jun',                      &
     &        'Jul','Aug','Sep','Oct','Nov','Dec'/                      
                                                                        
      IF(intmon.GE.1 .AND. intmon.LE.12) THEN 
          chmon=c3(intmon) 
      ELSE 
          chmon='???' 
      END IF 
                                                                        
      END                                           
! Copyright (C) 1999 by Mario Carpino (carpino@brera.mi.astro.it)       
! Version: November 5, 1999                                             
! --------------------------------------------------------------------- 
!                                                                       
! Besselian epoch as a function of Modified Julian Date                 
!                                                                       
      DOUBLE PRECISION FUNCTION bessep(tjme) 
      IMPLICIT NONE 
                                                                        
      DOUBLE PRECISION tjme 
! Copyright (C) 1999 by Mario Carpino (carpino@brera.mi.astro.it)       
! Version: November 11, 1999                                            
! --------------------------------------------------------------------- 
! Besselian year                                                        
!                                                                       
      DOUBLE PRECISION bessyr 
      PARAMETER (bessyr=365.242198781d0) 
                                                                        
      bessep=1900.d0+(tjme-15019.81352d0)/bessyr 
                                                                        
      END                                           
! Copyright (C) 1999 by Mario Carpino (carpino@brera.mi.astro.it)       
! Version: November 4, 1999                                             
! --------------------------------------------------------------------- 
!                                                                       
!  *****************************************************************    
!  *                                                               *    
!  *                         C H M O 2 I                           *    
!  *                                                               *    
!  *    Transforms 3-character month names into integer (1-12)     *    
!  *                   (0 = not a month code)                      *    
!  *                                                               *    
!  *****************************************************************    
!                                                                       
! INPUT:    INTMON    -  Integer month (1-12)                           
!                                                                       
      INTEGER FUNCTION chmo2i(chm) 
      IMPLICIT NONE 
                                                                        
      CHARACTER*(*) chm 
                                                                        
      CHARACTER*3 c3(12) 
      DATA c3/'JAN','FEB','MAR','APR','MAY','JUN',                      &
     &        'JUL','AUG','SEP','OCT','NOV','DEC'/                      
                                                                        
      INTEGER i 
      CHARACTER*3 chm1 
                                                                        
      chm1=chm 
      CALL upcase(chm1) 
                                                                        
      chmo2i=0 
      DO 1 i=1,12 
      IF(chm1.EQ.c3(i)) THEN 
          chmo2i=i 
          RETURN 
      END IF 
    1 END DO 
                                                                        
      END                                           
       subroutine julian(iy,imo,iday,ihr,imin,sec,jd) 
!                                                                       
!      compute julian date using the algorithm given by van flandern    
!      and pulkkinen, ap.j. suppl. vol 41, page 392.                    
!                                                                       
!      input arguments are iy...sec                                     
!                                                                       
!      output is the real julian date: jd                               
!                                                                       
      implicit none 
      double precision jd,sec 
      integer ihr,imin,iy,imo,iday,iaux 
      iaux=-7*(iy+(imo+9)/12)/4-3*((iy+(imo-9)/7)/100+1)/4+275*imo/9 
      jd=iaux+iday+367.d+0*iy+1721028.5d+0                              &
     &     +(ihr+imin/60.d+0+sec/3600.d+0)/24.d+0                       
      return 
      END                                           
!=========================================================              
! this is used only for rrdot                                           
!                                                                       
      SUBROUTINE TIMES (TDBJD,TTJD,SECDIF) 
!                                                                       
!     THIS ROUTINE COMPUTES THE TERRESTRIAL TIME (TT) JULIAN DATE       
!     CORRESPONDING TO A BARYCENTRIC DYNAMICAL TIME (TDB) JULIAN DATE.  
!     EXPRESSIONS USED IN THIS VERSION ARE APPROXIMATIONS RESULTING     
!     IN ACCURACIES OF ABOUT 20 MICROSECONDS.  SEE EXPLANATORY          
!     SUPPLEMENT TO THE ASTRONOMICAL ALMANAC, PP. 42-44 AND 316.        
!                                                                       
!          TDBJD  = TDB JULIAN DATE (IN)                                
!          TTJD   = TT JULIAN DATE (OUT)                                
!          SECDIF = DIFFERENCE TDBJD-TTJD, IN SECONDS (OUT)             
!                                                                       
!                                                                       
      DOUBLE PRECISION TDBJD,TTJD,SECDIF,SECCON,REV,T0,ECC,             &
     &     TDAYS,M,L,LJ,E,DSIN                                          
!                                                                       
      DATA SECCON / 206264.8062470964D0 /,   REV / 1296000.D0 / 
      DATA T0 / 2451545.00000000D0 / 
!     T0 = TDB JULIAN DATE OF EPOCH J2000.0                             
      DATA ECC / 0.01671022D0 / 
!     ECC = ECCENTRICITY OF EARTH-MOON BARYCENTER ORBIT                 
!                                                                       
      TDAYS = TDBJD - T0 
      M  = ( 357.51716D0 + 0.985599987D0 * TDAYS ) * 3600.D0 
      L  = ( 280.46435D0 + 0.985609100D0 * TDAYS ) * 3600.D0 
      LJ = (  34.40438D0 + 0.083086762D0 * TDAYS ) * 3600.D0 
      M  = DMOD (  M, REV ) / SECCON 
      L  = DMOD (  L, REV ) / SECCON 
      LJ = DMOD ( LJ, REV ) / SECCON 
      E  = M + ECC * DSIN ( M ) + 0.5D0 * ECC**2 * DSIN ( 2.D0 * M ) 
      SECDIF =   1.658D-3 * DSIN ( E )                                  &
     &         + 20.73D-6 * DSIN ( L - LJ )                             
      TTJD = TDBJD - SECDIF / 86400.D0 
!                                                                       
      RETURN 
!                                                                       
      END                                           
! ======================================================================
! calendwri                                                             
! composes a calendar date in the appropriate format                    
! for printing                                                          
! ======================================================================
SUBROUTINE calendwri(tcl,calend) 
  IMPLICIT NONE 
  DOUBLE PRECISION :: tcl             ! input date MJD   
  CHARACTER(LEN=14) :: calend         ! output string   
  INTEGER iyear,imonth,iday           ! calendar date variables 
  DOUBLE PRECISION :: hour,h24
! ======================================================================= 
! calendar date                                                         
  call mjddat(tcl,iday,imonth,iyear,hour) 
  h24=hour/24.d0 
  IF(h24.gt.0.999d0)h24=0.999d0 
  write(calend,'(i4,a1,i2.2,a1,i2.2,f4.3)')                         &
     &     iyear,'/',imonth,'/',iday,h24                                
END SUBROUTINE calendwri


SUBROUTINE convert_mjd_cal(t_mjd,t_str)
  !=============================================================================
  DOUBLE PRECISION,  INTENT(IN)  :: t_mjd  ! Time in MJD
  CHARACTER(LEN=16), INTENT(OUT) :: t_str  ! Time in the format YYYY/MM/DD HH:SS
  !=============================================================================
  CHARACTER(LEN=14)  :: t_aux
  INTEGER            :: year, month 
  DOUBLE PRECISION   :: day
  INTEGER            :: day_int, hour, minute                               
  CHARACTER(LEN=2)   :: month_ch,day_ch,hour_ch,minute_ch 
  !=============================================================================
  CALL calendwri(t_mjd,t_aux)
  READ(t_aux(1:14),'(I4,1X,I2,1X,F6.3)') year,month,day
  day_int = FLOOR(day)
  hour    = FLOOR((day-day_int)*24d0)
  minute  = NINT(((day-day_int)*24d0-hour)*60d0)
  IF(month.LT.10)THEN 
     WRITE(month_ch,'(A1,I1)') '0', month
  ELSE
     WRITE(month_ch,'(I2)') month
  ENDIF
  IF(day_int.LT.10)THEN 
     WRITE(day_ch,'(A1,I1)') '0',day_int
  ELSE
     WRITE(day_ch,'(I2)') day_int
  ENDIF
  IF(hour.LT.10)THEN 
     WRITE(hour_ch,'(A1,I1)') '0',hour
  ELSE
     WRITE(hour_ch,'(I2)') hour
  ENDIF
  IF(minute.lt.10)THEN 
     WRITE(minute_ch,'(A1,I1)') '0',minute
  ELSE
     WRITE(minute_ch,'(I2)') minute
  ENDIF
  WRITE(t_str,'(I4,2(A1,A2),1X,A2,A1,A2)') year,'/',month_ch,'/',day_ch,hour_ch,':',minute_ch

END SUBROUTINE convert_mjd_cal

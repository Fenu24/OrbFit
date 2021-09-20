! Copyright (C) 1998 by Steven Chesley (chesley@dm.unipi.it)            
! Version: Dec. 15, 1998                                                
! --------------------------------------------------------------------- 
!                                                                       
!  *****************************************************************    
!  *                                                               *    
!  *                         W R I R W O                           *    
!  *                                                               *    
!  *       Writes a-priori standard deviation of observations      *    
!  *                       and fit residuals                       *    
!  *                                                               *    
!  *****************************************************************    
!                                                                       
! INPUT:    FILE      -  Output file name                               
!           OBJID     -  IAU Identifier for each observation            
!           IOBS      -  Observation type 1=a,d 2=r,rdot 3=r 4=rdot 5=sa
! WARNING: iobs has bene replaced by type, tech!!!
!           TUTM      -  Time (MJD, UTM)                                
!           OBSCOD    -  Observatory code                               
!           ALPHA     -  Right Ascension (radians)                      
!           RMSA      -  A-priori RMS of right ascension (rad)          
!           RESA      -  Residuals in right ascension (rad) (aka O-C)   
!           DELTA     -  Declination (radians)                          
!           RMSD      -  A-priori RMS of declination (rad)              
!           RESD      -  Residuals in declination (rad)                 
!           SMAG      -  Magnitude observations (string)                
!           RMSMAG       A-priori RMS of magnitude                      
!           RESMAG    -  Residuals in magnitude                         
!           RMSH      -  RMS of residuals in magnitude                  
!           SEL       -  Selection indicator (0=don't use; 1=use for fit
!                              2=use for fit & Gauss method)            
!           CHI2      -  CHI**2 value for each observation              
!           N         -  Number of observations                         
!           RMSRES    -  RMS of the residuals                           
!                                                                       
      SUBROUTINE wrirwo(file,objid,type,tech,tutm,obscod,               &
     &     alpha,rmsa,resa,delta,rmsd,resd,                             &
     &     smag,rmsmag,resmag,rmsh,                                     &
     &     sel,chi2,n,rmsres)                                           
      USE station_coordinates 
      USE fund_const
      IMPLICIT NONE 
                                                                        
      CHARACTER*(*) file 
      INTEGER n,obscod(n),sel(n)
!,iobs(n)
      CHARACTER*1 type(n),tech(n) 
      CHARACTER*(*) objid(n) 
      DOUBLE PRECISION tutm(n),alpha(n),delta(n),rmsa(n),rmsd(n) 
      DOUBLE PRECISION resa(n),resd(n),chi2(n),rmsres 
! magnitudes (string), a priori rms'                                    
      CHARACTER*6 smag(n) 
      DOUBLE PRECISION rmsmag(n) 
! fit residuals, rms                                                    
      DOUBLE PRECISION resmag(n),rmsh 
! =================END INTERFACE=============================            
      INCLUDE 'parcmc.h90' 
      INCLUDE 'jplhdr.h90' 
      INCLUDE 'parobx.h90' 
                                                                        
      INTEGER unit,i,is(nobx),ln 
      INTEGER iyear,imonth,iday,ihour,imin,isec,ideg 
      DOUBLE PRECISION day,hour,minu,sec,resnor 
      CHARACTER*19 magstri,tmpstr 
      CHARACTER*30 rastri,rdstri 
      CHARACTER*37 rstri 
      CHARACTER*3 radtyp 
      CHARACTER*1 obstr,signo 
      LOGICAL radar 
      INTEGER truncat,iobcur,iotr,iore 
! fix for exadecimal station code                                       
      CHARACTER*3 obsstr 
                                                                        
      CALL filopn(unit,file,'UNKNOWN') 
                                                                        
! check data set first                                                  
      radar=.false. 
      DO i=1,n 
         IF(type(i).eq.'R'.or.type(i).eq.'V')THEN 
            radar=.true. 
         ELSEIF(type(i).ne.'O'.and.type(i).ne.'S')THEN 
            WRITE(*,*)'wrirwo: obs.type ',type(i),' unknown, rec.no=',i 
            STOP 
         ENDIF 
      ENDDO 
! ========= SORT OBSERVATIONS BY TIME ==============                    
      call heapsort(tutm,n,is) 
! ========= HANDLE ASTROMETRY OBSERVATIONS =========                    
! astrometry header                                                     
      IF(rmsres.GT.0.D0) WRITE(unit,110) comcha,rmsres 
  110 FORMAT(A1,'RMS of orbit residuals = ',F8.3) 
                                                                        
      IF(rmsh.gt.0.d0)WRITE(unit,111) comcha,rmsh 
  111 FORMAT(A1,'RMS of mag residuals = ',F5.2) 
                                                                        
      WRITE(unit,120) comcha 
  120 FORMAT(A1,'++ OBJ ++ OBS +++++ DATE ++++++',                      &
     &     '  +++++ RIGHT ASCENSION ++++++',                            &
     &     '  ++++++ DECLINATION ++++++++',                             &
     &     '  +++++ APP MAG +++++',                                     &
     &     '  ++++ QUALITY +++')                                        
      WRITE(unit,221) comcha 
  221 FORMAT(A1,'+ DESIG + TYP YYYY MM DD.dddddd',                      &
     &     '  HH MM SS.sss   rms     resid',                            &
     &     '  DD MM SS.ss   rms     resid',                             &
     &     '  MAG COL rms   resid',                                     &
     &     '  OBS     CHI  SEL')                                        
! ========= ASTROMETRY LOOP==========================                   
      DO   i=1,n 
      IF(type(is(i)).eq.'O'.or.type(is(i)).eq.'S')THEN 
         obstr=tech(is(i))
! convert time                                                          
         CALL mjddat(tutm(is(i)),iday,imonth,iyear,hour) 
         day=iday+hour/24.d0 
! convert RA                                                            
         CALL sessag(alpha(is(i))*degrad/15.d0,signo,ihour,imin,sec) 
         IF(signo.eq.'-')STOP 'wrirwo error: negative right ascension.' 
! prepare RA string                                                     
         resnor=resa(is(i))*secrad*cos(delta(is(i))) 
         WRITE(tmpstr,FMT='(F6.3)') sec 
         CALL rmsp(tmpstr,ln) 
         IF(ln.lt.6)tmpstr='0'//tmpstr 
         IF(abs(resnor).gt.999.d0)THEN 
            WRITE(rastri,131)ihour,imin,tmpstr,                         &
     &           rmsa(is(i))*secrad,resnor                              
  131       FORMAT(2x,I2.2,1x,I2.2,1x,A6,F7.2,1P,E9.1) 
         ELSE 
            WRITE(rastri,130)ihour,imin,tmpstr,                         &
     &           rmsa(is(i))*secrad,resnor                              
  130       FORMAT(2x,I2.2,1x,I2.2,1x,A6,F7.2,F9.3) 
         ENDIF 
! convert DEC                                                           
         CALL sessag(delta(is(i))*degrad,signo,ideg,imin,sec) 
! prepare DEC string                                                    
         WRITE(tmpstr,FMT='(F5.2)') sec 
         CALL rmsp(tmpstr,ln) 
         IF(ln.lt.5)tmpstr='0'//tmpstr 
         IF(abs(resd(is(i))*secrad).gt.999.d0)THEN 
            WRITE(rdstri,171)signo,ideg,imin,tmpstr,                    &
     &           rmsd(is(i))*secrad,resd(is(i))*secrad                  
  171       FORMAT(1x,A1,I2.2,1x,I2.2,1x,A5,F7.2,1P,E9.1) 
         ELSE 
            WRITE(rdstri,170)signo,ideg,imin,tmpstr,                    &
     &           rmsd(is(i))*secrad,resd(is(i))*secrad                  
  170       FORMAT(1x,A1,I2.2,1x,I2.2,1x,A5,F7.2,F9.3) 
         ENDIF 
! prepare MAG string                                                    
         IF(rmsmag(is(i)).lt.0)THEN 
            WRITE(magstri,121)smag(is(i)) 
  121       FORMAT(a6,13x) 
         ELSEIF(resmag(is(i)).gt.1.d6)THEN 
            WRITE(magstri,122)smag(is(i)),rmsmag(is(i)) 
  122       FORMAT(a6,1x,f5.2,7x) 
         ELSE 
            WRITE(magstri,123)smag(is(i)),rmsmag(is(i)),resmag(is(i)) 
  123       FORMAT(a6,1x,f5.2,1x,f6.2) 
         ENDIF 
! output  astrometry                                                    
         WRITE(tmpstr,FMT='(F9.6)') day 
         CALL rmsp(tmpstr,ln) 
         IF(ln.lt.9)tmpstr='0'//tmpstr 
! fix for exadecimal station code                                       
         CALL codestat(obscod(is(i)),obsstr) 
         WRITE(unit,101) objid(is(i)),obstr,iyear,imonth,tmpstr,        &
     &        rastri,rdstri,magstri,                                    &
     &        obsstr,sqrt(chi2(is(i))),sel(is(i))                       
  101    FORMAT(1x,A9,2x,a1,2x,I4,1x,I2.2,1x,A9,                        &
     &        A30,A29,2x,A19,                                           &
     &        2x,A3,f9.2,2X,I1,3x)                                      
      ENDIF 
      ENDDO 
! ========= HANDLE RADAR OBSERVATIONS ==========                        
      If(.not.radar) GOTO 99 
! radar header                                                          
      WRITE(unit,128) comcha 
  128 FORMAT(A1,'++ OBJ ++ OBS ++++++ DATE +++++++ ',                   &
     &     '++++++++ RADAR RANGE/RANGE RATE +++++++++ ',                &
     &     '++++++ QUALITY +++++')                                      
      WRITE(unit,228) comcha 
  228 FORMAT(A1,'+ DESIG + TYP YYYY MM DD hh:mm:ss ',                   &
     &     'TYP   KM or KM/DAY  a priori rms residual ',                &
     &     'TRX REC     CHI  SEL')                                      
! ========= RADAR LOOP==========================                        
      DO i=1,n 
         IF(type(is(i)).eq.'R'.or.type(is(i)).eq.'V')THEN 
            IF(tech(is(i)).eq.'s')THEN 
! surface bounce                                                        
               obstr='r' 
            ELSE 
! already corrected to center of mass 
               obstr='R' 
            ENDIF 
! convert time                                                          
            CALL mjddat(tutm(is(i)),iday,imonth,iyear,hour) 
! convert hour to 12:12:12                                              
            ihour=truncat(hour,1d-7) 
            minu=(hour-ihour)*60.d0 
            imin=truncat(minu,1d-5) 
            sec=(minu-imin)*60.d0 
            isec=truncat(sec,1d-3) 
            IF(type(is(i)).eq.'R')THEN 
! range observation                                                     
               radtyp='DEL' 
               IF(abs(resa(is(i))*au).gt.999.d0)THEN 
                  WRITE(rstri,143)alpha(is(i))*au,rmsa(is(i))*au,       &
     &                 resa(is(i))*au                                   
               ELSE 
                  WRITE(rstri,141)alpha(is(i))*au,rmsa(is(i))*au,       &
     &                 resa(is(i))*au                                   
               ENDIF 
            ELSEIF(type(is(i)).eq.'V')THEN 
! range-rate observation                                                
               radtyp='DOP' 
               IF(abs(resd(is(i))*au).gt.999.d0)THEN 
                  WRITE(rstri,143)delta(is(i))*au,rmsd(is(i))*au,       &
     &                 resd(is(i))*au                                   
               ELSE 
                  WRITE(rstri,141)delta(is(i))*au,rmsd(is(i))*au,       &
     &                 resd(is(i))*au                                   
               ENDIF 
            ELSE 
               STOP '*** wrirwo.f: internal error(1) ***' 
            ENDIF 
  143       FORMAT(f16.5,1x,f9.5,1x,1p,e10.4) 
  141       FORMAT(f16.5,1x,f9.5,1x,f10.5) 
! find codes of two observatories                                       
            iotr= obscod(is(i))/10000 
            iore= obscod(is(i))-iotr*10000 
! output radar data                                                     
            WRITE(unit,102) objid(is(i)),obstr,iyear,imonth,iday,       &
     &           ihour,imin,isec,radtyp,rstri,                          &
     &           iotr,iore,sqrt(chi2(is(i))),sel(is(i))                 
  102       FORMAT(1x,A9,2x,a1,2x,I4,1x,I2.2,1x,I2.2,1x,                &
     &           I2.2,':',I2.2,':',I2.2,1x,A3,1x,A37,                   &
     &           1X,I3.3,1X,I3.3,f9.2,2X,I1,3x)                         
         ENDIF 
      ENDDO 
! ============================================                          
   99 CALL filclo(unit,' ') 
      END SUBROUTINE wrirwo                                          
 ! Copyright (C) 1998 by Steven Chesley (chesley@dm.unipi.it)            
! Version: Jan. 22, 1998                                                
! Version: September 10, 1999                                           
! --------------------------------------------------------------------- 
!                                                                       
!  *****************************************************************    
!  *                                                               *    
!  *                         R E A R W O                           *    
!  *                                                               *    
!  *                 Reads observations, apriori rms,              *    
!  *                   and post fit redisuals                      *    
!  *                                                               *    
!  *****************************************************************    
!                                                                       
! INPUT:    FILE      -  Input file name                                
!           NLEF      -  max dimension for arrays                       
! OUTPUT:                                                               
!           OBJID     -  IAU Identifier for each observation            
!           IOBS      -  Observation type 1000+x astrometry 2000+x radar
!           TAU       -  Time(MJD, TDT)                                 
!           TUTM      -  Time (MJD, UTM)                                
!           OBSCOD    -  Observatory code                               
!           ALPHA     -  Right Ascension (radians)                      
!           RMSA      -  A-priori RMS of right ascension (rad)          
!           RESA      -  Residuals in right ascension (rad) (aka O-C)   
!           DELTA     -  Declination (radians)                          
!           RMSD      -  A-priori RMS of declination (rad)              
!           RESD      -  Residuals in declination (rad)                 
!           SMAG      -  Magnitude observations (string)                
!           RMSMAG       A-priori RMS of magnitude                      
!           RESMAG    -  Residuals in magnitude                         
!           SEL       -  Selection indicator (0=don't use; 1=use for fit
!                              2=use for fit & Gauss method)            
!           CHI2      -  CHI**2 value for each observation              
!           N         -  Number of observations                         
!                                                                       
      SUBROUTINE rearwo(file,objid,iobs,tau,tutm,obscod,                &
     &     alpha,rmsa,resa,                                             &
     &     delta,rmsd,resd,                                             &
     &     smag,rmsmag,resmag,                                          &
     &     sel,chi2,n,nlef)                                             
      USE fund_const
      USE station_coordinates 
      USE output_control
      IMPLICIT NONE 
! input file, space left in arrays                                      
      CHARACTER*(*) file 
      INTEGER nlef 
! number of obs, observatory, selection flag, obs. type, object identifi
      INTEGER n,obscod(nlef),sel(nlef),iobs(nlef) 
      CHARACTER*(*) objid(nlef) 
! time of observation                                                   
      DOUBLE PRECISION tutm(nlef),tau(nlef) 
! observations and residuals                                            
      DOUBLE PRECISION alpha(nlef),delta(nlef),rmsa(nlef),rmsd(nlef) 
      DOUBLE PRECISION resa(nlef),resd(nlef),chi2(nlef) 
! magnitudes (string), a priori rms'                                    
      CHARACTER*6 smag(nlef) 
      DOUBLE PRECISION rmsmag(nlef) 
! fit residuals, rms                                                    
      DOUBLE PRECISION resmag(nlef) 
! =================END INTERFACE=============================           
      INCLUDE 'parcmc.h90' 
                                                                        
      INTEGER unit,i,iday,month,year,deg,mindec,hr,minra,isec,imin,ihour 
      DOUBLE PRECISION day,secra,secdec 
                                                                        
      CHARACTER*140 rec,tmprec 
      CHARACTER*1 obstyp,sign 
      CHARACTER*3 scale 
      CHARACTER*37 rstri 
      CHARACTER*3 radtyp 
      CHARACTER*9 chistr 
      DOUBLE PRECISION chi,sec,sect 
      INTEGER mjd,mjdt,ll,iotr,iore 
                                                                        
      DOUBLE PRECISION tjm1 
      EXTERNAL tjm1 
! fix for exadecimal station code                                       
      CHARACTER*3 obsstr 
                                                                        
      CALL filopn(unit,file,'UNKNOWN') 
                                                                        
! ========= PROCESS RECORDS SEQUENTIALLY ===========                    
      n=0 
      DO 1  i=1,nlef+20 
         READ(unit,100,END=99,ERR=10)rec 
  100    FORMAT(A) 
! skip comments                                                         
         tmprec=rec 
         CALL rmsp(tmprec,ll) 
         IF(tmprec(1:1).eq.comcha) GOTO 1 
! otherwise get observation code                                        
         obstyp=rec(13:13) 
         IF(obstyp.ne.'R'.and.obstyp.ne.'r'.and.obstyp.ne.'S')THEN 
            n=n+1 
            IF(n.gt.nlef) STOP 'rearwo: nobs > nobx.' 
            iobs(n)=1000+ichar(obstyp) 
! Read astrometry observation                                           
            READ(rec,101,ERR=10) objid(n),year,month,day,               &
     &           hr, minra, secra, rmsa(n),resa(n),                     &
     &           sign,deg,mindec,secdec,rmsd(n),resd(n),                &
     &           smag(n),rmsmag(n),resmag(n),                           &
     &           obsstr,chistr,sel(n)                                   
  101       FORMAT(1x,A9,5x,I4,I3,F10.6,                                &
     &           2x,I2,1x,I2,F7.3,F7.2,G9.3,                            &
     &           1x,A1,I2,1x,I2,F6.2,F7.2,G9.1,2x,                      &
     &           a6,1x,f5.2,1x,f6.2,                                    &
     &           2X,A3,A9,2X,I1,3x)                                     
            READ(chistr,FMT='(F9.2)',ERR=11)chi 
            chi2(n)=chi**2 
! fix for exadecimal station code                                       
            CALL statcode(obsstr,obscod(n)) 
            GOTO 12 
   11       WRITE(ierrou,*)'rearwo: error in chi ',chistr 
            WRITE(*,*)'rearwo: error in chi ',chistr 
            chi2(n)=0.d0 
            numerr=numerr+1 
   12       CONTINUE 
            IF(rec(100:112).eq.'             ') rmsmag(n)=-1.d0 
!            IF(rmsmag(n).eq.0.d0) rmsmag(n)=-1.d0                      
! convert time                                                          
            IF(year.LT.1972) THEN 
               scale='UT1' 
            ELSE 
               scale='UTC' 
            ENDIF 
            iday=day 
            sec=(day-iday)*86400.d0 
            mjd=nint(tjm1(iday,month,year,0.d0)) 
            CALL cnvtim(mjd,sec,scale,mjdt,sect,'TDT') 
            tutm(n)=mjd+sec/86400.d0 
                                      !test                             
            tau(n)=mjdt+sect/86400.d0 
! convert DEC                                                           
            delta(n)=(deg*3600.d0+mindec*60.d0+secdec)/secrad 
            IF(sign.eq.'-')delta(n)=-delta(n) 
            rmsd(n)=rmsd(n)/secrad 
            resd(n)=resd(n)/secrad 
! convert RA (residuals conversion depends upon delta)                  
            alpha(n)=15.d0*(hr*3600.d0+minra*60.d0+secra)/secrad 
            rmsa(n)=rmsa(n)/secrad 
            resa(n)=resa(n)/cos(delta(n))/secrad 
! Radar Observation                                                     
         ELSEIF(obstyp.eq.'R'.or.obstyp.eq.'r')THEN 
            n=n+1 
            IF(n.gt.nlef) STOP 'rearwo: nobs > nobx.' 
! Read radar observation                                                
            READ(rec,102,ERR=10)objid(n),year,month,iday,ihour,         &
     &           imin,isec,radtyp,rstri,                                &
     &           iotr,iore,chi,sel(n)                                   
  102       FORMAT(1x,A9,2x,1x,2x,I4,1x,I2,1x,I2,1x,I2,':',I2,':',I2,1x,&
     &           A3,1x,A37,                                             &
     &           1X,I3.3,1X,I3.3,f9.2,2X,I1,3x)                         
! warning: if a radar station has  a non-numeric code, we are in a mess.
            obscod(n)=iotr*10000+iore 
            chi2(n)=chi**2 
!           WRITE(*,102) objid(n),year,month,iday,ihour,imin,isec,      
!    +           rstri,vstri,                                           
!    +           iotr,iore,chi,sel(n)                                   
! convert time                                                          
            sec=isec+60.d0*imin+3600.d0*ihour 
            mjd=nint(tjm1(iday,month,year,0.d0)) 
            tutm(n)=mjd+sec/86400.d0 
            IF(year.LT.1972) THEN 
               scale='UT1' 
            ELSE 
               scale='UTC' 
            ENDIF 
            CALL cnvtim(mjd,sec,scale,mjdt,sect,'TDT') 
            tau(n)=mjdt+sect/86400.d0 
            IF(radtyp.eq.'DEL')THEN 
! Prcess range string                                                   
               READ(rstri,141,ERR=10)alpha(n),rmsa(n),resa(n) 
  141          FORMAT(f16.5,1x,f9.5,1x,f10.5) 
               alpha(n)=alpha(n)/aukm 
               rmsa(n)=rmsa(n)/aukm 
               resa(n)=resa(n)/aukm 
               iobs(n)=2002 
               delta(n)=0.d0 
               rmsd(n)=-1.d0 
            ELSEIF(radtyp.eq.'DOP')THEN 
! Prcess range rate string                                              
               READ(rstri,141,ERR=10)delta(n),rmsd(n),resd(n) 
               delta(n)=delta(n)/aukm 
               rmsd(n)=rmsd(n)/aukm 
               resd(n)=resd(n)/aukm 
               iobs(n)=2003 
               alpha(n)=0.d0 
               rmsa(n)=-1.d0 
            ELSE 
               STOP '*** rearwo: internal error (1) ***' 
            ENDIF 
! photometry has no meaning                                             
            smag(n)='      ' 
            rmsmag(n)=-1.d0 
! surface bounce correction required                                    
            IF(obstyp.eq.'r')iobs(n)=iobs(n)+100 
         ELSE 
!           Skip record if unknown type                                 
            WRITE(*,*) 'Unknown obs type ',obstyp,' at line ',          &
     &           i,' in ',file                                          
         ENDIF 
! spaghetti code: only go to 10 on error, else skip line 10.            
         GOTO 1 
   10    WRITE(ierrou,*) 'ERROR while reading line ',i,' of ',file 
         WRITE(ierrou,*) 'skipping record: ',rec 
         WRITE(*,*) 'ERROR while reading line ',i,' of ',file 
         WRITE(*,*) 'skipping record: ',rec 
         n=n-1 
         numerr=numerr+1 
    1 ENDDO 
! ============================================                          
   99 CALL filclo(unit,' ') 
      RETURN 
      END SUBROUTINE rearwo                                        

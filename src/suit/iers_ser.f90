! =======MODULE iers_ser===========================                     
! IERS TIME SERIES: Carpino         
! CONTAINS        
! SUBROUTINES              
! delthp	high precision version of deltt     
! diurot	diurnal rotation matrix    
! dut1r		computation of DUT1R = UT1 - UT1R   
! dut1s		computation of DUT1S = UT1 - UT1S   
! equeqd	equation of the equinoxes (with der)         
! gmsnom	nominal Greenwich Mean Sidereal Time (linear with time)        
! gmstd		Greenwich Mean Sidereal Time (with der)      
! ierini	initialization of routine iersts    
! iersts	IERS time series for EOP            
! *isbadr	check records of IERS Bulletin A   
! nutarg	fundamental arguments of the IAU 1980 theory of nutation       
! nutnd		nutation angles (with der)          
! nutwhr	Wahr's nutation series (with der)     
! obleqd	mean obliquity of ecliptic (with der)      
! precd		precession matrix (with der)     
! *rdbula	read IERS Bulletin A  
! rnutd		nutation matrix (with der)       
! rotpv		rotation of position and velocity vectors  
! rotsys	rotation of reference system     
! xypol		coordinates of the terrestrial pole        
! 
! 
!   LINEAR ALGEBRA, but used only in iers_ser.f     
! pd1mat    
! pd2mat    
! rotmt1    
! rotmt2    
    
! HEADERS and MODULES  
!       iers_ser.o: \
!	comlib.h90 \
!	fund_const.o \
!       chebi_pol.o \
!	parlgi.h90   with chebi_pol 
MODULE iers_ser

IMPLICIT NONE

PUBLIC rotpv, rotsys

! data common to the module


! Besselian year
DOUBLE PRECISION, PARAMETER ::  bessyr=365.242198781d0

! ET (TDT) - TAI (s)
DOUBLE PRECISION, PARAMETER :: etmtai=32.184d0

! Copyright (C) 1999 by Mario Carpino (carpino@brera.mi.astro.it)
! Version: November 12, 1999
! ---------------------------------------------------------------------
! Max number of data points in IERS time series
INTEGER, PARAMETER :: niersx=30000
! Max number of separate interpolations (pcwlgi)
INTEGER, PARAMETER :: niplix=5000

! Time of J2000 in MJD (ET) (integer number of days + sec)    
INTEGER, PARAMETER, PUBLIC :: mj2000=51544 
DOUBLE PRECISION, PARAMETER, PUBLIC :: s2000=43200.D0 

! IERS time series
!
! ieidir      -  IERS input directory
! ieidil      -  Length of IERS input directory name
! ieilis      -  IERS input file list
! ieilil      -  Length of IERS input file list name
! isamp       -  Undersampling factor
! tiers       -  Time of data points (MJD, TDT)
! xiers       -  Value of data points:
!xiers(1) = X pole (arcsec)
!xiers(2) = Y pole (arcsec)
!xiers(3) = TDT-UT1 (s)
!xiers(4) = Dpsi (arcsec)
!xiers(5) = Depsilon (arcsec)
! niers       -  Number of data points
! iutsmo      -  Smoothing used for UT1 (xiers(i,3)) :
!         0 = none
!         1 = UT1R
!         2 = UT1S
! cciera      -  Apply consistency correction (flag)
! flcier      -  IERS consistency correction file
! cncor0(5)   -  Consistency correction (constant term)
! cncor1(5)   -  Consistency correction (linear term)
! cnep0       -  "Zero" epoch for consistency correction (Bess. year)
! nlpler      -  NL (interpolation length) for pcwlgi
! nvpler      -  NV (length of empty zone) for pcwlgi
! nspler      -  NS (length of superposition zone) for pcwlgi
! nsmopl      -  Order of smoothing for pcwlgi
! rmserx(5)   -  Max interpolation RMS error allowed for pcwlgi
! extra       -  Allow extrapolation
! blause      -  Use Bulletin A for extension in the future of real data
! blafil      -  Name of input file for Bulletin A
! tint1       -  Starting time for interpolation (end of extrapolation before data start)
! tint2       -  Ending time for interpolation (start of extrapolation after data end)
! dutd        -  Average time derivative of TDT-UT1
! utd1        -  Starting value for extrapolation of TDT-UT1 (tjme<tint1)
! utd2        -  Starting value for extrapolation of TDT-UT1 (tjme>tint2)
! iicier      -  Initialization check
!
CHARACTER*100 ieidir,ieilis,flcier,blafil
INTEGER ieidil,ieilil,isamp,niers,iutsmo,iicier
INTEGER nlpler,nvpler,nspler,nsmopl
DOUBLE PRECISION tiers(niersx),xiers(niersx,5)
DOUBLE PRECISION cncor0(5),cncor1(5),cnep0,rmserx(5),tint1,tint2
DOUBLE PRECISION dutd,utd1,utd2
LOGICAL cciera,extra,blause


CONTAINS
! 
! Copyright (C) 1999 by Mario Carpino (carpino@brera.mi.astro.it)       
! Version: November 11, 1999    
! --------------------------------------------------------------------- 
! 
!  ***************************************************************      
!  * *      
!  *     D E L T H P     *      
!  * *      
!  *        Computation of DT = ET-UT1 as a function of ET       *      
!  *         (with time derivatives)         *      
!  * *      
!  ***************************************************************      
! 
! INPUT:    MJDE      -  Modified Julian Time (ET): integer part        
! SECE      -  Modified Julian Time (ET): seconds within      
!    the day ( 0<= sece < 86400, so that  
!    TJME = mjde + sece/86400)  
! NDER      -  (=0,1,2) required order for derivatives        
! 
! OUTPUT:   DT        -  DT = ET - UT1  (s)         
! DT1       -  dDT/dET        (s/s)      (only if nder >= 1)  
! DT2       -  d^2 DT/dET^2   (s/s^2)    (only if nder  = 2)  
! 
  SUBROUTINE delthp(mjde,sece,dt,dt1,dt2,nder) 
    IMPLICIT NONE 
    INTEGER mjde,nder 
    DOUBLE PRECISION sece,dt,dt1,dt2 
    DOUBLE PRECISION c1,c2,eop(5),eopd(5),eopdd(5),tjme 
    LOGICAL first 
    SAVE first,c1,c2 
    DATA first/.true./ 
    IF(first) THEN 
       first=.false. 
       c1=1.d0/86400.d0 
       c2=c1/86400.d0 
    END IF
    tjme=mjde+sece/86400.d0 
    CALL iersts(tjme,eop,eopd,eopdd,nder) 
    dt=eop(3) 
    dt1=eopd(3)*c1 
    dt2=eopdd(3)*c2 
  END SUBROUTINE delthp
! Copyright (C) 1999 by Mario Carpino (carpino@brera.mi.astro.it)       
! Version: November 9, 1999     
! --------------------------------------------------------------------- 
! 
!  *****************************************************************    
!  *   *    
!  *     D I U R O T       *    
!  *   *    
!  *       Diurnal rotation matrix (with time derivatives)         *    
!  *   *    
!  *   *    
!  *****************************************************************    
! 
! INPUT:    MJDE      -  Modified Julian Time (ET): integer part        
! SECE      -  Modified Julian Time (ET): seconds within      
!    the day (0 <= SECU < 86400)
! NDER      -  Required order for derivatives (0/1/2)         
! 
! OUTPUT:   ROT       -  Diurnal rotation matrix D(t)         
! ROT1      -  First time derivative dD/dt (s^(-1)) 
!    (only if NDER >= 1)        
! ROT2      -  Second time derivative d^2D/dt^2 (s^(-2))      
!    (only if NDER = 2)         
! 
! The diurnal rotation matrix D(t) performs the transformation from     
! the cartesian coordinates  Xtrue referred to the true equator and     
! equinox of date (at time t) to the body-fixed coordinates XBF,        
! referred to the true equator and Greenwich meridian (at the same      
! time t)   
! 
!  XBF = D(t) Xtrue   
! 
  SUBROUTINE diurot(mjde,sece,rot,rot1,rot2,nder) 
    IMPLICIT NONE 
    INTEGER mjde,nder 
    DOUBLE PRECISION sece,rot(3,3),rot1(3,3),rot2(3,3) 
    INTEGER mjdu 
    DOUBLE PRECISION dt,dt1,dt2,tjme,eqeq,eqeq1,eqeq2,secu,utdet 
    DOUBLE PRECISION gmst,gmstd1,gmstd2,gast,gast1,gast2,utddet 
    DOUBLE PRECISION xpol(2),x1pol(2),x2pol(2) 
    DOUBLE PRECISION ra(3,3),rb(3,3),rc(3,3) 
    DOUBLE PRECISION ra1(3,3),rb1(3,3),rc1(3,3) 
    DOUBLE PRECISION ra2(3,3),rb2(3,3),rc2(3,3) 
    IF(nder.LT.0 .OR. nder.GT.2) STOP '**** diurot: nder = ? ****' 
! Computation of DeltaT = ET-UT1
    CALL delthp(mjde,sece,dt,dt1,dt2,nder) 
! Computation of the coordinates of the terrestrial pole Xpol and Ypol  
    tjme=mjde+sece/86400.d0 
    CALL xypol(tjme,xpol,x1pol,x2pol,nder) 
! Computation of the equation of the equinoxes Delta_psi*cos(eps)       
    CALL equeqd(tjme,eqeq,eqeq1,eqeq2,nder) 
! Computation of UT1 = ET-DeltaT
    mjdu=mjde 
    secu=sece-dt 
    CALL timnf(mjdu,secu,'UT1') 
! Computation of Greenwich Mean Sidereal Time Theta_mean      
    CALL gmstd(mjdu,secu,gmst,gmstd1,gmstd2,nder) 
! Computation of Greenwich Apparent Sidereal Time   
! Theta_app = Theta_mean + Delta_psi*cos(eps)       
    gast=gmst+eqeq 
! D(t) = R2(-Xpol) R1(-Ypol) R3(Theta_app)
    CALL rotmt( -xpol(1) , ra , 2) 
    CALL rotmt( -xpol(2) , rb , 1) 
    CALL rotmt(     gast , rc , 3) 
    IF(nder.GE.1) THEN 
!  dUT1/dET = 1 - dDeltaT/dET   
       utdet=1.d0-dt1 
!  dTheta_app/dET =   
!       = dTheta_mean/dET + d[Delta_psi*cos(eps)]/dET         
!       = dTheta_mean/dUT1*dUT1/dET + d[Delta_psi*cos(eps)]/dET         
       gast1=gmstd1*utdet+eqeq1 
       CALL rotmt1( -xpol(1) , ra1 , 2, -x1pol(1)) 
       CALL rotmt1( -xpol(2) , rb1 , 1, -x1pol(2)) 
       CALL rotmt1(     gast , rc1 , 3,     gast1) 
    END IF
    IF(nder.GE.2) THEN 
!  d^2 UT1/dET^2 = -d^2 DeltaT/dET^2      
       utddet=-dt2 
!  d^2 Theta_app/dET^2 =        
!       = d^2 Theta_mean/dET^2 + d^2[(Delta_psi*cos(eps)]/dET^2         
!       = d^2 Theta_mean/dUT1^2 * (dUT1/dET)^2 +    
!         + dTheta_ mean/dUT1 * d^2 UT1/dET^2 +     
!         + d^2 [Delta_psi*cos(eps)]/dET^2
       gast2=gmstd2*(utdet**2)+gmstd1*utddet+eqeq2 
       CALL rotmt2( -xpol(1) , ra2 , 2, -x1pol(1), -x2pol(1)) 
       CALL rotmt2( -xpol(2) , rb2 , 1, -x1pol(2), -x2pol(2)) 
       CALL rotmt2(     gast , rc2 , 3,     gast1,     gast2) 
    END IF
    IF(nder.GE.2)   rot2=rc2 !CALL assmat(rot2,rc2) 
    IF(nder.GE.1) rot1=rc1 ! CALL assmat(rot1,rc1) 
    rot=rc ! CALL assmat(rot,rc) 
    IF(nder.GE.2) CALL pd2mat(rb,rb1,rb2,rot,rot1,rot2) 
    IF(nder.GE.1) CALL pd1mat(rb,rb1,rot,rot1) 
    rot=MATMUL(rb,rot) ! CALL pdmat(rb,rot) 
    IF(nder.GE.2) CALL pd2mat(ra,ra1,ra2,rot,rot1,rot2) 
    IF(nder.GE.1) CALL pd1mat(ra,ra1,rot,rot1) 
    rot=MATMUL(ra,rot) ! CALL pdmat(ra,rot) 
  END SUBROUTINE diurot
! Copyright (C) 1999 by Mario Carpino (carpino@brera.mi.astro.it)       
! Version: November 11, 1999    
! --------------------------------------------------------------------- 
! 
!  *****************************************************************    
!  *   *    
!  *      D U T 1 R        *    
!  *   *    
!  *   Computation of DUT1R = UT1 - UT1R       *    
!  *   *    
!  *****************************************************************    
! 
! INPUT:    TJME      -  Modified Julian Day (TDT)  
! NTDER     -  Required order of time derivatives   
! 
! OUTPUT:   DU        -  DUT1R (s)        
! DUD       -  dUT1R/dt (s/d)   
! DUDD      -  d2UT1R/d2 (s/d**2)         
! 
! UT1R is a regularized form of UT1, which does not include terms       
! due to zonal tides with periods up to 35 days     
! 
  SUBROUTINE dut1r(tjme,du,dud,dudd,ntder) 
    IMPLICIT NONE 
    INTEGER ntder 
    DOUBLE PRECISION tjme,du,dud,dudd 
    INTEGER nterms 
    PARAMETER (nterms=41) 
    INTEGER unit,i,k 
    INTEGER kmul(5,nterms) 
    DOUBLE PRECISION p,cc,a,ad,add,coef(nterms) 
    DOUBLE PRECISION arg(5),argd(5),argdd(5),sina,cosa 
    CHARACTER*10 rec 
    LOGICAL first 
    SAVE first,kmul,coef 
    DATA first/.true./ 
    IF(ntder.LT.0 .OR. ntder.GT.2) STOP '**** dut1r: ntder = ? ****' 
    IF(first) THEN 
       first=.false. 
       CALL filopl(unit,'dut1r.coe') 
1      CONTINUE 
       READ(unit,100) rec 
       IF(rec.NE.'----------') GOTO 1 
       DO 2 k=1,nterms 
          READ(unit,*)(kmul(i,k),i=1,5),p,cc 
          coef(k)=cc*1.D-4 
2      ENDDO 
       READ(unit,*,END=3) p 
       STOP ' **** dut1r: error (01) ****' 
3      CONTINUE 
       CLOSE(unit) 
    END IF
100 FORMAT(a) 
    CALL nutarg(tjme,arg,argd,argdd,ntder) 
    du=0.d0 
    dud=0.d0 
    dudd=0.d0 
    DO 7 k=1,nterms 
       a=0.d0 
       DO 4 i=1,5 
          a=a+kmul(i,k)*arg(i) 
4      END DO
       sina=SIN(a) 
       du=du+coef(k)*sina 
       IF(ntder.GE.1) THEN 
          ad=0.d0 
          DO 5 i=1,5 
             ad=ad+kmul(i,k)*argd(i) 
5         ENDDO
          cosa=COS(a) 
          dud=dud+coef(k)*cosa*ad 
       END IF
       IF(ntder.GE.2) THEN 
          add=0.d0 
          DO 6 i=1,5 
             add=add+kmul(i,k)*argdd(i) 
6         ENDDO
          dudd=dudd+coef(k)*(COS(a)*add-sina*(ad**2)) 
       END IF
7   END DO
  END SUBROUTINE dut1r
! Copyright (C) 1999 by Mario Carpino (carpino@brera.mi.astro.it)       
! Version: November 11, 1999    
! --------------------------------------------------------------------- 
! 
!  *****************************************************************    
!  *   *    
!  *      D U T 1 S        *    
!  *   *    
!  *   Computation of DUT1S = UT1 - UT1S       *    
!  *   *    
!  *****************************************************************    
! 
! INPUT:    TJME      -  Modified Julian Day (TDT)  
! NTDER     -  Required order of time derivatives   
! 
! OUTPUT:   DU        -  DUT1S (s)        
! DUD       -  dUT1S/dt (s/d)   
! DUDD      -  d2UT1S/d2 (s/d**2)         
! 
! UT1S is a regularized form of UT1, which does not include terms       
! due to zonal tides with periods up to 35 days     
! 
      SUBROUTINE dut1s(tjme,du,dud,dudd,ntder) 
      IMPLICIT NONE 
  
      INTEGER ntder 
      DOUBLE PRECISION tjme,du,dud,dudd 
  
      INTEGER nterms 
      PARAMETER (nterms=62) 
  
      INTEGER unit,i,k 
      INTEGER kmul(5,nterms) 
      DOUBLE PRECISION cs,cc,p,arg(5),argd(5),argdd(5) 
      DOUBLE PRECISION a,ad,add,cosa,sina 
      DOUBLE PRECISION coefc(nterms),coefs(nterms) 
      CHARACTER*10 rec 
      LOGICAL first 
      SAVE first,kmul,coefc,coefs 
      DATA first/.true./ 
  
  
  
      IF(ntder.LT.0 .OR. ntder.GT.2) STOP '**** dut1s: ntder = ? ****' 
  
      IF(first) THEN 
first=.false. 
CALL filopl(unit,'dut1s.coe') 
    1     CONTINUE 
READ(unit,100) rec 
IF(rec.NE.'----------') GOTO 1 
DO 2 k=1,nterms 
READ(unit,101)(kmul(i,k),i=1,5),cs,cc 
coefs(k)=cs*1.d-4 
coefc(k)=cc*1.d-4 
    2     CONTINUE 
READ(unit,*,END=3) p 
STOP '**** dut1s: internal error (01) ****' 
    3     CONTINUE 
CLOSE(unit) 
      END IF 
  100 FORMAT(A) 
  101 FORMAT(I3,4I4,11X,F9.2,F7.2) 
  
      CALL nutarg(tjme,arg,argd,argdd,ntder) 
  
      du=0.d0 
      dud=0.d0 
      dudd=0.d0 
  
      DO 7 k=1,nterms 
      a=0.d0 
      DO 4 i=1,5 
      a=a+kmul(i,k)*arg(i) 
    4 END DO 
      cosa=COS(a) 
      sina=SIN(a) 
      du=du+coefc(k)*cosa+coefs(k)*sina 
  
      IF(ntder.GE.1) THEN 
         ad=0.d0 
         DO 5 i=1,5 
            ad=ad+kmul(i,k)*argd(i) 
    5     CONTINUE 
            dud=dud+ad*(-coefc(k)*sina+coefs(k)*cosa) 
      END IF 
      IF(ntder.GE.2) THEN 
         add=0.d0 
         DO 6 i=1,5 
            add=add+kmul(i,k)*argdd(i) 
6           CONTINUE 
            dudd=dudd-coefc(k)*(sina*add+cosa*(ad**2))&
                 &   +coefs(k)*(cosa*add-sina*(ad**2))
      END IF 
    7 END DO 
  
      END SUBROUTINE dut1s  
! Copyright (C) 1999 by Mario Carpino (carpino@brera.mi.astro.it)       
! Version: November 8, 1999     
! --------------------------------------------------------------------- 
! 
!  *****************************************************************    
!  *   *    
!  *      G M S T D        *    
!  *   *    
!  *       Equation of the equinoxes (with time derivatives)       *    
!  *   *    
!  *   *    
!  *****************************************************************    
! 
! INPUT:    TJME      -  Modified Julian Time (ET)  
! NDER      -  Required order for derivatives (0/1/2)         
! 
! OUTPUT:   EQEQ      -  Equation of the equinoxes (difference
!    between apparent sidereal time and mean        
!    sidereal time) in radians  
! EQEQ1     -  dEQEQ/dt (rad/s) (only if NDER >= 1) 
! EQEQ2     -  d^2EQEQ/dt^2 (rad/s^2) (only if NDER = 2)      
! 
      SUBROUTINE equeqd(tjme,eqeq,eqeq1,eqeq2,nder) 
      IMPLICIT NONE 
  
      INTEGER nder 
      DOUBLE PRECISION tjme,eqeq,eqeq1,eqeq2 
  
      DOUBLE PRECISION obl,obl1,obl2,dpsi,deps,dpsi1,deps1,dpsi2,deps2 
      DOUBLE PRECISION cose,sine 
  
      IF(nder.LT.0.OR.nder.GT.2) STOP '**** equeqd: nder = ? ****' 
      CALL obleqd(tjme,obl,obl1,obl2,nder) 
      CALL nutnd(tjme,dpsi,deps,dpsi1,deps1,dpsi2,deps2,nder) 
  
      cose=COS(obl) 
      eqeq=dpsi*cose 
      IF(nder.GE.1) THEN 
         sine=SIN(obl) 
         eqeq1=dpsi1*cose-dpsi*sine*obl1 
      END IF 
  
      IF(nder.GE.2) THEN 
         eqeq2=dpsi2*cose-2.d0*dpsi1*sine*obl1     &
     &-dpsi1*cose*obl1**2-dpsi*sine*obl2  
      END IF 
  
      END SUBROUTINE equeqd  
! Copyright (C) 1999 by Mario Carpino (carpino@brera.mi.astro.it)       
! Version: November 8, 1999     
! --------------------------------------------------------------------- 
! 
!  *****************************************************************    
!  *   *    
!  *     G M S N O M       *    
!  *   *    
!  * Nominal (strictly linear with time) Greenwich Mean Sid. Time  *    
!  * as a function of ET   *    
!  *   *    
!  *   *    
!  *****************************************************************    
! 
! INPUT:    MJD       -  Modified Julian Time (ET): integer part        
! SEC       -  Modified Julian Time (ET): seconds within      
!    the day (0 <= SECU < 86400)
! 
! OUTPUT:   GMST      -  Greenwich Mean Sidereal Time referred
!    to the mean equinox of date (rad)    
! GMSTD1    -  dGMST/dUT1 (rad/s)         
! 
      SUBROUTINE gmsnom(mjd,sec,gmst,gmstd1) 
      IMPLICIT NONE 
  
      INTEGER mjd 
      DOUBLE PRECISION sec,gmst,gmstd1 
  
      INTEGER mjd0,icheck 
      DOUBLE PRECISION gms0,gms1,sec0 
  
      DOUBLE PRECISION dt 
  
      IF(icheck.NE.33) STOP '**** gmsnom: COMMON not initialized ****' 
  
! GMSNOM(t) is linear with time (ET) by definition: 
! GMSNOM(t) = GMSNOM(t_0) + GMSTD1*(t-t_0)
      dt=(mjd-mjd0)*86400.d0+(sec-sec0) 
      gmst=gms0+gms1*dt 
      gmstd1=gms1 
  
      END SUBROUTINE gmsnom  
! Copyright (C) 1999 by Mario Carpino (carpino@brera.mi.astro.it)       
! Version: November 8, 1999     
! --------------------------------------------------------------------- 
! 
!  *****************************************************************    
!  *   *    
!  *      G M S T D        *    
!  *   *    
!  *       Greenwich Mean Sidereal Time        *    
!  *as a function of UT1 (with time derivatives)         *    
!  *   *    
!  *   *    
!  *****************************************************************    
! 
! INPUT:    MJDU      -  Modified Julian Time (UT1): integer part       
! SECU      -  Modified Julian Time (UT1): seconds within     
!    the day (0 <= SECU < 86400)
! NDER      -  Required order for derivatives (0/1/2)         
! 
! OUTPUT:   GMST      -  Greenwich Mean Sidereal Time referred
!    to the mean equinox of date (rad)    
! GMSTD1    -  dGMST/dUT1 (rad/s) (only if NDER >= 1)         
! GMSTD2    -  d^2GMST/dUT1^2 (rad/s^2) (only if NDER = 2)    
! 
      SUBROUTINE gmstd(mjdu,secu,gmst,gmstd1,gmstd2,nder)
      USE fund_const 
      IMPLICIT NONE 
  
      INTEGER mjdu,nder 
      DOUBLE PRECISION secu,gmst,gmstd1,gmstd2 
  
      DOUBLE PRECISION s2r,s2c,s2c2,a,b,c,d,fday,tu,gmst0 
      DOUBLE PRECISION ssum,dsum,ddsum,c0,c1,c2,c3,t 
      LOGICAL first 
      SAVE first,s2r,s2c,s2c2,a,b,c,d 
      DATA first/.true./ 
  
! Computation of a 3rd degree polynomial with derivatives     
      ssum  (c0,c1,c2,c3,t) = ((c3*t+c2)*t+c1)*t+c0 
      dsum (c0,c1,c2,c3,t) = (3.d0*c3*t+2.d0*c2)*t+c1 
      ddsum(c0,c1,c2,c3,t) = 6.d0*c3*t+2.d0*c2 
  
      IF(nder.LT.0.OR.nder.GT.2) STOP '**** gmstd: nder = ? ****' 
  
      IF(first) THEN 
         first=.false. 
         s2r = pig/(12*3600) 
         s2c = 1.d0/(36525.d0*86400.d0) 
         s2c2= s2c**2 
! Polynomial representation of GMST (see Supplement to the    
! Astronomical Almanac, 1984, S15)        
         a  =    24110.54841d0  * s2r 
         b  =  8640184.812866d0 * s2r 
         c  =        0.093104d0 * s2r 
         d  =       -6.2d-6     * s2r 
      END IF 
  
! Fraction of the day elapsed from 0h UT1 to SECU   
      fday=secu/86400.d0 
  
! Number of centuries elapsed from 2000 January 1, 12h UT1    
! to 0h UT1 of the current day (TJM=MJDU) 
      tu=(mjdu-51544.5d0)/36525.d0 
  
! GMST at 0h UT1 (TJM=MJDU)     
      gmst0=ssum(a,b,c,d,tu) 
  
! dGMST/dUT1 is computed at the midpoint between MJDU and MJDU+FDAY     
      tu=(mjdu+fday/2.d0-51544.5d0)/36525.d0 
      gmstd1=dsum(a,b,c,d,tu)*s2c+dpig/86400.d0 
  
! GMST at TJM=MJDU+FDAY         
      gmst=gmst0+gmstd1*secu 
      gmst=MOD(gmst,dpig) 
      IF(gmst.LT.0.d0) gmst=gmst+dpig 
  
! Derivatives are computed again at TJM=MJDU+FDAY   
! only when requested 
      tu=(mjdu+fday-51544.5d0)/36525.d0 
      IF(nder.GE.1) gmstd1=dsum(a,b,c,d,tu)*s2c+dpig/86400.d0 
      IF(nder.GE.2) gmstd2=ddsum(a,b,c,d,tu)*s2c2 
  
      END SUBROUTINE gmstd  
! Copyright (C) 1999 by Mario Carpino (carpino@brera.mi.astro.it)       
! Version: November 17, 1999    
! --------------------------------------------------------------------- 
! 
!  *****************************************************************    
!  *   *    
!  *     I E R I N I       *    
!  *   *    
!  *   Initialization of routine IERSTS        *    
!  *   *    
!  *****************************************************************    
! 
      SUBROUTINE ierini 
      USE chebi_pol
      IMPLICIT NONE 
  
  
! NEEDED common blocks:         
      INCLUDE 'comlib.h90' 
  
      DOUBLE PRECISION x,y,dt1,dlod,dpsi,deps,dt,tjme,du,dud,dudd,dtsec 
      DOUBLE PRECISION t0,tn,c0,c1,c2 
      INTEGER unilis,uniier,unicie,year,nr,lf,lr,day,mjd,mjdp,mjdc 
      INTEGER month,npt,i,k,itint1,itint2,niers1,k1 
      INTEGER nu,na,ni,ne1,ne2,mjd1,mjd2,kamin 
      CHARACTER ierfil*80,file*150,rec*80,cm3*3,cval*50 
      LOGICAL found,fail1,fail,first 
      INTEGER lench,itaiut 
!     INTEGER chmo2i 
      DOUBLE PRECISION tjm1 
      EXTERNAL lench,itaiut,tjm1 
!     EXTERNAL chmo2i
  
      IF(iiclib.NE.36) STOP '**** ierini: internal error (01) ****' 
  
      npt=0 
      fail=.false. 
  
! Input of options    
! Directory where input IERS EOPC04 files are located         
      ieidir=libdir(1:lenld)//'eopc04/' 
      CALL rdncha('IERS.eopc04.','dir',ieidir,.false.,found,fail1,fail) 
      ieidil=lench(ieidir) 
! File containing the list of IERS files to be used 
      ieilis=libdir(1:lenld)//'eopc04.lis' 
      CALL rdncha('IERS.eopc04.','list',ieilis,.false.,found,fail1,fail) 
      ieilil=lench(ieilis) 
! Smoothing to be used
      iutsmo=0 
      CALL sv2int('IERS.','smoothing',cval,iutsmo,.false.,    &
     &  found,fail1,fail)       
! Undersampling ratio 
      isamp=1 
      CALL rdnint('IERS.','sampling',isamp,.false.,found,fail1,fail) 
! Apply consistency correction  
      cciera=.true. 
      CALL rdnlog('IERS.ccor.','use',cciera,.false.,found,fail1,fail) 
! Name of consistency correction file     
      IF(cciera) THEN 
         flcier='eopc04.cnc' 
CALL rdncha('IERS.ccor.','file',flcier,.false.,     &
     &      found,fail1,fail)   
      ELSE 
         flcier=' ' 
      END IF 
! Interpolation length for pcwlgi         
      nlpler=15 
      CALL rdnint('IERS.','lint',nlpler,.false.,found,fail1,fail) 
! Length of empty zone for pcwlgi         
      nvpler=3 
      CALL rdnint('IERS.','lez',nvpler,.false.,found,fail1,fail) 
! Length of superposition zone for pcwlgi 
      nspler=3 
      CALL rdnint('IERS.','lsz',nspler,.false.,found,fail1,fail) 
! Order of smoothing for pcwlgi 
      nsmopl=6 
      CALL rdnint('IERS.','smord',nsmopl,.false.,found,fail1,fail) 
! Use of Bulletin A (EOP predictions into the future)         
      blause=.false. 
      CALL rdnlog('IERS.bulA.','use',blause,.false.,found,fail1,fail) 
      IF(blause) THEN 
         blafil='bulletinA' 
         CALL rdncha('IERS.bulA.','file',blafil,.false.,     &
     &      found,fail1,fail)   
      END IF 
! Allow extrapolation 
      extra=.true. 
      CALL rdnlog('IERS.','extrapolation',extra,.false.,      &
     &  found,fail1,fail)       
  
! Check on input values         
      nu=nlpler-2*(nvpler+nspler) 
      na=nspler+nu 
      IF(na.LT.1) THEN 
         WRITE(*,220) na 
         STOP '**** ierini: abnormal end ****' 
      END IF 
  220 FORMAT('ERROR: bad parameter definition (IERS.eopc04.)'/&
     &       '       Please supply values such that lint-2*lez-lsz > 0'/&
     &       '       (present value is lint-2*lez-lsz =',I3)  
  
      IF(fail) STOP '*** ierini: abnormal end ****' 
  
! Reserve space for extention records before beginning of data
      IF(extra) THEN 
         kamin=2+nvpler-nu 
         k=kamin/na 
         IF(k*na.LT.kamin) k=k+1 
         IF(k*na.LT.kamin) STOP '**** ierini: internal error (02) ****' 
         IF(k.LT.0) STOP '**** ierini: internal error (03) ****' 
         k=MAX(2,k) 
         ne1=(k+1)*na 
      ELSE 
         ne1=0 
      END IF 
      niers=ne1 
      first=.true. 
  
! Loop on IERS files  
      CALL filopn(unilis,ieilis,'OLD') 
    1 CONTINUE 
      READ(unilis,100,END=10) ierfil 
  100 FORMAT(A) 
      file=ieidir(1:ieidil)//ierfil 
      lf=lench(file) 
  
! Scan IERS file      
      CALL filopn(uniier,file,'OLD') 
      nr=0 
  
! Skip file header and read the year      
    2 CONTINUE 
      nr=nr+1 
      READ(uniier,100,ERR=20) rec 
!     IF(rec(1:10).NE.'  YEAR ==>') GOTO 2 
!     READ(rec(11:),*) year 
!     IF(year.LT.1900 .OR. year.GT.2100) GOTO 20 
!
! Skip file header and records until JAN-18-1972, when TAI-UTC.dat has data available
! and to avoid warning messages from itaiut if we start reading records from JAN-1-1972
!     IF(rec(18:22).NE.'41334') GOTO 2
      IF(rec(15:19).NE.'41334') GOTO 2
  
! Read and store data records   
    3 CONTINUE 
      nr=nr+1 
      READ(uniier,100,ERR=20,END=9) rec 
      lr=lench(rec) 
      IF(lr.LE.0) GOTO 3 
!     READ(rec,101,ERR=20) cm3,day,mjd,x,y,dt1,dlod,dpsi,deps 
! 101 FORMAT(2X,A3,1X,I3,2X,I5,2F9.5,F10.6,2X,F10.6,2X,2F9.5) 
!     READ(rec,101,ERR=20) year,cm3,day,mjd,x,y,dt1,dlod,dpsi,deps
! 101 FORMAT(2X,I4,2X,A3,1X,I3,2X,I5,2F9.6,F10.7,2X,F10.7,2X,2F9.6)
      READ(rec,101,ERR=20) year,month,day,mjd,x,y,dt1,dlod,dpsi,deps
  101 FORMAT(3(I4),I7,2F11.6,2F12.7,2F11.6)
      IF(year.LT.1900 .OR. year.GT.2100) GOTO 20
! Check on MJD value  
!     month=chmo2i(cm3) 
      IF(month.LE.0) GOTO 20 
      mjdc=NINT(tjm1(day,month,year,0.D0)) 
      IF(mjdc.NE.mjd) THEN 
WRITE(*,203) mjdc,mjd 
GOTO 20 
      END IF 
  203 FORMAT('ERROR: inconsistent MJD:',2I8) 
      IF(niers.GT.ne1) THEN 
IF(mjd-mjdp.NE.1) THEN 
    WRITE(*,202) mjd,mjdp 
    GOTO 20 
END IF 
      END IF 
  202 FORMAT('WRONG time sequence:',2I7) 
      mjdp=mjd 
      npt=npt+1 
      IF(MOD(npt,isamp).NE.0) GOTO 3 
! TDT-UT1 = (TDT-TAI) + (TAI-UTC) + (UTC-UT1)       
      dtsec=etmtai+itaiut(mjd) 
      dt=dtsec-dt1 
! Transformation TJM(UTC) -> TJM(TDT)     
! TDT ( = ET ) = (ET-TAI) + (TAI-UTC) + UTC         
      tjme=mjd+dtsec/86400.d0 
  
! Applying required smoothing to TDT-UT1  
      IF(iutsmo.eq.0) THEN 
         CONTINUE 
      ELSEIF(iutsmo.EQ.1) THEN 
         CALL dut1r(tjme,du,dud,dudd,0) 
         dt=dt+du 
      ELSEIF(iutsmo.EQ.2) THEN 
         CALL dut1s(tjme,du,dud,dudd,0) 
         dt=dt+du 
      ELSE 
         STOP '**** ierini: internal error (04) ****' 
      END IF 
! Storing values in arrays      
      niers=niers+1 
      IF(niers.GT.niersx) STOP '**** ierini: niers > niersx ****' 
      IF(first) THEN 
         first=.false. 
         mjd1=mjd 
      END IF 
      mjd2=mjd 
      tiers(niers)=tjme 
      xiers(niers,1)=x 
      xiers(niers,2)=y 
      xiers(niers,3)=dt 
      xiers(niers,4)=dpsi 
      xiers(niers,5)=deps 
      GOTO 3 
  
    9 CONTINUE 
      CALL filclo(uniier,' ') 
  
! End of loop on IERS files     
      GOTO 1 
   10 CONTINUE 
      CALL filclo(unilis,' ') 
  
      IF(niers.LE.0) STOP '**** ierini: empty EOPC04 files ****' 
      niers1=niers 
      WRITE(*,201) niers,mjd1,mjd2 
  201 FORMAT('Reading IERS time series:'/ &
     &       '        EOPC04 data:',I6,' points (from MJD=',  &
     &       I5,' to MJD=',I5,')')        
  
! Reading predictions from IERS Bulletin A
      IF(blause) THEN 
CALL rdbula(blafil,mjd2,tiers,xiers,niers,niersx,isamp,npt,   &
     &      iutsmo)   
      END IF 
  
! Number of extension records after end of data     
      IF(extra) THEN 
         ni=(niers-ne1-nlpler)/na 
         k=nlpler-nspler-nvpler+(ni-1)*na 
         ne2=mjd2-mjd1-k+2*na 
         ne2=MAX(2,ne2) 
      ELSE 
         ne2=0 
      END IF 
  
! Fill extension records        
      IF(extra) THEN 
         dutd=(xiers(niers,3)-xiers(ne1+1,3))/     &
     &    (tiers(niers)-tiers(ne1+1))     
         itint1=mjd1+(nvpler+nspler-ne1)*isamp 
         dtsec=etmtai+itaiut(itint1) 
         tint1=itint1+dtsec/86400.d0 
         dtsec=etmtai+itaiut(mjd1) 
         t0=mjd1+dtsec/86400.d0 
         utd1=xiers(ne1+1,3)+dutd*(tint1-tiers(ne1+1)) 
         DO 11 i=1,ne1 
            k=ne1+1-i 
            mjd=mjd1-k*isamp 
            dtsec=etmtai+itaiut(mjd) 
            tjme=mjd+dtsec/86400.d0 
            tiers(i)=tjme 
            IF(mjd.LT.itint1) THEN 
               c0=0.d0 
            ELSE 
               tn=(tjme-tint1)/(t0-tint1) 
               CALL smoocn(tn,c0,c1,c2,2) 
            END IF
            DO 13 k=1,5 
               xiers(i,k)=c0*xiers(1,k) 
13          ENDDO
            xiers(i,3)=utd1+dutd*(tjme-tint1) 
11       ENDDO
         IF(niers+ne2.GT.niersx)         &
              &        STOP '**** ierini: niers > niersx ****'         
         ni=(niers+ne2-ne1-nlpler)/na 
         k=nlpler-nspler-nvpler+(ni-1)*na 
         itint2=mjd1+k*isamp 
         dtsec=etmtai+itaiut(itint2) 
         tint2=itint2+dtsec/86400.d0 
         dtsec=etmtai+itaiut(mjd2) 
         t0=mjd2+dtsec/86400.d0 
         utd2=xiers(niers,3)+dutd*(tint2-tiers(niers)) 
         DO 12 i=1,ne2 
            mjd=mjd2+i*isamp 
            dtsec=etmtai+itaiut(mjd) 
            tjme=mjd+dtsec/86400.d0 
            k=niers+i 
            k1=niers1+i 
            IF(mjd.GT.itint2) THEN 
               c0=0.d0 
            ELSE 
               tn=(tjme-tint2)/(t0-tint2) 
               CALL smoocn(tn,c0,c1,c2,2) 
            END IF
            tiers(k)=tjme 
            xiers(k,1)=c0*xiers(niers,1) 
            xiers(k,2)=c0*xiers(niers,2) 
            xiers(k,3)=utd2+dutd*(tjme-tint2) 
            xiers(k1,4)=c0*xiers(niers1,4) 
            xiers(k1,5)=c0*xiers(niers1,5) 
12       ENDDO
         niers=niers+ne2 
         IF(ni.LT.0) STOP '**** ierini: internal error (05) ****' 
      END IF 
! Read consistency corrections  
      IF(cciera) THEN 
         CALL filopl(unicie,flcier) 
4        READ(unicie,100) rec 
         IF(rec(1:5).NE.'-----') GOTO 4 
         READ(unicie,*) cncor0 
         READ(unicie,*) cncor1 
         READ(unicie,*) cnep0 
         CALL filclo(unicie,' ') 
! Unit conversion     
         cncor0(1)=cncor0(1)/1.d3 
         cncor0(2)=cncor0(2)/1.d3 
         cncor0(3)=cncor0(3)/1.d4 
         cncor0(4)=cncor0(4)/1.d3 
         cncor0(5)=cncor0(5)/1.d3 
         cncor1(1)=cncor1(1)/1.d3 
         cncor1(2)=cncor1(2)/1.d3 
         cncor1(3)=cncor1(3)/1.d4 
         cncor1(4)=cncor1(4)/1.d3 
         cncor1(5)=cncor1(5)/1.d3 
      END IF 
  
! Definition of max interpolation error allowed for pcwlgi    
      DO 5 i=1,5 
      rmserx(i)=0.d0 
    5 END DO 
  
      DO 6 k=1,niers 
      DO 6 i=1,5 
      rmserx(i)=rmserx(i)+xiers(k,i)**2 
    6 CONTINUE 
  
      DO 7 i=1,5 
      rmserx(i)=1.D-11*SQRT(rmserx(i)/niers) 
    7 END DO 
  
      iicier=36 
  
      RETURN 
  
! Error termination   
   20 CONTINUE 
      WRITE(*,200) file(1:lf),nr 
  200 FORMAT('ERROR in reading file "',A,'" at record ',I5) 
      STOP '**** ierini: abnormal end ****' 
  
      END SUBROUTINE ierini  
! Copyright (C) 1999 by Mario Carpino (carpino@brera.mi.astro.it)       
! Version: November 12, 1999    
! --------------------------------------------------------------------- 
! 
!  *****************************************************************    
!  *   *    
!  *     I E R S T S       *    
!  *   *    
!  *IERS time series for polar motion, TDT-UT1 *    
!  *         and nutation correction *    
!  *   *    
!  *****************************************************************    
! 
! INPUT:    TJME      -  Modified Julian Day (ET)   
! NDER      -  Required order for derivatives (0/1/2)         
! 
! OUTPUT:   EOP(5)    -  Earth orientation parameters at TJME 
!        EOP(1) = X pole (arcsec)         
!        EOP(2) = Y pole (arcsec)         
!        EOP(3) = TDT-UT1 (s)   
!        EOP(4) = Dpsi (arcsec) 
!        EOP(5) = Depsilon (arcsec)       
! EOPD(5)   -  dEOP/dt (/d)     
! EOPDD(5)  -  d2EOP/dt2 (/d**2)
! 
      SUBROUTINE iersts(tjme,eop,eopd,eopdd,ntder) 
      USE chebi_pol
      IMPLICIT NONE 
  
      INTEGER ntder 
      DOUBLE PRECISION tjme 
      DOUBLE PRECISION eop(5),eopd(5),eopdd(5) 
  
      DOUBLE PRECISION eps 
      PARAMETER (eps=1.D-6) 
  
      INTEGER iret,i,iextr 
      DOUBLE PRECISION du,dud,dudd,dtbes 
  
! Work arrays for pcwlgi        
      DOUBLE PRECISION copler(0:lgintx,niplix,5),tlpler(2,niplix) 
      LOGICAL vdpler(niplix) 
  
      LOGICAL start 
      DOUBLE PRECISION bessep 
      EXTERNAL bessep 
  
      SAVE start,copler,tlpler,vdpler 
      DATA start/.true./ 
  
  
  
      IF(iicier.NE.36) THEN 
         CALL ierini 
         IF(iicier.NE.36) STOP ' **** iersts: internal error (01) ****' 
      END IF 
  
! Extrapolation       
      iextr=0 
      IF(extra) THEN 
         IF(tjme.LT.tint1) THEN 
            iextr=-1 
         ELSEIF(tjme.GT.tint2) THEN 
            iextr=1 
         END IF
      END IF 
  
      IF(iextr.EQ.0) THEN 
         CALL pcwlgi(xiers,tiers,niers,5,niersx,   &
     &      tjme,eop,eopd,eopdd,&
     &      nlpler,nvpler,nspler,nsmopl,ntder,start,rmserx,   &
     &      copler,tlpler,vdpler,lgintx,niplix,iret)
         IF(iret.NE.0) THEN 
            IF(extra) THEN 
               IF(ABS(tjme-tint1).LE.eps) THEN 
                  iextr=-1 
               ELSEIF(ABS(tjme-tint2).LE.eps) THEN 
                  iextr=1 
               END IF
               IF(iextr.EQ.0)STOP ' **** iersts: internal error (02) ****'      
            ELSE 
               WRITE(*,100) tjme 
               STOP '**** iersts: abnormal END ****' 
            END IF
         END IF
      END IF
  100 FORMAT(' iersts: time is out of bounds: TJME =',f12.5) 
  
      IF(iextr.NE.0) THEN 
         DO 2 i=1,5 
            eop(i)=0.d0 
            eopd(i)=0.d0 
            eopdd(i)=0.d0 
2        ENDDO
         eopd(3)=dutd 
         IF(iextr.EQ.-1) THEN 
            eop(3)=utd1+dutd*(tjme-tint1) 
         ELSEIF(iextr.EQ.1) THEN 
            eop(3)=utd2+dutd*(tjme-tint2) 
         ELSE 
            STOP '**** iersts: internal error (03) ****' 
         END IF
      END IF
  
! Subtracting applied smoothing 
      IF(iutsmo.EQ.0) THEN 
         CONTINUE 
      ELSEIF(iutsmo.EQ.1) THEN 
         CALL dut1r(tjme,du,dud,dudd,ntder) 
! (TDT-UT1) = (TDT-UT1R) - (UT1-UT1R)     
         eop(3)=eop(3)-du 
         IF(ntder.GE.1) eopd(3)=eopd(3)-dud 
         IF(ntder.GE.2) eopdd(3)=eopdd(3)-dudd 
      ELSEIF(iutsmo.EQ.2) THEN 
         CALL dut1s(tjme,du,dud,dudd,ntder) 
! (TDT-UT1) = (TDT-UT1S) - (UT1-UT1S)     
         eop(3)=eop(3)-du 
         IF(ntder.GE.1) eopd(3)=eopd(3)-dud 
         IF(ntder.GE.2) eopdd(3)=eopdd(3)-dudd 
      ELSE 
         STOP '**** iersts: internal error (04) ****' 
      END IF 
  
! Adding consistency correction 
      IF(cciera) THEN 
         dtbes=bessep(tjme)-cnep0 
         DO 1 i=1,5 
            eop(i)=eop(i)+cncor0(i)+cncor1(i)*dtbes 
            IF(ntder.GE.1) eopd(i)=eopd(i)+cncor1(i)/bessyr 
    1     CONTINUE 
      END IF 
  
      END SUBROUTINE iersts  
       
! Copyright (C) 1999 by Mario Carpino (carpino@brera.mi.astro.it)       
! Version: November 11, 1999    
! --------------------------------------------------------------------- 
! 
!  *****************************************************************    
!  *   *    
!  *     N U T A R G       *    
!  *   *    
!  *   Fundamental arguments of the IAU 1980 theory of nutation    *    
!  *   *    
!  *****************************************************************    
! 
! INPUT:    TJME      -  Modified Julian Day (TDT)  
! NTDER     -  Required order of time derivatives   
! 
! OUTPUT:   ARG(5)    -  Fundamental arguments (rad)
! ARGD(5)   -  dARG/dt (rad/d)  
! ARGDD(5)  -  d2ARG/dt2 (rad/d**2)       
! 
! The fundamental arguments are:
!   arg(1) = l     = mean anomaly of the Moon       
!   arg(2) = l'    = mean anomaly of the Sun        
!   arg(3) = F     = argument of latitude of the Moon         
!   arg(4) = D     = mean elongation of the Moon from the Sun 
!   arg(5) = Omega = mean longitude of the ascending node of the Moon   
! 
      SUBROUTINE nutarg(tjme,arg,argd,argdd,ntder)
      USE fund_const 
      IMPLICIT NONE 
  
      INTEGER ntder 
      DOUBLE PRECISION tjme,arg(5),argd(5),argdd(5) 
   
      INTEGER i 
      DOUBLE PRECISION s2r,d2c,d2c2,dl0,dl1,dl2,dl3,dp0,dp1,dp2,dp3 
      DOUBLE PRECISION df0,df1,df2,df3,dd0,dd1,dd2,dd3,dn0,dn1,dn2,dn3 
      DOUBLE PRECISION t1 
      LOGICAL first 
  
      SAVE first,d2c,d2c2,dl0,dl1,dl2,dl3,dp0,dp1,dp2,dp3 
      SAVE df0,df1,df2,df3,dd0,dd1,dd2,dd3,dn0,dn1,dn2,dn3 
  
      DOUBLE PRECISION ssum,dsum,ddsum,c0,c1,c2,c3,t 
      DATA first/.true./ 
  
! Computation of a 3rd degree polynomial and derivatives      
      ssum  (c0,c1,c2,c3,t) = ((c3*t+c2)*t+c1)*t+c0 
      dsum (c0,c1,c2,c3,t) = (3.d0*c3*t+2.d0*c2)*t+c1 
      ddsum(c0,c1,c2,c3,t) = 6.d0*c3*t+2.d0*c2 
  
      IF(ntder.LT.0 .OR. ntder.GT.2) STOP '**** nutarg: ntder = ? ****' 
  
      IF(first) THEN 
         first=.false. 
         s2r = pig/(180*3600) 
         d2c = 1.d0/36525.d0 
         d2c2= d2c**2 
! Coefficients of polynomial representation of fundamental arguments    
! from The Astronomical Almanac, 1984, page S26     
! 
! l = mean anomaly of the Moon  
         dl0 =     485866.733d0 * s2r 
         dl1 = 1717915922.633d0 * s2r 
         dl2 =         31.310d0 * s2r 
         dl3 =0.064d0 * s2r 
! l' = mean anomaly of the Sun  
         dp0 =    1287099.804d0 * s2r 
         dp1 =  129596581.224d0 * s2r 
         dp2 =        - 0.577d0 * s2r 
         dp3 =        - 0.012d0 * s2r 
! F = argument of latitude of the moon    
         df0 =     335778.877d0 * s2r 
         df1 = 1739527263.137d0 * s2r 
         df2 =       - 13.257d0 * s2r 
         df3 =0.011d0 * s2r 
! D = mean elongation of the Moon from the Sun      
         dd0 =    1072261.307d0 * s2r 
         dd1 = 1602961601.328d0 * s2r 
         dd2 =        - 6.891d0 * s2r 
         dd3 =0.019d0 * s2r 
! Omega = mean longitude of the ascending node of the Moon    
         dn0 =     450160.280d0 * s2r 
         dn1 =  - 6962890.539d0 * s2r 
         dn2 =7.455d0 * s2r 
         dn3 =0.008d0 * s2r 
      END IF 
  
! Time in Julian centuries since J2000.0  
      t1=(tjme-51544.5d0)/36525.d0 
  
      arg(1) = ssum(dl0,dl1,dl2,dl3,t1) 
      arg(2) = ssum(dp0,dp1,dp2,dp3,t1) 
      arg(3) = ssum(df0,df1,df2,df3,t1) 
      arg(4) = ssum(dd0,dd1,dd2,dd3,t1) 
      arg(5) = ssum(dn0,dn1,dn2,dn3,t1) 
  
      DO 1 i=1,5 
         arg(i)=MOD(arg(i),dpig) 
    1 END DO 
  
! First derivatives of the fundamental arguments    
      IF(ntder.GE.1) THEN 
         argd(1) = dsum(dl0,dl1,dl2,dl3,t1) * d2c 
         argd(2) = dsum(dp0,dp1,dp2,dp3,t1) * d2c 
         argd(3) = dsum(df0,df1,df2,df3,t1) * d2c 
         argd(4) = dsum(dd0,dd1,dd2,dd3,t1) * d2c 
         argd(5) = dsum(dn0,dn1,dn2,dn3,t1) * d2c 
      END IF 
  
! Second derivatives of the fundamental arguments   
      IF(ntder.GE.2) THEN 
         argdd(1) = ddsum(dl0,dl1,dl2,dl3,t1) * d2c2 
         argdd(2) = ddsum(dp0,dp1,dp2,dp3,t1) * d2c2 
         argdd(3) = ddsum(df0,df1,df2,df3,t1) * d2c2 
         argdd(4) = ddsum(dd0,dd1,dd2,dd3,t1) * d2c2 
         argdd(5) = ddsum(dn0,dn1,dn2,dn3,t1) * d2c2 
      END IF 
  
      END SUBROUTINE nutarg  
! Copyright (C) 1999 by Mario Carpino (carpino@brera.mi.astro.it)       
! Version: November 10, 1999    
! --------------------------------------------------------------------- 
! 
!  *****************************************************************    
!  *   *    
!  *      N U T N D        *    
!  *   *    
!  *  Nutation angles with time derivatives    *    
!  *   *    
!  *****************************************************************    
! 
! INPUT:    TJME      -  Modified Julian Day (ET)   
! NDER      -  Required order for derivatives (0/1/2)         
! 
! OUTPUT:   DPSI      -  Nutation in longitude Dpsi (rad)     
! DEPS      -  Nutation in obliquity Deps (rad)     
! DPSI1     -  dDpsi/dt (rad/s) (only if NDER >= 1) 
! DEPS1     -  dDeps/dt (rad/s) (only if NDER >= 1) 
! DPSI2     -  d^2 Dpsi/dt^2 (rad/s^2) (only if NDER = 2)     
! DEPS2     -  d^2 Deps/dt^2 (rad/s^2) (only if NDER = 2)     
! 
      SUBROUTINE nutnd(tjme,dpsi,deps,dpsi1,deps1,dpsi2,deps2,nder) 
      IMPLICIT NONE 
  
      INTEGER nder 
      DOUBLE PRECISION tjme,dpsi,deps,dpsi1,deps1,dpsi2,deps2 
  
! OPTIONS:  
! NUTCOR: correction to nutation theory   
!         0 - No correction     
!         1 - IERS correction   
      INTEGER nutcor 
      PARAMETER (nutcor=1) 
! Time tollerance (years)       
      DOUBLE PRECISION ttol 
      PARAMETER (ttol=1.D-8) 
  
      DOUBLE PRECISION fc0,fc1,fc2,eop(5),eopd(5),eopdd(5) 
      LOGICAL first 
  
  
      SAVE first,fc0,fc1,fc2 
      DATA first/.true./ 
  
      IF(first) THEN 
         first=.false. 
         fc0=ATAN(1.d0)/(45.d0*3600.d0) 
         fc1=fc0/86400.d0 
         fc2=fc1/86400.d0 
      END IF 
  
! Nominal model for nutation    
      CALL nutwhr(tjme,dpsi,deps,dpsi1,deps1,dpsi2,deps2,nder) 
  
! Corrections to the nominal model        
      IF(nutcor.LT.0.OR.nutcor.GT.1) STOP '**** nutnd: nutcor = ? ****' 
      IF(nutcor.EQ.1) THEN 
         CALL iersts(tjme,eop,eopd,eopdd,nder) 
         dpsi  = dpsi  + fc0 * eop(4) 
         deps  = deps  + fc0 * eop(5) 
         dpsi1 = dpsi1 + fc1 * eopd(4) 
         deps1 = deps1 + fc1 * eopd(5) 
         dpsi2 = dpsi2 + fc2 * eopdd(4) 
         deps2 = deps2 + fc2 * eopdd(5) 
      END IF 
  
      END SUBROUTINE nutnd  
! Copyright (C) 1999 by Mario Carpino (carpino@brera.mi.astro.it)       
! Version: November 8, 1999     
! --------------------------------------------------------------------- 
! 
!  *****************************************************************    
!  *   *    
!  *     N U T W H R       *    
!  *   *    
!  *    Nutation angles with derivatives       *    
!  *   (Wahr's nutation series for axis b for Earth model 1066a)   *    
!  *   *    
!  *   *    
!  *****************************************************************    
! 
! INPUT:    TJME      -  Modified Julian Time (ET)  
! NDER      -  Required order for derivatives (0/1/2)         
! 
! OUTPUT:   DPSI      -  Nutation in longitude (rad)
! DEPS      -  Nutation in obliquity (rad)
! DPSI1     -  dPSI/dt (rad/s) (only if NDER >= 1)  
! DEPS1     -  dEPS/dt (rad/s) (only if NDER >= 1)  
! DPSI2     -  d^2PSI/dt^2 (rad/s^2) (only if NDER = 2)       
! DEPS2     -  d^2EPS/dt^2 (rad/s^2) (only if NDER = 2)       
! 
      SUBROUTINE nutwhr(tjme,dpsi,deps,dpsi1,deps1,dpsi2,deps2,nder)
      USE fund_const 
      IMPLICIT NONE 
  
      INTEGER nder 
      DOUBLE PRECISION tjme,dpsi,deps,dpsi1,deps1,dpsi2,deps2 
  
      INTEGER i 
      DOUBLE PRECISION x(9,107) 
      DOUBLE PRECISION x1(180),x2(180),x3(180),x4(180),x5(180),x6(63) 
      DOUBLE PRECISION s2r,s2c,s2c2,dl0,dl1,dl2,dl3,dp0,dp1,dp2,dp3 
      DOUBLE PRECISION df0,df1,df2,df3,dd0,dd1,dd2,dd3,dn0,dn1,dn2,dn3 
      DOUBLE PRECISION t1,dl,dp,df,dd,dn,el,elp,f,d,omega 
      DOUBLE PRECISION el1,elp1,f1,d1,omega1,el2,elp2,f2,d2,omega2 
      DOUBLE PRECISION arg,coefp,coefe,sina,cosa,darg,dcoefp,dcoefe 
      DOUBLE PRECISION dsina,dcosa,ddarg,ddsina,ddcosa 
  
      LOGICAL first 
  
      EQUIVALENCE(x(1,  1),x1(1)) 
      EQUIVALENCE(x(1, 21),x2(1)) 
      EQUIVALENCE(x(1, 41),x3(1)) 
      EQUIVALENCE(x(1, 61),x4(1)) 
      EQUIVALENCE(x(1, 81),x5(1)) 
      EQUIVALENCE(x(1,101),x6(1)) 
  
      SAVE first,x,x1,x2,x3,x4,x5,x6,s2r,s2c,s2c2,dl0,dl1,dl2,dl3 
      SAVE dp0,dp1,dp2,dp3,df0,df1,df2,df3,dd0,dd1,dd2,dd3 
      SAVE dn0,dn1,dn2,dn3 
  
      DOUBLE PRECISION ssum,dsum,ddsum,c0,c1,c2,c3,t 
  
! Table of multiples of arguments (l,l',F,D,Omega) and coefficients     
! for the nutation in longitude and obliquity (unit =0.0001 arcsec);    
! see: The Astronomical Almanac, 1984, page S23     
      DATA x1/ 3.,  0.,  0.,  0.,  0.,       2.,    0.0,      0.,  0.0, &
     &         2.,  1.,  0., -2.,  0.,       1.,    0.0,      0.,  0.0, &
     &         2.,  0., -2.,  0.,  0.,      11.,    0.0,      0.,  0.0, &
     &         2.,  0.,  0., -2.,  0.,      48.,    0.0,      1.,  0.0, &
     &         2.,  0.,  0., -4.,  0.,      -1.,    0.0,      0.,  0.0, &
     &         2.,  0.,  0.,  2.,  0.,       1.,    0.0,      0.,  0.0, &
     &         2.,  0.,  0.,  0.,  0.,      29.,    0.0,     -1.,  0.0, &
     &         1., -1.,  0., -1.,  0.,      -3.,    0.0,      0.,  0.0, &
     &         1., -1.,  0., -2.,  0.,       1.,    0.0,      0.,  0.0, &
     &         1., -1.,  0.,  0.,  0.,       5.,    0.0,      0.,  0.0, &
     &         1.,  1.,  0., -2.,  0.,      -7.,    0.0,      0.,  0.0, &
     &         1.,  1.,  0.,  0.,  0.,      -3.,    0.0,      0.,  0.0, &
     &         1.,  0., -2., -2.,  0.,      -1.,    0.0,      0.,  0.0, &
     &         1.,  0., -2.,  2.,  0.,      -1.,    0.0,      0.,  0.0, &
     &         1.,  0., -2.,  0.,  0.,       4.,    0.0,      0.,  0.0, &
     &         1.,  0.,  2., -2.,  0.,      -1.,    0.0,      0.,  0.0, &
     &         1.,  0.,  2.,  0.,  0.,       3.,    0.0,      0.,  0.0, &
     &         1.,  0.,  0., -1.,  0.,      -4.,    0.0,      0.,  0.0, &
     &         1.,  0.,  0., -2.,  0.,    -158.,    0.0,     -1.,  0.0, &
     &         1.,  0.,  0., -4.,  0.,      -1.,    0.0,      0.,  0.0/ 
      DATA x2/ 1.,  0.,  0.,  2.,  0.,       6.,    0.0,      0.,  0.0, &
     &         1.,  0.,  0.,  0.,  0.,     712.,    0.1,     -7.,  0.0, &
     &         0.,  1., -2.,  2.,  0.,      -1.,    0.0,      0.,  0.0, &
     &         0.,  1.,  2., -2.,  0.,      -1.,    0.0,      0.,  0.0, &
     &         0.,  1.,  0., -2.,  0.,      -4.,    0.0,      0.,  0.0, &
     &         0.,  1.,  0.,  2.,  0.,      -1.,    0.0,      0.,  0.0, &
     &         0.,  1.,  0.,  1.,  0.,       1.,    0.0,      0.,  0.0, &
     &         0.,  0.,  2., -2.,  0.,     -22.,    0.0,      0.,  0.0, &
     &         0.,  0.,  2.,  0.,  0.,      26.,    0.0,     -1.,  0.0, &
     &         0.,  0.,  0.,  2.,  0.,      63.,    0.0,     -2.,  0.0, &
     &         0.,  0.,  0.,  1.,  0.,      -4.,    0.0,      0.,  0.0, &
     &        -1., -1.,  0.,  2.,  1.,       1.,    0.0,      0.,  0.0, &
     &        -1.,  0.,  2., -2.,  1.,      -2.,    0.0,      1.,  0.0, &
     &        -1.,  0.,  2.,  2.,  1.,     -10.,    0.0,      5.,  0.0, &
     &        -1.,  0.,  2.,  0.,  1.,      21.,    0.0,    -10.,  0.0, &
     &         0., -2.,  2., -2.,  1.,      -2.,    0.0,      1.,  0.0, &
     &        -1.,  0.,  0.,  2.,  1.,      16.,    0.0,     -8.,  0.0, &
     &        -1.,  0.,  0.,  1.,  1.,       1.,    0.0,      0.,  0.0, &
     &        -1.,  0.,  0.,  0.,  1.,     -58.,   -0.1,     32.,  0.0, &
     &        -2.,  0.,  2.,  0.,  1.,      46.,    0.0,    -24.,  0.0/ 
      DATA x3/-2.,  0.,  0.,  2.,  1.,      -6.,    0.0,      3.,  0.0, &
     &        -2.,  0.,  0.,  0.,  1.,      -2.,    0.0,      1.,  0.0, &
     &         2.,  0., -2.,  0.,  1.,       1.,    0.0,      0.,  0.0, &
     &         2.,  0.,  2., -2.,  1.,       1.,    0.0,     -1.,  0.0, &
     &         2.,  0.,  2.,  0.,  1.,      -5.,    0.0,      3.,  0.0, &
     &         2.,  0.,  0., -2.,  1.,       4.,    0.0,     -2.,  0.0, &
     &         2.,  0.,  0.,  0.,  1.,       2.,    0.0,     -1.,  0.0, &
     &         1.,  1.,  0., -2.,  1.,      -1.,    0.0,      0.,  0.0, &
     &         1.,  0.,  2., -2.,  1.,       6.,    0.0,     -3.,  0.0, &
     &         1.,  0.,  2.,  2.,  1.,      -1.,    0.0,      1.,  0.0, &
     &         1.,  0.,  2.,  0.,  1.,     -51.,    0.0,     27.,  0.0, &
     &         1.,  0.,  0., -2.,  1.,     -13.,    0.0,      7.,  0.0, &
     &         1.,  0.,  0.,  2.,  1.,      -1.,    0.0,      0.,  0.0, &
     &         1.,  0.,  0.,  0.,  1.,      63.,    0.1,    -33.,  0.0, &
     &         0., -1.,  2., -2.,  1.,      -5.,    0.0,      3.,  0.0, &
     &         0., -1.,  2.,  0.,  1.,      -1.,    0.0,      0.,  0.0, &
     &         0., -1.,  0.,  0.,  1.,     -12.,    0.0,      6.,  0.0, &
     &         0.,  1.,  2., -2.,  1.,       4.,    0.0,     -2.,  0.0, &
     &         0.,  1.,  2.,  0.,  1.,       1.,    0.0,      0.,  0.0, &
     &         0.,  1.,  0.,  0.,  1.,     -15.,    0.0,      9.,  0.0/ 
      DATA x4/ 0.,  0., -2.,  2.,  1.,       1.,    0.0,      0.,  0.0, &
     &         0.,  0., -2.,  0.,  1.,      -1.,    0.0,      0.,  0.0, &
     &         0.,  0.,  2., -2.,  1.,     129.,    0.1,    -70.,  0.0, &
     &         0.,  0.,  2.,  2.,  1.,      -7.,    0.0,      3.,  0.0, &
     &         0.,  0.,  2.,  0.,  1.,    -386.,   -0.4,    200.,  0.0, &
     &         0.,  0.,  0., -2.,  1.,      -5.,    0.0,      3.,  0.0, &
     &         0.,  0.,  0.,  2.,  1.,      -6.,    0.0,      3.,  0.0, &
     &         0.,  0.,  0.,  0.,  1.,       0.,    0.0,      0.,  0.0, &
     &        -1., -1.,  2.,  2.,  2.,      -3.,    0.0,      1.,  0.0, &
     &        -1.,  0.,  4.,  0.,  2.,       1.,    0.0,      0.,  0.0, &
     &        -1.,  0.,  2.,  4.,  2.,      -2.,    0.0,      1.,  0.0, &
     &        -1.,  0.,  2.,  2.,  2.,     -59.,    0.0,     26.,  0.0, &
     &        -1.,  0.,  2.,  0.,  2.,     123.,    0.0,    -53.,  0.0, &
     &        -1.,  0.,  0.,  0.,  2.,       1.,    0.0,     -1.,  0.0, &
     &        -2.,  0.,  2.,  4.,  2.,      -1.,    0.0,      1.,  0.0, &
     &        -2.,  0.,  2.,  2.,  2.,       1.,    0.0,     -1.,  0.0, &
     &        -2.,  0.,  2.,  0.,  2.,      -3.,    0.0,      1.,  0.0, &
     &         3.,  0.,  2., -2.,  2.,       1.,    0.0,      0.,  0.0, &
     &         3.,  0.,  2.,  0.,  2.,      -3.,    0.0,      1.,  0.0, &
     &         2.,  0.,  2., -2.,  2.,       6.,    0.0,     -3.,  0.0/ 
      DATA x5/ 2.,  0.,  2.,  2.,  2.,      -1.,    0.0,      0.,  0.0, &
     &         2.,  0.,  2.,  0.,  2.,     -31.,    0.0,     13.,  0.0, &
     &         1., -1.,  2.,  0.,  2.,      -3.,    0.0,      1.,  0.0, &
     &         1.,  1.,  2., -2.,  2.,       1.,    0.0,     -1.,  0.0, &
     &         1.,  1.,  2.,  0.,  2.,       2.,    0.0,     -1.,  0.0, &
     &         1.,  0.,  2., -2.,  2.,      29.,    0.0,    -12.,  0.0, &
     &         1.,  0.,  2.,  2.,  2.,      -8.,    0.0,      3.,  0.0, &
     &         1.,  0.,  2.,  0.,  2.,    -301.,    0.0,    129., -0.1, &
     &         1.,  0.,  0.,  0.,  2.,      -2.,    0.0,      1.,  0.0, &
     &         0., -1.,  2.,  2.,  2.,      -3.,    0.0,      1.,  0.0, &
     &         0., -1.,  2.,  0.,  2.,      -7.,    0.0,      3.,  0.0, &
     &         0.,  1.,  2.,  0.,  2.,       7.,    0.0,     -3.,  0.0, &
     &         0.,  1.,  0.,  0.,  2.,       1.,    0.0,      0.,  0.0, &
     &         0.,  0.,  4., -2.,  2.,       1.,    0.0,      0.,  0.0, &
     &         0.,  0.,  2., -1.,  2.,      -1.,    0.0,      0.,  0.0, &
     &         0.,  0.,  2.,  4.,  2.,      -1.,    0.0,      0.,  0.0, &
     &         0.,  0.,  2.,  2.,  2.,     -38.,    0.0,     16.,  0.0, &
     &         0.,  0.,  2.,  1.,  2.,       2.,    0.0,     -1.,  0.0, &
     &         0.,  0.,  2.,  0.,  2.,   -2274.,   -0.2,    977., -0.5, &
     &         0.,  0.,  0.,  0.,  2.,    2062.,    0.2,   -895.,  0.5/ 
      DATA x6/ 0.,  2.,  0.,  0.,  0.,      17.,   -0.1,      0.,  0.0, &
     &         0.,  1.,  0.,  0.,  0.,    1426.,   -3.4,     54., -0.1, &
     &         0., -1.,  2., -2.,  2.,     217.,   -0.5,    -95.,  0.3, &
     &         0.,  2.,  2., -2.,  2.,     -16.,    0.1,      7.,  0.0, &
     &         0.,  1.,  2., -2.,  2.,    -517.,    1.2,    224., -0.6, &
     &         0.,  0.,  2., -2.,  2.,  -13187.,   -1.6,   5736., -3.1, &
     &         0.,  0.,  0.,  0.,  1., -171996., -174.2,  92025.,  8.9/ 
  
      DATA first/.true./ 
  
! Computation of a 3rd degree polynomial with derivatives     
      ssum  (c0,c1,c2,c3,t) = ((c3*t+c2)*t+c1)*t+c0 
      dsum (c0,c1,c2,c3,t) = (3.d0*c3*t+2.d0*c2)*t+c1 
      ddsum(c0,c1,c2,c3,t) = 6.d0*c3*t+2.d0*c2 
  
      IF(nder.LT.0 .OR. nder.GT.2) STOP '**** nutwhr: nder = ? ****' 
  
      IF(first) THEN 
         first=.false. 
         s2r = pig/(180*3600) 
         s2c = 1.d0/(36525.d0*86400.d0) 
         s2c2= s2c**2 
! Polynomial representation of the fundamental arguments      
! see: The Astronomical Almanac, 1984, page S26     
         dl0 =     485866.733d0 * s2r 
         dl1 = 1717915922.633d0 * s2r 
         dl2 =         31.310d0 * s2r 
         dl3 =0.064d0 * s2r 
         dp0 =    1287099.804d0 * s2r 
         dp1 =  129596581.224d0 * s2r 
         dp2 =        - 0.577d0 * s2r 
         dp3 =        - 0.012d0 * s2r 
         df0 =     335778.877d0 * s2r 
         df1 = 1739527263.137d0 * s2r 
         df2 =       - 13.257d0 * s2r 
         df3 =0.011d0 * s2r 
         dd0 =    1072261.307d0 * s2r 
         dd1 = 1602961601.328d0 * s2r 
         dd2 =        - 6.891d0 * s2r 
         dd3 =0.019d0 * s2r 
         dn0 =     450160.280d0 * s2r 
         dn1 =  - 6962890.539d0 * s2r 
         dn2 =7.455d0 * s2r 
         dn3 =0.008d0 * s2r 
      END IF 
  
! Computation of fundamental arguments    
      t1=(tjme-51544.5d0)/36525.d0 
      dl = ssum(dl0,dl1,dl2,dl3,t1) 
      dp = ssum(dp0,dp1,dp2,dp3,t1) 
      df = ssum(df0,df1,df2,df3,t1) 
      dd = ssum(dd0,dd1,dd2,dd3,t1) 
      dn = ssum(dn0,dn1,dn2,dn3,t1) 
      el    = MOD(dl,dpig) 
      elp   = MOD(dp,dpig) 
      f     = MOD(df,dpig) 
      d     = MOD(dd,dpig) 
      omega = MOD(dn,dpig) 
  
! First derivatives of the fundamental arguments (rad/d)      
      IF(nder.GE.1) THEN 
         el1    = dsum(dl0,dl1,dl2,dl3,t1) * s2c 
         elp1   = dsum(dp0,dp1,dp2,dp3,t1) * s2c 
         f1     = dsum(df0,df1,df2,df3,t1) * s2c 
         d1     = dsum(dd0,dd1,dd2,dd3,t1) * s2c 
         omega1 = dsum(dn0,dn1,dn2,dn3,t1) * s2c 
      END IF 
  
! Second derivatives of the fundamental arguments (rad/d^2)   
      IF(nder.GE.2) THEN 
         el2    = ddsum(dl0,dl1,dl2,dl3,t1) * s2c2 
         elp2   = ddsum(dp0,dp1,dp2,dp3,t1) * s2c2 
         f2     = ddsum(df0,df1,df2,df3,t1) * s2c2 
         d2     = ddsum(dd0,dd1,dd2,dd3,t1) * s2c2 
         omega2 = ddsum(dn0,dn1,dn2,dn3,t1) * s2c2 
      END IF 
  
! Sum of the terms of the nutation series 
      dpsi  = 0.d0 
      deps  = 0.d0 
      dpsi1 = 0.d0 
      deps1 = 0.d0 
      dpsi2 = 0.d0 
      deps2 = 0.d0 
  
      DO 10 i = 1,107 
  
! Formation of multiples of arguments     
      arg    = x(1,i)*el  + x(2,i)*elp  + x(3,i)*f  + x(4,i)*d&
     &+ x(5,i)*omega  
      arg    = dmod(arg,dpig) 
  
! Formation of coefficients (first degree polynomials)        
      coefp  = x(6,i) + x(7,i)*t1 
      coefe  = x(8,i) + x(9,i)*t1 
  
! Evaluate nutation term        
      sina   =  SIN(arg) 
      cosa   =  COS(arg) 
      dpsi   = coefp*sina + dpsi 
      deps   = coefe*cosa + deps 
  
      IF(nder.GE.1) THEN 
! Formation of first derivatives of multiples of arguments    
         darg   = x(1,i)*el1 + x(2,i)*elp1 + x(3,i)*f1 + x(4,i)*d1     &
     &    + x(5,i)*omega1       
  
! Formation of first derivatives of coefficients    
         dcoefp = x(7,i) * s2c 
         dcoefe = x(9,i) * s2c 
  
! Evaluate first derivative of nutation term        
         dsina  =  cosa*darg 
         dcosa  = -sina*darg 
         dpsi1  = dcoefp*sina + coefp*dsina + dpsi1 
         deps1  = dcoefe*cosa + coefe*dcosa + deps1 
      END IF 
  
      IF(nder.GE.2) THEN 
! Formation of second derivatives of multiples of arguments   
         ddarg  = x(1,i)*el2 + x(2,i)*elp2 + x(3,i)*f2 + x(4,i)*d2     &
     &    + x(5,i)*omega2       
  
! Evaluate second derivative of nutation term       
         ddsina =  dcosa*darg + cosa*ddarg 
         ddcosa = -dsina*darg - sina*ddarg 
         dpsi2  = 2.d0*dcoefp*dsina + coefp*ddsina + dpsi2 
         deps2  = 2.d0*dcoefe*dcosa + coefe*ddcosa + deps2 
      END IF 
  
   10 END DO 
  
      dpsi  = dpsi  * 1.0d-4 * s2r 
      deps  = deps  * 1.0d-4 * s2r 
      dpsi1 = dpsi1 * 1.0d-4 * s2r 
      deps1 = deps1 * 1.0d-4 * s2r 
      dpsi2 = dpsi2 * 1.0d-4 * s2r 
      deps2 = deps2 * 1.0d-4 * s2r 
  
      END SUBROUTINE nutwhr  
! Copyright (C) 1999 by Mario Carpino (carpino@brera.mi.astro.it)       
! Version: November 8, 1999     
! --------------------------------------------------------------------- 
! 
!  *****************************************************************    
!  *   *    
!  *     O B L E Q D       *    
!  *   *    
!  *      Mean obliquity of ecliptic (with time derivatives)       *    
!  *  (see Astronomical Almanac 1987, B18)     *    
!  *   *    
!  *   *    
!  *****************************************************************    
! 
! INPUT:    TJME      -  Modified Julian Time (ET)  
! NDER      -  Required order for derivatives (0/1/2)         
! 
! OUTPUT:   OBL       -  Mean obliquity of ecliptic (rad)     
! OBL1      -  dOBL/dt (rad/s) (only if NDER >= 1)  
! OBL2      -  d^2OBL/dt^2 (rad/s^2) (only if NDER = 2)       
! 
      SUBROUTINE obleqd(tjme,obl,obl1,obl2,nder) 
      USE fund_const
      IMPLICIT NONE 
  
      INTEGER nder 
      DOUBLE PRECISION tjme,obl,obl1,obl2 
  
      DOUBLE PRECISION s2r,s2c,s2c2,ob0,ob1,ob2,ob3 
      DOUBLE PRECISION ts 
      LOGICAL first 
  
  
      SAVE first,s2r,s2c,s2c2,ob0,ob1,ob2,ob3 
  
      DOUBLE PRECISION ssum,dsum,ddsum,c0,c1,c2,c3,t 
      DATA first/.true./ 
  
! Computation of a 3rd degree polynomial with derivatives     
      ssum  (c0,c1,c2,c3,t) = ((c3*t+c2)*t+c1)*t+c0 
      dsum (c0,c1,c2,c3,t) = (3.d0*c3*t+2.d0*c2)*t+c1 
      ddsum(c0,c1,c2,c3,t) = 6.d0*c3*t+2.d0*c2 
  
      IF(nder.LT.0 .OR. nder.GT.2) STOP '**** obleqd: nder = ? ****' 
  
      IF(first) THEN 
         first=.false. 
         s2r=pig/(180.d0*3600.d0) 
         s2c = 1.d0/(36525.d0*86400.d0) 
         s2c2= s2c**2 
  
! Polynomial coefficients of the representation of OBL (see   
! Astronomical Almanac 1984, S21)         
         ob0 =  (23*3600+26*60+21.448d0) * s2r 
         ob1 =  -46.8150d0  * s2r 
         ob2 =  -0.00059d0  * s2r 
         ob3 =   0.001813d0 * s2r 
      END IF 
  
      ts  = (tjme-51544.5d0)/36525.d0 
      obl = ssum(ob0,ob1,ob2,ob3,ts) 
      IF(nder.GE.1) obl1 =  dsum(ob0,ob1,ob2,ob3,ts) * s2c 
      IF(nder.GE.2) obl2 = ddsum(ob0,ob1,ob2,ob3,ts) * s2c2 
  
      END SUBROUTINE obleqd  
! Copyright (C) 1999 by Mario Carpino (carpino@brera.mi.astro.it)       
! Version: November 9, 1999     
! --------------------------------------------------------------------- 
! 
!  *****************************************************************    
!  *   *    
!  *      P R E C D        *    
!  *   *    
!  *   Precession matrix (with time derivatives)         *    
!  *   *    
!  *****************************************************************    
! 
! INPUT:    TJME      -  Modified Julian Time (ET)  
! NDER      -  Required order for derivatives (0/1/2)         
! 
! OUTPUT:   ROT       -  Precession matrix P(t)     
! ROT1      -  First time derivative dP/dt (s^(-1)) 
!    (only if NDER >= 1)        
! ROT2      -  Second time derivative d^2P/dt^2 (s^(-2))      
!    (only if NDER = 2)         
! 
! The precession matrix P(t) performs the transformation from the       
! cartesian coordinates Xj2000 referred to the mean equator and equinox 
! of J2000 to the coordinates Xt referred to the mean equator and       
! equinox of date (corresponding to time t)         
! 
!  Xt = P(t) Xj2000   
! 
      SUBROUTINE precd(tjme,rot,rot1,rot2,nder)
      USE fund_const 
      IMPLICIT NONE 
  
      INTEGER nder 
      DOUBLE PRECISION tjme,rot(3,3),rot1(3,3),rot2(3,3) 
   
      DOUBLE PRECISION s2r,s2c,s2c2,tt,zeta,z,theta,zeta1,z1,theta1 
      DOUBLE PRECISION zeta2,z2,theta2 
      DOUBLE PRECISION zed,zd,thd,zedd,zdd,thdd,zeddd,zddd,thddd 
      DOUBLE PRECISION ra(3,3),rb(3,3),rc(3,3) 
      DOUBLE PRECISION ra1(3,3),rb1(3,3),rc1(3,3) 
      DOUBLE PRECISION ra2(3,3),rb2(3,3),rc2(3,3) 
      LOGICAL first 
  
  
      SAVE first,s2r,s2c,s2c2,zed,zd,thd,zedd,zdd,thdd,zeddd,zddd,thddd 
  
      DOUBLE PRECISION ssum,dsum,ddsum,c1,c2,c3,t 
      DATA first/.true./ 
  
! Computation of a 3rd  degree polynomial with derivatives    
      ssum  (c1,c2,c3,t) = ((c3*t+c2)*t+c1)*t 
      dsum (c1,c2,c3,t) = (3.d0*c3*t+2.d0*c2)*t+c1 
      ddsum(c1,c2,c3,t) = 6.d0*c3*t+2.d0*c2 
  
      IF(nder.LT.0 .OR. nder.GT.2) STOP '**** precd: nder = ? ****' 
  
! Time series for precession angles (zeta_A, z_A and theta_A):
! IAU 1980 coefficients (see Supplement to the Astronomical Almanac     
! 1984, page S19)     
      IF(first) THEN 
         first=.false. 
         s2r = pig/(180.d0*3600.d0) 
         s2c = 1.d0/(36525.d0*86400.d0) 
         s2c2= s2c**2 
! Linear term         
         zed   =   2306.2181d0 * s2r 
         zd    =   2306.2181d0 * s2r 
         thd   =   2004.3109d0 * s2r 
! Quadratic term      
         zedd  =   0.30188d0 * s2r 
         zdd   =   1.09468d0 * s2r 
         thdd  = - 0.42665d0 * s2r 
! Cubic term
         zeddd =   0.017998d0 * s2r 
         zddd  =   0.018203d0 * s2r 
         thddd = - 0.041833d0 * s2r 
      END IF 
  
! Computation of fundamental angles zeta_A, z_A and theta_A   
      tt    = ( tjme - 51544.5d0 ) / 36525.d0 
      zeta  = ssum(zed,zedd,zeddd,tt) 
      z     = ssum( zd, zdd, zddd,tt) 
      theta = ssum(thd,thdd,thddd,tt) 
  
! P(t) = R3 (-z_A) R2(theta_A) R3(-zeta_A)
      CALL rotmt(    -z , ra , 3) 
      CALL rotmt( theta , rb , 2) 
      CALL rotmt( -zeta , rc , 3) 
      IF(nder.GE.1) THEN 
         zeta1  = dsum(zed,zedd,zeddd,tt) * s2c 
         z1     = dsum( zd, zdd, zddd,tt) * s2c 
         theta1 = dsum(thd,thdd,thddd,tt) * s2c 
         CALL rotmt1(    -z , ra1 , 3,    -z1) 
         CALL rotmt1( theta , rb1 , 2, theta1) 
         CALL rotmt1( -zeta , rc1 , 3, -zeta1) 
      END IF 
      IF(nder.GE.2) THEN 
         zeta2  = ddsum(zed,zedd,zeddd,tt) * s2c2 
         z2     = ddsum( zd, zdd, zddd,tt) * s2c2 
         theta2 = ddsum(thd,thdd,thddd,tt) * s2c2 
         CALL rotmt2(    -z , ra2 , 3,    -z1,    -z2) 
         CALL rotmt2( theta , rb2 , 2, theta1, theta2) 
         CALL rotmt2( -zeta , rc2 , 3, -zeta1, -zeta2) 
      END IF 
  
      IF(nder.GE.2) rot2=rc2 !CALL assmat(rot2,rc2) 
      IF(nder.GE.1) rot1=rc1 !CALL assmat(rot1,rc1) 
      rot=rc !CALL assmat(rot,rc) 
  
      IF(nder.GE.2) CALL pd2mat(rb,rb1,rb2,rot,rot1,rot2) 
      IF(nder.GE.1) CALL pd1mat(rb,rb1,rot,rot1) 
      rot=MATMUL(rb,rot) ! CALL pdmat(rb,rot) 
  
      IF(nder.GE.2) CALL pd2mat(ra,ra1,ra2,rot,rot1,rot2) 
      IF(nder.GE.1) CALL pd1mat(ra,ra1,rot,rot1) 
      rot=MATMUL(ra,rot) ! CALL pdmat(ra,rot) 
  
      END SUBROUTINE precd  
! Copyright (C) 1999 by Mario Carpino (carpino@brera.mi.astro.it)       
! Version: November 16, 1999    
! --------------------------------------------------------------------- 
! 
!  *****************************************************************    
!  *   *    
!  *     R D B U L A       *    
!  *   *    
!  *  Read predictions of EOP time series      *    
!  *         from IERS Bulletin A    *    
!  *   *    
!  *****************************************************************    
! 
! INPUT:    FILE      -  Input file name  
! MJD2      -  Time (MJD,UTC) of last point read    
! TIERS     -  Time of data points (MJD, TDT)       
! XIERS     -  Value of data points:      
!        xiers(1) = X pole (arcsec)       
!        xiers(2) = Y pole (arcsec)       
!        xiers(3) = TDT-UT1 (s) 
!        xiers(4) = Dpsi (arcsec)         
!        xiers(5) = Depsilon (arcsec)     
! NIERS     -  No. of data points already contained in arrays 
! NIERSX    -  Physical DIMENSION of arrays         
! ISAMP     -  Undersampling factor       
! NPT       -  Number of data points read so far    
! IUTSMO    -  Smoothing used for UT1 (xiers(i,3)) :
!        0 = none     
!        1 = UT1R     
!        2 = UT1S     
! 
! OUTPUT:   MJD2      -  Time (MJD,UTC) of last point read    
! TIERS     -  (updated by adding values from Bulletin A)     
! XIERS     -  (updated by adding values from Bulletin A)     
! NIERS     -  (updated by adding values from Bulletin A)     
! 
      SUBROUTINE rdbula(file,mjd2,tiers,xiers,niers,niersx,isamp,npt,   &
     &        iutsmo) 
      IMPLICIT NONE 
  
      INTEGER mjd2,niers,niersx,isamp,npt,iutsmo 
      DOUBLE PRECISION tiers(niersx),xiers(niersx,5) 
      CHARACTER*(*) file 
  
      INTEGER unit,lf,nr,n,year,month,day,mjd,mjdc,mjdp,mjd1,mjdf,skip,i 
      DOUBLE PRECISION x,y,dt1,dtsec,dt,tjme,du,dud,dudd,dtsf 
      CHARACTER*100 rec 
      LOGICAL fill,strrec 
  
      INTEGER lench,itaiut 
      DOUBLE PRECISION tjm1  
      EXTERNAL lench,itaiut,tjm1
  
      lf=lench(file) 
      CALL filopl(unit,file) 
      nr=0 
      n=0 
  
! First scan of Bulletin A file: understand at which record data start  
      skip=0 
      strrec=.false. 
    1 CONTINUE 
      READ(unit,100,END=20) rec 
  100 FORMAT(A) 
      nr=nr+1 
      IF(isbadr(rec)) THEN 
         IF(strrec) THEN 
            IF(nr-skip.GE.15) GOTO 2 
         ELSE 
            skip=nr-1 
            strrec=.true. 
         END IF
      ELSE 
         skip=0 
         strrec=.false. 
      END IF
      GOTO 1 
  
    2 CONTINUE 
      IF(.NOT.strrec) GOTO 20 
  
! Skip non-data records         
      REWIND(unit) 
      DO 3 i=1,skip 
      READ(unit,100,END=20) rec 
    3 END DO 
      nr=skip 
  
    4 CONTINUE 
      nr=nr+1 
      READ(unit,100,END=10,ERR=30) rec 
      IF(.NOT.isbadr(rec)) GOTO 10 
      READ(rec,*,ERR=30) year,month,day,mjd,x,y,dt1 
      fill=.false. 
  
! Check on MJD value  
      mjdc=NINT(tjm1(day,month,year,0.D0)) 
      IF(mjdc.NE.mjd) THEN 
         WRITE(*,203) mjdc,mjd 
         GOTO 30 
      END IF 
  203 FORMAT('ERROR: inconsistent MJD:',2I8) 
  
! Check continuity with EOPC04 data       
      IF(n.EQ.0) THEN 
         mjd1=mjd 
! Why IERS issues predictions starting 2 days later the last record     
! present in their EOPC04 file? If it's so, here is the fix   
         IF(mjd.EQ.mjd2+2) THEN 
            fill=.true. 
         ELSEIF(mjd.NE.mjd2+1) THEN 
            WRITE(*,204) file(1:lf),mjd,mjd2 
         END IF
      ELSE 
! Check time consistency        
         IF(mjd-mjdp.NE.1) THEN 
            WRITE(*,205) mjd,mjdp 
            GOTO 30 
         END IF
      END IF 
  204 FORMAT('ERROR: beginning date of file "',A,'" (MJD=',I5,&
     &       ') is not consistent'/       &
     &       '       with ending date of EOPC04 data (MJD=',I5,')'/     &
     &       '       Try to get updated versions of IERS files')        
  205 FORMAT('WRONG time sequence:',2I7) 
  
      mjdp=mjd 
! TDT-UT1 = (TDT-TAI) + (TAI-UTC) + (UTC-UT1)       
      dtsec=etmtai+itaiut(mjd) 
      dt=dtsec-dt1 
! Transformation TJM(UTC) -> TJM(TDT)     
! TDT ( = ET ) = (ET-TAI) + (TAI-UTC) + UTC         
      tjme=mjd+dtsec/86400.d0 
! Applying required smoothing to TDT-UT1  
      IF(iutsmo.eq.0) THEN 
         CONTINUE 
      ELSEIF(iutsmo.EQ.1) THEN 
         CALL dut1r(tjme,du,dud,dudd,0) 
         dt=dt+du 
      ELSEIF(iutsmo.EQ.2) THEN 
         CALL dut1s(tjme,du,dud,dudd,0) 
         dt=dt+du 
      ELSE 
         STOP '**** rdbula: internal error (01) ****' 
      END IF 
! Interpolation of missing record         
      IF(fill) THEN 
         npt=npt+1 
         IF(MOD(npt,isamp).EQ.0) THEN 
            niers=niers+1 
            IF(niers.GT.niersx)         &
     &  STOP '**** rdbula: niers > niersx ****'     
            mjdf=mjd2+1 
            dtsf=etmtai+itaiut(mjdf) 
            tiers(niers)=mjdf+dtsf/86400.d0 
            xiers(niers,1)=(x+xiers(niers-1,1))/2 
            xiers(niers,2)=(y+xiers(niers-1,2))/2 
            xiers(niers,3)=(dt+xiers(niers-1,3))/2 
            xiers(niers,4)=0.D0 
            xiers(niers,5)=0.D0 
         END IF
      END IF
      npt=npt+1 
      IF(MOD(npt,isamp).NE.0) GOTO 4 
! Storing values in arrays      
      niers=niers+1 
      n=n+1 
      IF(niers.GT.niersx) STOP '**** rdbula: niers > niersx ****' 
      mjd2=mjd 
      tiers(niers)=tjme 
      xiers(niers,1)=x 
      xiers(niers,2)=y 
      xiers(niers,3)=dt 
      xiers(niers,4)=0.D0 
      xiers(niers,5)=0.D0 
      GOTO 4 
  
   10 CONTINUE 
      WRITE(*,206) n,mjd1,mjd2 
  206 FORMAT('    Bulletin A data:',I6,' points (from MJD=',  &
     &       I5,' to MJD=',I5,')')        
      CALL filclo(unit,' ') 
  
      RETURN 
  
   20 CONTINUE 
      WRITE(*,200) file(1:lf) 
  200 FORMAT('ERROR: sorry, unexpected format of file "',A,'"'/         &
     &       '       try running the program with ',&
     &       '"IERS.eopc04.bulA.use = .false."')    
      STOP '**** rdbula: abnormal end ****' 
  
   30 CONTINUE 
      WRITE(*,201) file(1:lf),nr 
  201 FORMAT('ERROR in reading file "',A,'" at record ',I5) 
     
      CONTAINS
     ! Copyright (C) 1999 by Mario Carpino (carpino@brera.mi.astro.it)       
! Version: November 16, 1999    
! --------------------------------------------------------------------- 
! 
!  *****************************************************************    
!  *   *    
!  *     I S B A D R       *    
!  *   *    
!  *Check whether a record is a data record    *    
!  *         from IERS Bulletin A    *    
!  *   *    
!  *****************************************************************    
! 
! INPUT:    REC       -  Input record     
! 
      LOGICAL FUNCTION isbadr(rec) 
      IMPLICIT NONE 
  
      CHARACTER*(*) rec 
  
      LOGICAL isnum,isbnum 
      EXTERNAL isnum,isbnum 
  
      isbadr=.false. 
  
! Example of record:  
! 
!       1999 11 12  51494        .0277       .3755       .41726         
!12345678901234567890123456789012345678901234567890123456789012         
!         1         2         3         4         5         6 
  
      IF(rec(1:7).NE.'       ') RETURN 
      IF(.NOT.isnum(rec(8:11))) RETURN 
      IF(rec(12:12).NE.' ') RETURN 
      IF(.NOT.isbnum(rec(13:14))) RETURN 
      IF(rec(15:15).NE.' ') RETURN 
      IF(.NOT.isbnum(rec(16:17))) RETURN 
      IF(rec(18:19).NE.'  ') RETURN 
      IF(.NOT.isnum(rec(20:24))) RETURN 
      IF(rec(25:28).NE.'    ') RETURN 
      IF(rec(33:33).NE.'.') RETURN 
      IF(.NOT.isbnum(rec(34:37))) RETURN 
      IF(rec(38:41).NE.'    ') RETURN 
      IF(rec(45:45).NE.'.') RETURN 
      IF(.NOT.isbnum(rec(46:49))) RETURN 
      IF(rec(50:53).NE.'    ') RETURN 
      IF(rec(57:57).NE.'.') RETURN 
      IF(.NOT.isbnum(rec(58:62))) RETURN 
  
      isbadr=.true. 
  
      END FUNCTION isbadr     
  
      END SUBROUTINE rdbula  
! Copyright (C) 1999 by Mario Carpino (carpino@brera.mi.astro.it)       
! Version: November 9, 1999     
! --------------------------------------------------------------------- 
! 
!  *****************************************************************    
!  *   *    
!  *      R N U T D        *    
!  *   *    
!  *     Nutation matrix according to Wahr (IAU 1980) theory       *    
!  *         (with time derivatives) *    
!  *   *    
!  *****************************************************************    
! 
! INPUT:    TJME      -  Modified Julian Time (ET)  
! NDER      -  Required order for derivatives (0/1/2)         
! 
! OUTPUT:   ROT       -  Nutation matrix N(Delta_psi,Delta_eps)         
! ROT1      -  First time derivative dN/dt (s^(-1)) 
!    (only if NDER >= 1)        
! ROT2      -  Second time derivative d^2N/dt^2 (s^(-2))      
!    (only if NDER = 2)         
! 
! The nutation matrix N(t) performs the transformation from the         
! cartesian coordinates referred to the mean equator and equinox        
! of date (at time t) Xmean to the coordinates referred to the true     
! equator and equinox of date (at the same time t) Xtrue      
! 
!  Xtrue = N(t) Xmean 
! 
      SUBROUTINE rnutd(tjme,rot,rot1,rot2,nder) 
      IMPLICIT NONE 
  
      INTEGER nder 
      DOUBLE PRECISION tjme,rot(3,3),rot1(3,3),rot2(3,3) 
  
      DOUBLE PRECISION epsm,epsm1,epsm2,epst,epst1,epst2 
      DOUBLE PRECISION dpsi,deps,dpsi1,deps1,dpsi2,deps2 
      DOUBLE PRECISION ra(3,3),rb(3,3),rc(3,3) 
      DOUBLE PRECISION ra1(3,3),rb1(3,3),rc1(3,3) 
      DOUBLE PRECISION ra2(3,3),rb2(3,3),rc2(3,3) 
  
      IF(nder.LT.0 .OR. nder.GT.2) STOP '**** rnutd: nder = ? ****' 
  
! Obliquity of the ecliptic eps(t)        
      CALL obleqd(tjme,epsm,epsm1,epsm2,nder) 
  
! Nutation angles Delta_psi(t) and Delta_eps(t)     
      CALL nutnd(tjme,dpsi,deps,dpsi1,deps1,dpsi2,deps2,nder) 
  
!  N(Delta_psi,Delta_eps) = R1(-eps-Delta_eps) R3(-Delta_psi) R1(eps)   
      epst = epsm + deps 
      CALL rotmt( -epst , ra , 1) 
      CALL rotmt( -dpsi , rb , 3) 
      CALL rotmt(  epsm , rc , 1) 
      IF(nder.GE.1) THEN 
         epst1 = epsm1 + deps1 
         CALL rotmt1( -epst , ra1 , 1, -epst1) 
         CALL rotmt1( -dpsi , rb1 , 3, -dpsi1) 
         CALL rotmt1(  epsm , rc1 , 1,  epsm1) 
      END IF 
      IF(nder.GE.2) THEN 
         epst2 = epsm2 + deps2 
         CALL rotmt2( -epst , ra2 , 1, -epst1, -epst2) 
         CALL rotmt2( -dpsi , rb2 , 3, -dpsi1, -dpsi2) 
         CALL rotmt2(  epsm , rc2 , 1,  epsm1,  epsm2) 
      END IF 
  
      IF(nder.GE.2) rot2=rc2 !CALL assmat(rot2,rc2) 
      IF(nder.GE.1) rot1=rc1 !CALL assmat(rot1,rc1) 
      rot=rc !CALL assmat(rot,rc) 
  
      IF(nder.GE.2) CALL pd2mat(rb,rb1,rb2,rot,rot1,rot2) 
      IF(nder.GE.1) CALL pd1mat(rb,rb1,rot,rot1) 
      rot=MATMUL(rb,rot) ! CALL pdmat(rb,rot) 
  
      IF(nder.GE.2) CALL pd2mat(ra,ra1,ra2,rot,rot1,rot2) 
      IF(nder.GE.1) CALL pd1mat(ra,ra1,rot,rot1) 
      rot=MATMUL(ra,rot) ! CALL pdmat(ra,rot) 
  
      END SUBROUTINE rnutd  
! Copyright (C) 1999 by Mario Carpino (carpino@brera.mi.astro.it)       
! Version: November 9, 1999     
! --------------------------------------------------------------------- 
! 
!  *****************************************************************    
!  *   *    
!  *      R O T P V        *    
!  *   *    
!  *        Computation of position and velocity vectors *    
!  *     in different reference frames         *    
!  *   *    
!  *   *    
!  *****************************************************************    
! 
! INPUT:    RSYS1     -  Initial reference system   
! DRAG1     -  Dragging of velocities in the initial
!    reference system (see below)         
! MJD1      -  Modified Julian Day (ET) of the initial        
!    reference system (integer part)      
! SEC1      -  Modified Julian Day (ET) of the initial        
!    reference system (seconds within the day)      
! P1        -  Position vector in RSYS1   
! V1        -  Velocity vector in RSYS1   
! RSYS2     -  Final reference system     
! DRAG2     -  Dragging of velocities in the final  
!    reference system (see below)         
! MJD2      -  Modified Julian Day (ET) of the final
!    reference system (integer part)      
! SEC2      -  Modified Julian Day (ET) of the final
!    reference system (seconds within the day)      
! 
! OUTPUT:   P2        -  Position vector in RSYS2   
! V2        -  Velocity vector in RSYS2   
! 
! "Initial" and "final" reference systems are specified through         
! the string (CHARACTER*4) variables RSYS1 and RSYS2, with the
! following conventions:        
!    'ECLM'  =  mean ecliptic and equinox 
!    'MEAN'  =  mean equator and equinox  
!    'TRUE'  =  true equator and equinox  
!    'PBF '  =  pseudo body-fixed (without pole motion and with uniform 
!     rotation rate)  
!    'BF  '  =  body-fixed (with pole motion and UT1 from IERS file)    
! 
! The flags DRAG1 and DRAG2 indicate whether the velocities in the      
! corresponding reference system are computed with (DRAGn=.true.)       
! or without (DRAGn=.false.) dragging, namely whether they are the real 
! velocities with respect to the considered frame or barely the inertial
! velocities rotated in the new frame.    
! 
! WARNING: if the final reference system is BF or PBF and the 
! transformation is performed with the aim of obtaining the initial     
! conditions for integrating the equations of motions in those systems, 
! the time of the reference system must always coincide with the time   
! at which position and velocity vectors are given. 
! 
      SUBROUTINE rotpv(rsys1,drag1,mjd1,sec1,p1,v1, &
     &       rsys2,drag2,mjd2,sec2,p2,v2) 
      IMPLICIT NONE 
  
      INTEGER mjd1,mjd2 
      DOUBLE PRECISION sec1,sec2,p1(3),v1(3),p2(3),v2(3) 
      CHARACTER*4 rsys1,rsys2 
      LOGICAL drag1,drag2 
  
      DOUBLE PRECISION rot(3,3),rot1(3,3),rot2(3,3),piner(3),viner(3) 
      DOUBLE PRECISION vtr(3) 
  
! SUBROUTINE ROTSYS rotates the axes, not the vectors; it gives         
! the derivatives of the rotation matrix with respect to the time of    
! the final frame     
! 
! STEP 1: from RSYS1 to J2000; computation of position and velocity     
! vectors in an intermediate inertial frame (Xin, Vin)        
! If the velocity vector in the initial reference frame V1 includes     
! dragging, then the transformation from the coordinates expressed      
! in the J2000 reference frame (Xin, Vin) to (X1, V1) is given by       
! 
!  X1 = R1 Xin        
!  V1 = dR1/dt Xin + R1 Vin     
! 
! where R1 is the rotation matrix which performs the transformation     
! J2000 -> RSYS1 ;the inverse transformation is therefore given by      
! 
!  Xin = R1' X1       
!  Vin = R1' V1 - R1' dR1/dt Xin
! 
! (where M' indicates the transpose of marix M)     
      IF(drag1) THEN 
         CALL rotsys('MEAN',mj2000,s2000,rsys1,mjd1,sec1,rot,rot1,     &
     &   rot2,1)      
         CALL trsp3(rot) 
         CALL prodmv(piner,rot,p1) 
         CALL prodmv(viner,rot,v1) 
         rot1=MATMUL(rot,rot1) !CALL pdmat(rot,rot1) 
         CALL prodmv(vtr,rot1,piner) 
         vtr=-1.d0*vtr ! CALL prodvs(vtr,-1.d0,vtr) 
         viner=vtr+viner ! CALL sumv(viner,vtr,viner) 
      ELSE 
! If the velocity in the initial reference frame does not include       
! dragging, the transformation equations are simply 
! 
!  X1 = R1 Xin        
!  V1 = R1 Vin        
! 
! and       
! 
!  Xin = R1' X1       
!  Vin = R1' V1       
         CALL rotsys('MEAN',mj2000,s2000,rsys1,mjd1,sec1,rot,rot1,     &
     &   rot2,0)      
         CALL trsp3(rot) 
         CALL prodmv(piner,rot,p1) 
         CALL prodmv(viner,rot,v1) 
      END IF 
  
! STEP 2: from J2000 to RSYS2; computation of position and velocity     
! vectors in the final reference frame (X2, V2) from the intermediate   
! inertial vectors (Xin, Vin)   
! 
! If the final reference frame includes dragging for the velocity       
! vector, the transformation is given by  
! 
!  X2 = R2 Xin        
!  V2 = dR2/dt Xin + R2 Vin     
! 
! where R2 is the rotation matrix which performs the transformation     
! J2000 -> RSYS2      
  
      IF(drag2) THEN 
         CALL rotsys('MEAN',mj2000,s2000,rsys2,mjd2,sec2,rot,rot1,     &
     &   rot2,1)      
         CALL prodmv(p2,rot,piner) 
         CALL prodmv(v2,rot,viner) 
         CALL prodmv(vtr,rot1,piner) 
         v2=vtr+v2 ! CALL sumv(v2,vtr,v2) 
      ELSE 
! otherwise, it is simply       
! 
!  X2 = R2 Xin        
!  V2 = R2 Vin        
         CALL rotsys('MEAN',mj2000,s2000,rsys2,mjd2,sec2,rot,rot1,     &
     &   rot2,0)      
         CALL prodmv(p2,rot,piner) 
         CALL prodmv(v2,rot,viner) 
      END IF 
  
      END SUBROUTINE rotpv  
! Copyright (C) 1999 by Mario Carpino (carpino@brera.mi.astro.it)       
! Version: November 9, 1999     
! --------------------------------------------------------------------- 
! 
!  *****************************************************************    
!  *   *    
!  *     R O T S Y S       *    
!  *   *    
!  *Rotation matrix for the transformation     *    
!  *  between different reference systems      *    
!  *   *    
!  *   *    
!  *****************************************************************    
! 
! INPUT:    RSYS1     -  Initial reference system   
! MJDE1     -  Date of the initial reference system t1        
!    (MJD, ET, integer part)    
! SECE1     -  Date of the initial reference system t1        
!    (MJD, ET, seconds within the day)    
! RSYS2     -  Final reference system     
! MJDE2     -  Date of the final reference system t2
!    (MJD, ET, integer part)    
! SECE2     -  Date of the final reference system t2
!    (MJD, ET, seconds within the day)    
! NDER      -  Required order for derivatives (0/1/2)         
! 
! OUTPUT:   ROT       -  Rotation matrix R
! ROT1      -  First time derivative dR/dt2 (s^(-1))
!    with respect to t2 and to an inertial reference
!    frame (only if NDER >= 1)  
! ROT2      -  Second time derivative d^2R/dt2^2 (s^(-2))     
!    with respect to t2 and to an inertial reference
!    frame(only if NDER = 2)    
! 
! The rotation matrix R performs the transformation from the cartesian  
! coordinates X1 referred to the "initial" reference system to the      
! cartesian coordinates X2 referred to the "final" reference system     
! 
!  X2 = R X1
! 
! The matrix R is computed as the product of two rotation matrices      
! R = R2(t2) R1(t1), where R1(t1) performs the transformation from      
! the "initial" reference system to J2000, and R2(t2) from J2000 to     
! the "final" reference system. The time derivatives of R are defined as
! 
!  dR/dt2     = R1(t1) dR2/dt2  
!  d^2R/dt2^2 = R1(t1) d^2R2/dt2^2        
! 
! Due to this definition, the derivatives of the rotation matrix can    
! be used directly in dynamic equations.  
! 
! "Initial" and "final" reference systems are specified through         
! the string (CHARACTER*4) variables RSYS1 and RSYS2, with the
! following conventions:        
!    'ECLM'  =  mean ecliptic and equinox 
!    'MEAN'  =  mean equator and equinox  
!    'TRUE'  =  true equator and equinox  
!    'PBF '  =  pseudo body-fixed (without pole motion and with uniform 
!     rotation rate)  
!    'BF  '  =  body-fixed (with pole motion and UT1 from IERS file)    
! 
      SUBROUTINE rotsys(rsys1,mjde1,sece1,rsys2,mjde2,sece2,  &
     &        rot,rot1,rot2,nder)         
      IMPLICIT NONE 
  
      INTEGER mjde1,mjde2,nder 
      DOUBLE PRECISION sece1,sece2,rot(3,3),rot1(3,3),rot2(3,3) 
      CHARACTER*4 rsys1,rsys2 
  
! Time tollerance     
      DOUBLE PRECISION eps 
      PARAMETER (eps=1.D-6) 
  
      INTEGER mjd1,mjd2,mjd,i,j 
      DOUBLE PRECISION sec1,sec2,sec,t1,t2,t,dt,obl,obl1,obl2,gmst,gmst1 
      DOUBLE PRECISION eqeq,eqeq1,eqeq2,gast,gast1,gast2 
      DOUBLE PRECISION r(3,3),r1(3,3),r2(3,3) 
      CHARACTER*4 rsys 
      LOGICAL gomean 
  
! Check of input parameters     
      IF(nder.LT.0 .OR. nder.GT.2) STOP '**** rotsys: nder = ? ****' 
      IF(rsys1.NE.'ECLM'.AND.rsys1.ne.'MEAN'.and.rsys1.ne.'TRUE'        &
     &        .AND.rsys1.NE.'PBF '.and.rsys1.ne.'BF  ')       &
     &         STOP '**** rotsys: unknown RSYS1 ****'         
      IF(rsys2.NE.'ECLM'.AND.rsys2.ne.'MEAN'.and.rsys2.ne.'TRUE'        &
     &        .AND.rsys2.NE.'PBF '.and.rsys2.ne.'BF  ')       &
     &         STOP '**** rotsys: unknown RSYS2 ****'         
  
! Time normalization  
      mjd1=mjde1 
      sec1=sece1 
      CALL timnf(mjd1,sec1,'ET ') 
      mjd2=mjde2 
      sec2=sece2 
      CALL timnf(mjd2,sec2,'ET ') 
  
! GENERAL CONCEPT: the complete rotation matrix R is built up by chain  
! multiplication of single rotation matrices which perform the total    
! transformation step by step; the route from the initial to the final  
! reference frames is computed by following the graph:        
! 
!     PBF  PBF        
!      |    |         
!      BF  ---  TRUE TRUE  ---  BF        
!      |    |         
!     MEAN   -----   MEAN       
!      |    |         
!     ECLM ECLM       
! 
!    (t=t1)         (t=t2)      
! 
! where the transformation between reference system which are adjacent  
! in the graph is standard (and is given by separate routines). As an   
! example, the transformation between BF(t1) and PBF(t2) is   
! accomplished with the following intermediate passages:      
! 
! BF(t1) -> TRUE(t1) -> MEAN(t1) -> MEAN(t_2) -> TRUE(t_2) -> PBF(t_2)  
! 
! however, if t2=t1 and no derivatives are required, the following      
! route is sufficient:
! 
! BF(t1) -> TRUE(t1) -> PBF(t1) 
! 
! the computation of the derivatives of the rotation matrix always      
! requires the passage through the inertial (MEAN) reference frame,     
! in order to be able to compute dR2/dt2. 
  
      t1=mjd1+sec1/86400.d0 
      t2=mjd2+sec2/86400.d0 
  
! The total rotation matrix from RSYS1 to RSYS2 is built up by going    
! along the trasformation tree step by step and multiplying the rotation
! matrix previously obtained by the matrix corresponding to the new     
! step. The variables RSYS1, MJD, SEC and T keep memory of the point    
! reached so far. At the beginning, this is coincident with the initial 
! reference system    
  
      rsys=rsys1 
      mjd=mjd1 
      sec=sec1 
      t=t1 
  
! The initial value of the total rotation matrix R is the unit matrix;  
! all its derivatives are zero  

      rot1=0.d0
      rot2=0.d0 !CALL assmat(rot2,rot1) 
      rot=0.d0 ! CALL assmat(rot,rot1) 
      DO  i=1,3 
         rot(i,i)=1.d0 
      ENDDO 
      dt=(mjd-mjd2)*86400.D0+sec-sec2 
  
! The flag GOMEAN indicates whether the direction of the trasformation  
! is approaching J2000 mean reference system from RSYS1 (phase 1),      
! or going away from J2000 towards RSYS2 (phase 2). 
! Passage through J2000 is necessary when:
! a) t1 is different from t2 and the application of the       
!    precession matrix is required;       
! b) the derivatives of the rotation matrix are required (theese        
!    must be always computed from the mean, inertial reference system)  
  
      gomean=(ABS(dt).GT.eps .OR. nder.GE.1) 
    3 CONTINUE 
  
! PHASE 1: path from RSYS1 to J2000. In this phase the time derivatives 
! of the rotation matrix are NOT updated at each step, since the total  
! derivative is computed with respect to J2000      
  
      IF(gomean) THEN 
IF(rsys.EQ.'ECLM') THEN 
! Trasformation ECLM -> MEAN: rotation around x1 axis by the negative   
! value of the mean obliquity of the ecliptic       
    CALL obleqd(t,obl,obl1,obl2,0) 
    CALL rotmt(-obl,r,1) 
    rot=MATMUL(r,rot) ! CALL pdmat(r,rot) 
    rsys='MEAN' 
ELSEIF(rsys.EQ.'TRUE') THEN 
! Trasformation TRUE -> MEAN: transpose of the nutation matrix
    CALL rnutd(t,r,r1,r2,0) 
    CALL trsp3(r) 
    rot=MATMUL(r,rot) ! CALL pdmat(r,rot) 
    rsys='MEAN' 
ELSEIF(rsys.EQ.'PBF ') THEN 
! Trasformation PBF -> TRUE: rotation around x3 axis by the negative    
! value of the nominal Greenwich Apparent Sidereal Time       
    CALL gmsnom(mjd,sec,gmst,gmst1) 
    CALL equeqd(t,eqeq,eqeq1,eqeq2,0) 
    gast=gmst+eqeq 
    CALL rotmt(-gast,r,3) 
    rot=MATMUL(r,rot) ! CALL pdmat(r,rot) 
    rsys='TRUE' 
ELSEIF(rsys.EQ.'BF  ') THEN 
! Trasformation BF -> TRUE: transpose of the diurnal rotation matrix    
    CALL diurot(mjd,sec,r,r1,r2,0) 
    CALL trsp3(r) 
    rot=MATMUL(r,rot) ! CALL pdmat(r,rot) 
    rsys='TRUE' 
ELSEIF(rsys.EQ.'MEAN') THEN 
! Trasformation MEAN(t1) -> MEAN(t2): precession matrix       
    CALL precd(t,r,r1,r2,0) 
    CALL trsp3(r) 
    rot=MATMUL(r,rot) ! CALL pdmat(r,rot) 
! In this point J2000 reference frame has been reached; here starts the 
! computation of the derivatives of the rotation matrix       
    CALL precd(t2,r,r1,r2,nder) 
    IF(nder.GE.2) CALL pd2mat(r,r1,r2,rot,rot1,rot2) 
    IF(nder.GE.1) CALL pd1mat(r,r1,rot,rot1) 
    rot=MATMUL(r,rot) ! CALL pdmat(r,rot) 
    mjd=mjd2 
    sec=sec2 
    t=t2 
    gomean=.false. 
ELSE 
    STOP '**** rotsys: internal error (01) ****' 
END IF 
  
! PHASE 2: path from J2000 to RSYS2. In this phase the time derivatives 
! of the rotation matrix are updated at each step   
      ELSEIF(rsys.NE.rsys2) THEN 
IF(rsys.EQ.'TRUE') THEN 
    IF(rsys2.EQ.'MEAN'.OR.rsys2.eq.'ECLM') THEN 
! Trasformation TRUE -> MEAN: transpose of the nutation matrix
        CALL rnutd(t,r,r1,r2,nder) 
        CALL trsp3(r) 
        IF(nder.GE.1) CALL trsp3(r1) 
        IF(nder.GE.2) CALL trsp3(r2) 
        IF(nder.GE.2) CALL pd2mat(r,r1,r2,rot,rot1,rot2) 
        IF(nder.GE.1) CALL pd1mat(r,r1,rot,rot1) 
        rot=MATMUL(r,rot) ! CALL pdmat(r,rot) 
        rsys='MEAN' 
    ELSEIF(rsys2.EQ.'PBF ') THEN 
! Trasformation TRUE -> PBF: rotation around x3 axis by the nominal     
! Greenwich Apparent Sidereal Time        
        CALL gmsnom(mjd,sec,gmst,gmst1) 
        CALL equeqd(t,eqeq,eqeq1,eqeq2,nder) 
        gast=gmst+eqeq 
        CALL rotmt(gast,r,3) 
        IF(nder.GE.1) THEN 
  gast1=gmst1+eqeq1 
  CALL rotmt1(gast,r1,3,gast1) 
        END IF 
        IF(nder.GE.2) THEN 
  gast2=eqeq2 
  CALL rotmt2(gast,r2,3,gast1,gast2) 
        END IF 
        IF(nder.GE.2) CALL pd2mat(r,r1,r2,rot,rot1,rot2) 
        IF(nder.GE.1) CALL pd1mat(r,r1,rot,rot1) 
        rot=MATMUL(r,rot) ! CALL pdmat(r,rot) 
        rsys='PBF ' 
    ELSEIF(rsys2.EQ.'BF  ') THEN 
! Trasformation TRUE -> BF: diurnal rotation matrix 
        CALL diurot(mjd,sec,r,r1,r2,nder) 
        IF(nder.GE.2) CALL pd2mat(r,r1,r2,rot,rot1,rot2) 
        IF(nder.GE.1) CALL pd1mat(r,r1,rot,rot1) 
        rot=MATMUL(r,rot) ! CALL pdmat(r,rot) 
        rsys='BF  ' 
    ELSE 
        STOP ' **** rotsys: internal error (02) ****' 
    END IF 
ELSEIF(rsys.EQ.'MEAN') THEN 
    IF(rsys2.EQ.'TRUE'.OR.rsys2.eq.'PBF ' &
     &      .OR.rsys2.EQ.'BF  ') THEN     
! Trasformation MEAN -> TRUE: nutation matrix       
        CALL rnutd(t,r,r1,r2,nder) 
        IF(nder.GE.2) CALL pd2mat(r,r1,r2,rot,rot1,rot2) 
        IF(nder.GE.1) CALL pd1mat(r,r1,rot,rot1) 
        rot=MATMUL(r,rot) ! CALL pdmat(r,rot) 
        rsys='TRUE' 
    ELSEIF(rsys2.EQ.'ECLM') THEN 
! Trasformation MEAN -> ECLM: rotation around x1 axis by the mean       
! obliquity of the ecliptic     
        CALL obleqd(t,obl,obl1,obl2,nder) 
        CALL rotmt(obl,r,1) 
        IF(nder.GE.1) CALL rotmt1(obl,r1,1,obl1) 
        IF(nder.GE.2) CALL rotmt2(obl,r2,1,obl1,obl2) 
        IF(nder.GE.2) CALL pd2mat(r,r1,r2,rot,rot1,rot2) 
        IF(nder.GE.1) CALL pd1mat(r,r1,rot,rot1) 
        rot=MATMUL(r,rot) ! CALL pdmat(r,rot) 
        rsys='ECLM' 
    ELSE 
        STOP ' **** rotsys: internal error (03) ****' 
    END IF 
ELSEIF(rsys.EQ.'ECLM') THEN 
! Trasformation ECLM -> MEAN: rotation around x1 axis by the negative   
! value of the mean obliquity of the ecliptic       
    CALL obleqd(t,obl,obl1,obl2,nder) 
    CALL rotmt(-obl,r,1) 
    IF(nder.GE.1) CALL rotmt1(-obl,r1,1,-obl1) 
    IF(nder.GE.2) CALL rotmt2(-obl,r2,1,-obl1,-obl2) 
    IF(nder.GE.2) CALL pd2mat(r,r1,r2,rot,rot1,rot2) 
    IF(nder.GE.1) CALL pd1mat(r,r1,rot,rot1) 
    rot=MATMUL(r,rot) ! CALL pdmat(r,rot) 
    rsys='MEAN' 
ELSEIF(rsys.EQ.'PBF ') THEN 
! Trasformation PBF -> TRUE: rotation around x3 axis by the negative    
! value of the nominal Greenwich Apparent Sidereal Time       
    CALL gmsnom(mjd,sec,gmst,gmst1) 
    CALL equeqd(t,eqeq,eqeq1,eqeq2,nder) 
    gast=gmst+eqeq 
    CALL rotmt(-gast,r,3) 
    IF(nder.GE.1) THEN 
        gast1=gmst1+eqeq1 
        CALL rotmt1(-gast,r1,3,-gast1) 
    END IF 
    IF(nder.GE.2) THEN 
        gast2=eqeq2 
        CALL rotmt2(-gast,r2,3,-gast1,-gast2) 
    END IF 
    IF(nder.GE.2) CALL pd2mat(r,r1,r2,rot,rot1,rot2) 
    IF(nder.GE.1) CALL pd1mat(r,r1,rot,rot1) 
    rot=MATMUL(r,rot) ! CALL pdmat(r,rot) 
    rsys='TRUE' 
ELSEIF(rsys.EQ.'BF  ') THEN 
! Trasformation BF -> TRUE: transpose of the diurnal rotation matrix    
    CALL diurot(mjd,sec,r,r1,r2,nder) 
    CALL trsp3(r) 
    IF(nder.GE.1) CALL trsp3(r1) 
    IF(nder.GE.2) CALL trsp3(r2) 
    IF(nder.GE.2) CALL pd2mat(r,r1,r2,rot,rot1,rot2) 
    IF(nder.GE.1) CALL pd1mat(r,r1,rot,rot1) 
    rot=MATMUL(r,rot) ! CALL pdmat(r,rot) 
    rsys='TRUE' 
ELSE 
    STOP ' **** rotsys: internal error (04) ****' 
END IF 
      ELSE 
RETURN 
      END IF 
      GOTO 3 
  
      END SUBROUTINE rotsys  
! Copyright (C) 1999 by Mario Carpino (carpino@brera.mi.astro.it)       
! Version: November 11, 1999    
! --------------------------------------------------------------------- 
! 
!  *****************************************************************    
!  *   *    
!  *      X Y P O L        *    
!  *   *    
!  *     Computation of the coordinates of the terrestrial pole    *    
!  *       with their time derivatives         *    
!  *   *    
!  *****************************************************************    
! 
! INPUT:    TJME      -  Modified Julian Day (TDT)  
! NTDER     -  Required order of time derivatives (0/1/2)     
! 
! OUTPUT:   XPOL      -  (x,y) pole coordinates (rad)         
! X1POL     -  d(x,y)/dET (rad/s) (only if NDER >= 1)         
! X2POL     -  d^2 (x,y)/dET^2 (rad/s^2) (only if NDER = 2)   
! 
      SUBROUTINE xypol(tjme,xpol,x1pol,x2pol,nder) 
      USE fund_const
      IMPLICIT NONE 
  
      INTEGER nder 
      DOUBLE PRECISION tjme,xpol(2),x1pol(2),x2pol(2) 
   
      DOUBLE PRECISION c0,c1,c2,eop(5),eopd(5),eopdd(5) 
      INTEGER i 
      LOGICAL first 
  
      SAVE first,c0,c1,c2 
      DATA first/.true./ 
  
      IF(first) THEN 
         first=.false. 
         c0=pig/(180.d0*3600.d0) 
         c1=c0/86400.d0 
         c2=c1/86400.d0 
      END IF 
  
! Interpolation of IERS data    
      CALL iersts(tjme,eop,eopd,eopdd,nder) 
  
! Transformation of units from  arcsec, arcsec/d, arcsec/d^2  
! to rad, rad/s, rad/s^2        
      DO 1 i=1,2 
      xpol(i)=eop(i)*c0 
      x1pol(i)=eopd(i)*c1 
      x2pol(i)=eopdd(i)*c2 
    1 END DO 
  
      END SUBROUTINE xypol  
 
      
! Copyright (C) 1999 by Mario Carpino (carpino@brera.mi.astro.it)       
! Version: November 10, 1999    
! --------------------------------------------------------------------- 
! 
!  *****************************************************************    
!  *   *    
!  *     P D 1 M A T       *    
!  *   *    
!  *      First time derivative of the product of two matrices     *    
!  *   *    
!  *****************************************************************    
! 
! INPUT:    R1        -  R1 matrix        
! R1D       -  dR1/dt 
! R2        -  R2 matrix        
! R2D       -  dR2/dt 
! 
! OUTPUT:   R2D       -  d(R1 R2)/dt      
! 
! The time derivative of the product R1 R2 is computed according        
! to the formula      
! 
! d(R1 R2)/dt = dR1/dt R2 + R1 dR2/dt     
! 
      SUBROUTINE pd1mat(r1,r1d,r2,r2d) 
      IMPLICIT NONE 
  
      DOUBLE PRECISION r1(3,3),r1d(3,3),r2(3,3),r2d(3,3) 
  
      DOUBLE PRECISION p1(3,3),p2(3,3) 
  
      CALL prodmm(p1,r1d,r2) 
      CALL prodmm(p2,r1,r2d) 
      r2d=p1+p2 ! CALL summat(r2d,p1,p2) 
  
      END SUBROUTINE pd1mat  
! Copyright (C) 1999 by Mario Carpino (carpino@brera.mi.astro.it)       
! Version: November 10, 1999    
! --------------------------------------------------------------------- 
! 
!  *****************************************************************    
!  *   *    
!  *     P D 2 M A T       *    
!  *   *    
!  *     Second time derivative of the product of two matrices     *    
!  *   *    
!  *****************************************************************    
! 
! INPUT:    R1        -  R1 matrix        
! R1D       -  dR1/dt 
! R1DD      -  d^2 R1/dt^2      
! R2        -  R2 matrix        
! R2D       -  dR2/dt 
! R2DD      -  d^2 R2/dt^2      
! 
! OUTPUT:   R2DD      -  d^2 (R1 R2)/dt^2 
! 
! The time derivative of the product R1 R2 is computed according        
! to the formula      
! 
! d^2 (R1 R2)/dt^2 =  
!       = d^2 R1/dt^2 R2 + 2 (dR1/dt)*(dR2/dt) + R1 d^2 R2/dt^2         
! 
      SUBROUTINE pd2mat(r1,r1d,r1dd,r2,r2d,r2dd) 
      IMPLICIT NONE 
  
      DOUBLE PRECISION r1(3,3),r1d(3,3),r1dd(3,3) 
      DOUBLE PRECISION r2(3,3),r2d(3,3),r2dd(3,3) 
  
      INTEGER i,j 
      DOUBLE PRECISION p1(3,3),p2(3,3),p3(3,3) 
  
      CALL prodmm(p1,r1dd,r2) 
      CALL prodmm(p2,r1d,r2d) 
      CALL prodmm(p3,r1,r2dd) 
  
      DO 1 i=1,3 
      DO 1 j=1,3 
      r2dd(i,j)=p1(i,j)+2.d0*p2(i,j)+p3(i,j) 
    1 CONTINUE 
  
      END SUBROUTINE pd2mat  

! Copyright (C) 1999 by Mario Carpino (carpino@brera.mi.astro.it)       
! Version: November 10, 1999    
! --------------------------------------------------------------------- 
! 
!  *****************************************************************    
!  *   *    
!  *     R O T M T 1       *    
!  *   *    
!  * First time derivative of a rotation matrix*    
!  *   *    
!  *****************************************************************    
! 
! INPUT:    A         -  Rotation angle alpha (rad) 
! K         -  (=1,2,3) Rotation axis (1=x, 2=y, 3=z)         
! ADOT      -  First time derivative of the rotation angle    
!    dalpha/dt (rad/arbitrary time unit)  
! 
! OUTPUT:   R         -  First time derivative of the rotation matrix   
!    dR_k/dt (rad/arbitrary time unit)    
! 
! The time derivative of the rotation matrix R_k(alpha) (for its        
! definition, see SUBROUTINE ROTMT) is computed according to the formula
! 
! dR_k/dt = dR_k/dalpha * dalpha/dt       
! 
      SUBROUTINE rotmt1(a,r,k,adot) 
      IMPLICIT NONE 
  
      INTEGER k 
      DOUBLE PRECISION a,r(3,3),adot 
  
      INTEGER i1,i2,i3 
      DOUBLE PRECISION cosa,sina 
  
      IF(k.LT.1 .OR. k.GT.3) STOP '**** rotmt1: k = ??? ****' 
  
      cosa=COS(a) 
      sina=SIN(a) 
      i1=k 
      IF(i1.GT.3) i1=i1-3 
      i2=i1+1 
      IF(i2.GT.3) i2=i2-3 
      i3=i2+1 
      IF(i3.GT.3) i3=i3-3 
  
      r(i1,i1) =  0.d0 
      r(i1,i2) =  0.d0 
      r(i1,i3) =  0.d0 
      r(i2,i1) =  0.d0 
      r(i2,i2) = -sina*adot 
      r(i2,i3) =  cosa*adot 
      r(i3,i1) =  0.d0 
      r(i3,i2) = -cosa*adot 
      r(i3,i3) = -sina*adot 
  
      END SUBROUTINE rotmt1  
! Copyright (C) 1999 by Mario Carpino (carpino@brera.mi.astro.it)       
! Version: November 10, 1999    
! --------------------------------------------------------------------- 
! 
!  *****************************************************************    
!  *   *    
!  *     R O T M T 2       *    
!  *   *    
!  * Second time derivative of a rotation matrix         *    
!  *   *    
!  *****************************************************************    
! 
! INPUT:    A         -  Rotation angle alpha (rad) 
! K         -  (=1,2,3) Rotation axis (1=x, 2=y, 3=z)         
! ADOT      -  First time derivative of the rotation angle    
!    dalpha/dt (rad/arbitrary time unit)  
! ADDOT     -  Second time derivative of the rotation angle   
!    d^2 alpha/dt^2 (rad/(arbitrary time unit)^2)   
! 
! OUTPUT:   R         -  Second time derivative of the rotation matrix  
!    d^2 R_k/dt^2 (rad/(arbitrary time unit)^2)     
! 
! The time derivative of the rotation matrix R_k(alpha) (for its        
! definition, see SUBROUTINE ROTMT) is computed according to the formula
! 
! dR^2_k/dt^2 =       
!  = d^2 R_k/dalpha^2 * (dalpha/dt)^2 + (dR_k/dalpha) * (d^2 alpha/dt^2)
! 
      SUBROUTINE rotmt2(a,r,k,adot,addot) 
      IMPLICIT NONE 
  
      INTEGER k 
      DOUBLE PRECISION a,r(3,3),adot,addot 
  
      INTEGER i1,i2,i3 
      DOUBLE PRECISION cosa,sina,adot2 
  
      IF(k.LT.1 .OR. k.GT.3) STOP '**** rotmt2: k = ??? ****' 
  
      cosa=COS(a) 
      sina=SIN(a) 
      adot2=adot**2 
  
      i1=k 
      IF(i1.GT.3) i1=i1-3 
      i2=i1+1 
      IF(i2.GT.3) i2=i2-3 
      i3=i2+1 
      IF(i3.GT.3) i3=i3-3 
  
      r(i1,i1) =  0.d0 
      r(i1,i2) =  0.d0 
      r(i1,i3) =  0.d0 
      r(i2,i1) =  0.d0 
      r(i2,i2) = -cosa*adot2 -sina*addot 
      r(i2,i3) = -sina*adot2 +cosa*addot 
      r(i3,i1) =  0.d0 
      r(i3,i2) =  sina*adot2 -cosa*addot 
      r(i3,i3) = -cosa*adot2 -sina*addot 
  
      END SUBROUTINE rotmt2         

END MODULE iers_ser

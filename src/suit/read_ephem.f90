!============library read_ephem====================
! read JPL ephemerides, using jpl_ephem
! converted for f90 by orbFit consortium
! version for asteroids by Carpino
! PUBLIC ROUTINES                                                       
! earth  orbital elelemts                                               
! earcar cartesian                                                      
! rdbep		read and interpolate binary ephemeris files                    
! jpllis	get JPL masses and IDs for a list of planets                   
! trange        sets the time limits of JPL ephem available, also sets
!                  fundamental-constants.mod  
!
!  note: use jpllis in masjpl (otherwise unstable for changes in DExxx) 
!   
! ==========================================                            
! EARTH                                                                 
! get Earth elements                                                    
SUBROUTINE earth(t0,eqp) 
  USE fund_const
  USE planet_masses 
  IMPLICIT NONE 
! input: epoch time                                                     
  DOUBLE PRECISION t0 
! output: elements of Earth (equinoctal, ecliptic)                      
  DOUBLE PRECISION eqp(6) 
! =============JPL EPHEM===============                                 
! data for masses                                                       
  INCLUDE 'jplhdr.h90' 
! output of JPL routine, Julian date, rotation matrix                   
  DOUBLE PRECISION et(2),rrd(6),xea(6),enne 
! integers for call to JPl routines                                     
  INTEGER ntarg,ncent,istate
! temporary arrays with planet names, pointers
  DOUBLE PRECISION gmp(13)
! fail flag for call to jpllis
  LOGICAL fail 
! ====================================                                  
! JPL Earth vector at observation time                                  
  et(1)=2400000.5d0 
  et(2)=t0 
  ntarg=3 
  ncent=11 
! first read header of current JPL ephemerides, storing masses in array gmp
  fail=.false.
  CALL jpllis2(gmp,fail) 
  IF(fail)THEN
     STOP ' **** earth: failed call to jpllis ********'
  ENDIF
! duplicate computation of gmse, in case masjpl has not been called yet 
  gmse=gmp(12)+gmp(3)
! ****** added on Sat Jun 14 1997 ******                                
! first istate need to be=2  (dpleph calculates also vel.)              
  istate=2 
  call dpleph(et,ntarg,ncent,rrd,istate) 
! Change of reference system EQUM00 ---> ECLM00                         
!      call rotpn(rot,'EQUM','J2000',0.d0,'ECLM','J2000',0.d0) 
  xea(1:3)=MATMUL(roteqec,rrd(1:3)) 
  xea(4:6)=MATMUL(roteqec,rrd(4:6)) 
! elements of Earth                                                     
  call coocha(xea,'CAR',gmse,eqp,'EQU',enne) 
END SUBROUTINE earth
! ======================================================================
! EARCAR - get Earth cartesian coordinates (ecliptic J2000)             
! ======================================================================
SUBROUTINE earcar(t0,xea,ifla) 
  USE fund_const
  USE planet_masses
  IMPLICIT NONE 
! input: epoch time, flag for getting Earth (heliocentric; ifla=1)      
!        or Sun (barycentric; ifla=2)                                   
  DOUBLE PRECISION, INTENT(IN) :: t0 
  INTEGER, INTENT(IN) :: ifla 
! output: heliocentric state vector of Earth (equinoctal, ecliptic)     
!      or barycentric state vector of Sun (equinoctal, ecliptic)        
  DOUBLE PRECISION, INTENT(OUT) :: xea(6) 
! =============JPL EPHEM===============                                 
! data for masses                                                       
  INCLUDE 'jplhdr.h90' 
! output of JPL routine, Julian date, rotation matrix                   
  double precision et(2),rrd(6) 
! integers for call to JPl routines                                     
  integer ntarg,ncent,istate 
! temporary arrays with planet names, pointers
  DOUBLE PRECISION gmp(13)
! fail flag for call to jpllis
  LOGICAL fail 
! ====================================                                  
! JPL Earth vector at observation time                                  
  et(1)=2400000.5d0 
  et(2)=t0 
  if (ifla.eq.1) then 
     ntarg=3 
     ncent=11 
  else 
     ntarg=11 
     ncent=12 
  endif
! first read header of current JPL ephemerides, storing masses in array gmp
  fail=.false.
  CALL jpllis2(gmp,fail) 
  IF(fail)THEN
     STOP ' **** earth: failed call to jpllis ********'
  ENDIF
! duplicate computation of gmse, in case masjpl has not been called yet 
  gmse=gmp(12)+gmp(3)
! ****** added on Sat Jun 14 1997 ******                                
! first istate need to be=2  (dpleph calculates also vel.)              
  istate=2 
  call dpleph(et,ntarg,ncent,rrd,istate) 
! Change of reference system EQUM00 ---> ECLM00                         
!      call rotpn(rot,'EQUM','J2000',0.d0,'ECLM','J2000',0.d0) 
  xea(1:3)=MATMUL(roteqec,rrd(1:3)) 
  xea(4:6)=MATMUL(roteqec,rrd(4:6)) 
END SUBROUTINE earcar
! Copyright (C) 1997 by Mario Carpino (carpino@brera.mi.astro.it)       
! Version: March 20, 1997                                               
! --------------------------------------------------------------------- 
!                                                                       
!  *****************************************************************    
!  *                                                               *    
!  *                          R D B E P                            *    
!  *                                                               *    
!  *         Reads and interpolates binary ephemeris files         *    
!  *                                                               *    
!  *****************************************************************    
!                                                                       
! INPUT:    T         -  Time (MJD, TDT)                                
!           NB        -  Number of bodies                               
!           ID        -  Identity (order number) of requested bodies    
!                                                                       
! OUTPUT:   X,V       -  Position and velocity vectors of requested     
!                        bodies (refsys: ECLM J2000)                    
!                                                                       
! EXAMPLE:                                                              
!     NB=2, ID(1)=3, ID(2)=1 means that the third and first bodies      
!           contained in the binary ephemeris files are requested       
!     on output, X(i,1) will contain the position vector of body        
!     number 3 and X(i,2) the position vector of body number 1          
!                                                                       
SUBROUTINE rdbep(t,nb,id,x,v) 
  USE planet_masses
  IMPLICIT NONE 
  INTEGER nb 
  INTEGER id(nb) 
  DOUBLE PRECISION t,x(3,nb),v(3,nb) 
  INTEGER unit,nhl,vers,n,recl,i,k,ib,lf 
  DOUBLE PRECISION t1,t2,dt,gma(nbepx),gma1(nbepx) 
  DOUBLE PRECISION c1,c2 
  DOUBLE PRECISION xv1(6),xv2(6),enne 
  LOGICAL first 
  DOUBLE PRECISION bomb(1)
  SAVE first,unit,nhl,n,t1,t2,dt,gma,gma1 
! Buffer                                                                
  INTEGER ipt1,ipt2 
  DOUBLE PRECISION tb1,tb2,elem1(6,nbepx),elem2(6,nbepx) 
  DOUBLE PRECISION enne1(nbepx),enne2(nbepx),m01(nbepx),m02(nbepx)   
  SAVE ipt1,ipt2,tb1,tb2,elem1,elem2,enne1,enne2,m01,m02 
  INTEGER lench 
  EXTERNAL lench 
  INTEGER ivel 
  DATA first/.true./ 
  IF(first) THEN 
     CALL filass(unit,filbep) 
! Open with a small record length (just to read the number of bodies)   
     OPEN(unit,FILE=filbep,ACCESS='DIRECT',RECL=24,STATUS='OLD') 
     READ(unit,REC=1) vers,nbep,n 
     IF(vers.NE.102) THEN 
        lf=lench(filbep) 
        WRITE(*,120) filbep(1:lf),vers 
        STOP '**** rdbep: abnormal end ****' 
     END IF
     CLOSE(unit) 
     CALL chkpdf(nbep,nbepx,'nbepx') 
     IF(n.LT.2) STOP '**** rdbep: n < 2 ****' 
     nhl=2+2*nbep 
! Open again with the correct record length                             
     recl=8*6*nbep 
     OPEN(unit,FILE=filbep,ACCESS='DIRECT',RECL=recl,STATUS='OLD') 
     READ(unit,REC=2) t1,t2,dt 
! Read masses                                                           
     DO  i=1,nbep 
       READ(unit,REC=2+i) masbep(i),gma(i),gma1(i)
     ENDDO
! Initialization of buffer                                              
     ipt1=1 
     ipt2=2 
     READ(unit,REC=nhl+ipt1) ((elem1(i,k),i=1,6),k=1,nbep) 
     READ(unit,REC=nhl+ipt2) ((elem2(i,k),i=1,6),k=1,nbep) 
     tb1=t1+dt*(ipt1-1) 
     tb2=t1+dt*(ipt2-1) 
     DO k=1,nbep 
        enne1(k)=SQRT(gma1(k)/elem1(1,k)**3) 
        enne2(k)=SQRT(gma1(k)/elem2(1,k)**3) 
        m01(k)=elem1(6,k) 
        m02(k)=elem2(6,k) 
     ENDDO 
     first=.false. 
  END IF
120 FORMAT(' ERROR: file ',A,' has unsupported version',I4) 
  IF(t.LT.tb1 .OR. t.GT.tb2) THEN 
     IF(t.LT.t1 .OR. t.GT.t2)THEN
        WRITE(*,*)' rdbep: for asteroid ephemerides from file ', filbep 
        WRITE(*,*)' time requested', t,' interval available',t1,t2
        i=2
        WRITE(*,*)bomb(i)           
        STOP '**** rdbep: T is out of bounds ****'                
     ENDIF
     ipt1=(t-t1)/dt+1 
     ipt2=ipt1+1 
     IF(ipt1.LT.1) STOP '**** rdbep: internal error (01) ****' 
     IF(ipt2.GT.n) STOP '**** rdbep: internal error (02) ****' 
     READ(unit,REC=nhl+ipt1) ((elem1(i,k),i=1,6),k=1,nbep) 
     READ(unit,REC=nhl+ipt2) ((elem2(i,k),i=1,6),k=1,nbep) 
     tb1=t1+dt*(ipt1-1) 
     tb2=t1+dt*(ipt2-1) 
     IF(t.LT.tb1 .OR. t.GT.tb2)                                    &
          &        STOP '**** rdbep: internal error (03) ****'               
     DO k=1,nbep 
        enne1(k)=SQRT(gma1(k)/elem1(1,k)**3) 
        enne2(k)=SQRT(gma1(k)/elem2(1,k)**3) 
        m01(k)=elem1(6,k) 
        m02(k)=elem2(6,k) 
     ENDDO 
  END IF
! Coefficients for interpolation                                        
  c1=(tb2-t)/dt 
  c2=(t-tb1)/dt 
  DO 2 ib=1,nb 
     k=id(ib) 
     elem1(6,k)=m01(k)+enne1(k)*(t-tb1) 
     elem2(6,k)=m02(k)+enne2(k)*(t-tb2) 
     ivel=1 
     CALL coocha(elem1(1,k),'KEP',gma1(k),xv1,'CAR',enne)              
     CALL coocha(elem2(1,k),'KEP',gma1(k),xv2,'CAR',enne)              
     DO i=1,3 
        x(i,ib)=c1*xv1(i)+c2*xv2(i) 
        v(i,ib)=c1*xv1(i+3)+c2*xv2(i+3) 
     END DO
2 END DO
                                                                        
END SUBROUTINE rdbep

!  *****************************************************************    
!  *                                                               *    
!  *                         J P L L I S 2                         *    
!  *                                                               *    
!  *      Get the list of masses and IDs from JPL DE header        *
!  *      Version robust w.r. to changes in the list cnam          *    
!  *                                                               *    
!  *****************************************************************    
!
! OUTPUT:   GMP       -  G*M(planet)                                    
!           FAIL      -  Error flag  
SUBROUTINE jpllis2(gmp,fail) 
  IMPLICIT NONE 
  DOUBLE PRECISION, INTENT(OUT) :: gmp(13) 
  LOGICAL, INTENT(OUT) :: fail
! for call to xstate 
  DOUBLE PRECISION et2(2),pv(6,12),pnut(4) 
  INTEGER list(12)
! loop index
  INTEGER i 
! JPLDE header 
  INCLUDE 'jplhdr.h90'
  INTEGER lench 
  EXTERNAL lench 
  DATA et2/2*0.d0/ 
  DATA list/12*0/
! begin execution 
  fail=.false.
  gmp=0.d0
! Dummy call to STATE for reading JPLDE header                          
  CALL state(et2,list,pv,pnut,1) 
!   'MERCURY','VENUS','EARTH_MOON','MARS','JUPITER','SATURN','URANUS','NEPTUNE','PLUTO','MOON','EARTH','SUN','EMRAT'/    
! loop on the cnam array of constant names
  DO i=1,30 ! maybe more???
     IF(cnam(i).eq.'GM1')THEN
        gmp(1)=cval(i)
     ELSEIF(cnam(i).eq.'GM2')THEN
        gmp(2)=cval(i)
     ELSEIF(cnam(i).eq.'GMB')THEN
        gmp(3)=cval(i)
     ELSEIF(cnam(i).eq.'GM4')THEN
        gmp(4)=cval(i)
     ELSEIF(cnam(i).eq.'GM5')THEN
        gmp(5)=cval(i)
     ELSEIF(cnam(i).eq.'GM6')THEN
        gmp(6)=cval(i)
     ELSEIF(cnam(i).eq.'GM7')THEN
        gmp(7)=cval(i)
     ELSEIF(cnam(i).eq.'GM8')THEN
        gmp(8)=cval(i)
     ELSEIF(cnam(i).eq.'GM9')THEN
        gmp(9)=cval(i) ! actually Pluto's mass is not anymore used
     ELSEIF(cnam(i).eq.'GMS')THEN
        gmp(12)=cval(i)
     ELSEIF(cnam(i).eq.'EMRAT')THEN
        gmp(13)=cval(i)
     ENDIF
  ENDDO
  gmp(10)=gmp(3)/(1+gmp(13))
  gmp(11)=gmp(3)*gmp(13)/(1+gmp(13)) 
  DO i=1,13
    IF(gmp(i).eq.0.d0) fail=.true.
  ENDDO
END SUBROUTINE jpllis2

! Copyright (C) 1997 by Mario Carpino (carpino@brera.mi.astro.it)       
! Version: October 13, 1997                                             
! --------------------------------------------------------------------- 
!                                                                       
!  *****************************************************************    
!  *                                                               *    
!  *                         J P L L I S                           *    
!  *                                                               *    
!  *      Get the list of masses and IDs from JPL DE header        *    
!  *                                                               *    
!  *****************************************************************    
!                                                                       
! INPUT:    NAMES     -  Planet names                                   
!           N         -  Number of planets                              
!                                                                       
! OUTPUT:   GMP       -  G*M(planet)                                    
!           ID        -  Planet ID number                               
!           FAIL      -  Error flag                                     
!                                                                       
SUBROUTINE jpllis(names,n,gmp,id,fail) 
  IMPLICIT NONE 
  INTEGER n 
  CHARACTER*(*) names(n) 
  DOUBLE PRECISION gmp(n) 
  INTEGER id(n) 
  LOGICAL fail 
  DOUBLE PRECISION et2(2),pv(6,12),pnut(4) 
  INTEGER list(12),i,ln,emopt 
  CHARACTER*12 n1 
  LOGICAL emwarn 
! JPLDE header 
  INCLUDE 'jplhdr.h90'
  INTEGER lench 
  EXTERNAL lench 
  DATA et2/2*0.d0/ 
  DATA list/12*0/ 
! Dummy call to STATE for reading JPLDE header                          
  CALL state(et2,list,pv,pnut,1) 
! Choice for Earth/Moon: 1=barycenter, 2=distinct bodies                
  emopt=0 
  emwarn=.false. 
  DO 1 i=1,n 
     n1=names(i) 
     CALL upcase(n1) 
     IF(n1.EQ.'MERCURY') THEN 
        id(i)=1 
        gmp(i)=cval(9) 
        IF(cnam(9).NE.'GM1')                                          &
 &        STOP '**** jpllis: internal error (C01) ****'             
     ELSEIF(n1.EQ.'VENUS') THEN 
        id(i)=2 
        gmp(i)=cval(10) 
        IF(cnam(10).NE.'GM2')                                         &
     &        STOP '**** jpllis: internal error (C02) ****'             
     ELSEIF(n1.EQ.'EARTH') THEN 
        id(i)=3 
        gmp(i)=cval(11)*emrat/(1+emrat) 
        IF(cnam(11).NE.'GMB')                                         &
     &        STOP '**** jpllis: internal error (C03) ****'             
        IF(emopt.NE.0 .AND. emopt.NE.2) emwarn=.true. 
        emopt=2 
     ELSEIF(n1.EQ.'MOON') THEN 
        id(i)=10 
        gmp(i)=cval(11)/(1+emrat) 
        IF(cnam(11).NE.'GMB')                                         &
     &        STOP '**** jpllis: internal error (C10) ****'             
        IF(emopt.NE.0 .AND. emopt.NE.2) emwarn=.true. 
        emopt=2 
     ELSEIF(n1.EQ.'EARTH+MOON') THEN 
        id(i)=13 
        gmp(i)=cval(11) 
        IF(cnam(11).NE.'GMB')                                         &
     &        STOP '**** jpllis: internal error (C13) ****'             
        IF(emopt.NE.0 .AND. emopt.NE.1) emwarn=.true. 
        emopt=1 
     ELSEIF(n1.EQ.'MARS') THEN 
        id(i)=4 
        gmp(i)=cval(12) 
        IF(cnam(12).NE.'GM4')                                         &
     &        STOP '**** jpllis: internal error (C04) ****'             
     ELSEIF(n1.EQ.'JUPITER') THEN 
        id(i)=5 
        gmp(i)=cval(13) 
        IF(cnam(13).NE.'GM5')                                         &
     &        STOP '**** jpllis: internal error (C05) ****'             
     ELSEIF(n1.EQ.'SATURN') THEN 
        id(i)=6 
        gmp(i)=cval(14) 
        IF(cnam(14).NE.'GM6')                                         &
     &        STOP '**** jpllis: internal error (C06) ****'             
     ELSEIF(n1.EQ.'URANUS') THEN 
        id(i)=7 
        gmp(i)=cval(15) 
        IF(cnam(15).NE.'GM7')                                         &
     &        STOP '**** jpllis: internal error (C07) ****'             
     ELSEIF(n1.EQ.'NEPTUNE') THEN 
        id(i)=8 
        gmp(i)=cval(16) 
        IF(cnam(16).NE.'GM8')                                         &
     &        STOP '**** jpllis: internal error (C08) ****'             
     ELSEIF(n1.EQ.'PLUTO') THEN 
        id(i)=9 
        gmp(i)=cval(17) 
        IF(cnam(17).NE.'GM9')                                         &
     &        STOP '**** jpllis: internal error (C09) ****'             
     ELSE 
        ln=lench(names(i)) 
        WRITE(*,100) names(i)(1:ln) 
        fail=.true. 
     END IF
100  FORMAT(' ERROR: ',A,' is unknown among JPL planets') 
1 END DO
  IF(emwarn) THEN 
     WRITE(*,101) 
     fail=.true. 
  END IF
101 FORMAT(' ERROR in the list of JPL planets: please DO NOT select'/ &
     &       '       "Earth+Moon" (Earth-Moon barycenter) together'/    &
     &       '       with "Earth" and/or "Moon"')                       
                                                                        
END SUBROUTINE jpllis
!                                                                       
subroutine trange
! for vlight, aukm, emratio
  USE fund_const 
  implicit none 
!                                                                       
  include 'timespan.h90' 
!                                                                       
  double precision et2(2),pv(6,12),pnut(4) 
  integer list(12) 
  double precision tt,deltt 
!                                                                       
! JPL  header                                                           
  include 'jplhdr.h90' 
!                                                                       
  data et2/2*0.d0/ 
  data list/12*0/ 
!                                                                       
! Dummy call to STATE for reading JPLDE header                          
  CALL state(et2,list,pv,pnut,1) 
! store in common timespan time span of jpleph                          
! transformation from JD to MJD                                         
  tejpl1=ss(1)-2400000.5d0 
  tejpl2=ss(2)-2400000.5d0 
!      write(*,*) tejpl1,tejpl2                                         
!                                                                       
! Dummy call to deltt to read the ET-UT data                            
  tt=deltt(50000.d0) 
! speed of light (IAU 1976) in km/s                                     
  ckm=299792.458d0 
! id in km/day                                                          
  ckm=ckm*8.64d4 
! conversion to au/day                                                  
  vlight=ckm/au 
! au, earth/mooon mass ratio moved to fundamental_constants.mod
  aukm=au
  reau=eradkm/aukm
  emratio=emrat
! add dummy call to rdbep to setup the range of values for that too
END subroutine trange

! ======================================================================
! MOONCAR - get Moon cartesian coordinates (ecliptic J2000)             
! FB October 2008
! ======================================================================
SUBROUTINE mooncar(t0,xmoon,ifla) 
  USE fund_const
  USE planet_masses
  IMPLICIT NONE 
! input: epoch time, flag for getting Moon (heliocentric; ifla=1)      
!        or Sun (barycentric; ifla=2)                                   
  DOUBLE PRECISION, INTENT(IN) :: t0 
  INTEGER, INTENT(IN) :: ifla 
! output: heliocentric state vector of Moon (equinoctal, ecliptic)     
!         or barycentric state vector of Sun (equinoctal, ecliptic)        
  DOUBLE PRECISION, INTENT(OUT) :: xmoon(6) 
! =============JPL EPHEM===============                                 
! data for masses                                                       
  INCLUDE 'jplhdr.h90' 
! output of JPL routine, Julian date, rotation matrix                   
  double precision et(2),rrd(6) 
! integers for call to JPl routines                                     
  integer ntarg,ncent,istate 
! temporary arrays with planet names, pointers
  DOUBLE PRECISION gmp(13)
! fail flag for call to jpllis
  LOGICAL fail 
! ====================================                                  
! JPL Moon vector at observation time                                  
  et(1)=2400000.5d0 
  et(2)=t0 
  if (ifla.eq.1) then 
     ntarg=10 
     ncent=11 
  else 
     ntarg=11 
     ncent=12 
  endif
! first read header of current JPL ephemerides, storing masses in array gmp
  fail=.false.
  CALL jpllis2(gmp,fail) 
  IF(fail)THEN
     STOP ' **** earth: failed call to jpllis ********'
  ENDIF
! duplicate computation of gmse, in case masjpl has not been called yet 
  gmse=gmp(12)+gmp(3)
! ****** added on Sat Jun 14 1997 ******                                
! first istate need to be=2  (dpleph calculates also vel.)              
  istate=2 
  call dpleph(et,ntarg,ncent,rrd,istate) 
! Change of reference system EQUM00 ---> ECLM00                         
!      call rotpn(rot,'EQUM','J2000',0.d0,'ECLM','J2000',0.d0) 
  xmoon(1:3)=MATMUL(roteqec,rrd(1:3)) 
  xmoon(4:6)=MATMUL(roteqec,rrd(4:6)) 
END SUBROUTINE mooncar

! ======================================================================
! SUNMOON_ CAR - get Sun+Moon cartesian coordinates 
!     either ecliptic or equatorial J2000 
! GT November 2008
! ======================================================================
SUBROUTINE sunmoon_car(t0,equ,xsun,xmoon) 
  USE fund_const
  USE planet_masses
  IMPLICIT NONE 
! SUNMOON_CAR - get geocentric cartesian coordinates
! for Sun and Moon, equatorial (equ=TRUE) or ecpliptic (equ=FALSE) 
! input: epoch time 
  DOUBLE PRECISION, INTENT(IN)                          :: t0
! equatorial or ecliptic
  LOGICAL, INTENT(IN)                                   :: equ 
! output: geocentric state vector of Sun and Moon (ecliptic)      
  DOUBLE PRECISION, DIMENSION(6), INTENT(OUT)           :: xsun
  DOUBLE PRECISION, DIMENSION(6), INTENT(OUT), OPTIONAL :: xmoon 
! =============JPL EPHEM===============                                 
! data for masses                                                       
  INCLUDE 'jplhdr.h90' 
! output of JPL routine, Julian date, rotation matrix                   
  DOUBLE PRECISION ::  et(2),rrd(6) 
! integers for call to JPl routines                                     
  INTEGER :: ntarg,ncent,istate 
! ====================================                                  
! JPL Earth vector at observation time                                  
  et(1)=2400000.5d0 
  et(2)=t0 
  ncent=3 
! first istate need to be=2  (dpleph calculates also vel.)              
  istate=2
  ntarg=11 
  CALL dpleph(et,ntarg,ncent,rrd,istate) 
  IF(equ)THEN
     xsun=rrd
  ELSE
! Change of reference system EQUM00 ---> ECLM00                         
     xsun(1:3)=MATMUL(roteqec,rrd(1:3)) 
     xsun(4:6)=MATMUL(roteqec,rrd(4:6))
  ENDIF
  IF(PRESENT(xmoon))THEN
    ntarg=10
    CALL dpleph(et,ntarg,ncent,rrd,istate)
    IF(equ)THEN
       xmoon=rrd
    ELSE   
    ! Change of reference system EQUM00 ---> ECLM00                         
       xmoon(1:3)=MATMUL(roteqec,rrd(1:3)) 
       xmoon(4:6)=MATMUL(roteqec,rrd(4:6)) 
    ENDIF
  ENDIF
END SUBROUTINE sunmoon_car

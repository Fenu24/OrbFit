! MODULE station_coordinates
! including Tranformations from and to geodetic coordinates

! SUBROUTINE geodetic_to_cartesian
!      transformation of geodetic coordinates into cartesian
! SUBROUTINE cartesian_to_geodetic
!      transformation of cartesian coordinates into geodetic

MODULE station_coordinates
USE output_control
IMPLICIT NONE
PRIVATE

! Internal PARAMETERs
DOUBLE PRECISION, PARAMETER  :: earth_radius = 6.3781363d6                           ! reference equatorial radius of the Earth (m)
DOUBLE PRECISION, PARAMETER  :: reciprocal_flattening = 298.257d0                    ! 1/f
DOUBLE PRECISION, PARAMETER  :: flattening = 1.d0/reciprocal_flattening              ! flattening f
DOUBLE PRECISION, PARAMETER  :: eccentricity_squared = flattening*(2.d0-flattening)  ! e2 = f(2-f)

! LIST OF PUBLIC ENTITIES
! SUBROUTINEs
PUBLIC :: geodetic_to_cartesian, cartesian_to_geodetic, obscoo, codestat, statcode

CONTAINS

! Transformation of geodetic coordinates into cartesian
SUBROUTINE geodetic_to_cartesian(longitude,latitude,altitude,x)

DOUBLE PRECISION,               INTENT(IN)  :: longitude   ! geodetic longitude (positive east, rad)
DOUBLE PRECISION,               INTENT(IN)  :: latitude    ! geodetic latitude (rad)
DOUBLE PRECISION,               INTENT(IN)  :: altitude    ! altitude on ellipsoid (m)
DOUBLE PRECISION, DIMENSION(3), INTENT(OUT) :: x           ! cartesian coordinates (m)

DOUBLE PRECISION :: enne

enne=earth_radius/SQRT(1.d0-eccentricity_squared*SIN(latitude)**2)
x(1)=(enne+altitude)*COS(latitude)*COS(longitude)
x(2)=(enne+altitude)*COS(latitude)*SIN(longitude)
x(3)=(enne*(1.d0-eccentricity_squared)+altitude)*SIN(latitude)

END SUBROUTINE geodetic_to_cartesian

! Transformation of cartesian coordinates into geodetic
SUBROUTINE cartesian_to_geodetic(x,longitude,latitude,altitude,vertical)
USE fund_const

DOUBLE PRECISION, DIMENSION(3), INTENT(IN)            :: x           ! cartesian coordinates (m)
DOUBLE PRECISION,               INTENT(OUT)           :: longitude   ! geodetic longitude (positive east, rad)
DOUBLE PRECISION,               INTENT(OUT)           :: latitude    ! geodetic latitude (rad)
DOUBLE PRECISION,               INTENT(OUT)           :: altitude    ! altitude on ellipsoid (m)
DOUBLE PRECISION, DIMENSION(3), INTENT(OUT), OPTIONAL :: vertical    ! local vertical direction

DOUBLE PRECISION :: t,h,h0,dn,dnh,z,zt,delth,sinfi,cosfi
INTEGER          :: i

t=0.d0
h0=0.d0
dn=0.d0
z=x(3)
i=1

10 CONTINUE
zt=z+t
dnh=SQRT((x(1)**2+x(2)**2+zt**2))
sinfi=zt/dnh
dn=earth_radius/SQRT(1.d0-eccentricity_squared*sinfi**2)
t=dn*eccentricity_squared*sinfi
h=dnh-dn
delth=ABS(h-h0)
IF(delth > 1.d-2) THEN
   i=i+1
   h0=h
   IF(i < 40) GOTO 10
   WRITE(ierrou,*) 'cartesian_to_geodetic: not converged; DELTH= ',delth
END IF
altitude=h
IF(sinfi.gt.1.d0)THEN
   write(*,*)'sinfi=',sinfi
   IF(sinfi.lt.1.d0+10.d0*epsilon(1.d0))THEN
      sinfi=1.d0
   ELSE
      write(ierrou,*)'oss_dif2: error! sinfi=',sinfi
      numerr=numerr+1
   ENDIF
ENDIF
IF(sinfi.lt.-1.d0)THEN
   write(*,*)'sinfi=',sinfi
   IF(sinfi.gt.-1.d0-10.d0*epsilon(1.d0))THEN
      sinfi=-1.d0
   ELSE
      write(ierrou,*)'oss_dif2: error! sinfi=',sinfi 
      numerr=numerr+1
   ENDIF
ENDIF
latitude=ASIN(sinfi)
longitude=ATAN2(x(2),x(1))
IF(longitude < 0.d0) longitude=dpig+longitude
IF(PRESENT(vertical)) THEN
   cosfi=COS(latitude)
   vertical(1)=cosfi*dcos(longitude)
   vertical(2)=cosfi*dsin(longitude)
   vertical(3)=sinfi
END IF

END SUBROUTINE cartesian_to_geodetic

! Copyright (C) 1997 by Mario Carpino (carpino@brera.mi.astro.it)
! Version: February 24, 1997
!
!  ***************************************************************
!  *                                                             *
!  *                         O B S C O O                         *
!  *                                                             *
!  *          Body-fixed coordinates of an observatory           *
!  *                                                             *
!  ***************************************************************
!
!
! INPUT:    IDOBS     -  Identifier of the observatory (0-999)
!
! OUTPUT:   XBF(3)    -  Body-fixed position of the observatory (au)
!           NAME      -  Name of the observatory
!           SHIFTST   -  Fraction of day after GMT
!
SUBROUTINE obscoo(idobs,xbf,name, shiftst)
  USE fund_const
! NEEDED common blocks:                                                 
  INCLUDE 'comlib.h90'
  integer lnobnx
  parameter (lnobnx=47)
! INPUT
  integer, intent(in):: idobs
! OUTPUT
  double precision, intent(out) :: xbf(3)
  character*(*), intent(out) :: name
! OPTIONAL
  double precision, intent(out), optional :: shiftst
! END INTERFACE
  integer ns1,ns2
  parameter (ns1=0,ns2=6200)
  double precision al1,pxy1,pz1,xbfv(3,ns1:ns2), princ, shiftst_ar(ns1:ns2)
  integer unit,i,k
  character*3 ocod
  character*(lnobnx) name1,namev(ns1:ns2)
  character*80 rec
  logical first,loaded(ns1:ns2),found
  save first,loaded,xbfv,namev!,eradau
  data first/.true./
! begin execution
  if(first)then
! Get parallax data from MPC
     DO i=ns1,ns2
        loaded(i)=.false.
     ENDDO
     call filopl(unit,'OBSCODE.dat')
2    continue
     read(unit,110,end=3)ocod,al1,pxy1,pz1,name1
110  format(a3,f10.5,f8.6,f9.6,a)
     call statcode(ocod,k)
     if(k.lt.ns1.or.k.gt.ns2)stop ' **** obscoo: obscode out of range ****'
     al1=princ(al1*radeg)
     shiftst_ar(k)=al1/dpig
     xbfv(1,k)=reau*pxy1*cos(al1)
     xbfv(2,k)=reau*pxy1*sin(al1)
     xbfv(3,k)=reau*pz1
     namev(k)=name1
     loaded(k)=.true.
     goto 2
3    continue
     call filclo(unit,' ')
! parallax data for simulations
     INQUIRE(FILE=libdir(1:lenld)//'OBSCODE_sim.dat',EXIST=found)
     IF(found)THEN
        CALL filopl(unit,'OBSCODE_sim.dat')
22      continue
        read(unit,110,end=23)ocod,al1,pxy1,pz1,name1
        call statcode(ocod,k)
        if(k.lt.ns1.or.k.gt.ns2)stop ' **** obscoo: obscode out of range ****'
        al1=princ(al1*radeg)
        shiftst_ar(k)=al1/dpig
        xbfv(1,k)=reau*pxy1*cos(al1)
        xbfv(2,k)=reau*pxy1*sin(al1)
        xbfv(3,k)=reau*pz1
        namev(k)=name1
        loaded(k)=.true.
        goto 22
23      continue
        CALL filclo(unit,' ')
     ENDIF
!***********************************
! Get parallax data for radar sites
     call filopl(unit,'RADCODE.dat')
4    continue
     read(unit,'(a)',end=5)rec
     if(rec(1:1).eq.'!'.or.rec(1:1).eq.' ')goto 4
     read(rec,*)k,al1,pxy1,pz1,name1
     if(k.lt.ns1.or.k.gt.ns2)stop ' **** obscoo: internal error (11) ****'
     al1=al1*radeg
     xbfv(1,k)=1d-10*pxy1*cos(al1)
     xbfv(2,k)=1d-10*pxy1*sin(al1)
     xbfv(3,k)=1d-10*pz1
     namev(k)=name1
     loaded(k)=.true.
     goto 4
5    call filclo(unit,' ')
!***********************************
     first=.false.
  end if
  IF(PRESENT(shiftst)) THEN
     IF(loaded(idobs)) THEN
        shiftst=shiftst_ar(idobs)
     ENDIF
  ENDIF
  if(loaded(idobs)) then
     do i=1,3
        xbf(i)=xbfv(i,idobs)
     enddo
     name=namev(idobs)
  else
     write(*,101)idobs
101  format(' obscoo: observatory',i4,' is not listed in file "OBSCODE.dat"')
!          stop ' **** obscoo: abnormal end ****'
     write(ierrou,101)idobs
     numerr=numerr+1
     do i=1,3
        xbf(i)=0.d0
     enddo
     name='UNKNOWN'
  end if

END SUBROUTINE obscoo
! ==============================================                        
!  statcode, codestat                                                   
! conversion from/to exadecimal to/from numeric code for observing stati
! answer to mess done by MPC in April 2002                              
! ==============================================                        
      SUBROUTINE statcode(obsstr,iobs) 
      IMPLICIT NONE 
! input                                                                 
      CHARACTER*3 obsstr 
! output                                                                
      INTEGER iobs 
! end interface                                                         
      CHARACTER*1 alfanum 
      INTEGER hundreds, temp 
      LOGICAL isnum 
      READ(obsstr,100)iobs 
  100 FORMAT(1x,i2) 
      READ(obsstr,'(A1)')alfanum 
      IF(isnum(alfanum))THEN 
         READ(alfanum,'(I1)')hundreds 
         iobs=iobs+100*hundreds 
      ELSEIF(alfanum.eq.' ')THEN 
         iobs=iobs 
      ELSE 
         temp=ichar(alfanum) - 55 
         IF(temp.gt.35) temp = temp - 6 
         iobs=iobs+100*temp 
      ENDIF 
      RETURN 
      END SUBROUTINE statcode                                           
!================================================                       
      SUBROUTINE codestat(iobs,obsstr) 
      IMPLICIT NONE 
! input                                                                 
      INTEGER iobs 
! output                                                                
      CHARACTER*3 obsstr 
! end interface                                                         
      CHARACTER*1 alfanum 
      INTEGER hundreds, units 
      units=mod(iobs,100) 
      hundreds=(iobs-units)/100 
      IF(hundreds.le.9)THEN 
         WRITE(obsstr,'(I3.3)')iobs 
      ELSE 
         alfanum=char(55+hundreds) 
         WRITE(obsstr,101)alfanum,units 
  101    FORMAT(A1,I2.2) 
      ENDIF 
      RETURN 
      END SUBROUTINE codestat  

END MODULE station_coordinates





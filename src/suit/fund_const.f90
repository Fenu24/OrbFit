! ==================MODULE fund_const========================
! contains only data; should replace 4 headers:
!	../include/trig.h90 \
!	../include/jplhdr.h90 \ only for au in km and Earth-Moon mass ratio
!	../include/sunmass.h90  
!	../include/vlight.h90
MODULE fund_const
IMPLICIT NONE
PUBLIC !all the fundamental constants are public and available anywhere in the code,
! by including the USE fund_const statement
! Variable precision facility (double or quadruple precision)
INTEGER, PARAMETER :: dkind=KIND(1.0d0)
INTEGER, PARAMETER :: qkind=KIND(1.0q0)
! Copyright (C) 1997 by A.Milani and M.Carpino
! Version: October 8, 1997
! gravitational  constant: units, mass of the Sun, day, au
      double precision gk,gms,eradkm
      parameter (gk=0.01720209895d0)
      parameter (gms=gk*gk)
      parameter (eradkm=6.3781363d3)
! Trigonometric constants
DOUBLE PRECISION, PARAMETER :: dpig=6.28318530717958648d0    ! 2 pi
DOUBLE PRECISION, PARAMETER :: pig=dpig/2.0d0                ! pi
DOUBLE PRECISION, PARAMETER :: radeg=dpig/360.0d0            ! Radians from degrees
DOUBLE PRECISION, PARAMETER :: radsec=radeg/3600.0d0         ! Radians from arcseconds
DOUBLE PRECISION, PARAMETER :: degrad=360.0d0/dpig           ! Degrees from radians
DOUBLE PRECISION, PARAMETER :: secrad=3600.0d0/radeg         ! Arcseconds from radians
DOUBLE PRECISION, PARAMETER :: radh=dpig/24.0d0              ! Radians from hours
DOUBLE PRECISION, PARAMETER :: hrad=24.0d0/dpig              ! Hours from radians
DOUBLE PRECISION, PARAMETER :: r_sun=6.96d5                  ! Radius of Sun in Km
!the ones below need to be initialised from JPL ephemerides
!  parameters initialized in read_ephem/trange
      double precision vlight,ckm !speed of light in current units, in au/day and in km/s 
! part of JPL ephemerides header common; name needs to be changed, initialization
! also to be performed in trange aukm=au, emratio=emrat
      double precision aukm,emratio,reau ! au in km, Earth-Moon mass ratio, Earth radius in current units

! DATA: matrices of standard rotations
DOUBLE PRECISION, DIMENSION(3,3), PUBLIC :: roteqec,roteceq

! selection of right hand side
INTEGER, PUBLIC :: rhs ! 1=force_model 2=force_sat 3=orbit9

END MODULE fund_const

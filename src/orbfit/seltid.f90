! Copyright (C) 1997 by Mario Carpino (carpino@brera.mi.astro.it)
! Version: November 24, 1997
! ---------------------------------------------------------------------
!
!  *****************************************************************
!  *                                                               *
!  *                         S E L T I D                           *
!  *                                                               *
!  *   Selection of more appropriate times for initial conditions  *
!  *          of two object for orbit identification               *
!  *           (using average of initial conditions)               *
!  *                                                               *
!  *****************************************************************
!
! INPUT:    TDT1      -  Time (MJD, TDT)           \
!           RMSA1     -  A-priori RMS of RA (rad)   |
!           RMSD1     -  A-priori RMS of DEC (rad)  |  object 1
!           SEL1      -  Selection indicator        |
!           N1        -  Number of observations    /
!           TDT2      -  Time (MJD, TDT)           \
!           RMSA2     -  A-priori RMS of RA (rad)   |
!           RMSD2     -  A-priori RMS of DEC (rad)  |  object 2
!           SEL2      -  Selection indicator        |
!           N2        -  Number of observations    /
!
! OUTPUT:   TBE(2)    -  Best epochs for initial conditions (MJD, TDT)
!
      SUBROUTINE seltid(tdt1,rmsa1,rmsd1,sel1,n1,                       &
     &                  tdt2,rmsa2,rmsd2,sel2,n2,tbe)
      USE fund_const
      IMPLICIT NONE

      INTEGER n1,sel1(n1)
      INTEGER n2,sel2(n2)
      DOUBLE PRECISION tdt1(n1),rmsa1(n1),rmsd1(n1)
      DOUBLE PRECISION tdt2(n2),rmsa2(n2),rmsd2(n2)
      DOUBLE PRECISION tbe(2)

! TUNING PARAMETERS
! Max RMS of a "good" observation (arcsec)
      DOUBLE PRECISION maxrms
      PARAMETER (maxrms=5.d0)

      DOUBLE PRECISION rms1,ts(2),te(2)
      INTEGER i1,i2

      rms1=maxrms*radsec
      CALL ustsp(tdt1,rmsa1,rmsd1,sel1,n1,rms1,ts(1),te(1))
      CALL ustsp(tdt2,rmsa2,rmsd2,sel2,n2,rms1,ts(2),te(2))

! I1 is the index of the object whose observation timespan starts
! earlier
      i1=1
      IF(ts(2).LT.ts(1)) i1=2
      i2=3-i1

! Disjoint intervals
      IF(ts(i2).GE.te(i1)) THEN
          tbe(i1)=te(i1)
          tbe(i2)=ts(i2)
! Timespan I2 totally contained in timespan I1
      ELSEIF(te(i2).LE.te(i1)) THEN
          tbe(i1)=(ts(i2)+te(i2))/2
          tbe(i2)=tbe(i1)
! Timespans partially overlapping
      ELSE
          tbe(i1)=(ts(i2)+te(i1))/2
          tbe(i2)=tbe(i1)
      END IF

      END

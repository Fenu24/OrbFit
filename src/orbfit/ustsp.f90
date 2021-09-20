! Copyright (C) 1997 by Mario Carpino (carpino@brera.mi.astro.it)
! Version: November 24, 1997
! ---------------------------------------------------------------------
!
!  *****************************************************************
!  *                                                               *
!  *                          U S T S P                            *
!  *                                                               *
!  *              Useful timespan of observations                  *
!  *                                                               *
!  *****************************************************************
!
! INPUT:    T         -  Observation times (MJD, TDT)
!           RMSA      -  A-priori RMS of RA (rad)
!           RMSD      -  A-priori RMS of DEC (rad)
!           SEL       -  Selection indicator
!           N         -  Number of observations
!           RMS       -  Max RMS of a "good" observation (rad)
!           T1,T2     -  Time of first and last "good" observation (MJD,
!
      SUBROUTINE ustsp(t,rmsa,rmsd,sel,n,rms,t1,t2)
      IMPLICIT NONE

      INTEGER n,sel(n)
      DOUBLE PRECISION t(n),rmsa(n),rmsd(n)
      DOUBLE PRECISION rms,t1,t2

      INTEGER i
      LOGICAL prec,prec1,first

! Determine whether there are "precise" observations
      prec=.false.
      DO 1 i=1,n
      IF(sel(i).LE.0) GOTO 1
      prec1=((rmsa(i).LE.rms).AND.(rmsd(i).LE.rms))
      IF(prec1) THEN
          prec=.true.
          GOTO 2
      END IF
    1 END DO
    2 CONTINUE

! Determine useful timespan of observations
      first=.true.
      IF(prec) THEN
          DO 3 i=1,n
          IF(sel(i).LE.0) GOTO 3
          prec1=((rmsa(i).LE.rms).AND.(rmsd(i).LE.rms))
          IF(.NOT.prec1) GOTO 3
          IF(first) THEN
              t1=t(i)
              t2=t(i)
              first=.false.
          ELSE
              t1=MIN(t1,t(i))
              t2=MAX(t2,t(i))
          END IF
    3     CONTINUE
      ELSE
          DO 4 i=1,n
          IF(sel(i).LE.0) GOTO 4
          IF(first) THEN
              t1=t(i)
              t2=t(i)
              first=.false.
          ELSE
              t1=MIN(t1,t(i))
              t2=MAX(t2,t(i))
          END IF
    4     CONTINUE
      END IF

! Emergency action
      IF(first) THEN
          t1=t(1)
          t2=t(n)
      END IF

      END

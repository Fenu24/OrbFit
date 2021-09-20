! Copyright (C) 1997-1999 by Mario Carpino (carpino@brera.mi.astro.it)
! Version: February 9, 1999
! ---------------------------------------------------------------------
!
!  *****************************************************************
!  *                                                               *
!  *                         S R T O B S                           *
!  *                                                               *
!  *     Sorting of observations in order of increasing time       *
!  *                                                               *
!  *****************************************************************
!
! INPUT:    OBJID     -  Object identifier
!           TDT       -  Time (MJD, TDT)
!           TUTM      -  Time (MJD, UTM)
!           ALPHA     -  Right ascension (rad)
!           DELTA     -  Declination (rad)
!           IOBS      -  Observation type
!           SMAG      -  Measured magnitude (string)
!           OBSCOD    -  Observatory code
!           SEL       -  Selection index
!           RMSA      -  A-priori RMS of RA (rad)
!           RMSD      -  A-priori RMS of DEC (rad)
!           RMSMAG    -  A-priori RMS of magnitude
!           N         -  Number of observations
!
      SUBROUTINE srtobs(objid,tdt,tutm,alpha,delta,iobs,smag,obscod,    &
     &                  sel,rmsa,rmsd,rmsmag,n)
      IMPLICIT NONE

      INTEGER n
      INTEGER obscod(n),iobs(n),sel(n)
      DOUBLE PRECISION tdt(n),tutm(n),alpha(n),delta(n)
      DOUBLE PRECISION rmsa(n),rmsd(n),rmsmag(n)
      CHARACTER*(*) objid(n),smag(n)

      INTEGER i1,i2,itmp,nit
      CHARACTER tobjid*9,tsmag*6
      LOGICAL change
      DOUBLE PRECISION rtmp

      nit=0

    1 CONTINUE
      change=.false.
      nit=nit+1
      IF(nit.GT.n+2) STOP '**** srtobs: internal error (01) ****'

      DO 2 i1=1,n-1
      i2=i1+1
      IF(tdt(i1).GT.tdt(i2)) THEN
          tobjid=objid(i1)
          objid(i1)=objid(i2)
          objid(i2)=tobjid
          rtmp=tdt(i1)
          tdt(i1)=tdt(i2)
          tdt(i2)=rtmp
          rtmp=tutm(i1)
          tutm(i1)=tutm(i2)
          tutm(i2)=rtmp
          rtmp=alpha(i1)
          alpha(i1)=alpha(i2)
          alpha(i2)=rtmp
          rtmp=delta(i1)
          delta(i1)=delta(i2)
          delta(i2)=rtmp
          itmp=iobs(i1)
          iobs(i1)=iobs(i2)
          iobs(i2)=itmp
          tsmag=smag(i1)
          smag(i1)=smag(i2)
          smag(i2)=tsmag
          itmp=obscod(i1)
          obscod(i1)=obscod(i2)
          obscod(i2)=itmp
          itmp=sel(i1)
          sel(i1)=sel(i2)
          sel(i2)=itmp
          rtmp=rmsa(i1)
          rmsa(i1)=rmsa(i2)
          rmsa(i2)=rtmp
          rtmp=rmsd(i1)
          rmsd(i1)=rmsd(i2)
          rmsd(i2)=rtmp
          rtmp=rmsmag(i1)
          rmsmag(i1)=rmsmag(i2)
          rmsmag(i2)=rtmp
          change=.true.
      END IF
    2 END DO

      IF(change) GOTO 1

      END

!
!  *****************************************************************
!  *                                                               *
!  *                         O F P R O P                           *
!  *                                                               *
!  *                Auxiliary routine for ORBFIT:                  *
!  *               propagation of orbital elements                 *
!  *                                                               *
!  *****************************************************************
!
! INPUT:    UNIREP    -  FORTRAN unit for report
!           UNIELE    -  FORTRAN unit for output of orbital elements
!           OPELE     -  Flag (need to open UNIELE)
!           ELEOUT    -  Name of file for output of orbital elements
!           NAME      -  Object names
!           DEFORB    -  Orbit definition flag
!           DEFCN     -  Tells whether covariance/normal matrices
!                            are defined
!           ELEM      -  Orbital elements
!           ELEM_UNC  -  Orbital element uncertainty
!           MASS      -  Object masses
!           COMELE    -  Comment on orbital elements
!           NOBJ      -  Number of objects
!           GMSUN     -  G*Mass(Sun)
!           OEPTIM    -  Epoch of output elements (MJD, TDT)
!           OEPSET    -  Flag stating that an output epoch is requested
!           OETYPE    -  Type of output elements (CAR/EQU/KEP/EQP)
!
      SUBROUTINE ofprop(unirep,uniele,opele,eleout,name,deforb,defcn,elem,elem_unc,mass,comele,nobj,gmsun,oeptim,oepset,oetype)
      USE orbit_elements
      USE propag_state
      IMPLICIT NONE

      INTEGER unirep,uniele,nobj
      TYPE(orbit_elem),DIMENSION(nobj) :: elem
      TYPE(orb_uncert),DIMENSION(nobj) :: elem_unc
      DOUBLE PRECISION gmsun
      DOUBLE PRECISION mass(nobj),oeptim
      LOGICAL opele,deforb(nobj),defcn(nobj),oepset
      CHARACTER*(*) eleout,name(nobj),comele(nobj),oetype

      INCLUDE 'parcmc.h90'

      INTEGER i,ln,lc,fail
      DOUBLE PRECISION gma1,enne,dxde(6,6)
      DOUBLE PRECISION elemt(6),covet(6,6),noret(6,6)
      DOUBLE PRECISION elemo(6),coveo(6,6),noreo(6,6)
      LOGICAL error,defnro
      TYPE(orbit_elem) :: elem1,elem2
      TYPE(orb_uncert) :: unc1,unc2

      INTEGER lench
      EXTERNAL lench

      IF(.NOT.oepset) RETURN
      WRITE(unirep,100)
  100 FORMAT('Propagation of orbital elements:')

      DO 1 i=1,nobj
      IF(.NOT.deforb(i)) GOTO 1
      gma1=gmsun*(1+mass(i))

      IF(defcn(i)) THEN
         CALL pro_ele(elem(i),oeptim,elem1,elem_unc(i),unc1)
         CALL coo_cha(elem1,oetype,elem2,fail,dxde)
!===     IF(fail /= 0) ???
         CALL convertunc(unc1,dxde,unc2)
      ELSE
         CALL pro_ele(elem(i),oeptim,elem1)
         CALL coo_cha(elem1,oetype,elem2,fail)
      END IF
      ln=lench(name(i))
      lc=lench(comele(i))
      IF(lc.LE.0) THEN
          WRITE(unirep,101) name(i)(1:ln)
      ELSE
          WRITE(unirep,102) name(i)(1:ln),comele(i)(1:lc)
      END IF
      CALL outele(unirep,elem2%coord,elem2%coo,elem2%t,' ',.true.,.false.)
      IF(opele) THEN
          OPEN(uniele,FILE=eleout,STATUS='UNKNOWN')
          opele=.false.
          WRITE(uniele,301) comcha
          CALL wromlh(uniele,'ECLM','J2000')
      END IF
      IF(lc.LE.0) THEN
          WRITE(uniele,103) comcha,name(i)(1:ln)
      ELSE
          WRITE(uniele,104) comcha,name(i)(1:ln),comele(i)(1:lc)
      END IF
      CALL wromlr(uniele,name(i),elem2%coord,elem2%coo,elem2%t,unc2%g,defcn(i),unc2%c,defcn(i),elem2%h_mag,elem2%g_mag,mass(i),6)
    1 END DO
  101 FORMAT(5X,'Propagated orbital elements for object ',A,':')
  102 FORMAT(5X,'Propagated orbital elements for object ',A,' (',A,'):')
  103 FORMAT(A,' Orbital elements of object ',A)
  104 FORMAT(A,' Orbital elements of object ',A,' (',A,')')
  301 FORMAT(A,' Propagated orbits')

      END

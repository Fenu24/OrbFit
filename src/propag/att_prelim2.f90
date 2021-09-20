! ================================================
! ATT_PRELIM
! attributable preliminary orbit
! ================================================
SUBROUTINE att_prelim2(name0,attr,elk,uncatt,fail)
  USE attributable
  USE orbit_elements
  USE output_control
!  USE triangles
  CHARACTER*(*), INTENT(IN) :: name0 ! asteroid designation
  TYPE(attrib), INTENT(IN) :: attr  ! attributable as computed 
  TYPE(orbit_elem), INTENT(OUT) :: elk ! new elements
  TYPE(orb_uncert), INTENT(OUT) :: uncatt ! its covariance
  LOGICAL, INTENT(OUT) :: fail
! END INTERFACE
  INTEGER iunmat,nroots,nrootsd
  DOUBLE PRECISION ::  c(0:5)! polynomial coefficients
  DOUBLE PRECISION, DIMENSION(3) :: roots ! real roots
  DOUBLE PRECISION, DIMENSION(8) :: rootsd ! roots of the derivative
  DOUBLE PRECISION :: a_orb,E_bound ! bound for the energy
! matrix of the r eps theta ref. system, topocentric equatorial
  DOUBLE PRECISION :: refs(3,3)
 ! observer pos/vel, asteroid pos/vel, equatorial heliocentric
  DOUBLE PRECISION :: xo(3), vo(3)
  DOUBLE PRECISION rr, rdot,ecc

  LOGICAL bizarre
  a_orb = 100.d0 ! comet boundary
!   a_orb = 2.5d0 ! main belt
  fail=.true.
  iunmat=0
  CALL admis_reg(name0,iunmat,attr%tdtobs,attr%angles,attr%obscod,  &
   &  nroots,roots,c,refs,xo,vo,a_orb,E_bound,nrootsd,rootsd)
! handle rare cases
  IF(nroots.eq.3)THEN
     write(iun_log,*)'att_prelim: ',name0,' two c.c.',roots(1:nroots)
     rr=(roots(3)+roots(2))/2.d0  ! privilege given to farther component
  ELSEIF(nroots.eq.2)THEN
     WRITE(iun_log,*)'att_prelim: ',name0, ' rare, nroots= ',nroots,roots(1:nroots)
     rr=roots(1)*0.8d0 ! somewhat inside the admissible region
  ELSEIF(nroots.eq.1)THEN
     rr=roots(1)*0.8d0 ! somewhat inside the admissible region
  ENDIF
  rdot=-c(1)/2.d0
  CALL attelements(attr,rr,rdot,elk,uncatt)
  
  fail=bizarre(elk,ecc)
  WRITE(iun_log,100)name0,rr,ecc
100 FORMAT(A,1X,'rho=',f8.3,1X,'ecc=',f8.4)
END SUBROUTINE att_prelim2

! ORBIT_ELEMENTS
! contains:
!
! OUT OF MOD:
!      appmag     correction to absolute magnitude to get apparent


MODULE orbit_elements
USE output_control
USE fund_const
USE dyn_param
IMPLICIT NONE
PRIVATE

! public routines
PUBLIC coo_cha, equal_orbels, write_elems_old, read_elems_old, ecc_peri
PUBLIC convertunc, propagunc
PUBLIC cootyp, attelements, eldiff, tpcar_mtpcar, opik_tpcar
PUBLIC centerbody_mass, att_prelim2 ! ??? used by bizarre
PUBLIC find_elems, wrikep, wromlr,rd_orb !old I/O
PUBLIC write_elems, read_elems           ! new I/O
PUBLIC wro1lh2_ngr, wromlh2, wro1lh2, wro1lr     ! new I/O
PUBLIC wro1lr_matlab
PUBLIC write_mtp_header, write_mtp_record
PUBLIC write_movsample_header, write_movsample_record
PUBLIC write_arsample_header, write_arsample_record
PUBLIC read_movsample_header, read_movsample_record
PUBLIC undefined_orb_uncert
PUBLIC rdelem
PUBLIC v_infty0

DOUBLE PRECISION :: eccbou ! boundary in e between use of Kepler equation
                           ! and use of f-g series

! Generic orbital elements set
TYPE orbit_elem
CHARACTER(LEN=3)       :: coo     ! elements type; known types are:
!    'CAR' : cartesian positions and velocities
!    'EQU' : equinoctal elements
!    'KEP' : classical keplerian elements, singular for
!             0 eccentricity, 0 and 180 degrees inclination;
!    'COM' : cometary elements, valid for ecc.ge.1
!    'ATT' : attributable plus r rdot
!    'OPI' : Opik type elements not yet implemented

!INTEGER :: ndim  ! actual dimension of elements vector

DOUBLE PRECISION, DIMENSION(6) :: coord ! elements vector
              ! DANGER: some of the coordinates may actually be angles...
DOUBLE PRECISION :: t  ! epoch time, MJD, TDT

LOGICAL :: mag_set ! true if some absolute magnitude is available
DOUBLE PRECISION :: h_mag, g_mag ! Absolute magnitude H, opposition effect G

INTEGER :: center ! 0 for Sun, code for planet as in JPL ephemerides

INTEGER :: obscode ! observatory code, for ATT type only

END TYPE orbit_elem

! default value when undefined
DOUBLE PRECISION, DIMENSION(6), PARAMETER :: zero_6d_vect = &
&    (/ 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0 /)

! undefined orbital elements
TYPE(orbit_elem), PARAMETER :: undefined_orbit_elem = ORBIT_ELEM( &
&  'UNK',           & ! no coordynate defined
& zero_6d_vect,     & ! null elements
&   1.d99,          & ! null epoch
&   .false.,        & ! default no magnitudes
&   -9.99d0,        & ! null magnitude
&    0.15d0,        & ! default opposition effect
&    0,             & ! default is heliocentric
&    500            & ! default is not topocentric correction
! &  ,gms            & ! if needed
&     )  !

TYPE orb_uncert
! normal and covariance matrix (warning: covariance not full if not
! all solved)
DOUBLE PRECISION, DIMENSION(ndimx,ndimx):: c,g
! norm of residuals, of last correction, RMS of magnitudes
LOGICAL :: succ ! logical for success (convergent differential corrections)
INTEGER :: ndim  ! actual dimension of solve for vector
END TYPE orb_uncert

! default value when undefined
DOUBLE PRECISION, DIMENSION(6,6), PARAMETER :: zero_6x6_matrix = &
& RESHAPE((/ 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, &
& 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, &
& 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, &
& 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0 /),(/ 6, 6 /))

TYPE orb_with_uncert
TYPE(orbit_elem) :: orbit ! nominal orbit, solution of least squares
TYPE(orb_uncert) :: uncertainty ! covariance matrix etc.
END TYPE orb_with_uncert

! public entities
PUBLIC orbit_elem,undefined_orbit_elem,orb_uncert
PUBLIC eccbou

CONTAINS

SUBROUTINE undefined_orb_uncert(nd,unc)
 INTEGER, INTENT(IN) :: nd ! dimension of matrices
 TYPE(orb_uncert), INTENT(OUT) :: unc
 unc%ndim=nd
 unc%succ=.false.
 unc%g(1:ndimx,1:ndimx)=0.d0
 unc%c(1:ndimx,1:ndimx)=0.d0
END SUBROUTINE undefined_orb_uncert

! ===================================================================
! COO_CHA
! ===================================================================
!   general purpose elements change
!
!  fail _flag = 0 conversion OK
!             = 1 zero eccentricity (conv. to KEP/COM fails)
!             = 2 zero inclination               "
!             = 3 zero ecc and zero inclination   "
!             = 5 eccentricity 1 and more (conv to EQU/KEP fails)
!             = 9 180 deg inclination (conv to KEP/COM/EQU fails)
!             = 10 zero angular momentum (conv to EQU/KEP/COM fails)
!             = 20 inconsistent calling sequence
!  if fail_flag =1,2,3 the jacobian is undefined for kep_equ, replaced
!                by an arbitrary matrix (to avoid undefined on exit)
!  if fail_flag =5 conversion to KEP, EQU is refused, something else
!                  is provided in output (should be COM)
!  if fail_flag > 8 there is no conversion provided....
!               we should set up the messages in a suitable way,
!               providing for a smooth failure
! ===================================================================
SUBROUTINE coo_cha(x,cooy,y,fail_flag,del,obscode)
! MODIFICATION 2009/debris: not allowed to change center, by passing argument iplanet
TYPE(orbit_elem), INTENT(IN) :: x
CHARACTER*3, INTENT(IN) :: cooy
TYPE(orbit_elem), INTENT(INOUT) :: y ! needed to allow for x=y
INTEGER, INTENT(OUT) :: fail_flag ! see above
! partial derivatives \partial y/\partial x
DOUBLE PRECISION, INTENT(OUT), OPTIONAL, DIMENSION(6,6) :: del
INTEGER, INTENT(IN), OPTIONAL :: obscode ! output observer
! END INTERFACE
TYPE(orbit_elem) z
LOGICAL kepgr_in,kepgr_ou,cargr_in,cargr_ou,opposite, done
CHARACTER*3 coox
! target_center is the output gravitational center
INTEGER obscod,target_center
DOUBLE PRECISION, DIMENSION(6,6) :: del1,del2

done=.false.
fail_flag=0
coox=x%coo
IF(rhs.eq.1.or.rhs.eq.3)THEN
   target_center=0
ELSEIF(rhs.eq.2)THEN
   target_center=3
ELSE
   WRITE(*,*)' coo_cha: rhs not initialized, rhs=', rhs
   STOP
ENDIF

IF(coox.eq.cooy)THEN
   y=x
   IF(PRESENT(del)) CALL eye(6,del)
   done=.true.
   RETURN
ENDIF
IF(.not.PRESENT(obscode).and.(cooy.eq.'ATT'))THEN
   obscod=500
ELSEIF(PRESENT(obscode))THEN
   obscod=obscode
ENDIF
kepgr_in=(coox.eq.'KEP'.or.coox.eq.'EQU'.or.coox.eq.'COM'.or.coox.eq.'COT')
kepgr_ou=(cooy.eq.'KEP'.or.cooy.eq.'EQU'.or.cooy.eq.'COM'.or.cooy.eq.'COT')
cargr_in=(coox.eq.'CAR'.or.coox.eq.'ATT')
cargr_ou=(cooy.eq.'CAR'.or.cooy.eq.'ATT')
!opposite=((kepgr_in.and.cargr_ou).or.(cargr_in.and.kepgr_ou))

IF(cargr_in)THEN
   IF(coox.eq.'ATT')THEN
!      IF(planet.eq.0)THEN
! attributable case must be converted to cartesian heliocentric
      IF(PRESENT(del))THEN
         z=car_att(x,target_center,del1) ! convert to CAR
      ELSE
         z=car_att(x,target_center)      ! convert to CAR
      ENDIF
   ELSEIF(coox.eq.'CAR')THEN
      z=x                                 ! already cartesian
      IF(PRESENT(del))THEN
         CALl eye(6,del1)     ! initialize by identity
      ENDIF
   !   ELSEIF(coox.eq.'OPI')THEN
   !      WRITE(*,*)' coo_cha: this conversion not ready ', coox, ' ',cooy
   !      STOP
   ENDIF
   IF(cooy.eq.'CAR')THEN
      y=z ! done
      IF(PRESENT(del))del=del1
   ELSEIF(cooy.eq.'ATT')THEN
      IF(PRESENT(del))THEN
         y=att_car(z,obscod,target_center,del2)
         del=MATMUL(del2,del1)
      ELSE
         y=att_car(z,obscod,target_center)
      ENDIF
   ELSEIF(kepgr_ou)THEN
! convert to elements
      IF(PRESENT(del))THEN
         y=to_elems(z,cooy,fail_flag,del2)
         del=MATMUL(del2,del1)
      ELSE
         y=to_elems(z,cooy,fail_flag)
      ENDIF
! handle fail_flag
      IF(fail_flag.gt.0.and.fail_flag.le.3)THEN
! zero eccentricity/inclination
! EQU are differentiable
         IF(cooy.eq.'EQU')fail_flag=0
! KEP, COM are defined but not diff.
         IF(.not.PRESENT(del))fail_flag=0
      ELSEIF(fail_flag.eq.5)THEN
! force to COT, leaving fail_flag=5 as warning
         IF(PRESENT(del))THEN
 ! should not be necessary
            IF(y%coo.ne.'COT')THEN
               WRITE(ierrou,*)' coo_cha: fail_flag, outcoo ', fail_flag, y%coo
               numerr=numerr+1
               y=to_elems(z,'COT',fail_flag,del2)
            ENDIF
            del=MATMUL(del2,del1)
         ELSE
            IF(y%coo.ne.'COT')THEN ! should not be necessary
               WRITE(ierrou,*)' coo_cha: fail_flag, outcoo ', fail_flag, y%coo
               numerr=numerr+1
               y=to_elems(z,'COT',fail_flag)
            ENDIF
         ENDIF
      ELSE
! fail_flag=9,10,20: nothing to do
      ENDIF
      done=.true.
      RETURN
   ELSE
      PRINT *, ' coo_cha: this should not happen 1, coox=',coox,' cooy=', cooy
      WRITE(*,*) x
      STOP
   ENDIF
ENDIF
IF(kepgr_in.and.cargr_ou)THEN
   IF(PRESENT(del))THEN
      z=car_elems(x,fail_flag,del1)
   ELSE
      z=car_elems(x,fail_flag)
   ENDIF
   IF(cooy.eq.'ATT')THEN
      IF(PRESENT(del))THEN
         y=att_car(z,obscod,target_center,del2)
         del=MATMUL(del2,del1)
      ELSE
         y=att_car(z,obscod,target_center)
      ENDIF
   ELSEIF(cooy.eq.'CAR')THEN
      y=z ! done
      IF(PRESENT(del))del=del1
   ENDIF
! handle fail_flag
   IF(fail_flag.gt.0.and.fail_flag.le.3)THEN
      IF(coox.eq.'EQU')fail_flag=0 ! EQU are differentiable
      IF(.not.PRESENT(del))fail_flag=0 ! KEP, COM are defined but not diff.
   ELSE
! fail_flag=5 cannot happen, 9,10,20 cannot be fixed
   ENDIF
   done=.true.
   RETURN
ENDIF
IF(kepgr_in.and.kepgr_ou)THEN ! twelfe cases
! starting from KEP
   IF(coox.eq.'KEP')THEN
      IF(cooy.eq.'EQU')THEN
         IF(PRESENT(del))THEN
            y=equ_kep(x,del)
         ELSE
            y=equ_kep(x)
         ENDIF
      ELSEIF(cooy.eq.'COM')THEN
         IF(PRESENT(del))THEN
            y=com_kep(x,fail_flag,del)
         ELSE
            y=com_kep(x,fail_flag)
         ENDIF
      ELSEIF(cooy.eq.'COT')THEN
         IF(PRESENT(del))THEN
            y=cot_kep(x,fail_flag,del)
         ELSE
            y=cot_kep(x,fail_flag)
         ENDIF
      ENDIF
      done=.true.
      RETURN
   ENDIF
! starting from EQU
   IF(coox.eq.'EQU')THEN
      IF(cooy.eq.'KEP')THEN
         IF(PRESENT(del))THEN
            y=kep_equ(x,fail_flag,del)
            IF(fail_flag.gt.0) CALl eye(6,del)! del non exi, set to identity
         ELSE
            y=kep_equ(x,fail_flag)
         ENDIF
      ELSEIF(cooy.eq.'COM')THEN
         IF(PRESENT(del))THEN
            z=kep_equ(x,fail_flag,del1)
            IF(fail_flag.eq.0)THEN
               y=com_kep(z,fail_flag,del2)
               del=MATMUL(del2,del1)
            ELSE
               y=com_kep(z,fail_flag)
               del=del1 ! del2 not existent, set to identity
            ENDIF
         ELSE
            z=kep_equ(x,fail_flag) ! always existent, maybe not C1
            y=com_kep(z,fail_flag)
         ENDIF
      ELSEIF(cooy.eq.'COT')THEN
         IF(PRESENT(del))THEN
            z=kep_equ(x,fail_flag,del1)
            IF(fail_flag.eq.0)THEN
               y=cot_kep(z,fail_flag,del2)
               del=MATMUL(del2,del1)
            ELSE
               y=cot_kep(z,fail_flag)
               del=del1 ! del2 not existent, set to identity
            ENDIF
         ELSE
            z=kep_equ(x,fail_flag) ! always existent, maybe not C1
            y=cot_kep(z,fail_flag)
         ENDIF
      ENDIF
      done=.true.
      RETURN
   ENDIF
! starting from COM
   IF(coox.eq.'COM')THEN
      IF(cooy.eq.'KEP')THEN
         IF(PRESENT(del))THEN
            y=kep_com(x,fail_flag,del)
         ELSE
            y=kep_com(x,fail_flag)
         ENDIF
      ELSEIF(cooy.eq.'EQU')THEN
         IF(PRESENT(del))THEN
            z=kep_com(x,fail_flag,del1)
            IF(fail_flag.lt.4)THEN
               y=equ_kep(z,del2)
               del=MATMUL(del2,del1)
            ELSE
               y=z
               del=del1 ! del 2 non existent, set to identity
            ENDIF
         ELSE
            z=kep_com(x,fail_flag)
            IF(fail_flag.lt.4)THEN
               y=equ_kep(z)
            ELSE
               y=z
            ENDIF
         ENDIF
      ELSEIF(cooy.eq.'COT')THEN
         IF(PRESENT(del))THEN
            y=cot_com(x,fail_flag,del)
         ELSE
            y=cot_com(x,fail_flag)
         ENDIF
      ENDIF
      done=.true.
      RETURN
   ENDIF
! starting from COT
   IF(coox.eq.'COT')THEN
      IF(cooy.eq.'KEP')THEN
         IF(PRESENT(del))THEN
            y=kep_cot(x,fail_flag,del)
         ELSE
            y=kep_cot(x,fail_flag)
         ENDIF
      ELSEIF(cooy.eq.'COM')THEN
         IF(PRESENT(del))THEN
            y=com_cot(x,fail_flag,del)
            IF(fail_flag.gt.0) CALl eye(6,del)! del non exi, set to identity
         ELSE
            y=com_cot(x,fail_flag)
         ENDIF
      ELSEIF(cooy.eq.'EQU')THEN
         IF(PRESENT(del))THEN
            z=kep_cot(x,fail_flag,del1)
            IF(fail_flag.lt.4)THEN
               y=equ_kep(z,del2)
               del=MATMUL(del2,del1)
            ELSE
               y=z
               del=del1 ! del 2 non existent, set to identity
            ENDIF
         ELSE
            z=kep_cot(x,fail_flag)
            IF(fail_flag.lt.4)THEN
               y=equ_kep(z)
            ELSE
               y=z
            ENDIF
         ENDIF
      ENDIF
      done=.true.
      RETURN
   ENDIF
   IF(.not.done)THEN
      PRINT *, ' coo_cha: this should not happen 2, coox=',coox,' cooy=', cooy
      WRITE(*,*) x
      STOP
   ENDIF
ENDIF
END SUBROUTINE coo_cha
! ====================================================================
! UTILITIES
! ====================================================================
LOGICAL FUNCTION equal_orbels(el1,el2)
  TYPE(orbit_elem), INTENT(IN) :: el1,el2
  INTEGER i
  equal_orbels=.false.
  IF(el1%coo.ne.el2%coo.or.el1%center.ne.el2%center)RETURN
  DO i=1,6
     IF(el1%coord(i).ne.el2%coord(i))RETURN
  ENDDO
  equal_orbels=.true.
!magnitude does not matter!
END FUNCTION equal_orbels

SUBROUTINE eldiff(el1,el0,del)
  TYPE(orbit_elem), INTENT(IN) :: el1,el0
  DOUBLE PRECISION, DIMENSION(6) :: del
  DOUBLE PRECISION pridif
  IF(el1%coo.ne.el0%coo)THEN
     WRITE(ierrou,*)'eldiff: different coordinates ',el1%coo, el0%coo
     numerr=numerr+1
     del=0.d0
     RETURN
  ENDIF
  IF(el0%coo.eq.'CAR')THEN
     del=el1%coord-el0%coord
  ELSEIF(el0%coo.eq.'EQU')THEN
     del(1:5)=el1%coord(1:5)-el0%coord(1:5)
     del(6)=pridif(el1%coord(6),el0%coord(6))
  ELSEIF(el0%coo.eq.'KEP')THEN
     del(1:3)=el1%coord(1:3)-el0%coord(1:3)
     del(4)=pridif(el1%coord(4),el0%coord(4))
     del(5)=pridif(el1%coord(5),el0%coord(5))
     del(6)=pridif(el1%coord(6),el0%coord(6))
  ELSEIF(el0%coo.eq.'COM')THEN
     del(1:3)=el1%coord(1:3)-el0%coord(1:3)
     del(4)=pridif(el1%coord(4),el0%coord(4))
     del(5)=pridif(el1%coord(5),el0%coord(5))
     del(6)=el1%coord(6)-el0%coord(6)
  ELSEIF(el0%coo.eq.'ATT')THEN
     del(2:6)=el1%coord(2:6)-el0%coord(2:6)
     del(1)=pridif(el1%coord(1),el0%coord(1))
  ELSE
     WRITE(*,*)'eldiff: unknown coordinates ',el1%coo, el0%coo
     STOP
  ENDIF
END SUBROUTINE eldiff

! ====================================================================
! ECC_PERI
! computes eccentricty, pericenter, apocenter, mean motion
! if energy>0 then enne=0, qg=1.d+50
SUBROUTINE ecc_peri(el,ecc,q,qg,enne)
TYPE(orbit_elem), INTENT(IN) :: el
DOUBLE PRECISION, INTENT(OUT) :: ecc,q,qg,enne
! ===================================================================
DOUBLE PRECISION gm,eps,a
DOUBLE precision x(3),y(3),ang(3),vlenz(3),gei,gei2,r0,vel2,alpha,prscal,vsize
TYPE(orbit_elem) el1
INTEGER fail_flag
!
eps=100*epsilon(1.d0)
IF(rhs.eq.1)THEN
   gm=centerbody_mass(0) ! we want always eccentricity heliocentric
ELSEIF(rhs.eq.2)THEN
   gm=centerbody_mass(3) ! we want always eccentricity geocentric
ELSE
! WARNING: automatic choice of method/dynamical model not supported in ORBIT9
   WRITE(*,*)'ecc_peri: ORBIT9 not supported here '
   STOP
ENDIF
IF(el%coo.eq.'EQU')THEN
   ecc=sqrt(el%coord(2)**2+el%coord(3)**2)
   q=el%coord(1)*(1.d0-ecc)
   qg=el%coord(1)*(1.d0+ecc)
   enne=sqrt(gm/el%coord(1)**3)
ELSEIF(el%coo.eq.'KEP')THEN
   ecc=el%coord(2)
   q=el%coord(1)*(1.d0-ecc)
   qg=el%coord(1)*(1.d0+ecc)
   enne=sqrt(gm/el%coord(1)**3)
ELSEIF(el%coo.eq.'COM'.or.el%coo.eq.'COT')THEN
   ecc=el%coord(2)
   q=el%coord(1)
   IF(ecc.gt.1.d0-eps)THEN
      qg=1.d+50
      enne=0.d0
   ELSE
      a=el%coord(1)/(1.d0-ecc)
      qg=a*(1.d0+ecc)
      enne=sqrt(gm/a**3)
   ENDIF
ELSE
   IF(el%coo.eq.'ATT')THEN
! convert to cartesian if it is ATT
      CALL coo_cha(el,'CAR',el1,fail_flag)
   ELSEIF(el%coo.eq.'CAR')THEN
      el1=el
   ELSE
      WRITE(*,*)' ecc_peri: unknown coord type, coo= ',el%coo
      WRITE(*,*) el
      STOP
   ENDIF
! compute eccentricity
   x=el1%coord(1:3)
   y=el1%coord(4:6)
!  radius and velocity squared
   vel2=prscal(y,y)          ! velocity
   r0=vsize(x)               ! distance from center
   call prvec(x,y,ang)       !  angular momentum
! non singular first element
   gei2=prscal(ang,ang)
   gei=sqrt(gei2)
   IF(gei.eq.0.d0)THEN
      WRITE(*,*) ' ecc_peri: zero angular momentum ',el
   ENDIF
   call prvec(y,ang,vlenz) ! Lenz vector
   vlenz=vlenz*(1.d0/gm)-x(1:3)*(1.d0/r0)
   ecc=vsize(vlenz)
   q=gei2/(gm*(1.d0+ecc))
   alpha=vel2-2*gm/r0 ! 2* energy
   IF(ecc.gt.1.d0-eps)THEN
      qg=1.d+50
      enne=0.d0
   ELSE
      qg=gei2/(gm*(1.d0-ecc))
      a=-gm/alpha
      enne=sqrt(gm/a**3)
   ENDIF
ENDIF
END SUBROUTINE ecc_peri

DOUBLE PRECISION FUNCTION centerbody_mass(iplanet)
USE planet_masses
INTEGER,INTENT(IN) :: iplanet
IF(iplanet.eq.0)THEN
! WARNING: does not work in ORBIT9
   IF(rhs.gt.2)THEN
      WRITE(*,*)' centerbody_mass: not ready for rhs=', rhs
      STOP
   ENDIF
   centerbody_mass=gms
ELSEIF(iplanet.eq.3)THEN
   centerbody_mass=gmearth
ELSEIF(iplanet.le.npla+iatrue.and.iplanet.ne.3)THEN
   centerbody_mass=gm(iplanet)
ELSE
   WRITE(*,*)' centerbody_mass: center not available ', iplanet
   STOP
ENDIF
END FUNCTION centerbody_mass
!
TYPE(orbit_elem) FUNCTION coord_translate(el,iplanet)
TYPE(orbit_elem), INTENT(IN) :: el
INTEGER, INTENT(IN) :: iplanet ! planet; JPL ephem codes
DOUBLE PRECISION :: xpla(6)
IF(el%coo.ne.'CAR')THEN
   WRITE(*,*)' coord_translate: wrong input coord ',el
   STOP
ENDIF
coord_translate=el
! check planet, observers are on Earth!
IF(iplanet.eq.3.and.el%center.eq.0)THEN
   CALL earcar(el%t,xpla,1)
   coord_translate%coord=el%coord-xpla
   coord_translate%center=3
ELSEIF(iplanet.eq.0.and.el%center.eq.3)THEN
   CALL earcar(el%t,xpla,1)
   coord_translate%coord=el%coord+xpla
   coord_translate%center=0
ELSE
   WRITE(*,*)' coord_translate: not yet available for ', iplanet, el%center
   STOP
ENDIF
END FUNCTION coord_translate
! ====================================================================
! END UTILITIES
! ====================================================================
! CARTESIAN TO/FROM ELEMENTS
! ====================================================================
! TO_ELEMS
! ====================================================================
TYPE(orbit_elem) FUNCTION to_elems(el,cooy,fail_flag,del)
USE ever_pitkin
TYPE(orbit_elem),INTENT(IN):: el ! must be cartesian
CHARACTER*3, INTENT(IN) :: cooy ! can be 'KEP', 'EQU', 'COM', 'COT'
INTEGER, INTENT(OUT) :: fail_flag ! error flag
! optional jacobian matrix. in output: derivatives d(to_elems)/d(el)
DOUBLE PRECISION, DIMENSION(6,6),INTENT(OUT),OPTIONAL :: del
DOUBLE PRECISION gm
DOUBLE PRECISION, DIMENSION(3) :: x0, y0
DOUBLE PRECISION, DIMENSION(3) :: x, y, ang ! position, velocity, ang. momentum
DOUBLE PRECISION, DIMENSION(3) :: vlenz, fb, gb, ww ! Lenz vector, broucke f and g, pos x vlenz
DOUBLE PRECISION :: gei2, gei, ecc2, ecc, varpi, div, tgi2, omnod
DOUBLE PRECISION :: prscal, vsize, princ, eps
DOUBLE PRECISION :: r0, vel2, alpha, sig0, peri, psi, dt
DOUBLE PRECISION :: enne, rad, beta, chk, ch, ck, fe, x2, y2, cosf, sinf
DOUBLE PRECISION :: eq(6), cosv,sinv,dd
TYPE(orbit_elem) :: elback ! for inversion
DOUBLE PRECISION, DIMENSION(6,6) :: del1
INTEGER :: fail_flag1, ising
DOUBLE PRECISION :: det
! paranoia check of the calling arguments
IF(el%coo.ne.'CAR')THEN
   WRITE(*,*)' to_elems: wrong coordinates in input',x
   fail_flag=20
   RETURN
ENDIF
IF(cooy.ne.'COM'.and.cooy.ne.'KEP'.and.cooy.ne.'EQU'.and.cooy.ne.'COT')THEN
   WRITE(*,*)' to_elems: wrong output coordinates', cooy
   fail_flag=20
   RETURN
ENDIF
! initializations
fail_flag=0
gm=centerbody_mass(el%center) ! mass of central body
eps=100*epsilon(1.d0)         ! small number for rounding off
to_elems=el                   ! copy defaults
! =====================================================================
! portion of the algorithm used for all output coordinates
! =====================================================================
x=el%coord(1:3)               ! cartesian position
y=el%coord(4:6)               ! cartesian velocity
vel2=prscal(y,y)              ! velocity squared
r0=vsize(x)                   ! distance from center
CALL prvec(x,y,ang)           ! angular momentum vector
gei2=prscal(ang,ang)
gei=sqrt(gei2)                ! angular momentum scalar
IF(gei.eq.0.d0)THEN
   WRITE(*,*) ' to_elems: zero angular momentum ',el
   fail_flag=10
   RETURN
ENDIF
CALL prvec(y,ang,vlenz)       ! Lenz vector
vlenz=vlenz*(1.d0/gm)-x(1:3)*(1.d0/r0)
ang=ang*(1.d0)/gei            ! orbit normal unit vector
!   zero divide occurs for inclination of 180 degrees
div=1.d0+ang(3)
IF(div.eq.0.d0)THEN
   WRITE(*,*) ' to_elems: 180 deg. inclination ',el
   fail_flag=9
   RETURN
ENDIF
!  unit vectors of the equinoctal reference system (Broucke and
!  Cefola 1972, CM 5, 303--310) are fb, gb, ang
fb(1)=1.d0-ang(1)**2/div
fb(2)=-ang(1)*ang(2)/div
fb(3)=-ang(1)
CALL prvec(ang,fb,gb)
!  elements related to eccentricity and inclination
eq(2)=prscal(vlenz,gb)
eq(3)=prscal(vlenz,fb)
eq(4)=ang(1)/div
eq(5)=-ang(2)/div
ecc2=eq(2)**2+eq(3)**2
ecc=sqrt(ecc2)                ! eccentricity
IF(ecc.lt.eps)THEN
   varpi=0.d0
   fail_flag=1
ELSE
   varpi=atan2(eq(2),eq(3))   ! longitude of pericenter
ENDIF
tgi2=sqrt(eq(4)**2+eq(5)**2)  ! tanget of inclination/2
if(tgi2.lt.eps)then
   omnod=0.d0
   fail_flag=fail_flag+2
else
   omnod=princ(atan2(eq(4),eq(5))) ! Omega (longitude of node)
endif
alpha=vel2-2*gm/r0            ! 2* energy
IF(alpha.gt.-eps)THEN
  IF(cooy.eq.'EQU'.or.cooy.eq.'KEP') fail_flag=5 ! hyperbolic orbit; is a failure
! if requested coordinates in output are KEP, EQU
ENDIF
sig0=prscal(x,y)              ! scalar product as used in ever_pitkin
! =====================================================================
! part depending upon output coordinates
! =====================================================================
IF(cooy.eq.'COM'.or.cooy.eq.'COT'.or.fail_flag.eq.5)THEN
!         coord(1)=q=a(1-e)
!         coord(2)=e
!         coord(3)=I
!         coord(4)=Omega
!         coord(5)=omega
! if cooy=COM
!         coord(6)=t_0 (time at pericenter)
! elseif cooy=COT
!         coord(6)=v (true anomaly)
   eq(1)=gei2/(gm*(1.d0+ecc))
   eq(2)=ecc
   eq(5)=princ(varpi-omnod)
   eq(3)=2*atan(tgi2)
   eq(4)=omnod
   peri=eq(1)
   IF(cooy.eq.'COM')THEN
      IF(ecc.gt.eps)THEN
! WARNING: use of solve_peri may fail, add error flag which affects
!          the output fail_flag
         CALL solve_peri(r0,sig0,peri,gm,alpha,psi,dt)
         eq(6)=el%t+dt
         fail_flag=0
      ELSE
         fail_flag=1  ! numerically circular orbit
         eq(6)=el%t   ! always at perihelion
      ENDIF
      to_elems%coo='COM'
   ELSEIF(cooy.eq.'COT'.or.fail_flag.eq.5)THEN
      dd=r0*ecc
      cosv=prscal(x,vlenz)/dd  ! cosine of true anomaly
      CALL prvec(vlenz,x,ww)
      sinv=prscal(ang,ww)/dd   ! sine of true anomaly
      eq(6)=atan2(sinv,cosv)
      IF(eq(6).lt.0.d0) eq(6)=eq(6)+dpig
      to_elems%coo='COT'
   ENDIF
ELSE
   eq(1)=-gm/alpha
! selection of different methods
   IF(ecc.gt.eccbou)THEN
      peri=eq(1)*(1-ecc)
      enne=dsqrt(gm/eq(1)**3)
      CALL solve_peri(r0,sig0,peri,gm,alpha,psi,dt)
      eq(6)=princ(-dt*enne+varpi)
   ELSE
! mean longitude from non--singular Kepler equation
      rad=dsqrt(1.d0-ecc2)
      beta=1.d0/(1.d0+rad)
      chk=eq(2)*eq(3)*beta
      ch=1.d0-eq(2)**2*beta
      ck=1.d0-eq(3)**2*beta
      x2=prscal(x,fb)
      y2=prscal(x,gb)
      cosf=eq(3)+(ck*x2-chk*y2)/(eq(1)*rad)
      sinf=eq(2)+(ch*y2-chk*x2)/(eq(1)*rad)
! eccentric longitude
      fe=datan2(sinf,cosf)
      eq(6)=fe+eq(2)*cosf-eq(3)*sinf
! reduction to principal value
      eq(6)=princ(eq(6))
   ENDIF
   IF(cooy.eq.'KEP'.and.fail_flag.eq.0)THEN
      eq(2)=ecc
      eq(4)=omnod
      eq(3)=2*atan(tgi2)
      eq(5)=princ(varpi-omnod)
      eq(6)=princ(eq(6)-varpi)
      to_elems%coo='KEP'
   ELSEIF(cooy.eq.'EQU'.or.(fail_flag.ge.1.and.fail_flag.le.3))THEN
!          coord(2)= h = e sin(varpi) ; varpi = Omega + omega
!          coord(3)= k = e cos(varpi)
!          coord(4)= p = tg(I/2) sin(Omega)
!          coord(5)= q = tg(I/2) cos(Omega)
!          coord(6)= mean longitude, measured from fb vector (mean an + varpi)
      fail_flag=0
      to_elems%coo='EQU'
   ELSE
      PRINT *, ' to_elems: this should not happen ', cooy, fail_flag
   ENDIF
ENDIF
to_elems%coord=eq
IF(PRESENT(DEL))THEN
! compute jacobian matrix by using the inverse transformation
! this is dirty, but the derivatives of elems w.r. to cartesian
! are not as important as the inverse
elback=car_elems(to_elems,fail_flag1,del1)
! WARNING: how to use fail_flag1??????
! inversion
CALL matin(del1,det,6,0,6,ising,1)
IF(ising.eq.0)THEN
   del=del1
ELSE
   del=0.d0
! WARNING: how to set fail_flag??? as above????
ENDIF
ENDIF
END FUNCTION to_elems
! ====================================================================
! CAR_ELEMS
! ====================================================================
TYPE(orbit_elem) FUNCTION car_elems(el,fail_flag,del)
USE ever_pitkin
TYPE(orbit_elem),INTENT(IN):: el ! can be 'KEP', 'EQU', 'COM', 'COT'
INTEGER, INTENT(OUT) :: fail_flag
! optional jacobian matrix. in output: derivatives d(car)/d(el)
DOUBLE PRECISION, DIMENSION(6,6),INTENT(OUT),OPTIONAL :: del
! END INTERFACE
DOUBLE PRECISION :: gm ! center body mass
CHARACTER*(3) :: coox
TYPE(orbit_elem) :: equ ! used for equinoctal, to remove the KEP case
DOUBLE PRECISION, DIMENSION(6,6) :: del1,del0,del2,del3
DOUBLE PRECISION :: eq(6) ! coordinate 6d-vectors
DOUBLE PRECISION, DIMENSION(3) :: x, y, x0, y0, acc
DOUBLE PRECISION :: princ, vsize, eps ! functions, small no
DOUBLE PRECISION :: r0, v0, dt
DOUBLE PRECISION, DIMENSION(3) :: fb,gb,dfdp,dfdq,dgdp,dgdq ! Brouke's vectors
DOUBLE PRECISION :: tgim2, tgim, ecc, upq, ecc2, varpi, enne, mean_an
! vectors used in computing derivatives (not the complete derivatives)
DOUBLE PRECISION :: dedh, dedk, ddigdh, ddigdk, fact, rsvperi,q, upecv2
DOUBLE PRECISION cosn,sinn,coso,sino,cosi,sini, cosv,sinv,rsq,fgfact,eye3(3,3)
! paranoia check of input arguments
coox=el%coo
IF(coox.ne.'COM'.and.coox.ne.'KEP'.and.coox.ne.'EQU'.and.coox.ne.'COT')THEN
   WRITE(*,*)' car_elems: wrong input coordinates', coox
   fail_flag=20
   RETURN
ENDIF
! initializations
fail_flag=0
gm=centerbody_mass(el%center)
eps=100*epsilon(1.d0) ! small number
car_elems=el ! copy other data (e.g., magnitude)
! ================================================================
! the position and velocity at pericenter are computed for all
!  output coordiantes
! ================================================================
IF(coox.eq.'COM'.or.coox.eq.'COT')THEN
   ecc=el%coord(2)
   ecc2=ecc**2
   IF(ecc.lt.eps)THEN
      varpi=0.d0
      fail_flag=1
   ELSE
      varpi=el%coord(4)+el%coord(5)
   ENDIF
   IF(el%coord(3).lt.eps) fail_flag=fail_flag+2
! pericenter
   r0=el%coord(1)
   v0=sqrt(gm*(1.d0+ecc)/el%coord(1))
!   WRITE(*,*) ' r0,v0 ',r0,v0
   cosn=cos(el%coord(4))
   sinn=sin(el%coord(4))
   coso=cos(el%coord(5))
   sino=sin(el%coord(5))
   cosi=cos(el%coord(3))
   sini=sin(el%coord(3))
   x0=(/coso*cosn-sino*sinn*cosi, coso*sinn+sino*cosn*cosi, sino*sini/)*r0
   y0=(/-sino*cosn-coso*sinn*cosi, -sino*sinn+coso*cosn*cosi, coso*sini/)*v0
   IF(coox.eq.'COM')THEN
      dt=el%t-el%coord(6) ! time between peihelion and epoch, for 2-body propagation
   ELSEIF(coox.eq.'COT')THEN
      cosv=cos(el%coord(6)) ! cosine and sine of true anomaly, for 2-body formulae
      sinv=sin(el%coord(6))
   ENDIF
ELSE
   IF(el%coo.eq.'KEP')THEN
      IF(PRESENT(del))THEN
         equ=equ_kep(el,del1)
      ELSE
         equ=equ_kep(el)
      ENDIF
   ELSEIF(el%coo.eq.'EQU')THEN
      equ=el
      IF(PRESENT(del))CALL eye(6,del1)
   ENDIF
   eq=equ%coord
   ecc2=eq(2)**2+eq(3)**2
   ecc=sqrt(ecc2)
! pericenter
   r0=eq(1)*(1.d0-ecc)
   v0=sqrt(gm*(1.d0+ecc)/(eq(1)*(1.d0-ecc)))
   enne=sqrt(gm/eq(1)**3)
!   WRITE(*,*) ' r0,v0, enne ',r0,v0,enne
! varpi= Omega+ omega
   IF(ecc.lt.eps)THEN
      varpi=0.d0
      fail_flag=1
   ELSE
      varpi=atan2(eq(2),eq(3))
   ENDIF
   mean_an=princ(eq(6)-varpi)
   IF(mean_an.gt.pig)mean_an=mean_an-dpig
   dt=mean_an/enne
   tgim2=eq(4)**2+eq(5)**2
   tgim=dsqrt(tgim2)
   IF(tgim.lt.eps) fail_flag=fail_flag+2
   upq=1.d0+tgim2
!  unit vectors of the equinoctal reference system (Broucke and
!  Cefola 1972, CM 5, 303--310) are f, g and the angular momentum
!  unit vector
   fb(1)=(1.d0-eq(4)**2+eq(5)**2)/upq
   fb(2)=2*eq(4)*eq(5)/upq
   fb(3)=-2*eq(4)/upq
   gb(1)=2*eq(4)*eq(5)/upq
   gb(2)=(1.d0+eq(4)**2-eq(5)**2)/upq
   gb(3)=2*eq(5)/upq
! pericenter vector, parallel to Lenz vector
   x0=r0*(fb*cos(varpi)+gb*sin(varpi))
! velocity at pericenter
   y0=v0*(-fb*sin(varpi)+gb*cos(varpi))
ENDIF
! WARNING: for now the use of f-g series is avoided only for coox='COT'
IF(coox.eq.'COT')THEN
   q=el%coord(1)
   rsq=(1.d0+ecc)/(1.d0+ecc*cosv) ! r/q
   rsvperi=sqrt((1.d0+ecc)*q**3/gm)/(1+ecc*cosv)
   x=x0*(rsq*cosv)+y0*(rsvperi*sinv)
   fgfact=sqrt(gm/(q**3*(1.d0+ecc)))
   y=(-sinv*fgfact)*x0+((cosv+ecc)/(1.d0+ecc))*y0
ELSE
   IF(PRESENT(del))THEN
      CALL fser_propag_der(x0,y0,0.d0,dt,gm,x,y,del3)
   ELSE
      CALL fser_propag(x0,y0,0.d0,dt,gm,x,y)
   ENDIF
ENDIF
car_elems%coord(1:3)=x
car_elems%coord(4:6)=y
car_elems%coo='CAR'
IF(PRESENT(del))THEN

! compute partials of x0, y0 w.r. to elements eq(1:5)

   IF(coox.eq.'COM'.or.coox.eq.'COT')THEN
      del0=0.d0
      del0(1:3,1)=  x0*(1/el%coord(1))
      del0(4:6,1)=  y0*(-0.5/el%coord(1))
      del0(1:3,2)=  0.d0
      del0(4:6,2)=  y0*0.5/(1.d0+ecc)
      del0(1:3,3)=  (/ sinn*sino*sini, -cosn*sino*sini, sino*cosi/)*r0
      del0(4:6,3)=  (/ sinn*coso*sini, -cosn*coso*sini, coso*cosi/)*v0
      del0(1:3,4)=  (/-sinn*coso-cosn*sino*cosi, cosn*coso-sinn*sino*cosi, 0.d0/)*r0
      del0(4:6,4)=  (/sinn*sino-cosn*coso*cosi, -cosn*sino-sinn*coso*cosi, 0.d0/)*v0
      fact= r0/v0
      del0(1:3,5)=  fact*y0
      del0(4:6,5)= -(1.d0/fact)*x0
      IF(coox.eq.'COM')THEN
! now use Goodyear to compute jacobian of propagation from t=eq(6) to t=eqn%t
         del=MATMUL(del3,del0)
! then fix derivatives with respect to time of pericenter
         acc=-x*gm/vsize(x)**3
         del(1:3,6)=  -y
         del(4:6,6)=  -acc
      ELSEIF(coox.eq.'COT')THEN
         del3=0.d0 ! del3=\partial (x,y)/\partial (x0,y0)
         CALL eye(3,eye3)
         del3(1:3,1:3)=rsq*cosv*eye3
         del3(1:3,4:6)=rsvperi*sinv*eye3
         del3(4:6,1:3)=-fgfact*sinv*eye3
         del3(4:6,4:6)=((cosv+ecc)/(1.d0+ecc))*eye3
         del2=0.d0 ! del2= \partial (x,y)/\partial el%coord (not through x0,y0)
         upecv2=(1.d0+ecc*cosv)**(-2)
         del2(1:3,1)=1.5d0*sqrt(q*(1.d0+ecc)/gm)*sinv/(1+ecc*cosv) * y0 ! \partial x/\partial q
         del2(1:3,2)=cosv*(1.d0-cosv)*upecv2 * x0    &! \partial x/\partial ecc
&            +0.5d0*sqrt(q**3/((1.d0+ecc)*gm))*sinv*(1.d0-ecc*cosv-2.d0*cosv)*upecv2 * y0
         del2(1:3,6)=-sinv*(1.d0+ecc)*upecv2 * x0  & ! \partial x/\partial v
&             +sqrt(q**3*(1.d0+ecc)/gm)*(cosv+ecc)*upecv2* y0
         del2(4:6,1)=1.5d0*sqrt(gm/(q**5*(1.d0+ecc)))*sinv * x0 ! \partial y/\partial q
         del2(4:6,2)=0.5d0*sqrt(gm/(q**3*(1.d0+ecc)**3))*sinv * x0  & ! \partial y/\partial e
&             + (1.d0-cosv)/(1.d0+ecc)**2 * y0
         del2(4:6,6)= -sqrt(gm/(q**3*(1.d0+ecc)))*cosv * x0 - sinv/(1.d0+ecc) * y0 ! \partial y/\partial v
         del=MATMUL(del3,del0)+del2 ! chain rule
      ENDIF
   ELSEIF(coox.eq.'EQU'.or.coox.eq.'KEP')THEN
      del2=0.d0
      del2(1:3,1)=  x0*(1/eq(1))
      del2(4:6,1)=  y0*(-0.5/eq(1))
      dedh=eq(2)/ecc
      dedk=eq(3)/ecc
      ddigdh=eq(3)/ecc2
      ddigdk=-eq(2)/ecc2
      fact=r0/v0
      del2(1:3,2)= -dedh/(1.d0-ecc)*x0+fact*ddigdh*y0
      del2(1:3,3)= -dedk/(1.d0-ecc)*x0+fact*ddigdk*y0
!      del2(1:3,2)=  eq(1)*(-gb+ddigdh*(-sin(varpi)*fb+cos(varpi)*gb))
      del2(1:3,3)= -dedk/(1.d0-ecc)*x0+fact*ddigdk*y0
!      del2(1:3,3)=  eq(1)*(-fb+ddigdk*(-sin(varpi)*fb+cos(varpi)*gb))
      del2(4:6,2)=  dedh/(1.d0-ecc2)*y0-(1.d0/fact)*ddigdh*x0
      del2(4:6,3)=  dedk/(1.d0-ecc2)*y0-(1.d0/fact)*ddigdk*x0
      dfdp(1)=-2*eq(4)
      dfdp(2)=2*eq(5)
      dfdp(3)=-2.d0
      dfdp=(dfdp-fb*(2*eq(4)))/upq
      dfdq(1)=2*eq(5)
      dfdq(2)=2*eq(4)
      dfdq(3)=0.d0
      dfdq=(dfdq-fb*(2*eq(5)))/upq
      dgdp(1)=2*eq(5)
      dgdp(2)=2*eq(4)
      dgdp(3)=0.d0
      dgdp=(dgdp-gb*(2*eq(4)))/upq
      dgdq(1)=2*eq(4)
      dgdq(2)=-2*eq(5)
      dgdq(3)=2.d0
      dgdq=(dgdq-gb*(2*eq(5)))/upq
      del2(1:3,4)=  r0*(cos(varpi)*dfdp+sin(varpi)*dgdp)
      del2(1:3,5)=  r0*(cos(varpi)*dfdq+sin(varpi)*dgdq)
      del2(4:6,4)=  v0*(-sin(varpi)*dfdp+cos(varpi)*dgdp)
      del2(4:6,5)=  v0*(-sin(varpi)*dfdq+cos(varpi)*dgdq)
      del2(1:6,6)=  0.d0
      del0=MATMUL(del2,del1)
! now use Goodyear to compute jacobian of propagation from t=eq(6) to t=eqn%t
      del=MATMUL(del3,del0)
! then fix derivatives with respect to mean anomaly
      acc=-x*gm/vsize(x)**3
      del(1:3,6)=y/enne
      del(4:6,6)=acc/enne
! then correct derivatives with respect to a to account for enne in dt
      del(1:3,1)= del(1:3,1)+(1.5d0*mean_an/(enne*eq(1)))*y
      del(4:6,1)= del(4:6,1)+(1.5d0*mean_an/(enne*eq(1)))*acc
! then correct derivatives with respedct to h,k to account for varpi in dt
      IF(coox.eq.'EQU')THEN
         del(1:3,2)=del(1:3,2)-(ddigdh/enne)*y
         del(1:3,3)=del(1:3,3)-(ddigdk/enne)*y
         del(4:6,2)=del(4:6,2)-(ddigdh/enne)*acc
         del(4:6,3)=del(4:6,3)-(ddigdk/enne)*acc
      ENDIF
   ENDIF
ENDIF

END FUNCTION car_elems
! ====================================================================
! ELEMENTS INTERNAL CONVERSIONS
! ====================================================================
! KEP_EQU
! keplerian elements from equinoctal
! =======================================================
TYPE(orbit_elem) FUNCTION kep_equ(el,fail_flag,del)
TYPE(orbit_elem),INTENT(IN):: el ! must be equinoctal
INTEGER,INTENT(OUT):: fail_flag
! optional jacobian matrix. in output: derivatives d(el)/d(eqn)
DOUBLE PRECISION, DIMENSION(6,6),INTENT(OUT),OPTIONAL :: del
! end interface
DOUBLE PRECISION dig,ecc,ecc2,tgi2,tgi2s,princ,eps
DOUBLE PRECISION,DIMENSION(6) :: eq
!
IF(el%coo.ne.'EQU')THEN
   WRITE(*,*)' kep_equ: wrong input coordinates',el
   STOP
ENDIF
fail_flag=0
eps=100*epsilon(1.d0) ! negligible number
eq=el%coord
kep_equ=el  !default
kep_equ%coord(1)=eq(1) !semimajor axis
!  test on eccentricity
ecc2=eq(2)**2+eq(3)**2
ecc=sqrt(ecc2)
if(ecc.lt.eps)then
   fail_flag=1 ! but computation is completed
   dig=0.d0
else
   dig=atan2(eq(2),eq(3))
endif
kep_equ%coord(2)=ecc
!   test on tangent of half inclination
tgi2s=eq(4)**2+eq(5)**2
tgi2=sqrt(tgi2s)
if(tgi2.lt.eps)then
   fail_flag=fail_flag+2 ! but computation is completed
   kep_equ%coord(4)=0.d0
else
   kep_equ%coord(4)=atan2(eq(4),eq(5))
endif
kep_equ%coord(3)=2.d0*atan(tgi2)
!   angular variables
kep_equ%coord(5)=dig-kep_equ%coord(4)
kep_equ%coord(6)=eq(6)-dig
kep_equ%coord(4)=princ(kep_equ%coord(4))
kep_equ%coord(5)=princ(kep_equ%coord(5))
kep_equ%coord(6)=princ(kep_equ%coord(6))
kep_equ%coo='KEP'
IF(PRESENT(del).and.fail_flag.eq.0)THEN
   del=0.d0 ! jacobian is not usable in the nearly singular case
   del(1,1)=  1.d0
   del(2,2)=  eq(2)/ecc ! ecc=h^2+k^2
   del(2,3)=  eq(3)/ecc
   del(3,4)=  2*eq(4)/((1+tgi2**2)*tgi2)
   del(3,5)=  2*eq(5)/((1+tgi2**2)*tgi2)
   del(4,4)=  eq(5)/tgi2s ! Omega=atan2(p,q)
   del(4,5)= -eq(4)/tgi2s
   del(5,2)=  eq(3)/ecc2  ! dig=atan2(h,k)
   del(5,3)= -eq(2)/ecc2  ! omega=dig-Omega
   del(5,4)= -eq(5)/tgi2s
   del(5,5)=  eq(4)/tgi2s
   del(6,2)= -eq(3)/ecc2  ! ell=lambda-dig
   del(6,3)=  eq(2)/ecc2
   del(6,6)=  1.d0
! else the jacobian is not defined!!!
ENDIF

RETURN
END FUNCTION kep_equ

! =======================================================
! EQU_KEP
! equinoctal elements from keplerian
! =======================================================
TYPE(orbit_elem) FUNCTION equ_kep(el,del)
TYPE(orbit_elem),INTENT(IN):: el ! must be keplerian
! optional jacobian matrix. in output: derivatives d(el)/d(eqn)
DOUBLE PRECISION, DIMENSION(6,6),INTENT(OUT),OPTIONAL :: del
! end interface
DOUBLE PRECISION dig,ecc,tgim,princ, inc
DOUBLE PRECISION sindig,cosdig,sinom,cosom,dti2
!
IF(el%coo.ne.'KEP')THEN
   WRITE(*,*)' equ_kep: wrong input coordinates',el
   STOP
ENDIF
equ_kep=el  !default
equ_kep%coord(1)=el%coord(1)
dig=el%coord(4)+el%coord(5)
ecc=el%coord(2)
sindig=sin(dig)
cosdig=cos(dig)
equ_kep%coord(2)=ecc*sindig
equ_kep%coord(3)=ecc*cosdig
IF(el%coord(3).ge.pig)THEN
   inc=pig-1.d-4
ELSE
   inc= el%coord(3)
ENDIF
tgim=tan(inc/2.d0)
sinom=sin(el%coord(4))
cosom=cos(el%coord(4))
equ_kep%coord(4)=tgim*sinom
equ_kep%coord(5)=tgim*cosom
equ_kep%coord(6)=dig+el%coord(6)
equ_kep%coord(6)=princ(equ_kep%coord(6))
equ_kep%coo='EQU'
IF(PRESENT(del))THEN
   del=0.d0
   del(1,1)=1.d0
   del(2,2)=sindig
   del(2,5)=equ_kep%coord(3)
   del(2,4)=equ_kep%coord(3)
   del(3,2)=cosdig
   del(3,5)=-equ_kep%coord(2)
   del(3,4)=-equ_kep%coord(2)
   dti2=1.d0/(2*cos(el%coord(3)/2)**2)
   del(4,3)=dti2*sinom
   del(4,4)= equ_kep%coord(5)
   del(5,3)=dti2*cosom
   del(5,4)=- equ_kep%coord(4)
   del(6,4)=1.d0
   del(6,5)=1.d0
   del(6,6)=1.d0
ENDIF
RETURN
END FUNCTION equ_kep
! =======================================================
! COM_KEP
! cometary elelements from keplerian
! =======================================================
TYPE(orbit_elem) FUNCTION com_kep(el,fail_flag,del)
TYPE(orbit_elem),INTENT(IN):: el ! must be keplerian
INTEGER,INTENT(OUT):: fail_flag
! optional jacobian matrix. in output: derivatives d(el)/d(eqn)
DOUBLE PRECISION, DIMENSION(6,6),INTENT(OUT),OPTIONAL :: del
! end interface
DOUBLE PRECISION ecc, eps, gm, enne, princ
DOUBLE PRECISION,DIMENSION(6) :: eq
!
IF(el%coo.ne.'KEP')THEN
   WRITE(*,*)' com_kep: wrong input coordinates',el
   STOP
ENDIF
fail_flag=0
gm=centerbody_mass(el%center)
eps=100*epsilon(1.d0) ! negligible number
eq=el%coord
com_kep=el  !default
ecc=eq(2)
com_kep%coord(1)=eq(1)*(1.d0-ecc)
!  test on eccentricity
if(ecc.lt.eps) fail_flag= 1 ! but computation is completed
!   test on tangent of half inclination
if(eq(3).lt.eps) fail_flag=fail_flag+2 ! but computation is completed
! time of pericenter from mean anomaly
enne=sqrt(gm/eq(1)**3)
eq(6)=princ(eq(6))
com_kep%coord(6)=el%t-eq(6)/enne
com_kep%coo='COM'
IF(PRESENT(del).and.fail_flag.eq.0)THEN
   del=0.d0 ! jacobian is not usable in the nearly singular case
   del(1,1)=  1.d0-ecc
   del(1,2)=  -eq(1)
   del(2,2)=  1.d0
   del(3,3)=  1.d0
   del(4,4)=  1.d0
   del(5,5)=  1.d0
   del(6,1)= -1.5d0*eq(6)/(enne*eq(1))
   del(6,6)=  -1.d0/enne
ENDIF
END FUNCTION com_kep

! =======================================================
! COT_KEP
! cometary elements with true anomaly from keplerian
! =======================================================
TYPE(orbit_elem) FUNCTION cot_kep(el,fail_flag,del)
TYPE(orbit_elem),INTENT(IN):: el ! must be keplerian
INTEGER,INTENT(OUT):: fail_flag
! optional jacobian matrix. in output: derivatives d(el)/d(eqn)
DOUBLE PRECISION, DIMENSION(6,6),INTENT(OUT),OPTIONAL :: del
! end interface
DOUBLE PRECISION ecc, eps, gm, enne, princ
DOUBLE PRECISION, DIMENSION(6) :: eq
DOUBLE PRECISION ell,u,du,cosu,sinu,cosv,sinv,dvdu,dvde,dudell,dude
INTEGER, PARAMETER :: itmax=10
INTEGER i
!
IF(el%coo.ne.'KEP')THEN
   WRITE(*,*)' cot_kep: wrong input coordinates',el
   STOP
ENDIF
fail_flag=0
gm=centerbody_mass(el%center)
eps=100*epsilon(1.d0) ! negligible number
eq=el%coord
cot_kep=el  !default
ecc=eq(2)
cot_kep%coord(1)=eq(1)*(1.d0-ecc)
!  test on eccentricity
if(ecc.lt.eps) fail_flag= 1 ! but computation is completed
!   test on tangent of half inclination
if(eq(3).lt.eps) fail_flag=fail_flag+2 ! but computation is completed
! time of pericenter from mean anomaly
!enne=sqrt(gm/eq(1)**3)
! Kepler's equation to find eccentric anomaly
u=pig
DO i=1,itmax
   cosu=cos(u);sinu=sin(u)
   ell=u-ecc*sinu
   du=(eq(6)-ell)/(1.d0-ecc*cosu)
   u=u+du
   IF(abs(du).lt.eps*10) GOTO 3
ENDDO
WRITE(ierrou,*)' cot_kep: non convergent Kepler equation, u,e,du:', u,ecc,du
3 CONTINUE
cosv=(cosu-ecc)/(1.d0-ecc*cosu)
sinv=sqrt(1.d0-ecc**2)*sinu/(1.d0-ecc*cosu)
cot_kep%coord(6)=atan2(sinv,cosv)
cot_kep%coord(6)=princ(cot_kep%coord(6))
cot_kep%coo='COT'
IF(PRESENT(del).and.fail_flag.eq.0)THEN
   del=0.d0 ! jacobian is not usable in the nearly singular case
   del(1,1)=  1.d0-ecc
   del(1,2)=  -eq(1)
   del(2,2)=  1.d0
   del(3,3)=  1.d0
   del(4,4)=  1.d0
   del(5,5)=  1.d0
   dvdu=sqrt((1.d0+ecc)/(1.d0-ecc))*(cosv+1.d0)/(cosu+1.d0)
   dvde=tan(u/2.d0)*(cosv+1.d0)/sqrt((1.d0+ecc)*(1.d0-ecc)**3)
   dudell=1.d0/(1.d0-ecc*cosu)
   dude=sinu/(1.d0-ecc*cosu)
   del(6,6)=dvdu*dudell
   del(6,2)=dvde+dvdu*dude
ENDIF
END FUNCTION cot_kep
! =======================================================
! KEP_COT
! keplerian elelements from cometary with true anomaly
! =======================================================
TYPE(orbit_elem) FUNCTION kep_cot(el,fail_flag,del)
TYPE(orbit_elem),INTENT(IN):: el ! must be cometary
INTEGER,INTENT(OUT):: fail_flag
! optional jacobian matrix. in output: derivatives d(el)/d(eqn)
DOUBLE PRECISION, DIMENSION(6,6),INTENT(OUT),OPTIONAL :: del
! end interface
DOUBLE PRECISION ecc, eps, gm, enne, a, dldn
DOUBLE PRECISION,DIMENSION(6) :: eq
DOUBLE PRECISION r,x1,y1,sinu,cosu,u,delldu,dellde,dudv,dude,v2
!
IF(el%coo.ne.'COT')THEN
   WRITE(*,*)' kep_cot: wrong input coordinates',el
   STOP
ENDIF
fail_flag=0
gm=centerbody_mass(el%center)
eps=100*epsilon(1.d0) ! negligible number
eq=el%coord
kep_cot=el  !default
ecc=eq(2)
if(ecc.ge.1.d0-eps)THEN
! this conversion is impossible
   fail_flag=5
! in output the elements are left as 'COT', del= identity
   IF(PRESENT(del)) CALL eye(6,del)
   RETURN
ENDIF
a=eq(1)/(1.d0-ecc)
kep_cot%coord(1)=a
!  test on eccentricity
if(ecc.lt.eps) fail_flag= 1 ! but computation is completed
!   test on tangent of half inclination
if(eq(3).lt.eps) fail_flag=fail_flag+2 ! but computation is completed
! eccentric anomaly from true anomaly
r=a*(1.d0-ecc**2)/(1.d0+ecc*cos(eq(6)))
x1=r*cos(eq(6)); y1=r*sin(eq(6))
sinu=y1/(a*sqrt(1.d0-ecc**2)); cosu=x1/a+ecc
u=atan2(sinu,cosu)
! kepler's equation (easy way)
kep_cot%coord(6)=u-ecc*sinu
kep_cot%coo='KEP' ! successful conversion
IF(PRESENT(del).and.fail_flag.eq.0)THEN
   del=0.d0 ! jacobian is not usable in the nearly singular case
   del(1,1)=  1.d0/(1.d0-ecc)
   del(1,2)=  +eq(1)/(1.d0-ecc)**2
   del(2,2)=  1.d0
   del(3,3)=  1.d0
   del(4,4)=  1.d0
   del(5,5)=  1.d0
! derivatives of mean anomaly (composite through eccentric anomaly)
   delldu=1.d0-ecc*cosu
   dellde=-sinu
   v2=eq(6)/2.d0
   dudv=sqrt((1.d0-ecc)/(1.d0+ecc))*cos(u/2.d0)**2/cos(v2)**2
   dude=(-2*tan(v2)*cos(u/2.d0)**2)/sqrt((1.d0-ecc)*(1.d0+ecc)**3)
   del(6,6)=delldu*dudv
   del(6,2)=dellde + delldu*dude
ENDIF
END FUNCTION kep_cot
! =======================================================
! COT_COM
! cometary elements with true anomaly from cometary
! =======================================================
TYPE(orbit_elem) FUNCTION cot_com(el,fail_flag,del)
TYPE(orbit_elem),INTENT(IN):: el ! must be keplerian
INTEGER,INTENT(OUT):: fail_flag
! optional jacobian matrix. in output: derivatives d(el)/d(eqn)
DOUBLE PRECISION, DIMENSION(6,6),INTENT(OUT),OPTIONAL :: del
! end interface
DOUBLE PRECISION ecc, eps , gm
DOUBLE PRECISION delcar(6,6),dcarel(6,6)
TYPE(orbit_elem) elcar,elcot
DOUBLE PRECISION dt, alpha,v,dvdq,dvde,dvdt0, vel2
INTEGER fail_flag1, fail_flag2
! ======================================================
IF(el%coo.ne.'COM')THEN
   WRITE(*,*)' cot_com: wrong input coordinates',el
   STOP
ENDIF
fail_flag=0
gm=centerbody_mass(el%center)
eps=100*epsilon(1.d0) ! negligible number
! eq=el%coord
cot_com=el  !default
ecc=el%coord(2)
!  test on eccentricity
if(ecc.lt.eps) fail_flag= 1 ! but computation is completed
!   test on tangent of half inclination
if(el%coord(3).lt.eps) fail_flag=fail_flag+2 ! but computation is completed
! find true anomaly from dt
dt=el%t-el%coord(6)
alpha=gm*(el%coord(2)-1.d0)/el%coord(1) ! 2* energy
! CALL ta_from_t0(el%coord(1),ecc,dt,gm,alpha,v,dvdq,dvde,dvdt0)
! cot_com%coord(6)=v
cot_com%coo='COT'
! use brutal composition of two available changes
IF(PRESENT(del))THEN
   CALL coo_cha(el,'CAR',elcar,fail_flag1,delcar)
   CALL coo_cha(elcar,'COT',elcot, fail_flag2,dcarel)
   del=MATMUL(dcarel,delcar)
ELSE
   CALL coo_cha(el,'CAR',elcar,fail_flag1)
   CALL coo_cha(elcar,'COT',elcot, fail_flag2)
ENDIF
cot_com%coord(6)=elcot%coord(6)
! cot_com%coord(6)=v
IF(PRESENT(del))THEN
  CALL eye(5,del(1:5,1:5))
!  del(6,1)=dvdq ! NOT READY
!  del(6,2)=dvde
  del(6,3:5)=0.d0
!  del(6,6)=dvdt0
  del(1:5,6)=0.d0
ENDIF
RETURN
! WARNING: what to do with fail_flag1 and fail_flag2 ???
END FUNCTION cot_com
! ======================================================
! TA_FROM_T0 true anomaly from time of perihelion
! ======================================================
SUBROUTINE ta_from_t0(q,ecc,dt,mu,alpha,v,dvdq,dvde,dvdt0)
  USE ever_pitkin
  DOUBLE PRECISION,INTENT(IN) :: q,ecc ! perihelion distance, eccentricity
  DOUBLE PRECISION,INTENT(IN) :: dt ! dt=t-t0 difference between
                       ! present time and time of passage at perihelion
  DOUBLE PRECISION,INTENT(IN) :: mu, alpha ! grav. constant, 2*energy
  DOUBLE PRECISION,INTENT(OUT) :: v ! true anomaly
  DOUBLE PRECISION,INTENT(OUT) :: dvdq,dvde,dvdt0 ! derivatives of true anomaly
! end interface                                   ! w.r.t. q,e,t0
  DOUBLE PRECISION :: sig0, rdot, psi, s0,s1,s2,s3, r
  DOUBLE PRECISION :: cosv,dvdr, dpsidt0, drdt0,drdpsi, ds2da,ds0da, drda
  DOUBLE PRECISION s0v,s1v,s2v,s3v,dd, psiv, eccv,qv, vv,cosvv
! *******************************
! compute psi at time t (psi0=0)
! *******************************
  sig0 = 0.d0
! solve universal Kepler's equation and compute Stumpff's functions
  CALL solve_kepuniv2(dt,q,sig0,mu,alpha,ecc,psi,s0,s1,s2,s3,DS2DA=ds2da,DS0DA=ds0da)
  dd=1.d-12
!  CALL s_funct(psi,alpha+dd,s0v,s1v,s2v,s3v)
!  CALL solve_kepuniv2(dt,q,sig0,mu,alpha+dd,ecc,psiv,s0v,s1v,s2v,s3v)
!  WRITE(*,*)' ds2da=',(s2v-s2)/dd, ds2da, s2
!  WRITE(*,*)' ds0da=',(s0v-s0)/dd, ds0da, s0
!  WRITE(*,*)' psi=', psi, psiv,psi-psiv
!  ds2da=(s2v-s2)/dd
!  ds0da=(s0v-s0)/dd
! *******************************
! compute r at time t
! *******************************
  r = q*s0 + sig0*s1 + mu*s2
  drdpsi = sig0*s0 + (mu+alpha*q)*s1
  cosv =(-1.d0 + q/r*(1.d0+ecc))/ecc
  v=acos(cosv)*sign(1.d0,drdpsi)

  qv=q+dd
  cosvv =(-1.d0 + qv/r*(1.d0+ecc))/ecc
  vv=acos(cosvv)*sign(1.d0,drdpsi)
  dvdq=(vv-v)/dd
  WRITE(*,*)' dvdq direct inc=', dvdq
  eccv=ecc+dd
  cosvv =(-1.d0 + q/r*(1.d0+eccv))/eccv
  vv=acos(cosvv)*sign(1.d0,drdpsi)
  dvde=(vv-v)/dd
  WRITE(*,*)' dvde direct inc=', dvde

! derivatives of v w.r. to q,e
! ***********WARNING: WRONG***************
! missing dependence of r from q,e
! dr/dq=s0+q*d(s0)/dq+ mu*d(s2)/dq
! dr/d(alpha)=q*d(s0)/d(alpha)+ mu*d(s2)/d(alpha)
! d(sj)/dq=d(sj)/d(alpha) d(alpha)/dq ; d(sj)/de=d(sj)/d(alpha) d(alpha)/de
! d(alpha)/dq=-alpha/q ; d(alpha)/de= mu/q
  dvdq = -1.d0/sqrt(1.d0-cosv**2)*(1.d0+ecc)/(r*ecc) *sign(1.d0,drdpsi)
  WRITE(*,*)' dvdq direct=', dvdq
  dvdr =  1.d0/sqrt(1.d0-cosv**2)*q*(1.d0+ecc)/(r**2*ecc)*sign(1.d0,drdpsi)! OK
  drda= q*ds0da + mu*ds2da
!  drdq=s0+(q*ds0da+mu*ds2da)*(-alpha/q)
  dvdq=dvdq + dvdr*(drda*(-alpha/q)+s0)
  dvde = -1.d0/sqrt(1.d0-cosv**2)*(1.d0-q/r)/ecc**2*sign(1.d0,drdpsi)
  WRITE(*,*)' dvde direct=', dvde
  dvde= dvde + dvdr*drda*mu/q
! ************END WARNING*****************
! derivative of v w.r. to t0
  dpsidt0 = -1.d0/(q*s0 + sig0*s1 + mu*s2)
  drdt0 = drdpsi*dpsidt0
  dvdt0 = dvdr*drdt0
END SUBROUTINE ta_from_t0

! =======================================================
! COM_COT
! cometary elelements from cometary with true anomaly
! =======================================================
TYPE(orbit_elem) FUNCTION com_cot(el,fail_flag,del)
TYPE(orbit_elem),INTENT(IN):: el ! must be cometary
INTEGER,INTENT(OUT):: fail_flag
! optional jacobian matrix. in output: derivatives d(el)/d(eqn)
DOUBLE PRECISION, DIMENSION(6,6),INTENT(OUT),OPTIONAL :: del
! end interface
DOUBLE PRECISION ecc, eps, gm !, enne, a, dldn, dt
!DOUBLE PRECISION,DIMENSION(6) :: eq
!DOUBLE PRECISION r,x1,y1,sinu,cosu,u,delldu,dellde,dudv,dude,v2,ell,delldv
DOUBLE PRECISION delcar(6,6),dcarel(6,6)
TYPE(orbit_elem) elcar,elcom
INTEGER fail_flag1, fail_flag2
! =======================================================
IF(el%coo.ne.'COT')THEN
   WRITE(*,*)' com_cot: wrong input coordinates',el
   STOP
ENDIF
fail_flag=0
eps=100*epsilon(1.d0) ! negligible number
com_cot=el  !default
ecc=el%coord(2) ! eccentricity
!  test on eccentricity
IF(ecc.lt.eps) fail_flag= 1 ! but computation is completed
!   test on inclination
IF(el%coord(3).lt.eps) fail_flag=fail_flag+2 ! but computation is completed
! use brutal composition of two available changes
IF(PRESENT(del))THEN
   CALL coo_cha(el,'CAR',elcar,fail_flag1,delcar)
   CALL coo_cha(elcar,'COM',elcom, fail_flag2,dcarel)
   del=MATMUL(dcarel,delcar)
ELSE
   CALL coo_cha(el,'CAR',elcar,fail_flag1)
   CALL coo_cha(elcar,'COM',elcom, fail_flag2)
ENDIF
! WARNING: what to do with fail_flag1 and fail_flag2 ???
com_cot=elcom
END FUNCTION com_cot
! =======================================================
! KEP_COM
! keplerian elelements from cometary
! =======================================================
TYPE(orbit_elem) FUNCTION kep_com(el,fail_flag,del)
TYPE(orbit_elem),INTENT(IN):: el ! must be cometary
INTEGER,INTENT(OUT):: fail_flag
! optional jacobian matrix. in output: derivatives d(el)/d(eqn)
DOUBLE PRECISION, DIMENSION(6,6),INTENT(OUT),OPTIONAL :: del
! end interface
DOUBLE PRECISION ecc, eps, gm, enne, a, dldn
DOUBLE PRECISION,DIMENSION(6) :: eq
!
IF(el%coo.ne.'COM')THEN
   WRITE(*,*)' kep_com: wrong input coordinates',el
   STOP
ENDIF
fail_flag=0
gm=centerbody_mass(el%center)
eps=100*epsilon(1.d0) ! negligible number
eq=el%coord
kep_com=el  !default
ecc=eq(2)
if(ecc.ge.1.d0-eps)THEN
! this conversion is impossible
   fail_flag=5
! in output the elements are left as 'COM', del= identity
   IF(PRESENT(del)) CALL eye(6,del)
   RETURN
ENDIF
a=eq(1)/(1.d0-ecc)
kep_com%coord(1)=a
!  test on eccentricity
if(ecc.lt.eps) fail_flag= 1 ! but computation is completed
!   test on tangent of half inclination
if(eq(3).lt.eps) fail_flag=fail_flag+2 ! but computation is completed
! time of pericenter from mean anomaly
enne=sqrt(gm/a**3)
kep_com%coord(6)=(el%t-eq(6))*enne
kep_com%coo='KEP'
IF(PRESENT(del).and.fail_flag.eq.0)THEN
   del=0.d0 ! jacobian is not usable in the nearly singular case
   del(1,1)=  1.d0/(1.d0-ecc)
   del(1,2)=  +eq(1)/(1.d0-ecc)**2
   del(2,2)=  1.d0
   del(3,3)=  1.d0
   del(4,4)=  1.d0
   del(5,5)=  1.d0
   dldn=      el%t-eq(6)
   del(6,1)=  -1.5d0*enne/eq(1)*dldn
   del(6,2)=  -1.5d0*enne/(1.d0-ecc)*dldn
   del(6,6)=  -enne
ENDIF

RETURN
END FUNCTION kep_com

! ====================================================================
! END ELEMENTS INTERNAL CONVERSIONS
! ======================================================================
! ATT submodule
! ====================================================================
! =======================================================
! ATT_CAR
! attributable elements from cartesian position and velocity,
! planetocentric
! attributable is:
!     %coord(1) = alpha
!     %coord(2) = delta
!     %coord(3) = d(alpha)/dt
!     %coord(4) = d(delta)/dt
!     %coord(5) = range
!     %coord(6) = range rate
! WARNING: att in equatorial, car in ecliptic coordinates
! if the center is different from the Earth
! =======================================================
TYPE(orbit_elem) FUNCTION att_car(el,obscode,target_center,del)
  TYPE(orbit_elem),INTENT(IN):: el
  INTEGER, INTENT(IN) :: obscode
  INTEGER, INTENT(IN) ::target_center
! optional jacobian matrix. in output: derivatives d(el)/d(eqn)
  DOUBLE PRECISION, DIMENSION(6,6),INTENT(OUT),OPTIONAL :: del
! end interface
  DOUBLE PRECISION, DIMENSION(6) :: att ! attributable
  DOUBLE PRECISION, DIMENSION(3) :: dadx,dddx,d,dv !first derivatives,
  DOUBLE PRECISION, DIMENSION(3,3) :: ddadx,ddddx ! second derivatives
  DOUBLE PRECISION :: dz,dis0,vsize,prscal,tobs,princ ! distances, functions
  DOUBLE PRECISION, DIMENSION(6) :: xtop,xobs
! aux. var for computation of second derivatives
  DOUBLE PRECISION den,x2,y2,z2,x4,y4
  DOUBLE PRECISION x2y2,x2z2,y2z2, ddd(3)
! =======================================================
  IF(el%coo.ne.'CAR')THEN
     WRITE(*,*)' att_car: wrong source coordinates ', el%coo
     STOP
  ENDIF
  IF(target_center.ne.0.and.target_center.ne.3)THEN
     WRITE(*,*)'att_car: Center different from Earth or Sun'
     STOP
  END IF
! compute observer position, taking into account light time
  CALL iter_obs(el%t,el%coord,obscode,xtop,tobs,dis0,target_center)
  d=xtop(1:3)
  dv=xtop(4:6)
! Computation of observation: right ascension (radians)
  dz=d(1)**2+d(2)**2
  if (dz.le.100*epsilon(1.d0)) then
! remove singularity at poles
     att(1)=0.d0
  else
     att(1)=atan2(d(2),d(1))
  endif
! Computation of observation: declination (radians)
  att(2)=asin(d(3)/dis0)
! Computation of first derivatives of $\alpha$ and $\delta$ w.r. to posi
! (if required): we derived eq. (2.20)
  dadx(1)=-d(2)/dz
  dadx(2)=d(1)/dz
  dadx(3)=0.d0
  dddx(1)=-d(3)*(d(1)/(sqrt(dz)*dis0**2))
  dddx(2)=-d(3)*(d(2)/(sqrt(dz)*dis0**2))
  dddx(3)=sqrt(dz)/dis0**2
! Apparent motion:
  att(3)=prscal(dadx,dv)
  att(4)=prscal(dddx,dv)
! range, range rate
  att(5)=dis0
  att(6)=prscal(d,dv)/dis0
! assign attributable type elements
  att_car=el !defaults
  att_car%center=3
  att_car%coord=att
  att_car%coo='ATT'
! aberration correction NO!!!
  att_car%obscode=obscode
! h magnitude is already the absolute one
! derivatives, if required
  IF(PRESENT(del))THEN
     IF(target_center.eq.0)THEN
        del=0.d0
        del(1,1:3)= MATMUL(roteqec,dadx) ! = MATMUL(dadx,roteceq)
        del(2,1:3)=MATMUL(roteqec,dddx)
        del(3,4:6)=del(1,1:3) ! commutation of d/dt
        del(4,4:6)=del(2,1:3) ! commutation of d/dt
        del(5,1:3)=MATMUL(roteqec,d*(1.d0/dis0))
        del(6,1:3)=MATMUL(roteqec,d*(-att(6)/dis0**2)+dv*(1.d0/dis0))
        del(6,4:6)=del(5,1:3) ! commutation of d/dt
! =====================================================================
! partials of adot,ddot with respect to positions require the
! second derivatives of alpha, delta with respect to the
! equatorial reference system
! =====================================================================
! Computation of second derivatives of $\alpha$ w.r. to positions
        ddadx(1,1)=2.d0*d(1)*d(2)/dz**2
        ddadx(1,2)=(d(2)**2-d(1)**2)/dz**2
        ddadx(2,1)=ddadx(1,2)
        ddadx(2,2)=-ddadx(1,1)
        ddadx(3,1)=0.d0
        ddadx(3,2)=0.d0
        ddadx(3,3)=0.d0
        ddadx(2,3)=0.d0
        ddadx(1,3)=0.d0
! chain rule for derivatives of adot
        ddd=MATMUL(ddadx,dv)
! rotation to the equatorial reference system
        del(3,1:3)=MATMUL(roteqec,ddd)
! =======================================================
! Computation of second derivatives of $\delta$ w.r. to positions
        den=1.d0/(dis0**4*dz*sqrt(dz))
        x2=d(1)**2
        y2=d(2)**2
        z2=d(3)**2
        x4=x2*x2
        y4=y2*y2
        x2y2=x2*y2
        x2z2=x2*z2
        y2z2=y2*z2
        ddddx(1,1)=d(3)*(2.d0*x4+x2y2-y2z2-y4)*den
        ddddx(2,2)=d(3)*(2.d0*y4+x2y2-x2z2-x4)*den
        ddddx(1,2)=d(1)*d(2)*d(3)*(z2+3.d0*x2+3.d0*y2)*den
        ddddx(2,1)=ddddx(1,2)
        ddddx(3,3)=-2.d0*d(3)*dz**2*den
        ddddx(1,3)=d(1)*dz*(z2-x2-y2)*den
        ddddx(3,1)=ddddx(1,3)
        ddddx(2,3)=d(2)*dz*(z2-x2-y2)*den
        ddddx(3,2)=ddddx(2,3)
! chain rule for derivatives of ddot
        CALL prodmv(ddd,ddddx,dv)
! rotation to the equatorial reference system
        CALL prodmv(del(4,1:3),roteqec,ddd)
! the same for the earth, but without rotation
     ELSEIF(target_center.eq.3)THEN
        del=0.d0
        del(1,1:3)=dadx
        del(2,1:3)=dddx
        del(3,4:6)=del(1,1:3)
        del(4,4:6)=del(2,1:3)
        del(5,1:3)=d*(1.d0/dis0)
        del(6,1:3)=d*(-att(6)/dis0**2)+dv*(1.d0/dis0)
        del(6,4:6)=del(5,1:3) ! commutation of d/dt
        ddadx(1,1)=2.d0*d(1)*d(2)/dz**2
        ddadx(1,2)=(d(2)**2-d(1)**2)/dz**2
        ddadx(2,1)=ddadx(1,2)
        ddadx(2,2)=-ddadx(1,1)
        ddadx(3,1)=0.d0
        ddadx(3,2)=0.d0
        ddadx(3,3)=0.d0
        ddadx(2,3)=0.d0
        ddadx(1,3)=0.d0
        CALL prodmv(ddd,ddadx,dv)
        del(3,1:3)=ddd
        den=1.d0/(dis0**4*dz*sqrt(dz))
        x2=d(1)**2
        y2=d(2)**2
        z2=d(3)**2
        x4=x2*x2
        y4=y2*y2
        x2y2=x2*y2
        x2z2=x2*z2
        y2z2=y2*z2
        ddddx(1,1)=d(3)*(2.d0*x4+x2y2-y2z2-y4)*den
        ddddx(2,2)=d(3)*(2.d0*y4+x2y2-x2z2-x4)*den
        ddddx(1,2)=d(1)*d(2)*d(3)*(z2+3.d0*x2+3.d0*y2)*den
        ddddx(2,1)=ddddx(1,2)
        ddddx(3,3)=-2.d0*d(3)*dz**2*den
        ddddx(1,3)=d(1)*dz*(z2-x2-y2)*den
        ddddx(3,1)=ddddx(1,3)
        ddddx(2,3)=d(2)*dz*(z2-x2-y2)*den
        ddddx(3,2)=ddddx(2,3)
        CALL prodmv(ddd,ddddx,dv)
        del(4,1:3)=ddd
     ELSE
        WRITE(*,*)'att_car: Center different from Earth or Sun'
        STOP
     END IF
  ENDIF
END FUNCTION att_car

! =======================================================
! ITER-OBS
! =======================================================
SUBROUTINE iter_obs(t,x,obscode,xtop,tobs,r,target_center)
USE reference_systems, ONLY: observer_position
! x is heliocentric ecliptic (tc=0), geocentric equatorial (tc=3)
DOUBLE PRECISION, INTENT(IN) :: t, x(6)
INTEGER, INTENT(IN) :: obscode, target_center
DOUBLE PRECISION, DIMENSION(6), INTENT(OUT):: xtop ! topocentric
! pos/vel with respect to observer obscode at appropriate time
! also ecliptic heliocentric position of observer
DOUBLE PRECISION, INTENT(OUT) :: tobs, r !obs time, range
DOUBLE PRECISION, DIMENSION(6) :: xpla
DOUBLE PRECISION vsize,tc, xobs(6)
INTEGER, PARAMETER :: jmax=5
DOUBLE PRECISION, PARAMETER :: tcont=1.d-8
INTEGER j
! change center: get ecliptic coordinates of the observer
tc=t
DO j=1,jmax
   IF(target_center.eq.0)THEN
      CALL observer_position(tc,xobs(1:3),xobs(4:6),OBSCODE=obscode)
      CALL earcar(tc,xpla,1) ! xpla=heliocentric ecliptic of geocenter
      xobs=xpla+xobs ! now xobs= heliocentric ecliptic of topos
      xtop=x-xobs  ! xtop=topocentric ecliptic of asteroid
      CALL prodmv(xtop(1:3),roteceq,xtop(1:3)) ! xtop=topocentric equatorial
      CALL prodmv(xtop(4:6),roteceq,xtop(4:6))
      IF(rhs.gt.1) STOP '**** noniter_obs: target_center=0 and rhs>1 not allowed***'
   ELSEIF(target_center.eq.3)THEN
      CALL observer_position(tc,xobs(1:3),xobs(4:6),OBSCODE=obscode)
! rotation to equatorial is handled by observer_position if rhs=2
      IF(rhs.ne.2) STOP '**** noniter_obs: target_center=3 and rhs=1 not allowed***'
      xtop=x-xobs ! topocentric equatorial
   ELSE
      STOP
   ENDIF
! ecliptic topocentric vector of the asteroid
   r=vsize(xtop)
   tobs=t +r/vlight !aberration correction
   IF(abs(tc-tobs).lt.tcont)RETURN
   tc=tobs
ENDDO
WRITE(*,*) 'iter_obs: slow convergence ', j-1, tobs, t, r
END SUBROUTINE iter_obs

! =======================================================
! NONITER_OBS
! =======================================================
SUBROUTINE noniter_obs(tobs,obscode,xobs, target_center)
  USE reference_systems, ONLY: observer_position
  DOUBLE PRECISION, INTENT(IN) :: tobs
  INTEGER, INTENT(IN) :: obscode, target_center
  DOUBLE PRECISION, DIMENSION(6), INTENT(OUT):: xobs ! topocentric ecliptic
! pos/vel with respect to observer obscode at appropriate time
  DOUBLE PRECISION, DIMENSION(6) :: xpla
! change center: get ecliptic coordinates of the observer
  IF(target_center.eq.0)THEN
     CALL observer_position(tobs,xobs(1:3),xobs(4:6),OBSCODE=obscode)
     CALL earcar(tobs,xpla,1) ! xpla=heliocentric ecliptic of geocenter
     xobs=xpla+xobs ! now xobs= heliocentric ecliptic of topos
     IF(rhs.gt.1) STOP '**** noniter_obs: target_center=0 and rhs>1 not allowed***'
  ELSEIF(target_center.eq.3)THEN
     CALL observer_position(tobs,xobs(1:3),xobs(4:6),OBSCODE=obscode)
! rotation to equatorial is handled by observer_position if rhs=2
     IF(rhs.ne.2) STOP '**** noniter_obs: target_center=3 and rhs=1 not allowed***'
  ELSE
     STOP
  ENDIF
END SUBROUTINE noniter_obs

! =======================================================
! CAR_ATT
! cartesian position and velocity, planetocentric, from
! attributable elements
! WARNING:  car in ecliptic coordinates, att in equatorial
! =======================================================
TYPE(orbit_elem) FUNCTION car_att(el,target_center,del)
  TYPE(orbit_elem),INTENT(IN) :: el
! New variable to know the center
  INTEGER, INTENT(in) ::target_center
! optional jacobian matrix. in output: derivatives d(el)/d(eqn)
  DOUBLE PRECISION, DIMENSION(6,6),INTENT(OUT),OPTIONAL :: del
! end interface
  DOUBLE PRECISION, DIMENSION(3) :: rhat,ralphat,rdelhat,ralpalp, &
     & rdeldel,ralpdel
  DOUBLE PRECISION, DIMENSION(6) :: att,x
  DOUBLE PRECISION cosa,cosd,sina,sind,r,rdot,tobs
  DOUBLE PRECISION, DIMENSION(6) :: xobs
!
  IF(el%coo.ne.'ATT')THEN
     WRITE(*,*)' car_att: wrong source coordinates ', el%coo
     STOP
  ENDIF
  IF(target_center.ne.3.and.target_center.ne.0)THEN
     WRITE(*,*)'car_att: The center is different from Earth and Sun'
     STOP
  ENDIF
  att=el%coord
! UNIT vectors of the r alpha delta ref. system
  cosa=cos(att(1))
  cosd=cos(att(2))
  sina=sin(att(1))
  sind=sin(att(2))
  r=att(5)
  rhat(1)=cosa*cosd
  rhat(2)=sina*cosd
  rhat(3)=sind
  ralphat(1)=-sina*cosd ! note: contains cosd and is not an unit vector
  ralphat(2)= cosa*cosd ! that is d(rhat)/d(alpha)=ralphat
  ralphat(3)=0.d0
  rdelhat(1)=-cosa*sind
  rdelhat(2)=-sina*sind
  rdelhat(3)=cosd
! second derivatives
  ralpalp(1)=-cosa*cosd
  ralpalp(2)=-sina*cosd
  ralpalp(3)=0.d0
  ralpdel(1)=sina*sind
  ralpdel(2)=-cosa*sind
  ralpdel(3)=0.d0
  rdeldel(1)=-cosa*cosd
  rdeldel(2)=-sina*cosd
  rdeldel(3)=-sind
  IF(target_center.eq.0)THEN
! planetocentric position and velocity, equatorial, converted to ecliptic
     CALL prodmv(x(1:3),roteqec,att(5)*rhat)
     CALL prodmv(x(4:6),roteqec,att(6)*rhat+r*(att(3)*ralphat+att(4)*rdelhat))
  ELSEIF(target_center.eq.3)THEN
     x(1:3)=att(5)*rhat
     x(4:6)=att(6)*rhat+r*(att(3)*ralphat+att(4)*rdelhat)
  ELSE
     STOP
  END IF
! assign output
  car_att=el
! aberration correction
  tobs=el%t+r/vlight
  CALL noniter_obs(tobs,el%obscode,xobs,target_center)
! change center: heliocentric, ecliptic/ geocentric, equatorial pos. at epoch
  car_att%center=target_center
  car_att%coord=x+xobs !
  car_att%coo='CAR'
! Obscode does not change for the Earth case
  car_att%obscode=500
! h magnitude is already absolute
  IF(PRESENT(del))THEN
! for the earth case
     IF(target_center.eq.3)THEN
        del=0.d0
        del(1:3,1)=ralphat*r
        del(1:3,2)=rdelhat*r
        del(1:3,5)=rhat
        del(4:6,1)=att(6)*ralphat+att(5)*(att(3)*ralpalp+att(4)*ralpdel)
        del(4:6,2)=att(6)*rdelhat+att(5)*(att(3)*ralpdel+att(4)*rdeldel)
        del(4:6,3)=del(1:3,1) ! commutation of d/dt
        del(4:6,4)=del(1:3,2) ! commutation of d/dt
        del(4:6,5)=att(3)*ralphat+att(4)*rdelhat
        del(4:6,6)=del(1:3,5) ! commutation of d/dt
     ELSE
        del=0.d0
        CALL prodmv(del(1:3,1),roteqec,ralphat*r)
        CALL prodmv(del(1:3,2),roteqec,rdelhat*r)
        CALL prodmv(del(1:3,5),roteqec,rhat)
        CALL prodmv(del(4:6,1),roteqec,att(6)*ralphat+att(5)*(att(3)*ralpalp+att(4)*ralpdel))
        CALL prodmv(del(4:6,2),roteqec,att(6)*rdelhat+att(5)*(att(3)*ralpdel+att(4)*rdeldel))
        del(4:6,3)=del(1:3,1) ! commutation of d/dt
        del(4:6,4)=del(1:3,2) ! commutation of d/dt
        CALL prodmv(del(4:6,5),roteqec,att(3)*ralphat+att(4)*rdelhat)
        del(4:6,6)=del(1:3,5) ! commutation of d/dt
     END IF
  ENDIF
END FUNCTION car_att
!=============================
! ATTELEMENTS
! ATT elements from attributable and r,rdot
!=============================
SUBROUTINE attelements(att,r,rdot,elatt,unc)
  USE station_coordinates, ONLY: statcode
  USE attributable
  USE reference_systems, ONLY: observer_position
  TYPE(attrib), INTENT(IN):: att ! 4-dim attributable
  DOUBLE PRECISION, INTENT(IN) :: r,rdot ! range, range rate
  TYPE(orbit_elem), INTENT(OUT):: elatt ! elements of type ATT
  TYPE(orb_uncert), INTENT(OUT), OPTIONAL :: unc ! covar/normal matrices
! end interface
  TYPE(orbit_elem) xcar
  INTEGER fail_flag,indp
  DOUBLE PRECISION rsun,phase,cosph, xx(3),xpla(6),vsize,prscal,appmag,ws(4)
  DOUBLE PRECISION pos(3),vel(3),xx1(3),xpla1(6)
!=============================
  elatt=undefined_orbit_elem
  elatt%center=3 ! ???????????????????????????
  elatt%coord(1:4)=att%angles
  elatt%coord(5)=r
  elatt%coord(6)=rdot
  elatt%coo='ATT'
  elatt%t=att%tdtobs-r/vlight
  CALL statcode(att%obscod, elatt%obscode)
! cartesian coordinates needed to get phase right
  CALL coo_cha(elatt,'CAR',xcar,fail_flag)
  xx=xcar%coord(1:3)
  CALL earcar(att%tdtobs,xpla,1)

  IF(rhs.eq.1)THEN
! xx1 is heliocentric ecliptic position of asteroid?????
     xx1=xx
  ELSEIF(rhs.eq.2)THEN
! xx1 is heliocentric ecliptic position of debris??????????
     CALL prodmv(xx1,roteqec,xx)
     xx1=xx1+xpla(1:3)
  ELSE
     WRITE(*,*)'attelements: rhs=', rhs
     STOP
  END IF
  rsun=vsize(xx1)
! xpla1 is the heliocentric/ecliptic observer position
  IF(rhs.eq.1)THEN
     xpla1=xpla
  ELSEIF(rhs.eq.2)THEN
     CALL observer_position(att%tdtobs,pos,vel,OBSCODE=elatt%obscode)
     pos=MATMUL(roteqec,pos)
     vel=MATMUL(roteqec,vel)
     xpla1(1:3)=pos+xpla(1:3)
     xpla1(4:6)=vel+xpla(4:6)
  ELSE
     WRITE(*,*)'attelements: rhs=', rhs
     STOP
  END IF
  cosph=prscal(xx1,xx1-xpla1(1:3))/(rsun*r)
  IF(cosph.gt.1.d0)cosph=1.d0
  IF(cosph.lt.-1.d0)cosph=-1.d0
  phase=acos(cosph)
  IF(att%apm.gt.0.d0)THEN
     elatt%mag_set=.true.
     elatt%h_mag=att%apm-appmag(0.d0,elatt%g_mag,rsun,r,phase)
  ELSE
     elatt%mag_set=.false.
  ENDIF
  IF(PRESENT(unc))THEN
     CALL undefined_orb_uncert(6,unc)
     unc%g=0.d0
     unc%c=0.d0
     unc%g(1:4,1:4)=att%g
     unc%succ=.true.
     CALL tchinv(att%g,4,unc%c(1:4,1:4),ws,indp)
  ENDIF
END SUBROUTINE attelements

! ================================================
! ATT_PRELIM
! attributable preliminary orbit
! ================================================
SUBROUTINE att_prelim2(name0,attr,elk,uncatt,fail)
  USE attributable
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
  IF(verb_prelim.gt.10)THEN
     WRITE(iun_log,100)name0,rr,ecc
100  FORMAT(A,1X,'rho=',f8.3,1X,'ecc=',f8.4)
  ENDIF
END SUBROUTINE att_prelim2


! ====================================================================
! END ATT submodule
! ====================================================================
! ====================================================================
! OPIK submodule
! ====================================================================
! TPCAR_MTPCAR from perihelion position and velocity, obtain position on the
! asymptote with rotation by -gamma/2 and rescaling by b/d of position,
! d/b of velocity
! ====================================================================
TYPE(orbit_elem) FUNCTION tpcar_mtpcar(el,fail_flag,bd,del)
  TYPE(orbit_elem), INTENT(IN) :: el ! must be cartesian
  INTEGER, INTENT(OUT) :: fail_flag
  DOUBLE PRECISION, INTENT(OUT), OPTIONAL :: bd ! b/d rescaling factor
  DOUBLE PRECISION, DIMENSION(6,6), OPTIONAL :: del ! jacobian matrix
! =============END INTERFACE==========================================
  DOUBLE PRECISION gm, v, d, eps, sig0, vsize, prscal
  DOUBLE PRECISION bsd, bsd2, sinhalfgam,coshalfgam
  DOUBLE PRECISION dbsddd, dbsddv, dgdd, dgdv ! derivatives of b/d, sinhalfgamma w.r. to d,v
  DOUBLE PRECISION, DIMENSION(3) :: x,y, ang, xtp,ytp, xmtp,xhat,yhat, xp,yp, xphat,yphat
  DOUBLE PRECISION, DIMENSION(3,3) :: rottp,drotpdg2 ! R, dr/d(sin(gam/2)
  DOUBLE PRECISION, DIMENSION(6,6) :: dlbd, dxpdx ! derivatives of rescaling,
                          ! of rotation by -halfgamma
  DOUBLE PRECISION, DIMENSION(3,3,3) :: dtaudJ ! d(rotation)/d(ang)
  DOUBLE PRECISION, DIMENSION(3,3,3) :: ericci ! Ricci tensor
  DOUBLE PRECISION, DIMENSION(3,3)   :: dbJdX
  DOUBLE PRECISION, DIMENSION(3,3)   :: dbJdY
  DOUBLE PRECISION                   :: normJ
  INTEGER i,j,k,h ! =1,3
  IF(el%coo.ne.'CAR')THEN
     WRITE(*,*)' tpcar_mtpcar: wrong coords; center ',el%coo,' ', el%center
     STOP
  ENDIF
  fail_flag=0
  gm=centerbody_mass(el%center)
  tpcar_mtpcar=el ! copy defaults
  x=el%coord(1:3) ! position vector
  y=el%coord(4:6) ! velocity vector
!  radius and velocity
  v=vsize(y)      ! velocity
  xmtp=x-prscal(x,y)*y/v**2
  d=vsize(xmtp)      ! distance from center
  eps=2.d-6
! check on orthogonality
  sig0=prscal(x,y)/(v*d)
  IF(abs(sig0).gt.eps)THEN
     WRITE(*,199) sig0,eps
 199 FORMAT(' tpcar_mtpcar: pos, vel not orthogonal ',1P,2d10.3)
  ENDIF
  call prvec(x,y,ang)           !  angular momentum
! conformal change in the TP
  bsd2=(d*v**2)/(d*v**2-2*gm)
  IF(bsd2.gt.0.d0)THEN
     bsd=sqrt(bsd2)
  ELSE
     fail_flag=6 ! means the orbit is not hyperbolic
     RETURN
  ENDIF
  sinhalfgam=gm/(v**2*d-gm)
  coshalfgam=sqrt(1.d0-sinhalfgam**2)
! halfgamma=asin(sinhalfgam)
  normJ=vsize(ang)
  ang=ang/normJ
! coordinates rotated by -halfgamma around ang
  CALL rotmhalfgam(sinhalfgam,ang,rottp,drotpdg2,dtaudJ)
  xp=MATMUL(rottp,x)
  yp=MATMUL(rottp,y)
!  xp=x
!  yp=y
! rescaling to asymptote
  xtp=xp*bsd
  ytp=yp*(1.d0/bsd)
  tpcar_mtpcar%coord(1:3)=xtp
  tpcar_mtpcar%coord(4:6)=ytp
!  tpcar_mtpcar%coord(1:3)=xp
!  tpcar_mtpcar%coord(4:6)=yp
  IF(PRESENT(bd))THEN
    bd=bsd
  ENDIF
  tpcar_mtpcar%coo='TPC'
! partial derivatives
  IF(PRESENT(del))THEN
! derivatives d(xp,yp)/d(x,y)
! neglecting derivatives with respect to ang
     xhat=x*(1.d0/d)
     yhat=y*(1.d0/v)
     dgdd=-v**2/gm*sinhalfgam**2
     dgdv=-2.d0*v*d/gm*sinhalfgam**2
     dxpdx=0.d0
     dxpdx(1:3,1:3)=rottp
     dxpdx(4:6,4:6)=rottp
     DO i=1,3
       DO j=1,3
         DO k=1,3
             dxpdx(i,j)=dxpdx(i,j)+drotpdg2(i,k)*x(k)*dgdd*xhat(j)
             dxpdx(i+3,j+3)=dxpdx(i+3,j+3)+drotpdg2(i,k)*y(k)*dgdv*yhat(j)
             dxpdx(i,j+3)=dxpdx(i,j+3)+drotpdg2(i,k)*x(k)*dgdv*yhat(j)
             dxpdx(i+3,j)=dxpdx(i+3,j)+drotpdg2(i,k)*y(k)*dgdd*xhat(j)
         ENDDO
       ENDDO
     ENDDO
! derivatives of rottp w.r. to ang
     ericci(1:3,1:3,1:3)=0.d0
     ericci(1,2,3)=1.d0
     ericci(2,3,1)=1.d0
     ericci(3,1,2)=1.d0
     ericci(1,3,2)=-1.d0
     ericci(3,2,1)=-1.d0
     ericci(2,1,3)=-1.d0
     DO j=1,3
        DO h=1,3
         dbJdX(h,j)=1.d0/normJ*(DOT_PRODUCT(ericci(j,1:3,h),y)-v/d*x(j)*ang(h))
         dbJdY(h,j)=1.d0/normJ*(DOT_PRODUCT(ericci(1:3,j,h),x)-d/v*y(j)*ang(h))
      ENDDO
     ENDDO
! add term depending upon d(rottp0/d(ang)
      DO i=1,3
        DO j=1,3
          DO k=1,3
             dxpdx(i,j)=dxpdx(i,j)+DOT_PRODUCT(dtaudJ(i,k,1:3),dbJdX(1:3,j))*x(k)
             dxpdx(i,j+3)=dxpdx(i,j+3)+DOT_PRODUCT(dtaudJ(i,k,1:3),dbJdY(1:3,j))*x(k)
             dxpdx(i+3,j)=dxpdx(i+3,j)+DOT_PRODUCT(dtaudJ(i,k,1:3),dbJdX(1:3,j))*y(k)
             dxpdx(i+3,j+3)=dxpdx(i+3,j+3)+DOT_PRODUCT(dtaudJ(i,k,1:3),dbJdY(1:3,j))*y(k)
          ENDDO
        ENDDO
     ENDDO
! derivatives d(txp,typ)/d(xp,yp)
     xphat=xp*(1.d0/d)
     yphat=yp*(1.d0/v)
     dbsddd=1.d0/(2.d0*bsd)*(-2.d0*gm*v**2)/(v**2*d-2.d0*gm)**2
     dbsddv=1.d0/(2.d0*bsd)*(-4.d0*gm*v*d)/(v**2*d-2.d0*gm)**2
     dlbd=0.d0
     DO i=1,3
       dlbd(i,i)=bsd
       dlbd(3+i,3+i)=1.d0/bsd
       DO j=1,3
         dlbd(i,j)=dlbd(i,j)+dbsddd*xphat(j)*xp(i)
         dlbd(i,j+3)=dbsddv*yphat(j)*xp(i)
         dlbd(i+3,j)=(-1.d0/bsd**2)*dbsddd*xphat(j)*yp(i)
         dlbd(i+3,j+3)=dlbd(i+3,j+3)+(-1.d0/bsd**2)*dbsddv*yphat(j)*yp(i)
       ENDDO
     ENDDO
     del=MATMUL(dlbd,dxpdx)
!     del=dlbd
!     del=dxpdx
  ENDIF
CONTAINS
SUBROUTINE rotmhalfgam(sgmez,bJ,tau,dtaudsgm,dtaudJ)
  DOUBLE PRECISION,INTENT(IN) :: sgmez !sinus(gamma/2)
  DOUBLE PRECISION,DIMENSION(3),INTENT(IN) :: bJ !ang. momentum unit vector
! rotation of -gamma/2 around bJ and its derivative w.r.t. sinus(gamma/2)
  DOUBLE PRECISION,DIMENSION(3,3),INTENT(OUT) :: tau,dtaudsgm
! derivative of tau w.r.t. the components of bJ
! dtaudJ(*,*,k) is the derivative w.r.t. the k-th component
  DOUBLE PRECISION,DIMENSION(3,3,3),INTENT(OUT) :: dtaudJ
! ------------- end interface -------------------------------------------
  DOUBLE PRECISION :: normJ,HRoy
  DOUBLE PRECISION :: cI,sI,cO,sO, Incl, Omeg, gmez
  DOUBLE PRECISION :: cgmez,tgmez ! cosinus and tangent of gamma/2
! rotation matrices
  DOUBLE PRECISION, DIMENSION(3,3) :: ROmeg
  DOUBLE PRECISION, DIMENSION(3,3) :: RmOmeg
  DOUBLE PRECISION, DIMENSION(3,3) :: RIncl
  DOUBLE PRECISION, DIMENSION(3,3) :: RmIncl
  DOUBLE PRECISION, DIMENSION(3,3) :: Rmhgam
! for the derivatives
  DOUBLE PRECISION,DIMENSION(3) :: dHRoydJ
  DOUBLE PRECISION, DIMENSION(3,3) :: delta
  DOUBLE PRECISION, DIMENSION(3,3) :: dRmhgam
  DOUBLE PRECISION, DIMENSION(3,3) :: RORI,RmIRmO,tauauxr,tauauxl !auxiliary
  DOUBLE PRECISION, DIMENSION(3,3,3) :: dROmegdJ,dRmOmegdJ
  DOUBLE PRECISION, DIMENSION(3,3,3) :: dRIncldJ,dRmIncldJ
  DOUBLE PRECISION, DIMENSION(3,3,3) :: tensor1,tensor2,tensor3,tensor4
  DOUBLE PRECISION :: vsize !functions
  INTEGER :: h,k ! loop indexes
! ---------------------------------------------
  normJ=vsize(bJ)
  IF(abs(normJ-1.d0).gt.100*epsilon(1.d0))THEN
     WRITE(*,*) 'rotmhalfgam: angular momentum not normalized! normJ=',normJ
  ENDIF
  HRoy=sqrt(bJ(1)**2+bJ(2)**2)
  IF(HRoy/normJ.lt.epsilon(1.d0)*100)THEN
     cI=1.d0
     sI=0.d0
     cO=1.d0
     sO=0.d0
  ELSE
     cI = bJ(3) ! bJ(3)/normJ (but J is unitary)
     sI = HRoy  ! HRoy/normJ
     cO = -bJ(2)/HRoy ! singular for zero I
     sO = bJ(1)/HRoy ! singular for zero I
  ENDIF
  cgmez = sqrt(1.d0-sgmez**2)
  tgmez = sgmez/cgmez

  Incl=atan2(sI,cI)
  Omeg=atan2(sO,cO)
  gmez=atan2(sgmez,cgmez)

! *******************************************************
! rotation matrix tau = ROmeg*RIncl*Rmhgam*RmIncl*RmOmeg
! *******************************************************
  CALL rotmt(-Omeg,ROmeg,3)
  CALL rotmt(Omeg,RmOmeg,3)
  CALL rotmt(-Incl,RIncl,1)
  CALL rotmt(Incl,RmIncl,1)
  CALL rotmt(gmez,Rmhgam,3)
  RORI=matmul(ROmeg,RIncl)
  RmIRmO=matmul(RmIncl,RmOmeg)
  tauauxr = matmul(Rmhgam,RmIRmO)
  tauauxl = matmul(RORI,Rmhgam) ! will be used later
  tau = matmul(RORI,tauauxr)

! *********************************************
! derivative of Rmhgam w.r.t. sgmez
! dtaudsgm = ROmeg*RIncl*dRmhgam*RmIncl*RmOmeg
! *********************************************
  dRmhgam(1:3,1:3)=0.d0! initialization
  dRmhgam(1,1)=-tgmez
  dRmhgam(1,2)=1.d0
  dRmhgam(2,1)=-1.d0
  dRmhgam(2,2)=-tgmez
  dtaudsgm = matmul(dRmhgam, RmIRmO) ! temporary step
  dtaudsgm = matmul(RORI, dtaudsgm)

! *************************************************
! derivative of Rmhgam w.r.t. the components of bJ
! *************************************************
!  IF(HRoy/normJ.lt.epsilon(1.d0)*100)THEN
     dHRoydJ(1) = sO
     dHRoydJ(2) = -cO
!  ELSE
!     dHRoydJ(1) = bJ(1)/HRoy ! singular for zero I
!     dHRoydJ(2) = bJ(2)/HRoy ! singular for zero I
!  ENDIF
  dHRoydJ(3) = 0.d0

  CALL eye(3,delta) ! Kronecker's delta

  dROmegdJ(1:3,1:3,1:3)=0.d0 ! initialization
  IF(HRoy/normJ.ge.epsilon(1.d0)*100)THEN
     dROmegdJ(1,1,1:3)= -(1.d0/HRoy)*delta(2,1:3) + (bJ(2)/HRoy**2)*dHRoydJ(1:3) ! singular for zero I
     dROmegdJ(1,2,1:3)= -(1.d0/HRoy)*delta(1,1:3) + (bJ(1)/HRoy**2)*dHRoydJ(1:3) ! singular for zero I
  ELSE
!     dROmegdJ(1,1,1:3)= -(1.d0/HRoy)*delta(2,1:3) +  (-cO/HRoy)*dHRoydJ(1:3) ! singular for zero I
!     dROmegdJ(1,2,1:3)= -(1.d0/HRoy)*delta(1,1:3) +  (sO/HRoy)*dHRoydJ(1:3) ! singular for zero I
  ENDIF
  dROmegdJ(2,1,1:3)= -dROmegdJ(1,2,1:3)
  dROmegdJ(2,2,1:3)=  dROmegdJ(1,1,1:3)
  dRmOmegdJ(1:3,1:3,1:3)= dROmegdJ(1:3,1:3,1:3)
  dRmOmegdJ(1,2,1:3)= -dROmegdJ(1,2,1:3)
  dRmOmegdJ(2,1,1:3)= -dRmOmegdJ(1,2,1:3)
  dRIncldJ(1:3,1:3,1:3)=0.d0 ! initialization
  dRIncldJ(2,2,1:3)=  delta(3,1:3)
  dRIncldJ(2,3,1:3)= -dHRoydJ(1:3)
  dRIncldJ(3,2,1:3)= -dRIncldJ(2,3,1:3)
  dRIncldJ(3,3,1:3)=  dRIncldJ(2,2,1:3)

  dRmIncldJ(1:3,1:3,1:3)= dRIncldJ(1:3,1:3,1:3)
  dRmIncldJ(2,3,1:3)= -dRIncldJ(2,3,1:3)
  dRmIncldJ(3,2,1:3)= -dRmIncldJ(2,3,1:3)

! **********************************************
! tensor1 = dROmegdJ*RIncl*Rmhgam*RmIncl*RmOmeg
! **********************************************
  DO k=1,3
     tensor1(1:3,1:3,k) = matmul(RIncl,tauauxr)
     tensor1(1:3,1:3,k) = matmul(dROmegdJ(1:3,1:3,k),tensor1(1:3,1:3,k))
  ENDDO
! **********************************************
! tensor2 = ROmeg*dRIncldJ*Rmhgam*RmIncl*RmOmeg
! **********************************************
  DO k=1,3
     tensor2(1:3,1:3,k) = matmul(dRIncldJ(1:3,1:3,k),tauauxr)
     tensor2(1:3,1:3,k) = matmul(ROmeg,tensor2(1:3,1:3,k))
  ENDDO
! **********************************************
! tensor3 = ROmeg*RIncl*Rmhgam*dRmIncldJ*RmOmeg
! **********************************************
  DO k=1,3
     tensor3(1:3,1:3,k) = matmul(tauauxl,dRmIncldJ(1:3,1:3,k))
     tensor3(1:3,1:3,k) = matmul(tensor3(1:3,1:3,k),RmOmeg)
  ENDDO
! **********************************************
! tensor4 = ROmeg*RIncl*Rmhgam*RmIncl*dRmOmegdJ
! **********************************************
  DO k=1,3
     tensor4(1:3,1:3,k) = matmul(tauauxl,RmIncl)
     tensor4(1:3,1:3,k) = matmul(tensor4(1:3,1:3,k),dRmOmegdJ(1:3,1:3,k))
  ENDDO

  dtaudJ = tensor1 + tensor2 + tensor3 + tensor4

END SUBROUTINE rotmhalfgam

END FUNCTION tpcar_mtpcar
! ====================================================================
! OPIK_TPCAR from perihelion position and velocity, obtain Opik elements
! ====================================================================
! opik= U, long, lat, csi, zeta, eta
! U= unperturbed v_infty
! long, lat  = antiradians (in ecliptical coordinates)
! csi, zeta = TP coordinates with axis zeta along projection of normal to
!                   ecliptic
! eta = U(t-t0), t0 time of pericenter passage
! WARNING: in .clo file we need also the time dependent theta, phi
! to be computed using the Earth velocity at t0
! ====================================================================
TYPE(orbit_elem) FUNCTION opik_tpcar(el,fail_flag,del,vt3)
  TYPE(orbit_elem), INTENT(IN) :: el ! must be cartesian
  INTEGER, INTENT(OUT) :: fail_flag
  DOUBLE PRECISION, DIMENSION(6,6), OPTIONAL :: del ! jacobian matrix
  DOUBLE PRECISION, DIMENSION(3,3), OPTIONAL :: vt3 ! rotation of TP to plane eta=0
! =============END INTERFACE==========================================
  DOUBLE PRECISION opik(6), xytp(6), vsize, b, dz, prscal, cosa, dt
  DOUBLE PRECISION v1(3), v2(3), rv(3,3),xx(3),drv(3,3,3)!,rvt(3,3), rv1(3,3)
  INTEGER i,j
  IF(el%coo.ne.'TPC')THEN
     WRITE(*,*)' opik_tpcar: wrong coordinates ', el
     STOP
  ENDIF
  fail_flag=0
  opik_tpcar=el ! set default; note that epoch remains the one of pericenter
  xytp=el%coord
  opik(1)=vsize(xytp(4:6))
  v1=xytp(4:6)*(1.d0/opik(1))
  b=vsize(xytp(1:3))
  cosa=prscal(xytp(1:3),xytp(4:6))/(b*opik(1))
  dt=-b*cosa/opik(1)
  IF(abs(cosa).gt.2.d-6)THEN
!     xytp(1:3)=xytp(1:3)-b*cosa*v1
!     b=vsize(xytp(1:3))
     WRITE(*,199)cosa, b, opik(1)
 199 FORMAT('opik_tpcar: pos not orthogonal to vel ',1p,3D10.3)
!     fail_flag=7
  ENDIF
! polar coordinates
  dz=xytp(4)**2+xytp(5)**2
  if (dz.le.100*epsilon(1.d0)) then
! remove singularity at poles
     opik(2)=0.d0
  else
! ecliptic longitude of velocity (radians)
     opik(2)=atan2(xytp(5),xytp(4))
     if (opik(2).lt.0.d0) then
        opik(2)=opik(2)+dpig
     endif
  endif
! ecliptic latitude of velocity (radians)
  opik(3)=asin(xytp(6)/opik(1))
! target plane
  CALL drvdytp(xytp(4:6),drv,rv)
  IF(PRESENT(vt3))THEN
     vt3=rv
  ENDIF
  v1=xytp(4:6)*(1.d0/opik(1))
!  v2=(/0.d0, 0.d0, 1.d0/)
!  CALL mtp_ref(v1,v2,rvt,rv)
  xx=MATMUL(rv,xytp(1:3))
  opik(4)=xx(2)
  opik(5)=xx(3)
  IF(abs(xx(1)).gt.2.d-7)THEN
     WRITE(*,*)' opik_tpcar: pos, vel not orthogonal ', xx(1), xx(1)/b
  ENDIF
! last element is the time of pericenter passage
  opik(6)=opik(1)*dt
! resulting elements
  opik_tpcar%t=el%t
  opik_tpcar%coord=opik
  IF(PRESENT(del))THEN
! partial derivatives
     del(1:3,1:3)=0.d0
     del(1,4:6)=v1
     del(2,4)=-xytp(5)/dz
     del(2,5)=xytp(4)/dz
     del(2,6)=0.d0
     del(3,4)=-xytp(6)*(xytp(4)/(sqrt(dz)*opik(1)**2))
     del(3,5)=-xytp(6)*(xytp(5)/(sqrt(dz)*opik(1)**2))
     del(3,6)=sqrt(dz)/opik(1)**2
     del(4:5,4:6)=0.d0 ! neglecting derivatives of rv
     del(4:5,1:3)=rv(2:3,1:3)
     del(6,1:3)=-(1.d0/opik(1))*xytp(4:6)
     del(6,4:6)=-(1.d0/opik(1))*xytp(1:3)+b/opik(1)*cosa*v1
! contribution from derivatives of rv
     DO i=4,5
       DO j=4,5
         del(i,j)=del(i,j)+DOT_PRODUCT(drv(i-2,1:3,j-3),xytp(1:3))
       ENDDO
     ENDDO
  ENDIF
  CONTAINS

SUBROUTINE drvdytp(ytp,drv,rv)
DOUBLE PRECISION, DIMENSION(3), INTENT(IN) :: ytp
                  ! Velocity vector in TP reference frame
DOUBLE PRECISION, DIMENSION(3,3), INTENT(OUT) :: rv
                  ! rotation matrix
DOUBLE PRECISION, DIMENSION(3,3,3), INTENT(OUT) :: drv
                  ! Derivatives of rotation matrix rv
DOUBLE PRECISION :: U,sqrx2y2,trhalf,d1,d2, x,y,z
                  ! Auxiliary quantities
! ==================================================================
 U=sqrt(ytp(1)**2+ytp(2)**2+ytp(3)**2)          ! Module of velocity
 x=ytp(1)
 y=ytp(2)
 z=ytp(3)
! Matrix rv
 rv(1,1)=x/U
 rv(1,2)=y/U
 rv(1,3)=z/U
 rv(3,1)=-x*z/(U*sqrt(x**2+y**2))
 rv(3,2)=-y*z/(U*sqrt(x**2+y**2))
 rv(3,3)=(x**2+y**2)/(U*sqrt(x**2+y**2))
 rv(2,1)=y/sqrt(x**2+y**2)
 rv(2,2)=-x/sqrt(x**2+y**2)
 rv(2,3)=0.d0

! ==================================================================

sqrx2y2=sqrt(ytp(1)**2+ytp(2)**2)
trhalf=sqrx2y2*(ytp(1)**2+ytp(2)**2)
d1=U**3*trhalf
d2=U**3*sqrx2y2
! ===================================================
! Derivatives of first row of matrix R_v
drv(1,1,1)=(ytp(2)**2+ytp(3)**2)/U**3
drv(1,1,2)=-ytp(1)*ytp(2)/U**3
drv(1,1,3)=-ytp(1)*ytp(3)/U**3

drv(1,2,1)=-ytp(1)*ytp(2)/U**3
drv(1,2,2)=(ytp(1)**2+ytp(3)**2)/U**3
drv(1,2,3)=-ytp(2)*ytp(3)/U**3

drv(1,3,1)=-ytp(1)*ytp(3)/U**3
drv(1,3,2)=-ytp(2)*ytp(3)/U**3
drv(1,3,3)=(ytp(1)**2+ytp(2)**2)/U**3
! ===================================================
! Derivatives of third row of matrix R_v
drv(3,1,1)=ytp(3)*((ytp(1)**2-U**2)*(ytp(1)**2+ytp(2)**2)+U**2*ytp(1)**2)/d1
drv(3,1,2)=ytp(1)*ytp(2)*ytp(3)*(ytp(1)**2+ytp(2)**2+U**2)/d1
drv(3,1,3)=ytp(1)*(ytp(3)**2-U**2)/d2

drv(3,2,1)=ytp(1)*ytp(2)*ytp(3)*(ytp(1)**2+ytp(2)**2+U**2)/d1
drv(3,2,2)=ytp(3)*((ytp(2)**2-U**2)*(ytp(1)**2+ytp(2)**2)+U**2*ytp(2)**2)/d1
drv(3,2,3)=ytp(2)*(ytp(3)**2-U**2)/d2

drv(3,3,1)=ytp(1)*ytp(3)**2/d2
drv(3,3,2)=ytp(2)*ytp(3)**2/d2
drv(3,3,3)=-ytp(3)*sqrx2y2/U**3
! ====================================================
! Derivatives of second row of matrix R_v
drv(2,1,1)=-ytp(1)*ytp(2)/trhalf
drv(2,1,2)=ytp(1)**2/trhalf
drv(2,1,3)=0.d0

drv(2,2,1)=-ytp(2)**2/trhalf
drv(2,2,2)=ytp(1)*ytp(2)/trhalf
drv(2,2,3)=0.d0

drv(2,3,1)=0.d0
drv(2,3,2)=0.d0
drv(2,3,3)=0.d0
! =====================================================
END SUBROUTINE drvdytp


END FUNCTION opik_tpcar
! ====================================================================
! END OPIK submodule
! ====================================================================
! ====================================================================
! I/O submodule
! ====================================================================
SUBROUTINE read_elems_old(el,name,eof,file,unit,covar,dp,nd)
  USE io_elems, ONLY: obscod
  USE station_coordinates, ONLY: statcode
  TYPE(orbit_elem), INTENT(OUT) :: el
  CHARACTER*(*),INTENT(OUT) :: name
  LOGICAL, INTENT(OUT) :: eof
  CHARACTER*(*), INTENT(IN) :: file ! input file
  INTEGER, INTENT(IN), OPTIONAL :: unit ! input unit, only if file already opened
  DOUBLE PRECISION,INTENT(OUT),DIMENSION(ndyx),OPTIONAL :: dp ! non grav parameters
  INTEGER,INTENT(IN),OPTIONAL :: nd ! number of solve for parameters
  TYPE(orb_uncert), INTENT(OUT), OPTIONAL :: covar
! end interface
  DOUBLE PRECISION elem(6),t0,h,g,cove(ndimx,ndimx),nore(ndimx,ndimx),mass !,cnv(6)
  CHARACTER*3 eltype
  CHARACTER*10 rsys,epoch
  LOGICAL defcov,defnor,defcn
  INTEGER kr, iobs
  INTEGER nd_cur ! number of solve for parameters (either nd if present, or 6)

  IF(PRESENT(nd))THEN
     nd_cur=nd
  ELSE
     nd_cur=6
  ENDIF

! if present dp, then nd is required
  IF(PRESENT(dp).and..not.PRESENT(nd))THEN
     WRITE(*,*)'read_elems_old: nd not present'
     STOP
  ENDIF

  IF(.not.PRESENT(unit))THEN
     CALL oporbf(file,0)
  ENDIF
  el=undefined_orbit_elem
  IF(rhs.EQ.2)THEN
     el%center=3
  END IF
  IF(PRESENT(dp))THEN
     dp=0.d0
     CALL rd_orb(name,elem,eltype,t0,cove,defcov,nore,defnor,     &
          &                 h,g,mass,rsys,epoch,kr,eof,dp,nd)
  ELSE
     CALL rd_orb(name,elem,eltype,t0,cove,defcov,nore,defnor,     &
          &                 h,g,mass,rsys,epoch,kr,eof)
  ENDIF
  IF(eof)THEN
     IF(.not.PRESENT(unit)) CALL clorbf
     RETURN
  ENDIF
! unit conversions already done by rdorb
  IF(eltype.eq.'ATT')THEN
! input obscod
     CALL statcode(obscod,iobs)
  END IF
! what to do with rsys,epoch????
  IF(rsys.ne.'ECLM'.or.epoch.ne.'J2000')THEN
     WRITE(*,*)' read_elems_old: other reference systems not handled ', rsys,'  ', epoch
     STOP
  ENDIF
! copy into output elements
  el%coo=eltype
  el%coord=elem
  el%t=t0
  IF(h.lt.-10.d0)THEN
     el%mag_set=.false.
  ELSE
     el%mag_set=.true.
     el%h_mag=h
     el%g_mag=g
  ENDIF
  IF(el%coo.eq.'ATT')el%obscode=iobs

  IF(PRESENT(covar))THEN
     CALL undefined_orb_uncert(nd_cur,covar)
     IF(defcov) covar%g(1:nd_cur,1:nd_cur)=cove(1:nd_cur,1:nd_cur)
     IF(defnor) covar%c(1:nd_cur,1:nd_cur)=nore(1:nd_cur,1:nd_cur)
     CALL fixcnm(defcov,defnor,defcn,cove,nore)
     IF(defcn)THEN
        covar%succ=.true.
     ELSE
        covar%succ=.false.
     ENDIF
  ENDIF

  IF(.not.PRESENT(unit)) CALL clorbf
END SUBROUTINE read_elems_old

SUBROUTINE write_elems_old(el,name,form,file,unit,covar,incfit,dp)
  USE name_rules
  USE io_elems, ONLY: obscod
  USE station_coordinates, ONLY: codestat
  USE dyn_param
  TYPE(orbit_elem), INTENT(IN) :: el
  CHARACTER*(*),INTENT(IN) :: name
  CHARACTER*(*), INTENT(IN) :: form ! can be either '1L' or 'ML'
  CHARACTER*(*), INTENT(IN), OPTIONAL :: file ! output file
                       ! (not used if PRESENT(unit))
  INTEGER, INTENT(IN), OPTIONAL :: unit ! input unit, only if already opened
! WARNING: if file already opened, it must also have the header already written
  TYPE(orb_uncert), INTENT(IN), OPTIONAL :: covar
  INTEGER, INTENT(IN), OPTIONAL :: incfit ! result from incomplete fit;
    ! warn the users of this!!!!
  DOUBLE PRECISION, OPTIONAL, INTENT(IN) :: dp(ndyx)
! =========== end interface================
  INTEGER :: nd
  INTEGER uniout,le
  CHARACTER*(idnamvir_len) namloc
  CHARACTER*10 rsys,epoch
  LOGICAL defcov,defnor
  DOUBLE PRECISION mass, cove(ndimx,ndimx), nore(ndimx,ndimx)
! number of parameters
  nd=6+nls
! open output file if not opened already
  IF(.not.PRESENT(unit))THEN
     CALL filopn(uniout,file,'unknown')
! write file header
     rsys='ECLM'
     epoch='J2000'
     IF(form.eq.'1L')THEN
! what to do when el%coo=ATT ?????
        CALL codestat(el%obscode,obscod)
        CALL wro1lh(uniout,rsys,epoch,el%coo)
     ELSEIF(form.eq.'ML')THEN
        CALL wromlh(uniout,rsys,epoch)
     ENDIF
  ELSE
     uniout=unit
  ENDIF
  namloc=name
  CALL rmsp(namloc,le)
  IF(PRESENT(incfit))THEN
     IF(incfit.eq.5)THEN
        WRITE(uniout,200)namloc(1:le)
200     FORMAT('! CONSTRAINED SOLUTION for ',a)
     ELSEIF(incfit.eq.4)THEN
        WRITE(uniout,201)namloc(1:le)
201     FORMAT('! 4_PARAMETERS SOLUTION for ',a)
     ENDIF
  ENDIF
! COVARIANCE; HOWEVER NOT WRITTEN IF FORMAT= '1L'
  IF(PRESENT(covar))THEN
     IF(covar%succ)THEN
        defcov=.true.
        defnor=.true.
        cove(1:nd,1:nd)=covar%g(1:nd,1:nd)
        nore(1:nd,1:nd)=covar%c(1:nd,1:nd)
     ELSE
        defcov=.false.
        defnor=.false.
        cove=0.d0
        nore=0.d0
     ENDIF
  ELSE
     defcov=.false.
     defnor=.false.
     cove=0.d0
     nore=0.d0
  ENDIF
  mass=0.d0
  IF(form.eq.'1L')THEN
     CALL wro1lr(uniout,namloc(1:le),el%coord,el%coo,el%t,el%h_mag,el%g_mag)
  ELSEIF(form.eq.'ML')THEN
     CALL codestat(el%obscode,obscod)
     IF(PRESENT(dp))THEN
        CALL wromlr(uniout,namloc(1:le),el%coord,el%coo,el%t,cove(1:nd,1:nd), &
       &     defcov,nore(1:nd,1:nd),defnor,el%h_mag,el%g_mag,mass,nd,dp)
     ELSE
        CALL wromlr(uniout,namloc(1:le),el%coord,el%coo,el%t,cove(1:nd,1:nd), &
       &     defcov,nore(1:nd,1:nd),defnor,el%h_mag,el%g_mag,mass,nd)
     ENDIF
  ENDIF

  IF(.not.PRESENT(unit)) CALL filclo(uniout,' ')

END SUBROUTINE write_elems_old

! propagation to covariance/normal matrix of conversion of coordinates,
! with a 6x6 jacobian
SUBROUTINE convertunc(unc1,dee,unc2)
! unc1 is input, unc2 is output, but they can be on the same memory location
  TYPE(orb_uncert), INTENT(IN) :: unc1
  TYPE(orb_uncert), INTENT(INOUT) :: unc2
  DOUBLE PRECISION, DIMENSION(6,6), INTENT(IN) :: dee
!*******************************************************
  LOGICAl error
  INTEGER nd

! Dimension
  nd=unc1%ndim
! covariance matrix is propagated by similarity transformation
  CALL convertcovdp(nd,unc1%g(1:nd,1:nd),dee,unc2%g(1:nd,1:nd))
! normal matrix is propagated by similarity transformation
  CALL norprsdp(unc1%c(1:nd,1:nd),dee,nd,unc2%c(1:nd,1:nd),error)
  IF(error)THEN
     unc2%succ=.false.
     unc2%ndim=unc1%ndim
  ELSE
     unc2%succ=.true.
     unc2%ndim=unc1%ndim
  ENDIF
END SUBROUTINE convertunc

! propagation of covariance/normal to a different time,
! with a nd x nd jacobian
SUBROUTINE propagunc(nd,unc1,depdep,unc2)
! unc1 is input, unc2 is output, but they can be on the same memory location
  INTEGER, INTENT(IN) :: nd
  TYPE(orb_uncert), INTENT(IN) :: unc1
  TYPE(orb_uncert), INTENT(INOUT) :: unc2
  DOUBLE PRECISION, DIMENSION(nd,nd), INTENT(IN) :: depdep
  LOGICAl error
! covariance matrix is propagated by similarity transformation
  CALL convertcov(nd,unc1%g(1:nd,1:nd),depdep(1:nd,1:nd),unc2%g(1:nd,1:nd))
! normal matrix is propagated by similarity transformation
  CALL norprs(unc1%c(1:nd,1:nd),depdep(1:nd,1:nd),nd,unc2%c(1:nd,1:nd),error)
  IF(error)THEN
     unc2%succ=.false.
     unc2%ndim=unc1%ndim
  ELSE
     unc2%succ=.true.
     unc2%ndim=unc1%ndim
  ENDIF
END SUBROUTINE propagunc

CHARACTER*3 FUNCTION cootyp(ityp)
  INTEGER, INTENT(IN) :: ityp
  if(ityp.eq.1)then
     cootyp='KEP'
  elseif(ityp.eq.2)then
     cootyp='EQU'
  elseif(ityp.eq.3)then
     cootyp='CAR'
  elseif(ityp.eq.4)then
     cootyp='COM'
  elseif(ityp.eq.5)then
     cootyp='COT'
  elseif(ityp.eq.6)then
     cootyp='ATT'
  else
     WRITE(*,*)' cootyp: coord code not known ',ityp
     STOP
  endif
END FUNCTION cootyp

! Copyright (C) 1997-1998 by Orbfit Consortium
! Version: December 16, 1998 Steven Chesley chesley@dm.unipi.it
! Version 3.1 September 2003 A.Milani
! ---------------------------------------------------------------------
!
!  *****************************************************************
!  *                                                               *
!  *                         W R I K E P                           *
!  *                                                               *
!  *  Writes an orbital element record in an orbital element file  *
!  *                     (multi-line format)                       *
!  *                                                               *
!  *****************************************************************
!
! INTPUT:   FILE      -  Output file name
!           NAME      -  Name of planet/asteroid/comet
!           ELEM      -  Orbital elements
!           COVAR     -  Uncertainty of orbital elements
!           UNIT      -  If output file already opened
!
! OUPUT:    MOID      -  MOID
!           PHA       -  0 or 1 PHA flag
!
! WARNING: the routine does write the header of the file: this
!          must not be generated by calling subroutine wromlh
!
SUBROUTINE wrikep(file,name,elem,covar,moid,pha,unit0)
  CHARACTER*(*), INTENT(INOUT) :: name
  TYPE(orbit_elem), INTENT(IN) :: elem
  TYPE(orb_uncert), INTENT(IN) :: covar
  DOUBLE PRECISION, INTENT(OUT) ::  moid
  CHARACTER*(*), INTENT(INOUT) :: file
  LOGICAL, INTENT(OUT) :: pha
   INTEGER, OPTIONAL, INTENT(IN) :: unit0
  TYPE(orbit_elem) :: elk
  TYPE(orb_uncert) :: covk,cove, unc6, unc6k
  INCLUDE 'parcmc.h90'
  INTEGER lnn,lnf,i,j,k,unit, fail_flag
  DOUBLE PRECISION cnv(ndimx),std(ndimx)
  DOUBLE PRECISION enne, dkde(6,6), correl(ndimx,ndimx)
  INTEGER :: nd

! "Other useful data"
  DOUBLE PRECISION perihe,aphe,period,dnode,anode,vinfty
  INTEGER iconv

  nd=covar%ndim
  IF(PRESENT(unit0))THEN
     unit=unit0
  ELSE
! Open FIle
     call rmsp(file,lnf)
     call filopn(unit,file(1:lnf),'unknown')
     CALL wromlh (unit,'ECLM','J2000')
  ENDIF
! Object Name
  call rmsp(name,lnn)
  WRITE(unit,100) name(1:lnn)
100 FORMAT(A)
! Get Keplerian Elements
  CALL coo_cha(elem,'KEP',elk,fail_flag,dkde)
! WARNING: not working for cometary elements
! Write Orbital elements
  cnv(1:2)=1.d0
  cnv(3:6)=degrad
  IF(nd.GT.6)THEN
     cnv(7:6+dyn%ndp)=1.d0
  END IF
  WRITE(unit,201) comcha
201 FORMAT(A,' Keplerian elements: a, e, i, long. node,',  &
     &         ' arg. peric., mean anomaly')
  WRITE(unit,101) (elk%coord(i)*cnv(i),i=1,6)
101 FORMAT(' KEP ',2F12.6,4F12.3)
! Epoch
  WRITE(unit,104) elk%t
104 FORMAT(' MJD ',F12.4,' TDT')
! Mass
!      IF(mass.NE.0.d0) WRITE(unit,105) mass
!  105 FORMAT(' MAS ',1P,E20.12)

! Magnitudes
  IF(elk%mag_set) WRITE(unit,106) elk%h_mag,elk%g_mag
106 FORMAT(' MAG ',2F7.3)
! nongrav parameters
  IF(dyn%nmod.ge.0)THEN
! ndp.gt.0 implies PRESENT(nmod) nmod.gt.0
! record giving the matrix dimension and list of solved parameters
     IF(nls.GT.0)THEN
        WRITE(unit,149) comcha
149     FORMAT(A,' Non-grav parameters: model used, actual number in use, dimension, ', &
             & 'list of solve for parameters to be currently determined')
        WRITE(unit,150) dyn%nmod, dyn%ndp, nd, ls(1:nls)
150     FORMAT(' LSP ',1X,I2,1X,I2,3X,I2,2X,4(1X,I2))
     ELSE
        WRITE(unit,'(A)') '! Non-grav parameters: model used, actual number in use, dimension'
        WRITE(unit,155) dyn%nmod, dyn%ndp, nd
155     FORMAT(' LSP ',1X,I2,1X,I2,3X,I2)
     END IF
  END IF
! NGR
  IF(nd.GT.6)THEN
     WRITE(unit,112) ((dyn%dp(i))*cnv(6+i),i=1,dyn%ndp)
112  FORMAT(' NGR ',1P,4E12.3) !4E22.12
  END IF

! A-priori for non-gravitational parameter
! ========================================
! Dynamical model nmod=1
     IF(dyn%nmod.EQ.1)THEN
! A/M and A2 to be determined
        IF(ls(1).EQ.1.AND.ls(2).EQ.2)THEN
! A/M and A2 a-priori available
           IF(dyn%apriori(2).AND.dyn%apriori(1))THEN
              WRITE(unit,'(A)') '! A priori value and standard deviation (A/M and Yarkovsky effect)'
              WRITE(unit,352) dyn%apriori_nominal(1), dyn%apriori_std(1)
              WRITE(unit,351) dyn%apriori_nominal(2), dyn%apriori_std(2)
! only A2 a-priori available
           ELSE IF(dyn%apriori(2))THEN
              WRITE(unit,'(A)') '! A priori value and standard deviation (only Yarkovsky effect)'
              WRITE(unit,351) dyn%apriori_nominal(2), dyn%apriori_std(2)
! only A/M a-priori available
           ELSE IF(dyn%apriori(1))THEN
              WRITE(unit,'(A)') '! A priori value and standard deviation (only A/M)'
              WRITE(unit,352) dyn%apriori_nominal(1), dyn%apriori_std(1)
           END IF
! only A2 to be determined with a-priori
        ELSE IF(ls(1).EQ.2.AND.dyn%apriori(2))THEN
           WRITE(unit,'(A)') '! A priori value and standard deviation (only Yarkovsky effect)'
           WRITE(unit,351) dyn%apriori_nominal(2), dyn%apriori_std(2)
! only A/M to be determined with a-priori
        ELSE IF(ls(1).EQ.1.AND.dyn%apriori(1))THEN
           WRITE(unit,'(A)') '! A priori value and standard deviation (only A/M)'
           WRITE(unit,352) dyn%apriori_nominal(1), dyn%apriori_std(1)
        END IF

! Dynamical model nmod=2
     ELSEIF(dyn%nmod.EQ.2)THEN
        IF(ls(1).EQ.1.AND.ls(2).EQ.2.AND.ls(3).EQ.3.AND.ls(4).EQ.4)THEN
           IF(dyn%apriori(1).AND.dyn%apriori(2).AND.dyn%apriori(3).AND.dyn%apriori(4))THEN
              WRITE(unit,'(A)') '! A priori value and standard deviation (A1, A2, A3 and dtdelay)'
              WRITE(unit,353) dyn%apriori_nominal(1), dyn%apriori_std(1)
              WRITE(unit,354) dyn%apriori_nominal(2), dyn%apriori_std(2)
              WRITE(unit,355) dyn%apriori_nominal(3), dyn%apriori_std(3)
              WRITE(unit,356) dyn%apriori_nominal(4), dyn%apriori_std(4)
           END IF
        ELSE IF(ls(1).EQ.1.AND.ls(2).EQ.2.AND.ls(3).EQ.3)THEN
           IF(dyn%apriori(1).AND.dyn%apriori(2).AND.dyn%apriori(3))THEN
              WRITE(unit,'(A)') '! A priori value and standard deviation (A1, A2, A3)'
              WRITE(unit,353) dyn%apriori_nominal(1), dyn%apriori_std(1)
              WRITE(unit,354) dyn%apriori_nominal(2), dyn%apriori_std(2)
              WRITE(unit,355) dyn%apriori_nominal(3), dyn%apriori_std(3)
           END IF
        ELSE IF(ls(1).EQ.1.AND.ls(2).EQ.2)THEN
           IF(dyn%apriori(1).AND.dyn%apriori(2))THEN
              WRITE(unit,'(A)') '! A priori value and standard deviation (A1, A2)'
              WRITE(unit,353) dyn%apriori_nominal(1), dyn%apriori_std(1)
              WRITE(unit,354) dyn%apriori_nominal(2), dyn%apriori_std(2)
           END IF
        END IF
     END IF
351  FORMAT(' APR YAR ',1P,2E22.14)
352  FORMAT(' APR A/M ',1P,2E22.14)
353  FORMAT(' APR A1C ',1P,2E22.14)
354  FORMAT(' APR A2C ',1P,2E22.14)
355  FORMAT(' APR A3C ',1P,2E22.14)
356  FORMAT(' APR DTD ',1P,2E22.14)

! Other useful data:
! PERIHELION 1.2345
! APHELION 2.3456
! ANODE 1.01234
! DNODE 1.04321
! MOID 0.00001
! PERIOD 1.11111
! PHA F
! VINFTY 14.5678 (in km/s)
  perihe=elk%coord(1)*(1.d0-elk%coord(2))
  aphe=elK%coord(1)*(1.d0+elk%coord(2))
  enne=sqrt(gms/elk%coord(1)**3)
  period=dpig/enne
  CALL nomoid(elem%t,elem,moid,anode,dnode)
  pha = moid.le.0.05d0.and.elk%h_mag.lt.22d0.and.elk%mag_set
  vinfty= v_infty0(elk)
  WRITE(unit,102)comcha,perihe,comcha,aphe,comcha,                  &
     &     anode,comcha,dnode,comcha,moid,comcha,period,comcha,pha, &
     &     comcha,vinfty
  102 FORMAT(A1,' PERIHELION ',F9.4/                                    &
     &       A1,' APHELION ',F9.4/                                      &
     &       A1,' ANODE ',F10.5/                                        &
     &       A1,' DNODE ',F10.5/                                        &
     &       A1,' MOID ',F10.5/                                         &
     &       A1,' PERIOD ',F16.4/                                       &
     &       A1,' PHA ',L1/                                             &
     &       A1,' VINFTY ',F9.4)
! Get Keplerian Covariance
  IF(covar%succ) THEN
     cove=covar
     IF(nd.GT.6)THEN
        unc6%c(1:6,1:6)=cove%c(1:6,1:6)
        unc6%g(1:6,1:6)=cove%g(1:6,1:6)
        unc6%succ=cove%succ
        unc6%ndim=6
        CALL convertunc(unc6,dkde,unc6k)
        covk%c(1:6,1:6)=unc6k%c(1:6,1:6)
        covk%g(1:6,1:6)=unc6k%g(1:6,1:6)
        covk%succ=unc6k%succ
        covk%g(1:nd,7:nd)=cove%g(1:nd,7:nd)
        covk%g(7:nd,1:nd)=cove%g(7:nd,1:nd)
        covk%c(1:nd,7:nd)=cove%c(1:nd,7:nd)
        covk%c(7:nd,1:nd)=cove%c(7:nd,1:nd)
     ELSE
        CALL convertunc(cove,dkde,covk)
     END IF
! Get Correlation Matrix
     DO i=1,nd
        std(i) = sqrt(covk%g(i,i))
     ENDDO
     DO i=1,nd
        DO j=1,nd
           correl(i,j) = covk%g(i,j)/std(i)/std(j)
        ENDDO
     ENDDO
! Covariance matrix
     WRITE(unit,107) comcha,(std(i)*cnv(i),i=1,nd)
     WRITE(unit,108) ((covk%g(i,k)*cnv(i)*cnv(k),k=i,nd),i=1,nd)
107  FORMAT(A1,' RMS ',1P,10E12.3)
108  FORMAT(' COV ',1P,3E23.8)
! Correlation matrix
     WRITE(unit,109)  ((correl(i,k),k=i,nd),i=1,nd)
109  FORMAT(' COR ',3F23.8)
  ENDIF
  IF(.not.PRESENT(unit0))THEN
     call filclo(unit,' ')
  ENDIF

END SUBROUTINE wrikep

! ===================================================
! velocity at infinity with repect to Earth
! (circular approx) in au/day
! this is a duplicate of the one in tp_trace, to be used by wrikep
DOUBLE PRECISION FUNCTION v_infty0(el0)
! input elements
  TYPE(orbit_elem), INTENt(IN) :: el0
! END INTERFACE
  TYPE(orbit_elem) eleq
  DOUBLE PRECISION eq0(6)
  DOUBLE PRECISION cosi,v2
  CHARACTER*3 coo
  INTEGER fail_flag
  coo=el0%coo
  CALL coo_cha(el0,'EQU',eleq,fail_flag)
  eq0=eleq%coord
  IF(fail_flag.eq.0)THEN
! v_infinity computation in equinoctal
     cosi=(1.d0-eq0(4)**2-eq0(5)**2)/                                  &
     &     (1.d0+eq0(4)**2+eq0(5)**2)
     v2=3.d0-1.d0/eq0(1) -                                             &
     &     2.d0*sqrt(eq0(1)*(1.d0-eq0(2)**2-eq0(3)**2))*cosi
  ELSEIF(fail_flag.eq.5)THEN
! v_infinity computation in cometary
     cosi=cos(eq0(3))
     v2=3.d0-(1.d0-eq0(2))/eq0(1) - 2.d0*sqrt(eq0(1)*(1.d0+eq0(2)))*cosi
  ELSE
     WRITE(*,*)' v_infty: coordinate conversion failed ', el0,fail_flag
  ENDIF
! handle non hyperbolic (w.r. to Earth) case
  IF(v2.gt.0.d0)THEN
     v_infty0=sqrt(v2)
  ELSE
     v_infty0=0.d0
  ENDIF
! normalization in au/day by Gauss constant
  v_infty0=v_infty0*gk
! conversion to km/s required by wrikep
   v_infty0=v_infty0*aukm/86400.d0
END FUNCTION v_infty0

! ===================================
! FIND_ELEMS
! find elements in .eq0 file, if necessary convert into
! required coordinate system, and report on availability
!  ===================================
SUBROUTINE find_elems(nam0,eledir,buildpath,coox,iunout,el0,unc0,ini0,cov0)
  CHARACTER*(*),INTENT(IN) :: nam0 ! asteroid name, could be long
  CHARACTER*60, INTENT(IN) ::  eledir ! root directory for .eq0 files
  LOGICAL, INTENT(IN):: buildpath ! to compose path with eledir, if not use
                                  ! eledir as complete pathname
  CHARACTER*3, INTENT(IN) :: coox ! in which coords the elements are required
   ! if coox=' ' coords are not changed, but left as they are
  INTEGER,INTENT(IN) :: iunout ! log output unit
  TYPE(orbit_elem), INTENT(OUT) :: el0 ! epoch time (MJD), elements, abs. mag.
  TYPE(orb_uncert), INTENT(OUT) :: unc0 ! covariance and normal matrices
  LOGICAL,INTENT(OUT) :: ini0,cov0 ! successful input flag, covar. available
! === END INTERFACE =====================
  INCLUDE 'sysdep.h90'
  LOGICAL eof, error, err
  INTEGER nlsloc, lsloc(ndyx)
  INTEGER nd
  CHARACTER*60 elefi0 ! elements  file name
  INTEGER le, fail_flag
  CHARACTER*19 name
  CHARACTER*120 file
  TYPE(orbit_elem) :: elk
  TYPE(orb_uncert) :: unck
  DOUBLE PRECISION :: dee(6,6)

!
  ini0=.false.
  cov0=.false.
  IF(buildpath)THEN
     CALL fidinam(eledir,nam0,'eq0',elefi0,le)
  ELSE
     file=eledir//dircha//nam0//'.eq0'
     CALL rmsp(file,le)
     elefi0(1:le)=file(1:le)
  ENDIF
  INQUIRE(file=elefi0(1:le),exist=ini0)
  IF(ini0)THEN
     CALL read_elems(elk,name,eof,nlsloc,err,lsloc,UNC=unck,FILE=elefi0(1:le))
     IF(dyn%nmod.GT.0)THEN
        nls=nlsloc
        ls=lsloc
     END IF
     ini0=(.not.eof).and.(name.eq.nam0)
     IF(ini0)THEN
        cov0=unck%succ
        WRITE(iunout,*)nam0,' OK'
     ELSE
        WRITE(*,*)elefi0(1:le),'not read properly'
         CALL write_err(nam0,iunout,' ELEMENTS NOT FOUND')
     ENDIF
  ELSE
     WRITE(*,*)elefi0(1:le),'not found'
     CALL write_err(nam0,iunout,' ELEMENT FILE NOT FOUND')
  ENDIF
  nd=6+nls
  IF(coox.eq.' ')RETURN
  IF(ini0)THEN
! coordinate change to equinoctal
     IF(cov0)THEN
        CALL coo_cha(elk,coox,el0,fail_flag,dee)
        CALL convertunc(unck,dee,unc0)
        cov0=unc0%succ
     ELSE
        CALL coo_cha(elk,coox,el0,fail_flag)
     ENDIF
     IF(fail_flag.gt.5)THEN
        CALL write_err(nam0,iunout,' FAILED COORD CHANGE ')
        WRITE(ierrou,*) elk%coo,fail_flag,coox
        WRITE(iunout,*) elk%coo,fail_flag,coox
        ini0=.false.
        cov0=.false.
     ELSEIF(fail_flag.eq.5.and.(coox.eq.'EQU'.or.coox.eq.'KEP'))THEN
        CALL write_err(nam0,iunout,' HYPERBOLIC COORD CHANGE ')
        WRITE(ierrou,*) elk%coo,fail_flag,coox
        WRITE(iunout,*) elk%coo,fail_flag,coox
        ini0=.false.
        cov0=.false.
     ENDIF
!  no covariances available
  ENDIF
END SUBROUTINE find_elems

! Copyright (C) 1997-1998 by Mario Carpino (carpino@brera.mi.astro.it)
! Modified by Andrea Milani, vers. 1.8.3, January 1999
! additional changes for scaling, vers.3.3.1, August 2005
! ---------------------------------------------------------------------
!
!  *****************************************************************
!  *                                                               *
!  *                         W R O M L R                           *
!  *                                                               *
!  *  Writes an orbital element record in an orbital element file  *
!  *                     (multi-line format)                       *
!  *                                                               *
!  *****************************************************************
!
! INPUT:   UNIT      -  Output FORTRAN unit
!           NAME0      -  Name of planet/asteroid/comet
!           ELEM(6)   -  Orbital element vector
!           ELTYPE    -  Type of orbital elements (KEP/EQU/CAR)
!           T0        -  Epoch of orbital elements (MJD, TDT)
!           COVE      -  Covariance matrix of orbital elements
!           DEFCOV    -  Tells whether the covariance matrix is defined
!           NORE      -  Normal matrix of orbital elements
!           DEFNOR    -  Tells whether the normal matrix is defined
!           H         -  H absolute magnitude (if <-100, missing)
!           G         -  G slope parameter
!           MASS      -  Mass (solar masses)
!
! WARNING: the routine does not write the header of the file: this
!          must be generated by calling subroutine wromlh
!
SUBROUTINE wromlr(unit,name0,elem,eltype,t0,cove,defcov,nore,defnor,h,g,mass,n_par,dp)
  USE fund_const
  USE output_control
  USE name_rules
  USE io_elems, ONLY: obscod
  USE dyn_param
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: unit
  DOUBLE PRECISION, INTENT(IN) :: elem(6),t0,h,g,cove(n_par,n_par),nore(n_par,n_par),mass
  CHARACTER*(*), INTENT(IN):: name0,eltype
  LOGICAL, INTENT(IN) ::  defcov,defnor
  INTEGER, INTENT(IN) :: n_par
  DOUBLE PRECISION, INTENT(IN), OPTIONAL :: dp(ndyx)
! ------------END INTERFACE -------------------------------
  DOUBLE PRECISION gammas(ndimx,ndimx),scales(ndimx),units(ndimx) ! for scaling
  CHARACTER*(idnamvir_len) name
  INCLUDE 'parcmc.h90'
  INTEGER l1,ln,i,k,iwd,j
  DOUBLE PRECISION cnv(ndimx),std(ndimx), ele(6), princ
  INTEGER lench
  EXTERNAL lench
! eigenvalues, eigenvectors
  DOUBLE PRECISION eigvec(ndimx,ndimx),eigval(ndimx),fv1(ndimx),fv2(ndimx),wdir(ndimx),sdir
  INTEGER ierr
! Name
  name=' '
  name=name0
  CALL rmsp(name,ln)
  WRITE(unit,100) name(1:ln)
100 FORMAT(A)
! Orbital elements
  cnv=1.d0
  ele=elem
  IF(eltype.EQ.'KEP') THEN
     cnv(3:6)=degrad
     ele(6)=princ(ele(6))
!     if(ele(6).lt.0.d0)ele(6)=ele(6)+dpig
     if(ele(5).lt.0.d0)ele(5)=ele(5)+dpig
     if(ele(4).lt.0.d0)ele(4)=ele(4)+dpig
     WRITE(unit,201) comcha
  201 FORMAT(A,' Keplerian elements: a, e, i, long. node,', &
           &         ' arg. peric., mean anomaly')
     WRITE(unit,101) (ele(i)*cnv(i),i=1,6)
101  FORMAT(' KEP ',1P,E24.16,0P,F18.15,4F18.13)
  ELSEIF(eltype.EQ.'CAR') THEN
     WRITE(unit,202) comcha
202  FORMAT(A,' Cartesian position and velocity vectors')
     WRITE(unit,102) ele
102  FORMAT(' CAR ',1P,6E22.14)
  ELSEIF(eltype.EQ.'EQU') THEN
     cnv(6)=degrad
     ele(6)=princ(ele(6))
     WRITE(unit,203) comcha
  203 FORMAT(A,' Equinoctial elements: a, e*sin(LP), e*cos(LP),',       &
     &         ' tan(i/2)*sin(LN), tan(i/2)*cos(LN), mean long.')
     WRITE(unit,103) (ele(i)*cnv(i),i=1,6)
  103 FORMAT(' EQU ',1P,E22.14,0P,2(1x,f19.15),2(1x,f20.15),1x,F17.13)
  ELSEIF(eltype.EQ.'COM') THEN
     cnv(3:5)=degrad
     if(ele(5).lt.0.d0)ele(5)=ele(5)+dpig
     if(ele(4).lt.0.d0)ele(4)=ele(4)+dpig
     WRITE(unit,204) comcha
  204 FORMAT(A,' Cometary elements: q, e, i, long. node,', &
           &         ' arg. peric., pericenter time')
     IF(abs(ele(6)).lt.1.d6)THEN
        WRITE(unit,114) (ele(i)*cnv(i),i=1,6)
114     FORMAT(' COM ',1P,E22.14,0P,1x,F18.15,1x,3(F18.13,1x),F18.10)
     ELSE
        WRITE(unit,215) (ele(i)*cnv(i),i=1,6)
215     FORMAT(' COM ',1P,E22.14,0P,1x,F18.15,1x,3(F18.13,1x),1P,D18.11)
     ENDIF
  ELSEIF(eltype.EQ.'COT') THEN
     cnv(3:6)=degrad
     if(ele(5).lt.0.d0)ele(5)=ele(5)+dpig
     if(ele(4).lt.0.d0)ele(4)=ele(4)+dpig
     if(ele(6).lt.0.d0)ele(6)=ele(6)+dpig
     WRITE(unit,206) comcha
  206 FORMAT(A,' Cometary elements: q, e, i, long. node,', &
           &         ' arg. peric., true anomaly')
     WRITE(unit,117) (ele(i)*cnv(i),i=1,6)
117  FORMAT(' COT ',1P,E22.14,0P,1x,F18.15,1x,4(F18.13,1x))
  ELSEIF(eltype.eq.'ATT') THEN
     cnv(1:4)=degrad
     IF(ele(1).lt.0.d0)ele(1)=ele(1)+dpig
          WRITE(unit,205) comcha
  205 FORMAT(A,' Attributable elements: R.A., DEC, R.A.dot, DECdot, r, rdot,', &
           &         ' auxiliary information on observer needed')
     WRITE(unit,115) (ele(i)*cnv(i),i=1,6)
115  FORMAT(' ATT ',2(F18.13,1x),2(F18.13,1x),2(F18.13,1x))
     WRITE(unit,116) ' STA ',obscod
116  FORMAT(A5,A3)
  ELSE
     WRITE(*,*)  '**** wromlr: unsupported orbital element type ****', eltype
     STOP '**** wromlr: unsupported orbital element type ****'
  END IF

! Epoch
  WRITE(unit,104) t0
104 FORMAT(' MJD ',F19.9,' TDT')
! Mass
  IF(mass.NE.0.d0) WRITE(unit,105) mass
105 FORMAT(' MAS ',1P,E20.12)
! Magnitudes
  IF(h.GT.-100.d0) WRITE(unit,106) h,g
106 FORMAT(' MAG ',2F7.3)

! nongrav parameters
  IF(n_par.gt.6)THEN
     cnv(7:6+dyn%ndp)=1.d0
     IF(PRESENT(dp))THEN
        WRITE(unit,112) (dp(i)*cnv(6+i),i=1,dyn%ndp)
     ELSE
        WRITE(unit,112) (dp(i)*cnv(6+i),i=1,dyn%ndp)
     ENDIF
112  FORMAT(' NGR ',1P,4E22.14)
  ENDIF


! Covariance matrix
  IF(defcov) THEN
     DO i=1,n_par
        std(i)=SQRT(cove(i,i))
     ENDDO
! eigenvalues
     gammas(1:n_par,1:n_par)=cove(1:n_par,1:n_par)
     CALL weak_dir(gammas(1:n_par,1:n_par),wdir(1:n_par),sdir,-1,eltype,ele,units(1:n_par),n_par)
! Eigenvalues have not to be rescaled
!     scales(1:n_par)=1.d0/units(1:n_par)
!     DO i=1,n_par
!        DO j=1,n_par
!           gammas(i,j)=gammas(i,j)*(scales(i)*scales(j))
!        ENDDO
!     ENDDO
     CALL rs(n_par,n_par,gammas(1:n_par,1:n_par),eigval(1:n_par),1,eigvec(1:n_par,1:n_par),fv1,fv2,ierr)
     DO i=1,n_par
        IF(eigval(i).gt.0.d0)THEN
           eigval(i)=sqrt(eigval(i))
        ELSE
           IF(eigval(i).lt.-1.d-10.and.verb_io.gt.9)THEN
              WRITE(*,*)'wromlr: zero/negative eigenvalue', eigval(i),'  for asteroid ',name(1:ln)
           ENDIF
           eigval(i)=-sqrt(-eigval(i))
        ENDIF
     ENDDO
! RMS, eigenvalues and weak direction are so far commented
! format to be changed for ndimx>10
     WRITE(unit,107) comcha,(std(i)*cnv(i),i=1,n_par)
107  FORMAT(A1,' RMS ',1P,10E14.5)
     WRITE(unit,111) comcha,eigval(1:n_par)
111  FORMAT(A1,' EIG',1P,10E14.5)
     WRITE(unit,110) comcha,wdir(1:n_par)
110  FORMAT(A1,' WEA',10F10.5)
! covariance matrix is given uncommented, to be readable
     WRITE(unit,108) ((cove(i,k)*cnv(i)*cnv(k),k=i,n_par),i=1,n_par)
108  FORMAT(' COV ',1P,3E23.15)
  END IF
! Normal matrix
  IF(defnor) THEN
     WRITE(unit,109) ((nore(i,k)/(cnv(i)*cnv(k)),k=i,n_par),i=1,n_par)
109  FORMAT(' NOR ',1P,3E23.15)
  END IF
END SUBROUTINE wromlr

! ---------------------------------------------------------------------
!
!  *****************************************************************
!  *                                                               *
!  *                          R D _ O R B                            *
!  *                                                               *
!  *     Read orbital elements from a file opened with OPORBF      *
!  *                                                               *
!  *****************************************************************
!
! The routine operates in a sequential way, returning each time
! the orbital elements of the next object contained in the file
!
! OUTPUT:   NAME      -  Name of planet/asteroid/comet
!           ELEM(6)   -  Orbital element vector
!           ELTYPE    -  Type of orbital elements (KEP/EQU/CAR/COM/ATT)
!           T0        -  Epoch of orbital elements (MJD, TDT)
!           COVE      -  Covariance matrix of orbital elements
!           DEFCOV    -  Tells whether the covariance matrix is defined
!           NORE      -  Normal matrix of orbital elements
!           DEFNOR    -  Tells whether the normal matrix is defined
!           H         -  H absolute magnitude (if <-100, missing)
!           G         -  G slope parameter
!           MASS      -  Mass (solar masses)
!           RSYS      -  Reference system type (EQUM/EQUT/ECLM)
!           EPOCH     -  Epoch specification (J2000/OFDATE)
!           KR        -  Record number at which object is found
!           EOF       -  End-of-file flag
!
 SUBROUTINE rd_orb(name,elem,eltype,t0,cove,defcov,nore,defnor,     &
     &                 h,g,mass,rsys,epoch,kr,eof,dp,nd)
   USE fund_const
   USE io_elems
   IMPLICIT NONE
   DOUBLE PRECISION elem(6),t0,h,g,cove(ndimx,ndimx),nore(ndimx,ndimx),mass
   DOUBLE PRECISION,DIMENSION(ndyx),OPTIONAL :: dp ! gravitational parameters
   INTEGER,OPTIONAL :: nd
   CHARACTER*(*) name,eltype,rsys,epoch
   LOGICAL defcov,defnor,eof
   INTEGER kr
   CHARACTER rec*200,rest*200,scale*3
   INTEGER lf,nit,i,k,mjd,mjde,ik
   DOUBLE PRECISION sec,sece,tmp(ndimx*(ndimx+1)/2),cnv(ndimx)
   LOGICAL error,end1,noep
   INTEGER lench,nitchs
! number of elements in upper part of a symmetric matrix
! number of rows to be read
   INTEGER :: nel,nrows
! number of solve for parameters
   INTEGER :: nd_cur
! number of solve for parameters
   EXTERNAL lench,nitchs
   IF(iicorb.NE.36) STOP '**** rd_orb: internal error (01) ****'
   IF(orbunt.LE.0) STOP '**** rd_orb: internal error (02) ****'
! if present dp, then nd is required
  IF(PRESENT(dp).and..not.PRESENT(nd))THEN
     WRITE(*,*)'rd_orb: nd not present'
     STOP
  ENDIF
! number of rows for normal/covariance matrix
   IF(PRESENT(nd))THEN
      nd_cur=nd
   ELSE
      nd_cur=6
   ENDIF
   nel=nd_cur*(nd_cur+1)/2
   nrows=nel/3
   IF(MOD(nel,3).NE.0)THEN
      nrows=nrows+1
   ENDIF
   rsys=dfrsty
   epoch=dfrsep
   mass=0.d0
   DO 1 i=1,ndimx
      DO 2 k=1,ndimx
         cove(i,k)=0
         nore(i,k)=0
2     END DO
1  END DO
   defcov=.false.
   defnor=.false.
   h=-1.d9
   g=-1.d9
   IF(rectyp.EQ.'1L') THEN
      CALL getrsc(orbunt,rec,orbnr,eof)
      IF(eof) RETURN
      CALL strcnt(rec,name,rest,error)
      IF(error) GOTO 20
      kr=orbnr
      nit=nitchs(rest)
      IF(deft0) THEN
         t0=dept0
         IF(nit.LT.8) THEN
            READ(rest,*,ERR=20) elem
         ELSE
            READ(rest,*,ERR=20) elem,h,g
         END IF
      ELSE
         IF(nit.LT.9) THEN
            READ(rest,*,ERR=20) t0,elem
         ELSE
            READ(rest,*,ERR=20) t0,elem,h,g
         END IF
      END IF
      eltype=deltyp
   ELSEIF(rectyp.EQ.'ML') THEN
      IF(nxtend) THEN
         eof=.true.
         RETURN
      END IF
      noep=.true. ! no epoch
! Name
      CALL getrsc(orbunt,name,orbnr,eof)
      IF(eof) RETURN
      kr=orbnr
! Orbital elements (mandatory, immediately after the name)
      CALL getrsc(orbunt,rec,orbnr,end1)
      IF(end1) GOTO 20
      IF(rec(1:1).ne.' ') GOTO 20
      eltype=rec(2:4)
      READ(rec(5:),*,ERR=20) elem
! Other keywords
3     CONTINUE
      CALL getrsc(orbunt,rec,orbnr,end1)
      IF(end1) THEN
         nxtend=.true.
         GOTO 4
      END IF
      IF(rec(1:1).NE.' ') THEN
         BACKSPACE(orbunt)
         orbnr=orbnr-1
         GOTO 4
      END IF
! Epoch of elements
      IF(rec(1:4).EQ.' MJD' .OR. rec(1:4).EQ.' JD ' .OR.            &
           &       rec(1:4).EQ.' CAL') THEN
         CALL ch2tim(rec,mjd,sec,scale,error)
         IF(error) GOTO 20
         CALL cnvtim(mjd,sec,scale,mjde,sece,'TDT')
         t0=mjde+sece/86400.d0
         noep=.false.
      ELSEIF(rec(1:4).EQ.' MAG') THEN
         READ(rec(5:),*,ERR=20) h,g
      ELSEIF(rec(1:4).EQ.' STA')THEN ! station code for ATT only
         READ(rec(5:),*,ERR=20) obscod
      ELSEIF(rec(1:4).EQ.' MAS') THEN
         READ(rec(5:),*,ERR=20) mass
      ELSEIF(rec(1:4).EQ.' NGR') THEN
         IF(PRESENT(dp))THEN
            READ(rec(5:),*,ERR=20) dp(1:dyn%ndp)
         ELSE
            WRITE(*,*)'rd_orb: dp not present, cannot read non gravitational parameters'
            STOP
         ENDIF
      ELSEIF(rec(1:4).EQ.' COV') THEN
         READ(rec(5:),*,ERR=40) (tmp(i),i=1,3)
         DO 17 k=1,nrows-1
            CALL getrsc(orbunt,rec,orbnr,end1)
            IF(end1) GOTO 40
            IF(rec(1:4).NE.' COV') GOTO 40
! the last row can be different if the number of elements if MOD(nel,3) NE 0
            IF(k.EQ.nrows-1.AND.MOD(nel,3).NE.0)THEN
               READ(rec(5:),*,ERR=20) (tmp(i),i=3*k+1,3*k+MOD(nel,3))
            ELSE
               READ(rec(5:),*,ERR=20) (tmp(i),i=3*k+1,3*k+3)
            ENDIF
17       CONTINUE
         ik=0
         DO 8 i=1,nd_cur
            DO 7 k=i,nd_cur
               ik=ik+1
               cove(i,k)=tmp(ik)
7           ENDDO
8        ENDDO
         IF(ik.NE.nel) STOP '**** rd_orb: internal error (03) ****'
         defcov=.true.
      ELSEIF(rec(1:4).EQ.' NOR') THEN
         READ(rec(5:),*,ERR=40) (tmp(i),i=1,3)
         DO 27 k=1,nrows-1
            CALL getrsc(orbunt,rec,orbnr,end1)
            IF(end1) GOTO 40
            IF(rec(1:4).NE.' NOR') GOTO 40
            IF(k.EQ.nrows-1.AND.MOD(nel,3).NE.0)THEN
               READ(rec(5:),*,ERR=20) (tmp(i),i=3*k+1,3*k+MOD(nel,3))
            ELSE
               READ(rec(5:),*,ERR=40) (tmp(i),i=3*k+1,3*k+3)
            ENDIF
   27    ENDDO
         ik=0
         DO 38 i=1,nd_cur
            DO 37 k=i,nd_cur
               ik=ik+1
               nore(i,k)=tmp(ik)
37          ENDDO
38       ENDDO
         IF(ik.NE.nel) STOP '**** rd_orb: internal error (04) ****'
         defnor=.true.
      ELSE
         GOTO 40
      END IF
      GOTO 3

40    CONTINUE
      lf=lench(orbfn)
      WRITE(*,240) orbfn(1:lf),orbnr
240   FORMAT(' ERROR in covariance, file ',A,' at line',I6)

4     CONTINUE
      IF(noep) THEN
         IF(deft0) THEN
            t0=dept0
         ELSE
            GOTO 20
         END IF
      END IF
   ELSE
      STOP '**** rd_orb: internal error (05) ****'
   END IF
! Transformation of angles in orbital elements
! and covariance matrix
   cnv(1:nd_cur)=1
   IF(eltype.EQ.'KEP') THEN
      cnv(3:6)=radeg
   ELSEIF(eltype.EQ.'COM') THEN
      cnv(3:5)=radeg
   ELSEIF(eltype.EQ.'COT')THEN
      cnv(3:6)=radeg
   ELSEIF(eltype.EQ.'EQU') THEN
      cnv(6)=radeg
   ELSEIF(eltype.EQ.'CAR') THEN
      CONTINUE
   ELSEIF(eltype.EQ.'ATT') THEN
      cnv(1:4)=radeg
   ELSE
      STOP '**** rd_orb: internal error (06) ****'
   END IF
   DO  i=1,6
      elem(i)=elem(i)*cnv(i)
   END DO
   IF(defcov) THEN
      DO 13 i=1,nd_cur
         cove(i,i)=cove(i,i)*(cnv(i)**2)
         DO 12 k=i+1,nd_cur
            cove(i,k)=cove(i,k)*cnv(i)*cnv(k)
            cove(k,i)=cove(i,k)
12       ENDDO
13    ENDDO
   END IF
   IF(defnor) THEN
      DO 15 i=1,nd_cur
         nore(i,i)=nore(i,i)/(cnv(i)**2)
         DO 14 k=i+1,nd_cur
            nore(i,k)=nore(i,k)/(cnv(i)*cnv(k))
            nore(k,i)=nore(i,k)
14       ENDDO
15    ENDDO
   END IF
   eof=.false.
   RETURN
20 CONTINUE
   lf=lench(orbfn)
   WRITE(*,200) orbfn(1:lf),orbnr
200 FORMAT(' ERROR in file ',A,' at line',I6)
   STOP '**** rd_orb: abnormal end ****'
END SUBROUTINE rd_orb

! upgrade to a consistent convention for handling non-gravitational
! perturbations
! ==================================================================
! READ_ELEMS
! New version of read_elemns to handle non-gravitational parameters
! Reads elements and dynamical parameters, if present, by using
! only information from the file. The size of the covar/normal
! matrices is determined by the LSP record. The dynamical parameters dp are
! output if given in the file.
! ==================================================================
SUBROUTINE read_elems(el,name,eof,nlsloc,err,lsloc,form,unc,file,unit,mass,rsys,epoch)
  USE station_coordinates, ONLY: statcode
  USE io_elems
  USE planet_masses, ONLY: icrel
! ====BEGIN INTERFACE ==================================
  TYPE(orbit_elem),           INTENT(OUT)   :: el       ! elements
  CHARACTER*(*),              INTENT(OUT)   :: name     ! object name (read from input file)
  LOGICAL,                    INTENT(OUT)   :: eof      ! end of file
  INTEGER,                    INTENT(OUT)   :: nlsloc   ! number of effective solve-for par
  LOGICAL,                    INTENT(OUT)   :: err      ! error flag
  ! optional
  INTEGER,                    INTENT(OUT)   :: lsloc(ndyx) ! solve-for dynamical parameters
  CHARACTER*(*),    OPTIONAL, INTENT(INOUT) :: form        ! can be either '1L' or 'ML' (single/multi line)
  TYPE(orb_uncert), OPTIONAL, INTENT(OUT)   :: unc         ! covariance/normal matrices
  CHARACTER*(*),    OPTIONAL, INTENT(IN)    :: file        ! input file (not used if PRESENT(unit))
  INTEGER,          OPTIONAL, INTENT(IN)    :: unit        ! input unit, only if already opened
  DOUBLE PRECISION, OPTIONAL, INTENT(OUT)   :: mass        ! mass of the object
  CHARACTER*10,     OPTIONAL, INTENT(OUT)   :: rsys,epoch
  !  INTEGER, INTENT(OUT), OPTIONAL :: incfit ! result from incomplete fit;
! WARNING: if file already opened, it must also have the header already written
! ====END INTERFACE ==================================
  INTEGER          :: nd                                      ! dimension of covar/normal matrices
  INTEGER          :: inunit                                  ! input orbit file (unit)
  INTEGER          :: inunitngr                               ! input 1-line catalog with non-grav param
  INTEGER          :: le, i, k, j, ik, lench, lf              ! integers: loop indices and length of files
  CHARACTER*2      :: formcur
  CHARACTER*50     :: namloc
  CHARACTER*6      :: vers                                    ! version (we need to switch to the new version 2.0)
  LOGICAL          :: error, eop, noep, end1                  ! various case of failure
  CHARACTER        :: rec*200,rest*200,scale*3                ! record, full and left
  INTEGER          :: nit,nitchs
  INTEGER          :: kr
  INTEGER          :: mjd,mjde                                ! time (mjd)
  DOUBLE PRECISION :: sec,sece, tmp(ndimx*(ndimx+1)/2)
  DOUBLE PRECISION :: cnv(ndimx),std(ndimx), princ, mass_loc
  DOUBLE PRECISION :: cove(ndimx,ndimx),nore(ndimx,ndimx)     ! covariance and normal matrices
  TYPE(dyn_par)    :: dyn_loc
  INTEGER          :: ndloc
  INTEGER          :: sum_active                              ! sum of the active non-grav parameters
! number of elements in upper part of a symmetric matrix
! number of rows to be read
  INTEGER :: nel,nrows
  INCLUDE 'parcmc.h90'
!--------------------------------------------------------
! Initizialing
  err=.false.;eof=.false.
  dyn_loc = undefined_dyn_par
!  IF(PRESENT(lsloc))THEN
  lsloc=0
!  ENDIF
! to avoid unassigned variable
  IF(PRESENT(mass))THEN
     mass=0.d0
  ENDIF
! open output file if not opened already
  IF(.not.PRESENT(unit))THEN
     CALL oporbf(file,0)
     vers=oef_vers
     IF(PRESENT(form))THEN
        form=rectyp
     ENDIF
     inunit=orbunt
     IF(orbunt.LE.0)THEN
        WRITE(*,*) ' read_elems: failed oporbf, orbunt= ',orbunt
        STOP ' read_elems: failed oporbf'
     ENDIF
  ELSE
     inunit=unit
     vers=oef_vers
  ENDIF
  IF(PRESENT(rsys).AND.PRESENT(epoch))THEN
     rsys=dfrsty
     epoch=dfrsep
  END IF
  IF(vers.EQ.'OEF1.1')THEN
     nd=6
     nlsloc=nd-6
!     IF(PRESENT(lsloc))THEN
        lsloc=0
!     ENDIF
  END IF
! which record type?
  IF(rectyp.eq.'1L')THEN
     el=undefined_orbit_elem
     CALL getrsc(orbunt,rec,orbnr,eof)
     IF(eof) RETURN
     CALL strcnt(rec,name,rest,error)
     IF(error) GOTO 20
     kr=orbnr
     nit=nitchs(rest)
     IF(nit.LT.8) THEN
! no magnitude available
        READ(rest,*,ERR=20) el%t,el%coord
        el%mag_set=.false.
     ELSEIF(nit.EQ.9)THEN
! include magnitude H,G
        READ(rest,*,ERR=20) el%t,el%coord,el%h_mag,el%g_mag
        el%mag_set=.true.
! old version, no dyn-param
     ELSEIF(nit.EQ.10.AND.vers.EQ.'OEF1.1')THEN
        ! include magnitude H,G
        ! No dyn parameters
        READ(rest,*,ERR=20) el%t,el%coord,el%h_mag,el%g_mag
        el%mag_set=.TRUE.
! new version: non-parameters
     ELSEIF(nit.EQ.10.AND.vers.EQ.'OEF2.0')THEN
        READ(rest,*,ERR=20) el%t,el%coord,el%h_mag,el%g_mag,sum_active
        el%mag_set=.TRUE.
        nd=6+sum_active
     ELSE
        WRITE(*,*)'read_elems: one magnitude missing, nit=',nit
        READ(rest,*,ERR=20) el%t,el%coord,el%h_mag
        el%mag_set=.true.
     END IF
     el%coo=deltyp
  ELSEIF(rectyp.eq.'ML')THEN
      IF(nxtend)THEN
        eof=.true.
        RETURN
     ENDIF
! setup of number of parameters in the default way, without nongrav
     el=undefined_orbit_elem
! initialize unc
     IF(PRESENT(unc))CALL undefined_orb_uncert(ndimx,unc)
     noep=.true. ! no epoch
! Name
     CALL getrsc(orbunt,name,orbnr,eof)
     IF(eof) GOTO 4
     kr=orbnr
! Orbital elements (mandatory, immediately after the name)
     CALL getrsc(orbunt,rec,orbnr,end1)
     IF(end1) GOTO 20
     IF(rec(1:1).ne.' ') GOTO 20
     el%coo=rec(2:4)
     READ(rec(5:),*,ERR=20) el%coord
     ! Other keywords
     DO j=1,100
        CALL getrsc(orbunt,rec,orbnr,end1)
        IF(end1) THEN
           nxtend=.true.
           IF(noep) GOTO 20 ! epoch time missing
           GOTO 4 !end of file, but normalize the last
        END IF
        IF(rec(1:1).NE.' ') THEN
! next object, needs to be unread
           BACKSPACE(orbunt)
           orbnr=orbnr-1
           IF(noep) GOTO 20 ! epoch time missing
           GOTO 4  !end of one object, normalize
        END IF
! Epoch of elements
        IF(rec(1:4).EQ.' MJD' .OR. rec(1:4).EQ.' JD ' .OR.            &
             &       rec(1:4).EQ.' CAL') THEN
           CALL ch2tim(rec,mjd,sec,scale,error)
           IF(error) GOTO 20
           CALL cnvtim(mjd,sec,scale,mjde,sece,'TDT')
           el%t=mjde+sece/86400.d0
           noep=.false. ! epoch time available
        ELSEIF(rec(1:4).EQ.' MAG') THEN
           READ(rec(5:),*,ERR=20) el%h_mag,el%g_mag
           el%mag_set=.true.
        ELSEIF(rec(1:4).EQ.' STA')THEN ! station code for ATT only
           READ(rec(5:),*,ERR=20) obscod
           IF(el%coo.eq.'ATT') THEN
              CALL statcode(obscod, el%obscode)
           ENDIF
        ELSEIF(rec(1:4).EQ.' MAS') THEN
! problem: this mass is not used
           READ(rec(5:),*,ERR=20) mass_loc
           IF(PRESENT(mass))THEN
              mass=mass_loc
           END IF
! nongrav handling
        ELSEIF(rec(1:4).EQ.' LSP')THEN
           nit=nitchs(rec(5:))
           IF(nit-3.GE.0.AND.nit-3.LE.ndyx)THEN
              READ(rec(5:),*,ERR=20)dyn_loc%nmod, dyn_loc%ndp, nd, lsloc(1:nit-3)
              nlsloc=nd-6
            ELSE IF(nit-3.LT.0) THEN
              WRITE(*,*) 'read_elems: error in reading LSP field: ', rec(5:)
              STOP
           ENDIF
        ELSEIF(rec(1:4).EQ.' NGR') THEN
           READ(rec(5:),*,ERR=20) dyn_loc%dp(1:dyn_loc%ndp)
           IF(dyn_loc%dp(1).NE.0.d0)THEN
              dyn_loc%active(1)=.TRUE.
           END IF
           IF(dyn_loc%dp(2).NE.0.d0)THEN
              dyn_loc%active(2)=.TRUE.
              IF(dyn_loc%nmod.EQ.1)THEN
                 icrel=2
              END IF
           END IF
           IF(dyn_loc%dp(3).NE.0.d0)THEN
              dyn_loc%active(3)=.TRUE.
           END IF
           IF(dyn_loc%dp(4).NE.0.d0)THEN
              dyn_loc%active(4)=.TRUE.
           END IF
        ELSEIF(rec(1:9).EQ.' APR YAR ') THEN
           dyn_loc%apriori(2)=.TRUE.
           READ(rec(10:),*,ERR=20) dyn_loc%apriori_nominal(2), dyn_loc%apriori_std(2)
        ELSEIF(rec(1:9).EQ.' APR A/M ') THEN
           dyn_loc%apriori(1)=.TRUE.
           READ(rec(10:),*,ERR=20) dyn_loc%apriori_nominal(1), dyn_loc%apriori_std(1)
        ELSEIF(rec(1:9).EQ.' APR A1C ') THEN
           dyn_loc%apriori(1)=.TRUE.
           READ(rec(10:),*,ERR=20) dyn_loc%apriori_nominal(1), dyn_loc%apriori_std(1)
        ELSEIF(rec(1:9).EQ.' APR A2C ') THEN
           dyn_loc%apriori(2)=.TRUE.
           READ(rec(10:),*,ERR=20) dyn_loc%apriori_nominal(2), dyn_loc%apriori_std(2)
        ELSEIF(rec(1:9).EQ.' APR A3C ') THEN
           dyn_loc%apriori(3)=.TRUE.
           READ(rec(10:),*,ERR=20) dyn_loc%apriori_nominal(3), dyn_loc%apriori_std(3)
        ELSEIF(rec(1:9).EQ.' APR DTD ') THEN
           dyn_loc%apriori(4)=.TRUE.
           READ(rec(10:),*,ERR=20) dyn_loc%apriori_nominal(4), dyn_loc%apriori_std(4)
        ELSEIF(rec(1:4).EQ.' COV') THEN
! covariance and normal matrices: size
           nel=nd*(nd+1)/2
           nrows=nel/3
! initialize cove and nore
           DO  i=1,ndimx
              DO  k=1,ndimx
                 cove(i,k)=0
                 nore(i,k)=0
              END DO
           END DO
! initialize unc
           IF(PRESENT(unc))CALL undefined_orb_uncert(nd,unc)
           IF(MOD(nel,3).NE.0)THEN
              nrows=nrows+1
           ENDIF
           READ(rec(5:),*,ERR=40) (tmp(i),i=1,3)
           DO  k=1,nrows-1
              CALL getrsc(orbunt,rec,orbnr,end1)
              IF(end1) GOTO 40
              IF(rec(1:4).NE.' COV') GOTO 40
! the last row can be different if the number of elements if MOD(nel,3) NE 0
              IF(k.EQ.nrows-1.AND.MOD(nel,3).NE.0)THEN
                 READ(rec(5:),*,ERR=20) (tmp(i),i=3*k+1,3*k+MOD(nel,3))
              ELSE
                 READ(rec(5:),*,ERR=20) (tmp(i),i=3*k+1,3*k+3)
              ENDIF
           ENDDO
           ik=0
           DO i=1,nd
              DO  k=i,nd
                 ik=ik+1
                 cove(i,k)=tmp(ik)
              ENDDO
           ENDDO
           IF(ik.NE.nel) STOP 'read_elems: insufficient COV records'
        ELSEIF(rec(1:4).EQ.' NOR') THEN
! assumed NOR is after COV, thus size has been computed
           READ(rec(5:),*,ERR=40) (tmp(i),i=1,3)
           DO k=1,nrows-1
              CALL getrsc(orbunt,rec,orbnr,end1)
              IF(end1) GOTO 40
              IF(rec(1:4).NE.' NOR') GOTO 40
              IF(k.EQ.nrows-1.AND.MOD(nel,3).NE.0)THEN
                 READ(rec(5:),*,ERR=20) (tmp(i),i=3*k+1,3*k+MOD(nel,3))
              ELSE
                 READ(rec(5:),*,ERR=40) (tmp(i),i=3*k+1,3*k+3)
              ENDIF
           ENDDO
           ik=0
           DO i=1,nd
              DO k=i,nd
                 ik=ik+1
                 nore(i,k)=tmp(ik)
              ENDDO
           ENDDO
           IF(ik.NE.nel) STOP 'read_elems: insufficient NOR records'
! covariance is available, but copied lateafter rescaling
           IF(PRESENT(unc)) unc%succ=.true.
        ELSE
           WRITE(*,*) 'read_elems: unrecognized keyword =',rec(1:4)
           STOP ' unrecognized record'
        END IF
! potentially infinite loop removed
     ENDDO

40   CONTINUE
     lf=lench(orbfn)
     WRITE(*,240) orbfn(1:lf),orbnr
240  FORMAT(' ERROR in covariance, file ',A,' at line',I6)
  ELSE
     WRITE(*,*) ' read_elems: unknown format, form=',form
     STOP ' orbit format?? '
  ENDIF
4 CONTINUE
! elements and possibly matrices read, now rescale
! Transformation of angles in orbital elements
! and covariance matrix
  cnv(1:nd)=1
   IF(el%coo.EQ.'KEP') THEN
      cnv(3:6)=radeg
   ELSEIF(el%coo.EQ.'COM') THEN
      cnv(3:5)=radeg
   ELSEIF(el%coo.EQ.'COT')THEN
      cnv(3:6)=radeg
   ELSEIF(el%coo.EQ.'EQU') THEN
      cnv(6)=radeg
   ELSEIF(el%coo.EQ.'CAR') THEN
      CONTINUE
   ELSEIF(el%coo.EQ.'ATT') THEN
      cnv(1:4)=radeg
   ELSE
      WRITE(*,*)'read_elems: unknown coordinate type=',el%coo
      STOP 'unknown coordinate type'
   ENDIF
   DO  i=1,6
      el%coord(i)=el%coord(i)*cnv(i)
   END DO
   IF(PRESENT(unc))THEN
      IF(unc%succ)THEN
         DO i=1,nd
            cove(i,i)=cove(i,i)*(cnv(i)**2)
            DO k=i+1,nd
               cove(i,k)=cove(i,k)*cnv(i)*cnv(k)
               cove(k,i)=cove(i,k)
            ENDDO
         ENDDO
         DO i=1,nd
            nore(i,i)=nore(i,i)/(cnv(i)**2)
            DO k=i+1,nd
               nore(i,k)=nore(i,k)/(cnv(i)*cnv(k))
               nore(k,i)=nore(i,k)
            ENDDO
         ENDDO
! copy into output, if output array available
         unc%g(1:nd,1:nd)=cove(1:nd,1:nd)
         unc%c(1:nd,1:nd)=nore(1:nd,1:nd)
      ENDIF
   ENDIF
   IF(.NOT.ngr_opt)THEN
      dyn=dyn_loc
   ELSE
   END IF
   IF(.not.PRESENT(unit)) CALL clorbf
! succesful termination
   RETURN
! error condition, reading not succeded
20 CONTINUE
   WRITE(*,*)'read_elems: error condition '
   err=.true.
 END SUBROUTINE read_elems


! ==================================================================
! WRITE_ELEMS
! Writes elements and dynamical parameters, if present, by using
! only information from dummy arguments. The size of the covar/normal
! matrices is determined by unc%ndim.
! A record contains the list ls, the integer nls.
! ==================================================================
SUBROUTINE write_elems(el,name,form,dynam_par,nlsloc,lsloc,unc,incfit,file,unit,unit2)
  USE station_coordinates, ONLY: codestat
  USE io_elems, ONLY: obscod
! ====BEGIN INTERFACE ====================================================================================
  TYPE(orbit_elem),          INTENT(IN) :: el            ! orbital elements
  CHARACTER*(*),             INTENT(IN) :: name          ! name
  CHARACTER*(*),   OPTIONAL, INTENT(IN) :: form          ! can be either '1L' or 'ML',
  TYPE(dyn_par),   OPTIONAL, INTENT(IN) :: dynam_par     ! dynamical non-grav parameters
  INTEGER,         OPTIONAL, INTENT(IN) :: nlsloc        ! number of effective solve-for par
  INTEGER,         OPTIONAL, INTENT(IN) :: lsloc(ndyx)   ! solve-for dynamic parameters
  TYPE(orb_uncert),OPTIONAL, INTENT(IN) :: unc           ! covariance/normal matrices
  INTEGER,         OPTIONAL, INTENT(IN) :: incfit        ! result from incomplete fit;
  CHARACTER*(*),   OPTIONAL, INTENT(IN) :: file          ! output file
                                                         ! (not used if PRESENT(unit))
  INTEGER,         OPTIONAL, INTENT(IN) :: unit          ! input unit, only if already opened
                                                         ! WARNING: if file already opened,
                                                         ! it must also have the header already written
  INTEGER,         OPTIONAL, INTENT(IN) :: unit2         ! input unit, only if already opened for 1L catalog
                                                         ! in case of non-grav
! ====END INTERFACE =======================================================================================
  INTEGER           :: nd                                         ! dimension of covar/normal matrices
  INTEGER           :: uniout, uniout_ngr                         ! output orbit file and non-grav file
  CHARACTER*10      :: rsys,epoch
  INTEGER           :: le, i, k                                   ! loop indices
! Local version of the optional variables
  CHARACTER*2       :: formcur
  CHARACTER*50      :: namloc
  TYPE(dyn_par)     :: dyn_loc
  INTEGER           :: sum_active                                 ! sum of the active non-grav parameters
  CHARACTER*6       :: vers                                       ! version
  DOUBLE PRECISION  :: cnv(ndimx),std(ndimx), ele(6), princ, mass
  DOUBLE PRECISION  :: units(ndimx)                               ! for scaling
! eigenvalues, eigenvectors
  DOUBLE PRECISION  :: eigvec(ndimx,ndimx),eigval(ndimx),fv1(ndimx),fv2(ndimx),wdir(ndimx),sdir
  INTEGER           :: ierr
  INCLUDE 'parcmc.h90'
! ==========================================================================================================
! initialization
  sum_active=0
! form
  IF(PRESENT(form))THEN
     formcur=form
  ELSE
     formcur='ML'
  ENDIF
! version
  vers='OEF2.0'
! dynamical parameters
  IF(PRESENT(dynam_par))THEN
     dyn_loc=dynam_par
  ELSE
     dyn_loc=undefined_dyn_par
  END IF
! dyn%nmod must be ge 0
  IF(dyn_loc%nmod.lt.0)THEN
     WRITE(*,*)'write_elems: wrong nongrav model nmod=', dyn_loc%nmod
     STOP
  END IF
! ndyx is the maximum for ndp
  IF(dyn_loc%ndp.gt.ndyx)THEN
     WRITE(*,*)'write_elems: wrong dimension given for dp, ndp= ', dyn_loc%ndp
     STOP
  ENDIF
! lsloc and nlsloc must come together, and with meaningful dimensions
  IF(PRESENT(lsloc).AND..NOT.PRESENT(nlsloc))THEN
     WRITE(*,*)'write_elems: lsloc is present without nlsloc ', lsloc
     STOP
  END IF
  IF(PRESENT(lsloc).AND.nlsloc.GT.0)THEN
     IF(nlsloc.GT.ndyx)THEN
        WRITE(*,*)'write_elems: wrong dimension for lsloc, nlsloc= ',nlsloc
        STOP
     ELSEIF(dyn_loc%ndp.LT.nlsloc)THEN
        WRITE(*,*)'write_elems: matrix larger than state vector ',nlsloc, dyn_loc%ndp
        STOP
     ENDIF
  ELSE IF(.NOT.PRESENT(lsloc))THEN
     IF(PRESENT(nlsloc).AND.nlsloc.gt.0)THEN
        WRITE(*,*)'write_elems: dimension given for non existing lsloc, nlsloc=',nlsloc
        STOP
     ENDIF
  ENDIF
  IF(.NOT.PRESENT(unit))THEN
     IF(PRESENT(file))THEN
        CALL filopn(uniout,file,'unknown')
! write file header
        rsys='ECLM'
        epoch='J2000'
        IF(formcur.eq.'1L')THEN
           CALL wro1lh2(uniout,rsys,epoch,el%coo,vers)
        ELSEIF(formcur.eq.'ML')THEN
           CALL wromlh2(uniout,rsys,epoch,vers)
        ENDIF
     ENDIF
  ELSE
     uniout=unit
  ENDIF
  IF(PRESENT(unit2))THEN
     uniout_ngr=unit2
  ELSE
     IF(PRESENT(file).AND.dyn_loc%ndp.GT.0)THEN
        IF(formcur.EQ.'1L')THEN
           CALL filopn(uniout_ngr,file,'unknown')
        ! write file header
           rsys='ECLM'
           epoch='J2000'
           CALL wro1lh2_ngr(uniout_ngr,rsys,epoch,vers,dyn_loc%ndp)
        ENDIF
     END IF
  END IF
! name without blanks
  namloc=' '
  namloc=name
  CALL rmsp(namloc,le)
! solutions with less than 6 parameters
  IF(PRESENT(incfit))THEN
     IF(incfit.eq.5)THEN
        WRITE(uniout,200)namloc
200     FORMAT('! CONSTRAINED SOLUTION for ',a)
     ELSEIF(incfit.eq.4)THEN
        WRITE(uniout,291)name
291     FORMAT('! 4_PARAMETERS SOLUTION for ',a)
     ENDIF
  ENDIF

! covariance/normal matrices
  IF(PRESENT(nlsloc).AND.nlsloc.gt.0)THEN
     nd=6+nlsloc
  ELSE
     nd=6
  ENDIF

  IF(PRESENT(unc))THEN
     IF(unc%ndim.ne.nd)THEN
        WRITE(*,*)' write_elems: ndim.ne.nd ', unc%ndim, nd
        !STOP
     ENDIF
  ENDIF

! when el%coo=ATT, station numeric code is needed
  CALL codestat(el%obscode,obscod)

! write in format 1L/ML
  IF(formcur.eq.'1L')THEN
! Check active non-grav parameters
     IF(dyn_loc%nmod.GT.0)THEN
        DO i=1,nlsloc
           IF(dyn_loc%active(lsloc(i)))THEN
              sum_active=sum_active+1
           END IF
        END DO
     END IF
     IF(sum_active.GT.0)THEN
        CALL wro1lr(uniout,namloc(1:le),el%coord,el%coo,el%t,el%h_mag,el%g_mag,sum_active)
        CALL wro1lr_ngr(uniout_ngr,namloc(1:le),dyn_loc%nmod,dyn_loc%ndp,dyn_loc%dp)
     ELSE
        CALL wro1lr(uniout,namloc(1:le),el%coord,el%coo,el%t,el%h_mag,el%g_mag)
     END IF
  ELSEIF(formcur.eq.'ML')THEN
! object name
     WRITE(uniout,'(A)') namloc(1:le)
! Orbital elements
     cnv=1.d0
     ele=el%coord
     IF(el%coo.EQ.'KEP') THEN
        cnv(3:6)=degrad
        ele(6)=princ(ele(6))
        IF(ele(5).lt.0.d0)ele(5)=ele(5)+dpig
        IF(ele(4).lt.0.d0)ele(4)=ele(4)+dpig
        WRITE(uniout,201) comcha
201     FORMAT(A,' Keplerian elements: a, e, i, long. node, arg. peric., mean anomaly')
        WRITE(uniout,101) (ele(i)*cnv(i),i=1,6)
101     FORMAT(' KEP ',1P,E24.16,0P,F18.15,4F18.13)
     ELSEIF(el%coo.EQ.'CAR') THEN
        WRITE(uniout,202) comcha
202     FORMAT(A,' Cartesian position and velocity vectors')
        WRITE(uniout,102) ele
102     FORMAT(' CAR ',1P,6E22.14)
     ELSEIF(el%coo.EQ.'EQU') THEN
        cnv(6)=degrad
        ele(6)=princ(ele(6))
        WRITE(uniout,203) comcha
203     FORMAT(A,' Equinoctial elements: a, e*sin(LP), e*cos(LP),',       &
             &         ' tan(i/2)*sin(LN), tan(i/2)*cos(LN), mean long.')
        WRITE(uniout,103) (ele(i)*cnv(i),i=1,6)
103     FORMAT(' EQU ',1P,E24.16,0P,2(1x,f19.15),2(1x,f20.15),1x,F17.13)
     ELSEIF(el%coo.EQ.'COM') THEN
        cnv(3:5)=degrad
        if(ele(5).lt.0.d0)ele(5)=ele(5)+dpig
        if(ele(4).lt.0.d0)ele(4)=ele(4)+dpig
        WRITE(uniout,204) comcha
204     FORMAT(A,' Cometary elements: q, e, i, long. node,', &
             &         ' arg. peric., pericenter time')
        IF(abs(ele(6)).lt.1.d6)THEN
           WRITE(uniout,114) (ele(i)*cnv(i),i=1,6)
114        FORMAT(' COM ',1P,E22.14,0P,1x,F18.15,1x,3(F18.13,1x),F18.10)
        ELSE
           WRITE(uniout,215) (ele(i)*cnv(i),i=1,6)
215        FORMAT(' COM ',1P,E22.14,0P,1x,F18.15,1x,3(F18.13,1x),1P,D18.11)
        ENDIF
     ELSEIF(el%coo.EQ.'COT') THEN
        cnv(3:6)=degrad
        if(ele(5).lt.0.d0)ele(5)=ele(5)+dpig
        if(ele(4).lt.0.d0)ele(4)=ele(4)+dpig
        if(ele(6).lt.0.d0)ele(6)=ele(6)+dpig
        WRITE(uniout,206) comcha
206     FORMAT(A,' Cometary elements: q, e, i, long. node,', &
             &         ' arg. peric., true anomaly')
        WRITE(uniout,117) (ele(i)*cnv(i),i=1,6)
117     FORMAT(' COT ',1P,E22.14,0P,1x,F18.15,1x,4(F18.13,1x))
     ELSEIF(el%coo.eq.'ATT') THEN
        cnv(1:4)=degrad
        IF(ele(1).lt.0.d0)ele(1)=ele(1)+dpig
        WRITE(uniout,205) comcha
205     FORMAT(A,' Attributable elements: R.A., DEC, R.A.dot, DECdot, r, rdot,', &
             &         ' auxiliary information on observer needed')
        WRITE(uniout,115) (ele(i)*cnv(i),i=1,6)
115     FORMAT(' ATT ',2(F18.13,1x),2(F18.13,1x),2(F18.13,1x))
        WRITE(uniout,116) ' STA ',obscod
116     FORMAT(A5,A3)
     ELSE
        WRITE(*,*)  '**** write_elems: unsupported orbital element type ****', el%coo
        STOP '**** write_elems: unsupported orbital element type ****'
     END IF
! Epoch
     WRITE(uniout,104) el%t
104  FORMAT(' MJD ',F19.9,' TDT')
! Mass: not ready
!     IF(mass.NE.0.d0) WRITE(uniout,105) mass
!105  FORMAT(' MAS ',1P,E20.12)
! Magnitudes
     !IF(el%mag_set) WRITE(uniout,106)el%h_mag,el%g_mag
     IF(el%h_mag.gt.-100.d0) WRITE(uniout,106)el%h_mag,el%g_mag
106  FORMAT(' MAG ',2F7.3)

     !IF(formcur.EQ.'ML')THEN
! nongrav parameters
     IF(dyn_loc%nmod.ge.0)THEN
! ndp.gt.0 implies PRESENT(nmod) nmod.gt.0
! record giving the matrix dimension and list of solved parameters
        IF(PRESENT(nlsloc).AND.nlsloc.GT.0)THEN
           WRITE(uniout,149) comcha
149        FORMAT(A,' Non-grav parameters: model used, actual number in use, dimension, ', &
                & 'list of solve for parameters to be currently determined')
           WRITE(uniout,150) dyn_loc%nmod, dyn_loc%ndp, nd, lsloc(1:nlsloc)
150        FORMAT(' LSP ',1X,I2,1X,I2,3X,I2,2X,4(1X,I2))
        ELSEIF(PRESENT(incfit).AND.incfit.EQ.5) THEN
           WRITE(uniout,'(A)') '! Non-grav parameters: model used, actual number in use, dimension'
           WRITE(uniout,155) dyn_loc%nmod, dyn_loc%ndp, 5
        ELSE
           WRITE(uniout,'(A)') '! Non-grav parameters: model used, actual number in use, dimension'
           WRITE(uniout,155) dyn_loc%nmod, dyn_loc%ndp, nd
155        FORMAT(' LSP ',1X,I2,1X,I2,3X,I2)
        END IF
        IF(dyn_loc%nmod.gt.0)THEN
           ! because ndp.gt.nlsloc, weights are available for output matrices
           cnv(7:6+dyn_loc%ndp)=1.d0
           WRITE(uniout,'(A)') '! Vector of dynamical parameters'
           WRITE(uniout,112) ((dyn_loc%dp(i))*cnv(6+i),i=1,dyn_loc%ndp)
112        FORMAT(' NGR ',1P,4E22.14)
        END IF
     ENDIF
     !END IF

! A-priori for non-gravitational parameter
! ========================================
! Dynamical model nmod=1
     IF(dyn_loc%nmod.EQ.1)THEN
! A/M and A2 to be determined
        IF(lsloc(1).EQ.1.AND.lsloc(2).EQ.2)THEN
! A/M and A2 a-priori available
           IF(dyn_loc%apriori(2).AND.dyn_loc%apriori(1))THEN
              WRITE(uniout,'(A)') '! A priori value and standard deviation (A/M and Yarkovsky effect)'
              WRITE(uniout,352) dyn_loc%apriori_nominal(1), dyn_loc%apriori_std(1)
              WRITE(uniout,351) dyn_loc%apriori_nominal(2), dyn_loc%apriori_std(2)
! only A2 a-priori available
           ELSE IF(dyn_loc%apriori(2))THEN
              WRITE(uniout,'(A)') '! A priori value and standard deviation (only Yarkovsky effect)'
              WRITE(uniout,351) dyn_loc%apriori_nominal(2), dyn_loc%apriori_std(2)
! only A/M a-priori available
           ELSE IF(dyn_loc%apriori(1))THEN
              WRITE(uniout,'(A)') '! A priori value and standard deviation (only A/M)'
              WRITE(uniout,352) dyn_loc%apriori_nominal(1), dyn_loc%apriori_std(1)
           END IF
! only A2 to be determined with a-priori
        ELSE IF(lsloc(1).EQ.2.AND.dyn_loc%apriori(2))THEN
           WRITE(uniout,'(A)') '! A priori value and standard deviation (only Yarkovsky effect)'
           WRITE(uniout,351) dyn_loc%apriori_nominal(2), dyn_loc%apriori_std(2)
! only A/M to be determined with a-priori
        ELSE IF(lsloc(1).EQ.1.AND.dyn_loc%apriori(1))THEN
           WRITE(uniout,'(A)') '! A priori value and standard deviation (only A/M)'
           WRITE(uniout,352) dyn_loc%apriori_nominal(1), dyn_loc%apriori_std(1)
        END IF

! Dynamical model nmod=2
     ELSEIF(dyn_loc%nmod.EQ.2)THEN
        IF(lsloc(1).EQ.1.AND.lsloc(2).EQ.2.AND.lsloc(3).EQ.3.AND.lsloc(4).EQ.4)THEN
           IF(dyn_loc%apriori(1).AND.dyn_loc%apriori(2).AND.dyn_loc%apriori(3).AND.dyn_loc%apriori(4))THEN
              WRITE(uniout,'(A)') '! A priori value and standard deviation (A1, A2, A3 and dtdelay)'
              WRITE(uniout,353) dyn_loc%apriori_nominal(1), dyn_loc%apriori_std(1)
              WRITE(uniout,354) dyn_loc%apriori_nominal(2), dyn_loc%apriori_std(2)
              WRITE(uniout,355) dyn_loc%apriori_nominal(3), dyn_loc%apriori_std(3)
              WRITE(uniout,356) dyn_loc%apriori_nominal(4), dyn_loc%apriori_std(4)
           END IF
        ELSE IF(lsloc(1).EQ.1.AND.lsloc(2).EQ.2.AND.lsloc(3).EQ.3)THEN
           IF(dyn_loc%apriori(1).AND.dyn_loc%apriori(2).AND.dyn_loc%apriori(3))THEN
              WRITE(uniout,'(A)') '! A priori value and standard deviation (A1, A2, A3)'
              WRITE(uniout,353) dyn_loc%apriori_nominal(1), dyn_loc%apriori_std(1)
              WRITE(uniout,354) dyn_loc%apriori_nominal(2), dyn_loc%apriori_std(2)
              WRITE(uniout,355) dyn_loc%apriori_nominal(3), dyn_loc%apriori_std(3)
           END IF
        ELSE IF(lsloc(1).EQ.1.AND.lsloc(2).EQ.2)THEN
           IF(dyn_loc%apriori(1).AND.dyn_loc%apriori(2))THEN
              WRITE(uniout,'(A)') '! A priori value and standard deviation (A1, A2)'
              WRITE(uniout,353) dyn_loc%apriori_nominal(1), dyn_loc%apriori_std(1)
              WRITE(uniout,354) dyn_loc%apriori_nominal(2), dyn_loc%apriori_std(2)
           END IF
        END IF
     END IF
351  FORMAT(' APR YAR ',1P,2E22.14)
352  FORMAT(' APR A/M ',1P,2E22.14)
353  FORMAT(' APR A1C ',1P,2E22.14)
354  FORMAT(' APR A2C ',1P,2E22.14)
355  FORMAT(' APR A3C ',1P,2E22.14)
356  FORMAT(' APR DTD ',1P,2E22.14)

! covariance matrix
     IF(PRESENT(unc)) THEN
        IF(.not.unc%succ)THEN
           WRITE(*,*)' write_elems: covariance not given'
        ELSE
           DO i=1,nd
              std(i)=SQRT(unc%g(i,i))
           ENDDO
! eigenvalues
           CALL weak_dir(unc%g(1:nd,1:nd),wdir(1:nd),sdir,-1,el%coo,ele,units(1:nd),nd)
           CALL rs(nd,nd,unc%g(1:nd,1:nd),eigval(1:nd),1,eigvec(1:nd,1:nd),fv1,fv2,ierr)
           DO i=1,nd
              IF(eigval(i).gt.0.d0)THEN
                 eigval(i)=sqrt(eigval(i))
              ELSE
                 IF(eigval(i).lt.-1.d-10.and.verb_io.gt.9)THEN
                    WRITE(*,*)'wromlr: zero/negative eigenvalue', eigval(i),'  for asteroid ',namloc(1:le)
                 ENDIF
                 eigval(i)=-sqrt(-eigval(i))
              ENDIF
           ENDDO
! RMS, eigenvalues and weak direction are so far commented
! format to be changed for ndimx>10
           WRITE(uniout,107) comcha,(std(i)*cnv(i),i=1,nd)
107        FORMAT(A1,' RMS ',1P,10E14.5)
           WRITE(uniout,111) comcha,eigval(1:nd)
111        FORMAT(A1,' EIG',1P,10E14.5)
           WRITE(uniout,110) comcha,wdir(1:nd)
110        FORMAT(A1,' WEA',10F10.5)
! covariance matrix is given uncommented, to be readable
           WRITE(uniout,108) ((unc%g(i,k)*cnv(i)*cnv(k),k=i,nd),i=1,nd)
108        FORMAT(' COV ',1P,3E23.15)
           WRITE(uniout,109) ((unc%c(i,k)/(cnv(i)*cnv(k)),k=i,nd),i=1,nd)
109        FORMAT(' NOR ',1P,3E23.15)
        ENDIF
     END IF
  ENDIF
  IF(.not.PRESENT(unit)) CALL filclo(uniout,' ')
END SUBROUTINE write_elems

! Copyright (C) 1997-1998 by Mario Carpino (carpino@brera.mi.astro.it)
! Version: June 21, 1998
! Version December 18, 1998 Steven Chesley (chesley@dm.unipi.it)
! Added ELTYPE, Changed time field to F13.4 to handle jpleph406
! Version May 20, 2015 Federica Spoto (spoto@mail.dm.unipi.it)
! Now the 1L catalogs must contain a flag to indicate if the dynamical
! parameters are active
! ---------------------------------------------------------------------
!
!  *****************************************************************
!  *                                                               *
!  *                         W R O 1 L R                           *
!  *                                                               *
!  *  Writes an orbital element record in an orbital element file  *
!  *   (single-line format, different epochs, keplerian elements)  *
!  *                     (multi-line format)                       *
!  *                                                               *
!  *****************************************************************
!
! OUTPUT:   UNIT      -  Output FORTRAN unit
!           NAME      -  Name of planet/asteroid/comet
!           ELEM(6)   -  Orbital element vector (keplerian elements)
!           ELTYPE    -  Type of orbital elements (KEP/EQU/CAR)
!           T0        -  Epoch of orbital elements (MJD, TDT)
!           H         -  H absolute magnitude (if <-100, missing)
!           G         -  G slope parameter
!
! WARNING: the routine does not write the header of the file: this
!          must be generated by calling subroutine wro1lh
!
SUBROUTINE wro1lr(unit,name,elem,eltype,t0,h,g,ngr_act)
  USE fund_const
  IMPLICIT NONE
! Begin interface
  INTEGER,           INTENT(IN)           :: unit
  DOUBLE PRECISION,  INTENT(IN)           :: elem(6),t0,h,g
  CHARACTER*(*),     INTENT(IN)           :: name,eltype
  INTEGER,           INTENT(IN), OPTIONAL :: ngr_act
! End interface
! Expected max length of name
  INTEGER, PARAMETER :: namtl=12
  DOUBLE PRECISION   :: cnvele(6), princ
  INTEGER            ::ln,nb,i
  CHARACTER*(namtl)  :: blanks
  INTEGER            ::lench
  EXTERNAL lench
  CHARACTER*80       ::elefmt
  INTEGER            ::ngr_act_loc
!============================================================
! Name
  ln=MAX(1,lench(name))
  nb=MAX(1,namtl-ln)
  blanks=' '
! Check if there are non-grav parameters
  IF(PRESENT(ngr_act))THEN
     ngr_act_loc=ngr_act
  ELSE
     ngr_act_loc=0
  END IF
! Convert to degrees
  cnvele=elem
  IF(eltype .eq. 'CAR')THEN
     elefmt='('''''''',A,'''''''',A,1x,F13.6,'//                    &
     &        '6(1x,F17.13),'//                                     &
     &        '2F6.2,1X,I2)'
  ELSEIF(eltype .eq. 'KEP')THEN
     cnvele(3:6)=cnvele(3:6)*degrad
     elefmt='('''''''',A,'''''''',A,F13.6,'//                       &
     &        '1p,6(e25.16),'//                                     &
     &        '0p,2F6.2,1X,I2)'
!        elefmt='('''''''',A,'''''''',A,F13.6,'//
!    +        'F16.12,1x,F12.10,4(1x,F12.7),'//
!    +        '2F6.2)'
  ELSEIF(eltype .eq. 'COM')THEN
     cnvele(3:5)=cnvele(3:5)*degrad
     elefmt='('''''''',A,'''''''',A,F13.6,'//                       &
     &        '1p,5(e25.16),0p,f14.6,'//                            &
     &        '0p,2F6.2,1X,I2)'
  ELSEIF(eltype .eq. 'COT')THEN
     cnvele(3:6)=cnvele(3:6)*degrad
     elefmt='('''''''',A,'''''''',A,F13.6,'//                       &
     &        '1p,6(e25.16),0p,2F6.2,1X,I2)'
  ELSEIF(eltype .eq. 'EQU')THEN
     cnvele(6)=princ(cnvele(6))
     cnvele(6)=cnvele(6)*degrad
     elefmt='('''''''',A,'''''''',A,F13.6,'//                       &
     &        'F16.12,4(1x,f12.9),1x,F12.7,'//                      &
     &        '2F6.2,1X,I2)'
  ELSEIF(eltype.eq.'ATT')THEN
     WRITE(*,*)' wro1lr: warning, the ATT elelments in  this format have no station code '
     cnvele(1:4)=cnvele(1:4)*degrad
     elefmt='('''''''',A,'''''''',A,F13.6,'//                       &
     &        '1p,6(e25.16),'//                                     &
     &        '0p,2F6.2,1X,I2)'
  ELSE
     WRITE(*,*)'*** wro1lr: unsupported coordinate type', eltype
     STOP '*** wro1lr: unsupported coordinate type'
  ENDIF
! Write
  IF(h.GT.-100.d0) THEN
     WRITE(unit,elefmt) name(1:ln),blanks(1:nb),t0,cnvele,h,g,ngr_act_loc
  ELSE
     WRITE(unit,elefmt) name(1:ln),blanks(1:nb),t0,cnvele,ngr_act_loc
  ENDIF
END SUBROUTINE wro1lr

!====================================================================!
! WRO1LR_MATLAB                                                      !
!====================================================================!
! It writes the files with orbital elements for MATLAB (1L)          !
!====================================================================!
SUBROUTINE wro1lr_matlab(unit,name,elem,eltype,t0,h,ngr_act,chi,rms,succ)
  USE fund_const
  IMPLICIT NONE
  !=======================================================================================================
  INTEGER,           INTENT(IN)           :: unit     ! File unit
  CHARACTER*(*),     INTENT(IN)           :: name     ! Name of the object
  DOUBLE PRECISION,  INTENT(IN)           :: elem(6)  ! Orbital elements
  CHARACTER*(*),     INTENT(IN)           :: eltype   ! Coordinate of the orbital elements
  DOUBLE PRECISION,  INTENT(IN)           :: t0       ! Time of the orbital elements
  DOUBLE PRECISION,  INTENT(IN)           :: h        ! Absolute magnitude
  INTEGER,           INTENT(IN), OPTIONAL :: ngr_act  ! Non grav. paramerters activated
  DOUBLE PRECISION,  INTENT(IN), OPTIONAL :: chi      ! Chi value
  DOUBLE PRECISION,  INTENT(IN), OPTIONAL :: rms      ! RMS value
  INTEGER,           INTENT(IN), OPTIONAL :: succ     ! Succ flag for the AR connected components
  !=======================================================================================================
  INTEGER, PARAMETER :: namtl=12  ! Expected max length of name
  DOUBLE PRECISION   :: cnvele(6) ! Elements to be written
  DOUBLE PRECISION   :: sma       ! Semi-major axis, for COM write
  DOUBLE PRECISION   :: princ     ! Function for the principal part of an angle
  INTEGER            :: ln,nb,i
  CHARACTER*(namtl)  :: blanks
  CHARACTER*80       :: elefmt
  INTEGER            :: ngr_act_loc
  DOUBLE PRECISION   :: chi_loc
  DOUBLE PRECISION   :: rms_loc
  INTEGER            :: succ_loc
  INTEGER            :: lench
  EXTERNAL lench
  !=======================================================================================================
  ln = MAX(1,lench(name))
  nb = MAX(1,namtl-ln)
  blanks = ' '
  ! Fill the optional variables
  IF(PRESENT(ngr_act))THEN
     ngr_act_loc = ngr_act
  ELSE
     ngr_act_loc = 0
  END IF
  IF(PRESENT(chi))THEN
     chi_loc = chi
  ELSE
     chi_loc = 0.d0
  END IF
  IF(PRESENT(rms))THEN
     rms_loc = rms
  ELSE
     rms_loc = 0.d0
  END IF
  IF(PRESENT(succ))THEN
     succ_loc = succ
  END IF
  ! Convert radians to degrees and inizialization of the format
  cnvele = elem
  IF(eltype.EQ.'CAR')THEN
     elefmt = '(A,A,F13.6,1P,6(F17.13),0P,F6.2,1X,I2,1P,1X,D12.5,1X,D12.5)'
  ELSEIF(eltype.EQ.'KEP')THEN
     cnvele(3:6) = cnvele(3:6)*degrad
     IF(PRESENT(succ)) THEN
        elefmt = '(A,A,F13.6,1P,6(E25.16),0P,F6.2,1X,I2,1P,1X,D12.5,1X,D12.5,1X,I1)'
     ELSE
        elefmt = '(A,A,F13.6,1P,6(E25.16),0P,F6.2,1X,I2,1P,1X,D12.5,1X,D12.5)'
     END IF
  ELSEIF(eltype.EQ.'COM')THEN
     cnvele(3:5) = cnvele(3:5)*degrad
     ! In this case write also the semi-major axis
     sma = cnvele(1)/(1-cnvele(2))
     IF(PRESENT(succ))THEN
        elefmt = '(A,A,F13.6,1P,6(E25.16),0P,F14.6,0P,F6.2,2X,I2,1P,2X,D12.5,1X,D12.5,5X,I1)'
     ELSE
        elefmt = '(A,A,F13.6,1P,6(E25.16),0P,F14.6,0P,F6.2,1X,I2,1P,1X,D12.5,1X,D12.5)'
     END IF
  ELSEIF(eltype.EQ.'COT')THEN
     cnvele(3:6) = cnvele(3:6)*degrad
     elefmt = '(A,A,F13.6,1P,6(E25.16),0P,F6.2,1X,I2,1P,1X,D12.5,1X,D12.5,1X,I1)'
  ELSEIF(eltype.EQ.'EQU')THEN
     cnvele(6) = princ(cnvele(6))
     cnvele(6) = cnvele(6)*degrad
     elefmt = '(A,A,F13.6,F16.12,4(1X,F12.9),1X,F12.7,F6.2,1X,I2,1P,1X,D12.5,1X,D12.5)'
  ELSEIF(eltype.EQ.'ATT')THEN
     WRITE(*,*)' wro1lr_matlab: WARNING! The ATT elements in this format have no station code '
     cnvele(1:4) = cnvele(1:4)*degrad
     elefmt = '(A,A,F13.6,1P,6(E25.16),0P,F6.2,1X,I,1P,1X,D12.52,1X,D12.5)'
  ELSE
     WRITE(*,*) 'wro1lr_matlab: ERROR! Unsupported coordinate type ', eltype
     STOP
  ENDIF
  ! Write the orbital elements
  IF(PRESENT(succ))THEN
     IF(eltype.EQ.'COM')THEN
        WRITE(unit,elefmt) name(1:ln),blanks(1:nb),t0,sma,cnvele,h,ngr_act_loc, chi_loc, rms_loc, succ_loc
     ELSE
        WRITE(unit,elefmt) name(1:ln),blanks(1:nb),t0,cnvele,h,ngr_act_loc, chi_loc, rms_loc, succ_loc
     END IF
  ELSE
     IF(eltype.EQ.'COM')THEN
        WRITE(unit,elefmt) name(1:ln),blanks(1:nb),t0,sma,cnvele,h,ngr_act_loc, chi_loc, rms_loc
     ELSE
        WRITE(unit,elefmt) name(1:ln),blanks(1:nb),t0,cnvele,h,ngr_act_loc, chi_loc, rms_loc
     END IF
  END IF
END SUBROUTINE wro1lr_matlab

!------------------------------------------------!
! SUBROUTINE write_mtp_header                    !
! It writes the header of the file .mtp.         !
!------------------------------------------------!
! Author: D. Bracali Cioci                       !
! Last Update: Feb 2020                          !
!------------------------------------------------!
SUBROUTINE write_mtp_header(iun_mtp,name)

  INTEGER,       INTENT(IN) :: iun_mtp  ! unit of file .mtp
  CHARACTER*(*), INTENT(IN) :: name     ! name of the object
  ! Local variables
  INTEGER :: ln
  INTEGER :: lench
  EXTERNAL lench

  ! Length of object name
  ln = MAX(1,lench(name))

  ! Write header
  WRITE(iun_mtp,600) name(1:ln)

600 FORMAT('% Name = ',A/'%'/&
         & '% i   j   t_cla         csi           zeta          d_cla        H')

END SUBROUTINE write_mtp_header


!----------------------------------------------!
! SUBROUTINE write_mtp_record                  !
! It writes a record of the file .mtp.         !
!----------------------------------------------!
! Author: D. Bracali Cioci                     !
! Last Update: Feb 2020                        !
!----------------------------------------------!
SUBROUTINE write_mtp_record(iun_mtp,i_ind,j_ind,t_cla,csi,zeta,d_cla,abs_mag)

  INTEGER,          INTENT(IN) :: iun_mtp      ! unit of file .mtp
  INTEGER,          INTENT(IN) :: i_ind,j_ind  ! indices of orbit in AR sampling
  REAL(KIND=dkind), INTENT(IN) :: t_cla        ! time of closest approach
  REAL(KIND=dkind), INTENT(IN) :: csi,zeta     ! coordinates on MTP
  REAL(KIND=dkind), INTENT(IN) :: d_cla        ! distance from Earth centre of mass
  REAL(KIND=dkind), INTENT(IN) :: abs_mag      ! absolute magnitude

  ! Write record
  WRITE(iun_mtp,601) i_ind,j_ind,t_cla,csi,zeta,d_cla,abs_mag

601 FORMAT(I3,2X,I3,1X,F12.5,1P,2(2X,E12.5),2X,E12.5,2X,0P,F5.2)

END SUBROUTINE write_mtp_record


!---------------------------------------------------!
! SUBROUTINE write_arsample_header                  !
! It writes the header of the file .ar_sample.      !
!---------------------------------------------------!
! Author: D. Bracali Cioci                          !
! Last Update: Feb 2020                             !
!---------------------------------------------------!
SUBROUTINE write_arsample_header(iun_ars,name,attr_coo,coeff,roots,rhodot_min,rhodot_max,type_sample)

  INTEGER,           INTENT(IN) :: iun_ars     ! unit of file .ar_sample
  CHARACTER*(*),     INTENT(IN) :: name        ! name of the object
  REAL(KIND=dkind),  INTENT(IN) :: attr_coo(4) ! attributable vector
  REAL(KIND=dkind),  INTENT(IN) :: coeff(0:5)  ! coefficients for the construction of V(r) (see cobweb.f90)
  REAL(KIND=dkind),  INTENT(IN) :: roots(3)    ! positive roots of the 6-th degree polynomial V(r)
  REAL(KIND=dkind),  INTENT(IN) :: rhodot_min  ! minimum range-rate of AR
  REAL(KIND=dkind),  INTENT(IN) :: rhodot_max  ! maximum range-rate of AR
  CHARACTER(LEN=18), INTENT(IN) :: type_sample ! type of AR sampling and scale in rho axis
  ! Local variables
  INTEGER :: ln
  INTEGER :: lench
  EXTERNAL lench

  ! Length of object name
  ln = MAX(1,lench(name))

  ! Write header
  WRITE(iun_ars,604) name(1:ln),attr_coo,coeff,roots,rhodot_min,rhodot_max,type_sample
  WRITE(iun_ars,'(A1)') '%'
  WRITE(iun_ars,'(A116)') '% i   j   rho               rho_dot           cc    succ   '// &
       &                  'chi          Jac (AR-MOV)     E_Sun         E_Earth       H'

  604 FORMAT('% Name            = ',A/'% Attributable    =',2(F10.7,1X),1P,E12.5,1X,E12.5/   &
         & '% Coefficients    =',6(E12.5,1X)/'% Roots           =',3(E12.5,1X),0P,' [au]'/ &
         & '% Min/Max rho_dot =',2(F15.6,1X),' [au/day]'/'% Sample          = ',A)

END SUBROUTINE write_arsample_header


!-------------------------------------------------!
! SUBROUTINE write_arsample_record                !
! It writes one record of the file .ar_sample.    !
!-------------------------------------------------!
! Author: D. Bracali Cioci                        !
! Last Update: Feb 2020                           !
!-------------------------------------------------!
SUBROUTINE write_arsample_record(iun_ars,i_ind,j_ind,rho,rho_dot,conn_comp,succ,chi,jac,en_sun,en_earth,hmag)

  INTEGER,          INTENT(IN) :: iun_ars     ! unit of file .ar_sample
  INTEGER,          INTENT(IN) :: i_ind,j_ind ! indices of orbit in AR sampling
  REAL(KIND=dkind), INTENT(IN) :: rho,rho_dot ! range and range-rate
  INTEGER,          INTENT(IN) :: conn_comp   ! AR connected component
  LOGICAL,          INTENT(IN) :: succ        ! diff. cor. convergence flag
  REAL(KIND=dkind), INTENT(IN) :: chi         ! chi of MOV fit
  REAL(KIND=dkind), INTENT(IN) :: jac         ! jacobian of AR to MOV map
  REAL(KIND=dkind), INTENT(IN) :: en_sun      ! two-body heliocentric energy
  REAL(KIND=dkind), INTENT(IN) :: en_earth    ! two-body geocentric energy
  REAL(KIND=dkind), INTENT(IN) :: hmag        ! absolute magnitude
  ! Local variables
  INTEGER :: int_flag

  int_flag = 0
  IF(succ) int_flag = 1

  ! Write record
  WRITE(iun_ars,605) i_ind,j_ind,rho,rho_dot,conn_comp,int_flag,chi,jac,en_sun,en_earth,hmag

605 FORMAT(I3,2X,I3,1X,1P,2(E17.10,1X),1X,I1,5X,I1,4X,2(1X,1PE12.5),3X,2(2X,E12.5),0P,2X,F8.4)

END SUBROUTINE write_arsample_record



!--------------------------------------------------!
! SUBROUTINE write_movsample_header                !
! It writes the header of the file .mov_sample.    !
!--------------------------------------------------!
! Author: D. Bracali Cioci                         !
! Last Update: Feb 2020                            !
!--------------------------------------------------!
SUBROUTINE write_movsample_header(iun_movs,name,epoch,coo)

  INTEGER,          INTENT(IN) :: iun_movs  ! unit of file .mov_sample
  CHARACTER*(*),    INTENT(IN) :: name      ! Name of the object
  REAL(KIND=dkind), INTENT(IN) :: epoch     ! epoch of orbit
  CHARACTER(LEN=3), INTENT(IN) :: coo       ! type of coordinates
  ! Local variables
  INTEGER :: ln
  INTEGER :: lench
  EXTERNAL lench

  ! Length of object name
  ln = MAX(1,lench(name))

  ! Wrire header
  WRITE(iun_movs,602) name(1:ln),epoch,coo

602 FORMAT('% Name        =  ',A/'% Epoch [MJD] = ',F12.5,' TDT',/'% Coordinates =  ',A3/'% Ref. Sys.   =  ECLM J2000'/'%')

  WRITE(iun_movs,'(A)') '% i   j     q [au]                   e                        I [deg]'//    &
       &         '                  long. node [deg]         arg. peric. [deg]       peric. time'// &
       &         '    H      chi          RMS          cc   MOID [au]     V_inf [km/s]  imp  geoc'

END SUBROUTINE write_movsample_header

!----------------------------------------------------!
! SUBROUTINE write_movsample_record                  !
! It writes one record of the file .mov_sample.      !
!----------------------------------------------------!
! Author: D. Bracali Cioci                           !
! Last Update: Feb 2020                              !
!----------------------------------------------------!
SUBROUTINE write_movsample_record(iun_movs,i_ind,j_ind,el_orb,chi,rms,cc_flag,moid,v_infty,imp_flag,e_earth)

  INTEGER,          INTENT(IN) :: iun_movs     ! unit of file .mov_sample
  INTEGER,          INTENT(IN) :: i_ind,j_ind  ! indices of orbit in AR sampling
  TYPE(orbit_elem), INTENT(IN) :: el_orb       ! orbital elements
  REAL(KIND=dkind), INTENT(IN) :: chi          ! chi of MOV fit
  REAL(KIND=dkind), INTENT(IN) :: rms          ! rms of MOV fit
  INTEGER,          INTENT(IN) :: cc_flag      ! flag of connected component
  REAL(KIND=dkind), INTENT(IN) :: moid         ! MOID
  REAL(KIND=dkind), INTENT(IN) :: v_infty      ! velocity at infinity
  INTEGER,          INTENT(IN) :: imp_flag     ! impacting orbit
  REAL(KIND=dkind), INTENT(IN) :: e_earth      ! energy w.r.t. Earth
  ! Local variables
  TYPE(orbit_elem) :: elcom     ! Cometarian elements of the sampling points
  INTEGER          :: fail_flag ! failure of the coordinate change
  INTEGER          :: geoc      ! flag for geocentric orbit

  CALL coo_cha(el_orb,'COM',elcom,fail_flag)
  IF(fail_flag.GT.0)THEN
     WRITE(*,*) 'write_movsample_record: failed conversion to COM. Fail flag = ', fail_flag
  ELSE
     geoc = 0
     IF(e_earth.LE.0.d0) geoc = 1
     elcom%coord(3:5) = elcom%coord(3:5)*degrad
     ! Write record
     WRITE(iun_movs,603) i_ind,j_ind,elcom%coord(1:6),elcom%h_mag,chi,rms,cc_flag,moid,v_infty,imp_flag,geoc
  END IF

603 FORMAT(I3,2X,I3,1X,1P,5(E25.16),0P,F14.6,2X,F6.2,2(1X,1PE12.5),0P,2X,I1,3X,F10.7,4X,F7.3,8X,I1,4X,I1)

END SUBROUTINE write_movsample_record


!--------------------------------------------------!
! SUBROUTINE read_movsample_header                 !
! It reads the header of the file .mov_sample.     !
!--------------------------------------------------!
! Author: D. Bracali Cioci                         !
! Last Update: Mar 2020                            !
!--------------------------------------------------!
SUBROUTINE read_movsample_header(iun_movs,name,epoch,coo)

  INTEGER,          INTENT(IN)  :: iun_movs  ! unit of file .mov_sample
  CHARACTER*(*),    INTENT(OUT) :: name      ! name of the object
  REAL(KIND=dkind), INTENT(OUT) :: epoch     ! epoch of orbit
  CHARACTER(LEN=3), INTENT(OUT) :: coo       ! type of coordinates
  ! Local variables
  INTEGER            :: status
  CHARACTER(LEN=200) :: aux,aux1

  READ(iun_movs,'(A17,A7)',IOSTAT=status) aux,name
  READ(iun_movs,'(A17,F12.5,A4)',IOSTAT=status) aux,epoch,aux1
  READ(iun_movs,'(A17,A3)',IOSTAT=status) aux,coo
  IF(status.GT.0)THEN
     WRITE(*,*) 'read_movsample_record: error in reading header of file .mov_sample'
     WRITE(ierrou,*) 'read_movsample_record: error in reading header of file .mov_sample'
     STOP
  END IF
  READ(iun_movs,*)
  READ(iun_movs,*)
  READ(iun_movs,*)

END SUBROUTINE read_movsample_header


!---------------------------------------------------!
! SUBROUTINE read_movsample_record                  !
! It reads one record of the file .mov_sample.      !
!---------------------------------------------------!
! Author: D. Bracali Cioci                          !
! Last Update: Mar 2020                             !
!---------------------------------------------------!
SUBROUTINE read_movsample_record(iun_movs,i_ind,j_ind,elcom,chi,rms,cc_flag, &
     &                           moid,v_infty,imp,geoc,eof)

  INTEGER,          INTENT(IN)  :: iun_movs     ! unit of file .mov_sample
  INTEGER,          INTENT(OUT) :: i_ind,j_ind  ! indices of orbit in AR sampling
  TYPE(orbit_elem), INTENT(OUT) :: elcom        ! orbital elements
  REAL(KIND=dkind), INTENT(OUT) :: chi          ! chi of MOV fit
  REAL(KIND=dkind), INTENT(OUT) :: rms          ! rms of MOV fit
  INTEGER,          INTENT(OUT) :: cc_flag      ! flag of connected component
  REAL(KIND=dkind), INTENT(OUT) :: moid         ! MOID
  REAL(KIND=dkind), INTENT(OUT) :: v_infty      ! velocity at infinity
  INTEGER,          INTENT(OUT) :: imp          ! impacting orbit
  INTEGER,          INTENT(OUT) :: geoc         ! geocentric orbit
  LOGICAL,          INTENT(OUT) :: eof          ! flag for end of file
  ! Local variables
  INTEGER :: status

  eof = .FALSE.
  READ(iun_movs,605,IOSTAT=status,END=100) i_ind,j_ind,elcom%coord(1:6),elcom%h_mag,chi,rms,cc_flag,moid,v_infty,imp,geoc
  IF(status.GT.0)THEN
     WRITE(*,*) 'read_movsample_record: error in reading one record of file .mov_sample'
     WRITE(ierrou,*) 'read_movsample_record: error in reading one record of file .mov_sample'
     STOP
  END IF
  elcom%coord(3:5) = elcom%coord(3:5)*radeg
  elcom%mag_set    = .TRUE.
  RETURN
100 eof = .TRUE.
605 FORMAT(I3,2X,I3,1X,1P,5(E25.16),0P,F14.6,2X,F6.2,2(1X,1PE12.5),0P,2X,I1,3X,F10.7,4X,F7.3,8X,I1,4X,I1)

END SUBROUTINE read_movsample_record




!  *****************************************************************
!  *                                                               *
!  *                         W R O 1 L H 2                         *
!  *                                                               *
!  *       Writes the header of an orbital element file            *
!  *   (single-line format, different epochs, elements)            *
!  *                                                               *
!  *****************************************************************
!
! INPUT:    UNIT      -  Output FORTRAN unit
!           RSYS      -  Reference system type (EQUM/EQUT/ECLM)
!           EPOCH     -  Reference system epoch (J2000/OFDATE)
!           ELTYPE    -  Type of orbital elements (EQU/KEP/CAR)
!           VERS      -  Version (old: OEF1.1 or new: OEF2.0 )
!

SUBROUTINE wro1lh2(unit,rsys,epoch,eltype,version)
  IMPLICIT NONE
! Begin interface
  INTEGER,      INTENT(IN)           :: unit
  CHARACTER*(*),INTENT(IN)           :: rsys,epoch,eltype
  CHARACTER*6,  INTENT(IN), OPTIONAL :: version
! End interface
  INCLUDE 'parcmc.h90'
  INTEGER l1,l2,l3,nb
  CHARACTER*100 bl
  CHARACTER*6 :: vers
  INTEGER lench
  EXTERNAL lench
!========================================================
  ! Default is OEF2.0
  IF(PRESENT(version)) THEN
     vers=version
  ELSE
     vers='OEF2.0'
  END IF
  l1=lench(rsys)
  l2=lench(epoch)
  l3=lench(eltype)
  nb=MAX(14-l1-l2,1)
  bl=' '
  WRITE(unit,100)     vers(1:6),comcha,comcha,eltype(1:l3),comcha,      &
       &                rsys(1:l1),epoch(1:l2),                         &
       &                bl(1:nb),comcha
100 FORMAT("format  = '",A,"'       ",A,' file format'/                     &
         &       'rectype = ''1L''           ',A,' record type (1L/ML)'/    &
         &       'elem    = ''',A,'''          ',A,                         &
         &                   ' type of orbital elements'/                   &
         &       'refsys  = ',A,1X,A,A,A,' default reference system'/       &
         &       'END_OF_HEADER')

  IF(eltype.EQ.'KEP') THEN
     WRITE(unit,201) comcha
201  FORMAT(A,' Name, Epoch(MJD), a, e, i, long. node,',' arg. peric., mean anomaly, ' &
          & 'absolute magnitude, slope param, non-grav param')
  ELSEIF(eltype.EQ.'CAR') THEN
     WRITE(unit,202) comcha
202  FORMAT(A,' Name, Epoch(MJD), cartesian position and velocity',' vectors, ' &
          & 'absolute magnitude, slope param, non-grav param')
  ELSEIF(eltype.EQ.'EQU') THEN
     WRITE(unit,203) comcha
203  FORMAT(A,' Name, Epoch(MJD), a, e*sin(LP), e*cos(LP),',' tan(i/2)*sin(LN), tan(i/2)*cos(LN), mean long., ' &
          & 'absolute magnitude, slope param, non-grav param')
  ELSEIF(eltype.EQ.'COM') THEN
     WRITE(unit,204) comcha
204  FORMAT(A,' Name, Epoch(MJD), q, e, i, long. node,',' arg. peric., perihelion time, ' &
          & 'absolute magnitude, slope param, non-grav param')
  ELSEIF(eltype.EQ.'COT') THEN
     WRITE(unit,205) comcha
205  FORMAT(A,' Name, Epoch(MJD), q, e, i, long. node,',' arg. peric., true anomaly, ' &
          & 'absolute magnitude, slope param, non-grav param')
  ELSEIF(eltype.eq.'ATT')THEN
     WRITE(unit,206) comcha
206  FORMAT(A,' Name, Epoch(MJD), R.A., DEC, R.A.dot, DECdot, r, rdot, ' &
          & 'absolute magnitude, slope param, non-grav param')
  END IF
END SUBROUTINE wro1lh2

!  *****************************************************************
!  *                                                               *
!  *                   W R O 1 L H 2 _ N G R                       *
!  *                                                               *
!  *       Writes the header of a file containing non-grav param   *
!  *              in case of a single line catalog                 *
!  *                                                               *
!  *****************************************************************
!
! INPUT:    UNIT      -  Output FORTRAN unit
!           RSYS      -  Reference system type (EQUM/EQUT/ECLM)
!           EPOCH     -  Reference system epoch (J2000/OFDATE)
!           VERS      -  Version ( OEF2.0 )
!

SUBROUTINE wro1lh2_ngr(unit,rsys,epoch,vers,n_dp)
  IMPLICIT NONE
! Begin interface
  INTEGER,      INTENT(IN)            :: unit
  CHARACTER*(*),INTENT(IN)            :: rsys,epoch
  CHARACTER*6,  INTENT(IN)            :: vers
  INTEGER,      INTENT(IN), OPTIONAL  :: n_dp
! End interface
  INCLUDE 'parcmc.h90'
  INTEGER       :: l1,l2,nb
  CHARACTER*100 :: bl
  INTEGER       :: lench
  EXTERNAL lench
!======================================================
  l1=lench(rsys)
  l2=lench(epoch)
  nb=MAX(14-l1-l2,1)
  bl=' '
  WRITE(unit,500)     vers(1:6),comcha,comcha,                              &
       &              rsys(1:l1),epoch(1:l2),                               &
       &              bl(1:nb),comcha
500 FORMAT("format  = '",A,"'       ",A,' file format'/                     &
         &       'rectype = ''1L''           ',A,' record type (1L/ML)'/    &
         &       'refsys  = ',A,1X,A,A,A,' default reference system'/       &
         &       'END_OF_HEADER')

  IF(PRESENT(n_dp)) THEN
     IF(n_dp.EQ.2) THEN
        WRITE(unit,501) comcha
501     FORMAT(A,' Name, dynamical model, actual number of parameters in use, A/M, A2 ')
     ELSE IF(n_dp.EQ.3) THEN
        WRITE(unit,502) comcha
502     FORMAT(A,' Name, dynamical model, actual number of parameters in use, A1, A2, A3 ')
     ELSE IF(n_dp.EQ.4) THEN
        WRITE(unit,503) comcha
503     FORMAT(A,' Name, dynamical model, actual number of parameters in use, A1, A2, A3, delay ')
     ELSE IF(n_dp.EQ.1) THEN
        WRITE(unit,504) comcha
504     FORMAT(A,' Name, dynamical model, actual number of parameters in use, eta ')
     ELSE IF(n_dp.NE.0) THEN
        WRITE(*,*) 'wro1lh2_ngr: n_dp, wrong value ', n_dp
        STOP
     END IF
  ELSE
     WRITE(unit,505) comcha
505     FORMAT(A,' Name, dynamical model, actual number of parameters in use, A/M, A2' &
             & '            (asteroid case)'/                                          &
             & '!                                                             A1, A2'  &
             & ', A3, delay (comet case)')
  END IF
END SUBROUTINE wro1lh2_ngr


! ---------------------------------------------------------------------
!
!  *****************************************************************
!  *                                                               *
!  *                         W R O M L H 2                           *
!  *                                                               *
!  *       Writes the header of an orbital element file            *
!  *                     (multi-line format)                       *
!  *                                                               *
!  *****************************************************************
!
! OUTPUT:   UNIT      -  Output FORTRAN unit
!           RSYS      -  Reference system type (EQUM/EQUT/ECLM)
!           EPOCH     -  Reference system epoch (J2000/OFDATE)
!           VERS      -  Version (OEF2.0)
!
SUBROUTINE wromlh2(unit,rsys,epoch,version)
  IMPLICIT NONE
! Begin interface
  INTEGER,       INTENT(IN)           :: unit
  CHARACTER*(*), INTENT(IN)           :: rsys,epoch
  CHARACTER*6,   INTENT(IN), OPTIONAL :: version
! End interface
  INCLUDE 'parcmc.h90'
  INTEGER       :: l1,l2,nb
  CHARACTER*100 :: bl
  CHARACTER*6   :: vers
  INTEGER       :: lench
  EXTERNAL lench
!=============================================================
  IF(PRESENT(version)) THEN
     vers=version
  ELSE
     vers='OEF2.0'
  END IF
  l1=lench(rsys)
  l2=lench(epoch)
  nb=MAX(14-l1-l2,1)
  bl=' '
  WRITE(unit,100) vers(1:6),comcha,comcha,rsys(1:l1),epoch(1:l2),       &
       &                bl(1:nb),comcha
100 FORMAT("format  = '",A6,"'       ",A,' file format'/                    &
         &       'rectype = ''ML''           ',A,' record type (1L/ML)'/    &
         &       'refsys  = ',A,1X,A,A,A,' default reference system'/       &
         &       'END_OF_HEADER')
END SUBROUTINE wromlh2


! Copyright (C) 1997-2003 by Mario Carpino (carpino@brera.mi.astro.it),
!                            Andrea Milani (milani@dm.unipi.it),
!                            Zoran Knezevic (zoran@aob.aob.bg.ac.yu)
! Version: 12 August 2003
! ---------------------------------------------------------------------
!
!  *****************************************************************
!  *                                                               *
!  *                         R D E L E M                           *
!  *                                                               *
!  *          Read orbital elements for a list of objects          *
!  *                     from a list of files                      *
!  *                                                               *
!  *****************************************************************
!
! INPUT:    UNIT      -  Output FORTRAN unit (report file)
!           OBJNAM    -  Object names
!           NOBJ      -  Number of objects
!           INFILES   -  Input files
!           NFIL      -  Number of input files
!           DEFORB    -  Is the orbit already defined?
!
! OUTPUT:   DEFORB    -  Was the orbit found?
!           DEFCN     -  Tells whether covariance/normal matrices
!                            are defined
!           ELTYPE    -  Type of orbital elements (EQU/KEP/CAR)
!           TELEM     -  Epoch of orbital elements (MJD, TDT)
!           ELEM      -  Orbital elements (ECLM J2000)
!           COVE      -  Covariance matrix of orbital elements
!           NORE      -  Normal matrix of orbital elements
!           MASS      -  Mass (solar masses)
!           HMAG      -  H absolute magnitude (if <-100, missing)
!           GMAG      -  G slope parameter
!           COMELE    -  Comment on orbital elements
!
! WARNING: the routine is designed to allow multiple calls on
!          different lists of files/objects: for this reason,
!          orbital elements of objects which have already an
!          orbit defined (DEFORB=.true.) ARE NOT modified.
!          For this reason, the user must define suitably DEFORB
!          BEFORE calling RDELEM
!
SUBROUTINE rdelem(unit,objnam,nobj,infiles,nfil,deforb,defcn,      &
     &                  eltype,telem,elem,cove,nore,              &
     &                  mass,h,g,comele)

  INTEGER nfilx
  PARAMETER(nfilx=20)
  INTEGER,          INTENT(IN)    :: unit,nobj,nfil
  CHARACTER(LEN=*), INTENT(IN)    :: objnam(nobj),infiles(nfil)
  LOGICAL,          INTENT(INOUT) :: deforb(nobj)
  LOGICAL,          INTENT(OUT)   :: defcn(nobj)
  CHARACTER(LEN=*), INTENT(OUT)   :: eltype(nobj),comele(nobj)
  DOUBLE PRECISION, INTENT(OUT)   :: telem(nobj),elem(6,nobj)
  DOUBLE PRECISION, INTENT(OUT)   :: cove(ndimx,ndimx,nobj), nore(ndimx,ndimx,nobj)
  DOUBLE PRECISION, INTENT(OUT)   :: mass(nobj),h(nobj),g(nobj)
!*********************************************************************************
  INTEGER   :: i,lf,is1,lfo,uniin,n
  CHARACTER :: form*10,infil1*100
  LOGICAL   :: opened
  INTEGER   :: lench
  EXTERNAL lench

! Number of input files
  IF(nfil.GT.nfilx) STOP '**** rdelem: nfil > nfilx ****'

! Nothing to do
  IF(nfil.LE.0) RETURN

! Loop on files
  DO 2 i=1,nfil
     lf=lench(infiles(i))
     opened=.false.

! Understand file format, if given in the argument infils(i)
     is1=index(infiles(i)(1:lf),'[')
     IF(is1.GT.0 .AND. infiles(i)(lf:lf).EQ.']') THEN
        form=infiles(i)(is1+1:lf-1)
        lf=is1-1
        infil1=infiles(i)(1:lf)
     ELSE
        infil1=infiles(i)
        form=' '
        CALL filopn(uniin,infil1(1:lf),'OLD')
        opened=.true.
        CALL oefdet(uniin,infil1(1:lf),form) ! SELF-DETECT FILE FORMAT
        REWIND(uniin)
     END IF
     lfo=lench(form)
     IF(form.NE.' ') WRITE(*,102) infil1(1:lf),form(1:lfo)
102  FORMAT('Scanning file "',A,'" (format: ',A,')')

! Reading file
     IF(form.EQ.'OEF') THEN
        IF(.NOT.opened) CALL filopn(uniin,infil1,'OLD')
        CALL rdoef(uniin,infil1,objnam,nobj,deforb,defcn,eltype,telem,&
             &     elem,cove,nore,mass,h,g,comele)
        CALL clorbf
     ELSEIF(form.EQ.'BA1') THEN
        WRITE(*,*)' rdelem: obsolete Lowell format pre-1999'
        STOP
!        IF(.NOT.opened) CALL filopn(uniin,infil1,'OLD')
!        CALL rdast1(uniin,infil1,objnam,nobj,deforb,defcn,            &
!     &                eltype,telem,elem,cove,nore,mass,h,g,comele)
!        CALL filclo(uniin,' ')
     ELSEIF(form.EQ.'BA2') THEN
        IF(.NOT.opened) CALL filopn(uniin,infil1,'OLD')
        CALL rdast2(uniin,infil1,objnam,nobj,deforb,defcn,            &
             &                eltype,telem,elem,cove,nore,mass,h,g,comele)
        CALL filclo(uniin,' ')
     ELSEIF(form.EQ.'MPC-A') THEN
        IF(.NOT.opened) CALL filopn(uniin,infil1,'OLD')
        CALL rdmpca(uniin,infil1,objnam,nobj,deforb,defcn,            &
             &                eltype,telem,elem,cove,nore,mass,h,g,comele)
        CALL filclo(uniin,' ')
     ELSEIF(form.EQ.' ') THEN
        WRITE(*,120) infil1(1:lf)
        STOP '**** rdelem: abnormal end ****'
     ELSE
        WRITE(*,121) form(1:lfo),infil1(1:lf)
        STOP '**** rdelem: abnormal end ****'
     END IF
120  FORMAT('ERROR: unknown format type for file "',A,'"')
121  FORMAT('ERROR: unsupported format "',A,'" for file "',A,'"')
! set default for magnitude data
     DO n=1,nobj
        IF(deforb(n))THEN
           IF(h(n).lt.-1.d6)THEN
! leave it
           ENDIF
           IF(g(n).lt.-1.d6)THEN
              g(n)=0.15d0
           ENDIF
        ENDIF
     ENDDO
! End of loop on files
2 END DO

END SUBROUTINE rdelem

! Copyright (C) 1998 by Mario Carpino (carpino@brera.mi.astro.it),
! Version: June 19, 1998
! ---------------------------------------------------------------------
!
!  *****************************************************************
!  *                                                               *
!  *                          R D O E F                            *
!  *                                                               *
!  *          Read orbital elements for a list of objects          *
!  *         from a file written in internal ORBFIT format         *
!  *                                                               *
!  *****************************************************************
!
! INPUT:    UNIIN     -  Input unit (file is already opened)
!           FILE      -  Input file name (for error messages only)
!           OBJNAM    -  Object names (without embedded blanks)
!           NOBJ      -  Number of objects
!
! OUTPUT:   DEFORB    -  Tells whether orbital elements are defined
!           DEFCN     -  Tells whether covariance/normal matrices
!                            are defined
!           ELTYPE    -  Type of orbital elements (EQU/KEP/CAR/COM?ATT)
!           TELEM     -  Epoch of orbital elements (MJD, TDT)
!           ELEM      -  Orbital elements (ECLM J2000)
!           COVE      -  Covariance matrix of orbital elements
!           NORE      -  Normal matrix of orbital elements
!           MASS      -  Mass (solar masses)
!           H         -  H absolute magnitude (if <-100, missing)
!           G         -  G slope parameter
!           COMELE    -  Comment on orbital elements
!
! WARNING: the routine assumes that objects having DEFORB=.true.
!          have already orbital elements defined (possibly from another
!          input file) and does not overwrite them
!
SUBROUTINE rdoef(uniin,file,objnam,nobj,deforb,defcn,eltype,telem,&
     &                 elem,cove,nore,mass,h,g,comele)
  USE fund_const
  USE reference_systems
  USE name_rules
  IMPLICIT NONE
  ! BEGIN INTERFACE
  INTEGER, INTENT(IN) :: uniin
  INTEGER, INTENT(IN) :: nobj
  CHARACTER*(*), INTENT(IN) :: file,objnam(nobj)
  DOUBLE PRECISION, INTENT(OUT) :: telem(nobj),elem(6,nobj)
  DOUBLE PRECISION, INTENT(OUT) :: cove(ndimx,ndimx,nobj), nore(ndimx,ndimx,nobj)
  DOUBLE PRECISION, INTENT(OUT) :: mass(nobj),h(nobj),g(nobj)
  CHARACTER*(*), INTENT(OUT) :: eltype(nobj),comele(nobj)
  LOGICAL, INTENT(INOUT) :: deforb(nobj),defcn(nobj)
  ! END INTERFACE

  CHARACTER*(20) filedum
  INTEGER kr,k,j1,j2,lf,ln,lc,lc1,nrem
  DOUBLE PRECISION t1,gmsun,gma,gma1,enne,h1,g1,m1
  DOUBLE PRECISION rot(3,3)
  TYPE(orbit_elem) :: elem_loc
  TYPE(orb_uncert) :: unc, unc2
  DOUBLE PRECISION elem1(6),xv(6),cove1(ndimx,ndimx),cove2(ndimx,ndimx)
  DOUBLE PRECISION de(6,6),nore1(ndimx,ndimx),nore2(ndimx,ndimx)
  CHARACTER eltyp1*3,rsys*10,epoch*10,krc*10
  LOGICAL defcov,defnor,end,err
  CHARACTER*(idnamvir_len) name1,nc1
  INTEGER lench
  INTEGER :: nlsloc,lsloc(ndyx)
  !Dimensions
  INTEGER :: nd
  EXTERNAL lench

  gmsun=gms

! Number of remaining object (orbit not yet found)
  nrem=0
  DO 1 k=1,nobj
     IF(.NOT.deforb(k)) nrem=nrem+1
1 END DO
  IF(nrem.LE.0) RETURN

  CALL oporbf(file,uniin)
  lf=lench(file)

3 CONTINUE
!  CALL rdorb(name1,elem1,eltyp1,t1,cove1,defcov,nore1,defnor,       &
!     &           h1,g1,m1,rsys,epoch,kr,end)
  CALL read_elems(elem_loc,name1,end,nlsloc,err,lsloc,UNC=unc,FILE=file(1:lf),UNIT=uniin,MASS=m1,&
       & RSYS=rsys, EPOCH=epoch)
! Save dyn%dp, nls, and ls
  IF(dyn%ndp.GT.0.AND.nlsloc.NE.0.AND..NOT.ngr_opt)THEN
     nls=nlsloc
     ls=lsloc
  ENDIF
! Set dimension
  nd=6+nls

!  IF(end) GOTO 20
! Name match is performed disregarding embedded blanks
  nc1=name1
  CALL rmsp(nc1,lc1)
  IF(lc1.LE.0) GOTO 3
  DO 4 k=1,nobj
     IF(deforb(k)) GOTO 4
     IF(nc1(1:lc1).EQ.objnam(k)) THEN
        IF(rsys.EQ.'ECLM' .AND. epoch.EQ.'J2000') THEN
           eltype(k)=elem_loc%coo
           DO j1=1,nd
              DO  j2=1,nd
                 cove(j1,j2,k)=unc%g(j1,j2)
                 nore(j1,j2,k)=unc%c(j1,j2)
              ENDDO
           ENDDO
           DO j1=1,6
              elem(j1,k)=elem_loc%coord(j1)
           END DO
        ELSE
           gma=gmsun*m1
           gma1=gma+gmsun
           IF(unc%succ) THEN
! Transformation in cartesian coordinates
              CALL cooder(elem_loc%coord,elem_loc%coo,gma1,xv,'CAR',enne,de)
              !CALL convertcovdp(6,cove1,de,cove2)
              !CALL norprsdp(nore1,de,6,nore2,error)
              CALL convertunc(unc,de,unc2)
              !IF(error) THEN
              !   ln=lench(name1)
              !   WRITE(*,120) file(1:lf),name1(1:ln)
              !   defnor=.false.
              !END IF
! Transformation of reference system
              CALL rotpn(rot,rsys,epoch,elem_loc%t,'ECLM','J2000',0.d0)
              CALL prodmv(elem(1,k),rot,xv(1))
              CALL prodmv(elem(4,k),rot,xv(4))
              DO j1=1,3
                 DO j2=1,3
                    de(j1,j2)=rot(j1,j2)
                    de(j1+3,j2)=0
                    de(j1,j2+3)=0
                    de(j1+3,j2+3)=rot(j1,j2)
                 ENDDO
              ENDDO
               !IF(defcov) CALL convertcovdp(6,cove2,de,cove(1:6,1:6,k))
               !IF(defnor) THEN
               !   CALL norprsdp(nore2,de,6,nore(1,1,k),error)
               !   IF(error) THEN
               !      ln=lench(name1)
               !      WRITE(*,120) file(1:lf),name1(1:ln)
               !      defnor=.false.
               !   END IF
               !END IF
              IF(unc%succ)THEN
                 CALL convertunc(unc2,de,unc)
              END IF
           ELSE
              CALL coocha(elem_loc,eltyp1,gma1,xv,'CAR',enne)
              CALL rotpn(rot,rsys,epoch,t1,'ECLM','J2000',0.d0)
              CALL prodmv(elem(1,k),rot,xv(1))
              CALL prodmv(elem(4,k),rot,xv(4))
           END IF
           eltype(k)='CAR'
        END IF
        deforb(k)=.true.
        telem(k)=elem_loc%t
        IF(unc%succ)THEN
           defcov=.TRUE.
           defnor=.TRUE.
        END IF
        CALL fixcnm(defcov,defnor,defcn(k),cove(1,1,k),nore(1,1,k))
        mass(k)=m1
        h(k)=elem_loc%h_mag
        g(k)=elem_loc%g_mag
!        WRITE(krc,101) kr
 !       CALL rmsp(krc,lc)
        comele(k)='read from file "'//file(1:lf)//'"'
! &              '" at record '//krc(1:lc)
        nrem=nrem-1
        IF(nrem.LE.0) GOTO 20
     END IF
101  FORMAT(I6)
120  FORMAT(' rdoef: error in transforming normal matrix'/             &
          &       '        (file "',A,'", object "',A,'")')
4 END DO
  GOTO 3
20 CONTINUE
END SUBROUTINE rdoef


END MODULE orbit_elements

! =====================================================================
!
!  ***************************************************************
!  *                                                             *
!  *                         A P P M A G                         *
!  *                                                             *
!  *     Calcolo della magnitudine apparente di un asteroide     *
!  *                                                             *
!  ***************************************************************
!
!
! INPUT:    H         -  Magnitudine assoluta
!           G         -  Parametro di slope
!           DS        -  Distanza dal Sole (au)
!           DT        -  Distanza dall'osservatore (au)
!           BETA      -  Angolo di fase solare (rad): angolo tra il
!                        Sole e l'osservatore (Terra), visto dall'astero
!
DOUBLE PRECISION FUNCTION appmag(h,g,ds,dt,beta)
  IMPLICIT NONE
  DOUBLE PRECISION h,g,ds,dt,beta
  DOUBLE PRECISION a1,a2,b1,b2
  DOUBLE PRECISION tb2,phi1,phi2
  SAVE a1,a2,b1,b2
! Costanti per il calcolo della magnitudine
  DATA a1,a2,b1,b2/3.33d0,1.87d0,0.63d0,1.22d0/
  tb2=TAN(beta/2.d0)
  phi1=EXP(-a1*tb2**b1)
  phi2=EXP(-a2*tb2**b2)
  appmag=5.d0*LOG10(ds*dt)+h-2.5d0*LOG10((1.d0-g)*phi1+g*phi2)
END FUNCTION appmag

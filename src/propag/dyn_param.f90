MODULE dyn_param

USE fund_const

IMPLICIT NONE

PRIVATE

! public dimension parameters
INTEGER, PARAMETER, PUBLIC :: ndyx=4
INTEGER, PARAMETER, PUBLIC :: ndimx=6+ndyx

! data types

TYPE dyn_par
  REAL(KIND=dkind) :: dp(ndyx)              ! vector of dynamical parameters
  INTEGER          :: ndp                   ! actual number in use
  INTEGER          :: nmod                  ! nmod=0 no dyn. par, 
                                            ! nmod=1 dir.rad.p.+ Yarko; 
                                            ! nmod=2 A1,A2,A3, delay; 
                                            ! nmod=3 EP violation 
  LOGICAL          :: active(ndyx)          ! included in right hand side
  LOGICAL          :: apriori(ndyx)         ! flag for a priori
  DOUBLE PRECISION :: apriori_nominal(ndyx) ! a priori value 
  DOUBLE PRECISION :: apriori_std(ndyx)     ! a priori std 
END TYPE dyn_par

! default value when undefined
DOUBLE PRECISION, DIMENSION(4), PARAMETER :: zero_4d_vect = &
&    (/ 0.d0, 0.d0, 0.d0, 0.d0 /)

! undefined dynamical parameters
TYPE(dyn_par), PARAMETER :: undefined_dyn_par = DYN_PAR( &
& zero_4d_vect,                              & ! null parameters
&   0,                                       & ! no parameters
&   0,                                       & ! default no model
&  (/ .FALSE., .FALSE., .FALSE., .FALSE. /), & ! no active
&  (/ .FALSE., .FALSE., .FALSE., .FALSE. /), & ! no a priori
&  zero_4d_vect,                             & ! no value for a priori
&  zero_4d_vect                              & ! no std for a priori 
&     )   

PUBLIC dyn_par, undefined_dyn_par

! actual data

!options
LOGICAL, PUBLIC:: ngr_opt    ! if present, read options 
                             ! for non-gravitational perturbations 
                             ! from the option file

! parameters current value
TYPE(dyn_par), PUBLIC :: dyn
! list of solve for parameters to be currently determined (index in vector dyn)
! 0 means no parameter to be determined; note that number of indexes .ne.0 must not exceed ndp
! nls is the number of parameters actually being solved
INTEGER, PUBLIC :: ls(ndyx), nls 
! options for non gravitational perturbations: asteroids (dyn%nmod=1)
INTEGER,          PUBLIC :: iyark                   ! Yarkovsky model: 0=none, 1=diurnal, 2=seasonal, 3=sec acc 
INTEGER,          PUBLIC :: det_drp                 ! solve for: 0=none 1= direct, spherical 2=secular accel 3=both
INTEGER,          PUBLIC :: ipa2m                   ! direct radiation pressure activation
DOUBLE PRECISION, PUBLIC :: drpa2m                  ! A/M x CR coefficient of radiation pressure (m^2/ton)
DOUBLE PRECISION, PUBLIC :: A2                      ! secular perturbation due to nongrav (def=0)
! options for non gravitational perturbations: comets (dyn%nmod=2)
INTEGER,          PUBLIC :: ioutgas                 ! flag to select outgassing model, 0=no, 1=sym, 2=asym
INTEGER,          PUBLIC :: det_outgas              ! solve for: 
                                                    ! 0=none 1=a1ng,a2ng  2=a1ng,a2ng,a3ng 3=a1ng,a2ng,a3ng,dtdelay
DOUBLE PRECISION, PUBLIC :: a1ng,a2ng,a3ng,dtdelay  ! A1 term contribution (radial component)
                                                    ! A2 term contribution (transversal component)
                                                    ! A3 term contribution (orthogonal component)
                                                    ! Delta T delay for the asymmetric model
! options for non gravitational perturbations: SEP violation (dyn%nmod=3)
LOGICAL,          PUBLIC :: sep_viol                ! equivalence principle test
DOUBLE PRECISION, PUBLIC :: eta_sep                 ! amount of violation of SEP, eta (assuming Omega0=-3.52d-6)

! parameters control value
DOUBLE PRECISION, DIMENSION(2), PUBLIC, PARAMETER :: dpmin=(/-1.d0, -1.d0/), dpmax=(/2.d2, 2.d0/)

!public routines
PUBLIC out_ng, check_det

CONTAINS

 SUBROUTINE out_ng(iunng,debnac,indsol)
   USE name_rules
   INTEGER, INTENT(IN)       :: iunng  ! output unit
   INTEGER, INTENT(IN)       :: indsol ! index of solution for that name
   CHARACTER*(*), INTENT(IN) :: debnac ! object name
   INTEGER le
   CALL rmsp(debnac,le)
   WRITE(iunng,100)debnac(1:le),indsol,dyn%nmod,dyn%ndp,dyn%dp(1:dyn%ndp)
100 FORMAT(A,1X,I3,1X,I2,1X,I2,4(1X,F12.7))
 END SUBROUTINE out_ng

! Check that the determined parameters are active
 LOGICAL FUNCTION check_det(dyn,ls,nls)
! Begin interface
   TYPE(dyn_par), INTENT(IN) :: dyn
   INTEGER,       INTENT(IN) :: ls(ndyx), nls
! End interface
   INTEGER j
!========================================================
   ! check that no parameter is determined without being active
   check_det=.true.
   DO j=1,nls
      IF (dyn%active(ls(j)))THEN
         WRITE(*,*) 'determining parameter ', ls(j)
      ELSE 
         WRITE(*,*) 'determining inactive parameter ', ls(j)
         check_det=.false.
      ENDIF
   ENDDO
 END FUNCTION check_det


END MODULE dyn_param

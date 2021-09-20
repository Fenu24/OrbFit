!======================================================================!
! MODULE COBWEB                                                        !
!======================================================================!
! It contains the routines for the systematic ranging (Admissible      !
! Region sampling, MOV computation), for the computation of the        !
! probability density associated to the MOV and for the prediction     !
! of imminent impactors.                                               !
!                                                                      !
! Authors: A. Del Vigna, F. Spoto, G. Tommei, A. Milani                !
!                                                                      !
! Last Update: April 2020                                              !
!                                                                      !
! References:                                                          !
! [1] Spoto, Del Vigna (2018): "Short arc orbit determination and      !
!                               imminent impactors in the Gaia era"    !
!                                                                      !
! [2] Del Vigna (2018): "On Impact Monitoring of Near-Earth Asteroids" !
!                       (PhD Thesis)                                   !
!======================================================================!
MODULE cobweb

  USE output_control
  USE fund_const
  USE astrometric_observations
  USE orbit_elements
  USE least_squares
  USE name_rules

  IMPLICIT NONE

  PRIVATE

  TYPE mov_point
     TYPE(orbit_elem) :: el_mov     ! MOV orbit
     TYPE(orb_uncert) :: unc_mov    ! uncertainty matrices of the MOV orbit
     DOUBLE PRECISION :: chi_mov    ! Chi of the MOV orbit
     DOUBLE PRECISION :: rms_mov    ! RMS of the MOV orbit residuals
     LOGICAL          :: succ_mov   ! success flag for the MOV orbit computation
     DOUBLE PRECISION :: r_rdot(2)  ! (rho,rho_dot) of the MOV orbit
     DOUBLE PRECISION :: r_theta(2) ! (R,theta) of the MOV orbit (only spider web)
     DOUBLE PRECISION :: e_earth    ! two-body energy w.r.t. Earth
     DOUBLE PRECISION :: e_sun      ! two-body energy w.r.t. Sun
     DOUBLE PRECISION :: v_infty    ! velocity at infinity w.r.t. to Earth
     INTEGER          :: conn_comp  ! AR connected component of MOV orbit (0=outside AR)
     INTEGER          :: imp_flag   ! Flag for impacting orbit
     DOUBLE PRECISION :: tcla       ! time of closest approach of the MOV orbit
     DOUBLE PRECISION :: dcla       ! closest approach distance of the MOV orbit
     DOUBLE PRECISION :: csi        ! csi coordinate on the MTP
     DOUBLE PRECISION :: zeta       ! zeta coordinate on the MTP
     DOUBLE PRECISION :: moid       ! MOID
     DOUBLE PRECISION :: jac_sigma  ! Jacobian of the map from the sampling space to the AR
     DOUBLE PRECISION :: jac_mu     ! Jacobian of the map from the AR to the MOV
  END type mov_point

  !---------------------------------!
  ! Public parameters and variables !
  !---------------------------------!
  PUBLIC :: mov_point
  INTEGER, PARAMETER, PUBLIC :: nmax  = 250       ! Maximum number of points along one axis
  INTEGER, PARAMETER, PUBLIC :: nmax2 = nmax*nmax ! Maximum number of points in the grid
  !
  TYPE(mov_point),  PUBLIC :: mov_orbit(0:nmax,0:nmax) ! MOV orbits
  LOGICAL,          PUBLIC :: use_nominal              ! True when the nominal is used
  INTEGER,          PUBLIC :: i_chimin,j_chimin        ! Indices of MOV orbit with minimum chi
  DOUBLE PRECISION, PUBLIC :: sigx                     ! Max of chi for impact search
  DOUBLE PRECISION, PUBLIC :: hmax                     ! Shooting star limit
  INTEGER,          PUBLIC :: ndir                     ! Number of rho sample points
  INTEGER,          PUBLIC :: np                       ! Number of rho_dot sample points
  INTEGER,          PUBLIC :: ndir2                    ! Number of rho sample points (densified grid)
  INTEGER,          PUBLIC :: np2                      ! Number of rho-dot sample points (densified grid)
  LOGICAL,          PUBLIC :: grid_lev_curve           ! If true, a 50x50 grid is additionally computed in cobweb case
  LOGICAL,          PUBLIC :: propag_geoc_orbit        ! If true, the geocentric orbits are propagated for VI search
  LOGICAL,          PUBLIC :: mov_flag                 ! True if the MOV is available (used in fitobs)
  DOUBLE PRECISION, PUBLIC :: elong                    ! Elongation from the Sun (output of the AR computation)
  DOUBLE PRECISION, PUBLIC :: xogeo(3)                 ! Geocentric observer position (Equatorial)
  DOUBLE PRECISION, PUBLIC :: vogeo(3)                 ! Geocentric observer velocity (Equatorial)

  !---------------------------------!
  ! Common parameters and variables !
  !---------------------------------!
  DOUBLE PRECISION, PARAMETER :: rho_max_unif = 3.16227766016838d0 ! sqrt(10), max rho to switch to a unif. sampling
  DOUBLE PRECISION, PARAMETER :: chi_min_aux = 1.d8
  !
  DOUBLE PRECISION :: rdotmax_ar   ! Maximum value for rho_dot of the AR
  DOUBLE PRECISION :: rdotmin_ar   ! Minimum value for rho_dot of the AR
  DOUBLE PRECISION :: E_lim        ! Limit value for the energy
  DOUBLE PRECISION :: chi_min      ! Minimum chi among MOV orbits
  LOGICAL          :: significant  ! Significance flag (> 2 observations, arc > 0.5h and select_risk)
  LOGICAL          :: select_risk  ! True when an orbit with cov. is available and geoc_chi > 1
  LOGICAL          :: log_grid     ! Logarithmic grid in rho
  LOGICAL          :: good_first   ! Flag whith value TRUE if the first grid produces at least one point with chi<5

  !-----------------!
  ! Public routines !
  !-----------------!
  PUBLIC :: spider,mov,sample_web_grid,immediate_imp,mc_res_att

CONTAINS

  !============================================================!
  ! SAMPLE_WEB_GRID                                            !
  !============================================================!
  ! Sampling of admissible region, either by cobweb or by grid !
  !============================================================!
  SUBROUTINE sample_web_grid(obj_name,ini0,cov0,elc,unc,res_norm,m,obs,obsw,nd, &
       &                     arctype,geoc_chi,acce_chi,chi_curv)
    USE attributable
    USE arc_control
    USE triangles
    !=======================================================================================================
    CHARACTER*(name_len), INTENT(IN)    :: obj_name    ! Object name
    LOGICAL,              INTENT(IN)    :: ini0        ! Availability of orbital elements
    LOGICAL,              INTENT(IN)    :: cov0        ! Availability of covariance
    TYPE(orbit_elem),     INTENT(INOUT) :: elc         ! Nominal elements
    TYPE(orb_uncert),     INTENT(IN)    :: unc         ! Uncertainty matrices
    DOUBLE PRECISION,     INTENT(IN)    :: res_norm    ! Residuals norm
    INTEGER,              INTENT(IN)    :: m           ! Number of observations
    TYPE(ast_obs),        INTENT(IN)    :: obs(m)      ! Observations
    TYPE(ast_wbsr),       INTENT(INOUT) :: obsw(m)     ! Weights and residuals
    INTEGER,              INTENT(IN)    :: nd          ! Dimension of the orbital elements space
    INTEGER,              INTENT(OUT)   :: arctype    ! Arc type
    DOUBLE PRECISION,     INTENT(OUT)   :: geoc_chi    ! Normalised value of the geodesic curvature
    DOUBLE PRECISION,     INTENT(OUT)   :: acce_chi    ! Normalised value for the along-track acceleration
    DOUBLE PRECISION,     INTENT(OUT)   :: chi_curv    ! Chi-value of geodetic curvature and acceleration
    !===== Arc type computation ============================================================================
    TYPE(attrib)                :: attrc                 ! Attributable
    LOGICAL                     :: error                 ! Error in the attributable computation
    DOUBLE PRECISION            :: t_att_round           ! Rounded time for the attributable (MJD)
    INTEGER                     :: nigarc                ! Number of nights for the arc type
    INTEGER                     :: fail_arty             ! Arc type computation failure
    CHARACTER(LEN=3)            :: staz                  ! Observatory code (from attributable computation)
    DOUBLE PRECISION, PARAMETER :: sphx=2.d0             ! Maximum arc span (in degrees)
    !===== Admissible Region computation ===================================================================
    DOUBLE PRECISION, PARAMETER :: a_max=100.d0          ! Maximum value for the semimajor axis
    DOUBLE PRECISION            :: refs(3,3)             ! Matrix of the reference system
    DOUBLE PRECISION            :: E_bound               ! Bound for E_Sun
    DOUBLE PRECISION            :: c(0:5)                ! Coefficients
    INTEGER                     :: nroots_pol               ! Number of roots of the 6 deg. poly. (V(r), comments below)
    DOUBLE PRECISION            :: roots(3)              ! Roots of the polynomial V(r)
    INTEGER                     :: nroots_der               ! Number of roots of the pol. derivative V'(r)
    DOUBLE PRECISION            :: rootsd(8)             ! Roots of the pol. derivative V'(r)
    DOUBLE PRECISION            :: xo(3)                 ! Heliocentric observer position (Equatorial)
    DOUBLE PRECISION            :: vo(3)                 ! Heliocentric observer velocity (Equatorial)
    DOUBLE PRECISION            :: rmin                  ! Minimum value for rho
    INTEGER                     :: ind                   ! Index of a zero of V'(r) between rmin and the first pos. root
    DOUBLE PRECISION            :: rder0                 ! Zero of V'(r) between rmin and the first pos. root
    INTEGER                     :: ntot                  ! Total number of points for the boundary computation
    DOUBLE PRECISION            :: frv(npox)             ! Metric value on the boundary of the AR
    DOUBLE PRECISION            :: rhodot(npox)          ! Values of rho_dot in the boundary sampling
    DOUBLE PRECISION            :: rdw                   ! Scaling factor for range in rhodot
    !===== Cobweb/Grid computation =========================================================================
    DOUBLE PRECISION            :: rdotmax               ! Maximum value for rho_dot
    DOUBLE PRECISION            :: rdotmin               ! Minimum value for rho_dot
    DOUBLE PRECISION            :: rho                   ! Rho values in the grid
    DOUBLE PRECISION            :: rmax                  ! Maximum value for rho
    DOUBLE PRECISION            :: rms                   ! RMS for cobweb/grid computation
    DOUBLE PRECISION            :: jacthr(0:nmax,0:nmax) ! Jacobian of (R,theta) -> (rho,rho_dot)
    DOUBLE PRECISION            :: score_nea             ! Score point: NEA
    DOUBLE PRECISION            :: score_mba             ! Score point: MBA
    DOUBLE PRECISION            :: score_distant         ! Score point: distant
    DOUBLE PRECISION            :: score_scatt           ! Score point: scattered
    DOUBLE PRECISION            :: score_pha             ! Score point: PHA
    !===== File variables ==================================================================================
    CHARACTER(LEN=100)          :: file                  ! File name
    INTEGER                     :: le                    ! Lenght of file
    INTEGER                     :: iun_att                ! Unit for the att file
    INTEGER                     :: iunpol                ! Unit for the pol file (AR computation output)
    INTEGER                     :: iuncob                ! Unit for the web file (info on the sampling points)
    INTEGER                     :: iuncom                ! Unit for the com file (info on the sampling points)
    INTEGER                     :: iunele                ! Unit for the kepele file (KEP el. for the sampling points)
    INTEGER                     :: iunscore              ! Unit for the score file
    !===== Functions =======================================================================================
    DOUBLE PRECISION            :: prscal                ! Scalar product
    DOUBLE PRECISION            :: vsize                 ! Norm of a vector
    !===== Other variables =================================================================================
    DOUBLE PRECISION            :: crho(2,2)             ! Submatrix of C for (rho,rho_dot)
    DOUBLE PRECISION            :: grho(2,2)             ! Submatrix of G for (rho,rho_dot)
    DOUBLE PRECISION            :: att(4)                ! Attributable for the energy computation
    TYPE(orbit_elem)            :: elk                   ! Keplerian elements of the sampling points
    INTEGER                     :: fail_flag             ! Failure of the coordinate change
    INTEGER                     :: i, j, k, ii, jj, iii  ! Loop indexes
    DOUBLE PRECISION            :: dt                    ! Arc length
    LOGICAL                     :: curv_flag             ! Very significant curvature flag
    !=======================================================================================================
    IF(nd.NE.6)THEN
       WRITE(*,*) 'sample_web_grid: not ready for non-gravitational parameters, nd = ',nd
       STOP
    END IF
    ! Inizialization
    log_grid  = .FALSE.
    E_lim = -gms/(2.d0*a_max) ! energy value for the outer boundary of AR
    !**************************!
    ! Compute the attributable !
    !**************************!
    CALL tee(iun_log,'------------------------=')
    CALL tee(iun_log,'| COMPUTE ATTRIBUTABLE |=')
    CALL tee(iun_log,'------------------------=')
    CALL attri_comp(m,obs,obsw,attrc,error)
    IF(error)THEN
       WRITE(*,*) 'sample_web_grid: error during the attributable computation'
       STOP
    END IF
    ! Write the attributable computation output
    t_att_round = NINT(attrc%tdtobs)
    IF(attrc%sph*radeg.GT.sphx)THEN
       WRITE(*,*) 'sample_web_grid: too wide arc, arc span [deg] = ', attrc%sph*radeg
    ELSE
       CALL wri_attri(0,0,obj_name,attrc,t_att_round)
       CALL wri_attri(iun_log,iun_log,obj_name,attrc,t_att_round)
    ENDIF
    !**********************!
    ! Compute the arc type !
    !**********************!
    CALL tee(iun_log,'--------------------=')
    CALL tee(iun_log,'| COMPUTE ARC TYPE |=')
    CALL tee(iun_log,'--------------------=')
    arctype = arc_type(obs,obsw,m,geoc_chi,acce_chi,chi_curv,nigarc,fail_arty)
    ! Write the arc type
    WRITE(*,176) arctype
    WRITE(*,177) m,nigarc,fail_arty,geoc_chi,acce_chi,chi_curv
    WRITE(iun_log,176) arctype
    WRITE(iun_log,177) m,nigarc,fail_arty,geoc_chi,acce_chi,chi_curv
176 FORMAT(' The arc type is ',I3)
177 FORMAT(' nobs = ',I4,', nights = ',I3,', errcode = ',I4/ &
         & ' geoc_chi = ',1P,D9.2,', acce_chi = ',D9.2,', chi_curv = ',D9.2)
    !*******************************!
    ! Compute the Admissible Region !
    !*******************************!
    !------------------------------------------------------------------!
    ! The most important condition defining the AR is the              !
    ! non-positivity of the two body energy of the body w.r.t the      !
    ! Sun. This condition gives solutions for r_dot if and only if the !
    ! following inequality holds:                                      !
    !                                                                  !
    !                          V(r) <= 4k^4                            !
    !                                                                  !
    ! where V(r) is a 6-degree polynomial. We can have only three      !
    ! possibilities:                                                   !
    !                                                                  !
    ! (1) V(r) has 4 simple real roots -> the AR has 2 connected       !
    !     components;                                                  !
    !                                                                  !
    ! (2) V(r) has 3 distinct real roots, 2 simple and 1 with even     !
    !     multiplicity -> the second c.c. reduces to a point;          !
    !                                                                  !
    ! (3) V(r) has 2 distinct real roots, 1 simple and 1 with odd      !
    !     multiplicity -> the AR has 1 c.c.                            !
    !                                                                  !
    ! Moreover, V'(r) cannot have more than 3 distinct roots. If they  !
    ! are exactly three, then there cannot be any root with            !
    ! multiplicity 2.                                                  !
    !------------------------------------------------------------------!
    !--------------------------------!
    !  Output of the AR computation  !
    !--------------------------------!
    file = obj_name//'.att'
    CALL rmsp(file,le)
    CALL filopn(iun_att,file(1:le),'unknown')
    staz = attrc%obscod
    CALL tee(iun_log,'-----------------------------=')
    CALL tee(iun_log,'| COMPUTE ADMISSIBLE REGION |=')
    CALL tee(iun_log,'-----------------------------=')
    CALL admis_reg(obj_name,iun_att,attrc%tdtobs,attrc%angles,staz,nroots_pol,roots,c, &
         &         refs,xo,vo,a_max,E_bound,nroots_der,rootsd,attrc%apm,xogeo,vogeo,elong)
    CALL filclo(iun_att,' ')
    ! Output the roots of the polynomial V(r)
    WRITE(*,*) 'Roots = ', roots(1:nroots_pol)
    WRITE(iun_log,*) 'Roots = ', roots(1:nroots_pol)
    !----------------------------------------!
    !  Find the range of values for rho_dot  !
    !----------------------------------------!
    ! Arbitrary lower boundary for the AR
    rmin = 1.d-4
    ! Search for a zero (rder0) of the derivative V'(r) between rmin and roots(1)
    rder0 = -1.d0
    ind   = 0
    IF(nroots_der.GE.1)THEN
       IF(nroots_der.EQ.1)THEN
          ! Do nothing
       ELSEIF(nroots_der.EQ.2)THEN
          ! Do nothing
       ELSEIF(nroots_der.EQ.3)THEN
          DO iii=1,nroots_der
             IF(rootsd(iii).GT.rmin .AND. rootsd(iii).LT.roots(1))THEN
                ind=iii
             ENDIF
          ENDDO
          IF(ind.GT.2)THEN
             IF(rootsd(ind-2).GT.rmin)THEN
                rder0=rootsd(ind-2)
                WRITE(*,*) 'rmin = ',rmin,', rder0 = ',rder0,', roots(1) = ',roots(1)
             ENDIF
          ENDIF
          IF(rder0.GE.roots(1))THEN
             WRITE(*,*) 'sample_web_grid: ERROR! rder0 > roots(1)'
          ENDIF
       ELSE
          ! V'(r) cannot have more than 3 distinct roots
          WRITE(*,*) 'sample_web_grid: number of roots of the derivative of V(r) = ', nroots_der
          WRITE(*,*) '                 It cannot be possible!'
       ENDIF
    ENDIF
    ! Set AR boundary sampling metric (public variable)
    nfunc = 2 ! Logaritmic scale
    CALL sample_bound_ne1(E_bound,attrc%eta**2,rmin,rder0,roots(1),8,5, &
         &                ntot,c,frv,rhodot,npox,rdw)
    rdotmin = MINVAL(rhodot(1:ntot))
    rdotmax = MAXVAL(rhodot(1:ntot))
    rdotmin_ar = rdotmin
    rdotmax_ar = rdotmax

    !**********************************!
    ! Set flags for systematic ranging !
    !**********************************!
    !------------------------------------------------------------------!
    ! The sampling of the AR depends on the availability of a nominal  !
    ! solution. If there is a nominal solution, with its covariance,   !
    ! and if the SNR of the geodetic curvature is > 3, then we trust   !
    ! the nominal solution and we do a cobweb sampling.                !
    !------------------------------------------------------------------!
    use_nominal = (ini0 .AND. cov0 .AND. ABS(geoc_chi).GT.3.d0) ! reliability of nominal solution
    !-------------------------------------------------------------!
    ! We consider significant a case that has more than 3         !
    ! observations, and an arc of more than 0.5 hours (i.e., 30   !
    ! minutes), or a case which has been selected through the     !
    ! select_risk flag, which is true when geoc_chi > 1.          !
    !-------------------------------------------------------------!
    dt          = (obs(m)%time_tdt-obs(1)%time_tdt)*24.d0       ! Arc duration (hours)
    select_risk = (ini0 .AND. cov0 .AND. ABS(geoc_chi).GT.1.d0) ! Significance despite m and arc length
    significant = (m.GT.2 .AND. dt.GT.0.5d0) .OR. (select_risk) ! Significant flag
    !---------------------------------------------------------------!
    ! If the orbit quality is high geocentric orbits are propagated !
    ! for the impact search (which is false by default).            !
    !---------------------------------------------------------------!
    propag_geoc_orbit = .FALSE.
    IF(use_nominal)THEN
       mov_orbit(0,0)%r_rdot(1:2) = elc%coord(5:6)
       CALL energy_earth(attrc%angles,xogeo,vogeo,.FALSE.)
       curv_flag         = chi_curv**2.GT.10 .AND. ABS(geoc_chi).GT.5.d0
       propag_geoc_orbit = (mov_orbit(0,0)%e_earth.LT.0.d0 .AND. (curv_flag .OR. arctype.GT.2))
    END IF
    !*************************!
    ! Cobweb/Grid computation !
    !*************************!
    DO k=1,2
       IF(k.EQ.1)THEN
          ! Skip the first grid computation: spider web case and when grid_lev_curve is false
          IF(use_nominal .AND. .NOT.grid_lev_curve) CYCLE
          !-----------------------------------------------------------------!
          !                  First iteration (first grid)                   !
          !-----------------------------------------------------------------!
          good_first = .TRUE.
          ! Selection of logarithmic grid in rho if rmax is too big (only with one component)
          IF(roots(1).LT.rho_max_unif) log_grid = .TRUE.
          ! Selection of other grid parameters
          CALL select_grid_param(k,nroots_pol,roots(1:3),rmin,rmax,rdotmin,rdotmax)
          !--------------------!
          ! Systematic ranging !
          !--------------------!
          CALL syst_ranging(obj_name,m,obs,obsw,nd,attrc,elc,unc,res_norm,c(0:5),nroots_pol,roots(1:3), &
               &            rmin,rmax,rdotmin,rdotmax,k)
          IF(.NOT.use_nominal)THEN
             ! In the grid case, compute NEA score to use for the second iteration
             CALL compute_score(score_distant,score_mba,score_nea,score_scatt,score_pha)
          END IF
       ELSEIF(k.EQ.2)THEN
          IF(.NOT.use_nominal)THEN
             !-------------------------------------------------------------------!
             !                Second iteration (densified grid)                  !
             !-------------------------------------------------------------------!
             !-------------------------------------------------------------!
             ! Once we obtain a first preliminary grid, we densify for a   !
             ! higher resolution. We select the minimum and the maximum    !
             ! value of rho and rho_dot among all the values of the        !
             ! points for which we had convergence with chi < 5.           !
             !-------------------------------------------------------------!
             IF(nroots_pol.EQ.1)THEN
                log_grid = .FALSE.
                IF(score_nea.GT.0.5d0) log_grid = .TRUE.
             END IF
             ! Selection of other grid parameters
             CALL select_grid_param(k,nroots_pol,roots(1:3),rmin,rmax,rdotmin,rdotmax)
             IF(.NOT.good_first)THEN
                ! No good point in the first iteration
                CALL tee(iun_log,'sample_web_grid: No good points in the previous iteration=')
                CALL write_arsample_file(obj_name,attrc%angles,c(0:5),roots(1:3),k)
                CALL write_movsample_file(obj_name,attrc%tdtobs,k)
                EXIT
             END IF
             !--------------------!
             ! Systematic ranging !
             !--------------------!
             CALL syst_ranging(obj_name,m,obs,obsw,nd,attrc,elc,unc,res_norm,c(0:5),nroots_pol,roots(1:3), &
                  &            rmin,rmax,rdotmin,rdotmax,k)
          ELSE
             !-------------------------------------------------------------------!
             !                    Second iteration (cobweb)                      !
             !-------------------------------------------------------------------!
             ! Set number of directions and number of points per direction
             IF(propag_geoc_orbit)THEN
                CALL tee(iun_log,'Automatic choice for the propagation of geocentric orbits=')
                ndir       = ndir2/4
                np         = np2/4
             ELSE
                ndir       = ndir2
                np         = np2
             END IF
             good_first = .TRUE.  ! No first iteration
             CALL syst_ranging(obj_name,m,obs,obsw,nd,attrc,elc,unc,res_norm,c(0:5),nroots_pol,roots(1:3), &
                  &            rmin,rmax,rdotmin,rdotmax,k)
          END IF
       END IF
    END DO
    !***************!
    ! Compute score !
    !***************!
    CALL compute_score(score_distant,score_mba,score_nea,score_scatt,score_pha)
    ! Open output file
    file = obj_name//'.score'
    CALL rmsp(file,le)
    CALL filopn(iunscore,file(1:le),'unknown')
    ! Write header score file
    WRITE(iunscore,'(A)') '% Object name = '//obj_name
    WRITE(iunscore,'(A)') '%'
    WRITE(iunscore, 595) 'Score Near-Earth = ', score_nea
    WRITE(iunscore, 595) 'Score Main Belt  = ', score_mba
    WRITE(iunscore, 595) 'Score Distant    = ', score_distant
    WRITE(iunscore, 595) 'Score Scattered  = ', score_scatt
    WRITE(iunscore, 595) 'Score PHA        = ', score_pha
595 FORMAT(A, 3X, F4.2)
    CALL filclo(iunscore,' ')

  END SUBROUTINE sample_web_grid


  !============================================================!
  ! SELECT_GRID_PARAM                                          !
  !============================================================!
  ! Selection of grid end-points and other parameters.         !
  !============================================================!
  SUBROUTINE select_grid_param(k_ind,n_roots,roots,rmin,rmax,rdotmin,rdotmax)
    !=======================================================================================================
    INTEGER,          INTENT(IN)     :: k_ind     ! Iteration index
    INTEGER,          INTENT(IN)     :: n_roots   ! Number of roots of the 6 deg. poly.
    DOUBLE PRECISION, INTENT(IN)     :: roots(3)  ! Roots of the polynomial defining the AR
    DOUBLE PRECISION, INTENT(INOUT)  :: rmin      ! Starting point for the rho sampling
    DOUBLE PRECISION, INTENT(INOUT)  :: rmax      ! Ending point for the rho sampling
    DOUBLE PRECISION, INTENT(INOUT)  :: rdotmin   ! Starting point for the rhodot sampling
    DOUBLE PRECISION, INTENT(INOUT)  :: rdotmax   ! Ending point for the rhodot sampling
    !=======================================================================================================
    DOUBLE PRECISION, PARAMETER :: rmin_aux = 1000.d0     ! Auxiliary rmin
    DOUBLE PRECISION, PARAMETER :: rmax_aux = 0.d0        ! Auxiliary rmax
    DOUBLE PRECISION, PARAMETER :: rdotmin_aux = 1000.d0  ! Auxiliary rdotmin
    DOUBLE PRECISION, PARAMETER :: rdotmax_aux = -1000.d0 ! Auxiliary rdotmax
    LOGICAL                     :: same_side              ! TRUE when rmin and rmax are on the same side w.r.t. r_1
    INTEGER                     :: i_min,i_max            ! Index of min and max rho
    INTEGER                     :: j_min,j_max            ! Index of min and max rho_dot
    INTEGER                     :: i,j                    ! Loop indexes
    INTEGER                     :: ndir_first,np_first    ! Auxiliary grid dimensions
    !=======================================================================================================
    IF(k_ind.EQ.1)THEN
       !----------------------------------------------------------!
       ! First iteration: rho selected from the AR, rho_dot as is !
       !----------------------------------------------------------!
       IF(n_roots.EQ.3)THEN
          ! Two connected components
          rmin = 1.d-3
          rmax = roots(3)
       ELSE
          ! One connected component
          rmax = roots(1)
       END IF
       CALL tee(iun_log,'----------------=')
       CALL tee(iun_log,'| COMPUTE GRID |=')
       CALL tee(iun_log,'----------------=')
    ELSE
       ! Store first iteration grid dimension
       ndir_first = ndir
       np_first   = np
       !----------------!
       ! Densified grid !
       !----------------!
       ndir = ndir2
       np   = np2
       !---------------------------------------------------------------------------!
       ! Second iteration: search for successful points to increase the resolution !
       !---------------------------------------------------------------------------!
       ! Initialise AR sampling boundaries (variables *_aux defined above as parameters)
       rmin    = rmin_aux
       rmax    = rmax_aux
       rdotmin = rdotmin_aux
       rdotmax = rdotmax_aux
       DO i=1,ndir_first
          DO j=1,np_first
             IF(mov_orbit(i,j)%chi_mov.LT.5.d0 .AND. mov_orbit(i,j)%el_mov%h_mag.LT.hmax .AND. mov_orbit(i,j)%succ_mov)THEN
                IF(mov_orbit(i,j)%r_rdot(1).LT.rmin)THEN
                   rmin  = mov_orbit(i,j)%r_rdot(1)
                   i_min = i
                END IF
                IF(mov_orbit(i,j)%r_rdot(1).GT.rmax)THEN
                   rmax  = mov_orbit(i,j)%r_rdot(1)
                   i_max = i
                END IF
                IF(mov_orbit(i,j)%r_rdot(2).LT.rdotmin)THEN
                   rdotmin = mov_orbit(i,j)%r_rdot(2)
                   j_min   = j
                END IF
                IF(mov_orbit(i,j)%r_rdot(2).GT.rdotmax)THEN
                   rdotmax = mov_orbit(i,j)%r_rdot(2)
                   j_max   = j
                END IF
             END IF
          END DO
       END DO
       IF(rmin.EQ.rmin_aux .OR. rmax.EQ.rmax_aux .OR. rdotmin.EQ.rdotmin_aux .OR. rdotmax.EQ.rdotmax_aux)THEN
          IF(chi_min .EQ. chi_min_aux)THEN
             ! No successful points in the first iteration
             good_first = .FALSE.
             RETURN
          ELSE
             i_min = i_chimin
             i_max = i_chimin
             j_min = j_chimin
             j_max = j_chimin
             ! No need of a densified grid in a grid around a single point
             ndir = ndir_first
             np   = np_first
          END IF
       END IF
       !--------------------------------------------------------------!
       ! To avoid loss of good points due to the lower resolution of  !
       ! the first grid, enlarge to the nearby columns                !
       !--------------------------------------------------------------!
       IF(i_min.EQ.1) i_min = 2
       IF(j_min.EQ.1) j_min = 2
       IF(i_max.EQ.ndir_first) i_max = ndir_first-1
       IF(j_max.EQ.np_first)   j_max = np_first-1
       rmin    = mov_orbit(i_min-1,1)%r_rdot(1)
       rmax    = mov_orbit(i_max+1,1)%r_rdot(1)
       rdotmin = mov_orbit(1,j_min-1)%r_rdot(2)
       rdotmax = mov_orbit(1,j_max+1)%r_rdot(2)
       CALL tee(iun_log,'--------------------------=')
       CALL tee(iun_log,'| COMPUTE DENSIFIED GRID |=')
       CALL tee(iun_log,'--------------------------=')
    END IF
    !----------------------------------------------!
    ! Set number of sample points in rho if needed !
    !----------------------------------------------!
    ! Position of rmin and rmax w.r.t. r_1
    same_side = same_side_r1(rmin,rmax,roots(1))
    IF(n_roots.EQ.3 .AND. .NOT.same_side)THEN
       ! Two connected components to cover: double the resolution in rho
       ndir = 2*ndir
    END IF
    !-----------------!
    ! Output messages !
    !-----------------!
    WRITE(*,*) '   Grid kind selection: logarithmic grid = ', log_grid
    WRITE(*,101) ndir,np
    IF(same_side)THEN
       WRITE(*,102) rmin,rmax
    ELSE
       WRITE(*,103) rmin,roots(1),roots(2),rmax
    END IF
    WRITE(*,104) rdotmin,rdotmax
101 FORMAT('    Grid dimension = ',I4,' x 'I4)
102 FORMAT('    Sampling interval in rho = [',F12.6,','F12.6,']')
103 FORMAT('    Sampling interval in rho = [',F12.6,','F12.6,'] U [',F12.6,','F12.6,']')
104 FORMAT('    Sampling interval in rho_dot = [',F12.6,','F12.6,']')

  END SUBROUTINE select_grid_param


  !====================================================================!
  ! SYST_RANGING                                                       !
  !====================================================================!
  ! It performs the systematic ranging on the AR: it samples a grid or !
  ! a spider web, compute the Manifold Of Variations and write the     !
  ! output files.                                                      !
  !====================================================================!
  SUBROUTINE syst_ranging(obj_name,m,obs,obsw,nd,attrc,el_nom,unc_nom,res_norm,coeff,n_roots,roots, &
       &                  rmin,rmax,rhodot_min,rhodot_max,k_ind)
    USE attributable
    USE station_coordinates, ONLY: statcode
    !=======================================================================================================
    CHARACTER*(name_len), INTENT(IN)    :: obj_name      ! Asteroid name
    INTEGER,              INTENT(IN)    :: m           ! Number of observations
    TYPE(ast_obs),        INTENT(IN)    :: obs(m)      ! Observations
    TYPE(ast_wbsr),       INTENT(INOUT) :: obsw(m)     ! Weights, residuals
    INTEGER,              INTENT(IN)    :: nd          ! Dimension of the parameters space
    TYPE(attrib),         INTENT(IN)    :: attrc       ! Attributable
    TYPE(orbit_elem),     INTENT(INOUT) :: el_nom      ! Nominal elements
    TYPE(orb_uncert),     INTENT(IN)    :: unc_nom     ! Uncertainty of the nominal elements
    DOUBLE PRECISION,     INTENT(IN)    :: res_norm    ! Residuals norm
    DOUBLE PRECISION,     INTENT(IN)    :: coeff(0:5)  ! Coefficients for the AR computation
    INTEGER,              INTENT(IN)    :: n_roots     ! Number of roots of the 6 deg. poly.
    DOUBLE PRECISION,     INTENT(IN)    :: roots(3)    ! Roots of the polynomial defining the AR
    DOUBLE PRECISION,     INTENT(IN)    :: rmin        ! Minimum value for rho
    DOUBLE PRECISION,     INTENT(IN)    :: rmax        ! Maximum value for rho
    DOUBLE PRECISION,     INTENT(IN)    :: rhodot_max  ! Maximum value for rho_dot
    DOUBLE PRECISION,     INTENT(IN)    :: rhodot_min  ! Minimum value for rho_dot
    INTEGER,              INTENT(IN)    :: k_ind       ! Loop index, for output
    !=======================================================================================================
    DOUBLE PRECISION  :: res_norm_aux                   ! RMS for cobweb/grid computation
    DOUBLE PRECISION  :: jacthr(0:nmax,0:nmax) ! Jacobian of (R,theta) -> (rho,rho_dot)
    DOUBLE PRECISION  :: att(4)                ! Attributable for the energy computation
    INTEGER           :: i,j
    !=======================================================================================================
    ! Initialization
    mov_orbit%el_mov%h_mag = 0.d0
    mov_orbit%chi_mov      = 0.d0
    mov_orbit%succ_mov     = .FALSE.
    mov_orbit%conn_comp    = 0
    mov_orbit%imp_flag     = 0
    mov_orbit%r_rdot(1)    = 0.d0
    mov_orbit%r_rdot(2)    = 0.d0
    mov_orbit%moid         = 0.d0
    mov_orbit%v_infty      = 0.d0
    mov_orbit%jac_sigma    = 0.d0
    mov_orbit%jac_mu       = 0.d0
    ! Computation of the AR sampling
    IF(use_nominal .AND. k_ind.EQ.2)THEN
       !-----------------------------!
       ! Spider web case (100 x 100) !
       !-----------------------------!
       IF(el_nom%coo.NE.'ATT')THEN
          WRITE(*,*) ' syst_ranging: nominal orbit must be in ATT. It is in ', el_nom%coo
          STOP
       ENDIF
       CALL tee(iun_log,'------------------=')
       CALL tee(iun_log,'| COMPUTE COBWEB |=')
       CALL tee(iun_log,'------------------=')
       ! Residual norm
       res_norm_aux = res_norm
       ! Computation of spider web
       CALL spider(el_nom,unc_nom)
    ELSE
       !-----------!
       ! Grid case !
       !-----------!
       ! Residual norm
       res_norm_aux = -1.d0
       ! Grid sampling
       CALL sample_grid(rmin,rmax,rhodot_min,rhodot_max,n_roots,roots)
    END IF
    !-----------------------!
    ! Energy w.r.t. the Sun !
    !-----------------------!
    CALL energy_sun(coeff,roots,rmin)
    !************************!
    ! Manifold Of Variations !
    !************************!
    ! If a reliable nominal solution does not exist, use only the computed attributable
    IF(.NOT.use_nominal)THEN
       !el_nom     = undefined_orbit_elem
       el_nom%coo = 'ATT'                          ! ATT coordinates
       el_nom%t   = attrc%tdtobs                   ! Time of the attributable
       CALL statcode(attrc%obscod,el_nom%obscode)  ! Observatory code
       el_nom%coord(1:4) = attrc%angles            ! Attributable for the first 4 coordinates
       el_nom%coord(5:6) = (/0.d0, 0.d0/)          ! Coordinates (0,0) for (rho,rho_dot)
    END IF
    ! Computation of the Manifold Of Variations
    CALL mov(m,obs,obsw,el_nom,unc_nom,res_norm_aux,nd,k_ind)
    mov_flag = .TRUE.
    ! Energy w.r.t. the Earth
    att = attrc%angles ! Attributable (angles and derivatives)
    CALL energy_earth(att,xogeo,vogeo,.TRUE.)
    ! Write output
    CALL write_arsample_file(obj_name,att,coeff,roots,k_ind)
    IF(k_ind.EQ.1)THEN
       ! At first iteration the MOV sample can be written here because we do not search for impactors
       CALL write_movsample_file(obj_name,attrc%tdtobs,k_ind)
    END IF

  END SUBROUTINE syst_ranging


  !============================================================!
  ! SAMPLE_GRID                                                !
  !============================================================!
  ! Grid in the (rho,rho_dot) space.                           !
  !============================================================!
  SUBROUTINE sample_grid(rmin,rmax,rdotmin,rdotmax,n_roots,roots)
    !=======================================================================================================
    DOUBLE PRECISION, INTENT(IN) :: rmin      ! Starting point for the rho sampling
    DOUBLE PRECISION, INTENT(IN) :: rmax      ! Ending point for the rho sampling
    DOUBLE PRECISION, INTENT(IN) :: rdotmin   ! Starting point for the rhodot sampling
    DOUBLE PRECISION, INTENT(IN) :: rdotmax   ! Ending point for the rhodot sampling
    INTEGER,          INTENT(IN) :: n_roots   ! Number of roots of the 6 deg. poly.
    DOUBLE PRECISION, INTENT(IN) :: roots(3)  ! Roots of the polynomial defining the AR
    !=======================================================================================================
    DOUBLE PRECISION :: dr        ! Sampling step in rho
    DOUBLE PRECISION :: dr1       ! Sampling step in rho (first c.c.)
    DOUBLE PRECISION :: dr2       ! Sampling step in rho (second c.c.)
    DOUBLE PRECISION :: drdot     ! Sampling step in rho_dot
    LOGICAL          :: same_side ! TRUE when rmin and rmax are on the same side w.r.t. r_1
    INTEGER          :: i,j       ! Loop indexes
    INTEGER          :: ndir_mid  ! Half of ndir
    !=======================================================================================================
    ! Establish how many connected components of the AR are to sample
    same_side =  same_side_r1(rmin,rmax,roots(1))
    IF(same_side)THEN
       !***********************************!
       ! One connected component to sample !
       !***********************************!
       IF(log_grid)THEN
          dr = (LOG10(rmax) - LOG10(rmin))/(ndir - 1)
       ELSE
          dr = (rmax - rmin)/(ndir - 1)
       END IF
       !-----------------!
       ! Sampling of rho !
       !-----------------!
       DO i=1,ndir
          IF(log_grid)THEN
             mov_orbit(i,1:np)%r_rdot(1) = 10**(dr*(i-1) + LOG10(rmin))
          ELSE
             mov_orbit(i,1:np)%r_rdot(1) = dr*(i-1) + rmin
          END IF
       ENDDO
    ELSE
       !**************************!
       ! Two connected components !
       !**************************!
       IF(MOD(ndir,2).NE.0)THEN
          CALL tee(iun_log,' The value of ndir must be even!=')
          STOP
       END IF
       ndir_mid = ndir/2
       IF(log_grid)THEN
          dr1 = (LOG10(roots(1)) - LOG10(rmin))/(ndir_mid - 1)
          dr2 = (LOG10(rmax) - LOG10(roots(2)))/(ndir_mid - 1)
       ELSE
          dr1 = (roots(1) - rmin)/(ndir_mid - 1)
          dr2 = (rmax - roots(2))/(ndir_mid - 1)
       END IF
       !-----------------!
       ! Sampling of rho !
       !-----------------!
       DO i=1,ndir_mid
          IF(log_grid)THEN
             mov_orbit(i,1:np)%r_rdot(1)          = 10**(dr1*(i-1) + LOG10(rmin))
             mov_orbit(ndir_mid+i,1:np)%r_rdot(1) = 10**(dr2*(i-1) + LOG10(roots(2)))
          ELSE
             mov_orbit(i,1:np)%r_rdot(1)          = dr1*(i-1) + rmin
             mov_orbit(ndir_mid+i,1:np)%r_rdot(1) = dr2*(i-1) + roots(2)
          END IF
       END DO
    END IF
    !---------------------!
    ! Sampling of rho_dot !
    !---------------------!
    drdot = (rdotmax - rdotmin)/(np-1)
    DO j=1,np
        mov_orbit(1:ndir,j)%r_rdot(2) = rdotmin + drdot*(j-1)
    END DO
    !--------------------------------------------!
    ! Determinant of the sampling map (cfr. [1]) !
    !--------------------------------------------!
    IF(log_grid)THEN
       mov_orbit(1:ndir,1:np)%jac_sigma = LOG(10.d0)*mov_orbit(1:ndir,1:np)%r_rdot(1)
    ELSE
       mov_orbit(1:ndir,1:np)%jac_sigma = 1.d0
    END IF

  END SUBROUTINE sample_grid


  !====================================================================!
  ! SPIDER                                                             !
  !====================================================================!
  ! Computation of the sample points in the (R,theta) space of polar   !
  ! elliptic coordinates.                                              !
  !====================================================================!
  SUBROUTINE spider(el_nom,unc)
    !=======================================================================================================
    TYPE(orbit_elem), INTENT(IN)  :: el_nom    ! Nominal elements (ATT coord.)
    TYPE(orb_uncert), INTENT(IN)  :: unc       ! Uncertainty of el_nom
    !=======================================================================================================
    DOUBLE PRECISION :: gamma_rho(2,2) ! Covariance matrix for (rho,rho_dot)
    DOUBLE PRECISION :: eigenval(2)    ! Eigenvalues of gamma_rho
    DOUBLE PRECISION :: eigenvec(2,2)  ! Matrix of eigenvectors of gamma_rho
    DOUBLE PRECISION :: fv1(2),fv2(2)  ! Temporary storage arrays for rs
    DOUBLE PRECISION :: v1(2),v2(2)    ! Eigenvectors of gamma_rho
    DOUBLE PRECISION :: a,b            ! Ellipse's semiaxis
    INTEGER          :: ierr           ! Error flag
    INTEGER          :: i,j            ! Loop indexes
    !=======================================================================================================
    mov_orbit(0,0)%r_rdot(1:2) = el_nom%coord(5:6)
    ! Covariance matrix for (rho,rho_dot)
    gamma_rho(1,1) = unc%g(5,5)
    gamma_rho(1,2) = unc%g(5,6)
    gamma_rho(2,1) = unc%g(6,5)
    gamma_rho(2,2) = unc%g(6,6)
    CALL rs(2,2,gamma_rho,eigenval,1,eigenvec,fv1,fv2,ierr)
    IF(ierr.NE.0)THEN
       WRITE(*,*) 'spider: failed computation of eigenvalues and eigenvectors with flag = ',ierr
       STOP
    END IF
    ! Eigenvectors and eigenvalues of gamma_rho
    v1 = eigenvec(:,1)     ! Longest axis direction
    v2 = eigenvec(:,2)     ! Shortest axis direction
    a  = SQRT(eigenval(1)) ! Semi-major axis
    b  = SQRT(eigenval(2)) ! Semi-minor axis
    !***************************************************************!
    !  Regular grid in (R,theta) plane: elliptic polar coordinates  !
    !***************************************************************!
    !------------------------------------------------------------------!
    ! Create a regular grid of points in the space of polar elliptic   !
    ! coordinates (R,theta), where R is in [0,sigx] (sigx is a         !
    ! parameter defining the maximum value of chi we consider reliable !
    ! in our analysis), and theta is in [0,2*pi). sigx is currently    !
    ! 5.d0, and it is set into fitobs.def.                             !
    !------------------------------------------------------------------!
    ! Map from (R,theta) to (rho,rho_dot)
    DO i=1,ndir
       DO j=1,np
          ! Fill (R, theta)
          mov_orbit(i,j)%r_theta(1) = (sigx/np)*j
          mov_orbit(i,j)%r_theta(2) = (dpig/ndir)*i
          ! Fill (rho, rho_dot)
          mov_orbit(i,j)%r_rdot(1) = mov_orbit(i,j)%r_theta(1)*(a*COS(mov_orbit(i,j)%r_theta(2))*v1(1) - &
               &                     b*SIN(mov_orbit(i,j)%r_theta(2))*v1(2)) + el_nom%coord(5)
          mov_orbit(i,j)%r_rdot(2) = mov_orbit(i,j)%r_theta(1)*(a*COS(mov_orbit(i,j)%r_theta(2))*v1(2) + &
               &                     b*SIN(mov_orbit(i,j)%r_theta(2))*v1(1)) + el_nom%coord(6)
          !--------------------------------------------!
          ! Determinant of the sampling map (cfr. [1]) !
          !--------------------------------------------!
          mov_orbit(i,j)%jac_sigma = mov_orbit(i,j)%r_theta(1)*a*b
       END DO
    END DO

  END SUBROUTINE spider


  !====================================================================!
  ! MOV                                                                !
  !====================================================================!
  ! Computation of the Manifold Of Variations.                         !
  !====================================================================!
  SUBROUTINE mov(m,obs,obsw,el_nom,unc_nom,rms,nd,k_ind,unc4_mov)
    INCLUDE 'parobx.h90'
    !=======================================================================================================
    INTEGER,          INTENT(IN)            :: m                        ! Number of observations
    TYPE(ast_obs),    INTENT(IN)            :: obs(m)                   ! Observations
    TYPE(ast_wbsr),   INTENT(INOUT)         :: obsw(m)                  ! Weights and observations
    TYPE(orbit_elem), INTENT(IN)            :: el_nom                   ! Nominal elements
    TYPE(orb_uncert), INTENT(IN)            :: unc_nom                  ! Uncertainty of the nominal orbit
    DOUBLE PRECISION, INTENT(IN)            :: rms                      ! Residuals norm
    INTEGER,          INTENT(IN)            :: nd                       ! Dimension of the parameters space
    INTEGER,          INTENT(IN)            :: k_ind                    ! Iteration index
    TYPE(orb_uncert), INTENT(OUT), OPTIONAL :: unc4_mov(0:nmax,0:nmax)  ! Uncertainty mat. of the 4-dim fit
    !===== Computation of the Jacobian =====================================================================
    DOUBLE PRECISION, PARAMETER :: jac_mu_min=1.d-4 ! Minimum value for jac_mu
    TYPE(orb_uncert) :: unc4_mov_loc(0:nmax,0:nmax) ! Uncertainty mat. of the 4-dim fit
    INTEGER          :: scal_obs                    ! Number of scalar observations
    DOUBLE PRECISION :: des_mat(nob2x,nd)           ! Design matrix
    DOUBLE PRECISION :: b_rho(nob2x,2)              ! Submatrix of B relative to (rho,rho_dot)
    DOUBLE PRECISION :: b_A(nob2x,4)                ! Submatrix of B relative to A
    DOUBLE PRECISION :: b_A_rho(4,2)                ! Product of b_A and b_rho
    DOUBLE PRECISION :: dAdrho(4,2)                 ! Derivative of the corrected attributable w.r.t. (rho,rho_dot)
    DOUBLE PRECISION :: M_A(2,2)                    ! Matrix of the map AR -> MOV
    !=======================================================================================================
    TYPE(orbit_elem) :: eltemp(0:nmax,0:nmax) ! Temporary orbit elements
    INTEGER          :: nused(0:nmax,0:nmax)  ! Number of observations used in the differential corrections
    DOUBLE PRECISION :: rmsh(0:nmax,0:nmax)   ! Weighted RMS of the residuals
    INTEGER          :: i,j,ii                ! Loop indexes
    LOGICAL          :: succ4                 ! Success of fourdim_fit
    DOUBLE PRECISION :: rms_min               ! Minimum RMS among the succesful sampling points
    INTEGER          :: verb_dif_old          ! Temporary verbosity level
    INTEGER          :: num_zer               ! Number of null eigenvalues of C_A
    DOUBLE PRECISION :: out_mat(4,4)          ! Output matrix of QR inversion
    DOUBLE PRECISION :: eig_val(4)            ! Eigenvalues of C_A
    DOUBLE PRECISION :: det_aux               ! Auxiliary determinant
    DOUBLE PRECISION :: delnor                ! Correction norm
    DOUBLE PRECISION :: dnp,dnm               ! Nodal distances
    !=======================================================================================================
    !****************!
    ! Inizialization !
    !****************!
    verb_dif_old = verb_dif ! Anti-verbose for differential corrections
    verb_dif     = 1
    !****************************************************!
    ! 4-dim differential correction along each direction !
    !****************************************************!
    ! Initialization to determine the minimum value of the RMS
    rms_min = 10000.d0
    dir_loop: DO i=0,ndir
       eltemp(i,0) = el_nom ! For j=0, give the nominal elements or the attributable
       pt_loop: DO j=0,np
          IF(i.EQ.0 .AND. j.NE.0) CYCLE
          IF(i.NE.0 .AND. j.EQ.0) CYCLE
          IF(i.EQ.0 .AND. j.EQ.0 .AND. (.NOT.use_nominal)) CYCLE
          IF(mov_orbit(i,j)%conn_comp.GT.0)THEN
             IF(i.EQ.0 .AND. j.EQ.0)THEN
                mov_orbit(0,0)%el_mov     = el_nom
                mov_orbit(0,0)%unc_mov    = unc_nom
                mov_orbit(0,0)%rms_mov    = rms
                mov_orbit(0,0)%chi_mov    = 0.d0
                mov_orbit(0,0)%succ_mov   = .TRUE.
                ! Set minimum RMS as the RMS of the nominal
                rms_min = rms
             ELSE
                !---------------------------------------------------------!
                ! If the diff. corr. on the previous points were          !
                ! successful, we give the obtained orbit as a first       !
                ! guess. Otherwise we give the elements el. If j.NE.1 and !
                ! we have a nominal, we use a CYCLE to change direction.  !
                !---------------------------------------------------------!
                IF(j.EQ.1)THEN
                   eltemp(i,j) = el_nom
                ELSE
                   IF(mov_orbit(i,j-1)%succ_mov)THEN
                      eltemp(i,j) = eltemp(i,j-1)
                   ELSE
                      eltemp(i,j) = el_nom
                      IF(use_nominal .AND. k_ind.EQ.2 .AND. mov_orbit(0,0)%conn_comp.GT.0)THEN
                         !---------------------------------------------------------------------!
                         ! Remark: the nominal solution has to be inside the AR, otherwise the !
                         ! ray starting at the nominal intersect the AR only for j>0           !
                         !---------------------------------------------------------------------!
                         WRITE(iun_log,'(A,I3,A,I3,A)') 'mov: point i = ',i,', j = ',j, &
                              &           ' with prev. not conv. in cobweb sampling. Change direction.'
                         WRITE(*,'(A,I3,A,I3,A)') 'mov: point i = ',i,', j = ',j, &
                              &     ' with prev. not conv. in cobweb sampling. Change direction.'
                         CYCLE dir_loop
                      END IF
                   ENDIF
                END IF
                ! (rho_i,rhodot_j) in the last two coordinates
                eltemp(i,j)%coord(5:6) = mov_orbit(i,j)%r_rdot(1:2)
                !*********************************************!
                ! Doubly constrained differential corrections !
                !*********************************************!
                WRITE(iun_log,'(A24,I3,2X,I3)') ' Fit orbit with indices ',i,j
                CALL fourdim_fit(m,obs,obsw,eltemp(i,j),mov_orbit(i,j)%el_mov,nd,unc4_mov_loc(i,j), &
                     &           mov_orbit(i,j)%unc_mov,mov_orbit(i,j)%rms_mov,delnor,rmsh(i,j),    &
                     &           nused(i,j),succ4,des_mat,scal_obs)
                ! Successful 4-dim diff. cor: MOV orbit reached
                mov_orbit(i,j)%succ_mov = mov_orbit(i,j)%unc_mov%succ.AND.succ4
             END IF
             IF(mov_orbit(i,j)%succ_mov .AND. mov_orbit(i,j)%el_mov%h_mag.LE.hmax)THEN
                !---------------------------!
                ! Compute MOID of MOV orbit !
                !---------------------------!
                CALL nomoid(mov_orbit(i,j)%el_mov%t,mov_orbit(i,j)%el_mov,mov_orbit(i,j)%moid,dnp,dnm)
                !-------------------------------------------!
                ! Compute velocity at infinity of MOV orbit !
                !-------------------------------------------!
                mov_orbit(i,j)%v_infty = v_infty0(mov_orbit(i,j)%el_mov)
                !--------------------------------------------------!
                ! Find the minimum RMS among the successful points !
                !--------------------------------------------------!
                rms_min = MIN(mov_orbit(i,j)%rms_mov,rms_min)
                !***************************************************!
                ! Computation of the determinant of the AR->MOV map !
                !***************************************************!
                !-----------------------------------------------------!
                ! Full non-linear computation (cf. [2], Theorem 2.21) !
                !-----------------------------------------------------!
!!$                DO r=1,nd
!!$                   norm_fact(r) = NORM2(unc_mov_mov_loc(i,j)%c(1:nd,r))
!!$                END DO
!!$                DO r=1,nd
!!$                   DO s=1,nd
!!$                      C_norm(r,s) = unc_mov_mov_loc(i,j)%c(r,s)/(norm_fact(r)*norm_fact(s))
!!$                   END DO
!!$                END DO
!!$                err = 1.d2*EPSILON(1.d0)
!!$                P_mat(1:nd,1:nd) = C_norm(1:nd,1:nd)
!!$                CALL tchol(P_mat(1:nd,1:nd),nd,nd,flag,err)
!!$                IF(flag.NE.0)THEN
!!$                   WRITE(*,*) 'Non-zero flag for tchol, flag,i,j', flag,i,j
!!$                END IF
!!$                DO r=1,nd
!!$                   P_mat(r,1:nd) = P_mat(r,1:nd)*norm_fact(r)
!!$                END DO
!!$                det_P(i,j)= 1.d0
!!$                DO r=1,nd
!!$                   det_P(i,j) = det_P(i,j)*P_mat(r,r)
!!$                END DO
                !---------------------------------------------------------!
                ! Gramian determinant of f_mu : AR -> MOV ([2], Eq. 2.16) !
                !---------------------------------------------------------!
                b_rho(1:scal_obs,:) = des_mat(1:scal_obs,5:6)
                b_A(1:scal_obs,:)   = des_mat(1:scal_obs,1:4)
                b_A_rho             = MATMUL(TRANSPOSE(b_A(1:scal_obs,:)),b_rho(1:scal_obs,:))
                dAdrho              = MATMUL(unc4_mov_loc(i,j)%g(1:4,1:4),b_A_rho)
                M_A                 = MATMUL(TRANSPOSE(dAdrho),dAdrho)
                det_aux             = (M_A(1,1)+1.d0)*(M_A(2,2)+1.d0)-M_A(1,2)*M_A(2,1)
                IF(det_aux.LT.0.d0)THEN
                   WRITE(*,*) 'mov: det(M_mu)^2 is negative at sample point ',i,j
                   det_aux = 0.d0
                END IF
                mov_orbit(i,j)%jac_mu = SQRT(det_aux)
                IF(mov_orbit(i,j)%jac_mu.LE.jac_mu_min) WRITE(*,*) 'mov: quasi-singular matrix M_mu (AR to MOV)'
             END IF
          ELSE
             mov_orbit(i,j)%succ_mov = .FALSE.
          END IF
       END DO pt_loop
    END DO dir_loop
    !***********************************!
    !  Computation of the chi function  !
    !***********************************!
    nused(0,0) = 2*m
    IF(use_nominal .AND. rms_min.LT.rms)THEN
       WRITE(*,*) 'mov: WARNING! Found RMS lower than the nominal one. RMS = ', rms_min
    END IF
    chi_min  = chi_min_aux
    i_chimin = 0
    j_chimin = 0
    DO i=0,ndir
       DO j=0,np
          IF(i.EQ.0 .AND. j.NE.0) CYCLE
          IF(i.NE.0 .AND. j.EQ.0) CYCLE
          IF(i.EQ.0 .AND. j.EQ.0 .AND. (.NOT.use_nominal)) CYCLE
          IF(.NOT.mov_orbit(i,j)%succ_mov) CYCLE
          IF(mov_orbit(i,j)%el_mov%h_mag.GT.hmax) CYCLE
          mov_orbit(i,j)%chi_mov = SQRT(ABS(mov_orbit(i,j)%rms_mov**2-rms_min**2)*nused(i,j))
          ! Find minimum chi and corresponding indices
          IF(mov_orbit(i,j)%chi_mov.LT.chi_min)THEN
             chi_min = mov_orbit(i,j)%chi_mov
             i_chimin = i
             j_chimin = j
          END IF
       END DO
    END DO
    ! Output of unc4_mov (if present)
    IF(PRESENT(unc4_mov)) unc4_mov = unc4_mov_loc
    ! Restore the beginning verbosity level
    verb_dif = verb_dif_old

  END SUBROUTINE mov


  !============================================================!
  ! WRITE_MOVSAMPLE_FILE                                       !
  !============================================================!
  ! Write output file of MOV sampling.                         !
  !============================================================!
  SUBROUTINE write_movsample_file(obj_name,att_t,k_ind)
    !=======================================================================================================
    CHARACTER*(name_len), INTENT(IN) :: obj_name    ! Object name
    DOUBLE PRECISION,     INTENT(IN) :: att_t       ! Time of attributable
    INTEGER,              INTENT(IN) :: k_ind       ! Grid number (first or second)
    !=======================================================================================================
    INTEGER             :: i,j               ! Loop indices
    CHARACTER(LEN=100)  :: file_mov          ! Cometarian file
    INTEGER             :: iun_mov           ! Unit for cometarian file
    INTEGER             :: le,len_mov        ! File lenght
    CHARACTER(LEN=18)   :: samp_type         ! Sample type
    !=======================================================================================================

    ! Open output file
    IF(k_ind.EQ.1)THEN
       file_mov = obj_name//'_1.mov_sample'
    ELSEIF(k_ind.EQ.2)THEN
       file_mov = obj_name//'.mov_sample'
    END IF
    CALL rmsp(file_mov,len_mov)
    CALL filopn(iun_mov,file_mov(1:len_mov),'UNKNOWN')
    !---------------------------!
    ! Header of the output file !
    !---------------------------!
    CALL write_movsample_header(iun_mov,obj_name,att_t,'COM')
    IF(good_first)THEN
       !************************!
       ! Loop on the MOV orbits !
       !************************!
       ! MOV orbit with minimum chi as first orbit of the output files
       CALL write_movsample_record(iun_mov,i_chimin,j_chimin,mov_orbit(i_chimin,j_chimin)%el_mov,              &
            &                      mov_orbit(i_chimin,j_chimin)%chi_mov,mov_orbit(i_chimin,j_chimin)%rms_mov,  &
            &                      mov_orbit(i_chimin,j_chimin)%conn_comp,mov_orbit(i_chimin,j_chimin)%moid,   &
            &                      mov_orbit(i_chimin,j_chimin)%v_infty,mov_orbit(i_chimin,j_chimin)%imp_flag, &
            &                      mov_orbit(i_chimin,j_chimin)%e_earth)
       WRITE(iun_mov,'(A)') '%'
       DO i=0,ndir
          DO j=0,np
             IF(i.EQ.0 .AND. j.NE.0) CYCLE
             IF(i.NE.0 .AND. j.EQ.0) CYCLE
             IF(i.EQ.0 .AND. j.EQ.0 .AND. (.NOT.use_nominal)) CYCLE
             ! Orbit with minimum chi already written
             IF(i.EQ.i_chimin .AND. j.EQ.j_chimin) CYCLE
             ! Select MOV orbits actually computed
             IF(mov_orbit(i,j)%succ_mov)THEN
                !----------------------------------------------!
                ! Write record in the MOV sampling output file !
                !----------------------------------------------!
                CALL write_movsample_record(iun_mov,i,j,mov_orbit(i,j)%el_mov,mov_orbit(i,j)%chi_mov, &
                     &                      mov_orbit(i,j)%rms_mov,mov_orbit(i,j)%conn_comp,          &
                     &                      mov_orbit(i,j)%moid,mov_orbit(i,j)%v_infty,               &
                     &                      mov_orbit(i,j)%imp_flag,mov_orbit(i,j)%e_earth)
             END IF
          END DO
       END DO
    END IF
    CALL filclo(iun_mov,' ')

  END SUBROUTINE write_movsample_file


  !============================================================!
  ! WRITE_ARSAMPLE_FILE                                        !
  !============================================================!
  ! Write output file of AR sampling.                          !
  !============================================================!
  SUBROUTINE write_arsample_file(obj_name,angles,coeff,roots,k_ind)
    !=======================================================================================================
    CHARACTER*(name_len), INTENT(IN) :: obj_name    ! Object name
    DOUBLE PRECISION,     INTENT(IN) :: angles(4)   ! Attributable angles and derivative
    DOUBLE PRECISION,     INTENT(IN) :: coeff(0:5)  ! Coefficients to compute V(r)
    DOUBLE PRECISION,     INTENT(IN) :: roots(3)    ! Positive roots of V(r)
    INTEGER,              INTENT(IN) :: k_ind       ! Grid number (first or second)
    !=======================================================================================================
    INTEGER             :: i,j               ! Loop indices
    CHARACTER(LEN=100)  :: file_ar           ! Cobweb file
    INTEGER             :: iun_ars           ! Unit for web file
    INTEGER             :: le,len_ar         ! File lenght
    CHARACTER(LEN=18)   :: samp_type         ! Sample type
    !=======================================================================================================
    ! Define sampling type string for output
    IF(.NOT.use_nominal)THEN
       IF(log_grid)THEN
          samp_type = 'Grid - Logarithmic'
       ELSE
          samp_type = 'Grid - Uniform'
       END IF
    ELSE
       samp_type = 'Spider web'
    END IF
    ! Open output file
    IF(k_ind.EQ.1)THEN
       file_ar  = obj_name//'_1.ar_sample'
    ELSEIF(k_ind.EQ.2)THEN
       file_ar  = obj_name//'.ar_sample'
    END IF
    CALL rmsp(file_ar,len_ar)
    CALL filopn(iun_ars,file_ar(1:len_ar),'UNKNOWN')
    !---------------------------!
    ! Header of the output file !
    !---------------------------!
    CALL write_arsample_header(iun_ars,obj_name,angles,coeff,roots,rdotmin_ar,rdotmax_ar,samp_type)
    !************************!
    ! Loop on the MOV orbits !
    !************************!
    ! MOV orbit with minimum chi as first orbit of the output files
    CALL write_arsample_record(iun_ars,i_chimin,j_chimin,mov_orbit(i_chimin,j_chimin)%r_rdot(1),               &
         &                     mov_orbit(i_chimin,j_chimin)%r_rdot(2), mov_orbit(i_chimin,j_chimin)%conn_comp, &
         &                     mov_orbit(i_chimin,j_chimin)%succ_mov,mov_orbit(i_chimin,j_chimin)%chi_mov,     &
         &                     mov_orbit(i_chimin,j_chimin)%jac_mu,mov_orbit(i_chimin,j_chimin)%e_sun,         &
         &                     mov_orbit(i_chimin,j_chimin)%e_earth,mov_orbit(i_chimin,j_chimin)%el_mov%h_mag)
    WRITE(iun_ars,'(A)') '%'
    DO i=0,ndir
       DO j=0,np
          IF(i.EQ.0 .AND. j.NE.0) CYCLE
          IF(i.NE.0 .AND. j.EQ.0) CYCLE
          IF(i.EQ.0 .AND. j.EQ.0 .AND. (.NOT.use_nominal)) CYCLE
          ! Orbit with minimum chi already written
          IF(i.EQ.i_chimin .AND. j.EQ.j_chimin) CYCLE
          !---------------------------------------------!
          ! Write record in the AR sampling output file !
          !---------------------------------------------!
          CALL write_arsample_record(iun_ars,i,j,mov_orbit(i,j)%r_rdot(1),mov_orbit(i,j)%r_rdot(2),           &
               &                     mov_orbit(i,j)%conn_comp,mov_orbit(i,j)%succ_mov,mov_orbit(i,j)%chi_mov, &
               &                     mov_orbit(i,j)%jac_mu,mov_orbit(i,j)%e_sun,mov_orbit(i,j)%e_earth,       &
               &                     mov_orbit(i,j)%el_mov%h_mag)
       END DO
    END DO
    CALL filclo(iun_ars,' ')

  END SUBROUTINE write_arsample_file



  !====================================================================!
  ! COMPUTE_SCORE                                                      !
  !====================================================================!
  ! Compute the score of the object to be a NEO, a MBA, a distant      !
  ! object, or everything else.                                        !
  !====================================================================!
  SUBROUTINE compute_score(score_dist,score_mba,score_neo,score_scatt,score_pha)
    DOUBLE PRECISION, INTENT(OUT) :: score_neo     ! Score of the object: NEO
    DOUBLE PRECISION, INTENT(OUT) :: score_mba     ! Score of the object: MBA / Trojan
    DOUBLE PRECISION, INTENT(OUT) :: score_dist    ! Score of the object: Distant
    DOUBLE PRECISION, INTENT(OUT) :: score_scatt   ! Score of the object: Scattered
    DOUBLE PRECISION, INTENT(OUT) :: score_pha     ! Score of the object: PHA
    !=======================================================================================================
    DOUBLE PRECISION :: nposs                         ! Number of possible points
    DOUBLE PRECISION :: n_neo,n_mba,n_dis,n_sca,n_pha ! Number of possible points per class
    LOGICAL          :: mba                           ! Flag for MBA category
    INTEGER          :: i,j                           ! Loop indexes
    TYPE(orbit_elem) :: elcom                         ! Cometary orbital elements
    INTEGER          :: fail_flag                     ! Fail flag in changing coordinates
    DOUBLE PRECISION :: q_NEO, q_dist                 ! Perihelion for NEO and for distant objects
    DOUBLE PRECISION :: q,ecc,sma                     ! Perihelion, eccentricity, semimajor axis
    DOUBLE PRECISION :: sum                           ! Sum of the four scores
    !=======================================================================================================
    nposs = 0
    n_neo = 0
    n_mba = 0
    n_dis = 0
    n_sca = 0
    n_pha = 0
    !------------------------------------------------------------!
    ! CATEGORIES DEFINITION. Each class of objects has a proper  !
    ! definition.                                                !
    !    - NEO: q < 1.3 au                                       !
    !                                                            !
    !    - MBA: non-NEO and ((1.7 < a < 4.5 and e < 0.4) or      !
    !                       (4.5 < a < 5.5 and e < 0.3))         !
    !                                                            !
    !    - DISTANT: q > 28 au                                    !
    !                                                            !
    !    - SCATTERED: everything else                            !
    !                                                            !
    !    - PHA: MOID < 0.05 au and H < 22                        !
    !------------------------------------------------------------!
    q_NEO  = 1.3d0
    q_dist = 28.d0
    DO i=0,ndir
       DO j=0,np
          IF(i.EQ.0 .AND. j.NE.0) CYCLE
          IF(i.NE.0 .AND. j.EQ.0) CYCLE
          IF(i.EQ.0 .AND. j.EQ.0 .AND. (.NOT.use_nominal)) CYCLE
          IF(.NOT.mov_orbit(i,j)%succ_mov .OR. mov_orbit(i,j)%chi_mov.GT.5.d0) CYCLE
          IF(mov_orbit(i,j)%el_mov%h_mag.GT.hmax) CYCLE
          nposs = nposs + EXP(-mov_orbit(i,j)%chi_mov**2/2)*mov_orbit(i,j)%jac_sigma*mov_orbit(i,j)%jac_mu
          CALL coo_cha(mov_orbit(i,j)%el_mov,'COM',elcom,fail_flag)
          ! Perihelion and aphelion
          q   = elcom%coord(1)
          ecc = elcom%coord(2)
          sma = q/(1.d0-ecc)
          mba = (sma.LT.4.5d0 .AND. sma.GT.1.7d0 .AND. ecc.LT.0.4d0).OR.(sma.LT.5.5d0 .AND. sma.GT.4.5d0 .AND. ecc.LT.0.3d0)
          IF(q.LT.q_NEO)THEN
             n_neo = n_neo + EXP(-mov_orbit(i,j)%chi_mov**2/2)*mov_orbit(i,j)%jac_sigma*mov_orbit(i,j)%jac_mu
          ELSEIF(mba)THEN
             n_mba = n_mba + EXP(-mov_orbit(i,j)%chi_mov**2/2)*mov_orbit(i,j)%jac_sigma*mov_orbit(i,j)%jac_mu
          ELSEIF(q.GT.q_dist)THEN
             n_dis = n_dis + EXP(-mov_orbit(i,j)%chi_mov**2/2)*mov_orbit(i,j)%jac_sigma*mov_orbit(i,j)%jac_mu
          ELSE
             n_sca = n_sca + EXP(-mov_orbit(i,j)%chi_mov**2/2)*mov_orbit(i,j)%jac_sigma*mov_orbit(i,j)%jac_mu
          END IF
          IF(mov_orbit(i,j)%moid.LT.0.05d0 .AND. mov_orbit(i,j)%el_mov%h_mag.LT.22.d0)THEN
             n_pha = n_pha + EXP(-mov_orbit(i,j)%chi_mov**2/2)*mov_orbit(i,j)%jac_sigma*mov_orbit(i,j)%jac_mu
          END IF
       END DO
    END DO
    IF(nposs.NE.0)THEN
       score_neo   = n_neo/nposs
       score_mba   = n_mba/nposs
       score_dist  = n_dis/nposs
       score_scatt = n_sca/nposs
       score_pha   = n_pha/nposs
       !-------------------!
       ! Output adjustment !
       !-------------------!
       IF(score_neo.GT.1.d0)   score_neo   = 1.d0
       IF(score_mba.GT.1.d0)   score_mba   = 1.d0
       IF(score_dist.GT.1.d0)  score_dist  = 1.d0
       IF(score_scatt.GT.1.d0) score_scatt = 1.d0
       IF(score_pha.GT.1.d0)   score_pha = 1.d0
       sum = score_neo + score_mba + score_dist + score_scatt
       IF(sum.NE.1.d0)THEN
          score_neo   = score_neo/sum
          score_mba   = score_mba/sum
          score_dist  = score_dist/sum
          score_scatt = score_scatt/sum
       END IF
    ELSE
       WRITE(iun_log,*) 'WARNING: skipped score computation because nposs = 0'
       WRITE(*,*) 'WARNING: skipped score computation because nposs = 0'
    END IF

  END SUBROUTINE compute_score


  !====================================================================!
  ! IMMEDIATE_IMP                                                      !
  !====================================================================!
  ! Search for immediate impactors, given cobweb or grid.              !
  !====================================================================!
  SUBROUTINE immediate_imp(obj_name,t_0,dt_im,res_norm,arc_type,geoc_chi,acce_chi,chi_curv,m,obs)
    USE propag_state, ONLY: pro_ele
    USE tp_trace,     ONLY: wri_tppoint, dtpde, tpplane, mtp_store, njc
    INCLUDE 'vers.h90'
    !=======================================================================================================
    CHARACTER*(name_len), INTENT(IN) :: obj_name  ! Object name
    DOUBLE PRECISION,     INTENT(IN) :: t_0       ! Orbital elements epoch
    DOUBLE PRECISION,     INTENT(IN) :: dt_im     ! Time (days) for the scan
    DOUBLE PRECISION,     INTENT(IN) :: res_norm  ! Residuals norm
    INTEGER,              INTENT(IN) :: arc_type  ! Arc type
    DOUBLE PRECISION,     INTENT(IN) :: geoc_chi  ! Normalized value of the geodetic curvature
    DOUBLE PRECISION,     INTENT(IN) :: acce_chi  ! Normalized value for the along-track acceleration
    DOUBLE PRECISION,     INTENT(IN) :: chi_curv  ! Chi-value of geodetic curvature and acceleration
    INTEGER,              INTENT(IN) :: m         ! Number of observations
    TYPE(ast_obs),        INTENT(IN) :: obs(m)    ! Observations
    !===== Impact Probability computation ==================================================================
    DOUBLE PRECISION :: imppro(0:nmax,0:nmax) ! Impact probability
    INTEGER          :: nposs                 ! Number of points with chi <= 5, H <= H_max
    INTEGER          :: nimp                  ! Number of impacting points
    TYPE(orbit_elem) :: el1                   ! Propagation of a sampling orbit
    DOUBLE PRECISION :: pronorm               ! Normalization factor for the Impact Probability
    DOUBLE PRECISION :: iptot                 ! Total IP
    !===== Auxiliary data ==================================================================================
    DOUBLE PRECISION   :: p_earth_bound ! Probability to be an Earth-bounded object
    DOUBLE PRECISION   :: h_mag_min     ! Minimum H among MOV impactors
    DOUBLE PRECISION   :: h_mag_max     ! Maximum H among MOV impactors
    DOUBLE PRECISION   :: t_min_imp     ! Minimum impact time among MOV impactors
    DOUBLE PRECISION   :: t_max_imp     ! Maximum impact among MOV impactors
    DOUBLE PRECISION   :: v_inf_min     ! Minimum velocity at infinity
    DOUBLE PRECISION   :: v_inf_max     ! Maximum velocity at infinity
    DOUBLE PRECISION   :: hour_min      ! Hour of t_min_imp
    DOUBLE PRECISION   :: tf_search     ! Final time of VI search
    INTEGER            :: iday_min      ! Day of t_min_imp
    INTEGER            :: imonth_min    ! Month of t_min_imp
    INTEGER            :: iyear_min     ! Year of t_min_imp
    DOUBLE PRECISION   :: hour_max      ! Hour of t_max_imp
    INTEGER            :: iday_max      ! Day of t_max_imp
    INTEGER            :: imonth_max    ! Month of t_max_imp
    INTEGER            :: iyear_max     ! Year of t_max_imp
    CHARACTER(LEN=16)  :: t_min_str     ! Minimum closest approach distance among MOV impactors (CAL)
    CHARACTER(LEN=16)  :: t_max_str     ! Maximum closest approach distance among MOV impactors (CAL)
    !===== Other variables =================================================================================
    CHARACTER(LEN=8)   :: date               ! Date of computation
    CHARACTER(LEN=10)  :: time               ! Time of computation
    INTEGER            :: i,j,k              ! Loop indexes
    CHARACTER(LEN=100) :: file               ! File name
    INTEGER            :: le                 ! Length of file
    INTEGER            :: iun_mtp            ! File unit for the .mtp file (MTP analysis)
    INTEGER            :: iunrisk            ! File unit for the risk file
    CHARACTER(LEN=20)  :: mulname            ! Name of the Virtual Asteroid
    INTEGER            :: dataflag           ! Data policy flag (values from 0 to 4)
    LOGICAL            :: mtp(0:nmax,0:nmax) ! MTP crossing flag
    DOUBLE PRECISION   :: d_min
    INTEGER            :: jc
    DOUBLE PRECISION   :: dt                 ! Arc length
    !=======================================================================================================
    
    ! Select the MTP plane (not the TP)
    tpplane = .FALSE.
    ! Final time for VI search
    tf_search = t_0 + dt_im
    ! Clean up the previous close app. records
    REWIND(iuncla)

    !*****************************************!
    !  COMPUTATION OF THE IMPACT PROBABILITY  !
    !*****************************************!
    nposs         = 0      ! Counter of "good" points
    pronorm       = 0.d0   ! Normalization factor
    nimp          = 0
    imppro        = 0.d0
    iptot         = 0.d0
    p_earth_bound = 0.d0
    h_mag_min     = 100.d0
    h_mag_max     = 0.d0
    t_min_imp     = 100000.d0
    t_max_imp     = 0.d0
    v_inf_min     = 1.d8
    v_inf_max     = 0.d0
    mtp           = .FALSE.
    DO i=0,ndir
       DO j=0,np
          IF(i.EQ.0 .AND. j.NE.0) CYCLE
          IF(i.NE.0 .AND. j.EQ.0) CYCLE
          IF(i.EQ.0 .AND. j.EQ.0 .AND. (.NOT.use_nominal)) CYCLE
          IF(.NOT.mov_orbit(i,j)%succ_mov .OR. mov_orbit(i,j)%chi_mov.GT.5.d0) CYCLE
          IF(mov_orbit(i,j)%el_mov%h_mag.GT.hmax) CYCLE
          !---------------------------------------------------!
          ! Normalization factor for IP (cf. [2], Eq. (2.19)) !
          !---------------------------------------------------!
          nposs   = nposs + 1
          pronorm = pronorm + EXP(-mov_orbit(i,j)%chi_mov**2/2.d0)*mov_orbit(i,j)%jac_mu*mov_orbit(i,j)%jac_sigma
          IF(mov_orbit(i,j)%e_earth.LE.0.d0)THEN
             !----------------------!
             ! Earth-bounded orbits !
             !----------------------!
             p_earth_bound = p_earth_bound + EXP(-mov_orbit(i,j)%chi_mov**2/2.d0)*mov_orbit(i,j)%jac_mu*mov_orbit(i,j)%jac_sigma
             IF(.NOT.propag_geoc_orbit) CYCLE
          END IF
          WRITE(mulname,590) i,j
590       FORMAT('cob_',I2,'_',I2)
          CALL rmsp(mulname,le)
          WRITE(iuncla,*) mulname
          CALL cov_not_av ! Propagate without covariance
          njc = 0         ! Number of minima during a close approach
          WRITE(*,*) 'MOV orbit ',i,'  ',j,' propagation'
          WRITE(iun_log,*) 'MOV orbit ',i,'  ',j,' propagation'
          CALL pro_ele(mov_orbit(i,j)%el_mov,tf_search,el1)
          !****************!
          !  MTP analysis  !
          !****************!
          IF(njc.GT.0)THEN
             !----------------------------------------------------!
             ! In case of multiple minima, select the closest one !
             !----------------------------------------------------!
             IF(njc.EQ.1)THEN 
                jc = 1 
             ELSE 
                d_min = 1.d10 
                DO k=1,njc
                   IF(mtp_store(k)%rcla.LT.d_min)THEN 
                      jc = k 
                      d_min = mtp_store(k)%rcla 
                   END IF
                END DO
                WRITE(*,*) ' Multiple minima: njc = ',njc,', selected = ',jc
                WRITE(iun_log,*) ' Multiple minima: njc = ',njc,', selected = ',jc
             END IF
             IF(mtp_store(jc)%iplam.EQ.3)THEN
                mtp(i,j) = .TRUE.
                ! MTP data
                mov_orbit(i,j)%csi  = mtp_store(jc)%tp_coord(1)/reau
                mov_orbit(i,j)%zeta = mtp_store(jc)%tp_coord(2)/reau
                mov_orbit(i,j)%tcla = mtp_store(jc)%tcla
                mov_orbit(i,j)%dcla = mtp_store(jc)%rcla
                ! Contribution of impacting points to the IP computation (numerator)
                IF(mov_orbit(i,j)%csi**2+mov_orbit(i,j)%zeta**2.LT.1.d0)THEN
                   mov_orbit(i,j)%imp_flag = 1
                   nimp = nimp + 1
                   !------------------------------------------------!
                   ! Contribution from the VI (cf. [2], eq. (2.22)) !
                   !------------------------------------------------!
                   imppro(i,j) = EXP(-mov_orbit(i,j)%chi_mov**2/2.d0)*mov_orbit(i,j)%jac_mu*mov_orbit(i,j)%jac_sigma
                   iptot       = iptot + imppro(i,j)
                   !************************************************!
                   ! Auxiliary data: H and impact time of impactors !
                   !************************************************!
                   h_mag_max = MAX(h_mag_max,mov_orbit(i,j)%el_mov%h_mag)
                   h_mag_min = MIN(h_mag_min,mov_orbit(i,j)%el_mov%h_mag)
                   t_max_imp = MAX(t_min_imp,mov_orbit(i,j)%tcla)
                   t_min_imp = MIN(t_min_imp,mov_orbit(i,j)%tcla)
                   v_inf_max = MAX(v_inf_max,mov_orbit(i,j)%v_infty)
                   v_inf_min = MIN(v_inf_min,mov_orbit(i,j)%v_infty)
                END IF
             END IF
          END IF
       END DO
    END DO
    ! Normalisation of probabilities
    IF((use_nominal .AND. nposs.GE.2) .OR. (.NOT.use_nominal .AND. nposs.GE.1))THEN
       p_earth_bound = p_earth_bound/pronorm
       iptot         = iptot/pronorm
    END IF

    !**********************************!
    !  WRITE MOV SAMPLING OUTPUT FILE  !
    !**********************************!
    CALL write_movsample_file(obj_name,t_0,2)

    !*************************!
    !  WRITE MTP OUTPUT FILE  !
    !*************************!
    file=obj_name//'.mtp'
    CALL rmsp(file,le)
    CALL filopn(iun_mtp,file(1:le),'unknown')
    !--------------------------!
    ! Write header of MTP file !
    !--------------------------!
    CALL write_mtp_header(iun_mtp,obj_name)
    IF(mtp(i_chimin,j_chimin))THEN
       ! MTP data of MOV orbit with minimum chi
       CALL write_mtp_record(iun_mtp,i_chimin,j_chimin,mov_orbit(i_chimin,j_chimin)%tcla,        &
            &                mov_orbit(i_chimin,j_chimin)%csi,mov_orbit(i_chimin,j_chimin)%zeta, &
            &                mov_orbit(i_chimin,j_chimin)%dcla,mov_orbit(i_chimin,j_chimin)%el_mov%h_mag)
       WRITE(iun_mtp,'(A)') '%'
    END IF
    DO i=0,ndir
       DO j=0,np
          IF(i.EQ.i_chimin .AND. j.EQ.j_chimin) CYCLE
          IF(mtp(i,j))THEN
             !------------------------------!
             ! Write record in the MTP file !
             !------------------------------!
             CALL write_mtp_record(iun_mtp,i,j,mov_orbit(i,j)%tcla,mov_orbit(i,j)%csi,mov_orbit(i,j)%zeta, &
                  &                mov_orbit(i,j)%dcla,mov_orbit(i,j)%el_mov%h_mag)
          END IF
       END DO
    ENDDO
    CALL filclo(iun_mtp,' ')

    !***********************!
    !  WRITE THE RISK FILE  !
    !***********************!
    dt = (obs(m)%time_tdt-obs(1)%time_tdt)*24.d0  ! Arc duration (hours)
    !-----------------------------------------------------------------------!
    ! Data policy flag assignment. We give the following data policy flag:  !
    !                                                                       !
    ! - 0 for IP <= 1E-4                                                    !
    !                                                                       !
    ! - 1 for 1E-4 < IP <= 1E-3                                             !
    !                                                                       !
    ! - 2 for 1E-3 < IP <= 1E-2                                             !
    !                                                                       !
    ! - 3 for IP > 1E-2, and when a nominal does not exists or exists       !
    !   but curvature chi^2 <= 10                                           !
    !                                                                       !
    ! - 4 for IP > 1E-2, when a nominal exists and has curvature chi^2 > 10 !
    !-----------------------------------------------------------------------!
    IF(iptot.LE.1.d-4)THEN
       dataflag = 0
    ELSEIF(iptot.LE.1.d-3)THEN
       dataflag = 1
    ELSEIF(iptot.LE.1.d-2)THEN
       dataflag = 2
    ELSEIF(use_nominal .AND. chi_curv**2.GT.10)THEN
       dataflag = 4
    ELSE
       dataflag = 3
    ENDIF
    !------------------------------------!
    ! Open risk/norisk/nosig output file !
    !------------------------------------!
    IF(.NOT.significant)THEN
       file = obj_name//'.nosig'
    ELSEIF(dataflag.GE.1)THEN
       file = obj_name//'.risk'
    ELSE
       file = obj_name//'.norisk'
    ENDIF
    CALL rmsp(file,le)
    CALL filopn(iunrisk,file(1:le),'unknown')
    WRITE(iunrisk,201) 'NEOCP     # obs      dt    Arc        Significance of curvature       ', &
         &             'RMS      Total    MOV    Imp      IP      Data '
    WRITE(iunrisk,201) 'Name                 (h)   type    Geodesic     Accel.     Overall    ', &
         &             '         # pts    pts    pts              flag '
    WRITE(iunrisk,201) '----------------------------------------------------------------------', &
         &             '-----------------------------------------------'
201 FORMAT(A70,A47)
    
    ! Content (IP in exp. format if < 1E-3, otherwise in floating point format)
    IF(iptot.GT.0.d0 .AND. iptot.LT.1.d-3)THEN
       WRITE(iunrisk,555) obj_name,m,dt,arc_type,geoc_chi,acce_chi,chi_curv, &
            &             mov_orbit(i_chimin,j_chimin)%rms_mov,ndir*np,nposs,nimp,iptot,dataflag
    ELSE
       WRITE(iunrisk,556) obj_name,m,dt,arc_type,geoc_chi,acce_chi,chi_curv, &
            &             mov_orbit(i_chimin,j_chimin)%rms_mov,ndir*np,nposs,nimp,iptot,dataflag
    END IF
    
555 FORMAT(A9,2X,I3,3X,F8.2,2X,I2,4X,F9.2,2X,F9.2,3X,F9.2,3X,F7.3,2X,I7,2X,I5,2X,I5,2X,1P,E10.3,4X,I1)
556 FORMAT(A9,2X,I3,3X,F8.2,2X,I2,4X,F9.2,2X,F9.2,3X,F9.2,3X,F7.3,2X,I7,2X,I5,2X,I5,5X,F5.3,5X,I1)

    !***************!
    !  Write notes  !
    !***************!
    !-------------------------!
    ! Elongation from the Sun !
    !-------------------------!
    WRITE(iunrisk,595) elong*degrad
    WRITE(iunrisk,*)
    !--------------------------------------------------------!
    ! Absolute magnitude, closest approach time and distance !
    !--------------------------------------------------------!
    IF(nimp.GT.0)THEN
       CALL convert_mjd_cal(t_min_imp,t_min_str)
       CALL convert_mjd_cal(t_max_imp,t_max_str)
       WRITE(iunrisk,596) h_mag_min,h_mag_max
       WRITE(iunrisk,597) t_min_str,t_max_str
       WRITE(iunrisk,598) v_inf_min,v_inf_max
    ENDIF
    !-------------------------------------------!
    ! Probability to be an Earth-bounded object !
    !-------------------------------------------!
    WRITE(iunrisk,599) p_earth_bound
    !---------------------!
    ! Shooting star limit !
    !---------------------!
    WRITE(iunrisk,600) hmax
    WRITE(iunrisk,*)
    !-------------------------------------!
    ! OrbFit version, date of computation !
    !-------------------------------------!
    CALL date_and_time(date,time)
    WRITE(iunrisk,121) pvers, date, time
    ! Reason for which the case is not significant
    IF(.NOT.significant)THEN
       WRITE(iunrisk,*)
       IF(m.LE.2)THEN
          WRITE(iunrisk,122) m
       ELSEIF(dt.LE.0.5d0)THEN
          WRITE(iunrisk,123) dt
       ENDIF
    ENDIF
595 FORMAT(/'Elongation from the Sun = ',F5.1)
596 FORMAT(/'Absolute magnitude of impactors: H_min = ',F5.2,' and H_max = ',F5.2)
597 FORMAT(/'Closest approach time of impactors: t_min = ',A16,' TDB and t_max = ',A16,' TDB')
598 FORMAT(/'Velocity at infinity of impactors: v_inf_min = ',F7.3,' [km/s] and v_inf_max = ',F7.3,' [km/s]')
599 FORMAT(/'Probability to be an Earth-bounded object = ', F5.3)
600 FORMAT(/'Shooting star limit: ',F5.2, ' magnitudes')
121 FORMAT(/'OrbFit software version ',A17,' Date of computation = ',A10,1X,A10,' CET')
122 FORMAT(/'This case is not significant because the number of observations is ',I3)
123 FORMAT(/'This case is not significant because dt is ',F8.2,' hours')
    CALL filclo(iunrisk,' ')

  END SUBROUTINE immediate_imp


 !============================================================!
 ! ENERGY_SUN                                                 !
 !============================================================!
 ! Energy w.r.t. the Sun and connected component selection.   !
 !============================================================!
 SUBROUTINE energy_sun(c,roots,rmin)
   !=======================================================================================================
   DOUBLE PRECISION, INTENT(IN) :: c(0:5)        ! Coefficients to construct V(r)
   DOUBLE PRECISION, INTENT(IN) :: roots(3)      ! Positive roots of V(r)
   DOUBLE PRECISION, INTENT(IN) :: rmin          ! Minimum value for rho
   !=======================================================================================================
   INTEGER          :: i,j         ! Loop indices
   DOUBLE PRECISION :: W_rho,S_rho ! Computation of the energy
   !=======================================================================================================
   DO i=0,ndir
      DO j=0,np
         ! Initialization
         mov_orbit(i,j)%e_sun = 0.d0
         IF(i.EQ.0 .AND. j.NE.0) CYCLE
         IF(i.NE.0 .AND. j.EQ.0) CYCLE
         IF(i.EQ.0 .AND. j.EQ.0 .AND. (.NOT.use_nominal)) CYCLE
         ! Auxiliary quantities to compute the two-body energy w.r.t. the Sun
         W_rho = c(2)*mov_orbit(i,j)%r_rdot(1)**2 + c(3)*mov_orbit(i,j)%r_rdot(1) + c(4)
         S_rho = mov_orbit(i,j)%r_rdot(1)**2 + c(5)*mov_orbit(i,j)%r_rdot(1) + c(0)
         ! Two-body energy w.r.t. the Sun
         mov_orbit(i,j)%e_sun = (mov_orbit(i,j)%r_rdot(2)**2 + c(1)*mov_orbit(i,j)%r_rdot(2) + &
              &                  W_rho - 2*gms/SQRT(S_rho))*0.5d0
         ! Exclude from the AR the points with E > E_lim
         IF(mov_orbit(i,j)%r_rdot(1).GT.rmin .AND. mov_orbit(i,j)%e_sun.LE.E_lim)THEN
            IF(roots(2).GT.0 .AND. mov_orbit(i,j)%r_rdot(1).GE.roots(2))THEN
               mov_orbit(i,j)%conn_comp = 2
            ELSE
               mov_orbit(i,j)%conn_comp = 1
            END IF
         END IF
      END DO
   END DO

 END SUBROUTINE energy_sun


 !============================================================!
 ! ENERGY_EARTH                                               !
 !============================================================!
 ! Energy w.r.t. the Earth.                                   !
 !============================================================!
 SUBROUTINE energy_earth(att,xogeo,vogeo,allorbits)
   !=======================================================================================================
   DOUBLE PRECISION, INTENT(IN)  :: att(4)      ! Attributables
   DOUBLE PRECISION, INTENT(IN)  :: xogeo(3)    ! Geocentric observer coordinates, equatorial
   DOUBLE PRECISION, INTENT(IN)  :: vogeo(3)    ! Geocentric observer coordinates velocity, equatorial
   LOGICAL,          INTENT(IN)  :: allorbits   ! Compute energy for all the orbits or just for the nominal
   !=======================================================================================================
   INTEGER          :: i,j        ! Loop indices
   DOUBLE PRECISION :: rhat(3)    ! Unit vector in the obs direction
   DOUBLE PRECISION :: rhalpha(3) ! Derivative of rhat w.r.t. alpha
   DOUBLE PRECISION :: rhdelta(3) ! Derivative of rhat w.r.t. delta
   !=======================================================================================================

   IF(allorbits .OR. use_nominal)THEN
      ! Unit vector along observer (geocentric equatorial)
      rhat(1)    =  COS(att(1))*COS(att(2))
      rhat(2)    =  SIN(att(1))*COS(att(2))
      rhat(3)    =  SIN(att(2))
      rhalpha(1) = -SIN(att(1))*COS(att(2))
      rhalpha(2) =  COS(att(1))*COS(att(2))
      rhalpha(3) =  0.d0
      rhdelta(1) = -COS(att(1))*SIN(att(2))
      rhdelta(2) = -SIN(att(1))*SIN(att(2))
      rhdelta(3) =  COS(att(2))
      IF(use_nominal .AND. .NOT.allorbits)THEN
         CALL compute_energy_earth(0,0,att,rhat,rhalpha,rhdelta)
      ELSE
         ! Energy w.r.t. Earth
         DO i=0,ndir
            DO j=0,np
               mov_orbit(i,j)%e_earth = 0.d0
               IF(i.EQ.0 .AND. j.NE.0) CYCLE
               IF(i.NE.0 .AND. j.EQ.0) CYCLE
               IF(i.EQ.0 .AND. j.EQ.0 .AND. (.NOT.use_nominal)) CYCLE
               CALL compute_energy_earth(i,j,att,rhat,rhalpha,rhdelta)
            END DO
         END DO
      END IF
   END IF

 END SUBROUTINE energy_earth


 !============================================================!
 ! COMPUTE_ENERGY_EARTH                                       !
 !============================================================!
 ! Energy w.r.t. the Earth of a single MOV orbit.             !
 !============================================================!
 SUBROUTINE compute_energy_earth(i,j,att,rhat,rhalpha,rhdelta)
   USE planet_masses, ONLY: gmearth
   !=======================================================================================================
   INTEGER,          INTENT(IN) :: i,j        ! Loop indices
   DOUBLE PRECISION, INTENT(IN) :: att(4)     ! Attributables
   DOUBLE PRECISION, INTENT(IN) :: rhat(3)    ! Unit vector in the obs direction
   DOUBLE PRECISION, INTENT(IN) :: rhalpha(3) ! Derivative of rhat w.r.t. alpha
   DOUBLE PRECISION, INTENT(IN) :: rhdelta(3) ! Derivative of rhat w.r.t. delta
   !=======================================================================================================
   DOUBLE PRECISION :: rgeo(3)    ! Geocentric position of the small body
   DOUBLE PRECISION :: rdotgeo(3) ! Geocentric velocity of the small body
   DOUBLE PRECISION :: dgeo       ! Geocentric distance of the small body
   DOUBLE PRECISION :: vgeo2      ! Squared geocentric velocity magnitude
   DOUBLE PRECISION :: prscal     ! Scalar product
   DOUBLE PRECISION :: vsize      ! Norm of a vector
   !=======================================================================================================
   rgeo    = xogeo + rhat*mov_orbit(i,j)%r_rdot(1)
   dgeo    = vsize(rgeo)
   rdotgeo = vogeo + mov_orbit(i,j)%r_rdot(2)*rhat + mov_orbit(i,j)%r_rdot(1)* &
        &    (rhalpha*att(3) + rhdelta*att(4))
   vgeo2   = prscal(rdotgeo,rdotgeo)
   ! Two-body energy w.r.t. the Earth
   mov_orbit(i,j)%e_earth = vgeo2/2 -gmearth/dgeo;
   
 END SUBROUTINE compute_energy_earth



 !===================================================================!
 ! SAME_SIDE_R1                                                      !
 !===================================================================!
 ! Logical function to assess if rmin and rmax are on the same side  !
 ! with respect to r_1.                                              !
 !===================================================================!
 LOGICAL FUNCTION same_side_r1(rmin,rmax,r_1)
   !=======================================================================================================
   DOUBLE PRECISION, INTENT(IN)  :: rmin      ! Starting point for the rho sampling
   DOUBLE PRECISION, INTENT(IN)  :: rmax      ! Ending point for the rho sampling
   DOUBLE PRECISION, INTENT(IN)  :: r_1       ! First root of the polynomial defining the AR
   !=======================================================================================================

   same_side_r1 = (rmin.LE.r_1 .AND. rmax.LE.r_1*1.05d0) .OR. (rmin.GT.r_1 .AND. rmax.GT.r_1)

 END FUNCTION same_side_r1

 
 !====================================================================!
 ! MC_RES_ATT                                                         !
 !====================================================================!
 ! Search for immediate impactors with MC method on residuals         !
 ! assuming a nominal solution is available                           !
 !====================================================================!
 SUBROUTINE mc_res_att(astnac,nmc,tf_search,m,obs,obsw,el0,obscodn,nd,nimp,isucc)
   USE tp_trace,     ONLY: tpplane, mtp_store, njc
   USE close_app,    ONLY: fix_mole, kill_propag
   USE propag_state, ONLY: pro_ele
   INCLUDE 'parobx.h90'
   !=======================================================================================================
   CHARACTER*(name_len), INTENT(IN)  :: astnac     ! Asteroid name
   INTEGER,              INTENT(IN)  :: nmc        ! Number of MC sample points
   DOUBLE PRECISION,     INTENT(IN)  :: tf_search  ! Final time of propagation
   INTEGER,              INTENT(IN)  :: obscodn    ! Observatory code
   INTEGER,              INTENT(IN)  :: m          ! Number of observations
   TYPE(ast_obs),        INTENT(IN)  :: obs(nobx)  ! Observations
   TYPE(ast_wbsr),       INTENT(IN)  :: obsw(nobx) ! Residuals and weigths
   TYPE(orbit_elem),     INTENT(IN)  :: el0        ! Nominal elements
   INTEGER,              INTENT(IN)  :: nd         ! Dimension of the par space
   INTEGER,              INTENT(OUT) :: nimp       ! Number of impacts
   INTEGER,              INTENT(OUT) :: isucc      ! Number of convergent differential corrections
   !=======================================================================================================
   CHARACTER(LEN=30) :: file        ! File name
   !===== MC computation ==================================================================================
   INTEGER           :: k           ! Index for the MC points
   INTEGER           :: kk          ! For one istance of noise
   DOUBLE PRECISION  :: rvnorm      ! Random noise
   TYPE(ast_obs)     :: obs1(nobx)  ! Perturbed observations
   TYPE(ast_wbsr)    :: obsw1(nobx) ! Perturbed weights
   INTEGER           :: icor(6)     ! Parameters to be corrected
   TYPE(orbit_elem)  :: elk         ! Keplerian elements (after diff_cor)
   TYPE(orb_uncert)  :: unck        ! Uncertainty matrices for elk
   TYPE(orbit_elem)  :: elatt       ! Attributable orbital elements
   TYPE(orbit_elem)  :: el1         ! Propagated elements
   LOGICAL           :: succ        ! Successful differential corrections
   DOUBLE PRECISION  :: csinok      ! Norm of the residuals
   DOUBLE PRECISION  :: delnok      ! Norm of the last correction
   DOUBLE PRECISION  :: csi         ! Xi on the TP
   DOUBLE PRECISION  :: zeta        ! Zeta on the TP
   !===== Other variables =================================================================================
   INTEGER           :: fail_flag   ! Failure flag for the coordinate change
   INTEGER           :: le          ! File length
   INTEGER           :: iuncob      ! Output file for the MC points
   !=======================================================================================================
   IF(nd.NE.6)THEN
      WRITE(*,*) 'mc_res_att: not ready for non-grav, nd = ',nd
      STOP
   END IF
   !******************!
   !  Inizialization  !
   !******************!
   verb_dif=1
   tpplane=.FALSE. ! Use MTP plane, not TP
   !***************************!
   !  File .mc for the output  !
   !***************************!
   file=astnac//'.mc'
   CALL rmsp(file,le)
   CALL filopn(iuncob,file(1:le),'unknown')

   !***************************************!
   !  LOOP ON INSTANCES OF GAUSSIAN NOISE  !
   !***************************************!
   nimp=0  ! Impact counter
   isucc=0 ! Counter of nominal solutions found
   DO k=1,nmc
      !********************************!
      !  Create one instance of noise  !
      !********************************!
      DO kk=1,m
         IF(obs(kk)%type.EQ.'O')THEN
            obs1(kk)=obs(kk)
            obs1(kk)%coord(1) = obs1(kk)%coord(1)+rvnorm()*obsw(kk)%rms_coord(1)/COS(obs1(kk)%coord(2))
            obs1(kk)%coord(2) = obs1(kk)%coord(2)+rvnorm()*obsw(kk)%rms_coord(2)
            obsw1(kk)=obsw(kk)
         ENDIF
      ENDDO
      !****************************!
      !  Compute nominal solution  !
      !****************************!
      icor=1
      CALL diff_cor(m,obs1,obsw1,el0,icor,-1,elk,unck,csinok,delnok,succ,nd)
      ! Check convergence
      IF(succ)THEN
         isucc=isucc+1
         ! First output the data point on the (rho,rho_dot) plane
         CALL coo_cha(elk,'ATT',elatt,fail_flag,OBSCODE=obscodn)
         WRITE(iuncob,130) k,isucc,elatt%coord(5:6),0.d0,0.d0,csinok,0.d0,0.d0,0.d0, &
              &            1,delnok,elatt%h_mag,0.d0
130      FORMAT(I4,1X,I4,1X,4(F12.7,1X),1P,D12.5,1X,D12.5,1X,D12.5,1X,D12.5,0P,1X,I1,1X,1P,D10.3,0P,1X,F6.3,1X,1P,D12.5)
         !******************************!
         !  Propagate until final time  !
         !******************************!
         CALL cov_not_av
         njc=0
         CALL pro_ele(elk,tf_search,el1)
         kill_propag=.FALSE.
         !******************************************!
         !  Read close approach data, if available  !
         !******************************************!
         IF(njc.GT.0)THEN
            IF(mtp_store(1)%iplam.EQ.3)THEN
               csi=mtp_store(1)%tp_coord(1)/reau
               zeta=mtp_store(1)%tp_coord(2)/reau
               ! Impactor?
               IF(csi**2+zeta**2.LT.1.d0)THEN
                  nimp=nimp+1
               ENDIF
            ENDIF
         ENDIF
      ENDIF
      WRITE(*,*) k, isucc, nimp, csinok
   ENDDO
   CALL filclo(iuncob,' ')
   verb_dif=20

 END SUBROUTINE mc_res_att

END MODULE cobweb

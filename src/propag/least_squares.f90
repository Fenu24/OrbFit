!***********************************************!
!                                               !
!             MODULE LEAST_SQUARES              !
!                                               !
!***********************************************!
! PUBLIC ROUTINES
!                errmod_set
!                fdiff_cor
!                diff_cor
!                constr_fit
!                step_fit (to be removed, replaced by step_fit2.multiple_sol)
!                sort_obs
!                unsort_obs
!                min_sol
!                blockset
!                fit_weight
!                mag_est
!                write_photom
!                difini
!                rejini
!                rms_compute
! NOT IN MODULE
!                meanti
!                bizarre
!                bizset
!                rescov
!                outcov
!                weak_dir
!                  scale_coef
! PRIVATE ROUTINES
!                blockdef
!                blocomp
!                invmat
!                sin_cor
!                reject_obs
!                   entopp
!
! MODULES AND HEADERS
! least_squares.o: \
!	../include/parobx.h90 \
!	../suit/astrometric_observations.mod \
!	../suit/fund_const.mod \
!	../suit/name_rules.mod \
!	../suit/output_control.mod \
!	../suit/util_suit.mod \
!	close_app.o \
!	least_squares.o \
!	obs_correl.o \
!	orbit_elements.o \
!	pred_obs.o

MODULE least_squares
  USE fund_const, ONLY: degrad
  USE astrometric_observations
  USE orbit_elements
  USE output_control
  USE dyn_param
  USE ades

  IMPLICIT NONE
  PRIVATE

  ! PUBLIC SUBROUTINEs
  PUBLIC :: fdiff_cor, diff_cor, sort_obs, unsort_obs, fit_weight ,difini, rejini
  PUBLIC :: mag_est, min_sol, rms_compute, write_photom, constr_fit, fourdim_fit
  PUBLIC :: blockset, step_fit, sin_cor

  ! PARAMETERS
  CHARACTER*20     :: error_model_priv                          ! Error model private string
  LOGICAL          :: errmod_given=.false.                      ! Is available
  LOGICAL          :: output_failrwo                            ! If true, write the rwo even upon failure
  DOUBLE PRECISION :: divrat,delcr, del_constr, step_sig
  INTEGER          :: iicdif                                    ! Initialization check
  LOGICAL          :: scaling_lov                               ! For weakdir, true if scaling required
  LOGICAL          :: second_lov                                ! For weakdir, selecting the second eigenvalue
  LOGICAL          :: aprc_flag                                 ! Controls full a priori normal matrix
  LOGICAL          :: output_des                                ! To use DES in output from fdiff_cor
  INTEGER          :: itmax                                     ! Max no iterations
  INTEGER          :: itgmax                                    ! Max iter. w. paralysed Q
  INTEGER          :: interactive                               ! Interactive: 0=automatic, 1=int. choice of solve-for parameters
  ! Options for outlier rejection/recovery
  LOGICAL          :: autrej                                    ! Automatic rejection flag PUBLIC
  LOGICAL          :: rejopp                                    ! Permit the rejection of an entire opp. PUBLIC
  LOGICAL          :: rej_fudge                                 ! Fudge factor active
  ! A priori
  LOGICAL          :: drpa2m_apriori, yarko_apriori             ! Use apriori on area/mass, Yarkovsky (def. F)
  LOGICAL          :: a1ng_apriori, a2ng_apriori                ! Use apriori a1ng, a2ng
  LOGICAL          :: a3ng_apriori, dtdelay_apriori             ! Use apriori for a3ng, dtdelay
  DOUBLE PRECISION :: drpa2m_nominal, yarko_nominal             ! A priori value (def. 0)
  DOUBLE PRECISION :: a1ng_nominal, a2ng_nominal                ! A priori value (def. 0)
  DOUBLE PRECISION :: a3ng_nominal, dtdelay_nominal             ! A priori value (def. 0)
  DOUBLE PRECISION :: drpa2m_std, yarko_std                     ! A priori standard deviation (no default)
  DOUBLE PRECISION :: a1ng_std, a2ng_std, a3ng_std, dtdelay_std ! A priori standard deviation (no default)

  PUBLIC :: errmod_set
  PUBLIC :: output_failrwo,autrej,rejopp,rej_fudge,itmax,itgmax,interactive,divrat,delcr,del_constr,step_sig
  PUBLIC :: scaling_lov
  PUBLIC :: second_lov
  PUBLIC :: output_des
  PUBLIC :: yarko_apriori, yarko_nominal, yarko_std, drpa2m_apriori, drpa2m_nominal, drpa2m_std, &
       & a1ng_apriori, a1ng_nominal, a1ng_std, a2ng_apriori, a2ng_nominal, a2ng_std, &
       & a3ng_apriori, a3ng_nominal, a3ng_std, dtdelay_apriori, dtdelay_nominal, dtdelay_std

  ! LIST OF PRIVATE ENTITIES, COMMON TO THE MODULE
  ! =============================================
  ! former CODODE.H requires parobx.h
  INCLUDE 'parobx.h90'
  DOUBLE PRECISION g(nob2x,ndimx) ! first derivatives of observations with respect to elements
  ! INTERFACE WITH blocks etc.
  ! block variables: max no blocks, max no residuals in block
  ! block definition:
  INTEGER stablk(nblx),nblk ! station code (integer) of the block, no. blocks
  DOUBLE PRECISION tblock(2,nblx) ! time interval of the block
  LOGICAL fullblk(nblx) ! block full, make another
  ! block variables: max no residuals in block
  ! ==== large arrays made allocatable==================
  INTEGER, ALLOCATABLE ::  indblk(:,:),noblk(:) ! pointers, no obs in blk
  INTEGER, ALLOCATABLE :: indblk_s(:,:) ! pointers to observations sorted in propagator order
  !INTEGER indblk(nxinbl_large,nblx),noblk(nblx) ! pointers, no obs in blk
  !INTEGER indblk_s(nxinbl_large,nblx) ! pointers to observations sorted in propagator order
  DOUBLE PRECISION, ALLOCATABLE ::  wblk(:,:,:) ! weight matrix for the bloc
  DOUBLE PRECISION, ALLOCATABLE ::  wblk_large(:,:,:) ! same for large blocks
  !DOUBLE PRECISION wblk(nxinbl,nxinbl,nblx) ! weight matrix for the block
  !DOUBLE PRECISION wblk_large(nxinbl_large,nxinbl_large,nbl_large) ! same for large blocks
  INTEGER large_blk(nblx)

  ! former comrej.h
  ! OPTIONS FOR OUTLIER REJECTION/RECOVERY
  !
  ! x2rej       -  Chi-square threshold for rejection
  ! x2rec       -  Chi-square threshold for recovery
  ! x2frac      -  Fraction of max chi-square at which start rejection
  !                in the "large residual" regime
  ! delrej      -  Convergency control (correction norm)
  ! fomax       -  Max fraction of outliers (total)
  ! itmaxr      -  Max number of iterations
  !
  DOUBLE PRECISION x2rej,x2rec,x2frac,delrej,fomax,x2mrej,x2mrec
  INTEGER itmaxr,iicrej

  ! ===============================================================
  ! data for mag_est, former mag headerh
  ! ===============================================================
  ! difference apparent minus absolute magnitude (to compute magnitude)
  ! dmagn is in original observations order,
  ! dmagns is sorted as required by propagator PUBLIC
  ! warning this header requires parobx.h
  ! also slope parameter known for current object
  ! phav,phavs is phase
  ! dsunv,dsunvs distance from the sun
  ! disv,disvs distance from Earth
  ! adotmv,adots ddotmv,ddots proper motion
  DOUBLE PRECISION dmagn(nobx),dmagns(nobx),gmagc
  DOUBLE PRECISION dsuna(nobx),dsunas(nobx),disa(nobx),disas(nobx)    &
       &  ,phaa(nobx),phaas(nobx),adots(nobx),ddots(nobx),                  &
       &   adotmv(nobx),ddotmv(nobx),elos(nobx),elomv(nobx)
  PUBLIC :: dmagns,elomv,adotmv,ddotmv

  LOGICAL, PUBLIC:: ab_mag ! are the magnitudes, as reported, AB magnitudes?
  !
  !
CONTAINS
  ! ===============================================
  ! errmod_set
  SUBROUTINE errmod_set(error_model)
    CHARACTER*(20) :: error_model ! error model file name
    error_model_priv=error_model
    errmod_given=.true.
  END SUBROUTINE errmod_set
  !
  ! Copyright (C) 1998 by OrbFit Consortium
  ! Version: December 15, 1997 Steven Chesley
  ! =========================================
  ! FDIFF_COR
  ! =========================================
  ! Differential corrections interface routine for FITOBS
  ! =========INTERFACE=============================
  SUBROUTINE fdiff_cor(batch,iarc,obs0,ini0,ok,cov0,el0,m,obs,obsw,nused,    &
       &     rwofi0,elc,uncert,csino0,delno0,rmsh,succ,rmshout,nmag,h_read,g_read,ades_type)
    USE util_suit
    !====== INTERFACE ============================================================================
    LOGICAL,          INTENT(IN)           :: batch          ! batch/interactive
    INTEGER,          INTENT(IN)           :: iarc           ! arc: 1,2 arcs, 3 together
    !                                                               <0 to avoid rwo and ades output
    LOGICAL,          INTENT(IN)           :: obs0,ini0      ! data availability
    LOGICAL,          INTENT(OUT)          :: ok             ! data available
    LOGICAL,          INTENT(OUT)          :: cov0           ! matrices available
    TYPE(orbit_elem), INTENT(INOUT)        :: el0            ! initial and corrected orb. el.
    INTEGER,          INTENT(IN)           :: m              ! observation number
    TYPE(ast_obs),    INTENT(IN)           :: obs(m)         ! observations
    TYPE(ast_wbsr),   INTENT(INOUT)        :: obsw(m)        ! weights/bias
    INTEGER,          INTENT(INOUT)        :: nused          ! no. scalar observations used (2 for optical, 1 for radar)
    CHARACTER(LEN=60),INTENT(IN)           :: rwofi0         ! residuals and weights file name
    TYPE(orbit_elem), INTENT(INOUT)        :: elc            ! corrected el., only OUT, but could be el0=elc
    !                                                          cannot be INTENT(OUT): eqc could have some address as eq0
    !                                                          WARNING: if accepted, they are also copied in eq0
    TYPE(orb_uncert), INTENT(OUT)          :: uncert         ! normal and covariance matrix (WARNING: not full if ncor.LT.6)
    DOUBLE PRECISION, INTENT(OUT)          :: csino0         ! norm of residuals
    DOUBLE PRECISION, INTENT(OUT)          :: delno0         ! norm of last correction
    DOUBLE PRECISION, INTENT(OUT)          :: rmsh           ! RMS of magnitudes
    LOGICAL,          INTENT(OUT)          :: succ           ! differential correction convergent
    DOUBLE PRECISION, INTENT(OUT),OPTIONAL :: rmshout        ! non-normalized RMS residual
    INTEGER,          INTENT(OUT),OPTIONAL :: nmag           ! no. magnitude obs
    DOUBLE PRECISION, INTENT(IN), OPTIONAL :: h_read, g_read ! H and G for correction of mag_est
    INTEGER,          INTENT(IN), OPTIONAL :: ades_type      ! flag for type of output file
    !====== END INTERFACE ========================================================================
    ! no. solve for variables, list, pseudoNewton/Newton
    INTEGER ncor,icor(ndimx),inew, isum
    ! rms of residuals
    DOUBLE PRECISION rms,sigma, rmsh1
    ! for manual discard., currently disabled
    ! discarded observations?
    !     LOGICAL disc
    !     DOUBLE PRECISION ra,rd,recov,disca
    ! loop indexes, counters, modes, units
    INTEGER i,nobnew,idif,iunf
    ! characters for menu
    CHARACTER*20 menunam
    ! dimension of solve for parameters vector
    INTEGER nd
    INTEGER ades_type_loc

    IF(.NOT.PRESENT(ades_type)) THEN
       ades_type_loc = 0   ! if flag is not present, we set it to print only .rwo file
    ELSE
       ades_type_loc = ades_type
    ENDIF
    ! =====================================================================
    ! check availability of observations and initial condition
    CALL cheobs(obs0,ini0,ok)
    IF(.not.ok) RETURN
    ! check availability of JPL ephemerides and ET-UT table
    CALL chetim(obs(1)%time_tdt,obs(m)%time_tdt,ok)
    IF(.not.ok) RETURN
    ! check availability of error_model
    IF(.not.errmod_given)THEN
       WRITE(*,*)' fdiff_cor: missing error_model'
       STOP
    ENDIF
    ! =====================================================================
    ! ok, go on with arc
    IF(.not.batch.and.iarc.lt.0)THEN
       ! case to remove output; should be in batch
       WRITE(*,*)'fdiff_cor: this should not happen, iarc=',iarc
    ENDIF
    ! =====================================================================
    ! dimension of solve for parameters vector
    nd=6+nls
    ! choose method, how many elements to correct (interactively???)
    IF(batch)THEN
       ! leave autoreject control autrej as it is set
       ! solve for all elements anyway
       interactive=0
       CALL whicor(interactive,nd,icor,ncor,inew)
    ELSE
       CALL whicor(interactive,nd,icor,ncor,inew)
    ENDIF
    ! differential corrections
    IF(batch)THEN
       iunf=-abs(iun_log)
    ELSE
       iunf=abs(iun_log)
    ENDIF
    CALL diff_cor(m,obs,obsw,el0,icor,iunf,elc,uncert,csino0,delno0,succ,nd)
    IF(.not.batch) WRITE(*,*)' RMS of weighed residuals is ',csino0
    !**********************!
    !  Estimate magnitude  !
    !**********************!
    IF(succ)THEN
       CALL mag_est(m,obs,obsw,elc%h_mag,rmsh,rmsh1,isum)
       IF(PRESENT(rmshout).and.PRESENT(nmag))THEN
          rmshout=rmsh1
          nmag=isum
       ENDIF
       ! Check for corrections
       IF(elc%h_mag.EQ.-9.99d0)THEN
          IF(PRESENT(h_read).AND.PRESENT(g_read))THEN
             elc%h_mag=h_read
             elc%g_mag=g_read
          ELSE
             WRITE(*,*) ' fdiff_cor: WARNING! No corrections for magnitude available!'
          END IF
       END IF
    ELSE
       IF(verb_dif.GT.19) WRITE(*,*)' fdiff_cor: no magnitude w/o new orbit'
       rmsh=-1.d0
       DO i=1,m
          obsw(i)%resm_def=.FALSE.
       ENDDO
    ENDIF
    IF(iarc.ge.0.and.(output_failrwo.or.succ))THEN
       ! output of residuals, weights and observations file
       ! warning: when the orbit is hyperbolic, these are the residuals of
       ! iteration one.
       ! when the divergence is mild (e.g. target function continues to increase or
       ! too many iterations) should be (TO BE FIXED)
       IF(output_des)THEN
          !        CALL tra_res_out(iuntra,m,obs,s2n,iunres,obsw)
       ELSE
          ! this would be used if rmshout (non-normalized) has to be written in the .rwo file
          !        IF(PRESENT(rmshout))THEN
          !           CALL write_rwo(rwofi0,obs,obsw,m,error_model_priv,csino0,rmshout)
          !        ELSE
          IF(ades_type_loc.GE.0) CALL write_rwo(rwofi0,obs,obsw,m,error_model_priv,csino0,rmsh)
          !        ENDIF
          IF (ades_type_loc .NE. 0) THEN
             CALL convert_fitobs_ades(obs,obsw,.TRUE.,.TRUE.)
             CALL write_ades(psv_name(1:len_ades_name),error_model_priv,ades_type_loc,csino0,rmsh)
          ENDIF
       ENDIF
    ENDIF
    ! is covariance available?
    IF(.not.succ)THEN
       cov0=.false.
    ELSE
       ! covariance matrix
       IF(ncor.eq.nd)THEN
          cov0=.true.
       ELSE
          cov0=.false.
       ENDIF
       ! output covariance in .fga file, not in batch
       IF(.not.batch.and.verb_dif.gt.9.and.iun_covar.ge.0)THEN
          IF(iarc.eq.1)then
             WRITE(iun_covar,*) 'COVARIANCE MATRIX FOR FIRST ARC'
          ELSEIF(iarc.eq.2)then
             WRITE(iun_covar,*) 'COVARIANCE MATRIX FOR SECOND ARC'
          ELSEIF(iarc.eq.4)THEN
             WRITE(iun_covar,*)'COVARIANCE MATRIX FOR BOTH ARCS'
          ENDIF
          CALL outcov(iun_covar,nd,icor,uncert%g,uncert%c)
       ENDIF
       ! improved orbital elements are accepted
       el0=elc
       ! observations used
       nused=0
       DO i=1,m
          IF(obsw(i)%sel_coord.gt.0)THEN
             IF(obs(i)%type.eq.'O'.or.obs(i)%type.eq.'S')THEN
                nused=nused+2
             ELSEIF(obs(i)%type.eq.'R'.or.obs(i)%type.eq.'V')THEN
                nused=nused+1
             ENDIF
          ENDIF
       ENDDO
    ENDIF
  END SUBROUTINE fdiff_cor

  !====================================================================!
  ! FOURDIM_FIT                                                        !
  !====================================================================!
  ! Linearly constrained differential corrector, only for attributable !
  ! elements, fixed (r,rdot).                                          !
  !====================================================================!
  ! Copyright 1998-2004, The OrbFit Consortium                         !
  ! - Version 3.2: November 2004                                       !
  ! - Last Update: July 2016, A. Del Vigna                             !
  !--------------------------------------------------------------------!
  ! REMARK. The weights are constant in this routine.                  !
  !====================================================================!
  SUBROUTINE fourdim_fit(m,obs,obsw,el0,elc,nd,uncert4,uncert,csinor,delnor,rmsh,nused,succ,des_mat,scal_obs)
    USE pred_obs
    USE close_app, ONLY: kill_propag
    INCLUDE 'parobx.h90'
    !=======================================================================================================
    INTEGER,          INTENT(IN)            :: m                 ! Number of observations
    TYPE(ast_obs),    INTENT(IN)            :: obs(m)            ! Observations
    TYPE(ast_wbsr),   INTENT(INOUT)         :: obsw(m)           ! Weights and biases
    TYPE(orbit_elem), INTENT(INOUT)         :: el0               ! Epoch of initial elements (not changed)
    TYPE(orbit_elem), INTENT(INOUT)         :: elc               ! Corrected elements at the same time
    INTEGER,          INTENT(IN)            :: nd                ! Dimension of solve for parameters vector
    TYPE(orb_uncert), INTENT(OUT)           :: uncert4           ! Uncertainty matrices of the 4-dim fit
    TYPE(orb_uncert), INTENT(OUT)           :: uncert            ! Uncertainty matrices of the 6-dim fit
    DOUBLE PRECISION, INTENT(OUT)           :: delnor            ! Differential corrections norm
    DOUBLE PRECISION, INTENT(OUT)           :: csinor            ! Residuals norm
    DOUBLE PRECISION, INTENT(OUT)           :: rmsh              ! RMS of the phtometric fit
    INTEGER,          INTENT(OUT)           :: nused             ! Number of observations used
    LOGICAL,          INTENT(OUT)           :: succ              ! Success flag
    DOUBLE PRECISION, INTENT(OUT), OPTIONAL :: des_mat(nob2x,nd) ! Design matrix
    INTEGER,          INTENT(OUT), OPTIONAL :: scal_obs          ! Actual number of observations (number of residuals)
    !===== Computation of alpha and delta ==================================================================
    INTEGER          :: iposs(nobx)  ! Sorting indexes for obs_sort
    INTEGER          :: ider         ! Flag to control computation of derivatives
    DOUBLE PRECISION :: tauj         ! Observations time
    DOUBLE PRECISION :: alj          ! Computed alpha
    DOUBLE PRECISION :: dej          ! Computed delta
    INTEGER          :: iocj         ! Observatory code
    DOUBLE PRECISION :: dade(6)      ! Partial derivatives of alpha w.r.t. elements
    DOUBLE PRECISION :: ddde(6)      ! Partial derivatives of delta w.r.t. elements
    DOUBLE PRECISION :: drde(6)      ! Partial derivatives of the range w.r.t. elements
    DOUBLE PRECISION :: dvde(6)      ! Partial derivatives of the range rate w.r.t. elements
    DOUBLE PRECISION :: pos(3)       ! Geocentric position of the observer
    DOUBLE PRECISION :: vel(3)       ! Geocentric velocity of the observer
    DOUBLE PRECISION :: alobs        ! To preserve integrity of obs_s
    DOUBLE PRECISION :: adot         ! Time derivative of alpha
    DOUBLE PRECISION :: ddot         ! Time derivative of delta
    DOUBLE PRECISION :: pha          ! Phase angle
    DOUBLE PRECISION :: dis          ! Distance from the Earth
    DOUBLE PRECISION :: dsun         ! Distance from the Sun
    DOUBLE PRECISION :: elo          ! Elongation
    DOUBLE PRECISION :: gallat       ! Galactic latitude
    DOUBLE PRECISION :: appmag       ! Apparent magnitude
    TYPE(ast_obs)    :: obs_s(nobx)  ! Sorted observations
    TYPE(ast_wbsr)   :: obsw_s(nobx) ! Sorted weights/biases
    !===== Differential corrections ========================================================================
    INTEGER          :: it           ! Iteration of differential corrections
    INTEGER          :: itg          ! Iteration with paralysed target function
    INTEGER          :: no           ! Actual number of observations (number of residuals)
    DOUBLE PRECISION :: g(nob2x,6)   ! Opposite of the design matrix
    DOUBLE PRECISION :: gr(nob2x,6)  ! Reduced deisgn matrix for 4-dim. fit (gr(:,5:6)=0.d0)
    INTEGER          :: icor(ndimx)  ! Array for the 4-dim diff. corrections (icor(5:6)=0)
    INTEGER          :: icor6(ndimx) ! Full correction (but NOT applied) to obtain the final matrices
    INTEGER          :: nsolv        ! Number of solve-for parameters
    INTEGER          :: indp         ! Inversion failure flag
    DOUBLE PRECISION :: deq(6)       ! Correction
    DOUBLE PRECISION :: cond         ! Conditioning number
    DOUBLE PRECISION :: csino0       ! Norm of residuals at first iteration
    DOUBLE PRECISION :: csino1       ! Auxiliary variable for csino
    DOUBLE PRECISION :: mu           ! Scaling factor for the uncertainty matrices
    DOUBLE PRECISION :: ecc          ! Eccentricity for bizarre control
    DOUBLE PRECISION :: deqv(nd)     ! Correction from the full fit
    !===== Functions =======================================================================================
    DOUBLE PRECISION :: snorm        ! Norm induced by a matrix
    DOUBLE PRECISION :: pridif       ! Difference of two angles
    LOGICAL          :: bizarre      ! Control of bizarre orbits
    DOUBLE PRECISION :: rescov       ! Scaling factor for the uncertainty matrices
    !=======================================================================================================
    INTEGER          :: j            ! Loop index on the observations
    INTEGER          :: i            ! Loop index on the solve-for parameters
    LOGICAL          :: twobo        ! Control of two body approximation (must be false)
    INTEGER          :: iun          ! Unit for output
    !=======================================================================================================
    ! Unit for output
    iun = ABS(iun_log)
    ! Control on dimension
    IF(nd.NE.6)THEN
       WRITE(*,*) 'fourdim_fit: not ready for non-grav, nd = ',nd
       STOP
    ENDIF
    ! Only for attributable elements
    IF(el0%coo.NE.'ATT')THEN
       WRITE(*,*) 'fourdim_fit: applied to wrong coordinates, ',el0%coo
       STOP
    ENDIF
    !******************!
    !  Inizialization  !
    !******************!
    icor      = 0       ! Array for the 4-dim diff. corrections
    icor(1:4) = 1       ! Number of solve-for variables (only the first 4)
    icor6     = 1       ! Full correction (but NOT applied) to obtain the final matrices
    ider      = 1       ! Assign ider depending from inew
    twobo     = .FALSE. ! Full n-body
    ! Inizialization of the uncertainty matrices
    CALL undefined_orb_uncert(nd,uncert)
    CALL undefined_orb_uncert(nd,uncert4)
    !*************************************!
    !  Definition of blocks if necessary  !
    !*************************************!
    CALL blockset(m,obs,obsw)                          ! Sort of times and reordering of obs. and weights
    CALL sort_obs(el0%t,obs,obsw,m,iposs,obs_s,obsw_s) ! Initialization with starting value for elements
    elc=el0                                            ! Assign el0 to elc, not to change el0
    IF(verb_dif.GT.9)  WRITE(iun,*) 'starting values ', el0
    IF(verb_dif.GT.19) WRITE(*,*) 'starting values ', el0

    !*****************************************!
    !  LOOP ON ITERATIONS OF NEWTON'S METHOD  !
    !*****************************************!
    itg=0
    DO 60 it=1,itmax
       ! Compute observations and derivatives
       CALL set_restart(.TRUE.)
       no=0
       DO 61 j=1,m
          tauj=obs_s(j)%time_tdt
          iocj=obs_s(j)%obscod_i
          pos=obs_s(j)%obspos
          vel=obs_s(j)%obsvel
          IF(obs_s(j)%type.EQ.'O' .OR. obs_s(j)%type .EQ. 'S')THEN
             CALL alph_del(elc,tauj,iocj,pos,vel,ider,twobo,nd,alj,dej,dade,ddde, &
                  &           adot,ddot,pha,dis,dsun)
             IF(kill_propag)THEN
                IF(verb_dif.GT.9) WRITE(*,*) 'fourdim_fit: kill_propag'
                succ=.FALSE.
                kill_propag=.FALSE.
                RETURN
             ENDIF
             ! Compute magnitude difference (apparent - absolute)
             dmagns(j)=appmag(0.d0,elc%g_mag,dsun,dis,pha)
          ELSEIF(obs_s(j)%type.EQ.'R' .OR. obs_s(j)%type.EQ.'V')THEN
             CALL r_rdot(elc,tauj,iocj,obs_s(j)%tech,vel,pos,alj,dej,nd,drde,dvde,ider)
          ELSE
             WRITE(*,*) 'fourdim_fit: obs. type ',obs_s(j)%type, ' not known'
             STOP 'fourdim_fit: wrong obs. type'
          ENDIF
          CALL set_restart(.FALSE.)
          !****************************************************!
          !  Compute residuals, form matrix g = -d(csi)/d(eq)  !
          !****************************************************!
          IF(obs_s(j)%type.EQ.'O' .OR. obs_s(j)%type.EQ.'S')THEN
             alobs                  = obs_s(j)%coord(1)
             obsw_s(j)%res_coord(1) = pridif(alobs,alj)
             obsw_s(j)%res_coord(2) = obs_s(j)%coord(2)-dej
             g(2*j-1,1:nd)          = dade
             g(2*j,1:nd)            = ddde
             no=no+2
          ELSEIF(obs_s(j)%type.EQ.'R')THEN
             obsw_s(j)%res_coord(1) = obs_s(j)%coord(1)-alj
             obsw_s(j)%res_coord(2) = 0.d0
             g(2*j-1,1:nd)          = drde
             g(2*j,1:nd)            = 0.d0
             no                     = no+1
          ELSEIF(obs_s(j)%type.EQ.'V')THEN
             obsw_s(j)%res_coord(1) = 0.d0
             obsw_s(j)%res_coord(2) = obs_s(j)%coord(2)-dej
             g(2*j-1,1:nd)          = 0.d0
             g(2*j,1:nd)            = dvde
             no                     = no+1
          ELSE
             WRITE(*,*) 'fourdim_fit: obs%type= ',obs_s(j)%type, ' not known'
          ENDIF
          obsw_s(j)%resc_def = .TRUE.
61     ENDDO
        !no=2*m ! Only for optical
       !*****************************************************************!
       !  Reduced design matrix neglecting derivatives w.r.t. (r,r_dot)  !
       !*****************************************************************!
       gr(1:no,1:4) = g(1:no,1:4)
       gr(1:no,5:6) = 0.d0
       ! Store the design matrix and the number of scalar obs for output
       IF(PRESENT(des_mat) .AND. PRESENT(scal_obs))THEN
          scal_obs = no
          des_mat(1:no,1:nd)=-g(1:no,1:nd)
       END IF
       !*************************************!
       !  ONE DIFFERENTIAL CORRECTIONS STEP  !
       !*************************************!
       ! Compute solution of linear least squares
       CALL min_sol(obs_s,obsw_s,m,gr,icor,iun,uncert4%c(1:nd,1:nd),deq(1:nd),uncert4%g(1:nd,1:nd),&
            & csinor,indp,cond,nd)
       IF(indp.NE.0) GOTO 9
       IF(it.EQ.1) csino0=csinor
       !*******************!
       !  Update solution  !
       !*******************!
       elc%coord=elc%coord+deq
       ! Norm of the correction (the zeros do not matter)
       delnor=snorm(deq(1:nd),uncert4%c(1:nd,1:nd),nd,nd)
       IF(verb_dif.GE.9) WRITE(iun,200) it,csinor,delnor,elc%coord
       IF(verb_dif.GE.19) WRITE(*,200) it,csinor,delnor,elc%coord
200    FORMAT('>>> Iteration 4d-fit ',I3,', RMS residuals =',1P,D12.4,    &
            &        ', norm corr =',D12.4,', new elem values:'/0P,6F13.7/)
       ! Control against hyperbolic and bizarre orbits
       IF(bizarre(elc,ecc))THEN
          IF(verb_dif.GE.19) WRITE(*,*) 'fourdim_fit: iter. ',it,' bizarre; e =',ecc
          IF(verb_dif.GE.9) WRITE(iun,*) 'fourdim_fit: iter. ',it,' bizarre; e =',ecc
          succ=.FALSE.
          csinor=csino0
          RETURN
       ENDIF
       !**************************************!
       !  Check if we need another iteration  !
       !**************************************!
       IF(delnor.LT.delcr)THEN
          IF(verb_dif.GE.9) WRITE(iun,*) 'fourdim_fit: convergence corrections small'
          IF(verb_dif.GE.19) WRITE(*,*) 'fourdim_fit: convergence corrections small'
          succ=.TRUE.
          GOTO 70
       ENDIF
       IF(it.GT.1)THEN
          IF(csinor.GT.csino1*1.1d0)THEN
             itg=itg+1
             IF(verb_dif.GE.9)WRITE(iun,*) 'fourdim_fit: target function increasing '
             IF(verb_dif.GE.19)WRITE(*,*) 'fourdim_fit: target function increasing '
             succ=.FALSE.
          ELSEIF(csinor.GT.csino1*divrat)THEN
             succ=.TRUE.
             IF(verb_dif.GE.9)WRITE(iun,*) 'fourdim_fit: target function paralyzed '
             IF(verb_dif.GE.19)WRITE(*,*) 'fourdim_fit: target function paralyzed '
             itg=itg+1
          ENDIF
          IF(itg.GT.itgmax)THEN
             IF(succ)THEN
                GOTO 70 ! Store residuals, normalize covariance, etc.
             ELSE
                csinor=csino0
                RETURN  ! Leave residuals from first iteration
             ENDIF
          ENDIF
       ENDIF
       csino1=csinor
60  ENDDO
70  CONTINUE ! Success achieved
    IF(.NOT.succ)THEN
       IF(verb_dif.GE.9)THEN
          WRITE(*,*) 'fourdim_fit: non convergent ',it,itg
       ELSEIF(verb_dif.GT.2)THEN
          WRITE(iun,*) 'fourdim_fit: non convergent ',it,itg
       ENDIF
    ENDIF

    !############################################################!
    !  COVARIANCE (AND NORMAL MATRIX) COMPUTATION AND RESCALING  !
    !############################################################!
    CALL min_sol(obs_s,obsw_s,m,g,icor6,iun,uncert%c(1:nd,1:nd),deqv,uncert%g(1:nd,1:nd),csinor,indp,cond,nd)
    IF(indp.EQ.0) uncert%succ=succ
    uncert4%succ=succ
    uncert%ndim=nd ! If using this routine nd=6
    uncert4%ndim=nd
    nsolv=nd
    !****************************************************!
    !  Compute the number of observations actually used  !
    !****************************************************!
    nused=0
    DO i=1,m
       IF(obsw_s(i)%sel_coord.GT.0)THEN
          IF(obs_s(i)%type.EQ.'O' .OR. obs_s(i)%type.EQ.'S')THEN
             nused=nused+2
          ELSEIF(obs_s(i)%type.EQ.'R' .OR. obs_s(i)%type.EQ.'V')THEN
             nused=nused+1
          ENDIF
       ENDIF
    ENDDO
    !****************************!
    !  Rescal covariance matrix  !
    !****************************!
    mu=rescov(nsolv,nused,csinor)
    ! Apply rescaling to both covariance and normal matrix
    DO  i=1,6
       DO  j=1,6
          uncert%g(i,j)=uncert%g(i,j)*mu**2
          uncert%c(i,j)=uncert%c(i,j)/mu**2
          uncert4%g(i,j)=uncert4%g(i,j)*mu**2
          uncert4%c(i,j)=uncert4%c(i,j)/mu**2
       ENDDO
    ENDDO

    !#######################!
    !  OUTPUT FINAL RESULT  !
    !#######################!
    IF(verb_dif.GE.9 .AND. verb_dif.LT.19) WRITE(*,201) it, csinor, delnor
    IF(verb_dif.LT.9 .AND. verb_dif.GT.2) WRITE(iun,201) it, csinor, delnor
201 FORMAT('Done 4-dim fit iter. ',I3,'  RMS =',1P,D12.4,' last corr. norm =',D12.4)
    ! Re-ordering the residuals for output
    CALL unsort_obs(iposs,m,obsw_s,obsw)
    ! Unsort dmagn
    DO i=1,m
       j=iposs(i)
       dmagn(j)=dmagns(i)
    ENDDO
    iun_log=iun
    !******************************!
    !  Compute absolute magnitude  !
    !******************************!
    CALL mag_est(m,obs,obsw,elc%h_mag,rmsh)
    elc%mag_set=.TRUE.
    RETURN

    ! Error case: inversion not possible
9   IF(verb_dif.GT.9)THEN
       WRITE(*,*) 'fourdim_fit: inversion failed, at row ', indp
       WRITE(iun,*) 'fourdim_fit: inversion failed, at row ', indp
    ENDIF
    succ=.FALSE.
  END SUBROUTINE fourdim_fit


  ! Copyright 2006, The Orbfit Consortium
  ! ===================================================================
  ! STEP_FIT differential corrector slippping along the LOV by steps
  ! WARNING: this version abandoned, used step_fit2 in multiple_sol.mod
  ! ===================================================================
  ! version 3.3.2, 6 January 2006
  !
  ! Input: m observations number
  !        obs observations
  !        obsw weights/bias
  !        el0 asteroid elements (NOT CHANGED)
  ! Output elc corrected orbital elements at same time
  !        uncert includes covariance matrix and inverse
  !        csinor residuals norm
  !        delnor differential corrections norm
  !        succ success flag; the solution is meaningful if succ=.true.
  ! ============= REMARK ===============================================
  ! The weights are constant in this routine
  ! =============INTERFACE===== =========================================
  SUBROUTINE step_fit(m,obs,obsw,el0,elc,uncert,csinor,delnor,rmsh,nused,succ,nd)
    ! =====================================================================
    USE pred_obs
    USE close_app, ONLY: kill_propag
    ! ================input data==========================
    INTEGER, INTENT(IN) ::  m ! number of observations
    TYPE(ast_obs), DIMENSION(m), INTENT(IN)  :: obs
    TYPE(ast_wbsr), DIMENSION(m), INTENT(INOUT) :: obsw
    TYPE(orbit_elem), INTENT(INOUT) :: el0  ! epoch time, initial elements
    INTEGER, INTENT(IN) :: nd ! dimension of solve for parameters vector
    ! ================output ==========================
    TYPE(orbit_elem), INTENT(INOUT) :: elc  ! corrected elements
    TYPE(orb_uncert), INTENT(OUT)  :: uncert ! normal and covar. matrix
    DOUBLE PRECISION, INTENT(OUT) :: delnor,csinor,rmsh ! corr,res norm
    INTEGER nused ! no obs. used
    LOGICAL succ ! success flag
    ! =============END INTERFACE============================================
    DOUBLE PRECISION :: peq(6) !  orthogonal vector
    TYPE(orbit_elem) elb
    ! proper motion, data to compute magnitude
    DOUBLE PRECISION :: adot,ddot,pha,dis,dsun,elo,gallat, appmag
    ! input data sorted
    INCLUDE 'parobx.h90'
    TYPE(ast_obs),DIMENSION(nobx) :: obs_s
    TYPE(ast_wbsr),DIMENSION(nobx) :: obsw_s
    integer iposs(nobx),iocj,no
    DOUBLE PRECISION tauj,alj,dej
    ! condition number, sigmna along weak dir
    DOUBLE PRECISION cond, sdir
    INTEGER indp ! inversion failure flag
    ! ====================================================================
    ! partial derivatives of observables w.r. to elements
    DOUBLE PRECISION dade(6),ddde(6),drde(6),dvde(6)
    ! position and velocity (geocentric) of observer
    DOUBLE PRECISION pos(3), vel(3)
    ! to preserve integrity of obs_s
    DOUBLE PRECISION alobs
    ! first derivatives of observations
    DOUBLE PRECISION g(nob2x,6),gr(nob2x,6),gs(nob2x,6)
    !  corr. and residuals norm, their controls
    integer ider,icor(ndimx),icor6(ndimx),nsolv
    DOUBLE PRECISION csino0,csino1,mu,rescov,gam(6,6),c(6,6)
    ! differential correction, eccentricity
    DOUBLE PRECISION deq(6),deq6(6),deq6_lov(6),deq6_const(6),ecc
    LOGICAL lov_step
    ! iteration indexes
    integer it,itg,itm
    ! matrices for coordinate change
    DOUBLE PRECISION v(6,6),deqv(6)
    ! loop indexes: j=1,m; i,k=1,6 , ibk for blocks
    integer j,k,i, ibk
    ! DOUBLE PRECISION functions
    DOUBLE PRECISION snorm,snormd, pridif
    logical twobo ! control of two body approximation (must be false)
    LOGICAL bizarre ! function for control of bizarre orbit
    INTEGER ir ! index for relaxation loop, to avoid bizarre corrections
    ! unit for output
    integer iun
    ! scaling LOV
    DOUBLE PRECISION, DIMENSION(6) :: units
    ! ====================================================================
    ! control on dimension
    IF(nd.ne.6)THEN
       WRITE(*,*)'step-fit not to be used with dyn_param, nd=',nd
       STOP
    ENDIF
    ! number of solve-for variables
    icor(1)=0
    icor(2:nd)=1
    icor6=1 ! for final matrices
    ! ====================================================================
    ! assign ider depending from inew
    ider=1
    ! full n-body
    twobo=.false.
    ! definition of blocks if necessary
    CALL blockset(m,obs,obsw)
    ! sort of times and reordering of obs. and weights
    CALL sort_obs(el0%t,obs,obsw,m,iposs,obs_s,obsw_s)
    ! Initialisation with starting value for elements
    elc=el0
    ! ====================================================================
    iun=abs(iun_log)
    if(verb_dif.gt.9)write(iun,*)'starting values ', el0
    if(verb_dif.gt.19) write(*,*)'starting values ',  el0
    ! ================== main loop ==============================
    ! iteration control parameters: use default
    ! assign ider
    ider=1
    ! Loop on iterations of differential corrections, alternate
    itg=0
    itm=3*itmax
    do 60 it=1,itm
       ! ================== iteration ==============================
       ! Compute observations and derivatives
       CALL set_restart(.true.)
       do 61 j=1,m
          tauj=obs_s(j)%time_tdt
          iocj=obs_s(j)%obscod_i
          pos=obs_s(j)%obspos
          vel=obs_s(j)%obsvel
          IF(obs_s(j)%type.eq.'O'.or.obs_s(j)%type.eq.'S')THEN
             CALL alph_del(elc,tauj,iocj,pos,vel,ider,twobo,nd,alj,dej,dade,ddde, &
                  &           adot,ddot,pha,dis,dsun)
             IF(kill_propag)THEN
                IF(verb_dif.gt.9) WRITE(*,*)' step_fit: kill_propag'
                succ=.false.
                kill_propag=.false.
                RETURN
             ENDIF
             !  compute magnitude difference (apparent minus absolute)
             dmagns(j)=appmag(0.d0,elc%g_mag,dsun,dis,pha)
          ELSEIF(obs_s(j)%type.eq.'R'.or.obs_s(j)%type.eq.'V')THEN
             CALL r_rdot(elc,tauj,iocj,obs_s(j)%tech,vel,pos,alj,dej,nd,drde,dvde,ider)
          ELSE
             WRITE(*,*)'step_fit: obs. type ',obs_s(j)%type, ' not known'
             STOP 'step_fit: obstype'
          ENDIF
          CALL set_restart(.false.)
          ! Compute residuals, form matrix g=-d(csi)/d(eq)
          IF(obs_s(j)%type.eq.'O'.or.obs_s(j)%type.eq.'S')THEN
             alobs=obs_s(j)%coord(1)
             obsw_s(j)%res_coord(1)=pridif(alobs,alj)
             obsw_s(j)%res_coord(2)=obs_s(j)%coord(2)-dej
             g(2*j-1,1:6)=dade
             g(2*j,1:6)=ddde
          ELSEIF(obs_s(j)%type.eq.'R')THEN
             obsw_s(j)%res_coord(1)=obs_s(j)%coord(1)-alj
             obsw_s(j)%res_coord(2)=0.d0
             g(2*j-1,1:6)=drde
             g(2*j,1:6)=0.d0
          ELSEIF(obs_s(j)%type.eq.'V')THEN
             obsw_s(j)%res_coord(1)=0.d0
             obsw_s(j)%res_coord(2)=obs_s(j)%coord(2)-dej
             g(2*j-1,1:6)=0.d0
             g(2*j,1:6)=dvde
          ELSE
             WRITE(*,*)'step_fit: obs%type= ',obs_s(j)%type, ' not known'
          ENDIF
          obsw_s(j)%resc_def=.true.
61     ENDDO
       no=2*m
       ! ===========linear constraint=============================
       CALL min_sol(obs_s,obsw_s,m,g,icor6,iun,                         &
            &                uncert%c,deq6,uncert%g,csinor,indp,cond,nd)
       IF(indp.ne.0) GOTO 9
       CALL weak_dir(uncert%g,peq,sdir,-1,elc%coo,elc%coord,units,nd)
       ! scaling of matrix d (Xi)/dX
       IF(scaling_lov)THEN
          DO i=1,no
             DO j=1,6
                gs(i,j)=g(i,j)*units(j)
             ENDDO
          ENDDO
       ELSE
          units=1.d0
          gs(1:no,1:6)=g(1:no,1:6)
       ENDIF
       deq6=deq6/units ! diff. corr. in scaled coords
       IF(DOT_PRODUCT(peq,deq6).lt.0.d0) peq=-peq ! sign rule to decrease Q
       deq6_lov=DOT_PRODUCT(peq,deq6)*peq ! component along the LOV
       deq6_const=deq6-deq6_lov ! component orthogonal to the LOV
       deq6_lov=deq6_lov*units ! in unscaled coordinates
       deq6_const=deq6_const*units ! in unscaled coordinates
       !******can be removed after check************************
       ! find orthonormal basis with peq as first vector
       CALL graha_1(peq,6,v)
       ! convert derivatives
       gr(1:no,1:6)=MATMUL(gs(1:no,1:6),v)
       ! ===========one differential corrections step=================
       ! Compute solution of linear least squares
       call min_sol(obs_s,obsw_s,m,gr,icor,iun,c,deqv,gam,csinor,indp,cond,nd)
       IF(it.eq.1)csino0=csinor
       IF(indp.ne.0) GOTO 9
       ! norm of the constrained correction: the zeros do not matter!
       delnor=snorm(deqv,c,6,6)
       !***************end removal**********************************
       !     delnor=snorm(deq6_const,uncert%c,6,6)
       ! Check if we need another constrained iteration
       IF(delnor.ge.del_constr)THEN
          IF(verb_dif.ge.19)write(*,*)' constr.corr., delnor=',delnor
          lov_step=.false.
          ! use constrained corrections
          !        deq=deq6_const
          !******can be removed after check************************
          ! convert back correction to the previous base
          deq=MATMUL(v,deqv)
          ! convert correction to unscaled coordinates
          deq=deq*units
          !****************end removal*****************************
       ELSE
          succ=.true.
          ! use unconstrained corrections, of limited length in the weak direction
          delnor=snorm(deq6_lov,uncert%c,6,6)
          IF(verb_dif.ge.19)write(*,*)' full corr., delnor=',delnor
          lov_step=.true.
          IF(delnor.gt.step_sig)THEN
             deq=deq6_const+deq6_lov*(step_sig/delnor)
          ELSE
             deq=deq6_const+deq6_lov
          ENDIF
       ENDIF
       ! Update solution
       !     elc%coord=elc%coord+deq
       ! relaxation loop
       elb=elc
       DO ir=1,4
          ! Update solution
          elb%coord=elc%coord+deq
          IF(bizarre(elb,ecc))THEN
             IF(verb_dif.ge.10)WRITE(*,*)' step_fit: short step to avoid ecc=',ecc
             deq=deq/2.d0
             elb%coord=elc%coord+deq
          ELSE
             EXIT
          ENDIF
       ENDDO
       elc=elb
       ! norm of the correction: the zeros do not matter!
       IF(verb_dif.ge.9)write(iun,200)it,csinor,delnor,elc%coord
       IF(verb_dif.ge.19)write(*,200)it,csinor,delnor,elc%coord
200    format(' *** iteration ',i3,' RMS residuals =',1p,d12.4,    &
            &        '   norm corr =',d12.4,'  new elem values:'/0p,6f13.7/)
       ! control against hyperbolic and bizarre orbits
       IF(bizarre(elc,ecc))THEN
          IF(verb_dif.ge.19)write(*,*)' step_fit: iter. ',it,' bizarre; e=',ecc
          IF(verb_dif.ge.9) write(iun,*)' step_fit: iter. ',it,' bizarre; e=',ecc
          succ=.false.
          csinor=csino0
          RETURN
       ENDIF
       ! Check if we need another iteration
       IF(lov_step)THEN
          if(delnor.lt.delcr)then
             IF(verb_dif.ge.9)write(iun,*)' convergence corrections small'
             IF(verb_dif.ge.19)write(*,*)' convergence corrections small'
             succ=.true.
             goto 70
          endif
          if(it.gt.1)then
             if(csinor.gt.csino1*1.1d0)then
                itg=itg+1
                IF(verb_dif.ge.9)write(iun,*)' target function increasing '
                IF(verb_dif.ge.19)write(*,*)' target function increasing '
                succ=.false.
             elseif(csinor.gt.csino1*divrat)then
                succ=.true.
                IF(verb_dif.ge.9)write(iun,*)' target function paralyzed '
                IF(verb_dif.ge.19)write(*,*)' target funct. paralyzed '
                itg=itg+1
             endif
             if(itg.gt.itgmax)then
                IF(succ)THEN
                   GOTO 70 ! store residuals, normalize covariance, etc.
                ELSE
                   csinor=csino0
                   RETURN  ! leave residuals from first iteration
                ENDIF
             endif
          endif
       ENDIF
       csino1=csinor
60  ENDDO
    !  succ=.true.
70  continue
    IF(.not.succ)THEN
       IF(verb_dif.ge.9)THEN
          WRITE(*,*)' non convergent ',it,itg
       ELSEIF(verb_dif.gt.2)THEN
          WRITE(iun,*)' non convergent ',it,itg
       ENDIF
    ENDIF
    ! =========covariance (and normal matrix) rescaling==============
    nsolv=6
    nused=0
    DO i=1,m
       IF(obsw_s(i)%sel_coord.gt.0)THEN
          IF(obs_s(i)%type.eq.'O'.or.obs_s(i)%type.eq.'S')THEN
             nused=nused+2
          ELSEIF(obs_s(i)%type.eq.'R'.or.obs_s(i)%type.eq.'V')THEN
             nused=nused+1
          ENDIF
       ENDIF
    ENDDO
    mu=rescov(nsolv,nused,csinor)
    ! apply rescaling to both covariance and normal matrix
    DO  i=1,6
       DO  j=1,6
          uncert%g(i,j)=uncert%g(i,j)*mu**2
          uncert%c(i,j)=uncert%c(i,j)/mu**2
       ENDDO
    ENDDO
    uncert%succ=succ
    uncert%ndim=6 ! in step-fit no nongrav
    ! Output final result, with norm of the residuals
    IF(verb_dif.ge.9.and.verb_dif.lt.19)write(*,201)it,csinor,delnor
201 format(' done step_fit iter. ',i3,'  RMS=',1p,d12.4,' last corr. =',d12.4)
    !  IF(verb_dif.lt.9.and.verb_dif.gt.2)write(iun,201)it,csinor,delnor
    ! reordering the residuals for output
    CALL unsort_obs(iposs,m,obsw_s,obsw)
    ! unsort  dmagn, dsunv,disv,phav,adotv,ddotv
    ! (data for magnitude computation)
    DO i=1,m
       j=iposs(i)
       dmagn(j)=dmagns(i)
    ENDDO
    ! compute magnitude
    iun_log=iun
    CALL mag_est(m,obs,obsw,elc%h_mag,rmsh)
    elc%mag_set=.true.
    RETURN
    ! error case; inversion is not possible....
9   IF(verb_dif.gt.9)THEN
       WRITE(*,*)' inversion failed, at row ', indp
       WRITE(iun,*)' inversion failed, at row ', indp
    ENDIF
    succ=.false.
    ! =====================
  END SUBROUTINE step_fit
  !
  ! Copyright 1998-2003, The Orbfit Consortium
  ! ====================================================================
  ! CONSTR_FIT  linearly constrained differential corrector
  ! ====================================================================
  ! Version 3.1, 7 June 2004, with scaling_lov/second_lov
  !
  ! Input: m observations number
  !        obs observations
  !        obsw weights/bias
  !        el0 asteroid elements (NOT CHANGED)
  ! Output elc corrected orbital elements at same time
  !        peq vector to which deq must be orthogonal; in input
  !            it is used to give sign
  !        uncert includes covariance matrix and inverse
  !        csinor residuals norm
  !        delnor differential corrections norm
  !        succ success flag; the solution is meaningful if succ=.true.
  ! ============= REMARK ===============================================
  ! The weights are constant in this routine
  ! ====================================================================
  SUBROUTINE constr_fit(m,obs,obsw,el0,peq,elc,uncert,csinor,delnor,rmsh,nused,succ,nd,h_read,g_read)
    USE pred_obs
    USE close_app, ONLY: kill_propag
    !====== INTERFACE ============================================================================
    INTEGER,          INTENT(IN)           :: m                  ! number of observations
    TYPE(ast_obs),    INTENT(IN)           :: obs(m)             ! observations
    TYPE(ast_wbsr),   INTENT(INOUT)        :: obsw(m)            ! weights/bias
    TYPE(orbit_elem), INTENT(INOUT)        :: el0                ! epoch time, initial elements
    INTEGER,          INTENT(IN)           :: nd                 ! dim of solve for parameters vector
    DOUBLE PRECISION, INTENT(IN), OPTIONAL :: h_read, g_read     ! H and G for correction of mag_est
    DOUBLE PRECISION, INTENT(INOUT)        :: peq(nd)            ! orthogonal vector
    TYPE(orbit_elem), INTENT(INOUT)        :: elc                ! corrected elements
    TYPE(orb_uncert), INTENT(OUT)          :: uncert             ! normal and covar. matrix
    DOUBLE PRECISION, INTENT(OUT)          :: delnor,csinor,rmsh ! corr,res norm
    INTEGER,          INTENT(OUT)          :: nused              ! no obs. used
    LOGICAL,          INTENT(OUT)          :: succ               ! success flag
    !====== END INTERFACE ========================================================================
    TYPE(orbit_elem) elb
    ! proper motion, data to compute magnitude
    DOUBLE PRECISION :: adot,ddot,pha,dis,dsun,elo,gallat, appmag
    ! input data sorted
    INCLUDE 'parobx.h90'
    TYPE(ast_obs),DIMENSION(nobx) :: obs_s
    TYPE(ast_wbsr),DIMENSION(nobx) :: obsw_s
    integer iposs(nobx),iocj,no
    DOUBLE PRECISION tauj,alj,dej
    ! condition number, sigmna along weak dir, check weak-dir
    DOUBLE PRECISION cond, sdir, peq0(ndimx)
    INTEGER indp ! inversion failure flag
    ! ====================================================================
    ! first and second derivatives of alpha, delta w.r. to elements
    DOUBLE PRECISION dade(ndimx),ddde(ndimx),drde(ndimx),dvde(ndimx)
    ! position and velocity (geocentric) of observer
    DOUBLE PRECISION pos(3), vel(3)
    ! to preserve integrity of obs_s
    DOUBLE PRECISION alobs
    ! first derivatives of observations
    DOUBLE PRECISION g(nob2x,ndimx),gr(nob2x,ndimx),gs(nob2x,ndimx)
    !  corr. and residuals norm, their controls
    integer ider,icor(ndimx),icor6(ndimx),nsolv
    DOUBLE PRECISION csino0,csino1,mu,rescov,gam(ndimx,ndimx),c(ndimx,ndimx)
    ! differential correction, eccentricity
    DOUBLE PRECISION deq(ndimx),ecc
    ! iteration indexes
    integer it,itg,itm
    ! matrices for coordinate change
    DOUBLE PRECISION v(ndimx,ndimx),deqv(ndimx)
    ! loop indexes: j=1,m; i,k=1,6 , ibk for blocks
    integer j,k,i, ibk ,ii
    ! DOUBLE PRECISION functions
    DOUBLE PRECISION snorm,snormd, pridif
    logical twobo ! control of two body approximation (must be false)
    LOGICAL bizarre ! function for control of bizarre orbit
    INTEGER ir ! index for relaxation loop, to avoid bizarre corrections
    ! unit for output
    integer iun
    ! scaling LOV
    DOUBLE PRECISION, DIMENSION(ndimx) :: units
    ! ====================================================================
    ! number of solve-for variables
    icor(1)=0
    icor(2:ndimx)=1
    icor6=1 ! for final matrices
    ! ====================================================================
    ! assign ider depending from inew
    ider=1
    ! full n-body
    twobo=.false.
    peq0(1:nd)=peq
    ! definition of blocks if necessary
    CALL blockset(m,obs,obsw)
    ! sort of times and reordering of obs. and weights
    CALL sort_obs(el0%t,obs,obsw,m,iposs,obs_s,obsw_s)
    ! Initialisation with starting value for elements
    elc=el0
    ! ====================================================================
    iun=abs(iun_log)
    if(verb_dif.gt.9)write(iun,*)'starting values ', el0
    if(verb_dif.gt.19) write(*,*)'starting values ',  el0
    ! ================== main loop ==============================
    ! iteration control parameters: use default
    ! assign ider
    ider=1
    ! Loop on iterations of NEWTON's METHOD
    itg=0
    itm=itmax
    do 60 it=1,itmax
       ! ================== iteration ==============================
       ! Compute observations and derivatives
       CALL set_restart(.true.)
       no=0 ! counter of scalar observations
       do 61 j=1,m
          tauj=obs_s(j)%time_tdt
          iocj=obs_s(j)%obscod_i
          pos=obs_s(j)%obspos
          vel=obs_s(j)%obsvel
          IF(obs_s(j)%type.eq.'O'.or.obs_s(j)%type.eq.'S')THEN
             CALL alph_del(elc,tauj,iocj,pos,vel,ider,twobo,nd,alj,dej,dade,ddde, &
                  &           adot,ddot,pha,dis,dsun)
             IF(kill_propag)THEN
                IF(verb_dif.gt.9) WRITE(*,*)' constr_fit: kill_propag'
                succ=.false.
                kill_propag=.false.
                RETURN
             ENDIF
             !  compute magnitude difference (apparent minus absolute)
             dmagns(j)=appmag(0.d0,elc%g_mag,dsun,dis,pha)
          ELSEIF(obs_s(j)%type.eq.'R'.or.obs_s(j)%type.eq.'V')THEN
             CALL r_rdot(elc,tauj,iocj,obs_s(j)%tech,vel,pos,alj,dej,nd,drde,dvde,ider)
          ELSE
             WRITE(*,*)'constr_fit: obs. type ',obs_s(j)%type, ' not known'
             STOP 'constr_fit: obstype'
          ENDIF
          CALL set_restart(.false.)
          ! Compute residuals, form matrix g=-d(csi)/d(eq)
          IF(obs_s(j)%type.eq.'O'.or.obs_s(j)%type.eq.'S')THEN
             alobs=obs_s(j)%coord(1)-obsw_s(j)%bias_coord(1)
             obsw_s(j)%res_coord(1)=pridif(alobs,alj)
             obsw_s(j)%res_coord(2)=obs_s(j)%coord(2)-dej-obsw_s(j)%bias_coord(2)
             g(2*j-1,1:nd)=dade(1:nd)
             g(2*j,1:nd)=ddde(1:nd)
             no=no+2
          ELSEIF(obs_s(j)%type.eq.'R')THEN
             obsw_s(j)%res_coord(1)=obs_s(j)%coord(1)-alj-obsw_s(j)%bias_coord(1)
             obsw_s(j)%res_coord(2)=0.d0
             g(2*j-1,1:nd)=drde(1:nd)
             g(2*j,1:nd)=0.d0
             no=no+1
          ELSEIF(obs_s(j)%type.eq.'V')THEN
             obsw_s(j)%res_coord(1)=0.d0
             obsw_s(j)%res_coord(2)=obs_s(j)%coord(2)-dej-obsw_s(j)%bias_coord(2)
             g(2*j-1,1:nd)=0.d0
             g(2*j,1:nd)=dvde(1:nd)
             no=no+1
          ELSE
             WRITE(*,*)'constr_fit: obs%type= ',obs_s(j)%type, ' not known'
          ENDIF
          obsw_s(j)%resc_def=.true.
61     ENDDO
       ! ===========linear constraint=============================
       uncert%g=0.d0
       uncert%c=0.d0
       CALL min_sol(obs_s,obsw_s,m,g,icor6,iun,                         &
            &  uncert%c(1:nd,1:nd),deqv,uncert%g(1:nd,1:nd),csinor,indp,cond,nd)
       IF(indp.ne.0) GOTO 9
       CALL weak_dir(uncert%g(1:nd,1:nd),peq,sdir,-1,elc%coo,elc%coord,units,nd)
       IF(DOT_PRODUCT(peq,peq0(1:nd)).lt.0.d0)peq=-peq
       peq0(1:nd)=peq
       ! scaling of matrix d (Xi)/dX
       IF(scaling_lov)THEN
          DO i=1,2*m
             DO j=1,nd
                gs(i,j)=g(i,j)*units(j)
             ENDDO
          ENDDO
       ELSE
          units=1.d0
          gs(1:2*m,1:nd)=g(1:2*m,1:nd)
       ENDIF
       ! find orthonormal basis with peq as first vector
       CALL graha_1(peq,nd,v(1:nd,1:nd))
       ! convert derivatives
       gr(1:2*m,1:nd)=MATMUL(gs(1:2*m,1:nd),v(1:nd,1:nd))
       ! ===========one differential corrections step=================
       ! Compute solution of linear least squares
       call min_sol(obs_s,obsw_s,m,gr,icor,iun,c(1:nd,1:nd),deqv,gam(1:nd,1:nd),csinor,indp,cond,nd)
       IF(verb_dif.gt.10)THEN
          DO ii=7,nd
             WRITE(*,*) 'gamma(',ii,',',ii,')=',gam(ii,ii)
             WRITE(*,*) 'sigma=',SQRT(ABS(gam(ii,ii)))
          ENDDO
       ENDIF
       IF(indp.ne.0) GOTO 9
       IF(it.eq.1)csino0=csinor
       ! convert back correction to the previous base
       deq(1:nd)=MATMUL(v(1:nd,1:nd),deqv(1:nd))
       ! convert correction to unscaled coordinates
       deq(1:nd)=deq(1:nd)*units(1:nd)
       ! Update solution:  relaxation loop
       elb=elc
       DO ir=1,4
          ! Update solution
          elb%coord=elc%coord+deq(1:6)
          IF(bizarre(elb,ecc))THEN
             IF(verb_dif.ge.9)WRITE(iun_log,*)' constr_fit: short step to avoid ecc=',ecc
             deq=deq/2.d0
             elb%coord=elc%coord+deq(1:6)
          ELSE
             EXIT
          ENDIF
       ENDDO
       elc=elb
       ! if there are nongrav parameters, correct also them
       IF(nd.gt.6)THEN
          DO k=7,nd
             dyn%dp(ls(k-6))=dyn%dp(ls(k-6))+deq(k)
          ENDDO
       ENDIF
       ! norm of the correction: the zeros do not matter!
       delnor=snorm(deqv,c(1:nd,1:nd),nd,nd)
       IF(verb_dif.ge.9)write(iun,200)it,csinor,delnor,elc%coord
       IF(verb_dif.ge.19)write(*,200)it,csinor,delnor,elc%coord
200    format(' constr_fit iteration ',i3,' RMS residuals =',1p,d12.4,    &
            &        '   norm corr =',d12.4,'  new elem values:'/0p,6f13.7/)
       IF(verb_dif.ge.19.and.nd.gt.6)THEN
          WRITE(*,*)'nongrav ',dyn%dp(ls(1:nd-6))
       ENDIF
       ! control against hyperbolic and bizarre orbits
       IF(bizarre(elc,ecc))THEN
          IF(verb_dif.ge.19)write(*,*)' constr_fit: iter. ',it,' bizarre; e=',ecc
          IF(verb_dif.ge.9) write(iun,*)' constr_fit: iter. ',it,' bizarre; e=',ecc
          succ=.false.
          csinor=csino0
          RETURN
       ENDIF
       ! Check if we need another iteration
       if(delnor.lt.delcr)then
          IF(verb_dif.ge.9)write(iun,*)' convergence corrections small'
          IF(verb_dif.ge.19)write(*,*)' convergence corrections small'
          succ=.TRUE.
          GOTO 70
       ENDIF
       IF(it.GT.1)THEN
          IF(csinor.GT.csino1*1.1d0)then
             itg=itg+1
             IF(verb_dif.GE.9) WRITE(iun,*)' target function increasing '
             IF(verb_dif.GE.19) WRITE(*,*)' target function increasing '
             succ=.FALSE.
          ELSEIF(csinor.GT.csino1*divrat)THEN
             succ=.TRUE.
             IF(verb_dif.GE.9) WRITE(iun,*)' target function paralyzed '
             IF(verb_dif.GE.19) WRITE(*,*)' target funct. paralyzed '
             itg=itg+1
          ENDIF
          IF(itg.GT.itgmax)THEN
             IF(succ)THEN
                GOTO 70 ! store residuals, normalize covariance, etc.
             ELSE
                csinor=csino0
                RETURN  ! leave residuals from first iteration
             ENDIF
          ENDIF
       ENDIF
       csino1=csinor
60  ENDDO
    ! succ=.true.
70  CONTINUE
    IF(.NOT.succ)THEN
       IF(verb_dif.GE.9)THEN
          WRITE(*,*)' non convergent ',it,itg
       ELSEIF(verb_dif.GT.2)THEN
          WRITE(iun,*)' non convergent ',it,itg
       ENDIF
    ENDIF
    ! =========covariance (and normal matrix) rescaling==============
    nsolv=nd
    nused=0
    DO i=1,m
       IF(obsw_s(i)%sel_coord.gt.0)THEN
          IF(obs_s(i)%type.EQ.'O'.OR.obs_s(i)%type.EQ.'S')THEN
             nused=nused+2
          ELSEIF(obs_s(i)%type.EQ.'R'.OR.obs_s(i)%type.EQ.'V')THEN
             nused=nused+1
          ENDIF
       ENDIF
    ENDDO
    mu=rescov(nsolv,nused,csinor)
    ! apply rescaling to both covariance and normal matrix
    DO  i=1,nd
       DO  j=1,nd
          uncert%g(i,j)=uncert%g(i,j)*mu**2
          uncert%c(i,j)=uncert%c(i,j)/mu**2
       ENDDO
    ENDDO
    uncert%succ=succ
    uncert%ndim=nd
    ! Output final result, with norm of the residuals
    IF(verb_dif.GE.9.AND.verb_dif.LT.19) WRITE(*,201)it,csinor,delnor
201 FORMAT(' done constr. iter. ',i3,'  RMS=',1p,d12.4,' last corr. =',d12.4)
    !  IF(verb_dif.lt.9.and.verb_dif.gt.2)write(iun,201)it,csinor,delnor
    ! reordering the residuals for output
    CALL unsort_obs(iposs,m,obsw_s,obsw)
    ! unsort  dmagn, dsunv,disv,phav,adotv,ddotv
    ! (data for magnitude computation)
    DO i=1,m
       j=iposs(i)
       dmagn(j)=dmagns(i)
    ENDDO
    !*********************!
    !  Compute magnitude  !
    !*********************!
    iun_log=iun
    CALL mag_est(m,obs,obsw,elc%h_mag,rmsh)
    ! Check for corrections
    IF(elc%h_mag.EQ.-9.99d0)THEN
       IF(PRESENT(h_read).AND.PRESENT(g_read))THEN
          elc%h_mag=h_read
          elc%g_mag=g_read
       !ELSE
       !   WRITE(*,*) ' constr_fit: WARNING! No corrections for magnitude available!'
       END IF
    END IF
    elc%mag_set=.TRUE.
    RETURN
    !**************************************!
    !  Error case: inversion not possible  !
    !**************************************!
9   IF(verb_dif.gt.9)THEN
       WRITE(*,*) ' constr_fit: inversion failed, at row ', indp
       WRITE(iun,*) ' constr_fit: inversion failed, at row ', indp
    ENDIF
    succ=.FALSE.
  END SUBROUTINE constr_fit
  !
  ! Copyright (C) 1998 by OrbFit Consortium
  ! ===================================================================
  ! DIFF_COR  differential corrector
  ! ===================================================================
  ! version 3.0 Andrea Milani, November 2002
  ! version 1.8 Steven Chesley, Dec. 15, 1998
  !        works on arbitrary  elements:
  ! Input: m observations number
  !        obs, obsw observations, weights, bias, etc
  !        el0 orbital elements, including epoch time
  !        icor(nd) flags .ne.0 to correct this element, 0 to leave it as i
  !        iunf = unit file for output; if <0, no residuals output
  !        [divrat=control for paralyzed target function from comdif.h90]
  ! Output elc corrected orbital elements at same time
  !        uncert includes covariance matrix and inverse
  !           warning: if only some elements are corrected, only the
  !              corresponding entries in uncert%g and uncert%c are nonzero
  !        csinor residuals norm
  !        delnor differential corrections norm
  !        obsw contains residuals, chi, selection flags, etc.
  !        succ logical success flag
  ! =============INTERFACE===== =========================================
  SUBROUTINE diff_cor(m,obs,obsw,el0,icor,iunf,elc,uncert,csinor,delnor,succ,nd,c_pr)
    ! ================input data==========================
    INTEGER m ! no. observations
    ! new data types
    TYPE(ast_obs),DIMENSION(m), INTENT(IN) :: obs
    TYPE(ast_wbsr),DIMENSION(m) :: obsw
    INTEGER iunf ! unit file to write output
    INTEGER, INTENT(IN) :: nd ! dimension of solve for parameters vector
    INTEGER, INTENT(IN) :: icor(nd) ! controls
    ! initial equinoctal elements
    TYPE(orbit_elem), INTENT(IN) :: el0
    ! ================output ==========================
    ! corrected elements
    TYPE(orbit_elem), INTENT(INOUT) :: elc ! to allow for el0=elc
    ! normal and covar. matrix
    TYPE(orb_uncert), INTENT(OUT) :: uncert
    ! success flag
    LOGICAL, INTENT(OUT) :: succ
    !  corr. and residuals norm
    DOUBLE PRECISION,INTENT(OUT) :: delnor,csinor
    ! a priori normal matrix; if aprc_flag is true, has to be there to be used
    DOUBLE PRECISION,INTENT(IN), OPTIONAL :: c_pr(nd,nd)
    ! =============END INTERFACE============================================
    ! input data sorted,  new data types
    TYPE(ast_obs),DIMENSION(nobx) :: obs_s
    TYPE(ast_wbsr),DIMENSION(nobx) :: obsw_s
    INTEGER iposs(nobx)
    ! residuals norm: previous and first (zeroth) iterations
    DOUBLE PRECISION csino1,csinor0
    ! differential correction, eccentricity,rescaling
    INTEGER nused,nsolv
    DOUBLE PRECISION ecc,mu,rescov
    ! difference from a priori (better with max dim?)
    DOUBLE PRECISION diff(nd)
    ! loop variables: it iteration number, itg=iter. with increase of q,
    ! ittodo= itmax, but forced to 1 for itmax=0 (no correction)
    INTEGER it,itg,itrej
    ! loop indexes: j=1,m; i,k=1,6, ibk for blocks, ii for nongrav
    INTEGER j,k,i,ibk, ii
    ! Flag for only computing covariance, not changing orbit.
    LOGICAL matonly
    ! unit for output
    INTEGER iun
    ! outlier rejection
    INTEGER nmod
    ! function for control of bizarre orbit, deserving to give up
    LOGICAL bizarre
    ! ====================================================================
    ! ================== INITIALIZATION ==================================
    ! ====================================================================
    ! Are commons loaded?
    IF(iicdif.ne. 36) STOP 'difcor: internal error(1)'
    IF(iicrej.ne. 36) STOP 'difcor: internal error(2)'
    ! definition of blocks if necessary
    CALL blockset(m,obs,obsw)
    ! sort of times and reordering of obs. and weights
    CALL sort_obs(el0%t,obs,obsw,m,iposs,obs_s,obsw_s)
    ! count number of solve-for variables
    nsolv=0
    DO  j=1,nd
       IF(icor(j).ne.0)nsolv=nsolv+1
    ENDDO
    ! check on consistency of a-priori
    IF(aprc_flag.and..not.PRESENT(c_pr))THEN
       WRITE(*,*)' diff_cor: inconsistent choice of a-priori, no c_pr and aprc_flag= ', aprc_flag
       STOP
    ENDIF
    ! Initialisation with starting value for elements
    elc=el0
    ! output unit
    iun=abs(iunf)
    IF(verb_dif.gt.19.and.itmax.gt.0)write(*,220) el0%coo,el0%coord
    IF(iunf.gt.0.and.verb_dif.gt.9.and.itmax.gt.0)THEN
       IF(nd.eq.6)THEN
          write(iun,220) el0%coo,el0%coord
220       format(' starting values ', a3/6f15.9)
       ELSE
          write(iun,221) el0%coo,el0%coord, dyn%dp(1:dyn%ndp)
221       format(' starting values ', a3/10f15.9)
       ENDIF
    ENDIF
    ! norm of corrections for "zeroth" iteration
    csino1=1.0d25
    ! control for covariance only pass
    matonly=.false.
    IF(itmax.eq.0)THEN
       itmax=1
       matonly=.true.
    ENDIF
    ! ====================================================================
    ! ================== BEGIN OUTLIER LOOP ==============================
    ! ====================================================================
    DO itrej=1,itmaxr
       ! counter for increasing/paralyzed iterations
       itg=0
       ! +++++++++++++++++++++ Begin Inner Loop +++++++++++++++++++++++++++++
       DO it=1,itmax
          ! ================== single iteration ==============================
          IF(aprc_flag)THEN
             diff(1:6)=elc%coord-el0%coord
             IF(nd.gt.6)THEN
                WRITE(*,*)' diff_cor: case with C a-priori and nongrav not ready, nd= ',nd
                STOP
             ENDIF
             CALL sin_cor(m,obs_s,obsw_s,elc,icor,iun,matonly,delnor,csinor,uncert,nd,C_PR=c_pr,DIFF=diff)
          ELSE
             CALL sin_cor(m,obs_s,obsw_s,elc,icor,iun,matonly,delnor,csinor,uncert,nd)
          ENDIF
          ! save first pass (sorted) data
          ! These will be passed out in case of hyperbolicity
          IF(it.eq.1 .and. itrej.eq.1)THEN
             CALL unsort_obs(iposs,m,obsw_s,obsw)
             csinor0=csinor
             ! cleanup the chi2
             obsw_s(1:m)%chi=0.d0
          ENDIF
          ! case of matrix inversion failure
          IF(.not.uncert%succ)THEN
             succ=.false.
             ! return saved rms
             csinor=csinor0
             ! cleanup the chi2
             obsw%chi=0.d0
             return
          ENDIF
          ! for computation of covariance without correction, there is no rejectio
          ! and no failure possible
          IF(matonly)THEN
             succ=.true.
             uncert%succ=.true.
             IF(verb_dif.ge.9)write(iun,*)'One pass only. No corrections.'
             !           IF(verb_dif.ge.19)write(*,*) 'One pass only. No corrections.'
             GOTO 70
          ENDIF
          ! control against hyperbolic and bizarre orbits
          IF(bizarre(elc,ecc))THEN
             IF(verb_dif.ge.19)write(*,*)' diff_cor: iter. ',it,' bizarre; e=',ecc
             IF(verb_dif.ge.9)write(iun,*)' diff_cor: iter. ',it,' bizarre; e=',ecc
             !           WRITE(iun,*)elc%coo, elc%coord
             succ=.false.
             !     return saved rms
             csinor=csinor0
             ! cleanup the chi2
             obsw%chi=0.d0
             return
          ENDIF
          ! norm of the correction: the zeros DO not matter!
          IF(verb_dif.ge.9)write(iun,200)it,csinor,delnor,elc%coord
          IF(verb_dif.ge.19)write(*,200)it,csinor,delnor,elc%coord
200       format(' difcor iteration ',i3/'  RMS residuals =',1p,d12.4,&
               &           '   norm corr =',d12.4/,'  new elem values'/0p,6f13.7/)
          IF(verb_dif.gt.19.and.nd.gt.6)THEN
             WRITE(*,*)'nongrav ',dyn%dp(ls(1:nd-6))
             DO ii=7,nd
                WRITE(*,*) 'sigma(',ii,')=',SQRT(ABS(uncert%g(ii,ii)))
             ENDDO
          ENDIF
          IF(verb_dif.ge.9)THEN
             WRITE(iun,*)'nongrav ',dyn%dp(ls(1:nd-6))
             DO ii=7,nd
                WRITE(iun,*) 'sigma(',ii,')=',SQRT(ABS(uncert%g(ii,ii)))
             ENDDO
          ENDIF
          ! ====== DECISION #1 ===================================================
          ! Check if we need another iteration
          ! Small corrections to elements?
          IF(delnor.lt.delrej)THEN
             succ=.true.
             uncert%succ=.true.
             GOTO 77
          ENDIF
          ! Target function increasing/paralyzed?
          IF(csinor.gt.csino1*1.1d0)THEN
             itg=itg+1
             IF(verb_dif.ge.9)write(iun,*)' target function increasing '
             IF(verb_dif.ge.19)write(*,*)' target function increasing '
             succ=.false.
             !           succ=.true.
          ELSEIF(csinor.gt.csino1*divrat)THEN
             itg=itg+1
             IF(verb_dif.ge.9)write(iun,*)' target function paralyzed '
             IF(verb_dif.ge.19)write(*,*)' target function paralyzed '
             succ=.true.
             uncert%succ=.true.
          ELSE
             itg=0
          ENDIF
          IF(itg.gt.itgmax)THEN
             ! go to outlier rejection if paralyzed
             IF(succ) GOTO 77
             !  otherwise (increasing) exit with failure
             GOTO 70
          ENDIF
          csino1=csinor
          ! +++++++++++++++++++++ End Inner Loop +++++++++++++++++++++++++++++
       ENDDO
       !     succ=.false.
       IF(verb_dif.ge.19)write(*,*)' too many iterations'
       IF(verb_dif.ge.9)write(iun,*)' too many iterations'
       IF(succ) GOTO 77
       !           Must get computed orbit and selection flags to agree:
       IF(autrej)THEN
          IF(aprc_flag)THEN
             CALL sin_cor(m,obs_s,obsw_s,elc,icor,iun,matonly,delnor,csinor,uncert,nd,C_PR=c_pr,DIFF=diff)
          ELSE
             CALL sin_cor(m,obs_s,obsw_s,elc,icor,iun,matonly,delnor,csinor,uncert,nd)
          ENDIF
       ENDIF
       ! case of matrix inversion failure
       IF(.not.uncert%succ)THEN
          succ=.false.
          ! return saved rms
          csinor=csinor0
          ! cleanup the chi2
          obsw%chi=0.d0
          return
       ENDIF
       GOTO 70
       ! ================== AUTOMATIC OUTLIER REJECTION =====================
77     CONTINUE
       IF(autrej)THEN
          IF(verb_rej.ge.19)write(*,*)'REJECTING...'
          ! call with -iun to avoid logging rejection stats
          CALL reject_obs(iunf,csinor,obs_s,obsw_s,m,uncert%g(1:nd,1:nd),icor,nmod,nd)
          ! update weight vector with zeros for rejected obs
          ! not necessary
          IF(verb_rej.ge.19)write(*,201)itrej,it,nmod
          IF(verb_rej.ge.9)write(iun,201)itrej,it,nmod
201       format('Outlier rejection pass ',i3,                     &
               &         ' after ',i3,' single differential correction passes.'/  &
               &              'There were ',i3,' changes.')
       ELSE
          obsw_s%chi=0.d0
          nmod=0
          ! norm of the correction: the zeros DO not matter!
          IF(verb_dif.ge.9)write(*,200)it,csinor,delnor,elc%coord
          IF(verb_dif.gt.7)write(iun,222)it,csinor,delnor
222       format(' done difcor iter. ',i3,'  RMS =',1p,d12.4,&
               &           ' last corr. =',d12.4)
       ENDIF
       ! go to next loop if no modifications or if not rejecting anyway
       IF(nmod.eq.0 .or. .not.autrej)THEN
          GOTO 75
       ENDIF
       ! ================== END OF OUTLIER LOOP ==============================
    ENDDO
    ! ====================================================================
    ! ================== BEGIN FINAL LOOP ================================
    ! ====================================================================
75  CONTINUE
    !      WRITE(*,*)delnor,delcr,succ

    IF(delnor.gt.delcr)THEN
       IF(verb_dif.ge.9)write(iun,*)'Final convergence loop after ', &
            &           itrej,' outlier rejection loops.'
       IF(verb_dif.ge.19)write(*,*)'Final convergence loop after ', &
            &           itrej,' outlier rejection loops.'
       ! counter for increasing/paralyzed iterations
       itg=0
       DO it=1,itmax
          ! ================== single iteration ==================================
          IF(aprc_flag)THEN
             CALL sin_cor(m,obs_s,obsw_s,elc,icor,iun,matonly,delnor,csinor,uncert,nd,C_PR=c_pr, DIFF=diff)
          ELSE
             CALL sin_cor(m,obs_s,obsw_s,elc,icor,iun,matonly,delnor,csinor,uncert,nd)
          ENDIF
          ! case of matrix inversion failure
          IF(.not.uncert%succ)THEN
             succ=.false.
             ! return saved rms
             csinor=csinor0
             ! cleanup the chi2
             obsw%chi=0.d0
             return
          ENDIF
          ! control against hyperbolic and bizarre orbits
          IF(bizarre(elc,ecc))THEN
             IF(verb_dif.ge.19)write(*,*)'difcor: iter. ',it,' bizarre; e=',ecc
             IF(verb_dif.gt.9)write(iun,*)'difcor: iter. ',it,' bizarre; e=',ecc
             succ=.false.
             !     return saved rms
             csinor=csinor0
             ! cleanup the chi2
             obsw_s%chi=0.d0
             ! the residuals of the first iteration will be passed up
             RETURN
          ENDIF
          ! norm of the correction: the zeros DO not matter!
          IF(verb_dif.ge.19)write(*,200)it,csinor,delnor,elc%coord
          IF(verb_dif.gt.9)write(iun,222)it,csinor,delnor
          IF(verb_dif.ge.19.and.nd.gt.6)THEN
             WRITE(*,*)'nongrav ',dyn%dp(ls(1:nd-6))
          ENDIF
          ! ====== DECISION # 2 ===============================================
          ! Check if we need another iteration
          ! Small corrections to elements?
          IF(delnor.lt.delcr)THEN
             IF(verb_dif.ge.9)THEN
                write(iun,*) 'Done. Corrections small after ',       &
                     &                 it,' passes.'
             ENDIF
             IF(verb_dif.ge.19)THEN
                write(*,*) 'Done. Corrections small after ',       &
                     &                 it,' passes.'
             ENDIF
             succ=.true.
             uncert%succ=.true.
             goto 70
          ENDIF
          !     Target function increasing/paralyzed?
          IF(csinor.gt.csino1*1.1d0)THEN
             itg=itg+1
             IF(verb_dif.ge.19)THEN
                write(*,*)' target function increasing '
             ENDIF
             IF(verb_dif.ge.9)THEN
                write(iun,*)' target function increasing '
             ENDIF
             succ=.false.
          ELSEIF(csinor.gt.csino1*divrat)THEN
             itg=itg+1
             IF(verb_dif.ge.9)THEN
                write(iun,*)' target function paralyzed '
             ENDIF
             IF(verb_dif.ge.19)THEN
                write(*,*)' target function paralyzed '
             ENDIF
             succ=.true.
             uncert%succ=.true.
          ELSE
             itg=0
          ENDIF
          IF(itg.gt.itgmax)THEN
             IF(verb_dif.ge.19)WRITE(*,*)                           &
                  &           'Done. Target funct. not decreasing after ',it,' passes.'
             IF(verb_dif.ge.9)write(iun,*)                               &
                  &           'Done. Target funct. not decreasing after ',it,' passes.'
             goto 70
          ENDIF
          csino1=csinor
          ! ====================== END FINAL LOOP =============================
       ENDDO
       succ=.false.
       IF(verb_dif.ge.5)write(*,*)'difcor: too many iterations ', it
       IF(verb_dif.gt.2)write(iun,*)'difcor: too many iterations ',it
    ENDIF
    ! ====================== CLEAN UP AND RETURN ========================
70  CONTINUE
    ! covariance (and normal matrix) rescaling
    nused=0
    DO i=1,m
       IF(obsw_s(i)%sel_coord.gt.0)THEN
          IF(obs_s(i)%type.eq.'O'.or.obs_s(i)%type.eq.'S')THEN
             nused=nused+2
          ELSEIF(obs_s(i)%type.eq.'R'.or.obs_s(i)%type.eq.'V')THEN
             nused=nused+1
          ENDIF
       ENDIF
    ENDDO
    IF(nsolv .gt. nused)THEN
       STOP 'difcor: internal error, nsolv>nused'
    ELSEIF(nsolv.eq.nused)THEN
       WRITE(*,*) ' difcor: nsolv=nused ', nsolv, nused, m, csinor
    ENDIF

    ! To be done better .... maybe done
    mu=rescov(nsolv,nused,csinor)
    DO  i=1,nd
       DO  j=1,nd
          uncert%g(i,j)=uncert%g(i,j)*mu**2
          uncert%c(i,j)=uncert%c(i,j)/mu**2
       ENDDO
    ENDDO
    ! reordering the residuals for output
    CALL unsort_obs(iposs,m,obsw_s,obsw)
    ! unsort  dmagn, dsunv,disv,phav,adotv,ddotv
    ! (data for magnitude computation)
    do i=1,m
       j=iposs(i)
       dsuna(j)=dsunas(i)
       disa(j)=disas(i)
       phaa(j)=phaas(i)
       dmagn(j)=dmagns(i)
       adotmv(j)=adots(i)
       ddotmv(j)=ddots(i)
       elomv(j)=elos(i)
    enddo
    !=====================
  END SUBROUTINE diff_cor
  !
  ! Copyright (C) 1998 by OrbFit Consortium
  ! Version: December 11, 1997 Steven Chesley
  ! Single iteration of differential correction of orbital elements
  !
  ! Input: m      observations number
  !        iobsrt sorted obs. type 1000's=astrometry, 2000's=radar
  !        tsort  sorted observation times
  !        als    sorted right ascensions
  !        des    sorted declination
  !        iocos  sorted observatory codes
  !        icor(nd) flags .ne.0 to correct this element, 0 to leave it as i
  !        nd     dimension of covariance matrix
  !        ws     sorted weight vector
  !        iun    unit for output
  !        matonly = .true. if only need covariance; do not update orbit
  !
  ! Output elc     corrected orbital elements at time elc%t
  !        delnor norm of the corrections
  !        csi    residuals
  !        csinor norm of residuals
  !        uncert includes covariance matrix and inverse
  !           warning: if only some elements are corrected, only the
  !              corresponding entries in gamma and gtwg are nonzero
  ! ============================================================
  SUBROUTINE sin_cor(m,obs_s,obsw_s,elc,icor,  &
       &              iun,matonly,delnor,csinor,uncert,nd,deqout,c_pr,diff)
    USE pred_obs
    USE close_app, ONLY: kill_propag
    ! INPUT :
    INTEGER, INTENT(IN) :: m ! number of observations
    TYPE(ast_obs),DIMENSION(m),INTENT(IN) :: obs_s ! observations
    LOGICAL, INTENT(IN) :: matonly ! control for application of corrections
    INTEGER, INTENT(IN) :: iun ! log file unit
    INTEGER, INTENT(IN) :: nd ! dimension of solve for parameters vector
    INTEGER, INTENT(IN) :: icor(nd) ! variables to be corrected
    ! a priori normal matrix; if logical flag present and true, to be used
    DOUBLE PRECISION,INTENT(IN), OPTIONAL :: c_pr(nd,nd), diff(nd)
    ! INPUT AND OUTPUT
    TYPE(ast_wbsr),DIMENSION(m) :: obsw_s  ! weights and residuals
    TYPE(orbit_elem), INTENT(INOUT) :: elc ! first guess, new oelements
    ! OUTPUT :
    ! norm of residuals, of correction
    DOUBLE PRECISION, INTENT(OUT) :: csinor,delnor
    TYPE(orb_uncert), INTENT(OUT) :: uncert ! normal and covariance matrix
    DOUBLE PRECISION, INTENT(OUT), OPTIONAL :: deqout(nd) ! correction applied
    ! END INTERFACE=========================
    TYPE(orbit_elem) elb
    INTEGER j,iocj,ider,i,k,ii
    INTEGER nused
    CHARACTER*1 tech ! obs technology, used by r_rdot
    DOUBLE PRECISION tauj
    DOUBLE PRECISION alj,dej,cond,deq(ndimx)
    INTEGER indp ! inversion failure flag
    DOUBLE PRECISION alobs ! duplicate R.A.,to preserve integrity of obs_s
    ! first derivatives of alpha, delta w.r. to elements
    DOUBLE PRECISION dade(ndimx),ddde(ndimx),drde(ndimx),dvde(ndimx)
    DOUBLE PRECISION gamma(ndimx,ndimx), gtwg(ndimx,ndimx)
    LOGICAL twobo  ! control of two body approximation (must be false)
    DOUBLE PRECISION snormd,snorm,appmag,pridif ! DOUBLE PRECISION functions
    DOUBLE PRECISION pos(3), vel(3) ! observer position, velocity (geocentric)
    ! proper motion, data to compute magnitude
    DOUBLE PRECISION :: adot,ddot,pha,dis,dsun,elo,gallat
    LOGICAL bizarre ! function for control of bizarre orbit
    DOUBLE PRECISION ecc ! eccentricity computed by bizarre
    INTEGER ir ! index for relaxation loop, to avoid bizarre corrections
    DOUBLE PRECISION peq(ndimx),sdir ! weak dir
    ! ===========================================
    ! inizialization of uncert
    CALL undefined_orb_uncert(nd,uncert)
    ! assign ider
    ider=1
    twobo=.false.
    !  twobo=.true.
    ! Compute observations and derivatives
    CALL set_restart(.true.)
    DO 61 j=1,m
       tauj=obs_s(j)%time_tdt
       iocj=obs_s(j)%obscod_i
       pos=obs_s(j)%obspos
       vel=obs_s(j)%obsvel
       tech=obs_s(j)%tech
       IF(obs_s(j)%type.eq.'O'.or.obs_s(j)%type.eq.'S')THEN
          CALL alph_del(elc,tauj,iocj,pos,vel,ider,twobo,nd,alj,dej,dade,ddde, &
               &           adot,ddot,pha,dis,dsun,elo,gallat)
          IF(kill_propag)THEN
             IF(verb_dif.gt.9) WRITE(*,*)' sin_cor: kill_propag'
             uncert%succ=.false.
             kill_propag=.false.
             RETURN
          ENDIF
          !  compute magnitude difference (apparent minus absolute)
          dmagns(j)=appmag(0.d0,elc%g_mag,dsun,dis,pha)
          phaas(j)=pha
          dsunas(j)=dsun
          disas(j)=dis
          adots(j)=adot
          ddots(j)=ddot
          elos(j)=elo
       ELSEIF(obs_s(j)%type.eq.'R'.or.obs_s(j)%type.eq.'V')THEN
          ! note: in this case vel=position of receiver, pos=position of transmitter
          CALL r_rdot(elc,tauj,iocj,tech,vel,pos,alj,dej,nd,drde,dvde,ider)
          ! in output alj=range, dej=range-rate
       ELSE
          WRITE(*,*)'sin_cor: obs%type= ',obs_s(j)%type, ' not known'
       ENDIF
       CALL set_restart(.false.)
       ! Compute residuals, form matrix g =-d(csi)/d(eq)
       IF(obs_s(j)%type.eq.'O'.or.obs_s(j)%type.eq.'S')THEN
          alobs=obs_s(j)%coord(1)-obsw_s(j)%bias_coord(1)
          obsw_s(j)%res_coord(1)=pridif(alobs,alj)
          obsw_s(j)%res_coord(2)=obs_s(j)%coord(2)-dej-obsw_s(j)%bias_coord(2)
          g(2*j-1,1:nd)=dade(1:nd)
          g(2*j,1:nd)=ddde(1:nd)
       ELSEIF(obs_s(j)%type.eq.'R')THEN
          obsw_s(j)%res_coord(1)=obs_s(j)%coord(1)-alj-obsw_s(j)%bias_coord(1)
          obsw_s(j)%res_coord(2)=0.d0
          g(2*j-1,1:nd)=drde(1:nd)
          g(2*j,1:nd)=0.d0
       ELSEIF(obs_s(j)%type.eq.'V')THEN
          obsw_s(j)%res_coord(1)=0.d0
          obsw_s(j)%res_coord(2)=obs_s(j)%coord(2)-dej-obsw_s(j)%bias_coord(2)
          g(2*j-1,1:nd)=0.d0
          g(2*j,1:nd)=dvde(1:nd)
       ELSE
          WRITE(*,*)'sin_cor: obs%type= ',obs_s(j)%type, ' not handled'
       ENDIF
       obsw_s(j)%resc_def=.true.
61  ENDDO
    CALL set_restart(.true.)
    ! =========================================================
    ! Compute solution of linear least squares
    IF(aprc_flag.and.PRESENT(c_pr).and.PRESENT(diff))THEN
       CALL min_sol(obs_s,obsw_s,m,g(:,1:nd),icor,iun,gtwg(1:nd,1:nd),deq(1:nd),    &
            &     gamma(1:nd,1:nd),csinor,indp,cond,nd,c_pr,diff)
    ELSE
       CALL min_sol(obs_s,obsw_s,m,g(:,1:nd),icor,iun,gtwg(1:nd,1:nd),deq(1:nd),     &
            &    gamma(1:nd,1:nd),csinor,indp,cond,nd)
    ENDIF
    IF(verb_dif.gt.10)THEN
       DO ii=7,nd
          WRITE(*,*)'gamma(',ii,',',ii,')=',gamma(ii,ii)
          WRITE(*,*) 'sigma=',SQRT(ABS(gamma(ii,ii)))
       ENDDO
    ENDIF
    IF(PRESENT(deqout))  deqout=deq
    IF(indp.ne.0) GOTO 9
    ! Update solution (only solve-for variables)
    IF(.not.matonly)THEN
       DO k=1,nd
          IF(k.le.6)THEN
             ! correct one of the orbital elements
             IF(icor(k).ne.0) elc%coord(k)=elc%coord(k)+deq(k)
          ELSEIF(k.gt.6.and.k.le.6+dyn%ndp)THEN
             ! if nongrav are determined, there must be 6 elements to be determined and the nongrav paramaters
             ! from the vector dyn are the last parameters to be determined
             dyn%dp(ls(k-6))=dyn%dp(ls(k-6))+deq(k)
          ENDIF
       ENDDO
       ! relaxation loop
       !     elb=elc
       ! Update solution
       !     elb%coord=elc%coord+deq
       !     IF(bizarre(elb,ecc))THEN
       !        IF(verb_dif.ge.10)WRITE(*,*)'sincor: short step to avoid ecc=',ecc
       !        CALL weak_dir(gamma,peq,sdir,iun,elc%coo,elc%coord)
       !        deq=deq-0.8d0*DOT_PRODUCT(peq,deq)*peq
       !        elb%coord=elc%coord+deq
       !     ENDIF
       !     elc=elb
    END IF
    ! Norm of the corrections (not right for case with some icor(j)=0)
    delnor=snorm(deq(1:nd),gtwg(1:nd,1:nd),nd,nd)
    ! output matrices
    uncert%g(1:nd,1:nd)=gamma(1:nd,1:nd)
    uncert%c(1:nd,1:nd)=gtwg(1:nd,1:nd)
    uncert%ndim=nd
    uncert%succ=.true.
    RETURN
9   IF(verb_dif.gt.9)THEN
       WRITE(*,*)' inversion failed, at row ', indp
       WRITE(iun,*)' inversion failed, at row ', indp
    ENDIF
    uncert%succ=.false.
  END SUBROUTINE sin_cor
  ! ==========================================
  ! BLOCKSET
  ! defines blocks as a function of the control delta t=dtblock
  !
  SUBROUTINE blockset(m,obs,obsw)
    IMPLICIT NONE
    ! ============INPUT=======================
    ! no. observations
    INTEGER,INTENT(IN) :: m
    ! new data types
    TYPE(ast_obs),DIMENSION(m), INTENT(IN) :: obs
    TYPE(ast_wbsr),DIMENSION(m),INTENT(IN) :: obsw
    ! ==== END INTERFACE ======================
    INTEGER ibk, nolarge ! index of block, pointer to large block
    DOUBLE PRECISION cov_small(nxinbl,nxinbl) !workspace normal blocks
    DOUBLE PRECISION cov_large(nxinbl_large,nxinbl_large) ! id. large blocks
    ! definition of blocks while observations are still sorted by time
    IF(error_model_priv.ne.' '.and.error_model_priv.ne.'cbm10'&
         & .and.error_model_priv.ne.'fcct14'.and.error_model_priv.ne.'vfcc17') THEN
       CALL blockdef(obs,m)
       large_blk(1:nblk)=0
       nolarge=0
       DO ibk=1,nblk
          IF(2*noblk(ibk).gt.nxinbl)THEN
             WRITE(iun_log,*) ' block ',ibk, ' large, no. obs ', noblk(ibk), ' station ',stablk(ibk)
             nolarge=nolarge+1
             IF(nolarge.gt.nbl_large)THEN
                STOP ' **** increase nbl_large ********'
             ENDIF
             large_blk(ibk)=nolarge
             ! full weight matrix for the block, large case
             CALL blocomp(m,obs,obsw,stablk(ibk),indblk(1,ibk),  &
                  &            nxinbl_large,noblk(ibk),wblk_large(1,1,nolarge),cov_large)
          ELSE
             ! full weight matrix for the block, normal case
             CALL blocomp(m,obs,obsw,stablk(ibk),indblk(1,ibk),  &
                  &            nxinbl,noblk(ibk),wblk(1,1,ibk),cov_small)
          ENDIF
       ENDDO
    ENDIF
  END SUBROUTINE blockset
  ! ==========================================
  ! BLOCKDEF
  ! defines blocks as a function of the control delta t=dtblock
  !
  SUBROUTINE blockdef(obs,m)
    ! ============INPUT=======================
    ! observations
    INTEGER,INTENT(IN) :: m ! no. observations
    TYPE(ast_obs),DIMENSION(m), INTENT(IN) :: obs
    ! ===========INTERNAL OUTPUT=========================
    ! block variables: no blocks, station and time limites for each block
    ! ========================================
    DOUBLE PRECISION tt
    ! loop indexes
    INTEGER j,k,jj
    ! dtblock is zero for radar
    DOUBLE PRECISION dtb
    fullblk=.false. ! all block are empty and available
    nblk=0
    DO j=1,m
       tt=obs(j)%time_tdt
       IF(obs(j)%type.eq.'O'.or.obs(j)%type.eq.'S')THEN
          dtb=dtblock
       ELSEIF(obs(j)%type.eq.'R'.or.obs(j)%type.eq.'V')THEN
          dtb=0.d0
       ENDIF
       DO k=1,nblk
          IF(obs(j)%obscod_i.eq.stablk(k))THEN
             IF(noblk(k)*2.lt.nxinbl_large)THEN
                IF(tblock(1,k).le.tt.and.tt.le.tblock(2,k))THEN
                   ! belonging to a known block, no changes to boundaries
                   ! store index, increase counter
                   jj=noblk(k)+1
                   indblk(jj,k)=j
                   noblk(k)=jj
                   ! this should happen only when there are two obs. at the same time
                   IF(verb_dif.gt.19)WRITE(*,*)'blockdef: duplicate obs? ',obs(j)%time_tdt,obs(j)%type,tblock(1:2,k)
                   !                 numerr=numerr+1
                ELSEIF(tblock(1,k).ge.tt.and.tblock(1,k)-dtb.le.tt)THEN
                   WRITE(*,*)'blockdef: wrong ',nblk,tt,tblock(1:2,k)
                   STOP
                   ! in a known block, extension on the left  NO
                   !                tblock(1,k)=tt
                   !                GOTO 1
                ELSEIF(tt.ge.tblock(2,k).and.tt.le.tblock(2,k)+dtb)THEN
                   ! in a known block, extension on the right
                   tblock(2,k)=tt
                   ! store index, increase counter
                   jj=noblk(k)+1
                   indblk(jj,k)=j
                   noblk(k)=jj
                   GOTO 1
                ENDIF
             ELSE
                IF(.not.fullblk(k)) THEN
                   ! block becomes full
                   ! create new block to avoid array bound error
                   WRITE(ierrou,*)' blockdef: block too large, station ',obs(j)%obscod_s
                   numerr=numerr+1
                   WRITE(ierrou,*)' blockdef: block too large, station ',obs(j)%obscod_s
                   IF(nblk+1.gt.nblx)STOP '**** increase nblx *********'
                   fullblk(k)=.true.
                   ! if the block is already known to be full, another one has been created already
                ENDIF
             ENDIF
          ENDIF
          ! not in this block
       ENDDO
       ! block not found; create another one
       IF(nblk+1.gt.nblx)STOP '**** increase nblx *********'
       nblk=nblk+1
       tblock(1,nblk)=tt
       tblock(2,nblk)=tt
       stablk(nblk)=obs(j)%obscod_i
       noblk(nblk)=1
       indblk(1,nblk)=j
1      CONTINUE
    ENDDO
    !  WRITE(iun_log,*)' number of blocks ',nblk
    DO k=1,nblk
       IF(noblk(k)*2.le.nxinbl_large)THEN
          !        WRITE(*,*) 'blockdef: station ',stablk(k),' block of ',nxinbl_large
          !        WRITE(iun_log,*) 'blockdef: station ',stablk(k),' block of ', &
          ! &  noblk(k)*2
       ELSE
          WRITE(iun_log,*) 'blockdef: station ',stablk(k),' block of ', &
               &  noblk(k)*2,' max ',nxinbl_large
          WRITE(ierrou,*) 'blockdef: station ',stablk(k),' block of ', &
               &  noblk(k)*2,' max ',nxinbl_large
          numerr=numerr+1
       ENDIF
    ENDDO
  END SUBROUTINE blockdef
  ! ========================================
  ! BLOCOMP
  !  computes the non-diagonal weight matrix for a given block
  ! ========================================
  SUBROUTINE blocomp(m,obs,obsw,station,indblock,        &
       &       nxxbl,noinblock,wblock,cov)
    USE obs_correl
    IMPLICIT NONE
    ! ============INPUT=======================
    ! observations
    INTEGER,INTENT(IN) :: m ! no. observations
    TYPE(ast_obs),DIMENSION(m), INTENT(IN) :: obs
    TYPE(ast_wbsr),DIMENSION(m),INTENT(IN) :: obsw
    ! ============INPUT==================
    ! list of indexes, number of times (actal, max), station code
    INTEGER, INTENT(IN) :: noinblock
    INTEGER,INTENT(IN) :: indblock(noinblock),nxxbl,station
    ! ============OUTPUT=================
    ! weight matrix for the block,
    DOUBLE PRECISION, INTENT(OUT) :: wblock(nxxbl,nxxbl)
    ! workspace for covariance matrix
    DOUBLE PRECISION, INTENT(OUT) :: cov(nxxbl,nxxbl)
    ! ============END INTERFACE==========
    INTEGER j,k,jj,kk,iunf ! loop indexes, units
    ! off-diagonal covariances, diagonal weights
    DOUBLE PRECISION cova,covd, wj(3),wk(3),det
    DOUBLE PRECISION cond ! conditioning number
    INTEGER iun ! unit for input
    INTEGER indp ! flag indicating non positive definite if >0
    ! downsizing of non diagonal part of observations covariance
    DOUBLE PRECISION, PARAMETER :: cdec_ratio=0.9d0 ! decreasing factor
    DOUBLE PRECISION cov_decr ! accumulated decrease
    ! ==========================================================
    ! check bounds
    IF(noinblock.gt.nxinbl_large)THEN
       WRITE(*,*)' blocomp: block too big for station ',station,noinblock,nxinbl
       STOP
    ENDIF
    ! decision: use correlations whenever they are known
    ! covariance matrix
    DO jj=1,noinblock
       j=indblock(jj)
       ! diagonal weight
       CALL fit_weight(obs(j),obsw(j),.false.,wj)
       DO kk=1,noinblock
          k=indblock(kk)
          ! diagonal weight
          CALL fit_weight(obs(k),obsw(k),.false.,wk)
          IF(jj.eq.kk)THEN
             det=wj(1)*wj(2)-wj(3)**2
             ! variance of RA
             cov(2*jj-1,2*jj-1)=wj(2)/det
             ! variance of DEC
             cov(2*jj,2*jj)=wj(1)/det
             ! correlation between RA and DEC
             cov(2*jj-1,2*jj)=-wj(3)/det
             cov(2*jj,2*jj-1)=cov(2*jj-1,2*jj)
          ELSE
             ! is anyway zero where there are no correlations
             cov(2*jj-1,2*kk)=0.d0
             cov(2*jj,2*kk-1)=0.d0
             ! add the nondiagonal part,
             CALL obscor(error_model_priv,obs(j),wj,obs(k),wk,cova,covd)
             ! covariance on RA
             cov(2*jj-1,2*kk-1)=cova
             cov(2*kk-1,2*jj-1)=cova
             ! covariance on DEC
             cov(2*jj,2*kk)=covd
             cov(2*kk,2*jj)=covd
          ENDIF
       ENDDO
    ENDDO
    ! now we need the weight matrix...
    cov_decr=1.d0
1   CONTINUE
    CALL invmat(wblock,nxxbl,2*noinblock,cov,cond,indp,ierrou)
    IF(indp.gt.0)THEN
       cov_decr=cov_decr*cdec_ratio
       WRITE(ierrou,*)' blocomp: decreasing covar. by ',cov_decr,noinblock
       IF(verb_dif.gt.19)WRITE(*,*)' blocomp: decreasing covar. by ',cov_decr,noinblock
       numerr=numerr+1
       DO jj=1,2*noinblock
          DO kk=1,2*noinblock
             IF(jj.ne.kk)cov(jj,kk)=cdec_ratio*cov(jj,kk)
          ENDDO
       ENDDO
       GOTO 1
    ENDIF
    !      WRITE(*,*)'blocomp: station,noinblock, cond ', station,noinblock,cond
    !      IF(cond.gt.1e5)WRITE(*,*)station, noinblock, cond
    RETURN
  END SUBROUTINE blocomp

  !
  ! Copyright 1998, 2002 Orbfit Consortium
  ! Modified Dec. 2 1998 by Steve Chesley to include sorting of selection
  ! Modified 2001 by A. Milani to use quick sort
  ! Modified 2002 by A. Milani to use new data types
  ! SRTOSS: remake of srtarc
  ! Ordinamento dei tempi di osservazione di un arco rispetto a t0
  ! INPUT:    t0       -  Tempo centrale dell'arco
  !           obs      -  Observations
  !           obsw     -  observations weights, residuals, etc.
  !           noss     -  Numero di osservazioni
  ! OUTPUT:
  !           obs_s  -  Obs. sorted
  !           obsw_s  -  Obs. weights sorted
  !         iposs2(k) -  Mappa: I_{noss} --> I_{noss}

  SUBROUTINE sort_obs(t0,obs,obsw,noss,iposs2,obs_s,obsw_s)
    implicit none
    ! INPUT obs number and new data types
    INTEGER noss
    TYPE(ast_obs),DIMENSION(noss) :: obs
    TYPE(ast_wbsr),DIMENSION(noss) :: obsw
    !       time of initial conditions
    double precision t0
    ! OUPUT: new data types
    TYPE(ast_obs),DIMENSION(noss) :: obs_s
    TYPE(ast_wbsr),DIMENSION(noss) :: obsw_s
    !      sorting indexes
    integer  iposs2(noss)
    ! END INTERFACE
    ! indexes
    integer ii,kk,ibk,j,jj,jj2
    ! input times
    double precision toss(nobx)
    !      sorting indexes for increasing time, for inverse map
    integer iposs(nobx),iposs3(nobx)
    !****************
    !   static memory not required
    !****************
    toss(1:noss)=obs(1:noss)%time_tdt
    ! Sort: generating the indexes, toss is not changed
    CALL heapsort(toss,noss,iposs)
    ! Reordering of right ascensions declinations and weigths
    ! observatory codes, observation types and Selection Flags
    ! first the ones with time of observation > t0
    kk=0
    DO  ii=1,noss
       IF(toss(iposs(ii)).ge.t0)THEN
          kk=kk+1
          obs_s(kk)=obs(iposs(ii))
          obsw_s(kk)=obsw(iposs(ii))
          iposs2(kk)=iposs(ii)
          iposs3(iposs(ii))=kk
       ENDIF
    ENDDO
    ! now those with time of observation < t0
    DO  ii=noss,1,-1
       IF(toss(iposs(ii)).lt.t0)THEN
          kk=kk+1
          obs_s(kk)=obs(iposs(ii))
          obsw_s(kk)=obsw(iposs(ii))
          iposs2(kk)=iposs(ii)
          iposs3(iposs(ii))=kk
       ENDIF
    ENDDO
    IF(error_model_priv.ne.' '.and.error_model_priv.ne.'cbm10'&
         & .and.error_model_priv.ne.'fcct14'.and.error_model_priv.ne.'vfcc17') THEN
       ! now reordering of the indirect address tables for blocks
       DO ibk=1,nblk
          DO j=1,noblk(ibk)
             jj=indblk(j,ibk)
             jj2=iposs3(jj)
             indblk_s(j,ibk)=jj2
          ENDDO
       ENDDO
    ENDIF
    RETURN
  END SUBROUTINE sort_obs
  !==================================================
  !================     UNSORT     ==================
  !==================================================
  ! UNSORT: reordering of residuals in original order
  !
  ! Input:   iposs
  !          noss
  !          nos2
  !          csi
  !          tsort
  !          toss
  !          sels
  !
  ! Output:  csir
  !          sel
  !
  !==================================================
  SUBROUTINE unsort_obs(iposs,noss,obsw_s,obsw)
    implicit none
    ! INPUT: number of observations, sorting indexes
    INTEGER, INTENT(IN) ::  noss,iposs(noss)
    !    new data types
    TYPE(ast_wbsr),DIMENSION(noss), INTENT(IN) :: obsw_s
    ! OUPUT: new data types (but keeps some fields from input)
    TYPE(ast_wbsr),DIMENSION(noss), INTENT(INOUT) :: obsw
    ! END INTERFACE
    !      double precision tsort(noss),toss(noss),tmp
    integer sel(noss),sels(noss)
    !     loop indexes
    integer n,k
    !      toss(1:noss)=obs(1:noss)%time_tdt
    !      tsort(1:noss)=obs_s(1:noss)%time_tdt
    ! unsort checking that the times are back to the old order
    do  n=1,noss
       k=iposs(n)
       !         tmp=toss(k)-tsort(n)
       !         if(tmp.ne.0.d0)then
       !            write(*,*)'sorting problem ',tmp,n,tsort(n),k,toss(k)
       !         endif
       obsw(k)=obsw_s(n)
    enddo
    return
  END SUBROUTINE unsort_obs
  !
  ! ===========================================================
  !
  !                     M I N _ S O L
  ! version 1.3, A. Milani, 28/6/1997
  ! version 3.0, A. Milani November 2002, with correlated observations
  ! version 4.0, F. Spoto, April 2015: added non-gravitational parameters
  !
  ! least squares solver, including normal matrix normalisation
  !      INPUT :   obs_s, obsw_s   =  vector of observations and residuals
  !                m    =  observations number
  !                g     =  matrix of first derivatives of residuals
  !                          multiplied by -1
  !                icor(k)=  .ne.0 to correct this element, 0 to let it as
  !                iunf   = unit for output file
  !      OUTPUT:   gtwg  =  normal matrix GtWG
  !                dx0   =  corrections vector
  !                gamma =  inverse matrix (if all icor=1, else pseudo-inverse)
  !                csinor = norm of the residuals
  !                cond  =  conditioning number of gamma
  !                indp  = control; if not 0 computation not done
  !           warning: if only some elements are corrected, only the
  !              corresponding entries in gamma and gtwg are nonzero in ou
  ! ===========================================================
  SUBROUTINE min_sol(obs_s,obsw_s,m,g,icor,iunf,gtwg,dx0,gamma,csinor,indp,cond,nd,c_pr,diff)
    ! ===========================================================
    ! INPUT
    INTEGER, INTENT(IN) :: m ! number of  observations
    TYPE(ast_obs),DIMENSION(m),INTENT(IN) :: obs_s ! OBSERVATIONS
    TYPE(ast_wbsr),DIMENSION(m), INTENT(IN) :: obsw_s
    INTEGER, INTENT(IN) :: nd ! dimension of elements vector
    DOUBLE PRECISION, INTENT(IN) :: g(nob2x,nd) ! partial derivatives
    INTEGER, INTENT(IN) :: iunf ! unit for output file
    INTEGER, INTENT(IN) :: icor(nd) ! controls for corrections
    ! OUTPUT
    ! normal matrix, differential correction, covariance matrix, norm of residuals
    DOUBLE PRECISION, INTENT(OUT) :: gtwg(nd,nd),dx0(nd),gamma(nd,nd),csinor
    INTEGER, INTENT(OUT):: indp ! flag indicating non positive definite if >0
    DOUBLE PRECISION, INTENT(OUT) :: cond ! condition number
    ! a priori normal matrix; if logical flag is true nd present, to be used
    DOUBLE PRECISION,INTENT(IN), OPTIONAL :: c_pr(nd,nd),diff(nd)
    ! END INTERFACE ================================
    DOUBLE PRECISION w(3) ! weights for the current observation
    DOUBLE PRECISION cr(ndimx,ndimx),gr(ndimx,ndimx) ! reduced normal, covariance matrix
    DOUBLE PRECISION gtwcsi(ndimx) ! right hand side of the normal equation
    INTEGER ndc ! no. of solve for variables
    DOUBLE PRECISION wo(nxinbl_large,nxinbl_large) ! array for weights of current block
    INTEGER noinbl ! actual number of obs. in this block
    ! loop indexes: j,k=1,nd; i=1,no; jj,kk=1,ndc ,; ibk=1,nblk;
    !  jjo,kko=1,m, jji,kki=1,2*m
    INTEGER j,k,i,jj,kk,ibk,jjo,kko,jji,kki
    ! temporary variables and workspace
    INTEGER iun, nused
    ! for qrinv
    INTEGER izer
    DOUBLE PRECISION aval(ndimx),det
    ! ===========================================================
    ! control on dimensioning
    if (2*m.gt.nob2x) THEN
       write(*,*) '2*m > nob2x in min_sol',m,nob2x
       stop 'min_sol: 2*m > nob2x'
    ENDIF
    ! check availability of error_model
    IF(.not.errmod_given)THEN
       WRITE(*,*)' min_sol: missing error_model'
       STOP
    ENDIF
    ! compute number of scalar data
    nused=0
    DO  i=1,m
       IF(obsw_s(i)%sel_coord.gt.0)THEN
          IF(obs_s(i)%type.eq.'O'.or.obs_s(i)%type.eq.'S')THEN
             nused=nused+2
          ELSEIF(obs_s(i)%type.eq.'R'.or.obs_s(i)%type.eq.'V')THEN
             nused=nused+1
          ENDIF
       ENDIF
    ENDDO
    IF(nused.eq.0)THEN
       indp=0
       WRITE(iun_log,*)' no observations?????', nused,m
       RETURN
    ENDIF
    ! initialization of gtwg, gtwcsi, including a-priori
    gtwcsi=0.d0
    csinor=0.d0
    gtwg=0.d0
    ! =========================================
    ! A-priori for non-gravitational parameters
    ! =========================================
    IF(dyn%apriori(1).AND.nd.GT.6)THEN
       gtwg(7,7)=1.d0/dyn%apriori_std(1)**2
       gtwcsi(7)=(-dyn%dp(1)+dyn%apriori_nominal(1))/dyn%apriori_std(1)**2
       csinor=csinor+(-dyn%dp(1)+dyn%apriori_nominal(1))**2/dyn%apriori_std(1)**2
    ENDIF
    IF(dyn%apriori(2).AND.nd.GT.6)THEN
       ! Yarkovsky a-priori
       IF(dyn%nmod.EQ.1.AND.ls(1).EQ.2)THEN
          ! nd=7, parameter of Yarkovsky is number 7
          IF(nd.ne.7) STOP '**** min_sol: error nd.ne.7 ****'
          gtwg(7,7)=1.d0/dyn%apriori_std(2)**2
          gtwcsi(7)=(-dyn%dp(2)+dyn%apriori_nominal(2))/dyn%apriori_std(2)**2
       ELSEIF(dyn%nmod.EQ.1.AND.ls(1).EQ.1.AND.ls(2).EQ.2.)THEN
          ! nd=8, parameter of Yarkovsky is number 8
          IF(nd.ne.8) STOP '**** min_sol: error nd.ne.8 ****'
          gtwg(8,8)=1/dyn%apriori_std(2)**2
          gtwcsi(8)=(-dyn%dp(2)+dyn%apriori_nominal(2))/dyn%apriori_std(2)**2
       ENDIF
       IF(dyn%nmod.EQ.2)THEN
          gtwg(8,8)=1/dyn%apriori_std(2)**2
          gtwcsi(8)=(-dyn%dp(2)+dyn%apriori_nominal(2))/dyn%apriori_std(2)**2
       ENDIF
       csinor=csinor+(-dyn%dp(2)+dyn%apriori_nominal(2))**2/dyn%apriori_std(2)**2
    ENDIF
    IF(dyn%apriori(3).AND.nd.GT.6)THEN
       gtwg(9,9)=1.d0/dyn%apriori_std(3)**2
       gtwcsi(9)=(-dyn%dp(3)+dyn%apriori_nominal(3))/dyn%apriori_std(3)**2
       csinor=csinor+(-dyn%dp(3)+dyn%apriori_nominal(3))**2/dyn%apriori_std(3)**2
    ENDIF
    IF(dyn%apriori(4).AND.nd.GT.6)THEN
       gtwg(10,10)=1.d0/dyn%apriori_std(4)**2
       gtwcsi(10)=(-dyn%dp(4)+dyn%apriori_nominal(4))/dyn%apriori_std(4)**2
       csinor=csinor+(-dyn%dp(4)+dyn%apriori_nominal(4))**2/dyn%apriori_std(4)**2
    ENDIF

    ! all a priori matrix (from previous orbit)
    IF(aprc_flag)THEN
       IF(PRESENT(c_pr))THEN
          gtwg(1:nd,1:nd)=c_pr(1:nd,1:nd)
          gtwcsi(1:nd)=MATMUL(c_pr(1:nd,1:nd),diff(1:nd))
       ELSE
          WRITE(*,*)' min_sol: inconsistent choice on apriori, c_pr not present, aprc_flag= ', aprc_flag
       ENDIF
    ENDIF
    ! ===========================================================
    IF(error_model_priv.eq.' '.or.error_model_priv.eq.'cbm10'.or.error_model_priv.eq.'fcct14' &
         & .or.error_model_priv.eq.'vfcc17')THEN
       ! normal matrix GtWG of the pseudo-Newton method
       DO 1 j=1,nd
          DO 2 k=1,nd
             !           gtwg(j,k)=0.d0
             IF(icor(k).ne.0.and.icor(j).ne.0)THEN
                DO i=1,m
                   ! diagonal weight
                   CALL fit_weight(obs_s(i),obsw_s(i),.true.,w)
                   ! g^twg is the normal matrix
                   gtwg(j,k)=gtwg(j,k)+g(2*i-1,j)*w(1)*g(2*i-1,k)+   &
                        &              g(2*i,j)*w(2)*g(2*i,k)    +   &
                        &              w(3)*(g(2*i,j)*g(2*i-1,k)+g(2*i-1,j)*g(2*i,k))
                ENDDO
             ENDIF
2         ENDDO
1      ENDDO
       ! ===========================================================
       ! Computation of vectors GtWcsi; also norm of residuals
       !     gtwcsi=0.d0
       !     csinor=0.d0
       DO  i=1,m
          CALL fit_weight(obs_s(i),obsw_s(i),.true.,w)
          gtwcsi(1:nd)=gtwcsi(1:nd)+(g(2*i-1,1:nd)*w(1)+g(2*i,1:nd)*w(3))*obsw_s(i)%res_coord(1)   &
               &         +(g(2*i-1,1:nd)*w(3)+g(2*i,1:nd)*w(2))*obsw_s(i)%res_coord(2)
          csinor=csinor+ w(1)*obsw_s(i)%res_coord(1)**2 +    &
               &         w(2)*obsw_s(i)%res_coord(2)**2 +    &
               &         2.d0*w(3)*obsw_s(i)%res_coord(1)*obsw_s(i)%res_coord(2)
       ENDDO
       ! ===========================================================
       ! with off-diagonal weights (if error_model is given)
    ELSEIF(error_model_priv.ne.' '.and.error_model_priv.ne.'cbm10'&
         & .and.error_model_priv.ne.'fcct14'.and.error_model_priv.ne.'vfcc17')THEN
       ! start from scratch
       !     gtwcsi=0.d0
       !     csinor=0.d0
       !     gtwg=0.d0
       ! loop on blocks
       DO 3 ibk=1,nblk
          noinbl=noblk(ibk)
          ! copy current weights of the block
          IF(large_blk(ibk).gt.0)THEN
             wo(1:2*noinbl,1:2*noinbl)=wblk_large(1:2*noinbl,1:2*noinbl,large_blk(ibk))
          ELSE
             wo(1:2*noinbl,1:2*noinbl)=wblk(1:2*noinbl,1:2*noinbl,ibk)
          ENDIF
          ! full weight matrix for the block wblk is available
          DO 33 jj=1,noinbl
             jjo=indblk_s(jj,ibk)
             jji=jjo*2-1
             IF(obsw_s(jjo)%sel_coord.gt.0)THEN
                DO 34 kk=1,noinbl
                   kko=indblk_s(kk,ibk)
                   kki=kko*2-1
                   IF(obsw_s(kko)%sel_coord.gt.0)THEN
                      ! Q_jj target function of the block
                      csinor=csinor+obsw_s(jjo)%res_coord(1)*wo(2*jj-1,2*kk-1)* &
                           &      obsw_s(kko)%res_coord(1)+obsw_s(jjo)%res_coord(2)*wo(2*jj,2*kk)*  &
                           &                          obsw_s(kko)%res_coord(2)
                      ! alpha-delta correlations are not considered in this version
                      !                    +obsw_s(jjo)%res_coord(1)*wo(2*jj-1,2*kk)*       &
                      !     &                          obsw_s(kko)%res_coord(2)               &
                      !                    +obsw_s(jjo)%res_coord(2)*wo(2*jj,2*kk-1)*       &
                      !     &                          obsw_s(kko)%res_coord(1)               &
                      ! contribution of this couple of observations to the normal matrix
                      DO 4 j=1,nd
                         ! contribution of this couple of observations to the right hand side of normal equations
                         gtwcsi(j)=gtwcsi(j)                  &
                              &        +g(jji,j)*wo(2*jj-1,2*kk-1)*obsw_s(kko)%res_coord(1)     &
                              &        +g(jji+1,j)*wo(2*jj,2*kk)*obsw_s(kko)%res_coord(2)
                         ! alpha-delta correlations are not considered in this version
                         !    &         +g(jji,j)*wo(2*jj-1,2*kk)*obsw_s(kko)%res_coord(2)  &
                         !    &         +g(jji+1,j)*wo(2*jj,2*kk-1)*obsw_s(kko)%res_coord(1)  &
                         DO 5 k=1,nd
                            IF(icor(k).ne.0.and.icor(j).ne.0)THEN
                               gtwg(j,k)=gtwg(j,k)+g(jji,j)*wo(2*jj-1,2*kk-1)*g(kki,k)&
                                    &                 +g(jji+1,j)*wo(2*jj,2*kk)* g(kki+1,k)
                               ! alpha-delta correlations are not considered in this version
                               !    &                 +g(jji+1,j)*wo(2*jj,2*kk-1)* g(kki,k)         &
                               !    &              +g(jji,j)*wo(2*jj-1,2*kk)* g(kki+1,k)         &
                            ENDIF
5                        ENDDO
4                     ENDDO
                   ENDIF
34              ENDDO
             ENDIF
33        ENDDO
          ! end update due to the block
3      ENDDO
    ENDIF
    ! norm of the residuals
    csinor=sqrt(csinor/nused)
    ! ===========================================================
    ! number of parameters to be solved
    ndc=0
    DO  j=1,nd
       IF(icor(j).gt.0) ndc=ndc+1
    ENDDO
    IF(ndc.eq.nd)THEN
       ! simplyfied version if all parameters to be determined
       CALL invmat(gamma,nd,nd,gtwg,cond,indp,iunf)
       ! qr inversion
       IF(indp.ne.0)THEN
          CALL qr_inv(gtwg,gamma,nd,izer,det,aval(1:nd))
          IF(izer.eq.0)THEN
             indp=0
             cond=MAXVAL(aval(1:nd))/MINVAL(aval(1:nd))
          ELSE
             WRITE(iun_log,*)' qrinv failed, izer=',izer
             RETURN
          ENDIF
       ENDIF
       goto 78
    ENDIF
    ! ===========================================================
    ! squeeze normal matrix into a smaller one
    jj=0
    DO 7 j=1,nd
       IF(icor(j).ne.0)THEN
          jj=jj+1
          kk=0
          DO 8 k=1,nd
             IF(icor(k).ne.0)THEN
                kk=kk+1
                cr(jj,kk)=gtwg(j,k)
             ENDIF
8         ENDDO
       ENDIF
7   ENDDO
    IF(jj.eq.kk.and.jj.ne.0)THEN
       ndc=jj
    ELSE
       write(*,*)' min_sol, this should not happen ',jj,kk
       stop 'min_sol: jj,kk '
    ENDIF
    ! ===========================================================
    ! Cholewski inversion
    CALL invmat(gr,ndimx,ndc,cr,cond,indp,iunf)
    ! qr inversion
    IF(indp.ne.0)THEN
       CALL qr_inv(cr(1:ndc,1:ndc),gr(1:ndc,1:ndc),ndc,izer,det,aval(1:ndc))
       IF(izer.eq.0)THEN
          indp=0
          cond=MAXVAL(aval(1:ndc))/MINVAL(aval(1:ndc))
       ELSE
          WRITE(iun_log,*)' qrinv failed, izer=',izer
          RETURN
       ENDIF
    ENDIF
    ! ===========================================================
    !  covariance matrix
    !  warning: rows and columns not involved in the correction
    !  are set to zero (should be infinite!)
    kk=0
    DO 14 k=1,nd
       IF(icor(k).ne.0)THEN
          kk=kk+1
          jj=0
          DO 16 j=1,nd
             IF(icor(j).ne.0)THEN
                jj=jj+1
                gamma(k,j)=gr(kk,jj)
             ELSE
                gamma(k,j)=0.d0
             ENDIF
16        ENDDO
       ELSE
          DO  j=1,nd
             gamma(k,j)=0.d0
          ENDDO
       ENDIF
14  ENDDO
    ! =======================
78  continue
    IF(cond.gt.1.d12)THEN
       iun=abs(iunf)
       WRITE(iun,100)cond
100    FORMAT('Conditioning number: ',1p,d12.5)
       IF(verb_dif.ge.9)WRITE(*,100)cond
    ENDIF
    ! compute correction
    !   warning: correction is zero automatically for non solved-for variable
    dx0=matmul(gamma,gtwcsi(1:nd))
  END SUBROUTINE min_sol

  ! ================================================
  !  INVMAT  ndim x ndim  matrix inversion
  ! in this version it is assumed that the input matrix a is symmetric,
  ! definite positive, and in the calling program has
  ! dimension nx x ndim; the output c is of the same shape
  ! and the output matrix is symmetric, definite positive
  ! If this is not the case, a warning is issued on the standard ouput
  ! ===========INTERFACE=============================
  SUBROUTINE invmat(c,nx,ndim,a,cond,indp,iunf)
    implicit none
    INTEGER,INTENT(IN) :: nx,ndim,iunf
    DOUBLE PRECISION,INTENT(IN)::  a(nx,ndim)
    ! OUPUT: inverted matrix, conditioning number
    DOUBLE PRECISION,INTENT(OUT):: c(nx,ndim),cond
    INTEGER,INTENT(OUT) :: indp ! flag indicating failure of pos. def.
    ! ========END INTERFACE=============================
    INTEGER nxx
    ! warning: here it is assumed that ndim never exceeds 6;
    ! April 2003: changed for bigger matrices
    parameter(nxx=2*nxinbl_large)
    DOUBLE PRECISION an(nxx),v(nxx)
    INTEGER i,j,iun
    DOUBLE PRECISION err,omax,omin,da,eps
    logical sym
    ! ==================================================
    ! check that the matrix is indeed symmetric
    sym=.true.
    iun=abs(iunf)
    eps=1.d2*epsilon(1.d0)
    DO 1 i=1,ndim
       DO 2 j=1,i-1
          da=abs((a(i,j)-a(j,i))/(a(i,j)+a(j,i)))
          IF(da.gt.1000*eps)THEN
             IF(verb_dif.ge.9) write(iun,*)'invmat: ',i,j,a(i,j),a(j,i),da,100*eps
             IF(verb_dif.ge.19)                                          &
                  &            write(*,*)'invmat: ',i,j,a(i,j),a(j,i),da,100*eps
             sym=.false.
          ENDIF
2      ENDDO
1   END DO
    IF(.not.sym)THEN
       IF(verb_dif.ge.9)write(iun,*)'invmat: input matrix not symmetric'
       IF(verb_dif.ge.19)write(*,*)'invmat: input matrix not symmetric'
    ENDIF
    ! ==========================================================
    ! Cholewski
    ! ==========================================================
    ! normalisation of columns of matrix to be inverted
    DO  i=1,ndim
       an(i)=sqrt(abs(a(i,i)))
    ENDDO
    DO i=1,ndim
       DO  j=1,ndim
          c(i,j)=a(i,j)/(an(i)*an(j))
       ENDDO
    ENDDO
    DO i=ndim+1,nx
       DO  j=1,ndim
          c(i,j)=0.d0
       ENDDO
    ENDDO
    ! first Cholewsky factorisation
    err=eps
    CALL tchol(c,nx,ndim,indp,err)
    IF(indp.eq.0)THEN
       ! Control of conditioning number of the inverted matrix
       omax=c(1,1)
       omin=c(1,1)
       DO i=2,ndim
          if (c(i,i).gt.omax) THEN
             omax=c(i,i)
          ENDIF
          if (c(i,i).lt.omin) THEN
             omin=c(i,i)
          ENDIF
       ENDDO
       cond=(omax/omin)**2
       CALL inver(c,v,nx,ndim)
       ! unnormalize the matrix by norm of columns
       DO  i=1,ndim
          DO  j=1,ndim
             c(i,j)=c(i,j)/(an(i)*an(j))
          ENDDO
       ENDDO
    ELSE
       IF(verb_dif.gt.9)THEN
          write(iun,*)' matrix is not positive definite'
          write(iun,*)' pivot number ',indp,' is ',c(indp,indp)
       ELSEIF(verb_dif.gt.19)THEN
          write(*,*)' matrix is not positive definite'
          write(*,*)' pivot number ',indp,' is ',c(indp,indp)
       ENDIF
    ENDIF
  END SUBROUTINE invmat

  ! computation of rms of residuals
  DOUBLE PRECISION FUNCTION rms_compute(obs,obsw,n)
    IMPLICIT NONE
    ! INPUT observations: new data types
    INTEGER, INTENT(IN) :: n  !number of observations
    TYPE(ast_obs),DIMENSION(n),INTENT(IN) :: obs
    TYPE(ast_wbsr),DIMENSION(n),INTENT(IN) :: obsw
    ! END INTERFACE
    DOUBLE PRECISION csinor, w(3)
    INTEGER nused,i
    ! Computation norm of residuals
    csinor=0.d0
    nused=0
    DO  i=1,n
       CALL fit_weight(obs(i),obsw(i),.true.,w)
       IF(obsw(i)%sel_coord.gt.0)THEN
          csinor=csinor+ w(1)*obsw(i)%res_coord(1)**2 + &
               &         w(2)*obsw(i)%res_coord(2)**2 + &
               &         2.d0*w(3)*obsw(i)%res_coord(1)*obsw(i)%res_coord(2)
          ! note that later the correlations corresponding to the
          ! error model error_model_priv is NOT taken into account
          IF(obs(i)%type.eq.'O'.or.obs(i)%type.eq.'S')THEN
             nused=nused+2
          ELSEIF(obs(i)%type.eq.'R'.or.obs(i)%type.eq.'V')THEN
             nused=nused+1
          ENDIF
       ENDIF
    ENDDO
    ! missing: effect of correlations!
    csinor=sqrt(csinor/nused)
    rms_compute=csinor
  END FUNCTION rms_compute

  !
  ! Copyright (C) 1998 by Mario Carpino (carpino@brera.mi.astro.it)
  ! Version: January 15, 1998; corrected AM/ZK March 16, 1998
  ! Modified 2 Dec. 1998 by Steven Chesley
  ! deeply modified by A. Milani for Fortran 90 version, November 2002
  ! corrected for spherical metric,
  ! assuming residuals are delta alpha, delta delta
  ! ---------------------------------------------------------------------
  !
  !  *****************************************************************
  !  *                                                               *
  !  *                         F I T _ W E I G T                     *
  !  *                                                               *
  !  *      Computation of observation weights for orbital fit       *
  !  *            for one observation, from the a-priori RMS         *
  !  *                                                               *
  !  *****************************************************************
  !
  ! INPUT:  from OBS1, OBSW1
  !           RMSA      -  A-priori RMS of right ascension (rad)
  !           RMSD      -  A-priori RMS of declination (rad)
  !           DEL       -  Declination
  !           SEL       -  Selection index
  !           TYPE      -  Observation type: R,V are 1-dim
  !         RES0       - set weight to zero for rejected
  ! OUTPUT:   W         -  Weights
  !
  SUBROUTINE fit_weight(obs1,obsw1,res0,w)
    IMPLICIT NONE
    ! INPUT 1 observation, new data types
    TYPE(ast_obs),INTENT(IN) :: obs1
    TYPE(ast_wbsr), INTENT(IN) :: obsw1
    LOGICAl,INTENT(IN):: res0
    ! OUTPUT 1 couple of weights
    DOUBLE PRECISION, INTENT(OUT) :: w(3)
    ! END INTERFACE
    INTEGER sel
    DOUBLE PRECISION rmsa,rmsd,del,corr,det,rmsra2,rmsd2,cov
    CHARACTER*1 type
    !
    rmsa=obsw1%rms_coord(1)
    rmsd=obsw1%rms_coord(2)
    rmsra2=rmsa**2
    rmsd2=rmsd**2
    sel=obsw1%sel_coord
    del=obs1%coord(2)
    type=obs1%type
    !
    IF(sel.EQ.0.and.res0) THEN
       w(1:3)=0.d0
    ELSEIF(type.eq.'O'.or.type.eq.'S')THEN
       corr=obsw1%correl
       cov=corr*rmsa*rmsd
       IF(corr.GT.0.d0)THEN
          det=rmsra2*rmsd2-cov**2
          w(1)=rmsd2/det
          w(2)=rmsra2/det
          w(3)=-cov/det
       ELSE
          w(1)=1.d0/rmsra2
          w(2)=1.d0/rmsd2
          w(3)=0.d0
       END IF
    ELSEIF(type.eq.'R')THEN
       ! for radar obs, residuals, weights and partials are all zero for the other component
       w(1)=1.d0/rmsra2
       w(2)=1.d0
       w(3)=0.d0
    ELSEIF(type.eq.'V')THEN
       w(1)=1.d0
       w(2)=1.d0/rmsd2
       w(3)=0.d0
    END IF
  END SUBROUTINE fit_weight
  !
  ! Copyright (C) 1998 by Mario Carpino (carpino@brera.mi.astro.it)
  !                    and Steve Chesley (chesley@@dm.unipi.it)
  ! Version: December 15, 1998 Steven Chesley
  ! hacked for radar: A. Milani, January 16, 1998
  ! fortran 90 version 3.0 A. Milani November 2002
  ! ---------------------------------------------------------------------
  !
  !  *****************************************************************
  !  *                                                               *
  !  *                      R E J E C T _ O B S                      *
  !  *                                                               *
  !  *            Rejection and recovery of outliers                 *
  !  *                                                               *
  !  *****************************************************************
  !
  ! INPUT:    UNILOG    -  File unit for logging rejections (<0 for no log
  !                     (overridden by verb_rej
  !           CSINOR    -  Norm of residuals
  !           SEL       -  Selection flag (0 = not selected)
  !           CSI       -  O-C residuals
  !           W         -  Observation weights
  !           NOBS      -  Number of observations
  !           GAMMA     -  Covariance matrix of orbital elements
  !           ICOR      -  List of solve-for parameters
  !
  ! OUTPUT:   SEL       -  Updated selection flags
  !           X2        -  CHI**2 residual
  !           NMOD      -  Number of modifications to the outlier list
  !                        (additions or deletions)
  !
  ! WARNING: in some cases the weights are zero, because a scalar observat
  !          has been squeezed in a two-dimensional observation
  !          This is the case for obs. types 3,4 (radar)
  SUBROUTINE reject_obs(unilog,csinor,obs_s,obsw_s,nobs,gamma,icor,nmod,nd)
    IMPLICIT NONE
    ! number of observations
    INTEGER, INTENT(IN) :: nobs
    ! observations, new data types
    TYPE(ast_obs),DIMENSION(nobs), INTENT(IN) :: obs_s
    ! obs. weights and residduals, new data types
    TYPE(ast_wbsr),DIMENSION(nobs), INTENT(INOUT) :: obsw_s
    ! norm of residuals, covariance matrix
    DOUBLE PRECISION, INTENT(IN) ::  csinor,gamma(nd,nd)
    ! flags indicating the variables being solved for
    INTEGER, INTENT(IN) :: icor(nd)
    ! unit for output
    INTEGER, INTENT(IN) :: unilog
    ! no. solve for parameters
    INTEGER, INTENT(IN) :: nd
    ! hidden input
    ! ===========units for err,pro,clo files ========
    ! OUTPUT
    ! rejection mode
    INTEGER, INTENT(OUT):: nmod
    ! END INTERFACE  =======================
    ! Output of debugging information
    LOGICAL, PARAMETER :: debug=.false.
    ! variables to copy observation fields
    INTEGER sel(nobs)
    CHARACTER*1 type(nobs)
    DOUBLE PRECISION csi(2),w(3),x2(nobs),det_w,w_inv(3)
    DOUBLE PRECISION mjd(nobs)
    !
    INTEGER iob,ipa,ipd,i,k,nsel,minobs,selp,nrej,nrec
    DOUBLE PRECISION sa,sd,sad,gamga(ndimx),gamgd(ndimx),covres(2,2),wres(2,2)
    DOUBLE PRECISION det,x2max
    DOUBLE PRECISION covr,wrr
    ! calendar date variables
    INTEGER iyear,imonth,iday
    DOUBLE PRECISION hour
    CHARACTER*16 date
    ! fudge factor, stuck loops
    DOUBLE PRECISION fudge
    INTEGER nit,oldmod,oldnob,iunlog
    SAVE nit,oldmod,oldnob
    DATA nit/0/
    DATA oldmod/-1000000/
    DATA oldnob/-1000000/
    IF(iicrej.NE.36) STOP '**** reject: internal error (01) ****'
    !     log=.true.
    !     IF(unilog.le.0) log=.false.
    nrej=0
    nrec=0
    iunlog=abs(unilog)
    ! copy selection flags, types
    sel=obsw_s%sel_coord
    type=obs_s%type
    mjd=obs_s%time_utc
    !****************************************************
    ! First compute chi**2 for each obs
    !****************************************************
    ! Minimum number of observations
    !     First count solved for params
    minobs=nd
    DO i=1,nd
       IF(icor(i).NE.1) minobs=minobs-1
    END DO
    !     Next set minobs
    minobs=nint(0.5d0*minobs)
    IF(nobs.le.minobs)THEN
       IF(verb_rej.gt.9)WRITE(iunlog,*)                                 &
            &        'No autorejection with only ',nobs,' observations.'
       IF(verb_rej.ge.5)WRITE(*,*)                                    &
            &        'No autorejection with only ',nobs,' observations.'
       nmod=0
       RETURN
    ENDIF
    !
    x2max=0
    nsel=0
    IF(verb_rej.ge.9)WRITE(iunlog,300) csinor
    IF(verb_rej.ge.19) WRITE(*,300) csinor
300 FORMAT('==== OUTLIER REJECTION ======='/                          &
         &        'Residual norm =',1P,E13.5,0P)
    DO 1 iob=1,nobs
       ! Pointers to RA and DEC observations
       ipa=2*iob-1
       ipd=2*iob
       ! Expected variance of fit residuals
       DO 2 i=1,nd
          sa=0
          sd=0
          DO 3 k=1,nd
             sa=sa+gamma(i,k)*g(ipa,k)
             sd=sd+gamma(i,k)*g(ipd,k)
3         ENDDO
          gamga(i)=sa
          gamgd(i)=sd
2      ENDDO
       ! find weights and residuals for the two components of this obs.
       ! they have to be available even when sel(iob)=0
       CALL fit_weight(obs_s(iob),obsw_s(iob),.false.,w)
       csi=obsw_s(iob)%res_coord
       ! check obs. type
       IF(type(iob).eq.'O'.or.type(iob).eq.'S')THEN
          ! both components to be used
          sa=0
          sd=0
          sad=0
          DO  i=1,nd
             sa =sa +g(ipa,i)*gamga(i)
             sd =sd +g(ipd,i)*gamgd(i)
             sad=sad+g(ipa,i)*gamgd(i)
          ENDDO

          det_w=w(1)*w(2)-w(3)**2
          w_inv(1)=w(2)/det_w
          w_inv(2)=w(1)/det_w
          w_inv(3)=-w(3)/det_w
          ! separate case with zero weight in one of the two components
          IF(sel(iob).eq.0)THEN
             covres(1,1)=w_inv(1)+sa
             covres(2,2)=w_inv(2)+sd
             covres(1,2)=w_inv(3)+sad
          ELSE
             covres(1,1)=w_inv(1)-sa
             covres(2,2)=w_inv(2)-sd
             covres(1,2)=w_inv(3)-sad
          END IF
          covres(2,1)=covres(1,2)
          ! Chi-square value
          CALL inv22(covres,wres,det)
          IF (det.eq.0.d0)THEN
             WRITE(ierrou,*) 'WARNING: reject.f, det=0.'
             IF(verb_rej.gt.9)WRITE(iunlog,*) 'reject.f: det=0.'
             IF(verb_rej.gt.19)WRITE(*,*) 'reject.f: det=0.'
             numerr=numerr+1
             x2(iob)=0.d0
          ELSE
             x2(iob)=csi(1)*wres(1,1)*csi(1)+csi(2)*wres(2,2)*csi(2)+2*csi(1)*wres(1,2)*csi(2)
          !   IF(iob.eq.48) THEN
          !      WRITE(*,*) x2(iob), csi, wres
          !   ENDIF
          ENDIF
       ELSEIF(type(iob).eq.'V')THEN
          ! only second component to be used
          sd=0
          DO  i=1,nd
             sa =sa +g(ipa,i)*gamga(i)
             sd =sd +g(ipd,i)*gamgd(i)
             sad=sad+g(ipa,i)*gamgd(i)
          ENDDO
          IF(sel(iob).EQ.0) THEN
             covr=1/w(2)+sd
          ELSE
             covr=1/w(2)-sd
          END IF
          ! Chi-square value
          IF (covr.eq.0.d0)THEN
             WRITE(ierrou,*) 'WARNING: reject.f, covr=0. (1)'
             numerr=numerr+1
          ENDIF
          wrr=1.d0/covr
          x2(iob)=csi(2)*wrr*csi(2)
       ELSEIF(type(iob).eq.'R')THEN
          ! only first component to be used
          sa=0
          DO  i=1,nd
             sa =sa +g(ipa,i)*gamga(i)
          ENDDO
          IF(sel(iob).EQ.0) THEN
             covr=1/w(1)+sa
          ELSE
             covr=1/w(1)-sa
          END IF
          ! Chi-square value
          IF (covr.eq.0.d0)THEN
             WRITE(ierrou,*) 'WARNING: reject.f, covr=0. (2)'
             numerr=numerr+1
          ENDIF
          wrr=1.d0/covr
          x2(iob)=csi(1)*wrr*csi(1)
       ELSE
          WRITE(*,*)'reject: this should not happen ',                &
               &           w(ipa),' ',w(ipd)
          STOP 'reject: w '
       ENDIF
       IF (x2(iob).lt.0.d0)THEN
          IF(ierrou.le.0)WRITE(*,*) 'reject_obs: ierrou ', ierrou,&
               &       ' object name ', obs_s(1)%objdes
          WRITE(ierrou,*) 'WARNING: reject.f, x2=',x2(iob)
          If(verb_rej.ge.15)WRITE(*,*) 'WARNING: reject.f, x2=',x2(iob)
          numerr=numerr+1
          x2(iob)=0.d0
       ENDIF
       ! X tech observations are not serious ones, thus it does not matter
       ! if they have a high chi^2
       IF(sel(iob).NE.0.and.obs_s(iob)%tech.ne.'X') THEN
          x2max=MAX(x2(iob),x2max)
          nsel=nsel+1
       END IF
1   END DO
    !****************************************************
    ! Second perform tests
    !****************************************************
    IF(nsel.LE.0) STOP '**** reject: internal error (02) ****'
    !     For <=18 obs. or for inf. loop -> reject one at a time
    IF(nsel .le. minobs*6) THEN
       IF(unilog.gt.0.and.verb_rej.gt.9) WRITE(unilog,*) 'Rejecting no more than one.'
    ELSEIF(nit.ge.4) THEN
       IF(unilog.gt.0.and.verb_rej.gt.9) WRITE(unilog,*) 'Trying to get unstuck.'
    ELSE
       x2max=x2max*x2frac
    ENDIF
    IF(verb_rej.ge.9) WRITE(iunlog,302) x2max
    IF(verb_rej.ge.19) WRITE(*,302) x2max
302 FORMAT('Chi2 threshold =',1P,E13.5)
    ! This is a "fudge" factor to make it more difficult to reject when
    ! there are very few observations. If there are more than ~50 then
    ! this has little effect.
    IF(rej_fudge)THEN
       IF(error_model_priv.eq.' ')THEN
          fudge=400.d0 * (1.2)**(-nsel)
       ELSE
          ! if there is an error model, then the fudge has to be used
          ! only for less than ~8 observations.
          fudge=400.d0 * (3)**(-nsel)
       ENDIF
    ELSE
       fudge=0.d0
    ENDIF
    DO 5 iob=1,nobs
       selp=sel(iob)
       IF(obs_s(iob)%tech.eq.'X')THEN
          IF(sel(iob).gt.0)THEN
             call mjddat(mjd(iob),iday,imonth,iyear,hour)
             write(date,'(i4,a1,i2.2,a1,f8.5)')                       &
                  &              iyear,'/',imonth,'/',iday+hour/24d0
             WRITE(ierrou,*)'rejecting X observation',date,x2(iob),obs_s(iob)%obscod_s
             numerr=numerr+1
          ENDIF
          sel(iob)=0
       ELSEIF(selp.EQ.0) THEN
          IF(x2(iob).LE.x2rec + 0.75*fudge) THEN
             call mjddat(mjd(iob),iday,imonth,iyear,hour)
             write(date,'(i4,a1,i2.2,a1,f8.5)')                       &
                  &              iyear,'/',imonth,'/',iday+hour/24d0
             IF(verb_rej.ge.9)write(iunlog,*) 'Recover: ',date,x2(iob),type(iob)
             IF(verb_rej.ge.29)write(*,*) 'Recover: ',date,x2(iob),type(iob)
             sel(iob)=1
             nrec=nrec+1
          ENDIF
       ELSE
          IF(x2(iob).GE.x2max .AND. x2(iob).GT.x2rej + fudge) THEN
             ! check to see if this is the only remaining obs in an opposition
             call mjddat(mjd(iob),iday,imonth,iyear,hour)
             write(date,'(i4,a1,i2.2,a1,f8.5)')                       &
                  &              iyear,'/',imonth,'/',iday+hour/24d0
             IF((type(iob).eq.'R'.or.type(iob).eq.'V').and.         &
                  &             entopp_radar(iob,sel,mjd,type,nobs))THEN
                write(*,*)'*** WARNING: ',name_obj(1:lobjnam),       &
                     &              ' Did not reject last radar obs. in opp. Date:',date
                write(ierrou,*)'*** WARNING: ',name_obj(1:lobjnam),  &
                     &              ' Did not reject last radar obs. in opp. Date:',date
                numerr=numerr+1
             ELSEIF(entopp(iob,sel,mjd,nobs))THEN
                if(rejopp)then
                   write(*,*)                                         &
                        &                '*** WARNING: ',name_obj(1:lobjnam),              &
                        &               ' Rejecting last obs. in opp. Date:', date
                   write(ierrou,*)                                    &
                        &                '*** WARNING: ',name_obj(1:lobjnam),              &
                        &               ' Rejecting last obs. in opp. Date:', date
                   numerr=numerr+1
                   sel(iob)=0
                   nrej=nrej+1
                else
                   write(*,*)                                         &
                        &            '*** WARNING: ',name_obj(1:lobjnam),                   &
                        &            ' Did not reject last obs. in opp. Date:',date
                   write(ierrou,*)                                    &
                        &            '*** WARNING: ',name_obj(1:lobjnam),                   &
                        &            ' Did not reject last obs. in opp. Date:',date
                   numerr=numerr+1
                endif
             ELSE
                IF(verb_rej.ge.9)write(iunlog,*) 'Reject : ',date,x2(iob),type(iob)
                IF(verb_rej.ge.29)write(*,*) 'Reject : ',date,x2(iob),type(iob)
                sel(iob)=0
                nrej=nrej+1
             ENDIF
          END IF
       END IF

       IF(debug.and.unilog.gt.0)                                      &
            &        WRITE(unilog,301) iob,x2(iob),selp,sel(iob)
5   END DO
301 FORMAT(I5,1P,E13.5,0P,2I3)
    IF(verb_rej.ge.9) WRITE(iunlog,303) nobs,nsel,nrej,nrec
    IF(verb_rej.ge.19) WRITE(*,303) nobs,nsel,nrej,nrec
303 FORMAT('No. of observations:'/                                    &
         &       10X,'total     =',I5/                                      &
         &       10X,'selected  =',I5/                                      &
         &       10X,'rejected  =',I5/                                      &
         &       10X,'recovered =',I5)

    nmod=nrej+nrec
    !     Is it the same no. of modifications with the same number of obs?
    IF(oldmod.EQ.nmod .AND. nobs.EQ.oldnob)THEN
       nit=nit+1
    ELSE
       nit=0
    ENDIF
    oldmod=nmod
    oldnob=nobs
    ! copy x2, sel into chi of obsw_s
    obsw_s(1:nobs)%chi=sqrt(x2(1:nobs))
    obsw_s%sel_coord=sel

  CONTAINS

    ! END SUBROUTINE reject_obs

    ! decide if index is the last selected obs in an opposition
    LOGICAL FUNCTION entopp(index,sel,mjd,nobs)
      integer index,nobs,sel(nobs),i
      double precision mjd(nobs)
      entopp=.false.
      do i=1,nobs
         ! return if we find a different selected obs that is close enough
         if(i.ne.index.and.                                             &
              &        abs(mjd(i)-mjd(index)).le.180d0.and.                      &
              &        sel(i).eq.1) return
      enddo
      entopp=.true.
    END FUNCTION entopp
    ! decide if index is the last selected radar obs in an opposition
    LOGICAL FUNCTION entopp_radar(index,sel,mjd,type,nobs)
      integer index,nobs,sel(nobs),i
      double precision mjd(nobs)
      CHARACTER*1 type(nobs)
      entopp_radar=.false.
      do i=1,nobs
         ! return if we find a different selected obs that is close enough
         if(i.ne.index.and.                                             &
              &        abs(mjd(i)-mjd(index)).le.180d0.and.                      &
              &        sel(i).eq.1.and.(type(i).eq.'R'.or.type(i).eq.'V')) return
      enddo
      entopp_radar=.true.
    END FUNCTION entopp_radar

  END SUBROUTINE reject_obs
  !
  ! Copyright 1999 Orbfit Consortium
  ! Written by Milani and Chesley
  ! Last Modified 20/1/99 by S. Chesley
  ! AM, FS 2/7/2013 output unormalized rms
  !***********************************************************************
  ! ROUTINE MAGEST
  !***********************************************************************
  ! INPUTS:  from obs, obsw
  !          smag -     magnitude strings e.g., '15.55V'
  !          rmsmag -   a priori rms of the observation (NOTE: rmsmag will
  !                        be changed to 99.99 for outlier rejection)
  !          sel -      selection flag for the astrometric observation
  !          m -        total number of observations
  ! OUTPUTS: into obsw
  !          resmag -   residual of the mag observation (NOTE: If resmag =
  !                        1.d9 then the mag obs was not used in the fit.)
  !          h0 -       absolute magnitude
  !          rmsh -     weighted RMS of the residuals
  !          rmshout -  unweighted RMS of the residuals in magnitudes
  !          nmag    -  number of magnitude observations
  !***********************************************************************
  SUBROUTINE mag_est(m,obs,obsw,h0,rmsh,rmshout,nmag)
    IMPLICIT NONE
    ! INPUT
    ! number of observations
    INTEGER m
    ! new data types
    TYPE(ast_obs),DIMENSION(m), INTENT(IN) :: obs
    ! INPUT/OUTPUT
    TYPE(ast_wbsr),DIMENSION(m), INTENT(INOUT) :: obsw
    ! OUPUT (but may exist before)
    ! magnitudes: estimated values, fit rms
    DOUBLE PRECISION, INTENT(INOUT) :: h0, rmsh
    DOUBLE PRECISION, OPTIONAL, INTENT(OUT) :: rmshout
    INTEGER, OPTIONAL, INTENT(OUT) :: nmag
    ! END INTERFACE =======================
    ! magnitudes: observed (string), a priori rms
    CHARACTER*6 smag
    ! observed magnitudes, mean
    DOUBLE PRECISION obsm(nobx),hmean,hsum,wsum,wrsum
    DOUBLE PRECISION hnew,chi2, rmsh1
    INTEGER isum
    ! fit residuals, a priori rms, selection flag
    DOUBLE PRECISION resmag(nobx),rmsmag(nobx)
    INTEGER selmag(nobx)
    LOGICAL avail(nobx)
    CHARACTER*1 col
    CHARACTER*5 smag5
    INTEGER numrej,numchg
    ! debugging:
    LOGICAL verbose
    ! common magnitude data
    ! function to compute default mag RMS
    DOUBLE PRECISION magrms
    ! loop indexes
    INTEGER j,icount,nrej
    ! copy from observation and obs. weights
    ! ======================================================================
    verbose = .false.
    ! First create a vector of appmags and rms's
    DO j=1,m
       smag=obs(j)%mag_str
       READ(smag,101)col
101    FORMAT(5x,a1)
!!! maybe this is done in input.....
       !         obs(j)%mag_band=col
       avail(j)=smag.ne.'      '.and.                              &
            & (col.eq.'B'.or.col.eq.'V'.or.col.eq.'R'.or.col.eq.'U'.or.col.eq.'J'.or. &
            &  col.eq.'I'.or.col.eq.'C'.or.col.eq.'z'.or.col.eq.'i'.or.col.eq.'g'.or. &
            &  col.eq.'G'.or.col.eq.'H'.or.col.eq.' '.or.col.eq.'r'.or.col.eq.'y'.or. &
            &  col.eq.'w'.or.col.eq.'K'.or.col.eq.'L'.or.col.eq.'Y'.or.col.eq.'c'.or.col.eq.'o'.or.col.eq.'u')
       IF(.not.(col.eq.'B'.or.col.eq.'V'.or.col.eq.'R'.or.col.eq.'U'.or.col.eq.'J'    &
            &   .or.col.eq.'I'.or.col.eq.'C'.or.col.eq.'z'.or.col.eq.'i'.or.col.eq.'g' &
            &   .or.col.eq.'G'.or.col.eq.'H'.or.col.eq.' '.or.col.eq.'r'.or.col.eq.'y' &
            &   .or.col.eq.'w'.or.col.eq.'K'.or.col.eq.'L'.or.col.eq.'Y'.or.col.eq.'c'.or.col.eq.'o'.or.col.eq.'u'))THEN
          IF(ierrou.gt.0)THEN
             WRITE(ierrou,*) 'Unknown Color:',col,' Obs #',j
             numerr=numerr+1
          ELSE
             WRITE(*,*) 'Unknown Color:',col,' Obs #',j
          ENDIF
       ENDIF
!!! maybe this is done in input.....
       !         obs(j)%mag_def=avail(j)
!!!! fix
       ! prepare vector of rms
       rmsmag(j)=obsw(j)%rms_mag
       IF(avail(j))THEN
          ! Get the default rms (to be used for recovery)
          ! The last three arguments are dummies.
          !            tmprms(j)=magrms(smag,0d0,0,'0')
          ! Read the magnitude and make color corrections.
          smag5=smag(1:5)
          READ(smag5,*,ERR=3)obsm(j)
          ! to remove correction to AB magnitude
          IF(ab_mag)THEN
             IF(col.eq.'g')THEN
                obsm(j)=obsm(j)+0.093d0
             ELSEIF(col.eq.'r')THEN
                obsm(j)=obsm(j)-0.166d0
             ELSEIF(col.eq.'i')THEN
                obsm(j)=obsm(j)-0.397d0
             ELSEIF(col.eq.'z')THEN
                obsm(j)=obsm(j)-0.572d0
             ENDIF
          ENDIF
          ! conversion to V magnitude
          IF(col.eq.'V'.or.col.eq.'H')THEN
             CONTINUE
          ELSEIF(col.eq.' ')THEN
             obsm(j)=obsm(j)-0.8d0
          ELSEIF(col.eq.'B')THEN
             obsm(j)=obsm(j)-0.8d0
          ELSEIF(col.eq.'R')THEN
             obsm(j)=obsm(j)+0.4d0
          ELSEIF(col.eq.'I')THEN
             obsm(j)=obsm(j)+0.8d0
          ELSEIF(col.eq.'C')THEN
             obsm(j)=obsm(j)+0.4d0
             ! the following values for g,r,i,z,y,w are provided by
             ! Fitzsimmons November 2011 report
          ELSEIF(col.eq.'z')THEN
             obsm(j)=obsm(j)+0.37d0 ! -0.1d0
          ELSEIF(col.eq.'i')THEN
             obsm(j)=obsm(j)+0.39d0 ! +0.3d0
          ELSEIF(col.eq.'g')THEN
             obsm(j)=obsm(j)-0.28d0 ! -0.5d0
          ELSEIF(col.eq.'U'.or.col.eq.'J')THEN
             obsm(j)=obsm(j)+1.1d0
          ELSEIF(col.eq.'r')THEN
             obsm(j)=obsm(j)+0.23d0 ! -0.1d0
          ELSEIF(col.eq.'y')THEN !
             obsm(j)=obsm(j)+0.36d0 ! +0.0d0
          ELSEIF(col.eq.'w')THEN ! this is a conjecture, no data
             obsm(j)=obsm(j)+0.16d0
          ELSEIF(col.eq.'G')THEN
             obsm(j)=obsm(j)+0.28d0 ! Gaia correction by Tholen
          ELSEIF(col.eq.'c')THEN
             obsm(j)=obsm(j)-0.05d0
          ELSEIF(col.eq.'o')THEN
             obsm(j)=obsm(j)+0.33d0
          ELSEIF(col.eq.'L')THEN
             obsm(j)=obsm(j)+0.2d0
         ! ELSEIF(col.eq.'H')THEN
         !    obsm(j)=obsm(j)+1.4d0
          ELSEIF(col.eq.'K')THEN
             obsm(j)=obsm(j)+1.7d0
          ELSEIF(col.eq.'Y')THEN
             obsm(j)=obsm(j)+0.7d0
          ELSEIF(col.eq.'u')THEN
             obsm(j)=obsm(j)+2.5d0
          ENDIF
       ENDIF
    ENDDO
    ! should be 1 by default
    selmag(1:m)=obsw%sel_mag
    ! Start Outlier loop
    DO nrej=1,10
       hsum=0.d0
       icount=0
       wsum=0.d0
       ! main loop
       DO j=1,m
          IF(avail(j))THEN
             icount=icount+1
             hsum=hsum+(obsm(j)-dmagn(j))/rmsmag(j)
             wsum=wsum+1.d0/rmsmag(j)
          ENDIF
       ENDDO
       ! jump out if no good observations
       IF(icount.eq.0)THEN
          IF(verb_dif.ge.15)WRITE(*,*)' magnitude not estimated w/o obs'
          IF(PRESENT(rmshout))THEN
             rmshout=-1.d0
          ENDIF
          h0=-9.99d0
          RETURN
       ENDIF
       ! Compute H, etc.
       hmean=hsum/wsum
       rmsh=0.d0
       rmsh1=0.d0
       wrsum=0.d0
       isum=0
       DO j=1,m
          IF(avail(j))THEN
             resmag(j)=obsm(j)-dmagn(j)-hmean
             rmsh1=rmsh1+resmag(j)**2
             rmsh=rmsh+(resmag(j)/rmsmag(j))**2
             isum=isum+1
             wrsum=wrsum+(1.d0/rmsmag(j))**2
             selmag(j)=1
          ELSE
             !              if data not used
             selmag(j)=0
          ENDIF
       ENDDO
       ! Divide by icount for rmshout, by wsum for rmsh
       rmsh=sqrt(rmsh/wrsum)
       rmsh1=sqrt(rmsh1/isum)
       IF(PRESENT(rmshout).and.PRESENT(nmag))THEN
          rmshout=rmsh1
          nmag=isum
       ENDIF
       IF(verbose) write(*,*)'H,RMS',hmean,rmsh
       ! Handle outliers in a separate loop
       numrej=0
       numchg=0
       ! Never reject if there are less than 5 obs.
       IF(icount.le.4) GOTO 10
       DO j=1,m
          IF(avail(j))THEN
             IF(selmag(j).eq.0)THEN
                !              test for recovery
                hnew=(hsum+(obsm(j)-dmagn(j))/rmsmag(j))/             &
                     &                 (wsum+1.d0/rmsmag(j))
                chi2=(((obsm(j)-dmagn(j))-hnew)/rmsmag(j))**2
                IF(chi2.lt.x2mrec)THEN
                   numchg=numchg+1
                   selmag(j)=1
                   !                     rmsmag(j)=tmprms(j)
                ENDIF
                IF(verbose) write(*,*)'REC:Hnew,chi2,rms,res',        &
                     &                 hnew,chi2,rmsmag(j),resmag(j)
             ELSE
                !              test for rejection
                hnew=(hsum-(obsm(j)-dmagn(j))/rmsmag(j))/             &
                     &                 (wsum-1.d0/rmsmag(j))
                chi2=(((obsm(j)-dmagn(j))-hnew)/rmsmag(j))**2
                IF(chi2.gt.x2mrej)THEN
                   numchg=numchg+1
                   selmag(j)=0
                ENDIF
                IF(verbose) write(*,*)'REJ:Hnew,chi2,rms,res',        &
                     &                 hnew,chi2,rmsmag(j),resmag(j)
             ENDIF
             !              count rejection
             IF(selmag(j).eq.0) numrej=numrej+1
          ENDIF
       ENDDO
       IF(numchg.eq.0) GOTO 10
    ENDDO
    If(verb_dif.ge.5)                                                 &
         &     WRITE(iun_log,*)'WARNING: Did not finish rejecting photometry...'
    !***********************************************************************
10  IF(verb_dif.ge.9)WRITE(iun_log,190)h0,hmean,icount,rmsh,numrej,nrej
    IF(verb_dif.ge.19)WRITE(*,190)h0,hmean,icount,rmsh,numrej,nrej
190 FORMAT('Absol Magnitude: Old=',f5.2,' New=',f5.2,/             &
         &        'no. magnitude obs.=',i4, ' with RMS ',f5.2,/             &
         &        'no. rejected obs.=',i4,' in ',i4,' passes.')
    IF(verb_dif.ge.19) WRITE(*,*)
    ! output data: to calling program
    h0=hmean
    ! to obsw array
    DO j=1,m
       IF(avail(j))THEN
          obsw(j)%res_mag=resmag(j)
          obsw(j)%resm_def=.true.
          obsw(j)%sel_mag=selmag(j)
       ELSE
          obsw(j)%resm_def=.false.
       ENDIF
    ENDDO
    RETURN
    ! parsing error
3   If(verb_dif.ge.5)                                                 &
         &     WRITE(*,*)'magest: error in parsing magn. no ',j,' ',smag
    STOP
  END SUBROUTINE mag_est
  ! ========================================================
  ! WRITE_PHOTOM output of data for statistics of photometry
  ! ========================================================
  SUBROUTINE write_photom(iupho,name0,h0,gma0,obs,obsw,m)
    USE name_rules, ONLY: name_len
    INTEGER, INTENT(IN) :: iupho ! output unit
    CHARACTER*(name_len), INTENT(IN) :: name0 ! asteroid name
    DOUBLE PRECISION, INTENT(IN) :: h0,gma0 ! H and G magnitudes as estimated
    ! ===== observational data ===========================
    INTEGER m ! observation number
    ! observations new data type
    TYPE(ast_obs), INTENT(IN), DIMENSION(m) :: obs
    TYPE(ast_wbsr), INTENT(IN),DIMENSION(m) :: obsw
    ! end interface
    INTEGER ii
    WRITE(iupho,107)name0,h0,gma0
107 FORMAT('!',a18,1x,f7.4,1x,f6.3)
    DO ii=1,m
       IF(obs(ii)%mag_str.ne.'      ')WRITE(iupho,117)        &
            &                obs(ii)%obscod_s,obs(ii)%time_utc,                 &
            &          obs(ii)%mag_str,dmagn(ii)+h0,phaa(ii)*degrad,dsuna(ii),  &
            &      disa(ii),obs(ii)%coord(1)*degrad,obs(ii)%coord(2)*degrad,   &
            &                adotmv(ii)*degrad,ddotmv(ii)*degrad
117    FORMAT(i3.3,1x,f13.6,1x,a6,1x,f8.4,1x,f8.4,1x,f7.4,    &
            &                1x,f7.4,1x,f9.5,1x,f9.5,1x,f8.4,1x,f8.4)
    ENDDO
  END SUBROUTINE write_photom
  ! Copyright (C) 1998 by Mario Carpino (carpino@brera.mi.astro.it)
  ! Version: December 14, 1998 Steve Chesley
  ! ---------------------------------------------------------------------
  !
  !  *****************************************************************
  !  *                                                               *
  !  *                         R E J I N I                           *
  !  *                                                               *
  !  *        Initialization of options for outlier rejection        *
  !  *                                                               *
  !  *****************************************************************
  !
  SUBROUTINE rejini
    LOGICAL found,fail1,fail
    fail=.false.
    autrej=.true.
    CALL rdnlog('reject.','auto',autrej,.false.,found,                &
         &            fail1,fail)
    rejopp=.false.
    CALL rdnlog('reject.','rejopp',rejopp,.false.,found,              &
         &            fail1,fail)
    rej_fudge=.true.
    CALL rdnlog('reject.','rej_fudge',rej_fudge,.false.,found,  &
         &            fail1,fail)
    x2rej=10.0D0
    CALL rdnrea('reject.','chi2_reject',x2rej,.false.,found,          &
         &            fail1,fail)
    x2rec=9.21D0
    CALL rdnrea('reject.','chi2_recover',x2rec,.false.,found,         &
         &            fail1,fail)
    x2frac=0.25D0
    CALL rdnrea('reject.','chi2_frac',x2frac,.false.,found,           &
         &            fail1,fail)
    delrej=1.D-1
    CALL rdnrea('reject.','conv_cntr',delrej,.false.,found,           &
         &            fail1,fail)
    itmaxr=15
    CALL rdnint('reject.','nit_max',itmaxr,.false.,found,             &
         &            fail1,fail)
    fomax=50.D0
    CALL rdnrea('reject.','max_perc',fomax,.false.,found,             &
         &            fail1,fail)
    fomax=fomax/100
    ! Magnitude chi^2 values
    x2mrej=8.0D0
    CALL rdnrea('reject.','chi2_mag_rej',x2mrej,.false.,found,        &
         &            fail1,fail)
    x2mrec=7.0D0
    CALL rdnrea('reject.','chi2_mag_rec',x2mrec,.false.,found,        &
         &            fail1,fail)
    IF(fail) STOP '**** rejini: abnormal end ****'

    iicrej=36
  END SUBROUTINE rejini
  ! Copyright (C) 1998 by Mario Carpino (carpino@brera.mi.astro.it)
  ! Version: December 3, 1998
  ! ---------------------------------------------------------------------
  !
  !  *****************************************************************
  !  *                                                               *
  !  *                         D I F I N I                           *
  !  *                                                               *
  !  *        Initialization of options for iteration control        *
  !  *                     in ROUTINE DIFCOR                         *
  !  *                                                               *
  !  *****************************************************************
  !
  SUBROUTINE difini
    LOGICAL found,fail1,fail
    CHARACTER*80 :: comment
    fail=.false.
    itmax=20
    CALL rdnint('difcor.','nit_max',itmax,.false.,found,fail1,fail)
    itgmax=5
    CALL rdnint('difcor.','nitg_max',itgmax,.false.,found,fail1,fail)
    divrat=0.999d0
    CALL rdnrea('difcor.','div_cntr',divrat,.false.,found,fail1,fail)
    delcr=1.d-3
    CALL rdnrea('difcor.','conv_contr',delcr,.false.,found,fail1,fail)
    del_constr=1.d-5
    CALL rdnrea('difcor.','del_constr',del_constr,.false.,found,fail1,fail)
    step_sig=0.5d0
    CALL rdnrea('difcor.','step_sig',step_sig,.false.,found,fail1,fail)
    scaling_lov=.false.
    CALL rdnlog('difcor.','scaling_lov',scaling_lov,.false.,found,fail1,fail)
    second_lov=.false.
    CALL rdnlog('difcor.','second_lov',second_lov,.false.,found,fail1,fail)
    IF(fail) STOP '**** difini: abnormal end ****'
    iicdif=36
    ! default not in key/def file ???? there is output_des
    output_des=.false.
    output_failrwo=.true.
    ! a priori for full matrix; to be set in known_fast main (added AM 4/1/2016)
    aprc_flag=.false.
    ! =========================================
    ! A priori for non-gravitational parameters
    ! =========================================
    ! A priori for Yarkovsky (added 11/4/2014)
    ! A priori for all (added 3/2/2016, A. Del Vigna)
    IF(ngr_opt) THEN
       ! Inizialization
       yarko_apriori=.FALSE.
       drpa2m_apriori=.FALSE.
       yarko_nominal=0.d0
       drpa2m_nominal=0.d0
       a1ng_apriori=.FALSE.
       a1ng_nominal=0.d0
       a2ng_apriori=.FALSE.
       a2ng_nominal=0.d0
       a3ng_apriori=.FALSE.
       a3ng_nominal=0.d0
       dtdelay_apriori=.FALSE.
       dtdelay_nominal=0.d0
       IF(dyn%nmod.EQ.1)THEN
          ! Area/mass
          comment='direct radiation pressure a-priori'
          CALL input_log_opt('difcor.','drpa2m_apriori',drpa2m_apriori,.FALSE.,found,comment,iun_log)
          dyn%apriori(1)=drpa2m_apriori
          IF(dyn%apriori(1))THEN
             comment='a-priori nominal value for drp'
             CALL input_rea_opt('difcor.','drpa2m_nominal',drpa2m_nominal,.FALSE.,found,comment,iun_log)
             dyn%apriori_nominal(1)=drpa2m_nominal
             CALL rdnrea('difcor.','drpa2m_std',drpa2m_std,.false.,found,fail1,fail)
             IF(.not.found.and.drpa2m_apriori) THEN
                WRITE(*,*) '******* difini: abnormal end due to drpa2m_std *******'
                WRITE(*,*) 'ERROR: drpa2m_apriori is true, but drpa2m_std is missing'
                IF(iun_log.NE.0) WRITE(iun_log,*) 'ERROR: drpa2m_apriori is true, but drpa2m_std is missing'
                STOP
             END IF
             dyn%apriori_std(1)=drpa2m_std
          END IF
          ! Yarkovsky
          comment='Yarkovsky a-priori'
          CALL input_log_opt('difcor.','yarko_apriori',yarko_apriori,.FALSE.,found,comment,iun_log)
          dyn%apriori(2)=yarko_apriori
          IF(dyn%apriori(2))THEN
             comment='a-priori nominal value for Yarkovsky parameter'
             CALL input_rea_opt('difcor.','yarko_nominal',yarko_nominal,.FALSE.,found,comment,iun_log)
             dyn%apriori_nominal(2)=yarko_nominal
             CALL rdnrea('difcor.','yarko_std',yarko_std,.false.,found,fail1,fail)
             IF(.not.found.and.yarko_apriori) THEN
                WRITE(*,*) '******* difini: abnormal end due to yarko_std *******'
                WRITE(*,*) 'ERROR: yarko_apriori is true, but yarko_std is missing'
                IF(iun_log.NE.0) WRITE(iun_log,*) 'ERROR: yarko_apriori is true, but yarko_std is missing';
                STOP
             END IF
             dyn%apriori_std(2)=yarko_std
          END IF
       ELSE IF(dyn%nmod.EQ.2) THEN
          ! Outgassing a1ng
          comment='outgassing A1 a-priori'
          CALL input_log_opt('difcor.','a1ng_apriori',a1ng_apriori,.FALSE.,found,comment,iun_log)
          dyn%apriori(1)=a1ng_apriori
          IF(dyn%apriori(1))THEN
             comment='a-priori nominal value for the A1 parameter'
             CALL input_rea_opt('difcor.','a1ng_nominal',a1ng_nominal,.FALSE.,found,comment,iun_log)
             dyn%apriori_nominal(1)=a1ng_nominal
             CALL rdnrea('difcor.','a1ng_std',a1ng_std,.false.,found,fail1,fail)
             IF(.not.found.and.a1ng_apriori) THEN
                WRITE(*,*) '******* difini: abnormal end due to a1ng_std *******'
                WRITE(*,*) 'ERROR: a1ng_apriori is true, but a1ng_std is missing'
                IF(iun_log.NE.0) WRITE(iun_log,*) 'ERROR: a1ng_apriori is true, but a1ng_std is missing'
                STOP
             END IF
             dyn%apriori_std(1)=a1ng_std
          END IF
          ! Outgassing a2ng
          comment='outgassing A2 a-priori'
          CALL input_log_opt('difcor.','a2ng_apriori',a2ng_apriori,.FALSE.,found,comment,iun_log)
          dyn%apriori(2)=a2ng_apriori
          IF(dyn%apriori(2))THEN
             comment='a-priori nominal value for the A2 parameter'
             CALL input_rea_opt('difcor.','a2ng_nominal',a2ng_nominal,.FALSE.,found,comment,iun_log)
             dyn%apriori_nominal(2)=a2ng_nominal
             CALL rdnrea('difcor.','a2ng_std',a2ng_std,.false.,found,fail1,fail)
             IF(.not.found.and.a2ng_apriori) THEN
                WRITE(*,*) '******* difini: abnormal end due to a2ng_std *******'
                WRITE(*,*) 'ERROR: a2ng_apriori is true, but a2ng_std is missing'
                IF(iun_log.NE.0) WRITE(iun_log,*) 'ERROR: a2ng_apriori is true, but a2ng_std is missing'
                STOP
             END IF
             dyn%apriori_std(2)=a2ng_std
          END IF
          ! Check: if A1 a-priori is given, it must be given also A2 a-priori (and viceversa)
          IF(a1ng_apriori.AND..NOT.a2ng_apriori)THEN
             WRITE(*,*) '******* difini: abnormal end due to a1ng_apriori and a2ng_apriori *******'
             WRITE(*,*) 'ERROR: a1ng_apriori is true, but a2ng_apriori is false'
             STOP
          END IF
          IF(.NOT.a1ng_apriori.AND.a2ng_apriori)THEN
             WRITE(*,*) '******* difini: abnormal end due to a1ng_apriori and a2ng_apriori *******'
             WRITE(*,*) 'ERROR: a2ng_apriori is true, but a1ng_apriori is false'
             STOP
          END IF
          ! Outgassing a3ng
          comment='outgassing A3 a-priori'
          CALL input_log_opt('difcor.','a3ng_apriori',a3ng_apriori,.FALSE.,found,comment,iun_log)
          dyn%apriori(3)=a3ng_apriori
          IF(dyn%apriori(3))THEN
             comment='a-priori nominal value for the A3 parameter'
             CALL input_rea_opt('difcor.','a3ng_nominal',a3ng_nominal,.FALSE.,found,comment,iun_log)
             dyn%apriori_nominal(3)=a3ng_nominal
             CALL rdnrea('difcor.','a3ng_std',a3ng_std,.false.,found,fail1,fail)
             IF(.not.found.and.a3ng_apriori) THEN
                WRITE(*,*) '******* difini: abnormal end due to a3ng_std *******'
                WRITE(*,*) 'ERROR: a3ng_apriori is true, but a3ng_std is missing'
                IF(iun_log.NE.0) WRITE(iun_log,*) 'ERROR: a3ng_apriori is true, but a3ng_std is missing'
                STOP
             END IF
             dyn%apriori_std(3)=a3ng_std
          END IF
          ! Outgassing dtdelay
          comment='outgassing delay a-priori'
          CALL input_log_opt('difcor.','dtdelay_apriori',dtdelay_apriori,.FALSE.,found,comment,iun_log)
          dyn%apriori(4)=dtdelay_apriori
          IF(dyn%apriori(4))THEN
             comment='a-priori nominal value for the delay parameter'
             CALL input_rea_opt('difcor.','dtdelay_nominal',dtdelay_nominal,.FALSE.,found,comment,iun_log)
             dyn%apriori_nominal(4)=dtdelay_nominal
             CALL rdnrea('difcor.','dtdelay_std',dtdelay_std,.false.,found,fail1,fail)
             IF(.not.found.and.dtdelay_apriori) THEN
                WRITE(*,*) '******* difini: abnormal end due to dtdelay_std *******'
                WRITE(*,*) 'ERROR: dtdelay_apriori is true, but dtdelay_std is missing'
                IF(iun_log.NE.0) WRITE(iun_log,*) 'ERROR: dtdelay_apriori is true, but dtdelay_std is missing'
                STOP
             END IF
             dyn%apriori_std(4)=dtdelay_std
          END IF
       END IF
    END IF
    ! ==== large arrays allocated ==================
    ! total allocatable size:
    ! nblx*(2*nxinbl_large+1)*4 + 8*nxinbl^2*nblx + 8*nxinbl_large^2*nbl_large
    ALLOCATE (indblk(nxinbl_large,nblx),noblk(nblx),indblk_s(nxinbl_large,nblx))
    ALLOCATE (wblk(nxinbl,nxinbl,nblx),wblk_large(nxinbl_large,nxinbl_large,nbl_large))
  END SUBROUTINE difini

END MODULE least_squares


! =========================================================
! RESCOV: rescaling factor of covariance matrix
! to be used as follows:
!     gamma ---> gamma *rescov**2
!     normatr--> normatr/rescov**2
DOUBLE PRECISION FUNCTION rescov(nsolv,nused,csinor)
  IMPLICIT NONE
  ! input: number of parameters to be solved, nomber of scalar obs,
  INTEGER nsolv,nused
  ! norm of residuals (RMS relative to the given observatory weights)
  DOUBLE PRECISION csinor
  IF(nsolv.lt.nused)THEN
     IF(csinor.gt.1.d0)THEN
        rescov=csinor*sqrt(float(nused)/float(nused-nsolv))
     ELSE
        rescov=sqrt(float(nused)/float(nused-nsolv))
     ENDIF
  ELSE
     rescov=1.d0
  ENDIF
END FUNCTION rescov
! =========================================================
! MEANTI time in the middle of the observed arc
! =====================================================
DOUBLE PRECISION FUNCTION meanti(tau,rmsa,rmsd,del,m)
  IMPLICIT NONE
  INTEGER m
  DOUBLE PRECISION tau(m),rmsa(m),rmsd(m),del(m)
  INTEGER i
  DOUBLE PRECISION tw,w
  meanti=0.d0
  IF(m.le.0)THEN
     WRITE(*,*)'meanti: no data '
     STOP
  ENDIF
  tw=0.d0
  DO i=1,m
     w=1.d0/((rmsa(i)*COS(del(i)))**2+rmsd(i)**2)
     meanti=meanti+tau(i)*w
     tw=tw+w
  ENDDO
  meanti=meanti/tw
  RETURN
END FUNCTION meanti
! =============================
! BIZARRE
! logical function deciding when it is a good idea to give up
! differential corrections, both for difcor and difvin.
! uses parameters set in bizset
LOGICAL FUNCTION bizarre(el,ecc)
  USE orbit_elements
  USE output_control
  USE fund_const
  USE dyn_param
  IMPLICIT NONE
  ! input elements
  TYPE(orbit_elem), INTENT(IN) :: el
  DOUBLE PRECISION, INTENT(OUT) :: ecc
  INTEGER i
  ! controls
  DOUBLE PRECISION ecclim0, samin0,samax0,phmin0,ahmax0,qmax0
  COMMON / bizcon /ecclim0, samin0,samax0,phmin0,ahmax0,qmax0
  ! local variables
  DOUBLE PRECISION a, eps, vsize, prscal, r, gm, alpha, q, qg, enne
  DOUBLE PRECISION eq(6),x(3),y(3),ang(3),vlenz(3)
  TYPE(orbit_elem) elco ! for coordinate change
  INTEGER fail_flag
  ! ==========================
  eps=100*epsilon(1.d0)
  bizarre=.false.
  IF(el%coo.eq.'EQU')THEN
     eq=el%coord
     ecc=sqrt(eq(2)**2+eq(3)**2)
     a=eq(1)
     IF(ecc.ge.ecclim0.or.a.ge.samax0.or.a.le.samin0)bizarre=.true.
     IF(a*(1.d0-ecc).le.phmin0.or.a*(1.d0+ecc).ge.ahmax0)bizarre=.true.
     IF(a*(1.d0-ecc).ge.qmax0)bizarre=.true.
  ELSEIF(el%coo.eq.'KEP')THEN
     eq=el%coord
     ecc=eq(2)
     a=eq(1)
     IF(ecc.ge.ecclim0.or.a.ge.samax0.or.a.le.samin0)bizarre=.true.
     IF(a*(1.d0-ecc).le.phmin0.or.a*(1.d0+ecc).ge.ahmax0)bizarre=.true.
     IF(a*(1.d0-ecc).ge.qmax0)bizarre=.true.
  ELSEIF(el%coo.eq.'COM'.or.el%coo.eq.'COT')THEN
     eq=el%coord
     ecc=eq(2)
     IF(eq(1).ge.qmax0) bizarre=.true.
     IF(ecc.ge.1.d0-eps)THEN
        If(ecc.ge.ecclim0) bizarre=.true.
        IF(eq(1).le.phmin0) bizarre=.true.
     ELSE
        a=eq(1)/(1.d0-ecc)
        IF(ecc.ge.ecclim0.or.a.ge.samax0.or.a.le.samin0)bizarre=.true.
        IF(a*(1.d0-ecc).le.phmin0.or.a*(1.d0+ecc).ge.ahmax0)bizarre=.true.
     ENDIF
  ELSE
     IF(el%coo.eq.'CAR')THEN
        r=vsize(el%coord(1:3))
        IF(r.ge.1000.d0)THEN
           bizarre=.true.
           ecc=1.d0
           !           WRITE(ierrou,*)' bizarre: range too large=',r
           !           numerr=numerr+1
           RETURN
        ENDIF
     ELSEIF(el%coo.eq.'ATT')THEN
        IF(el%coord(5).ge.1000.d0)THEN
           bizarre=.true.
           ecc=1.d0
           RETURN
        ENDIF
     ENDIF
     CALL ecc_peri(el,ecc,q,qg,enne)
     IF(ecc.ge.1.d0-eps)THEN
        If(ecc.ge.ecclim0) bizarre=.true.
        IF(q.le.phmin0) bizarre=.true.
     ELSE
        a=q/(1.d0-ecc)
        IF(ecc.ge.ecclim0.or.a.ge.samax0.or.a.le.samin0)bizarre=.true.
        IF(q.le.phmin0.or.qg.ge.ahmax0)bizarre=.true.
     ENDIF
  ENDIF
  IF(rhs.eq.2.and.nls.gt.0)THEN
     DO i=1,2
        IF(dyn%dp(i).gt.dpmax(i).or.dyn%dp(i).lt.dpmin(i))THEN
           bizarre=.true.
        ENDIF
     ENDDO
  ENDIF

END FUNCTION bizarre
!=====================================================================!
! BIZSET                                                              !
! Assign values of bizarre orbit control parameters. If all the input !
! are zero, then default value are assigned                           !
!=====================================================================!
SUBROUTINE bizset(ecclim,samin,samax,phmin,ahmax,qmax)
  USE fund_const
  IMPLICIT NONE
  !=========================================================================================================
  DOUBLE PRECISION, INTENT(IN) :: ecclim ! Limit for the eccentricity
  DOUBLE PRECISION, INTENT(IN) :: samin  ! Minimum semimajor axis
  DOUBLE PRECISION, INTENT(IN) :: samax  ! Maximum semimajor axus
  DOUBLE PRECISION, INTENT(IN) :: phmin  ! Minimum perihelion
  DOUBLE PRECISION, INTENT(IN) :: ahmax  ! Maximum aphelion
  DOUBLE PRECISION, INTENT(IN) :: qmax   ! Maximum perihelion
  !=========================================================================================================
  ! Hidden output, available only to bizarre
  DOUBLE PRECISION :: ecclim0, samin0,samax0,phmin0,ahmax0,qmax0
  COMMON / bizcon /ecclim0, samin0,samax0,phmin0,ahmax0,qmax0
  !=========================================================================================================
  IF(ecclim.EQ.0.d0)THEN
     ecclim0=0.99d0  ! Exclude only almost hyperbolic
  ELSE
     ecclim0=ecclim
  ENDIF
  IF(samin.EQ.0.d0)THEN
     ! Exclude inner Vulcanoids
     IF(rhs.EQ.1)THEN
        samin0=0.1d0
     ELSEIF(rhs.EQ.2)THEN
        samin0=4.3d-5
     ELSEIF(rhs.EQ.3)THEN
        STOP '*** bizset: ERROR! rhs=3 ***'
     ENDIF
  ELSE
     samin0=samin
  ENDIF
  IF(samax.EQ.0.d0)THEN
     IF(rhs.EQ.1)THEN
        samax0=200.d0 ! Exclude long period comets
     ELSEIF(rhs.EQ.2)THEN
        samax0=0.02d0 ! Exclude interplanetary
     ENDIF
  ELSE
     samax0=samax
  ENDIF
  IF(phmin.EQ.0.d0)THEN
     IF(rhs.EQ.1)THEN
        phmin0=0.001d0  ! Exclude Sun dipping comets
     ELSEIF(rhs.EQ.2)THEN
        phmin0=4.24d-5  ! Exclude reentry
     ENDIF
  ELSE
     phmin0=phmin
  ENDIF
  IF(ahmax.eq.0.d0)THEN
     IF(rhs.eq.1)THEN
        ahmax0=400.d0  ! Exclude long periodic comets
     ELSEIF(rhs.eq.2)THEN
        ahmax0=0.03d0  ! Exclude interplanetary launches
     ENDIF
  ELSE
     ahmax0=ahmax
  ENDIF
  IF(qmax.EQ.0.d0)THEN
     IF(rhs.EQ.1)THEN
        qmax0=1000.d0  ! Exclude objects too far to be seen
     ELSEIF(rhs.EQ.2)THEN
        qmax0=0.01d0  ! Exclude out of sphere of influence
     ENDIF
  ELSE
     qmax0=qmax
  ENDIF
END SUBROUTINE bizset
! =====================================================================
! OUTCOV
! =====================================================================
!  output of covariance and normal matrix, computation of eigenvalues
!   input: iun   = output unit
!          nd    = no. parameters solved
!          icor  = flags >0 solved for varaiable
!          gamma = covariance matrix
!          c     = normal matrix
! WARNING: if not all the variables are solved for, the corresponding
!          rows and columns of both gamma and c are set to zero.
!          These null lines are removed in the output.
SUBROUTINE outcov(iun,nd,icor,gamma,c)
  USE dyn_param
  implicit none
  ! output unit, number of parameters
  integer iun,nd
  ! which variables have been solved for
  integer icor(nd)
  ! covariance and normal matrices
  double precision gamma(ndimx,ndimx),c(ndimx,ndimx)
  ! end interface
  ! loop indexes, error flag
  integer i,j,ii,jj,ierr,ncor
  ! covariance and normal matrices
  double precision gam(ndimx*ndimx),cc(ndimx*ndimx),eigv(ndimx*ndimx)
  ! eigenvalues, eigenvectors
  double precision eigval(ndimx),fv1(ndimx),fv2(ndimx)
  ! packing of matrices: no. lines and columns
  ncor=0
  do  i=1,nd
     if(icor(i).gt.0)ncor=ncor+1
  enddo
  if(ncor.eq.0)then
     write(*,*)' no correction, icor=',icor
     return
  endif
  write(iun,*)
  write(iun,*)'no. solve for=',ncor
  write(iun,*)icor
  ! packing of matrices:
  ii=0
  do 10 i=1,nd
     jj=0
     if(icor(i).gt.0)then
        ii=ii+1
        do 11 j=1,nd
           if(icor(j).gt.0)then
              jj=jj+1
              gam(ii+(jj-1)*ncor)=gamma(i,j)
              cc(ii+(jj-1)*ncor)=c(i,j)
           endif
11      enddo
     endif
10 enddo
  ! output covariance
  write(iun,*)
  write(iun,*) 'COVARIANCE MATRIX'
  do 2 j=1,ncor
     write(iun,109) (gam(i+(j-1)*ncor),i=1,ncor)
109  format(6e24.16)
2 enddo
  ! eigenvalues
  call rs(ncor,ncor,gam,eigval,1,eigv,fv1,fv2,ierr)
  write(iun,*)
  write(iun,*) 'EIGENVALUES '
  write(iun,109) (eigval(i),i=1,ncor)
  write(iun,*)
  write(iun,*) 'EIGENVECTORS'
  do 3 j=1,ncor
     write(iun,109) (eigv(i+(j-1)*ncor),i=1,ncor)
3 enddo
  ! normal matrix
  write(iun,*)
  write(iun,*) 'NORMAL MATRIX'
  do 4 j=1,ncor
     write(iun,109) (cc(i+(j-1)*ncor),i=1,ncor)
4 enddo
  write(iun,*)
end subroutine outcov

!
! =====================================================================
! WEAK_DIR
! =====================================================================
!  weak direction and sigma
!   input:       gamma = covariance matrix
!                iun8  = output unit (if positive)
!   output:      wdir  = weak direction vector
!                sdir  = sigma
! ============ INTERFACE====================
SUBROUTINE weak_dir(gamma,wdir,sdir,iun8,coo,coord,units,nd)
  USE fund_const, ONLY: dpig, pig, gk, rhs,gms,aukm
  USE least_squares
  USE dyn_param
  IMPLICIT NONE
  DOUBLE PRECISION, DIMENSION(nd,nd), INTENT(IN) :: gamma
  DOUBLE PRECISION, DIMENSION(nd), INTENT(OUT)  :: wdir
  DOUBLE PRECISION, INTENT(OUT)                :: sdir
  INTEGER, INTENT(IN) :: iun8
  CHARACTER(LEN=3), INTENT(IN)    :: coo
  DOUBLE PRECISION, DIMENSION(6), INTENT(IN) :: coord  ! orbital elements
  DOUBLE PRECISION, DIMENSION(nd),INTENT(OUT) :: units  ! for scaling
  INTEGER, INTENT(IN) :: nd  ! no solve for parameters
  ! ========END INTERFACE====================
  INTEGER :: ierr ! error flag
  INTEGER :: i,j ! loop index
  ! rescaled covariance, scales
  DOUBLE PRECISION, DIMENSION(ndimx,ndimx) ::  gammas
  DOUBLE PRECISION, DIMENSION(ndimx)   ::  scales
  ! eigenvalues, eigenvectors, workspace
  DOUBLE PRECISION, DIMENSION(ndimx,ndimx) :: eigvec
  DOUBLE PRECISION, DIMENSION(ndimx)   :: eigval,fv1,fv2
  DOUBLE PRECISION prscal ! function
  ! =========================================
  ! scaling control
  IF(scaling_lov)THEN
     CALL scale_coef!(coo,coord,scales)
     DO i=1,nd
        DO j=1,nd
           gammas(i,j)=gamma(i,j)*(scales(i)*scales(j))
        ENDDO
     ENDDO
     units(1:nd)=1.d0/scales(1:nd)
  ELSE
     gammas(1:nd,1:nd)=gamma(1:nd,1:nd)
     units=1.d0
  ENDIF
  ! eigenvalues
  call rs(nd,nd,gammas(1:nd,1:nd),eigval(1:nd),1,eigvec(1:nd,1:nd),fv1,fv2,ierr)
  IF(second_lov)THEN
     wdir=eigvec(1:nd,nd-1)
     sdir=sqrt(eigval(nd-1))
  ELSE
     wdir=eigvec(1:nd,nd)
     sdir=sqrt(eigval(nd))
  ENDIF
  ! axial vector continued as true vector
  IF(coo.eq.'EQU'.or.coo.eq.'KEP'.or.coo.eq.'COM'.or.coo.eq.'COT')THEN
     ! if coo='EQU', 'KEP' a growing with sigma, if 'COM' q growing with sigma
     IF(wdir(1).lt.0.d0) wdir=-wdir
  ELSEIF(coo.eq.'ATT')THEN
     ! if coo='ATT' r growing with sigma
     IF(wdir(5).lt.0.d0) wdir=-wdir
  ELSEIF(coo.eq.'CAR')THEN
     ! if coo='CAR' r growing with sigma
     IF(prscal(wdir(1:3),coord(1:3)).lt.0.d0) wdir=-wdir
  ENDIF

  IF(iun8.gt.0)THEN
     WRITE(iun8,*)
     CALL tee(iun8,'WEAK DIRECTION =')
     WRITE(iun8,109) (wdir(i),i=1,nd)
109  FORMAT(6e24.16)
     WRITE(*,110) (wdir(i),i=1,nd)
110  FORMAT(1p,6e12.4)
     WRITE(iun8,*)
     WRITE(iun8,111) sdir
     WRITE(*,111) sdir
111  FORMAT(' WEAK SIGMA ',1p,e12.4)
     WRITE(iun8,*)
  ENDIF


CONTAINS
  ! ==================INTERFACE======================================
  SUBROUTINE scale_coef !(coo,coord,scales)
    USE fund_const
    !  USE least_squares
    !  IMPLICIT NONE
    !  CHARACTER(LEN=3), INTENT(IN)    :: coo
    !  DOUBLE PRECISION, DIMENSION(6), INTENT(IN) :: coord
    !  DOUBLE PRECISION, DIMENSION(6), INTENT(OUT) ::  scales
    !  DOUBLE PRECISION  units(6),
    !==================END INTERFACE====================================
    DOUBLE PRECISION :: r,v, newgk
    DOUBLE PRECISION :: phi ! solar constant
    INTEGER :: fail_flag
    DOUBLE PRECISION :: angmom(3),lenz(3),ecc,p,enne,vsize,d
    IF(rhs.eq.1)THEN
       newgk=gk
    ELSE
       WRITE(*,*)'scale_coeff: not supported rhs=', rhs
       STOP
    ENDIF
    IF(coo.eq.'EQU')THEN
       units(1)=coord(1)   ! a:semimajor axis
       units(2:5)=1
       units(6)=dpig
    ELSEIF(coo.eq.'CAR')THEN
       r=sqrt(coord(1)**2+coord(2)**2+coord(3)**2)
       v=sqrt(coord(4)**2+coord(5)**2+coord(6)**2)
       units(1:3)=r
       units(4:6)=v
    ELSEIF(coo.eq.'KEP')THEN
       !       units(1)=coord(1)    ! a:semimajor axis
       units(1)=1.d-3
       units(2)=1.d0
       units(3)=pig
       units(4:6)=dpig
    ELSEIF(coo.eq.'COM'.or.coo.eq.'COT')THEN
       units(1)=coord(1)  ! q:perihelion distance
       units(2)=1.d0
       units(3)=pig
       units(4)=dpig
       units(5)=dpig
       IF(coo.eq.'COM')THEN
          ! T: period T=2pi*sqrt(a^3/gm)....a=q/(1-e) SINGULAR for e=1
          !     units(6)=dpig*sqrt(((coord(1)**3)/(1-coord(2))**3)/gms)
          units(6)=dpig/newgk*sqrt(coord(1)**3)/sqrt(1+coord(2)) !WARNING changed sign
       ELSEIF(coo.eq.'COT')THEN
          units(6)=dpig
       ENDIF
    ELSEIF(coo.eq.'ATT')THEN
       units(1)=dpig
       units(2)=pig
       units(3)=newgk
       units(4)=newgk
       units(5)=1.d0
       units(6)=newgk
    ELSE
       WRITE(*,*)' scale_coef: coordinates ', coo, ' not known '
       STOP
    ENDIF
    ! nongrav parameters, if any
    IF(nd.gt.6)THEN
       DO j=1,nd-6
          IF(dyn%nmod.eq.1.and.ls(j).eq.1)THEN
             ! radiation pressure ???
             !          units(6+j)=abs(dyn%dp(ls(j)))
             units(6+j)=1.d0
          ELSEIF(dyn%nmod.eq.1.and.ls(j).eq.2)THEN
             IF(coo.EQ.'ATT') THEN
                WRITE(*,*) 'scale_coeff: yark with eletype ATT not possible'
                STOP
             ENDIF
             ! WARNING: assumes unit is 10^-10 AU/d^2
             units(6+j)=1.d0
          ELSE
             ! not knowing what....
             units(6+j)=1.d0
          ENDIF
       ENDDO
    ENDIF
    scales(1:nd)=1.d0/units(1:nd) ! dY/dX has the inverse of the scale units
  END SUBROUTINE scale_coef

END SUBROUTINE weak_dir

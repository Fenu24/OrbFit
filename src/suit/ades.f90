! ----------------------------------------------------------------- !
! MODULE ades                                                       !
! Defines ADES data types and contains subroutines that read and    !
! write ADES data format as .xml and .psv files.                    !
! ----------------------------------------------------------------- !
! ------------------------   WARNING:   --------------------------- !
! All subroutines and functions dealing with reading/writing data   !
! from/to .xml files are provisionally neglected (the corresponding !
! parts of code are commented): being developed using modules       !
! from an external library ( FoX - Fortran XML library:             !
! https://github.com/andreww/fox )                                  !
! ----------------------------------------------------------------- !
! Author: Alessia Bertolucci  - SpaceDyS - 2020                     !
! ----------------------------------------------------------------- !
! Modules used:                                                     !
! fund_const                                                        !
! option_input                                                      !
! file_oper                                                         !
! char_str                                                          !
! output_control                                                    !
!                                                                   !
! Functions contained:                                              !
! get_prec                                                          !
! isreal                                                            !
! find_obs_ades                                                     !
! changed_obs_ades                                                  !
!                                                                   !
! Subroutines contained:                                            !
! read_ades_psv                                                     !
! read_ades_xml                                                     !
! read_xml_elements                                                 !
! read_rec_ades                                                     !
! check_fields                                                      !
! write_xml                                                         !
! write_ades                                                        !
! write_header_ades                                                 !
! write_segment_ades                                                !
! clean_residuals_group                                             !
! clean_remarks                                                     !
! get_residuals_group                                               !
! catcod_mpc                                                        !
! techn_mpc                                                         !
! rot_IAU_to_ECLJ2000                                               !
! geoIAU_to_cart                                                    !
! RADec_accuracy                                                    !
! distPA_to_RADec_cov_mat                                           !
! precision_group_control                                           !
! fill_psv_rec                                                      !
! addobs_psv_rec                                                    !
! sort_ades                                                         !
! update_ades_mpc                                                   !
! update_ades_rwo                                                   !
! ----------------------------------------------------------------- !

MODULE ades

  USE fund_const
  USE option_input
  USE file_oper
  USE char_str
  USE output_control
  !USE m_dom_parse
  !USE m_dom_dom
  !USE m_dom_utils
  !USE FoX_dom

  IMPLICIT NONE
  PRIVATE

  ! Define ADES data type
  ! Variables taken from ADES_Description_v13Jul2018 and ades_master_13Jul.pdf (13 Jul 2018)
  TYPE ades_observ

     ! Identification group (we do not consider artSat)
     CHARACTER(LEN=20)    :: permID            ! IAU permanent designation
     CHARACTER(LEN=30)    :: provID            ! MPC provisional designation (unpacked form) for unnumbered object
     CHARACTER(LEN=8)     :: trkSub            ! Observer-assigned tracklet identifier

     ! Other elements not belonging to any particular group
     CHARACTER(LEN=19)    :: obsID             ! Globally unique observation identifier assigned by the MPC
     CHARACTER(LEN=12)    :: trkID             ! Globally unique tracklet identifier assigned by the MPC
     CHARACTER(LEN=3)     :: mode              ! Mode of instrumentation
     CHARACTER(LEN=4)     :: stn               ! MPC observatory code (3-4 alphanum chars)
     CHARACTER(LEN=4)     :: trx               ! Code of transmitting antenna (3-4 alphanum chars)
     CHARACTER(LEN=4)     :: rcv               ! Code of receveing antenna (3-4 alphanum chars)

     ! Location group
     CHARACTER(LEN=7)     :: sys               ! Reference system for roving and satellite optical obs
     INTEGER              :: ctr               ! Origin of reference system, one int
     ! (https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/naif_ids.html#NAIF%20Object%20ID%20numbers)
     REAL(KIND=dkind)     :: pos(3)            ! Position of observer
     REAL(KIND=dkind)     :: posCov(6)         ! Upper triangular part of pos1, pos2, pos3 covariance matrix in same units as pos1-2-3

     ! Other elements not belonging to any particular group
     CHARACTER(LEN=2)     :: prog              ! Program code assigned by the MPC
     CHARACTER(LEN=40)    :: obsTime           ! UTC date and time of observation
     REAL(KIND=dkind)     :: rmsTime           ! Random uncertainty in obsTime in seconds as estimated by the observer

     ! Observations group
     REAL(KIND=dkind)     :: ra                ! RA for optical observations
     REAL(KIND=dkind)     :: dec               ! Dec for optical observations
     REAL(KIND=dkind)     :: raStar            ! For occultations, RA of the occulted star in degrees (used for conversion to ra)
     REAL(KIND=dkind)     :: decStar           ! For occultations, Dec of the occulted star in degrees (used for conversion to dec)
     CHARACTER(LEN=30)    :: obsCenter         ! Origin of offset observation (used for conversion to ra and dec)
     REAL(KIND=dkind)     :: deltaRA           ! Delta(RA)*cos(Delta(Dec)) in arcsec (used for conversion to ra)
     REAL(KIND=dkind)     :: deltaDec          ! Delta(Dec) in arcsec (used for conversion to dec)
     REAL(KIND=dkind)     :: dist              ! Distance in arcsec for offset and occultation obs (used for conversion to ra and dec)
     REAL(KIND=dkind)     :: pa                ! Position angle in degree for offset and occultation obs (used for conversion to ra and dec)
     REAL(KIND=dkind)     :: rmsRA             ! Random component of RA*cos(Dec) uncertainty (1 sigma) in arcsec from astrometric reduction
     REAL(KIND=dkind)     :: rmsDec            ! Random component of Dec uncertainty (1 sigma) in arcsec from astrometric reduction
     REAL(KIND=dkind)     :: rmsDist           ! Random component of distance in arcsec (1 sigma) for offset and occultation obs
     ! (used for conversion to rmsRA and rmsDec)
     REAL(KIND=dkind)     :: rmsPA             ! Random component of position angle in degree (1 sigma) for offset and occultation obs
     ! (used for conversion to rmsRA and rmsDec)
     REAL(KIND=dkind)     :: rmsCorr           ! Correlation between ra and dec or dist and pa
     REAL(KIND=dkind)     :: delay             ! Delay for radar observations in s
     REAL(KIND=dkind)     :: rmsDelay          ! Delay uncertainty (1 sigma) in mus
     REAL(KIND=dkind)     :: doppler           ! Doppler shift for radar observations in Hz
     REAL(KIND=dkind)     :: rmsDoppler        ! Doppler uncertainty (1 sigma) in Hz

     ! Other elements not belonging to any particular group
     CHARACTER(LEN=8)     :: astCat            ! Star catalog (up to 8 chars)

     ! Photometry group
     REAL(KIND=dkind)     :: mag               ! Apparent magnitude
     REAL(KIND=dkind)     :: rmsMag            ! Apparent magnitude uncertainty (1 sigma) in magnitudes
     CHARACTER(LEN=3)     :: band              ! Band of mag (up to 3 chars)
     CHARACTER(LEN=8)     :: photCat           ! Star catalog used for photometric reduction
     REAL(KIND=dkind)     :: photAp            ! Photometric aperture radius in arcsec
     INTEGER              :: nucMag            ! Nuclear magnitude flag for comets (0 = FALSE, 1 = TRUE)

     ! Other elements not belonging to any particular group
     REAL(KIND=dkind)     :: logSNR            ! log10 of signal-to-noise ratio
     REAL(KIND=dkind)     :: seeing            ! Size of seeing disc in arcsec
     REAL(KIND=dkind)     :: exp               ! Exposure time in seconds
     REAL(KIND=dkind)     :: rmsFit            ! RMS of fit of astrometric comparison stars in arcsec
     INTEGER              :: nStars            ! Number of reference stars in astrometric fit
     INTEGER              :: com               ! Flag to indicate whether the observation is reduced to the center of mass (0=FALSE, 1=TRUE)
     REAL(KIND=dkind)     :: frq               ! Carrier reference frequency in MHz
     CHARACTER(LEN=5)     :: ref               ! Standard reference field used for citations
     CHARACTER(LEN=1)     :: disc              ! Discovery flag, two allowed values: "*" = new discovery record,
     ! "+" = first measurement of a previously observed object
     CHARACTER(LEN=20)    :: subFrm            ! Originally reported reference frame for angular measurements
     CHARACTER(LEN=4)     :: subFmt            ! Format in which the observation was originally submitted to the MPC

     ! Precision group
     INTEGER              :: precTime          ! Precision of time in 1.d-6 days
     REAL(KIND=dkind)     :: precRA            ! Precision of RA in sec
     REAL(KIND=dkind)     :: precDec           ! Precision of Dec in arcsec

     ! Other elements not belonging to any particular group
     REAL(KIND=dkind)     :: uncTime           ! Estimated (presumed systematic) time uncertainty in seconds
     CHARACTER(LEN=6)     :: notes             ! Set of one char note flags (up to 6 chars). first letter
     CHARACTER(LEN=300)   :: remarks           ! Comments provided by the observer

     ! Residuals group
     CHARACTER(LEN=40)    :: orbProd           ! Orbit producer
     CHARACTER(LEN=40)    :: orbID             ! Local orbit ID
     REAL(KIND=dkind)     :: resRA             ! Residuals of RA*cos(Dec) in arcsec
     REAL(KIND=dkind)     :: resDec            ! Residuals of Dec in arcsec
     CHARACTER(LEN=1)     :: selAst            ! Selection flag for optical astrometry (1 char: "A/a", "D/d")
     REAL(KIND=dkind)     :: sigRA             ! Adopted uncertainty of RA*cos(Dec) (1 sigma) in arcsec
     REAL(KIND=dkind)     :: sigDec            ! Adopted uncertainty of Dec (1 sigma) in arcsec
     REAL(KIND=dkind)     :: sigCorr           ! Adopted uncertainty between RA*cos(Dec) and Dec
     REAL(KIND=dkind)     :: sigTime           ! Adopted time uncertainty (1 sigma) in second
     REAL(KIND=dkind)     :: biasRA            ! Adopted bias of RA*cos(Dec) in arcsec
     REAL(KIND=dkind)     :: biasDec           ! Adopted bias of Dec in arcsec
     REAL(KIND=dkind)     :: biasTime          ! Adopted time bias in seconds
     CHARACTER(LEN=40)    :: photProd          ! Producer of photometric residuals (if not present then = orbProd)
     REAL(KIND=dkind)     :: resMag            ! Photometric residuals in magnitudes
     CHARACTER(LEN=1)     :: selPhot           ! Selection flag for photometry (1 char: "A/a", "D/d")
     REAL(KIND=dkind)     :: sigMag            ! Adopted uncertainty of mag (1 sigma) in magnitudes
     REAL(KIND=dkind)     :: biasMag           ! Adopted photometric bias in mag
     CHARACTER(LEN=8)     :: photMod           ! Description of the photometric model used in obtaining the photometric residuals
     ! (up to 8 chars)
     REAL(KIND=dkind)     :: resDelay          ! Residuals of delay in mus
     REAL(KIND=dkind)     :: resDoppler        ! Residuals of Doppler shift in Hz
     CHARACTER(LEN=1)     :: selDelay          ! Selection flag for radar astrometry (1 char: "A/a", "D/d")
     CHARACTER(LEN=1)     :: selDoppler        ! Selection flag for radar astrometry (1 char: "A/a", "D/d")
     REAL(KIND=dkind)     :: sigDelay          ! Adopted uncertainty of delay in mus
     REAL(KIND=dkind)     :: sigDoppler        ! Adopted uncertainty of doppler in Hz

     ! Other elements not belonging to any particular group
     CHARACTER(LEN=1)     :: deprecated        ! Char ("X" only allowed value) to flag the observation preserved for
     ! historical purposes as deprecated not to be used for orbit fitting
     CHARACTER(LEN=40)    :: localUse          ! Only for xml format. There are no restrictions on what can be included under this element
  END TYPE ades_observ

  ! Define ADES flags. The variables are the same as in ades_observ DT (+ 4 variables, see end of DT)
  ! If a variable is TRUE, it means it has been read from the ADES input file
  TYPE ades_flags

     ! Identification group (we do not consider artSat)
     LOGICAL              :: permID
     LOGICAL              :: provID
     LOGICAL              :: trkSub

     ! Other elements not belonging to any particular group
     LOGICAL              :: obsID
     LOGICAL              :: trkID
     LOGICAL              :: mode
     LOGICAL              :: stn
     LOGICAL              :: trx
     LOGICAL              :: rcv

     ! Location group
     LOGICAL              :: sys
     LOGICAL              :: ctr
     LOGICAL              :: pos(3)
     LOGICAL              :: posCov(6)

     ! Other elements not belonging to any particular group
     LOGICAL              :: prog
     LOGICAL              :: obsTime
     LOGICAL              :: rmsTime

     ! Observations group
     LOGICAL              :: ra
     LOGICAL              :: dec
     LOGICAL              :: raStar
     LOGICAL              :: decStar
     LOGICAL              :: obsCenter
     LOGICAL              :: deltaRA
     LOGICAL              :: deltaDec
     LOGICAL              :: dist
     LOGICAL              :: pa
     LOGICAL              :: rmsRA
     LOGICAL              :: rmsDec
     LOGICAL              :: rmsDist
     LOGICAL              :: rmsPA
     LOGICAL              :: rmsCorr
     LOGICAL              :: delay
     LOGICAL              :: rmsDelay
     LOGICAL              :: doppler
     LOGICAL              :: rmsDoppler

     ! Other elements not belonging to any particular group
     LOGICAL              :: astCat

     ! Photometry group
     LOGICAL              :: mag
     LOGICAL              :: rmsMag
     LOGICAL              :: band
     LOGICAL              :: photCat
     LOGICAL              :: photAp
     LOGICAL              :: nucMag

     ! Other elements not belonging to any particular group
     LOGICAL              :: logSNR
     LOGICAL              :: seeing
     LOGICAL              :: exp
     LOGICAL              :: rmsFit
     LOGICAL              :: nStars
     LOGICAL              :: com
     LOGICAL              :: frq
     LOGICAL              :: ref
     LOGICAL              :: disc
     LOGICAL              :: subFrm
     LOGICAL              :: subFmt

     ! Precision group
     LOGICAL              :: precTime
     LOGICAL              :: precRA
     LOGICAL              :: precDec

     ! Other elements not belonging to any particular group
     LOGICAL              :: uncTime
     LOGICAL              :: notes
     LOGICAL              :: remarks

     ! Residuals group
     LOGICAL              :: orbProd
     LOGICAL              :: orbID
     LOGICAL              :: resRA
     LOGICAL              :: resDec
     LOGICAL              :: selAst
     LOGICAL              :: sigTime
     LOGICAL              :: sigRA
     LOGICAL              :: sigDec
     LOGICAL              :: sigCorr
     LOGICAL              :: biasRA
     LOGICAL              :: biasDec
     LOGICAL              :: biasTime
     LOGICAL              :: photProd
     LOGICAL              :: resMag
     LOGICAL              :: selPhot
     LOGICAL              :: sigMag
     LOGICAL              :: biasMag
     LOGICAL              :: photMod
     LOGICAL              :: resDelay
     LOGICAL              :: resDoppler
     LOGICAL              :: selDelay
     LOGICAL              :: selDoppler
     LOGICAL              :: sigDelay
     LOGICAL              :: sigDoppler

     ! Other elements not belonging to any particular group
     LOGICAL              :: deprecated
     LOGICAL              :: localUse

     ! four additional flags to distinguish between optical, offset, occultation and radar obs
     LOGICAL              :: opt
     LOGICAL              :: off
     LOGICAL              :: occ
     LOGICAL              :: rad

  END TYPE ades_flags

  ! Define an additional data type that stores the precision (number of digits after dot)
  ! of different ADES input entries. This will be used for computing different quantities, such as
  ! weights, but also to define the format of the residuals group to be written to the ADES output file
  ! The computation of these precisions is rather rough, but it is necessary to ensure that biases and
  ! weights read in from the output ADES file have enough decimal numbers in order to be consistent
  ! with the their computation within the code in case they are not taken directly from the file
  TYPE ades_precision

     INTEGER     :: rmsTime          ! converted to prec%sigTime
     INTEGER     :: ra               ! always used for resRA, biasRA
     INTEGER     :: dec              ! always used for resDec, biasDec
     INTEGER     :: raStar           ! converted to prec%ra
     INTEGER     :: decStar          ! converted to prec%dec
     INTEGER     :: deltaRA          ! converted to prec%ra
     INTEGER     :: deltaDec         ! converted to prec%dec
     INTEGER     :: dist             ! converted to prec%ra and prec%dec
     INTEGER     :: pa               ! converted to prec%ra and prec%dec
     INTEGER     :: rmsRA            ! converted to prec%sigRA
     INTEGER     :: rmsDec           ! converted to prec%sigDec
     INTEGER     :: rmsDist          ! converted to prec%sigRA and prec%sigDec
     INTEGER     :: rmsPA            ! converted to prec%sigRA and prec%sigDec
     INTEGER     :: rmsCorr          ! converted to prec%sigCorr
     INTEGER     :: precRA           ! converted to prec%sigRA
     INTEGER     :: precDec          ! converted to prec%sigDec
     INTEGER     :: delay            ! always used for resDelay
     INTEGER     :: rmsDelay         ! converted to prec%sigDelay
     INTEGER     :: doppler          ! always used for resDoppler
     INTEGER     :: rmsDoppler       ! converted to prec%sigDoppler
     INTEGER     :: frq              ! converted to prec%sigDoppler
     INTEGER     :: mag              ! always used for resMag
     INTEGER     :: rmsMag           ! converted to prec%sigMag
     INTEGER     :: sigRA            ! always used for sigRA
     INTEGER     :: sigDec           ! always used for sigDec
     INTEGER     :: sigCorr          ! always used for sigCorr
     INTEGER     :: sigDelay         ! always used for sigDelay
     INTEGER     :: sigDoppler       ! always used for sigDoppler
     INTEGER     :: sigMag           ! always used for sigMag

  END TYPE ades_precision

  ! Define an additional data type containing auxiliary variables with extra information about the fit needed for output
  TYPE ades_fit_var

     CHARACTER(LEN=1)    :: ast_flag            ! flag of astrometric fit for each observation
     CHARACTER(LEN=1)    :: radar_flag          ! flag of astrometric fit for each observation
     CHARACTER(LEN=1)    :: mag_flag            ! flag of photometric fit for each observation
     LOGICAL             :: force_w_flags(2)    ! flags for forcing weights of coordinates
     LOGICAL             :: ades_resc_def       ! astrometric residuals flag
     LOGICAL             :: ades_resm_def       ! photometric residuals flag
     REAL(KIND=dkind)    :: chi                 ! chi value of fit

  END TYPE ades_fit_var

  ! Number of ades fields
  ! psv case
  INTEGER, PARAMETER  ::  nfields_ades = 92

  ! Array of all ADES field names in psv format
  CHARACTER(LEN=15), PARAMETER :: field_names(nfields_ades) = &
       & (/"type          ","permID        ","provID        ","artSat        ", &
       &   "trkSub        ","obsID         ","trkID         ","mode          ", &
       &   "stn           ","trx           ","rcv           ","sys           ", &
       &   "ctr           ","pos1          ","pos2          ","pos3          ", &
       &   "posCov11      ","posCov12      ","posCov13      ","posCov22      ", &
       &   "posCov23      ","posCov33      ","prog          ","obsTime       ", &
       &   "rmsTime       ","ra            ","dec           ","raStar        ", &
       &   "decStar       ","obsCenter     ","deltaRA       ","deltaDec      ", &
       &   "dist          ","pa            ","rmsRA         ","rmsDec        ", &
       &   "rmsDist       ","rmsPA         ","rmsCorr       ","delay         ", &
       &   "rmsDelay      ","doppler       ","rmsDoppler    ","astCat        ", &
       &   "mag           ","rmsMag        ","band          ","photCat       ", &
       &   "photAp        ","nucMag        ","logSNR        ","seeing        ", &
       &   "rmsFit        ","nStars        ","com           ","frq           ", &
       &   "exp           ","ref           ","disc          ","subFrm        ", &
       &   "subFmt        ","precTime      ","precRA        ","precDec       ", &
       &   "uncTime       ","notes         ","remarks       ","orbProd       ", &
       &   "orbID         ","resRA         ","resDec        ","selAst        ", &
       &   "sigRA         ","sigDec        ","sigCorr       ","sigTime       ", &
       &   "biasRA        ","biasDec       ","biasTime      ","photProd      ", &
       &   "selPhot       ","sigMag        ","biasMag       ","photMod       ", &
       &   "resMag        ","resDelay      ","selDelay      ","sigDelay      ", &
       &   "resDoppler    ","selDoppler    ","sigDoppler    ","deprecated    "/)

  ! Array of all ADES field names in xml format
  CHARACTER(LEN=15), PARAMETER :: field_names_xml(nfields_ades) = &
       & (/"permID        ","provID        ","artSat        ","trkSub        ", &
       &   "obsID         ","trkID         ","mode          ","stn           ", &
       &   "trx           ","rcv           ","sys           ","ctr           ", &
       &   "pos1          ","pos2          ","pos3          ","posCov11      ", &
       &   "posCov12      ","posCov13      ","posCov22      ","posCov23      ", &
       &   "posCov33      ","prog          ","obsTime       ","rmsTime       ", &
       &   "ra            ","dec           ","raStar        ","decStar       ", &
       &   "obsCenter     ","deltaRA       ","deltaDec      ","dist          ", &
       &   "pa            ","rmsRA         ","rmsDec        ","rmsDist       ", &
       &   "rmsPA         ","rmsCorr       ","delay         ","rmsDelay      ", &
       &   "doppler       ","rmsDoppler    ","astCat        ","mag           ", &
       &   "rmsMag        ","band          ","photCat       ","photAp        ", &
       &   "nucMag        ","logSNR        ","seeing        ","exp           ", &
       &   "rmsFit        ","nStars        ","com           ","frq           ", &
       &   "ref           ","disc          ","subFrm        ","subFmt        ", &
       &   "precTime      ","precRA        ","precDec       ","uncTime       ", &
       &   "notes         ","remarks       ","orbProd       ","orbID         ", &
       &   "resRA         ","resDec        ","selAst        ","sigRA         ", &
       &   "sigDec        ","sigCorr       ","sigTime       ","biasRA        ", &
       &   "biasDec       ","biasTime      ","photProd      ","resMag        ", &
       &   "selPhot       ","sigMag        ","biasMag       ","photMod       ", &
       &   "resDelay      ","selDelay      ","sigDelay      ","resDoppler    ", &
       &   "selDoppler    ","sigDoppler    ","deprecated    ","localUse      "/)

  ! undefined ADES obs data type
  ! organized (new lines) as the blocks in the data type
  TYPE(ades_observ), PARAMETER :: undefined_ades_observ = &
       & ADES_OBSERV('','','',                                                                                     &
       & '','','','','','',                                                                                        &
       & '',0,(/0.d0,0.d0,0.d0/),(/0.d0,0.d0,0.d0,0.d0,0.d0,0.d0/),                                                &
       & '','',0.d0,                                                                                               &
       & 0.d0,0.d0,0.d0,0.d0,'',0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,                  &
       & '',                                                                                                       &
       & 0.d0,0.d0,'','',0.d0,0,                                                                                   &
       & 0.d0,0.d0,0.d0,0.d0,0,0,0.d0,'','','','',                                                                 &
       & 0,0.d0,0.d0,                                                                                              &
       & 0.d0,'','',                                                                                               &
       & '','',0.d0,0.d0,'',0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,'',0.d0,'',0.d0,0.d0,'',0.d0,0.d0,'','',0.d0,0.d0,  &
       & '','')

  ! undefined ADES flag data type
  ! organized (new lines) as the blocks in the data type
  TYPE(ades_flags), PARAMETER :: undefined_ades_flags = &
       & ADES_FLAGS(.FALSE.,.FALSE.,.FALSE.,                                                                       &
       & .FALSE.,.FALSE.,.FALSE.,.FALSE.,.FALSE.,.FALSE.,                                                          &
       & .FALSE.,.FALSE.,(/.FALSE.,.FALSE.,.FALSE./),(/.FALSE.,.FALSE.,.FALSE.,.FALSE.,.FALSE.,.FALSE./),          &
       & .FALSE.,.FALSE.,.FALSE.,                                                                                  &
       & .FALSE.,.FALSE.,.FALSE.,.FALSE.,.FALSE.,.FALSE.,.FALSE.,.FALSE.,.FALSE.,.FALSE.,.FALSE.,.FALSE.,.FALSE.,  & ! continues 
       & .FALSE.,.FALSE.,.FALSE.,.FALSE.,.FALSE.,                                                                  & ! end of block
       & .FALSE.,                                                                                                  &
       & .FALSE.,.FALSE.,.FALSE.,.FALSE.,.FALSE.,.FALSE.,                                                          &
       & .FALSE.,.FALSE.,.FALSE.,.FALSE.,.FALSE.,.FALSE.,.FALSE.,.FALSE.,.FALSE.,.FALSE.,.FALSE.,                  &
       & .FALSE.,.FALSE.,.FALSE.,                                                                                  &
       & .FALSE.,.FALSE.,.FALSE.,                                                                                  &
       & .FALSE.,.FALSE.,.FALSE.,.FALSE.,.FALSE.,.FALSE.,.FALSE.,.FALSE.,.FALSE.,.FALSE.,.FALSE.,.FALSE.,.FALSE.,  & ! continues 
       & .FALSE.,.FALSE.,.FALSE.,.FALSE.,.FALSE.,.FALSE.,.FALSE.,.FALSE.,.FALSE.,.FALSE.,.FALSE.,                  & ! end of block
       & .FALSE.,                                                                                                  &
       & .FALSE.,.FALSE.,.FALSE.,.FALSE.,.FALSE.)                                                                    ! 4 new variables

  ! undefined ADES precision data type
  TYPE(ades_precision), PARAMETER :: undefined_ades_precision = &
       & ADES_PRECISION(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)

  ! undefined ADES fit_var data type
  TYPE(ades_fit_var), PARAMETER :: undefined_ades_fit_var = &
       & ADES_FIT_VAR('','','',(/.FALSE.,.FALSE./),.FALSE.,.FALSE.,0.d0)

  INCLUDE 'parobx.h90' ! observation numbers: maximum

  ! public ades_obs and ades_flg arrays, along with the total number of observations
  CHARACTER(LEN=200)    :: ades_filename             ! ADES filename (.psv or .xml)
  TYPE(ades_observ)     :: ades_obs(nobx)            ! ADES observations
  TYPE(ades_flags)      :: ades_flg(nobx)            ! ADES keyword reading flags
  TYPE(ades_precision)  :: ades_prec(nobx)           ! ADES entry precisions
  TYPE(ades_fit_var)    :: ades_fit(nobx)            ! ADES fit info
  INTEGER               :: nobs_ades                 ! number of selected ADES observations
  INTEGER               :: ntot_ades                 ! total number of lines (including headers and key records) of ADES file
  CHARACTER(LEN=200)    :: ades_header(nobx)         ! ADES input file headers
  CHARACTER(LEN=600)    :: ades_key_rec(nobx)        ! ADES input file key records
  CHARACTER(LEN=600)    :: ades_data_rec(nobx)       ! ADES input file data records
  INTEGER               :: ntot_to_nobs_idx(nobx)    ! mapping array between total and observation indexes (ades_data_rec -> ades_obs/flg)
  INTEGER               :: nobs_to_ntot_idx(nobx)    ! mapping array between observation and total indexes (ades_obs/flg -> ades_data_rec)
  INTEGER               :: header_length             ! number of ObsContext lines 
  CHARACTER(LEN=10)     :: ades_error_model          ! error model used for ADES weights and bias computation
  LOGICAL               :: new_ades                  ! false if .psv input file exists and data from it have not been changed
  CHARACTER(LEN=200)    :: psv_name                  ! ades file name without file extension
  INTEGER               :: len_ades_name             ! length of psv_name
  LOGICAL               :: new_key_rec(nobx)         ! TRUE if an obs type is different wrt to the previous one and thus need a new keys line

  ! allowed values for different ADES keys (see also http://www.minorplanetcenter.net/iau/info/ADESFieldValues.html)
  ! observatory codes are checked during conversion
  CHARACTER(LEN=3), PARAMETER   :: mode_values(9) = (/"CCD","VID","PHO","ENC","PMT","MIC","MER","TDI","UNK"/)
  CHARACTER(LEN=5), PARAMETER   :: sys_rov_values(3) = (/"WGS84","ITRF ","IAU  "/)
  CHARACTER(LEN=7), PARAMETER   :: sys_sat_values(2) = (/"ICRF_AU","ICRF_KM"/)
  INTEGER,          PARAMETER   :: ctr_values(13) = (/0,1,2,3,4,5,6,7,8,9,10,301,399/)   ! ctr only planets + sun + moon (SPICE codes)
  CHARACTER(LEN=8), PARAMETER   :: astCat_values(47) = (/"Gaia2   ","Gaia1   ","URAT1   ","UCAC4   ","PPMXL   ","NOMAD   ", &   
       & "2MASS   ","UCAC3   ","UCAC2   ","UCAC1   ","USNOB1  ","USNOA2  ","USNOSA2 ","USNOA1  ","USNOSA1 ","Tyc2    ","Tyc1    ", &    
       & "Hip2    ","Hip1    ","ACT     ","GSCACT  ","GSC2.3  ","GSC2.2  ","GSC1.2  ","GSC1.1  ","GSC1.0  ","GSC     ","SDSS8   ", &
       & "SDSS7   ","CMC15   ","CMC14   ","SSTRC4  ","SSTRC1  ","MPOSC3  ","PPM     ","AC      ","SAO1984 ","SAO     ","AGK3    ", &
       & "FK4     ","ACRS    ","LickGas ","Ida93   ","Perth70 ","COSMOS  ","Yale    ","UNK     "/)
  CHARACTER(LEN=8), PARAMETER   :: photCat_values(47) = astCat_values                                
  CHARACTER(LEN=3), PARAMETER   :: band_values(38) = (/"Vj ","Rc ","Ic ","Bj ","Uj ","Sg ","Sr ","Si ","Sz ","Pg ","Pr ","Pi ", &
       & "Pz ","Pw ","Ao ","Ac ","B  ","V  ","H  ","R  ","I  ","C  ","U  ","J  ","g  ","r  ","y  ","c  ","o  ","i  ","z  ","w  ", &
       & "G  ","L  ","K  ","Y  ","u  ","   "/)
  INTEGER, PARAMETER            :: precTime_values(6) = (/1,10,100,1000,10000,100000/)
  REAL(KIND=dkind), PARAMETER   :: precRADec_values(7) = (/1.d-3,1.d-2,1.d-1,1.d0,6.d-1,6.d0,6.d1/)
  CHARACTER(LEN=1), PARAMETER   :: notes_values(40) = (/"A","a","B","b","c","C","D","d","E","F","f","G","g","H","h","I","i","K", &
       & "k","M","m","N","n","O","o","P","p","Q","R","r","S","s","T","t","U","u","V","W","w","X"/)
  CHARACTER(LEN=6), PARAMETER   :: photMod_values(7) = (/"HG    ","H     ","HG1G2 ","HG12  ","NEATM ","FRM   ","      "/)

  ! public types and variables
  PUBLIC :: ades_filename, ades_observ, ades_flags, ades_precision, ades_fit_var, ades_obs, ades_flg, ades_prec, ades_fit,  &
       &    nobs_ades, ades_header, ades_key_rec, ades_data_rec, ntot_to_nobs_idx, nobs_to_ntot_idx, ntot_ades, ades_error_model, &
       &    new_ades, psv_name, new_key_rec, header_length, len_ades_name

  ! public parameters
  PUBLIC :: undefined_ades_flags, undefined_ades_observ, undefined_ades_precision, undefined_ades_fit_var, notes_values

  ! public routines
  PUBLIC :: read_ades_psv, write_ades, get_astCat, catcod_mpc, techn_mpc, rot_IAU_to_ECLJ2000, &
       &       geoIAU_to_cart, RADec_accuracy, distPA_to_RADec_cov_mat,  precision_group_control, &
       &       fill_psv_rec, addobs_psv_rec, sort_ades, update_ades_mpc, update_ades_rwo ! ,read_ades_xml

  ! public functions
  PUBLIC :: get_prec, isreal

CONTAINS

  ! READING ADES

  ! --------------------------------------------- !
  ! SUBROUTINE read_ades_psv                      !
  ! It reads observations from ADES .psv file     !
  ! --------------------------------------------- !

  SUBROUTINE read_ades_psv(filename,ades_obs_out,ades_flg_out,ades_prec_out,ades_fit_out, &
       & nobs_ades_out,ntot_ades_out,ades_header_out,ades_key_rec_out,ades_data_rec_out,ntot_to_nobs_out, &
       & nobs_to_ntot_out,header_length_out,ades_error_model_out,new_key_rec_out)

    CHARACTER(LEN=*),     INTENT(IN)  :: filename                  ! name of input ADES file (either .obs_psv or .rwo_psv)

    TYPE(ades_observ),    INTENT(OUT) :: ades_obs_out(nobx)        ! ADES observations
    TYPE(ades_flags),     INTENT(OUT) :: ades_flg_out(nobx)        ! ADES keyword reading flags
    TYPE(ades_precision), INTENT(OUT) :: ades_prec_out(nobx)       ! ADES entry precisions
    TYPE(ades_fit_var),   INTENT(OUT) :: ades_fit_out(nobx)        ! ADES fit info
    INTEGER,              INTENT(OUT) :: nobs_ades_out             ! number of selected ADES observations
    INTEGER,              INTENT(OUT) :: ntot_ades_out             ! total number of lines (including headers and key records) of ADES file
    CHARACTER(LEN=200),   INTENT(OUT) :: ades_header_out(nobx)     ! ADES input file headers
    CHARACTER(LEN=600),   INTENT(OUT) :: ades_key_rec_out(nobx)    ! ADES input file key records
    CHARACTER(LEN=600),   INTENT(OUT) :: ades_data_rec_out(nobx)   ! ADES input file data records
    INTEGER,              INTENT(OUT) :: ntot_to_nobs_out(nobx)    ! mapping array between total and observation indexes 
    INTEGER,              INTENT(OUT) :: nobs_to_ntot_out(nobx)    ! mapping array between observation and total indexes 
    INTEGER,              INTENT(OUT) :: header_length_out         ! number of ObsContext lines 
    CHARACTER(LEN=10),    INTENT(OUT) :: ades_error_model_out      ! error model used for ADES weights and bias computation
    LOGICAL,              INTENT(OUT) :: new_key_rec_out(nobx)

    INTEGER               :: unit_ades           ! ADES file unit
    TYPE(ades_observ)     :: ades_loc_obs        ! local observation from ADES file
    TYPE(ades_flags)      :: ades_loc_flg        ! local reading flags from ADES file
    TYPE(ades_precision)  :: ades_loc_prec       ! local entry precisions form ADES file
    CHARACTER(LEN=200)    :: ades_loc_header     ! local header from ADES file
    CHARACTER(LEN=500)    :: ades_loc_key_rec    ! local key record from ADES file
    CHARACTER(LEN=500)    :: ades_loc_data_rec   ! local data record from ADES file
    TYPE(ades_fit_var)    :: ades_loc_fit        ! local fit info from ADES file
    CHARACTER(LEN=30)     :: keys(nfields_ades)  ! array of keys read from line
    LOGICAL               :: skipped_obs_loc     ! true if observation has been skipped because of format errors

    INTEGER           :: iline               ! i-th record
    LOGICAL           :: value_found         ! flag used to store data records
    INTEGER           :: old_count_keys      ! stored count_keys of previous key record
    LOGICAL           :: eof                 ! end of file flag
    INTEGER           :: lf,idx_start,idx_end
    LOGICAL           :: error
    CHARACTER(LEN=65) :: error_msg

    INTEGER, EXTERNAL :: lench
    EXTERNAL          :: filopn,filclo

    ! initialization of local variables 
    iline = 0
    idx_start = 0
    idx_end = 0
    nobs_ades_out = 0
    ntot_ades_out = 0
    ades_obs_out = undefined_ades_observ
    ades_flg_out = undefined_ades_flags
    ades_prec_out = undefined_ades_precision
    ades_fit_out = undefined_ades_fit_var
    ades_header_out = ""
    ades_key_rec_out = ""
    ades_data_rec_out = ""
    header_length_out = 0
    ntot_to_nobs_out = 0
    nobs_to_ntot_out = 0
    ades_error_model_out = ""
    new_key_rec_out = .FALSE.

    lf=lench(filename)

    ! open ADES file
    CALL filopn(unit_ades,filename,'OLD')

    ! loop over file lines
    DO

       ! number of line
       iline=iline+1

       ! read one record from ADES .psv file
       CALL read_rec_ades(unit_ades,iline,keys,filename,ades_loc_obs,ades_loc_flg,ades_loc_prec,ades_loc_fit, &
            & ades_loc_header,ades_loc_key_rec,ades_loc_data_rec,value_found,old_count_keys,eof,error,skipped_obs_loc)

       IF(eof)THEN
          EXIT
       ELSEIF(error)THEN 
          error_msg = 'ERROR: read_ades_psv: input error in .psv file: '
          WRITE(*,*) error_msg//filename(1:lf)
          WRITE(ierrou,*) error_msg//filename(1:lf)
          numerr=numerr+1
          RETURN
       ENDIF

       ! read error model, if present, in the header
       IF (INDEX(ades_loc_header,"Errormodel") .GT. 0) THEN
          idx_end = LEN(ades_loc_header)
          idx_start = INDEX(ades_loc_header,"=") + 1
          ades_error_model_out = trim(adjustl(ades_loc_header(idx_start:idx_end)))
       ENDIF

       IF (skipped_obs_loc) THEN
          IF (ntot_ades_out .GT. 0) ades_key_rec_out(ntot_ades-1) = ""  ! delete also previous line of keys
       ELSE
          ntot_ades_out = ntot_ades_out + 1
          ! check ntot_ades upper limit
          IF (ntot_ades_out .GT. nobx) THEN
             error_msg = 'ERROR: read_ades_psv: too many lines in input .psv file: '
             WRITE(*,*) error_msg//filename(1:lf)
             WRITE(ierrou,*) error_msg//filename(1:lf)
             error = .TRUE.
             numerr=numerr+1
             RETURN
          END IF

          ! store headers
          ades_header_out(ntot_ades_out) = ades_loc_header
          IF (ades_loc_header .NE. "") header_length_out = header_length_out + 1
          ! store key records
          ades_key_rec_out(ntot_ades_out) = ades_loc_key_rec
          ! store data records
          ades_data_rec_out(ntot_ades_out) = ades_loc_data_rec
          ! store ades_obs/flg index corresponding to ntot_ades index
          IF (value_found) THEN
             ntot_to_nobs_out(ntot_ades_out) = nobs_ades_out + 1
             nobs_to_ntot_out(nobs_ades_out+1) = ntot_ades_out
          ELSE
             ntot_to_nobs_out(ntot_ades_out) = 0    ! if it is a line of keys, the next is a data line related to obs nobs_ades+1 
             IF (ades_loc_key_rec .NE. "") new_key_rec_out(nobs_ades_out+1) = .TRUE.   ! new_key_rec depends on obs index
          ENDIF
          ! store observations and related flags
          IF (value_found) THEN
             nobs_ades_out = nobs_ades_out + 1
             ! fill obs and flags
             ades_obs_out(nobs_ades_out) = ades_loc_obs
             ades_flg_out(nobs_ades_out) = ades_loc_flg
             ades_prec_out(nobs_ades_out) = ades_loc_prec
             ades_fit_out(nobs_ades_out) = ades_loc_fit
          ENDIF
       ENDIF

    END DO

    IF (nobs_ades_out == 0) THEN
       error_msg = 'ERROR: read_ades_psv: no valid observation in input ADES file: '
       WRITE(*,*) error_msg//filename(1:lf)
       WRITE(ierrou,*) error_msg//filename(1:lf)
       numerr=numerr+1
       RETURN
    ENDIF

    ! close ADES file
    CALL filclo(unit_ades,' ')

  END SUBROUTINE read_ades_psv

  ! --------------------------------------------- !
  ! SUBROUTINE read_ades_xml                      !
  ! It reads an ADES xml file, storing data in    !
  ! ades public variables. Uses functions from    !
  ! FoX library.                                  !
  ! --------------------------------------------- !
  ! Lines of code related to FoX types, functions !
  ! and subroutines are commented                 !

  SUBROUTINE read_ades_xml(file)

    CHARACTER(LEN=*), INTENT(IN) :: file        ! name of input ADES file

    !  TYPE(Node),     POINTER   :: myDoc,np
    !  TYPE(NodeList), POINTER   :: adesList
    INTEGER                   :: line,obs_n,idx_start,idx_end,i,lf
    LOGICAL                   :: error_flg,val_found
    CHARACTER(LEN=40)         :: aux_str
    CHARACTER(LEN=90)         :: error_msg

    INTEGER,         EXTERNAL :: lench

    ! set public variables
    ades_filename = file
    nobs_ades = 0
    ntot_ades = 0
    header_length = 0
    ades_obs = undefined_ades_observ
    ades_flg = undefined_ades_flags
    ades_prec = undefined_ades_precision
    ades_fit = undefined_ades_fit_var
    ades_header = ""
    ades_key_rec = ""
    ades_data_rec = ""
    ades_error_model = ""
    ntot_to_nobs_idx = 0
    nobs_to_ntot_idx = 0

    lf=lench(file)

    ! initialize variables
    idx_start = 0
    idx_end = 0
    ades_header(1) = "# version=2017"    ! intialize header
    line = 1                             ! update header line number
    obs_n = 0
    error_flg = .FALSE.
    error_msg = ""

    ! initialize variables for read_xml_elements
    val_found = .FALSE.

    ! Load in the document
    !  myDoc => parseFile(file)

    ! Find the ades element
    !  adesList => getElementsByTagName(myDoc, "ades")

    !  np => item(adesList, 0)   
    !  CALL read_xml_elements(np,obs_n,line,error_flg,val_found)

    IF (error_flg) THEN      ! if error is found, writes a message and returns
       error_msg = 'ERROR: read_ades_xml: input error in ADES file'
       WRITE(*,*) trim(error_msg)," ",trim(adjustl(file(1:lf))),' at obs ',obs_n
       WRITE(ierrou,*) trim(error_msg)," ",trim(adjustl(file(1:lf))),' at obs ',obs_n
       numerr=numerr+1
       RETURN                
    ENDIF

    ! check obs number upper and lower limits
    IF (obs_n .GT. nobx) THEN
       error_msg = 'ERROR: read_ades_xml: too many observations in input ADES file: '
       WRITE(*,*) trim(error_msg)," ",trim(adjustl(file(1:lf)))
       WRITE(ierrou,*) trim(error_msg)," ",trim(adjustl(file(1:lf)))
       error_flg = .TRUE.
       numerr=numerr+1
    END IF
    IF (obs_n .LE. 1 .AND. .NOT. val_found) THEN
       error_msg = 'ERROR: read_ades_xml: no valid observations in input ADES file: '
       WRITE(*,*) trim(error_msg)," ",trim(adjustl(file(1:lf)))
       WRITE(ierrou,*) trim(error_msg)," ",trim(adjustl(file(1:lf)))
       error_flg = .TRUE.
       numerr=numerr+1
    ENDIF
    IF (error_flg) RETURN

    ! set public variables
    nobs_ades = obs_n
    header_length = line          ! store total number of header lines in .psv format

    ! read error model, if present, in the header
    DO i = 1,line
       IF (INDEX(ades_header(i),"Errormodel") .GT. 0) THEN
          aux_str = ades_header(i)
          idx_end = LEN(aux_str)
          idx_start = INDEX(aux_str,"=") + 1
          ades_error_model = trim(adjustl(aux_str(idx_start:idx_end)))
          EXIT
       ENDIF
    ENDDO

    CALL fill_psv_rec(.FALSE.,.FALSE.)  ! we need to store public variables and strings in psv format 

    ! Clear up all allocated memory 
    !  CALL destroy(myDoc)

  END SUBROUTINE read_ades_xml

  ! ---------------------------------------- !
  ! FUNCTION has_subelements(p)              !
  ! For .xml files, it uses FoX library.     !
  ! It is .true. if p has any ELEMENT_NODE   !
  ! children and .false. if not.             !
  ! ---------------------------------------- !
  ! LOGICAL FUNCTION has_subelements(p)
  !    TYPE(Node), POINTER :: p
  !    TYPE(NodeList), POINTER :: children
  !    TYPE(Node), POINTER :: cp
  !    INTEGER :: i
  !    has_subelements = .false.
  !    children => getChildNodes(p)
  !    DO i = 0, getLength(children) - 1
  !       cp => item(children, i)
  !       IF (getNodeTYPE(cp) == ELEMENT_NODE) THEN
  !          has_subelements=.true.
  !          RETURN
  !       ENDIF
  !    ENDDO
  ! END FUNCTION has_subelements

  ! ------------------------------------------------- !
  ! SUBROUTINE read_xml_elements                      !
  ! Reads xml elements walking recursively through a  !
  ! node and storing all the elements as it goes.     !
  ! If a node has no sub-elements, it stores          !
  ! the content as well as the tag name.              !
  ! ------------------------------------------------- !
  ! Lines of code related to FoX types, functions and !
  ! subroutines are commented                         !

  RECURSIVE SUBROUTINE read_xml_elements(obs_idx,line,error,value_found)

    !TYPE(Node), INTENT(IN), POINTER   :: p
    INTEGER,    INTENT(INOUT)         :: obs_idx                      ! observation index
    INTEGER,    INTENT(INOUT)         :: line                         ! obsContext line index
    LOGICAL,    INTENT(INOUT)         :: error                        ! error
    LOGICAL,    INTENT(INOUT)         :: value_found                  ! TRUE if we find and store data value

    LOGICAL             :: key_found                                  ! TRUE if we recognize key value
    CHARACTER(LEN=90)   :: warn_msg,error_msg                         ! warning and error messages
    CHARACTER(LEN=17)   :: rtn                                        ! this routine
    CHARACTER(LEN=20)   :: int_err                                    ! auxiliary string
    INTEGER             :: len_ie,len_chi,len_str,idx_chi,idx_AS      ! auxiliary length and position integers
    CHARACTER(LEN=40)   :: aux_str,time_str
    CHARACTER(LEN=15)   :: chi_str
    CHARACTER(LEN=1)    :: flg_str
    INTEGER             :: i,j,lf
    CHARACTER(LEN=30)   :: key,branch   
    CHARACTER(LEN=40)   :: ades_entry

    !TYPE(NodeList), POINTER :: children
    !TYPE(Node),     POINTER :: cp,np

    INTEGER,      EXTERNAL  :: lench
    LOGICAL,      EXTERNAL  :: isnum

    LOGICAL, SAVE :: data_flag            ! TRUE if we are in obsData block
    LOGICAL, SAVE :: warn_numb,warn_str   ! warning flags

    DATA  data_flag,warn_numb,warn_str/.FALSE.,.FALSE.,.FALSE./

    rtn = 'read_xml_elements'
    chi_str = ''
    flg_str = ''
    key_found = .FALSE.
    lf = lench(ades_filename)
    !np => null()

    !IF (has_subelements(p)) THEN     ! If the node has sub-elements it does not have ades data
    ! Sub-elements only in ObsContext block and at opening of ObsData
    !children => getChildNodes(p)
    !branch = getLocalName(p)

    !IF (branch == 'obsData') then
    data_flag = .TRUE.
    !ELSEIF (branch == 'optical' .OR. branch == 'offset' .OR. branch == 'occultation' .OR. branch == 'radar') THEN  ! there is a new observation block
    data_flag = .TRUE.  ! assignes value again, in case there is not obsData tag
    value_found = .TRUE.
    IF (obs_idx .GE. 1) THEN  ! control on observation just stored: if some inconsistency has been found, write warning message and skip obs
       IF (warn_numb .OR. warn_str) THEN
          WRITE(int_err,*) obs_idx
          CALL rmsp(int_err,len_ie)
          IF (warn_numb) THEN
             warn_msg = trim(adjustl(key))//' value not a number &
                  &in input ADES file: '//trim(ades_filename(1:lf))
          ELSE IF (warn_str) THEN
             warn_msg = trim(adjustl(key))//' string too long &
                  &in input ADES file: '//trim(ades_filename(1:lf))
          ENDIF
          ! write warning message
          WRITE(*,*) 'WARNING: ',rtn,': ',warn_msg
          WRITE(*,*) 'Skipping observation: '//int_err(1:len_ie)
          ! set value found to FALSE to skip element storage
          value_found = .FALSE.
       ENDIF
       ! check consistency of fields
       !CALL check_fields(filename,ades_obs(obs_idx),ades_flg(obs_idx),value_found,error,OBS_IDX=obs_idx)

       ! reset warning flags to false for next observation
       warn_numb = .FALSE.
       warn_str = .FALSE.

    ENDIF

    IF (value_found) obs_idx = obs_idx+1   ! update obs_idx at current number of observation
    ! if NOT value_found, obs_idx is not updated: obs must be skipped thus ades_obs(obs_idx) is overwritten

    ! store flag about type of current observation
    !IF (branch == 'optical') THEN
    ades_flg(obs_idx)%opt = .TRUE.
    !ELSEIF (branch == 'offset') THEN
    ades_flg(obs_idx)%off = .TRUE.
    !ELSEIF (branch == 'occultation') THEN
    ades_flg(obs_idx)%occ = .TRUE.
    !ELSEIF (branch == 'radar') THEN
    ades_flg(obs_idx)%rad = .TRUE.
    !ENDIF

    !ELSEIF (branch == 'obsContext') THEN
    data_flag = .FALSE.

    !ENDIF
    !DO i = 0, getLength(children) - 1
    ! IF (error) EXIT   ! if an error is found in previous observation, exit now from iteration
    !cp => item(children, i)
    !IF (branch == 'obsContext' .AND. has_subelements(cp)) THEN
    line = line + 1       ! is the corresponding line number as it would be in psv file
    !key = getLocalName(cp)
    ades_header(line) = "# "//trim(adjustl(key))     ! there is the need to store public ades_header with psv format
    !ENDIF
    !IF (getNodeTYPE(cp) == ELEMENT_NODE) THEN   
    !CALL read_xml_elements(cp,obs_idx,line,error,value_found)
    !ENDIF
    !ENDDO
    !ELSE
    ! nodes with no sub-elements have data
    ! key = getLocalName(p)    
    ! ades_entry = getTextContent(p)
    IF (.NOT. data_flag) THEN
       line = line + 1
       aux_str = trim(adjustl(ades_entry))
       len_str = len(trim(adjustl(ades_entry)))
       ades_header(line) = "! "//trim(adjustl(key))//" "//aux_str(1:len_str)
    ELSE
       DO j=1,size(field_names_xml)

          ! compares key entry to elements of field_names_xml array to check if key is recognized
          IF (trim(adjustl(key)) .EQ. trim(adjustl(field_names_xml(j)))) THEN
             key_found = .TRUE.
             EXIT
          ENDIF
       ENDDO
       ! error if key has not been recognized
       IF(.NOT. key_found) THEN
          WRITE(int_err,*) obs_idx
          CALL rmsp(int_err,len_ie)
          error_msg='key entry/ies not recognized &
               &in input ADES file: '//trim(ades_filename(1:lf))//' observation: '//int_err(1:len_ie)
          WRITE(ierrou,*) 'ERROR: ',rtn,': ',error_msg
          error=.TRUE.
       ENDIF
       SELECT CASE (trim(adjustl(key)))
       CASE('permID')
          IF (len(trim(adjustl(ades_entry))) .LE. 20) THEN
             aux_str = ades_entry
             CALL rmsp(aux_str,len_str)
             ades_obs(obs_idx)%permID = aux_str(1:len_str)
             ades_flg(obs_idx)%permID = .TRUE.
          ELSE
             warn_str = .TRUE.

          ENDIF
       CASE('provID')
          IF (len(trim(adjustl(ades_entry))) .LE. 30) THEN
             len_str = len(trim(adjustl(ades_entry)))
             aux_str = trim(adjustl(ades_entry))
             !CALL rmsp(aux_str,len_str)
             ades_obs(obs_idx)%provID = aux_str(1:len_str)
             ades_flg(obs_idx)%provID = .TRUE.
          ELSE
             warn_str = .TRUE.

          ENDIF
       CASE('trkSub')
          IF (len(trim(adjustl(ades_entry))) .LE. 8) THEN
             ades_obs(obs_idx)%trkSub = ades_entry
             ades_flg(obs_idx)%trkSub = .TRUE.
          ELSE
             warn_str = .TRUE.

          ENDIF

          ! other elements not belonging to any particular group
       CASE('obsID')
          IF (len(trim(adjustl(ades_entry))) .LE. 19) THEN
             ades_obs(obs_idx)%obsID = ades_entry
             ades_flg(obs_idx)%obsID = .TRUE.
          ELSE
             warn_str = .TRUE.

          ENDIF
       CASE('trkID')
          IF (len(trim(adjustl(ades_entry))) .LE. 12) THEN
             ades_obs(obs_idx)%trkID = ades_entry
             ades_flg(obs_idx)%trkID = .TRUE.
          ELSE
             warn_str = .TRUE.

          ENDIF
       CASE('mode')
          IF (len(trim(adjustl(ades_entry))) .LE. 3) THEN
             ades_obs(obs_idx)%mode = ades_entry
             ades_flg(obs_idx)%mode = .TRUE.
          ELSE
             warn_str = .TRUE.

          ENDIF
       CASE('stn')
          IF (len(trim(adjustl(ades_entry))) .LE. 4) THEN
             ades_obs(obs_idx)%stn = ades_entry
             ades_flg(obs_idx)%stn = .TRUE.
          ELSE
             warn_str = .TRUE.

          ENDIF
       CASE('trx')
          IF (len(trim(adjustl(ades_entry))) .LE. 4) THEN
             ades_obs(obs_idx)%trx = ades_entry
             ades_flg(obs_idx)%trx = .TRUE.
          ELSE
             warn_str = .TRUE.

          ENDIF
       CASE('rcv')
          IF (len(trim(adjustl(ades_entry))) .LE. 4) THEN
             ades_obs(obs_idx)%rcv = ades_entry
             ades_flg(obs_idx)%rcv = .TRUE.
          ELSE
             warn_str = .TRUE.

          ENDIF

          ! Location group
       CASE('sys')
          IF (len(trim(adjustl(ades_entry))) .LE. 7) THEN
             ades_obs(obs_idx)%sys = ades_entry
             ades_flg(obs_idx)%sys = .TRUE.
          ELSE
             warn_str = .TRUE.

          ENDIF
       CASE('ctr')
          IF (trim(adjustl(ades_entry)) .EQ. "") THEN
             ades_flg(obs_idx)%ctr = .TRUE.
          ELSEIF (isnum(trim(ades_entry))) THEN
             READ(ades_entry,*)ades_obs(obs_idx)%ctr
             ades_flg(obs_idx)%ctr = .TRUE.
          ELSE
             warn_numb = .TRUE.

          ENDIF
       CASE('pos1')
          IF (trim(adjustl(ades_entry)) .EQ. "") THEN
             ades_flg(obs_idx)%pos(1) = .TRUE.
          ELSEIF (isreal(ades_entry)) THEN
             READ(ades_entry,*)ades_obs(obs_idx)%pos(1)
             ades_flg(obs_idx)%pos(1) = .TRUE.
          ELSE
             warn_numb = .TRUE.

          ENDIF
       CASE('pos2')
          IF (trim(adjustl(ades_entry)) .EQ. "") THEN
             ades_flg(obs_idx)%pos(2) = .TRUE.
          ELSEIF (isreal(ades_entry)) THEN
             READ(ades_entry,*)ades_obs(obs_idx)%pos(2)
             ades_flg(obs_idx)%pos(2) = .TRUE.
          ELSE
             warn_numb = .TRUE.

          ENDIF
       CASE('pos3')
          IF (trim(adjustl(ades_entry)) .EQ. "") THEN
             ades_flg(obs_idx)%pos(3) = .TRUE.
          ELSEIF (isreal(ades_entry)) THEN
             READ(ades_entry,*)ades_obs(obs_idx)%pos(3)
             ades_flg(obs_idx)%pos(3) = .TRUE.
          ELSE
             warn_numb = .TRUE.

          ENDIF
       CASE('posCov11')
          IF (trim(adjustl(ades_entry)) .EQ. "") THEN
             ades_flg(obs_idx)%posCov(1) = .TRUE.
          ELSEIF (isreal(ades_entry)) THEN
             READ(ades_entry,*)ades_obs(obs_idx)%posCov(1)
             ades_flg(obs_idx)%posCov(1) = .TRUE.
          ELSE
             warn_numb = .TRUE.

          ENDIF
       CASE('posCov12')
          IF (trim(adjustl(ades_entry)) .EQ. "") THEN
             ades_flg(obs_idx)%posCov(2) = .TRUE.
          ELSEIF (isreal(ades_entry)) THEN
             READ(ades_entry,*)ades_obs(obs_idx)%posCov(2)
             ades_flg(obs_idx)%posCov(2) = .TRUE.
          ELSE
             warn_numb = .TRUE.

          ENDIF
       CASE('posCov13')
          IF (trim(adjustl(ades_entry)) .EQ. "") THEN
             ades_flg(obs_idx)%posCov(3) = .TRUE.
          ELSEIF (isreal(ades_entry)) THEN
             READ(ades_entry,*)ades_obs(obs_idx)%posCov(3)
             ades_flg(obs_idx)%posCov(3) = .TRUE.
          ELSE
             warn_numb = .TRUE.

          ENDIF
       CASE('posCov22')
          IF (trim(adjustl(ades_entry)) .EQ. "") THEN
             ades_flg(obs_idx)%posCov(4) = .TRUE.
          ELSEIF (isreal(ades_entry)) THEN
             READ(ades_entry,*)ades_obs(obs_idx)%posCov(4)
             ades_flg(obs_idx)%posCov(4) = .TRUE.
          ELSE
             warn_numb = .TRUE.

          ENDIF
       CASE('posCov23')
          IF (trim(adjustl(ades_entry)) .EQ. "") THEN
             ades_flg(obs_idx)%posCov(5) = .TRUE.
          ELSEIF (isreal(ades_entry)) THEN
             READ(ades_entry,*)ades_obs(obs_idx)%posCov(5)
             ades_flg(obs_idx)%posCov(5) = .TRUE.
          ELSE
             warn_numb = .TRUE.

          ENDIF
       CASE('posCov33')
          IF (trim(adjustl(ades_entry)) .EQ. "") THEN
             ades_flg(obs_idx)%posCov(6) = .TRUE.
          ELSEIF (isreal(ades_entry)) THEN
             READ(ades_entry,*)ades_obs(obs_idx)%posCov(6)
             ades_flg(obs_idx)%posCov(6) = .TRUE.
          ELSE
             warn_numb = .TRUE.

          ENDIF

          ! other elements not belonging to any particular group
       CASE('prog')
          IF (len(trim(adjustl(ades_entry))) .LE. 2) THEN
             ades_obs(obs_idx)%prog = ades_entry
             ades_flg(obs_idx)%prog = .TRUE.
          ELSE
             warn_str = .TRUE.

          ENDIF
       CASE('obsTime')
          time_str = ades_entry
          CALL rmsp(time_str,len_str)
          ! lower threshold should be len_str = 23, but we relax it to 20
          IF (len_str .GE. 20 .AND. len_str .LE. 40) THEN
             IF (isreal(time_str(1:4)) .AND. isreal(time_str(6:7)) .AND. isreal(time_str(9:10)) &
                  & .AND. isreal(time_str(12:13)) .AND. isreal(time_str(15:16)) .AND. &
                  & isreal(time_str(18:len_str-1))) THEN
                ades_obs(obs_idx)%obsTime = ades_entry
                ades_flg(obs_idx)%obsTime = .TRUE.
             ELSE
                warn_numb = .TRUE.

             ENDIF
          ELSE
             ! we write the warning message here since it may be triggered by a too short,
             ! rather than too long, string
             WRITE(int_err,*) obs_idx
             warn_msg = trim(adjustl(key))//' string too long or too short&
                  &in input ADES file: '//trim(ades_filename(1:lf))//' Skipping observation '//int_err
             WRITE(*,*) 'WARNING: ',rtn,': ',warn_msg
             value_found = .FALSE.
             RETURN
          ENDIF
       CASE('rmsTime')
          IF (trim(adjustl(ades_entry)) .EQ. "") THEN
             ades_flg(obs_idx)%rmsTime = .TRUE.
          ELSEIF (isreal(ades_entry)) THEN
             READ(ades_entry,*)ades_obs(obs_idx)%rmsTime
             ades_flg(obs_idx)%rmsTime = .TRUE.
             ades_prec(obs_idx)%rmsTime = get_prec(ades_entry)
          ELSE
             warn_numb = .TRUE.

          ENDIF
          ! observation group
       CASE('ra')
          IF (trim(adjustl(ades_entry)) .EQ. "") THEN
             ades_flg(obs_idx)%ra = .TRUE.
          ELSEIF (isreal(ades_entry)) THEN
             READ(ades_entry,*)ades_obs(obs_idx)%ra
             ades_flg(obs_idx)%ra = .TRUE.
             ades_prec(obs_idx)%ra = get_prec(ades_entry)
             ! flag observation as optical
             ades_flg(obs_idx)%opt = .TRUE.
          ELSE
             warn_numb = .TRUE.

          ENDIF
       CASE('dec')
          IF (trim(adjustl(ades_entry)) .EQ. "") THEN
             ades_flg(obs_idx)%dec = .TRUE.
          ELSEIF (isreal(ades_entry)) THEN
             READ(ades_entry,*)ades_obs(obs_idx)%dec
             ades_flg(obs_idx)%dec = .TRUE.
             ades_prec(obs_idx)%dec = get_prec(ades_entry)
             ! flag observation as optical
             ades_flg(obs_idx)%opt = .TRUE.
          ELSE
             warn_numb = .TRUE.

          ENDIF
       CASE('raStar')
          IF (trim(adjustl(ades_entry)) .EQ. "") THEN
             ades_flg(obs_idx)%raStar = .TRUE.
          ELSEIF (isreal(ades_entry)) THEN
             READ(ades_entry,*)ades_obs(obs_idx)%raStar
             ades_flg(obs_idx)%raStar = .TRUE.
             ades_prec(obs_idx)%raStar = get_prec(ades_entry)
             ! flag observation as occultation
             ades_flg(obs_idx)%occ = .TRUE.
          ELSE
             warn_numb = .TRUE.

          ENDIF
       CASE('decStar')
          IF (trim(adjustl(ades_entry)) .EQ. "") THEN
             ades_flg(obs_idx)%decStar = .TRUE.
          ELSEIF (isreal(ades_entry)) THEN
             READ(ades_entry,*)ades_obs(obs_idx)%decStar
             ades_flg(obs_idx)%decStar = .TRUE.
             ades_prec(obs_idx)%decStar = get_prec(ades_entry)
             ! flag observation as occultation
             ades_flg(obs_idx)%occ = .TRUE.
          ELSE
             warn_numb = .TRUE.

          ENDIF
       CASE('obsCenter')
          IF (len(trim(adjustl(ades_entry))) .LE. 20) THEN
             ades_obs(obs_idx)%obsCenter = ades_entry
             ades_flg(obs_idx)%obsCenter = .TRUE.
             ! flag observation as offset
             ades_flg(obs_idx)%off = .TRUE.
          ELSE
             warn_str = .TRUE.

          ENDIF
       CASE('deltaRA')
          IF (trim(adjustl(ades_entry)) .EQ. "") THEN
             ades_flg(obs_idx)%deltaRA = .TRUE.
          ELSEIF (isreal(ades_entry)) THEN
             READ(ades_entry,*)ades_obs(obs_idx)%deltaRA
             ades_flg(obs_idx)%deltaRA = .TRUE.
             ades_prec(obs_idx)%deltaRA = get_prec(ades_entry)
          ELSE
             warn_numb = .TRUE.

          ENDIF
       CASE('deltaDec')
          IF (trim(adjustl(ades_entry)) .EQ. "") THEN
             ades_flg(obs_idx)%deltaDec = .TRUE.
          ELSEIF (isreal(ades_entry)) THEN
             READ(ades_entry,*)ades_obs(obs_idx)%deltaDec
             ades_flg(obs_idx)%deltaDec = .TRUE.
             ades_prec(obs_idx)%deltaDec = get_prec(ades_entry)
          ELSE
             warn_numb = .TRUE.

          ENDIF
       CASE('dist')
          IF (trim(adjustl(ades_entry)) .EQ. "") THEN
             ades_flg(obs_idx)%dist = .TRUE.
          ELSEIF (isreal(ades_entry)) THEN
             READ(ades_entry,*)ades_obs(obs_idx)%dist
             ades_flg(obs_idx)%dist = .TRUE.
             ades_prec(obs_idx)%dist = get_prec(ades_entry)
          ELSE
             warn_numb = .TRUE.

          ENDIF
       CASE('pa')
          IF (trim(adjustl(ades_entry)) .EQ. "") THEN
             ades_flg(obs_idx)%pa = .TRUE.
          ELSEIF (isreal(ades_entry)) THEN
             READ(ades_entry,*)ades_obs(obs_idx)%pa
             ades_flg(obs_idx)%pa = .TRUE.
             ades_prec(obs_idx)%pa = get_prec(ades_entry)
          ELSE
             warn_numb = .TRUE.

          ENDIF
       CASE('rmsRA')
          IF (trim(adjustl(ades_entry)) .EQ. "") THEN
             ades_flg(obs_idx)%rmsRA = .TRUE.
          ELSEIF (isreal(ades_entry)) THEN
             READ(ades_entry,*)ades_obs(obs_idx)%rmsRA
             ades_flg(obs_idx)%rmsRA = .TRUE.
             ades_prec(obs_idx)%rmsRA = get_prec(ades_entry)
          ELSE
             warn_numb = .TRUE.

          ENDIF
       CASE('rmsDec')
          IF (trim(adjustl(ades_entry)) .EQ. "") THEN
             ades_flg(obs_idx)%rmsDec = .TRUE.
          ELSEIF (isreal(ades_entry)) THEN
             READ(ades_entry,*)ades_obs(obs_idx)%rmsDec
             ades_flg(obs_idx)%rmsDec = .TRUE.
             ades_prec(obs_idx)%rmsDec = get_prec(ades_entry)
          ELSE
             warn_numb = .TRUE.

          ENDIF
       CASE('rmsDist')
          IF (trim(adjustl(ades_entry)) .EQ. "") THEN
             ades_flg(obs_idx)%rmsDist = .TRUE.
          ELSEIF (isreal(ades_entry)) THEN
             READ(ades_entry,*)ades_obs(obs_idx)%rmsDist
             ades_flg(obs_idx)%rmsDist = .TRUE.
             ades_prec(obs_idx)%rmsDist = get_prec(ades_entry)
          ELSE
             warn_numb = .TRUE.

          ENDIF
       CASE('rmsPA')
          IF (trim(adjustl(ades_entry)) .EQ. "") THEN
             ades_flg(obs_idx)%rmsPA = .TRUE.
          ELSEIF (isreal(ades_entry)) THEN
             READ(ades_entry,*)ades_obs(obs_idx)%rmsPA
             ades_flg(obs_idx)%rmsPA = .TRUE.
             ades_prec(obs_idx)%rmsPA = get_prec(ades_entry)
          ELSE
             warn_numb = .TRUE.

          ENDIF
       CASE('rmsCorr')
          IF (trim(adjustl(ades_entry)) .EQ. "") THEN
             ades_flg(obs_idx)%rmsCorr = .TRUE.
          ELSEIF (isreal(ades_entry)) THEN
             READ(ades_entry,*)ades_obs(obs_idx)%rmsCorr
             ades_flg(obs_idx)%rmsCorr = .TRUE.
             ades_prec(obs_idx)%rmsCorr = get_prec(ades_entry)
          ELSE
             warn_numb = .TRUE.

          ENDIF
       CASE('delay')
          IF (trim(adjustl(ades_entry)) .EQ. "") THEN
             ades_flg(obs_idx)%delay = .TRUE.
          ELSEIF (isreal(ades_entry)) THEN
             READ(ades_entry,*)ades_obs(obs_idx)%delay
             ades_flg(obs_idx)%delay = .TRUE.
             ades_prec(obs_idx)%delay = get_prec(ades_entry)
             ! flag observation as radar
             ades_flg(obs_idx)%rad = .TRUE.
          ELSE
             warn_numb = .TRUE.

          ENDIF
       CASE('rmsDelay')
          IF (trim(adjustl(ades_entry)) .EQ. "") THEN
             ades_flg(obs_idx)%rmsDelay = .TRUE.
          ELSEIF (isreal(ades_entry)) THEN
             READ(ades_entry,*)ades_obs(obs_idx)%rmsDelay
             ades_flg(obs_idx)%rmsDelay = .TRUE.
             ades_prec(obs_idx)%rmsDelay = get_prec(ades_entry)
          ELSE
             warn_numb = .TRUE.

          ENDIF
       CASE('doppler')
          IF (trim(adjustl(ades_entry)) .EQ. "") THEN
             ades_flg(obs_idx)%doppler = .TRUE.
          ELSEIF (isreal(ades_entry)) THEN
             READ(ades_entry,*)ades_obs(obs_idx)%doppler
             ades_flg(obs_idx)%doppler = .TRUE.
             ades_prec(obs_idx)%doppler = get_prec(ades_entry)
             ! flag observation as radar
             ades_flg(obs_idx)%rad = .TRUE.
          ELSE
             warn_numb = .TRUE.

          ENDIF
       CASE('rmsDoppler')
          IF (trim(adjustl(ades_entry)) .EQ. "") THEN
             ades_flg(obs_idx)%rmsDoppler = .TRUE.
          ELSEIF (isreal(ades_entry)) THEN
             READ(ades_entry,*)ades_obs(obs_idx)%rmsDoppler
             ades_flg(obs_idx)%rmsDoppler = .TRUE.
             ades_prec(obs_idx)%rmsDoppler = get_prec(ades_entry)
          ELSE
             warn_numb = .TRUE.

          ENDIF

          ! other elements not belonging to any particular group
       CASE('astCat')
          IF (len(trim(adjustl(ades_entry))) .LE. 8) THEN
             ades_obs(obs_idx)%astCat = ades_entry
             ades_flg(obs_idx)%astCat = .TRUE.
          ELSE
             warn_str = .TRUE.

          ENDIF

          ! photometry group
       CASE('mag')
          IF (trim(adjustl(ades_entry)) .EQ. "") THEN
             ades_flg(obs_idx)%mag = .TRUE.
             ades_obs(obs_idx)%mag = 9.9d9
          ELSEIF (isreal(ades_entry)) THEN
             READ(ades_entry,*)ades_obs(obs_idx)%mag
             ades_flg(obs_idx)%mag = .TRUE.
             ades_prec(obs_idx)%mag = get_prec(ades_entry)
          ELSE
             warn_numb = .TRUE.

          ENDIF
       CASE('rmsMag')
          IF (trim(adjustl(ades_entry)) .EQ. "") THEN
             ades_flg(obs_idx)%rmsMag = .TRUE.
          ELSEIF (isreal(ades_entry)) THEN
             READ(ades_entry,*)ades_obs(obs_idx)%rmsMag
             ades_flg(obs_idx)%rmsMag = .TRUE.
             ades_prec(obs_idx)%rmsMag = get_prec(ades_entry)
          ELSE
             warn_numb = .TRUE.

          ENDIF
       CASE('band')
          IF (len(trim(adjustl(ades_entry))) .LE. 3) THEN
             ades_obs(obs_idx)%band = ades_entry
             ades_flg(obs_idx)%band = .TRUE.
          ELSE
             warn_str = .TRUE.

          ENDIF
       CASE('photCat')
          IF (len(trim(adjustl(ades_entry))) .LE. 8) THEN
             ades_obs(obs_idx)%photCat = ades_entry
             ades_flg(obs_idx)%photCat = .TRUE.
          ELSE
             warn_str = .TRUE.

          ENDIF
       CASE('photAp')
          IF (trim(adjustl(ades_entry)) .EQ. "") THEN
             ades_flg(obs_idx)%photAp = .TRUE.
          ELSEIF (isreal(ades_entry)) THEN
             READ(ades_entry,*)ades_obs(obs_idx)%photAp
             ades_flg(obs_idx)%photAp = .TRUE.
          ELSE
             warn_numb = .TRUE.

          ENDIF
       CASE('nucMag')
          IF (trim(adjustl(ades_entry)) .EQ. "") THEN
             ades_flg(obs_idx)%nucMag = .TRUE.
          ELSEIF (isnum(trim(ades_entry))) THEN
             READ(ades_entry,*)ades_obs(obs_idx)%nucMag
             ades_flg(obs_idx)%nucMag = .TRUE.
          ELSE
             warn_numb = .TRUE.

          ENDIF

          ! other elements not belonging to any particular group
       CASE('logSNR')
          IF (trim(adjustl(ades_entry)) .EQ. "") THEN
             ades_flg(obs_idx)%logSNR = .TRUE.
          ELSEIF (isreal(ades_entry)) THEN
             READ(ades_entry,*)ades_obs(obs_idx)%logSNR
             ades_flg(obs_idx)%logSNR = .TRUE.
          ELSE
             warn_numb = .TRUE.

          ENDIF
       CASE('seeing')
          IF (trim(adjustl(ades_entry)) .EQ. "") THEN
             ades_flg(obs_idx)%seeing = .TRUE.
          ELSEIF (isreal(ades_entry)) THEN
             READ(ades_entry,*)ades_obs(obs_idx)%seeing
             ades_flg(obs_idx)%seeing = .TRUE.
          ELSE
             warn_numb = .TRUE.

          ENDIF
       CASE('exp')
          IF (trim(adjustl(ades_entry)) .EQ. "") THEN
             ades_flg(obs_idx)%exp = .TRUE.
          ELSEIF (isreal(ades_entry)) THEN
             READ(ades_entry,*)ades_obs(obs_idx)%exp
             ades_flg(obs_idx)%exp = .TRUE.
          ELSE
             warn_numb = .TRUE.

          ENDIF
       CASE('rmsFit')
          IF (trim(adjustl(ades_entry)) .EQ. "") THEN
             ades_flg(obs_idx)%rmsFit = .TRUE.
          ELSEIF (isreal(ades_entry)) THEN
             READ(ades_entry,*)ades_obs(obs_idx)%rmsFit
             ades_flg(obs_idx)%rmsFit = .TRUE.
          ELSE
             warn_numb = .TRUE.

          ENDIF
       CASE('nStars')
          IF (trim(adjustl(ades_entry)) .EQ. "") THEN
             ades_flg(obs_idx)%nStars = .TRUE.
          ELSEIF (isnum(trim(ades_entry))) THEN
             READ(ades_entry,*)ades_obs(obs_idx)%nStars
             ades_flg(obs_idx)%nStars = .TRUE.
          ELSE
             warn_numb = .TRUE.

          ENDIF
       CASE('com')
          IF (trim(adjustl(ades_entry)) .EQ. "") THEN
             ades_flg(obs_idx)%com = .TRUE.
          ELSEIF (isnum(trim(ades_entry))) THEN
             READ(ades_entry,*)ades_obs(obs_idx)%com
             ades_flg(obs_idx)%com = .TRUE.
          ELSE
             warn_numb = .TRUE.

          ENDIF
       CASE('frq')
          IF (trim(adjustl(ades_entry)) .EQ. "") THEN
             ades_flg(obs_idx)%frq = .TRUE.
          ELSEIF (isreal(ades_entry)) THEN
             READ(ades_entry,*)ades_obs(obs_idx)%frq
             ades_flg(obs_idx)%frq = .TRUE.
             ades_prec(obs_idx)%frq = get_prec(ades_entry)
          ELSE
             warn_numb = .TRUE.

          ENDIF
       CASE('ref')
          IF (len(trim(adjustl(ades_entry))) .LE. 16) THEN
             ades_obs(obs_idx)%ref = ades_entry
             ades_flg(obs_idx)%ref = .TRUE.
          ELSE
             warn_str = .TRUE.

          ENDIF
       CASE('disc')
          IF (len(trim(adjustl(ades_entry))) .LE. 1) THEN
             ades_obs(obs_idx)%disc = ades_entry
             ades_flg(obs_idx)%disc = .TRUE.
          ELSE
             warn_str = .TRUE.

          ENDIF
       CASE('subFrm')
          IF (len(trim(adjustl(ades_entry))) .LE. 20) THEN
             ades_obs(obs_idx)%subFrm = ades_entry
             ades_flg(obs_idx)%subFrm = .TRUE.
          ELSE
             warn_str = .TRUE.

          ENDIF
       CASE('subFmt')
          IF (len(trim(adjustl(ades_entry))) .LE. 4) THEN
             ades_obs(obs_idx)%subFmt = ades_entry
             ades_flg(obs_idx)%subFmt = .TRUE.
          ELSE
             warn_str = .TRUE.

          ENDIF

          ! precision group
       CASE('precTime')
          IF (trim(adjustl(ades_entry)) .EQ. "") THEN
             ades_flg(obs_idx)%precTime = .TRUE.
          ELSEIF (isnum(trim(ades_entry))) THEN
             READ(ades_entry,*)ades_obs(obs_idx)%precTime
             ades_flg(obs_idx)%precTime = .TRUE.
          ELSE
             warn_numb = .TRUE.

          ENDIF
       CASE('precRA')
          IF (trim(adjustl(ades_entry)) .EQ. "") THEN
             ades_flg(obs_idx)%precRA = .TRUE.
          ELSEIF (isreal(ades_entry)) THEN
             READ(ades_entry,*)ades_obs(obs_idx)%precRA
             ades_flg(obs_idx)%precRA = .TRUE.
             ades_prec(obs_idx)%precRA = get_prec(ades_entry)
          ELSE
             warn_numb = .TRUE.

          ENDIF
       CASE('precDec')
          IF (trim(adjustl(ades_entry)) .EQ. "") THEN
             ades_flg(obs_idx)%precDec = .TRUE.
          ELSEIF (isreal(ades_entry)) THEN
             READ(ades_entry,*)ades_obs(obs_idx)%precDec
             ades_flg(obs_idx)%precDec = .TRUE.
             ades_prec(obs_idx)%precDec = get_prec(ades_entry)
          ELSE
             warn_numb = .TRUE.
          ENDIF

          ! other elements not belonging to any particular group
       CASE('uncTime')
          IF (trim(adjustl(ades_entry)) .EQ. "") THEN
             ades_flg(obs_idx)%uncTime = .TRUE.
          ELSEIF (isreal(ades_entry)) THEN
             READ(ades_entry,*)ades_obs(obs_idx)%uncTime
             ades_flg(obs_idx)%uncTime = .TRUE.
          ELSE
             warn_numb = .TRUE.
          ENDIF
       CASE('notes')
          IF (len(trim(adjustl(ades_entry))) .LE. 6) THEN
             ades_obs(obs_idx)%notes = ades_entry
             ades_flg(obs_idx)%notes = .TRUE.
          ELSE
             warn_str = .TRUE.
          ENDIF

       CASE('remarks')
          IF (len(trim(adjustl(ades_entry))) .LE. 200) THEN
             ades_obs(obs_idx)%remarks = ades_entry
             ades_flg(obs_idx)%remarks = .TRUE.
             idx_chi = INDEX(ades_entry,'chi =')
             idx_AS = INDEX(ades_entry,'A = ') + INDEX(ades_entry,'S = ')
             IF (idx_chi .GT. 0) THEN
                IF (idx_AS .GT. 0) THEN
                   chi_str = ades_entry(idx_chi+5:idx_AS-1)
                   CALL rmsp(chi_str,len_chi)
                   !IF (len_chi .GT. 9) THEN
                   !READ(chi_str,'(1P,E9.2)') chi(obs_idx)
                   !ELSE
                   !READ(chi_str,'(F9.2)') chi(obs_idx)
                   !ENDIF
                ENDIF
             ENDIF
             IF (idx_AS .GT. 0) THEN
                flg_str = ades_entry(idx_AS:idx_AS)
                !IF (flg_str .EQ. "S") radar_flag(obs_idx) = ades_entry(idx_AS+4:idx_AS+4)
                !IF (flg_str .EQ. "A") THEN
                !ast_flag(obs_idx) = ades_entry(idx_AS+4:idx_AS+4)
                !mag_flag(obs_idx) = ades_entry(idx_AS+10:idx_AS+10)
                !ENDIF
             ENDIF
             IF (idx_chi .GT. 0 .AND. idx_AS .LE. 0) THEN
                warn_msg = 'Missing astrometric fit flag in input ADES file: '//trim(ades_filename(1:lf))
                ! write warning message
                WRITE(*,*) 'WARNING: ',rtn,': ',warn_msg
                ! set value found to FALSE to skip record storage
                value_found = .FALSE.
                RETURN
             ENDIF
          ELSE
             warn_str = .TRUE.
          ENDIF

          ! Residuals group
       CASE('orbProd')
          IF (len(trim(adjustl(ades_entry))) .LE. 40) THEN
             ades_obs(obs_idx)%orbProd = ades_entry
             ades_flg(obs_idx)%orbProd = .TRUE.
          ELSE
             warn_str = .TRUE.

          ENDIF
       CASE('orbID')
          IF (len(trim(adjustl(ades_entry))) .LE. 40) THEN
             ades_obs(obs_idx)%orbID = ades_entry
             ades_flg(obs_idx)%orbID = .TRUE.
          ELSE
             warn_str = .TRUE.

          ENDIF
       CASE('resRA')
          IF (trim(adjustl(ades_entry)) .EQ. "") THEN
             ades_flg(obs_idx)%resRA = .TRUE.
             !ades_resc_def(obs_idx) = .FALSE.
          ELSEIF (isreal(ades_entry)) THEN
             READ(ades_entry,*)ades_obs(obs_idx)%resRA
             ades_flg(obs_idx)%resRA = .TRUE.
             !ades_resc_def(obs_idx) = .TRUE.
          ELSE
             warn_numb = .TRUE.

          ENDIF
       CASE('resDec')
          IF (trim(adjustl(ades_entry)) .EQ. "") THEN
             ades_flg(obs_idx)%resDec = .TRUE.
             !ades_resc_def(obs_idx) = .FALSE.
          ELSEIF (isreal(ades_entry)) THEN
             READ(ades_entry,*)ades_obs(obs_idx)%resDec
             ades_flg(obs_idx)%resDec = .TRUE.
          ELSE
             warn_numb = .TRUE.

          ENDIF
       CASE('selAst')
          IF (len(trim(adjustl(ades_entry))) .LE. 1) THEN
             ades_obs(obs_idx)%selAst = ades_entry
             ades_flg(obs_idx)%selAst = .TRUE.
          ELSE
             warn_str = .TRUE.

          ENDIF
       CASE('sigRA')
          IF (trim(adjustl(ades_entry)) .EQ. "") THEN
             ades_flg(obs_idx)%sigRA = .TRUE.
          ELSEIF (isreal(ades_entry)) THEN
             READ(ades_entry,*)ades_obs(obs_idx)%sigRA
             ades_flg(obs_idx)%sigRA = .TRUE.
             ades_prec(obs_idx)%sigRA = get_prec(ades_entry)
          ELSE
             warn_numb = .TRUE.

          ENDIF
       CASE('sigDec')
          IF (trim(adjustl(ades_entry)) .EQ. "") THEN
             ades_flg(obs_idx)%sigDec = .TRUE.
          ELSEIF (isreal(ades_entry)) THEN
             READ(ades_entry,*)ades_obs(obs_idx)%sigDec
             ades_flg(obs_idx)%sigDec = .TRUE.
             ades_prec(obs_idx)%sigDec = get_prec(ades_entry)
          ELSE
             warn_numb = .TRUE.

          ENDIF
       CASE('sigCorr')
          IF (trim(adjustl(ades_entry)) .EQ. "") THEN
             ades_flg(obs_idx)%sigCorr = .TRUE.
          ELSEIF (isreal(ades_entry)) THEN
             READ(ades_entry,*)ades_obs(obs_idx)%sigCorr
             ades_flg(obs_idx)%sigCorr = .TRUE.
             ades_prec(obs_idx)%sigCorr = get_prec(ades_entry)
          ELSE
             warn_numb = .TRUE.

          ENDIF
       CASE('sigTime')
          IF (trim(adjustl(ades_entry)) .EQ. "") THEN
             ades_flg(obs_idx)%sigTime = .TRUE.
          ELSEIF (isreal(ades_entry)) THEN
             READ(ades_entry,*)ades_obs(obs_idx)%sigTime
             ades_flg(obs_idx)%sigTime = .TRUE.
          ELSE
             warn_numb = .TRUE.

          ENDIF
       CASE('biasRA')
          IF (trim(adjustl(ades_entry)) .EQ. "") THEN
             ades_flg(obs_idx)%biasRA = .TRUE.
          ELSEIF (isreal(ades_entry)) THEN
             READ(ades_entry,*)ades_obs(obs_idx)%biasRA
             ades_flg(obs_idx)%biasRA = .TRUE.
          ELSE
             warn_numb = .TRUE.

          ENDIF
       CASE('biasDec')
          IF (trim(adjustl(ades_entry)) .EQ. "") THEN
             ades_flg(obs_idx)%biasDec = .TRUE.
          ELSEIF (isreal(ades_entry)) THEN
             READ(ades_entry,*)ades_obs(obs_idx)%biasDec
             ades_flg(obs_idx)%biasDec = .TRUE.
          ELSE
             warn_numb = .TRUE.

          ENDIF
       CASE('biasTime')
          IF (trim(adjustl(ades_entry)) .EQ. "") THEN
             ades_flg(obs_idx)%biasTime = .TRUE.
          ELSEIF (isreal(ades_entry)) THEN
             READ(ades_entry,*)ades_obs(obs_idx)%biasTime
             ades_flg(obs_idx)%biasTime = .TRUE.
          ELSE
             warn_numb = .TRUE.

          ENDIF
       CASE('photProd')
          IF (len(trim(adjustl(ades_entry))) .LE. 40) THEN
             ades_obs(obs_idx)%photProd = ades_entry
             ades_flg(obs_idx)%photProd = .TRUE.
          ELSE
             warn_str = .TRUE.

          ENDIF
       CASE('resMag')
          IF (trim(adjustl(ades_entry)) .EQ. "") THEN
             ades_flg(obs_idx)%resMag = .TRUE.
             !ades_resm_def(obs_idx) = .FALSE.
          ELSEIF (isreal(ades_entry)) THEN
             READ(ades_entry,*)ades_obs(obs_idx)%resMag
             ades_flg(obs_idx)%resMag = .TRUE.
             !ades_resm_def(obs_idx) = .TRUE.
          ELSE
             warn_numb = .TRUE.

          ENDIF
       CASE('selPhot')
          IF (len(trim(adjustl(ades_entry))) .LE. 1) THEN
             ades_obs(obs_idx)%selPhot = ades_entry
             ades_flg(obs_idx)%selPhot = .TRUE.
          ELSE
             warn_str = .TRUE.

          ENDIF
       CASE('sigMag')
          IF (trim(adjustl(ades_entry)) .EQ. "") THEN
             ades_flg(obs_idx)%sigMag = .TRUE.
          ELSEIF (isreal(ades_entry)) THEN
             READ(ades_entry,*)ades_obs(obs_idx)%sigMag
             ades_flg(obs_idx)%sigMag = .TRUE.
             ades_prec(obs_idx)%sigMag = get_prec(ades_entry)
          ELSE
             warn_numb = .TRUE.

          ENDIF
       CASE('biasMag')
          IF (trim(adjustl(ades_entry)) .EQ. "") THEN
             ades_flg(obs_idx)%biasMag = .TRUE.
          ELSEIF (isreal(ades_entry)) THEN
             READ(ades_entry,*)ades_obs(obs_idx)%biasMag
             ades_flg(obs_idx)%biasMag = .TRUE.
          ELSE
             warn_numb = .TRUE.

          ENDIF
       CASE('photMod')
          IF (len(trim(adjustl(ades_entry))) .LE. 8) THEN
             ades_obs(obs_idx)%photMod = ades_entry
             ades_flg(obs_idx)%photMod = .TRUE.
          ELSE
             warn_str = .TRUE.

          ENDIF
       CASE('resDelay')
          IF (trim(adjustl(ades_entry)) .EQ. "") THEN
             ades_flg(obs_idx)%resDelay = .TRUE.
             !ades_resc_def(obs_idx) = .FALSE.
          ELSEIF (isreal(ades_entry)) THEN
             READ(ades_entry,*)ades_obs(obs_idx)%resDelay
             ades_flg(obs_idx)%resDelay = .TRUE.
             !ades_resc_def(obs_idx) = .TRUE.
          ELSE
             warn_numb = .TRUE.

          ENDIF
       CASE('resDoppler')
          IF (trim(adjustl(ades_entry)) .EQ. "") THEN
             ades_flg(obs_idx)%resDoppler = .TRUE.
             !ades_resc_def(obs_idx) = .FALSE.
          ELSEIF (isreal(ades_entry)) THEN
             READ(ades_entry,*)ades_obs(obs_idx)%resDoppler
             ades_flg(obs_idx)%resDoppler = .TRUE.
             !ades_resc_def(obs_idx) = .TRUE.
          ELSE
             warn_numb = .TRUE.

          ENDIF
       CASE('selDelay')
          IF (len(trim(adjustl(ades_entry))) .LE. 1) THEN
             ades_obs(obs_idx)%selDelay = ades_entry
             ades_flg(obs_idx)%selDelay = .TRUE.
          ELSE
             warn_str = .TRUE.

          ENDIF
       CASE('selDoppler')
          IF (len(trim(adjustl(ades_entry))) .LE. 1) THEN
             ades_obs(obs_idx)%selDoppler = ades_entry
             ades_flg(obs_idx)%selDoppler = .TRUE.
          ELSE
             warn_str = .TRUE.

          ENDIF
       CASE('sigDelay')
          IF (trim(adjustl(ades_entry)) .EQ. "") THEN
             ades_flg(obs_idx)%sigDelay = .TRUE.
          ELSEIF (isreal(ades_entry)) THEN
             READ(ades_entry,*)ades_obs(obs_idx)%sigDelay
             ades_flg(obs_idx)%sigDelay = .TRUE.
             ades_prec(obs_idx)%sigDelay = get_prec(ades_entry)
          ELSE
             warn_numb = .TRUE.

          ENDIF
       CASE('sigDoppler')
          IF (trim(adjustl(ades_entry)) .EQ. "") THEN
             ades_flg(obs_idx)%sigDoppler = .TRUE.
          ELSEIF (isreal(ades_entry)) THEN
             READ(ades_entry,*)ades_obs(obs_idx)%sigDoppler
             ades_flg(obs_idx)%sigDoppler = .TRUE.
             ades_prec(obs_idx)%sigDoppler = get_prec(ades_entry)
          ELSE
             warn_numb = .TRUE.

          ENDIF

          ! other elements not belonging to any particular group
       CASE('deprecated')
          IF (len(trim(adjustl(ades_entry))) .LE. 1) THEN
             ades_obs(obs_idx)%deprecated = ades_entry
             ades_flg(obs_idx)%deprecated = .TRUE.
          ELSE
             warn_str = .TRUE.

          ENDIF
       CASE('localUse')
          IF (trim(adjustl(ades_entry)) .EQ. "") THEN  ! there are no restrictions on localUse 
             ades_flg(obs_idx)%localUse = .TRUE.
          ELSE
             ades_obs(obs_idx)%localUse = ades_entry
             ades_flg(obs_idx)%localUse = .TRUE.
          ENDIF
       END SELECT
    ENDIF
    !ENDIF

  END SUBROUTINE read_xml_elements


  ! ---------------------------------------- !
  ! SUBROUTINE read_rec_ades                 !
  ! It reads one record from ADES .psv file. !
  ! It skips the line if value is not a      !
  ! number or if string is too long          !
  ! ---------------------------------------- !

  SUBROUTINE read_rec_ades(unit,iline,keys,filename,ades_loc_obs,ades_loc_flg,ades_loc_prec,ades_loc_fit,ades_loc_header, &
       &                  ades_loc_key_rec,ades_loc_data_rec,value_found,old_count_keys,eof,error,skipped_obs_loc)

    INTEGER,              INTENT(IN)    :: unit                         ! ADES file unit
    INTEGER,              INTENT(INOUT) :: iline                        ! i-th record
    CHARACTER(LEN=30),    INTENT(INOUT) :: keys(nfields_ades)           ! array of keys read from line
    CHARACTER(LEN=200),   INTENT(IN)    :: filename                     ! ADES filename (either .obs_psv or .rwo_psv)
    TYPE(ades_observ),    INTENT(OUT)   :: ades_loc_obs                 ! observations from ADES file
    TYPE(ades_flags),     INTENT(OUT)   :: ades_loc_flg                 ! reading flags from ADES file
    TYPE(ades_precision), INTENT(OUT)   :: ades_loc_prec                ! entry precisions form ADES file
    CHARACTER(LEN=200),   INTENT(OUT)   :: ades_loc_header              ! header of ADES file
    CHARACTER(LEN=500),   INTENT(OUT)   :: ades_loc_key_rec             ! key record of ADES file
    CHARACTER(LEN=500),   INTENT(OUT)   :: ades_loc_data_rec            ! data record of ADES file
    TYPE(ades_fit_var),   INTENT(OUT)   :: ades_loc_fit                 ! local fit info from ADES file
    LOGICAL,              INTENT(OUT)   :: value_found                  ! TRUE if we find and store data record
    INTEGER,              INTENT(INOUT) :: old_count_keys               ! stored count_keys of previous key record
    LOGICAL,              INTENT(OUT)   :: eof                          ! end of file flag
    LOGICAL,              INTENT(OUT)   :: error                        ! error flag
    LOGICAL,              INTENT(OUT)   :: skipped_obs_loc

    INTEGER             :: status                    ! status of record
    CHARACTER(LEN=1000) :: record                    ! line record
    CHARACTER(LEN=90)   :: error_msg,warn_msg        ! error and warning messages
    CHARACTER(LEN=16)   :: rtn                       ! this routine
    CHARACTER(LEN=20)   :: int_err,FMT               ! auxiliary strings
    INTEGER             :: len_ie,len_chi,lf         ! auxiliary integer
    CHARACTER(LEN=40)   :: ades_entry(nfields_ades)  ! PSV entries in record
    LOGICAL             :: warn_numb,warn_str        ! if .TRUE. warning message is issued and the line is skipped

    CHARACTER(LEN=40)  :: time_str,aux_str
    INTEGER            :: i,j,old_rec_idx,len_str,pos_idx,idx_chi,idx_AS,idx_fw
    INTEGER            :: count_entries
    INTEGER            :: count_keys
    CHARACTER(LEN=15)  :: chi_str
    CHARACTER(LEN=1)   :: flg_str

    LOGICAL, EXTERNAL  :: isnum
    INTEGER, EXTERNAL  :: lench

    ! initialization
    ades_loc_obs=undefined_ades_observ
    ades_loc_flg=undefined_ades_flags
    ades_loc_prec=undefined_ades_precision
    ades_loc_fit=undefined_ades_fit_var
    ades_loc_header = ""
    ades_loc_key_rec = ""
    ades_loc_data_rec = ""
    value_found=.FALSE.
    eof=.FALSE.
    error=.FALSE.
    count_entries = 0
    count_keys = 0
    rtn = 'read_rec_ades'
    skipped_obs_loc = .FALSE.
    chi_str = ""
    flg_str = ""

    lf = lench(filename)

    ! read the line
    READ(unit,'(A)',IOSTAT=status) record

    ! check the status
    IF(status.GT.0)THEN
       WRITE(int_err,*) iline
       CALL rmsp(int_err,len_ie)
       error_msg='wrong format at line '//int_err(1:len_ie)
       WRITE(ierrou,*) 'ERROR: ',rtn,': ',trim(error_msg)
       error=.TRUE.
    ELSEIF(status.LT.0)THEN
       eof=.TRUE.
       RETURN
    END IF

    ! split fields in record and store them if record!=header
    ! if header store it in corresponding field
    IF (record(1:1) .EQ. "!" .OR. record(1:1) .EQ. "#") THEN
       ades_loc_header = record
       RETURN
       ! if record != header
    ELSE
       ! loop over each char in record
       old_rec_idx = 0
       DO i=1,len(record)

          ! find and store entry between pipes
          IF (record(i:i) .EQ. '|' .OR. i .EQ. len(record)) THEN
             count_entries = count_entries + 1
             ades_entry(count_entries)=trim(adjustl(record(old_rec_idx+1:i-1))) 
             IF (i .EQ. len(record)) ades_entry(count_entries)=trim(adjustl(record(old_rec_idx+1:i)))
             old_rec_idx = i

             ! if the forbidden localUse keyword is found, then issue an error 
             IF (ades_entry(count_entries) .EQ. "localUse")THEN
                WRITE(int_err,*) iline
                CALL rmsp(int_err,len_ie)
                error_msg=' localUse key found in input ADES file: '//filename(1:lf)//' line: '//int_err(1:len_ie)
                WRITE(ierrou,*) 'ERROR: ',rtn,': ',trim(error_msg)
                error=.TRUE.
             ENDIF

             ! check if record is a keyword or data record
             DO j=1,size(field_names)

                ! compare entry to elements of field_names array
                ! if key is recognized
                IF (ades_entry(count_entries) .EQ. trim(adjustl(field_names(j)))) THEN

                   ! store key record
                   ades_loc_key_rec = record

                   ! if first loop reset the keys array
                   IF (count_keys .EQ. 0) THEN
                      keys = ""
                      old_count_keys = 0
                   ENDIF

                   ! add entry to array of keys
                   count_keys = count_keys + 1
                   keys(count_keys) = ades_entry(count_entries)
                   EXIT
                ENDIF
             END DO
          ENDIF
       END DO

       ! if no key has been recognized it means we are in a data record
       IF (count_keys .EQ. 0) THEN

          ! store data record
          ades_loc_data_rec = record
          value_found = .TRUE.

          ! switches to write warning messages and skip line if entry is
          ! not a number or if the string is too long
          warn_numb = .FALSE.
          warn_str = .FALSE.

          ! loop over keys stored from latest keyword record
          ! store the corresponding value (and change flag and set precision) if key is useful
          DO i=1,size(keys)

             ! select key
             SELECT CASE (trim(adjustl(keys(i))))

                ! identification group (we do not consider artSat)
             CASE('permID')
                IF (len(trim(adjustl(ades_entry(i)))) .LE. 20) THEN
                   aux_str = ades_entry(i)
                   CALL rmsp(aux_str,len_str)
                   ades_loc_obs%permID = aux_str(1:len_str)
                   ades_loc_flg%permID = .TRUE.
                ELSE
                   warn_str = .TRUE.
                   EXIT
                ENDIF
             CASE('provID')
                IF (len(trim(adjustl(ades_entry(i)))) .LE. 30) THEN
                   len_str = len(trim(adjustl(ades_entry(i))))
                   aux_str = trim(adjustl(ades_entry(i)))
                   ! CALL rmsp(aux_str,len_str)
                   ades_loc_obs%provID = aux_str(1:len_str)
                   ades_loc_flg%provID = .TRUE.
                ELSE
                   warn_str = .TRUE.
                   EXIT
                ENDIF
             CASE('trkSub')
                IF (len(trim(adjustl(ades_entry(i)))) .LE. 8) THEN
                   ades_loc_obs%trkSub = ades_entry(i)
                   ades_loc_flg%trkSub = .TRUE.
                ELSE
                   warn_str = .TRUE.
                   EXIT
                ENDIF

                ! other elements not belonging to any particular group
             CASE('obsID')
                IF (len(trim(adjustl(ades_entry(i)))) .LE. 19) THEN
                   ades_loc_obs%obsID = ades_entry(i)
                   ades_loc_flg%obsID = .TRUE.
                ELSE
                   warn_str = .TRUE.
                   EXIT
                ENDIF
             CASE('trkID')
                IF (len(trim(adjustl(ades_entry(i)))) .LE. 12) THEN
                   ades_loc_obs%trkID = ades_entry(i)
                   ades_loc_flg%trkID = .TRUE.
                ELSE
                   warn_str = .TRUE.
                   EXIT
                ENDIF
             CASE('mode')
                IF (len(trim(adjustl(ades_entry(i)))) .LE. 3) THEN
                   ades_loc_obs%mode = ades_entry(i)
                   ades_loc_flg%mode = .TRUE.
                ELSE
                   warn_str = .TRUE.
                   EXIT
                ENDIF
             CASE('stn')
                IF (len(trim(adjustl(ades_entry(i)))) .LE. 4) THEN
                   ades_loc_obs%stn = ades_entry(i)
                   ades_loc_flg%stn = .TRUE.
                ELSE
                   warn_str = .TRUE.
                   EXIT
                ENDIF
             CASE('trx')
                IF (len(trim(adjustl(ades_entry(i)))) .LE. 4) THEN
                   ades_loc_obs%trx = ades_entry(i)
                   ades_loc_flg%trx = .TRUE.
                ELSE
                   warn_str = .TRUE.
                   EXIT
                ENDIF
             CASE('rcv')
                IF (len(trim(adjustl(ades_entry(i)))) .LE. 4) THEN
                   ades_loc_obs%rcv = ades_entry(i)
                   ades_loc_flg%rcv = .TRUE.
                ELSE
                   warn_str = .TRUE.
                   EXIT
                ENDIF

                ! Location group
             CASE('sys')
                IF (len(trim(adjustl(ades_entry(i)))) .LE. 7) THEN
                   ades_loc_obs%sys = ades_entry(i)
                   ades_loc_flg%sys = .TRUE.
                ELSE
                   warn_str = .TRUE.
                   EXIT
                ENDIF
             CASE('ctr')
                IF (trim(adjustl(ades_entry(i))) .EQ. "") THEN
                   ades_loc_flg%ctr = .TRUE.
                ELSEIF (isnum(trim(ades_entry(i)))) THEN
                   READ(ades_entry(i),*)ades_loc_obs%ctr
                   ades_loc_flg%ctr = .TRUE.
                ELSE
                   warn_numb = .TRUE.
                   EXIT
                ENDIF
             CASE('pos1')
                IF (trim(adjustl(ades_entry(i))) .EQ. "") THEN
                   ades_loc_flg%pos(1) = .TRUE.
                ELSEIF (isreal(ades_entry(i))) THEN
                   READ(ades_entry(i),*)ades_loc_obs%pos(1)
                   ades_loc_flg%pos(1) = .TRUE.
                ELSE
                   warn_numb = .TRUE.
                   EXIT
                ENDIF
             CASE('pos2')
                IF (trim(adjustl(ades_entry(i))) .EQ. "") THEN
                   ades_loc_flg%pos(2) = .TRUE.
                ELSEIF (isreal(ades_entry(i))) THEN
                   READ(ades_entry(i),*)ades_loc_obs%pos(2)
                   ades_loc_flg%pos(2) = .TRUE.
                ELSE
                   warn_numb = .TRUE.
                   EXIT
                ENDIF
             CASE('pos3')
                IF (trim(adjustl(ades_entry(i))) .EQ. "") THEN
                   ades_loc_flg%pos(3) = .TRUE.
                ELSEIF (isreal(ades_entry(i))) THEN
                   READ(ades_entry(i),*)ades_loc_obs%pos(3)
                   ades_loc_flg%pos(3) = .TRUE.
                ELSE
                   warn_numb = .TRUE.
                   EXIT
                ENDIF
             CASE('posCov11')
                IF (trim(adjustl(ades_entry(i))) .EQ. "") THEN
                   ades_loc_flg%posCov(1) = .TRUE.
                ELSEIF (isreal(ades_entry(i))) THEN
                   READ(ades_entry(i),*)ades_loc_obs%posCov(1)
                   ades_loc_flg%posCov(1) = .TRUE.
                ELSE
                   warn_numb = .TRUE.
                   EXIT
                ENDIF
             CASE('posCov12')
                IF (trim(adjustl(ades_entry(i))) .EQ. "") THEN
                   ades_loc_flg%posCov(2) = .TRUE.
                ELSEIF (isreal(ades_entry(i))) THEN
                   READ(ades_entry(i),*)ades_loc_obs%posCov(2)
                   ades_loc_flg%posCov(2) = .TRUE.
                ELSE
                   warn_numb = .TRUE.
                   EXIT
                ENDIF
             CASE('posCov13')
                IF (trim(adjustl(ades_entry(i))) .EQ. "") THEN
                   ades_loc_flg%posCov(3) = .TRUE.
                ELSEIF (isreal(ades_entry(i))) THEN
                   READ(ades_entry(i),*)ades_loc_obs%posCov(3)
                   ades_loc_flg%posCov(3) = .TRUE.
                ELSE
                   warn_numb = .TRUE.
                   EXIT
                ENDIF
             CASE('posCov22')
                IF (trim(adjustl(ades_entry(i))) .EQ. "") THEN
                   ades_loc_flg%posCov(4) = .TRUE.
                ELSEIF (isreal(ades_entry(i))) THEN
                   READ(ades_entry(i),*)ades_loc_obs%posCov(4)
                   ades_loc_flg%posCov(4) = .TRUE.
                ELSE
                   warn_numb = .TRUE.
                   EXIT
                ENDIF
             CASE('posCov23')
                IF (trim(adjustl(ades_entry(i))) .EQ. "") THEN
                   ades_loc_flg%posCov(5) = .TRUE.
                ELSEIF (isreal(ades_entry(i))) THEN
                   READ(ades_entry(i),*)ades_loc_obs%posCov(5)
                   ades_loc_flg%posCov(5) = .TRUE.
                ELSE
                   warn_numb = .TRUE.
                   EXIT
                ENDIF
             CASE('posCov33')
                IF (trim(adjustl(ades_entry(i))) .EQ. "") THEN
                   ades_loc_flg%posCov(6) = .TRUE.
                ELSEIF (isreal(ades_entry(i))) THEN
                   READ(ades_entry(i),*)ades_loc_obs%posCov(6)
                   ades_loc_flg%posCov(6) = .TRUE.
                ELSE
                   warn_numb = .TRUE.
                   EXIT
                ENDIF

                ! other elements not belonging to any particular group
             CASE('prog')
                IF (len(trim(adjustl(ades_entry(i)))) .LE. 2) THEN
                   ades_loc_obs%prog = ades_entry(i)
                   ades_loc_flg%prog = .TRUE.
                ELSE
                   warn_str = .TRUE.
                   EXIT
                ENDIF
             CASE('obsTime')
                time_str = ades_entry(i)
                CALL rmsp(time_str,len_str)
                ! lower threshold should be len_str = 23, but we relax it to 20
                IF (len_str .GE. 20 .AND. len_str .LE. 40) THEN
                   IF (isreal(time_str(1:4)) .AND. isreal(time_str(6:7)) .AND. isreal(time_str(9:10)) &
                        & .AND. isreal(time_str(12:13)) .AND. isreal(time_str(15:16)) .AND. &
                        & isreal(time_str(18:len_str-1))) THEN
                      ades_loc_obs%obsTime = ades_entry(i)
                      ades_loc_flg%obsTime = .TRUE.
                   ELSE
                      warn_numb = .TRUE.
                      EXIT
                   ENDIF
                ELSE
                   ! we write the warning message here since it may be triggered by a too short,
                   ! rather than too long, string
                   WRITE(int_err,*) iline
                   CALL rmsp(int_err,len_ie)
                   warn_msg = trim(adjustl(keys(i)))//' string too long or too short&
                        &in input ADES file: '//filename(1:lf)//' Skipping line '//int_err(1:len_ie)
                   WRITE(*,*) 'WARNING: ',rtn,': ',trim(warn_msg)
                   value_found = .FALSE.
                   RETURN
                ENDIF
             CASE('rmsTime')
                IF (trim(adjustl(ades_entry(i))) .EQ. "") THEN
                   ades_loc_flg%rmsTime = .TRUE.
                ELSEIF (isreal(ades_entry(i))) THEN
                   READ(ades_entry(i),*)ades_loc_obs%rmsTime
                   ades_loc_flg%rmsTime = .TRUE.
                   ades_loc_prec%rmsTime = get_prec(ades_entry(i))
                ELSE
                   warn_numb = .TRUE.
                   EXIT
                ENDIF
                ! observation group
             CASE('ra')
                IF (trim(adjustl(ades_entry(i))) .EQ. "") THEN
                   ades_loc_flg%ra = .TRUE.
                ELSEIF (isreal(ades_entry(i))) THEN
                   READ(ades_entry(i),*)ades_loc_obs%ra
                   ades_loc_flg%ra = .TRUE.
                   ades_loc_prec%ra = get_prec(ades_entry(i))
                   ! flag observation as optical
                   ades_loc_flg%opt = .TRUE.
                ELSE
                   warn_numb = .TRUE.
                   EXIT
                ENDIF
             CASE('dec')
                IF (trim(adjustl(ades_entry(i))) .EQ. "") THEN
                   ades_loc_flg%dec = .TRUE.
                ELSEIF (isreal(ades_entry(i))) THEN
                   READ(ades_entry(i),*)ades_loc_obs%dec
                   ades_loc_flg%dec = .TRUE.
                   ades_loc_prec%dec = get_prec(ades_entry(i))
                   ! flag observation as optical
                   ades_loc_flg%opt = .TRUE.
                ELSE
                   warn_numb = .TRUE.
                   EXIT
                ENDIF
             CASE('raStar')
                IF (trim(adjustl(ades_entry(i))) .EQ. "") THEN
                   ades_loc_flg%raStar = .TRUE.
                ELSEIF (isreal(ades_entry(i))) THEN
                   READ(ades_entry(i),*)ades_loc_obs%raStar
                   ades_loc_flg%raStar = .TRUE.
                   ades_loc_prec%raStar = get_prec(ades_entry(i))
                   ! flag observation as occultation
                   ades_loc_flg%occ = .TRUE.
                ELSE
                   warn_numb = .TRUE.
                   EXIT
                ENDIF
             CASE('decStar')
                IF (trim(adjustl(ades_entry(i))) .EQ. "") THEN
                   ades_loc_flg%decStar = .TRUE.
                ELSEIF (isreal(ades_entry(i))) THEN
                   READ(ades_entry(i),*)ades_loc_obs%decStar
                   ades_loc_flg%decStar = .TRUE.
                   ades_loc_prec%decStar = get_prec(ades_entry(i))
                   ! flag observation as occultation
                   ades_loc_flg%occ = .TRUE.
                ELSE
                   warn_numb = .TRUE.
                   EXIT
                ENDIF
             CASE('obsCenter')
                IF (len(trim(adjustl(ades_entry(i)))) .LE. 20) THEN
                   ades_loc_obs%obsCenter = ades_entry(i)
                   ades_loc_flg%obsCenter = .TRUE.
                   ! flag observation as offset
                   ades_loc_flg%off = .TRUE.
                ELSE
                   warn_str = .TRUE.
                   EXIT
                ENDIF
             CASE('deltaRA')
                IF (trim(adjustl(ades_entry(i))) .EQ. "") THEN
                   ades_loc_flg%deltaRA = .TRUE.
                ELSEIF (isreal(ades_entry(i))) THEN
                   READ(ades_entry(i),*)ades_loc_obs%deltaRA
                   ades_loc_flg%deltaRA = .TRUE.
                   ades_loc_prec%deltaRA = get_prec(ades_entry(i))
                ELSE
                   warn_numb = .TRUE.
                   EXIT
                ENDIF
             CASE('deltaDec')
                IF (trim(adjustl(ades_entry(i))) .EQ. "") THEN
                   ades_loc_flg%deltaDec = .TRUE.
                ELSEIF (isreal(ades_entry(i))) THEN
                   READ(ades_entry(i),*)ades_loc_obs%deltaDec
                   ades_loc_flg%deltaDec = .TRUE.
                   ades_loc_prec%deltaDec = get_prec(ades_entry(i))
                ELSE
                   warn_numb = .TRUE.
                   EXIT
                ENDIF
             CASE('dist')
                IF (trim(adjustl(ades_entry(i))) .EQ. "") THEN
                   ades_loc_flg%dist = .TRUE.
                ELSEIF (isreal(ades_entry(i))) THEN
                   READ(ades_entry(i),*)ades_loc_obs%dist
                   ades_loc_flg%dist = .TRUE.
                   ades_loc_prec%dist = get_prec(ades_entry(i))
                ELSE
                   warn_numb = .TRUE.
                   EXIT
                ENDIF
             CASE('pa')
                IF (trim(adjustl(ades_entry(i))) .EQ. "") THEN
                   ades_loc_flg%pa = .TRUE.
                ELSEIF (isreal(ades_entry(i))) THEN
                   READ(ades_entry(i),*)ades_loc_obs%pa
                   ades_loc_flg%pa = .TRUE.
                   ades_loc_prec%pa = get_prec(ades_entry(i))
                ELSE
                   warn_numb = .TRUE.
                   EXIT
                ENDIF
             CASE('rmsRA')
                IF (trim(adjustl(ades_entry(i))) .EQ. "") THEN
                   ades_loc_flg%rmsRA = .TRUE.
                ELSEIF (isreal(ades_entry(i))) THEN
                   READ(ades_entry(i),*)ades_loc_obs%rmsRA
                   ades_loc_flg%rmsRA = .TRUE.
                   ades_loc_prec%rmsRA = get_prec(ades_entry(i))
                ELSE
                   warn_numb = .TRUE.
                   EXIT
                ENDIF
             CASE('rmsDec')
                IF (trim(adjustl(ades_entry(i))) .EQ. "") THEN
                   ades_loc_flg%rmsDec = .TRUE.
                ELSEIF (isreal(ades_entry(i))) THEN
                   READ(ades_entry(i),*)ades_loc_obs%rmsDec
                   ades_loc_flg%rmsDec = .TRUE.
                   ades_loc_prec%rmsDec = get_prec(ades_entry(i))
                ELSE
                   warn_numb = .TRUE.
                   EXIT
                ENDIF
             CASE('rmsDist')
                IF (trim(adjustl(ades_entry(i))) .EQ. "") THEN
                   ades_loc_flg%rmsDist = .TRUE.
                ELSEIF (isreal(ades_entry(i))) THEN
                   READ(ades_entry(i),*)ades_loc_obs%rmsDist
                   ades_loc_flg%rmsDist = .TRUE.
                   ades_loc_prec%rmsDist = get_prec(ades_entry(i))
                ELSE
                   warn_numb = .TRUE.
                   EXIT
                ENDIF
             CASE('rmsPA')
                IF (trim(adjustl(ades_entry(i))) .EQ. "") THEN
                   ades_loc_flg%rmsPA = .TRUE.
                ELSEIF (isreal(ades_entry(i))) THEN
                   READ(ades_entry(i),*)ades_loc_obs%rmsPA
                   ades_loc_flg%rmsPA = .TRUE.
                   ades_loc_prec%rmsPA = get_prec(ades_entry(i))
                ELSE
                   warn_numb = .TRUE.
                   EXIT
                ENDIF
             CASE('rmsCorr')
                IF (trim(adjustl(ades_entry(i))) .EQ. "") THEN
                   ades_loc_flg%rmsCorr = .TRUE.
                ELSEIF (isreal(ades_entry(i))) THEN
                   READ(ades_entry(i),*)ades_loc_obs%rmsCorr
                   ades_loc_flg%rmsCorr = .TRUE.
                   ades_loc_prec%rmsCorr = get_prec(ades_entry(i))
                ELSE
                   warn_numb = .TRUE.
                   EXIT
                ENDIF
             CASE('delay')
                IF (trim(adjustl(ades_entry(i))) .EQ. "") THEN
                   ades_loc_flg%delay = .TRUE.
                ELSEIF (isreal(ades_entry(i))) THEN
                   READ(ades_entry(i),*)ades_loc_obs%delay
                   ades_loc_flg%delay = .TRUE.
                   ades_loc_prec%delay = get_prec(ades_entry(i))
                   ! flag observation as radar
                   ades_loc_flg%rad = .TRUE.
                ELSE
                   warn_numb = .TRUE.
                   EXIT
                ENDIF
             CASE('rmsDelay')
                IF (trim(adjustl(ades_entry(i))) .EQ. "") THEN
                   ades_loc_flg%rmsDelay = .TRUE.
                ELSEIF (isreal(ades_entry(i))) THEN
                   READ(ades_entry(i),*)ades_loc_obs%rmsDelay
                   ades_loc_flg%rmsDelay = .TRUE.
                   ades_loc_prec%rmsDelay = get_prec(ades_entry(i))
                ELSE
                   warn_numb = .TRUE.
                   EXIT
                ENDIF
             CASE('doppler')
                IF (trim(adjustl(ades_entry(i))) .EQ. "") THEN
                   ades_loc_flg%doppler = .TRUE.
                ELSEIF (isreal(ades_entry(i))) THEN
                   READ(ades_entry(i),*)ades_loc_obs%doppler
                   ades_loc_flg%doppler = .TRUE.
                   ades_loc_prec%doppler = get_prec(ades_entry(i))
                   ! flag observation as radar
                   ades_loc_flg%rad = .TRUE.
                ELSE
                   warn_numb = .TRUE.
                   EXIT
                ENDIF
             CASE('rmsDoppler')
                IF (trim(adjustl(ades_entry(i))) .EQ. "") THEN
                   ades_loc_flg%rmsDoppler = .TRUE.
                ELSEIF (isreal(ades_entry(i))) THEN
                   READ(ades_entry(i),*)ades_loc_obs%rmsDoppler
                   ades_loc_flg%rmsDoppler = .TRUE.
                   ades_loc_prec%rmsDoppler = get_prec(ades_entry(i))
                ELSE
                   warn_numb = .TRUE.
                   EXIT
                ENDIF

                ! other elements not belonging to any particular group
             CASE('astCat')
                IF (len(trim(adjustl(ades_entry(i)))) .LE. 8) THEN
                   ades_loc_obs%astCat = ades_entry(i)
                   ades_loc_flg%astCat = .TRUE.
                ELSE
                   warn_str = .TRUE.
                   EXIT
                ENDIF

                ! photometry group
             CASE('mag')
                IF (trim(adjustl(ades_entry(i))) .EQ. "") THEN
                   ades_loc_flg%mag = .TRUE.
                   !ades_loc_obs%mag = 9.9d9
                ELSEIF (isreal(ades_entry(i))) THEN
                   READ(ades_entry(i),*)ades_loc_obs%mag
                   ades_loc_flg%mag = .TRUE.
                   ades_loc_prec%mag = get_prec(ades_entry(i))
                ELSE
                   warn_numb = .TRUE.
                   EXIT
                ENDIF
             CASE('rmsMag')
                IF (trim(adjustl(ades_entry(i))) .EQ. "") THEN
                   ades_loc_flg%rmsMag = .TRUE.
                ELSEIF (isreal(ades_entry(i))) THEN
                   READ(ades_entry(i),*)ades_loc_obs%rmsMag
                   ades_loc_flg%rmsMag = .TRUE.
                   ades_loc_prec%rmsMag = get_prec(ades_entry(i))
                ELSE
                   warn_numb = .TRUE.
                   EXIT
                ENDIF
             CASE('band')
                IF (len(trim(adjustl(ades_entry(i)))) .LE. 3) THEN
                   ades_loc_obs%band = ades_entry(i)
                   ades_loc_flg%band = .TRUE.
                ELSE
                   warn_str = .TRUE.
                   EXIT
                ENDIF
             CASE('photCat')
                IF (len(trim(adjustl(ades_entry(i)))) .LE. 8) THEN
                   ades_loc_obs%photCat = ades_entry(i)
                   ades_loc_flg%photCat = .TRUE.
                ELSE
                   warn_str = .TRUE.
                   EXIT
                ENDIF
             CASE('photAp')
                IF (trim(adjustl(ades_entry(i))) .EQ. "") THEN
                   ades_loc_flg%photAp = .TRUE.
                ELSEIF (isreal(ades_entry(i))) THEN
                   READ(ades_entry(i),*)ades_loc_obs%photAp
                   ades_loc_flg%photAp = .TRUE.
                ELSE
                   warn_numb = .TRUE.
                   EXIT
                ENDIF
             CASE('nucMag')
                IF (trim(adjustl(ades_entry(i))) .EQ. "") THEN
                   ades_loc_flg%nucMag = .TRUE.
                ELSEIF (isnum(trim(ades_entry(i)))) THEN
                   READ(ades_entry(i),*)ades_loc_obs%nucMag
                   ades_loc_flg%nucMag = .TRUE.
                ELSE
                   warn_numb = .TRUE.
                   EXIT
                ENDIF

                ! other elements not belonging to any particular group
             CASE('logSNR')
                IF (trim(adjustl(ades_entry(i))) .EQ. "") THEN
                   ades_loc_flg%logSNR = .TRUE.
                ELSEIF (isreal(ades_entry(i))) THEN
                   READ(ades_entry(i),*)ades_loc_obs%logSNR
                   ades_loc_flg%logSNR = .TRUE.
                ELSE
                   warn_numb = .TRUE.
                   EXIT
                ENDIF
             CASE('seeing')
                IF (trim(adjustl(ades_entry(i))) .EQ. "") THEN
                   ades_loc_flg%seeing = .TRUE.
                ELSEIF (isreal(ades_entry(i))) THEN
                   READ(ades_entry(i),*)ades_loc_obs%seeing
                   ades_loc_flg%seeing = .TRUE.
                ELSE
                   warn_numb = .TRUE.
                   EXIT
                ENDIF
             CASE('exp')
                IF (trim(adjustl(ades_entry(i))) .EQ. "") THEN
                   ades_loc_flg%exp = .TRUE.
                ELSEIF (isreal(ades_entry(i))) THEN
                   READ(ades_entry(i),*)ades_loc_obs%exp
                   ades_loc_flg%exp = .TRUE.
                ELSE
                   warn_numb = .TRUE.
                   EXIT
                ENDIF
             CASE('rmsFit')
                IF (trim(adjustl(ades_entry(i))) .EQ. "") THEN
                   ades_loc_flg%rmsFit = .TRUE.
                ELSEIF (isreal(ades_entry(i))) THEN
                   READ(ades_entry(i),*)ades_loc_obs%rmsFit
                   ades_loc_flg%rmsFit = .TRUE.
                ELSE
                   warn_numb = .TRUE.
                   EXIT
                ENDIF
             CASE('nStars')
                IF (trim(adjustl(ades_entry(i))) .EQ. "") THEN
                   ades_loc_flg%nStars = .TRUE.
                ELSEIF (isnum(trim(ades_entry(i)))) THEN
                   READ(ades_entry(i),*)ades_loc_obs%nStars
                   ades_loc_flg%nStars = .TRUE.
                ELSE
                   warn_numb = .TRUE.
                   EXIT
                ENDIF
             CASE('com')
                IF (trim(adjustl(ades_entry(i))) .EQ. "") THEN
                   ades_loc_flg%com = .TRUE.
                ELSEIF (isnum(trim(ades_entry(i)))) THEN
                   READ(ades_entry(i),*)ades_loc_obs%com
                   ades_loc_flg%com = .TRUE.
                ELSE
                   warn_numb = .TRUE.
                   EXIT
                ENDIF
             CASE('frq')
                IF (trim(adjustl(ades_entry(i))) .EQ. "") THEN
                   ades_loc_flg%frq = .TRUE.
                ELSEIF (isreal(ades_entry(i))) THEN
                   READ(ades_entry(i),*)ades_loc_obs%frq
                   ades_loc_flg%frq = .TRUE.
                   ades_loc_prec%frq = get_prec(ades_entry(i))
                ELSE
                   warn_numb = .TRUE.
                   EXIT
                ENDIF
             CASE('ref')
                IF (len(trim(adjustl(ades_entry(i)))) .LE. 16) THEN
                   ades_loc_obs%ref = ades_entry(i)
                   ades_loc_flg%ref = .TRUE.
                ELSE
                   warn_str = .TRUE.
                   EXIT
                ENDIF
             CASE('disc')
                IF (len(trim(adjustl(ades_entry(i)))) .LE. 1) THEN
                   ades_loc_obs%disc = ades_entry(i)
                   ades_loc_flg%disc = .TRUE.
                ELSE
                   warn_str = .TRUE.
                   EXIT
                ENDIF
             CASE('subFrm')
                IF (len(trim(adjustl(ades_entry(i)))) .LE. 20) THEN
                   ades_loc_obs%subFrm = ades_entry(i)
                   ades_loc_flg%subFrm = .TRUE.
                ELSE
                   warn_str = .TRUE.
                   EXIT
                ENDIF
             CASE('subFmt')
                IF (len(trim(adjustl(ades_entry(i)))) .LE. 4) THEN
                   ades_loc_obs%subFmt = ades_entry(i)
                   ades_loc_flg%subFmt = .TRUE.
                ELSE
                   warn_str = .TRUE.
                   EXIT
                ENDIF

                ! precision group
             CASE('precTime')
                IF (trim(adjustl(ades_entry(i))) .EQ. "") THEN
                   ades_loc_flg%precTime = .TRUE.
                ELSEIF (isnum(trim(ades_entry(i)))) THEN
                   READ(ades_entry(i),*)ades_loc_obs%precTime
                   ades_loc_flg%precTime = .TRUE.
                ELSE
                   warn_numb = .TRUE.
                   EXIT
                ENDIF
             CASE('precRA')
                IF (trim(adjustl(ades_entry(i))) .EQ. "") THEN
                   ades_loc_flg%precRA = .TRUE.
                ELSEIF (isreal(ades_entry(i))) THEN
                   READ(ades_entry(i),*)ades_loc_obs%precRA
                   ades_loc_flg%precRA = .TRUE.
                   ades_loc_prec%precRA = get_prec(ades_entry(i))
                ELSE
                   warn_numb = .TRUE.
                   EXIT
                ENDIF
             CASE('precDec')
                IF (trim(adjustl(ades_entry(i))) .EQ. "") THEN
                   ades_loc_flg%precDec = .TRUE.
                ELSEIF (isreal(ades_entry(i))) THEN
                   READ(ades_entry(i),*)ades_loc_obs%precDec
                   ades_loc_flg%precDec = .TRUE.
                   ades_loc_prec%precDec = get_prec(ades_entry(i))
                ELSE
                   warn_numb = .TRUE.
                   EXIT
                ENDIF

                ! other elements not belonging to any particular group
             CASE('uncTime')
                IF (trim(adjustl(ades_entry(i))) .EQ. "") THEN
                   ades_loc_flg%uncTime = .TRUE.
                ELSEIF (isreal(ades_entry(i))) THEN
                   READ(ades_entry(i),*)ades_loc_obs%uncTime
                   ades_loc_flg%uncTime = .TRUE.
                ELSE
                   warn_numb = .TRUE.
                   EXIT
                ENDIF
             CASE('notes')
                IF (len(trim(adjustl(ades_entry(i)))) .LE. 6) THEN
                   ades_loc_obs%notes = ades_entry(i)
                   ades_loc_flg%notes = .TRUE.
                ELSE
                   warn_str = .TRUE.
                   EXIT
                ENDIF
             CASE('remarks')
                IF (len(trim(adjustl(ades_entry(i)))) .LE. 200) THEN
                   ades_loc_obs%remarks = ades_entry(i)
                   ades_loc_flg%remarks = .TRUE.
                   idx_chi = INDEX(ades_entry(i),'chi =')
                   idx_AS = INDEX(ades_entry(i),'A = ') + INDEX(ades_entry(i),'S = ')
                   idx_fw = INDEX(ades_entry(i),'RMS_F = ')
                   IF (idx_chi .GT. 0) THEN
                      IF (idx_AS .GT. 0) THEN
                         chi_str = ades_entry(i)(idx_chi+5:idx_AS-1)
                         CALL rmsp(chi_str,len_chi)
                         IF (len_chi .GT. 9) THEN
                            READ(chi_str,'(1P,E9.2)') ades_loc_fit%chi
                         ELSE
                            WRITE(FMT,'(i20)') len_chi
                            READ(chi_str,"(F"//trim(adjustl(FMT))//".2)") ades_loc_fit%chi
                         ENDIF
                      ELSE
                         warn_msg = 'Missing astrometric fit flag in input ADES file: '//filename(1:lf)
                         ! write warning message
                         WRITE(*,*) 'WARNING: ',rtn,': ',trim(warn_msg)
                      ENDIF
                      IF (idx_fw .LE. 0) THEN
                         warn_msg = 'Missing flags for forced weights in input ADES file: '//filename(1:lf)
                         ! write warning message
                         WRITE(*,*) 'WARNING: ',rtn,': ',trim(warn_msg)
                      ENDIF
                   ENDIF
                   IF (idx_AS .GT. 0) THEN
                      flg_str = ades_entry(i)(idx_AS:idx_AS)
                      IF (flg_str .EQ. "S") ades_loc_fit%radar_flag = ades_entry(i)(idx_AS+4:idx_AS+4)
                      IF (flg_str .EQ. "A") THEN
                         ades_loc_fit%ast_flag = ades_entry(i)(idx_AS+4:idx_AS+4)
                         ades_loc_fit%mag_flag = ades_entry(i)(idx_AS+10:idx_AS+10)
                      ENDIF
                   ENDIF
                   IF (idx_fw .GT. 0) THEN
                      IF (.NOT. ades_loc_flg%rad) THEN
                         READ(ades_entry(i)(idx_fw+8:idx_fw+8),'(L1)') ades_loc_fit%force_w_flags(1)
                         READ(ades_entry(i)(idx_fw+10:idx_fw+10),'(L1)') ades_loc_fit%force_w_flags(2)
                      ELSE
                         IF (ades_loc_flg%delay) READ(ades_entry(i)(idx_fw+8:idx_fw+8),'(L1)') ades_loc_fit%force_w_flags(1)
                         IF (ades_loc_flg%doppler) READ(ades_entry(i)(idx_fw+8:idx_fw+8),'(L1)') ades_loc_fit%force_w_flags(2)
                      ENDIF
                   ENDIF
                ELSE
                   warn_str = .TRUE.
                   EXIT
                ENDIF

                ! Residuals group
             CASE('orbProd')
                IF (len(trim(adjustl(ades_entry(i)))) .LE. 40) THEN
                   ades_loc_obs%orbProd = ades_entry(i)
                   ades_loc_flg%orbProd = .TRUE.
                ELSE
                   warn_str = .TRUE.
                   EXIT
                ENDIF
             CASE('orbID')
                IF (len(trim(adjustl(ades_entry(i)))) .LE. 40) THEN
                   ades_loc_obs%orbID = ades_entry(i)
                   ades_loc_flg%orbID = .TRUE.
                ELSE
                   warn_str = .TRUE.
                   EXIT
                ENDIF
             CASE('resRA')
                IF (trim(adjustl(ades_entry(i))) .EQ. "") THEN
                   ades_loc_flg%resRA = .TRUE.
                   ades_loc_fit%ades_resc_def = .FALSE.
                ELSEIF (isreal(ades_entry(i))) THEN
                   READ(ades_entry(i),*)ades_loc_obs%resRA
                   ades_loc_flg%resRA = .TRUE.
                   ades_loc_fit%ades_resc_def = .TRUE.
                ELSE
                   warn_numb = .TRUE.
                   EXIT
                ENDIF
             CASE('resDec')
                IF (trim(adjustl(ades_entry(i))) .EQ. "") THEN
                   ades_loc_flg%resDec = .TRUE.
                   ades_loc_fit%ades_resc_def = .FALSE.
                ELSEIF (isreal(ades_entry(i))) THEN
                   READ(ades_entry(i),*)ades_loc_obs%resDec
                   ades_loc_flg%resDec = .TRUE.
                ELSE
                   warn_numb = .TRUE.
                   EXIT
                ENDIF
             CASE('selAst')
                IF (len(trim(adjustl(ades_entry(i)))) .LE. 1) THEN
                   ades_loc_obs%selAst = ades_entry(i)
                   ades_loc_flg%selAst = .TRUE.
                ELSE
                   warn_str = .TRUE.
                   EXIT
                ENDIF
             CASE('sigRA')
                IF (trim(adjustl(ades_entry(i))) .EQ. "") THEN
                   ades_loc_flg%sigRA = .TRUE.
                ELSEIF (isreal(ades_entry(i))) THEN
                   READ(ades_entry(i),*)ades_loc_obs%sigRA
                   ades_loc_flg%sigRA = .TRUE.
                   ades_loc_prec%sigRA = get_prec(ades_entry(i))
                ELSE
                   warn_numb = .TRUE.
                   EXIT
                ENDIF
             CASE('sigDec')
                IF (trim(adjustl(ades_entry(i))) .EQ. "") THEN
                   ades_loc_flg%sigDec = .TRUE.
                ELSEIF (isreal(ades_entry(i))) THEN
                   READ(ades_entry(i),*)ades_loc_obs%sigDec
                   ades_loc_flg%sigDec = .TRUE.
                   ades_loc_prec%sigDec = get_prec(ades_entry(i))
                ELSE
                   warn_numb = .TRUE.
                   EXIT
                ENDIF
             CASE('sigCorr')
                IF (trim(adjustl(ades_entry(i))) .EQ. "") THEN
                   ades_loc_flg%sigCorr = .TRUE.
                ELSEIF (isreal(ades_entry(i))) THEN
                   READ(ades_entry(i),*)ades_loc_obs%sigCorr
                   ades_loc_flg%sigCorr = .TRUE.
                   ades_loc_prec%sigCorr = get_prec(ades_entry(i))
                ELSE
                   warn_numb = .TRUE.
                   EXIT
                ENDIF
             CASE('sigTime')
                IF (trim(adjustl(ades_entry(i))) .EQ. "") THEN
                   ades_loc_flg%sigTime = .TRUE.
                ELSEIF (isreal(ades_entry(i))) THEN
                   READ(ades_entry(i),*)ades_loc_obs%sigTime
                   ades_loc_flg%sigTime = .TRUE.
                ELSE
                   warn_numb = .TRUE.
                   EXIT
                ENDIF
             CASE('biasRA')
                IF (trim(adjustl(ades_entry(i))) .EQ. "") THEN
                   ades_loc_flg%biasRA = .TRUE.
                ELSEIF (isreal(ades_entry(i))) THEN
                   READ(ades_entry(i),*)ades_loc_obs%biasRA
                   ades_loc_flg%biasRA = .TRUE.
                ELSE
                   warn_numb = .TRUE.
                   EXIT
                ENDIF
             CASE('biasDec')
                IF (trim(adjustl(ades_entry(i))) .EQ. "") THEN
                   ades_loc_flg%biasDec = .TRUE.
                ELSEIF (isreal(ades_entry(i))) THEN
                   READ(ades_entry(i),*)ades_loc_obs%biasDec
                   ades_loc_flg%biasDec = .TRUE.
                ELSE
                   warn_numb = .TRUE.
                   EXIT
                ENDIF
             CASE('biasTime')
                IF (trim(adjustl(ades_entry(i))) .EQ. "") THEN
                   ades_loc_flg%biasTime = .TRUE.
                ELSEIF (isreal(ades_entry(i))) THEN
                   READ(ades_entry(i),*)ades_loc_obs%biasTime
                   ades_loc_flg%biasTime = .TRUE.
                ELSE
                   warn_numb = .TRUE.
                   EXIT
                ENDIF
             CASE('photProd')
                IF (len(trim(adjustl(ades_entry(i)))) .LE. 40) THEN
                   ades_loc_obs%photProd = ades_entry(i)
                   ades_loc_flg%photProd = .TRUE.
                ELSE
                   warn_str = .TRUE.
                   EXIT
                ENDIF
             CASE('resMag')
                IF (trim(adjustl(ades_entry(i))) .EQ. "") THEN
                   ades_loc_flg%resMag = .TRUE.
                   ades_loc_fit%ades_resm_def = .FALSE.
                ELSEIF (isreal(ades_entry(i))) THEN
                   READ(ades_entry(i),*)ades_loc_obs%resMag
                   ades_loc_flg%resMag = .TRUE.
                   ades_loc_fit%ades_resm_def = .TRUE.
                ELSE
                   warn_numb = .TRUE.
                   EXIT
                ENDIF
             CASE('selPhot')
                IF (len(trim(adjustl(ades_entry(i)))) .LE. 1) THEN
                   ades_loc_obs%selPhot = ades_entry(i)
                   ades_loc_flg%selPhot = .TRUE.
                ELSE
                   warn_str = .TRUE.
                   EXIT
                ENDIF
             CASE('sigMag')
                IF (trim(adjustl(ades_entry(i))) .EQ. "") THEN
                   ades_loc_flg%sigMag = .TRUE.
                ELSEIF (isreal(ades_entry(i))) THEN
                   READ(ades_entry(i),*)ades_loc_obs%sigMag
                   ades_loc_flg%sigMag = .TRUE.
                   ades_loc_prec%sigMag = get_prec(ades_entry(i))
                ELSE
                   warn_numb = .TRUE.
                   EXIT
                ENDIF
             CASE('biasMag')
                IF (trim(adjustl(ades_entry(i))) .EQ. "") THEN
                   ades_loc_flg%biasMag = .TRUE.
                ELSEIF (isreal(ades_entry(i))) THEN
                   READ(ades_entry(i),*)ades_loc_obs%biasMag
                   ades_loc_flg%biasMag = .TRUE.
                ELSE
                   warn_numb = .TRUE.
                   EXIT
                ENDIF
             CASE('photMod')
                IF (len(trim(adjustl(ades_entry(i)))) .LE. 8) THEN
                   ades_loc_obs%photMod = ades_entry(i)
                   ades_loc_flg%photMod = .TRUE.
                ELSE
                   warn_str = .TRUE.
                   EXIT
                ENDIF
             CASE('resDelay')
                IF (trim(adjustl(ades_entry(i))) .EQ. "") THEN
                   ades_loc_flg%resDelay = .TRUE.
                   ades_loc_fit%ades_resc_def = .FALSE.
                ELSEIF (isreal(ades_entry(i))) THEN
                   READ(ades_entry(i),*)ades_loc_obs%resDelay
                   ades_loc_flg%resDelay = .TRUE.
                   ades_loc_fit%ades_resc_def = .TRUE.
                ELSE
                   warn_numb = .TRUE.
                   EXIT
                ENDIF
             CASE('resDoppler')
                IF (trim(adjustl(ades_entry(i))) .EQ. "") THEN
                   ades_loc_flg%resDoppler = .TRUE.
                   ades_loc_fit%ades_resc_def = .FALSE.
                ELSEIF (isreal(ades_entry(i))) THEN
                   READ(ades_entry(i),*)ades_loc_obs%resDoppler
                   ades_loc_flg%resDoppler = .TRUE.
                   ades_loc_fit%ades_resc_def = .TRUE.
                ELSE
                   warn_numb = .TRUE.
                   EXIT
                ENDIF
             CASE('selDelay')
                IF (len(trim(adjustl(ades_entry(i)))) .LE. 1) THEN
                   ades_loc_obs%selDelay = ades_entry(i)
                   ades_loc_flg%selDelay = .TRUE.
                ELSE
                   warn_str = .TRUE.
                   EXIT
                ENDIF
             CASE('selDoppler')
                IF (len(trim(adjustl(ades_entry(i)))) .LE. 1) THEN
                   ades_loc_obs%selDoppler = ades_entry(i)
                   ades_loc_flg%selDoppler = .TRUE.
                ELSE
                   warn_str = .TRUE.
                   EXIT
                ENDIF
             CASE('sigDelay')
                IF (trim(adjustl(ades_entry(i))) .EQ. "") THEN
                   ades_loc_flg%sigDelay = .TRUE.
                ELSEIF (isreal(ades_entry(i))) THEN
                   READ(ades_entry(i),*)ades_loc_obs%sigDelay
                   ades_loc_flg%sigDelay = .TRUE.
                   ades_loc_prec%sigDelay = get_prec(ades_entry(i))
                ELSE
                   warn_numb = .TRUE.
                   EXIT
                ENDIF
             CASE('sigDoppler')
                IF (trim(adjustl(ades_entry(i))) .EQ. "") THEN
                   ades_loc_flg%sigDoppler = .TRUE.
                ELSEIF (isreal(ades_entry(i))) THEN
                   READ(ades_entry(i),*)ades_loc_obs%sigDoppler
                   ades_loc_flg%sigDoppler = .TRUE.
                   ades_loc_prec%sigDoppler = get_prec(ades_entry(i))
                ELSE
                   warn_numb = .TRUE.
                   EXIT
                ENDIF

                ! other elements not belonging to any particular group
             CASE('deprecated')
                IF (len(trim(adjustl(ades_entry(i)))) .LE. 1) THEN
                   ades_loc_obs%deprecated = ades_entry(i)
                   ades_loc_flg%deprecated = .TRUE.
                ELSE
                   warn_str = .TRUE.
                   EXIT
                ENDIF

             END SELECT

          END DO

          ! if some inconsistency has been found, write warning message and skip line
          IF (warn_numb .OR. warn_str) THEN
             WRITE(int_err,*) iline
             CALL rmsp(int_err,len_ie)
             IF (warn_numb) THEN
                warn_msg = trim(adjustl(keys(i)))//' value not a number &
                     &in input ADES file: '//filename(1:lf)//' Skipping line '//int_err(1:len_ie)
             ELSE IF (warn_str) THEN
                warn_msg = trim(adjustl(keys(i)))//' string too long &
                     &in input ADES file: '//filename(1:lf)//' Skipping line '//int_err(1:len_ie)
             ENDIF
             ! write warning message
             WRITE(*,*) 'WARNING: ',rtn,': ',trim(warn_msg)
             ! set value found to FALSE to skip record storage
             value_found = .FALSE.
             RETURN
          ENDIF

          ! check consistency of fields
          CALL check_fields(filename,ades_loc_obs,ades_loc_flg,value_found,error,ILINE=iline)

          IF (.NOT. value_found) THEN
             skipped_obs_loc = .TRUE.
             ades_loc_data_rec = ""
          ENDIF

          ! store old_count_keys for later check
       ELSEIF (count_keys .GT. 0) THEN
          old_count_keys = count_keys

          ! error if not all keys have been recognized
       ELSEIF (count_keys .GT. 0 .AND. count_entries .NE. count_keys) THEN
          WRITE(int_err,*) iline
          CALL rmsp(int_err,len_ie)
          error_msg='key entry/ies not recognized &
               &in input ADES file: '//filename(1:lf)//' line: '//int_err(1:len_ie)
          WRITE(ierrou,*) 'ERROR: ',rtn,': ',trim(error_msg)
          error=.TRUE.

       ENDIF

       ! error if there are different numbers of keys and values
       IF (count_entries .NE. old_count_keys) THEN
          WRITE(int_err,*) iline
          CALL rmsp(int_err,len_ie)
          error_msg='different number of keys and values'//&
               &' in input ADES file: '//filename(1:lf)//' line: '//int_err(1:len_ie)
          WRITE(ierrou,*) 'ERROR: ',rtn,': ',trim(error_msg)
          error=.TRUE.
       ENDIF

    ENDIF

  END SUBROUTINE read_rec_ades

  ! ---------------------------------------------- !
  ! SUBROUTINE check_fields                        !
  ! It checks that the keys and the values in      !
  ! input file are consistent with the ADES format !
  ! specifications. RETURN if required keys are    !
  ! not present and skip line if values are not    !
  ! consistent with format.                        !
  ! Do not check for keys that are not allowed     !
  ! to be present at the same time, but only for   !
  ! required keys.                                 !
  ! ---------------------------------------------- !

  SUBROUTINE check_fields(filename,ades_loc_obs,ades_loc_flg,value_found,error,iline,obs_idx)

    CHARACTER(LEN=200), INTENT(IN)     :: filename            ! ADES filename (either .obs_psv or .rwo_psv)
    TYPE(ades_observ),  INTENT(IN)     :: ades_loc_obs        ! ADES observations
    TYPE(ades_flags),   INTENT(IN)     :: ades_loc_flg        ! ADES reading flags
    LOGICAL,            INTENT(INOUT)  :: value_found         ! TRUE if we find and store data record
    LOGICAL,            INTENT(INOUT)  :: error               ! error flag
    INTEGER,     INTENT(IN), OPTIONAL  :: iline               ! i-th record (psv case)
    INTEGER,     INTENT(IN), OPTIONAL  :: obs_idx             ! observation number (xml case)

    CHARACTER(LEN=20)    :: int_err
    INTEGER              :: len_ie,i,len_str,lf
    CHARACTER(LEN=13)    :: rtn = "check_fields"
    CHARACTER(LEN=80)    :: err_strk,warn_strv
    CHARACTER(LEN=90)    :: error_msg,warn_msg,warn_skip      ! error and warning messages
    CHARACTER(LEN=40)    :: time_str
    INTEGER              :: month,day,hour,minute
    LOGICAL              :: location_flg(5),optastr_res_flg(7),optphot_res_flg(5),rad_res_flg(5),opt_prec_flg(3)
    REAL(KIND=dkind)     :: second

    ! prepare error messages
    IF (PRESENT(iline)) THEN   ! psv case: use line number
       WRITE(int_err,*) iline
       CALL rmsp(int_err,len_ie)
       err_strk = 'in input ADES file: '//trim(filename)//' line: '//int_err(1:len_ie)
       warn_strv = 'in input ADES file: '//trim(filename)//' line: '//int_err(1:len_ie)
       warn_skip = 'Skipping line: '//int_err(1:len_ie)
    ELSEIF (PRESENT(obs_idx)) THEN    ! xml case: use observation number
       WRITE(int_err,*) obs_idx
       CALL rmsp(int_err,len_ie)
       err_strk = 'in input ADES file: '//trim(filename)//' observation: '//int_err(1:len_ie)
       warn_strv = 'in input ADES file: '//trim(filename)//' observation: '//int_err(1:len_ie)
       warn_skip = ' Skipping observation: '//int_err(1:len_ie)
    ELSE
       WRITE(ierrou,*) 'ERROR: ',rtn,': Missing both iline and obs_idx input parameters.'
       error=.TRUE.
       RETURN
    ENDIF

    ! check whether the type of observation has been recognized
    IF (.NOT. ades_loc_flg%opt .AND. .NOT. ades_loc_flg%off .AND. .NOT. ades_loc_flg%occ .AND. .NOT. ades_loc_flg%rad) THEN
       IF (PRESENT(iline)) THEN 
          error_msg='Observation type (optical, offset, occultation, radar) not recognized &
               &in input ADES file: '//trim(filename)// ' line: '//int_err(1:len_ie)
       ELSEIF (PRESENT(obs_idx)) THEN
          error_msg='Observation type (optical, offset, occultation, radar) not recognized &
               &in input ADES file: '//trim(filename)// ' observation: '//int_err(1:len_ie)
       ENDIF
       WRITE(ierrou,*) 'ERROR: ',rtn,': ',error_msg
       error=.TRUE.
       RETURN
    ENDIF
    ! check whether the two or more types of observations have been recognized
    IF ((ades_loc_flg%opt .AND. (ades_loc_flg%off .OR. ades_loc_flg%occ .OR. ades_loc_flg%rad)) .OR. &
         & (ades_loc_flg%off .AND. (ades_loc_flg%occ .OR. ades_loc_flg%rad)) .OR. &
         & (ades_loc_flg%occ .AND. ades_loc_flg%rad)) THEN
       IF (PRESENT(iline)) THEN 
          error_msg='Multiple observation types (optical, offset, occultation, radar) recognized &
               &in input ADES file: '//trim(filename)// ' line: '//int_err(1:len_ie)
       ELSEIF (PRESENT(obs_idx)) THEN
          error_msg='Multiple observation types (optical, offset, occultation, radar) recognized &
               &in input ADES file: '//trim(filename)// ' observation: '//int_err(1:len_ie)
       ENDIF
       WRITE(ierrou,*) 'ERROR: ',rtn,': ',error_msg
       error=.TRUE.
       RETURN
    ENDIF

    ! skip offset observations (not supported)
    IF (ades_loc_flg%off) THEN
       ! write warning message
       warn_msg='Offset observation not supported '//warn_strv
       WRITE(*,*) 'WARNING: ',rtn,': ',warn_msg
       ! set value found to FALSE to skip record storage
       value_found = .FALSE.
       WRITE(*,*) warn_skip
       RETURN
    ENDIF

    ! deprecated observations: only warning, not skipping obs (deprecated obs are underweighted)
    IF (.NOT. ades_loc_flg%rad .AND. ades_loc_flg%deprecated .AND. (ades_loc_obs%deprecated .EQ. "X")) THEN
       warn_msg='Deprecated observation found '//warn_strv
       WRITE(*,*) 'WARNING: ',rtn,': ',warn_msg
       !value_found = .FALSE.
       !WRITE(*,*) warn_skip
       !RETURN
    ENDIF
    ! check whether observations is radar and has forbidden deprecated entry
    IF (ades_loc_flg%rad .AND. ades_loc_flg%deprecated) THEN
       error_msg='Radar observation has forbidden deprecated entry '//err_strk
       WRITE(ierrou,*) 'ERROR: ',rtn,': ',error_msg
       error=.TRUE.
       RETURN
    ENDIF

    ! check for permID provID and trkSUB keys (we do not consider artSat)
    ! if radar obs, then only permID and/or provID are permitted
    IF (ades_loc_flg%rad) THEN
       IF (.NOT. ades_loc_flg%permID .AND. .NOT. ades_loc_flg%provID) THEN
          error_msg='permID/provID key not found '//err_strk
          WRITE(ierrou,*) 'ERROR: ',rtn,': ',error_msg
          error=.TRUE.
          RETURN
       ENDIF
    ELSE
       ! if optical obs
       IF (.NOT. ades_loc_flg%permID .AND. .NOT. ades_loc_flg%provID .AND. .NOT. ades_loc_flg%trkSub) THEN
          error_msg='permID/provID/trkSub key not found '//err_strk
          WRITE(ierrou,*) 'ERROR: ',rtn,': ',error_msg
          error=.TRUE.
          RETURN
       ENDIF
    ENDIF

    ! mode must be present for optical and offset obs
    IF (.NOT. ades_loc_flg%rad .AND. .NOT. ades_loc_flg%occ .AND. .NOT. ades_loc_flg%mode) THEN
       error_msg='mode key not found '//err_strk
       WRITE(ierrou,*) 'ERROR: ',rtn,': ',error_msg
       error=.TRUE.
       RETURN
    ENDIF
    ! check mode values
    IF (.NOT. ades_loc_flg%rad .AND. .NOT. ades_loc_flg%occ .AND. ades_loc_flg%mode) THEN
       IF (ALL(mode_values .NE. ades_loc_obs%mode)) THEN
          warn_msg='Forbidden mode value '//warn_strv
          WRITE(*,*) 'WARNING: ',rtn,': ',warn_msg
          !value_found = .FALSE.
          !WRITE(*,*) warn_skip
          !RETURN
       ELSEIF (ades_loc_obs%mode .EQ. "UNK" .AND. ades_loc_obs%deprecated .NE. "X") THEN
          warn_msg='Unknown mode value '//warn_strv
          WRITE(*,*) 'WARNING: ',rtn,': ',warn_msg
       ENDIF
    ENDIF

    ! stn must be present for all optical obs
    IF (.NOT. ades_loc_flg%rad .AND. .NOT. ades_loc_flg%stn) THEN
       error_msg='stn key not found '//err_strk
       WRITE(ierrou,*) 'ERROR: ',rtn,': ',error_msg
       error=.TRUE.
       RETURN
    ENDIF

    ! trx must be present for radar obs
    IF (ades_loc_flg%rad .AND. .NOT. ades_loc_flg%trx) THEN
       error_msg='trx key not found '//err_strk
       WRITE(ierrou,*) 'ERROR: ',rtn,': ',error_msg
       error=.TRUE.
       RETURN
    ENDIF

    ! rcv must be present for radar obs
    IF (ades_loc_flg%rad .AND. .NOT. ades_loc_flg%rcv) THEN
       error_msg='rcv key not found '//err_strk
       WRITE(ierrou,*) 'ERROR: ',rtn,': ',error_msg
       error=.TRUE.
       RETURN
    ENDIF

    ! location group must be present for roving or satellite all optical obs
    ! but in ground-based obs
    location_flg = (/ades_loc_flg%pos(1),ades_loc_flg%pos(2),ades_loc_flg%pos(3),ades_loc_flg%sys,ades_loc_flg%ctr/)
    IF (.NOT. ades_loc_flg%rad .AND. .NOT. (ALL(location_flg .EQV. .TRUE.) .OR. ALL(location_flg .EQV. .FALSE.))) THEN
       error_msg='sys/ctr/pos keys must all be present or not present at all '//err_strk
       WRITE(ierrou,*) 'ERROR: ',rtn,': ',error_msg
       error=.TRUE.
       RETURN
    ENDIF
    ! check sys ctr and pos values
    ! for roving obs
    IF (.NOT. ades_loc_flg%rad .AND. ades_loc_flg%sys .AND. trim(adjustl(ades_loc_obs%stn)) .EQ. '247') THEN
       ! sys
       IF (ALL(sys_rov_values .NE. ades_loc_obs%sys)) THEN
          warn_msg='Forbidden sys value for roving observation '//warn_strv
          WRITE(*,*) 'WARNING: ',rtn,': ',warn_msg  
          value_found = .FALSE.
          WRITE(*,*) warn_skip
          RETURN
       ENDIF
       ! ctr
       IF (ades_loc_obs%ctr .NE. 399) THEN
          warn_msg='Forbidden ctr value for roving observation (only geocentric-ctr=399 allowed) '//warn_strv
          WRITE(*,*) 'WARNING: ',rtn,': ',warn_msg
          value_found = .FALSE.
          WRITE(*,*) warn_skip
          RETURN
       ENDIF
       ! check pos1,pos2 (long,lat) values
       IF (.NOT. ades_loc_flg%rad .AND. ANY(sys_rov_values .EQ. ades_loc_obs%sys) .AND. &
            & (ades_loc_obs%pos(1) .LT. 0.0 .OR. ades_loc_obs%pos(1) .GT. 360.0)) THEN
          warn_msg='Forbidden pos1 value '//warn_strv
          WRITE(*,*) 'WARNING: ',rtn,': ',warn_msg
          value_found = .FALSE.
          WRITE(*,*) warn_skip
          RETURN
       ENDIF
       IF (.NOT. ades_loc_flg%rad .AND. (trim(adjustl(ades_loc_obs%sys)) .EQ. "WGS84" .OR. trim(adjustl(ades_loc_obs%sys)) &
            & .EQ. "IAU") .AND. (ades_loc_obs%pos(2) .LT. -90.0 .OR. ades_loc_obs%pos(2) .GT. 90.0)) THEN
          warn_msg='Forbidden pos2 value '//warn_strv
          WRITE(*,*) 'WARNING: ',rtn,': ',warn_msg
          value_found = .FALSE.
          WRITE(*,*) warn_skip
          RETURN
       ENDIF
    ENDIF
    ! for satellite obs
    IF (.NOT. ades_loc_flg%rad .AND. ades_loc_flg%sys .AND. trim(adjustl(ades_loc_obs%stn)) .NE. '247') THEN
       ! sys
       IF (ALL(sys_sat_values .NE. ades_loc_obs%sys)) THEN
          warn_msg='Forbidden sys value for satellite observation '//warn_strv
          WRITE(*,*) 'WARNING: ',rtn,': ',warn_msg
          value_found = .FALSE.
          WRITE(*,*) warn_skip
          RETURN
       ENDIF
       ! ctr
       IF (ALL(ctr_values .NE. ades_loc_obs%ctr)) THEN
          warn_msg='Forbidden ctr value (only planets/Sun/Moon allowed) '//warn_strv
          WRITE(*,*) 'WARNING: ',rtn,': ',warn_msg
          value_found = .FALSE.
          WRITE(*,*) warn_skip
          RETURN
       ENDIF
    ENDIF

    ! obsTime must be always present
    IF (.NOT. ades_loc_flg%obsTime) THEN
       error_msg='obsTime key not found '//err_strk
       WRITE(ierrou,*) 'ERROR: ',rtn,': ',error_msg
       error=.TRUE.
       RETURN
       ! check obsTime values
    ELSE
       time_str = ades_loc_obs%obsTime
       CALL rmsp(time_str,len_str)
       READ(time_str(6:7),*)month
       READ(time_str(9:10),*)day
       READ(time_str(12:13),*)hour
       READ(time_str(15:16),*)minute
       READ(time_str(18:len_str-1),*)second
       IF ((month .LT. 1 .OR. month .GT. 12) .OR. (day .LT. 1 .OR. day .GT. 31) .OR. (hour .LT. 0 .OR. hour .GT. 24) &
            & .OR. (minute .LT. 0 .OR. minute .GT. 60) .OR. (second .LT. 0.d0 .OR. second .GT. 60.d0)) THEN
          warn_msg='Forbidden obsTime value '//warn_strv
          WRITE(*,*) 'WARNING: ',rtn,': ',warn_msg
          value_found = .FALSE.
          WRITE(*,*) warn_skip
          RETURN
       ENDIF
    ENDIF

    ! rmsTime must be positive
    IF (ades_loc_flg%rmsTime .AND. ades_loc_obs%rmsTime .LT. 0.d0) THEN
       warn_msg='Forbidden rmsTime value '//warn_strv
       WRITE(*,*) 'WARNING: ',rtn,': ',warn_msg
       !value_found = .FALSE.
       !WRITE(*,*) warn_skip
       !RETURN
    ENDIF

    ! ra and dec must be present for optical obs
    IF (ades_loc_flg%opt .AND. (.NOT. ades_loc_flg%ra .OR. .NOT. ades_loc_flg%dec)) THEN
       error_msg='ra and/or dec key not found '//err_strk
       WRITE(ierrou,*) 'ERROR: ',rtn,': ',error_msg
       error=.TRUE.
       RETURN
       ! check ra and dec values
    ELSEIF (ades_loc_flg%opt) THEN
       ! ra
       IF (ades_loc_obs%ra .LT. 0.d0 .OR. ades_loc_obs%ra .GE. 360.d0) THEN
          warn_msg='Forbidden ra value '//warn_strv
          WRITE(*,*) 'WARNING: ',rtn,': ',warn_msg
          value_found = .FALSE.
          WRITE(*,*) warn_skip
          RETURN
       ENDIF
       ! dec
       IF (ades_loc_obs%dec .LT. -90.d0 .OR. ades_loc_obs%dec .GT. 90.d0) THEN
          warn_msg='Forbidden dec value '//warn_strv
          WRITE(*,*) 'WARNING: ',rtn,': ',warn_msg
          value_found = .FALSE.
          WRITE(*,*) warn_skip
          RETURN
       ENDIF
    ENDIF

    ! raStar and decStar must be present for occultation obs
    IF (ades_loc_flg%occ .AND. (.NOT. ades_loc_flg%raStar .OR. .NOT. ades_loc_flg%decStar)) THEN
       error_msg='raStar and/or decStar key not found '//err_strk
       WRITE(ierrou,*) 'ERROR: ',rtn,': ',error_msg
       error=.TRUE.
       RETURN
       ! check raStar and decStar values
    ELSEIF (ades_loc_flg%occ) THEN
       ! raStar
       IF (ades_loc_obs%raStar .LT. 0.d0 .OR. ades_loc_obs%raStar .GE. 360.d0) THEN
          warn_msg='Forbidden raStar value '//warn_strv
          WRITE(*,*) 'WARNING: ',rtn,': ',warn_msg
          value_found = .FALSE.
          WRITE(*,*) warn_skip
          RETURN
       ENDIF
       ! decStar
       IF (ades_loc_obs%decStar .LT. -90.d0 .OR. ades_loc_obs%decStar .GT. 90.d0) THEN
          warn_msg='Forbidden decStar value '//warn_strv
          WRITE(*,*) 'WARNING: ',rtn,': ',warn_msg
          value_found = .FALSE.
          WRITE(*,*) warn_skip
          RETURN
       ENDIF
    ENDIF

    ! obsCenter must be present for offset obs
    IF (ades_loc_flg%off .AND. .NOT. ades_loc_flg%obsCenter) THEN
       error_msg='obsCenter key not found '//err_strk
       WRITE(ierrou,*) 'ERROR: ',rtn,': ',error_msg
       error=.TRUE.
       RETURN
    ENDIF

    ! if deltaRA is present for off/occ obs then deltaDec must be present as well and viceversa
    IF ((ades_loc_flg%off .OR. ades_loc_flg%occ) .AND. ((ades_loc_flg%deltaRA .AND. .NOT. ades_loc_flg%deltaDec) &
         & .OR. (.NOT. ades_loc_flg%deltaRA .AND. ades_loc_flg%deltaDec))) THEN
       error_msg='deltaRA/deltaDec must all be present or not present at all '//err_strk
       WRITE(ierrou,*) 'ERROR: ',rtn,': ',error_msg
       error=.TRUE.
       RETURN
    ENDIF
    ! deltaRA/deltaDec OR dist/pa must be present for off/occ obs
    IF ((ades_loc_flg%off .OR. ades_loc_flg%occ) .AND. (ades_loc_flg%deltaRA .AND. ades_loc_flg%dist)) THEN
       error_msg='both deltaRA/deltaDec and dist/pa keys found '//err_strk
       WRITE(ierrou,*) 'ERROR: ',rtn,': ',error_msg
       error=.TRUE.
       RETURN
    ENDIF

    ! if dist is present for off/occ obs then pa must be present as well and viceversa
    IF ((ades_loc_flg%off .OR. ades_loc_flg%occ) .AND. ((ades_loc_flg%dist .AND. .NOT. ades_loc_flg%pa) &
         & .OR. (.NOT. ades_loc_flg%dist .AND. ades_loc_flg%pa))) THEN
       error_msg='dist/pa must all be present or not present at all '//err_strk
       WRITE(ierrou,*) 'ERROR: ',rtn,': ',error_msg
       error=.TRUE.
       RETURN
    ELSEIF (ades_loc_flg%dist .AND. (ades_loc_flg%off .OR. ades_loc_flg%occ)) THEN
       ! check dist and pa values
       IF (ades_loc_flg%dist .AND. ades_loc_obs%dist .LT. 0.d0) THEN
          warn_msg='Forbidden dist value '//warn_strv
          WRITE(*,*) 'WARNING: ',rtn,': ',warn_msg
          value_found = .FALSE.
          WRITE(*,*) warn_skip
          RETURN
       ENDIF
       IF (ades_loc_flg%pa .AND. (ades_loc_obs%pa .LT. 0.d0 .OR. ades_loc_obs%pa .GE. 360.d0)) THEN
          warn_msg='Forbidden pa value '//warn_strv
          WRITE(*,*) 'WARNING: ',rtn,': ',warn_msg       
          value_found = .FALSE.
          WRITE(*,*) warn_skip
          RETURN
       ENDIF
    ENDIF

    ! if rmsRA is present for all optical obs then rmsDec must be present as well and viceversa
    IF (.NOT. ades_loc_flg%rad .AND. (ades_loc_flg%rmsRA .AND. .NOT. ades_loc_flg%rmsDec) .OR. &
         & (.NOT. ades_loc_flg%rmsRA .AND. ades_loc_flg%rmsDec)) THEN
       error_msg='rmsRA/rmsDec must all be present or not present at all '//err_strk
       WRITE(ierrou,*) 'ERROR: ',rtn,': ',error_msg
       error=.TRUE.
       RETURN
    ENDIF
    ! check rmsRA and rmsDec values
    IF (.NOT. ades_loc_flg%rad .AND. ades_loc_flg%rmsRA .AND. ades_loc_obs%rmsRA .LT. 0.d0) THEN
       warn_msg='Forbidden rmsRA value '//warn_strv
       WRITE(*,*) 'WARNING: ',rtn,': ',warn_msg
       !value_found = .FALSE.
       !RETURN
    ENDIF
    IF (.NOT. ades_loc_flg%rad .AND. ades_loc_flg%rmsDec .AND. ades_loc_obs%rmsDec .LT. 0.d0) THEN
       warn_msg='Forbidden rmsDec value '//warn_strv
       WRITE(*,*) 'WARNING: ',rtn,': ',warn_msg
       !value_found = .FALSE.
       !RETURN
    ENDIF

    ! if rmsDist is present for off/occ obs then rmsPA must be present as well and viceversa
    IF ((ades_loc_flg%off .OR. ades_loc_flg%occ) .AND. (ades_loc_flg%rmsDist .AND. .NOT. ades_loc_flg%rmsPA) .OR. &
         & (.NOT. ades_loc_flg%rmsDist .AND. ades_loc_flg%rmsPA)) THEN
       error_msg='rmsDist/rmsPA must all be present or not present at all '//err_strk
       WRITE(ierrou,*) 'ERROR: ',rtn,': ',error_msg
       error=.TRUE.
       RETURN
    ENDIF
    ! check rmsDist and rmsPA values
    IF ((ades_loc_flg%off .OR. ades_loc_flg%occ) .AND. ades_loc_flg%rmsDist .AND. ades_loc_obs%rmsDist .LT. 0.d0) THEN
       warn_msg='Forbidden rmsDist value '//warn_strv
       WRITE(*,*) 'WARNING: ',rtn,': ',warn_msg
       !value_found = .FALSE.
       !RETURN
    ENDIF
    IF ((ades_loc_flg%off .OR. ades_loc_flg%occ) .AND. ades_loc_flg%rmsPA .AND. ades_loc_obs%rmsPA .LT. 0.d0) THEN
       warn_msg='Forbidden rmsPA value '//warn_strv
       WRITE(*,*) 'WARNING: ',rtn,': ',warn_msg
       !value_found = .FALSE.
       !RETURN
    ENDIF

    ! check rmsCorr value
    IF (.NOT. ades_loc_flg%rad .AND. ades_loc_flg%rmsCorr .AND. &
         & (ades_loc_obs%rmsCorr .LT. -1.d0 .OR. ades_loc_obs%rmsCorr .GT. 1.d0)) THEN
       warn_msg='Forbidden rmsCorr value '//warn_strv
       WRITE(*,*) 'WARNING: ',rtn,': ',warn_msg
       !value_found = .FALSE.
       !RETURN
    ENDIF

    ! doppler OR delay must be present for radar obs, along with their rmsDoppler and rmsDelay
    IF (ades_loc_flg%rad .AND. ades_loc_flg%doppler .AND. ades_loc_flg%delay) THEN
       error_msg='both doppler and delay keys found '//err_strk
       WRITE(ierrou,*) 'ERROR: ',rtn,': ',error_msg
       error=.TRUE.
       RETURN
    ENDIF
    IF (ades_loc_flg%rad .AND. .NOT. ades_loc_flg%doppler .AND. .NOT. ades_loc_flg%delay) THEN
       error_msg='doppler or delay key not found '//err_strk
       WRITE(ierrou,*) 'ERROR: ',rtn,': ',error_msg
       error=.TRUE.
       RETURN
    ENDIF
    IF ((ades_loc_flg%doppler .AND. .NOT. ades_loc_flg%rmsDoppler) .OR. &
         &   (.NOT. ades_loc_flg%doppler .AND. ades_loc_flg%rmsDoppler)) THEN
       error_msg='doppler/rmsDoppler must all be present or not present at all '//err_strk
       WRITE(ierrou,*) 'ERROR: ',rtn,': ',error_msg
       error=.TRUE.
       RETURN
    ENDIF
    IF ((ades_loc_flg%delay .AND. .NOT. ades_loc_flg%rmsDelay) .OR. (.NOT. ades_loc_flg%delay .AND. ades_loc_flg%rmsDelay)) THEN
       error_msg='delay/rmsDelay must all be present or not present at all '//err_strk
       WRITE(ierrou,*) 'ERROR: ',rtn,': ',error_msg
       error=.TRUE.
       RETURN
    ENDIF
    ! check delay, rmsDelay, rmsDoppler values
    IF (ades_loc_flg%delay .AND. ades_loc_obs%delay .LT. 0.d0) THEN
       warn_msg='Forbidden delay value '//warn_strv
       WRITE(*,*) 'WARNING: ',rtn,': ',warn_msg
       value_found = .FALSE.
       WRITE(*,*) warn_skip
       RETURN
    ENDIF
    IF (ades_loc_flg%rad .AND. ades_loc_flg%rmsDelay .AND. ades_loc_obs%rmsDelay .LT. 0.d0) THEN
       warn_msg='Forbidden rmsDelay value '//warn_strv
       WRITE(*,*) 'WARNING: ',rtn,': ',warn_msg
       value_found = .FALSE.
       WRITE(*,*) warn_skip
       RETURN
    ENDIF
    IF (ades_loc_flg%rad .AND. ades_loc_flg%rmsDoppler .AND. ades_loc_obs%rmsDoppler .LT. 0.d0) THEN
       warn_msg='Forbidden rmsDoppler value '//warn_strv
       WRITE(*,*) 'WARNING: ',rtn,': ',warn_msg
       value_found = .FALSE.
       WRITE(*,*) warn_skip
       RETURN
    ENDIF

    ! astCat must be present for optical and occultation obs
    IF (.NOT. ades_loc_flg%rad .AND. .NOT.  ades_loc_flg%off .AND. .NOT. ades_loc_flg%astCat) THEN
       error_msg='astCat key not found '//err_strk
       WRITE(ierrou,*) 'ERROR: ',rtn,': ',error_msg
       error=.TRUE.
       RETURN
    ENDIF
    ! check astCat value
    ! if field is empty, obs is accepted but under-weighted
    ! if astCat value is a name that does not exist, set an error and return
    IF (.NOT. ades_loc_flg%rad .AND. .NOT.  ades_loc_flg%off .AND. ades_loc_flg%astCat .AND. &
         & ALL(astCat_values .NE. ades_loc_obs%astCat) .AND. ades_loc_obs%astCat .NE. "") THEN
       error_msg='Forbidden astCat value '//err_strk
       WRITE(ierrou,*) 'ERROR: ',rtn,': ',error_msg
       error=.TRUE.
       RETURN
    ENDIF

    ! if mag is present for optical obs then band must be present as well and viceversa
    IF (.NOT. ades_loc_flg%rad .AND. (ades_loc_flg%mag .AND. .NOT. ades_loc_flg%band) .OR. &
         &   (.NOT. ades_loc_flg%mag .AND. ades_loc_flg%band)) THEN
       error_msg='mag/band must all be present or not present at all '//err_strk
       WRITE(ierrou,*) 'ERROR: ',rtn,': ',error_msg
       error=.TRUE.
       RETURN
    ENDIF
    ! check rmsMag value
    IF (.NOT. ades_loc_flg%rad .AND. ades_loc_flg%rmsMag .AND. ades_loc_obs%rmsMag .LT. 0.d0) THEN
       warn_msg='Forbidden rmsMag value '//warn_strv
       WRITE(*,*) 'WARNING: ',rtn,': ',warn_msg
       !value_found = .FALSE.
       !RETURN
    ENDIF
    ! check band value
    IF (.NOT. ades_loc_flg%rad .AND. ades_loc_obs%mag .NE. 9.9d9 .AND. ades_loc_flg%band .AND. &
         & ALL(band_values .NE. ades_loc_obs%band)) THEN
       warn_msg='Forbidden band value '//warn_strv
       WRITE(*,*) 'WARNING: ',rtn,': ',warn_msg
       value_found = .FALSE.
       WRITE(*,*) warn_skip
       RETURN
    ENDIF
    ! check photCat value
    IF (.NOT. ades_loc_flg%rad .AND. ades_loc_obs%mag .NE. 9.9d9 .AND. ades_loc_flg%photCat .AND. &
         & ALL(photCat_values .NE. ades_loc_obs%photCat) .AND. ades_loc_obs%photCat .NE. "") THEN
       warn_msg='Forbidden photCat value '//warn_strv
       WRITE(*,*) 'WARNING: ',rtn,': ',warn_msg
       !value_found = .FALSE.
       !WRITE(*,*) warn_skip
       !RETURN
    ENDIF
    ! check photAp value
    IF (.NOT. ades_loc_flg%rad .AND. ades_loc_flg%photAp .AND. ades_loc_obs%photAp .LT. 0.d0) THEN
       warn_msg='Forbidden photAp value '//warn_strv
       WRITE(*,*) 'WARNING: ',rtn,': ',warn_msg
       !value_found = .FALSE.
       !RETURN
    ENDIF
    ! check nucMag value
    IF (.NOT. ades_loc_flg%rad .AND. ades_loc_flg%nucMag .AND. (ades_loc_obs%nucMag .NE. 0 .AND. ades_loc_obs%nucMag .NE. 1)) THEN
       warn_msg='Forbidden nucMag value '//warn_strv
       WRITE(*,*) 'WARNING: ',rtn,': ',warn_msg
       !value_found = .FALSE.
       !RETURN
    ENDIF

    ! check seeing value
    IF (.NOT. ades_loc_flg%rad .AND. ades_loc_flg%seeing .AND. ades_loc_obs%seeing .LT. 0.d0) THEN
       warn_msg='Forbidden seeing value '//warn_strv
       WRITE(*,*) 'WARNING: ',rtn,': ',warn_msg
       !value_found = .FALSE.
       !RETURN
    ENDIF
    ! check exp value
    IF ((ades_loc_flg%opt .OR. ades_loc_flg%off) .AND. ades_loc_flg%exp .AND. ades_loc_obs%exp .LT. 0.d0) THEN
       warn_msg='Forbidden exp value '//warn_strv
       WRITE(*,*) 'WARNING: ',rtn,': ',warn_msg
       !value_found = .FALSE.
       !RETURN
    ENDIF
    ! check rmsFit value
    IF ((ades_loc_flg%opt .OR. ades_loc_flg%off) .AND. ades_loc_flg%rmsFit .AND. ades_loc_obs%rmsFit .LT. 0.d0) THEN
       warn_msg='Forbidden rmsFit value '//warn_strv
       WRITE(*,*) 'WARNING: ',rtn,': ',warn_msg
       !value_found = .FALSE.
       !RETURN
    ENDIF
    ! check nStars value
    IF ((ades_loc_flg%opt .OR. ades_loc_flg%off) .AND. ades_loc_flg%nStars .AND. ades_loc_obs%nStars .LT. 0) THEN
       warn_msg='Forbidden nStars value '//warn_strv
       WRITE(*,*) 'WARNING: ',rtn,': ',warn_msg
       !value_found = .FALSE.
       !RETURN
    ENDIF
    ! check com value
    IF (ades_loc_flg%rad .AND. ades_loc_flg%com .AND. (ades_loc_obs%com .NE. 0 .AND. ades_loc_obs%com .NE. 1)) THEN
       warn_msg='Forbidden com value '//warn_strv
       WRITE(*,*) 'WARNING: ',rtn,': ',warn_msg
       value_found = .FALSE.
       WRITE(*,*) warn_skip
       RETURN
    ENDIF
    ! frq must be present for radar obs
    IF (ades_loc_flg%rad .AND. .NOT. ades_loc_flg%frq) THEN
       error_msg='frq key not found '//err_strk
       WRITE(ierrou,*) 'ERROR: ',rtn,': ',error_msg
       error=.TRUE.
       RETURN
    ENDIF
    ! check frq value
    IF (ades_loc_flg%rad .AND. ades_loc_flg%frq .AND. ades_loc_obs%frq .LT. 0.d0) THEN
       warn_msg='Forbidden frq value '//warn_strv
       WRITE(*,*) 'WARNING: ',rtn,': ',warn_msg
       value_found = .FALSE.
       WRITE(*,*) warn_skip
       RETURN
    ENDIF
    ! check disc value
    IF (.NOT. ades_loc_flg%rad .AND. ades_loc_flg%disc .AND. &
         &   (trim(adjustl(ades_loc_obs%disc)) .NE. "*" .AND. trim(adjustl(ades_loc_obs%disc)) .NE. "+" .AND. &
         &    trim(adjustl(ades_loc_obs%disc)) .NE. "")) THEN
       warn_msg='Forbidden disc value '//warn_strv
       WRITE(*,*) 'WARNING: ',rtn,': ',warn_msg
       !value_found = .FALSE.
       !WRITE(*,*) warn_skip
       !RETURN
    ENDIF

    ! if precision group is present, check that all the required fields are present
    opt_prec_flg = (/ades_loc_flg%precTime,ades_loc_flg%precRA,ades_loc_flg%precDec/)
    IF (.NOT. ades_loc_flg%rad .AND. .NOT. (ALL(opt_prec_flg .EQV. .TRUE.) .OR. ALL(opt_prec_flg .EQV. .FALSE.))) THEN
       error_msg='precTime,precRA,precDec must all be present or not present at all '//err_strk
       WRITE(ierrou,*) 'ERROR: ',rtn,': ',error_msg
       error=.TRUE.
       RETURN
    ENDIF
    ! check precTime value
    IF (.NOT. ades_loc_flg%rad .AND. ades_loc_flg%precTime .AND. ALL(precTime_values .NE. ades_loc_obs%precTime) &
         & .AND. ades_loc_obs%precTime .NE. 0) THEN
       warn_msg='Forbidden precTime value '//warn_strv
       WRITE(*,*) 'WARNING: ',rtn,': ',warn_msg
       !value_found = .FALSE.
       !RETURN
    ENDIF
    ! check precRA value
    IF (.NOT. ades_loc_flg%rad .AND. ades_loc_flg%precRA .AND. ALL(precRADec_values .NE. ades_loc_obs%precRA) &
         & .AND. ades_loc_obs%precRA .NE. 0.d0) THEN
       warn_msg='Forbidden precRA value '//warn_strv
       WRITE(*,*) 'WARNING: ',rtn,': ',warn_msg
       !value_found = .FALSE.
       !RETURN
    ENDIF
    ! check precDec value
    IF (.NOT. ades_loc_flg%rad .AND. ades_loc_flg%precDec .AND. ALL(precRADec_values .NE. ades_loc_obs%precDec) &
         & .AND. ades_loc_obs%precDec .NE. 0.d0) THEN
       warn_msg='Forbidden precDec value '//warn_strv
       WRITE(*,*) 'WARNING: ',rtn,': ',warn_msg
       !value_found = .FALSE.
       !RETURN
    ENDIF

    ! check uncTime value
    IF (.NOT. ades_loc_flg%rad .AND. ades_loc_flg%uncTime .AND. ades_loc_obs%uncTime .LT. 0.d0) THEN
       warn_msg='Forbidden uncTime value '//warn_strv
       WRITE(*,*) 'WARNING: ',rtn,': ',warn_msg
       !value_found = .FALSE.
       !RETURN
    ENDIF

    ! check notes values
    IF (.NOT. ades_loc_flg%rad .AND. ades_loc_flg%notes) THEN
       DO i=1,6
          IF (ALL(notes_values .NE. ades_loc_obs%notes(i:i)) .AND. ades_loc_obs%notes(i:i) .NE. " ") THEN
             warn_msg='Forbidden notes char value '//warn_strv
             WRITE(*,*) 'WARNING: ',rtn,': ',warn_msg
             !value_found = .FALSE.
             !RETURN
          ENDIF
       END DO
    ENDIF

    ! if residuals group is present, check that all the required fields are present
    ! optical astrometry residuals
    optastr_res_flg = (/ades_loc_flg%orbProd,ades_loc_flg%orbID,ades_loc_flg%resRA,ades_loc_flg%resDec,&
         &                 ades_loc_flg%selAst,ades_loc_flg%sigRA,ades_loc_flg%sigDec/)
    IF (.NOT. ades_loc_flg%rad .AND. .NOT. (ALL(optastr_res_flg .EQV. .TRUE.) .OR. ALL(optastr_res_flg .EQV. .FALSE.))) THEN
       error_msg='orbProd/orbID/resRA/resDec/selAst/sigRA/sigDec must all be present &
            &                 or not present at all '//err_strk
       WRITE(ierrou,*) 'ERROR: ',rtn,': ',error_msg
       error=.TRUE.
       RETURN
    ENDIF
    ! check selAst value
    IF (.NOT. ades_loc_flg%rad .AND. ades_loc_flg%selAst .AND. ades_loc_obs%selAst .NE. "A" &
         &   .AND. ades_loc_obs%selAst .NE. "a" .AND. ades_loc_obs%selAst .NE. "D" .AND. ades_loc_obs%selAst .NE. "d") THEN
       warn_msg='Forbidden selAst value '//warn_strv
       WRITE(*,*) 'WARNING: ',rtn,': ',warn_msg
       value_found = .FALSE.
       WRITE(*,*) warn_skip
       RETURN
    ENDIF
    ! check sigRA value
    IF (.NOT. ades_loc_flg%rad .AND. ades_loc_flg%sigRA .AND. ades_loc_obs%sigRA .LT. 0.d0) THEN
       warn_msg='Forbidden sigRA value '//warn_strv
       WRITE(*,*) 'WARNING: ',rtn,': ',warn_msg
       !value_found = .FALSE.
       !RETURN
    ENDIF
    ! check sigDec value
    IF (.NOT. ades_loc_flg%rad .AND. ades_loc_flg%sigDec .AND. ades_loc_obs%sigDec .LT. 0.d0) THEN
       warn_msg='Forbidden sigDec value '//warn_strv
       WRITE(*,*) 'WARNING: ',rtn,': ',warn_msg
       !value_found = .FALSE.
       !RETURN
    ENDIF
    ! check sigCorr value
    IF (.NOT. ades_loc_flg%rad .AND. ades_loc_flg%sigCorr .AND. &
         &   (ades_loc_obs%sigCorr .LT. -1.d0 .OR. ades_loc_obs%sigCorr .GT. 1.d0)) THEN
       warn_msg='Forbidden sigCorr value '//warn_strv
       WRITE(*,*) 'WARNING: ',rtn,': ',warn_msg
       !value_found = .FALSE.
       !RETURN
    ENDIF
    ! check sigTime value
    IF (.NOT. ades_loc_flg%rad .AND. ades_loc_flg%sigTime .AND. ades_loc_obs%sigTime .LT. 0.d0) THEN
       warn_msg='Forbidden sigTime value '//warn_strv
       WRITE(*,*) 'WARNING: ',rtn,': ',warn_msg
       !value_found = .FALSE.
       !RETURN
    ENDIF

    ! optical photometry residuals
    optphot_res_flg = (/ades_loc_flg%orbProd,ades_loc_flg%orbID,ades_loc_flg%resMag,ades_loc_flg%selPhot,ades_loc_flg%sigMag/)
    IF (.NOT. ades_loc_flg%rad .AND. ades_loc_flg%mag .AND. .NOT. (ALL(optphot_res_flg .EQV. .TRUE.) & 
         &   .OR. ALL(optphot_res_flg .EQV. .FALSE.))) THEN
       error_msg='orbProd/orbID/resMag/selPhot/sigMag must all be present or not present at all '//err_strk
       WRITE(ierrou,*) 'ERROR: ',rtn,': ',error_msg
       error=.TRUE.
       RETURN
    ENDIF
    ! check selPhot value
    IF (.NOT. ades_loc_flg%rad .AND. ades_loc_flg%selPhot .AND. ades_loc_obs%selPhot .NE. "A" &
         &   .AND. ades_loc_obs%selPhot .NE. "a" .AND. ades_loc_obs%selPhot .NE. "D" .AND. ades_loc_obs%selPhot .NE. "d") THEN
       warn_msg='Forbidden selPhot value '//warn_strv
       WRITE(*,*) 'WARNING: ',rtn,': ',warn_msg
       value_found = .FALSE.
       WRITE(*,*) warn_skip
       RETURN
    ENDIF
    ! check sigMag value
    IF (.NOT. ades_loc_flg%rad .AND. ades_loc_flg%sigMag .AND. ades_loc_obs%sigMag .LT. 0.d0) THEN
       warn_msg='Forbidden sigMag value '//warn_strv
       WRITE(*,*) 'WARNING: ',rtn,': ',warn_msg
       !value_found = .FALSE.
       !RETURN
    ENDIF
    ! check photMod value
    IF (.NOT. ades_loc_flg%rad .AND. ades_loc_obs%mag .NE. 9.9d9 .AND. ades_loc_flg%photMod &
         & .AND. ALL(photMod_values .NE. ades_loc_obs%photMod)) THEN
       warn_msg='Forbidden photMod value '//warn_strv
       WRITE(*,*) 'WARNING: ',rtn,': ',warn_msg
       !value_found = .FALSE.
       !RETURN
    ENDIF

    ! radar residuals
    rad_res_flg = (/ades_loc_flg%orbProd,ades_loc_flg%orbID,ades_loc_flg%resDelay,ades_loc_flg%selDelay,ades_loc_flg%sigDelay/)
    IF (ades_loc_flg%rad .AND. ades_loc_flg%delay .AND. &
         & .NOT. (ALL(rad_res_flg .EQV. .TRUE.) .OR. ALL(rad_res_flg .EQV. .FALSE.))) THEN
       error_msg='orbProd/orbID/resDelay/selDelay/sigDelay must all be present or not present at all '//err_strk
       WRITE(ierrou,*) 'ERROR: ',rtn,': ',error_msg
       error=.TRUE.
       RETURN
    ENDIF
    ! check selDelay value
    IF (ades_loc_flg%rad .AND. ades_loc_flg%selDelay .AND. ades_loc_obs%selDelay .NE. "A" .AND. ades_loc_obs%selDelay .NE. "a" &
         &   .AND. ades_loc_obs%selDelay .NE. "D" .AND. ades_loc_obs%selDelay .NE. "d") THEN
       warn_msg='Forbidden selDelay value '//warn_strv
       WRITE(*,*) 'WARNING: ',rtn,': ',warn_msg
       value_found = .FALSE.
       WRITE(*,*) warn_skip
       RETURN
    ENDIF
    ! check sigDelay value
    IF (ades_loc_flg%rad .AND. ades_loc_flg%sigDelay .AND. ades_loc_obs%sigDelay .LT. 0.d0) THEN
       warn_msg='Forbidden sigDelay value '//warn_strv
       WRITE(*,*) 'WARNING: ',rtn,': ',warn_msg
       !value_found = .FALSE.
       !RETURN
    ENDIF

    rad_res_flg = (/ades_loc_flg%orbProd,ades_loc_flg%orbID,ades_loc_flg%resDoppler,ades_loc_flg%selDoppler,&
         &               ades_loc_flg%sigDoppler/)
    IF (ades_loc_flg%rad .AND. ades_loc_flg%doppler .AND. & 
         & .NOT. (ALL(rad_res_flg .EQV. .TRUE.) .OR. ALL(rad_res_flg .EQV. .FALSE.))) THEN
       error_msg='orbProd/orbID/resDoppler/selDoppler/sigDoppler must all be present &
            &                 or not present at all '//err_strk
       WRITE(ierrou,*) 'ERROR: ',rtn,': ',error_msg
       error=.TRUE.
       RETURN
    ENDIF
    ! check selDoppler value
    IF (ades_loc_flg%rad .AND. ades_loc_flg%selDoppler .AND. ades_loc_obs%selDoppler .NE. "A" &
         & .AND. ades_loc_obs%selDoppler .NE. "a" .AND. ades_loc_obs%selDoppler .NE. "D" &
         & .AND. ades_loc_obs%selDoppler .NE. "d") THEN
       warn_msg='Forbidden selDoppler value '//warn_strv
       WRITE(*,*) 'WARNING: ',rtn,': ',warn_msg
       value_found = .FALSE.
       WRITE(*,*) warn_skip
       RETURN
    ENDIF
    ! check sigDoppler value
    IF (ades_loc_flg%rad .AND. ades_loc_flg%sigDoppler .AND. ades_loc_obs%sigDoppler .LT. 0.d0) THEN
       warn_msg='Forbidden sigDoppler value '//warn_strv
       WRITE(*,*) 'WARNING: ',rtn,': ',warn_msg
       !value_found =.FALSE.
       !RETURN
    ENDIF

    ! doppler OR delay residuals must be present for radar obs
    IF (ades_loc_flg%rad .AND. ades_loc_flg%resDoppler .AND. ades_loc_flg%resDelay) THEN
       error_msg='both resDoppler and resDelay keys found '//err_strk
       WRITE(ierrou,*) 'ERROR: ',rtn,': ',error_msg
       error=.TRUE.
       RETURN
    ENDIF

    ! check deprecated value
    IF (.NOT. ades_loc_flg%rad .AND. ades_loc_flg%deprecated .AND. (ades_loc_obs%deprecated .NE. "X" &
         & .AND. trim(adjustl(ades_loc_obs%deprecated)) .NE. "")) THEN
       warn_msg='Forbidden deprecated value '//warn_strv
       WRITE(*,*) 'WARNING: ',rtn,': ',warn_msg
       !value_found =.FALSE.
       !WRITE(*,*) warn_skip
       !RETURN
    ENDIF

  END SUBROUTINE check_fields


  ! WRITING output ades file

  ! ---------------------------------------------- !
  ! SUBROUTINE add_data_element                    !
  ! Add an element tag with text val to a parent   !
  ! node in ADES xml format. Uses functions and    !
  ! subroutines from FoX library.                  !
  ! ---------------------------------------------- !
  !  SUBROUTINE add_data_element(file,parent,tag,val)
  !    
  !    TYPE(Node), POINTER, INTENT(IN) :: file,parent 
  !    CHARACTER(len=*),    INTENT(IN) :: tag,val
  !    TYPE(Node), POINTER             :: temp
  !    
  !    temp => createElementNS(file,"",tag)  
  !    CALL setTextContent(temp,val)
  !    temp => appENDChild(parent,temp)
  !
  !  END SUBROUTINE add_data_element

  ! ------------------------------------------------ !
  ! FUNCTION add_element                             !
  ! Add an element tag to a parent node in ADES xml  !
  ! format, returning the added element.             !
  ! Uses functions and subroutines from FoX library. !
  ! ------------------------------------------------ !
  !  FUNCTION addelement(file,parent,tag) RESULT (el)
  !   
  !    TYPE(Node),      POINTER, INTENT(IN) :: file,parent
  !    CHARACTER(len=*),         INTENT(IN) :: tag
  !    TYPE(Node), POINTER :: temp,el
  !
  !    el => createElementNS(file,"",tag)
  !    temp => appENDChild(parent,el)
  !
  !  END FUNCTION addelement

  ! --------------------------------------------------- !
  ! SUBROUTINE write_xml                                !
  ! Writes a xml output ADES file using functions and   !
  ! subroutines from FoX library. Creates an xml tree   !
  ! in memory and writes it to a file.                  !
  ! --------------------------------------------------- !
  ! Lines of code related to FoX types, functions and   !
  ! subroutines are commented                           !

  SUBROUTINE  write_xml(filename,error_model,comment,RMSast,RMSmag)

    CHARACTER(LEN=*),   INTENT(IN)           :: filename                      ! name of output ADES file
    CHARACTER(LEN=*),   INTENT(IN)           :: error_model
    LOGICAL,            INTENT(IN)           :: comment                       ! true if the header contains the comment field
    REAL(KIND=dkind),   INTENT(IN)           :: RMSast                        ! global RMS of astrometric fit
    REAL(KIND=dkind),   INTENT(IN)           :: RMSmag                        ! global RMS of photometric fit

    CHARACTER(LEN=20)   :: rmsast_str,rmsmag_str,hstr1
    CHARACTER(LEN=100)  :: hstr2,hstr3,content,element,aux_hstr,hstr,aux_str
    LOGICAL             :: errmod_found,rms_found                             ! if false, we need to write Errormodel/RMS
    INTEGER             :: l,j,i,k,keys_idx,len_rec,tot_len,ele_len,cont_len
    !  TYPE(Node), POINTER :: ades_output,root,np,obsBlock,obsContext,obsData
    !  TYPE(Node), POINTER :: subContext,optical,radar,occultation,offset

    errmod_found = .FALSE.
    rms_found = .FALSE.
    aux_hstr = ""
    hstr = ""
    aux_str = ""
    element = ""
    content = ""
    !  subContext => null()
    !  optical => null()
    !  radar => null()
    !  occultation => null()
    !  offset => null()
    l = header_length    ! last line of header

    ! first, update ades_header with Errormodel and RMS
    ! store string hstr3 for writing rms
    IF (RMSast .GT. 0.d0 .OR. RMSmag .GT. 0.d0) THEN
       WRITE(rmsast_str,"(1P,E13.5)") RMSast
       WRITE(rmsmag_str,"(1P,E13.5)") RMSmag
       IF (RMSast .GT. 0.d0 .AND. RMSmag .GT. 0.d0) THEN
          hstr3 = "! line RMSast = "//trim(adjustl(rmsast_str))//" RMSmag = "//trim(adjustl(rmsmag_str))
       ELSEIF(RMSast .GT. 0.d0) THEN
          hstr3 = "! line RMSast = "//trim(adjustl(rmsast_str))
       ELSEIF(RMSmag .GT. 0.d0) THEN
          hstr3 = "! line RMSmag = "//trim(adjustl(rmsmag_str))
       ENDIF
    ELSE
       rms_found = .TRUE.    ! we do not need to write RMS in # comment
    END IF

    ! if there is not # comment in ades_header, need to add it at the end of the header
    IF (.NOT. comment) THEN
       ! add new lines to end of header
       hstr1 = "# comment"
       hstr2 = "! line Errormodel = "//trim(adjustl(error_model))
       ades_header(l+1) = trim(adjustl(hstr1))
       ades_header(l+2) = trim(adjustl(hstr2))
       ! add global rms only if at least one is greater than 0
       IF (RMSast .GT. 0.d0 .OR. RMSmag .GT. 0.d0) THEN
          ades_header(l+3) = trim(adjustl(hstr3))
          header_length = header_length + 3      ! update header_length
       ELSE
          header_length = header_length + 2
       ENDIF
       l = header_length
    ELSE ! if there is # comment in ades_header, replace or add Errormodel and rms
       DO j = 2,l
          ! if comment contains different Errormodel, replace it
          IF (INDEX(ades_header(j),"Errormodel") .GT. 0) THEN
             errmod_found = .TRUE.
             IF (ades_error_model.NE.trim(adjustl(error_model))) THEN
                ades_header(j) = "! line Errormodel = "//trim(adjustl(error_model))
             ENDIF
          ENDIF
          ! if comment contains RMSast and/or RMSmag, replace the line
          IF ((INDEX(ades_header(j),"RMSast") .GT. 0) .OR.(INDEX(ades_header(j),"RMSmag") .GT. 0)) THEN
             IF (RMSast .GT. 0.d0 .OR. RMSmag .GT. 0.d0) THEN
                ades_header(j) = trim(adjustl(hstr3))
                rms_found = .TRUE.
             ENDIF
          ENDIF
       ENDDO

       IF (.NOT.errmod_found) THEN  ! add Errormodel at the end of #comment (that is the end of header)
          hstr2 = "! line Errormodel = "//trim(adjustl(error_model))
          ades_header(l+1) = trim(adjustl(hstr2))
          l = l + 1
          header_length = l
       ENDIF

       IF (.NOT.rms_found) THEN  ! add RMS at the end of #comment (that is the end of header)
          ades_header(l+1) = trim(adjustl(hstr3))
          l = l + 1
          header_length = l
       ENDIF
    ENDIF

    ! now write xml tree    
    ! Create a new document and get a POINTER to the root element
    ! ades_output => createDocument(getImplementation(), "", "ades", null())  
    ! root => getDocumentElement(ades_output)

    ! Create new elements and append them to the document root element 
    ! obsBlock => addelement(ades_output,root,"obsBlock")
    ! obsContext => addelement(ades_output,root,"obsContext")     ! always write obsContext, although it is not given in input

    ! append elements to obsContext
    ! for each psv-format line of ades_header, create element or child element
    DO k = 2,l  ! l = header_length
       hstr = ades_header(k)
       IF (hstr(1:1) == '#') THEN   ! is an Element: in psv it is a row starting with #, in xml is a start-tag
          ele_len = len(trim(adjustl(hstr(2:))))
          element = trim(adjustl(hstr(2:)))
          !  subContext => addelement(ades_output,obsContext,element(1:ele_len))  
       ELSEIF (hstr(1:1) == '!')  THEN ! is a Child Element: in psv it is a row starting with !, in xml is an element (start-tag + content + end-tag)
          aux_hstr = trim(adjustl(hstr(2:)))
          i = INDEX(aux_hstr," ")    ! in psv format, space separates child element name from child element content
          element = trim(adjustl(aux_hstr(:i-1)))
          ele_len = len(trim(adjustl(aux_hstr(:i-1))))
          content = trim(adjustl(aux_hstr(i+1:)))
          cont_len = len(trim(adjustl(aux_hstr(i+1:))))
          ! IF (.NOT. ASSOCIATED(subContext)) THEN    ! control on consistency of element - child element structure
          ! write warning message
          ! WRITE(*,*) 'WARNING: write_xml: child elements without parent element in obsContext block &
          !      & of ADES file '//filename//'. Skipping line',k,'of header.'
          ! CYCLE
          !ENDIF
          !  CALL add_data_element(ades_output,subContext,element(1:ele_len),content(1:cont_len))   ! appending child element to parent element 
       ENDIF
    ENDDO

    !  obsData => addelement(ades_output, root, "obsData")
    DO j = 1,nobs_ades  ! loop over each observation (ades_obs(j))

       ! first, append observation element to obsData 
       IF(ades_flg(j)%opt) THEN
          !  optical => addelement(ades_output,obsData,"optical")
       ELSEIF(ades_flg(j)%rad) THEN
          !  radar => addelement(ades_output,obsData,"radar")
       ELSEIF(ades_flg(j)%occ) THEN
          !  occultation => addelement(ades_output,obsData,"occultation")
       ELSEIF(ades_flg(j)%off) THEN
          !  offset => addelement(ades_output,obsData,"offset")
       ENDIF

       ! then, append child elements
       ! use ades_data_rec and ades_key_rec that contain formatted data
       k = nobs_to_ntot_idx(j)  ! find data line number
       IF (new_key_rec(j)) keys_idx = k-1   ! find keys line number
       tot_len = MAX(len(trim(adjustl(ades_data_rec(k)))),len(trim(adjustl(ades_key_rec(keys_idx)))))
       i = INDEX(ades_data_rec(k),'|')
       len_rec = 0
       DO WHILE (i .GT. 0) 
          aux_str = ades_data_rec(k)(len_rec+1 : len_rec+i-1)
          content = trim(adjustl(aux_str))
          cont_len = len(trim(adjustl(aux_str)))
          aux_hstr = ades_key_rec(keys_idx)(len_rec+1 : len_rec+i-1)
          element = trim(adjustl(aux_hstr))
          ele_len = len(trim(adjustl(aux_hstr)))
          IF(ades_flg(j)%opt) THEN
             !  CALL add_data_element(ades_output,optical,element(1:ele_len),content(1:cont_len))
          ELSEIF(ades_flg(j)%rad) THEN
             !  CALL add_data_element(ades_output,radar,element(1:ele_len),content(1:cont_len))
          ELSEIF(ades_flg(j)%occ) THEN
             !  CALL add_data_element(ades_output,occultation,element(1:ele_len),content(1:cont_len))
          ELSEIF(ades_flg(j)%off) THEN
             !  CALL add_data_element(ades_output,offset,element(1:ele_len),content(1:cont_len))
          ENDIF
          len_rec = len_rec + i
          i = INDEX(ades_data_rec(k)(len_rec+1 :),'|')
       ENDDO
       ! add last key and value
       aux_str = ades_data_rec(k)(len_rec+1 : tot_len)
       content = trim(adjustl(aux_str))
       cont_len = len(trim(adjustl(aux_str)))
       aux_hstr = ades_key_rec(keys_idx)(len_rec+1 : tot_len)
       element = trim(adjustl(aux_hstr))
       ele_len = len(trim(adjustl(aux_hstr)))
       IF(ades_flg(j)%opt) THEN
          !  CALL add_data_element(ades_output,optical,element(1:ele_len),content(1:cont_len))
       ELSEIF(ades_flg(j)%rad) THEN
          !  CALL add_data_element(ades_output,radar,element(1:ele_len),content(1:cont_len))
       ELSEIF(ades_flg(j)%occ) THEN
          !  CALL add_data_element(ades_output,occultation,element(1:ele_len),content(1:cont_len))
       ELSEIF(ades_flg(j)%off) THEN
          !  CALL add_data_element(ades_output,offset,element(1:ele_len),content(1:cont_len))
       ENDIF

       ! add localUse, if present
       IF (ades_flg(j)%localUse) THEN
          aux_str = ades_obs(j)%localUse
          content = trim(adjustl(aux_str))
          cont_len = len(trim(adjustl(aux_str)))
          aux_hstr = "localUse"
          element = trim(adjustl(aux_hstr))
          ele_len = len(trim(adjustl(aux_hstr)))
          IF(ades_flg(j)%opt) THEN
             !  CALL add_data_element(ades_output,optical,element(1:ele_len),content(1:cont_len))
          ELSEIF(ades_flg(j)%rad) THEN
             !  CALL add_data_element(ades_output,radar,element(1:ele_len),content(1:cont_len))
          ELSEIF(ades_flg(j)%occ) THEN
             !  CALL add_data_element(ades_output,occultation,element(1:ele_len),content(1:cont_len))
          ELSEIF(ades_flg(j)%off) THEN
             !  CALL add_data_element(ades_output,offset,element(1:ele_len),content(1:cont_len))
          ENDIF
       ENDIF
    ENDDO

    ! routine storing xml tree on ades_output as standard xml file format
    !CALL setParameter(getDomConfig(ades_output), "invalid-pretty-print", .true.)

    ! CALL serialize to produce the document
    !CALL serialize(ades_output,filename)

    ! clean up the memory
    !CALL destroy(ades_output)

  END SUBROUTINE write_xml

  ! ---------------------------------------------------- !
  ! SUBROUTINE write_ades                                !
  ! It writes an output ADES file: .psv, .xml or both,   !
  ! according to the value of the input flag ades_type.  !
  ! ---------------------------------------------------- !
  ! Provisionally neglect the case of writing .xml files !

  SUBROUTINE write_ades(filename,error_model,ades_type,RMSast,RMSmag)

    CHARACTER(LEN=*),   INTENT(IN)           :: filename                ! name of output ADES file (without file extension)
    CHARACTER(LEN=*),   INTENT(IN)           :: error_model
    INTEGER,            INTENT(IN)           :: ades_type               ! flag for ADES output file (.psv/.xml)
    REAL(KIND=dkind),   INTENT(IN), OPTIONAL :: RMSast                  ! global RMS of astrometric fit
    REAL(KIND=dkind),   INTENT(IN), OPTIONAL :: RMSmag                  ! global RMS of photometric fit

    INTEGER             :: iline,unit_ades,k,i,j,ntot_ades_old,header_length_old,diff_idx
    REAL(KIND=dkind)    :: RMSast_loc,RMSmag_loc
    LOGICAL             :: comment   ! true if the header contains the comment field
    CHARACTER(LEN=200)  :: aux_hstr
    CHARACTER(LEN=20)   :: FMT,rmsast_str,rmsmag_str
    CHARACTER(LEN=100)  :: hstr

    CHARACTER(LEN=500), DIMENSION(nobx) :: ades_key_rec_tmp,ades_data_rec_tmp

    EXTERNAL            :: filopn,filclo

    ! get optional arguments    
    IF(PRESENT(RMSast)) THEN
       RMSast_loc = RMSast
    ELSE
       RMSast_loc = -1.d0
    ENDIF
    IF(PRESENT(RMSmag)) THEN
       RMSmag_loc = RMSmag
    ELSE
       RMSmag_loc = -1.d0
    ENDIF

    ! initialization
    comment = .FALSE.
    DO i = 1,header_length
       IF ((INDEX(ades_header(i),"# comment") .GT. 0) .OR. (INDEX(ades_header(i),"#comment") .GT. 0)) comment = .TRUE.
    END DO

    IF (.NOT. new_ades .AND. header_length .GT. 1) THEN 
       ! output ADES file is the same as the input ADES file, thus write again the ObsContext block

       !IF (ades_type == -1 .OR. ades_type == 1 .OR. ades_type == -3 .OR. ades_type == 3) THEN

       ! write records on .rwo_psv       
       ! open output .rwo_psv file
       CALL filopn(unit_ades,filename//'.rwo_psv','REPLACE')

       DO iline=1,ntot_ades
          ! write header segment
          aux_hstr = trim(ades_header(iline))
          IF (aux_hstr(1:1) .EQ. "#") THEN
             CALL write_header_ades(unit_ades,iline,RMSast_loc,RMSmag_loc,comment,error_model)
             ! write obsData segments
             ! if key record write the corresponding segment
          ELSEIF (trim(adjustl(ades_key_rec(iline))) .NE. "") THEN
             CALL write_segment_ades(iline,unit_ades)
          ENDIF
       END DO
       
       ! close file
       CALL filclo(unit_ades,' ')

       ! write .xml
       !IF (ades_type == -3 .OR. ades_type == 3) CALL write_xml(filename//'.xml',error_model,comment,RMSast_loc,RMSmag_loc)        
       !ELSEIF (ades_type == -2 .OR. ades_type == 2) THEN
       ! store remarks and Residuals group in ades_data_rec lines
       !DO iline=1,ntot_ades    
       !   IF (trim(adjustl(ades_key_rec(iline))) .NE. "") CALL write_segment_ades(iline)
       !END DO
       ! write .xml
       !CALL write_xml(filename//'.xml',error_model,comment,RMSast_loc,RMSmag_loc)       
       !ENDIF

    ELSE  ! there is no ADES input file or, if it exists, data from it have been changed
       ! (thus the output file cannot have the same ObsContext block from input)
       ! initialize header 
       ades_header = ""

       ! store first line and comment lines in ades_header
       header_length_old = header_length
       ades_header(1) = "# version=2017"
       ades_header(2) = "# comment"
       ades_header(3) = "! line Errormodel = "//trim(adjustl(error_model))
       comment = .TRUE.
       header_length = 3
          
       !IF (ades_type == -1 .OR. ades_type == 1 .OR. ades_type == -3 .OR. ades_type == 3) THEN
       ! write .psv
       ! open output .psv file
       CALL filopn(unit_ades,filename//'.rwo_psv','REPLACE')

       ! write ades_header
       DO k = 1,3
          WRITE(FMT,"(i20)") len(trim(adjustl(ades_header(k)))) 
          WRITE(unit_ades,"(A"//trim(adjustl(FMT))//")") trim(adjustl(ades_header(k)))
       END DO

       ! write RMS, if present
       IF (RMSast_loc .GT. 0.d0 .OR. RMSmag_loc .GT. 0.d0) THEN
          WRITE(rmsast_str,"(1P,E13.5)") RMSast_loc
          WRITE(rmsmag_str,"(1P,E13.5)") RMSmag_loc
          ! add rms 
          IF (RMSast_loc .GT. 0.d0 .AND. RMSmag_loc .GT. 0.d0) THEN
             hstr = "! line RMSast = "//trim(adjustl(rmsast_str))//" RMSmag = "//trim(adjustl(rmsmag_str))
          ELSEIF(RMSast_loc .GT. 0.d0) THEN
             hstr = "! line RMSast = "//trim(adjustl(rmsast_str))
          ELSEIF(RMSmag_loc .GT. 0.d0) THEN
             hstr = "! line RMSmag = "//trim(adjustl(rmsmag_str))
          ENDIF
          header_length = header_length + 1
          ! write new line to header
          WRITE(FMT,"(i20)") len(trim(adjustl(hstr)))
          WRITE(unit_ades,"(A"//trim(adjustl(FMT))//")") trim(adjustl(hstr))
       ENDIF

       ! written a new header, thus the total number of lines is changed
       diff_idx = header_length - header_length_old
       ntot_ades_old = ntot_ades
       ntot_ades = ntot_ades + diff_idx
       ! re-arranging ADES arrays related to lines of text
       ! observations number is unchanged
       DO j = 1,nobs_ades
          nobs_to_ntot_idx(j) = nobs_to_ntot_idx(j) + diff_idx
          ntot_to_nobs_idx(nobs_to_ntot_idx(j)) = j
       END DO
       ! shift data and keys lines
       ades_data_rec_tmp = ""
       ades_key_rec_tmp = ""
       ades_data_rec_tmp(header_length+1:ntot_ades) = ades_data_rec(header_length_old+1:ntot_ades_old)
       ades_key_rec_tmp(header_length+1:ntot_ades) = ades_key_rec(header_length_old+1:ntot_ades_old)
       ades_data_rec = ades_data_rec_tmp
       ades_key_rec = ades_key_rec_tmp
       
       ! write obsData segments
       DO iline=1,ntot_ades
          IF (trim(adjustl(ades_key_rec(iline))) .NE. "") CALL write_segment_ades(iline,unit_ades)
       END DO

       ! close file
       CALL filclo(unit_ades,' ')

       ! write .xml
       ! IF (ades_type == -3 .OR. ades_type == 3) CALL write_xml(filename//'.xml',error_model,comment,RMSast_loc,RMSmag_loc)

       !ELSEIF (ades_type == -2 .OR. ades_type == 2) THEN
       ! store remarks and Residuals group in ades_data_rec lines
       ! DO iline=1,ntot_ades    
       !    IF (trim(adjustl(ades_key_rec(iline))) .NE. "") CALL write_segment_ades(iline)
       ! END DO

       ! write .xml
       ! CALL write_xml(filename//'.xml',error_model,comment,RMSast_loc,RMSmag_loc)

       !ENDIF
    ENDIF

  END SUBROUTINE write_ades


  ! -------------------------------------------------- !
  ! SUBROUTINE write_header_ades                       !
  ! It writes the header (segment) of ADES psv file.   !
  ! -------------------------------------------------- !

  SUBROUTINE write_header_ades(unit,i,rms_ast,rms_mag,comment,error_model)

    INTEGER,            INTENT(IN)    :: unit            ! ADES file unit
    INTEGER,            INTENT(IN)    :: i               ! line number
    REAL(KIND=dkind),   INTENT(IN)    :: rms_ast         ! global RMS of astrometric fit
    REAL(KIND=dkind),   INTENT(IN)    :: rms_mag         ! global RMS of photometric fit
    LOGICAL,            INTENT(INOUT) :: comment         ! if FALSE, there is no # comment field in input ADES file ObsContext field
    CHARACTER(LEN=*),   INTENT(IN)    :: error_model

    CHARACTER(LEN=20)  :: FMT,rmsast_str,rmsmag_str,hstr1
    CHARACTER(LEN=100) :: hstr2,hstr3
    CHARACTER(LEN=200) :: aux_hstr
    INTEGER :: j,k

    ! if first line set the version and write
    IF (i .EQ. 1) THEN

       ades_header(i) = "# version=2017"
       WRITE(FMT,"(i20)") len(trim(adjustl(ades_header(i)))) 
       WRITE(unit,"(A"//trim(adjustl(FMT))//")") trim(adjustl(ades_header(i))) 

       ! if line is not first line and not related to comment write the segment
    ELSEIF (trim(adjustl(ades_header(i))) .NE. "# comment" .AND. trim(adjustl(ades_header(i))) .NE. "#comment") THEN

       ! write current line
       WRITE(FMT,"(i20)") len(trim(adjustl(ades_header(i))))
       WRITE(unit,"(A"//trim(adjustl(FMT))//")") trim(adjustl(ades_header(i)))

       ! write other lines belonging to segment
       j = i+1
       aux_hstr = trim(ades_header(j))
       DO WHILE (aux_hstr(1:1) .EQ. "!")
          WRITE(FMT,"(i20)") len(trim(adjustl(ades_header(j))))
          WRITE(unit,"(A"//trim(adjustl(FMT))//")") trim(adjustl(ades_header(j)))
          j = j+1
          aux_hstr = trim(ades_header(j))
       END DO

       ! if there is not # comment in ades_header, need to add it at the end of the header
       IF (.NOT. comment) THEN
          ! check if this is the last line of first header using j computed above
          IF (trim(adjustl(ades_header(j))) .EQ. "") THEN
             ! add new lines to end of header
             hstr1 = "# comment"
             hstr2 = "! line Errormodel = "//trim(adjustl(error_model))
             WRITE(FMT,"(i20)") len(trim(adjustl(hstr1)))
             WRITE(unit,"(A"//trim(adjustl(FMT))//")") trim(adjustl(hstr1))
             WRITE(FMT,"(i20)") len(trim(adjustl(hstr2)))
             WRITE(unit,"(A"//trim(adjustl(FMT))//")") trim(adjustl(hstr2))

             ! add rms
             IF (rms_ast .GT. 0.d0 .OR. rms_mag .GT. 0.d0) THEN
                WRITE(rmsast_str,"(1P,E13.5)") rms_ast
                WRITE(rmsmag_str,"(1P,E13.5)") rms_mag
                IF (rms_ast .GT. 0.d0 .AND. rms_mag .GT. 0.d0) THEN
                   hstr3 = "! line RMSast = "//trim(adjustl(rmsast_str))//" RMSmag = "//trim(adjustl(rmsmag_str))
                ELSEIF(rms_ast .GT. 0.d0) THEN
                   hstr3 = "! line RMSast = "//trim(adjustl(rmsast_str))
                ELSEIF(rms_mag .GT. 0.d0) THEN
                   hstr3 = "! line RMSmag = "//trim(adjustl(rmsmag_str))
                ENDIF
             ENDIF
          ENDIF
       ENDIF

       ! if line is comment and not first line
    ELSE

       ! write current line
       WRITE(FMT,"(i20)") len(trim(adjustl(ades_header(i))))
       WRITE(unit,"(A"//trim(adjustl(FMT))//")") trim(adjustl(ades_header(i)))
       ! update line index
       j = i+1   

       ! if next line contains Errormodel, replace it and write the line; otherwise, only write the line
       IF (INDEX(ades_header(i+1),"Errormodel") .GT. 0) THEN
          IF (ades_error_model .NE. trim(adjustl(error_model))) &
               & ades_header(i+1) = "! line Errormodel = "//trim(adjustl(error_model))
          WRITE(FMT,"(i20)") len(trim(adjustl(ades_header(i+1))))
          WRITE(unit,"(A"//trim(adjustl(FMT))//")") trim(adjustl(ades_header(i+1)))
          k = i + 2
          j = k
       ELSE
          hstr2 = "! line Errormodel = "//trim(adjustl(error_model))
          WRITE(FMT,"(i20)") len(trim(adjustl(hstr2)))
          WRITE(unit,"(A"//trim(adjustl(FMT))//")") trim(adjustl(hstr2))
          k = i + 1
       ENDIF

       ! if we need to write new RMSast and RMSmag then replace or write next line
       IF (rms_ast .GT. 0.d0 .OR. rms_mag .GT. 0.d0) THEN
          WRITE(rmsast_str,"(1P,E13.5)") rms_ast
          WRITE(rmsmag_str,"(1P,E13.5)") rms_mag

          IF (rms_ast .GT. 0.d0 .AND. rms_mag .GT. 0.d0) THEN
             hstr3 = "! line RMSast = "//trim(adjustl(rmsast_str))//" RMSmag = "//trim(adjustl(rmsmag_str))
          ELSEIF(rms_ast .GT. 0.d0) THEN
             hstr3 = "! line RMSast = "//trim(adjustl(rmsast_str))
          ELSEIF(rms_mag .GT. 0.d0) THEN
             hstr3 = "! line RMSmag = "//trim(adjustl(rmsmag_str))
          ENDIF
          ! if next line contains RMSast and RMSmag, replace it and write the line; otherwise, only write the line
          IF ((INDEX(ades_header(k),"RMSast") .GT. 0) .OR. (INDEX(ades_header(k),"RMSmag") .GT. 0)) THEN
             ades_header(k) = hstr3
             WRITE(FMT,"(i20)") len(trim(adjustl(ades_header(k))))
             WRITE(unit,"(A"//trim(adjustl(FMT))//")") trim(adjustl(ades_header(k)))
             j = j+1
          ELSE
             WRITE(FMT,"(i20)") len(trim(adjustl(hstr3)))
             WRITE(unit,"(A"//trim(adjustl(FMT))//")") trim(adjustl(hstr3))
          ENDIF

       ENDIF

       ! write other lines belonging to segment
       aux_hstr = trim(ades_header(j))
       DO WHILE (aux_hstr(1:1) .EQ. "!")
          WRITE(FMT,"(i20)") len(trim(adjustl(ades_header(j))))
          WRITE(unit,"(A"//trim(adjustl(FMT))//")") trim(adjustl(ades_header(j)))
          j = j+1
          aux_hstr = trim(ades_header(j))
       ENDDO
    ENDIF

  END SUBROUTINE write_header_ades

  ! ----------------------------------------------------- !
  ! SUBROUTINE write_segment_ades                         !
  ! It writes one segment (key + corresponding data) to   !
  ! ADES psv file.                                        !
  ! ----------------------------------------------------- !

  SUBROUTINE  write_segment_ades(i,unit)

    INTEGER,   INTENT(IN)            :: i                             ! line number
    INTEGER,   INTENT(IN), OPTIONAL  :: unit                          ! ADES file unit

    CHARACTER(LEN=600)                            :: key_res_group    ! residuals group string to be added to key record
    CHARACTER(LEN=600), DIMENSION(:), ALLOCATABLE :: data_res_group   ! residuals group string array to be added to data record

    CHARACTER(LEN=20) :: FMT
    CHARACTER(LEN=7) :: line
    CHARACTER(LEN=600) :: aux_kstr,aux_dstr
    INTEGER :: ndata_segment,idx_obs,j,k,idx_start,write_start_idx,len_rec,npipes_right,count_pipes
    LOGICAL :: pipe_found,skipped_segment,deprecated_found,remarks_found
    INTEGER :: status

    ! check that the key record has at least one corresponding data record
    IF (i .LT. ntot_ades .AND. trim(ades_data_rec(i+1)) .EQ. "") THEN
       WRITE(line,'(i7)')i
       WRITE(ierrou,*) 'write_ades_segment: key record does not have any corresponding data record &
            &for input ADES file: '//trim(ades_filename)//' line: '//trim(adjustl(line))
       numerr=numerr+1 
       RETURN
    ENDIF

    ! check if deprecated and remarks keys are present
    deprecated_found = .FALSE.
    IF (INDEX(ades_key_rec(i),"deprecated") .GT. 0) deprecated_found = .TRUE.
    remarks_found = .FALSE.
    IF (INDEX(ades_key_rec(i),"remarks") .GT. 0) remarks_found = .TRUE.

    ! if key record contains residuals group then remove the group from the key record and from the corresponding data records
    ! we apply this cleaning to all lines, including offset and deprecated observations
    CALL clean_residuals_group(i)

    ! if remarks entry in data record contains chi then remove this chi from data records belonging to the segment
    ! and do the same for the corresponding element in ades_obs
    ! we apply this cleaning to all lines, including offset and deprecated observations

    CALL clean_remarks(i)

    ! 1) check whether segment has at least one non-skipped obs
    ! 2) find segment's output writing start index
    ! 3) find number of data records in segment
    skipped_segment = .TRUE.
    write_start_idx = 0
    len_rec = 0
    ndata_segment = 0
    j = i+1
    DO WHILE (ades_data_rec(j) .NE. "")

       idx_obs = ntot_to_nobs_idx(j)
       ! if non-skipped line
       IF (idx_obs .GT. 0) THEN

          skipped_segment = .FALSE.

          ! find writing starting index using remarks field in data records
          ! in this case remarks corresponds to the last pipe of the line
          IF (ades_flg(idx_obs)%remarks .AND. trim(adjustl(ades_obs(idx_obs)%remarks)) .NE. "") THEN
             idx_start = len(ades_data_rec(j))
             pipe_found = .FALSE.
             DO WHILE (.NOT. pipe_found)
                IF (ades_data_rec(j)(idx_start:idx_start) .EQ. "|") pipe_found = .TRUE.
                idx_start = idx_start - 1
             END DO
             write_start_idx = MAX(write_start_idx,idx_start-1)
          ENDIF

          len_rec = MAX(len_rec,len(trim(ades_data_rec(j))))

       ENDIF

       ! get total number of data records in segment
       ndata_segment = ndata_segment +1
       j = j+1       
    END DO

    ! find writing starting index using remarks field in key records
    IF ((.NOT. skipped_segment) .AND. remarks_found) THEN
       idx_start = INDEX(ades_key_rec(i),"remarks")
       pipe_found = .FALSE.
       DO WHILE (.NOT. pipe_found)
          idx_start = idx_start - 1
          IF (ades_key_rec(i)(idx_start:idx_start) .EQ. "|") pipe_found = .TRUE.
       END DO
       write_start_idx = MAX(write_start_idx,idx_start-1)
    ENDIF

    ! if deprecated is present change write starting index
    IF ((.NOT. skipped_segment) .AND. deprecated_found) THEN
       ! data records
       write_start_idx = 0
       DO j=i+1,i+ndata_segment 
          idx_start = len(ades_data_rec(j))
          ! remarks are not removed for deprecated data lines, so we need to count two pipes
          IF (remarks_found) THEN
             npipes_right = 2
             ! otherwise deprecated corresponds to the last pipe
          ELSE
             npipes_right = 1
          ENDIF
          count_pipes = 0
          DO WHILE (count_pipes .LT. npipes_right)
             IF (ades_data_rec(j)(idx_start:idx_start) .EQ. "|") count_pipes = count_pipes + 1
             idx_start = idx_start - 1
          END DO
          write_start_idx = MAX(write_start_idx,idx_start-1)
       ENDDO
       ! key record
       idx_start = INDEX(ades_key_rec(i),"deprecated")
       pipe_found = .FALSE.
       DO WHILE (.NOT. pipe_found)
          idx_start = idx_start - 1
          IF (ades_key_rec(i)(idx_start:idx_start) .EQ. "|") pipe_found = .TRUE.
       END DO
       write_start_idx = MAX(write_start_idx,idx_start-1)
    ENDIF

    len_rec = MAX(len_rec,len(trim(ades_key_rec(i))))

    ! if no deprecated and/or remarks entry has been found, set segment's output writing start index to end of segment's longest line
    IF (write_start_idx .EQ. 0) write_start_idx = len_rec

    ! allocate array
    IF(ALLOCATED(data_res_group)) DEALLOCATE(data_res_group)
    ALLOCATE(data_res_group(ndata_segment),STAT=status)
    IF(status.GT.0)THEN
       WRITE(ierrou,*) 'write_ades_segment: insufficient memory to allocate data_res_group'
       numerr=numerr+1
    END IF
    ! get residuals + remarks formatted entries
    IF (.NOT. skipped_segment) CALL get_residuals_group(i,ndata_segment,key_res_group,data_res_group)

    ! loop over next data records belonging to the same segment defined by ades_key_rec(i)
    DO j=i+1,i+ndata_segment

       ! get index of observations
       idx_obs = ntot_to_nobs_idx(j)

       ! write key record if we are in first loop
       IF (j .EQ. i+1) THEN

          ! we add residuals group and remarks at the end of the line
          IF (.NOT. skipped_segment) ades_key_rec(i) = ades_key_rec(i)(1:write_start_idx)//key_res_group

          IF (PRESENT(unit)) THEN   
             ! write key record to psv file
             WRITE(FMT,"(i20)") len(trim(adjustl(ades_key_rec(i))))
             WRITE(unit,"(A"//trim(adjustl(FMT))//")") trim(adjustl(ades_key_rec(i)))
          ENDIF
       ENDIF

       ! add residuals group (depending whether the segment has been skipped)

       ! we add residuals group and remarks at the end of the line
       IF (.NOT. skipped_segment) ades_data_rec(j) = ades_data_rec(j)(1:write_start_idx)//data_res_group(j-i)

       IF (PRESENT(unit)) THEN
          ! write data record to psv file
          WRITE(FMT,"(i20)") len(trim(adjustl(ades_data_rec(j))))
          WRITE(unit,"(A"//trim(adjustl(FMT))//")") trim(adjustl(ades_data_rec(j)))
       ENDIF
    END DO

    ! deallocate array
    DEALLOCATE(data_res_group)

  END SUBROUTINE  write_segment_ades

  ! ----------------------------------------------------- !
  ! SUBROUTINE clean_residuals_group                      !
  ! For a given key record it removes the residuals group !
  ! from the key record and from the data records         !
  ! belonging to that specific key record.                !
  ! ----------------------------------------------------- !

  SUBROUTINE clean_residuals_group(i)

    INTEGER,            INTENT(IN)    :: i               ! line number of key record

    CHARACTER(LEN=7) :: line
    LOGICAL :: pipe_found,resgroup_found,deprecated_found,remarks_found
    INTEGER :: j,idx_start,idx_stop,len_str
    INTEGER :: npipes_resgroup,npipes_right,count_pipes

    ! initialization
    idx_start = 0
    idx_stop = 0
    resgroup_found = .FALSE.
    IF (INDEX(ades_key_rec(i),"orbProd") .GT. 0) resgroup_found = .TRUE.
    deprecated_found = .FALSE.
    IF (INDEX(ades_key_rec(i),"deprecated") .GT. 0) deprecated_found = .TRUE.
    remarks_found = .FALSE.
    IF (INDEX(ades_key_rec(i),"remarks") .GT. 0) remarks_found = .TRUE.

    ! if residuals group is not found then return
    IF (.NOT. resgroup_found) RETURN

    ! find start and stop indexes of residuals group
    ! length of string without leading and trailing spaces
    len_str = len(trim(adjustl(ades_key_rec(i))))

    ! find starting index
    idx_start = INDEX(ades_key_rec(i),"orbProd")
    ! find closest pipe towards left
    pipe_found = .FALSE.
    DO WHILE (.NOT. pipe_found)
       idx_start = idx_start - 1
       IF (ades_key_rec(i)(idx_start:idx_start) .EQ. "|") pipe_found = .TRUE.
    END DO

    ! find stop index
    ! if deprecated is present we need to clean up to that entry
    IF (deprecated_found) THEN
       idx_stop = INDEX(ades_key_rec(i),"deprecated")
       pipe_found = .FALSE.
       DO WHILE (.NOT. pipe_found)
          idx_stop = idx_stop - 1
          IF (ades_key_rec(i)(idx_stop:idx_stop) .EQ. "|") pipe_found = .TRUE.
       END DO
       ! if deprecated is missing but remarks is present we need to clean up to that entry
    ELSEIF (remarks_found) THEN
       idx_stop = INDEX(ades_key_rec(i),"remarks")
       pipe_found = .FALSE.
       DO WHILE (.NOT. pipe_found)
          idx_stop = idx_stop - 1
          IF (ades_key_rec(i)(idx_stop:idx_stop) .EQ. "|") pipe_found = .TRUE.
       END DO
       ! otherwise if remarks is missing it means that residuals group is at the end of the line, therefore clean up to the EOL
    ELSE
       idx_stop = len_str
    ENDIF

    ! if indexes are zero it means something went wrong
    IF (idx_start .EQ. 0 .OR. idx_stop .EQ. 0)THEN
       WRITE(line,'(i7)')i
       WRITE(ierrou,*) 'clean_ades_segment_resgroup: unable to find pipe corresponding to &
            & key record residuals group for input ADES file: '//trim(ades_filename)// ' line: '//trim(adjustl(line))
       numerr=numerr+1
       RETURN
    ENDIF

    ! count number of pipes (entries) of residuals group
    npipes_resgroup = 1 ! we already consider the one at idx_start
    DO j=idx_start+1,idx_stop
       IF (ades_key_rec(i)(j:j) .EQ. '|') npipes_resgroup = npipes_resgroup + 1
    END DO
    ! find the number of pipes from the right of the line to identify where the residuals group starts
    IF (deprecated_found .AND. remarks_found) THEN
       npipes_right = 2
    ELSEIF ((.NOT. deprecated_found .AND. remarks_found) .OR. (deprecated_found .AND. .NOT. remarks_found)) THEN
       npipes_right = 1
    ELSE
       npipes_right = 0
    ENDIF

    ! clean key record from residuals group
    IF (idx_stop .NE. len_str) THEN
       ades_key_rec(i) = trim(adjustl(ades_key_rec(i)(1:idx_start-1)//ades_key_rec(i)(idx_stop:len_str)))
    ELSE
       ades_key_rec(i) = ades_key_rec(i)(1:idx_start-1)
    ENDIF

    ! loop over next data records belonging to the same segment defined by ades_key_rec(i)
    j = i+1
    DO WHILE (ades_data_rec(j) .NE. "")

       !idx_start = 0  idx_start is the same as the key_record
       idx_stop = 0

       ! find index of npipes_right to define where the residuals group starts from the right in data record (idx_stop)
       ! if npipes_right = 0 then residuals group is up to the end of the string
       len_str = len(ades_data_rec(j))
       idx_stop = len_str
       IF (npipes_right .GT. 0) THEN
          count_pipes = 0
          pipe_found = .FALSE.
          DO WHILE (.NOT. pipe_found)
             IF (ades_data_rec(j)(idx_stop:idx_stop) .EQ. "|") count_pipes = count_pipes + 1
             IF (count_pipes .EQ. npipes_right) pipe_found = .TRUE.
             idx_stop = idx_stop - 1
          END DO
       ENDIF

       ! if indexes are zero it means something went wrong
       IF (idx_start .EQ. 0 .OR. idx_stop .EQ. 0)THEN
          WRITE(line,'(i7)')j
          WRITE(ierrou,*) 'clean_ades_segment_resgroup: unable to find pipe corresponding to &
               & data record residuals group for input ADES file: '//trim(ades_filename)// &
               & ' line: '//trim(adjustl(line))
          numerr=numerr+1 
          RETURN
       ENDIF

       ! clean data record from residuals group
       IF (idx_stop .NE. len_str) THEN
          ades_data_rec(j) = trim(adjustl(ades_data_rec(j)(1:idx_start-1)//ades_data_rec(j)(idx_stop:len_str)))
       ELSE
          ades_data_rec(j) = ades_data_rec(j)(1:idx_start-1)
       ENDIF

       j = j+1

    END DO

  END SUBROUTINE clean_residuals_group

  ! ----------------------------------------------------- !
  ! SUBROUTINE clean_remarks                              !
  ! For a given key record it removes chi and the flags   !
  ! of used fits from the remarks entry of segment        !
  ! data_record, if present. It also updates the          !
  ! corresponding value in ades_obs (ades_obs%remarks).   !
  ! ----------------------------------------------------- !

  SUBROUTINE clean_remarks(i)

    INTEGER,            INTENT(IN)    :: i               ! line number of key record
    CHARACTER(LEN=7) :: line
    LOGICAL :: pipe_found,remarks_found
    INTEGER :: j,idx_start,len_str,idx_obs,idx_chi,idx_A,idx_S
    CHARACTER(LEN=200) :: substr

    ! initialization
    remarks_found = .FALSE.
    IF (INDEX(ades_key_rec(i),"remarks") .GT. 0) remarks_found = .TRUE.

    ! return if remarks has not been found
    IF (.NOT. remarks_found) RETURN

    ! loop over next data records belonging to the same segment defined by ades_key_rec(i)
    j = i+1
    DO WHILE (ades_data_rec(j) .NE. "")

       idx_start = 0
       idx_obs = ntot_to_nobs_idx(j)

       ! remove chi and flags from remarks if present
       ! must be "chi =" or "A =" or "S =" at the end of entry
       len_str = len(ades_data_rec(j))

       ! find index of chi or flags
       idx_chi = INDEX(ades_data_rec(j),"chi =")
       idx_A = INDEX(ades_data_rec(j),"A =")
       idx_S = INDEX(ades_data_rec(j),"S =")

       ! if chi is present
       IF (idx_chi .GT. 0) THEN

          idx_start = idx_chi

          ! remove chi and flags from remarks entry (located at the end of line)
          ades_data_rec(j) = ades_data_rec(j)(1:idx_start-1)

          ! update ades_obs accordingly
          IF (idx_obs .GT. 0) THEN

             IF (ades_flg(idx_obs)%remarks) THEN

                ! remove chi and flags
                ! if only chi and flags in remarks
                IF (ades_obs(idx_obs)%remarks(1:5) .EQ. "chi =") THEN

                   ades_obs(idx_obs)%remarks = ""

                   ! if there's something else in remarks other than chi
                ELSE

                   idx_start = INDEX(ades_obs(idx_obs)%remarks,"chi =")
                   ! clean entry
                   IF (idx_start .GT. 0) THEN
                      substr = ades_obs(idx_obs)%remarks
                      substr = substr(1:idx_start-1)
                      ades_obs(idx_obs)%remarks = substr
                   ENDIF

                ENDIF

             ENDIF

          ENDIF

          ! if there is not chi, find index of flags and remove them  
       ELSE IF (idx_A .GT. 0) THEN
          ! case of optical obs
          idx_start = idx_A

          ! remove flags from remarks entry (located at the end of line)
          ades_data_rec(j) = ades_data_rec(j)(1:idx_start-1)

          ! update ades_obs accordingly
          IF (idx_obs .GT. 0) THEN

             IF (ades_flg(idx_obs)%remarks) THEN

                ! remove flags
                ! if only flags in remarks
                IF (ades_obs(idx_obs)%remarks(1:5) .EQ. "A =") THEN

                   ades_obs(idx_obs)%remarks = ""

                   ! if there's something else in remarks other than flags
                ELSE

                   idx_start = INDEX(ades_obs(idx_obs)%remarks,"A =")
                   ! clean entry
                   IF (idx_start .GT. 0) THEN
                      substr = ades_obs(idx_obs)%remarks
                      substr = substr(1:idx_start-1)
                      ades_obs(idx_obs)%remarks = substr
                   ENDIF

                ENDIF

             ENDIF

          ENDIF

       ELSE IF (idx_S .GT. 0) THEN
          ! case of radar obs
          idx_start = idx_S

          ! remove flag from remarks entry (located at the end of line)
          ades_data_rec(j) = ades_data_rec(j)(1:idx_start-1)

          ! update ades_obs accordingly
          IF (idx_obs .GT. 0) THEN

             IF (ades_flg(idx_obs)%remarks) THEN

                ! remove flag S
                ! if only flag S in remarks
                IF (ades_obs(idx_obs)%remarks(1:5) .EQ. "S =") THEN

                   ades_obs(idx_obs)%remarks = ""

                   ! if there's something else in remarks other than S
                ELSE

                   idx_start = INDEX(ades_obs(idx_obs)%remarks,"S =")
                   ! clean entry
                   IF (idx_start .GT. 0) THEN
                      substr = ades_obs(idx_obs)%remarks
                      substr = substr(1:idx_start-1)
                      ades_obs(idx_obs)%remarks = substr
                   ENDIF

                ENDIF

             ENDIF

          ENDIF

       ENDIF

       j = j+1

    END DO

  END SUBROUTINE clean_remarks

  ! ----------------------------------------------------- !
  ! SUBROUTINE get_residuals_group                        !
  ! For a given segment it finds the corresponding        !
  ! key and data residuals group and remarks, and writes  !
  ! chi value and flags of used fit at the end of remarks !
  ! ----------------------------------------------------- !

  SUBROUTINE get_residuals_group(i,ndata_segment,key_res_group,data_res_group)

    INTEGER,            INTENT(IN)    :: i                             ! line number of key record
    INTEGER,            INTENT(IN)    :: ndata_segment                 ! number of data records in segment
    CHARACTER(LEN=600), INTENT(OUT)   :: key_res_group                 ! residuals group string to be added to key record
    CHARACTER(LEN=600), INTENT(OUT)   :: data_res_group(ndata_segment) ! residuals group string array to be added to data record

    CHARACTER(LEN=40), DIMENSION(ndata_segment,18) :: sfmt ! array containing format strings for integer/real conversions to strings
    INTEGER :: length(20)                                  ! array of maximum lengths for segment

    INTEGER :: j,idx_obs,len_remarks,totlen_data,prec,len_res_group,idx_stop,len_str,count_pipes,npipes_right,k,k_end
    CHARACTER(LEN=9)   :: chi_str
    CHARACTER(LEN=24)  :: flags_str
    CHARACTER(LEN=1)   :: force_fl(2)
    CHARACTER(LEN=300) :: aux_str,data_str,empty_str,remarks_str
    CHARACTER(LEN=20)  :: FMT,prec_str,totlen_str
    LOGICAL :: is_opt,is_del,is_dopp,no_astrofit,no_magfit,deprecated_found,is_deprecated,remarks_found
    REAL(KIND=dkind) :: resid_coord

    INTEGER, PARAMETER :: min_prec = 1    ! minimum precision for data values

    ! initialization
    sfmt = ""
    length = 0
    key_res_group = ""
    data_res_group = ""
    is_opt = .FALSE.
    is_del = .FALSE.
    is_dopp = .FALSE.
    empty_str = ""

    ! check if deprecated and remarks keys are present
    deprecated_found = .FALSE.
    IF (INDEX(ades_key_rec(i),"deprecated") .GT. 0) deprecated_found = .TRUE.
    remarks_found = .FALSE.
    IF (INDEX(ades_key_rec(i),"remarks") .GT. 0) remarks_found = .TRUE.

    ! loop over segment (only data records)
    DO j=i+1,i+ndata_segment

       ! get index of observations
       idx_obs = ntot_to_nobs_idx(j)
       ! we get length only from lines not skipped, notice that we add a +1 to account for initial pipe in data
       IF (idx_obs .GT. 0) THEN
          ! compute len_remarks: remarks (if present) + chi (if defined) + flag for used fit = value (format F9.2)
          IF (ABS(ades_fit(idx_obs)%chi) .LT. 9.9d5)THEN
             WRITE(chi_str,'(F9.2)') ades_fit(idx_obs)%chi
          ELSE
             WRITE(chi_str,'(1P,E9.2)') ades_fit(idx_obs)%chi
          ENDIF
          ! write force_w_flags as strings
          WRITE(force_fl(1),'(L1)') ades_fit(idx_obs)%force_w_flags(1)
          WRITE(force_fl(2),'(L1)') ades_fit(idx_obs)%force_w_flags(2)
          ! flags_str if optical obs
          IF (.NOT. ades_flg(idx_obs)%rad) THEN
             WRITE(flags_str,'(A24)') " A = "//ades_fit(idx_obs)%ast_flag//" M = "//ades_fit(idx_obs)%mag_flag// &
                  & " RMS_F = "//force_fl(1)//" "//force_fl(2)
          ELSE
             ! flags_str if radar obs
             IF(ades_flg(idx_obs)%delay) THEN
                WRITE(flags_str,'(A16)') " S = "//ades_fit(idx_obs)%radar_flag//" RMS_F = "//force_fl(1)
             ELSE
                WRITE(flags_str,'(A16)') " S = "//ades_fit(idx_obs)%radar_flag//" RMS_F = "//force_fl(2)
             ENDIF
          ENDIF
          IF (.NOT. ades_flg(idx_obs)%remarks .OR. trim(adjustl(ades_obs(idx_obs)%remarks)) .EQ. "") THEN
             IF (ades_fit(idx_obs)%ades_resc_def) THEN
                len_remarks = 6+len(trim(adjustl(chi_str)))+len(trim(flags_str))
             ELSE
                len_remarks = len(trim(flags_str))
             ENDIF
          ELSE
             IF (ades_fit(idx_obs)%ades_resc_def) THEN
                len_remarks = len(trim(adjustl((ades_obs(idx_obs)%remarks))))+7+len(trim(adjustl(chi_str)))+len(trim(flags_str))
             ELSE
                len_remarks = len(trim(adjustl((ades_obs(idx_obs)%remarks))))+len(trim(flags_str))
             ENDIF
          ENDIF
          ! remarks - 1
          length(1) = MAX(1+len_remarks,len("|remarks"),length(1))

          ! orbProd - 2
          length(2) = MAX(1+len(trim(adjustl(ades_obs(idx_obs)%orbProd))),len("|orbProd"),length(2))

          ! orbID - 3
          length(3) = MAX(1+len(trim(adjustl(ades_obs(idx_obs)%orbID))),len("|orbID"),length(3))

          ! if optical obs
          IF (.NOT. ades_flg(idx_obs)%rad) THEN

             ! segment is optical
             is_opt = .TRUE.

             ! resRA - 4: conversion degree to arcsec 3600 -> prec ~ -4 + 1 for safety
             ! precision given by sigRA
             resid_coord = ades_obs(idx_obs)%resRA
             IF (abs(resid_coord) > 999.d0) THEN
                sfmt(j-i,4) = "(1P,E9.1)"
                totlen_data = 9
             ELSE
                prec = 3
                totlen_data = 4 + prec + 2  ! +2 for "-" and "."
                WRITE(totlen_str,"(i20)") totlen_data
                WRITE(prec_str,"(i20)") prec
                sfmt(j-i,4) = "(F"//trim(adjustl(totlen_str))//"."//trim(adjustl(prec_str))//")"
             ENDIF
             length(4) = MAX(1+totlen_data,len("|resRA"),length(4))

             ! resDec - 5: conversion degree to arcsec 3600 -> prec ~ 4 + 1 for safety
             ! number of digits before (including) dot = 4 (+/-NN.)
             ! precision given by sigDec
             resid_coord = ades_obs(idx_obs)%resDec
             IF (abs(resid_coord) > 999.d0) THEN
                sfmt(j-i,5) = "(1P,E9.1)"
                totlen_data = 9
             ELSE
                prec = 3
                totlen_data = 4 + prec + 2  ! +2 for "-" and "."
                WRITE(totlen_str,"(i20)") totlen_data
                WRITE(prec_str,"(i20)") prec
                sfmt(j-i,5) = "(F"//trim(adjustl(totlen_str))//"."//trim(adjustl(prec_str))//")"
             ENDIF
             length(5) = MAX(1+totlen_data,len("|resDec"),length(5))

             ! selAst - 6:
             length(6) = len("|selAst")

             ! sigRA - 7
             prec = 3 
             totlen_data = 4 + prec + 1  ! +1 for "."
             WRITE(totlen_str,"(i20)") totlen_data
             WRITE(prec_str,"(i20)") prec
             sfmt(j-i,7) = "(SS,F"//trim(adjustl(totlen_str))//"."//trim(adjustl(prec_str))//")"
             length(7) = 1 + totlen_data

             ! sigDec - 8
             prec = 3 
             totlen_data = 4 + prec + 1  ! +1 for "."
             WRITE(totlen_str,"(i20)") totlen_data
             WRITE(prec_str,"(i20)") prec
             sfmt(j-i,8) = "(SS,F"//trim(adjustl(totlen_str))//"."//trim(adjustl(prec_str))//")"
             length(8) = 1 + totlen_data

             ! sigCorr - 9
             prec = 3
             totlen_data = 1 + prec + 2  ! +2 for "-" and "."
             WRITE(totlen_str,"(i20)") totlen_data
             WRITE(prec_str,"(i20)") prec
             sfmt(j-i,9) = "(F"//trim(adjustl(totlen_str))//"."//trim(adjustl(prec_str))//")"
             length(9) = MAX(1+totlen_data,len("|sigCorr"),length(9))

             ! sigTime - 10
             prec = 3 
             totlen_data = 5 + prec + 1  ! +1 for "."
             WRITE(totlen_str,"(i20)") totlen_data
             WRITE(prec_str,"(i20)") prec
             sfmt(j-i,10) = "(SS,F"//trim(adjustl(totlen_str))//"."//trim(adjustl(prec_str))//")"
             length(10) = MAX(1+totlen_data,len("|sigTime"),length(10))

             ! biasRA - 11: conversion degree to arcsec 3600 -> prec ~ 4 + 1 for safety
             prec = 3 
             totlen_data = 4 + prec + 2  ! +2 for "-" and "."
             WRITE(totlen_str,"(i20)") totlen_data
             WRITE(prec_str,"(i20)") prec
             sfmt(j-i,11) = "(F"//trim(adjustl(totlen_str))//"."//trim(adjustl(prec_str))//")"
             !length(11) = MAX(1+totlen_data,len("|biasRA"),length(11))
             length(11) = 1 + totlen_data

             ! biasDec - 12: conversion degree to arcsec 3600 -> prec ~ 4 + 1 for safety
             prec = 3 
             totlen_data = 4 + prec + 2  ! +2 for "-" and "."
             WRITE(totlen_str,"(i20)") totlen_data
             WRITE(prec_str,"(i20)") prec
             sfmt(j-i,12) = "(F"//trim(adjustl(totlen_str))//"."//trim(adjustl(prec_str))//")"
             length(12) = 1 + totlen_data

             ! biasTime - 13
             prec = 3 
             totlen_data = 5 + prec + 2  ! +2 for "-" and "."
             WRITE(totlen_str,"(i20)") totlen_data
             WRITE(prec_str,"(i20)") prec
             sfmt(j-i,13) = "(F"//trim(adjustl(totlen_str))//"."//trim(adjustl(prec_str))//")"
             length(13) = MAX(1+totlen_data,len("|biasTime"),length(13))

             ! photProd - 14
             length(14) = MAX(1+len(trim(adjustl(ades_obs(idx_obs)%photProd))),len("|photProd"),length(14))

             ! resMag - 15
             ! precision given by sigMag
             prec = 2
             totlen_data = 2 + prec + 2  ! +2 for "-" and "."
             WRITE(totlen_str,"(i20)") totlen_data
             WRITE(prec_str,"(i20)") prec
             sfmt(j-i,15) = "(F"//trim(adjustl(totlen_str))//"."//trim(adjustl(prec_str))//")"
             length(15) = MAX(1+totlen_data,len("|resMag"),length(15))

             ! selPhot - 16:
             length(16) = len("|selPhot")

             ! sigMag - 17
             prec = 3 
             totlen_data = 2 + prec + 1  ! +1 for "."
             WRITE(totlen_str,"(i20)") totlen_data
             WRITE(prec_str,"(i20)") prec
             sfmt(j-i,17) = "(SS,F"//trim(adjustl(totlen_str))//"."//trim(adjustl(prec_str))//")"
             length(17) = MAX(1+totlen_data,len("|sigMag"),length(17))

             ! biasMag - 18:
             prec = 3 
             totlen_data = 2 + prec + 2  ! +2 for "-" and "."
             WRITE(totlen_str,"(i20)") totlen_data
             WRITE(prec_str,"(i20)") prec
             sfmt(j-i,18) = "(SS,F"//trim(adjustl(totlen_str))//"."//trim(adjustl(prec_str))//")"
             length(18) = MAX(1+totlen_data,len("|biasMag"),length(18))

             ! photMod - 19:
             length(19) = len("|photMod")

             ! if radar obs
          ELSE

             ! if delay
             IF (ades_flg(idx_obs)%delay) THEN

                ! segment is radar with delay obs
                is_del = .TRUE.

                ! resDelay - 4: conversion s to mus 1.e6 -> prec ~ -6
                ! precision given by sigDelay
                resid_coord = ades_obs(idx_obs)%resDelay
                IF (resid_coord > 99999.d0 .OR. resid_coord < -9999.d0) THEN
                   sfmt(j-i,4) = "(1P,E11.4)"
                   totlen_data = 11
                ELSE
                   prec = 5
                   totlen_data = 4 + prec + 2  ! +2 for "-" and "."
                   WRITE(totlen_str,"(i20)") totlen_data
                   WRITE(prec_str,"(i20)") prec
                   sfmt(j-i,4) = "(F"//trim(adjustl(totlen_str))//"."//trim(adjustl(prec_str))//")"
                ENDIF
                length(4) = MAX(1+totlen_data,len("|resDelay"),length(4))

                ! selDelay - 5
                length(5) = len("|selDelay")

                ! sigDelay - 6
                prec = 3 
                totlen_data = 4 + prec + 1  ! +1 for "."
                WRITE(totlen_str,"(i20)") totlen_data
                WRITE(prec_str,"(i20)") prec
                sfmt(j-i,6) = "(SS,F"//trim(adjustl(totlen_str))//"."//trim(adjustl(prec_str))//")"
                length(6) = MAX(1+totlen_data,len("|sigDelay"),length(6))

                ! if doppler
             ELSEIF (ades_flg(idx_obs)%doppler) THEN

                ! segment is radar with doppler obs
                is_dopp = .TRUE.

                ! resDoppler - 4
                ! precision given by sigDoppler
                resid_coord = ades_obs(idx_obs)%resDoppler
                IF (resid_coord > 99999.d0 .OR. resid_coord < -9999.d0) THEN
                   sfmt(j-i,4) = "(1P,E11.4)"
                   totlen_data = 11
                ELSE
                   prec = 5
                   totlen_data = 4 + prec + 2  ! +2 for "-" and "."
                   WRITE(totlen_str,"(i20)") totlen_data
                   WRITE(prec_str,"(i20)") prec
                   sfmt(j-i,4) = "(F"//trim(adjustl(totlen_str))//"."//trim(adjustl(prec_str))//")"
                ENDIF
                length(4) = MAX(1+totlen_data,len("|resDoppler"),length(4))

                ! selDoppler - 5
                length(5) = len("|selDoppler")

                ! sigDoppler - 6
                prec = 3
                totlen_data = 4 + prec + 1  ! +1 for "."
                WRITE(totlen_str,"(i20)") totlen_data
                WRITE(prec_str,"(i20)") prec
                sfmt(j-i,6) = "(SS,F"//trim(adjustl(totlen_str))//"."//trim(adjustl(prec_str))//")"
                length(6) = MAX(1+totlen_data,len("|sigDoppler"),length(6))

             ENDIF

          ENDIF
       ENDIF

       ! deprecated - 20
       length(20) = len("|deprecated")

    END DO

    ! format key entries
    ! orbProd - 2
    WRITE(FMT,"(i20)") length(2)
    WRITE(aux_str,"(A"//trim(adjustl(FMT))//")") "|orbProd"
    key_res_group = adjustl(aux_str(1:length(2)))
    len_res_group = length(2)
    ! orbID - 3
    WRITE(FMT,"(i20)") length(3)
    WRITE(aux_str,"(A"//trim(adjustl(FMT))//")") "|orbID"
    key_res_group = key_res_group(1:len_res_group)//adjustl(aux_str(1:length(3)))
    len_res_group = len_res_group + length(3)
    ! if optical obs
    IF (is_opt) THEN
       ! resRA - 4
       WRITE(FMT,"(i20)") length(4)
       WRITE(aux_str,"(A"//trim(adjustl(FMT))//")") "|resRA"
       key_res_group = key_res_group(1:len_res_group)//adjustl(aux_str(1:length(4)))
       len_res_group = len_res_group + length(4)
       ! resDec - 5
       WRITE(FMT,"(i20)") length(5)
       WRITE(aux_str,"(A"//trim(adjustl(FMT))//")") "|resDec"
       key_res_group = key_res_group(1:len_res_group)//adjustl(aux_str(1:length(5)))
       len_res_group = len_res_group + length(5)
       ! selAst - 6
       WRITE(FMT,"(i20)") length(6)
       WRITE(aux_str,"(A"//trim(adjustl(FMT))//")") "|selAst"
       key_res_group = key_res_group(1:len_res_group)//adjustl(aux_str(1:length(6)))
       len_res_group = len_res_group + length(6)
       ! sigRA - 7
       WRITE(FMT,"(i20)") length(7)
       WRITE(aux_str,"(A"//trim(adjustl(FMT))//")") "|sigRA"
       key_res_group = key_res_group(1:len_res_group)//adjustl(aux_str(1:length(7)))
       len_res_group = len_res_group + length(7)
       ! sigDec - 8
       WRITE(FMT,"(i20)") length(8)
       WRITE(aux_str,"(A"//trim(adjustl(FMT))//")") "|sigDec"
       key_res_group = key_res_group(1:len_res_group)//adjustl(aux_str(1:length(8)))
       len_res_group = len_res_group + length(8)
       ! sigCorr - 9
       WRITE(FMT,"(i20)") length(9)
       WRITE(aux_str,"(A"//trim(adjustl(FMT))//")") "|sigCorr"
       key_res_group = key_res_group(1:len_res_group)//adjustl(aux_str(1:length(9)))
       len_res_group = len_res_group + length(9)
       ! sigTime - 10
       WRITE(FMT,"(i20)") length(10)
       WRITE(aux_str,"(A"//trim(adjustl(FMT))//")") "|sigTime"
       key_res_group = key_res_group(1:len_res_group)//adjustl(aux_str(1:length(10)))
       len_res_group = len_res_group + length(10)
       ! biasRA - 11
       WRITE(FMT,"(i20)") length(11)
       WRITE(aux_str,"(A"//trim(adjustl(FMT))//")") "|biasRA"
       key_res_group = key_res_group(1:len_res_group)//adjustl(aux_str(1:length(11)))
       len_res_group = len_res_group + length(11)
       ! biasDec - 12
       WRITE(FMT,"(i20)") length(12)
       WRITE(aux_str,"(A"//trim(adjustl(FMT))//")") "|biasDec"
       key_res_group = key_res_group(1:len_res_group)//adjustl(aux_str(1:length(12)))
       len_res_group = len_res_group + length(12)
       ! biasTime - 13
       WRITE(FMT,"(i20)") length(13)
       WRITE(aux_str,"(A"//trim(adjustl(FMT))//")") "|biasTime"
       key_res_group = key_res_group(1:len_res_group)//adjustl(aux_str(1:length(13)))
       len_res_group = len_res_group + length(13)

       ! add Photometry subgroup
       ! photProd - 14
       WRITE(FMT,"(i20)") length(14)
       WRITE(aux_str,"(A"//trim(adjustl(FMT))//")") "|photProd"
       key_res_group = key_res_group(1:len_res_group)//adjustl(aux_str(1:length(14)))
       len_res_group = len_res_group + length(14)
       ! resMag - 15
       WRITE(FMT,"(i20)") length(15)
       WRITE(aux_str,"(A"//trim(adjustl(FMT))//")") "|resMag"
       key_res_group = key_res_group(1:len_res_group)//adjustl(aux_str(1:length(15)))
       len_res_group = len_res_group + length(15)
       ! selPhot - 16
       WRITE(FMT,"(i20)") length(16)
       WRITE(aux_str,"(A"//trim(adjustl(FMT))//")") "|selPhot"
       key_res_group = key_res_group(1:len_res_group)//adjustl(aux_str(1:length(16)))
       len_res_group = len_res_group + length(16)
       ! sigMag - 17
       WRITE(FMT,"(i20)") length(17)
       WRITE(aux_str,"(A"//trim(adjustl(FMT))//")") "|sigMag"
       key_res_group = key_res_group(1:len_res_group)//adjustl(aux_str(1:length(17)))
       len_res_group = len_res_group + length(17)
       ! biasMag - 18
       WRITE(FMT,"(i20)") length(18)
       WRITE(aux_str,"(A"//trim(adjustl(FMT))//")") "|biasMag"
       key_res_group = key_res_group(1:len_res_group)//adjustl(aux_str(1:length(18)))
       len_res_group = len_res_group + length(18)
       ! photMod - 19
       WRITE(FMT,"(i20)") length(19)
       WRITE(aux_str,"(A"//trim(adjustl(FMT))//")") "|photMod"
       key_res_group = key_res_group(1:len_res_group)//adjustl(aux_str(1:length(19)))
       len_res_group = len_res_group + length(19)

       ! if radar delay obs
    ELSE IF (is_del) THEN
       ! resDelay - 4
       WRITE(FMT,"(i20)") length(4)
       WRITE(aux_str,"(A"//trim(adjustl(FMT))//")") "|resDelay"
       key_res_group = key_res_group(1:len_res_group)//adjustl(aux_str(1:length(4)))
       len_res_group = len_res_group + length(4)
       ! selDelay - 5
       WRITE(FMT,"(i20)") length(5)
       WRITE(aux_str,"(A"//trim(adjustl(FMT))//")") "|selDelay"
       key_res_group = key_res_group(1:len_res_group)//adjustl(aux_str(1:length(5)))
       len_res_group = len_res_group + length(5)
       ! sigDelay - 6
       WRITE(FMT,"(i20)") length(6)
       WRITE(aux_str,"(A"//trim(adjustl(FMT))//")") "|sigDelay"
       key_res_group = key_res_group(1:len_res_group)//adjustl(aux_str(1:length(6)))
       len_res_group = len_res_group + length(6)
       ! if radar doppler obs
    ELSE IF (is_dopp) THEN
       ! resDoppler - 4
       WRITE(FMT,"(i20)") length(4)
       WRITE(aux_str,"(A"//trim(adjustl(FMT))//")") "|resDoppler"
       key_res_group = key_res_group(1:len_res_group)//adjustl(aux_str(1:length(4)))
       len_res_group = len_res_group + length(4)
       ! selDoppler - 5
       WRITE(FMT,"(i20)") length(5)
       WRITE(aux_str,"(A"//trim(adjustl(FMT))//")") "|selDoppler"
       key_res_group = key_res_group(1:len_res_group)//adjustl(aux_str(1:length(5)))
       len_res_group = len_res_group + length(5)
       ! sigDoppler - 6
       WRITE(FMT,"(i20)") length(6)
       WRITE(aux_str,"(A"//trim(adjustl(FMT))//")") "|sigDoppler"
       key_res_group = key_res_group(1:len_res_group)//adjustl(aux_str(1:length(6)))
       len_res_group = len_res_group + length(6)
    ENDIF
    ! if deprecated - 20
    IF (deprecated_found) THEN
       WRITE(FMT,"(i20)") length(20)
       WRITE(aux_str,"(A"//trim(adjustl(FMT))//")") "|deprecated"
       key_res_group = key_res_group(1:len_res_group)//adjustl(aux_str(1:length(20)))
       len_res_group = len_res_group + length(20)
    ENDIF
    ! remarks - 1
    WRITE(FMT,"(i20)") length(1)
    WRITE(aux_str,"(A"//trim(adjustl(FMT))//")") "|remarks"
    key_res_group = key_res_group(1:len_res_group)//adjustl(aux_str(1:length(1)))
    len_res_group = len_res_group + length(1)

    ! loop over and format data entries
    DO j=i+1,i+ndata_segment

       ! get index of observations
       idx_obs = ntot_to_nobs_idx(j)

       ! if skipped line add empty (+ X for deprecated if present) residuals group
       IF (idx_obs .EQ. 0) THEN

          ! set number of empty fields depending on whether the observation is optical or not
          IF (is_opt) THEN
             k_end = 19
          ELSE
             k_end = 6
          ENDIF

          ! fill empty residuals group
          data_str = empty_str
          WRITE(FMT,"(i20)") length(2)
          WRITE(aux_str,"(A"//trim(adjustl(FMT))//")") "|"//data_str
          data_res_group(j-i) = adjustl(aux_str(1:length(2)))
          len_res_group = length(2)
          DO k=3,k_end
             WRITE(FMT,"(i20)") length(k)
             WRITE(aux_str,"(A"//trim(adjustl(FMT))//")") "|"//data_str
             data_res_group(j-i) = data_res_group(j-i)(1:len_res_group)//adjustl(aux_str(1:length(k)))
             len_res_group = len_res_group + length(k)
          ENDDO

          ! add deprecated
          ! if there is a deprecated entry check if it is "X" or blank
          is_deprecated = .FALSE.
          IF (deprecated_found) THEN
             len_str = len(ades_data_rec(j))
             idx_stop = len_str
             count_pipes = 0
             ! deprecated is either the last or the second to last entry depending on the presence of remarks
             IF (remarks_found) THEN
                npipes_right = 2
             ELSE
                npipes_right = 1
             ENDIF
             remarks_str = ""
             DO WHILE (count_pipes .LT. npipes_right)
                idx_stop = idx_stop - 1
                IF (count_pipes .EQ. npipes_right-1 .AND. ades_data_rec(j)(idx_stop:idx_stop) .EQ. "X") THEN
                   is_deprecated = .TRUE.
                   EXIT
                ENDIF
                IF (ades_data_rec(j)(idx_stop:idx_stop) .EQ. "|") count_pipes = count_pipes + 1            
             ENDDO
             WRITE(FMT,"(i20)") length(20)
             IF (is_deprecated) THEN
                WRITE(aux_str,"(A"//trim(adjustl(FMT))//")") "|X"
             ELSE
                WRITE(aux_str,"(A"//trim(adjustl(FMT))//")") "| "
             ENDIF
             data_res_group(j-i) = data_res_group(j-i)(1:len_res_group)//adjustl(aux_str(1:length(20)))
             len_res_group = len_res_group + length(20)
          ENDIF

          ! add remarks if present
          IF (remarks_found) THEN
             ! save remarks entry
             remarks_str = ades_data_rec(j)(INDEX(ades_data_rec(j),"|",BACK=.TRUE.)+1:len(trim(ades_data_rec(j))))
             ! write remarks
             WRITE(FMT,"(i20)") length(1)
             WRITE(aux_str,"(A"//trim(adjustl(FMT))//")") "|"//remarks_str
             data_res_group(j-i) = data_res_group(j-i)(1:len_res_group)//adjustl(aux_str(1:length(1)))
             len_res_group = len_res_group + length(1)
          ENDIF

          ! only non-skipped lines
       ELSEIF (idx_obs .GT. 0) THEN

          ! find whether fit has been performed
          no_astrofit = .FALSE.
          IF (is_opt) THEN
             IF (trim(adjustl(ades_obs(idx_obs)%selAst)) .EQ. "D" .OR. trim(adjustl(ades_obs(idx_obs)%selAst)) .EQ. "d") THEN
                no_astrofit = .TRUE.
             ENDIF
             no_magfit = .FALSE.
             ! IF (trim(adjustl(ades_obs(idx_obs)%selPhot)) .EQ. "D" .OR. trim(adjustl(ades_obs(idx_obs)%selPhot)) .EQ. "d" .OR. &
             !      & ades_obs(idx_obs)%mag .EQ. 9.9d9) THEN
             IF(ades_obs(idx_obs)%mag .EQ. 9.9d9 .OR. ades_obs(idx_obs)%mag .EQ. 0.d0) THEN
                no_magfit = .TRUE.
             ENDIF
          ELSEIF (is_del) THEN
             IF (trim(adjustl(ades_obs(idx_obs)%selDelay)) .EQ. "D" .OR. trim(adjustl(ades_obs(idx_obs)%selDelay)) .EQ. "d") THEN
                no_astrofit = .TRUE.
             ENDIF
          ELSEIF (is_dopp) THEN
             IF (trim(adjustl(ades_obs(idx_obs)%selDoppler)) .EQ. "D" .OR. &
                  & trim(adjustl(ades_obs(idx_obs)%selDoppler)) .EQ. "d") THEN
                no_astrofit = .TRUE.
             ENDIF
          ENDIF

          ! orbProd - 2
          WRITE(FMT,"(i20)") length(2)
          WRITE(aux_str,"(A"//trim(adjustl(FMT))//")") "|"//trim(adjustl(ades_obs(idx_obs)%orbProd))
          data_res_group(j-i) = adjustl(aux_str(1:length(2)))
          len_res_group = length(2)
          ! orbID - 3
          WRITE(FMT,"(i20)") length(3)
          WRITE(aux_str,"(A"//trim(adjustl(FMT))//")") "|"//trim(adjustl(ades_obs(idx_obs)%orbID))
          data_res_group(j-i) = data_res_group(j-i)(1:len_res_group)//adjustl(aux_str(1:length(3)))
          len_res_group = len_res_group + length(3)

          ! if optical obs
          IF (is_opt) THEN
             ! resRA - 4
             WRITE(data_str,sfmt(j-i,4)) ades_obs(idx_obs)%resRA
             IF (.NOT. ades_fit(idx_obs)%ades_resc_def) data_str = empty_str
             WRITE(FMT,"(i20)") length(4)
             WRITE(aux_str,"(A"//trim(adjustl(FMT))//")") "|"//data_str
             data_res_group(j-i) = data_res_group(j-i)(1:len_res_group)//adjustl(aux_str(1:length(4)))
             len_res_group = len_res_group + length(4)
             ! resDec - 5
             WRITE(data_str,sfmt(j-i,5)) ades_obs(idx_obs)%resDec
             IF (.NOT. ades_fit(idx_obs)%ades_resc_def) data_str = empty_str
             WRITE(FMT,"(i20)") length(5)
             WRITE(aux_str,"(A"//trim(adjustl(FMT))//")") "|"//data_str
             data_res_group(j-i) = data_res_group(j-i)(1:len_res_group)//adjustl(aux_str(1:length(5)))
             len_res_group = len_res_group + length(5)
             ! selAst - 6
             WRITE(FMT,"(i20)") length(6)
             WRITE(aux_str,"(A"//trim(adjustl(FMT))//")") "|"//trim(adjustl(ades_obs(idx_obs)%selAst))
             data_res_group(j-i) = data_res_group(j-i)(1:len_res_group)//adjustl(aux_str(1:length(6)))
             len_res_group = len_res_group + length(6)
             ! sigRA - 7
             WRITE(data_str,sfmt(j-i,7)) ades_obs(idx_obs)%sigRA
             !IF (no_astrofit) data_str = empty_str
             WRITE(FMT,"(i20)") length(7)
             WRITE(aux_str,"(A"//trim(adjustl(FMT))//")") "|"//data_str
             data_res_group(j-i) = data_res_group(j-i)(1:len_res_group)//adjustl(aux_str(1:length(7)))
             len_res_group = len_res_group + length(7)
             ! sigDec - 8
             WRITE(data_str,sfmt(j-i,8)) ades_obs(idx_obs)%sigDec
             !IF (no_astrofit) data_str = empty_str
             WRITE(FMT,"(i20)") length(8)
             WRITE(aux_str,"(A"//trim(adjustl(FMT))//")") "|"//data_str
             data_res_group(j-i) = data_res_group(j-i)(1:len_res_group)//adjustl(aux_str(1:length(8)))
             len_res_group = len_res_group + length(8)
             ! sigCorr - 9
             WRITE(data_str,sfmt(j-i,9)) ades_obs(idx_obs)%sigCorr
             !IF (no_astrofit) data_str = empty_str
             WRITE(FMT,"(i20)") length(9)
             WRITE(aux_str,"(A"//trim(adjustl(FMT))//")") "|"//data_str
             data_res_group(j-i) = data_res_group(j-i)(1:len_res_group)//adjustl(aux_str(1:length(9)))
             len_res_group = len_res_group + length(9)
             ! sigTime - 10
             WRITE(data_str,sfmt(j-i,10)) ades_obs(idx_obs)%sigTime
             !IF (no_astrofit) data_str = empty_str
             WRITE(FMT,"(i20)") length(10)
             WRITE(aux_str,"(A"//trim(adjustl(FMT))//")") "|"//data_str
             data_res_group(j-i) = data_res_group(j-i)(1:len_res_group)//adjustl(aux_str(1:length(10)))
             len_res_group = len_res_group + length(10)
             ! biasRA - 11
             WRITE(data_str,sfmt(j-i,11)) ades_obs(idx_obs)%biasRA
             !IF (no_astrofit) data_str = empty_str
             WRITE(FMT,"(i20)") length(11)
             WRITE(aux_str,"(A"//trim(adjustl(FMT))//")") "|"//data_str
             data_res_group(j-i) = data_res_group(j-i)(1:len_res_group)//adjustl(aux_str(1:length(11)))
             len_res_group = len_res_group + length(11)
             ! biasDec - 12
             WRITE(data_str,sfmt(j-i,12)) ades_obs(idx_obs)%biasDec
             !IF (no_astrofit) data_str = empty_str
             WRITE(FMT,"(i20)") length(12)
             WRITE(aux_str,"(A"//trim(adjustl(FMT))//")") "|"//data_str
             data_res_group(j-i) = data_res_group(j-i)(1:len_res_group)//adjustl(aux_str(1:length(12)))
             len_res_group = len_res_group + length(12)
             ! biasTime - 13
             WRITE(data_str,sfmt(j-i,13)) ades_obs(idx_obs)%biasTime
             !IF (no_astrofit) data_str = empty_str
             WRITE(FMT,"(i20)") length(13)
             WRITE(aux_str,"(A"//trim(adjustl(FMT))//")") "|"//data_str
             data_res_group(j-i) = data_res_group(j-i)(1:len_res_group)//adjustl(aux_str(1:length(13)))
             len_res_group = len_res_group + length(13)

             ! photProd - 14
             WRITE(FMT,"(i20)") length(14)
             WRITE(aux_str,"(A"//trim(adjustl(FMT))//")") "|"//trim(adjustl(ades_obs(idx_obs)%photProd))
             data_res_group(j-i) = data_res_group(j-i)(1:len_res_group)//adjustl(aux_str(1:length(14)))
             len_res_group = len_res_group + length(14)
             ! resMag - 15
             WRITE(data_str,sfmt(j-i,15)) ades_obs(idx_obs)%resMag
             IF (no_magfit .OR. .NOT. ades_fit(idx_obs)%ades_resm_def) data_str = empty_str
             WRITE(FMT,"(i20)") length(15)
             WRITE(aux_str,"(A"//trim(adjustl(FMT))//")") "|"//data_str
             data_res_group(j-i) = data_res_group(j-i)(1:len_res_group)//adjustl(aux_str(1:length(15)))
             len_res_group = len_res_group + length(15)
             ! selPhot - 16
             WRITE(FMT,"(i20)") length(16)
             WRITE(aux_str,"(A"//trim(adjustl(FMT))//")") "|"//trim(adjustl(ades_obs(idx_obs)%selPhot))
             data_res_group(j-i) = data_res_group(j-i)(1:len_res_group)//adjustl(aux_str(1:length(16)))
             len_res_group = len_res_group + length(16)
             ! sigMag - 17
             WRITE(data_str,sfmt(j-i,17)) ades_obs(idx_obs)%sigMag
             IF (no_magfit) data_str = empty_str
             WRITE(FMT,"(i20)") length(17)
             WRITE(aux_str,"(A"//trim(adjustl(FMT))//")") "|"//data_str
             data_res_group(j-i) = data_res_group(j-i)(1:len_res_group)//adjustl(aux_str(1:length(17)))
             len_res_group = len_res_group + length(17)
             ! biasMag - 18
             WRITE(data_str,sfmt(j-i,18)) ades_obs(idx_obs)%biasMag
             IF (no_magfit) data_str = empty_str
             WRITE(FMT,"(i20)") length(18)
             WRITE(aux_str,"(A"//trim(adjustl(FMT))//")") "|"//data_str
             data_res_group(j-i) = data_res_group(j-i)(1:len_res_group)//adjustl(aux_str(1:length(18)))
             len_res_group = len_res_group + length(18)
             ! photMod - 19
             data_str = ades_obs(idx_obs)%photMod
             IF (no_magfit) data_str = empty_str
             WRITE(FMT,"(i20)") length(19)
             WRITE(aux_str,"(A"//trim(adjustl(FMT))//")") "|"//data_str
             data_res_group(j-i) = data_res_group(j-i)(1:len_res_group)//adjustl(aux_str(1:length(19)))
             len_res_group = len_res_group + length(19)

             ! if radar delay obs
          ELSE IF (is_del) THEN
             ! resDelay - 4
             WRITE(data_str,sfmt(j-i,4)) ades_obs(idx_obs)%resDelay
             IF (.NOT. ades_fit(idx_obs)%ades_resc_def) data_str = empty_str 
             WRITE(FMT,"(i20)") length(4)
             WRITE(aux_str,"(A"//trim(adjustl(FMT))//")") "|"//data_str
             data_res_group(j-i) = data_res_group(j-i)(1:len_res_group)//adjustl(aux_str(1:length(4)))
             len_res_group = len_res_group + length(4)
             ! selDelay - 5
             WRITE(FMT,"(i20)") length(5)
             WRITE(aux_str,"(A"//trim(adjustl(FMT))//")") "|"//trim(adjustl(ades_obs(idx_obs)%selDelay))
             data_res_group(j-i) = data_res_group(j-i)(1:len_res_group)//adjustl(aux_str(1:length(5)))
             len_res_group = len_res_group + length(5)
             ! sigDelay - 6
             WRITE(data_str,sfmt(j-i,6)) ades_obs(idx_obs)%sigDelay
             !IF (no_astrofit) data_str = empty_str
             WRITE(FMT,"(i20)") length(6)
             WRITE(aux_str,"(A"//trim(adjustl(FMT))//")") "|"//data_str
             data_res_group(j-i) = data_res_group(j-i)(1:len_res_group)//adjustl(aux_str(1:length(6)))
             len_res_group = len_res_group + length(6)

             ! if radar doppler obs
          ELSE IF (is_dopp) THEN
             ! resDoppler - 4
             WRITE(data_str,sfmt(j-i,4)) ades_obs(idx_obs)%resDoppler
             IF (.NOT. ades_fit(idx_obs)%ades_resc_def) data_str = empty_str
             WRITE(FMT,"(i20)") length(4)
             WRITE(aux_str,"(A"//trim(adjustl(FMT))//")") "|"//data_str
             data_res_group(j-i) = data_res_group(j-i)(1:len_res_group)//adjustl(aux_str(1:length(4)))
             len_res_group = len_res_group + length(4)
             ! selDoppler - 5
             WRITE(FMT,"(i20)") length(5)
             WRITE(aux_str,"(A"//trim(adjustl(FMT))//")") "|"//trim(adjustl(ades_obs(idx_obs)%selDoppler))
             data_res_group(j-i) = data_res_group(j-i)(1:len_res_group)//adjustl(aux_str(1:length(5)))
             len_res_group = len_res_group + length(5)
             ! sigDoppler - 6
             WRITE(data_str,sfmt(j-i,6)) ades_obs(idx_obs)%sigDoppler
             !IF (no_astrofit) data_str = empty_str
             WRITE(FMT,"(i20)") length(6)
             WRITE(aux_str,"(A"//trim(adjustl(FMT))//")") "|"//data_str
             data_res_group(j-i) = data_res_group(j-i)(1:len_res_group)//adjustl(aux_str(1:length(6)))
             len_res_group = len_res_group + length(6)

          ENDIF

          ! add deprecated empty entry - 20
          IF (deprecated_found) THEN
             WRITE(FMT,"(i20)") length(20)
             WRITE(aux_str,"(A"//trim(adjustl(FMT))//")") "|"//trim(adjustl(ades_obs(idx_obs)%deprecated))
             data_res_group(j-i) = data_res_group(j-i)(1:len_res_group)//adjustl(aux_str(1:length(20)))
             len_res_group = len_res_group + length(20)
          ENDIF

          ! remarks + chi + flag for fit - 1
          ! write force_w_flags as strings
          WRITE(force_fl(1),'(L1)') ades_fit(idx_obs)%force_w_flags(1)
          WRITE(force_fl(2),'(L1)') ades_fit(idx_obs)%force_w_flags(2)
          ! if optical obs
          IF (is_opt) THEN
             WRITE(flags_str,'(A24)') " A = "//ades_fit(idx_obs)%ast_flag//" M = "//ades_fit(idx_obs)%mag_flag// &
                  & " RMS_F = "//force_fl(1)//" "//force_fl(2)
             WRITE(FMT,"(i20)") length(1)
             IF (ades_fit(idx_obs)%ades_resc_def) THEN
                IF (ABS(ades_fit(idx_obs)%chi) .LT. 9.9d5)THEN
                   WRITE(chi_str,'(F9.2)') ades_fit(idx_obs)%chi
                ELSE
                   WRITE(chi_str,'(1P,E9.2)') ades_fit(idx_obs)%chi
                ENDIF
                IF (.NOT. ades_flg(idx_obs)%remarks .OR. trim(adjustl(ades_obs(idx_obs)%remarks)) .EQ. "") THEN
                   WRITE(aux_str,"(A"//trim(adjustl(FMT))//")") &
                        & "|chi = "//trim(adjustl(chi_str))//trim(flags_str)
                ELSE
                   WRITE(aux_str,"(A"//trim(adjustl(FMT))//")") &
                        & "|"//trim(adjustl(ades_obs(idx_obs)%remarks))//" chi = "//trim(adjustl(chi_str))//trim(flags_str)
                ENDIF
                data_res_group(j-i) = data_res_group(j-i)(1:len_res_group)//adjustl(aux_str(1:length(1)))
                len_res_group = len_res_group + length(1)

             ELSE
                IF (.NOT. ades_flg(idx_obs)%remarks .OR. trim(adjustl(ades_obs(idx_obs)%remarks)) .EQ. "") THEN
                   WRITE(aux_str,"(A"//trim(adjustl(FMT))//")") &
                        & "|"//trim(adjustl(flags_str))
                ELSE
                   WRITE(aux_str,"(A"//trim(adjustl(FMT))//")") "|"//trim(adjustl(ades_obs(idx_obs)%remarks))//trim(flags_str)
                ENDIF
                data_res_group(j-i) = data_res_group(j-i)(1:len_res_group)//adjustl(aux_str(1:length(1)))
                len_res_group = len_res_group + length(1)
             ENDIF
             ! if radar obs  
          ELSE IF (is_del .OR. is_dopp) THEN
             IF (is_del) THEN
                WRITE(flags_str,'(A16)') " S = "//ades_fit(idx_obs)%radar_flag//" RMS_F = "//force_fl(1)
             ELSE
                WRITE(flags_str,'(A16)') " S = "//ades_fit(idx_obs)%radar_flag//" RMS_F = "//force_fl(2)
             ENDIF
             WRITE(FMT,"(i20)") length(1)
             IF (ades_fit(idx_obs)%ades_resc_def) THEN
                IF (ABS(ades_fit(idx_obs)%chi) .LT. 9.9d5)THEN
                   WRITE(chi_str,'(F9.2)') ades_fit(idx_obs)%chi
                ELSE
                   WRITE(chi_str,'(1P,E9.2)') ades_fit(idx_obs)%chi
                ENDIF
                IF (.NOT. ades_flg(idx_obs)%remarks .OR. trim(adjustl(ades_obs(idx_obs)%remarks)) .EQ. "") THEN
                   WRITE(aux_str,"(A"//trim(adjustl(FMT))//")") &
                        & "|chi = "//trim(adjustl(chi_str))//trim(flags_str)
                ELSE
                   WRITE(aux_str,"(A"//trim(adjustl(FMT))//")") &
                        & "|"//trim(adjustl(ades_obs(idx_obs)%remarks))//" chi = "//trim(adjustl(chi_str))//trim(flags_str)
                ENDIF
                data_res_group(j-i) = data_res_group(j-i)(1:len_res_group)//adjustl(aux_str(1:length(1)))
                len_res_group = len_res_group + length(1)
             ELSE
                IF (.NOT. ades_flg(idx_obs)%remarks .OR. trim(adjustl(ades_obs(idx_obs)%remarks)) .EQ. "") THEN
                   WRITE(aux_str,"(A"//trim(adjustl(FMT))//")") &
                        & "|"//trim(adjustl(flags_str))
                ELSE
                   WRITE(aux_str,"(A"//trim(adjustl(FMT))//")") "|"//trim(adjustl(ades_obs(idx_obs)%remarks))//trim(flags_str)
                ENDIF
                data_res_group(j-i) = data_res_group(j-i)(1:len_res_group)//adjustl(aux_str(1:length(1)))
                len_res_group = len_res_group + length(1)
             ENDIF
          ENDIF

       ENDIF

    END DO

  END SUBROUTINE get_residuals_group

  ! ----------------------------------------------------- !
  ! SUBROUTINE get_astCat                                 !
  ! It converts the astrometric catalogue code used in    !
  ! error model from MPC catcod to ADES astCat            !
  ! ----------------------------------------------------- !

  SUBROUTINE get_astCat(icode,ocode)

    CHARACTER(LEN=*), INTENT(IN)  :: icode ! MPC code 
    CHARACTER(LEN=8), INTENT(OUT) :: ocode ! ADES code for astrometric catalogue (astCat)

    ocode = " "

    SELECT CASE (trim(adjustl(icode)))
    CASE("a")
       ocode = "USNOA1"
    CASE("b")
       ocode = "USNOSA1"
    CASE("c")
       ocode = "USNOA2"
    CASE("d")
       ocode = "USNOSA2"
    CASE("e")
       ocode = "UCAC1"
    CASE("f")
       ocode = "Tyc1"
    CASE("g")
       ocode = "Tyc2"
    CASE("h")
       ocode = "GSC1.0"
    CASE("i")
       ocode = "GSC1.1"
    CASE("j")
       ocode = "GSC1.2"
    CASE("k")
       ocode = "GSC2.2"
    CASE("l")
       ocode = "ACT"
    CASE("m")
       ocode = "GSCACT"
    CASE("n")
       ocode = "SDSS8"
    CASE("o")
       ocode = "USNOB1"
    CASE("p")
       ocode = "PPM"
    CASE("q")
       ocode = "UCAC4"
    CASE("r")
       ocode = "UCAC2"
    CASE("s")  
       ocode = "USNOB2"
    CASE("t")
       ocode = "PPMXL"
    CASE("u")
       ocode = "UCAC3"
    CASE("v")
       ocode = "NOMAD"
    CASE("w")
       ocode = "CMC14"
    CASE("x")
       ocode = "Hip2"
    CASE("y")
       ocode = "Hip1"
    CASE("z")
       ocode = "GSC"
    CASE("A")
       ocode = "AC"
    CASE("B")
       ocode = "SAO1984"
    CASE("C")
       ocode = "SAO"
    CASE("D")
       ocode = "AGK3"
    CASE("E")
       ocode = "FK4"
    CASE("F")
       ocode = "ACRS"
    CASE("G")
       ocode = "LickGas"
    CASE("H")
       ocode = "Ida93"
    CASE("I")
       ocode = "Perth70"
    CASE("J")
       ocode = "COSMOS"
    CASE("K")
       ocode = "Yale"
    CASE("L")
       ocode = "2MASS"
    CASE("M")
       ocode = "GSC2.3"
    CASE("N")
       ocode = "SDSS7"
    CASE("O")
       ocode = "SSTRC1"
    CASE("P")
       ocode = "MPOSC3"
    CASE("Q")
       ocode = "CMC15"
    CASE("R")
       ocode = "SSTRC4"
    CASE("S")
       ocode = "URAT1"
    CASE("T")
       ocode = "URAT2"
    CASE("U")
       ocode = "Gaia1"
    CASE("V")
       ocode = "Gaia2"
    CASE("W")
       ocode = "UCAC5"
    CASE(" ")
       ocode = "UNK"
    END SELECT

    IF (ocode .EQ. " ") THEN
       WRITE(ierrou,*) 'get_astCat: unable to convert from catcodmpc to astCat for output ADES file '//trim(ades_filename)// &
            & " Input code: "//trim(adjustl(icode))
       numerr=numerr+1 
       RETURN
    ENDIF

  END SUBROUTINE get_astCat

  ! ----------------------------------------------- !
  ! SUBROUTINE catcod_mpc                           !
  ! It finds the MPC astrometric catalogue code     !
  ! used in error model from ADES astCat            !
  ! ----------------------------------------------- !

  SUBROUTINE catcod_mpc(icode,ocode)

    CHARACTER(LEN=8), INTENT(IN)  :: icode ! ADES code for astrometric catalogue (astCat)
    CHARACTER(LEN=2), INTENT(OUT) :: ocode ! output (MPC) code for error model

    ocode = "  "

    SELECT CASE (trim(adjustl(icode)))
    CASE("USNOA1")
       ocode = " a"
    CASE("USNOSA1")
       ocode = " b"
    CASE("USNOA2")
       ocode = " c"
    CASE("USNOSA2")
       ocode = " d"
    CASE("UCAC1")
       ocode = " e"
    CASE("Tyc1")
       ocode = " f"
    CASE("Tyc2")
       ocode = " g"
    CASE("GSC1.0")
       ocode = " h"
    CASE("GSC1.1")
       ocode = " i"
    CASE("GSC1.2")
       ocode = " j"
    CASE("GSC2.2")
       ocode = " k"
    CASE("ACT")
       ocode = " l"
    CASE("GSCACT")
       ocode = " m"
    CASE("SDSS8")
       ocode = " n"
    CASE("USNOB1")
       ocode = " o"
    CASE("PPM")
       ocode = " p"
    CASE("UCAC4")
       ocode = " q"
    CASE("UCAC2")
       ocode = " r"
    CASE("USNOB2")  
       ocode = " s"
    CASE("PPMXL")
       ocode = " t"
    CASE("UCAC3")
       ocode = " u"
    CASE("NOMAD")
       ocode = " v"
    CASE("CMC14")
       ocode = " w"
    CASE("Hip2")
       ocode = " x"
    CASE("Hip1")
       ocode = " y"
    CASE("GSC")
       ocode = " z"
    CASE("AC")
       ocode = " A"
    CASE("SAO1984")
       ocode = " B"
    CASE("SAO")
       ocode = " C"
    CASE("AGK3")
       ocode = " D"
    CASE("FK4")
       ocode = " E"
    CASE("ACRS")
       ocode = " F"
    CASE("LickGas")
       ocode = " G"
    CASE("Ida93")
       ocode = " H"
    CASE("Perth70")
       ocode = " I"
    CASE("COSMOS")
       ocode = " J"
    CASE("Yale")
       ocode = " K"
    CASE("2MASS")
       ocode = " L"
    CASE("GSC2.3")
       ocode = " M"
    CASE("SDSS7")
       ocode = " N"
    CASE("SSTRC1")
       ocode = " O"
    CASE("MPOSC3")
       ocode = " P"
    CASE("CMC15")
       ocode = " Q"
    CASE("SSTRC4")
       ocode = " R"
    CASE("URAT1")
       ocode = " S"
    CASE("URAT2")
       ocode = " T"
    CASE("Gaia1")
       ocode = " U"
    CASE("Gaia2")
       ocode = " V"
    CASE("UCAC5")
       ocode = " W"
    CASE("UNK")
       ocode = "  "
    END SELECT

  END SUBROUTINE catcod_mpc

  ! ----------------------------------------------------- !
  ! SUBROUTINE techn_mpc                                  !
  ! It finds the MPC techn code from ADES mode            !
  ! and notes (unless it is an occultation observation).  !
  ! ----------------------------------------------------- !

  SUBROUTINE techn_mpc(mode,notes,ocode)

    CHARACTER(LEN=3), INTENT(IN)  :: mode  ! ADES mode
    CHARACTER(LEN=6), INTENT(IN)  :: notes ! ADES notes
    CHARACTER(LEN=1), INTENT(OUT) :: ocode ! output techn

    ocode = " "

    ! normal place
    IF (INDEX(notes,'n') .NE. 0) THEN
       ocode = "N"
    ENDIF

    SELECT CASE (mode)
    CASE("CCD")
       ocode = "C"
       ! worst rms in error model since we can't map it
    CASE("VID")
       ocode = "n"
    CASE("PHO")
       ocode = "P"
    CASE("ENC")
       ocode = "e"
       ! worst rms in error model since we can't map it
    CASE("PMT")
       !ocode = "M"
       ocode = " "
    CASE("MIC")
       ocode = "M"
    CASE("MER")
       ocode = "T"
    CASE("TDI")
       ocode = "C"
       !CASE("OCC")
       !ocode = "E"
    CASE("UNK")
       ocode = " "
    END SELECT

    IF (ocode .EQ. " ") THEN
       WRITE(ierrou,*) 'techn_mpc: unable to convert to techn for input ADES file '//trim(ades_filename)// &
            & " Input mode and notes: "//trim(adjustl(mode))//" , "//trim(adjustl(notes))
       numerr=numerr+1 
       !RETURN
    ENDIF

  END SUBROUTINE techn_mpc

  ! ----------------------------------------------------- !
  ! SUBROUTINE rot_IAU_to_ECLJ2000                        !
  ! It computes the rotation matrix for conversion from   !
  ! IAU Earth BF to ECLJ2000 reference frames             !
  ! ----------------------------------------------------- !

  SUBROUTINE rot_IAU_to_ECLJ2000(mjd,rot_mat)

    REAL(KIND=dkind),                 INTENT(IN)  :: mjd             ! modified julian date of observation
    REAL(KIND=dkind), DIMENSION(3,3), INTENT(OUT) :: rot_mat         ! rotation matrix

    REAL(KIND=dkind)                 :: dt,ww,asc_ret,decli
    REAL(KIND=dkind), DIMENSION(3,3) :: rot,r1,r3

    ! find difference in days between mjd and modified julian date of J2000 standard epoch (51544.5)
    dt = mjd - 51544.5d0

    ww=(190.147d0 + 360.9856235d0*dt)*pig/180.d0
    dt=dt/36525.d0                 ! in centuries
    asc_ret=(0.d0 - 0.641d0*dt)*pig/180.d0
    decli=(90.d0 - 0.557d0*dt)*pig/180.d0
    ! computation of rotation matrix from ECLJ2000 to Earth BF
    CALL rotmt(ww,r3,3)
    CALL rotmt(pig/2.d0-decli,r1,1)
    rot=MATMUL(r3,r1)
    CALL rotmt(pig/2.d0+asc_ret,r3,3)
    rot=MATMUL(rot,r3)
    rot=MATMUL(rot,TRANSPOSE(roteqec))           ! rotation matrix from ECLJ2000 to Earth BF
    rot_mat=TRANSPOSE(rot)                       ! rotation matrix from BF to ECLJ2000

  END SUBROUTINE rot_IAU_to_ECLJ2000

  ! -------------------------------------------------------------------- !
  ! SUBROUTINE geoIAU_to_cart                                            !
  ! It converts geocentric geodetic coordinates to cartesian following   !
  ! Earth size and shape parameters defined in IAU frame standard        !
  ! (see note 11 in                                                      !
  ! http://aa.usno.navy.mil/publications/reports/Archinaletal2011a.pdf)  !
  ! -------------------------------------------------------------------- !

  SUBROUTINE geoIAU_to_cart(geo_pos,cart_pos)

    REAL(KIND=dkind), INTENT(IN)  :: geo_pos(3)   ! input geodetic 3D position ([long]=deg,[lat]=deg,[alt]=m)
    REAL(KIND=dkind), INTENT(OUT) :: cart_pos(3)  ! output cartesian 3D position (in km)

    REAL(KIND=dkind) :: f,a,b,e_sq,N

    ! Earth size and shape parameters (see description)
    f = 1.d0/298.25642
    a = 6378136.6*1.d-3 ! -> km
    b = a - a*f
    e_sq = (a**2 -b**2)/a**2

    N = a/SQRT(1.d0 - e_sq * (SIN(geo_pos(2)*radeg)**2))

    cart_pos(1) = (N + geo_pos(3)*1.d-3) * COS(geo_pos(2)*radeg) * COS(geo_pos(1)*radeg)
    cart_pos(2) = (N + geo_pos(3)*1.d-3) * COS(geo_pos(2)*radeg) * SIN(geo_pos(1)*radeg)
    cart_pos(3) = ((1.d0 - e_sq)*N + geo_pos(3)*1.d-3) * SIN(geo_pos(2)*radeg)

  END SUBROUTINE geoIAU_to_cart

  ! -------------------------------------------------------------------- !
  ! SUBROUTINE RADec_accuracy                                            !
  ! It computes the accuracies (and precisions of rms) of RA and Dec     !
  ! entries.                                                             !
  ! -------------------------------------------------------------------- !

  SUBROUTINE RADec_accuracy(i,dec,ra_acc,dec_acc,prec)

    INTEGER,          INTENT(IN)    :: i            ! observation index in ades_obs and ades_flg arrays
    REAL(KIND=dkind), INTENT(IN)    :: dec          ! declination in radians
    REAL(KIND=dkind), INTENT(OUT)   :: ra_acc       ! RA accuracy in radians
    REAL(KIND=dkind), INTENT(OUT)   :: dec_acc      ! Dec accuracy in radians
    INTEGER,          INTENT(INOUT) :: prec(2)      ! output precisions of RA, Dec 

    REAL(KIND=dkind), DIMENSION(2,2) :: DP_cov_mat,RD_cov_mat

    ! for RA-Dec obs
    IF (ades_flg(i)%ra .OR. ades_flg(i)%deltaRA) THEN

       ! initialization, for dist and pa we get precisions in input (see calling routine above)
       prec = 0

       ! compute RA and Dec accuracy from ra and dec entries and "convert" them to radians
       IF (ades_flg(i)%opt) THEN
          ra_acc = 10.d0**(-ades_prec(i)%ra*1.d0)*radeg/COS(dec)
          ! this is a precision for sigRA, so include conversion from degree to arcsec 3600 -> prec ~ -3 + 1 for safety
          prec(1) = MAX(0,ades_prec(i)%ra -3 + 1)
          dec_acc = 10.d0**(-ades_prec(i)%dec*1.d0)*radeg
          ! this is a precision for sigDec, so include conversion from degree to arcsec 3600 -> prec ~ -3 + 1 for safety
          prec(2) = MAX(0,ades_prec(i)%dec -3 + 1)
       ELSEIF (ades_flg(i)%occ) THEN
          ra_acc = 10.d0**(-ades_prec(i)%deltaRA*1.d0)/(secrad*COS(dec))
          prec(1) = ades_prec(i)%deltaRA
          dec_acc = 10.d0**(-ades_prec(i)%deltaDec*1.d0)/secrad
          prec(2) = ades_prec(i)%deltaDec
       ENDIF

       ! if precRA and precDec are present, then compare them with ra_acc and dec_acc and take the maximum (worst)
       IF (ades_flg(i)%precRA) THEN
          IF (.NOT. ALL(precRADec_values .NE. ades_obs(i)%precRA)) &
               & ra_acc = MAX(ra_acc,ades_obs(i)%precRA*15.d0/(secrad*COS(dec)))
          IF (.NOT. ALL(precRADec_values .NE. ades_obs(i)%precDec)) &
               & dec_acc = MAX(dec_acc,ades_obs(i)%precDec/secrad)
          ! update prec if ra_acc/dec_acc have been updated
          ! (including conversion for precRA sec to arcsec 15 -> prec ~ -1 + 1 for safety = 0)
          IF (ra_acc .EQ. ades_obs(i)%precRA*15.d0/(secrad*COS(dec))) prec(1) = ades_prec(i)%precRA
          IF (dec_acc .EQ. ades_obs(i)%precDec/secrad) prec(2) = ades_prec(i)%precDec
       ENDIF

       ! for dist-PA obs
    ELSEIF (ades_flg(i)%dist) THEN

       ! get RA-Dec covariance matrix from dist and pa accuracies
       DP_cov_mat(1,1) = (10.d0**(-ades_prec(i)%dist*1.d0)/secrad)**2
       DP_cov_mat(2,2) = (10.d0**(-ades_prec(i)%pa*1.d0)*radeg)**2
       DP_cov_mat(1,2) = 0.d0
       DP_cov_mat(2,1) = 0.d0
       CALL distPA_to_RADec_cov_mat(DP_cov_mat,ades_obs(i)%dist/secrad,ades_obs(i)%pa*radeg, &
            &                                  ades_obs(i)%decStar*radeg,RD_cov_mat)
       ! find ra and dec accuracies
       ra_acc = SQRT(RD_cov_mat(1,1))
       dec_acc = SQRT(RD_cov_mat(2,2))

       ! if precRA and precDec are present, then compare them with ra_acc and dec_acc and take the maximum (worst)
       IF (ades_flg(i)%precRA) THEN
          IF (.NOT. ALL(precRADec_values .NE. ades_obs(i)%precRA)) &
               & ra_acc = MAX(ra_acc,ades_obs(i)%precRA*15.d0/(secrad*COS(dec)))
          IF (.NOT. ALL(precRADec_values .NE. ades_obs(i)%precDec)) &
               & dec_acc = MAX(dec_acc,ades_obs(i)%precDec/secrad)
          ! update prec if ra_acc/dec_acc have been updated
          ! (including conversion for precRA sec to arcsec 15 -> prec ~ -1 + 1 for safety = 0)
          IF (ra_acc .EQ. ades_obs(i)%precRA*15.d0/(secrad*COS(dec))) prec(1) = ades_prec(i)%precRA
          IF (dec_acc .EQ. ades_obs(i)%precDec/secrad) prec(2) = ades_prec(i)%precDec
       ENDIF

    ENDIF

  END SUBROUTINE RADec_accuracy

  ! ------------------------------------------------------------- !
  ! FUNCTION get_prec                                             !
  ! It computes the precision (number of digits after dot) of a   !
  ! "numeric" string.                                             !
  ! ------------------------------------------------------------- !

  INTEGER FUNCTION get_prec(string)

    CHARACTER(LEN=*),INTENT(IN) :: string ! input string
    CHARACTER(LEN=1) :: last_char
    INTEGER :: len_str,i

    LOGICAL, EXTERNAL  :: isnum

    ! check that the string is "numeric"
    IF (.NOT. isreal(trim(adjustl(string))) .AND. .NOT. isnum(trim(adjustl(string)))) THEN
       WRITE(ierrou,*) 'prec_str: Input string is not "numeric" for input ADES file '//trim(ades_filename)// &
            & ' Input string: '//trim(adjustl(string))
       numerr=numerr+1 
       RETURN
    ENDIF

    get_prec = 0
    CALL rmsp(string,len_str)

    IF (INDEX(string,'.') .NE. 0) THEN

       get_prec = len_str - INDEX(string,'.')

    ENDIF

  END FUNCTION get_prec

  ! ------------------------------------------------------------------------------------------------ !
  ! SUBROUTINE distPA_to_RADec_cov_mat                                                               !
  ! It converts the distance-position angle covariance matrix to RA-Dec covariance matrix            !
  ! ------------------------------------------------------------------------------------------------ !
  ! Method:                                                                                          !
  ! D = dist, PA = pos. angle, Y = (RA,Dec), X = (D,PA)                                              !
  ! A = sinD*sinPA, B = sinDec*cosD, C = cosD*cosPA                                                  !
  ! (see also conversion in convert_ades_fitobs)                                                     !
  ! Y = F(X), with F1: RA = asin(A/sqrt(1-(B+C)**2)), F2: Dec = asin(B+C) - Dec_ref                  !
  ! cov_mat_Y = CM * cov_mat_X * CM^T, with CM = dF/dX (partial derivatives)                         !
  ! ------------------------------------------------------------------------------------------------ !

  SUBROUTINE distPA_to_RADec_cov_mat(DP_cov_mat,dist,PA,delta0,RD_cov_mat)

    REAL(KIND=dkind), DIMENSION(2,2), INTENT(IN)  :: DP_cov_mat  ! input distance-position angle covariance matrix
    REAL(KIND=dkind), INTENT(IN)                  :: dist,PA     ! input distance-position coordinates in radians
    REAL(KIND=dkind), INTENT(IN)                  :: delta0      ! input Dec of reference point in radians
    REAL(KIND=dkind), DIMENSION(2,2), INTENT(OUT) :: RD_cov_mat  ! output RA-Dec covariance matrix

    REAL(KIND=dkind), DIMENSION(2,2)  :: CM
    REAL(KIND=dkind)                  :: A,B,C,D,E ! these are not the A,B,... defined in the routine description above

    ! CM11 = dF1/dD
    A = SIN(PA)
    B = SIN(delta0)**2+COS(PA)**2+2.d0*SIN(delta0)*COS(PA)
    C = 1.d0 - B*COS(dist)**2
    CM(1,1) = (A*COS(dist) - (A*B*SIN(dist)**2*COS(dist))/C) / SQRT(C - A**2*SIN(dist)**2)

    ! CM12 = dF1/dPA
    A = SIN(dist)
    B = SIN(delta0)**2*COS(dist)**2
    C = COS(dist)**2
    D = 2.d0*SIN(delta0)*COS(dist)**2
    E = 1.d0 - B - C*COS(PA)**2 - D*COS(PA)
    CM(1,2) = (A*COS(PA) - (A*SIN(PA)**2*(2.d0*C*COS(PA)+D))/(2.d0*E)) / SQRT(E - A**2*SIN(PA)**2)

    ! CM21 = dF2/dD
    CM(2,1) = -(SIN(dist)*(SIN(delta0)+COS(PA))) / SQRT(1.d0 - COS(dist)**2*(SIN(delta0)+COS(PA))**2)

    ! CM22 = dF2/dPA
    CM(2,2) = -(SIN(PA)*COS(dist)) / SQRT(1.d0 - (COS(PA)*COS(dist)+SIN(delta0)*COS(dist))**2)

    ! find RA-Dec cov mat via cov_mat_Y = CM * cov_mat_X * CM^T 
    CALL convertcov(2,DP_cov_mat,CM,RD_cov_mat)

  END SUBROUTINE distPA_to_RADec_cov_mat

  ! ------------------------------------------------- !
  ! SUBROUTINE precision_group_control                !
  ! It controls that values in ADES fields precRA and !
  ! precDec follow the required rules.                !
  ! ------------------------------------------------- !

  SUBROUTINE precision_group_control(precision,string)
    REAL(KIND=dkind),  INTENT(INOUT)  :: precision
    CHARACTER(LEN=40), INTENT(OUT)    :: string
    INTEGER                           :: l,indx

    ! using required default psv formatting
    WRITE(string,'(SS,F6.3)') precision
    CALL rmsp(string,l)
    indx = INDEX(string,'.')
    IF (indx .NE. 0) THEN
       ! non-integer values are printed with no trailing zeros
       DO WHILE (string(l:l) .EQ. '0')
          string = string(1:l-1)
          l = l-1
       ENDDO
       IF  (string(l:l) .EQ. '.') THEN
          ! integer values are printed with no decimal point
          string = string(1:l-1)  
       ENDIF
    ENDIF

    READ(string,*) precision

  END SUBROUTINE precision_group_control

  ! --------------------------------------------------- !
  ! SUBROUTINE fill_psv_rec                             !
  ! Stores arrays and variables to be used for writing  !
  ! a whole .psv file when it is not given as input.    !
  ! --------------------------------------------------- !

  SUBROUTINE  fill_psv_rec(new_header,rms)

    LOGICAL, INTENT(IN) :: new_header  ! flag for writing new obsContext
    LOGICAL, INTENT(IN) :: rms         ! TRUE if there is global RMS of astrometric/photometric fit

    INTEGER            :: obs_line
    INTEGER            :: i,j

    ! take into account the number of ObsContext lines 
    ! if data have not changed (but re-sorted), write the same ObsContext Block and header_length is already stored
    ! if data have changed (or if there is no obsContext block in input file), define new ObsContext and header_length
    IF(new_header .OR. header_length .LE. 1) header_length = 3  ! define number of ObsContext lines: version line + 1 comment line + 1 fields line
    IF (rms) header_length = header_length + 1  ! add rms line in header
    ades_key_rec = ""
    ades_data_rec = ""
    ntot_to_nobs_idx = 0
    nobs_to_ntot_idx = 0
    new_key_rec(1) = .TRUE.
    ntot_ades = header_length + nobs_ades + 1  ! consider 1 line for each observation (values lines)
    obs_line = ntot_ades - nobs_ades + 1  ! line number corresponding to first line of observational data
    ntot_to_nobs_idx(obs_line) = 1
    nobs_to_ntot_idx(1) = obs_line
    ! add a fields line if type of observation has changed wrt previous observation
    ! define arrays ntot_to_nobs_idx and nobs_to_ntot_idx
    DO j = 2,nobs_ades
       IF ((ades_flg(j)%opt .NEQV. ades_flg(j-1)%opt) .OR. (ades_flg(j)%occ .NEQV. ades_flg(j-1)%occ) .OR. &
            & (ades_flg(j)%delay .NEQV. ades_flg(j-1)%delay) .OR. (ades_flg(j)%doppler .NEQV. ades_flg(j-1)%doppler) .OR. &
            & (ades_flg(j)%sys .NEQV. ades_flg(j-1)%sys)) THEN
          ntot_ades = ntot_ades + 1  ! add a line
          new_key_rec(j) = .TRUE.
          obs_line = obs_line + 2  ! obs_line + 1 is a fields line
          ntot_to_nobs_idx(obs_line) = j
          nobs_to_ntot_idx(j) = obs_line
       ELSE
          new_key_rec(j) = .FALSE.
          obs_line = obs_line + 1
          ntot_to_nobs_idx(obs_line) = j
          nobs_to_ntot_idx(j) = obs_line
       ENDIF
    ENDDO

    ! for each observations fill ades_key_rec and ades_data_rec
    DO i = 1,nobs_ades
       CALL addobs_psv_rec(i)
    ENDDO

  END SUBROUTINE fill_psv_rec

  ! --------------------------------------------- !
  ! SUBROUTINE addobs_psv_rec                     !
  ! For a given observation entry i in ADES data  ! 
  ! type, stores the corresponding line in arrays !
  ! ades_data_rec and ades_key_rec, recording     !
  ! lines to be written on .psv file.             !
  ! --------------------------------------------- !

  SUBROUTINE addobs_psv_rec(i)

    INTEGER, INTENT(IN) :: i        ! index of observation

    INTEGER             :: lin      ! line number of the observation
    INTEGER             :: tot_len,prec,len_rec
    INTEGER             :: len_field(36)
    CHARACTER(LEN=100)  :: aux_str1,aux_str2,data
    CHARACTER(LEN=20)   :: FMT,prec_str,totlen_str,FMT_1,sfmt(36)

    ! variables initialization
    lin = nobs_to_ntot_idx(i)
    len_field = 0
    len_rec = 0

    ! store lengths and formats
    ! Identification Group
    IF (ades_flg(i)%permID) THEN
       len_field(1) = 8 !MAX(len("|permID"),len(trim(adjustl(ades_obs(i)%permID)))+1) 
    ELSEIF (ades_flg(i)%provID) THEN
       len_field(1) = 12 !MAX(len("|provID"),len(trim(adjustl(ades_obs(i)%provID)))+1) 
    ELSEIF (ades_flg(i)%trkSub) THEN 
       len_field(1) = 9 !MAX(len("|trkSub"),len(trim(adjustl(ades_obs(i)%trkSub)))+1) 
    ENDIF
    ! obsID
    !len_field(2) = 20  ! according to ADES formatting rules
    !trkID
    ! IF (.NOT. ades_flg(i)%rad) THEN
    !   len_field(3) = 13
    !ENDIF
    ! mode
    IF (.NOT. ades_flg(i)%rad .AND. .NOT. ades_flg(i)%occ) THEN
       len_field(2) = 5   ! according to ADES formatting rules
    ENDIF
    ! stn
    IF (.NOT. ades_flg(i)%rad) THEN
       len_field(3) = 5   ! according to ADES formatting rules
    ENDIF
    ! trx and rcv
    IF (ades_flg(i)%rad) THEN
       len_field(4) = 5   ! according to ADES formatting rules
       len_field(5) = 5
    ENDIF
    ! Location Group
    ! must not be present if stn is associated with a stationary MPC observatory code
    IF (ades_flg(i)%sys) THEN
       ! sys
       len_field(6) = 8 !MAX(len("|sys"),len(trim(adjustl(ades_obs(i)%sys)))+1)
       ! ctr
       WRITE(data,"(i20)") ades_obs(i)%ctr
       len_field(7) = 4 !MAX(len("|ctr"),len(trim(adjustl(data)))+1) 
       ! pos
       prec = 4
       tot_len = 6 + prec + 2 ! +2 to account for decimal point and the sign
       WRITE(totlen_str,"(i20)") tot_len
       WRITE(prec_str,"(i20)") prec
       sfmt(8) = "(F"//trim(adjustl(totlen_str))//"."//trim(adjustl(prec_str))//")"
       len_field(8) = tot_len+1
       ! posCov
       !prec = 4
       !tot_len = 6 + prec + 2 ! +2 to account for decimal point and the sign
       !WRITE(totlen_str,"(i20)") tot_len
       !WRITE(prec_str,"(i20)") prec
       !sfmt(9) = "(F"//trim(adjustl(totlen_str))//"."//trim(adjustl(prec_str))//")"
       !len_field(9) = tot_len+1
    ENDIF

    ! type of observation
    IF (ades_flg(i)%opt) THEN
       len_field(10) = 7 
    ELSEIF (ades_flg(i)%off) THEN
       len_field(10) = 6 
    ELSEIF (ades_flg(i)%occ) THEN 
       len_field(10) = 11 
    ELSEIF (ades_flg(i)%rad) THEN
       len_field(10) = 5
    ENDIF

    ! obsTime (always required)
    len_field(11) = 25  ! according to ADES formatting rules

    ! Observation Group
    IF (ades_flg(i)%opt) THEN   
       ! ra
       !prec = ades_prec(i)%ra  
       !prec = 7                 ! using required default psv formatting width tot_len = 11
       prec = 9
       tot_len = 3 + prec + 1   ! +1 for '.' 
       WRITE(totlen_str,"(i20)") tot_len
       WRITE(prec_str,"(i20)") prec
       sfmt(12) = "(F"//trim(adjustl(totlen_str))//"."//trim(adjustl(prec_str))//")"
       len_field(12) = tot_len+1  ! +1 for |
       ! dec
       !prec = ades_prec(i)%dec
       !prec = 7                 ! using required default psv formatting width tot_len =11
       prec = 9
       tot_len = 3 + prec + 1   ! +1 for '.' 
       WRITE(totlen_str,"(i20)") tot_len
       WRITE(prec_str,"(i20)") prec
       sfmt(13) = "(F"//trim(adjustl(totlen_str))//"."//trim(adjustl(prec_str))//")"
       len_field(13) = tot_len+1  ! +1 for |
    ELSEIF (ades_flg(i)%occ) THEN
       ! raStar
       !prec = ades_prec(i)%raStar
       !prec = 7                 ! using required default psv formatting width tot_len = 11
       prec = 9
       tot_len = 3 + prec + 1   ! +1 for '.' 
       WRITE(totlen_str,"(i20)") tot_len
       WRITE(prec_str,"(i20)") prec
       sfmt(14) = "(F"//trim(adjustl(totlen_str))//"."//trim(adjustl(prec_str))//")"
       len_field(14) = tot_len+1
       ! decStar
       !prec = ades_prec(i)%decStar
       !prec = 7                 ! using required default psv formatting width tot_len = 11
       prec = 9
       tot_len = 3 + prec + 1   ! +1 for '.' 
       WRITE(totlen_str,"(i20)") tot_len
       WRITE(prec_str,"(i20)") prec
       sfmt(15) = "(F"//trim(adjustl(totlen_str))//"."//trim(adjustl(prec_str))//")"
       len_field(15) = tot_len+1 
       ! deltaRA
       prec = 3
       tot_len = 7 + prec + 1
       WRITE(totlen_str,"(i20)") tot_len
       WRITE(prec_str,"(i20)") prec
       sfmt(16) = "(F"//trim(adjustl(totlen_str))//"."//trim(adjustl(prec_str))//")"
       len_field(16) = MAX(len("|deltaRA"),tot_len+1) 
       ! deltaDec
       prec = 3
       tot_len = 7 + prec + 1
       WRITE(totlen_str,"(i20)") tot_len
       WRITE(prec_str,"(i20)") prec
       sfmt(17) = "(F"//trim(adjustl(totlen_str))//"."//trim(adjustl(prec_str))//")"
       len_field(17) = MAX(len("|deltaDec"),tot_len+1)
    ELSEIF (ades_flg(i)%delay) THEN
       ! delay
       !prec = ades_prec(i)%delay
       prec = 8                
       tot_len = 4 + prec + 1  ! +1 for "."
       WRITE(totlen_str,"(i20)") tot_len
       WRITE(prec_str,"(i20)") prec
       sfmt(18) = "(SS,F"//trim(adjustl(totlen_str))//"."//trim(adjustl(prec_str))//")"
       len_field(18) = MAX(len("|delay"),tot_len+1) 
       ! rmsDelay
       !prec = ades_prec(i)%rmsDelay
       prec = 3
       tot_len = 2 + prec + 1  ! +1 for "."
       WRITE(totlen_str,"(i20)") tot_len
       WRITE(prec_str,"(i20)") prec
       sfmt(19) = "(SS,F"//trim(adjustl(totlen_str))//"."//trim(adjustl(prec_str))//")"
       len_field(19) = MAX(len("|rmsDelay"),tot_len+1) 
    ELSEIF (ades_flg(i)%doppler) THEN
       ! doppler
       !prec = ades_prec(i)%doppler
       prec = 2
       tot_len = 7 + prec + 2  ! +2 for sign and "."
       WRITE(totlen_str,"(i20)") tot_len
       WRITE(prec_str,"(i20)") prec
       sfmt(20) = "(F"//trim(adjustl(totlen_str))//"."//trim(adjustl(prec_str))//")"
       len_field(20) = MAX(len("|doppler"),tot_len+1) 
       ! rmsDoppler
       !prec = ades_prec(i)%rmsDoppler
       prec = 3
       tot_len = 2 + prec + 1  ! +1 for "."
       WRITE(totlen_str,"(i20)") tot_len
       WRITE(prec_str,"(i20)") prec
       sfmt(21) = "(SS,F"//trim(adjustl(totlen_str))//"."//trim(adjustl(prec_str))//")"
       len_field(21) = MAX(len("|rmsDoppler"),tot_len+1) 
    ENDIF
    ! astCat
    IF (ades_flg(i)%opt .OR. ades_flg(i)%occ) THEN
       len_field(22) = 9 ! 8 + 1 to account for |
    ENDIF
    ! Photometry Group
    IF (.NOT. ades_flg(i)%rad) THEN
       ! mag
       IF (ades_obs(i)%mag .NE. 9.9d9) THEN
          ! prec = ades_prec(i)%mag  ! required default psv formatting width = 5
          prec = 2
          tot_len = 5  ! 2 + prec + 1 for '.'
          WRITE(totlen_str,"(i20)") tot_len
          WRITE(prec_str,"(i20)") prec
          sfmt(23) = "(F"//trim(adjustl(totlen_str))//"."//trim(adjustl(prec_str))//")"
       ENDIF
       len_field(23) = 6 ! + 1 to account for |
       ! band
       len_field(24) = 5 ! MAX(3,len(key)) + 1 to account for |
       ! photCat
       len_field(25) = 9 ! 8 + 1 to account for |
    ENDIF
    ! com and frq
    IF (ades_flg(i)%rad) THEN
       ! com
       len_field(26) = 4 ! LEN(key) + 1 to account for | 
       ! frq
       prec = 1 !prec = ades_prec(i)%frq
       tot_len = 4 + prec + 1  ! +1 for "."   ! as read from .rad file
       WRITE(totlen_str,"(i20)") tot_len
       WRITE(prec_str,"(i20)") prec
       sfmt(27) = "(SS,F"//trim(adjustl(totlen_str))//"."//trim(adjustl(prec_str))//")"
       len_field(27) = MAX(len("|frq"),tot_len+1)
    ENDIF

    !prog
    len_field(28) = 5  ! it is a 2-characters string, thus we take length of |prog

    !ref
    !len_field(28) = 6  !MAX(len("|ref"),len(trim(adjustl(ades_obs(i)%ref)))+1)

    ! disc, subFmt, Precision Group and notes
    IF (.NOT. ades_flg(i)%rad)  THEN
       ! disc
       len_field(29) = 5 ! LEN(key) + 1 to account for |

       ! subFmt
       len_field(30) = 7 !LEN(key) + 1 to account for | (subFmt is up to 4 aplhanumeric characters)

       ! Precision Group
       ! precTime
       WRITE(data,"(i20)")  ades_obs(i)%precTime
       len_field(31) = MAX(len("|precTime"),len(trim(adjustl(data)))+1)

       ! precRA
       prec = 3    ! required default psv formatting
       tot_len = 2 + prec + 1 ! 1 for '.'
       !WRITE(totlen_str,"(i20)") tot_len
       !WRITE(prec_str,"(i20)") prec
       !sfmt(32) = "(SS,F"//trim(adjustl(totlen_str))//"."//trim(adjustl(prec_str))//")"
       len_field(32) = MAX(len("|precRA"),tot_len+1)

       ! precDec
       prec = 3    ! required default psv formatting
       tot_len = 2 + prec + 1  ! 1 for '.'
       !WRITE(totlen_str,"(i20)") tot_len
       !WRITE(prec_str,"(i20)") prec
       !sfmt(33) = "(SS,F"//trim(adjustl(totlen_str))//"."//trim(adjustl(prec_str))//")"
       len_field(33) = MAX(len("|precDec"),tot_len+1)

       ! notes
       len_field(34) = 7 ! 6 + 1 to account for |
    ENDIF
    ! remarks
    len_field(35) = MAX(len("|remarks"),len(trim(adjustl(ades_obs(i)%remarks)))+1)

    ! deprecated
    IF (.NOT. ades_flg(i)%rad) len_field(36) = 11

    ! fill ades_key_rec
    IF (new_key_rec(i)) THEN
       ! type of observation
       WRITE(FMT,"(i20)") len_field(10)
       WRITE(aux_str1,"(A"//trim(adjustl(FMT))//")") adjustl("type")
       ades_key_rec(lin-1) = adjustl(aux_str1(1:len_field(10)))
       len_rec = len_field(10)
       ! Identification Group
       WRITE(FMT,"(i20)") len_field(1)
       IF (ades_flg(i)%permID) THEN
          WRITE(aux_str1,"(A"//trim(adjustl(FMT))//")") "|permID"
       ELSEIF (ades_flg(i)%provID) THEN
          WRITE(aux_str1,"(A"//trim(adjustl(FMT))//")") "|provID"
       ELSEIF (ades_flg(i)%trkSub) THEN 
          WRITE(aux_str1,"(A"//trim(adjustl(FMT))//")") "|trkSub"
       ENDIF
       ades_key_rec(lin-1) = ades_key_rec(lin-1)(1:len_rec)//adjustl(aux_str1(1:len_field(1)))
       len_rec =  len_rec + len_field(1)
       ! obsID
       !WRITE(FMT,"(i20)") len_field(2)
       !WRITE(aux_str1,"(A"//trim(adjustl(FMT))//")") "|obsID"
       !ades_key_rec(lin-1) = ades_key_rec(lin-1)(1:len_rec)//adjustl(aux_str1(1:len_field(2)))
       !len_rec = len_rec + len_field(2)
       !trkID
       !IF (.NOT. ades_flg(i)%rad) THEN
       !  WRITE(FMT,"(i20)") len_field(3)
       !  WRITE(aux_str1,"(A"//trim(adjustl(FMT))//")") "|trkID"
       !  ades_key_rec(lin-1) = ades_key_rec(lin-1)(1:len_rec)//adjustl(aux_str1(1:len_field(3)))
       !  len_rec = len_rec + len_field(3)
       !ENDIF
       ! mode
       IF (.NOT. ades_flg(i)%rad .AND. .NOT. ades_flg(i)%occ) THEN
          WRITE(FMT,"(i20)") len_field(2)
          WRITE(aux_str1,"(A"//trim(adjustl(FMT))//")") "|mode"
          ades_key_rec(lin-1) = ades_key_rec(lin-1)(1:len_rec)//adjustl(aux_str1(1:len_field(2)))
          len_rec = len_rec + len_field(2)
       ENDIF
       ! stn
       IF (.NOT. ades_flg(i)%rad) THEN
          WRITE(FMT,"(i20)") len_field(3)
          WRITE(aux_str1,"(A"//trim(adjustl(FMT))//")") "|stn"
          ades_key_rec(lin-1) = ades_key_rec(lin-1)(1:len_rec)//adjustl(aux_str1(1:len_field(3)))
          len_rec = len_rec + len_field(3)
       ENDIF
       ! trx and rcv
       IF (ades_flg(i)%rad) THEN
          WRITE(FMT,"(i20)") len_field(4)
          WRITE(aux_str1,"(A"//trim(adjustl(FMT))//")") "|trx"
          ades_key_rec(lin-1) = ades_key_rec(lin-1)(1:len_rec)//adjustl(aux_str1(1:len_field(4)))
          len_rec = len_rec + len_field(4)
          WRITE(FMT,"(i20)") len_field(5)
          WRITE(aux_str1,"(A"//trim(adjustl(FMT))//")") "|rcv"
          ades_key_rec(lin-1) = ades_key_rec(lin-1)(1:len_rec)//adjustl(aux_str1(1:len_field(5)))
          len_rec = len_rec + len_field(5)
       ENDIF
       ! Location Group
       ! must not be present if stn is associated with a stationary MPC observatory code
       IF (ades_flg(i)%sys) THEN
          ! sys
          WRITE(FMT,"(i20)") len_field(6)
          WRITE(aux_str1,"(A"//trim(adjustl(FMT))//")") "|sys"
          ades_key_rec(lin-1) = ades_key_rec(lin-1)(1:len_rec)//adjustl(aux_str1(1:len_field(6)))
          len_rec = len_rec + len_field(6)
          ! ctr
          WRITE(FMT,"(i20)") len_field(7)
          WRITE(aux_str1,"(A"//trim(adjustl(FMT))//")") "|ctr"
          ades_key_rec(lin-1) = ades_key_rec(lin-1)(1:len_rec)//adjustl(aux_str1(1:len_field(7)))
          len_rec = len_rec + len_field(7)
          ! pos1
          WRITE(FMT,"(i20)") len_field(8)
          WRITE(aux_str1,"(A"//trim(adjustl(FMT))//")") "|pos1"
          ades_key_rec(lin-1) = ades_key_rec(lin-1)(1:len_rec)//adjustl(aux_str1(1:len_field(8)))
          len_rec = len_rec + len_field(8)
          ! pos2
          WRITE(FMT,"(i20)") len_field(8)
          WRITE(aux_str1,"(A"//trim(adjustl(FMT))//")") "|pos2"
          ades_key_rec(lin-1) = ades_key_rec(lin-1)(1:len_rec)//adjustl(aux_str1(1:len_field(8)))
          len_rec = len_rec + len_field(8)
          ! pos3
          WRITE(FMT,"(i20)") len_field(8)
          WRITE(aux_str1,"(A"//trim(adjustl(FMT))//")") "|pos3"
          ades_key_rec(lin-1) = ades_key_rec(lin-1)(1:len_rec)//adjustl(aux_str1(1:len_field(8)))
          len_rec = len_rec + len_field(8)
          ! posCov11
          !WRITE(FMT,"(i20)") len_field(9)
          !WRITE(aux_str1,"(A"//trim(adjustl(FMT))//")") "|posCov11"
          !ades_key_rec(lin-1) = ades_key_rec(lin-1)(1:len_rec)//adjustl(aux_str1(1:len_field(9)))
          !len_rec = len_rec + len_field(9)
          ! posCov12
          !WRITE(FMT,"(i20)") len_field(9)
          !WRITE(aux_str1,"(A"//trim(adjustl(FMT))//")") "|posCov12"
          !ades_key_rec(lin-1) = ades_key_rec(lin-1)(1:len_rec)//adjustl(aux_str1(1:len_field(9)))
          !len_rec = len_rec + len_field(9)
          ! posCov13
          !WRITE(FMT,"(i20)") len_field(9)
          !WRITE(aux_str1,"(A"//trim(adjustl(FMT))//")") "|posCov13"
          !ades_key_rec(lin-1) = ades_key_rec(lin-1)(1:len_rec)//adjustl(aux_str1(1:len_field(9)))
          !len_rec = len_rec + len_field(9)
          ! posCov22
          !WRITE(FMT,"(i20)") len_field(9)
          !WRITE(aux_str1,"(A"//trim(adjustl(FMT))//")") "|posCov22"
          !ades_key_rec(lin-1) = ades_key_rec(lin-1)(1:len_rec)//adjustl(aux_str1(1:len_field(9)))
          !len_rec = len_rec + len_field(9)
          ! posCov23
          !WRITE(FMT,"(i20)") len_field(9)
          !WRITE(aux_str1,"(A"//trim(adjustl(FMT))//")") "|posCov23"
          !ades_key_rec(lin-1) = ades_key_rec(lin-1)(1:len_rec)//adjustl(aux_str1(1:len_field(9)))
          !len_rec = len_rec + len_field(9)
          ! posCov33
          !WRITE(FMT,"(i20)") len_field(9)
          !WRITE(aux_str1,"(A"//trim(adjustl(FMT))//")") "|posCov33"
          !ades_key_rec(lin-1) = ades_key_rec(lin-1)(1:len_rec)//adjustl(aux_str1(1:len_field(9)))
          !len_rec = len_rec + len_field(9)
       ENDIF

       ! prog
       WRITE(FMT,"(i20)") len_field(28)
       WRITE(aux_str1,"(A"//trim(adjustl(FMT))//")") "|prog"
       ades_key_rec(lin-1) = ades_key_rec(lin-1)(1:len_rec)//adjustl(aux_str1(1:len_field(28)))
       len_rec = len_rec + len_field(28)

       ! obsTime (always required)
       WRITE(FMT,"(i20)") len_field(11)
       WRITE(aux_str1,"(A"//trim(adjustl(FMT))//")") "|obsTime"
       ades_key_rec(lin-1) = ades_key_rec(lin-1)(1:len_rec)//adjustl(aux_str1(1:len_field(11)))
       len_rec = len_rec + len_field(11)

       ! Observation Group
       IF (ades_flg(i)%opt) THEN
          ! ra
          WRITE(FMT,"(i20)") len_field(12)
          WRITE(aux_str1,"(A"//trim(adjustl(FMT))//")") "|ra"
          ades_key_rec(lin-1) = ades_key_rec(lin-1)(1:len_rec)//adjustl(aux_str1(1:len_field(12)))
          len_rec = len_rec + len_field(12)
          ! dec
          WRITE(FMT,"(i20)") len_field(13)
          WRITE(aux_str1,"(A"//trim(adjustl(FMT))//")") "|dec"
          ades_key_rec(lin-1) = ades_key_rec(lin-1)(1:len_rec)//adjustl(aux_str1(1:len_field(13)))
          len_rec = len_rec + len_field(13)
       ELSEIF (ades_flg(i)%occ) THEN
          ! raStar
          WRITE(FMT,"(i20)") len_field(14)
          WRITE(aux_str1,"(A"//trim(adjustl(FMT))//")") "|raStar"
          ades_key_rec(lin-1) = ades_key_rec(lin-1)(1:len_rec)//adjustl(aux_str1(1:len_field(14)))
          len_rec = len_rec + len_field(14)
          ! decStar
          WRITE(FMT,"(i20)") len_field(15)
          WRITE(aux_str1,"(A"//trim(adjustl(FMT))//")") "|decStar"
          ades_key_rec(lin-1) = ades_key_rec(lin-1)(1:len_rec)//adjustl(aux_str1(1:len_field(15)))
          len_rec = len_rec + len_field(15)
          ! deltaRA
          WRITE(FMT,"(i20)") len_field(16)
          WRITE(aux_str1,"(A"//trim(adjustl(FMT))//")") "|deltaRA"
          ades_key_rec(lin-1) = ades_key_rec(lin-1)(1:len_rec)//adjustl(aux_str1(1:len_field(16)))
          len_rec = len_rec + len_field(16)
          ! deltaDec
          WRITE(FMT,"(i20)") len_field(17)
          WRITE(aux_str1,"(A"//trim(adjustl(FMT))//")") "|deltaDec"
          ades_key_rec(lin-1) = ades_key_rec(lin-1)(1:len_rec)//adjustl(aux_str1(1:len_field(17)))
          len_rec = len_rec + len_field(17)
       ELSEIF (ades_flg(i)%delay) THEN
          ! delay
          WRITE(FMT,"(i20)") len_field(18)
          WRITE(aux_str1,"(A"//trim(adjustl(FMT))//")") "|delay"
          ades_key_rec(lin-1) = ades_key_rec(lin-1)(1:len_rec)//adjustl(aux_str1(1:len_field(18)))
          len_rec = len_rec + len_field(18)
          ! rmsDelay
          WRITE(FMT,"(i20)") len_field(19)
          WRITE(aux_str1,"(A"//trim(adjustl(FMT))//")") "|rmsDelay"
          ades_key_rec(lin-1) = ades_key_rec(lin-1)(1:len_rec)//adjustl(aux_str1(1:len_field(19)))
          len_rec = len_rec + len_field(19)
       ELSEIF (ades_flg(i)%doppler) THEN
          ! doppler
          WRITE(FMT,"(i20)") len_field(20)
          WRITE(aux_str1,"(A"//trim(adjustl(FMT))//")") "|doppler"
          ades_key_rec(lin-1) = ades_key_rec(lin-1)(1:len_rec)//adjustl(aux_str1(1:len_field(20)))
          len_rec = len_rec + len_field(20)
          ! rmsDoppler
          WRITE(FMT,"(i20)") len_field(21)
          WRITE(aux_str1,"(A"//trim(adjustl(FMT))//")") "|rmsDoppler"
          ades_key_rec(lin-1) = ades_key_rec(lin-1)(1:len_rec)//adjustl(aux_str1(1:len_field(21)))
          len_rec = len_rec + len_field(21)
       ENDIF

       ! astCat
       IF (ades_flg(i)%opt .OR. ades_flg(i)%occ) THEN
          WRITE(FMT,"(i20)") len_field(22)
          WRITE(aux_str1,"(A"//trim(adjustl(FMT))//")") "|astCat"
          ades_key_rec(lin-1) = ades_key_rec(lin-1)(1:len_rec)//adjustl(aux_str1(1:len_field(22)))
          len_rec = len_rec + len_field(22)
       ENDIF

       ! Photometry Group
       IF (.NOT. ades_flg(i)%rad) THEN
          ! mag
          WRITE(FMT,"(i20)") len_field(23)
          WRITE(aux_str1,"(A"//trim(adjustl(FMT))//")") "|mag"
          ades_key_rec(lin-1) = ades_key_rec(lin-1)(1:len_rec)//adjustl(aux_str1(1:len_field(23)))
          len_rec = len_rec + len_field(23)
          ! band
          WRITE(FMT,"(i20)") len_field(24)
          WRITE(aux_str1,"(A"//trim(adjustl(FMT))//")") "|band"
          ades_key_rec(lin-1) = ades_key_rec(lin-1)(1:len_rec)//adjustl(aux_str1(1:len_field(24)))
          len_rec = len_rec + len_field(24)
          ! photCat
          WRITE(FMT,"(i20)") len_field(25)
          WRITE(aux_str1,"(A"//trim(adjustl(FMT))//")") "|photCat"
          ades_key_rec(lin-1) = ades_key_rec(lin-1)(1:len_rec)//adjustl(aux_str1(1:len_field(25)))
          len_rec = len_rec + len_field(25)
       ENDIF

       ! com and frq
       IF (ades_flg(i)%rad) THEN
          ! com
          WRITE(FMT,"(i20)") len_field(26)
          WRITE(aux_str1,"(A"//trim(adjustl(FMT))//")") "|com"
          ades_key_rec(lin-1) = ades_key_rec(lin-1)(1:len_rec)//adjustl(aux_str1(1:len_field(26)))
          len_rec = len_rec + len_field(26)
          ! frq
          WRITE(FMT,"(i20)") len_field(27)
          WRITE(aux_str1,"(A"//trim(adjustl(FMT))//")") "|frq"
          ades_key_rec(lin-1) = ades_key_rec(lin-1)(1:len_rec)//adjustl(aux_str1(1:len_field(27)))
          len_rec = len_rec + len_field(27)
       ENDIF

       ! ref
       !WRITE(FMT,"(i20)") len_field(28)
       !WRITE(aux_str1,"(A"//trim(adjustl(FMT))//")") "|ref"
       !ades_key_rec(lin-1) = ades_key_rec(lin-1)(1:len_rec)//adjustl(aux_str1(1:len_field(28)))
       !len_rec = len_rec + len_field(28)

       ! disc, Precision Group and notes
       IF (.NOT. ades_flg(i)%rad)  THEN
          ! disc
          WRITE(FMT,"(i20)") len_field(29)
          WRITE(aux_str1,"(A"//trim(adjustl(FMT))//")") "|disc"
          ades_key_rec(lin-1) = ades_key_rec(lin-1)(1:len_rec)//adjustl(aux_str1(1:len_field(29)))
          len_rec = len_rec + len_field(29)

          ! subFmt
          WRITE(FMT,"(i20)") len_field(30)
          WRITE(aux_str1,"(A"//trim(adjustl(FMT))//")") "|subFmt"
          ades_key_rec(lin-1) = ades_key_rec(lin-1)(1:len_rec)//adjustl(aux_str1(1:len_field(30)))
          len_rec = len_rec + len_field(30)

          ! Precision Group
          ! precTime
          WRITE(FMT,"(i20)") len_field(31)
          WRITE(aux_str1,"(A"//trim(adjustl(FMT))//")") "|precTime"
          ades_key_rec(lin-1) = ades_key_rec(lin-1)(1:len_rec)//adjustl(aux_str1(1:len_field(31)))
          len_rec = len_rec + len_field(31)
          ! precRA
          WRITE(FMT,"(i20)") len_field(32)
          WRITE(aux_str1,"(A"//trim(adjustl(FMT))//")") "|precRA"
          ades_key_rec(lin-1) = ades_key_rec(lin-1)(1:len_rec)//adjustl(aux_str1(1:len_field(32)))
          len_rec = len_rec + len_field(32)
          ! precDec
          WRITE(FMT,"(i20)") len_field(33)
          WRITE(aux_str1,"(A"//trim(adjustl(FMT))//")") "|precDec"
          ades_key_rec(lin-1) = ades_key_rec(lin-1)(1:len_rec)//adjustl(aux_str1(1:len_field(33)))
          len_rec = len_rec + len_field(33)

          ! notes
          WRITE(FMT,"(i20)") len_field(34)
          WRITE(aux_str1,"(A"//trim(adjustl(FMT))//")") "|notes"
          ades_key_rec(lin-1) = ades_key_rec(lin-1)(1:len_rec)//adjustl(aux_str1(1:len_field(34)))
          len_rec = len_rec + len_field(34)
       ENDIF

       ! remarks
       WRITE(FMT,"(i20)") len_field(35)
       WRITE(aux_str1,"(A"//trim(adjustl(FMT))//")") "|remarks"
       ades_key_rec(lin-1) = ades_key_rec(lin-1)(1:len_rec)//adjustl(aux_str1(1:len_field(35)))
       len_rec = len_rec + len_field(35)

       ! deprecated
       IF (.NOT. ades_flg(i)%rad) THEN
          WRITE(FMT,"(i20)") len_field(36)
          WRITE(aux_str1,"(A"//trim(adjustl(FMT))//")") "|deprecated"
          ades_key_rec(lin-1) = ades_key_rec(lin-1)(1:len_rec)//adjustl(aux_str1(1:len_field(36)))
          len_rec = len_rec + len_field(36)
       ENDIF
    ELSE
       ades_key_rec(lin-1) = ""
    ENDIF

    ! fill ades_data_rec
    ! type of observation
    WRITE(FMT,"(i20)") len_field(10)
    IF (ades_flg(i)%opt) THEN
       WRITE(aux_str2,"(A"//trim(adjustl(FMT))//")") "optical"
    ELSEIF (ades_flg(i)%off) THEN
       WRITE(aux_str2,"(A"//trim(adjustl(FMT))//")") "offset"
    ELSEIF (ades_flg(i)%occ) THEN 
       WRITE(aux_str2,"(A"//trim(adjustl(FMT))//")") "occultation"
    ELSEIF (ades_flg(i)%rad) THEN
       WRITE(aux_str2,"(A"//trim(adjustl(FMT))//")") "radar"
    ENDIF
    ades_data_rec(lin) = adjustl(aux_str2(1:len_field(10)))
    len_rec = len_field(10)
    ! Identification Group
    WRITE(FMT,"(i20)") len_field(1)
    IF (ades_flg(i)%permID) THEN
       WRITE(aux_str2,"(A"//trim(adjustl(FMT))//")") "|"//trim(adjustl(ades_obs(i)%permID))
       ades_data_rec(lin) =  ades_data_rec(lin)(1:len_rec)//adjustl(aux_str2(1:len_field(1)))
    ELSEIF (ades_flg(i)%provID) THEN
       WRITE(aux_str2,"(A"//trim(adjustl(FMT))//")") "|"//trim(adjustl(ades_obs(i)%provID))
       ades_data_rec(lin) =  ades_data_rec(lin)(1:len_rec)//adjustl(aux_str2(1:len_field(1)))
    ELSEIF (ades_flg(i)%trkSub) THEN
       WRITE(aux_str2,"(A"//trim(adjustl(FMT))//")") "|"//trim(adjustl(ades_obs(i)%trkSub))
       ades_data_rec(lin) =  ades_data_rec(lin)(1:len_rec)//adjustl(aux_str2(1:len_field(1)))
    ENDIF
    len_rec = len_rec + len_field(1)
    ! obsID and trkID left empty
    ! obsID
    !WRITE(FMT,"(i20)") len_field(2)
    !WRITE(aux_str2,"(A"//trim(adjustl(FMT))//")") "|"
    !ades_data_rec(lin) = ades_data_rec(lin)(1:len_rec)//adjustl(aux_str2(1:len_field(2)))
    !len_rec = len_rec + len_field(2)
    !trkID
    !IF (.NOT. ades_flg(i)%rad) THEN
    !   WRITE(FMT,"(i20)") len_field(3)
    !   WRITE(aux_str2,"(A"//trim(adjustl(FMT))//")") "|"
    !   ades_data_rec(lin) = ades_data_rec(lin)(1:len_rec)//adjustl(aux_str2(1:len_field(3)))
    !   len_rec = len_rec + len_field(3)
    !ENDIF

    ! mode
    IF (.NOT. ades_flg(i)%rad .AND. .NOT. ades_flg(i)%occ) THEN
       WRITE(FMT_1,"(i20)") len_field(2) - 1
       WRITE(aux_str2,"(A"//trim(adjustl(FMT_1))//")") ades_obs(i)%mode
       aux_str2 = "|"//aux_str2
       ades_data_rec(lin) = ades_data_rec(lin)(1:len_rec)//adjustl(aux_str2(1:len_field(2)))
       len_rec = len_rec + len_field(2)
    ENDIF
    ! stn
    IF (.NOT. ades_flg(i)%rad) THEN
       WRITE(FMT,"(i20)") len_field(3)
       WRITE(aux_str2,"(A"//trim(adjustl(FMT))//")") "|"//trim(adjustl(ades_obs(i)%stn))
       ades_data_rec(lin) = ades_data_rec(lin)(1:len_rec)//adjustl(aux_str2(1:len_field(3)))
       len_rec = len_rec + len_field(3)
    ENDIF
    ! trx and rcv
    IF (ades_flg(i)%rad) THEN
       WRITE(FMT,"(i20)") len_field(4)
       WRITE(aux_str2,"(A"//trim(adjustl(FMT))//")") "|"//trim(adjustl(ades_obs(i)%trx))
       ades_data_rec(lin) = ades_data_rec(lin)(1:len_rec)//adjustl(aux_str2(1:len_field(4)))
       len_rec = len_rec + len_field(4)
       WRITE(FMT,"(i20)") len_field(5)
       WRITE(aux_str2,"(A"//trim(adjustl(FMT))//")") "|"//trim(adjustl(ades_obs(i)%rcv))
       ades_data_rec(lin) = ades_data_rec(lin)(1:len_rec)//adjustl(aux_str2(1:len_field(5)))
       len_rec = len_rec + len_field(5)
    ENDIF
    ! Location Group
    ! must not be present if stn is associated with a stationary MPC observatory code
    IF (ades_flg(i)%sys) THEN
       ! sys
       WRITE(FMT,"(i20)") len_field(6)
       WRITE(aux_str2,"(A"//trim(adjustl(FMT))//")") "|"//trim(adjustl(ades_obs(i)%sys))
       ades_data_rec(lin) = ades_data_rec(lin)(1:len_rec)//adjustl(aux_str2(1:len_field(6)))
       len_rec = len_rec + len_field(6)
       ! ctr
       WRITE(data,"(i20)") ades_obs(i)%ctr
       WRITE(FMT,"(i20)") len_field(7)
       WRITE(aux_str2,"(A"//trim(adjustl(FMT))//")") "|"//trim(adjustl(data))
       ades_data_rec(lin) = ades_data_rec(lin)(1:len_rec)//adjustl(aux_str2(1:len_field(7)))
       len_rec = len_rec + len_field(7)
       ! pos1
       WRITE(data,sfmt(8)) ades_obs(i)%pos(1)
       WRITE(FMT,"(i20)") len_field(8)
       WRITE(aux_str2,"(A"//trim(adjustl(FMT))//")") "|"//trim(adjustl(data))
       ades_data_rec(lin) = ades_data_rec(lin)(1:len_rec)//adjustl(aux_str2(1:len_field(8)))
       len_rec = len_rec + len_field(8)
       ! pos2
       WRITE(data,sfmt(8)) ades_obs(i)%pos(2)
       WRITE(FMT,"(i20)") len_field(8)
       WRITE(aux_str2,"(A"//trim(adjustl(FMT))//")") "|"//trim(adjustl(data))
       ades_data_rec(lin) = ades_data_rec(lin)(1:len_rec)//adjustl(aux_str2(1:len_field(8)))
       len_rec = len_rec + len_field(8)
       ! pos3
       WRITE(data,sfmt(8)) ades_obs(i)%pos(3)
       WRITE(FMT,"(i20)") len_field(8)
       WRITE(aux_str2,"(A"//trim(adjustl(FMT))//")") "|"//trim(adjustl(data))
       ades_data_rec(lin) = ades_data_rec(lin)(1:len_rec)//adjustl(aux_str2(1:len_field(8)))
       len_rec = len_rec + len_field(8)
       ! posCov11
       !WRITE(data,sfmt(9)) ades_obs(i)%posCov(1)
       !IF (ALL(ades_obs(i)%posCov .EQ. 0.d0)) data = ""
       !WRITE(FMT,"(i20)") len_field(9)
       !WRITE(aux_str2,"(A"//trim(adjustl(FMT))//")") "|"//trim(adjustl(data))
       !ades_data_rec(lin) = ades_data_rec(lin)(1:len_rec)//adjustl(aux_str2(1:len_field(9)))
       !len_rec = len_rec + len_field(9)
       ! posCov12
       !WRITE(data,sfmt(9)) ades_obs(i)%posCov(2)
       !IF (ALL(ades_obs(i)%posCov .EQ. 0.d0)) data = ""
       !WRITE(FMT,"(i20)") len_field(9)
       !WRITE(aux_str2,"(A"//trim(adjustl(FMT))//")") "|"//trim(adjustl(data))
       !ades_data_rec(lin) = ades_data_rec(lin)(1:len_rec)//adjustl(aux_str2(1:len_field(9)))
       !len_rec = len_rec + len_field(9)
       ! posCov13
       !WRITE(data,sfmt(9)) ades_obs(i)%posCov(3)
       !IF (ALL(ades_obs(i)%posCov .EQ. 0.d0)) data = ""
       !WRITE(FMT,"(i20)") len_field(9)
       !WRITE(aux_str2,"(A"//trim(adjustl(FMT))//")") "|"//trim(adjustl(data))
       !ades_data_rec(lin) = ades_data_rec(lin)(1:len_rec)//adjustl(aux_str2(1:len_field(9)))
       !len_rec = len_rec + len_field(9)
       ! posCov22
       !WRITE(data,sfmt(9)) ades_obs(i)%posCov(4)
       !IF (ALL(ades_obs(i)%posCov .EQ. 0.d0)) data = ""
       !WRITE(FMT,"(i20)") len_field(9)
       !WRITE(aux_str2,"(A"//trim(adjustl(FMT))//")") "|"//trim(adjustl(data))
       !ades_data_rec(lin) = ades_data_rec(lin)(1:len_rec)//adjustl(aux_str2(1:len_field(9)))
       !len_rec = len_rec + len_field(9)
       ! posCov23
       !WRITE(data,sfmt(9)) ades_obs(i)%posCov(5)
       !IF (ALL(ades_obs(i)%posCov .EQ. 0.d0)) data = ""
       !WRITE(FMT,"(i20)") len_field(9)
       !WRITE(aux_str2,"(A"//trim(adjustl(FMT))//")") "|"//trim(adjustl(data))
       !ades_data_rec(lin) = ades_data_rec(lin)(1:len_rec)//adjustl(aux_str2(1:len_field(9)))
       !len_rec = len_rec + len_field(9)
       ! posCov33
       !WRITE(data,sfmt(9)) ades_obs(i)%posCov(6)
       !IF (ALL(ades_obs(i)%posCov .EQ. 0.d0)) data = ""
       !WRITE(FMT,"(i20)") len_field(9)
       !WRITE(aux_str2,"(A"//trim(adjustl(FMT))//")") "|"//trim(adjustl(data))
       !ades_data_rec(lin) = ades_data_rec(lin)(1:len_rec)//adjustl(aux_str2(1:len_field(9)))
       !len_rec = len_rec + len_field(9)
    ENDIF

    !prog
    WRITE(FMT,"(i20)") len_field(28)
    WRITE(aux_str2,"(A"//trim(adjustl(FMT))//")") "|"//trim(adjustl(ades_obs(i)%prog))
    ades_data_rec(lin) = ades_data_rec(lin)(1:len_rec)//adjustl(aux_str2(1:len_field(28)))
    len_rec = len_rec + len_field(28)

    ! obsTime (always required)
    WRITE(FMT,"(i20)") len_field(11)
    WRITE(aux_str2,"(A"//trim(adjustl(FMT))//")") "|"//trim(adjustl(ades_obs(i)%obsTime))
    ades_data_rec(lin) = ades_data_rec(lin)(1:len_rec)//adjustl(aux_str2(1:len_field(11)))
    len_rec = len_rec + len_field(11)

    ! Observation Group
    IF (ades_flg(i)%opt) THEN
       ! ra
       WRITE(data,sfmt(12)) ades_obs(i)%ra
       WRITE(FMT,"(i20)") len_field(12)
       WRITE(aux_str2,"(A"//trim(adjustl(FMT))//")") "|"//trim(adjustl(data))
       ades_data_rec(lin) = ades_data_rec(lin)(1:len_rec)//adjustl(aux_str2(1:len_field(12)))
       len_rec = len_rec + len_field(12)
       ! dec
       WRITE(data,sfmt(13)) ades_obs(i)%dec
       WRITE(FMT,"(i20)") len_field(13)
       WRITE(aux_str2,"(A"//trim(adjustl(FMT))//")") "|"//trim(adjustl(data))
       ades_data_rec(lin) = ades_data_rec(lin)(1:len_rec)//adjustl(aux_str2(1:len_field(13)))
       len_rec = len_rec + len_field(13)
    ELSEIF (ades_flg(i)%occ) THEN
       ! raStar
       WRITE(data,sfmt(14)) ades_obs(i)%raStar
       WRITE(FMT,"(i20)") len_field(14)
       WRITE(aux_str2,"(A"//trim(adjustl(FMT))//")") "|"//trim(adjustl(data))
       ades_data_rec(lin) = ades_data_rec(lin)(1:len_rec)//adjustl(aux_str2(1:len_field(14)))
       len_rec = len_rec + len_field(14)
       ! decStar
       WRITE(data,sfmt(15)) ades_obs(i)%decStar
       WRITE(FMT,"(i20)") len_field(15)
       WRITE(aux_str2,"(A"//trim(adjustl(FMT))//")") "|"//trim(adjustl(data))
       ades_data_rec(lin) = ades_data_rec(lin)(1:len_rec)//adjustl(aux_str2(1:len_field(15)))
       len_rec = len_rec + len_field(15)
       ! deltaRA
       WRITE(data,sfmt(16)) ades_obs(i)%deltaRA
       WRITE(FMT,"(i20)") len_field(16)
       WRITE(aux_str2,"(A"//trim(adjustl(FMT))//")") "|"//trim(adjustl(data))
       ades_data_rec(lin) = ades_data_rec(lin)(1:len_rec)//adjustl(aux_str2(1:len_field(16)))
       len_rec = len_rec + len_field(16)
       ! deltaDec
       WRITE(data,sfmt(17)) ades_obs(i)%deltaDec
       WRITE(FMT,"(i20)") len_field(17)
       WRITE(aux_str2,"(A"//trim(adjustl(FMT))//")") "|"//trim(adjustl(data))
       ades_data_rec(lin) = ades_data_rec(lin)(1:len_rec)//adjustl(aux_str2(1:len_field(17)))
       len_rec = len_rec + len_field(17)
    ELSEIF (ades_flg(i)%delay) THEN
       ! delay
       WRITE(data,sfmt(18)) ades_obs(i)%delay
       WRITE(FMT,"(i20)") len_field(18)
       WRITE(aux_str2,"(A"//trim(adjustl(FMT))//")") "|"//trim(adjustl(data))
       ades_data_rec(lin) = ades_data_rec(lin)(1:len_rec)//adjustl(aux_str2(1:len_field(18)))
       len_rec = len_rec + len_field(18)
       ! rmsDelay
       WRITE(data,sfmt(19)) ades_obs(i)%rmsDelay
       WRITE(FMT,"(i20)") len_field(19)
       WRITE(aux_str2,"(A"//trim(adjustl(FMT))//")") "|"//trim(adjustl(data))
       ades_data_rec(lin) = ades_data_rec(lin)(1:len_rec)//adjustl(aux_str2(1:len_field(19)))
       len_rec = len_rec + len_field(19)
    ELSEIF (ades_flg(i)%doppler) THEN
       ! doppler
       WRITE(data,sfmt(20)) ades_obs(i)%doppler
       WRITE(FMT,"(i20)") len_field(20)
       WRITE(aux_str2,"(A"//trim(adjustl(FMT))//")") "|"//trim(adjustl(data))
       ades_data_rec(lin) = ades_data_rec(lin)(1:len_rec)//adjustl(aux_str2(1:len_field(20)))
       len_rec = len_rec + len_field(20)
       ! rmsDoppler
       WRITE(data,sfmt(21)) ades_obs(i)%rmsDoppler
       WRITE(FMT,"(i20)") len_field(21)
       WRITE(aux_str2,"(A"//trim(adjustl(FMT))//")") "|"//trim(adjustl(data))
       ades_data_rec(lin) = ades_data_rec(lin)(1:len_rec)//adjustl(aux_str2(1:len_field(21)))
       len_rec = len_rec + len_field(21)
    ENDIF

    ! astCat
    IF (ades_flg(i)%opt .OR. ades_flg(i)%occ) THEN
       WRITE(FMT_1,"(i20)") len_field(22) - 1
       WRITE(aux_str2,"(A"//trim(adjustl(FMT_1))//")") adjustr(ades_obs(i)%astCat)
       aux_str2 = "|"//aux_str2
       ades_data_rec(lin) = ades_data_rec(lin)(1:len_rec)//adjustl(aux_str2(1:len_field(22)))
       len_rec = len_rec + len_field(22)
    ENDIF

    ! Photometry Group
    IF (.NOT. ades_flg(i)%rad) THEN
       ! mag
       WRITE(FMT,"(i20)") len_field(23)
       IF (ades_obs(i)%mag .GT. 0d0 .AND. ades_obs(i)%mag .NE. 9.9d9) THEN
          WRITE(data,sfmt(23)) ades_obs(i)%mag
          WRITE(aux_str2,"(A"//trim(adjustl(FMT))//")") "|"//trim(adjustl(data))
       ELSE
          data = ""
          WRITE(aux_str2,"(A"//trim(adjustl(FMT))//")") "|"//data
       ENDIF
       ades_data_rec(lin) = ades_data_rec(lin)(1:len_rec)//adjustl(aux_str2(1:len_field(23)))
       len_rec = len_rec + len_field(23)
       ! band
       WRITE(FMT_1,"(i20)") len_field(24) - 1
       IF (ades_obs(i)%mag .GT. 0d0 .AND. ades_obs(i)%mag .NE. 9.9d9) THEN
          WRITE(aux_str2,"(A"//trim(adjustl(FMT_1))//")") adjustr(ades_obs(i)%band)
       ELSE
          aux_str2 = ""
       ENDIF
       aux_str2 = "|"//aux_str2
       ades_data_rec(lin) = ades_data_rec(lin)(1:len_rec)//adjustl(aux_str2(1:len_field(24)))
       len_rec = len_rec + len_field(24)
       ! photCat
       WRITE(FMT_1,"(i20)") len_field(25) - 1
       IF (ades_obs(i)%mag .GT. 0d0 .AND. ades_obs(i)%mag .NE. 9.9d9) THEN
          WRITE(aux_str2,"(A"//trim(adjustl(FMT_1))//")") adjustr(ades_obs(i)% photCat)
       ELSE
          aux_str2 = ""
       ENDIF
       aux_str2 = "|"//aux_str2
       ades_data_rec(lin) = ades_data_rec(lin)(1:len_rec)//adjustl(aux_str2(1:len_field(25)))
       len_rec = len_rec + len_field(25)
    ENDIF

    ! com and frq
    IF (ades_flg(i)%rad) THEN
       ! com
       WRITE(data,"(i20)")  ades_obs(i)%com
       WRITE(FMT,"(i20)") len_field(26)
       WRITE(aux_str2,"(A"//trim(adjustl(FMT))//")") "|"//trim(adjustl(data))
       ades_data_rec(lin) = ades_data_rec(lin)(1:len_rec)//adjustl(aux_str2(1:len_field(26)))
       len_rec = len_rec + len_field(26)
       ! frq
       WRITE(data,sfmt(27)) ades_obs(i)%frq
       WRITE(FMT,"(i20)") len_field(27)
       WRITE(aux_str2,"(A"//trim(adjustl(FMT))//")") "|"//trim(adjustl(data))
       ades_data_rec(lin) = ades_data_rec(lin)(1:len_rec)//adjustl(aux_str2(1:len_field(27)))
       len_rec = len_rec + len_field(27)
    ENDIF

    ! ref
    !WRITE(FMT,"(i20)") len_field(28)
    !WRITE(aux_str2,"(A"//trim(adjustl(FMT))//")") "|"//trim(adjustl(ades_obs(i)%ref))
    !ades_data_rec(lin) = ades_data_rec(lin)(1:len_rec)//adjustl(aux_str2(1:len_field(28)))
    !len_rec = len_rec + len_field(28)

    ! disc, Precision Group and notes
    IF (.NOT. ades_flg(i)%rad)  THEN
       ! disc
       WRITE(FMT,"(i20)") len_field(29)
       WRITE(aux_str2,"(A"//trim(adjustl(FMT))//")") "| "//trim(adjustl(ades_obs(i)%disc))
       ades_data_rec(lin) = ades_data_rec(lin)(1:len_rec)//adjustl(aux_str2(1:len_field(29)))
       len_rec = len_rec + len_field(29)

       ! subFmt
       WRITE(FMT,"(i20)") len_field(30)
       WRITE(aux_str2,"(A"//trim(adjustl(FMT))//")") "|"//trim(adjustl(ades_obs(i)%subFmt))
       ades_data_rec(lin) = ades_data_rec(lin)(1:len_rec)//adjustl(aux_str2(1:len_field(30)))
       len_rec = len_rec + len_field(30)

       ! Precision Group
       ! precTime
       WRITE(data,"(i20)")  ades_obs(i)%precTime
       IF (ades_obs(i)%precTime .EQ. 0) data = ""
       WRITE(FMT,"(i20)") len_field(31)
       WRITE(aux_str2,"(A"//trim(adjustl(FMT))//")") "|"//trim(adjustl(data))
       ades_data_rec(lin) = ades_data_rec(lin)(1:len_rec)//adjustl(aux_str2(1:len_field(31)))
       len_rec = len_rec + len_field(31)
       ! precRA
       CALL precision_group_control(ades_obs(i)%precRA,data)
       IF (ades_obs(i)%precRA .EQ. 0.d0) data = ""
       !WRITE(data,sfmt(32)) ades_obs(i)%precRA
       WRITE(FMT,"(i20)") len_field(32)
       WRITE(aux_str2,"(A"//trim(adjustl(FMT))//")") "|"//trim(adjustl(data))
       ades_data_rec(lin) = ades_data_rec(lin)(1:len_rec)//adjustl(aux_str2(1:len_field(32)))
       len_rec = len_rec + len_field(32)
       ! precDec
       CALL precision_group_control(ades_obs(i)%precDec,data)
       IF (ades_obs(i)%precDec .EQ. 0.d0) data = ""
       !WRITE(data,sfmt(33)) ades_obs(i)%precDec
       WRITE(FMT,"(i20)") len_field(33)
       WRITE(aux_str2,"(A"//trim(adjustl(FMT))//")") "|"//trim(adjustl(data))
       ades_data_rec(lin) = ades_data_rec(lin)(1:len_rec)//adjustl(aux_str2(1:len_field(33)))
       len_rec = len_rec + len_field(33)

       ! notes
       WRITE(FMT,"(i20)") len_field(34)
       WRITE(aux_str2,"(A"//trim(adjustl(FMT))//")") "|"//trim(adjustl(ades_obs(i)%notes))
       ades_data_rec(lin) = ades_data_rec(lin)(1:len_rec)//adjustl(aux_str2(1:len_field(34)))
       len_rec = len_rec + len_field(34)
    ENDIF
    ! remarks
    WRITE(FMT,"(i20)") len_field(35)
    WRITE(aux_str2,"(A"//trim(adjustl(FMT))//")") "|"//trim(adjustl(ades_obs(i)%remarks))
    ades_data_rec(lin) = ades_data_rec(lin)(1:len_rec)//adjustl(aux_str2(1:len_field(35)))
    len_rec = len_rec + len_field(35)

    ! deprecated
    IF (.NOT. ades_flg(i)%rad) THEN
       WRITE(FMT,"(i20)") len_field(36)
       WRITE(aux_str2,"(A"//trim(adjustl(FMT))//")") "|"//trim(adjustl(ades_obs(i)%deprecated))
       ades_data_rec(lin) = ades_data_rec(lin)(1:len_rec)//adjustl(aux_str2(1:len_field(36)))
       len_rec = len_rec + len_field(36)
    ENDIF
    ! Residuals Groups fields are handled by subroutine get_residuals_group

  END SUBROUTINE addobs_psv_rec

  ! ---------------------------------------------------------------- !
  ! FUNCTION isreal                                                  !
  ! It tells whether a character string contains only digits         !
  ! It uses the isnum function after having stripped the "." char    !
  ! (and "d"/"e"/+"/"-") from the string.                            !
  ! WARNING: assumes ASCII internal representation of characters     !
  ! ---------------------------------------------------------------- !

  LOGICAL FUNCTION isreal(string)

    CHARACTER(LEN=*),INTENT(IN) :: string   ! input string

    CHARACTER(LEN=len(string)) :: aux_str ! auxiliary string
    INTEGER :: idx,idx_dot ! index of "d"/"e"/"."/"+"/"-"

    LOGICAL, EXTERNAL  :: isnum

    aux_str = string

    ! remove "d" or "e" chars if after dot
    idx_dot = INDEX(aux_str,'.')
    idx = INDEX(aux_str,'d')
    IF ((idx_dot .GT. 0 .AND. idx .GT. idx_dot) .OR. (idx_dot .EQ. 0)) THEN
       IF (idx .EQ. 1) THEN
          aux_str = trim(aux_str(idx+1:len(aux_str)))
       ELSEIF (idx .GT. 1) THEN
          aux_str = trim(aux_str(1:idx-1))//trim(aux_str(idx+1:len(aux_str)))
       ENDIF
    ENDIF
    idx = INDEX(aux_str,'D')
    IF ((idx_dot .GT. 0 .AND. idx .GT. idx_dot) .OR. (idx_dot .EQ. 0)) THEN
       IF (idx .EQ. 1) THEN
          aux_str = trim(aux_str(idx+1:len(aux_str)))
       ELSEIF (idx .GT. 1) THEN
          aux_str = trim(aux_str(1:idx-1))//trim(aux_str(idx+1:len(aux_str)))
       ENDIF
    ENDIF
    idx = INDEX(aux_str,'e')
    IF ((idx_dot .GT. 0 .AND. idx .GT. idx_dot) .OR. (idx_dot .EQ. 0)) THEN
       IF (idx .EQ. 1) THEN
          aux_str = trim(aux_str(idx+1:len(aux_str)))
       ELSEIF (idx .GT. 1) THEN
          aux_str = trim(aux_str(1:idx-1))//trim(aux_str(idx+1:len(aux_str)))
       ENDIF
    ENDIF
    idx = INDEX(aux_str,'E')
    IF ((idx_dot .GT. 0 .AND. idx .GT. idx_dot) .OR. (idx_dot .EQ. 0)) THEN
       IF (idx .EQ. 1) THEN
          aux_str = trim(aux_str(idx+1:len(aux_str)))
       ELSEIF (idx .GT. 1) THEN
          aux_str = trim(aux_str(1:idx-1))//trim(aux_str(idx+1:len(aux_str)))
       ENDIF
    ENDIF

    ! remove "." 
    IF (idx_dot .EQ. 1) THEN
       aux_str = trim(aux_str(idx_dot+1:len(aux_str)))
    ELSEIF (idx_dot .GT. 1) THEN
       aux_str = trim(aux_str(1:idx_dot-1))//trim(aux_str(idx_dot+1:len(aux_str)))
    ENDIF

    ! remove "+"
    idx = INDEX(aux_str,'+')
    IF (idx .EQ. 1) THEN
       aux_str = trim(aux_str(idx+1:len(aux_str)))
    ELSEIF (idx .GT. 1) THEN
       aux_str = trim(aux_str(1:idx-1))//trim(aux_str(idx+1:len(aux_str)))
    ENDIF

    ! remove "-"
    idx = INDEX(aux_str,'-')
    IF (idx .EQ. 1) THEN
       aux_str = trim(aux_str(idx+1:len(aux_str)))
    ELSEIF (idx .GT. 1) THEN
       aux_str = trim(aux_str(1:idx-1))//trim(aux_str(idx+1:len(aux_str)))
    ENDIF

    isreal = isnum(trim(aux_str))

  END FUNCTION isreal

  ! -------------------------------------------------------- !
  ! SUBROUTINE sort_ades                                     !
  ! After a re-sorting of observations with respect to time, !
  ! sorts in the same way all ADES arrays associated with    !
  ! observations and .psv lines.                             !
  ! -------------------------------------------------------- !

  SUBROUTINE sort_ades(nobs,iperm)

    INTEGER, INTENT(IN)               :: nobs     ! number  of observations
    INTEGER, INTENT(IN), DIMENSION(:) :: iperm    ! array of indexes of re-sorted obs

    ! observation numbers: maximum, space left
    INCLUDE 'parobx.h90'

    TYPE(ades_observ),     DIMENSION(nobx) :: ades_obs_tmp
    TYPE(ades_flags),      DIMENSION(nobx) :: ades_flg_tmp
    TYPE(ades_precision),  DIMENSION(nobx) :: ades_prec_tmp
    TYPE(ades_fit_var),    DIMENSION(nobx) :: ades_fit_tmp
    LOGICAL,               DIMENSION(nobx) :: new_key_rec_tmp !ades_resm_def_tmp,ades_resc_def_tmp,
    INTEGER,               DIMENSION(nobx) :: nobs_to_ntot_idx_tmp
    CHARACTER(LEN=500),    DIMENSION(nobx) :: ades_key_rec_tmp,ades_data_rec_tmp

    INTEGER  :: iline,permline,i,j,obs_line

    ! first, sort arrays associated with observations data
    ! copy public arrays into temp vectors
    ades_obs_tmp(1:nobs) = ades_obs(1:nobs)
    ades_flg_tmp(1:nobs) = ades_flg(1:nobs)
    ades_prec_tmp(1:nobs) = ades_prec(1:nobs)
    ades_fit_tmp(1:nobs) = ades_fit(1:nobs)

    ! copy back public vectors in sorted order
    DO i = 1,nobs
       ades_obs(i) = ades_obs_tmp(iperm(i)) 
       ades_flg(i) = ades_flg_tmp(iperm(i))
       ades_prec(i) = ades_prec_tmp(iperm(i))
       ades_fit(i) = ades_fit_tmp(iperm(i)) 
    ENDDO

    ! sort arrays associated with lines of .psv file
    ! re-initialize arrays
    ades_data_rec_tmp(1:ntot_ades) = ades_data_rec(1:ntot_ades)
    ades_key_rec_tmp(1:ntot_ades) = ades_key_rec(1:ntot_ades)
    ades_data_rec = ""
    ades_key_rec = ""
    nobs_to_ntot_idx_tmp(1:nobs) = nobs_to_ntot_idx(1:nobs)
    ntot_to_nobs_idx = 0
    nobs_to_ntot_idx = 0
    new_key_rec_tmp(1:nobs) = new_key_rec(1:nobs)
    new_key_rec = .FALSE.
    new_key_rec(1) = .TRUE.
    ntot_ades = header_length + nobs + 1 ! consider 1 line for each observation + first keys line
    obs_line = ntot_ades - nobs + 1  ! line number corresponding to first line of observational data
    ntot_to_nobs_idx(obs_line) = 1
    nobs_to_ntot_idx(1) = obs_line
    ! add a fields line if type of observation has changed wrt previous observation
    ! define arrays ntot_to_nobs_idx and  nobs_to_ntot_idx
    DO j = 2,nobs
       IF ((ades_flg(j)%opt .NEQV. ades_flg(j-1)%opt) .OR. (ades_flg(j)%occ .NEQV. ades_flg(j-1)%occ) .OR. &
            & (ades_flg(j)%delay .NEQV. ades_flg(j-1)%delay) .OR. (ades_flg(j)%doppler .NEQV. ades_flg(j-1)%doppler) .OR. &
            & (ades_flg(j)%sys .NEQV. ades_flg(j-1)%sys)) THEN
          ntot_ades = ntot_ades + 1  ! add a line
          new_key_rec(j) = .TRUE.
          obs_line = obs_line + 2  ! obs_line + 1 is a fields line
          ntot_to_nobs_idx(obs_line) = j
          nobs_to_ntot_idx(j) = obs_line
       ELSE
          new_key_rec(j) = .FALSE.
          obs_line = obs_line + 1
          ntot_to_nobs_idx(obs_line) = j
          nobs_to_ntot_idx(j) = obs_line
       ENDIF
    ENDDO
    ! sort arrays of keys and data fields
    DO i = 1,nobs
       iline = nobs_to_ntot_idx(i)
       permline = nobs_to_ntot_idx_tmp(iperm(i))
       ades_data_rec(iline) = ades_data_rec_tmp(permline)
       IF (new_key_rec(i)) THEN
          DO j = 1,nobs
             IF (new_key_rec_tmp(j) .AND. (ades_flg_tmp(j)%opt .EQV. ades_flg(i)%opt) .AND. &
                  & (ades_flg_tmp(j)%occ .EQV. ades_flg(i)%occ) .AND. (ades_flg_tmp(j)%delay .EQV. ades_flg(i)%delay) .AND. &
                  & (ades_flg_tmp(j)%doppler .EQV. ades_flg(i)%doppler) .AND. (ades_flg_tmp(j)%sys .EQV. ades_flg(i)%sys)) THEN
                ades_key_rec(iline-1) = ades_key_rec_tmp(nobs_to_ntot_idx_tmp(j)-1)
                EXIT
             ENDIF
          ENDDO
       ENDIF
    ENDDO
  END SUBROUTINE sort_ades

  ! ------------------------------------------------------- !
  ! FUNCTION find_obs_ades                                  !
  ! Find an observation from a list, matching the given one !
  ! ------------------------------------------------------- !
  
  INTEGER FUNCTION find_obs_ades(obs,obs_t,obs_fl,obs_flt,m,double_obs)

    TYPE(ades_observ), INTENT(IN), DIMENSION(nobx) :: obs          ! array of ADES observations
    TYPE(ades_observ), INTENT(IN)                  :: obs_t        ! one new ADES observation
    TYPE(ades_flags),  INTENT(IN), DIMENSION(nobx) :: obs_fl       ! array of ADES flags
    TYPE(ades_flags),  INTENT(IN)                  :: obs_flt      ! ADES flags of new observation
    INTEGER,           INTENT(IN)                  :: m            ! number of obs
    INTEGER,           INTENT(OUT)                 :: double_obs   ! location of the duplicate observation

    INTEGER           :: j,len_str,len_str_t
    CHARACTER(LEN=30) :: time_str,time_str_t
    LOGICAL           :: same_type
    CHARACTER(LEN=16) :: ades_time,ades_time_t
    REAL(KIND=dkind)  :: seconds,seconds_t

    REAL(KIND=dkind), PARAMETER :: epst = 1.d-3 ! time (seconds) control

    find_obs_ades = 0
    double_obs = 0
    DO j = 1,m     ! loop on observations
       ! compare type of observations
       same_type = (obs_fl(j)%opt .EQV. obs_flt%opt) .AND. (obs_fl(j)%occ .EQV. obs_flt%occ) .AND. &
            & (obs_fl(j)%delay .EQV. obs_flt%delay) .AND. (obs_fl(j)%doppler .EQV. obs_flt%doppler) &
            & .AND. (obs_fl(j)%sys .EQV. obs_flt%sys)
       time_str = trim(adjustl(obs(j)%obsTime))
       time_str_t = trim(adjustl(obs_t%obsTime))
       CALL rmsp(time_str,len_str)
       CALL rmsp(time_str_t,len_str_t)
       ! compare strings yyyy-mm-ddThh:mm
       ades_time = time_str(1:16)
       ades_time_t = time_str_t(1:16)
       ! store seconds
       READ(time_str(18:len_str-1),*) seconds
       READ(time_str_t(18:len_str_t-1),*) seconds_t
       IF (ades_time .EQ. ades_time_t .AND. ABS(seconds - seconds_t) .LT. epst .AND. obs(j)%stn .EQ. obs_t%stn .AND. &
            & obs(j)%trx .EQ. obs_t%trx .AND. obs(j)%rcv .EQ. obs_t%rcv .AND. same_type) THEN
          ! observation found; is it a double?
           IF (find_obs_ades .NE. 0) THEN
             IF (double_obs .EQ. 0) THEN
                double_obs = j
                IF (ierrou .GT. 0) THEN
                   ! excluding deprecated case (X observation is often a duplicate)
                   IF(obs_t%deprecated.NE.'X' .AND. obs(find_obs_ades)%deprecated.NE.'X' .AND. obs(j)%deprecated.NE.'X')THEN
                      WRITE(ierrou,*)'find_obs_ades: two same time',find_obs_ades,j,obs_t%obsTime 
                      numerr=numerr+1
                   ENDIF
                ELSE
                   WRITE(*,*)'find_obs_ades: two same time', find_obs_ades,j,obs_t%obsTime
                ENDIF
             ELSE
                IF(ierrou .GT. 0)THEN ! disaster case: triple observation
                   WRITE(ierrou,*)'find_obs_ades: three same time', &
                        &                     find_obs_ades,double_obs,j,obs_t%obsTime
                ELSE
                   WRITE(*,*)'find_obs_ades: three same time',                &
                        &                     find_obs_ades,double_obs,j,obs_t%obsTime
                ENDIF
             ENDIF
          ELSE
             find_obs_ades = j ! normal case, no double
          ENDIF
       ENDIF
    ENDDO
  END FUNCTION find_obs_ades

  ! -------------------------------------- !
  ! FUNCTION changed_obs_ades              !
  ! Compare two ADES observations: TRUE if !
  ! astrometric data have changed          !
  ! -------------------------------------- !
  
  LOGICAL FUNCTION changed_obs_ades(obs,obs_t,obs_fl)

    TYPE(ades_observ), INTENT(IN) :: obs          ! old ADES observation
    TYPE(ades_observ), INTENT(IN) :: obs_t        ! new ADES observation
    TYPE(ades_flags),  INTENT(IN) :: obs_fl       ! ADES flags of old obs
    
    REAL(KIND=dkind),  PARAMETER  ::  epsd = 1.d-9, epsa = 1.d-3  ! control on angles (decimal degrees and arcsec)
    REAL(KIND=dkind),  PARAMETER  ::  epsdel = 1.d-8, epsdop = 1.d-2  ! control on radar (seconds and Hz)

    changed_obs_ades = .FALSE.
    IF (.NOT. obs_fl%rad) THEN
       IF (obs_fl%opt) THEN
          IF (ABS(obs%ra-obs_t%ra) .GT. epsd*ABS(obs%ra) .OR. ABS(obs%dec-obs_t%dec) .GT. epsd*ABS(obs%dec)) &
               & changed_obs_ades = .TRUE.
       ELSE ! ignored offset case (obsCenter value)
          IF (ABS(obs%raStar-obs_t%raStar) .GT. epsd*ABS(obs%raStar) .OR. &
               & ABS(obs%decStar-obs_t%decStar) .GT. epsd*ABS(obs%decStar) &
               & .OR. ABS(obs%deltaRA-obs_t%deltaRA) .GT. epsa*ABS(obs%deltaRA) &
               & .OR. ABS(obs%deltaDec-obs_t%deltaDec) .GT. epsa*ABS(obs%deltaDec) &
               & .OR. ABS(obs%dist-obs_t%dist) .GT. epsa*ABS(obs%dist) &
               & .OR. ABS(obs%pa-obs_t%pa) .GT. epsd*ABS(obs%pa)) changed_obs_ades = .TRUE.
       ENDIF
    ELSE
       IF (obs_fl%delay) THEN
          IF (ABS(obs%delay-obs_t%delay) .GT. epsdel*ABS(obs%delay)) changed_obs_ades = .TRUE.
       ELSE
          IF (ABS(obs%doppler-obs_t%doppler) .GT. epsdop*ABS(obs%doppler)) changed_obs_ades = .TRUE.
       ENDIF
    ENDIF
         
  END FUNCTION changed_obs_ades

  ! ---------------------------------------------------- !
  ! SUBROUTINE update_ades_mpc                           !
  ! Add information from .obs_psv to the one available   !
  ! from a file .rwo_psv. The number of observations can !
  ! increase, from m to a maximum m+mt.                  !
  ! ---------------------------------------------------- !

  SUBROUTINE update_ades_mpc(ades_obs_out,ades_flg_out,ades_prec_out,ades_fit_out,m,ntot_ades_out,ades_key_rec_out, &
       & ades_data_rec_out,ntot_to_nobs_out,nobs_to_ntot_out,new_key_rec_out,ades_obs_tmp,ades_flg_tmp, &
       & ades_prec_tmp,ades_fit_tmp,mt,ades_key_rec_tmp,ades_data_rec_tmp,nobs_to_ntot_tmp,change,mnew)

    TYPE(ades_observ),    INTENT(INOUT), DIMENSION(nobx) :: ades_obs_out
    TYPE(ades_flags),     INTENT(INOUT), DIMENSION(nobx) :: ades_flg_out
    TYPE(ades_precision), INTENT(INOUT), DIMENSION(nobx) :: ades_prec_out
    TYPE(ades_fit_var),   INTENT(INOUT), DIMENSION(nobx) :: ades_fit_out
    INTEGER,              INTENT(INOUT)                  :: m
    INTEGER,              INTENT(INOUT)                  :: ntot_ades_out
    CHARACTER(LEN=600),   INTENT(INOUT), DIMENSION(nobx) :: ades_key_rec_out
    CHARACTER(LEN=600),   INTENT(INOUT), DIMENSION(nobx) :: ades_data_rec_out
    INTEGER,              INTENT(INOUT), DIMENSION(nobx) :: ntot_to_nobs_out
    INTEGER,              INTENT(INOUT), DIMENSION(nobx) :: nobs_to_ntot_out
    LOGICAL,              INTENT(INOUT), DIMENSION(nobx) :: new_key_rec_out

    TYPE(ades_observ),    INTENT(IN), DIMENSION(nobx) :: ades_obs_tmp
    TYPE(ades_flags),     INTENT(IN), DIMENSION(nobx) :: ades_flg_tmp
    TYPE(ades_precision), INTENT(IN), DIMENSION(nobx) :: ades_prec_tmp
    TYPE(ades_fit_var),   INTENT(IN), DIMENSION(nobx) :: ades_fit_tmp
    INTEGER,              INTENT(IN)                  :: mt
    CHARACTER(LEN=600),   INTENT(IN), DIMENSION(nobx) :: ades_key_rec_tmp
    CHARACTER(LEN=600),   INTENT(IN), DIMENSION(nobx) :: ades_data_rec_tmp
    INTEGER,              INTENT(IN), DIMENSION(nobx) :: nobs_to_ntot_tmp

    LOGICAL,              INTENT(OUT)  :: change      ! true if data in .obs_psv have changed wrt .rwo_psv
    INTEGER,              INTENT(OUT)  :: mnew

    INTEGER  :: j,mj,double,k

    ! monitor changes
    mnew = m
    change = .FALSE.
    ! scan supposedly new observations
    DO j = 1,mt
       mj = find_obs_ades(ades_obs_out,ades_obs_tmp(j),ades_flg_out,ades_flg_tmp(j),m,double)
       IF (mj .NE. 0 .AND. double .EQ. 0) THEN
          ! the observation was already there
          IF (changed_obs_ades(ades_obs_out(mj),ades_obs_tmp(j),ades_flg_out(mj))) THEN 
             ! the observation has changed
             change = .TRUE.
             ades_obs_out(mj) = ades_obs_tmp(j)
             ades_flg_out(mj) = ades_flg_tmp(j)
             ades_prec_out(mj) = ades_prec_tmp(j)
             ades_fit_out(mj) = ades_fit_tmp(j)
             ! substitution of line of formatted data to be written on output .rwo_psv file
             ades_data_rec_out(nobs_to_ntot_out(mj)) = ades_data_rec_tmp(nobs_to_ntot_tmp(j))
          ELSEIF (ades_obs_out(mj)%astCat .NE. ades_obs_tmp(j)%astCat) THEN
             change = .TRUE.
             ades_obs_out(mj)%astCat = ades_obs_tmp(j)%astCat
             ades_obs_out(mj)%biasRA = ades_obs_tmp(j)%biasRA
             ades_obs_out(mj)%biasDec = ades_obs_tmp(j)%biasDec
          ELSEIF (ades_obs_out(mj)%photCat .NE. ades_obs_tmp(j)%photCat) THEN
             change = .TRUE.
             ades_obs_out(mj)%photCat = ades_obs_tmp(j)%photCat
             ades_obs_out(mj)%biasMag = ades_obs_tmp(j)%biasMag
          ENDIF
       ELSEIF (mj .NE. 0 .AND. double .NE. 0) THEN
          ! the observation was already there, in double copy
          IF (changed_obs_ades(ades_obs_out(mj),ades_obs_tmp(j),ades_flg_out(mj))) THEN
             IF (changed_obs_ades(ades_obs_out(double),ades_obs_tmp(j),ades_flg_out(double))) THEN
                change = .TRUE.
                ! double and changed
                WRITE(*,*)'update_ades_mpc: double and changed'
                WRITE(*,*)' records ',mj,' and ',double,' in .rwo_psv'
                WRITE(*,*)' record ',j,' in .obs_psv'
             ENDIF
          ENDIF
       ELSEIF (mj .EQ. 0) THEN
          ! the observation is new: add it
          change = .TRUE.
          mnew = mnew+1
          ades_obs_out(mnew) = ades_obs_tmp(j)
          ades_flg_out(mnew) = ades_flg_tmp(j)
          ades_prec_out(mnew) = ades_prec_tmp(j)
          ades_fit_out(mnew) = ades_fit_tmp(j)
          IF ((ades_flg_out(mnew)%opt .NEQV. ades_flg_out(mnew-1)%opt) .OR. &
               & (ades_flg_out(mnew)%occ .NEQV. ades_flg_out(mnew-1)%occ) .OR. &
               & (ades_flg_out(mnew)%delay .NEQV. ades_flg_out(mnew-1)%delay) .OR. &
               & (ades_flg_out(mnew)%doppler .NEQV. ades_flg_out(mnew-1)%doppler) .OR. &
               & (ades_flg_out(mnew)%sys .NEQV. ades_flg_out(mnew-1)%sys)) THEN
             new_key_rec_out(mnew) = .TRUE.
             ntot_ades_out = ntot_ades_out + 2  ! add one data line and one key line
          ELSE
             new_key_rec_out(mnew) = .FALSE.
             ntot_ades_out = ntot_ades_out + 1   ! add only one data line
          ENDIF
          ntot_to_nobs_out(ntot_ades_out) = mnew
          nobs_to_ntot_out(mnew) = ntot_ades_out
          ! addition of the lines of formatted data to be written on output .rwo_psv file
          IF (new_key_rec_out(mnew)) THEN
             DO k = j,1,-1
                IF (ades_key_rec_tmp(nobs_to_ntot_tmp(k)-1) .NE. "") THEN
                   ades_key_rec_out(nobs_to_ntot_out(mnew)-1) = ades_key_rec_tmp(nobs_to_ntot_tmp(k)-1)
                   EXIT
                ENDIF
             ENDDO
          ENDIF
          ades_data_rec_out(nobs_to_ntot_out(mnew)) = ades_data_rec_tmp(nobs_to_ntot_tmp(j))
       ENDIF
    ENDDO
  END SUBROUTINE update_ades_mpc

  ! ---------------------------------------------------- !
  ! SUBROUTINE update_ades_rwo                           !
  ! Add information from .rwo_psv (selection flags and   !
  ! manually fixed weights with their flags force_w) to  !
  ! the one available from ADES data. The number of      !
  ! observations cannot increase, remains m. Arrays      !
  ! ades_data_rec and ades_key_rec do not need updates,  !
  ! since only Residuals Group changes, which is handled !
  ! during writing process of .rwo_psv file.             !
  ! ---------------------------------------------------- !
  
  SUBROUTINE update_ades_rwo(ades_obs_out,ades_flg_out,ades_prec_out,ades_fit_out,ades_obs_tmp, &
       & ades_flg_tmp,ades_prec_tmp,ades_fit_tmp,m,mt,change_rwo)

    TYPE(ades_observ),    INTENT(INOUT), DIMENSION(nobx) :: ades_obs_out
    TYPE(ades_flags),     INTENT(INOUT), DIMENSION(nobx) :: ades_flg_out
    TYPE(ades_precision), INTENT(INOUT), DIMENSION(nobx) :: ades_prec_out
    TYPE(ades_fit_var),   INTENT(INOUT), DIMENSION(nobx) :: ades_fit_out
    TYPE(ades_observ),    INTENT(IN),    DIMENSION(nobx) :: ades_obs_tmp
    TYPE(ades_flags),     INTENT(IN),    DIMENSION(nobx) :: ades_flg_tmp
    TYPE(ades_precision), INTENT(IN),    DIMENSION(nobx) :: ades_prec_tmp
    TYPE(ades_fit_var),   INTENT(IN),    DIMENSION(nobx) :: ades_fit_tmp
    INTEGER,              INTENT(IN)                     :: m,mt
    LOGICAL,              INTENT(OUT)                    :: change_rwo ! true if data have changed wrt rwo

    INTEGER  :: j,mj,double

    change_rwo = .FALSE.
    DO j=1,mt
       mj=find_obs_ades(ades_obs_out,ades_obs_tmp(j),ades_flg_out,ades_flg_tmp(j),m,double)
       IF (mj .NE. 0 .AND. double .EQ. 0) THEN
          ! this observation was already there
          IF (changed_obs_ades(ades_obs_out(mj),ades_obs_tmp(j),ades_flg_out(mj))) THEN  
             ! the observation has changed
             change_rwo = .TRUE.
          ELSE
             ! if no change to observation then preserve the
             ! manually fixed weights and selection flags
             IF(ades_fit_tmp(j)%force_w_flags(1))THEN
                ades_fit_out(mj)%force_w_flags(1) = ades_fit_tmp(j)%force_w_flags(1)
                ades_obs_out(mj)%sigRA = ades_obs_tmp(j)%sigRA
                ades_flg_out(mj)%sigRA = ades_flg_tmp(j)%sigRA 
                ades_prec_out(mj)%sigRA = ades_prec_tmp(j)%sigRA
                ades_obs_out(mj)%sigDelay = ades_obs_tmp(j)%sigDelay
                ades_flg_out(mj)%sigDelay = ades_flg_tmp(j)%sigDelay
                ades_prec_out(mj)%sigDelay = ades_prec_tmp(j)%sigDelay
             ENDIF
             IF(ades_fit_tmp(j)%force_w_flags(2))THEN
                ades_fit_out(mj)%force_w_flags(2) = ades_fit_tmp(j)%force_w_flags(2)
                ades_obs_out(mj)%sigDec = ades_obs_tmp(j)%sigDec
                ades_flg_out(mj)%sigDec = ades_flg_tmp(j)%sigDec
                ades_prec_out(mj)%sigDec = ades_prec_tmp(j)%sigDec
                ades_obs_out(mj)%sigDoppler = ades_obs_tmp(j)%sigDoppler
                ades_flg_out(mj)%sigDoppler = ades_flg_tmp(j)%sigDoppler
                ades_prec_out(mj)%sigDoppler = ades_prec_tmp(j)%sigDoppler
             ENDIF
             IF(ades_fit_tmp(j)%force_w_flags(1) .AND. ades_fit_tmp(j)%force_w_flags(2))THEN
                ades_obs_out(mj)%sigCorr = ades_obs_tmp(j)%sigCorr
                ades_prec_out(mj)%sigCorr = ades_prec_tmp(j)%sigCorr
                ades_flg_out(mj)%sigCorr = ades_flg_tmp(j)%sigCorr
             ENDIF 
             IF(ades_flg_tmp(j)%mag .AND. ades_obs_tmp(j)%mag .EQ. 9.9d9) THEN
                ades_obs_out(mj)%sigMag = ades_obs_tmp(j)%sigMag
                ades_flg_out(mj)%sigMag = ades_flg_tmp(j)%sigMag
                ades_prec_out(mj)%sigMag = ades_prec_tmp(j)%sigMag
             ENDIF
             ! selection flags are preserved anyway
             ades_obs_out(mj)%selAst = ades_obs_tmp(j)%selAst
             ades_obs_out(mj)%selPhot = ades_obs_tmp(j)%selPhot
             ades_obs_out(mj)%selDelay = ades_obs_tmp(j)%selDelay
             ades_obs_out(mj)%selDoppler = ades_obs_tmp(j)%selDoppler
             ades_fit_out(mj)%ast_flag = ades_fit_tmp(j)%ast_flag
             ades_fit_out(mj)%mag_flag = ades_fit_tmp(j)%mag_flag
             ades_fit_out(mj)%radar_flag = ades_fit_tmp(j)%radar_flag
          ENDIF
       ELSEIF (mj .NE. 0 .AND. double .NE. 0) THEN
          ! the observation was already there, in double copy
          IF (changed_obs_ades(ades_obs_out(mj),ades_obs_tmp(j),ades_flg_out(mj))) THEN
             IF (changed_obs_ades(ades_obs_out(double),ades_obs_tmp(j),ades_flg_out(double))) THEN
                change_rwo = .TRUE.
                ! double and changed
                WRITE(*,*)'update_ades_rwo: double and changed'
                WRITE(ierrou,*)'update_ades_rwo: double and changed'
                numerr=numerr+1
                WRITE(*,*)' records ',mj,' and ',double,' in .obs_psv'
                WRITE(*,*)' record ',j,' in .rwo_psv'
                WRITE(ierrou,*)' records ',mj,' and ',double,' in .obs_psv'
                WRITE(ierrou,*)' record ',j,' in .rwo_psv'
             ELSE
                ! OK, it is the double
                ! it is the same, so preserve the selection flags
                ! manually fixed weights and selection flags
                IF(ades_fit_tmp(j)%force_w_flags(1))THEN
                   ades_fit_out(double)%force_w_flags(1) = ades_fit_tmp(j)%force_w_flags(1)
                   ades_obs_out(double)%sigRA = ades_obs_tmp(j)%sigRA
                   ades_flg_out(double)%sigRA = ades_flg_tmp(j)%sigRA
                   ades_prec_out(double)%sigRA = ades_prec_tmp(j)%sigRA
                   ades_obs_out(double)%sigDelay = ades_obs_tmp(j)%sigDelay
                   ades_flg_out(double)%sigDelay = ades_flg_tmp(j)%sigDelay
                   ades_prec_out(double)%sigDelay = ades_prec_tmp(j)%sigDelay
                ENDIF
                IF(ades_fit_tmp(j)%force_w_flags(2))THEN
                   ades_fit_out(double)%force_w_flags(2) = ades_fit_tmp(j)%force_w_flags(2)
                   ades_obs_out(double)%sigDec = ades_obs_tmp(j)%sigDec
                    ades_flg_out(double)%sigDec = ades_flg_tmp(j)%sigDec
                   ades_prec_out(double)%sigDec = ades_prec_tmp(j)%sigDec
                   ades_obs_out(double)%sigDoppler = ades_obs_tmp(j)%sigDoppler
                   ades_flg_out(double)%sigDoppler = ades_flg_tmp(j)%sigDoppler
                   ades_prec_out(double)%sigDoppler = ades_prec_tmp(j)%sigDoppler
                ENDIF
                IF(ades_fit_tmp(j)%force_w_flags(1) .AND. ades_fit_tmp(j)%force_w_flags(2))THEN
                   ades_obs_out(double)%sigCorr = ades_obs_tmp(j)%sigCorr
                   ades_prec_out(double)%sigCorr = ades_prec_tmp(j)%sigCorr
                   ades_flg_out(double)%sigCorr = ades_flg_tmp(j)%sigCorr
                ENDIF
                IF(ades_flg_tmp(j)%mag .AND. ades_obs_tmp(j)%mag .EQ. 9.9d9) THEN
                   ades_obs_out(double)%sigMag = ades_obs_tmp(j)%sigMag
                   ades_flg_out(double)%sigMag = ades_flg_tmp(j)%sigMag
                   ades_prec_out(double)%sigMag = ades_prec_tmp(j)%sigMag
                ENDIF
                ! selection flags are preserved anyway
                ades_obs_out(double)%selAst = ades_obs_tmp(j)%selAst
                ades_obs_out(double)%selPhot = ades_obs_tmp(j)%selPhot
                ades_obs_out(double)%selDelay = ades_obs_tmp(j)%selDelay
                ades_obs_out(double)%selDoppler = ades_obs_tmp(j)%selDoppler
                ades_fit_out(double)%ast_flag = ades_fit_tmp(j)%ast_flag
                ades_fit_out(double)%mag_flag = ades_fit_tmp(j)%mag_flag
                ades_fit_out(double)%radar_flag = ades_fit_tmp(j)%radar_flag
             ENDIF
          ELSE
             ! it is the same, so preserve the selection flags
             ! manually fixed weights and selection flags
              IF(ades_fit_tmp(j)%force_w_flags(1))THEN
                ades_fit_out(mj)%force_w_flags(1) = ades_fit_tmp(j)%force_w_flags(1)
                ades_obs_out(mj)%sigRA = ades_obs_tmp(j)%sigRA
                ades_flg_out(mj)%sigRA = ades_flg_tmp(j)%sigRA 
                ades_prec_out(mj)%sigRA = ades_prec_tmp(j)%sigRA
                ades_obs_out(mj)%sigDelay = ades_obs_tmp(j)%sigDelay
                ades_flg_out(mj)%sigDelay = ades_flg_tmp(j)%sigDelay
                ades_prec_out(mj)%sigDelay = ades_prec_tmp(j)%sigDelay
             ENDIF
             IF(ades_fit_tmp(j)%force_w_flags(2))THEN
                ades_fit_out(mj)%force_w_flags(2) = ades_fit_tmp(j)%force_w_flags(2)
                ades_obs_out(mj)%sigDec = ades_obs_tmp(j)%sigDec
                ades_flg_out(mj)%sigDec = ades_flg_tmp(j)%sigDec
                ades_prec_out(mj)%sigDec = ades_prec_tmp(j)%sigDec
                ades_obs_out(mj)%sigDoppler = ades_obs_tmp(j)%sigDoppler
                ades_flg_out(mj)%sigDoppler = ades_flg_tmp(j)%sigDoppler
                ades_prec_out(mj)%sigDoppler = ades_prec_tmp(j)%sigDoppler
             ENDIF
             IF(ades_fit_tmp(j)%force_w_flags(1) .AND. ades_fit_tmp(j)%force_w_flags(2))THEN
                ades_obs_out(mj)%sigCorr = ades_obs_tmp(j)%sigCorr
                ades_prec_out(mj)%sigCorr = ades_prec_tmp(j)%sigCorr
                ades_flg_out(mj)%sigCorr = ades_flg_tmp(j)%sigCorr
             ENDIF                
             IF(ades_flg_tmp(j)%mag .AND. ades_obs_tmp(j)%mag .EQ. 9.9d9) THEN
                ades_obs_out(mj)%sigMag = ades_obs_tmp(j)%sigMag
                ades_flg_out(mj)%sigMag = ades_flg_tmp(j)%sigMag
                ades_prec_out(mj)%sigMag = ades_prec_tmp(j)%sigMag
             ENDIF
             ! selection flags are preserved anyway
             ades_obs_out(mj)%selAst = ades_obs_tmp(j)%selAst
             ades_obs_out(mj)%selPhot = ades_obs_tmp(j)%selPhot
             ades_obs_out(mj)%selDelay = ades_obs_tmp(j)%selDelay
             ades_obs_out(mj)%selDoppler = ades_obs_tmp(j)%selDoppler
             ades_fit_out(mj)%ast_flag = ades_fit_tmp(j)%ast_flag
             ades_fit_out(mj)%mag_flag = ades_fit_tmp(j)%mag_flag
             ades_fit_out(mj)%radar_flag = ades_fit_tmp(j)%radar_flag
          ENDIF
       ELSEIF (mj .EQ. 0) THEN
          ! if it is not found in .rwo_psv, leave the input weights and selection flags
          change_rwo = .TRUE.
       ENDIF
    ENDDO
    ! check if there are extra (added) observations in .obs_psv file
    IF (mt .LT. m) change_rwo = .TRUE.
    ! erase all residuals and chi (for safety: .obs_psv should not contain Residuals group)
    ades_obs_out%resRA = 0.d0
    ades_obs_out%resDec = 0.d0
    ades_obs_out%resMag = 0.d0
    ades_obs_out%resDelay = 0.d0
    ades_obs_out%resDoppler = 0.d0
    ades_fit_out%ades_resc_def = .FALSE.
    ades_fit_out%ades_resm_def = .FALSE.
    ades_fit_out%chi = 0.d0
    
  END SUBROUTINE update_ades_rwo

END MODULE ades

! MODULE astrometric_observations
! Handling of astrometric observations

! TYPE ast_obs
! 
!     astrometric observations (optical or radar)
! TYPE ast_wbsr
!     weighting, bias, selection flag and residuals of astrometric observation
! CONTAINS

! SUBROUTINE write_rwo(file,obs,obsw,n,error_model,rms,rmsmag,sort,unit0)
! output observations in .rwo format
! SUBROUTINE write_rwo_rad(unit, obs, obsw)
! SUBROUTINE write_rwo_opt(unit,obs,obsw,error,nrec)
! SUBROUTINE write_rwo_pos(unit,obs,obsw)
! SUBROUTINE write_rwo_header(unit, error_model, nrec, rms, rmsmag)

! SUBROUTINE read_rwo(file,obs,obsw,n,error_model,rms,rmsag,unit0,eof0)
! input observations in .rwo format
! SUBROUTINE read_rwo_rec(unit,nr,obs1,obsw1,error, eof)
! SUBROUTINE read_rwo_header(unit, error_model, error, rms, rmsmag)

! SUBROUTINE mpc_obs_input(obs_found,obs,nobs,unit,filnam,eof,ons_mode)
! SUBROUTINE mpcrec_transform(mpcrec,obs,complete,error,skip_rad,ons_mode)
!     transformation of a MPC record into an observation (TYPE ast_obs)
! SUBROUTINE mpcrec_add_obspos(obs,error,mpcrec1,mpcrec2)

! SUBROUTINE jpl_radarobs(obs_found,file,obs,nobs)
! SUBROUTINE jplradar_transform(rec,obs,error)
! CHARACTER(LEN=11) FUNCTION radar_resid_string(value,defined)
! INTEGER FUNCTION radar_station(stastr)

! SUBROUTINE input_obs(obsdir,astna0,precob,error_model,obs0,obs,obsw,m,iun20,change,rms,rmsmag)
!     input of observations: reads sequentially the .rwo, .obs, .rad file and combines
!     the data according to the precedence rule specified by precob

! SUBROUTINE addobs_mpc(obs,obsw,m,obst,obswt,mt,mnew,change)
! SUBROUTINE addobs_rwo(obs,obsw,m,obst,obswt,mt,change)
! SUBROUTINE update_ades_mpc
! SUBROUTINE update_ades_rwo
! SUBROUTINE new_ades_obs
! INTEGER FUNCTION find_obs(obst,obs,m,double)
! LOGICAL FUNCTION changed_obs(obst,obswt,obs,obsw)
! INTEGER FUNCTION find_obs_ades(obst,obs,m,double)
! LOGICAL FUNCTION changed_obs_ades(obst,obswt,obs,obsw)

! SUBROUTINE aster_radius(objid,obstype,m)

! SUBROUTINE rdanga(string,angle,acc,error)
! SUBROUTINE observ_rms(obs,error_model,init,obsw,n)
! SUBROUTINE observ_rms_singobs(obs,error_model,init,obsw)
! SUBROUTINE astrow_bias(error_model,obs,obsw)
! SUBROUTINE astrow_nobias(error_model,obs,obsw)
! DOUBLE PRECISION FUNCTION magrms_new(magstr,tdt,idsta,typ)
! SUBROUTINE convert_ades_fitobs(error_model,observ,wei_bia)
! SUBROUTINE convert_fitobs_ades(observ,wei_bia,chi)

! HEADERS AND MODULES
! astrometric_observations.o: \
!	../include/parcmc.h90 \
!	../include/parobx.h90 \
!	../include/sysdep.h90 \
!	fund_const.o \
!	output_control.o \
!	reference_systems.o \
!	station_coordinates.o 


MODULE astrometric_observations
  USE fund_const
  USE output_control
  USE name_rules
  USE file_oper
  USE char_str
  USE reference_systems
  USE station_coordinates
  USE ades        
  IMPLICIT NONE
  PRIVATE
  ! SUBROUTINEs
  PUBLIC :: input_obs,mpc_obs_input, addobs_rwo, addobs_mpc, aster_radius,jpl_radarobs,find_obs
  PUBLIC :: read_rwo,write_rwo,read_rwo_rec,read_rwo_header, write_rwo_opt,write_rwo_header
  PUBLIC :: observ_rms,observ_rms_singobs
  PUBLIC :: convert_ades_fitobs,convert_fitobs_ades
  ! OPERATOR OVERLOADING
  PUBLIC :: drop_geopos, add_geopos
  ! CHANGE OF RWO FORMAT

  ! common data for sub-module observ_rms_bias
  ! to read cbm10.sta
  INTEGER,          PARAMETER      :: nstarmsx=100
  DOUBLE PRECISION, SAVE           :: rmsstabin(2,nstarmsx)     ! rms(a*cosd), rmsd, in arcsec
  CHARACTER*3,      SAVE           :: obscodrms(nstarmsx)       ! star catalog code (' '=not known)
  INTEGER,          SAVE           :: indstarms(nstarmsx)       ! lexicographic sort index
  INTEGER,          SAVE           :: nstarms                   ! actual number
  ! to read cbm10.rules, fcct14.rules, vfcc17.rules
  INTEGER,          PARAMETER      :: nrulesx=300
  DOUBLE PRECISION, SAVE           :: rmsbin(2,nrulesx)         ! rms(a*cosd), rmsd, in arcsec
  CHARACTER*4,      SAVE           :: obscatcod(nrulesx)        ! star catalog code (' '=not known)
  CHARACTER*6,      SAVE           :: errmod_rules(nrulesx)     ! star catalog code (' '=not known)
  INTEGER,          SAVE           :: indrules(nrulesx)         ! lexicographic sort index
  DOUBLE PRECISION, SAVE           :: mjd_date_rules(2,nrulesx) ! lower_date, upper_date (range date for rule)
  INTEGER,          SAVE           :: nrules                    ! actual number
 
  CHARACTER*(6), PUBLIC :: change_format ! INPUT,OUTPUT, ' ' 

  LOGICAL, PUBLIC :: ons_name  ! if true, the identifier in an .obs file is not a designation

  !INTEGER, PARAMETER :: len_objdes = 9   ! length of the object designation
  !PUBLIC len_objdes

  ! Generic astrometric observation
  TYPE ast_obs
     CHARACTER(LEN=name_len)         :: objdes     ! object designation
     CHARACTER(LEN=1)                :: type       ! observation type:
                                                   !    O = optical
                                                   !    R = radar range
                                                   !    V = radar range rate
                                                   !    S = satellite/space
     CHARACTER(LEN=1)                :: tech       ! observation technology (column 15 of MPC record):
                                                   !    P = Photographic
                                                   !    E = Encoder
                                                   !    C = CCD
                                                   !    T = Transit or meridian circle
                                                   !    M = Micrometer
                                                   !    c = corrected for center of mass (radar only)
                                                   !    s = surface bounce (radar only)
     CHARACTER(LEN=1)                :: par_unit   ! parallax unit for satellite observations:
                                                   !    1 = km
                                                   !    2 = au
     CHARACTER(LEN=1)                :: note       ! note on observation circumstances (column 14 of MPC record)
     CHARACTER(LEN=2)                :: catcodmpc  ! MPC catalog flag:
                                                   ! a = USNO-A1.0
                                                   ! b = USNO-SA1.0
                                                   ! c = USNO-A2.0
                                                   ! d = USNO-SA2.0
                                                   ! e = UCAC-1
                                                   ! f = Tycho-1
                                                   ! g = Tycho-2
                                                   ! h = GSC-1.0
                                                   ! i = GSC-1.1
                                                   ! j = GSC-1.2
                                                   ! k = GSC-2.2
                                                   ! l = ACT
                                                   ! m = GSC-ACT
                                                   ! n = TRC
                                                   ! o = USNO-B1.0
                                                   ! p = PPM
                                                   ! q = UCAC-2-beta
                                                   ! r = UCAC-2
                                                   ! s = USNO-B2.0
                                                   ! t = UCAC-3-beta
                                                   ! u = UCAC-3
                                                   ! v = NOMAD
                                                   ! w = CMC
                                                   ! x = Hip-2
                                                   ! z = GSC (unspecified version)
                                                   ! ' ' = unknown

     CHARACTER(LEN=7)                :: obscod_s   ! observer's MPC code (string, official); for radar, TX(A3)//' '//RX(A3)
     INTEGER                         :: obscod_i   ! observer's MPC code (integer, internal representation); for radar,
                                                   ! encoding of both transmitter and receiver (see jplradar_transform)
     DOUBLE PRECISION                :: time_utc   ! observation UTC time (MJD)
     DOUBLE PRECISION                :: time_tdt   ! observation TDT time (MJD)
     DOUBLE PRECISION, DIMENSION(2)  :: coord      ! observation coordinates:
                                                   !    RA,DEC (optical & satellite) (rad)
                                                   !    range,range rate (radar)
     DOUBLE PRECISION                :: mag        ! object apparent magnitude
     CHARACTER(LEN=6)                :: mag_str    ! magnitude original string in the record
     CHARACTER(LEN=1)                :: mag_band   ! magnitude color band
     LOGICAL                         :: mag_def    ! TRUE if magnitude is defined
     DOUBLE PRECISION                :: acc_time   ! accuracy of time (as given in input format) (d)
     DOUBLE PRECISION, DIMENSION(2)  :: acc_coord  ! accuracy of coordinates (as given in input format)
     DOUBLE PRECISION                :: acc_mag    ! accuracy of magnitude (as given in input format)
     ! Observer's position and velocity are given with respect to the mean ecliptic and equinox J2000
     ! WARNING: in case of radar observations:
     !         obspos = body-fixed position of transmitting station
     !         obsvel = body-fixed position of receiving station
     DOUBLE PRECISION, DIMENSION(3)  :: obspos     ! observer's geocentric position (au)
     DOUBLE PRECISION, DIMENSION(3)  :: obsvel     ! observer's geocentric velocity (au/d)
     ! The following array is used only for roving observatory observations (obs%type='O' and obs%obscod_s='247')
     ! and is ALLOCATEd and initialized by SUBROUTINE mpcrec_add_obspos or SUBROUTINE read_rwo
     !     geopos(1) = E longitude of observing site (rad)
     !     geopos(2) = N latitude of observing site (rad)
     !     geopos(3) = Altitude of observing site (m)
     DOUBLE PRECISION, DIMENSION(:), POINTER :: geopos => NULL()
  END TYPE ast_obs
  ! Generic astrometric observation
  TYPE ast_obs_direct
     CHARACTER(LEN=name_len)         :: objdes     ! object designation
     CHARACTER(LEN=1)                :: type       ! observation type:
                                                   !    O = optical
                                                   !    R = radar range
                                                   !    V = radar range rate
                                                   !    S = satellite/space
     CHARACTER(LEN=1)                :: tech       ! observation technology (column 15 of MPC record):
                                                   !    P = Photographic
                                                   !    E = Encoder
                                                   !    C = CCD
                                                   !    T = Transit or meridian circle
                                                   !    M = Micrometer
                                                   !    c = corrected for center of mass (radar only)
                                                   !    s = surface bounce (radar only)
     CHARACTER(LEN=1)                :: par_unit   ! parallax unit for satellite observations:
                                                   !    1 = km
                                                   !    2 = au
     CHARACTER(LEN=1)                :: note       ! note on observation circumstances (column 14 of MPC record)
     CHARACTER(LEN=2)                :: catcodmpc  ! MPC catalog flag:
                                                   ! a = USNO-A1.0                               
                                                   ! b = USNO-SA1.0
                                                   ! c = USNO-A2.0
                                                   ! d = USNO-SA2.0
                                                   ! e = UCAC-1
                                                   ! f = Tycho-1
                                                   ! g = Tycho-2
                                                   ! h = GSC-1.0
                                                   ! i = GSC-1.1
                                                   ! j = GSC-1.2
                                                   ! k = GSC-2.2
                                                   ! l = ACT
                                                   ! m = GSC-ACT
                                                   ! n = TRC
                                                   ! o = USNO-B1.0
                                                   ! p = PPM
                                                   ! q = UCAC-2-beta
                                                   ! r = UCAC-2
                                                   ! s = USNO-B2.0
                                                   ! t = UCAC-3-beta
                                                   ! u = UCAC-3
                                                   ! v = NOMAD
                                                   ! w = CMC
                                                   ! x = Hip-2
                                                   ! z = GSC (unspecified version)
                                                   ! ' ' = unknown
     CHARACTER(LEN=7)                :: obscod_s   ! observer's MPC code (string, official); for radar, TX(A3)//' '//RX(A3)
     INTEGER                         :: obscod_i   ! observer's MPC code (integer, internal representation); for radar,
                                                   !  encoding of both transmitter ans receiver (see jplradar_transform)
     DOUBLE PRECISION                :: time_utc   ! observation UTC time (MJD)
     DOUBLE PRECISION                :: time_tdt   ! observation TDT time (MJD)
     DOUBLE PRECISION, DIMENSION(2)  :: coord      ! observation coordinates:
                                                   !    RA,DEC (optical & satellite) (rad)
                                                   !    range,range rate (radar)
     DOUBLE PRECISION                :: mag        ! object apparent magnitude
     CHARACTER(LEN=6)                :: mag_str    ! magnitude original string in the record
     CHARACTER(LEN=1)                :: mag_band   ! magnitude color band
     LOGICAL                         :: mag_def    ! TRUE if magnitude is defined
     DOUBLE PRECISION                :: acc_time   ! accuracy of time (as given in input format) (d)
     DOUBLE PRECISION, DIMENSION(2)  :: acc_coord  ! accuracy of coordinates (as given in input format)
     DOUBLE PRECISION                :: acc_mag    ! accuracy of magnitude (as given in input format)
     ! Observer's position and velocity are given with respect to the mean ecliptic and equinox J2000
     ! WARNING: in case of radar observations:
     !         obspos = body-fixed position of transmitting station
     !         obsvel = body-fixed position of receiving station
     DOUBLE PRECISION, DIMENSION(3)  :: obspos     ! observer's geocentric position (au)
     DOUBLE PRECISION, DIMENSION(3)  :: obsvel     ! observer's geocentric velocity (au/d)
  END TYPE ast_obs_direct

  !INTERFACE OPERATOR(=)
  !  MODULE PROCEDURE drop_geopos
  !END INTERFACE

  DOUBLE PRECISION, DIMENSION(3), PARAMETER :: zero_3d_vect = (/ 0.d0, 0.d0, 0.d0 /)

  ! Undefined astrometric observation (useful for initializing variables before assigning values)
  TYPE(ast_obs), PARAMETER :: undefined_ast_obs = AST_OBS( &
       &   ' ' ,                &  ! object designation
       &   '?' ,                &  ! observation type
       &   '?' ,                &  ! observation technology
       &   ' ' ,                &  ! parallax unit (only for satellite observations)
       &   ' ' ,                &  ! note on observation circumstances (column 14 of MPC record)
       &   '  ',                &  ! astrometric catalog code
       &   '       ' ,          &  ! observer's MPC code (official)
       &   -999999 ,            &  ! observer's MPC code (internal)
       &   1.d99 ,              &  ! observation UTC time
       &   1.d99 ,              &  ! observation TDT time (MJD)
       &   (/ 0.d0, 0.d0 /) ,   &  ! observation coordinates
       &   1.d99 ,              &  ! object apparent magnitude
       &   '      ' ,           &  ! magnitude original string in the record
       &   '?' ,                &  ! magnitude color band
       &   .false. ,            &  ! magnitude is not defined
       &   1.d99 ,              &  ! accuracy of time (as given in input format) (d)
       &   (/ 1.d99, 1.d99 /) , &  ! accuracy of coordinates (as given in input format)
       &   1.d99 ,              &  ! accuracy of magnitude (as given in input format)
       &   zero_3d_vect ,       &  ! observer's geocentric position (au)
       &   zero_3d_vect ,       &  ! observer's geocentric velocity (au/d)
       &   NULL()               )  ! observer's geodetic coordinates (longitude, latitude, altitude)

  ! Weighting, bias, selection flag and residuals of astrometric observation
  TYPE ast_wbsr
     INTEGER                        :: sel_coord  ! selection flag for coordinates:
                                                  !   0 = not used
                                                  !   1 = used in LS fit
                                                  !   2 = used in initial orbit determination + LS fit
     DOUBLE PRECISION, DIMENSION(2) :: rms_coord  ! RMS of coordinates to be used in weighting
     DOUBLE PRECISION               :: correl     ! correlation of (alpha*cos(delta),delta)
     DOUBLE PRECISION, DIMENSION(2) :: bias_coord ! biases of coordinates
     LOGICAL,          DIMENSION(2) :: force_w    ! TRUE if weight of coordinates has been forced
     DOUBLE PRECISION, DIMENSION(2) :: res_coord  ! residuals of coordinates
     LOGICAL                        :: resc_def   ! TRUE if residuals of coordinates are defined
     DOUBLE PRECISION               :: chi        ! 2D chi (SQRT(chisquare)) of residuals
     DOUBLE PRECISION               :: rms_mag    ! RMS of magnitude to be used in weighting
     INTEGER                        :: sel_mag    ! selection flag for magnitude:
                                                  !   0 = not used
                                                  !   1 = used in LS fit
     DOUBLE PRECISION               :: res_mag    ! residual of magnitude
     LOGICAL                        :: resm_def   ! TRUE if residual of magnitude is defined
  END TYPE ast_wbsr


  ! Undefined weight, bias, selection flag and residuals of astrometric observation
  ! (useful for initializing variables before assigning values)
  TYPE(ast_wbsr), PARAMETER :: undefined_ast_wbsr = AST_WBSR( &
       1 ,                      &  ! selection flag for coordinates
       (/ 9.d9, 9.d9 /) ,       &  ! RMS of coordinates to be used in weighting
       0.d0,                    &  ! correlation of (alpha*cos(delta),delta)
       (/ 0.d0, 0.d0 /) ,       &  ! biases of coordinates
       (/ .false., .false. /) , &  ! TRUE if weight of coordinates has been forced
       (/ 0.d0, 0.d0 /) ,       &  ! residuals of coordinates
       .false. ,                &  ! TRUE if residuals of coordinates are defined
       0.d0 ,                   &  ! 2D chi (SQRT(chisquare)) of residuals
       9.d9 ,                   &  ! RMS of magnitude to be used in weighting
       0 ,                      &  ! selection flag for magnitude
       0.d0 ,                   &  ! residual of magnitude
       .false.                  )  ! TRUE if residual of magnitude is define

  INTEGER, PARAMETER :: rwo_version = 2           !  current version of RWO file format

  INTEGER, PARAMETER, PUBLIC :: nradobsx=10000 ! max total number of radar observations
  CHARACTER*100, PUBLIC, DIMENSION(nradobsx) :: radar_rec  ! JPL formatted records

  ! list of valid technology descriptors according to MPC
  CHARACTER(LEN=1), PARAMETER :: tech_values(16)=(/"C","c","P","A","N","M","T","E","e","V","S","O","H","n","x","X"/)
  
  ! LIST OF PUBLIC ENTITIES
  ! Derived TYPEs
  PUBLIC :: ast_obs, ast_wbsr, ast_obs_direct
  ! PARAMETERs
  PUBLIC :: undefined_ast_obs, undefined_ast_wbsr
  ! DATA: only asteroid radius (it is considered part of the observation data)
  DOUBLE PRECISION radius
  PUBLIC :: radius

CONTAINS

  TYPE(ast_obs_direct) FUNCTION drop_geopos(obs)
    TYPE(ast_obs), INTENT(IN) :: obs
    drop_geopos%objdes=obs%objdes
    drop_geopos%type=obs%type
    drop_geopos%tech=obs%tech
    drop_geopos%par_unit=obs%par_unit
    drop_geopos%note=obs%note
    drop_geopos%obscod_s=obs%obscod_s
    drop_geopos%obscod_i=obs%obscod_i
    drop_geopos%time_utc=obs%time_utc
    drop_geopos%time_tdt=obs%time_tdt
    drop_geopos%coord=obs%coord
    drop_geopos%mag=obs%mag
    drop_geopos%mag_str=obs%mag_str
    drop_geopos%mag_band=obs%mag_band
    drop_geopos%mag_def=obs%mag_def
    drop_geopos%acc_time=obs%acc_time
    drop_geopos%acc_coord=obs%acc_coord
    drop_geopos%acc_mag=obs%acc_mag
    drop_geopos%obspos=obs%obspos
    drop_geopos%obsvel=obs%obsvel
  END FUNCTION drop_geopos
  TYPE(ast_obs) FUNCTION add_geopos(obs)
    TYPE(ast_obs_direct), INTENT(IN) :: obs
    add_geopos%objdes=obs%objdes
    add_geopos%type=obs%type
    add_geopos%tech=obs%tech
    add_geopos%par_unit=obs%par_unit
    add_geopos%note=obs%note
    add_geopos%obscod_s=obs%obscod_s
    add_geopos%obscod_i=obs%obscod_i
    add_geopos%time_utc=obs%time_utc
    add_geopos%time_tdt=obs%time_tdt
    add_geopos%coord=obs%coord
    add_geopos%mag=obs%mag
    add_geopos%mag_str=obs%mag_str
    add_geopos%mag_band=obs%mag_band
    add_geopos%mag_def=obs%mag_def
    add_geopos%acc_time=obs%acc_time
    add_geopos%acc_coord=obs%acc_coord
    add_geopos%acc_mag=obs%acc_mag
    add_geopos%obspos=obs%obspos
    add_geopos%obsvel=obs%obsvel
    add_geopos%geopos => NULL()
  END FUNCTION add_geopos


  ! SUBMODULE iorwo
  ! ==========================================
  ! Write a .rwo file- structured
  ! allows to select sorting or not sorting by observation time
  ! however, radar observations are always after optical
  ! ==========================================
  SUBROUTINE write_rwo(file,obs,obsw,n,error_model,rms,rmsmag,sort,unit0)

    CHARACTER(LEN=*),             INTENT(IN)           :: file        ! output file
    INTEGER,                      INTENT(IN)           :: n           ! number of observations to write
    TYPE(ast_obs),  DIMENSION(n), INTENT(IN)           :: obs         ! set of n observations
    TYPE(ast_wbsr), DIMENSION(n), INTENT(IN)           :: obsw        ! set of n weights and residuals
    CHARACTER(LEN=*),             INTENT(IN)           :: error_model ! error model file name 
    DOUBLE PRECISION,             INTENT(IN), OPTIONAL :: rms         ! RMS of astrometric fit
    DOUBLE PRECISION,             INTENT(IN), OPTIONAL :: rmsmag      ! RMS of photometric fit
    LOGICAL,                      INTENT(IN), OPTIONAL :: sort        ! sort times?
    ! note default is to sort
    INTEGER,                      INTENT(IN), OPTIONAL :: unit0 ! output unit if file is already opened
    INCLUDE 'parcmc.h90'     ! comment character
    INTEGER, DIMENSION(n) :: time_ord
    INTEGER :: unit,i,k,lf,nrec
    DOUBLE PRECISION :: rms0, rmsmag0
    LOGICAL :: error,radar
    INTEGER, EXTERNAL :: lench
    EXTERNAL :: filopn,filclo,mjddat,real2string,sessag
    DATA nrec /0/
    ! ===============================================
    ! check operating mode
    IF(PRESENT(unit0))THEN
       unit=unit0
    ELSE! open output file
       CALL filopn(unit,file,'UNKNOWN')
       nrec=0
       IF(PRESENT(rms)) THEN
          rms0=rms
       ELSE
          rms0=-1.
       ENDIF
       IF(PRESENT(rmsmag)) THEN
          rmsmag0=rmsmag
       ELSE
          rmsmag0=-1.
       ENDIF
       ! write header
       CALL write_rwo_header(unit, error_model, nrec, rms0, rmsmag0)
    ENDIF
    IF(PRESENT(sort))THEN
       IF(sort)CALL heapsort(obs(1:n)%time_tdt,n,time_ord)
    ELSE
       CALL heapsort(obs(1:n)%time_tdt,n,time_ord) !default is to sort
    ENDIF
    ! First pass: only optical observations (including roving observatory and satellite) are processed
    radar=.false.
    DO k=1,n
       IF(PRESENT(sort))THEN
          IF(sort)THEN
             i=time_ord(k)
          ELSE
             i=k
          ENDIF
       ELSE
          i=time_ord(k)
       ENDIF
       SELECT CASE (obs(i)%type)
       CASE ('O', 'S')                  ! OPTICAL/SATELLITE OBSERVATION
          CALL write_rwo_opt(unit,obs(i),obsw(i),error,nrec)
          IF(error)THEN
             lf=lench(file)
             WRITE(ierrou,121) file(1:lf)
121          FORMAT('write_rwo: error in DAY field: file ',A)
             ! numerr has already been incremeted inside write_rwo_opt
          ENDIF
          ! Second record for satellite observations
          IF(obs(i)%type == 'S'.or.obs(i)%obscod_s == '247') THEN
             CALL write_rwo_pos(unit,obs(i),obsw(i))
             nrec=nrec+1
          ENDIF
       CASE ('R', 'V')             ! RADAR OBSERVATION (delayed to next pass)
          radar=.TRUE.
       CASE DEFAULT
          lf=lench(file)
          WRITE(*,130) file(1:lf),obs(i)%type
          WRITE(ierrou,130) file(1:lf),obs(i)%type
130       FORMAT('write_rwo: cannot write in file ',A,' observation of unknown type ',A)
          numerr=numerr+1
       END SELECT
    END DO

    ! Second pass: radar observations
    IF(radar) THEN
       WRITE(unit,107) comcha,comcha
107    FORMAT(A1,' Object   Obser ====== Date =======  ============ Radar range/range rate (km or km/d) ============= ',  &
            'Station    Residual'/  &
            A1,' Design   K T N YYYY MM DD hh:mm:ss        Measure     Accuracy    rms    F      Bias       ',          &
            'Resid   TRX RCX     Chi   S')
       DO k=1,n
          IF(PRESENT(sort))THEN
             IF(sort)THEN
                i=time_ord(k)
             ELSE
                i=k
             ENDIF
          ELSE
             i=time_ord(k)
          ENDIF
          IF(obs(i)%type.eq.'R'.or.obs(i)%type.eq.'V')THEN
             CALL write_rwo_rad(unit,obs(i), obsw(i))
             nrec=nrec+1
          ENDIF
       END DO
    END IF
    IF(.not.PRESENT(unit0))CALL filclo(unit,' ')
  END SUBROUTINE write_rwo

  ! Writes one record of radar observations
  SUBROUTINE write_rwo_rad(unit, obs, obsw)
    INTEGER,        INTENT(IN) :: unit 
    TYPE(ast_obs),  INTENT(IN) :: obs
    TYPE(ast_wbsr), INTENT(IN) :: obsw
    INTEGER                    :: iday,month,year,iss,hh,mm 
    DOUBLE PRECISION           :: day,hour,rss
    CHARACTER(LEN=9)  :: chi_str
    CHARACTER(LEN=11) :: resid_rad,bias_rad
    CHARACTER(LEN=36) :: tmp2
    SELECT CASE (obs%type)
    CASE ('R', 'V')             ! RADAR OBSERVATION
       CALL mjddat(obs%time_utc,iday,month,year,hour)
       rss=hour*3600
       iss=NINT(rss)
       IF(iss<0)THEN
          WRITE(*,*) '**** write_rwo_rad: negative hour ****', hour, rss, iss
          STOP
       ELSEIF(ABS(rss-iss)>1.D-3)THEN
          WRITE(*,*) '**** write_rwo_rad: problem in rounding time ****',hour,rss,iss
          STOP
       ENDIF
       hh=iss/3600
       iss=iss-hh*3600
       mm=iss/60
       iss=iss-mm*60
       WRITE(tmp2,140) obs%objdes,obs%type,obs%tech,obs%note,  &
            &            year,month,iday,hh,mm,iss
140    FORMAT(1X,A9,1X,A1,1X,A1,1X,A1,1X,I4.4,1X,I2.2,1X,I2.2,1X,I2.2,':',I2.2,':',I2.2)
    END SELECT
    SELECT CASE (obs%type)
    CASE ('R')                  ! RADAR RANGE OBSERVATION
       bias_rad=radar_resid_string(obsw%bias_coord(1)*aukm,.TRUE.)
       resid_rad=radar_resid_string(obsw%res_coord(1)*aukm,obsw%resc_def)
       IF(obsw%resc_def) THEN
          IF(abs(obsw%chi).lt.9.9e5)THEN
             WRITE(chi_str,128) obsw%chi
          ELSE
             WRITE(chi_str,1128) obsw%chi
          ENDIF
       ELSE
          chi_str=' '
       END IF
       WRITE(unit,202) tmp2,obs%coord(1)*aukm,obs%acc_coord(1)*aukm,    &
            &                obsw%rms_coord(1)*aukm,obsw%force_w(1),   &
            &               bias_rad,resid_rad,obs%obscod_s,chi_str,obsw%sel_coord
    CASE ('V')                  ! RADAR RANGE RATE OBSERVATION
       bias_rad=radar_resid_string(obsw%bias_coord(2)*aukm,.TRUE.)
       resid_rad=radar_resid_string(obsw%res_coord(2)*aukm,obsw%resc_def)
       IF(obsw%resc_def) THEN
          IF(abs(obsw%chi).lt.9.9e5)THEN
             WRITE(chi_str,128) obsw%chi
128          FORMAT(F9.2)
          ELSE
             WRITE(chi_str,1128) obsw%chi
1128         FORMAT(1P,D9.2)
          ENDIF
       ELSE
          chi_str=' '
       END IF
       WRITE(unit,202) tmp2,obs%coord(2)*aukm,obs%acc_coord(2)*aukm,     &
            &                obsw%rms_coord(2)*aukm,obsw%force_w(2),   &
            &               bias_rad,resid_rad,obs%obscod_s,chi_str,obsw%sel_coord
202    FORMAT(A,F16.5,1X,F9.5,1X,F9.5,1X,L1,1X,A,1X,A,1X,A7,1X,A,1X,I1)
    END SELECT
  END SUBROUTINE write_rwo_rad

  ! writes one record of a .rwo file with an optical observation
  ! If satellite, roving only the first line
  SUBROUTINE write_rwo_opt(unit,obs,obsw,error,nrec)
    INTEGER,        INTENT(IN)    :: unit 
    TYPE(ast_obs),  INTENT(IN)    :: obs
    TYPE(ast_wbsr), INTENT(IN)    :: obsw
    LOGICAL,        INTENT(OUT)   :: error
    INTEGER,        INTENT(INOUT) :: nrec
    INTEGER                       :: iday,month,year,nd2,ra_h,ra_min,dec_deg,dec_min,hh,mm,iss 
    DOUBLE PRECISION              :: day,hour,ra_sec,resnor,dec_sec,rss
    CHARACTER(LEN=1)              :: ra_sign,dec_sign
    CHARACTER(LEN=5)              :: dec_sec_str
    CHARACTER(LEN=6)              :: ra_sec_str
    CHARACTER(LEN=8)              :: rmsa_str,biasa_str,rmsd_str,biasd_str
    CHARACTER(LEN=9)              :: resida_str,residd_str,chi_str
    CHARACTER(LEN=13)             :: sday
    CHARACTER(LEN=19)             :: mag_str
    CHARACTER(LEN=49)             :: tmp1
    INTEGER, EXTERNAL             :: lench
    ! check obs.type 
    IF(obs%type.ne.'O'.and.obs%type.ne.'S')THEN
       error=.true.
       WRITE(ierrou,223) obs%type, nrec+1
223    FORMAT('write_rwo: attempt to write as optical type ',A1,'  : record',I5)
       numerr=numerr+1
    ENDIF
    ! Conversion of time
    CALL mjddat(obs%time_utc,iday,month,year,hour)
    ! Compose common part of output string (designation, observation type, time)
    day=iday+hour/24.d0
    nd2=-NINT(LOG10(obs%acc_time))
    CALL real2string(day,2,nd2,sday,error)
    IF(error)THEN
       WRITE(ierrou,121) nrec+1
121    FORMAT('write_rwo: error in DAY field: record',I5)
       numerr=numerr+1
    ENDIF
    WRITE(tmp1,120) obs%objdes,obs%type,obs%tech,obs%note,year,month,sday,obs%acc_time
120 FORMAT(1X,A9,1X,A1,1X,A1,1X,A1,1X,I4.4,1X,I2.2,1X,A,1X,1P,E10.3,0P)
    ! Conversion of RA
    CALL sessag(obs%coord(1)*hrad,ra_sign,ra_h,ra_min,ra_sec)
    IF(ra_sign == '-')THEN
       WRITE(*,*) '**** write_rwo_opt: negative right ascension **** ',obs%objdes,' ',ra_sign,ra_h,ra_min,ra_sec
       WRITE(ierrou,*)   '**** write_rwo_opt: negative right ascension **** ',obs%objdes,' ',ra_sign,ra_h,ra_min,ra_sec
       numerr=numerr+1 
    ENDIF
    WRITE(ra_sec_str,122) ra_sec
122 FORMAT(F6.3)
    IF(ra_sec_str(1:1) == ' ') ra_sec_str(1:1)='0'
    WRITE(rmsa_str,125)  obsw%rms_coord(1)*secrad*COS(obs%coord(2))
    WRITE(biasa_str,125) obsw%bias_coord(1)*secrad*COS(obs%coord(2))
    !  WRITE(biasa_str,125) obsw%bias_coord(1)*secrad*COS(obs%coord(2)) wrong correction for cos(dec)
125 FORMAT(F8.3)
    IF(obsw%resc_def) THEN
       resnor=obsw%res_coord(1)*secrad*COS(obs%coord(2))
       IF(ABS(resnor) > 999.d0)THEN
          WRITE(resida_str,123) resnor
123       FORMAT(1P,E9.1,0P)
       ELSE
          WRITE(resida_str,124) resnor
124       FORMAT(F9.3)
       END IF
    ELSE
       resida_str=' '
    END IF
    ! Conversion of DEC
    CALL sessag(obs%coord(2)*degrad,dec_sign,dec_deg,dec_min,dec_sec)
    WRITE(dec_sec_str,126) dec_sec
    IF(dec_sec_str(1:1) == ' ') dec_sec_str(1:1)='0'
    WRITE(rmsd_str,125)  obsw%rms_coord(2)*secrad
    WRITE(biasd_str,125) obsw%bias_coord(2)*secrad
    IF(obsw%resc_def) THEN
       resnor=obsw%res_coord(2)*secrad
       IF(ABS(resnor) > 999.d0)THEN
          WRITE(residd_str,123) resnor
       ELSE
          WRITE(residd_str,124) resnor
       END IF
    ELSE
       residd_str=' '
    END IF
    ! Chi
    IF(obsw%resc_def) THEN
       IF(abs(obsw%chi).lt.9.9e5)THEN
          WRITE(chi_str,128) obsw%chi
128       FORMAT(F9.2)
       ELSE
          WRITE(chi_str,1128) obsw%chi
1128      FORMAT(1P,D9.2)
       ENDIF
    ELSE
       chi_str=' '
    END IF
    ! Conversion of MAG
    mag_str=' '
    IF(obs%mag_def) THEN
       mag_str(1:6)=obs%mag_str
       ! This is not necessary since accuracy can be determined from mag_str
       ! IF(obs%acc_mag >= 0.d0) WRITE(mag_str(8:12),126) obs%acc_mag
       IF(obsw%rms_mag >= 0.d0) WRITE(mag_str(8:12),126) obsw%rms_mag
126    FORMAT(F5.2)
       IF(obsw%resm_def) WRITE(mag_str(14:19),127) obsw%res_mag
127    FORMAT(F6.2)
    END IF
    IF(lench(obs%obscod_s)>3) THEN
       !     lf=lench(file)
       WRITE(ierrou,131) obs%obscod_s !,nrec+1,file(1:lf)
131    FORMAT('write_rwo: abnormal observatory name "',A)
       numerr=numerr+1
    END IF
    WRITE(unit,201) tmp1,ra_h,ra_min,ra_sec_str,obs%acc_coord(1)*secrad*COS(obs%coord(2)),   &
         &               rmsa_str,obsw%force_w(1),biasa_str,resida_str,            &
         &               dec_sign,dec_deg,dec_min,dec_sec_str,obs%acc_coord(2)*secrad, &
         &               rmsd_str,obsw%force_w(2),biasd_str,residd_str,            &
         &               mag_str,obs%catcodmpc,obs%obscod_s,chi_str,obsw%sel_coord,obsw%sel_mag

201 FORMAT(A,1X,I2.2,1X,I2.2,1X,A6,1X,1P,E10.3,0P,1X,A8,1X,L1,1X,A8,A9, &
         1X,A1,I2.2,1X,I2.2,1X,A5,1X,1P,E10.3,0P,1X,A8,1X,L1,1X,A8,A9, &
         1X,A,2X,A2,1X,A3,1X,A,1X,I1,1X,I1)
    nrec=nrec+1
  END SUBROUTINE write_rwo_opt

  ! Writes second record with position for satellite/roving observations
  SUBROUTINE write_rwo_pos(unit,obs,obsw)
    INTEGER,        INTENT(IN) :: unit 
    TYPE(ast_obs),  INTENT(IN) :: obs
    TYPE(ast_wbsr), INTENT(IN) :: obsw
    DOUBLE PRECISION, DIMENSION(3)   :: equpos
    CHARACTER(LEN=33) :: tmp3
    CHARACTER(LEN=49) :: tmp1
    DOUBLE PRECISION  :: conv
    CHARACTER(LEN=13) :: sday
    INTEGER           :: iday,month,year,nd2
    DOUBLE PRECISION  :: day,hour
    LOGICAL           :: error
    ! Conversion of time
    CALL mjddat(obs%time_utc,iday,month,year,hour)
    ! Compose common part of output string (designation, observation type, time)
    day=iday+hour/24.d0
    nd2=-NINT(LOG10(obs%acc_time))
    CALL real2string(day,2,nd2,sday,error)
    WRITE(tmp1,120) obs%objdes,obs%type,obs%tech,obs%note,year,month,sday
120 FORMAT(1X,A9,1X,A1,1X,A1,1X,A1,1X,I4.4,1X,I2.2,1X,A)
    IF(obs%type == 'S') THEN
       tmp3=tmp1(1:33)
       tmp3(14:14)='s'
       SELECT CASE (obs%par_unit)
       CASE ('1')
          conv=aukm
       CASE ('2')
          conv=1.d0
       CASE DEFAULT
          STOP '**** write_rwo: internal error (04) ****'
       END SELECT
       ! Transformation of ecliptical coordinates into equatorial
       !         CALL rotpn(rot,'ECLM','J2000',0.d0,'EQUM','J2000',0.d0)
       equpos=MATMUL(roteceq,obs%obspos)
       equpos=equpos*conv
       ! Write parallax unit and satellite position
       WRITE(unit,205) tmp3,obs%par_unit,equpos,obs%obscod_s(1:3)
205    FORMAT(A,1X,A1,3F24.12,1X,A)
    END IF
    ! Second record for roving observatory observations
    IF(obs%obscod_s == '247') THEN
       IF(obs%tech /= 'V') STOP '**** write_rwo: internal error (05) ****'
       IF(.NOT.ASSOCIATED(obs%geopos)) STOP '**** write_rwo: internal error (06) ****'
       tmp3=tmp1(1:33)
       tmp3(14:14)='v'
       WRITE(unit,206) tmp3,obs%geopos(1:2)*degrad,obs%geopos(3),obs%obscod_s(1:3)
206    FORMAT(A,1X,F10.6,1X,F10.6,1X,F8.1,1X,A)
    END IF
  END SUBROUTINE write_rwo_pos

  CHARACTER(LEN=11) FUNCTION radar_resid_string(value,defined)
    DOUBLE PRECISION, INTENT(IN) :: value
    LOGICAL,          INTENT(IN) :: defined
    LOGICAL                      :: fixed
    radar_resid_string=' '
    IF(.NOT.defined) RETURN
    fixed=.true.
    IF(value>99999.D0) fixed=.false.
    IF(value<-9999.d0) fixed=.false.
    IF(fixed) THEN
       WRITE(radar_resid_string,101) value
101    FORMAT(F11.5)
    ELSE
       WRITE(radar_resid_string,102) value
102    FORMAT(1P,E11.4)
    END IF
  END FUNCTION radar_resid_string

  SUBROUTINE write_rwo_header(unit, error_model, nrec, rms, rmsmag)
    INTEGER,          INTENT(IN) :: unit !output unit
    CHARACTER(LEN=*), INTENT(IN) :: error_model ! error model file name
    INTEGER,       INTENT(INOUT) :: nrec ! record counter 
    DOUBLE PRECISION, INTENT(IN), OPTIONAL :: rms         ! RMS of astrometric fit
    DOUBLE PRECISION, INTENT(IN), OPTIONAL :: rmsmag      ! RMS of photometric fit
    INTEGER lench, ll
    INCLUDE 'parcmc.h90'     ! comment character
    ! write header

    WRITE(unit,101) rwo_version
101 FORMAT('version = ',I3)
    nrec=nrec+1
    ll=lench(error_model)
    IF(ll > 0) THEN
       WRITE(unit,102) error_model(1:ll)
102    FORMAT('errmod  = ''',A,'''')
    ELSE
       WRITE(unit,102) ' '
    END IF
    nrec=nrec+1
    IF(PRESENT(rms))THEN
       IF(rms >= 0.d0) THEN
          WRITE(unit,103) rms
103       FORMAT('RMSast  = ',1P,E13.5)
          nrec=nrec+1
       ENDIF
    ENDIF
    IF(PRESENT(rmsmag))THEN
       IF(rmsmag >= 0.d0) THEN
          WRITE(unit,104) rmsmag
104       FORMAT('RMSmag  = ',1P,E13.5)
          nrec=nrec+1
       ENDIF
    ENDIF
    WRITE(unit,105)
105 FORMAT('END_OF_HEADER')
    WRITE(unit,106) comcha,comcha
    nrec=nrec+3
106 FORMAT(A1,' Object   Obser ============= Date ============= ================== Right Ascension =================', &
         '  ================= Declination ===================== ==== Magnitude ==== Ast Obs  Residual SEL'/  &
         A1,' Design   K T N YYYY MM DD.dddddddddd   Accuracy HH MM SS.sss  Accuracy      RMS  F     Bias    Resid sDD MM SS.ss', &
         '  Accuracy      RMS  F     Bias    Resid Val  B   RMS  Resid Cat Cod       Chi A M')

  END SUBROUTINE write_rwo_header

  ! input from rwo files - structured
  SUBROUTINE read_rwo(file,obs,obsw,n,error_model,rms,rmsmag,unit0,eof0)
    ! INTERFACE
    CHARACTER(LEN=*),             INTENT(IN)            :: file        ! input file, to be opened
    ! unless the optional variable unit0 is present
    TYPE(ast_obs),  DIMENSION(:), INTENT(OUT)           :: obs         ! observations
    TYPE(ast_wbsr), DIMENSION(:), INTENT(OUT)           :: obsw        ! observations weights/residuals
    INTEGER,                      INTENT(OUT)           :: n           ! number of observations found
    CHARACTER(LEN=*),             INTENT(OUT)           :: error_model ! error model file name 
    DOUBLE PRECISION,             INTENT(OUT), OPTIONAL :: rms         ! RMS of astrometric fit
    DOUBLE PRECISION,             INTENT(OUT), OPTIONAL :: rmsmag      ! RMS of photometric fit
    INTEGER,                      INTENT(IN) , OPTIONAL :: unit0 ! input unit if file is already opened
    ! then the routine must read until the obj.designation changes, then leave open
    LOGICAL,                      INTENT(OUT), OPTIONAL :: eof0 ! end of file if file was opened
    ! note eof0 is required if unit0 is present
    INCLUDE 'parobx.h90' ! observation numbers: maximum
    ! END INTERFACE
    LOGICAL error, eof, stored
    INTEGER :: unit,nr,kr,lf,j
    TYPE(ast_obs)  :: obs1
    TYPE(ast_wbsr) :: obsw1
    INTEGER, EXTERNAL :: lench
    CHARACTER*(name_len) ::  name0, name 
    SAVE obs1,obsw1,stored
    DATA stored /.false./
    ! ==========================================================
    ! check operating mode
    lf=lench(file)
    IF(PRESENT(unit0))THEN
       IF(.not.PRESENT(eof0))THEN
          WRITE(*,302) file(1:lf)
302       FORMAT('read_rwo: missing eof0 to warn of end of file ',A)
          STOP
       ENDIF
       IF(PRESENT(rms))THEN
          rms=-99.d0
       ENDIF
       IF(PRESENT(rmsmag))THEN
          rmsmag=-99.d0
       ENDIF
       unit=unit0
    ELSE
       ! open file
       CALL filopn(unit,file,'OLD')
       n=0 ! observations counter
       CALL rdfnam(unit,file,nr)
       !read header; should count records...
       CALL read_rwo_header(unit, error_model, error, nr, rms, rmsmag)
       IF(error)THEN
          WRITE(ierrou,301) file(1:lf)
301       FORMAT('read_rwo: input error in header of file ',A)
       ENDIF
    ENDIF
    ! main loop on observations
    DO j=1,nobx
       IF(.not.stored)THEN
          CALL read_rwo_rec(unit,nr,obs1,obsw1,error,eof)
          ! redifine rejection when tech is equal 'x'
          IF(obs1%tech.EQ.'x'.OR.obs1%tech.EQ.'X') THEN
             obsw1%sel_coord=0
          END IF
          IF(PRESENT(eof0))THEN
             eof0=eof
          ENDIF
          IF(eof)THEN
             stored=.false.
             EXIT
          ELSEIF(error)THEN 
             WRITE(*,202) nr,file(1:lf)
             IF(ierrou > 0) WRITE(ierrou,202) nr,file(1:lf)
202          FORMAT('read_rwo: input error at record',I6,' of file ',A)
             numerr=numerr+1
             stored=.false.
          ENDIF
       ELSE
          IF(PRESENT(eof0))THEN
             eof0=.false.
          ENDIF
       ENDIF
       IF(PRESENT(unit0))THEN
          IF(j.eq.1)THEN
             name0=obs1%objdes
             stored=.false.
          ELSE
             ! check if the designation has changed
             name=obs1%objdes
             IF(name.ne.name0)THEN
                ! store for next call
                stored=.true.                 
                EXIT
             ELSE
                stored=.false.
             ENDIF
          ENDIF
       ENDIF
       ! observation has been read
       n=n+1
       obs(n)=obs1
       obsw(n)=obsw1
    END DO
    IF(.not.PRESENT(unit0))CALL filclo(unit,' ')
  END SUBROUTINE read_rwo

  ! ======================================
  ! read_rwo_rec
  ! ======================================
  SUBROUTINE read_rwo_rec(unit,nr,obs1,obsw1,error, eof)
    USE station_coordinates
    USE reference_systems
    INTEGER,         INTENT(IN)   :: unit      ! input unit
    INTEGER,         INTENT(INOUT):: nr        ! record counter
    TYPE(ast_obs),   INTENT(OUT)  :: obs1      ! observation
    TYPE(ast_wbsr),  INTENT(OUT)  :: obsw1     ! weights residuals etc.
    LOGICAL,         INTENT(OUT)  :: error,eof ! error flag, end of file
    INCLUDE 'parcmc.h90'                       ! comment character
    CHARACTER(LEN=200)            :: record
    DOUBLE PRECISION              :: day,sec,sect,ra_sec,dec_sec,coord1,acc_coord1,rms_coord1,conv
    DOUBLE PRECISION, DIMENSION(3):: equpos,geopos,bfpos0
    CHARACTER(LEN=1)              :: dec_sign,col1,col2
    CHARACTER(LEN=3)              :: scale,obscod_s1
    CHARACTER(LEN=5)              :: mag_field
    CHARACTER(LEN=8)              :: rmsa_str,rmsd_str
    CHARACTER(LEN=9)              :: resida_str,residd_str,chi_str
    CHARACTER(LEN=11)             :: resid_rad,bias_rad
    CHARACTER(LEN=19)             :: mag_str
    CHARACTER(LEN=33)             :: tmp3
    LOGICAL                       :: force1
    INTEGER                       :: year,month,iday,mjd,mjdt,ra_h,ra_min,dec_deg,dec_min
    INTEGER                       :: len_mag_field,pos_point,n_dec_mag
    INTEGER                       :: hh,mm,iss,ipr0,ipr1,code_tx,code_rx
    DOUBLE PRECISION, EXTERNAL    :: tjm1
    CHARACTER(LEN=16)             :: stname
    CHARACTER(LEN=200)            :: tem_rec
    INTEGER                       :: len_comp,vec_comp,i0,len_aux,i  ! variables for splitting string of satellite position

    ! ======================================
    ! setup of flags
    error=.FALSE.
    eof=.FALSE.
    ! read one record
3   CONTINUE
    READ(unit,201,END=2) record
201 FORMAT(A)
    nr=nr+1
    IF(record(1:1) == comcha) GOTO 3
    obs1=undefined_ast_obs
    obsw1=undefined_ast_wbsr
    READ(record,110, ERR=1, END=2) obs1%objdes,obs1%type,obs1%tech,obs1%note
110 FORMAT(1X,A9,1X,A1,1X,A1,1X,A1)
    SELECT CASE (obs1%type)
    CASE ('O', 'S')                  ! OPTICAL/SATELLITE OBSERVATION
       ! time of observations
       READ(record(18:),120,ERR=1) year,month,day,obs1%acc_time
120    FORMAT(I4,1X,I2,1X,F13.10,1X,E10.3)
       iday=day
       sec=(day-iday)*86400.d0
       mjd=NINT(tjm1(iday,month,year,0.d0))
       IF(year < 1972) THEN
          scale='UT1'
       ELSE
          scale='UTC'
       ENDIF
       CALL cnvtim(mjd,sec,scale,mjdt,sect,'TDT')
       obs1%time_utc=mjd+sec/86400.d0
       obs1%time_tdt=mjdt+sect/86400.d0
       ! change_format='INPUT'
       IF(change_format.ne.'INPUT')THEN
          ! read alpha, delta with weights, bias, residuals... also magnitude, chi, selection flags
          READ(record(51:),521,ERR=1) ra_h,ra_min,ra_sec,obs1%acc_coord(1),rmsa_str,obsw1%force_w(1),&
               &              obsw1%bias_coord(1),resida_str,                                              &
               &              dec_sign,dec_deg,dec_min,dec_sec,obs1%acc_coord(2),rmsd_str,obsw1%force_w(2),&
               &              obsw1%bias_coord(2),residd_str,                                              &
               &              mag_str,obs1%catcodmpc,obs1%obscod_s,chi_str,obsw1%sel_coord,obsw1%sel_mag
521       FORMAT(I2,1X,I2,1X,F6.6,1X,E10.3,1X,A8,1X,L1,1X,F8.3,A9, &
               &           1X,A1,I2,1X,I2,1X,F5.2,1X,E10.3,1X,A8,1X,L1,1X,F8.3,A9, &
               &           1X,A19,2X,A2,1X,A3,1X,A9,1X,I1,1X,I1)
       ELSE
          ! read alpha, delta with weights, bias, residuals... also magnitude, chi, selection flags, no catcodmpc
          READ(record(51:),121,ERR=1) ra_h,ra_min,ra_sec,obs1%acc_coord(1),rmsa_str,obsw1%force_w(1),&
               &              obsw1%bias_coord(1),resida_str,                                              &
               &              dec_sign,dec_deg,dec_min,dec_sec,obs1%acc_coord(2),rmsd_str,obsw1%force_w(2),&
               &              obsw1%bias_coord(2),residd_str,                                              &
               &              mag_str,obs1%obscod_s,chi_str,obsw1%sel_coord,obsw1%sel_mag
121       FORMAT(I2,1X,I2,1X,F6.6,1X,E10.3,1X,A8,1X,L1,1X,F8.3,A9, &
               &           1X,A1,I2,1X,I2,1X,F5.2,1X,E10.3,1X,A8,1X,L1,1X,F8.3,A9, &
               &           1X,A19,1X,A3,1X,A9,1X,I1,1X,I1)
       ENDIF
       obs1%acc_coord=obs1%acc_coord/secrad
       IF(rmsa_str == ' ') THEN
          obsw1%rms_coord(1)=-1.d0
       ELSE
          READ(rmsa_str,*,ERR=1) obsw1%rms_coord(1)
          obsw1%rms_coord(1)=obsw1%rms_coord(1)/secrad
       END IF
       IF(rmsd_str == ' ') THEN
          obsw1%rms_coord(2)=-1.d0
       ELSE
          READ(rmsd_str,*,ERR=1) obsw1%rms_coord(2)
          obsw1%rms_coord(2)=obsw1%rms_coord(2)/secrad
       END IF
       obsw1%bias_coord=obsw1%bias_coord/secrad
       obs1%coord(1)=15.d0*(ra_h*3600.d0+ra_min*60.d0+ra_sec)/secrad
       obs1%coord(2)=(dec_deg*3600.d0+dec_min*60.d0+dec_sec)/secrad
       ! correct for cos(DEC)
       obs1%acc_coord(1)=obs1%acc_coord(1)/COS(obs1%coord(2))
       obsw1%bias_coord(1)=obsw1%bias_coord(1)/COS(obs1%coord(2))
       obsw1%rms_coord(1)=obsw1%rms_coord(1)/COS(obs1%coord(2))
       IF(dec_sign == '-') obs1%coord(2)=-obs1%coord(2)
       obsw1%resc_def=(resida_str /= ' ' .AND. residd_str /= ' ')
       IF(obsw1%resc_def) THEN
          READ(resida_str,*,ERR=1) obsw1%res_coord(1)
          READ(residd_str,*,ERR=1) obsw1%res_coord(2)
          obsw1%res_coord=obsw1%res_coord/secrad
          obsw1%res_coord(1)=obsw1%res_coord(1)/COS(obs1%coord(2))
          READ(chi_str,*,ERR=1) obsw1%chi
       ELSE
          obsw1%chi=-1.d0
       END IF
       ! station
       CALL statcode(obs1%obscod_s,obs1%obscod_i)
       ! magnitude
       obs1%mag_str=mag_str(1:6)
       mag_field=mag_str(1:5)
       CALL rmsp(mag_field,len_mag_field)
       mag_field=mag_str(1:5) !WARNING: to fix 644 mistake
       IF(len_mag_field>0) THEN
          READ(mag_field,*,ERR=1) obs1%mag
          obs1%mag_def=.true.
          pos_point=INDEX(mag_field,'.')
          IF(pos_point == 0) THEN
             obs1%acc_mag=1
          ELSE
             n_dec_mag=len_mag_field-pos_point
             obs1%acc_mag=10.0d0**(-n_dec_mag)
          END IF
       ELSE
          obs1%mag_def=.false.
          obs1%mag=0.d0
          obs1%acc_mag=99.d9
       END IF
       obs1%mag_band=mag_str(6:6)
       IF(mag_str(8:12) == ' ') THEN
          obsw1%rms_mag=-1.d0
       ELSE
          READ(mag_str(8:12),*,ERR=1) obsw1%rms_mag
       END IF
       IF(mag_str(14:19) == ' ') THEN
          obsw1%resm_def=.false.
          obsw1%res_mag=0.d0
       ELSE
          READ(mag_str(14:19),*,ERR=1) obsw1%res_mag
          obsw1%resm_def=.true.
       END IF
       ! Observer's position is computed for ground-based observations and read from the following record
       ! for roving observatory and satellite observations
       SELECT CASE (obs1%type)
       CASE ('O')
          SELECT CASE (obs1%obscod_s)
          CASE ('247')
             IF(obs1%tech /= 'V') GOTO 1
             READ(unit,206,ERR=1) tmp3,geopos,obscod_s1
206          FORMAT(A33,1X,F10.6,1X,F10.6,1X,F8.1,1X,A)
             nr=nr+1
             IF(tmp3(14:14) /= 'v') GOTO 1
             tmp3(14:14)='V'
             IF(tmp3 /= record(1:33)) GOTO 1
             IF(obscod_s1 /= obs1%obscod_s(1:3)) GOTO 1
             geopos(1:2)=geopos(1:2)*radeg
             IF(ASSOCIATED(obs1%geopos)) DEALLOCATE(obs1%geopos)
             ALLOCATE(obs1%geopos(3))
             obs1%geopos=geopos
             CALL geodetic_to_cartesian(geopos(1),geopos(2),geopos(3),bfpos0)
             bfpos0=bfpos0/(aukm*1d3)
             CALL observer_position(obs1%time_tdt,obs1%obspos,obs1%obsvel,BFPOS=bfpos0)
          CASE DEFAULT
             CALL observer_position(obs1%time_tdt,obs1%obspos,obs1%obsvel,OBSCODE=obs1%obscod_i)
          END SELECT
       CASE ('S')
          ! Read parallax unit and satellite position
          READ(unit,205,ERR=1) tmp3,obs1%par_unit,tem_rec
205       FORMAT(A33,1X,A1,A)
          ! Split string of satellite position with respect to blanks
          len_comp=0
          vec_comp=0
          len_aux=LEN_TRIM(tem_rec)
          DO i=1,len_aux-3
             IF(tem_rec(i:i).EQ.'-'.OR.tem_rec(i:i).EQ.'+')THEN
                i0=i
                len_comp=len_comp+1
             ELSEIF(isnumber(tem_rec(i:i)).OR.tem_rec(i:i).EQ.'.')THEN
                len_comp=len_comp+1
                IF(len_comp.EQ.1) i0=i
                IF(i.EQ.len_aux-3)THEN
                   vec_comp=vec_comp+1
                   READ(tem_rec(i0:i0+len_comp-1),*) equpos(vec_comp)
                   len_comp=0
                END IF
             ELSEIF(tem_rec(i:i).EQ.' ')THEN
                IF(len_comp.GT.0)THEN
                   vec_comp=vec_comp+1
                   READ(tem_rec(i0:i0+len_comp-1),*) equpos(vec_comp)
                   len_comp=0
                END IF
             END IF
          END DO
          obscod_s1 = tem_rec(len_aux-2:len_aux)
          nr=nr+1
          tmp3(14:14)='S'
          IF(tmp3 /= record(1:33)) GOTO 1
          IF(obscod_s1 /= obs1%obscod_s(1:3)) GOTO 1
          SELECT CASE (obs1%par_unit)
          CASE ('1')
             conv=1.d0/aukm
          CASE ('2')
             conv=1.d0
          CASE DEFAULT
             GOTO 1
          END SELECT
          equpos=equpos*conv
          ! Transformation of equatorial coordinates into ecliptical
          obs1%obspos=MATMUL(roteqec,equpos)
          obs1%obsvel=0.d0
       CASE DEFAULT
          STOP '**** read_rwo: internal error (01) ****'
       END SELECT

    CASE ('R', 'V')                  ! RADAR OBSERVATION (RANGE OR RANGE RATE)
       READ(record(18:),122,ERR=1) year,month,iday,hh,col1,mm,col2,iss
122    FORMAT(I4,1X,I2,1X,I2,1X,I2,A1,I2,A1,I2)
       IF(col1 /= ':') GOTO 1
       IF(col2 /= ':') GOTO 1
       mjd=NINT(tjm1(iday,month,year,0.d0))
       sec=iss+mm*60+hh*3600
       IF(year < 1972) THEN
          scale='UT1'
       ELSE
          scale='UTC'
       ENDIF
       CALL cnvtim(mjd,sec,scale,mjdt,sect,'TDT')
       obs1%time_utc=mjd+sec/86400.d0
       obs1%time_tdt=mjdt+sect/86400.d0
       READ(record(37:),123,ERR=1) coord1,acc_coord1,rms_coord1,force1, &
            &       bias_rad,resid_rad,obs1%obscod_s,chi_str,obsw1%sel_coord
123    FORMAT(F16.5,1X,F9.5,1X,F9.5,1X,L1,1X,A,1X,A,1X,A7,1X,A,1X,I1)
       CALL statcode(obs1%obscod_s(1:3),code_tx)
       CALL statcode(obs1%obscod_s(5:7),code_rx)
       obs1%obscod_i=code_tx*10000+code_rx
       CALL obscoo(code_tx,obs1%obspos,stname)
       CALL obscoo(code_rx,obs1%obsvel,stname)
       obs1%mag_def=.false.
       SELECT CASE (obs1%type)
       CASE ('R')
          ipr1=1
       CASE ('V')
          ipr1=2
       END SELECT
       ipr0=3-ipr1
       obs1%coord(ipr1)=coord1/aukm
       obs1%coord(ipr0)=0
       obs1%acc_coord(ipr1)=acc_coord1/aukm
       obsw1%rms_coord(ipr1)=rms_coord1/aukm
       obsw1%force_w(ipr1)=force1
       obsw1%force_w(ipr0)=.false.
       IF(bias_rad == ' ') THEN
          obsw1%bias_coord(ipr1)=0.d0
       ELSE
          READ(bias_rad,*,ERR=1) obsw1%bias_coord(ipr1)
          obsw1%bias_coord(ipr1)=obsw1%bias_coord(ipr1)/aukm
       END IF
       obsw1%bias_coord(ipr0)=0.d0
       obsw1%resc_def=(resid_rad /= ' ')
       IF(obsw1%resc_def) THEN
          READ(resid_rad,*,ERR=1) obsw1%res_coord(ipr1)
          obsw1%res_coord(ipr1)=obsw1%res_coord(ipr1)/aukm
       ELSE
          obsw1%res_coord(ipr1)=0.d0
       END IF
       obsw1%res_coord(ipr0)=0.d0
       IF(chi_str /= ' ') THEN
          READ(chi_str,*,ERR=1) obsw1%chi
       ELSE
          obsw1%chi=-9.d99
       END IF
       obs1%acc_time=10d-10
    CASE DEFAULT
       WRITE(*,203) obs1%type,nr
       WRITE(ierrou,203) obs1%type,nr
203    FORMAT('read_rwo: unknown observation type ',A,' at record',I6,' of file ',A)
       numerr=numerr+1
       error=.true.
    END SELECT
    RETURN
1   CONTINUE ! error cases
    error=.true.
    RETURN
2   CONTINUE ! regular end of file
    eof=.true.
  END SUBROUTINE read_rwo_rec

  ! ========================================
  !  read_rwo_header
  ! ========================================
  SUBROUTINE  read_rwo_header(unit, error_model, error, nr, rms, rmsmag)
    INTEGER,          INTENT(IN)  :: unit               !output unit
    CHARACTER(LEN=*), INTENT(OUT) :: error_model        ! error model file name
    LOGICAL,          INTENT(OUT) :: error              ! error flag
    INTEGER,          INTENT(OUT) :: nr               ! record counter
    DOUBLE PRECISION, INTENT(OUT), OPTIONAL :: rms      ! RMS of astrometric fit
    DOUBLE PRECISION, INTENT(OUT), OPTIONAL :: rmsmag   ! RMS of photometric fit
    EXTERNAL :: rdfnam,rdfint,rdfcha,rdfrea,cnvtim,rmsp
    INTEGER :: version,kr
    LOGICAL :: found
    ! ========================================
    error=.false.
    nr=0
    CALL rdfint(unit,'version',.false.,version,found,kr)
    IF(found) THEN
       nr=nr+1
       IF(version == 1) THEN
          change_format='INPUT'
       ELSEIF(version == rwo_version) THEN
          change_format=' '
       ELSE
          WRITE(*,200) version
200       FORMAT('Unknown version ',I2)
          WRITE(ierrou,200) version
          numerr=numerr+1
          error=.true.
          RETURN
       END IF
    ELSE
       WRITE(*,*) 'no rwo version  in header'
       WRITE(*,*) 'may be old f77 format '
       WRITE(ierrou,*) 'no rwo version  in header'
       numerr=numerr+1
       error=.true.
       RETURN
    END IF
    error_model=' '
    CALL rdfcha(unit,'errmod',.false.,error_model,found,kr)
    IF(.not.found)THEN
       WRITE(*,222) version
222    FORMAT('No error_model, version ',I2)
       WRITE(ierrou,200) version
       numerr=numerr+1
       error=.true.
       RETURN
    ENDIF
    nr=nr+1
    IF(PRESENT(rms)) THEN
       rms=-99.d0
       CALL rdfrea(unit,'RMSast',.false.,rms,found,kr)
       If(found) nr=nr+1
    END IF
    IF(PRESENT(rmsmag)) THEN
       rmsmag=-99.d0
       CALL rdfrea(unit,'RMSmag',.false.,rmsmag,found,kr)
       If(found) nr=nr+1
    END IF
  END SUBROUTINE read_rwo_header
  ! END SUBMODULE iorwo

  ! SUBMODULE iompc
  ! read file of MPC observation data

  SUBROUTINE mpc_obs_input(obs_found,obs,nobs,unit,filnam,eof,ons_mode,disc)

    LOGICAL, INTENT(OUT) :: obs_found ! true if there is something
    INTEGER, INTENT(OUT) :: nobs ! number of observations
    TYPE(ast_obs), INTENT(OUT), DIMENSION(:) :: obs ! vector of observations
    INTEGER, INTENT(IN), OPTIONAL :: unit ! input unit if file is already opened
    ! then the routine must read until the obj.designation changes, then leave open
    CHARACTER(*), INTENT(IN), OPTIONAL :: filnam ! input file name if it is to be opened

    ! then the routine must open, read until end of file, then close
    LOGICAL, INTENT(OUT), OPTIONAL :: eof ! end of file flag in case the file was open
    LOGICAL, INTENT(IN), OPTIONAL :: ons_mode ! if present, reading an ONS
    CHARACTER*1, INTENT(OUT), DIMENSION(:), OPTIONAL :: disc
    ! requiring a new designations manufactured by the perl routine
    !end interface

    INTEGER  iunit, nobsx, n, lf
    INTEGER, EXTERNAL :: lench 

    ! record and logical controls for mpcrec_transform
    CHARACTER*(80) mpcrec,mpcrec2
    LOGICAL skip_rad, error, complete, rad_in_mpcfile

    Character*9 name ! for message

    ! check that one or the other operations mode is selected (replacing either mpcin or mpcin3)

    obs_found=.FALSE.
    rad_in_mpcfile=.FALSE.

    nobs=0
    IF(PRESENT(unit).and.PRESENT(filnam))THEN
       WRITE(*,*) ' mpc_obs_input: both options present'
       STOP
    ENDIF

    nobsx=SIZE(obs,1)

    ! if necessary open file
    IF(PRESENT(filnam))THEN
       CALL filopn(iunit,filnam,'old')
    ELSEIF(PRESENT(unit))THEN
       iunit=unit
       IF(.not.PRESENT(eof))THEN
          WRITE(*,*)' mpc_obs_input: inconsistent optional arguments'
       ENDIF
    ELSE
       WRITE(*,*) ' mpc_obs_input: no options present'
       STOP
    ENDIF

    nobs=0
    DO n=1,nobsx
       nobs=nobs+1
       READ(iunit,'(a)',end=10, err=11) mpcrec
       ! blank line is like file termination
       IF(lench(mpcrec).eq.0) GOTO 10
       ! interpret one record
       IF(PRESENT(disc)) disc(nobs) = mpcrec(13:13)
       IF(PRESENT(ons_mode))THEN 
          CALL mpcrec_transform(mpcrec,obs(nobs),complete,error,skip_rad,ons_mode)
       ELSE  
          CALL mpcrec_transform(mpcrec,obs(nobs),complete,error,skip_rad)
       ENDIF
       IF(n.eq.1)name=obs(1)%objdes
       ! if error skip one record and hope it can go ahead
       IF(error) THEN
          WRITE(ierrou,*)'object ',name,'error at MPC observation record no ',n
          numerr=numerr+1
          nobs=nobs-1
          CYCLE
       ENDIF
       ! skipping MPC radar data and issue a warning
       IF(skip_rad)THEN
          nobs=nobs-1
          rad_in_mpcfile=.TRUE.
          CYCLE
       ENDIF
       ! if asteroid is changed, go back of the only record read so far
       IF(PRESENT(unit))THEN
          IF(obs(nobs)%objdes.ne.obs(1)%objdes)THEN
             eof=.FALSE.
             nobs=nobs-1
             BACKSPACE(iunit)
             RETURN
          ENDIF
       ENDIF
       ! handle cases not including satellite and roving observer observations
       IF(complete)THEN
          ! find observer position
          CALL mpcrec_add_obspos(obs(nobs),error)
          ! if error skip the two records and hope it can go ahead
          IF(error) THEN
             WRITE(ierrou,*)'object ',name,'error in position at MPC record no ',n
             numerr=numerr+1
             nobs=nobs-1
             CYCLE
          ENDIF
       ELSE
          ! handle roving observer/satellite position
          IF(obs(nobs)%type.eq.'S'.or.obs(nobs)%obscod_s.eq.'247')THEN
             ! read second line with satellite/roving obs. position
             READ(iunit,'(a)',end=10) mpcrec2
             IF(mpcrec2(15:15).eq.'X')THEN
                ! ignore X record if inserted between S and s
                READ(iunit,'(a)',end=10) mpcrec2
                CALL mpcrec_add_obspos(obs(nobs),error,mpcrec,mpcrec2)
             ELSEIF(mpcrec2(15:15).eq.'s'.or.mpcrec2(15:15).eq.'v')THEN
                !            WRITE(*,*)' processing satellite/roving observer data, t=',obs(nobs)%time_utc
                CALL mpcrec_add_obspos(obs(nobs),error,mpcrec,mpcrec2)
             ELSE
                error=.true.
             ENDIF
             ! if error skip the two records and hope it can go ahead
             IF(error) THEN
                IF(obs(nobs)%type.eq.'S')THEN
                   WRITE(ierrou,*)'object ',name,'error in satellite position at MPC record no ',n
                ELSEIF(obs(nobs)%obscod_s.eq.'247')THEN
                   WRITE(ierrou,*)'object ',name,'error in roving observer pos. at MPC record no ',n
                ENDIF
                numerr=numerr+1
                nobs=nobs-1
                CYCLE
             ENDIF
          ELSE
             ! logical error
             WRITE(ierrou,*)' object ', name,' not complete but not roving/sat, rec= ',n
          ENDIF
       ENDIF
       obs_found=.TRUE.
    ENDDO
    ! too many observations
    WRITE(ierrou, *)' too many MPC observations, record no ',n
    numerr=numerr+1

    ! error
11  WRITE(*,*) 'mpc_obs_input: error in reading at nobs =',nobs
    ! end of file
10  nobs=nobs-1

    ! If necessary close file and issue warning for presence of radar observations
    IF(PRESENT(filnam))THEN
       IF(rad_in_mpcfile) THEN
          lf=lench(filnam)
          !      WRITE(ierrou,120) filnam(1:lf)
          !      WRITE(ierrou,122)
          !      numerr=numerr+1
       END IF
       CALL filclo(iunit,' ')
    ELSE
       IF(rad_in_mpcfile) THEN
          !      WRITE(ierrou,121)
          !      WRITE(ierrou,122)
          !      numerr=numerr+1
       END IF
       eof=.TRUE.
    END IF
120 FORMAT('MPC file ',A,' contains radar observations')
121 FORMAT('Some input MPC file contains radar observations')
122 FORMAT(5X,'Radar observations cannot be handled correctly by MPC format:'/   &
         5X,'please use appropriate .rad file instead')

  END SUBROUTINE mpc_obs_input
  ! Transformation of a MPC record into an observation (TYPE ast_obs)
  ! This routine works on a single record and therefore cannot handle correctly satellite and radar observations
  SUBROUTINE mpcrec_transform(mpcrec,obs,complete,error,skip_rad,ons_mode)
    USE station_coordinates
    USE output_control

    CHARACTER(LEN=*),                INTENT(IN)  :: mpcrec      ! MPC record
    TYPE(ast_obs),                   INTENT(OUT) :: obs         ! astrometric observations
    LOGICAL,                         INTENT(OUT) :: error       ! error conversion flag,
    LOGICAL,                         INTENT(OUT) :: complete    ! input complete flag: if(.not.complete) the next record has to be read
    LOGICAL,                         INTENT(OUT) :: skip_rad    ! radar skip flag: if(skip_rad) a radar data has been skipped
    LOGICAL, INTENT(IN), OPTIONAL :: ons_mode ! if present, reading an ONS
    ! requiring a new designations manufactured by the perl routine
    !end interface
    CHARACTER(LEN=80) :: error_code          ! internal error code
    LOGICAL           :: err_tmp             ! error codes for called subroutines
    INTEGER           :: year,month,day
    DOUBLE PRECISION  :: seconds_utc,seconds_tdt
    CHARACTER(LEN=9)  :: chdate
    INTEGER           :: len_chdate,pos_point,n_dec_day,mjd_utc,mjd_tdt,len_error_code,len_mag_field,n_dec_mag
    CHARACTER(LEN=3)  :: timescale
    CHARACTER(LEN=5)  :: mag_field

    ! Obsolete CALLs (to be substituted by modules)
    DOUBLE PRECISION, EXTERNAL :: tjm1
    INTEGER,          EXTERNAL :: lench
    INTEGER le
    error=.true.
    error_code='undefined'

    ! Default values
    obs=undefined_ast_obs

    ! Observation type and technology
    obs%tech=mpcrec(15:15)
    ! error if it is a second record with observer position
    IF(obs%tech.eq.'v'.or.obs%tech.eq.'s')THEN
       error_code='second record with obs. pos. given as first record'
       GOTO 10
    ELSEIF(obs%tech.eq."'")THEN  ! anti-quote safety
       obs%tech=' '
    ELSEIF(obs%tech.eq.'"')THEN
       obs%tech=' '
    ELSEIF(obs%note.eq."\")THEN !" anti-backslash safety
       obs%note=' ' 
    ENDIF
    ! Radar and satellite observations are not handled at once, another record has to be read
    obs%type='O'
    skip_rad=.FALSE.
    IF(obs%tech.eq.'S') THEN
       complete=.FALSE.
       obs%type='S'
    ELSE
       complete=.TRUE.
    END IF

    IF(PRESENT(ons_mode))THEN
       IF(ons_mode)THEN
          ! Object Name: read designation created by the perl script, given by the observer, etc.
          obs%objdes=mpcrec(4:12)
          CALL rmsp(obs%objdes,le)
       ELSE
          ! Object Name
          CALL iaucod(mpcrec(1:12),obs%objdes,err_tmp)
          IF (err_tmp) THEN
             error_code='conversion of object name'
             GOTO 10
          END IF
       ENDIF
    ELSE
       ! Object Name
       CALL iaucod(mpcrec(1:12),obs%objdes,err_tmp)
       IF (err_tmp) THEN
          error_code='conversion of object name'
          GOTO 10
       END IF
    ENDIF
    ! Radar observations are not read from MPC files
    IF(obs%tech.eq.'R' .OR. obs%tech.eq.'r') THEN
       error=.FALSE.
       skip_rad=.TRUE.
       RETURN
    END IF
    ! Default value for technology descriptor is 'P' (Photographic)
    IF(obs%tech.eq.' ') obs%tech='P'
    ! Check validity of technology descriptor
    IF(ALL(tech_values.ne.obs%tech)) THEN
       WRITE(*,*) 'WARNING mpcrec_transform: technology descriptor ',obs%tech,' does not exist'
       WRITE(iun_log,*) 'WARNING mpcrec_transform: technology descriptor ',obs%tech,' does not exist'
    ENDIF

    ! Note on observation circumstances
    obs%note=mpcrec(14:14)
    IF(obs%note.eq."'")THEN ! anti-quote safety
       obs%note=' '
    ELSEIF(obs%note.eq.'"')THEN
       obs%note=' '
    ELSEIF(obs%note.eq.'\')THEN !' anti-backslash safety
       obs%note=' '
    ENDIF
    ! Observatory code
    error_code='input field: observatory code'

    READ(mpcrec(78:80),*,ERR=10) obs%obscod_s
    CALL statcode(obs%obscod_s,obs%obscod_i)
    IF(obs%obscod_i.eq.247)THEN
       complete=.false.
    ENDIF

    ! Read time fields
    error_code='input field: year, month and date'
    READ(mpcrec,100,ERR=10) year,month,chdate
100 FORMAT(15X,I4,1X,I2,1X,A9)
    len_chdate=lench(chdate)
    IF(len_chdate.LE.0) GOTO 10

    ! Time accuracy
    pos_point=INDEX(chdate,'.')
    IF(pos_point == 0) THEN
       error_code='input field: day'
       READ(chdate,*,ERR=10) day
       seconds_utc=0
       obs%acc_time=1
    ELSE
       error_code='input field: day'
       READ(chdate(1:pos_point-1),*,ERR=10) day
       error_code='input field: seconds'
       READ(chdate(pos_point:),*,ERR=10) seconds_utc
       seconds_utc=seconds_utc*86400
       n_dec_day=len_chdate-pos_point
       obs%acc_time=10.0d0**(-n_dec_day)
    END IF

    ! Conversion of time fields into MJD
    IF(year.LT.1972) THEN
       timescale='UT1'
    ELSE
       timescale='UTC'
    END IF
    mjd_utc=NINT(tjm1(day,month,year,0.d0))
    CALL cnvtim(mjd_utc,seconds_utc,timescale,mjd_tdt,seconds_tdt,'TDT')
    obs%time_utc=mjd_utc+seconds_utc/86400.d0
    obs%time_tdt=mjd_tdt+seconds_tdt/86400.d0

    ! Conversion of RA and DEC
    CALL rdanga(mpcrec(45:56),obs%coord(2),obs%acc_coord(2),err_tmp)
    IF (err_tmp) THEN
       error_code='conversion of DEC'
       GOTO 10
    ELSE
       obs%coord(2)=obs%coord(2)*radeg
       obs%acc_coord(2)=obs%acc_coord(2)*radeg
    END IF
    CALL rdanga(mpcrec(33:44),obs%coord(1),obs%acc_coord(1),err_tmp)
    IF (err_tmp) THEN
       error_code='conversion of RA'
       GOTO 10
    ELSE
       obs%coord(1)=obs%coord(1)*radh
       obs%acc_coord(1)=obs%acc_coord(1)*radh
       obs%acc_coord(1)=obs%acc_coord(1)/COS(obs%coord(2))
    END IF

    ! Magnitude
    obs%mag_str=mpcrec(66:71)
    obs%mag_band=mpcrec(71:71)
    mag_field=mpcrec(66:70) !WARNING: to fix 644 mistake
    CALL rmsp(mag_field,len_mag_field)
    mag_field=mpcrec(66:70)
    IF(len_mag_field>0) THEN
       error_code='input field: magnitude'
       READ(mag_field,*,ERR=10) obs%mag
       obs%mag_def=.true.
       pos_point=INDEX(mag_field,'.')
       IF(pos_point == 0) THEN
          obs%acc_mag=1
       ELSE
          n_dec_mag=len_mag_field-pos_point
          obs%acc_mag=10.0d0**(-n_dec_mag)
       END IF
    END IF
    ! Astrometric catalog code from MPC
    obs%catcodmpc(2:2)=mpcrec(72:72)

    ! discovery records to be removed
    !IF (obs%tech.eq.'X') THEN
    !   error_code=' corrected discovery obs. removed'
    !   GOTO 10
    !ENDIF
    ! Normal termination
    error=.false.
    RETURN

    ! Error termination
10  CONTINUE
    len_error_code=lench(error_code)
    WRITE(ierrou,101)mpcrec(1:12),error_code(1:len_error_code)
101 FORMAT(A12,'  ERROR in mpcrec_transform: ',A)
    WRITE(ierrou,'(A80)')mpcrec
  END SUBROUTINE mpcrec_transform

  ! Add the geocentric position of the observer to an astrometric observation
  ! In case of a satellite or "roving observatory" observation, the two MPC records must be passed
  SUBROUTINE mpcrec_add_obspos(obs,error,mpcrec1,mpcrec2)
    USE station_coordinates
    USE reference_systems

    TYPE(ast_obs),      INTENT(INOUT)          :: obs        ! astrometric observation to be completed
    LOGICAL,            INTENT(OUT)            :: error      ! error conversion flag
    CHARACTER(LEN=*),   INTENT(IN),   OPTIONAL :: mpcrec1    ! first  MPC record
    CHARACTER(LEN=*),   INTENT(IN),   OPTIONAL :: mpcrec2    ! second MPC record

    CHARACTER(LEN=80)                :: error_code             ! internal error code
    CHARACTER*3                      :: obscods
    INTEGER                          :: len_error_code,sign,obscod,le
    DOUBLE PRECISION                 :: longitude,latitude,altitude,conv
    DOUBLE PRECISION, DIMENSION(3)   :: bfpos0,equpos
    ! DOUBLE PRECISION, DIMENSION(3,3) :: rot
    CHARACTER*5 :: altstr

    ! astronomical unit au in km from fund_const

    ! Obsolete CALLs (to be substituted by modules)
    INTEGER, EXTERNAL :: lench
    Character*9 name ! for message

    error=.false.
    name=obs%objdes

    SELECT CASE (obs%type)
    CASE ('O')
       IF(obs%obscod_s == '247') THEN
          ! Roving observatory
          IF(.NOT.PRESENT(mpcrec1)) STOP '**** mpcrec_add_obspos: internal error (01) ****'
          IF(.NOT.PRESENT(mpcrec2)) STOP '**** mpcrec_add_obspos: internal error (02) ****'
          IF(mpcrec1(1:14) /= mpcrec2(1:14)) THEN
             error_code='mismatch in col 1-14'
             GOTO 10
          END IF
          IF(mpcrec2(15:15) /= 'v') THEN
             error_code='missing code "v" in column 15'
             GOTO 10
          END IF
          IF(mpcrec1(16:32) /= mpcrec2(16:32)) THEN
             error_code='mismatch in col 16-32'
             GOTO 10
          END IF
          IF(mpcrec2(33:34) /= '1 ') THEN
             error_code='abnormal content of columns 33-34'
             GOTO 10
          END IF
          IF(mpcrec2(78:80) /= '247') THEN
             error_code='abnormal content of columns 78-80'
             GOTO 10
          END IF
          ! Input and transformation of coordinates
          error_code='input field: longitude'
          READ(mpcrec2(35:44),*,ERR=10) longitude
          SELECT CASE (mpcrec2(46:46))
          CASE (' ', '+')
             sign=+1
          CASE ('-')
             sign=-1
          CASE DEFAULT
             error_code='input field: sign of latitude'
             GOTO 10
          END SELECT
          error_code='input field: latitude'
          READ(mpcrec2(47:55),*,ERR=10) latitude
          latitude=sign*latitude
          error_code='input field: altitude'
          altstr=mpcrec2(57:61)
          CALL rmsp(altstr,le)
          IF(le.gt.0)THEN
             READ(mpcrec2(57:61),*,ERR=10) altitude
          ELSE
             ! altitude error case; does not matter
             WRITE(ierrou,102)mpcrec2(57:61),obs%objdes
102          FORMAT('mpcrec_add_obspos: altitude is ',a5,' for ',a9)
             numerr=numerr+1
             altitude=0.d0
          ENDIF
          longitude=longitude*radeg
          latitude=latitude*radeg
          ! Store information on the geodetic position of the observatory (needed only for output)
          IF(ASSOCIATED(obs%geopos)) DEALLOCATE(obs%geopos)
          ALLOCATE(obs%geopos(3))
          obs%geopos(1)=longitude
          obs%geopos(2)=latitude
          obs%geopos(3)=altitude
          CALL geodetic_to_cartesian(longitude,latitude,altitude,bfpos0)
          bfpos0=bfpos0/(aukm*1d3)
          CALL observer_position(obs%time_tdt,obs%obspos,obs%obsvel,BFPOS=bfpos0)
       ELSE
           ! Fixed astronomical observatory
          CALL observer_position(obs%time_tdt,obs%obspos,obs%obsvel,OBSCODE=obs%obscod_i)
       END IF

    CASE ('S')
       ! Satellite observation
       IF(.NOT.PRESENT(mpcrec1)) STOP '**** mpcrec_add_obspos: internal error (03) ****'
       IF(.NOT.PRESENT(mpcrec2)) STOP '**** mpcrec_add_obspos: internal error (04) ****'
       IF(mpcrec1(1:12) /= mpcrec2(1:12)) THEN
          error_code='mismatch in col 1-12 '//name
          GOTO 10
       END IF
       IF(mpcrec2(15:15) /= 's') THEN
          error_code='missing code "s" in column 15 '//name
          GOTO 10
       END IF
       IF(mpcrec1(16:32) /= mpcrec2(16:32)) THEN
          error_code='mismatch in col 16-32 '//name
          GOTO 10
       END IF
       ! Conversion of unit of length
       SELECT CASE (mpcrec2(33:33))
       CASE ('1')
          conv=1.d0/aukm
       CASE ('2')
          conv=1.d0
       CASE DEFAULT
          error_code='input field: unit of parallax indicator '//name
          GOTO 10
       END SELECT
       obs%par_unit=mpcrec2(33:33)
       error_code='input field: parallax vector'
       READ(mpcrec2(36:45),*,ERR=10) equpos(1)
       IF(mpcrec2(35:35) == '-') equpos(1)=-equpos(1)
       READ(mpcrec2(48:57),*,ERR=10) equpos(2)
       IF(mpcrec2(47:47) == '-') equpos(2)=-equpos(2)
       READ(mpcrec2(60:69),*,ERR=10) equpos(3)
       IF(mpcrec2(59:59) == '-') equpos(3)=-equpos(3)
       ! Transformation of equatorial coordinates into ecliptical
       !      CALL rotpn(rot,'EQUM','J2000',0.d0,'ECLM','J2000',0.d0)
       obs%obspos=MATMUL(roteqec,equpos)
       obs%obspos=conv*obs%obspos
       obs%obsvel=0.d0
       error_code='input field: observatory code '//name
       READ(mpcrec1(78:80),*,ERR=10) obscods
       IF(obscods /= obs%obscod_s) THEN
          error_code='abnormal content of columns 78-80 '//name
          GOTO 10
       END IF

    CASE ('R', 'V')
       ! Radar observations
       WRITE(*,*)' mpcrec_add_obspos: MPC radar not alllowed '
       STOP
       CALL observer_position(obs%time_tdt,obs%obspos,obs%obsvel,OBSCODE=obs%obscod_i,PRECISION=2)
    CASE DEFAULT
       ! Unknown type of observation
       error_code='observation type: ' // obs%type
       GOTO 10

    END SELECT

    RETURN

10  CONTINUE
    len_error_code=lench(error_code)
    WRITE(ierrou,101) error_code(1:len_error_code)
101 FORMAT('*** ERROR in mpcrec_add_obspos: ',A)
    numerr=numerr+1
    error=.true.
  END SUBROUTINE mpcrec_add_obspos
  ! END SUBMODULE iompc

  ! SUBMODULE ioradar
  SUBROUTINE  jpl_radarobs(obs_found,file,obs,nobs,frq)

    CHARACTER*(*), INTENT(IN) :: file
    LOGICAL, INTENT(OUT) :: obs_found ! true if there is something
    INTEGER, INTENT(OUT) :: nobs ! number of observations
    TYPE(ast_obs), INTENT(OUT), DIMENSION(:) :: obs ! vector of observations
    DOUBLE PRECISION, INTENT(OUT), DIMENSION(:), OPTIONAL :: frq
    ! end interface
    LOGICAL error
    CHARACTER*100 rec
    INTEGER i,l,unit,nlef
    external lench
    integer lench
    nlef=SIZE(obs)
    call filopn(unit,file,'old')
    nobs=0
    DO  i=1,nlef
       read(unit,101,end=10) rec
101    format(a)
       l=lench(rec)
       IF(l.eq.0) CYCLE
       nobs=nobs+1
       radar_rec(i)=rec
       ! parse one jpl radar record
       CALL jplradar_transform(rec,obs(i),error)
       IF(error)THEN
          l=lench(file)
          write(*,102) file(1:l),nobs
102       format(' **** jpl_radarobs: input conversion error ****'/          &
               &        ' **** file = "',a,'", line',i5,' ****')
          nobs=0
          obs_found=.false.
          call filclo(unit,' ')
          RETURN
       ENDIF
       IF(PRESENT(frq)) READ(rec,103) frq(i)
103    FORMAT(73X,F5.0)
    ENDDO
    ! too many records
    STOP ' **** jpl_radarobs: nobs > nlef ****'

    ! regular ending
10  continue
    call filclo(unit,' ')
    IF(nobs.eq.0)THEN
       obs_found=.false.
       write(*,*) 'Warning: Found ',nobs,' obs in ',file
    ELSE
       obs_found=.true.
       IF(verb_io.gt.9)write(*,*) 'Found ',nobs,' obs in ',file
    ENDIF

    RETURN
  END SUBROUTINE jpl_radarobs

  ! Revised by Genny on 5 November 2005 to deal with the new format of 
  ! radar_data_ast.txt to include numbered objects > 99999
  SUBROUTINE jplradar_transform(rec,obs,error)
    USE station_coordinates

    CHARACTER*(*), INTENT(IN)  :: rec
    TYPE(ast_obs), INTENT(OUT) :: obs
    LOGICAL,       INTENT(OUT) :: error

    INTEGER errcod
    ! name
    character*6 number
    character*17 nametmp,stname
    integer lnum,lnam
    ! time
    integer year,month,day,hour,min,isec,mjd,mjdt
    double precision sec,sect
    character*3 scale
    ! funcs
    INTEGER ix
    DOUBLE PRECISION tjm1
    EXTERNAL tjm1
    ! measurements
    double precision observ,rms,freq
    DOUBLE PRECISION r,v, accr,accv
    character*2 unit
    character*3 surf
    logical range
    ! obscode
    character*(name_len) trxstr,recstr
    integer iotr,iore
    CHARACTER(LEN=3) :: code_trx,code_rcx

    error=.true.
    errcod=1
    obs=undefined_ast_obs
    ! ========== HANDLE OBJECT NAME =================
    READ(rec,101,ERR=10) number,nametmp
101 FORMAT(a6,1x,a16)
    call rmsp(number,lnum)
    call rmsp(nametmp,lnam)
    if(lnam.gt.9)lnam=9
    if(lnum.eq.0)then
       ! remove identified object if present ('1991AQ=1994RD')
       ix=index(nametmp,'=')
       if(ix.gt.0)lnam=ix-1
       obs%objdes=nametmp(1:lnam)
    else
       obs%objdes=number(1:lnum)
    endif
    ! ========== HANDLE DATE AND TIME =================
    READ(rec,102,ERR=10) year,month,day,hour,min,isec
102 FORMAT(24X,I4,5(1X,I2))
    IF(year.LT.1972) THEN
       scale='UT1'
    ELSE
       scale='UTC'
    END IF
    sec=(hour*60d0+min)*60d0+isec
    mjd=nint(tjm1(day,month,year,0.d0))
    CALL cnvtim(mjd,sec,scale,mjdt,sect,'TDT')
    obs%time_utc=mjd+sec/86400.d0
    obs%time_tdt=mjdt+sect/86400.d0
    obs%acc_time=10d-10
    ! ========== READ MEASUREMENT =================
    READ(rec,103,ERR=10) observ,rms,unit,surf,freq
103 FORMAT(44x,f13.2,1x,f7.3,1x,a2,1x,a3,1x,f5.0)
    freq=freq*1d6
    if(unit.eq.'us')then
       obs%type='R'
       range=.true.
    elseif(unit.eq.'Hz')then
       obs%type='V'
       range=.false.
    else
       errcod=37
       WRITE(*,*)rec
       goto 10
    endif
    if(surf.eq.'COM')then
       obs%tech='c'
    elseif(surf.eq.'PP ')then
       obs%tech='s'
    else
       errcod=47
       goto 10
    endif
    if(range)then
       r=observ*1d-6/86400d0*vlight/2.d0
       accr=rms*1d-6/86400d0*vlight/2.d0
       obs%coord(1)=r
       obs%acc_coord(1)=accr
       obs%coord(2)=0.d0
       obs%acc_coord(2)=-1.d0
    else
       v=-observ*vlight/(freq*2.d0)
       accv=rms*vlight/(freq*2.d0)
       obs%coord(2)=v
       obs%acc_coord(2)=accv
       obs%coord(1)=0.d0
       obs%acc_coord(1)=-1.d0
    endif
    ! ======================= HANDLE OBSERVATORY CODES ====================
    READ(rec,104,ERR=10) trxstr,recstr
104 FORMAT(79x,a9,1x,a9)
    errcod=5
    IF(trxstr.eq.'-73')THEN
       trxstr='EISCAT'
       recstr='EISCAT'
    END IF
    iotr=radar_station(trxstr)
    iore=radar_station (recstr)
    obs%obscod_i=iotr*10000+iore
    CALL codestat(iotr,code_trx)
    CALL codestat(iore,code_rcx)
    obs%obscod_s=code_trx//' '//code_rcx
    ! obs%obspos = body-fixed position of transmitting station
    ! obs%obsvel = body-fixed position of receiving station
    CALL obscoo(iotr,obs%obspos,stname)
    CALL obscoo(iore,obs%obsvel,stname)
    ! ======================= FILL DUMMY FIELDS ===========================
    obs%note=' '
    obs%mag=0.d0
    obs%mag_str='      '
    obs%mag_band=' '
    obs%acc_mag=0.d0
    obs%mag_def=.FALSE.
    ! ====================== CLEAN UP =====================================
    error=.false.
10  CONTINUE
    IF(error) THEN
       WRITE(*,105) errcod,rec
105    FORMAT(' jplrad: error code',I3,'rec:',A)
    ENDIF
    RETURN

  END SUBROUTINE jplradar_transform

  ! ====================== station func =============================
  ! compute observatory code from string
  INTEGER FUNCTION radar_station(stastr)
    character*(*) stastr
    ! We have the following stations:
    !     'Arecibo  ' = 251
    !     'DSS 13   ' = 252
    !     'DSS 14   ' = 253
    !     'DSS 25   ' = 257
    !     'Haystack ' = 254
    !     'Evpatoria' = 255 (not in MPC file...)
    !     'Greenbank' = 256
    !     'EISCAT'    = 273 (not in MPC file...)
    if(stastr.eq.'Arecibo  ')then
       radar_station=251
    elseif(stastr.eq.'DSS 13   ')then
       radar_station=252
    elseif(stastr.eq.'DSS 14   ')then
       radar_station=253
    elseif(stastr.eq.'DSS 25   ')then
       radar_station=257
    elseif(stastr.eq.'Haystack ')then
       radar_station=254
    elseif(stastr.eq.'Evpatoria')then
       radar_station=255
    elseif(stastr.eq.'Greenbank')then
       radar_station=256
    elseif(stastr.eq.'EISCAT   ')then
       !     radar_station=273
       radar_station=259
    else
       radar_station=500
       write(*,*)'jplrad: unknown station: ',stastr
       STOP '***unknown radar station****'
    endif
  END FUNCTION radar_station
  ! END SUBMODULE ioradar
  ! Input of observations: reads sequentially the .rwo, .obs, .rad file and combines
  ! the data according to the precedence rule specified by precob

  SUBROUTINE input_obs(obsdir,astna0,precob,error_model,obs0,obs,obsw,m,iun20,change,rms,rmsmag,ades_type)

    CHARACTER(LEN=*),             INTENT(IN)            :: obsdir      ! input dir. where files astna0.rwo, astna0.psv/xml, astna0.obs and astna0.rad are
    CHARACTER(LEN=*),             INTENT(IN)            :: astna0      ! asteroid name
    INTEGER,                      INTENT(IN)            :: iun20       ! messages unit
    LOGICAL,                      INTENT(IN)            :: precob      ! .TRUE. for overwriting .rwo, .FALSE. for updating .rwo
    CHARACTER(LEN=*),             INTENT(IN)            :: error_model ! error model file name
    INTEGER,                      INTENT(IN), OPTIONAL  :: ades_type   ! flag for type of ADES output file (.psv/.xml)
    LOGICAL,                      INTENT(OUT)           :: obs0        ! successful input
    LOGICAL,                      INTENT(OUT)           :: change      ! existence of new data
    INTEGER,                      INTENT(OUT)           :: m           ! number of observations
    TYPE(ast_obs),  DIMENSION(:), INTENT(OUT)           :: obs         ! observations
    TYPE(ast_wbsr), DIMENSION(:), INTENT(OUT)           :: obsw        ! weights and possibly residuals
    DOUBLE PRECISION,             INTENT(OUT), OPTIONAL :: rms         ! RMS of astrometric fit (-1 = undefined)
    DOUBLE PRECISION,             INTENT(OUT), OPTIONAL :: rmsmag      ! RMS of photometric fit (-1 = undefined)

    ! observation numbers: maximum, space left
    INCLUDE 'parobx.h90'
    INTEGER :: nlef
    INTEGER :: ades_type_loc
    ! file names
    CHARACTER*77 file
    CHARACTER*8 ext
    INTEGER lfile
    LOGICAL rwo,mpc,rad,mpc_psv,rwo_psv  ! flags for existence of input files
    LOGICAL ades                 ! flag for existence of ADES (either psv or xml) input file 
    ! new obs. number
    INTEGER mnew,mr,nlefm,mnewt,mj,double
    ! ===== observational data: temporary copy===========
    INTEGER mt ! observation number
    TYPE(ast_obs),         DIMENSION(nobx) :: obst  ! observations
    TYPE(ast_wbsr),        DIMENSION(nobx) :: obswt ! weights/residuals
    
    CHARACTER*(20) :: error_modelt ! error model file name
    LOGICAL        :: init ! for observ_rms: if not initialized, obsw should be set to default

    ! observational ADES data: temporary copy
    CHARACTER(LEN=200) :: ades_filename_out,ades_filename_tmp
    INTEGER            :: nobs_ades_out,nobs_ades_tmp
    INTEGER            :: ntot_ades_out,ntot_ades_tmp
    INTEGER            :: header_length_out,header_length_tmp
    CHARACTER(LEN=10)  :: ades_error_model_out,ades_error_model_tmp 
    TYPE(ades_observ),    DIMENSION(nobx) :: ades_obs_out,ades_obs_tmp
    TYPE(ades_flags),     DIMENSION(nobx) :: ades_flg_out,ades_flg_tmp
    TYPE(ades_precision), DIMENSION(nobx) :: ades_prec_out,ades_prec_tmp
    TYPE(ades_fit_var),   DIMENSION(nobx) :: ades_fit_out,ades_fit_tmp
     CHARACTER(LEN=200),  DIMENSION(nobx) :: ades_header_out,ades_header_tmp
    CHARACTER(LEN=600),   DIMENSION(nobx) :: ades_key_rec_out,ades_key_rec_tmp
    CHARACTER(LEN=600),   DIMENSION(nobx) :: ades_data_rec_out,ades_data_rec_tmp
    INTEGER,              DIMENSION(nobx) :: ntot_to_nobs_out,ntot_to_nobs_tmp
    INTEGER,              DIMENSION(nobx) :: nobs_to_ntot_out,nobs_to_ntot_tmp
    LOGICAL,              DIMENSION(nobx) :: new_key_rec_out,new_key_rec_tmp
    ! ===============================================
    ! sorting
    INTEGER iperm(nobx)
    ! directory char
    INCLUDE 'sysdep.h90'
    ! loop index
    INTEGER i,j
    INTEGER ld
    INTEGER, EXTERNAL :: lench
    LOGICAL changemod       ! error model changed with respect to previous .rwo file
    LOGICAL change_ades     ! true if ADES data has been changed
    LOGICAL change_ades_errmod    ! true if error_model has changed wrt previous ADES file
    CHARACTER(LEN=1), DIMENSION(nobx) :: disc_flag,disct,disc_flag_loc  ! discovery flags (from/for ADES file, from .obs)
    REAL(KIND=dkind), DIMENSION(nobx) :: freq,frqt,frq_loc              ! carrier reference frequency (for ADES file, from .rad)
    CHARACTER(LEN=4), DIMENSION(nobx) :: subFmt,subFmt_loc              ! Format in which the observation was originally submitted to the MPC
    ! =============EXECUTION BEGINS======================
    IF(PRESENT(rms)) rms=-1.d0
    IF(PRESENT(rmsmag)) rmsmag=-1.d0
    IF(.NOT.PRESENT(ades_type)) THEN
       ades_type_loc = 0   ! if flag is not present, we set it to print only .rwo file
    ELSE
       ades_type_loc = ades_type
    ENDIF
    
    changemod = .FALSE.
    new_ades = .TRUE.
    change_ades = .FALSE.
    change_ades_errmod = .FALSE.
    disc_flag = ""
    disct = ""
    freq = 0.d0
    frqt = 0.d0
    subFmt = ""
    subFmt_loc = ""

    !  compute file name
    ld=lench(obsdir)
    IF(ld.GT.0) THEN
       IF(obsdir(ld:ld).EQ.dircha) THEN
          file=obsdir(1:ld)//astna0
       ELSE
          file=obsdir(1:ld)//dircha//astna0
       END IF
    ELSE
       file=astna0
    END IF
    CALL rmsp(file,lfile)

    ! initialize ades file name 
    psv_name = file(1:lfile)
    CALL rmsp(psv_name,len_ades_name)

    ! WARNING: cases dealing with ADES data from/to .xml file are provisionally neglected
    ! thus, the corresponding parts of code are commented
    
    ! existence of .rwo, .obs, .rad, .obs_psv, .rwo_psv
    INQUIRE(file=file(1:lfile)//'.rwo',exist=rwo)
    INQUIRE(file=file(1:lfile)//'.obs',exist=mpc)
    INQUIRE(file=file(1:lfile)//'.rad',exist=rad)
    INQUIRE(file=file(1:lfile)//'.rwo_psv',exist=rwo_psv)
    INQUIRE(file=file(1:lfile)//'.obs_psv',exist=mpc_psv)
    !INQUIRE(file=file(1:lfile)//'.psv',exist=psv)
    !INQUIRE(file=file(1:lfile)//'.xml',exist=xml) 
    !ades = psv .OR. xml
    ades = rwo_psv .OR. mpc_psv

    ! check whether there is not coexistence of .obs and .obs_psv
    ! in any case, give precedence to ADES
    IF (mpc .AND. mpc_psv) THEN
       WRITE(*,*) ' input_obs: ',file(1:lfile)//'.obs and ',file(1:lfile)//'.obs_psv should not coexist.'
       WRITE(iun20,*) ' input_obs: ',file(1:lfile)//'.obs and ',file(1:lfile)//'.obs_psv should not coexist.'
    ENDIF

    ! control about consistence of precob and the existence of obs_psv
    IF (precob .AND. .NOT. mpc_psv) ades = .FALSE.  ! use data from .obs/.rad + .rwo files, if they exist; otherwise, return

    ! check sizes of arrays
    nlef=SIZE(obs)
    IF(nlef.ne.SIZE(obsw))THEN
       WRITE(*,*)' input_obs: problems in dimensions of arrays obs, obsw ',nlef,SIZE(obsw)
       STOP
    ENDIF
    m=0
    ! select operations mode
    IF(.NOT.rwo.AND..NOT.ades)THEN
       ! For now we will not do radar only orbits
       ! optical observations file is required
       IF(.NOT.mpc)THEN
          !WRITE(*,*)'You must provide at least one of .obs, .rwo, .psv and .xml file in directory ',obsdir
          WRITE(*,*)'You must provide at least one of .obs, .rwo or .obs_psv, .rwo_psv files in directory ',obsdir
          obs0=.false.
          RETURN
       ENDIF
       IF(verb_io.gt.9)WRITE(*,*) 'Neither .rwo nor ADES files, reading .obs and/or .rad files.'
       ! there is neither .rwo nor ades, so read .obs and .rad
       IF(mpc)THEN
          ! Input of astrometric observations from a file (MPC format)
          CALL mpc_obs_input(mpc,obs,m,FILNAM=file(1:lfile)//'.obs',ONS_MODE=ons_name,DISC=disc_flag_loc)
          disc_flag(1:m) = disc_flag_loc(1:m)
          subFmt(1:m) = "M92" 
          IF(verb_io.gt.9)WRITE(*,*)'input_obs: ',m,' obs in ',file(1:lfile)//'.obs'
          WRITE(iun20,*)'input_obs: ',m,' obs in ',file(1:lfile)//'.obs'
       ENDIF
       ! read radar jpl data
       IF(rad)THEN
          nlefm=nlef-m
          CALL jpl_radarobs(rad,file(1:lfile)//'.rad',obs(m+1:),mr,FRQ=frq_loc)
          IF(verb_io.gt.9)WRITE(*,*)'input_obs:',mr,' obs from  ',file(1:lfile)//'.rad'
          WRITE(iun20,*)'input_obs:',mr,' obs  ',file(1:lfile)//'.rad'
          ! adjust values positions in frequency array
          freq(m+1:m+mr) = frq_loc(1:mr)
          m=mr+m
       ENDIF
       obs0=mpc.or.rad
       IF(.not.obs0)RETURN
       ! find weights for these; a priori RMS of astrometric observations
       init=.true.
       CALL observ_rms(obs,error_model,init,obsw,m)
       ! output data for possible manual fixing: create weights file; no rms available
       IF (ades_type_loc .GE. 0) THEN
          CALL write_rwo(file(1:lfile)//'.rwo',obs,obsw,m,error_model)
          change = .TRUE.
       ENDIF
       IF (ades_type_loc .NE. 0) THEN
          ! convert data to ADES format
          ades_filename = psv_name(1:len_ades_name)//'.rwo_psv'
          nobs_ades = m
          CALL convert_fitobs_ades(obs,obsw,ades,.FALSE.,disc_flag,freq,subFmt)
          ! create output ADES file
          CALL write_ades(psv_name(1:len_ades_name),error_model,ades_type_loc)
          change_ades = .TRUE.
       ENDIF

    ELSEIF (.NOT.ades) THEN
       IF(verb_io.gt.9)WRITE(*,*) file(1:lfile),' ADES files not found, ALL obs will come from either .obs/.rad ', & 
            & 'files or .rwo file.'
       ! select between update and overwrite of .rwo using precob flag: .psv file is written using data from either
       ! rwo file in the former case or obs/rad files in the latter case  
       IF(precob)THEN
          IF(verb_io.gt.9)WRITE(*,*) file(1:lfile),'.rwo found but ALL obs will come from .obs/.rad files.'
          ! give the precedence to the observation files .obs and .rad
          ! with respect to .rwo, which is overwritten
          IF(mpc)THEN
             ! Input of astrometric observations from a file (MPC format)
             CALL mpc_obs_input(mpc,obs,m,FILNAM=file(1:lfile)//'.obs',ONS_MODE=ons_name,DISC=disc_flag_loc)
             disc_flag(1:m) = disc_flag_loc(1:m)
             subFmt(1:m) = "M92"
             IF(verb_io.gt.9)WRITE(*,*)'input_obs: ',m,' obs in ',file(1:lfile)//'.obs'
             WRITE(iun20,*)'input_obs: ',m,' obs in ',file(1:lfile)//'.obs'
          ELSE
             m=0
          ENDIF
          ! read radar jpl data
          IF(rad)THEN
             nlefm=nlef-m
             CALL jpl_radarobs(rad,file(1:lfile)//'.rad',obs(m+1:),mr,FRQ=frq_loc)
             IF(verb_io.gt.9)WRITE(*,*)'input_obs:',mr,' obs from  ',file(1:lfile)//'.rad'
             WRITE(iun20,*)'input_obs:',mr,' obs  ',file(1:lfile)//'.rad'
             ! adjust values positions in frequency array
             freq(m+1:m+mr) = frq_loc(1:mr)
             m=mr+m
          ENDIF
          obs0=rad.or.mpc
          ! If no obs then object does not "exist", so rwo should not be read:
          IF(.not.obs0)return
          ! find weights for these; a priori RMS of astrometric observations
          init=.true.
          CALL observ_rms(obs,error_model,init,obsw,m)
          ! give default selection flag of 1
          obsw(1:m)%sel_coord=1
          ! format transition:
          !  CALL read_rwo_trans(file(1:lfile)//'.rwo',obst,obswt,mt,error_modelt,rms,rmsmag)
          ! read .rwo  anyway, but store in temporary array the data
          CALL read_rwo(file(1:lfile)//'.rwo',obst,obswt,mt,error_modelt,rms,rmsmag)
          IF(error_model.ne.error_modelt)THEN
             init=.false. ! structure is initialized already
             CALL observ_rms(obst,error_model,init,obswt,mt)
             changemod=.true.
          ELSE
             changemod=.false.
          ENDIF
          ! recover informations from .rwo (only selection flags)
          CALL addobs_rwo(obs,obsw,m,obst,obswt,mt,change)
          ! if error model is changed, the data have to be considered changed
          change=(change.or.changemod) 
          IF(change)THEN
             IF(verb_io.gt.9)WRITE(*,*)'There are new/changed obs. New numobs=',m
             ! output updated .rwo file, if there are new observations (erasing residuals)
             IF (ades_type .GE. 0) CALL write_rwo(file(1:lfile)//'.rwo',obs,obsw,m,error_model)
          ELSE
             IF(verb_io.gt.9)WRITE(*,*)'There are no updates in .obs or .rad files.'
          ENDIF
       ELSE
          ! give the precedence to .rwo, the .obs and .rad files are intended
          ! as additional observations only; data in .rwo are not erased, can only be
          ! changed
          ! read .rwo, and store data in final array
          IF(verb_io.gt.9)WRITE(*,*)'Using .rwo file, but checking .obs,.rad for update.'
          CALL read_rwo(file(1:lfile)//'.rwo',obs,obsw,m,error_modelt,rms,rmsmag)
          IF(verb_io.gt.9)WRITE(*,*)'input_obs: ',m,' obs from  ',file(1:lfile)//'.rwo'
          WRITE(iun20,*)'input_obs: ',m,' obs from ',file(1:lfile)//'.rwo'
          IF(m.eq.0)THEN
             obs0=.false.
          ELSE
             obs0=.true.
          ENDIF
          IF(error_model.ne.error_modelt)THEN
             init=.false. ! structure is initialized already
             CALL observ_rms(obs,error_model,init,obsw,m)
             changemod=.true.
          ELSE
             changemod=.false.
          ENDIF
          ! if there are input data
          IF(mpc)THEN
             ! Input of astrometric observations into temporary from a file (MPC form)
             CALL mpc_obs_input(mpc,obst,mt,FILNAM=file(1:lfile)//'.obs',ONS_MODE=ons_name,DISC=disc_flag_loc)
             disct(1:mt) = disc_flag_loc(1:mt)
             subFmt_loc(1:mt) = "M92"
             !      IF(error_model.EQ.'cbm10') THEN
             !         obs%catcodmpc=obst%catcodmpc
             !      ENDIF
             IF(.not.mpc)THEN
                WRITE(*,*) file(1:lfile)//'.obs is possibly corrupt. Not using any data from this file.'
             ELSE
                IF(verb_io.gt.9)WRITE(*,*)'input_obs:',mt,' obs from  ',file(1:lfile)//'.obs'
                WRITE(iun20,*)'input_obs:',mt,' from ',file(1:lfile)//'.obs'
             ENDIF
          ENDIF
          ! read radar jpl data
          IF(rad)THEN
             nlefm=nobx-mt
             CALL jpl_radarobs(rad,file(1:lfile)//'.rad',obst(mt+1:),mr,FRQ=frq_loc)
             IF(.not.rad)THEN
                WRITE(*,*) file(1:lfile)//'.rad is possibly corrupt. Not using any data from this file.'
             ELSE
                IF(verb_io.gt.9)WRITE(*,*)'input_obs:',mr,' obs from  ',file(1:lfile)//'.rad'
                WRITE(iun20,*)'input_obs:',mr,' obs  ',file(1:lfile)//'.rad'
                ! adjust values positions in frequency array
                frqt(mt+1:mt+mr) = frq_loc(1:mr)
                mt=mr+mt
             ENDIF
          ENDIF
          ! add information from .obs and .rad
          IF(mpc.or.rad)THEN
             obs0=.true.
             ! find weights for these; a priori RMS of astrometric observations
             init=.true.
             CALL observ_rms(obst,error_model,init,obswt,mt)
             CALL addobs_mpc(obs,obsw,m,obst,obswt,mt,mnew,change)
             ! re-assign positions in frequency and disc_flag array according to obs order given by .rwo and addobs_mpc
             DO j=1,mt
                mj=find_obs(obst(j),obs,mnew,double)
                IF (mj .GT. 0) THEN 
                   IF(rad) freq(mj) = frqt(j)
                   disc_flag(mj) = disct(j)
                   subFmt(mj) = subFmt_loc(j)
                ENDIF
             ENDDO
             m=mnew
          ELSE
             change=.true.
          ENDIF
          IF (.NOT.obs0) RETURN
          ! if error model is changed, the data have to be considered changed
          change=(change.or.changemod)   
          IF(change)THEN
             IF(verb_io.gt.9)WRITE(*,*)'There are new/changed obs. New numobs=',m
             ! output updated .rwo file, if there are new observations (erasing residuals)
             IF (ades_type_loc .GE. 0) CALL write_rwo(file(1:lfile)//'.rwo',obs,obsw,m,error_model)
          ELSE
             IF(verb_io.gt.9)WRITE(*,*)'There are no updates in .obs or .rad files.'
          ENDIF
          !change=.true.
       ENDIF
       IF (ades_type_loc .NE. 0) THEN
          ! write output ADES file
          ades_filename = psv_name(1:len_ades_name)//'.rwo_psv'
          nobs_ades = m
          CALL convert_fitobs_ades(obs,obsw,ades,.FALSE.,disc_flag,freq,subFmt)
          CALL write_ades(psv_name(1:len_ades_name),error_model,ades_type_loc)
          change_ades = .TRUE.
       ENDIF

    ELSEIF (ades) THEN
       ! neglect .xml input file
       !IF (xml) THEN   
       !   IF(verb_io.gt.9) WRITE(*,*) file(1:lfile),'.xml file found. Using .xml file.'
       !   CALL read_ades_xml(psv_name(1:len_ades_name)//'.xml')
       !   ext = '.xml'
       !ELSE
       ! initialization of public ADES variables
       nobs_ades = 0
       ntot_ades = 0
       ades_obs = undefined_ades_observ
       ades_flg = undefined_ades_flags
       ades_prec = undefined_ades_precision
       ades_fit = undefined_ades_fit_var
       ades_header = ""
       ades_key_rec = ""
       ades_data_rec = ""
       header_length = 0
       ntot_to_nobs_idx = 0
       nobs_to_ntot_idx = 0
       ades_error_model = ""
       new_key_rec = .FALSE.
       ! read .obs_psv
       IF(.NOT.rwo_psv)THEN
          IF(verb_io.gt.9)WRITE(*,*) 'No .rwo_psv file, reading .obs_psv file.'
          ! there is no .rwo_psv, so read .obs_psv
          IF(verb_io.gt.9) WRITE(*,*) file(1:lfile),'.obs_psv file found. Using .obs_psv file'
          ext = '.obs_psv'
          ades_filename_out = psv_name(1:len_ades_name)//ext
          CALL read_ades_psv(ades_filename_out,ades_obs_out,ades_flg_out,ades_prec_out, &
               & ades_fit_out,nobs_ades_out,ntot_ades_out,ades_header_out,ades_key_rec_out,ades_data_rec_out, &
               & ntot_to_nobs_out,nobs_to_ntot_out,header_length_out,ades_error_model_out,new_key_rec_out)
          m = nobs_ades_out
          IF(m .EQ. 0) THEN
             obs0 = .FALSE.
          ELSE
             obs0 = .TRUE.
          ENDIF
          IF (.NOT. obs0) RETURN
          IF(verb_io.gt.9)WRITE(*,*)'input_obs: ',m,' obs in ',trim(adjustl(ades_filename_out)) 
          WRITE(iun20,*)'input_obs: ',m,' obs in ',trim(adjustl(ades_filename_out))
          change_ades = .TRUE.
          new_ades = .FALSE.
       ELSE
          ! select between update and overwrite of .rwo_psv using precob flag: .rwo_psv file is written using data from either
          ! rwo_psv file in the former case or .obs_psv file in the latter case  
          IF(precob)THEN
             IF(verb_io.gt.9)WRITE(*,*) file(1:lfile),'.rwo_psv found but ALL obs will come from .obs_psv file.'
             ! give the precedence to the observation file .obs_psv
             ! with respect to .rwo_psv, which is overwritten
             IF(mpc_psv)THEN
                ! Input of astrometric observations from a file (psv from MPC database)
                ext = '.obs_psv'
                ades_filename_out = psv_name(1:len_ades_name)//ext
                CALL read_ades_psv(ades_filename_out,ades_obs_out,ades_flg_out,ades_prec_out,ades_fit_out, &
                     & nobs_ades_out,ntot_ades_out, ades_header_out,ades_key_rec_out,ades_data_rec_out, &
                     & ntot_to_nobs_out,nobs_to_ntot_out,header_length_out,ades_error_model_out,new_key_rec_out)
                m = nobs_ades_out
             ENDIF
             IF(m .EQ. 0) THEN
                obs0 = .FALSE.
             ELSE
                obs0 = .TRUE.
             ENDIF
             IF (.NOT. obs0) RETURN
             IF(verb_io.gt.9)WRITE(*,*)'input_obs: ',m,' obs in ',trim(adjustl(ades_filename_out)) 
             WRITE(iun20,*)'input_obs: ',m,' obs in ',trim(adjustl(ades_filename_out))
             ! find weights for these; a priori RMS of astrometric observations
             !CALL convert_ades_fitobs(error_model,obs,obsw)
             ! read .rwo_psv  anyway, but store in temporary array the data
             IF(verb_io.gt.9)WRITE(*,*)'Recovering information from .rwo_psv'
             ades_filename_tmp = psv_name(1:len_ades_name)//'.rwo_psv'
             CALL read_ades_psv(ades_filename_tmp,ades_obs_tmp,ades_flg_tmp,ades_prec_tmp, &
               & ades_fit_tmp,nobs_ades_tmp,ntot_ades_tmp,ades_header_tmp,ades_key_rec_tmp,ades_data_rec_tmp, &
               & ntot_to_nobs_tmp,nobs_to_ntot_tmp,header_length_tmp,ades_error_model_tmp,new_key_rec_tmp)
             mt = nobs_ades_tmp
             IF(verb_io.gt.9)WRITE(*,*)'input_obs:',mt,' obs in ',trim(adjustl(ades_filename_tmp))
             WRITE(iun20,*)'input_obs:',mt,' in ',trim(adjustl(ades_filename_tmp))
             !CALL convert_ades_fitobs(error_model,obst,obswt)
             ! recover information from .rwo_psv (forced weights and selection flags)
             CALL update_ades_rwo(ades_obs_out,ades_flg_out,ades_prec_out,ades_fit_out,ades_obs_tmp, &
                           & ades_flg_tmp,ades_prec_tmp,ades_fit_tmp,m,mt,change_ades)
             ades_error_model_out = ades_error_model_tmp
             IF(ades_error_model_out .NE. error_model) change_ades_errmod = .TRUE.
             change_ades = (change_ades .OR. change_ades_errmod)
             IF (change_ades) THEN
                IF(verb_io.gt.9)WRITE(*,*)'There are new/changed obs in ',ext,' file. New numobs=',m
             ELSE
                IF(verb_io.gt.9)WRITE(*,*)'There are no updates in ',ext,' file.'
                new_ades = .FALSE.
             ENDIF
          ELSE
             ! give the precedence to .rwo_psv, the .obs_psv file is intended as additional
             ! observations only; data in .rwo_psv are not erased, can only be changed
             ! read .rwo_psv, and store data in final array
             IF(verb_io.gt.9)WRITE(*,*)'Using .rwo_psv file, but checking .obs_psv for update.'
             ext = '.rwo_psv'
             ades_filename_out = psv_name(1:len_ades_name)//ext
             CALL read_ades_psv(ades_filename_out,ades_obs_out,ades_flg_out,ades_prec_out,ades_fit_out, &
               & nobs_ades_out,ntot_ades_out,ades_header_out,ades_key_rec_out,ades_data_rec_out, &
               & ntot_to_nobs_out,nobs_to_ntot_out,header_length_out,ades_error_model_out,new_key_rec_out)
             m = nobs_ades_out
             IF(m .EQ. 0) THEN
                obs0 = .FALSE.
             ELSE
                obs0 = .TRUE.
                IF(verb_io.gt.9)WRITE(*,*)'input_obs: ',m,' obs in ',trim(adjustl(ades_filename_out)) 
                WRITE(iun20,*)'input_obs: ',m,' obs in ',trim(adjustl(ades_filename_out))
             ENDIF
             !CALL convert_ades_fitobs(error_model,obs,obsw)
             IF(mpc_psv)THEN
                ! read .obs_psv and store in temporary array the data
                ades_filename_tmp = psv_name(1:len_ades_name)//'.obs_psv'
                CALL read_ades_psv(ades_filename_tmp,ades_obs_tmp,ades_flg_tmp,ades_prec_tmp,ades_fit_tmp, &
                     & nobs_ades_tmp,ntot_ades_tmp, ades_header_tmp,ades_key_rec_tmp,ades_data_rec_tmp, &
                     & ntot_to_nobs_tmp,nobs_to_ntot_tmp,header_length_tmp,ades_error_model_tmp,new_key_rec_tmp)
                mt = nobs_ades_tmp
                IF(mt .GT. 0) THEN
                   obs0 = .TRUE.
                ELSE
                   WRITE(*,*) trim(adjustl(ades_filename_tmp)),' is possibly corrupt. Not using any data from this file.'
                ENDIF
             ENDIF
             IF (.NOT. obs0) RETURN
             IF (mpc_psv) THEN
                IF(verb_io.gt.9)WRITE(*,*)'input_obs: ',mt,' obs in ',trim(adjustl(ades_filename_tmp)) 
                WRITE(iun20,*)'input_obs: ',mt,' obs in ',trim(adjustl(ades_filename_tmp))
                ! find weights for these; a priori RMS of astrometric observations
                !CALL convert_ades_fitobs(error_model,obst,obswt)
                ! add information from .obs_psv
                CALL update_ades_mpc(ades_obs_out,ades_flg_out,ades_prec_out,ades_fit_out,m,ntot_ades_out,ades_key_rec_out, &
                           & ades_data_rec_out,ntot_to_nobs_out,nobs_to_ntot_out,new_key_rec_out,ades_obs_tmp,ades_flg_tmp, &
                           & ades_prec_tmp,ades_fit_tmp,mt,ades_key_rec_tmp,ades_data_rec_tmp,nobs_to_ntot_tmp,change_ades,mnew)
                m = mnew
                nobs_ades_out = mnew
             ELSE
                change_ades = .TRUE.
             ENDIF
             ! if error model is changed, the data have to be considered changed
             IF(ades_error_model_out .NE. error_model) change_ades_errmod = .TRUE.
             change_ades = (change_ades .OR. change_ades_errmod)
             IF(change_ades)THEN
                IF(verb_io.gt.9)WRITE(*,*)'There are new/changed obs. New numobs=',m
             ELSE
                IF(verb_io.gt.9)WRITE(*,*)'There are no updates in .obs_psv file.'
                new_ades = .FALSE.
             ENDIF
          ENDIF
       ENDIF
       ! store public ADES arrays
       ades_filename = ades_filename_out
       nobs_ades = nobs_ades_out
       ntot_ades = ntot_ades_out
       ades_obs = ades_obs_out
       ades_flg = ades_flg_out
       ades_prec = ades_prec_out 
       ades_fit = ades_fit_out
       ades_header = ades_header_out
       ades_key_rec = ades_key_rec_out
       ades_data_rec = ades_data_rec_out
       header_length = header_length_out
       ntot_to_nobs_idx = ntot_to_nobs_out
       nobs_to_ntot_idx = nobs_to_ntot_out
       ades_error_model = ades_error_model_out
       new_key_rec = new_key_rec_out
       ! conversion of ADES data to internal data types (and computation of a priori RMS)
       CALL convert_ades_fitobs(error_model,obs,obsw)
       ! write output file according to input flag ades_type
       IF (ades_type_loc .GE. 0) THEN  ! write .rwo
          CALL write_rwo(file(1:lfile)//'.rwo',obs,obsw,m,error_model)
          !change = .TRUE.
       ENDIF
       IF (ades_type_loc .NE. 0 .AND. change_ades) THEN
          ! write output ADES file .rwo_psv
          ! convert data to ADES format (getting Residuals Group)
          CALL convert_fitobs_ades(obs,obsw,ades,.FALSE.)
          ! create output ADES file
          CALL write_ades(psv_name(1:len_ades_name),error_model,ades_type_loc)
       ENDIF
       change = change_ades
    ENDIF

    !==================================
    ! get asteroid radius (if necessary) before returning
    CALL aster_radius(obs%objdes,obs%type,m)
    ! =======  sort data before returning =========
    CALL heapsort(obs%time_tdt,m,iperm)
    ! copy output arrays into temp vectors
    obst(1:m) = obs(1:m)
    obswt(1:m) = obsw(1:m)
    ! copy back input output vectors in sorted order
    DO i = 1,m
       obs(i) = obst(iperm(i))
       obsw(i) = obswt(iperm(i))
    ENDDO
    ! sort ades arrays 
    IF (nobs_ades.GT.0) THEN
       CALL sort_ades(m,iperm)
    ENDIF
    ! re-initialize ades public variable to be used next
    new_ades = .FALSE.

  END SUBROUTINE input_obs

  ! SUBMODULE addobs
  ! =========================================
  !  A D D O B S _ M P C
  !
  ! add information from .obs and .rad to the one available from a file .rwo
  ! The number of observations can increase, from m to a maximum
  ! which is m+mt; however, the vectors have a maximum dimension nlef
  SUBROUTINE addobs_mpc(obs,obsw,m,obst,obswt,mt,mnew,change)
    ! ===============================================
    ! change flag, new obs. number
    LOGICAL, INTENT(OUT) :: change
    INTEGER, INTENT(OUT) ::  mnew
    ! ===== observational data ===========================
    INTEGER m ! observation number
    TYPE(ast_obs), INTENT(INOUT), DIMENSION(:) :: obs ! observations
    TYPE(ast_wbsr), INTENT(INOUT), DIMENSION(:) :: obsw ! weights and possibly residuals
    ! ===== observational data: temporary copy===========
    INTEGER mt ! observation number
    TYPE(ast_obs), INTENT(IN), DIMENSION(:) :: obst ! observations
    TYPE(ast_wbsr), INTENT(IN), DIMENSION(:) :: obswt ! weights and possibly residuals
    ! ===========================================
    INTEGER nlef ! length of obs,obsw vectors
    INTEGER j,mj,double
    ! ==================================================
    nlef=SIZE(obs)
    IF(nlef.ne.SIZE(obsw))THEN
       WRITE(*,*)' addobs_mpc: problems in dimensions of arrays obs, obsw ',nlef,SIZE(obsw)
       STOP
    ENDIF
    ! monitor changes
    mnew=m
    change=.false.
    ! scan supposedly new observations
    DO 1 j=1,mt
       mj=find_obs(obst(j),obs,m,double)
       IF(mj.ne.0.and.double.eq.0)THEN
          !     the observation was already there
          IF(changed_obs(obst(j),obswt(j),obs(mj),obsw(mj)))THEN
             ! ... but it is changed
             change=.true.
             obs(mj)=obst(j)
             obsw(mj)=obswt(j)
          ELSEIF(change_format.EQ.'INPUT') THEN
             change=.true.
             obs(mj)%catcodmpc=obst(j)%catcodmpc
             obsw(mj)%bias_coord=obswt(j)%bias_coord
             IF(.not.obsw(mj)%force_w(1)) THEN
                obsw(mj)%rms_coord(1)=obswt(j)%rms_coord(1)
             ENDIF
             IF(.not.obsw(mj)%force_w(2)) THEN
                obsw(mj)%rms_coord(2)=obswt(j)%rms_coord(2)
             ENDIF
          ELSEIF(obs(mj)%catcodmpc.NE.obst(j)%catcodmpc) THEN
             change=.true.
             obs(mj)%catcodmpc=obst(j)%catcodmpc
             obsw(mj)%bias_coord=obswt(j)%bias_coord
          ENDIF
       ELSEIF(mj.ne.0.and.double.ne.0)THEN
          ! the observation was already there, in double copy!
          IF(changed_obs(obst(j),obswt(j),obs(mj),obsw(mj)))THEN
             IF(changed_obs(obst(j),obswt(j),obs(double),obsw(double)))THEN
                change=.true.
                ! double, and changed! human intervention required
                WRITE(*,*)'addobs_mpc: double and changed'
                WRITE(*,*)' records ',mj,' and ',double,' in .rwo'
                WRITE(*,*)' record ',j,' in .obs'
                !             STOP
             ELSE
                ! OK, it is the double
             ENDIF
          ELSE
             ! OK, it is the first one
          ENDIF
       ELSEIF(mj.eq.0)THEN
          ! the observation is new: add it
          change=.true.
          !       WRITE(*,*)'addobs: new observation at record ',j
          mnew=mnew+1
          obs(mnew)=obst(j)
          obsw(mnew)=obswt(j)
       ENDIF
1   ENDDO
  END SUBROUTINE addobs_mpc
  ! =======================================================
  ! F I N D _ O B S
  ! find an observation from a list, matching the given one
  INTEGER FUNCTION find_obs(obst,obs,m,double)
    ! ============= INPUT =====================
    ! number of obs. record to be scanned, time, time to be found
    INTEGER m ! number of obs
    TYPE(ast_obs), INTENT(IN) :: obst ! one new observation
    TYPE(ast_obs),  INTENT(IN), DIMENSION(m) :: obs ! observations
    ! ============OUTPUT (impure func!) ===========
    INTEGER, INTENT(OUT) :: double ! location of the duplicate observation
    ! =========END INTERFACE=====================
    INTEGER j
    DOUBLE PRECISION, PARAMETER ::  epst= 1.d-8 ! time control
    find_obs=0
    double=0
    DO 1 j=1,m
       IF(abs(obst%time_utc-obs(j)%time_utc).lt.epst.and.                                &
            &      obst%obscod_i.eq.obs(j)%obscod_i.and.obst%type.eq.obs(j)%type)THEN
          ! observation found; is it a double?
          IF(find_obs.ne.0)THEN
             IF(double.eq.0)THEN
                double=j
                IF(ierrou.gt.0)THEN
                   IF(obst%tech.eq.'X'.or.obs(find_obs)%tech.eq.'X'.or.obs(j)%tech.eq.'X')THEN
                      !                   X obs. is often a duplicate, no problem
                   ELSE
                      WRITE(ierrou,*)'findob: two same time',find_obs,j,obst%time_utc,obst%obscod_s,obst%tech 
                      numerr=numerr+1
                   ENDIF
                ELSE
                   IF(verb_io.gt.9)WRITE(*,*)'findob: two same time', find_obs,j,obst%time_utc, obst%obscod_s
                ENDIF
             ELSE
                IF(ierrou.gt.0)THEN ! disaster case: triple observation!!!
                   WRITE(ierrou,*)'findob: three same time', &
                        &                     find_obs,double, j,obst%time_utc, obst%obscod_s, obst%tech 
                ELSE
                   WRITE(*,*)'findob: three same time',                &
                        &                     find_obs,double, j,obst%time_utc, obst%obscod_s, obst%tech 
                ENDIF
                !             STOP
             ENDIF
          ELSE
             find_obs=j !normal case, no double
          ENDIF
       ENDIF
1   ENDDO
  END FUNCTION find_obs
  ! =======================================================
  ! C H A N G E D _ O B S
  ! is an observation changed?
  LOGICAL FUNCTION changed_obs(obst,obswt,obs,obsw)
    TYPE(ast_obs), INTENT(IN) :: obst ! one new observation
    TYPE(ast_wbsr), INTENT(IN) :: obswt ! one new observation weight
    TYPE(ast_obs),  INTENT(IN) :: obs ! one old observations
    TYPE(ast_wbsr), INTENT(IN) :: obsw ! one old observation weight
    DOUBLE PRECISION::  epsa, epsw, epsb ! control on angle, weight, bias
    changed_obs=.false.
    epsa=1d-8
    IF(abs(obs%coord(1)-obst%coord(1)).gt.epsa*abs(obs%coord(1)).or.             &
         &      abs(obs%coord(2)-obst%coord(2)).gt.epsa*abs(obs%coord(2)))changed_obs=.true.
    ! temporary fix
    !      IF(obs%mag_str.ne.obst%mag_str.or.obs%tech.ne.obst%tech)changed_obs=.true.
    ! problem: if the weights are changed, should this be considered a changed observation??
    ! obsw, obswt are included to make this possible: but no
    !       epsw=1.d-2
    !       IF(abs(obsw%rms_coord(1)-obswt%rms_coord(1)).gt.epsw*abs(obsw%rms_coord(1)).or.      &
    !      &   abs(obsw%rms_coord(2)-obswt%rms_coord(2)).gt.epsw*abs(obsw%rms_coord(2)))       &
    !      & changed_obs=.true. 
    !       epsb=1.d-8 ! control on bias is absolute, because bias can be small!!
    !      IF(abs(obsw%bias_coord(1)-obswt%bias_coord(1)).gt.epsb.or.     &
    !      &   abs(obsw%bias_coord(2)-obswt%bias_coord(2)).gt.epsb) changed_obs=.true. 
  END FUNCTION changed_obs
  ! =========================================
  !  A D D O B S _ R W O
  !
  ! add information from .rwo file (only selection flags and manually fixed weights)
  ! to the observations from .obs and .rad
  ! The number of observations cannot increase, remains m
  SUBROUTINE addobs_rwo(obs,obsw,m,obst,obswt,mt,change)
    USE output_control
    ! logical change flag
    LOGICAl change
    ! ===== observational data ===========================
    INTEGER m ! observation number
    TYPE(ast_obs), INTENT(INOUT), DIMENSION(:) :: obs ! observations
    TYPE(ast_wbsr), INTENT(INOUT), DIMENSION(:) :: obsw ! weights and possibly residuals
    ! ===== observational data: temporary copy===========
    INTEGER mt ! observation number
    TYPE(ast_obs), INTENT(IN), DIMENSION(:) :: obst ! observations
    TYPE(ast_wbsr), INTENT(IN), DIMENSION(:) :: obswt ! weights and possibly residuals
    ! ===========================================
    !     INTEGER nlef ! length of obs,obsw vectors
    ! ===========================================
    INTEGER j,mj,double
    ! if all the same...
    change=.false.
    ! scan old weigts and selection flags, see if they match  observations
    DO 1 j=1,mt
       mj=find_obs(obst(j),obs,m,double)
       IF(mj.ne.0.and.double.eq.0)THEN
          ! this observation is still present in .obs, .rad files
          IF(changed_obs(obst(j),obswt(j),obs(mj),obsw(mj)))THEN
             ! ... but it is changed, thus leave selection flag=1 and default weight
             change=.true.
          ELSE
             ! if no change to observation then preserve the
             ! manually fixed weights and selection flags
             IF(obswt(j)%force_w(1))THEN
                obsw(mj)%rms_coord(1)=obswt(j)%rms_coord(1)
                obsw(mj)%force_w(1)=obswt(j)%force_w(1)
             ENDIF
             IF(obswt(j)%force_w(2))THEN
                obsw(mj)%rms_coord(2)=obswt(j)%rms_coord(2)
                obsw(mj)%force_w(2)=obswt(j)%force_w(2)
             ENDIF
             IF(obswt(j)%rms_mag.lt.0.d0)obsw(mj)%rms_mag=obswt(j)%rms_mag
             ! selection flags are preserved anyway
             obsw(mj)%sel_coord= obswt(j)%sel_coord
             obsw(mj)%sel_mag= obswt(j)%sel_mag
          ENDIF
       ELSEIF(mj.ne.0.and.double.ne.0)THEN
          ! the observation was already there, in double copy!
          IF(changed_obs(obst(j),obswt(j),obs(mj),obsw(mj)))THEN
             IF(changed_obs(obst(j),obswt(j),obs(double),obsw(double)))THEN
                change=.true.
                ! double, and changed! human intervention required
                WRITE(*,*)'addobs_rwo: double and changed'
                WRITE(ierrou,*)'addobs_rwo: double and changed'
                numerr=numerr+1
                WRITE(*,*)' records ',mj,' and ',double,' in .obs'
                WRITE(*,*)' record ',j,' in .rwo'
                WRITE(ierrou,*)' records ',mj,' and ',double,' in .obs'
                WRITE(ierrou,*)' record ',j,' in .rwo'
                !             STOP
             ELSE
                ! OK, it is the double
                ! it is the same, so preserve the selection flags
                ! manually fixed weights and selection flags
                IF(obswt(j)%force_w(1))THEN
                   obsw(double)%rms_coord(1)=obswt(j)%rms_coord(1)
                   obsw(double)%force_w(1)=obswt(j)%force_w(1)
                ENDIF
                IF(obswt(j)%force_w(2))THEN
                   obsw(double)%rms_coord(2)=obswt(j)%rms_coord(2)
                   obsw(double)%force_w(2)=obswt(j)%force_w(2)
                ENDIF
                IF(obswt(j)%rms_mag.lt.0.d0)obsw(mj)%rms_mag=obswt(j)%rms_mag
                obsw(double)%sel_coord=obswt(j)%sel_coord
                obsw(double)%sel_mag=obswt(j)%sel_mag
             ENDIF
          ELSE
             ! OK, it is the first one
             ! it is the same, so preserve the selection flags
             ! manually fixed weights and selection flags
             IF(obswt(j)%force_w(1))THEN
                obsw(mj)%rms_coord(1)=obswt(j)%rms_coord(1)
                obsw(mj)%force_w(1)=obswt(j)%force_w(1)
             ENDIF
             IF(obswt(j)%force_w(2))THEN
                obsw(mj)%rms_coord(2)=obswt(j)%rms_coord(2)
                obsw(mj)%force_w(2)=obswt(j)%force_w(2)
             ENDIF
             IF(obswt(j)%rms_mag.lt.0.d0)obsw(mj)%rms_mag=obswt(j)%rms_mag
             obsw(mj)%sel_coord=obswt(j)%sel_coord
             obsw(mj)%sel_mag=obswt(j)%sel_mag
          ENDIF
       ELSEIF(mj.eq.0)THEN
          ! if it is not found in .rwo, leave the default weights and selection flags
          change=.true.
       ENDIF
1   ENDDO
    ! check if there are extra (added) observations in .obs file
    ! it might be better to loop on the .obs rather than the .rwo data
    IF (mt.lt.m) change=.true.
  END SUBROUTINE addobs_rwo

  ! Copyright Orbfit Consortium 1999
  ! this routine determines an asteroid radius if needed for radar
  ! Revised by Genny on 5 November 2005 to deal with the new format of 
  ! astorb.dat to include numbered objects > 99999
  SUBROUTINE aster_radius(objid,obstype,m)
    INTEGER,INTENT(IN) :: m
    CHARACTER*(name_len), INTENT(IN), DIMENSION(m) ::  objid  ! designation
    CHARACTER*1, INTENT(iN), DIMENSION(m) ::  obstype ! find if radar
    ! OUTPUT: radius of asteroid is public radius
    ! ====================================
    integer unit,ln,lnam,lnum,isav,i
    character*18 oid,number,name
    double precision hmag,diam,exponent
    logical needed

    ! First find out if we need the radius at all
    needed=.false.
    do i=1,m
       if(obstype(i).eq.'R'.or.obstype(i).eq.'V')then
          needed=.true.
          isav=i
          exit
       endif
    enddo
    if(.not.needed)then
       ! bail out, leaving as it is.
       !         write (*,*) 'Asteroid radius not required.'
       !         radius=-1d0
       return
    else
       ! get radius
       call filopl(unit,'astorb.rad')
       oid=objid(isav)
       call rmsp(oid,ln)
       ! read records in an infinite loop (it's F(UGLY)77)
1      continue
       read(unit,101,end=109) number,name,hmag,diam
       call rmsp(number,lnum)
       call rmsp(name,lnam)
       if(oid(1:ln).eq.number(1:lnum) .or.                       &
            &           oid(1:ln).eq.name(1:lnam)    )THEN
          !              found object in file
          if(diam.gt.0d0)then
             radius=diam/2d0/1.4998d8
             IF(verb_io.gt.9)write(*,102)'RADAR: Using radius = ',                 &
                  &                 radius*1.4998d8,' km from IRAS diameter.'
          else
             exponent=0.5d0*(6.3d0-log10(0.2d0)-0.4d0*hmag)
             radius=10**exponent/2d0/1.4998d8
             IF(verb_io.gt.9)write(*,102)'RADAR: Using radius = ',                 &
                  &                 radius*1.4998d8,' km from absolute magnitude.'
          endif
          ! exit infinite loop
          goto 109
       endif
       goto 1
    endif
    write(*,*)'*** astrad warning: ',name,' not found in astorb.dat.'
    radius=1d0/1.4998d8
    write(*,*)'RADAR: Using radius = ',radius*1.4998d8,' km.'
109 call filclo(unit,' ')
    ! ==============================================================
    !                  num   nam    comp   Hmag  Gmag    col   diam
101 format(a6,1X,A18,1X,15x,1X,f5.2,1X,5x,1X,4x,1X,f5.1)
102 format(a,f6.2,a)
    ! ==============================================================
  END SUBROUTINE aster_radius
  !
  ! Copyright (C) 1997 by Mario Carpino
  ! Version: February 24, 1997
  !
  !  *****************************************************************
  !  *                                                               *
  !  *                         R D A N G A                           *
  !  *                                                               *
  !  *    Read an angle with its accuracy from a character string    *
  !  *                                                               *
  !  *****************************************************************
  !
  ! INPUT:    STRING    -  Character string
  !
  ! OUTPUT:   ANGLE     -  Angle (no unit conversion is performed)
  !           ACC       -  Angle accuracy
  !           ERROR     -  Conversion error
  !
  SUBROUTINE rdanga(string,angle,acc,error)
    CHARACTER*(*) string
    DOUBLE PRECISION angle,acc
    LOGICAL error

    ! Max string length
    INTEGER lx
    PARAMETER (lx=200)

    CHARACTER*(lx) c1,c,field
    INTEGER l,isig,nf,i,pp,iv,ll
    LOGICAL nospli
    DOUBLE PRECISION fact,rv


    INTEGER lench,nitchs
    EXTERNAL lench,nitchs

    error=.true.
    IF(lench(string).GT.lx) STOP '**** rdanga: LEN(string) > lx ****'
    c1=string
    CALL norstr(c1,l)

    ! The sign may be separated from the value by blank space
    isig=1
    IF(c1(1:1).EQ.'+') THEN
       c=c1(2:)
       CALL norstr(c,l)
    ELSEIF(c1(1:1).EQ.'-') THEN
       isig=-1
       c=c1(2:)
       CALL norstr(c,l)
    ELSE
       c=c1
    END IF

    nf=nitchs(c)
    IF(nf.LT.1.OR.nf.GT.3) RETURN
    angle=0
    fact=1
    DO 1 i=1,nf
       CALL stspli(c,' ',field,nospli)
       IF(nospli) RETURN
       pp=INDEX(field,'.')
       IF(pp.GT.0) THEN
          IF(i.NE.nf) RETURN
          READ(field,*,err=10,end=10) rv
       ELSE
          READ(field,*,err=10,end=10) iv
          rv=iv
       END IF
       angle=angle+fact*rv
       IF(i.EQ.nf) THEN
          CALL norstr(field,ll)
          pp=INDEX(field,'.')
          IF(pp.EQ.0) THEN
             acc=1
          ELSE
             ll=lench(field)
             acc=10.0d0**(pp-ll)
          END IF
          acc=acc*fact
       ELSE
          fact=fact/60.d0
       END IF
1   END DO
    angle=angle*isig
    error=.false.

10  CONTINUE
  END SUBROUTINE rdanga

  !SUB-MODULE observ_bias_weight
  ! Version 5.0 Copyright OrbFit consortium 2018
  ! ==========================================================================
  ! OBSERV_RMS
  !     Computation of bias from catalogs and of a-priori RMS of observations
  !
  ! ==========================================================================
  !
  ! INPUT:    OBS          - observations
  !           ERROR_MODEL  - name of weighing scheme to be used
  !           N            - number of observations
  !           INIT         - reinitialize all the structure (when init=.false.,
  !                          information different from weight/bias is preserved)
  !
  ! OUTPUT:   OBSW         - observation weights
  !
  SUBROUTINE observ_rms(obs,error_model,init,obsw,n)
    INTEGER       ,               INTENT(IN)    :: n           ! number of obs
    LOGICAL       ,               INTENT(IN)    :: init        ! if true reset obsw 
    TYPE(ast_obs) , DIMENSION(n), INTENT(INOUT) :: obs         ! observations
    CHARACTER*(*) ,               INTENT(IN)    :: error_model ! code for errmod
    TYPE(ast_wbsr), DIMENSION(n), INTENT(INOUT) :: obsw        ! when init=.false.,
    ! ============================== END INTERFACE =====================================
    INTEGER i,j,k,jj,iline,it,ic,ip
    INTEGER unit ! for input of RMS values
    INTEGER low_year,low_month,low_day,up_year,up_month,up_day
    INTEGER ntech, ncat, nprog, lench, len_ld, len_ud
    DOUBLE PRECISION rmsrmi,rmsvmi,rmstmp(2),t2,lower_mjd,upper_mjd
    ! DOUBLE PRECISION, DIMENSION(2,nrulesx)    :: mjd_date_rules
    ! Get default magnitude weights from a function magrms_new in same module
    CHARACTER q*4, catcod*1, qq*3
    CHARACTER obscod*3, tech_type*5, catcodstr*10, program_code*10, lower_date*10, upper_date*10
    INTEGER ind(n),nob(n), nob1(n) !sorting by time, number obs in batch
    INTEGER ios ! for errors and end of file management
    TYPE(ast_obs), DIMENSION(n) :: obs1
    DOUBLE PRECISION, PARAMETER :: gapmax=8.d0/24.d0 ! 8 hours max gap
    DOUBLE PRECISION, EXTERNAL :: tjm1
    ! ==================================================================================
    LOGICAL, SAVE :: first
    DATA first/.true./
    ! ======================= initialize arrays of RMS rules ===========================
    IF(first)THEN
       IF(error_model.EQ.'vfcc17')THEN
          CALL filopl(unit,'weights/vfcc17.rules')
          nrules=0
          iline=0
          ERR:DO
             READ(unit, 301, IOSTAT=ios)obscod,tech_type,catcodstr,program_code,lower_date,upper_date,rmstmp
301          FORMAT(A3,3X,A5,3X,A10,3X,A3,3X,A10,3X,A10,3X,F5.2,2X,F5.2)
             IF(ios > 0) THEN
                WRITE(*,*) "Error in reading weights/vfcc17.rules"
             ELSE IF(ios < 0) THEN
                EXIT ERR
             ELSE
                CALL rmsp(tech_type,ntech)
                CALL rmsp(catcodstr,ncat)
                CALL rmsp(program_code,nprog)
                IF(lower_date(1:4).EQ.'    ')THEN
                   lower_mjd=-99999.d0
                ELSE
                   READ(lower_date(9:10),'(I2)') low_day
                   READ(lower_date(6:7) ,'(I2)') low_month
                   READ(lower_date(1:4) ,'(I4)') low_year
                   lower_mjd=tjm1(low_day,low_month,low_year,0.d0)
                ENDIF
                IF(upper_date(1:4).EQ.'    ')THEN
                   upper_mjd=99999.d0
                ELSE
                   READ(upper_date(9:10),'(I2)') up_day
                   READ(upper_date(6:7) ,'(I2)') up_month
                   READ(upper_date(1:4) ,'(I4)') up_year
                   upper_mjd=tjm1(up_day,up_month,up_year,0.d0)
                ENDIF
                DO it=1, ntech
                   DO ic=1, ncat
                      IF(nprog.EQ.0) THEN
                         iline=iline+1
                         nprog=1
                         errmod_rules(iline)=obscod//tech_type(it:it)//catcodstr(ic:ic)//' '
                         mjd_date_rules(1:2,iline)=(/lower_mjd,upper_mjd/)
                         rmsbin(1:2,iline)=rmstmp
                      ELSE
                         DO ip=1, nprog
                            iline=iline+1
                            errmod_rules(iline)=obscod//tech_type(it:it)//catcodstr(ic:ic)//program_code(ip:ip)
                            mjd_date_rules(1:2,iline)=(/lower_mjd,upper_mjd/)
                            rmsbin(1:2,iline)=rmstmp
                         ENDDO
                      ENDIF
                   ENDDO
                ENDDO
                nrules=iline
             ENDIF
          ENDDO ERR
97        CALL filclo(unit,' ')
          !       CALL heapsortnachlo(obscatcod,nrules,4,q,indrules)       
       ELSEIF(error_model.EQ.'fcct14')THEN
          CALL filopl(unit,'weights/fcct14.rules')
          nrules=0
          DO j=1, nrulesx
             READ(unit, 101, END=89)obscod,catcodstr,rmstmp
101          FORMAT(A3,3X,A6,2X,F5.2,2X,F5.2)
             ncat=lench(catcodstr)-2
             DO i=1,ncat
                obscatcod(nrules+i)=obscod//catcodstr(2+i:2+i)
                rmsbin(1:2,nrules+i)=rmstmp
             ENDDO
             nrules=nrules+ncat
          ENDDO
89        CALL filclo(unit,' ')
          CALL heapsortnachlo(obscatcod,nrules,4,q,indrules)
       ELSEIF(error_model.EQ.'cbm10')THEN
          ! rules for observatory, whatever star catalog
          CALL filopl(unit,'weights/cbm10.sta')
          READ(unit,*) ! skip header
          nstarms=0
          DO j=1, nstarmsx
             READ(unit, *, END=98)obscod,rmsstabin(1:2,j)
             ! example: 244        0.20       0.20
             obscodrms(j)=obscod
             nstarms=nstarms+1
          ENDDO
98        CALL filclo(unit,' ')
          CALL heapsortnachlo(obscodrms,nstarms,3,qq,indstarms)
          ! rules for observatory and star catalog
          CALL filopl(unit,'weights/cbm10.rules')
          READ(unit,*) ! skip header
          nrules=0
          DO j=1, nrulesx
             READ(unit, 100, END=99)obscod,catcod,rmsbin(1:2,j)
             obscatcod(j)=obscod//catcod
             ! example: 704:cc @ 1.23, 1.19
100          FORMAT(A3,2X,A1,2X,F5.2,1X,F5.2)
             nrules=nrules+1
          ENDDO
99        CALL filclo(unit,' ')
          CALL heapsortnachlo(obscatcod,nrules,4,q,indrules)
       ENDIF
       first=.false.
    ENDIF
    ! for each observations build the corresponding weights and bias record
    DO  i=1,n
       IF(init) obsw(i)=undefined_ast_wbsr
       IF((obs(i)%type.EQ.'O'.OR.obs(i)%type.EQ.'S').AND.&
            & (error_model.EQ.'cbm10'.OR.error_model.EQ.'fcct14'.OR.error_model.EQ.'vfcc17'))THEN
          CALL astrow_bias(error_model,obs(i),obsw(i))
          ! magnitude weigthing
          obsw(i)%rms_mag=magrms_new(obs(i)%mag_str,obs(i)%time_tdt,obs(i)%obscod_i,obs(i)%tech)
       ELSEIF(obs(i)%type.EQ.'O'.OR.obs(i)%type.EQ.'S')THEN
          ! astrometric (telescope) observations; get weight
          CALL astrow_nobias(error_model,obs(i),obsw(i))
          ! magnitude weigthing
          obsw(i)%rms_mag=magrms_new(obs(i)%mag_str,obs(i)%time_tdt,obs(i)%obscod_i,obs(i)%tech)
       ELSEIF(obs(i)%type.EQ.'R')THEN
          ! weighting of radar data is given with observations, but there is
          ! minimum credible (for a given acccuracy of the models)
          rmsrmi=1.d-3/aukm
          obsw(i)%rms_coord(1)=MAX(obs(i)%acc_coord(1),rmsrmi)
          obsw(i)%rms_coord(2)=9.d9
          obsw(i)%rms_mag=9.d9
          obsw(i)%bias_coord(2)=0.d0
       ELSEIF(obs(i)%type.EQ.'V')THEN
          rmsvmi=1.d-3/aukm
          obsw(i)%rms_coord(2)=MAX(obs(i)%acc_coord(2),rmsvmi)
          obsw(i)%rms_coord(1)=9.d9
          obsw(i)%rms_mag=9.d9
          obsw(i)%bias_coord(1)=0.d0
       ELSE
          WRITE(*,*)'observ_rms: obs. type', obs(i)%type,' not known ', &
               &           'at time ',obs(i)%time_tdt,' number ',i
          STOP
       ENDIF
       IF(obs(i)%tech.EQ.'x'.OR.obs(i)%tech.EQ.'X') THEN
          obsw(i)%sel_coord=0
       END IF
    END DO
    ! apply sqrt(nobs) factor within one batch (ending with an 8 hour gap).
    ! sort by time first
    CALL heapsort(obs(1:n)%time_tdt,n,ind)
    DO i=1,n
       obs1(i)=obs(ind(i))
    ENDDO
    ! find number of observations in batch from same station
    nob1=0  
    DO  j=1,n
       IF(obs1(j)%type.EQ.'O'.AND.(error_model.EQ.'fcct14'.OR.error_model.EQ.'vfcc17'))THEN
          !IF(obs1(j)%type.EQ.'O'.AND.error_model.EQ.'fcct14')THEN
          IF(nob1(j).GT.0)CYCLE ! already done
          t2=obs1(j)%time_tdt
          nob1(j)=1
          IF(obs1(j)%tech.EQ.'X'.OR.obs1(j)%tech.EQ.'x')CYCLE
          ! skip also obs having nonexistent technology descriptor 
          IF(ALL(tech_values.NE.obs1(j)%tech)) CYCLE
          ! find batch from same station
          DO k=j+1,n
             IF(obs1(k)%time_tdt-t2.GT.gapmax) EXIT ! batch finished
             IF(obs1(k)%obscod_s.EQ.obs1(j)%obscod_s)THEN
                IF (obs1(k)%tech.EQ.'X'.OR.obs1(k)%tech.EQ.'x') CYCLE
                IF(ALL(tech_values.NE.obs1(k)%tech)) CYCLE
                nob1(j)=nob1(j)+1
                t2=obs1(k)%time_tdt
             ENDIF
          ENDDO
          ! assign number of obs in batch to each obs in batch
          DO k=j+1,n
             IF(obs1(k)%time_tdt-t2.GT.gapmax) EXIT ! batch finished
             IF(obs1(k)%obscod_s.EQ.obs1(j)%obscod_s)THEN
                IF (obs1(k)%tech.EQ.'X'.OR.obs1(k)%tech.EQ.'x') CYCLE
                IF(ALL(tech_values.NE.obs1(k)%tech)) CYCLE
                nob1(k)=nob1(j)
             ENDIF
          ENDDO
       ENDIF
    ENDDO
    !  DO i=1,n
    !     WRITE(*,*) nob1(i)
    !  END DO
    !  STOP
    ! go back to original array of obs., apply sqrt
    DO k=1,n
       jj=ind(k)
       nob(jj)=nob1(k)
       !IF(nob(jj).GE.3)THEN
       IF(nob(jj).NE.0)THEN
          IF(error_model.EQ.'vfcc17') THEN
             IF(nob(jj).GE.5) THEN
                IF(.NOT.obsw(jj)%force_w(1))THEN
                   obsw(jj)%rms_coord(1)=obsw(jj)%rms_coord(1)*SQRT(nob(jj)*0.25d0)
                END IF
                IF(.NOT.obsw(jj)%force_w(2))THEN
                   obsw(jj)%rms_coord(2)=obsw(jj)%rms_coord(2)*SQRT(nob(jj)*0.25d0)
                END IF
             ENDIF
          ELSE
             IF(.NOT.obsw(jj)%force_w(1))THEN
                obsw(jj)%rms_coord(1)=obsw(jj)%rms_coord(1)*sqrt(1.d0*nob(jj))
             END IF
             IF(.NOT.obsw(jj)%force_w(2))THEN
                obsw(jj)%rms_coord(2)=obsw(jj)%rms_coord(2)*sqrt(1.d0*nob(jj))
             END IF
          ENDIF
       ENDIF
    ENDDO

  END SUBROUTINE observ_rms

  ! ------------------------------------------------------ !
  ! SUBROUTINE observ_rms_singobs                          !
  ! It retrieves a-priori RMS and bias from catalogs       !
  ! of a single observation.                               !
  ! ------------------------------------------------------ !
  ! Author: Alessia Bertolucci, SpaceDyS 2020              !
  ! ------------------------------------------------------ !

   SUBROUTINE observ_rms_singobs(obs,error_model,init,obsw)
    LOGICAL,               INTENT(IN)    :: init        ! if true reset obsw 
    TYPE(ast_obs),         INTENT(INOUT) :: obs         ! observation
    CHARACTER*(*),         INTENT(IN)    :: error_model ! code for errmod
    TYPE(ast_wbsr),        INTENT(INOUT) :: obsw        ! when init=.false.
   
    INTEGER  :: i,j,k,jj,iline,it,ic,ip
    INTEGER  :: unit                        ! for input of RMS values
    INTEGER  :: ios                         ! for errors and end of file management
    INTEGER  :: low_year,low_month,low_day,up_year,up_month,up_day
    INTEGER  :: ntech, ncat, nprog, lench, len_ld, len_ud

    REAL(KIND=dkind) :: rmsrmi,rmsvmi,rmstmp(2),t2,lower_mjd,upper_mjd
    CHARACTER        :: q*4, catcod*1, qq*3
    CHARACTER        :: obscod*3, tech_type*5, catcodstr*10, program_code*10, lower_date*10, upper_date*10
    
    REAL(KIND=dkind), EXTERNAL  :: tjm1
    ! ==================================================================================
    LOGICAL, SAVE :: first
    DATA first/.true./
    ! ======================= initialize arrays of RMS rules ===========================
    IF(first)THEN
       IF(error_model.EQ.'vfcc17')THEN
          CALL filopl(unit,'weights/vfcc17.rules')
          nrules=0
          iline=0
          ERR:DO
             READ(unit, 301, IOSTAT=ios)obscod,tech_type,catcodstr,program_code,lower_date,upper_date,rmstmp
301          FORMAT(A3,3X,A5,3X,A10,3X,A3,3X,A10,3X,A10,3X,F5.2,2X,F5.2)
             IF(ios > 0) THEN
                WRITE(*,*) "Error in reading weights/vfcc17.rules"
             ELSE IF(ios < 0) THEN
                EXIT ERR
             ELSE
                CALL rmsp(tech_type,ntech)
                CALL rmsp(catcodstr,ncat)
                CALL rmsp(program_code,nprog)
                IF(lower_date(1:4).EQ.'    ')THEN
                   lower_mjd=-99999.d0
                ELSE
                   READ(lower_date(9:10),'(I2)') low_day
                   READ(lower_date(6:7) ,'(I2)') low_month
                   READ(lower_date(1:4) ,'(I4)') low_year
                   lower_mjd=tjm1(low_day,low_month,low_year,0.d0)
                ENDIF
                IF(upper_date(1:4).EQ.'    ')THEN
                   upper_mjd=99999.d0
                ELSE
                   READ(upper_date(9:10),'(I2)') up_day
                   READ(upper_date(6:7) ,'(I2)') up_month
                   READ(upper_date(1:4) ,'(I4)') up_year
                   upper_mjd=tjm1(up_day,up_month,up_year,0.d0)
                ENDIF
                DO it=1, ntech
                   DO ic=1, ncat
                      IF(nprog.EQ.0) THEN
                         iline=iline+1
                         nprog=1
                         errmod_rules(iline)=obscod//tech_type(it:it)//catcodstr(ic:ic)//' '
                         mjd_date_rules(1:2,iline)=(/lower_mjd,upper_mjd/)
                         rmsbin(1:2,iline)=rmstmp
                      ELSE
                         DO ip=1, nprog
                            iline=iline+1
                            errmod_rules(iline)=obscod//tech_type(it:it)//catcodstr(ic:ic)//program_code(ip:ip)
                            mjd_date_rules(1:2,iline)=(/lower_mjd,upper_mjd/)
                            rmsbin(1:2,iline)=rmstmp
                         ENDDO
                      ENDIF
                   ENDDO
                ENDDO
                nrules=iline
             ENDIF
          ENDDO ERR
97        CALL filclo(unit,' ')
           
       ELSEIF(error_model.EQ.'fcct14')THEN
          CALL filopl(unit,'weights/fcct14.rules')
          nrules=0
          DO j=1, nrulesx
             READ(unit, 101, END=89)obscod,catcodstr,rmstmp
101          FORMAT(A3,3X,A6,2X,F5.2,2X,F5.2)
             ncat=lench(catcodstr)-2
             DO i=1,ncat
                obscatcod(nrules+i)=obscod//catcodstr(2+i:2+i)
                rmsbin(1:2,nrules+i)=rmstmp
             ENDDO
             nrules=nrules+ncat
          ENDDO
89        CALL filclo(unit,' ')
          CALL heapsortnachlo(obscatcod,nrules,4,q,indrules)
       ELSEIF(error_model.EQ.'cbm10')THEN
          ! rules for observatory, whatever star catalog
          CALL filopl(unit,'weights/cbm10.sta')
          READ(unit,*) ! skip header
          nstarms=0
          DO j=1, nstarmsx
             READ(unit, *, END=98)obscod,rmsstabin(1:2,j)
             ! example: 244        0.20       0.20
             obscodrms(j)=obscod
             nstarms=nstarms+1
          ENDDO
98        CALL filclo(unit,' ')
          CALL heapsortnachlo(obscodrms,nstarms,3,qq,indstarms)
          ! rules for observatory and star catalog
          CALL filopl(unit,'weights/cbm10.rules')
          READ(unit,*) ! skip header
          nrules=0
          DO j=1, nrulesx
             READ(unit, 100, END=99)obscod,catcod,rmsbin(1:2,j)
             obscatcod(j)=obscod//catcod
             ! example: 704:cc @ 1.23, 1.19
100          FORMAT(A3,2X,A1,2X,F5.2,1X,F5.2)
             nrules=nrules+1
          ENDDO
99        CALL filclo(unit,' ')
          CALL heapsortnachlo(obscatcod,nrules,4,q,indrules)
       ENDIF
       first=.false.
    ENDIF
    ! for the observation build the corresponding weights and bias record
    IF(init) obsw = undefined_ast_wbsr
    IF((obs%type.EQ.'O'.OR.obs%type.EQ.'S').AND.&
         & (error_model.EQ.'cbm10'.OR.error_model.EQ.'fcct14'.OR.error_model.EQ.'vfcc17'))THEN
       CALL astrow_bias(error_model,obs,obsw)
    ELSEIF(obs%type.EQ.'O'.OR.obs%type.EQ.'S')THEN
       ! astrometric (telescope) observations; get weight
       CALL astrow_nobias(error_model,obs,obsw)
    ELSEIF(obs%type.EQ.'R')THEN
       ! weighting of radar data is given with observations, but there is
       ! minimum credible (for a given acccuracy of the models)
       rmsrmi=1.d-3/aukm
       obsw%rms_coord(1)=MAX(obs%acc_coord(1),rmsrmi)
       obsw%rms_coord(2)=9.d9
       obsw%rms_mag=9.d9
       obsw%bias_coord(2)=0.d0
    ELSEIF(obs%type.EQ.'V')THEN
       rmsvmi=1.d-3/aukm
       obsw%rms_coord(2)=MAX(obs%acc_coord(2),rmsvmi)
       obsw%rms_coord(1)=9.d9
       obsw%rms_mag=9.d9
       obsw%bias_coord(1)=0.d0
    ELSE
       WRITE(*,*)'observ_rms_singobs: obs. type', obs%type,' not known ', &
            &           'at time ',obs%time_tdt
       STOP
    ENDIF
    IF(obs%tech.EQ.'x'.OR.obs%tech.EQ.'X') THEN
       obsw%sel_coord=0
    END IF
   
  END SUBROUTINE observ_rms_singobs

  ! ============================================================================
  ! ASTROW_BIAS
  !     Computation of bias from catalogs, and of a-priori RMS of observations.
  !     Error models:
  !          CBM10 from Chesley, Baer and Monet (2010)
  !          FCCT14 from Farnocchia, Chesley, Chamberlin, Tholen (2014)
  !          VFCC17 from Veres, Farnocchia, Chesley, Chamberlin (2017)
  ! Last version introduces VFCC17 with different format in weight/vfcc17.rules
  ! wrt the previous error model.
  !
  ! Written by Andrea Chessa @SpaceDyS 2017
  ! Updated by Fabrizio Bernardi @SpaceDyS 2018
  ! ============================================================================
  !
  ! INPUT:    OBS          - observations
  !           ERROR_MODEL  - name of weighing scheme to be used
  !           N            - number of observations
  !           INIT         - reinitialize all the structure (when init=.false.,
  !                          information different from weight/bias is preserved)
  !
  ! OUTPUT:   OBSW         - observation weights
  ! ============================================================================
  SUBROUTINE astrow_bias(error_model,obs,obsw)
    USE cat_debias
    IMPLICIT NONE
    TYPE(ast_obs)   , INTENT(INOUT) :: obs
    CHARACTER*(*)   , INTENT(IN)    :: error_model
    TYPE(ast_wbsr)  , INTENT(INOUT) :: obsw ! when init=.false.,
    ! information different from weight/bias is preserved
    ! =============================end interface==================================
    CHARACTER*2 catcode
    DOUBLE PRECISION                :: bias_ra, bias_dec, pmRA, pmDEC, rmsmin, tdt, rmsa, rmsd, deltat
    CHARACTER                       :: observ*3, obscat*4, final_rule*6
    CHARACTER                       :: tech, catalog, prog_code
    LOGICAL                         :: ermrms, biasflag, first_err
    INTEGER                         :: i, j, ioc, iocall, indobs, num_non_zero, idx_non_zero
    INTEGER                         :: num_star, num_blank, final_star, final_blank, nloop
    INTEGER                         :: idx_selected(nrules), idx_selected_all(nrules)         
    DOUBLE PRECISION, PARAMETER     :: mjd2000=51544.d0
    ! ============================================================================

    IF(error_model.NE.'cbm10'.AND.error_model.NE.'fcct14'.AND.error_model.NE.'vfcc17')THEN
       WRITE(*,*)  '**** WRONG ERROR MODEL ***** ',error_model
       STOP
    ENDIF
    first_err=.TRUE.
    tdt=obs%time_utc
    observ=obs%obscod_s
    tech=obs%tech
    catalog=obs%catcodmpc(2:2)
    prog_code=obs%note
    ! years since 2000
    deltat=(tdt-mjd2000)/365.25d0
    ! find bias
    CALL cat_bias(obs%coord(1),obs%coord(2),obs%catcodmpc(2:2),bias_ra,bias_dec,pmRA,pmDEC,biasflag,error_model)
    ! record bias in ast_wbsr data type
    IF(biasflag)THEN
       ! pmRA and pmDEC [mas/y], bias_ra, bias_dec [arcsec]
       obsw%bias_coord(1)=(bias_ra+pmRA*deltat/1000.d0)*radsec/cos(obs%coord(2)) !!!Checked with Steve
       obsw%bias_coord(2)=(bias_dec+pmDEC*deltat/1000.d0)*radsec
    ELSE
       obsw%bias_coord(1:2)=0.d0
    ENDIF
    !####################!
    ! ERROR MODEL VFCC17 !
    !####################!
    IF(error_model.eq.'vfcc17')THEN
       IF(obs%tech.eq.'x'.or.obs%tech.eq.'X')THEN
          rmsa=999.d0
          rmsd=999.d0
          GOTO 330
       ENDIF
       ! ioc = index obs code
       ! idx_selected index containig the obscode
       ioc=0
       iocall=0
       idx_selected=0
       idx_non_zero=1
       ! selection of all records containing the obscode
       DO j=1, nrules
          IF(index(errmod_rules(j),observ).NE.0) THEN
             ioc=ioc+1
             idx_selected(ioc)=j ! store the selected index containig the obscode
          ENDIF
          IF(index(errmod_rules(j),'ALL').NE.0) THEN
             iocall=iocall+1
             idx_selected_all(iocall)=j ! store the selected index containig the 'ALL' obscode
             !IF(index(errmod_rules(j),'C').NE.0) THEN
             !   IF((index(errmod_rules(j),'U').LE.0).AND.(index(errmod_rules(j),'V').LE.0)) THEN
             !      gen_ind_all=j ! This is for the generic astrometric data
             !   ELSE
             !      gen_ind_gaia=j  ! This is for the generic astrometric Gaia catalog data
             !   ENDIF
             !ENDIF
          ENDIF
       ENDDO
       nloop=0
       ! if the record is not found, we select the standard for ALL
300    IF(ioc.EQ.0) THEN
          observ='ALL'
          ioc=iocall
          DO i=1, iocall
             idx_selected(i)=idx_selected_all(i)
          ENDDO
       ENDIF
       ! select gen_ind for generic astrometric data
       !IF((index(catalog,'U').NE.0).OR.(index(catalog,'V').NE.0)) THEN
       !   gen_ind=gen_ind_gaia
       !ELSE
       !   gen_ind=gen_ind_all
       !ENDIF
       !idx_non_zero=COUNT(idx_selected.NE.0)
       ! read the selected rules and check the specific field (tech, cat, prog), line by line.
       ! if not the index becomes zero
       ! the mj_date comparison occurs after the selection of the rules
       idx_non_zero=ioc
       DO i=1, idx_non_zero
          IF(tdt.GT.mjd_date_rules(1,idx_selected(i)).AND.tdt.LE.mjd_date_rules(2,idx_selected(i))) THEN
             IF(index(errmod_rules(idx_selected(i)),observ//tech//catalog//prog_code).NE.0) THEN
                rmsa=rmsbin(1,idx_selected(i))
                rmsd=rmsbin(2,idx_selected(i))
                GOTO 330
             ENDIF
          ENDIF
       ENDDO
       DO i=1, idx_non_zero
          IF(tdt.GT.mjd_date_rules(1,idx_selected(i)).AND.tdt.LE.mjd_date_rules(2,idx_selected(i))) THEN
             IF(index(errmod_rules(idx_selected(i)),observ//tech//catalog//' ').NE.0) THEN
                rmsa=rmsbin(1,idx_selected(i))
                rmsd=rmsbin(2,idx_selected(i))
                GOTO 330
             ENDIF
          ENDIF
       ENDDO
       DO i=1, idx_non_zero
          IF(tdt.GT.mjd_date_rules(1,idx_selected(i)).AND.tdt.LE.mjd_date_rules(2,idx_selected(i))) THEN
             IF(index(errmod_rules(idx_selected(i)),observ//tech//'*'//prog_code).NE.0) THEN
                rmsa=rmsbin(1,idx_selected(i))
                rmsd=rmsbin(2,idx_selected(i))
                GOTO 330
             ENDIF
          ENDIF
       ENDDO
       DO i=1, idx_non_zero
          IF(tdt.GT.mjd_date_rules(1,idx_selected(i)).AND.tdt.LE.mjd_date_rules(2,idx_selected(i))) THEN
             IF(index(errmod_rules(idx_selected(i)),observ//tech//'* ').NE.0) THEN
                rmsa=rmsbin(1,idx_selected(i))
                rmsd=rmsbin(2,idx_selected(i))
                GOTO 330
             ENDIF
          ENDIF
       ENDDO
       DO i=1, idx_non_zero
          IF(tdt.GT.mjd_date_rules(1,idx_selected(i)).AND.tdt.LE.mjd_date_rules(2,idx_selected(i))) THEN
             IF(index(errmod_rules(idx_selected(i)),observ//'**'//prog_code).NE.0) THEN
                rmsa=rmsbin(1,idx_selected(i))
                rmsd=rmsbin(2,idx_selected(i))
                GOTO 330
             !ELSE
                !rmsa=rmsbin(1,gen_ind)
                !rmsd=rmsbin(2,gen_ind)
                !ioc=0
                !GOTO 300
             ENDIF
          !ELSE
             !rmsa=rmsbin(1,gen_ind)
             !rmsd=rmsbin(2,gen_ind)
             !ioc=0
             !GOTO 300
          ENDIF
       ENDDO
       nloop = nloop+1
       IF(nloop.LT.2) THEN
          ! the record is not found, set observ='ALL' and repeat the research
          ioc=0
          GOTO 300
       ELSE
          ! the record is still not found, assign rms for generic astrometric data
          rmsa=1.d0
          rmsd=1.d0
          GOTO 330
       ENDIF

       !####################!
       ! ERROR MODEL FCCT14 !
       !####################!
    ELSEIF(error_model.eq.'fcct14')THEN
       ! first set default based on observation technique
       IF(obs%tech.eq.'P'.or.obs%tech.eq.'A'.or.obs%tech.eq.'N')THEN
          ! Photografic: Default value of RMS (time dependent, in arcsec)
          ! Before 1890
          IF(tdt.LT.11368.d0) THEN
             rmsa=3.d0
             rmsd=3.d0
             ! From 1890 to 1950
          ELSE IF(tdt.LT.33282.d0) THEN
             rmsa=2.d0
             rmsd=2.d0
             ! After 1950
          ELSE
             rmsa=1.5d0
             rmsd=1.5d0
          ENDIF
          !Occultations
       ELSEIF(obs%tech.eq.'E')THEN
          rmsa=0.2d0
          rmsd=0.2d0
          !Hypparcos
       ELSEIF(obs%tech.eq.'H')THEN
          rmsa=0.4d0
          rmsd=0.4d0
          !Transit circle
       ELSEIF(obs%tech.eq.'T')THEN
          rmsa=0.5d0
          rmsd=0.5d0
          !Encoder
       ELSEIF(obs%tech.eq.'e')THEN
          rmsa=0.75d0
          rmsd=0.75d0
          !Micrometer
       ELSEIF(obs%tech.eq.'M')THEN
          rmsa=3.0d0
          rmsd=3.0d0
          !CCD - normal place - roving - satellites
       ELSEIF(obs%tech.eq.'c'.or.obs%tech.eq.'C'.or.obs%tech.eq.'n'.or.obs%tech.eq.'V'.or.obs%tech.eq.'S')THEN
          rmsa=1.0d0
          rmsd=1.0d0
          !Discovery observations - normally duplicated or very poor
       ELSEIF(obs%tech.eq.'x'.or.obs%tech.eq.'X')THEN
          rmsa=999.d0
          rmsd=999.d0
       ELSE
          WRITE(*,*)' obs. type unknown', obs%type,' ', obs%tech,' ',obs%obscod_s
          STOP 
       ENDIF
       ! use rules to overwrite default, but only for CCD observations 
       IF(obs%tech.eq.'c'.or.obs%tech.eq.'C'.or.obs%tech.eq.'n'.or.obs%tech.eq.'V'.or.obs%tech.eq.'S')THEN
          obscat=observ//obs%catcodmpc(2:2) ! four character string for a rule
          CALL bin_search(obscat,nrules,obscatcod,indrules,indobs)
          ! use rms appropriate for observatory, star catalog (in arcsec)
          IF(indobs.gt.0.and.indobs.le.nrules)THEN
             rmsa=rmsbin(1,indobs)
             rmsd=rmsbin(2,indobs)
          ELSE
             obscat=observ//'*'
             CALL bin_search(obscat,nrules,obscatcod,indrules,indobs)
             ! use rms appropriate for observatory, any catalog
             IF(indobs.gt.0.and.indobs.le.nrules)THEN
                rmsa=rmsbin(1,indobs)
                rmsd=rmsbin(2,indobs)
             ELSE
                obscat='ALL'//obs%catcodmpc(2:2)
                CALL bin_search(obscat,nrules,obscatcod,indrules,indobs)
                ! use rms appropriate for catalog, any station
                IF(indobs.gt.0.and.indobs.le.nrules)THEN
                   rmsa=rmsbin(1,indobs)
                   rmsd=rmsbin(2,indobs)
                ENDIF
             ENDIF
          ENDIF
       ENDIF
       !####################!
       ! ERROR MODEL CBM10  !
       !####################!
    ELSEIF(error_model.eq.'cbm10')THEN
       IF(biasflag)THEN
          ! use rms appropriate for observatory, whatever the star catalog
          CALL bin_search(observ,nstarms,obscodrms,indstarms,indobs)
          IF(indobs.gt.0.and.indobs.le.nstarms)THEN
             rmsa=rmsstabin(1,indobs)
             rmsd=rmsstabin(2,indobs)
             ermrms=.true.
          ELSE
             ! use rms appropriate for observatory, star catalog (in arcsec)
             obscat=observ//obs%catcodmpc(2:2) ! four character string for a rule
             CALL bin_search(obscat,nrules,obscatcod,indrules,indobs)
             IF(indobs.gt.0.and.indobs.le.nrules)THEN
                rmsa=rmsbin(1,indobs)
                rmsd=rmsbin(2,indobs)
                ermrms=.true.
             ELSE
                obscat='ALL'//obs%catcodmpc(2:2) ! rule for all other obs with given catalog
                CALL bin_search(obscat,nrules,obscatcod,indrules,indobs)
                IF(indobs.gt.0.and.indobs.le.nrules)THEN
                   rmsa=rmsbin(1,indobs)
                   rmsd=rmsbin(2,indobs)
                   ermrms=.true.
                ELSE
                   !              WRITE(iwarou,'(A,A,A,A,A,F11.5)')'astrow_bias: case without weights ',obscat, &
                   !    &     ' Observatory: ',observ, ' Time UTC: ',obs%time_utc
                   numwar=numwar+1
                   ermrms=.false.
                   ! Default value of RMS (time dependent, in arcsec)
                   rmsa=1.5d0
                   rmsd=1.5d0
                ENDIF
             ENDIF
          ENDIF
!!! for non-CCD default RMS
          !!IF(obs%tech.ne.'c'.or.obs%tech.ne.'C')THEN
          !!  ! Default value of RMS (time dependent, in arcsec)
          !!  rmsa=1.5d0
          !!  rmsd=1.5d0
          !!ENDIF
       ELSEIF(obs%tech.eq.'c'.or.obs%tech.eq.'C'.OR.&
!!! DF when the catalog is known, use proper weights
            obs%catcodmpc.EQ.' h'.OR.obs%catcodmpc.EQ.' i'&
            .OR.obs%catcodmpc.EQ.' j'.OR.obs%catcodmpc.EQ.' k'&
            .OR.obs%catcodmpc.EQ.' l'.OR.obs%catcodmpc.EQ.' m'&
            .OR.obs%catcodmpc.EQ.' n'.OR.obs%catcodmpc.EQ.' p'&
            .OR.obs%catcodmpc.EQ.' v'.OR.obs%catcodmpc.EQ.' w'&
            .OR.obs%catcodmpc.EQ.' x'.OR.obs%catcodmpc.EQ.' z')THEN
          ! only CCD observations have special model
          ! is there a rule for this observatory?
          CALL bin_search(observ,nstarms,obscodrms,indstarms,indobs)
          IF(indobs.gt.0.and.indobs.le.nrules)THEN
             rmsa=rmsstabin(1,indobs)
             rmsd=rmsstabin(2,indobs)
             ermrms=.true.
          ELSE
             ! use rms appropriate for observatory, star catalog (in arcsec)
             obscat=observ//obs%catcodmpc(2:2) ! four character string for a rule
             CALL bin_search(obscat,nrules,obscatcod,indrules,indobs)
             IF(indobs.gt.0.and.indobs.le.nrules)THEN
                rmsa=rmsbin(1,indobs)
                rmsd=rmsbin(2,indobs)
                ermrms=.true.
             ELSE
                obscat='ALL'//obs%catcodmpc(2:2) ! rule for all other obs with given catalog
                CALL bin_search(obscat,nrules,obscatcod,indrules,indobs)
                IF(indobs.gt.0.and.indobs.le.nrules)THEN
                   rmsa=rmsbin(1,indobs)
                   rmsd=rmsbin(2,indobs)
                   ermrms=.true.
                ELSE
                   !              WRITE(iwarou,'(A,A,A,A,A,F11.5)')'astrow_bias: case without weights ',obscat,  &
                   !&         ' Observatory: ',observ, ' Time UTC: ',obs%time_utc
                   numwar=numwar+1
                   ermrms=.false.
                   ! Default value of RMS (time dependent, in arcsec)
                   rmsa=1.5d0
                   rmsd=1.5d0
                ENDIF
             ENDIF
          ENDIF
       ELSE
          ermrms=.false.
          ! non-CCD observations
!!! DF transit type
          IF(obs%tech.eq.'t'.or.obs%tech.eq.'T')THEN
             rmsa=0.5d0
             rmsd=0.5d0
             ! Default value of RMS (time dependent, in arcsec)
             ! Before 1890
          ELSEIF(tdt.LT.11368.d0) THEN
             rmsa=3.d0
             rmsd=3.d0
             ! From 1890 to 1950
          ELSE IF(tdt.LT.33282.d0) THEN
             rmsa=2.d0
             rmsd=2.d0
             ! After 1950
          ELSE
             rmsa=1.2d0
             rmsd=1.2d0
          ENDIF
       ENDIF
    ENDIF
    ! change rms value only if it is not forced
330 IF(.not.obsw%force_w(1))THEN 
       obsw%rms_coord(1)=MAX(obs%acc_coord(1),rmsa*radsec/cos(obs%coord(2)))
    ENDIF
    IF(.not.obsw%force_w(2))THEN
       obsw%rms_coord(2)=MAX(obs%acc_coord(2),rmsd*radsec)
    ENDIF
  END SUBROUTINE astrow_bias
  ! ====================================================================
  ! ASTROW_NOBIAS
  ! No bias from catalogs, assign RMS
  ! ====================================================================
  ! INPUT:    MPCTYP    -  Observation type (column 15 of MPC record)
  !           TDT       -  Time of observation (MJD, TDT)
  !           IDSTA     -  Observatory code 9numeric)
  !           ACCA      -  Accuracy of right ascension (rad)
  !           ACCD      -  Accuracy of declination (rad)
  !
  ! OUTPUT:   RMSA      -  A-priori RMS of right ascension (rad)
  !           RMSD      -  A-priori RMS of declination (rad)
  !           BIASA     -  Bias in  right ascension (rad)
  !           BIASD     -  Bias in  declination (rad)
  !           DECSTEP   -  Decision steps of RMS assignment
  !
  SUBROUTINE astrow_nobias(error_model,obs,obsw)
    TYPE(ast_obs) , INTENT(IN)    :: obs         ! observation
    CHARACTER*(*) , INTENT(IN)    :: error_model ! code for errmod
    TYPE(ast_wbsr), INTENT(INOUT) :: obsw        ! assigned weight

    ! ==============END INTERFACE==================================
    !mpctyp,tdt,idsta,acca,accd,        &
    !     &        rmsa,rmsd,biasa,biasd,decstep)
    DOUBLE PRECISION              :: tdt,acca,accd
    CHARACTER*(1)                 :: mpctyp
    DOUBLE PRECISION              :: rmsa,rmsd,biasa,biasd, rmsmin
    LOGICAL                       :: error
    ! =============================================================
    ! no bias anyway
    obsw%bias_coord(1:2)=0.d0
    ! default, rough rule-of-thumb error model
    ! Default value of RMS (time dependent)
    tdt=obs%time_tdt
    ! Before 1890
    IF(tdt.LT.11368.d0) THEN
       rmsmin=3.d0
       ! From 1890 to 1950
    ELSE IF(tdt.LT.33282.d0) THEN
       rmsmin=2.d0
       ! After 1950
    ELSE
       rmsmin=1.d0
    END IF
    rmsmin=rmsmin*radsec
    acca=obs%acc_coord(1)
    accd=obs%acc_coord(2)
    rmsa=MAX(rmsmin/cos(obs%coord(2)),acca)
    rmsd=MAX(rmsmin,accd)
    obsw%rms_coord(1)=rmsa
    obsw%rms_coord(2)=rmsd
  END SUBROUTINE astrow_nobias
  !***********************************************************************
  ! 'magrms' returns the default magnitude rms based on the MPC obs string
  !***********************************************************************
  DOUBLE PRECISION FUNCTION magrms_new(magstr,tdt,idsta,typ)
    CHARACTER*6 magstr
    ! obs. type (from column 15 of .obs format), station code, time MJD
    CHARACTER*1 typ
    INTEGER idsta
    DOUBLE PRECISION tdt
    INTEGER ll,lench
    ! magnitude weighting for now is simple;
    ! should be based on digits given,color
    ll=lench(magstr)
    IF(ll.le.0)THEN
       magrms_new=-1.d0
       RETURN
    ENDIF
    IF(magstr(3:5).eq.'   ')THEN
       magrms_new=1.0
    ELSEIF(magstr(5:5).eq. ' ')THEN
       magrms_new=0.7
    ELSE
       magrms_new=0.5
    ENDIF
  END  FUNCTION magrms_new

  ! END SUB_MODULE observ_bias_weight

  ! ------------------------------------------------------ !
  ! SUBROUTINE convert_ades_fitobs                         !
  ! It fills the general data type used by differential    !
  ! corrections using the values read from the ADES file   !
  ! Finally, it updates ades_obs with sigmas and biases    !
  ! computed in this routine. These will be written in the !
  ! output file as part of the residuals group             !
  ! ------------------------------------------------------ !
  ! Author: Alessia Bertolucci, SpaceDyS 2020              !
  ! ------------------------------------------------------ !

  SUBROUTINE convert_ades_fitobs(error_model,observ,wei_bia)
    
    CHARACTER(LEN=*),   INTENT(IN)  :: error_model            ! error model file name 
    TYPE(ast_obs),      INTENT(OUT) :: observ(nobs_ades)      ! observed observables for diffcor routine
    TYPE(ast_wbsr),     INTENT(OUT) :: wei_bia(nobs_ades)     ! weights and biases for diffcor routine

    TYPE(ast_obs)       :: obs1(nobs_ades),obs_loc
    TYPE(ast_wbsr)      :: obsw_loc
    INTEGER             :: i,len_str,lenna,k,aux_i1,aux_i2,prec,j,jj,lendes
    CHARACTER(LEN=4)    :: aux_str1,aux_str2
    CHARACTER(LEN=5)    :: mag_string
    CHARACTER(LEN=7)    :: obsnum
    CHARACTER(LEN=30)   :: time_str

    ! Times: time_ut1utc,tdb,tdtj
    INTEGER             :: mjd_ut1utc,mjd_out,mjd_out_tdt     ! integer MJD date (1 day=86400.0 s)
    REAL(KIND=dkind)    :: sec_ut1utc,sec_out,sec_out_tdt     ! time from midnight, in seconds
    CHARACTER(LEN=3)    :: scale_ut1utc                       ! Time scale 
    INTEGER             :: year,month,day,hour,minute,mjdd
    REAL(KIND=dkind)    :: tjmd,second,secs,t2,day_frac
    REAL(KIND=dkind)    :: ra,dec,ra_ref,dec_ref,dist,pa,deltaRa,deltaDec,range,rrate
    REAL(KIND=dkind)    :: deltaRa_arcsec,deltaDec_arcsec,deltaRa_deg,deltaDec_deg
    REAL(KIND=dkind)    :: biasRA,biasDec,rms_radf
    REAL(KIND=dkind)    :: A,B,C
    REAL(KIND=dkind)    :: xsun(6)                ! Sun ECLJ geoc. state vector
    REAL(KIND=dkind)    :: obs_pos(3)             ! observer's geocentric position (au)
    REAL(KIND=dkind)    :: geo_pos(3)             ! roving observatory observations (observ%type='O' and observ%obscod_s='247')
                                                  ! geo_pos(1) = E longitude of observing site (rad)
                                                  ! geo_pos(2) = N latitude of observing site (rad)
                                                  ! geo_pos(3) = Altitude of observing site (m)

    REAL(KIND=dkind), DIMENSION(3,3)         :: rot                     ! rotation matrix
    REAL(KIND=dkind), PARAMETER              :: gapmax=8.d0/24.d0       ! 8 hours max gap
    REAL(KIND=dkind), DIMENSION(2,nobs_ades) :: min_rms
    REAL(KIND=dkind), DIMENSION(2,2)         :: DP_cov_mat,RD_cov_mat
    
    CHARACTER(LEN=80)   :: name_em_loc                          ! local name of error model
    CHARACTER(LEN=16)   :: stname
    CHARACTER(LEN=11)   :: objdes
    CHARACTER(LEN=6)    :: notes
    CHARACTER(LEN=2)    :: prog
    LOGICAL             :: compute_errormodel_rms,compute_errormodel_bias    ! if TRUE we use the error model to compute either rms or bias
    LOGICAL             :: compute_all_errormodel,use_errmod(nobs_ades)      ! if TRUE we use the error model to compute weights and biases 
    LOGICAL             :: init          ! if true, reinitialize wei_bia, else information different from weight/bias is preserved
    INTEGER             :: new_prec(2),len_notes,len_prog,l
    INTEGER             :: ind(nobs_ades),nob(nobs_ades),nob1(nobs_ades) !sorting by time, number obs in batch
    REAL(KIND=dkind)    :: bias_rad,ra_acc,dec_acc,rms_rad,cov_radec,detg2
    REAL(KIND=dkind)    :: corr_radec(nobs_ades)
    REAL(KIND=dkind)    :: rmsRA,rmsDec,factor,varRA,varDec,corr2,eigval(2)
    CHARACTER(LEN=100)  :: deltara_arcsec_str,deltadec_arcsec_str,deltara_deg_str,deltadec_deg_str
    CHARACTER(LEN=15)   :: day_fstr
    
    INTEGER,          PARAMETER  :: precTime_values(6) = (/1,10,100,1000,10000,100000/)
    CHARACTER(LEN=3), PARAMETER   :: station_rules(16) = (/"F51","F52","703","G96","I52","V06","V00", &
         & "T05","T08","I41","C51","H01","568","T09","T12","T14"/)
    
    REAL(KIND=dkind), EXTERNAL :: tjm1
    
    ! Initialization
    name_em_loc = error_model          
    CALL rmsp(name_em_loc,lenna)
    observ  = undefined_ast_obs
    wei_bia = undefined_ast_wbsr
    compute_all_errormodel = .FALSE.
    use_errmod(1:nobs_ades) = .FALSE.
    min_rms(1:2,1:nobs_ades) = 0.d0
    
    ! if error model has changed, we use it to re-compute weights and biases of all observations from ADES file
    IF (ades_error_model .NE. "" .AND. ades_error_model .NE. error_model) compute_all_errormodel = .TRUE.

    ! loop over ADES observations
    DO i=1, nobs_ades
       
       !----- objdes (observ) -----!

       IF (ades_flg(i)%permID) THEN
          observ(i)%objdes = trim(adjustl(ades_obs(i)%permID))
       ELSEIF (ades_flg(i)%provID) THEN
          objdes = ades_obs(i)%provID
          CALL rmsp(objdes,lendes)
          observ(i)%objdes = objdes
       ELSE
          objdes = ades_obs(i)%trkSub
          CALL rmsp(objdes,lendes)
          observ(i)%objdes = objdes
       ENDIF
       
       !----- type (observ) -----!
              
       ! if optical obs
       IF (.NOT. ades_flg(i)%rad) THEN
          ! if location group is present and not roving then satellite/space
          IF (ades_flg(i)%sys .AND. (trim(adjustl(ades_obs(i)%stn)) .NE. '247')) THEN
             observ(i)%type = 'S'
             ! otherwise optical for ground-based or roving
          ELSE
             observ(i)%type = 'O'
          ENDIF
          ! if radar obs
       ELSE
          ! if delay then radar range
          IF (ades_flg(i)%delay) THEN
             observ(i)%type = 'R'
             ! if doppler then radar range-rate
          ELSE
             observ(i)%type = 'V'
          ENDIF
       ENDIF
       
       !----- tech (observ) -----!

       IF (.NOT. ades_flg(i)%rad) THEN
          IF (trim(adjustl(ades_obs(i)%deprecated)) .EQ. 'X') THEN
             observ(i)%tech = "X"
          ELSE
             ! if occultation obs
             IF (ades_flg(i)%occ) THEN
                observ(i)%tech = "E"
                ! if roving
             ELSEIF (trim(adjustl(ades_obs(i)%stn)) .EQ. '247') THEN
                observ(i)%tech = "V"                                                    
                ! if satellite 
             ELSEIF (ades_flg(i)%sys) THEN
                observ(i)%tech = "S" 
                ! otherwise ground-based
             ELSE
                CALL techn_mpc(ades_obs(i)%mode,ades_obs(i)%notes,observ(i)%tech)             
             ENDIF
          ENDIF
       ELSE
          IF (ades_flg(i)%com .AND. ades_obs(i)%com == 0) THEN
             observ(i)%tech = "s"
          ELSE  ! if %com is not present, we assume it true
             observ(i)%tech = "c"
          ENDIF
       ENDIF

       !----- note (observ) -----!
      
       IF (ades_flg(i)%notes) THEN
          notes = ades_obs(i)%notes
          CALL rmsp(notes,len_notes)
       ENDIF
       IF (ades_flg(i)%prog) THEN
          prog = ades_obs(i)%prog
          CALL rmsp(prog,len_prog)
       ENDIF
       IF (ades_flg(i)%notes .AND. len_notes == 1 .AND. notes .NE. "n") THEN
          observ(i)%note = notes
       ELSEIF (ades_flg(i)%prog .AND. len_prog == 1) THEN
          observ(i)%note = prog
       ELSEIF (ades_flg(i)%prog .AND. prog == "pi") THEN
          observ(i)%note = "|"
       ENDIF

       !----- catcodmpc (observ) -----!
           
       ! astrometric catalog code 
       IF (.NOT. ades_flg(i)%rad) THEN
          CALL catcod_mpc(ades_obs(i)%astCat,observ(i)%catcodmpc)
          ! if astCat is not specified by the observer, use the error model to get weights and biases 
          IF (trim(adjustl(ades_obs(i)%astCat)) .EQ. "" .OR. &
               & trim(adjustl(ades_obs(i)%astCat)) .EQ. "UNK") use_errmod(i) = .TRUE.
       ENDIF
       
       !----- obscod (observ) -----!

       ! if optical obs then stn
       IF (.NOT. ades_flg(i)%rad) THEN
          ! we take the first tree chars of stn (which is up to 4 chars)
          aux_str1 = trim(adjustl(ades_obs(i)%stn))
          observ(i)%obscod_s = aux_str1(1:3)
          CALL statcode(observ(i)%obscod_s,observ(i)%obscod_i)
       ELSE  ! if radar obs then trx+" "+rcv
          ! obscod_s is at most 7 chars, so we restrict rcv and trx to 3 chars each
          aux_str1 = trim(adjustl(ades_obs(i)%trx))
          aux_str2 = trim(adjustl(ades_obs(i)%rcv))
          observ(i)%obscod_s = aux_str1(1:3)//' '//aux_str2(1:3)
          CALL statcode(aux_str1(1:3),aux_i1)
          CALL statcode(aux_str2(1:3),aux_i2)
          observ(i)%obscod_i = aux_i1*10000 + aux_i2                                             
       ENDIF
       ! If observation is optical and has not dedicated rule in neocp.rules,
       ! the weight is assigned according to error model
       IF (.NOT.ades_flg(i)%rad .AND. ons_name .AND. &
            & ALL(station_rules .NE. trim(adjustl(ades_obs(i)%stn)))) &
            & use_errmod(i) = .TRUE.

       !--- time_tdt,time_utc,acc_time (observ) ---!
       
       time_str = trim(adjustl(ades_obs(i)%obsTime))
       CALL rmsp(time_str,len_str)
       READ(time_str(1:4),*)year
       READ(time_str(6:7),*)month
       READ(time_str(9:10),*)day
       READ(time_str(12:13),*)hour
       READ(time_str(15:16),*)minute
       READ(time_str(18:len_str-1),*)second
       ! computation of Modified Julian Date
       tjmd = tjm1(day,month,year,0.d0)
       mjd_ut1utc = FLOOR(tjmd)
       sec_ut1utc = second+minute*60.d0+hour*3600.d0

       ! correction for bias in time        
       secs = sec_ut1utc-ades_obs(i)%biasTime  ! TO BE CHECKED after real data processing
       mjdd = FLOOR(secs/86400.d0)
       mjd_ut1utc = mjd_ut1utc + mjdd
       sec_ut1utc = secs-mjdd*86400.d0  !seconds within a day

       IF (year.LT.1972) THEN
          scale_ut1utc='UT1'
       ELSE
          scale_ut1utc='UTC'
       END IF
       mjd_out = mjd_ut1utc
       sec_out = sec_ut1utc
       day_frac = sec_out/86400.d0
       observ(i)%time_utc = mjd_out + day_frac
       CALL cnvtim(mjd_out,sec_out,scale_ut1utc,mjd_out_tdt,sec_out_tdt,'TDT')
       observ(i)%time_tdt = mjd_out_tdt + sec_out_tdt/86400.d0

       ! fill acc_time
       IF (.NOT. ades_flg(i)%rad) THEN
          ! get precision of seconds to day conversion (1/86400) 
          WRITE(day_fstr,'(F9.6)') day_frac     ! max precision is 10**-6
          CALL rmsp(day_fstr,l)
          IF (INDEX(day_fstr,'.') .NE. 0) THEN
             DO WHILE (day_fstr(l:l) .EQ. '0')  ! delete trailing zeros
                IF (day_fstr(l-1:l-1) .EQ. '.') EXIT  ! in case of only zeros, leave one zero after decimal point
                day_fstr = day_fstr(1:l-1)
                l = l-1
             ENDDO
          ENDIF
          IF (INDEX(day_fstr,'.') .GT. 0) THEN
             prec = get_prec(day_fstr)
          ELSE
             prec = 0
          ENDIF
          !prec = 5
          !IF (INDEX(time_str,'.') .NE. 0) prec = prec + len_str - INDEX(time_str,'.') - 1
          observ(i)%acc_time = 10.d0**(-prec*1.d0)
          ! compare with precTime if present
          IF (ades_flg(i)%precTime .AND. (.NOT. ALL(precTime_values .NE. ades_obs(i)%precTime))) THEN
             observ(i)%acc_time = MAX(observ(i)%acc_time,ades_obs(i)%precTime/1.d6)
          ENDIF
       ELSE
          ! for radar observations always 10d-10
          observ(i)%acc_time = 10d-10
       ENDIF

       !----- coord (observ) -----!
       
       ! if all optical obs then ra, dec
       IF (.NOT. ades_flg(i)%rad) THEN
          ! if optical we already have ra and dec
          IF (ades_flg(i)%opt) THEN
             ! convert to radians
             ra = ades_obs(i)%ra*radeg
             dec = ades_obs(i)%dec*radeg
             ! if occultation we need to compute ra and dec
          ELSEIF (ades_flg(i)%occ) THEN
             ! find ra and dec of the origin (occulted star)
             ! convert to radians
             ra_ref = ades_obs(i)%raStar*radeg
             dec_ref = ades_obs(i)%decStar*radeg
             ! update precision for output format
             ades_prec(i)%ra = ades_prec(i)%raStar
             ades_prec(i)%dec = ades_prec(i)%decStar
             ! if we have dist and pa convert them in deltaRA and deltaDec
             IF (ades_flg(i)%dist) THEN
                ! convert to radians
                dist = ades_obs(i)%dist*radsec
                pa = ades_obs(i)%pa*radeg
                ! convert dist and pa in deltaRA and deltaDec following equations 1-2 in
                ! http://www.physicspages.com/2015/04/11/angular-distances-on-the-celestial-sphere/
                ! Carroll, Bradley W. & Ostlie, Dale A. (2007), An Introduction to Modern Astrophysics
                A = SIN(dist)*SIN(pa)
                B = SIN(dec_ref)*COS(dist)
                C = COS(dec_ref)*COS(pa)
                ! deltaRA and deltaDec in radians
                deltaRA = ASIN(A / SQRT(1 - (B+C)**2))
                deltaDec = ASIN(B+C) - dec_ref
                ! convert deltaRA and deltaDec to arcsec in order to compute their precisions for output format
                deltaRa_arcsec = deltaRA*secrad
                deltaDec_arcsec = deltaDec*secrad
                WRITE(deltara_arcsec_str,*) deltaRa_arcsec
                WRITE(deltadec_arcsec_str,*) deltaDec_arcsec
                ades_prec(i)%deltaRA = get_prec(deltara_arcsec_str)
                ades_prec(i)%deltaDec = get_prec(deltadec_arcsec_str)

                ! otherwise we already have deltaRA and deltaDec
             ELSE
                ! convert to radians
                deltaRA = ades_obs(i)%deltaRA*radsec
                deltaDec = ades_obs(i)%deltaDec*radsec
             ENDIF
             ! update precision of ra and dec for output format
             ! conversion of deltaRa and deltaDec from radians to degrees
             deltaRa_deg = deltaRA*degrad
             deltaDec_deg = deltaDec*degrad
             ! precision of deltaRa and deltaDec converted to degrees
             WRITE(deltara_deg_str,*) deltaRa_deg
             WRITE(deltadec_deg_str,*) deltaDec_deg
             new_prec(1) = get_prec(deltara_deg_str)
             new_prec(2) = get_prec(deltadec_deg_str)
             IF (new_prec(1) .GT. ades_prec(i)%raStar) ades_prec(i)%ra = new_prec(1)
             IF (new_prec(2) .GT. ades_prec(i)%decStar) ades_prec(i)%dec = new_prec(2) 
             ! find ra and dec by adding the corresponding deltas to the ref values
             dec = dec_ref + deltaDec
             ! scale deltaRA and add it to ra_ref
             deltaRA = deltaRA/COS(dec)
             ra = ra_ref + deltaRA
          ENDIF
          ! check values of ra and dec after conversion(s)
          IF (ra/radeg .LT. 0.d0 .OR. ra/radeg .GT. 360.d0) THEN
             WRITE(obsnum,'(i7)')i
             WRITE(ierrou,*) 'convert_ades_fitobs: Computed RA value not in [0,360] degrees &
                  &for input ADES file '//trim(ades_filename)//' obs number: '//trim(adjustl(obsnum))
             numerr=numerr+1 
          ENDIF
          IF (dec/radeg .LT. -90.d0 .OR. dec/radeg .GT. 90.d0) THEN
             WRITE(obsnum,'(i7)')i
             WRITE(ierrou,*) 'convert_ades_fitobs: Computed Dec value not in [-90,90] degrees &
                  &for input ADES file '//trim(ades_filename)//' obs number: '//trim(adjustl(obsnum))
             numerr=numerr+1 
          ENDIF
          ! store coordinates
          observ(i)%coord(1) = ra
          observ(i)%coord(2) = dec
          ! if radar obs then range (au) or range-rate (au/d)
       ELSE
          ! if delay, then range is the distance = (vlight/8.64d4)*delay/2 (half of back-and-forth travel time)
          IF (ades_flg(i)%delay) THEN
             range = 0.5d0*ades_obs(i)%delay*vlight/8.64d4
             observ(i)%coord(1) = range
             ! if doppler, then velocity = -doppler shift * carrier wavelength (again divided 2)
          ELSE
             rrate = -0.5d0*ades_obs(i)%doppler*(1.d-6*vlight/ades_obs(i)%frq)
             observ(i)%coord(2) = rrate
          ENDIF
       ENDIF
       
       ! --- acc_coord (observ) --- !
       ! fill observ(i)%acc_coord in case of optical observations (case of radar observations below)
       IF (.NOT. ades_flg(i)%rad) THEN
          new_prec(1) = 0
          new_prec(2) = 0
          CALL RADec_accuracy(i,dec,observ(i)%acc_coord(1),observ(i)%acc_coord(2),new_prec) 
       ENDIF
 
       ! --- magnitude stuff for observ --- !
 
       IF (ades_flg(i)%mag) THEN
          ! fill observ 
          ! mag: apparent magnitude
          observ(i)%mag = ades_obs(i)%mag
          ! mag_band
          observ(i)%mag_band = ades_obs(i)%band
          ! mag_str
          IF (ades_obs(i)%mag .LE. 0.0d0) THEN
             mag_string = "     "
             ! acc_mag
             observ(i)%acc_mag = 99.d9
          ELSE
             IF (ades_prec(i)%mag .LE. 0) THEN
                WRITE(mag_string,"(I2)") NINT(ades_obs(i)%mag)
             ELSEIF (ades_prec(i)%mag .EQ. 1) THEN
                WRITE(mag_string,"(F4.1)") ades_obs(i)%mag
                mag_string = adjustl(mag_string)//" "
             ELSEIF (ades_prec(i)%mag .GE. 2) THEN
                WRITE(mag_string,"(F5.2)") ades_obs(i)%mag   
             ENDIF
             ! acc_mag
             observ(i)%acc_mag = 10.d0**(-ades_prec(i)%mag*1.d0)
          ENDIF
          observ(i)%mag_str = mag_string//ades_obs(i)%band
          ! mag_def
          IF (trim(adjustl(mag_string)) .NE. "") THEN
             observ(i)%mag_def = .TRUE.
          ELSE
             observ(i)%mag_def = .FALSE.
          ENDIF

       ELSEIF (.NOT. ades_flg(i)%rad) THEN                       
          observ(i)%mag_def=.FALSE.
          observ(i)%mag=0.d0
          observ(i)%acc_mag=99.d9
       ENDIF

       !--- obspos/obsvel, par_unit (observ) ---!
       
       ! for optical obs: ECLJ2000 geocentric position of the observer
       ! ground-based, not roving
       IF (observ(i)%type .EQ. 'O' .AND. observ(i)%obscod_s .NE. '247') THEN
          ! find positions
          ! get state vector from internal code
          CALL observer_position(observ(i)%time_tdt,observ(i)%obspos,observ(i)%obsvel,OBSCODE=observ(i)%obscod_i)

          ! satellite
       ELSEIF (observ(i)%type .EQ. 'S') THEN
          ! find positions
          ! we have already positions in au 
          IF (ades_obs(i)%sys .EQ. "ICRF_AU") THEN
             obs_pos = ades_obs(i)%pos
             ! convert positions in au from km (sys = ICRF_KM)
          ELSE
             obs_pos = ades_obs(i)%pos/aukm
          ENDIF
          ! apply rotation to Ecliptic J2000.0 reference system
          obs_pos = MATMUL(roteqec,obs_pos)

          ! translation to geocentric position (unless ctr = 399, i.e. we are already geocentric)
          IF (ades_obs(i)%ctr .EQ. 399) THEN
             observ(i)%obspos = obs_pos
             observ(i)%obsvel = 0.d0
          ELSE
             CALL sunmoon_car(observ(i)%time_tdt,.FALSE.,xsun)
             ! translate observer position to geocentric
             ! xsun is an array containing geocentric cartesian coordinates ECLJ2000 of Sun. The units are au and au/day.
             observ(i)%obspos = obs_pos + xsun(1:3)
             observ(i)%obsvel = 0.d0
          ENDIF

          ! fill par_unit
          IF (ades_obs(i)%sys .EQ. "ICRF_AU") THEN
             observ(i)%par_unit = "2"
          ELSE
             observ(i)%par_unit = "1"
          ENDIF

          ! roving
       ELSEIF (observ(i)%obscod_s .EQ. '247') THEN
          IF(ASSOCIATED(observ(i)%geopos)) DEALLOCATE(observ(i)%geopos)
          ALLOCATE(observ(i)%geopos(3))
          IF (trim(adjustl(ades_obs(i)%sys)) .EQ. "WGS84") THEN
             geo_pos(1) = ades_obs(i)%pos(1)*radeg
             geo_pos(2) = ades_obs(i)%pos(2)*radeg
             geo_pos(3) = ades_obs(i)%pos(1)
             observ(i)%geopos = geo_pos
             ! convert to cartesian coordinates (in m)
             CALL geodetic_to_cartesian(geo_pos(1),geo_pos(2),geo_pos(3),obs_pos)
             ! convert in au
             obs_pos=obs_pos/(aukm*1d3)
             ! compute position and velocity in Ecliptic J2000
             CALL observer_position(observ(i)%time_tdt,observ(i)%obspos,observ(i)%obsvel,BFPOS=obs_pos)

          ELSEIF (trim(adjustl(ades_obs(i)%sys)) .EQ. "ITRF") THEN
             ! convert from cylindrical to cartesian coordinates (in km)
             obs_pos(1) = ades_obs(i)%pos(2)*COS(ades_obs(i)%pos(1)*radeg)
             obs_pos(2) = ades_obs(i)%pos(2)*SIN(ades_obs(i)%pos(1)*radeg)
             obs_pos(3) = ades_obs(i)%pos(3)
             ! convert in au
             obs_pos=obs_pos/(aukm)
             ! compute position and velocity in Ecliptic J2000
             CALL observer_position(observ(i)%time_tdt,observ(i)%obspos,observ(i)%obsvel,BFPOS=obs_pos)
             ! fill observ(i)%geopos
             observ(i)%geopos(1) = ades_obs(i)%pos(1)*radeg
             observ(i)%geopos(2) = ATAN(ades_obs(i)%pos(3)/ades_obs(i)%pos(2))
             observ(i)%geopos(3) = SQRT(ades_obs(i)%pos(2)**2+ades_obs(i)%pos(3)**2)*1d3

          ELSEIF (trim(adjustl(ades_obs(i)%sys)) .EQ. "IAU") THEN

             ! convert to cartesian coordinates (in km)
             CALL geoIAU_to_cart(ades_obs(i)%pos,obs_pos)
             ! get rotation to ECLJ2000 reference system
             CALL rot_IAU_to_ECLJ2000(observ(i)%time_tdt,rot)
             ! apply rotation to ECLJ2000 reference system
             obs_pos = MATMUL(rot,obs_pos)
             ! get equatorial coordinates (used in call to observer_position for computing %obsvel)
             obs_pos = MATMUL(roteceq,obs_pos)
             ! convert in au
             obs_pos=obs_pos/(aukm)
             ! compute position and velocity in Ecliptic J2000
             CALL observer_position(observ(i)%time_tdt,observ(i)%obspos,observ(i)%obsvel,BFPOS=obs_pos)

             ! fill observ(i)%geopos
             observ(i)%geopos(1) = ades_obs(i)%pos(1)*radeg
             observ(i)%geopos(2) = ades_obs(i)%pos(2)*radeg
             observ(i)%geopos(3) = ades_obs(i)%pos(3)

          ENDIF

          ! radar observations
          ! obspos = body-fixed position of transmitting station
          ! obsvel = body-fixed position of receiving station
       ELSE IF (observ(i)%type .EQ. 'R' .OR. observ(i)%type .EQ. 'V') THEN
          CALL statcode(observ(i)%obscod_s(1:3),aux_i1)
          CALL obscoo(aux_i1,observ(i)%obspos,stname)    
          CALL statcode(observ(i)%obscod_s(5:7),aux_i2)
          CALL obscoo(aux_i2,observ(i)%obsvel,stname)
       ENDIF

       ! ----- sel_coord (wei_bia) ----- !
       
       IF (.NOT. ades_flg(i)%rad .AND. ades_fit(i)%ast_flag .NE. "") THEN
          READ(ades_fit(i)%ast_flag,'(i1)') wei_bia(i)%sel_coord
       ELSEIF (ades_flg(i)%rad .AND.  ades_fit(i)%radar_flag .NE. "") THEN
          READ(ades_fit(i)%radar_flag,'(i1)') wei_bia(i)%sel_coord
       ELSE
          ! initialize to 1 as default
          ! we neglect wei_bia(i)%sel_coord = 2 value
          wei_bia(i)%sel_coord = 1
          ! optical obs
          IF (.NOT. ades_flg(i)%rad) THEN
             ! if selAst forced rejected or automatic rejected 
             IF (ades_flg(i)%selAst .AND. (ades_obs(i)%selAst .EQ. "d" .OR. ades_obs(i)%selAst .EQ. "D")) wei_bia(i)%sel_coord = 0
             ! for radar obs
          ELSE
             ! selDelay
             IF (ades_flg(i)%selDelay) THEN
                ! if forced rejected or automatic rejected
                IF (ades_obs(i)%selDelay .EQ. "d" .OR. ades_obs(i)%selDelay .EQ. "D") wei_bia(i)%sel_coord = 0
                ! selDoppler
             ELSEIF (ades_flg(i)%selDoppler) THEN
                ! if forced rejected or automatic rejected
                IF (ades_obs(i)%selDoppler .EQ. "d" .OR. ades_obs(i)%selDoppler .EQ. "D") wei_bia(i)%sel_coord = 0
             ENDIF
          ENDIF
       ENDIF
       
       !--- force_w (wei_bia) ---!

       wei_bia(i)%force_w(1) = ades_fit(i)%force_w_flags(1)
       wei_bia(i)%force_w(2) = ades_fit(i)%force_w_flags(2)
       
       !--- chi (wei_bia) ---!

       wei_bia(i)%chi = ades_fit(i)%chi

       ! ----- sel_mag (wei_bia) ----- !
              
       ! we also set sel_mag here, given that selPhot might be neglected later if we use error model
       IF (.NOT. ades_flg(i)%rad) THEN
          IF (ades_fit(i)%mag_flag .NE. "") THEN
             READ(ades_fit(i)%mag_flag,'(i1)') wei_bia(i)%sel_mag
          ELSE
             ! initialize to 0 as default
             wei_bia(i)%sel_mag = 0
             IF (ades_flg(i)%selPhot) THEN 
                ! if selPhot forced rejected or automatic rejected
                IF (ades_obs(i)%selPhot .EQ. "a" .OR. ades_obs(i)%selPhot .EQ. "A") wei_bia(i)%sel_mag = 1
             ENDIF
          ENDIF
       ELSE
          wei_bia(i)%sel_mag = 0
       ENDIF
       
       ! ---- bias and weights (wei_bia) ----- !
       
       ! initialize flags for computation of bias and weights via error model as default
       compute_errormodel_bias = .TRUE.
       compute_errormodel_rms = .TRUE.

       ! erase non-forced weigths
       IF (.NOT. ades_flg(i)%rad) THEN
          IF(.NOT.ades_fit(i)%force_w_flags(1)) ades_obs(i)%sigRA = 0.d0
          IF(.NOT.ades_fit(i)%force_w_flags(2)) ades_obs(i)%sigDec = 0.d0
          IF(.NOT.(ALL(ades_fit(i)%force_w_flags))) ades_obs(i)%sigCorr = 0.d0
       ELSE
          IF(ades_flg(i)%delay .AND. .NOT. ades_fit(i)%force_w_flags(1)) ades_obs(i)%sigDelay = 0.d0
          IF(ades_flg(i)%doppler .AND. .NOT. ades_fit(i)%force_w_flags(2)) ades_obs(i)%sigDoppler = 0.d0
       ENDIF

       IF (.NOT. ades_flg(i)%rad) THEN

          ! store forced weights: sigRA and sigDec, if present in input ADES file
          IF (ades_flg(i)%sigRA .AND. ades_obs(i)%sigRA .GT. 0.d0 .OR. ades_obs(i)%sigDec .GT. 0.d0) THEN
             
             ! convert to radians
             IF (ades_obs(i)%sigRA .GT. 0.d0) THEN
                rmsRA = ades_obs(i)%sigRA/(secrad*COS(dec))
             ENDIF
             IF (ades_obs(i)%sigDec .GT. 0.d0) THEN 
                rmsDec = ades_obs(i)%sigDec/secrad
             ENDIF

             ! store correlation
             ! sigCorr is an optional element
             IF (ades_flg(i)%sigCorr .AND. ades_obs(i)%sigRA .GT. 0.d0 .AND. ades_obs(i)%sigDec .GT. 0.d0 &
                  & .AND.(ades_obs(i)%sigCorr .GE. -1.d0 .AND. ades_obs(i)%sigCorr .LE. 1.d0)) THEN  
                corr_radec(i) = ades_obs(i)%sigCorr
             ELSE
                corr_radec(i) = 0.d0
             ENDIF
          
             ! set rms computation flag to FALSE
             IF (ades_obs(i)%sigRA .GT. 0.d0 .AND. ades_obs(i)%sigDec .GT. 0.d0) compute_errormodel_rms = .FALSE.

             ! store (provisionally) weights
             wei_bia(i)%rms_coord(1) = rmsRA
             wei_bia(i)%rms_coord(2) = rmsDec
             wei_bia(i)%correl = corr_radec(i)
          ENDIF

          IF (.NOT. compute_all_errormodel) THEN
             
             ! *** BIAS *** for optical obs
             ! we store the biases directly from input ADES file if present
             IF (.NOT.ades_flg(i)%rad.AND.ades_flg(i)%biasRA .AND. (ades_obs(i)%biasRA.NE.0.d0 .OR. &
                  & (ades_obs(i)%biasRA .EQ. 0.d0 .AND. (ades_obs(i)%astCat .EQ. "Gaia1" .OR. &
                  & ades_obs(i)%astCat .EQ. "Gaia2")))) THEN
                
                ! convert to radians
                biasRA = ades_obs(i)%biasRA/(secrad*COS(dec))
                biasDec = ades_obs(i)%biasDec/secrad

                ! set bias computation flag to FALSE
                compute_errormodel_bias = .FALSE.

             ENDIF

             ! *** WEIGHTS *** for optical obs
             ! we store the rms (ra-dec) directly from input ADES file if present
             ! after this loop over ADES obs we will store the rms into wei_bia 

             ! RA and Dec if rmsRA-rmsDec present in input ADES file
             IF (ades_flg(i)%rmsRA .AND. ades_obs(i)%rmsRA .GT. 0.d0 .AND. .NOT. use_errmod(i) .AND. &
                  & .NOT.(ALL(ades_fit(i)%force_w_flags))) THEN

                ! store correlation and precision for output format
                ! rmsCorr is an optional element
                ! if only one weight is forced, set correlation to 0
                IF (ades_flg(i)%rmsCorr .AND. .NOT.(ANY(ades_fit(i)%force_w_flags)) &
                     & .AND. (ades_obs(i)%rmsCorr .GE. -1.d0 .AND. ades_obs(i)%rmsCorr .LE. 1.d0)) THEN
                   corr_radec(i) = ades_obs(i)%rmsCorr
                   ades_prec(i)%sigCorr = ades_prec(i)%rmsCorr
                ELSE
                   corr_radec(i) = 0.d0
                ENDIF

                ! convert to radians and update precision for output format
                IF (.NOT.ades_fit(i)%force_w_flags(1)) THEN
                   rmsRA = ades_obs(i)%rmsRA/(secrad*COS(dec))
                   ades_prec(i)%sigRA = ades_prec(i)%rmsRA
                ENDIF
                IF (.NOT.ades_fit(i)%force_w_flags(2)) THEN
                   rmsDec = ades_obs(i)%rmsDec/secrad
                   ades_prec(i)%sigDec = ades_prec(i)%rmsDec
                ENDIF
                
                ! compute RA and Dec accuracy
                CALL RADec_accuracy(i,dec,ra_acc,dec_acc,new_prec)

                ! rmsRA and rmsDec are the maximum (worst) between rms entries and corresponding accuracies
                IF (ra_acc .GT. rmsRA .AND. dec_acc .GT. rmsDec) THEN
                   rmsRA = ra_acc
                   rmsDec = dec_acc
                   corr_radec(i) = 0.d0
                   ! update precisions from RADec_accuracy
                   IF (.NOT.ades_fit(i)%force_w_flags(1)) ades_prec(i)%sigRA = new_prec(1)
                   IF (.NOT.ades_fit(i)%force_w_flags(2)) ades_prec(i)%sigDec = new_prec(2)
                   !ades_prec(i)%sigCorr = 0
                ELSEIF ((ra_acc .GT. rmsRA .AND. dec_acc .LT. rmsDec) .OR. &
                     &     (ra_acc .LT. rmsRA .AND. dec_acc .GT. rmsDec)) THEN
                   WRITE(obsnum,'(i7)')i
                   WRITE(ierrou,*) 'convert_ades_fitobs: ra_acc and dec_acc must be both GT or LT rmsRA and rmsDec &
                        &for input ADES file '//trim(ades_filename)//' obs number: '//trim(adjustl(obsnum))
                   numerr=numerr+1 
                ENDIF

                ! set rms computation flag to FALSE
                compute_errormodel_rms = .FALSE.

                ! RA and Dec if rmsDist-rmsPA present in input ADES file
             ELSEIF (ades_flg(i)%rmsDist .AND. ades_obs(i)%rmsDist .GT. 0.d0 .AND. .NOT. use_errmod(i) .AND. &
                  & .NOT.(ALL(ades_fit(i)%force_w_flags))) THEN

                ! get RA-Dec covariance matrix from rmsDist and rmsPA (and eventually from their correlation)
                DP_cov_mat(1,1) = (ades_obs(i)%rmsDist/secrad)**2
                DP_cov_mat(2,2) = (ades_obs(i)%rmsPA*radeg)**2
                IF (ades_flg(i)%rmsCorr .AND. .NOT.(ANY(ades_fit(i)%force_w_flags)) &
                     & .AND. (ades_obs(i)%rmsCorr .GE. -1.d0 .AND. ades_obs(i)%rmsCorr .LE. 1.d0)) THEN
                   DP_cov_mat(1,2) = ades_obs(i)%rmsCorr*ades_obs(i)%rmsDist/secrad*ades_obs(i)%rmsPA*radeg
                   DP_cov_mat(2,1) = DP_cov_mat(1,2)
                ELSE
                   DP_cov_mat(1,2) = 0.d0
                   DP_cov_mat(2,1) = 0.d0
                ENDIF
                CALL distPA_to_RADec_cov_mat(DP_cov_mat,dist,pa,dec_ref,RD_cov_mat)
                ! find rmsRA-rmsDec and eventually the correlation
                ! update precision for output format, in this case we set them to 6 decimal numbers (1.e-6 arcsec) by hand
                IF (.NOT.ades_fit(i)%force_w_flags(1)) THEN
                   rmsRA = SQRT(RD_cov_mat(1,1))
                   new_prec(1) = 6
                   ades_prec(i)%sigRA = new_prec(1)
                ENDIF
                IF (.NOT.ades_fit(i)%force_w_flags(2)) THEN
                   rmsDec = SQRT(RD_cov_mat(2,2))
                   new_prec(2) = 6
                   ades_prec(i)%sigDec = new_prec(2)
                ENDIF
                IF (ades_flg(i)%rmsCorr .AND. .NOT.(ANY(ades_fit(i)%force_w_flags)) &
                     & .AND. (ades_obs(i)%rmsCorr .GE. -1.d0 .AND. ades_obs(i)%rmsCorr .LE. 1.d0)) THEN
                   corr_radec(i) = RD_cov_mat(1,2)/(rmsRA*rmsDec)
                   ades_prec(i)%sigCorr = 12 !new_prec(1)+new_prec(2)
                ELSE
                   corr_radec(i) = 0.d0
                ENDIF

                ! compute RA and Dec accuracy
                CALL RADec_accuracy(i,dec,ra_acc,dec_acc,new_prec)

                ! rmsRA and rmsDec are the maximum (worst) between rms entries and corresponding accuracies
                IF (ra_acc .GT. rmsRA .AND. dec_acc .GT. rmsDec) THEN
                   rmsRA = ra_acc
                   rmsDec = dec_acc
                   corr_radec(i) = 0.d0
                   ! update precisions
                   IF (.NOT.ades_fit(i)%force_w_flags(1)) ades_prec(i)%sigRA = new_prec(1)
                   IF (.NOT.ades_fit(i)%force_w_flags(2)) ades_prec(i)%sigDec = new_prec(2)
                   !ades_prec(i)%sigCorr = 0
                ELSEIF ((ra_acc .GT. rmsRA .AND. dec_acc .LT. rmsDec) .OR. &
                     &     (ra_acc .LT. rmsRA .AND. dec_acc .GT. rmsDec)) THEN
                   WRITE(obsnum,'(i7)')i
                   WRITE(ierrou,*) 'convert_ades_fitobs: ra_acc and dec_acc must be both GT or LT rmsRA and rmsDec &
                        &for input ADES file '//trim(ades_filename)//' obs number: '//trim(adjustl(obsnum))
                   numerr=numerr+1 
                ENDIF

                ! set rms computation flag to FALSE
                compute_errormodel_rms = .FALSE.
             ENDIF
             
             ! update value for use_errmod(i) according to rms computation flag
             use_errmod(i) = compute_errormodel_rms

             ! computation of (optical) bias and rms via error model
             IF (.NOT. ades_flg(i)%rad .AND. (compute_errormodel_bias .OR. compute_errormodel_rms)) THEN
                ! error models
                init = .FALSE.
                CALL observ_rms_singobs(observ(i),name_em_loc(1:lenna),init,wei_bia(i))
                ! retrieve biases
                IF (compute_errormodel_bias) THEN
                   biasRA = wei_bia(i)%bias_coord(1)
                   biasDec = wei_bia(i)%bias_coord(2)
                ENDIF
                ! retrieve rms
                IF (compute_errormodel_rms) THEN
                   IF(.NOT.ades_fit(i)%force_w_flags(1)) THEN
                      rmsRA = wei_bia(i)%rms_coord(1)
                      ! update precision if rmsRA has changed (prec = 2 from error model)
                      IF (wei_bia(i)%rms_coord(1) .NE. observ(i)%acc_coord(1)) ades_prec(i)%sigRA = 2
                   ENDIF
                   IF(.NOT.ades_fit(i)%force_w_flags(2)) THEN
                      rmsDec = wei_bia(i)%rms_coord(2)
                      ! update precision if rmsDec has changed (prec = 2 from error model)
                      IF (wei_bia(i)%rms_coord(2) .NE. observ(i)%acc_coord(2)) ades_prec(i)%sigDec = 2
                   ENDIF
                   ! set correlation to 0
                   corr_radec(i) = 0.d0
                   ! update precision
                   ades_prec(i)%sigCorr = 0
                ENDIF
             ENDIF
          
          ELSE
             ! in psv input file are present sigmas and biases computed with a different error model
             ! thus compute all again using the requested error model
             ! forced weights (and their associated correlation) are left as input
             init = .FALSE.
             CALL observ_rms_singobs(observ(i),name_em_loc(1:lenna),init,wei_bia(i))
             ! retrieve biases
             biasRA = wei_bia(i)%bias_coord(1)
             biasDec = wei_bia(i)%bias_coord(2)
             ! retrieve rms
             rmsRA = wei_bia(i)%rms_coord(1)
             rmsDec = wei_bia(i)%rms_coord(2)
             ! update precisions if rms have changed (prec = 2 from error model)
             IF (wei_bia(i)%rms_coord(1) .NE. observ(i)%acc_coord(1)) ades_prec(i)%sigRA = 2
             IF (wei_bia(i)%rms_coord(2) .NE. observ(i)%acc_coord(2)) ades_prec(i)%sigDec = 2
          ENDIF

          ! store optical bias and weights 
          wei_bia(i)%bias_coord(1) = biasRA
          wei_bia(i)%bias_coord(2) = biasDec
          wei_bia(i)%rms_coord(1) = rmsRA
          wei_bia(i)%rms_coord(2) = rmsDec
          wei_bia(i)%correl = corr_radec(i)  
       ENDIF
       
       ! if ADES data contain optical residuals, store them
       IF (ades_flg(i)%resRA .AND. ades_flg(i)%resDec .AND. ades_fit(i)%ades_resc_def) THEN
          wei_bia(i)%resc_def=.TRUE.
          wei_bia(i)%res_coord(1) = ades_obs(i)%resRA/(secrad*COS(dec))
          wei_bia(i)%res_coord(2) = ades_obs(i)%resDec/secrad
       ENDIF

       ! radar obs
       IF (ades_flg(i)%rad) THEN

          ! *** WEIGHTS ***          
          rms_radf = 0.d0
          ! as for optical obs, rms is taken from sigma entry; if the latter is not present it is taken from rms entry (REQ)
          ! range
          IF (ades_flg(i)%sigDelay .AND. ades_obs(i)%sigDelay .GT. 0.d0) THEN
             rms_radf = 0.5d0*ades_obs(i)%sigDelay*1.d-6*vlight/8.64d4
          ELSEIF (ades_flg(i)%delay) THEN
             rms_radf = 0.5d0*ades_obs(i)%rmsDelay*1.d-6*vlight/8.64d4
             ades_prec(i)%sigDelay = ades_prec(i)%rmsDelay
          ENDIF
          ! range-rate
          IF (ades_flg(i)%sigDoppler .AND. ades_obs(i)%sigDoppler .GT. 0.d0) THEN
             rms_radf = 0.5d0*ades_obs(i)%sigDoppler*(1.d-6*vlight/ades_obs(i)%frq)
          ELSEIF (ades_flg(i)%doppler) THEN
             rms_radf = 0.5d0*ades_obs(i)%rmsDoppler*(1.d-6*vlight/ades_obs(i)%frq)
             ades_prec(i)%sigDoppler = ades_prec(i)%rmsDoppler
          ENDIF
          rms_rad = MAX(rms_radf,1.d-3/aukm)

          ! *** BIAS *** 
          ! bias for radar obs always 0
          bias_rad = 0.d0        
       ENDIF
      
       ! store radar bias and rms
       ! store observ(i)%acc_coord for radar observations
       IF (ades_flg(i)%delay) THEN
          wei_bia(i)%bias_coord(1) = bias_rad
          wei_bia(i)%rms_coord(1) = rms_rad
          observ(i)%acc_coord(1) = rms_radf
       ELSEIF (ades_flg(i)%doppler) THEN
          wei_bia(i)%bias_coord(2) = bias_rad
          wei_bia(i)%rms_coord(2) = rms_rad
          observ(i)%acc_coord(2) = rms_radf
       ENDIF
       ! if ADES data contain radar residuals, store them
       IF (ades_flg(i)%resDelay .AND. ades_fit(i)%ades_resc_def) THEN
          wei_bia(i)%resc_def=.TRUE.
          wei_bia(i)%res_coord(1) = ades_obs(i)%resDelay*vlight*1.d-6*0.5d0/8.64d4
       ELSEIF (ades_flg(i)%resDoppler .AND. ades_fit(i)%ades_resc_def) THEN
          wei_bia(i)%resc_def=.TRUE.
          wei_bia(i)%res_coord(2) = ades_obs(i)%resDoppler*vlight/(2.d6*ades_obs(i)%frq)
       ENDIF

       ! --- magnitude fields for wei_bia --- !

       IF (ades_flg(i)%mag) THEN
          ! fill wei_bia
          ! rms_mag
          wei_bia(i)%rms_mag = -1.d0 
          ! we use sigma mag if given in input ADES file
          IF (ades_flg(i)%sigMag .AND. ades_obs(i)%sigMag .GT. 0) THEN
             wei_bia(i)%rms_mag = ades_obs(i)%sigMag
             ! otherwise we use rms mag (if given)
          ELSE IF (ades_flg(i)%rmsMag .AND. ades_obs(i)%rmsMag .GT. 0) THEN
             wei_bia(i)%rms_mag = ades_obs(i)%rmsMag
             ades_prec(i)%sigMag = ades_prec(i)%rmsMag             
             ! otherwise we recompute
             ! notice that for number of digits after dot > 2 we still use the same rms as for 2 digits
          ELSE
             wei_bia(i)%rms_mag = magrms_new(observ(i)%mag_str,observ(i)%time_tdt,observ(i)%obscod_i,observ(i)%tech)   
             ! update precision
             ades_prec(i)%sigMag = 1
          ENDIF
          ! if ADES data contain photometric residuals, store them
          IF (ades_flg(i)%resMag .AND. ades_fit(i)%ades_resm_def) THEN
             wei_bia(i)%resm_def = .TRUE.
             wei_bia(i)%res_mag = ades_obs(i)%resMag
          ENDIF

       ELSE IF (.NOT. ades_flg(i)%rad) THEN                       
          wei_bia(i)%rms_mag = -1.d0
          wei_bia(i)%resm_def=.FALSE.
          wei_bia(i)%res_mag=0.d0
       ENDIF

    END DO
    
    ! Read .rules file dedicated to NEOCP observations with ADES RMS and assign weights
    IF (ons_name .AND. .NOT.(compute_all_errormodel) .AND. .NOT.(ALL(use_errmod))) THEN
       CALL stat_min_rms(nobs_ades,observ,min_rms)
       DO j = 1,nobs_ades
          IF (ades_obs(j)%rmsRA .NE. 0.d0 .AND. .NOT. use_errmod(j)) THEN  ! rmsRA and rmsDec were selected as rms_coord
             ! sum variance with credible minimum value read in neocp.rules only for NEOCP observations
             IF (.NOT.(ades_fit(j)%force_w_flags(1))) &
                  & wei_bia(j)%rms_coord(1) = SQRT(wei_bia(j)%rms_coord(1)**2 + min_rms(1,j)**2)
             IF (.NOT.(ades_fit(j)%force_w_flags(2))) &
                  & wei_bia(j)%rms_coord(2) = SQRT(wei_bia(j)%rms_coord(2)**2 + min_rms(2,j)**2)
          ENDIF
       END DO
    END IF
             
    ! find number of optical obs within the same batch (ending with an 8 hours gap) and apply sqrt(N) factor 
    ! also, update ades_obs with computed sigmas and biases
    
    ! sort by time first
    CALL heapsort(observ(1:nobs_ades)%time_tdt,nobs_ades,ind)
    DO i=1, nobs_ades
       obs1(i)=observ(ind(i))
    END DO
    
    ! find number of observations in batch from same station
    nob1=0  
    DO j=1, nobs_ades
       IF (obs1(j)%type .EQ. 'O' .AND. obs1(j)%obscod_s .NE. '247') THEN
          IF (nob1(j).GT.0) CYCLE ! already done
          t2=obs1(j)%time_tdt
          nob1(j)=1
          IF (obs1(j)%tech.EQ.'X'.OR.obs1(j)%tech.EQ.'x') CYCLE
          ! skip also obs having nonexistent technology descriptor 
          IF (ALL(tech_values.NE.obs1(j)%tech)) CYCLE
          ! find batch from same station
          DO k=j+1, nobs_ades
             IF (obs1(k)%time_tdt-t2.GT.gapmax) EXIT ! batch finished
             IF (obs1(k)%obscod_s.EQ.obs1(j)%obscod_s) THEN
                IF (obs1(k)%tech.EQ.'X'.OR.obs1(k)%tech.EQ.'x') CYCLE
                IF (ALL(tech_values.NE.obs1(k)%tech)) CYCLE
                nob1(j)=nob1(j)+1
                t2=obs1(k)%time_tdt
             ENDIF
          END DO
          ! assign number of obs in batch to each obs in batch
          DO k=j+1, nobs_ades
             IF(obs1(k)%time_tdt-t2.GT.gapmax) EXIT ! batch finished
             IF(obs1(k)%obscod_s.EQ.obs1(j)%obscod_s)THEN
                IF (obs1(k)%tech.EQ.'X'.OR.obs1(k)%tech.EQ.'x') CYCLE
                IF (ALL(tech_values.NE.obs1(k)%tech)) CYCLE
                nob1(k)=nob1(j)
             ENDIF
          END DO
       ENDIF
    END DO
    
    ! go back to original array of obs, apply sqrt(N) or sqrt(N/4) factor to RMS from error model
    DO i = 1, nobs_ades
       jj=ind(i)
       nob(jj)=nob1(i)
       factor = 1.d0
       IF (.NOT. ades_flg(jj)%rad .AND. .NOT.(ALL(ades_fit(jj)%force_w_flags)) .AND.  &
            &  (use_errmod(jj) .OR. compute_all_errormodel)) THEN
          ! rescale rms
          IF (nob(jj).NE.0) THEN
             IF(error_model.EQ.'vfcc17') THEN
                IF(nob(jj).GE.5) THEN
                   factor = SQRT(nob(jj)*0.25d0)
                ENDIF
             ELSE
                factor = SQRT(1.d0*nob(jj))
             ENDIF
             !IF (use_errmod(jj) .OR. compute_all_errormodel) THEN  ! non-forced weights have been assigned according to error model
             ! rescale with factor
             IF (.NOT.(ades_fit(jj)%force_w_flags(1))) wei_bia(jj)%rms_coord(1) = wei_bia(jj)%rms_coord(1)*factor
             IF (.NOT.(ades_fit(jj)%force_w_flags(2))) wei_bia(jj)%rms_coord(2) = wei_bia(jj)%rms_coord(2)*factor
          ENDIF
       ENDIF
    ENDDO

    DO i = 1, nobs_ades

       IF (.NOT. ades_flg(i)%rad) THEN
          ! update ades_obs with new sigmas and biases (optical case, already updated in radar case)
          dec = ades_obs(i)%dec*radeg
          ades_obs(i)%sigRA = wei_bia(i)%rms_coord(1)*secrad*COS(dec)
          ades_obs(i)%sigDec = wei_bia(i)%rms_coord(2)*secrad
          ades_obs(i)%sigCorr = corr_radec(i)  
          ades_obs(i)%biasRA = wei_bia(i)%bias_coord(1)*secrad*COS(dec)
          ades_obs(i)%biasDec = wei_bia(i)%bias_coord(2)*secrad
          IF (wei_bia(i)%rms_mag .NE. -1.d0) ades_obs(i)%sigMag = wei_bia(i)%rms_mag
          
          ! perform checks on sigmas
          IF (ades_obs(i)%sigRA .LT. 0.d0) THEN
             WRITE(obsnum,'(i7)')i
             WRITE(ierrou,*) 'convert_ades_fitobs: Computed sigRA negative &
                  &for input ADES file '//trim(ades_filename)//' obs number: '//trim(adjustl(obsnum))
             numerr=numerr+1
          ENDIF
          IF (ades_obs(i)%sigDec .LT. 0.d0) THEN
             WRITE(obsnum,'(i7)')i
             WRITE(ierrou,*) 'convert_ades_fitobs: Computed sigDec negative &
                  &for input ADES file '//trim(ades_filename)//' obs number: '//trim(adjustl(obsnum))
             numerr=numerr+1
          ENDIF
          IF (ades_obs(i)%sigCorr .GT. 1.d0 .OR. ades_obs(i)%sigCorr .LT. -1.d0) THEN
             WRITE(ierrou,*) 'convert_ades_fitobs: Computed sigCorr is not in [-1,1] interval &
                  &for input ADES file '//trim(ades_filename)//' obs number: '//trim(adjustl(obsnum))
             numerr = numerr+1
          ENDIF

          IF (ades_flg(i)%rmsRA) THEN
             IF (ades_obs(i)%sigRA*secrad*COS(dec) .LT. ades_obs(i)%rmsRA) THEN
                WRITE(obsnum,'(i7)')i
                WRITE(ierrou,*) 'convert_ades_fitobs: Computed sigRA less than rmsRA &
                     &for input ADES file '//trim(ades_filename)//' obs number: '//trim(adjustl(obsnum))
                numerr=numerr+1
             ENDIF
             IF (ades_obs(i)%sigDec*secrad .LT. ades_obs(i)%rmsDec) THEN
                WRITE(obsnum,'(i7)')i
                WRITE(ierrou,*) 'convert_ades_fitobs: Computed sigDec less than rmsDec &
                     &for input ADES file '//trim(ades_filename)//' obs number: '//trim(adjustl(obsnum))
                numerr=numerr+1
             ENDIF
          ENDIF
          IF (ades_obs(i)%sigMag .LT. 0.d0) THEN
             WRITE(obsnum,'(i7)')i
             WRITE(ierrou,*) 'convert_ades_fitobs: Computed sigMag negative &
                  &for input ADES file '//trim(ades_filename)//' obs number: '//trim(adjustl(obsnum))
             numerr=numerr+1
          ENDIF
          IF (ades_flg(i)%rmsMag) THEN
             IF (ades_obs(i)%sigMag .LT. ades_obs(i)%rmsMag) THEN
                WRITE(obsnum,'(i7)')i
                WRITE(ierrou,*) 'convert_ades_fitobs: Computed sigMag less than rmsMag &
                     &for input ADES file '//trim(ades_filename)//' obs number: '//trim(adjustl(obsnum))
                numerr=numerr+1
             ENDIF
          ENDIF
       ELSE         
          IF (ades_flg(i)%delay) THEN

             ades_obs(i)%sigDelay = wei_bia(i)%rms_coord(1)*2.d6*8.64d4/vlight

             IF (ades_flg(i)%rmsDelay) THEN
                IF (ades_obs(i)%sigDelay .LT. ades_obs(i)%rmsDelay) THEN  
                   WRITE(obsnum,'(i7)')i
                   WRITE(ierrou,*) 'convert_ades_fitobs: Computed sigDelay less than rmsDelay &
                        &for input ADES file '//trim(ades_filename)//' obs number: '//trim(adjustl(obsnum))
                   numerr = numerr+1
                ENDIF
             ENDIF

          ELSE

             ades_obs(i)%sigDoppler = wei_bia(i)%rms_coord(2)*2.d6*ades_obs(i)%frq/vlight

             IF (ades_flg(i)%rmsDoppler) THEN
                IF (ades_obs(i)%sigDoppler .LT. ades_obs(i)%rmsDoppler) THEN
                   WRITE(obsnum,'(i7)')i
                   WRITE(ierrou,*) 'convert_ades_fitobs: Computed sigDoppler less than rmsDoppler &
                        &for input ADES file '//trim(ades_filename)//' obs number: '//trim(adjustl(obsnum))
                   numerr = numerr+1
                ENDIF
             ENDIF
          ENDIF
       ENDIF   
    END DO
  END SUBROUTINE convert_ades_fitobs

  ! ------------------------------------------------------ !
  ! SUBROUTINE convert_fitobs_ades                         !
  ! It fills the ADES ObsData fields in ades data types.   !
  ! If a .psv/.xml file was given in input, it only        !
  ! converts and prepares the output (residuals group +    !
  ! chi + selection flags) to be written to output ADES    !
  ! file, otherwise it fills also all the suitable and     !
  ! required ADES fields. It also performs checks on the   !
  ! output values to ensure they are consistent with the   !
  ! ADES format rules.                                     !
  ! ------------------------------------------------------ !
  ! Author: Alessia Bertolucci, SpaceDyS 2020              !
  ! ------------------------------------------------------ !

  SUBROUTINE convert_fitobs_ades(observ,wei_bia,ades_found,RMS,disc_flag,freq,subFmt)
    
    TYPE(ast_obs),     INTENT(IN)    :: observ(nobs_ades)               ! observations
    TYPE(ast_wbsr),    INTENT(IN)    :: wei_bia(nobs_ades)              ! rms and biases
    LOGICAL,           INTENT(IN)    :: ades_found                      ! TRUE if ADES file has been found
    LOGICAL,           INTENT(IN)    :: RMS                             ! TRUE if there is global RMS of astrometric/photometric fit
    CHARACTER(LEN=1),  INTENT(IN), OPTIONAL :: disc_flag(nobs_ades)     ! discovery flag
    REAL(KIND=dkind),  INTENT(IN), OPTIONAL :: freq(nobs_ades)          ! reference frequency
    CHARACTER(LEN=4),  INTENT(IN), OPTIONAL :: subFmt(nobs_ades)        ! Format in which the observation was originally submitted to the MPC

    INTEGER            :: i,j,count_objdes(nobs_ades),lendes,l
    CHARACTER(LEN=20)  :: objdes,temp_objdes1,temp_objdes2,count_str,des
    REAL(KIND=dkind)   :: cosd
    CHARACTER(LEN=7)   :: obsnum,auxstr
    CHARACTER(LEN=500) :: err_str,data_str

    ! initialization
    ades_fit = undefined_ades_fit_var

    ! if ades file has not been found, fills all the ades fields
    IF (.NOT.ades_found) THEN
       ades_obs = undefined_ades_observ
       ades_flg = undefined_ades_flags
       ades_prec = undefined_ades_precision
            
       ! loop over ADES observations
       DO i=1, nobs_ades
          CALL new_ades_obs(i,observ,wei_bia,disc_flag,freq,subFmt)
       END DO
       
    ELSE

       ! if an input psv file exists, fill only residuals group fields
       ! loop over observations
       DO i=1,nobs_ades
          ! prepare end of error string
          WRITE(obsnum,'(i7)')i
          err_str = 'for ADES file: '//trim(psv_name(1:len_ades_name))//' obs number: '//trim(adjustl(obsnum))
          
          ! compute COS(Dec)
          cosd = COS(observ(i)%coord(2))

          ! if designation is provisional, it has been modified in convert_ades_fitobs to fit observ%objdes value
          ! then re-write it according to ADES default format (space between year and letters)
          IF (ades_flg(i)%provID) THEN
             des = trim(adjustl(ades_obs(i)%provID))
             CALL rmsp(des,l)
             des = des(1:4)//" "//des(5:l)   ! according to ADES formatting guidelines
             ades_obs(i)%provID = des
          ENDIF

          ! fill chi
          IF (wei_bia(i)%resc_def) ades_fit(i)%chi = wei_bia(i)%chi
          IF (ades_fit(i)%chi .LT. 0.d0) THEN
             WRITE(ierrou,*) 'ERROR: convert_fitobs_ades: computed chi is negative '//err_str
             numerr = numerr+1
          ENDIF

          ! fill force_w_flags
          ades_fit(i)%force_w_flags(1) = wei_bia(i)%force_w(1)
          ades_fit(i)%force_w_flags(2) = wei_bia(i)%force_w(2)
                    
          ! fill orbProd
          ades_flg(i)%orbProd = .TRUE.
          ades_obs(i)%orbProd = "OrbFit5.0"

          ! optical obs
          IF (.NOT.ades_flg(i)%rad) THEN

             ! fill ast_flag and mag_flag
             WRITE(ades_fit(i)%ast_flag,'(i1)') wei_bia(i)%sel_coord
             WRITE(ades_fit(i)%mag_flag,'(i1)') wei_bia(i)%sel_mag

             ! fill resRA and resDec (convert to arcsec and apply COS(Dec) factor)
             ades_obs(i)%resRA = wei_bia(i)%res_coord(1)*secrad*cosd
             ades_obs(i)%resDec = wei_bia(i)%res_coord(2)*secrad
             ades_flg(i)%resRA = .TRUE.
             ades_flg(i)%resDec = .TRUE.
             ades_fit(i)%ades_resc_def = wei_bia(i)%resc_def

             ! fill selAst
             ades_flg(i)%selAst = .TRUE.
             IF (wei_bia(i)%sel_coord .EQ. 0) THEN
                ades_obs(i)%selAst = "D"
             ELSE
                ades_obs(i)%selAst = "A" 
             ENDIF

             ! fill sigRA, sigDec, sigCorr: already stored
             ades_flg(i)%sigRA = .TRUE.
             WRITE(data_str,*) ades_obs(i)%sigRA
             ades_prec(i)%sigRA = get_prec(data_str)
             ades_flg(i)%sigDec = .TRUE.
             WRITE(data_str,*) ades_obs(i)%sigDec
             ades_prec(i)%sigDec = get_prec(data_str)
             ades_flg(i)%sigCorr = .TRUE.
             WRITE(data_str,*) ades_obs(i)%sigCorr
             ades_prec(i)%sigCorr = get_prec(data_str)
                       
             ! control on sigmas
             IF (ades_obs(i)%sigRA .LT. 0.d0) THEN
                WRITE(ierrou,*) 'ERROR: convert_fitobs_ades: computed sigRA is negative '//err_str
                numerr = numerr+1
             ENDIF
             IF (ades_obs(i)%sigDec .LT. 0.d0) THEN
                WRITE(ierrou,*) 'ERROR: convert_fitobs_ades: computed sigDec is negative '//err_str
                numerr = numerr+1
             ENDIF
             IF (ades_obs(i)%sigCorr .GT. 1.d0 .OR. ades_obs(i)%sigCorr .LT. -1.d0) THEN
                WRITE(ierrou,*) 'ERROR: convert_fitobs_ades: computed sigCorr is not in [-1,1] interval'//err_str
                numerr = numerr+1
             ENDIF

             ! biasRA and biasDec already computed
             ades_flg(i)%biasRA = .TRUE.
             ades_flg(i)%biasDec = .TRUE.

             ! sigTime and biasTime obs stored with their default values, if not present  
             ades_flg(i)%sigTime = .TRUE.
             ades_flg(i)%biasTime= .TRUE.

             ! fill photProd if not present (otherwise it is already stored)
             IF (.NOT. ades_flg(i)%photProd) THEN
                ades_obs(i)%photProd = ades_obs(i)%orbProd
                ades_flg(i)%photProd = .TRUE.
             ENDIF

             ! fill resMag
             ades_flg(i)%resMag = .TRUE.
             ades_obs(i)%resMag = wei_bia(i)%res_mag
             ades_fit(i)%ades_resm_def = wei_bia(i)%resm_def

             ! fill selPhot, if not stored yet
             ades_flg(i)%selPhot = .TRUE.
             IF (wei_bia(i)%sel_mag .EQ. 0) THEN
                ades_obs(i)%selPhot = "D"
             ELSE
                ades_obs(i)%selPhot = "A" 
             ENDIF

             ! sigMag is already stored
             ades_flg(i)%sigMag = .TRUE.
             WRITE(data_str,*) ades_obs(i)%sigMag
             ades_prec(i)%sigMag = get_prec(data_str)
                          
             ! control on sigma
             IF (ades_obs(i)%sigMag .LT. 0.d0) THEN
                WRITE(ierrou,*) 'ERROR: convert_fitobs_ades: computed sigMag is negative '//err_str
                numerr = numerr+1
             ENDIF

             ! biasMag
             ades_flg(i)%biasMag = .TRUE.
             ! default is zero
             ades_obs(i)%biasMag = 0.0
             ! assign bias corresponding to band (values for z,i,g,r,y,w are provided by Fitzsimmons November 2011 report)
             IF(trim(adjustl(ades_obs(i)%band)) .EQ. 'V' .OR. trim(adjustl(ades_obs(i)%band)) .EQ. 'H') THEN
                ades_obs(i)%biasMag = +0.00d0
             ELSE IF(trim(adjustl(ades_obs(i)%band)) .EQ. ' ' .OR. trim(adjustl(ades_obs(i)%band)) .EQ. 'B') THEN
                ades_obs(i)%biasMag = -0.80d0
             ELSE IF(trim(adjustl(ades_obs(i)%band)) .EQ. 'I') THEN
                ades_obs(i)%biasMag = +0.80d0
             ELSE IF(trim(adjustl(ades_obs(i)%band)) .EQ. 'R' .OR. trim(adjustl(ades_obs(i)%band)) .EQ. 'C') THEN
                ades_obs(i)%biasMag = +0.40d0
             ELSE IF(trim(adjustl(ades_obs(i)%band)) .EQ. 'z') THEN
                ades_obs(i)%biasMag = +0.37d0
             ELSE IF(trim(adjustl(ades_obs(i)%band)) .EQ. 'i') THEN
                ades_obs(i)%biasMag = +0.39d0
             ELSE IF(trim(adjustl(ades_obs(i)%band)) .EQ. 'g') THEN
                ades_obs(i)%biasMag = -0.28d0
             ELSE IF(trim(adjustl(ades_obs(i)%band)) .EQ. 'U' .OR. trim(adjustl(ades_obs(i)%band)) .EQ. 'J') THEN
                ades_obs(i)%biasMag = +1.10d0
             ELSE IF(trim(adjustl(ades_obs(i)%band)) .EQ. 'r') THEN
                ades_obs(i)%biasMag = +0.23d0
             ELSE IF(trim(adjustl(ades_obs(i)%band)) .EQ. 'y') THEN
                ades_obs(i)%biasMag = +0.36d0
             ELSE IF(trim(adjustl(ades_obs(i)%band)) .EQ. 'w') THEN
                ades_obs(i)%biasMag = +0.16d0
             ELSE IF(trim(adjustl(ades_obs(i)%band)) .EQ. 'G') THEN
                ades_obs(i)%biasMag = +0.28d0                       ! Gaia correction by Tholen
             ELSE IF(trim(adjustl(ades_obs(i)%band)) .EQ. 'c') THEN
                ades_obs(i)%biasMag = -0.05d0
             ELSE IF(trim(adjustl(ades_obs(i)%band)) .EQ. 'o') THEN
                ades_obs(i)%biasMag = +0.33d0
             ELSE IF(trim(adjustl(ades_obs(i)%band)) .EQ. 'L') THEN
                ades_obs(i)%biasMag =  +0.20d0
             ELSE IF(trim(adjustl(ades_obs(i)%band)) .EQ. 'K') THEN
                ades_obs(i)%biasMag = +1.70d0
             ELSE IF(trim(adjustl(ades_obs(i)%band)) .EQ. 'Y') THEN
                ades_obs(i)%biasMag = +0.70d0
             ELSE IF(trim(adjustl(ades_obs(i)%band)) .EQ. 'u') THEN
                ades_obs(i)%biasMag = +2.50d0
             ELSE
                WRITE(ierrou,*) 'ERROR: convert_fitobs_ades: unknown color, unable to find corresponding biasMag '//err_str 
                numerr=numerr+1
             END IF

             ! set photometric model photMod
             ades_flg(i)%photMod = .TRUE.
             ades_obs(i)%photMod = "H"

             ! radar obs
          ELSE

             ! fill radar_flag
             WRITE(ades_fit(i)%radar_flag,'(i1)') wei_bia(i)%sel_coord

             ! fill residuals for delay
             IF (ades_flg(i)%delay) THEN

                ! fill resDelay
                ades_flg(i)%resDelay = .TRUE.
                ades_obs(i)%resDelay = 2.d0*wei_bia(i)%res_coord(1)*1.d6*8.64d4/vlight
                ades_fit(i)%ades_resc_def = wei_bia(i)%resc_def

                ! fill selDelay, if not stored yet
                ades_flg(i)%selDelay = .TRUE.
                IF (wei_bia(i)%sel_coord .EQ. 0) THEN
                   ades_obs(i)%selDelay = "D"
                ELSE
                   ades_obs(i)%selDelay = "A" 
                ENDIF

                ! sigDelay is already stored
                ades_flg(i)%sigDelay = .TRUE.
                WRITE(data_str,*) ades_obs(i)%sigDelay
                ades_prec(i)%sigDelay = get_prec(data_str)
                
                ! control on sigma
                IF (ades_obs(i)%sigDelay .LE. 0.d0) THEN
                   WRITE(ierrou,*) 'ERROR: convert_fitobs_ades: computed sigDelay is negative  '//err_str
                   numerr = numerr+1
                ENDIF

                ! fill residuals for Doppler
             ELSE
                ! fill resDoppler
                ades_flg(i)%resDoppler = .TRUE.
                ades_obs(i)%resDoppler = 2.d6*wei_bia(i)%res_coord(2)*ades_obs(i)%frq/vlight
                ades_fit(i)%ades_resc_def = wei_bia(i)%resc_def
                
                ! fill selDoppler
                ades_flg(i)%selDoppler = .TRUE.
                IF (wei_bia(i)%sel_coord .EQ. 0) THEN
                   ades_obs(i)%selDoppler = "D"
                ELSE
                   ades_obs(i)%selDoppler = "A" 
                ENDIF

                ! sigDoppler is already stored
                ades_flg(i)%sigDoppler = .TRUE.
                WRITE(data_str,*) ades_obs(i)%sigDoppler
                ades_prec(i)%sigDoppler = get_prec(data_str)
                
                ! control on sigma
                IF (ades_obs(i)%sigDoppler .LE. 0.d0) THEN
                   WRITE(ierrou,*) 'ERROR: convert_fitobs_ades: computed sigDoppler is negative '//err_str
                   numerr = numerr+1
                ENDIF

             ENDIF
          ENDIF
       END DO
    ENDIF
    
    ! fills the other Residuals Group fields which are handled in the same way in both cases (existence or non existence of input psv file)
    ! count how many times an object appears in ades_obs (and then set objdes for filling orbID entry)
    count_objdes = 1
    DO i=1,nobs_ades

       IF (.NOT. ades_flg(i)%rad) THEN
          IF (ades_flg(i)%permID) THEN
             temp_objdes1 = ades_obs(i)%permID
          ELSEIF (ades_flg(i)%provID) THEN
             temp_objdes1 = ades_obs(i)%provID
             CALL rmsp(temp_objdes1,lendes)
          ELSE
             temp_objdes1 = ades_obs(i)%trkSub
          ENDIF
       ELSE
          IF (ades_flg(i)%permID) THEN
             temp_objdes1 = ades_obs(i)%permID
          ELSE
             temp_objdes1 = ades_obs(i)%provID
             CALL rmsp(temp_objdes1,lendes)
          ENDIF
       ENDIF

       IF (i .GT. 1) THEN
          j = i-1
          DO WHILE (j .GE. 1)

             IF (.NOT. ades_flg(j)%rad) THEN
                IF (ades_flg(j)%permID) THEN
                   temp_objdes2 = ades_obs(j)%permID
                ELSEIF (ades_flg(j)%provID) THEN
                   temp_objdes2 = ades_obs(j)%provID
                   CALL rmsp(temp_objdes2,lendes)
                ELSE
                   temp_objdes2 = ades_obs(j)%trkSub
                ENDIF
             ELSE
                IF (ades_flg(j)%permID) THEN
                   temp_objdes2 = ades_obs(j)%permID
                ELSE
                   temp_objdes2 = ades_obs(j)%provID
                   CALL rmsp(temp_objdes2,lendes)
                ENDIF
             ENDIF

             IF ((temp_objdes2 .EQ. temp_objdes1) .AND. (ades_flg(j)%rad .EQV. ades_flg(i)%rad) .AND. &
                  & (ades_flg(j)%opt .EQV. ades_flg(i)%opt) .AND. (ades_flg(j)%occ .EQV. ades_flg(i)%occ)) THEN
                count_objdes(i) = count_objdes(i) + 1
             ENDIF

             j = j-1

          ENDDO

       ENDIF

    ENDDO
    ! loop over observations
    DO i=1,nobs_ades
       ! prepare end of error string
       WRITE(obsnum,'(i7)')i
       err_str = 'for ADES file: '//trim(psv_name(1:len_ades_name))//' obs number: '//trim(adjustl(obsnum))

       ! fill orbID
       IF (.NOT. ades_flg(i)%rad) THEN
          IF (ades_flg(i)%permID) THEN
             objdes = ades_obs(i)%permID
          ELSEIF (ades_flg(i)%provID) THEN
             objdes = ades_obs(i)%provID
             CALL rmsp(objdes,lendes)
          ELSE
             objdes = ades_obs(i)%trkSub
             CALL rmsp(objdes,lendes)
          ENDIF
       ELSE
          IF (ades_flg(i)%permID) THEN
             objdes = ades_obs(i)%permID
          ELSE
             objdes = ades_obs(i)%provID
             CALL rmsp(objdes,lendes)
          ENDIF
       ENDIF
       ades_flg(i)%orbID = .TRUE.
       WRITE(count_str,"(i20)") count_objdes(i)
       IF (ades_flg(i)%opt) THEN
          ades_obs(i)%orbID = "OrbFit5.0_"//trim(adjustl(objdes))//"_#"//trim(adjustl(count_str))//"opt"
       ELSEIF (ades_flg(i)%occ) THEN
          ades_obs(i)%orbID = "OrbFit5.0_"//trim(adjustl(objdes))//"_#"//trim(adjustl(count_str))//"occ"
       ELSEIF (ades_flg(i)%rad) THEN
          ades_obs(i)%orbID = "OrbFit5.0_"//trim(adjustl(objdes))//"_#"//trim(adjustl(count_str))//"rad"
       ENDIF

    ENDDO

    ! if there is not an input ADES file, store arrays ades_key_rec and ades_data_rec for writing .psv
    IF (.NOT. ades_found) THEN
       CALL fill_psv_rec(new_ades,RMS)
    ENDIF
    
  END SUBROUTINE convert_fitobs_ades

  !----------------------------------------------------------- !
  ! SUBROUTINE new_ades_obs                                    !
  ! For a single observation, it fills the residuals group     !
  ! (except for orbID, that is handled in convert_fitobs_ades) !
  ! and also all the suitable and required ADES fields.        !
  ! It also performs checks on the output values to ensure     !
  ! they are consistent with the ADES format rules.            !
  ! ---------------------------------------------------------- !
  ! Author: Alessia Bertolucci - SpaceDyS - 2020               !
  ! ---------------------------------------------------------- !

  SUBROUTINE new_ades_obs(i,observ,wei_bia,disc,freq,subfmt)

    INTEGER,           INTENT(IN)    :: i                              ! index of observation
    TYPE(ast_obs),     INTENT(IN)    :: observ(nobs_ades)              ! observations
    TYPE(ast_wbsr),    INTENT(IN)    :: wei_bia(nobs_ades)             ! rms and biases
    
    CHARACTER(LEN=1),  INTENT(IN), OPTIONAL :: disc(nobs_ades)         ! ADES field 'disc' (discovery flag: * or +)
    REAL(KIND=dkind),  INTENT(IN), OPTIONAL :: freq(nobs_ades)         ! reference frequency
    CHARACTER(LEN=4),  INTENT(IN), OPTIONAL :: subfmt(nobs_ades)       ! Format in which the observation was originally submitted to the MPC
     
    INTEGER            :: j,count_objdes(nobs_ades),l,lenfile,unit,ind_point,len_month,len_d,len_h,len_min,len_sec
    CHARACTER(LEN=20)  :: objdes,temp_objdes1,temp_objdes2,count_str
    REAL(KIND=dkind)   :: cosd
    REAL(KIND=dkind)   :: obs_pos(3)              ! observer's geocentric position (au)
    CHARACTER(LEN=7)   :: obsnum,auxstr
    CHARACTER(LEN=1)   :: mode
    CHARACTER(LEN=200) :: rec
    CHARACTER(LEN=20)  :: des,file
    CHARACTER(LEN=500) :: err_str
    CHARACTER(LEN=40)  :: time_str,data_str,date,time
    CHARACTER(LEN=3)   :: auxstr1,auxstr2,ades_mode
    CHARACTER(LEN=6)   :: notes
    LOGICAL            :: perm_des,prov_des 

    ! Times
    INTEGER            :: mjd_utc,mjd_out,mjdd   ! integer MJD date (1 day=86400.0 s)
    REAL(KIND=dkind)   :: tjm_utc,tjmd,tjm_ut1   ! real MJD date
    REAL(KIND=dkind)   :: secs,sec_out,seconds,hour,minutes    
    INTEGER            :: year,month,day
    CHARACTER(LEN=4)   :: year_str
    CHARACTER(LEN=2)   :: month_str,day_str
    CHARACTER(LEN=10)  :: sec_str,hour_str,min_str

    LOGICAL, EXTERNAL  :: isnum,islett

    ! compute COS(Dec)
    cosd = COS(observ(i)%coord(2))

    notes = ""
    
    ! prepare end of error string
    WRITE(obsnum,'(i7)')i
    err_str = 'for ADES file: '//trim(psv_name(1:len_ades_name))//' obs number: '//trim(adjustl(obsnum))
    
    ! Identify type of observation and fill corresponding ades flags
    ! optical obs
    IF (trim(adjustl(observ(i)%type)) .EQ. 'O') THEN
       IF (trim(adjustl(observ(i)%tech)) .EQ. 'E') THEN
          ades_flg(i)%occ = .TRUE.   ! it is an occultation observation
       ELSE
          ades_flg(i)%opt = .TRUE.
          IF (trim(adjustl(observ(i)%obscod_s)) .EQ. '247') ades_flg(i)%sys = .TRUE.
       ENDIF
    ELSEIF (trim(adjustl(observ(i)%type)) .EQ. 'S') THEN
       ades_flg(i)%opt = .TRUE.
       ades_flg(i)%sys = .TRUE.
    ELSEIF (trim(adjustl(observ(i)%type)) .EQ. 'R'.OR. trim(adjustl(observ(i)%type)) .EQ. 'V') THEN
       ades_flg(i)%rad = .TRUE.
    ENDIF

    ! fill Identification Group using observ%objdes 
    des = trim(adjustl(observ(i)%objdes))
    CALL rmsp(des,l)
    perm_des = isnum(des(1:l))
    !com_des = isnum(des(1:l-1))
    prov_des = isnum(des(1:4)) .AND. islett(des(5:6))
    ! IF (perm_des .OR. (com_des .AND. (des(l:l) .EQ. 'P' .OR. des(l:l) .EQ. 'D'))) THEN
    ! objdes is a number (minor planet) or a number ending with D or P (comet)
    IF (perm_des) THEN
       ! objdes is a number (minor planet) 
       ades_obs(i)%permID = des                                   ! this means the designation is permanent 
       ades_flg(i)%permID = .TRUE.
       ! ELSEIF (.NOT. perm_des .AND. .NOT. com_des .AND. (prov_des .OR. ades_flg(i)%rad)) THEN
    ELSEIF (.NOT. perm_des .AND. (prov_des .OR. ades_flg(i)%rad)) THEN
       ! objdes is a string containing first 4 numbers (year) and then letters thus the designation is provisional
       ! radar observations must always include at least one of permID or provID
       des = des(1:4)//" "//des(5:l)   ! according to ADES formatting guidelines
       ades_obs(i)%provID = des
       ades_flg(i)%provID = .TRUE.
    ELSE
       ades_obs(i)%trkSub = des
       ades_flg(i)%trkSub = .TRUE.
    ENDIF
    
    ! fill mode (this field is forbidden in case of occultation or radar observations)
    IF (.NOT. ades_flg(i)%rad) THEN 
       mode = trim(adjustl(observ(i)%tech))
       IF (mode .EQ. "X" .OR. mode .EQ. "x") THEN
          ades_obs(i)%deprecated = "X"
          ades_flg(i)%deprecated = .TRUE.
       ENDIF
       IF (.NOT. ades_flg(i)%occ) THEN 
          ades_mode = "UNK"
          SELECT CASE (mode)
          CASE("P")
             ades_mode = "PHO"
          CASE("e")
             ades_mode = "ENC"
          CASE("C")
             ades_mode = "CCD"
          CASE("T")
             ades_mode = "MER"
          CASE("M")
             ades_mode = "MIC"
          CASE("c")
             ades_mode = "CCD"
          CASE("N")
             !ades_mode = "NOR"
             notes = "n"
          CASE("n")
             ades_mode = "VID"
          CASE("A")
             ades_mode = "PHO"
             !CASE("E")
             !ades_mode = "OCC"
          CASE("S")
             ades_mode = "CCD"   ! satellite obs set as CCD
          CASE("V")
             ades_mode = "CCD"   ! roving obs set as CCD
          END SELECT
          IF (ades_mode .EQ. "UNK" .AND. .NOT. ades_flg(i)%deprecated) THEN
             WRITE(ierrou,*) 'ERROR: convert_fitobs_ades: unable to convert to mode: tech = '//trim(adjustl(mode))
             numerr=numerr+1 
          ENDIF
          ades_obs(i)%mode = ades_mode 
          ades_flg(i)%mode = .TRUE.
       ENDIF
    ENDIF

    ! fill notes (and prog)
    IF (.NOT. ades_flg(i)%rad) THEN
       ades_flg(i)%notes = .TRUE.
       IF (ALL(notes_values .NE. observ(i)%note)) THEN
          notes = trim(adjustl(notes))
          ! the discarded value is related to the program used, thus it is put in %prog field
          ades_flg(i)%prog = .TRUE.
          IF (observ(i)%note == "|") THEN
             ades_obs(i)%prog = "pi"
          ELSE
             ades_obs(i)%prog = observ(i)%note
          ENDIF
       ELSE
          notes = trim(adjustl(notes))//observ(i)%note
       ENDIF
       ades_obs(i)%notes = notes
    ENDIF
    
    ! fill stn, trx/rcv
    auxstr = trim(adjustl(observ(i)%obscod_s))
    IF (.NOT. ades_flg(i)%rad) THEN
       ades_obs(i)%stn = auxstr
       ades_flg(i)%stn = .TRUE.
    ELSE
       auxstr1 = auxstr(1:3)
       auxstr2 = auxstr(5:7)
       ades_obs(i)%trx = auxstr1
       ades_obs(i)%rcv = auxstr2
       ades_flg(i)%trx = .TRUE.
       ades_flg(i)%rcv = .TRUE.
    ENDIF

    ! fill Location Group
    IF (ades_flg(i)%sys) THEN
       IF (auxstr .EQ. '247') THEN
          ades_obs(i)%sys = 'WGS84'
          ades_obs(i)%pos(1) = observ(i)%geopos(1)/radeg
          ades_obs(i)%pos(2) = observ(i)%geopos(2)/radeg
          ades_obs(i)%pos(3) = observ(i)%geopos(3)
       ELSE
          ades_obs(i)%sys = 'ICRF_KM'
          obs_pos = observ(i)%obspos*aukm  ! convert from au to km
          ! apply rotation from Ecliptic J2000.0 reference system into equatorial coordinates
          obs_pos = MATMUL(roteceq,obs_pos)
          ades_obs(i)%pos = obs_pos
       ENDIF
       ades_flg(i)%pos = .TRUE.
       ades_obs(i)%ctr = 399
       ades_flg(i)%ctr = .TRUE.
       ! Covariance matrix of (pos1,pos2,pos3) is left to default value
       !ades_flg(i)%posCov = .TRUE.
    ENDIF

    ! fill obsTime
    tjm_utc = observ(i)%time_utc
    ! computation of UTC calendar date and time of observation
    CALL mjddat(tjm_utc,day,month,year,hour)
    secs = hour*3600     ! seconds within a day
    mjd_utc = FLOOR(tjm_utc)
    ! IF (year .LT. 1972) THEN  ! if observation was taken before 1972, scale was UT1 and then we need to convert UTC calendar date and time
    !    CALL cnvtim(mjd_utc,secs,'UTC',mjd_out,sec_out,'UT1')
    !    mjdd = FLOOR(sec_out/86400.d0)
    !    mjd_out = mjd_out + mjdd
    !    sec_out = sec_out-mjdd*86400.d0  !seconds within a day
    !    tjm_ut1 = mjd_out + sec_out/86400.d0
    !    CALL mjddat(tjm_ut1,day,month,year,hour) 
    ! ELSE
    !    mjd_out = mjd_utc
    !    sec_out = secs
    ! ENDIF
    mjd_out = mjd_utc
    sec_out = secs
    hour = FLOOR(hour)
    minutes = FLOOR((sec_out-hour*3600)/60.d0)
    seconds = sec_out-hour*3600-minutes*60
    WRITE(year_str,'(I4)') year
    WRITE(month_str,'(I2)') month
    WRITE(day_str,'(I2)') day
    WRITE(hour_str,'(F3.0)') hour
    WRITE(min_str,'(F3.0)') minutes
    WRITE(sec_str,'(F6.3)') seconds
   
    ! adjust length of time strings to fit ADES format hh:mm:ss.sss (h and m integer without decimal point)
    ind_point = INDEX(hour_str,'.')
    IF (ind_point .GT. 0) hour_str = hour_str(1:ind_point-1)
    ind_point = INDEX(min_str,'.')
    IF (ind_point .GT. 0) min_str = min_str(1:ind_point-1)

    year_str = trim(adjustl(year_str))
    month_str = trim(adjustl(month_str))
    len_month = len(trim(adjustl(month_str)))
    day_str = trim(adjustl(day_str))
    len_d = len(trim(adjustl(day_str)))
    hour_str = trim(adjustl(hour_str))
    len_h = len(trim(adjustl(hour_str)))
    min_str = trim(adjustl(min_str))
    len_min = len(trim(adjustl(min_str)))
    sec_str = trim(adjustl(sec_str))
    len_sec = len(trim(adjustl(sec_str)))

    !control on length, according to default ADES formatting rules: yyyy-mm-ddThh:mm:ss.sssZ
    IF (len_month .EQ. 1) THEN
       month_str = "0"//month_str
       len_month = 2
    ENDIF
    IF (len_d .EQ. 1) THEN
       day_str = "0"//day_str
       len_d = 2
    ENDIF
    IF (len_h .EQ. 1) THEN
       hour_str = "0"//hour_str
       len_h = 2
    ENDIF
    IF (len_min .EQ. 1) THEN
       min_str = "0"//min_str
       len_min = 2
    ENDIF
    IF (len_sec .LT. 6) THEN   
       sec_str = "0"//sec_str
       len_sec = len_sec + 1
       IF (len_sec .NE. 6) WRITE(ierrou,*) 'ERROR: new_ades_obs: wrong format in time string for ObsTime field'
    ENDIF

    date = year_str//'-'//month_str(1:len_month)//'-'//day_str(1:len_d)
    time = 'T'//hour_str(1:len_h)//':'//min_str(1:len_min)//':'//sec_str(1:len_sec)//'Z'
    ades_obs(i)%obsTime = trim(adjustl(date))//trim(adjustl(time))
    ades_flg(i)%obsTime = .TRUE.

    IF (ades_flg(i)%rad) THEN
       ! fill frq
       ! freq is in MHz
       IF(PRESENT(freq) .AND. freq(i) .GT. 0.d0) THEN
          ades_obs(i)%frq = freq(i)
          ades_flg(i)%frq = .TRUE.
          WRITE(data_str,*) ades_obs(i)%frq
          ades_prec(i)%frq = get_prec(data_str)
       ENDIF

       ! fill com
       IF (trim(adjustl(observ(i)%tech)) .EQ. 'c') THEN
          ades_obs(i)%com = 1
          ades_flg(i)%com = .TRUE.
       ENDIF
    ENDIF

    ! fill Observation Group 
    ! fill RA and Dec if observation is optical, raStar, decStar and deltaRA, deltaDec in occultation case
    IF (ades_flg(i)%opt) THEN
       ! ra
       ades_obs(i)%ra = observ(i)%coord(1)*degrad
       ades_flg(i)%ra = .TRUE.
       WRITE(data_str,*) ades_obs(i)%ra
       ades_prec(i)%ra = get_prec(data_str)
       ! dec
       ades_obs(i)%dec = observ(i)%coord(2)*degrad
       ades_flg(i)%dec = .TRUE.
       WRITE(data_str,*) ades_obs(i)%dec
       ades_prec(i)%dec = get_prec(data_str)
    ELSEIF (ades_flg(i)%occ) THEN
       ! raStar
       ades_obs(i)%raStar = observ(i)%coord(1)*degrad
       ades_flg(i)%raStar = .TRUE.
       ades_prec(i)%raStar = ades_prec(i)%ra
       ! decStar
       ades_obs(i)%decStar =  observ(i)%coord(2)*degrad
       ades_flg(i)%decStar = .TRUE.
       ades_prec(i)%decStar = ades_prec(i)%dec
       ! deltaRA
       ades_obs(i)%deltaRA = 0.d0
       ades_flg(i)%deltaRA = .TRUE.
       ! deltaDec
       ades_obs(i)%deltaDec =  0.d0
       ades_flg(i)%deltaDec = .TRUE.
    ELSEIF (ades_flg(i)%rad) THEN
       IF (trim(adjustl(observ(i)%type)) .EQ. 'R') THEN
          ! delay
          ades_obs(i)%delay = observ(i)%coord(1)*8.64d4*2.d0/vlight
          ades_flg(i)%delay = .TRUE.
          WRITE(data_str,*) ades_obs(i)%delay
          ades_prec(i)%delay = get_prec(data_str)
          ! rmsDelay
          ades_obs(i)%rmsDelay = wei_bia(i)%rms_coord(1)*8.64d4*2.d0*1.d6/vlight
          ades_flg(i)%rmsDelay = .TRUE.
          WRITE(data_str,*) ades_obs(i)%rmsDelay
          ades_prec(i)%rmsDelay = get_prec(data_str)
       ELSEIF(trim(adjustl(observ(i)%type)) .EQ. 'V') THEN
          ! doppler
          ades_obs(i)%doppler = - observ(i)%coord(2)*ades_obs(i)%frq*1.d6*2.d0/vlight
          ades_flg(i)%doppler = .TRUE.
          WRITE(data_str,*) ades_obs(i)%doppler
          ades_prec(i)%doppler = get_prec(data_str)
          ! rmsDoppler
          ades_obs(i)%rmsDoppler = wei_bia(i)%rms_coord(2)*ades_obs(i)%frq/(0.5d0*1.d-6*vlight)
          ades_flg(i)%rmsDoppler = .TRUE.
          WRITE(data_str,*) ades_obs(i)%rmsDoppler
          ades_prec(i)%rmsDoppler = get_prec(data_str)
       ENDIF
    ENDIF

    ! fill astCat
    IF (ades_flg(i)%opt .OR. ades_flg(i)%occ) THEN
       CALL get_astCat(observ(i)%catcodmpc,ades_obs(i)%astCat)
       ades_flg(i)%astCat = .TRUE.
    ENDIF

    ! fill Photometry Group
    IF (.NOT.ades_flg(i)%rad) THEN
       ades_flg(i)%mag = .TRUE.
       IF (observ(i)%mag_def) THEN
          ades_obs(i)%mag = observ(i)%mag
          WRITE(data_str,*) ades_obs(i)%mag
          ades_prec(i)%mag = get_prec(data_str)
          ades_obs(i)%photCat = ades_obs(i)%astCat
       !ELSE
          !ades_obs(i)%mag = 9.9d9
       ENDIF
       ades_obs(i)%band = observ(i)%mag_band
       ades_flg(i)%band = .TRUE.
       ades_flg(i)%photCat = .TRUE.
    ENDIF

    ! fill disc and Precision Group 
    IF (.NOT.ades_flg(i)%rad) THEN
       ! fill disc 
       IF(PRESENT(disc) .AND. disc(i) .NE. "") THEN
          ades_flg(i)%disc = .TRUE.
          ades_obs(i)%disc = disc(i)
       ENDIF
       ! fill subfmt
       IF(PRESENT(subfmt) .AND. subfmt(i) .NE. "") THEN
          ades_flg(i)%subFmt = .TRUE.
          ades_obs(i)%subFmt = subfmt(i)
       ENDIF
       
       ! fill Precision Group
       ades_obs(i)%precTime = FLOOR(observ(i)%acc_time*1.d6)
       ades_flg(i)%precTime = .TRUE.
       ades_obs(i)%precRA = observ(i)%acc_coord(1)*cosd*secrad/15.d0 ! RA accuracy converted in seconds: 1 second is equivalent to 15 arcsecs
       ades_flg(i)%precRA = .TRUE. 
       ades_obs(i)%precDec = observ(i)%acc_coord(2)*secrad          ! Dec accuracy converted in arcsec
       ades_flg(i)%precDec = .TRUE.
       ! control on required format
       CALL precision_group_control(ades_obs(i)%precRA,data_str)
       ades_prec(i)%precRA = get_prec(data_str)
       CALL precision_group_control(ades_obs(i)%precDec,data_str)
       ades_prec(i)%precDec = get_prec(data_str)
    ENDIF

    !fill Residuals Group

    ! fill orbProd
    ades_flg(i)%orbProd = .TRUE.
    ades_obs(i)%orbProd = "OrbFit5.0"

    ! optical obs
    IF (.NOT. ades_flg(i)%rad) THEN

       ! fill ast_flag and mag_flag
       WRITE(ades_fit(i)%ast_flag,'(i1)') wei_bia(i)%sel_coord
       WRITE(ades_fit(i)%mag_flag,'(i1)') wei_bia(i)%sel_mag
       
       ! fill resRA and resDec (convert to arcsec and apply COS(Dec) factor)
       ades_flg(i)%resRA = .TRUE.
       ades_obs(i)%resRA = wei_bia(i)%res_coord(1)*secrad*cosd                                                            
       ades_flg(i)%resDec = .TRUE.
       ades_obs(i)%resDec = wei_bia(i)%res_coord(2)*secrad
       ades_fit(i)%ades_resc_def = wei_bia(i)%resc_def
                          
       ! fill selAst
       ades_flg(i)%selAst = .TRUE.
       IF (wei_bia(i)%sel_coord .EQ. 0) THEN
          ades_obs(i)%selAst = "D"
       ELSE
          ades_obs(i)%selAst = "A" 
       ENDIF

       ! fill sigRA, sigDec, sigCorr 
       ades_flg(i)%sigRA = .TRUE.
       ades_obs(i)%sigRA = wei_bia(i)%rms_coord(1)*secrad*cosd
       WRITE(data_str,*) ades_obs(i)%sigRA
       ades_prec(i)%sigRA = get_prec(data_str)
       IF (ades_obs(i)%sigRA .LT. 0.d0) THEN
          WRITE(ierrou,*) 'ERROR: convert_fitobs_ades: computed sigRA is negative '//err_str
          numerr = numerr+1
       ENDIF
       ades_flg(i)%sigDec = .TRUE.
       ades_obs(i)%sigDec = wei_bia(i)%rms_coord(2)*secrad
       WRITE(data_str,*) ades_obs(i)%sigDec
       ades_prec(i)%sigDec = get_prec(data_str)
       IF (ades_obs(i)%sigDec .LT. 0.d0) THEN
          WRITE(ierrou,*) 'ERROR: convert_fitobs_ades: computed sigDec is negative '//err_str
          numerr = numerr+1
       ENDIF
       ades_flg(i)%sigCorr = .TRUE.
       ades_obs(i)%sigCorr = wei_bia(i)%correl
       WRITE(data_str,*) ades_obs(i)%sigCorr
       ades_prec(i)%sigCorr = get_prec(data_str)
       IF (ades_obs(i)%sigCorr .GT. 1.d0 .OR. ades_obs(i)%sigCorr .LT. -1.d0) THEN
          WRITE(ierrou,*) 'ERROR: convert_fitobs_ades: computed sigCorr is not in [-1,1] interval'//err_str
          numerr = numerr+1
       ENDIF
       
       ! sigTime and biasTime obs stored with their default values  
       ades_flg(i)%sigTime = .TRUE.
       ades_flg(i)%biasTime= .TRUE.

       ! fill biasRA and biasDec
       ades_flg(i)%biasRA = .TRUE.
       ades_obs(i)%biasRA = wei_bia(i)%bias_coord(1)*secrad*cosd
       ades_flg(i)%biasDec = .TRUE.
       ades_obs(i)%biasDec = wei_bia(i)%bias_coord(2)*secrad

       ! fill photometry subgroup
       ! fill photProd 
       ades_obs(i)%photProd = ades_obs(i)%orbProd
       ades_flg(i)%photProd = .TRUE.

       ! fill resMag
       ades_flg(i)%resMag = .TRUE.
       ades_obs(i)%resMag = wei_bia(i)%res_mag
       ades_fit(i)%ades_resm_def = wei_bia(i)%resm_def

       ! fill selPhot
       ades_flg(i)%selPhot = .TRUE.
       IF (wei_bia(i)%sel_mag .EQ. 0) THEN
          ades_obs(i)%selPhot = "D"
       ELSE
          ades_obs(i)%selPhot = "A" 
       ENDIF

       ! sigMag 
       IF (wei_bia(i)%rms_mag .GE. 0) ades_obs(i)%sigMag = wei_bia(i)%rms_mag
       ades_flg(i)%sigMag = .TRUE.
       WRITE(data_str,*) ades_obs(i)%sigMag
       ades_prec(i)%sigMag = get_prec(data_str)
       IF (ades_obs(i)%sigMag .LT. 0.d0) THEN
          WRITE(ierrou,*) 'ERROR: convert_fitobs_ades: computed sigMag is negative '//err_str
          numerr = numerr+1
       ENDIF

       ! biasMag
       ades_flg(i)%biasMag = .TRUE.
       ! default is zero
       ades_obs(i)%biasMag = 0.d0
       ! assign bias corresponding to band (values for z,i,g,r,y,w are provided by Fitzsimmons November 2011 report)
       IF(trim(adjustl(ades_obs(i)%band)) .EQ. 'V' .OR. trim(adjustl(ades_obs(i)%band)) .EQ. 'H') THEN
          ades_obs(i)%biasMag = +0.00d0
       ELSE IF(trim(adjustl(ades_obs(i)%band)) .EQ. ' ' .OR. trim(adjustl(ades_obs(i)%band)) .EQ. 'B') THEN
          ades_obs(i)%biasMag = -0.80d0
       ELSE IF(trim(adjustl(ades_obs(i)%band)) .EQ. 'I') THEN
          ades_obs(i)%biasMag = +0.80d0
       ELSE IF(trim(adjustl(ades_obs(i)%band)) .EQ. 'R' .OR. trim(adjustl(ades_obs(i)%band)) .EQ. 'C') THEN
          ades_obs(i)%biasMag = +0.40d0
       ELSE IF(trim(adjustl(ades_obs(i)%band)) .EQ. 'z') THEN
          ades_obs(i)%biasMag = +0.37d0
       ELSE IF(trim(adjustl(ades_obs(i)%band)) .EQ. 'i') THEN
          ades_obs(i)%biasMag = +0.39d0
       ELSE IF(trim(adjustl(ades_obs(i)%band)) .EQ. 'g') THEN
          ades_obs(i)%biasMag = -0.28d0
       ELSE IF(trim(adjustl(ades_obs(i)%band)) .EQ. 'U' .OR. trim(adjustl(ades_obs(i)%band)) .EQ. 'J') THEN
          ades_obs(i)%biasMag = +1.10d0
       ELSE IF(trim(adjustl(ades_obs(i)%band)) .EQ. 'r') THEN
          ades_obs(i)%biasMag = +0.23d0
       ELSE IF(trim(adjustl(ades_obs(i)%band)) .EQ. 'y') THEN
          ades_obs(i)%biasMag = +0.36d0
       ELSE IF(trim(adjustl(ades_obs(i)%band)) .EQ. 'w') THEN
          ades_obs(i)%biasMag = +0.16d0
       ELSE IF(trim(adjustl(ades_obs(i)%band)) .EQ. 'G') THEN
          ades_obs(i)%biasMag = +0.28d0                       ! Gaia correction by Tholen
       ELSE IF(trim(adjustl(ades_obs(i)%band)) .EQ. 'c') THEN
          ades_obs(i)%biasMag = -0.05d0
       ELSE IF(trim(adjustl(ades_obs(i)%band)) .EQ. 'o') THEN
          ades_obs(i)%biasMag = +0.33d0
       ELSE IF(trim(adjustl(ades_obs(i)%band)) .EQ. 'L') THEN
          ades_obs(i)%biasMag =  +0.20d0
       ELSE IF(trim(adjustl(ades_obs(i)%band)) .EQ. 'K') THEN
          ades_obs(i)%biasMag = +1.70d0
       ELSE IF(trim(adjustl(ades_obs(i)%band)) .EQ. 'Y') THEN
          ades_obs(i)%biasMag = +0.70d0
       ELSE IF(trim(adjustl(ades_obs(i)%band)) .EQ. 'u') THEN
          ades_obs(i)%biasMag = +2.50d0
       ELSE
          WRITE(ierrou,*) 'ERROR: convert_fitobs_ades: unknown color, unable to find corresponding biasMag '//err_str 
          numerr=numerr+1
       END IF

       ! set photometric model photMod
       ades_flg(i)%photMod = .TRUE.
       ades_obs(i)%photMod = "H"

       ! radar obs
    ELSE

       ! fill radar_flag
       WRITE(ades_fit(i)%radar_flag,'(i1)') wei_bia(i)%sel_coord

       IF (ades_flg(i)%delay) THEN

          ! fill resDelay
          ades_flg(i)%resDelay = .TRUE.
          ades_obs(i)%resDelay = 2.d0*wei_bia(i)%res_coord(1)*1.d6*8.64d4/vlight
          ades_fit(i)%ades_resc_def = wei_bia(i)%resc_def

          ! fill selDelay
          ades_flg(i)%selDelay = .TRUE.
          IF (wei_bia(i)%sel_coord .EQ. 0) THEN
             ades_obs(i)%selDelay = "D"
          ELSE
             ades_obs(i)%selDelay = "A" 
          ENDIF

          ! fill sigDelay 
          ades_flg(i)%sigDelay = .TRUE.
          ades_obs(i)%sigDelay = 2.d0*1.d6*wei_bia(i)%rms_coord(1)*8.64d4/vlight
          WRITE(data_str,*) ades_obs(i)%sigDelay
          ades_prec(i)%sigDelay = get_prec(data_str)
          IF (ades_obs(i)%sigDelay .LE. 0.d0) THEN
             WRITE(ierrou,*) 'ERROR: convert_fitobs_ades: computed sigDelay is negative  '//err_str
             numerr = numerr+1
          ENDIF

       ELSE
          ! fill resDoppler
          ades_flg(i)%resDoppler = .TRUE.
          ades_obs(i)%resDoppler = 2.d6*wei_bia(i)%res_coord(2)*ades_obs(i)%frq/vlight
          ades_fit(i)%ades_resc_def = wei_bia(i)%resc_def

          ! fill selDoppler
          ades_flg(i)%selDoppler = .TRUE.
          IF (wei_bia(i)%sel_coord .EQ. 0) THEN
             ades_obs(i)%selDoppler = "D"
          ELSE
             ades_obs(i)%selDoppler = "A" 
          ENDIF

          ! fill sigDoppler 
          ades_flg(i)%sigDoppler = .TRUE.
          ades_obs(i)%sigDoppler = wei_bia(i)%rms_coord(2)*ades_obs(i)%frq*2.d0*1.d6/vlight
          WRITE(data_str,*) ades_obs(i)%sigDoppler
          ades_prec(i)%sigDoppler = get_prec(data_str)
          IF (ades_obs(i)%sigDoppler .LE. 0.d0) THEN
             WRITE(ierrou,*) 'ERROR: convert_fitobs_ades: computed sigDoppler is negative '//err_str
             numerr = numerr+1
          ENDIF
       ENDIF
    ENDIF

    ! fill chi
    IF (wei_bia(i)%resc_def) ades_fit(i)%chi = wei_bia(i)%chi
    IF (ades_fit(i)%chi .LT. 0.d0) THEN
       WRITE(ierrou,*) 'ERROR: convert_fitobs_ades: computed chi is negative '//err_str
       numerr = numerr+1
    ENDIF

    ! fill force_w_flags
    ades_fit(i)%force_w_flags(1) = wei_bia(i)%force_w(1)
    ades_fit(i)%force_w_flags(2) = wei_bia(i)%force_w(2)
    
  END SUBROUTINE new_ades_obs

  ! --------------------------------------------------------------------- !
  ! SUBROUTINE isnumber                                                   !
  ! It tells whether a character string contains only digits              !
  ! WARNING: assumes ASCII internal representation of characters          !
  !                                                                       !
  !  0 nul   16 dle   32 sp    48 0     64 @     80 P     96 `    112 p   !
  !  1 soh   17 dc1   33 !     49 1     65 A     81 Q     97 a    113 q   !
  !  2 stx   18 dc2   34 "     50 2     66 B     82 R     98 b    114 r   !
  !  3 etx   19 dc3   35 #     51 3     67 C     83 S     99 c    115 s   !
  !  4 eot   20 dc4   36 $     52 4     68 D     84 T    100 d    116 t   !
  !  5 enq   21 nak   37 %     53 5     69 E     85 U    101 e    117 u   !
  !  6 ack   22 syn   38 &     54 6     70 F     86 V    102 f    118 v   !
  !  7 bel   23 etb   39 '     55 7     71 G     87 W    103 g    119 w   !
  !  8 bs    24 can   40 (     56 8     72 H     88 X    104 h    120 x   !
  !  9 ht    25 em    41 )     57 9     73 I     89 Y    105 i    121 y   !
  ! 10 nl    26 sub   42 *     58 :     74 J     90 Z    106 j    122 z   !
  ! 11 vt    27 esc   43 +     59 ;     75 K     91 [    107 k    123 {   !
  ! 12 np    28 fs    44 ,     60 <     76 L     92 \    108 l    124 |   !
  ! 13 cr    29 gs    45 -     61 =     77 M     93 ]    109 m    125 }   !
  ! 14 so    30 rs    46 .     62 >     78 N     94 ^    110 n    126 ~   !
  ! 15 si    31 us    47 /     63 ?     79 O     95 _    111 o    127 del !
  ! --------------------------------------------------------------------- !
  ! Author: D. Bracali Cioci                                              !
  ! Last Update: 2013 Sep                                                 !
  ! --------------------------------------------------------------------- !
  LOGICAL FUNCTION isnumber(string) 

    CHARACTER(LEN=*),INTENT(IN) :: string   ! input string

    INTEGER :: i,icc 

    isnumber = .FALSE. 
    DO i=1,LEN(string) 
       icc = ICHAR(string(i:i)) 
       IF(icc.LT.48.OR.icc.GT.57) RETURN 
    END DO

    isnumber = .TRUE.           
  END FUNCTION isnumber

  ! --------------------------------------------------- !
  ! SUBROUTINE stat_min_rms                             !
  ! Read file neocp.rules containing minimum credible   !
  ! uncertainties for observations of the most          !
  ! productive stations (to be used in association with !
  ! ADES RMS for weighting observations)                !
  ! --------------------------------------------------- !
  ! Author: A. Bertolucci                               !
  ! Last Update: Nov 2020                               !
  ! --------------------------------------------------- !

  SUBROUTINE stat_min_rms(n,obs,min_wei)

    INTEGER,          INTENT(IN)                    :: n
    TYPE(ast_obs),    INTENT(IN),  DIMENSION(n)     :: obs
    REAL(KIND=dkind), INTENT(OUT), DIMENSION(2,n)   :: min_wei
  
    INTEGER          :: j,k
    INTEGER          :: unit ! for input of RMS values
    INTEGER          :: n_rules
    CHARACTER(LEN=3) :: obscatcod(nrulesx), observ, obscod
    REAL(KIND=dkind) :: rmsa, rmsd
    LOGICAL          :: code_found
    
    CALL filopl(unit,'weights/neocp.rules')
    READ(unit,*) ! skip header
    n_rules = 0
    DO j = 1, nrulesx
       READ(unit, 500, END=501) obscod,rmsbin(1:2,j)
       obscatcod(j) = obscod
500    FORMAT(A3,2X,F5.2,1X,F5.2)
       n_rules = n_rules + 1
    ENDDO
501 CALL filclo(unit,' ')
    DO k = 1, n
       observ = obs(k)%obscod_s
       code_found = .FALSE.
       ! selection of the record containing the obscode
       DO j = 1, n_rules
          IF (INDEX(obscatcod(j),observ) .NE. 0) THEN 
             rmsa = rmsbin(1,j)
             rmsd = rmsbin(2,j)
             code_found = .TRUE.
             EXIT
          ENDIF
       ENDDO
       IF (.NOT. code_found) THEN ! If the record matching obs_cod was not found, assign default value
          rmsa = 0.d0
          rmsd = 0.d0
          WRITE(iun_log,*) 'convert_ades_fitobs: station code ',observ,' not found in file neocp.rules. & 
               & Using error model values for weighting.'
       ENDIF
       ! store values
       min_wei(1,k) = rmsa*radsec/COS(obs(k)%coord(2))
       min_wei(2,k) = rmsd*radsec
    END DO
    
  END SUBROUTINE stat_min_rms

END MODULE astrometric_observations

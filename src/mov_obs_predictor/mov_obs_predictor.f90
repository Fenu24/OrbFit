!=============================================================!
! PROGRAM mov_obs_predictor                                   !
! It computes the prediction corresponding to a certain date  !
! and observatory for each orbit of MOV sampling              !
!=============================================================!
! Author: D. Bracali Cioci                                    !
! Last Update: Mar 2020                                       !
!=============================================================!
PROGRAM mov_obs_predictor

  USE fund_const
  USE output_control
  USE name_rules
  USE dyn_param
  USE propag_state, ONLY: pro_ele
  USE close_app,    ONLY: kill_propag
  USE orbit_elements
  USE pred_obs
  USE mov_obs_predictor_util

  IMPLICIT NONE

  TYPE(orbit_elem)   :: el_orb           ! orbital elements
  CHARACTER(LEN=3)   :: coo              ! coordinates of MOV sampling
  CHARACTER(LEN=6)   :: progna           ! program name
  CHARACTER(LEN=15)  :: obj_name         ! object name in .mov_sample
  CHARACTER(LEN=80)  :: run              ! run name
  CHARACTER(LEN=100) :: filnam           ! name of generic file
  REAL(KIND=dkind)   :: epoch            ! epoch of initial orbit
  REAL(KIND=dkind)   :: chi              ! chi of MOV fit
  REAL(KIND=dkind)   :: rms              ! rms of MOV fit
  REAL(KIND=dkind)   :: alpha,delta      ! output of prediction
  REAL(KIND=dkind)   :: vis_mag          ! visual magnitude
  REAL(KIND=dkind)   :: alpha_chimin     ! ra of prediction of orbit with minimum chi
  REAL(KIND=dkind)   :: delta_chimin     ! dec of prediction of orbit with minimum chi
  REAL(KIND=dkind)   :: d_alpha,d_delta  ! difference between ra and dec of current
  !                                        prediction and prediction of orbit with minimum chi
  INTEGER            :: cc_flag          ! flag of connected component
  INTEGER            :: imp_flag         ! flag for impact
  INTEGER            :: len_filna        ! length of filnam
  INTEGER            :: len1,len2        ! auxiliary lengths
  INTEGER            :: idone            ! unit of file .done
  INTEGER            :: iun_orb          ! unit of file .eq0
  INTEGER            :: iun_movsample    ! unit of file .mov_sample
  INTEGER            :: iun_movpred      ! unit of file .mov_pred
  INTEGER            :: iun_nompred      ! unit of file .nom_pred
  INTEGER            :: i_ind,j_ind      ! indices of sampling
  INTEGER            :: inl              ! flag for nonlinearity
  INTEGER            :: imp_orb          ! flag for impacting orbits
  INTEGER            :: geoc_orb         ! flag for geocentric orbits
  REAL(KIND=dkind)   :: moid             ! MOID
  REAL(KIND=dkind)   :: v_infty          ! velocity at infinity
  LOGICAL            :: eof              ! flag of end of file
  LOGICAL            :: first_pred       ! flag of first prediction
  REAL(KIND=dkind)   :: princ            ! principal value
  INTEGER            :: deg              ! DEC degrees
  INTEGER            :: hour             ! hour for RA
  INTEGER            :: min              ! minutes for RA
  REAL(KIND=dkind)   :: sec              ! seconds for RA
  INTEGER            :: ln               ! auxiliary length
  CHARACTER(LEN=19)  :: str              ! auxiliary string
  LOGICAL            :: fail             ! auxiliary flag
  CHARACTER(LEN=12)  :: alpha_str        ! string for RA output
  CHARACTER(LEN=12)  :: delta_str        ! string for DEC output
  CHARACTER(LEN=1)   :: signo            ! signo of RA
  REAL(KIND=dkind)   :: appm_a           ! apparent motion in RA
  REAL(KIND=dkind)   :: appm_d           ! apparent motion in DEC
  REAL(KIND=dkind)   :: appm             ! apparent motion modulus
  REAL(KIND=dkind)   :: pa               ! position angle
  REAL(KIND=dkind)   :: airmass          ! airmass
  REAL(KIND=dkind)   :: cosaux           ! auxiliary variable
  REAL(KIND=dkind)   :: adot             ! time derivative of RA
  REAL(KIND=dkind)   :: ddot             ! time derivative of DEC
  REAL(KIND=dkind)   :: pha              ! phase angle
  REAL(KIND=dkind)   :: dis              ! Earth distance
  REAL(KIND=dkind)   :: dsun             ! Sun distance
  REAL(KIND=dkind)   :: elo              ! Solar elongation
  REAL(KIND=dkind)   :: gallat           ! galactic latitude
  REAL(KIND=dkind)   :: gallon           ! galactic longitude
  REAL(KIND=dkind)   :: elmoon           ! Lunar elongation
  REAL(KIND=dkind)   :: elsun            ! Sun ele
  REAL(KIND=dkind)   :: elev             ! object elevation
  REAL(KIND=dkind)   :: cvf              ! auxiliary conversion factor

  ! Run name
  WRITE(*,*) 'Run name ='
  READ(*,100) run
100 FORMAT(A)

  IF(run.EQ.'') STOP 'No run specified.'
  ! input options
  progna = 'movpre'

  ! Output files for control, errors and warnings
  filnam = run//'.log'
  CALL rmsp(filnam,len_filna)
  CALL filopn(iun_log,filnam(1:len_filna),'UNKNOWN')
  filnam = run//'.err'
  CALL rmsp(filnam,len_filna)
  CALL filopn(ierrou,filnam(1:len_filna),'UNKNOWN')
  numerr = 0
  filnam = run//'.warn'
  CALL rmsp(filnam,len_filna)
  CALL filopn(iwarou,filnam(1:len_filna),'UNKNOWN')
  numwar = 0
  filnam = run//'.done'
  CALL rmsp(filnam,len_filna)
  CALL filopn(idone,filnam(1:len_filna),'UNKNOWN')

  ! Input options, for the propagator and the specific main program
  CALL mov_pred_opt(progna,run)
  CALL rmsp(name_obj,len2)

  !----------------------------------------------!
  ! Output file for MOV observational prediction !
  !----------------------------------------------!
  filnam = run//'.mov_pred'
  CALL rmsp(filnam,len_filna)
  CALL filopn(iun_movpred,filnam(1:len_filna),'UNKNOWN')
  ! Nominal observation prediction data
  filnam = run//'.nom_pred'
  CALL rmsp(filnam,len_filna)
  CALL filopn(iun_nompred,filnam(1:len_filna),'UNKNOWN')

  ! Open MOV sample file
  filnam = movdir(1:len_movdir)//name_obj(1:len2)//'.mov_sample'
  CALL rmsp(filnam,len_filna)
  CALL filopn(iun_movsample,filnam(1:len_filna),'UNKNOWN')

  ! Read header of file .mov_sample
  CALL read_movsample_header(iun_movsample,obj_name,epoch,coo)
  CALL rmsp(name_obj,len1)
  IF(name_obj(1:len1).NE.obj_name(1:len2))THEN
     WRITE(*,*) 'mov_obs_predictor: inconsistent object name in option file and .mov_sample'
     WRITE(ierrou,*) 'mov_obs_predictor: inconsistent object name in option file and .mov_sample'
     STOP
  END IF

  !----------------------------------------!
  ! PREDICTION COMPUTATION OF MOV SAMPLING !
  !----------------------------------------!
  WRITE(iun_movpred,'(A)') '% i    j   RA*COS(DELTA) [deg]      DELTA [deg]              V        '// &
       &                   'chi         q [au]                   e                         '//        &
       &                   'MOID [au]     V_inf [km/s]  imp  geoc'
  first_pred = .TRUE.
  inl = 0
  DO
     ! Read one record of file .mov_sample and propagate
     el_orb     = undefined_orbit_elem
     el_orb%coo = coo
     el_orb%t   = epoch
     CALL read_movsample_record(iun_movsample,i_ind,j_ind,el_orb,chi,rms,cc_flag, &
          &                     moid,v_infty,imp_orb,geoc_orb,eof)
     IF(eof) EXIT
     ! Prediction for points within the value set in the options file
     IF(chi.GT.op_sigma .AND. .NOT.first_pred) CYCLE
     !------------------------!
     ! Observation prediction !
     !------------------------!
     IF(first_pred)THEN
        CALL predic_obs(el_orb,id_obscode,t_pred,'O',alpha,delta,vis_mag,inl, &
             &    ADOT0=adot,DDOT0=ddot,PHA0=pha,DIS0=dis,DSUN0=dsun,ELO0=elo,GALLAT0=gallat,GALLON0=gallon,     &
             &    ELEV0=elev,ELSUN0=elsun,ELMOON0=elmoon)
        ! Orbit with minimum chi
        alpha = princ(alpha)
        alpha_chimin = alpha
        delta_chimin = delta
        !-------------------------------------------------!
        ! Output messages for the nominal obs. prediction !
        !-------------------------------------------------!
        ! Output RA
        CALL sessag(alpha_chimin*degrad/15.d0,signo,hour,min,sec)
        WRITE(str,'(F6.3)') sec
        CALL rmsp(str,ln)
        IF(ln.LT.6) str = '0'//str
        WRITE(alpha_str,200) hour,min,str
        ! Output DEC
        CALL sessag(delta_chimin*degrad,signo,deg,min,sec)
        WRITE(str,'(F5.2)') sec
        CALL rmsp(str,ln)
        IF(ln.LT.5) str = '0'//str
        WRITE(delta_str,210) signo,deg,min,str
        WRITE(*,'(A)') '::::::::::::::::::::::::::::::::::::'
        WRITE(*,'(A)') 'Nominal observation prediction'
        WRITE(*,'(A)') '::::::::::::::::::::::::::::::::::::'
        ! Airmass (Young, 1994)
        cosaux  = COS(pig/2-elev)
        airmass = 1.002432d0*cosaux**2+0.148386d0*cosaux+0.0096467d0
        airmass = airmass/(cosaux**3+0.149864d0*cosaux**2+0.0102963d0*cosaux+0.000303978d0)
        CALL angvcf('"/min',cvf,fail)
        ! Apparent motion
        appm_a = secrad*adot*COS(delta)/1440.d0
        appm_d = secrad*ddot/1440.d0
        appm   = SQRT(appm_a**2+appm_d**2)   ! apparent motion
        pa     = ATAN2(adot*COS(delta),ddot) ! position angle
        IF(pa.LT.0.d0) pa = pa+dpig
        pa=pa*degrad
        !
        WRITE(iun_nompred,220) cal_date(1:19),t_pred,id_obscode,op_sigma,alpha_str,alpha*degrad,delta_str,&
             &                 delta*degrad,appm_a,appm_d,appm,pa,vis_mag,elo*degrad,elmoon*degrad,       &
             &                 gallat*degrad,gallon*degrad
        IF(elev.GT.0.d0)THEN
           WRITE(iun_nompred,221) elev*degrad,airmass,pha*degrad
        ELSE
           WRITE(iun_nompred,222) elev*degrad,pha*degrad
        END IF
        !
        first_pred = .FALSE.
        READ(iun_movsample,*)
        CYCLE
     ELSE
        CALL predic_obs(el_orb,id_obscode,t_pred,'O',alpha,delta,vis_mag,inl)
        ! Relative data with respect to the orbit with minimum chi
        d_alpha = alpha-alpha_chimin+pig
        IF(d_alpha.LT.0) d_alpha = d_alpha+dpig
        d_alpha = (d_alpha-pig)*COS(delta_chimin)
        d_delta = delta-delta_chimin
     END IF
     WRITE(iun_movpred,230) i_ind,j_ind,d_alpha*degrad,d_delta*degrad,vis_mag,chi,el_orb%coord(1:2),moid,v_infty,imp_orb,geoc_orb
  END DO
  CALL filclo(iun_movpred,' ')


200 FORMAT(I2.2,':',I2.2,':',A6)
210 FORMAT(A1,I2.2,1X,I2.2,1X,A5)
220 FORMAT('Prediction time = ',A19,' UTC; ',F12.5,' MJD,   Observatory = ',I4.4,',   Sigma = ',F3.1/ &
         & ''/&
         & 'RA  = ',A12,' [HH:MM:SS]    ',F11.5,' [deg]'/    &
         & 'DEC = ',A12,' [deg min sec] ',F11.5,' [deg]'/ &
         & ''/&
         & 'RA/DEC Apparent motion = ',F9.3,' [arcsec/min]',F9.3,' [arcsec/min]'/&
         & 'Apparent motion = ',F9.3,' [arcsec/min], Position angle = ',F7.3,' [deg]'/&
         & 'Visual magnitude = ',F5.2,/&
         & 'Solar elongation = ',F7.2,' [deg], Lunar elongation = ',F7.2, &
         & ' [deg], Galactic latitude and longitude = ',F7.2,' [deg] ',F7.2,' [deg]')
221 FORMAT('Elevation = ',F7.2,' [deg], Airmass = ',F8.3,', Phase angle = ',F7.2,' [deg]')
222 FORMAT('Elevation = ',F7.2,' [deg], Airmass =    INF  , Phase angle = ',F7.2,' [deg]')
230 FORMAT(I3,2X,I3,1X,1P,2(E24.16,1X),0P,F9.4,1X,1P,E10.3,2(E25.16),0P,3X,F10.7,4X,F7.3,8X,I1,4X,I1)

END PROGRAM mov_obs_predictor

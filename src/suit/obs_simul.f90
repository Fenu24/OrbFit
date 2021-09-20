!==========================================================
! OBS_SIMUL 
! generation of observation record
! Copyright OrbFit consortium, 2005
! A.Milani, May 2005
! =========================================================
SUBROUTINE obs_simul(type1,t1,tut,astnam,ids,rmssec,rmsmag,alpha,delta,hmagn,obs,obsw)
  USE astrometric_observations
  USE name_rules
  USE station_coordinates, ONLY: codestat,obscoo
  USE reference_systems, ONLY: observer_position
  USE fund_const
  IMPLICIT NONE   
! ================input=============================
  CHARACTER*1, INTENT(IN) :: type1 ! O= optical R= radar range V= radar range rate
  DOUBLE PRECISION, INTENT(IN) :: t1, tut ! times, MJD, TDT and UTC
  CHARACTER*(name_len), INTENT(IN) :: astnam ! object designation (written in objdes field)
  INTEGER, INTENT(IN) :: ids ! station code (integer)
  DOUBLE PRECISION, INTENT(IN) :: rmssec, rmsmag ! RMS of random noise in astrometry (arcsec),
                                           ! in photometry (magnitudes)
  DOUBLE PRECISION, INTENT(IN) :: alpha, delta, hmagn ! observation: RA, DEC, app. magnitude
! ===============output=============================
  TYPE(ast_obs),INTENT(OUT) :: obs
  TYPE(ast_wbsr),INTENT(OUT) :: obsw
! ========end interface=========================
  CHARACTER(LEN=16)  :: stname ! station name
  DOUBLE PRECISION   :: rvnorm ! Gaussian generator
! ==============================================
  obs=undefined_ast_obs
  obs%objdes=astnam
  obs%type=type1  
! technology and station codes
  IF(type1.eq.'O')THEN
     obs%tech='C'
     obs%obscod_i=ids
     WRITE(obs%obscod_s(3:7),'(I4)')ids
     CALL codestat(ids,obs%obscod_s)
     IF(rhs.ne.1.and.rhs.ne.2)THEN
        WRITE(*,*)'obs_simul: rhs=', rhs
        STOP
     ELSE
        CALL observer_position(t1,obs%obspos,obs%obsvel,OBSCODE=ids)
     END IF
  ELSEIF(type1.eq.'R'.or.type1.eq.'V')THEN
     obs%tech='c'
     obs%obscod_i=ids*10000+ids
     WRITE(obs%obscod_s,'(I3,1X,I3)')ids,ids
     CALL obscoo(ids,obs%obspos,stname)
     CALL obscoo(ids,obs%obsvel,stname)
  ENDIF
! observation time 
  obs%time_utc=tut
  obs%time_tdt=t1
  obs%acc_time=0.000001d0
! photometry
  IF(type1.eq.'O')THEN
     obs%mag=hmagn
     obs%mag_def=.true.
     WRITE(obs%mag_str,110)hmagn
110  FORMAT(F4.1,' V')
     obs%acc_mag=1.e-2
     obs%mag_band='V'
  ENDIF 
! weights
  obsw=undefined_ast_wbsr
  obsw%force_w=.true.
  IF(type1.eq.'O')THEN
     obsw%rms_coord(1)=rmssec*radsec  ! noise in alpha, delta
     obsw%rms_coord(2)=rmssec*radsec
     obs%coord(2)=delta+rvnorm()*obsw%rms_coord(2)
     obs%coord(1)=alpha+rvnorm()*obsw%rms_coord(1)/cos(delta)
     obs%acc_coord=1e-7
     obsw%rms_mag=rmsmag ! noise in photometry
     obsw%sel_mag=1
     obs%mag=hmagn+rvnorm()*rmsmag
  ELSEIF(type1.eq.'R')THEN  
     obsw%rms_coord(1)=rmssec/aukm ! noise in range
     obs%acc_coord(1)=rmssec/aukm
     obs%coord(1)=alpha+rvnorm()*obsw%rms_coord(1)
  ELSEIF(type1.eq.'V')THEN
     obsw%rms_coord(2)=rmssec/aukm ! noise in range rate
     obs%acc_coord(2)=rmssec/aukm
     obs%coord(2)=delta+rvnorm()*obsw%rms_coord(2)
  ENDIF
END SUBROUTINE obs_simul

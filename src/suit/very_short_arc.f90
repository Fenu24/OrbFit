MODULE very_short_arc
USE name_rules
USE output_control
USE fund_const
USE astrometric_observations
USE attributable
IMPLICIT NONE
PRIVATE

INCLUDE 'par_sharc.h90'

TYPE vsa

  INTEGER NOBS ! number of observations
  TYPE(ast_obs) :: obs(nobs_shx) ! observations
  TYPE(ast_wbsr) :: obsw(nobs_shx) ! observation weights
  TYPE(attrib) :: att    ! attributable
  INTEGER nights ! number of nights of observations
  CHARACTER*19 :: name ! name, also in A=B form
END TYPE vsa

PUBLIC vsa, nobs_shx ! type, parameter

PUBLIC vsa_read, in_ons ! routines  
CONTAINS 

SUBROUTINE vsa_read(iunobs,error_model,timestep,shift,v,obscodi,qua, &
&   error,eof,iundup,iunlst,iunatt,iunrat,mmin,rwoout,obsdir)  
  INTEGER, INTENT(IN) :: iunobs ! input logical unit
  INTEGER, INTENT(IN) :: iundup,iunlst,iunatt,iunrat ! output units
  INTEGER, INTENT(IN) :: mmin ! min no obs to write the .rwo file
  DOUBLE PRECISION, INTENT(IN) :: timestep,shift ! for rounded times 
  TYPE(vsa), INTENT(OUT) :: v ! very short arc 
  INTEGER, INTENT(OUT) :: qua ! quality code 1-10
  INTEGER, INTENT(OUT) :: obscodi ! integer observatory code
  CHARACTER*(*), INTENT(IN) :: error_model ! error model file
  LOGICAL, INTENT(OUT) :: eof, error ! if error, try next; if eof, end
  LOGICAL, INTENT(IN) :: rwoout ! flag to control .rwo file output
  CHARACTER*(*), INTENT(IN), OPTIONAL :: obsdir ! directories for output .rwo 
     ! required is rwoout=.true.
  logical obs0 ! successful input flag
  INTEGER nights, quality ! nights counting, quality functions
  INTEGER m, minoss
  CHARACTER*(name_len) namshort
  LOGICAL change !always true, that is .rwo files not consulted
  DOUBLE PRECISION trou, arc
! =============================================
  IF(rwoout.and..not.PRESENT(obsdir))THEN
     WRITE(*,*)'vsa_read: error in calling sequence ', rwoout
     STOP
  ENDIF
  error=.false.
  IF(rwoout)THEN
     CALL in_ons_rwo(iunobs,iundup,error_model,obsdir,namshort,obs0 &
&      ,m,v%obs,v%obsw,obscodi,nobs_shx,iun_log,mmin)
  ELSE
     CALL in_ons(iunobs,error_model,namshort,obs0 &
&      ,m,v%obs,v%obsw,obscodi,nobs_shx,iun_log,mmin)
  ENDIF
  IF(.not.obs0)THEN 
! input observation file ended
     eof=.true.
     RETURN
  ENDIF
  eof=.false.
! quality WARNING: should we remove the degraded observations?
  v%name=' '
  v%name=namshort
  v%nobs=m
  v%nights=nights(v%nobs,v%obs,v%obsw)
  arc=v%obs(m)%time_utc-v%obs(1)%time_utc
  qua=quality(arc,m,v%nights,minoss)
  change=.true.
! list of ONS for input to other programs
  CALL obssta(iunlst,v%name,v%name,v%obs(1:m)%time_utc          &
       &          ,v%obs(1:m)%type,m,v%nights,change)
  IF(m.lt.mmin) RETURN
! compute attributable
  CALL attri_comp(v%nobs,v%obs,v%obsw,v%att,error)
  IF(error) RETURN
! copy into vsa record already done
! find rounded time:                                                    
  trou=timestep*nint((v%att%tdtobs-shift)/timestep)+shift 
! output .att and .rat files
  CALL wri_attri(iunatt,iunrat,v%name,v%att,trou)
END SUBROUTINE vsa_read
! ====================================================
! IN_ONS                  
! Input of astrometric observations from a file (MPC format)
! without output on .rwo file
! anti duplication check should be done in the calling module            
! ==================================================== 
SUBROUTINE in_ons(iunobs,error_model,name0,obs0&
&      ,m,obs,obsw,obscodi,nobshx,iunout,mmin) 
  INTEGER, INTENT(IN) :: nobshx ! max numbe rof observations
  INTEGER, INTENT(IN) :: iunobs, iunout ! input/output units
  INTEGER, INTENT(IN) :: mmin ! min no obs to write the .rwo file
  CHARACTER*(*), INTENT(IN) :: error_model ! error model
  logical, intent(out) :: obs0 ! successful input flag
  CHARACTER*(name_len), INTENT(OUT) :: name0 ! ons unique fingerprint
  integer, intent(out) ::  m ! observation number
! observations new data type
  TYPE(ast_obs),DIMENSION(nobshx), INTENT(OUT) :: obs
  TYPE(ast_wbsr),DIMENSION(nobshx), INTENT(OUT) :: obsw
  INTEGER obscodi !observatory numeric code
  INTEGER, DIMENSION(nobshx) :: time_ord(nobshx)
! end interface
  TYPE(ast_obs),DIMENSION(nobshx) :: obs_ns ! not sorted
  INTEGER j,i,k ! char. counters,loop indexes,duplicate counter
  LOGICAl eof,init
  CHARACTER*4 n_enc
  CHARACTER*(name_len) name1
! ====================================================                  
! Input of astrometric observations from a file (MPC format)            
! ==================================================== 
  CALL mpc_obs_input(obs0,obs_ns,m,UNIT=iunobs,EOF=eof,ONS_MODE=.true.)
  name0=obs(1)%objdes ! object name
! problem: is the check on max size of obs done in the routine?
! yes, but a more flexible error management would be desirable
  IF(m.eq.0)THEN 
     obs0=.false. 
     RETURN
  ELSEIF(m.gt.nobshx)THEN
     obs0=.false.
     WRITE(ierrou,*)obs(1)%objdes,' too many ',m,' observations'
     numerr=numerr+1
     RETURN
  ELSE 
     obs0=.true. 
  ENDIF
! sort times first
   CALL heapsort(obs_ns(1:m)%time_tdt,m,time_ord)
   DO k=1,m
      i=time_ord(k)
      obs(k)=obs_ns(i)
  ENDDO
! find weights for these; a priori RMS of astrometric observations
  init=.true. ! nothing is set of the wsbr structure
  CALL observ_rms(obs,error_model,init,obsw,m)
! observatory code
  obscodi=obs(1)%obscod_i
  DO j=2,m
    IF(obs(j)%obscod_i.ne.obscodi)THEN
       WRITE(ierrou,*)name0,' changing obscode ', obs(1)%obscod_s,' to ',obs(j)%obscod_s
       numerr=numerr+1
    ENDIF
 ENDDO
 name0=obs(1)%objdes 
END SUBROUTINE in_ons
! ====================================================
! IN_ONS_RWO                  
! Input of astrometric observations from a file (MPC format)
! output on .rwo file, with anti duplication check            
! ==================================================== 
SUBROUTINE in_ons_rwo(iunobs,iundup,error_model,obsdir,name0,obs0&
&      ,m,obs,obsw,obscodi,nobshx,iunout,mmin) 
  INTEGER, INTENT(IN) :: nobshx ! max numbe rof observations
  INTEGER, INTENT(IN) :: iunobs, iundup, iunout ! input/output units
  INTEGER, INTENT(IN) :: mmin ! min no obs to write the .rwo file
  CHARACTER*(*), INTENT(IN) :: error_model ! error model
  CHARACTER*(*), INTENT(IN) :: obsdir ! directories for data 
  logical, intent(out) :: obs0 ! successful input flag
  CHARACTER*(name_len), INTENT(OUT) :: name0 ! ons unique fingerprint
  integer, intent(out) ::  m ! observation number
! observations new data type
  TYPE(ast_obs),DIMENSION(nobshx), INTENT(OUT) :: obs
  TYPE(ast_wbsr),DIMENSION(nobshx), INTENT(OUT) :: obsw
  INTEGER obscodi !observatory numeric code
  INTEGER, DIMENSION(nobshx) :: time_ord(nobshx)
! end interface
  TYPE(ast_obs),DIMENSION(nobshx) :: obs_ns ! not sorted
  CHARACTER*80 rwofiout ! file names 
  INTEGER lrwout,j,i,k,dupcou! char. counters,loop indexes,duplicate counter
  LOGICAl eof,init,exi,first
  CHARACTER*4 n_enc
  CHARACTER*(name_len) name1
  data first /.true./
  SAVE first, dupcou
  IF(first)THEN
     dupcou=0
     first=.false.
  ENDIF
! ====================================================                  
! Input of astrometric observations from a file (MPC format)            
! ==================================================== 
  CALL mpc_obs_input(obs0,obs_ns,m,UNIT=iunobs,EOF=eof,ONS_MODE=.true.)
! problem: is the check on max size of obs done in the routine?
! yes, but a more flexible error management would be desirable
  IF(m.eq.0)THEN 
     obs0=.false. 
     RETURN
  ELSEIF(m.gt.nobshx)THEN
     obs0=.false.
     WRITE(ierrou,*)obs(1)%objdes,' too many ',m,' observations'
     numerr=numerr+1
     RETURN
  ELSE 
     obs0=.true. 
  ENDIF
! sort times first
   CALL heapsort(obs_ns(1:m)%time_tdt,m,time_ord)
   DO k=1,m
      i=time_ord(k)
      obs(k)=obs_ns(i)
  ENDDO
! find weights for these; a priori RMS of astrometric observations
  init=.true. ! nothing is set of the wsbr structure
  CALL observ_rms(obs,error_model,init,obsw,m)
! compute file name                                                     
2 name0=obs(1)%objdes 
  CALL fidinam(obsdir,name0,'rwo',rwofiout,lrwout)
  INQUIRE(file=rwofiout(1:lrwout),exist=exi) 
  IF(exi)THEN
     dupcou=dupcou+1
     WRITE(*,*)'in_ons: duplicate rwo file for ',name0,dupcou
     CALL base_10_to_64(15581440+dupcou+10000,n_enc)
     name1=n_enc//name0(5:9)
     WRITE(*,*)'in_ons: replaced with new name ', name1
     obs%objdes=name1
     WRITE(iundup,100)name0,name1
100  FORMAT(' duplicate ',a9,2x,a9)
     GOTO 2
  ENDIF
! output observations and weights file
  IF(m.ge.mmin)THEN
     CALL write_rwo(rwofiout(1:lrwout),obs,obsw,m,error_model)
  ENDIF
! observatory code
  obscodi=obs(1)%obscod_i
  DO j=2,m
    IF(obs(j)%obscod_i.ne.obscodi)THEN
       WRITE(ierrou,*)name0,' changing obscode ', obs(1)%obscod_s,' to ',obs(j)%obscod_s
       numerr=numerr+1
    ENDIF
 ENDDO
END SUBROUTINE in_ons_rwo

END MODULE very_short_arc

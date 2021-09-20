MODULE fitobs_mod

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: finopt, finobs, f_gauss

CONTAINS

  ! ===================================================================   
  ! FINOPT                                                                
  ! ===================================================================   
  ! input options, for the propagator and the specific main program       
  ! input: progna = program name (6 characters)                           
  !        run    = run identifier (80 characters, up to 76 non-blank)    
  SUBROUTINE finopt(progna,run,astna0,astnap,error_model,ades_type) 
    USE output_control
    USE tp_trace, ONLY: tpplane
    USE eval_risk, ONLY: givenmass, massgiven
    USE fund_const
    USE astrometric_observations, ONLY: radius, ons_name
    USE offlov_checktp, ONLY: shrinkea
    USE cobweb, ONLY: hmax,sigx,ndir,np,ndir2,np2,grid_lev_curve,propag_geoc_orbit
    USE dyn_param, ONLY : ngr_opt
    USE multi_store, ONLY : prob_sampl
    USE cla_store, ONLY : dmeacontr
    implicit none 
    character*6,INTENT(IN) :: progna 
    character*80,INTENT(IN) :: run 
    CHARACTER*(*),INTENT(OUT) :: astna0,astnap 
    CHARACTER*(20),INTENT(OUT) :: error_model ! error model file name
    INTEGER, INTENT(OUT), OPTIONAL :: ades_type  ! flag about format of ADES output file
    ! hidden output through bizset: controls for bizarre orbits
    DOUBLE PRECISION ecclim,samin,samax,phmin,ahmax,qmax 
    ! ==========END INTERFACE============================================   
    integer le
    character*90 file 
    LOGICAL ireq,found
    CHARACTER*60 comment 
    CHARACTER*100 filnam
    INTEGER iunout, iuncovar ! now only local names
    ! modification to input mass from physical observations
    INCLUDE 'parlib.h90' 
    CHARACTER*200 filmass
    LOGICAL ok
    INTEGER iunmass
    ! ===================================================================== 
    CALL initopt(progna,run,'fop')                                         
    ! ===================================================================== 
    ! Output files: for control and results, for covariance  
    ! === WARNING === 
    ! File fou to be opened here because rmodel could write in it
    filnam=run//'.fou' 
    CALL rmsp(filnam,le) 
    CALL filopn(iunout,filnam,'UNKNOWN')
    iun_log=iunout ! later iunout to be abolished
    filnam=run//'.fga' 
    CALL rmsp(filnam,le) 
    CALL filopn(iuncovar,filnam,'UNKNOWN')
    iun_covar=iuncovar
    ! ===================================================================== 
    ! read option for physical model and integration method                 
    CALL rmodel(1) 
    ! ===================================================================== 
    ! initializations for Gauss method                                      
    CALL iodini 
    ! ===================================================================== 
    ! Output files: for errors, close approaches, for propagator parameters 
    filnam=run//'.err' 
    CALL rmsp(filnam,le)
    CALL filopn(ierrou,filnam,'UNKNOWN') 
    numerr=0 
    filnam=run//'.clo' 
    CALL rmsp(filnam,le)
    CALL filopn(iuncla,filnam,'UNKNOWN') 
    numcla=0
    filnam=run//'.pro' 
    CALL rmsp(filnam,le)
    CALL filopn(ipirip,filnam,'UNKNOWN')
    ! =============================
    ! asteroid name(s)                                                      
    ! =============SOME ASTEROID NAME NEEDED===================
    ireq=.false.
    astna0=run
    comment='first arc asteroid name'
    CALL input_cha_opt(progna,'astna0',astna0,ireq,found,comment,iunout)
    CALL rmsp(astna0,le)
    ons_name=.false.
    comment='not a designation'
    CALL input_log_opt(progna,'ons_name',ons_name,ireq,found,comment,iunout)
    ireq=.false. 
    astnap=' '
    comment='second arc asteroid name'
    CALL input_cha_opt(progna,'astnap',astnap,ireq,found,comment,iunout)
    CALL rmsp(astnap,le)
    ! =============SELECTED ERROR MODEL===================   
    error_model=' '
    ireq=.false.
    comment='error model file'
    CALL input_cha_opt(progna,'error_model',error_model,ireq,found,comment,iunout)
    ! =============FOR COBWEB ============================
    comment='compute a grid for the level curves in case of cobweb'
    CALL input_log_opt(progna,'grid_lev_curve',grid_lev_curve,ireq,found,comment,iunout)
    comment='number of rho values at first iteration'
    CALL input_int_opt(progna,'cob_ndir',ndir,ireq,found,comment,iunout)
    comment='number of rho_dot values at first iteration'
    CALL input_int_opt(progna,'cob_np',np,ireq,found,comment,iunout)
    comment='number of rho values at second iteration (also for cobweb)'
    CALL input_int_opt(progna,'cob_ndir2',ndir2,ireq,found,comment,iunout)
    comment='number of rho_dot values at second iteration (also for cobweb)'
    CALL input_int_opt(progna,'cob_np2',np2,ireq,found,comment,iunout)
    comment='sigma max for cobweb'
    CALL input_rea_opt(progna,'cob_sigx',sigx,ireq,found,comment,iunout)
    comment='max absolute magnitude'
    CALL input_rea_opt(progna,'cob_hmax',hmax,ireq,found,comment,iunout)
    comment='propag geocentric orbits'
    CALL input_log_opt(progna,'propag_geoc_orbit',propag_geoc_orbit,ireq,found,comment,iunout)
    ! ======selection of target plane=====================
    !  tpplane is set to false in rmodel; true is default in fitobs
    tpplane=.true.
    comment='use b-plane as target plane' 
    CALL input_log_opt(progna,'tpplane',tpplane,ireq,found,comment,iunout)
    ! ======is mass given in a separate file?=====================
    massgiven=.false.
    comment='mass given in a separate file'
    CALL input_log_opt(progna,'massgiven',massgiven,ireq,found,comment,iunout)
    IF(massgiven)THEN
       ! get mass known from other sources, if available
       filmass=dlibd//'/'//astna0//'.mass'
       CALL rmsp(filmass,le)
       INQUIRE(file=filmass,exist=ok)
       IF(ok)THEN
          CALL filopn(iunmass,filmass(1:le),'old')
          READ(iunmass,*) givenmass, radius
          massgiven=.true.
          radius=radius/aukm ! convert from km to AU
          CALL filclo(iunmass,' ')
       ELSE
          WRITE(*,*) filmass(1:le),'not found, no mass data'
          massgiven=.false.
          STOP
       ENDIF
    ENDIF
    ! =============bizarre control=======================   
    ! ================limit eccentricity============================            
    ireq=.false. 
    ecclim=0.d0  ! leave at the bizset default
    comment='limit eccentricity'
    CALL input_rea_opt(progna,'ecclim',ecclim,ireq,found,comment,iunout)
    ! ================limit semimajor axis============================            
    ireq=.false. 
    samax=0.d0 !leave at the bizset default
    comment='max semimajor axis'
    CALL input_rea_opt(progna,'samax',samax,ireq,found,comment,iunout)
    ireq=.false. 
    samin=0.d0 !leave at the bizset default
    comment='min semimajor axis'
    CALL input_rea_opt(progna,'samin',samin,ireq,found,comment,iunout)
    ! ================limit perihelion, aphelion============================
    ireq=.false. 
    ahmax=0.d0 !leave at the bizset default
    comment='max aphelion'
    CALL input_rea_opt(progna,'ahmax',ahmax,ireq,found,comment,iunout)
    ireq=.false. 
    phmin=0.d0 !leave at the bizset default
    comment='min perihelion'
    CALL input_rea_opt(progna,'phmin',phmin,ireq,found,comment,iunout)
    qmax=0.d0 !leave at the bizset default
    comment='max perihelion'
    CALL input_rea_opt(progna,'qmax',qmax,ireq,found,comment,iunout)
    ! ================control on bizarre orbits================ 
    CALL bizset(ecclim,samin,samax,phmin,ahmax,qmax) 
    ! shrink of Earth cross section
    ireq=.false. 
    shrinkea=1.d0
    comment='shrink of Earth cross section'
    CALL input_rea_opt(progna,'shrinkea',shrinkea,ireq,found,comment,iunout)
    IF(shrinkea.GT.1.d0)THEN
       WRITE(*,*)'resintpopt: shrinkea must be smaller than 1, shrinkea=',shrinkea
       STOP
    ENDIF
    ! =================use of IP sampling/sigma_lov sampling===================
    ireq=.false.
    prob_sampl=.false.
    comment='use of constant steps in IP'
    CALL input_log_opt(progna,'prob_sampl',prob_sampl,ireq,found,comment,iunout)
    ! ================selection of smaller disk (for total ignore of close app)=====
    ireq=.false.  
    comment='max dist. for inclolinctp'
    !  dmeacontr=dmea
    CALL input_rea_opt(progna,'dmeacontr',dmeacontr,ireq,found,comment,iunout)
    ! ==============selection of type of ADES output file===================   
    IF(PRESENT(ades_type))THEN
       ades_type=0
       ireq=.false.
       comment='flag for ADES output file'
       CALL input_int_opt(progna,'ades_type',ades_type,ireq,found,comment,iunout)
    ENDIF
  END SUBROUTINE finopt

  ! ==================================================================    
  ! FINOBS                                                                
  ! ==================================================================    
  ! Arc observations input                                               
  subroutine finobs(progna,iar,astna0,obs0,nlef,m,obs,obsw,rwofi0  &
       &     ,error_model,ades_type) 
    USE astrometric_observations
    USE output_control
    implicit none 
    ! =======INPUT==================================                        
    ! ========= file names, i/o control ============ 
    character*(6), INTENT(IN) :: progna ! program name
    character*(*), INTENT(IN) :: astna0 ! asteroid name (9 characters)
    integer, INTENT(IN) :: iar ! flag for arcs 1-2 
    CHARACTER*20, INTENT(IN) ::  error_model ! weighing model
    INTEGER, INTENT(IN), OPTIONAL :: ades_type ! flag for type of output file
    ! =======OUTPUT================================== 
    ! ===== main output: observational data =========================== 
    integer, INTENT(IN) :: nlef ! observation number
    integer, INTENT(OUT) :: m ! observation number
    ! new data types
    TYPE(ast_obs),DIMENSION(nlef), INTENT(OUT) :: obs
    TYPE(ast_wbsr),DIMENSION(nlef), INTENT(OUT) :: obsw
    ! auxiliary output                     
    character*(60), INTENT(OUT) :: rwofi0  ! residuals and weights file name  
    logical obs0 ! successful input flag
    logical change ! change flag
    ! ============END INTERFACE========================== 
    CHARACTER*60 comment ! for input-options
    logical found,ireq ! input flags                            
    character*60 obsdir0 ! file names, , for observations   
    CHARACTER*18 nam0                                                          
    INTEGER lnam, le, lench  ! string length (after blank removal)  
    integer ll ! loop indexes                                              
    INCLUDE 'sysdep.h90' ! directory char  
    ! precedence to .obs file with respect to .rwo file: false for use in fitobs
    LOGICAL precob
    precob=.false.                                             
    ! ====================================================   
    obs0=.false. ! nothing found yet
    m=0 
    nam0=astna0
    ! find observations file name
    ireq=.false. 
    IF(iar.eq.1)THEN
       CALL rmsp(nam0,le)
       !     CALL rmsp(astna0,le)
       nam0=astna0(1:le) 
       comment='first object obs. file name' 
       CALL input_cha_opt(progna,'nam0',nam0,ireq,found,comment,iun_log)
    ELSEIF(iar.eq.2 )THEN
       CALL rmsp(nam0,le)
       nam0=astna0(1:le)
       comment='second object obs. file name' 
       CALL input_cha_opt(progna,'namp',nam0,ireq,found,comment,iun_log)
    ENDIF
    CALL rmsp(nam0,lnam)
    IF(lnam.eq.0) RETURN ! blank object does not exist
    ! find directory of observations    
    ireq=.false.
    obsdir0='obsdata' ! locally called obsdir0 even for second arc    
    if(iar.eq.1)then 
       comment='first object obs. directory name'
       CALL input_cha_opt(progna,'obsdir0',obsdir0,ireq,found,comment,iun_log)
    elseif(iar.eq.2)then 
       comment='second object obs. directory name'
       CALL input_cha_opt(progna,'obsdirp',obsdir0,ireq,found,comment,iun_log)
    endif
    CALL rmsp(obsdir0,le) 
    ! input data 
    IF(PRESENT(ades_type))THEN
       CALL input_obs(obsdir0,nam0,precob,error_model,obs0,obs,obsw,m,   &
            &    iun_log,change,ADES_TYPE=ades_type)
    ELSE
       CALL input_obs(obsdir0,nam0,precob,error_model,obs0,obs,obsw,m,   &
            &    iun_log,change)
    ENDIF
    ! compute name of .rwo file                                             
    ll=lench(obsdir0) 
    rwofi0=obsdir0(1:ll)//dircha 
    ll=lench(rwofi0) 
    rwofi0=rwofi0(1:ll)//nam0//'.rwo' 
    CALL rmsp(rwofi0,ll) 
  END SUBROUTINE finobs

  ! Copyright (C) 1997 by Mario Carpino (carpino@brera.mi.astro.it)       
  ! Version OFINOD: December 1, 1997                                      
  ! Adapted for FITOBS, AM/ZK March 1998                                  
  ! Adapted for use of VAISALA, ZK November 1998                          
  ! IOBS added (Feb 10, 1999) MC   
  ! Fortran 90 version 3.0 A. Milani November 2002   
  ! --------------------------------------------------------------------- 
  !                                                                       
  !  *****************************************************************    
  !  *                                                               *    
  !  *                         F G A U S S                           *    
  !  *                                                               *    
  !  *                Auxiliary routine for FITOBS:                  *    
  !  *                 initial orbit determination                   *    
  !  *                                                               *    
  !  *****************************************************************    
  !                                                                       
  ! In the present version, only Gauss' method is supported               
  !                                                                       
  ! INPUT:    UNIELE    -  FORTRAN unit for elements                      
  !           NAME      -  Object name                                    
  !           DEFORB    -  Orbit definition flag                          
  !           DEFCOV    -  Orbit covariance definition flag               
  !           RWFIL     -  Names of residual/weight files
  !           OBS, OBSW -  Input observations and weights                 
  !           OBSW%SEL_COORD  -  Selection index (0=don't use; 1=use for fit;   
  !                        2=use for fit & Gauss method)                  
  !           N         -  Number of observations for each object         
  !           IMETH     -  Method to be used: 1=Gauss, 2=Vaisala          
  !                                                                       
  ! OUTPUT:   EL        -  Orbital elements                    
  ! NOT OUTPUT, but available: residulas inside obsw
  !                                                                       
  ! Input variables DEFORB, DEFCOV, OBSW%SEL_COORD are modified by the routine
  !                                                                       
  SUBROUTINE f_gauss(uniele,name,deforb,defcov,          &
       &     rwfil,obs,obsw,n,error_model,imeth,el,ades_type)
    USE fund_const, ONLY: gms 
    USE astrometric_observations 
    USE least_squares, ONLY: rms_compute
    USE orbit_elements
    USE output_control
    USE iodet
    IMPLICIT NONE 
    ! INPUT observations: new data types
    INTEGER, INTENT(IN) :: n  !number of observations
    TYPE(ast_obs),DIMENSION(n),INTENT(IN) :: obs
    TYPE(ast_wbsr),DIMENSION(n),INTENT(INOUT) :: obsw
    CHARACTER*20, INTENT(IN) ::  error_model ! weighing model
    ! other input   
    INTEGER,INTENT(IN) :: uniele, imeth ! unit for output, method
    CHARACTER*18, INTENT(IN) :: name !asteroid name 
    CHARACTER*60, INTENT(IN) :: rwfil  ! output file for residuals
    ! OUTPUT
    LOGICAL, INTENT(OUT) :: deforb,defcov ! orbit state on output 
    TYPE(orbit_elem), INTENT(OUT) :: el ! elements
    INTEGER, INTENT(IN), OPTIONAL :: ades_type ! flag for type of output file
    ! END INTERFACE

    ! HEADERS (to be revised)
    INCLUDE 'parcmc.h90'
    INCLUDE 'pariod.h90' ! initial orbit determination parameters
    INCLUDE 'comiod.h90' ! initial orbit determination common
    !
    INTEGER ln,j, nused
    INTEGER, EXTERNAL :: lench 
    LOGICAL fail 
    CHARACTER*3 eltype 
    CHARACTER*60 comele 
    DOUBLE PRECISION enne,rms 
    DOUBLE PRECISION, EXTERNAL :: snormd 
    DOUBLE PRECISION elem(6),h,g,gg(6,6),cc(6,6) 
    DOUBLE PRECISION eq(6), telem
    ln=lench(name) 
    IF(imeth.eq.1)THEN 
       WRITE(iun_log,203) name(1:ln) 
203    FORMAT('Automatic method for object ',A,':') 
       iodnm=2 
       iodmet(1)=1 
       iodmen(1)='GAUSS' 
       iodmet(2)=2 
       iodmen(2)='VAISALA' 
    ELSEIF(imeth.eq.2)THEN 
       WRITE(iun_log,201) name(1:ln) 
201    FORMAT('Gauss method for object ',A,':') 
       iodnm=1 
       iodmet(1)=1 
       iodmen(1)='GAUSS' 
    ELSEIF(imeth.eq.3)THEN 
       WRITE(iun_log,202) name(1:ln) 
202    FORMAT('Vaisala method for object ',A,':') 
       iodnm=1 
       iodmet(1)=2 
       iodmen(1)='VAISALA' 
    ENDIF
    IF(PRESENT(ades_type))THEN
       CALL io_det_ades(iun_log,rwfil,name,obs,obsw,n,error_model,1,elem,   &
            &     telem,eltype,comele,fail,ades_type)                          
    ELSE
       CALL io_det(iun_log,rwfil,name,obs,obsw,n,error_model,1,elem,   &
            &     telem,eltype,comele,fail)
    ENDIF
    ! **********************************************************************
    deforb=(.NOT.fail) 
    defcov=.false. 
    IF(fail) THEN 
       WRITE(*,204) 
204    FORMAT(2X,'INITIAL ORBIT DETERMINATION FAILED') 
    ELSE 
       ! change to equinoctal                                                  
       CALL coocha(elem,eltype,gms,eq,'EQU',enne) 
       WRITE(*,220) telem,eq 
220    FORMAT(' preliminary orbit elements for epoch=',f12.4/6f13.7)
       el=undefined_orbit_elem
       el%coo='EQU'
       el%coord=eq
       el%t=telem 
       h=el%h_mag !default undefined values
       g=el%g_mag
       ! output new elements, multiline format, no header, in .fel file.       
       !CALL wro1lh(uniele,'ECLM','J2000','KEP') 
       CALL wromlr(uniele,name,elem,'KEP',telem,gg,.false.,cc,.false. &
            &        ,h,g,0.d0,6)  
       rms=rms_compute(obs,obsw,n)
       WRITE(iun_log,222)rms 
       WRITE(*,222)rms 
222    FORMAT(' RMS of residuals, preliminary orbit=',f10.4) 

    END IF
  END SUBROUTINE f_gauss


END MODULE fitobs_mod

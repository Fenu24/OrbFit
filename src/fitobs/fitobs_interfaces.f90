!                                                                       
! ====================================================                  
! FCLOMON2 impact monitoring from FITOBS
! ====================================================  
SUBROUTINE fclomon2(progna,astname,m,obs,obsw,mm1,mm2,tmcla,sigma,nd,limit_stretch,sigma0)
  USE dyn_param
  USE fund_const
  USE propag_state
  USE ephem_prop
  USE tp_trace
  USE cla_store
  USE output_control
  USE orbit_elements
  USE ret_analysistp
  USE multiple_sol
  USE multi_store
  USE obssto
  USE name_rules, ONLY: name_len
  USE force_model, ONLY: masjpl
  USE planet_masses, ONLY: dmea
  USE close_app, ONLY: fix_mole, kill_propag
  USE eval_risk, ONLY: massgiven,givenmass
  USE offlov_checktp, ONLY: bsdmin, bsdmax, header_risk,shrinkea, header_esa_risk
  IMPLICIT NONE
! =================INPUT=========================================
  character*(6), INTENT(IN) :: progna ! program name
  CHARACTER*(*), INTENT(IN) :: astname ! asteroid name 
  integer, INTENT(IN) :: m ! observation number
! new data types
  TYPE(ast_obs),DIMENSION(m), INTENT(IN) :: obs
  TYPE(ast_wbsr),DIMENSION(m), INTENT(IN) :: obsw
! min and max index of alternate orbits                                
  INTEGER, INTENT(IN) :: mm1
  INTEGER, INTENT(IN) :: mm2
! search until time (MJD)
  DOUBLE PRECISION, INTENT(IN) :: tmcla 
  DOUBLE PRECISION, INTENT(IN) :: sigma  ! max excursion in sigma
  INTEGER, INTENT(IN) :: nd ! dimension of the parameter space
  DOUBLE PRECISION, INTENT(IN) :: limit_stretch ! Limit in the stretching
  DOUBLE PRECISION, INTENT(IN), OPTIONAL :: sigma0 ! LOV parameter of center
! workspace por pro_ele 
  TYPE(orbit_elem) el1
  TYPE(orb_uncert) unm1
  CHARACTER*9 astna0
  INTEGER nm,j,le  ! number of virtual objects, index, string length
  INTEGER no                                    ! number of close 
                                                ! approaches found
!  INTEGER, PARAMETER :: nox=10000               ! max no close app
!  TYPE(tp_point), DIMENSION(nox) :: vas_trace      ! array of tp_traces(global)
!  TYPE(tp_point), DIMENSION(nox) :: vas_traceloc     ! array of tp_traces(of return)
!  DOUBLE PRECISION :: dist(nox)
!  DOUBLE PRECISION :: dminpos(nox) 
  LOGICAL :: is_min, newflag                    ! local minimum, low  MOID
  INTEGER :: no_risk                   ! number of Virtual impactors found
!  INTEGER, PARAMETER :: nshx=10000              ! maximum number of showers
!  INTEGER, PARAMETER :: nretx=50000             ! maximum number of returns 
  INTEGER            :: nsho                    ! showers number
 ! INTEGER            :: isho(nshx)              ! showers index
 ! INTEGER            :: iret(nretx)             ! returns (trails) index
  INTEGER            :: nret                    ! returns (trails) number
  INTEGER            :: lre,ire,nr,ns
  DOUBLE PRECISION   :: tcat                    ! time of initial conditions
  CHARACTER(LEN=15)  :: despla  ! desired planet 
  DOUBLE PRECISION :: stretchmin
  INTEGER            :: nmin            ! Minimum counter
!=======================IN/OUT FILES===================================
  CHARACTER(LEN=80) :: clodir,repdir,catname,obsdir,mulsodir,eledir,elefil
  CHARACTER(LEN=100):: riskfile,riskesafile,repfile,newfile,outfile,warfile, clo2file
  CHARACTER(LEN=100) :: wfile             ! for vvv_tp
! ============================UNITS FOR IN/OUT======================================
  INTEGER iunnew,iunwarn,iunout, iunrisk,iunrep, nrisk, nnew, iwdir
!=============================OPTIONS=====================================  
  DOUBLE PRECISION :: dt             ! length of showers 
  DOUBLE PRECISION :: tgap           ! time gap desired
  DOUBLE PRECISION :: dmin_used      ! radius of TP
  DOUBLE PRECISION :: dnewton        ! control on distance
  LOGICAL ireq,found
  INTEGER vdifold, vmultold
  CHARACTER*60 comment 
  DOUBLE PRECISION sigma00 ! for sigma0 optional
! Scatter plane
  REAL(KIND=qkind) :: vvv_q(ndimx)  ! weak direction (quadruple precision)
! =======================================================================
! modification August 16, 2005 to input mass from physical observations
  INCLUDE 'parlib.h90' 
  CHARACTER*200 filmass
  LOGICAL ok
  INTEGER iunmass
! memory management
  nox=30000
  nshx=1000
  nretx=5000
  ALLOCATE(vas_trace(nox),vas_trasrt(nox),isho(nshx),iret(nretx),tclo(nox),tclos(nox))
  noxloc=nox/10
  ALLOCATE(vas_traceloc(noxloc),vas_trlocsrt(noxloc),dist(noxloc),dminpos(noxloc),stretch(noxloc))
! ===============================================================
! input options
! ================shower length time============================  
  ireq=.false.
  dt=90.d0  
  comment='shower length time'
  CALL input_rea_opt(progna,'dt',dt,ireq,found,comment,iun_log)
! ================shower gap time============================           
  ireq=.false.
  tgap=30  
  comment='shower gap time'
  CALL input_rea_opt(progna,'tgap',tgap,ireq,found,comment,iun_log)
! ================radius of target plane used===========================
  dmin_used=dmea ! must be equal to the one used in the propagation
! ================control for startup of Newton's method================
! in earth radii                                                        
  ireq=.false.
  dnewton=60.d0 
  comment='newton control value, in Earth radii'
  CALL input_rea_opt(progna,'dnewton',dnewton,ireq,found,comment,iun_log)
! ================control factor on beta================
! to identify entangled                                
  ireq=.false. 
  beta_factor=20.d0
  comment='control factor on beta to identify entangled'
  CALL input_rea_opt(progna,'beta_factor',beta_factor,ireq,found,comment,iun_log)
! ================fix mole by using 2-body inside Earth?==============
  ireq=.false.
  fix_mole=.true.            ! default for fix_mole     
  comment='fix mole by using 2-body inside Earth'
  CALL input_log_opt(progna,'fix_mole',fix_mole,ireq,found,comment,iun_log)
! verbosity control: lower than for the rest of fitobs
  vdifold=verb_dif
  verb_dif=1
  vmultold=verb_mul
  verb_mul=1
! realignement of sigma0
  IF(PRESENT(sigma0))THEN
     sigma00=sigma0
  ELSE
     sigma00=0.d0
  ENDIF
! ============get mass known from other sources, if available==================
  filmass=dlibd//'/'//name_obj//'.mass'
  CALL rmsp(filmass,le)
  INQUIRE(file=filmass,exist=ok)
  IF(ok)THEN
     CALL filopn(iunmass,filmass(1:le),'old')
     READ(iunmass,*) givenmass
     massgiven=.true.
     CALL filclo(iunmass,' ')
  ELSE
     massgiven=.false.
  ENDIF
! END modification August 16, 2005 to input mass from physical observations
! ================OUTPUT ALSO ON FILES===============================
! open .out, .rep, .new, .warn, .risk
!----------------------------FILE .out-----------------------------------------
  iunout=0
! OPEN FILE .out
  outfile=astname//'.out' 
  CALL rmsp(outfile,le) 
  CALL filopn(iunout,outfile(1:le),'unknown') 
 !----------------------------FILE .rep-----------------------------------------
  repfile=astname//'.rep' 
  CALL rmsp(repfile,le) 
  CALL filopn(iunrep,repfile(1:le),'unknown') 
  CALL header_rep(iunrep)                           
!---------------------------FILE .new-----------------------------------------
  newfile=astname//'.new' 
  CALL rmsp(newfile,le) 
  CALL filopn(iunnew,newfile(1:le),'unknown') 
  CALL header_new(iunnew)                              
!---------------------------FILE .risk-----------------------------------------
  nrisk=0 
  bsdmin=1.d10 
  bsdmax=.5d0 
  shrinkea=1.d0
  riskfile=astname//'.risk' 
  CALL rmsp(riskfile,le) 
  CALL filopn(iunrisk,riskfile(1:le),'unknown') 
  WRITE(iunrisk,300) astname 
300 FORMAT('Object: ',a9,' '/) 
  CALL header_risk(iunrisk)
!--------------------------FILE .warn----------------------------------
  warfile=astname//'.warn' 
  CALL rmsp(warfile,le) 
  CALL filopn(iunwarn,warfile(1:le),'unknown') 
  CALL filclo(iunwarn,'DELETE') 
  CALL filopn(iunwarn,warfile(1:le),'unknown') 
  WRITE(*,*)' units iunout,iunnew,iunwarn,iunrisk ',iunout,iunnew,iunwarn,iunrisk
! ===================================================================
  IF(prob_sampl)THEN
! for densification when IP sampling this needs to be modified
     imul0=imi0
! sigma_max 
  ELSE
                       ! this is ok also for densification ????
! WARNING: begins modification to allow for densification
! modification 28/10/2013 to use in densification
! imul0, sigma_max must be recomputed in such a way that vas_tr(1:no)%sigma
! is aligned with sigma_q
! ???????????
!     imul0=imi0-sigma00/delta_sigma
     imi0=imi0-sigma00/delta_sigma
     sigma_max=abs(sigma00)+sigma ! range of LOV parameter is restricted, 
                          ! ignore close app with minima ouside range
  ENDIF

! AS AN ALTERNATIVE, CORRECT VALUE OF vas_tr(1:no)%sigma after input 
! deltasig, imul0, sigma_max must be recomputed in such a way that vas_tr(1:no)%sigma
! is aligned with sigma_q
! if imi0 corresponds to sigma0, sigma(imi)=sigma0+(imi-imi0)*delta_sigma
! on the other hand rescaltp uses sigma(imi)=deltasig*(imi-imul0)
! END WARNING 

  CALL masjpl  ! should not matter, unless fclomon2 is executed too early....
! ======================================================================== 
! copy observations 
  m_m=m
  obs_m(1:m)=obs(1:m)
  obsw_m(1:m)=obsw(1:m)
! setup close app. output
  tpplane=.true. ! use TP plane, not MTP
  REWIND(iuncla) ! cleanup previous close app. records
  IF(scatterplane) THEN
     WRITE(*,*) 'FITOBS: scatter_stretch = TRUE' 
     scatter_stretch=.TRUE.
     wfile=astname//'.wdir'
     CALL rmsp(wfile,le)
     CALL filopn(iwdir,wfile(1:le),'old')
! read quadruple precision, convert in double 
     READ(iwdir,*) vvv_q(1:nd)
     vvv_tp(1:nd)=vvv_q(1:nd)
     CALL filclo(iwdir,' ')   
  END IF
! propagation and output close approach files
  nm=mm2-mm1+1
  DO j=mm1,mm2
     astna0=' '
     IF(j.le.9)THEN              
        WRITE(astna0,111) j 
111     FORMAT('va_00000',I1)
     ELSEIF(j.le.99)THEN
        WRITE(astna0,112) j
112     FORMAT('va_0000',I2)
     ELSEIF(j.le.999)THEN
        WRITE(astna0,113) j
113     FORMAT('va_000',I3)
     ELSEIF(j.le.9999)THEN
        WRITE(astna0,114) j
114     FORMAT('va_00',I4)
     ELSEIF(j.le.99999)THEN
        WRITE(astna0,115) j
115     FORMAT('va_0',I5)
     ELSEIF(j.le.999999)THEN
        WRITE(astna0,116) j
116     FORMAT('va_',I6)
     ELSE
        WRITE(*,*)' nmulti: change format for imult>499999'
        STOP
     ENDIF
     WRITE(iuncla,*)astna0
     CALL cov_avai(unm(j),elm(j)%coo,elm(j)%coord)
     IF(dyn%nmod.GT.0)THEN
        dyn%dp(1:dyn%ndp)=dpm(1:dyn%ndp,j)
     END IF
     WRITE(*,*) ' VA number ',j
     CALL pro_ele(elm(j),tmcla,el1,unm(j),unm1)
  ENDDO
! Rewind
  REWIND(iuncla)
  kill_propag=.FALSE.
! now input tp records
  despla='EARTH'
  CALL inclolinctp(iun_log,iuncla,despla,no)

  IF(no.le.0)THEN
     WRITE(*,*) ' no close approach to ', despla, ' found up to time MJD ',tmcla
     RETURN
  ENDIF
! close .clo file and open .clo2 file
  CALL filclo(iuncla,' ')
! close approach file (but for the falsi runs, to be deleted)           
  clo2file=name_obj//'.clo2' 
  CALL rmsp(clo2file,le) 
  CALL filopn(iuncla,clo2file(1:le),'unknown') 
! shower analysis
  CALL showret3tp(iunout,no,dt,tgap,nsho,nret)
! main loop on returns                                                  
  ns=1 
  CALL header_rep(iun_log)
  no_risk=0
  nnew=0
  DO 1 nr=1,nret 
! shower counter                                                        
     if(iret(nr).ge.isho(ns+1))ns=ns+1 
! reopen files for filament report                                      
     ire=iret(nr) 
! warning: remember to define iret(nr+1)                                
     lre=iret(nr+1)-iret(nr) 
! analyse return, finding minimum distance etc.                         
     WRITE(iunout,197)ns,nr,lre,vas_trace(ire)%rindex,vas_trace(ire+lre-1)%rindex 
197  FORMAT('shower no. ', i4,' return no. ',i5,' length ',i5,' from ',f6.1,' to ',f6.1)
! vas_trace  copied in a "return record' vas_traceloc                                  
     CALL arrcut(vas_trace,ire,lre,nox,iunout,vas_traceloc) 
     dist(1:lre)=vas_traceloc(1:lre)%b
     dminpos(1:lre)=vas_traceloc(1:lre)%minposs
     newflag=.false. 
     nmin=0
     IF(lre.gt.1)THEN 
! finding local minima                                                  
        DO 2 j=1,lre 
           is_min=.false. 
           IF(j.eq.1)THEN 
              IF(dist(j).lt.dist(j+1))is_min=.true. 
           ELSEIF(j.eq.lre)THEN 
              IF(dist(j).lt.dist(j-1))is_min=.true. 
           ELSE 
              IF(dist(j).lt.dist(j+1).and.dist(j).lt.dist(j-1))     &
     &                 is_min=.true.                                    
           ENDIF
! if local minimum, check minimum possible                              
           IF(is_min)THEN
              CALL header_rep(iun_log)
              CALL header_rep(0)
              CALL wrireptp(vas_traceloc(j),vas_traceloc(1)%rindex,vas_traceloc(lre)%rindex,iunrep)   
              IF(dminpos(j).lt.dnewton)THEN 
                 newflag=.true. 
                 nmin=nmin+1
                 stretch(nmin)= vas_traceloc(j)%stretch
              ENDIF
           ENDIF
2       ENDDO
     ELSE 
! singletons: check if moid is small 
        CALL header_rep(iunrep) 
        CALL wrireptp(vas_traceloc(1),vas_traceloc(1)%rindex,vas_traceloc(lre)%rindex,iunrep)
        DO j=1,lre 
           IF(dminpos(j).lt.dnewton)THEN
              newflag=.true.
              nmin=nmin+1
              stretch(nmin)= vas_traceloc(1)%stretch
           END IF
        ENDDO
     ENDIF
! selected for further analysis?                                        
     IF(newflag)THEN 
        IF(limit_stretch.GT.0.d0)THEN
           stretchmin=MINVAL(stretch(1:nmin))
           IF(stretchmin.GT.limit_stretch)THEN
              CLOSE(iunrep) 
              CLOSE(iunout)
              CYCLE   
           ENDIF
        ENDIF
! reopen files for newton report                                        
        WRITE(*,*)' return number nr=',nr, ' falsi/newton ' 
        nnew=nnew+1 
        CALL ret_min(lre,vas_traceloc,tdt_cat,dnewton,no_risk,iunnew,iunwarn,iunrisk,limit_stretch)
! test output                                                           
        nrisk=nrisk+no_risk 
     ELSE 
        WRITE(*,*)' return number nr=',nr, 'NO  falsi/newton ' 
     ENDIF
! end loop on returns 
1 ENDDO
  verb_dif=vdifold
  verb_mul=vmultold
! close output files
  IF(nrisk.gt.0)THEN
     CALL filclo(iunrisk,' ')
  ELSE
     CALL filclo(iunrisk,'DELETE')
  ENDIF
  CALL filclo(iunout,' ')
  CALL filclo(iunrep,' ')
  IF(nnew.gt.0)THEN
     CALL filclo(iunnew,' ')
     CALL filclo(iunwarn,' ')
  ELSE
     CALL filclo(iunnew,'DELETE')
     CALL filclo(iunwarn,'DELETE')
  ENDIF
END SUBROUTINE fclomon2

! ============ FINOBS =========

SUBROUTINE f_gaussdeg8(uniele,name,deforb,defcov,          &
     &     rwfil,obs,obsw,n,error_model,el)
  USE fund_const, ONLY: gms 
  USE astrometric_observations 
  USE least_squares, ONLY: rms_compute
  USE orbit_elements
  USE output_control
  IMPLICIT NONE 
! INPUT observations: new data types
  INTEGER, INTENT(IN) :: n  !number of observations
  TYPE(ast_obs),DIMENSION(n),INTENT(IN) :: obs
  TYPE(ast_wbsr),DIMENSION(n),INTENT(INOUT) :: obsw
  CHARACTER*20, INTENT(IN) ::  error_model ! weighing model
! other input   
  INTEGER,INTENT(IN) :: uniele ! unit for output
  CHARACTER*18, INTENT(IN) :: name !asteroid name 
  CHARACTER*60, INTENT(IN) :: rwfil ! output file for residuals
! OUTPUT
  LOGICAL, INTENT(OUT) :: deforb,defcov ! orbit state on output 
  TYPE(orbit_elem), INTENT(OUT) :: el ! elements  
! END INTERFACE
  DOUBLE PRECISION, DIMENSION(3) :: tobs,alpha3,delta3
  DOUBLE PRECISION tr,dt,dt1 ! mean time, half interval, dist. to center
  INTEGER :: obscod(3), isel(3),j,nselect
! for call to gaussdeg8
  INTEGER nroots,nsol
  LOGICAL fail, debug
  CHARACTER*(20) msg 
  DOUBLE PRECISION, PARAMETER :: ecc_max=2.d0, q_max=500.d0
  TYPE(orbit_elem), DIMENSION(3) :: elv
  DOUBLE PRECISION :: rr(3) ! topocentric distance
  DOUBLE PRECISION rms
! select obs: first, last closest to mean time
     defcov=.false.
  IF(n.le.2)THEN
     deforb=.false.
     RETURN
  ENDIF
  isel(1)=1
  isel(3)=n
  tr= (obs(1)%time_tdt+obs(n)%time_tdt)/2.d0
  dt=(obs(n)%time_tdt-obs(1)%time_tdt)/2.d0
  isel(2)=1
  DO j=2,n
    dt1=abs(obs(j)%time_tdt-tr)
    IF(dt1.lt.dt)THEN
       isel(2)=j
       dt=dt1
    ENDIF
  ENDDO
  IF(isel(2).eq.1.or.isel(2).eq.n)THEN
     WRITE(*,*) 'f_gaussdeg8: logical error in times'
     WRITE(*,*) obs(1:n)%time_tdt
     STOP
  ENDIF     
! copy in array
  DO j=1,3
    tobs(j)=obs(isel(j))%time_tdt
    alpha3(j)=obs(isel(j))%coord(1)
    delta3(j)=obs(isel(j))%coord(2)
    obscod(j)=obs(isel(j))%obscod_i
    obsw(isel(j))%sel_coord=2
  ENDDO
! find solution(s)
  debug=.true.
  CALL gaussdeg8(tobs,alpha3,delta3,obscod,ecc_max,q_max,elv,nroots,nsol,rr,fail,msg,debug)
  IF(msg.ne.' ')WRITE(*,*)' error message from gaussdeg8=',msg
! assess solutions
  IF(nsol.eq.0)THEN
     deforb=.false.
     RETURN 
  ENDIF
!  DO j=1,nsol
!     rms=rms_compute(obs,obsw,n)
!     WRITE(iun_log,222)j, rms 
!     WRITE(*,222)j, rms 
!222  FORMAT(' preliminary orbit no. ',i3,' RMS of residuals=',f10.4)  
!  ENDDO
! select solution, if more than one
  IF(nsol.eq.1)THEN
     nselect=1
  ELSE
22   WRITE(*,*)' select one solution between 1 and ', nsol
     READ(*,*) nselect
     IF(nselect.lt.1.or.nselect.gt.nsol) GOTO 22
  ENDIF
  deforb=.true.
  el=elv(nselect)  
  WRITE(*,*) ' selected preliminary orbit ', el 

END SUBROUTINE f_gaussdeg8

! ====================== F_GAUSS ==========

! ====================================================                  
! FOBPRE predict observations                                           
! ====================================================                  
SUBROUTINE fobpre(icov,ini00,cov00,ok,titnam,filnam,  &
     &   el00,unc00,ids,type1,t1,tut,aobs0,dobs0,t2,dt,astnam)
  USE multiple_sol, ONLY: outmul 
  USE pred_obs
  USE orbit_elements
  USE util_suit
  USE output_control
  USE astrometric_observations
  USE two_states
  USE fund_const
  USE station_coordinates, ONLY: codestat,obscoo
  USE reference_systems, ONLY: observer_position
  USE ephem_prop
  IMPLICIT NONE 
! =================INPUT=========================================       
  INTEGER,INTENT(IN) :: icov ! requirements on covariance
  LOGICAL ini00,cov00 ! availability of initial conditions, covariance
  DOUBLE PRECISION, INTENT(IN) :: t1,t2,dt ! observation time, 
! also beginning of ephemerides time, end of ephemerides time, step
  DOUBLE PRECISION tut ! UTC of observation 
  INTEGER ids ! station code
  CHARACTER*(1) type1 ! observation type  
  TYPE(orbit_elem), INTENT(IN) :: el00 ! elements
  TYPE(orb_uncert), INTENT(IN) :: unc00 !covariance and normal matrix 
  CHARACTER*80, INTENT(INOUT) :: titnam ! asteroid name etc. 
  CHARACTER*60, INTENT(IN) :: filnam                     
  DOUBLE PRECISION, INTENT(IN) :: aobs0,dobs0 ! actual obs. for comparison 
  CHARACTER*(*), INTENT(IN) :: astnam ! asteroid name
! ONLY DIRECT OUTPUT
  LOGICAL, INTENT(OUT) :: ok ! necessary data available
! ===== predicted observations ===========                              
! angles (best fit prediction), apparent magnitude
  DOUBLE PRECISION alpha,delta,hmagn
! covariance of the observations                                        
  DOUBLE PRECISION gamad(2,2),axes(2,2),sig(2),gamad1(2,2)
! noise in observations
  DOUBLE PRECISION rmssec, rmsmag
! confidence boundary, line of max variation 
  INTEGER npo, ibv, npo1, npop 
  DOUBLE PRECISION sigma
  DOUBLE PRECISION :: aobs,dobs,adot,ddot ! alpha, delta, proper motion
  DOUBLE PRECISION :: pha,dis,dsun,elo,gallat,gallon,elmoon,elev,elsun ! phase, dist. Earth, dist. Sun
  INTEGER  inl ! menu: handling of nonlinearity
  CHARACTER*20 menunam ! menu                                     
  CHARACTER*100 file,fields ! ephemerides output 
  CHARACTER*3 scale 
  INTEGER ln,iuneph,le 
! ===================================================================== 
! options                                                               
! ===================================================================   
! ===================================================================   
! chose handling of nonlinearity                                        
7 IF(icov.ge.3.and.icov.lt.6)THEN 
     menunam='prednonl' 
     CALL menu(inl,menunam,3,'How to handle nonlinearity?=',        &
     &         'linear map=',                                           &
     &         '2-body nonlinearity=',                                  &
     &         'full n-body nonlinearity=')
     IF(inl.eq.0)GOTO 7 
  ELSEIF(icov.eq.6)THEN
     inl=1
  ELSE ! icov=1 for simple obs, icov=2 for use simulated obs
     inl=-1
  ENDIF
  ibv=0
! ===================================================================== 
! check availability of initial conditions (also covariance for icov=2) 
  CALL chereq(icov,ini00,cov00,el00%t,iun_log,ok) 
  IF(.not.ok)RETURN 
! ===================================================================== 
! check availability of JPL ephemerides and ET-UT table                 
  CALL chetim(t1,t1,ok) 
  IF(.not.ok)RETURN 
! ===================================================================== 
! compute prediction; without and with covariance                       
! ===================================================================== 
  IF(icov.eq.1.or.icov.eq.2)THEN 
! ===================================================================== 
! only alpha, delta, magnitude
     inl=1                                          
     CALL predic_obs(el00,ids,t1,type1,             &
     &        alpha,delta,hmagn,inl,                                      &
     &        ADOT0=adot,DDOT0=ddot,PHA0=pha,DIS0=dis,                  &
     &        DSUN0=dsun,ELO0=elo,GALLAT0=gallat,ELMOON0=elmoon,GALLON0=gallon,ELSUN0=elsun,ELEV0=elev)
    IF(icov.eq.1)THEN
     CALL outobc(iun_log,type1,ids,tut,alpha,delta,hmagn,adot,ddot,    &
     &     elo,dis,icov,gamad,sig,axes,elmoon,gallat,gallon,elev,dsun,elsun,pha)

     ELSE
! add second arc, if not there
        IF(.not.obsp)THEN
           obsp=.true.
           astnap='sim'
           rwofip='sim.rwo'
           elefip='sim.fel'
           CALL rmsp(elefip,le) 
           CALL filopn(iunelp,elefip(1:le),'unknown')
        ENDIF
        IF(obs0)obstwo=.true.
! weights and elements files for identification                         
        IF(obstwo.and.iunelt.eq.0)THEN 
           CALL titast(3,astna0,astnap,titnam,rwofil,le) 
           CALL rmsp(rwofil,le) 
           rwotwo=rwofil(1:le)//'.rwo' 
! compose elements file full name                                       
           eletwo=rwofil(1:le)//'.fel' 
           CALL rmsp(eletwo,le) 
           CALL filopn(iunelt,eletwo(1:le),'unknown') 
! output header
           CALL wromlh (iunelt,'ECLM','J2000') 
        ENDIF
! select noise level
         WRITE(*,*)' give RMS of simulated observation'
         IF(type1.eq.'O')THEN
            WRITE(*,*)' for both alpha, delta in arcsec'
            READ(*,*) rmssec
         ELSEIF(type1.eq.'R')THEN
            WRITE(*,*)' for distance, in km'
            READ(*,*) rmssec
         ELSEIF(type1.eq.'V')THEN
            WRITE(*,*)' for range rate, in km/day'
            READ(*,*) rmssec
         ENDIF
         rmsmag=0.7d0
! add observation to second arc
         mall=mall+1
         mp=mp+1   
         CALL obs_simul(type1,t1,tut,astnam,ids,rmssec,rmsmag,alpha,delta,   &
      &                hmagn,obs(mall),obsw(mall))
         IF(type1.eq.'R'.or.type1.eq.'V')CALL radar_ob(obs(1:mall)%type,mall)
      ENDIF
   ELSEIF(icov.eq.3)THEN 
! ===================================================================== 
! alpha, delta, magnitude, covariance and ellipse of confidence         
     CALL predic_obs(el00,ids,t1,type1,             &
     &        alpha,delta,hmagn,inl,                                    &
     &        UNCERT=unc00,GAMAD=gamad,SIG=sig,AXES=axes,          &
     &        ADOT0=adot,DDOT0=ddot,PHA0=pha,DIS0=dis,                  &
     &        DSUN0=dsun,ELO0=elo,GALLAT0=gallat,GALLON0=gallon,ELMOON0=elmoon,ELSUN0=elsun,ELEV0=elev)
     CALL outobc(iun_log,type1,ids,tut,alpha,delta,hmagn,adot,ddot,    &
     &     elo,dis,icov,gamad,sig,axes,elmoon,gallat,gallon,elev,dsun,elsun,pha)                                 
! ===================================================================   
! generation of sky epehemrides                                         
  ELSEIF(icov.eq.6)THEN 
! check availability of JPL ephemerides and ET-UT table for entire time 
     CALL chetim(t1,t2,ok) 
     IF(.not.ok)RETURN 
     IF(nint(abs(t2-t1)/dt).gt.6000)THEN 
        write(*,*)' Too many ephemerides points:',                  &
     &           nint(abs(t2-t1)/dt)                                    
        write(*,*)'Select a time interval and span to ',            &
     &           'ensure that there are fewer than 6000 points.'         
     ELSE 
! open ephemerides file in current directory                            
        file=astnam//'.eph' 
        CALL rmsp(file,ln) 
        CALL filopn(iuneph,file(1:ln),'unknown') 
        fields='cal,mjd,coord,mag,elev,airm,elsun,elong,mooel,glat,glon,r,delta,appmot,skyerr' 
        scale='UTC' 
        CALL ephemc(iuneph,el00,unc00,.true.,t1,t2,dt,ids,scale,fields)
        CALL filclo(iuneph,' ') 
        WRITE(*,*)' Generated ephemeris in file: ',file(1:ln) 
     ENDIF
! ===================================================================== 
  ELSEIF(icov.eq.4.or.icov.eq.5)THEN 
! ===================================================================== 
! alpha, delta, magnitude, covariance and confidence boundary;          
! input specification of set of points                                  
     CALL asscbd(iun_log,npoinx,npo,sigma,ibv) 
! ===================================================================== 
! compute prediction, boundary                                          
     CALL predic_obs(el00,ids,t1,type1,             &
     &        alpha,delta,hmagn,inl,                                    &
     &        unc00,sigma,npo,ibv,gamad,sig,axes,npo1,                   &
     &        adot,ddot,pha,dis,dsun,elo,gallat,ELMOON0=elmoon,GALLON0=gallon,ELSUN0=elsun,ELEV0=elev)
     CALL outobc(iun_log,type1,ids,tut,alpha,delta,hmagn,adot,ddot,    &
     &        elo,dis,icov,gamad,sig,axes,elmoon,gallat,gallon,elev,dsun,elsun,pha)                             
     IF(npo1.le.0)THEN 
        WRITE(*,*)'fobpre: no elliptic orbits ',npo1 
        RETURN 
     ENDIF
     IF(ibv.eq.1)THEN 
! confidence boundary; one point added to close line                    
        al_m(npo1+1)=al_m(1) 
        de_m(npo1+1)=de_m(1) 
        hmag_m(npo1+1)=hmag_m(1) 
        el_m(1:6,npo1+1)=el_m(1:6,1) 
!            IF(inl.eq.3)disv(npo+1)=disv(1) 
        npop=npo1+1 
     ELSEIF(ibv.eq.2)THEN 
! line of variations                                                    
        npop=npo1 
     ENDIF
! if no observation is given, use the nominal marked with a cross       
     IF(icov.eq.4)THEN 
        aobs=alpha 
        dobs=delta 
     ELSE
        aobs=aobs0
        dobs=dobs0
     ENDIF
! ===================================================================== 
! output observation, apparent motion, confidence boundary              
     CALL outmul(titnam,filnam,tut,sigma,alpha,delta,               &
     &              al_m,de_m,hmag_m,1,npop,1,icov-3,aobs,dobs,type1)
  ENDIF
END SUBROUTINE fobpre                                      
! ===================================================================   
! FINELE                                                                
! ===================================================================   
! input of initial conditions                                           
SUBROUTINE finele(progna,iar,astna,el,ini,cov,unc)
  USE fund_const                                
  USE orbit_elements
  USE output_control 
  USE dyn_param          
  IMPLICIT NONE 
! ==========INPUT===================                                    
! file names, i/o control                                      
  CHARACTER*(*), INTENT(IN) :: astna  ! asteroid full name                               
  CHARACTER*6, INTENT(IN)  :: progna  ! program name 
  INTEGER, INTENT(IN)      :: iar     ! flag for arcs 1-2 
! ==========OUTPUT==================                                    
! epoch time (MJD), elements (equinoctal), absolute magnitude, opp.effec
  TYPE(orbit_elem), INTENT(INOUT) :: el !remains defined!
  LOGICAL, INTENT(OUT) ::  ini ! successful input flag 
  TYPE(orb_uncert), INTENT(OUT) :: unc ! covariance matrices 
  LOGICAL, INTENT(OUT) :: cov ! covariance available
! =========END INTERFACE=============                                   
! asteroid name for elements (18 characters)                            
  CHARACTER*50 namel,name
  INTEGER lnam, lnam1                                                   
  LOGICAL ireq, found ! logical input flags 
  CHARACTER*60 comment                     
  INTEGER le,lench ! length of names
! variables for reading routines  
  CHARACTER*100 elefi   ! file name for elements    
! output from read_elems
  INTEGER ndp, nlsloc, nmod, lsloc(ndyx) ! LSP record
  DOUBLE PRECISION dp(ndyx) ! nongrav parameters
  TYPE(dyn_par) :: dynam_par
  CHARACTER*4 form ! 1l, ML
  LOGICAL err,eof ! error, end of file (without error)
! ====================================  
  lnam=lench(astna) 
! this arc does not exist                                               
  IF(lnam.eq.0)RETURN 
! ====================================================                  
! inizialization of ini
  ini=.FALSE.
! check that the asteroid name has been read 
  ireq=.false.
  namel=' ' 
  namel=astna
  IF(iar.eq.1)THEN 
     comment='first arc asteroid name'
     CALL input_cha_opt(progna,'namel0',namel,ireq,found,comment,iun_log)
  ELSEIF(iar.eq.2)THEN 
     comment='second arc asteroid name'
     CALL input_cha_opt(progna,'namelp',namel,ireq,found,comment,iun_log)
  ENDIF
  CALL rmsp(namel,lnam) 
! read the name of the elements file, inquire                           
  ireq=.false. 
  elefi='ast.cat' ! default not used; only to fail graciously when no orbit given
  IF(iar.eq.1)THEN 
     comment='first arc elements file'
     CALL input_cha_opt(progna,'elefi0',elefi,ireq,found,comment,iun_log)
  ELSEIF(iar.eq.2)THEN 
     comment='second arc elements file'
     CALL input_cha_opt(progna,'elefip',elefi,ireq,found,comment,iun_log)
  ENDIF
  CALL rmsp(elefi,le)
!============================================= 
  INQUIRE(file=elefi,exist=found) 
  IF(found)THEN 
     CALL read_elems(el,name,eof,nlsloc,err,lsloc,form,unc,elefi)
     CALL rmsp(name,lnam1)
     IF(name.ne.namel)THEN
! WARNING: elefi must contain only the orbit of the one we are looking for; scanning
!          of the file for finding one object has been disabled.
        WRITE(*,*)' finele: asteroid ',namel(1:lnam),' not found in file ',elefi(1:le),' found ',name(1:lnam1)
        RETURN
     ENDIF
! control of errors
     IF(err) RETURN
! initial conditions found
     ini=.true.     
     IF(unc%succ)THEN
        cov=.true.
     ELSE
        cov=.false. ! even if covariance was available before, it is not applicable now
     ENDIF
     WRITE(*,*) 'input orbit file correctly read'
     WRITE(iun_log,*) 'input orbit file correctly read'
     WRITE(*,*) el
     WRITE(iun_log,*) el
! if consistent with the model being used, store dp in dyn_param.dyn
     IF(.NOT.ngr_opt.AND.dyn%ndp.GT.0.AND.nlsloc.NE.0)THEN 
        nls=nlsloc
        ls=lsloc
     ENDIF
! lsloc, nlsloc are not used, that is the options in the .fop file prevail on the options used
!        in the computation of orbit contained in elefi
  ELSE 
     WRITE(*,*) 'File ',elefi(1:le),' not found!' 
     WRITE(iun_log,*) 'File ',elefi(1:le),' not found!' 
  ENDIF
END SUBROUTINE finele
                                         
! ===================================================================== 
! WRIEQU (write initial conditions, equinoctal)                         
! ===================================================================== 
subroutine wriequ(iun,astna0,t0,eq0) 
  implicit none 
  double precision t0,eq0(6) 
  integer iun 
  character*18 astna0 
! initial conditions found                                              
  write(*,108)astna0,t0 
108 format(1x,a18,' initial elem (a,h,k,p,q,lam), epoch=',f8.1) 
  write(iun,104) eq0 
  write(*,104) eq0 
104 format(6f13.7) 
  write(iun,*)' ' 
  return 
END SUBROUTINE wriequ

! ========================================                              
! FIDENT compute identification norm                                    
! =======================================                               
SUBROUTINE fident(id_dim,cov0,covp,el0,elp,unc0,uncp,elid,ff)  
  USE orbit_elements 
  USE output_control
  IMPLICIT NONE 
  LOGICAL, INTENT(OUT):: ff !failure of identification proposal
  LOGICAL fail 
  INTEGER,INTENT(IN) :: id_dim ! method control (normally 6) 
  TYPE(orbit_elem), INTENT(IN) :: el0,elp
  LOGICAL,INTENT(IN) :: cov0,covp
  TYPE(orb_uncert), INTENT(IN) :: unc0,uncp
  TYPE(orbit_elem), INTENT(INOUT) :: elid
! two input elements, their uncertainties
  DOUBLE PRECISION eq0(6),eqp(6)
  DOUBLE PRECISION c0(6,6),cp(6,6),g0(6,6),gp(6,6) 
! identification elements                                          
  DOUBLE PRECISION eqf(6) 
! similarity norms   
  DOUBLE PRECISION d2,da2,dista,dq,dqalt 
! determinants and eigenvalues                                          
  DOUBLE PRECISION detc2,detc5,detc6,eigen5(5),eigen6(6) 
! ========================================                              
! check requirements 
  ff=.false.                                                   
  IF(.not.cov0)THEN 
     WRITE(*,*)' covariance matrix for arc 1 not ready' 
     RETURN 
  ENDIF
  IF(.not.covp)THEN 
     WRITE(*,*)' covariance matrix for arc 2 not ready' 
     RETURN 
  ENDIF
  IF(el0%t.ne.elp%t)THEN 
     WRITE(*,*)' to=',el0%t,' tp=', elp%t
     WRITE(*,*)' elements and covariances must be available' 
     WRITE(*,*)' for the same time, use propagation first' 
     RETURN 
  ENDIF
!  IF(el0%coo.ne.'EQU'.or.elp%coo.ne.'EQU')THEN
!     WRITE(*,*)' identification guess only in EQU, not in ', el0%coo, elp%coo
!     RETURN
!  ENDIF
  eq0=el0%coord
  eqp=elp%coord
  g0(1:6,1:6)=unc0%g(1:6,1:6)
  c0(1:6,1:6)=unc0%c(1:6,1:6)
  gp(1:6,1:6)=uncp%g(1:6,1:6)
  cp(1:6,1:6)=uncp%c(1:6,1:6)
! selection of algorithm                                                
  IF(id_dim.eq.5)THEN 
! identification by 5x5 matrix                                          
     CALL idno5(eq0,eqp,g0,c0,gp,cp,                                &
     &    d2,da2,dq,dqalt,dista,detc2,detc5,eqf,fail,eigen5)            
  ELSEIF(id_dim.eq.6)THEN 
! identification by 6x6 matrix                                          
     CALL idno6(eq0,eqp,g0,c0,gp,cp,                                &
     &    d2,da2,dq,dqalt,dista,detc2,detc6,eqf,fail,eigen6)            
  ELSE 
     WRITE(*,*)' fident: id_dim=',id_dim,' not understood' 
     RETURN 
  ENDIF
! output                                                                
  IF(fail)THEN 
     CALL tee(iun_log,' FAILED IDENTIFICATION ALGORITHM=') 
     ff=.false. 
  ELSE 
     WRITE(*,*)' INITIAL COND. AVAILABLE, NOW USE DIFF. CORR.' 
     ff=.true. 
! store as proposed identification 
     elid=el0 ! magnitude from the first...
     elid%coord=eqf
  ENDIF
! penalties                                                             
  WRITE(*,*)' d2, d2alt,  dq,   dqalt,    dista' 
  WRITE(*,194)d2,da2,dq,dqalt,dista 
194 FORMAT(4(f13.4,1x),f10.6) 
  WRITE(iun_log,*)' d2, d2alt,  dq,   dqalt,    dista' 
  WRITE(iun_log,194)d2,da2,dq,dqalt,dista 
! proposed elements                                                     
  WRITE(iun_log,208)elid%t 
208 FORMAT(' ident. elem (a,h,k,p,q,lam), epoch=',f8.1) 
  WRITE(*,208)elid%t 
  WRITE(*,104)eqf 
  WRITE(iun_log,104)eqf 
104 FORMAT(6f13.7) 
! store as proposed identification                                      
END SUBROUTINE fident








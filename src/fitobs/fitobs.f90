! Copyright (C) 1997,2003 ORBFIT consortium
! version 3.1: A. Milani December 2003
! ===================================================================
! PROGRAM FITOBS
! ==================================================================
! Two arc differential corrector, state propagator,
!     observation predictor
! ===================================================================
PROGRAM fitobs
  USE fund_const
  USE astrometric_observations
  USE orbit_elements
  USE output_control
  USE least_squares
  USE multiple_sol
  USE multi_store
  USE yark_pert
  USE util_suit
  USE force_model
  USE two_states
  USE attributable
  USE arc_control
  USE virtual_impactor
  USE tp_trace, ONLY:  wri_tppoint,tpplane,dtpde,mtp_store,njc, vvv_tp, scatter_stretch
  USE close_app, ONLY: fix_mole, kill_propag
  USE dyn_param
  USE cobweb
  USE station_coordinates, ONLY : statcode
  USE triangles
  USE ephem_prop
  USE fitobs_mod
  IMPLICIT NONE
! ===== input observations ======================
  INTEGER ifobs ! menu index
! ===== input orbital elements===================
  INTEGER imeth,iele,nselect ! to select prelim orbit metod, which arc
  TYPE(orbit_elem) el3(3) ! alternate preliminary orbits
  LOGICAL failpre
  CHARACTER*(50) msg
  INTEGER nsolp,nroots ! no preliminary orbits, no roots deg. 8
  DOUBLE PRECISION rr3(3) ! alternate topocentric distances
! ===== differential corrections ==============
  LOGICAL succ ! success flag
  double precision meanti ! function for time of difcor
  INTEGER idif, itsav, ist_fit, incfit ! menu index, save itmax when changed, step fit flag, No. param
  DOUBLE PRECISION peq(6) ! unit vector in weak direction
  DOUBLE PRECISION   :: h_ini, g_ini ! H and G magnitudes for correction of mag_est
! ======== proposed identifications =======
! equinoctal orbital elements (guess from identif), epoch time
  TYPE(orbit_elem) :: elide
  DOUBLE PRECISION enm ! mean motion for joint orbits
  INTEGER ng, igue !  number of revolutions, menu index
  LOGICAL ff  ! identification routine success flag
  LOGICAL inide, iniboth !  successful input flags, for identifications
! ============ propagation ============================
  DOUBLE PRECISION tr ! target time
  DOUBLE PRECISION tf,step,interv ! generate ephemerides
  INTEGER numsav, iprop,iyes ! ephem records, menu index, yes/no interrogation
  CHARACTER*(8) ans
! ===== predicted observations ===========
  DOUBLE PRECISION  tut1,tut2,t1,t2,dt ! time (TDT, UTC)
  CHARACTER*(1) type1  ! obs. type
  INTEGER ids, iprob ! station code, menu index
  INTEGER  im  ! confidence boundary, line of max variation
! ======== multiple solutions ==========================
  DOUBLE PRECISION sigma, sigma0, csino00  ! max value, central value,
!  residual norm at nominal value
  logical mult_adopted
  INTEGER imult,ifff,iff,iffat,imi,nmult,iarm,marc,iscal,isampl,iscat,iwdir
  INTEGER m1,m2 ! interval
  CHARACTER*160 catname ! input multiple solution catalog
! propagation times, close appr. analysis time
  DOUBLE PRECISION trmult,tcmult,tmcla, dtclo
! minimum v_infty with respect to Earth (circular approx.)
  DOUBLE PRECISION vel_inf
! max value of sigma change allowed for newton step
  DOUBLE PRECISION siglim
  INTEGER nvi, iunvi
  TYPE(orbit_elem) elop
  INTEGER, PARAMETER :: nj=20
  DOUBLE PRECISION :: del(2,nj)
! scatter plane
  REAL(KIND=qkind) :: vvv_q(ndimx)  ! weak direction (quadruple precision)
  CHARACTER(LEN=100) :: wfile             ! for vvv_tp
  DOUBLE PRECISION :: limit_stretch ! Limit in the stretching
! ===========close approach analysis======================
  INTEGER iclan
  DOUBLE PRECISION tlim
  CHARACTER*120 helpfi,ddocd1
  INTEGER ll,iunit,lench
  INCLUDE 'doclib.h90'
! ===============to change coord=================
  DOUBLE PRECISION, DIMENSION(6,6) :: dee !partials
  CHARACTER*3 cooy  ! target type
  CHARACTER choice  ! choice if you want to propagate even if nd.NE.unc0%ndim
  INTEGER fail_flag, icoord ! control coo_cha, menu index
  !Wrikep
  LOGICAL :: pha
  DOUBLE PRECISION :: moid
  CHARACTER*60 :: elefil
! ===============attributables ===================
  TYPE(attrib) attr0,attrp,attr,attrc !attributables
  DOUBLE PRECISION trou ! rounded time (at the integer MJD)
  DOUBLE PRECISION, PARAMETER :: sphx=2.d0 ! max arc span in degrees
  INTEGER iatt ! attributable menu selection
  LOGICAL error ! only 1 obs and so on
  DOUBLE PRECISION r, rdot ! to complete an attributable to an ATT elem.
  INTEGER artyp, nigarc, fail_arty ! arc type
  INTEGER obscodi ! for ATT elements
  INTEGER :: iunpol, iunfla ! file unit for admissibile region
  DOUBLE PRECISION rmsmax,geoc_chi, acce_chi, chi !arc type
  DOUBLE PRECISION dx2n,dx2v ! for Poincare' topocentric correction
  INTEGER iii,ind, ntot, nimp, nposs
! =======================================================================
! output of admis_reg
  INTEGER nrootsf,nrootsd
  DOUBLE PRECISION roots(3),c(0:5),refs(3,3),xo(3),vo(3),E_bound,rootsd(8)
  DOUBLE PRECISION r1,r2,rder0, rdw                 ! for range in rhodot
  DOUBLE PRECISION,DIMENSION(npox) :: frv,rhodot
  DOUBLE PRECISION rmax, rdotmax,rdotmin, dr, drdot ! for grid
! =======================================================================
! MonteCarlo impacts
  INTEGER  obscodn,isucc, nmc
  DOUBLE PRECISION :: simppro,impp
! === multiple minima=================
  INTEGER jmin
  DOUBLE PRECISION chi_min
! ======== output moid =====================
  DOUBLE PRECISION moid0, dnp0, dnm0
! ========= calendar to julian =================
  DOUBLE PRECISION jd,sec
  INTEGER ihr,imin,iy,imo,iday
! ========= input control ======================
  LOGICAL ok ! available data
! file names depending upon run identifier
  CHARACTER*80 run
  CHARACTER*80 titnam
  CHARACTER*100 filname,dummyfile, file,file_orb
  INTEGER le,lnam
  CHARACTER*6 progna
! logical units
  INTEGER iunrms,nrej,nscobs ! for statistical report
! ======== controls and flags ===============
  LOGICAL batch ! batch control
! main menus, choice of arc,covariance required?,copy to/from,test deriv.
  INTEGER ifun,iarc,icov,icop,ider2
  CHARACTER*20 menunam ! characters for menu
  INTEGER, PARAMETER:: iope=0 ! iope=1 experimental; iope=0 distribution vers.
! short circuit to force input of obs and orbit
  LOGICAL init,init2
! asteroids with mass
  LOGICAL found,chi_min_orb
  INTEGER nfound
  INTEGER nd ! number of solve for parameters
  INTEGER ades_type
  CHARACTER*6 vers ! OEF version for output orbit file
!Partials for nd parameters
  DOUBLE PRECISION, DIMENSION(ndimx,ndimx) :: depdep !partials
! ======== loop indexes =====================                           
!  ii=subset of 1,6 etc,j=1,6,i=1,m                                     
  INTEGER j,ii,i 
! ===================================================================== 
! Run name                                                              
  WRITE(*,*) 'Run name =' 
  READ(*,100) run 
100 FORMAT(a) 
  IF(run.eq.'')stop 'No run specified.' 
! input options                                                         
  progna='fitobs' 
  CALL finopt(progna,run,astna0,astnap,error_model,ades_type) 
! error model
  CALL errmod_set(error_model)
! dimension of matrices determined by ls,nls
  nd=6+nls
! version of orbit files determined by ndp
!  IF(det_drp.gt.0)THEN
  vers='OEF2.0'
!  ELSE
!     vers='OEF1.1'
!  ENDIF
! ================fix mole by using 2-body inside Earth?==============
  fix_mole=.true.            ! default for fix_mole
 ! setup close app. output
  tpplane=.true. ! use TP plane, not MTP
! initialization of flag for MOV orbit with minimum chi
  chi_min_orb = .false.
! check for asteroid masses
  IF(lench(astnap).ne.0)THEN
     CALL selpert2(astna0,astnap,nfound)
  ELSE
     CALL selpert(astna0,found)
  ENDIF
! =====================================================================
! initializations
  batch=.false. ! batch control
  sigma0=0.d0  ! central solution is nominal
  mult_adopted=.false.
! asteroid name for messages
  name_obj=astna0
  CALL rmsp(name_obj,lobjnam)
! =====================================================================
! verbosity levels for an interactive program
  verb_pro=10
  verb_clo=10
  verb_dif=20
  verb_mul=10
  verb_rej=20
  verb_io=10
  verb_moid=1
! setting of logical flags: nothing is available at the beginning
  CALL set_state_def
  inide=.false.
  elide=undefined_orbit_elem
  elide%t=0.d0
! intiialization of times, even for non defined elements
  el0=undefined_orbit_elem
  el0%t=0.d0
  elp=undefined_orbit_elem
  elp%t=0.d0
! multiple solutions not yet computed
  imip=0
  imim=0
! unit numbers for orbital elements;zero if not opened
  iunel0=0
  iunelp=0
  iunelt=0
! ================SHORT CIRCUIT ==============================
  init=.true.
  init2=.false.
! ================MAIN MENU===================================
! Choice of function
50 CONTINUE
  IF(init2)THEN
     ifun=2
     GOTO 60
  ENDIF
  IF(init)THEN
     ifun=1
     GOTO 60
  ENDIF
  menunam='mainmenu'
  CALL menu(ifun,menunam,12,'What would you like?=',                    &
     &   'input of observational data=',                                &
     &   'acquire orbital elements=',                                   &
     &   'differential corrections=',                                   &
     &   'first guess for identification=',                             &
     &   'state propagation=',                                          &
     &   'predictions of observations=',                                &
     &   'multiple solutions=',                                         &
     &   'coordinate change=',                                          &
     &   'attributables=',                                              &
     &   'status=',                                                     &
     &   'date conversion=',                                            &
     &   'test derivatives=')
60 IF(ifun.eq.0)THEN
! ==========TERMINATE CLEANLY=========================
     CALL filclo(iun_log,' ')
     CALL filclo(iun_covar,' ')
! close close approach file
     IF(numcla.gt.0)THEN
        CALL filclo(iuncla,' ')
     ELSE
        CALL filclo(iuncla,'DELETE')
     ENDIF
! close propagator parameters file
     CALL filclo(ipirip,'  ')
! close error file
     IF(numerr.gt.0)THEN
        CALL filclo(ierrou,' ')
     ELSE
        CALL filclo(ierrou,'DELETE')
     ENDIF
     IF(obs0)THEN
        CALL filclo(iunrms,' ')
     END IF
     STOP
  ELSEIF(ifun.eq.1)THEN
! ================MENU 1: INPUT OBS============================
     IF(init)THEN
        ifobs=3
        GOTO 61
     ENDIF
     WRITE(*,*)' INPUT OF OBSERVATIONAL DATA'
51   menunam='inputobs'
     CALL menu(ifobs,menunam,5,' which data to input?=',            &
     &      'first arc=','second arc=','both=',                         &
     &      'undiscard outliers, arc 1=',                               &
     &      'undiscard outliers, arc 2=')
61   IF(ifobs.eq.0) GOTO 50
! =====================================================================
! input data, according to request
     IF(ifobs.eq.1.or.ifobs.eq.3)THEN
! =====================================================================
! Arc 1 input
        CALL tee(iun_log,' INPUT OF OBSERVATIONAL DATA, ARC 1=')
        CALL finobs(progna,1,astna0,obs0,nobx,m,obs,obsw, &
      &          rwofi0,error_model,ades_type)
        IF(obs0)THEN
! compose elements  file full name
           IF(iunel0.eq.0)THEN
              elefi0=astna0//'.fel'
              CALL rmsp(elefi0,le)
              CALL filopn(iunel0,elefi0(1:le),'unknown')
! output header
              CALL wromlh2 (iunel0,'ECLM','J2000',vers)
              file=astna0//'.rms'
              CALL rmsp(file,le)
              CALL filopn(iunrms,file(1:le),'UNKNOWN') ! statistical report
           ENDIF
        ELSE
           m=0
        ENDIF
     ENDIF
     IF(ifobs.eq.2.or.ifobs.eq.3)THEN
! =====================================================================
! Arc 2 input
        CALL tee(iun_log,' INPUT OF OBSERVATIONAL DATA, ARC 2=')
        CALL finobs(progna,2,astnap,obsp,nobx-m,mp,obs(m+1:nobx), &
     &   obsw(m+1:nobx),rwofip,error_model,ades_type) 
        IF(obsp)THEN 
! compose elements  file full name  
           IF(iunelp.eq.0)THEN
              elefip=astnap//'.fel'
              CALL rmsp(elefip,le)
              CALL filopn(iunelp,elefip(1:le),'unknown')
! output header
              CALL wromlh2 (iunelp,'ECLM','J2000',vers)
           ENDIF
        ENDIF
     ENDIF
! =====================================================================
! reintroduce outliers
! Arc 1
     IF(ifobs.eq.4)THEN
        obsw(1:m)%sel_coord=1
! Arc 2
     ELSEIF(ifobs.eq.5)THEN
        obsw(m+1:m+mp)%sel_coord=1
     ENDIF
! =====================================================================
!  end input observational data: summary of input status
     obstwo=obs0.and.obsp
     mall=m+mp
     IF(mp.eq.0)obs(m+1)%time_tdt=0.d0
! find if there are radar data
     CALL radar_ob(obs(1:mall)%type,mall)
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
        CALL wromlh2 (iunelt,'ECLM','J2000',vers)
     ENDIF
! =====================================================================
  ELSEIF(ifun.eq.2)THEN
! =====================================================================
! input orbital elements
     IF(init2)THEN
        iele=3
        GOTO 62
     ENDIF
     WRITE(*,*)' INPUT OF ORBITAL ELEMENTS'
! ================MENU 2: INPUT ELEMENTS======================
52   menunam='inputele'
     CALL menu(iele,menunam,11,                                      &
     &      ' Which orbital elements to input/compute?=',               &
     &      ' input arc 1=',' input arc 2=',                            &
     &      ' input both arcs=',                                        &
     &      ' compute arc 1 by Gauss/Vaisala method=',                  &
     &      ' compute arc 2 by Gauss/Vaisala method=',                  &
     &      ' compute both arcs by Gauss/Vaisala=',                     &
     &      ' find multiple solutions of deg. 8 equation, arc 1=',      &
     &      ' find multiple solutions of deg. 8 equation, arc 2=',      &
     &      ' compute arc 1 by Laplace-Poincare=',                      &
     &      ' compute arc 2 by Laplace-Poincare=',                      &
     &      ' give to arc 2 ele of arc 1=')
62   IF(iele.eq.0)GOTO 50
     IF(iele.eq.1.or.iele.eq.3)THEN
! ====================================================================
! initial conditions for arc 1
        CALL tee(iun_log,' INPUT OF ORBITAL ELEMENTS, ARC 1=')
        CALL finele(progna,1,astna0,el0,ini0,cov0,unc0)
        ! Store H and G magnitude from catalog or from eq0 (to be passed to fdiff cor)
        h_ini=el0%h_mag
        g_ini=el0%g_mag
        nd=6+nls
        IF(.not.ini0)THEN
           IF(.not.init2)GOTO 52
        ENDIF
     ENDIF
     IF(iele.eq.2.or.iele.eq.3)THEN
! ===================================================================
! initial conditions for arc 2
        CALL tee(iun_log,' INPUT OF ORBITAL ELEMENTS, ARC 2=')
        CALL finele(progna,2,astnap,elp,inip,covp,uncp)
        IF(.not.inip)THEN
           IF(.not.init2)GOTO 52
        ENDIF
     ENDIF
     IF(iele.ge.4.and.iele.le.6)THEN
! ===================================================================
! use Gauss/Vaisala method for preliminary orbit
        menunam='prelimet'
        CALL menu(imeth,menunam,3,' Which method to use?=',         &
     &            ' Automatic=',                                        &
     &            ' Gauss=',                                            &
     &            ' Vaisala=')
        IF(imeth.eq.1)THEN
           CALL tee(iun_log,' AUTO SELECT METHOD=')
        ELSEIF(imeth.eq.2)THEN
           CALL tee(iun_log,' GAUSS METHOD=')
        ELSEIF(imeth.eq.3)THEN
           CALL tee(iun_log,' VAISALA METHOD=')
        ENDIF
        IF(iele.eq.4.or.iele.eq.6)THEN
! use Gauss/Vaisala method for preliminary orbit, arc1
           IF(.not.obs0)THEN
              WRITE(*,*)'missing observations for arc 1'
              GOTO 52
           ENDIF
           CALL tee(iun_log,' PRELIM. ORB. ARC 1=')
           CALL f_gauss(iunel0,astna0,ini0,cov0,         &
     &           rwofi0,obs,obsw,m,error_model,imeth,el0,ades_type)
        ENDIF
        IF(iele.eq.5.or.iele.eq.6)THEN
! use Gauss/Vaisala method for preliminary orbit, arc 2
           IF(.not.obsp)THEN
              WRITE(*,*)'missing observations for arc 2'
              GOTO 52
           ENDIF
           CALL tee(iun_log,' PRELIM. ORB. ARC 2=')
           CALL f_gauss(iunelp,astnap,inip,covp,rwofip,   &
     &        obs(m+1:m+mp),obsw(m+1:m+mp),mp,error_model,imeth,elp,ades_type)
        ENDIF
     ENDIF
! =====================================================================
! find multiple solutions of deg. 8 equation, and select one
     IF(iele.eq.7)THEN
        IF(.not.obs0)THEN
           WRITE(*,*)'missing observations for arc 1'
           GOTO 52
        ENDIF
        CALL tee(iun_log,' PRELIM. ORB.DEG 8,  ARC 1=')
        CALL f_gaussdeg8(iunel0,astna0,ini0,cov0,          &
     &     rwofi0,obs,obsw,m,error_model,el0)
     ELSEIF(iele.eq.8)THEN
  ! use Gauss/Vaisala method for preliminary orbit, arc 2
        IF(.not.obsp)THEN
           WRITE(*,*)'missing observations for arc 2'
           GOTO 52
        ENDIF
        CALL tee(iun_log,' PRELIM. ORB.DEG 8,  ARC 2=')
        CALL f_gaussdeg8(iunelp,astnap,inip,covp,          &
     &     rwofip,obs(m+1:mall),obsw(m+1:mall),mp,error_model,elp)
     ENDIF
     IF(iele.eq.9)THEN
! use Laplace method for preliminary orbit, arc1
        IF(.not.obs0)THEN
           WRITE(*,*)'missing observations for arc 1'
           GOTO 52
        ENDIF
        CALL tee(iun_log,' PRELIM. ORB. LAPLACE,  ARC 1=')
        CALL laplace_poincare(m,obs,obsw,el3,nroots,nsolp,rr3,failpre,msg,.true.)
        IF(nsolp.eq.1)THEN
           el0=el3(1)
           ini0=.true.
        ELSEIF(nsolp.gt.1)THEN
21         WRITE(*,*)' select one solution between 1 and ', nsolp
           READ(*,*) nselect
           IF(nselect.lt.1.or.nselect.gt.nsolp) GOTO 21
           el0=el3(nselect)
           ini0=.true.
        ELSE
           WRITE(*,*)' no solution '
        ENDIF
     ELSEIF(iele.eq.10)THEN
! use Laplace method for preliminary orbit, arc 2
        IF(.not.obsp)THEN
           WRITE(*,*)'missing observations for arc 2'
           GOTO 52
        ENDIF
        CALL tee(iun_log,' PRELIM. ORB. LAPLACE,  ARC 2=')
        CALL laplace_poincare(mp,obs(m+1:mall),obsw(m+1:mall),el3,nroots,nsolp,rr3,failpre,msg,.true.)
        IF(nsolp.eq.1)THEN
           elp=el3(1)
           inip=.true.
        ELSEiF(nsolp.gt.1)THEN
22         WRITE(*,*)' select one solution between 1 and ', nsolp
           READ(*,*) nselect
           IF(nselect.lt.1.or.nselect.gt.nsolp) GOTO 22
           elp=el3(nselect)
           inip=.true.
        ELSE
           WRITE(*,*)' no solution '
        ENDIF
     ENDIF
     IF(iele.eq.11)THEN
! =====================================================================
! copy elements of arc 1 into elements of arc 2
        IF(ini0)THEN
           elp=el0
           inip=.true.
! covariance is not copied
           covp=.false.
        ELSE
           WRITE(*,*)' initial conditions for arc 1 not available'
        ENDIF
     ENDIF
! initialisation of Yarkovsky after acquiring elements
! (the physical model of the first asteroid is assumed)
     IF(ini0) CALL yarkinit(astna0,el0)
! =====================================================================
  ELSEIF(ifun.eq.3)THEN
     CALL tee(iun_log,' DIFFERENTIAL CORRECTIONS=')
! =================MENU 3: DIFFERENTIAL CORRECTIONS====================
53   CALL orb_sel2(.false.,iarc)
     IF(iarc.eq.0)GOTO 553
     CALL obs_cop(1,iarc) ! copy observations to obsc, obswc
     IF(iarc.eq.1)THEN
        rwofic=rwofi0 ! ARC 1
        iunelc=iunel0
     ELSEIF(iarc.eq.2)THEN
        rwofic=rwofip ! ARC 2
        iunelc=iunelp
     ELSEIF(iarc.eq.3)THEN
        rwofic=rwotwo ! ARCS IDENTIFIED
        iunelc=iunelt
     ENDIF
! =========== select mode======================
     menunam='difcomod'
553  CALL menu(idif,menunam,5,' select correction and reject mode=',&
     &      'correct all, autoreject=',                                 &
     &      'correct all, no rejections=',                              &
     &      'constrained solution on LOV (no rejections)=',             &
     &      'correct only some elements (no rejections)=',              &
     &      'compute residuals and covariance (no correction)=')
     IF(idif.eq.0) GOTO 50
     itsav=itmax
     IF(idif.eq.5)itmax=0
! auto/manual reject
     IF(idif.eq.1)THEN
        autrej=.true.
     ELSE
        autrej=.false.
     ENDIF
! which elements to correct
     IF(idif.ne.4)THEN
        interactive=0
     ELSE
        interactive=1
     ENDIF
! ====================constrained solution===========================
     IF(idif.eq.3)THEN
        incfit=5
! check availability of observations and initial condition
        CALL cheobs(obsflag,inic,ok)
        IF(.not.ok) GOTO 53
! check availability of JPL ephemerides and ET-UT table
        CALL chetim(obsc(1)%time_tdt,obsc(m)%time_tdt,ok)
        IF(.not.ok) GOTO 53
681     WRITE(*,*)' use scaling, 1=yes, 0=no?'
        READ(*,*)iscal
        IF(iscal.eq.0)THEN
           scaling_lov=.false.
        ELSEIF(iscal.eq.1)THEN
           scaling_lov=.true.
        ELSE
           WRITE(*,*)' answer ', iscal, ' not understood'
           GOTO 681
        ENDIF
682     WRITE(*,*)' which LOV, 1=largest eigenv., 2=second'
        READ(*,*)iscal
        IF(iscal.eq.1)THEN
           second_lov=.false.
        ELSEIF(iscal.eq.2)THEN
           second_lov=.true.
        ELSE
           WRITE(*,*)' answer ', iscal, ' not understood'
           GOTO 682
        ENDIF
        CALL tee(iun_log,' CONSTRAINED DIFFERENTIAL CORRECTIONS=')
        peq=0.d0
        CALL constr_fit(mc,obsc,obswc,elc,peq,elc,uncc,csinoc,delnoc,rmshc,iobc,succ,nd,h_ini,g_ini)
        covc=succ
! =========================step-fit=================================
        IF(succ)THEN
           WRITE(*,*)' step_fit along the LOV? 1=yes 0=no'
           READ(*,*) ist_fit
           IF(ist_fit.eq.1)THEN
              CALL step_fit2(mc,obsc,obswc,elc,uncc,csinoc,delnoc,iobc,succ)
           ENDIF
        ENDIF
! =========================full solution=============================
     ELSE
        incfit=6
        CALL tee(iun_log,' FULL DIFFERENTIAL CORRECTIONS=')
        CALL fdiff_cor(batch,1,obsflag,inic,ok,covc,elc,mc,obsc,obswc,iobc,    &
             &         rwofic,elc,uncc,csinoc,delnoc,rmshc,succ,H_READ=h_ini,G_READ=g_ini,ADES_TYPE=ades_type)
! output for global parameter determination
        nrej=0
        nscobs=0
        DO j=1,mc
           IF(obsc(j)%type.eq.'O')THEN
              nscobs=nscobs+2
              IF(obswc(j)%sel_coord.eq.0)nrej=nrej+2
           ELSEIF(obsc(j)%type.eq.'R'.or.obsc(j)%type.eq.'V')THEN
              nscobs=nscobs+1
              IF(obswc(j)%sel_coord.eq.0)nrej=nrej+1
           ENDIF
        ENDDO
        WRITE(iunrms,222)astnac,nscobs,nrej,nscobs-nrej,csinoc,csinoc**2*(nscobs-nrej)
222     FORMAT(A9,1X,I5,1X,I4,1X,I5,1X,F10.7,1X,F14.9)
     ENDIF
! availability of observations, initial condition, JPL and ET-UT data
     IF(.not.ok.or..not.succ) GOTO 50
! output new elements
     CALL write_elems(elc,astnac,'ML',dyn,nls,ls,UNC=uncc,UNIT=iunelc)
     CALL nomoid(elc%t,elc,moid0,dnp0,dnm0)
     write(*,199)moid0,0,dnp0,dnm0
199  format('orb.dist.      dist.n+  dist.n-'/                &
          &              f8.5,1x,i4,1x,f8.5,1x,f8.5)
     write(*,*)
     write(iunelc,198)moid0,0,dnp0,dnm0
198  format('!MOID ',f8.5,1x,i4/'!NODES ',f8.5,1x,f8.5)
! restore state
     icop=2
     CALL sta_cop(icop,iarc)
     CALL obs_cop(icop,iarc)
     itmax=itsav
! =====================================================================
  ELSEIF(ifun.eq.4)THEN
! =====================================================================
! check availability of JPL ephemerides and ET-UT table
     CALL chetim(obs(1)%time_tdt,obs(mall)%time_tdt,ok)
     IF(.not.ok) GOTO 50
! ok, go on with arc
     CALL tee(iun_log,' FIRST GUESS FOR IDENTIFICATION=')
! ==================MENU 4: FIRST GUESS=======================
54   menunam='firstgue'
     CALL menu(igue,menunam,4,' which initial guess?=',             &
     &      'use averages, fit longitudes=',                        &
     &      'use elements of arc 1=',                               &
     &      'use elements of arc 2=',                               &
     &      'recompute ident=')
     IF(igue.eq.0)GOTO 50
! =====================================================================
     IF(igue.eq.1)THEN
! =====================================================================
! first guess of a unique solution for both arcs, with averages
!   and fit to longitudes
!  =====================================================================
! check availability of observations and initial condition
        iniboth=ini0.and.inip
        CALL cheobs(obstwo,iniboth,ok)
        IF(.not.ok) GOTO 54
! can be done
        CALL tee(iun_log,' GUESS FROM LONGITUDE FIT=')
! =====================================================================
! epoch time in the middle, unless they are too close
        IF(abs(el0%t-elp%t).lt.1.d0)THEN
           WRITE(*,*)' initial times ',el0%t,elp%t,' too close'
           GOTO 54
        ENDIF
! =====================================================================
! magnitude is the one of the first (not knowing any better)
        el=el0
        el%t=(el0%t+elp%t)/2.d0
        WRITE(iun_log,123) el%t
        WRITE(*,123) el%t
123     FORMAT(1x,'tm =',f8.1)
        WRITE(iun_log,*)
! =====================================================================
! use as first guess linear interpolation for h,k,p,q
        el%coord(2:5)=(el0%coord(2:5)+elp%coord(2:5))/2.d0
! Estimate of the mean motion (hence semimajor axis)
! allowing to combine two arcs
        CALL start(el0,elp,1,ng,enm,el%coord(1),el%coord(6),ok)
        IF(ok)THEN
           initwo=.true.
           CALL wriequ(iun_log,astna0,el%t,el%coord)
           WRITE(*,*)' INITIAL COND. AVAILABLE, NOW USE DIFF. CORR.'
        ENDIF
! =====================================================================
     ELSEIF(igue.eq.2)THEN
! =====================================================================
        CALL cheobs(obstwo,ini0,ok)
        IF(.not.ok) GOTO 54
! can be done
        CALL tee(iun_log,' USE ELEMENTS OF ARC 1=')
        el=el0 ! magnitude from first arc
        initwo=.true.
        WRITE(*,*)' INITIAL COND. AVAILABLE, NOW USE DIFF. CORR.'
! =====================================================================
     ELSEIF(igue.eq.3)THEN
! =====================================================================
        CALL cheobs(obstwo,inip,ok)
        IF(.not.ok) GOTO 54
! can be done
        CALL tee(iun_log,' USE ELEMENTS OF ARC 2=')
        el=elp ! magnitude from second arc
        initwo=.true.
        WRITE(*,*)' INITIAL COND. AVAILABLE, NOW USE DIFF. CORR.'
! =====================================================================
     ELSEIF(igue.eq.4)THEN
! =====================================================================
! recompute first guess with the identif algorithm
! =====================================================================
!        IF(el0%coo.ne.'EQU'.or.elp%coo.ne.'EQU')THEN
!           WRITE(*,*)'not possible unless EQU, given ', el0%coo, elp%coo
!           GOTO 54
!        ENDIF
        CALL fident(6,cov0,covp,el0,elp,unc0,uncp,elide,ff)
        IF(ff)THEN
           el=elide
           initwo=.true.
           WRITE(*,*)' INITIAL COND. AVAILABLE, NOW USE DIFF. CORR.'
        ELSE
           GOTO 54
        ENDIF
! =====================================================================
     ENDIF
! =====================================================================
  ELSEIF(ifun.eq.5)THEN
! =====================================================================
! Computation of asteroid elements at a required time
     CALL tee(iun_log,' PROPAGATION OF ELEMENTS=')
! ================= MENU 5: PROPAGATION =======================
55   menunam='propagat'
     CALL menu(iprop,menunam,8,'what to propagate?=',               &
     &      'propagate arc 1, variable time=',                      &
     &      'propagate arc 2, variable time=',                      &
     &      'propagate identified orbit, variable time=',           &
     &      'propagate arc 1, weighted center of observations=',    &
     &      'propagate arc 2, weighted center of observations=',    &
     &      'generate time history of orbital elements, arc 1=',    &
     &      'generate time history of orbital elements, arc 2=',    &
     &   'generate time history of orbital elements, identified orbit=')
     IF(iprop.eq.0) GOTO 50
! =====================================================================
! state vector only? also covariance? (not meaningful for ephemerides)
! =====================================================================
     IF(iprop .le. 5)THEN
        menunam='null'
        CALL menu(icov,menunam,2,'What is required?=',              &
     &           'orbital elements only=',                          &
     &           'also covariance matrix=')
        IF(icov.eq.0) GOTO 55
     ELSE
! default no covariance
        icov=1
     ENDIF
     IF(iprop.le.5)THEN
! =====================================================================
! propagation
! =====================================================================
        IF(iprop.eq.4)THEN ! selection of target epoch
           tr= meanti(obs(1:m)%time_tdt,obsw(1:m)%rms_coord(1),        &
     &           obsw(1:m)%rms_coord(2),obs(1:m)%coord(2),m)
        ELSEIF(iprop.eq.5)THEN
           tr= meanti(obs(m+1:m+mp)%time_tdt,obsw(m+1:m+mp)%rms_coord(1),&
     &           obsw(m+1:m+mp)%rms_coord(2),obs(m+1:m+mp)%coord(2),mp)
! warning: funny result if mp=0; check obsp?
        ELSE
           WRITE(*,*)' propagate to epoch (MJD)?   '
           READ(*,*)tr
        ENDIF
        IF(iprop.eq.1.or.iprop.eq.4)THEN
           iarc=1
           iunelc=iunel0
        ELSEIF(iprop.eq.2.or.iprop.eq.5)THEN
           iarc=2
           iunelc=iunelp
        ELSEIF(iprop.eq.3)THEN
           iarc=3
           iunelc=iunelt
        ENDIF
        CALL orb_sel2(.true.,iarc)
        IF(icov.eq.2)THEN
           CALL cov_avai(uncc,elc%coo,elc%coord)
        ENDIF
        CALL fstpro(.false.,icov,inic,covc,iun_log,iun_covar,ok,       &
     &         elc,uncc,tr,elc,uncc)
        IF(ok.and.(.not.kill_propag))THEN
! problem with name for identification; used arc 1, but...
           IF(dyn%ndp.gt.0)THEN
              CALL write_elems(elc,astnac,'ML',dyn,nls,ls,UNC=uncc,FILE=dummyfile,UNIT=iunelc)
           ELSE
              CALL write_elems(elc,astnac,'ML',UNC=uncc,FILE=dummyfile,UNIT=iunelc)
           END IF
           CALL nomoid(elc%t,elc,moid0,dnp0,dnm0)
           write(*,199)moid0,0,dnp0,dnm0
           write(iunelc,198)moid0,0,dnp0,dnm0
! copy into the right state
           CALL sta_cop(2,iarc)
        ENDIF
        GOTO 50
     ELSEIF(iprop.gt.5.and.iprop.le.8)THEN
! =====================================================================
! ephemerides (orbital elements) generation
! =====================================================================
        WRITE(*,*)' Current time is : ',el0%t,'(MJD).'
        WRITE(*,*)' begin ephemerides from epoch (MJD)?   '
        READ(*,*)tr
        WRITE(*,*)' end ephemerides from epoch (MJD)?   '
        READ(*,*)tf
        WRITE(*,*)' time step in days?'
        READ(*,*)step
        WRITE(*,*) 'Is data correct? (y/n)'
        READ(*,*)ans
        IF(ans(1:1).eq.'n' .or. ans(1:1).eq.'N') GOTO 55
!   determine total number of steps (with margin)
        interv=tf-tr
        numsav=interv/step+10
! determine type of elements output
        menunam='coord'
        call menu(icoord,menunam,5,'What type of elements?=',         &
     &           'KEPlerian=', 'EQUinoctial=',                        &
     &           'CARtesian=', 'COMetary=','ATTributable=')
        IF(icoord.eq.0)THEN
           goto 55
        ELSE
           cooy=cootyp(icoord)
        ENDIF
        IF(iprop.eq.6)THEN
           iarc=1
        ELSEIF(iprop.eq.7)THEN
           iarc=2
        ELSEIF(iprop.eq.8)THEN
           astnaj=astna0//'joint'
           iarc=3
        ENDIF
        CALL orb_sel2(.true.,iarc)
! avoiding self perturbations
        CALL selpert(astna0,found)
        IF(found)THEN
           WRITE(*,*)astna0,' massive asteroid, no selfperturbation'
        ENDIF
        CALL fsteph(astnac,'.',inic,ok,elc,                &
     &           tr,tf,step,numsav,.true.,cooy,.true.)
! if some data are not available, this cannot be done
        IF(.not.ok)THEN
           WRITE(iun_log,*)'    DATA NOT AVAILABLE'
           GOTO 55
        ENDIF
     ENDIF
! =====================================================================
  ELSEIF(ifun.eq.6)THEN
! =====================================================================
! prediction of observations
     CALL tee(iun_log,' PREDICTION OF OBSERVATIONS=')
! MENU 6: PREDICTIONS
56   CALL orb_sel2(.false.,iarc)
     IF(iarc.eq.0) GOTO 50
! setup title string for graphics output
     CALL titast(iarc,astna0,astnap,titnam,filname,lnam)
! =====================================================================
! observations vector only? also covariance?
! ==================================================================
! move to fobpre???
     menunam='predicbd'
     CALL menu(iprob,menunam,6,'What is required?=',                &
     &      'observations (alpha, delta) only=',                    &
     &      'use simulated observation=',                           &
     &      'also covariance matrix=',                              &
     &      'confidence boundary=',                                 &
     &      'compare CB with observations=',                        &
     &      'ephemerides (on the sky)=')
     ok=.true.
     IF(iprob.eq.0)GOTO 50
! ===================================================================
! assign observation time
556  IF(iprob.le.5)THEN
        CALL asstim(iprob,obs(1:mall)%type,obs(1:mall)%time_tdt,    &
     &       obs(1:mall)%time_utc,obs(1:mall)%obscod_i,m,mall,im,   &
     &       type1,t1,tut1,ids)
     ELSE
        CALL seleph(tut1,t1,tut2,t2,dt,ids)
        im=1
     ENDIF
! ===================================================================
! predict
! =====================================================================
     CALL fobpre(iprob,inic,covc,ok,                &
     &           titnam,filname,elc,uncc,ids,type1,t1,        &
     &           tut1,obs(im)%coord(1),obs(im)%coord(2),t2,dt,astnac)
! if some data are not available, this cannot be done
     IF(.not.ok)THEN
        WRITE(iun_log,*)'    DATA NOT AVAILABLE'
        GOTO 56
     ENDIF
     IF(iprob.eq.2)THEN
        WRITE(*,*)'more simulated observations? 0=no'
        READ(*,*)iyes
        IF(iyes.ne.0)GOTO 556
     ENDIF
! =====================================================================
  ELSEIF(ifun.eq.7)THEN
! =====================================================================
! search for alternate solutions
     CALL tee(iun_log,'MULTIPLE SOLUTIONS=')
! MENU 7: MULTIPLE SOLUTIONS
! choice of arc
58   menunam='multisol'
     CALL menu(marc,menunam,5,'which orbit?=',                   &
     &      'arc 1=','arc 2=',                                   &
     &      'joint computed orbit=',                             &
     &      'use already computed=',                             &
     &      'input from file=')
     IF(marc.eq.0) GOTO 50
! choice of LOV
581  WRITE(*,*)' use scaling, 1=yes, 0=no, 2=linear LOV?'
     READ(*,*)iscal
     IF(iscal.eq.0)THEN
        scaling_lov=.false.
     ELSEIF(iscal.eq.1)THEN
        scaling_lov=.true.
     ELSEIF(iscal.eq.2)THEN
! note that scatter plane has no meaning if not lin_lov
        lin_lov=.true.
583     WRITE(*,*)' use scatterplane? 0=no 1=yes'
        READ(*,*)iscat
        IF(iscat.eq.0)THEN
           scatterplane=.false.
        ELSEIF(iscat.eq.1)THEN
           scatterplane=.true.
           WRITE(*,*) 'date MJD before scattering encounter'
           READ(*,*) before_scatter
           WRITE(*,*) 'date MJD after scattering encounter'
           READ(*,*) after_scatter
        ELSE
           WRITE(*,*)' answer ', iscat, ' not understood'
           GOTO 583
        ENDIF
     ELSE
        WRITE(*,*)' answer ', iscal, ' not understood'
        GOTO 581
     ENDIF
582  IF(iscal.lt.2)THEN
        WRITE(*,*)' which LOV, 1=largest eigenv., 2=second'
        READ(*,*)iscal
        IF(iscal.eq.1)THEN
           second_lov=.false.
        ELSEIF(iscal.eq.2)THEN
           second_lov=.true.
        ELSE
           WRITE(*,*)' answer ', iscal, ' not understood'
           GOTO 582
        ENDIF
     END IF
584  WRITE(*,*)' which sampling? 1=const step in sigma 2=const step in IP'
     READ(*,*)isampl
     IF(isampl.eq.1)THEN
! Constant step in sigma
        prob_sampl=.false.
        WRITE(*,*)' sigma_max (extent on LoV)?'
        READ(*,*) sigma
559     WRITE(*,*)' steps on each side'
        READ(*,*) imult
        WRITE(iun_log,*)' sigma=',sigma,'  in ',imult,' steps'
        IF(2*imult+1.gt.mulx)THEN
           WRITE(*,*)' too many; max is ',mulx
           GOTO 559
        ENDIF
     ELSEIF(isampl.eq.2)THEN
! Constant step in IP
        prob_sampl=.true.
        WRITE(*,*)' prob_step (completeness limit)?'
        READ(*,*) prob_step
        WRITE(*,*)' sigma_step_max (control against too long steps)?'
        READ(*,*) sigma_step_max
        WRITE(*,*)' sigma_max (extent on LoV)?'
        READ(*,*) sigma_max
        WRITE(iun_log,*) ' sigma_max=',sigma_max,'  prob_step ',prob_step,' and sigma_step_max', sigma_step_max
     ELSE
        WRITE(*,*)' answer ', isampl, ' not understood'
        GOTO 584
     ENDIF
! =====================================================================
! compute multiple solutions
! =====================================================================
     IF(marc.ge.1.and.marc.le.3)THEN
        CALL orb_sel2(.true.,marc)
        CALL obs_cop(1,marc) ! copy observations to obsc, obswc
        CALL tee(iun_log,'COMPUTE MULTIPLE SOL=')
        IF(lin_lov)THEN
           IF(scatterplane)THEN
              CALL filnam('.',astnac,'wdir',file,le)
              CALL filopn(iwdir,file(1:le),'unknown')
           ENDIF
! use linear lov
           IF(prob_sampl)THEN
              CALL lin_multi_gen(batch,elc,uncc,prob_step,sigma_step_max,sigma_max,nd,imult,iwdir,mc,obs,obsw,csinoc)
           ELSE
              CALL lin_multi(batch,elc,uncc,sigma,imult,nd,iwdir,mc,obsc,obswc,csinoc)
           END IF
           IF(scatterplane) CALL filclo(iwdir,' ')
        ELSEIF(mult_adopted)THEN
! use alternate nominal
           IF(prob_sampl)THEN
              CALL multi_gen(batch,elc,uncc,csinoc,delnoc,mc,obsc,obswc,   &
                   &          prob_step,sigma_step_max,sigma_max,nd,imult,CSINO0=csino00)
           ELSE
              CALL f_multi(batch,obsflag,inic,ok,covc, &
                   &   elc,uncc,csinoc,delnoc,mc,obsc,obswc,sigma,imult,nd,CSINO0=csino00)
           ENDIF
        ELSE
! use of true nominal
           IF(prob_sampl)THEN
              CALL multi_gen(batch,elc,uncc,csinoc,delnoc,mc,obsc,obswc,   &
                   &          prob_step,sigma_step_max,sigma_max,nd,imult)
           ELSE
              CALL f_multi(batch,obsflag,inic,ok,covc, &
                   &   elc,uncc,csinoc,delnoc,mc,obsc,obswc,sigma,imult,nd)
           ENDIF
        ENDIF
     ELSEIF(marc.eq.4)THEN
        CALL tee(iun_log,'USE THE ALREADY COMPUTED ONES=')
        ok=imip-imim.gt.0
! input from file
     ELSEIF(marc.eq.5)THEN
        CALL tee(iun_log,'INPUT MULTIPLE SOL. FROM FILE=')
        WRITE(*,*)' File name? no suffix'
        READ(*,'(A)')catname
        CALL mult_input(catname,ok,sigma,nd,sigma_step_max)
        CALL orb_sel2(.false.,iarc) ! need to know which arc they are for
        CALL obs_cop(1,iarc) ! copy observations to obsc, obswc
     ENDIF
     IF(.not.ok)GOTO 58
     nmult=imip-imim+1
     WRITE(*,*)' number of multiple sol. available', nmult
     IF(nmult.le.0)THEN
        CALL tee(iun_log,'FAILED MULTIPLE SOLUTIONS')
        GOTO 58
     ENDIF
! ======================================================
! how to use multiple solutions?
! ======================================================
145  menunam='multiuse'
     CALL menu(ifff,menunam,5,'what to do?=',                       &
     &      'plot of multiple solutions=',                          &
     &      'multiple predicted observations=',                         &
     &      'adopt one of the above solution=',                         &
     &      'propagate multiple solutions=',                            &
     &      'close approach analysys=')
     IF(ifff.eq.0) GOTO 50
! =================================================================
     CALL orb_sel2(.false.,iarc)
! setup title string for graphics output
     CALL titast(iarc,astna0,astnap,titnam,filname,lnam)
! ======================================================
     IF(ifff.eq.1)THEN
! ======================================================
! a-e plot of multiple solutions
        IF(marc.ge.4) elc=elm(imi0)
        CALL fmuplo(titnam,sigma)
! ======================================================
     ELSEIF(ifff.eq.2)THEN
! ======================================================
! multiple predicted observation
        menunam='null'
        CALL menu(iff,menunam,2,'What is required?=',               &
     &      'predicted observations (alpha, delta) only=',          &
     &      'compare with actual observations=')
        IF(iff.eq.0) GOTO 145
! =====================================================================
! assign observation time
        IF(iff.eq.2)iffat=5
        IF(iff.eq.1)iffat=3
        CALL asstim(iff+2,obs%type,obs%time_tdt,obs%time_utc,           &
     &           obs%obscod_i,m,mall,im,type1,t1,tut1,ids)
! =====================================================================
! check availability of JPL ephemerides and ET-UT table
        CALL chetim(t1,t1,ok)
        IF(.not.ok) GOTO 145
! can be done
        CALL tee(iun_log,'MULTIPLE PREDICTED OBSERVATIONS=')
! =================================================
! only alpha, delta, magnitude
        CALL fmuobs(type1,ids,t1,tut1,sigma,     &
     &      obs(im)%coord(1),obs(im)%coord(2),iff, &
     &            titnam,filname,iun_log)
! ======================================================
     ELSEIF(ifff.eq.3)THEN
! ======================================================
! adopt alternate solution
        IF(prob_sampl)THEN
           WRITE(*,*)' densification not available for IP sampling'
           GOTO 145
        ENDIF
        WRITE(*,*)' which one to keep? 0=none'
        READ(*,*) imi
        IF(imi.ge.imim.and.imi.le.imip)THEN
           CALL tee(iun_log,'ALTERNATE SOLUTION ADOPTED=')
           IF(mult_adopted)THEN
              WRITE(*,*)' WARNING: multiple densification is'
              WRITE(*,*) ' incompatible with impact monitoring'
           ENDIF
           mult_adopted=.true.
           sigma0=(imi-imi0)*delta_sigma
           WRITE(*,*)' NUMBER ',imi
           WRITE(iun_log,*)' NUMBER ',imi
           WRITE(iun_log,*)
           WRITE(*,144)imi,elm(imi)%coord,csinom(imi)
144        FORMAT(i5,1P,6d15.8,e13.5,e12.3)
! copy state vector and matrices, norms
           icop=2
! handle case in which the multiple solutions do not come from arc
! CORRECTION AMC 28/10/2013 iarc-> marc
           IF(marc.eq.4)THEN
147           WRITE(*,*)' for which arc? 1,2=arcs, 3=joint'
              READ(*,*)iarm
              IF(iarm.ge.1.and.iarm.le.3)THEN
                 iarc=iarm
              ELSE
                 WRITE(*,*)' must be 1,2,3'
                 GOTO 147
              ENDIF
              csino00=csinom(imi0)
           ELSEIF(marc.eq.5)THEN
              csino00=csinom(imi0)
           ENDIF
           IF(iarc.eq.1)THEN
              CALL stacop(icop,el0,unc0,csino0,delno0,             &
     &               elm(imi),unm(imi),csinom(imi),delnom(imi))
           ELSEIF(iarc.eq.2)THEN
              CALL stacop(icop,elp,uncp,csinop,delnop,             &
     &               elm(imi),unm(imi),csinom(imi),delnom(imi))
           ELSEIF(iarc.eq.3)THEN
              CALL stacop(icop,el,unc,csinor,delnor,                &
     &               elm(imi),unm(imi),csinom(imi),delnom(imi))
           ELSE
              WRITE(*,*)' iarm=',iarm,' not understood'
              GOTO 145
           ENDIF
! handle dynamical parameters
! WARNING: only one set of dyn.par. for arcs 1,2,3
           IF(dyn%ndp.gt.0)THEN
              dyn%dp(1:ndyx)=dpm(1:ndyx,imi)
           ENDIF
        ELSE
           CALL tee(iun_log,'BACK TO THE ORIGINAL SOLUTION=')
        ENDIF
     ELSEIF(ifff.eq.4)THEN
! =====================================================================
! select propagation time
        CALL tee(iun_log,'MULTIPLE PROPAGATION=')
        tcmult=elm(imi0)%t
        WRITE(*,*)' Current time is : ',tcmult,'(MJD).'
        WRITE(*,*)' propagate to epoch (MJD)?   '
        READ(*,*)trmult
! check availability of JPL ephemerides
        CALL chetim(tcmult,trmult,ok)
        IF(.not.ok)GOTO 145
! propagation
        IF(scatterplane) THEN
           WRITE(*,*) 'FITOBS: scatter_stretch = TRUE'
           scatter_stretch=.TRUE.
           wfile=astnac//'.wdir'
           CALL rmsp(wfile,le)
           CALL filopn(iwdir,wfile(1:le),'old')
! read quadruple precision, convert in double
           READ(iwdir,*) vvv_q(1:nd)
           vvv_tp(1:nd)=vvv_q(1:nd)
           CALL filclo(iwdir,' ')
        END IF
        CALL fmupro(iun_log,trmult)
! close approach analysis on multiple solutions (as in CLOMON2)
     ELSEIF(ifff.eq.5)THEN
        CALL tee(iun_log,'CLOSE APPROACH ANALYSIS=')
        ddocd1=ddocd
        ll=lench(ddocd1)
        helpfi=ddocd1(1:ll)//'/'//'closapan'
        CALL rmsp(helpfi,ll)
        helpfi=helpfi(1:ll)//'.help'
        CALL filopn(iunit,helpfi,'OLD')
        CALL filcat(iunit)
        CALL filclo(iunit,' ')
        filname=run//'.vis'
        CALL rmsp(filname,le)
        CALL filopn(iunvi,filname(1:le),'UNKNOWN')
! select interval in index (and sigma) space
146     WRITE(*,*)' SELECT INTERVAL, between ',imim,' and ',imip
        READ(*,*)m1,m2
        IF(m1.lt.imim.or.m2.gt.imip.or.m1.ge.m2)THEN
           WRITE(*,*)'must be ',imim,' <= ',m1,' < ',m2,' <= ',imip
           WRITE(*,*)' try again'
           GOTO 146
        ENDIF
        WRITE(iun_log,*)' VA interval between m1=',m1,' m2=',m2
! select final time
        WRITE(*,*)' search for close approaches until time (MJD)?'
        READ(*,*) tmcla
        WRITE(iun_log,*)' search until time',tmcla,' MJD'
! Control on stretching?
        WRITE(*,*)' Stretching limit (0 for not using it)?'
        READ(*,*) limit_stretch
! modified for densification
        IF(mult_adopted)THEN
           CALL fclomon2(progna,astnac,mc,obsc,obswc,m1,m2,tmcla,sigma,nd,limit_stretch,sigma0)
        ELSE
           CALL fclomon2(progna,astnac,mc,obsc,obswc,m1,m2,tmcla,sigma,nd,limit_stretch,0.d0)
        ENDIF
        IF(num_vi.le.0)THEN
           CALL filclo(iunvi,'DELETE')
           GOTO 145
        ENDIF
        DO nvi=1,num_vi
           CALL wri_tppoint(vis(nvi)%tp,iunvi,.true.)
           IF(dyn%ndp.gt.0)THEN
              CALL write_elems(vis(nvi)%ele,astnac,'ML',dyn,nls,ls,UNC=vis(nvi)%unc,FILE=dummyfile,UNIT=iunvi)
           ELSE
              CALL write_elems(vis(nvi)%ele,astnac,'ML',UNC=vis(nvi)%unc,FILE=dummyfile,UNIT=iunvi)
           END IF
        ENDDO
        CALL filclo(iunvi,' ')
! what to do with VIs
        IF(iope.eq.0) GOTO 145
345     CONTINUE
        DO nvi=1,num_vi
           WRITE(*,*)nvi,vis(nvi)%tp%tcla, vis(nvi)%tp%b, vis(nvi)%tp%txi, vis(nvi)%tp%tze
        ENDDO
        WRITE(*,*)' which VI to analyse? 0=none'
        READ(*,*) nvi
        IF(nvi.le.0) GOTO 145
        curr_vi=vis(nvi)
        CALL vi_draw(del)
        GOTO 345
     ENDIF
! stay inside the multiple orbits case for repeated use of the data
     GOTO 145
! =====================================================================
  ELSEIF(ifun.eq.8)THEN
! =====================================================================
! coordinate changes
     CALL tee(iun_log,'COORDINATE CHANGES=')
     menunam='coord'
     CALL menu(icoord,menunam,6,'Coordinates?=',          &
     &      'KEPlerian=',                                   &
     &      'EQUinoctal=',                                  &
     &      'CARtesian=',                                   &
     &      'COMetary=',                                    &
     &      'COmetary True anomaly=',                       &
     &      'ATTributables=')
     ok=.true.
     IF(icoord.eq.0)GOTO 50
     cooy=cootyp(icoord)
     menunam='cooarc'
     CALL menu(iarc,menunam,4,                               &
     &      ' Which orbital elements to convert?=',          &
     &      ' arc 1=',' arc 2=',' both arcs=',               &
     &      ' identification=')
     IF(iarc.eq.1.or.iarc.eq.3)THEN
        IF(ini0)THEN
! handle special case of attributable elements
           IF(same_station(obs,m))THEN
              obscodi=obs(1)%obscod_i
           ELSE
              obscodi=500
           ENDIF
           IF(cov0)THEN
              CALL coo_cha(el0,cooy,el0,fail_flag,dee)
              depdep(1:6,1:6)=dee
              IF(nd.gt.6)THEN
                 depdep(7:nd,1:nd)=0.d0
                 depdep(1:nd,7:nd)=0.d0
                 DO j=7,nd
                    depdep(j,j)=1.d0
                 ENDDO
              ENDIF
              IF(nd.ne.unc0%ndim)THEN
                 WRITE(*,*)' Warning '
                 WRITE(*,*)' propagunc: inconsistent dimension of matrices, nd',nd
                 WRITE(*,*)' unc0%ndim: ', unc0%ndim
260              WRITE(*,*)' do you want to continue? (y=yes or n=no)'
                 READ(*,*) choice
                 IF(choice.EQ.'y')THEN
                    CALL propagunc(unc0%ndim,unc0,depdep(1:unc0%ndim,1:unc0%ndim),unc0)
                 ELSE IF(choice.EQ.'n') THEN
                    WRITE(*,*) 'covariance non propagated'
                    unc0%succ=.FALSE.
                 ELSE
                    WRITE(*,*) 'Wrong choice, please select only y or n.'
                    GO TO 260
                 END IF
              ELSE
                 CALL propagunc(nd,unc0,depdep(1:nd,1:nd),unc0)
              END IF
           ELSE
              CALL coo_cha(el0,cooy,el0,fail_flag,OBSCODE=obscodi)
              unc0%succ=.false.
           ENDIF
           IF(mov_flag.AND.chi_min_orb)THEN
              CALL write_elems(el0,astna0,'ML',UNC=unc0,FILE=dummyfile)
              chi_min_orb = .FALSE.
           ELSE
              IF(dyn%ndp.gt.0.AND.nd.EQ.unc0%ndim)THEN
                 CALL write_elems(el0,astna0,'ML',dyn,nls,ls,UNC=unc0,FILE=dummyfile,UNIT=iunel0)
              ELSE
                 CALL write_elems(el0,astna0,'ML',UNC=unc0,FILE=dummyfile,UNIT=iunel0)
              END IF
           END IF
! write .ke0 file as in neofit
           IF(cooy.EQ.'KEP') THEN
              CALL filnam('.',astna0,'ke',elefil,le)
              el0%mag_set=.TRUE.
              CALL wrikep(elefil(1:le),astna0,el0,unc0,moid,pha)
           END IF
           WRITE(*,*)' elements for arc 1', el0
        ELSE
           WRITE(*,*)' initial conditions not available for ',astna0
        ENDIF
     ENDIF
     IF(iarc.eq.2.or.iarc.eq.3)THEN
        IF(inip)THEN
! handle special case of attributable elements
           IF(same_station(obs(m+1:m+mp),mp))THEN
              obscodi=obs(m+1)%obscod_i
           ELSE
              obscodi=500
           ENDIF
           IF(covp)THEN
              CALL coo_cha(elp,cooy,elp,fail_flag,dee)
              depdep(1:6,1:6)=dee
              IF(nd.gt.6)THEN
                 depdep(7:nd,1:nd)=0.d0
                 depdep(1:nd,7:nd)=0.d0
                 DO j=7,nd
                    depdep(j,j)=1.d0
                 ENDDO
              ENDIF
              CALL propagunc(nd,uncp,depdep(1:nd,1:nd),uncp)
           ELSE
              CALL coo_cha(elp,cooy,elp,fail_flag,OBSCODE=obscodi)
              uncp%succ=.false.
           ENDIF
           IF(dyn%ndp.gt.0)THEN
              CALL write_elems(elp,astnap,'ML',dyn,nls,ls,UNC=uncp,FILE=dummyfile,UNIT=iunelp)
           ELSE
              CALL write_elems(elp,astnap,'ML',UNC=uncp,FILE=dummyfile,UNIT=iunelp)
           END IF
           WRITE(*,*)' elements for arc 2', elp
        ELSE
           WRITE(*,*)' initial conditions not available for ',astnap
        ENDIF
     ENDIF
     IF(iarc.eq.4)THEN
        IF(initwo)THEN
! handle special case of attributable elements
           IF(same_station(obs,mall))THEN
              obscodi=obs(1)%obscod_i
           ELSE
              obscodi=500
           ENDIF
           IF(covtwo)THEN
              CALL coo_cha(el,cooy,el,fail_flag,dee,obscodi)
              depdep(1:6,1:6)=dee
              IF(nd.gt.6)THEN
                 depdep(7:nd,1:nd)=0.d0
                 depdep(1:nd,7:nd)=0.d0
                 DO j=7,nd
                    depdep(j,j)=1.d0
                 ENDDO
              ENDIF
              CALL propagunc(nd,unc,depdep(1:nd,1:nd),unc)
           ELSE
              CALL coo_cha(el,cooy,el,fail_flag,OBSCODE=obscodi)
              unc%succ=.false.
           ENDIF
           IF(dyn%ndp.gt.0)THEN
              CALL write_elems(el,astna0,'ML',dyn,nls,ls,UNC=unc,FILE=dummyfile,UNIT=iunelt)
           ELSE
              CALL write_elems(el,astna0,'ML',UNC=unc,FILE=dummyfile,UNIT=iunelt)
           END IF
           WRITE(*,*)' elements for both arcs', el
        ELSE
           WRITE(*,*)' initial conditions not available for '  &
                &                ,astna0//'='//astnap
        ENDIF
     ENDIF
! =====================================================================
  ELSEIF(ifun.eq.9)THEN
! =====================================================================
! compute attributables
     CALL tee(iun_log,'COMPUTE ATTRIBUTABLES=')
     CALL orb_sel(.false.,iarc)
     IF(iarc.eq.0) GOTO 50
     CALL obs_cop(1,iarc) ! copy observations to obsc, obswc
     IF(iarc.eq.1)THEN
        astnac=astna0
     ELSEIF(iarc.eq.2)THEN
        astnac=astnap
     ELSEIF(iarc.eq.3)THEN
        astnac=astna0
     ENDIF
     CALL attri_comp(mc,obsc,obswc,attrc,error)
     IF(error) GOTO 50
     trou=nint(attrc%tdtobs)
     IF(attrc%sph*radeg.gt.sphx)THEN
        WRITE(*,*)' arc too wide ', attrc%sph*radeg
     ELSE
        CALL wri_attri(0,0,astnac,attrc,trou)
     ENDIF
     IF(iarc.eq.1)THEN
        attr0=attrc
        elc=el0
     ELSEIF(iarc.eq.2)THEN
        attrp=attrc
        elc=elp
     ELSEIF(iarc.eq.3)THEN
        attr=attrc
        elc=el
     ENDIF
! set controls to detect missing steps
     artyp=0
     menunam='dummy'
 59  CONTINUE
     CALL menu(iatt,menunam,7,'What to do with attributable?=',     &
 &                    'Assign r, rdot to form ATT elements=',       &
 &                    'Compute arc type=',                          &
 &                    'Compute Admissible Region=',                 &
 &                    'Systematic ranging (AR and MOV sampling)=',  &
 &                    'Imminent Impacts=',                          &
 &                    'MC Imminent Impacts=',                       &
 &                    'Select MOV orbit with minimum chi=')

     IF(iatt.eq.0) GOTO 50
     IF(iatt.eq.1)THEN
599     WRITE(*,*) 'assign r (au) '
        READ(*,*,ERR=599) r
598     WRITE(*,*) 'assign rdot (au/day) '
        READ(*,*,ERR=598) rdot
        CALL attelements(attrc,r,rdot,elc,uncc)
        IF(iarc.eq.1)THEN
           el0=elc
           ini0=.true.
           cov0=.true.
           unc0=uncc
           WRITE(*,*)' ATT elements for arc 1'
           WRITE(*,*) el0
        ELSEIF(iarc.eq.2)THEN
           elp=elc
           inip=.true.
           covp=.true.
           uncp=uncc
           WRITE(*,*)' ATT elements for arc 2'
           WRITE(*,*) elp
        ELSEIF(iarc.eq.3)THEN
           el=elc
           initwo=.true.
           covtwo=.true.
           unc=uncc
           WRITE(*,*)' ATT elements for identification'
           WRITE(*,*) elc
        ENDIF
     ELSEIF(iatt.eq.2)THEN
        CALL obs_cop(1,iarc)
        artyp=arc_type(obsc,obswc,mc,geoc_chi,acce_chi,chi,nigarc,fail_arty)
        WRITE(*,*)' arc type of arc ', iarc,' is ', artyp
        WRITE(*,177)mc,nigarc,fail_arty,geoc_chi,acce_chi,chi
 177    FORMAT('nobs=',I4,' nights=',i3,' errcode=',i4/         &
&              'geoc_chi=',1p,d9.2,' acce_chi=',d9.2,' chi=',d9.2)
        !WRITE(iun_log,*)' arc type of arc ', iarc,' is ', artyp      ! a lot of iun_log WRITEs are removed
        !WRITE(iun_log,177)mc,nigarc,fail_arty,geoc_chi,acce_chi,chi  ! a lot of iun_log WRITEs are removed
     ELSEIF(iatt.eq.3)THEN
        CALL tee(iun_log,'COMPUTE ADMISSIBLE REGION=')
        file=astnac//'.att'
        CALL rmsp(file,le)
        CALL filopn(iunfla,file(1:le),'unknown')
        CALL admis_reg(astnac,iunfla,attrc%tdtobs,attrc%angles,attrc%obscod,nrootsf, &
       & roots,c,refs,xo,vo,100.d0,E_bound,nrootsd,rootsd,attrc%apm,xogeo,vogeo,elong)
        CALL filclo(iunfla,' ')
        WRITE(*,*) 'roots= ', roots(1:nrootsf)
        !WRITE(iun_log,*) 'roots= ', roots(1:nrootsf)  ! a lot of iun_log WRITEs are removed
! =============================================================
! WARNING: this may be not useful, now all in sample_web_grid
! find the range of values for rhodot
        r2=roots(1)
        r1=1.d-4 ! arbitrary lower boundary
! search for a zero of the derivative between r1,r2
        rder0 = -1.d0
        ind = 0
        IF(nrootsd.ge.1) THEN
           IF(nrootsd.eq.1) THEN
              ! do nothing
           ELSEIF(nrootsd.eq.2) THEN
              ! do nothing
           ELSEIF(nrootsd.eq.3) THEN
              DO iii=1,nrootsd
                 IF((rootsd(iii).gt.r1).and.(rootsd(iii).lt.r2)) THEN
                    ind=iii
                 ENDIF
              ENDDO
              IF(ind.gt.2)THEN
                 IF(rootsd(ind-2).gt.r1) THEN
                    rder0=rootsd(ind-2)
                    write(*,*)'r1,rder0,r2',r1,rder0,r2
                 ENDIF
              ENDIF
              IF(rder0.ge.r2) THEN
                 WRITE(*,*)'ERROR!!: rder0 > r2'
              ENDIF
              !            write(*,*)'r1,rder0,r2',r1,rder0,r2
           ELSE
              WRITE(*,*)'not included case: r1,r2 nrootsd=',nrootsd
           ENDIF
        ENDIF
! use logaritmic scale
        nfunc=2
        CALL sample_bound_ne1(0.d0,attrc%eta**2,r1,rder0,r2,8,5,&
             &ntot,c,frv,rhodot,npox,rdw)
        rdotmin=MINVAL(rhodot(1:ntot))
        rdotmax=MAXVAL(rhodot(1:ntot))
        file=astnac//'.pol'
        CALL rmsp(file,le)
        CALL filopn(iunpol,file(1:le),'unknown')
        WRITE(iunpol,178) roots, c, rdotmin,rdotmax
178     FORMAT(11(1x,d14.6))
        CALL filclo(iunpol,' ')
! END WARNING
     ELSEIF(iatt.eq.4)THEN
!  ndir=50;np=50;sigx=5.d0; hmax=34.5d0 defaults assigned in finopt.f90
        CALL sample_web_grid(astnac,ini0,cov0,elc,uncc,     &
             &       csino0,m,obs(1:m),obsw(1:m),nd,artyp,geoc_chi,acce_chi,chi)
     ELSEIF(iatt.eq.5)THEN
        CALL tee(iun_log,'SEARCH FOR IMMEDIATE IMPACTORS=')
        IF(.NOT.mov_flag)THEN
           WRITE(*,*) ' Compute the MOV first!'
           GOTO 59
        ENDIF
! select final time
        WRITE(*,*)' Initial conditions at ',elc%t,'. Select VI search time span: '
        READ(*,*) dtclo
        WRITE(*,*) ' VI search until final time = ', elc%t+dtclo
! search and write output
        CALL immediate_imp(astnac,elc%t,dtclo,csino0,artyp,geoc_chi,acce_chi,chi,m,obs(1:m))
     ELSEIF(iatt.eq.6)THEN
! MonteCarlo computation of impact probability
        IF(.not.(ini0.and.cov0))THEN
           WRITE(*,*)' Montecarlo not ready without nominal solution'
           GOTO 59
        ENDIF
        CALL statcode(attrc%obscod,obscodn)
        WRITE(*,*)' number of MC sample points? 0=exit'
        READ(*,*) nmc
        IF(nmc.le.0) GOTO 59
! select final time
        WRITE(*,*)' Initial conditions at ',el0%t,' select delta t'
        READ(*,*)dtclo
        tmcla=el0%t+dtclo
        CALL mc_res_att(astnac,nmc,tmcla,m,obs,obsw,el0,obscodn,nd,nimp,isucc)
! compute probability and report
        impp=(nimp*1.d0)/(isucc*1.d0)
        simppro=sqrt(impp*(1.d0-impp)/isucc)
        WRITE(*,778)nimp, isucc,impp,simppro
778     FORMAT('impactors= ',I4,' out of nominal= ',I4,' IP= ',F6.4,' sigma(IP)= ', F6.4)
     ELSEIF(iatt.eq.7)THEN
        IF(.NOT.mov_flag)THEN
           WRITE(*,*) ' Compute the MOV first!'
           GOTO 59
        END IF
        !---------------------------------!
        ! Find MOV orbit with minimum chi !
        !---------------------------------!
        WRITE(*,*) ' Orbit with minimum chi: i = ', i_chimin,', j = ', j_chimin,', chi = ', mov_orbit(i_chimin,j_chimin)%chi_mov
        ini0 = .TRUE.
        cov0 = unc0%succ
        el0  = mov_orbit(i_chimin,j_chimin)%el_mov
        unc0 = mov_orbit(i_chimin,j_chimin)%unc_mov
        ! Set variables for output
        chi_min_orb = .TRUE.
        file_orb = astnac//'.orb_chimin'
        CALL rmsp(file_orb,le)
        dummyfile = file_orb(1:le)
     ELSE
        WRITE(*,*) ' Option ', iatt, ' not operational'
     ENDIF
     GOTO 59
! =====================================================================
  ELSEIF(ifun.eq.10)THEN
! =====================================================================
! show all the status flags
     WRITE(*,180)obs0,obsp,ini0,inip,inide,initwo,                  &
   &         cov0,covp,covtwo
180  FORMAT('   obs0',' obsp ',' ini0 ',' inip ','inide ','initwo',    &
     &       ' cov0 ',' covp ','covtwo'                                 &
     &           /10L6)
! check that there are data
     IF(m.eq.0)THEN
        WRITE(*,*)' no observation available '
        GOTO 50
     ENDIF
! give the epoch times for all the elements
     WRITE(*,181)el0%t,elp%t,el%t,elide%t, &
   & obs(1)%time_tdt,obs(m)%time_tdt,obs(m+1)%time_tdt,obs(mall)%time_tdt
181  FORMAT('   t0   ','   tp   ','   tm   ','  tide  ',            &
     &          '  tin0  ','  tfi0  ','  tinp  ','  tfip  '/            &
     &          8f8.1)
! observational data available
     IF(obs0)THEN
        WRITE(*,*)' no obs =',m,' of asteroid ',astna0
        IF(cov0)THEN
           WRITE(*,*)'      obs used in fit=',iob0, ' RMS =',csino0
        ELSE
           WRITE(*,*)'      obs used in fit=',iob0
        ENDIF
     ENDIF
     IF(obsp)THEN
        WRITE(*,*)' no obs =',mp,' of asteroid ',astnap
        IF(covp)THEN
           WRITE(*,*)'      obs used in fit=',iobp, ' RMS =',csinop
        ELSE
           WRITE(*,*)'      observations used in fit=',iobp
        ENDIF
     ENDIF
     IF(obs0.and.obsp)THEN
        WRITE(*,*)' no obs =',mall,' of both asteroids '
        IF(covtwo)THEN
           WRITE(*,*)'      obs used in fit=',iobtwo, ' RMS =',csinor
        ELSE
           WRITE(*,*)'      observations used in fit=',iobtwo
        ENDIF
     ENDIF
     IF(ini0)WRITE(*,*)' coord. orbit of ',astna0,' are ',el0%coo
     IF(inip)WRITE(*,*)' coord. orbit of ',astnap,' are ',elp%coo
     IF(initwo)WRITE(*,*)' coord. orbit of ',astna0//'='//astnap,' are ',el%coo
! =====================================================================
  ELSEIF(ifun.eq.11)THEN
! =====================================================================
! calendar to MJD conversion
     WRITE(*,*)'calendar date, year, month, day, hour?'
     READ(*,*)iy,imo,iday,ihr
     sec=0.d0
     imin=0
     CALL julian(iy,imo,iday,ihr,imin,sec,jd)
     WRITE(*,777)jd-2400000.5d0
777  FORMAT('MJD=',f13.6)
! =====================================================================
  ELSEIF(ifun.eq.12)THEN
! =====================================================================
! Check first and possibly second derivatives at the starting points
     CALL twotes(1,obs(1)%time_tdt,obs(1)%obscod_i,el0,obs(1)%obspos,obs(1)%obsvel,nd)
! =====================================================================
  ELSE
! =====================================================================
! non existing option
     WRITE(*,*)' This I cannot do'
  ENDIF
  IF(init)THEN
     init=.false.
     init2=.true.
  ELSEIF(init2)THEN
     init2=.false.
  ENDIF
  GOTO 50
END PROGRAM fitobs

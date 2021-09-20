! Copyright (C) 1998 by OrbFit Consortium                               
! Version: December 14, 1998 Mario Carpino                              
! ================================================================      
! RMODEL                                                                
! ================================================================      
! This subroutine reads from file 'propag.def' all the propagator option
!   WARNING: we advise the user against changing the file propag.def;    
!          do it at your risk !....                               
! Input options for phisical model! version 1.2 of propagator  
!   MOON/MERCURY/PLUTO         
!       ilun=1  ! 0=no moon    1=yes                                      
!       imerc=1 ! 0=no mercury 1=yes (recommended =1 for asteroids)      
!       iplut=0 ! 0=no pluto   1=yes (recommended =0 for asteroids, =1 for 
!                                   transneptunian)                     
!   RELATIVITY
!       irel=0  ! 0=newtonian 1=gen. relativity (recommended =0 for main 
!                                                =1 for Earth crossing)      
!               ! 2=eliocentric acceleration of the asteroid on the (1/c^2)
!               !   (= PN) level
!   MAJOR PERTURBERS (ASTEROIDS)
!       iast=0      !  0=no asteroids with mass n=no. of massive asteroids    
!       filbe='CPV' ! name of the asteroid ephemerides file 
!                   ! CPV = Ceres, Pallas, and Vesta
!                   ! AST17 = 16 massive asteroids+Pluto
!   CLOSE APPROACHES/ABERRATION/TOPOCENTRIC CORRECTION
!       iclap=1 ! 0=no close approach control 1=yes (recommended =1) 
!               - iclap=1 means for Radau alg. to calculate time pos. and vel.  
!                 at closest approach; for multistep alg. the propagator 
!                 detects the time of close-appr.                        
!               - iclap =0 the subroutine force does not check for possible close  
!                 approaches to major planets andor massive asteroids.  
!       iaber=1 ! aberration 0=no 1=yes (recommended =1 always)          
!       istat=1 ! 0=no topocentric corr 1= yes (recommended =1 always)   
!                 iclap.h90                                                               
!   CONTROLS FOR MINIMUM DISTANCE
!   The minimal distance to define a close-approch is controlled by:    
!       dmea (for the Earth only)
!       dter (for Mercury, Venus)
!       dmoon (for the Moon)
!       dmjup (for giant planets)         
!       dmast (for massive asteroids)
!   YARKOVSKY PARAMETERS   
!   In module yark_pert   
!       iyark=1  ! Yarkovsky force                                       
!       iyarpt=0 ! Partials of Yarkovsky are generally not computed      
!                ! for rhs=2
! ---------------------------------------------------------------
! RHS.EQ.2 : DEBRIS
! options applicable only to SATELLITE case
!       ites=2      ! max harmonic degree
!       irad=0      ! radiation pressure 1=spher.sat
!       itide=0     ! tidal perturbation 1=k2 no lag
!       ipla=2      ! 0=2-body 2=Sun+Moon
! ================================================================      
SUBROUTINE rmodel(rhs0)
  USE dyn_param
  USE least_squares, ONLY : difini, rejini, output_des, ab_mag
  USE propag_state,  ONLY: inipro ! calls inipro
  USE tp_trace,      ONLY: tpplane,scatter_stretch
  USE close_app ! sets former closapl options
  USE yark_pert
  USE perturbations
  USE planet_masses 
  USE force_model
  USE force_sat
  USE spher_harm
  USE astrometric_observations, ONLY: radius ! to set default value
  USE fund_const
  USE output_control
  USE reference_systems
  USE multiple_sol, ONLY: lin_lov, scatterplane,no_outcov
  IMPLICIT NONE 
! ========BEGIN INTERFACE======================                         
  INTEGER,INTENT(IN) :: rhs0
! ========END INTERFACE========================
! close approach options in force_model.mod
! radar quality flag in force_model.mod
! default asteroid radius in astrometric_observations.mod
  CHARACTER*80     :: filbe                                               
  LOGICAL          :: fail,fail1,found 
  INTEGER          :: ll,ia 
! controls for bizarre orbits                                           
  DOUBLE PRECISION :: ecclim, samin,samax,phmin,ahmax,qmax 
  CHARACTER*60 :: comment 
! initialization of rhs
  rhs=rhs0
!****************                                                       
! static memory not required (used only once)                         
!****************                                                       
! read time-range (initialize JPL ephemerides) 
  CALL trange 
! setup the rotation from equatorial to ecliptic and back
  CALL rotpn(roteqec,'EQUM','J2000',0.d0,'ECLM','J2000',0.d0) 
  roteceq=TRANSPOSE(roteqec) 
! restart flag set to default value                                     
  CALL set_restart(.TRUE.) 
! flag for AB magnitudes
  CALL rdnlog('propag.','ab_mag',ab_mag,.true.,found,fail1,fail)
!
! Options for dynamical parameters
  dyn=undefined_dyn_par ! default non nongrav
  nls=0 ! default no dynamical parameters to solve
  ls=0
  IF(rhs.EQ.1)THEN
! Force model flags for asteroids OrbFit
     fail=.FALSE. 
     CALL rdnint('propag.','ilun',ilun,.TRUE.,found,fail1,fail) 
     CALL rdnint('propag.','imerc',imerc,.TRUE.,found,fail1,fail) 
     CALL rdnint('propag.','iplut',iplut,.TRUE.,found,fail1,fail) 
     CALL rdnint('propag.','irel',irel,.TRUE.,found,fail1,fail) 
! In case selmet is not called                                          
     icrel=irel 
     comment='number of perturbers '
     CALL input_int_opt('propag.','iast',iast,.TRUE.,found,comment,iun_log)
     comment='binary ephemerides file '
     CALL input_cha_opt('propag.','filbe',filbe,.TRUE.,found,comment,iun_log)
!     CALL rdnint('propag.','iast',iast,.TRUE.,found,fail1,fail) 
!     CALL rdncha('propag.','filbe',filbe,.TRUE.,found,fail1,fail) 
     CALL rmsp(filbe,ll) 
! Name of ephemerides file
     filbep=filbe(1:ll)//'.bep' 
     filbec=filbe(1:ll)//'.bai' 
     DO ia=1,iast 
        astid(ia)=ia 
     ENDDO
     iatrue=iast
! Check and controls on iast and filbe
     IF(iast.EQ.0)THEN
        IF(iun_log.ne.0) WRITE(iun_log,*) 'perturbing asteroids used: none'
        WRITE(*,*) 'perturbing asteroids used: none'
     ELSEIF(iast.EQ.3)THEN
        IF(filbe.EQ.'')THEN
           WRITE(*,*) 'ERROR! rmodel: iast=', iast,', while filbe=',filbe(1:ll)
           STOP
        ELSEIF(filbe.EQ.'CPV')THEN
           IF(iun_log.ne.0) WRITE(iun_log,*) 'perturbing asteroids used: Ceres, Pallas, and Vesta'
           WRITE(*,*) 'perturbing asteroids used: Ceres, Pallas, and Vesta'
        ELSEIF(filbe.EQ.'AST17')THEN
           WRITE(*,*) 'WARNING! rmodel: iast=', iast,', while filbe=',filbe(1:ll)
        END IF
     ELSEIF(iast.EQ.17)THEN
        IF(filbe.EQ.'' .OR. filbe.EQ.'CPV')THEN
           WRITE(*,*) 'ERROR! rmodel: iast=', iast,', while filbe=',filbe(1:ll)
           STOP
        ELSEIF(filbe.EQ.'AST17')THEN
           IF(iun_log.ne.0) WRITE(iun_log,*) 'perturbing asteroids used: 16 major perturbing asteroids plus Pluto'
           WRITE(*,*) 'perturbing asteroids used: 16 major perturbing asteroids plus Pluto'
        END IF
     ELSEIF(iast.GT.0.AND.iast.LT.3)THEN
        IF(filbe.EQ.'')THEN
           WRITE(*,*) 'ERROR! rmodel: iast=', iast,', while filbe=',filbe(1:ll)
           STOP
        ELSEIF(filbe.EQ.'CPV'.OR.filbe.EQ.'AST17')THEN
           WRITE(*,*) 'WARNING! rmodel: iast=', iast,', while filbe=',filbe(1:ll)
        END IF
     ELSEIF(iast.GT.3.AND.iast.LT.17)THEN
        IF(filbe.EQ.''.OR.filbe.EQ.'CPV')THEN
           WRITE(*,*) 'ERROR! rmodel: iast=', iast,', while filbe=',filbe(1:ll)
           STOP
        ELSEIF(filbe.EQ.'AST17')THEN
           WRITE(*,*) 'WARNING! rmodel: iast=', iast,', while filbe=',filbe(1:ll)
        END IF
     ELSEIF(iast.GT.17)THEN
        WRITE(*,*) 'ERROR! rmodel: iast=', iast,', while filbe=',filbe(1:ll)
        STOP
     ELSE
        WRITE(*,*) 'ERROR! rmodel: meaningless value for iast, iast=',iast
     END IF
! ================WARNING=============================
! yark_exp (if different from default=2) has to be read from option file
! because it is not written in the orbit file; e.g Bennu (Chesley et al. 2014)
     CALL rdnrea('propag.','yark_exp',yark_exp,.TRUE.,found,fail1,fail)
! ================END WARNING=========================
!
! ====================================================
! NON-GRAVITATIONAL FORCES OPTIONS
! ====================================================
     CALL rdnlog('propag.','ngr_opt',ngr_opt,.FALSE.,found,fail1,fail)
     IF(iun_log.gt.0)WRITE(iun_log,*) 'ngr_opt = ',ngr_opt
     WRITE(*,*) 'ngr_opt = ',ngr_opt
     IF(ngr_opt)THEN
        WRITE(*,*) '   rmodel: overriding the default for ngr_opt, now ', ngr_opt
     END IF
     IF(ngr_opt) THEN
! SEP violation
        CALL rdnlog('propag.','sep_viol',sep_viol,.TRUE.,found,fail1,fail)
! dyn%nmod=3 (SEP violation)
        IF(sep_viol)THEN
           dyn%nmod=3
           CALL rdnrea('propag.','eta_sep',eta_sep,.TRUE.,found,fail1,fail)
           dyn%ndp=1
           dyn%dp(1)=eta_sep
           dyn%active(1)=.TRUE.
        ENDIF
! ====================================================
! Solar direct radiation pressure and Yarkovsky effect
! ====================================================
! Direct radiation pressure
        comment='direct radiation pressure in the dynamical model'
        CALL input_int_opt('propag.','ipa2m',ipa2m,.TRUE.,found,comment,iun_log)
        IF(ipa2m.GT.0)THEN
           IF(dyn%nmod.EQ.3)THEN
              WRITE(*,*) '   rmodel: contradictory selection: sep_viol incompatible with ipa2m'
              STOP
           ENDIF
           dyn%nmod=1
        ENDIF
        comment='direct radiation pressure coefficient'
        CALL input_rea_opt('propag.','drpa2m',drpa2m,.TRUE.,found,comment,iun_log)
! Yarkovsky effect
        comment='Yarkovsky effect in the dynamical model'
        CALL input_int_opt('propag.','iyark',iyark,.TRUE.,found,comment,iun_log)
        IF(iyark.GT.0 .AND. dyn%nmod.EQ.3)THEN
           WRITE(*,*) '   rmodel: contradictory selection: sep_viol incompatible with iyark=',iyark
           STOP
        END IF
        IF(iyark.EQ.1.or.iyark.EQ.2)THEN ! from model, no parameters 
           CALL rdnint('propag.','iyarpt',iyarpt,.TRUE.,found,fail1,fail) 
           CALL rdncha('propag.','yardir',yardir,.TRUE.,found,fail1,fail) 
           yarfil=.FALSE. 
           yarini=.FALSE. 
        ELSEIF(iyark.EQ.3)THEN
           dyn%nmod=1
        END IF
        comment='A2 Yarkovsky parameter'
        CALL input_rea_opt('propag.','A2',A2,.TRUE.,found,comment,iun_log)
! Setup of dynamical model nmod=1
        IF(iyark.EQ.3.OR.ipa2m.GT.0)THEN
           dyn%ndp=2
           dyn%dp(1)=drpa2m
           dyn%dp(2)=A2
           IF(iyark.EQ.3) dyn%active(2)=.TRUE.
           IF(ipa2m.GT.0) dyn%active(1)=.TRUE.
        END IF
! Determination of dynamical parameters in model nmod=1
        comment='solve-for parameters selector (drpa2m, A2, or both)'
        CALL input_int_opt('propag.','det_drp',det_drp,.TRUE.,found,comment,iun_log)
        IF(ipa2m.EQ.0)THEN
           IF(det_drp.EQ.1 .OR. det_drp.EQ.3)THEN
              WRITE(*,*) '   rmodel: contradictory selection'
              WRITE(*,*) '   ERROR: det_drp=',det_drp,', but ipa2m=',ipa2m
              STOP
           END IF
        END IF
        IF(iyark.EQ.0)THEN
           IF(det_drp.EQ.2 .OR. det_drp.EQ.3)THEN
              WRITE(*,*) '   rmodel: contradictory selection'
              WRITE(*,*) '   ERROR: det_drp=',det_drp,', but iyark=',iyark
              STOP
           END IF
        END IF
! Setup non-grav parameter determination for asteroids (yarko or drpa2m)
        IF(det_drp.EQ.1)THEN
           nls=1
           ls(1)=1
        ELSEIF(det_drp.EQ.2)THEN
           icrel=2
           nls=1
           ls(1)=2
        ELSEIF(det_drp.EQ.3)THEN
           icrel=2
           nls=2
           ls(1)=1
           ls(2)=2
        ELSEIF(det_drp.NE.0)THEN
           WRITE(*,*) '   rmodel: choice not possible, det_drp ',det_drp
           STOP
        ENDIF
! ====================================================
! Outgassing model, cometary style
! ====================================================
! Selection of the model and input of the parameters
        comment='outgassing model'
        CALL input_int_opt('propag.','ioutgas',ioutgas,.TRUE.,found,comment,iun_log)
        comment='A1 outgassing parameter'
        CALL input_rea_opt('propag.','a1ng',a1ng,.TRUE.,found,comment,iun_log)
        comment='A2 outgassing parameter'
        CALL input_rea_opt('propag.','a2ng',a2ng,.TRUE.,found,comment,iun_log)
        comment='A3 outgassing parameter'
        CALL input_rea_opt('propag.','a3ng',a3ng,.TRUE.,found,comment,iun_log)
        comment='delay outgassing parameter'
        CALL input_rea_opt('propag.','dtdelay',dtdelay,.TRUE.,found,comment,iun_log)
! Compatibility check
        IF(ioutgas.GT.0)THEN
           IF(dyn%nmod.EQ.3)THEN
              WRITE(*,*) '   rmodel: contradictory selection: sep_viol incompatible with ioutgas'
              STOP
           ELSEIF(dyn%nmod.EQ.1.OR.iyark.GT.0)THEN
              WRITE(*,*) '   rmodel: contradictory selection: ipa2m/iyark incompatible with ioutgas'
              STOP
           ENDIF
           dyn%nmod=2
! Setup of outgassing dynamical model nmod=2
           dyn%ndp=4
           dyn%dp(1)=a1ng
           dyn%dp(2)=a2ng
           dyn%dp(3)=a3ng
           dyn%active(1:3)=.TRUE.
           IF(ioutgas.GT.1)THEN
              dyn%dp(4)=dtdelay
              dyn%active(4)=.TRUE.
           END IF
        ENDIF
! Determination of outgassing parameters in nmod=2
        comment='outgassing parameter to solve-for'
        CALL input_int_opt('propag.','det_outgas',det_outgas,.TRUE.,found,comment,iun_log)
        IF(ioutgas.EQ.0 .AND. det_outgas.GT.0)THEN    
           WRITE(*,*) '   rmodel: contradictory selection'
           WRITE(*,*) '   ERROR: det_outgas=',det_outgas,', but ioutgas=',ioutgas
           STOP
        END IF
        IF(ioutgas.EQ.1 .AND. det_outgas.EQ.3)THEN
           WRITE(*,*) '   rmodel: contradictory selection'
           WRITE(*,*) '   ERROR: det_outgas=',det_outgas,', but ioutgas=',ioutgas
           STOP
        END IF
! Setup non-grav parameter determination for outgassing models
        IF(det_outgas.EQ.1)THEN
           nls=2
           ls(1)=1
           ls(2)=2
        ELSE IF(det_outgas.EQ.2)THEN
           nls=3
           ls(1)=1
           ls(2)=2
           ls(3)=3
        ELSE IF(det_outgas.EQ.3)THEN
           IF(dyn%dp(1).EQ.0 .AND. dyn%dp(2).EQ.0 .AND. dyn%dp(3).EQ.0)THEN
              WRITE(*,*) '   rmodel: invalid selection. If det_outgas=3, at least one among a1ng, a2ng and a3ng must be non zero'
              STOP
           ELSE
              nls=4
              ls(1)=1
              ls(2)=2
              ls(3)=3
              ls(4)=4
           END IF
        ELSE IF(det_outgas.NE.0)THEN
           WRITE(*,*) '   rmodel: wrong value of det_outgas ', det_outgas
           STOP
        END IF
! ====================================================     
! Logical control, setup of dyn_param
        IF((iyark.GT.0 .OR. ipa2m.GT.0 .OR. ioutgas.GT.0).AND.sep_viol)THEN
           WRITE(*,*) '   rmodel: incompatible non-grav models and SEP violation'
           WRITE(*,*) 'iyark=',iyark,', ipa2m=',ipa2m,', ioutgas=',ioutgas,' and sep_viol=',sep_viol
           STOP
        ELSEIF(iyark.GT.0.AND.ioutgas.GT.0)THEN
           WRITE(*,*) '   rmodel: incompatible non-grav models, iyark=', iyark,', ioutgas=', ioutgas
           STOP
        END IF
        IF(.NOT.check_det(dyn,ls,nls))THEN
           WRITE(*,*) '   rmodel: check_det has failed. Check ls=',ls,' nls=',nls
           STOP
        END IF
     END IF
  ELSEIF(rhs.EQ.2)THEN
     STOP '***** rmodel: rhs=2 not used in this version *****'
  ELSEIF(rhs.EQ.3)THEN
     WRITE(*,*) '***** rmodel: for rhs=3 not used *****'
     STOP
  ENDIF
! Close approach control                                                
  CALL rdnint('propag.','iclap',iclap,.TRUE.,found,fail1,fail) 
  fix_mole=.FALSE.
  kill_propag=.FALSE.
  tpplane=.FALSE. ! use MTP unless otherwise stated
  min_dist=.TRUE.
  surf_stop=.FALSE.
  eprdot=1.d-10
! Linear LOV and scatter plane 
  lin_lov=.false.
  scatterplane=.false.
  scatter_stretch=.false.
  no_outcov=.false.
!                                                                       
  IF(iclap.NE.0)THEN 
! close approach control if requeststed                                 
     CALL rdnint('propag.','npoint',npoint,.TRUE.,found,fail1,fail) 
     CALL rdnrea('propag.','dmea',dmea,.TRUE.,found,fail1,fail) 
     CALL rdnrea('propag.','dmoon',dmoon,.TRUE.,found,fail1,fail) 
     CALL rdnrea('propag.','dmjup',dmjup,.TRUE.,found,fail1,fail) 
     CALL rdnrea('propag.','dmast',dmast,.TRUE.,found,fail1,fail) 
     CALL rdnrea('propag.','dter',dter,.TRUE.,found,fail1,fail) 
  ENDIF
! optical observations
  CALL rdnint('propag.','iaber',iaber,.TRUE.,found,fail1,fail) 
  CALL rdnint('propag.','istat',istat,.TRUE.,found,fail1,fail) 
! radar flag defaults to false          
  radar=.FALSE. 
  radius=-1.d0
! numerical integrator options                                          
  CALL inipro 

! Options for difcor (including outlier rejection)                      
  CALL difini 
  CALL rejini
! output new/old format for rwo files (only applies in fdiff_cor)
  CALL rdnlog('propag.','output_des',output_des,.TRUE.,found,fail1,fail)
! Options for stopping difcor at bizarre orbits; use default            
  IF(rhs.eq.1)THEN
     ecclim=0.d0 
     samin=0.d0 
     samax=0.d0 
     phmin=0.0d0 
     ahmax=0.d0 
     qmax=0.d0
     CALL bizset(ecclim,samin,samax,phmin,ahmax,qmax)
  ELSEIF(rhs.eq.2)THEN
     WRITE(iun_log,*)' rmodel: bizarre controls set by fitdeb.def'
  ELSE
     WRITE(*,*) 'rmodel: not supported rhs=', rhs
     STOP
  ENDIF 
! availaility of covarinace matrix is false at start                   
  CALL cov_not_av 
! verbosity is set at the minimum level by default                      
  verb_clo=1 
  verb_pro=1 
  verb_obs=1
  verb_dif=1 
  verb_mul=1 
  verb_rej=1 
  verb_io=1
  verb_moid=1
  verb_prelim=1
  verb_covariance=1
  verb_matrix=1
! ====================================================
  IF(fail)STOP '**** rmodel: abnormal end ****' 
  RETURN 
END SUBROUTINE rmodel

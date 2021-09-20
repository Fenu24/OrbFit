! ==============MODULE PRED_OBS======================
! replaces both predict_obs.f90 and obs_compute.f90 
! (but also alph_del.f90 is required, may be added here)
! headers npoint and phase are eliminated 
! OUT OF MODULE
!                  aber1          aberration 
!                  aber2          aberration second order
!                  outobc         observation output
!                  rise_and_set   times of rise and set 
!                  rise_set_iter  iteration of rise_and_set
! ROUTINES   
!                  predic_obs  replaces preobX
!                  predic_obs2 for attributables
!                  pre_obs_att for attribute only
!                  alph_del2    optical obs. with proper motion   public       
!                    oss_dif2       "   
!                  alph-del     optical obs.                      public
!                    oss_dif
!                  r_rdot      radar observations                 public       
!                    deltau        "                                    
!                    deldop1       "                                    
!                    deldop2       "   
! HEADERS and MODULES
!     pred_obs.o: \
!	../suit/ASTROMETRIC_OBSERVATIONS.mod \
!	../suit/FUND_CONST.mod \
!	../suit/OUTPUT_CONTROL.mod \
!	../suit/REFERENCE_SYSTEMS.mod \
!	../suit/STATION_COORDINATES.mod \
!	force_model.o \
!	propag_state.o 

MODULE pred_obs
USE fund_const
USE output_control
USE dyn_param
IMPLICIT NONE
PRIVATE

! PUBLIC ROUTINES
PUBLIC predic_obs, alph_del2, r_rdot, alph_del, pre_obs_att, predic_obs2, outobc

! PUBLIC DATA

! ===============================================================
! former phase header 
! galactic pole alpha=12h 51.3m, delta=27d 7.0min
double precision algal,degal,pangGalCen,gax,gay,gaz
parameter (algal=3.36603348393342d0)
parameter (degal=4.73478785724624d-1)
parameter (pangGalCen= 2.14556815606167d0)
parameter (gax=-0.86787505968543d0)
parameter (gay=-0.19757464383054d0)
parameter (gaz=0.45580384036475d0)
! from "The Proper Motion of Sagittarius A*", Reid, M. J., and Brunthaler, 
! A., The Astrophysical Journal, v616, pg. 883
! The galactic north pole is at RA = 12:51:26.282, Dec = +27:07:42.01 (J2000.0)
! the galactic center at RA = 17:45:37.224, Dec = -28:56:10.23 (J2000.0). 
! The Position angle of the Galactic Centre is 122.932 deg. 

! ===============================================================
! former npoint header
! common to handle data on multiple observations
integer,parameter :: npoinx=4000
DOUBLE PRECISION :: al_m(npoinx),de_m(npoinx),hmag_m(npoinx) ! observation
DOUBLE PRECISION :: el_m(6,npoinx),dynm(ndyx,npoinx) ! line of elements, of dyn.par 
PUBLIC npoinx,al_m,de_m,hmag_m,el_m
! phase, distance to Earth, distance to Sun (to compute magnitude)
! elongation, galactic latitude, apparent motion
double precision phav(npoinx),disv(npoinx),dsunv(npoinx)
double precision elov(npoinx)
double precision gallav(npoinx),adotv(npoinx),ddotv(npoinx)
! but the above need to be public???
PUBLIC disv

CONTAINS
! ===================================================================== 
! PREDIC_OBS   fortran90 version, A. Milani, January 2003
! ===================================================================== 
! method: call with minimum set of arguments results in simple
!         nominal prediction; with optional arguments for linear
!         confidence boundary output in optional args.
!         for nonlinear boundary, the multiple data structure is
!         in the public data of the module, but is not accessed if the 
!         USE statement is restricted to the routine pred_obs  
! problems:
!          1: adot,ddot and other data from former phase header are 
!             part of the output, supplied by alph_del2 or similar
!  input: minimum
!          el    = orbital elements, including
!                     coordinate type EQU, KEP, COM, CAR, ATT
!                     epoch time (MJD,TDT)
!                     orbital elements vector 
!                     absolute magnitude                                   
!                     opposition effect coefficient
!          idsta = station code (integer)                                      
!          tobs  = prediction time (MJD,TDT)                                
!          type  = observation type (O for optical, R=V for radar; S not handled)
!  input optional
!          uncert = orbital uncertainty, including 
!                   covariance at epoch time 
!                    normal matrix at epoch time (should be the inverse)
!          sigma = level of the confidence boundary in RMS values       
!          npo   = number of points in the boundary                     
!          ibv   = Type of depiction                                    
!                       =0 for automatic selection                      
!                       =1 for confidence boundary                      
!                       =2 for line of maximum variation                
!          inl   = handling of nonlinearity                             
!                       =0 for automatic selection                      
!                       =1 for linear ellipse                           
!                       =2 for 2-Body nonlinear propagation of covarianc
!                       =3 for n-body nonlinear propagation of covarianc
!                                                                       
!  output: minimum
!          alpha = right ascension (equatorial J2000), radians          
!          delta = declination (equatorial J2000), radians              
!          hmagn = apparent magnitude, as predicted, from h and g given 
!  output optional :
!          gamad = covariance matrix of observations alpha, delta       
!          sig   = sqrt(eigenvalues) of gamad                           
!          axes  = the eigenvectors of gamad are the columns of this mat
!  hidden output :
!          npo1  = number of output dta points (could be less than npo) 
!          al(npo1),de(npo1) points on the confidence boundary          
!                   (difference with respect to best prediciton, radians
!          elm(npo1) alternate elements for observation time            
!                                                                       
!  In the linear approximation, the ellipse of confidence has principal axes 
!          along axes; the semiaxes lenghts are sig                     
!  In the nonlinear approximation, the boundary is a map of the         
!          confidence ellipse in the elements space                     
!                                                                       
!  WARNING: the magnitudes are often very poorly predictable            
! ============INTERFACE=================================================
SUBROUTINE predic_obs(el,idsta,tobs,type,       &
     &    alpha,delta,hmagn,inl,                                    &
     &    uncert,sigma,npo,ibv,gamad,sig,axes,npo1,               &  
     &    adot0,ddot0,pha0,dis0,dsun0,elo0,gallat0,ecllat0,twobo, &
     &    elev0,azimuth0,elsun0,elmoon0,gallon0,sub_ast_station_light0, umbra0, &
     &    penumbra0,phamoon0)
  USE station_coordinates 
  USE orbit_elements                 
  USE semi_linear
  USE close_app, ONLY: fix_mole,kill_propag
! ============= input ==================================================
  TYPE(orbit_elem), INTENT(IN) :: el ! elements
  DOUBLE PRECISION, INTENT(IN) :: tobs ! elements, epoch time, obs.time
  INTEGER,INTENT(IN) ::  idsta ! station code   
  CHARACTER*1, INTENT(IN) :: type ! observation type
! flag for 2-body approximation; must be .false. for full n-body computa
  LOGICAL,INTENT(IN),OPTIONAL :: twobo 
! ============= output =================================================
  DOUBLE PRECISION, INTENT(OUT) :: alpha,delta,hmagn ! best fit obs., apparent magnitude
! ======optional input ================================
  TYPE(orb_uncert), INTENT(IN),OPTIONAL :: uncert ! covariance, normal matrices
  INTEGER, INTENT(INOUT) :: inl ! nonlinearity; used in output when 0 in input 
  INTEGER, INTENT(IN),OPTIONAL :: npo ! no points
  INTEGER, INTENT(INOUT),OPTIONAL ::  ibv ! conf.bd/LOV; used in output when 0 in input
  DOUBLE PRECISION, INTENT(IN),OPTIONAL :: sigma ! sigma value for the boundary
! ======optional output =================================================
  DOUBLE PRECISION,INTENT(OUT),OPTIONAL :: adot0,ddot0 ! proper motion
  DOUBLE PRECISION,INTENT(OUT),OPTIONAL :: gamad(2,2),axes(2,2),sig(2) ! covariance on sky plane
  DOUBLE PRECISION,INTENT(OUT),OPTIONAL :: pha0,dis0,dsun0,elo0,gallat0,gallon0,ecllat0,elmoon0,phamoon0 ! phase, distance to Earth, distance to Sun 
                                             ! elongation, galactic latitude, Moon phase
  DOUBLE PRECISION,INTENT(OUT),OPTIONAl :: elev0,azimuth0,elsun0 ! elevation and azimuth of obj, elevation of Sun
! points on the confidence boundary (difference w.r. to alpha,delta)    
! WARNING! the output number of points is npo1.le.npo; this beacuse hyperbolic points are discarded 
  INTEGER,INTENT(OUT),OPTIONAL :: npo1 
  LOGICAL,INTENT(OUT),OPTIONAL :: sub_ast_station_light0, umbra0, penumbra0
! ============END INTERFACE=============================================
  DOUBLE PRECISION :: gameq(ndimx,ndimx),ceq(ndimx,ndimx) ! covariance, normal matrices 
! partial derivatives of alpha, delta, w.r. to elements (by columns)    
  DOUBLE PRECISION daddet(ndimx,2),dummy(ndimx) 
! second derivatives of alpha, delta, w.r. to elements (not used)       
!   DOUBLE PRECISION ddade(6,6),dddde(6,6) 
! ===================================================================   
! orthonormal basis, matrix defining the plane of the ellipse           
  DOUBLE PRECISION v(ndimx,ndimx),ceicel(ndimx-2,2) 
! transformation matrix between the two planes                          
  DOUBLE PRECISION b(2,2)
  TYPE(orbit_elem) :: elv
! number of full revolutions around the sky                             
  INTEGER ng,nrev 
! functions                                                             
  DOUBLE PRECISION appmag
! for astronomical unit in km from fund_const
! elongation,distance to Earth, distance to Sun (to compute magnitude)  
  DOUBLE PRECISION adot,ddot,pha,dis,rdot,dsun,elo,gallat,gallon, ecllat, elev,azimuth, elsun, elmoon, phamoon
  LOGICAL sub_ast_station_light, umbra, penumbra
  double precision obs4(4) ! observations: alpha, delta in RAD  alphadot, deltadot in RAD/day
  double precision dobde(4,ndimx) ! partial derivatives of obs, w.r. to asteroid elements and parameters
! for r_rdot
  DOUBLE PRECISION posr(3) ! B-F position of reciever, assumed to be aldo transmitter
  CHARACTER*1 tech ! assumed center of mass correction applied
  CHARACTER(LEN=16)  :: stname 
  INTEGER idstarad  ! with encoded trasnmitter and receiver   
  INTEGER nd ! number of paremeters to be solved, including incond and nongrav 
! ===================================================================   
! constant of gravitation, trigonometric constants from fund_const.mod 
! temporaries, indexes                                                  
  DOUBLE PRECISION dal,ddl,maxsig,minsig,allin,delin 
  INTEGER n,inlu,ibvu,ider,k
  LOGICAL cov ! true if covariance computations are possible
  LOGICAL twobo1 
!****************                                                       
!   static memory not required                                          
!****************                                                       
  IF(PRESENT(twobo))THEN
     twobo1=twobo
  ELSE
     twobo1=.false.
  ENDIF
! ===================================================================== 
  if((type.eq.'R'.or.type.eq.'V').and.(inl.eq.2.or.twobo1))then 
     WRITE(*,*)' predic_obs: mixing of radar and two-body '//      &
          &        'approximation not permitted'                             
     RETURN 
  endif
! simple prediction without derivatives
  cov=PRESENT(uncert).and.PRESENT(sig).and.PRESENT(axes)
  IF(cov)cov=(cov.and.uncert%succ)
  IF(cov)THEN
    gameq=uncert%g
    ceq=uncert%c
    ider=1
  ELSE
    ider=0
  ENDIF
! ===================================================================== 
! decide dimension of vector of parameters to be solved
  nd=6+nls
! compute observation; derivatives (of order 1) if required                
  IF(type.eq.'O')THEN 
     CALL alph_del2 (el,tobs,idsta,obs4,ider,nd,dobde,      &
            &   pha,dis,rdot,dsun,elo,gallat,ecllat,twobo1,elev,azimuth,&
            &   elsun,elmoon,gallon,sub_ast_station_light, umbra,penumbra,phamoon) 
     IF(fix_mole.and.kill_propag) RETURN
     alpha=obs4(1)
     delta=obs4(2)
     IF(ider.gt.0)THEN
        daddet(1:nd,1)=dobde(1,1:nd)
        daddet(1:nd,2)=dobde(2,1:nd)
     ENDIF
! store true apparent motion, etc.  
     IF(PRESENT(adot0))adot0=obs4(3)
     IF(PRESENT(ddot0))ddot0=obs4(4)    
     IF(PRESENT(pha0))pha0=pha 
     IF(PRESENT(dis0))dis0=dis 
     IF(PRESENT(dsun0))dsun0=dsun 
     IF(PRESENT(elo0))elo0=elo 
     IF(PRESENT(gallat0))gallat0=gallat 
     IF(PRESENT(gallon0))gallon0=gallon 
     IF(PRESENT(ecllat0))ecllat0=ecllat
     IF(PRESENT(elev0))elev0=elev
     IF(PRESENT(azimuth0))azimuth0=azimuth
     IF(PRESENT(elsun0))elsun0=elsun
     IF(PRESENT(elmoon0))elmoon0=elmoon
     IF(PRESENT(sub_ast_station_light0))sub_ast_station_light0=sub_ast_station_light
     IF(PRESENT(umbra0)) umbra0=umbra
     IF(PRESENT(penumbra0)) penumbra0=penumbra
     IF(PRESENT(phamoon0)) phamoon0=phamoon
! compute apparent magnitude at time of observation (only if it is not at default value -9.99)                   
     IF(el%h_mag.eq.-9.99d0)THEN
        hmagn=el%h_mag
     ELSE
        hmagn=appmag(el%h_mag,el%g_mag,dsun,dis,pha)
     END IF
  ELSEIF(type.eq.'R'.or.type.eq.'V')THEN 
     tech='c'  ! center of mass correction assumed
     CALL obscoo(idsta,posr,stname) ! transmitter and receiver assumed both idsta
     IF(rhs.eq.2)THEN
        STOP '***** rhs=2 with radar *****'
     ENDIF
     idstarad=idsta*10000+idsta
     CALL r_rdot (el,tobs,idstarad,tech,posr,posr,alpha,delta,       &
     &       nd,daddet(1,1), daddet(1,2),ider) 
     IF(fix_mole.and.kill_propag) RETURN
     hmagn=0.d0 
  ELSE
! output adot...... 
     WRITE(*,*)' predic_obs: this observation type not supported ',type 
     RETURN                                                         
  ENDIF
  IF(.not.cov) RETURN
! ===================================================================== 
! *******************fix for infamous polar bug*********************
  IF(type.eq.'O'.or.type.eq.'S')THEN
     daddet(1:nd,1)=daddet(1:nd,1)*cos(delta)
  ENDIF
! compute ellipse of covariance of alpha*cos(delta),delta
! *******************end fix polar bug******************************
! compute ellipse of covariance of alpha,delta                          
  CALL ellips(daddet(1:nd,1:2),gameq(1:nd,1:nd),sig,axes,gamad,nd)
! confidence boundary?
  IF(.not.PRESENT(ibv).or..not.PRESENT(npo))RETURN
  IF(.not.PRESENT(uncert))THEN
     WRITE(*,*)' predic_obs: error in arguments, c0 missing, inl=',inl
     RETURN
  ENDIF
! If inl=0 then use automatic selection method                          
  IF(inl.eq.0)THEN                                                  
     maxsig=max(sig(1),sig(2))*degrad 
     if(maxsig.le.1.0d0)then                                        
        inlu=1                                                       
! Is it safe to use two body if we may have close approach? NO           
!         elseif(maxsig .le. 5.d0)then                                  
!            inl=2                                                      
     else                                                           
        inlu=3                                                       
     endif
  else
     inlu=inl                                                        
  endif
  inl=inlu ! to output choice done when in auto mode
! If ibv=0 then use automatic selection method                          
  IF(ibv.eq.0)THEN                                                  
     maxsig=max(sig(1),sig(2))                                      
     minsig=min(sig(1),sig(2))                                      
     if(maxsig/minsig.le.200.d0)then                                
        ibvu=1                                                       
     else                                                           
        ibvu=2                                                       
     endif
  else
     ibvu=ibv! left as it was
  endif
  ibv=ibvu ! to output choice done when in auto mode
  if(inlu.eq.2)then                                                  
! 2-body aproximation for the central point, no derivatives   
     CALL alph_del2 (el,tobs,idsta,obs4,ider,nd,dobde,TWOBO=.true.,ELOMOON0=elmoon,PHAMOON0=phamoon) 
     allin=obs4(1)
     delin=obs4(2)                                   
  endif
! ===================================================================== 
! compute ellipse in the elements space                                 
  CALL slinel(daddet(1:nd,1:2),gameq(1:nd,1:nd),ceq(1:nd,1:nd),ceicel(1:nd-2,1:2),b,v(1:nd,1:nd),nd)
! ===========================================================           
! compute line of orbital elements                                      
  CALL linobs(ibvu,npo,el,axes,sig,b,v(1:nd,1:nd),sigma,ceicel(1:nd-2,1:2),el_m,dynm,npo1,nd)
! ===========================================================           
  ng=0                                                              
  DO 7 n=1,npo1   
! chose method to handle nonlinearity                                   
     IF(inlu.eq.1)THEN                                                
! linear map from ellipse
        dal=DOT_PRODUCT(el_m(1:6,n),daddet(1:6,1))
        ddl=DOT_PRODUCT(el_m(1:6,n),daddet(1:6,2))
        al_m(n)=dal          
!       Correction for cos(delta) on sky not required, included in daddet(:,1)
!        al_m(n)=dal*cos(delta)  
        de_m(n)=ddl                                                    
! apparent magnitude is the one of the nominal orbit                    
        hmag_m(n)=hmagn  
! compute elements in plane corresponding to sky plane
        el_m(1:6,n)=el%coord+el_m(1:6,n)                   
     ELSEIF(inlu.eq.2)THEN  
! compute elements in plane corresponding to sky plane
        el_m(1:6,n)=el%coord+el_m(1:6,n) 
        elv=el
        elv%coord=el_m(1:6,n)
! 2-body propagation from ellipse, no derivatives      
        CALL alph_del2 (elv,tobs,idsta,obs4,ider,nd,dobde,&
     &   pha,dis,rdot,dsun,elo,gallat,TWOBO=.true.,ELOMOON0=elmoon,PHAMOON0=phamoon) 
! compute apparent magnitude at time of observation                     
        hmag_m(n)=appmag(elv%h_mag,elv%g_mag,dsun,dis,pha)
! difference is with respect to 2-body approx., used w.r. to true orbit 
!        al_m(n)=obs4(1)-allin                                            
!        Correction for cos(delta) on sky             
        al_m(n)=(obs4(1)-allin)*cos(delin)                                            
        de_m(n)=obs4(2)-delin  
! should we store elongation etc????                                          
     ELSEIF(inlu.eq.3)THEN   
! compute elements in plane corresponding to sky plane
        el_m(1:6,n)=el%coord+el_m(1:6,n)
        elv=el
        elv%coord=el_m(1:6,n)                       
        IF(nd.gt.6)THEN
           DO k=7,nd
              dynm(1:ndyx,n)=dynm(1:ndyx,n)+dyn%dp(1:ndyx)  
           ENDDO       
        ENDIF

! full n-body propagation from ellipse, no derivatives   
        if(type.eq.'O')then  
           CALL alph_del2 (elv,tobs,idsta,obs4,ider,nd,dobde, &
     &   pha,dis,rdot,dsun,elo,gallat,ELOMOON0=elmoon,PHAMOON0=phamoon) 
           al_m(n)=obs4(1)
           de_m(n)=obs4(2)                     
! other prediction data stored in common                                
           phav(n)=pha                                               
           disv(n)=dis                                               
           dsunv(n)=dsun                                             
           elov(n)=elo                                               
           gallav(n)=gallat                                          
           adotv(n)=obs4(3)                                             
           ddotv(n)=obs4(4)                                            
! compute apparent magnitude at time of observation                     
           hmag_m(n)=appmag(elv%h_mag,elv%g_mag,dsun,dis,pha)
        elseif(type.eq.'R'.or.type.eq.'V')then
! tech, posr,idstarad already set           
           CALL r_rdot (elv,tobs,idstarad,tech,posr,posr,al_m(n),de_m(n), & 
                  &             nd,dummy,dummy,0)  
           hmag_m(n)=0.d0                                             
        ELSE                                                         
           WRITE(*,*) 'predic_obs: type ',type,' not handled'  
           RETURN                         
        ENDIF
!        al_m(n)=al_m(n)-alpha                                            
!        de_m(n)=de_m(n)-delta                                            
!       Correction for cos(delta)  -  Projection on sky
        al_m(n)=(al_m(n)-alpha)*cos(delta)                                            
        de_m(n)=de_m(n)-delta                                            
     ELSE                                                            
        WRITE(*,*)' predic_obs: this we have not invented yet ', inl     
        RETURN                                                       
     ENDIF
! keep count of lost revolutions
     IF(type.eq.'O')THEN                                        
        IF(n.eq.1)THEN                                                  
           IF(al_m(n).gt.pig)al_m(n)=al_m(n)-dpig          
        ELSE                                                            
           CALL angupd(al_m(n),al_m(n-1),ng)                                
        ENDIF
     ENDIF
! temporary output
     if(type.eq.'O'.or.type.eq.'S')then         
        write(*,*)n,', RA*cos(delta)/DEC (deg)',al_m(n)*degrad,de_m(n)*degrad,ng    
     elseif(type.eq.'R'.or.type.eq.'V')then        
        write(*,*)n,', R/RDOT (km,km/day)',al_m(n)*aukm,de_m(n)*aukm,ng
     endif
7 ENDDO
! ===================================================================== 
! ensure that LOV is consistent with nominal point                      
! first find midpoint of LOV, assume npo is even                        
  if(ibvu.eq.2)then                                                  
     nrev=nint((al_m(npo/2)+al_m(npo/2+1))/2.d0/dpig)                   
     write(*,*)'debug: nrev:',nrev                                  
     if(nrev.ne.0)then                                              
        do n=1,npo1                                                 
           al_m(n)=al_m(n)-nrev*dpig                                    
        enddo
     endif
  endif
END SUBROUTINE predic_obs
! ======================================================================
! PREDIC_OBS2
! generates attributable with its covariance
! ONLY to be used for ./src/triang directory
! =====================================================================
SUBROUTINE predic_obs2(el,idsta,tobs,att,nd,uncert,rr,pha,dsun,twobo,dobde,elmoon)
  USE station_coordinates 
  USE orbit_elements                 
  USE attributable
  USE close_app, ONLY: fix_mole,kill_propag    
  USE semi_linear
! ============= input ==================================================
  TYPE(orbit_elem), INTENT(IN) :: el ! elements
  DOUBLE PRECISION, INTENT(IN) :: tobs ! obs.time TDT
  INTEGER,INTENT(IN) ::  idsta ! station code   
! ============= output =================================================
  TYPE(attrib), INTENT(OUT) :: att ! best fit obs., apparent magnitude
! observations: alpha, delta in RAD  alphadot, deltadot in RAD/day
! ======optional input ================================
  INTEGER,INTENT(IN):: nd! number of parameters to be solved, 
                                  ! including incond and nongrav 
  TYPE(orb_uncert), INTENT(IN),OPTIONAL :: uncert ! covariance, normal matrices
! flag for 2-body approximation; must be .false. for full n-body computa
  LOGICAL,INTENT(IN),OPTIONAL :: twobo
! ======optional output =================================================
  DOUBLE PRECISION, INTENT(OUT),OPTIONAL :: rr(2) ! predicted r,rdot
  DOUBLE PRECISION, INTENT(OUT),OPTIONAL :: pha,dsun,elmoon !phase, dist. to sun
! partial derivatives of obs, w.r. to asteroid elements (= dadde)
  DOUBLE PRECISION, INTENT(OUT), OPTIONAL:: dobde(4,ndimx)
! ============END INTERFACE=============================================
  DOUBLE PRECISION :: gameq(ndimx,ndimx) ! covariance of elements
  DOUBLE PRECISION :: att4(4), gamad(4,4) ! angles,covariance on att. 4-space
! partial derivatives of alpha, delta, w.r. to elements (by columns, by rows)    
  DOUBLE PRECISION daddet(ndimx,4),dadde(4,ndimx)
  DOUBLE PRECISION axes(4,4),sig(4) ! axes, eigenvalues (NOT to be used:polar bug)
  DOUBLE PRECISION :: pha0,dis0,dsun0,rdot0,elmoon0
! ===================================================================   
! functions                                                             
  DOUBLE PRECISION appmag
! elongation, galactic latitude (not really used)
  DOUBLE PRECISION elo,gallat,ecllat 
! TIME CONVERSION
  integer MJD1,MJD2
  double precision SEC1,SEC2
! ===================================================================   
! temporaries, indexes
  CHARACTER*3 obscod                                                  
  LOGICAL cov ! true if covariance computations are possible
  LOGICAL twobo1,ttwo 
  INTEGER ider,i,j
!****************                                                       
!   static memory not required                                          
!****************   
! simple prediction without derivatives
  cov=PRESENT(uncert)
  IF(cov)cov=(cov.and.uncert%succ)
  IF(cov)THEN
    gameq=uncert%g
    ider=1
  ELSE
    ider=0
  ENDIF
  ttwo=PRESENT(twobo)                              
  IF(ttwo)THEN
     twobo1=twobo
  ELSE
     twobo1=.false.
  ENDIF

! paranoia check
  IF(nd.ne.6+nls)THEN
     write(*,*)'wrong number of parameters! nd,6+nls:',nd,6+nls
     stop
  ENDIF
! ===================================================================== 
! compute observation; derivatives (of order 1) if required                
  CALL alph_del2 (el,tobs,idsta,att4,ider,nd,dadde,      &
            &   pha0,dis0,rdot0,dsun0,elo,gallat,ecllat,twobo1,ELOMOON0=elmoon0)
  IF(fix_mole.and.kill_propag) RETURN
! optional arguments for magnitude
  IF(PRESENT(pha))pha=pha0
  IF(PRESENT(dsun))dsun=dsun0 
  IF(PRESENT(elmoon))elmoon=elmoon0
! also computation of r, rdot
  IF(PRESENT(rr))THEN
     rr(1)=dis0
     rr(2)=rdot0
  ENDIF
  att=undefined_attrib
  DO i=1,4
     att%angles(i)=att4(i)
  ENDDO
  att%tdtobs=tobs
! find UT of observation
  mjd1=FLOOR(att%tdtobs)
  sec1=(att%tdtobs-mjd1)*86400.d0
  CALL cnvtim(mjd1,sec1,'TDT',mjd2,sec2,'UTC')
  att%tutobs=mjd2+sec2/86400.d0
! station code
  CALL codestat(idsta,obscod)
  att%obscod=obscod
! compute apparent magnitude at time of observation                     
  att%apm=appmag(el%h_mag,el%g_mag,dsun0,dis0,pha0) 
  IF(.not.cov) RETURN
  daddet=TRANSPOSE(dadde)
  IF(PRESENT(dobde))THEN
     dobde=dadde
  ENDIF
! *******************fix for infamous polar bug*********************
! NOT DONE 
!  daddet(1:6,1)=daddet(1:6,1)*cos(att(2))
!  daddet(1:6,3)=daddet(1:6,3)*cos(att(2))
! *******************end fix polar bug******************************
! compute ellipsoid of covariance of alpha,delta,alphadot,deltadot
  CALL ellipsoid(nd,daddet(1:nd,1:4),gameq(1:nd,1:nd),sig,axes,gamad)
  att%g=gamad
END SUBROUTINE predic_obs2

! ===================================================================== 
! ALPH_DEL2
! ===================================================================== 
! Computation of alpha, delta, alphadot, deltadot  and their derivatives
! ===================================================================== 
!
! Input                                                                 
!    t0: epoch time                                                     
!    tobs: observation time                                             
!    iobscod: station code (integer) 
!    el: orbital elements 
!    ider: flag for derivatives options:                                
!        0 no derivatives                                               
!        1 only first deriv. of alpha,delta w.r. to coord
! Optional input      
!    twobo: logical flag for 2-body approximation; if .true., 2-body    
!        approximation (but Earth position comes from JPL ephem); if .fa
!        all orbit propagations are full n-body                         
! Output      
!   obs4: attributable, including
!           alpha,delta  computed at time tobs                        
!           adot,ddot their time derivatives
!   dobde: their partial derivatives with respect to 
!         the coordinates of el
! Optional output: pha0,dis0,rdot0,dsun0,elo0,gallat0
!
! ==============INTERFACE============================================   
SUBROUTINE alph_del2 (el,tobs,iobscod,obs4,ider,nd,  &
     &   dobde,pha0,dis0,rdot0,dsun0,elo0,gallat0,eclat0,&
     &   twobo,elev0,azimuth0,elsun0,elomoon0,gallon0,sub_ast_station_light0, &
     &   umbra0,penumbra0,phamoon0) 
  USE propag_state
  USE orbit_elements                 
  USE close_app, ONLY: fix_mole,kill_propag                
! ==============INPUT==========================
  TYPE(orbit_elem), INTENT(IN) :: el ! elements
  double precision, intent(in) :: tobs ! observation time (MJD, TDT)
  integer,intent(in) :: iobscod ! observatory code
  integer, intent(IN) :: ider ! flag to control computation of derivatives
! ======optional input =================================================
  logical,intent(in), optional :: twobo ! two-body approx: if .false, full n-body  
! ============OUTPUT==============================                      
  DOUBLE PRECISION obs4(4) ! observations: alpha, delta in RAD  
! alphadot, deltadot in RAD/day
  INTEGER, INTENT(IN) :: nd ! number of partial derivatives
  double precision, intent(OUT), OPTIONAL :: dobde(4,nd) ! partial derivatives 
! of obs, w.r. to asteroid elements
! ======optional output =================================================
! phase, distance to Earth, distance to Sun, solar elongation, galactic and ecliptic latitude
  DOUBLE PRECISION,INTENT(OUT),OPTIONAL :: pha0,dis0,dsun0,rdot0,elo0,gallat0,gallon0,eclat0,elev0,azimuth0,elsun0,elomoon0,phamoon0
  LOGICAL, INTENT(OUT),OPTIONAL :: sub_ast_station_light0, umbra0,penumbra0
! =============END INTERFACE=========================================   
  LOGICAL twobo1 ! control of 2-body approximation
  DOUBLE PRECISION xast(6) ! asteroid cartesian coordinates 
  DOUBLE PRECISION xea(6) ! cartesian coordinates of the Earth 
  DOUBLE PRECISION dobdx(4,6) ! partial derivatives of alpha, delta, 
                        !  adot, ddot  w.r. elements and parameters
! derivatives of cart. coord. w. r. elements 
  DOUBLE PRECISION dxdpar(6,nd)!,ddxde(3,6,6) 
! elongation,distance to Earth, distance to Sun (to compute magnitude)  
  DOUBLE PRECISION pha,dis,rdot,dsun,elo,gallat,eclat,elev,azimuth,elsun, elomoon,gallon,phamoon
  LOGICAL sub_ast_station_light, umbra,penumbra
! **********************************************************************
!****************                                                       
!   static memory not required                                          
!****************             
! Orbit propagation:                                                    
  IF(PRESENT(twobo))THEN
     twobo1=twobo
  ELSE
     twobo1=.false. ! default is full n-body
  ENDIF
! propagation to time t2
  CALL propag(el,tobs,xast,xea,ider,nd,dxdpar,twobo1) 
  IF(fix_mole.and.kill_propag) RETURN
! Computation of observations                                           
  CALL oss_dif2(xast,xea,tobs,iobscod,obs4,ider,dobdx,      &
     &   pha,dis,rdot,dsun,elo,gallat,eclat,elev,azimuth,elsun,elomoon,gallon,sub_ast_station_light, umbra,penumbra,phamoon) 
 ! store true phase, etc.  
  IF(PRESENT(pha0))pha0=pha 
  IF(PRESENT(dis0))dis0=dis 
  IF(PRESENT(dsun0))dsun0=dsun
  IF(PRESENT(rdot0))rdot0=rdot 
  IF(PRESENT(elo0))elo0=elo 
  IF(PRESENT(gallat0))gallat0=gallat 
  IF(PRESENT(gallon0))gallon0=gallon 
  IF(PRESENT(eclat0))eclat0=eclat
  IF(PRESENT(elev0))elev0=elev
  IF(PRESENT(azimuth0))azimuth0=azimuth
  IF(PRESENT(elsun0))elsun0=elsun
  IF(PRESENT(elomoon0))elomoon0=elomoon
  IF(PRESENT(phamoon0))phamoon0=phamoon
  IF(PRESENT(sub_ast_station_light0))sub_ast_station_light0=sub_ast_station_light
  IF(PRESENT(umbra0)) umbra0=umbra
  IF(PRESENT(penumbra0)) penumbra0=penumbra
  IF(ider.lt.1) RETURN 
! derivatives with respect to equinoctal elements                       
  IF(PRESENT(dobde))THEN
     dobde=MATMUL(dobdx,dxdpar(1:6,1:nd))
  ENDIF 
END SUBROUTINE alph_del2
! ===================================================================== 
! OSS_DIF2                                                               
! ===================================================================== 
! Corrections to observations                                           
! ===================================================================== 
! Input                                                                 
! xast asteroid cartesian coordinates at time tauj                      
! xea Earth cartesian coordinates at time tauj                          
!     both in the ecliptic system                                       
! tauj observation time                                                 
! ioc station code                                                      
! ider flag for derivatives options:                                    
!       =0 no derivatives                                               
!       <2 only first deriv. of $\alpha,\delta$ w.r.                    
! to positions                                                          
!       great equal 2 also second partial derivatives                   
! (approximation with 2-body case)                                      
!                                                                       
! Output                                                                
! obs=alpha,delta,adot,ddot computed at time tauj (in the equatorial sys
! (if required)                                                         
! dobdx matrix of first derivatives w. r. to ecliptic positions         
! (if required)                                                         
! ===================================================================== 
SUBROUTINE oss_dif2(xast,xea,tobs,iobscod,obs4,ider,dobdx,      &
     &   pha0,dis0,rdot0,dsun0,elo0,gallat0,eclat0,elev,azimuth,elsun,elomoon,gallon0,sub_ast_station_light,umbra,penumbra,phamoon) 
  USE reference_systems, ONLY: observer_position,obs2ecl 
  USE force_model
  USE fund_const
  USE station_coordinates
  USE planet_masses, ONLY:istat, iaber
! INPUT
  DOUBLE PRECISION, INTENT(IN) :: tobs ! observation time (MJD, TDT)
  DOUBLE PRECISION, INTENT(IN) :: xast(6),xea(6) ! ast and Earth cartesian coord
  INTEGER, INTENT(IN) :: iobscod,ider ! obs.code(integer), control on no. derivatives
! OUTPUT
  DOUBLE PRECISION, INTENT(OUT) :: obs4(4) ! observations alpha, delta, aodt,ddot
  DOUBLE PRECISION, INTENT(OUT) :: dobdx(4,6) ! partials of obs  w.r.t. cartesian ecliptic
! phase, distance to Earth, distance to Sun elongation, galactic latitude, elevation, elev. Sun
  DOUBLE PRECISION,INTENT(OUT) :: pha0,dis0,rdot0,dsun0,elo0,gallat0,gallon0,eclat0, elev,elsun,elomoon,azimuth 
  DOUBLE PRECISION,INTENT(OUT),OPTIONAL :: phamoon
  LOGICAL, INTENT(OUT),OPTIONAL :: sub_ast_station_light, umbra, penumbra
! ===================================================================== 
  double precision d(6),dc(3) ! vector difference of cartesian coordinates, ecliptic 
! WARNING: even if the velocities are not always propagated, they are available
! at the time of the observations                                       
!     note that alpha, delta are equatorialWRITE(*,*) 'Test-1',elmoon

  DOUBLE PRECISION dcor(6) ! after aberration correction
  DOUBLE PRECISION xo(3),vo(3),xobsea(6) ! topocentric position of the observatory, helyocentric 
  DOUBLE PRECISION xo2(3),vo2(3),xast1(6),xea1(6),d1(6),xo1(3),vo1(3),xobsea1(6),d2(6)
  DOUBLE PRECISION cospha,coselo,cosphamoon ! phase and elongation cosine
  DOUBLE PRECISION xmoon(6),xobmoon(3),coselomoon,sinelomoon(3) ! heliocentric, obs-centric of Moon,cos and sin of Moon elongation
  DOUBLE PRECISION sinelo,vvv(3)
  DOUBLE PRECISION dz,dea,dc_s ! auxiliary var
  DOUBLE PRECISION vsize,prscal ! real functions 
  DOUBLE PRECISION deq(6),tmp(3),ddd(3) ! rotation matr, rotated vectors
  DOUBLE PRECISION alpha,delta,adot,ddot ! scalar variables to preserve old code 
! partials of equatorial alpha, delta w.r. to cartesian ecliptic coord. 
  DOUBLE PRECISION dadx(3),dddx(3),ddadx(3,3),ddddx(3,3) 
! aux. var for computation of second derivatives                        
  DOUBLE PRECISION x,y,z
  DOUBLE PRECISION h_umbra, z_ast, d_cone, h_penumbra, z_ast_penumbra, d_cone_penumbra,par_rad,out_penumbra
  DOUBLE PRECISION den,x2,y2,z2,x4,y4 
  DOUBLE PRECISION x2y2,x2z2,y2z2 
  DOUBLE PRECISION coscoelev, coscoelsun,singallon0,cosgallon0,gast
  DOUBLE PRECISION xolon,xolat,xoalt,xovert(3),xovert2(3),xoequ(3),angH
  CHARACTER(LEN=16)  :: stname 
! ===================================================================== 
! Difference vector in ecliptic/equatorial 
! coordinates gives a geocentric vector
  d=xast-xea
  dobdx=0 ! initialization
! ===================================================================== 
! Displacement of the station with respect to the center of the Earth   
  IF(istat.gt.0.and.iobscod.ne.500)THEN 
     IF(rhs.ne.1.and.rhs.ne.2)THEN
        WRITE(*,*)'oss_dif2: rhs=', rhs
        STOP
     ENDIF
     CALL observer_position(tobs,xo,vo,OBSCODE=iobscod,GAST=gast)
!     call pvobs(tobs,iobscod,xo,vo) 
! topocentric correction in ecliptic(rhs=1)/equatorial(rhs=2) coordinates
     d(1:3)=d(1:3)-xo  
     d(4:6)=d(4:6)-vo
     d2(1:3)=MATMUL(roteceq,d(1:3))*149597870691.d0
     d2(4:6)=MATMUL(roteceq,d(4:6))*149597870691.d0
! Elevation over horizon 
     CALL obscoo(iobscod,xoequ,stname)
     xoequ=xoequ*149597870691.d0
     CALL cartesian_to_geodetic(xoequ,xolon,xolat,xoalt,xovert)
     CALL obs2ecl(tobs,iobscod,xovert,xovert2)
     coscoelev=prscal(xovert2,d)/(vsize(d)*vsize(xovert2))
     IF(coscoelev.ge.-1.d0.and.coscoelev.le.1.d0)THEN
        elev=pig/2.d0-acos(coscoelev)
     ELSEIF(coscoelev.gt.1.d0)THEN
        elev=0.d0
        IF(verb_obs.gt.9)THEN
!          WRITE(ierrou,*) 'oss_dif2: xo,d,dis0, coscoelev',xo,d,dis0,coscoelev
           WRITE(ierrou,*) 'oss_dif2: coscoelev',coscoelev
           numerr=numerr+1
        ENDIF
     ELSE
        elev=-pig/2.d0
        IF(verb_obs.gt.9)THEN
!          WRITE(ierrou,*) 'oss_dif2: xo,d,dis0, coscoelev',xo,d,dis0,coscoelev
           WRITE(ierrou,*) 'oss_dif2: coscoelev',coscoelev
           numerr=numerr+1
        ENDIF
     ENDIF
  ELSE
     xo=0.d0
     vo=0.d0
     IF (iobscod.eq.500) THEN
        elev=0.d0
     ENDIF
  ENDIF
  xobsea(1:3)=xo+xea(1:3)
  xobsea(4:6)=vo+xea(4:6)
! ===================================================================== 
! Aberration (only time delay)                                    
  IF(iaber.eq.1)THEN 
     CALL aber1(d,xast(4:6),dcor)
     dcor(4:6)=d(4:6)
  ELSEIF(iaber.eq.2)THEN
     CALL aber2(d,xast,xobsea,dcor) 
  ENDIF
  d=dcor
  dis0=vsize(d) 
  dc=d(1:3)+xo(1:3)
  dc_s=vsize(dc)
! ====================================================================
! Computation of Earth Shadow logical and sub-ast location on light logical
! Written by Fabrizio Bernardi 25 Sep 2013
! Equation of cone from (x^2+y^2)/c^2=(z-z0)^2, where c=r/h and z0=h
! r= radius of cone base (Earth Radius,ER), h= height of cone
! h= ER/sin(alpha)=ER*d/(Rsun-ER)
umbra=.FALSE.
sub_ast_station_light=.FALSE.
dea=vsize(xea)*aukm
h_umbra=eradkm*dea/(r_sun-eradkm)
z_ast=(dot_product(dc(1:3),xea(1:3))/vsize(xea))*aukm
IF(z_ast.GT.0.AND.z_ast.LE.h_umbra) THEN
   d_cone=dsqrt(z_ast**2+((eradkm/h_umbra*(z_ast-h_umbra))**2)) ! distance of cone point from Earth center at specific z_ast
   IF(dc_s*aukm.LE.d_cone) umbra=.TRUE. ! object in umbra
ENDIF
IF(z_ast.LE.0) sub_ast_station_light=.TRUE. ! sub-asteroid location on light
! ====================================================================
! Computation of Earth Penumbra logical 
! Written by Fabrizio Bernardi 26 Sep 2013
! Cone equation as for umbra, but h=ER*d/(Rsun+ER)
penumbra=.FALSE.
h_penumbra=eradkm*dea/(r_sun+eradkm)
! max distance considered for penumbra. Out of this the dimming of brightness is meaningless
! depends on par_rad that is the angular ratio between sun and earth as seen from object
par_rad=4.d0
out_penumbra=-par_rad*eradkm*dea/(r_sun-par_rad*eradkm)
z_ast_penumbra=-z_ast
IF(z_ast_penumbra.LT.0.AND.z_ast_penumbra.GE.out_penumbra) THEN ! signs depend on conventions
   ! distance of penumbra cone point from Earth center at specific z_ast
   d_cone_penumbra=dsqrt(z_ast_penumbra**2+((eradkm/h_penumbra*(z_ast_penumbra-h_penumbra))**2)) ! distance of cone point from Earth center at specific z_ast
   IF(dc_s*aukm.LE.d_cone_penumbra) penumbra=.TRUE. ! object in penumbra
ENDIF
! ===================================================================== 
! Computation of solar distance, earth distance, phase, solar elongation, moon elongation
! heliocentric/ecliptic position of asteroid/debris and earth, topocentric 
! ecliptic position of asteroid/debris and
! geocentric/ecliptic position of observer
  IF(rhs.eq.1)THEN
     xast1=xast ! xast1 is the distance from the sun
     xea1=xea
     d1=d
     xo1=xo
     vo1=vo
  ELSEIF(rhs.eq.2)THEN
     xast1(1:3)=MATMUL(roteqec,xast(1:3))
     xast1(4:6)=MATMUL(roteqec,xast(4:6))
     CALL earcar(tobs,xea1,1)
     xast1=xast1+xea1
     d1(1:3)=MATMUL(roteqec,d(1:3))
     d1(4:6)=MATMUL(roteqec,d(4:6))
     xo1=MATMUL(roteqec,xo)
     vo1=MATMUL(roteqec,vo)
  ELSE
     STOP
  END IF
  dsun0=vsize(xast1) 
  xobsea1(1:3)=xea1(1:3)+xo1
  xobsea1(4:6)=xea1(4:6)+vo1
  rdot0=prscal(d1(1:3),d1(4:6))/dis0
  cospha=prscal(d1,xast1)/(dis0*dsun0) 
  IF(cospha.gt.1.d0)THEN
     IF(cospha.lt.1.d0+10.d0*epsilon(1.d0))THEN
        cospha=1.d0
     ELSE
        write(ierrou,*)'oss_dif2: error! cospha=',cospha 
        numerr=numerr+1
     ENDIF
  ENDIF
  IF(cospha.lt.-1.d0)THEN
     IF(cospha.gt.-1.d0-10.d0*epsilon(1.d0))THEN
        cospha=-1.d0
     ELSE
        write(ierrou,*)'oss_dif2: error! cospha=',cospha 
        numerr=numerr+1
     ENDIF
  ENDIF

  pha0=acos(cospha) 
  coselo=-prscal(d1,xobsea1(1:3))/(dis0*vsize(xobsea1)) 
  CALL prvec(d1,xobsea1(1:3),vvv)
  sinelo=-vvv(3)
  elo0=acos(coselo)
  IF(sinelo.lt.0.d0)elo0=-elo0
  eclat0=asin(d1(3)/dis0)
  CALL mooncar(tobs,xmoon,1)
  xobmoon=xmoon(1:3)-xobsea1(1:3)
  coselomoon=prscal(d1,xobmoon)/(dis0*vsize(xobmoon))
  cosphamoon=prscal(xmoon,xobmoon)/(vsize(xmoon)*vsize(xobmoon))
  IF(cosphamoon.gt.1.d0)THEN
     IF(cosphamoon.lt.1.d0+10.d0*epsilon(1.d0))THEN
        cosphamoon=1.d0
     ELSE
        write(ierrou,*)'oss_dif2: error! cosphamoon=',cosphamoon
        numerr=numerr+1
     ENDIF
  ENDIF
  IF(cosphamoon.lt.-1.d0)THEN
     IF(cosphamoon.gt.-1.d0-10.d0*epsilon(1.d0))THEN
        cosphamoon=-1.d0
     ELSE
        write(ierrou,*)'oss_dif2: error! cosphamoon=',cosphamoon
        numerr=numerr+1
     ENDIF
  ENDIF
  phamoon=acos(cosphamoon)
  CALL prvec(d1,xobmoon,vvv)
  sinelomoon=MATMUL(roteceq,vvv)
  elomoon=acos(coselomoon)
  IF(sinelomoon(3).lt.0.d0)elomoon=-elomoon
! =====================================================================
! illumination angle at the station
!************************* VERSION 4.3
!  coscoelsun=prscal(d1,-xea1(1:3))/(dis0*dsun0)
!*************************************
! VERSION 4.2
  coscoelsun=prscal(xo1,-xobsea1(1:3))/(vsize(xo1)*vsize(xobsea1)) 
  IF(coscoelsun.ge.1.d0)THEN
     elsun=pig/2.d0
  ELSEIF(coscoelsun.le.-1.d0)THEN
     elsun=-pig/2.d0
  ELSE
     elsun=pig/2.d0-acos(coscoelsun)
  ENDIF
! ===================================================================== 
! rotation to the equatorial reference system   
  IF(rhs.eq.1)THEN
     deq(1:3)=MATMUL(roteceq,d(1:3)) 
     deq(4:6)=MATMUL(roteceq,d(4:6)) 
     d=deq                 
  END IF
! ===================================================================== 
! Computation of observation: right ascension (radians)                 
  dz=d(1)**2+d(2)**2 
  IF (dz.le.100.d0*epsilon(1.d0)) THEN
! remove singularity at poles 
     alpha=0.d0 
  ELSE 
     alpha=atan2(d(2),d(1)) 
     IF (alpha.lt.0.d0) THEN 
        alpha=alpha+dpig 
     ENDIF
  ENDIF
! Computation of observation: declination (radians)                     
  delta=asin(d(3)/dis0) 
! Computation of Azimuth
  IF(istat.gt.0.and.iobscod.ne.500)THEN 
     angH= gast + xolon - alpha
     azimuth=atan2((-cos(delta)*sin(angH)),(cos(xolat)*sin(delta)-sin(xolat)*cos(delta)*cos(angH)))
     IF(azimuth.LT.0) THEN 
        azimuth = azimuth + dpig
     ENDIF
  ELSE
     azimuth = 0
  ENDIF


! ===================================================================== 
! galactic latitude. Corrected 13/10/2008 by F.Bernardi  
  gallat0=asin(sin(degal)*sin(delta)+cos(degal)*cos(delta)*cos(alpha-algal))
  singallon0=cos(delta)*sin(alpha-algal)/cos(gallat0)
  IF(singallon0.gt.1.d0) singallon0=1.d0
  IF(singallon0.lt.-1.d0) singallon0=-1.d0
  cosgallon0=(cos(degal)*sin(delta)-sin(degal)*cos(delta)*cos(alpha-algal))/cos(gallat0)
  IF(cosgallon0.gt.1.d0) cosgallon0=1.d0
  IF(cosgallon0.lt.-1.d0) cosgallon0=-1.d0
  if(singallon0.ge.0.and.cosgallon0.ge.0) gallon0=pangGalCen-asin(singallon0)
  if(singallon0.ge.0.and.cosgallon0.lt.0) gallon0=pangGalCen-acos(cosgallon0)
  if(singallon0.lt.0.and.cosgallon0.lt.0) gallon0=pangGalCen+acos(cosgallon0)-2*pig
  if(singallon0.lt.0.and.cosgallon0.ge.0) gallon0=pangGalCen+acos(cosgallon0)-2*pig
  if(gallon0.lt.0) gallon0=gallon0+2*pig
! Old Version
!  gallat0=pig/2d0-acos((d(1)*gax+d(2)*gay+d(3)*gaz)/dis0) 

! ===================================================================== 
! Computation of first derivatives of $\alpha$ and $\delta$ w.r. to posi
! (if required): we derived eq. (2.20)                                  
  dadx(1)=-d(2)/dz 
  dadx(2)=d(1)/dz 
  dadx(3)=0.d0 
  dddx(1)=-d(3)*(d(1)/(sqrt(dz)*dis0**2)) 
  dddx(2)=-d(3)*(d(2)/(sqrt(dz)*dis0**2)) 
  dddx(3)=sqrt(dz)/dis0**2 
! ===================================================================== 
! Apparent motion:                                                      
  adot=prscal(dadx,d(4)) 
  ddot=prscal(dddx,d(4)) 
! store into obs vector                                                 
  obs4(1)=alpha 
  obs4(2)=delta 
  obs4(3)=adot 
  obs4(4)=ddot 
! check if observation partials are required                            
  if(ider.eq.0)RETURN 
! ===================================================================== 
! partials of alpha, delta have already been computed                   
! partials of adot,ddot with respect to velocities are the same         
! rotation to the ecliptic reference system
! dobdx is the derivatives of observations w.r.t. cartesian ecliptic     
  IF(rhs.eq.1)THEN                      
     tmp=MATMUL(roteqec,dadx)
  ELSEIF(rhs.eq.2)THEN
     tmp=dadx
  ELSE
     STOP
  END IF
  dobdx(1,1:3)=tmp 
  dobdx(1,4:6)=0.d0 
  dobdx(3,4:6)=tmp
  IF(rhs.eq.1)THEN
     tmp=MATMUL(roteqec,dddx)
  ELSEIF(rhs.eq.2)THEN
     tmp=dddx
  ELSE
     STOP
  END IF
  dobdx(2,1:3)=tmp 
  dobdx(2,4:6)=0.d0 
  dobdx(4,4:6)=tmp
! ===================================================================== 
! partials of adot,ddot with respect to positions require the           
! second derivatives of alpha, delta with respect to the                
! equatorial reference system                                           
! ===================================================================== 
! Computation of second derivatives of $\alpha$ w.r. to positions       
  ddadx(1,1)=2.d0*d(1)*d(2)/dz**2 
  ddadx(1,2)=(d(2)**2-d(1)**2)/dz**2 
  ddadx(2,1)=ddadx(1,2) 
  ddadx(2,2)=-ddadx(1,1) 
  ddadx(3,1)=0.d0 
  ddadx(3,2)=0.d0 
  ddadx(3,3)=0.d0 
  ddadx(2,3)=0.d0 
  ddadx(1,3)=0.d0 
! Computation of second derivatives of $\delta$ w.r. to positions       
  den=1.d0/(dis0**4*dz*sqrt(dz)) 
  x=d(1) 
  y=d(2) 
  z=d(3) 
!                                                                       
  x2=x*x 
  y2=y*y 
  z2=z*z 
  x4=x2*x2 
  y4=y2*y2 
!                                                                       
  x2y2=x2*y2 
  x2z2=x2*z2 
  y2z2=y2*z2 
!                                                                       
  ddddx(1,1)=z*(2.d0*x4+x2y2-y2z2-y4)*den 
  ddddx(2,2)=z*(2.d0*y4+x2y2-x2z2-x4)*den 
  ddddx(1,2)=x*y*z*(z2+3.d0*x2+3.d0*y2)*den 
  ddddx(2,1)=ddddx(1,2) 
  ddddx(3,3)=-2.d0*z*dz**2*den 
  ddddx(1,3)=x*dz*(z2-x2-y2)*den 
  ddddx(3,1)=ddddx(1,3) 
  ddddx(2,3)=y*dz*(z2-x2-y2)*den 
  ddddx(3,2)=ddddx(2,3) 
! =======================================================               
! chain rule for derivatives of adot                                    
  ddd=MATMUL(ddadx,d(4:6))
! =======================================================               
! rotation to the equatorial reference system                           
  IF(rhs.ne.2)THEN
     dobdx(3,1:3)=MATMUL(roteqec,ddd)
  END IF
! =======================================================               
! chain rule for derivatives of adot                                    
  ddd=MATMUL(ddddx,d(4:6))
! =======================================================               
! rotation to the equatorial reference system                           
  IF(rhs.ne.2)THEN
     dobdx(4,1:3)=MATMUL(roteqec,ddd)
  END IF
END SUBROUTINE oss_dif2
! ===================================================================== 
! PRE_OBS_ATT   fortran90 version, A. Milani, May 2003
! ===================================================================== 
! method: call with minimum set of arguments for attribute only
!  input: minimum
!          el    = orbital elements, including
!                  epoch time (MJD,TDT)                                     
!          tobs  = prediction time (MJD,TDT)                                
!          eq    = equinoctial orbital elements vector at time t0
!  output: minimum
!          alpha = right ascension (equatorial J2000), radians          
!          delta = declination (equatorial J2000), radians
!          optional adot,ddot  proper motion            
! ============INTERFACE=================================================
SUBROUTINE pre_obs_att(el,tobs,alpha,delta,adot,ddot,twobo)
!  USE reference_systems 
  USE propag_state
  USE orbit_elements
  USE close_app, ONLY:fix_mole, kill_propag
! ============= input ==================================================
  TYPE(orbit_elem), INTENT(IN) :: el
  DOUBLE PRECISION, INTENT(IN) :: tobs ! obs.time
!  INTEGER, INTENT(IN) ::  idsta ! station code   
! ============= output =================================================
  DOUBLE PRECISION, INTENT(OUT) :: alpha,delta ! best fit obs.
  DOUBLE PRECISION, INTENT(OUT), OPTIONAL :: adot,ddot ! proper motion
  LOGICAL, INTENT(IN), OPTIONAL :: twobo
! ============END INTERFACE============================================= 
! constant of gravitation, trigonometric constants from fund_const.mod 
  double precision xast(6) ! asteroid cartesian coordinates 
  double precision xea(6) ! cartesian coordinates of the Earth 
  double precision dxdpar(6,ndimx) ! derivatives of cart. coord. w.r. to elements 
  INTEGER ider
  double precision dadx(3),dddx(3) ! partial of angles w.r. cart. coord
  double precision d(6),deq(6) ! diff. of cartesian coord., ecliptic, equatorial 
  double precision dis0,vsize,dz,prscal ! distance from earth, from z axis
  double precision rot(3,3),rotinv(3,3) ! rotation matr,inverse
  LOGICAL twobo1
! ===================================================================== 
  IF(rhs.gt.1)THEN
     STOP '****** pre_obs_att: rhs =2 not allowed *****' 
  ENDIF
! Orbit propagation:                                                    
  IF(PRESENT(twobo))then
     twobo1=twobo
  ELSE
     twobo1=.false. ! default is full n-body
  ENDIF
! compute observation; no derivatives
  ider=0
! full n-body numerical integration                                     
  call propag(el,tobs,xast,xea,ider,6,TWOBO=twobo1) 
  IF(fix_mole.and.kill_propag) RETURN
! Difference vector
  d=xast-xea 
! rotation to the equatorial reference system   
  deq(1:3)=MATMUL(roteceq,d(1:3)) 
  deq(4:6)=MATMUL(roteceq,d(4:6))
  d=deq                         
! Computation of observation: right ascension (radians)                 
  dz=d(1)**2+d(2)**2 
  if (dz.le.100*epsilon(1.d0)) then
! remove singularity at poles 
     alpha=0.d0 
  else 
     alpha=atan2(d(2),d(1)) 
     if (alpha.lt.0.d0) then 
        alpha=alpha+dpig 
     endif
  endif
! Computation of observation: declination (radians) 
  dis0=vsize(d)
  delta=asin(d(3)/dis0) 
  IF(.not.(PRESENT(adot).or.PRESENT(ddot)))RETURN
! Computation of first derivatives of $\alpha$ and $\delta$ w.r. to posi
! (if required): we derived eq. (2.20)                                  
  dadx(1)=-d(2)/dz 
  dadx(2)=d(1)/dz 
  dadx(3)=0.d0 
  dddx(1)=-d(3)*(d(1)/(sqrt(dz)*dis0**2)) 
  dddx(2)=-d(3)*(d(2)/(sqrt(dz)*dis0**2)) 
  dddx(3)=sqrt(dz)/dis0**2 
! Apparent motion:
  IF(PRESENT(adot))THEN                                                      
     adot=prscal(dadx,d(4))
  ENDIF
  IF(PRESENT(ddot))THEN  
     ddot=prscal(dddx,d(4))
  ENDIF  
END SUBROUTINE pre_obs_att
! ===================================================================== 
! ALPH_DEL  vers. 3.0
! A. Milani, November 2002                    
! ===================================================================== 
! Computation of alpha and delta and their derivatives                  
! ===================================================================== 
!
! Input
!    el : orbital elements, including epoch time
!    tauj: observation time
!    ioj: station code; pos, vel position and velocity as computed in input    
!    ider: flag for derivatives options:                                
!        0 no derivatives  1 derivatives  
!    twobo: logical flag for 2-body approximation; if .true., 2-body    
!        approximation (but Earth position comes from JPL ephem); if .fa
!        all orbit propagations are full n-body                         
! Output
!   alj,dej: alpha,delta  computed at time tauj                         
!   dade,ddde:  matrices of first derivatives (if required)             
!
! ==============INTERFACE============================================   
SUBROUTINE alph_del (el,tauj,iocj,pos,vel,ider,twobo,nd,alj,dej,dade,ddde, &
     &       adot0,ddot0,pha0,dis0,dsun0,elo0,gallat0,elomoon0,gallon0)
  USE propag_state
  USE orbit_elements
  USE close_app, ONLY : fix_mole,kill_propag  
! kill_propag e' gia' nelle definizioni.....
!=========================================================
! INPUT
! flag to control computation of derivatives                            
  INTEGER, INTENT(IN) ::  ider
! flag to control two-body approximation: if .false., full n-body       
  LOGICAL, INTENT(IN) :: twobo 
! times: epoch time for asteroid elements, observation time (MJD)       
  DOUBLE PRECISION, INTENT(IN) :: tauj 
! asteroid elements                                          
  TYPE(orbit_elem), INTENT(IN) ::  el 
! observatory code and position, velocity
  INTEGER, INTENT(in) :: iocj 
  DOUBLE PRECISION, INTENT(IN) :: pos(3),vel(3)
! OUTPUT
! observations: alpha (right ascension) delta (eclination), in RAD      
  DOUBLE PRECISION, INTENT(OUT) ::  alj,dej
! partial derivatives of alpha, delta, w.r. to asteroid coordinates   
  INTEGER, INTENT(IN) :: nd   ! number of parameters, including init. cond and dynamical
  DOUBLE PRECISION, INTENT(OUT) :: dade(nd),ddde(nd) 
! ======optional output =================================================
  DOUBLE PRECISION,INTENT(OUT),OPTIONAL :: adot0,ddot0 ! proper motion
  DOUBLE PRECISION,INTENT(OUT),OPTIONAL :: pha0,dis0,dsun0,elo0,gallat0,gallon0,elomoon0 ! phase, 
             ! distance to Earth, distance to Sun, elongation, galactic latitude
! =============END INTERFACE=========================================   
  double precision xea(6),xobs(6) ! cartesian coordinates of the Earth, of the observer   
  double precision xast(6) ! asteroid cartesian coordinates
  double precision dadx(3),dddx(3) ! first derivatives of alpha, delta, w.r. to coord
  DOUBLE PRECISION :: adot,ddot,pha,dis,dsun,elo,gallat,gallon,elomoon ! proper motion, data to compute magnitude
  double precision dxdpar(6,ndimx) ! first derivatives of cartesian coordinates with respect to all parameters
  integer j ! loop variables j=1,6;
  double precision prscal ! double precision functions
! **********************************************************************
!****************                                                       
!   static memory not required                                          
!****************                                                       
! Orbit propagation:                                                    
! full n-body numerical integration                                     
  call propag(el,tauj,xast,xea,ider,nd,dxdpar,twobo)
  xobs(1:3)=xea(1:3)+pos
  xobs(4:6)=xea(4:6)+vel
!  WRITE(11,899)tauj,xast,xobs
!899 FORMAT(f15.9,3(1x,f15.12),3(1x,f15.12),3(1x,f15.12),3(1x,f15.12))
!  IF(kill_propag) STOP '**** alph_del: kill propag *****'
  IF(fix_mole.and.kill_propag) RETURN
! Computation of observations                                           
  CALL oss_dif(xast,xea,tauj,iocj,pos,vel,alj,dej,ider,dadx,dddx,&
     &       adot,ddot,pha,dis,dsun,elo,gallat,elomoon,gallon)
! store true phase, etc.
  IF(PRESENT(adot0))adot0=adot
  IF(PRESENT(ddot0))ddot0=ddot
  IF(PRESENT(pha0))pha0=pha 
  IF(PRESENT(dis0))dis0=dis 
  IF(PRESENT(dsun0))dsun0=dsun 
  IF(PRESENT(elo0))elo0=elo 
  IF(PRESENT(gallat0))gallat0=gallat 
  IF(PRESENT(gallon0))gallon0=gallon 
  IF(PRESENT(elomoon0))elomoon0=elomoon 
  if(ider.lt.1)return                                               
! derivatives with respect to equinoctal elements
  DO  j=1,nd
!     dade(j)=prscal(dadx,dxdpar(1,j)) 
!     ddde(j)=prscal(dddx,dxdpar(1,j))
     dade(j)=prscal(dadx,dxdpar(1:6,j)) 
     ddde(j)=prscal(dddx,dxdpar(1:6,j)) 
  ENDDO
END SUBROUTINE alph_del
! ===================================================================== 
! OSS_DIF vers. 3.0 
! vers. 3.0 A. Milani, November 2002 (fortran 90)
! vers. 1.3 A. Milani, June 14, 1997 (equatorial coordinates for observations)
! ===================================================================== 
! Corrections to observations                                           
! ===================================================================== 
! Input                                                                 
! xast asteroid cartesian coordinates at time tauj                      
! xea Earth cartesian coordinates at time tauj                          
!     both in the ecliptic system                                       
! tauj observation time                                                 
! idst station code, pos, vel geocentric position and velocity
! ider flag for derivatives options:                                    
!       =0 no derivatives                                               
!       =1 first deriv. of $\alpha,\delta$ w.r. to positions 
! Output                                                                
! alj,dej,alpha,delta computed at time tauj (in the equatorial system)  
! (if required)                                                         
! dadx,dddx matrices of first derivatives w. r. to ecliptic positions   
! (if required)                                                         
! ddadx,ddddx matrices of second deriv., 2-b approximation              
! ===================================================================== 
 SUBROUTINE oss_dif(xast,xea,tauj,idst,pos,vel,alj,dej,ider,dadx,dddx, &
     &       adot,ddot,pha,dis,dsun,elo,gallat,elomoon,gallon)
   USE fund_const
   USE force_model
   USE planet_masses, ONLY:istat, iaber
! ===================================================================== 
! INPUT
! control derivatives                        
   integer, intent(in) :: ider
! cartesian coordinates, ecliptic, of asteroid, of earth, 
! WARNING: even if the velocities are not always propagated, they are 
! available at the time of the observations                                       
   DOUBLE PRECISION, INTENT(IN) :: xast(6),xea(6)
! time tdtd of the observation
   DOUBLE PRECISION, INTENT(IN) :: tauj
! observatory code, geocentric position, velocity 
   INTEGER, INTENT(in) :: idst
   DOUBLE PRECISION, INTENT(IN) :: pos(3),vel(3)
! OUTPUT
! observations alpha, delta (equatorial, radians)              
   DOUBLE PRECISION, INTENT(OUT) :: alj,dej
! partials of equatorial alpha, delta w.r. to cartesian ecliptic heliocentric coord. (rhs=1)
!                                                 "     equatorial geocentric cood. (rhs=2) 
   DOUBLE PRECISION, INTENT(OUT) :: dadx(3),dddx(3)
! proper motion, data to compute magnitude
   DOUBLE PRECISION, INTENT(OUT) :: adot,ddot,pha,dis,dsun,elo,gallat,gallon,elomoon
! END INTERFACE
! difference vector
   DOUBLE PRECISION d(6)
   double precision dcor(6) ! after aberration correction
! topocentric position of the observatory ***** temporary for check
   double precision xo(3),vo(3),xobsea(6),xast1(6),xea1(6),d1(6),pos1(3),vel1(3),xobsea1(6)
! phase and solar elongation cosine                                           
   double precision cospha,coselo,sinelo,vvv(3)
! heliocentric, obs-centric of Moon,cos and sin of Moon elongation,sin and cos galact.longitude
   DOUBLE PRECISION xmoon(6),xobmoon(3),coselomoon,sinelomoon(3),singallon,cosgallon 
! auxiliary var.                                          
   double precision dz 
! real functions                                                        
   double precision vsize,prscal 
! integer loop index                                                    
   integer i,ii,ij,j 
! rotation matrices from reference_systems.mod, rotated vectors
   double precision deq(6),tmp(3) 
! ===================================================================== 
! Difference vector 
! difference in ecliptic(rhs=1)/equatorial(rhs=2) coordinates gives a geocentric vector
   d=xast-xea 
! ===================================================================== 
! Displacement of the station with respect to the center of the Earth   
   if(istat.gt.0)then 
      IF(idst.ne.500)THEN 
! check *********************
!           call pvobs(tauj,idst,xo,vo) 
!           write(*,*)'pos -xo ', pos-xo
!           write(*,*)'vel -vo ', vel-vo
! end check *********************
!            call vdiff(d,xo,d) 
!            call vdiff(d(4),vo,d(4))
         IF(rhs.ne.1.and.rhs.ne.2)THEN
            WRITE(*,*)'oss_dif2: rhs=', rhs
            STOP
         ENDIF
         d(1:3)=d(1:3)-pos
         d(4:6)=d(4:6)-vel
      ENDIF
   ENDIF
! ===================================================================== 
! Aberration (only time delay)                                          
   xobsea(1:3)=xea(1:3)+pos
   xobsea(4:6)=xea(4:6)+vel
   IF(iaber.eq.1)THEN 
      CALL aber1(d,xast(4:6),dcor) 
      dcor(4:6)=d(4:6)
   ELSEIF(iaber.eq.2)THEN
      CALL aber2(d,xast,xobsea,dcor) 
   ENDIF
   d=dcor
   dis=vsize(d) 
! ===================================================================== 
! Computation of solar distance, earth distance, phase, elongation  
! heliocentric/ecliptic position of asteroid/debris and earth, topocentric 
! ecliptic position of asteroid/debris and
! geocentric/ecliptic position of observer
   IF(rhs.eq.1)THEN
      xast1=xast ! xast1 is the distance from the sun
      xea1=xea
      d1=d
      pos1=pos
      vel1=vel
   ELSEIF(rhs.eq.2)THEN
      xast1(1:3)=MATMUL(roteqec,xast(1:3))
      xast1(4:6)=MATMUL(roteqec,xast(4:6))
      CALL earcar(tauj,xea1,1)
      xast1=xast1+xea1
      d1(1:3)=MATMUL(roteqec,d(1:3))
      d1(4:6)=MATMUL(roteqec,d(4:6))
      pos1=MATMUL(roteqec,pos)
      vel1=MATMUL(roteqec,vel)
   ELSE
      STOP
   END IF    
   dsun=vsize(xast1) 
   xobsea1(1:3)=xea1(1:3)+pos1
   xobsea1(4:6)=xea1(4:6)+vel1
   cospha=prscal(d1,xast1)/(dis*dsun) 
   pha=acos(cospha) 
   coselo=-prscal(d1,xobsea1)/(dis*vsize(xobsea1)) 
   elo=acos(coselo) 
   CALL prvec(d1,xobsea1,vvv)
   sinelo=-vvv(3)
   If(sinelo.lt.0.d0)elo=-elo
   CALL mooncar(tauj,xmoon,1)
   xobmoon=xmoon(1:3)-xobsea1(1:3)
   coselomoon=prscal(d1,xobmoon)/(dis*vsize(xobmoon))
   CALL prvec(d1,xobmoon,vvv)
   sinelomoon=MATMUL(roteceq,vvv)
   elomoon=acos(coselomoon)
   IF(sinelomoon(3).lt.0.d0)elomoon=-elomoon
! ===================================================================== 
! rotation to the equatorial reference system
   IF(rhs.eq.1)THEN                           
      call prodmv(deq,roteceq,d) 
      call prodmv(deq(4),roteceq,d(4))
! trick to change as little as possible from vers. 1.2 to 1.3           
      d(1:6)=deq(1:6)
   ELSEIF(rhs.eq.3)THEN
      WRITE(*,*)'oss_dif: rhs= ', rhs
      STOP 
   ENDIF
! ===================================================================== 
! Computation of observation: right ascension (radians)                 
   dz=d(1)**2+d(2)**2 
   if (dz.le.100.d0*epsilon(1.d0)) then 
      alj=0.d0 
   else 
      alj=atan2(d(2),d(1)) 
      if (alj.lt.0.d0) then 
         alj=alj+dpig 
      endif
   endif
! Computation of observation: declination (radians)                     
   dej=asin(d(3)/dis) 
! ===================================================================== 
! galactic latitude. Corrected 13/10/2008 by F.Bernardi  
  gallat=asin(sin(degal)*sin(dej)+cos(degal)*cos(dej)*cos(alj-algal))
  singallon=cos(dej)*sin(alj-algal)/cos(gallat)
  cosgallon=(cos(degal)*sin(dej)-sin(degal)*cos(dej)*cos(alj-algal))/cos(gallat)
  if(singallon.ge.0.and.cosgallon.ge.0) gallon=pangGalCen-asin(singallon)
  if(singallon.ge.0.and.cosgallon.lt.0) gallon=pangGalCen-acos(cosgallon)
  if(singallon.lt.0.and.cosgallon.lt.0) gallon=pangGalCen+acos(cosgallon)-2*pig
  if(singallon.lt.0.and.cosgallon.ge.0) gallon=pangGalCen+acos(cosgallon)-2*pig
  if(gallon.lt.0) gallon=gallon+2*pig
! Old Version
!  gallat0=pig/2d0-acos((d(1)*gax+d(2)*gay+d(3)*gaz)/dis0) 
! ===================================================================== 
! Computation of first derivatives of $\alpha$ and $\delta$ w.r. to posi
! (if required): we derived eq. (2.20)                                  
   dadx(1)=-d(2)/dz 
   dadx(2)=d(1)/dz 
   dadx(3)=0.d0 
   dddx(1)=-d(3)*(d(1)/(sqrt(dz)*dis**2)) 
   dddx(2)=-d(3)*(d(2)/(sqrt(dz)*dis**2)) 
   dddx(3)=sqrt(dz)/dis**2
! ===================================================================== 
! Apparent motion:                                                      
   adot=prscal(dadx,d(4)) 
   ddot=prscal(dddx,d(4)) 
   if(ider.eq.0)return 
! ===================================================================== 
! rotation to the equatorial reference system 
   IF(rhs.eq.1)THEN  
      dadx=matmul(roteqec,dadx)
      dddx=matmul(roteqec,dddx) 
   ENDIF
 END SUBROUTINE oss_dif
! ==================================
!
! R _ R D O T
!
! radar observations
SUBROUTINE r_rdot (el,tr,ioc,tech,posr,post,r,v,nd,drde,dvde,ider) 
  USE reference_systems, ONLY: observer_position 
  USE propag_state
  USE astrometric_observations, ONLY: radius ! to set default value
  USE orbit_elements
  USE fund_const
  USE ever_pitkin
! =============INPUT====================                                
  TYPE(orbit_elem),intent(IN) :: el ! asteroid equinoctal elements 
! observation technology (used for surface bounce correction)                 
  CHARACTER*1, INTENT(IN):: tech 
! observation time (MJD)       
  double precision, intent(in) :: tr 
! observatory code (integer) with encoded transmitter and receiver
  integer, intent(in) :: ioc 
! flag to control computation of derivatives                            
  integer, intent(in):: ider 
  DOUBLE PRECISION, INTENT(IN) :: posr(3), post(3) 
! normally position and velocity of observer
! but in this special case, body fixed geocentric position 
! of receiver and of transmitterer
! ============OUTPUT====================                                
! observations: range and range rate in AU, au/day                      
  double precision, intent(out) :: r,v 
! partial derivatives of range and range w.r. to asteroid coordinates
  INTEGER, INTENT(IN) :: nd ! no. parameters including incond nongrav   
  double precision, intent(out) ::  drde(nd),dvde(nd) 
! =============END INTERFACE=========================================   
! cartesian coordinates of the Earth, of the asteroid, id. at receive time; 
! barycentric cood. of Sun
  double precision xea(6),xast(6),xastr(6),xsun(6) 
! first partial derivatives of r, rdot, w.r. to ast. coordinates, vel.  
  double precision drdx(6),dvdx(6) 
! asteroid radius (for surface bounce correction)                       
  double precision rb 
! station codes                                                         
  INTEGER iotr,iore 
! station positions: geocentric, heliocentric                           
  double precision xre(3),yre(3),xtr(3),ytr(3)
  double precision xre1(3),yre1(3),xtr1(3),ytr1(3) 
  double precision rre(3),vre(3),rtr(3),vtr(3) 
! difference vector, distance, velocity difference, size                
  double precision rhorv(3),rhor,rhordv(3),rhord 
  double precision rhotv(3),rhot,rhotdv(3),rhotd 
  double precision vsize,prscal !functions 
! solar system barycentric velocities (for 1/c^2 Doppler level)         
  double precision vressb(3),vastssb(3),vtrssb(3) 
  double precision vressb2,vtrssb2 
! down leg time, bounce time                                            
  double precision taud,taudold,tb,tbold 
! up leg time, transmit time                                            
  double precision tauu,tt,ttold,tauuold 
! delta tau correction function is part of module 
! speed of light from fund-const.mod
! radius of asteroid from astrometric_observations.mod 
! iteration index, control and flag for the position/velocity subroutine 
  INTEGER i,ifla 
  INTEGER, PARAMETER :: itmax=10 
! control on convergence set at 0.05 microseconds (Yeomans et al. AJ 199
  DOUBLE PRECISION ept 
!     PARAMETER (ept=6.d-13)                                            
! correction 9/1/2002: control on taud, not on time in MJD              
  PARAMETER (ept=1.d-16) 
  DOUBLE PRECISION tbdif(itmax),taudif(itmax) 
! time scale                                                            
  double precision tretdb,temp,tdiffr,tdifft 
! 2-body step                                                           
  DOUBLE PRECISION eqast(6),enne 
! first derivatives of cartesian coordinates with respect to elements   
  double precision dxdpar(6,ndimx) 
! auxiliary scalars for Doppler GR & tropospheric corrections           
  double precision rgeo,rgeod,deldoptr,deldopre,rast,rastd 
  double precision rsta,rstad,scal1,scal2,scal3,levelc 
  double precision dtaud,dtauu
! =================================================                     
! deal with surface vs. mass center return:                             
  if(tech.eq.'s')then 
     rb=radius 
     if(rb.le.0)then 
        write(*,*)  '**** rrdot: internal error', radius 
        stop 
     endif
  elseif(tech.eq.'c')then 
! no correction is needed                                               
     rb=0 
  else
     WRITE(*,*)' r_rdot; error in bs%tech', tech,' at time ', tr 
     WRITE(ierrou,*)' r_rdot; error in bs%tech', tech,' at time ', tr 
     numerr=numerr+1
     rb=0        
  endif
! ======================================                                
! Displacement of the stations with respect to the center of the Earth  
! find codes of two observatories                                       
  iotr=ioc/10000 
  iore=ioc-iotr*10000 
  ifla=1 
! compute position of receiver at time tr   
  CALL observer_position(tr,xre,yre,BFPOS=posr,PRECISION=2)
!  call pvobs3(tr,posr,xre,yre)
! =======================================                               
! get receive time in TDB                                               
  call times(tr+2400000.5,temp,tdiffr) 
  tretdb=tr+tdiffr/86400.d0 
! Compute down leg time                                                 
! Orbit propagation at time tr                                          
  call propag(el,tretdb,xastr,xea,ider,nd,dxdpar) 
! ... new initial conditions are xastr                                  
! receive station position at time tr                                   
  rre=xea(1:3) + xre  ! CALL vsumg(3,xea,xre,rre) 
  vre=xea(4:6) + yre  ! CALL vsumg(3,xea(4),yre,vre) 
! initial guess is bounce time(position) equal to receive time(position)
  tb=tretdb 
  taud=0.d0 
  xast=xastr 
 ! Loop on down leg                                                      
  DO i=1,itmax 
     rhorv=xast(1:3)-rre ! CALL vdiff(xast,rre,rhorv) 
     rhor=vsize(rhorv) 
     dtaud=deltau(xast,xre,rhorv,rre) 
! correction                                                            
     taudold=taud 
     taud=(rhor-rb)/vlight + dtaud 
     tbold=tb 
     tb=tretdb-taud 
! convergence control                                                   
     tbdif(i)=tb-tbold 
     taudif(i)=taud-taudold 
!         WRITE(*,*)tb,tb-tbold,taud,taud-taudold                       
     IF(abs(taud-taudold).lt.ept) GOTO 9 
! Orbit propagation at time tb    
!     CALL coocha(xastr,'CAR',gms,eqast,'EQU',enne) 
     CALL fser_propag(xastr(1:3),xastr(4:6),tretdb,tb,gms,xast(1:3),xast(4:6))
!     CALL prop2b(tretdb,eqast,tb,xast,gms,0,dxde,ddxde) 
  ENDDO
! too many iterations                                                   
  WRITE(*,*)' slow conv. on down leg time ',(tbdif(i),i=1,itmax),   &
     &  (taudif(i),i=1,itmax)                                           
! compute relative velocity between asteroid and receiver               
9 CONTINUE 
  rhordv=xast(4:6)-vre  ! CALL vdiff(xast(4),vre,rhordv) 
  rhord=prscal(rhorv,rhordv)/rhor 
! solar velocity at the receive and bounce epochs to get the asteroid   
! barycentric velocity (vastssb) and the receiver barycentric velocity  
! (vressb)                                                              
  ifla=2 
! - receive                                                             
  CALL earcar(tretdb,xsun,ifla) 
  vressb=xsun(4:6)+vre ! CALL vsumg(3,vre,xsun(4),vressb) 
  vressb2=vressb(1)*vressb(1)+vressb(2)*vressb(2)+vressb(3)*vressb(3)
! - bounce                                                              
  CALL earcar(tb,xsun,ifla) 
  vastssb=xast(4:6)+xsun(4:6)  ! CALL vsumg(3,xast(4),xsun(4),vastssb) 
! =======================================                               
! compute upleg time                                                    
! fist guess is up leg time = down leg time                             
  ifla=1 
  tt=tb-taud 
  tauu=taud 
! Loop on upleg time                                                    
  DO i=1,itmax 
! compute transmitter position at estimated transmit time tt            
! call the Earth-rotation model in TDT                                  
     CALL times(tt+2400000.5,temp,tdifft) 
     CALL earcar(tt,xea,ifla) 
     CALL observer_position(tt-(tdifft/86400.d0),xtr,ytr,BFPOS=post,PRECISION=2)
!     CALL pvobs3(tt,post,xtr,ytr)
     rtr=xea(1:3) + xtr   
! upleg time                                                            
     rhotv=xast(1:3)-rtr 
     rhot=vsize(rhotv) 
     dtauu=deltau(xast,xtr,rhotv,rtr) 
     tauuold=tauu 
     tauu=(rhot-rb)/vlight + dtauu 
     ttold=tt 
     tt=tb-tauu 
! convergence control                                                   
!        WRITE(*,*)tt,tt-ttold,tauu, tauu-tauuold                       
     IF(abs(tauu-tauuold).lt.ept) GOTO 19 
  ENDDO
! too many iterations                                                   
  WRITE(*,*)' slow convergence up leg time ',tt-ttold, tauu-tauuold 
19 CONTINUE 
!      WRITE(*,*)' tauu at convergence ', tauu
!      WRITE(*,*)' rtr ', rtr
! compute relative velocity between asteroid and transmitter            
  vtr=xea(4:6)+ytr    
  rhotdv=xast(4:6)-vtr   
  rhotd=prscal(rhotv,rhotdv)/rhot 
! ==========================================================            
! compute distance                                                      
  r=0.5d0*(tauu+taud+(tdifft-tdiffr)/86400.d0)*vlight 
! compute relative frequency shift (up to 1/c^3 level);                 
! solar velocity at the bounce epochs and the asteroid barycentric      
! velocity (vastssb)                                                    
  ifla=2 
  CALL earcar(tt,xsun,ifla) 
  vtrssb=xsun(4:6) + vtr 
  vtrssb2=vtrssb(1)*vtrssb(1)+vtrssb(2)*vtrssb(2)+ vtrssb(3)*vtrssb(3)
  levelc=rhotd+rhord 
  scal1=prscal(rhotv,vtrssb)/rhot/vlight 
  scal2=prscal(rhorv,vastssb)/rhor/vlight 
  scal3=0.5d0*(vtrssb2-vressb2)+                                    &
     &      gms*((1.d0/vsize(rtr))-(1.d0/vsize(rre)))                   
  v=0.5d0*(levelc+rhotd*scal1*(1.d0+scal1)                          &
     &               -rhord*scal2*(1.d0-scal2)                          &
     &               -(rhotd*rhord*(1.d0+scal1-scal2)                   &
     &               -scal3*(1.d0-(levelc/vlight)))/vlight)             
! - get GR and tropospheric corrections to Doppler                      
! -- GR stuff                                                           
  rast=vsize(xast) 
  rastd=prscal(xast,xast(4))/rast 
! a) upleg piece                                                        
  rsta=vsize(rtr) 
  rstad=prscal(rtr,vtr)/rsta 
  call deldop1(rast,rastd,rhot,rhotd,rsta,rstad,deldoptr) 
! b) downleg piece                                                      
  rsta=vsize(rre) 
  rstad=prscal(rre,vre)/rsta 
  call deldop1(rast,rastd,rhor,rhord,rsta,rstad,deldopre) 
  v=v+0.5d0*(deldoptr+deldopre) 
! -- troposheric stuff                                                  
! a) at transmit passage -->                                            
  rgeo=vsize(xtr) 
  rgeod=prscal(xtr,ytr)/rgeo 
  call deldop2(xast,xtr,ytr,rgeo,rgeod,rhotv,rhotdv,rhot,rhotd,deldoptr)
! b) at receive passage <--                                             
  rgeo=vsize(xre) 
  rgeod=prscal(xre,yre)/rgeo 
  call deldop2(xast,xre,yre,rgeo,rgeod,rhorv,rhordv,rhor,rhord,deldopre)
  v=v+0.5d0*(deldoptr+deldopre) 
! rem interplanetary environment effects (e^- plasma) neglected         
! ==========================================================            
! Derivatives                                                           
  IF(ider.eq.0)RETURN 
! derivs of r,rdot wrt cartesian                                        
  do i=1,3 
! d(range)/d(r)                                                         
     drdx(i)=(rhotv(i)/rhot +rhorv(i)/rhor)/2.d0 
! d(range)/d(v) = 0                                                     
     drdx(i+3)=0 
! d(range-rate)/d(r)                                                    
     dvdx(i)=-((rhotd*rhotv(i)/rhot-rhotdv(i))/rhot +               &
     &             (rhord*rhorv(i)/rhor-rhordv(i))/rhor)/2.d0           
! d(range-rate)/d(v)= c * d(range)/d(r)                                 
     dvdx(i+3)=drdx(i) 
  enddo
! derivs of cartesian with respect to elements               
  do i=1,nd 
     drde(i)=DOT_PRODUCT(drdx,dxdpar(1:6,i))  
     dvde(i)=DOT_PRODUCT(dvdx,dxdpar(1:6,i))  
  enddo
  RETURN 
CONTAINS
! ======================================================================
! DELTAU - "small" corrections to radar time of flight                  
! ======================================================================
!     xast - asteroid position, heliocentric                            
!     xsta - station position, relative to Earth center                 
!     rho - asteroid position, relative to station                      
!     r - station position, heliocentric                                
  DOUBLE PRECISION FUNCTION deltau(xast,xsta,rho,r) 
    double precision xast(3),xsta(3),rho(3),r(3) 
    double precision vsize 
    double precision rsta,e,p,q,cosz,cotz,deltau1,deltau2 
!     double precision sinha,sin2ha,fghz,ampli,finte1,fun1,fun2         
!     double precision aprim,bprim,deltau3,omeg1,alpha                  
! ================================                                      
! Relativistic delay                                                    
    e=vsize(r) 
    p=vsize(xast) 
    q=vsize(rho) 
    deltau1=2d0*gms*log(abs((e+p+q)/(e+p-q)))/vlight**3 
! Earth ionospheric/tropospheric delay                                  
! ref. EM Standish, A&A 233, 252 (1990)                                 
    rsta=vsize(xsta) 
    if(rsta.lt.1d-12)then 
       write(*,*)'deltau: radar station at geocenter!' 
       cosz=1d0-1d-8 
    else 
       cosz=prscal(xsta,rho)/rsta/q ! cosz=prscag(3,xsta,rho)/rsta/q 
    endif
    if(cosz.eq.1d0)cosz=cosz-1d-8 
    cotz=cosz/sqrt(1d0-cosz**2) 
    deltau2=(7d-9/86400d0)/(cosz+1.4d-3/(4.5d-2+cotz)) 
! Interplanetary medium propagation effect                              
! ref. EM Standish, A&A 233, 252 (1990) with constants of DE118         
! rem. a more precise model might be needed here; e.g.                  
!      Muhleman & Anderson, ApJ 247, 1093 (1981) or newer               
!      alpha=q/e                                                        
!      cosha=(r(1)*rho(1)+r(2)*rho(2)+r(3)*rho(3))/e/q                  
!      sin2ha=1.d0-cosha*cosha                                          
!      sinha=dsqrt(sin2ha)                                              
! X- or S-band;                                                         
!c !!! Information about the frequency is not passed here at the        
!c     moment; it should be decided manually !!!                        
!c      fghz=2.38d0                                                     
!      fghz=8.51d0                                                      
!      ampli=(2.01094d-8)/fghz/fghz/e/86400.d0                          
!      aprim=1.237265d-6                                                
!      bprim=9.524021d0                                                 
!      finte1=(datan((alpha+cosha)/sinha)-datan(cosha/sinha))/sinha     
!      fun1=bprim*finte1                                                
!      omeg1=1.d0+alpha*(2.d0*cosh+alpha)                               
!      fun2=aprim*(0.25d0*((alpha+cosha)*((1.d0/omeg1)+(1.5d0/sin2ha))  
!     .           /omeg1-cosha*(1.d0+(1.5d0/sin2ha)))/sin2h             
!     .           +0.375d0*finte1/sin2ha/sin2ha)/(e**4)                 
!      deltau3=ampli*(fun1+fun2)                                        
! Add 'em up                                                            
!c      deltau=deltau1+deltau2+deltau3                                  
    deltau=deltau1+deltau2 
    return 
  END FUNCTION deltau
! ======================================================================
! DELDOP1 - "small" corrections to radar-rate measurements              
! ======================================================================
  SUBROUTINE deldop1(p,pdot,q,qdot,e,edot,deldop) 
    double precision p,pdot,q,qdot,e,edot,deldop 
    double precision brac1,brac2 
! ================================                                      
! relativistic range-rate correction                                    
    brac1=-q*(edot+pdot)+qdot*(e+p) 
    brac2=((e+p)**2)-(q**2) 
    deldop=4.d0*gms*brac1/brac2/(vlight**2) 
  END SUBROUTINE deldop1
! ======================================================================
! DELDOP2 - "small" corrections to radar-rate measurements              
! ======================================================================
!     xast(6) - asteroid r & v heliocentric                             
!     xtr(3),ytr(3) - station r & v relative to Earth center            
!     rsta,drsta - |xtr| & d|xtr|/dt                                    
!     rhov(3),drhov(3) - asteroid r & v relative to station             
!     rho,drho - |rhov| & d|rhov|/dt                                    
  SUBROUTINE deldop2(xast,xsta,vsta,rsta,drsta,rhov,drhov,rho,drho,deldop)
    double precision xast(6),xsta(3),vsta(3),rhov(3),drhov(3) 
    double precision rsta,drsta,rho,drho,deldop 
    double precision cosz,sinz,cotz,dcoszdt,phiz 
    double precision brac,brac1,brac2,scal1,scal2 
! ================================                                      
! rate of change of the Earth ionospheric/tropospheric ~ Doppler shift  
    if(rsta.lt.1d-12)then 
       write(*,*)'deltau: radar station at geocenter!' 
       cosz=1d0-1d-8 
    else 
       cosz=(rhov(1)*xsta(1)+rhov(2)*xsta(2)+rhov(3)*xsta(3))/rsta/rho 
    endif
    sinz=sqrt(1.d0-cosz**2) 
    cotz=cosz/sinz 
    brac=0.045d0+cotz 
    brac1=1.d0-(0.0014d0/brac/brac/(sinz**3)) 
    brac2=cosz+(0.0014d0/brac) 
    phiz=-vlight*(7.d-9/86400.d0)*brac1/(brac2**2) 
    scal1=drhov(1)*xsta(1)+drhov(2)*xsta(2)+drhov(3)*xsta(3) 
    scal2=rhov(1)*vsta(1)+rhov(2)*vsta(2)+rhov(3)*vsta(3) 
    dcoszdt=((scal1+scal2)/rho/rsta)-cosz*((drho/rho)+(drsta/rsta)) 
    deldop=phiz*dcoszdt 
  END SUBROUTINE deldop2

END SUBROUTINE r_rdot


! ===================================================================== 
! OUTOBC                                                                
! ===================================================================== 
!  output of predicted observation, possibly with confidence ellipse    
!   input: iun   = output unit                                          
!          type  = observation type                                     
!          ids = station code                                           
!          t1 = time of observation (UTC)                               
!          alpha, delta, hmagn = observation                            
!          adot,ddot = proper motion                                    
!          elo,dis = elongation, distance from Earth                    
!          elmoon = lunar elongation 
!          elsun = Sun elevation above horizon
!          icov  = 1 for observations only, 2 to add confidence ellipse 
!          gamad,sig,axes = covariance matrix, sigmas along axes        
!                      (only for icov=2, otherwise dummy)               
! ===================================================================== 
SUBROUTINE outobc(iun,type,ids,t1,alpha,delta,hmagn,adot,ddot,    &
     &     elo,dis,icov,gamad,sig,axes,elmoon,gallat,gallon,elev,dsun,elsun,pha)  
  USE fund_const                               
  implicit none 
! needs AU value in km: from fund_cons
! output unit, station code, obs. type                                  
      integer iun,ids
      CHARACTER*(1) type
! observations                                                          
      DOUBLE PRECISION t1,alpha,delta,hmagn,adot,ddot,elo,dis,cosangzen, airmass 
      DOUBLE PRECISION, INTENT(IN), OPTIONAL :: elmoon, gallat,gallon
      DOUBLE PRECISION, INTENT(IN), OPTIONAL :: elev,dsun,elsun,pha
! covariance                                                            
      integer icov 
      DOUBLE PRECISION gamad(2,2),axes(2,2),sig(2), velsiz, pa, pad, cvf
! ================end interface===============================          
      DOUBLE PRECISION princ 
      integer i,j 
! time variables                                                        
      integer ideg,iday,imonth,iyear,ihour,imin,isec,ln,truncat 
      double precision hour,minu,sec 
      CHARACTER*22 timstr 
      CHARACTER*19 tmpstr 
      CHARACTER*12 rastri,rdstri 
      CHARACTER*1 signo
      LOGICAL fail
! convert time                                                          
      CALL mjddat(t1,iday,imonth,iyear,hour) 
! convert hour to 12:12:12                                              
      ihour=truncat(hour,1d-7) 
      minu=(hour-ihour)*60.d0 
      imin=truncat(minu,1d-5) 
      sec=(minu-imin)*60.d0 
      isec=truncat(sec,1d-3) 
      WRITE(timstr,192) iyear,imonth,iday,ihour,imin,isec,sec-isec 
  192 FORMAT(I4,'/',I2.2,'/',I2.2,1x,I2.2,':',I2.2,':',I2.2,f3.2) 
! =================== select by observation type ===================    
      IF(type.eq.'O'.or.type.eq.'S')THEN 
! %%%%%%%%%%%% ASTROMETRY %%%%%%%%%%%%%%%%                              
! convert RA                                                            
         alpha=princ(alpha) 
         CALL sessag(alpha*degrad/15.d0,signo,ihour,imin,sec) 
         IF(signo.eq.'-')STOP 'wrirms error: negative right ascension.' 
! prepare RA string                                                     
         WRITE(tmpstr,FMT='(F6.3)') sec 
         CALL rmsp(tmpstr,ln) 
         IF(ln.lt.6)tmpstr='0'//tmpstr 
         WRITE(rastri,130)ihour,imin,tmpstr 
  130    FORMAT(I2.2,':',I2.2,':',a6) 
! convert DEC                                                           
         CALL sessag(delta*degrad,signo,ideg,imin,sec) 
! prepare DEC string                                                    
         WRITE(tmpstr,FMT='(F5.2)') sec 
         CALL rmsp(tmpstr,ln) 
         IF(ln.lt.5)tmpstr='0'//tmpstr 
         WRITE(rdstri,170)signo,ideg,imin,tmpstr 
  170    FORMAT(A1,I2.2,1x,I2.2,1x,a5)
! Airmass calculation according Young (1994)
         IF(PRESENT(elev))THEN 
            cosangzen=cos(pig/2-elev)
            airmass=1.002432d0*cosangzen**2+0.148386d0*cosangzen+0.0096467d0
            airmass=airmass/(cosangzen**3+0.149864d0*cosangzen**2+0.0102963d0*cosangzen+0.000303978)
         ENDIF
         CALL angvcf('"/min',cvf,fail) 
         velsiz=SQRT((secrad*adot*cos(delta)/1440.d0)**2+(secrad*ddot/1440.d0)**2)
         pa=ATAN2(adot*cos(delta),ddot) 
         IF(pa.LT.0.D0) pa=pa+dpig
         pa=pa*degrad 
         IF(PRESENT(elmoon).and.PRESENT(gallat).and.PRESENT(gallon).and.PRESENT(elev).and.PRESENT(dsun).and.PRESENT(elsun))THEN 
            IF(elev.le.0) THEN
               write(iun,101)timstr,t1,ids,                                &
                    &        rastri,alpha*degrad,                                      &
                    &        rdstri,delta*degrad,                                      &
                    &        secrad*adot*cos(delta)/1440.d0,secrad*ddot/1440.d0,  &
                    &        secrad*adot*cos(delta)*50.d0/(3.d0*1440.d0),secrad*ddot*50.d0/(3.d0*1440.d0),  &
                    &        velsiz,pa,dsun,                                                &
                    &        dis,elsun*degrad,elo*degrad,elmoon*degrad,gallat*degrad,         &
                    &        gallon*degrad,hmagn,pha*degrad,elev*degrad
               write(*,101)timstr,t1,ids,                                  &
                    &        rastri,alpha*degrad,                                      &
                    &        rdstri,delta*degrad,                                      &
                    &        secrad*adot*cos(delta)/1440.d0,secrad*ddot/1440.d0,  &
                    &        secrad*adot*cos(delta)*50.d0/(3.d0*1440.d0),secrad*ddot*50.d0/(3.d0*1440.d0),  &
                    &        velsiz,pa,dsun,                                                &
                    &        dis,elsun*degrad,elo*degrad,elmoon*degrad,gallat*degrad,         &
                    &        gallon*degrad,hmagn,pha*degrad,elev*degrad
            ELSE IF(ids.eq.500.or.ids.eq.245.or.ids.eq.247.or.ids.eq.248.or.ids.eq.249.or.ids.eq.250) THEN
               write(iun,104)timstr,t1,ids,                                &
                    &        rastri,alpha*degrad,                                      &
                    &        rdstri,delta*degrad,                                      &
                    &        secrad*adot*cos(delta)/1440.d0,secrad*ddot/1440.d0,  &
                    &        secrad*adot*cos(delta)*50.d0/(3.d0*1440.d0),secrad*ddot*50.d0/(3.d0*1440.d0),  &
                    &        velsiz,pa,dsun,                                                &
                    &        dis,elo*degrad,elmoon*degrad,gallat*degrad,gallon*degrad, &
                    &        hmagn,pha*degrad
               write(*,104)timstr,t1,ids,                                  &
                    &        rastri,alpha*degrad,                                      &
                    &        rdstri,delta*degrad,                                      &
                    &        secrad*adot*cos(delta)/1440.d0,secrad*ddot/1440.d0,  &
                    &        secrad*adot*cos(delta)*50.d0/(3.d0*1440.d0),secrad*ddot*50.d0/(3.d0*1440.d0),  &
                    &        velsiz,pa,dsun,                                                &
                    &        dis,elo*degrad,elmoon*degrad,gallat*degrad,gallon*degrad, &
                    &        hmagn
            ELSE
               write(iun,105)timstr,t1,ids,                                &
                    &        rastri,alpha*degrad,                                      &
                    &        rdstri,delta*degrad,                                      &
                    &        secrad*adot*cos(delta)/1440.d0,secrad*ddot/1440.d0,  &
                    &        secrad*adot*cos(delta)*50.d0/(3.d0*1440.d0),secrad*ddot*50.d0/(3.d0*1440.d0),  &
                    &        velsiz,pa,dsun,                                                &
                    &        dis,elsun*degrad,elo*degrad,elmoon*degrad,gallat*degrad,         &
                    !     &        dis,elo*degrad,elmoon*degrad,gallat*degrad,         &
                    &        gallon*degrad,hmagn,pha*degrad,elev*degrad,airmass
               write(*,105)timstr,t1,ids,                                  &
                    &        rastri,alpha*degrad,                                      &
                    &        rdstri,delta*degrad,                                      &
                    &        secrad*adot*cos(delta)/1440.d0,secrad*ddot/1440.d0,  &
                    &        secrad*adot*cos(delta)*50.d0/(3.d0*1440.d0),secrad*ddot*50.d0/(3.d0*1440.d0),  &
                    &        velsiz,pa,dsun,                                                &
                    &        dis,elsun*degrad,elo*degrad,elmoon*degrad,gallat*degrad,         &
                    !     &        dis,elo*degrad,elmoon*degrad,gallat*degrad,         &
                    &        gallon*degrad,hmagn,pha*degrad,elev*degrad,airmass
101            format('Astrometric Observation Prediction'/                   &
                    &        'For ',a19,' (UTC); ',f12.5,'(MJD)'/                      &
                    &        'Observatory code = ',i4.4/                               &
                    &        'RA = ',a12,' (HH:MM:SS); ',f11.5,' (deg)'/               &
                    &        'DEC = ',a12,' (deg min sec); ',f11.5,' (deg)'/           &
                    &        'RA*cos(DEC)/DEC Apparent motion =',2(2x,f9.3),' (arcsec/min)',2(2x,f10.1),' (marcsec/s)'/    &
                    &        'Apparent motion velocity = ',f9.3,' arcsec/min)', ' and Position Angle = ',f7.3, ' (deg)'/ &
                    &        'Sun distance = ',f8.4,' (au)'/                           &
                    &        'Earth distance = ',f8.4,' (au)'/                         &
                    &        'Sun elevation above horizon = ',f7.2,' (deg)'/           &
                    &        'Solar elongation = ',f7.2,' (deg)','    Lunar elongation = ',f7.2,' (deg)'/                       &
                    &        'Galactic latitude = ',f7.2,' (deg)','   Galactic longitude = ',f7.2,' (deg)'/                     &
                    &        'Apparent magnitude = ',f5.2/                             &
                    &        'Phase angle = ',f7.2,' (deg)'/                           &
                    &        'Altitude = ',f7.2,' (deg)','    Airmass =   INF ')
104            format('Astrometric Observation Prediction'/                   &
                    &        'For ',a19,' (UTC); ',f12.5,'(MJD)'/                      &
                    &        'Observatory code = ',i4.4/                               &
                    &        'RA = ',a12,' (HH:MM:SS); ',f11.5,' (deg)'/               &
                    &        'DEC = ',a12,' (deg min sec); ',f11.5,' (deg)'/           &
                    &        'RA*cos(DEC)/DEC Apparent motion =',2(2x,f9.3),' (arcsec/min)',2(2x,f10.1),' (marcsec/s)'/    &
                    &        'Apparent motion velocity = ',f9.3,' arcsec/min)', ' and Position Angle = ',f7.3, ' (deg)'/ &
                    &        'Sun distance = ',f8.4,' (au)'/                           &
                    &        'Earth distance = ',f8.4,' (au)'/                         &
                    &        'Solar elongation = ',f7.2,' (deg)','    Lunar elongation = ',f7.2,' (deg)'/                       &
                    &        'Galactic latitude = ',f7.2,' (deg)','   Galactic longitude = ',f7.2,' (deg)'/                     &
                    &        'Apparent magnitude = ',f5.2/                             &
                    &        'Phase angle = ',f7.2,' (deg)'/                           &
                    &        'Altitude = N.A.             Airmass =   N.A.')
105            format('Astrometric Observation Prediction'/                   &
                    &        'For ',a19,' (UTC); ',f12.5,'(MJD)'/                      &
                    &        'Observatory code = ',i4.4/                               &
                    &        'RA = ',a12,' (HH:MM:SS); ',f11.5,' (deg)'/               &
                    &        'DEC = ',a12,' (deg min sec); ',f11.5,' (deg)'/           &
                    &        'RA*cos(DEC)/DEC Apparent motion =',2(2x,f9.3),' (arcsec/min)',2(2x,f10.1),' (marcsec/s)'/    &
                    &        'Apparent motion velocity = ',f9.3,' arcsec/min)', ' and Position Angle = ',f7.3, ' (deg)'/ &
                    &        'Sun distance = ',f8.4,' (au)'/                           &
                    &        'Earth distance = ',f8.4,' (au)'/                         &
                    &        'Sun elevation above horizon = ',f7.2,' (deg)'/           &
                    &        'Solar elongation = ',f7.2,' (deg)','    Lunar elongation = ',f7.2,' (deg)'/                       &
                    &        'Galactic latitude = ',f7.2,' (deg)','   Galactic longitude = ',f7.2,' (deg)'/                     &
                    &        'Apparent magnitude = ',f5.2/                             &
                    &        'Phase angle = ',f7.2,' (deg)'/                           &
                    &        'Altitude = ',f7.2,' (deg)','    Airmass =',f8.3)
            ENDIF
         ELSE IF(PRESENT(elmoon).and.PRESENT(dsun))THEN 
            write(iun,103)timstr,t1,ids,                                &
                 &        rastri,alpha*degrad,                                      &
                 &        rdstri,delta*degrad,                                      &
                 &        secrad*adot*cos(delta)/1440.d0,secrad*ddot/1440.d0,  &
                 &        secrad*adot*cos(delta)*50.d0/(3.d0*1440.d0),secrad*ddot*50.d0/(3.d0*1440.d0),  &
                 &        velsiz,pa,dsun,                                                &
                 &        dis,elo*degrad,elmoon*degrad,hmagn,pha*degrad
            write(*,103)timstr,t1,ids,                                  &
                 &        rastri,alpha*degrad,                                      &
                 &        rdstri,delta*degrad,                                      &
                 &        secrad*adot*cos(delta)/1440.d0,secrad*ddot/1440.d0,  &
                 &        secrad*adot*cos(delta)*50.d0/(3.d0*1440.d0),secrad*ddot*50.d0/(3.d0*1440.d0),  &
                 &        velsiz,pa,dsun,                                                &
                 &        dis,elo*degrad,elmoon*degrad,hmagn,pha*degrad
103         format('Astrometric Observation Prediction'/                   &
                 &        'For ',a19,' (UTC); ',f12.5,'(MJD)'/                      &
                 &        'Observatory code = ',i4.4/                               &
                 &        'RA = ',a12,' (HH:MM:SS); ',f11.5,' (deg)'/               &
                 &        'DEC = ',a12,' (deg min sec); ',f11.5,' (deg)'/           &
                 &        'RA*cos(DEC)/DEC Apparent motion =',2(2x,f9.3),' (arcsec/min)',2(2x,f10.1),' (marcsec/s)'/    &
                 &        'Apparent motion velocity = ',f9.3,' arcsec/min)', ' and Position Angle = ',f7.3, ' (deg)'/ &
                 &        'Sun distance = ',f8.4,' (au)'/                           &
                 &        'Earth distance = ',f8.4,' (au)'/                         &
                 &        'Solar elongation = ',f7.2,' (deg)'/                      &
                 &        'Lunar elongation = ',f7.2,' (deg)'/                      &
                 &        'Apparent magnitude = ',f5.2/                             &
                 &        'Phase angle = ',f7.2,' (deg)')              
         ELSE
            write(iun,221)timstr,t1,ids,                                &
                 &        rastri,alpha*degrad,                                      &
                 &        rdstri,delta*degrad,                                      &
                 &        secrad*adot/24.d0,secrad*ddot/24.d0,                      &
                 &        velsiz,pa,dsun,                                                &
                 &        dis,elo*degrad,hmagn,pha*degrad
            write(*,221)timstr,t1,ids,                                  &
                 &        rastri,alpha*degrad,                                      &
                 &        rdstri,delta*degrad,                                      &
                 &        secrad*adot/24.d0,secrad*ddot/24.d0,                      &
                 &        velsiz,pa,dsun,                                                &
                 &        dis,elo*degrad,hmagn,pha*degrad
221         format('Astrometric Observation Prediction'/                   &
                 &        'For ',a19,' (UTC); ',f12.5,'(MJD)'/                      &
                 &        'Observatory code = ',i4.4/                               &
                 &        'RA = ',a12,' (HH:MM:SS); ',f11.5,' (deg)'/               &
                 &        'DEC = ',a12,' (deg min sec); ',f11.5,' (deg)'/           &
                 &        'RA/DEC Apparent motion =',2(2x,f9.2),' (arcsec/hour)'/   &
                 &        'Apparent motion velocity = ',f9.3,' arcsec/min)', ' and Position Angle = ',f7.3, ' (deg)'/ &
                 &        'Sun distance = ',f8.4,' (au)'/                           &
                 &        'Earth distance = ',f8.4,' (au)'/                         &
                 &        'Solar elongation = ',f7.2,' (deg)'/                      &
                 &        'Apparent magnitude = ',f5.2/                            &
                 &        'Phase angle = ',f7.2,' (deg)')                           
         ENDIF
         IF(icov.eq.1)RETURN 
! rescaling in arcsec                                                   
         do  i=1,2 
            sig(i)=sig(i)*secrad 
            do  j=1,2 
               gamad(i,j)=gamad(i,j)*secrad**2 
            enddo 
         enddo 
         write(iun,201)(sig(j),(axes(i,j),i=1,2),j=1,2) 
         write(*,201)(sig(j),(axes(i,j),i=1,2),j=1,2) 
  201    format(                                                        &
     &'Size and orientation of 1-sigma uncertainty ellipse'/            &
     &'Short axis : Size = ',1p,g12.6 ,' (arcsec); Direction = ',       &
     & 0p,2(1x,f8.5)/                                                   &
     &'Long axis : Size = ',1p,g12.6 ,' (arcsec); Direction = ',        &
     & 0p,2(1x,f8.5))                                                   
      ELSEIF(type.eq.'R'.or.type.eq.'V')THEN 
! %%%%%%%%%%%% RADAR %%%%%%%%%%%%%%%%                                   
         write(iun,102)t1,ids,alpha*aukm,delta*aukm 
         write(*,102)t1,ids,alpha*aukm,delta*aukm 
  102    format('time, MJD=',f13.6,'  station=',i4/                     &
     &       ' range (KM)         = ',f16.5/                            &
     &       ' range rate (KM/DAY)=  ',f15.5)                           
         IF(icov.eq.1)RETURN 
! rescaling in km, km/day                                               
         do  i=1,2 
            sig(i)=sig(i)*aukm 
            do  j=1,2 
               gamad(i,j)=gamad(i,j)*aukm**2 
            enddo 
         enddo 
         write(iun,202)(sig(j),(axes(i,j),i=1,2),j=1,2) 
         write(*,202)(sig(j),(axes(i,j),i=1,2),j=1,2) 
  202    format(' in the range (KM), range-rate (KM/DAY) plane'/        &
     &          ' sigma1 = ',1p,g14.7 ,' axis1 = ',2(1x,g12.5)/          &
     &          ' sigma2 = ',1p,g14.7 ,' axis2 = ',2(1x,g12.5))          
      ELSEIF(type.eq.'P')THEN 
! %%%%%%%%%%%% PROPER MOTION %%%%%%%%%%%%%%%%                           
         write(iun,109)t1,ids,alpha*secrad/24.d0,delta*secrad/24.d0 
         write(*,109)t1,ids,alpha*secrad/24.d0,delta*secrad/24.d0 
  109    format('time, MJD=',f13.6,'  station=',i4/                     &
     &       ' RA motion (arcsec/hour)     = ',f9.2/                    &
     &       ' DEC motion (arcsec/hour)    = ',f9.2)                    
         IF(icov.eq.1)RETURN 
! rescaling in arcsec/hour                                              
         do  i=1,2 
            sig(i)=sig(i)*secrad/24.d0 
            do  j=1,2 
               gamad(i,j)=gamad(i,j)*(secrad/24.d0)**2 
            enddo 
         enddo 
         write(iun,209)(sig(j),(axes(i,j),i=1,2),j=1,2) 
         write(*,209)(sig(j),(axes(i,j),i=1,2),j=1,2) 
  209    format('sigma1 (arcsec/hr) = ',1p,g14.7 ,' axis1 = ',2(1x,g12.5)/&
     &          'sigma2 (arcsec/hr) = ',1p,g14.7 ,' axis2 = ',2(1x,g12.5))
      ELSE 
         WRITE(*,*)'outobc: type=',type,' not understood' 
      ENDIF 
      return 
      END SUBROUTINE outobc




END MODULE pred_obs  

! Copyright (C) 1997 by Mario Carpino (carpino@brera.mi.astro.it)       
! Version: February 24, 1997                                            
!                                                                       
!  ***************************************************************      
!  *                                                             *      
!  *                          A B E R 1                          *      
!  *                                                             *      
!  *          Correzione approssimata per aberrazione            *      
!  *                  stellare e/o planetaria                    *      
!  *                                                             *      
!  ***************************************************************      
!                                                                       
!                                                                       
! INPUT:    XREL(3)   -  Posizione relativa vera del corpo osservato (UA
!           VREL(3)   -  Velocita` relativa (UA/d)                      
!                                                                       
! OUTPUT:   XCOR(3)   -  Posizione relativa apparente (tenendo conto del
!                        aberrazione)                                   
!                                                                       
! NOTA: in generale per ottenere la correzione completa (comprendente le
!       cosiddette aberrazioni "stellare" + "planetaria") bisogna che VR
!       sia la velocita` relativa del corpo osservato rispetto all'osser
!       tore:                                                           
!                 VREL  =   V(pianeta) - V(osservatore)                 
!                                                                       
!       Se si vuole ottenere solo la correzione per l'aberrazione "stell
!       bisogna porre VREL =  - V(osservatore).                         
! 
SUBROUTINE aber1(xrel,vrel,xcor) 
  USE fund_const
  IMPLICIT NONE
  DOUBLE PRECISION :: xrel(3),vrel(3)
  DOUBLE PRECISION :: xcor(3) ! removed INTENT to allow xrel in the same location as xcor
  DOUBLE PRECISION ro,vsize,dt
!  distanza
  ro=vsize(xrel) 
!  effetto di ritardo                                                   
  dt=ro/vlight
  xcor=xrel-dt*vrel 

END SUBROUTINE aber1

! Copyright (C) 20087 by Fabrizio Bernardi (bernardi@adams.dm.unipi.it)  
! Version: September 9, 2008                                            
!                                                                       
!  ***************************************************************      
!  *                                                             *      
!  *                          A B E R 2                          *      
!  *                                                             *      
!  *          Correzione approssimata per aberrazione            *      
!  *       stellare e/o planetaria al secondo ordine in t        *
!  *                                                             *      
!  ***************************************************************      
!                                                                       
!                                                                       
! INPUT:    XREL(3)   -  Posizione relativa vera del corpo osservato (UA
!           VREL(3)   -  Velocita` relativa (UA/d)                      
!                                                                       
! OUTPUT:   XCOR(3)   -  Posizione relativa apparente (tenendo conto del
!                        aberrazione)                                   
!                                                                       
! NOTA: in generale per ottenere la correzione completa (comprendente le
!       cosiddette aberrazioni "stellare" + "planetaria") bisogna che VR
!       sia la velocita` relativa del corpo osservato rispetto all'osser
!       tore:                                                           
!                 VREL  =   V(pianeta) - V(osservatore)                 
!                                                                       
!       Se si vuole ottenere solo la correzione per l'aberrazione "stell
!       bisogna porre VREL =  - V(osservatore).                         
! 

SUBROUTINE aber2(xrel,xast,xoss,xcor) 
  USE fund_const
  USE ever_pitkin
  IMPLICIT NONE
  DOUBLE PRECISION :: xrel(3),xast(6),xoss(6)
  DOUBLE PRECISION :: xcor(6),xastcor(6)    ! removed INTENT to allow xrel in the same location as xcor
  DOUBLE PRECISION r,ro,rdotdot(3),vsize,dt,t0
  t0=0
!  distance
  ro=vsize(xrel) 
!  delay effect    
                                               
  dt=ro/vlight
! first call of fser_propag. Determine new ro and dt
  CALL fser_propag(xast(1:3),xast(4:6),t0,-dt,gms,xastcor(1:3),xastcor(4:6))
  ro=vsize(xastcor(1:3)-xoss(1:3)) 
  dt=ro/vlight
! second call. Determine acceleration rdotdot and final corrected coordinates.
  CALL fser_propag(xast(1:3),xast(4:6),t0,-dt,gms,xastcor(1:3),xastcor(4:6))
  xcor(1:3)=xastcor(1:3)-xoss(1:3)
  xcor(4:6)=xastcor(4:6)-xoss(4:6)
END SUBROUTINE aber2


! Copyright (C) 2012 by Davide Bracali Cioci (bracalicioci@spacedys.com)  
! Version: February 27, 2012                                                                                                                  
!  ***************************************************************      
!  *                                                             *      
!  *                          RISE_SET_ITER                      *      
!  *                                                             *      
!  *         Computes times of rise and set of an object         *      
!  *              with respect to a fixed altitude               *
!  *        through linear propagation of right ascension        *
!  *        and declination to the times of rise and set         *
!  *        computed in first approximation by rise_and_set      *
!  *         (taking account of atmospheric refraction)          *
!  *                                                             *      
!  ***************************************************************      
!                                                                       
!                                                                       
! INPUT:    alpha,delta,adot,ddot   -  right ascension, declination and
!                                      angular velocities 
!           longitude,latitude      -  longitude and latitude of the 
!                                      station (rad)
!           altitude                -  fixed altitude for the computation 
!                                      of times of rise and set (rad) 
!           tref                    -  mjd                                   
! OUTPUT:   t_rise,t_set,t_transit     (day)
!           visib                   -  flag for the visibility of the object
!                                      visib = 0: object never rises (not visible)
!                                      visib = 1: object never sets
!                                      visib = 2: object rises and sets
! OPTIONAL OUTPUT: alphavec,deltavec  -  bidimensional vectors such that
!                                        alphavec(1)=alpha_rise 
!                                        alphavec(2)=alpha_set
!                                        deltavec(1)=delta_rise
!                                        deltavec(2)=delta_set
!                                        after linear propagation

SUBROUTINE rise_set_iter(alpha,delta,adot,ddot,tref,longitude,latitude,altitude,&
     & t_rise,t_set,t_transit,visib,alphavec,deltavec)
  USE iers_ser
  USE fund_const
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT(IN) :: alpha,delta,adot,ddot       ! right ascension,declination and
                                                              ! angular velocities
  DOUBLE PRECISION, INTENT(IN) :: longitude,latitude          ! longitudine, latitude of the station
  DOUBLE PRECISION, INTENT(IN) :: altitude                    ! altitude of object
  DOUBLE PRECISION, INTENT(IN) :: tref                        ! local time
  DOUBLE PRECISION, INTENT(OUT) :: t_rise,t_set,t_transit                    ! rise, set and transit of object
                                                                             ! corresponding to the altitude in input
  DOUBLE PRECISION, DIMENSION(2), INTENT(OUT), OPTIONAL :: alphavec,deltavec ! r.a. and dec. at rise and set
  INTEGER, INTENT(OUT) :: visib                                              ! flag for the visibility of the object
    
  !======================END INTERFACE========================
  
  ! variables used inside two loops
  DOUBLE PRECISION :: alpha_loop1,alpha_loop2,delta_loop1,delta_loop2   
  DOUBLE PRECISION :: t_rise_aux,t_set_aux,t_transit_aux
  DOUBLE PRECISION :: sec_loop1,sec_loop2,tref_loop1,tref_loop2
  DOUBLE PRECISION, PARAMETER :: min=1.d0/1440.d0,sec=1.d0/86400.d0
  INTEGER :: visib_aux 
  
  LOGICAL :: first ! to check that at least one iteration is performed
  
  ! compute starting guess for t_rise and t_set to be used inside two loops
  CALL rise_and_set(alpha,delta,tref,longitude,latitude,altitude,t_rise_aux,t_set_aux,t_transit,visib)
  
  ! set initial values for loop variables
  alpha_loop1=alpha
  alpha_loop2=alpha
  delta_loop1=delta
  delta_loop2=delta
  sec_loop1=tref-FLOOR(tref) 
  sec_loop2=sec_loop1
  tref_loop1=tref
  tref_loop2=tref
  
  first=.TRUE. ! set logical variable to start the loop 
  DO WHILE(first.OR.(ABS(t_rise_aux-sec_loop1).GT.sec))
     alpha_loop1=alpha_loop1+adot*(t_rise_aux-sec_loop1) ! linear propagation of the right ascension to the t_rise
     delta_loop1=delta_loop1+ddot*(t_rise_aux-sec_loop1) ! linear propagation of the declination to the t_rise
     sec_loop1=t_rise_aux
     ! if time of rise is negative we recompute it at the previous day
     IF(t_rise_aux.LT.0)THEN
        tref_loop1=tref_loop1-1.d0  
        sec_loop1=t_rise_aux+1.d0
     END IF
     ! compute new time of rise
     CALL rise_and_set(alpha_loop1,delta_loop1,tref_loop1,longitude,latitude,altitude,t_rise_aux,t_set_aux,t_transit_aux,visib_aux)
     first=.FALSE.
     ! computation at the previous day is already done
     IF(t_rise_aux.LT.0)t_rise_aux=t_rise_aux+1.d0 
     !!!! exit loop if new time differ from the previous by one second  
  ENDDO
  
  t_rise=t_rise_aux
  
  first=.TRUE. ! set logical variable to start the loop
  DO WHILE(first.OR.(ABS(t_set_aux-sec_loop2).GT.sec))
     alpha_loop2=alpha_loop2+adot*(t_set_aux-sec_loop2) ! linear propagation of the right ascension to the t_set
     delta_loop2=delta_loop2+ddot*(t_set_aux-sec_loop2) ! linear propagation of the declination to the t_set
     sec_loop2=t_set_aux
     ! if time of set is greater than 1 we recompute it at the next day
     IF(t_set_aux.GT.1.d0)THEN
        tref_loop1=tref_loop1+1.d0
        sec_loop2=t_set_aux-1.d0
     END IF
     ! compute new time of set
     CALL rise_and_set(alpha_loop2,delta_loop2,tref_loop2,longitude,latitude,altitude,t_rise_aux,t_set_aux,t_transit_aux,visib_aux)
     first=.FALSE.
     ! computation at the next day is already done
     IF(t_set_aux.GT.1.d0)t_set_aux=t_set_aux-1.d0  
     !!!! exit loop if new time differ from the previous by one second
  END DO
  
  t_set=t_set_aux
  ! store right ascensions and declinations of rise and set after linear propagation 
  IF(PRESENT(alphavec))THEN
     alphavec(1)=alpha_loop1
     alphavec(2)=alpha_loop2
  END IF
  IF(PRESENT(deltavec))THEN
     deltavec(1)=delta_loop1
     deltavec(2)=delta_loop2
  END IF

END SUBROUTINE rise_set_iter



! Copyright (C) 2012 by Davide Bracali Cioci (bracalicioci@spacedys.com)  
! Version: February 27, 2012
!  ***************************************************************      
!  *                                                             *      
!  *                          RISE_AND_SET                       *      
!  *                                                             *      
!  *         Computes times of rise and set of an object         *      
!  *              with respect to a fixed altitude               *
!  *         (taking account of atmospheric refraction)          *
!  *           without rescaling times of rise and set           *
!  *            between 0 and 1 day (i.e. 0 and 2*pi)            *
!  *                                                             *      
!  ***************************************************************      
!                                                                       
!                                                                       
! INPUT:    alpha,delta             -  right ascension and declination 
!           longitude,latitude      -  longitude and latitude of the 
!                                      station (rad)
!           altitude                -  fixed altitude for the computation 
!                                      of times of rise and set (rad) 
!           tref                    -  mjd                                   
! OUTPUT:   t_rise,t_set,t_transit     (day)
!           visib                   -  flag for the visibility of the object
!                                      visib = 0: object never rises (not visible)
!                                      visib = 1: object never sets
!                                      visib = 2: object rises and sets


SUBROUTINE rise_and_set(alpha,delta,tref,longitude,latitude,altitude,t_rise,t_set,t_transit,visib)
  USE iers_ser
  USE fund_const
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT(IN) :: alpha,delta                 ! right ascension,declination
  DOUBLE PRECISION, INTENT(IN) :: longitude,latitude          ! longitudine, latitude of the station
  DOUBLE PRECISION, INTENT(IN) :: altitude                    ! altitude of object
  DOUBLE PRECISION, INTENT(IN) :: tref                        ! local time
  DOUBLE PRECISION, INTENT(OUT) :: t_rise,t_set,t_transit     ! rise, set and transit of object
                                                              ! corresponding to the altitude in input
  INTEGER, INTENT(OUT) :: visib                               ! flag for the visibility of the object
    
  !======================END INTERFACE========================

  INTEGER :: mjd_aux                       ! mjd corresponding to the tref in input
  DOUBLE PRECISION :: ut                   ! universal time
  DOUBLE PRECISION :: lha,cos_lha          ! local hour angle and its cosine
  DOUBLE PRECISION :: gmst,gmstd1,gmstd2   ! Greenwich Mean Sidereal Time and derivatives
  DOUBLE PRECISION :: atm_ref,atm_ref_delta,atm_ref_alpha           ! taking account atmospheric refraction
  DOUBLE PRECISION :: true_alt,true_delta,true_alpha
  !DOUBLE PRECISION :: cos_parall,sin_parall ! cosine and sine of the parallactic angle between 
                                            ! celestial pole and zenith as seen from the object position
  DOUBLE PRECISION, PARAMETER :: ben1=1.351520851d-3,ben2=0.6069468169
  
  mjd_aux=FLOOR(tref)
  !Bennet formula
  atm_ref=1.d0/TAN(altitude+ben1/(altitude*ben2+1-0.5*pig*ben2))
  !cos_parall=(SIN(latitude)-SIN(delta)*SIN(altitude))/(COS(delta)*COS(altitude))
  true_alt=altitude-(atm_ref/60.d0)*radeg
  !Correction for R.A. and Dec. from the article "An Accurate Method for 
  !Computing Atmospheric Refraction",1996,Stone, R. C
  true_delta=delta!-(atm_ref/60.d0)*radeg*cos_parall
  true_alpha=alpha!-(atm_ref_alpha/60.d0)*radeg
  CALL gmstd(mjd_aux,0.d0,gmst,gmstd1,gmstd2,0)       ! compute Greenwich Mean Sidereal Time
  ut=MOD(true_alpha-gmst-longitude,dpig)                   ! compute universal time
  IF(ut.lt.0)THEN
     ut=ut+dpig                ! corresponding value of UT in [0,dpig]
  END IF
  ! compute cosine of local hour angle
  cos_lha=(SIN(true_alt)-SIN(latitude)*SIN(true_delta))/(COS(latitude)*COS(true_delta)) 
  IF(cos_lha.ge.1)THEN
     lha=0.d0      ! object always below the altitude
     visib=0       ! object never rises
  ELSEIF(cos_lha.le.-1)THEN
     lha=pig       ! object always above the altitude
     visib=1       ! object never sets
  ELSE
     lha=acos(cos_lha)
     visib=2
  END IF
  !sin_parall=COS(latitude)*SIN(lha)/COS(altitude)
  !true_alpha=alpha-(atm_ref/60.d0)*radeg*sin_parall/COS(delta)
  !CALL gmstd(mjd_aux,0.d0,gmst,gmstd1,gmstd2,0)       ! compute Greenwich Mean Sidereal Time
  !ut=MOD(true_alpha-gmst-longitude,dpig)                   ! compute universal time
  !IF(ut.lt.0)THEN
  !   ut=ut+dpig                ! corresponding value of UT in [0,dpig]
  !END IF
  ! compute times of transit, rise and set
  t_transit=ut             
  t_rise=ut-lha
  !IF(t_rise.lt.0)THEN
  !   t_rise=t_rise+dpig
  !END IF
  t_set=ut+lha
  !IF(t_set.lt.0)THEN
  !   t_set=t_set+dpig
  !END IF
  ! convert times from radians to days
  t_rise=t_rise/dpig
  t_transit=t_transit/dpig
  t_set=t_set/dpig
END SUBROUTINE rise_and_set



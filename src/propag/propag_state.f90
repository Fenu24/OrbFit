!==========MODULE PROPAGATOR==============================

! PUBLIC ROUTINES  
!              pro_ele       propagator of elements
!              propag        general purpose propagator (to cartesian)  
!              inipro        initialises propagator options

! MODULE CONTAINS:                                                      
! ROUTINES         
!              sel_met                                        
!              propin  
!              bdnstev 
!              clocms 
!              bdintrp
!                compco_intrp                                                
!              catst                                                          
!              invaxv                                                 
!              inivar        
!              selste                                                   
!              compco                   
!              zed  
!              bessel  
!              clo_test  
! OUT OF MODULE
!              set_restart   to operate on restart control 
!              varwra  used by  tp_trace
!              varunw  used by  close_app (tp_fser)
!              vawrxv  used by  close_app (tp_fser) 
! 
!  HEADERS   
! propag_state.o: \
!	../include/nvarx.h90 \
!	../suit/FUND_CONST.mod \
!	../suit/OUTPUT_CONTROL.mod \
!	../suit/PLANET_MASSES.mod \
!	../suit/REFERENCE_SYSTEMS.mod \
!	close_app.o \
!	force_model.o \
!	propag_state.o \
!	ra15_mod.o \
!	runge_kutta_gauss.o \
!	tp_trace.o 

MODULE propag_state

USE fund_const
USE output_control
USE dyn_param
IMPLICIT NONE

PRIVATE

! public routines
PUBLIC pro_ele, propag, inipro, clo_test, intein, propin

! PUBLIC DATA
! to force restart anyway
LOGICAL, PUBLIC :: restar 

! fix to be able to force step submultiple of output interval
LOGICAL, PUBLIC :: fix_mul_step 
DOUBLE PRECISION, PUBLIC :: integration_step
    
! common to the module, but private

! former parint  header
INTEGER,PARAMETER :: mmax=20
INTEGER,PARAMETER :: m2max=2*mmax+2
INTEGER, PARAMETER :: itmaxx = 50 ! must be larger than all itmax's
! former comint  header
! controls for numerical integration: multistep only
DOUBLE PRECISION hms,deltos,error,epms,hmax_me
PUBLIC hms
INTEGER mms,imet,iusci,icmet
!   multistep coefficients: values computed from the beginning
DOUBLE PRECISION c_ms0(mmax+1),f_ms0(mmax+1),b_ms0(mmax+1),a_ms0(mmax+1)
!   multistep coefficients: duplicated array, for the current choice of mms
DOUBLE PRECISION c_ms(m2max),f_ms(m2max),b_ms(m2max),a_ms(m2max),cerr

CONTAINS
! ===================================================================== 
! PRO_ELE version with elements and optional covariance propagation
! ===================================================================== 
!                                                                       
!  input: el0 elements, including
!                  coordinate type EQU, KEP, CAR, COM                        
!                  epoch time (MJD) 
!                  orbital elements vector 
!         t1 prediction time (MJD)                                      
!         unc0  uncertainty of el0, including covariance and normal matrix
!  output:                                                              
!         el1 orbital elements for the asteroid at time t1            
!         unc1 uncertainty of el1, including  covariance and normal matrix
! ============INTERFACE=================================================
SUBROUTINE pro_ele(el0,t1,el1,unc0,unc1,obscod0,twobo,dede0) 
  USE orbit_elements
  USE close_app, ONLY: kill_propag
! equinoctal elements and epoch times
  TYPE(orbit_elem), INTENT(IN) :: el0  
  double precision, INTENT(IN) :: t1 
  TYPE(orbit_elem), INTENT(INOUT) :: el1
  INTEGER, INTENT(IN), OPTIONAL :: obscod0
  LOGICAl, INTENT(IN), OPTIONAL :: twobo
! normal, covarianvce matrices    
  TYPE(orb_uncert), INTENT(IN), OPTIONAL ::  unc0                       
  TYPE(orb_uncert), INTENT(INOUT), OPTIONAL ::  unc1
! state transition matrix
  DOUBLE PRECISION, INTENT(OUT),OPTIONAL :: dede0(6,6)
! ============END INTERFACE=============================================
! interface with propag
  double precision xastr(6),xear(6),dxdpar(6,ndimx)
  INTEGER ider 
! interface with coo_cha
  TYPE(orbit_elem) :: el2
  INTEGER fail_flag, obscode
! derivatives of cartesian w.r. elelments, new elements w.r. to cartesian,
!  new elements w.r. old
  double precision de1dx(6,6),de1de0(ndimx,ndimx)
  LOGICAL error, twobo1
  INTEGER j, nd ! for variable dimension covariance
! static memory not required                                          
! ===================================================================== 
  IF((PRESENT(unc0).and..not.PRESENT(unc1)).or.     &
&          (PRESENT(unc1).and..not.PRESENT(unc0)))THEN
     WRITE(*,*)' pro_ele: missing argument'
     STOP
  ELSEIF(PRESENT(unc0))THEN
     ider=1
  ELSE
     ider=0
  ENDIF    
  IF(PRESENT(obscod0))THEN
     obscode=obscod0
  ELSE
     obscode=el0%obscode
  ENDIF
  IF(PRESENT(twobo))THEN
     twobo1=twobo
  ELSE
     twobo1=.false.
  ENDIF
  IF(PRESENT(dede0)) ider=1
! ===================================================================== 
! call private propagation routine, requiring derivatives               
  IF(PRESENT(unc0))THEN
     nd=unc0%ndim
  ELSE
     nd=6+nls
  ENDIF
  CALL propag(el0,t1,xastr,xear,ider,nd,dxdpar,TWOBO=twobo1)
  IF(kill_propag)RETURN
! compute new elements, with partial derivatives  
  el2=el0
  el2%t=t1
  el2%coord=xastr
  el2%coo='CAR'
  IF(PRESENT(dede0))THEN
     CALL coo_cha(el2,el0%coo,el1,fail_flag,de1dx,OBSCODE=obscode)
     IF(fail_flag.ne.0)THEN
        WRITE(*,*)'pro_ele: coord ', el0%coo, ' fail_flag=',fail_flag
     ENDIF
! chain rule to obtain d(east1)/d(east0)                                
     dede0=MATMUL(de1dx,dxdpar(1:6,1:6)) 
  END IF
  IF(PRESENT(unc0))THEN
     CALL coo_cha(el2,el0%coo,el1,fail_flag,de1dx,OBSCODE=obscode)
     IF(fail_flag.ne.0)THEN
        WRITE(*,*)'pro_ele: coord ', el0%coo, ' fail_flag=',fail_flag
     ENDIF
! chain rule to obtain d(east1)/d(east0)
     de1de0(1:6,1:nd)=MATMUL(de1dx,dxdpar(1:6,1:nd))
     IF(nd.GT.6)THEN
        de1de0(7:nd,1:nd)=0.d0
        DO j=7,nd
           de1de0(j,j)=1.d0
        ENDDO
     ENDIF 
! covariance/normal matrix propagated by similarity transformation   
     CALL propagunc(nd,unc0,de1de0(1:nd,1:nd),unc1)
!     CALL convertunc(unc0,de1de0,unc1)
     IF(.not.unc1%succ)THEN 
        WRITE(*,*)'pro_ele: this should not happen' 
     ENDIF
  ELSE
     CALL coo_cha(el2,el0%coo,el1,fail_flag,OBSCODE=obscode)
     IF(fail_flag.ne.0)THEN
        WRITE(*,*)'pro_ele: coord ', el0%coo, ' fail_flag=',fail_flag
     ENDIF
  ENDIF       
END SUBROUTINE pro_ele
! ===================================================================== 
! PROPAG vers. 3.1
! A. Milani August 2003 - version with orbit_elem data type   
! ===================================================================== 
! N+1-body problem propagator with automatic restart control
!    if the orbital elements at epoch are unchanged, integration continues
!    from the state reached at thew previous call rather than restarting.
!    This has enormous efficiency gain if the observation times are passed
!    to the propagator in a suitable order, see least-squares.sort_obs 
!    Anyway the propagator can handle times passed in an arbitrary order,
!    and restart can be forced through the variable restar 
! 
! Uses JPL ephemerides as source for the planetary position, binary
!   ephemerides file generated by bineph for massive asteroids
! 
!  WARNING: the input elements and the output cartesian coordinates   
!           and derivatives are ecliptic (mean of J2000.0)
! ===================================================================== 
! 
!  input:        
!         t2 prediction time    
!         el orbital element record
!  output:  
!        xast position and velocity vector in heliocentric cartesian    
!               coordinates for the asteroid at time t2       
!        xea  position and velocity vector in heliocentric cartesian    
!               coordinates for the Earth at time t2
!        dxde first derivatives of position vector with respect to elements;
!             this is optional, variational equation not computed if 
!             this argument not present
!
! ================INTERFACE===========================================  
SUBROUTINE propag(el,t2,xast,xea,ider,nd,dxdpar,twobo)
  USE force_model
  USE force_sat
  USE orbit_elements
  USE fund_const
  USE ever_pitkin, ONLY: fser_propag,fser_propag_der
  USE close_app, ONLY: kill_propag
  USE planet_masses
! ===============INPUT==============================
! elements, including epoch, target epoch (MJD)      
  TYPE(orbit_elem),INTENT(IN) :: el
  DOUBLE PRECISION, INTENT(IN) :: t2
  INTEGER, INTENT(IN):: ider ! 1= derivatives required, 0=no
! OPTIONAL INPUT
  LOGICAL, INTENT(IN), OPTIONAL :: twobo
! ===============OUTPUT=============================
! cartesian coord., at time t2, of the asteroid, of the Earth, derivative
  DOUBLE PRECISION, INTENT(OUT) :: xast(6), xea(6)
  INTEGER, INTENT(IN) :: nd ! number of parameters in variational equation, .ge.6
  DOUBLE PRECISION, INTENT(OUT), OPTIONAL :: dxdpar(6,nd) ! d(pos,vel)/d(incond,param)
! ===========END INTERFACE===========================================   
  INCLUDE 'nvarx.h90' ! max dimensions for propagated state vector
! =============STATE VECTOR AND DERIVATIVES======================       
! main state vector for only one asteroid + variational eq.   
  double precision y1(nvarx),y2(nvarx) 
  TYPE(orbit_elem) elcar
! matrices of partial derivatives         
! derivatives of initial cartesian with respect to elements
! cartesian with respect to initial cartesian
  DOUBLE PRECISION :: dx0de(6,6), dxdxp0(6,ndimx)
! jacobian of all solve for parameters with respect to initial state, in elements and in cartesian
  DOUBLE PRECISION, DIMENSION(ndimx,ndimx):: dxp0dep0,dxpdxp0 
! other used in 2-body case
  DOUBLE PRECISION, DIMENSION(6,6) :: dxda,dxdx
  DOUBLE PRECISION det
  INTEGER ising
! =============RESTART CONTROL=========================       
! times: current  
  double precision t1 
! asteroid elements at the previous call  
  TYPE(orbit_elem) :: elsave
!asteroid equinoctal elements for 2-body propagation
  TYPE(orbit_elem) :: eleq
  INTEGER fail_flag ! for coo_cha call
  LOGICAL equal 
! restart control     
  integer nfl 
! close approach at initial time monitoring         
  INTEGER iplamloc 
!  DOUBLE PRECISION texit 
! ===================================================         
! stepsize: as given by selste, maximum alowed, current, previous  
  DOUBLE PRECISION :: hgiv,hmax,h,hs 
  DOUBLE PRECISION :: ecc,q,qg,enne ! eccentricity, perielion, aphelion, m.mot 
  DOUBLE PRECISION :: gmcur ! current central body mass (depends upon rhs)
! =============       
! integers for dimensions       
  INTEGER :: nv,nv1,nvar,nvar2 
! time step for attributables file is fixed submultiple  
  INTEGER :: nn 
! store previous time for Earth coordinates
  DOUBLE PRECISION :: t2old,xea1(6)
! initialisation control (to have dummy initial conditions), old ider
  INTEGER :: lflag,ider0
  LOGICAL :: twobo1
  DOUBLE PRECISION :: ddxde(3,6,6) ! 2nd derivatives of cart. w. r. elements
  INTEGER :: j ! loop index
! ************************************** 
! static memory allocation      
  save 
! initialization control    
  data lflag/0/ 
! **************************************  
!  dummy initial conditions to force restart the first time   
  IF(lflag.eq.0)THEN 
     lflag=1 
! Store fictitious epoch time and elements for the asteroid   
! (to be sure they are different the first time)    
     elsave=undefined_orbit_elem
     IF(rhs.EQ.2)THEN ! maybe not necessary
        elsave%center=3
     END IF
     t2old=-1d+55 
     t1=el%t
     ider0=-1
! masses, distance control and so on needed the first time    
     IF(rhs.eq.1)THEN
        CALL masjpl
     ELSEIF(rhs.eq.2)THEN
        CALL eamoon_mass
     ELSEIF(rhs.eq.3)THEN
! should have been done by input9
     END IF
  ENDIF
! ===================================================================== 
! JPL Earth vector at observation time   
  IF(t2.ne.t2old)THEN 
     IF(rhs.eq.1)THEN
        CALL earcar(t2,xea1,1)
     ELSEIF(rhs.eq.2)THEN
        xea1=0.d0
     ELSEIF(rhs.eq.3)THEN
        xea1=0.d0 ! not used 
     ENDIF
     t2old=t2
  ENDIF
  xea=xea1  
! ===================================================================== 
  IF(PRESENT(twobo))THEN
     twobo1=twobo
  ELSE
     twobo1=.false.
  ENDIF
  IF(twobo1)THEN
     IF(rhs.eq.1)THEN
        gmcur=gms
     ELSEIF(rhs.eq.2)THEN
        gmcur=gmearth
     ELSEIF(rhs.eq.3)THEN
        gmcur=gms ! not used
     ENDIF
! two body propagation  
     IF(el%coo.eq.'EQU')THEN
! use old style 2-body propagator with Kepler's equation
        CALL prop2b(el%t,el%coord,t2,xast,gmcur,ider,dxdpar,ddxde) 
        RETURN
     ELSEIF(el%coo.eq.'CAR')THEN
! ready for f-g series propagation
        IF(ider.eq.1) CALL eye(6,dxda)
     ELSE
        IF(nd.gt.6)THEN
           WRITE(*,*)' propag: non-grav with two body: ', nd,twobo1
           STOP
        ENDIF               
! convert to cartesian
        IF(ider.eq.1)THEN
           CALL coo_cha(el,'CAR',eleq,fail_flag,dxda) ! the matrix is d(car)/de 
           IF(fail_flag.ge.4)THEN
              WRITE(*,*)' propag: conversion to CAR for 2-body failed ',fail_flag,el
              STOP
           ENDIF
           CALL matin(dxda,det,6,0,6,ising,1)  ! converted to de/d(car)
        ELSE
           CALL coo_cha(el,'CAR',eleq,fail_flag)
           IF(fail_flag.ge.4)THEN
              WRITE(*,*)' propag: conversion to CAR for 2-body failed ',fail_flag,el
              STOP
           ENDIF
        ENDIF
     ENDIF
! f-g series propagation
     IF(ider.eq.0)THEN
        CALL fser_propag(el%coord(1:3),el%coord(4:6),el%t,t2,gmcur,xast(1:3),xast(4:6))
     ELSE
        CALL fser_propag_der(el%coord(1:3),el%coord(4:6),el%t,t2,gmcur,xast(1:3),xast(4:6),dxdx)
        dxdpar(1:6,1:6)=MATMUL(dxdx,dxda)
     ENDIF
     RETURN
  ENDIF
! ===================================================================== 
! Check if time and/or elements changed: we may not need to compute   
! initial conditions for the asteroid        
  equal=equal_orbels(el,elsave)    
  IF(.not.equal.or.abs(t2-t1).gt.abs(t2-el%t).or.ider.ne.ider0.or.restar)THEN  
!  new arc; propin needs to be informed   
     nfl=0         
! also cloapp needs to be informed        
     CALL set_clost(.true.)
! ===================================================================== 
! Compute asteroid cartesian elements at epoch time
     CALL coo_cha(el,'CAR',elcar,fail_flag,dx0de)
     dxp0dep0(1:6,1:6)=dx0de
     IF(nd.gt.6)THEN
        dxp0dep0(7:nd,1:nd)=0.d0
        dxp0dep0(1:6,7:nd)=0.d0
        DO j=7,nd
          dxp0dep0(j,j)=1.d0
        ENDDO     
     ENDIF
!  need new clotest using coo_cha
     CALL clo_test(elcar,iplamloc) 
     elsave=el 
     ider0=ider
! ===================================================================== 
! choices about integration method and stepsize
     IF(rhs.eq.1)THEN  
        CALL ecc_peri(el,ecc,q,qg,enne)
     ELSEIF(rhs.eq.2)THEN
        CALL ecc_peri(el,ecc,q,qg,enne) 
     ELSEIF(rhs.eq.3)THEN
! WARNING: automatic choice of method/dynamical model not supported in ORBIT9
        WRITE(*,*)'propag: ORBIT9 not supported here '
        STOP
     ENDIF
     if(imet.eq.0) then 
! automatical choice of numerical integration method 
! logic must be: if e<<1 we can use multistep, otherwise always ra15
! but we use q and Q and e and a in selmet, thus new selmet is needed
        IF(rhs.eq.1)THEN
           call sel_met(ecc,q,qg,hmax)
        ELSEIF(rhs.eq.2)THEN
           call sel_met_sat(ecc,q,qg,hmax)
        ELSEIF(rhs.eq.3)THEN
! WARNING: automatic choice of method/dynamical model not supported in ORBIT9
           WRITE(*,*)'propag: ORBIT9 not supported here '
           STOP
        END IF
     else 
        icmet=imet 
        icrel=irel
        iclun=ilun
     endif
! read masses from JPL header (stored in a module)         
! alignement of masses might have changed, because of different list of 
! asteroids (to avoid self perturbation) 
     IF(rhs.eq.1)THEN 
        CALL masjpl    
     ELSEIF(rhs.eq.2)THEN
        CALL eamoon_mass
     ELSEIF(rhs.eq.3)THEN
! WARNING: automatic choice of method/dynamical model not supported in ORBIT9
        WRITE(*,*)'propag: ORBIT9 not supported here '
        STOP 
     ENDIF
! control if both position and vel. are needed in dpleph      
     if(icmet.eq.3.and.iclap.eq.1) then 
        velo_req=.true. 
     elseif(icrel.ge.1)then 
        velo_req=.true.
     elseif(dyn%active(2))THEN
        velo_req=.true.
     else 
        velo_req=.false.
        if(icrel.ge.1)then 
           velo_req=.true.
        elseif(dyn%active(2))THEN
           velo_req=.true.
        endif
     endif
! *************************************** 
! selection of an appropriate stepsize    
     if(icmet.ne.3)then 
! for the multistep (the RK starter uses the same stepsize):
        call sel_ste(ecc,enne,error,mms,hmax,hgiv) 
! TO BE FIXED *********************
        IF(fix_mul_step)THEN
! then use the logic of former propag2, to force stepsize to be
! a sub-multiple of the interval between attributables rounded times
           If(hgiv.eq.integration_step)THEN
              h=hgiv
           ELSE
              nn=integration_step/hgiv+1 
              h=integration_step/nn
              hgiv=h
           ENDIF
        ELSE
           h=hgiv
        ENDIF
     elseif(icmet.eq.3)then 
! for Everhart: step selection is automatic, only sign matters        
        h=hms
     endif
! sign control:       
     if(t2-el%t.lt.0.d0)then 
        h=-abs(h) 
     endif
     hs=h
! ===================================================================== 
! position vector dimension     
     nv=3 
     if(ider.ge.1)then 
! Variational equations: 3 X 6 components 
        nv1=18+3*nls 
     else 
        nv1=0 
     endif
! total number of variables (and total number of positions)   
     nvar2=nv+nv1 
     nvar=nvar2*2 
! ===================================================================== 
! Vector y1, length 6, contains at least  
! the asteroid position and velocity      
     y1(1:3)=elcar%coord(1:3)
     y1(nvar2+1:nvar2+3)=elcar%coord(4:6) 
     if(ider.ge.1)then 
! in this case the vector y1 contain also the matrices        
! of partial derivatives; then we need to 
! initialise the  variation matrix as the 6 X 6 identity plus zeros
! for derivatives w.r. to dynamical parameters     
        CALL inivar(y1,nvar2,nvar) 
     endif
! ===================================================================== 
! Restart from epoch  
     t1=el%t 
  else 
! old arc:  
! stepsize was already decided; do not change!!!!
!     if(icmet.eq.3)then 
!        h=hms 
!     else 
!        h=hgiv 
!     endif
! sign control: 
     IF(t2-t1.ne.0.d0)THEN 
        if(t2-t1.lt.0.d0)then 
           h=-abs(h) 
        elseif(t2-t1.gt.0.d0)then
           h=abs(h)
        endif
! need to check that direction is not changed       
        if(h*hs.lt.0.d0)then 
! restart is necessary because of U-turn 
           IF(nfl.eq.1)THEN
              WRITE(*,*)'propag: reverse ',h,hs,t1,t2,t2-t1,nfl 
              nfl=0
           ENDIF
        else 
! restart not necessary; propin can go on 
           IF(t1.ne.el%t)THEN
              nfl=1 
           ELSE
              nfl=0
           ENDIF
        endif
        hs=h
     ENDIF
  endif
! ===================================================================== 
! propagator to compute asteroid and planets orbits
  IF(t1.ne.t2)THEN 
     CALL propin(nfl,y1,t1,t2,y2,h,nvar,nd,dxp0dep0(1:nd,1:nd))
  ELSE
     y2=y1
  ENDIF
  IF(kill_propag)RETURN
!================================================================== 
! Asteroid coordinates
  xast(1:3)=y2(1:3) 
  xast(4:6)=y2(nvar2+1:nvar2+3) 
  if(ider.ge.1)then 
! if partial derivative are required      
! rewrap vector into (6 x nd) matrix    
     CALL varwra(y2,dxdxp0(1:6,1:nd),nd,nvar,nvar2) 
! ===================================================================== 
! Chain rule: we compute dxdpar=dxdxp0*dxp0dep0 6 x nd=(6 x nd) (nd x nd) 
     dxdpar=MATMUL(dxdxp0(1:6,1:nd),dxp0dep0(1:nd,1:nd))
  endif
  return 
END SUBROUTINE propag         

! ================================================================      
! INIPRO    
! ================================================================      
! this subroutine reads from file 'propag.def' all the propagator option
!  WARNING: we advise the user against changing the file propag.def;    
!                do it at your risk !.... 
!  options: 
!  control of propagation methods:                  
!           imet =1 (multistep) =2 (runge-kutta) =3 (everhart)
!             =0 automatic (multistep for main belt, Everhart for high  
!                 eccentricity and/or planet crossing)        
!    Runge-Kutta-Gauss: (imet=2; also used as starter for imet=1)       
!           isrk = order of the method is isrk*2    
!           h = integration step
!           eprk = convergence control in the solution of the implicit e
!           lit1,lit2 = no. of Gauss-Seidel iterations (first step, afte
! 
!    multistep:       
!           mms = multistep number of previous steps
!         WARNING: order is mms+2; e.g in Milani and Nobili, 1988, m=mms
! 
!    everhart:        
!           h = stepsize (fixed if  llev.le.0; initial if llev.gt.0))   
!           llev = control is 10**(-llev) 
! 
!           iusci = control of output numerical parameters    
! ================================================================      
SUBROUTINE inipro 
! no interface    
  USE runge_kutta_gauss
  USE ra15_mod
  logical fail,fail1,found 
  integer iork,iork_c,iord,isfl 
  integer j 
!****************     
!   static memory not required (used only once)     
!****************     
  fail=.false. 
! choice of propagator method (0=auto)         
  call rdnint('propag.','imet',imet,.true.,found,fail1,fail) 
! accuracy of time 
  call rdnrea('propag.','deltos',deltos,.true.,found,fail1,fail)
! use of fixed multiple of stepsize: default no
  fix_mul_step=.false.
! options for multistep
  call rdnrea('propag.','error',error,.true.,found,fail1,fail) 
  call rdnint('propag.','iord',iord,.true.,found,fail1,fail) 
  mms=iord-2 ! iord= multistep order; mms=no. back steps required to start
  call rdnrea('propag.','hms',hms,.true.,found,fail1,fail) 
  hms=abs(hms) ! max stepsize (only for RKG/multistep)
  call rdnrea('propag.','hmax_me',hmax_me,.true.,found,fail1,fail) 
  hmax_me=abs(hmax_me)   ! max stepsize if there is Mercury (for RKG/multistep)
  call rdnrea('propag.','epms',epms,.true.,found,fail1,fail)
! note epms can be negative/zero, meaning no corrector 
! options for RKGauss  
  call rdnint('propag.','iork',iork,.true.,found,fail1,fail) 
  isrk=iork/2  !number of intermediate steps
  call rdnrea('propag.','eprk',eprk,.true.,found,fail1,fail) 
  call rdnint('propag.','lit1',lit1,.true.,found,fail1,fail) 
  call rdnint('propag.','lit2',lit2,.true.,found,fail1,fail) 
! options for Radau      
  call rdnint('propag.','llev',llev,.true.,found,fail1,fail) 
  call rdnrea('propag.','hev',hev,.true.,found,fail1,fail) 
  hev=abs(hev) ! this is the initial stepsize, then it is automatic
  call rdnrea('propag.','eprk_r',eprk_r,.true.,found,fail1,fail)
  call rdnint('propag.','lit1_r',lit1_r,.true.,found,fail1,fail) 
  call rdnint('propag.','lit2_r',lit2_r,.true.,found,fail1,fail)
! options for Radau when used with fixed stepsize
  call rdnint('propag.','lit1_rc',lit1_rc,.true.,found,fail1,fail) 
  call rdnint('propag.','lit2_rc',lit2_rc,.true.,found,fail1,fail)
  deriv_safe=.true.
! options for  RKGauss during close approach 
  iork_c=iork   
  isrk_c=iork_c/2 !number of intermediate steps
  call rdnrea('propag.','eprk_c',eprk_c,.true.,found,fail1,fail) 
  call rdnint('propag.','lit1_c',lit1_c,.true.,found,fail1,fail) 
! verbosity options, obsolescent        
  call rdnint('propag.','iusci',iusci,.true.,found,fail1,fail) 
  if(fail)stop '**** inipro: abnormal end ****' 
!   initialize multistep and implicit Runge-Kutta-Gauss       
! ********************************************************************  
!  input coefficients for Runge-Kutta     
!  if(imet.ne.3)then 
     call legnum(isrk,isfl) 
     if(isfl.ne.0)then 
        write(*,997)isrk,isfl 
997     format(' required isrk=',i4,'  found only up to ',i4) 
        stop 
     endif
 ! endif
! ********************************************************************  
!   calcolo coefficienti multistep        
!   cambio notazione-nell'input iord=m nell'art. cel.mech.    
!   d'ora in poi mms come in revtst, orbit8a        
  if(mms.gt.mmax)then 
     write(ipirip,998)mms,mmax 
998  format(' chiesto mms=',i4,'  spazio solo per ',i4) 
     stop  
  endif
!   calcolo coefficienti predittore       
!   c=Cow pred f=Cow corr b=Ad pred a=Ad corr       
!   warning: one order more than used, because c(m+2) is required       
!   by the error formula        
  call compco(mms+2,c_ms0,f_ms0,b_ms0,a_ms0) 
  cerr=c_ms0(mms+2) 
!  duplicazione per ridurre i calcoli di indici nel multistep 
! this is useful only if sel_ste is not used
  do  j=1,mms+1 
     c_ms(j)=c_ms0(j) 
     f_ms(j)=f_ms0(j)
     a_ms(j)=a_ms0(j) 
     b_ms(j)=b_ms0(j)
     c_ms(j+mms+1)=c_ms0(j) 
     f_ms(j+mms+1)=f_ms0(j)
     a_ms(j+mms+1)=a_ms0(j) 
     b_ms(j+mms+1)=b_ms0(j)    
  enddo
END SUBROUTINE inipro
! *****************************************************                 
!  INTEIN initialisation of the  propagator                             
!  version for propin9c, with automatic stepsize control                
! *****************************************************                 
SUBROUTINE intein(iun,dtin,aa,enne,ecc,norb,norbx)
  USE runge_kutta_gauss
  USE ra15_mod
  INTEGER, INTENT(IN) :: iun, norb, norbx ! option unit, no. bodies, max
! filter input delta t, mean motion, semimajor axis, eccentricity 
  DOUBLE PRECISION, INTENT(IN):: dtin,enne(norbx),aa(norbx),ecc(norbx)
! END INTERFACE 
  double precision step(norbx),z(norbx),err(norbx) 
! output units 
! include 'proout.h' 
  character*1 cc 
!                                                                  
  integer iauto,nn,j,nb,iork,isfl,iord,m 
  double precision hh,hhh,emin,econv,eps,sign,ee, hmax 
  double precision roff 
  INTEGER icha ! now dummy
! *****************************************************                 
!  options for numerical integration: input from options file           
  read(iun,108)cc 
108 format(a1) 
! stepsize (either fixed or maximum)                                    
  call reaflo(iun,'h',hmax) 
!  h must go in the right direction - not here                          
  if(hmax*dtin.lt.0.d0)hmax=-hmax 
  hh=abs(hmax) 
  sign=hmax/hh 
!  precision for output times                                           
  deltos=abs(hmax)/1.d5 
!  options for automatic stepsize control                               
  call reaint(iun,'iauto',iauto) 
  call reaflo(iun,'error',error) 
!  options:                                                             
!  controllo metodi propagazione:                                       
!           imet=1 (multistep) =2 (runge-kutta) =3 (everhart)           
!    runge-kutta:                                                       
!           isrk=ordine del runge-kutta (diviso 2)                      
!           h=passo di integrazione                                     
!           eprk=controllo di convergenza nell'equaz. implicita del rk  
!           lit1,lit2=iterazioni di gauss-seidel al primo passo e dopo  
!           iusci=flag uscita controlli numerici                        
!    multistep:                                                         
!           mms=ordine del multistep ('m-2' nell'art. Cel.Mech)         
! nnnn      ipc=iteraz correttore; epms=controllo;                      
!           iusci=flag uscita controlli numerici                        
!    everhart:                                                          
! nnnn      isrk=ordine (per ora solo 15)                               
!           h=passo (fisso se ll.le.0)                                  
!           ll=controllo e' 10**(-ll); se ll.gt.0, scelta automatica del
! nnnn      iusci=flag uscita contr. num.; se iusci.ge.0, uscita passo c
  call reaint(iun,'iord',iord) 
  mms=iord-2 
  call reaflo(iun,'epms',epms) 
  call reaint(iun,'iork',iork) 
  isrk=iork/2 
  call reaflo(iun,'eprk',eprk) 
  call reaint(iun,'lit1',lit1) 
  call reaint(iun,'lit2',lit2) 
  call reaint(iun,'imet',imet) 
  IF(imet.le.0.or.imet.gt.3)THEN
     WRITE(*,*)' intein: option not supported imet=', imet
     STOP
  ENDIF
  icmet=imet
!  IF(imet.eq.3)THEN
     hev=hmax
     lit1_r=lit1
     lit2_r=lit2
     eprk_r=eprk
!  ENDIF
  call reaint(iun,'iusci',iusci) 
  call reaint(iun,'icha',icha) 
  call reaint(iun,'ll',llev) 
! ********************************************************************  
!  find machine error; if the machine has statistically even (IEEE 754) 
!  rounding off, the rounding off error will be of the order of eps     
!  per revolution squared, or even less; otherwise, the rounding off    
!  will be the dominant part of the error                               
  eps=EPSILON(1.d0) 
! ********************************************************************  
!  input coefficients for Runge-Kutta                                   
  if(imet.ne.3)then 
     call legnum(isrk,isfl) 
     if(isfl.ne.0)then 
        write(9,997)isrk,isfl 
997     format(' required isrk=',i4,'  found only up to ',i4) 
        stop 
     endif
  endif
! ********************************************************************  
!   calcolo coefficienti multistep                                      
!   cambio notazione-nell'input iord=m nell'art. cel.mech.              
!   d'ora in poi m come in revtst, orbit8a                              
!     write(9,*)iord                                                    
  if(imet.eq.1)then 
     m=mms 
     if(m.gt.mmax)then 
        write(9,998)m,mmax 
998     format(' required m=',i4,' space only for ',i4) 
        stop 998 
     endif
!   calcolo coefficienti predittore                                     
!   c=Cow pred f=Cow corr b=Ad pred a=Ad corr                           
!   warning: one order more than used, because c(m+2) is required       
!   by the error formula                                                
     call compco(m+2,c_ms0,f_ms0,b_ms0,a_ms0)
 !  duplicazione per ridurre i calcoli di indici nel multistep 
! this is useful only if sel_ste is not used
     do  j=1,m+1 
        c_ms(j)=c_ms0(j) 
        f_ms(j)=f_ms0(j)
        a_ms(j)=a_ms0(j) 
        b_ms(j)=b_ms0(j)
        c_ms(j+mms+1)=c_ms0(j) 
        f_ms(j+mms+1)=f_ms0(j)
        a_ms(j+mms+1)=a_ms0(j) 
        b_ms(j+mms+1)=b_ms0(j)    
     enddo
     cerr=c_ms0(m+2) 
! ***************************************************************       
!  file header for                                                      
!  numeric controls file (unit 10)                                      
     ipirip=10 
     open(ipirip,file='orb9.num',status='unknown') 
     if(iusci.gt.0)then 
        write(ipirip,113)hmax,iord,epms,iork,eprk,                   &
     &           lit1,lit2,imet,iusci,icha                              
113     format(' hmax=',f10.5,' iord=',i3,' epms=',1P,d12.4,0p,'iork=',i3,' eprk=',1P,d12.4/5i4)
     endif
! ********************************************************************  
!  for each orbit, compute the truncation error in longitude            
!  according to Milani and Nobili, 1988, Celest. Mech. 43, 1--34        
!  and select a stepsize giving a truncation  error equivalent to       
!  one roundoff per revolution squared (or to the value given in input) 
! (the maximum is however held fixed at hh)                             
     emin=0.05 
     econv=1.d0 
     write(9,*)' Selection of the best stepsize, maximum=',hmax
     DO 10 j=1,norb 
        ee=max(ecc(j),emin) 
        call sel_ste(ee,enne(j),error,mms,hh,step(j)) 
        write(9,109)j,step(j),aa(j),ecc(j) 
109     format(' orbit no.',i3,' sugg.step=',f7.4,' a,e=',f6.3,f8.4) 
10   ENDDO
! ********************************************************************  
!  now find the minimum for all orbits                                  
     hhh=hh  ! =abs(hmax)
     DO j=1,norb 
        hhh=min(hhh,step(j)) 
     ENDDO
     if(iauto.gt.0)then 
        hms=hhh*sign 
        nn=dtin/hms+1 
        hms=dtin/nn 
        write(9,110)hhh,hms,nn 
110     format(' min. step required ',f9.6,' selected ',f9.6,' nn=',i6) 
     else 
        write(9,111)hms 
111     format(' fixed stepsize h=',d12.5) 
     endif
  endif
! ********************************************************************  
END SUBROUTINE intein



! version 3.6.1 A. Milani D. Farnocchia Oct 2008  
! =====================================================================      
! *** SEL_MET_SAT ***
! Choice of numerical integration method for satellite  
! 
! This routine is called from propag only if  imet=0
! in 'namerun'.top.   
!            INPUT: eccentricity, pericenter, apocenter 
!            OUTPUT : icmet etc. stored in model.h  
! ===================================================================== 
SUBROUTINE sel_met_sat(ecc,q,qg,hmax)
  USE force_model
  USE planet_masses
! satellite  eccentricity, perielion, aphelion 
  DOUBLE PRECISION, INTENT(IN) :: ecc,q,qg 
  DOUBLE PRECISION, INTENT(OUT) :: hmax
  hmax=1.d0 
!  IF(ecc.gt.0.3d0) icmet=3
END SUBROUTINE sel_met_sat

!=====================================================================!
! SEL_MET                                                             !
!=====================================================================!
! Choice of numerical integration method and of the dynamical model.  !
!---------------------------------------------------------------------!
! Version 3.1: A. Milani, Aug. 2003                                   !
!=====================================================================!
SUBROUTINE sel_met(ecc,q,qg,hmax) 
  USE runge_kutta_gauss
  USE ra15_mod
  USE force_model
  USE planet_masses
  !==============================================================================================
  DOUBLE PRECISION, INTENT(IN)  :: ecc  ! Eccentricity
  DOUBLE PRECISION, INTENT(IN)  :: q    ! Perihelion
  DOUBLE PRECISION, INTENT(IN)  :: qg   ! Aphelion
  DOUBLE PRECISION, INTENT(OUT) :: hmax ! Maximum stepsize
  !==============================================================================================
  INTEGER          :: iord   ! Multistep order
  INTEGER          :: iork   ! Runge-Kutta order
  DOUBLE PRECISION :: qmin   ! Minimum value for the perihelion
  DOUBLE PRECISION :: qgmax  ! Maximum value for the aphelion
  DOUBLE PRECISION :: eccmax ! Maximum for the eccentricity
  INTEGER          :: j      ! Loop index
  !==============================================================================================
  !----------------!
  ! MBA definition !
  !----------------!
  qmin   = 1.4d0  ! All NEO
  qgmax  = 4.3d0  ! Almost Jupiter crossing
  eccmax = 0.25d0 ! For the multistep activation
  !--------------------------!
  ! Storage from the options !
  !--------------------------!
  icrel=irel ! Relativity flag
  iclun=ilun ! Lunar flag

  !****************************************************************!
  ! SELECTION OF THE DYNAMICAL MODEL AND OF THE INTEGRATION METHOD !
  !****************************************************************!
  IF(radar)THEN
     !----------------------------------------------!
     ! Radar data: same dyn. model and Radau anyway !
     !----------------------------------------------!
     icrel = 2 
     iclun = 1 
     imerc = 1 
     llev  = MAX(12,llev) 
     icmet = 3 
     !iast   = MAX(iast,3) 
     !iatrue = iast
  ELSEIF(qg.LT.5.d0)THEN 
     !------------------------!
     ! Asteroids (NEA or MBA) !
     !------------------------!
     iplut = 0 ! Pluto is irrelevant anyway
     imerc = 1 ! Mercury is required
     IF(q.GT.qmin .AND. qg.LT.qgmax .AND. ecc.LT.eccmax)THEN 
        !----------------!
        ! Main Belt case !
        !----------------!
        icmet = 1 ! Multistep
        iord  = 8 
        mms   = iord-2 
        iork  = 8 
        isrk  = iork/2 
     ELSE 
        !-----------------------------------------------------------------!
        ! Non MBA: NEA, Mars crosser, high eccentricity, Jupiter crossing !
        !-----------------------------------------------------------------!
        IF(q.LT.qmin)THEN
           icrel = 2 ! Before it was MAX(irel,1), now always set to 2, ADV 25/01/2017 
           iclun = 1 ! Moon as separate body
        END IF
        icmet=3 ! Radau method
     ENDIF
  ELSEIF(q.LT.qmin)THEN
     !-------------------------------------------!
     ! Short periodic comet/comet like asteroids !
     !-------------------------------------------!
     icrel = 1 
     iclun = 1 
     icmet = 3 
  ELSEIF(qg.GT.32.d0 .AND. q.GT.20.d0 .AND. ecc.LT.eccmax)THEN 
     !-----!
     ! EKO !
     !-----!
     iplut = 0 
     icrel = 1 
     imerc = 1 
     ! Multistep integration method  
     icmet = 1
     iord  = 8 
     mms   = iord-2 
     iork  = 8 
     isrk  = iork/2 
     ! Trojans
  ELSEIF(q.GT.4.d0 .AND. qg.LT.7.d0)THEN
     !---------!
     ! Trojans !
     !---------!
     iplut = 0 
     icrel = 1 
     imerc = 1 
     ! Multistep integration method
     icmet = 1 
     iord  = 8 
     mms   = iord-2 
     iork  = 8 
     isrk  = iork/2 
  ELSE
     !----------------------------------------------!
     ! Everything else: Centaurus, comets (Hildas?) !
     !----------------------------------------------!
     icmet = 3 
     icrel = 1 
     iplut = 0
     imerc = 1 
  ENDIF
  !---------------------------------------------------------!
  ! Control on stepsize to avoid instability due to Mercury !
  !---------------------------------------------------------!
  IF(imerc.EQ.1)THEN
     hmax = hmax_me
  ELSE
     hmax = hms
  ENDIF
  DO  j=1,mms+1 
     c_ms(j)       = c_ms0(j) 
     f_ms(j)       = f_ms0(j)
     a_ms(j)       = a_ms0(j) 
     b_ms(j)       = b_ms0(j)
     c_ms(j+mms+1) = c_ms0(j) 
     f_ms(j+mms+1) = f_ms0(j)
     a_ms(j+mms+1) = a_ms0(j) 
     b_ms(j+mms+1) = b_ms0(j)    
  ENDDO
  ! Write the selected/forced option
  IF(verb_pro.GT.10)THEN 
     WRITE(ipirip,*) 'e,q,Q = ',ecc,q,qg,' icmet =',icmet 
     WRITE(ipirip,*) 'hms, epms, eprk, deltos, error, mms, isrk, ',&
          & 'lit1, lit2, iusci, llev, hev'       
     WRITE(ipirip,*) hms,epms,eprk,deltos,error,mms,isrk, &
          & lit1,lit2,iusci,llev,hev        
     WRITE(ipirip,*) 'Forced: ilun, imerc, iplut, irel, iast, iaber, istat' 
     WRITE(ipirip,*) ilun,imerc,iplut,irel,iast,iaber,istat 
  ENDIF
END SUBROUTINE sel_met
   
! **********************************************************  
!  PROPIN   ORBFIT version, with multistep interpolator, July 2003 
!  purpose: propagator/interpolator; the state vector y1      
!           at time t1 is propagated to y2 at time t2         
!  input:   
!      nfl: if 0, integration needs to be restarted 
!      t1: initial time         
!      y1: state vector at t1   
!      nfl: flag $>0$ if continuing, $=0$ if it has to restart
!      t2: final time 
!      h: stepsize (only if it is fixed; ra15 uses hev as initial step)  
!      nvar: length(y1)     !      nvar2: nvar/2     
!      dx0pdep: derivatives of initial cartesian with respect to elements,
!               including dynamical parameters  
!  output:  
!      t2: time to which y2 refers (different by .lt.deltos from input) 
!      y2: state vector at time t2        
!      t1,y1: time and state vector after last integration step         
! 
!  osservazione: se nelle chiamate successive alla prima t1 ed y1 non   
!                sono stati modificati dal programma chiamante,l'integra
!                zione prosegue a partire da t1,y1 usando i dati dei pas
!                si precedenti accumulati nelle cataste ck (per rkimp)  
!                e dd, delta (per il multistep)     
!                se invece t1,y1 sono stati modificati, il che risulta  
!                da nfl=0,l'integrazione riprende con un passo iniziale.
!                se pero' t1,y1 non sono modificati,i tempi t2 devono es
!                sere in successione monotona consistente con h         
!  subroutine chiamate:         
!       ra15: metodo di everhart
!       force: external per chiamata secondo membro 
!       legnum :lettura coefficienti r-k  
!       rkimp :propagatore runge-kutta    
!       kintrp :interpolatore di ck       
!       bdnste :multistep BDS che fa nstep passi
!       bdinterp: interpolatore multistep    
!       catst :caricatore catasta differenze e somme
!       coeff :calcolo coefficienti multistep       
!       rkstep :variatore automatico del passo      
!       fct :secondo membro equazioni differenziali (fisso; riduce      
!  al primo ordine quanto fornito da force)         
!   codice indici : i=1,nvar, j=1,isrk    
! **********************************************************  
SUBROUTINE propin(nfl,y1,t1,t2,y2,h,nvar,nd,dxp0dep0) 
  USE runge_kutta_gauss
  USE ra15_mod
  USE tp_trace, ONLY: str_clan
  USE surf_trace, ONLY: str_surfan
  USE force_model
  USE planet_masses
  USE force_sat
  USE force9d
  USE close_app, ONLY: kill_propag,min_dist,surf_stop
! INPUT
  integer, intent(in) ::  nfl ! restart control 
  integer, intent(in) ::  nvar! actual dimension of y1,y2
  integer, intent(in) ::  nd  ! number of parameters being solved
  double precision, intent(inout) :: h ! times, stepsize
! derivatives of state vector with respect to initial state vector,
! including dynamical parameters (if any)
  DOUBLE PRECISION, INTENT(IN), OPTIONAL :: dxp0dep0(nd,nd)
! INPUT/OUTPUT
  INCLUDE 'nvarx.h90'
  double precision, intent(in) :: t2 ! final time
  double precision,intent(inout) :: t1,y1(nvarx) ! current state vector at t1
! OUTPUT
  double precision, intent(out) :: y2(nvarx) ! state vector at t2
! END INTERFACE
!  workspace per rkimp
  double precision ck(ismax,nvarx),ck1(ismax,nvarx),ck2(ismax,nvarx) 
  double precision ep(itmaxx),dery(nvarx) 
!  workspace for bdstep         
  double precision dd(nvar2x,4),delta(nvar2x,m2max)
  double precision :: y3(nvarx) ! test state vector at t2
!  close approaching planet     
  double precision xxpla(6) 
! internal times, stepsize squared   
  double precision t0,tint,din,h2
! controls  
  double precision sdelto,epin 
  integer ndim,nvar2 ! dimensions
! counters  
  integer nstep,nrk,npas 
! integers  
  integer lflag, idc,j,j1,n,lf,it,lit,i,m 
  integer nclass ! type of equation for ra15
! arrays to store state transition matrix (initial conditions only)
  DOUBLE PRECISION stm(6,nd),stmout(6,nd),stm0(6,nd),tcur
! test for need of velocity at each step
!  LOGICAL, EXTERNAL ::  velocity_req 
! name of right hand side routine         
!      external force 
! **************************************  
! static memory allocation      
  save 
! **********************************************************  
!  controllo metodi propagazione:         
!           imet=1 (multistep) =2 (runge-kutta) =3 (everhart) 
!    runge-kutta:     
!           isrk=ordine del runge-kutta (diviso 2)            
!           h=passo di integrazione       
!           eprk=controllo di convergenza nell'equaz. implicita del rk  
!           lit1,lit2=iterazioni di gauss-seidel al primo passo e dopo  
!           iusci=flag uscita controlli numerici    
!    multistep:       
!           mms=ordine del multistep ('m' nell'art. Cel.Mech) 
! nnnn      ipc=iteraz correttore; epms=controllo   
! nnnn      iusms=flag uscita controlli numerici    
!    everhart:        
! nnnn      isrk=ordine (per ora solo 15) 
!           h=passo (fisso se ll.le.0)    
!           ll=controllo 10**(-ll); se ll.gt.0,     
!    scelta automatica del passo
!           iusci=flag uscita contr. num.;
!                  se iusci.ge.0, uscita passo cambiato       
  data lflag/0 / 
! *******************************************************   
!  dimensione del sec membro eq ordine 2  
  nvar2=nvar/2 
  IF(abs(t1-t2).lt.deltos)THEN 
     y2(1:nvar)=y1(1:nvar) 
     RETURN 
  ENDIF
! Radau 15
  if(icmet.eq.3)then 
!  everhart method (propagates to exact time):      
     if(velo_req)then 
        nclass=2 
     else 
        nclass=-2 
     endif
! This IF separates the case of OrbFit, call from propag_state.mod, from
! The one of orbit9, call fromn orbit9.f90 and similar
     IF(PRESENT(dxp0dep0))THEN
! accumulate state transition matrix (for initial conditions, but also
! for dynamical parameters if any)
        IF(nvar.gt.6)THEN 
           CALL vawrxv(y1,y1(nvar2+1),stm0,nvar2,nd) 
           CALL invaxv(y1,y1(nvar2+1),nvar2) 
        ENDIF
     ENDIF
666  CONTINUE 
     CALL ra15(y1,y1(nvar2+1),t1,t2,tcur,nvar2,nclass,idc) 
     IF(tcur.eq.t2)THEN 
!  current time and state is the final one;         
        IF(PRESENT(dxp0dep0))THEN
           IF(nvar.gt.6)THEN 
!  but the accumulated state    
!  transition matrix stm0 has to be used:     
              CALL vawrxv(y1,y1(nvar2+1),stm,nvar2,nd) 
! for derivatives w.r. to intiial conditions
              stmout(1:6,1:6)=MATMUL(stm(1:6,1:6),stm0(1:6,1:6))
! for derivatives w.r. to parameters
              IF(nd.gt.6)THEN
                 DO j=7,nd
                    stmout(1:6,j)=MATMUL(stm(1:6,1:6),stm0(1:6,j))+stm(1:6,j)
                 ENDDO
              ENDIF
! the partials with respect to the dyn.parameters
! are already renormalized; they are left as they are
              CALL varunw(stmout,y1,y1(nvar2+1),nd,nvar2) 
           ENDIF
           ! If the propagation ends but a close approach is still
           ! going on, the output mtp/tp trace must be available anyway
           IF(min_dist .AND. idc.NE.0) CALL str_clan(stm0,nd,dxp0dep0)
        ENDIF
        y2(1:nvar)=y1(1:nvar) 
        t1=t2 
! the propagation with Radau15 is terminated
        RETURN
     ELSE 
! write message      
!       WRITE(*,*)'end close approach to planet',idc,tcur       
! propagation has been interrupted because of a close approach
        IF(PRESENT(dxp0dep0))THEN
           IF(nvar.gt.6)THEN 
! setup the close approach record with derivatives  
              IF(min_dist) CALL str_clan(stm0,nd,dxp0dep0)
! Store information for IC center
              IF(surf_stop.AND.kill_propag) CALL str_surfan(stm0,nd,dxp0dep0) 
! multiply the accumulated state transition matrix          
              CALL vawrxv(y1,y1(nvar2+1),stm,nvar2,nd) 
! for derivatives w.r. to intiial conditions
              stmout(1:6,1:6)=MATMUL(stm(1:6,1:6),stm0(1:6,1:6))
! for derivatives w.r. to parameters
              IF(nd.gt.6)THEN
                 DO j=7,nd
                    stmout(1:6,j)=MATMUL(stm(1:6,1:6),stm0(1:6,j))+stm(1:6,j)
                 ENDDO
              ENDIF
              stm0=stmout
              CALL invaxv(y1,y1(nvar2+1),nvar2)
           ELSE
              ! setup the close approach record without derivatives
              IF(min_dist) CALL str_clan(stm0,nd,dxp0dep0)
           ENDIF
           IF(kill_propag)THEN
              IF(verb_pro.gt.9.AND.(.NOT.surf_stop)) WRITE(*,*)' propin: returning because of kill_propag, t=',tcur 
              RETURN
           ENDIF
        ENDIF
        t1=tcur 
        GOTO 666 
     ENDIF
  endif
! **********************************************************
! set m (order of multistep method) 
  m=mms 
  if(lflag.eq.0)then 
!  inizializzazione     
!  contapassi:npas caricamento delta
!  nrk caricamento ck   
     npas=0 
     t0=t1 
     nrk=0 
!  dimensione usata per controllo di convergenza (no variational eq.)   
     ndim=3 
!  inizializzazioni per il multistep
     j1=0 
     h2=h*h 
!  fine inizializzazioni
     lflag=1 
  else 
!  controllo se nuovo arco          
     if(nfl.eq.0)then 
        npas=0 
        t0=t1 
        nrk=0 
!  inizializzazioni per il multistep
        j1=0 
        h2=h*h 
     endif
  endif
!******************************************************     
!  controllo direzione del tempo    
  if((t2-t1)*h.lt.0.d0) then 
     write(*,999)t1,t2,h 
999  format('propin: from t1=',f12.4,' to t2=',f12.4,  &
     &         ' with step h=',d12.4)           
     write(*,*)'propin: this should not happen' 
     STOP 
!          h=-h         
  endif
! **********************************************************
!  main loop:           
!  controllo di funzione: propagatore o interpolatore       
5 if(dabs(t2-t1).lt.dabs(h)-deltos)goto 90 
!  scelta propagatore   
  if(npas.lt.m.or.icmet.ge.2)then 
! **********************************************************
!  runge kutta implicito fa un solo passo       
24   continue 
!  passo iniziale?      
     IF(nrk.le.0)THEN 
!  passo iniziale: store secondo membro  
        IF(rhs.eq.1)THEN       
           CALL force(y1,y1(nvar2+1),t1,delta(1,m+1),nvar2,idc,xxpla,0,1)
! close approach control
           CALL clocms(idc,t1,y1(1:3),y1(nvar2+1:nvar2+3),xxpla) 
        ELSEIF(rhs.eq.2)THEN
           CALL forcesat(y1,y1(nvar2+1),t1,delta(1,m+1),nvar2,idc,xxpla,0,1)
        ELSEIF(rhs.eq.3)THEN
           CALL force9(y1,y1(nvar2+1),t1,delta(1,m+1),nvar2,idc,xxpla,0,1) 
           CALL clocms9(idc,t1)
        ENDIF
!  passo iniziale: inizializzazione di ck a zero
        ck=0.d0 
        lit=lit1 
     ELSE 
!  passo non iniziale : interpolazione dei ck   
        ck1(1:isrk,1:nvar)=ck(1:isrk,1:nvar) 
        call kintrp(ck1,ck,isrk,nvar) 
        lit=lit2 
     ENDIF
!  un passo del rk      
22   call rkimp(t1,h,y1,dery,ck,isrk,y2,lit,nvar,eprk,ep,lf,ndim)
!     WRITE(*,*)t1,h,npas,nrk,m
!  controllo di avvenuta convergenza
     if(lf.le.0)then 
!  caso di non convergenza          
!           call camrk(ep,npas,nrk,lf)
! WARNING: quick fix for strange orbits resulting in RKG divergence  
        IF(rhs.eq.2)THEN      
           kill_propag=.true.
           RETURN
        ENDIF
        CALL rkstep(ep,npas,nrk,lf,h) 
        h2=h*h 
        if(lf.eq.2)goto 24 
        goto 22 
     endif
!  passo del rk adottato
     npas=npas+1 
     nrk=nrk+1 
     t1=t0+h*npas 
     y1(1:nvar)=y2(1:nvar) 
     if(iusci.gt.0)then 
        it=iabs(lf) 
        if(npas.eq.iusci*(npas/iusci))write(ipirip,1000)npas,       &
             &    (ep(i),i=1,it)
1000    format(' npas',i6,' ep ',5d12.3/(5d12.3)) 
     endif
     if(icmet.lt.2)then 
!  preparazione per il multistep   
        IF(rhs.eq.1)THEN      
           CALL force(y1,y1(nvar2+1),t1,delta(1,m+1-npas),nvar2,idc,xxpla,0,1)
! close approach control
           CALL clocms(idc,t1,y1(1:3),y1(nvar2+1:nvar2+3),xxpla) 
        ELSEIF(rhs.eq.2)THEN
           CALL forcesat(y1,y1(nvar2+1),t1,delta(1,m+1-npas),nvar2,idc,xxpla,0,1)
        ELSEIF(rhs.eq.3)THEN
           CALL force9(y1,y1(nvar2+1),t1,delta(1,m+1-npas),nvar2,idc,xxpla,0,1) 
           CALL clocms9(idc,t1)
        ENDIF
! store in differences array        
        if(npas.eq.m)call catst(m,m+1,nvar2,nvar2x,nvar,delta,dd,y1,h,h2)
     endif
     goto 5 
! ***************************************************************       
!  propagatore multistep fa nstep passi, l'ultimo con la velocita'      
  else 
!  occorrono altri passi; quanti?   
      sdelto=deltos*dabs(h)/h 
      din=(t2-t1+sdelto)/h 
      nstep=din 
!  propagazione per nstep passi     
32    call bdnste(t1,y1,h,h2,nstep,m,j1,dd,delta,nvar2,nvar2x,nvar)
!  passo del ms adottato
      npas=npas+nstep 
      t1=t0+h*npas 
   endif
!  nuovo passo          
   goto 5 
! *************************************************************** 
   90 continue 
!  controllo se occorre interpolare 
   if(dabs(t2-t1).le.deltos)then 
!  t2 e' circa t1 : non occorre interpolare     
      y2(1:nvar)=y1(1:nvar) 
   else 
!  occorre interpolare  
      IF(npas.ge.m.and.icmet.eq.1)THEN
! interpolatore multistep
         tint=t2-t1 
         CALL bdintrp(t1,h,tint,m,dd,delta,j1,nvar2,nvar2x,nvar,y2)
      ELSE   
!  con il propagatore runge-kutta usato con passo variato    
         ck2(1:isrk,1:nvar)=0.d0
         lit=lit1 
         epin=eprk 
         tint=t2-t1 
94       call rkimp(t1,tint,y1,dery,ck2,isrk,y2,lit,nvar,epin,ep,lf,ndim)
         if(lf.le.0)then 
!  interpolatore impazzito          
            it=iabs(lf) 
            write(ipirip,1002)tint,epin,(ep(j),j=1,it) 
1002        format( ' interpolatore impazzito,passo h=', &
         &      f10.7,' controllo =',d12.3/' ep= ',5 d12.3/(5d12.3/))       
            epin=epin*10.d0 
            goto 94 
         endif
      ENDIF
   endif
!************************************************        
 END SUBROUTINE propin      
! ******************************************************************    
!  b  d  n  s  t  e        
!           
!  propagatore multistep   
!  bds ordine m+2          
!  esegue nstep passi, all'ultimo calcola anche le velocita'            
!  con aggiornamento catasta a flip-flop  
!           
!  input:   
!           h=passo; h2=h*h
!           m=max ord diff da tenere (predittore ha ord.m+3)
!           nvar=dim di y1; nvar2=nvar/2; nvar2x= maximum
!           fct2=nome secondo membro (forma con solo der. 2e)           
!  input e output:         
!           t1=tempo (esce invariato)     
!           y1=vett stato pos+vel (esce aggiornato se converge)         
!           j1=0,1 indice flip=flop per l'indirizzamento 
!           delta=differenze; dd=somme    
!           
 SUBROUTINE bdnste(t1,y1,h,h2,nstep,m,j1,dd,delta,nvar2,nvar2x,nvar)
   USE force_model, ONLY: velo_req, force
   USE force_sat
   USE force9d
   integer, intent(IN) :: nstep,m ! no. steps, order-2
   integer, intent(inout) :: j1 ! flipflop control, 
   double precision, intent(in) :: h,h2 ! stepsize
!   external fct2 ! name of right hand side
! workspace
   integer, intent(IN) ::  nvar,nvar2,nvar2x 
   double precision, intent(inout) ::  delta(nvar2x,m2max),dd(nvar2x,4)
! result
   double precision, intent(inout) :: y1(nvar), t1 ! state, epoch 
! end interface
   DOUBLE PRECISION yc(nvar2x*2),ycmax,yvmax,tt
!   close approaching planet              
   integer idc 
   double precision xxpla(6) 
! test for need of velocity at each step
!   LOGICAL, EXTERNAL ::  velocity_req
!   indexes for address computations      
   integer id1,id2,in1,in2,kmax,kmin,kj1,kj2,kj2p 
!   loop indexes           
   integer i,k,n 
!   scalar temporary       
   double precision dtemp 
! **************************************  
! static memory not required              
! **************************************  
!           
!  loop sul numero di passi
   tt=t1 
   do 1 n=1,nstep 
!           
!  indirizzamento controllato dal flip-flop j1           
      if(j1.eq.0)then 
         id2=2 
         id1=1 
         in2=4 
         in1=3 
         kmax=m+1 
         kmin=1 
         kj1=0 
         kj2=m+1 
      else 
         id2=4 
         id1=3 
         in2=2 
         in1=1 
         kmax=2*m+2 
         kmin=m+2 
         kj1=m+1 
         kj2=0 
      endif 
!           
!  predittore  ! c=Cow pred b=Cow corr f=Ad pred a=Ad corr 
!           
      if(m.eq.10)then 
         do  i=1,nvar2 
           y1(i)=h2*(delta(i,kmax)*c_ms(kmax)+delta(i,kmax-1)*c_ms(kmax-1)+   &
     &         delta(i,kmax-2)*c_ms(kmax-2)+delta(i,kmax-3)*c_ms(kmax-3)+     &
     &         delta(i,kmax-4)*c_ms(kmax-4)+delta(i,kmax-5)*c_ms(kmax-5)+     &
     &         delta(i,kmax-6)*c_ms(kmax-6)+delta(i,kmax-7)*c_ms(kmax-7)+     &
     &         delta(i,kmax-8)*c_ms(kmax-8)+delta(i,kmax-9)*c_ms(kmax-9)+     &
     &         delta(i,kmax-10)*c_ms(kmax-10)+dd(i,id2))    
         enddo 
      elseif(m.eq.6)then 
         do  i=1,nvar2 
           y1(i)=h2*(delta(i,kmax)*c_ms(kmax)+delta(i,kmax-1)*c_ms(kmax-1)+   &
     &         delta(i,kmax-2)*c_ms(kmax-2)+delta(i,kmax-3)*c_ms(kmax-3)+     &
     &         delta(i,kmax-4)*c_ms(kmax-4)+delta(i,kmax-5)*c_ms(kmax-5)+     &
     &         delta(i,kmax-6)*c_ms(kmax-6)+dd(i,id2))      
         enddo 
      else 
         DO i=1,nvar2 
            dtemp=0.d0 
            DO k=kmax,kmin,-1 
               dtemp=dtemp+delta(i,k)*c_ms(k)
            ENDDO
            y1(i)=(dtemp+dd(i,id2))*h2
         ENDDO
      endif 
      tt=tt+h 
      IF(velo_req)THEN 
         do  i=1,nvar2 
            dtemp=0.d0 
            do k=kmax,kmin,-1 
               dtemp=dtemp+delta(i,k)*f_ms(k) 
            enddo 
            y1(i+nvar2)=h*(dtemp+dd(i,id1)) 
         enddo 
      ENDIF                                           
!   accelerazione nel nuovo punto       
      kj2p=kj2+1 
      IF(rhs.eq.1)THEN       
         CALL force(y1,y1(nvar2+1),tt,delta(1,kj2p),nvar2,idc,xxpla,0,1)
! close approach control
         CALL clocms(idc,tt,y1(1:3),y1(nvar2+1:nvar2+3),xxpla) 
      ELSEIF(rhs.eq.2)THEN
         CALL forcesat(y1,y1(nvar2+1),tt,delta(1,kj2p),nvar2,idc,xxpla,0,1)
      ELSEIF(rhs.eq.3)THEN
         CALL force9(y1,y1(nvar2+1),tt,delta(1,kj2p),nvar2,idc,xxpla,0,1) 
         CALL clocms9(idc,tt)
      ENDIF
! ===============================       
!   aggiornamento catasta               
      DO i=1,nvar2 
         dd(i,in1)=dd(i,id1)+delta(i,kj2p) 
         dd(i,in2)=dd(i,id2)+dd(i,in1) 
      ENDDO 
      DO k=1,m 
         delta(1:nvar2,k+1+kj2)=delta(1:nvar2,k+kj2)-delta(1:nvar2,k+kj1)
      ENDDO 
! ============================================          
!   aggiunto un correttore per controllo
      IF(epms.gt.0.d0)THEN
         ycmax=0.d0             
         DO i=1,nvar2 
            dtemp=0.d0 
            DO k=kmax+kj2-kj1,kmin+kj2-kj1,-1 
               dtemp=dtemp+delta(i,k)*b_ms(k)
            ENDDO
            yc(i)=(dtemp+dd(i,in2)-dd(i,in1))*h2
            IF(i.le.3.and.abs(yc(i)-y1(i)).gt.ycmax)ycmax=abs(yc(i)-y1(i)) 
         ENDDO
         yvmax=0.d0
         IF(velo_req)THEN 
            DO  i=1,nvar2 
               dtemp=0.d0 
               DO k=kmax+kj2-kj1,kmin+kj2-kj1,-1 
                  dtemp=dtemp+delta(i,k)*a_ms(k) 
               ENDDO
               yc(i+nvar2)=h*(dtemp+dd(i,in1)) 
               IF(i.le.3.and.abs(yc(i+nvar2)-y1(i+nvar2)).gt.yvmax)        &
     &                yvmax=abs(yc(i+nvar2)-y1(i+nvar2))
            ENDDO
         ENDIF
         IF(ycmax.gt.epms.or.yvmax.gt.epms)THEN
!           WRITE(*,199)'pred ',tt,ycmax,yvmax,h
!   accelerazione nel nuovo punto       
            kj2p=kj2+1 
            IF(rhs.eq.1)THEN       
               CALL force(yc,yc(nvar2+1),tt,delta(1,kj2p),nvar2,idc,xxpla,1,1)
            ELSEIF(rhs.eq.2)THEN
               CALL forcesat(yc,yc(nvar2+1),tt,delta(1,kj2p),nvar2,idc,xxpla,1,1)
            ELSEIF(rhs.eq.3)THEN
               CALL force9(yc,yc(nvar2+1),tt,delta(1,kj2p),nvar2,idc,xxpla,1,1) 
            ENDIF
! nuovo  aggiornamento catasta               
            DO i=1,nvar2 
               dd(i,in1)=dd(i,id1)+delta(i,kj2p) 
               dd(i,in2)=dd(i,id2)+dd(i,in1) 
            ENDDO
            DO  k=1,m 
               DO i=1,nvar2
                  delta(i,k+1+kj2)=delta(i,k+kj2)-delta(i,k+kj1)
               ENDDO
            ENDDO
!   aggiunto un correttore per convergenza
            ycmax=0.d0             
            DO i=1,nvar2 
               dtemp=0.d0 
               DO k=kmax+kj2-kj1,kmin+kj2-kj1,-1 
                  dtemp=dtemp+delta(i,k)*b_ms(k)
               ENDDO
               y1(i)=(dtemp+dd(i,in2)-dd(i,in1))*h2
               IF(i.le.3.and.abs(yc(i)-y1(i)).gt.ycmax)ycmax=abs(yc(i)-y1(i)) 
            ENDDO
            yvmax=0.d0
            IF(velo_req)THEN 
               DO  i=1,nvar2 
                  dtemp=0.d0 
                  DO k=kmax+kj2-kj1,kmin+kj2-kj1,-1 
                     dtemp=dtemp+delta(i,k)*a_ms(k) 
                  ENDDO
                  y1(i+nvar2)=h*(dtemp+dd(i,in1)) 
                  IF(i.le.3.and.abs(yc(i+nvar2)-y1(i+nvar2)).gt.yvmax)        &
        &   yvmax=abs(yc(i+nvar2)-y1(i+nvar2))
               ENDDO
            ENDIF
            IF(ycmax.gt.epms.or.yvmax.gt.epms)     &
        &      WRITE(ipirip,199) 'corr ',tt,ycmax,yvmax,h
199         FORMAT(a5,f12.4,1p,2d10.2,0p,1x,f6.2)
         ENDIF
      ENDIF
!   passo concluso  
      j1=1-j1 
!   fine loop sul numero di passi
1  ENDDO
!   calcolo velocita', se non gia' fatto      
   IF(.not.velo_req)THEN 
      do  i=1,nvar2 
         dtemp=0.d0 
         do k=kmax,kmin,-1 
            dtemp=dtemp+delta(i,k)*f_ms(k) 
         enddo
         y1(i+nvar2)=h*(dtemp+dd(i,id1)) 
      enddo
   ENDIF
 END SUBROUTINE bdnste
! ========================================================  
 SUBROUTINE bdintrp(t1,h,s,m,dd,delta,j1,nvar2,nvar2x,nvar,y)   
   integer, intent(in) :: m,j1 ! order-2,flipflop control
!   time, fixed stepsize, interpolation step
   double precision, intent(in) :: t1,h,s 
!   workspace       
   integer, intent(in) :: nvar,nvar2,nvar2x 
   double precision,intent(in) :: delta(nvar2x,m2max),dd(nvar2x,4)
   double precision, intent(out) :: y(nvar)
! end interface
   INTEGER i,k,j ! loop indexes
   INTEGER id1,id2,kmin,kmax ! flipflop adresses, address in rotating stack
   DOUBLE PRECISION dtemp !scalar temporary
   DOUBLE PRECISION st ! relative step
   DOUBLE PRECISION c_s(m2max),f_s(m2max)
! relative step
   st=s/h
   IF(st.lt.0.d0.or.st.gt.1.d0)THEN
      STOP '*** bdintrp: st.lt.0.d0.or.st.gt.1.d0 ****'
   ENDIF
!  indirizzamento controllato dal flip-flop j1
   if(j1.eq.0)then 
      id2=2 
      id1=1 
      kmax=m+1 
      kmin=1 
   else 
      id2=4 
      id1=3 
      kmax=2*m+2 
      kmin=m+2
   endif
! computation of interpolation coefficients
   CALL compco_intrp()
!  duplicazione per ridurre i calcoli di indici            
   do  j=1,m+1 
      c_s(j+m+1)=c_s(j) 
      f_s(j+m+1)=f_s(j) 
   enddo
! interpolation : position           
   do i=1,nvar2 
      dtemp=0.d0 
      do  k=kmax,kmin,-1 
         dtemp=dtemp+delta(i,k)*c_s(k) 
      enddo
      y(i)=(dtemp+(st-1.d0)*dd(i,id1)+dd(i,id2))*h*h 
   enddo
! interpolation : velocity
   do  i=1,nvar2 
      dtemp=0.d0 
      do k=kmax,kmin,-1 
         dtemp=dtemp+delta(i,k)*f_s(k) 
      enddo
      y(i+nvar2)=h*(dtemp+dd(i,id1)) 
   enddo
   RETURN
 CONTAINS
   SUBROUTINE compco_intrp()
     DOUBLE PRECISION e_s(0:m2max)
     DOUBLE PRECISION temp,temp2 !scalar temporary
     e_s(0)=1.d0 ! e_s(j)=eta_j(s)
     DO k=1,m+2
        e_s(k)=e_s(k-1)*(st+k-1)/k
     ENDDO
     c_s(1)=1.d0 ! c_s(j)= gamma_{j-1}(s) j=1,m+2
     f_s(1)=1.d0 ! f_s(j)= alpha_{j-1}(s)
     DO j=1,m+1
        temp=e_s(j) -e_s(j-1)/2! e_s(j)*beta(0)+ e_s(j)*beta(1)
        temp2=e_s(j)-e_s(j-1) ! e_s(j)*delta(0)+ e_s(j-1)*delta(1)
        DO k=2,j 
!     sum [eta_{j-k}(s)*beta_{k}]=alpha_{j}(s)=f_s(j+1)
!     sum [eta_{j-k}(s)*delta_{k}]=gamma_{j}(s)=c_s(j+1)
           temp=temp+e_s(j-k)*a_ms(k) ! e_s(j-k)*beta(k)
           temp2=temp2+e_s(j-k)*b_ms(k-1) ! e_s(j-k)*delta(k)
        ENDDO
        f_s(j+1)=temp
        c_s(j+1)=temp2
     ENDDO
     temp2=e_s(m+2)-e_s(m+1)
     DO k=2,m+2
        temp2=temp2+e_s(m+2-k)*b_ms0(k-1)
     ENDDO
     c_s(m+3)=temp2
! che casino questi indici, vanno buttati i primi per Adams,
! i primi due per Cowell
     DO j=1,m+1
        f_s(j)=f_s(j+1)
        c_s(j)=c_s(j+2)
     ENDDO
   END SUBROUTINE compco_intrp
 END SUBROUTINE bdintrp
!=============================== 
! CLOCMS9
! close approach control for orbit9
 SUBROUTINE clocms9(idc,tt)
   USE massmod
   INTEGER, INTENT(IN) :: idc
   DOUBLE PRECISION, INTENT(IN) :: tt
   INCLUDE 'comnbo.h90'
   INTEGER np,n
   np=MOD(idc,nbod-1)
   IF(np.eq.0) np=nbod-1
   n=(idc-np)/(nbod-1)
   if(idc.ne.0)then
!      write(*,*)'t =',tt,' close approach of ',ida(n),' to planet ',nompla(np)
      write(iuncla,*)'t =',tt,' close approach of ',ida(n),' to planet ',nompla(np)
   endif
 END SUBROUTINE clocms9
!=============================== 
! CLOCMS            
! close approach control         
 SUBROUTINE clocms(idc,tt,xa,va,xpla) 
   USE planet_masses
   USE close_app, ONLY : clost
   USE force_model, ONLY: iclap,velo_req ! only for close app. control
   INTEGER,INTENT(IN):: idc ! current planet approached, if any
   DOUBLE PRECISION,INTENT(IN) :: tt,xpla(6),xa(3),va(3) ! time, 
! heliocentric position and velocity, planet pos-vel
   LOGICAL velo_reqsave ! to save state before close app.
! =======end interface=================
   INTEGER ic ! planet currently approached
   DOUBLE PRECISION x(3),v(3) ! planetocentric position and velocity
   DOUBLE PRECISION dcla,vsize,dist,vcla ! distance from planet,velocity
   DOUBLE PRECISION tcla ! time of close app.
! calendar date variables        
   INTEGER iyear,imonth,iday 
   DOUBLE PRECISION hour 
   CHARACTER*16 date 
! memory of close approach state
   SAVE dcla,vcla,tcla,ic
! startup
   IF(clost)THEN 
      ic=0 
      clost=.false.
      v=0.d0 
      velo_reqsave=velo_req
   ENDIF
   IF(iclap.eq.0)RETURN
! separation of orb9/orb11 cases: WARNING: use rhs
!      IF(iorbfit.eq.11)THEN 
!      IF(idc.ne.0)THEN
! to be improved with a real safety feature   
!         write(*,*)'t =',tt,' close approach to=',ordnam(idc) 
!         write(iuncla,*)'t =',tt,' close approach to =',ordnam(idc)
!      ENDIF
! control on close approaches  
   IF(ic.eq.0)THEN 
! close approach not going on; check if it is beginning    
      IF(idc.ne.0)THEN 
! close approach detected (first step inside influence region)          
! relative position and velocity
         x=xa-xpla(1:3)
         IF(velo_req)THEN
            v=va-xpla(4:6)
         ELSE
            velo_reqsave=velo_req
            velo_req=.true.
         ENDIF
! planet with which close approach is taking place         
         ic=idc  
         dcla=vsize(x)
         vcla=vsize(v)
         tcla=tt
      ENDIF
   ELSE
! close approach taking place    
!    control of inconsistences 
      IF(idc.ne.0)THEN 
         IF(idc.ne.ic.or.idc.gt.nmass.or.idc.lt.0)THEN 
            WRITE(*,*)' closap: this should not happen',idc,ic 
         ENDIF
! relative position and velocity, distance 
         x=xa-xpla(1:3)
         v=va-xpla(4:6)
         dist=vsize(x) ! 
         IF(dist.lt.dcla)THEN
            dcla=dist
            vcla=vsize(v)
            tcla=tt
         ENDIF
      ELSE
!   close approach ended; find date, write record
         call mjddat(tcla,iday,imonth,iyear,hour) 
         write(date,'(i4,a1,i2.2,a1,i2.2,f6.5)')   &
 &        iyear,'/',imonth,'/',iday,hour/24d0 
         write(*,199)date,dcla,vcla,ordnam(ic)
199      FORMAT(a16,' dist= ',f7.4,' vel= ',f9.6,' to ',a15) 
         write(iuncla,199)date,dcla,vcla,ordnam(ic)
         numcla=numcla+1 ! to preserve clo file
         ic=0 ! reset state, ready for next close approach
         velo_req=velo_reqsave ! reset velocity computation
      ENDIF
   ENDIF
!      ELSEIF(iorbfit.eq.9)THEN  WARNING: use rhs
!         if(idc.ne.0)then 
!            write(*,*)'t =',tt,' close approach code=',idc 
!            write(iuncla,*)'t =',tt,' close approach code =',idc 
!         endif 
!      ENDIF 
 END SUBROUTINE clocms
! *******************************************************************   
!  {\bf catst} ORB8V
!   riempimento tavola differenze
 SUBROUTINE catst(m,m1,nvar2,nvar2x,nvar,delta,dd,y1,h,h2) 
! order, number of variables (positions only) 
   integer, intent(in) :: m,m1,nvar2,nvar2x,nvar 
! table of differences, sums, state vector      
   double precision, intent(inout) :: delta(nvar2x,m1),dd(nvar2x,2)
   double precision, intent(in) :: y1(nvar)
   double precision, intent(in) :: h,h2 ! stepsize, squared
! loop indexes      
   integer i,k,l 
   double precision dtemp 
!****************   
!   static memory not required   
!****************   
   DO k=1,m 
      DO l=1,m-k+1 
         delta(1:nvar2,m+2-l)=delta(1:nvar2,m+1-l)-delta(1:nvar2,m+2-l)
      ENDDO
   ENDDO
!   prima e seconda somma determinate dalle condizioni iniziali         
   DO 3 i=1,nvar2 
      dtemp=y1(i+nvar2)/h 
      DO k=0,m 
         dtemp=dtemp-a_ms(k+1)*delta(i,k+1)
      ENDDO
      dd(i,1)=dtemp 
      dtemp=y1(i)/h2+dtemp 
      DO k=0,m 
         dtemp=dtemp-b_ms(k+1)*delta(i,k+1)
      ENDDO
      dd(i,2)=dtemp 
3  ENDDO
 END SUBROUTINE catst
! version 3.1 A. Milani Aug 2003   
! ===================================================================== 
! SEL_STE: selection of the stepsize for multistep          
! ===================================================================== 
!  input:           
!         ecc eccentricity       
!         enne mean motion    
!         error allowed error/rev^2
!         mms order-2 for multistep   
!  output:          
!         h recommended stepsize for multistep
! ===================================================================== 
SUBROUTINE sel_ste(ecc,enne,error,mms,hmax,h) 
! input/output      
  DOUBLE PRECISION, INTENT(IN) :: ecc,enne,error,hmax
  DOUBLE PRECISION, INTENT(OUT) :: h 
  integer mms 
! scalar temporaries
  double precision eps,econv,hh,z,emul,step,err 
  integer ila,igr,nb 
!   static memory not required   
! ===================================================================== 
! the truncation error will be of the order of error       
! per revolution squared; however, error cannot be less than            
! machine accuracy, otherwise the rounding off would be the dominant    
! source of integration error.   
  eps=epsilon(1.d0)
  eps=max(eps,error) 
  hh=hmax 
! ===================================================================== 
! for each orbit, compute the truncation error in longitude
! according to Milani and Nobili, 1988, Celest. Mech. 43, 1--34         
! and select a stepsize giving a truncation  error equivalent to        
! one roundoff per revolution squared         
  econv=1.d0 
!      write(ipirip,*)' Selection of the best stepsize'    
  call  zed(ecc,mms+2,z,econv,ila,igr) 
! if the series with Bessel functions is divergent, e.g. for e=0.2      
! with iord=12, some wild guess is used because the error  
! estimation formula is no good anyway; hope in your luck to have       
! a not too bad orbit            
!     if(ila.ge.20)z=1.d6        
  if(mod(mms,2).eq.0)then 
     emul=dpig**2*3/2*cerr*z 
     step=(eps/emul)**(1.d0/(mms+3))/enne 
     err=emul*(enne*hh)**(mms+3) 
  else 
     emul=dpig**2*3/4*mms*cerr*z 
     step=(eps/emul)**(1.d0/(mms+4))/enne 
     err=emul*(enne*hh)**(mms+4) 
  endif
!     write(ipirip,*)enne,ecc,z  
! ===================================================================== 
! now compare with h given in input           
  h=min(step,hmax)
  IF(verb_pro.gt.9)THEN 
     write(ipirip,110)hmax,h 
110  format(' max. step required ',f9.6,' selected ',f9.6)
  ENDIF 
END SUBROUTINE sel_ste
         
! **********************************************************            
!   {\bf compco} ORB8V           
!      
!   calcolo coefficienti multistep            
!   cs=Cow pred bs=Cow corr fs=Ad pred as=Ad corr          
!   versione con calcoli real*8 (trasportabile)            
!   per calcoli in quadrupla precisione cambiare:          
!      - tutti i .d0 in .q0      
!      - real*16 a(mmax) etc     
!   m1=m+1          
!   i coefficienti per le differenze fino     
!   all'ordine iord-2=m sono calcolati e restituiti        
!   negli array con indice da 1 a m+1         
      SUBROUTINE compco(m1,cs,fs,bs,as) 
      integer m,m1,i,j,i1,i2,k 
      double precision a(mmax),b(mmax),c(mmax),d(mmax) 
      double precision as(m1),bs(m1),cs(m1),fs(m1) 
!****************   
!   static memory not required (used only once)            
!****************   
! formulas refer to Milani 7 Nobili, CM 43(1988), 1-34 
      m=m1-1 
! a(j)=beta(j-1) (A.5)
      a(1)=1.d0 
      do 10 i=2,m+3 
        k=i-1 
        a(i)=0.d0 
        do 10 j=1,k 
   10     a(i)=a(i)-a(j)/(k-j+2.d0) 
! b(j)=delta(j-1) (A.11)
      do 11 j=1,m+3 
        b(j)=0.d0 
        do 11 k=1,j 
   11     b(j)=b(j)+a(j+1-k)*a(k) 
      do 12 i=1,m+3 
! c(j)=alpha(j-1) (A.8)
! d(j)=gamma(j-1) (A.13)
        c(i)=0.d0 
        d(i)=0.d0 
        do 12 j=1,i 
          c(i)=c(i)+a(j) 
   12     d(i)=d(i)+b(j) 
! now indexes have to be shifted down, by 1 for Adams, by 2 for Cowell
! to be used in the first and second sum formulae
      do 17 i=2,m+2 
        i1=i-1 
        as(i1)=a(i) !as(j)=a(j+1)=beta(j) j=1,m+1
   17   fs(i1)=c(i) !fs(j)=c(j+1)=alpha(j) j=1,m+1
      do 18 i=3,m+3 
        i2=i-2 
        bs(i2)=b(i) !bs(j)=b(j+2)=delta(j+1) j=1,m+1
   18   cs(i2)=d(i) !cs(j)=d(j+2)=gamma(j+1) j=1,m+1
!     write(6,200)(cs(i),fs(i),bs(i),as(i),i=1,m+1)        
!200  format(4f19.16)            
      return 
      END SUBROUTINE compco   
!********************************************************************** 
!  {\bf zed}  ORB8V 
! Calcolo di Z(m,e) -- correzione all'errore di troncamento
! dovuta alla eccentricita'. m e' l'ordine del multistep nella          
! forma Stormer predictor; in LONGSTOP m=12.  
! eps controllo di convergenza.  
!      
! Routine che calcola Z(m,e);imax numero massimo iterazioni (dato 20)   
!********************************************************************** 
SUBROUTINE zed(e,m,f,eps,i,igr)
  INTEGER, INTENT(IN) ::  m ! order -2 
  INTEGER, INTENT(OUT) :: i,igr ! last term, largest term
  DOUBLE PRECISION, INTENT(IN) :: e,eps ! ecc, control convergence of series
  DOUBLE PRECISION, INTENT(OUT) :: f ! Z(e,m)
  INTEGER, PARAMETER :: imax=20
  DOUBLE PRECISION x,adf,ri,cie,sie,c2ie,s2ie,sum,df,vadf,ff
  f=0.d0 
  adf=0.d0 
 ! IF(e.gt.0.27d0.and.e.lt.0.582d0)THEN 
 !    f=1.d6 
 !    i=imax+1
 !    igr=imax
 !    IF(verb_pro.gt.10)THEN
 !       WRITE(ipirip,*)' with e>0.27 flat stepsize'
 !    ENDIF
 !    RETURN 
 ! ELSEIF(e.ge.0.582d0.and.e.lt.0.8d0)THEN 
 !    ff=sqrt(1.d0-e**2)/(1.d0-e)**2
 !    f=ff**(m+1)
 !    i=imax+1
 !    igr=imax
 !    WRITE(ipirip,*)' with e>0.582 true anomaly at pericenter' 
 !    RETURN 
 ! ELSEIF(e.ge.0.8d0)THEN
 !    ff=sqrt(1.d0-0.8d0**2)/(1.d0-0.8d0)**2
 !    f=ff**(m+1)
 !    i=imax+1
 !    igr=imax
 !    IF(verb_pro.gt.10)THEN
 !       WRITE(ipirip,*)' with e>0.8 flat stepsize'
 !    ENDIF
 !    RETURN 
 ! ENDIF
  ff=sqrt(1.d0-e**2)/(1.d0-e)**2
  f=ff**(m+1)
  IF(f.gt.1.d6)THEN
     i=imax+1
     igr=imax
     RETURN
  ENDIF
  do 1 i=1,imax 
     x=i*e 
     ri=i 
     cie=(bessel_jn(i-1,x)-bessel_jn(i+1,x))/i 
     sie=(bessel_jn(i-1,x)+bessel_jn(i+1,x))/i 
     c2ie=cie*cie 
     s2ie=sie*sie 
     sum=(c2ie+s2ie)/2.d0 
     df=sum*ri**(m+4) 
     f=f+df 
!c      write(*,100)i,df         
!c100   format(i5,d18.6)         
     vadf=adf 
     adf=dabs(df) 
     if(adf.gt.vadf)igr=i 
     if(adf.lt.eps.and.i.gt.1)return 
1 enddo
  write(ipirip,101)i-1,df,igr 
101 format(' non convergence in zed; last term=',i5,d18.6/            &
         &       ' max was for ',i4) 
  f=1.d6 
END SUBROUTINE zed   
!********************************************************************** 
!  {\bf bessel_jn}  ORB8V           
! Function che calcola la funzione di Bessel J(i,x).       
! imax numero massimo di iterazioni (dato 20) 
!********************************************************************** 
double precision function bessel_jn(i,x) 
  INTEGER, INTENT(IN) :: i
  DOUBLE PRECISION, INTENT(IN) :: x
  INTEGER, PARAMETER :: imax=20
  DOUBLE PRECISION, PARAMETER  :: epbs=1.d-10
  DOUBLE PRECISION x2,x2p,dbess
  INTEGER j,jpifat, jfat ,ifl,ifat
  x2=x/2.d0 
  x2p=x2 
  ifat=1 
  DO j=2,i 
     x2p=x2p*x2 
     ifat=ifat*j 
  ENDDO
  if(i.eq.0) x2p=1.d0 
  jfat=1 
  jpifat=1 
  ifl=1 
  bessel_jn=0.d0 
  do 2 j=1,imax 
     dbess=x2p*ifl/ifat 
     dbess=dbess/jfat 
     dbess=dbess/jpifat 
!c      write(*,*)dbess          
     bessel_jn=bessel_jn+dbess 
!c      write(*,*)i,j,ifat,jfat,jpifat        
     if(dabs(dbess).lt.epbs)return 
     jpifat=jpifat*(i+j) 
     jfat=jfat*j 
     ifl=-ifl 
     x2p=x2p*x2*x2 
2 ENDDO
  write(ipirip,101)j-1,dbess 
101 format('bessel_jn: non convergence; last term=',i5,d18.6) 
END function bessel_jn
! ================================            
! invaxv            
!      
! initialise the  variation matrix as the 6 X 6 identity  
! plus columns of zeros for the derivatives w.r. to dynamical parameters 
! ================================            
SUBROUTINE invaxv(x,v,nvar2) 
  INTEGER, INTENT(IN) :: nvar2 
  DOUBLE PRECISION, INTENT(INOUT) ::  x(nvar2),v(nvar2) 
! end interfface    
  INTEGER i,j,iii,ij,nd
  iii=3 
  do 7 j=1,6 
     do  i=1,3 
        ij=i+3*(j-1) 
        x(iii+ij)=0.d0 
        v(iii+ij)=0.d0 
        if(i.eq.j)then 
           x(iii+ij)=1.d0 
        elseif(j.eq.i+3)then 
           v(iii+ij)=1.d0 
        endif
     enddo
7 ENDDO
! case of additional parameters to be solved
  IF(nvar2.gt.21)THEN
! initialize array of zeros as derivatives of state with respect to
! dynamical parameters
     nd=(nvar2-3)/3
     DO j=7,nd
       DO i=1,3
         ij=i+3*(j-1)
         x(iii+ij)=0.d0 
         v(iii+ij)=0.d0 
       ENDDO
     ENDDO
  ENDIF
ENDSUBROUTINE invaxv
!======================================
! inivar            
!      
! initialise the  variation matrix as the 6 X 6 identity
! plus columns of zeros for the derivatives w.r. to dynamical parameters   
! ================================            
SUBROUTINE inivar(y1,nvar2,nvar) 
  INTEGER,INTENT(IN) ::  nvar,nvar2 
  DOUBLE PRECISION, INTENT(INOUT) ::  y1(nvar) 
! end interfface    
  INTEGER i,j,iii,ij,nd 
! case of initial conditions as only parameters to be solved
  iii=3 
  DO j=1,6 
     DO i=1,3 
        ij=i+3*(j-1) 
        y1(iii+ij)=0.d0 
        y1(iii+ij+nvar2)=0.d0 
        IF(i.eq.j)THEN 
           y1(iii+ij)=1.d0 
        ELSEIF(j.eq.i+3)THEN 
           y1(iii+ij+nvar2)=1.d0 
        ENDIF
     ENDDO
  ENDDO
! case of additional parameters to be solved
  IF(nvar2.gt.21)THEN
! initialize array of zeros as derivatives of state with respect to
! dynamical parameters
     nd=(nvar2-3)/3
     DO j=7,nd
       DO i=1,3
         ij=i+3*(j-1)
         y1(iii+ij)=0.d0 
         y1(iii+ij+nvar2)=0.d0 
       ENDDO
     ENDDO
  ENDIF
END SUBROUTINE inivar
! ====================================================                  
! CLOTEST                                                               
! tests for presence of a close approach at the time t0 for             
! asteroid with equinoctal elements east. On output iplnet is the number
! of the approaching planet, 0 if none.                                 
! ********************************************************              
! WARNING: this version actually checks only for close approaches to the Earth
! should be fixed to test for all planets.                              
!=====================================================                  
SUBROUTINE clo_test(el0,iplanet,dist0) 
  USE planet_masses
  USE orbit_elements
! input: time, asteroid elements
  TYPE(orbit_elem), INTENT(IN) :: el0  
! output: planet being approached (0 if none), time to get clear        
  INTEGER, INTENT(OUT) ::  iplanet 
  DOUBLE PRECISION, OPTIONAL, INTENT(OUT) :: dist0 ! distance
! end interface                                                         
  DOUBLE PRECISION xea(6),xast(6),dist
  INTEGER fail_flag
  TYPE(orbit_elem) :: elcar
! from planet_masses.mod mass of the sun from fund_const
! radius of the close approach sphere for the planets                   
! JPL Earth vector at observation time                                  
  IF(rhs.eq.1)THEN
     CALL earcar(el0%t,xea,1)
  ELSEIF(rhs.eq.2)THEN
     xea=0.d0
  ELSEIF(rhs.eq.3)THEN
! should not happen
     STOP
  ENDIF 
! cartesian coordinates of the asteroid                                 
  CALL coo_cha(el0,'CAR',elcar,fail_flag) 
  xast=elcar%coord
! distance                                                              
  dist=sqrt((xea(1)-xast(1))**2+(xea(2)-xast(2))**2+(xea(3)-xast(3))**2)
  IF(PRESENT(dist0))THEN
     dist0=dist
  ENDIF
! test                                                                  
  IF(dist.le.dmea)THEN 
!        WRITE(*,*)' clotest: warning! close approach to Earth'         
!        WRITE(*,*)' at the initial epoch, dist=',dist                  
     iplanet=3
! change 5 January 2004: time to exit from close app too short,
! this is a fix waiting to use aftclov 
!     texit=20.d0 
!     texit=50.d0 
! fixed by using aftclov
  ELSE 
     iplanet=0 
!     texit=0.d0 
  ENDIF
END SUBROUTINE clo_test

END MODULE propag_state          
! ==========================================================            
! SET_RESTART       
! with argument logical, true or false        
! if true forces restart anyway  
! if false, restart is avoided it it is a continuation orbit            
! as controlled in propag        
! ======================================================== 
SUBROUTINE set_restart(res_log)
  USE propag_state 
  IMPLICIT NONE 
  LOGICAL res_log 
  restar=res_log 
END SUBROUTINE set_restart
! ====================================================     
! varwra            
!      
! First parzial derivatives: rewrap vector into 6xnd matrix 
! ====================================================     
SUBROUTINE varwra(y2,dxdx0,nd,nvar,nvar2) 
  IMPLICIT NONE 
  INTEGER,INTENT(IN) :: nvar,nvar2,nd 
  DOUBLE PRECISION, INTENT(INOUT) ::  y2(nvar) 
  DOUBLE PRECISION, INTENT(OUT) ::  dxdx0(6,nd) 
  INTEGER i,j,ij
! ==================================================== 
! check dimensions
  IF(nvar2.ne.3+3*nd)THEN
     WRITE(*,*)' varwra: nd,nvar2 ', nd,nvar2
     STOP
  ENDIF
! portion of the state transition matrix about initial conditions
  DO j=1,3 
     DO  i=1,3 
        ij=i+3*(j-1)+3 
        dxdx0(i,j)=y2(ij) 
        dxdx0(i,j+3)=y2(ij+9) 
        dxdx0(i+3,j)=y2(ij+nvar2) 
        dxdx0(i+3,j+3)=y2(ij+9+nvar2) 
     ENDDO
  ENDDO
! portion about dynamical parameters
  IF(nd.gt.6)THEN
     DO j=7,nd
        DO i=1,3
          ij=i+3*(j-1)+3       
          dxdx0(i,j)=y2(ij)
          dxdx0(i+3,j)=y2(ij+nvar2)
        ENDDO
     ENDDO
  ENDIF
END  SUBROUTINE varwra
! ====================================================     
! varunw            
!      
! First parzial derivatives: unwrap 6x6 matrix into vector 
! ====================================================     
SUBROUTINE varunw(dxdx0,x,v,nd,nvar2) 
  INTEGER, INTENT(IN) :: nvar2,nd 
  DOUBLE PRECISION, INTENT(IN) :: dxdx0(6,nd)
  DOUBLE PRECISION, INTENT(INOUT) :: x(nvar2),v(nvar2) 
  INTEGER i,j,ij
! ==================================================== 
! check dimensions
  IF(nvar2.ne.3+3*nd)THEN
     WRITE(*,*)' varunw: nd,nvar2 ', nd,nvar2
     STOP
  ENDIF
  DO j=1,3 
     DO  i=1,3 
        ij=i+3*(j-1)+3 
        x(ij)=dxdx0(i,j) 
        x(ij+9)=dxdx0(i,j+3) 
        v(ij)=dxdx0(i+3,j) 
        v(ij+9)=dxdx0(i+3,j+3) 
     ENDDO
  ENDDO
! portion about dynamical parameters
  IF(nd.gt.6)THEN
     DO j=7,nd
        DO i=1,3
          ij=i+3*(j-1)+3       
          x(ij)=dxdx0(i,j)
          v(ij)=dxdx0(i+3,j)
        ENDDO
     ENDDO
  ENDIF
END SUBROUTINE varunw

! ====================================================     
! vawrxv
!            
! First parzial derivatives: rewrap vector into 6x6 matrix 
! ====================================================     
SUBROUTINE vawrxv(x,v,dxdx0,nvar2,nd) 
  INTEGER, INTENT(IN) ::  nvar2, nd
  DOUBLE PRECISION, INTENT(IN) :: x(nvar2),v(nvar2) 
  DOUBLE PRECISION, INTENT(OUT) ::  dxdx0(6,nd)
!  DOUBLE PRECISION, INTENT(OUT) ::  dxdx0(6,6)
  INTEGER i,j,ij 
! ====================================================     
! check dimensions
  IF(nvar2.ne.3+3*nd)THEN
     WRITE(*,*)' vawrxv: nd,nvar2 ', nd,nvar2
     STOP
  ENDIF
  DO j=1,3 
     DO  i=1,3 
        ij=i+3*(j-1)+3 
        dxdx0(i,j)=x(ij) 
        dxdx0(i,j+3)=x(ij+9) 
        dxdx0(i+3,j)=v(ij) 
        dxdx0(i+3,j+3)=v(ij+9) 
     ENDDO
  ENDDO
! portion about dynamical parameters
  IF(nd.gt.6)THEN
     DO j=7,nd
        DO i=1,3
           ij=i+3*(j-1)+3       
           dxdx0(i,j)=x(ij)
           dxdx0(i+3,j)=v(ij)
        ENDDO
     ENDDO
  ENDIF
END SUBROUTINE vawrxv


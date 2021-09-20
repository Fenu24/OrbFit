! ===========MODULE close_app=============   
                                    
! MODULE CONTAINS  
! PUBLIC ROUTINES  
!            cloapp                                 
! PRIVATE ROUTINES 
!               stepcon                                                 
!               falsi                                                   
!               strclo  
!               falsi_surf
!               strsurf
!               geod_hsurf
! 
! OUT OF MODULE
!            npoint_set             
!            cov_avai/cov_not_av      
!            set_clost 
!
! HEADERS AND MODULES  
! close_app.o: \
!	../include/nvarx.h90 \
!	../suit/FUND_CONST.mod \
!	close_app.o \
!	tp_trace.o \
!	runge_kutta_gauss.o 


MODULE close_app

USE fund_const
USE output_control
IMPLICIT NONE
PRIVATE

! PUBLIC SUBROUTINES

PUBLIC :: cloapp

! PUBLIC DATA


! parameters with controls for close app.
DOUBLE PRECISION eprdot
! Minimum number of data points for a deep close appr
INTEGER npoint
! startup from a non close approaching state at each integration restart
LOGICAL clost
! fix mole flag... to be used only by resret2...find close app min (vsa_attx)
LOGICAL fix_mole, kill_propag, min_dist
! logical to activate impact corridor computation
LOGICAL surf_stop

PUBLIC eprdot,npoint,clost,fix_mole, kill_propag, min_dist, surf_stop

! PRIVATE COMMON DATA   

CONTAINS
!                                                                       
! ==============================================================        
! CLOAPP                                                                
! close approach control; driver for falsi                              
! vers. 1.9.2, A.Milani, August 25, 1999                                
! with decreasing stepsize near closest approach                        
! also with end of close approach flag                                  
! ==============================================================        
  SUBROUTINE cloapp(tcur,th,xa,va,nv,idc,xpla,xldir,dir,nes,cloend) 
    USE tp_trace
    USE planet_masses
! close approach control: which planet is close (0=none)                
    INTEGER idc 
! stepsize control: flag for fixed step, stepsize, direction            
    LOGICAL nes 
    DOUBLE PRECISION xldir,dir 
! current time, previous stepsize                                       
    DOUBLE PRECISION tcur,th 
! positon, velocity of asteroid, of close encounter planet              
    INTEGER nv 
    DOUBLE PRECISION xa(nv),va(nv) 
    DOUBLE PRECISION xpla(6) 
! logical flag for close approach termination                           
! to allow for storage of state transition matrix                       
    LOGICAL cloend 
! =========END INTERFACE========================================        
! stepsize tricks                                                       
    LOGICAL nesold 
    DOUBLE PRECISION xlold,tmin,dt,thfirst 
!               ,tleft                                                  
! which planet is being approached, vector difference                   
    INTEGER ic 
    DOUBLE PRECISION x(3),v(3) 
!                                                                       
    DOUBLE PRECISION vsize 
! counters for arrays of multiple minima                      
    INTEGER jc 
! loop indexes i=1,3                                                    
    INTEGER i 
! logical flags for regulae falsi being initiaited                      
    LOGICAL first 
! static memory model required                                          
    SAVE 
! ===========================================                           
    IF(clost)THEN 
       ic=0 
       clost=.false. 
    ENDIF
! close approach is not ending, in most cases                           
    cloend=.false. 
! ========================================================              
! control on close approaches                                         
    IF(ic.eq.0)THEN 
! close approach not going on; check if it is beginning                 
       IF(idc.ne.0)THEN 
! close approach detected (first step inside influence region)          
! relative position and velocity  
          x=xa(1:3)-xpla(1:3)
          v=va(1:3)-xpla(4:6)         
! planet with which close approach is taking place                      
          ic=idc 
          iplam=idc 
!           WRITE(iuncla,*)' approach to planet', iplam                 
! close approach minima counter                                         
          jc=0 
! first call to falsi                                                   
          first=.true. 
          IF(min_dist)THEN
             CALL falsi(tcur,xa,va,nv,xpla,jc,first,iplam)
             njc=jc
          END IF
          IF(surf_stop.AND.iplam.EQ.3) CALL falsi_surf(tcur,xa,va,nv,xpla,first,iplam)
! setup of fixed stepsize                                               
          nesold=nes 
          nes=.true. 
          xlold=xldir 
!  stepsize based upon angular velocity                                 
          CALL stepcon(x,v,npoint,dmin(idc),dir,dt) 
! not allowed to be longer than the one used before the close approach  
          xldir=dir*min(abs(th),abs(dt)) 
! which is stored for alter use                                         
          thfirst=th 
       ENDIF
    ELSE 
! close approach taking place                                           
!    control of inconsistences                                          
       IF(idc.ne.0)THEN 
          IF(idc.ne.ic.or.idc.gt.nmass.or.idc.lt.0)THEN 
             WRITE(*,*)' cloapp: this should not happen',idc,ic 
          ENDIF
! relative position and velocity                                        
          x(1:3)=xa(1:3)-xpla(1:3) 
          v(1:3)=va(1:3)-xpla(4:6) 
          first=.false. 
          IF(min_dist)THEN
             CALL falsi(tcur,xa,va,nv,xpla,jc,first,iplam)
             njc=jc
          END IF
          IF(surf_stop.AND.iplam.EQ.3) CALL falsi_surf(tcur,xa,va,nv,xpla,first,iplam)
! mole fix case: check kill_propag
          IF(kill_propag)THEN
             cloend=.true. 
             xldir=xlold 
             nes=nesold 
             ic=0 
             njc=jc 
             jc=0 
             RETURN
          ENDIF 
!  stepsize based upon angular velocity                                 
          CALL stepcon(x,v,npoint,dmin(idc),dir,dt) 
! not allowed to be longer than the one used before the close approach  
          xldir=dir*min(abs(thfirst),abs(dt)) 
       ELSE 
!   close approach ended; reset flags etc.                              
          cloend=.true. 
          xldir=xlold 
          nes=nesold 
          ic=0 
          njc=jc 
          jc=0 
       ENDIF
    ENDIF
! =================================================                     
  END SUBROUTINE cloapp
! =================================================                     
! STEPCON                                                               
! numerical stepsize control based on true anomaly                      
! =================================================                     
      SUBROUTINE stepcon(x,v,npoin,disc,dir,dt) 
! input: geocentric state                                               
      DOUBLE PRECISION x(3),v(3) 
! input: min no points, size of MTP disc, direction of time             
      INTEGER npoin 
      DOUBLE PRECISION disc,dir 
! outpput: suggested stepsize                                           
      DOUBLE PRECISION dt 
! end interface                                                         
      DOUBLE PRECISION prscal,vsize 
      DOUBLE PRECISION dtheta,cosdt1,dtheta1,tmin,dtold 
      DOUBLE PRECISION x1(3),cgeo(3),jgeo 
      dtheta=pig/npoin 
      CALL prvec(x,v,cgeo) 
      jgeo=vsize(cgeo) 
      dt=prscal(x,x)*dtheta/jgeo*dir 
      CALL lincom(x,1.d0,v,dt,x1) 
      cosdt1=prscal(x,x1)/(vsize(x)*vsize(x1)) 
      dtheta1=acos(cosdt1) 
      IF(dtheta1.gt.dtheta)THEN 
         dt=(dtheta/dtheta1)*dt 
      ENDIF 
      tmin=2*disc/vsize(v) 
      dtold=tmin/npoin*dir
      IF(abs(dt).lt.1.d-5)THEN 
         WRITE(iun_log,*)'stepcon: ',dt, dtold,vsize(x),vsize(v),dtheta1,jgeo  
         dt=1.d-5
         IF(.not.min_dist.and.(.not.surf_stop))kill_propag=.true.
      ENDIF     
      END SUBROUTINE stepcon                                          
! =====================================                                 
! FALSI  regula falsi                                                   
! =====================================                                 
 SUBROUTINE falsi(t,xa,va,nv,xpla,jc,first,iplam) 
   USE output_control
   USE runge_kutta_gauss
   USE planet_masses
   USE critical_points
   USE orbit_elements
! INPUT                                                                 
! current position and velocity                                         
   INTEGER nv 
   DOUBLE PRECISION t,xa(nv),va(nv),xpla(6) 
! is this the beginning of an integration?                              
   LOGICAL first 
! planet being approached                                               
   INTEGER iplam 
! INPUT AND OUTPUT                                                      
! close approach counter                              
   INTEGER jc
! =====================================                                 
! fixed stepsize time interval                                          
   DOUBLE PRECISION tt1,tt2 
! data for regula falsi                                                 
   DOUBLE PRECISION r1,r2,r0,rdot1,rdot2,rdot0,t1,t2,t0 
   DOUBLE PRECISION z1,z2,z0,zdot1,zdot2,zdot0 
! functions                                                             
   DOUBLE PRECISION vsize,prscal 
! state vectors: with partials 
   INCLUDE 'nvarx.h90' 
   DOUBLE PRECISION x(nvar2x),v(nvar2x),xt(nvar2x),vt(nvar2x),xplat(6) 
   DOUBLE PRECISION xat(nvar2x),vat(nvar2x) 
! time steps                                                            
   DOUBLE PRECISION dt,tt,di,hh 
! iterations                                                            
   INTEGER it,i 
! ======for ID========                                                
!   DOUBLE PRECISION moid0 
! variables for call to dmintil_rms                                  
   DOUBLE PRECISION x6(6) 
   TYPE(orbit_elem) el1,el2
! cartesian coordinates asteroid and planet at minimum                  
   DOUBLE PRECISION c1min(3,nminx),c2min(3,nminx) 
!     SQUARED DISTANCE function                                         
   DOUBLE PRECISION dmintil(nminx) 
!     number of relative minima found                                   
   INTEGER nummin 
   INTEGER, PARAMETER :: itmax_f=50 !iterations of falsi
   DOUBLE PRECISION :: peri
!
! ====================                                                
! memory model is static                                                
   SAVE 
! default for flag stopping propagation due to collision
   kill_propag=.false.
! planetocentric position with derivatives
   x(1:nv)=xa(1:nv)
   v(1:nv)=va(1:nv)
   x(1:3)=xa(1:3)-xpla(1:3) 
   v(1:3)=va(1:3)-xpla(4:6) 
! initialisation at the beginning of each close approach                
   IF(first)THEN 
      tt1=t 
      r1=vsize(x) 
      rdot1=prscal(x,v)/r1 
! MOID needed also when not propagating covariance
!      IF(nv.ge.21)THEN
! MOID at the beginning of each encounter                               
      x6(1:3)=xa(1:3) 
      x6(4:6)=va(1:3)      
      el2=undefined_orbit_elem
!      IF(rhs.EQ.2)THEN
!         el2%center=3
!      END IF
      el2%t=t
      el2%coo='CAR'
      el1=el2
      el2%coord=x6 ! asteroid cartesian elements
      el1%coord=xpla
      CALL dmintil_rms(el1,el2,nummin,dmintil,c1min,c2min)
!      ENDIF 
! end initialisation of close approach                                  
      RETURN 
   ENDIF
! compute r, rdot                                                       
   r2=vsize(x) 
   rdot2=prscal(x,v)/r2 
   tt2=t 
! direction of time propagation                                         
   IF(tt2.gt.tt1)THEN 
      di=1.d0 
   ELSEIF(tt2.lt.tt1)THEN 
      di=-1.d0 
   ELSE 
      WRITE(iun_log,199)tt1,tt2,x(1:3) !,v,xpla
 199  FORMAT(' falsi: zero step ',f8.2,1x,f8.2,1P,3(1x,d10.3))
      RETURN 
   ENDIF
! if inside Earth, use 2-body propagation
   IF(iplam.eq.3.and.r2.lt.reau.and.fix_mole)THEN
! WARNING: in this way the minimum is not found
! let this be handled by strclo       
      CALL strclo(iplam,t,xpla,xa,va,nv,jc,r2,rdot2,         &
     &        dmintil,c1min,nummin)                                         
      rdot2=0.d0 
! END WARNING
! check for minimum distance                                            
   ELSEIF(abs(rdot2).lt.eprdot)THEN 
! already at stationary point; WARNING: is it a minimum?                
      CALL strclo(iplam,t,xpla,xa,va,nv,jc,r2,rdot2,         &
     &        dmintil,c1min,nummin)                                         
      rdot2=0.d0 
   ELSEIF(rdot1*di.lt.0.d0.and.rdot2*di.gt.0.d0)THEN 
! rdot changes sign, in a way appropriate for a minimum                 
      t1=tt1 
      t2=tt2 
!        WRITE(*,*) ' r. f. begins, t,r,rdot=',t1,r1,rdot1,t2,r2,rdot2  
!        WRITE(iuncla,*) ' r. f. begins ',t1,r1,rdot1,t2,r2,rdot2       
      dt=-rdot2*(t1-t2)/(rdot1-rdot2) 
      tt=dt+t2 
      hh=tt-t 
! iterate regula falsi                                                  
      DO it=1,itmax_f 
         CALL rkg(t,xa,va,nv,hh,xat,vat,xplat) 
! planetocentric position                                               
         xt(1:3)=xat(1:3)-xplat(1:3) 
         vt(1:3)=vat(1:3)-xplat(4:6) 
         r0=vsize(xt) 
         rdot0=prscal(xt,vt)/r0 
         t0=tt 
! selection of next couple of points with opposite sign                 
         IF(abs(rdot0*dt).lt.eprdot)THEN 
!  already at stationary point; WARNING: is it a minimum? 
            CALL strclo(iplam,tt,xplat,xat,vat,nv,jc,r0,rdot0,&
     &             dmintil,c1min,nummin)                                    
            GOTO 2 
         ELSEIF(rdot0*rdot1.lt.0.d0)THEN 
            r2=r0 
            rdot2=rdot0 
            t2=t0 
         ELSEIF(rdot0*rdot2.lt.0.d0)THEN 
            r1=r0 
            rdot1=rdot0 
            t1=t0 
         ENDIF
         dt=-rdot2*(t1-t2)/(rdot1-rdot2) 
         tt=dt+t2 
         hh=tt-t 
      ENDDO
! failed convergence                                                    
      CALL strclo(iplam,t0,xplat,xat,vat,nv,jc,r0,rdot0,     &
     &        dmintil,c1min,nummin)                                         
      WRITE(ierrou,*)' falsi: failed convergence for ',ordnam(iplam) 
      WRITE(ierrou,*)'t1,r1,rdot1,t2,r2,rdot2,dt,tt' 
      WRITE(ierrou,111)t1,r1,rdot1,t2,r2,rdot2,dt,tt 
!        WRITE(*,*)' falsi: failed convergence for ',ordnam(iplam)      
!        WRITE(*,*)'t1,r1,rdot1,t2,r2,rdot2,dt,tt'                      
!        WRITE(*,111)t1,r1,rdot1,t2,r2,rdot2,dt,tt                      
111   format(8(1x,f16.9)) 
      numerr=numerr+1 
   ENDIF
! found zero                                                            
2  CONTINUE 
! save current state to be previous state next time                     
   tt1=tt2 
   r1=r2 
   rdot1=rdot2 
 END SUBROUTINE falsi
! ===================================================                    
! TP_FSER
!  xt,vt = geocentric coordinates and partials at time t0
!  xat,vat= geocentric cood and partials at time t0+dt of perigee
! =================================================== 
  SUBROUTINE  tp_fser(nv,t0,xt,vt,iplam,xpla,xat,vat,dt)
    USE ever_pitkin
    USE planet_masses
    USE dyn_param, ONLY: ndimx
    INTEGER, INTENT(IN) :: nv,iplam
    DOUBLE PRECISION, INTENT(IN) :: t0,xt(nv),vt(nv)
    DOUBLE PRECISION, INTENT(OUT) :: dt, xpla(6),xat(nv),vat(nv)
! ===================== 
    DOUBLE PRECISION, DIMENSION(6,6) :: dxvdxtvt ! only initial conditions
    DOUBLE PRECISION, DIMENSION(6,ndimx) :: dxvdx0,dxdx0 ! also possibly dyn.par
    DOUBLE PRECISION, DIMENSION(3) :: vlenz
    DOUBLE PRECISION :: vel2,gei2,gei,ecc,peri,sig0,alpha,psi,r0,conv_contr
!  ang era definito scalare....
    DOUBLE PRECISION, DIMENSION(3) :: ang
! functions                                                             
    DOUBLE PRECISION vsize,prscal
    INTEGER nd 
    r0=vsize(xt)
    vel2=prscal(vt,vt)
    call prvec(xt,vt,ang)
    gei2=prscal(ang,ang)
    gei=sqrt(gei2)
    IF(gei.eq.0.d0)THEN 
       WRITE(*,*) ' tp_fser: zero angular momentum ',xt,vt,t0 
    ENDIF
    call prvec(vt,ang,vlenz)                    ! Lenz vector
    vlenz=vlenz/gm(iplam)-xt(1:3)/r0
    ecc=vsize(vlenz)
    peri=gei2/(gm(iplam)*(1.d0+ecc))
    sig0=prscal(xt(1:3),vt(1:3))
    alpha=vel2-2*gm(iplam)/r0 ! 2* energy
!    WRITE(*,*) 'before ', ecc,peri,gei,alpha
    conv_contr=10*epsilon(1.d0)
    CALL solve_peri(r0,sig0,peri,gm(iplam),alpha,psi,dt,conv_contr)
    IF(nv.eq.3)THEN
       CALL fser_propag(xt,vt,0.d0,dt,gm(iplam),xat(1:3),vat(1:3))
    ELSE
       nd=nv/3-1
       CALL vawrxv(xt,vt,dxdx0(:,1:nd),nv,nd)
       CALL fser_propag_der(xt,vt,0.d0,dt,gm(iplam),xat(1:3),vat(1:3),dxvdxtvt)
       dxvdx0(:,1:nd)=MATMUL(dxvdxtvt,dxdx0(:,1:nd))
       CALL varunw(dxvdx0,xat,vat,nd,nv)
       r0=vsize(xat)
       vel2=prscal(vat,vat)
       call prvec(xat,vat,ang)
       gei2=prscal(ang,ang)
       gei=sqrt(gei2)
       IF(gei.eq.0.d0)THEN 
          WRITE(*,*) ' falsi: zero angular momentum ',xat,vat,t0+dt 
       ENDIF
       call prvec(vat,ang,vlenz)                    ! Lenz vector
       vlenz=vlenz/gm(iplam)-xat(1:3)/r0
       ecc=vsize(vlenz)
       peri=gei2/(gm(iplam)*(1.d0+ecc))
       alpha=vel2-2*gm(iplam)/r0 ! 2* energy
!       WRITE(*,*) 'after ', ecc,peri,gei,alpha,gm(iplam)
    ENDIF
    CALL earcar(t0+dt,xpla,1)
END SUBROUTINE tp_fser
! ===================================================                   
! STRCLO                                                                
! store close approach                                                  
! ===================================================                   
SUBROUTINE strclo(iplam0,tcur,xpla,xa,va,nv,jc,r,rdot,            &
     &             dmintil,c1min,nummin) 
  USE tp_trace
  USE planet_masses
  USE dyn_param, ONLY:nls
! input                                                                 
! planet number
  INTEGER, INTENT(IN) :: iplam0
! time current, distance, radial velocity (should be small)             
  DOUBLE PRECISION,INTENT(IN) :: tcur,r,rdot 
! cartesian coordinates of planet at local MOID                           
  DOUBLE PRECISION, INTENT(IN) :: c1min(3,nummin) 
! SQUARED DISTANCE function at local MOIDs
  DOUBLE PRECISION, INTENT(IN):: dmintil(nummin) 
! number of relative minima (local MOID) found
  INTEGER, INTENT(IN) :: nummin 
!  planet position and velocity                                         
  DOUBLE PRECISION, INTENT(IN):: xpla(6) 
! heliocentrci position and velocity                                    
!  with partial derivatives if nv=21, without if nv=3                   
  INTEGER, INTENT(IN) :: nv 
  DOUBLE PRECISION, INTENT(IN) :: xa(nv),va(nv) 
!  input/output : close approach minima conter                          
  INTEGER, INTENT(INOUT) :: jc 
! =========end interface==================== 
! state vectors: without partials ! with partials 
  INCLUDE 'nvarx.h90' 
  DOUBLE PRECISION x(nvar2x),v(nvar2x)
  DOUBLE PRECISION xat(nvar2x),vat(nvar2x) 
  CHARACTER*30 :: planam 
  DOUBLE PRECISION vl,vsize,dt,r2,rdot2,tcur2, xpla2(6)
  INTEGER i 
  INTEGER lpla,lench 
! calendar date variables                                               
  INTEGER iyear,imonth,iday 
  DOUBLE PRECISION hour 
  CHARACTER*16 date 
! multiple minima analysis                                              
  DOUBLE PRECISION cosa,angle,prscal,angmin 
  INTEGER iun
  INTEGER j
  iun=abs(iun_log)
  iplam=iplam0
  planam=ordnam(iplam)
! store current local moid (just before encounter)                      
  IF(nummin.le.0)THEN
     moid0=-0.9999d0
     angmoid=0.d0
  ELSE   
     moid0=abs(dmintil(1)) 
     angmoid=4.d2 
     angmin=4.d2
! control to avoid computation of moid in strange orbits 
     IF(nv.ge.21)THEN
        DO j=1,nummin 
           cosa=prscal(c1min(1,j),xpla)/(vsize(c1min(1,j))*vsize(xpla)) 
           angle=acos(cosa) 
           IF(angle.lt.angx*radeg)THEN
              IF(angle.lt.angmin*radeg)THEN 
                 moid0=abs(dmintil(j)) 
                 angmoid=angle*degrad
                 angmin=angle*degrad 
              ENDIF
           ENDIF
        ENDDO
     ENDIF
  ENDIF
! planetocentric position with derivatives
  x(1:nv)=xa(1:nv)
  v(1:nv)=va(1:nv)
  x(1:3)=xa(1:3)-xpla(1:3) 
  v(1:3)=va(1:3)-xpla(4:6) 
! handle near miss
! WARNING: why 20*reau? should be less???????????????????????????
   IF(iplam.eq.3.and.r.lt.reau)THEN
! WARNING ????????????????????????????????????????
! if inside Earth, use 2-body propagation
      CALL tp_fser(nv,tcur,x,v,iplam,xpla2,xat,vat,dt)
! xat,vat are planetocentric a periplanet
      r2=vsize(xat(1:3)) 
! minimum forced by 2-body solution, store
      rdot2=prscal(xat(1:3),vat(1:3))/r2
      tcur2=tcur+dt
! only for resret, if passing inside Earth, returns are not real
      IF(fix_mole.and.r.lt.reau) kill_propag=.true.
! the following setting doesn't work as expected
!     IF(fix_mole.and.r.lt.reau/3.d0) kill_propag=.true. 
   ELSE
      xat(1:nv)=x(1:nv)
      vat(1:nv)=v(1:nv)
      r2=r
      rdot2=rdot
      tcur2=tcur
      xpla2=xpla
   ENDIF    
! store close approach point                                            
  jc=jc+1 
  IF(jc.gt.njcx)THEN 
     WRITE(*,*)' strclo: jc>njcx ',jc,njcx,' at ', tcur
     WRITE(ierrou,*)' strclo: jc>njcx ',jc,njcx,' at ', tcur
     numerr=numerr+1
     IF(fix_mole)THEN
        kill_propag=.true.
     ENDIF
     jc=jc-1 
     RETURN
  ENDIF
  tcla(jc)=tcur2 
  rmin(jc)=r2 
  DO i=1,3 
     xcla(i,jc)=xat(i)
     vcla(i,jc)=vat(i) 
     xplaj(i,jc)=xpla(i) 
     xplaj(i+3,jc)=xpla(3+i) 
  ENDDO
  lpla=lench(planam) 
  numcla=numcla+1 
  call mjddat(tcur2,iday,imonth,iyear,hour) 
  IF(verb_clo.gt.9)THEN
     write(date,'(i4,a1,i2.2,a1,i2.2,f6.5)') iyear,'/',imonth,'/',iday,hour/24d0
     WRITE(*,97)planam(1:lpla),date,tcur2,r2 
     WRITE(iun,97)planam(1:lpla),date,tcur2,r2 
97   FORMAT(' Close approach to ',a,' on ',a16,f12.5,' MJD at ',f10.8,' au')
  ENDIF
  IF(nv.eq.3)THEN 
! nothing to do
  ELSEIF(nv.ge.21)THEN
     IF(nv.ne.21+3*nls)THEN
        WRITE(*,*)'strclo: error in dimensions, nv, nls ', nv,nls
        STOP
     ENDIF
! store partials at close approach time                                 
     DO i=4,nv 
        xcla(i,jc)=xat(i) 
        vcla(i,jc)=vat(i) 
     ENDDO
! planet velocity is needed [to define reference system] 
     xplaj(1:6,jc)=xpla(1:6) 
  ELSE 
     write(*,*)' strclo: nv=',nv 
     STOP 
  ENDIF
END SUBROUTINE strclo
  
 ! =================================================
 ! FALSI_SURF
 ! regula falsi for impact corridor computation
 ! =================================================
 SUBROUTINE falsi_surf(t,xa,va,nv,xpla,first,iplam) 
   USE output_control
   USE runge_kutta_gauss
   USE planet_masses
   USE critical_points
   USE orbit_elements
   USE surf_trace, ONLY: fail_surf
   ! current position and velocity of asteroid and planet                                         
   INTEGER,          INTENT(IN) :: nv 
   REAL(KIND=dkind), INTENT(IN) :: t,xa(nv),va(nv),xpla(6) 
   ! is this the beginning of an integration?                              
   LOGICAL,          INTENT(IN) :: first 
   ! planet being approached                                               
   INTEGER,          INTENT(IN) :: iplam 

   ! fixed stepsize time interval                                          
   REAL(KIND=dkind)   :: tt1,tt2 
   ! data for regula falsi        
   REAL(KIND=dkind)   :: ell_dist1,ell_dist2,ell_dist0,t1,t2,t0,hei
   ! state vectors: with partials 
   INCLUDE 'nvarx.h90' 
   REAL(KIND=dkind)   :: x(nvar2x),v(nvar2x),xt(nvar2x),vt(nvar2x)
   REAL(KIND=dkind)   :: xat(nvar2x),vat(nvar2x),xplat(6) 
   ! time steps                                                            
   REAL(KIND=dkind)   :: dt,tt,di,hh
   ! loop indices 
   INTEGER            :: it,i  
   ! functions                                                             
   !REAL(KIND=dkind)   :: vsize,prscal
   ! max no. iterations of falsi 
   INTEGER, PARAMETER :: itmax_f=50 
   
   ! memory model is static                                                
   SAVE 
   ! default for flag stopping propagation due to collision
   kill_propag=.false.
   ! planetocentric position with derivatives
   x(1:nv)=xa(1:nv)
   v(1:nv)=va(1:nv)
   x(1:3)=xa(1:3)-xpla(1:3) 
   v(1:3)=va(1:3)-xpla(4:6) 
   ! initialisation at the beginning of each close approach   
  
   IF(first)THEN
      fail_surf=.FALSE.
      tt1=t 
      ell_dist1=geod_hsurf(iplam,x(1:3))
      RETURN 
   ENDIF
   ! compute distance to ellipsoid surface           
   ell_dist2=geod_hsurf(iplam,x(1:3))
   tt2=t 
   ! direction of time propagation                                         
   IF(tt2.gt.tt1)THEN 
      di=1.d0 
   ELSEIF(tt2.lt.tt1)THEN 
      di=-1.d0 
   ELSE 
      WRITE(iun_log,199) tt1,tt2,x(1:3)
199   FORMAT('falsi_surf: zero step ',f8.2,1x,f8.2,1P,3(1x,d10.3))
      RETURN 
   ENDIF
   IF(abs(ell_dist2).lt.eprdot)THEN 
      ! already at ellipsoid surface: store result
      CALL strsurf(iplam,t,xpla,xa,va,nv)                                         
      ell_dist2=0.d0
   ELSEIF(ell_dist1*di.gt.0.d0.and.ell_dist2*di.lt.0.d0)THEN 
      ! ell_dist changes sign, in a way appropriate for an object approaching the ellipsoid surface         
      t1=tt1 
      t2=tt2 
      dt=-ell_dist2*(t1-t2)/(ell_dist1-ell_dist2)
      tt=dt+t2 
      hh=tt-t 
      ! iterate regula falsi                                                  
      DO it=1,itmax_f 
         CALL rkg(t,xa,va,nv,hh,xat,vat,xplat)
         ! planetocentric position                                        
         xt(1:3)=xat(1:3)-xplat(1:3) 
         vt(1:3)=vat(1:3)-xplat(4:6)
         ell_dist0=geod_hsurf(iplam,xt(1:3))
         t0=tt 
         ! selection of next couple of points with opposite sign                 
         IF(abs(ell_dist0).lt.eprdot)THEN 
            ! already at ellipsoid surface: store result 
            CALL strsurf(iplam,tt,xplat,xat,vat,nv)     
            GOTO 2 
         ELSEIF(ell_dist0*ell_dist1.lt.0.d0)THEN 
            ell_dist2=ell_dist0
            t2=t0 
         ELSEIF(ell_dist0*ell_dist2.lt.0.d0)THEN 
            ell_dist1=ell_dist0
            t1=t0 
         ENDIF
         dt=-ell_dist2*(t1-t2)/(ell_dist1-ell_dist2)
         tt=dt+t2 
         hh=tt-t 
      ENDDO
      ! failed convergence              
      CALL strsurf(iplam,t0,xplat,xat,vat,nv)
      WRITE(ierrou,*) 'falsi_surf: failed convergence for ', ordnam(iplam) 
      WRITE(ierrou,*) 't1,ell_dist1,t2,ell_dist2,dt,tt' 
      WRITE(ierrou,111) t1,ell_dist1,t2,ell_dist2,dt,tt
111   FORMAT(6(1x,E24.16)) 
      numerr=numerr+1 
   ELSEIF(ell_dist1*di.lt.0.d0.and.ell_dist2*di.gt.0.d0)THEN
      numerr=numerr+1
      WRITE(ierrou,*) 'falsi_surf: object comes from interior of surface!'
      WRITE(*,*) 'falsi_surf: object comes from interior of surface!'
      kill_propag=.TRUE.
      fail_surf=.TRUE.
      RETURN
      !STOP
   ENDIF
   ! found zero                                                            
2  CONTINUE 
   ! save current state to be previous state next time                     
   tt1=tt2 
   ell_dist1=ell_dist2
 END SUBROUTINE falsi_surf

 ! ===================================================                   
 ! GEOD_HSURF
 ! It computes the difference between the altitude over the 
 ! geoid surface of the planet and the altitude height_surf
 ! Zero value corresponds to a point at height height_surf 
 ! over the geoid surface of the planet.
 ! ===================================================                   
 REAL(KIND=dkind) FUNCTION geod_hsurf(iplam0,pos_pla_ecl)
   USE surf_trace
   INTEGER,                        INTENT(IN) :: iplam0       ! index of planet
   REAL(KIND=dkind), DIMENSION(3), INTENT(IN) :: pos_pla_ecl  ! geocentric ecliptic position
   REAL(KIND=dkind), DIMENSION(3) :: x                        ! geocentric equatorial position
   INTEGER                        :: i                        ! loop index
   REAL(KIND=dkind)               :: ecc_squared              ! eccentricity of geoid squared
   REAL(KIND=dkind)               :: z,t,zt,sinfi,h,dn,dnh,delth,h0

   IF(iplam0.eq.3)THEN
      ! this is an approximation, the altitude computed will not be exactly the altitude on Earth
      x=MATMUL(roteceq,pos_pla_ecl)  
   ELSE
      numerr=numerr+1
      WRITE(ierrou,*) 'geod_hsurf: not ready for impact corridor on planet different from Earth'
      WRITE(*,*) 'geod_hsurf: not ready for impact corridor on planet different from Earth'
      STOP
   ENDIF
   ecc_squared=flatt_surf*(2.d0-flatt_surf)
   ! initialization
   t=0.d0
   h0=0.d0
   dn=0.d0
   ! conversion to km for internal loop
   x=x*aukm
   z=x(3)
   DO i=1,40
      zt=z+t
      dnh=SQRT((x(1)**2+x(2)**2+zt**2))
      sinfi=zt/dnh
      !dn=eradkm/SQRT(1.d0-ecc_squared*sinfi**2)
      dn=rad_surf/SQRT(1.d0-ecc_squared*sinfi**2)
      t=dn*ecc_squared*sinfi
      h=dnh-dn
      delth=ABS(h-h0)
      IF(delth.lt.1.d-5) EXIT
      h0=h
   ENDDO
   IF(i.ge.40)THEN
      numwar=numwar+1
      WRITE(iwarou,*) 'geod_hsurf: convergence not attained, i= ',i
   ENDIF
   geod_hsurf=(h-height_surf*aukm)/aukm
 END FUNCTION geod_hsurf

 ! ===================================================                   
 ! STRSURF                                                             
 ! store intersection point with surface at height h
 ! ===================================================                   
 SUBROUTINE strsurf(ipl,tcur,xpla,xa,va,nv) 
   USE surf_trace
   USE dyn_param, ONLY:nls
   ! input             
   INTEGER, INTENT(IN) :: ipl   ! planet index 
   ! time current
   DOUBLE PRECISION,INTENT(IN) :: tcur
   !  planet position and velocity                                         
   DOUBLE PRECISION, INTENT(IN):: xpla(6) 
   ! heliocentrci position and velocity                                    
   !  with partial derivatives if nv=21, without if nv=3
   !  (also dynamical paramenters if nv>21)                   
   INTEGER, INTENT(IN) :: nv 
   DOUBLE PRECISION, INTENT(IN) :: xa(nv),va(nv) 
   ! =========end interface==================== 
   ! state vectors: without partials ! with partials 
   INCLUDE 'nvarx.h90' 
   DOUBLE PRECISION x(nvar2x),v(nvar2x)
   INTEGER i 

   ! planetocentric position with derivatives
   x(1:nv)=xa(1:nv)
   v(1:nv)=va(1:nv)
   x(1:3)=xa(1:3)-xpla(1:3) 
   v(1:3)=va(1:3)-xpla(4:6) 

   ! store surface intersection point                                            
   tsurf=tcur 
   ipla_ic=ipl   ! planet index   
   DO i=1,3
      xsurf(i)=x(i)
      vsurf(i)=v(i) 
      xplajs(i)=xpla(i) 
      xplajs(i+3)=xpla(3+i) 
   ENDDO
      
   IF(nv.GE.21)THEN
      IF(nv.ne.21+3*nls)THEN
         numerr=numerr+1
         WRITE(ierrou,*)'strsurf: error in dimensions, nv, nls ', nv,nls
         WRITE(*,*)'strsurf: error in dimensions, nv, nls ', nv,nls
         STOP
      ENDIF
      ! store partials at close approach time                                 
      DO i=4,nv 
         xsurf(i)=x(i) 
         vsurf(i)=v(i) 
      ENDDO
      ! planet velocity is needed [to define reference system] 
      xplajs(1:6)=xpla(1:6) 
   ELSEIF(nv.NE.3)THEN
      numerr=numerr+1
      WRITE(ierrou,*)' strsurf: nv=',nv 
      WRITE(*,*)' strsurf: nv=',nv 
      STOP 
   ENDIF
   ! STOP PROPAGATION
   kill_propag=.TRUE.
 END SUBROUTINE strsurf


END MODULE close_app  

! routines not belonging to the module
! ==============================================================  
! NPOINT_SET
! ==============================================================  
SUBROUTINE npoint_set(npmult)
  USE close_app
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: npmult
  INTEGER npoint0
  DOUBLE PRECISION eprdot0
  LOGICAL :: init = .true.
  IF(init)THEN
     npoint0=npoint
     eprdot0=eprdot
     init=.false.
  ENDIF
  npoint=npoint0*npmult
  IF(npmult.eq.1)THEN
     eprdot=eprdot0
  ELSE
     eprdot=1d-3*eprdot0
  ENDIF
END SUBROUTINE npoint_set
! ==============================================================  
! SET_CLOST
SUBROUTINE set_clost(clostin)
  USE close_app
  IMPLICIT NONE
  LOGICAL clostin
  clost=clostin
END SUBROUTINE set_clost
! ==============================================================        
! COV_AVAI/COV_NOT_AV                                                   
! routines to manipulate the covariance common used                     
! for online target plane analysys                                      
! =============================================================         
SUBROUTINE cov_avai(unc,coo,coord) 
  USE tp_trace
  USE orbit_elements
  IMPLICIT NONE 
! INPUT: covariance matrix, elements, coordinates
  TYPE(orb_uncert), INTENT(IN) :: unc
  DOUBLE PRECISION, INTENT(IN) :: coord(6)
  CHARACTER*3, INTENT(IN) :: coo
! HIDDEN OUTPUT: covariance matrix and flag of availability             
  unc_store=unc 
  coord_store=coord
  coo_store=coo
  covava=.true. 
END SUBROUTINE cov_avai
SUBROUTINE cov_not_av 
  USE tp_trace
  IMPLICIT NONE 
! HIDDEN OUTPUT: flag of availability                                   
  covava=.false. 
END SUBROUTINE cov_not_av



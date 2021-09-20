! ============================================================!
! MODULE surf_trace                                           !
! contains the routines and public variables to store impact  !
! at a fixed height over the geoid surface of the planet      !
! being impacted.                                             !
! ============================================================!
!                                                             !
! author:  L. Dimare, D. Bracali Cioci                        !
! date of creation: 9 November 2015                           !
! date of last update: 25 November 2016                       !
! ============================================================!
MODULE surf_trace
  USE fund_const
  USE output_control
  USE dyn_param, ONLY: ndyx,ndimx
  IMPLICIT NONE
  PRIVATE
  
  ! height of the ellipsoid surface wrt the planet surface (AU)
  REAL(KIND=dkind), PUBLIC             :: height_surf
  ! flattening and semi-major axis of the planet reference ellipsoid
  ! rad_surf is in KM
  REAL(KIND=dkind), PUBLIC             :: flatt_surf, rad_surf

  ! time current (ET)
  REAL(KIND=dkind), PUBLIC :: tsurf
  !  planet position and velocity                                         
  REAL(KIND=dkind), PUBLIC :: xplajs(6) 
  ! asteroid planetocentric pos/vel
  REAL(KIND=dkind), PUBLIC :: xsurf(21+3*ndyx),vsurf(21+3*ndyx) 
  ! currently (last) approached planet
  INTEGER, PUBLIC :: ipla_ic
  ! failure of falsi_surf: initial conditions are inside Earth
  LOGICAL, PUBLIC :: fail_surf 

  ! surface point data
  TYPE surf_point
     ! time of intersection with surface (ET)
     REAL(KIND=dkind) :: t         
     ! time of intersection with surface (UTC)
     REAL(KIND=dkind) :: tutc         
     ! planetocentric ecliptic pos/vel at intersection with surface
     REAL(KIND=dkind) :: pgeo(3),vgeo(3)
     ! heliocentric ecliptic pos/vel of planet
     REAL(KIND=dkind) :: pear(3),vear(3)
     ! geodetic (body-fixed) planetocentric coordinates of intersection point (rad)
     REAL(KIND=dkind) :: long,lat,h
     ! planetocentric body-fixed pos/vel at intersection with surface (au, au/d)
     REAL(KIND=dkind) :: pbf(3),vbf(3)
     INTEGER          :: iplam   ! planet index
     ! partials of pos/vel wrt initial parameters 
     REAL(KIND=dkind) :: dxdx0(6,ndimx)
     ! partials of impact time wrt initial parameters
     REAL(KIND=dkind) :: dtdx0(ndimx)
     ! differential of the map from initial elements to surface
     REAL(KIND=dkind) :: dsurfdele(2,ndimx)
     ! differential of the map from initial elements to body-fixed surface
     REAL(KIND=dkind) :: dsurfbfdele(2,ndimx)
  END TYPE surf_point
    
  TYPE(surf_point), PUBLIC :: pt_surf
  
  ! LIST OF PUBLIC ENTITIES
  ! derived TYPEs
  PUBLIC surf_point
  ! subroutines
  PUBLIC undefined_surf_point, str_surfan, surf_geodetic, unc_semiax, dist_geoid

CONTAINS
  SUBROUTINE undefined_surf_point(sfpoint)
    TYPE(surf_point), INTENT(OUT) :: sfpoint

    sfpoint%t=0.d0
    sfpoint%tutc=0.d0
    sfpoint%pgeo=0.d0
    sfpoint%vgeo=0.d0
    sfpoint%pear=0.d0
    sfpoint%vear=0.d0
    sfpoint%long=0.d0
    sfpoint%lat=0.d0
    sfpoint%h=0.d0
    sfpoint%pbf=0.d0
    sfpoint%vbf=0.d0
    sfpoint%iplam=0.d0
    sfpoint%dxdx0=0.d0
    sfpoint%dtdx0=0.d0
    sfpoint%dsurfdele=0.d0
    sfpoint%dsurfbfdele=0.d0

  END SUBROUTINE undefined_surf_point

  ! =============================================================!
  ! Store information on IC surface intersection point and       !
  ! tangent plane                                                !
  !                                                              !
  ! author: L. Dimare                                            !
  ! date of creation: 05 February 2016                           !
  ! date of last update: 30 November 2016                        !
  ! =============================================================!
  SUBROUTINE str_surfan(dx1dx0,nd,dx0de)
    USE planet_masses
    USE iers_ser, ONLY: rotpv
    INTEGER,          INTENT(IN) :: nd           ! number of parameters
    REAL(KIND=dkind), INTENT(IN) :: dx1dx0(6,nd) ! accumulated state transition matrix, for initial conditions
                                                 ! and for dynamical parametrs
    REAL(KIND=dkind), INTENT(IN) :: dx0de(nd,nd) ! derivatives of initial state wrt parameters, including dyn.parameters
    TYPE(surf_point)  :: sfpoint
    ! planet names                                                          
    CHARACTER(LEN=30) :: planam 
    INTEGER           :: lpla
    ! calendar date variables                                               
    INTEGER           :: iyear,imonth,iday 
    REAL(KIND=dkind)  :: hour 
    CHARACTER(LEN=16) :: date 
    ! number of propagated equations (pos+partials of pos)
    INTEGER           :: nvar2   
    ! state transition matrix including dyn. parameters
    REAL(KIND=dkind)  :: dx1pdx0p(nd,nd)
    ! derivatives with respect state vector (incond+dyn.par)
    REAL(KIND=dkind)  :: dxdx1(6,nd),dxdx0(6,nd)
    ! derivatives with respect state vector (incond+dyn.par)
    REAL(KIND=dkind)  :: dxde(6,nd)
    ! differential of function to impact surface
    REAL(KIND=dkind)  :: dadtde(3,nd),dadde(2,nd)
    ! equatorial body-fixed coordinates
    REAL(KIND=dkind)  :: deq(3),veq(3)
    ! geodetic coordinates
    REAL(KIND=dkind)  :: long,lat,alt
    ! time: days,seconds
    INTEGER           :: timed,timed2
    REAL(KIND=dkind)  :: timesec,timesec2
    INTEGER           :: j         ! loop index
    ! functions
    INTEGER           :: lench 

    CALL undefined_surf_point(sfpoint)
    ! planet
    planam=ordnam(ipla_ic) 
    lpla=lench(planam)
    ! Gregorian calendar date and hour
    CALL mjddat(tsurf,iday,imonth,iyear,hour) 
    WRITE(date,'(i4,a1,i2.2,a1,i2.2,f6.5)') iyear,'/',imonth,'/',iday,hour/24d0 
    ! new state transition matrix including dynamical parameters 
    dx1pdx0p(1:6,1:nd)=dx1dx0
    IF(nd.gt.6)THEN
       dx1pdx0p(7:nd,1:nd)=0.d0
       DO j=7,nd
          dx1pdx0p(j,j)=1.d0
       ENDDO
    ENDIF
!!! Do we want to store information also for boundary???
!!! No. It is enough to store only the center of IC with derivatives. 
!!! In this case we know that nvar>6, because we are propagating the variational eq. 
    sfpoint%t=tsurf
    sfpoint%pgeo=xsurf(1:3)
    sfpoint%vgeo=vsurf(1:3)
    sfpoint%pear=xplajs(1:3)
    sfpoint%vear=xplajs(4:6)
    sfpoint%iplam=ipla_ic
    ! WARNING: the state transition matrix stored in xsurf and vsurf is the result
    ! of the RKG propagation inside the routine falsi_surf (from current state inside RA15 
    ! to x=state at intersection time). 
    ! Then it does not give the global state transition matrix and must be multiplied by dx1dx0,
    ! obtained at a previous call of RA15. It gives dxdx1. 
    ! Actually for IC computation, RA15 is called once and the propagation stops with kill-propag
    ! after the close approach, so that dx1dx0 should be the identity, unless the propagation
    ! is a continuation of a previous one.
    nvar2=3*(nd+1)
    ! partial state transition matrix 6xnd: dxdx1
    ! extrapolate matrix dxdx1 from vectors xsurf, vsurf 
    ! x1=initial state of current propagation
    CALL vawrxv(xsurf,vsurf,dxdx1,nvar2,nd) 
    ! Chain rule: compute dxde=dxdx1*dx1pdx0p*dx0de (6xnd)=(6xnd)x(ndxnd)x(ndxnd)
    dxdx0=MATMUL(dxdx1,dx1pdx0p)
    dxde=MATMUL(dxdx0,dx0de)
    ! differential of function to IC surface:
    ! WARNING: IC surface is not defined in the body-fixed frame, 
    ! it is in 'MEAN' not in 'BF  '
    dadtde(1:3,1:nd)=diff_icfun(dxde(1:6,1:nd),nd)
    sfpoint%dxdx0(1:6,1:nd)=dxde
    sfpoint%dsurfdele(1:2,1:nd)=dadtde(1:2,1:nd)
    sfpoint%dtdx0(1:nd)=dadtde(3,1:nd)

    ! ADDITIONAL INFORMATION IN BODY-FIXED REFERENCE
    ! rotation to the equatorial body-fixed reference system   
    dadde(1:2,1:nd)=diff_icfun2(dxde(1:6,1:nd),nd)
    sfpoint%dsurfbfdele(1:2,1:nd)=dadde
    timed=FLOOR(tsurf)
    timesec=(tsurf-timed)*86400.d0
    CALL rotpv('ECLM',.true.,51544,43200.d0,xsurf(1:3),vsurf(1:3),'BF  ',.true.,timed,timesec,deq,veq)
    ! geodetic coordinates 
    CALL surf_geodetic(deq,long,lat,alt)
    sfpoint%long=long
    sfpoint%lat=lat
    sfpoint%h=alt
    sfpoint%pbf=deq
    sfpoint%vbf=veq
    ! convert time to UTC
    CALL cnvtim(timed,timesec,'ET ',timed2,timesec2,'UTC')
    sfpoint%tutc=timed2+timesec2/86400.d0
    ! store information in public variable TYPE(surf_point)
    pt_surf=sfpoint    

  END SUBROUTINE str_surfan

  ! =============================================================!
  ! Compute the derivatives of geodetic angular cordinates       ! 
  ! (long,lat) and impact time wrt the initial state vector,     !
  ! where (long,lat) are the coordinates on the surface          !
  ! at altitude height_surf wrt the geoid.                       !
  ! This function uses public data xsurf,vsurf.                  !
  !                                                              !
  ! author: L. Dimare, D. Bracali Cioci                          !
  ! date of creation: 2 December 2015                            !
  ! date of last update: 30 November 2016                        !
  ! =============================================================!
  FUNCTION diff_icfun(dphi,nd)
    INTEGER, INTENT(IN) :: nd
    REAL(KIND=dkind), DIMENSION(6,nd), INTENT(IN) :: dphi
    REAL(KIND=dkind), DIMENSION(3,nd) :: diff_icfun
    
    REAL(KIND=dkind) :: d(3),deq(3)
    REAL(KIND=dkind) :: alpha,delta,alt,dadx(3),dddx(3),dgdx(3),der(3,3)
    REAL(KIND=dkind) :: dis0
    REAL(KIND=dkind) :: adot,ddot,gdot
    REAL(KIND=dkind), ALLOCATABLE :: dgdx0(:),dadx0(:),dddx0(:),dtdx0(:)

    ALLOCATE(dgdx0(nd),dadx0(nd),dddx0(nd),dtdx0(nd))

    diff_icfun=0.d0
    ! Planetocentric position
    d=xsurf(1:3)     
    ! rotation to the equatorial reference system   
    deq=MATMUL(roteceq,d)
    CALL surf_geodetic(deq,alpha,delta,alt,DER=der)
    dadx=der(1,1:3)
    dddx=der(2,1:3)
    dgdx=der(3,1:3)
    ! Come back to Ecliptical
    dadx=MATMUL(dadx,roteceq)
    dddx=MATMUL(dddx,roteceq)
    dgdx=MATMUL(dgdx,roteceq)
    ! Computation of time derivatives of (long,lat,alt)
    adot=DOT_PRODUCT(dadx,vsurf(1:3))
    ddot=DOT_PRODUCT(dddx,vsurf(1:3))
    gdot=DOT_PRODUCT(dgdx,vsurf(1:3))

    dgdx0=MATMUL(dgdx,dphi(1:3,1:nd))

    dtdx0=-dgdx0/gdot
    dadx0=MATMUL(dadx,dphi(1:3,1:nd))+adot*dtdx0
    dddx0=MATMUL(dddx,dphi(1:3,1:nd))+ddot*dtdx0

    diff_icfun(1,1:nd)=dadx0
    diff_icfun(2,1:nd)=dddx0
    diff_icfun(3,1:nd)=dtdx0

    DEALLOCATE(dgdx0,dadx0,dddx0,dtdx0)
  END FUNCTION  diff_icfun

  ! =============================================================!
  ! Compute the derivatives of geodetic angular cordinates       ! 
  ! (long,lat) wrt the initial state vector,                     !
  ! where (long,lat) are the coordinates on the surface          !
  ! at altitude height_surf wrt the geoid, in the body-fixed     ! 
  ! reference frame.                                             !
  ! This function uses public data tsurf,xsurf,vsurf.            !
  !                                                              !
  ! author: L. Dimare                                            !
  ! date of creation: 14 November 2016                           !
  ! date of last update: 14 November 2016                        !
  ! =============================================================!
  FUNCTION diff_icfun2(dphi,nd)
    USE iers_ser, ONLY: rotpv,rotsys,mj2000,s2000
    INTEGER, INTENT(IN) :: nd
    REAL(KIND=dkind), DIMENSION(6,nd), INTENT(IN) :: dphi
    REAL(KIND=dkind), DIMENSION(2,nd) :: diff_icfun2
    
    REAL(KIND=dkind) :: d(3),dv(3),deq(3),veq(3)
    REAL(KIND=dkind) :: alpha,delta,alt,dadx(3),dddx(3),dgdx(3),der(3,3)
    REAL(KIND=dkind) :: adot,ddot,gdot
    REAL(KIND=dkind), ALLOCATABLE :: dgdx0(:),dadx0(:),dddx0(:),dtdx0(:)
    REAL(KIND=dkind) :: rot(3,3),rot1(3,3),rot2(3,3)
    INTEGER          :: timed       ! days of time of surface intersection (ET)
    REAL(KIND=dkind) :: timesec     ! seconds of time of surface intersection (ET)
    INTEGER          :: i

    ALLOCATE(dgdx0(nd),dadx0(nd),dddx0(nd),dtdx0(nd))
    
    diff_icfun2=0.d0
    
    ! DERIVATIVE OF IMPACT TIME (if falsi uses roteceq, i.e. the 
    ! impact surface is defined in the equatorial J2000 reference)
    d=xsurf(1:3)            ! Planetocentric position 
    deq=MATMUL(roteceq,d)
    CALL surf_geodetic(deq,alpha,delta,alt,DER=der)
    dgdx=MATMUL(der(3,1:3),roteceq)
    dgdx0=MATMUL(dgdx,dphi(1:3,1:nd))
    gdot=DOT_PRODUCT(dgdx,vsurf(1:3))
    dtdx0=-dgdx0/gdot
   
    ! rotation to the body-fixed reference system
    timed=FLOOR(tsurf)
    timesec=(tsurf-timed)*86400.d0
    ! transformation of velocity from AU/d to AU/s
    dv=vsurf(1:3)/86400.d0
    CALL rotpv('ECLM',.true.,mj2000,s2000,d,dv,'BF  ',.true.,timed,timesec,deq,veq)
    ! transformation of velocity from AU/s to AU/d
    veq=veq*86400.d0
    CALL rotsys('ECLM',mj2000,s2000,'BF  ',timed,timesec,rot,rot1,rot2,0)
    CALL surf_geodetic(deq,alpha,delta,alt,DER=der)
    dadx=der(1,1:3)
    dddx=der(2,1:3)
    ! Come back to Ecliptical
    dadx=MATMUL(dadx,rot)
    dddx=MATMUL(dddx,rot)
    
    dadx0(1:nd)=MATMUL(dadx,dphi(1:3,1:nd))+DOT_PRODUCT(der(1,1:3),veq)*dtdx0
    dddx0(1:nd)=MATMUL(dddx,dphi(1:3,1:nd))+DOT_PRODUCT(der(2,1:3),veq)*dtdx0
    
    diff_icfun2(1,1:nd)=dadx0
    diff_icfun2(2,1:nd)=dddx0
        
    DEALLOCATE(dgdx0,dadx0,dddx0,dtdx0)
  END FUNCTION  diff_icfun2

  ! =============================================================!
  ! Compute geodetic angular coordinates of a given point on the !
  ! IC surface                                                   !
  !                                                              !
  ! author: L. Dimare, D. Bracali Cioci                          !
  ! date of creation: 10 December 2015                           !
  ! date of last update: 5 February 2016                         !
  ! =============================================================!
  SUBROUTINE surf_geodetic(pos_car,long,lat,alt,vert,der)
    REAL(KIND=dkind), DIMENSION(3),   INTENT(IN)            :: pos_car     ! cartesian coordinates 
    REAL(KIND=dkind),                 INTENT(OUT)           :: long        ! geodetic longitude (positive east, rad)
    REAL(KIND=dkind),                 INTENT(OUT)           :: lat         ! geodetic latitude (rad)
    REAL(KIND=dkind),                 INTENT(OUT)           :: alt         ! altitude on ellipsoid 
    REAL(KIND=dkind), DIMENSION(3),   INTENT(OUT), OPTIONAL :: vert        ! local vertical direction
    REAL(KIND=dkind), DIMENSION(3,3), INTENT(OUT), OPTIONAL :: der         ! matrix of derivatives wrt x
    REAL(KIND=dkind) :: ecc_squared
    REAL(KIND=dkind) :: x(3)                                               ! cartesian coordinates 
    REAL(KIND=dkind) :: t,h,h0,dn,dnh,z,zt,delth,sinfi,cosfi,der_aux(3,3)
    REAL(KIND=dkind) :: det,dd(9),mat(9)
    INTEGER          :: i,j,ising

    ecc_squared=flatt_surf*(2.d0-flatt_surf)
   
    t=0.d0
    h0=0.d0
    dn=0.d0
    ! conversion to km
    x=pos_car*aukm
    z=x(3)

    DO i=1,40
       zt=z+t
       dnh=SQRT((x(1)**2+x(2)**2+zt**2))
       sinfi=zt/dnh
       dn=rad_surf/SQRT(1.d0-ecc_squared*sinfi**2)
!       dn=eradkm/SQRT(1.d0-ecc_squared*sinfi**2)
       t=dn*ecc_squared*sinfi
       h=dnh-dn
       delth=ABS(h-h0)
       IF(delth .lt. 1.d-5) EXIT
       h0=h
    ENDDO
    IF (i.ge.40) THEN
       WRITE(ierrou,*) 'surf_geodetic: not converged; DELTH= ',delth
       numerr=numerr+1
    ENDIF
    alt=h   ! altitude in km
    IF(sinfi.gt.1.d0)THEN
       IF(sinfi.lt.1.d0+10.d0*epsilon(1.d0))THEN
          sinfi=1.d0
       ELSE
          WRITE(ierrou,*)'surf_geodetic: error! sinfi=',sinfi
          numerr=numerr+1
       ENDIF
    ENDIF
    IF(sinfi.lt.-1.d0)THEN
       IF(sinfi.gt.-1.d0-10.d0*epsilon(1.d0))THEN
          sinfi=-1.d0
       ELSE
          WRITE(ierrou,*) 'surf_geodetic: error! sinfi=',sinfi 
          numerr=numerr+1
       ENDIF
    ENDIF
    lat=ASIN(sinfi)
    long=ATAN2(x(2),x(1))
    IF(long .lt. 0.d0) long=dpig+long
    IF(PRESENT(vert)) THEN
       cosfi=COS(lat)
       vert(1)=cosfi*dcos(long)
       vert(2)=cosfi*dsin(long)
       vert(3)=sinfi
    END IF
    ! Derivatives
    IF(PRESENT(der))THEN
       der_aux(1:3,1)=dnh*COS(lat)*(/-SIN(long),COS(long),0.d0/)
       der_aux(1:3,2)=(dn*(1.d0-ecc_squared)/(1.d0-ecc_squared*sinfi**2)+alt)* &
            & (/-SIN(lat)*COS(long),-SIN(lat)*SIN(long),COS(lat)/)
       der_aux(1:3,3)=(/COS(lat)*COS(long),COS(lat)*SIN(long),SIN(lat)/)
       !CALL matin(der_aux,det,3,3,3,ising,1) ! it is not precise enough
       dd=0.d0   ! initialization of identity matrix
       DO i=1,3
          DO j=1,3
             mat((i-1)*3+j)=der_aux(i,j)
          ENDDO
          dd((i-1)*3+i)=1.d0
       ENDDO
!       CALL dgelg(dd,mat,3,3,1.d-30,ising) 
       CALL dgelg(dd,mat,3,3,1.d-15,ising) 
       IF(ising.NE.0)THEN
          WRITE(*,*) 'surf_geodetic: inversion of der_aux failed. Derivatives not computed'
          STOP
       ENDIF
       DO i=1,3
          DO j=1,3
             der_aux(i,j)=dd((i-1)*3+j)
          ENDDO
       ENDDO
       der=der_aux
       ! Conversion to AU
       der(1:2,1:3)=der(1:2,1:3)*aukm
    END IF
    ! Conversion to AU
    alt=alt/aukm
  END SUBROUTINE surf_geodetic
  
  ! Solve linear equations with Gauss elimination method      !
  ! using total maximum pivot (double precision).             !
  ! n linear systems with the same coefficients matrix A      !
  ! are solved simultaneosly.                                 !
  !                                                           !
  ! See Bini D., Capovani M., Menchi O., Metodi numerici      ! 
  ! per l'algebra lineare,  Zanichelli, Chapter 4.            !
  SUBROUTINE dgelg(r,a,m,n,eps,ier) 
    implicit none 
    
    integer :: i,ier,ii,ist,j,k,lend,l,ll,lst,n,nm,m,mm 
    double precision :: eps,piv,pivi,tb,tol 
    double precision :: a(m*m),r(m*n) 
    
    if(m.le.0) GOTO 23 
    ier=0 
    piv=0.d0 
    mm=m*m 
    nm=n*m 
    DO l=1,mm 
       tb=dabs(a(l)) 
       IF(tb-piv.gt.0.d0)THEN 
          piv=tb 
          i=l
       ENDIF
    ENDDO
    tol=eps*piv 
    lst=1 
    do 17 k=1,m 
       IF(piv.le.0.d0) GOTO 23
       IF(ier.eq.0)THEN
          IF(piv-tol.le.0.d0)THEN
             ier=k-1
          ENDIF
       ENDIF
       pivi=1.d0/a(i) 
       j=(i-1)/m 
       i=i-j*m-k 
       j=j+1-k 
       DO l=k,nm,m 
          ll=l+i 
          tb=pivi*r(ll) 
          r(ll)=r(l) 
          r(l)=tb
       ENDDO
       IF(k-m.ge.0) GOTO 18 
       lend=lst+m-k 
       IF(j.le.0) GOTO 12 
       ii=j*m 
       DO l=lst,lend 
          tb=a(l) 
          ll=l+ii 
          a(l)=a(ll) 
          a(ll)=tb
       ENDDO
12     CONTINUE
       DO l=lst,mm,m 
          ll=l+i 
          tb=pivi*a(ll) 
          a(ll)=a(l) 
          a(l)=tb
       ENDDO
       a(lst)=j 
       piv=0.d0 
       lst=lst+1 
       j=0 
       do 16 ii=lst,lend 
          pivi=-a(ii) 
          ist=ii+m 
          j=j+1 
          DO l=ist,mm,m 
             ll=l-j 
             a(l)=a(l)+pivi*a(ll) 
             tb=dabs(a(l)) 
             IF(tb-piv.gt.0.d0)THEN 
                piv=tb 
                i=l 
             ENDIF
          ENDDO
          DO l=k,nm,m 
             ll=l+j 
             r(ll)=r(ll)+pivi*r(l)
          ENDDO
16     ENDDO
       lst=lst+m
17  ENDDO
18  IF(m-1.lt.0) GOTO 23
    IF(m.eq.1) GOTO 22 
    ist=mm+m 
    lst=m+1 
    DO i=2,m 
       ii=lst-i 
       ist=ist-lst 
       l=ist-m 
       l=a(l)+0.5d0 
       DO j=ii,nm,m 
          tb=r(j) 
          ll=j 
          DO k=ist,mm,m 
             ll=ll+1 
             tb=tb-a(k)*r(ll)
          ENDDO
          k=j+l 
          r(j)=r(k) 
          r(k)=tb 
       ENDDO
    ENDDO
22  return 
23  ier=-1 
    return 
  END SUBROUTINE dgelg

  ! =============================================================!
  ! Compute the lengths of the axes of the linear uncertainty    !
  ! ellipse on the plane tangent to the IC surface at altitude   !
  ! hsurf fixed above reference geoid.                           !
  ! This subroutine uses public data rad_surf and flatt_surf.    !
  !                                                              !
  ! author: L. Dimare                                            !
  ! date of creation: 30 May 2016                                !
  ! date of last update: 27 Oct 2016                             !
  ! =============================================================!
  SUBROUTINE unc_semiax(lat,hsurf,sig,axes,unc_ax,azimuth,metric,spher_unc)
    REAL(KIND=dkind), INTENT(IN)  :: lat                  ! latitude (rad)
    REAL(KIND=dkind), INTENT(IN)  :: hsurf                ! altitude of surface on geoid in km    
    REAL(KIND=dkind), INTENT(IN)  :: sig(2)               ! eigenvalues of covariance on tangent plane
    REAL(KIND=dkind), INTENT(IN)  :: axes(2,2)            ! eigenvectors of covariance matrix
    REAL(KIND=dkind), INTENT(OUT) :: unc_ax(2)            ! lengths of semi-axes (km)
    REAL(KIND=dkind), INTENT(OUT), OPTIONAL :: azimuth    ! semi-major axis azimuth (rad)
    REAL(KIND=dkind), INTENT(OUT), OPTIONAL :: metric(2)  ! coeff. ee,gg of first fund. form (km^2)
    LOGICAL,          INTENT(IN), OPTIONAL  :: spher_unc  ! activate spherical approximation 
   
    REAL(KIND=dkind) :: ecc_squared             ! squared eccentricity of reference geoid
    REAL(KIND=dkind) :: cos_p,sin_p             ! cosine and sine of latitude phi
    REAL(KIND=dkind) :: dn,dnh                  ! auxiliary quantities for geodetic coordinates
    REAL(KIND=dkind) :: ee,gg                   ! metric on tangent plane
    REAL(KIND=dkind) :: vec(2),vv               ! vector on tangent plane and its squared norm 

    ecc_squared=flatt_surf*(2.d0-flatt_surf)
    cos_p=COS(lat)
    sin_p=SIN(lat)
    dn=rad_surf/SQRT(1.d0-ecc_squared*sin_p**2) ! in Km
    dnh=dn+hsurf                                ! in Km  

    ! first fundamental form on geoid/ metric on tangent plane
    ee=(dnh*cos_p)**2
    gg=(dn*(1.d0-ecc_squared)/(1.d0-ecc_squared*sin_p**2)+hsurf)**2

    IF(PRESENT(spher_unc).and.spher_unc)THEN
       ee=((rad_surf+hsurf)*cos_p)**2
       gg=(rad_surf+hsurf)**2
    ENDIF
    
    ! semi-minor axis    
    vec=axes(1:2,1)
    vv=ee*vec(1)**2+gg*vec(2)**2
    unc_ax(1)=SQRT(vv)*sig(1)
    ! semi-major axis
    vec=axes(1:2,2)
    vv=ee*vec(1)**2+gg*vec(2)**2
    unc_ax(2)=SQRT(vv)*sig(2)
    ! semi-major axis azimuth
    IF(axes(1,2).LT.0.d0)THEN
       azimuth=ACOS(-axes(2,2)*SQRT(gg)/SQRT(vv))
    ELSE
       azimuth=ACOS(axes(2,2)*SQRT(gg)/SQRT(vv))
    ENDIF
    IF(PRESENT(metric))THEN
       metric(1)=ee
       metric(2)=gg
    ENDIF

  END SUBROUTINE unc_semiax

  ! =============================================================!
  ! Compute the distance between two points on the geoid.        !
  ! We compute the distance in radians, considering the          ! 
  ! geoid semi-major axis equal to 1.                            !
  ! This subroutine uses the public datum flatt_surf.            !
  !                                                              !
  ! author: L. Dimare                                            !
  ! date of creation: 31 May 2016                                !
  ! date of last update: 25 November 2016                        !
  ! =============================================================!
  SUBROUTINE dist_geoid(pt1,pt2,dist)
    REAL(KIND=dkind), INTENT(IN) :: pt1(2)  ! longitude and latitude of first point (rad)
    REAL(KIND=dkind), INTENT(IN) :: pt2(2)  ! longitude and latitude of second point (rad)
    REAL(KIND=dkind), INTENT(OUT):: dist    ! ellipsoidal distance between points (rad)
    
    INTEGER, PARAMETER :: itmax=50
    REAL(KIND=dkind) :: ecc_squared         ! squared eccentricity of ellipsoid
    REAL(KIND=dkind) :: uu1,uu2             ! reduced longitudes
    REAL(KIND=dkind) :: ldif                ! difference between longitudes
    REAL(KIND=dkind) :: lambda,lpre
    REAL(KIND=dkind) :: snu1,snu2,csu1,csu2 ! sine and cosine of reduced latitudes
    REAL(KIND=dkind) :: snl,csl             ! sine and cosine of lambda
    REAL(KIND=dkind) :: sigma,sn,cs         ! arc length and its sine, cosine
    REAL(KIND=dkind) :: sna,cs2a            ! sine and squared cosine of azimuth
    REAL(KIND=dkind) :: cs_m
    REAL(KIND=dkind) :: aa,bb,cc,k1
    REAL(KIND=dkind) :: u_squared, delsig
    INTEGER          :: i
    ! functions
    REAL(KIND=dkind) :: pridif

    !!! INSERT CONTROL ON ANTIPODAL POINTS !!!

    ecc_squared=flatt_surf*(2.d0-flatt_surf)
    ! First approach: Vincenty's formulae
    ! Reduced latitude 
    uu1=ATAN((1-flatt_surf)*TAN(pt1(2)))
    uu2=ATAN((1-flatt_surf)*TAN(pt2(2)))
    snu1=SIN(uu1)
    csu1=COS(uu1)
    snu2=SIN(uu2)
    csu2=COS(uu2)
    ! Difference in longitude: long2 - long1 in (-pi,pi)
    ldif=pridif(pt2(1),pt1(1))
    lambda=ldif ! initialization
    ! loop to calculate lambda
    DO i=1,itmax
       lpre=lambda
       snl=SIN(lambda)
       csl=COS(lambda)
       ! sigma=arc length between points
       ! sin(sigma) and cos(sigma)
       sn=SQRT((csu2*snl)**2+(csu1*snu2-snu1*csu2*csl)**2)
       cs=snu1*snu2+csu1*csu2*csl
       sigma=ATAN2(sn,cs)
       ! alpha=azimuth at the equator
       ! sin(alpha) and cos(alpha)
       sna=csu1*csu2*snl/sn
       cs2a=1-sna**2
       cs_m=cs-(2.d0*snu1*snu2)/cs2a
       cc=flatt_surf*cs2a*(4.d0+flatt_surf*(4.d0-3.d0*cs2a))/16.d0
       ! update lambda
       lambda=ldif+(1.d0-cc)*flatt_surf*sna*(sigma+cc*sn*(cs_m+cc*cs*(-1.d0+2.d0*cs_m)))
       ! convergence control: 1.d-12 corresponds to about 0.06mm on Earth
       IF(ABS(lambda-lpre).LE.1.d-12) EXIT
    END DO
    IF(i.GE.itmax)THEN
       numerr=numerr+1
       WRITE(ierrou,*) 'dist_geoid: failed in computing distance on geoid. Spherical approximation is used instead.'
       snl=SIN(ABS(ldif))
       csl=COS(ABS(ldif))
       sn=SQRT((csu2*snl)**2+(csu1*snu2-snu1*csu2*csl)**2)
       cs=snu1*snu2+csu1*csu2*csl
       sigma=ATAN2(sn,cs)
       dist=sigma   ! angular distance in radians
       RETURN
    ENDIF
    u_squared=cs2a*ecc_squared/(1-ecc_squared)
    k1=(SQRT(1+u_squared)-1.d0)/(SQRT(1+u_squared)+1.d0)
    aa=(1.d0+k1**2/4.d0)/(1-k1)
    bb=k1*(1.d0-3.d0*k1**2/8.d0)
    delsig=bb*sn*(cs_m+bb/4.d0*(cs*(-1.d0+2.d0*cs_m**2)-&
         & bb/6.d0*cs_m*(-3.d0+4.d0*sn**2)*(-3.d0+4.d0*cs_m**2)))
    dist=(1-flatt_surf)*aa*(sigma-delsig)  ! angular distance in radians
!!$ FOR TESTING
!!$    WRITE(iun_log,*) ' === SUBROUTINE dist_geoid === '
!!$    WRITE(iun_log,*) 'Geodetic coordinates of points (long,lat) in deg: '
!!$    WRITE(iun_log,*) pt1*degrad
!!$    WRITE(iun_log,*) pt2*degrad
!!$    WRITE(iun_log,*) 'spherical approximation values: '
!!$    snl=SIN(ABS(ldif))
!!$    csl=COS(ABS(ldif))
!!$    sn=SQRT((csu2*snl)**2+(csu1*snu2-snu1*csu2*csl)**2)
!!$    cs=snu1*snu2+csu1*csu2*csl
!!$    sigma=ATAN2(sn,cs)
!!$    WRITE(iun_log,*) 'distance = ',sigma*degrad,' deg'
!!$    WRITE(iun_log,*) 'distance = ',sigma*rad_surf,' km'
!!$    WRITE(iun_log,*) ' ============================ '
  END SUBROUTINE dist_geoid

END MODULE surf_trace

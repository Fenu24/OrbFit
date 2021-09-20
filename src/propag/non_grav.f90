! ==========MODULE non_grav============
! CONTAINS
! comet_non_grav

MODULE non_grav
  USE fund_const
  USE dyn_param
  IMPLICIT NONE
  PRIVATE

  PUBLIC radpast, comet_non_grav_symm, comet_non_grav_asymm
  
CONTAINS 
! ======================================================
! Direct solar radiation pressure subroutine for 
! the purpose of the asteroid work
  SUBROUTINE radpast(x,accrad) 
! x(3) = helio-centric position vector of the asteroid\par              
! accrad(3) = geo-centric acceleration of the satellite             
!   due to the radiation pressure, for unit A/M m^2/ton\par
! \interface
!  INPUT                                                                
    REAL(KIND=dkind), INTENT (IN) :: x(3)
!  OUTPUT                                                               
    REAL(KIND=dkind), INTENT (OUT) :: accrad(3)
! \endint
! local variables                                                       
    REAL(KIND=dkind) :: rsun,rsun2
    REAL(KIND=dkind) :: au2,dsrpma,scc  
! functions                                                             
    REAL(KIND=dkind) :: prscal
! ======================================================
! distance from Sun
    rsun2=prscal(x,x) 
    rsun=dsqrt(rsun2) 
! get radiation force-factor; note the coefficient is
! drpa2m=A/M x C_R coefficient of radiation pressure (m^2/ton)
    au2=(aukm*1.d3)**2 ! AU in m
! Solar Constant at 1AU (ton/day^3) divided by velocity of light
    scc=1.37d3*(86400.d0**3)*1.d-3/vlight
! Solar flux scaled, in such a way that coefficient can be in $m^2/ton$
    dsrpma=scc/au2
! Spherical satellite assumed                                           
! $ F=\frac{\Phi}{c}\cdot\frac{C_R*A}{M}*\frac{\mathbf{x}}{|\mathbf{x}|^3}$
    accrad=dsrpma*x/(rsun2*rsun)
  END SUBROUTINE radpast
  

! ====================================================================
! COMET NON GRAV. TERMS added 29/10/2008 - Symmetric Model
! Marsden and Sekanina 1973 - "Comets and nongravitational forces. V."
! version 2.0, A. Del Vigna, Jan 2016: added derivatives
! ====================================================================
  SUBROUTINE comet_non_grav_symm(x,v,nongrav)
! interface: INPUT
    DOUBLE PRECISION, INTENT(IN), DIMENSION(3) :: x,v ! position and velocity, heliocentric
! interface: OUTPUT
    DOUBLE PRECISION, INTENT(OUT), DIMENSION(3,6) :: nongrav ! each column is one component, order R,T,W 
! end interface
!   DOUBLE PRECISION :: vvec(3),vv2, rr, vsize, prscal, factor, conv
    DOUBLE PRECISION :: rr,rr2,rv,trans_size,trans(3),xu(3),tu(3),hu(3)
    DOUBLE PRECISION :: vsize, prscal
    DOUBLE PRECISION :: al1,rr0,mesp,nesp,kesp,grfac
! ====================================================================
    rr=vsize(x)
    rr0=2.808d0         ! in AU - scale distance for water ice(frost line)
    al1=0.111262d0      ! normalizing constant
    mesp=2.15d0         ! esponent m of the formula
    nesp=5.093d0        ! esponent n of the formula
    kesp=4.6142d0       ! esponent k of the formula
    grfac=al1*((rr/rr0)**(-mesp))*(1+(rr/rr0)**nesp)**(-kesp)
    xu=x/rr         ! radial
    rr2=prscal(x,x)
    rv =prscal(x,v) 
    CALL prvec(x,v,hu)
    hu=hu/vsize(hu) ! out of plane
    trans=v-(x*rv/rr2)
    trans_size=vsize(trans)
    tu=trans/trans_size ! transversal in plane
    nongrav(1:3,1)=grfac*xu
    nongrav(1:3,2)=grfac*tu
    nongrav(1:3,3)=grfac*hu
    nongrav(1:3,4:6)=0.d0
  END SUBROUTINE comet_non_grav_symm

! ====================================================================
! COMET NON GRAV. TERMS added 07/11/2008 - Asymmetric Model
! Yeomans and Chodas 1989 - "An Asymmetric Outgassing Model for
!                            Cometary Nongravitational Accelerations"
! version 2.0, A. Del Vigna, Jan 2016: added derivatives
! ====================================================================
! Outgassing force A1*g(r')*r + A2*g(r')*t + A3*g(r')*w, where
! 1) r,t,w are the three unit vectors along the radial, transversal 
!    and out of plane directions (at time t)
! 2) r' is r computed at time t' = t - DT (DT is the delay parameter)
! 3) g(r)=alpha*(r/r0)^(-m)*[1+(r/r0)^n]^(-k)
! ====================================================================
  SUBROUTINE comet_non_grav_asymm(x,v,nongrav)
    USE ever_pitkin
! === INPUT ===========================================================================
    DOUBLE PRECISION, INTENT(IN), DIMENSION(3) :: x, v ! heliocentric position and velocity
! === OUTPUT ==========================================================================
    DOUBLE PRECISION, INTENT(OUT), DIMENSION(3,6) :: nongrav ! each column is one component, order R,T,W
                                                             ! columns 4-6 are the partials of col. 1,2,3 w.r. to
                                                             ! the delay parameter 
! =====================================================================================
    DOUBLE PRECISION               :: al1, rr0, mesp, nesp, kesp, grfac ! computation of g(r)
    DOUBLE PRECISION               :: xdelay(6)                         ! state vector at time t'=t-DT
    DOUBLE PRECISION               :: dgddt                             ! derivative of g wrt time delay dt
    DOUBLE PRECISION, DIMENSION(3) :: xu, tu, hu                        ! unit vectors at time t
    DOUBLE PRECISION               :: rrd                               ! norm of xdelay(1:3)
    DOUBLE PRECISION               :: rr2                               ! norm of x and its square
    DOUBLE PRECISION               :: t0                                ! time for fser_propag (set to 0)
! other variables
    DOUBLE PRECISION :: trans_size,trans(3)
    DOUBLE PRECISION :: vsize, prscal, rv
! =====================================================================================
    t0=0.d0
    CALL fser_propag(x,v,t0,-dyn%dp(4),gms,xdelay(1:3),xdelay(4:6))
    rrd=vsize(xdelay(1:3))
    !********************************!
    !      Computation of g(r')      !
    !********************************!
    rr0=2.808d0         ! in AU - scale distance for water ice(frost line)
    al1=0.111262d0      ! normalizing constant alpha
    mesp=2.15d0         ! esponent m of the formula
    nesp=5.093d0        ! esponent n of the formula
    kesp=4.6142d0       ! esponent k of the formula
    grfac=al1*((rrd/rr0)**(-mesp))*(1+(rrd/rr0)**nesp)**(-kesp)
    
    !*********************************!
    !   Computation of unit vectors   !
    !*********************************!
    xu=x/vsize(x)              ! radial unit vector at time t
    rr2=prscal(x,x)
    rv =prscal(x,v)
    CALL prvec(x,v,hu)
    hu=hu/vsize(hu)            ! out of plane unit vector at time t
    trans=v-(x*rv/rr2)
    trans_size=vsize(trans)
    tu=trans/trans_size        ! transversal unit vector at time t
    
    !********************************!
    !    Computation of the force    !
    !********************************!
    nongrav(1:3,1)=grfac*xu
    nongrav(1:3,2)=grfac*tu
    nongrav(1:3,3)=grfac*hu
    
    !**************************************!
    !  Derivatives of g(r) wrt time delay  !
    !**************************************!
    dgddt=grfac/rrd*(mesp+kesp*nesp/(1+(rrd/rr0)**(-nesp)))*DOT_PRODUCT(xdelay(4:6),xu)
    nongrav(1:3,4)=dgddt*xu
    nongrav(1:3,5)=dgddt*tu
    nongrav(1:3,6)=dgddt*hu
  END SUBROUTINE comet_non_grav_asymm
END MODULE non_grav

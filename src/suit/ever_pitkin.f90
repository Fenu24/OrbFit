MODULE ever_pitkin
USE output_control
USE fund_const  
IMPLICIT NONE

PRIVATE

PUBLIC fser_propag, fser_propag_der, solve_peri, solve_kepuniv, solve_kepuniv2, pdtime, s_funct

! PUBLIC s_funct, r_of_psi !maybe to be made private later

! References: Everhart, E. and Pitkin, E.T., Am. J. Phys vol 51, 712-717
!             Goodyear, W.H., Astron. J. vol. 70, 189--192
!               with Errata published on Astron. J.  446

CONTAINS
! ==========================================================
! FSER_PROPAG f,g series 2-body propagator
! ==========================================================
  SUBROUTINE fser_propag(x0,y0,t0,t,mu,x,y)
    DOUBLE PRECISION, DIMENSION(3), INTENT(IN):: x0,y0 ! initial conditions
    DOUBLE PRECISION, INTENT(IN):: t,t0,mu ! time of propag, epoch, G*mass
    DOUBLE PRECISION, DIMENSION(3), INTENT(OUT):: x,y ! propagated pos and vel
! end interface
    DOUBLE PRECISION :: s0,s1,s2,s3 !Stumpff's functions
    DOUBLE PRECISION r0, v0, vsize, prscal, alpha, sig0
    DOUBLE PRECISION psi, r,f,g,fdot,gdot 
    DOUBLE PRECISION period, dt
    LOGICAL accept, ecc_cont
    DOUBLE PRECISION ecc,q,energy
    DOUBLE PRECISION, PARAMETER :: eps=-3.d-6 ! corresponds to a period of 1000y
! compute initial range, velocity, scalar product, 2*energy 
    r0=vsize(x0)
    sig0=prscal(x0,y0)
    accept=ecc_cont(x0,y0,mu,1.d0,1.d0,ecc,q,energy)
    v0=vsize(y0)
    IF(r0.eq.0.d0.or.v0.eq.0.d0)THEN
       WRITE(*,*)' fser_propag: zero radius/velocity '
       WRITE(*,*)x0,y0
       STOP
    ENDIF
!    alpha=v0**2-2*mu/r0 ! 2*E
    alpha=2*energy
! remove integer periods
    IF(alpha.lt.eps)THEN
       period=dpig*mu/(-alpha)**(1.5d0)
!      write(*,*)' period ', period
       dt=t-t0-anint((t-t0)/period)*period
    ELSE
       dt=t-t0
    ENDIF
! solve universal Kepler's equation and compute Stumpff's functions
    CALL solve_kepuniv2(dt,r0,sig0,mu,alpha,ecc,psi,s0,s1,s2,s3)
! compute f,g series
    f=1.d0-mu*s2/r0
    g=dt-mu*s3
! compute current range, derivatives of f,g series
    r=r0*s0+sig0*s1+mu*s2
    fdot=-mu*s1/(r0*r)
    gdot=1.d0-mu*s2/r
! current position and velocity
    x=x0*f+y0*g
    y=x0*fdot+y0*gdot
    RETURN
  END SUBROUTINE fser_propag
! ==========================================================
! FSER_PROPAG_DER f,g series 2-body propagator with partials
! ==========================================================
  SUBROUTINE fser_propag_der(x0,y0,t0,t,mu,x,y,dxydxy0)
    DOUBLE PRECISION, DIMENSION(3), INTENT(IN) :: x0,y0 ! initial conditions
    DOUBLE PRECISION, INTENT(IN) :: t,t0,mu ! time of propag, epoch, G*mass
    DOUBLE PRECISION, DIMENSION(3), INTENT(OUT) :: x,y ! propagated pos and vel
    DOUBLE PRECISION, DIMENSION(6,6), INTENT(OUT) :: dxydxy0 ! state transition matrix
! end interface
    DOUBLE PRECISION :: s0,s1,s2,s3,s4,s5 !Stumpff's functions
    DOUBLE PRECISION r0, v0, vsize, prscal, alpha, sig0
    DOUBLE PRECISION psi, r,f,g,fdot,gdot
    DOUBLE PRECISION period, dt
    DOUBLE PRECISION :: u_good
    DOUBLE PRECISION, DIMENSION(3,3) :: unit3,dxdx0,dxdy0,dydx0,dydy0
    DOUBLE PRECISION, DIMENSION(3,1) :: xdd, xd
    DOUBLE PRECISION, DIMENSION(3,2) :: x_xd, tmp 
    DOUBLE PRECISION, DIMENSION(1,3) :: xdd0, xd0
    DOUBLE PRECISION, DIMENSION(2,3) :: x0_xd0
    DOUBLE PRECISION, DIMENSION(2,2) :: ff
    LOGICAL accept, ecc_cont
    DOUBLE PRECISION ecc,q,energy
    DOUBLE PRECISION, PARAMETER :: eps=-3.d-6 ! corresponds to a period of 1000y
! compute initial range, velocity, scalar product, 2*energy 
    r0=vsize(x0)
    sig0=prscal(x0,y0)
    accept=ecc_cont(x0,y0,mu,1.d0,1.d0,ecc,q,energy)
    v0=vsize(y0)
    IF(r0.eq.0.d0.or.v0.eq.0.d0)THEN
       WRITE(*,*)' fser_propag: zero radius/velocity '
       WRITE(*,*)x0,y0
       STOP
    ENDIF
!    alpha=v0**2-2*mu/r0 ! 2*E
    alpha=2*energy
! remove integer periods
    IF(alpha.lt.eps)THEN
       period=dpig*mu/(-alpha)**(1.5d0)
!      write(*,*)' period ', period
       dt=t-t0-anint((t-t0)/period)*period
    ELSE
       dt=t-t0
    ENDIF
! solve universal Kepler's equation and compute Stumpff's functions
    CALL solve_kepuniv2(dt,r0,sig0,mu,alpha,ecc,psi,s0,s1,s2,s3)
! compute f,g series
    f=1.d0-mu*s2/r0
    g=dt-mu*s3
! compute current range, derivatives of f,g series
    r=r0*s0+sig0*s1+mu*s2
    fdot=-mu*s1/(r0*r)
    gdot=1.d0-mu*s2/r
! current position and velocity
    x=x0*f+y0*g
    y=x0*fdot+y0*gdot
! additional functions for derivatives
    s4=(s2-psi**2/2)/alpha
    s5=(s3-psi**3/6)/alpha
! U function (Goodyear p. 190)
    u_good=s2*dt+mu*(psi*s4-3*s5)
    CALL eye(3,unit3)
    xd(1:3,1)=y
    xd0(1,1:3)=y0
    xdd(1:3,1)=-x*mu/r**3
    xdd0(1,1:3)=-x0*mu/r0**3
    x_xd(1:3,1)=x
    x_xd(1:3,2)=y
    x0_xd0(1,1:3)=x0
    x0_xd0(2,1:3)=y0
! two instances of fdot (errata)
    ff(1,1)= -(fdot*s1+(f-1.d0)/r0)/r0 
    ff(1,2)= -fdot*s2
    ff(2,1)= ((f-1.d0)*s1/r0)
    ff(2,2)= (f-1.d0)*s2 
    tmp= MATMUL(x_xd, ff)
    dxdx0=unit3*f + MATMUL(xd,xdd0)*u_good + MATMUL(tmp,x0_xd0)
! one instance of fdot (errata)
    ff(1,1)= -fdot*s2 
    ff(1,2)= -(gdot-1.d0)*s2 
    ff(2,1)= (f-1.d0)*s2
    ff(2,2)= g*s2 
    tmp= MATMUL(x_xd, ff)
    dxdy0=unit3*g - MATMUL(xd,xd0)*u_good +  MATMUL(tmp,x0_xd0)
! three instaces of fdot (errata)
    ff(1,1)= -fdot*(s0/(r*r0)+ 1/r**2 +1/r0**2)
    ff(1,2)= -(fdot*s1+(gdot-1.d0)/r)/r
    ff(2,1)= (fdot*s1+(f-1.d0)/r0)/r0  
    ff(2,2)= fdot*s2
    tmp= MATMUL(x_xd, ff)
! note unit3*fdot (errata)
    dydx0=unit3*fdot + MATMUL(xdd,xdd0)*u_good + MATMUL(tmp,x0_xd0)
! two instanced of fdot, not one (errata)
    ff(1,1)=-(fdot*s1+(gdot-1.d0)/r)/r
    ff(1,2)= -(gdot-1.d0)*s1/r
    ff(2,1)=fdot*s2 
    ff(2,2)=(gdot-1.d0)*s2
    tmp= MATMUL(x_xd, ff)
    dydy0=unit3*gdot - MATMUL(xdd,xd0)*u_good + MATMUL(tmp,x0_xd0)
    dxydxy0(1:3,1:3)=dxdx0
    dxydxy0(1:3,4:6)=dxdy0
    dxydxy0(4:6,1:3)=dydx0
    dxydxy0(4:6,4:6)=dydy0
    RETURN
  END SUBROUTINE fser_propag_der
! =======================================================================
! SOLVE_KEPUNIV solves universal Kepler's equation
! in inout, r0,sig0 are values at reference time t0, 
!           dt=t-t0 (possibly with integer periods removed for elliptic case)
!           mu=G*mass alpha=2*energy
!           conv_control (optional) is the control on Delta psi
! in output, psi is the universal anomaly corresponding to t
! =======================================================================
  SUBROUTINE solve_kepuniv(dt,r0,sig0,mu,alpha,psi,s0,s1,s2,s3,conv_contr) 
    DOUBLE PRECISION, INTENT(IN):: dt,r0,sig0,mu,alpha 
    DOUBLE PRECISION, INTENT(IN), OPTIONAL :: conv_contr ! convergence control
    DOUBLE PRECISION, INTENT(OUT):: psi, s0, s1, s2, s3
! end interface
    DOUBLE PRECISION contr,dpsi,psi1,fun, funp
    DOUBLE PRECISION :: a0, u0, e0, cosu, enne, du
    INTEGER j
    INTEGER, PARAMETER :: jmax=100, itx=10
! control
    IF(PRESENT(conv_contr))THEN
       contr=conv_contr
    ELSE
       contr=100*epsilon(1.d0)
    ENDIF
!initial guess: rough value
    psi=dt/r0 
    IF(alpha.lt.0.d0)THEN
       IF(abs(sig0).lt.contr)THEN
          a0=-mu/alpha
          IF(a0.gt.1.d3)THEN
             WRITE(ierrou,*)' alpha, mu, a0 ', alpha,mu, a0
             numerr=numerr+1
          ELSE
! use Kepler's equation to compute ecc. anomaly
             e0=1.d0-r0/a0
             enne=sqrt(-alpha**3)/mu
             u0=pig
             if(enne*dt.lt.0.d0)u0=-pig
             DO j=1,itx
                du=-(u0-e0*sin(u0)-enne*dt)/(1.d0-e0*cos(u0))
                u0=u0+du
                IF(abs(du).lt.contr)EXIT
             ENDDO
! if ecc. anomaly found, use relationship whith psi
             IF(j.le.itx)psi=u0/sqrt(-alpha)
          ENDIF
       ENDIF
    ELSE
       psi=psi+dt**2*sig0/(2*r0**3) ! improved initial guess
! in all other cases, stick to psi=dt/r0
    ENDIF
    DO j=1,jmax
      IF(abs(psi).gt.1.d5)THEN
         WRITE(ierrou,*)' psi, alpha, r0', psi, alpha, r0
         numerr=numerr+1
      ENDIF
      CALL s_funct(psi,alpha,s0,s1,s2,s3)
      fun=r0*s1+sig0*s2+mu*s3-dt
      funp=r0*s0+sig0*s1+mu*s2
      dpsi=-fun/funp
      psi1=psi+dpsi
      IF(psi1*psi.lt.0.d0)THEN
         psi=psi/2.d0
      ELSE
         psi=psi1
      ENDIF
      IF(abs(dpsi).lt.contr.or.abs(dpsi).lt.contr*10*abs(psi)) RETURN
    ENDDO
    IF(abs(dpsi).gt.100*contr.or.abs(dpsi).gt.contr*1000*abs(psi))THEN
       WRITE(ierrou,*)' solve_kepuniv: poor convergence ', jmax,dpsi,psi,contr
       numerr=numerr+1
    ENDIF
  END SUBROUTINE solve_kepuniv
! =======================================================================
! SOLVE_KEPUNIV2 solves universal Kepler's equation
! using elliptic and hyperbolic Kepler's equation to provide
! starting guess for universal Kepler's equation
! in inout, r0,sig0 are values at reference time t0, 
!           dt=t-t0 (possibly with integer periods removed for elliptic case)
!           mu=G*mass alpha=2*energy e0=eccentricity
!           conv_control (optional) is the control on Delta psi
! in output, psi is the universal anomaly corresponding to t
! =======================================================================
  SUBROUTINE solve_kepuniv2(dt,r0,sig0,mu,alpha,e0,psi,s0,s1,s2,s3,conv_contr,ds2da,ds0da) 
    DOUBLE PRECISION, INTENT(IN):: dt,r0,sig0,mu,alpha,e0 
    DOUBLE PRECISION, INTENT(IN), OPTIONAL :: conv_contr ! convergence control
    DOUBLE PRECISION, INTENT(OUT):: psi, s0, s1, s2, s3
    DOUBLE PRECISION, INTENT(OUT), OPTIONAL :: ds2da,ds0da
! end interface
    DOUBLE PRECISION contr,dpsi,psi1, fun, funp, pridif !for Newton's method
    DOUBLE PRECISION :: a0, u0, cosu0, u, cosu, enne, du, ell, ell0, princ! elliptic case
    DOUBLE PRECISION :: f0, f, coshf0, coshf, df, ff, psi0 ! hyperbolic case
    INTEGER j
    INTEGER, PARAMETER :: jmax=100, itx=20
! control
    IF(PRESENT(conv_contr))THEN
       contr=conv_contr
    ELSE
       contr=100*epsilon(1.d0)
    ENDIF
!initial guess: rough value
    psi=dt/r0 
    IF(alpha.lt.0.d0)THEN
! elliptic case
       a0=-mu/alpha
       enne=sqrt(-alpha**3)/mu
       IF(e0.lt.contr)THEN
          u0=0.d0
          u=enne*dt
       ELSE
! eccentic anomaly when r=r0
          cosu0=(1.d0-r0/a0)/e0
! avoiding rounding off near pericienter/apocenter
          IF(cosu0.lt.1.d0.and.cosu0.gt.-1.d0)THEN
             u0=acos(cosu0)
          ELSEIF(cosu0.ge.1.d0)THEN
             u0=0.d0
          ELSEIF(cosu0.le.-1.d0)THEN
             u0=pig
          ENDIF
          IF(sig0.lt.0.d0)u0=-u0
          u0=princ(u0)
          ell0=princ(u0-e0*sin(u0))
! eccentric anomaly from Kepler's equation        
          u=pig !starting from inflection point
          ell=princ(ell0+enne*dt)
          DO j=1,itx
             du=-(u-e0*sin(u)-ell)/(1.d0-e0*cos(u))
             u=u+du
             IF(abs(du).lt.contr*1.d3)EXIT
          ENDDO
       ENDIF
! if ecc. anomaly found, use relationship whith psi
       psi0=pridif(u,u0)/sqrt(-alpha)
    ELSEIF(alpha.gt.0.d0)THEN
! hyperbolic case
       a0=-mu/alpha
       enne=sqrt(alpha**3)/mu
! hyperbolic eccentic anomaly when r=r0
       coshf0=(1.d0-r0/a0)/e0
! avoiding rounding off near pericenter
       IF(coshf0.gt.1.d0)THEN
          f0=log(coshf0+sqrt(coshf0**2-1))
       ELSE
          f0=0.d0
       ENDIF
       IF(sig0.lt.0.d0)f0=-f0
       ell0=e0*sinh(f0)-f0
       f=0.d0 !starting from inflection point
!       f=psi*sqrt(alpha) ! starting from crude  approximation
       ell=ell0+enne*dt
! hyp. eccentric anomaly from hyperbolic Kepler's equation
       DO j=1,itx
          IF(abs(f).lt.1.5d1)THEN 
             df=-(e0*sinh(f)-f-ell)/(e0*cosh(f)-1.d0)
             ff=f+df
             IF(f*ff.lt.0.d0)THEN
                df=-f/2.0
                f=f+df
             ELSE
                f=ff
             ENDIF
          ELSE
             df=-f/2.d0
             f=f+df             
          ENDIF
          
          IF(abs(df).lt.contr*1.d3) GOTO 3
       ENDDO
       WRITE(ierrou,*)' solve_kepuniv2: non convergent hyperb. kepler eq. '
       WRITE(ierrou,*)f,df,e0,j
! if ecc. anomaly found, use relationship whith psi
3      psi0=(f-f0)/sqrt(alpha)
    ENDIF
    psi=psi0
! Newton's method in universal Kepler's equation
    DO j=1,jmax
      IF(abs(psi).gt.1.d5)THEN
         WRITE(ierrou,*)' psi, alpha, r0', psi, alpha, r0
         numerr=numerr+1
      ENDIF
      IF(PRESENT(ds2da).and.PRESENT(ds0da))THEN
         CALL s_funct(psi,alpha,s0,s1,s2,s3,DS2DA=ds2da,DS0DA=ds0da)
      ELSE
         CALL s_funct(psi,alpha,s0,s1,s2,s3)
      ENDIF
      fun=r0*s1+sig0*s2+mu*s3-dt
      funp=r0*s0+sig0*s1+mu*s2
      dpsi=-fun/funp
      IF(abs(s3).gt.1.d-2/epsilon(1.d0)) THEN
         WRITE(ierrou,*) psi, alpha
         WRITE(ierrou,*) s0,s1,s2,s3
         numerr=numerr+1
         EXIT
      ENDIF
      psi1=psi+dpsi
      IF(psi1*psi.lt.0.d0)THEN
         psi=psi/2.d0
      ELSE
         psi=psi1
      ENDIF
      IF(abs(dpsi).lt.contr.or.abs(dpsi).lt.contr*10*abs(psi)) RETURN
    ENDDO
    IF(abs(dpsi).gt.10*contr.or.abs(dpsi).gt.contr*100*abs(psi))THEN
       WRITE(ierrou,*)' solve_kepuniv2: poor convergence ', jmax,dpsi,psi,contr
       numerr=numerr+1
    ENDIF
  END SUBROUTINE solve_kepuniv2
! ==========================================================
! SOLVE_PERI 
!     finds psi (universal anomaly) and time dt at perihelion 
!       from current r0,sig0 and perihelion value peri
!             given mu=G*mass, alpha=2*Energy 
! ==========================================================
  SUBROUTINE solve_peri(r0,sig0,peri,mu,alpha,psi,dt,conv_contr)
    DOUBLE PRECISION, INTENT(IN):: r0,sig0,mu,alpha,peri
    DOUBLE PRECISION, INTENT(IN), OPTIONAL :: conv_contr ! convergence control
    DOUBLE PRECISION, INTENT(OUT):: psi, dt
    DOUBLE PRECISION :: s0, s1, s2, s3
    DOUBLE PRECISION contr,dpsi, rp, rpp, psi0, psiper, r, a0, u0, e0, cosu
    DOUBLE PRECISION coshf0, f0, enne, ell0,ell,dt0, princ
    INTEGER j, icount
    INTEGER, PARAMETER ::jmax=60
! control
    IF(PRESENT(conv_contr))THEN
       contr=conv_contr
    ELSE
       contr=100*epsilon(1.d0)
    ENDIF
    IF(alpha.le.0.d0)THEN
! use eccentric anomaly as first guess for psi
       a0=-mu/alpha
       e0=1.d0-peri/a0
       cosu=(1.d0-r0/a0)/e0
       IF(cosu.le.-1.d0)THEN
          u0=pig
       ELSEIF(cosu.ge.1.d0)THEN
          u0=0.d0
       ELSE
          u0=acos(cosu)
       ENDIF
       psi=-sign(1.d0,sig0)*u0/sqrt(-alpha)
! use classical Kepler's equation as first guess for dt
       enne=sqrt(-alpha**3)/mu
       ell0=princ(u0-e0*sin(u0))
       IF(ell0.gt.pig)ell0=ell0-dpig
       dt0=ell0/enne
    ELSE
! hyperbolic case
       a0=-mu/alpha
       e0=1.d0-peri/a0
! hyperbolic eccentic anomaly when r=r0
       coshf0=(1.d0-r0/a0)/e0
! avoiding rounding off near pericenter
       IF(coshf0.gt.1.d0)THEN
          f0=log(coshf0+sqrt(coshf0**2-1))
       ELSE
          f0=0.d0
       ENDIF
! use hyperbolic eccentric anomaly as first guess for psi
       psi=-sign(1.d0,sig0)*f0/sqrt(alpha)
! use hyperbolic Kepler's equation as first guess for dt
       enne=sqrt(alpha**3)/mu
       ell0=e0*sinh(f0)-f0
       dt0=ell0/enne
    ENDIF
    icount=1  
    dt=dt0
    DO j=1,jmax
      CALL s_funct(psi,alpha,s0,s1,s2,s3)
      IF(abs(s3).gt.1.d-2/epsilon(1.d0)) THEN
         WRITE(*,*) psi, alpha
         WRITE(*,*) s0,s1,s2,s3
         EXIT
      ENDIF
      r=r0*s0+sig0*s1+mu*s2
      rp=(r0*alpha+mu)*s1+sig0*s0
      rpp=(r0*alpha+mu)*s0+sig0*alpha*s1
      IF(ABS(rp).lt.contr)THEN
         GOTO 1
      ENDIF
      dpsi=-rp/rpp
      psi=psi+dpsi
      IF(ABS(dpsi).lt.contr.or.abs(rp).lt.contr) GOTO 1
      IF(r.lt.peri)GOTO 1      
!      IF(ABS(psi).gt.psiper.or.psi*psi0.lt.0.d0)THEN
!         WRITE(*,*)' solve_peri: unstable newton ', icount, j, psi
!         icount=icount+1
!         psi=psi0/icount
!      ENDIF
    ENDDO
    IF(verb_io.gt.9)THEN
       WRITE(*,*)' solve_peri: poor convergence ', jmax,dpsi,contr
    ENDIF
1   CONTINUE
    dt=r0*S1+sig0*s2+mu*s3
    RETURN
  END SUBROUTINE solve_peri
! =======================================================================
! R_OF_PSI  range from universal anomaly, energy, initial value, angle, mass
! =======================================================================
  DOUBLE PRECISION FUNCTION r_of_psi(psi,alpha,r0,sig0,mu)
    DOUBLE PRECISION, INTENT(IN) :: psi, alpha, r0,sig0,mu
    DOUBLE PRECISION :: s0,s1,s2,s3 ! Stumpff's functions
    CALL s_funct(psi,alpha,s0,s1,s2,s3)
    r_of_psi=r0*s0+sig0*s1+mu*s2
  END FUNCTION r_of_psi
! =======================================================================
! S_FUNCT computation of Stumpff's functions s0, s1 s2, s3
! =======================================================================
  SUBROUTINE s_funct(psi,alpha,s0,s1,s2,s3,conv_contr,ds2da,ds0da)
    DOUBLE PRECISION, INTENT(IN) :: psi, alpha ! univ. ecc. anom., 2*energ 
    DOUBLE PRECISION, INTENT(OUT) :: s0,s1,s2,s3 ! Stumpff functions
    DOUBLE PRECISION, INTENT(IN), OPTIONAL :: conv_contr ! covergence control
    DOUBLE PRECISION, INTENT(OUT), OPTIONAL :: ds0da,ds2da
    DOUBLE PRECISION beta, term2, term3, term0, term1, psi2, s02, s12
    DOUBLE PRECISION contr, dif1, dif2, sval, overfl, termda
    INTEGER j,jj,nhalf
    INTEGER, PARAMETER:: jmax=70, halfmax=30
    DOUBLE PRECISION, PARAMETER :: betacontr=1.d2 !to make it not operational
! control
    IF(PRESENT(conv_contr))THEN
       contr=conv_contr
    ELSE
       contr=100*epsilon(1.d0)
    ENDIF
!    overfl=HUGE(1.d0)/1.d3
    overfl=1.d0/epsilon(1.d0)
! compute s2,s3
    beta=alpha*psi**2 
! computation of d(s2)/d(alpha)
    IF(PRESENT(ds2da))THEN
       termda=psi**4/24.d0
       ds2da=termda
       DO j=1,jmax
         termda=termda*beta/((2.d0*j+4.d0)*(2.d0*j+3.d0))
         ds2da=ds2da+(j+1)*termda
         IF(abs(termda*(j+1)).lt.contr)THEN
!            write(*,*)'j,term_j=',j,termda*(j+1)
            EXIT
         ENDIF
         IF(abs(termda*(j+1)).gt.overfl) EXIT
       ENDDO
       IF(j.gt.jmax.or.abs(termda*(j+1)).gt.overfl)THEN
          WRITE(*,*)' ds2da, termda', ds2da, termda
       ENDIF
    ENDIF
! if beta>1000 then...
    IF(abs(beta).lt.betacontr)THEN
       term2=psi**2/2.d0
       term3=term2*psi/3.d0
       s2=term2
       s3=term3
       DO j=1,jmax
          term2=term2*beta/((2.d0*j+1.d0)*(2.d0*j+2.d0))
          s2=s2+term2
          IF(abs(term2).gt.overfl)EXIT
          IF(abs(term2).lt.contr)EXIT
       ENDDO
       IF(j.gt.jmax.or.abs(term2).gt.overfl)THEN
          WRITE(*,*)' psi,alpha, j, term2 ', psi, alpha, j-1, term2
       ENDIF
       DO j=1,jmax
          term3=term3*beta/((2*j+2)*(2*j+3))
          s3=s3+term3
          IF(abs(term3).gt.overfl)EXIT
          IF(abs(term3).lt.contr)EXIT
       ENDDO
       IF(j.gt.jmax.or.abs(term3).gt.overfl)THEN
          WRITE(*,*)' psi,alpha,j, term3', psi,alpha,j-1, term3
       ENDIF
! recursive formulae
       s1=psi+alpha*s3
       s0=1.d0+alpha*s2
       IF(PRESENT(ds0da))THEN
          ds0da=s2+alpha*ds2da
       ENDIF
    ELSE
       psi2=psi
       nhalf=0
       DO jj=1,halfmax
         psi2=psi2*0.5d0
         nhalf=nhalf+1
         beta=alpha*psi2**2
         IF(abs(beta).lt.betacontr)EXIT
       ENDDO
!       WRITE(*,*)' halving', psi,nhalf,beta
       term0=1.d0
       term1=psi2
       s0=1.d0
       s1=psi2
       DO j=1,jmax
          term0=term0*beta/((2*j-1)*(2*j))
          s0=s0+term0
          IF(abs(term0).gt.overfl)EXIT
          IF(abs(term0).lt.contr)EXIT
       ENDDO
       IF(j.gt.jmax.or.abs(term0).gt.overfl)THEN
           WRITE(*,*)' psi,alpha, j, term0 ', psi,alpha,j-1, term0
       ENDIF
       DO j=1,jmax
          term1=term1*beta/((2*j)*(2*j+1))
          s1=s1+term1          
          IF(abs(term1).gt.overfl)EXIT
          IF(abs(term1).lt.contr)EXIT
       ENDDO
       IF(j.gt.jmax.or.abs(term1).gt.overfl)THEN
          WRITE(*,*)'psi,alpha, j, term1 ', psi,alpha,j-1, term1
       ENDIF
! duplication formula
       DO jj=1,nhalf
         s02=2*s0**2-1.d0
         s12=2*s0*s1
         s0=s02
         s1=s12
       ENDDO
       IF(abs(s0).gt.overfl.or.abs(s1).gt.overfl)THEN
          WRITE(ierrou,*)' s_funct: halving failed, psi, s0,s1 ', psi,s0,s1 
          numerr=numerr+1
       ENDIF
! recursive formulae
       s3=(s1-psi)/alpha
       s2=(s0-1.d0)/alpha
!       WRITE(*,*)s0,s1,s2,s3
    ENDIF
! controls
    dif1=2*s2+alpha*s2**2-s1**2
    dif2=s2+s0*s2-s1**2
    sval=s1**2+s2**2
    IF(abs(dif1).gt.100*contr*sval.or.abs(dif2).gt.100*contr*sval)THEN
    IF(abs(dif1).gt.1000*contr.or.abs(dif2).gt.1000*contr)THEN
       WRITE(ierrou,100)psi,alpha,s0,s1,s2,s3
  100  FORMAT('psi, alpha, s functions ', 1p, 6d13.5)
       WRITE(ierrou,101)dif1,dif2,contr*sval
  101  FORMAT('check id1, check id2, contr ', 1p, 3d13.5)
       numerr=numerr+1
    ENDIF
    ENDIF
  END SUBROUTINE s_funct
! ==========================================================
! PDTIME 
!     finds psi (universal anomaly) at current time and time dt=t-t0 
!     where t0 is the time of passage at perihelion
!     from the perihelion distance r0 and the current value of r.
!     mu=G*mass, alpha=2*Energy 
! ==========================================================
  SUBROUTINE pdtime(r0,r,mu,alpha,psi_ini,psi,dt,conv_contr)
    DOUBLE PRECISION, INTENT(IN):: r0,r,mu,alpha
    DOUBLE PRECISION, INTENT(IN) :: psi_ini !starting guess
    DOUBLE PRECISION, INTENT(IN), OPTIONAL :: conv_contr ! covergence control
    DOUBLE PRECISION, INTENT(OUT):: psi,dt
    DOUBLE PRECISION :: s0,s1,s2,s3
    DOUBLE PRECISION contr,dpsi,fpsi,dfpsi
    INTEGER j, icount
    INTEGER, PARAMETER ::jmax=60
! control
    IF(PRESENT(conv_contr))THEN
       contr=conv_contr
    ELSE
       contr=100*epsilon(1.d0)
    ENDIF
! starting guess for Newton's method
    psi=psi_ini
    icount=1  
    DO j=1,jmax
      CALL s_funct(psi,alpha,s0,s1,s2,s3)
      fpsi=r0*s0+mu*s2-r
      dfpsi=(r0*alpha+mu)*s1
      IF(ABS(fpsi).lt.contr)THEN
!         WRITE(*,*)j,psi,fpsi,dfpsi,' no correction'
         GOTO 1
      ENDIF
      dpsi=-fpsi/dfpsi
      WRITE(*,*)j,psi,fpsi,dfpsi
      psi=psi+dpsi
      IF(ABS(dpsi).lt.contr.or.abs(fpsi).lt.contr) GOTO 1

      IF(r.lt.r0) THEN
         WRITE(*,*)'pdtime: error! r,r0',r,r0
         STOP
      ENDIF

    ENDDO
    IF(verb_io.gt.9)THEN
       WRITE(*,*)' pdtime: poor convergence ', jmax,dpsi,contr
    ENDIF
1   CONTINUE
    dt=r0*s1+mu*s3
!    WRITE(*,*)' dt ',dt
    RETURN
  END SUBROUTINE pdtime
END MODULE ever_pitkin
! =======================================================================
! ECC_CONT  to estimate eccentricity etc. and avoid going too far
! is true if solution is "acceptable", that is ecc.le.eccmax.and.q.le.qmax
! =======================================================================
LOGICAL FUNCTION ecc_cont(x,y,mu,eccmax,qmax,ecc,q,energy)
   IMPLICIT NONE
   DOUBLE PRECISION, INTENT(IN) :: x(3),y(3) ! position, velocity
   DOUBLE PRECISION, INTENT(IN) :: mu ! G*mass
   DOUBLE PRECISION, INTENT(IN) :: eccmax,qmax ! controls
   DOUBLE PRECISION, INTENT(OUT) :: ecc,q,energy ! eccentricity, perihelion, energy
   DOUBLE PRECISION vel2, prscal, vsize, r0, ang(3), vlenz(3), gei,gei2
! compute eccentricity
!  radius and velocity squared                                          
   vel2=prscal(y,y)          ! velocity
   r0=vsize(x)               ! distance from center
   call prvec(x,y,ang)       !  angular momentum 
! non singular first element
   gei2=prscal(ang,ang)
   gei=sqrt(gei2)
   IF(gei.eq.0.d0)THEN
      WRITE(*,*) ' ecc_cont: zero angular momentum ',x,y
      STOP
   ENDIF
   call prvec(y,ang,vlenz) ! Lenz vector
   vlenz=vlenz*(1.d0/mu)-x(1:3)*(1.d0/r0)
   ecc=vsize(vlenz)
   q=gei2/(mu*(1.d0+ecc))
   energy=vel2/2.d0-mu/r0
   IF(ecc.le.eccmax.and.q.le.qmax)THEN
      ecc_cont=.true. 
   ELSE
      ecc_cont=.false.
   ENDIF
 END FUNCTION ecc_cont

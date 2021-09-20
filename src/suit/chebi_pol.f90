MODULE chebi_pol

PUBLIC smoocn, pcwlgi
! common data


! Max degree of polynomial interpolation using Legendre, Chebychev
! or similar polynomials
INTEGER, PARAMETER :: lgintx=30
INTEGER, PARAMETER :: lgin1x=lgintx+1

PUBLIC lgintx,lginlx

CONTAINS
! Copyright (C) 1999 by Mario Carpino (carpino@brera.mi.astro.it)       
! Version: November 11, 1999                                            
! --------------------------------------------------------------------- 
!                                                                       
!  ***************************************************************      
!  *                                                             *      
!  *                         C H E B Y M                         *      
!  *                                                             *      
!  *       Recursive computation of Chebyschev polynomials       *      
!  *                                                             *      
!  ***************************************************************      
!                                                                       
! INPUT:    X         -  Independent variable                           
!           N         -  Degree                                         
!                                                                       
! OUTPUT:   P(i)      -  Chebyschev polynomials (i=0,N)                 
!                                                                       
      SUBROUTINE chebym(x,n,p) 
      IMPLICIT NONE 
                                                                        
      INTEGER n 
      DOUBLE PRECISION x,p(0:n) 
                                                                        
      INTEGER j 
                                                                        
      p(0)=1.d0 
      IF(n.EQ.0) RETURN 
      p(1)=x 
      IF(n.EQ.1) RETURN 
      DO 1 j=2,n 
      p(j)=2.d0*x*p(j-1)-p(j-2) 
    1 END DO 
                                                                        
      END SUBROUTINE chebym                                          
! Copyright (C) 1999 by Mario Carpino (carpino@brera.mi.astro.it)       
! Version: November 11, 1999                                            
! --------------------------------------------------------------------- 
!                                                                       
!  *****************************************************************    
!  *                                                               *    
!  *                         C H E B Y D                           *    
!  *                                                               *    
!  *         Recursive computation of Chebychev polynomials        *    
!  *             with derivatives up to second order               *    
!  *                                                               *    
!  *****************************************************************    
!                                                                       
! INPUT:    X         -  Independent variable                           
!           N         -  Degree                                         
!           NTDER     -  Required time derivative (=0, 1 or 2)          
!                                                                       
! OUTPUT:   P(i)      -  (i=0,N) Chebychev polynomials of degree i P(x) 
!           PD(i)     -  First derivative dP/dx (if NTDER >= 1)         
!           PDD(i)    -  Second derivative d2P/dx2 (if NTDER >= 2)      
!                                                                       
      SUBROUTINE chebyd(x,n,ntder,p,pd,pdd) 
      IMPLICIT NONE 
                                                                        
      INTEGER n,ntder 
      DOUBLE PRECISION x,p(0:n),pd(0:n),pdd(0:n) 
                                                                        
      INTEGER j 
      DOUBLE PRECISION x2 
                                                                        
      IF(ntder.LT.0 .OR. ntder.GT.2) STOP '**** chebyd: ntder = ? ****' 
                                                                        
      x2=2.d0*x 
                                                                        
! Chebychev polynomials                                                 
      p(0)=1.d0 
      IF(n.EQ.0) GOTO 2 
      p(1)=x 
      IF(n.EQ.1) GOTO 2 
      DO 1 j=2,n 
      p(j)=x2*p(j-1)-p(j-2) 
    1 END DO 
                                                                        
! First derivatives                                                     
    2 CONTINUE 
      IF(ntder.LT.1) RETURN 
      pd(0)=0.d0 
      IF(n.EQ.0) GOTO 4 
      pd(1)=1.d0 
      IF(n.EQ.1) GOTO 4 
      DO 3 j=2,n 
      pd(j)=x2*pd(j-1)+2*p(j-1)-pd(j-2) 
    3 END DO 
                                                                        
! Second derivatives                                                    
    4 CONTINUE 
      IF(ntder.LT.2) RETURN 
      pdd(0)=0.d0 
      IF(n.EQ.0) RETURN 
      pdd(1)=0.d0 
      IF(n.EQ.1) RETURN 
      DO 5 j=2,n 
      pdd(j)=x2*pdd(j-1)+4*pd(j-1)-pdd(j-2) 
    5 END DO 
                                                                        
      END SUBROUTINE chebyd                                          
! Copyright (C) 1999 by Mario Carpino (carpino@brera.mi.astro.it)       
! Version: November 11, 1999                                            
! --------------------------------------------------------------------- 
!                                                                       
!  *****************************************************************    
!  *                                                               *    
!  *                         L G N I N T                           *    
!  *                                                               *    
!  *    Least squares interpolation using Chebychev polynomials    *    
!  *                                                               *    
!  *****************************************************************    
!                                                                       
! INPUT:    X(i)      -  (i=1,NPT) Dependent variable                   
!           T(i)      -  (i=1,NPT) Independent variable (NORMALIZED)    
!           NPT       -  Number of data points                          
!           NGR       -  Polynomial degree                              
!                                                                       
! OUTPUT:   COEF(i)   -  (i=0,NGR) Coefficients of Legendre polynomials 
!           SIGMA     -  Fit RMS error                                  
!                                                                       
      SUBROUTINE lgnint(x,t,npt,ngr,coef,sigma) 
      IMPLICIT NONE 
                                                                        
      INTEGER npt,ngr 
      DOUBLE PRECISION x(npt),t(npt),coef(0:ngr),sigma 
                                                                        
      INTEGER i,j,k,indp 
      DOUBLE PRECISION atx(0:lgintx),covc(0:lgintx,0:lgintx) 
      DOUBLE PRECISION poly(0:lgintx),v(0:lgintx),xk,xt,resk,cx                        
      DOUBLE PRECISION eps ! round off control
                                                                        
      IF(ngr.GT.lgintx) STOP '**** lgnint: ngr > lgintx ****' 
                                                                        
      DO 1 i=0,ngr 
      atx(i)=0.d0 
      DO 1 j=i,ngr 
      covc(i,j)=0.d0 
    1 CONTINUE 
                                                                        
! Computation of the scalar products with Legendre polynomials          
      DO 2 k=1,npt 
      xk=x(k) 
      CALL chebym(t(k),ngr,poly) 
      DO 2 i=0,ngr 
      atx(i)=atx(i)+poly(i)*xk 
      DO 2 j=i,ngr 
      covc(i,j)=covc(i,j)+poly(i)*poly(j) 
    2 CONTINUE 
                                                                        
! Least squares fit                                                     
      eps=epsilon(1.d0)*100 
      CALL tchol(covc,lgintx+1,ngr+1,indp,eps) 
      IF(indp.NE.0) THEN 
          WRITE(6,100) indp 
          STOP '**** lgnint: abnormal end ****' 
      END IF 
  100 FORMAT(' **** lgnint: INDP =',i4,' ****') 
      CALL inver(covc,v,lgintx+1,ngr+1) 
                                                                        
! Polynomial coefficients                                               
      DO 4 i=0,ngr 
      cx=0.d0 
      DO 3 j=0,ngr 
      cx=cx+covc(i,j)*atx(j) 
    3 END DO 
      coef(i)=cx 
    4 END DO 
                                                                        
! RMS error                                                             
      sigma=0.d0 
      DO 6 k=1,npt 
      CALL chebym(t(k),ngr,poly) 
      xt=0.d0 
      DO 5 j=0,ngr 
      xt=xt+coef(j)*poly(j) 
    5 END DO 
      resk=x(k)-xt 
      sigma=sigma+resk**2 
    6 END DO 
      sigma=sqrt(sigma/npt) 
                                                                        
      END SUBROUTINE lgnint                                          
! Copyright (C) 1999 by Mario Carpino (carpino@brera.mi.astro.it)       
! Version: November 11, 1999                                            
! --------------------------------------------------------------------- 
!                                                                       
!  *****************************************************************    
!  *                                                               *    
!  *                         P C W L G I                           *    
!  *                                                               *    
!  *      Piece-wise interpolation using Chebychev polynomials     *    
!  *                                                               *    
!  *****************************************************************    
!                                                                       
! ARGUMENTS:                                                            
!    X(i,k)  (IN)   -  (i=1,N; k=1,ND) Dependent variable               
!    T(i)    (IN)   -  (i=1,N) Independent variable                     
!    N       (IN)   -  Number of data points                            
!    ND      (IN)   -  Actual dimension of dependent variable           
!    NX      (IN)   -  First dimension of X and T arrays                
!                                                                       
!    TI      (IN)   -  Interpolation time                               
!    XI(k)   (OUT)  -  Interpolated data: value                         
!    XID(k)  (OUT)  -  Interpolated data: first derivative              
!    XIDD(k) (OUT)  -  Interpolated data: second derivative             
!                                                                       
!    NL      (IN)   -  Length of interpolation interval                 
!    NV      (IN)   -  Length of unused part of interpolation           
!    NS      (IN)   -  Length of superposition                          
!    NSMORD  (IN)   -  Order of smoothness for SMOOCN                   
!    NTDER   (IN)   -  Requested order of time derivatives              
!    START   (IN)   -  Startup flag (must be a variable)                
!    RMSX(k) (IN)   -  Max allowed interpolation RMS                    
!                                                                       
!    COEF(0:NGX,NIX,ND)  (WORK)  -  Polynomial coefficients             
!    TLIM(2,NIX)         (WORK)  -  Limits of interpolation intervals   
!    VOID(NIX)           (WORK)  -  (logical) interpolation flag        
!    NGX     (IN)   -  Max polynomial degree (dim of COEF)              
!    NIX     (IN)   -  Max number of intervals (dim of COEF,VOID)       
!    IRET    (OUT)  -  Return code:                                     
!                        0   = OK                                       
!                        1-4 = time is out of limits                    
!                                                                       
      SUBROUTINE pcwlgi(x,t,n,nd,nx,                                    &
     &                  ti,xi,xid,xidd,                                 &
     &                  nl,nv,ns,nsmord,ntder,start,rmsx,               &
     &                  coef,tlim,void,ngx,nix,iret)                    
      IMPLICIT NONE 
                                                                        
      INTEGER n,nd,nx,nl,nv,ns,nsmord,ntder,ngx,nix,iret 
      DOUBLE PRECISION x(nx,nd),t(nx) 
      DOUBLE PRECISION ti,xi(nd),xid(nd),xidd(nd),rmsx(nd) 
      DOUBLE PRECISION coef(0:ngx,nix,nd),tlim(2,nix) 
      LOGICAL start,void(nix) 
                                                                        
      INTEGER ndx 
      PARAMETER (ndx=10) 
                                                                        
      INTEGER nu,na,ni,i,k,kint,kkn,kip,ksip,kofp,kis,ksis,kofs 
      DOUBLE PRECISION dtmed,t1,t2,tin,ft,ft2,sigma,ss,ssd,ssdd 
      DOUBLE PRECISION ts1,ts2,fc,c1,c1d,c1dd,c2,c2d,c2dd 
      DOUBLE PRECISION tn(lgin1x) 
      DOUBLE PRECISION pol(0:lgintx),pold(0:lgintx),poldd(0:lgintx) 
      DOUBLE PRECISION xi1(ndx),xi1d(ndx),xi1dd(ndx) 
      DOUBLE PRECISION xi2(ndx),xi2d(ndx),xi2dd(ndx) 
      LOGICAL super 
                                                                        
      IF(ntder.LT.0 .OR. ntder.GT.2) STOP '**** pcwlgi: ntder = ? ****' 
                                                                        
      nu=nl-2*(nv+ns) 
      na=ns+nu 
      ni=(n-nl)/na 
      dtmed=(t(n)-t(1))/(n-1) 
                                                                        
      IF(start) THEN 
          IF(nu.LT.0) STOP '**** pcwlgi: nu < 0  ****' 
          IF(ni.GT.nix) STOP '**** pcwlgi: ni > nix  ****' 
          DO 12 i=1,ni 
          void(i)=.true. 
   12     CONTINUE 
          IF(nl.GT.ngx) STOP '**** pcwlgi: nl > ngx  ****' 
          IF(nl.GT.lgintx) STOP '**** pcwlgi: nl > lgintx  ****' 
          IF(nd.GT.ndx) STOP '**** pcwlgi: nd > ndx  ****' 
          start=.false. 
      END IF 
                                                                        
      IF(ti.LT.t(1).OR.ti.GT.t(n)) THEN 
          iret=1 
          RETURN 
      END IF 
                                                                        
! KINT = number of the interval t(kint) <= ti < t(kint+1)               
      kint=(ti-t(1))/dtmed+1 
    1 CONTINUE 
      IF(ti.LT.t(kint)) THEN 
          kint=kint-1 
      ELSEIF(ti.GE.t(kint+1)) THEN 
          kint=kint+1 
      ELSE 
          GOTO 2 
      END IF 
      IF(kint.LT.1.OR.kint.GT.n-1) THEN 
          iret=2 
          RETURN 
      END IF 
      GOTO 1 
    2 CONTINUE 
                                                                        
! KIP = number of "primary" interpolation                               
      kkn=kint-nv-1 
      kip=kkn/(ns+nu)+1 
      IF(kkn.LT.0) kip=kip-1 
      IF(kip.LT.1.OR.kip.GT.ni) THEN 
          iret=3 
          RETURN 
      END IF 
                                                                        
! KSIP = starting interval of primary interpolation                     
      ksip=(ns+nu)*(kip-1)+1 
      IF(ksip.LT.1.OR.ksip.GT.n-nl)                                     &
     &    STOP '**** pcwlgi: internal error (01) ****'                  
                                                                        
! KOFP = offset with respect to starting interval (interval number)     
      kofp=kint-ksip+1 
      IF(kofp.LE.nv) STOP '**** pcwlgi: internal error (02) ****' 
      IF(kofp.GT.nv+ns+nu) STOP '**** pcwlgi: internal error (03) ****' 
                                                                        
! SUPER = superposition (with previous interpolation)                   
      super=(kofp.LE.nv+ns) 
                                                                        
! Primary interpolation                                                 
      IF(void(kip)) THEN 
! Normalization of time axis                                            
          t1=t(ksip) 
          t2=t(ksip+nl) 
          tlim(1,kip)=t1 
          tlim(2,kip)=t2 
          ft=2/(t2-t1) 
          DO 3 i=1,nl+1 
          tn(i)=ft*(t(i+ksip-1)-t1)-1 
    3     CONTINUE 
          DO 4 i=1,nd 
          CALL lgnint(x(ksip,i),tn,nl+1,nl,coef(0,kip,i),sigma) 
          IF(sigma.GT.rmsx(i)) STOP '**** pcwlgi: rms > rmsx ****' 
    4     CONTINUE 
          void(kip)=.false. 
      END IF 
      t1=tlim(1,kip) 
      t2=tlim(2,kip) 
      ft=2/(t2-t1) 
      ft2=ft**2 
      tin=ft*(ti-t1)-1 
      CALL chebyd(tin,nl,ntder,pol,pold,poldd) 
      DO 6 i=1,nd 
         ss  =0.d0 
         ssd =0.d0 
         ssdd=0.d0 
         DO 5 k=0,nl 
            ss  =ss  +coef(k,kip,i)*pol(k) 
            IF(ntder.ge.1)ssd =ssd +coef(k,kip,i)*pold(k) 
            IF(ntder.ge.2)ssdd=ssdd+coef(k,kip,i)*poldd(k) 
    5    CONTINUE 
         xi1(i)  =ss 
         xi1d(i) =ssd*ft 
         xi1dd(i)=ssdd*ft2 
    6 END DO 
                                                                        
! Auxiliary interpolation (which is before the primary)                 
      IF(super) THEN 
          kis=kip-1 
! KSIS = starting interval of auxiliary interpolation                   
          ksis=(ns+nu)*(kis-1)+1 
          IF(kis.LE.0) THEN 
              iret=4 
              RETURN 
          END IF 
          IF(ksis.LT.1 .OR. ksis.GT.n-nl)                               &
     &           STOP '**** pcwlgi: internal error (04) ****'           
! KOFS = offset with respect to starting interval (interval number)     
          kofs=kint-ksis+1 
          IF(kofs.LE.nv+ns+nu)                                          &
     &        STOP '**** pcwlgi: internal error (05) ****'              
          IF(kofs.GT.nv+2*ns+nu)                                        &
     &        STOP '**** pcwlgi: internal error (06) ****'              
          IF(void(kis)) THEN 
              t1=t(ksis) 
              t2=t(ksis+nl) 
              tlim(1,kis)=t1 
              tlim(2,kis)=t2 
              ft=2/(t2-t1) 
              DO 13 i=1,nl+1 
              tn(i)=ft*(t(i+ksis-1)-t1)-1 
   13         CONTINUE 
              DO 14 i=1,nd 
              CALL lgnint(x(ksis,i),tn,nl+1,nl,coef(0,kis,i),sigma) 
   14         CONTINUE 
              void(kis)=.false. 
          END IF 
          t1=tlim(1,kis) 
          t2=tlim(2,kis) 
          ft=2/(t2-t1) 
          ft2=ft**2 
          tin=ft*(ti-t1)-1 
          CALL chebyd(tin,nl,ntder,pol,pold,poldd) 
          DO 16 i=1,nd 
             ss  =0.d0 
             ssd =0.d0 
             ssdd=0.d0 
             DO 15 k=0,nl 
                ss  =ss  +coef(k,kis,i)*pol(k) 
                IF(ntder.ge.1)ssd =ssd +coef(k,kis,i)*pold(k) 
                IF(ntder.ge.2)ssdd=ssdd+coef(k,kis,i)*poldd(k) 
   15        CONTINUE 
             xi2(i)  =ss 
             xi2d(i) =ssd*ft 
             xi2dd(i)=ssdd*ft2 
   16     CONTINUE 
! Superposition of the two interpolations                               
          ts1=t(ksip+nv) 
          ts2=t(ksip+nv+ns) 
          ft=1.d0/(ts2-ts1) 
          ft2=ft**2 
          fc=ft*(ti-ts1) 
          IF(fc.LT.0.d0.OR.fc.GT.1.d0)                                  &
     &        STOP '**** pcwlgi: internal error (07) ****'              
          CALL smoocn(fc,c1,c1d,c1dd,nsmord) 
          c1d =c1d*ft 
          c1dd=c1dd*ft2 
          c2  = 1-c1 
          c2d =  -c1d 
          c2dd=  -c1dd 
          DO 7 i=1,nd 
          xi(i)=c1*xi1(i)+c2*xi2(i) 
          IF(ntder.GE.1) THEN 
              xid(i)=c1d*xi1(i)+c1*xi1d(i)+c2d*xi2(i)+c2*xi2d(i) 
          ELSE 
              xid(i)=0.d0 
          END IF 
          IF(ntder.GE.2) THEN 
              xidd(i)=c1dd*xi1(i)+2*c1d*xi1d(i)+c1*xi1dd(i)             &
     &               +c2dd*xi2(i)+2*c2d*xi2d(i)+c2*xi2dd(i)             
          ELSE 
              xidd(i)=0.d0 
          END IF 
    7     CONTINUE 
      ELSE 
          DO 8 i=1,nd 
          xi(i)=xi1(i) 
          IF(ntder.GE.1) THEN 
              xid(i)=xi1d(i) 
          ELSE 
              xid(i)=0.d0 
          END IF 
          IF(ntder.GE.2) THEN 
              xidd(i)=xi1dd(i) 
          ELSE 
              xidd(i)=0.d0 
          END IF 
    8     CONTINUE 
      END IF 
                                                                        
      iret=0 
                                                                        
      END SUBROUTINE pcwlgi                                          
! Copyright (C) 1999 by Mario Carpino (carpino@brera.mi.astro.it)       
! Version: November 11, 1999                                            
! --------------------------------------------------------------------- 
!                                                                       
!  *****************************************************************    
!  *                                                               *    
!  *                         S M O O C N                           *    
!  *                                                               *    
!  *               Smooth Cn link between 0 and 1                  *    
!  *                                                               *    
!  *****************************************************************    
!                                                                       
! INPUT:    X         -  Independent variable                           
!           NSMORD      -  Order of smoothness                          
!                                                                       
! OUTPUT:   Y         -  y(x)                                           
!           Y1        -  dy/dx                                          
!           Y2        -  d2y/dx2                                        
!                                                                       
! The function y(x) satisfies the following conditions:                 
!   1) y(x)=0 for x <= 0                                                
!   2) y(x)=1 for x >= 1                                                
!   3) y(x) is continuous with its derivatives up to order NSMORD       
!           for any x                                                   
!                                                                       
! In the interval 0 <= x <= 1, y(x) is defined as the minimum           
! degree polynomial P(x) satisfying the conditions:                     
!   1) P(0) = 0                                                         
!   2) P(1) = 1                                                         
!   3) d^n P/dx^n = 0 for x=0 and x=1 (n=1,NSMORD)                      
!                                                                       
! Since these are 2*NSMORD+2 independent equations, they can be         
! satisfied by a polynomial of degree 2*NSMORD+1                        
!                                                                       
      SUBROUTINE smoocn(x,y,y1,y2,nsmord) 
      IMPLICIT NONE 
                                                                        
      INTEGER nsmord 
      DOUBLE PRECISION x,y,y1,y2 
! Max order of smoothness for SMOOCN routine
      INTEGER nd1x,nd2x
      PARAMETER (nd1x=16)
      PARAMETER (nd2x=nd1x+1)
                                                                        
      INTEGER i,k,ng,nd1,nd2,ider,ising 
      INTEGER ic(nd1x),id(nd1x) 
      DOUBLE PRECISION a(nd1x,nd2x),det,xn,xn2 
      DOUBLE PRECISION cp(nd1x,nd1x),cp1(nd1x,nd1x),cp2(nd1x,nd1x) 
      LOGICAL first,cpcmp(nd1x) 
      SAVE first,cpcmp,cp,cp1,cp2 
      DATA first/.true./ 
                                                                        
!     SAVE first,cpcmp,cp,cp1,cp2                                       
                                                                        
      IF(nsmord.LT.0) STOP '**** smoocn: internal error (01) ****' 
                                                                        
      IF(first) THEN 
          DO 10 k=1,nd1x 
          cpcmp(k)=.true. 
   10     CONTINUE 
          first=.false. 
      END IF 
                                                                        
! Polynomial degree                                                     
      ng=2*nsmord+1 
      nd1=nsmord+1 
      nd2=nsmord+2 
      IF(nd1.GT.nd1x) STOP '**** smoocn: nd1 > nd1x ****' 
                                                                        
      IF(cpcmp(nd1)) THEN 
                                                                        
          DO 1 k=1,nd1 
          i=2*k-1 
          ic(k)=1 
          id(k)=i 
          a(1,k)=ic(k) 
          a(1,nd2)=0.5d0 
    1     CONTINUE 
                                                                        
          DO 3 ider=1,nsmord 
          DO 2 k=1,nd1 
          i=2*k-1 
          IF(id(k).GT.0) THEN 
              ic(k)=ic(k)*id(k) 
          ELSE 
              ic(k)=0 
          END IF 
          id(k)=id(k)-1 
          a(ider+1,k)=ic(k) 
    2     CONTINUE 
          a(ider+1,nd2)=0.d0 
    3     CONTINUE 
                                                                        
          CALL matin(a,det,nd1,1,nd1x,ising,0) 
          IF(ising.NE.0) STOP '**** smoocn: internal error (02) ****' 
                                                                        
          DO 4 k=1,nd1 
          i=2*k-1 
          ic(k)=1 
          id(k)=i 
    4     cp(nd1,k)=a(k,nd2) 
                                                                        
          DO 5 k=1,nd1 
          i=2*k-1 
          IF(id(k).GT.0) THEN 
              ic(k)=ic(k)*id(k) 
          ELSE 
              ic(k)=0 
          END IF 
          id(k)=id(k)-1 
    5     cp1(nd1,k)=cp(nd1,k)*ic(k) 
                                                                        
          DO 6 k=1,nd1 
          i=2*k-1 
          IF(id(k).GT.0) THEN 
              ic(k)=ic(k)*id(k) 
          ELSE 
              ic(k)=0 
          END IF 
    6     cp2(nd1,k)=cp(nd1,k)*ic(k) 
                                                                        
          cpcmp(nd1)=.false. 
      END IF 
                                                                        
      xn=2*x-1 
      y=0.d0 
      y1=0.d0 
      y2=0.d0 
                                                                        
      IF(x.LE.0.d0) RETURN 
      IF(x.GE.1.d0) THEN 
          y=1.d0 
          RETURN 
      END IF 
                                                                        
      xn2=xn**2 
      DO 11 k=nd1,1,-1 
      y=y*xn2+cp(nd1,k) 
      y1=y1*xn2+cp1(nd1,k) 
   11 END DO 
      y=0.5d0+y*xn 
      y1=2*y1 
      DO 12 k=nd1,2,-1 
      y2=y2*xn2+cp2(nd1,k) 
   12 END DO 
      y2=4*y2*xn 
                                                                        
      END SUBROUTINE smoocn                                          

END MODULE chebi_pol

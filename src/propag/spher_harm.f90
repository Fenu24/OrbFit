! SPHER_HARM: Gravitational potential routines for satellites
MODULE spher_harm
USE fund_const
USE output_control
USE planet_masses
IMPLICIT NONE 
PRIVATE
! ROUTINES 
!    geopot : Input of geopotential coefficients                   
!    namcoe : Computes the l, m and i of an harmonic with a given ndex 
!    parm10 : Computes accelerations due to spherical harmonics    
!             (Not normalized)                                      
!    legit  : Recursive computation of Legendre polynomials        
!             and Legendre associated functions                   
!                                                                       
!  FUNCTIONS                                                       
!    indcoe : Address of harmonics (l,m,i)                        
!    numcoe : Number of coefficients up to degree ites             
!    fnspar : Computes the factor for the normalization            
!             of spherical harmonic coefficients                   
!
!   narmx=max.no.harmonic coefficients, itesx=max.degree and order
INTEGER, PARAMETER, PUBLIC :: itesx=50
INTEGER, PARAMETER, PUBLIC :: narmx=(itesx+1)**2
INTEGER, PUBLIC :: ites ! max harmonic degree used for Earth field
! name of harmonic coef/const files
CHARACTER*(50), PUBLIC :: modfile
! data on the planet we are orbiting around
! GM_{planet}, R_{planet}; 
REAL(KIND=dkind), PUBLIC :: gm_planet, r_planet
! spherical harmonics coefficients, non normalized
REAL(KIND=dkind), PUBLIC  :: harmco(narmx)

! public routines

PUBLIC :: parm10, geopot, numcoe


CONTAINS

  SUBROUTINE geopot(ites) 
!  GEOPOT: input of geopotential coefficients
! name of potential model file: modfile (shared data)
! max degree, normalization
    INTEGER, INTENT(IN) :: ites
! harmonic coefficients, normalized if armno=T
! are placed in the vector harmco
! indexes, scalar temporaries 
    INTEGER :: lll, ind, n 
! harmonic indexes and coefficients
    INTEGER :: l,m,i, ncoef 
    REAL(KIND=dkind) :: c, s
 ! unit for coefficients file
    INTEGER :: iun 
! name of potential model file
    CHARACTER (LEN=60) :: modelfile 
! dimensions check         
    ncoef=numcoe(ites)                                             
    IF(ncoef.gt.narmx)THEN 
       write(*,*) 'geopot: ncoef =',ncoef,' narmx=',narmx 
       STOP 
    ENDIF
! name of the geopotential coefficients files                            
    modelfile=modfile
    call rmsp(modelfile,lll) 
! assign planet whose harmonic expansion is used
! place radius, mass of planet in spher\_harm.mod
! should be consistent with gravity model
    r_planet=reau
    gm_planet=gmearth
! read coefficients in [modfile].coe file
    CALL filopn(iun,modelfile(1:lll)//'.coe','old') 
! default coefficients are zero
    DO  n=1,narmx 
       read(iun,*,end=2)l,m,c, s
! selection of the harmonics up to degree l=ites                        
! that is the first numcoe(ites) coefficients 
       if(l.gt.ites) CYCLE
       ind=indcoe(l,m,1)
       harmco(ind)=c/fnspar(l,m)
       IF(m.ne.0)THEN
          ind=indcoe(l,m,0)
          harmco(ind)=s/fnspar(l,m)
       ENDIF
     ENDDO 
2    CALL filclo(iun,' ')
   END SUBROUTINE geopot

  SUBROUTINE parm10 (x,ites,arm,partials) 
!  parm: vers. orbit12, revised October 1999                    
!  computes accelerations due to Mercury spherical                      
!  harmonics up to degree $l=ites > 1$.                                 
!  x= body fixed cartesian coordinates                                  
!  ites = maximum degree l                                              
!  arm(k,1,lmi)= acceleration due to harmonic lmi, k=1,3                
!  arm(k,j+1,lmi)= partial derivative of component k of accel. with resp
!      component j of x, only if partials=T}                             
!   Warning: to be optimized for performance                            
!  INPUT                                                                
    REAL(KIND=dkind), DIMENSION(3), INTENT (IN) :: x 
    INTEGER, INTENT (IN) :: ites
    LOGICAL, INTENT(IN) :: partials 
!  OUTPUT                                                               
    REAL(KIND=dkind), DIMENSION(3,4,narmx), INTENT(OUT) :: arm 
!  INCLUDE 'planet.h' moved to module shared data
    REAL(KIND=dkind), DIMENSION(narmx) :: pt,p 
    REAL(KIND=dkind), DIMENSION(itesx) :: s,c,alfa
    REAL(KIND=dkind), DIMENSION(3) :: dvy
    REAL(KIND=dkind), DIMENSION(3,3) :: ddv, dy
    REAL(KIND=dkind), DIMENSION(3,3,3) :: t 
    REAL(KIND=dkind) :: x12,x22,x32,r,r2,r3,rl,rt,rxy,rxy2,rxz2,ryz2, &
     &                    den2,rxyq                                     
    REAL(KIND=dkind) :: d,sf,cf,tgf,cl,sl,pp 
    INTEGER :: ii,lmi 
! indexes                                                               
    INTEGER :: i,j,k,l,m,n,nn                               
    SAVE 

! Warning: dipole terms do not exist                                    
    if(ites.le.1)return 
! Polar coordinates $(r, \phi, \lambda)$, radius, latitude, longitude:  
! \[ \sin(\phi)=\frac{z}{r} \quad \cos(\phi)=\frac{x^2+y^2}{r} \]          
! \[ \cos(\lambda)=\frac{x}{x^2+y^2} \quad \sin(\lambda)=\frac{y}{x^2+y^2} \]
    rt=r_planet 
    rxy2=x(1)**2+x(2)**2 
    r2=rxy2+x(3)**2 
    r=dsqrt(r2) 
    rxy=dsqrt(rxy2) 
    sf=x(3)/r 
    cf=rxy/r 
    tgf=x(3)/rxy 
    cl=x(1)/rxy 
    sl=x(2)/rxy 

! If partials have to be computed, the tensor t(i,j,k), j.ge.k,        
! of the second derivatives of $(r, \phi, \lambda)$ with respect       
! to $(x,y,z)$ has to be computed                                      
    if (partials)then 
       r3=r2*r 
       x12=x(1)**2 
       x22=x(2)**2 
       x32=x(3)**2 
       rxz2=r2-x22 
       ryz2=r2-x12 
       den2=r2*r2*rxy**3 
       rxyq=rxy2**2 
! \[ \der2parz{r}{(xyz)}=\frac{1}{r^3}\begin{pmatrix}
! y^2+z^2 & -xy & -xz \cr   
! -xy & x^2 + z^2 & -yz \cr             
! -xz &  -yz  &  x^2+y^2 \cr 
! \end{pmatrix} \]         
       t(1,1,1)=ryz2/r3 
       t(1,2,1)=-x(1)*x(2)/r3 
       t(1,3,1)=-x(1)*x(3)/r3 
       t(1,2,2)=rxz2/r3 
       t(1,3,2)=-x(2)*x(3)/r3 
       t(1,3,3)=rxy2/r3 
! \[ \der2parz{\phi}{(xyz)}=\frac{1}{r^4(x^2+y^2)^{3/2}} \times \]       
! \[ \times \begin{pmatrix} 
! zx^2(2x^2+y^2)-zy^2(z^2+y^2) &                     
! idem & idem\cr xyz[3x^2+3y^2+z^2] &                             
!                 zy^2(2y^2+x^2)-zx^2(z^2+x^2) & idem\cr                
!         x[z^2-(x^2+y^2)](x^2+y^2) &                                   
!         y[z^2-(x^2+y^2)](x^2+y^2) & -2z(x^2+y^2)^2 \cr
! \end{pmatrix} \]             
       t(2,1,1)=x(3)*(x12*(x12+rxy2)-x22*ryz2)/den2 
       t(2,2,1)=x(1)*x(2)*x(3)*(x12+2.d0*rxy2+ryz2)/den2 
       t(2,3,1)=x(1)*(x32-rxy2)*rxy2/den2 
       t(2,2,2)=x(3)*(x22*(x22+rxy2)-x12*rxz2)/den2 
       t(2,3,2)=x(2)*(x32-rxy2)*rxy2/den2 
       t(2,3,3)=-2*x(3)*rxyq/den2 
! \[ \der2parz{\lambda}{(xyz)}= \frac{1}{(x^2+y^2)^2}\begin{pmatrix}         
!    2xy & y^2-x^2 & 0 \cr y^2-x^2 & -2xy & 0 \cr 0 & 0 & 0 \cr
! \end{pmatrix} \]     
       t(3,1,1)=2*x(1)*x(2)/rxyq 
       t(3,2,1)=(x22-x12)/rxyq 
       t(3,3,1)=0.d0 
       t(3,2,2)=-2*x(1)*x(2)/rxyq 
       t(3,3,2)=0.d0 
       t(3,3,3)=0.d0 
    endif

!     $dy(i,j)=$ jacobian matrix of the polar coordinates               
!     with respect to the cartesian ones                                
! \[ dy=\derparz{(r,\phi,\lambda)}{(xyz)} \quad                        
!     \derparz{r}{x}=\frac{x}{r} \, \, etc. \]                                 
! \[ \derparz{\phi}{x}=\frac{-xz}{r^2 (x^2+y^2)} \quad                   
!    \derparz{\phi}{y}=\frac{-yz}{r^2 (x^2+y^2)} \quad                   
!    \derparz{\phi}{z}=\frac{x^2+y^2}{r^2} \]                              
! \[ \derparz{\lambda}{x}=\frac{-y}{x^2+y^2} \quad                        
!    \derparz \lambda y=\frac{x}{x^2+y^2} \]                             
    dy(1,1)=x(1)/r 
    dy(1,2)=x(2)/r 
    dy(1,3)=x(3)/r 
    dy(2,1)=-x(1)*x(3)/(r2*rxy) 
    dy(2,2)=-x(2)*x(3)/(r2*rxy) 
    dy(2,3)=rxy/r2 
    dy(3,1)=-x(2)/rxy2 
    dy(3,2)=x(1)/rxy2 
    dy(3,3)=0.d0 

!  Legendre associated functions                                        
! \[ p(lmi)=P_{l,m}(\sin\phi) \ \ for\ i=1 \]                           
! \[ pt(lmi)=P_{l,m}(\sin\phi)\begin{cases} 
! \cos(m\lambda) \ \ for\ i=1\cr     
!      \sin(m\lambda) \ \ for\ i=0,m\neq 0
! \end{cases} \]                           
    ii=(ites+1)**2 
    call legit(sf,cf,sl,cl,p,pt,s,c,ii,ites)                           
! \[ \alpha(l)=\frac {GM}{r} \frac{R_\oplus^l}{r^l} \]                   
    alfa(1)=gm_planet*rt/r2 
! \[ U=\sum_{l,m,i} J_{l,m,i} U_{l,m,i} \quad                          
!    U_{l,m,i}=\alpha(l)P_{lm}(\sin\phi)                                
!   \begin{cases} \cos(m\lambda)\, \, i=1\cr \sin(m\lambda) \, \, i=0, m\neq 0
! \end{cases} \]   
!  loop on degree l                                                     
    do 31 l=2,ites 
       alfa(l)=alfa(l-1)*rt/r 
       rl=-(l+1)/r 
!  loop on order m                                                      
       do 5 m=0,l 
!  sine and cosine terms                                                
          do 51 i=0,1 
             lmi=l*l+2*m+i 
             if(m.eq.0.and.i.eq.0)then 
!  for $m=0$ cosine only                                                  
                CYCLE
             endif

! gradient of the potential in polar coordinates                       
! \[ \derparz{U_{l,m,i}}{r}= \frac{-(l+1)}{r} \alpha(l)\, P_{lm}            
!     \begin{cases} \cos(m\lambda)\, \, i=1\cr \sin(m\lambda\, \, i=0
! \end{cases} \]          
             dvy(1)=rl*alfa(l)*pt(lmi) 
! \[ \derparz{U_{l,m,i}}{\lambda}= m \alpha(l)\, P_{lm}                 
!   \begin{cases} -\sin(m\lambda)\, \, i=1 \cr \cos(m\lambda)\, \, i=0
! \end{cases} \]          
             dvy(3)=(-1)**i*m*alfa(l)*pt(lmi+1-2*i) 
! \[ \derparz{U_{l,m,i}}{\phi}= \alpha(l) P'_{lm} \cos\phi             
!     \begin{cases} \cos(m\lambda)\, \, i=1\cr \sin(m\lambda)\, \, i=0
! \end{cases} \]         
! \[ P'_{lm}\cos\phi=P_{l,m+1}-m \tan\phi\, P_{lm} \]                    
! Warning: $P_{l,l+1}=0$                                                
             if(m.eq.l)then 
                dvy(2)=alfa(l)*(-m*tgf*pt(lmi)) 
             elseif(m.eq.0)then 
                dvy(2)=alfa(l)*p(lmi+2) 
             else 
                pp=p(l*l+2*m+3)*(s(m)*(1-i)+c(m)*i) 
                dvy(2)=alfa(l)*(pp-m*tgf*pt(lmi)) 
             endif

!  gradient of the potential in cartesian coordinates                   
! \[ \derparz{U}{(xyz)}=                                                
!       \derparz{U}{(r\phi\lambda)} \derparz{(r\phi\lambda)}{(xyz)} \]   
             DO j=1,3 
                d=0.d0 
                DO k=1,3 
                   d=d+dvy(k)*dy(k,j)
                ENDDO 
                arm(j,1,lmi)=d
             ENDDO
             if(partials)then 
!  ddv(i,j): matrix of second derivatives of U                          
! \[ \der2parz{U_{l,m,i}}{r} = \frac{(l+1)(l+2)}{r^2}\alpha(l)\,          
!       P_{lm} \begin{cases} \cos(m\lambda)\cr                                 
!     \sin(m\lambda)
! \end{cases} \]                                                 
                ddv(1,1)=-(l+2)*dvy(1)/r 
! \[ \frac{\partial^2 U_{l,m,i}}{\partial r \partial\phi}=              
!     \frac{-(l+1)}r \alpha(l) P'_{lm} \begin{cases} \cos(m\lambda)\cr         
!     \sin(m\lambda)
! \end{cases}= \frac{-(l+1)}{r}\derparz{U_{l,m,i}}{\phi} \]    
                ddv(1,2)=-(l+1)*dvy(2)/r 
                ddv(2,1)=ddv(1,2) 
! \[ \frac{\partial^2 U_{l,m,i}}{\partial r \partial\lambda}=           
!    \frac{-(l+1)}{r} \alpha(l) m\, P_{lm}\begin{cases}                        
!    -\sin(m\lambda) \cr \cos(m\lambda)
! \end{cases}                                
!     = \frac{-(l+1)}r \derparz {U_{l,m,i}}{\lambda} \]    
                ddv(1,3)=-(l+1)*dvy(3)/r 
                ddv(3,1)=ddv(1,3) 
! \[ \der2parz{U_{l,m,i}}{\phi}= \alpha(l) \der2parz{P_{lm}(\sin\phi)}{
!     \begin{cases} \cos(m\lambda) \cr \sin(m\lambda)
! \end{cases}} \]                         
! \[ \der2parz{P_{lm}(\sin\phi)}{\phi}=\tan\phi P_{l,m+1} +             
!  [m^2 \sec^2\phi -m\tan^2\phi-l(l+1)] P_{lm}(\sin\phi) \]              
                ddv(2,2)=tgf*dvy(2)+alfa(l)*pt(lmi)*((m/cf)**2-l*(l+1)) 
! \[ \frac{\partial^2 U_{l,m,i}}{\partial\lambda \partial\phi}=         
!    m \alpha(l) P'_{lm} \cos\phi \begin{cases} 
! -\sin(m\lambda) \, \, i=1\cr     
!      \cos(m\lambda) \, \, i=0
! \end{cases} \]                                       
! \[ P'_{lm}\cos\phi=P_{l,m+1}-m \tan\phi\ P_{lm} \]                    
! Warning: $P_{l,l+1}=0$                                                
                if(m.eq.l)then 
                   ddv(2,3)=m*alfa(l)*(-tgf*m*pt(lmi+1-2*i))*(-1)**i 
                elseif(m.eq.0)then 
                   ddv(2,3)=0.d0 
                else 
                   pp=p(l*l+2*m+3)*((1-i)*c(m)+i*s(m)) 
                   ddv(2,3)=m*alfa(l)*(pp-                               &
     &                     tgf*m*pt(lmi+1-2*i))*(-1)**i                 
                endif
                ddv(3,2)=ddv(2,3) 
! \[ \der2parz U\lambda= -m^2 \alpha(l)\  P_{lm} \begin{cases} \cos(m\lambda)\,
!     \sin(m\lambda)
! \end{cases} \]                                                 
                ddv(3,3)=-m*m*alfa(l)*pt(lmi) 
!  matrix of partials of the accelerations with respect to coordinates  
! \[ \der2parz{U}{(xyz)}=\bigl[\derparz{(r\phi\lambda)}{(xyz)}\bigl]^T  
!   \der2parz{U}{(r\phi\lambda)} \derparz{(r\phi\lambda)}{(xyz)} +      
!   \derparz{U}{(r\phi\lambda)} \der2parz{(r\phi\lambda)}{(xyz)} \]     
                DO k=1,3 
                   DO j=1,k 
                      d=0.d0 
                      DO n=1,3 
                         d=d+t(n,k,j)*dvy(n) 
                         DO nn=1,3 
                           d=d+ddv(n,nn)*dy(n,k)*dy(nn,j)
                         ENDDO 
                      ENDDO
                      arm(j,k+1,lmi)=d 
                      arm(k,j+1,lmi)=d
                   ENDDO
                ENDDO
             endif
!  end loop on i=0,1                                                    
51        ENDDO
!  end loop on m=0,l                                                    
5      ENDDO
!  end loop on l=2,ites                                                 
31  ENDDO
  END subroutine parm10
      
! \bighorline
  SUBROUTINE legit(sf,cf,sl,cl,p,pt,s,c,ii,i) 
!\comments
!  {\bf legit}                                                          
!   recursive computation of Legendre polynomials                       
!   and Legendre associated functions, unnormalised,                    
!   up to degree $l=i$                                                  
! {\obeylines                                                           
!     INPUT:                                                            
!     $sf=\sin\phi\quad cf=\cos\phi$                                   
!     $sl=\sin\lambda\quad cf=\cos\lambda $                            
!     OUTPUT:                                                           
!     p,pt=vectors of length $ii=(i+1)^2$:                              
!     $p(l^2+2m+1)=P_{l,m}(\sin\phi)$ Legendre associated function      
!     $pt(l^2+2m+1)=P_{l,m}(\sin\phi)\cos(m\lambda)$                    
!     $pt(l^2+2m)=P_{l,m}(\sin\phi)\sin(m\lambda) \ \ \ m\neq 0$} 
! \interface
!  INPUT                                                                
    REAL(KIND=dkind), INTENT (IN) :: sf,sl,cf,cl 
    INTEGER, INTENT (IN) :: i,ii 
!  OUTPUT              
    REAL(KIND=dkind), DIMENSION(i), INTENT (OUT) :: c,s 
    REAL(KIND=dkind), DIMENSION(ii), INTENT (OUT) :: p,pt
! \endint
    INTEGER :: j 
!  indexes                                                              
    INTEGER :: l,m 
    SAVE

!                                                                       
!  initialization of recursion                                          
! \[ P_{00}(\sin\phi)=1 \quad P_{10}(\sin\phi)=\sin\phi \quad         
!    P_{11}(\sin\phi)=\cos\phi \]                                       
    p(1)=1.d0 
    p(2)=sf 
    p(4)=cf 
! \[ c(m)=\cos(m\lambda) \quad s(m)=\sin(m\lambda) \]                  
    c(1)=cl 
    s(1)=sl 
!  recursion formulae:                                                  
    DO l=2,i 
! for Legendre polynomials ($m=0$):                                     
! \[ P_{l}(\sin\phi)=\frac 1l [(2l-1)\sin\phi P_{l-1}-(l-1)P_{l-2} \]
       p(l*l+1)=((2*l-1)*sf*p((l-1)*(l-1)+1)-(l-1)*p((l-2)*(l-2)+1))/l 
!  for $m=l$:                                                           
! \[ P_{ll}(\sin\phi)= (2l-1)\cos\phi P_{l-1,l-1} \]                   
       p((l+1)*(l+1))=(2*l-1)*cf*p(l*l) 
!  for $m=l-1$:                                                         
! \[ P_{l,l-1}(\sin\phi)=(2l-1)\cos\phi P_{l-1,l-2} \]                  
       p((l+1)*(l+1)-2)=(2*l-1)*cf*p(l*l-2) 
! for $1\leq m \leq l-2\ ,\ l\geq 2$:                                   
! \[ P_{l,m}(\sin\phi)=P_{l-2,m}+(2l-1)P_{l-1,m-1}\cos\phi \]           
       j=l-2 
       IF(j.ge.0)THEN 
          DO m=1,j 
             p(l*l+2*m+1)=p((l-2)*(l-2)+2*m+1)+(2*l-1)*cf                  &
     &                 *p((l-1)*(l-1)+2*m-1)
          ENDDO
       ENDIF
!  addition formulae for sines and cosines of multiples of longitude    
       s(l)=s(l-1)*cl+c(l-1)*sl 
       c(l)=c(l-1)*cl-s(l-1)*sl 
    ENDDO
    DO l=1,i 
!  zonal case                                                           
       pt(l*l+1)=p(l*l+1) 
!  sectorial and tesseral cases                                         
       DO m=1,l 
          j=l*l+2*m 
          pt(j)=p(j+1)*s(m) 
          pt(j+1)=p(j+1)*c(m)
       ENDDO
    ENDDO
  END SUBROUTINE legit
! \bighorline                                               
  INTEGER FUNCTION indcoe(l,m,i) 
! \comments
! \begin{verbatim}
! Harmonics addressing routines: $lmi=l^2+2m+i$                     
! Warning: list of coefficients begins with $\Delta C_{00}$;          
! space is left even for the dipole terms, never used
! \end{verbatim}                
! {\bf indcoe}: address of harmonics (l,m,i)      
! \interface
    INTEGER, INTENT (IN) :: i,l,m 
! \endint
    if(l.lt.0)stop ' **** indcoe: l < 0 ****' 
    if(m.lt.0)stop ' **** indcoe: m < 0 ****' 
    if(m.gt.l)stop ' **** indcoe: m > l ****' 
    if(m.eq.0)then 
       if(i.ne.1)stop ' **** indcoe: i .ne. 1 (m=0) ****' 
    else 
       if(i.ne.1.and.i.ne.0)stop ' **** indcoe: i = ? (m>0) ****' 
    endif
    indcoe=l*l+2*m+i 
  END FUNCTION indcoe

! \bighorline
  SUBROUTINE namcoe(lmi,l,m,i) 
! \comments
!  {\bf namcoe} inverse of indcoe 
! \interface
    INTEGER, INTENT (IN) :: lmi 
    INTEGER, INTENT (OUT) :: l,m,i 
! \endint
    REAL(KIND=dkind) :: sl 
    sl=lmi 
    l=dsqrt(sl) 
! fix for case of $S_l0$ not existing                                   
    IF(l*l-lmi.eq.0)l=l-1 
    m=(lmi-l*l)/2 
    i=lmi-l*l-2*m 
  END SUBROUTINE namcoe

! \bighorline
    INTEGER FUNCTION numcoe(ites) 
! \comments
! {\bf numcoe} number of coefficients up to degree ites   
! \interface
    INTEGER, INTENT (IN) :: ites
! \endint 
    numcoe=(ites+1)**2 
  END FUNCTION numcoe

! \bighorline        
  REAL(KIND=dkind) FUNCTION fnspar(l,m)
! \comments
! \begin{verbatim}
! F N S P A R : Factor for the normalization of spherical harmonic 
! coefficients 
!                                                                       
! INPUT:   
! L : Degree                                        
! M : Order                                         
! OUTPUT:
! FNSPAR: Ratio between normalized and unnormalized coefficient:
! Cnorm(l,m)=fnspar(l,m)*C(l,m)
! \end{verbatim}                  
! \interface                                                              
    INTEGER, INTENT (IN) :: l,m ! degree, order
! \endint
    INTEGER :: ne 
    REAL(KIND=dkind) :: f 
! indexes                                                               
    INTEGER :: j 
!                                                                       
    fnspar=1.d0 
    ne=0 
    if(m.ne.0)then 
       do j=l-m+1,l+m 
          fnspar=fnspar*j 
          if(fnspar.gt.1.d10)then 
             fnspar=fnspar/1.d10 
             ne=ne+1 
          end if
       enddo
       fnspar=fnspar/2.d0 
    end if
    fnspar=fnspar/(2*l+1) 
    fnspar=sqrt(fnspar) 
    if(ne.ne.0)then 
       f=1.d5**ne 
       fnspar=f*fnspar 
    end if
  END FUNCTION fnspar

    
END MODULE spher_harm

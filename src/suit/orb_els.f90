! ================ LIBRARY orb_els===============               
! ORBITAL ELEMENTS:                                                     
! CONTAINS                                                              
! SUBROUTINES                                                           
!   COORDINATE CHANGE (to be replaced by orbit_elements)                                                   
!   coocha      general purpose coordinate/element change (with and w/o 
! prop2b        2-body propagator, with derivatives (2nd derivatives??) 
! *carequ	cartesian to equinoctal transformation                        
! *equcar	equinoctal to cartesian transformation                        
! *kepequ	keplerian to equinoctal transformation                        
! *equkep	equinoctal to keplerian transformation                        
! *eqpequ	equinoctal polar to equinoctal transformation                 
! *equeqp	equinoctal to equinoctal polar transformation                 
! *ekensd	keplerian to equinoctal transformation                        
! *kepcar       fast conversion directly from keplerian to cartesian 

! HEADERS and MODULES                                                   
! orb_els.o: \                                                          
!	../include/parcmc.h \ public                                                   
!	comorb.h \        private                                                     
!	reference_systems.o                                                   
!                                                                               
!                                                                       
! ===================================================================   
! COOCHA-COODER                                                         
! ===================================================================   
!   general purpose coordinate to/from                                  
!   elements change                                                     
!                                                                       
!  specifically designed to handle correctly low                        
!  eccentricity and low inclination orbits;                             
!  polar orbits are supported without derivatives                       
!                                                                       
!  limitations: only elliptic orbits are handled in this                
!  version; singularity for 180 degrees inclination                     
!                                                                       
!  coordinate types implememted:                                        
!    'CAR' : cartesian positions and velocities in a vector             
!            of length 6                                                
!    'EQU' : equinoctal elements, see carequ, in a vector of            
!            length 6; the mean motion enne is also provided            
!    'EQP' : equinoctal elements for highly inclined orbit              
!    'KEP' : classical keplerian elements,  see kepequ, in a            
!             vector of length 6; singular for                          
!             0 eccentricity, 0 and 180 degrees inclination;            
!                                                                       
!    INPUT:  x(6), coox (string defining the coord. type of x)          
!            gm=G x Mass of the Sun (anyway the central mass, e.g.      
!                 mass of Sun = mass of planet for massive palnet)      
!            cooy (string defining the coord. type)                     
!                                                                       
!    OUTPUT: y(6)                                                       
!            enne=mean motion; WARNING: enne is not computed            
!            for identical transformation                               
! ====================================================================  
! COOCHA                                                                
! ====================================================================  
! in this simpler version the elements are computed without derivatives 
! ===============INTERFACE========================================      
      subroutine coocha(x,coox,gm,y,cooy,enne) 
      implicit none 
! ================= input/output ================                       
      double precision x(6),y(6),gm,enne 
      character*3 coox,cooy 
! =============END INTERFACE====================                        
! ================ workspace ============                               
      double precision z(6) 
! ================ loop indexes ==============                          
      integer j 
! ================ rounding off problem ===========                     
      double precision eps 
!  required in Kepler equation and controls for singular cases
!*********************************                                      
!  static memory no                                              
!*********************************                                      
! error return with enne=0                                              
      enne=0.d0 
!  dummy case                                                           
      if(coox.eq.cooy)then 
         do 1 j=1,6 
    1      y(j)=x(j) 
         return 
      endif 
!  input is interpreted according to input type coox,                   
!  and anyway transformed to equinoctal elements z; enne is             
!  also computed in any case apart from identical transformation        
      if(coox.eq.'EQU')then 
         do 2 j=1,6 
    2      z(j)=x(j) 
         if(x(1).gt.0.d0)then 
            enne=sqrt(gm/x(1)**3) 
         else 
            enne=0.d0 
         endif 
      elseif(coox.eq.'KEP')then 
         call kepequ(x,z) 
         if(x(1).gt.0.d0)then 
            enne=sqrt(gm/x(1)**3) 
         else 
            enne=0.d0 
         endif 
      elseif(coox.eq.'CAR')then 
         call carequ(x,gm,z,enne) 
         IF(enne.eq.0.d0)RETURN 
      elseif(coox.eq.'EQP')then 
         call eqpequ(x,z) 
         if(x(1).gt.0.d0)then 
            enne=sqrt(gm/x(1)**3) 
         else 
            enne=0.d0 
         endif 
      else 
         write(*,*)'**** coocha: in coord.type ',coox,' unknown ****' 
         stop 
      endif 
!  transformation to the output type cooy 
      eps=1.d2*epsilon(1.d0) 
      if(cooy.eq.'EQU')then 
         do 3 j=1,6 
    3      y(j)=z(j) 
      elseif(cooy.eq.'KEP')then 
!  this is a potentially singular case; the control for negligible      
!  eccentricity and inclination is set to $100 \times$ rounding off     
         call equkep(z,eps,y) 
      elseif(cooy.eq.'EQP')then 
!  this case is singular only for zero inclination                      
         call equeqp(z,eps,y) 
      elseif(cooy.eq.'CAR')then 
!  control for convergence in kepler equation is set to $100 \times$    
!  rounding off                                                         
         call equcar(z,gm,eps,y) 
      else 
         write(*,*)'**** coocha: out coord.type ',cooy,' unknown ****' 
         stop 
      endif 
      return 
      END SUBROUTINE coocha                                         
! ======================================================                
! COODER                                                                
! ======================================================                
!   general purpose coordinate to/from                                  
!   elements change: version with partial derivatives                   
!   definitions and input as for coocha (but no equinoctal polar)       
!   in output, derpar(6,6) = matrix of partial derivatives              
! ==============INTERFACE============================                   
      subroutine cooder(x,coox,gm,y,cooy,enne,derpar) 
      implicit none 
! ================= input/output ================                       
      double precision x(6),y(6),gm,enne,derpar(6,6) 
      character*3 coox,cooy 
! =============END INTERFACE====================                        
! ================ workspace ============                               
      double precision z(6),derws(6,6),derws2(6,6),w(6),ddxde(3,6,6) 
! =============== scalars ===============                               
      double precision t0,det 
! ================ loop indexes ==============                          
      integer j 
! ================ derivatives, inversion control =======               
      integer ider,ising 
! ================ rounding off problem ===========                     
      double precision eps 
!  required in Kepler equation and controls for singular cases   
!*********************************                                      
!  static memory no 
!*********************************                                      
! error return with enne=0                                              
      enne=0.d0 
!  dummy case                                                           
      if(coox.eq.cooy)then 
         do 1 j=1,6 
    1      y(j)=x(j) 
         call eye(6,derpar) 
         return 
      endif 
!  input is interpreted according to input type coox,                   
!  and anyway transformed to equinoctal elements z; enne is             
!  also computed in any case apart from identical transformation        
      if(coox.eq.'EQU')then 
         do 2 j=1,6 
    2      z(j)=x(j) 
         if(x(1).gt.0.d0)then 
            enne=sqrt(gm/x(1)**3) 
         else 
            enne=0.d0 
         endif 
         call eye(6,derws) 
      elseif(coox.eq.'KEP')then 
         call ekensd(x,z,derws) 
         if(x(1).gt.0.d0)then 
            enne=sqrt(gm/x(1)**3) 
         else 
            enne=0.d0 
         endif 
      elseif(coox.eq.'CAR')then 
         call carequ(x,gm,z,enne) 
         IF(enne.eq.0.d0)RETURN 
         t0=0.d0 
         ider=1 
         call prop2b(t0,z,t0,w,gm,ider,derws,ddxde) 
         call matin(derws,det,6,0,6,ising,1) 
      elseif(coox.eq.'EQP')then 
         write(*,*)' partial derivatives for EQP not implemented' 
         stop 
      else 
         write(*,*)'**** coocha: in coord.type ',coox,' unknown ****' 
         stop 
      endif 
!  transformation to the output type cooy                               
         eps=1.d2*epsilon(1.d0) 
      if(cooy.eq.'EQU')then 
         do 3 j=1,6 
    3      y(j)=z(j) 
         call eye(6,derws2) 
      elseif(cooy.eq.'KEP')then 
!  this is a potentially singular case; the control for negligible      
!  eccentricity and inclination is set to $100 \times$ rounding off     
         call equkep(z,eps,y) 
         call ekensd(y,w,derws2) 
! ***    write(*,*)(w(ii)-z(ii),ii=1,6)                                 
         call matin(derws2,det,6,0,6,ising,1) 
      elseif(cooy.eq.'EQP')then 
         write(*,*)' partial derivatives for EQP not implemented' 
         stop 
      elseif(cooy.eq.'CAR')then 
!  control for convergence in kepler equation is set to $100 \times$    
!  rounding off                                                         
         t0=0.d0 
         ider=1 
         call prop2b(t0,z,t0,y,gm,ider,derws2,ddxde) 
         call equcar(z,gm,eps,w) 
! ***    write(*,*)(w(ii)-y(ii),ii=1,6)                                 
      else 
         write(*,*)'**** coocha: out coord.type ',cooy,' unknown ****' 
         stop 
      endif 
!  multiplication of jacobian matrices to get the jacobian of the compos
      derpar=MATMUL(derws,derws2) ! call mulmat(derws,6,6,derws2,6,6,derpar) 
      return 
   END SUBROUTINE cooder
!                                                                       
!  PROP2B                                                               
!                                                                       
!  Task:       Solves the two body problem in equinoctal coordinates    
!                                                                       
!  Features:   the Newton's method is exploited in the interval         
!              [W, 2PI + W] where W is the peri-planet longitude:       
!              W = w_p + RAAN                                           
!                                                                       
!  Input arguments:                                                     
!              E         =     equinoctal orbital elements              
!                               e=(a,h,k,p,q,lambda)                    
!              T0        =     epoch time                               
!              T1        =     prediction time                          
!              GM        =     gravitational constant of the system GM  
!              IDER       =     flag for derivatives option :           
!                                0 = only position and velocites        
!                                1 = first derivatives                  
!                                2 = first and second partial derivative
!  Output arguments:                                                    
!              X         =     position and velocity components         
!                               in absolute cartesian coordinates       
!              DXDE      =     first derivatives of position vector x   
!                               with respect to elements                
!              DDXDE     =     second derivatives of position vector x  
!                               with respect to elements                
!                                                                       
!****************                                                       
!   static memory not required                                          
!****************                                                       
      subroutine prop2b(t0,e,t1,x,gm,ider,dxde,ddxde)
      USE fund_const 
      implicit none 
      double precision e(6),x(6) 
      double precision f(3),g(3),w(3) 
      double precision dxde(6,6),ddxde(3,6,6) 
! scalars                                                               
      double precision t0,t1,gm 
      integer ider 
! iteration maximum                                                     
      integer iter 
      parameter (iter=25) 
! scalar temporaries                                                    
      double precision enne, pml,ecc2,eps,pol,princ 
      double precision sinel,cosel,rf,rdf,del,el,beta,r,upq,            &
     &  xe,ye,coe,xpe,ype,tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7,tmp8,      &
     &  tmp9,tmp10,tmp11,tmp12,dx1de2,dx2de2,dx1de3,dx2de3,dx4de2,      &
     &  dx5de2,dx4de3,dx5de3                                            
      integer j
!  Mean motion                                                          
      enne=sqrt(gm/e(1)**3) 
!  Mean longitude PML at time T1                                        
      pml=e(6)+enne*(t1-t0) 
!  Absolute peri-planet longitude POL at epoch T0                       
      ecc2=e(2)**2+e(3)**2 
      eps=epsilon(1.d0)*1.d2 
      if(ecc2.lt.eps)then 
         pol=0.d0 
      elseif(ecc2.ge.1.d0)then 
         write(*,*)' dcc.ge.1, ecc**2=',ecc2 
         write(*,*)e 
         stop 
      else 
         pol=atan2(e(2),e(3)) 
         pol=princ(pol) 
      endif 
!  Mean longitude restriction to [POL, POL + 2*PIGR]                    
!  (Newton's method continuity conditions)                              
      pml = princ(pml) 
      if (pml.lt.pol)then 
         pml = pml + dpig 
      endif 
!  Newton's method for motion equation solution:                        
!  R(F,lambda) = F - ksinF + hcosF - lambda = 0                         
!  search for F for a given lambda within [POL, POL + 2PIGR]            
      el = pig+pol 
      do 70 j=1,iter 
             sinel=sin(el) 
             cosel=cos(el) 
             rf = el - e(3)*sinel + e(2)*cosel - pml 
             rdf = 1.d0 - e(3)*cosel - e(2)*sinel 
             del = - rf/rdf 
             el = el + del 
             if (abs(del).lt.eps) goto 100 
   70 continue 
      write(*,*)' Too many iter. in newton, iter=',iter,' del=',del 
      write(*,*) ' eq,eps ',e, eps 
!      stop                                                             
!  Computation of position and velocity on the orbit plane              
!  in equinoctal cartesian coordinates (f,g,w)                          
  100 beta=1.d0/(1.d0+sqrt(1.d0-ecc2)) 
      xe=e(1)*((1.d0-beta*e(2)**2)*cosel+e(2)*e(3)*beta*sinel-e(3)) 
      ye=e(1)*((1.d0-beta*e(3)**2)*sinel+e(2)*e(3)*beta*cosel-e(2)) 
!  Equinoctal reference frame                                           
      upq=1.d0+e(4)**2+e(5)**2 
      f(1)=(1.d0-e(4)**2+e(5)**2)/upq 
      f(2)=2.d0*e(4)*e(5)/upq 
      f(3)=-2.d0*e(4)/upq 
      g(1)=2.d0*e(4)*e(5)/upq 
      g(2)=(1.d0+e(4)**2-e(5)**2)/upq 
      g(3)=2.d0*e(5)/upq 
!  Conversion from equinoctal to absolute coordinates                   
      call lincom(f,xe,g,ye,x) 
!  Computation of velocities                                            
      coe=enne*e(1)**2/sqrt(xe**2+ye**2) 
      xpe=coe*(e(2)*e(3)*beta*cosel-(1.d0-beta*e(2)**2)*sinel) 
      ype=coe*((1.d0-beta*e(3)**2)*cosel-e(2)*e(3)*beta*sinel) 
      call lincom(f,xpe,g,ype,x(4)) 
!  Computation of partials if required                                  
      if(ider.lt.1)return 
!  Equinoctal reference frame, third vector                             
      w(1)=2.d0*e(4)/upq 
      w(2)=-2.d0*e(5)/upq 
      w(3)=(1.d0-e(4)**2-e(5)**2)/upq 
!  Computation of same temporary variables                              
      r=sqrt(xe**2+ye**2) 
      tmp1=pml-el 
      tmp2=beta+e(2)**2*beta**3/(1.d0-beta) 
      tmp3=e(2)*e(3)*beta**3/(1.d0-beta) 
      tmp4=beta*e(2)-sinel 
      tmp5=beta*e(3)-cosel 
      tmp6=beta+e(3)**2*beta**3/(1.d0-beta) 
      tmp7=1.d0-r/e(1) 
      tmp8=sinel-e(2) 
      tmp9=cosel-e(3) 
      tmp10=e(1)*cosel/r 
      tmp11=e(1)*sinel/r 
      tmp12=enne*e(1)**2/r 
!  Computation of derivatives of position vector w. r. to the elements  
      dxde(1,1)=(x(1)-3.d0*x(4)*(t1-t0)/2.d0)/e(1) 
      dxde(2,1)=(x(2)-3.d0*x(5)*(t1-t0)/2.d0)/e(1) 
      dxde(3,1)=(x(3)-3.d0*x(6)*(t1-t0)/2.d0)/e(1) 
      dx1de2=-e(1)*(tmp1*tmp2+e(1)*cosel*tmp4/r) 
      dx2de2=e(1)*(tmp1*tmp3-1.d0+e(1)*cosel*tmp5/r) 
      call lincom(f,dx1de2,g,dx2de2,dxde(1,2)) 
      dx1de3=-e(1)*(tmp1*tmp3+1.d0-e(1)*sinel*tmp4/r) 
      dx2de3=e(1)*(tmp1*tmp6-e(1)*sinel*tmp5/r) 
      call lincom(f,dx1de3,g,dx2de3,dxde(1,3)) 
      dxde(1,4)=2.d0*(e(5)*(ye*f(1)-xe*g(1))-xe*w(1))/upq 
      dxde(2,4)=2.d0*(e(5)*(ye*f(2)-xe*g(2))-xe*w(2))/upq 
      dxde(3,4)=2.d0*(e(5)*(ye*f(3)-xe*g(3))-xe*w(3))/upq 
      dxde(1,5)=2.d0*(e(4)*(-ye*f(1)+xe*g(1))+ye*w(1))/upq 
      dxde(2,5)=2.d0*(e(4)*(-ye*f(2)+xe*g(2))+ye*w(2))/upq 
      dxde(3,5)=2.d0*(e(4)*(-ye*f(3)+xe*g(3))+ye*w(3))/upq 
      dxde(1,6)=x(4)/enne 
      dxde(2,6)=x(5)/enne 
      dxde(3,6)=x(6)/enne 
!  Computation of derivatives of velocity vector w. r. to the elements  
      dxde(4,1)=-(x(4)-3.d0*gm*x(1)*(t1-t0)/r**3)/(2.d0*e(1)) 
      dxde(5,1)=-(x(5)-3.d0*gm*x(2)*(t1-t0)/r**3)/(2.d0*e(1)) 
      dxde(6,1)=-(x(6)-3.d0*gm*x(3)*(t1-t0)/r**3)/(2.d0*e(1)) 
      dx4de2=tmp12*(tmp7*tmp2+e(1)**2*tmp8*tmp4/r**2+tmp10*cosel) 
      dx5de2=-tmp12*(tmp7*tmp3+e(1)**2*tmp8*tmp5/r**2-tmp10*sinel) 
      call lincom(f,dx4de2,g,dx5de2,dxde(4,2)) 
      dx4de3=tmp12*(tmp7*tmp3+e(1)**2*tmp9*tmp4/r**2-tmp11*cosel) 
      dx5de3=-tmp12*(tmp7*tmp6+e(1)**2*tmp9*tmp5/r**2+tmp11*sinel) 
      call lincom(f,dx4de3,g,dx5de3,dxde(4,3)) 
      dxde(4,4)=2.d0*(e(5)*(ype*f(1)-xpe*g(1))-xpe*w(1))/upq 
      dxde(5,4)=2.d0*(e(5)*(ype*f(2)-xpe*g(2))-xpe*w(2))/upq 
      dxde(6,4)=2.d0*(e(5)*(ype*f(3)-xpe*g(3))-xpe*w(3))/upq 
      dxde(4,5)=2.d0*(e(4)*(-ype*f(1)+xpe*g(1))+ype*w(1))/upq 
      dxde(5,5)=2.d0*(e(4)*(-ype*f(2)+xpe*g(2))+ype*w(2))/upq 
      dxde(6,5)=2.d0*(e(4)*(-ype*f(3)+xpe*g(3))+ype*w(3))/upq 
      dxde(4,6)=-enne*e(1)**3*x(1)/r**3 
      dxde(5,6)=-enne*e(1)**3*x(2)/r**3 
      dxde(6,6)=-enne*e(1)**3*x(3)/r**3 
!  Computation of second derivatives if required                        
      if(ider.lt.2)return 
      write(*,*) 'second derivatives not supported in this version' 
      stop 
   END SUBROUTINE prop2b
! ======================================================                
!   {\bf carequ}: coordinate change from cartesian to                   
!   equinoctal elements                                                 
!                                                                       
!     x= position and velocity, cartesian coordinates                   
!                                                                       
!     gm= $G \times Mass$                                               
!                                                                       
!         eq(1)= a                                                      
!                                                                       
!         eq(2)=ecc*sin(dig)                                            
!                                                                       
!         eq(3)=ecc*cos(dig) ;  dig=longitude pericentre                
!                                                                       
!         eq(4)=tgim*cos(omega) ; tgim=tg(i/2)                          
!                                                                       
!         eq(5)=tgim*sin(omega) ; omega=longitude ascending node        
!                                                                       
!         eq(6)=mean longitude                                          
!                                                                       
!         enne= mean motion                                             
!                                                                       
!  if the energy is not negative, eq(6) is set equal to true longitude, 
!  eq(1)=negative a and enne=0                                          
!                                                                       
!  bugs: zero angular momentum and/or energy are handled in such        
!  a way that zero divide are avoided, but overflows can occur          
!                                                                       
!  this definition of elements is singular only for a planar            
!  retrograde orbit and parabolic orbits                                
!=======================================================                
      subroutine carequ(x,gm,eq,enne) 
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) ::  x(6),gm
      DOUBLE PRECISION, INTENT(OUT) :: eq(6),enne 
! =========END INTERFACE=====================
      DOUBLE PRECISION, DIMENSION(3) :: ang,f,g,vlenz 
      DOUBLE PRECISION vel2,r,gei,d,ainv,cosf,sinf,ecc2,rad,beta,chk,ch,ck,fe,x2,y2
      DOUBLE PRECISION princ,vsize,prscal ! functions
      INTEGER i
! error return when enne is zero                                        
      enne=0.d0 
!  radius and velocity squared                                          
      vel2=prscal(x(4),x(4)) 
      r=vsize(x) 
!  angular momentum                                                     
      call prvec(x(1),x(4),ang) 
!   zero divide occurs for zero angular momentum                        
      gei=vsize(ang) 
      IF(gei.eq.0.d0)THEN 
         write(*,*) '****** carequ: zero angular momentum ******' 
         return 
      ENDIF 
!  angular momentum unit vector, Lenz vector                            
      call prvec(x(4),ang,vlenz) 
      do 1 i=1,3 
        vlenz(i)=vlenz(i)/gm-x(i)/r 
    1   ang(i)=ang(i)/gei 
!   zero divide occurs for inclination of 180 degrees                   
      d=1.d0+ang(3) 
      IF(d.eq.0.d0)THEN 
         WRITE(*,*) '****** carequ: 180 deg. inclination ******' 
         RETURN 
      ENDIF 
!  unit vectors of the equinoctal reference system (Broucke and         
!  Cefola 1972, CM 5, 303--310) are f, g, ang                           
      f(1)=1.d0-ang(1)**2/d 
      f(2)=-ang(1)*ang(2)/d 
      f(3)=-ang(1) 
      call prvec(ang,f,g) 
!  elements related to eccentricity and inclination                     
      eq(2)=prscal(vlenz,g) 
      eq(3)=prscal(vlenz,f) 
      eq(4)=ang(1)/d 
      eq(5)=-ang(2)/d 
!     tgim1=dsqrt(eq(4)**2+eq(5)**2)                                    
!  test on energy                                                       
      ainv=2/r-vel2/gm 
      if(ainv.eq.0.d0)then 
! eq(1) is q                                                            
         eq(1)=(r+prscal(vlenz,x))/(1+vsize(vlenz)) 
         enne=0 
         cosf=prscal(x,f) 
         sinf=prscal(x,g) 
         eq(6)=atan2(sinf,cosf) 
         eq(6)=princ(eq(6)) 
         return 
!         stop '****** carequ: parabolic orbit ******'                  
      elseif(ainv.lt.0.d0)then 
         eq(1)=1.d0/ainv 
         enne=0 
         cosf=prscal(x,f) 
         sinf=prscal(x,g) 
         eq(6)=atan2(sinf,cosf) 
         eq(6)=princ(eq(6)) 
         return 
      endif 
!   semimajor axis and mean motion                                      
      eq(1)=1.d0/ainv 
      enne=dsqrt(gm/eq(1)**3) 
!   mean longitude from non--singular Kepler equation                   
      ecc2=eq(2)**2+eq(3)**2 
      rad=dsqrt(1.d0-ecc2) 
      beta=1.d0/(1.d0+rad) 
      chk=eq(2)*eq(3)*beta 
      ch=1.d0-eq(2)**2*beta 
      ck=1.d0-eq(3)**2*beta 
      x2=prscal(x,f) 
      y2=prscal(x,g) 
      cosf=eq(3)+(ck*x2-chk*y2)/(eq(1)*rad) 
      sinf=eq(2)+(ch*y2-chk*x2)/(eq(1)*rad) 
!   eccentric longitude                                                 
      fe=datan2(sinf,cosf) 
      eq(6)=fe+eq(2)*cosf-eq(3)*sinf 
!   reduction to principal value                                        
      eq(6)=princ(eq(6)) 
      return 
      END subroutine carequ                                          
! ======================================================                
!   {\bf equcar}: coordinate change from                                
!   equinoctal elements to cartesian                                    
!                                                                       
!     x= position and veloctiy, cartesian coordinates                   
!                                                                       
!     gm= $G \times Mass$                                               
!                                                                       
!     eps= convergence control for non--singular Kepler equation        
!                                                                       
!     eq: equinoctal elements, see carequ                               
!                                                                       
!         enne= mean motion                                             
!                                                                       
!  bugs: if  eq(1) is negative, the hyperbolic case is not handled in   
!  this version                                                         
!                                                                       
!  this definition of elements is singular only for a planar            
!  retrograde orbit and parabolic orbits                                
!=======================================================                
      subroutine equcar(eq,gm,eps,x) 
!  ==================================================================== 
      USE fund_const
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) ::  eq(6),gm,eps
      DOUBLE PRECISION, INTENT(OUT) :: x(6)
! =========END INTERFACE=====================
      DOUBLE PRECISION, DIMENSION(3) :: f,g 
      DOUBLE PRECISION cosf,sinf,ecc2,rad,beta,chk,ch,ck,fe,x2,y2,df,xp2,yp2
      DOUBLE PRECISION tgim2,ecc,opwz,dig, tlong,princ,enne,de
      INTEGER i
!   test for hyperbolic orbit                                           
      if(eq(1).le.0.d0)then 
         stop '****** equcar: hyperbolic/parabolic orbit ******' 
      endif 
!  non--singular intermediate variables                                 
      ecc2=eq(2)**2+eq(3)**2 
      rad=dsqrt(1.d0-ecc2) 
      beta=1.d0/(1.d0+rad) 
      chk=eq(2)*eq(3)*beta 
      ch=1.d0-eq(2)**2*beta 
      ck=1.d0-eq(3)**2*beta 
      tgim2=eq(4)**2+eq(5)**2 
!     tgim=dsqrt(tgim2)                                                 
      opwz=1.d0+tgim2 
!   mean motion                                                         
      enne=dsqrt(gm/eq(1)**3) 
!  unit vectors of the equinoctal reference system (Broucke and         
!  Cefola 1972, CM 5, 303--310) are f, g and the angular momentum       
!  unit vector                                                          
      f(1)=(1.d0-eq(4)**2+eq(5)**2)/opwz 
      f(2)=2*eq(4)*eq(5)/opwz 
      f(3)=-2*eq(4)/opwz 
      g(1)=2*eq(4)*eq(5)/opwz 
      g(2)=(1.d0+eq(4)**2-eq(5)**2)/opwz 
      g(3)=2*eq(5)/opwz 
!  Non singular Kepler equation                                         
      ecc=dsqrt(ecc2) 
      if(ecc.lt.eps)then 
!  for negligible eccentricity, the eccentric longitude is              
!  set equal to the mean longitude                                      
         fe=eq(6) 
         cosf=cos(fe) 
         sinf=sin(fe) 
      else 
!  mean longitude is reduced to the interval (dig, dig+2*pig)           
!  with dig=longitude of pericentre                                     
         tlong=eq(6) 
         dig=atan2(eq(2),eq(3)) 
         tlong=princ(tlong-dig)+dig 
!  initial condition for Newton's method is always                      
!  the longitude of apocentre; this ensures the right                   
!  convexity for secure convergence                                     
         fe=dig+pig 
         do 16 i=1,100 
           cosf=cos(fe) 
           sinf=sin(fe) 
           df=(fe-tlong+eq(2)*cosf-eq(3)*sinf)/                         &
     &        (1.d0-eq(2)*sinf-eq(3)*cosf)                              
           if(dabs(df).lt.eps)goto 17 
   16      fe=fe-df 
!  convergence problems -- this should happen only for                  
!  extremely high eccentricity                                          
         WRITE(*,*)' eqcar : df=',df, ' elements ',eq 
!         stop '****** equcar: 100 iterations of Newton ******'          
      endif 
!  cartesian coordinates in the equinoctal frame x2, y2, 0              
   17 x2=eq(1)*(ch*cosf+chk*sinf-eq(3)) 
      y2=eq(1)*(ck*sinf+chk*cosf-eq(2)) 
      do 18 i=1,3 
   18   x(i)=x2*f(i)+y2*g(i) 
!  cartesian velocities in the equinoctal frame xp2, yp2, 0             
      de=enne*eq(1)**2/dsqrt(x2**2+y2**2) 
      xp2=de*(chk*cosf-ch*sinf) 
      yp2=de*(ck*cosf-chk*sinf) 
      do 19 i=1,3 
   19   x(i+3)=xp2*f(i)+yp2*g(i) 
      return 
      END subroutine equcar                                          
!=======================================================                
!   {\bf kepequ}: coordinate change from keplerian to                   
!   equinoctal elements                                                 
!     eq: equinoctal elements, see carequ                               
!                                                                       
!     el(1)=a                                                           
!                                                                       
!     el(2)=ecc                                                         
!                                                                       
!     el(3)=inclination (radians)                                       
!                                                                       
!     el(4)=longitude asc. node (radians)                               
!                                                                       
!     el(5)=argument of pericentre (radians)                            
!                                                                       
!     el(6)=mean anomaly (radians)                                      
!                                                                       
!   bugs: no test for hyperbolic/parabolic case; element a              
!   is copied whatever its meaning                                      
!=======================================================                
      subroutine kepequ(el,eq) 
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) :: el(6)
      DOUBLE PRECISION, INTENT(OUT) :: eq(6)
      DOUBLE PRECISION dig,ecc,tgim,princ
      eq(1)=el(1) 
      dig=el(4)+el(5) 
      ecc=el(2) 
      eq(2)=ecc*sin(dig) 
      eq(3)=ecc*cos(dig) 
      tgim=tan(el(3)/2.d0) 
      eq(4)=tgim*dsin(el(4)) 
      eq(5)=tgim*dcos(el(4)) 
      eq(6)=dig+el(6) 
      eq(6)=princ(eq(6)) 
      return 
      END subroutine kepequ                                          
!=======================================================                
!   {\bf equkep}: coordinate change from equinoctal elements            
!   to keplerian                                                        
!                                                                       
!     eq: equinoctal elements, see equcar/carequ                        
!                                                                       
!     el: keplerian elements, see kepequ                                
!                                                                       
!     eps=control on eccentricity/inclination; for smaller values,      
!     the angles eq(4), eq(5) are set to arbitrary value 0              
!                                                                       
!   bugs: no test for hyperbolic/parabolic case; element a              
!   is copied whatever its meaning                                      
!=======================================================                
      subroutine equkep(eq,eps,el)
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) :: eq(6),eps
      DOUBLE PRECISION, INTENT(OUT) :: el(6)
      DOUBLE PRECISION dig,ecc,tgi2,princ
! 
      el(1)=eq(1) 
!  test on eccentricity                                                 
      ecc=sqrt(eq(2)**2+eq(3)**2) 
      if(ecc.lt.eps)then 
         dig=0.d0 
      else 
         dig=atan2(eq(2),eq(3)) 
      endif 
      el(2)=ecc 
!   test on tangent of half inclination                                 
      tgi2=sqrt(eq(4)**2+eq(5)**2) 
      if(tgi2.lt.eps)then 
         el(4)=0.d0 
      else 
         el(4)=atan2(eq(4),eq(5)) 
      endif 
      el(3)=2.d0*atan(tgi2) 
!   angular variables                                                   
      el(5)=dig-el(4) 
      el(6)=eq(6)-dig 
      el(4)=princ(el(4)) 
      el(5)=princ(el(5)) 
      el(6)=princ(el(6)) 
      return 
      END subroutine equkep                                          
!=======================================================                
!   {\bf eqpequ}: coordinate change from equinoctal                     
!   polar to equinoctal elements                                        
!     eq: equinoctal elements, see carequ                               
!                                                                       
!     el(1)=a                                                           
!                                                                       
!     el(2)=e sin (omega)                                               
!                                                                       
!     el(3)=e cos (omega)                                               
!                                                                       
!     el(4)=tg (I/2) sin (Omega)                                        
!                                                                       
!     el(5)=tg (I/2) cos (Omega)                                        
!                                                                       
!     el(6)=mean argument of latitude (=mean anom+ omega)               
!                                                                       
!       omega= argument of pericentre                                   
!       Omega= longitude of asc. node                                   
!                                                                       
!   limitations: no test for hyperbolic/parabolic case; element a       
!   is copied whatever its meaning                                      
!=======================================================                
      subroutine eqpequ(el,eq)
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) :: el(6)
      DOUBLE PRECISION, INTENT(OUT) :: eq(6)
      DOUBLE PRECISION omnod,co,so,princ
!
      eq(1)=el(1) 
      omnod=atan2(el(4),el(5)) 
      co=cos(omnod) 
      so=sin(omnod) 
      eq(2)=el(2)*so+el(3)*co 
      eq(3)=el(3)*co-el(2)*so 
      eq(4)=el(4) 
      eq(5)=el(5) 
      eq(6)=princ(el(6)+omnod) 
      return 
      END subroutine eqpequ                                          
!=======================================================                
!   {\bf equeqp}: coordinate change from equinoctal elements            
!   to equinoctal polar                                                 
!                                                                       
!     eq: equinoctal elements, see equcar/carequ                        
!                                                                       
!     el: equinoctal polar elements, see eqpequ                         
!                                                                       
!     eps=control on inclination; for smaller values,                   
!       the computation stops                                           
!   limitation: no test for hyperbolic/parabolic case; element a        
!   is copied whatever its meaning                                      
!=======================================================                
      subroutine equeqp(eq,eps,el) 
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) :: eq(6),eps
      DOUBLE PRECISION, INTENT(OUT) :: el(6)
      DOUBLE PRECISION tgi2,omnod,co,so,princ
!
      el(1)=eq(1) 
!  test on inclination                                                  
      tgi2=sqrt(eq(4)**2+eq(5)**2) 
      if(tgi2.lt.eps)then 
         write(*,*)' inclination zero cannot be handled in eqp' 
         stop 
      else 
         omnod=atan2(eq(4),eq(5)) 
      endif 
      co=cos(omnod) 
      so=sin(omnod) 
      el(2)=eq(2)*co-eq(3)*so 
      el(3)=eq(3)*co+eq(2)*so 
      el(4)=eq(4) 
      el(5)=eq(5) 
      el(6)=princ(eq(6)-omnod) 
      return 
      END subroutine equeqp                                          
! Copyright (C) 1997 by Mario Carpino (carpino@brera.mi.astro.it)       
! Version: February 24, 1997                                            
!                                                                       
! Transformation between keplerian and non-singular (equinoctal) element
! with jacobian matrix                                                  
! (all angles in radians)                                               
!                                                                       
! INPUT:    elem      -  keplerian elements                             
!                            elem(1) = semimajor axis                   
!                            elem(2) = eccentricity                     
!                            elem(3) = inclination                      
!                            elem(4) = longitude of ascending node      
!                            elem(5) = argument of pericenter           
!                            elem(6) = mean anomaly                     
!                                                                       
! OUTPUT:   ekns      -  non-singular elements                          
!                            ekns(1) = semimajor axis                   
!                            ekns(2) = h = e*sin(varpi)                 
!                            ekns(3) = k = e*cos(varpi)                 
!                            ekns(4) = P = tan(i/2)*sin(Omega)          
!                            ekns(5) = Q = tan(i/2)*cos(Omega)          
!                            ekns(6) = mean longitude                   
!           dkdns     -  jacobian matrix d(ekns)/d(elem)                
!                                                                       
      subroutine ekensd(elem,ekns,dkdns) 
      USE fund_const
      implicit none 
                                                                        
      double precision elem(6),ekns(6),dkdns(6,6) 
      double precision alp,ti2,sinalp,cosalp,dti2,sinom,cosom 
      integer i,k 
                                                                        
      ekns(1)=elem(1) 
! Longitude of pericenter                                               
      alp=elem(5)+elem(4) 
      sinalp=sin(alp) 
      cosalp=cos(alp) 
      ekns(2)=elem(2)*sinalp 
      ekns(3)=elem(2)*cosalp 
      ti2=tan(elem(3)/2) 
      sinom=sin(elem(4)) 
      cosom=cos(elem(4)) 
      ekns(4)=ti2*sinom 
      ekns(5)=ti2*cosom 
      ekns(6)=alp+elem(6) 
      ekns(6)=mod(ekns(6),dpig) 
                                                                        
      do 1 i=1,6 
      do 1 k=1,6 
    1 dkdns(i,k)=0.d0 
                                                                        
      dkdns(1,1)=1.d0 
      dkdns(2,2)=sinalp 
      dkdns(2,5)=ekns(3) 
      dkdns(2,4)=ekns(3) 
      dkdns(3,2)=cosalp 
      dkdns(3,5)=-ekns(2) 
      dkdns(3,4)=-ekns(2) 
      dti2=1.d0/(2*cos(elem(3)/2)**2) 
      dkdns(4,3)=dti2*sinom 
      dkdns(4,4)=ekns(5) 
      dkdns(5,3)=dti2*cosom 
      dkdns(5,4)=-ekns(4) 
      dkdns(6,4)=1 
      dkdns(6,5)=1 
      dkdns(6,6)=1 
                                                                        
   END SUBROUTINE ekensd
!  ==================================================================== 
!  KEPCAR - fast coordinate change                                      
!  ==================================================================== 
      SUBROUTINE kepcar(ek,gm,ivel,x) 
!  ====================================================================
      USE fund_const 
      IMPLICIT NONE 
! INPUT                                                                 
! keplerian elements, mass of central body                              
      DOUBLE PRECISION ek(6),gm 
! ivel=1 also velocities, ivel=0 no                                     
      INTEGER ivel 
! OUTPUT                                                                
      DOUBLE PRECISION x(6) 
! END INTERFACE                                                         
      DOUBLE PRECISION f(3),g(3),ecc,h,k,rad,beta,dig,chk,ch,ck,p,q,ecc2 
      DOUBLE PRECISION tgim2,tlong,fe,sinf,cosf,princ,x2,y2,xp2,yp2,de 
      DOUBLE PRECISION opwz,enne,df 
! controls: convergence of kepler's equation                            
      DOUBLE PRECISION eps 
      INTEGER i 
! ======================================================================
!   test for hyperbolic orbit                                           
      if(ek(1).le.0.d0)then 
         WRITE(*,*)'****** kepcar: hyperbolic/parabolic orbit ******' 
      endif 
!  non--singular intermediate variables                                 
      ecc=ek(2) 
      ecc2=ecc**2 
      rad=dsqrt(1.d0-ecc2) 
      beta=1.d0/(1.d0+rad) 
      dig=ek(4)+ek(5) 
      k=ecc*cos(dig) 
      h=ecc*sin(dig) 
      chk=h*k*beta 
      ch=1.d0-h**2*beta 
      ck=1.d0-k**2*beta 
      tgim2=tan(ek(3)/2.d0) 
!     tgim=dsqrt(tgim2)                                                 
      opwz=1.d0+tgim2 
      p=tgim2*sin(ek(4)) 
      q=tgim2*cos(ek(4)) 
!   mean motion                                                         
      enne=dsqrt(gm/ek(1)**3) 
!  unit vectors of the equinoctal reference system (Broucke and         
!  Cefola 1972, CM 5, 303--310) are f, g and the angular momentum       
!  unit vector                                                          
      f(1)=(1.d0-p**2+q**2)/opwz 
      f(2)=2*p*q/opwz 
      f(3)=-2*p/opwz 
      g(1)=2*p*q/opwz 
      g(2)=(1.d0+p**2-q**2)/opwz 
      g(3)=2*q/opwz 
!  Non singular Kepler equation  
      eps=100*epsilon(1.d0)                                       
      IF(ecc.lt.eps)THEN 
!  for negligible eccentricity, the eccentric longitude is              
!  set equal to the mean longitude                                      
         fe=princ(dig+ek(6)) 
         cosf=cos(fe) 
         sinf=sin(fe) 
      ELSE 
!  mean longitude is reduced to the interval (dig, dig+2*pig)           
!  with dig=longitude of pericentre                                     
         tlong=dig+ek(6) 
         tlong=princ(tlong-dig)+dig 
!  initial condition for Newton's method is always                      
!  the longitude of apocentre; this ensures the right                   
!  convexity for secure convergence                                     
         fe=dig+pig 
         DO i=1,10 
            cosf=cos(fe) 
            sinf=sin(fe) 
            df=(fe-tlong+h*cosf-k*sinf)/                                &
     &           (1.d0-h*sinf-k*cosf)                                   
            fe=fe-df 
            IF(dabs(df).lt.eps)GOTO 17 
         ENDDO 
!  convergence problems -- this should happen only for                  
!  extremely high eccentricity                                          
         WRITE(*,*)'****** kepcar: 10 iterations of Newton ******' 
      ENDIF 
!  cartesian coordinates in the equinoctal frame x2, y2, 0              
   17 x2=ek(1)*(ch*cosf+chk*sinf-k) 
      y2=ek(1)*(ck*sinf+chk*cosf-h) 
      DO i=1,3 
        x(i)=x2*f(i)+y2*g(i) 
      ENDDO 
      IF(ivel.eq.0)RETURN 
!  cartesian velocities in the equinoctal frame xp2, yp2, 0             
      de=enne*ek(1)**2/dsqrt(x2**2+y2**2) 
      xp2=de*(chk*cosf-ch*sinf) 
      yp2=de*(ck*cosf-chk*sinf) 
      DO  i=1,3 
        x(i+3)=xp2*f(i)+yp2*g(i) 
      ENDDO 
      return 
    END SUBROUTINE kepcar
                                                             

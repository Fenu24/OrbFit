! ==========MODULE yark_pert============
! CONTAINS
!                 yarkdi   yarkovsky diurnal                            
!                 yarkse   yarkovsky seasonal                           
!               yarkinit   initialisation of yarkovsky 
! MODULES AND HEADERS
! yark_pert.o: \
!	../include/sysdep.h90 \
!	../suit/FUND_CONST.mod 

MODULE yark_pert
USE fund_const
USE dyn_param
IMPLICIT NONE
PRIVATE

DOUBLE PRECISION :: yark_exp

PUBLIC yarkdi, yarkse, yarkinit, secular_nongrav,yark_exp, sec_nong9

! private common data
! former yarkom.h
! controls and directory for Yarkovski force
      INTEGER iyarpt
      CHARACTER*80 yardir

! fromer yarkov.h 
! common containing all physical information
! on the current asteroid needed to compute Yarkovsky acceleration
      DOUBLE PRECISION yarkp(10),beya(7),alya(7)
      DOUBLE PRECISION spya,sqya,etaya75,fmeaya,thfacya,radfluya
! logical flggs: availability of physical data, 
!  has yarkinit routine been called
      LOGICAL yarfil,yarini
! secular perturbation in semiajor axis 
!      DOUBLE PRECISION dadt
      PUBLIC iyarpt,yardir,yarfil,yarini !,dadt
      
CONTAINS 
! =========================================================
! SEC_NONG9   added by AM 9/9/2009
! secular perturbation on semimajor axis, presumably due to
! non gravitational perturbations (including Yarkovsky)
! implemented as acceleration along the transverse direction
! 22/3/2018 added radial component 
SUBROUTINE sec_nong9(xb,vb,s,sv,secacc,iyark,dadt)
!  interface: INPUT
  DOUBLE PRECISION, INTENT(IN), DIMENSION(3) :: xb,vb ! position and velocity, barycentric
  DOUBLE PRECISION, INTENT(IN), DIMENSION(3) :: s,sv !  position and velocity of the sun, barycentric
  INTEGER,INTENT(IN) :: iyark ! 3= A2 only; 4=A2 and A1
  DOUBLE PRECISION, INTENT(IN) :: dadt ! in au/My
!  interface: OUTPUT
  DOUBLE PRECISION, INTENT(OUT), DIMENSION(3) :: secacc
! end interface
!        DOUBLE PRECISION :: vvec(3),vv2, rr, vsize, prscal, factor, conv
  DOUBLE PRECISION :: rr,rr2,vv2,rv
  DOUBLE PRECISION :: factor,conv,yark_acc,trans_size, angm2, gmsy, sqarg
  DOUBLE PRECISION :: vsize, prscal
  DOUBLE PRECISION, DIMENSION(3) :: x, v, trans, angm
! =========================================================
! conversion to heliocentric coordinates
  x = xb-s
  v = vb-sv
! Steve Chesley's implementation, assumes heliocentric position and velocity
  rr = vsize(x)
  rr2=prscal(x,x)
  vv2=prscal(v,v)
  rv =prscal(x,v)
  CALL prvec(x,v,angm)
  angm2=prscal(angm,angm)
  gmsy=gms*365.25d0**2 ! here units are AU, y
  sqarg=2*gmsy/rr-vv2
  IF(sqarg.lt.0.d0)sqarg=0
  factor=angm2*dsqrt(sqarg)/gmsy !fails for heliocentric hyperbolic
! warning: dadt is in au/My, thus it has to be converted in au/y
  conv=1.d6 
  yark_acc=0.5*factor/conv/rr2
  trans=v-(x*rv/rr2)
  trans_size=vsize(trans)
  secacc=dadt*yark_acc*trans/trans_size
  IF(iyark.eq.4)THEN
     secacc=secacc+(30.d0*yark_acc*abs(dadt)/rr)*x
  ENDIF
END SUBROUTINE sec_nong9
! =========================================================
! SECULAR_NONGRAV   added 17/9/2008
! secular perturbation on semimajor axis, presumably due to
! non gravitational perturbations (including Yarkovsky)
! implemented as acceleration along the velocity
! 22.11.2016 - Added derivatives, LDimare
      SUBROUTINE secular_nongrav(x,v,secacc,der)
!  interface: INPUT
        DOUBLE PRECISION, INTENT(IN), DIMENSION(3) :: x,v ! position and velocity, heliocentric
!  interface: OUTPUT
        DOUBLE PRECISION, INTENT(OUT), DIMENSION(3) :: secacc
        DOUBLE PRECISION, INTENT(OUT), OPTIONAL :: der(3,6) ! derivatives of acc wrt pos and vel
! end interface
        DOUBLE PRECISION :: rr,rr2,vv2,rv,angm_x,angm_y,angm_z,angm2
        DOUBLE PRECISION :: factor,conv,yark_acc,trans_size,trans(3)
        DOUBLE PRECISION :: vsize, prscal
        DOUBLE PRECISION,DIMENSION(3) :: ang_mom
!!! 1/r^3 model
        DOUBLE PRECISION :: p,e,a,n,a2,integral
! derivatives
        DOUBLE PRECISION, DIMENSION(3,6) :: dtrans,dtvers
        DOUBLE PRECISION, DIMENSION(3,3) :: idmat ! id. matrix
        INTEGER          :: i,j
! =========================================================
! secacc=A2/r^yark_exp
        rr=vsize(x)
        trans=v-(x*DOT_PRODUCT(x,v)/rr**2)
        trans_size=vsize(trans)
!        secacc=1.d-15/rr**yark_exp*trans/trans_size
! Good for Apophis!
        secacc=1.d-10/rr**yark_exp*trans/trans_size
        IF(PRESENT(der))THEN
           ! derivatives of vector trans
           dtrans=0.d0
           CALL eye(3,idmat)
           dtrans(1:3,1:3)=-DOT_PRODUCT(x,v)/rr**2*idmat
           dtrans(1:3,4:6)=idmat
           DO i=1,3
              DO j=1,3
                 dtrans(i,j)=dtrans(i,j)-x(i)*v(j)/rr**2+2.d0*DOT_PRODUCT(x,v)/rr**4*x(i)*x(j)
              ENDDO
              DO j=4,6
                 dtrans(i,j)=dtrans(i,j)-x(j-3)*x(i)/rr**2
              ENDDO
           ENDDO
           ! derivatives of trans/trans_size
           dtvers=1.d0/trans_size*dtrans
           DO i=1,3
              DO j=1,6
                 dtvers(i,j)=dtvers(i,j)-1.d0/trans_size**3*DOT_PRODUCT(trans,dtrans(1:3,j))*trans(i)
              ENDDO
           ENDDO
           der=0.d0
           DO i=1,3    ! index of acceleration component to be derived
              DO j=1,3 ! index of position component wrt derive
                 ! first piece: derivative of 1/rr**(yark_exp)
                 der(i,j)=-yark_exp/rr**(yark_exp+2)*x(j)*trans(i)/trans_size
              ENDDO
              DO j=1,6 
                 ! second piece: derivative of trans/trans_size
                 der(i,j)=der(i,j)+dtvers(i,j)/rr**yark_exp
              ENDDO
           ENDDO
           der=der*1.d-10
        ENDIF

        RETURN
!! below we have obsolete formulations to directly compute
!! dadt
! useful quantities
        rr=vsize(x)
        rr2=rr**2
        CALL prvec(x,v,ang_mom)
        a=-gms/2/(DOT_PRODUCT(v,v)/2-gms/rr)
        n=SQRT(gms/a**3)
        p=DOT_PRODUCT(ang_mom,ang_mom)/gms
        e=SQRT(ABS(1-p/a))
        trans=v-(x*DOT_PRODUCT(x,v)/rr2)
        trans_size=vsize(trans)
! warning: dadt is in AU/My, thus it has to be converted in AU/d
        conv=1.d6*365.25d0
! 1/r^2 model
        secacc=0.5d0/conv*n*p*a/rr2*trans/trans_size
!        return
! 1/r^3 model
        secacc=1.d0/conv*n*p**2*a/(2+e**2)/rr**3*trans/trans_size
      END SUBROUTINE secular_nongrav
! ******************************************************************    
      SUBROUTINE yarkdi(xast,a,iparti) 
! ******************************************************************    
!                                                                       
! This subroutine computes the heliocentric components of the           
! Yarkovsky thermal acceleration -- the diurnal variant only.           
! If the flag (iparti) in the common block is set to 1, one gets at     
! the output also partials wrt to some parameters of the thermal        
! model (if these are desired to be adjusted).                          
!                                                                       
! Input parameters:                                                     
! -----------------                                                     
!                                                                       
! - via header   xast(3) ... heliocentric coordinates of the body (in AU
! - via common   yarkp(1-3) ... sx, sy, sz (unit vector of the          
!                               body's spin axis orientation)           
!                yarkp(4-5) ... k_0 and k_1 parameters of the           
!                               surface thermal conductivity            
!                               [K(T) = k_0 + k_1 T_av^3]               
!                yarkp(6) ... density of the surface layer              
!                yarkp(7) ... radius of the body                        
!                yarkp(8) ... rotation frequency                        
!                yarkp(9) ... surface absorptivity   
!                yarkp(10)... bulk density   
!                iparti   ... partials (yes=1/no=0)                     
!                             [presently only partials listed below,    
!                              a(4) - a(21) are available]              
!                                                                       
!                                                                       
! Output parameters: a(1-3) ... diurnal acceleration                    
! ------------------ a(4-6) ... partials wrt the radius of the body     
!                    a(7-9) ... partials wrt the thermal conductivity   
!                               parameter k_0                           
!                    a(10-12) ... partials wrt the thermal conductivity 
!                                 parameter k_1                         
!                    a(13-15) ... partials wrt the x-component of the   
!                                 spin axis unit vector                 
!                    a(16-18) ... partials wrt the y-component of the   
!                                 spin axis unit vector                 
!                    a(19-21) ... partials wrt the z-component of the   
!                                 spin axis unit vector                 
!                                                                       
! SI units are assumed internally in the subroutine, but the results    
! (e.g. accelerations) are given in AU and days.                        
!                                                                       
! Written by: D. Vokrouhlicky, Oct 99                                   
! (queries to vokrouhl@mbox.cesnet.cz)                                  
! ..................................................................    
      implicit double precision (a-h,o-z) 
! here we specify two more parameters that eventually might be changed: 
! -- the average bulk density of the body (densityb) which is now set   
!    to 2 g/cm^3                                                        
! -- the heat capacity of the surface layer (capacity) which is now set 
!    to 680 J/kg/K                                                      
      DOUBLE PRECISION, parameter :: capacity=680.d0,solcon=1371.d0
      DOUBLE PRECISION, parameter :: emiss=0.9d0,stefboltz=5.66962d-8,clight3=8.99377374d8 
      DOUBLE PRECISION, parameter :: dsqrt2=1.414213562373d0,dsqrt23=1414.213562373d0 
      DOUBLE PRECISION, parameter :: aceuni=0.049900176d0
      DOUBLE PRECISION :: densityb 
! input: asteroid position, flag for partials                           
      double precision xast(3) 
      integer iparti 
! output: acceleration and partials                                     
      double precision a(21) 
! internal variables                                                    
      double precision vprod1(3),vprod2(3)
      double precision rau2,rau,xn,yn,zn,radfku,tstar, tav1000,surcon,bgama,theta,diudepth,rp 
! physical data on the current asteroid                                 
! ----------------------------------------------------------------------
      densityb=yarkp(10)
      rau2=xast(1)*xast(1)+xast(2)*xast(2)+xast(3)*xast(3) 
      rau=dsqrt(rau2) 
      xn=xast(1)/rau 
      yn=xast(2)/rau 
      zn=xast(3)/rau 
! initializations & constants                                           
      radflu=solcon/rau2 
! - subsolar temperature                                                
      tstar=(yarkp(9)*radflu/emiss/stefboltz)**0.25d0 
      tav1000=tstar/dsqrt23 
! - surface conductivity                                                
      surcon=yarkp(4)+yarkp(5)*(tav1000**3) 
! - thermal inertia & diurnal thermal parameter                         
      bgama=dsqrt(surcon*yarkp(6)*capacity) 
      theta=bgama*dsqrt(yarkp(8))/emiss/stefboltz/(tstar**3) 
      diudepth=dsqrt(surcon/yarkp(6)/capacity/yarkp(8)) 
! - radius of the body scaled by the depth of the diurnal wave          
      rp=yarkp(7)/diudepth 
      al=dsqrt2*rp 
      tau=theta/al 
      tau1=1.d0+tau 
! - the auxiliary functions A-D, a,b                                    
      cal=dcos(al) 
      sal=dsin(al) 
      if (al.lt.90.d0) then 
       ealm=dexp(-al) 
      else 
       ealm=0.d0 
      endif 
      af=3.d0*(al+2.d0)*ealm+(3.d0*(al-2.d0)*cal+al*(al-3.d0)*sal) 
      bf=al*(al+3.d0)*ealm+(-al*(al-3.d0)*cal+3.d0*(al-2.d0)*sal) 
      caf=-(al+2.d0)*ealm+(-(al-2.d0)*cal+al*sal) 
      cbf=-al*ealm-(al*cal+(al-2.d0)*sal) 
      ccf=caf+tau*af/tau1 
      cdf=cbf+tau*bf/tau1 
! - G exp(i delta) & amplitude computed                                 
      facp=aceuni*yarkp(9)*radflu/yarkp(7)/densityb/clight3 
      deno=ccf*ccf+cdf*cdf 
      deno1=deno*tau1 
      gcosd=(caf*ccf+cbf*cdf)/deno1 
      gsind=(cbf*ccf-caf*cdf)/deno1 
! geometric products                                                    
! - r x s                                                               
      vprod1(1)=yn*yarkp(3)-zn*yarkp(2) 
      vprod1(2)=zn*yarkp(1)-xn*yarkp(3) 
      vprod1(3)=xn*yarkp(2)-yn*yarkp(1) 
! - s x (r x s) = r - (r.s) s                                           
      scalar=xn*yarkp(1)+yn*yarkp(2)+zn*yarkp(3) 
      vprod2(1)=xn-scalar*yarkp(1) 
      vprod2(2)=yn-scalar*yarkp(2) 
      vprod2(3)=zn-scalar*yarkp(3) 
! diurnal acceleration                                                  
      a(1)=facp*(gsind*vprod1(1)+gcosd*vprod2(1)) 
      a(2)=facp*(gsind*vprod1(2)+gcosd*vprod2(2)) 
      a(3)=facp*(gsind*vprod1(3)+gcosd*vprod2(3)) 
! Partials?                                                             
      if (iparti.eq.0) return 
! - general                                                             
      cafp=-ealm+cal+(2.d0*al-1.d0)*sal 
      cbfp=-ealm-(2.d0*al-1.d0)*cal+sal 
      afp=3.d0*ealm+(al*al-3.d0)*cal+(al*(al-4.d0)+3.d0)*sal 
      bfp=(2.d0*al+3.d0)*ealm-(al*(al-4.d0)+3.d0)*cal                   &
     &     +(al*al-3.d0)*sal                                            
! - thermal conductivity parameters (k_0,k_1)                           
      xi1r=caf*ccf-cbf*cdf 
      xi1i=cbf*ccf+caf*cdf 
      xi2r=cafp*af-cbfp*bf 
      xi2i=cbfp*af+cafp*bf 
      xi2r=xi2r-caf*afp+cbf*bfp 
      xi2i=xi2i-cbf*afp-caf*bfp 
      deno=xi1r*xi1r+xi1i*xi1i 
      facr=1.d0+0.5d0*al*(xi2r*xi1r+xi2i*xi1i)/deno 
      faci=     0.5d0*al*(xi2i*xi1r-xi2r*xi1i)/deno 
      derikr=-tau*(gcosd*facr-gsind*faci)/tau1 
      deriki=-tau*(gsind*facr+gcosd*faci)/tau1 
      a(7)=facp*(deriki*vprod1(1)+derikr*vprod2(1)) 
      a(8)=facp*(deriki*vprod1(2)+derikr*vprod2(2)) 
      a(9)=facp*(deriki*vprod1(3)+derikr*vprod2(3)) 
      a(10)=a(7)*(tav1000**3) 
      a(11)=a(8)*(tav1000**3) 
      a(12)=a(9)*(tav1000**3) 
! - radius of the body                                                  
      rfac=(tau+tau1)/tau1 
      a(4)=-a(1)*rfac-2.d0*a(7) 
      a(5)=-a(2)*rfac-2.d0*a(8) 
      a(6)=-a(3)*rfac-2.d0*a(9) 
! - partials d_K (a), d_R (a) ...                                       
      a(4)=a(4)/yarkp(7) 
      a(5)=a(5)/yarkp(7) 
      a(6)=a(6)/yarkp(7) 
      a(7)=a(7)/surcon 
      a(8)=a(8)/surcon 
      a(9)=a(9)/surcon 
      a(10)=a(10)/surcon 
      a(11)=a(11)/surcon 
      a(12)=a(12)/surcon 
! - spin axis components                                                
! ... sx                                                                
      a(13)=-facp*gcosd*(xn*yarkp(1)+scalar) 
      a(14)=facp*(gsind*zn-gcosd*xn*yarkp(2)) 
      a(15)=-facp*(gsind*yn+gcosd*xn*yarkp(3)) 
! ... sy                                                                
      a(16)=-facp*(gsind*zn+gcosd*yn*yarkp(1)) 
      a(17)=-facp*gcosd*(yn*yarkp(2)+scalar) 
      a(18)=facp*(gsind*xn-gcosd*yn*yarkp(3)) 
! ... sz                                                                
      a(19)=facp*(gsind*yn-gcosd*zn*yarkp(1)) 
      a(20)=-facp*(gsind*xn+gcosd*zn*yarkp(2)) 
      a(21)=-facp*gcosd*(zn*yarkp(3)+scalar) 
      return 
      END SUBROUTINE yarkdi                                          
! ******************************************************************    
      SUBROUTINE yarkse(xast,vast,a,iparti) 
! ******************************************************************    
!                                                                       
! This subroutine computes the heliocentric components of the           
! Yarkovsky thermal acceleration -- the seasonal variant only.          
! If the flag (iparti) is set to 1, one gets at the output also         
! partials wrt to some parameters of the thermal model (if these        
! are desired to be adjusted).                                          
!                                                                       
! Input parameters:                                                     
! -----------------                                                     
!                                                                       
! - via header   iparti  ... partials (yes=1/no=0)                      
!                (xast,vast) ... state vector of the asteroid           
! - via common   yarkp(1-3) ... sx, sy, sz (unit vector of the          
!                               body's spin axis orientation)           
!                yarkp(4-5) ... k_0 and k_1 parameters of the           
!                               surface thermal conductivity            
!                               [K(T) = k_0 + k_1 T_av^3]               
!                yarkp(6) ... density of the surface layer              
!                yarkp(7) ... radius of the body                        
!                yarkp(8) ... rotation frequency                        
!                yarkp(9) ...  surface absorptivity  
!                yarkp(10)... bulk density   
!                + some more precomputed useful variables               
!                                                                       
! Output parameters: a(1-3) ... seasonal acceleration                   
! ------------------ a(4-6) ... partials wrt the radius of the body     
!                    a(7-9) ... partials wrt the thermal conductivity   
!                                                                       
! REM. PARTIALS ARE DISABLED AT THIS MOMENT                             
!                                                                       
! SI units are assumed throughout the subroutine, but the results       
! (e.g. accelerations) are given in AU and days.                        
!                                                                       
! Written by: D. Vokrouhlicky, Oct 99                                   
! (queries to vokrouhl@mbox.cesnet.cz)                                  
! ..................................................................    
      implicit double precision (a-h,o-z) 
! modified 17/4/2007, to comply with the Golevka code at JPL.
     integer, parameter:: napprox=7
!     integer, parameter:: napprox=12
      parameter (capacity=680.d0,dsqrt2=1.414213562373d0) 
      parameter (emiss=0.9d0,clight3=8.99377374d8,aceuni=0.049900176d0) 
      dimension xast(3),vast(3) 
      dimension brac(napprox),bras(napprox),gcosd(napprox),gsind(napprox),a(21)
      double precision densityb 
      integer iparti
      integer k
! ----------------------------------------------------------------------
       densityb=yarkp(10)
! - thermal inertia & seasonal thermal parameter                        
       bgama=dsqrt(yarkp(4)*yarkp(6)*capacity) 
       theta=bgama*thfacya/emiss 
       seadepth=dsqrt(yarkp(4)/yarkp(6)/capacity/fmeaya) 
! - radius of the body scaled by the depth of the seasonal wave         
       rp=yarkp(7)/seadepth 
       rp2=dsqrt2*rp 
       tau=theta*etaya75/rp2 
       tau1=1.d0+tau 
! - amplitude of the effect                                             
       fac=aceuni*yarkp(9)*radfluya/yarkp(7)/densityb/clight3/tau1 
! - G_k cos(d_k) & G_K sin(d_k) functions computed                      
       do 10 k=1,napprox 
        fk=k 
        alk=dsqrt(fk)*rp2 
! - the auxiliary functions A-D, a,b                                    
        cal=dcos(alk) 
        sal=dsin(alk) 
        if (alk.lt.90.d0) then 
         ealm=dexp(-alk) 
        else 
         ealm=0.d0 
        endif 
        af=3.d0*(alk+2.d0)*ealm+(3.d0*(alk-2.d0)*cal+alk*(alk-3.d0)*sal) 
        bf=alk*(alk+3.d0)*ealm+(-alk*(alk-3.d0)*cal+3.d0*(alk-2.d0)*sal) 
        caf=-(alk+2.d0)*ealm+(-(alk-2.d0)*cal+alk*sal) 
        cbf=-alk*ealm-(alk*cal+(alk-2.d0)*sal) 
        ccf=caf+tau*af/tau1 
        cdf=cbf+tau*bf/tau1 
! - G exp(i delta)                                                      
        deno=ccf*ccf+cdf*cdf 
        gcosd(k)=(caf*ccf+cbf*cdf)/deno 
        gsind(k)=(cbf*ccf-caf*cdf)/deno 
! compute cos- & sin-related brackets                                   
        brac(k)=spya*alya(k)*gcosd(k)+sqya*beya(k)*gsind(k) 
        bras(k)=sqya*beya(k)*gcosd(k)-spya*alya(k)*gsind(k) 
   10  continue 
! mean anomaly determined                                               
! - computed from the state vector                                      
      r2=xast(1)*xast(1)+xast(2)*xast(2)+xast(3)*xast(3) 
      v2=vast(1)*vast(1)+vast(2)*vast(2)+vast(3)*vast(3) 
      rdot=xast(1)*vast(1)+xast(2)*vast(2)+xast(3)*vast(3) 
      r=dsqrt(r2) 
      aaxi=1.d0/(2.d0/r-(v2/gms)) 
      esinu=rdot/dsqrt(aaxi*gms) 
      ecosu=(r*v2/gms)-1.d0 
      uano=datan2(esinu,ecosu) 
      anomaly=uano-esinu 
      if (anomaly.lt.0.d0) anomaly=anomaly+dpig 
! compute the sum...                                                    
      fact=0.d0 
      do 100 k=napprox,1,-1 
       fk=k 
       canomaly=dcos(fk*anomaly) 
       sanomaly=dsin(fk*anomaly) 
       fact=fact+(brac(k)*canomaly+bras(k)*sanomaly) 
  100 continue 
      fact=fact*fac 
! seasonal acceleration (~ factor * {\bf s})                            
      a(1)=fact*yarkp(1) 
      a(2)=fact*yarkp(2) 
      a(3)=fact*yarkp(3) 
! Partials? -- DISABLED AT THE MOMENT                                   
!      if (iparti.eq.0) return                                          
! - general                                                             
!      cafp=-ealm+cal+(2.d0*al-1.d0)*sal                                
!      cbfp=-ealm-(2.d0*al-1.d0)*cal+sal                                
!      afp=3.d0*ealm+(al*al-3.d0)*cal+(al*(al-4.d0)+3.d0)*sal           
!      bfp=(2.d0*al+3.d0)*ealm-(al*(al-4.d0)+3.d0)*cal                  
!     .     +(al*al-3.d0)*sal                                           
! - thermal conductivity parameters (k_0,k_1)                           
!      xi1r=caf*ccf-cbf*cdf                                             
!      xi1i=cbf*ccf+caf*cdf                                             
!      xi2r=cafp*af-cbfp*bf                                             
!      xi2i=cbfp*af+cafp*bf                                             
!      xi2r=xi2r-caf*afp+cbf*bfp                                        
!      xi2i=xi2i-cbf*afp-caf*bfp                                        
!      deno=xi1r*xi1r+xi1i*xi1i                                         
!      facr=1.d0+0.5d0*al*(xi2r*xi1r+xi2i*xi1i)/deno                    
!      faci=     0.5d0*al*(xi2i*xi1r-xi2r*xi1i)/deno                    
!      derikr=-tau*(gcosd*facr-gsind*faci)/tau1                         
!      deriki=-tau*(gsind*facr+gcosd*faci)/tau1                         
!      a(7)=fac*(deriki*vprod1(1)+derikr*vprod2(1))                     
!      a(8)=fac*(deriki*vprod1(2)+derikr*vprod2(2))                     
!      a(9)=fac*(deriki*vprod1(3)+derikr*vprod2(3))                     
!      a(10)=a(7)*(tav1000**3)                                          
!      a(11)=a(8)*(tav1000**3)                                          
!      a(12)=a(9)*(tav1000**3)                                          
! - radius of the body                                                  
!      rfac=(tau+tau1)/tau1                                             
!      a(4)=-a(1)*rfac-2.d0*a(7)                                        
!      a(5)=-a(2)*rfac-2.d0*a(8)                                        
!      a(6)=-a(3)*rfac-2.d0*a(9)                                        
! - partials d_K (a), d_R (a) ...                                       
!      a(4)=a(4)/yarkp(7)                                               
!      a(5)=a(5)/yarkp(7)                                               
!      a(6)=a(6)/yarkp(7)                                               
!      a(7)=a(7)/surcon!!!!!! --> yarkp(4)                              
!      a(8)=a(8)/surcon                                                 
!      a(9)=a(9)/surcon                                                 
!      a(10)=a(10)/surcon                                               
!      a(11)=a(11)/surcon                                               
!      a(12)=a(12)/surcon                                               
! - spin axis components                                                
      return 
      END SUBROUTINE yarkse                                          
! ==================================================================    
! yarkinit: initialisation of the yarkovsky force model 
!           for a given asteroid
!           written by A. Milani & D. Vokrouhlicky, Oct 99              
SUBROUTINE yarkinit(astnam,elem)
  USE orbit_elements 
  USE dyn_param
  IMPLICIT NONE 
! ============ BEGIN INTERFACE ===================
  CHARACTER*(*),    INTENT(IN) ::  astnam 
  TYPE(orbit_elem), INTENT(IN) :: elem
! ============ END INTERFACE =====================
  CHARACTER*80 :: file 
  DOUBLE PRECISION :: lat,long,emiss,stefboltz,argu,argu2,tstarya,eta 
  TYPE(orbit_elem) :: elekep
  DOUBLE PRECISION :: elkep(6),pvya(3),qvya(3),nvya(3),enne,cgam,obli 
  INTEGER          :: unit,le, fail_flag
  INCLUDE 'sysdep.h90' 
! yar is the logical flag for the existence of the physical data        
! allowing computation of Yarkovsky; otherwise, the non gravitational   
! force is set to zero                                                  
!                                                                       
  IF(iyark.EQ.0.OR.iyark.GE.3)RETURN 
  yarini=.true. 
! convert elements to keplerian                                         
  call coo_cha(elem,'KEP',elekep,fail_flag)
  IF(fail_flag.ge.4)THEN
     WRITE(*,*)' yarkinit: not possible with comet ', fail_flag, elekep
     STOP
  ELSE
     elkep=elekep%coord
  ENDIF 
! compute the name of the file which could contain the yarkovsky data   
  CALL filnam(yardir,astnam,'yar',file,le) 
  INQUIRE(file=file(1:le),exist=yarfil) 
  IF(yarfil)THEN 
     call filopn(unit,file(1:le),'old') 
     read(unit,*,end=111) 
! equatorial RA and DEC of the spin axis; beware of old names
     read(unit,*,end=111)long 
     read(unit,*,end=111)lat 
! - via common   yarkp(1-3) ... sx, sy, sz (unit vector of the          
!                               body's spin axis orientation)           
!                yarkp(4-5) ... k_0 and k_1 parameters of the           
!                               surface thermal conductivity            
!                               [K(T) = k_0 + k_1 T_av^3]               
!                yarkp(6) ... density of the surface layer              
!                yarkp(7) ... radius of the body                        
!                yarkp(8) ... rotation frequency                        
!                yarkp(9) ... surface absorptivity
!                yarkp(10)... bulk density                      
     yarkp(1)=dcos(lat*radeg)*dcos(long*radeg) 
     yarkp(2)=dcos(lat*radeg)*dsin(long*radeg) 
     yarkp(3)=dsin(lat*radeg) 
! rotate to ecliptic
     yarkp(1:3)=MATMUL(roteqec,yarkp(1:3))
     read(unit,*,end=111)yarkp(4) 
     read(unit,*,end=111)yarkp(5) 
     read(unit,*,end=111)yarkp(6) 
     read(unit,*,end=111)yarkp(7) 
     read(unit,*,end=111)yarkp(8) 
     read(unit,*,end=111)yarkp(9) 
     read(unit,*,end=111)yarkp(10) 
! precompute some variables for the seasonal variant of the Yarkovsky   
! effect:                                                               
! - constants                                                           
     emiss=0.9d0 
     stefboltz=5.66962d-8 
! - mean motion & solar radiation flux at r=a                           
     fmeaya=(1.9909837d-7)/elkep(1)/dsqrt(elkep(1)) 
     radfluya=1371.d0/elkep(1)/elkep(1) 
! - subsolar temperature                                                
     tstarya=(yarkp(9)*radfluya/emiss/stefboltz)**0.25d0 
     thfacya=dsqrt(fmeaya)/stefboltz/(tstarya**3) 
! - projections s_P and s_Q of the spin axis computed                   
     pvya(1)=dcos(elkep(4))*dcos(elkep(5))-                         &
     &           dcos(elkep(3))*dsin(elkep(4))*dsin(elkep(5))           
     pvya(2)=dsin(elkep(4))*dcos(elkep(5))+                         &
     &           dcos(elkep(3))*dcos(elkep(4))*dsin(elkep(5))           
     pvya(3)=dsin(elkep(3))*dsin(elkep(5)) 
     qvya(1)=-dcos(elkep(4))*dsin(elkep(5))-                        &
     &           dcos(elkep(3))*dsin(elkep(4))*dcos(elkep(5))           
     qvya(2)=-dsin(elkep(4))*dsin(elkep(5))+                        &
     &         dcos(elkep(3))*dcos(elkep(4))*dcos(elkep(5))             
     qvya(3)=dsin(elkep(3))*dcos(elkep(5)) 
     nvya(1)=dsin(elkep(3))*dsin(elkep(4)) 
     nvya(2)=-dsin(elkep(3))*dcos(elkep(4)) 
     nvya(3)=dcos(elkep(3)) 
     spya=yarkp(1)*pvya(1)+yarkp(2)*pvya(2)+yarkp(3)*pvya(3) 
     sqya=yarkp(1)*qvya(1)+yarkp(2)*qvya(2)+yarkp(3)*qvya(3) 
     cgam=yarkp(1)*nvya(1)+yarkp(2)*nvya(2)+yarkp(3)*nvya(3) 
     obli=dacos(cgam)/radeg 
     write(*,*)'Incl, Omega, omega (deg) ',elkep(3)/radeg,elkep(4)/radeg,elkep(5)/radeg
     write(*,*)' Obliquity of the spin axis; Yarkovsky: ',obli 
! - compute the \alpha(k) and \beta(k) coefficients                     
     eta=dsqrt(1.d0-elkep(2)*elkep(2)) 
     etaya75=eta**0.75d0 
! -- \beta_1(x) ... \beta_7(x) functions                                
     argu=elkep(2) 
     argu2=argu*argu 
     beya(1)=eta*(1.d0+argu2*(-1152.d0+argu2*(48.d0-argu2))/9216.d0) 
     argu=2.d0*elkep(2) 
     argu2=argu*argu 
     beya(2)=eta*argu*(1.d0+argu2*(-1920.d0+argu2*(60.d0-argu2))    &
     &           /23040.d0)                                             
     argu=3.d0*elkep(2) 
     argu2=argu*argu 
     beya(3)=3.d0*eta*argu2*(1.d0+argu2*(-40.d0+argu2)/640.d0)/8.d0 
     argu=4.d0*elkep(2) 
     argu2=argu*argu 
     beya(4)=eta*argu2*argu*(1.d0+argu2*(-48.d0+argu2)/960.d0)/12.d0 
     argu=5.d0*elkep(2) 
     argu2=argu*argu 
     beya(5)=5.d0*eta*argu2*argu2*(1.d0-argu2/24.d0)/384.d0 
     argu=6.d0*elkep(2) 
     argu2=argu*argu 
     beya(6)=eta*argu2*argu2*argu*(1.d0-argu2/28.d0)/640.d0 
     argu=7.d0*elkep(2) 
     argu2=argu*argu 
     beya(7)=7.d0*eta*argu2*argu2*argu2/46080.d0 
! -- \alpha_1(x) ... \alpha_7(x) functions                              
     argu=elkep(2) 
     argu2=argu*argu 
     alya(1)=1.d0+argu2*(-3456.d0+argu2*(240.d0-7.d0*argu2))/9216.d0 
     argu=2.d0*elkep(2) 
     argu2=argu*argu 
     alya(2)=argu*(1.d0+argu2*(-960.d0+argu2*(45.d0-argu2))/5760.d0) 
     argu=3.d0*elkep(2) 
     argu2=argu*argu 
     alya(3)=3.d0*argu2*(1.d0+argu2*(-200.d0+7.d0*argu2)/1920.d0)/  &
     &           8.d0                                                   
     argu=4.d0*elkep(2) 
     argu2=argu*argu 
     alya(4)=argu*argu2*(1.d0+argu2*(-36.d0+argu2)/480.d0)/12.d0 
     argu=5.d0*elkep(2) 
     argu2=argu*argu 
     alya(5)=argu2*argu2*(1.d0-7.d0*argu2/120.d0)/76.8d0 
     argu=6.d0*elkep(2) 
     argu2=argu*argu 
     alya(6)=argu*argu2*argu2*(1.d0-argu2/21.d0)/640.d0 
     argu=7.d0*elkep(2) 
     argu2=argu*argu 
     alya(7)=7.d0*argu2*argu2*argu2/46080.d0 
! close the input file                                                  
     call filclo(unit,' ') 
  ELSE 
     WRITE(*,*)' Yarkovsky datafile not found:',file(1:le) 
     stop 
  ENDIF
  WRITE(*,*)' Yarkovsky data loaded for asteroid ', astnam 
  RETURN 
111 yarfil=.false. 
  WRITE(*,*)' incomplete yarkovsky file for asteroid ', astnam 
END SUBROUTINE yarkinit

END MODULE yark_pert

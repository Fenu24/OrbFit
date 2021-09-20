! Copyright (C) 1998 by Zoran Knezevic (zoran@aob.bg.ac.yu)             
! Version: November, 1998                                               
!                                                                       
! Subroutine VAISUB                                                     
!                                                                       
! Contains the logic and computation of auxiliaries; calls two          
! included auxiliary subroutines: compdif and zeroin                    
!                                                                       
! Input is mostly contained in common (comvais.h)                       
!         alpha    - right ascension (rad)                              
!         delta    - declination (rad)                                  
!         xt,yt,zt - rectangular equatorial topocentric                 
!                    coordinates of the Sun                             
!         tjd0,tjd - uncorrected and light-time corrected obs. times    
!         xi,eta,zeta - direction cosines                               
!         ifo      - interactivity control 1=yes(fitobs), 2=no(orbfit)  
!                                                                       
! Output: geocentric distance at the second instant: rho                
!         difference of the z-components of the velocity: zddif         
!         coordinates at t2: x2,y2,z2                                   
!         velocity components at t2: x2p,y2p,z2p                        
!         time vector in MJD (TDT): tjd                                 
!         logical control for the succes of computation: fail           
!                                                                       
! ==============================================================        
      subroutine vaisub(ifo,rho,zddif,x2,y2,z2,x2p,y2p,z2p,fail) 
      implicit none 
      include 'comvais.h90' 
      include 'parrho.h90' 
!                                                                       
      integer j,ic,ifo 
!                                                                       
      double precision d(3),rho2(3),rho,zddif 
      double precision cb1,cb3,cc,x2,y2,z2,x2p,y2p,z2p 
      character*1 ans 
!                                                                       
      logical fail 
! ======================================================                
! direction cosines                                                     
      do 1 j=1,3 
        xi(j)=cos(alpha(j))*cos(delta(j)) 
        eta(j)=sin(alpha(j))*cos(delta(j)) 
        zeta(j)=sin(delta(j)) 
        if(abs(xi(j)**2+eta(j)**2+zeta(j)**2-1.d0).gt.1.d-7)then 
          write(*,*)'VAISUB: Observations too near to great circle.' 
!!!!!!!!!!!!!!!???????????????????!!!!!!!!!!!!!!!!!!
          fail=.true. 
          return 
        endif 
    1 continue 
! auxiliaries                                                           
      cb1=xi(1)*yt(1)-eta(1)*xt(1) 
      cb3=xi(3)*yt(3)-eta(3)*xt(3) 
      cc=eta(1)*xi(3)-eta(3)*xi(1) 
! ======================================================                
! initializations (see parrho.h)                                        
! ======================================================                
! counter                                                               
      ic=1 
! initial interval initialization test (can be removed in final vers.)  
!      if(rho2ini1.eq.rho2ini2) stop 'VAISUB: wrong initialization.'    
      if(rho2ini1.lt.rho2ini2)then 
        rho2(1)=rho2ini1 
        rho2(2)=rho2ini2 
      else 
        rho2(1)=rho2ini2 
        rho2(2)=rho2ini1 
      endif 
! ============================================================          
! Vaisala starts here; approximation 1 for both initial                 
! values rho2(1), rho2(2)                                               
! ============================================================          
   99 continue 
      call compdif(ic,rho2,cb1,cb3,cc,                                  &
     &             d,x2,y2,z2,x2p,y2p,z2p,fail)                         
! something horrible happened in compdif?                               
      IF(fail)THEN 
         RETURN 
      ENDIF 
! set counter ic                                                        
      if(ic.eq.1)then 
        ic=2 
        goto 99 
      endif 
      if(ic.eq.2)then 
! check for the root in the initial interval; if not ...                
        if(d(2)*d(1).gt.0.d0)then 
! ... find where is it                                                  
          if(abs(d(2)).lt.abs(d(1)))then 
             write(*,*) 'VAISUB: Object geocen. dist. > ', rho2(2),' AU' 
             if(ifo.eq.1)then 
               write(*,*) 'Do you want me to continue (y/n)? ' 
               read(*,*) ans 
               if(ans.eq.'y')then 
! extend the initial interval by 10 AU, and try again                   
                 rho2(2)=rho2(2)+10.d0 
                 goto 99 
               else 
                 fail=.true. 
                 return 
               endif 
             elseif(ifo.eq.2)then 
               fail=.true. 
               return 
             endif 
          else 
             if(ifo.eq.1)then 
!             write(*,*) 'VAISUB: Object geocen. dist. < ', rho2(1),' AU
!             write(*,*) 'INSIDE THE SPHERE OF INFLUENCE OF THE EARTH'  
!               write(*,*) 'Do you want me to continue (y/n)? '         
!               read(*,*) ans                                           
!               if(ans.eq.'y')then                                      
! extend the initial interval by 0.0060 AU, and try again               
!                 rho2(1)=rho2(1)-6.d-3                                 
!                 if(rho2(1).lt.2.d-4)then                              
!                   write(*,*)'OBJECT HAS ALREADY HIT THE EARTH!!'      
!                   write(*,*)'Check the input data, and try again'     
!                   fail=.true.                                         
!                   return                                              
!                 endif                                                 
!                 ic=1                                                  
!                 goto 99                                               
!               else                                                    
                 fail=.true. 
                 return 
!               endif                                                   
             elseif(ifo.eq.2)then 
               fail=.true. 
               return 
             endif 
          endif 
        endif 
! root in the interval: find it                                         
        call zeroin(rho2(1),rho2(2),d(1),d(2),cb1,cb3,cc,zdacc,         &
     &     rho2(3),fail)                                                
        IF(fail)RETURN 
! final approximation                                                   
        ic=3 
        goto 99 
      endif 
      rho=rho2(3) 
      zddif=d(3) 
      return 
      END                                           
! ==========================================================            
!                                                                       
! Subroutine COMPDIF                                                    
!                                                                       
! All the computations, except for the auxiliaries, are                 
! performed in this subroutine                                          
!                                                                       
! Input:                                                                
!     ic         -  iterations control                                  
!     rho2       -  trial geocentric distances                          
!     cb1,cb2,cc -  auxiliaries for the first and third instant         
!     d          -  difference of the z-components of velocity          
!     x2,y2,z2,x2p,y2p,z2p - see above                                  
! ==========================================================            
      subroutine compdif(ic,rho2,cb1,cb3,cc,                            &
     &             d,x2,y2,z2,x2p,y2p,z2p,fail)                         
      implicit none 
      LOGICAl fail 
      include 'comvais.h90' 
!                                                                       
      integer jc,j,ic 
!                                                                       
      double precision tau(3),rho2(3),rho1,rho3 
      double precision d(3),f(3),g(3) 
      double precision x2,y2,z2,r2,x2p,y2p,z2p 
      double precision sv2,s2,c1,c2,a2,a3,a4,a5,b3,b4,b5 
      double precision ca1,ca3,cb1,cb3,cc,z1,z3,z2p1,z2p3 
! ======================================================                
! parameters                                                            
      double precision gk,tlt 
      parameter(gk=1.720209895d-2) 
      parameter(tlt=5.78d-3) 
      SAVE 
! initializations                                                       
      jc=1 
   99 continue 
! time intervals; light time only when rho1/3 are known                 
      if(jc.eq.1)then 
        tjd(1)=tjd0(1) 
        tjd(2)=tjd0(2) 
        tjd(3)=tjd0(3) 
      elseif(jc.eq.2)then 
        tjd(1)=tjd0(1)-tlt*rho1 
        tjd(2)=tjd0(2)-tlt*rho2(ic) 
        tjd(3)=tjd0(3)-tlt*rho3 
      endif 
      do 4 j=1,3 
        tau(j)=gk*(tjd(j)-tjd(2)) 
    4 continue 
      IF(tau(3)-tau(2).eq.0.d0.or.tau(2)-tau(1).eq.0.d0.or.             &
     &  tau(3)-tau(2).eq.0.d0)THEN                                      
        fail=.true. 
        WRITE(*,*)'compdif: tau ',tau 
        RETURN 
      ELSE 
        fail=.false. 
      ENDIF 
! radius vector at the second instant                                   
      x2=rho2(ic)*xi(2)-xt(2) 
      y2=rho2(ic)*eta(2)-yt(2) 
      z2=rho2(ic)*zeta(2)-zt(2) 
      r2=sqrt(x2**2+y2**2+z2**2) 
! f,g series; to degree 5 only in second passage                        
      a2=0.5d0/r2**3 
      b3=a2/3.d0 
      if(jc.eq.2)then 
        sv2=x2p**2+y2p**2+z2p**2 
        s2=x2*x2p+y2*y2p+z2*z2p 
        c1=s2/r2**2 
        c2=0.25d0*(sv2-1.d0/r2-(s2/r2)**2)/r2**2 
        a3=c1*a2 
        a4=a2*(b3/2.d0+c2-c1**2) 
        a5=-c1*(2.d0*c1*a3+3.d0*a4) 
        b4=a3/2.d0 
        b5=(3.d0*a4-a2*b3)/5.d0 
      endif 
      do 5 j=1,3 
        f(j)=1.d0-a2*tau(j)**2 
        g(j)=tau(j)-b3*tau(j)**3 
        if(jc.eq.2)then 
          f(j)=f(j)+a3*tau(j)**3+a4*tau(j)**4+a5*tau(j)**5 
          g(j)=g(j)+b4*tau(j)**4+b5*tau(j)**5 
        endif 
    5 continue 
! velocity components for the second instant                            
      ca1=((xi(1)*y2-eta(1)*x2)*f(1)+cb1)/g(1) 
      ca3=((xi(3)*y2-eta(3)*x2)*f(3)+cb3)/g(3) 
      x2p=(ca1*xi(3)-ca3*xi(1))/cc 
! topocentric distances for the first and third instants                
      rho1=(f(1)*x2+g(1)*x2p+xt(1))/xi(1) 
      rho3=(f(3)*x2+g(3)*x2p+xt(3))/xi(3) 
! z-components of velocity and their difference                         
      z1=rho1*zeta(1)-zt(1) 
      z3=rho3*zeta(3)-zt(3) 
      z2p1=(z1-f(1)*z2)/g(1) 
      z2p3=(z3-f(3)*z2)/g(3) 
      d(ic)=z2p3-z2p1 
! velocity components in f,g series from the third approximation only   
      y2p=(ca1*eta(3)-ca3*eta(1))/cc 
      z2p=z2p1+g(3)*d(ic)/(g(3)-g(1)) 
! recompute with more accurate formulae                                 
      if(jc.eq.1)then 
        jc=2 
        goto 99 
      endif 
      return 
      END                                           
! ===================================================================   
!                                                                       
! Subroutine ZEROIN                                                     
!                                                                       
! Finds the root contained in the interval, with a given accuracy       
!                                                                       
! Input: boundaries of the initial interval: ax,bx                      
!        values of the function at interval boundaries: f1,f2           
!        auxiliaries for the computation of function: cb1,cb3,cc        
!        the size of the final interval containing the root: tol        
!                                                                       
! Output: rho- the root, i.e. geocentric dist. at the second instant    
!                                                                       
! Courtesy: Gojko Djurasevic, Astronomical Observatory, Belgrade        
! Adjusted by Zoran Knezevic, 1997                                      
!====================================================================   
      subroutine zeroin(ax,bx,f1,f2,cb1,cb3,cc,tol,                     &
     &      rho,fail)                                                   
      implicit none 
      LOGICAL fail 
      include 'comvais.h90' 
!                                                                       
      integer ic 
!                                                                       
      double precision ax,bx,f1,f2,cb1,cb3,cc,tol,rho 
      double precision dz(3),rho2(3),x2,y2,z2,x2p,y2p,z2p 
      double precision eps,tol1,a,b,c,d,e,fa,fb,fc,xm,p,q,r,s 
      eps=1.d0 
   10 eps=eps/2.d0 
      tol1=1.d0+eps 
      if(tol1.gt.1.d0) goto 10 
      a=ax 
      b=bx 
      fa=f1 
      fb=f2 
   20 c=a 
      fc=fa 
      d=b-a 
      e=d 
   30 if(abs(fc).ge.abs(fb)) goto 40 
      a=b 
      b=c 
      c=a 
      fa=fb 
      fb=fc 
      fc=fa 
   40 tol1=2.d0*eps*abs(b)+0.5d0*tol 
      xm=0.5d0*(c-b) 
      if(abs(xm).le.tol1) goto 90 
      if(fb.eq.0.d0) goto 90 
      if(abs(e).lt.tol1) goto 70 
      if(abs(fa).le.abs(fb)) goto 70 
      if(a.ne.c) goto 50 
      s=fb/fa 
      p=2.d0*xm*s 
      q=1.d0-s 
      goto 60 
   50 q=fa/fc 
      r=fb/fc 
      s=fb/fa 
      p=s*(2.d0*xm*q*(q-r)-(b-a)*(r-1.d0)) 
      q=(q-1.d0)*(r-1.d0)*(s-1.d0) 
   60 if(p.gt.0.d0) q=-q 
      p=abs(p) 
      if((2.d0*p).ge.(3.d0*xm*q-abs(tol1*q))) goto 70 
      if(p.ge.abs(0.5d0*e*q)) goto 70 
      e=d 
      d=p/q 
      goto 80 
   70 d=xm 
      e=d 
   80 a=b 
      fa=fb 
      if(abs(d).gt.tol1) then 
        b=b+d 
      else 
        b=b+sign(tol1,xm) 
      endif 
      ic=3 
      rho2(3)=b 
      call compdif(ic,rho2,cb1,cb3,cc,                                  &
     &             dz,x2,y2,z2,x2p,y2p,z2p,fail)                        
      fb=dz(3) 
      if((fb*(fc/abs(fc))).gt.0.d0) goto 20 
      goto 30 
   90 rho=b 
      return 
      END                                           

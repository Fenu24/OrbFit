MODULE terms9
 USE controlmod
 IMPLICIT NONE
 PRIVATE

! public routines
 PUBLIC terms, rgener, trigfu

CONTAINS

! ******************************                                        
!  terms - subroutines used to compute long periodic                  
!  perturbations; copyright Milani and Knezevic 1992                    
!  version with z2; April 1993                                          
! =================================================================**** 
!                                                                       
!   t e r m s (vers. 7.0 and later)                                     
!   computation of secular perturbations resulting from degree 4 terms  
!                                                                       
! =================================================================**** 
SUBROUTINE terms(al,r1,r2,s1,s2,n5,n6,m5,m6,x5,x6,e5,e6,          &
     &    sld,cld,ya,yb,de,pm,ff)                                       
  IMPLICIT NONE
! INPUT
  DOUBLE PRECISION, INTENT(IN) :: al,r1,r2,s1,s2,n5,n6,m5,m6,x5,x6,e5,e6
! =================================================================**** 
!     trigonometric functions.                                          
!  SLD(i,j,k,m)=\sin(i\lambda_5+j\lambda_6+k\dlta_5+m\dlta_6)           
!                                                                       
!  CLD(i,j,k,m)=\cos(i\lambda_5+j\lambda_6+k\dlta_5+m\dlta_6)           
  DOUBLE PRECISION, INTENT(IN) :: cld(-2:2,-2:2,-2:2,-2:2),sld(-2:2,-2:2,-2:2,-2:2) 
!  laplace, leverrier coefficients, divisors                            
  DOUBLE PRECISION, INTENT(IN) :: ya(4),yb(28),de(ndex) 
!OUTPUT
!  divisors, amplitudes
  DOUBLE PRECISION, INTENT(OUT) :: ff(4,ndex) 
! END INTERFACE
! =================================================================**** 
! yuasa coefficients, planetary masses
  DOUBLE PRECISION yu(28),pm
!  arrays for precomputation of powers                                  
  DOUBLE PRECISION x5v(0:3),x6v(0:3),e5v(0:3),e6v(0:3) 
  DOUBLE PRECISION n5v(0:3),n6v(0:3),m5v(0:3),m6v(0:3) 
  DOUBLE PRECISION r1v(0:3),r2v(0:3),s1v(0:3),s2v(0:3) 
! scalar temporaries
  DOUBLE PRECISION tf, trm
! loop indexes
  INTEGER i,j1
!   Adjustment of some coefficients appearing in generating function    
  yu(1:27)=yb(1:27)      
  yu(1)=yb(1)-ya(1)/4.d0 
  yu(3)=yb(3)-ya(2)/4.d0 
  yu(5)=yb(5)+ya(2)/2.d0 
  yu(8)=yb(8)-ya(3)/8.d0 
  yu(10)=yb(10)+ya(3)/4.d0 
  yu(12)=yb(12)-ya(4)/8.d0 
!    z2 term from second order                                          
  yu(28)=(yu(5)*2.d0*x6*e6+yu(10)*x6*m6+yb(14)*e6*n6+yb(16)*m6*n6)* &
     &   (((2.d0*yu(1)+yu(5))*x6+(yu(12)+yb(14))*n6/2.d0)*(1.d0/de(10)  &
     &   -1.d0/de(2))-(2*pm/al)*((yu(1)*4.d0*x6+yu(12)*n6)*(yu(1)*2.d0+ &
     &   yu(5))*(r1**2+r2**2)+(yu(5)*2.d0*x6+yb(14)*n6)*(yu(5)+yu(3)*   &
     &   2.d0)*(s1**2+s2**2))*(1.d0/de(10)**2)+(pm/al)*(yu(1)*4.d0*x6*  &
     &   (r1**2+r2**2+2.d0*x5**2+x6**2)+yu(5)*2.d0*x6*(s1**2+s2**2+     &
     &   e5**2+e6**2)+yu(12)*(n6*(r1**2+r2**2+2.d0*x5**2+3.d0*x6**2)+   &
     &   4.d0*n5*x5*x6)+yb(14)*n6*(s1**2+s2**2+e5**2+e6**2))*(yu(1)+    &
     &   yu(5))*(1.d0/de(2)**2))                                        
  yu(28)=yu(28)+(yu(1)*2.d0*x6*x6+yu(12)*x6*n6+yb(19)*n6*n6)*       &
     &   ((yu(5)*2.d0*e6+yu(10)*m6)*(1.d0/de(57)-1.d0/de(5))-(12.d0*pm/ &
     &   al)*yu(1)*(r1**2+r2**2)*(1.d0/de(57)**2)+(pm/al)*(yu(3)*       &
     &   4.d0*e6*(s1**2+s2**2+2.d0*e5**2+e6**2)+yu(5)*2.d0*e6**2*(r1**2 &
     &   +r2**2+x5**2+x6**2)+yu(8)*(m6*(s1**2+s2**2+2.d0*e5**2+3.d0*    &
     &   e6**2)+4.d0*m5*e5*e6)+yu(10)*m6*(r1**2+r2**2+x5**2+x6**2))*    &
     &   yu(5)*(1.d0/de(5)**2))                                         
  yu(28)=yu(28)+(yb(20)*2.d0*e6+yb(22)*m6)*(yb(20)*x6*x6+yb(23)*    &
     &   x6*n6+yb(26)*n6*n6)*((1.d0/de(48)-1.d0/de(24))-(4.d0*pm/al)*   &
     &   (s1**2+s2**2)*(yu(3)*3.d0*(1.d0/de(48)**2)-(yu(5)-yu(3))*2.d0  &
     &   *(1.d0/de(24)**2)))                                            
  yu(28)=-yu(28)*pm 
! precomputation of all useful powers; all up to degree 3               
  x5v(0)=1.d0 
  x6v(0)=1.d0 
  e5v(0)=1.d0 
  e6v(0)=1.d0 
  n5v(0)=1.d0 
  n6v(0)=1.d0 
  m5v(0)=1.d0 
  m6v(0)=1.d0 
  r1v(0)=1.d0 
  r2v(0)=1.d0 
  s1v(0)=1.d0 
  s2v(0)=1.d0 
  DO i=1,3 
     x5v(i)=x5v(i-1)*x5 
     x6v(i)=x6v(i-1)*x6 
     e5v(i)=e5v(i-1)*e5 
     e6v(i)=e6v(i-1)*e6 
     m5v(i)=m5v(i-1)*m5 
     m6v(i)=m6v(i-1)*m6 
     n5v(i)=n5v(i-1)*n5 
     n6v(i)=n6v(i-1)*n6 
     r1v(i)=r1v(i-1)*r1 
     r2v(i)=r2v(i-1)*r2 
     s1v(i)=s1v(i-1)*s1 
     s2v(i)=s2v(i-1)*s2 
  ENDDO
!  for dynamic variables, degree 4 cannnot occur                        
!        r1v(4)=r1v(3)*r1                                               
!        r2v(4)=r2v(3)*r2                                               
!        s1v(4)=s1v(3)*s1                                               
!        s2v(4)=s2v(3)*s2                                               
! perturbation accumulators are set to zero                             
   ff=0.d0 
! computation begins here                                               
   do 2 j1=1,nig 
      if(ig(15,j1).gt.1) then 
         tf=cld(ig(16,j1),ig(17,j1),ig(18,j1),ig(19,j1)) 
      else 
         tf=sld(ig(16,j1),ig(17,j1),ig(18,j1),ig(19,j1)) 
      endif
      trm=yu(ig(1,j1))*ig(2,j1)                                       &
     &  *x5v(ig(3,j1))*x6v(ig(4,j1))                                    &
     &  *e5v(ig(5,j1))*e6v(ig(6,j1))                                    &
     &  *n5v(ig(7,j1))*n6v(ig(8,j1))                                    &
     &  *m5v(ig(9,j1))*m6v(ig(10,j1))                                   &
     &  *r1v(ig(11,j1))*r2v(ig(12,j1))                                  &
     &  *s1v(ig(13,j1))*s2v(ig(14,j1))*tf                               
      if(ig(11,j1).gt.0) then 
         ff(1,ig(20,j1))=ff(1,ig(20,j1))+ig(11,j1)/r1*trm 
      endif
      if(ig(12,j1).gt.0) then 
         ff(2,ig(20,j1))=ff(2,ig(20,j1))+ig(12,j1)/r2*trm 
      endif
      if(ig(13,j1).gt.0) then 
         ff(3,ig(20,j1))=ff(3,ig(20,j1))+ig(13,j1)/s1*trm 
      endif
      if(ig(14,j1).gt.0) then 
         ff(4,ig(20,j1))=ff(4,ig(20,j1))+ig(14,j1)/s2*trm 
      endif
2  ENDDO
 END SUBROUTINE terms
! =============================                                         
!  trigfu (vers. 6.8 and later)                                         
! =================================================================**** 
!     trigonometric functions; FFT version                              
! =================================================================**** 
 SUBROUTINE trigfu(l5,l6,d5,d6,sld,cld) 
   IMPLICIT NONE
   DOUBLE PRECISION, INTENT(IN) :: l5,l6,d5,d6 ! secular angles J+S
! ************************************************************          
! SLD(i,j,k,m)=\sin(i\lambda_5+j\lambda_6+k\delta_5+m\delta_6)          
!                                                                       
! CLD(i,j,k,m)=\cos(i\lambda_5+j\lambda_6+k\delta_5+m\delta_6)          
! ************************************************************          
   DOUBLE PRECISION, INTENT(OUT) :: cld(-2:2,-2:2,-2:2,-2:2),sld(-2:2,-2:2,-2:2,-2:2) 
! END INTERFACE
   DOUBLE PRECISION sl5,cl5,sl6,cl6,sd5,cd5,sd6,cd6 ! sine/cosine
   sl5=sin(l5) 
   cl5=cos(l5) 
   sl6=sin(l6) 
   cl6=cos(l6) 
   cld(0,0,0,0)=1.d0 
   sld(1,0,0,0)=sl5 
   cld(1,0,0,0)=cl5 
   sld(0,1,0,0)=sl6 
   cld(0,1,0,0)=cl6 
   sld(1,-1,0,0)=sl5*cl6-cl5*sl6 
   cld(1,-1,0,0)=cl5*cl6+sl5*sl6 
   sld(2,0,0,0)=2*cl5*sl5 
   cld(2,0,0,0)=cl5*cl5-sl5*sl5 
   sld(0,2,0,0)=2*cl6*sl6 
   cld(0,2,0,0)=cl6*cl6-sl6*sl6 
   sld(2,-1,0,0)=sl5*cld(1,-1,0,0)+cl5*sld(1,-1,0,0) 
   cld(2,-1,0,0)=cl5*cld(1,-1,0,0)-sl5*sld(1,-1,0,0) 
   sld(-1,2,0,0)=-sl5*cld(0,2,0,0)+cl5*sld(0,2,0,0) 
   cld(-1,2,0,0)=cl5*cld(0,2,0,0)+sl5*sld(0,2,0,0) 
   sld(1,1,0,0)=sl5*cl6+cl5*sl6 
   cld(1,1,0,0)=cl5*cl6-sl5*sl6 
   sd5=sin(d5) 
   cd5=cos(d5) 
   sd6=sin(d6) 
   cd6=cos(d6) 
   sld(0,0,1,0)=sd5 
   cld(0,0,1,0)=cd5 
   sld(0,0,0,1)=sd6 
   cld(0,0,0,1)=cd6 
   sld(0,0,1,-1)=sd5*cd6-sd6*cd5 
   cld(0,0,1,-1)=cd5*cd6+sd5*sd6 
   sld(0,0,2,0)=2*sd5*cd5 
   cld(0,0,2,0)=cd5*cd5-sd5*sd5 
   sld(0,0,0,2)=2*sd6*cd6 
   cld(0,0,0,2)=cd6*cd6-sd6*sd6 
   sld(0,0,1,1)=sd5*cd6+sd6*cd5 
   cld(0,0,1,1)=cd5*cd6-sd5*sd6 
   sld(0,0,2,-1)=sd5*cld(0,0,1,-1)+cd5*sld(0,0,1,-1) 
   cld(0,0,2,-1)=cd5*cld(0,0,1,-1)-sd5*sld(0,0,1,-1) 
   sld(0,0,-1,2)=-sd5*cld(0,0,0,2)+cd5*sld(0,0,0,2) 
   cld(0,0,-1,2)=cd5*cld(0,0,0,2)+sd5*sld(0,0,0,2) 
   sld(1,0,1,0)=sl5*cd5+cl5*sd5 
   cld(1,0,1,0)=cl5*cd5-sl5*sd5 
   sld(1,0,-1,0)=sl5*cd5-cl5*sd5 
   cld(1,0,-1,0)=cl5*cd5+sl5*sd5 
   sld(-1,0,1,0)=-sld(1,0,-1,0) 
   cld(-1,0,1,0)=cld(1,0,-1,0) 
   sld(0,1,1,0)=sl6*cd5+cl6*sd5 
   cld(0,1,1,0)=cl6*cd5-sl6*sd5 
   sld(0,1,-1,0)=sl6*cd5-cl6*sd5 
   cld(0,1,-1,0)=cl6*cd5+sl6*sd5 
   sld(0,-1,1,0)=-sld(0,1,-1,0) 
   cld(0,-1,1,0)=cld(0,1,-1,0) 
   sld(1,0,0,1)=sl5*cd6+cl5*sd6 
   cld(1,0,0,1)=cl5*cd6-sl5*sd6 
   sld(1,0,0,-1)=sl5*cd6-cl5*sd6 
   cld(1,0,0,-1)=cl5*cd6+sl5*sd6 
   sld(-1,0,0,1)=-sld(1,0,0,-1) 
   cld(-1,0,0,1)=cld(1,0,0,-1) 
   sld(0,1,0,1)=sl6*cd6+cl6*sd6 
   cld(0,1,0,1)=cl6*cd6-sl6*sd6 
   sld(0,1,0,-1)=sl6*cd6-cl6*sd6 
   cld(0,1,0,-1)=cl6*cd6+sl6*sd6 
   sld(0,-1,0,1)=-sld(0,1,0,-1) 
   cld(0,-1,0,1)=cld(0,1,0,-1) 
   sld(2,0,-1,0)=sl5*cld(1,0,-1,0)+cl5*sld(1,0,-1,0) 
   cld(2,0,-1,0)=cl5*cld(1,0,-1,0)-sl5*sld(1,0,-1,0) 
   sld(2,0,0,-1)=sl5*cld(1,0,0,-1)+cl5*sld(1,0,0,-1) 
   cld(2,0,0,-1)=cl5*cld(1,0,0,-1)-sl5*sld(1,0,0,-1) 
   sld(0,2,-1,0)=sl6*cld(0,1,-1,0)+cl6*sld(0,1,-1,0) 
   cld(0,2,-1,0)=cl6*cld(0,1,-1,0)-sl6*sld(0,1,-1,0) 
   sld(0,2,0,-1)=sl6*cld(0,1,0,-1)+cl6*sld(0,1,0,-1) 
   cld(0,2,0,-1)=cl6*cld(0,1,0,-1)-sl6*sld(0,1,0,-1) 
   sld(-1,0,2,0)=sd5*cld(-1,0,1,0)+cd5*sld(-1,0,1,0) 
   cld(-1,0,2,0)=cd5*cld(-1,0,1,0)-sd5*sld(-1,0,1,0) 
   sld(-1,0,0,2)=sd6*cld(-1,0,0,1)+cd6*sld(-1,0,0,1) 
   cld(-1,0,0,2)=cd6*cld(-1,0,0,1)-sd6*sld(-1,0,0,1) 
   sld(0,-1,2,0)=sd5*cld(0,-1,1,0)+cd5*sld(0,-1,1,0) 
   cld(0,-1,2,0)=cd5*cld(0,-1,1,0)-sd5*sld(0,-1,1,0) 
   sld(0,-1,0,2)=sd6*cld(0,-1,0,1)+cd6*sld(0,-1,0,1) 
   cld(0,-1,0,2)=cd6*cld(0,-1,0,1)-sd6*sld(0,-1,0,1) 
   sld(1,0,-1,1)=sd6*cld(1,0,-1,0)+cd6*sld(1,0,-1,0) 
   cld(1,0,-1,1)=cd6*cld(1,0,-1,0)-sd6*sld(1,0,-1,0) 
   sld(1,0,1,-1)=sd5*cld(1,0,0,-1)+cd5*sld(1,0,0,-1) 
   cld(1,0,1,-1)=cd5*cld(1,0,0,-1)-sd5*sld(1,0,0,-1) 
   sld(0,1,-1,1)=sd6*cld(0,1,-1,0)+cd6*sld(0,1,-1,0) 
   cld(0,1,-1,1)=cd6*cld(0,1,-1,0)-sd6*sld(0,1,-1,0) 
   sld(0,1,1,-1)=sd5*cld(0,1,0,-1)+cd5*sld(0,1,0,-1) 
   cld(0,1,1,-1)=cd5*cld(0,1,0,-1)-sd5*sld(0,1,0,-1) 
   sld(-1,0,1,1)=sd6*cld(-1,0,1,0)+cd6*sld(-1,0,1,0) 
   cld(-1,0,1,1)=cd6*cld(-1,0,1,0)-sd6*sld(-1,0,1,0) 
   sld(0,-1,1,1)=sd6*cld(0,-1,1,0)+cd6*sld(0,-1,1,0) 
   cld(0,-1,1,1)=cd6*cld(0,-1,1,0)-sd6*sld(0,-1,1,0) 
   sld(-1,1,1,0)=sl6*cld(-1,0,1,0)+cl6*sld(-1,0,1,0) 
   cld(-1,1,1,0)=cl6*cld(-1,0,1,0)-sl6*sld(-1,0,1,0) 
   sld(1,-1,1,0)=sl5*cld(0,-1,1,0)+cl5*sld(0,-1,1,0) 
   cld(1,-1,1,0)=cl5*cld(0,-1,1,0)-sl5*sld(0,-1,1,0) 
   sld(-1,1,0,1)=sl6*cld(-1,0,0,1)+cl6*sld(-1,0,0,1) 
   cld(-1,1,0,1)=cl6*cld(-1,0,0,1)-sl6*sld(-1,0,0,1) 
   sld(1,-1,0,1)=sl5*cld(0,-1,0,1)+cl5*sld(0,-1,0,1) 
   cld(1,-1,0,1)=cl5*cld(0,-1,0,1)-sl5*sld(0,-1,0,1) 
   sld(1,1,-1,0)=sl5*cld(0,1,-1,0)+cl5*sld(0,1,-1,0) 
   cld(1,1,-1,0)=cl5*cld(0,1,-1,0)-sl5*sld(0,1,-1,0) 
   sld(1,1,0,-1)=sl5*cld(0,1,0,-1)+cl5*sld(0,1,0,-1) 
   cld(1,1,0,-1)=cl5*cld(0,1,0,-1)-sl5*sld(0,1,0,-1) 
!  for z2 terms, some deg. 6 arguments                                  
   sld(0,2,0,1)=sld(0,2,0,0)*cd6+cld(0,2,0,0)*sd6 
   cld(0,2,0,1)=cld(0,2,0,0)*cd6-sld(0,2,0,0)*sd6 
 END SUBROUTINE trigfu
! =======================================================               
!  rgener                                                               
!  read integer lists (generic term representation) of the              
!  determining function                                                 
 SUBROUTINE rgener(iun)
   INTEGER, INTENT(IN) :: iun ! input unit    
   INTEGER jg,k 
   read(iun,199) 
199 format(//////////) 
   nig=0 
   DO jg=1,igx 
      read(iun,202,end=25)(ig(k,jg),k=1,20) 
202   format(20i3) 
      nig=nig+1 
   ENDDO
   write(*,*)' too many terms in gener.lng; max was ', igx 
   stop 
25 continue 
 END SUBROUTINE rgener

END MODULE terms9

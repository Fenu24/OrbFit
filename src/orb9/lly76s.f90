! ********************************************************************* 
!   llyco                                                               
!   interface rotine for lly                                            
!   Laplace's, Leverrier's and Yuasa's coefficients                     
! ********************************************************************* 
SUBROUTINE llyco(a,pa,an,pn,sm,ya,yb,pm,iord)
  IMPLICIT NONE
! INPUT
  INTEGER, INTENT(IN) :: iord ! order of computation, set to 1 for first order
! semimajor axis: asteroid, planet; mean motion: asteroid, planet, mass: sun, planet
  DOUBLE PRECISION, INTENT(IN) :: a, pa, an, pn, sm, pm
! OUTPUT
!  yuasa coefficients  a and b                                          
  DOUBLE PRECISION, INTENT(OUT) :: ya(4),yb(27) 
! END INTERFACE
  INTEGER nalfa
  DOUBLE PRECISION apr, alfa, sa
  INTEGER np
!  order control                                                        
  if(iord.eq.0)return 
!  inner/outer perturbation cases handled by Williams (1969) trick      
  if(a.lt.pa) then 
     alfa=a/pa 
     apr=pa 
  else 
     alfa=pa/a 
     apr=a 
  endif
  sa=sqrt(a) 
  np=nalfa(alfa) 
  call lly(a,alfa,apr,sa,sm,an,pn,np,ya,yb,pm,iord)
END SUBROUTINE llyco
! **************************************************************        
!     (l)aplace,(l)everrier,(y)uasa (co)efficients                      
SUBROUTINE lly(aa,a,apr,sa,sm,an,pn,np,ya,yb,pm,io)
  IMPLICIT NONE
! INPUT
  INTEGER, INTENT(IN) :: io ! order of computation (io=1 to avoid second
! order computation
! semimajor axis: value and ratios, square root 
! mean motion: asteroid, planet, mass: sun, planet
  DOUBLE PRECISION, INTENT(IN) :: aa, a, apr, sa, an, pn, sm, pm
  INTEGER, INTENT(IN) :: np ! number of terms to be computed
! OUPUT
!  Yuasa coefficients  a and b                                          
  DOUBLE PRECISION, INTENT(OUT) :: ya(4),yb(27) 
! END INTERFACE
  DOUBLE PRECISION b(180),c(180),e(5) 
  DOUBLE PRECISION b1(180),b2(180),b3(180),b4(5),c1(180),c2(5) 
  DOUBLE PRECISION cl(38,361), ci, d11, d12, d13, d14, d15, d31, d32, d33, d34 
  DOUBLE PRECISION yy(27), d61, d62, d63, d64, d65, d66 
  DOUBLE PRECISION di11, di12, di13, di14, di15, di16, di17, di31, di32
  DOUBLE PRECISION di33, di34, di35, di61, di62, di63, di64, di65, di66
  DOUBLE PRECISION di67, di68, di69, di610
  DOUBLE PRECISION yd1, yd3, yd5, yd6
  INTEGER k, k5, k6, nnp, i, j, k7, k8, j1, j2, j3, ll1, ll2
  INTEGER ll3, ll4, ll5, ll6
!     laplace's coefficients.                                           
  call lclong (a,np,b,c,e) 
  call lcderlon(a,np,b1,b2,b3,b4,c1,c2) 
!     yuasa's coefficients of the first order (yuasa, 1973; table 1.)   
!     ya(1)=a1; 2=3; 3=5; 4=6;                                          
!     yb(1)=b1; 2=2; 3=4; 4=6; 5=7; 6=8; 7=9; 8=11; 9=12; 10=13;        
!     11=14; 12=15; 13=16; 14=17; 15=18; 16=19; 17=20; 18=21; 19=22;    
!     20=23; 21=24; 22=25; 23=26; 24=27; 25=28; 26=29; 27=31;           
!                                                                       
  ya(1)=(b1(1)+b2(1)/2.d0)/4.d0 
  ya(2)=-a*c(2)/8.d0 
  ya(3)=-2.d0*ya(2) 
  ya(4)=(2.d0*b(2)-2.d0*b1(2)-b2(2))/4.d0 
  yb(1)=(4.d0*b3(1)+b4(1))/128.d0 
  yb(2)=(4.d0*b1(1)+14.d0*b2(1)+8.d0*b3(1)+b4(1))/32.d0 
  yb(3)=(-c(2)+3.d0*a/2.d0*(e(3)/2.d0+e(1)))*a/32.d0 
  yb(4)=(c(2)+3.d0*a*(e(3)/2.d0+e(1)))*a/16.d0 
  yb(5)=-(2.d0*c(2)+4.d0*c1(2)+c2(2))*a/32.d0 
  yb(6)=yb(5) 
  yb(7)=yb(5) 
  yb(8)=-(e(1)+e(3)/2.d0)*3.d0*a*a/16.d0 
  yb(9)=yb(8) 
  yb(10)=-2.d0*yb(5) 
  yb(11)=yb(10) 
  yb(12)=-(4.d0*b2(2)+6.d0*b3(2)+b4(2))/32.d0 
  yb(13)=(4.d0*b(2)-4.d0*b1(2)-22.d0*b2(2)-10.d0*b3(2)-b4(2))/32.d0 
  yb(14)=(4.d0*c1(1)+4.d0*c1(3)+c2(1)+c2(3))*a/32.d0 
  yb(15)=yb(14) 
  yb(16)=-yb(14)+ya(4)/4.d0 
  yb(17)=-yb(14)-ya(4)/4.d0 
  yb(18)=-yb(8)/2.d0 
  yb(19)=(b(3)-b1(3)+b2(3)/2.d0+2.d0*b3(3)/3.d0+                    &
     &b4(3)/12.d0)*3.d0/16.d0                                           
  yb(20)=(6.d0*c(2)+4.d0*c1(2)+c2(2)/2.d0)*a/32.d0 
  yb(21)=yb(20) 
  yb(22)=-2.d0*yb(20) 
  yb(23)=(2.d0*c1(1)+c2(1)/2.d0)*(-a/16.d0) 
  yb(24)=yb(23) 
  yb(25)=-yb(23)*2.d0 
  yb(26)=a*c2(2)/64.d0 
  yb(27)=-2.d0*yb(26) 
  DO k5=1,4 
     ya(k5)=ya(k5)/apr
  ENDDO 
  DO k6=1,27 
     yb(k6)=yb(k6)/apr
  ENDDO
! if first order computation, set to zero Yuasa's D coefficients        
  if(io.lt.2)then 
     return 
  endif
!     leverrier's coefficients.                                         
!     cl(1,k)=(1)^(i); 2=2^i; 3=2^i+1; 4=2^-i-1; 5=11^i; 6=21^i;        
!     7=21^i-1; 8=21^-i-1; 9=50^i; 10=50^i+1; 11=51^i; 12=60^i;         
!     13=70^i; 14=70^i-1; 15=71^i-1; 16=90^-i-1; 17=172^i;              
!     18=172^i+1; 19=182^i; 20=212^i;                                   
!     derivatives of leverrier's coefficients.                          
!     cl(21,k)=(1)_a^(i); 22=1_a^i-1; 23=1_a^-i-1; 24=2_a^i;            
!     25=11_a^i; 26=21_a^i-1; 27=21_a^-i-1; 28=50_a^i; 29=50_a^i+1;     
!     30=70_a^i-1; 31=1_a^i+1;'                                         
!     additional coefficients for fourth degree corrections             
!     32=52; 33=72; 34=80; 35=100; 36=120; 37=130; 38=192;              
  nnp=2*np-1 
  DO 80 k=1,nnp 
     i=k-np 
     ci=i*1.d0 
     j=iabs(i)+1 
     j1=iabs(i+1)+1 
     j2=iabs(-i-1)+1 
     j3=iabs(i-1)+1 
     cl(1,k)=b(j)/2.d0 
     cl(2,k)=-2.d0*ci*ci*b(j)+b1(j)+b2(j)/2.d0 
     cl(3,k)=-2.d0*(ci+1.d0)*(ci+1.d0)*b(j1)+b1(j1)+b2(j1)/2.d0 
     cl(4,k)=-2.d0*(-ci-1.d0)*(-ci-1.d0)*b(j2)+b1(j2)+b2(j2)/2.d0 
     cl(5,k)=-a/4.d0*(c(j3)+c(j1)) 
     cl(6,k)=(2.d0*ci+4.d0*ci*ci)*b(j)-2.d0*b1(j)-b2(j) 
     cl(7,k)=(2.d0-6.d0*ci+4.d0*ci*ci)*b(j3)-2.d0*b1(j3)-b2(j3) 
     cl(8,k)=(2.d0+6.d0*ci+4.d0*ci*ci)*b(j2)-2.d0*b1(j2)-b2(j2) 
     cl(9,k)=-2.d0*ci*b(j)-b1(j) 
     cl(10,k)=-2.d0*(ci+1.d0)*b(j1)-b1(j1) 
     cl(11,k)=(ci-5.d0*ci*ci+4.d0*ci*ci*ci)*b(j)+.5d0*                 &
     &(3.d0-7.d0*ci+4.d0*ci*ci)*b1(j)-(1.d0+ci)*b2(j)-b3(j)/2.d0        
     cl(12,k)=a*((.5d0+ci)*(c(j3)+c(j1))+.5d0*(c1(j3)+c1(j1))) 
     cl(13,k)=(1.d0+2.d0*ci)*b(j)+b1(j) 
     cl(14,k)=(2.d0*ci-1.d0)*b(j3)+b1(j3) 
     cl(15,k)=(4.d0-16.d0*ci+20.d0*ci*ci-8.d0*ci*ci*ci)*b(j3)+         &
     &(-4.d0+12.d0*ci-4.d0*ci*ci)*b1(j3)+(3.d0+2.d0*ci)*b2(j3)+b3(j3)   
     cl(16,k)=.5d0*((1.d0-ci-10.d0*ci*ci-8.d0*ci*ci*ci)*b(j2)-         &
     &(1.d0+ci+4.d0*ci*ci)*b1(j2)+(3.d0+2.d0*ci)*b2(j2)+b3(j2))         
     cl(17,k)=.5d0*((-5.d0*ci+4.d0*ci*ci)*b(j)+                        &
     &(-2.d0+4.d0*ci)*b1(j)+b2(j))                                      
     cl(18,k)=.5d0*((-1.d0+3.d0*ci+4.d0*ci*ci)*b(j1)+                  &
     &(2.d0+4.d0*ci)*b1(j1)+b2(j1))                                     
     cl(19,k)=(-2.d0*ci-4.d0*ci*ci)*b(j)-(2.d0+4.d0*ci)*b1(j)-b2(j) 
     cl(20,k)=.5d0*a*c(j3) 
     cl(21,k)=.5d0/aa*b1(j) 
     cl(22,k)=.5d0/aa*b1(j3) 
     cl(23,k)=.5d0/aa*b1(j2) 
     cl(24,k)=1.d0/aa*((1.d0-2.d0*ci*ci)*b1(j)+2.d0*b2(j)+b3(j)/2.d0) 
     cl(25,k)=-.25d0*a/aa*(c1(j3)+c(j3)+c1(j1)+c(j1)) 
     cl(26,k)=1.d0/aa*((-6.d0*ci+4.d0*ci*ci)*b1(j3)-4.d0*b2(j3)-b3(j3)) 
     cl(27,k)=1.d0/aa*((6.d0*ci+4.d0*ci*ci)*b1(j2)-4.d0*b2(j2)-b3(j2)) 
     cl(28,k)=-1.d0/aa*((1.d0+2.d0*ci)*b1(j)+b2(j)) 
     cl(29,k)=-1.d0/aa*((3.d0+2.d0*ci)*b1(j1)+b2(j1)) 
     cl(30,k)=1.d0/aa*(2.d0*ci*b1(j3)+b2(j3)) 
     cl(31,k)=.5d0/aa*b1(j1) 
     cl(32,k)=8.d0*ci*ci*ci*b(j)-(2.d0+4.d0*ci-4.d0*ci*ci)             &
     &*b1(j)-(4.d0+2.d0*ci)*b2(j)-b3(j)                                 
     cl(33,k)=.5d0*((-1.d0-5.d0*ci-14.d0*ci*ci-8.d0*ci*ci*ci)*         &
     &b(j)+(7.d0+ci-4.d0*ci*ci)*b1(j)+(7.d0+2.d0*ci)*b2(j)+b3(j))       
     cl(34,k)=.5d0*a*((-2.d0-2.d0*ci)*(c(j3)+c(j1))-(c1(j3)+           &
     &c1(j1)))                                                          
     cl(35,k)=.5d0*((8.d0*ci+18.d0*ci*ci+8.d0*ci*ci*ci)*b(j)-          &
     &(10.d0+ci-4.d0*ci*ci)*b1(j)-(8.d0+2.d0*ci)*b2(j)-b3(j))           
     cl(36,k)=.5d0*a*((2.d0*ci-5.d0)*c(j3)-c1(j3)) 
     cl(37,k)=.5d0*a*(2.d0*(1.d0-ci)*c(j3)+c1(j3)) 
     cl(38,k)=.5d0*((4.d0+9.d0*ci+4.d0*ci*ci)*b(j)+(6.d0+              &
     &4.d0*ci)*b1(j)+b2(j))                                             
     DO k7=1,38 
        cl(k7,k)=cl(k7,k)/apr
     ENDDO
80 ENDDO
!     yuasa's coefficients of the second order.                         
  yd1=0.d0 
  yd3=0.d0 
  yd5=0.d0 
  yd6=0.d0 
  DO 90 k=1,nnp 
     j=nnp+1-k 
     i=k-np 
     ci=i*1.d0 
     if (k .eq. np ) goto 100 
     d11=-1.d0/(pn-an)*(.25d0*sa/sm*(cl(1,k)*(cl(24,k)+cl(24,j))       &
     &+cl(2,k)*(cl(21,k)+cl(21,j)))-                                    &
     &1.d0/(8.d0*sa*sm)*cl(1,k)*(cl(2,k)+cl(2,j)))                      
     d14=3.d0/(2.d0*aa*aa)*(1/(8.d0*(pn-an)*(pn-an))*(cl(1,k)          &
     &*(cl(2,k)+cl(2,j))+cl(2,k)*(cl(1,k)+cl(1,j))))                    
     goto 110 
100  d11=0.d0 
     d14=0.d0 
110  d12=1.d0/(ci*pn-(ci-1.d0)*an)*(-(ci-1.d0)*sa/sm/4.d0*cl(28,k)*    &
     &cl(9,k)+1.d0/(8.d0*sa*sm)*cl(9,k)*cl(11,k)+                       &
     &(ci-2.d0)/(16.d0*sa*sm)*cl(9,k)*cl(9,k))                          
     d13=1.d0/((ci*pn-(ci-2.d0)*an)*8.d0*sa*sm)*cl(17,k)*cl(17,k) 
     d15=3.d0/(2.d0*aa*aa)*(ci-1.d0)*(ci-1.d0)/((ci*pn-(ci-1.d0)*an)*  &
     &(ci*pn-(ci-1.d0)*an)*8.d0)*cl(9,k)*cl(9,k)                        
     yd1=yd1+d11+d12+d13+d14+d15 
     if (k .eq. np) goto 120 
     d31=-1.d0/(pn-an)*(sa/sm/4.d0*(cl(1,k)*(cl(25,k)+                 &
     &cl(25,j))+cl(5,k)*(cl(21,k)+cl(21,j)))                            &
     &-1.d0/(8.d0*sa*sm)*cl(1,k)*(cl(5,k)+cl(5,j)))                     
     d34=3.d0/(2.d0*aa*aa*(pn-an)*(pn-an)*8.d0)*(cl(1,k)*(cl(5,k)      &
     &+cl(5,j))+cl(5,k)*(cl(1,k)+cl(1,j)))                              
     goto 130 
120  d31=0.d0 
     d34=0.d0 
130  d32=1.d0/(ci*pn-(ci-1.d0)*an)*1.d0/(16.d0*sa*sm)*cl(9,k)*cl(12,k) 
     d33=1.d0/((ci*pn-(ci-2.d0)*an)*8.d0*sa*sm)*cl(20,k)*cl(20,k) 
     yd3=yd3+d31+d32+d33+d34 
     if (k .eq. np) goto 140 
     d61=-1.d0/(pn-an)*(sa/4.d0/sm*cl(1,k)*(cl(26,k)+cl(27,k))-        &
     &1.d0/(16.d0*sa*sm)*cl(1,k)*(cl(7,k)+cl(8,k)))                     
     d65=3.d0/(2.d0*aa*aa*(pn-an)*(pn-an)*8.d0)*(3.d0*cl(1,k)          &
     &*cl(7,k)+cl(1,k)*cl(8,k))                                         
     goto 150 
140  d61=0.0 
     d65=0.0 
150  if (k .eq. (nnp-1)/2) goto 160 
     d62=-1.d0/((ci+1.d0)*(pn-an))*((ci+1.d0)*sa/4.d0/sm*cl(6,k)*      &
     &(cl(31,k)+cl(23,k))+1.d0/(16.d0*sa*sm)*cl(6,k)*(cl(3,k)+cl(4,k))) 
     goto 170 
160  d62=0.d0 
170  d63=1.d0/(ci*pn-(ci-1.d0)*an)*(-(ci-1.d0)*sa/4.d0/sm*(cl(30,k)*   &
     &cl(9,k)+cl(28,k)*cl(14,k))+1.d0/sa/sm*((ci-1.d0)/16.d0*cl(9,k)*   &
     &cl(14,k)+1.d0/16.d0*cl(9,k)*cl(15,k)+1.d0/8.d0*cl(9,k)*cl(16,k))) 
     d64=1.d0/(((ci+1.d0)*pn-(ci-1.d0)*an)*8.d0*sa*sm)*cl(18,k)*       &
     &cl(19,k)                                                          
     d66=3.d0*(1.d0-ci)*(1.d0-ci)/(2.d0*aa*aa*4.d0                     &
     &  *(ci*pn-(ci-1.d0)*an)*                                          &
     &(ci*pn-(ci-1.d0)*an))*cl(9,k)*cl(14,k)                            
     yd6=yd6+d61+d62+d63+d64+d65+d66
90 ENDDO
  ll1=(nnp-1)/2-1 
  ll2=ll1+1 
  ll3=ll1+2 
  ll4=ll1+3 
  ll5=ll1+4 
  ll6=ll1+5 
  di11=-1.d0/(pn-an)*(sa/sm/apr/apr*(.5d0*(cl(1,ll2)+cl(1,ll4))-    &
     &.25d0*(cl(2,ll2)+cl(2,ll4)))-a/sa/sm/apr*(.25d0*(cl(1,ll2)+       &
     &cl(1,ll4))-.125d0*(cl(2,ll2)+cl(2,ll4)))+sa*a/sm/apr*(.5d0*       &
     &(cl(21,ll2)+cl(21,ll4))-.25d0*(cl(24,ll2)+cl(24,ll4)))-           &
     &sa*a/sm/apr**3+.25d0*(a/apr)**2/sa/sm)                            
  di12=1.d0/(pn-2.d0*an)*(sa/sm/apr/apr/2.d0*cl(9,ll2)+             &
     &a/sa/sm/apr*(-.75d0*cl(9,ll2)+.125d0*cl(11,ll2))+sa*a/sm/apr/2.d0 &
     &*cl(28,ll2)-sa*a/2.d0/sm/apr**3+9.d0*a*a/16.d0/sa/sm/apr/apr)     
  di13=(3.d0*a/8.d0/sm/sa/apr*(-cl(9,ll4)+cl(11,ll4))-              &
     &9.d0*a*a/16.d0/sm/sa/apr/apr)/pn                                  
  di14=-1.d0/(pn+an)*(a/sm/sa/apr/8.d0*cl(17,ll4)-                  &
     &a*a/32.d0/sa/sm/apr/apr)                                          
  di15=1.d0/(pn-3.d0*an)*(3.d0*a/8.d0/sa/sm/apr*cl(17,ll2)-         &
     &9.d0*a*a/32.d0/sa/sm/apr/apr)                                     
  di16=3.d0/(2.d0*aa*aa)*(1.d0/((pn-an)*(pn-an))*(a/apr*(.5d0*      &
     &(cl(1,ll2)+cl(1,ll4))-.25d0*(cl(2,ll2)+cl(2,ll4)))-               &
     &.5d0*(a/apr)**2))                                                 
  di17=3.d0/(2.d0*aa*aa)*(-1.d0/((pn-2.d0*an)*(pn-2.d0*an))*        &
     &(a/apr*cl(9,ll2)-.5d0*(a/apr)**2))                                
  yd1=yd1+di11+di12+di13+di14+di15+di16+di17 
  di31=-1.d0/(pn-an)*(.25d0*sa/sm/apr**2*(cl(1,ll2)+cl(1,ll4)-      &
     &cl(5,ll2)-cl(5,ll4))+a/sa/sm/apr*(-.125d0*(cl(1,ll2)+cl(1,ll4))   &
     &+.125d0*(cl(5,ll2)+cl(5,ll4)))+.25d0*a*sa/sm/apr*(cl(21,ll2)+     &
     &cl(21,ll4)-cl(25,ll2)-cl(25,ll4))-.5d0*sa*a/sm/apr**3             &
     &+.125d0*(a/apr)**2/sa/sm)                                         
  di32=-1.d0/(pn-2.d0*an)*(a/16.d0/sa/sm/apr*(cl(9,ll2)-            &
     &cl(12,ll2))-(a/apr)**2/16.d0/sa/sm)                               
  di33=(3.d0*a/16.d0/sa/sm/apr*(-cl(9,ll4)+cl(12,ll4))-             &
     &9.d0*(a/apr)**2/16.d0/sa/sm)/pn                                   
  di34=-1.d0/(pn+an)*(.25d0*a/sa/sm/apr*cl(20,ll4)                  &
     &-.125d0*(a/apr)**2/sa/sm)                                         
  di35=3.d0/(2.d0*aa*aa)*(1.d0/((pn-an)*(pn-an))*(a/apr/4.d0*       &
     &(cl(1,ll2)+cl(1,ll4)-cl(5,ll2)-cl(5,ll4))-.25d0*(a/apr)**2))      
  yd3=yd3+di31+di32+di33+di34+di35 
  yd5=-2.d0*yd3 
  di61=1.d0/(pn-an)*(sa/sm/apr**2*(cl(1,ll1)+cl(1,ll5)+             &
     &.25d0*(cl(6,ll1)+cl(6,ll3)))+a/sa/sm/apr*(-.25d0*(cl(1,ll1)       &
     &+cl(1,ll5))+.125d0*(cl(2,ll1)+cl(2,ll5))-3.d0/16.d0*cl(6,ll3)     &
     &+cl(6,ll1)/16.d0)+a*sa/sm/apr*(cl(21,ll1)+cl(21,ll5)+.25d0*       &
     &(cl(26,ll3)+cl(26,ll1))))                                         
  di62=1.d0/(2.d0*pn-an)*(sa/sm/apr**2*cl(9,ll5)+a/sa/sm/apr/4.d0*  &
     &cl(9,ll5)+sa*a/sm/apr*cl(28,ll5))                                 
  di63=1.d0/(2.d0*pn-3.d0*an)*.75d0*a/sa/sm/apr*cl(9,ll1) 
  di64=1.d0/(pn-2.d0*an)*(.5d0*sa/sm/apr**2*cl(13,ll1)+             &
     &.5d0*sa*a/sm/apr*cl(30,ll1)+a/sm/sa/apr*(-.125d0*cl(13,ll1)       &
     &+cl(15,ll1)/16.d0+.125d0*cl(16,ll3)))                             
  di65=(a/sa/sm/apr*(3.d0*cl(15,ll3)/16.d0+3.d0*cl(16,ll1)/8.d0     &
     &+.75d0*cl(17,ll5)))/pn                                            
  di66=-1.d0/(pn+an)*a/16.d0/sa/sm/apr*cl(19,ll3) 
  di67=1.d0/(pn-3.d0*an)*3.d0*a/16.d0/sa/sm/apr*cl(19,ll1) 
  di68=3.d0/(2.d0*aa*aa)*(-1.d0/((pn-an)*(pn-an))*a/apr*            &
     &(cl(1,ll1)+cl(1,ll5)+.25d0*(cl(6,ll3)+cl(6,ll1))))                
  di69=3.d0/(2.d0*aa*aa)*(-1.d0/((2.d0*pn-an)*(2.d0*pn-an))*        &
     &a/apr*cl(9,ll5))                                                  
  di610=3.d0/(2.d0*aa*aa)*(-1.d0/((pn-2.d0*an)*(pn-2.d0*an))        &
     &*a/apr*cl(13,ll1))                                                
  yd6=yd6+di61+di62+di63+di64+di65+di66+di67+di68+di69+di610 
  ya(1)=ya(1)+yd1*pm 
  ya(2)=ya(2)+yd3*pm 
  ya(3)=ya(3)+yd5*pm 
  ya(4)=ya(4)+yd6*pm 
! Corrections of the second order and fourth degree (only 2:1, 3:1)     
  yy(1)=1/((3.d0*pn-an)*(3.d0*pn-an))*cl(17,ll6)*cl(17,ll6)/32.d0 
  yy(2)=1/((3.d0*pn-an)*(3.d0*pn-an))*cl(19,ll5)*cl(19,ll5)/32.d0 
  yy(3)=1/((3.d0*pn-an)*(3.d0*pn-an))*cl(20,ll6)*cl(20,ll6)/32.d0 
  yy(4)=4.d0*yy(3) 
  yy(5)=0.d0 
  yy(6)=yy(5) 
  yy(7)=0.d0 
  yy(8)=-yy(4) 
  yy(9)=yy(8) 
  yy(10)=-2.d0*yy(5) 
  yy(11)=-2.d0*yy(7) 
  yy(12)=1/((3.d0*pn-an)*(3.d0*pn-an))*cl(17,ll6)*cl(19,ll5)/16.d0 
  yy(13)=1/((3.d0*pn-an)*(3.d0*pn                                   &
     &-an))*(cl(19,ll5)*cl(38,ll4)/16.d0-a/apr*27.d0*cl(19,ll5)/32.d0)  
  yy(14)=0.d0 
  yy(15)=yy(14) 
  yy(16)=0.d0 
  yy(17)=0.d0 
  yy(18)=2.d0*yy(3) 
  yy(19)=1/((3.d0*pn-an)*(3.d0*pn-an))*(cl(17,ll6)*cl(38,ll4)/16.d0-&
     &a/apr*27.d0*cl(17,ll6)/32.d0)                                     
  yy(20)=1/((3.d0*pn-an)*(3.d0*pn-an))*cl(17,ll6)*cl(20,ll6)/16.d0 
  yy(21)=yy(20) 
  yy(22)=-2.d0*yy(20) 
  yy(23)=1/((3.d0*pn-an)*(3.d0*pn-an))*cl(19,ll5)*cl(20,ll6)/16.d0 
  yy(24)=yy(23) 
  yy(25)=-2.d0*yy(23) 
  yy(26)=1/((3.d0*pn-an)*(3.d0*pn-an))                              &
     &*(cl(38,ll4)*cl(20,ll6)/16.d0-a/apr*27.d0*cl(20,ll6)/32.d0)       
  yy(27)=-2.d0*yy(26) 
  DO k8=1,27 
     yb(k8)=yb(k8)+pm*3.d0/(2.d0*aa*aa)*yy(k8) 
  ENDDO
END SUBROUTINE lly
! *********************************************************             
SUBROUTINE lclong(a,np,b,c,e) 
  USE controlmod
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT(IN) :: a
  INTEGER, INTENT(IN) :: np 
  DOUBLE PRECISION, INTENT(OUT) :: b(180),c(180),e(5) 
! END INTERFACE
  DOUBLE PRECISION v(180,3), s, cj, aku, f
  INTEGER i, np0, k, j
! computation of the initial Laplace coefficients (i=0)                 
  do 30 i=1,3 
     s=.5d0+(i-1)*1.d0 
     v(1,i)=2.d0 
     aku=1.d0 
     cj=0.d0 
5    f=((s+cj)/(cj+1.d0)) 
     aku=aku*a*a*f*f 
     v(1,i)=v(1,i)+2.d0*aku 
     if(2.d0*aku.gt.excrit) then 
        cj=cj+1.d0 
        goto 5 
     endif
! interface for computation of other-than-initial Laplace coefficients  
     if(abs(s-2.5d0).lt.1d-3)then 
        np0=5 
     else 
        np0=np+2 
     endif
     DO 20 k=2,np0 
        j=k-1 
        call laplon(s,j,a,v(k,i)) 
20   ENDDO
30 ENDDO
  DO k=1,np+2 
     b(k)=v(k,1) 
     c(k)=v(k,2) 
  ENDDO
  DO  k=1,5 
     e(k)=v(k,3) 
  ENDDO   
END SUBROUTINE lclong
!  *****************************************************************    
!  computation of other-than-initial individual Laplace coefficients    
!                                                                       
SUBROUTINE laplon(s,j,alpha,v) 
  USE controlmod
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT(IN) :: alpha,s
  INTEGER, INTENT(IN) :: j
  DOUBLE PRECISION, INTENT(OUT) :: v
! END INTERFACE
  DOUBLE PRECISION z, z1, z2, count, term
  INTEGER l, l1, iexp
  v=0.d0 
  z=1.d0 
  DO 10 l1=1,j 
     l=l1-1 
     z=z*(s+l*1.d0)/(l1*1.d0) 
10 ENDDO
  z =2.d0*z*(alpha**j) 
  z1 = 1.d0 
  z2 = 1.d0 
  count = 0.d0 
  iexp = 2 
  term = 1.d0 
20 v = v + z*term 
  z1 = z1*((s + count)/(count + 1.d0)) 
  z2 = z2*((s + j*1.d0+count)/((j + 1)*1.d0+count)) 
  term = z1*z2*(alpha**iexp) 
  if (z*term.le.excrit) goto 30 
  count = count + 1.d0 
  iexp = iexp + 2 
  goto 20 
30 continue 
END SUBROUTINE laplon
! *********************************************************             
!  computation of derivatives                                           
!                                                                       
SUBROUTINE lcderlon(a,np,b1,b2,b3,b4,c1,c2) 
  USE controlmod
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT(IN) :: a ! ratio of semimajor axes
  INTEGER, INTENT(IN) :: np ! no terms 
  DOUBLE PRECISION, INTENT(OUT) :: b1(180),b2(180),b3(180),b4(5),c1(180),c2(5)
! END INTERFACE
  DOUBLE PRECISION dv(180,6), s, aku, cj, f
  INTEGER i, nn, j, k, ider, np1
!     first derivatives                                                 
!     computation of the initial first derivatives                      
  do 996 i=1,2 
     s=5.d-1+(i-1)*1.d0 
     dv(1,i)=0.d0 
     aku=1.d0 
     cj=0.d0 
     nn=2 
1    f=((s+cj)/(cj+1.d0)) 
     aku=aku*a*a*f*f 
     dv(1,i)=dv(1,i)+2.d0*aku*nn 
     if(2.d0*aku*nn.gt.excrit) then 
        cj=cj+1.d0 
        nn=nn+2 
        goto 1 
     endif
!                                                                       
! inteface for computation of other first derivatives                   
     ider=1 
     DO 10 k=2,np+2 
        j=k-1 
        call derivlon(s,j,dv(k,i),a,ider) 
10   ENDDO
996 ENDDO
  DO 19 k=1,np+2 
     b1(k)=dv(k,1) 
     c1(k)=dv(k,2) 
19 ENDDO
!                                                                       
! second derivatives                                                    
!                                                                       
!     computation of the initial second derivatives                     
  DO 997 i=3,4 
     s=5.d-1+(i-3)*1.d0 
     dv(1,i)=0.d0 
     aku=1.d0 
     cj=0.d0 
     nn=2 
2    f=((s+cj)/(cj+1.d0)) 
     aku=aku*a*a*f*f 
     dv(1,i)=dv(1,i)+2.d0*aku*nn*(nn-1) 
     if(2.d0*aku*nn*(nn-1).gt.excrit) then 
        cj=cj+1.d0 
        nn=nn+2 
        goto 2 
     endif
!                                                                       
! inteface for computation of other second derivatives                  
     ider=2 
     if(s.eq.1.5d0) then 
        np1=5 
     else 
        np1=np+2 
     endif
     DO 11 k=2,np1 
        j=k-1 
        call derivlon(s,j,dv(k,i),a,ider) 
11   ENDDO
997 ENDDO
  DO 20 k=1,np+2 
     b2(k)=dv(k,3) 
20 ENDDO
  DO 210 k=1,5 
     c2(k)=dv(k,4) 
210 ENDDO
!                                                                       
! third derivatives                                                     
!                                                                       
!     computation of the initial third derivative                       
  s=5.d-1 
  dv(1,5)=0.d0 
  aku=s*s*a*a 
  cj=1.d0 
  nn=4 
3 f=((s+cj)/(cj+1.d0)) 
  aku=aku*a*a*f*f 
  dv(1,5)=dv(1,5)+2.d0*aku*nn*(nn-1)*(nn-2) 
  if(2.d0*aku*nn*(nn-1)*(nn-2).gt.excrit) then 
     cj=cj+1.d0 
     nn=nn+2 
     goto 3 
  endif
!                                                                       
! inteface for computation of other third derivatives                   
  ider=3 
  DO 12 k=2,np+2 
     j=k-1 
     call derivlon(s,j,dv(k,5),a,ider) 
12 ENDDO
998 CONTINUE 
  DO 21 k=1,np+2 
     b3(k)=dv(k,5) 
21 ENDDO
!                                                                       
! fourth derivative                                                     
!                                                                       
!     computation of the initial fourth derivative                      
  s=5.d-1 
  dv(1,6)=0.d0 
  aku=s*s*a*a 
  cj=1.d0 
  nn=4 
4 f=((s+cj)/(cj+1.d0)) 
  aku=aku*a*a*f*f 
  dv(1,6)=dv(1,6)+2.d0*aku*nn*(nn-1)*(nn-2)*(nn-3) 
  if(2.d0*aku*nn*(nn-1)*(nn-2)*(nn-3).gt.excrit) then 
     cj=cj+1.d0 
     nn=nn+2 
     goto 4 
  endif
!                                                                       
! inteface for computation of other fourth derivatives                  
  ider=4 
  DO 13 k=2,5 
     j=k-1 
     call derivlon(s,j,dv(k,6),a,ider) 
13 ENDDO
  Do 22 k=1,5 
     b4(k)=dv(k,6) 
22 ENDDO
!                                                                       
END SUBROUTINE lcderlon
! *********************************************************             
!  computation of individual derivatives                                
!                                                                       
SUBROUTINE derivlon(s,j,dv,alpha,ider)
  USE controlmod, ONLY: excrit
  IMPLICIT NONE 
  DOUBLE PRECISION, INTENT(IN) :: s, alpha
  INTEGER, INTENT(IN) :: ider,j 
  DOUBLE PRECISION, INTENT(OUT) :: dv
! END INTERFACE
  DOUBLE PRECISION z, z1, z2, count, term
  INTEGER l1, l, iexp
  dv  = 0.d0 
  z = 1.d0 
  DO 10 l1=1,j 
     l=l1-1 
     z = z*(s + l *1.d0)/(l1*1.d0) 
10 ENDDO
  z=2.d0*z 
  if(ider.eq.1) then 
     z1 = 1.d0 
     z2 = 1.d0 
     count = 0.d0 
     iexp = 0 
     term = alpha**j 
20   dv = dv + z*(j+iexp)*term 
     z1 = z1*((s + count)/(count + 1.d0)) 
     z2 = z2*((s + j*1.d0+count)/((j + 1)*1.d0+count)) 
     iexp = iexp + 2 
     term = z1*z2*(alpha**(j+iexp)) 
     if (z*(j+iexp)*term.lt.excrit) goto 30 
     count = count + 1.d0 
     goto 20 
30   continue 
  elseif(ider.eq.2) then 
     z1 = 1.d0 
     z2 = 1.d0 
     count = 0.d0 
     iexp = 0 
     term = alpha**j 
40   dv = dv + z*(j+iexp)*(j+iexp-1)*term 
     z1 = z1*((s + count)/(count + 1.d0)) 
     z2 = z2*((s + j*1.d0+count)/((j + 1)*1.d0+count)) 
     iexp=iexp+2 
     term = z1*z2*(alpha**(j+iexp)) 
     if (z*(j+iexp)*(j+iexp-1)*term.lt.excrit) goto 50 
     count = count + 1.d0 
     goto 40 
50   continue 
  elseif(ider.eq.3) then 
     z1 = 1.d0 
     z2 = 1.d0 
     count = 0.d0 
     iexp = 0 
     term = alpha**j 
60   dv = dv + z*(j+iexp)*(j+iexp-1)*(j+iexp-2)*term 
     z1 = z1*((s + count)/(count + 1.d0)) 
     z2 = z2*((s + j*1.d0+count)/((j + 1)*1.d0+count)) 
     iexp=iexp+2 
     term = z1*z2*(alpha**(j+iexp)) 
     if (z*(j+iexp)*(j+iexp-1)*(j+iexp-2)*term.lt.excrit) goto 70 
     count = count + 1.d0 
     goto 60 
70   continue 
  elseif(ider.eq.4) then 
     z1 = 1.d0 
     z2 = 1.d0 
     count = 0.d0 
     iexp = 0 
     term = alpha**j 
80   dv = dv + z*(j+iexp)*(j+iexp-1)*(j+iexp-2)*(j+iexp-3)*term 
     z1 = z1*((s + count)/(count + 1.d0)) 
     z2 = z2*((s + j*1.d0+count)/((j + 1)*1.d0+count)) 
     iexp=iexp+2 
     term = z1*z2*(alpha**(j+iexp)) 
     if((j.eq.1).and.(iexp.eq.2)) goto 85 
     if (z*(j+iexp)*(j+iexp-1)*(j+iexp-2)*(j+iexp-3)*term.lt.excrit)   &
          &goto 90                                                           
85   count = count + 1.d0 
     goto 80 
90   continue 
  endif
END SUBROUTINE derivlon

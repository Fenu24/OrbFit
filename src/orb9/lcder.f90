! *********************************************************             
!  computation of derivatives                                           
!                                                                       
SUBROUTINE lcder(a,np,excrit,b1,b2,b3,b4,b5,c1,c2,c3,e1) 
  IMPLICIT NONE 
  DOUBLE PRECISION, INTENT(IN) :: a ! semimajor axis
  INTEGER, INTENT(IN) :: np ! no. planets
  DOUBLE PRECISION, INTENT(IN) :: excrit ! control on convergence of coefficient
  DOUBLE PRECISION, INTENT(OUT) :: b1(180),b2(180),b3(180),b4(180),b5(180), &
     &     c1(180),c2(180),c3(180),e1(180) 
! END INTERFACE
  DOUBLE PRECISION dv(180,9)
! scalar temporaries
  DOUBLE PRECISION s, aku, cj,f
! loop indexes
  INTEGER i, nn, ider, k , j
!                                                                       
!     first derivatives                                                 
!                                                                       
!     computation of the initial first derivatives                      
  do 996 i=1,3 
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
     DO k=2,np 
        j=k-1 
        call deriv(s,j,dv(k,i),a,excrit,ider) 
     ENDDO
996 ENDDO
  DO k=1,np 
     b1(k)=dv(k,1) 
     c1(k)=dv(k,2) 
     e1(k)=dv(k,3) 
  ENDDO
!                                                                       
! second derivatives                                                    
!                                                                       
!     computation of the initial second derivatives                     
  DO 997 i=4,5 
     s=5.d-1+(i-4)*1.d0 
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
     DO k=2,np 
        j=k-1 
        call deriv(s,j,dv(k,i),a,excrit,ider) 
     ENDDO
997 ENDDO
  DO k=1,np 
     b2(k)=dv(k,4) 
     c2(k)=dv(k,5) 
  ENDDO
!                                                                       
! third derivatives                                                     
!                                                                       
!     computation of the initial third derivatives                      
  DO 998 i=6,7 
     s=5.d-1+(i-6)*1.d0 
     dv(1,i)=0.d0 
     aku=s*s*a*a 
     cj=1.d0 
     nn=4 
3    f=((s+cj)/(cj+1.d0)) 
     aku=aku*a*a*f*f 
     dv(1,i)=dv(1,i)+2.d0*aku*nn*(nn-1)*(nn-2) 
     if(2.d0*aku*nn*(nn-1)*(nn-2).gt.excrit) then 
        cj=cj+1.d0 
        nn=nn+2 
        goto 3 
     endif
!                                                                       
! inteface for computation of other third derivatives                   
     ider=3 
     DO k=2,np 
        j=k-1 
        call deriv(s,j,dv(k,i),a,excrit,ider) 
     ENDDO
998 ENDDO
  DO k=1,np 
     b3(k)=dv(k,6) 
     c3(k)=dv(k,7) 
  ENDDO
!                                                                       
! fourth and fifth derivative                                           
!                                                                       
!     computation of the initial fourth derivative                      
  s=5.d-1 
  dv(1,8)=0.d0 
  aku=s*s*a*a 
  cj=1.d0 
  nn=4 
4 f=((s+cj)/(cj+1.d0)) 
  aku=aku*a*a*f*f 
  dv(1,8)=dv(1,8)+2.d0*aku*nn*(nn-1)*(nn-2)*(nn-3) 
  if(2.d0*aku*nn*(nn-1)*(nn-2)*(nn-3).gt.excrit) then 
     cj=cj+1.d0 
     nn=nn+2 
     goto 4 
  endif
!                                                                       
! inteface for computation of other fourth derivatives                  
  ider=4 
  DO k=2,np 
     j=k-1 
     call deriv(s,j,dv(k,8),a,excrit,ider) 
  ENDDO 
  DO k=1,np 
     b4(k)=dv(k,8) 
  ENDDO 
!                                                                       
!     computation of the initial fifth derivative                       
  dv(1,9)=0.d0 
  aku=s*s*(s+1.d0)/2.d0*(s+1.d0)/2.d0*a*a*a*a 
  cj=2.d0 
  nn=6 
5 f=((s+cj)/(cj+1.d0)) 
  aku=aku*a*a*f*f 
  dv(1,9)=dv(1,9)+2.d0*aku*nn*(nn-1)*(nn-2)*(nn-3)*(nn-4) 
  if(2.d0*aku*nn*(nn-1)*(nn-2)*(nn-3)*(nn-4).gt.excrit) then 
     cj=cj+1.d0 
     nn=nn+2 
     goto 5 
  endif
!                                                                       
! inteface for computation of other fifth derivatives                   
  ider=5 
  DO k=2,np 
     j=k-1 
     call deriv(s,j,dv(k,9),a,excrit,ider) 
  ENDDO 
  DO k=1,np 
     b5(k)=dv(k,9) 
  ENDDO 
END SUBROUTINE lcder

! *********************************************************             
!  computation of individual derivatives                                
!                                                                       
SUBROUTINE deriv(s,j,dv,alpha,excrit,ider) 
  IMPLICIT NONE
! INPUT  
  DOUBLE PRECISION, INTENT(IN) :: s, alpha
  INTEGER, INTENT(IN) :: j, ider  
  DOUBLE PRECISION, INTENT(IN) :: excrit ! convergence control
! OUTPUT
  DOUBLE PRECISION, INTENT(OUT) :: dv
! END INTERFACE  
  DOUBLE PRECISION z, z1, z2, term, count
  INTEGER l, l1, iexp
  dv  = 0.d0 
  z = 1.d0 
  DO l1=1,j 
     l=l1-1 
     z = z*(s + l *1.d0)/(l1*1.d0) 
  ENDDO
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
     if (z*(j+iexp)*(j+iexp-1)*(j+iexp-2)*(j+iexp-3)*term.lt.excrit) goto 90
85   count = count + 1.d0 
     goto 80 
90   continue 
  elseif(ider.eq.5) then 
     z1 = 1.d0 
     z2 = 1.d0 
     count = 0.d0 
     iexp = 0 
     term = alpha**j 
100  dv = dv + z*(j+iexp)*(j+iexp-1)*(j+iexp-2)*(j+iexp-3)*(j+iexp-4)*term
     z1 = z1*((s + count)/(count + 1.d0)) 
     z2 = z2*((s + j*1.d0+count)/((j + 1)*1.d0+count)) 
     iexp=iexp+2 
     term = z1*z2*(alpha**(j+iexp)) 
     if((j.le.2).and.(iexp.eq.2)) goto 105 
     if (z*(j+iexp)*(j+iexp-1)*(j+iexp-2)*(j+iexp-3)*(j+iexp-4)        &
     &    *term.lt.excrit) goto 110                                     
105  count = count + 1.d0 
     goto 100 
110  continue 
  endif
END SUBROUTINE deriv

!                                                                       
! *******************************************************************   
!  interface subroutine to compute Laplace coefficients and their       
!  derivatives                                                          
SUBROUTINE lapco (a,np,excrit,b,b1,b2,b3,b4,b5,c,c1,c2,c3,e,e1) 
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT(IN) :: a, excrit ! semimajor axis ratio, control
  INTEGER, INTENT(IN) :: np ! number of Laplace coefficients is 2*np+1
  DOUBLE PRECISION, INTENT(OUT) :: b(180),c(180),e(180)
  DOUBLE PRECISION, INTENT(OUT) :: b1(180),b2(180),b3(180),b4(180),b5(180)
  DOUBLE PRECISION, INTENT(OUT) :: c1(180),c2(180),c3(180),e1(180)
!  computation of Laplace coefficients                                  
  call lc(a,np,excrit,b,c,e) 
!  computation of derivatives of Laplace coefficients                   
  call lcder(a,np,excrit,b1,b2,b3,b4,b5,c1,c2,c3,e1) 
END SUBROUTINE lapco
! *********************************************************             
SUBROUTINE lc(a,np,excrit,b,c,e) 
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT(IN) :: a, excrit ! semimajor axis ratio, control
  INTEGER, INTENT(IN) :: np ! number of Laplace coefficients is 2*np+1
  DOUBLE PRECISION, INTENT(OUT) :: b(180),c(180),e(180)
! END INTERFACE 
  DOUBLE PRECISION  v(180,3), s, aku, cj, f
  INTEGER i, j, k 
! computation of the initial Laplace coefficients (i=0)                 
  DO 30 i=1,3 
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
     DO k=2,np 
        j=k-1 
        call lap(s,j,a,excrit,v(k,i)) 
     ENDDO
30 ENDDO
  DO k=1,np 
     b(k)=v(k,1) 
     c(k)=v(k,2) 
     e(k)=v(k,3) 
  ENDDO
END SUBROUTINE lc
!  *****************************************************************    
!  computation of other-than-initial individual Laplace coefficients    
!                                                                       
SUBROUTINE lap(s,j,alpha,excrit,v) 
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT(IN) :: s, alpha, excrit
  INTEGER, INTENT(IN) :: j
  DOUBLE PRECISION, INTENT(OUT) :: v
! END INTERFACE
  DOUBLE PRECISION z,z1, z2, count, term
  INTEGER l, l1, iexp
  v=0.d0 
  z=1.d0 
  DO l1=1,j 
     l=l1-1 
     z=z*(s+l*1.d0)/(l1*1.d0) 
  ENDDO
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
END SUBROUTINE lap

!                                                                       
! *******************************************************************   
!  subroutine to compute Leverrier coefficients                         
SUBROUTINE levco(i,b,b1,b2,b3,b4,c,c1,c2,e,a,pai,cl) 
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT(IN) :: a, pai ! ratio of semimaj axis, inverse of plan. semimaj axis
!  laplace coefficients 
  DOUBLE PRECISION, INTENT(IN) ::  b(180),c(180),e(180) 
  DOUBLE PRECISION, INTENT(IN) ::  b1(180),b2(180),b3(180),b4(180),c1(180),c2(180) 
  INTEGER, INTENT(IN) :: i ! order of term 
!  leverrier coefficients                                               
  DOUBLE PRECISION, INTENT(OUT) :: cl(380) 
! END INTERFACE
  INTEGER j,j1,j2,j3,j4
  DOUBLE PRECISION di
!                                                                       
  j=iabs(i)+1 
  j1=iabs(i-1)+1 
  j2=iabs(i+1)+1 
  j3=iabs(i-2)+1 
  j4=iabs(i+2)+1 
  di=i*1.d0 
  cl(1)=pai*b(j)/2.d0 
  cl(2)=pai*(-2.d0*di*di*b(j)+b1(j)+b2(j)/2.d0) 
  cl(3)=cl(2) 
  cl(4)=pai*((-9.d0*di*di+16.d0*di*di*di*di)*b(j)/8.d0-             &
     &di*di*(b1(j)+b2(j))+b3(j)/2.d0+b4(j)/8.d0)                        
  cl(5)=pai*(8.d0*di*di*di*di*b(j)+(2.d0-8.d0*di*di)*b1(j)          &
     &+(7.d0-4.d0*di*di)*b2(j)+4.d0*b3(j)+b4(j)/2.d0)                   
  cl(6)=pai*((-17.d0*di*di+16.d0*di*di*di*di)*b(j)/8.d0+            &
     &3.d0*(1.d0-di*di)*b1(j)+(9.d0-2.d0*di*di)*b2(j)/2.d0              &
     &+1.5d0*b3(j)+b4(j)/8.d0)                                          
  cl(11)=-a*pai/4.d0*(c(j1)+c(j2)) 
  cl(12)=a*pai*((di*di-.5d0)*(c(j1)+c(j2))-(c1(j1)+c1(j2))          &
     &-.25d0*(c2(j1)+c2(j2)))                                           
  cl(13)=cl(12) 
  cl(17)=3.d0*a*a*pai/16.d0*(e(j3)+4.d0*e(j)+e(j4)) 
  cl(21)=pai*((2.d0*di+4.d0*di*di)*b(j)-2.d0*b1(j)-b2(j)) 
  cl(22)=pai*((-di-7.d0*di*di-14.d0*di*di*di-8.d0*di*di*di*di)      &
     &*b(j)+(3.d0+9.d0*di+6.d0*di*di)*b1(j)+.5d0*(-3.d0+9.d0*di         &
     &+8.d0*di*di)*b2(j)-3.d0*b3(j)-b4(j)/2.d0)                         
  cl(23)=pai*((-di-5.d0*di*di-14.d0*di*di*di-8.d0*di*di*di*di)      &
     &*b(j)+(-3.d0+9.d0*di+10.d0*di*di)*b1(j)+.5d0*(-21.d0+9.d0         &
     &*di+8.d0*di*di)*b2(j)-5.d0*b3(j)-b4(j)/2.d0)                      
  cl(27)=a*pai*((1.d0-di-2.d0*di*di)*(c(j1)+c(j2))+2.d0*            &
     &(c1(j1)+c1(j2))+.5d0*(c2(j1)+c2(j2)))                             
  cl(31)=.25d0*pai*((20.d0*di+61.d0*di*di+56.d0*di*di*di            &
     &+16.d0*di*di*di*di)*b(j)-(20.d0+36.d0*di+16.d0*di*di)*b1(j)+      &
     &(2.d0-18.d0*di-8.d0*di*di)*b2(j)+8.d0*b3(j)+b4(j))                
  cl(36)=.25d0*a*pai*((12.d0-15.d0*di+4.d0*di*di)*c(j1)+            &
     &(8.d0-4.d0*di)*c1(j1)+c2(j1))                                     
  cl(40)=a*pai*((-5.d0+7.d0*di-2.d0*di*di)*c(j1)+(-4.d0             &
     &+2.d0*di)*c1(j1)-c2(j1)/2.d0)                                     
  cl(44)=.25d0*a*pai*((10.d0-13.d0*di+4.d0*di*di)*c(j1)+            &
     &(8.d0-4.d0*di)*c1(j1)+c2(j1))                                     
  cl(50)=pai*(-2.d0*di*b(j)-b1(j)) 
  cl(51)=pai*((di-5.d0*di*di+4.d0*di*di*di)*b(j)+.5d0*              &
     &(3.d0-7.d0*di+4.d0*di*di)*b1(j)-(1.d0+di)*b2(j)-b3(j)/2.d0)       
  cl(52)=pai*(8.d0*di*di*di*b(j)-(2.d0+4.d0*di-4.d0*di*di)          &
     &*b1(j)-(4.d0+2.d0*di)*b2(j)-b3(j))                                
  cl(60)=a*pai*((.5d0+di)*(c(j1)+c(j2))+.5d0*(c1(j1)+c1(j2))) 
  cl(70)=pai*((1.d0+2.d0*di)*b(j)+b1(j)) 
  cl(71)=pai*((-4.d0*di*di-8.d0*di*di*di)*b(j)+4.d0*(1.d0+          &
     &di-di*di)*b1(j)+(5.d0+2.d0*di)*b2(j)+b3(j))                       
  cl(72)=.5d0*pai*((-1.d0-5.d0*di-14.d0*di*di-8.d0*di*di*di)*       &
     &b(j)+(7.d0+di-4.d0*di*di)*b1(j)+(7.d0+2.d0*di)*b2(j)+b3(j))       
  cl(80)=.5d0*a*pai*((-2.d0-2.d0*di)*(c(j1)+c(j2))-(c1(j1)+         &
     &c1(j2)))                                                          
  cl(90)=.5d0*pai*((5.d0*di+14.d0*di*di+8.d0*di*di*di)*b(j)-        &
     &(4.d0+7.d0*di+4.d0*di*di)*b1(j)+(1.d0-2.d0*di)*b2(j)+b3(j))       
  cl(100)=.5d0*pai*((8.d0*di+18.d0*di*di+8.d0*di*di*di)*b(j)-       &
     &(10.d0+di-4.d0*di*di)*b1(j)-(8.d0+2.d0*di)*b2(j)-b3(j))           
  cl(120)=.5d0*a*pai*((2.d0*di-5.d0)*c(j1)-c1(j1)) 
  cl(130)=.5d0*a*pai*(2.d0*(1.d0-di)*c(j1)+c1(j1)) 
  cl(172)=.5d0*pai*((-5.d0*di+4.d0*di*di)*b(j)+(-2.d0+4.d0*di)      &
     &*b1(j)+b2(j))                                                     
  cl(173)=pai/3.d0*((11.d0*di-32.d0*di*di+30.d0*di*di*di-           &
     &8.d0*di*di*di*di)*b(j)+(8.d0-23.d0*di+24.d0*di*di-8.d0*di*        &
     &di*di)*b1(j)+1.5d0*(-4.d0+3.d0*di)*b2(j)+2.d0*di*b3(j)+           &
     &b4(j)/2.d0)                                                       
  cl(174)=pai*((10.d0*di*di*di-8.d0*di*di*di*di)*b(j)-(2.d0+        &
     &di-8.d0*di*di+8.d0*di*di*di)*b1(j)+.5d0*(-2.d0+11.d0*di)*         &
     &b2(j)+(2.d0+2.d0*di)*b3(j)+b4(j)/2.d0)                            
  cl(178)=.25d0*a*pai*((2.d0+di-4.d0*di*di)*(c(j1)+c(j2))-          &
     &4.d0*di*(c1(j1)+c1(j2))-(c2(j1)+c2(j2)))                          
  cl(182)=pai*((-2.d0*di-4.d0*di*di)*b(j)-(2.d0+4.d0*di)*           &
     &b1(j)-b2(j))                                                      
  cl(183)=pai*((di-3.d0*di*di-6.d0*di*di*di+8.d0*di*di*di*di)       &
     &*b(j)+(3.d0-3.d0*di-8.d0*di*di+8.d0*di*di*di)*b1(j)-.5d0*         &
     &(3.d0+17.d0*di)*b2(j)-(3.d0+2.d0*di)*b3(j)-b4(j)/2.d0)            
  cl(184)=pai*((di+5.d0*di*di+14.d0*di*di*di+8.d0*di*di*di*di)      &
     &*b(j)-(3.d0+5.d0*di-8.d0*di*di-8.d0*di*di*di)*b1(j)-.5d0*         &
     &(21.d0+19.d0*di)*b2(j)-(5.d0+2.d0*di)*b3(j)-b4(j)/2.d0)           
  cl(188)=a*pai*((1.d0+3.d0*di+2.d0*di*di)*(c(j1)+c(j2))+           &
     &(2.d0+2.d0*di)*(c1(j1)+c1(j2))+.5d0*(c2(j1)+c2(j2)))              
  cl(192)=.5d0*pai*((4.d0+9.d0*di+4.d0*di*di)*b(j)+(6.d0+           &
     &4.d0*di)*b1(j)+b2(j))                                             
  cl(193)=pai*((-8.d0*di*di-18.d0*di*di*di-8.d0*di*di*di*di)*       &
     &b(j)+(10.d0+13.d0*di-8.d0*di*di-8.d0*di*di*di)*b1(j)+.5d0*        &
     &(34.d0+25.d0*di)*b2(j)+(6.d0+2.d0*di)*b3(j)+b4(j)/2.d0)           
  cl(194)=pai/3.d0*((-8.d0-31.d0*di-56.d0*di*di-38.d0*di*di*di      &
     &-8.d0*di*di*di*di)*b(j)+(16.d0-5.d0*di-24.d0*di*di-8.d0*di*       &
     &di*di)*b1(j)+1.5d0*(20.d0+9.d0*di)*b2(j)+(8.d0+2.d0*di)*          &
     &b3(j)+b4(j)/2.d0)                                                 
  cl(198)=-.25d0*a*pai*((10.d0+13.d0*di+4.d0*di*di)*(c(j1)+         &
     &c(j2))+(8.d0+4.d0*di)*(c1(j1)+c1(j2))+c2(j1)+c2(j2))              
  cl(202)=pai/3.d0*((13.d0*di+41.d0*di*di+34.d0*di*di*di+8.d0*      &
     &di*di*di*di)*b(j)-(9.d0+23.d0*di+24.d0*di*di+8.d0*di*di*di)       &
     &*b1(j)+1.5d0*(3.d0+di)*b2(j)+(1.d0+2.d0*di)*b3(j)-b4(j)/2.d0)     
  cl(206)=pai/3.d0*((27.d0*di+65.d0*di*di+42.d0*di*di*di+8.d0*di    &
     &*di*di*di)*b(j)-(39.d0+7.d0*di-24.d0*di*di-8.d0*di*di*di)*        &
     &b1(j)-1.5d0*(27.d0+11.d0*di)*b2(j)-(9.d0+2.d0*di)*b3(j)-          &
     &b4(j)/2.d0)                                                       
  cl(212)=.5d0*a*pai*c(j1) 
  cl(213)=a*pai*((-7.d0+8.d0*di-2.d0*di*di)*c(j1)+2.d0*             &
     &c1(j1)+c2(j1)/2.d0)                                               
  cl(214)=a*pai*((1.d0-2.d0*di*di)*c(j1)+2.d0*c1(j1)+c2(j1)/2.d0) 
  cl(218)=-3.d0*a*a*pai/4.d0*(e(j3)+e(j)) 
  cl(222)=a*pai*((-5.d0-3.d0*di+2.d0*di*di)*c(j1)-4.d0*c1(j1)-      &
     &c2(j1)/2.d0)                                                      
  cl(226)=a*pai*((3.d0-5.d0*di+2.d0*di*di)*c(j1)-c2(j1)/2.d0) 
  cl(240)=pai/3.d0*((-13.d0*di+15.d0*di*di-4.d0*di*di*              &
     &di)*b(j)-1.5d0*(3.d0-9.d0*di+4.d0*di*di)*b1(j)+                   &
     &(3.d0-3.d0*di)*b2(j)-b3(j)/2.d0)                                  
  cl(250)=pai/2.d0*((-5.d0*di-6.d0*di*di+8.d0*di*di*di)*            &
     &b(j)-(4.d0+di-12.d0*di*di)*b1(j)+(1.d0+6.d0*di)*b2(j)+b3(j))      
  cl(260)=pai*((-4.d0*di-9.d0*di*di-4.d0*di*di*di)*b(j)             &
     &-0.5d0*(10.d0+25.d0*di+12.d0*di*di)*b1(j)-(4.d0+3.d0*             &
     &di)*b2(j)-b3(j)/2.d0)                                             
  cl(270)=pai/2.d0*((27.d0+65.d0*di+42.d0*di*di+8.d0*               &
     &di*di*di)*b(j)/3.d0+(17.d0+17.d0*di+4.d0*di*di)*b1(j)             &
     &+(5.d0+2.d0*di)*b2(j)+b3(j)/3.d0)                                 
  cl(290)=a*pai/2.d0*((3.d0-2.d0*di)*c(j1)-c1(j1)) 
  cl(300)=a*pai*((1.d0+di)*c(j1)+c1(j1)/2.d0) 
  cl(336)=pai/6.d0*((-206.d0*di+283.d0*di*di-120.d0*                &
     &di*di*di+16.d0*di*di*di*di)*b(j)/4.d0-(16.d0-59.d0*di+42.d0*      &
     &di*di-8.d0*di*di*di)*b1(j)+1.5d0*(8.d0-13.d0*di+4.d0*             &
     &di*di)*b2(j)-(3.d0-2.d0*di)*b3(j)+b4(j)/4.d0)                     
  cl(340)=pai/3.d0*((-13.d0*di-11.d0*di*di+26.d0*di*di*di           &
     &-8.d0*di*di*di*di)*b(j)-(9.d0-5.d0*di-30.d0*di*di+16.d0*          &
     &di*di*di)*b1(j)+1.5d0*(3.d0+7.d0*di-8.d0*di*di)*b2(j)+            &
     &(1.d0-4.d0*di)*b3(j)-b4(j)/2.d0)                                  
  cl(344)=pai/4.d0*((-20.d0*di-29.d0*di*di+16.d0*di*di*di           &
     &+16.d0*di*di*di*di)*b(j)-(20.d0+16.d0*di-48.d0*di*di-32.d0        &
     &*di*di*di)*b1(j)+(2.d0+36.d0*di+24.d0*di*di)*b2(j)+(8.d0          &
     &+8.d0*di)*b3(j)+b4(j))                                            
  cl(348)=pai/3.d0*((-27.d0*di-65.d0*di*di-42.d0*di*di*di           &
     &-8.d0*di*di*di*di)*b(j)-(39.d0+109.d0*di+78.d0*di*di+16.d0        &
     &*di*di*di)*b1(j)-1.5d0*(27.d0+31.d0*di+8.d0*di*di)*b2(j)          &
     &-(9.d0+4.d0*di)*b3(j)-b4(j)/2.d0)                                 
  cl(352)=pai/6.d0*((256.d0+646.d0*di+499.d0*di*di+152.d0*          &
     &di*di*di+16.d0*di*di*di*di)*b(j)/4.d0+(142.d0+173.d0*di+66.d0*    &
     &di*di+8.d0*di*di*di)*b1(j)+1.5d0*(38.d0+25.d0*di+4.d0*            &
     &di*di)*b2(j)+(7.d0+2.d0*di)*b3(j)+b4(j)/4.d0)                     
  cl(358)=a*pai/4.d0*((16.d0-17.d0*di+4.d0*di*di)*c(j1)-            &
     &(8.d0-4.d0*di)*c1(j1)+c2(j1))                                     
  cl(362)=a*pai*((3.d0+di-2.d0*di*di)*c(j1)-2.d0*di*c1(j1)          &
     &-.5d0*c2(j1))                                                     
  cl(366)=a*pai/4.d0*((10.d0+13.d0*di+4.d0*di*di)*c(j1)+            &
     &(8.d0+4.d0*di)*c1(j1)+c2(j1))                                     
  cl(372)=3.d0*a**2*pai*e(j3)/8.d0 
END SUBROUTINE levco

!                                                                       
! ****************************************************************      
!  setting the delimiters for calculation of                            
!  Laplace, LeVerrier and Yuasa coefficients                            
!                                                                       
INTEGER FUNCTION nalfa(a) 
  DOUBLE PRECISION, INTENT(IN) :: a ! ratio of semimajor axes
  if(a.lt.0.3d0) then 
     nalfa=11 
  elseif(a.lt.0.4d0) then 
     nalfa=21 
  elseif(a.lt.0.5d0) then 
     nalfa=31 
  elseif(a.lt.0.6d0)then 
     nalfa=43 
  elseif(a.lt.0.7d0)then 
     nalfa=61 
  else 
     nalfa=81 
  endif
END FUNCTION nalfa

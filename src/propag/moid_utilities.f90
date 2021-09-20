! SUBROUTINES: rvfft, irvfft, choosedeg, code_input, decode_out

!------------------------------------------------------------c          
!  A real-valued, in-place, split-radix FFT program          c          
! Real input and output data in array x                      c          
! Length is N=2**M                                           c          
! Decimation-in-time                                         c          
! Output in order                                            c          
!   [Re(0),Re(1),...,Re(N/2),Im(N/2-1),...,Im(1)]            c          
! H.V. Sorensen, Rice University, 1985                       c          
! Modified by D. Bini 1993                                   c          
!------------------------------------------------------------c          
!------------------------------------------------------------c          
! For convolution (polynomial multiplication)                c          
! X must contain the coefficients of the first polynomial    c          
! Y the coefficients of the second polynomial                c          
! Then compute u=irvfft(x), v=irfft(y) and set               c          
!     w(0)=u(0)*v(0), w(i)=u(i)*v(i)-u(n-i+1)*v(n-i+1)       c          
!     w(n-i+1)=u(i)*v(n-i+1)+u(n-i+1)*v(i)                   c          
! Finally output rvfft(w)/n                                  c          
!------------------------------------------------------------c          
                                                                        
                                                                        
      subroutine rvfft(x,n,m)
      implicit double precision (a-h,o-z)
      dimension x(n) 
!----DIGIT REVERSE COUNTER NOT REMOVED -----------------c               
  100  j=1 
       n1=n-1 
       do 104 i=1,n1 
           if(i.ge.j)goto 101 
           xt=x(j) 
           x(j)=x(i) 
           x(i)=xt 
  101      k=n/2 
  102      if (k.ge.j)goto 103 
               j=j-k 
               k=k/2 
               goto 102 
  103      j=j+k 
  104  continue 
!----LENGTH TWO BUTTERFLIES-------------------------c                   
      is=1 
      id=4 
   70 do 60 i0=is,n,id 
        i1=i0+1 
        r1=x(i0) 
        x(i0)=r1+x(i1) 
        x(i1)=r1-x(i1) 
   60 continue 
      is=2*id-1 
      id=4*id 
      if(is.lt.n) goto 70 
!-----L SHAPED BUTTERFLIES--------------------------c                   
      n2=2 
      do 10 k=2,m 
        n2=n2*2 
        n4=n2/4 
        n8=n2/8 
        e=6.283185307179586D0/n2 
        is=0 
        id=n2*2 
   40   do 38 i=is,n-1,id 
          i1=i+1 
          i2=i1+n4 
          i3=i2+n4 
          i4=i3+n4 
          t1=x(i4)+x(i3) 
          x(i4)=x(i4)-x(i3) 
          x(i3)=x(i1)-t1 
          x(i1)=x(i1)+t1 
          if(n4.eq.1)goto 38 
          i1=i1+n8 
          i2=i2+n8 
          i3=i3+n8 
          i4=i4+n8 
          t1=(x(i3)+x(i4))/dsqrt(2.0d0) 
          t2=(x(i3)-x(i4))/dsqrt(2.0d0) 
          x(i4)=x(i2)-t1 
          x(i3)=-x(i2)-t1 
          x(i2)=x(i1)-t2 
          x(i1)=x(i1)+t2 
   38 continue 
         is=2*id-n2 
         id=4*id 
      if(is.lt.n)goto 40 
      a=e 
      do 32 j=2,n8 
         a3=3*a 
         cc1=dcos(a) 
         ss1=dsin(a) 
         cc3=dcos(a3) 
         ss3=dsin(a3) 
         a=j*e 
         is=0 
         id=2*n2 
   36    do 30 i=is,n-1,id 
            i1=i+j 
            i2=i1+n4 
            i3=i2+n4 
            i4=i3+n4 
            i5=i+n4-j+2 
            i6=i5+n4 
            i7=i6+n4 
            i8=i7+n4 
            t1=x(i3)*cc1+x(i7)*ss1 
            t2=x(i7)*cc1-x(i3)*ss1 
            t3=x(i4)*cc3+x(i8)*ss3 
            t4=x(i8)*cc3-x(i4)*ss3 
            t5=t1+t3 
            t6=t2+t4 
            t3=t1-t3 
            t4=t2-t4 
            t2=x(i6)+t6 
            x(i3)=t6-x(i6) 
            x(i8)=t2 
            t2=x(i2)-t3 
            x(i7)=-x(i2)-t3 
            x(i4)=t2 
            t1=x(i1)+t5 
            x(i6)=x(i1)-t5 
            x(i1)=t1 
            t1=x(i5)+t4 
            x(i5)=x(i5)-t4 
            x(i2)=t1 
   30    continue 
            is=2*id-n2 
            id=4*id 
         if(is.lt.n)goto 36 
   32    continue 
   10 continue 
      return 
   END SUBROUTINE rvfft

!-------------------------------------------------------c               
!  A real-valued, in-place, split-radix IFFT program    c               
! Hermitian symmetric input and real output in array X  c               
! Length is N=2**M                                      c               
! Decimation-in-frequency                               c               
! Input  order                                          c               
!     [Re(0),Re(1),...,Re(N/2),Im(N/2-1),...,Im(1)]     c               
! H.V. Sorensen, Rice University, 1985                  c               
! Modified by D. Bini 1993                              c               
!-------------------------------------------------------c               
      subroutine irvfft(x,n,m) 
      implicit double precision (a-h,o-z) 
      dimension x(n) 
!----L SHAPED BUTTERFLIES-------------------------------c               
      n2=2*n 
      do 10 k=1,m-1 
         is=0 
         id=n2 
         n2=n2/2 
         n4=n2/4 
         n8=n4/2 
        e=6.283185307179586D0/n2 
   17    do 15 i=is,n-1,id 
            i1=i+1 
            i2=i1+n4 
            i3=i2+n4 
            i4=i3+n4 
            t1=x(i1)-x(i3) 
            x(i1)=x(i1)+x(i3) 
            x(i2)=2*x(i2) 
            x(i3)=t1-2*x(i4) 
            x(i4)=t1+2*x(i4) 
            if (n4.eq.1)goto 15 
            i1=i1+n8 
            i2=i2+n8 
            i3=i3+n8 
            i4=i4+n8 
            t1=(x(i2)-x(i1))/dsqrt(2.0d0) 
            t2=(x(i4)+x(i3))/dsqrt(2.0d0) 
            x(i1)=x(i1)+x(i2) 
            x(i2)=x(i4)-x(i3) 
            x(i3)=2*(-t2-t1) 
            x(i4)=2*(-t2+t1) 
   15    continue 
            is=2*id-n2 
            id=4*id 
         if(is.lt.n-1)goto 17 
         a=e 
         do 20 j=2,n8 
            a3=3*a 
            cc1=dcos(a) 
            ss1=dsin(a) 
            cc3=dcos(a3) 
            ss3=dsin(a3) 
            a=j*e 
            is=0 
            id=2*n2 
   40       do 30 i=is,n-1,id 
               i1=i+j 
               i2=i1+n4 
               i3=i2+n4 
               i4=i3+n4 
               i5=i+n4-j+2 
               i6=i5+n4 
               i7=i6+n4 
               i8=i7+n4 
               t1=x(i1)-x(i6) 
               x(i1)=x(i1)+x(i6) 
               t2=x(i5)-x(i2) 
               x(i5)=x(i2)+x(i5) 
               t3=x(i8)+x(i3) 
               x(i6)=x(i8)-x(i3) 
               t4=x(i4)+x(i7) 
               x(i2)=x(i4)-x(i7) 
               t5=t1-t4 
               t1=t1+t4 
               t4=t2-t3 
               t2=t2+t3 
               x(i3)=t5*cc1+t4*ss1 
               x(i7)=-t4*cc1+t5*ss1 
               x(i4)=t1*cc3-t2*ss3 
               x(i8)=t2*cc3+t1*ss3 
   30      continue 
               is=2*id-n2 
               id=4*id 
           if(is.lt.n-1)goto 40 
   20    continue 
   10 continue 
!-----LENGTH TWO BUTTERFLIES--------------------------------c           
      is=1 
      id=4 
   70 do 60 i0=is,n,id 
      i1=i0+1 
      r1=x(i0) 
      x(i0)=r1+x(i1) 
      x(i1)=r1-x(i1) 
   60 continue 
      is=2*id-1 
      id=4*id 
      if(is.lt.n)goto 70 
!----DIGIT REVERSE COUNTER NOT REMOVED -----------------c               
  100  j=1 
       n1=n-1 
       do 104 i=1,n1 
           if(i.ge.j)goto 101 
           xt=x(j) 
           x(j)=x(i) 
           x(i)=xt 
  101      k=n/2 
  102      if (k.ge.j)goto 103 
               j=j-k 
               k=k/2 
               goto 102 
  103      j=j+k 
  104  continue 
                                                                        
!-----DIVISION BY N REMOVED---------------------------------c           
       do 99 i=1,n 
          x(i)=x(i)/n 
   99  continue 
       return 
      END SUBROUTINE irvfft

! ****************************************
! subroutine to obtain angles in [0,360]
! input and output in degrees
! written by G.F.Gronchi Nov.2004
! last modified 21/10/2005 GFG
! ****************************************
 SUBROUTINE choosedeg(lambda,chl)
   USE output_control
   USE fund_const
   IMPLICIT NONE
   REAL(KIND=8), INTENT(IN) :: lambda
   REAL(KIND=8), INTENT(OUT) :: chl
! ------------- end interface ------------
   REAL(KIND=8) :: x,y !sin(lambda),cos(lambda)
! ===================================================   
   x=sin(radeg*lambda)
   y=cos(radeg*lambda)
! chl in [-180,180]                                                 
   chl=degrad*atan2(x,y)
! chl in [0,360]                                                    
   IF((chl.ge.0.d0).and.(chl.le.180.d0)) THEN 
      chl = chl 
   ELSEIF((chl.lt.0.d0).and.(chl.ge.-180.d0)) THEN 
      chl = 360.d0 + chl 
   ELSE 
      WRITE(iun_log,*)'choosedeg: input angle outside its established range!'&
           & ,chl 
!      STOP
   ENDIF
!     control                                                           
   if((chl.gt.360.d0).or.(chl.lt.0.d0)) then 
      WRITE(iun_log,*)'choosedeg: output angle outside its established range!'&
           & ,chl 
!      STOP
   endif
 END SUBROUTINE choosedeg
 
! ============================================================
! ******* CODING of the n evaluations of a polynomial ********
! ********* p: p(om(1)), p(om(n)) to apply irvfft ************
! ******** written by GIOVANNI F. GRONCHI (2004) *************
! ============================================================
  SUBROUTINE code_input(N,eval,codeval)
    IMPLICIT NONE
    INTEGER :: N
! evaluation of poly p(x) in the N-th roots of unity
    COMPLEX(KIND=8),DIMENSION(N) :: eval
! -----------------------------------------------------------
! if omega is a primitive complex root of unity then we have
! Re(p(omega^s)) = Re(p(omeg^{-s}))
! Im(p(omega^s)) = -Im(p(omeg^{-s}))   for s=1,..,N/2
! -----------------------------------------------------------
    REAL(KIND=8),DIMENSION(N) :: codeval ! coded evaluations
! ------ end interface ----------------
    INTEGER :: i ! loop indexes
! ============================================================

    DO i = 1,N/2+1
       codeval(i) = DBLE(eval(i)) ! from 1 to N/2+1 take the real part
    ENDDO
! (roots of unity are taken clockwise)
    DO i = 1,N/2-1 
       codeval(N/2+1+i) = DIMAG(eval(N/2-i+1))
    ENDDO
  END SUBROUTINE code_input

!     ******************************************************************
!     OUTPUT CONVERSION FROM THE STANDARD FORM of DTF TO THE USUAL ONE  
!     ******************************************************************
!     *********** written by GIOVANNI F. GRONCHI (2001) ****************
!     ********** Department of Mathematics, UNIVERSITY of PISA *********
!     ==================================================================
      SUBROUTINE decode_out(N,dftout1,dftout2,dfteval) 
      IMPLICIT NONE 
      INTEGER N 
      DOUBLE PRECISION dftout1( * ),dftout2( * ) 
      COMPLEX*16 dfteval( * ) 
!     ----------------------------------- end interface ----------------
!     loop indexes                                                      
      INTEGER j 
!     ==================================================================
      dfteval(1) = DCMPLX(dftout1(1),0.d0) 
      DO j = 2,N/2 
         dfteval(j) = DCMPLX(dftout1(j),dftout1(N+2-j)) 
         dfteval(N/2+j) = DCMPLX(dftout2(N/2+j),dftout2(N/2+2-j)) 
      ENDDO 
      dfteval(N/2+1) = DCMPLX(dftout1(N/2+1),0.d0) 
                                                                        
      RETURN 
      END                                           


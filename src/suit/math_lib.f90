! ===========LIBRARY math_lib=====================                      
! note: many routines to be removed in fortran90                        
! MATHEMATICAL ROUTINES:                                      
!
!   3D LINEAR ALGEBRA OK
! prscal	scalar product of two 3D vectors                               
! prvec		vector product of two 3D vectors                               
! vsize		size of a 3D vector                                            
! rotmt         rotation matrix around some coordinate axis
! prodmv	product of a 3x3 matrix by a 3D vector  (could be removed)  
! prodmm        product of two  3x3 matrices            (could be removed)
! trsp3         transposed 3x3 matrix                   (could be removed) 
! lincom	linear combination of two 3D vectors    (could be removed)  
! mtp_ref       Graham_Schmidt with two assigned vectors
!   GENERAL LINEAR ALGEBRA 
! matin		matrix inversion (Gauss method) (WARNING: old fortran loops)
! eye		identity matrix 
! inv22         inversion of a 2x2 matrix (WARNING: not safe for det=0)  
! bilin         bilinear form of two general vectors 
!   LEAST SQUARES OK                                                       
! tchinv        standard interface for tchol, inver  
! tchol		Tcholesky factorization of a positive-defined matrix           
! inver         Tcholesky inversion                                     
! qr_inv        inversion by eigenvector/eigenvalues
! snorm		weighed RMS norm sqrt(xTWx/n)                                  
! snormd        RMS norm sqrt(xTx/n)   
! convertcov    conversion of covariance matrix
! convertnor    conversion of normal matrix  
! norprs        normal matrix propagation for a NxN matrix (uses convertnor)  
! linfi3        linear regression
! quadratic_fit quadratic least squares fit (also linear)
! cubic_fit     cubic least squares fit
!   SORTING  OK                                                           
! heapsort      sort by real, leaving in place and storing indexes      
! heapsorti     sort by integer, leaving in place and storing indexes   
!   ROUNDING OK
! truncat       truncation accounting for rounding off (to be removed) 
!                                                                       
! HEADERS math_lib.o: none                                   
! MODULES: fund_const.o    
!                       
! =========================================================== 
! 3D LINEAR ALGEBRA
! =========================================================== 
!  ***************************************************************      
!  *                                                             *      
!  *                          P R S C A L                        *      
!  *                                                             *      
!  *                 Scalar product:   C = <A, B>                *      
!  *                                                             *      
!  ***************************************************************      
double precision function prscal(a,b) 
  double precision a(3),b(3) 
  prscal=a(1)*b(1)+a(2)*b(2)+a(3)*b(3) 
  return 
END function prscal
!  ***************************************************************      
!  *                                                             *      
!  *                          P R V E C                          *      
!  *                                                             *      
!  *                 Vector product:   C = A x B                 *      
!  *                                                             *      
!  ***************************************************************      

subroutine prvec(a,b,c) 
  IMPLICIT NONE
  DOUBLE PRECISION, DIMENSION(3) :: a,b,c 
  c(1)=a(2)*b(3)-a(3)*b(2) 
  c(2)=a(3)*b(1)-a(1)*b(3) 
  c(3)=a(1)*b(2)-a(2)*b(1) 
  return 
END subroutine prvec
!  ***************************************************************      
!  *                                                             *      
!  *                          V S I Z E                          *      
!  *                                                             *      
!  *                     Size of a 3-D vector                    *      
!  *                                                             *      
!  ***************************************************************      
!                                                                       
DOUBLE PRECISION FUNCTION vsize(x) 
  IMPLICIT NONE 
  DOUBLE PRECISION x(3),s 
  s=x(1)*x(1)+x(2)*x(2)+x(3)*x(3) 
  vsize=SQRT(s) 
END FUNCTION vsize

! Author: Mario Carpino (carpino@brera.mi.astro.it)                     
! Version: June 5, 1996                                                 
!                                                                       
!  ***************************************************************      
!  *                                                             *      
!  *                         P R O D M V                         *      
!  *                                                             *      
!  *            Product of a 3x3 matrix by a 3-D vector          *      
!  *                           y = Ax                            *      
!  *                                                             *      
!  ***************************************************************      
!                                                                       
SUBROUTINE prodmv(y,a,x) 
  IMPLICIT NONE 
  DOUBLE PRECISION, INTENT(IN) :: a(3,3)
  DOUBLE PRECISION :: x(3),y(3)
  DOUBLE PRECISION s, z(3)
  INTEGER j,l 
  z=x                                                                     
  DO  j=1,3 
     s=0.d0 
     DO  l=1,3 
        s=s+a(j,l)*z(l) 
     ENDDO
     y(j)=s 
   ENDDO
 END SUBROUTINE prodmv
! Author: Mario Carpino (carpino@brera.mi.astro.it)                     
! Version: June 5, 1996                                                 
!                                                                       
!  ***************************************************************      
!  *                                                             *      
!  *                          R O T M T                          *      
!  *                                                             *      
!  *              Rotation matrix around k axis	                 *      
!  *                                                             *      
!  ***************************************************************      
!                                                                       
!  If X are "old" coordinates and X' are "new" coordinates (referred    
!  to a frame which is rotated by an angle alpha around k axis in       
!  direct sense), then X' = R X                                         
!                                                                       
subroutine rotmt(alpha,r,k) 
  implicit none 
  double precision alpha,r(3,3) 
  double precision cosa,sina 
  integer k,i1,i2,i3 
  if(k.lt.1.or.k.gt.3)stop ' **** ROTMT: k = ??? ****' 
  cosa=cos(alpha) 
  sina=sin(alpha) 
  i1=k 
  if(i1.gt.3)i1=i1-3 
  i2=i1+1 
  if(i2.gt.3)i2=i2-3 
  i3=i2+1 
  if(i3.gt.3)i3=i3-3 
  r(i1,i1)=1.d0 
  r(i1,i2)=0.d0 
  r(i1,i3)=0.d0 
  r(i2,i1)=0.d0 
  r(i2,i2)=cosa 
  r(i2,i3)=sina 
  r(i3,i1)=0.d0 
  r(i3,i2)=-sina 
  r(i3,i3)=cosa
END subroutine rotmt
 
                                
! Author: Mario Carpino (carpino@brera.mi.astro.it)                     
! Version: June 5, 1996                                                 
!  ***************************************************************      
!  *                                                             *      
!  *                         P R O D M M                         *      
!  *                                                             *      
!  *                 Product of two 3x3 matrices                 *      
!  *                           A = BC                            *      
!  *                                                             *      
!  ***************************************************************      
SUBROUTINE prodmm(a,b,c) 
  IMPLICIT NONE 
  DOUBLE PRECISION :: b(3,3),c(3,3)
  DOUBLE PRECISION :: a(3,3)
  DOUBLE PRECISION w(3,3),z(3,3),s 
  INTEGER j,k,l
  z=c
  w=b 
  DO j=1,3 
     DO  k=1,3 
        s=0.d0 
        DO l=1,3 
          s=s+w(j,l)*z(l,k)
        ENDDO 
        a(j,k)=s 
     ENDDO
   ENDDO
END SUBROUTINE prodmm
! Copyright (C) 1999 by Mario Carpino (carpino@brera.mi.astro.it)       
! Version: November 10, 1999                                            
! --------------------------------------------------------------------- 
!                                                                       
!  *****************************************************************    
!  *                                                               *    
!  *                          T R S P 3                            *    
!  *                                                               *    
!  *                     Transpose of a matrix                     *    
!  *                                                               *    
!  *****************************************************************    
!                                                                       
! INPUT:    R         -  Input matrix                                   
!                                                                       
! OUTPUT:   R         -  R' (transpose of input matrix)                 
!                                                                       
SUBROUTINE trsp3(r) 
  IMPLICIT NONE                                                                       
  DOUBLE PRECISION r(3,3)                                                             
  INTEGER i,j 
  DOUBLE PRECISION rt                                                               
  DO 1 i=1,2 
     DO 1 j=i+1,3 
        rt=r(i,j) 
        r(i,j)=r(j,i) 
        r(j,i)=rt 
1 CONTINUE                                                                         
END SUBROUTINE trsp3  
!   
subroutine lincom(v1,a,v2,b,v3) 
!                                                                       
!    Computes the linear combination vector v3 of the                   
!    three-elements vectors v1,v2 with coefficients a,b                 
!                                                                       
  IMPLICIT NONE 
  double precision v1(3),v2(3),v3(3),a,b 
  INTEGER i
  do i=1,3 
     v3(i)=a*v1(i)+b*v2(i) 
  enddo
END subroutine lincom
! =================================================================     
! reference system with tpno as first axis, +dx as third axis           
! to be used with dx= normal to ecliptic                                
      SUBROUTINE mtp_ref(tpno,dx,v3,vt3) 
      IMPLICIT NONE 
      DOUBLE PRECISION tpno(3),dx(3),v3(3,3),vt3(3,3) 
! end interface                                                         
      DOUBLE PRECISION vl,vsize,vv,prscal 
      INTEGER i 
! the first vector is normal to the target plane                         
      v3(1:3,1)=tpno ! CALL vcopy(3,tpno,v3(1,1)) 
! the third vector is the projection of -dx on the plane orthogonal to t
      vv=prscal(tpno,dx) 
      DO i=1,3 
!         v3(i,3)=-(dx(i)-vv*tpno(i)) ! obsolete to have - dx (e.g. Earth velocity...)
         v3(i,3)=(dx(i)-vv*tpno(i)) 
      ENDDO 
! reduced to unit length                                                
      vl=vsize(v3(1,3)) 
      DO i=1,3 
         v3(i,3)=v3(i,3)/vl 
      ENDDO 
! the second vector must be orthogonal to both                          
      CALL prvec(tpno,v3(1,3),v3(1,2)) 
! paranoia check of unit length                                         
      vl=vsize(v3(1,2)) 
      DO i=1,3 
         v3(i,2)=v3(i,2)/vl 
      ENDDO 
! the inverse is the transposed matrix                                  
      vt3=TRANSPOSE(v3) !CALL transp(v3,3,3,vt3) 
      RETURN 
      END SUBROUTINE mtp_ref 
! =========================================================== 
! GENERAL LINEAR ALGEBRA
! =========================================================== 
!  VERSOR
  SUBROUTINE versor(n,a,epsi,b,vl,ize)  
    integer n,i 
    logical ize 
    double precision a(n),b(n),epsi,vl 
    vl=sqrt(DOT_PRODUCT(a,a)) 
    if(vl.lt.epsi)then 
       ize=.true. 
    else 
       ize=.false. 
       do  i=1,n 
          b(i)=a(i)/vl 
       enddo
    endif
  END  SUBROUTINE versor
!  VERSOR_Q
  SUBROUTINE versor_q(n,a,epsi,b,vl,ize)  
    USE fund_const, ONLY: qkind
    IMPLICIT NONE
    INTEGER n,i 
    LOGICAL ize 
    REAL(KIND=qkind) a(n),b(n),epsi,vl 
    vl=sqrt(DOT_PRODUCT(a,a)) 
    if(vl.lt.epsi)then 
       ize=.true. 
    else 
       ize=.false. 
       do  i=1,n 
          b(i)=a(i)/vl 
       enddo
    endif
  END  SUBROUTINE versor_q

! ================================================                      
!                                                                       
!  ***************************************************************      
!  *                                                             *      
!  *                          M A T I N                          *      
!  *                                                             *      
!  *              Gauss' method for matrix inversion             *      
!  *         and solution of systems of linear equations         *      
!  *                                                             *      
!  ***************************************************************      
!                                                                       
! INPUT:    A(i,j)    -   Matrix of coefficients of the system          
!           N         -   Order of A(i,j)                               
!           L         -   Number of linear systems to be solved (the    
!                         right hand  sides are stored in A(i,N+j),     
!                         j=1,L)                                        
!           NL        -   First dimension of A(i,j)                     
!           INVOP     -   If =1 the inverse of A(i,j) is computed       
!                         explicitly,                                   
!                         if =0 only the linear systems are solved      
!                                                                       
! OUTPUT:   A(i,j)    -   If INVOP=1, contains the inverse of the input 
!                         if INVOP=0, contains the triangular           
!                         factorization of the input matrix.            
!                         In both cases A(i,N+j), j=1,L contain the     
!                         solutions of the input systems.               
!           DET       -   Determinant of the input matrix A(i,j)        
!           ISING     -   If =0 the algorithm was completed             
!                         successfully,                                 
!                         if =1 the input matrix was singular           
!                                                                       
 subroutine matin(a,det,n,l,nl,ising,invop) 
   IMPLICIT NONE
   INTEGER, INTENT(IN) :: n,l,nl,invop
   INTEGER, INTENT(OUT) :: ising
   DOUBLE PRECISION, INTENT(INOUT) ::  a(nl,n+l)
   DOUBLE PRECISION, INTENT(OUT) :: det
   INTEGER, PARAMETER :: nmax=1000 
   DOUBLE PRECISION vet(nmax) ! workspace 
   INTEGER  p(nmax),j,i,k,nc,npiv,nsc
   DOUBLE PRECISION amax,sc
   if(n.gt.nmax)stop ' **** matin: n > nmax ****' 
   do  j=1,n 
      p(j)=j 
   enddo
   do 4 j=1,n-1 
      amax=0.d0 
      DO i=j,n 
         if(amax.ge.dabs(a(i,j))) CYCLE
         npiv=i 
         amax=dabs(a(i,j)) 
      ENDDO
! modification 2006: error flag must be set
      IF(amax.eq.0.d0)THEN
         ising=1
         RETURN 
      ENDIF
! end modification
      if(npiv.eq.j)goto 6 
      nsc=p(npiv) 
      p(npiv)=p(j) 
      p(j)=nsc 
      do  i=1,n+l 
         sc=a(npiv,i) 
         a(npiv,i)=a(j,i) 
         a(j,i)=sc
      ENDDO 
    6 DO i=j+1,n 
         a(i,j)=a(i,j)/a(j,j)
      ENDDO 
      do 4 i=j+1,n 
      do 4 k=j+1,n+l 
4  a(i,k)=a(i,k)-a(i,j)*a(j,k) 
   if(l.eq.0)goto 7 
   do 9 i=1,l 
   do 9 j=n,1,-1 
      if(j.eq.n)goto 9 
      do  k=j+1,n 
         a(j,n+i)=a(j,n+i)-a(j,k)*a(k,n+i)
      enddo 
9  a(j,n+i)=a(j,n+i)/a(j,j) 
7  det=1.d0 
   do j=1,n 
      det=det*a(j,j) 
   enddo
   ising=0 
   if(invop.eq.0)return 
   if(n.ne.1)goto 20 
   a(1,1)=1.d0/a(1,1) 
   return 
20 do 12 i=2,n 
   do 12 j=1,i-1 
      a(i,j)=-a(i,j) 
      if(j+2.gt.i)goto 12 
      do  k=j+1,i-1 
         a(i,j)=a(i,j)-a(i,k)*a(k,j)
      enddo 
12 continue 
   do 15 k=n,1,-1 
      do 14 i=1,n 
         vet(i)=0.d0 
         if(i.eq.k)vet(i)=1.d0 
         if(k.gt.i)vet(i)=a(k,i) 
         if(k.eq.n)goto 14 
         do  j=k+1,n 
            vet(i)=vet(i)-a(k,j)*a(j,i)
         enddo
   14 continue 
      sc=a(k,k) 
      do 15  i=1,n 
         a(k,i)=vet(i)/sc
15 continue
18 nc=0 
   do 16 j=1,n 
      if(j.eq.p(j))goto 16 
      nsc=p(j) 
      p(j)=p(nsc) 
      p(nsc)=nsc 
      do 17 i=1,n 
         sc=a(i,nsc) 
         a(i,nsc)=a(i,j) 
17       a(i,j)=sc 
         nc=1 
16 continue 
   if(nc.eq.1)goto 18 
   return 
END subroutine matin                                          
!==============================================                         
! identity matrix ndim x ndim in output a                                
subroutine eye(ndim,a) 
  implicit none 
  integer, intent(in) :: ndim
  double precision, intent(out) :: a(ndim,ndim) 
  integer i,j 
  do  i=1,ndim 
     do  j=1,ndim 
        if(i.eq.j)then 
           a(i,j)=1.d0 
        else 
           a(i,j)=0.d0 
        endif
     enddo
  enddo
END subroutine eye
! =======================================================               
! computation of a bilinear form                                        
double precision function bilin(v1,v2,nv,w,nx,n) 
  implicit none 
  integer, intent(in) :: nv,nx,n 
  double precision, intent(in) :: v1(nv),v2(nv),w(nx,n) 
  integer i,j 
!  control on dimensions                                                
  if(nv.ne.n)write(*,*)' dimension quadr, nv=',nv,' ne n=',n 
  if(n.gt.nx)write(*,*)' dimension quadr, n=',n,' gt nx=',nx 
  bilin=0.d0 
  do 10 i=1,n 
     do 10 j=1,n 
10      bilin=bilin+v1(i)*w(i,j)*v2(j) 
END function bilin                                               

! ===================================================================== 
! inverse of a 2x2 matrix
! warning: this can give NaN if deta=0
subroutine inv22(a,b,deta) 
  implicit none 
  double precision, intent(in) :: a(2,2)
  double precision, intent(out) :: b(2,2),deta 
  deta=a(1,1)*a(2,2)-a(1,2)*a(2,1) 
  if(deta.eq.0.d0)then 
     write(*,*)' deta is zero ' 
  endif
  b(1,1)=a(2,2)/deta 
  b(2,2)=a(1,1)/deta 
  b(1,2)=-a(1,2)/deta 
  b(2,1)=-a(2,1)/deta 
END subroutine inv22

! ===================================================================== 
! inverse of a 2x2 matrix
! warning: this can give NaN if deta=0
SUBROUTINE inv22_q(a,b,deta) 
  USE fund_const, ONLY: qkind
  IMPLICIT NONE
  REAL(KIND=qkind), INTENT(in) :: a(2,2)
  REAL(KIND=qkind), INTENT(out) :: b(2,2),deta 
  deta=a(1,1)*a(2,2)-a(1,2)*a(2,1) 
  if(deta.eq.0.d0)then 
     write(*,*)' deta is zero ' 
  endif
  b(1,1)=a(2,2)/deta 
  b(2,2)=a(1,1)/deta 
  b(1,2)=-a(1,2)/deta 
  b(2,1)=-a(2,1)/deta 
END SUBROUTINE inv22_q
!
! =========================================================== 
! LEAST SQUARES
! =========================================================== 
! ===========================================================           
! T C H I N V
! Cholewski method for inversion                                        
! requires a workspace vector, but does not have dimension limits       
! ======================================================                
subroutine tchinv(c,n,cinv,ws,indp) 
  USE output_control
  implicit none 
! input: dimension of input matrix c                                    
  integer, intent(in) :: n 
  double precision, intent(in) :: c(n,n) 
! output: possible degeneracy (0=no, else row where pivot is too small) 
  integer, intent(out)::  indp 
  double precision, intent(out) :: cinv(n,n) 
! workspace                                                             
  double precision, intent(out) :: ws(n) 
! end interface                                                         
  integer i
  double precision err,omax,omin,cond 
! =========================================                             
  cinv=c ! call mcopy(n,n,c,cinv) 
  err=epsilon(1.d0)*100 
! matrix factorized by Tcholesky method:                                
  call tchol(cinv,n,n,indp,err) 
  if(indp.ne.0.and.verb_matrix.gt.9)then 
     write(*,*)' indp=',indp,' in tchol' 
  endif
! ===========================================================           
! Control of conditioning number and inversion of matrix                
  omax=cinv(1,1) 
  omin=cinv(1,1) 
  DO i=2,n 
     if (cinv(i,i).gt.omax) then 
        omax=cinv(i,i) 
     endif
     if (cinv(i,i).lt.omin) then 
        omin=cinv(i,i) 
     endif
  ENDDO
  cond=(omax/omin)**2 
  IF(cond.gt.1.d12.and.verb_matrix.gt.9)write(*,111)n,n,cond 
111 format(' Conditioning of ',i2,'x',i2,'matrix=',d12.4) 
! ===========================================================           
! inversion                                                             
  call inver(cinv,ws,n,n) 
END subroutine tchinv
! ===========================================================           
!                                                                       
!                            t c h o l                                  
!                                                                       
!      Tcholesky factorization of a positive-defined matrix             
!                                                                       
!             originally written by prof. Luigi Mussio                  
!         Istituto di Topografia del Politecnico di Milano              
!                                                                       
! ===========================================================           
!                                                                       
!    warning: only the upper triangle of the matrix is processed        
!                                                                       
!                                                                       
! input:    a(nmax,nmax)   - matrix to be processed                     
!           nmax           - first dimension of a as declared in the dim
!                            statement of the calling program           
!           n              - actual dimension of a                      
!           err            - control on pivot (if <0, control is automat
!                                                                       
!                                                                       
! output:   a              - contains the tcholesky factorization of the
!                            matrix (to be supplied to subroutine inver)
!           indp           - if it is different from zero, the input mat
!                            was not positive defined and the algorithm 
!                            does not work                              
!                                                                       
 subroutine tchol(a,nmax,n,indp,err) 
   implicit none 
   integer, intent(in) :: n, nmax
   double precision, intent(in) :: err
   integer, intent(out) :: indp 
   double precision, intent(inout) :: a(nmax,n)
   double precision :: errrel
! ========END INTERFACE=========================                        
   double precision as,errloc,s 
   integer i,k,ii,l,j 
   as=0.d0 
   do  i=1,n 
      do k=1,i 
         as=max(as,abs(a(k,i)))
      enddo
   enddo 
! in this version we are not checking the roundoff value                
! of the machine, because the control err is assumed to hve
! already BEEN compared with the rounding off; 
! thus err must be positive                                             
   IF(err.le.0.d0)THEN
      indp=1
      RETURN
   ELSE
      errrel=err*as
   ENDIF 
   do 40 i=1,n 
      l=i-1 
      s=0.d0 
      if(l.eq.0) goto 15 
      do  k=1,l 
         s=s+a(k,i)*a(k,i) 
      enddo
15    a(i,i)=a(i,i)-s 
      if(a(i,i).le.errrel)then 
         indp=i 
         return 
      endif
      a(i,i)=sqrt(a(i,i)) 
      if(i.eq.n) goto 40 
      ii=i+1 
      DO 30 j=ii,n 
         s=0.d0 
         if(l.eq.0)goto 25 
         do  k=1,l 
            s=s+a(k,i)*a(k,j) 
         enddo
25       a(i,j)=(a(i,j)-s)/a(i,i) 
30    ENDDO
40 ENDDO
   indp=0 
   return 
 END subroutine tchol
! ===========================================================           
!                                                                       
!                            i n v e r                                  
!                                                                       
!              inversion of a positive-defined matrix                   
!                                                                       
!                  written by prof. luigi mussio                        
!         istituto di topografia del politecnico di milano              
!                                                                       
!                                                                       
! input:    a(nmax,nmax)   - input matrix, previously factorized by     
!                            subroutine tchol                           
!           v(nmax)        - work area                                  
!           nmax           - first dimension of a as declared in the dim
!                            statement of the calling program           
!           n              - actual dimension of a                      
!                                                                       
!                                                                       
! output:   a              - inverse of input matrix                    
!                                                                       
!                                                                       
 subroutine inver(a,v,nmax,n) 
   implicit none 
   integer, INTENT(IN) :: nmax,n 
   double precision, INTENT(INOUT) :: a(nmax,n),v(n) 
! ==========END INTERFACE=========                                      
   integer i,l,k,j 
   double precision s 
! ================================                                      
   do 8 i=n,1,-1 
      l=i+1 
      if(i.eq.n)goto 4 
      do 3 j=n,l,-1 
         s=0.d0 
         do 2 k=n,l,-1 
            if(j.lt.k)goto 1 
            s=s+a(i,k)*a(k,j) 
            goto 2 
1           s=s+a(i,k)*a(j,k) 
2        ENDDO
         v(j)=-s/a(i,i) 
3     enddo
4     s=0.d0 
      if(i.eq.n)goto 7 
      do  k=n,l,-1 
         s=s+a(i,k)*v(k) 
      enddo
      do  j=n,l,-1 
         a(i,j)=v(j) 
      enddo
7     a(i,i)=(1.d0/a(i,i)-s)/a(i,i) 
8  ENDDO
   do  i=2,n 
      do  j=1,i-1 
        a(i,j)=a(j,i)
      enddo
   enddo 
 END subroutine inver      


! ===========================================================           
! T C H I N V _Q
! Cholewski method for inversion                                        
! requires a workspace vector, but does not have dimension limits       
! ======================================================                
subroutine tchinv_q(c,n,cinv,ws,indp) 
  USE fund_const, ONLY: qkind
  USE output_control
  implicit none 
! input: dimension of input matrix c                                    
  integer, intent(in) :: n 
  REAL(KIND=qkind), intent(in) :: c(n,n) 
! output: possible degeneracy (0=no, else row where pivot is too small) 
  integer, intent(out)::  indp 
  REAL(KIND=qkind), intent(out) :: cinv(n,n) 
! workspace                                                             
  REAL(KIND=qkind), intent(out) :: ws(n) 
! end interface                                                         
  integer i
  REAL(KIND=qkind) err,omax,omin,cond 
! =========================================                             
  cinv=c ! call mcopy(n,n,c,cinv) 
  err=epsilon(1.q0)*1.q10 
! matrix factorized by Tcholesky method:                                
  call tchol_q(cinv,n,n,indp,err) 
  if(indp.ne.0.and.verb_matrix.gt.9)then 
     write(*,*)' indp=',indp,' in tchol' 
  endif
! ===========================================================           
! Control of conditioning number and inversion of matrix                
  omax=cinv(1,1) 
  omin=cinv(1,1) 
  DO i=2,n 
     if (cinv(i,i).gt.omax) then 
        omax=cinv(i,i) 
     endif
     if (cinv(i,i).lt.omin) then 
        omin=cinv(i,i) 
     endif
  ENDDO
  cond=(omax/omin)**2 
  IF(cond.gt.1.d12.and.verb_matrix.gt.9)write(*,111)n,n,cond 
111 format(' Conditioning of ',i2,'x',i2,'matrix=',d12.4) 
! ===========================================================           
! inversion                                                             
  call inver_q(cinv,ws,n,n) 
END subroutine tchinv_q
! ===========================================================           
!                                                                       
!                            t c h o l _ q 
!                                                                       
!      Tcholesky factorization of a positive-defined matrix             
!                                                                       
!             originally written by prof. Luigi Mussio                  
!         Istituto di Topografia del Politecnico di Milano              
!                                                                       
! ===========================================================           
!                                                                       
!    warning: only the upper triangle of the matrix is processed        
!                                                                       
!                                                                       
! input:    a(nmax,nmax)   - matrix to be processed                     
!           nmax           - first dimension of a as declared in the dim
!                            statement of the calling program           
!           n              - actual dimension of a                      
!           err            - control on pivot (if <0, control is automat
!                                                                       
!                                                                       
! output:   a              - contains the tcholesky factorization of the
!                            matrix (to be supplied to subroutine inver)
!           indp           - if it is different from zero, the input mat
!                            was not positive defined and the algorithm 
!                            does not work                              
!                                                                       
 subroutine tchol_q(a,nmax,n,indp,err) 
   USE fund_const, ONLY: qkind
   implicit none 
   integer, intent(in) :: n, nmax
   REAL(KIND=qkind), intent(in) :: err
   integer, intent(out) :: indp 
   REAL(KIND=qkind), intent(inout) :: a(nmax,n)
   REAL(KIND=qkind) :: errrel
! ========END INTERFACE=========================                        
   REAL(KIND=qkind) as,errloc,s 
   integer i,k,ii,l,j 
   as=0.d0 
   do  i=1,n 
      do k=1,i 
         as=max(as,abs(a(k,i)))
      enddo
   enddo 
! in this version we are not checking the roundoff value                
! of the machine, because the control err is assumed to hve
! already BEEN compared with the rounding off; 
! thus err must be positive                                             
   IF(err.le.0.d0)THEN
      indp=1
      RETURN
   ELSE
      errrel=err*as
   ENDIF 
   do 40 i=1,n 
      l=i-1 
      s=0.d0 
      if(l.eq.0) goto 15 
      do  k=1,l 
         s=s+a(k,i)*a(k,i) 
      enddo
15    a(i,i)=a(i,i)-s 
      if(a(i,i).le.errrel)then 
         indp=i 
         return 
      endif
      a(i,i)=sqrt(a(i,i)) 
      if(i.eq.n) goto 40 
      ii=i+1 
      DO 30 j=ii,n 
         s=0.d0 
         if(l.eq.0)goto 25 
         do  k=1,l 
            s=s+a(k,i)*a(k,j) 
         enddo
25       a(i,j)=(a(i,j)-s)/a(i,i) 
30    ENDDO
40 ENDDO
   indp=0 
   return 
 END subroutine tchol_q
! ===========================================================           
!                                                                       
!                            i n v e r _ q                                 
!                                                                       
!              inversion of a positive-defined matrix                   
!                                                                       
!                  written by prof. luigi mussio                        
!         istituto di topografia del politecnico di milano              
!                                                                       
!                                                                       
! input:    a(nmax,nmax)   - input matrix, previously factorized by     
!                            subroutine tchol                           
!           v(nmax)        - work area                                  
!           nmax           - first dimension of a as declared in the dim
!                            statement of the calling program           
!           n              - actual dimension of a                      
!                                                                       
!                                                                       
! output:   a              - inverse of input matrix                    
!                                                                       
!                                                                       
 subroutine inver_q(a,v,nmax,n)
   USE fund_const, ONLY: qkind
   implicit none 
   integer, INTENT(IN) :: nmax,n 
   REAL(KIND=qkind), INTENT(INOUT) :: a(nmax,n),v(n) 
! ==========END INTERFACE=========                                      
   integer i,l,k,j 
   REAL(KIND=qkind) s 
! ================================                                      
   do 8 i=n,1,-1 
      l=i+1 
      if(i.eq.n)goto 4 
      do 3 j=n,l,-1 
         s=0.d0 
         do 2 k=n,l,-1 
            if(j.lt.k)goto 1 
            s=s+a(i,k)*a(k,j) 
            goto 2 
1           s=s+a(i,k)*a(j,k) 
2        ENDDO
         v(j)=-s/a(i,i) 
3     enddo
4     s=0.d0 
      if(i.eq.n)goto 7 
      do  k=n,l,-1 
         s=s+a(i,k)*v(k) 
      enddo
      do  j=n,l,-1 
         a(i,j)=v(j) 
      enddo
7     a(i,i)=(1.d0/a(i,i)-s)/a(i,i) 
8  ENDDO
   do  i=2,n 
      do  j=1,i-1 
        a(i,j)=a(j,i)
      enddo
   enddo 
 END subroutine inver_q      

!====================================================                   
 subroutine qr_inv(a,b,n,izer,det,aval) 
   implicit none 
!*****************************************************                  
! a = input matrix                                                      
! b = output matrix    
! they are n by n; note n has a maximum of nmax
! aval= eigenvalues vector                                              
! q = eigenvectors matrix                                               
!*****************************************************
   integer, parameter :: nmax=10
   integer, intent(in) :: n
   double precision, intent(in) :: a(n,n)
   double precision, intent(out) :: b(n,n)
   integer,intent(out) :: izer
   double precision q(nmax,nmax),dstar(nmax,nmax) 
   double precision qt(nmax,nmax),tmp(nmax,nmax) 
   double precision aval(nmax),w1(nmax),w2(nmax) 
   double precision det 
   integer i,j,ierr 
   call rs(n,n,a,aval(1:n),1,q(1:n,1:n),w1,w2,ierr)
   dstar=0.d0 
   det=1.d0 
   izer=0 
   do 70 i=1,n 
      det=det*aval(i) 
      if (aval(i).gt.0.d0) then 
         dstar(i,i)=1.d0/aval(i) 
      else 
         izer=izer+1 
      endif
70 enddo
   qt(1:n,1:n)=TRANSPOSE(q(1:n,1:n)) 
   tmp(1:n,1:n)=MATMUL(dstar(1:n,1:n),qt(1:n,1:n))  
   b=MATMUL(q(1:n,1:n),tmp(1:n,1:n)) 
 END subroutine qr_inv
! ==============================================                        
! SNORM norm of vector v according to metric a                          
DOUBLE PRECISION FUNCTION snorm(v,a,n,nnx) 
  IMPLICIT NONE 
  INTEGER, INTENT(IN) :: nnx,n
  INTEGER i,k 
  DOUBLE PRECISION, INTENT(IN) ::  v(nnx),a(nnx,nnx) 
  snorm=0.d0 
  DO  i=1,n 
     DO  k=1,n 
        snorm=snorm+v(i)*v(k)*a(i,k) 
     ENDDO
  ENDDO
  IF(snorm.gt.0.d0)THEN
     snorm=sqrt(snorm/n) 
  ELSE
     snorm=0.d0
  ENDIF
END FUNCTION snorm
! ===========================================                           
! SNORMD RMS norm of vector v with diagonal metric a                    
! warning: if a weight is zero, the component is not                    
! counted in the vector length                                          
DOUBLE PRECISION FUNCTION snormd(v,a,n,nused) 
  IMPLICIT NONE 
  INTEGER,INTENT(IN) :: n
  DOUBLE PRECISION, INTENT(IN) :: v(n),a(n) 
  INTEGER, INTENT(OUT) :: nused 
  INTEGER i
  snormd=0.d0 
  nused=0 
  DO  i=1,n 
     snormd=snormd+(v(i)**2)*a(i) 
     IF(a(i).gt.0.d0)nused=nused+1 
  ENDDO
  IF(snormd.gt.0.d0.and.nused.gt.0)THEN
     snormd=sqrt(snormd/nused)
  ELSE
     snormd=0.d0
  ENDIF 
END FUNCTION snormd
! =========================================================== 
! CONVERTCOV, CONVERTNOR
! conversion by similarity of covariance and normal matrix
! ===========================================================  
! conversion of covariance matrix, version nd x nd
                            
SUBROUTINE convertcov(nd,g0,derpar,gk) 
  IMPLICIT NONE 
! input: g0=covariance matrix in elements E, derpar=\partial X/\partial E
! output: gk covariance matrix in elements X
  INTEGER, INTENT(IN) :: nd                            
  DOUBLE PRECISION, INTENT(IN) :: g0(nd,nd),derpar(nd,nd)
  DOUBLE PRECISION, INTENT(INOUT) :: gk(nd,nd) ! to allow for g0=gk
! workspace                                                             
  DOUBLE PRECISION tmp6(nd,nd),derpart(nd,nd) 
! compute normal and covariance matrix for keplerian elements           
  tmp6=MATMUL(derpar,g0) 
  derpart=TRANSPOSE(derpar) 
  gk=MATMUL(tmp6,derpart) 
END SUBROUTINE convertcov

! version for nd x nd matrices with 6 x 6 jacobian
SUBROUTINE convertcovdp(nd,g0,derpar,gk) 
  IMPLICIT NONE 
! input: g0=covariance matrix in elements E, derpar=\partial X/\partial E
! output: gk covariance matrix in elements X 
  INTEGER, INTENT(IN) :: nd                           
  DOUBLE PRECISION, INTENT(IN) :: g0(nd,nd),derpar(6,6)
  DOUBLE PRECISION, INTENT(INOUT) :: gk(nd,nd) ! to allow for g0=gk
! workspace                                                             
  DOUBLE PRECISION tmp6(6,6),derpart(6,6) 
! compute normal and covariance matrix for keplerian elements           
  tmp6=MATMUL(derpar,g0(1:6,1:6)) 
  derpart=TRANSPOSE(derpar) 
  gk(1:6,1:6)=MATMUL(tmp6,derpart)
  IF(nd.eq.6) RETURN
  gk(7:nd,7:nd)=g0(7:nd,7:nd)
  gk(7:nd,1:6)=MATMUL(g0(7:nd,1:6),derpart)
  gk(1:6,7:nd)=TRANSPOSE(gk(7:nd,1:6)) 
END SUBROUTINE convertcovdp
                                                                        
! conversion of normal matrix nd x nd                      
SUBROUTINE convertnor(nd,c0,derpar,ck) 
  IMPLICIT NONE 
! input: c0=normal matrix in elements E, derpar=\partial E/\partial X   
! output: ck= normal matrix in elements X
  INTEGER, INTENT(IN) :: nd                               
  DOUBLE PRECISION, INTENT(IN) ::  c0(nd,nd),derpar(nd,nd)
  DOUBLE PRECISION, INTENT(INOUT) :: ck(nd,nd) 
! workspace                                                             
  DOUBLE PRECISION tmp6(nd,nd),derpart(nd,nd) 
! compute normal and covariance matrix for keplerian elements           
  derpart=TRANSPOSE(derpar) 
  tmp6=MATMUL(derpart, c0)  
  ck=MATMUL(tmp6,derpar)    
END SUBROUTINE convertnor

! conversion of normal matrix  nd x nd with jacobian 6 x 6                    
SUBROUTINE convertnordp(nd,c0,derpar,ck) 
  IMPLICIT NONE 
! input: c0=normal matrix in elements E, derpar=\partial E/\partial X   
! output: ck= normal matrix in elements X
  INTEGER, INTENT(IN) :: nd                                
  DOUBLE PRECISION, INTENT(IN) ::  c0(nd,nd),derpar(6,6)
  DOUBLE PRECISION, INTENT(INOUT) :: ck(nd,nd) 
! workspace                                                             
  DOUBLE PRECISION tmp6(6,6),derpart(6,6) 
! compute normal and covariance matrix for keplerian elements           
  derpart=TRANSPOSE(derpar) 
  tmp6=MATMUL(derpart, c0(1:6,1:6))  
  ck(1:6,1:6)=MATMUL(tmp6,derpar)    
  IF(nd.eq.6) RETURN
  ck(7:nd,7:nd)=c0(7:nd,7:nd)
  ck(7:nd,1:6)=MATMUL(c0(7:nd,1:6),derpar)
  ck(1:6,7:nd)=TRANSPOSE(ck(7:nd,1:6))   
END SUBROUTINE convertnordp

! Copyright (C) 1998 by Mario Carpino (carpino@brera.mi.astro.it)       
! Version: November 2009 AMC, GT, DF                                           
! --------------------------------------------------------------------- 
!                                                                       
!  *****************************************************************    
!  *                                                               *    
!  *                         N O R P R S                           *    
!  *                                                               *    
!  *              Propagation of a normal matrix (NxN)             *    
!  *                                                               *    
!  *****************************************************************    
!                                                                       
! INPUT:    NOR1      -  Normal matrix (X1 coord.)                      
!           JAC       -  Jacobian matrix dX2/dX1                        
!           N         -  Dimension (logical AND physical) of matrices   
!                                                                       
! OUTPUT:   NOR2      -  Normal matrix (X2 coord.) =                    
!                           = JAC'^(-1) NOR1 JAC^(-1)                   
!           ERROR     -  Error flag (singular Jacobian matrix)          
!                

! version for nd x nd matrices with 6 x 6 jacobian
SUBROUTINE norprsdp(nor1,jac,nd,nor2,error) 
  IMPLICIT NONE
  INTEGER, INTENT(IN) ::  nd 
  DOUBLE PRECISION, INTENT(IN) :: nor1(nd,nd),jac(6,6) 
  DOUBLE PRECISION, INTENT(INOUT) :: nor2(nd,nd) ! to allow for nor1=nor2
  LOGICAL, INTENT(OUT) ::  error 
! END INTERFACE
  INTEGER ising 
  DOUBLE PRECISION jacinv(6,6),det 
! Inversion of Jacobian matrix                                          
  jacinv=jac
  CALL matin(jacinv,det,6,0,6,ising,1) 
  error=(ising.NE.0) 
  IF(error) THEN 
     nor2(1:nd,1:nd)=0.d0 
     RETURN 
  END IF
  CALL convertnordp(nd,nor1,jacinv(1:6,1:6),nor2)
END SUBROUTINE norprsdp

! version for nd x nd matrices        
SUBROUTINE norprs(nor1,jac,nd,nor2,error) 
  IMPLICIT NONE
  INTEGER, INTENT(IN) ::  nd 
  DOUBLE PRECISION, INTENT(IN) :: nor1(nd,nd),jac(nd,nd) 
  DOUBLE PRECISION, INTENT(INOUT) :: nor2(nd,nd) ! to allow for nor1=nor2
  LOGICAL, INTENT(OUT) ::  error 
! END INTERFACE
  INTEGER ising 
  DOUBLE PRECISION jacinv(nd,nd),det 
! Inversion of Jacobian matrix                                          
  jacinv=jac 
  CALL matin(jacinv,det,nd,0,nd,ising,1) 
  error=(ising.NE.0) 
  IF(error) THEN 
     nor2(1:nd,1:nd)=0.d0 
     RETURN 
  END IF
  CALL convertnor(nd,nor1,jacinv,nor2)
END SUBROUTINE norprs

                             
!  ***************************************************************      
!  LINFI3                                                               
!      linear fit to a real function, with residuals                    
!  ***************************************************************      
subroutine linfi3(tv,angv,rate,rms,cost,res,npo) 
  implicit none 
!  input                                                                
  integer npo 
  double precision tv(npo),angv(npo) 
!  output                                                               
  double precision res(npo),rate,rms,cost 
!  workspace                                                            
  double precision b(2,2),c(2),det,tmp,t,ang1,tm 
  integer i,j,n 
! ***************************************************************       
! compute the mean of tv                                                
  tm=0.d0 
  do n=1,npo 
     tm=tm+tv(n) 
  enddo
  tm=tm/npo 
! reduce to tm as origin                                             
  do 1 i=1,2 
     do 1 j=1,2 
1       b(i,j)=0.d0 
  do 2 i=1,2 
2    c(i)=0.d0 
  do 10 n=1,npo 
     t=tv(n)-tm 
     b(1,2)=b(1,2)+t 
     b(2,2)=b(2,2)+t**2 
     ang1=angv(n) 
     c(1)=c(1)+ang1 
     c(2)=c(2)+ang1*t 
10 continue 
! fit minimi quadrati: inversione matrice 2x2                           
  b(1,1)=npo 
  det=b(1,1)*b(2,2)-b(1,2)**2 
  b(2,1)=-b(1,2)/det 
  b(1,2)=b(2,1) 
  tmp=b(1,1)/det 
  b(1,1)=b(2,2)/det 
  b(2,2)=tmp 
  cost=b(1,1)*c(1)+b(1,2)*c(2) 
  rate=b(2,1)*c(1)+b(2,2)*c(2) 
! calcolo rms e residui                                                 
  rms=0.d0 
  do 5 i=1,npo 
     res(i)=angv(i)-rate*(tv(i)-tm)-cost 
     rms=rms+res(i)**2 
 5 continue 
   cost=cost-tm*rate 
   rms=dsqrt(rms/npo) 

 END SUBROUTINE linfi3

 SUBROUTINE cubic_fit(x,y,sy,m,ntime,g4,s4,rms_4,ising)
    IMPLICIT NONE
!INPUT:
    INTEGER,INTENT(IN) :: m, ntime ! number of data, 
            ! separate observation times
    DOUBLE PRECISION, INTENT(IN) :: x(m),y(m),sy(m) ! fit done
        ! on y, with RMS sy, as polynomial in x
!OUTPUT
    DOUBLE PRECISION, INTENT(OUT), DIMENSION(4) :: s4  ! cubic sol.
! covariance matrices
    DOUBLE PRECISION, INTENT(OUT) :: g4(4,4) !of cubic fit
    DOUBLE PRECISION, INTENT(OUT) :: rms_4 ! RMS of cubic fit
    INTEGER, INTENT(OUT) :: ising ! row of degeneracy
!END INTERFACE
    DOUBLE PRECISION :: sw,swxy(0:3),swy,sr,ww,detc
            !sums, weight, determinant, right hand side
    DOUBLE PRECISION :: c4(4,4) ! normal matr. of cubic fit
    DOUBLE PRECISION :: aa(4,5),det,swx(0:6),tmp 
            ! linear system for degree 2
    INTEGER j,k
! total weight and central time
    sw=0.d0
    swx=0.d0
! normal matrix and right hand side
    swxy=0.d0
    swy=0.d0 
    DO j=1,m
       ww=1/sy(j)**2
       sw=sw+ww
       DO k=0,6
          swx(k)=swx(k)+ww*x(j)**k
       ENDDO
       swy=swy+y(j)*ww
       DO k=0,3
          swxy(k)=swxy(k)+(x(j)**k)*y(j)*ww
       ENDDO
    ENDDO
    IF(m.gt.3.and.ntime.gt.3)THEN
       DO j=1,4
         DO k=1,4
           aa(j,k)=swx(8-k-j)
         ENDDO
         aa(j,5)=swxy(4-j)
       ENDDO
    ELSE
       STOP
    ENDIF
    CALL matin(aa,det,4,1,4,ising,1)
    g4=aa(1:4,1:4)
    s4=aa(1:4,5)
! post fit residuals
    sr=0.d0
    DO j=1,m
       ww=1/sy(j)**2
       tmp=ww*(y(j)-s4(1)*x(j)**3-s4(2)*x(j)**2-s4(3)*x(j)-s4(4))**2
       sr=sr+tmp
    ENDDO
    rms_4=sqrt(sr/m)
 END SUBROUTINE cubic_fit

 SUBROUTINE quadratic_fit(x,y,sy,m,ntime,g2,g3,s2,s3,rms_2,rms_3,ising)
    IMPLICIT NONE
!INPUT:
    INTEGER,INTENT(IN) :: m, ntime ! number of data, 
            ! separate observation times
    DOUBLE PRECISION, INTENT(IN) :: x(m),y(m),sy(m) ! fit done
        ! on y, with RMS sy, as polynomial in x
!OUTPUT
    DOUBLE PRECISION, INTENT(OUT), DIMENSION(2) :: s2  ! linear sol.
    DOUBLE PRECISION, INTENT(OUT), DIMENSION(3) :: s3  ! quadratic sol.
! covariance matrices
    DOUBLE PRECISION, INTENT(OUT) :: g2(2,2), g3(3,3) !of linear, quadratic fit
    DOUBLE PRECISION, INTENT(OUT) :: rms_2, rms_3 ! RMS of linear/quadr fit
    INTEGER, INTENT(OUT) :: ising ! row of degeneracy
!END INTERFACE
    DOUBLE PRECISION :: sw,swx,swx2,swxy,swy,sr,ww,detc,d2(2) 
            !sums, weight, determinant, right hand side
    DOUBLE PRECISION :: c2(2,2), c3(3,3) ! normal matr. of lin/quadr fit
    DOUBLE PRECISION :: aa(3,4),det,swx4,swx3,swx2y,tmp 
            ! linear system for degree 2
    INTEGER j
! total weight and central time
    sw=0.d0
    swx=0.d0
! normal matrix and right hand side
    swx2=0.d0
    swxy=0.d0
    swy=0.d0 
! for curvature
    swx3=0.d0
    swx4=0.d0
    swx2y=0.d0
    DO j=1,m
       ww=1/sy(j)**2
       sw=sw+ww
       swx=swx+x(j)*ww
       swx2=swx2+x(j)**2*ww
       swy=swy+y(j)*ww
       swxy=swxy+x(j)*y(j)*ww
! for curvature
       swx3=swx3+x(j)**3*ww
       swx4=swx4+x(j)**4*ww
       swx2y=swx2y+x(j)**2*y(j)*ww      
    ENDDO
!    IF(abs(swx).gt.100*epsilon(1.d0))WRITE(*,*)' central time trouble ',swx
! normal matrix for linear fit
    c2(1,1)=swx2
    c2(1,2)=swx
    c2(2,1)=swx
    c2(2,2)=sw
    d2(1)=swxy !right hand side for linear fit
    d2(2)=swy
    detc=c2(1,1)*c2(2,2)-c2(1,2)*c2(2,1)
    IF(abs(detc).lt.epsilon(detc)*c2(1,1)*100)THEN
       WRITE(*,*)' lin_reg: var(x) too small', detc,c2
    ENDIF
    g2(1,1)=c2(2,2)/detc ! covariance matrix for linear fit 
    g2(2,2)=c2(1,1)/detc
    g2(2,1)=-c2(2,1)/detc
    g2(1,2)=-c2(1,2)/detc
! solution
    s2=MATMUL(g2,d2)
! rms of residuals
    sr=0.d0
    DO j=1,m
       ww=1/sy(j)**2
       sr=sr+ww*(y(j)-s2(1)*x(j)-s2(2))**2
    ENDDO
    rms_2=sqrt(sr/m)
! check if quadratic fit is possible
    IF(m.gt.2.and.ntime.gt.2)THEN
       aa(1,1)=swx4
       aa(1,2)=swx3
       aa(2,1)=swx3
       aa(2,2)=swx2
       aa(3,1)=swx2
       aa(1,3)=swx2
       aa(3,2)=swx
       aa(2,3)=swx
       aa(3,3)=sw
       aa(1,4)=swx2y
       aa(2,4)=swxy
       aa(3,4)=swy
       c3=aa(1:3,1:3)
       CALL matin(aa,det,3,1,3,ising,1)
       g3=aa(1:3,1:3)
       s3=aa(1:3,4)
! post fit residuals
        sr=0.
        DO j=1,m
           ww=1/sy(j)**2
           tmp=ww*(y(j)-s3(1)*x(j)**2-s3(2)*x(j)-s3(3))**2
!              write(*,*)j,tmp
           sr=sr+tmp
        ENDDO
        rms_3=sqrt(sr/m)
    ELSE
       g3=0.d0
       s3=0.d0
       rms_3=0.d0
    ENDIF
  END SUBROUTINE quadratic_fit            
! =========================================================== 
! SORTING
! =========================================================== 
! ===================================                                   
! HEAPSORT                                                              
! sorting by heapsort of a double precision array a of lenght n         
! the input array a is not changed, in output                           
! ind is the indirect address of the sorted array                       
! ===================================                                   
SUBROUTINE heapsort(a,n,ind) 
  IMPLICIT NONE 
  INTEGER,INTENT(IN) ::  n 
  DOUBLE PRECISION,INTENT(IN) :: a(n) 
  INTEGER, INTENT(OUT) :: ind(n) 
! =====end interface========                                            
  integer j,l,ir,indt,i 
  double precision q 
! initialise indexes to natural order                                   
  DO j=1,n 
     ind(j)=j 
  ENDDO
  IF(n.eq.1)RETURN 
! counters for recursion, length of heap                                
  l=n/2+1 
  ir=n 
! recursive loop                                                        
1 CONTINUE 
  IF(l.gt.1)THEN 
     l=l-1 
     indt=ind(l) 
     q=a(indt) 
  ELSE 
     indt=ind(ir) 
     q=a(indt) 
     ind(ir)=ind(1) 
     ir=ir-1 
     IF(ir.eq.1)THEN 
        ind(1)=indt 
        RETURN 
     ENDIF
  ENDIF
  i=l 
  j=l+l 
2 IF(j.le.ir)THEN 
     IF(j.lt.ir)THEN 
        IF(a(ind(j)).lt.a(ind(j+1)))j=j+1 
     ENDIF
     IF(q.lt.a(ind(j)))THEN 
        ind(i)=ind(j) 
        i=j 
        j=j+j 
     ELSE 
        j=ir+1 
     ENDIF
     GOTO 2 
  ENDIF
  ind(i)=indt 
  GOTO 1 
END SUBROUTINE heapsort
! ===================================                                   
! HEAPSORTI                                                             
! sorting by heapsort of an integer array a of lenght n                 
! the input array a is not changed, in output                           
! ind is the indirect address of the sorted array                       
! ===================================                                   
SUBROUTINE heapsorti(a,n,ind) 
  IMPLICIT NONE 
  INTEGER, INTENT(IN) :: n 
  INTEGER, INTENT(IN) :: a(n) 
  INTEGER, INTENT(OUT) :: ind(n) 
! =====end interface========                                            
  integer j,l,ir,indt,i 
  double precision q 
! initialise indexes to natural order                                   
  DO j=1,n 
     ind(j)=j 
  ENDDO
  IF(n.eq.1)RETURN 
! counters for recursion, length of heap                                
  l=n/2+1 
  ir=n 
! recursive loop                                                        
1 CONTINUE 
  IF(l.gt.1)THEN 
     l=l-1 
     indt=ind(l) 
     q=a(indt) 
  ELSE 
     indt=ind(ir) 
     q=a(indt) 
     ind(ir)=ind(1) 
     ir=ir-1 
     IF(ir.eq.1)THEN 
        ind(1)=indt 
        RETURN 
     ENDIF
  ENDIF
  i=l 
  j=l+l 
2 IF(j.le.ir)THEN 
     IF(j.lt.ir)THEN 
        IF(a(ind(j)).lt.a(ind(j+1)))j=j+1 
     ENDIF
     IF(q.lt.a(ind(j)))THEN 
        ind(i)=ind(j) 
        i=j 
        j=j+j 
     ELSE 
        j=ir+1 
     ENDIF
     GOTO 2 
  ENDIF
  ind(i)=indt 
  GOTO 1 
END SUBROUTINE heapsorti
! =========================================================== 
! ROUNDING/TRUNCATING
! ===========================================================
! This function returns the truncated integer of the input.             
! Importantly it rounds to the nearest integer if the input             
! is within the given precision of an integer value                       
integer function truncat(flt,eps) 
  implicit none                                                              
  double precision, intent(in) :: flt,eps
  double precision one 
  one=1.d0-eps 
  truncat=flt 
  if(abs(truncat-flt).ge.one) truncat=nint(flt) 
END function truncat

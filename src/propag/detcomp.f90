! ==================================================================
! ************** M O D U L E   D E T C O M P ***********************
!           to compute the determinant of a matrix
! written by G.F.Gronchi, February 2007
! ==================================================================
! this module contains the following routines:
!  SUBROUTINE ZGETF2
!  SUBROUTINE ZGERU
!  subroutine  zscal
!  subroutine  zswap
!  SUBROUTINE XERBLA
! and the following functions:
!  integer function izamax
!  double precision function dcabs1
! ==================================================================
!MODULE detcomp
!  IMPLICIT NONE 
!  PRIVATE
!! subroutines
!  PUBLIC :: cdetcomp
!  PUBLIC :: ZGETF2,ZGERU,zscal,zswap,XERBLA 
!! functions
!  PUBLIC :: izamax,dcabs1
!CONTAINS

! ******************************************************************
! ******** computes the determinant of matrix A(6,6) employing *****
! ******** LU factorization (A(6,6) has complex coefficients) ******
! ******************************************************************
! *********** written by GIOVANNI F. GRONCHI (2001) ****************
! ********** Department of Mathematics, UNIVERSITY of PISA *********
! ==================================================================
  SUBROUTINE cdetcomp(A,detA) 
    IMPLICIT NONE 
! A is a matrix 6x6                                                 
    COMPLEX*16 :: A(6,6),detA 
! ============== end interface =================
    INTEGER :: INFO, LDA, M, N 
! IPIV must have the same dimension of N
    INTEGER :: IPIV(6)
    INTEGER :: signdetA 
! loop indexes                                                      
    INTEGER :: i,j 
! ==================================================================
    LDA= 6 
! number of columns of A                                            
    N=6 
! number of rows of A                                               
    M=6 
! initialization                                                    
    detA=1.d0 
    signdetA=1 
! compute LU factorization                                          
    CALL ZGETF2( M, N, A, LDA, IPIV, INFO ) 
    DO i=1,N 
       IF (IPIV(i).ne.i) THEN 
          signdetA = -signdetA 
       ENDIF
    ENDDO
! compute determinant of A                                          
    DO i=1,N 
       detA=detA*A(i,i) 
    ENDDO
    detA = signdetA*detA 
  END SUBROUTINE cdetcomp

! ==============================================================
  SUBROUTINE ZGETF2( M, N, A, LDA, IPIV, INFO ) 
!                                                                       
!  -- LAPACK routine (version 2.0) --                                   
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,       
!     Courant Institute, Argonne National Lab, and Rice University      
!     September 30, 1994                                                
!                                                                       
! .. Scalar Arguments ..                                            
    INTEGER            INFO, LDA, M, N 
! ..                                                                
! .. Array Arguments ..                                             
    INTEGER            IPIV( * ) 
    COMPLEX*16         A( LDA, * ) 
!  ..                                                                
!                                                                       
!  Purpose                                                              
!  =======                                                              
!                                                                       
!  ZGETF2 computes an LU factorization of a general m-by-n matrix A     
!  using partial pivoting with row interchanges.                        
!                                                                       
!  The factorization has the form                                       
!     A = P * L * U                                                     
!  where P is a permutation matrix, L is lower triangular with unit     
!  diagonal elements (lower trapezoidal if m > n), and U is upper       
!  triangular (upper trapezoidal if m < n).                             
!                                                                       
!  This is the right-looking Level 2 BLAS version of the algorithm.     
!                                                                       
!  Arguments                                                            
!  =========                                                            
!                                                                       
!  M       (input) INTEGER                                              
!          The number of rows of the matrix A.  M >= 0.                 
!                                                                       
!  N       (input) INTEGER                                              
!          The number of columns of the matrix A.  N >= 0.              
!                                                                       
!  A       (input/output) COMPLEX*16 array, dimension (LDA,N)           
!          On entry, the m by n matrix to be factored.                  
!          On exit, the factors L and U from the factorization          
!          A = P*L*U; the unit diagonal elements of L are not stored.   
!                                                                       
!  LDA     (input) INTEGER                                              
!          The leading dimension of the array A.  LDA >= max(1,M).      
!                                                                       
!  IPIV    (output) INTEGER array, dimension (min(M,N))                 
!          The pivot indices; for 1 <= i <= min(M,N), row i of the      
!          matrix was interchanged with row IPIV(i).                    
!                                                                       
!  INFO    (output) INTEGER                                             
!          = 0: successful exit                                         
!          < 0: if INFO = -k, the k-th argument had an illegal value    
!          > 0: if INFO = k, U(k,k) is exactly zero. The factorization  
!               has been completed, but the factor U is exactly         
!               singular, and division by zero will occur if it is used 
!               to solve a system of equations.                         
!                                                                       
!  =====================================================================
!                                                                       
!     .. Parameters ..                                                  
    COMPLEX*16         ONE, ZERO 
    PARAMETER          ( ONE = ( 1.0D+0, 0.0D+0 ),                    &
         &                   ZERO = ( 0.0D+0, 0.0D+0 ) )                    
!     ..                                                                
!     .. Local Scalars ..                                               
    INTEGER            J, JP 
!     ..                                                                
!     .. External Functions ..                                          
    INTEGER            IZAMAX 
    EXTERNAL           IZAMAX 
!     ..                                                                
!     .. External Subroutines ..                                        
    EXTERNAL           XERBLA, ZGERU, ZSCAL, ZSWAP 
!     ..                                                                
!     .. Intrinsic Functions ..                                         
    INTRINSIC          MAX, MIN 
!     ..                                                                
!     .. Executable Statements ..                                       
!                                                                       
!     Test the input parameters.                                        
!                                                                       
    INFO = 0 
    IF( M.LT.0 ) THEN 
       INFO = -1 
    ELSE IF( N.LT.0 ) THEN 
       INFO = -2 
    ELSE IF( LDA.LT.MAX( 1, M ) ) THEN 
       INFO = -4 
    END IF
    IF( INFO.NE.0 ) THEN 
       CALL XERBLA( 'ZGETF2', -INFO ) 
       RETURN 
    END IF
!                                                                       
!     Quick return if possible                                          
!                                                                       
    IF( M.EQ.0 .OR. N.EQ.0 )                                          &
         &   RETURN                                                         
!                                                                       
    DO 10 J = 1, MIN( M, N ) 
!                                                                       
!        Find pivot and test for singularity.                           
!                                                                       
       JP = J - 1 + IZAMAX( M-J+1, A( J, J ), 1 ) 
       IPIV( J ) = JP 
       IF( A( JP, J ).NE.ZERO ) THEN 
!                                                                       
!           Apply the interchange to columns 1:N.                       
!                                                                       
          IF( JP.NE.J )                                               &
               &         CALL ZSWAP( N, A( J, 1 ), LDA, A( JP, 1 ), LDA )         
!                                                                       
!           Compute elements J+1:M of J-th column.                      
!                                                                       
          IF( J.LT.M )                                                &
               &         CALL ZSCAL( M-J, ONE / A( J, J ), A( J+1, J ), 1 )       
!                                                                       
       ELSE IF( INFO.EQ.0 ) THEN 
!                                                                       
          INFO = J 
       END IF
!                                                                       
       IF( J.LT.MIN( M, N ) ) THEN 
!                                                                       
!           Update trailing submatrix.                                  
!                                                                       
          CALL ZGERU( M-J, N-J, -ONE, A( J+1, J ), 1, A( J, J+1 ),    &
               &                  LDA, A( J+1, J+1 ), LDA )                       
       END IF
10  END DO
    RETURN 
!                                                                       
!     End of ZGETF2                                                     
!                                                                       
  END SUBROUTINE ZGETF2

! ==========================================================
      SUBROUTINE ZGERU ( M, N, ALPHA, X, INCX, Y, INCY, A, LDA ) 
!     .. Scalar Arguments ..                                            
      COMPLEX*16         ALPHA 
      INTEGER            INCX, INCY, LDA, M, N 
!     .. Array Arguments ..                                             
      COMPLEX*16         A( LDA, * ), X( * ), Y( * ) 
!     ..                                                                
!                                                                       
!  Purpose                                                              
!  =======                                                              
!                                                                       
!  ZGERU  performs the rank 1 operation                                 
!                                                                       
!     A := alpha*x*y' + A,                                              
!                                                                       
!  where alpha is a scalar, x is an m element vector, y is an n element 
!  vector and A is an m by n matrix.                                    
!                                                                       
!  Parameters                                                           
!  ==========                                                           
!                                                                       
!  M      - INTEGER.                                                    
!           On entry, M specifies the number of rows of the matrix A.   
!           M must be at least zero.                                    
!           Unchanged on exit.                                          
!                                                                       
!  N      - INTEGER.                                                    
!           On entry, N specifies the number of columns of the matrix A.
!           N must be at least zero.                                    
!           Unchanged on exit.                                          
!                                                                       
!  ALPHA  - COMPLEX*16      .                                           
!           On entry, ALPHA specifies the scalar alpha.                 
!           Unchanged on exit.                                          
!                                                                       
!  X      - COMPLEX*16       array of dimension at least                
!           ( 1 + ( m - 1 )*abs( INCX ) ).                              
!           Before entry, the incremented array X must contain the m    
!           element vector x.                                           
!           Unchanged on exit.                                          
!                                                                       
!  INCX   - INTEGER.                                                    
!           On entry, INCX specifies the increment for the elements of  
!           X. INCX must not be zero.                                   
!           Unchanged on exit.                                          
!                                                                       
!  Y      - COMPLEX*16       array of dimension at least                
!           ( 1 + ( n - 1 )*abs( INCY ) ).                              
!           Before entry, the incremented array Y must contain the n    
!           element vector y.                                           
!           Unchanged on exit.                                          
!                                                                       
!  INCY   - INTEGER.                                                    
!           On entry, INCY specifies the increment for the elements of  
!           Y. INCY must not be zero.                                   
!           Unchanged on exit.                                          
!                                                                       
!  A      - COMPLEX*16       array of DIMENSION ( LDA, n ).             
!           Before entry, the leading m by n part of the array A must   
!           contain the matrix of coefficients. On exit, A is           
!           overwritten by the updated matrix.                          
!                                                                       
!  LDA    - INTEGER.                                                    
!           On entry, LDA specifies the first dimension of A as declared
!           in the calling (sub) program. LDA must be at least          
!           max( 1, m ).                                                
!           Unchanged on exit.                                          
!                                                                       
!                                                                       
!  Level 2 Blas routine.                                                
!                                                                       
!  -- Written on 22-October-1986.                                       
!     Jack Dongarra, Argonne National Lab.                              
!     Jeremy Du Croz, Nag Central Office.                               
!     Sven Hammarling, Nag Central Office.                              
!     Richard Hanson, Sandia National Labs.                             
!                                                                       
!                                                                       
!     .. Parameters ..                                                  
      COMPLEX*16         ZERO 
      PARAMETER        ( ZERO = ( 0.0D+0, 0.0D+0 ) ) 
!     .. Local Scalars ..                                               
      COMPLEX*16         TEMP 
      INTEGER            I, INFO, IX, J, JY, KX 
!     .. External Subroutines ..                                        
      EXTERNAL           XERBLA 
!     .. Intrinsic Functions ..                                         
      INTRINSIC          MAX 
!     ..                                                                
!     .. Executable Statements ..                                       
!                                                                       
!     Test the input parameters.                                        
!                                                                       
      INFO = 0 
      IF     ( M.LT.0 )THEN 
         INFO = 1 
      ELSE IF( N.LT.0 )THEN 
         INFO = 2 
      ELSE IF( INCX.EQ.0 )THEN 
         INFO = 5 
      ELSE IF( INCY.EQ.0 )THEN 
         INFO = 7 
      ELSE IF( LDA.LT.MAX( 1, M ) )THEN 
         INFO = 9 
      END IF 
      IF( INFO.NE.0 )THEN 
         CALL XERBLA( 'ZGERU ', INFO ) 
         RETURN 
      END IF 
!                                                                       
!     Quick return if possible.                                         
!                                                                       
      IF( ( M.EQ.0 ).OR.( N.EQ.0 ).OR.( ALPHA.EQ.ZERO ) )               &
     &   RETURN                                                         
!                                                                       
!     Start the operations. In this version the elements of A are       
!     accessed sequentially with one pass through A.                    
!                                                                       
      IF( INCY.GT.0 )THEN 
         JY = 1 
      ELSE 
         JY = 1 - ( N - 1 )*INCY 
      END IF 
      IF( INCX.EQ.1 )THEN 
         DO 20, J = 1, N 
            IF( Y( JY ).NE.ZERO )THEN 
               TEMP = ALPHA*Y( JY ) 
               DO 10, I = 1, M 
                  A( I, J ) = A( I, J ) + X( I )*TEMP 
   10          CONTINUE 
            END IF 
            JY = JY + INCY 
   20    CONTINUE 
      ELSE 
         IF( INCX.GT.0 )THEN 
            KX = 1 
         ELSE 
            KX = 1 - ( M - 1 )*INCX 
         END IF 
         DO 40, J = 1, N 
            IF( Y( JY ).NE.ZERO )THEN 
               TEMP = ALPHA*Y( JY ) 
               IX   = KX 
               DO 30, I = 1, M 
                  A( I, J ) = A( I, J ) + X( IX )*TEMP 
                  IX        = IX        + INCX 
   30          CONTINUE 
            END IF 
            JY = JY + INCY 
   40    CONTINUE 
      END IF 
!                                                                       
      RETURN 
!                                                                       
!     End of ZGERU .                                                    
!                                                                       
      END SUBROUTINE ZGERU 

! ===============================================================
      subroutine  zscal(n,za,zx,incx) 
!                                                                       
!     scales a vector by a constant.                                    
!     jack dongarra, 3/11/78.                                           
!     modified 3/93 to return if incx .le. 0.                           
!     modified 12/3/93, array(1) declarations changed to array(*)       
!                                                                       
      complex*16 za,zx(*) 
      integer i,incx,ix,n 
!                                                                       
      if( n.le.0 .or. incx.le.0 )return 
      if(incx.eq.1)go to 20 
!                                                                       
!        code for increment not equal to 1                              
!                                                                       
      ix = 1 
      do 10 i = 1,n 
        zx(ix) = za*zx(ix) 
        ix = ix + incx 
   10 continue 
      return 
!                                                                       
!        code for increment equal to 1                                  
!                                                                       
   20 do 30 i = 1,n 
        zx(i) = za*zx(i) 
   30 continue 
      return 
      END subroutine zscal

! ===============================================================
      subroutine  zswap (n,zx,incx,zy,incy) 
!                                                                       
!     interchanges two vectors.                                         
!     jack dongarra, 3/11/78.                                           
!     modified 12/3/93, array(1) declarations changed to array(*)       
!                                                                       
      complex*16 zx(*),zy(*),ztemp 
      integer i,incx,incy,ix,iy,n 
!                                                                       
      if(n.le.0)return 
      if(incx.eq.1.and.incy.eq.1)go to 20 
!                                                                       
!       code for unequal increments or equal increments not equal       
!         to 1                                                          
!                                                                       
      ix = 1 
      iy = 1 
      if(incx.lt.0)ix = (-n+1)*incx + 1 
      if(incy.lt.0)iy = (-n+1)*incy + 1 
      do 10 i = 1,n 
        ztemp = zx(ix) 
        zx(ix) = zy(iy) 
        zy(iy) = ztemp 
        ix = ix + incx 
        iy = iy + incy 
   10 continue 
      return 
!                                                                       
!       code for both increments equal to 1                             
   20 do 30 i = 1,n 
        ztemp = zx(i) 
        zx(i) = zy(i) 
        zy(i) = ztemp 
   30 continue 
      return 
      END subroutine zswap                                         

! ===============================================================
      SUBROUTINE XERBLA( SRNAME, INFO ) 
!                                                                       
!  -- LAPACK auxiliary routine (version 2.0) --                         
!     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,       
!     Courant Institute, Argonne National Lab, and Rice University      
!     September 30, 1994                                                
!                                                                       
!     .. Scalar Arguments ..                                            
      CHARACTER*6        SRNAME 
      INTEGER            INFO 
!     ..                                                                
!                                                                       
!  Purpose                                                              
!  =======                                                              
!                                                                       
!  XERBLA  is an error handler for the LAPACK routines.                 
!  It is called by an LAPACK routine if an input parameter has an       
!  invalid value.  A message is printed and execution stops.            
!                                                                       
!  Installers may consider modifying the STOP statement in order to     
!  call system-specific exception-handling facilities.                  
!                                                                       
!  Arguments                                                            
!  =========                                                            
!                                                                       
!  SRNAME  (input) CHARACTER*6                                          
!          The name of the routine which called XERBLA.                 
!                                                                       
!  INFO    (input) INTEGER                                              
!          The position of the invalid parameter in the parameter list  
!          of the calling routine.                                      
!                                                                       
! ===================================================================== 
!                                                                       
!     .. Executable Statements ..                                       
!                                                                       
      WRITE( *, FMT = 9999 )SRNAME, INFO 
!                                                                       
      STOP 
!                                                                       
 9999 FORMAT( ' ** On entry to ', A6, ' parameter number ', I2, ' had ',&
     &      'an illegal value' )                                        
!                                                                       
!     End of XERBLA                                                     
!                                                                       
      END SUBROUTINE XERBLA

! ===============================================================
      integer function izamax(n,zx,incx) 
!                                                                       
!     finds the index of element having max. absolute value.            
!     jack dongarra, 1/15/85.                                           
!     modified 3/93 to return if incx .le. 0.                           
!     modified 12/3/93, array(1) declarations changed to array(*)       
!                                                                       
      complex*16 zx(*) 
      double precision smax 
      integer i,incx,ix,n 
      double precision dcabs1 
!                                                                       
      izamax = 0 
      if( n.lt.1 .or. incx.le.0 )return 
      izamax = 1 
      if(n.eq.1)return 
      if(incx.eq.1)go to 20 
!                                                                       
!        code for increment not equal to 1                              
!                                                                       
      ix = 1 
      smax = dcabs1(zx(1)) 
      ix = ix + incx 
      do 10 i = 2,n 
         if(dcabs1(zx(ix)).le.smax) go to 5 
         izamax = i 
         smax = dcabs1(zx(ix)) 
    5    ix = ix + incx 
   10 continue 
      return 
!                                                                       
!        code for increment equal to 1                                  
!                                                                       
   20 smax = dcabs1(zx(1)) 
      do 30 i = 2,n 
         if(dcabs1(zx(i)).le.smax) go to 30 
         izamax = i 
         smax = dcabs1(zx(i)) 
   30 continue 
      return 
      END function izamax

! ===============================================================
      double precision function dcabs1(z) 
       complex*16 z 
!      complex*16 z,zz                                                  
!      double precision t(2)                                            
!      equivalence (zz,t(1))                                            
!      zz = z                                                           
!      dcabs1 = dabs(t(1)) + dabs(t(2))                                 
!       dcabs1=abs(real(z))+abs(real(z*(0,1)))                          
       dcabs1=abs(DBLE(z))+abs(DBLE(z*(0,1))) 
                                                                        
       return 
      END function dcabs1                                         

!END MODULE detcomp

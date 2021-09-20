      subroutine rs(nm,n,a,w,matz,z,fv1,fv2,ierr) 
!                                                                       
      integer n,nm,ierr,matz 
      double precision a(nm,n),w(n),z(nm,n),fv1(n),fv2(n) 
!                                                                       
!     this subroutine calls the recommended sequence of                 
!     subroutines from the eigensystem subroutine package (eispack)     
!     to find the eigenvalues and eigenvectors (if desired)             
!     of a real symmetric matrix.                                       
!                                                                       
!     on input:                                                         
!                                                                       
!        nm  must be set to the row dimension of the two-dimensional    
!        array parameters as declared in the calling program            
!        dimension statement;                                           
!                                                                       
!        n  is the order of the matrix  a;                              
!                                                                       
!        a  contains the real symmetric matrix;                         
!                                                                       
!        matz  is an integer variable set equal to zero if              
!        only eigenvalues are desired;  otherwise it is set to          
!        any non-zero integer for both eigenvalues and eigenvectors.    
!                                                                       
!     on output:                                                        
!                                                                       
!        w  contains the eigenvalues in ascending order;                
!                                                                       
!        z  contains the eigenvectors if matz is not zero;              
!                                                                       
!        ierr  is an integer output variable set equal to an            
!        error completion code described in section 2b of the           
!        documentation.  the normal completion code is zero;            
!                                                                       
!        fv1  and  fv2  are temporary storage arrays.                   
!                                                                       
!     questions and comments should be directed to b. s. garbow,        
!     applied mathematics division, argonne national laboratory         
!                                                                       
!     ------------------------------------------------------------------
!                                                                       
      if (n .le. nm) go to 10 
      ierr = 10 * n 
      go to 50 
!                                                                       
   10 if (matz .ne. 0) go to 20 
!     :::::::::: find eigenvalues only ::::::::::                       
      call  tred1(nm,n,a,w,fv1,fv2) 
      call  tqlrat(n,w,fv2,ierr) 
      go to 50 
!     :::::::::: find both eigenvalues and eigenvectors ::::::::::      
   20 call  tred2(nm,n,a,w,fv1,z) 
      call  tql2(nm,n,w,fv1,z,ierr) 
   50 return 
!     :::::::::: last card of rs ::::::::::                             
      END                                           
      subroutine tred1(nm,n,a,d,e,e2) 
!                                                                       
      integer i,j,k,l,n,ii,nm,jp1 
      double precision a(nm,n),d(n),e(n),e2(n) 
      double precision f,g,h,scale 
!                                                                       
!     this subroutine is a translation of the algol procedure tred1,    
!     num. math. 11, 181-195(1968) by martin, reinsch, and wilkinson.   
!     handbook for auto. comp., vol.ii-linear algebra, 212-226(1971).   
!                                                                       
!     this subroutine reduces a real symmetric matrix                   
!     to a symmetric tridiagonal matrix using                           
!     orthogonal similarity transformations.                            
!                                                                       
!     on input:                                                         
!                                                                       
!        nm must be set to the row dimension of two-dimensional         
!          array parameters as declared in the calling program          
!          dimension statement;                                         
!                                                                       
!        n is the order of the matrix;                                  
!                                                                       
!        a contains the real symmetric input matrix.  only the          
!          lower triangle of the matrix need be supplied.               
!                                                                       
!     on output:                                                        
!                                                                       
!        a contains information about the orthogonal trans-             
!          formations used in the reduction in its strict lower         
!          triangle.  the full upper triangle of a is unaltered;        
!                                                                       
!        d contains the diagonal elements of the tridiagonal matrix;    
!                                                                       
!        e contains the subdiagonal elements of the tridiagonal         
!          matrix in its last n-1 positions.  e(1) is set to zero;      
!                                                                       
!        e2 contains the squares of the corresponding elements of e.    
!          e2 may coincide with e if the squares are not needed.        
!                                                                       
!     questions and comments should be directed to b. s. garbow,        
!     applied mathematics division, argonne national laboratory         
!                                                                       
!     ------------------------------------------------------------------
!                                                                       
      do 100 i = 1, n 
  100 d(i) = a(i,i) 
!     :::::::::: for i=n step -1 until 1 do -- ::::::::::               
      do 300 ii = 1, n 
         i = n + 1 - ii 
         l = i - 1 
         h = 0.0d0 
         scale = 0.0d0 
         if (l .lt. 1) go to 130 
!     :::::::::: scale row (algol tol then not needed) ::::::::::       
         do 120 k = 1, l 
  120    scale = scale + abs(a(i,k)) 
!                                                                       
         if (scale .ne. 0.0d0) go to 140 
  130    e(i) = 0.0d0 
         e2(i) = 0.0d0 
         go to 290 
!                                                                       
  140    do 150 k = 1, l 
            a(i,k) = a(i,k) / scale 
            h = h + a(i,k) * a(i,k) 
  150    continue 
!                                                                       
         e2(i) = scale * scale * h 
         f = a(i,l) 
         g = -sign(sqrt(h),f) 
         e(i) = scale * g 
         h = h - f * g 
         a(i,l) = f - g 
         if (l .eq. 1) go to 270 
         f = 0.0d0 
!                                                                       
         do 240 j = 1, l 
            g = 0.0d0 
!     :::::::::: form element of a*u ::::::::::                         
            do 180 k = 1, j 
  180       g = g + a(j,k) * a(i,k) 
!                                                                       
            jp1 = j + 1 
            if (l .lt. jp1) go to 220 
!                                                                       
            do 200 k = jp1, l 
  200       g = g + a(k,j) * a(i,k) 
!     :::::::::: form element of p ::::::::::                           
  220       e(j) = g / h 
            f = f + e(j) * a(i,j) 
  240    continue 
!                                                                       
         h = f / (h + h) 
!     :::::::::: form reduced a ::::::::::                              
         do 260 j = 1, l 
            f = a(i,j) 
            g = e(j) - h * f 
            e(j) = g 
!                                                                       
            do 260 k = 1, j 
               a(j,k) = a(j,k) - f * e(k) - g * a(i,k) 
  260    continue 
!                                                                       
  270    do 280 k = 1, l 
  280    a(i,k) = scale * a(i,k) 
!                                                                       
  290    h = d(i) 
         d(i) = a(i,i) 
         a(i,i) = h 
  300 continue 
!                                                                       
      return 
!     :::::::::: last card of tred1 ::::::::::                          
      END                                           
      subroutine tqlrat(n,d,e2,ierr) 
!                                                                       
      integer i,j,l,m,n,ii,l1,mml,ierr 
      double precision d(n),e2(n) 
      double precision b,c,f,g,h,p,r,s,machep 
!      real*8 dsqrt,dabs,dsign                                          
!                                                                       
!     this subroutine is a translation of the algol procedure tqlrat,   
!     algorithm 464, comm. acm 16, 689(1973) by reinsch.                
!                                                                       
!     this subroutine finds the eigenvalues of a symmetric              
!     tridiagonal matrix by the rational ql method.                     
!                                                                       
!     on input:                                                         
!                                                                       
!        n is the order of the matrix;                                  
!                                                                       
!        d contains the diagonal elements of the input matrix;          
!                                                                       
!        e2 contains the squares of the subdiagonal elements of the     
!          input matrix in its last n-1 positions.  e2(1) is arbitrary. 
!                                                                       
!      on output:                                                       
!                                                                       
!        d contains the eigenvalues in ascending order.  if an          
!          error exit is made, the eigenvalues are correct and          
!          ordered for indices 1,2,...ierr-1, but may not be            
!          the smallest eigenvalues;                                    
!                                                                       
!        e2 has been destroyed;                                         
!                                                                       
!        ierr is set to                                                 
!          zero       for normal return,                                
!          j          if the j-th eigenvalue has not been               
!                     determined after 30 iterations.                   
!                                                                       
!     questions and comments should be directed to b. s. garbow,        
!     applied mathematics division, argonne national laboratory         
!                                                                       
!     ------------------------------------------------------------------
!                                                                       
!     :::::::::: machep is a machine dependent parameter specifying     
!                the relative precision of floating point arithmetic.   
!                machep = 16.0d0**(-13) for long form arithmetic        
!                on s360 ::::::::::                                     
!      data machep/z'3410000000000000'/                                 
      data machep/1.d-16/ 
!                                                                       
      ierr = 0 
      if (n .eq. 1) go to 1001 
!                                                                       
      do 100 i = 2, n 
  100 e2(i-1) = e2(i) 
!                                                                       
      f = 0.0d0 
      b = 0.0d0 
      e2(n) = 0.0d0 
!                                                                       
      do 290 l = 1, n 
         j = 0 
         h = machep * (abs(d(l)) + sqrt(e2(l))) 
         if (b .gt. h) go to 105 
         b = h 
         c = b * b 
!     :::::::::: look for small squared sub-diagonal element :::::::::: 
  105    do 110 m = l, n 
            if (e2(m) .le. c) go to 120 
!     :::::::::: e2(n) is always zero, so there is no exit              
!                through the bottom of the loop ::::::::::              
  110    continue 
!                                                                       
  120    if (m .eq. l) go to 210 
  130    if (j .eq. 30) go to 1000 
         j = j + 1 
!     :::::::::: form shift ::::::::::                                  
         l1 = l + 1 
         s = sqrt(e2(l)) 
         g = d(l) 
         p = (d(l1) - g) / (2.0d0 * s) 
         r = sqrt(p*p+1.0d0) 
         d(l) = s / (p + sign(r,p)) 
         h = g - d(l) 
!                                                                       
         do 140 i = l1, n 
  140    d(i) = d(i) - h 
!                                                                       
         f = f + h 
!     :::::::::: rational ql transformation ::::::::::                  
         g = d(m) 
         if (g .eq. 0.0d0) g = b 
         h = g 
         s = 0.0d0 
         mml = m - l 
!     :::::::::: for i=m-1 step -1 until l do -- ::::::::::             
         do 200 ii = 1, mml 
            i = m - ii 
            p = g * h 
            r = p + e2(i) 
            e2(i+1) = s * r 
            s = e2(i) / r 
            d(i+1) = h + s * (h + d(i)) 
            g = d(i) - e2(i) / g 
            if (g .eq. 0.0d0) g = b 
            h = g * p / r 
  200    continue 
!                                                                       
         e2(l) = s * g 
         d(l) = h 
!     :::::::::: guard against underflow in convergence test :::::::::: 
         if (h .eq. 0.0d0) go to 210 
         if (abs(e2(l)) .le. abs(c/h)) go to 210 
         e2(l) = h * e2(l) 
         if (e2(l) .ne. 0.0d0) go to 130 
  210    p = d(l) + f 
!     :::::::::: order eigenvalues ::::::::::                           
         if (l .eq. 1) go to 250 
!     :::::::::: for i=l step -1 until 2 do -- ::::::::::               
         do 230 ii = 2, l 
            i = l + 2 - ii 
            if (p .ge. d(i-1)) go to 270 
            d(i) = d(i-1) 
  230    continue 
!                                                                       
  250    i = 1 
  270    d(i) = p 
  290 continue 
!                                                                       
      go to 1001 
!     :::::::::: set error -- no convergence to an                      
!                eigenvalue after 30 iterations ::::::::::              
 1000 ierr = l 
 1001 return 
!     :::::::::: last card of tqlrat ::::::::::                         
      END                                           
      subroutine tred2(nm,n,a,d,e,z) 
!                                                                       
      integer i,j,k,l,n,ii,nm,jp1 
      double precision a(nm,n),d(n),e(n),z(nm,n) 
      double precision f,g,h,hh,scale 
!      real*8 dsqrt,dabs,dsign                                          
!                                                                       
!     this subroutine is a translation of the algol procedure tred2,    
!     num. math. 11, 181-195(1968) by martin, reinsch, and wilkinson.   
!     handbook for auto. comp., vol.ii-linear algebra, 212-226(1971).   
!                                                                       
!     this subroutine reduces a real symmetric matrix to a              
!     symmetric tridiagonal matrix using and accumulating               
!     orthogonal similarity transformations.                            
!                                                                       
!     on input:                                                         
!                                                                       
!        nm must be set to the row dimension of two-dimensional         
!          array parameters as declared in the calling program          
!          dimension statement;                                         
!                                                                       
!        n is the order of the matrix;                                  
!                                                                       
!        a contains the real symmetric input matrix.  only the          
!          lower triangle of the matrix need be supplied.               
!                                                                       
!     on output:                                                        
!                                                                       
!        d contains the diagonal elements of the tridiagonal matrix;    
!                                                                       
!        e contains the subdiagonal elements of the tridiagonal         
!          matrix in its last n-1 positions.  e(1) is set to zero;      
!                                                                       
!        z contains the orthogonal transformation matrix                
!          produced in the reduction;                                   
!                                                                       
!        a and z may coincide.  if distinct, a is unaltered.            
!                                                                       
!     questions and comments should be directed to b. s. garbow,        
!     applied mathematics division, argonne national laboratory         
!                                                                       
!     ------------------------------------------------------------------
!                                                                       
      do 100 i = 1, n 
!                                                                       
         do 100 j = 1, i 
            z(i,j) = a(i,j) 
  100 continue 
!                                                                       
      if (n .eq. 1) go to 320 
!     :::::::::: for i=n step -1 until 2 do -- ::::::::::               
      do 300 ii = 2, n 
         i = n + 2 - ii 
         l = i - 1 
         h = 0.0d0 
         scale = 0.0d0 
         if (l .lt. 2) go to 130 
!     :::::::::: scale row (algol tol then not needed) ::::::::::       
         do 120 k = 1, l 
  120    scale = scale + abs(z(i,k)) 
!                                                                       
         if (scale .ne. 0.0d0) go to 140 
  130    e(i) = z(i,l) 
         go to 290 
!                                                                       
  140    do 150 k = 1, l 
            z(i,k) = z(i,k) / scale 
            h = h + z(i,k) * z(i,k) 
  150    continue 
!                                                                       
         f = z(i,l) 
         g = -sign(sqrt(h),f) 
         e(i) = scale * g 
         h = h - f * g 
         z(i,l) = f - g 
         f = 0.0d0 
!                                                                       
         do 240 j = 1, l 
            z(j,i) = z(i,j) / h 
            g = 0.0d0 
!     :::::::::: form element of a*u ::::::::::                         
            do 180 k = 1, j 
  180       g = g + z(j,k) * z(i,k) 
!                                                                       
            jp1 = j + 1 
            if (l .lt. jp1) go to 220 
!                                                                       
            do 200 k = jp1, l 
  200       g = g + z(k,j) * z(i,k) 
!     :::::::::: form element of p ::::::::::                           
  220       e(j) = g / h 
            f = f + e(j) * z(i,j) 
  240    continue 
!                                                                       
         hh = f / (h + h) 
!     :::::::::: form reduced a ::::::::::                              
         do 260 j = 1, l 
            f = z(i,j) 
            g = e(j) - hh * f 
            e(j) = g 
!                                                                       
            do 260 k = 1, j 
               z(j,k) = z(j,k) - f * e(k) - g * z(i,k) 
  260    continue 
!                                                                       
  290    d(i) = h 
  300 continue 
!                                                                       
  320 d(1) = 0.0d0 
      e(1) = 0.0d0 
!     :::::::::: accumulation of transformation matrices ::::::::::     
      do 500 i = 1, n 
         l = i - 1 
         if (d(i) .eq. 0.0d0) go to 380 
!                                                                       
         do 360 j = 1, l 
            g = 0.0d0 
!                                                                       
            do 340 k = 1, l 
  340       g = g + z(i,k) * z(k,j) 
!                                                                       
            do 360 k = 1, l 
               z(k,j) = z(k,j) - g * z(k,i) 
  360    continue 
!                                                                       
  380    d(i) = z(i,i) 
         z(i,i) = 1.0d0 
         if (l .lt. 1) go to 500 
!                                                                       
         do 400 j = 1, l 
            z(i,j) = 0.0d0 
            z(j,i) = 0.0d0 
  400    continue 
!                                                                       
  500 continue 
!                                                                       
      return 
!     :::::::::: last card of tred2 ::::::::::                          
      END                                           
      subroutine tql2(nm,n,d,e,z,ierr) 
!                                                                       
      integer i,j,k,l,m,n,ii,l1,nm,mml,ierr 
      double precision d(n),e(n),z(nm,n) 
      double precision b,c,f,g,h,p,r,s,machep 
!      real*8 dsqrt,dabs,dsign                                          
!                                                                       
!     this subroutine is a translation of the algol procedure tql2,     
!     num. math. 11, 293-306(1968) by bowdler, martin, reinsch, and     
!     wilkinson.                                                        
!     handbook for auto. comp., vol.ii-linear algebra, 227-240(1971).   
!                                                                       
!     this subroutine finds the eigenvalues and eigenvectors            
!     of a symmetric tridiagonal matrix by the ql method.               
!     the eigenvectors of a full symmetric matrix can also              
!     be found if  tred2  has been used to reduce this                  
!     full matrix to tridiagonal form.                                  
!                                                                       
!     on input:                                                         
!                                                                       
!        nm must be set to the row dimension of two-dimensional         
!          array parameters as declared in the calling program          
!          dimension statement;                                         
!                                                                       
!        n is the order of the matrix;                                  
!                                                                       
!        d contains the diagonal elements of the input matrix;          
!                                                                       
!        e contains the subdiagonal elements of the input matrix        
!          in its last n-1 positions.  e(1) is arbitrary;               
!                                                                       
!        z contains the transformation matrix produced in the           
!          reduction by  tred2, if performed.  if the eigenvectors      
!          of the tridiagonal matrix are desired, z must contain        
!          the identity matrix.                                         
!                                                                       
!      on output:                                                       
!                                                                       
!        d contains the eigenvalues in ascending order.  if an          
!          error exit is made, the eigenvalues are correct but          
!          unordered for indices 1,2,...,ierr-1;                        
!                                                                       
!        e has been destroyed;                                          
!                                                                       
!        z contains orthonormal eigenvectors of the symmetric           
!          tridiagonal (or full) matrix.  if an error exit is made,     
!          z contains the eigenvectors associated with the stored       
!          eigenvalues;                                                 
!                                                                       
!        ierr is set to                                                 
!          zero       for normal return,                                
!          j          if the j-th eigenvalue has not been               
!                     determined after 30 iterations.                   
!                                                                       
!     questions and comments should be directed to b. s. garbow,        
!     applied mathematics division, argonne national laboratory         
!                                                                       
!     ------------------------------------------------------------------
!                                                                       
!     :::::::::: machep is a machine dependent parameter specifying     
!                the relative precision of floating point arithmetic.   
!                machep = 16.0d0**(-13) for long form arithmetic        
!                on s360 ::::::::::                                     
!      data machep/z'3410000000000000'/                                 
!      data machep/1.d-16/ 
      machep=epsilon(1.d0)
!
                                                                             
      ierr = 0 
      if (n .eq. 1) go to 1001 
!                                                                       
      do 100 i = 2, n 
  100 e(i-1) = e(i) 
!                                                                       
      f = 0.0d0 
      b = 0.0d0 
      e(n) = 0.0d0 
!                                                                       
      do 240 l = 1, n 
         j = 0 
         h = machep * (dabs(d(l)) + dabs(e(l))) 
         if (b .lt. h) b = h 
!     :::::::::: look for small sub-diagonal element ::::::::::         
         do 110 m = l, n 
            if (dabs(e(m)) .le. b) go to 120 
!     :::::::::: e(n) is always zero, so there is no exit               
!                through the bottom of the loop ::::::::::              
  110    continue 
!
         WRITE(*,*) ' tql2: WRONG! m,n=',m,n
  120    if (m .eq. l) go to 220 
  130    if (j .eq. 30) go to 1000 
         j = j + 1 
!     :::::::::: form shift ::::::::::                                  
         l1 = l + 1 
         g = d(l) 
         p = (d(l1) - g) / (2.0d0 * e(l)) 
         r = dsqrt(p*p+1.0d0) 
         d(l) = e(l) / (p + dsign(r,p)) 
         h = g - d(l) 
!                                                                       
         do 140 i = l1, n 
  140    d(i) = d(i) - h 
!                                                                       
         f = f + h 
!     :::::::::: ql transformation ::::::::::                           
         IF(m.eq.0)THEN
            WRITE(*,*) 'tql2: ',nm,n,m,l
         ENDIF
         p = d(m) 
         c = 1.0d0 
         s = 0.0d0 
         mml = m - l 
!     :::::::::: for i=m-1 step -1 until l do -- ::::::::::             
         do 200 ii = 1, mml 
            i = m - ii 
            g = c * e(i) 
            h = c * p 
            if (dabs(p) .lt. dabs(e(i))) go to 150 
            c = e(i) / p 
            r = dsqrt(c*c+1.0d0) 
            e(i+1) = s * p * r 
            s = c / r 
            c = 1.0d0 / r 
            go to 160 
  150       c = p / e(i) 
            r = dsqrt(c*c+1.0d0) 
            e(i+1) = s * e(i) * r 
            s = 1.0d0 / r 
            c = c * s 
  160       p = c * d(i) - s * g 
            d(i+1) = h + s * (c * g + s * d(i)) 
!     :::::::::: form vector ::::::::::                                 
            do 180 k = 1, n 
               h = z(k,i+1) 
               z(k,i+1) = s * z(k,i) + c * h 
               z(k,i) = c * z(k,i) - s * h 
  180       continue 
!                                                                       
  200    continue 
!                                                                       
         e(l) = s * p 
         d(l) = c * p 
         if (dabs(e(l)) .gt. b) go to 130 
  220    d(l) = d(l) + f 
  240 continue 
!     :::::::::: order eigenvalues and eigenvectors ::::::::::          
      do 300 ii = 2, n 
         i = ii - 1 
         k = i 
         p = d(i) 
!                                                                       
         do 260 j = ii, n 
            if (d(j) .ge. p) go to 260 
            k = j 
            p = d(j) 
  260    continue 
!                                                                       
         if (k .eq. i) go to 300 
         d(k) = d(i) 
         d(i) = p 
!                                                                       
         do 280 j = 1, n 
            p = z(j,i) 
            z(j,i) = z(j,k) 
            z(j,k) = p 
  280    continue 
!                                                                       
  300 continue 
!                                                                       
      go to 1001 
!     :::::::::: set error -- no convergence to an                      
!                eigenvalue after 30 iterations ::::::::::              
 1000 ierr = l 
 1001 return 
!     :::::::::: last card of tql2 ::::::::::                           
      END                                           

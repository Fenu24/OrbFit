! Slightly modified by M.Carpino (December 14, 1998): I have only change
! the name of SUBROUTINE START into SUBROUTINE PZSTART, in order to     
! avoid a name clash with propag/start.f                                
!***********************************************************************
!    NUMERICAL COMPUTATION OF THE ROOTS OF A POLYNOMIAL HAVING          
!        COMPLEX COEFFICIENTS, BASED ON ABERTH'S METHOD.                
!                      Version 1.4, June   1996                         
!    (D. Bini, Dipartimento di Matematica, Universita' di Pisa)         
!                         (bini@dm.unipi.it)                            
!***********************************************************************
! Work performed under the support of the ESPRIT BRA project 6846 POSSO 
!***********************************************************************
!**********         SUBROUTINES AND FUNCTIONS                 **********
!***********************************************************************
!  The following modules are listed:                                    
!  POLZEROS  :  computes polynomial roots by means of Aberth's method   
!    ABERTH  :  computes the Aberth correction                          
!    NEWTON  :  computes p(x)/p'(x) by means of Ruffini-Horner's rule   
!    START   :  Selects N starting points by means of Rouche's theorem  
!    CNVEX   :  Computes the convex hull, used by START                 
!    CMERGE  :  Used by CNVEX                                           
!    LEFT    :  Used by CMERGE                                          
!    RIGHT   :  Used by CMERGE                                          
!    CTEST   :  Convexity test, Used by CMERGE                          
!***********************************************************************
!                                                                       
!                                                                       
!***********************************************************************
!********************** SUBROUTINE POLZEROS ****************************
!***********************************************************************
!                        GENERAL COMMENTS                               
!***********************************************************************
!  This routine approximates the roots of   the  polynomial             
!  p(x)=a(n+1)x^n+a(n)x^(n-1)+...+a(1), a(j)=cr(j)+I ci(j), I**2=-1,    
!  where a(1) and a(n+1) are nonzero.                                   
!  The coefficients are COMPLEX*16 numbers. The routine is fast, robust 
!  against overflow, and allows to deal with polynomials of any degree. 
!  Overflow situations are very unlikely and may occurr if there exist  
!  simultaneously coefficients of moduli close to BIG and close to      
!  SMALL, i.e., the greatest and the smallest positive DOUBLE PRECISION 
!  respectively. In this limit situation the program outputs a warning  
!  message. The computation can be speeded up by performing some side   
!  computations in single precision, thus slightly reducing the         
!  robustness of the program (see the comments in the routine ABERTH).  
!  Besides a set of approximations to the roots, the program delivers a 
!  set of a-posteriori error bounds which are guaranteed in the most    
!  part of cases. In the situation where underflow does not allow to    
!  compute a guaranteed bound, the program outputs a warning message    
!  and sets the bound to 0. In the situation where the root cannot be   
!  represented as a COMPLEX*16 number the error bound is set to -1.     
!***********************************************************************
!  The computation is performed by means of Aberth's method             
!  according to the formula                                             
!           x(i)=x(i)-newt/(1-newt*abcorr), i=1,...,n             (1)   
!  where newt=p(x(i))/p'(x(i)) is the Newton correction and abcorr=     
!  =1/(x(i)-x(1))+...+1/(x(i)-x(i-1))+1/(x(i)-x(i+1))+...+1/(x(i)-x(n)) 
!  is the Aberth correction to the Newton method.                       
!***********************************************************************
!  The value of the Newton correction is computed by means of the       
!  synthetic division algorithm (Ruffini-Horner's rule) if |x|<=1,      
!  otherwise the following more robust (with respect to overflow)       
!  formula is applied:                                                  
!                    newt=1/(n*y-y**2 R'(y)/R(y))                 (2)   
!  where                                                                
!                    y=1/x                                              
!                    R(y)=a(1)*y**n+...+a(n)*y+a(n+1)            (2')   
!  This computation is performed by the routine NEWTON.                 
!***********************************************************************
!  The starting approximations are complex numbers that are             
!  equispaced on circles of suitable radii. The radius of each          
!  circle, as well as the number of roots on each circle and the        
!  number of circles, is determined by applying Rouche's theorem        
!  to the functions a(k+1)*x**k and p(x)-a(k+1)*x**k, k=0,...,n.        
!  This computation is performed by the routine START.                  
!***********************************************************************
!                              STOP CONDITION                           
!***********************************************************************
! If the condition                                                      
!                     |p(x(j))|<EPS s(|x(j)|)                      (3)  
! is satisfied,    where      s(x)=s(1)+x*s(2)+...+x**n * s(n+1),       
! s(i)=|a(i)|*(1+3.8*(i-1)),  EPS is the machine precision (EPS=2**-53  
! for the IEEE arithmetic), then the approximation x(j) is not updated  
! and the subsequent iterations (1)  for i=j are skipped.               
! The program stops if the condition (3) is satisfied for j=1,...,n,    
! or if the maximum number NITMAX of  iterations   has   been reached.  
! The condition (3) is motivated by a backward rounding error analysis  
! of the Ruffini-Horner rule, moreover the condition (3) guarantees     
! that the computed approximation x(j) is an exact root of a slightly   
! perturbed polynomial.                                                 
!***********************************************************************
!             INCLUSION DISKS, A-POSTERIORI ERROR BOUNDS                
!***********************************************************************
! For each approximation x of a root, an a-posteriori absolute error    
! bound r is computed according to the formula                          
!                   r=n(|p(x)|+EPS s(|x|))/|p'(x)|                 (4)  
! This provides an inclusion disk of center x and radius r containing a 
! root.                                                                 
!***********************************************************************
!***********************************************************************
!*************       MEANING OF THE INPUT VARIABLES         ************
!***********************************************************************
!***********************************************************************
!                                                                       
!  -- N     : degree of the polynomial.                                 
!  -- POLY  : complex vector of N+1 components, POLY(i) is the          
!           coefficient of x**(i-1), i=1,...,N+1 of the polynomial p(x) 
!  -- EPS   : machine precision of the floating point arithmetic used   
!            by the computer, EPS=2**(-53)  for the IEEE standard.      
!  -- BIG   : the max DOUBLE PRECISION, BIG=2**1023 for the IEEE standar
!  -- SMALL : the min positive DOUBLE PRECISION, SMALL=2**(-1074) for th
!  -- NITMAX: the max number of allowed iterations.                     
!***********************************************************************
!***********************************************************************
!*************      MEANING OF THE OUTPUT VARIABLES         ************
!***********************************************************************
!***********************************************************************
!  ROOT   : complex vector of N components, containing the              
!           approximations to the roots of p(x).                        
!  RADIUS : real vector of N components, containing the error bounds to 
!           the approximations of the roots, i.e. the disk of center    
!           ROOT(i) and radius RADIUS(i) contains a root of p(x), for   
!           i=1,...,N. RADIUS(i) is set to -1 if the corresponding root 
!           cannot be represented as floating point due to overflow or  
!           underflow.                                                  
!  ERR    : vector of N components detecting an error condition;        
!           ERR(j)=.TRUE. if after NITMAX iterations the stop condition 
!                         (3) is not satisfied for x(j)=ROOT(j);        
!           ERR(j)=.FALSE.  otherwise, i.e., the root is reliable,      
!                         i.e., it can be viewed as an exact root of a  
!                         slightly perturbed polynomial.                
!           The vector ERR is used also in the routine convex hull for  
!           storing the abscissae of the vertices of the convex hull.   
!  ITER   : number of iterations peformed.                              
!***********************************************************************
!***********************************************************************
!************    MEANING OF THE AUXILIARY VARIABLES         ************
!***********************************************************************
!***********************************************************************
!  APOLY  : real vector of N+1 components used to store the moduli of   
!           the coefficients of p(x) and the coefficients of s(x) used  
!           to test the stop condition (3).                             
!  APOLYR : real vector of N+1 components used to test the stop         
!           condition                                                   
!***********************************************************************
!*****         WARNING:   2 is the output unit                    ******
!***********************************************************************
      SUBROUTINE POLZEROS (N, POLY, EPS, BIG, SMALL, NITMAX,            &
     & ROOT, RADIUS, ERR, ITER, APOLY,APOLYR)                           
      INTEGER N,I,NZEROS,ITER,NITMAX 
      COMPLEX *16  POLY(N+1),ROOT(N) 
      COMPLEX *16 CORR,ABCORR 
      DOUBLE PRECISION RADIUS(N),APOLY(N+1),APOLYR(N+1) 
      DOUBLE PRECISION AMAX,EPS,SMALL,BIG 
      LOGICAL ERR(N+1) 
! Check consistency of data                                             
      IF (ABS(POLY(N+1)) .EQ. 0.0D0)THEN 
        WRITE(*,*)'Inconsistent data: the leading coefficient is zero' 
        STOP 
      ENDIF 
      IF (ABS(POLY(1)) .EQ. 0.0D0)THEN 
        WRITE(*,*)'The constant term is zero: deflate the polynomial' 
        STOP 
      ENDIF 
! Compute the moduli of the coefficients                                
      AMAX=0 
      DO 10 I=1,N+1 
        APOLY(I)=ABS(POLY(I)) 
        AMAX=MAX(AMAX,APOLY(I)) 
        APOLYR(I)=APOLY(I) 
   10 END DO 
      IF((AMAX).GE.(BIG/(N+1)))THEN 
        WRITE(*,*)'WARNING: COEFFICIENTS TOO BIG, OVERFLOW IS LIKELY' 
        WRITE(2,*)'WARNING: COEFFICIENTS TOO BIG, OVERFLOW IS LIKELY' 
      ENDIF 
! Initialize                                                            
      DO 20 I=1,N 
        RADIUS(I)=0 
        ERR(I)=.TRUE. 
   20 END DO 
! Select the starting points                                            
      CALL PZSTART(N,APOLYR,ROOT,RADIUS,NZEROS,SMALL,BIG,ERR) 
! Compute the coefficients of the backward-error polynomial             
      DO 30 I=1,N+1 
        APOLYR(N-I+2)=EPS*APOLY(I)*(3.8*(N-I+1)+1) 
        APOLY(I)=EPS*APOLY(I)*(3.8*(I-1)+1) 
   30 END DO 
      IF((APOLY(1).EQ.0).OR.(APOLY(N+1).EQ.0))THEN 
        WRITE(*,*)'WARNING: THE COMPUTATION OF SOME INCLUSION RADIUS' 
        WRITE(*,*)'MAY FAIL. THIS IS REPORTED BY RADIUS=0' 
        WRITE(2,*)'WARNING: THE COMPUTATION OF SOME INCLUSION RADIUS' 
        WRITE(2,*)'MAY FAIL. THIS IS REPORTED BY RADIUS=0' 
      ENDIF 
      DO 40 I=1,N 
        ERR(I)=.TRUE. 
        IF(RADIUS(I).EQ.-1)ERR(I)=.FALSE. 
   40 END DO 
! Starts Aberth's iterations                                            
      DO 50 ITER=1,NITMAX 
        DO 50 I=1,N 
          IF (ERR(I)) THEN 
            CALL NEWTON(N,POLY,APOLY,APOLYR,ROOT(I),SMALL,RADIUS(I),    &
     &           CORR,ERR(I))                                           
            IF (ERR(I)) THEN 
              CALL ABERTH(N,I,ROOT,ABCORR) 
              ROOT(I)=ROOT(I)-CORR/(1-CORR*ABCORR) 
            ELSE 
              NZEROS=NZEROS+1 
              IF (NZEROS.EQ.N) RETURN 
            ENDIF 
          ENDIF 
   50 CONTINUE 
      END                                           
                                                                        
                                                                        
                                                                        
                                                                        
!***********************************************************************
!                             SUBROUTINE NEWTON                         
!***********************************************************************
! Compute  the Newton's correction, the inclusion radius (4) and checks 
! the stop condition (3)                                                
!***********************************************************************
! Input variables:                                                      
!     N     : degree of the polynomial p(x)                             
!     POLY  : coefficients of the polynomial p(x)                       
!     APOLY : upper bounds on the backward perturbations on the         
!             coefficients of p(x) when applying Ruffini-Horner's rule  
!     APOLYR: upper bounds on the backward perturbations on the         
!             coefficients of p(x) when applying (2), (2')              
!     Z     : value at which the Newton correction is computed          
!     SMALL : the min positive DOUBLE PRECISION, SMALL=2**(-1074) for th
!***********************************************************************
! Output variables:                                                     
!     RADIUS: upper bound to the distance of Z from the closest root of 
!             the polynomial computed according to (4).                 
!     CORR  : Newton's correction                                       
!     AGAIN : this variable is .true. if the computed value p(z) is     
!             reliable, i.e., (3) is not satisfied in Z. AGAIN is       
!             .false., otherwise.                                       
!***********************************************************************
      SUBROUTINE NEWTON(N,POLY,APOLY,APOLYR,Z,SMALL,RADIUS,CORR,AGAIN) 
      INTEGER N,I 
      COMPLEX *16 POLY(N+1) 
      COMPLEX *16 Z,P,P1,CORR,ZI,DEN,PPSP 
      DOUBLE PRECISION RADIUS,APOLY(N+1),APOLYR(N+1) 
      DOUBLE PRECISION AP,AZ,AZI,ABSP 
      DOUBLE PRECISION SMALL 
      LOGICAL AGAIN 
      AZ=ABS(Z) 
! If |z|<=1 then apply Ruffini-Horner's rule for p(z)/p'(z)             
! and for the computation of the inclusion radius                       
      IF(AZ.LE.1)THEN 
        P=POLY(N+1) 
        AP=APOLY(N+1) 
        P1=P 
        DO 10 I=N,2,-1 
          P=P*Z+POLY(I) 
          P1=P1*Z+P 
          AP=AP*AZ+APOLY(I) 
   10   CONTINUE 
        P=P*Z+POLY(1) 
        AP=AP*AZ+APOLY(1) 
        CORR=P/P1 
        ABSP=ABS(P) 
        AP=AP 
        AGAIN=(ABSP.GT.(SMALL+AP)) 
        IF(.NOT.AGAIN) RADIUS=N*(ABSP+AP)/ABS(P1) 
        RETURN 
      ELSE 
! If |z|>1 then apply Ruffini-Horner's rule to the reversed polynomial  
! and use formula (2) for p(z)/p'(z). Analogously do for the inclusion  
! radius.                                                               
        ZI=1/Z 
        AZI=1/AZ 
        P=POLY(1) 
        P1=P 
        AP=APOLYR(N+1) 
        DO 20 I=N,2,-1 
          P=P*ZI+POLY(N-I+2) 
          P1=P1*ZI+P 
          AP=AP*AZI+APOLYR(I) 
   20   CONTINUE 
        P=P*ZI+POLY(N+1) 
        AP=AP*AZI+APOLYR(1) 
        ABSP=ABS(P) 
        AGAIN=(ABSP.GT.(SMALL+AP)) 
        PPSP=(P*Z)/P1 
        DEN=N*PPSP-1 
        CORR=Z*(PPSP/DEN) 
        IF(AGAIN)RETURN 
        RADIUS=ABS(PPSP)+(AP*AZ)/ABS(P1) 
        RADIUS=N*RADIUS/ABS(DEN) 
        RADIUS=RADIUS*AZ 
      ENDIF 
      END                                           
                                                                        
                                                                        
                                                                        
!***********************************************************************
!                             SUBROUTINE ABERTH                         
!***********************************************************************
! Compute  the Aberth correction. To save time, the reciprocation of    
! ROOT(J)-ROOT(I) could be performed in single precision (complex*8)    
! In principle this might cause overflow if both ROOT(J) and ROOT(I)    
! have too small moduli.                                                
!***********************************************************************
! Input variables:                                                      
!     N     : degree of the polynomial                                  
!     ROOT  : vector containing the current approximations to the roots 
!     J     : index of the component of ROOT with respect to which the  
!             Aberth correction is computed                             
!***********************************************************************
! Output variable:                                                      
!     ABCORR: Aberth's correction (compare (1))                         
!***********************************************************************
      SUBROUTINE ABERTH(N,J,ROOT,ABCORR) 
      INTEGER N,I,J 
      COMPLEX *16 ABCORR,ZJ 
      COMPLEX *16 ROOT(N) 
! The next variable Z could be defined as complex*8 to speed up the     
! computation, this slightly reduces the robustness of the program      
      COMPLEX *16 Z 
      ABCORR=0 
      ZJ=ROOT(J) 
      DO 10 I=1,J-1 
        Z=ZJ-ROOT(I) 
        ABCORR=ABCORR+1/Z 
   10 END DO 
      DO 20 I=J+1,N 
        Z=ZJ-ROOT(I) 
        ABCORR=ABCORR+1/Z 
   20 END DO 
      END                                           
                                                                        
!***********************************************************************
!                             SUBROUTINE START                          
!***********************************************************************
! Compute  the starting approximations of the roots                     
!***********************************************************************
! Input variables:                                                      
!     N     :  number of the coefficients of the polynomial             
!     A     :  moduli of the coefficients of the polynomial             
!     SMALL : the min positive DOUBLE PRECISION, SMALL=2**(-1074) for th
!     BIG   : the max DOUBLE PRECISION, BIG=2**1023 for the IEEE standar
! Output variables:                                                     
!     Y     :  starting approximations                                  
!     RADIUS:  if a component is -1 then the corresponding root has a   
!              too big or too small modulus in order to be represented  
!              as double float with no overflow/underflow               
!     NZ    :  number of roots which cannot be represented without      
!              overflow/underflow                                       
! Auxiliary variables:                                                  
!     H     :  needed for the computation of the convex hull            
!***********************************************************************
! This routines selects starting approximations along circles center at 
! 0 and having suitable radii. The computation of the number of circles 
! and of the corresponding radii is performed by computing the upper    
! convex hull of the set (i,log(A(i))), i=1,...,n+1.                    
!***********************************************************************
      SUBROUTINE PZSTART(N,A,Y,RADIUS,NZ,SMALL,BIG,H) 
      INTEGER N,NZ,I,IOLD,NZEROS,J,JJ 
      LOGICAL H(N+1) 
      COMPLEX*16 Y(N) 
      DOUBLE PRECISION A(N+1),RADIUS(N) 
      DOUBLE PRECISION R,TH,ANG,TEMP,SIGMA,PI2 
      DOUBLE PRECISION SMALL,BIG,XSMALL,XBIG 
      PARAMETER(PI2=6.2831853071796,SIGMA=0.7) 
      XSMALL=LOG(SMALL) 
      XBIG=LOG(BIG) 
      NZ=0 
! Compute the logarithm A(I) of the moduli of the coefficients of       
! the polynomial and then the upper covex hull of the set (A(I),I)      
      DO 10 I=1,N+1 
        IF(A(I).NE.0)THEN 
          A(I)=LOG(A(I)) 
          ELSE 
          A(I)=-1.D30 
        ENDIF 
   10 END DO 
      CALL CNVEX(N+1,A,H) 
! Given the upper convex hull of the set (A(I),I) compute the moduli    
! of the starting approximations by means of Rouche's theorem           
      IOLD=1 
      TH=PI2/N 
      DO 20 I=2,N+1 
        IF (H(I)) THEN 
          NZEROS=I-IOLD 
          TEMP=(A(IOLD)-A(I))/NZEROS 
! Check if the modulus is too small                                     
          IF((TEMP.LT.-XBIG).AND.(TEMP.GE.XSMALL))THEN 
            WRITE(*,*)'WARNING:',NZEROS,' ZERO(S) ARE TOO SMALL TO' 
            WRITE(*,*)'REPRESENT THEIR INVERSES AS COMPLEX*16, THEY' 
            WRITE(*,*)'ARE REPLACED BY SMALL NUMBERS, THE CORRESPONDING' 
            WRITE(*,*)'RADII ARE SET TO -1' 
            WRITE(2,*)'WARNING:',NZEROS,' ZERO(S) ARE TOO SMALL TO ' 
            WRITE(2,*)'REPRESENT THEIR INVERSES AS COMPLEX*16, THEY' 
            WRITE(2,*)'ARE REPLACED BY SMALL NUMBERS, THE CORRESPONDING' 
            WRITE(2,*)'RADII ARE SET TO -1' 
            NZ=NZ+NZEROS 
            R=1.0D0/BIG 
          ENDIF 
          IF(TEMP.LT.XSMALL)THEN 
            NZ=NZ+NZEROS 
            WRITE(*,*)'WARNING: ',NZEROS,' ZERO(S) ARE TOO SMALL TO BE' 
            WRITE(*,*)'REPRESENTED AS COMPLEX*16, THEY ARE SET TO 0' 
            WRITE(*,*)'THE CORRESPONDING RADII ARE SET TO -1' 
            WRITE(2,*)'WARNING: ',NZEROS,' ZERO(S) ARE TOO SMALL TO BE' 
            WRITE(2,*)'REPRESENTED AS COMPLEX*16, THEY ARE SET 0' 
            WRITE(2,*)'THE CORRESPONDING RADII ARE SET TO -1' 
          ENDIF 
! Check if the modulus is too big                                       
          IF(TEMP.GT.XBIG)THEN 
            R=BIG 
            NZ=NZ+NZEROS 
            WRITE(*,*)'WARNING: ',NZEROS,' ZEROS(S) ARE TOO BIG TO BE' 
            WRITE(*,*)'REPRESENTED AS COMPLEX*16,' 
            WRITE(*,*)'THE CORRESPONDING RADII ARE SET TO -1' 
            WRITE(2,*)'WARNING: ',NZEROS,' ZERO(S) ARE TOO BIG TO BE' 
            WRITE(2,*)'REPRESENTED AS COMPLEX*16,' 
            WRITE(2,*)'THE CORRESPONDING RADII ARE SET TO -1' 
          ENDIF 
          IF((TEMP.LE.XBIG).AND.(TEMP.GT.MAX(-XBIG,XSMALL)))THEN 
            R=EXP(TEMP) 
          ENDIF 
! Compute NZEROS approximations equally distributed in the disk of      
! radius R                                                              
          ANG=PI2/NZEROS 
          DO 30 J=IOLD,I-1 
            JJ=J-IOLD+1 
            IF((R.LE.(1.0D0/BIG)).OR.(R.EQ.BIG))RADIUS(J)=-1 
            Y(J)=R*(COS(ANG*JJ+TH*I+SIGMA)+(0,1)*SIN(ANG*JJ+TH*I+SIGMA)) 
   30     CONTINUE 
          IOLD=I 
        ENDIF 
   20 END DO 
      END                                           
                                                                        
                                                                        
!***********************************************************************
!                             SUBROUTINE CNVEX                          
!***********************************************************************
! Compute  the upper convex hull of the set (i,a(i)), i.e., the set of  
! vertices (i_k,a(i_k)), k=1,2,...,m, such that the points (i,a(i)) lie 
! below the straight lines passing through two consecutive vertices.    
! The abscissae of the vertices of the convex hull equal the indices of 
! the TRUE  components of the logical output vector H.                  
! The used method requires O(nlog n) comparisons and is based on a      
! divide-and-conquer technique. Once the upper convex hull of two       
! contiguous sets  (say, {(1,a(1)),(2,a(2)),...,(k,a(k))} and           
! {(k,a(k)), (k+1,a(k+1)),...,(q,a(q))}) have been computed, then       
! the upper convex hull of their union is provided by the subroutine    
! CMERGE. The program starts with sets made up by two consecutive       
! points, which trivially constitute a convex hull, then obtains sets   
! of 3,5,9... points,  up to  arrive at the entire set.                 
! The program uses the subroutine  CMERGE; the subroutine CMERGE uses   
! the subroutines LEFT, RIGHT and CTEST. The latter tests the convexity 
! of the angle formed by the points (i,a(i)), (j,a(j)), (k,a(k)) in the 
! vertex (j,a(j)) up to within a given tolerance TOLER, where i<j<k.    
!***********************************************************************
      SUBROUTINE CNVEX(N,A,H) 
      INTEGER N,I,J,K,M,NJ,JC 
      LOGICAL H(N) 
      DOUBLE PRECISION A(N) 
      DO 10 I=1,N 
        H(I)=.TRUE. 
   10 END DO 
! compute K such that N-2<=2**K<N-1                                     
      K=INT(LOG(N-2.0D0)/LOG(2.0D0)) 
      IF(2**(K+1).LE.(N-2)) K=K+1 
! For each M=1,2,4,8,...,2**K, consider the NJ pairs of consecutive     
! sets made up by M+1 points having the common vertex                   
! (JC,A(JC)), where JC=M*(2*J+1)+1 and J=0,...,NJ,                      
! NJ=MAX(0,INT((N-2-M)/(M+M))).                                         
! Compute the upper convex hull of their union by means of the          
! subroutine CMERGE                                                     
      M=1 
      DO 20 I=0,K 
        NJ=MAX(0,INT((N-2-M)/(M+M))) 
        DO 30 J=0,NJ 
          JC=(J+J+1)*M+1 
          CALL CMERGE(N,A,JC,M,H) 
   30   CONTINUE 
        M=M+M 
   20 END DO 
      END                                           
                                                                        
!***********************************************************************
!                             SUBROUTINE LEFT                           
!***********************************************************************
! Given as input the integer I and the vector H of logical, compute the 
! the maximum integer IL such that IL<I and H(IL) is TRUE.              
!***********************************************************************
! Input variables:                                                      
!     N   : length of the vector H                                      
!     H   : vector of logical                                           
!     I   : integer                                                     
!***********************************************************************
! Output variable:                                                      
!     IL  : maximum integer such that IL<I, H(IL)=.TRUE.                
!***********************************************************************
      SUBROUTINE LEFT(N,H,I,IL) 
      INTEGER N,I,IL 
      LOGICAL H(N) 
      DO 10 IL=I-1,0,-1 
        IF (H(IL)) RETURN 
   10 END DO 
      END                                           
                                                                        
!***********************************************************************
!                             SUBROUTINE RIGHT                          
!***********************************************************************
!***********************************************************************
! Given as input the integer I and the vector H of logical, compute the 
! the minimum integer IR such that IR>I and H(IL) is TRUE.              
!***********************************************************************
!***********************************************************************
! Input variables:                                                      
!     N   : length of the vector H                                      
!     H   : vector of logical                                           
!     I   : integer                                                     
!***********************************************************************
! Output variable:                                                      
!     IR  : minimum integer such that IR>I, H(IR)=.TRUE.                
!***********************************************************************
      SUBROUTINE RIGHT(N,H,I,IR) 
      INTEGER N,I,IR 
      LOGICAL H(N) 
      DO 10 IR=I+1,N 
        IF (H(IR)) RETURN 
   10 END DO 
      END                                           
                                                                        
!***********************************************************************
!                             SUBROUTINE CMERGE                         
!***********************************************************************
! Given the upper convex hulls of two consecutive sets of pairs         
! (j,A(j)), compute the upper convex hull of their union                
!***********************************************************************
! Input variables:                                                      
!     N    : length of the vector A                                     
!     A    : vector defining the points (j,A(j))                        
!     I    : abscissa of the common vertex of the two sets              
!     M    : the number of elements of each set is M+1                  
!***********************************************************************
! Input/Output variable:                                                
!     H    : vector defining the vertices of the convex hull, i.e.,     
!            H(j) is .TRUE. if (j,A(j)) is a vertex of the convex hull  
!            This vector is used also as output.                        
!***********************************************************************
      SUBROUTINE CMERGE(N,A,I,M,H) 
      INTEGER N,M,I,IR,IL,IRR,ILL 
      LOGICAL H(N) 
      LOGICAL TSTL,TSTR,CTEST 
      DOUBLE PRECISION A(N) 
! at the left and the right of the common vertex (I,A(I)) determine     
! the abscissae IL,IR, of the closest vertices of the upper convex      
! hull of the left and right sets, respectively                         
      CALL LEFT(N,H,I,IL) 
      CALL RIGHT(N,H,I,IR) 
! check the convexity of the angle formed by IL,I,IR                    
      IF (CTEST(N,A,IL,I,IR)) THEN 
        RETURN 
      ELSE 
! continue the search of a pair of vertices in the left and right       
! sets which yield the upper convex hull                                
        H(I)=.FALSE. 
   10   CONTINUE 
        IF (IL.EQ.(I-M)) THEN 
          TSTL=.TRUE. 
          ELSE 
          CALL LEFT(N,H,IL,ILL) 
          TSTL=CTEST(N,A,ILL,IL,IR) 
        ENDIF 
        IF (IR.EQ.MIN(N,I+M))THEN 
          TSTR=.TRUE. 
        ELSE 
          CALL RIGHT(N,H,IR,IRR) 
          TSTR=CTEST(N,A,IL,IR,IRR) 
        ENDIF 
        H(IL)=TSTL 
        H(IR)=TSTR 
        IF (TSTL.AND.TSTR) RETURN 
        IF(.NOT.TSTL) IL=ILL 
        IF(.NOT.TSTR) IR=IRR 
        GOTO 10 
      ENDIF 
      END                                           
                                                                        
!***********************************************************************
!                             FUNCTION CTEST                            
!***********************************************************************
! Test the convexity of the angle formed by (IL,A(IL)), (I,A(I)),       
! (IR,A(IR)) at the vertex (I,A(I)), up to within the tolerance         
! TOLER. If convexity holds then the function is set to .TRUE.,         
! otherwise CTEST=.FALSE. The parameter TOLER is set to 0.4 by default. 
!***********************************************************************
! Input variables:                                                      
!     N       : length of the vector A                                  
!     A       : vector of double                                        
!     IL,I,IR : integers such that IL<I<IR                              
!***********************************************************************
! Output:                                                               
!     .TRUE. if the angle formed by (IL,A(IL)), (I,A(I)), (IR,A(IR)) at 
!            the vertex (I,A(I)), is convex up to within the tolerance  
!            TOLER, i.e., if                                            
!            (A(I)-A(IL))*(IR-I)-(A(IR)-A(I))*(I-IL)>TOLER.             
!     .FALSE.,  otherwise.                                              
!***********************************************************************
      FUNCTION CTEST(N,A,IL,I,IR) 
      INTEGER N,I,IL,IR 
      DOUBLE PRECISION A(N) 
      DOUBLE PRECISION S1,S2,TOLER 
      LOGICAL CTEST 
      PARAMETER(TOLER=0.4) 
      S1=A(I)-A(IL) 
      S2=A(IR)-A(I) 
      S1=S1*(IR-I) 
      S2=S2*(I-IL) 
      CTEST=.FALSE. 
      IF(S1.GT.(S2+TOLER)) CTEST=.TRUE. 
      END                                           

!                                                                       
!  ***************************************************************      
!  *                                                             *      
!  *                        R V N O R M                          *      
!  *                                                             *      
!  *            Random variable normally distributed             *      
!  *                    (mean=0, variance=1)                     *      
!  *                                                             *      
!  ***************************************************************      
!                                                                       
      DOUBLE PRECISION FUNCTION rvnorm() 
      IMPLICIT NONE 
                                                                        
      DOUBLE PRECISION v1,v2,s,x 
                                                                        
      LOGICAL stor 
      DATA stor/.false./ 
                                                                        
      SAVE stor,v2,x 
                                                                        
      DOUBLE PRECISION urand 
      EXTERNAL urand 
                                                                        
      IF(stor) THEN 
          rvnorm=v2*x 
          stor=.false. 
      ELSE 
    1     v1=2.d0*urand()-1.d0 
          v2=2.d0*urand()-1.d0 
          s=v1**2+v2**2 
          IF(s.GE.1.d0) GOTO 1 
          x=SQRT(-(2.d0*LOG(s)/s)) 
          rvnorm=v1*x 
          stor=.true. 
      END IF 
                                                                        
      END                                           
!                                                                       
!  ***************************************************************      
!  *                                                             *      
!  *                           R A N D                           *      
!  *                                                             *      
!  *            Random variable uniformly distributed            *      
!  *                       between 0 and 1                       *      
!  *                                                             *      
!  ***************************************************************      
!                                                                       
      DOUBLE PRECISION FUNCTION urand() 
      IMPLICIT NONE 
                                                                        
      INTEGER savrnd 
                                                                        
      INTEGER ia(55),jrand,i,unit 
      DATA jrand/0/ 
                                                                        
      SAVE jrand,ia 
                                                                        
      INTEGER irn55 
      EXTERNAL irn55 
                                                                        
      IF(jrand.EQ.0) THEN 
          CALL filopl(unit,'random.dat') 
          DO 1 i=1,55 
          READ(unit,100) ia(i) 
    1     CONTINUE 
          CALL filclo(unit,' ') 
      END IF 
                                                                        
      jrand=jrand+1 
      IF(jrand.GT.55) jrand=irn55(ia) 
      urand=ia(jrand)*1.D-9 
                                                                        
      RETURN 
!                                                                       
!  ***************************************************************      
!  *                                                             *      
!  *                         S A V R N D                         *      
!  *                                                             *      
!  *             Generates a new version of the file             *      
!  *      containing the seeds for random number generation      *      
!  *                                                             *      
!  ***************************************************************      
!                                                                       
      ENTRY savrnd() 
                                                                        
      savrnd=irn55(ia) 
      CALL filopn(unit,'random.dat','UNKNOWN') 
      DO 2 i=1,55 
      WRITE(unit,100) ia(i) 
  100 FORMAT(I10) 
    2 END DO 
      CALL filclo(unit,' ') 
                                                                        
      END                                           
!                                                                       
!  ***************************************************************      
!  *                                                             *      
!  *                         I N T R N D                         *      
!  *                                                             *      
!  *                 Initialization of the file                  *      
!  *      containing the seeds for random number generation      *      
!  *                                                             *      
!  ***************************************************************      
!                                                                       
      SUBROUTINE intrnd 
      IMPLICIT NONE 
                                                                        
      INTEGER ia(55),ix,unit,i 
                                                                        
      WRITE(*,*)' ENTER SEED' 
      READ(*,*) ix 
      CALL in55(ia,ix) 
      CALL filopn(unit,'random.dat','UNKNOWN') 
      DO 1 i=1,55 
      WRITE(unit,100) ia(i) 
  100 FORMAT(I10) 
    1 END DO 
      CALL filclo(unit,' ') 
                                                                        
      END                                           
!                                                                       
!  ***************************************************************      
!  *                                                             *      
!  *                          I R N 5 5                          *      
!  *                                                             *      
!  ***************************************************************      
!                                                                       
      INTEGER FUNCTION irn55(ia) 
      IMPLICIT NONE 
                                                                        
      INTEGER ia(*) 
                                                                        
      INTEGER i,j 
                                                                        
      DO 1 i=1,24 
      j=ia(i)-ia(i+31) 
      IF(j.LT.0) j=j+1000000000 
      ia(i)=j 
    1 END DO 
                                                                        
      DO 2 i=25,55 
      j=ia(i)-ia(i-24) 
      IF(j.LT.0) j=j+1000000000 
      ia(i)=j 
    2 END DO 
                                                                        
      irn55=1 
                                                                        
      END                                           
!                                                                       
!  ***************************************************************      
!  *                                                             *      
!  *                           I N 5 5                           *      
!  *                                                             *      
!  ***************************************************************      
!                                                                       
      SUBROUTINE in55(ia,ix) 
      IMPLICIT NONE 
                                                                        
      INTEGER ia(*),ix 
                                                                        
      INTEGER i,j,k,ii 
                                                                        
      INTEGER irn55 
      EXTERNAL irn55 
                                                                        
      ia(55)=ix 
      j=ix 
      k=1 
      DO 1 i=1,54 
      ii=MOD(21*i,55) 
      ia(ii)=k 
      k=j-k 
      IF(k.LT.0) k=k+1000000000 
      j=ia(ii) 
    1 END DO 
                                                                        
      DO 2 i=1,10 
      ii=irn55(ia) 
    2 END DO 
                                                                        
      END                                           

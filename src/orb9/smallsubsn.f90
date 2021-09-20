! =====================================                                 
!  small routines                                                       
! =====================================
! RMS with respect to mean                                              
DOUBLE PRECISION FUNCTION sigma(x,n) 
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT(IN) :: x(n)
  INTEGER, INTENT(IN) :: n
! END INTERFACE
  INTEGER i
  DOUBLE PRECISION av 
  av=sum(x(1:n))/n 
  sigma=0.d0 
  DO i=1,n 
     sigma=sigma+(x(i)-av)**2
  ENDDO
  sigma=SQRT(sigma/n)  
END FUNCTION sigma
! =====================================                                 
! RMS with respect to fixed value                                       
DOUBLE PRECISION FUNCTION sigmf(x,av,n) 
  IMPLICIT NONE 
  INTEGER, INTENT(IN) :: n 
  DOUBLE PRECISION, INTENT(IN) ::  x(n),av
  INTEGER i
  sigmf=0.d0 
  DO i=1,n 
     sigmf=sigmf+(x(i)-av)**2 
  ENDDO
  sigmf=sqrt(sigmf/n) 
END FUNCTION sigmf

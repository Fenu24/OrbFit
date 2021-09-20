!                                                                       
! Computation of the right hand side of the equations of motion         
!                                                                       
      SUBROUTINE forceb(x,v,tm,f) 
      IMPLICIT NONE 
                                                                        
      INCLUDE 'dim.h90' 
                                                                        
      DOUBLE PRECISION x(3,nbtx),v(3,nbtx),f(3,nbtx),tm,et2(2),r6(6) 
                                                                        
      DOUBLE PRECISION r3(nbtx),xjpl(3,nbjx),r2,fi(3),r3j(nbjx),rr3 
      INTEGER j,k 
                                                                        
      INCLUDE 'comast.h90' 
      INCLUDE 'comjpl.h90' 
                                                                        
      et2(1)=2400000.5d0 
      et2(2)=tm 
                                                                        
! Distance from the Sun (asteroids)                                     
      DO 1 k=1,nbt 
      r2=x(1,k)**2+x(2,k)**2+x(3,k)**2 
      r3(k)=SQRT(r2)*r2 
    1 END DO 
! Distance from the Sun (JPL planets)                                   
      DO 10 k=1,nbj 
! Call DPLEPH, computing only positions (ISTATE=1)                      
      CALL dpleph(et2,jplid(k),11,r6,1) 
      xjpl(1,k)=r6(1) 
      xjpl(2,k)=r6(2) 
      xjpl(3,k)=r6(3) 
      r2=xjpl(1,k)**2+xjpl(2,k)**2+xjpl(3,k)**2 
      r3j(k)=SQRT(r2)*r2 
   10 END DO 
                                                                        
      DO 2 k=1,nbt 
      fi(1)=0 
      fi(2)=0 
      fi(3)=0 
      DO 3 j=1,nbm 
      IF(j.NE.k) THEN 
          r2=(x(1,k)-x(1,j))**2+(x(2,k)-x(2,j))**2+(x(3,k)-x(3,j))**2 
          rr3=SQRT(r2)*r2 
          fi(1)=fi(1)+gma(j)*((x(1,j)-x(1,k))/rr3-x(1,j)/r3(j)) 
          fi(2)=fi(2)+gma(j)*((x(2,j)-x(2,k))/rr3-x(2,j)/r3(j)) 
          fi(3)=fi(3)+gma(j)*((x(3,j)-x(3,k))/rr3-x(3,j)/r3(j)) 
      END IF 
    3 END DO 
      DO 4 j=1,nbj 
      r2=(x(1,k)-xjpl(1,j))**2                                          &
     &  +(x(2,k)-xjpl(2,j))**2                                          &
     &  +(x(3,k)-xjpl(3,j))**2                                          
      rr3=SQRT(r2)*r2 
      fi(1)=fi(1)+gmj(j)*((xjpl(1,j)-x(1,k))/rr3-xjpl(1,j)/r3j(j)) 
      fi(2)=fi(2)+gmj(j)*((xjpl(2,j)-x(2,k))/rr3-xjpl(2,j)/r3j(j)) 
      fi(3)=fi(3)+gmj(j)*((xjpl(3,j)-x(3,k))/rr3-xjpl(3,j)/r3j(j)) 
    4 END DO 
      IF(k.LE.nbm) THEN 
          f(1,k)=fi(1)-gma1(k)*x(1,k)/r3(k) 
          f(2,k)=fi(2)-gma1(k)*x(2,k)/r3(k) 
          f(3,k)=fi(3)-gma1(k)*x(3,k)/r3(k) 
      ELSE 
          f(1,k)=fi(1)-gmsun*x(1,k)/r3(k) 
          f(2,k)=fi(2)-gmsun*x(2,k)/r3(k) 
          f(3,k)=fi(3)-gmsun*x(3,k)/r3(k) 
      END IF 
    2 END DO 
                                                                        
      END SUBROUTINE forceb                                          

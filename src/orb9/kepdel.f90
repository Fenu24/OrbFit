! {\bf kepdel} subroutine computes transformation of keplerian          
!              a,e,i to Delaunay's L,G,H                                
!                                                                       
SUBROUTINE kepdel(el,del)
  USE massmod, ONLY: gm 
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT(IN) ::   el(6)
  DOUBLE PRECISION, INTENT(OUT) ::  del(6) 
  del(1)=sqrt(gm(1))*sqrt(el(1)) 
  del(2)=del(1)*sqrt(1.d0-el(2)*el(2)) 
  del(3)=del(2)*cos(el(3)) 
  del(4:6)=el(4:6) 
END SUBROUTINE kepdel
! ===========================================                           
! {\bf delkep} subroutine computes transformation of Delaunay's         
!              L,G,H to keplerian a,e,i                                 
!                                                                       
SUBROUTINE delkep(del,el) 
  USE massmod, ONLY: gm
  USE fund_const, ONLY: dpig
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT(IN) ::   del(6)
  DOUBLE PRECISION, INTENT(OUT) ::  el(6)
  DOUBLE PRECISION ecsq, cosi
  el(1)=(del(1))**2/gm(1) 
  ecsq=1.d0-(del(2)/del(1))**2 
  if (ecsq .lt. 0.d0)then 
     el(2)=0.d0 
  else 
     el(2)=sqrt(ecsq) 
  endif
  cosi=del(3)/del(2) 
  if (cosi .lt. 0.d0) cosi=-cosi 
  if(cosi.gt.1.d0)cosi=1.d0 
  el(3)=dpig/4.d0-asin(cosi) 
  el(4:6)=del(4:6) 
END SUBROUTINE delkep

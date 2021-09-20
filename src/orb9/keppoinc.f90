! ==================================================                    
!  keppoinc                                                             
SUBROUTINE keppoinc(em,ap0,aq0,ar0,as0,al,a,aa,aunit,sm) 
  USE fund_const
  IMPLICIT NONE 
! input                                                                 
  DOUBLE PRECISION em(6),aunit,sm 
! output                                                                
  DOUBLE PRECISION ap0,aq0,ar0,as0,al,a,aa 
  DOUBLE PRECISION aw,ao,ai,eh,ek,tgi2,ep,eq,eps,ae,ati2,ag,ah,awt 
!   change from degrees to radians                                      
  aw=em(5)*radeg 
  ao=em(4)*radeg 
  ai=em(3)*radeg 
! semimajor axis                                                        
  aa=em(1) 
!   equinoctal elements                                                 
  eh=em(2)*dsin(aw+ao) 
  ek=em(2)*dcos(aw+ao) 
  tgi2=dtan(ai*0.5d0) 
  ep=tgi2*dsin(ao) 
  eq=tgi2*dcos(ao) 
!     delaunay's and poincare's variables                               
!     change from equinoctal to delaunay                                
!     equinoctal elements are: a, h=e sin w, k=e cos w,                 
!     p=tg(i/2) sin omega, q=tg(i/2) cos omega, mean longitude          
!     eps=control for essentially zero mean                             
!     eccentricity and/or inclination                                   
  eps=1.0d-8 
  a=aa/aunit 
  al=sm*dsqrt(a) 
  ae=dsqrt(eh**2+ek**2) 
  ati2=dsqrt(ep**2+eq**2) 
  ai=datan(ati2)*2.d0 
  ag=al*dsqrt(1.d0-ae*ae) 
  ah=ag*dcos(ai) 
  if(ai.gt.eps)then 
     ao=datan2(ep,eq) 
  else 
     ao=0.d0 
  endif
  if(ae.gt.eps)then 
     awt=datan2(eh,ek) 
  else 
     awt=0.d0 
  endif
  aw=awt-ao 
  ap0=dsqrt(2.d0*(al-ag))*dcos(aw+ao) 
  aq0=-dsqrt(2.d0*(al-ag))*dsin(aw+ao) 
  ar0=dsqrt(2.d0*(ag-ah))*dcos(ao) 
  as0=-dsqrt(2.d0*(ag-ah))*dsin(ao) 
 
END SUBROUTINE keppoinc

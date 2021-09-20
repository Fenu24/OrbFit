! =========moid_compute==============
! ROUTINES                                                              
!             nomoid                                                    
!             nodedi                                                    
!                                                                       
!  HEADERS                                                             
! moid_compute.o:                                                       
!            parbep.h masses.h masses Earth and Moon                    
!            fund_const.mod for mass of the Sun

! =========================================                             
! NOMOID - compute moid using Gronchi's routine                         
! nodal distances and Minimum Orbital Intersection Distance             
! with respect to the Earth                                             
! for an orbit, given in any elements, at a given time           
! =========================================                             
SUBROUTINE nomoid(t0,el0,moid,dnp,dnm) 
  USE fund_const
  USE planet_masses
  USE orbit_elements
  USE critical_points
  IMPLICIT NONE 
! ===========INPUT=====================                                 
! epoch (MJD)                       
  DOUBLE PRECISION, INTENT(IN) :: t0 
  TYPE(orbit_elem), INTENT(IN) :: el0
! ===========OUTPUT=================================                    
! MOID, iteration count and succ. flag                                  
! ascending node distance,  descending node distance                    
  DOUBLE PRECISION moid,dnp,dnm 
! ==========END INTERFACE================================               
! elements of Earth (equinoctal*, cartesian coordinates   
  TYPE(orbit_elem) :: elea,elcar
  INTEGER fail_flag,fail_flag1 
! minimal values of distance, signed                  
  DOUBLE PRECISION dmintil(nminx)
! number of relative minima found                                   
  INTEGER nummin 
! ======================================================                
! Earth's elements at the same time                                    
  elea=undefined_orbit_elem
  elea%t=t0
  elea%coo='CAR'
  CALL earcar(t0,elea%coord,1)
! compute moid  
!  WRITE(*,*) elea,el0
  CALL dmintil_rms(elea,el0,nummin,dmintil)
  moid=abs(dmintil(1)) 
! compute nodal distances
  CALL coo_cha(el0,'CAR',elcar, fail_flag1)
  CALL nodedi(elcar%coord,elea%coord,dnp,dnm) 
END SUBROUTINE nomoid
! ==========================================                            
! NODEDI                                                                
! nodal distances of two orbits                                
! ==========================================                            
SUBROUTINE nodedi(x,xpl,dnp,dnm) 
  USE fund_const
  USE planet_masses
  IMPLICIT NONE 
! ===========INPUT=====================                                 
! cartesian coordinates of asteroid, of Earth 
  DOUBLE PRECISION x(6),xpl(6) 
! ===========OUTPUT=================================                    
! output ascending node distance,  descending node distance             
  DOUBLE PRECISION dnp,dnm 
! ==========END INTERFACE================================               
  DOUBLE PRECISION c(3),cpl(3),vlenz(3),vlenzpl(3) 
  DOUBLE PRECISION vnod(3),vnl,ome,omepl,qpl,q, c2,c2pl 
  DOUBLE PRECISION ecc,eccpl,cosf,cosfpl,rp,rm,rplp,rplm 
  INTEGER i 
  DOUBLE PRECISION prscal,vsize 
!  angular momentum                                                     
  CALL prvec(x,x(4),c) 
  CALL prvec(xpl,xpl(4),cpl) 
! ascending node                                                        
  CALL prvec(cpl,c,vnod) 
  vnl=vsize(vnod) 
!  angular momentum unit vector, Lenz vector                            
  CALL  prvec(x(4),c,vlenz) 
  DO i=1,3 
     vlenz(i)=vlenz(i)/gms-x(i)/vsize(x) 
  ENDDO
  ecc=vsize(vlenz) 
  CALL  prvec(xpl(4),cpl,vlenzpl) 
  DO i=1,3 
     vlenzpl(i)=vlenzpl(i)/gmse-xpl(i)/vsize(xpl) 
  ENDDO
  eccpl=vsize(vlenzpl) 
! true anomaly at mutual node= - arg. of perihelion                     
  cosf=prscal(vnod,vlenz)/(vnl*ecc) 
  ome=acos(cosf) 
  cosfpl=prscal(vnod,vlenzpl)/(vnl*eccpl) 
  omepl=acos(cosfpl) 
! nodal points and distances
  c2=prscal(c,c)
  q=(c2/gms)/(1.d0+ecc)
  c2pl=prscal(cpl,cpl)
  qpl=(c2pl/gmse)/(1.d0+eccpl)
  rp=q*(1.d0+ecc)/(1.d0+ecc*cosf)
  rm=q*(1.d0+ecc)/(1.d0-ecc*cosf)
  rplp=qpl*(1.d0+eccpl)/(1.d0+eccpl*cosfpl)
  rplm=qpl*(1.d0+eccpl)/(1.d0-eccpl*cosfpl)
  dnp=rp-rplp
  dnm=rm-rplm
  IF(rp.lt.0.d0.or.rplp.lt.0.d0)THEN
     dnp=dnm ! one of the two must be real
  ELSEIF(rm.lt.0.d0.or.rplm.lt.0.d0)THEN
     dnm=dnp
  ENDIF
!  WRITE(*,*)ecc,cosf,rp,rm,dnm,dnp     
END SUBROUTINE nodedi

! =============================================                         
! VAREQI                                                                
! variational vectors input                                             
SUBROUTINE vareqi(ilca,semim,jjin,xv,v0,v1,vvel,vpos,vslog,yv) 
  USE massmod  
  IMPLICIT NONE 
! =============INPUT=================                                   
  INTEGER ilca 
  DOUBLE PRECISION semim 
! arrays for variational equations                                      
  DOUBLE PRECISION xv(6,nvzx),jjin(nastx) 
! =============OUTPUT================                                   
! normalisation controls                                                
  DOUBLE PRECISION v0,v1,vvel,vpos,vslog(nvzx) 
! arrays for variational equations                                      
  DOUBLE PRECISION yv(6,nvzx) 
  INTEGER jv 
! ============END INTERFACE===========                                  
! headers                                                               
  INCLUDE 'comnbo.h90' 
! random number seed                                                    
!  INTEGER ix
! random number function
  DOUBLE PRECISION urand 
! arrays for variational equations                                      
  DOUBLE PRECISION dy(6) 
! scalar temporaries                                                    
  DOUBLE PRECISION sv,vv,vnor,renorm,sd,sn,snp,snv 
! functions                                                             
  DOUBLE PRECISION drandu 
! loop indexes                                                          
  INTEGER j,i 
! ********************************************************************  
!  variational equations for the first nvz asteroids: initial           
!  conditions are chosen at random, unless there is continuation        
! **********************************************************            
!  ix=1234567 
  v0=1.d0 
!  coefficients of the metric to be used for LCE                        
!  at semimajor axis semim, position and velocity vector                
!  have both norm 1 for a circular orbit                                
  vpos=1.d0/semim**2 
  vvel=semim/gm(1) 
!  initial conditions for variational equations                         
  do 40 j=1,nvz 
     jv=jjin(j) 
     if(jv.gt.ilca)then 
!  new variational equation:                                            
!  create a random initial vector for variational equation              
        sd=0.d0 
        DO i=1,6 
           dy(i)=urand() 
           sd=sd+dy(i)**2 
        ENDDO
        sn=v0/dsqrt(sd) 
        snp=sn*semim 
        snv=sn*sqrt(gm(1)/semim) 
        DO i=1,3 
           yv(i,j)=dy(i)*snp 
           yv(3+i,j)=dy(i+3)*snv 
        ENDDO
     elseif(jv.le.ilca)then 
!  continuation of a previous computation:                              
!  renormalize available variation vector (its norm is already          
!  accounted for in gamma)
        yv(1:6,j)=xv(1:6,jv) 
        sv=0.d0 
        vv=0.d0 
        DO i=1,3 
           sv=sv+yv(i,j)**2 
           vv=vv+yv(i+3,j)**2 
        ENDDO
        vnor=sqrt(sv*vpos+vv*vvel) 
        renorm=v0/vnor           
        yv(1:6,j)=yv(1:6,j)*renorm 
     endif
40 ENDDO
END SUBROUTINE vareqi

! =======================================================               
! subs8: various subroutines used by long8;                             
! copyright Milani and Knezevic, 1992, 1997                             
! ********************************************************************  
!   convun                                                              
!   conversion of units for longit/secres                               
! ********************************************************************* 
SUBROUTINE convun(pm,pa,pg,ps,pn,aunit,tunit,pmunit, &
     &   nplin,nplou,fnl,nnl)                                           
  USE fund_const, ONLY: gk
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: nplin, nplou, nnl
! units
  DOUBLE PRECISION, INTENT(OUT) :: aunit, tunit, pmunit
!  secular theory frequencies                                           
  DOUBLE PRECISION, INTENT(INOUT) ::  pg(8),ps(8),fnl(nnl) 
!  masses, semimajor axes and mean motions of the planets               
  DOUBLE PRECISION, INTENT(INOUT) ::  pm(8),pn(8),pa(8) 
! END INTERFACE
  DOUBLE PRECISION pnj, pmj, pmsun
  INTEGER n,j 
!  Mean motion of Jupiter
  pnj=gk*3.6525d2*sqrt((1.d0+pm(5))/pa(5)**3) 
!   unit of time=period of jupiter/dpig                                 
  tunit=1.d0/pnj
!     write(*,*)' tunit=',tunit                                         
!  conversion of frequencies to new time unit                           
  DO j=nplin,nplou 
     pg(j)=pg(j)*tunit 
     ps(j)=ps(j)*tunit 
  ENDDO
  DO n=1,nnl 
     fnl(n)=fnl(n)*tunit 
  ENDDO
!   unit mass=mass of sun+jupiter                                       
  pmj=pm(5) 
  pmunit=1.d0+pmj 
  DO j=nplin,nplou 
     pm(j)=pm(j)/pmunit
  ENDDO
  pmsun=1.d0/pmunit 
!   unit of length=semimajor axis of jupiter                            
  aunit=pa(5) 
  DO j=nplin,nplou 
     pa(j)=pa(j)/aunit
  ENDDO
!  mean motions of planets - n'(jup)=1.d0                               
  DO j=nplin,nplou 
     pn(j)=sqrt((pmsun+pm(j))/(pa(j))**3.d0)
  ENDDO
  if(abs(pn(5)-1.d0).gt.1.d-15)then 
     write(*,*)'pn(5)',pn(5) 
     pn(5)=1.d0
  endif 
END SUBROUTINE convun
!                                                                       
! **********************************************************************
!  inplda                                                               
!  inner planets data (used to compute g0, s0 in longit/secres)         
! **********************************************************************
SUBROUTINE inplda(pa,pmv) 
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT(INOUT) :: pa(8),pmv(8) 
!  semimajor axis from Bretagnon (1974)                                 
  pa(1)=0.387099d0 
  pa(2)=0.723332d0 
  pa(3)=1.d0 
  pa(4)=1.523691d0 
!  masses from Standish and Campbell (1984), see Nobili et al.(1989)    
  pmv(1)=1.d0/6.0236d6 
  pmv(2)=1.d0/4.085235d5 
  pmv(3)=1.d0/3.289055d5 
  pmv(4)=1.d0/3.09871d6 
END SUBROUTINE inplda
! ****************************************************************      
!  gsfour                                                               
! subroutine which calculates the degree four corrections of            
! asteroid free oscillation frequencies g0,s0. output delta g,s.        
! ****************************************************************      
SUBROUTINE gsfour(pm,ya,yb,xi,et,pni,pmi,ani,ami,cg,cs,nplin,nplou)
  IMPLICIT NONE 
! INPUT
  DOUBLE PRECISION, INTENT(IN) :: ya(4),yb(27),xi(8),et(8),pni(8),pmi(8), ani, ami, pm
  INTEGER, INTENT(IN) :: nplin, nplou
! OUTPUT
  DOUBLE PRECISION, INTENT(OUT) :: cg, cs
! END INTERFACE
  INTEGER k
  DOUBLE PRECISION sumxi, sumet,xini,etmi,sumni,summi
  sumxi=0.d0 
  sumet=0.d0 
  xini=0.d0 
  etmi=0.d0 
  sumni=0.d0 
  summi=0.d0 
  DO k=nplin,nplou 
     sumxi=sumxi+xi(k)**2 
     sumet=sumet+et(k)**2 
     xini=xini+xi(k)*pni(k) 
     etmi=etmi+et(k)*pmi(k) 
     sumni=sumni+pni(k)**2 
     summi=summi+pmi(k)**2 
  ENDDO
  cg=pm*((yb(1)-ya(1)/4.d0)*(2.d0*ani+4.d0*sumxi)                   &
     &+(yb(5)+ya(2)/2.d0)*(ami+sumet)+                                  &
     &(yb(10)+ya(3)/4.d0)*etmi+(yb(12)-ya(4)/8.d0)*2.d0*xini+           &
     &yb(2)*sumni+yb(6)*summi)                                          
  cs=pm*((yb(5)+ya(2)/2.d0)*(ani+sumxi)                             &
     &+(yb(3)-ya(2)/4.d0)*(2.d0*ami+4.d0*sumet)+                        &
     &(yb(8)-ya(3)/8.d0)*2.d0*etmi+yb(4)*summi                          &
     &+yb(7)*sumni+yb(14)*xini)                                         
END SUBROUTINE gsfour
! ********************************************************************* 
!   g0s0                                                                
!   second order/degree frequencies for perihelion and node             
!   of the asteroid; planets from nplin to nplou are                    
!   included in the computation                                         
!   second order effects are included in yd; however, if                
!   for a given planet the second order effect is not required,         
!   yd is set to zero in lly                                            
! ********************************************************************* 
SUBROUTINE g0s0(al,ya,pm,g0,s0,nplin,nplou) 
  IMPLICIT NONE
! INPUT
  INTEGER, INTENT(IN) :: nplin, nplou !list of planets
!  masses of the planets                                                
  DOUBLE PRECISION, INTENT(IN) :: pm(8) 
!  yuasa coefficients--warning!! firt index is yuasa index,             
!    second index is planet index; Delaunay L
  DOUBLE PRECISION, INTENT(IN) :: ya(4,8), al 
! OUTPUT
  DOUBLE PRECISION, INTENT(OUT) :: g0, s0
! END INTERFACE
  DOUBLE PRECISION dsl
  INTEGER j
  dsl=2.d0/al 
  g0=0.d0 
  s0=0.d0 
  DO j=nplin,nplou
     g0=g0+dsl*pm(j)*ya(1,j)
     s0=s0+dsl*pm(j)*ya(2,j)
  ENDDO
END SUBROUTINE g0s0
! ********************************************************************  
!  forced                                                               
!  forced oscillation amplitudes.                                       
!  second order effects are accounted for; if they are not              
!  required, yd is set to zero in lly                                   
!  the forced terms included are the ones with frequencies              
!  indexed from nplin to nplou for planets nplin to nplou               
! ********************************************************************  
SUBROUTINE forcedgs(g0,s0,pg,ps,pm,pni,pmi,ya,al,xi,et,nplin,nplou)
  IMPLICIT NONE
! proper frequencies, linear, Delaunay L
   DOUBLE PRECISION, INTENT(IN) :: g0, s0, al
!  linear secular perturbation theory-- in two dim arrays,              
!     first index is planet index, second index is frequency index      
!  secular theory constants                                             
  DOUBLE PRECISION, INTENT(IN) :: pni(8,8),pg(8),pmi(8,8),ps(8) 
  INTEGER, INTENT(IN) :: nplin, nplou !list of planets
!  masses of the planets                                                
  DOUBLE PRECISION, INTENT(IN) :: pm(8) 
!  yuasa coefficients--warning!! firt index is yuasa index,             
!    second index is planet index                                       
!  yuasa coefficients  a                                                
  DOUBLE PRECISION, INTENT(IN) :: ya(4,8) 
!  OUTPUT: forced terms
  DOUBLE PRECISION, INTENT(OUT) :: xi(8),et(8)
! END INTERFACE
  DOUBLE PRECISION  pxi(8,8),pet(8,8), divs, divg
  INTEGER j, k
!  for eccentricities                                                   
  DO 51 k=nplin,nplou 
     divg=1.d0/(g0-pg(k))/al 
     DO j=nplin,nplou 
        pxi(j,k)=-pni(j,k)*pm(j)*(ya(4,j))*divg 
     ENDDO
51 ENDDO
!  for inclinations                                                     
  DO 50 k=nplin,nplou 
!         if(k.eq.5)goto 50                                             
     divs=1.d0/(s0-ps(k))/al 
     DO j=nplin,nplou 
        pet(j,k)=-pmi(j,k)*pm(j)*(ya(3,j))*divs 
     ENDDO
50 ENDDO
!  the forced term at a given frequency is the sum of the contributions 
!  from all the planets; that is, mixed terms are included              
  DO 67 k=nplin,nplou 
     xi(k)=0.d0 
     et(k)=0.d0 
     DO j=nplin,nplou 
        xi(k)=xi(k)+pxi(j,k) 
        et(k)=et(k)+pet(j,k)
     ENDDO
67 ENDDO
END SUBROUTINE forcedgs
! **************************************************                    
SUBROUTINE nlforc(nnl,pninl,xinl,nnlx,ya) 
  IMPLICIT NONE
! INPUT
  INTEGER, INTENT(IN) :: nnl,nnlx 
!  nonlinear secular perturbation theory for Jupiter, Saturn            
  DOUBLE PRECISION, INTENT(IN) :: pninl(2,nnl) 
!  yuasa coefficients--warning!! firt index is yuasa index,             
!    second index is planet index                                       
!  yuasa coefficients  a                                                
  DOUBLE PRECISION, INTENT(IN) :: ya(4,8)
! OUTPUT
  DOUBLE PRECISION, INTENT(OUT) :: xinl(nnlx,8)
! END INTERFACE
  INTEGER k
!  nonlinear forced terms                                               
  DO k=1,nnl 
     xinl(k,5)=-pninl(1,k)*ya(4,5) 
     xinl(k,6)=-pninl(2,k)*ya(4,6) 
   ENDDO 
 END SUBROUTINE nlforc
!                                                                       
! *******************************************************************   
!  propmo version 8 (Apr. 1997)                                         
!  proper modes and proper actions                                      
!  that is removal of the forced terms with frequencies indexed         
!  from nplin to nplou                                                  
! *******************************************************************   
SUBROUTINE propmo(ap,aq,ar,as,al                                  &
     &            ,xi,et,pl,pd,nplin,nplou                              &
     &           ,rho)                                                  
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: nplin, nplou !list of planets
  DOUBLE PRECISION, INTENT(IN) :: al, ap, aq, ar, as
  DOUBLE PRECISION, INTENT(IN) :: xi(8),et(8),pl(8),pd(8) 
  DOUBLE PRECISION, INTENT(OUT) :: rho(4)
! END INTERFACE 
  DOUBLE PRECISION sl, fro1, fro2, fsi1, fsi2
  INTEGER j 
  sl=dsqrt(al) 
!  forced eccentricities and inclinations                               
  fro1=0.d0 
  fro2=0.d0 
  fsi1=0.d0 
  fsi2=0.d0 
  DO j=nplin,nplou 
     fro1=fro1+xi(j)*cos(pl(j)) 
     fro2=fro2+xi(j)*sin(pl(j)) 
     fsi1=fsi1+et(j)*cos(pd(j)) 
     fsi2=fsi2+et(j)*sin(pd(j))
  ENDDO
!  proper Poincare' variables                                           
  rho(1)=ap/sl-fro1 
  rho(2)=-aq/sl-fro2 
  rho(3)=ar/sl-fsi1 
  rho(4)=-as/sl-fsi2 
END SUBROUTINE propmo

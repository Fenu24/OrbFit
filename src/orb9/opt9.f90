! =====================================================                 
!  version-dependent options                                            
SUBROUTINE lonopt(inp,ive)
  USE controlmod
  IMPLICIT NONE
  INTEGER, INTENT(IN) ::  inp,ive
!==================================
  INTEGER j                                     
!  version
  if(ive.eq.-1)then 
     vers='9-f90  23/11/2008'
  elseif(ive.eq.-2)then
     vers='9  June 1, 2000' 
  elseif(ive.eq.31)then 
     vers='8.3.1 15/5/97' 
  else 
     write(*,*)' version ',ive,' unknown' 
     stop 
  endif
!  control of the secular perturbation theory to be used:               
!  isec=1 LONGSTOP 1B; isec=2 ORBIT8V                                   
  if(inp.eq.1.or.inp.eq.4)then 
     isec=1 
  elseif(inp.eq.2.or.inp.eq.3)then 
     isec=2 
  else 
     write(*,*)' unknown option inp=',inp 
     stop 
  endif
!  controls for planets taken into account                              
!                                                                       
!  to be accounted for in the forced terms                              
!  (both frequency index and planet index)                              
!  outer planets from Jupiter to Neptune in all versions                
  nplin=5 
  nplou=8 
!  planets to be accounted for in the computation of g0,s0              
!  including inner planets for real objects, excluding for output of num
  if(inp.eq.1.or.inp.eq.4)then 
     nplin0=2 
!     nplin0=5
  elseif(inp.eq.2.or.inp.eq.3)then 
     nplin0=5 
  endif
  nplou0=8 
!  warning; nplin0.le.nplin.le.nplou0.le.nplou                          
!  flags to control first/second order computation                      
  ipla=0
  DO j=nplin0,nplou 
     ipla(j)=1 
  ENDDO
!  in this version second order only Jupiter                            
  ipla(5)=2 
!  planets for which degree four corrections to elements                
!  and frequencies are computed                                         
!  in this version Jupiter and Saturn                                   
  npli4=5 
  nplo4=6 
! flags to remove terms depending upon some divisors                    
!  no divisors removed in this version; see denom5.f for j
  dd(1:ndex)=1.d0 
!  control for divisor rated as small (for automatic elimination        
!  of near resonant terms), in arcsec/yr                                
  divm=0.5d0 
!  critical value of QCE for attempting to switch to adapted theory     
  iqcr=10 
!  controls for convergence (in arcsec/yr)                              
  fconte=1.d-3 
  fconti=1.d-3 
  nitmax=20 
  divrat=1.d0 
END SUBROUTINE lonopt
! **********************************************                        
!  selection of discarded cases to be written elsewhere                 
SUBROUTINE select(inp,aa,pre,prsini,gsy,ssy,iqce,iqcf,iqcm,discard)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: inp
  DOUBLE PRECISION, INTENT(IN):: aa, pre,prsini,gsy,ssy
  INTEGER, INTENT(IN) :: iqce, iqcf, iqcm
  LOGICAL, INTENT(OUT) :: discard 
  if(inp.eq.1.or.inp.eq.4)then 
!  totally wrong cases                                                  
     if(iqcm.gt.13.or.aa.lt.2.d0)then 
        discard=.true. 
!               ipro=3                                                  
!               iang=11                                                 
!  too high inclination or eccentricity                                 
     elseif(prsini.gt.0.3d0.or.pre.gt.0.3.or.gsy.lt.28.24)then
        discard=.true. 
!              ipro=3                                                   
!              iang=11                                                  
     else 
        discard=.false. 
!              ipro=2                                                   
!              iang=9                                                   
     endif
  else 
!  if input from numerical integration, no separate files               
     discard=.false. 
!           ipro=2                                                      
!           iang=9                                                      
  endif
END SUBROUTINE select

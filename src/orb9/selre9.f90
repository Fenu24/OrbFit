! =================================                                     
!   selres (version 8.0)                                                
!   automatic selection of terms to be discarded                        
!     input:                                                            
!         ipas: last iteration                                          
!         nde: number of divisors                                       
!         amp(1:nde,1:ipas): amplitude (total e and I), per div. and per
!         de(nde): divisor                                              
!         da: size of last jump in proper elements                      
!         divm: control on divisors rated small (from controlmod)
!         conv: conversion of frequencies to arcsec/y                   
!     output:                                                           
!         ires=1 removal applied; ires=0 go on; ires=-1 stop iterations 
!         irfl= code for the first resonance to be removed              
!         irfl2= code for the second                                     
! =================================                                     
SUBROUTINE selres(amp,de,da,ipas,irfl,irfl2,conv,ires,       &
     &      nde,nam0,t,inp,nada,nada2,rho0,rho,iun10)                   
  USE controlmod
  IMPLICIT NONE
! INPUT
  INTEGER, INTENT(IN) :: nde ! number of divisors
  DOUBLE PRECISION, INTENT(IN) :: da, t ! amplitude change, time
  CHARACTER*9, INTENT(IN) ::  nam0 ! asteroid name
  INTEGER, INTENT(IN) ::  iun10, inp ! output unit, input mode
! amplitudes, divisors for each combination frequency
  DOUBLE PRECISION, INTENT(IN) :: amp(ndex,maxit),de(ndex,maxit)
  DOUBLE PRECISION, INTENT(IN) :: conv
! linear proper elements
  DOUBLE PRECISION, INTENT(IN) :: rho0(4)
! OUTPUT
  INTEGER, INTENT(INOUT) :: irfl, irfl2, ires !resonance flags
  INTEGER, INTENT(INOUT) :: ipas ! iteration number
  DOUBLE PRECISION, INTENT(INOUT) :: rho(4) ! proper elements
! END INTERFACE
  INTEGER i, jmax, jmax1, nada,nada2
  DOUBLE PRECISION ampx, ampx1
!  safety factor
  DOUBLE PRECISION saf                                                        
  data saf/3.d0/ 
!  if resonances already selected, what do you want?                    
  if(irfl2.ne.0)then 
     ires=-1 
     return 
  endif
!   look for the largest amplitude change, only among the ones with     
!   not too big a divisor (control divm), now                           
  call fmaxi(amp(1,ipas),1,nde,de(1,ipas),conv                     &
     &          ,divm,irfl,jmax,ampx)                                   
  if(ipas.gt.1)then 
     call fmaxi(amp(1,ipas-1),1,nde,de(1,ipas-1),conv              &
     &         ,divm,irfl,jmax1,ampx1)                                  
  endif
!   check whether the largest amplitude is large enough to be           
!   responsible for the jump                                            
  if(ampx*saf.gt.da)then 
     if(abs(de(jmax,ipas)*conv).lt.divm)then 
        if(irfl.eq.0)then 
           irfl=jmax 
        else 
           irfl2=jmax 
        endif
        ires=jmax 
     endif
  else 
!   try with the previous iteration                                     
     if(ipas.gt.1)then 
        if(ampx1*saf.gt.da)then 
           if(abs(de(jmax1,ipas-1)*conv).lt.divm)then 
              if(irfl.eq.0)then 
                 irfl=jmax1 
              else 
                 irfl2=jmax1 
              endif
              ires=jmax1 
           endif
        else 
           ires=0 
        endif
     else 
        ires=0 
     endif
  endif
!  added operations                                                     
  if(ires.gt.0)then 
!  write warning in resonance file                                      
     if(inp.eq.1.or.inp.eq.4)then 
        write(iun10,110)nam0,ipas,ires,amp(ires,ipas)              &
     &               ,de(ires,ipas)*conv                                
110     format(' removal',a9,i3,i3,1p,e12.4,e12.4) 
     else 
        write(iun10,120)t,ipas,ires,amp(ires,ipas)                 &
     &               ,de(ires,ipas)*conv                                
120     format(' removal',f12.2,i3,i3,1p,e12.4,e12.4) 
     endif
!  counters of adapted cases                                            
     if(irfl2.eq.0)then 
        nada=nada+1 
     else 
        nada2=nada2+1 
     endif
!  restart from linear proper elements                                  
     do  i=1,4 
        rho(i)=rho0(i) 
     enddo
!  iteration restarts; divergence counter reset to zero                 
     ipas=0 
  endif
END SUBROUTINE selres
! =============                                                         
!   fmaxi                                                               
SUBROUTINE fmaxi(a,i1,i2,de,conv,divm,j1,jx,ax) 
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: i1,i2,j1
  DOUBLE PRECISION, INTENT(IN) :: a(i2),de(i2), conv,divm 
  DOUBLE PRECISION, INTENT(OUT) :: ax
  INTEGER, INTENT(OUT) :: jx
  INTEGER j
  ax=0.d0 
  jx=0 
  DO 1 j=i1,i2 
     if(abs(a(j)).gt.ax.and.j.ne.j1)then 
        if(abs(de(j))*conv.le.divm)then 
           jx=j 
           ax=a(j) 
        endif
     endif
1 ENDDO
END SUBROUTINE fmaxi
! *************************************************************         
!                                                                       
!   d e n o m (vers. 6.8 and later)                                     
!   computes divisors for given frequencies;                            
!   includes doubles, without the factor 2                              
!   includes divisors 60, 61, 62, and cross reference for               
!   forced nonlinear frequencies                                        
!                                                                       
! *************************************************************         
SUBROUTINE denom (g,s,pg,ps,de,nde,inl) 
  USE controlmod
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT(IN) :: g,s ! proper frequencies
!  secular theory frequencies for the planets g_j, s_j                  
  DOUBLE PRECISION, INTENT(IN) :: pg(8),ps(8) 
!  divisors, actual number                                           
  DOUBLE PRECISION, INTENT(OUT) :: de(ndex) 
  INTEGER, INTENT(OUT) :: nde
!  cross reference indexes for forced nonlinear terms, actual number
  INTEGER, INTENT(OUT) :: inl(nnlx)
! END INTERFACE
  DOUBLE PRECISION g5,g6,s7,g7, s6  
  g5=pg(5) 
  g6=pg(6) 
  s6=ps(6) 
  s7=ps(7) 
  g7=pg(7) 
  de(1)=g-g5 
  de(2)=g-g6 
! this one cannot give rise to a resonance                              
  de(3)=g5-g6 
  de(4)=s-s7 
  de(5)=s-s6 
! this one cannot give rise to a resonance                              
  de(6)=s7-s6 
  de(7)=g+s-s7-g5 
  de(8)=g+s-s7-g6 
  de(9)=g+s-s6-g5 
  de(10)=g+s-s6-g6 
  de(11)=2.d0*g-2.d0*s 
  de(12)=g-2.d0*g5+g6 
  de(13)=g+g5-2.d0*g6 
  inl(1)=13 
  de(14)=2.d0*g-g5-g6 
  de(15)=-g+s+g5-s7 
  de(16)=-g+s+g6-s7 
  de(17)=-g+s+g5-s6 
  de(18)=-g+s+g6-s6 
  de(19)=g-g5+s7-s6 
  de(20)=g-g5-s7+s6 
  de(21)=g-g6+s7-s6 
  de(22)=g-g6-s7+s6 
  de(23)=2.d0*g-s-s7 
  de(24)=2.d0*g-s-s6 
  de(25)=-g+2.d0*s-g5 
  de(26)=-g+2.d0*s-g6 
  de(27)=2.d0*g-2.d0*s7 
  de(28)=2.d0*g-2.d0*s6 
  de(29)=2.d0*g-s7-s6 
  de(30)=g-s+g5-s7 
  de(31)=g-s+g5-s6 
  de(32)=g-s+g6-s7 
  de(33)=g-s+g6-s6 
  de(34)=g+g5-2.d0*s7 
  de(35)=g+g6-2.d0*s7 
  de(36)=g+g5-2.d0*s6 
  de(37)=g+g6-2.d0*s6 
  de(38)=g+g5-s7-s6 
  de(39)=g+g6-s7-s6 
  de(40)=s-2.d0*s7+s6 
  de(41)=s+s7-2.d0*s6 
  de(42)=2.d0*s-s7-s6 
  de(43)=s+g5-g6-s7 
  de(44)=s-g5+g6-s7 
  de(45)=s+g5-g6-s6 
  de(46)=s-g5+g6-s6 
  de(47)=2.d0*s-2.d0*g5 
  de(48)=2.d0*s-2.d0*g6 
  de(49)=2.d0*s-g5-g6 
  de(50)=s-2.d0*g5+s7 
  de(51)=s-2.d0*g5+s6 
  de(52)=s-2.d0*g6+s7 
  de(53)=s-2.d0*g6+s6 
  de(54)=s-g5-g6+s7 
  de(55)=s-g5-g6+s6 
  de(56)=2.d0*g-2.d0*g5 
  de(57)=2.d0*g-2.d0*g6 
  de(58)=2.d0*s-2.d0*s7 
  de(59)=2.d0*s-2.d0*s6 
!  divisors appearing only in forced terms                              
  de(60)=g-2.d0*g6+g7 
  inl(2)=60 
  de(61)=g-3.d0*g6+2.d0*g5 
  inl(3)=61 
!  degree six divisor z2                                                
  de(62)=2.d0*(g-g6)+(s-s6) 
!  other nonlinear forced terms                                         
  de(63)=g+g5-g6-g7 
  inl(4)=63 
  de(64)=g-g5-g6+g7 
  inl(5)=64 
  de(65)=g+g5-2*g6-s6+s7 
  inl(6)=65 
  inl(7)=20 
  inl(8)=12 
!  total number of divisors appearing in the theory                     
  nde=65 
!  divisors not appearing at all in the theory (wait for next vers.)    
  de(71)=3.d0*(g-g6)+(s-s6)   
END SUBROUTINE denom

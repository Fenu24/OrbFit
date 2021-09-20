MODULE short9
  USE fund_const
  USE massmod
  USE controlmod
  IMPLICIT NONE

  PRIVATE
! shared data

! controls                                                              
  DOUBLE PRECISION eps,epdiv,epsdel
!  control on eccentricity/inclination considered zero                  
  DOUBLE PRECISION, PARAMETER :: epsil=1.d-6 
! generic terms of the hamiltonian                                      
  INTEGER ntd,nti,ih(ntx,17),ihi(ntx,15) 
! routines, public
 PUBLIC short

CONTAINS


! *********************************************************             
!  {\bf short} version 8 march 27, 1997                                 
!   short calculates the short-periodic perturbations                   
!     in a, e, i, o, wt.                                                
!   input from file oscel. output to file meanel.                       
!   revision history:                                                   
!     version 6.0 with generic terms                                    
!     version 6.0.1 trigonometric computations replaced by addition form
!     version 6.1.0 with input routines and orbit8v formats;            
!     version 6.1.1 with variable number of terms in the hamiltonian fil
!     version 6.2.1 with Laplace and Leverrier coefficients as in longit
!     version 6.3.1 with computation of mean mean anomaly; also format  
!          in input changed to allow for the asteroid name              
!     version 6.3.3 to be used only for asteroids with up to date       
!           orbit                                                       
!     version 6.3.4 to be used for asteroids with elements at variable  
!           epoch                                                       
!     version 7.0 unifyes 6.3.3 and 6.3.4 (inp=3 for 6.3.4, inp=1,2     
!                 for 6.3.3)                                            
!                 also provides for asteroid name as 6 char string      
!                 as an alternative to 8 digits code number (only in out
!     version 7.7 as 7.0 but with new software                          
!     version 8.0 adds second order corrections with RKG order 2        
!     versions 8.0-8.1-8.2-8.3-8.4-8.5-8.6                              
!   input:                                                              
!     any reference and coordinate system, as specified in the input hea
!   output:                                                             
!     heliocentric keplerian elements, invariable plane LONGSTOP 1b     
! *********************************************************             
SUBROUTINE short(el,pel,nam0,tt,inp,em,iqc,it,delta,ierr,        &
     &     iun19,iun21)                                                 
! =============================================================         
! INPUT                                                                 
  INTEGER inp 
! asteroid number (inp=1,3) or time (inp=2)                             
  DOUBLE PRECISION tt 
  CHARACTER*9 nam0 
! el,em: 1=semimaj. ax. 2= ecc. 3=inclin.                               
!        4=node 5=arg. peri 6=mean anomaly                              
!  planetary elements; second index is 1=jup 2=sat                      
  DOUBLE PRECISION pel(6,2) 
!  asteroid keplerian elements: osculating                              
  DOUBLE PRECISION el(6) 
!  unit for output of warnings and headers                              
  INTEGER iun21,iun19 
! =============================================================         
! OUTPUT                                                                
!  asteroid keplerian elements: mean, quality code                      
  DOUBLE PRECISION em(6) 
  INTEGER iqc 
! iteration results                                                     
  DOUBLE PRECISION delta 
  INTEGER it 
! error control                                                         
  INTEGER ierr 
! END INTERFACE                                                         
! =============================================================         
!  intermediate values                                                  
  DOUBLE PRECISION elz(6) 
!  asteroid names (alphanumeric)                                        
  character*6 astna 
! masses etc.                                                           
  DOUBLE PRECISION h,h2 
! control flags                                                         
  INTEGER integr 
! iteration control                                                     
  INTEGER itmax 
!  version identifier                                                   
  character*60 vers 
! loop indexes                                                          
  INTEGER i 
! functions                                                             
  DOUBLE PRECISION princ 
! =============================================================         
! initialisation                                                        
  INTEGER flag 
  SAVE 
  DATA flag/0/ 
  IF(flag.eq.0)THEN 
!  files with  Hamiltonian: opening and reading                         
     call iniham 
!  step is mass ratio Jupiter/sun                                       
     h=gm(2) 
! private options file, not known to the user                           
     open(1,file='short8s.opt',status='old') 
     call skip (1,1) 
!  controls for resonances                                              
     call reaflo(1,'eps',eps) 
     call reaflo(1,'epdiv',epdiv) 
! *************************************************************         
!  choice of method to integrate Lie series equation                    
     call reaint(1,'integr',integr) 
     call reaflo(1,'epsdel',epsdel) 
     if(integr.eq.-1)then 
        vers=' 7.7, (Euler=7.0) Pisa  Nov 1996' 
     elseif(integr.eq.0)then 
        vers=' 8.0 (RKG2) Pisa, Dec 1996' 
     elseif(integr.eq.1)then 
        vers=' 8.1 (RK2) Pisa, Apr 1997' 
     elseif(integr.eq.2)then 
        vers=' 8.2 (RKG2+reverse) Pisa, Dec 1996' 
     elseif(integr.eq.3)then 
        vers=' 8.3 (RK2+reverse) Pisa, Dec 1996' 
     elseif(integr.eq.4)then 
        vers=' 8.4 (RK2-Aitken) Pisa, Apr 1997' 
     elseif(integr.eq.5)then 
        vers=' 8.5 (RKG2A-2step)Pisa, Apr 1997' 
     elseif(integr.eq.6)then 
        vers=' 8.6 (RKG4) Pisa, Dec 1996' 
     elseif(integr.eq.31)then 
        vers=' 8.3.1 (Euler+reverse) Pisa, May 1997' 
     else 
        write(*,*)' unknown option integr=',integr 
        stop 
     endif
!  header resonance warning file                                        
     write(iun21,211)vers 
211  format(                                                           &
     &' Mean motion res, zero ecc/incl, iterat failure warnings'/       &
     &   ' vers. ',a60/                                                 &
     & ' ast.no, i+b13, i+b12, b15, b17, b14, b16, cfac, per(yr) - or'/ &
     & ' zero e/sI warning,procedure, no(tt), mean elements - or'/      &
     &        ' ast.no,  iter,  delta a - or'/                          &
     &        ' iteration failure warning'/)                            
! *************************************************************         
!  output header                                                        
     if(inp.eq.1)then 
        write(iun19,232)vers 
232     format('  vers.',a60/                                       &
     &   'ast.no.   m.ano.    arg.per. long.nod.   inclin.    ecc.'     &
     &   ,'     sem.ax.  QCM IER QCO  magn.')                           
     elseif(inp.eq.2)then 
        write(iun19,233)vers 
233     format('  vers.',a60/                                       &
     &   'time,y  m.ano. arg.per. long.nod.  inclin.    ecc. ',         &
     &   '   sem.ax.  QCM IER')                                         
     else 
        write(*,*)' wrong input type in ios63, inp=',inp 
        stop 
     endif
     close(1) 
! end initialisation                                                    
     flag=1 
  ENDIF
! =============================================================         
!  choice of integration method for Lie series                          
! *********************************************************             
  if(integr.eq.-1)then 
! *********************************************************             
! version 7.7 (=7.0) Euler method:                                      
     call euler(el,pel,nam0,tt,inp,h,em,iqc,iun21) 
! *********************************************************             
  elseif(integr.eq.0)then 
! *********************************************************             
!  version 8.0: RKG2 integration step                                   
     itmax=10 
     call rkg2(el,pel,nam0,tt,inp,h,itmax,epsdel,em,iqc,         &
     &        it,delta,iun21)                                           
! *********************************************************             
  elseif(integr.eq.1)then 
! *********************************************************             
! version 8.1: RK explicit  integration step                            
     call rk2(el,pel,nam0,tt,inp,h,em,iqc,ierr,iun21) 
  elseif(integr.eq.2)then 
! *********************************************************             
!  version 8.2: RKG2 integration step is starting point for             
!      iteration of reverse map                                         
! RKG method:                                                           
     itmax=10 
     call rkg2(el,pel,nam0,tt,inp,h,itmax,epsdel,elz,iqc,        &
     &        it,delta,iun21)                                           
! *********************************************************             
!  eccentricity/inclination control                                     
     if(elz(2).lt.epsil.or.elz(3).lt.epsil)then 
        DO i=1,6 
           em(i)=elz(i) 
        ENDDO
        ierr=-2 
        goto 33 
     endif
! *********************************************************             
!  iteration of reverse map                                             
     itmax=10 
     call revers(elz,el,pel,nam0,tt,inp,h,itmax,epsdel,          &
     &            em,iqc,it,delta,ierr,iun21)                           
! ********************************************************              
  elseif(integr.eq.3)then 
! *********************************************************             
! version 8.3: RK explicit  integration step is starting point          
! for reverse map method:                                               
     call rk2(el,pel,nam0,tt,inp,h,elz,iqc,ierr,iun21) 
! *********************************************************             
!  eccentricity/inclination control                                     
     if(elz(2).lt.epsil.or.elz(3).lt.epsil)then 
        DO i=1,6 
           em(i)=elz(i) 
        ENDDO
        ierr=-2 
        goto 33 
     endif
! *********************************************************             
! iteration of reverse map                                              
     itmax=10 
     call revers(elz,el,pel,nam0,tt,inp,h,itmax,epsdel,          &
     &         em,iqc,it,delta,ierr,iun21)                              
! ********************************************************              
  elseif(integr.eq.31)then 
! *********************************************************             
! version 8.3.1: Euler integration step is starting point               
! for reverse map method:                                               
! *********************************************************             
     call euler(el,pel,nam0,tt,inp,h,elz,iqc,iun21) 
!  eccentricity/inclination control                                     
     if(elz(2).lt.epsil.or.elz(3).lt.epsil)then 
        DO i=1,6 
           em(i)=elz(i) 
        ENDDO
        ierr=-1 
        goto 33 
     endif
! *********************************************************             
! iteration of reverse map                                              
     itmax=10 
     call revers(elz,el,pel,nam0,tt,inp,h,itmax,epsdel,          &
     &         em,iqc,it,delta,ierr,iun21)                              
! *********************************************************             
  elseif(integr.eq.4)then 
! *********************************************************             
! version 8.4 RKG method with Aitken accelerator                        
     call rkg2ai(el,pel,nam0,tt,inp,h,em,iqc,it,delta,iun21) 
! *********************************************************             
  elseif(integr.eq.5)then 
! *********************************************************             
! version 8.5: two half steps with RKG-Aitken                           
     h2=h/2.d0 
     call rkg2ai(el,pel,nam0,tt,inp,h2,elz,iqc,it,delta,iun21) 
     call rkg2ai(elz,pel,nam0,tt,inp,h2,em,iqc,it,delta,iun21) 
! ********************************************************              
  elseif(integr.eq.6)then 
! *********************************************************             
! version 8.6: RKG4 implicit method                                     
     call rkg4(el,pel,nam0,tt,inp,h,itmax,epsdel,em,iqc,it,delta,&
     &        iun21)                                                    
! ********************************************************              
  endif
! *********************************************************             
!  control for zero eccentricity/inclination                            
! *********************************************************             
33 continue 
!  eccentricity/inclination control                                     
  if(em(2).lt.epsil.or.em(3).lt.epsil)then 
     if(inp.eq.2)then 
        write(iun21,*)'zero ecc/incl main',tt, em 
     elseif(inp.eq.1.or.inp.eq.3)then 
        write(iun21,*)'zero ecc.incl main',nam0, em 
     endif
     if(em(2).lt.epsil)then 
        em(2)=0.d0 
     endif
     if(em(3).lt.epsil)then 
        em(3)=0.d0 
     endif
  endif
! *********************************************************             
!  conversion to principal value of angular variables                   
  em(4)=princ(em(4)) 
  em(5)=princ(em(5)) 
  em(6)=princ(em(6)) 
! ***************************************************************       
!  angles in degrees                                                    
  DO i=3,6 
     el(i)=el(i)*degrad 
     em(i)=em(i)*degrad 
  ENDDO
! go back to osculating value in all cases in which the short           
! periodic perturbation is unreasonably huge;                           
! zero eccentricity cases and divergent cases (ierr.ne.0) are not reset 
! to osculating values                                                  
  if(iqc.ge.15)then 
     em(1:6)=el(1:6) 
  endif

END SUBROUTINE short

! ===========================================                           
!   euler method of order 1                                             
! -------------------------------------------                           
subroutine euler(el,pel,nam0,tt,inp,h,em,iqc,iun21) 
  CHARACTER*9 nam0 
  INTEGER iun21, inp,iqc 
  DOUBLE PRECISION :: tt, h
!  planetary elements; second index is 1=jup 2=sat                      
  DOUBLE PRECISION ::  pel(6,2) 
!  asteroid elements (osculating, mean; keplerian, Delaunay)            
  DOUBLE PRECISION ::  el(6),em(6),del(6),dem(6) 
!  corrections to Delaunay elements                                     
  DOUBLE PRECISION ::  dd(6) 
! -------------------------------------------                           
!  delaunay variables l,g,h - osculating values                         
  call kepdel(el,del) 
! computation of corrections in lghwo;                                  
  call dercor(dd,iqc,pel,el,inp,nam0,tt,iun21) 
! -------------------------------------------                           
  dem(1:6)=del(1:6)+h*dd(1:6)
!  conversion to mean keplerian                                         
  call delkep(dem,em) 
END subroutine euler
! ===========================================                           
!   Runge Kutta Gauss  method of order 2                                
!   with Aitken's accelerated convergence (3 steps only)                
! -------------------------------------------                           
subroutine rkg2ai(el,pel,nam0,tt,inp,h,em,iqc,it,delta,iun21) 
  CHARACTER*9 nam0 
  INTEGER iun21, inp,iqc 
  DOUBLE PRECISION :: tt, h
!  planetary elements; second index is 1=jup 2=sat                      
  DOUBLE PRECISION ::  pel(6,2) 
!  asteroid elements (osculating, mean; keplerian, Delaunay)            
  DOUBLE PRECISION ::  el(6),em(6),del(6),dem(6) 
!  workspace:                                                           
!  intermediate values for implicit method                              
  DOUBLE PRECISION ::  dez(6,5),elz(6,5) 
!  corrections to Delaunay elements                                     
  DOUBLE PRECISION ::  dd(6,4) 
!  control on eccentricity/inclination considered zero                  
  DOUBLE PRECISION epsil, delta
  data epsil/1.d-6/ 
  INTEGER it,i
! -------------------------------------------                           
!  delaunay variables l,g,h - osculating values                         
  call kepdel(el,del) 
!  starting point is osculating value                                   
  do 5 i=1,6 
     elz(i,1)=el(i) 
     dez(i,1)=del(i) 
5 enddo
!  iteration loop                                                       
  do 1 it=1,3 
! computation of corrections in lghwo;                                  
     call dercor(dd(1,it),iqc,pel,elz(1,it),inp,nam0,tt,iun21) 
! new intermediate point                                                
     do 6 i=1,6 
        dez(i,it+1)=del(i)+(h/2.d0)*dd(i,it) 
6    enddo
     call delkep(dez(1,it+1),elz(1,it+1)) 
     if(elz(2,it+1).lt.epsil.or.elz(3,it+1).lt.epsil)then 
        if(inp.eq.2)then 
           write(iun21,*)'zero ecc/incl rkg2ai',tt,elz 
        elseif(inp.eq.1.or.inp.eq.3)then 
           write(iun21,*)'zero ecc/incl rkg2ai',nam0,elz 
        endif
        goto 9 
     endif
1 enddo
! Aitken's method to accelerate convergence                             
  do 7 i=1,6 
     dd(i,4)=dd(i,1)-                                                &
     &   (dd(i,2)-dd(i,1))**2/(dd(i,3)-2.d0*dd(i,2)+dd(i,1))            
     dez(i,5)=del(i)+(h/2.d0)*dd(i,4) 
7 enddo
  call delkep(dez(1,5),elz(1,5)) 
! control on last change                                                
9 delta=abs(elz(1,it+1)-elz(1,it)) 
! computation of mean value                                             
  do 8 i=1,6 
     dem(i)=del(i)+h*dd(i,it) 
8 enddo
!  conversion to mean keplerian                                         
  call delkep(dem,em) 
! -------------------------------------------                           
END subroutine rkg2ai
! ===========================================                           
!   Runge Kutta Gauss  method of order 2                                
!   with convergence control on semimajor axis                          
! -------------------------------------------                           
subroutine rkg2(el,pel,nam0,tt,inp,h,itmax,epsdel,em,iqc,      &
     &    it,delta,iun21)                                               
  CHARACTER*9 nam0 
  INTEGER iun21, inp,iqc, itmax 
  DOUBLE PRECISION :: tt, h, epsdel
!  planetary elements; second index is 1=jup 2=sat                      
  DOUBLE PRECISION ::  pel(6,2) 
!  asteroid elements (osculating, mean; keplerian, Delaunay)            
  DOUBLE PRECISION ::  el(6),em(6),del(6),dem(6) 
!  workspace:                                                           
  INTEGER, PARAMETER :: nitmax=20 
!  intermediate values for implicit method                              
  DOUBLE PRECISION ::  dez(6,nitmax),elz(6,nitmax) 
!  corrections to Delaunay elements                                     
  DOUBLE PRECISION ::  dd(6,nitmax) 
!  control on eccentricity/inclination considered zero                  
  DOUBLE PRECISION epsil, delta
  data epsil/1.d-6/ 
  INTEGER it,i
! -------------------------------------------                           
!  delaunay variables l,g,h - osculating values                         
  call kepdel(el,del) 
!  starting point is osculating value                                   
  do 5 i=1,6 
     elz(i,1)=el(i) 
     dez(i,1)=del(i) 
5 enddo
!  iteration loop                                                       
  do 1 it=1,itmax 
! computation of corrections in lghwo;                                  
     call dercor(dd(1,it),iqc,pel,elz(1,it),inp,nam0,tt,iun21) 
! new intermediate point                                                
     do 6 i=1,6 
        dez(i,it+1)=del(i)+(h/2.d0)*dd(i,it) 
6    enddo
     call delkep(dez(1,it+1),elz(1,it+1)) 
     delta=abs(elz(1,it+1)-elz(1,it)) 
     if(delta.lt.epsdel)goto 9 
     if(elz(2,it+1).lt.epsil.or.elz(3,it+1).lt.epsil)then 
        if(inp.eq.2)then 
           write(iun21,*)'zero ecc/incl rkg2',tt,elz 
        elseif(inp.eq.1.or.inp.eq.3)then 
           write(iun21,*)'zero ecc/incl rkg2',nam0,elz 
        endif
        goto 9 
     endif
1 enddo
  write(*,*)' RKG iteration failure ', delta 
  write(iun21,*)' RKG iteration failure ' 
9 continue 
! computation of mean value                                             
  do 8 i=1,6 
     dem(i)=del(i)+h*dd(i,it) 
8 enddo
!  conversion to mean keplerian                                         
  call delkep(dem,em) 
! -------------------------------------------                           
END subroutine rkg2
! ===========================================                           
!   Runge Kutta explicit method of order 2                              
! -------------------------------------------                           
subroutine rk2(el,pel,nam0,tt,inp,h,em,iqc,ierr,iun21) 
  CHARACTER*9 nam0 
  INTEGER iun21, inp,iqc, ierr 
  DOUBLE PRECISION :: tt, h
!  planetary elements; second index is 1=jup 2=sat                      
  DOUBLE PRECISION ::  pel(6,2) 
!  asteroid elements (osculating, mean; keplerian, Delaunay)            
  DOUBLE PRECISION ::  el(6),em(6),del(6),dem(6) 
!  intermediate values                                                  
  DOUBLE PRECISION ::  dez(6),elz(6) 
!  corrections to Delaunay elements                                     
  DOUBLE PRECISION ::  dd(6),ddz(6) 
!  control on eccentricity/inclination considered zero                  
  DOUBLE PRECISION epsil
  data epsil/1.d-6/ 
  INTEGER i
! -------------------------------------------                           
!  delaunay variables l,g,h - osculating values                         
  call kepdel(el,del) 
! computation of corrections in lghwo;                                  
  call dercor(dd,iqc,pel,el,inp,nam0,tt,iun21) 
! emergency exit                                                        
  ierr=0 
  if(iqc.gt.15)then 
     write(iun21,*)' too close to resonance, iqc= ',iqc 
     write(iun21,*)dd 
     ierr=1 
  endif
!  calculation of intermediate value                                    
  DO i=1,6 
     dez(i)=del(i)+h*dd(i)
  ENDDO 
  call delkep(dez,elz) 
  if(elz(2).lt.epsil.or.elz(3).lt.epsil)then 
!  zero ecc/incl: exit with euler                                       
     if(inp.eq.2)then 
        write(iun21,*)'zero ecc/incl rk2',tt,elz 
     elseif(inp.eq.1.or.inp.eq.3)then 
        write(iun21,*)'zero ecc/incl rk2',nam0,elz 
     endif
     ierr=-1 
     DO i=1,6 
        em(i)=elz(i) 
     ENDDO
     return 
  endif
8 continue 
! computation of corrections at intermediate value;                     
  call dercor(ddz,iqc,pel,elz,inp,nam0,tt,iun21) 
! emergency exit                                                        
  if(iqc.gt.15)then 
     write(iun21,*)' too close to resonance, iqc= ',iqc 
     write(iun21,*)dd 
     ierr=2 
  endif
!  starting point for iteration of reverse map                          
  DO i=1,6 
     dem(i)=del(i)+0.5d0*h*(dd(i)+ddz(i)) 
  ENDDO
  call delkep(dem,em) 
END subroutine rk2
! ===========================================                           
!   iteration of inverse map                                            
!   with convergence control on semimajor axis                          
! -------------------------------------------                           
subroutine revers(elz,el,pel,nam0,tt,inp,h,itmax,epsdel,       &
     &         em,iqc,it,delta,ierr,iun21) 
  CHARACTER*9 nam0 
  INTEGER iun21, inp,iqc, itmax, ierr 
  DOUBLE PRECISION :: tt, h, epsdel
!  planetary elements; second index is 1=jup 2=sat                      
  DOUBLE PRECISION ::  pel(6,2) 
!  asteroid elements (osculating, mean; keplerian, Delaunay)            
  DOUBLE PRECISION ::  el(6),em(6),del(6),dem(6) 
!  workspace:                                                           
  INTEGER, PARAMETER :: nitmax=20 
!  intermediate values for implicit method                              
  DOUBLE PRECISION ::  elz(6),elzold(6) 
!  corrections to Delaunay elements                                     
  DOUBLE PRECISION ::  dd(6) 
!  control on eccentricity/inclination considered zero
  DOUBLE PRECISION epsil, delta
  data epsil/1.d-6/ 
  INTEGER i, it
! -------------------------------------------                           
!  delaunay variables l,g,h - osculating values                         
  call kepdel(el,del) 
!  iteration loop                                                       
  ierr=0 
  do 7 it=1,itmax 
! computation of corrections in lghwo;                                  
     call dercor(dd,iqc,pel,elz,inp,nam0,tt,iun21) 
! emergency exit                                                        
     if(iqc.gt.15)then 
        write(iun21,*)' too close to resonance, iqc= ',iqc 
        write(iun21,*)dd 
        ierr=it+2 
        goto 8 
     endif
!  calculation of mean value next iteration                             
     DO i=1,6 
        elzold(i)=elz(i) 
        dem(i)=del(i)+h*dd(i) 
     ENDDO
     call delkep(dem,elz) 
! iteration control                                                     
     delta=abs(elz(1)-elzold(1)) 
!            write(iun21,*)it,delta                                     
     if(delta.lt.epsdel)goto 8 
     if(elz(2).lt.epsil.or.elz(3).lt.epsil)then 
        if(inp.eq.2)then 
           write(iun21,*)'zero ecc/incl revers',tt,elz 
        elseif(inp.eq.1.or.inp.eq.3)then 
           write(iun21,*)'zero ecc/incl revers',nam0,elz 
        endif
        ierr=-it-2 
        goto 8 
     endif
7 ENDDO
  write(*,*)' inverse map iteration failure ', delta 
  write(iun21,*)' inverse map iteration failure ' 
  ierr=it+2 
8 continue 
!  conversion to mean keplerian 
  em(1:6)=elz(1:6) 
END subroutine revers
! ===========================================                           
!   Runge Kutta Gauss  method of order 4                                
!   with convergence control on semimajor axis                          
! -------------------------------------------                           
subroutine rkg4(el,pel,nam0,tt,inp,h,itmax,epsdel,em,iqc,      &
     &   it,delta,iun21)  
  CHARACTER*9 nam0 
  INTEGER iun21, inp,iqc, itmax
  DOUBLE PRECISION :: tt, h, epsdel
!  planetary elements; second index is 1=jup 2=sat                      
  DOUBLE PRECISION ::  pel(6,2) 
!  asteroid elements (osculating, mean; keplerian, Delaunay)            
  DOUBLE PRECISION ::  el(6),em(6),del(6),dem(6) 
!  workspace:                                                           
  INTEGER,PARAMETER ::  nitmax=20 
!  intermediate values for implicit method                              
  DOUBLE PRECISION ::  dez(6,2,nitmax),elz(6,2,nitmax) 
!  corrections to Delaunay elements                                     
  DOUBLE PRECISION ::  dd(6,2,nitmax) 
!  coeffients for RKG4                                                  
  DOUBLE  PRECISION ::  a(2,2),b(2),c(2) 
!  control on eccentricity/inclination considered zero
  DOUBLE  PRECISION epsil, delta       
  data epsil/1.d-6/ 
  INTEGER it,i,k
! -------------------------------------------                           
! RKG coefficients                                                      
  c(1)=(3.d0-sqrt(3.d0))/6.d0 
  c(2)=(3.d0+sqrt(3.d0))/6.d0 
  b(1)=0.5d0 
  b(2)=0.5d0 
  a(1,1)=0.25d0 
  a(2,2)=0.25d0 
  a(1,2)=(3.d0-sqrt(1.2d1))/1.2d1 
  a(2,1)=(3.d0+sqrt(1.2d1))/1.2d1 
! -------------------------------------------                           
!  delaunay variables l,g,h - osculating values                         
  call kepdel(el,del) 
!  starting from osculating value                                       
  DO k=1,2 
     elz(1:6,k,1)=el(1:6) 
     dez(1:6,k,1)=del(1:6)
  ENDDO
! computation of corrections at osculating value                        
  call dercor(dd(1,1,1),iqc,pel,el,inp,nam0,tt,iun21) 
! emergency exit                                                        
  if(iqc.gt.20)then 
     write(iun21,*)' too close to resonance, iqc= ',iqc 
     write(iun21,*)dd(6,1,1) 
     goto 88 
  endif
! -------------------------------------------                           
!  first guess of intermediate value (by Euler method)                  
  DO k=1,2 
     DO  i=1,6 
        dez(i,k,1)=del(i)+c(2)*h*dd(i,1,1) 
     ENDDO
     call delkep(dez(1,k,1),elz(1,k,1)) 
  ENDDO
! -------------------------------------------                           
!  iteration loop                                                       
  itmax=10 
  do 8 it=1,itmax 
! computation of corrections in lghwo;                                  
     DO k=1,2 
        call dercor(dd(1,k,it+1),iqc,pel,elz(1,k,it),inp,nam0,tt,iun21) 
     ENDDO
! emergency exit                                                        
     if(iqc.gt.20)then 
        write(iun21,*)' too close to resonance, iqc= ',iqc 
        write(iun21,*)dd(1,2,it+1) 
        em(1:6)=el(1:6) 
        return 
     endif
!  calculation of intermediate value                                    
     DO k=1,2 
        DO i=1,6 
           dez(i,k,it+1)=                                          &
     &           del(i)+h*(a(k,1)*dd(i,k,it+1)+a(k,2)*dd(i,k,it+1))     
        ENDDO
        call delkep(dez(1,k,it+1),elz(1,k,it+1)) 
     ENDDO
! -------------------------------------------                           
! iteration control                                                     
     delta=(abs(elz(1,1,it+1)-elz(1,1,it))+                       &
     &         abs(elz(1,2,it+1)-elz(1,2,it)))*0.5d0                    
!          write(iun21,*)it,delta                                       
     if(delta.lt.epsdel)goto 88 
     if(elz(2,1,it+1).lt.epsil.or.elz(3,1,it+1).lt.epsil)then 
        if(inp.eq.2)then 
           write(iun21,*)'zero ecc/incl, point 2 ',tt,elz 
        elseif(inp.eq.1.or.inp.eq.3)then 
           write(iun21,*)'zero ecc/incl, point 2 ',nam0,elz 
        endif
        goto 88 
     elseif(elz(2,2,it+1).lt.epsil.or.elz(3,2,it+1).lt.epsil)then 
        write(iun21,*)'zero ecc/incl, point 1', elz 
        if(inp.eq.2)then 
           write(iun21,*)'zero ecc/incl, point 1 ',tt,elz 
        elseif(inp.eq.1.or.inp.eq.3)then 
           write(iun21,*)'zero ecc/incl, point 1 ',nam0,elz 
        endif
        goto 88 
     endif
8 ENDDO
! -------------------------------------------                           
!  exit from main loop: failure                                         
  write(*,*)' RKG iteration failure ', delta 
  write(iun21,*)' RKG iteration failure ' 
! -------------------------------------------                           
!  exit from main loop: success                                         
88 continue 
  DO i=1,6 
     dem(i)=del(i)+h*(b(1)*dd(i,1,it+1)+b(2)*dd(i,2,it+1)) 
  ENDDO
!  conversion to mean keplerian                                         
  call delkep(dem,em) 
END subroutine rkg4

! *********************************************************             
!   subroutine iniham                                                   
!  initialises the integer table of hamiltonian terms                   
! max number of terms in generic Hamiltonian (degree 4)                 
SUBROUTINE iniham 
! input units
  INTEGER iun60,iun61
! terms counters
  INTEGER  m,nt, n 
!  files with  Hamiltonian: opening and reading                         
  call filopn(iun60,'hamil.dir','old') 
  call filopn(iun61,'hamil.ind','old') 
  read(iun60,399) ntd 
399 format(/i6) 
  DO m=1,ntd 
     read(iun60,400)(ih(m,n),n=1,17)
400  format(17i4) 
  ENDDO
  read(iun61,399) nti 
  DO m=1,nti 
     read(iun61,405)(ihi(m,n),n=1,15) 
405  format(15i4) 
  ENDDO
  call filclo(iun60,' ') 
  call filclo(iun61,' ') 
END subroutine iniham

! {\bf dercor} version 8; December 1996                                 
! dercor calculates partial derivatives of generating                   
! function and summs over planets accounted for                         
!                                                                       
!   input: el asteroid keplerian elements                               
!          pel perturbing planet keplerian elements                     
!          no ast. name (encoded) (only for error messages)             
!          tt time (idem)                                               
!          inp choice between tt,no                                     
!   output: iqc  quality code for mean elements                         
!           dd(6) corrections to Delaunay elements (to be multiplied by 
SUBROUTINE dercor(dd,iqc,pel,el,inp,nam0,tt,iun21)
! INPUT
  CHARACTER*9, INTENT(IN) ::  nam0 ! asteroid name
  INTEGER, INTENT(IN) :: iun21  ! resonance warning output unit
! =============================================================         
! el: 1=semimaj. ax. 2= ecc. 3=inclin.                                  
!        4=node 5=arg. peri 6=mean anomaly                              
! aul: 1=long. peri 2=mean long. 3=mean mot.                            
!  planetary elements; second index is 1=jup 2=sat                      
  DOUBLE PRECISION, INTENT(IN) :: pel(6,2)
!  asteroid elements (osculating)                                       
  DOUBLE PRECISION, INTENT(IN) :: el(6) 
! OUTPUT
!  corrections to Delaunay elements
  DOUBLE PRECISION, INTENT(OUT) ::  dd(6)
! quality control flag
  INTEGER, INTENT(OUT) :: iqc
! END INTERFACE
!  trigonometric functions to be used many times                        
  DOUBLE PRECISION ccg(ntx),scg(ntx),cck(350),sck(350) 
! planetary auxiliary variables, asteroid auxiliary
  DOUBLE PRECISION pau(5), aul(7) 
!  corrections to Delaunay elements: direct, indirect
  DOUBLE PRECISION sd(6),sdi(6) 
!  laplace coefficients                                                 
  DOUBLE PRECISION b(180),c(180),e(180) 
  DOUBLE PRECISION bm(180),cm(180),e0m(180) 
!  derivatives of laplace coefficients                                  
  DOUBLE PRECISION b1(180),b2(180),b3(180),b4(180),b5(180),                &
     &      c1(180),c2(180),c3(180),e1(180)                             
  DOUBLE PRECISION b1m(180),b2m(180),b3m(180),b4m(180),                    &
     &      c1m(180),c2m(180)                                           
!  leverrier coefficients, hamiltonian integer lists                    
  DOUBLE PRECISION cl(380),cla(380) 
! auxiliary variable for mass transformation, scalar temporaries 
  DOUBLE PRECISION alph(4), pa, pe,pin,po, a ,pai, tt,pery, dlcon,sc
  DOUBLE PRECISION cfaci,cargi, dgc0, tmp1, efac, cfac,sifac,si2fac,div,faca,cc,afac, carg, cargk
! functions
  DOUBLE PRECISION princ
  INTEGER nalfa
!  zero eccentricity/inclination flag                                   
  logical ezfl 
! loop indexes 
  INTEGER i,k,i1,i2, iiqc,m, inp, ipl, iplan, np, np1,nnp, iab
!  quality control                                                      
  iqc=0 
!                                                                       
!  Poincare' trick for masses depending on only one small parameter     
  DO ipl=1,2 
     alph(ipl)=gm(ipl+1)/gm(2) 
  ENDDO
!   asteroid mean motion rad/day.                                       
  aul(3)=sqrt(gm(1))/sqrt(el(1)**3) 
!  osculating longitude of perihelion, mean longitude,                  
!  also sin/cos of inclination                                          
  aul(1)=princ(el(5)+el(4)) 
  aul(2)=princ(el(6)+el(5)+el(4)) 
  aul(4)=sin(el(3)) 
  aul(5)=sin(el(3)/2.d0) 
  aul(6)=cos(el(3)) 
  aul(7)=cos(el(3)/2.d0)/2.d0 
!     initialization of corrections to l,g,h,w,o.                       
  dd=0.d0 
! *********************************************************             
!   loop on number of perturbing planets                                
  DO 1 iplan=1,2 
!     assignment of jupiter/saturn values to programme variables.       
     pa=pel(1,iplan) 
     pe=pel(2,iplan) 
     pin=pel(3,iplan) 
     po=pel(4,iplan) 
     pau(1)=princ(pel(4,iplan)+pel(5,iplan)) 
     pau(2)=princ(pel(6,iplan)+pau(1)) 
     pau(4)=sin(pin) 
     pau(5)=sin(pin/2.d0) 
!  planet mean motion rad/day                                           
     pau(3)=sqrt(smp(iplan)/pa**3) 
! **********************************************************            
!  interface for Laplace and Leverrier coefficients                     
!  inner/outer perturbation cases handled by Williams trick             
     if(el(1).lt.pa) then 
        a=el(1)/pa 
        pai=1.d0/pa 
     else 
        a=pa/el(1) 
        pai=1.d0/el(1) 
     endif
!  setting delimiters for calculation of Lap. Lev. coeff.               
     np=nalfa(a) 
!  two more Laplace coefficients needed to compute Lev. coeff.          
     np1=np+2 
     nnp=2*np-1 
!  calculation of the laplace coefficients.                             
     CALL lapco (a,np1,excrit,b,b1,b2,b3,b4,b5,c,c1,c2,c3,e,e1) 
     DO iab=1,np1 
        bm(iab)=b1(iab)/el(1) 
        b1m(iab)=(b2(iab)+b1(iab))/el(1) 
        b2m(iab)=(b3(iab)+2*b2(iab))/el(1) 
        b3m(iab)=(b4(iab)+3*b3(iab))/el(1) 
        b4m(iab)=(b5(iab)+4*b4(iab))/el(1) 
        cm(iab)=c1(iab)/el(1) 
        c1m(iab)=(c2(iab)+c1(iab))/el(1) 
        c2m(iab)=(c3(iab)+2*c2(iab))/el(1) 
        e0m(iab)=e1(iab)/el(1) 
     ENDDO
!  *********************************************************            
!  initialization of the correction to l,g,h,w,o                        
!  due to a single planet                                               
     sd=0.d0 
     sdi=0.d0 
! ***************************                                           
!   precomputation of the trigonometric functions (the part             
!   with arguments not depending upon k)                                
     DO m=1,ntd 
        carg = princ(ih(m,12)*pau(2)-ih(m,13)*aul(2)+ih(m,14)       &
    &      *pau(1)+ih(m,15)*aul(1)+ih(m,16)*po+ih(m,17)*el(4))         
        ccg(m)=cos(carg) 
        scg(m)=sin(carg)
     ENDDO
!   precomputation of trigonometric function, part with                 
!   arguments depending on k (but not upon m)                           
     DO k=1,nnp 
        i=k-np 
        cargk=princ(i*(pau(2)-aul(2))) 
        cck(k)=cos(cargk) 
        sck(k)=sin(cargk) 
     ENDDO
!*********************************************                          
!  main iteration on k                                                  
     DO 3 k=1,nnp 
        i=k-np 
        iab=abs(i) 
!  calculation of the leverrier coefficients. #i# < 21,51,101,173.      
        call levco(i,b,b1,b2,b3,b4,c,c1,c2,e,a,pai,cl) 
        call levco(i,bm,b1m,b2m,b3m,b4m,cm,c1m,c2m,e0m,a,pai,cla) 
! *********************************************************             
! computation of perturbations (main part) by means of the              
! generic terms.                                                        
        DO 7 m=1,ntd 
           if((i.eq.-ih(m,12)).and.(ih(m,12).eq.ih(m,13))) CYCLE 
           ezfl=.false. 
           if(el(2).lt.epsil)then 
              if(ih(m,6).eq.0)then 
                 efac=1.d0 
              else 
                 efac=0.d0 
                 ezfl=.true. 
              endif
           else 
              efac=el(2)**ih(m,6) 
           endif
           if(aul(4).lt.epsil)then 
              if(ih(m,8).eq.0)then 
                 sifac=1.d0 
              else 
                 sifac=0.d0 
                 ezfl=.true. 
              endif
              if(ih(m,10).eq.0)then 
                 si2fac=1.d0 
              else 
                 si2fac=0.d0 
                 ezfl=.true. 
              endif
           else 
              sifac=aul(4)**ih(m,8) 
              si2fac=aul(5)**ih(m,10) 
           endif
           div =(i+ih(m,12))*pau(3)-(i+ih(m,13))*aul(3) 
           faca = 1.d0/((i+ih(m,12))*pau(3)-(i+ih(m,13))*aul(3))           &
     &      *ih(m,1)/ih(m,2)*(-1)**ih(m,4)*efac                         &
     &      *pe**ih(m,7)*sifac*pau(4)**ih(m,9)                          &
     &      *si2fac*pau(5)**ih(m,11)                                    
           cfac = faca*cl(ih(m,3)) 
           afac = faca*cla(ih(m,3)) 
           if(ih(m,5).ne.0)then 
              cfac=cfac*i**ih(m,5) 
              afac=afac*i**ih(m,5) 
           endif
           cc=ccg(m)*cck(k)-scg(m)*sck(k) 
           sc=scg(m)*cck(k)+ccg(m)*sck(k) 
           sd(1)=sd(1)+(-i-ih(m,13))*cc*cfac 
           sd(2)=sd(2)+(-i-ih(m,13)+ih(m,15))*cc*cfac 
           sd(3)=sd(3)+(-i-ih(m,13)+ih(m,15)+ih(m,17))*cc*cfac 
           sd(4)=sd(4)+ih(m,6)/el(2)*sc*cfac 
           sd(5)=sd(5)+(ih(m,8)*aul(6)/aul(4)+ih(m,10)*                    &
     &           aul(7)/aul(5))*sc*cfac                                          
           sd(6)=sd(6)+sc*afac 
!  resonance warning                                                    
           if(abs(cfac).gt.eps.and.abs(div).lt.epdiv)then 
              i1=i+ih(m,12) 
              i2=i+ih(m,13) 
              pery=dpig/div 
              if(inp.eq.1)then 
                 write(iun21,110)nam0,i2,i1,ih(m,15),ih(m,17),ih(m,14)     &
     &              ,ih(m,16),cfac,pery                                 
110              format(a9,2i5,2x,4i3,2x,f10.3,2x,f7.1) 
              elseif(inp.eq.2)then 
                 write(iun21,111)tt,i2,i1,ih(m,15),ih(m,17),ih(m,14)       &
     &              ,ih(m,16),cfac,pery                                 
111              format(f12.2,2i5,2x,4i3,2x,f10.3,2x,f7.1) 
              endif
           endif
!  quality code                                                         
           dlcon=dlog(1.d1) 
           if(cfac.gt.0.d0)then 
              iiqc=10 +10*dlog(cfac)/dlcon 
           elseif(cfac.lt.0.d0)then 
              iiqc=10 +10*dlog(-cfac)/dlcon 
           elseif(.not.ezfl)then 
!  coefficient is zero, while it should not be                          
              if(i.ne.0.or.ih(m,5).eq.0)then 
                 write(*,*)el 
                 i1=i+ih(m,12) 
                 i2=i+ih(m,13) 
                 write(*,*)nam0,i,ih(m,12),ih(m,13),ih(m,15),ih(m,17),ih(m,14)&
     &              ,ih(m,16),faca,cl(ih(m,3)),ih(m,3)                  
                 goto 3 
              else 
                 iiqc=0 
              endif
           endif
           iqc=max(iqc,iiqc) 
7       ENDDO
3    ENDDO
! computation of perturbations (indirect part) by means of the          
! generic terms.                                                        
     DO m=1,nti 
        cfaci = 1.d0/(ihi(m,10)*pau(3)-ihi(m,11)*aul(3))                &
     &      *el(1)/pa/pa*ihi(m,1)/ihi(m,2)*(-1)**ihi(m,3)*el(2)**       &
     &      ihi(m,4)*pe**ihi(m,5)*aul(4)**ihi(m,6)*pau(4)**             &
     &      ihi(m,7)*aul(5)**ihi(m,8)*pau(5)**ihi(m,9)                  
        cargi = princ(ihi(m,10)*pau(2)-ihi(m,11)*aul(2)+ihi(m,12)*pau(1)  &
     &      +ihi(m,13)*aul(1)+ihi(m,14)*po+ihi(m,15)*el(4))             
        sdi(1)=sdi(1)-ihi(m,11)*cfaci*cos(cargi) 
        sdi(2)=sdi(2)+(-ihi(m,11)+ihi(m,13))*cfaci*cos(cargi) 
        sdi(3)=sdi(3)+(-ihi(m,11)+ihi(m,13)+ihi(m,15))*cfaci*cos(cargi) 
        sdi(4)=sdi(4)+ihi(m,4)/el(2)*cfaci*sin(cargi) 
        sdi(5)=sdi(5)+(ihi(m,6)*aul(6)/aul(4)+ihi(m,8)*                 &
     &      aul(7)/aul(5))*cfaci*sin(cargi)                            
        sdi(6)=sdi(6)+cfaci/el(1)*sin(cargi) 
     ENDDO
! *********************************************************             
!     summing up the corrections due to single planets.                 
     DO i=1,6 
        dd(i)=dd(i)+(sd(i)+sdi(i))*alph(iplan) 
     ENDDO
!  end of loop on planets                                               
1 ENDDO
! *********************************************************             
! Delaunay equations have minus sign for action variables               
  dd(1:3)=-dd(1:3)
! Delaunay equations are different from Lagrange equations              
  dgc0=sqrt(gm(1)*el(1)*(1-el(2)*el(2))) 
  tmp1=-dd(5)/(sin(el(3))*dgc0) 
  dd(6)=2.d0*sqrt(el(1))/sqrt(gm(1))*dd(6)+(1.d0-el(2)*el(2))/              &
     &       (sqrt(gm(1))*sqrt(el(1))*el(2))*dd(4) 
  dd(5)=-sqrt(1.d0-el(2)*el(2))/                                   &
     & (sqrt(gm(1))*sqrt(el(1))*el(2))*dd(4)-cos(el(3))*tmp1
  dd(4)=tmp1 
! *********************************************************             
END SUBROUTINE dercor


END MODULE short9

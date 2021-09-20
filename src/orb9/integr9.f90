! ===========================================                           
!   euler method of order 1                                             
! -------------------------------------------                           
subroutine euler(el,pel,nam0,tt,inp,h,em,iqc,iun21) 
  USE fund_const
  IMPLICIT NONE
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
  USE fund_const
  IMPLICIT NONE
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
  USE fund_const
  IMPLICIT NONE
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
  USE fund_const
  IMPLICIT NONE
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
  USE fund_const
  IMPLICIT NONE
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
  USE fund_const
  IMPLICIT NONE                      
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

! *********************************************************             
! mainpro vers. 9                                                       
! May 2000                                                              
! this program computes at once mean and proper elements                
! version with OEF (OrbFit compatible) input                            
! *********************************************************             
PROGRAM prop9
  USE massmod 
  USE short9
!  masses
  IMPLICIT NONE 
! max no. of entries                                                    
  INTEGER, PARAMETER :: numx=20000 
! =============================================================         
! el,em: 1=semimaj. ax. 2= ecc. 3=inclin.                               
!        4=node 5=arg. peri 6=mean anomaly                              
!  planetary elements; second index is 1=jup 2=sat                      
  DOUBLE PRECISION pel(6,2) 
!  asteroid keplerian elements (osculating, mean)                       
  DOUBLE PRECISION el(6),em(6) 
!  end of file                                                          
  logical eof 
!  file names                                                           
! these are for input planets and baric                                 
  character*72 filpla,filbar 
! real input/output file names                                          
  character*72 filast,filfil,filpro 
  INTEGER len 
! options                                                               
  INTEGER inp,nmax 
! input to short8s subroutine                                           
  DOUBLE PRECISION tt 
  CHARACTER*9 nam0 
! reference time for long9                                               
  DOUBLE PRECISION t0 
! outputs from short8 subroutine                                        
  INTEGER iqcm,it,ierr 
  DOUBLE PRECISION delta 
! osculating elements quality code, magnitude just copied               
  INTEGER iqco 
  DOUBLE PRECISION hm 
! counters                                                              
  INTEGER kr,krdon 
! output proper elements:frequencies, actions, angles                   
  DOUBLE PRECISION gp,sp,g0,s0,pre,prsini,apip,atep,apiv,atev 
  INTEGER ngper,ngnod 
! quality codes                                                         
  INTEGER iqce,iqcf,irfl,irfl2,ipas 
! unit conversion                                                       
  DOUBLE PRECISION conv 
! unit numbers                                                          
  INTEGER iun1,iun2,iun3,iun8,iun9,iun10,iun11,iun19,iun21,iun18    &
     &     ,iun7,iun17, iundone                                                 
! ********************************************************              
! use orbfit library                                                    
  CALL libini 
!  choice of input format, options:                                     
  CALL filopn(iun1,'prop.opt','old') 
  call skip (iun1,1) 
!  inp=1 from catalogue; inp=2 numerical integration; inp=3 filtered dat
  call reaint(iun1,'inp',inp) 
  if(inp.le.0.or.inp.gt.3)then 
     write(*,*)' unknown input option, inp= ',inp 
     stop 
  endif
!  control on number of asteroids                                       
  call reaint(iun1,'nmax',nmax) 
! ----------------------------------------------------------------------
! input/output file names                                               
!  input file: asteroid osculating elements                             
  call reastr(iun1,'filast',filast) 
! Warning: INP=3 IS CURRENTLY NOT SUPPORTED
!  input file for mean (filtered) elements:                             
! (only for inp=3, otherwise just to skip line in option file
  call reastr(iun1,'filfil',filfil) 
  IF(INP.eq.3)THEN
     WRITE(*,*)' inp=3 NOT SUPPORTED'
     STOP
  ENDIF
! End Warning AM, ZK 14/6/2017
!  output file: beginning of name                                       
  call reastr(iun1,'filpro',filpro) 
! ----------------------------------------------------------------------
! planetary data for short periodic perturbations and barycentric correc
  call skip(iun1,1) 
!  input file: planets                                                  
  call reastr(iun1,'filpla',filpla) 
!  input file: barycenter of the inner solar system; used only for inp=3
!  ibar=1 barycentric correction, ibar=0 no correction                  
  call reastr(iun1,'filbar',filbar) ! only to skip line in .opt file
!  end input options                                                    
  close(iun1) 
! ********************************************************************* 
!  initialisations, file opening, for short and long                    
! **********************************************************************
  call openfi(inp,t0,filast,filfil,filpro,filpla,filbar,            &
     &    iun2,iun3,iun8,iun9,iun10,iun11,iun19,iun21,              &
! iun18,
     &    iun7,iun17)  
! *********************************************************             
!   loop on input asteroid data begins here                             
!  counter for mean/proper elements                                     
  krdon=0 
  DO 10 kr=1,nmax ! WARNING: update nmax in prop.opt.allnum, now 600,000
!   input osculating elements                                           
     if(inp.eq.1.or.inp.eq.2)then 
        call iosho(kr,inp,iun8,iun7,iun17,nam0,t0,tt,el,            &
   &           pel,eof,iqco,hm,iun19)                                 
        if(kr.eq.1000*(kr/1000))then 
           if(inp.eq.1)then 
              write(*,*)nam0,'  read no ',kr 
           elseif(inp.eq.2)then 
              write(*,*)tt,'  read no ',kr 
           endif
        endif
        if(eof)goto 99 
! *********************************************************             
!  excluding the non-main-belt asteroids.                               
        if(inp.eq.1)then 
           if ((el(1)*(1.d0-el(2))).lt.1.1d0.or.el(1).gt.3.8d0      &
     &              .or.el(1).lt.1.8d0)then                             
              goto 10 
           endif
        endif
! counter for mean elements                                             
        krdon=krdon+1 
! ***************************************************                   
        CALL short(el,pel,nam0,tt,inp,em,iqcm,it,delta,ierr,iun19,iun21)
! ********************************************************              
!  output of the mean elements                                          
        if(inp.eq.1)then 
           write(iun21,911)nam0,it,delta 
911        format(a9,i3,d12.4) 
           write (iun19,410) nam0,em(6),em(5),em(4),em(3),em(2),    &
     &             em(1),iqcm,ierr,iqco,hm                              
410        format(a9,4f10.5,f11.7,f11.7,i3,i3,i3,1x,f5.2) 
        elseif(inp.eq.2)then 
           write(iun21,912)tt,it,delta 
912        format(f12.2,i3,d12.4) 
           write (iun19,411)tt,em(6),em(5),em(4),em(3),em(2),em(1), &
     &              iqcm,ierr                                           
411        format(f12.2,4f10.5,f10.7,f12.7,i3,i3) 
        endif
! ********************************************************              
     else 
!  filtered mean elements to be read for inp=3                          
!   input from output of orbit8v, filtered                              
!   only one orbit in the input, equinoctal
        WRITE(*,*)' inp=3 not supported, long period'
        STOP
        call inpfil(iun18,tt,em,eof) 
!   iqcm is not given, it is assumed that filtering has                 
!   eliminated all the short periods anyway                             
        iqcm=0 
!  end job and counters                                                 
        if(eof)goto 99 
        if(kr.eq.(kr/1000)*1000)then 
           write(*,*)tt,'  read  ',kr 
        endif
! ********************************************************              
!  mean elements are available                                          
     endif
! compute proper elements                                               
     call long9(em,nam0,tt,t0,inp,                                  &
      &        gp,sp,pre,prsini,iqce,                                    &
      &        apip,atep,g0,s0,ipas,iqcf,                                &
      &        irfl,irfl2,                                               &
      &        ngper,apiv,ngnod,atev,conv,                               &
      &        iun2,iun3,iun9,iun10,iun11)                               
! output proper elements                                                
     call selout(inp,tt,nam0,                                       &
      &        gp,sp,em(1),pre,prsini,iqce,iqcm,iqco,                    &
      &        hm,apip,atep,g0,s0,ipas,iqcf,                             &
      &        irfl,irfl2,                                               &
      &        ngper,apiv,ngnod,atev,conv,iun2,iun3,iun9,iun11)          
! counter                                                               
     if(krdon.eq.1000*(krdon/1000))then 
        if(inp.eq.1)then 
           write(*,*)nam0,'  mean el. done no ',krdon 
        elseif(inp.eq.2.or.inp.eq.3)then 
           write(*,*)tt,'  mean el. done no ',krdon 
        endif
     endif
! ***************************************************                   
!  end loop on asteroid catalog                                         
10 ENDDO
!     displaying the information on job termination.                    
99 continue 
  write(*,*)krdon,'  mean and proper elements generated ' 
! proper elements output files (always)                                 
  call filclo(iun2,' ') 
  call filclo(iun9,' ') 
  call filclo(iun10,' ') 
! discarded proper elements (only for catalog input)                    
  IF(inp.eq.1)call filclo(iun3,' ') 
  IF(inp.eq.1)call filclo(iun11,' ') 
! mean elements output files (if mean elements are computed)            
  call filclo(iun19,' ') 
  call filclo(iun21,' ') 
! input planetary data (if mean elements are computed)                  
  call filclo(iun17,' ') 
  call filclo(iun7,' ') 
! input filtered data (if mean elements are not computed)               
!  IF(inp.eq.3)call filclo(iun18,' ') 
  CALL filopn(iundone,'prop9.done','unknown') 
  CALL filclo(iundone,' ')
END PROGRAM prop9

! ********************************************************************* 
!   version 8 (1997):                                                   
!   improvements with respect to previous versions:                     
!   z2 second order terms added (6.8)                                   
!   nonlinear forced terms reduced to 2g6-g5 (6.7.7)                    
!   selection of secular perturbation theory (6.7)                      
!   with correction to frequencies before iteration (6.5.3)             
!   with second order degree four terms 3:1 (6.5.1)                     
!   with nonlinear forced terms for eccentricity (6.1)                  
!   s7 terms in degree four long periodic terms (5.7)                   
!   generic terms (6.3.0)                                               
!   input/output for 1994 format elements with magnitudes, QCO          
!   improved logic for divergent cases, with automatic switch to adapted
!   theory, using two consecutive iterations (6.5.2)                    
! ********************************************************************* 
SUBROUTINE long9(em,nam0,tt,t0,inp,                               &
     &        gp,sp,pre,prsini,iqce,                                    &
     &        apip,atep,g0,s0,npas,iqcf,                                &
     &        irfl,irfl2,                                               &
     &        ngper,apiv,ngnod,atev,conv,                               &
     &        iun2,iun3,iun9,iun10,iun11)                               
  USE fund_const  
  USE controlmod  
  USE terms9
  IMPLICIT NONE 
! INPUT                                                                 
! elements, epoch, reference epoch                                       
  DOUBLE PRECISION, INTENT(IN) :: em(6),tt,t0 
! asteroid name, input type 
  CHARACTER*9, INTENT(IN) :: nam0 
  INTEGER, INTENT(IN) :: inp 
! unit numbers                                                          
  INTEGER, INTENT(IN) :: iun2,iun3,iun9,iun10,iun11 
! OUTPUT                                                                
! frequencies (nonlinear, linear)                                       
  DOUBLE PRECISION, INTENT(out) :: gp,sp,g0,s0 
! proper elements                                                       
  DOUBLE PRECISION, INTENT(out) :: pre,prsini,apip,atep 
! proper angles with number revolutions                                 
  DOUBLE PRECISION, INTENT(out) :: apiv,atev 
  INTEGER, INTENT(out) :: ngper,ngnod 
! quality flags                                                     
  INTEGER, INTENT(out) :: iqce,iqcf,irfl,irfl2 
! END INTERFACE     
! ====================== proper elems. per iteration====================
!  vector of Poincare' elements, its derivative with respect to mass,   
!  linear proper elements                                               
  DOUBLE PRECISION rho(4),rho0(4) 
!  frequencies (internal units) and distances (in arcsec/year)          
  DOUBLE PRECISION gv(maxit),sv(maxit) 
!  differences of frequencies (arcsec/yr), of proper amplitude          
  DOUBLE PRECISION ds(maxit),dg(maxit),d(maxit),da(maxit) 
!  previous values of  mu, nu (proper amplitudes for the asteroid)      
!  and proper elements (including angles)                               
  DOUBLE PRECISION amiv(maxit),aniv(maxit) 
  DOUBLE PRECISION api(maxit),ate(maxit),prev(maxit),prsinv(maxit) 
!  divisors, amplitudes of largest terms                                
  DOUBLE PRECISION dev(ndex,maxit),amp(ndex,maxit),                 &
     &       ampi(ndex,maxit),ampe(ndex,maxit)                          
! ======================                                                
!  forced terms for the asteroid xi_j, eta_j; j frequency index         
  DOUBLE PRECISION xi(8),et(8) 
!  nonlinear forced terms: frequencies, phases, amplitudes, forced ecc. 
  DOUBLE PRECISION fnl(nnlx),phanl(nnlx),pninl(2,nnlx) 
!  cross reference index for nonlinear forced terms                     
  INTEGER inl(nnlx) 
! ======================                                                
!  secular theory constants for the planets                             
!  one-dim arrays; index is frequency index lambda_j, delta_j           
!  nu_j, mu_j; j is frequency index for the current planet              
  DOUBLE PRECISION pl(8),pd(8) 
!     trigonometric functions                                           
  DOUBLE PRECISION cld(-2:2,-2:2,-2:2,-2:2),sld(-2:2,-2:2,-2:2,-2:2) 
!  secular theory frequencies for the planets g_j, s_j                  
  DOUBLE PRECISION pg(8),ps(8) 
!  dummy variable for protection of constants not to be rescaled twice  
  DOUBLE PRECISION dum(8) 
!  two-dim arrays; first index is planet index, second is freq. index   
! nu_{j,k}, mu_{j,k} with j planet ind, k freq. index                   
  DOUBLE PRECISION pmi(8,8),pni(8,8) 
!  masses, mean motions and semimajor axes of the planets               
  DOUBLE PRECISION pmv(8),pn(8),pa(8) 
! ======================                                                
!  yuasa coefficients  a and b (first order effects)                    
!  warning: here planet index is the second one, first is Yuasa subscrip
  DOUBLE PRECISION yya(4,8),yyb(28,8) 
! ======================                                                
! mean elements to be used for discarded cases                          
!  DOUBLE PRECISION ae,ai0,aw0,ao0 
! constants and units                                                   
  DOUBLE PRECISION sm,conv,aunit,tunit,pmunit 
! Poincare' variables                                                   
  DOUBLE PRECISION ap0,aq0,ar0,as0,ap,aq,ar,as,al,aa,a,an 
! frequencies                                                           
  DOUBLE PRECISION g,s 
! proper elements                                                       
  DOUBLE PRECISION ani,ami,api1,ate1 
! other variables                                                       
  DOUBLE PRECISION awtrad,ao0rad,sini0,dej,da1,dani,dami 
  INTEGER ipas,nnl,irescl,kk,ive1,nde,ires,ip1,npas
  INTEGER iun1 ! options file unit
! =================================================================**** 
!  asteroid name alphanumeric                                           
  CHARACTER*6 astna 
! options from long8s.opt file                                          
  DOUBLE PRECISION eps,deps 
  INTEGER ive 
! loop indexes                                                          
  INTEGER i,j 
! counters                                                              
  INTEGER ndisc,ndiver,nada,nada2 
! flags                                                                 
  INTEGER idisc 
! initialisation                                                        
  INTEGER flag 
  SAVE 
  DATA flag /0/ 
  IF(flag.eq.0)THEN 
! ********************************************************************* 
!  read controls from file long8s.opt                                   
! ********************************************************************* 
     call filopn(iun1,'long8s.opt','old') 
     call skip(iun1,1) 
!  control for amplitude of term to be considered resonant              
     call reaflo(iun1,'eps',eps) 
!  control for divisor rated small                                      
     call reaflo(iun1,'deps',deps) 
!  version flag                                                         
     call reaint(iun1,'ive',ive) 
!  version-dependent options, defining the theory used (see optxx)      
     call lonopt(inp,ive) 
     call filclo(iun1,' ') 
! complete headers of output files                                      
     IF(inp.eq.1)THEN 
        write(iun2,503)vers 
        write(iun3,503)vers 
503     format(' version ',a30/                                     &
     &           '   no; a (au); ecc ; sinI;  g ("/yr); s ("/yr);'      &
     &           ,'RFL;QCM;QCO;magn.')                                  
        write(iun9,504)vers 
        write(iun11,504)vers 
504     format(' version ',a30/'   no; l.per, l.node(deg); g0 ("/yr) ;'&
     &           ,'s0 ("/yr);npas;QCF;QCE;IRE2')                        
        write(iun10,506)vers,eps,deps 
506     format(' version ',a30/' only amplitudes > ',1p,d10.2,      &
     &           '  with divisors <',d10.2)                             
        write(iun10,505) 
505     format('   no; iter; no.div; divisor ("/yr);'               &
     &           ,' amplitude e; amplitude I')                          
     ELSEIF(inp.eq.2.or.inp.eq.3)THEN 
        write(iun2,513)vers 
  513       format(' version ',a30/                                     &
     &         ' time; a (au); ecc ; sinI;  g ("/yr); s ("/yr);QCE;QCM')
        write(iun9,514)vers 
  514       format(' version ',a30/                                     &
     &           ' time; l.per, l.node(deg); g0 ("/yr) ;'               &
     &           ,' s0 ("/yr);npas;RFL;QCE;IRE2')                       
        write(iun10,506)vers,eps,deps 
        write(iun10,515) 
  515       format(                                                     &
     &           ' time; iter; no.div; total amplitude; divisor ("/yr)')
     ENDIF
!   counters are used to keep track of the number of revolutions of     
!   both api1 and ate1 (proper perihelion and node)                     
     ngper=0 
     apiv=0.d0 
     ngnod=0 
     atev=0.d0 
!  read integer lists (generic term representation) of the              
!  determining function                                                 
! *********************************************                         
! open file with generic integer records                                
     open(50,file='gener.lng') 
     call rgener(50) 
     close(50) 
!  counters for discarded, non-convergent and adapted cases             
     ndisc=0 
     ndiver=0 
     nada=0 
     nada2=0 
! ********************************************************************  
!  if the time is fixed,                                                
!  the secular perturbation theory has to be consulted to find the      
!  appropriate semimajor axis of the planets, and later the secular     
!  eccentricities and inclinations will be used as well as the          
!  fundamental secular angles (see secthe5)                             
!  the fixed time is defined by t0 (in JD), with 0 shift                
! ********************************************************************* 
     call secthe(0.d0,t0,pmv,pa,pni,pg,pl,pmi,ps,pd     &
     &        ,nnl,fnl,pninl,phanl)                                
!  trigonometric functions of the fundamental angles (see trigfu)       
     call trigfu(pl(5),pl(6),pd(7),pd(6),sld,cld) 
!  data on the inner planets for computation of g0, s0                  
     call inplda(pa,pmv) ! these are constant data, even for inp=2
! ********************************************************************* 
!  constants and units                                                  
! *******************************************************************   
! change of scale, units, fundamental constants                         
     call convun(pmv,pa,pg,ps,pn,aunit,tunit,pmunit,nplin0,nplou,fnl,nnl)
! square root of solar mass                                             
     sm=dsqrt(1.d0/pmunit) 
!  conversion factor to get frequencies in arcsec/year                  
     conv=1.d0/tunit*360.d0/dpig*3600.d0 
! ********************************************************************* 
! end initialisation                                                    
     flag=1 
  ENDIF
! ======================================================================
! actual computations to be done at each call                           
! ======================================================================
! conversion from keplerian to Poincare' variables                      
  call keppoinc(em,ap,aq,ar,as,al,a,aa,aunit,sm) 
! ********************************************************************* 
! the secular perturbation theory has to be consulted every time if the 
! time is changing; rescaling is not necessary since only the phases    
! are recomputed                                                        
  if(inp.gt.1)then 
     call secthe(tt,t0,pmv,pa,pni,pg,pl,pmi,ps,pd    &
     &       ,nnl,fnl,pninl,phanl)                                 
     call trigfu(pl(5),pl(6),pd(7),pd(6),sld,cld) 
!  data on the inner planets for computation of g0, s0                  
     call inplda(pa,pmv) 
! ********************************************************************* 
!  constants and units                                                  
! *******************************************************************   
! change of scale, units, fundamental constants                         
     call convun(pmv,pa,pg,ps,pn,aunit,tunit,pmunit,nplin0,nplou,fnl,nnl)
  endif
! ********************************************************************* 
!  functions of the asteroid semimajor axis                             
! ********************************************************************* 
!  asteroid mean motion                                                 
  an=sm/al**3 
! ********************************************************************* 
!   laplace's, leverrier's and yuasa's coefficients.                    
! ********************************************************************* 
  do j=nplin0,nplou 
     call llyco(a,pa(j),an,pn(j),sm,yya(1,j),yyb(1,j),pmv(j),ipla(j))
  enddo
! ********************************************************************* 
!   second order/degree frequencies and forced oscillation amplitudes.  
! ********************************************************************* 
  call g0s0(al,yya,pmv,g0,s0,nplin0,nplou0) 
!   forced oscillation amplitudes; from outer planets only              
  call forcedgs(g0,s0,pg,ps,pmv,pni,pmi,yya,al,xi,et,nplin,nplou) 
! ********************************************************************* 
!                                                                       
!  introducing the running variables:                                   
!  only the proper modes are involved in the iteration                  
!                                                                       
! ********************************************************************* 
!  linear proper modes and linear proper actions                        
  call propmo(ap,aq,ar,as,al,xi,et,pl,pd,nplin,nplou,rho0)
! ********************************************************************* 
!   ready for the first iteration: initial frequencies are the          
!   ones obtained by theory of order 2, degree 2; initial               
!   elements are the linear proper ones                                 
! ********************************************************************* 
  rho(1:4)=rho0(1:4) 
  gv(1)=g0 
  sv(1)=s0 
! ********************************************************              
!  check for the possible blow up of the proper elements since          
!  the linear theory (this can occur e.g. near g0=g6); otherwise        
!  continue with the nonlinear theory                                   
! ********************************************************              
  call elems8(rho,ani,ami,api1,ate1,pre,prsini,idisc) 
  amiv(1)=ami 
  aniv(1)=ani 
!  test flag idisc                                                      
  if(idisc.eq.0)then 
!  OK, store proper elements as first iteration                         
     prev(1)=pre 
     api(1)=api1 
     prsinv(1)=prsini 
     ate(1)=ate1 
  else 
!  ni>2 and/or mi>2:in the case discarded since linear theory,          
!  the mean values and g0,s0 are output                                 
     sini0=sin(em(3)*dpig/3.6d2) 
     awtrad=(em(4)+em(5))*dpig/3.6d2 
     ao0rad=em(4)*dpig/3.6d2 
!  with  quality codes 200 (idisc=1,3) or 201 (idisc=2)                 
     if(idisc.eq.2)then 
        iqcf=201 
        iqce=201 
     else 
        iqcf=200 
        iqce=200 
     endif
     ndisc=ndisc+1 
     irfl=iqce 
     gp=g0 
     sp=s0 
     pre=em(2) 
     prsini=sini0 
     apip=awtrad 
     atep=ao0rad 
     npas=0 
!         call selout(inp,namfor,t,no,                                  
!    +        g0,s0,aa,ae,sini0,iqce,iqcm,iqco,                         
!    +        hm,awtrad,ao0rad,g0,s0,ipas,iqcf,                         
!    +        irfl,irfl2,                                               
!    +        ngper,apiv,ngnod,atev,conv)                               
     return 
  endif
! ********************************************************************* 
!  no divergence detected so far                                        
  irescl=0 
!  no removal so far                                                    
  irfl=0 
  irfl2=0 
! ********************************************************************* 
2 continue 
  ipas=0 
! ********************************************************************* 
!                                                                       
!  computation of a starting point (other than linear proper)           
!                                                                       
! ********************************************************************* 
  if(ive.eq.3)then 
! ===================================================================== 
!   breiter method:                                                     
!                                                                       
!   first call euler, then rest of RK; handled (from point of view      
!   of divergence, removal, etc.) as two separate iterations            
! ===================================================================== 
     do 9 kk=1,2 
        if(kk.eq.1)then 
           ive1=-1 
        else 
           ive1=ive 
        endif
        ipas=ipas+1 
! ===================================================================== 
!   first line: functions of a only                                     
!   second line: functions of planetary orbits                          
!   third  line: controls for resonance removal                         
!   fourth line: linear and current elements,                           
!   fifth line: output: freq., divisors, amplitudes                     
!   sixth line:          new elements                                   
        call inte8(ive1,yya,yyb,xi,et,al,aa,g0,s0,                   &
     &        pni,pmi,pg,ps,pmv,pninl,nnl,inl,cld,sld,phanl,            &
     &        irfl,irfl2,                                               &
     &        rho0,rho,                                                 &
     &        g,s,dev(1,ipas),nde,amp(1,ipas),ampi(1,ipas),ampe(1,ipas),&
     &        rho)                                                      
! ********************************************************************* 
!    test for resonance                                                 
! ********************************************************************* 
        DO 41 j=1,nde 
           if(amp(j,ipas).gt.eps)then 
              dej=dev(j,ipas)*conv 
              if(abs(dej).lt.deps)then 
                 if(inp.eq.1)then 
                    write(iun10,109)nam0,ipas,j,dej,ampe(j,ipas),ampi(j,ipas) 
                 else 
                    write(iun10,119)tt,ipas,j,ampe(j,ipas),ampi(j,ipas),dej 
                 endif
              endif
           endif
41      ENDDO
! ********************************************************************* 
        gv(ipas+1)=g 
        sv(ipas+1)=s 
        call elems8(rho,ani,ami,api1,ate1,pre,prsini,idisc) 
        aniv(ipas+1)=ani 
        amiv(ipas+1)=ami 
!   test blow-up flag idisc                                             
        if(idisc.eq.0)then 
!   OK, store proper elements from this iteration                       
           prev(ipas+1)=pre 
           api(ipas+1)=api1 
           prsinv(ipas+1)=prsini 
           ate(ipas+1)=ate1 
        else 
!  in the discarded case, if ni>2 and/or mi>2                           
!  the output is the proper elements and frequencies                    
!  of the previous iteration,                                           
!  with  quality codes 150/160 (idisc=1,3) or 151/161 (idisc=2)         
           if(idisc.eq.2)then 
              iqcf=141+10*kk 
              iqce=141+10*kk 
           else 
              iqcf=140+10*kk 
              iqce=140+10*kk 
           endif
!  if divergent, try to eliminate the small divisor term                
!  and go back  to the first iteration                                  
           da1=1.d-1 
           call selres(amp,dev,da1,ipas,irfl,irfl2,conv,ires    &
     &               ,nde,nam0,tt,inp,nada,nada2,rho0,rho,iun10)        
           if(ires.gt.0)then 
              irescl=0 
              goto 2 
           else 
!  if removal impossible, give up with irfl>99                          
              ndisc=ndisc+1 
              irfl=iqce 
              gp=gv(ipas) 
              sp=sv(ipas) 
              pre=prev(ipas) 
              prsini=prsinv(ipas) 
              apip=api(ipas) 
              atep=ate(ipas) 
              npas=ipas 
!                call selout(inp,namfor,t,no,                           
!    +             gv(ipas),sv(ipas),aa,prev(ipas),prsinv(ipas),        
!    +             iqce,iqcm,iqco,hm,                                   
!    +             api(ipas),ate(ipas),g0,s0,ipas,iqcf,                 
!    +             irfl,irfl2,                                          
!    +             ngper,apiv,ngnod,atev,conv)                          
              return 
           endif
        endif
9    enddo
! ********************************************************************* 
  elseif(ive.gt.0)then 
     write(*,*)' option not yet supported, ive=',ive 
     stop 
  endif
! ********************************************************************* 
!                                                                       
!  iteration starts here                                                
!                                                                       
! ********************************************************************* 
1 continue 
!  no do loop counter, because I want to                                
!  be able to reset it back! however, this is the beginning of the      
!  iteration loop                                                       
  ipas=ipas+1 
! ===================================================================== 
!   euler method (ive<0), fixed point iteration                         
! ===================================================================== 
!   first line: functions of a only                                     
!   second line: functions of planetary orbits                          
!   third  line: controls for resonance removal                         
!   fourth line: linear and current elements,                           
!   fifth line: output: freq., divisors, amplitudes                     
!   sixth line:          new elements (here at the same address,        
!                        because of the fixed point iteration)          
  ive1=-1
  call inte8(ive1,yya,yyb,xi,et,al,aa,g0,s0,                        &
     &   pni,pmi,pg,ps,pmv,pninl,nnl,inl,cld,sld,phanl,                 &
     &   irfl,irfl2,                                                    &
     &   rho0,rho,                                                      &
     &   g,s,dev(1,ipas),nde,amp(1,ipas),ampi(1,ipas),ampe(1,ipas),     &
     &   rho)                                                           
! ********************************************************************* 
!    test for resonance                                                 
! ********************************************************************* 
  DO 40 j=1,nde 
     if(amp(j,ipas).gt.eps)then 
        dej=dev(j,ipas)*conv 
        if(abs(dej).lt.deps)then 
           if(inp.eq.1)then 
              write(iun10,109)nam0,ipas,j,dej,ampe(j,ipas),ampi(j,ipas) 
109           format(a9,2i4,f9.4,1p,2e12.4) 
           else 
              write(iun10,119)tt,ipas,j,ampe(j,ipas),ampi(j,ipas),dej 
119           format(f12.2,2i4,1p,e12.4,e12.4,0p,f9.4) 
           endif
        endif
     endif
40 ENDDO
! ===================================================================== 
  gv(ipas+1)=g 
  sv(ipas+1)=s 
! ********************************************************************* 
!   computation of proper elements nu, mu, ecc, sini and angles;        
!   test for cases to be discarded:                                     
!   this can occur near nonlinear resonances.                           
!   idisc=0 OK idisc=1 ni>2 idisc=2 mi>2 idisc=3 mi and ni >1           
! ********************************************************************* 
  call elems8(rho,ani,ami,api1,ate1,pre,prsini,idisc) 
  aniv(ipas+1)=ani 
  amiv(ipas+1)=ami 
!   test blow-up flag idisc                                             
  if(idisc.eq.0)then 
!   OK, store proper elements from this iteration                       
     prev(ipas+1)=pre 
     api(ipas+1)=api1 
     prsinv(ipas+1)=prsini 
     ate(ipas+1)=ate1 
  else 
!  in the discarded case, if ni>2 and/or mi>2                           
!  the output is the proper elements and frequencies                    
!  of the previous iteration,                                           
!  with  quality codes 100 (idisc=1,3) or 101 (idisc=2)                 
     if(idisc.eq.2)then 
        iqcf=101 
        iqce=101 
     else 
        iqcf=100 
        iqce=100 
     endif
!  if divergent, try to eliminate the small divisor term                
!  and go back  to the first iteration                                  
     da1=1.d-1 
     call selres(amp,dev,da1,ipas,irfl,irfl2,conv,ires         &
     &               ,nde,nam0,tt,inp,nada,nada2,rho0,rho,iun10)        
     if(ires.gt.0)then 
        irescl=0 
        goto 2 
     else 
!  if removal impossible, give up with irfl>99                          
        ndisc=ndisc+1 
        irfl=iqce 
        gp=gv(ipas) 
        sp=sv(ipas) 
        pre=prev(ipas) 
        prsini=prsinv(ipas) 
        apip=api(ipas) 
        atep=ate(ipas) 
        npas=ipas 
!            call selout(inp,namfor,t,no,                               
!    +        gv(ipas),sv(ipas),aa,prev(ipas),prsinv(ipas),             
!    +        iqce,iqcm,iqco,hm,                                        
!    +        api(ipas),ate(ipas),g0,s0,ipas,iqcf,                      
!    +        irfl,irfl2,                                               
!    +        ngper,apiv,ngnod,atev,conv)                               
        return 
     endif
  endif
! ********************************************************************* 
! is the next iteration necessary? compare with convergence control     
! for the frequencies g, s                                              
! ********************************************************************* 
  dg(ipas)=(gv(ipas+1)-gv(ipas))*conv 
  ds(ipas)=(sv(ipas+1)-sv(ipas))*conv 
  d(ipas)=dsqrt(dg(ipas)**2+ds(ipas)**2) 
  dani=dsqrt(aniv(ipas+1))-dsqrt(aniv(ipas)) 
  dami=dsqrt(amiv(ipas+1))-dsqrt(amiv(ipas)) 
  da(ipas)=dsqrt(dani**2+dami**2) 
!  quality codes based upon elements and frequencies convergence        
  iqce=40+10*dlog(da(ipas))/dlog(10.d0) 
  iqcf=30+10*dlog(d(ipas))/dlog(10.d0) 
!  convergence control                                                  
  if((dabs(dg(ipas)).lt.fconte).and.(dabs(ds(ipas)).lt.fconti))then
! ****************************************************************      
! if convergence occurred, write output.                                
! number of iteration is the current one. remark: the number of         
! iteration is 1 if the first degree 4 correction is adopted            
     ip1=ipas+1 
! check for divisor z3, not appearing in theory but sometimes           
! important in the main belt                                            
     if(irfl.eq.0)then 
        if(abs(dev(71,ipas))*conv.lt.0.5d0)irfl=71 
     elseif(irfl2.eq.0)then 
        if(abs(dev(71,ipas))*conv.lt.0.5d0)irfl2=71 
     endif
     gp=gv(ip1) 
     sp=sv(ip1) 
     pre=prev(ip1) 
     prsini=prsinv(ip1) 
     apip=api(ip1) 
     atep=ate(ip1) 
     npas=ipas 
!        call selout(inp,namfor,t,no,                                   
!    +        gv(ip1),sv(ip1),aa,pre,prsini,iqce,iqcm,iqco,             
!    +        hm,api(ip1),ate(ip1),g0,s0,ipas,iqcf,                     
!    +        irfl,irfl2,                                               
!    +        ngper,apiv,ngnod,atev,conv)                               
     return 
! *****************************************************************     
! if divergence is already showing up, stop before it is too late!      
  elseif(ipas.ge.4)then 
     if(da(ipas).gt.divrat*da(ipas-1))then 
        if(iqce.gt.iqcr)then 
!  if strongly divergent, try to eliminate the small divisor term       
!  and go back to the first iteration                                   
           da1=da(ipas) 
           call selres(amp,dev,da1,ipas,irfl,irfl2,conv,ires   &
    &              ,nde,nam0,tt,inp,nada,nada2,rho0,rho,iun10)         
           if(ires.gt.0)then 
!  successful removal, restart iteration                                
              irescl=0 
              goto 2 
           elseif(ires.lt.0)then 
!   already two terms removed, give up                                  
              ndiver=ndiver+1 
              npas=-ipas 
              if(irfl.eq.0)then 
                 if(abs(dev(71,ipas))*conv.lt.0.5d0)irfl=71 
              elseif(irfl2.eq.0)then 
                 if(abs(dev(71,ipas))*conv.lt.0.5d0)irfl2=71 
              endif
              gp=gv(ipas) 
              sp=sv(ipas) 
              pre=prev(ipas) 
              prsini=prsinv(ipas) 
              apip=api(ipas) 
              atep=ate(ipas) 
!           call selout(inp,namfor,t,no,                                
!    +             gv(ipas),sv(ipas),aa,prev(ipas),prsinv(ipas),        
!    +             iqce,iqcm,iqco,hm,                                   
!    +             api(ipas),ate(ipas),g0,s0,npas,iqcf,                 
!    +             irfl,irfl2,                                          
!    +             ngper,apiv,ngnod,atev,conv)                          
              return 
           elseif(ires.eq.0)then 
              irescl=irescl+1 
              if(irescl.ge.3)then 
!   unsuccessful removal after two divergences, give up
                 ndiver=ndiver+1 
                 npas=-ipas 
                 if(irfl.eq.0)then 
                    if(abs(dev(71,ipas))*conv.lt.0.5d0)irfl=71 
                 elseif(irfl2.eq.0)then 
                    if(abs(dev(71,ipas))*conv.lt.0.5d0)irfl2=71 
                 endif
                 gp=gv(ipas) 
                 sp=sv(ipas) 
                 pre=prev(ipas) 
                 prsini=prsinv(ipas) 
                 apip=api(ipas) 
                 atep=ate(ipas) 
!               call selout(inp,namfor,t,no,                            
!    +             gv(ipas),sv(ipas),aa,prev(ipas),prsinv(ipas),        
!    +             iqce,iqcm,iqco,hm,                                   
!    +             api(ipas),ate(ipas),g0,s0,npas,iqcf,                 
!    +             irfl,irfl2,                                          
!    +             ngper,apiv,ngnod,atev,conv)                          
                 return 
              endif
           endif
        else 
!   divergence not too strong, write as is                              
!  the elements of the previous iteration are printed                   
!  the current frequencies are not used (because they are "worse")      
!  however, the quality code is based on the discarded iteration(?)     
!  the number of iterations npas is set to a negative value             
           ndiver=ndiver+1 
           npas=-ipas 
           if(irfl.eq.0)then 
              if(abs(dev(71,ipas))*conv.lt.0.5d0)irfl=71 
           elseif(irfl2.eq.0)then 
              if(abs(dev(71,ipas))*conv.lt.0.5d0)irfl2=71 
           endif
           gp=gv(ipas) 
           sp=sv(ipas) 
           pre=prev(ipas) 
           prsini=prsinv(ipas) 
           apip=api(ipas) 
           atep=ate(ipas) 
!          call selout(inp,namfor,t,no,                                 
!    +             gv(ipas),sv(ipas),aa,prev(ipas),prsinv(ipas),        
!    +             iqce,iqcm,iqco,hm,                                   
!    +             api(ipas),ate(ipas),g0,s0,npas,iqcf,                 
!    +             irfl,irfl2,                                          
!    +             ngper,apiv,ngnod,atev,conv)                          
           return 
        endif
     endif
  endif
! ********************************************************              
!  test the number of iterations                                        
! ********************************************************              
  if(ipas.gt.nitmax)then 
! ********************************************************              
!  too many iterations; look for resonances                             
! ********************************************************              
     da1=da(ipas) 
     if(irfl.eq.0)then 
        if(abs(dev(71,ipas))*conv.lt.0.5d0)irfl=71 
     elseif(irfl2.eq.0)then 
        if(abs(dev(71,ipas))*conv.lt.0.5d0)irfl2=71 
     endif
     call selres(amp,dev,da1,ipas,irfl,irfl2,conv,ires         &
     &               ,nde,nam0,tt,inp,nada,nada2,rho0,rho,iun10)        
     if(ires.gt.0)then 
!  successful removal, restart iteration                                
        irescl=0 
        goto 2 
     else 
!  if too many iterations, unsuccessful removal                         
!  print the last one anyway                                            
        ndiver=ndiver+1 
        gp=gv(ipas) 
        sp=sv(ipas) 
        apip=api(ipas) 
        atep=ate(ipas) 
        npas=ipas 
!           call selout(inp,namfor,t,no,                                
!    +        gv(ipas),sv(ipas),aa,pre,prsini,iqce,iqcm,iqco,           
!    +        hm,api(ipas),ate(ipas),g0,s0,ipas,iqcf,                   
!    +        irfl,irfl2,                                               
!    +        ngper,apiv,ngnod,atev,conv)                               
        return 
     endif
  else 
! ********************************************************************* 
!  otherwise, a new iteration is necessary                              
     goto 1 
  endif
! **********************************************************************
!  this is the end of the iteration loop                                
! ********************************************************************* 
! ********************************************************************* 
! 999  write(*,998) it-1,ndisc,ndiver,nada,nada2                        
!998  format(' done, no. input=',i6,/'  discarded=',i5,                 
!    +    '  non conv=',i5,' adapted=',i5,' twice ad.=',i5)             
END SUBROUTINE long9

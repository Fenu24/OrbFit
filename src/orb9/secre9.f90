! ********************************************************************* 
!   version 8.0:  with second order degree four terms                   
!   with nonlinear forced terms for eccentricity                        
!    s7 terms in degree four long periodic terms                        
! ********************************************************************* 
PROGRAM secre8 
  USE controlmod
  USE fund_const
  IMPLICIT NONE
!  forced terms for the asteroid xi_j, eta_j; j frequency index         
  DOUBLE PRECISION xi(8),et(8),xip(8),etp(8) 
!  corrections to the frequencies (for each planet)                     
  DOUBLE PRECISION cg(8),cs(8) 
!  secular theory constants for the planets                             
!  one-dim arrays; index is frequency index lambda_j, delta_j           
!  nu_j, mu_j; j is frequency index for the current planet              
  DOUBLE PRECISION pl(8),pd(8),ppni(8),ppmi(8) 
!  secular theory frequencies for the planets g_j, s_j                  
  DOUBLE PRECISION pg(8),ps(8) 
!  two-dim arrays; first index is planet index, second is freq. index   
! nu_{j,k}, mu_{j,k} with j planet ind, k freq. index                   
  DOUBLE PRECISION pmi(8,8),pni(8,8) 
!  masses, mean motions and semimajor axes of the planets               
  DOUBLE PRECISION pmv(8),pn(8),pa(8) 
!  yuasa coefficients  a and b (first order effects)                    
!  warning: here planet index is the second one, first is Yuasa subscrip
  DOUBLE PRECISION yya(4,8),yyb(27,8) 
!  nonlinear forced terms: frequencies, phases, amplitudes, forced ecc. 
  DOUBLE PRECISION fnl(nnlx),phanl(nnlx),pninl(2,nnlx),xinl(nnlx) 
!  proper elements to be used to compute frequencies
  DOUBLE PRECISION ae,ai0,aw0,ao0, ao, aw, g, s, gp, sp, g0, s0
! scalar temporaries
  DOUBLE PRECISION t,t0,tunit,aunit,pmunit,sm,conv,delti,delte,delta,deltir
  DOUBLE PRECISION e0,a0,ai0r,ro1,ro2,sig1,sig2,ani,ami,amp1,amp2,ai,ap,aq,ar,as
  DOUBLE PRECISION aa, a, sl, an, al, ag, ah
! loop indexes and other integer scalar
  INTEGER inp, ive, it1,it2,it3, j, ite, ita, nnl, iti, jp, nfr, n, ngrid, iter
!                                                                       
! ********************************************************************* 
!                                                                       
!  read controls from file secre8.opt                                   
!                                                                       
! ********************************************************************* 
  open(1,file='secre8.opt') 
  call gridsel(ngrid,ita,a0,delta,ite,e0,delte,iti,ai0,delti,aw,ao)
  deltir=delti*radeg 
  ai0r= ai0*radeg 
  open(3,file='datam.m') 
  write(3,555)ngrid,a0,delta,ai0,delti,e0,delte,aw,ao 
555 format(' data=[',i5,1x,2(f10.7,1x),4(f10.7,1x),2(f10.4,1x),']'/) 
!  version-dependent options, defining the theory used (see optxx)      
!  ive < 0 Euler; ive=3 RK2 explicit (see inte8)                        
  inp=1 
  ive=-1 
  call lonopt(inp,ive) 
! ********************************************************************* 
!                                                                       
!  initialisations, file opening, etc.                                  
!                                                                       
! **********************************************************************
  open(2,file='gs.fla',status='unknown') 
  open(4,file='forced.fla',status='unknown') 
!                                                                       
! ********************************************************************  
!  the secular perturbation theory has to be consulted to find the      
!  appropriate semimajor axis of the planets, and later the secular     
!  eccentricities and inclinations will be used as well as the          
!  fundamental secular angles (see secthe67)                            
!  the fixed time is defined by t0 (in JD), with 0 shift                
  t=0.d0                                         
  t0=2.4404005d6
  call secthe(t,t0,pmv,pa,pni,pg,pl,pmi,ps,pd,nnl,fnl,pninl,phanl)
!  data on the inner planets for computation of g0, s0                  
  call inplda(pa,pmv) 
! *******************************************************************   
! change of scale, units, fundamental constants                         
  call convun(pmv,pa,pg,ps,pn,aunit,tunit,pmunit,nplin0,nplou,fnl,nnl)
! square root of solar mass                                             
  sm=sqrt(1.d0/pmunit) 
!  conversion factor to get frequencies in arcsec/year                  
  conv=1.d0/tunit*secrad 
! ********************************************************************* 
!                                                                       
!  computations depending on semimajor axis only                        
!                                                                       
! **********************************************************************
  DO 300 it3=1,ita 
     aa=a0+delta*(it3-1) 
     a=aa/aunit 
     al=sm*sqrt(a) 
     sl=sqrt(al) 
!  asteroid mean motion                                                 
     an=sm/al**3 
! ********************************************************************* 
!   laplace's, leverrier's and yuasa's coefficients.                    
! ********************************************************************* 
     DO  j=nplin0,nplou 
        CALL llyco(a,pa(j),an,pn(j),sm,yya(1,j),yyb(1,j),pmv(j),ipla(j))
     ENDDO
! ********************************************************************* 
!   second order/degree frequencies and forced oscillation amplitudes.  
! ********************************************************************* 
     CALL g0s0(al,yya,pmv,g0,s0,nplin0,nplou0) 
!   forced oscillation amplitudes; from outer planets only              
     CALL forcedgs(g0,s0,pg,ps,pmv,pni,pmi,yya,al,xi,et,                 &
     &    nplin,nplou)                                                  
! ********************************************************************* 
!   computations depending on eccentricity and inclinations             
! *********************************************************************
      iter=0 
      DO 400 it2=1,iti 
         ai=ai0r+(it2-1)*deltir 
         DO 500 it1=1,ite 
            ae=e0+(it1-1)*delte 
!  counter of grid points                                               
            iter=iter+1 
!                                                                       
!  delaunay's and poincare's variables                                  
!                                                                       
            ag=al*sqrt(1.d0-ae*ae) 
            ah=ag*cos(ai) 
            amp1=max(2.d0*(al-ag),0.d0) 
            ap=sqrt(amp1)*cos(aw+ao) 
            aq=-sqrt(amp1)*sin(aw+ao) 
            amp2=max(2.d0*(ag-ah),0.d0) 
            ar=sqrt(amp2)*cos(ao) 
            as=-sqrt(amp2)*sin(ao) 
!                                                                       
!  proper modes.                                                        
!                                                                       
            ro1=ap/sl 
            ro2=-aq/sl 
            sig1=ar/sl 
            sig2=-as/sl 
!                                                                       
! proper actions (doubled)                                              
!                                                                       
            ani=ro1**2+ro2**2 
            ami=sig1**2+sig2**2 
! ********************************************************************* 
!                                                                       
!  Degree 4 corrections to the frequencies, a function of the proper    
!  amplitudes ami, ani and of the forced terms xi, eta.                 
!  In this version the perturbing function for planets from npli4 to npl
!  is taken into account to degree 4; however both for the planets and f
!  the asteroid the forced terms with all the frequencies               
!  (from nplin to nplou) are accounted for.                             
!                                                                       
! ********************************************************************* 
            g=g0 
            s=s0 
!  loop on planets                                                      
            DO 51 jp=npli4,nplo4 
!  proper modes of current perturbing planet                            
               DO nfr=nplin,nplou 
                  ppni(nfr)=pni(jp,nfr) 
                  ppmi(nfr)=pmi(jp,nfr) 
               ENDDO
               CALL gsfour(pmv(jp),yya(1,jp),yyb(1,jp),xi,et,ppni,ppmi, &
     &             ani,ami,cg(jp),cs(jp),nplin,nplou)                       
               g=g+2.d0/al*cg(jp) 
               s=s+2.d0/al*cs(jp) 
51          ENDDO
!   forced oscillation amplitudes; now with new frequencies g,s         
            CALL forcedgs(g,s,pg,ps,pmv,pni,pmi,yya,al,xip,etp,nplin,nplou)
! ********************************************************************* 
!                                                                       
!  conversion and output                                                
!                                                                       
            gp=g*conv 
            sp=s*conv 
            if(gp.ge.1.d3)gp=999.999d0 
            if(gp.le.-1.d3)gp=-999.999d0 
            if(sp.ge.1.d3)sp=999.999d0 
            if(sp.le.-1.d3)sp=-999.999d0 
            write(2,601)gp,sp 
601         format(2f12.4) 
!                                                                       
            if(xip(6).gt.1.d0)then 
               xip(6)=1.d0 
            elseif(xip(6).lt.-1.d0)then 
               xip(6)=-1.d0 
            endif
            if(etp(6).gt.1.d0)then 
               etp(6)=1.d0 
            elseif(etp(6).lt.-1.d0)then 
               etp(6)=-1.d0 
            endif
            if(xip(5).gt.1.d0)then 
               xip(5)=1.d0 
            elseif(xip(5).lt.-1.d0)then 
               xip(5)=-1.d0 
            endif
            if(etp(5).gt.1.d0)then 
               etp(5)=1.d0 
            elseif(etp(5).lt.-1.d0)then 
               etp(5)=-1.d0 
            endif
            write(4,602)xip(5),etp(5),xip(6),etp(6) 
602         format(2f10.6) 
! ********************************************************************* 
!                                                                       
!   end of grid loops                                                   
!                                                                       
500      ENDDO
400   ENDDO
!     write(*,*)it3                                                     
300 ENDDO
   n=iter 
   close(2) 
   close(4) 
 END PROGRAM secre8

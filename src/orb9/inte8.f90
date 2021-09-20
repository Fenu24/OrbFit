! ===================================================================== 
!   inte8                                                               
!   euler method                                                        
! ===================================================================== 
!   first line: functions of a only                                     
!   second line: functions of planetary orbits                          
!   third  line: controls for resonance removal                         
!   fourth line: linear and current elements,                           
!   fifth line: output: freq., divisors, amplitudes                     
!   sixth line: new elements                                            
SUBROUTINE inte8(ive,yya,yyb,xi,et,al,aa,g0,s0,                   &
     &   pni,pmi,pg,ps,pmv,pninl,nnl,inl,cld,sld,phanl,                 &
     &   irfl,irfl2,                                                    &
     &   rho0,rho,                                                      &
     &   g,s,dev,nde,amp,ampi,ampe,                                     &
     &   rho1) 
  USE controlmod
  IMPLICIT NONE  
! ===================================================================== 
! INPUT
  INTEGER, INTENT(IN) :: ive 
  DOUBLE PRECISION, INTENT(IN) :: al,aa,g0,s0 ! functions of a only
! number of nonlinear terms
  INTEGER, INTENT(IN) :: nnl
!  nonlinear forced terms: phases, amplitudes,                          
  DOUBLE PRECISION, INTENT(IN) :: phanl(nnlx),pninl(2,nnlx) 
!  secular theory constants for the planets                             
!  two-dim arrays; first index is planet index, second is freq. index   
! nu_{j,k}, mu_{j,k} with j planet ind, k freq. index                   
  DOUBLE PRECISION, INTENT(IN) :: pmi(8,8),pni(8,8) 
!  secular theory frequencies for the planets g_j, s_j                  
  DOUBLE PRECISION, INTENT(IN) :: pg(8),ps(8) 
!  masses of the planets                                                
  DOUBLE PRECISION, INTENT(IN) :: pmv(8) 
! ===================================================================== 
!  yuasa coefficients  a and b (first order effects)                    
!  warning: here planet index is the second one, first is Yuasa subscript
  DOUBLE PRECISION, INTENT(IN) ::  yya(4,8),yyb(28,8) 
!     trigonometric functions                                           
  DOUBLE PRECISION, INTENT(IN) :: cld(-2:2,-2:2,-2:2,-2:2),sld(-2:2,-2:2,-2:2,-2:2) 
! ======================                                                
!  forced terms for the asteroid xi_j, eta_j; j frequency index         
  DOUBLE PRECISION, INTENT(IN) :: xi(8),et(8) 
!  current linear proper elements 
  DOUBLE PRECISION, INTENT(IN) :: rho(4),rho0(4)
! OUTPUT
  INTEGER, INTENT(IN) :: irfl,irfl2 ! resonance removal flags
  DOUBLE PRECISION, INTENT(OUT) :: g,s ! updated frequencies
!  cross reference index for nonlinear forced terms                     
  INTEGER, INTENT(OUT) :: inl(nnlx) 
!  divisors, amplitudes of largest terms                                
  DOUBLE PRECISION, INTENT(OUT) :: dev(ndex),amp(ndex),ampe(ndex),ampi(ndex)
! output new iteration elements
  DOUBLE PRECISION, INTENT(INOUT) ::  rho1(4) 
! ??????
  INTEGER, INTENT(OUT) :: nde
! ====================== proper elems. per iteration====================
!  vector of Poincare' elements: derivative with respect to mass,        
  DOUBLE PRECISION drh1(4),drh2(4)
  DOUBLE PRECISION g1,s1, h ! corrected frequencies, integration step
  INTEGER i ! loop index
! ===================================================================== 
!  euler step                                                           
  if(ive.lt.0)then 
!   computation of the right hand side of Lie series diff. euqation     
!   first line: functions of a only                                     
!   second line: functions of planetary orbits, time independent        
!   third line: planetary functions dependent upon time                 
!   fourth line: frequencies from linear theory                         
!   fifth line: controls for resonance removal                          
!   sixth line: current elements, freq., divisors, amplitudes           
!   seventh line: derivatives of elements w.r. to mass of planets       
     call iterat(yya,yyb,xi,et,al,aa,                               &
     &      pni,pmi,pg,ps,pmv,pninl,nnl,inl,                            &
     &      cld,sld,phanl,                                              &
     &      g0,s0,                                                      &
     &      irfl,irfl2,                                                 &
     &      rho,g,s,dev,nde,ampe,ampi,amp,                              &
     &      drh1)                                                       
! ********************************************************************* 
!   proper modes corrected directly (without going back to non-proper   
!   variables): euler method                                            
! ********************************************************************* 
     h=pmv(5) 
     DO i=1,4 
        rho1(i)=rho0(i)+h*drh1(i) 
     ENDDO
! ===================================================================== 
  else 
! ===================================================================== 
!  Runge-Kutta-2 (explicit)                                             
     if(ive.eq.3)then 
!  compute right hand side at intermediate (euler) value of rho         
        call iterat(yya,yyb,xi,et,al,aa,                            &
     &         pni,pmi,pg,ps,pmv,pninl,nnl,inl,                         &
     &         cld,sld,phanl,                                           &
     &         g0,s0,                                                   &
     &         irfl,irfl2,                                              &
     &         rho1,g1,s1,dev,nde,ampe,ampi,amp,                        &
     &         drh2)                                                    
        h=pmv(5) 
        do i=1,4 
           rho1(i)=rho0(i)+(h*0.5d0)*(drh1(i)+drh2(i)) 
        enddo
!  output frequencies are last computed (consistent with divisors dev)  
        g=g1 
        s=s1 
     else 
        write(*,*)' option not yet supported, ive=',ive 
        stop 
     endif
  endif
END SUBROUTINE inte8

! ===================================================================== 
!   iterat: one iteration of fixed point map                            
!   first line: functions of a only                                     
!   second line: functions of planetary orbits, time independent        
!   third line: planetary functions dependent upon time                 
!   fourth line: frequencies from linear theory                         
!   fifth line: controls for resonance removal                          
!   sixth line: current elements, freq., divisors, amplitudes           
!   seventh line: derivatives of elements w.r. to mass of planets       
! ===================================================================== 
SUBROUTINE iterat(yya,yyb,xi,et,al,aa,                            &
     &   pni,pmi,pg,ps,pmv,pninl,nnl,inl,                               &
     &   cld,sld,phanl,                                                 &
     &   g0,s0,                                                         &
     &   irfl,irfl2,                                                    &
     &   rho,g,s,de,nde,ampe,ampi,amp,                                  &
     &   drh)                                                           
! ===================================================================== 
  USE controlmod
  USE terms9
  IMPLICIT NONE
! ===================================================================== 
! INPUT
! semimajor axis etc.
   DOUBLE PRECISION, INTENT(IN) :: al, aa, g0, s0
! ===================================================================== 
!  yuasa coefficients  a and b (first order effects)                    
!  warning: here planet index is the second one, first is Yuasa subscrip
   DOUBLE PRECISION, INTENT(IN) :: yya(4,8),yyb(28,8) 
!  forced terms for the asteroid xi_j, eta_j; j frequency index         
   DOUBLE PRECISION, INTENT(IN) :: xi(8),et(8) 
! ===================================================================== 
!  secular theory constants for the planets                             
!  two-dim arrays; first index is planet index, second is freq. index   
! nu_{j,k}, mu_{j,k} with j planet ind, k freq. index                   
   DOUBLE PRECISION, INTENT(IN) :: pmi(8,8),pni(8,8) 
!  one-dim arrays; index is frequency index lambda_j, delta_j           
!  secular theory frequencies for the planets g_j, s_j                  
   DOUBLE PRECISION, INTENT(IN) :: pg(8),ps(8) 
!  masses of the planets                                                
   DOUBLE PRECISION, INTENT(IN) :: pmv(8) 
!  trigonometric functions                                           
   DOUBLE PRECISION, INTENT(IN) :: cld(-2:2,-2:2,-2:2,-2:2),sld(-2:2,-2:2,-2:2,-2:2) 
! ===================================================================== 
! removed resonance flags, no. divisors
   INTEGER, INTENT(IN) :: irfl,irfl2
!  nonlinear forced terms, actual number
   INTEGER, INTENT(IN) :: nnl 
!  nonlinear forced terms: phases, amplitudes, forced ecc. 
   DOUBLE PRECISION, INTENT(IN) :: phanl(nnlx),pninl(2,nnlx) 
   DOUBLE PRECISION, INTENT(IN) :: rho(4) ! proper elements
! ===================================================================== 
! OUTPUT
!  divisors, amplitudes of largest terms                 
   DOUBLE PRECISION, INTENT(OUT) :: de(ndex),ampe(ndex),ampi(ndex),amp(ndex)
   DOUBLE PRECISION, INTENT(OUT) :: g,s ! frequencies
   DOUBLE PRECISION, INTENT(OUT) :: drh(4) ! correction to proper elements
! no. divisors, cross reference index for nonlinear forced terms
   INTEGER, INTENT(OUT) :: nde, inl(nnlx)
! END INTERFACE 
!  nonlinear forced terms
   DOUBLE PRECISION xinl(nnlx,8)
! ===================================================================== 
!  degree four corrections                                              
   DOUBLE PRECISION fp(4,ndex),f(4,ndex),sr(4) 
!  nu_j, mu_j; j is frequency index for the current planet              
   DOUBLE PRECISION ppni(8),ppmi(8) 
!  corrections to the frequencies (for each planet)                     
   DOUBLE PRECISION cg(8),cs(8)
! scalar temporaries
   DOUBLE PRECISION h, ro1,ro2,sig1,sig2,ani,ami
! loop indexes
   INTEGER nfr, jp , i, j,n 
! =================================================================**** 
   save 
! ===================================================================== 
!   proper actions (doubled)                                            
   ro1=rho(1) 
   ro2=rho(2) 
   sig1=rho(3) 
   sig2=rho(4) 
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
   DO  jp=npli4,nplo4 
!  proper modes of current perturbing planet                            
      DO nfr=nplin,nplou 
         ppni(nfr)=pni(jp,nfr) 
         ppmi(nfr)=pmi(jp,nfr) 
      ENDDO
      call gsfour(pmv(jp),yya(1,jp),yyb(1,jp),xi,et,ppni,ppmi,        &
     &         ani,ami,cg(jp),cs(jp),nplin,nplou)                       
      g=g+2.d0/al*cg(jp) 
      s=s+2.d0/al*cs(jp) 
   ENDDO
! ********************************************************************* 
!                                                                       
!   Computation of secular perturbations resulting from degree 4 terms  
!   in the perturbing function from planets from npli4 to nplo4;        
!   they are function of the proper elements and of the forcing         
!   terms; however only four planetary frequencies (two for the periheli
!   and two for the nodes) are taken into account; in this version they 
!   are g5, g6, s7, and s6.                                             
!   The frequencies appear only in the denominators de.                 
!                                                                       
! ********************************************************************* 
   CALL denom(g,s,pg,ps,de,nde,inl) 
!  perturbation accumulators are set to zero  
   fp=0.d0
   sr=0.d0                          
!  loop on perturbing planets                                           
   DO 31 jp=npli4,nplo4 
      call terms(al,ro1,ro2,sig1,sig2,pni(jp,5),pni(jp,6),            &
     &    pmi(jp,7),pmi(jp,6),xi(5),xi(6),et(7),et(6),                  &
     &    sld,cld,yya(1,jp),yyb(1,jp),de,pmv(jp),f)                     
! ********************************************************************* 
!    nonlinear forced terms (for each planet) are added here            
!    they are computed with current value of g                          
! ********************************************************************* 
      call nlforc(nnl,pninl,xinl,nnlx,yya) 
      DO n=1,nnl 
         f(2,inl(n))=f(2,inl(n))-xinl(n,jp)*cos(phanl(n)) 
         f(1,inl(n))=f(1,inl(n))+xinl(n,jp)*sin(phanl(n)) 
      ENDDO
! ********************************************************************* 
!                                                                       
!     corrections-final equations. (with divisors removed if required)  
!                                                                       
! ********************************************************************* 
      DO 20 j=1,nde 
         DO 10 i=1,4 
            fp(i,j)=fp(i,j)+pmv(jp)*f(i,j)/de(j) 
!  forced removal of some term(s)                                       
!  this is due to the relationship of the z2 term with other divisors   
            if(j.eq.62)then 
               if(irfl.ne.10.and.irfl2.ne.10.and.                       &
     &         irfl.ne.2.and.irfl2.ne.2.and.                            &
     &         irfl.ne.57.and.irfl2.ne.57.and.                          &
     &         irfl.ne.5.and.irfl2.ne.5.and.                            &
     &         irfl.ne.48.and.irfl2.ne.48.and.                          &
     &         irfl.ne.24.and.irfl2.ne.24.and.                          &
     &         irfl.ne.62.and.irfl2.ne.62)then                          
                  if(aa.lt.2.48d0)then 
                     sr(i)=sr(i)+pmv(jp)*dd(j)*f(i,j)/de(j) 
                  endif
               endif
            elseif(j.ne.irfl.and.j.ne.irfl2)then 
               sr(i)=sr(i)+pmv(jp)*dd(j)*f(i,j)/de(j) 
            endif
10       ENDDO
20    ENDDO
!  perturbation from all planets have been accumulated in sr, fp        
31 ENDDO
! ================================================================      
!  corrections without mass of Jupiter as factor                        
   h=pmv(5) 
   drh(2)=-sr(1)/al/h 
   drh(1)=sr(2)/al/h 
   drh(4)=-sr(3)/al/h 
   drh(3)=sr(4)/al/h 
! ********************************************************************* 
!    amplitudes stored for resonance test                               
! ********************************************************************* 
   DO 40 j=1,nde 
      ampe(j)=0.d0 
      ampi(j)=0.d0 
      DO i=1,2 
         ampi(j)=ampi(j)+fp(i+2,j)**2 
         ampe(j)=ampe(j)+fp(i,j)**2 
      ENDDO
      amp(j)=sqrt(ampe(j)+ampi(j))/al 
      ampe(j)=sqrt(ampe(j))/al 
      ampi(j)=sqrt(ampi(j))/al*dd(j) 
40 ENDDO
!====================================================================== 
!     return from iterat                                                
 END SUBROUTINE iterat

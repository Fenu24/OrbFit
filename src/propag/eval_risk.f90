MODULE eval_risk
!============================================================
! To compute impact probability and Palermo Scale value
! transformed into module, A. Milani August 2005          
! Created by Giacomo Tommei, 25 November 2002 
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: palermo, torino ! public functions
  LOGICAL, PUBLIC :: massgiven ! is mass available from physical observations?
  DOUBLE PRECISION, PUBLIC ::  givenmass ! then this is the value in 
CONTAINS

! ======================================                                
! computation of impact risk scales                                     
DOUBLE PRECISION FUNCTION palermo(U,hmagn,prob,dt,          &
     &   mass,v_imp,energy,e_tilde,fb,rel_prob,iunwarn)
  USE fund_const              
  IMPLICIT NONE 

! INPUT                                                                 
! velocity at infinity, H magnitude, probability, difference in time)   
  DOUBLE PRECISION, INTENT(IN) :: U, hmagn, prob, dt 
! output warnign file                                                   
  INTEGER, INTENT(IN) :: iunwarn 
! OUTPUT                                                                
! mass, impact velocity, energy, expected energy                        
  DOUBLE PRECISION, INTENT(OUT) :: mass,v_imp,energy,e_tilde 
! frequency of impacts with given energy, relative probability          
  DOUBLE PRECISION, INTENT(OUT) :: fb, rel_prob 
! ==================end interface===========================            
! default values: rho*p_v**(-3/2) for 'standard' NEA (to be replaced if 
! the units are kg/km**3                                                
  DOUBLE PRECISION rhopv32 
  DATA rhopv32 / 43.d12/ 
  DOUBLE PRECISION v_earth,v_escape 
  DATA v_earth /29.78d0/ 
  DATA v_escape /11.18d0/ 
! velocity of impact squared, difference in time years, conversion au/d 
  DOUBLE PRECISION v_imp2,dyear,convfac 
  INTEGER iunwarn0
! =======================================================  
  iunwarn0=abs(iunwarn)      
  IF(massgiven)THEN
     mass=givenmass
  ELSE       
! mass in kg; 1329 is a mistery number, where is it from? reference needed
     mass=pig/6.d0 *rhopv32 * 1.329d3**3 * 10**(-3.d0*hmagn/5.d0) 
  ENDIF
! impact velocity in km/s                                               
  convfac=aukm/86400.d0 
  v_imp2=(U*convfac)**2+v_escape**2 
  v_imp=sqrt(v_imp2) 
! energy of impact in Joule, in megaton                                 
  energy=mass*v_imp2/2.d0*1.d6 
  energy=energy/4.2e15 
! expected energy                                                       
  e_tilde=energy*prob 
! probability density of impact at the given energy                     
  fb=3e-2*energy**(-0.8d0) 
! relative probability                                                  
  dyear=dt/365.25d0 
  rel_prob=prob/(fb*dyear) 
! Palermo logarithmic scale                                             
  palermo=log10(rel_prob) 
! output: debug
  IF(iunwarn.lt.0)THEN
     WRITE(*,*)' risk case analysed for Palermo scale' 
     WRITE(*,*)'mass, impact velocity, energy' 
     WRITE(*,100) mass, v_imp, energy 
100  FORMAT(1p,d9.2,' kg ',0p,f5.2,' km/s ', 1p,d9.2,' MT') 
     WRITE(*,*)' expected E, background prob, ratio, scale value' 
     WRITE(*,101) e_tilde, fb, rel_prob, palermo 
101  FORMAT(1p,d9.2,' MT ',d9.2,' backgr ',d9.2,' ratio ',             &
     &     0p,f6.2,' scale')  
  ENDIF
  WRITE(iunwarn0,*)' risk case analysed for Palermo scale' 
  WRITE(iunwarn0,*)'mass, impact velocity, energy' 
  WRITE(iunwarn0,100) mass, v_imp, energy                              
  WRITE(iunwarn0,*)' expected E, background prob, ratio, scale value' 
  WRITE(iunwarn0,101) e_tilde, fb, rel_prob, palermo           
END FUNCTION palermo
!#######################################################################
!                     TORINO_SCALE
!#######################################################################
! torino_scale:  input-> impact probability, energy of impact
!#######################################################################

INTEGER FUNCTION torino(IP,energy)
  DOUBLE PRECISION, INTENT(IN) :: IP, energy

  IF (IP.ge.0d0.and.IP.le.1d0.and.energy.ge.0d0.and.energy.le.10**8) THEN
     IF ((log10(IP)).le.-2d0.and.(log10(energy)).gt.0d0.and.(log10(energy)).gt.&
            &(-4d0-3d0/2d0*log10(IP)).and.(log10(energy)).le.(-1d0-3d0/2d0*log10(IP))) THEN 
        torino=1
     ELSEIF ((log10(IP)).le.-2d0.and.(log10(energy)).gt.(-1d0-3d0/2d0*log10(IP)).and.&
    &(log10(energy)).le.(2d0-3d0/2d0*log10(IP))) THEN
        torino=2
     ELSEIF ((log10(IP)).gt.-2d0.and.IP.le.0.99d0.and.(log10(energy)).gt.0d0.and.(log10(energy)).le.2d0) THEN
        torino=3
     ELSEIF ((log10(IP)).gt.-2d0.and.(log10(energy)).ge.2d0.and.(log10(energy)).le.(2d0-3d0/2d0*log10(IP))) THEN
        torino=4
     ELSEIF (IP.le.0.99d0.and.(log10(energy)).le.5d0.and.(log10(energy)).gt.(2d0-3d0/2d0*log10(IP))) THEN
        torino=5
     ELSEIF ((log10(IP)).le.-2d0.and.(log10(energy)).gt.(2d0-3d0/2d0*log10(IP))) THEN
        torino=6
     ELSEIF ((log10(IP)).gt.-2d0.and.IP.le.0.99d0.and.(log10(energy)).gt.5d0) THEN
        torino=7
     ELSEIF (IP.gt.0.99d0.and.(log10(energy)).gt.0d0.and.(log10(energy)).le.2d0) THEN
        torino=8
     ELSEIF (IP.gt.0.99d0.and.(log10(energy)).gt.2d0.and.(log10(energy)).le.5d0) THEN
        torino=9
     ELSEIF (IP.gt.0.99d0.and.(log10(energy)).gt.5d0) THEN
        torino=10
     ELSE    
        torino = 0    
     ENDIF
  ELSE
     torino= -1
  ENDIF
END FUNCTION torino



END MODULE eval_risk

DOUBLE PRECISION FUNCTION prob_1dim(be,s0,r0)
  IMPLICIT NONE 
!     input data                                                        
  DOUBLE PRECISION, INTENT(IN) :: be,s0,r0 ! impact cross
! section radius, stretching, minimum LOV distance 
!     output: impact probability prob, assuming width=0  
  DOUBLE PRECISION chord
  DOUBLE PRECISION, PARAMETER :: span_sigma=70.d0
  IF(r0.gt.be)THEN
     prob_1dim=0.d0
     RETURN
  ENDIF
  chord=2.d0*sqrt(be**2-r0**2)
  prob_1dim=chord/(span_sigma*s0)
END FUNCTION prob_1dim

!============================================================          
DOUBLE PRECISION FUNCTION prob(csi0,zeta0,sigma0,w0,s0,r0) 
  IMPLICIT NONE 
!     input data                                                        
  DOUBLE PRECISION, INTENT(IN) :: csi0,zeta0,sigma0,w0,s0,r0 
!     output: impact probability prob                                   
! ============ end interface =============                          
!     common per passare ad f                                           
  DOUBLE PRECISION csi,zeta,sigma,w,s,r 
  COMMON/probabil/csi,zeta,sigma,w,s,r 
  DOUBLE PRECISION eee 
  INTEGER nnn 
!     per dqags                                                         
  INTEGER limx,limx4 
  PARAMETER (limx=100,limx4=4*limx) 
  INTEGER ier,limit,iwork(limx),lenw,last 
!     function evaluations                                              
  INTEGER neval 
  DOUBLE PRECISION epsabs,epsrel,abserr,work(limx4) 
!     integration interval in the csi variable                          
  DOUBLE PRECISION rmax,rmin 
!     output of dqags                                                   
  DOUBLE PRECISION rm 
  EXTERNAL ff_prob 
!     =======================================================           
!     copy into common                                                  
  csi=csi0 
  zeta=zeta0 
  sigma=sigma0 
  w=w0 
  s=s0 
  rmax=min(r0,csi0+8*w0) 
  rmin=max(-r0,csi0-8*w0) 
  IF(rmin.gt.rmax)THEN 
     prob=0.0d0 
     RETURN 
  ENDIF
  r=r0 
! check for essentially 1-dimensional case                              
                                                                        
!     set accumulators to zero                                          
  prob=0.d0 
  eee=0.d0 
  nnn=0 
                                                                        
!     preparativi chiamata dqags                                        
  epsabs=1.d-9 
  epsrel=1.d-5 
  limit=limx 
  lenw=limx4 
                                                                        
!     perturbing function                                               
  CALL dqags(ff_prob,rmin,rmax,epsabs,epsrel,rm,abserr,neval,ier,limit   &
     &     ,lenw,last,iwork,work)                                       
                                                                        
  prob = rm 
  eee = abserr 
  nnn = neval 
                                                                        
END FUNCTION prob
                                                                        
!     ===========================================================       
DOUBLE PRECISION FUNCTION ff_prob(u) 
  IMPLICIT NONE 
!     u variable                                                        
  double precision u 
!     ========== end interface ===========                              
  DOUBLE PRECISION csi,zeta,sigma,w,s,r 
!     common per passare ad f                                           
  COMMON/probabil/csi,zeta,sigma,w,s,r 
  DOUBLE PRECISION espo1 
!     common per passare ad f                                           
  COMMON/proble/espo1 
!     estremi di integrazione                                           
  DOUBLE PRECISION fstext,sndext 
!     per dqags                                                         
  INTEGER limx,limx4 
  PARAMETER (limx=100,limx4=4*limx) 
  INTEGER neval,ier,limit,iwork(limx),lenw,last 
  DOUBLE PRECISION epsabs,epsrel,abserr,work(limx4),result 
  EXTERNAL f_prob 
!     ===========================================================       
                                                                        
!     calcoli dipendenti solo da u                                      
  espo1 = -(1.d0/2.d0)*((u-csi)/w)**2 
                                                                        
!     estremi integrale                                                 
  fstext = -dsqrt(r**2-u**2) 
  sndext = dsqrt(r**2-u**2) 
                                                                        
!     preparativi per dqags                                             
  epsabs=1.d-9 
  epsrel=1.d-5 
  limit=limx 
  lenw=limx4 
                                                                        
!     integrale                                                         
  call dqagsc(f_prob,fstext,sndext,epsabs,epsrel,result,abserr,neval,ier,&
     &   limit,lenw,last,iwork,work)                                    
                                                                        
!     risultato                                                         
  ff_prob=result 

END FUNCTION ff_prob
                                                                        
!     ===========================================================       
DOUBLE PRECISION FUNCTION f_prob(uu)
  USE fund_const 
  IMPLICIT NONE 
!     uu variable                                                       
  double precision uu 
!     ========== end interface ===========                              
  DOUBLE PRECISION csi,zeta,sigma,w,s,r 
  COMMON/probabil/csi,zeta,sigma,w,s,r 
  DOUBLE PRECISION espo1 
  COMMON/proble/espo1 
  DOUBLE PRECISION espo2 
!     ===========================================================       
                                                                        
  espo2 = -(1.d0/2.d0)*((uu-zeta)/s+sigma)**2 
                                                                        
!     FUNZIONE INTEGRANDA                                               
  f_prob = (1.d0/(2.d0*pig*s*w))*exp(espo1 + espo2) 
                                                                        
END FUNCTION f_prob
                                                                        

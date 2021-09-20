! ============================================                          
!  FITSUBS                                                              
! ============================================                          
! routines used by orbit determination interactive programs, such as fit
!                                                                       
!  CONTAINS:                                                            
!                                                                       
!   whichcor   selection of parameters to be solved for                 
!   chereq     check availability of data required                      
!   cheobs     id for observations                                      
!   chetim     id for time series such as JPL ephem.)                   
!   seleph     interactive selection of interval                        
!   asstim        "          "       of time                            
!   asscbd        "          "       of confidence boundary             
!   orbsel        "          "       of orbit                           
!   start      initial guess for identification (by interpolation)      
! ===================================================================   
! WHICOR                                                                
! ===================================================================   
! interrogation routine for inew, icor                                  
SUBROUTINE whicor(inter,nd,icor,ncor,inew) 
  implicit none 
  INTEGER, INTENT(IN) :: inter, nd 
  INTEGER, INTENT(OUT) :: icor(nd),ncor,inew
! end interface
  integer iansw,i 
! Choose method (NOW FIXED AT PSEUDO-NEWTON)                            
  inew=2 
!60   write(*,*)' 1=true Newton 2=pseudo Newton'                        
!     read(*,*)inew                                                     
!     if(inew.ne.1.and.inew.ne.2)then                                   
!           write(*,*)'This we have not invented yet'                   
!           goto 60                                                     
!     endif                                                             
! interactive/automatic default version                                 
  if(inter.eq.0)then 
     icor(1:nd)=1 
     ncor=nd 
     return 
  else 
! interactive version                                                   
!                                                                       
! Component of orbital element vector that need to be corrected         
     ncor=0 
     do 63 i=1,nd 
        write(*,*)'Element:  ',i,   ' 1=correct, 0=no' 
        read(*,*) iansw 
        if(iansw.ne.0)then 
           ncor=ncor+1 
           icor(i)=iansw 
        else 
           icor(i)=0 
        endif
63   enddo
  endif
END SUBROUTINE whicor
! ====================================================================  
! CHEREQ                                                                
! ====================================================================  
! checking the availability of the required data, for propagation       
! and preciction of observations                                        
! ====================================================================  
SUBROUTINE chereq(icov,ini,cov,t,iunout,ok) 
  implicit none 
  integer icov 
  logical ini,cov 
  integer iunout
  double precision t 
  logical ok 
! ===================================================================   
! logical check: initial conditions available                           
! ===================================================================   
  if(ini)then 
     IF(iunout.gt.0)write(iunout,222) t 
     ok=.true. 
  else 
     write(*,*)' initial conditions not available' 
     ok=.false. 
  endif
! ===================================================================   
! logical check: if covariance has to be propagated,                    
! it needs to be available for the initial time                         
! ===================================================================   
  if(icov.gt.1)then 
     if(.not.cov)then 
        write(*,*)' initial covariance not available' 
        ok=.false. 
     endif
  endif
  222 format('  Initial epoch (MJD): ',f8.1) 
END SUBROUTINE chereq
! ====================================================================  
! CHEOBS                                                                
! ====================================================================  
! check availability of observations and initial condition              
! ===================================================================   
SUBROUTINE cheobs(obs,ini,ok) 
  logical obs,ini,ok 
  if(ini)then 
     if(obs)then 
        ok=.true. 
     else 
        write(*,*)'missing observations for this arc' 
        ok=.false. 
     endif
  else 
     write(*,*)'missing initial conditions for this arc' 
     ok=.false. 
  endif
END SUBROUTINE cheobs
! ====================================================                  
! CHETIM check availability of required time-dependent data             
!     input: t1,t2 dates MJD; assumed t1.le.t2                          
!     ok: .true. if available                                           
! ====================================================                  
SUBROUTINE chetim(t1,t2,ok) 
  implicit none 
  double precision t1,t2 
  logical ok 
! ======== time spans for JPL data etc. =========                       
  include 'timespan.h90' 
  ok=.true. 
! ==================================================================    
! check availability of JPL ephemerides                                 
  if(t1.lt.tejpl1.or.t2.gt.tejpl2)then 
     write(*,*)' JPL ephemerides not available for =',t1,t2 
     write(*,*)' but only for interval ',tejpl1,tejpl2 
     ok=.false. 
  endif
! ===================================================================   
! check availability of ET-UT table                                     
  if(t1.lt.temut1.or.t2.gt.temut2)then 
     if(temute) then 
!           write(*,*)' ET-UT not available for ',t1,t2                 
!           write(*,*)' but only for interval ',temut1,temut2           
!           write(*,*)' however, extrapolation will be used'            
     else 
        write(*,*)' ET-UT not available for ',t1,t2 
        write(*,*)' but only for interval ',temut1,temut2 
        ok=.false. 
     endif
  endif
END SUBROUTINE chetim
! ===================================================                   
! SELEPH                                                                
! select time interval, step                                            
SUBROUTINE seleph(tut1,tdt1,tut2,tdt2,dt,idsta)
  USE station_coordinates
  IMPLICIT NONE 
! output                                                                
  DOUBLE PRECISION  tut1,tdt1,tut2,tdt2,dt 
  INTEGER idsta 
! times in various formats used internally                              
  CHARACTER*3 scale 
  INTEGER mjd,mjdtdt 
  DOUBLE PRECISION sec,sectdt
 ! output                                                                
  CHARACTER*3 obsstr 
  WRITE(*,*)' Initial time (MJD UTC)?' 
  READ(*,*)tut1 
  WRITE(*,*)' Final time (MJD UTC)?' 
  READ(*,*)tut2 
  WRITE(*,*)' Time interval (days)?' 
  READ(*,*)dt 
  WRITE(*,*)' Observatory code?' 
  READ(*,*)obsstr 
  CALL statcode(obsstr,idsta)
  scale='UTC' 
! =========== TIME CONVERSION ================                          
! starting time                                                         
  mjd=tut1 
  sec=(tut1-mjd)*86400 
  call cnvtim(mjd,sec,scale,mjdtdt,sectdt,'TDT') 
  tdt1=sectdt/86400.d0+float(mjdtdt) 
! stopping time                                                         
  mjd=tut2 
  sec=(tut2-mjd)*86400 
  call cnvtim(mjd,sec,scale,mjdtdt,sectdt,'TDT') 
  tdt2=sectdt/86400.d0+float(mjdtdt) 
END SUBROUTINE seleph
! ====================================================                  
! ASSTIM assign time for predictions                                    
! operations controlled by flag icov                                    
!  ICOV=1,2,3  ask the user to assign prediciton time t1                
!  ICOV=4      ask the user to assign observation number im             
!                      then t1=tau(im)                                  
!  for ICOV.ne.4 also station code ids and obs. type iob1               
!        have to be assigned by the user                                
! ====================================================                  
SUBROUTINE asstim(icov,typ,tau,tut,idsta,m,mall,im,type1,t1,tut1,ids)
  USE util_suit
  USE station_coordinates, ONLY: statcode
  IMPLICIT NONE 
  INTEGER m,mall
! ==============input=============  
  INTEGER icov,mp,idsta(mall)
  CHARACTER*(1) typ(mall)                                   
  DOUBLE PRECISION tau(mall),tut(mall) 
! ==============output=============                                     
  INTEGER im,ids
  CHARACTER*(1) type1
  CHARACTER*3 ids_aux ! observatory alpha-numeric code
! end interface
  INTEGER iob1
  DOUBLE PRECISION t1,tut1 
! time conversion                                                       
  INTEGER mjd1,mjd2
  DOUBLE PRECISION sec1,sec2 
! characters for menu                                                   
  CHARACTER*20 menunam 
! calendar to MJD conversion
  INTEGER iy,imo,iday 
  DOUBLE PRECISION hr,tjm1
! ===================================================================   
  if(mall.lt.m)then 
     write(*,*)'asstim: this should not happen, m,mall ',m,mall 
     mp=0 
  else 
     mp=mall-m 
  endif
! assign observation time                                               
  if(icov.eq.5)then 
! compare confidence boundary with observations                         
1    write(*,*)' observed arcs: from, to, no. obs' 
     write(*,182)tau(1),tau(m),m 
182  format('arc 1: ',2f8.1,i6) 
     IF(mall.gt.m)THEN 
        write(*,183)tau(m+1),tau(mall),mp 
183     format('arc 2: ',2f8.1,i6) 
     ENDIF
     write(*,*)' observation number?   ' 
     read(*,*)im 
     if(im.lt.0.or.im.gt.mall)then 
        write(*,*)' observation no. im=',im,' not available' 
        goto 1 
     endif
     t1=tau(im) 
     tut1=tut(im) 
     ids=idsta(im)
     type1=typ(im) 
  else 
! dummy obs. number (not used, to avoid out of bounds)                  
     im=1 
! assign arbitrary time   
     WRITE(*,*)'simulated obs. UTC calendar date, year, month, day, hour (with fraction)?' 
     READ(*,*)iy,imo,iday,hr 
     tut1= tjm1(iday,imo,iy,hr)
!     write(*,*)' give time of prediction (MJD)   ' 
!     read(*,*)t1 
! universal time of the required observation                            
     mjd1=FLOOR(tut1) 
     sec1=(tut1-float(mjd1))*86400.d0 
     CALL cnvtim(mjd1,sec1,'UTC',mjd2,sec2,'TDT') 
     t1=sec2/86400.d0+float(mjd2) 
!         write(*,*)t1,tut1                                             
! assign observation type     
     menunam='predtype'                                          
2    CALL menu(iob1,menunam,5,'Observation type?=',                 &
     &      'optical (alpha, delta)=',                                  &
     &      'radar range=',                                     &
     &      'radar range rate=',                               &
     &      'satellite (alpha,delta)=',                                 &
     &      'proper motion (adot, ddot)=')
     IF(iob1.eq.0)THEN
        WRITE(*,*) 'please select observation type'
        GOTO 2
     ENDIF
! convert code to type
     IF(iob1.eq.1)THEN
        type1='O'
     ELSEIF(iob1.eq.2)THEN
        type1='R'
     ELSEIF(iob1.eq.3)THEN
        type1='V'
     ELSEIF(iob1.eq.4)THEN
        type1='S'
     ELSEIF(iob1.eq.5)THEN
! type P does not exist, but is udes by fobpre to predict proper motion
        type1='P'
     ENDIF
! assign observatory code                                               
     write(*,*) 
     write(*,*)' observatory code (geocenter=500)?   ' 
     read(*,*) ids_aux
     CALL statcode(ids_aux,ids)
! secret nationalistic feature                                          
     if(ids.lt.0)then 
        ids=599 
     endif
  endif
END SUBROUTINE asstim
                                                                        
! ===================================================================   
!  ASSCBD                                                               
! ===================================================================   
! alpha, delta, magnitude, covariance and confidence boundary;          
! input specification of set of points                                  
! ===================================================================   
      SUBROUTINE asscbd(iun20,npox,npo,sigma,ibv) 
      IMPLICIT NONE 
      INTEGER npox,npo,ibv,iun20 
      DOUBLE PRECISION sigma 
      WRITE(*,*)' How many sigmas?   ' 
      READ(*,*) sigma 
    1 WRITE(*,*)' 1=confidence boundary 2=line of max variation 0=auto' 
      READ(*,*) ibv 
      IF(ibv.ne.1.and.ibv.ne.2.and.ibv.ne.0)THEN 
         WRITE(*,*)' option not understood ',ibv 
         GOTO 1 
      ENDIF 
      WRITE(*,*)' how many points (even, please)?   ' 
      READ(*,*) npo 
      IF(npo+2.gt.npox)THEN 
         WRITE(*,*)' npo=',npo,' too large for npox=',npox 
         npo=npox-2 
         WRITE(*,*)' used npox-2' 
      ENDIF 
      WRITE(iun20,*)'no. points ',npo 
      RETURN 
    END SUBROUTINE asscbd
! ================================================                      
! START   Estimate of the mean motion and semimajor axis                
!         allowing to combine two arcs                                  
! Input   el0,elp elements of both arcs  
!         MUST BE EQU!!!
!         iorb select orbit 1,2                                         
! Output  ng number of revolutions                                      
!         enm mean motion value for central time                        
!         am  semimajor axis value                                      
!         plm mean longitude value                                      
! ===============INTERFACE========================                      
SUBROUTINE start(el0,elp,iorb,ng,enm,am,plm,ok)
  USE fund_const 
  USE orbit_elements
  implicit none 
! ===========INPUT==================
  TYPE(orbit_elem), INTENT(IN) :: el0,elp 
  integer,intent(in) :: iorb 
! ===========OUTPUT=================                                    
  double precision,intent(out):: enm,am,plm 
  integer, intent(out) :: ng 
  logical, intent(out) :: ok
! ==========END INTERFACE===========  
  double precision eq0(6),eqp(6),t0,tp 
  double precision en1,pl12,pl12p,en2,pl21,pl21p,enm0,enmp 
  double precision tm,pll1,pll2 
  integer ng0,ngp 
! functions                                                             
  DOUBLE PRECISION primea 
! =================================================                     
  IF(el0%coo.ne.'EQU'.or.elp%coo.ne.'EQU')THEN
     WRITE(*,*)'start: not possible unless EQU, given ', el0%coo, elp%coo
     ok=.false.
     RETURN
  ENDIF  
  eq0=el0%coord
  eqp=elp%coord
  t0=el0%t
  tp=elp%t
! Prediction of $\lambda$ at time tp, starting from elements at time t0 
  en1=sqrt(gms/eq0(1)**3) 
  pl12=eq0(6)+en1*(tp-t0) 
  ng0=FLOOR(pl12/dpig) 
  pl12p=pl12-ng0*dpig 
! The number of revolutions is computed to reduce the discrepancies in  
! lambda between the two arcs                                           
  if(eqp(6)-pl12p.gt.pig)then 
     ng0=ng0-1 
  elseif(pl12p-eqp(6).gt.pig)then 
     ng0=ng0+1 
  endif
! Prediction of $\lambda$ at time t0, starting from elements at time tp 
  en2=sqrt(gms/eqp(1)**3) 
  pl21=eqp(6)+en2*(t0-tp) 
  ngp=FLOOR(pl21/dpig) 
  pl21p=pl21-ngp*dpig 
! The number of revolutions is computed to reduce the discrepancies in  
! $\lambda$ \hfil\break                                                 
! between the two arcs                                                  
  if(eq0(6)-pl21p.gt.pig)then 
     ngp=ngp-1 
  elseif(pl21p-eq0(6).gt.pig)then 
     ngp=ngp+1 
  endif
! Computation of mean mean motion using the 2-body model, starting      
! firstly from arc 1 and then from arc 2                                
  enm0=(eqp(6)-eq0(6)+ng0*dpig)/(tp-t0) 
  enmp=(eq0(6)-eqp(6)+ngp*dpig)/(t0-tp) 
! Choose  the best orbit to start from                                  
  write(*, 177)ng0,ngp 
177 format(' number of rev. ',i4, i4) 
  if(iorb.eq.1)then 
     enm=enm0 
     ng=ng0 
  elseif(iorb.eq.2)then 
     enm=enmp 
     ng=-ngp 
  else 
     write(*,*)' start: iorb=',iorb 
     stop 
  endif
  am=(gms/enm**2)**(1.d0/3.d0) 
! 1997 change:use estimated mean motion to compute mean                 
! longitude at time in the middle                                       
  tm=(t0+tp)/2.d0 
  pll1=eq0(6)+enm*(tm-t0) 
  pll2=eqp(6)+enm*(tm-tp) 
  plm=primea(pll1,pll2)
  ok=.true.
END SUBROUTINE start

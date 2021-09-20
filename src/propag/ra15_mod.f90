! ================ MODULE ra15_mod =======
! CONTAINS 
!PUBLIC ROUTINES
!                ra15  Everhart Radau
! PRIVATE ROUTINES
!                rasust
!                rapred
!                rabeta
!                radcon
!                bintrp
! HEADERS, MODULES used
! ra15_mod.o: \
!	../include/nvarx.h90 \
!	../suit/OUTPUT_CONTROL.mod \
!	close_app.o \
!	force_model.o 
!
MODULE ra15_mod
USE output_control
USE fund_const
IMPLICIT NONE
PRIVATE
!PUBLIC ROUTINES
PUBLIC ra15
! PUBLIC DATA
      INTEGER iusci ! output control
      PUBLIC iusci
      LOGICAL deriv_safe ! to get control on variational eqns too
      PUBLIC deriv_safe
! controls for numerical integrator: Everhart
      INTEGER llev,lit1_r,lit2_r,lit1_rc,lit2_rc
      DOUBLE PRECISION hev,eprk_r
      PUBLIC llev,lit1_r,lit2_r,lit1_rc,lit2_rc,hev,eprk_r
!PRIVATE COMMON DATA

      INTEGER nwrite
      INTEGER, PARAMETER :: nwritx=10
      INTEGER, PARAMETER::  itmaxr=20 ! max no iterations
! constants of the integration method 
      DOUBLE PRECISION h(8),w(7),u(7),c(21),d(21),r(21),w1 
CONTAINS
! ========================================================              
! RA15V                                                                 
! Vers 2.29 last update 9 January 2002 A. Milani                        
! with partial recomputation of right hand side.                        
!                                                                       
! structured version; subroutines are in rasubs.f                       
!                                                                       
!  integrator radau by e. everhart, physics department, university of de
!  this 15th-order version, called ra15, is written out for faster execu
!  y'=f(y,t) is  nclass=1,  y"=f(y,t) is nclass= -2,  y"=f(y',y,t) is nc
!  tfin is t(final); tini is t(initial).                                
!  nv = the number of simultaneous differential equations.              
!  the dimensioning below assumes nv will not be larger than 60.        
!  ll controls sequence size. thus ss=10**(-ll) controls the size of a t
!  a typical ll-value is in the range 6 to 12 for this order 15 program.
!  however, if ll.lt.0 then xl is the constant sequence size used.      
!  x and v enter as the starting position-velocity vector, at time tini,
!  and are output as the final position-velocity vector at time tfin.   
! =========================INTERFACE===========================         
      SUBROUTINE ra15(x,v,tini,tfin,tcur,nv,nclass,idcend) 
      USE close_app
      USE force_model
      USE force_sat
      USE force9d
! =================INPUT============================                    
! initial and final time, current time at output                        
! (in case propagation is interrupted)                                  
      DOUBLE PRECISION tini,tfin,tcur 
! dimension of state variables                                          
      INTEGER nv 
! error control                                                         
      INTEGER ll 
! initial stepsize                                                      
      DOUBLE PRECISION xl 
! equation type                                                         
      INTEGER nclass 
! =================OUPUT============================                    
! state variables                                                       
      INCLUDE 'nvarx.h90' 
      DOUBLE PRECISION x(nvar2x),v(nvar2x) 
! close approach end flag                                               
      INTEGER idcend 
! =====================END INTERFACE===========================         
! total number of function calls                                        
      INTEGER nf 
! state variables (WARNING: there is a fixed limit here!)                
      INTEGER, PARAMETER :: ndx=nvar2x,nvx=nvarx
! right hand side at begin/end of step                                  
      DOUBLE PRECISION f1(nvx) 
! memory control flag                                                   
      INTEGER ips,imem 
! storage arrays                                                        
      DOUBLE PRECISION b(7,nvx),g(7,nvx),e(7,nvx),bd(7,nvx) 
! ===============CLOSE APPROACH CONTROL====================             
! no headers: close app. options from force_model.mod 
      INTEGER idc 
! positon and velocities of planet which has close-encounter            
! whith the asteroid                                                    
      DOUBLE PRECISION xpla(6) 
      LOGICAL cloend 
! logical flag for variational eq., state transition matrix accumulation
      LOGICAL variaz 
! ========================================================              
! logical control                                                       
      LOGICAL npq,nsf,nper,ncl,nes 
! ========================================================              
! time and stepsize                                                     
! initial stepsize with sign,time direction, final and current time from
      DOUBLE PRECISION xldir,dir,tf,tm 
! stepsize control: first guess, current stepsize, t or t**2            
      DOUBLE PRECISION tp,t,t2 
! error control, abs. val. stepsize, truncation error parameter         
      DOUBLE PRECISION ss,tval,hv 
      DOUBLE PRECISION ep(itmaxr) 
! intermediate times                                                    
      DOUBLE PRECISION s,q 
! scalar temporaries                                                    
      DOUBLE PRECISION g_k,temp 
! loop indexes k=1,nv; l=1,7; j=2,8                                     
      INTEGER k,l,j 
! iteration count for implicit equations m=1,ni                         
      INTEGER m 
! iteration number control                                              
      INTEGER ni 
! count of steps done                                                   
      INTEGER ns 
! count of step reductions                                              
      INTEGER ncount 
! ========================================================              
! scalar variables                                                      
! constants                                                             
      DOUBLE PRECISION pw 
! static memory allocation                                              
      SAVE 
! ===============to remove ORBFIT commons========================       
!     INTEGER lit1,lit2,itmax,iusci,ipirip                              
!     DOUBLE PRECISION eprk                                             
!     parameter (itmax=20)                                              
!     lit1=10                                                           
!     lit2=4                                                            
!     ipirip=10                                                         
!     iusci=0                                                           
!     eprk=1.d-12                                                       
!     iclap=0                                                           
! ========================================================              
! begin execution                                                       
! ========================================================              
! control against excessive write
      nwrite=0
! setup of logical controls                                             
!  y'=f(y,t)  ncl=.true.    y"=f(y,t)  ncl=.false.   y"=f(y',y,t) ncl=.f
!  nclass=1   npq=.true.    nclass= -2 npq=.true.    nclass= 2    npq=.f
      npq=nclass.lt.2 
      ncl=nclass.eq.1 
! set to zero storage arrays (for first order equations only)           
      DO  k=1,nv 
        IF(ncl) v(k)=0.d0 
      ENDDO 
!  nper is .true. only on last sequence of the integration.             
      nper=.false. 
!  nsf is .false. on starting sequence, otherwise .true.                
      nsf=.false. 
!  nes is .true. only if ll is negative. then the sequence size is xl.  
      ll=llev 
      xl=hev 
      nes=ll.lt.0 
! variaz is true if the state transition matrix is being propagated;    
      variaz=nv.gt.3 
! ===============================================================       
!  evaluate the constants in the h-, w-, u-, c-, d-, and r-vectors      
! ===============================================================       
      CALL radcon(ncl) 
      pw=1.d0/9.d0 
! ==========================================================            
! initialisation                                                        
! ==========================================================            
! direction of time                                                     
      dir=1.d0 
      tf=tfin-tini 
      if(tf.lt.0.d0) dir=-1.d0 
      xldir=dir*abs(xl) 
      ss=10.**(-ll) 
! ===============================================================       
!  the statements above are used only once in an integration to set up t
!  constants.  the  next set in                                         
!  a reasonable estimate to tp based on experience. same sign as dir.   
!  an initial first sequence size can be set with xl even with ll positi
! ===============================================================       
      IF(xldir.ne.0.d0)then 
        tp=xldir 
      ELSE 
         IF(nes)then 
            write(*,*)' ra15v: fixed stepsize xl=0; ll=',ll 
            stop 
         ELSE 
            tp=0.1d0*dir 
         ENDIF 
      ENDIF 
      IF(tp/tf.gt.1.d0)then 
         tp=tf 
      ELSEIF(tp/tf.gt.0.5d0)then 
         tp=0.5d0*tf 
      ENDIF 
      ncount=0 
! ============================================================          
! information on initial state of the integrator                        
      IF(iusci.gt.100)write(ipirip,999)tini,tfin,tp*dir 
  999 format(' ra15: from tini=',f12.4,' to tfin=',f12.4,               &
     &   ' with max. step h=',1p,d12.4)                                 
! ===============================================================       
!  line 4000 is the starting place of the first sequence.               
! ===============================================================       
 4000 ns=0 
      nf=0 
! ===============================================================       
! force initial value for first step                                    
! ================================================================      
      tm=0.d0 
      ips=0 
      imem=1 
      IF(rhs.eq.1)THEN
         CALL force(x, v, tm+tini, f1,nv,idc,xpla,ips,imem)
      ELSEIF(rhs.eq.2)THEN
         CALL forcesat(x, v, tm+tini, f1,nv,idc,xpla,ips,imem)
      ELSEIF(rhs.eq.3)THEN
         CALL force9(x,v,tm+tini,f1,nv,idc,xpla,ips,imem) 
      ENDIF 
      nf=nf+1 
!                                                                       
 3000 CONTINUE 
! set to zero storage arrays                                            
! when extrapolation from previous step is not to be used               
      do  k=1,nv 
        do  l=1,7 
          e(l,k)=0.d0 
          b(l,k)=0.d0 
          bd(l,k)=0.d0 
        enddo 
      enddo 
! ===============================================================       
! line 722  begins every step after the first.                          
  722 CONTINUE 
! ===============================================================       
! compute array g from b                                                
      CALL rabeta(nv,b,d,g) 
! ===============================================================       
! set time variables                                                    
      t=tp 
      t2=t*t 
      IF(ncl) t2=t 
      tval=abs(t) 
! ===============================================================       
!  loop 175 is lit1_r iterations on first step and lit2_r iterations the
      IF(nes)THEN 
! problem: if the use of fixed step size is forced by setting in input l
! (in the option file, propag.llev ,0), then the number of iterations   
! used for close approaches extends to the entire integration.
         IF(nsf)THEN 
            ni=lit1_rc 
         ELSE 
            ni=lit2_rc 
         ENDIF 
      ELSE 
         IF(nsf)THEN 
            ni=lit2_r 
         ELSE 
            ni=lit1_r 
         ENDIF 
      ENDIF 
      do 175 m=1,ni 
! =========do on substep=========================                       
        CALL rasust(m,t,t2,tm,tini,x,v,b,f1,nv,ncl,npq,g,ep(m),nf) 
! =====================================================                 
!  iteration of step is over.                                           
! =====================================================                 
! iteration control based on norm of differences                        
        IF(m.eq.1)goto 175 
        IF(ep(m)/ep(1).lt.eprk_r)then 
! prepare for stepsize determination of next step                       
           IF(.not.nes)then 
! if variable stepsize, step size control                               
              hv=0.d0 
              do 635 k=1,nv 
  635           hv=dmax1(hv,abs(b(7,k))) 
              hv=hv*w(7)/tval**7 
           ENDIF 
           IF(verb_pro.gt.39)write(ipirip,996)ni,tini+tm+t,ep(1)         &
     &            ,(ep(j)/ep(1),j=2,m)                                  
  996      format('ra15v: good convergence iter ',i3,' time=',f10.2,    &
     &    ' controls:'/(5d12.4/))                                       
           goto 176 
        ENDIF 
! =========end do on iterations=========================                
  175 continue 
! ======================================================                
! bad convergence of iterations                                         
      IF(verb_pro.gt.19)write(ipirip,997)ni,eprk_r,ep(1)                      &
     &            ,(ep(j)/ep(1),j=2,ni)                                 
  997 format('ra15v: bad convergence iter ',i3,' eprk_r=',d10.2,          &
     &    ' controls:'/(5d12.4/))                                       
      IF(nes)then 
         IF(nwrite.lt.2*nwritx)THEN
            write(*,*)tm,ep 
            write(*,*)' ra15v: non convergence with fixed step ',t 
         ENDIF                                       
      ELSE 
         tp=.8d0*tp 
         goto 3000 
      ENDIF 
! ======================================================                
! satisfactory convergence; continue as in old version                  
  176 IF (.not.nsf)then 
! ======================================================                
! block executed only for first step                                    
! ======================================================                
         IF(nes)then 
! fixed stepsize                                                        
            tp=xldir 
         ELSE 
! variable stepsize: compute initial stepsize                           
            tp=(ss/hv)**pw*dir 
! not too quick increase                                                
            IF(tp/t.gt.1.4d0) tp=t*1.4d0 
! not too quick decrease                                                
            IF(tp/t.lt.1.d0)then 
               tp=.8d0*tp 
               ncount=ncount+1 
               IF(ncount.gt.20) then 
! too many shortenings of the stepsize                                  
                  write(ipirip,*)' ra15v: ncount.gt.10',xl,tp 
                  stop 
               ENDIF 
        IF(ncount.gt.1.and.iusci.gt.100)write(ipirip,888)ncount,t,tp 
  888          format (2x,2i3,2d18.10) 
!  restart with 0.8x sequence size if new size called for is smaller tha
!  originally chosen starting sequence size on first sequence.          
               go to 4000 
            ELSE 
!  initial step accepted                                                
               IF(iusci.gt.100)write(ipirip,*)'ra15v: second step=',tp 
            ENDIF 
         ENDIF 
         nsf=.true. 
! ======================================================                
         IF(iusci.gt.100)write(ipirip,*)' tp ',tp 
      ENDIF 
! =====================================================                 
      CALL rapred(ncl,x,v,t,t2,f1,b,nv) 
      DO k=1,nv 
        IF((abs(x(k)).gt.1.d15.or.abs(v(k)).gt.1.d15).and.nwrite.lt.nwritx)THEN 
            WRITE(ierrou,*)' rapred ',k,x(k),v(k),f1(k) 
            nwrite=nwrite+1
            numerr=numerr+1
        ENDIF 
      ENDDO 
! =====================================================                 
! current time update                                                   
      tm=tm+t 
      ns=ns+1 
!  return if done.                                                      
      IF(nper.or.abs(t-tf).lt.1.d-10)then 
         tcur=tfin 
         IF(iusci.gt.100)write(ipirip,*) 'nf,ns,tfin ',nf,ns,tfin 
         return 
      ENDIF 
!  control on size of next sequence and adjust last sequence to exactly 
!  cover the integration span. nper=.true. set on last sequence.        
! ===============================================================       
! force initial value for next step                                     
! ================================================================      
      idcend=idc 
      ips=0 
      imem=1 
      IF(rhs.eq.1)THEN
         CALL force(x,v,tm+tini,f1,nv,idc,xpla,ips,imem)
      ELSEIF(rhs.eq.2)THEN
         CALL forcesat(x,v,tm+tini,f1,nv,idc,xpla,ips,imem)
      ELSEIF(rhs.eq.3)THEN
         CALL force9(x,v,tm+tini,f1,nv,idc,xpla,ips,imem)
      ENDIF
      nf=nf+1 
! ================================================================      
! close approach control                                                
      IF(iclap.eq.1)THEN 
! problem: if the integration begins inside a clsoe approach, then the f
! is done with variable stepsize; will the stepsize adjust automatically
! to a value clsoe to the true-anomaly controlled one???   
         IF(rhs.eq.1)THEN             
            CALL cloapp(tm+tini,t,x,v,nv,idc,xpla,xldir,dir,nes,cloend) 
         ELSEIF(rhs.eq.2)THEN
! ???????????????????????????
            STOP
         ELSEIF(rhs.eq.3)THEN
            CALL cloapp9(idc,tm+tini)
         ENDIF
         IF(cloend)THEN
            tcur=tm+tini 
            return 
         ELSE 
!            idcend=0                                                   
         ENDIF 
      ENDIF 
! ================================================================      
! check if the integration required is finished                         
      IF(nes)THEN 
         tp=xldir 
      ELSE 
         tp=dir*(ss/hv)**pw 
         IF(tp/t.gt.1.4d0) tp=t*1.4d0 
!         WRITE(*,*)' next step ',tp                                    
         IF(iusci.gt.100)write(ipirip,*)'step=',tp 
      ENDIF 
      IF(dir*(tm+tp).gt.dir*tf-1.d-8)THEN 
         tp=tf-tm 
         nper=.true. 
      ENDIF 
! ================================================================      
!  new step, no extrapolation, recompute right hand side                
      IF(.not.nsf) GOTO 4000 
!  if the integration continues without discontinuities,                
!  extrapolate b-values for next step.                                  
      q=tp/t 
      CALL bintrp(q,b,e,bd,nv,ns) 
!  new step                                                             
!     IF(.not.nsf) GOTO 3000                                            
      go to 722 
      END SUBROUTINE ra15                                          
! ================================================                      
! rasust                                                                
!                                                                       
! iterative convergence of the right hand sides at intermediate points  
! ================================================                      
      SUBROUTINE rasust(m,t,t2,tm,tini,x,v,b,f1,nv,ncl,npq,g,epsi,nf)
      USE force_model                                               
      USE force_sat
      USE force9d
      INTEGER,INTENT(IN) :: m ! iteration number 
      DOUBLE PRECISION epsi ! control to be used for convergence 
      INTEGER,INTENT(IN) :: nv ! dimension of state vector 
! logical flag for first order diff. eq., for second order diff.eq.     
      LOGICAL, INTENT(IN) :: ncl,npq 
      DOUBLE PRECISION x(nv),v(nv) ! state vector, pos and vel 
      DOUBLE PRECISION f1(nv) ! right hand side at begining of step 
      DOUBLE PRECISION b(7,nv),g(7,nv) ! storage arrays 
      DOUBLE PRECISION t,t2,tm,tini ! times
! =====================END INTERFACE===========================         
! memory control                                                        
      INTEGER ips,imem 
! workspace state variables
      INTEGER ndx,nvx 
      PARAMETER (ndx=1221,nvx=2*ndx) 
! x,v temporary                                                         
      DOUBLE PRECISION y(nvx),z(nvx) 
      DOUBLE PRECISION ck(nvx,8),fj(nvx) 
! total number of function calls, number components for convergence test
      INTEGER nf, nvv 
! loop indexes                                                          
      INTEGER j,k,jd 
! scalar temporaries                                                    
      DOUBLE PRECISION q,s,temp,g_k 
! close approach flags: detected, which planet                          
      INTEGER idc 
! positon and velocities of planet which has close-encounter            
! whith the asteroid                                                    
      DOUBLE PRECISION xpla(6) 
! memory model static                                                   
      SAVE 
! initialize previous iteration array to zero if first iteration        
      IF(m.eq.1)THEN 
         DO  k=1,nv 
            DO j=2,8 
               ck(k,j)=0.d0 
            ENDDO 
         ENDDO 
      ENDIF 
      epsi=0.d0 
! =========do on substep=========================                       
!  loop 174 is for each substep within a sequence.                      
      do 174 j=2,8 

!         jdm=j-2                                                       
          s=h(j) 
          q=s 
          IF(ncl) q=1.d0 
! ===============================================================       
!  use eqs. (2.9) and (2.10) of text to predict positions at each subste
!  these collapsed series are broken into two parts because an otherwise
!  excellent  compiler could not handle the complicated expression.     
! ===============================================================       
          do 130 k=1,nv 
            temp=w(3)*b(3,k)+s*(w(4)*b(4,k)+s*(w(5)*b(5,k)              &
     &        +s*(w(6)*b(6,k)+s*w(7)*b(7,k))))                          
            y(k)=x(k)+q*(t*v(k)+t2*s*(f1(k)*w1+s*(w(1)*b(1,k)           &
     &          +s*(w(2)*b(2,k)+s*temp))))                              
            IF(.not.npq)then 
!  next are calculated the velocity predictors need for general class ii
               temp=u(3)*b(3,k)+s*(u(4)*b(4,k)+s*(u(5)*b(5,k)           &
     &              +s*(u(6)*b(6,k)+s*u(7)*b(7,k))))                    
               z(k)=v(k)+s*t*(f1(k)+s*(u(1)*b(1,k)+s*(u(2)*b(2,k)       &
     &             +s*temp)))                                           
            ENDIF 
  130     continue 
! ==================================================================    
!  find forces at each substep.                                         
! ==================================================================    
          IF(m.eq.1)THEN 
             ips=0 
          ELSE 
             ips=-1 
          ENDIF 
          imem=j 
          IF(rhs.eq.1)THEN
             CALL force(y,z,tm+s*t+tini,fj,nv,idc,xpla,ips,imem)
          ELSEIF(rhs.eq.2)THEN
             CALL forcesat(y,z,tm+s*t+tini,fj,nv,idc,xpla,ips,imem)
          ELSEIF(rhs.eq.3)THEN
             CALL force9(y,z,tm+s*t+tini,fj,nv,idc,xpla,ips,imem)
          ENDIF
          IF(deriv_safe)THEN
             nvv=nv
          ELSE
             nvv=3
          ENDIF
          do k=1,nvv
            epsi=epsi+dabs(fj(k)-ck(k,j)) 
          enddo 
          do  k=1,nv 
            ck(k,j)=fj(k) 
          enddo 
          nf=nf+1 
! ================do on components===================                   
          do 171 k=1,nv 
! ===============================================================       
!  find g-value for the force fj found at the current substep. this     
!  section, including the many-branched goto, uses eq. (2.4) of text.   
! ===============================================================       
!		write(*,*)j-1,k
            temp=g(j-1,k) 
            g_k=(fj(k)-f1(k))/s 
            IF(j.le.2)then 
               g(1,k)=g_k 
            ELSEIF(j.eq.3)then 
                g(2,k)=(g_k-g(1,k))*r(1) 
            ELSEIF(j.eq.4)then 
                g(3,k)=((g_k-g(1,k))*r(2)-g(2,k))*r(3) 
            ELSEIF(j.eq.5)then 
                g(4,k)=(((g_k-g(1,k))*r(4)-g(2,k))*r(5)-g(3,k))*r(6) 
            ELSEIF(j.eq.6)then 
                g(5,k)=((((g_k-g(1,k))*r(7)-g(2,k))*r(8)-g(3,k))*r(9)-   &
     &                    g(4,k))*r(10)                                 
            ELSEIF(j.eq.7)then 
                g(6,k)=(((((g_k-g(1,k))*r(11)-g(2,k))*r(12)              &
     &                -g(3,k))*r(13)-g(4,k))*r(14)-g(5,k))*r(15)        
            ELSEIF(j.eq.8)then 
                g(7,k)=((((((g_k-g(1,k))*r(16)-g(2,k))*r(17)             &
     &                -g(3,k))*r(18)-g(4,k))*r(19)                      &
     &                -g(5,k))*r(20)-g(6,k))*r(21)                      
            ENDIF 
! ===============================================================       
!  upgrade all b-values                                                 
! ===============================================================       
            temp=g(j-1,k)-temp 
            b(j-1,k)=b(j-1,k)+temp 
! ===============================================================       
!  temp is now the improvement on g(jd,k) over its former value.        
!  now we upgrade the b-value using this dfference in the one term.     
!  this section is based on eq. (2.5).                                  
! ===============================================================       
            IF(j.eq.3)then 
                b(1,k)=b(1,k)+c(1)*temp 
            ELSEIF(j.eq.4)then 
                b(1,k)=b(1,k)+c(2)*temp 
                b(2,k)=b(2,k)+c(3)*temp 
            ELSEIF(j.eq.5)then 
                b(1,k)=b(1,k)+c(4)*temp 
                b(2,k)=b(2,k)+c(5)*temp 
                b(3,k)=b(3,k)+c(6)*temp 
            ELSEIF(j.eq.6)then 
                b(1,k)=b(1,k)+c(7)*temp 
                b(2,k)=b(2,k)+c(8)*temp 
                b(3,k)=b(3,k)+c(9)*temp 
                b(4,k)=b(4,k)+c(10)*temp 
            ELSEIF(j.eq.7)then 
                b(1,k)=b(1,k)+c(11)*temp 
                b(2,k)=b(2,k)+c(12)*temp 
                b(3,k)=b(3,k)+c(13)*temp 
                b(4,k)=b(4,k)+c(14)*temp 
                b(5,k)=b(5,k)+c(15)*temp 
            ELSEIF(j.eq.8)then 
                b(1,k)=b(1,k)+c(16)*temp 
                b(2,k)=b(2,k)+c(17)*temp 
                b(3,k)=b(3,k)+c(18)*temp 
                b(4,k)=b(4,k)+c(19)*temp 
                b(5,k)=b(5,k)+c(20)*temp 
                b(6,k)=b(6,k)+c(21)*temp 
            ENDIF 
! =========end do on components=========================                
  171     continue 
! =========end do on substeps=========================                  
  174   continue 
        RETURN 
      END SUBROUTINE rasust                                          
                                                                        
! =====================================================                 
! rapred                                                                
!                                                                       
! finds new x and v values at end of sequence using eqs. (2.11),(2.12)  
! =====================================================                 
      SUBROUTINE rapred(ncl,x,v,t,t2,f1,b,nv) 
      LOGICAL  ncl ! flag for first order equations
      INTEGER nv ! actual dimension                    
      DOUBLE PRECISION,INTENT(IN) :: t,t2 ! stepsize, squared for 2nd order eq.
      DOUBLE PRECISION x(nv),v(nv) ! state variables         
      DOUBLE PRECISION, INTENT(IN) :: f1(nv) ! right hand side, begin step
      DOUBLE PRECISION, INTENT(IN) :: b(7,nv) ! storage arrays 
! end interface                                                         
      INTEGER k ! loop indexes 
! ============================================                          
      DO k=1,nv 
        x(k)=x(k)+v(k)*t+t2*(f1(k)*w1+b(1,k)*w(1)+b(2,k)*w(2)           &
     &   +b(3,k)*w(3)+b(4,k)*w(4)+b(5,k)*w(5)+b(6,k)*w(6)+b(7,k)*w(7))  
        IF(.not.ncl)THEN 
           v(k)=v(k)+t*(f1(k)+b(1,k)*u(1)+b(2,k)*u(2)+b(3,k)*u(3)       &
     &       +b(4,k)*u(4)+b(5,k)*u(5)+b(6,k)*u(6)+b(7,k)*u(7))          
        ENDIF 
      ENDDO 
      RETURN 
      END SUBROUTINE rapred                                          
! ===============================================================       
! rabeta                                                                
!                                                                       
!  find new beta-values from the predicted b-values,                    
!  following eq. (2.7) in text.                                         
! ===============================================================       
      SUBROUTINE rabeta(nv,b,d,g) 
! number of equations                                                   
      INTEGER nv 
! storage arrays                                                        
      DOUBLE PRECISION b(7,nv),g(7,nv),d(21) 
! loop indexes                                                          
      INTEGER k 
! =========================                                             
      DO k=1,nv 
        g(1,k)=b(1,k)+d(1)*b(2,k)+d(2)*b(3,k)+                          &
     &             d(4)*b(4,k)+d( 7)*b(5,k)+d(11)*b(6,k)+d(16)*b(7,k)   
        g(2,k)=            b(2,k)+d(3)*b(3,k)+                          &
     &             d(5)*b(4,k)+d( 8)*b(5,k)+d(12)*b(6,k)+d(17)*b(7,k)   
        g(3,k)=            b(3,k)+                                      &
     &             d(6)*b(4,k)+d( 9)*b(5,k)+d(13)*b(6,k)+d(18)*b(7,k)   
        g(4,k)=         b(4,k)+d(10)*b(5,k)+d(14)*b(6,k)+d(19)*b(7,k) 
        g(5,k)=                      b(5,k)+d(15)*b(6,k)+d(20)*b(7,k) 
        g(6,k)=                                   b(6,k)+d(21)*b(7,k) 
        g(7,k)=                                                b(7,k) 
      ENDDO 
      RETURN 
      END SUBROUTINE rabeta                                          
! ===============================================================       
!  bintrp                                                               
!                                                                       
!  now predict b-values for next step. the predicted values from the pre
!  sequence were saved in the e-matrix. te correction bd between the act
!  b-values found and these predicted values is applied in advance to th
!  next sequence. the gain in accuracy is significant. using eqs. (2.13)
! ===============================================================       
      SUBROUTINE bintrp(q,b,e,bd,nv,ns) 
! =============INPUT===============                                     
! rescaled time                                                         
      DOUBLE PRECISION q 
! dimension of state vector, number of steps completed                  
      INTEGER nv,ns 
! =============INPUT AND OUTPUT=============                            
! storage arrays                                                        
      DOUBLE PRECISION e(7,nv),b(7,nv),bd(7,nv) 
! =============END INTERFACE================                            
      INTEGER k,j,l 
! ===============================================================       
      DO 39 k=1,nv 
        IF(ns.ne.1)then 
           DO  j=1,7 
             bd(j,k)=b(j,k)-e(j,k) 
           ENDDO 
                                                                        
        ENDIF 
        e(1,k)=      q*(b(1,k)+ 2.d0*b(2,k)+ 3.d0*b(3,k)+               &
     &           4.d0*b(4,k)+ 5.d0*b(5,k)+ 6.d0*b(6,k)+ 7.d0*b(7,k))    
        e(2,k)=                q**2*(b(2,k)+ 3.d0*b(3,k)+               &
     &           6.d0*b(4,k)+10.d0*b(5,k)+15.d0*b(6,k)+21.d0*b(7,k))    
        e(3,k)=                             q**3*(b(3,k)+               &
     &           4.d0*b(4,k)+10.d0*b(5,k)+20.d0*b(6,k)+35.d0*b(7,k))    
        e(4,k)= q**4*(b(4,k)+ 5.d0*b(5,k)+15.d0*b(6,k)+35.d0*b(7,k)) 
        e(5,k)=              q**5*(b(5,k)+ 6.d0*b(6,k)+21.d0*b(7,k)) 
        e(6,k)=                           q**6*(b(6,k)+ 7.d0*b(7,k)) 
        e(7,k)=                                         q**7*b(7,k) 
        DO  l=1,7 
          b(l,k)=e(l,k)+bd(l,k) 
        ENDDO 
   39 ENDDO 
      RETURN 
      END SUBROUTINE bintrp                                          
! ===============================================================       
!  RADCON assigns h spacings and                                        
!  evaluates the constants in the w-, u-, c-, d-, and r-vectors         
! ===========begin interface=========================================   
      SUBROUTINE radcon(ncl) 
      LOGICAL ncl 
! ===========end interface===========================================   
! loop indexes                                                          
      INTEGER n,k,la,lc,lb,ld,le,l 
! scalar real                                                           
      DOUBLE PRECISION ww 
      DOUBLE PRECISION half,one 
! integers                                                              
      INTEGER nw(8) 
      DATA half, one/0.5d0, 1.0d0/ 
      DATA nw/0,0,1,3,6,10,15,21/ 
! ===================================================================   
!  these h values are the gauss-radau spacings, scaled to the range 0 to
!  for integrating to order 15.                                         
      h(1)=0.d0 
      h(2)= .05626256053692215d0 
      h(3)= .18024069173689236d0 
      h(4)=.35262471711316964d0 
      h(5)=.54715362633055538d0 
      h(6)= .73421017721541053d0 
      h(7)=.88532094683909577d0 
      h(8)= .97752061356128750d0 
!  the sum of the h-values should be 3.73333333333333333                
! ====================================================================  
!  evaluate the constants in the w-, u-, c-, d-, and r-vectors          
      do 14 n=2,8 
        ww=n+n*n 
        IF(ncl) ww=n 
        w(n-1)=one/ww 
        ww=n 
   14   u(n-1)=one/ww 
      w1=half 
      IF(ncl) w1=one 
      c(1)=-h(2) 
      d(1)=h(2) 
      r(1)=one/(h(3)-h(2)) 
      la=1 
      lc=1 
      do 73 k=3,7 
        lb=la 
        la=lc+1 
        lc=nw(k+1) 
        c(la)=-h(k)*c(lb) 
        c(lc)=c(la-1)-h(k) 
        d(la)=h(2)*d(lb) 
        d(lc)=-c(lc) 
        r(la)=one/(h(k+1)-h(2)) 
        r(lc)=one/(h(k+1)-h(k)) 
        IF(k.eq.3) go to 73 
        do 72 l=4,k 
          ld=la+l-3 
          le=lb+l-4 
          c(ld)=c(le)-h(k)*c(le+1) 
          d(ld)=d(le)+h(l-1)*d(le+1) 
   72     r(ld)=one/(h(k+1)-h(l-1)) 
   73 continue 
      return 
      END SUBROUTINE radcon            
!=============================== 
! CLOAPP9
! close approach control for orbit9
 SUBROUTINE cloapp9(idc,tt)
   USE massmod
   INTEGER, INTENT(INOUT) :: idc
   DOUBLE PRECISION, INTENT(IN) :: tt
   INCLUDE 'comnbo.h90'
   INTEGER np,n
   np=MOD(idc,nbod-1)
   IF(np.eq.0) np=nbod-1
   n=(idc-np)/(nbod-1)
   if(idc.ne.0)then
      write(iuncla,*)'t =',tt,' close approach of ',ida(n),' to planet ',nompla(np)
      idc=0
   endif
 END SUBROUTINE cloapp9

END MODULE ra15_mod

! Copyright (C) 2000 by Mario Carpino (carpino@brera.mi.astro.it)       
! Version: May 25, 2000                                                 
! --------------------------------------------------------------------- 
!                                                                       
!  *****************************************************************    
!  *                                                               *    
!  *                         M O D G E Q                           *    
!  *                                                               *    
!  *            Modified solution of Gauss' equation               *    
!  *                                                               *    
!  *****************************************************************    
!                                                                       
! Solution of Gauss' equation for finding v(t1) from r(t1) and r(t2).   
! The method uses closed form formulae which are valid for elliptic,    
! parabolic and hyperbolic orbits                                       
!                                                                       
! INPUT:    GK        -  Gauss' k constant (possibly multiplied by      
!                        the square root of the sum of the masses)      
!           X1        -  Position vector at time t1                     
!           X2        -  Position vector at time t2                     
!           TAU       -  Normalized time difference ( = k*(t2-t1) )     
!           ELEM      -  Orbital elements at time t1                    
!           ELTYPE    -  Type of orbital elements ('KEP'/'COM')         
!                                                                       
! OUTPUT:   FSER      -  F series                                       
!           GSER      -  G series                                       
!           VEL1      -  Orbital velocity at time t1                    
!                                                                       
      SUBROUTINE modgeq(gk,x1,x2,tau,elem,eltype,fser,gser,vel1) 
      IMPLICIT NONE 
                                                                        
      DOUBLE PRECISION gk,x1(3),x2(3),tau,elem(6),fser,gser,vel1(3) 
      CHARACTER*(*) eltype 
                                                                        
      INTEGER i 
      DOUBLE PRECISION ae1,ae2,emme2,ndk,p,dcap,tau0,esse,u3 
      DOUBLE PRECISION sinthe,costhe,theta,sinf,cosf,th2,ecth2,sma 
      DOUBLE PRECISION beta2,beta,b,t0,r,coshf,sinhf,expx,expmx,effe 
      DOUBLE PRECISION emme,anec 
                                                                        
      DOUBLE PRECISION anecc,vsize,hypan 
      EXTERNAL anecc,vsize,hypan 
                                                                        
      IF(eltype.EQ.'KEP') THEN 
          ae1=anecc(elem(6),elem(2)) 
          ndk=elem(1)**(-1.5D0) 
          emme2=elem(6)+tau*ndk 
          ae2=anecc(emme2,elem(2)) 
          fser=(1-COS(ae2-ae1))/(1-elem(2)*COS(ae1)) 
          fser=1-fser 
          gser=ae2-ae1-SIN(ae2-ae1) 
          gser=(tau-gser/ndk)/gk 
      ELSEIF(eltype.EQ.'COM') THEN 
          u3=1.d0/3.d0 
! Parabolic orbit                                                       
          IF(elem(2).EQ.1.D0) THEN 
              p=2.d0*elem(1) 
              ndk=SQRT(1.D0/p**3) 
              dcap=TAN(elem(6)/2.d0) 
              tau0=-dcap*(1.d0+dcap**2/3.d0)/(2.d0*ndk) 
              esse=3.d0*ndk*tau 
              esse=ATAN(1.d0/esse) 
              sinthe=SIN(esse/2.d0) 
              costhe=COS(esse/2.d0) 
              IF(sinthe.GT.0.d0) THEN 
                  sinthe=sinthe**u3 
              ELSEIF(sinthe.LT.0.d0) THEN 
                  sinthe=-(ABS(sinthe)**u3) 
              END IF 
              IF(costhe.GT.0.d0) THEN 
                  costhe=costhe**u3 
              ELSEIF(costhe.LT.0.d0) THEN 
                  costhe=-(ABS(costhe)**u3) 
              END IF 
              theta=ATAN2(sinthe,costhe) 
              sinf=2.d0*COS(2.d0*theta) 
              cosf=SIN(2.d0*theta) 
              th2=2.d0*ATAN2(sinf,cosf) 
! Hyperbolic orbit                                                      
          ELSE 
              p=(1+elem(2))*elem(1) 
              sma=elem(1)/(elem(2)-1) 
              beta2=elem(2)**2-1 
              beta=SQRT(beta2) 
              ndk=SQRT(1/sma**3) 
              b=SQRT(sma*p) 
              IF(elem(6).EQ.0.d0) THEN 
                  t0=0.d0 
              ELSE 
                  cosf=COS(elem(6)) 
                  sinf=SIN(elem(6)) 
                  r=p/(1.d0+elem(2)*cosf) 
                  coshf=elem(2)-r*cosf/sma 
                  sinhf=r*sinf/b 
                  expx=coshf+sinhf 
                  expmx=coshf-sinhf 
                  IF(expx.GT.expmx) THEN 
                      effe=LOG(expx) 
                  ELSE 
                      effe=-LOG(expmx) 
                  END IF 
                  t0=elem(2)*SINH(effe)-effe 
                  t0=-t0/ndk 
              END IF 
              emme=ndk*(tau-t0) 
! Hyperbolic eccentric anomaly                                          
              anec=hypan(emme,elem(2)) 
              th2=SQRT((elem(2)+1)/(elem(2)-1))*TANH(anec/2) 
              th2=2*ATAN(th2) 
          END IF 
          ecth2=elem(2)*COS(th2) 
          fser=(COS(th2-elem(6))+ecth2)/(1+ecth2) 
          gser=vsize(x1)*vsize(x2)*SIN(th2-elem(6)) 
      ELSE 
          STOP '**** modgeq: ELTYPE = ??? ****' 
      END IF 
                                                                        
      DO 1 i=1,3 
      vel1(i)=(x2(i)-fser*x1(i))/gser 
    1 END DO 
                                                                        
      END                                           

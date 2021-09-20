! Copyright (C) 2000 by Mario Carpino (carpino@brera.mi.astro.it)       
! Version: May 29, 2000                                                 
! --------------------------------------------------------------------- 
!                                                                       
!  *****************************************************************    
!  *                                                               *    
!  *                          C C E K 1                            *    
!  *                                                               *    
!  *     Transformation from cartesian coordinates (pos/vel)       *    
!  *                    to keplerian elements                      *    
!  *                                                               *    
!  *****************************************************************    
!                                                                       
!     see: A.E.Roy, Orbital Motion, Adam Hilger Ltd, Chapter 4          
!                                                                       
! INPUT:    XV(6)     -  Position/velocity vector                       
!           GM        -  G * ( M + m )                                  
!                                                                       
! OUTPUT:   ELEM(6)   -  Orbital elements                               
!                            ELEM(1) = a/q                              
!                            ELEM(2) = e                                
!                            ELEM(3) = i                                
!                            ELEM(4) = long. node                       
!                            ELEM(5) = arg.pericenter                   
!                            ELEM(6) = M/f                              
!           TYPE      -  Type of orbital elements (KEP/COM)             
!                                                                       
      SUBROUTINE ccek1(elem,type,xv,gm) 
      USE fund_const
      IMPLICIT NONE 
                                                                        
      DOUBLE PRECISION elem(6),xv(6),gm 
      CHARACTER*(*) type 
                                                                        
      DOUBLE PRECISION eps 
      PARAMETER (eps=5.d-15) 
                                                                        
      INTEGER j 
      DOUBLE PRECISION elle(3),elv(3),xorb(3),vorb(3),el2,elmod,rv,rs 
      DOUBLE PRECISION rot(3,3),r1(3,3),r2(3,3),sini,ainc,anod,v2 
      DOUBLE PRECISION reca,sma,enne,esine,ecose,ecc,anec,emme,x1,x2 
      DOUBLE PRECISION rad,xm,sinper,cosper,argper,p,q,cosf,sinf,effe 
      DOUBLE PRECISION ecosf,esinf 
                                                                        
! Angular momentum vector                                               
      elle(1)=xv(2)*xv(6)-xv(3)*xv(5) 
      elle(2)=xv(3)*xv(4)-xv(1)*xv(6) 
      elle(3)=xv(1)*xv(5)-xv(2)*xv(4) 
                                                                        
! Angular momentum unit vector                                          
      el2=elle(1)**2+elle(2)**2+elle(3)**2 
      elmod=SQRT(el2) 
      DO 1 j=1,3 
      elv(j)=elle(j)/elmod 
    1 END DO 
                                                                        
! Orbital inclination and longitude of the node                         
      sini=SQRT(elv(1)**2+elv(2)**2) 
      ainc=ATAN2(sini,elv(3)) 
      IF(ainc.GT.dpig) ainc=ainc-dpig 
      IF(ainc.LT.0.d0) ainc=ainc+dpig 
      IF(sini.EQ.0.d0) THEN 
          ainc=0.d0 
          anod=0.d0 
      ELSE 
          anod=ATAN2(elv(1),-elv(2)) 
          IF(anod.GT.dpig) anod=anod-dpig 
          IF(anod.LT.0.d0) anod=anod+dpig 
      END IF 
                                                                        
! Cartesian coordinates with respect to the "orbital frame" (with X     
! axis directed along the line of nodes)                                
      CALL rotmt(ainc,r1,1) 
      CALL rotmt(anod,r2,3) 
      CALL prodmm(rot,r1,r2) 
      CALL prodmv(xorb,rot,xv(1)) 
      CALL prodmv(vorb,rot,xv(4)) 
      rv=xorb(1)*vorb(1)+xorb(2)*vorb(2) 
                                                                        
! Reciprocal semimajor axis                                             
      rs=SQRT(xorb(1)**2+xorb(2)**2) 
      v2=vorb(1)**2+vorb(2)**2 
      reca=2.d0/rs-v2/gm 
    2 CONTINUE 
                                                                        
! CASE 1: RECA > 0 (elliptic orbit)                                     
      IF(reca.GT.0.d0) THEN 
          type='KEP' 
          sma=1.d0/reca 
          enne=SQRT(gm/(sma**3)) 
                                                                        
! Eccentricity                                                          
          esine=rv/(enne*sma**2) 
          ecose=v2*rs/gm-1.d0 
          ecc=SQRT(esine**2+ecose**2) 
          IF(ABS(ecc-1.d0).LT.eps) THEN 
              reca=0.d0 
              GOTO 2 
          END IF 
                                                                        
! Eccentric and mean anomalies                                          
          anec=ATAN2(esine,ecose) 
          emme=anec-ecc*SIN(anec) 
          IF(emme.LT.0.d0) emme=emme+dpig 
          IF(emme.GT.dpig) emme=emme-dpig 
                                                                        
! Argument of pericenter                                                
          x1=COS(anec)-ecc 
          rad=1.d0-ecc**2 
          x2=SQRT(rad)*SIN(anec) 
          xm=SQRT(x1**2+x2**2) 
          x1=x1/xm 
          x2=x2/xm 
          sinper=x1*xorb(2)-x2*xorb(1) 
          cosper=x1*xorb(1)+x2*xorb(2) 
          argper=ATAN2(sinper,cosper) 
                                                                        
! CASE 2: RECA = 0 (parabolic orbit)                                    
      ELSEIF(reca.EQ.0.d0) THEN 
          type='COM' 
          ecc=1.d0 
                                                                        
! Semilatus rectum and pericenter distance                              
          p=el2/gm 
          q=p/2.d0 
                                                                        
! True anomaly                                                          
          cosf=p/rs-1.d0 
          sinf=rv*p/(rs*elmod) 
          effe=ATAN2(sinf,cosf) 
                                                                        
! Argument of pericenter                                                
          argper=ATAN2(xorb(2),xorb(1))-effe 
                                                                        
! CASE 3: RECA < 0 (hyperbolic orbit)                                   
      ELSE 
          type='COM' 
                                                                        
! Semilatus rectum, true anomaly and eccentricity                       
          p=el2/gm 
          ecosf=p/rs-1.d0 
          esinf=rv*p/(elmod*rs) 
          effe=ATAN2(esinf,ecosf) 
          ecc=SQRT(ecosf**2+esinf**2) 
          IF(ABS(ecc-1.d0).LT.eps) THEN 
              reca=0.d0 
              GOTO 2 
          END IF 
                                                                        
! Pericenter distance                                                   
          q=p/(1.d0+ecc) 
                                                                        
! Argument of pericenter                                                
          argper=ATAN2(xorb(2),xorb(1))-effe 
      END IF 
                                                                        
      IF(argper.LT.0.d0) argper=argper+dpig 
      IF(argper.GT.dpig) argper=argper-dpig 
      IF(type.EQ.'KEP') THEN 
          elem(1)=sma 
          elem(6)=emme 
      ELSEIF(type.EQ.'COM') THEN 
          elem(1)=q 
          elem(6)=effe 
      ELSE 
          STOP '**** ccek1: internal error (01) ****' 
      END IF 
      elem(2)=ecc 
      elem(3)=ainc 
      elem(5)=argper 
      elem(4)=anod 
                                                                        
      END                                           

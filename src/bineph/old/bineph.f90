! Copyright (C) 1997 by Mario Carpino (carpino@brera.mi.astro.it)       
! Version: December 17, 1997                                            
! --------------------------------------------------------------------- 
!                                                                       
! Generation of a binary (random access) ephemeris file                 
!                                                                       
      PROGRAM bineph 
      USE reference_systems 
      IMPLICIT NONE 
                                                                        
      INCLUDE 'dim.h90' 
                                                                        
      CHARACTER*80 run,optfil,auxfil,binfil 
      INTEGER nbout,lr,nout,ln,i,k,recl,nhl,ll,nv,iptf1,it,iptb1 
      INTEGER unit 
      DOUBLE PRECISION t0,x0(3,nbmx),v0(3,nbmx),teph1,teph2,dteph 
      DOUBLE PRECISION t,x(3,nbmx),v(3,nbmx),xvt(6) 
      DOUBLE PRECISION tmp,xl,enne,elem(6,nbmx),rot(3,3) 
      LOGICAL found 
                                                                        
      INCLUDE 'comast.h90' 
                                                                        
      INTEGER lench 
      EXTERNAL lench 
                                                                        
      WRITE(*,110) 
  110 FORMAT(' Run name =') 
      READ(*,100) run 
  100 FORMAT(A) 
      lr=lench(run) 
      WRITE(*,111) run(1:lr) 
  111 FORMAT(' Run name = ',A) 
                                                                        
! Input of options                                                      
      CALL namini 
      CALL libini 
                                                                        
      CALL filopf(unit,'bineph.def',found) 
      IF(found) THEN 
          CALL rdnam(unit) 
          CALL filclo(unit,' ') 
      ELSE 
          WRITE(*,210) 
          STOP '**** bineph: abnormal end ****' 
      END IF 
  210 FORMAT(' ERROR: no default option file (bineph.def) found') 
      optfil=run(1:lr)//'.bop' 
      CALL filopn(unit,optfil,'OLD') 
      CALL rdnam(unit) 
      CALL filclo(unit,' ') 
                                                                        
      CALL rdklst('bineph.key') 
      CALL chkkey 
                                                                        
      CALL rdoptb(nbout,t0,x0,v0,teph1,teph2,dteph) 
      IF(nbout.LE.0) STOP '**** bineph: object list is empty ****' 
                                                                        
! Initializations                                                       
      CALL rotpn(rot,'EQUM','J2000',0.d0,'ECLM','J2000',0.d0) 
      dteph=ABS(dteph) 
      IF(teph1.GT.teph2) THEN 
          tmp=teph1 
          teph1=teph2 
          teph2=tmp 
      END IF 
      nout=(teph2-teph1)/dteph+1 
      IF(teph2.GT.teph1+dteph*(nout-1)) nout=nout+1 
      teph2=teph1+dteph*(nout-1) 
                                                                        
! Deleting old copies of output files                                   
      auxfil=run(1:lr)//'.bai' 
      binfil=run(1:lr)//'.bep' 
      CALL dlifex(auxfil) 
      CALL dlifex(binfil) 
                                                                        
! Auxiliary file                                                        
      OPEN(1,FILE=auxfil,STATUS='UNKNOWN') 
      WRITE(1,200) nbout,teph1,teph2,dteph 
  200 FORMAT(I5/3F13.5) 
      DO 1 i=1,nbout 
      ln=lench(astnam(i)) 
      WRITE(1,201) astmas(i),astnam(i)(1:ln) 
    1 END DO 
  201 FORMAT(1P,E18.10,1X,A) 
      CLOSE(1) 
                                                                        
! Binary ephemeris file                                                 
! Assumes RECL given in single precision words                          
      recl=8*6*nbout 
      OPEN(1,FILE=binfil,ACCESS='DIRECT',RECL=recl,STATUS='NEW') 
! First record: version number, number of bodies, number of records     
      WRITE(1,REC=1) 102,nbout,nout 
! Second record: first and last epoch (MJD, TDT), stepsize (d)          
      WRITE(1,REC=2) teph1,teph2,dteph 
! Records from 3 to 2+nbout: masses                                     
      DO 8 i=1,nbout 
      WRITE(1,REC=2+i) astmas(i),gma(i),gma1(i) 
    8 END DO 
! Records from 3+nbout to 2+2*nbout: names                              
      DO 9 i=1,nbout 
      WRITE(1,REC=2+nbout+i) astnam(i) 
    9 END DO 
      nhl=2+2*nbout 
                                                                        
! Numerical integration                                                 
      ll=12 
      xl=1.d0 
      nv=3*nbt 
! First point to be computed in the forward integration                 
      tmp=(t0-teph1)/dteph+1 
      iptf1=tmp 
      IF(tmp-iptf1.GT.0.d0) iptf1=iptf1+1 
                                                                        
! Forward integration                                                   
      IF(iptf1.LE.nout) THEN 
          DO 2 k=1,nbt 
          DO 3 i=1,3 
          x(i,k)=x0(i,k) 
          v(i,k)=v0(i,k) 
    3     CONTINUE 
    2     CONTINUE 
          t=teph1+dteph*(iptf1-1) 
          CALL ra15b(x,v,t0,t-t0,xl,ll,nv,-2) 
          DO 4 i=1,nbout 
          CALL prodmv(xvt(1),rot,x(1,i)) 
          CALL prodmv(xvt(4),rot,v(1,i)) 
          CALL coocha(xvt,'CAR',gma1(i),elem(1,i),'KEP',enne) 
    4     CONTINUE 
! Records from 3+nbout to end: keplerian elements (a,e,i,N,w,M);        
! units: AU, radians; refsys: ECLM J2000                                
          WRITE(1,REC=nhl+iptf1) ((elem(i,k),i=1,6),k=1,nbout) 
          DO 5 it=iptf1+1,nout 
          CALL ra15b(x,v,t,dteph,xl,ll,nv,-2) 
          DO 6 i=1,nbout 
          CALL prodmv(xvt(1),rot,x(1,i)) 
          CALL prodmv(xvt(4),rot,v(1,i)) 
          CALL coocha(xvt,'CAR',gma1(i),elem(1,i),'KEP',enne) 
    6     CONTINUE 
          WRITE(1,REC=nhl+it) ((elem(i,k),i=1,6),k=1,nbout) 
          t=teph1+dteph*(it-1) 
    5     CONTINUE 
          iptb1=iptf1-1 
      ELSE 
          iptb1=nout 
      END IF 
                                                                        
! Backward integration                                                  
      IF(iptb1.GE.1) THEN 
          DO 12 k=1,nbt 
          DO 13 i=1,3 
          x(i,k)=x0(i,k) 
          v(i,k)=v0(i,k) 
   13     CONTINUE 
   12     CONTINUE 
          t=teph1+dteph*(iptb1-1) 
          CALL ra15b(x,v,t0,t-t0,xl,ll,nv,-2) 
          DO 14 i=1,nbout 
          CALL prodmv(xvt(1),rot,x(1,i)) 
          CALL prodmv(xvt(4),rot,v(1,i)) 
          CALL coocha(xvt,'CAR',gma1(i),elem(1,i),'KEP',enne) 
   14     CONTINUE 
          WRITE(1,REC=nhl+iptb1) ((elem(i,k),i=1,6),k=1,nbout) 
          DO 15 it=iptb1-1,1,-1 
          CALL ra15b(x,v,t,-dteph,xl,ll,nv,-2) 
          DO 16 i=1,nbout 
          CALL prodmv(xvt(1),rot,x(1,i)) 
          CALL prodmv(xvt(4),rot,v(1,i)) 
          CALL coocha(xvt,'CAR',gma1(i),elem(1,i),'KEP',enne) 
   16     CONTINUE 
          WRITE(1,REC=nhl+it) ((elem(i,k),i=1,6),k=1,nbout) 
          t=teph1+dteph*(it-1) 
   15     CONTINUE 
      END IF 
                                                                        
      CLOSE(1) 
                                                                        
      END PROGRAM bineph                                          

!*********************************************************************  
!  {\bf secth} PROTRO                                                   
!  Routine which provides the fundamental frequencies and the semimajor 
!  axis of Jupiter according to the LONGSTOP 1B synthetic theory.       
!  See Nobili et al.  A&A 210,313, 1989                                 
!  g(9) is used to accomodate the most important combination frequency  
!  UNITS: semimaj axis= AU; freq= arcsec/yr                             
SUBROUTINE secth(inflag,g,s,aj) 
  IMPLICIT NONE 
  INTEGER, INTENT(IN) :: inflag ! 1= 5My 2=50My
  DOUBLE PRECISION, INTENT(OUT) :: g(9),s(9),aj
  DOUBLE PRECISION gv(9,2),sv(9,2)
! frequencies for 5 Myr integrations                                    
  gv(5,1)=4.233030463d0                                               
  gv(6,1)=28.234843727d0                                              
  gv(7,1)=3.088221248d0                                               
  gv(8,1)=0.673627325d0                                               
  sv(6,1)=-26.333553362d0                                             
  sv(7,1)=-2.992590375d0                                              
  sv(8,1)=-0.689319899d0                                              
! frequencies for 50 Myr integrations                                   
  gv(5,2)=4.232814411d0 
  gv(6,2)=28.234770846d0 
  gv(7,2)=3.088071817d0 
  gv(8,2)=0.673266012d0 
  sv(6,2)=-26.333616452d0 
  sv(7,2)=-2.993420737d0 
  sv(8,2)=-0.690906565d0
! choice depending upon time span
  g(1:4)=0.d0
  g(5:8)=gv(5:8,inflag)
  s(1:5)=0.d0
  s(6:8)=sv(6:8,inflag) 
! combination frequency                                                 
  g(9)=2.d0*g(6)-g(5) 
  s(9)=0.d0 
! This has some problems...                                             
  s(5)=0.d0 
! Semimajor axes (au) derived from LONGSTOP 1B (average of filtered     
!  --periods  < 100,000 yr wiped out)                                   
  aj=5.2025696d0 
END SUBROUTINE secth
! ========================================================              
SUBROUTINE propa(sy,sc,sc2,ss,ss2,scs,scy,ssy,sy2,nt,             &
     &              ym,sig,d0,d1,d2,sp,n,i1,i2,ns,i)  
  USE synthtro                  
  IMPLICIT NONE 
  double precision sy(nsx,nastx,3,2),sc(nsx,nastx,3,2),ss(nsx,nastx,3,2)
  double precision ss2(nsx,nastx,3,2),scs(nsx,nastx,3,2) 
  double precision scy(nsx,nastx,3,2),sy2(nsx,nastx,3,2) 
  double precision ssy(nsx,nastx,3,2),sc2(nsx,nastx,3,2) 
  double precision ym(nsx,3,2),sig(nsx,3,2),sp(nsx,3,2) 
  double precision d0(nsx,3,2),d1(nsx,3,2),d2(nsx,3,2) 
  double precision swy,swc,swc2,sws,sws2,swcs,swsy,swy2,swcy 
  integer j,m,i1,i2,ns,i,n,ntt 
  integer nt(nsx) 
!c    dimension amp(nsx,3,2),ph(nsx,3,2),tm(nsx)                        
  DO 21 j=1,2 
     DO m=1,3 
        swy=SUM(sy(i1:i2,n,m,j)) 
        swc=SUM(sc(i1:i2,n,m,j)) 
        swc2=SUM(sc2(i1:i2,n,m,j)) 
        sws=SUM(ss(i1:i2,n,m,j)) 
        sws2=SUM(ss2(i1:i2,n,m,j)) 
        swcs=SUM(scs(i1:i2,n,m,j)) 
        swcy=SUM(scy(i1:i2,n,m,j)) 
        swsy=SUM(ssy(i1:i2,n,m,j)) 
        swy2=SUM(sy2(i1:i2,n,m,j))
        ntt=SUM(nt(i1:i2)) 
        CALL ferraz(ntt,swy,swc,swc2,sws,sws2,swcs,swcy,swsy,swy2,   &
     &     ym(i,m,j),sig(i,m,j),d0(i,m,j),d1(i,m,j),d2(i,m,j),sp(i,m,j))
     ENDDO
21 ENDDO
END SUBROUTINE propa
! *********************************************                         
SUBROUTINE ferraz(nt,swy,swc,swc2,sws,sws2,swcs,swcy,swsy,swy2,   &
     &  ym,sig,d0,d1,d2,sp)                                             
  IMPLICIT NONE 
  DOUBLE PRECISION sw,swy,swcy,swc,swcs,ym,swsy,sws,swy2,siig 
  DOUBLE PRECISION a0,a1,a2,c1,c2,sp,d0,d1,d2,swc2,sws2,sig,aa1 
  INTEGER nt 
  sw=dfloat(nt) 
  ym=swy/sw 
  swcy=swcy-ym*swc 
  swsy=swsy-ym*sws 
  swy2=swy2-(sw-2.d0)*ym**2 
  sig=sqrt(swy2/sw) 
!   calcolo densita' spettrale                                          
  IF(sw.gt.0.d0)THEN 
     a0=1/sqrt(sw) 
  ELSE 
     WRITE(*,*)' ferraz: sw<0 ', sw 
     a0=0.d0 
  ENDIF
  aa1=swc2-a0**2*swc**2 
  IF(aa1.gt.0.d0)THEN 
     a1=1/sqrt(aa1) 
  ELSE 
     WRITE(*,*)' ferraz: aa1<0 ', aa1 
     a1=0.d0 
  ENDIF
  a2=sws2-a0**2*sws**2-a1**2*swcs**2-a0**4*a1**2*(swc*sws)**2 
  a2=a2+2*a0**2*a1**2*swc*sws*swcs 
  IF(a2.gt.0.d0)THEN 
     a2=1./sqrt(a2) 
  ELSE 
     WRITE(*,*)' ferraz: a2<0 ', a2 
     a2=0.d0 
  ENDIF
  c1=a1*swcy 
  c2=a2*swsy-a1*a2*c1*(swcs-a0**2*swc*sws) 
  sp=c1**2+c2**2 
  sp=sp/swy2 
!    calcolo elemento periodico                                         
  d2=a2*c2 
  d1=a1*c1 +a2*a1**2*c2*(a0**2*swc*sws-swcs) 
  d0=-a0**2*(d1*swc+d2*sws)+ym 
END SUBROUTINE ferraz

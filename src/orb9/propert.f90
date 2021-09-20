! =========================================                             
!  {\bf propert} Proper elements for trojans                            
!   Vers. 1.1-1.2 Meudon,December 91/June 92                            
!    1.2 with adaptive handling of forced terms 
! f95 version November 2008                         
! =========================================                             
PROGRAM propert 
  USE fund_const
  USE synthtro
  IMPLICIT NONE 
! dimensions
  
  DOUBLE PRECISION, ALLOCATABLE :: sy(:,:,:,:),sc(:,:,:,:),ss(:,:,:,:)
  DOUBLE PRECISION, ALLOCATABLE :: ss2(:,:,:,:),scs(:,:,:,:),scy(:,:,:,:)
  DOUBLE PRECISION, ALLOCATABLE :: ssy(:,:,:,:),sc2(:,:,:,:),sy2(:,:,:,:)
  DOUBLE PRECISION ym(nsx,3,2),sig(nsx,3,2),sp(nsx,3,2) 
  DOUBLE PRECISION d0(nsx,3,2),d1(nsx,3,2),d2(nsx,3,2),clib(nastx) 
!    double precision amp(nsx,3,2),ph(nsx,3,2),tm(nsx)                 
  DOUBLE PRECISION ft1(nsx),ft2(nsx) 
  INTEGER nt(nsx),klisg(nforc),kliss(nforc) 
  DOUBLE precision gp(nforc),sn(nforc),tf(ntx)
!  dimensions: elements, sampled and filtered                           
  DOUBLE PRECISION, ALLOCATABLE :: elf(:,:,:) 
!  workspaces                                                           
  DOUBLE PRECISION x(ntx),y(ntx),th(ntx),dy(ntx),gam(ntx) 
  INTEGER inflag
  INTEGER ntt,nbb
  DOUBLE PRECISION ae(nbx),de(nbx),fe(nbx),dfe(nbx),pt(nbx) 
  DOUBLE PRECISION ai(nbx),di(nbx),fi(nbx),dfi(nbx),fref(nbx) 
!  file names                                                           
  CHARACTER*60 initia,filnam 
!  asteroid numbers                                                     
  CHARACTER*9 number(nastx) 
! loop indexes                                                          
  INTEGER jp 
! scalar integers                                                       
  INTEGER len                                          
  INTEGER ntf,naf,na,i,n,m,mm,nn,ns,j,i1,i2,iii,itot 
  INTEGER iwri,ij,i0,nb,l45,ny,nda,ndl,ndpnu,nde,ndg,ndsi,nds 
! various scalar real
  DOUBLE PRECISION aj,rate,rms,cost,perth,dperth,sd1
  DOUBLE PRECISION sigma ! function 
  DOUBLE PRECISION sd2,dd1,dd2,vd1,vd2,da,sda,dda 
  DOUBLE PRECISION cr,scr,dcs,sp100,ampe,dampe,freg,dfreg,amp,damp 
  DOUBLE PRECISION fres,dfres,dcr,sperth,sampe,samp 
  DOUBLE PRECISION sfreg,sfres,sfref,dfref 
  DOUBLE PRECISION carex, tlyap,pnu,g,s,princ,si,dpnu,dg,dsi,s1 
  DOUBLE PRECISION ds,plce,t0 
! added December 2003 to handle shift of libration center               
! corrected September 2005 to avoid shift of 360 deg                    
  DOUBLE PRECISION dcen,ympr 
  LOGICAL bad 
! ================================================
! allocate arrays
  ALLOCATE(sy(nsx,nastx,3,2),sc(nsx,nastx,3,2),ss(nsx,nastx,3,2))
  ALLOCATE(ss2(nsx,nastx,3,2),scs(nsx,nastx,3,2),scy(nsx,nastx,3,2))
  ALLOCATE(ssy(nsx,nastx,3,2),sc2(nsx,nastx,3,2),sy2(nsx,nastx,3,2))
  ALLOCATE(elf(8,nastx,ntx))
! ================================================                      
!  chose forced terms                                                   
  klisg=0 
  kliss=0
  klisg(5)=1 
  klisg(6)=1 
  klisg(7)=1 
  klisg(8)=0 
  klisg(9)=0 
  kliss(5)=0 
  kliss(6)=0 
  kliss(7)=0 
  kliss(8)=0 
  kliss(9)=0 
! =================================================                     
!   1: input data                                                       
! ==================================================                    
!   1b: filtered data                                                   
  write(*,*)' initial part of file names?' 
  read(*,*)initia 
!      open(22,file='propert.idx',status='old')                         
!      read(22,122)initia                                               
! 122  format(a60)                                                      
!  flag for the running box
  WRITE(*,*)' running box flag 1 = 5000 recs;&
& 2 = 50000 recs extended'
  READ(*,517)inflag
517 FORMAT(i1)
!  inflag=2
  CALL parprtt(inflag,nbb,ntt,ndap,nshi,isum,iste)
  call cmpstr(initia,len) 
  filnam=initia(1:len)//'ast.fil' 
  open(1,file=filnam,status='old') 
!   main output files                                                   
  filnam=initia(1:len)//'prt.out' 
  open(9,file=filnam,status='unknown') 
  filnam=initia(1:len)//'prt.pro' 
  open(10,file=filnam,status='unknown') 
  filnam=initia(1:len)//'prt.del' 
  open(11,file=filnam,status='unknown') 
  filnam=initia(1:len)//'prt.sig' 
  open(12,file=filnam,status='unknown') 
  filnam=initia(1:len)//'prt.libc' 
  open(62,file=filnam,status='unknown') 
  filnam=initia(1:len)//'prt.err' 
  open(63,file=filnam,status='unknown')
! data input 
  CALL  inptro(1,ntf,naf) 
  close(1) 
!     write(*,177)tf(1),tf(ntf),ntf                                     
!177  format(' data from t=',f14.4,' to t=',f14.4,' no datapoints=',i5/ 
!    +  ' how many datapoints together?')                               
  write(9,179)tf(1),tf(ntf),ntf 
179 format(' data from t=',f14.4,' to t=',f14.4,' no datapoints=',i5) 
!        read(*,*)ndap                                                  
!c      ndap=499                                                        
!      write(*,178)ndap                                                 
!178  format(' batches of ',i5,' datapoints, shift by how many?')       
!        read(*,*)nshi                                                  
!c      nshi=50                                                         
  write(9,*)' together ',ndap,'  step',nshi 
! ====================================================                  
!   1a: partial sums for libr. amplitude                                
! ====================================================                  
  if(len.eq.1)then 
     open(1,file='orb9t.sum',status='old') 
  else 
     filnam=initia(1:len)//'o9t.sum' 
     open(1,file=filnam,status='old') 
  endif
  read(1,101)na 
  if(na.ne.naf)then 
     write(*,*)' discrepant input,na=',na,'  naf=',naf 
     stop 
  endif
101 format(i3) 
  DO 76 i=1,nsx 
     read(1,119,end=77)nt(i),ft1(i),ft2(i) 
119  format(i4,f14.4,f14.4) 
!  read partial sums, for each asteroid, 3 harmonics                    
     DO 75 n=1,na 
        DO m=1,3 
           read(1,175,end=77)nn,mm,sy(i,n,m,1),sc(i,n,m,1),             &
     &     sc2(i,n,m,1),ss(i,n,m,1),ss2(i,n,m,1),                       &
     &     scs(i,n,m,1),scy(i,n,m,1),ssy(i,n,m,1),sy2(i,n,m,1),         &
     &     sy(i,n,m,2),sc(i,n,m,2),                                     &
     &     sc2(i,n,m,2),ss(i,n,m,2),ss2(i,n,m,2),                       &
     &     scs(i,n,m,2),scy(i,n,m,2),ssy(i,n,m,2),sy2(i,n,m,2)          
175        format(2i3,3d24.16,5(/3d24.16)) 
           if(nn.ne.n.or.mm.ne.m)then 
              write(*,*)nn,n,mm,m 
              stop 
           endif
        ENDDO
75   ENDDO
76 ENDDO
77 CONTINUE 
  ns=i-1 
  close(1) 
!  input complete                                                       
  write(9,192)ns,ft1(1),ft2(ns) 
192 format('part.sums, intervals=',i3,' from t=',f14.4,' to t=',f14.4)
!   chose intervals                                                     
!      write(*,*)' no sums=',ns                                         
!        read(*,*)isum                                                  
!c       isum=20                                                        
!       write(*,*)' how many in each step?'                             
!         read(*,*)iste                                                 
!c       iste=2                                                         
  write(9,*)' together ',isum,' step',iste 
! =====================================================                 
!   1c: constant data                                                   
  call secth(inflag,gp,sn,aj) 
! =====================================================                 
!  loop on the naf asteroids 
  DO 1 n=1,naf 
! =====================================================                 
!   2: determination of frequencies by fit to arguments                 
!   2a: libration argument theta                                        
! =====================================================                 
     DO j=1,ntf 
!   quick fix for error in input                                        
        if(elf(1,n,j).eq.0.) goto 1 
        th(j)=elf(7,n,j)
     ENDDO
     call linfi3(tf,th,rate,rms,cost,dy,ntf) 
     perth=dpig/rate 
     dperth=perth-dpig/(rate+rms/(tf(ntf)-tf(1))) 
     write(9,120)number(n),rate,perth,rms,cost*degrad 
120  format(' ============================================= '/         &
     &   ' Asteroid no. ',a9/                                           &
     &    ' Libration ph:',                                             &
     &  ' freq ',1p,d13.6,0p,' per ',f11.6,                             &
     &  ' rms ',1p,d12.4,0p,' ph ',f9.4)                                
!======================================================                 
!   3: find the proper amplitude and phase                              
!   3a: semimajor axis and critical argument                            
! =====================================================                 
!  for each asteroid, compute proper elements da,cr                     
!  compute on all                                                       
     i1=1 
     i2=ns 
     i=ns+1 
     CALL propa(sy,sc,sc2,ss,ss2,scs,scy,ssy,sy2,nt,                &
     &              ym,sig,d0,d1,d2,sp,n,i1,i2,ns,i)                    
     clib(n)=d0(ns+1,1,2) 
!  running box                                                          
     iii=0 
     DO 10 i=isum,ns,iste 
        iii=iii+1 
        i1=i+1-isum 
        i2=i 
        call propa(sy,sc,sc2,ss,ss2,scs,scy,ssy,sy2,nt,              &
     &              ym,sig,d0,d1,d2,sp,n,i1,i2,ns,iii)                  
10   ENDDO
     itot=ns+1 
     do 32 j=1,2 
        do 31 m=1,3 
           sd1=sigma(d1(1,m,j),iii) 
           sd2=sigma(d2(1,m,j),iii) 
           dd1=maxval(d1(1:iii,m,j))-minval(d1(1:iii,m,j)) 
           dd2=maxval(d2(1:iii,m,j))-minval(d2(1:iii,m,j)) 
           vd1=sum(d1(1:iii,m,j))/iii 
           vd2=sum(d2(1:iii,m,j))/iii 
           if(m.eq.1)then 
              if(j.eq.1)then 
                 write(9,103)ym(itot,m,j),sig(itot,m,j) 
103              format(' =========== delta a,AU ',                     &
     &                  ' mean=',1p,d12.4,' sig=',d12.4)                
                 da=d2(itot,m,j) 
                 sda=sd2 
                 dda=dd2 
              elseif(j.eq.2)then 
                 write(9,104)ym(itot,m,j),sig(itot,m,j) 
104              format('=========== libr ampl,rad ',                   &
     &                  ' mean=',1p,d12.4,' sig=',d12.4)                
                 cr=d1(itot,m,j)*degrad 
                 scr=sd1*degrad 
                 dcr=dd1*degrad 
                 pnu=360/perth 
! handle shift of center with respect to \pm 60 deg                     
                 ympr=princ(ym(itot,m,j)) 
                 if(ympr.gt.pig)ympr=ympr-dpig 
                 IF(ympr.gt.0.d0)THEN 
                    dcen=(ympr-dpig/6.d0)*degrad 
                 ELSE 
                    dcen=(ympr+dpig/6.d0)*degrad 
                 ENDIF
                 bad=(abs(dcen)+2.d0.gt.cr) 
                 IF(bad)THEN 
                    WRITE(62,299)number(n),dcen,cr,pnu 
299                 FORMAT(a9,1x,f10.6,1x,f10.6,1x,f7.4) 
                 ENDIF
              endif
           endif
           sp100=sp(itot,m,j)*100 
           write(9,102)m,sp100,d0(itot,m,j)                             &
     &            ,d1(itot,m,j),vd1,sd1,dd1,                            &
     &            d2(itot,m,j),vd2,sd2,dd2                              
102        format(' harm=',i2,' S% ',f6.2,' cost=',1p,d12.4/            &
     &       ' cos ',d13.5,d13.5,d12.4,d12.4/                           &
     &       ' sin ',d13.5,d13.5,d12.4,d12.4)
31      ENDDO
32   ENDDO
! =====================================================                 
     write(9,121) 
121  format(' =========== eccentricity ') 
!   compute proper eccentricity and inclination from all data           
     DO j=1,ntf 
        x(j)=elf(3,n,j) 
        y(j)=elf(2,n,j) 
     ENDDO
     iwri=1 
     call doinc(x,y,tf,ntf,gp,klisg,ampe,dampe,freg,dfreg,iwri) 
! =====================================================                 
     write(9,126) 
126  format(' =========== inclination ') 
     DO j=1,ntf 
        x(j)=elf(5,n,j) 
        y(j)=elf(4,n,j)
     ENDDO
     iwri=1 
     call doinc(x,y,tf,ntf,sn,kliss,amp,damp,fres,dfres,iwri) 
     amp=2.d0*atan(amp) 
! =====================================================                 
!   running box test for accuracy                                       
     iii=0 
     iwri=0 
     DO 50 ij=ndap,ntf,nshi 
        i0=ij-ndap 
        i1=ij-ndap+1 
        i2=ij 
        iii=iii+1 
!   libration peri                                                      
        call linfi3(tf(i1),th(i1),fref(iii),rms,cost,dy,ndap) 
        pt(iii)=dpig/fref(iii) 
!   compute proper eccentricity and inclination from box data           
        DO j=1,ndap 
           x(j)=elf(3,n,j+i0) 
           y(j)=elf(2,n,j+i0)
        ENDDO
        CALL doinc(x,y,tf,ndap,gp,klisg,ae(iii),de(iii),fe(iii),dfe(iii),iwri) 
        DO j=1,ndap 
           x(j)=elf(5,n,j+i0) 
           y(j)=elf(4,n,j+i0)
        ENDDO
        CALL doinc(x,y,tf,ndap,sn,kliss,ai(iii),di(iii),fi(iii),dfi(iii),iwri)
50   ENDDO
     nb=iii 
!   compute max deviation, standard deviation                           
     sperth=sigma(pt,nb) 
     sfref=sigma(fref,nb)*degrad 
     dfref=(maxval(fref(1:nb))-minval(fref(1:nb)))*degrad 
     sampe=sigma(ae,nb) 
     samp=sigma(ai,nb) 
     sfreg=sigma(fe,nb) 
     sfres=sigma(fi,nb) 
     dperth=maxval(pt(1:nb))-minval(pt(1:nb)) 
     dampe=maxval(ae,nb)-minval(ae(1:nb)) 
     damp=maxval(ai,nb)-minval(ai(1:nb)) 
     dfreg=maxval(fe(1:nb))-minval(fe(1:nb)) 
     dfres=maxval(fi(1:nb))-minval(fi(1:nb)) 
     write(9,250)sperth,dperth,sampe,dampe,sfreg,dfreg,                &
     &              samp,damp,sfres,dfres                               
250  format(' ======== Libr. per : sig=',1p,d12.4,' del=',d12.4/       &
     &       ' Proper e: sig',d12.4,' del',d12.4,                       &
     &       ' Fr.peri: sig',d12.4,' del',d12.4/                        &
     &       ' Proper I: sig',d12.4,' del',d12.4,                       &
     &       ' Fr.node: sig',d12.4,' del',d12.4)                        
! ==========================================================            
!   2d: Lyapounov characteristic exponents                              
! ==========================================================            
     DO  j=1,ntf 
        gam(j)=elf(8,n,j)
     ENDDO
     call linfi3(tf,gam,carex,rms,cost,dy,ntf) 
     tlyap=1.d0/carex 
     write(9,150)carex,rms,cost,tlyap 
150  format(' ======== max. LCE',                                      &
     &  1p,d12.4,' rms.res',d12.4,' cost',d12.4,' tlyap',d12.4)         
!      write(12,151)number(n),carex,tlyap,clib(n)                       
! 151  format(1x,a9,3x,1p,d12.4,3x,d12.4,1x,0p,f10.6)                   
! =====================================================                 
!  new output format                                                    
     g=freg 
     s=fres 
     clib(n)=princ(clib(n)) 
     if(clib(n).lt.dpig/2)then 
        l45=4 
     else 
        l45=5 
     endif
     ny=tf(ntf)/1.d6+1 
     si=sin(amp) 
! file with proper elements                                             
     write(10,220)number(n),da,cr,pnu,ampe,g,si,s,l45,ny 
220  format(1x,a9,1x,f5.4,1x,f6.2,1x,f6.3,1x,                          &
     &       f5.4,1x,f7.2,1x,f5.4,1x,f7.2,1x,I2,1x,I2)                  
! lyapounov exponent for 10^5 years                                     
     plce=carex*1.d5 
! file with delta                                                       
     write(11,221)number(n),dda,dcr,dfref,dampe,dfreg,damp,dfres,plce 
221  format(1x,a9,1x,f6.5,1x,f7.3,1x,f7.4,1x,                          &
     &       f6.5,1x,f8.3,1x,f6.5,1x,f8.3,1x,f5.2)                      
! file with RMS                                                         
     write(12,221)number(n),sda,scr,sfref,sampe,sfreg,samp,sfres,plce 
! ndpnu,nde,ndg,ndsi,                                                   
!    +       nds,plce                                                   
! ==================================================                    
!  end loop on asteroids                                                
! ==================================================                    
1 ENDDO

CONTAINS
! =====================================================                 
!  {\bf inptro} PROTRO data input                                       
 SUBROUTINE inptro(iun,nt,na)
! SUBROUTINE inptro(iun,t,el,nt,na,number)
!  USE fund_const
!  USE synthtro, ONLY: nastx, ntx 
!  IMPLICIT NONE 
  INTEGER, INTENT(IN) ::  iun ! input unit
  INTEGER, INTENT(OUT) :: nt,na ! no times, no asteroids
!  DOUBLE PRECISION, INTENT(OUT) :: el(8,nastx0,ntx0),t(ntx) 
!  CHARACTER*9, INTENT(OUT) :: number(nastx)
! END INTERFACE   
  double precision x(8),t0 
  integer ng(2),ilce,i,j,n,nerr,jj,le 
  character coo*3,sys*3,ref*6,unit*3,comme*60,colhea*100 
! ===========================================================           
!  read header of file on unit iun (Trojans)                            
  call reahea(iun,number,comme,colhea,coo,sys,ref,unit,na,ilce,t0) 
  write(*,100)(number(i),i=1,na) 
  write(63,100)(number(i),i=1,na) 
100 format(10(2x,a9)) 
  DO jj=1,na
    CALL rmsp(number(jj),le)
  ENDDO
! ===========================================================           
  nerr=0 
!  read loop                                                            
  DO 1 j=1,ntx
     READ(iun,*,end=3,err=4)tf(j) 
     GOTO 6
!  error cases                                                          
4    nerr=nerr+1 
     WRITE(63,*)' error in time input, record ',j,' ast ',n 
     IF(nerr.lt.100)THEN 
        write(*,*)' error in time input, record ',j,' ast ',n 
     ENDIF
     DO i=1,7 
        x(i)=0.d0 
     ENDDO
     ng(1)=0 
     CYCLE ! next asteroid
6    CONTINUE
     DO 10 n=1,na 
        read(iun,*,end=5,err=5)(x(i),i=1,6),ng(1),x(7),ng(2),x(8)
        GOTO 7
5       WRITE(63,*)' error in elements input, record ',j,' ast ',n 
        IF(nerr.lt.100)THEN 
           write(*,*)' error in elements input, record ',j,' ast ',n 
        ENDIF
        DO i=1,7 
           x(i)=0.d0 
        ENDDO
        ng(1)=0 
        CYCLE ! next asteroid
7       CONTINUE 
        elf(1:5,n,j)=x(1:5) 
        elf(6,n,j)=x(6)+dpig*ng(1) 
        elf(7,n,j)=x(7)+dpig*ng(2) 
        elf(8,n,j)=x(8) 
10   ENDDO
1 ENDDO
  WRITE(*,*)' too many records, max was',ntx 
3 nt=j-1 
END SUBROUTINE inptro

END PROGRAM propert

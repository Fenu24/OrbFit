! **********************************************************            
!                                                                       
!     ORBIT9T  Pisa, May 2002                                           
!                                                                       
! ***********************************************************           
PROGRAM orbit9t
  USE massmod 
  USE fund_const
  USE fir_filters
  USE force9d
  USE propag_state, ONLY: propin, hms
  IMPLICIT NONE 
!  parameters                                                           
  INTEGER, PARAMETER ::  nto=5000
  INCLUDE 'comnbo.h90' 
  INCLUDE 'libcen.h90' 
! arrays for initial conditions and elements, output                    
  DOUBLE PRECISION yp(6,nbox),ennep(nbox),elp(6,nbox),plv(nbox) 
  DOUBLE PRECISION astll(nastx),allv(nastx) 
  INTEGER ngir(nastx) 
!  arrays for proper elements; cyclic storage                           
  DOUBLE PRECISION thef(nastx,nto),da(nastx,nto),cr(nastx,nto),     &
     &          ft(nto),ft2(nto)                                        
  DOUBLE PRECISION sy(nastx,3,2),sc(nastx,3,2),ss(nastx,3,2),       &
     &          sc2(nastx,3,2)                                          
  DOUBLE PRECISION ss2(nastx,3,2),scs(nastx,3,2),scy(nastx,3,2) 
  DOUBLE PRECISION ssy(nastx,3,2),sy2(nastx,3,2) 
! center of mass of inner solar system                                  
  DOUBLE PRECISION bar(6) 
  DOUBLE PRECISION ya(6,nastx),ennea(nastx),ela(6,nastx),alv(nastx) 
  INTEGER ngip(nbox),ngia(nastx),ngpp(nbox),ngpa(nastx) 
!  cartesian position and velocity vectors                              
  DOUBLE PRECISION y1(nvarx),y2(nvarx) 
! filter parameters 
  INTEGER nskk,nskk2,nskip,nskipa,nskip2,nskia2
  DOUBLE PRECISION hlen, hlen2
!  array for filter input/output                                        
  DOUBLE PRECISION eq(nvarx),eqf(nvarx),angl(norbx),angf(norbx) 
  DOUBLE PRECISION eqf2(nvarx),angf2(norbx) 
  INTEGER ng(norbx),ngf(norbx),ngf2(norbx) 
!  for variational equations                                            
  DOUBLE PRECISION gamma(nvzx),vslog(nvzx),vnor(nvzx),yv(6,nastx) 
!  identifiers of coordinate systems for output                         
  CHARACTER*3 sysz 
  CHARACTER*6 refy,refz 
  LOGICAL outf,outfa,reno,outf2,outfa2 
! dimensions                                                            
  INTEGER ndim,ndim2,nvar2,nvar,npla,norb,ncof 
! controls                                                              
  INTEGER nout,idump,iprqua,nin,nfl,iout0,ioucou 
  INTEGER njump,njump2 
  DOUBLE PRECISION v0,v1,vpos,vvel 
! time                                                                  
  DOUBLE PRECISION dt,t,t1,t2,tf,dinf,dinf2,tf2, tf2a 
! loop indexes                                                          
  INTEGER j,ja,jq,jj,iout,i 
! scalar temporaries                                                    
  DOUBLE PRECISION sv,vv,renorm 
! *********************************                                     
! only to test derivatives                                              
  DOUBLE PRECISION rinc,step,f0(nvar2x),f1(nvar2x) 
  INTEGER ii 
! close approach data                                                   
  INTEGER idc 
  DOUBLE PRECISION xxpla(6) 
! functions                                                             
  DOUBLE PRECISION princ 
  INTEGER iround 
! renormalization of libration ellipse, libration angle 
  DOUBLE PRECISION ren,the 
! other stuff                                                           
  DOUBLE PRECISION xx,res,yy 
  INTEGER nangf,nf,nf2,iis,n,nf2c,iouf2,lfp,nd,jshift 
  INTEGER inf2,nt,inf2c,m,iic,nf1,iouf 
! *********************************************************             
  CALL libini 
!   input                                                               
  CALL inpo9t(dt,nout,njump,njump2,idump,                    &
     &             sysz,refz,refy,                           &
     &             iprqua,                                   &
     &             v0,v1,vpos,vvel,vslog,                    &
     &             t1,yp,ya,ngip,ngia,yv)
!  renormalisation for amplitude of libration of critical arg.          
  ren=1/0.2783 
!  number of components                                                 
  ndim=3 
  ndim2=ndim*2 
  nvar2=ndim*(nbod-1+na+nvz) 
  nvar=2*nvar2 
  npla=nbod-1 
  norb=nbod-1+na 
  nangf=npla+2*na 
  lfp=0 
!  number of real components to be filtered                             
  ncof=norb*(2*ndim-1)+nvz 
  GOTO 99 
! ==================================================                    
! test derivatives; d(i)/d(ii)                                          
  step=1.d-8 
  idc=0 
  DO 90 j=1,na 
     write(*,*)(ya(i,j),i=1,3) 
     DO 91 ii=1,3 
! -----------------------------------                                   
! set variation vector to  E_ii                                         
        DO i=1,3 
           yv(i,j)=0.d0 
           yv(i+3,j)=0.d0 
        ENDDO
        yv(ii,j)=1.d0 
! reorder                                                               
        call stack(y1,nvar,yp,npla,ya,na,yv,nvz,ndim,ndim2) 
! compute derivatives                                                   
        CALL force9(y1,y1(nvar2+1),t1,f0,nvar2,idc,xxpla,0,1) 
! increment variable ii by step                                         
        ya(ii,j)=ya(ii,j)+step 
! reorder                                                               
        call stack(y1,nvar,yp,npla,ya,na,yv,nvz,ndim,ndim2) 
! compute increment                                                     
        CALL force9(y1,y1(nvar2+1),t1,f1,nvar2,idc,xxpla,0,1) 
        DO 94 i=1,3 
! compute incremental ratio                                             
           rinc=(f1(3*npla+3*j-3+i)-f0(3*npla+3*j-3+i))/step 
           write(*,*)i,ii,f0(3*norb+3*j-3+i),rinc,                     &
     &            (f0(3*norb+3*j-3+i)-rinc)/rinc                        
94      ENDDO
91   ENDDO
! ----------------------------------                                    
90 ENDDO
  STOP 
99 CONTINUE 
! *********************************************************             
!  initial conditions: all arrays copied in the state                   
!  vector                                                               
  CALL stack(y1,nvar,yp,npla,ya,na,yv,nvz,ndim,ndim2) 
! ****************************************************                  
!  output at initial time:                                              
!                                                                       
!  planets                                                              
  CALL coord(yp,'BAR',nbod,'CAR',refy,elp,sysz,'EQU',refz,ennep,bar) 
!  revolution counters:                                                 
  DO j=1,nbod-1 
     plv(j)=elp(6,j) 
  ENDDO
! arrays for filtering and output                                       
  jq=0 
  DO 58 j=1,nbod-1 
     angl(j)=plv(j) 
     ng(j)=ngip(j) 
     DO i=1,5 
        eq(i+jq)=elp(i,j) 
     ENDDO
     jq=jq+5 
58 ENDDO
!  asteroids                                                            
  CALL cooast(ya,'BAR',na,'CAR',refy,bar,nbod,ela,sysz,'EQU',refz,ennea)
  DO 6 ja=1,na 
!  resonance variables: since the revolution counters are not yet       
!  effective, some trick must be played to reconvert to the             
!  right interval                                                       
     xx=princ(ela(6,ja)-elp(6,1)) 
     if(xx.gt.dpig/2)xx=xx-dpig 
     if(xx.gt.0)then 
!          xx=(xx-dpig/6.d0)/ren                                        
        xx=(xx-dpig/6.d0-dcen(ja))/ren 
     else 
!          xx=(xx+dpig/6.d0)/ren                                        
        xx=(xx+dpig/6.d0-dcen(ja))/ren 
     endif
     yy=ela(1,ja)-elp(1,1) 
     res=atan2(yy,xx) 
     astll(ja)=princ(res) 
     angl(2*ja+nbod-1)=astll(ja)
     ng(2*ja+nbod-1)=0
     ngir(ja)=0 
     allv(ja)=astll(ja) 
!  revolution counters:
     alv(ja)=ela(6,ja) 
! arrays for filtering and output
     angl(2*ja+nbod-2)=alv(ja)
     ng(2*ja+nbod-2)=ngia(ja) 
     DO i=1,5 
        eq(i+jq)=ela(i,ja)
     ENDDO
     jq=jq+5 
6 ENDDO
  DO j=1,nvz 
     eq(jq+j)=vslog(j) 
  ENDDO
  CALL outsubt(13,14,t1,eq,angl,ng,ncof,nangf,na,nvz,.false.) 
  IF(nvz.ne.0)THEN 
     DO j=1,nvz 
        write(14,*)(yv(i,j),i=1,ndim) 
        write(14,*)(yv(i+ndim,j),i=1,ndim) 
     ENDDO
  ENDIF
! ****************************************************                  
!  initialisations                                                      
!                                                                       
!  interval to filter input                                             
  dinf=dt/nsamp 
  nin=nout*nsamp 
!  initial time is t1; how many filter inputs are behind?               
  iout0=iround(t1/dinf) 
!  starter is needed to start                                           
  nfl=0 
!  output counters                                                      
  ioucou=0 
!  pointers in the proper elements stack: current position              
  nf=0 
  nf2=0 
!  output file for partial sums                                         
  OPEN(19,file='orb9t.sum',status='unknown') 
  WRITE(19,219)na 
219 FORMAT(i3,'  number of Trojans in this file') 
!  control on cyclic storage                                            
  if(nto.lt.idump)then 
     write(9,*)' cyclic array insufficient, nto=',nto,' idump=',idump
     stop 
  endif
!  filter skip (.fil remains in phase)                                  
  nskk=mod((nfil+1)/2,nsamp) 
  if(nskk.eq.0)then 
     nskip=0 
  else 
     nskip=nsamp-nskk 
  endif
  nskipa=nskip 
  hlen=dinf*(nfil-1)/2 
  nskk2=mod((nfil2+1)/2,nsamp2) 
  if(nskk2.eq.0)then 
     nskip2=0 
  else 
     nskip2=nsamp2-nskk2 
  endif
  nskia2=nskip2 
  dinf2=dt/nsamp
  hlen2=dinf2*(nfil2-1)/2 
  write(9,*)' filters parameters'
  write(9,*)nskip,hlen,dinf,nskip2,hlen2,dinf2 
!                                                                       
! *********************************************************             
!   main loop begins here                                               
!                                                                       
!  propagation to next filter input time                                
  DO 10 iout=iout0+1,nin 
     t2=dinf*iout 
     CALL propin(nfl,y1,t1,t2,y2,hms,nvar,6) 
     nfl=1 
     t=t2 
! *********************************************************             
!  output: reordering of the arrays                                     
     CALL reord(y2,nvar,yp,npla,ya,na,yv,nvz,ndim,ndim2) 
!  planets: coordinate change                                           
     CALL coord(yp,'BAR',nbod,'CAR',refy,elp,sysz,'EQU',refz,ennep,bar) 
!  asteroids: coordinate change                                         
     CALL cooast(ya,'BAR',na,'CAR',refy,bar,nbod,ela,sysz,'EQU',refz,ennea)
!  planets                                                              
     jq=0 
     DO 17 j=1,nbod-1 
!  planets: revolution counters update                                  
        CALL predicl(plv(j),ngip(j),dinf,ennep(j),ngpp(j)) 
        ngip(j)=ngpp(j) 
        CALL prngup(elp(6,j),plv(j),ngip(j)) 
!  planets: arrays for filtering and output                             
        angl(j)=plv(j) 
        ng(j)=ngip(j) 
        DO i=1,5 
           eq(i+jq)=elp(i,j) 
        ENDDO
        jq=jq+5 
17   ENDDO
!  asteroids                                                            
     DO 16 ja=1,na 
!  asteroids: revolution counters update                                
        CALL predicl(alv(ja),ngia(ja),dinf,ennea(ja),ngpa(ja)) 
        ngia(ja)=ngpa(ja) 
        CALL prngup(ela(6,ja),alv(ja),ngia(ja)) 
!  resonance variables: force                                           
!  the difference in the right interval.                                
        xx=princ(ela(6,ja)-elp(6,1)) 
        if(xx.gt.dpig/2)xx=xx-dpig 
        if(xx.gt.0)then 
!          xx=(xx-dpig/6.d0)/ren                                        
           xx=(xx-dpig/6.d0-dcen(ja))/ren 
        else 
!          xx=(xx+dpig/6.d0)/ren                                        
           xx=(xx+dpig/6.d0-dcen(ja))/ren 
        endif 
        yy=ela(1,ja)-elp(1,1) 
        res=atan2(yy,xx) 
        astll(ja)=princ(res) 
        angl(2*ja+nbod-1)=astll(ja) 
        call prngup(astll(ja),allv(ja),ngir(ja)) 
        ng(2*ja+nbod-1)=ngir(ja) 
!  asteroids: arrays for filtering and output                           
        angl(2*ja+nbod-2)=alv(ja) 
        ng(2*ja+nbod-2)=ngia(ja) 
        DO i=1,5 
           eq(i+jq)=ela(i,ja)
        ENDDO 
        jq=jq+5 
16   ENDDO
!  variational equations                                                
     DO 20 j=1,nvz 
!  norm of the variation vector; here the norm**2 is defined by         
!  weight vpos for positions and vvel for velocities                    
        sv=0.d0 
        vv=0.d0 
        DO i=1,ndim 
           sv=sv+yv(i,j)**2 
           vv=vv+yv(i+ndim,j)**2 
        ENDDO
        vnor(j)=dsqrt(sv*vpos+vv*vvel) 
!  log of divergence ratio; no /(t-t0) used here                        
        gamma(j)=(dlog(vnor(j))-dlog(v0)+vslog(j)) 
        jq=jq+1 
        eq(jq)=gamma(j) 
20   ENDDO
! ***************************************************************       
!  input filters                                                        
     CALL filter(t,eq,ncof,nskip,dinf,tf,eqf,outf) 
     CALL filtan(t,angl,ng,nangf,nskipa,dinf,tf,angf,ngf,outfa)
     CALL filter2(t,eq,ncof,nskip2,dinf2,tf2,eqf2,outf2)
     CALL filtan2(t,angl,ng,nangf,nskia2,dinf2,tf2a,angf2,ngf2,outfa2)
! ****************************************************************      
! output control  
!                                                                       
! filtered output:                                                      
     if(outf2.neqv.outfa2)then 
        write(9,*)' wrong 2nd filter control',tf2,outf2,outfa2 
     endif
     if(outf.neqv.outfa)then 
        write(9,*)' wrong filter control',tf,outf,outfa 
     endif
     if(outf)then 
        nf=nf+1 
        nf1=mod(nf,nto) 
        if(nf1.eq.0)nf1=nto 
        ft(nf1)=tf 
        jq=npla*5 
        ja=npla 
        DO 70 n=1,na 
           da(n,nf1)=eqf(jq+1)-eqf(1) 
           jq=jq+5 
           cr(n,nf1)=angf(ja+1)-angf(1)+(ngf(ja+1)-ngf(1))*dpig 
           ja=ja+2 
70      ENDDO
        iouf=iround(tf/dt) 
        if(mod(iouf,njump).eq.0.and.iprqua.gt.1)then 
           CALL outsubt(11,21,tf,eqf,angf,ngf,ncof,nangf,na,nvz,.true.) 
        endif
     endif
     if(outf2)then 
!     WRITE(*,*) 'tf2, tf2a ',tf2,tf2a 
        nf2=nf2+1 
        nf2c=mod(nf2,nto) 
        if(nf2c.eq.0)nf2c=nto 
        ft2(nf2c)=tf2 
        ja=npla 
        DO n=1,na 
           thef(n,nf2c)=angf2(ja+2)+ngf2(ja+2)*dpig 
           ja=ja+2 
        ENDDO
        iouf2=iround(tf2/dt) 
        if(mod(iouf2,njump2).eq.0)then 
           CALL outsubt(12,22,tf2,eqf2,angf2,ngf2,ncof,nangf,na,nvz,.true.)
        endif
     endif
! unfiltered output and periodic operations                             
     if(mod(iout,nsamp).eq.0)then 
        ioucou=ioucou+1 
!  periodic operations: dump and renormalisation                        
!                                                                       
!  dump                                                                 
        if(mod(ioucou,idump).eq.0)then 
           call outsubt(13,14,t,eq,angl,ng,ncof,nangf,na,nvz,.false.) 
           if(nvz.ne.0)then 
              DO j=1,nvz 
                 write(14,*)(yv(i,j),i=1,ndim) 
                 write(14,*)(yv(i+ndim,j),i=1,ndim) 
              ENDDO
           endif
!  write partial sums for proper elements                               
           if(lfp.eq.0)then 
              lfp=1 
              nd=nsamp2/nsamp 
!  check alignement of ft,ft2                                           
              DO 72 jj=1,nf 
                 if(abs(ft2(1)-ft(jj)).lt.1.d-7)then 
                    jshift=jj-1 
                    goto 73 
                 endif
72            ENDDO
!  this should not happen!                                              
              jshift=999
              write(9,*)' unmatched filter times, ft(1)=',ft(1),       &
     &                     ' ft2(1)=',ft2(1),nf,jshift                  
              stop 
73            continue 
              inf2=1 
           endif
           nt=nf2-inf2+1 
           nf2c=mod(nf2,nto) 
           if(nf2c.eq.0)nf2c=nto 
           inf2c=mod(inf2,nto) 
           if(inf2c.eq.0)inf2c=nto 
           write(19,119)nt,ft2(inf2c),ft2(nf2c) 
119        format(i4,f14.4,f14.4) 
!  compute partial sums, for each asteroid, 3 harmonics                 
           DO 750 n=1,na 
              DO 75 m=1,3 
                 call sumzer(sy(n,m,1),sc(n,m,1),sc2(n,m,1),ss(n,m,1),     &
     &        ss2(n,m,1),scs(n,m,1),scy(n,m,1),ssy(n,m,1),sy2(n,m,1))   
                 call sumzer(sy(n,m,2),sc(n,m,2),sc2(n,m,2),ss(n,m,2),     &
     &        ss2(n,m,2),scs(n,m,2),scy(n,m,2),ssy(n,m,2),sy2(n,m,2))   
                 DO 74 ii=inf2,nf2 
                    iic=mod(ii,nto) 
                    if(iic.eq.0)iic=nto 
                    the=thef(n,iic)*m 
                    iis=jshift+1+(ii-1)*nd 
                    iis=mod(iis,nto) 
                    if(iis.eq.0)iis=nto 
                    if(abs(ft(iis)-ft2(iic)).gt.1.d-7)then 
                       write(9,*)' unmatched times, f1=',iis,ft(iis),       &
     &                       ' filter2=',iic,ft2(iic)                   
                       write(9,*)'inf2=',inf2,' nf2=',nf2,' nf=',nf,        &
     &                       'nf1=',nf1,' t=',t                         
                       stop 
                    endif
                    call sumper(the,da(n,iis),sy(n,m,1),sc(n,m,1),          &
     &          sc2(n,m,1),ss(n,m,1),ss2(n,m,1),                        &
     &          scs(n,m,1),scy(n,m,1),ssy(n,m,1),sy2(n,m,1))            
                    call sumper(the,cr(n,iis),sy(n,m,2),sc(n,m,2),          &
     &          sc2(n,m,2),ss(n,m,2),ss2(n,m,2),                        &
     &          scs(n,m,2),scy(n,m,2),ssy(n,m,2),sy2(n,m,2))            
74               ENDDO
                 write(19,120)n,m,sy(n,m,1),sc(n,m,1),                     &
     &          sc2(n,m,1),ss(n,m,1),ss2(n,m,1),                        &
     &          scs(n,m,1),scy(n,m,1),ssy(n,m,1),sy2(n,m,1),            &
     &          sy(n,m,2),sc(n,m,2),                                    &
     &          sc2(n,m,2),ss(n,m,2),ss2(n,m,2),                        &
     &          scs(n,m,2),scy(n,m,2),ssy(n,m,2),sy2(n,m,2)             
120              format(2i3,3d24.16,5(/3d24.16)) 
75            ENDDO
750        ENDDO
           inf2=nf2+1 
! *********************************************************             
!  renormalisation                                                      
           reno=.false. 
           DO 50 j=1,nvz 
              if(vnor(j).gt.v1)then 
                 reno=.true. 
!  renormalisation of  variation vector; orbit is unchanged             
                 renorm=v0/vnor(j) 
                 vslog(j)=vslog(j)-dlog(renorm) 
                 jj=3*(nbod-1+na+j)-3 
                 DO i=1,ndim 
                   y2(nvar2+jj+i)=y2(nvar2+jj+i)*renorm 
                   y2(i+jj)=y2(i+jj)*renorm 
                 ENDDO
!                write(*,301)j,t,vnor(j),gamma(j)                       
                 write(9,301)j,t,vnor(j),gamma(j) 
301              format(' renormalisation var.eq.no ',i4/               &
     &             ' time=',f18.6,' norm=',1p,d15.5,                    &
     &             ' gamma=',d15.5)                                     
              endif
50         ENDDO
           if(reno)then 
!  the starter must be used again                                       
              t1=t 
              DO j=1,nvar 
                 y1(j)=y2(j)
              ENDDO
              nfl=0 
           endif
        endif
     endif
10 ENDDO
!                                                                       
!   end main loop                                                       
!                                                                       
! *********************************************************             
!   end                                                                 
  write(9,199)t 
199 format(' end orbit9t, time = ',f18.4) 
END PROGRAM orbit9t
!                                                                       
!  *******************************************                          
SUBROUTINE sumper(the,y,swy,swc,swc2,sws,sws2,swcs,swcy,swsy,swy2) 
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT(IN) :: the, y
  DOUBLE PRECISION, INTENT(OUT) :: swy, swc,swc2,sws,sws2,swcs,swcy,swsy,swy2
  DOUBLE PRECISION c,s 
  swy=swy+y 
  c=cos(the) 
  s=sin(the) 
  swc=swc+c 
  swc2=swc2+c**2 
  sws=sws+s 
  sws2=sws2+s**2 
  swcs=swcs+c*s 
  swcy=swcy+c*y 
  swsy=swsy+s*y 
  swy2=swy2+y**2 
END SUBROUTINE sumper
!  *******************************************                          
SUBROUTINE sumzer(swy,swc,swc2,sws,sws2,swcs,swcy,swsy,swy2) 
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT(OUT) :: swy, swc,swc2,sws,sws2,swcs,swcy,swsy,swy2
  swy=0.d0 
  swc=0.d0 
  swc2=0.d0 
  sws=0.d0 
  sws2=0.d0 
  swcs=0.d0 
  swcy=0.d0 
  swsy=0.d0 
  swy2=0.d0 
END SUBROUTINE sumzer

! **********************************************************            
!  iosho9 vers. with input OEF single line; iqco is assumed 0           
! this version can be used only with numbered and qual1 orbits          
!  elements are read and  converted to coordinates specified by         
!  cooy and  sysy                                                       
! **********************************************************            
SUBROUTINE iosho(kr,inp,iun8,iun7,iun17,nam0,t0,tt,el,pel,        &
     &     eof,iqco,hm,iun19)                                           
  USE fund_const    
  USE massmod
  IMPLICIT NONE 
! unit numbers                                                          
  INTEGER iun19,iun7,iun17,iun8 
! only remaining include in prop9.x
  INCLUDE 'comnbo.h90' 
! arrays for orbital elements                                           
  DOUBLE PRECISION xp(6,nbox),yp(6,nbox),ennep(nbox),bar(6) 
  DOUBLE PRECISION zzp(6,nbox),zbp(6,nbox) 
  INTEGER nga(6,nastx),ngp(6,nbox) 
  DOUBLE PRECISION xa(6,nastx),ya(6,nastx),ennea(nastx) 
  DOUBLE PRECISION zza(6,nastx),zba(6,nastx) 
  DOUBLE PRECISION barin(6),pm(10),el(6),pel(6,2) 
!  identifiers of coordinate systems                                    
  character*3 coox,cooy,sysx,sysy,unitx,cooxp,sysxp,unitxp 
  character*6 refx,refy,refb,refxp 
!  names and numbers of input planets and asteroids, comments           
  character*9 numbin(nastx)
  character*60 compla,comast 
  character*100 colhep,colhea 
!  end of file                                                          
  logical eof 
! controls                                                              
  INTEGER inp,ibar 
! output of rdorb                                                       
  character*5 epoch 
  character*4 rsys 
! name of asteroid being computed: long. short                          
  CHARACTER*19 namid 
  CHARACTER*9 nam0 
! times                                                                 
  DOUBLE PRECISION t0,t0p,t0a,t0b,ta,tp,tt 
! definition flags                                                      
  LOGICAL cov0,defnor 
! magnitude, opposition effect, mass, matrices                          
  DOUBLE PRECISION hm,gma0,mass,g0(6,6),c0(6,6) 
! loop indexes                                                          
  INTEGER i,j,kr,nrecord,ll,jj,ii 
! number of bodies etc.                                                 
  INTEGER ninp,nast,npla,ilce,ilc,ilca,nang,nact,numang,numact 
! lce                                                                   
  DOUBLE PRECISION gamma 
! quality codes                                                         
  INTEGER iqco 
  SAVE 
! **********************************************************            
!  if this is the first call, read headers                              
!  (and barycentric correction if required);                            
!  do initialisations                                                   
  if(kr.eq.1)then 
! open input file for barycenter of inner solar system                  
     if(inp.eq.2)then 
        ibar=0 
     elseif(inp.eq.1)then 
        ibar=1 
     elseif(inp.eq.3)then 
        write(*,*)' wrong call to iosho with inp=',inp 
        stop 
     else 
        write(*,*)' unknown input type, inp=',inp 
        stop 
     endif
! *******************************************************               
!  elimination of short periodic terms is done in keplerian             
!  heliocentric elements, invariable plane                              
     cooy='KEP' 
     sysy='HEL' 
     refy='INVL1B' 
! **********************************************************            
!  barycenter correction, if required, from unit 17                     
! **********************************************************            
     ninp=0 
     if(ibar.gt.0)then 
        call reahea(iun17,nompla,compla,colhep,coox,sysx,refb,      &
     &      unitx,npla,ilc,t0b)                                         
        if(coox.ne.'CAR')then 
           write(*,*)'wrong baryc.file',coox 
           stop 
        endif
        call reaint(iun17,'ninp',ninp) 
        read(iun17,*)(pm(ii),ii=1,ninp) 
! read J2 and relativistic correction constants                         
        call skip(iun17,1) 
        read(iun17,*) cj2,scw 
        call skip(iun17,1) 
        read(iun17,*)barin 
     else 
!  do not worry, ninp is dummy but is needed for dimension.ne.0         
        ninp=1 
     endif
! **********************************************************            
!  planetary elements header, from unit 7:                              
! **********************************************************            
!  (1) header                                                           
     call reahea(iun7,nompla,compla,colhep,cooxp,sysxp,refxp,       &
     &     unitxp,npla,ilc,t0p)                                         
     nang=numang(cooxp) 
     nact=numact(cooxp) 
!  Jupiter and Saturn have to be available                              
     if(nompla(1).ne.'JUPITER  ')then 
        write(*,*)' Jupiter not available, first is',nompla(1) 
        stop 
     endif
     if(nompla(2).ne.'SATURN  ')then 
        write(*,*)' Saturn not available, second is',nompla(1) 
        stop 
     endif
!  check consistency of barycenter data (if required)                   
     if(ibar.eq.1)then 
        if(abs(t0b-t0p).gt.1.d-7)then 
           write(*,*)' refer.time for baryc. incompatible',t0b,t0p 
           stop 
        endif
        if(refb.ne.refxp)then 
           write(*,*)' ref. system of baryc.incompatible',refb,refx 
           stop 
        endif
     endif
     nbod=npla+1 
!  (2) input masses, sun=1 (warning: if baryc. corr. applied, mass      
!  of the sun is not 1 but 1+inner planets); all the mass ratios,       
!  gravitational constant, etc. are computed and stored in massmod
     call reamas(iun7,nbod,ibar,pm,ninp) 
     write(*,421)cooxp,sysxp,refxp,unitxp 
421  format(' input: planets, coox=',a3,                            &
     &          ' sysx=',a3,' refx=',a6,' unitx=',a3)                   
     write(*,*)' barycenter corr.flag ',ibar 
!c         call wrimas(9,nbod)                                          
! **********************************************************            
!  asteroids elements header, from unit 8:                              
! **********************************************************            
     if(inp.eq.1)then 
!  (1a) input header from catalogue in OEF format:                      
!   already done by oporbf                                              
     elseif(inp.eq.2)then 
!  (1b) input header from orbit9 format:                                
        call reahea(iun8,numbin,comast,colhea,coox,sysx,refx,       &
     &        unitx,nast,ilca,t0a)                                      
!  only one asteroid assumed in the file (output of conv9)              
        if(nast.gt.1)then 
           write(*,*)' nast=',nast,' use conv9 to select one' 
           stop 
        endif
        write(*,422)numbin(1),coox,sysx,refx,unitx,t0a 
422     format(' asteroid: ',a4,' input from orbit8v output'/       &
     &        ' coox=',a3,' sysx=',a3,' refx=',a6,' unitx=',a3/         &
     &        ' JD=',f13.4)                                             
!  (2) check consistency of reference times                             
        if(abs(t0a-t0p).gt.1.d-7)then 
           write(*,*)' reference times incompatible',t0a,t0p 
           stop 
        else 
           t0=t0a 
        endif
     endif
! *************************************************************         
!  output header                                                        
     if(inp.eq.1)then 
        write(iun19,232) 
232     format('  Mean elements Milani & Knezevic') 
     elseif(inp.eq.2)then 
        write(iun19,233)t0a,refy,numbin(1) 
233     format(' t0=',f10.1,'; ref=',a6,' asteroid ',a4,            &
     &   '  Mean elements Milani & Knezevic')                           
     else 
        write(*,*)' wrong input type in iosho, inp=',inp 
        stop 
     endif
! *************************************************************         
!  end of initialisations                                               
  endif
! *************************************************************         
!  begins input of data, to be performed at each call for the           
!  asteroids; only at the first call for planets if input is            
!  from catalogue                                                       
! *************************************************************         
  if(kr.eq.1.or.inp.eq.2)then 
! *************************************************************         
!  (3) input planetary elements from unit 7                             
     j=0 
     read(iun7,*,end=7,err=8)tp 
     if(inp.eq.1)then 
!  (3a) check consistency of input time: all at the reference epoch     
        if(abs(tp).gt.1.d-7)then 
           write(*,*)' initial times incompatible',tp 
           stop 
        endif
     endif
     DO j=1,nbod-1 
        read(iun7,*,end=9,err=8)(xp(i,j),i=1,nact)                 &
     &          ,(xp(i+nact,j),ngp(i,j),i=1,nang)                       
!c           write(9,*)tp,(xp(i,j),i=1,nact)                            
!c     +        ,(xp(i+nact,j),ngp(i,j),i=1,nang)                       
!  (4) conversion of angles to radiants                                 
        call radiant(xp(1,j),cooxp,unitxp,dpig) 
     ENDDO 
! *******************************************************               
!  coordinate change to elimination system                              
!  also barycentric correction, if required                             
     if(ibar.eq.0)then 
!  change to elimination coordinate system yp                           
        call coord(xp,sysxp,nbod,cooxp,refxp,yp,                    &
     &                sysy,cooy,refy,ennep,bar)                         
     else 
!  (5) barycenter correction: zzp is cartesian,                         
!   zbp has the bar. correction added                                   
        call coord(xp,sysxp,nbod,cooxp,refxp,zzp,                   &
     &                'HEL','CAR',refxp,ennep,bar)                      
        DO  j=1,nbod-1 
           zbp(1:6,j)=zzp(1:6,j)+barin(1:6) 
        ENDDO
!  change to elimination coordinate system yp                           
        call coord(zbp,'HEL',nbod,'CAR',refxp,yp,                   &
     &            sysy,cooy,refy,ennep,bar)                             
     endif
     pel(1:6,1)=yp(1:6,1) 
     pel(1:6,2)=yp(1:6,2)    
! quality code for orbit is assumed 0                                   
     iqco=0 
!c         write(9,*)pel                                                
  endif
! ***********************************************************           
!  input asteroid elements from OEF file                                
  if(inp.eq.1)then 
!   from catalogue:                                                     
!  (4a) input elements:                                                 
90   CONTINUE 
     CALL rdorb(namid,xa,coox,t0,g0,cov0,c0,defnor,                 &
     &           hm,gma0,mass,rsys,epoch,nrecord,eof)                   
     IF(eof)return 
! check reference system                                                
     IF(rsys.ne.'ECLM'.or.epoch.ne.'J2000')THEN 
        WRITE(*,*)' incompatible reference systems ' 
        WRITE(*,*)' ref=',rsys,' epoch=',epoch 
        WRITE(*,*)' should be ECLM, J2000' 
        STOP 
     ELSE 
        refx='ECLM00' 
        sysx='HEL' 
     ENDIF
!  if barycenter correction is required, it must all be consistent      
!  (the restriction on ref. syst. could be eliminated, see orbit8t      
     if(refx.ne.refb)then 
        write(*,*)' reference system of baryc. inconsistent' 
        write(*,*)' with asteroid input ',refb,refx 
        stop 
     endif
! handle funny identification names                                     
! convert name in case it is of the form nam0=namp                      
     ll=index(namid,'=')-1 
     IF(ll.lt.0)THEN 
        nam0=namid(1:9) 
     ELSE 
        nam0=namid(1:ll) 
     ENDIF
!90        call reaele(iun8,eof,no,xa(1,1),iq,iqco,iy,im,iday,hm)       
!  date in file; date has to be always the same as                      
!  t0p given in the header of the planets                               
     t0=t0+2400000.5d0 
     if(abs(t0-t0p).gt.0.4d0)then 
        write(*,*)namid,' wrong date ',t0 
        goto 90 
     endif
! change to elimination coordinate system ya                            
     na=1 
!  (8) adding the barycenter correction: zza is cartesian, zba is correc
     call cooast(xa,sysx,na,coox,refx,bar,nbod,zza,'HEL','CAR',     &
     &              refx,ennea)                                         
     do i=1,6 
        zba(i,1)=zza(i,1)+barin(i) 
     enddo
!  change to elimination coordinate system ya                           
     call cooast(zba,'HEL',na,'CAR',refx,bar,nbod,ya,sysy,cooy,     &
     &              refy,ennea)                                         
  elseif(inp.eq.2)then 
!  input from output file; it is assumed to contain only one            
!  asteroid (that is, output from conv8v)                               
!  (4b) input elements                                                  
     jj=1 
     if(ilca.gt.0)then 
        read(iun8,*,end=7,err=78)ta,(xa(i,jj),i=1,nact)             &
     &        ,(xa(i+nact,jj),nga(i,jj),i=1,nang),gamma                 
     else 
        read(iun8,*,end=7,err=78)ta,(xa(i,jj),i=1,nact)             &
     &        ,(xa(i+nact,jj),nga(i,jj),i=1,nang)                       
     endif
!  (3b) check consistency of reference times                            
     if(abs(tp-ta).gt.1.d-7)then 
        write(*,*)' record out of phase, planets ',tp,              &
     &       ' asteroids ',ta                                           
     endif
     na=1 
! (7) conversion of angles to radiants                                  
     call radiant(xa,coox,unitx,dpig) 
! *******************************************************               
! change to the elimination system;                                     
     call cooast(xa,sysx,na,coox,refx,bar,nbod,ya,sysy,             &
     &    cooy,refy,ennea)                                              
  endif
! ********************************************************************  
!  copy in the output arrays                                            
  tt=tp 
  el(1:6)=ya(1:6,1)       
!c     write(9,*)no,el                                                  
!  input of one record successfully completed                           
  eof=.false. 
  return 
! *************************************************************         
!  all the input error cases                                            
!  planet reading error                                                 
9 write(*,*)tp,' not enough planets, required ',npla,'  last is ',j 
  stop 
8 write(*,*)tp,' error in reading planet ',j,' required ',npla 
  stop 
!  asteroid reading error                                               
78 write(*,*)ta,' error in reading asteroid elements' 
  stop 
!  regular end of file (either planets or asteroids)                    
7 eof=.true. 
END SUBROUTINE iosho

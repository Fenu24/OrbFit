! **********************************************************            
!  REAPLA                                                               
!  planetary initial conditions, from unit 2:                           
! **********************************************************            
SUBROUTINE reapla(iunin,iunbar,iunou,ibar,iprqua,                 &
     &    coox,sysx,refx,sysz,refz,nbod,t0p,tp,xp,ngp,barin)            
  USE massmod
  USE fund_const
  IMPLICIT NONE 
! ====================INPUT==========================                   
! unit for input pla, input bar, for output, baricenter required, outp. 
  INTEGER, INTENT(IN) :: iunin,iunbar,iunou,ibar,iprqua 
! reference systems for output                                          
  CHARACTER*3, INTENT(IN) :: sysz 
  CHARACTER*6, INTENT(IN) :: refz 
! ====================OUTPUT=========================                   
! current number of bodies with mass (incl. the sun)                    
  INTEGER, INTENT(OUT) :: nbod 
! reference systems                                                     
  CHARACTER*3, INTENT(OUT) :: sysx,coox 
  CHARACTER*6, INTENT(OUT) :: refx 
! initial conditions for the planets                                    
  DOUBLE PRECISION, INTENT(OUT) :: xp(6,nbox) 
  INTEGER, INTENT(OUT) :: ngp(nangx,nbox) 
! epoch times: absolute, relative                                       
  DOUBLE PRECISION, INTENT(OUT) :: t0p,tp 
! barycenter of inner solar system                                      
  DOUBLE PRECISION, INTENT(OUT) :: barin(6) 
! ====================END INTERFACE==================                   
!  names and numbers of input planets and asteroids, comments           
  character*60 compla 
  character*100 colhep,colhea 
! inverse masses of the inner planets                                   
  DOUBLE PRECISION pm(10) 
! times, reference systems of baryc                                     
  DOUBLE PRECISION t0b 
  CHARACTER*6 refb 
  CHARACTER*3 unitx 
! dimensions                                                            
  INTEGER npla,nact,nang,ilc,ninp 
! functions                                                             
  INTEGER numang,numact 
! loop indexes                                                          
  INTEGER j,i,ii 
! =================BARYCENTER===========================                
!  barycenter correction, if required                                   
  ninp=0 
  if(ibar.gt.0)then 
     call reahea(iunbar,nompla,compla,colhep,coox,sysx,refb,unitx,  &
     &   npla,ilc,t0b)                                                  
     if(coox.ne.'CAR')then 
        write(9,*)'wrong baryc.file',coox 
        stop 
     endif
! read inner planets inverse masses (not used)                          
     call reaint(iunbar,'ninp',ninp) 
     read(iunbar,*)(pm(ii),ii=1,ninp) 
! read J2 and relativistic correction constants                         
     call skip(iunbar,1) 
     read(iunbar,*) cj2,scw 
     call skip(iunbar,1) 
! read center of mass of inner solar system                             
     read(iunbar,*)barin 
     close(iunbar) 
  else 
!  do not worry, ninp is dummy but is needed for dimension.ne.0         
     ninp=1 
  endif
! ================PLANETS========================                       
!  (1) input header                                                     
  call reahea(iunin,nompla,compla,colhep,coox,sysx,refx,unitx,      &
     &   npla,ilc,t0p)                                                  
!     WRITE(*,*)(nompla(j),j=1,npla)                                    
  if(ilc.ne.0)then 
     write(9,*)'LCE for planets not supported',ilc 
     stop 
  endif
  nbod=npla+1 
!  (2) input masses, sun=1                                              
  call reamas(iunin,nbod,ibar,pm,ninp) 
  write(9,*)'in input ',npla,' planets: ',coox,' ',sysx,' ',refx 
  write(9,409)ibar,ninp 
409 format('barycenter corr.flag ',i2,' no. inner pla. ',i3) 
  call wrimas(iunou,nbod) 
!  (3) input initial conditions                                         
  nang=numang(coox) 
  nact=numact(coox) 
  read(iunin,*,end=8,err=8)tp 
!     write(*,*)coox,nang,nact                                          
  write(iunou,*)' planetary initial conditions as read' 
  DO 2 j=1,nbod-1 
     read(iunin,*,end=9,err=8)(xp(i,j),i=1,nact)                     &
     &    ,(xp(i+nact,j),ngp(i,j),i=1,nang)                             
!  (4) conversion of angles to radiants                                 
     call radiant(xp(1,j),coox,unitx,dpig) 
2 ENDDO
  close(iunin) 
! check consistency of times, ref.system                                
  if(ibar.eq.1)then 
     if(abs(t0b-t0p).gt.1.d-7.or.tp.ne.0.d0)then 
        write(9,*)' refer.time for baryc. incompatible',t0b,t0p,tp 
        stop 
     endif
     if(refb.ne.refx)then 
        write(9,*)' ref. system of baryc.incompatible',refb,refx 
        stop 
     endif
  endif
! =======================output headers======================           
!  column headers for equinoctal elements and numbers                   
  colhea='  a(au)      h      k      p      q     lambda  nrev  gam' 
!  sampled elements, only for iprq even, unit 11 and 21                 
  if(mod(iprqua,2).eq.0)then 
     open(11,file='vpla.dat',status='unknown') 
     call wrihea(11,nompla,compla,colhea,'EQU',sysz,refz,'RAD',npla,0,t0p)
     call wrimas(11,nbod) 
  endif
!  filtered elements, always, unit 12 and 22                            
  open(12,file='vpla.fil',status='unknown') 
  call wrihea(12,nompla,compla,colhea,'EQU',sysz,refz,'RAD',npla,0,t0p)
  call wrimas(12,nbod) 
!  dump files                                                           
  open(13,file='orb9.dmp',status='unknown') 
  call wrihea(13,nompla,compla,colhea,'EQU',sysz,refz,'RAD',npla,0,t0p)
  call wrimas(13,nbod) 
  return 
!  planet not found                                                     
9 write(9,*)' not enough planets, required ',npla,'  last was ',j 
  stop 
8 write(9,*)' error in reading planet no. ',j,' required ',npla 
  stop 
END SUBROUTINE reapla

! ======================================================                
!  {\bf jobinp} input and extraction job definition                     
!  to be used by conv9                                                  
SUBROUTINE jobinp(inter,nast,imas,nbod,ifo,number,iast,nia,coox,  &
     &     sysx,refx,unitx,ilce,cooy,sysy,refy,uniy,iall,ifil,iorfl)    
  USE massmod   
  IMPLICIT NONE 
  INTEGER iast(nastx) 
  DOUBLE PRECISION dum(4) 
!  character arrays                                                     
  character*40 infile,outfile(nastx),plafil 
  character*3 coox,cooy,sysx,sysy,unitx,uniy,coox1,sysx1,unitx1 
  character*6 refx,refy,refx1 
  character*60 commen,compla,comme1 
  character*9 number(nastx),numbin(nastx),nmbs(nastx) 
  character*9 num,nnum,num1 
  character*6 postfi 
  character*100 colhea 
!  integers                                                             
  integer ifil, iorfl(nastx), iorflp, icnt,isa,ich,ilce0 
  INTEGER iall, ilce, nia, ifo, nbod,imas,nast,inter, ihead,iia,i,len,j,ihou
  INTEGER iconve, ilce1,npla,len1, iun, ill
  DOUBLE PRECISION t0, t0p
!  interactive/batch mode                                               
  write(*,*)' interactive=1 from conv9.inp=0' 
  read(*,*)inter 
!  open and input file description                                      
  if(inter.eq.1)then 
     write(*,*)' input  file ?' 
     read(*,140)infile 
140  format(a40) 
  else 
     open(4,file='conv9.inp',status='old') 
     read(4,140)infile 
     write(*,*)'input file=',infile 
  endif
  open(1,file=infile, status='old') 
  if(inter.eq.1)then 
     write(*,*)' header provided ? 1=yes 0=no' 
     read(*,*)ihead 
  else 
     read(4,*)ihead 
  endif
  if(ihead.gt.0)then 
     call reahea(1,numbin,commen,colhea,coox,                       &
     &      sysx,refx,unitx,nast,ilce,t0)                               
     write(*,*)'coox=',coox,' sysx=',sysx,' refx=',refx 
     write(*,*)'unitx=',unitx,' nast=',nast,' ilce=',ilce,' t0=',t0 
  else 
     write(*,*)' please supply header' 
     call askhea(numbin,commen,colhea,coox,sysx,refx,unitx,nast,ilce,t0) 
  endif
!  asteroid names in input                                              
  if(inter.eq.1)then 
!         do 45 iia=1,nast                                              
!           write(*,*)numbin(iia),' available'                          
! 45      continue                                                      
     write(*,*)' with mass ? 1=yes(planets) 0=no(asteroids)' 
     read(*,*)imas 
  else 
     read(4,*)imas 
     write(*,*)'imas=',imas 
  endif
  if(imas.eq.1)then 
     nbod=nast+1 
     call reamas(1,nbod,0,dum,1) 
  endif
!  asteroid names in input                                              
  DO iia=1,nast 
     write(*,*)numbin(iia),' available' 
  ENDDO
!  choice between two formats: time,record if one orbit (anyway with    
!  free format reading the newline does not matter),                    
!  time/records for more orbits                                         
  if(nast.gt.1)then 
     ifo=1 
  elseif(nast.eq.1)then 
     ifo=0 
  elseif(nast.lt.1.and.nast.gt.nastx)then 
     write(*,*)' no. orbits ',nast,' not allowed' 
     stop 
  endif
!  output selection                                                     
  iia=1 
43 continue 
  if(inter.eq.1)then 
     write(*,*)' which body among the above ?' 
     read(*,143)num 
  else 
     read(4,143)num 
  endif
143 format(a9) 
  if(num.eq.'all')then 
     DO i=1,nast 
        number(i)=numbin(i) 
        iast(i)=i 
     ENDDO
     if(inter.eq.0) read(4,*) 
     iia=nast+1 
     goto 47 
  endif
  num1=num 
  call cmpstr(num1,len) 
  if(len.le.0)goto 47 
  do 46 i=1,nast 
     nnum=numbin(i) 
     if(num.eq.nnum)then 
        do 48 j=1,iia 
           nnum=number(j) 
           if(num.eq.nnum)then 
              write(*,*)num,' already selected' 
              goto 43 
           endif
48      ENDDO
        write(*,*)num,' OK' 
        iast(iia)=i 
        number(iia)=num 
        iia=iia+1 
        goto 43 
     endif
46 ENDDO
  write(*,*)num,' not available' 
  goto 43 
!  files opening and headers                                            
47 nia=iia-1 
  if(inter.eq.1)then 
     write(*,*)' postfix for output files ?' 
     read(*,145)postfi 
     write(*,*)' header in output files? 1=yes 0=no' 
     read(*,*)ihou 
!  choice of output system                                              
     write(*,*)' conversion? 1=yes 0=no' 
     read(*,*)iconve 
     if(iconve.eq.0)then 
        cooy=coox 
        sysy=sysx 
        refy=refx 
        uniy=unitx 
     else 
        call asksys(cooy,sysy,refy,uniy) 
     endif
  else 
     read(4,145)postfi 
     read(4,*)ihou 
     write(*,*)'ihou=',ihou,' postfix for output=',postfi 
     read(4,144)cooy 
     write(*,*)'cooy=',cooy 
     read(4,144)sysy 
     write(*,*)'sysy=',sysy 
     read(4,145)refy 
     write(*,*)'refy=',refy 
     read(4,144)uniy 
     write(*,*)'uniy=',uniy 
144  format(a3) 
145  format(a6) 
  endif
!  planets needed for coordinate transformation?                        
  if(imas.eq.0)then 
     if(sysx.eq.sysy.and.coox.eq.cooy.and.refx.eq.refy)then 
!  no coordinate transformation is really used                          
        goto 146 
     else 
        if(inter.eq.1)then 
           write(*,*)' file with planets? (must have headers)' 
           read(*,140)plafil 
        else 
           read(4,140)plafil 
           write(*,*)'file with planets=',plafil 
        endif
        open(2,file=plafil, status='old') 
        call reahea(2,nompla,compla,colhea,coox1,                      &
     &        sysx1,refx1,unitx1,npla,ilce1,t0p)                        
        nbod=npla+1 
        call reamas(2,nbod,0,dum,1) 
!  this restriction could be removed if needed                          
!  (but in practice the output of orbit9 is always the same             
!  for planets and asteroids)                                           
        if(coox1.ne.coox.or.sysx1.ne.sysx.or.refx1.ne.refx.or.unitx1.ne.unitx)then                                              
           write(*,*)' incompatible coord. systems ' 
           stop 
        elseif(abs(t0-t0p).gt.1.d-3)then 
           write(*,*)' incompatible times t0=',t0,'  t0p=',t0p 
           stop 
        endif
     endif
146  continue 
  endif
  if(ihou.gt.0)then 
     call cmpstr(commen,len) 
     comme1='from file '//infile//' ' 
     call cmpstr(comme1,len1) 
     comme1=comme1(1:len1)//commen(1:len) 
  endif
  if(cooy.eq.'EQU')then 
     colhea='  a(au)   h    k    p   q   lambda  nrev gamma' 
  elseif(cooy.eq.'KEP')then 
     colhea='  a(au)   e    I    Omega  omega  Mean.an.  nrev gamma' 
  elseif(cooy.eq.'CAR')then 
     colhea='   x     y     z     vx     vy     vz    gamma' 
  else 
     write(*,*)' output coordinates ',cooy,'  unknown' 
     stop 
  endif
  if(inter.eq.1)then 
     write(*,*)' 0=one orb/file; 1=more orbs/file; 2=all together' 
     read(*,*)iall 
  else 
     read(4,*)iall 
     write(*,*)'iall= ',iall 
  endif
! file(s) opening and headers                                           
! one orbit per file                                                    
  if(iall.eq.0)then 
     DO 49 iia=1,nia 
        outfile(iia)='v'//number(iia)//'.'//postfi 
        call rmbl(outfile(iia),len) 
        iun=10+iia 
        open(iun,file=outfile(iia), status='unknown') 
        if(ilce.ge.iia)then 
           ill=1 
        else 
           ill=0 
        endif
        if(ihou.gt.0)then 
           call wrihea(iun,number(iia),comme1,colhea,cooy,           &
     &          sysy,refy,uniy,1,ill,t0)                                
        endif
49   ENDDO
! more orbits per file                                                  
  elseif(iall.eq.1)then 
     icnt=0 
     ifil=1 
71   if(inter.eq.1)then 
        write(*,*)' How many orbits in file',ifil,' ?' 
        read(*,*)iorfl(ifil) 
     else 
        read(4,*)iorfl(ifil) 
     endif
! REGULAR EXIT                                                          
     if(iorfl(ifil).eq.0)return 
     write(*,*)'File',ifil,' with',iorfl(ifil),' bodies' 
     iorflp=icnt+iorfl(ifil) 
! IRREGULAR EXIT - this should not happen                               
     if(iorflp.gt.nia)then 
        if(inter.eq.1)then 
           write(*,*)' More than',nia,' bodies selected.' 
           write(*,*)' Repeat the last entry!' 
           goto 71 
        else 
           write(*,*)' More than',nia,' bodies selected.' 
           write(*,*)' Check the entries for no. of orb/file!' 
           stop 
        endif
     endif
     outfile(ifil)='vm'//achar(48+ifil)//'.'//postfi 
     call rmbl(outfile(ifil),len) 
     iun=10+ifil 
     open(iun,file=outfile(ifil), status='unknown') 
     if(ilce.eq.0)then 
        ill=0 
     elseif(ilce.lt.iorflp)then 
        ill=ilce-icnt 
     else 
        ill=iorfl(ifil) 
     endif
     if(ihou.gt.0)then 
        DO isa=1,iorfl(ifil) 
           nmbs(isa)=number(icnt+isa) 
        ENDDO
        call wrihea(iun,nmbs,comme1,colhea,cooy,sysy,refy,uniy,     &
     &                iorfl(ifil),ill,t0)                               
     endif
     icnt=icnt+iorfl(ifil) 
     ifil=ifil+1 
     goto 71 
! all together                                                          
  else 
     outfile(1)='vall'//'.'//postfi 
     call rmbl(outfile(1),len) 
     iun=11 
     open(iun,file=outfile(1), status='unknown') 
     if(ihou.gt.0)then 
        if(ilce.gt.nia)ilce=nia 
        call wrihea(iun,number,comme1,colhea,cooy,sysy,refy,uniy,nia,ilce,t0)
        if(imas.eq.1)call wrimas(iun,nia+1) 
     endif
  endif
END SUBROUTINE jobinp

! ========================================================              
!  {\bf reahea} ORB8V CONV8V GIFFV read header of data file on unit iun 
subroutine reahea(iun,numbin,commen,colhea,coo,sys,ref,           &
     &           uni,nast,ilc,t0)                                       
  IMPLICIT NONE 
  double precision t0 
  character*3 coo,sys,uni 
  character*6 ref 
  integer iun,nast,ilc 
  character*9 numbin(*) 
  character*70 nome 
  character*100 colhea 
  character*60 commen 
  INTEGER j,nr,nn 
!  read name of the object, data on integration (first line)            
  read(iun,100)nast,nome 
100 format(i4,a70) 
!     rewind iun                                                        
  IF(nast.lt.7)THEN 
     read(nome,103)(numbin(j),j=1,nast) 
103  FORMAT(7(1x,a9)) 
  ELSE 
     read(nome,103)(numbin(j),j=1,7) 
  ENDIF
!  read number and names of the objects,following lines                 
  IF(nast.gt.7)THEN 
     READ(iun,101)(numbin(j),j=8,nast) 
101  FORMAT((3x,7(1x,a9))) 
  ENDIF
!      WRITE(*,101)(numbin(j),j=1,nast)                                 
!  reference time and comment                                           
  call reaflc(iun,'t0',t0,commen) 
!  read definition of reference system,no. LCE, etc (second line)       
  read(iun,110)coo,sys,ref,uni,ilc 
110 format(8x,a3,1x,a3,1x,a6,1x,a3,7x,i4) 
!  read column headers                                                  
  read(iun,102)colhea 
102 format(a100) 
END subroutine reahea
! ========================================================              
!  {\bf wrihea} ORB8V CONV8V write header of data file on unit iun      
subroutine wrihea(iun,numbin,commen,colhea,                       &
     &      coo,sys,ref,uni,nast,ilc,t0)                                
  IMPLICIT NONE 
  double precision t0 
  character*3 coo,sys,uni 
  character*6 ref 
  integer iun,ilc,nast 
  character*9 numbin(*) 
  character*100 colhea 
  character*60 commen 
  INTEGER j 
!  write number and names of the objects                                
  IF(nast.le.7)THEN 
     WRITE(iun,100)nast,(numbin(j),j=1,nast) 
100  FORMAT(i4,7(1x,a9)) 
  ELSE 
     WRITE(iun,101)nast,(numbin(j),j=1,nast) 
101  FORMAT(i4,7(1x,a9)/(3x,7(1x,a9))) 
  ENDIF
!  reference time and comments                                          
  call wriflo(iun,'t0',t0,commen) 
!  write definition of reference system, no. LCE, etc (second line)     
  write(iun,110)coo,sys,ref,uni,ilc 
110 format('coosys=''',a3,1x,a3,1x,a6,1x,a3,''';ilce=',i4) 
!  write column headers                                                 
  write(iun,102)colhea 
102 format(a100) 
END subroutine wrihea
! ======================================================================
!   {\bf askhea} CONV8V GIFFV interactive input of header               
subroutine askhea(numbin,commen,colhea,coox,                      &
     &         sysx,refx,unitx,nast,ilce,t0)                            
  USE parameters_giff
  INTEGER, PARAMETER :: nastx=nplx 
  double precision t0 
  character*3 coox,sysx,unitx 
  character*6 refx 
  character*60 commen 
  character*9 numbin(nastx) 
  character*100 colhea 
!  number and name of bodies in input file                              
  write(*,*)' number of orbits in input file' 
  read(*,*)nast 
  if(nast.gt.nastx)then 
     write(*,*)' too many: ',nast,' > ',nastx 
     stop 
  endif
  write(*,*)' give ',nast,' identifiers, 1 per line, format a4' 
  do  n=1,nast 
     read(*,100)numbin(n)
  enddo
100 format(a9) 
  write(*,*)' out of ',nast,'  orbits, how many with max LCE?' 
  read(*,*)ilce 
  write(*,*)' reference time in JD (e.g. 2448600.5 12 Dec. 91) ?' 
  read(*,*)t0 
  write(*,*)' comment for output files?' 
  read(*,*)commen 
!  coordinate and reference system                                      
  write(*,*)' coordinates, EQU, KEP, CAR ?' 
  read(*,101)coox 
101 format(a3) 
  write(*,*)' coordinate system HEL, JAC, BAR, HEC ?' 
  read(*,101)sysx 
  write(*,*)' reference sys. ECLM50,EQUM00,EQUM50,INVL1B,ECLM00 ?' 
  read(*,102)refx 
102 format(a5) 
  write(*,*)' units for angles, RAD DEG REV ?' 
  read(*,101)unitx 
  write(*,*)' column headers line ?' 
  read(*,*)colhea 
END subroutine askhea
! ============================================================          
!  {\bf asksys}  CONV8V                                                 
!  ask for output coordinates and reference                             
subroutine asksys(cooy,sysy,refy,uniy) 
  character*3 cooy,sysy,uniy 
  character*6 refy 
  write(*,*)' coordinates, EQU, KEP, CAR ?' 
  read(*,101)cooy 
101 format(a3) 
  write(*,*)' coordinate system HEL, JAC, BAR, HEC ?' 
  read(*,101)sysy 
  write(*,*)' reference sys,ECLM50,EQUM50,INVL1B,ECLM00,EQUM00 ?' 
  read(*,102)refy 
102 format(a6) 
  write(*,*)' units for angles RAD, DEG, REV ?' 
  read(*,101)uniy 
END subroutine asksys
! ====================================================                  
!  {\bf radiant} ORB8V converts angles of elements in radiants          
subroutine radiant(x,coo,unit,dpig) 
  IMPLICIT NONE 
  DOUBLE PRECISION, INTENT(INOUT) :: x(6) 
  DOUBLE PRECISION, INTENT(IN) :: dpig
  CHARACTER*3, INTENT(IN) ::  unit,coo
! END INTERFACE
  DOUBLE PRECISION undeg, unday 
  INTEGER i
  if(unit.eq.'DEG')then 
     undeg=dpig/3.6d2 
  elseif(unit.eq.'REV')then 
     undeg=dpig 
  elseif(unit.eq.'RAD')then 
     undeg=1.d0 
  elseif(unit.eq.'DAY')then 
     if(coo.eq.'CAR')then 
        unday=3.6525d2 
     else 
        write(*,*)' unit DAY used only for coord. CAR' 
        stop 
     endif
  else 
     write(9,*)' units for angles ',unit,'  not supported' 
     stop 
  endif
!  warning: do not use numang: in KEP, inclination is                   
!  an angle but not an angle variable                                   
  if(coo.eq.'KEP')then 
     do  i=3,6 
        x(i)=x(i)*undeg
     enddo 
  elseif(coo.eq.'EQU'.or.coo.eq.'EQP')then 
     x(6)=x(6)*undeg 
  elseif(coo.eq.'CAR')then 
     if(unit.eq.'DAY')then 
        do  i=4,6 
           x(i)=x(i)*unday 
        enddo
     endif
  endif
END subroutine radiant
! ====================================================                  
!  {\bf radang} CONV8V converts angles of elements in radiants,         
!   vector version                                                      
subroutine radang(x,n,coo,unit,dpig) 
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: n
  DOUBLE PRECISION, INTENT(INOUT) :: x(6,n) 
  DOUBLE PRECISION, INTENT(IN) :: dpig
  CHARACTER*3, INTENT(IN) ::  unit,coo
! END INTERFACE
  DOUBLE PRECISION undeg, unday 
  INTEGER j,i
  if(unit.eq.'DEG')then 
     undeg=dpig/3.6d2 
  elseif(unit.eq.'REV')then 
     undeg=dpig 
  elseif(unit.eq.'RAD')then 
     undeg=1.d0 
  elseif(unit.eq.'DAY')then 
     if(coo.eq.'CAR')then 
        unday=3.6525d2 
     else 
        write(*,*)' unit DAY used only for coord. CAR' 
        stop 
     endif
  else 
     write(9,*)' units for angles ',unit,'  not supported' 
     stop 
  endif
  if(coo.eq.'KEP')then 
     do  j=1,n 
        do i=3,6 
           x(i,j)=x(i,j)*undeg
        enddo
     enddo 
  elseif(coo.eq.'EQU'.or.coo.eq.'EQP')then 
     do  j=1,n 
        x(6,j)=x(6,j)*undeg
     enddo
  elseif(coo.eq.'CAR')then 
     if(unit.eq.'DAY')then 
        do  j=1,n 
           do i=4,6 
              x(i,j)=x(i,j)*unday
           enddo
        enddo
     endif
  endif
END subroutine radang
! ====================================================                  
!  {\bf unrad} CONV8V converts angles of elements from radiants,        
!   vector version                                                      
subroutine unrad(x,n,coo,unit,dpig) 
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: n
  DOUBLE PRECISION, INTENT(INOUT) :: x(6,n) 
  DOUBLE PRECISION, INTENT(IN) :: dpig
  CHARACTER*3, INTENT(IN) ::  unit,coo
! END INTERFACE
  DOUBLE PRECISION undeg, unday 
  INTEGER j,i
  if(unit.eq.'DEG')then 
     undeg=dpig/3.6d2 
  elseif(unit.eq.'REV')then 
     undeg=dpig 
  elseif(unit.eq.'RAD')then 
     undeg=1.d0 
  elseif(unit.eq.'DAY')then 
     if(coo.eq.'CAR')then 
        unday=3.6525d2 
     else 
        write(*,*)' unit DAY used only for coord. CAR' 
        stop 
     endif
  else 
     write(9,*)' units for angles ',unit,'  not supported' 
     stop 
  endif
  if(coo.eq.'KEP')then 
     do  j=1,n 
        do i=3,6 
           x(i,j)=x(i,j)/undeg 
        enddo
     enddo
  elseif(coo.eq.'EQU')then 
     do  j=1,n 
        x(6,j)=x(6,j)/undeg
     enddo
  elseif(coo.eq.'CAR')then 
     if(unit.eq.'DAY')then 
        do  j=1,n 
           do  i=4,6 
              x(i,j)=x(i,j)/unday
           enddo
        enddo
     endif
  endif
END subroutine unrad
! ==============================================                        
!  {\bf numang} ORB8V CONV8V GIFFV                                      
!  number of angle variables in coordinate systems                      
integer function numang(coo) 
  character*3 coo 
  if(coo.eq.'KEP')then 
     numang=3 
  elseif(coo.eq.'EQU')then 
     numang=1 
  elseif(coo.eq.'CAR')then 
     numang=0 
  else 
     write(9,*)' numang: this should not happen ',coo 
     stop 
  endif
END function numang
! ==============================================                        
!  {\bf numact} ORB8V CONV8V GIFFV                                      
!   number of real variables in coordinate systems                      
integer function numact(coo) 
  character*3 coo 
  if(coo.eq.'KEP')then 
     numact=3 
  elseif(coo.eq.'EQU')then 
     numact=5 
  elseif(coo.eq.'CAR')then 
     numact=6 
  else 
     write(9,*)' numact: this should not happen ',coo 
     stop 
  endif
END function numact

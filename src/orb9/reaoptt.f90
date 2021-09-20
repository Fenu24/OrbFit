! **********************************************************            
!  {\bf reaopt} input options and open files                            
! **********************************************************            
SUBROUTINE reaoptt(iun,ibar,catalo,iprqua,                        &
     &  dt,nout,idump,nsamp,nsamp2,njump,njump2,                        &
     &  sysz,refz,v1,semim) 
  USE fund_const
  USE massmod   
  USE dyn_param, ONLY: iyark                                         
  IMPLICIT NONE 
! input                                                                 
  INTEGER, INTENT(IN) :: iun 
! output                                                                
  INTEGER, INTENT(OUT) :: ibar,idump,nsamp,nsamp2,njump,njump2,nout,iprqua 
  CHARACTER, INTENT(OUT) :: sysz*3,refz*6 
  DOUBLE PRECISION, INTENT(OUT) :: dt,v1,semim 
! end interface                                                         
!  file names                                                           
  CHARACTER inplan*72,inast*72,inbar*72 
! record with asteroid name                                             
  CHARACTER*20 reco 
  INTEGER len,lplan,lbar,last,le 
!  flag for catalog input/input from output                             
  LOGICAL catalo 
! dimensions, units, loops                                              
  INTEGER j,jj,k 
  INCLUDE 'comnbo.h90' 
  INCLUDE 'libcen.h90' 
! asteroid names for libration center                                   
  CHARACTER*9 idac(nastx) 
  DOUBLE PRECISION dcent(nastx) 
! **********************************************************            
!  skip record                                                          
  CALL skip(iun,1) 
!  input files:                                                         
!  planetary initial conditions                                         
  call skip(iun,1) 
  call reastr(iun,'inplan',inplan) 
  call rmbl(inplan,lplan) 
  open(2,file=inplan(1:lplan),status='old') 
!  barycenter correction, if required                                   
  call reaint(iun,'ibar',ibar) 
  if(ibar.gt.0)then 
     call reastr(iun,'inbar',inbar) 
     call rmbl(inbar,lbar) 
     open(7,file=inbar(1:lbar),status='old') 
  else 
     call skip(iun,1) 
  endif
!  input file: asteroid initial conditions                              
  call reastr(iun,'inast',inast) 
  call rmbl(inast,last) 
  CALL oporbf(inast(1:last), -1) 
!  two cases need to be separated: input from catalogue                 
!  (with  file name containing the string 'oscel') and                  
!  input from output of a previous computation                          
  if((index(inast,'.dat').ne.0.and.inast.ne.'astorb.dat')           &
     &     .or.index(inast,'.dma').ne.0)then                            
!  input from output of numerical integration                           
     catalo=.false. 
  else 
!  input from OEF format for orbit catalog                              
     catalo=.true. 
  endif
!  variational equations included?                                      
  call skip(iun,1) 
  call reaint(iun,'nvz',nvz) 
!  input list of asteroids required                                     
  na=0 
  DO 3 j=1,nastx 
     read(iun,100,end=4)reco 
100  format(a20) 
     call rmbl(reco,len) 
     if(len.eq.0)goto 30 
     ida(j)=reco(1:len) 
     na=na+1 
3 ENDDO
!  error in input of the options                                        
  write(9,*)na,' too many asteroids at once, max was', nastx 
  stop 
!  control on no. var.eq.                                               
30 if(nvz.gt.na)nvz=na 
  if(nvz.gt.nvzx)then 
     write(9,*)nvz,' too many variational equations, max was ',nvzx 
     stop 
  endif
  if(nvz.ne.0.and.nvz.lt.na)then 
     write(9,*)' warning: nvz =',nvz,' less than na=',na 
     write(9,*)' not allowed in all versions; in this case' 
     write(9,*)' nvz has been forced to be equal to na' 
     nvz=na 
  endif
! input of libration center data                                        
! note: the file libcenter.dat contains the absolute value of dcen      
! the sign must be adjusted taking into account that in L4 (lambda-lambd
! dcen >0, in l5 dcen<0 (see in orbit9t.f)                              
  open(99,file='libcenter.dat',status='unknown') 
  jj=0 
  DO j=1,nastx 
     READ(99,*,end=5) idac(j),dcent(j) 
     CALL rmsp(idac(j),le) 
     WRITE(*,*) idac(j),dcent(j) 
     jj=jj+1 
  ENDDO
5 CONTINUE 
  DO j=1,na 
     dcen(j)=0.d0 
     DO k=1,jj 
        IF(idac(k).eq.ida(j))THEN 
           dcen(j)=dcent(k)*radeg 
        ENDIF
     ENDDO
  ENDDO
  close (99) 
! Yarkovsky effect not used for Trojans
  iyark=0
!  output  control                                                      
  call skip(iun,1) 
  call reaflo(iun,'dt',dt) 
  call reaint(iun,'nout',nout) 
  call reaint(iun,'idump',idump) 
  call reaint(iun,'nsamp',nsamp) 
  call reaint(iun,'nsamp2',nsamp2) 
  call reaint(iun,'njump',njump) 
  njump2=njump*nsamp/nsamp2 
  if(njump*nsamp.ne.njump2*nsamp2)then 
     write(9,*)' bad output synchro.',njump,nsamp,njump2,nsamp2 
     stop 
  else 
     write(9,*)' nsamp=',nsamp,' njump=',njump,' nsamp2=',nsamp2,   &
     &     ' njump2=',njump2                                            
  endif
  call reaint(iun,'iprqua',iprqua) 
  call reastr(iun,'sysz',sysz) 
  call reastr(iun,'refz',refz) 
!  options for computation of Lyapounov exponents                       
  call skip(iun,1) 
  call reaflo(iun,'v1',v1) 
  call reaflo(iun,'semim',semim) 
  RETURN 
!  option file incomplete                                               
4 write(9,*)' error in options file ' 
  STOP
END SUBROUTINE reaoptt

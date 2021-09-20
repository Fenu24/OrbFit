! **********************************************************            
!  {\bf reaopt} input options and open files                            
! **********************************************************            
SUBROUTINE reaopt(iun,ibar,catalo,                                &
     &  dt,nout,idump,nsamp,iprqua,                                     &
     &  sysz,refz,v1,semim)                                             
  USE massmod
  IMPLICIT NONE 
! input                                                                 
  INTEGER, INTENT(IN) ::  iun 
! output                                                                
  INTEGER, INTENT(OUT) :: ibar,idump,nsamp,iprqua,nout 
  CHARACTER, INTENT(OUT) :: sysz*3,refz*6 
  DOUBLE PRECISION, INTENT(OUT) :: dt,v1,semim 
! end interface                                                         
!  file names                                                           
  character inplan*72,inast*72,inbar*72 
! record with asteroid name                                             
  character*20 reco 
  integer len,lplan,lbar,last 
!  flag for catalog input/input from output                             
  logical catalo 
! dimensions, units, loops                                              
  integer j 
! commons                                                               
  INCLUDE 'comnbo.h90' 
! **********************************************************            
!  skip record                                                          
  call skip(iun,1) 
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
  CALL oporbf(inast(1:last),-1) 
!  two cases need to be separated: input from catalogue                 
!  (with  file name containing the string 'oscel') and                  
!  input from output of a previous computation                          
  if((index(inast,'.dat').ne.0.and.inast.ne.'astorb.dat')   &
&           .or.index(inast,'.dma').ne.0)then                            
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
!  output  control                                                      
  call skip(iun,1) 
  call reaflo(iun,'dt',dt) 
  call reaint(iun,'nout',nout) 
  call reaint(iun,'idump',idump) 
  call reaint(iun,'nsamp',nsamp) 
  call reaint(iun,'iprqua',iprqua) 
  call reastr(iun,'sysz',sysz) 
  call reastr(iun,'refz',refz) 
!  options for computation of Lyapounov exponents                       
  call skip(iun,1) 
  call reaflo(iun,'v1',v1) 
  call reaflo(iun,'semim',semim) 
  return 
!  option file incomplete                                               
4 write(9,*)' error in options file ' 
  stop 
END SUBROUTINE reaopt

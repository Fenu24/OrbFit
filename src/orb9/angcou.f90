! ==============================================                        
!  {\bf numrev} preserves revolution counters across a change of        
!  coordinate type and system                                           
subroutine numrev(y,nast,x,ng,yold,ngold,coox,cooy) 
   USE fund_const
   IMPLICIT NONE
   INTEGER, parameter :: nelex=6,nangx=3 
   INTEGER :: nast
   DOUBLE PRECISION ::  y(nelex,nast),x(nelex,nast),yold(nelex,nast) 
   INTEGER ::  ng(nangx,nast),ngold(nangx,nast) 
   character*3, INTENT(IN) ::  coox,cooy 
! end interface
!  same coordinates: angles should never change by pi                   
!  input: x,ng old coordinates coox                                     
!         y new coordinates cooy=coox                                   
!  output: ng number rev. in new coordinates                            
!  warning: in this case x would be spoiled, has to be restored         
   INTEGER nang, n, j, i, numang, iround, k, nangj
   DOUBLE PRECISION tmp, dl
   if(coox.eq.cooy)then 
      nang=numang(coox) 
      DO 1 n=1,nast 
         DO  j=1,nang 
            tmp=x(7-j,n) 
            call prngup(y(7-j,n),x(7-j,n),ng(nang+1-j,n)) 
            x(7-j,n)=tmp 
            ngold(nang+1-j,n)=ng(nang+1-j,n) 
            yold(7-j,n)=y(7-j,n) 
         ENDDO
1     ENDDO
      return 
!  from keplerian to equinoctal: mean longitude is the sum              
!  of the three keplerian angles, including number of rev.              
!  input: x,ng keplerian coord; y equinoctal coord                      
!  output: ng equinoctal no. rev (one angle only)                       
   elseif(coox.eq.'KEP'.and.cooy.eq.'EQU')then 
      do 2 n=1,nast 
         dl=x(6,n)+x(5,n)+x(4,n)-y(6,n) 
         k=iround(dl/dpig) 
         ng(1,n)=ng(1,n)+ng(2,n)+ng(3,n)+k 
2     ENDDO
      return 
!  from equinoctal to keplerian: the assumption is that                 
!  the argument of perihelion and the longitude of node do              
!  not perform more than half a revolution between two outputs          
!  input: x,ng equinoctal coord                                         
!         yold,ngold keplerian coordinates of previous output           
!  output: ng=ngold keplerian no rev (3 angles)                         
   elseif(coox.eq.'EQU'.and.cooy.eq.'KEP')then 
      DO 3 n=1,nast 
         dl=y(6,n)+y(5,n)+y(4,n)-x(6,n) 
         k=iround(dl/dpig) 
!  update slow angles                                                   
         call prngup(y(4,n),yold(4,n),ngold(1,n)) 
         call prngup(y(5,n),yold(5,n),ngold(2,n)) 
!  for mean anomaly                                                     
         ngold(3,n)=ng(1,n)-k-ngold(1,n)-ngold(2,n) 
         DO i=1,3 
            ng(i,n)=ngold(i,n) 
         ENDDO
3     ENDDO
      return 
!  if output is cartesian, nothing to do                                
   elseif(cooy.eq.'CAR')then 
      return 
!  if input is cartesian, no way to do it                               
   elseif(coox.eq.'CAR')then 
      nang=numang(cooy) 
      ng(1:nangj,1:nast)=0 
      return 
   else 
      write(9,*)' numrev: this should not happen',coox,cooy 
   endif
   return 
 END subroutine numrev
! **********************************************************            
!  {\bf predicl} number of revolutions: prediction of likely value      
 subroutine predicl(pl,ngi,tout,enne,ngip) 
   USE fund_const
   IMPLICIT NONE
   DOUBLE PRECISION :: pl, tout, enne
   INTEGER ngi, ngip, iround
! end interface
   INTEGER ndl
   DOUBLE PRECISION :: dl
   dl=enne*tout 
   ndl=iround(dl/dpig) 
   dl=dl-dpig*ndl 
   pl=pl+dl 
   ngip=ngi+ndl 
   if(pl.gt.dpig)then 
      ngip=ngip+1 
      pl=pl-dpig 
   elseif(pl.lt.0.d0)then 
      ngip=ngip-1 
      pl=pl+dpig 
   endif
   return 
 END subroutine predicl
! **********************************************************            
!  {\bf iround} rounding of a real to an integer, in an even way        
 integer function iround(x) 
   double precision x,d 
   i=x 
   d=x-i 
   if(d.gt.0.5d0)i=i+1 
   if(d.le.-0.5d0)i=i-1 
   iround=i 
 END function iround
! **********************************************************            
!  {\bf prngup}                                                         
!   given an angle ang and an angle vang with ng revolutions            
!   finds a new new ng in the assumption that                           
!   less than half a revolution has occurred; the new ang               
!   replaces vang in the output.                                        
!   in case ang is not between 0 and dpig (eg between                   
!   -pig and pig, as in the output of datan2), it is reduced            
!   to principal value                                                  
 subroutine prngup(ang,vang,ng) 
   USE fund_const
   IMPLICIT NONE
   DOUBLE PRECISION :: ang, vang
   INTEGER :: ng
! end interface
   DOUBLE PRECISION :: d
   if(ang.lt.0.d0)ang=ang+dpig 
   d=ang-vang 
   if(d.gt.dpig*0.5d0)ng=ng-1 
   if(d.lt.-dpig*0.5d0)ng=ng+1 
   vang=ang 
  END subroutine prngup

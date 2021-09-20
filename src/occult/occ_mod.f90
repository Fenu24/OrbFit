! module with routines and functions for star occultation

MODULE occ_mod



IMPLICIT NONE



DOUBLE PRECISION ::  occ_a,occ_b,occ_c,occ_d, occ_v, occ_t0, occ_dt
! occulting shape is ellipse $a x^2+2b x(y-d)+ c(y-d)^2\le 1$ 
! center of occulting body passes from x=0 at time t0, then moves with 
! velocity v; sampling of time with step dt

! limb darkening parameters: limb_dark=1-limb1*(1-cosf) +limb2*(1-cosf)**2
DOUBLE PRECISION :: limb1, limb2

INTEGER, PARAMETER :: nsam_x= 1000000 ! maximum sample tested 

INTEGER, PARAMETER :: ntim_x= 10000 ! maximum number of times

DOUBLE PRECISION, DIMENSION(nsam_x):: x_rect,y_rect 
! sampling rectangle containing occulter at t0

DOUBLE PRECISION, DIMENSION(nsam_x):: x_occ, y_occ
! sampling occulter shape at t0

DOUBLE PRECISION, DIMENSION(ntim_x):: t_occ, rel_lum

!functions/routines sample_occult, belong_occult, belong_star


CONTAINS

!selects rectangle containing the occulting cross section
! given the shape defined by public module data
! also computed time interval of occultation and eigenvalues of normal matrix
 SUBROUTINE rectangle_time_select(xmin,xmax,ymin,ymax,area,tmin,tmax,eig)
  USE fund_const

  DOUBLE PRECISION, INTENT(OUT):: xmin,xmax,ymin,ymax,area,tmin,tmax, eig(2)
  
  DOUBLE PRECISION det, tra, angle, stdx,stdy
! invariants
  det=occ_a*occ_c-occ_b**2
  tra=occ_a+occ_c
  IF(det.lt.0.d0)THEN
     WRITE(*,*) ' negative determinant, not an ellipse: a,b,c,det,tra ='
     WRITE(*,*)      occ_a,' ',occ_b,' ',occ_c,' ',tra,' ',det
     STOP
  ENDIF
! eigenvalues (to be output)
  eig(1)=(tra+sqrt(tra**2-4*det))/2
  eig(2)=(tra-sqrt(tra**2-4*det))/2
! note that det is also the product of the two eigenvalues, thus pig/sqrt(det)
! = pig/sqrt(eig(1)*eig(2)) is the area of the ellipse; we omit pig
  area=pig/sqrt(det)
  angle=0 ! to be done
! marginal standard deviations, from covariance matrix  
  stdx=sqrt(occ_c/det)
  stdy=sqrt(occ_a/det)
! rectangle
  xmin=-stdx; xmax=stdx
  ymin=occ_d-stdy;ymax=occ_d+stdy
! transit time contained in interval
  tmin=occ_t0-(1.d0+xmax)/occ_v
  tmax=occ_t0+(1.d0-xmin)/occ_v

 END SUBROUTINE rectangle_time_select

! sample of npo random points inside the occulter stationary cross section
 SUBROUTINE sample_occult(xmin,xmax,ymin,ymax,npo,npr)

 DOUBLE PRECISION, INTENT(IN):: xmin,xmax,ymin,ymax
  INTEGER, INTENT(IN):: npo ! number of points to be provided
! inside occulting cross section, with coordinates in x_occ, y_occ
  INTEGER, INTENT(OUT):: npr ! number of points tested in rectangle
! with coordinates in vectors x_rect, y_rect
  INTEGER i,n
  DOUBLE PRECISION xspan,yspan
  DOUBLE PRECISION, EXTERNAL ::  urand
! sides of rectangle  
  xspan=xmax-xmin
  yspan=ymax-ymin
! counter of point in occulter cross section 
  n=0
  DO i=1,nsam_x
    x_rect(i)= xmin+ urand()*xspan
    y_rect(i)= ymin+ urand()*yspan 
    IF(belong_occult(x_rect(i),y_rect(i)))THEN
       n=n+1    
       x_occ(n)=x_rect(i)
       y_occ(n)=y_rect(i)
       IF(n.eq.npo)THEN
          npr=i
          EXIT
       ENDIF
    ENDIF   
  ENDDO
  IF(i.gt.nsam_x.or.n.lt.npo)THEN
     WRITE(*,*)' sampling failed, npo= ',npo,' n= ',n,' i= ',i
     STOP
  ENDIF

 END SUBROUTINE sample_occult

 DOUBLE PRECISION FUNCTION limb_dark(sinf)
  DOUBLE PRECISION, INTENT(IN):: sinf ! sine of angle from center 
  DOUBLE PRECISION cosf
  cosf=sqrt(1.d0-sinf**2)
  limb_dark=1.d0+limb1*(1.d0-cosf) +limb2*(1.d0-cosf)**2
 END FUNCTION limb_dark
 
 DOUBLE PRECISION FUNCTION lum_normalize(nn)
  INTEGER, INTENT(IN) :: nn ! number of points for integral
  INTEGER i
  DOUBLE PRECISION sinf, step, lum

  step=1.d0/float(nn)
  sinf=0.d0
  lum=0.d0
  DO i=1,nn
    sinf=step*i
    lum=2.d0*sinf*limb_dark(sinf)*step+lum
  ENDDO
  lum_normalize=lum
 END FUNCTION lum_normalize

 LOGICAL FUNCTION belong_occult(x,y)

  DOUBLE PRECISION, INTENT(IN) :: x,y ! cartesian coordinates of tested point 
  DOUBLE PRECISION val

  val=occ_a*x**2 + 2*occ_b*x*(y-occ_d) + occ_c*(y-occ_d)**2
  belong_occult=val.le.1.d0
 END FUNCTION belong_occult

 LOGICAL FUNCTION belong_star(x,y,t,val)

  DOUBLE PRECISION, INTENT(IN) :: x,y,t ! cartesian coordinates of point, time 
  DOUBLE PRECISION, INTENT(OUT)::  val ! in output distance from 
                                       !  center of star disk
  DOUBLE PRECISION xp
  
  xp=x+occ_v*(t-occ_t0)
  val=xp**2+y**2
! unit of length is star radius
  belong_star=val.le.1.d0
  val=sqrt(val)
 END FUNCTION belong_star



END MODULE occ_mod

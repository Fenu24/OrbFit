PROGRAM occult_sim

 USE occ_mod
 USE fund_const

 IMPLICIT NONE

! declarations

DOUBLE PRECISION :: xmin,xmax,ymin,ymax,area,tmin,tmax, eig(2)
DOUBLE PRECISION :: npor, nr, relerror, square, lum, nor, sinf

INTEGER npo, npr, ntim
INTEGER j, i,n  ! loop indexes, counters
INTEGER ilum    ! units

! input options: shape, velocity, time, time step

 CALL occinop


! initialize random numbers (not necessary if random.dat is available)
!  CALL intrnd()

! compute rectangle
 CALL rectangle_time_select(xmin,xmax,ymin,ymax,area,tmin,tmax,eig)
 WRITE(*,100) xmin,xmax,ymin,ymax
100 FORMAT('min X, max X ',F6.3,1X,F6.3,' min Y, Max Y ',F6.3,1X,F6.3)
 WRITE(*,101) area,tmin,tmax, 1.d0/sqrt(eig(1)),1.d0/sqrt(eig(2))
101 FORMAT('area ',F6.3,' min t, max t ',F10.3,1X,F10.3,' axes ',F6.3,1X,F6.3)

! select poins inside occulting cross section at t=t0
 CALL libini
 CALL sample_occult(xmin,xmax,ymin,ymax,npo,npr)
 WRITE(*,*)' generated sample of ',npo,' points, tested ',npr
 npor=float(npo)
 square=(xmax-xmin)*(ymax-ymin)
 relerror=(npor/npr)/(area/square)-1.d0
 WRITE(*,102) relerror
102 FORMAT(' relative area error ',1P,D12.4)

! compute total luminosity before occultation
 nor=lum_normalize(10000)

! start loop on times
 DO j=-10,ntim_x
    t_occ(j+11)=tmin+j*occ_dt
    IF(t_occ(j+11).gt.tmax+10*occ_dt) EXIT
    n=0 !counter of random points in front of star
    lum=0.d0 ! initialize integral of luminosity
    DO i=1,npo
       IF(belong_star(x_occ(i),y_occ(i),t_occ(j+11),sinf))THEN
          n=n+1
! contribution to integral
          lum=lum+limb_dark(sinf)
! to make image define flag...          
       ENDIF
    ENDDO
    nr=float(n) !number of points occulting

    rel_lum(j+11)= (nor-(lum/npor)*(area/pig))/nor
! below would be withoyt limb darkening
!    rel_lum(j+11)=1.d0-(nr/npor)*(area/pig)
 ENDDO
 ntim=j+10

! output table of data
 CALL filopn(ilum,'occultation.fla','unknown')
  DO j=1,ntim
     WRITE(ilum,*) t_occ(j),rel_lum(j)
  ENDDO
 CALL filclo(ilum,'  ')

CONTAINS

! options input 
 SUBROUTINE occinop
  USE fund_const
! USE occ_mod
! IMPLICIT NONE
  INTEGER iunout,iunit
  CHARACTER*60 comment
  CHARACTER*6 progna
  CHARACTER*80 run, file
  CHARACTER*3 suffix
  LOGICAL ireq,found
  INTEGER lerun,lepro
  DOUBLE PRECISION longax,shortax,anglong,ang
  DOUBLE PRECISION, DIMENSION (2,2) :: rot,diag,normal

! initialize option readout
  progna='occult'; run='occult'; suffix='opt'
  CALL libini 
  CALL namini 
! default options for propag                                   
!  CALL filopl(iunit,'propag.def') 
!  CALL rdnam(iunit) 
!  CALL filclo(iunit,' ') 
! read specific options; currently no default
  CALL rmsp(run,lerun) 
  file=run(1:lerun)//'.'//suffix 
  INQUIRE(FILE=file,EXIST=found) 
  IF(found) THEN 
     CALL filopn(iunit,file,'OLD') 
     CALL rdnam(iunit) 
     CALL filclo(iunit,' ') 
  ELSE 
     write(*,*)'**** file not found: ',file 
     write(*,*)'******* ',progna,' using defaults ****' 
  ENDIF

! possible options for progna
  CALL rmsp(progna,lepro)
  file=progna(1:lepro)//'.key' 
  CALL rdklst(file) 
! check for non-existing options 
  CALL chkkey 
! open file containing options used
  CALL filopn(iunout,'occultation.log','UNKNOWN') 

! numerical accuracy/noise
  ireq=.false.
  npo=50000
  comment='number of sample points inside the disk'
  CALL input_int_opt('occult','npo',npo,ireq,found,comment,iunout)
  IF(.not.found)STOP '****missing npo*********'
  ireq=.true.
! shape and position of the occulting disk
!  comment='coefficient of X^2'
!  CALL input_rea_opt('occult','occ_a',occ_a,ireq,found,comment,iunout)
!  IF(.not.found)STOP '****missing a*********'
!  comment='coefficient of XY'
!  CALL input_rea_opt('occult','occ_b',occ_b,ireq,found,comment,iunout)
!  IF(.not.found)STOP '****missing b*********'
!  comment='coefficient of Y^2'
!  CALL input_rea_opt('occult','occ_c',occ_c,ireq,found,comment,iunout)
!  IF(.not.found)STOP '****missing c*********'
  comment='long semiaxis of ellipse'
  CALL input_rea_opt('occult','longax',longax,ireq,found,comment,iunout)
  IF(.not.found)STOP '****missing longax*********'
  comment='short semiaxis of ellipse'
  CALL input_rea_opt('occult','shortax',shortax,ireq,found,comment,iunout)
  IF(.not.found)STOP '****missing shortax*********'
  comment='angle between X axis and long axis (degrees)'
  CALL input_rea_opt('occult','anglong',anglong,ireq,found,comment,iunout)
  IF(.not.found)STOP '****missing anglong*********'
  diag=0.d0
  diag(1,1)=1.d0/longax**2
  diag(2,2)=1.d0/shortax**2
  ang=anglong*radeg
  rot(1,1)=cos(ang);rot(2,2)=cos(ang)
  rot(1,2)=-sin(ang);rot(2,1)=sin(ang)
  normal=MATMUL(rot,MATMUL(diag,TRANSPOSE(rot)))
  occ_a=normal(1,1);occ_b=normal(1,2);occ_c=normal(2,2)
! position of the ellipse center w.r. to equator  
  comment='displacement in Y'
  CALL input_rea_opt('occult','occ_d',occ_d,ireq,found,comment,iunout)
  IF(.not.found)STOP '****missing d*********'
! times of the occultations
  comment='relative velocity in star radius/day'
  CALL input_rea_opt('occult','occ_v',occ_v,ireq,found,comment,iunout)
  IF(.not.found)STOP '****missing v*********'
  comment='reference time (occulter at center)'
  CALL input_rea_opt('occult','occ_t0',occ_t0,ireq,found,comment,iunout)
  IF(.not.found)STOP '****missing t0*********'
  comment='time step for luminosity measurement'
  CALL input_rea_opt('occult','occ_dt',occ_dt,ireq,found,comment,iunout)
  IF(.not.found)STOP '****missing dt*********'
! limb darkening parameters: limb_dark=1-limb1*(1-cosf) +limb2*(1-cosf)**2
  comment='limb darkening, linear coef.'
  CALL input_rea_opt('occult','limb1',limb1,ireq,found,comment,iunout)
  IF(.not.found)STOP '****missing limb1*********'
  comment='limb darkening, quadratic coef.'
  CALL input_rea_opt('occult','limb2',limb2,ireq,found,comment,iunout)
  IF(.not.found)STOP '****missing limb2*********'
  CALL filclo(iunout,' ')

 END SUBROUTINE occinop
END PROGRAM occult_sim

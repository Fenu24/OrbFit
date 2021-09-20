MODULE cat_debias
  USE healpix_types
  IMPLICIT NONE
  PRIVATE 
! common data: biases per catalog
  DOUBLE PRECISION, ALLOCATABLE:: catrecord(:,:)   !catrecord(49152,26)
  INTEGER(KIND=i4b), PARAMETER :: ns_max=8192 ! 2^13 : largest nside available

  !initialise array x2pix, y2pix and pix2x, pix2y used in several routines
  integer(KIND=i4b), save, dimension(128) :: x2pix=0,y2pix=0

  integer(KIND=i4b), save, dimension(0:1023) :: pix2x=0, pix2y=0

! public routines
  PUBLIC cat_bias
  
  CONTAINS

  SUBROUTINE cat_bias_init(error_model)
    INTEGER                     :: error, i,iun
    CHARACTER*(*),INTENT(IN)    :: error_model
! opening file and allocate catrecord, 3 catalogs times 4 field
    IF(error_model.eq.'cbm10')THEN
       CALL filopl(iun,'xcat.bias_cbm10')
       ALLOCATE(catrecord(49152,3*4))
! strip header, two lines
       READ(iun,*)
       READ(iun,*)
! 19 catalogs times 4 fields
    ELSEIF(error_model.EQ.'fcct14'.OR.error_model.EQ.'vfcc17')THEN
       CALL filopl(iun,'xcat.bias_fcct14')
       ALLOCATE(catrecord(49152,19*4))
       ! strip header, six lines
       DO i=1,6
          READ(iun,*)
       END DO
    ELSE
       WRITE(*,*) 'cat_bias_init : **** WRONG ERROR MODEL ***** ',error_model
       STOP
    END IF

! read bias table      
    DO i=1,49152
       READ(iun,*,IOSTAT=error) catrecord(i,:)
       IF (error/=0)THEN
         WRITE(*,*)' cat_bias_init: error in reading, record=',i 
         WRITE(*,*) catrecord(i,:)
       ENDIF
    ENDDO
   END SUBROUTINE cat_bias_init

  !------------------------------------------------------------------
  ! SUBROUTINE cat_bias
  ! The purpose of this routine is to provide the astrometric catalog
  ! biases from the xcat.bias file. Biases provided by S. Chesley and
  ! Baer.
  ! Outputs are bias in RA and DEC, given the coordinates of the 
  ! observations
  ! Creation on 21 May 2009 by F. Bernardi
  !-----------------------------------------------------------------

  SUBROUTINE cat_bias(phi,del,catcodmpc,bias_RA,bias_DEC,pmRA,pmDEC,biasflag,error_model)
    USE healpix_types
    DOUBLE PRECISION,      INTENT(OUT) :: bias_DEC,pmDEC
    DOUBLE PRECISION,      INTENT(OUT) :: bias_RA,pmRA
    LOGICAL,               INTENT(OUT) :: biasflag        
    CHARACTER*1,           INTENT(IN)  :: catcodmpc               ! this is MPC catalogue
!    DOUBLE PRECISION,      INTENT(IN)  :: catrecord(49152,26)   ! This array contains the cat biases
    DOUBLE PRECISION,      INTENT(IN)  :: del           ! DEC of the observations in radians
    DOUBLE PRECISION                   :: dpig=6.28318530717958648d0
    CHARACTER*(*),INTENT(IN)           :: error_model
    INTEGER(KIND=4)                    :: ipring
    INTEGER(KIND=4),       PARAMETER   :: nside=64      ! This is the SIDE number for the HEALpix pixellization
    DOUBLE PRECISION,      INTENT(IN)  :: phi           ! RA of the observations in radians
    DOUBLE PRECISION                   :: theta         ! colatitude in equatorial coo.
    LOGICAL,                      SAVE :: first=.true.
    INTEGER                            :: nskip         ! index for debias table
    IF(first)THEN
       CALL cat_bias_init(error_model)
       first=.false.
    ENDIF
! default if no success in finding debiasing
    bias_RA=0.d0
    bias_DEC=0.d0
    pmRA=0.d0
    pmDEC=0.d0
    biasflag=.false.
! find appropriate sky pixel
    theta=dpig/4-del         
    CALL ang2pix_ring(nside, theta, phi, ipring)
! use debiasing data, available for some catalogs
    IF(error_model.eq.'cbm10')THEN
       IF(catcodmpc.EQ.'a'.OR.catcodmpc.EQ.'b') THEN         ! Case for USNO-A1 catalog
          bias_RA=catrecord((ipring+1),1)
          bias_DEC=catrecord((ipring+1),2)
          pmRA=catrecord((ipring+1),3)
          pmDEC=catrecord((ipring+1),4)
          biasflag=.true.
       ELSEIF(catcodmpc.EQ.'c'.OR.catcodmpc.EQ.'d') THEN     ! Case for USNO-A2 catalog
          bias_RA=catrecord((ipring+1),5)
          bias_DEC=catrecord((ipring+1),6)
          pmRA=catrecord((ipring+1),7)
          pmDEC=catrecord((ipring+1),8)
          biasflag=.true.
       ELSEIF(catcodmpc.EQ.'o') THEN                         ! Case for USNO-B1 catalog
          bias_RA=catrecord((ipring+1),9)
          bias_DEC=catrecord((ipring+1),10)
          pmRA=catrecord((ipring+1),11)
          pmDEC=catrecord((ipring+1),12)
          biasflag=.true.
       ENDIF
    ELSEIF(error_model.EQ.'fcct14'.OR.error_model.EQ.'vfcc17')THEN
       SELECT CASE (catcodmpc)
       CASE('a')
! USNO-A1.0
          nskip=0
       CASE('b')
! USNO-SA1.0
          nskip=1
       CASE('c')
! USNO-A2.0
          nskip=2
       CASE('d')
! USNO-SA2.0
          nskip=3
       CASE('e')
! UCAC-1
          nskip=4
       CASE('g')
! Tycho-2 not to be debiased
          nskip=-1
       CASE('i')
! GSC-1.1
          nskip=6
       CASE('j')
! GSC-1.2
          nskip=7
       CASE('l')
! ACT not to be debiased
          nskip=-1
       CASE('m')
! GSC-ACT
          nskip=9
       CASE('o')
! USNO-B1.0
          nskip=10
       CASE('p')
! PPM
          nskip=11
       CASE('q')
! UCAC-4
          nskip=12
       CASE('r')
! UCAC-2
          nskip=13
       CASE('u')
! UCAC-3
          nskip=14
       CASE('v')
! NOMAD
          nskip=15
       CASE('w')
! CMC-14
          nskip=16
       CASE('L')
! 2MASS
          nskip=17
       CASE('N')
! SDSS-DRT7
          nskip=18
       CASE DEFAULT
          nskip=-1
       END SELECT
! Debiading of the selected catalog
       IF(nskip.GE.0)THEN
          bias_RA=catrecord((ipring+1),4*nskip+1)
          bias_DEC=catrecord((ipring+1),4*nskip+2)
          pmRA=catrecord((ipring+1),4*nskip+3)
          pmDEC=catrecord((ipring+1),4*nskip+4)
          biasflag=.true.
       END IF
    END IF
    
! add part for new debiasing scheme
  END SUBROUTINE cat_bias

  !=======================================================================
  subroutine ang2pix_ring(nside, theta, phi, ipix)
    !=======================================================================
    !     renders the pixel number ipix (RING scheme) for a pixel which contains
    !     a point on a sphere at coordinates theta and phi, given the map
    !     resolution parameter nside
    !=======================================================================
    INTEGER(KIND=i4b), INTENT(IN) :: nside
    INTEGER(KIND=i4b), INTENT(OUT) :: ipix
    REAL (KIND=dp), INTENT(IN) ::  theta, phi

    INTEGER(KIND=i4b) ::  nl4, jp, jm
    REAL (KIND=dp)  ::  z, za, tt, tp, tmp, temp1, temp2
    INTEGER(KIND=i4b) ::  ir, ip, kshift

    !-----------------------------------------------------------------------
    interface fatal_error
       module procedure fatal_error_womsg, fatal_error_msg
    end interface
    if (nside<1 .or. nside>ns_max) call fatal_error ("nside out of range")
    if (theta<0.0_dp .or. theta>pi)  then
       print*,"ANG2PIX_RING: theta : ",theta," is out of range [0, Pi]"
       call fatal_error
    endif

    z = COS(theta)
    za = ABS(z)
    tt = MODULO( phi, twopi) / halfpi  ! in [0,4)


    if ( za <= twothird ) then ! Equatorial region ------------------
       temp1 = nside*(.5_dp+tt)
       temp2 = nside*.75_dp*z
       jp = int(temp1-temp2) ! index of  ascending edge line
       jm = int(temp1+temp2) ! index of descending edge line

       ir = nside + 1 + jp - jm ! in {1,2n+1} (ring number counted from z=2/3)
       kshift = 1 - modulo(ir,2) ! kshift=1 if ir even, 0 otherwise

       nl4 = 4*nside
       ip = INT( ( jp+jm - nside + kshift + 1 ) / 2 ) ! in {0,4n-1}
       if (ip >= nl4) ip = ip - nl4

       ipix = 2*nside*(nside-1) + nl4*(ir-1) + ip

    else ! North & South polar caps -----------------------------

       tp = tt - INT(tt)      !MODULO(tt,1.0_dp)
       tmp = nside * SQRT( 3.0_dp*(1.0_dp - za) )

       jp = INT(tp          * tmp ) ! increasing edge line index
       jm = INT((1.0_dp - tp) * tmp ) ! decreasing edge line index

       ir = jp + jm + 1        ! ring number counted from the closest pole
       ip = INT( tt * ir )     ! in {0,4*ir-1}
       if (ip >= 4*ir) ip = ip - 4*ir

       if (z>0._dp) then
          ipix = 2*ir*(ir-1) + ip
       else
          ipix = 12*nside**2 - 2*ir*(ir+1) + ip
       endif

    endif

    return
  end subroutine ang2pix_ring

!-----------------------------------------------------
!  subroutine fatal_error (msg)
!    character(len=*), intent(in), optional :: msg
!
!    if (present(msg)) then
!       print *,'Fatal error: ', trim(msg)
!    else
!       print *,'Fatal error'
!    endif
!    call exit_with_status(1)
!  end subroutine fatal_error

  subroutine fatal_error_msg (msg)
    character(len=*), intent(in) :: msg
       print *,'Fatal error: ', trim(msg)
    call exit_with_status(1)
  end subroutine fatal_error_msg

  subroutine fatal_error_womsg
      print *,'Fatal error'
    call exit_with_status(1)
  end subroutine fatal_error_womsg

  subroutine exit_with_status (code, msg)
    ! ===========================================================
    integer, intent(in) :: code
    character (len=*), intent(in), optional :: msg
    ! ===========================================================
!  USE f90_unix, ONLY : iargc, getarg, exit
  
    if (present(msg)) print *,trim(msg)
    print *,'program exits with exit code ', code
    
!#if (defined (RS6000))
!    call exit_ (code)
!#else
!    call exit (code)
!#endif
  end subroutine exit_with_status


END MODULE cat_debias

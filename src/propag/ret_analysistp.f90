!===========================================================================!
! MODULE RET_ANALYSIS                                                       !
!===========================================================================!
! This module contains the routines for the resonant returns analysis       !
!===========================================================================!
! Authors: G. Tommei, A. Milani (Jan 2005)                                  !
! Last Update: July 2016 (A. Del Vigna, L. Dimare)                          !
!===========================================================================!
! Other modules used:                                                       !
!    fund_const                                                             !
!    tp_trace                                                               !
!    dyn_param                                                              !
!    multi_store                                                            !
! Subroutines contained:                                                    ! 
!    large_vis, risklarge, ret_min, ret_minctp, falslog4tp, achillestp,     !
!    findminctp, propclosapp, lovclosapptp, aftclo2v, wrireptp, header_rep, !
!    header_out, wriouttp, wriwarn, header_new, arrcut                      !
!===========================================================================!
MODULE ret_analysistp
  USE fund_const
  USE tp_trace
  USE dyn_param
  USE multi_store

  IMPLICIT NONE
  PRIVATE

  !**********************!
  !  PUBLIC SUBROUTINES  !
  !**********************!
  PUBLIC :: ret_minctp, ret_min, arrcut
  PUBLIC :: wrireptp, header_rep, wriwarn, header_new, header_out, wriouttp
  PUBLIC :: check_largevi, risklarge,  lovclosapptp!, propclosapp
  !********************!
  !  PUBLIC VARIABLES  !
  !********************!
  DOUBLE PRECISION, PUBLIC :: beta_factor    ! Beta factor (entangled case)
  LOGICAL,          PUBLIC :: segment_method ! Controls for methods   
  
CONTAINS 
  !=====================================================================!
  ! CHECK_LARGEVI (function)                                            !
  !=====================================================================!
  ! WARNING: this subroutine is now obsolete, and has been substituted  !
  ! by large_vis. This routine is currently used only by                !
  ! panst/resret_des.f90.                                               !
  !=====================================================================!
  LOGICAL FUNCTION check_largevi(vas,lre)
    USE output_control 
    !=======================================================================================================
    TYPE(tp_point), INTENT(IN) :: vas(lre) ! Return
    INTEGER,        INTENT(IN) :: lre      ! Length of the return
    !=======================================================================================================
    INTEGER, PARAMETER :: nvimin=10 ! Minimum number of VAs inside the Earth 
    INTEGER            :: j         ! Loop index
    INTEGER            :: nin       ! Number of points inside the Earth section
    !=======================================================================================================
    ! Count VAs inside the Earth
    !--------------------------------------------------------------------! 
    ! It is used the field d, which represents the actual distance (that !
    ! is, the distance on the MTP) expressed in Earth radii              !
    !--------------------------------------------------------------------!
    nin=0
    DO j=1,lre
       IF(vas(j)%d.LT.1.d0)THEN
          nin=nin+1
       ENDIF
    ENDDO
    IF(nin.GE.nvimin)THEN
       check_largevi=.TRUE.
    ELSE
       check_largevi=.FALSE.
    ENDIF
  END FUNCTION check_largevi

  !====================================================================!
  ! LARGE_VIS                                                          !
  !====================================================================!
  ! This subroutine takes as input a return vas, and count the number  !
  ! of segments of consecutive points inside the Earth section. For    !
  ! each one of these, the first and the last indexes are stored.      !
  !====================================================================!
  ! Authors: A. Del Vigna, L. Dimare                                   !
  ! Date of creation: July 2016                                        ! 
  !====================================================================!
  SUBROUTINE large_vis(vas,lre,num_large,ind1_large,ind2_large)
    USE output_control 
    !=======================================================================================================
    TYPE(tp_point), INTENT(IN)  :: vas(lre)       ! Return
    INTEGER,        INTENT(IN)  :: lre            ! Length of the return
    INTEGER,        INTENT(OUT) :: num_large      ! Number of large VIs
    INTEGER,        INTENT(OUT) :: ind1_large(10) ! Array with the first index of each large VI
    INTEGER,        INTENT(OUT) :: ind2_large(10) ! Array with the last index of eachlarge VI
    !                                             !----------------------------------------------!
    !                                             ! Maximum number of large VIs in a return = 10 !
    !                                             !----------------------------------------------!
    !=======================================================================================================
    INTEGER, PARAMETER :: nvimin=10 ! Minimum number of VAs inside the Earth 
    INTEGER            :: j         ! Loop index
    INTEGER            :: nin       ! Number of points inside the Earth section
    INTEGER            :: ind_start ! Temporary storage of the start of the large VI
    !=======================================================================================================
    nin=0
    num_large=0
    ind1_large=0
    ind2_large=0
    !******************************!
    !  Count VAs inside the Earth  !
    !******************************!
    !--------------------------------------------------------------------! 
    ! It is used the field d, which represents the actual distance (that !
    ! is, the distance on the MTP) expressed in Earth radii              !
    !--------------------------------------------------------------------!
    IF(vas(1)%d.LT.1.d0)THEN
       nin=1
       ind_start=1
    END IF
    DO j=2,lre
       IF(vas(j)%d.LT.1.d0 .AND. vas(j-1)%d.LT.1.d0)THEN
          ! The possible large VI is continuing
          nin=nin+1
       ELSEIF(vas(j)%d.GT.1.d0 .AND. vas(j-1)%d.LT.1.d0)THEN
          ! The possible large VI is ended
          !************************!
          !  Check if it is large  !
          !************************!
          IF(nin.GE.nvimin)THEN
             num_large=num_large+1
             ind1_large(num_large)=ind_start
             ind2_large(num_large)=j-1
          END IF
       ELSEIF(vas(j)%d.LT.1.d0 .AND. vas(j-1)%d.GT.1.d0)THEN
          ! A new possible large VI begins
          ind_start=j
          nin=1
       ENDIF
    ENDDO
    ! last point inside Earth
    IF(vas(lre)%d.LT.1.d0)THEN
       IF(nin.GE.nvimin)THEN
          num_large=num_large+1
          ind1_large(num_large)=ind_start
          ind2_large(num_large)=lre
       END IF
    ENDIF
  END SUBROUTINE large_vis

  !======================================================================!
  ! RISKLARGE                                                            !
  !======================================================================!
  ! This subroutine is called in case a large VI occurs. It computes     !
  ! the IP by means of a 1-dim formula, and store as VI                  !
  ! representative, the VA (inside the Earth) with the minimum           !
  ! distance from the Earth centre.                                      !
  !======================================================================!
  ! Authors: G. Tommei, A. Milani                                        !
  ! Last update: September 2016 (A. Del Vigna, L. Dimare)                !
  !======================================================================!
  SUBROUTINE risklarge(vas,lre,lmin,no_risk,t0,iunnew,iunwarn,iunrisk,riskfile,riskesafile)
    USE eval_risk
    USE multiple_sol,    ONLY: lovmagn, elm, unm, dpm, imip, imim, imul
    USE offlov_checktp , ONLY: header_risk, bsdmin, bsdmax
    USE planet_masses,   ONLY: gmearth
    USE virtual_impactor
    !=======================================================================================================
    TYPE(tp_point) ,  INTENT(IN)           :: vas(lre)    ! Return
    INTEGER,          INTENT(IN)           :: lre         ! Length of the return
    INTEGER,          INTENT(IN)           :: lmin        ! Location of the minimum b inside the return
    INTEGER,          INTENT(INOUT)        :: no_risk     ! Risk counter
    DOUBLE PRECISION, INTENT(IN)           :: t0          ! Elements initial time
    INTEGER,          INTENT(IN)           :: iunnew      ! Unit for the new file
    INTEGER,          INTENT(IN)           :: iunwarn     ! Unit for the warn file
    INTEGER,          INTENT(IN), OPTIONAL :: iunrisk     ! Unit for the risk file
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: riskfile    ! Risk file name
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: riskesafile ! ESA risk file name
    !===== Risk assessment =================================================================================
    DOUBLE PRECISION  :: tcl                               ! Time of closest approach
    DOUBLE PRECISION  :: deltat                            ! Interval time from t0 to tcl
    DOUBLE PRECISION  :: sigimp                            ! Sigma on the second LOV
    DOUBLE PRECISION  :: sigma                             ! Sigma LOV
    DOUBLE PRECISION  :: width                             ! Width
    DOUBLE PRECISION  :: p_imp                             ! Impact Probability
    DOUBLE PRECISION  :: dcur                              ! Distance of the minimum
    DOUBLE PRECISION  :: mu                                ! Gravitational constant
    DOUBLE PRECISION  :: b_e                               ! Earth impact cross section
    DOUBLE PRECISION  :: tpdist                            ! Distance from the center
    DOUBLE PRECISION  :: stretch                           ! Mean stretching in the large VI
    INTEGER           :: i1,i2                             ! First and last index of the large VI
    INTEGER           :: ii1,ii2                           ! First and last index for mean quantities computation
    DOUBLE PRECISION  :: deltasig                          ! Step size on the LOV
    INTEGER           :: ind                               ! Real index of the VI representative
    !===== Palermo Scale computation =======================================================================
    DOUBLE PRECISION  :: rrrr                              ! Real index of the minimum
    DOUBLE PRECISION  :: vvvv                              ! Velocity at infinity as output of lovmagn
    DOUBLE PRECISION  :: h                                 ! Magnitude as computed from lovmagn
    DOUBLE PRECISION  :: U                                 ! Velocity at infinity
    DOUBLE PRECISION  :: mass                              ! Mass
    DOUBLE PRECISION  :: v_imp                             ! Impact velocity
    DOUBLE PRECISION  :: energy                            ! Energy
    DOUBLE PRECISION  :: e_tilde                           ! Expected energy
    DOUBLE PRECISION  :: rel_prob                          ! Relative probability
    DOUBLE PRECISION  :: fb                                ! Frequency of impacts with given energy
    DOUBLE PRECISION  :: ps                                ! Palermo Scale
    INTEGER           :: ts                                ! Torino Scale
    !===== File units ======================================================================================
    INTEGER           :: iunrisk0                          ! Risk file
    INTEGER           :: iunrisk1                          ! ESA risk file
    INTEGER           :: iunnew0                           ! New file
    INTEGER           :: iunwarn0                          ! Warn file
    INTEGER           :: iunesarisk                        ! Unit for the ESA risk file
    !===== Other variables =================================================================================
    INTEGER           :: i,j                               ! Loop indexes
    INTEGER           :: le                                ! Length of riskfile
    INTEGER           :: lee                               ! Length of riskesafile
    CHARACTER(LEN=14) :: calend                            ! Calendar date
    INTEGER           :: year, month, dayint, hour, minute ! Date computation
    DOUBLE PRECISION  :: day                               ! Day of the date
    CHARACTER(LEN=2)  :: monthch,daych,hourch,minutech     ! Names for the date
    CHARACTER(LEN=17) :: impactdate                        ! Impact date
    DOUBLE PRECISION  :: diff_ts,cur_time
    INTEGER           :: curr_hour
    INTEGER           :: date_values(8)
    DOUBLE PRECISION, EXTERNAL :: tjm1
    !=======================================================================================================
    !********************************************************!
    !  Compute IP by finite approximation to 1-dim integral  !
    !********************************************************!
    p_imp=0
    i1=lre+1
    i2=0
    DO j=1,lre
       IF(vas(j)%d.LT.1.d0)THEN
          IF(prob_sampl)THEN
! With variable step size, the step is just the difference of the the sigma values
             deltasig=sigmavalv(NINT(vas(j)%rindex))-sigmavalv(NINT(vas(j)%rindex)-1)
          ELSE
! With fixed step in sigma_lov
             deltasig=delta_sigma
          ENDIF
          p_imp=p_imp+deltasig*EXP(-vas(j)%sigma**2/2.d0)/SQRT(dpig)
          i1=MIN(i1,j)
          i2=MAX(i2,j)
       ENDIF
    ENDDO
    
    IF(p_imp.GT.1.d0) p_imp=1.d0 ! Security against IP > 1

    !****************************************************************************!
    !  Choose as VI representative the one with minimum b from the Earth centre  !
    !****************************************************************************!
    ind=NINT(vas(lmin)%rindex)             ! It should be always an integer 
    sigma=vas(lmin)%sigma
    curr_vi%ele=elm(ind)                   ! Elements of the VI repr.  
    curr_vi%unc=unm(ind)                   ! Uncertainty matrices of the VI repr.
    curr_vi%tp=vas(lmin)                   ! TP point
    curr_vi%dyn=dyn                        ! Assignment to the dyn field, inizialization (25/07/2016)
    curr_vi%dyn%dp(1:ndyx)=dpm(1:ndyx,ind) ! The actual non-grav parameters are the ones of the ind-th VA
    stored_vi=.TRUE.
    !*****************************!
    !  Mean width and stretching  !
    !*****************************!
    ! Computed using points outside the impact cross section
    iunwarn0=ABS(iunwarn)
    IF(i1.GT.1)THEN
       ii1=i1-1
    ELSE
       ii1=i1
       IF(iunwarn.LT.0)WRITE(*,*)' VI ends at beginning of return'
       WRITE(iunwarn0,*)' VI ends at beginning of return'
       WRITE(*,*)' VI ends at beginning of return'
    ENDIF
    IF(i2.LT.lre)THEN
       ii2=i2+1
    ELSE
       ii2=i2
       IF(iunwarn.LT.0) WRITE(*,*)' VI ends at end of return'
       WRITE(iunwarn0,*) ' VI ends at end of return'
    ENDIF
    tpdist=SQRT((vas(ii1)%opik%coord(4)-vas(ii2)%opik%coord(4))**2+         &
         &      (vas(ii1)%opik%coord(5)-vas(ii2)%opik%coord(5))**2)
    stretch=tpdist/ABS(vas(ii1)%sigma-vas(ii2)%sigma)
    width=(vas(ii1)%width+vas(ii2)%width)*0.5d0
    !*************!
    !  Warn file  !
    !*************!
    IF(iunwarn.LT.0) WRITE(*,*) ' Computing IP for risk of type large'
    IF(iunwarn.LT.0) WRITE(*,*) ' VI-filament from ',FLOOR(vas(i1)%rindex),' to ', FLOOR(vas(i2)%rindex)
    IF(iunwarn.LT.0) WRITE(*,*) ' Impact Probability ',p_imp
    WRITE(iunwarn0,*) ' Computing IP for risk of type large' 
    WRITE(iunwarn0,*) ' VI-filament from ',FLOOR(vas(i1)%rindex),' to ', FLOOR(vas(i2)%rindex)
    WRITE(iunwarn0,*) ' Impact Probability ',p_imp

    CALL  wriwarn(vas(lmin)%tcla,vas(lmin)%rindex,vas(lmin)%b, &
         &           vas(lmin)%dd2_ds,vas(lmin)%stretch,vas(lmin)%width,vas(lmin)%moid,       &
         &           0.d0,0.d0,0,'LARGEVI','RLCP',iunnew) 

    !*************!
    !  Risk file  !
    !*************!
    IF(p_imp.GT.1e-11)THEN 
       ! Palermo scale computation
       U=vas(ii1)%opik%coord(1)  ! Velocity at infinity (using ii1, to avoid too much vicinity to the Earth centre)
       tcl=vas(lmin)%tcla        ! Time of close approach
       deltat=tcl-t0             ! Distance from t0
       rrrr=vas(lmin)%rindex     ! Real indexd of the minimum
       CALL lovmagn(rrrr,vvvv,h) ! Magnitude computation
       ps=palermo(U*reau,h,p_imp,deltat,mass,v_imp,energy,e_tilde,fb,rel_prob,iunwarn)     
       ! calendar date                                                         
       CALL calendwri(tcl,calend)
       ! Write calendar date in ESA format
       READ(calend(1:14),'(I4,1X,I2,1X,F6.3)') year, month, day
       dayint=FLOOR(day)
       hour=FLOOR((day-dayint)*24d0)
       minute=NINT(((day-dayint)*24d0-hour)*60d0)
       IF(month.LT.10)THEN 
          WRITE(monthch,'(A1,I1)') '0',month
       ELSE
          WRITE(monthch,'(I2)') month
       ENDIF
       IF(dayint.LT.10)THEN 
          WRITE(daych,'(A1,I1)') '0',dayint
       ELSE
          WRITE(daych,'(I2)') dayint
       ENDIF
       IF(hour.LT.10)THEN 
          WRITE(hourch,'(A1,I1)') '0',hour
       ELSE
          WRITE(hourch,'(I2)') hour
       ENDIF
       IF(minute.lt.10)THEN 
          WRITE(minutech,'(A1,I1)') '0',minute
       ELSE
          WRITE(minutech,'(I2)') minute
       ENDIF
       WRITE(impactdate,'(1X,I4,2(A1,A2),2(1X,A2))') year,'-',monthch,'-',daych,hourch,minutech
       ! Other stuff for the output
       sigimp=0.d0
       dcur=vas(lmin)%b
       ! Store information on size of Earth impact cross section 
       mu=gmearth/reau**3
       b_e=SQRT(1.d0+(2.d0*mu)/U**2)
       IF(b_e.GT.bsdmax) bsdmax=b_e 
       IF(b_e.LT.bsdmin) bsdmin=b_e 
       no_risk=no_risk+1 ! Risk counter
       IF(PRESENT(riskfile))THEN
          CALL rmsp(riskfile,le) 
          iunrisk0=44 
          iunrisk1=iunrisk0
          OPEN(UNIT=iunrisk0, FILE=riskfile(1:le),POSITION='APPEND')
       ELSEIF(PRESENT(iunrisk))THEN
          iunrisk0=abs(iunrisk)
          iunrisk1=iunrisk
       ELSE
          WRITE(*,*) 'risklarge: missing output optional arguments'
          STOP
       ENDIF
       !*****************!
       !  ESA risk file  !
       !*****************!
       IF(PRESENT(riskesafile))THEN
          CALL rmsp(riskesafile,lee) 
          iunesarisk=49 
          OPEN(UNIT=iunesarisk, FILE=riskesafile(1:lee),POSITION='APPEND')
          !------------!
          ! Compute TS !
          !------------!
          CALL DATE_AND_TIME(VALUES=date_values)
          curr_hour = date_values(5)+(date_values(6)+date_values(7)/60.d0)/60.d0
          cur_time  = tjm1(date_values(3),date_values(2),date_values(1),curr_hour)
          diff_ts = ABS(tcl-cur_time)/365.25d0
          IF(diff_ts.LT.100.d0)THEN
             ts=torino(p_imp,(e_tilde/p_imp))
          ELSE
             ts=-1
          END IF
          WRITE(iunesarisk,50) impactdate,p_imp,ps,ts, v_imp
50        FORMAT(A17,2X,1P,E9.2,4X,0P,F6.2,5X,I2,6X,F7.2)
       END IF
       IF(iunrisk1.LE.0) CALL header_risk(0)
       !****************!
       !  Store the VI  !
       !****************!
       WRITE(*,*) ' enter store_vi'
       CALL store_vi()
       !****************!
       !  Final output  !
       !****************!
       IF(width.GE.1.0d4)THEN 
          IF(iunrisk1.LT.0) WRITE(*,200)calend,tcl,sigma,sigimp,dcur,          &
               &              width,stretch,p_imp,e_tilde,ps
          WRITE(iunrisk0,200)calend,tcl,sigma,sigimp,dcur,          &
               &              width,stretch,p_imp,e_tilde,ps                     
200       FORMAT(A14,1X,F10.3,1X,F7.3,1X,F5.3,1X,F7.2,' +/- ',      &
               &              F8.2,1X,1P,E9.2,1X,E9.2,1X,E9.2,1X,0P,F6.2)         
       ELSE
          IF(iunrisk1.LT.0) WRITE(*,100)calend,tcl,sigma,sigimp,dcur,          &
               &              width,stretch,p_imp,e_tilde,ps
          WRITE(iunrisk0,100)calend,tcl,sigma,sigimp,dcur,          &
               &              width,stretch,p_imp,e_tilde,ps                     
100       FORMAT(A14,1X,F10.3,1X,F7.3,1X,F5.3,1X,F7.2,' +/- ',      &
               &              F8.3,1X,1P,E9.2,1X,E9.2,1X,E9.2,1X,0P,F6.2)         
       ENDIF
       IF(PRESENT(riskfile)) CLOSE(iunrisk0)
    ENDIF
  END SUBROUTINE risklarge
  
  !============================================================================
  ! RET_MIN (public subroutine) FITOBS
  !============================================================================
  SUBROUTINE ret_min(lre,vas_tr,t0,dnewton,no_risk,iunnew,iunwarn,iunrisk,limit_stretch)  
    USE offlov_checktp 
    USE output_control    
    USE virtual_impactor
    !=======================INPUT================================================
    INTEGER, INTENT(IN) :: lre  ! no close approach records in filament  
    TYPE(tp_point), DIMENSION(lre),INTENT(IN) :: vas_tr 
    DOUBLE PRECISION, INTENT(IN) :: t0           ! epoch 
    DOUBLE PRECISION, INTENT(IN)   :: dnewton    ! control for falsi/newton 
    INTEGER, INTENT(IN):: iunnew,iunwarn,iunrisk ! output units
    DOUBLE PRECISION, INTENT(IN) :: limit_stretch
    !=======================OUTPUT===============================================
    INTEGER, INTENT(OUT) :: no_risk   ! number of risk cases found
    !=================== END INTERFACE===========================================
    DOUBLE PRECISION, DIMENSION(lre) :: dminpos
    INTEGER :: i,j 
    DOUBLE PRECISION, DIMENSION(lre) :: dr2ds 
    DOUBLE PRECISION :: dminpo
    DOUBLE PRECISION :: x1,x2 
    CHARACTER(LEN=7) :: type          ! output from falslog               
    CHARACTER(LEN=100) :: riskfile
    TYPE(tp_point) :: va_tracemin      ! output from findminctp   
    DOUBLE PRECISION :: tmin,smin, distmin 
    LOGICAL :: fals_conv,fals_notp 
    INTEGER :: nvai 
    INTEGER :: niter !,iunwarn0,iunnew0                                                         
    DOUBLE PRECISION :: siglim,dalpha1,dalpha2,pridif,dfalsi, deltasig
    INTEGER :: i1,i2, lmin, int1,int2
    DOUBLE PRECISION, PARAMETER :: del_fal=2.3d-3    ! convergence control 
    ! this is the limit change in distance in RE (=1e-7 AU)  
    INTEGER :: num_large, ind1_large(10), ind2_large(10), ind1, ind2
    !=============================================================================
    dminpos=vas_tr%minposs
    no_risk=0    
    ! derivative of distance squared with respect to sigma                  
    dr2ds(1:lre)=vas_tr(1:lre)%dd2_ds 
    ! scan filament                                                         
    i1=vas_tr(1)%rindex 
    i2=vas_tr(lre)%rindex 
    !-----------------------------------------------!
    ! WARNING: the use of check_largevi is obsolete !
    !-----------------------------------------------!
!!$    IF(check_largevi(vas_tr(1:lre),lre))THEN
!!$       ! finding minimum distance and index lmin
!!$       distmin=MINVAL(vas_tr(1:lre)%b)
!!$       DO i=1,lre
!!$          IF(vas_tr(i)%b.eq.distmin)lmin=i
!!$       ENDDO
!!$       WRITE(*,*)' enter risklarge'
!!$       CALL risklarge(vas_tr(1:lre),lre,lmin,no_risk,t0,iunnew,iunwarn,iunrisk)
!!$       WRITE(*,*)' exit risklarge'
    CALL large_vis(vas_tr(1:lre),lre,num_large,ind1_large,ind2_large)
    IF(num_large.GT.0)THEN
       DO j=1,num_large
          ind1=MAX(ind1_large(j)-1,1)
          ind2=MIN(lre,ind2_large(j)+1)
          lmin=MINLOC(vas_tr(ind1:ind2)%b,DIM=1)
          CALL risklarge(vas_tr(ind1:ind2),ind2-ind1+1,lmin,no_risk,t0,&
               & iunnew,iunwarn,RISKFILE=riskfile)
       END DO
    ELSE
       IF(.NOT.prob_sampl) deltasig=delta_sigma
       DO 1 j=1,lre+1
          IF(j.eq.1)THEN 
             ! handle heads (including singletons)                                   
             IF(dr2ds(1).gt.0.d0)THEN 
                IF(lre.eq.1)THEN 
                   type='SINGLET' 
                ELSE 
                   type='ASCHEAD' 
                ENDIF
                ! increasing head                                                       
                IF(dminpos(1).lt.dnewton)THEN
                   IF(limit_stretch.GT.0.d0)THEN
                      IF(vas_tr(1)%stretch.GT.limit_stretch .OR. &
                           vas_tr(1)%stretch*vas_tr(1)%width.GT.vas_tr(1)%b**2/2.d0*1.d11)THEN
                      CYCLE
                      END IF
                   END IF
                   IF(prob_sampl)THEN
                      deltasig=sigmavalv(i1+1)-sigmavalv(i1)
                   ENDIF
                   siglim=deltasig/2.d0 
                   CALL achillestp(vas_tr(1),siglim,type,        &
                        &                iunwarn,iunnew,va_tracemin,niter,fals_conv,fals_notp,deltasig,limit_stretch)
                   IF(fals_notp)THEN 
                      WRITE(iunnew,*)' achillestp no TP ',type 
                   ELSE 
                      CALL riskchecktp(va_tracemin,t0,type,no_risk,iunnew,iunwarn,iunrisk)
                   ENDIF
                ENDIF
             ENDIF
          ELSEIF(j.eq.lre+1)THEN 
             ! handle tails (including singletons)                                   
             IF(dr2ds(lre).lt.0.d0)THEN 
                IF(lre.eq.1)THEN 
                   type='SINGLET' 
                ELSE 
                   type='DESTAIL' 
                ENDIF
                ! decreasing tail                                                       
                IF(dminpos(lre).lt.dnewton)THEN 
                   IF(limit_stretch.GT.0.d0)THEN
                      IF(vas_tr(lre)%stretch.GT.limit_stretch .OR. &
                           & vas_tr(lre)%stretch*vas_tr(lre)%width.GT.vas_tr(lre)%b**2/2.d0*1.d11)THEN
                         CYCLE
                      END IF
                   END IF
                   IF(prob_sampl)THEN
                      deltasig=sigmavalv(i2)-sigmavalv(i2-1)
                   ENDIF
                   siglim=deltasig/2.d0
                   CALL achillestp(vas_tr(lre),siglim,type,      &
                        &                iunwarn,iunnew,va_tracemin,niter,fals_conv,fals_notp,deltasig,limit_stretch)
                   ! if failure write failure record, else possibly risk record            
                   IF(fals_notp)THEN 
                      WRITE(iunnew,*)' achillestp no TP ',type 
                   ELSE 
                      CALL riskchecktp(va_tracemin,t0,type,no_risk,iunnew,iunwarn,iunrisk)
                   ENDIF
                ENDIF
             ENDIF
          ELSE 
             x2=vas_tr(j)%rindex 
             x1=vas_tr(j-1)%rindex 
             IF(prob_sampl)THEN
                CALL find_interval(x1,x2,int1,int2)
                deltasig=sigmavalv(int2)-sigmavalv(int1)
             ENDIF
             CALL falslog4tp(vas_tr(j-1),iunwarn,type) 
             ! check for interesting cases                                           
             dminpo=min(dminpos(j),dminpos(j-1)) 
             IF(dminpo.lt.dnewton)THEN 
                IF(limit_stretch.GT.0.d0)THEN
                   IF(vas_tr(j)%stretch.GT.limit_stretch .OR. &
                        & vas_tr(j)%stretch*vas_tr(j)%width.GT.vas_tr(j)%b**2/2.d0*1.d11)THEN
                   CYCLE
                   END IF
                END IF
!  if falslog proposes to try regula falsi, and the MOID is low, do it  
                IF(type.eq.'SIMPMIN'.or.type.eq.'INTEMIN')THEN 
                   IF(limit_stretch.GT.0.d0)THEN
                      IF(vas_tr(j-1)%stretch.GT.limit_stretch .OR. &
                           & vas_tr(j-1)%stretch*vas_tr(j-1)%width.GT.vas_tr(j-1)%b**2/2.d0*1.d11)THEN
                         CYCLE
                      END IF
                   END IF
                   CALL findminctp(vas_tr(j-1),x1,x2,type,                  &
                        &                iunwarn,iunnew,va_tracemin,                  &
                        &                fals_conv,niter,fals_notp,deltasig)  
                   ! if failure write failure record, else possibly risk record            
                   IF(fals_notp)THEN 
                      WRITE(iunnew,*)' falsi no TP ' 
                      siglim=deltasig/2.d0 
                      CALL achillestp(vas_tr(j-1),siglim,     &
                           &                      type,iunwarn,iunnew,va_tracemin,               &
                           &                      niter,fals_conv,fals_notp,deltasig,limit_stretch)                  
                      IF(fals_notp.or..not.fals_conv)THEN 
                         WRITE(iunnew,*)' achillestp no TP/conv ',type 
                      ELSE 
                         CALL riskchecktp(va_tracemin,t0,type,no_risk,iunnew,iunwarn,iunrisk)
                      ENDIF
                      WRITE(iunnew,*)' falsi no TP ' 
                      siglim=delta_sigma/2.d0
                      CALL achillestp(vas_tr(j),siglim,type,iunwarn,iunnew,va_tracemin,    &
                           &                      niter,fals_conv,fals_notp,deltasig,limit_stretch)                  
                      IF(fals_notp.or..not.fals_conv)THEN 
                         WRITE(iunnew,*)' achillestp no TP/conv ',type 
                      ELSE 
                         CALL riskchecktp(va_tracemin,t0,type,no_risk,iunnew,iunwarn,iunrisk) 
                      ENDIF
                   ELSE 
                      CALL riskchecktp(va_tracemin,t0,type,no_risk,iunnew,iunwarn,iunrisk)
                   ENDIF
                   ! but for the interrupted minima, there is always the suspicion that the
                   ! composite, that is two minima and one maximum                         
                   IF(type.eq.'INTEMIN')THEN 
                      ! change in direction of the LOV                                        
                      dalpha1=pridif(vas_tr(j-1)%alpha_lov,va_tracemin%alpha_lov) 
                      dalpha2=pridif(vas_tr(j)%alpha_lov,va_tracemin%alpha_lov)
                      ! store the minimum distance output from falsi                          
                      dfalsi=va_tracemin%b 
                      WRITE(iunwarn,*)'check for almost interrupted',     &
                           &                   dalpha1*degrad,dalpha2*degrad   
                      IF(cos(dalpha2).gt.0.70d0)THEN 
                         siglim=deltasig*abs(vas_tr(j-1)%rindex-va_tracemin%rindex)/2.d0 
                         CALL achillestp(vas_tr(j-1),siglim,     &
                              &                      type,iunwarn,iunnew,va_tracemin,               &
                              &                      niter,fals_conv,fals_notp,deltasig,limit_stretch)                  
                         IF(fals_notp.or..not.fals_conv)THEN 
                            WRITE(iunnew,*)' achillestp no TP/conv ',type 
                         ELSE 
                            IF(abs(dfalsi-va_tracemin%b).gt.del_fal*dfalsi)THEN
                               CALL riskchecktp(va_tracemin,t0,type,no_risk,iunnew,iunwarn,iunrisk) 
                            ENDIF
                         ENDIF
                      ELSEIF(cos(dalpha1).gt.0.7d0)THEN 
                         siglim=deltasig*abs(vas_tr(j)%rindex-va_tracemin%rindex)/2.d0 
                         CALL achillestp(vas_tr(j),siglim,        &
                              &                      type,iunwarn,iunnew,va_tracemin,               &
                              &                      niter,fals_conv,fals_notp,deltasig,limit_stretch)                  
                         IF(fals_notp.or..not.fals_conv)THEN 
                            WRITE(iunnew,*)' achillestp no TP/conv ', type 
                         ELSE 
                            IF(abs(dfalsi-va_tracemin%b).gt.del_fal*dfalsi)   THEN
                               CALL riskchecktp(va_tracemin,t0,type,no_risk,iunnew,iunwarn,iunrisk) 
                            ENDIF
                         ENDIF
                      ENDIF
                   ENDIF
                ELSEIF(type.eq.'INTEUNK')THEN 
                   siglim=deltasig/2.d0 
                   IF(dr2ds(j-1).lt.0.d0)THEN 
                      IF(limit_stretch.GT.0.d0)THEN
                         IF(vas_tr(j-1)%stretch.GT.limit_stretch .OR. &
                              & vas_tr(j-1)%stretch*vas_tr(j-1)%width.GT.vas_tr(j-1)%b**2/2.d0*1.d11)THEN
                            CYCLE
                         END IF
                      END IF
                      CALL achillestp(vas_tr(j-1),siglim,type,   &
                           &                iunwarn,iunnew,va_tracemin,niter,fals_conv,fals_notp,deltasig,limit_stretch)
                   ELSEIF(dr2ds(j).gt.0.d0)THEN 
                      CALL achillestp(vas_tr(j),siglim,type,     &
                           &                iunwarn,iunnew,va_tracemin,niter,fals_conv,fals_notp,deltasig,limit_stretch)
                   ELSE 
                      WRITE(iunnew,*)' case not understood ', type 
                   ENDIF
                   IF(fals_notp)THEN 
                      WRITE(iunnew,*)' achillestp no TP ',type 
                   ELSE 
                      CALL riskchecktp(va_tracemin,t0,type,no_risk,iunnew,iunwarn,iunrisk)
                   ENDIF
                ELSEIF(type.eq.'ENTADES')THEN 
                   IF(limit_stretch.GT.0.d0)THEN
                      IF(vas_tr(j-1)%stretch.GT.limit_stretch .OR. &
                           & vas_tr(j-1)%stretch*vas_tr(j-1)%width.GT.vas_tr(j-1)%b**2/2.d0*1.d11)THEN
                         CYCLE
                      END IF
                   END IF
                   ! look for first minimum
                   siglim=deltasig/2.d0 
                   CALL achillestp(vas_tr(j-1),siglim,type,      &
                        &                iunwarn,iunnew,va_tracemin,niter,fals_conv,fals_notp,deltasig,limit_stretch)
                   IF(fals_notp)THEN 
                      WRITE(iunnew,*)' achillestp no TP ',type 
                   ELSE 
                      CALL riskchecktp(va_tracemin,t0,type,no_risk,iunnew,iunwarn,iunrisk)
                   ENDIF
                   ! missing second half of the entagled case                              
                   WRITE(iunnew,*)type, ' second minimum not handled' 
                ELSEIF(type.eq.'ENTAASC')THEN 
                   ! look for last minimum    
                   siglim=deltasig/2.d0 
                   CALL achillestp(vas_tr(j),siglim,type,        &
                        &                iunwarn,iunnew,va_tracemin,niter,fals_conv,fals_notp,deltasig,limit_stretch)
                   IF(fals_notp)THEN 
                      WRITE(iunnew,*)' achillestp no TP ',type 
                   ELSE 
                      CALL riskchecktp(va_tracemin,t0,type,no_risk,iunnew,iunwarn,iunrisk)    
                   ENDIF
                   ! missing second half of the entagled case                              
                   WRITE(iunnew,*)type, ' second minimum not handled' 
                ELSE 
                   IF(type.eq.'NOWORRY'.or.type.eq.'INTEMAX'.or.          &
                        &                   type.eq.'SIMPMAX')THEN                         
                      !  ok not to be handled                                                 
                   ELSE 
                      WRITE(iunnew,*)' case not handled, type=',type,     &
                           &                   vas_tr(j-1)%tcla, vas_tr(j-1)%rindex    
                      WRITE(iunwarn,*)' case not handled, type=',type 
                      CALL wriouttp(vas_tr(j-1),iunwarn) 
                      CALL wriouttp(vas_tr(j),iunwarn) 
                   ENDIF
                ENDIF
             ENDIF
          ENDIF

1      ENDDO
    ENDIF

  END SUBROUTINE ret_min
  !============================================================================
  ! RET_MINCTP (public subroutine) RESRET2TP
  !============================================================================
  SUBROUTINE ret_minctp(lre,vas_tr,t0,iunwarn,iunnew,riskfile,dnewton,no_risk,limit_stretch,riskesafile)
    USE offlov_checktp  
    USE virtual_impactor
    !=======================INPUT================================================
    INTEGER, INTENT(IN) :: lre                   ! number of close approach
    ! records in filament  
    TYPE(tp_point), DIMENSION(lre), INTENT(IN) :: vas_tr 
    DOUBLE PRECISION, INTENT(IN) :: t0           ! epoch 
    DOUBLE PRECISION, INTENT(IN)   :: dnewton     ! control for falsi/newton 
    !=====================OUTPUT UNITS===========================================
    INTEGER            :: iunwarn,iunnew 
    CHARACTER(LEN=100) :: riskfile  
    DOUBLE PRECISION, INTENT(IN) :: limit_stretch
    CHARACTER(LEN=100), INTENT(IN), OPTIONAL :: riskesafile  
    !=======================OUTPUT===============================================
    INTEGER, INTENT(OUT) :: no_risk 
    !=================== END INTERFACE===========================================
    DOUBLE PRECISION, DIMENSION(lre) :: dminpos
    INTEGER :: i,j 
    DOUBLE PRECISION, DIMENSION(lre) :: dr2ds 
    DOUBLE PRECISION :: dminpo
    DOUBLE PRECISION :: x1,x2 
    CHARACTER(LEN=7) :: typ          ! output from falslog               
    TYPE(tp_point) :: va_tracemin      ! output from findminctp   
    DOUBLE PRECISION :: tmin,smin,distmin 
    LOGICAL :: fals_conv,fals_notp 
    INTEGER :: nbigrisk, nsegrisk
    INTEGER :: nvai 
    INTEGER :: niter                                                         
    DOUBLE PRECISION :: siglim,dalpha1,dalpha2,pridif,dfalsi, deltasig
    INTEGER :: i1,i2, lmin, int1, int2
    DOUBLE PRECISION, PARAMETER :: del_fal=2.3d-3    ! convergence control 
    ! this is the limit change in distance in RE (=1e-7 AU)  
    !============================================================================
    INTEGER :: num_large, ind1_large(10), ind2_large(10), ind1, ind2
    !============================================================================
    dminpos=vas_tr%minposs
    no_risk=0
! derivative of distance squared with respect to sigma
    dr2ds(1:lre)=vas_tr(1:lre)%dd2_ds 
! scan filament                                                         
    i1=vas_tr(1)%rindex 
    i2=vas_tr(lre)%rindex
    !-----------------------------------------------!
    ! WARNING: the use of check_largevi is obsolete !
    !-----------------------------------------------!
!!$    IF(check_largevi(vas_tr(1:lre),lre))THEN
!!$       ! finding minimum distance and index lmin
!!$       distmin=MINVAL(vas_tr(1:lre)%b)
!!$       DO i=1,lre
!!$          IF(vas_tr(i)%b.eq.distmin)lmin=i
!!$       ENDDO
!!$       CALL risklarge(vas_tr(1:lre),lre,lmin,no_risk,t0,iunnew,iunwarn,RISKFILE=riskfile,RISKESAFILE=riskesafile)
    CALL large_vis(vas_tr(1:lre),lre,num_large,ind1_large,ind2_large)
    IF(num_large.GT.0)THEN
       DO j=1,num_large
          ind1=MAX(ind1_large(j)-1,1)
          ind2=MIN(lre,ind2_large(j)+1)
          lmin=MINLOC(vas_tr(ind1:ind2)%b,DIM=1)
          CALL risklarge(vas_tr(ind1:ind2),ind2-ind1+1,lmin,no_risk,t0,&
                & iunnew,iunwarn,RISKFILE=riskfile,RISKESAFILE=riskesafile)
       END DO
    ELSE
       IF(.NOT.prob_sampl) deltasig=delta_sigma
       DO 1 j=1,lre+1
          IF(j.EQ.1)THEN 
! handle heads (including singletons)                                   
             IF(dr2ds(1).GT.0.d0)THEN 
                IF(lre.EQ.1)THEN 
                   typ='SINGLET' 
                ELSE 
                   typ='ASCHEAD'
                ENDIF
! increasing head                                                       
                IF(dminpos(1).LT.dnewton)THEN 
                   IF(limit_stretch.GT.0.d0)THEN
                      IF(vas_tr(1)%stretch.GT.limit_stretch .OR. &
                           vas_tr(1)%stretch*vas_tr(1)%width.GT.vas_tr(1)%b**2/2.d0*1.d11)THEN
                      CYCLE
                      END IF
                   END IF
                   IF(prob_sampl)THEN
                      deltasig=sigmavalv(i1+1)-sigmavalv(i1)
                   ENDIF
                   siglim=deltasig/2.d0 
                   CALL achillestp(vas_tr(1),siglim,typ,        &
                        &                iunwarn,iunnew,va_tracemin,niter,fals_conv,fals_notp,deltasig,limit_stretch)
                   IF(fals_notp)THEN 
                      !                      WRITE(iunnew,*)' achillestp no TP ',typ 
                   ELSE 
                      CALL riskchecktp(va_tracemin,t0,typ,no_risk,    &
                           &                   iunnew,iunwarn,RISKFILE=riskfile,RISKESAFILE=riskesafile)
                   ENDIF
                ENDIF
             ENDIF
          ELSEIF(j.EQ.lre+1)THEN 
! handle tails (including singletons)                                   
             IF(dr2ds(lre).LT.0.d0)THEN 
                IF(lre.EQ.1)THEN 
                   typ='SINGLET' 
                ELSE 
                   typ='DESTAIL' 
                ENDIF
! decreasing tail                                                       
                IF(dminpos(lre).LT.dnewton)THEN 
                   IF(limit_stretch.GT.0.d0)THEN
                      IF(vas_tr(lre)%stretch.GT.limit_stretch .OR. &
                           & vas_tr(lre)%stretch*vas_tr(lre)%width.GT.vas_tr(lre)%b**2/2.d0*1.d11)THEN
                         CYCLE
                      END IF
                   END IF
                   IF(prob_sampl)THEN
                      deltasig=sigmavalv(i2)-sigmavalv(i2-1)
                   ENDIF
                   siglim=deltasig/2.d0 
                   CALL achillestp(vas_tr(lre),siglim,typ,      &
                        &                iunwarn,iunnew,va_tracemin,niter,fals_conv,fals_notp,deltasig,limit_stretch)
                   ! if failure write failure record, else possibly risk record            
                   IF(fals_notp)THEN 
                      !                      WRITE(iunnew,*)' achillestp no TP ',typ 
                   ELSE 
                      CALL riskchecktp(va_tracemin,t0,typ,no_risk,    &
                           &                   iunnew,iunwarn,RISKFILE=riskfile,RISKESAFILE=riskesafile)
                   ENDIF
                ENDIF
             ENDIF
          ELSE 
             x2=vas_tr(j)%rindex 
             x1=vas_tr(j-1)%rindex 
             IF(prob_sampl)THEN
                CALL find_interval(x1,x2,int1,int2)
                deltasig=sigmavalv(int2)-sigmavalv(int1)
             ENDIF
             CALL falslog4tp(vas_tr(j-1),iunwarn,typ) 
! check for interesting cases                                           
             dminpo=MIN(dminpos(j),dminpos(j-1)) 
             IF(dminpo.LT.dnewton)THEN 
                IF(limit_stretch.GT.0.d0)THEN
                   IF(vas_tr(j)%stretch.GT.limit_stretch .OR. &
                      & vas_tr(j)%stretch*vas_tr(j)%width.GT.vas_tr(j)%b**2/2.d0*1.d11)THEN
                      CYCLE
                   END IF
                END IF
!  if falslog proposes to try regula falsi, and the MOID is low, do it  
                IF(typ.EQ.'SIMPMIN'.OR.typ.EQ.'INTEMIN')THEN 
                   IF(limit_stretch.GT.0.d0)THEN
                      IF(vas_tr(j-1)%stretch.GT.limit_stretch .OR. &
                           & vas_tr(j-1)%stretch*vas_tr(j-1)%width.GT.vas_tr(j-1)%b**2/2.d0*1.d11)THEN
                         CYCLE
                      END IF
                   END IF
                   CALL findminctp(vas_tr(j-1),x1,x2,typ,                  &
                        &                iunwarn,iunnew,va_tracemin,                  &
                        &                fals_conv,niter,fals_notp,deltasig)                        
                   ! if failure write failure record, else possibly risk record            
                   IF(fals_notp)THEN 
                      WRITE(iunnew,*)' falsi no TP '
                      siglim=deltasig/2.d0 
                      CALL achillestp(vas_tr(j-1),siglim,     &
                           &                      typ,iunwarn,iunnew,va_tracemin,               &
                           &                      niter,fals_conv,fals_notp,deltasig,limit_stretch)                  
                      IF(fals_notp.or..not.fals_conv)THEN 
                         !                         WRITE(iunnew,*)' achillestp no TP/conv ',typ 
                      ELSE 
                         CALL riskchecktp(va_tracemin,t0,typ,      &
                              &                         no_risk,iunnew,iunwarn,RISKFILE=riskfile,RISKESAFILE=riskesafile)
                      ENDIF
                      siglim=deltasig/2.d0 
                      CALL achillestp(vas_tr(j),siglim,       &
                           &                      typ,iunwarn,iunnew,va_tracemin,               &
                           &                      niter,fals_conv,fals_notp,deltasig,limit_stretch)                  
                      IF(fals_notp.or..not.fals_conv)THEN 
                         !                         WRITE(iunnew,*)' achillestp no TP/conv ',typ 
                      ELSE 
                         CALL riskchecktp(va_tracemin,t0,typ,      &
                              &                         no_risk,iunnew,iunwarn,RISKFILE=riskfile,RISKESAFILE=riskesafile)
                      ENDIF
                   ELSE 
                      nbigrisk=0
                      CALL riskchecktp(va_tracemin,t0,typ,no_risk,    &
                           &                      iunnew,iunwarn,RISKFILE=riskfile,RISKESAFILE=riskesafile)
                   ENDIF
                   ! but for the interrupted minima, there is always the suspicion that the
                   ! composite, that is two minima and one maximum                         
                   IF(typ.eq.'INTEMIN')THEN 
                      ! change in direction of the LOV                                        
                      dalpha1=pridif(vas_tr(j-1)%alpha_lov,va_tracemin%alpha_lov) 
                      dalpha2=pridif(vas_tr(j)%alpha_lov,va_tracemin%alpha_lov)
                      ! store the minimum distance output from falsi                          
                      dfalsi=va_tracemin%b 
                      WRITE(iunwarn,*)'check for almost interrupted',     &
                           &                   dalpha1*degrad,dalpha2*degrad 
                      IF(cos(dalpha2).gt.0.70d0)THEN
                         siglim=deltasig*abs(vas_tr(j-1)%rindex-va_tracemin%rindex)/2.d0 
                         CALL achillestp(vas_tr(j-1),siglim,     &
                              &                      typ,iunwarn,iunnew,va_tracemin,               &
                              &                      niter,fals_conv,fals_notp,deltasig,limit_stretch)                  
                         IF(fals_notp.or..not.fals_conv)THEN 
                            !                            WRITE(iunnew,*)' achillestp no TP/conv ',typ 
                         ELSE 
                            IF(abs(dfalsi-va_tracemin%b).gt.del_fal*dfalsi)   &
                                 &                         CALL riskchecktp(va_tracemin,t0,typ, &
                                 &                         no_risk,iunnew,iunwarn,RISKFILE=riskfile,RISKESAFILE=riskesafile)
                         ENDIF
                      ELSEIF(cos(dalpha1).gt.0.7d0)THEN 
                         siglim=deltasig*abs(vas_tr(j)%rindex-va_tracemin%rindex)/2.d0 
                         CALL achillestp(vas_tr(j),siglim,        &
                              &                      typ,iunwarn,iunnew,va_tracemin,               &
                              &                      niter,fals_conv,fals_notp,deltasig,limit_stretch)                  
                         IF(fals_notp.or..not.fals_conv)THEN 
                            !                            WRITE(iunnew,*)' achillestp no TP/conv ', typ 
                         ELSE 
                            IF(abs(dfalsi-va_tracemin%b).gt.del_fal*dfalsi)   &
                                 &                         CALL riskchecktp(va_tracemin,t0,typ, &
                                 &                         no_risk,iunnew,iunwarn,RISKFILE=riskfile,RISKESAFILE=riskesafile)
                         ENDIF
                      ENDIF
                   ENDIF
                ELSEIF(typ.eq.'INTEUNK')THEN 
                   siglim=deltasig/2.d0 
                   IF(dr2ds(j-1).lt.0.d0)THEN 
                      IF(limit_stretch.GT.0.d0)THEN
                         IF(vas_tr(j-1)%stretch.GT.limit_stretch .OR. &
                              & vas_tr(j-1)%stretch*vas_tr(j-1)%width.GT.vas_tr(j-1)%b**2/2.d0*1.d11)THEN
                            CYCLE
                         END IF
                      END IF
                      CALL achillestp(vas_tr(j-1),siglim,typ,   &
                           &                iunwarn,iunnew,va_tracemin,niter,fals_conv,fals_notp,deltasig,limit_stretch)
                   ELSEIF(dr2ds(j).gt.0.d0)THEN
                      CALL achillestp(vas_tr(j),siglim,typ,     &
                           &                iunwarn,iunnew,va_tracemin,niter,fals_conv,fals_notp,deltasig,limit_stretch)
                   ELSE 
                      WRITE(iunnew,*)' case not understood ', typ 
                   ENDIF
                   IF(fals_notp)THEN 
                      !                      WRITE(iunnew,*)' achillestp no TP ',typ 
                   ELSE 
                      CALL riskchecktp(va_tracemin,t0,typ,            &
                           &                   no_risk,iunnew,iunwarn,RISKFILE=riskfile,RISKESAFILE=riskesafile)
                   ENDIF
                ELSEIF(typ.eq.'ENTADES')THEN 
                   IF(limit_stretch.GT.0.d0)THEN
                      IF(vas_tr(j-1)%stretch.GT.limit_stretch .OR. &
                           & vas_tr(j-1)%stretch*vas_tr(j-1)%width.GT.vas_tr(j-1)%b**2/2.d0*1.d11)THEN
                         CYCLE
                      END IF
                   END IF
                   ! look for first minimum
                   siglim=deltasig/2.d0 
                   CALL achillestp(vas_tr(j-1),siglim,typ,      &
                        &                iunwarn,iunnew,va_tracemin,niter,fals_conv,fals_notp,deltasig,limit_stretch)
                   IF(fals_notp)THEN 
                      !                      WRITE(iunnew,*)' achillestp no TP ',typ 
                   ELSE 
                      CALL riskchecktp(va_tracemin,t0,typ,            &
                           &                   no_risk,iunnew,iunwarn,RISKFILE=riskfile,RISKESAFILE=riskesafile)
                   ENDIF
                   ! missing second half of the entagled case                              
                   WRITE(iunnew,*)typ, ' second minimum not handled' 
                ELSEIF(typ.eq.'ENTAASC')THEN 
                   ! look for last minimum
                   siglim=deltasig/2.d0 
                   CALL achillestp(vas_tr(j),siglim,typ,        &
                        &                iunwarn,iunnew,va_tracemin,niter,fals_conv,fals_notp,deltasig,limit_stretch)
                   IF(fals_notp)THEN 
                      !                      WRITE(iunnew,*)' achillestp no TP ',typ 
                   ELSE 
                      CALL riskchecktp(va_tracemin,t0,typ,            &
                           &                   no_risk,iunnew,iunwarn,RISKFILE=riskfile,RISKESAFILE=riskesafile)
                   ENDIF
                   ! missing second half of the entagled case                              
                   WRITE(iunnew,*)typ, ' second minimum not handled' 
                ELSE 
                   IF(typ.eq.'NOWORRY'.or.typ.eq.'INTEMAX'.or.          &
                        &                   typ.eq.'SIMPMAX')THEN                         
                      !  ok not to be handled                                                 
                   ELSE 
                      WRITE(iunnew,*)' case not handled, typ=',typ,     &
                           &                   vas_tr(j-1)%tcla, vas_tr(j-1)%rindex    
                      WRITE(iunwarn,*)' case not handled, typ=',typ 
                      CALL wriouttp(vas_tr(j-1),iunwarn) 
                      CALL wriouttp(vas_tr(j),iunwarn) 
                   ENDIF
                ENDIF
             ENDIF
          ENDIF
1      ENDDO
    ENDIF

  END SUBROUTINE ret_minctp
  !============================================================================
  SUBROUTINE falslog4tp(vas_tracec,iunwar,type) 
    !===================== INPUT=========================================
    TYPE(tp_point), DIMENSION(2), INTENT(IN) :: vas_tracec    
    ! 2 close approach records
    !====================OUTPUT UNIT=====================================
    INTEGER :: iunwar 
    !===================== OUTPUT========================================
    !  ENCOUNTER TYPE                                                        
    !  SIMPMIN= simple minimum (not interrupted, presumably single minimum) 
    !  SIMPMAX= simple maximum (not interrupted, presumably single maximum) 
    !  NOWORRY= increasing or decreasing (not interrupted) NO               
    !  ENTAASC= entangled, multiple stationary points, first stationary poin
    !  ENTADES= entangled, multiple stationary points, first stationary poin
    !  INTEMIN= simple interrupted; at least one minimum in between, possibl
    !  INTEMAX= simple maximum (interrupted, presumably single maximum) NO  
    !  INTEUNK= interrupted, but not expected  DENSIFY                      
    !  UNKNOWN= LOV turning (change of direction between 60 and 90 deg) DENS
    CHARACTER(LEN=7) :: type 
    !================== END INTERFACE=======================================
    DOUBLE PRECISION :: dalpha,dbeta1,dbeta2,num,den,pridif 
    DOUBLE PRECISION :: betalim1,betalim2,gk,factor,sin20
    !======================================================================
    ! change in direction of the LOV 
    IF(((vas_tracec(2)%opik%coord(4)-vas_tracec(1)%opik%coord(4)))==0.d0.and.  &
         &  (vas_tracec(2)%opik%coord(5)-vas_tracec(1)%opik%coord(5))== 0.d0 )THEN
       type='NOWORRY'
    ELSE                                
       dalpha=pridif(vas_tracec(1)%alpha_lov,vas_tracec(2)%alpha_lov) 
       ! angle between tangent to LOV and vector between two TP points         
       num=(vas_tracec(2)%opik%coord(4)-vas_tracec(1)%opik%coord(4))*         &
            &    (-vas_tracec(1)%stretch_lov*sin(vas_tracec(1)%alpha_lov))+        &
            &     (vas_tracec(2)%opik%coord(5)-vas_tracec(1)%opik%coord(5))         &
            &    *(vas_tracec(1)%stretch_lov*cos(vas_tracec(1)%alpha_lov))         
       den=vas_tracec(1)%stretch_lov*           &
            &     sqrt((vas_tracec(2)%opik%coord(4)-vas_tracec(1)%opik%coord(4))**2+                   &
            &              (vas_tracec(2)%opik%coord(5)-vas_tracec(1)%opik%coord(5))**2)  
       dbeta1=acos(num/den) 
       num=(vas_tracec(2)%opik%coord(4)-vas_tracec(1)%opik%coord(4))      &
            &    *(-vas_tracec(2)%stretch_lov*sin(vas_tracec(2)%alpha_lov))+        &
            &     (vas_tracec(2)%opik%coord(5)-vas_tracec(1)%opik%coord(5))       &
            &    *(vas_tracec(2)%stretch_lov*cos(vas_tracec(2)%alpha_lov))         
       den=vas_tracec(2)%stretch_lov*           &
            &     sqrt((vas_tracec(2)%opik%coord(4)-vas_tracec(1)%opik%coord(4))**2+                   &
            &     (vas_tracec(2)%opik%coord(5)-vas_tracec(1)%opik%coord(5))**2)       
       dbeta2=acos(num/den) 
       ! selection of type                                                     
       factor=beta_factor
       gk=0.01720209895d0 
       sin20=0.34202014d0
       betalim1=MAX(factor*gk*vas_tracec(1)%b/vas_tracec(1)%opik%coord(1),sin20) 
       betalim2=MAX(factor*gk*vas_tracec(2)%b/vas_tracec(2)%opik%coord(1),sin20)
       IF(cos(dalpha).gt.0.d0)THEN 
          IF(abs(sin(dbeta1)).lt.betalim1.or.abs(sin(dbeta2)).lt.betalim2)THEN
             IF(vas_tracec(1)%dd2_ds.lt.0.d0.and.vas_tracec(2)%dd2_ds.gt.0.0d0)THEN 
                ! simple minimum                                                        
                type='SIMPMIN' 
             ELSEIF(vas_tracec(1)%dd2_ds.gt.0.d0.and.vas_tracec(2)%dd2_ds.lt.0.0d0)THEN 
                type='SIMPMAX' 
             ELSE 
                ! could be increasing or decreasing, we do not care                     
                type='NOWORRY' 
             ENDIF
          ELSE 
             ! this is a complicated case, to be further discriminated               
             IF(vas_tracec(2)%dd2_ds.gt.0.d0)THEN 
                type='ENTAASC' 
                IF(vas_tracec(1)%dd2_ds.lt.0.d0)THEN 
                   WRITE(iunwar,*)' ENTAASC but first decreasing' 
                ENDIF
             ELSE 
                type='ENTADES' 
                IF(vas_tracec(1)%dd2_ds.gt.0.d0)THEN 
                   WRITE(iunwar,*)' ENTADES but second increasing' 
                ENDIF
             ENDIF
          ENDIF
       ELSEIF(cos(dalpha).lt.0.d0)THEN 
          IF(vas_tracec(1)%dd2_ds.lt.0.d0.and.vas_tracec(2)%dd2_ds.gt.0.0d0)THEN 
             ! there is at least one minimum of the interrupted type in between      
             type='INTEMIN' 
          ELSEIF(vas_tracec(1)%dd2_ds.gt.0.d0.and.vas_tracec(2)%dd2_ds.lt.0.0d0)THEN 
             ! there is at least one maximum of the almost interrupted type          
             type='INTEMAX' 
          ELSE 
             ! interrputed, but distance derivative does not change???               
             ! could have one maximum and one minimum                                
             type='INTEUNK' 
          ENDIF
       ELSE 
          ! case abolished                                                        
          ! dalpha is between 90 and 90 degrees, can this happen?                 
          type='UNKNOWN' 
       ENDIF
    ENDIF
  END SUBROUTINE falslog4tp
  !===================================================================================
  SUBROUTINE achillestp(va_trace,siglim,type,                  &
       &     iunwar0,iunnew0,va_tracemin,niter,fals_conv,fals_notp,deltasig,limit_stretch) 
    USE output_control  
    USE planet_masses,   ONLY: gmearth
    !========================== INPUT================================
    TYPE(tp_point), INTENT(IN) :: va_trace
    DOUBLE PRECISION, INTENT(IN) :: siglim 
    CHARACTER(LEN=7), INTENT(INOUT) :: type 
    INTEGER, INTENT(IN) :: iunwar0,iunnew0 
    DOUBLE PRECISION, INTENT(IN) :: deltasig
    DOUBLE PRECISION, INTENT(IN) :: limit_stretch
    !========================== OUTPUT===============================     
    TYPE(tp_point), INTENT(OUT) :: va_tracemin                
    LOGICAL, INTENT(OUT) :: fals_notp,fals_conv 
    INTEGER, INTENT(OUT) :: niter 
    !========================END INTERFACE===========================
    CHARACTER(LEN=4) :: method 
    DOUBLE PRECISION :: strlov,alphalov,csi,zeta,tc 
    INTEGER :: j,it,nn, iunwar,iunnew 
    INTEGER, PARAMETER :: nit=8 
    DOUBLE PRECISION, DIMENSION(nit+1) :: xx
    TYPE(tp_point), DIMENSION(nit+1) :: ar 
    DOUBLE PRECISION, DIMENSION(nit+1) ::f,g,ddx,xlim,d 
    TYPE(tp_point), DIMENSION(2) :: va2_trace       ! for call to findminctp 
    DOUBLE PRECISION :: dddx,cont_dist,cont_prob,fm,dds,           &
         &     ddw,cont_dist1,cont_dist2               ! convergence control      
    DOUBLE PRECISION, PARAMETER :: eps_fal=1.d-10    ! this is the limit probability 
    ! (with uniform prob.density)     
    DOUBLE PRECISION, PARAMETER :: del_fal=2.3d-3    ! this is the limit change 
    ! in distance in RE (=1e-7 AU)
    DOUBLE PRECISION  :: b_e,v,mu
    !=================================================================
    WRITE(*,*)' achillestp entry ', type,va_trace%tcla,va_trace%rindex 
    iunwar=abs(iunwar0)
    iunnew=abs(iunnew0)
    WRITE(iunwar,*)' beginning Achillestp: ' 
    WRITE(iunwar,*)' x,t,f',va_trace%rindex,va_trace%tcla,va_trace%dd2_ds 
    ar(1)=va_trace
    IF(.not.prob_sampl.and.deltasig.ne.delta_sigma) THEN
       WRITE(*,*) 'achillestp : deltasig ', deltasig, ' not equal to delta_sigma', delta_sigma
       STOP
    END IF
    xlim(1)=siglim/deltasig
    xx(1)=va_trace%rindex 
    ! iterations                                                            
    nn=0 
    DO 1 it=1,nit 
       nn=nn+1 
       ! computation of f,g; uses the stretching as given in input, the LOV one
       strlov=ar(nn)%stretch_lov 
       alphalov=ar(nn)%alpha_lov 
       csi=ar(nn)%opik%coord(4) 
       zeta=ar(nn)%opik%coord(5) 
       tc=ar(nn)%tcla 
       d(nn)=ar(nn)%b
       f(nn)=ar(nn)%dd2_ds 
       g(nn)=2.d0*strlov**2 
       ! newton method                                                         
       dds=-f(nn)/g(nn) 
       ddx(nn)=dds/deltasig
       IF(abs(ddx(nn)).gt.xlim(nn))THEN 
          ddx(nn)=ddx(nn)*abs(xlim(nn)/ddx(nn)) 
       ENDIF
       xx(nn+1)=xx(nn)+ddx(nn) 
       WRITE(iunwar,*)nn,ddx(nn),xlim(nn) 
       ! find and propagate to close approach the interpolated orbit
       CALL lovclosapptp(xx(nn+1),tc,tc,iunwar,fals_notp,ar(nn+1))
       IF(fals_notp)THEN 
          WRITE(iunwar,*)' achillestp: no TP',ddx(nn),xx(nn),it 
          !          WRITE(*,*)' achillestp: no TP',ddx(nn),xx(nn),it 
          nn=nn-1 
          method='ANTP' 
          niter=it 
          CALL  wriwarn(ar(nn+1)%tcla,ar(nn+1)%rindex,ar(nn+1)%b, &
               &           ar(nn+1)%dd2_ds,                                       &
               &           ar(nn+1)%stretch,ar(nn+1)%width,ar(nn+1)%moid,       &
               &           0.d0,0.d0,niter,type,method,iunwar)          
          fals_conv=.false. 
          xlim(nn+1)=min(xlim(nn+1)/2.d0,abs(ddx(nn+1)/2.d0)) 
       ELSE 
          ! convergence criteria                                                  
          ! intervals for control                                                 
          fm=abs(ar(nn+1)%dd2_ds) 
          d(nn+1)=ar(nn+1)%b 
          dddx=abs(xx(nn+1)-xx(nn)) 
          cont_dist1=fm*deltasig*dddx/(d(nn+1)+d(nn))
          cont_dist2=abs(d(nn+1)-d(nn)) 
          IF(it.eq.1)THEN 
             cont_dist=cont_dist1 
          ELSE 
             cont_dist=min(cont_dist1,cont_dist2) 
          ENDIF
          cont_prob=dddx*deltasig/6 
          ! write warn record                                                     
          method='Ac' 
          niter=nn 
          CALL  wriwarn(ar(nn+1)%tcla,ar(nn+1)%rindex,ar(nn+1)%b, &
               &           ar(nn+1)%dd2_ds,                                       &
               &           ar(nn+1)%stretch,ar(nn+1)%width,ar(nn+1)%moid,       &
               &           cont_dist,cont_prob,niter,type,method,iunwar0)            
          ddw=d(nn+1) 
          ! Convergence control
          v   = ar(nn+1)%opik%coord(1) 
          mu  = gmearth/reau**3
          b_e = sqrt(1.d0+(2.d0*mu)/v**2)
          IF(cont_dist.LT.del_fal .OR. cont_dist.LT.del_fal*ddw            &
               &           .OR. cont_prob.LT.eps_fal .OR. ddw.LT.b_e)THEN           
             ! convergence achieved                                                  
             WRITE(iunwar,*)' minimum found, iter= ',it,nn 
             WRITE(*,*)' minimum found, iter= ',it,nn 
             ! write new record                                                      
             CALL  wriwarn(ar(nn+1)%tcla,ar(nn+1)%rindex,ar(nn+1)%b, &
                  &           ar(nn+1)%dd2_ds,                                       &
                  &           ar(nn+1)%stretch,ar(nn+1)%width,ar(nn+1)%moid,       &
                  &           cont_dist,cont_prob,niter,type,method,iunnew0)      
             ! copy output close approach array                                      
             va_tracemin=ar(nn+1) 
             fals_conv=.true. 
             RETURN 
          ELSE 
             !------------------------------------!
             ! See if possible to switch to falsi !
             !------------------------------------!
             ! If limit_stretch is not 0, check ar(nn+1)%stretch
             IF(limit_stretch.GT.0.d0)THEN
                IF(ar(nn+1)%stretch.GT.limit_stretch .OR. &
                     ar(nn+1)%stretch*ar(nn+1)%width.GT.ar(nn+1)%b**2/2.d0*1.d11)THEN
                   niter       = nn
                   fals_conv   = .FALSE.
                   va_tracemin = ar(nn+1) 
                   WRITE(*,*) ' achillestp: stretching exceeds max value'
                   RETURN
                END IF
             END IF
             IF(xx(nn).lt.xx(nn+1))THEN 
                IF(ar(nn)%dd2_ds.lt.0.d0.and.ar(nn+1)%dd2_ds.gt.0.d0)THEN 
                   WRITE(iunwar,*)' achillestp: minimum identified' 
                   ! copy arrays in the order nn, nn+1]
                   va2_trace(1)=ar(nn)
                   va2_trace(2)=ar(nn+1)                                     
                   CALL falslog4tp(va2_trace,iunwar,type) 
                   ! new type                                                              
                   IF(type.eq.'SIMPMIN')THEN 
                      type='ASIMPMI' 
                   ELSEIF(type.eq.'INTEMIN')THEN 
                      type='AINTEMI' 
                      ! use falsi only for interrupted                                        
                      WRITE(iunwar,*)' achillestp: switch to falsi' 
                      CALL findminctp(va2_trace,xx(nn),xx(nn+1),type,        &
                           &                       iunwar,iunnew,va_tracemin,                      &
                           &                       fals_conv,niter,fals_notp,deltasig)                 
                      RETURN 
                   ENDIF
                ENDIF
             ELSE 
                IF(ar(nn+1)%dd2_ds.lt.0.d0.and.ar(nn)%dd2_ds.gt.0.d0)THEN 
                   WRITE(iunwar,*)' achillestp: minimum identified' 
                   ! copy arrays in the order nn+1,nn                                      
                   va2_trace(1)=ar(nn+1)
                   va2_trace(2)=ar(nn)
                   CALL falslog4tp(va2_trace,iunwar,type) 
                   ! new type                                                              
                   IF(type.eq.'SIMPMIN')THEN 
                      type='ASIMPMI' 
                   ELSEIF(type.eq.'INTEMIN')THEN 
                      type='AINTEMI' 
                      ! use falsi only for interrupted                                        
                      WRITE(iunwar,*)' achillestp: switch to falsi' 
                      CALL findminctp(va2_trace,xx(nn+1),xx(nn),type,        &
                           &                       iunwar,iunnew,va_tracemin,             &
                           &                       fals_conv,niter,fals_notp,deltasig)                 
                      RETURN 
                   ENDIF
                ENDIF
             ENDIF
          ENDIF
          xlim(nn+1)=xlim(nn)/2.d0 
       ENDIF
1   ENDDO
    ! case not convergent                                                   
    IF(fals_notp)THEN 
       method='ANTP' 
       CALL  wriwarn(ar(nn+1)%tcla,xx(nn+1),ar(nn+1)%b,ar(nn+1)%dd2_ds,         &
            &     ar(nn+1)%stretch,ar(nn+1)%width,ar(nn+1)%moid,                 &
            &     0.d0,0.d0,niter,type,method,iunwar)
    ELSE 
       method='AcNC' 
       CALL  wriwarn(ar(nn+1)%tcla,xx(nn+1),ar(nn+1)%b,ar(nn+1)%dd2_ds,         &
            &     ar(nn+1)%stretch,ar(nn+1)%width,ar(nn+1)%moid,                 &
            &     cont_dist,cont_prob,niter,type,method,iunnew0)                
       CALL  wriwarn(ar(nn+1)%tcla,xx(nn+1),ar(nn+1)%b,ar(nn+1)%dd2_ds,         &
            &     ar(nn+1)%stretch,ar(nn+1)%width,ar(nn+1)%moid,                 &
            &     cont_dist,cont_prob,niter,type,method,iunwar) 
    ENDIF
    niter=nn

    ! this case not convergent, but the output can be used                  
    fals_conv=.false.
    va_tracemin=ar(nn+1) 
    RETURN 
  END SUBROUTINE achillestp
  !===============================================================
  ! ==============================================================        
  ! FINDMINCTP  Regula falsi, to be used only when the                      
  !           two previous points have opposite sign                      
  ! ============================================================          
  !                                                                       
  SUBROUTINE findminctp(arrc,a0,b0,type,iunwar0,iunnew0,arrmin,         &
       &            fals_conv,niter,fals_notp,deltasig) 
    USE planet_masses
    !========================== INPUT===============================
    DOUBLE PRECISION, INTENT(IN)             :: a0,b0 
    TYPE(tp_point), DIMENSION(2), INTENT(IN) :: arrc   ! the last 
    !two used in the return  
    CHARACTER(LEN=7), INTENT(IN)             :: type
    DOUBLE PRECISION, INTENT(IN)             :: deltasig
    INTEGER, INTENT(IN)                      :: iunwar0,iunnew0 
    !======================= INPUT/OUTPUT===========================
    !    TYPE(tp_point), DIMENSION(2) :: arrf       ! array contains the last two 
    ! used in this routine   
    !========================== OUTPUT==============================
    TYPE(tp_point), INTENT(OUT)  :: arrmin
    INTEGER, INTENT(OUT)         :: niter 
    LOGICAL, INTENT(OUT)         :: fals_conv,fals_notp 
    !======================= END INTERFACE==========================
    INTEGER :: it,nit 
    !    INTEGER, PARAMETER :: itmax=200,itma1=itmax+1,itma2=itmax+2 
    INTEGER, PARAMETER :: itmax=30,itma1=itmax+1,itma2=itmax+2 
    DOUBLE PRECISION, DIMENSION(itma1) :: a,b,ta,tb
    DOUBLE PRECISION, DIMENSION(itma1) :: da,db,fa,fb 
    DOUBLE PRECISION, DIMENSION(itma1) :: fw,tw,dw
    DOUBLE PRECISION, DIMENSION(itma2) :: w
    DOUBLE PRECISION, DIMENSION(itma1) :: wstr,wwid,wmoidgf 
    DOUBLE PRECISION :: ffa,ffb 
    CHARACTER(LEN=7) :: type2 
    DOUBLE PRECISION :: ddx,dold,cont_dist,cont_prob,fm  &
         &     ,ddw,cont_dist1,cont_dist2       ! convergence control       
    DOUBLE PRECISION, PARAMETER :: eps_fal=1.d-10 ! limit probability 
    ! (with uniform prob.density)
    DOUBLE PRECISION, PARAMETER :: del_fal=2.3d-3 ! limit change in 
    ! distance in RE (=1e-7 AU) 
    DOUBLE PRECISION :: tafter 
    LOGICAl          :: falsok
    CHARACTER(LEN=4) :: method 
    INTEGER          :: i, npset, iunwar, iunnew
    DOUBLE PRECISION  :: b_e,vsize,v,mu
    REAL(KIND=qkind) :: fbb,faa,aa,bb,xx
    !=========================================================================
    iunwar=abs(iunwar0)
    iunnew=abs(iunnew0) 
    WRITE(iunwar,*)'findminctp ',type,arrc(1)%tcla,a0,b0 
    IF(a0.eq.b0)THEN
       WRITE(iunwar,*) ' supect multiple minimum '
       fals_conv=.false.
       niter=0       
       fals_notp=.false.
       arrmin=arrc(1)
       RETURN
    ENDIF
    !    arrf(1)=arrc(1)
    !    arrf(2)=arrc(2)
    ! consult arrays for initial two points                                 
    a(1)=a0 
    b(1)=b0 
    ta(1)=arrc(1)%tcla 
    tb(1)=arrc(2)%tcla 
    da(1)=arrc(1)%b 
    db(1)=arrc(2)%b 
    fa(1)=arrc(1)%dd2_ds 
    fb(1)=arrc(2)%dd2_ds 
    ! store a as w(1)                                                       
    w(1)=a0 
    fw(1)=fa(1) 
    dw(1)=da(1) 
    tw(1)=ta(1) 
    wstr(1)=arrc(1)%stretch 
    wwid(1)=arrc(1)%width 
    wmoidgf(1)=arrc(1)%moid 
    ! beginning record                                                      
    WRITE(iunwar,120)ta(1),a(1),da(1),fa(1),tb(1),b(1),db(1),fb(1) 
120 FORMAT('beginning falsi:'/'t1,x1,d1,f1 ',f10.2,1x,f11.5,1p,2d11.3,&
         &  0p/'t2,x2,d2,f2 ',f10.2,1x,f11.5,1p,2d11.3)                     
    ! control on sign                                                       
    IF(fa(1)*fb(1).gt.0.d0)THEN 
       WRITE(iunwar,*)' this case not suitable for regula falsi' 
    ENDIF
    ! first point to be tested                                              
    w(2)=(fb(1)*a(1)-fa(1)*b(1))/(fb(1)-fa(1)) 
    ! iteartion loop                                                        
    DO 1 it=1,itmax 
       npset=0
       ! find and propagate to close approach the interpolated orbit
11     CALL lovclosapptp(w(it+1),ta(it),tb(it),iunwar,          &
            &        fals_notp,arrmin)                                      
       ! target plane found?                                                   
       IF(fals_notp)THEN 
          fals_conv=.false. 
          IF(it.eq.1)THEN 
             cont_dist=-1.d0 
             cont_prob=-1.d0 
          ENDIF
          niter=it 
          write(iunwar,*)' falsi failed: no target plane', w(it+1) 
          method='FNTP' 
          CALL wriwarn(tw(it),w(it),dw(it),fw(it),wstr(it),wwid(it),  &
               &        wmoidgf(it),cont_dist,cont_prob,niter,type,method,iunwar) 
          CALL wriwarn(tw(it),w(it),dw(it),fw(it),wstr(it),wwid(it),  &
               &        wmoidgf(it),cont_dist,cont_prob,niter,type,method,iunnew) 
          RETURN 
       ELSEIF(.not.arrmin%tp_conv)THEN
          ! set a larger value for npoint
          npset=npset+1
          IF(npset.gt.0)THEN
             niter=it 
             write(iunwar,*)' falsi failed: singularity', w(it+1) 
             write(iunnew,*)' falsi failed: singularity', w(it+1) 
             method='FSIN' 
             tw(it+1)=arrmin%tcla
             dw(it+1)=arrmin%b
             wstr(it+1)=arrmin%stretch 
             wwid(it+1)=arrmin%width
             fw(it+1)=arrmin%dd2_ds
             wmoidgf(it+1)=arrmin%moid
             CALL wriwarn(tw(it+1),w(it+1),dw(it+1),fw(it+1),wstr(it+1),wwid(it+1),  &
                  &             wmoidgf(it+1),cont_dist,cont_prob,niter,type,method,iunwar) 
             CALL wriwarn(tw(it+1),w(it+1),dw(it+1),fw(it+1),wstr(it+1),wwid(it+1),  &
                  &                wmoidgf(it+1),cont_dist,cont_prob,niter,type,method,iunnew)
             arrmin%rindex=w(it+1)
             RETURN
          ELSE
             write(iunwar,*)' findminctp: npoint increased by ',2**npset
             write(iunnew,*)' findminctp: npoint increased by ',2**npset
             !             CALL npoint_set(2**npset) ! try increasing npoint by a factor 2
             GOTO 11
          ENDIF
       ENDIF
       ! data on the new point                                                 
       tw(it+1)=arrmin%tcla
       dw(it+1)=arrmin%b
       wstr(it+1)=arrmin%stretch 
       wwid(it+1)=arrmin%width
       fw(it+1)=arrmin%dd2_ds
       wmoidgf(it+1)=arrmin%moid 
       ! falsi logic                                                           
       IF(fa(it)*fw(it+1).le.0.d0.and.fb(it)*fw(it+1).ge.0.d0)THEN 
          ! standard regula falsi                                                 
          b(it+1)=w(it+1) 
          fb(it+1)=fw(it+1) 
          db(it+1)=dw(it+1) 
          tb(it+1)=tw(it+1) 
          !...........................
          a(it+1)=a(it) 
          fa(it+1)=fa(it) 
          da(it+1)=da(it) 
          ta(it+1)=ta(it) 
          ! improved regula falsi                                                 
          IF(fw(it)*fw(it+1).gt.0.d0.and.it.gt.1)THEN 
             IF(it.gt.1.and.fw(it-1)*fw(it).gt.0.d0)THEN 
                ! superaccelerator                                                      
                ffa=fa(it)/4.d0 
             ELSE 
                ! accelerator                                                           
                ffa=fa(it)/2.d0 
             ENDIF
          ELSE 
             ! ordinary regula falsi                                                 
             ffa=fa(it) 
          ENDIF
          ffb=fb(it) 
       ELSEIF(fb(it)*fw(it+1).le.0.d0.and.fa(it)*fw(it+1).ge.0.d0)THEN 
          ! standard regula falsi                                                 
          a(it+1)=w(it+1) 
          fa(it+1)=fw(it+1) 
          da(it+1)=dw(it+1) 
          ta(it+1)=tw(it+1) 
          b(it+1)=b(it) 
          fb(it+1)=fb(it) 
          db(it+1)=db(it) 
          tb(it+1)=tb(it) 
          ! improved regula falsi                                                 
          IF(fw(it)*fw(it+1).gt.0.d0.and.it.gt.1)THEN 
             IF(it.gt.1.and.fw(it-1)*fw(it).gt.0.d0)THEN 
                ! superaccelerator                                                      
                ffb=fb(it)/4.d0 
             ELSE 
                ! accelerator                                                           
                ffb=fb(it)/2.d0 
             ENDIF
          ELSE 
             ! ordinary regula falsi                                                 
             ffb=fb(it) 
          ENDIF
          ffa=fa(it) 
       ELSE 
          ! these are not regula falsi, but the secant method                     
          WRITE(iunwar,*)' findminctp should not be used with same sign' 
          WRITE(*,*) it, a(it),fa(it),b(it),fb(it),w(it+1),fw(it+1) 
       ENDIF
       ! compute next step                                                     
       w(it+2)=(ffb*a(it+1)-ffa*b(it+1))/(ffb-ffa) 
       ! convergence criteria                                                  
       ! intervals for control                                                 
       fm=abs(fw(it+1)) 
       ddx=abs(b(it+1)-a(it+1)) 
       IF(.not.prob_sampl.and.deltasig.ne.delta_sigma) THEN
          WRITE(*,*) 'findminctp : deltasig ', deltasig, ' not equal to delta_sigma', delta_sigma
          STOP
       END IF
       cont_dist1=fm*deltasig*ddx/(da(it+1)+db(it+1)) 
       cont_dist2=abs(dw(it+1)-dw(it)) 
       IF(it.eq.1)THEN 
          cont_dist=cont_dist1 
       ELSE 
          cont_dist=min(cont_dist1,cont_dist2) 
       ENDIF
       cont_prob=ddx*deltasig/6 
       niter=it 
       method='Fa' 
       CALL wriwarn(tw(it+1),w(it+1),dw(it+1),fw(it+1),wstr(it+1),    &
            &        wwid(it+1),wmoidgf(it+1),cont_dist,cont_prob,niter,       &
            &        type,method,iunwar)                                       
       ddw=dw(it+1) 
       ! check type           
       CALL falslog4tp(arrc,iunwar,type2) 
       IF(type2.ne.type)THEN 
          IF(type.eq.'AINTEMI'.and.type2.eq.'INTEMIN')THEN 
          ELSEIF(type.eq.'ASIMPMI'.and.type2.eq.'SIMPMIN')THEN 
          ELSE 
             WRITE(iunwar,*)' type change, from ',type,' to ',type2 
          ENDIF
       ENDIF
       ! convergence control
       !         tp_flag_cur=arrmin%tp_conv
       !         IF(tp_flag_cur)THEN 
       ! on TP, gravitational focusing 
       v=arrmin%opik%coord(1) 
       mu=gmearth/reau**3
       b_e=sqrt(1.d0+(2.d0*mu)/v**2)
       IF(cont_dist.lt.del_fal  & ! .or.cont_dist.lt.del_fal*ddw      &
            &           .or.cont_prob.lt.eps_fal.or.ddw.lt.b_e)THEN           
          WRITE(iunwar,*)' minimum found, iter= ',niter 
          WRITE(*,*)' minimum found, iter= ',niter 
          fals_conv=.true. 
          ! write record in newton file                                           
          CALL wriwarn(tw(it+1),w(it+1),dw(it+1),fw(it+1),wstr(it+1), &
               &           wwid(it+1),wmoidgf(it+1),cont_dist,cont_prob,niter,    &
               &           type,method,iunnew)                                    
          arrmin%rindex=w(it+1) 
          RETURN 
       ENDIF
1   ENDDO
    ! non convergent case                                                   
    niter=itmax 
    WRITE(iunwar,*)' too many iter, ',cont_dist,cont_prob,itmax 
    IF(dw(it).lt.1.d0)THEN 
       fals_conv=.true. 
    ELSE 
       fals_conv=.false. 
    ENDIF
    ! write record in newton file                                           
    nit=niter+1 
    method='FaNC' 
    CALL wriwarn(tw(nit),w(nit),dw(nit),fw(nit),wstr(nit),wwid(nit),  &
         &   wmoidgf(nit),cont_dist,cont_prob,niter,type,method,iunwar)     
    CALL wriwarn(tw(nit),w(nit),dw(nit),fw(nit),wstr(nit),wwid(nit),  &
         &   wmoidgf(nit),cont_dist,cont_prob,nit,type,method,iunnew)     
    arrmin%rindex=w(it+1) 
  END SUBROUTINE findminctp
  ! ==========================================================================
  ! NOT IN USE
  ! =========================================================================
  SUBROUTINE propclosapp(jlov,t1,t2,iunwar,rindex,fals_notp,va_tracemin) 
    USE close_app, ONLY: kill_propag
    USE orbit_elements
    USE propag_state
    USE multiple_sol
    USE dyn_param
    ! ========================= INPUT========================================
    INTEGER, INTENT(IN)          :: jlov   ! index in array elm_loc
    DOUBLE PRECISION, INTENT(IN) ::  t1,t2 ! min, max time for propagation
    INTEGER, INTENT(IN)          :: iunwar 
    DOUBLE PRECISION, INTENT(IN) :: rindex ! real LOV index
    ! ==========================OUTPUT==========================================
    LOGICAL, INTENT(OUT) :: fals_notp 
    TYPE(tp_point), INTENT(OUT) :: va_tracemin        
    ! =======================END INTERFACE======================================
    LOGICAL :: falsok 
    DOUBLE PRECISION :: tafter,tbefore           ! epoch time
    TYPE(orbit_elem) :: el0,el
    TYPE(orb_uncert) :: unc0,unc
    INTEGER :: ipla, nd 
    DOUBLE PRECISION :: v_infty, v_inf0             ! velocity at the boundary 
    ! of Earth's sphere 
    ! of influence
    LOGICAL :: batch
    ! =========================================================================== 
    nd=6+nls
    fals_notp=.false. 
    ! find orbit with interpolated value of sigma   
    el0=elm_loc(jlov)
    ! dynamical parameters
    if(nd.gt.6) dyn%dp=dpm_loc(:,jlov)
    ! select propagation time                                               
    ipla=3 
    v_inf0=v_infty(el0) 
    CALL aftclo2v(ipla,el0%t,t1,t2,v_inf0,tbefore,tafter) 
    ! availability of covariance for TP analysis online: NO                     
    CALL cov_not_av
    ! reset close approach storage                                          
    njc=0 
    iplam=0 
    ! propagate                                                             
    batch=.true. 
    CALL pro_ele(el0,tafter,el) 
    !  WRITE(*,*)'lovclosapptp: t1,t2,tbefore,tafter,tclas: ',t1,t2,tbefore,tafter,tcla(1:njc)
    IF(njc.eq.0)THEN 
       WRITE(iunwar,*)' lovclosapp failed; no close approach' 
       falsok=.false. 
    ELSEIF(iplam.ne.3)THEN 
       WRITE(iunwar,*)' lovclosapp failed; close app. planet ', iplam
       falsok=.false. 
       kill_propag=.false.
    ELSEIF(tp_store(njc)%tcla.gt.max(tbefore,tafter).or.tp_store(njc)%tcla.lt.min(tbefore,tafter))THEN   
       IF(kill_propag)THEN
          WRITE(iunwar,*)' lovclosapp failed: previous collision at ',tcla(njc)
       ELSE
          WRITE(iunwar,*)' lovclosapp failed; tcla too far ',tcla(njc) 
       ENDIF
       falsok=.false.
       kill_propag=.false.
    ELSE 
       IF(kill_propag)THEN
          WRITE(iunwar,*)' lovclosapp stopped, collision found '
       ENDIF
       kill_propag=.false.
       falsok=.true.
    ENDIF
    IF(falsok)THEN      
       CALL arrloadtp(va_tracemin,rindex)
       IF(.not.va_tracemin%tp_conv)THEN
          WRITE(iunwar,*)' lovclosapp failed: close approach elliptic ',va_tracemin%tcla
          fals_notp=.true.
       ENDIF
    ELSE 
       fals_notp=.true. 
    ENDIF
  END SUBROUTINE propclosapp
  ! =========================================================================
  SUBROUTINE lovclosapptp(rindex,t1,t2,iunwar,fals_notp,va_tracemin) 
    USE close_app, ONLY: kill_propag
    USE orbit_elements
    USE propag_state
    USE multiple_sol, ONLY: lovinterp, lovinterp3
    ! ========================= INPUT========================================
    DOUBLE PRECISION, INTENT(IN) :: rindex 
    DOUBLE PRECISION, INTENT(IN) ::  t1,t2
    INTEGER, INTENT(IN) :: iunwar 
    ! ==========================OUTPUT==========================================
    LOGICAL, INTENT(OUT) :: fals_notp 
    TYPE(tp_point), INTENT(OUT) :: va_tracemin        
    ! =======================END INTERFACE======================================
    LOGICAL :: falsok 
    DOUBLE PRECISION :: tafter,tbefore           ! epoch time
    TYPE(orbit_elem) :: el0,el
    TYPE(orb_uncert) :: unc0,unc
    INTEGER :: ipla 
    DOUBLE PRECISION :: v_infty, v_inf0             ! velocity at the boundary 
    ! of Earth's sphere 
    ! of influence
    DOUBLE PRECISION :: sigmaout   ! value of sigma_lov at the interpolated point
    ! =========================================================================== 
    fals_notp=.false. 
    ! find orbit with interpolated value of sigma                           
    IF(prob_sampl)THEN
       CALL lovinterp3(rindex,sigmaout,el0,unc0,falsok)
    ELSE
       CALL lovinterp(rindex,delta_sigma,el0,unc0,falsok) 
    ENDIF
    IF(.not.falsok) THEN 
       WRITE(iunwar,*)' lovinterp fails', rindex 
       fals_notp=.true. 
       RETURN 
    ENDIF
    ! select propagation time                                               
    ipla=3 
    v_inf0=v_infty(el0) 
    CALL aftclo2v(ipla,el0%t,t1,t2,v_inf0,tbefore,tafter) 
    ! availability of covariance for TP analysis online                     
    CALL cov_avai(unc0,el0%coo,el0%coord) 
    ! reset close approach storage                                          
    njc=0 
    iplam=0 
    ! propagate                                                              
    CALL pro_ele(el0,tafter,el,unc0,unc) 
    !  WRITE(*,*)'lovclosapptp: t1,t2,tbefore,tafter,tclas: ',t1,t2,tbefore,tafter,tcla(1:njc)
    IF(njc.eq.0)THEN 
       WRITE(iunwar,*)' lovclosapp failed; no close approach' 
       falsok=.false. 
    ELSEIF(iplam.ne.3)THEN 
       WRITE(iunwar,*)' lovclosapp failed; close app. planet ', iplam
       falsok=.false. 
       kill_propag=.false.
    ELSEIF(tp_store(njc)%tcla.gt.max(tbefore,tafter).or.tp_store(njc)%tcla.lt.min(tbefore,tafter))THEN   
       IF(kill_propag)THEN
          WRITE(iunwar,*)' lovclosapp failed: previous collision at ',tcla(njc)
       ELSE
          WRITE(iunwar,*)' lovclosapp failed; tcla too far ',tcla(njc) 
       ENDIF
       falsok=.false.
       kill_propag=.false.
    ELSE 
       IF(kill_propag)THEN
          WRITE(iunwar,*)' lovclosapp stopped, collision found '
       ENDIF
       kill_propag=.false.
       falsok=.true.
    ENDIF
    IF(falsok)THEN 
       IF(prob_sampl)THEN
          CALL arrloadtp(va_tracemin,rindex,sigmaout)
       ELSE
          CALL arrloadtp(va_tracemin,rindex)
       ENDIF
       IF(.not.va_tracemin%tp_conv)THEN
          WRITE(iunwar,*)' lovclosapp failed: close approach elliptic ',va_tracemin%tcla
          fals_notp=.true.
       ENDIF
    ELSE 
       fals_notp=.true. 
    ENDIF
  END SUBROUTINE lovclosapptp
  ! ================================================================      
  ! SUBROUTINE aftclo2v                                                  
  ! select time interval to get after the close approach                  
  ! taking into account the relative velocity w.r. to Earth
  ! ================================================================        
  SUBROUTINE aftclo2v(iplam,t0,t1,t2,v_inf0,tbefore,tafter) 
    USE planet_masses
    ! ======================INPUT===================================== 
    INTEGER, INTENT(IN) :: iplam             ! planet number
    DOUBLE PRECISION, INTENT(IN) :: t0       ! time of initial 
    ! conditions
    DOUBLE PRECISION, INTENT(IN) :: t1,t2    ! time of two known 
    ! close approaches
    DOUBLE PRECISION, INTENT(IN) :: v_inf0   ! velocity 
    ! =====================OUTPUT=====================================
    DOUBLE PRECISION, INTENT(OUT) :: tafter    ! time "after"
    DOUBLE PRECISION, INTENT(OUT) :: tbefore   ! time "before"
    !  depending upon sense of propagation  
    ! ===================END INTERFACE=================================      
    DOUBLE PRECISION :: delt_tp       ! time interval to 
    ! exit from TP disk    
    ! ===================================================================   
    ! warning: really done only for Earth                                   
    IF(ordnam(iplam).ne.'EARTH') THEN 
       WRITE(*,*)' aftclo2v: not to be used for planet ',ordnam(iplam) 
    ENDIF
    ! time interval to exit from TP disk   
    IF(v_inf0.gt.0.d0)THEN                                 
       delt_tp=2*dmin(iplam)/v_inf0
    ELSE
       delt_tp=365.25d0
    ENDIF
    ! forced to avoid infinte intervals for v-inf=0
    delt_tp=MIN(delt_tp,365.25d0)
    ! forced to avoid short intervals for fast encounters                   
    IF(delt_tp.lt.20.d0)delt_tp=20.d0 
    ! time interval to be clear out of the TP disk is given as deltat       
    IF(t0.lt.t1)THEN 
       ! future close approaches                                               
       tafter=max(t1,t2)+delt_tp 
       tbefore=min(t1,t2)-delt_tp 
    ELSE 
       ! past close approaches                                                 
       tafter=min(t1,t2)-delt_tp 
       tbefore=max(t1,t2)+delt_tp 
    ENDIF
  END SUBROUTINE aftclo2v
  ! =========================================================================







  ! ===========================================================
  ! SMALL I/O ROUTINES
  ! ===========================================================

  SUBROUTINE wrireptp(vatr,imul1,imul2,iunrep,no_outcov)
    ! =================================================================
    INTEGER, INTENT(IN)  :: iunrep ! output unit
    DOUBLE PRECISION, INTENT(IN) :: imul1,imul2 ! return interval 
    TYPE(tp_point), INTENT(IN) :: vatr ! tp point with minimum distance
    LOGICAL, INTENT(IN), OPTIONAL :: no_outcov
    CHARACTER(LEN=14) :: calend                ! calendar date
    INTEGER iun
    LOGICAL nooutcov
    ! =================================================================
    IF(PRESENT(no_outcov))THEN
       nooutcov=no_outcov
    ELSE
       nooutcov=.false.
    ENDIF
    CALL calendwri(vatr%tcla,calend) 
    ! write .rep file
    IF(nooutcov)THEN
       IF(iunrep.ge.0)THEN
          WRITE(iunrep,101)calend,vatr%tcla,vatr%rindex,imul1,imul2,vatr%b,  &
               &     vatr%opik%coord(4),vatr%opik%coord(5),vatr%minposs, &
               &     vatr%opik%coord(1)
       ELSEIF(iunrep.lt.0)THEN
          iun=abs(iunrep)
          WRITE(iun,101)calend,vatr%tcla,vatr%rindex,imul1,imul2,vatr%b,  &
               &     vatr%opik%coord(4),vatr%opik%coord(5),vatr%minposs, &
               &     vatr%opik%coord(1)
          WRITE(*,101)calend,vatr%tcla,vatr%rindex,imul1,imul2,vatr%b,  &
               &     vatr%opik%coord(4),vatr%opik%coord(5),vatr%minposs,&
               &     vatr%opik%coord(1)
       ENDIF
    ELSE
       IF(iunrep.ge.0)THEN
          WRITE(iunrep,100)calend,vatr%tcla,vatr%rindex,imul1,imul2,vatr%b,  &
               &     vatr%opik%coord(4),vatr%opik%coord(5),vatr%stretch,  &
               &     vatr%width,vatr%alpha*degrad,vatr%minposs,      &
               &     vatr%stretch_lov,vatr%alpha_lov*degrad,vatr%opik%coord(1)
       ELSEIF(iunrep.lt.0)THEN
          iun=abs(iunrep)
          WRITE(iun,100)calend,vatr%tcla,vatr%rindex,imul1,imul2,vatr%b,  &
               &     vatr%opik%coord(4),vatr%opik%coord(5),vatr%stretch,  &
               &     vatr%width,vatr%alpha*degrad,vatr%minposs,      &
               &     vatr%stretch_lov,vatr%alpha_lov*degrad,vatr%opik%coord(1)
          WRITE(*,100)calend,vatr%tcla,vatr%rindex,imul1,imul2,vatr%b,  &
               &     vatr%opik%coord(4),vatr%opik%coord(5),vatr%stretch,  &
               &     vatr%width,vatr%alpha*degrad,vatr%minposs,      &
               &     vatr%stretch_lov,vatr%alpha_lov*degrad,vatr%opik%coord(1)
       ENDIF
100    FORMAT(a14,1x,f12.5,1x,f6.1,1x,f6.1,1x,f6.1,1x,3(1x,f12.6),1p,          &
            &     2(1x,e10.3),0p,1x,f8.3,1x,f12.6,1x,1p,e11.4,0p,1x,f10.5,     &
            &     1p,e12.4)
101    FORMAT(a14,1x,f12.5,1x,f6.1,1x,f6.1,1x,f6.1,1x,3(1x,f12.6),             &
            &     1x,f12.6,1x,1p,e12.4)               
    ENDIF
  END SUBROUTINE wrireptp

  !=============================================================================
  ! SUBROUTINE HEADER_REP
  !=============================================================================

  SUBROUTINE header_rep(iunrep,no_outcov)
    INTEGER, INTENT(IN) :: iunrep ! output unit
    LOGICAL, INTENT(IN), OPTIONAL :: no_outcov
    LOGICAL nooutcov
    !=============================================================================
    IF(PRESENT(no_outcov))THEN
       nooutcov=no_outcov
    ELSE
       nooutcov=.false.
    ENDIF
    IF(nooutcov)THEN
       WRITE(iunrep,198) 
198    FORMAT(' CALENDAR DATE ',1x,'  MJD DATE  ',1x,'IMIN'                   &
            &     ,1x,'   IMU1 ',1x,' IMU2',1x,'  MIN_FOUND ',1x                     &
            &           ,1x,'  CSI (RE) ',1x,'  ZETA (RE) ',1x,'      MOID     ',1x    &
            &                ,1x,'  ')     
    ELSE
       WRITE(iunrep,199) 
199    FORMAT('  CALENDAR DATE ',1x,'  MJD DATE  ',1x,'IMIN'             &
            &     ,1x,'IMU1',1x,'IMU2',1x,'   MIN_FOUND ',1x                     &
            &           ,1x,'  CSI (RE) ',1x,'    ZETA (RE) ',1x,'STRETCHING', 1x,&
            &  '   WIDTH  ',1x,' ANG_ELL',1x,'  MOID     ',                    &
            &        1x,'STRETC_LOV',1x,' ANG_LOV')     
    ENDIF
  END SUBROUTINE header_rep

  !=============================================================================
  ! SUBROUTINE HEADER_OUT
  !=============================================================================

  SUBROUTINE header_out(iunout,no_outcov)
    INTEGER, INTENT(IN) :: iunout ! output unit
    LOGICAL, INTENT(IN), OPTIONAL :: no_outcov
    LOGICAL nooutcov
    !=============================================================================
    IF(PRESENT(no_outcov))THEN
       nooutcov=no_outcov
    ELSE
       nooutcov=.false.
    ENDIF
    IF(nooutcov)THEN
       WRITE(iunout,188) 
188    FORMAT('RINDEX ',1x,' CALENDAR DATE ',1x,' MJD DATE  ',1x,' CSI(OPIK) ',1x,' ZETA(OPIK) ',1x,     &
            & ' CSI ',1x,' ZETA ',1x,' D ',1x,' TETA ',1x,' PHI ',1x,' MOID ')
    ELSE
       WRITE(iunout,189) 
189    FORMAT('RINDEX ',1x,' CALENDAR DATE ',1x,' MJD DATE ',1x,' CSI(OPIK) ',1x,' ZETA(OPIK) ',1x,       &
            & ' CSI ',1x,' ZETA ',1x,' STRETCH ',1x,' WIDTH ',1x,' ALPHA ',1x,' STRETCH_LOV ',1x,         &
            & ' ALP_LOV ',1x,' D ',1x,' TETA ',1x,' PHI ',1x,' MOID')
    ENDIF
  END SUBROUTINE header_out

  SUBROUTINE wriouttp(vatr,iunout,no_outcov) 
    ! =========================================================================
    INTEGER, INTENT(IN) :: iunout
    TYPE(tp_point), INTENT(IN) :: vatr
    LOGICAL, INTENT(IN), OPTIONAL :: no_outcov
    LOGICAL nooutcov  
    CHARACTER(LEN=14) :: calend                ! calendar date
    ! ====================================================================
    IF(PRESENT(no_outcov))THEN
       nooutcov=no_outcov
    ELSE
       nooutcov=.false.
    ENDIF
    CALL calendwri(vatr%tcla,calend)
    ! write .out file  
    IF(nooutcov)THEN
       WRITE(iunout,101)vatr%rindex,calend,vatr%tcla,vatr%txi,vatr%tze,                   &
            &           vatr%opik%coord(4),vatr%opik%coord(5),vatr%d,                     &
            &           vatr%theta*degrad,vatr%phi*degrad,vatr%moid
101    FORMAT(f6.1,1x,a14,1x,f12.5,2(1x,f18.12),                                          & 
            &     2(1x,f18.12),1x,f15.10,                                                 &
            &     1x,f10.6,1x,f11.6,1x,f12.6)
    ELSE                                                     
       WRITE(iunout,100)vatr%rindex,calend,vatr%tcla,vatr%txi,vatr%tze,                   &
            &           vatr%opik%coord(4),vatr%opik%coord(5),                            &
            &           vatr%stretch,vatr%width,vatr%alpha*degrad,                        &
            &           vatr%stretch_lov,vatr%alpha_lov*degrad,                           &
            &           vatr%d,vatr%theta*degrad,vatr%phi*degrad,vatr%moid 
100    FORMAT(f6.1,1x,a14,1x,f12.5,1x,2(1x,f18.12),1x,                                    &
            &     2(1x,f18.12),1x,1p,e14.6,2x,0p,f10.5,3x,f8.3,                           &
            &     1x,1p,e14.6,3x,0p,f8.3,                                                 &
            &     3x,f18.12,1x,f10.6,1x,f11.6,1x,f12.6)
    ENDIF
  END SUBROUTINE wriouttp

  !========================================================================= 
  SUBROUTINE wriwarn(tw,w,dw,fw,wstr,wwid,wmoidgf,                  &
       &     cont_dist,cont_prob,niter,type,method,iunwar)                
    ! ===========================INPUT=========================================
    DOUBLE PRECISION, INTENT(IN) :: tw,w,dw,fw,wstr,wwid,wmoidgf ! data 
    ! on current close approach 
    DOUBLE PRECISION, INTENT(IN) ::  cont_dist,cont_prob  ! convergence control 
    INTEGER, INTENT(IN) :: niter               ! number of iterations so far
    CHARACTER(LEN=7), INTENT(IN) :: type       ! string comments 
    CHARACTER(LEN=4), INTENT(IN) :: method 
    INTEGER,INTENT(IN) :: iunwar               ! output unit
    ! ==========================================================================
    CHARACTER(LEN=14) :: calend                ! calendar date 
    INTEGER iun
    ! =========================================================================
    CALL calendwri(tw,calend) 
    IF(iunwar.gt.0)THEN
       WRITE(iunwar,119)calend,tw,w,dw,fw,wstr                           &
            &     ,wwid,wmoidgf,cont_dist,cont_prob,niter,type,method
    ELSEIF(iunwar.eq.0)THEN
       CALL header_new(iunwar)
       WRITE(*,119)calend,tw,w,dw,fw,wstr                           &
            &     ,wwid,wmoidgf,cont_dist,cont_prob,niter,type,method
    ELSEIF(iunwar.lt.0)THEN
       iun=abs(iunwar)
       CALL header_new(0)
       WRITE(*,119)calend,tw,w,dw,fw,wstr                           &
            &     ,wwid,wmoidgf,cont_dist,cont_prob,niter,type,method
       WRITE(iun,119)calend,tw,w,dw,fw,wstr                           &
            &     ,wwid,wmoidgf,cont_dist,cont_prob,niter,type,method
    ENDIF
119 FORMAT(a14,1x,f12.5,1x,f10.5,1x,f7.2,1x,1p,d9.2,1x,d9.2,1x,d9.2,  &
         &           1x,d9.2,1x,d9.2,1x,d9.2,1x,i3,1x,a7,1x,a4)             
  END SUBROUTINE wriwarn

  SUBROUTINE header_new(iunnew)
    INTEGER, INTENT(IN) :: iunnew ! output unit
    WRITE(iunnew,198) 
198 FORMAT('  CALENDAR DATE ',1x,'      MJD  ',1x,'  INDEX  ',1x,     &
         &       'LOV_DIST ',1x,'DR2/DSIG',1x,'  STRETCH.', 1x,'  WIDTH ',  &
         &       '   MOIDGF',1x,'  DIST.CON',1x,' PROB.CON',1x,' IT',1x,    &
         & 'TYPE',1x,'METH')  
  END SUBROUTINE header_new

!===========================================================================
! SUBROUTINE arrcut 
!===========================================================================
!SUBROUTINE arrcut(vas_trace,ire,lre,nox,iunout,vas_traceloc,no_outcov,dmeacontr2,dnewton) 
SUBROUTINE arrcut(vas_trace,ire,lre,nox,iunout,vas_traceloc,no_outcov) 
  INTEGER, INTENT(IN) :: iunout   ! unit for out record
  INTEGER, INTENT(IN) :: ire,nox  ! address of filament, dimension
  INTEGER, INTENT(INOUT) :: lre   ! length of filament, possibly reduced in output
  TYPE(tp_point), DIMENSION(nox), INTENT(IN) :: vas_trace 
                                      ! global close approach records,
  LOGICAL, INTENT(IN), OPTIONAL :: no_outcov  ! no covariance 
!  DOUBLE PRECISION, INTENT(IN), OPTIONAL:: dmeacontr2, dnewton ! radius of second disk, min(minposs)
! -------------------LOCAL ARRAYS (for a single trail)--------------------
  TYPE(tp_point), DIMENSION(lre), INTENT(OUT) :: vas_traceloc     
  INTEGER j,k,jj !,nmoid,nr2
  LOGICAL nooutcov 
!-------------------------------------------------------------------------
  IF(PRESENT(no_outcov))THEN
     nooutcov=no_outcov
  ELSE
     nooutcov=.false.
  ENDIF
!  IF(PRESENT(dmeacontr2).and.PRESENT(dnewton))THEN
!     nmoid=0
!     nr2=0
!     DO j=1,lre
!        IF(vas_trace(j+ire-1)%minposs.lt.dnewton)THEN
!           nmoid=nmoid+1
!           IF(vas_trace(j+ire-1)%b.lt.dmeacontr2)THEN
!              nr2=nr2+1
!           ENDIF
!        ENDIF
!     ENDDO
!  ELSE
! copy all and traslate, removing double VA
     jj=0
     DO j=1,lre
        jj=jj+1
!        IF(jj.gt.1)THEN
!           IF(vas_trace(j+ire-2)%rindex.eq.vas_trace(j+ire-1)%rindex)THEN
!              IF(vas_trace(j+ire-2)%b.gt.vas_trace(j+ire-1)%b)THEN
!                 vas_traceloc(jj-1)=vas_trace(j+ire-1)
!              ENDIF
!              jj=jj-1
!           ELSE
!              vas_traceloc(jj)=vas_trace(j+ire-1)   
!           ENDIF
!        ELSE
           vas_traceloc(jj)=vas_trace(j+ire-1)  
!        ENDIF
     ENDDO
     DO j=1,jj
         CALL wriouttp(vas_traceloc(j),iunout,nooutcov)
     ENDDO
     lre=jj
 ! ENDIF
END SUBROUTINE arrcut

!===========================================================================!
! FIND_INTERVAL                                                             !
!===========================================================================!
! Author: A. Milani, A. Del Vigna                                           !
!===========================================================================!
! This subroutine takes as input two real numbers x < y, and (if all        !
! works) gives as output two consecutive integer numbers intmin and         !
! intmax such that intmin <= x < y <= intmax                                !
!===========================================================================!
SUBROUTINE find_interval(x,y,intmin,intmax) 
  USE output_control
  !=========================================================================================================
  DOUBLE PRECISION, INTENT(IN)  :: x,y            ! Real numbers
  INTEGER,          INTENT(OUT) :: intmin, intmax ! Integer numbers such that intmin <= x < y <= intmax
  !=========================================================================================================
  DOUBLE PRECISION, PARAMETER :: eps=1.d-10
  INTEGER          :: ix,iy           ! Nearest integer for x and y
  DOUBLE PRECISION :: sx, sy          ! x-NINT(x) and the same for y
  INTEGER          :: intminx,intmaxx ! Consecutive integers around x
  INTEGER          :: intminy,intmaxy ! Consecutive integers around y
  !=========================================================================================================
  ix=NINT(x)
  iy=NINT(y)
  sx=x-ix
  sy=y-iy
  IF(ABS(sx).LE.eps)THEN ! x practically an integer
     intminx=ix
     intmaxx=ix
  ELSEIF(sx.GT.0.d0)THEN ! ix < x
     intminx=ix
     intmaxx=ix+1
  ELSEIF(sx.LT.0.d0)THEN ! x < ix
     intminx=ix-1
     intmaxx=ix
  END IF
  IF(ABS(sy).LE.eps)THEN ! y practically an integer
     intminy=iy
     intmaxy=iy
  ELSEIF(sy.GT.0.d0)THEN ! iy < y
     intminy=iy
     intmaxy=iy+1
  ELSEIF(sy.LT.0.d0)THEN ! y < iy
     intminy=iy-1
     intmaxy=iy
  END IF   
  ! Paranoia check: x>y must not occur
  IF(x.GE.y)THEN
     WRITE(ierrou,*) ' find_interval: x ge y ',x,y
     numerr=numerr+1
     STOP '***** find_interval: wrong order *****'
  END IF
  IF(intmaxx.EQ.intmaxy .AND. intminx.EQ.intminy .AND. intmaxx.NE.intminx)THEN
     ! Two real values with same integer part
     intmax=intmaxx
     intmin=intminx
  ELSEIF(intminx.EQ.intmaxx .AND. intminy.EQ.intminx)THEN
     ! x integer and y with the same floor
     intmax=intmaxy
     intmin=intminx
  ELSEIF(intminy.EQ.intmaxy .AND. intmaxy.EQ.intmaxx)THEN
     ! y integer and x with the same ceiling
     intmax=intmaxy
     intmin=intminx
  ELSEIF(intminy.EQ.intmaxy .AND. intminx.EQ.intmaxx .AND. intminy.EQ.intminx+1)THEN
     ! x and y consecutive integers
     intmax=intmaxy
     intmin=intminx
  ELSE
     WRITE(ierrou,*) 'find_interval: not understood ',x,y
     numerr=numerr+1
     intmin=intminx
     intmax=intmaxy
  END IF
END SUBROUTINE find_interval

END MODULE ret_analysistp





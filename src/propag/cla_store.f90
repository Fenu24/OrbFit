!=======================================================================!
! MODULE CLA_STORE                                                      !
!=======================================================================!
! This module contains the routines for the close approach data reading !
! from clo file, and subsequent split into showers and returns. It also !
! contains a routine for the handling of duplications in a return and   !
! the sorting routines needed for the various operations.               !
!=======================================================================!
! Author: G. Tommei, A. Milani                                          !
! Last Update: July 2016 (A. Del Vigna, A. Chessa)                      !
!=======================================================================!
! Other modules used:                                                   !
!    tp_trace                                                           !
!    multi_store                                                        !
!                                                                       !
! Subroutines contained:                                                ! 
!    inclolinctp, showret3tp, arrcut_dupret, split_dupret, sort_time2   !   
!    sort_index2, sort_time2_loc, sort_index2_loc, sort_rindex2_loc,    !
!    showretloc                                                         !
!=======================================================================!
MODULE cla_store
  USE tp_trace
  USE multi_store
  
  IMPLICIT NONE
  
  PUBLIC
  !********************!
  !  PUBLIC VARIABLES  !
  !********************!
  INTEGER                                   :: nox          ! Maximum no. of close app.
  INTEGER                                   :: noxloc       ! Maximum no. of close app. in a return
  INTEGER                                   :: nshx         ! Maximum no. of showers
  INTEGER                                   :: nretx        ! Maximum no. of returns
  TYPE(tp_point),  DIMENSION(:),ALLOCATABLE :: vas_trace    ! Array of Target Plane points
  TYPE(tp_point),  DIMENSION(:),ALLOCATABLE :: vas_trasrt   ! The same, for sorting routines
  TYPE(tp_point),  DIMENSION(:),ALLOCATABLE :: vas_traceloc ! Array of Target Plane points for local storage
  TYPE(tp_point),  DIMENSION(:),ALLOCATABLE :: vas_trlocsrt ! The same, for sorting routines
  INTEGER,         DIMENSION(:),ALLOCATABLE :: isho         ! Index of beginning of the showers
  INTEGER,         DIMENSION(:),ALLOCATABLE :: iret         ! Index of beginning of the returns
  DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE :: tclo         ! Time of closest approach (for each VA)
  DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE :: tclos        ! The same, for sorting routines
  DOUBLE PRECISION                          :: dmeacontr    ! Max dist for cl. app. to be considered (in AU)
  ! Variables for each close app. in return
  DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE :: dist         ! Distance
  DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE :: dminpos      ! Minimum possible distance
  DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE :: stretch      ! Stretching
  DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE :: rindex       ! Real index of the orbit

CONTAINS
  !=============================================================================!
  ! INCLOLINCTP                                                                 !
  !=============================================================================!
  ! - Input close approach data, with confidence ellipse                        !
  ! - Version 2tp, 6 Jan. 2005                                                  !
  !=============================================================================!
  SUBROUTINE inclolinctp(iunout,iunclo,despla,no,no_outcov)
    !=============================================================================================
    INTEGER,           INTENT(IN)           :: iunout    ! Output log unit
    INTEGER,           INTENT(IN)           :: iunclo    ! Input unit for the clo file
    CHARACTER(LEN=15), INTENT(IN)           :: despla    ! Desired planet
    LOGICAL,           INTENT(IN), OPTIONAL :: no_outcov ! No covariance info 
    INTEGER,           INTENT(OUT)          :: no        ! Total number close app 
    !=============================================================================================
    TYPE(tp_point)    :: va_trace ! TP point read from file
    CHARACTER(LEN=15) :: curpla   ! Current planet
    CHARACTER(LEN=15) :: planam   ! Planet name
    INTEGER           :: imultot  ! Counter for records read
    INTEGER           :: imulcur  ! VA index read from clo file
    LOGICAL           :: error    ! Error in reading clo file
    LOGICAL           :: eof      ! End of file
    LOGICAL           :: nooutcov ! The same as no_outcov, input for rea_clorectp
    !=============================================================================================
    curpla=despla 
    no=0 
    imultot=0 
    IF(PRESENT(no_outcov))THEN
       nooutcov=no_outcov
    ELSE
       nooutcov=.FALSE.
    ENDIF

    !###################!
    !  LOOP ON RECORDS  !
    !###################!
1   CONTINUE 
    !**************************************!
    !  Read next record on desired planet  !
    !**************************************!
    CALL rea_clorectp(iunclo,curpla,imulcur,va_trace,planam,error,eof,nooutcov) 
    IF(error)THEN 
       WRITE(*,*) ' inclolinctp error: cloapp ',no,' orbit ', imulcur 
       STOP
    ELSEIF(eof)THEN 
       GOTO 2 
    ELSE 
       ! Discard to avoid complex behaviour
       IF(va_trace%d.GT.dmeacontr) GOTO 1
       no=no+1 
       va_trace%rindex=imulcur
       IF(prob_sampl)THEN
          va_trace%sigma=sigmavalv(imulcur)
       ELSE
          va_trace%sigma=delta_sigma*(va_trace%rindex-imi0)
       ENDIF
       imultot=imultot+1 
       ! Conversion of angles already done in rea_tppoint
       ! Rescaling in Re, find minimum possible
       CALL rescaltp(va_trace)
       IF(no.GT.nox)THEN
          WRITE(*,*) ' too many .clo records, max was ', nox, ' to be increased'
          STOP
       ENDIF
       vas_trace(no)=va_trace
    ENDIF
    GOTO 1 
    !###############!
    !  END OF LOOP  !
    !###############!

    ! End of file exit point
2   WRITE(iunout,*) ' read ',no,' closapp. from ',imultot,' orbits' 
  END SUBROUTINE inclolinctp


  !=============================================================================!
  ! SHOWRET3TP                                                                  !
  !=============================================================================!
  ! - Decomposition of list of close approaches (to the same planet)            !
  !   into showers and returns                                                  !
  ! - Version 4.0: split of returns containing duplicated points                !
  !===================FILE shower_analysis.f====================================!
  ! Created by Giacomo Tommei, 25 November 2002                                 !
  ! TP version 3.2, Milani and Tommei, December 2004                            !
  !=============================================================================!
  SUBROUTINE showret3tp(iunlog,no,dt,tgap,nsho,nret)
    !==============================================================================================
    INTEGER,          INTENT(IN)  :: iunlog ! Log output unit
    INTEGER,          INTENT(IN)  :: no     ! Total number of close approaches
    DOUBLE PRECISION, INTENT(IN)  :: dt     ! Max time span defining a shower
    DOUBLE PRECISION, INTENT(IN)  :: tgap   ! Desired time gap
    INTEGER,          INTENT(OUT) :: nsho   ! Showers number
    INTEGER,          INTENT(OUT) :: nret   ! Returns (trails) number
    !==============================================================================================
    INTEGER          :: lsho          ! Length of shower
    DOUBLE PRECISION :: ts1,ts2       ! Time span of showers
    INTEGER          :: j,js          ! Loop indexes
    INTEGER          :: first_retsho  ! First return of this shower 
    INTEGER          :: nshoret       ! Number of returns in shower 
    LOGICAL          :: dupret(nretx) ! Return with duplicates
    INTEGER          :: lret          ! Current return length
    !*************************************************************!
    !  Variables for the handling of duplicate points in returns  !
    !*************************************************************!
    INTEGER          :: nret_loc      ! Number of sub-returns in the current return
    INTEGER          :: iret_loc(nox) ! List of sub-returns 
    !                                   (iret_loc(1) is the beginning of the current return)
    !==============================================================================================
    !######################!
    !  CUTTING IN SHOWERS  !
    !######################!
    !****************!
    !  Sort by time  !
    !****************!
    tclo(1:no)=vas_trace(1:no)%tcla 
    CALL sort_time2(tclo,no) 
    !********************!
    !  Cut by time span  !
    !********************!
    nsho=1 
    isho(1)=1
    DO j=2,no 
       IF(tclo(j).GT.tclo(isho(nsho))+dt.OR.tclo(j).GT.tclo(j-1)+tgap)THEN
          ! New shower found
          IF(nsho+1.GT.nshx)THEN 
             WRITE(*,*) ' showret3tp: increase nshx, was ',nshx 
             STOP ' ***** showret3tp: FAILED *****' 
          ENDIF
          nsho=nsho+1 
          isho(nsho)=j 
          ! Check for time gap                                                    
          IF(tclo(j).GT.tclo(j-1)+tgap)THEN 
             WRITE(iunlog,178) tclo(j)-tclo(j-1),j,nsho,vas_trace(j)%rindex,tclo(j),vas_trace(j-1)%rindex,tclo(j-1) 
178          FORMAT('showret3tp: shower gap ',1X,F12.5,2(1X,I6),4(1X,F12.5))
          ENDIF
       ENDIF
    ENDDO
    isho(nsho+1)=no+1 
 
    !######################!
    !  CUTTING IN RETURNS  !
    !######################!
    !*******************!
    !  Loop on showers  !
    !*******************!
    ! Initialization of return variables
    nret=0 
    iret(1)=1 
    dupret(1)=.FALSE.
    DO 1 js=1,nsho
       first_retsho=nret        ! First return of this shower
       nshoret=1                ! Number of returns in this shower
       ts1=tclo(isho(js))       ! Beginning time of the shower
       ts2=tclo(isho(js+1)-1)   ! End time of the shower
       lsho=isho(js+1)-isho(js) ! Duration of the shower
       !*************************************************!
       !  Sort by multiple sol. index inside the shower  !
       !*************************************************!
       CALL sort_index2(tclo,no,isho(js),lsho) 
       nret=nret+1            ! First return begins at beginning of shower                         
       iret(nret)=isho(js)    ! Index of the beginning of the first return
       dupret(nret)=.FALSE.   ! Flag for duplicate indexes in the return     
       ! Loop on the current shower
       DO j=isho(js)+1,isho(js+1)-1 
          !**********************************!
          !  Cut by non consecutive indexes  !
          !**********************************!
          IF(vas_trace(j)%rindex.EQ.vas_trace(j-1)%rindex)THEN 
             WRITE(iunlog,188) js, nret, vas_trace(j)%rindex, tclo(j), tclo(j-1)
188          FORMAT('showret3tp: duplicate point, shower, return ',1X,I3,1X,I6,3(1X,F12.5))
             dupret(nret)=.TRUE.
          ELSEIF(vas_trace(j)%rindex.GT.vas_trace(j-1)%rindex+1)THEN 
             lret=j-iret(nret)  ! Return length
             ! New return found, but check if it has duplications
             IF(dupret(nret))THEN
                !********************************************!
                !  Handling of duplicate points in filament  !
                !********************************************!
                CALL arrcut_dupret(iret(nret),lret,nox)                      ! Copy current return in vas_traceloc
                CALL split_dupret(nret,lret,nox,iunlog,nret_loc,iret_loc)    ! Split current return to remove duplications
                vas_trace(iret(nret):iret(nret)+lret-1)=vas_traceloc(1:lret) ! Copy the return after the split into many
                iret(nret:nret+nret_loc-1)=iret_loc(1:nret_loc)              ! Update the list of returns' initial indexes
                dupret(nret:nret+nret_loc-1)=.FALSE.                         ! Update the flag for duplications in a return
                nret=nret+nret_loc       ! It is nret+nret_loc-1 to take into account the sub-ret, then +1 to add the next ret
                iret(nret)=j             ! Beginning of the next return
                nshoret=nshoret+nret_loc ! Count of the sub-returns as returns of this shower
             ELSE
                nret=nret+1          ! Count the finished return
                nshoret=nshoret+1    ! Number of returns in this shower 
                iret(nret)=j         ! Beginning of the following return
                dupret(nret)=.FALSE. ! Inizialization for the following return dupret
             END IF
          END IF
       ENDDO
       !*************************************************!
       !  Last return of a shower can have duplications  !
       !*************************************************!
       IF(j.EQ.isho(js+1))THEN
          lret=j-iret(nret)
          IF(dupret(nret))THEN
             !********************************************!
             !  Handling of duplicate points in filament  !
             !********************************************!
             CALL arrcut_dupret(iret(nret),lret,nox)                      ! Copy current return in vas_traceloc
             CALL split_dupret(nret,lret,nox,iunlog,nret_loc,iret_loc)    ! Split current return to remove duplications
             vas_trace(iret(nret):iret(nret)+lret-1)=vas_traceloc(1:lret) ! Copy the return after the split into many
             iret(nret:nret+nret_loc-1)=iret_loc(1:nret_loc)              ! Update the list of returns' initial indexes
             dupret(nret:nret+nret_loc-1)=.FALSE.                         ! Update the flag for duplications in a return
             nret=nret+nret_loc-1                                         ! The current return was already added
             nshoret=nshoret+nret_loc-1                                   ! Count of the sub-returns as returns of this shower
          END IF
       END IF
       !**************************!
       !  Print in the .out file  !
       !**************************!
       WRITE(iunlog,199) js,nshoret,ts1,ts2 
199    FORMAT(' Shower ',i4,' split in ',i4,' returns, time from ',f13.5 , ' to ', f13.5)                           
1   ENDDO
    !***********************!
    !  End loop on showers  !
    !***********************!
    iret(nret+1)=no+1 
  END SUBROUTINE showret3tp


  !=============================================================================!
  ! ARRCUT_DUPRET                                                               !
  !=============================================================================!
  ! Storage of a filament belonging to vas_trace in vas_traceloc                !
  !=============================================================================!
  SUBROUTINE arrcut_dupret(ire,lre,nox)
    !=======================================================================================
    INTEGER, INTENT(IN)    :: ire ! Address of filament
    INTEGER, INTENT(IN)    :: nox ! Maximum dimension
    INTEGER, INTENT(INOUT) :: lre ! Length of filament
    !=======================================================================================
    IF(lre.GT.nox)THEN
       WRITE(*,*) ' arrcut_dupret: lre =', lre,' nox =', nox
       STOP
    END IF
    vas_traceloc(1:lre)=vas_trace(ire:lre+ire-1)
  END SUBROUTINE arrcut_dupret


  !=============================================================================!
  ! SPLIT_DUPRET                                                                !
  !=============================================================================!
  ! Split a return containing duplications into sub-returns.                    !
  !   (1) The return is sorted by time, then we apply a recursive procedure:    !
  !       we cut the return when a duplicated index is found among the ones     !
  !       from the beginning of the previous cut. This is called a sub-shower.  !
  !   (2) Each sub-shower is sorted by time and is cut into sub-returns         !
  !       following the usual return definition, i.e. a sub-return is a subset  !
  !       of the sub-shower with consecutive indexes.                           !
  !=============================================================================!
  ! Authors: A. Del Vigna, A. Chessa                                            !
  !=============================================================================!
  SUBROUTINE split_dupret(nret,lret,nox,iunlog,nret_loc,iret_loc)
    !====================================================================================================
    INTEGER, INTENT(IN)  :: nret          ! Number of current returns
    INTEGER, INTENT(IN)  :: lret          ! Length of filament
    INTEGER, INTENT(IN)  :: nox           ! Dimension
    INTEGER, INTENT(IN)  :: iunlog        ! Unit for log output
    INTEGER, INTENT(OUT) :: nret_loc      ! Number of sub-returns in the current return
    INTEGER, INTENT(OUT) :: iret_loc(nox) ! List of sub-returns (iret_loc(1) is the beginning the return)
    !====================================================================================================
    DOUBLE PRECISION :: tclo_loc(lret) ! Time of close approach for vas_traceloc
    INTEGER          :: nsho_loc       ! Number of sub-showers
    INTEGER          :: isho_loc(nox)  ! List of sub-showers beginning index
    INTEGER          :: lsho_loc       ! Length of sub-shower
    DOUBLE PRECISION :: ts1            ! Beginning time of the sub-shower
    DOUBLE PRECISION :: ts2            ! End time of the sub-shower
    INTEGER          :: nshoret_loc    ! Number of sub-returns in sub-shower 
    INTEGER          :: j,jj,js        ! Loop indexes
    !====================================================================================================
    ! Initializations of outputs
    nret_loc=0
    iret_loc=0

    !##########################!
    !  CUTTING IN SUB-SHOWERS  !
    !##########################!
    !-------------------------------------------------------------!
    ! The return is sorted by time, then we cut the return when a !
    ! duplicated index is found among the ones from the beginning !
    ! of the previous cut. This is called a sub-shower            !
    !-------------------------------------------------------------!
    !***************************!
    !   Sort tclo_loc by time   !
    !***************************!
    ! Use the routine sort_time2_loc to order by time and to order vas_traceloc at the same times
    tclo_loc(1:lret)=vas_traceloc(1:lret)%tcla 
    CALL sort_time2_loc(tclo_loc,lret)
    !***************************!
    !  Cut by repeated indexes  !
    !***************************!
    nsho_loc=1 
    isho_loc(1)=1
    DO j=2,lret
       DO jj=isho_loc(nsho_loc),j-1
          ! Check if the current index is among the previous ones belonging to the current sub-shower
          IF(vas_traceloc(j)%rindex.EQ.vas_traceloc(jj)%rindex)THEN
             ! New sub-shower beginning at index j
             nsho_loc=nsho_loc+1
             isho_loc(nsho_loc)=j
             ! Stop checking duplications
             EXIT
          END IF
       END DO
    END DO
    isho_loc(nsho_loc+1)=lret+1 
    !***********************!
    !  Loop on sub-showers  !
    !***********************!
    nret_loc=0 
    iret_loc(1)=iret(nret) ! iret_loc has to start from the index of the current return
    DO 1 js=1,nsho_loc
       nshoret_loc=1                        ! Number of sub-returns in this sub-shower
       ts1=tclo_loc(isho_loc(js))           ! Beginning time of the sub-shower
       ts2=tclo_loc(isho_loc(js+1)-1)       ! End time of the sub-shower
       lsho_loc=isho_loc(js+1)-isho_loc(js) ! Time span of the sub-shower
       !*****************************************************!
       !  Sort by multiple sol. index inside the sub-shower  !
       !*****************************************************!
       ! Use the routine sort_index2_loc to order by index and to order vas_traceloc at the same times
       CALL sort_index2_loc(lret,isho_loc(js),lsho_loc) 
       ! First return begins at beginning of sub-shower                         
       nret_loc=nret_loc+1 
       iret_loc(nret_loc)=isho_loc(js)+iret(nret)-1   ! iret_loc has to be shifted by iret(nret)
       ! Cut by non consecutive indexes      
       DO j=isho_loc(js)+1,isho_loc(js+1)-1 
          IF(vas_traceloc(j)%rindex.EQ.vas_traceloc(j-1)%rindex)THEN 
             WRITE(*,*) ' split_dupret: paranoia check, duplication still present!!!'
             STOP
          ELSEIF(vas_traceloc(j)%rindex.GT.vas_traceloc(j-1)%rindex+1)THEN 
             nret_loc=nret_loc+1                ! Add the sub-return just finished
             nshoret_loc=nshoret_loc+1          ! Counter of sub-returns
             iret_loc(nret_loc)=j+iret(nret)-1  ! Beginning of the next return
          END IF
       ENDDO
       !**************************!
       !  Print in the .out file  !
       !**************************!
       WRITE(iunlog,799)js,nshoret_loc,ts1,ts2 
799    FORMAT(' Sub-shower ',i4,' split in ',i4,' sub-returns',               &
            &   ' time from ',f13.5 , ' to ', f13.5)                                                                      
1   ENDDO
    !***********************!
    !  End loop on showers  !
    !***********************!
  END SUBROUTINE split_dupret


  !======================================================================!
  ! SORT_TIME2                                                           !
  !======================================================================!
  ! Sorts a list of close approaches by time                             !
  !======================================================================!
  SUBROUTINE sort_time2(t,no) 
    !============================================================================
    INTEGER,          INTENT(IN)    :: no       ! Actual no of records  
    DOUBLE PRECISION, INTENT(INOUT) :: t(no)    ! Sort key
    !=============END INTERFACE==================================================
    INTEGER                         :: ipo(nox) ! Pointers 
    INTEGER                         :: j        ! Loop indexes 
    !============================================================================
    ! Sort by time, get pointer ipo with the indexes defining the order
    CALL heapsort(t,no,ipo(1:no))
    ! Re-order using the pointer ipo                                            
    DO j=1,no 
       tclos(j)=t(ipo(j)) 
       vas_trasrt(j)=vas_trace(ipo(j))
    ENDDO
    ! Output in input arrays                                                
    DO j=1,no 
       t(j)=tclos(j) 
       vas_trace(j)=vas_trasrt(j) 
    ENDDO
  END SUBROUTINE sort_time2


  !======================================================================!
  ! SORT_INDEX2                                                          !
  !======================================================================!
  ! Sorts a list of close approaches by index of multiple solution       !
  ! but only starting from record number isho                            !
  !======================================================================!
  SUBROUTINE sort_index2(t,no,startsho,lsho) 
    !==============================================================================             
    INTEGER,          INTENT(IN)    :: no       ! Total no of records
    INTEGER,          INTENT(IN)    :: startsho ! Starting row
    INTEGER,          INTENT(IN)    :: lsho     ! Number of rows
    DOUBLE PRECISION, INTENT(INOUT) :: t(no)    ! Time
    !==============================================================================
    INTEGER                         :: i(nox)   ! Sort key 
    INTEGER                         :: ipo(nox) ! Pointers 
    INTEGER                         :: j,jj     ! Loop indexes
    !==============================================================================
    ! Sort by index, get pointer ipo with the indexes defining the order
    DO j=1,no
       i(j)=vas_trace(j)%rindex
    END DO
    CALL heapsorti(i(startsho:startsho+no-1),lsho,ipo)               
    ! Re-order using the pointer ipo
    DO j=1,lsho 
       jj=j+startsho-1 
       tclos(jj)=t(ipo(j)+startsho-1) 
       vas_trasrt(jj)=vas_trace(ipo(j)+startsho-1) 
    ENDDO
    ! Output in input arrays                                                
    DO j=1,lsho 
       jj=j+startsho-1 
       t(jj)=tclos(jj)  
       vas_trace(jj)=vas_trasrt(jj) 
    ENDDO
  END SUBROUTINE sort_index2


  !========================================================================!
  ! SORT_TIME2_LOC                                                         !
  !========================================================================!
  ! Sorts a list of close approaches by time                               !
  !========================================================================!
  SUBROUTINE sort_time2_loc(t,no) 
    !==============================================================================        
    INTEGER,          INTENT(IN)    :: no       ! Actual no of records  
    DOUBLE PRECISION, INTENT(INOUT) :: t(no)    ! Sort key
    !==============================================================================
    INTEGER                         :: ipo(nox) ! Pointer
    INTEGER                         :: j,k      ! Loop indexes 
    !==============================================================================
    ALLOCATE(vas_trlocsrt(no))
    ! Sort by time, get pointer ipo with the indexes defining the order
    CALL heapsort(t,no,ipo(1:no))         
    ! Re-order using the pointer ipo
    DO j=1,no 
       tclos(j)=t(ipo(j)) 
       vas_trlocsrt(j)=vas_traceloc(ipo(j))
    ENDDO
    ! Output in input arrays                                                
    DO j=1,no 
       t(j)=tclos(j) 
       vas_traceloc(j)=vas_trlocsrt(j) 
    ENDDO
    DEALLOCATE(vas_trlocsrt)
  END SUBROUTINE sort_time2_loc


  !======================================================================!
  ! SORT_INDEX2_LOC                                                      !
  !======================================================================!
  ! Sorts a list of close approaches by index of multiple solution       !
  ! but only starting from record number startsho                        !
  !======================================================================!
  SUBROUTINE sort_index2_loc(no,startsho,lsho) 
    !===========================================================================
    INTEGER, INTENT(IN) :: no       ! Total no of records
    INTEGER, INTENT(IN) :: startsho ! Starting row
    INTEGER, INTENT(IN) :: lsho     ! Number of rows
    !===========================================================================
    INTEGER             :: i(nox)   ! Sort key 
    INTEGER             :: ipo(nox) ! Pointers 
    INTEGER             :: j,jj     ! Loop indexes
    !===========================================================================
    ALLOCATE(vas_trlocsrt(no))
    ! Sort by index, get pointer ipo with the indexes defining the order
    DO j=1,no
       i(j)=vas_traceloc(j)%rindex
    END DO
    CALL heapsorti(i(startsho:startsho+no-1),lsho,ipo)               
    ! Re-order using the pointer ipo
    DO j=1,lsho 
       jj=j+startsho-1 
       vas_trlocsrt(jj)=vas_traceloc(ipo(j)+startsho-1) 
    ENDDO
    ! Output in input arrays                                                
    DO j=1,lsho 
       jj=j+startsho-1 
       vas_traceloc(jj)=vas_trlocsrt(jj) 
    ENDDO
    DEALLOCATE(vas_trlocsrt)
  END SUBROUTINE sort_index2_loc


  !======================================================================!
  ! SORT_RINDEX2_LOC                                                     !
  !======================================================================!
  ! Sorts a list of close approaches by real index of multiple solution  !      
  ! but only starting from record number startsho                        !
  !======================================================================!
  SUBROUTINE sort_rindex2_loc(no,startsho,lsho) 
    !=========================================================================          
    INTEGER, INTENT(IN) :: no       ! Total no of records
    INTEGER, INTENT(IN) :: startsho ! Starting row
    INTEGER, INTENT(IN) :: lsho     ! Number of rows
    !=========================================================================
    DOUBLE PRECISION    :: r(nox)   ! Sort key 
    INTEGER             :: ipo(nox) ! Pointers 
    INTEGER             :: j,jj     ! Loop indexes
    !=========================================================================
    ! Sort by index, get pointer ipo with the indexes defining the order
    DO j=1,no
       r(j)=vas_traceloc(j)%rindex
    END DO
    CALL heapsort(r(startsho:startsho+no-1),lsho,ipo)               
    ! Re-order using the pointer ipo
    DO j=1,lsho 
       jj=j+startsho-1 
       vas_trlocsrt(jj)=vas_traceloc(ipo(j)+startsho-1) 
    ENDDO
    ! Output in input arrays                                                
    DO j=1,lsho 
       jj=j+startsho-1 
       vas_traceloc(jj)=vas_trlocsrt(jj) 
    ENDDO
  END SUBROUTINE sort_rindex2_loc


  !================================================================!
  ! SHOWRETLOC                                                     !    
  !================================================================!
  ! Decomposition of list of close approaches (to the same planet) !
  ! into showers and returns                                       !      
  !================================================================!
  SUBROUTINE showretloc(iunlog,no,dt,tgap,dindex,nsho,nret)
    ! ==============================================================================
    INTEGER,          INTENT(IN)  :: iunlog  ! Log output unit
    INTEGER,          INTENT(IN)  :: no      ! Total number of  close approaches
    DOUBLE PRECISION, INTENT(IN)  :: dt      ! Max time span defining a shower
    DOUBLE PRECISION, INTENT(IN)  :: dindex  ! Max gap in LOV index to cut a return
    DOUBLE PRECISION, INTENT(IN)  :: tgap    ! Desired time gap 
    INTEGER,          INTENT(OUT) :: nsho    ! Showers number
    INTEGER,          INTENT(OUT) :: nret    ! Returns (trails) number
    !===============================================================================
    INTEGER                       :: lsho    ! Length of shower
    DOUBLE PRECISION              :: ts1     ! Geninning time of the return
    DOUBLE PRECISION              :: ts2     ! End time of the return
    INTEGER                       :: nshoret ! Number of returns in shower
    INTEGER                       :: j,js    ! Loop indexes 
    !===============================================================================
    !######################!
    !  CUTTING IN SHOWERS  !
    !######################!
    !****************!
    !  Sort by time  !
    !****************!
    tclo(1:no)=vas_traceloc(1:no)%tcla 
    CALL sort_time2_loc(tclo,no) 
    !********************!
    !  Cut by time span  !
    !********************!
    nsho=1 
    isho(1)=1 
    DO j=2,no 
       IF(tclo(j).gt.tclo(isho(nsho))+dt.or.tclo(j).gt.tclo(j-1)+tgap)THEN
          ! New shower found
          IF(nsho+1.gt.nshx)THEN 
             WRITE(*,*)' showretloc: increase nshx, was ',nshx 
             STOP ' ***** showretloc: FAILED *****' 
          ENDIF
          nsho=nsho+1 
          isho(nsho)=j 
          ! Check for time gap
          IF(tclo(j).gt.tclo(j-1)+tgap)THEN 
             WRITE(iunlog,*)'showretloc: shower gap ', tclo(j)-tclo(j-1),j,nsho                            
             WRITE(iunlog,*) vas_traceloc(j)%rindex,tclo(j),vas_traceloc(j-1)%rindex,tclo(j-1) 
          ENDIF
       ENDIF
    ENDDO
    isho(nsho+1)=no+1 

    !######################!
    !  CUTTING IN RETURNS  !
    !######################!
    !*******************!
    !  Loop on showers  !
    !*******************!
    ! Initialization of return variables
    nret=0 
    iret(1)=1 
    DO 1 js=1,nsho 
       nshoret=1                ! Number of returns in this shower
       ts1=tclo(isho(js))       ! Beginning time of the shower
       ts2=tclo(isho(js+1)-1)   ! End time of the shower
       lsho=isho(js+1)-isho(js) ! Duration of the shower
       !*************************************************!
       !  Sort by multiple sol. index inside the shower  !
       !*************************************************!
       CALL sort_index2_loc(no,isho(js),lsho)            
       nret=nret+1            ! First return begins at beginning of shower                         
       iret(nret)=isho(js)    ! Index of the beginning of the first return
       DO j=isho(js)+1,isho(js+1)-1 
          !**********************************!
          !  Cut by non consecutive indexes  !
          !**********************************!
          IF(vas_traceloc(j)%rindex.eq.vas_traceloc(j-1)%rindex)THEN 
             WRITE(iunlog,*)'duplicate point, shower ', js,' return ',nret 
             WRITE(iunlog,*) vas_traceloc(j)%rindex,vas_traceloc(j)%tcla,vas_traceloc(j-1)%tcla 
          ELSEIF(vas_traceloc(j)%rindex.gt.vas_traceloc(j-1)%rindex+dindex)THEN 
             ! New return found
             !***************************************************!
             ! WARNING! Duplication not handled in this version! !
             !***************************************************!
             nret=nret+1 
             nshoret=nshoret+1 
             iret(nret)=j 
          ENDIF
       ENDDO
       !**************************!
       !  Print in the .out file  !
       !**************************!
       WRITE(iunlog,199)js,nshoret,ts1,ts2 
199    FORMAT(' Shower ',i4,' split in ',i4,' returns',               &
            &   ' time from ',f13.5 , ' to ', f13.5)                                                                            
1   ENDDO
    !***********************!
    !  End loop on showers  !
    !***********************!
    iret(nret+1)=no+1 
  END SUBROUTINE showretloc

END MODULE cla_store

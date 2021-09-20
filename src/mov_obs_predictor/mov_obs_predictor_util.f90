!=======================================================!
! MODULE mov_obs_predictor_util                         !
! Utilities for MOV observation predictor main program  !
!=======================================================!
! Author: D. Bracali Cioci                              !
! Last Update: Mar 2020                                 !
!=======================================================!
MODULE mov_obs_predictor_util

  USE fund_const
  USE char_str
  USE output_control
  USE station_coordinates, ONLY: statcode

  IMPLICIT NONE

  PRIVATE

  CHARACTER(LEN=80), PUBLIC :: movdir      ! input directory for MOV sample
  INTEGER,           PUBLIC :: len_movdir  ! length of movdir
  CHARACTER(LEN=50), PUBLIC :: cal_date    ! prediction time (CAL)
  REAL(KIND=dkind),  PUBLIC :: t_pred      ! prediction time (MJD)
  REAL(KIND=dkind),  PUBLIC :: op_sigma    ! obs. prediction maximum sigma
  INTEGER,           PUBLIC :: id_obscode  ! MPC observatory code

  PUBLIC :: mov_pred_opt

CONTAINS

  ! ======================================================!
  ! MODULE mov_pred_opt                                   !
  ! Options ingestion for main mov_obs_predictor          !
  ! ======================================================!
  ! author: D.Bracali Cioci                               !
  ! last update: Mar 2020                                 !
  ! ======================================================!

  SUBROUTINE mov_pred_opt(progna,run)

    CHARACTER(LEN=6), INTENT(IN) :: progna   ! program name
    CHARACTER(LEN=80),INTENT(IN) :: run      ! run name
    ! Local variables
    CHARACTER(LEN=100)   :: comment
    CHARACTER(LEN=3)     :: obs_code
    LOGICAL              :: ireq,found

    ! Read options from file <run name>.opt
    CALL initopt(progna,run,'opt')

    ! Read option for physical model and integration method
    CALL rmodel(1)

    WRITE(iun_log,*) 'Reading options:'
    ! ====================================================================
    ! read options for this program from stored list of character arrays
    ireq    = .TRUE.
    comment = 'asteroid name '
    CALL input_cha_opt(progna,'astna',name_obj,ireq,found,comment,iun_log)
    CALL rmsp(name_obj,lobjnam)
    ! ====================================================================
    ireq    = .FALSE.
    movdir  = './movsample'
    comment = 'input directory for MOV sample '
    CALL input_cha_opt(progna,'movdir',movdir,ireq,found,comment,iun_log)
    CALL rmsp(movdir,len_movdir)
    ! ====================================================================
    ireq    = .TRUE.
    comment = 'Calendar date of prediction '
    CALL input_cha_opt(progna,'cal_date',cal_date,ireq,found,comment,iun_log)
    CALL caldate_to_time(cal_date,t_pred)
    ! ====================================================================
    ireq    = .TRUE.
    comment = 'MPC Observatory code'
    CALL input_cha_opt(progna,'obs_code',obs_code,ireq,found,comment,iun_log)
    CALL statcode(obs_code,id_obscode)
    ! ====================================================================
    ireq     = .FALSE.
    op_sigma = 1.d0
    comment  = 'obs. prediction maximum sigma'
    CALL input_rea_opt(progna,'op_sigma',op_sigma,ireq,found,comment,iun_log)

  END SUBROUTINE mov_pred_opt

  !==============================================================!
  ! SUBROUTINE caldate_to_time                                   !
  ! Conversion from a calendar date string to the corresponding  !
  ! type time in UTC scale                                       !
  !==============================================================!
  ! Author: D. Bracali Cioci                                     !
  ! Last Update: Mar 2020                                        !
  !==============================================================!
  SUBROUTINE caldate_to_time(calend_cha,time_out)

    CHARACTER(LEN=50), INTENT(IN)  :: calend_cha ! calendar date string
    REAL(KIND=dkind),  INTENT(OUT) :: time_out   ! corresponding time
    ! Local variables
    REAL(KIND=dkind) :: tjm1
    REAL(KIND=dkind) :: hour                           ! hours + fractions of hours
    INTEGER          :: year,month,day_int             ! year, month and day
    INTEGER          :: hour_int,minute_int,second_int ! hours, minutes and seconds
    LOGICAL          :: isnum
    EXTERNAL         :: tjm1,isnum

    ! Find the first / separator in date
    IF(calend_cha(5:5).NE.'/')THEN
       WRITE(*,*) 'caldate_to_time'//&
            & ' missing the first / separator in the input calendar date'
       WRITE(ierrou,*) 'caldate_to_time'//&
            & ' missing the first / separator in the input calendar date'
       STOP
    END IF
    ! Reading year of date
    IF(isnum(calend_cha(1:4)))THEN
       READ(calend_cha(1:4),'(I4)') year
    ELSE
       WRITE(*,*) 'caldate_to_time'//&
            & ' field "year" is not a number in the input calendar date'
       WRITE(ierrou,*) 'caldate_to_time'//&
            & ' field "year" is not a number in the input calendar date'
       STOP
    END IF
    ! Find the second / separator in date
    IF(calend_cha(8:8).NE.'/')THEN
       WRITE(*,*) 'caldate_to_time'//&
            & ' missing the second / separator in the input calendar date'
       WRITE(ierrou,*) 'caldate_to_time'//&
            & ' missing the second / separator in the input calendar date'
       STOP
    END IF
    ! Reading month of date
    IF(isnum(calend_cha(6:7)))THEN
       READ(calend_cha(6:7),'(I2)') month
    ELSE
       WRITE(*,*) 'caldate_to_time'//&
            & ' field "month" is not a number in the input calendar date'
       WRITE(ierrou,*) 'caldate_to_time'//&
            & ' field "month" is not a number in the input calendar date'
       STOP
    END IF
    ! Find comma separator in date
    IF(calend_cha(11:11).NE.',')THEN
       WRITE(*,*) 'caldate_to_time'//&
            & ' missing comma separator in the input calendar date'
       WRITE(ierrou,*) 'caldate_to_time'//&
            & ' missing comma separator in the input calendar date'
       STOP
    END IF
    ! Reading day of date
    IF(isnum(calend_cha(9:10)))THEN
       READ(calend_cha(9:10),'(I2)') day_int
    ELSE
       WRITE(*,*) 'caldate_to_time'//&
            & ' field "day" is not a number in the input calendar date'
       WRITE(ierrou,*) 'caldate_to_time'//&
            & ' field "day" is not a number in the input calendar date'
       STOP
    END IF
    ! Find first colon separator
    IF(calend_cha(14:14).NE.':')THEN
       WRITE(*,*) 'caldate_to_time'//&
            & ' missing first colon separator in the input calendar date'
       WRITE(ierrou,*) 'caldate_to_time'//&
            & ' missing first colon separator in the input calendar date'
       STOP
    END IF
    ! Reading hours of date
    IF(isnum(calend_cha(12:13)))THEN
       READ(calend_cha(12:13),'(I2)') hour_int
    ELSE
       WRITE(*,*) 'caldate_to_time'//&
            & 'field "hours" is not a number in the input calendar date'
       WRITE(ierrou,*) 'caldate_to_time'//&
            & 'field "hours" is not a number in the input calendar date'
       STOP
    END IF
    ! Find second colon separator
    IF(calend_cha(17:17).NE.':')THEN
       WRITE(*,*) 'caldate_to_time'//&
            & 'missing second colon separator in the input calendar date'
       WRITE(ierrou,*) 'caldate_to_time'//&
            & 'missing second colon separator in the input calendar date'
       STOP
    END IF
    ! Reading minutes of date
    IF(isnum(calend_cha(15:16)))THEN
       READ(calend_cha(15:16),'(I2)') minute_int
    ELSE
       WRITE(*,*) 'caldate_to_time'//&
            & 'field "minutes" is not a number in the input calendar date'
       WRITE(ierrou,*) 'caldate_to_time'//&
            & 'field "minutes" is not a number in the input calendar date'
       STOP
    END IF
    ! Reading seconds of date
    IF(isnum(calend_cha(18:19)))THEN
       READ(calend_cha(18:19),'(I2)') second_int
    ELSE
       WRITE(*,*) 'caldate_to_time'//&
            & 'field "seconds" is not a number in the input calendar date'
       WRITE(ierrou,*) 'caldate_to_time'//&
            & 'field "seconds" is not a number in the input calendar date'
       STOP
    END IF
    ! Final conversion into type time
    hour  = hour_int+(minute_int+second_int/60.d0)/60.d0
    time_out = tjm1(day_int,month,year,hour)

  END SUBROUTINE caldate_to_time

END MODULE mov_obs_predictor_util

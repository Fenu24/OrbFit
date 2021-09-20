! Copyright (C) 1999 by Mario Carpino (carpino@brera.mi.astro.it)
! Version: February 8, 1999
! ---------------------------------------------------------------------
!
!  *****************************************************************
!  *                                                               *
!  *                         O F C L R F                           *
!  *                                                               *
!  *                      Close report files                       *
!  *                                                               *
!  *****************************************************************
!
SUBROUTINE ofclrf
  USE output_control
  IMPLICIT NONE
  
  CHARACTER(LEN=12):: status

! Close report file for errors
  IF(numerr.LE.0) THEN
     status='DELETE'
  ELSE
     status='KEEP'
  END IF
  CALL filclo(ierrou,status)
  
! Close report file for close approaches
  IF(numcla.LE.0) THEN
     status='DELETE'
  ELSE
     status='KEEP'
  END IF
  CALL filclo(iuncla,status)

! Close report file for propagator parameters
  CALL filclo(ipirip,' ')
  
END SUBROUTINE ofclrf

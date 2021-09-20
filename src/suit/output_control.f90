MODULE output_control

! output control 
!   propagator param. file unit, 
integer ipirip
!   close app. file unit, number close app, error file unit, number of errors
!   convention: if units are negative, the close app/errors are output
!   to standard output and the errors are not counted
integer iuncla,numcla,ierrou,numerr,iwarou,numwar
! asteroid name for messages
CHARACTER*30 name_obj
INTEGER lobjnam
! unit for log report (replace the old iun20, iun8, iuncovar, iunlog, iundummy, and so on)
INTEGER iun_log,iun_covar
! =========verbosity control==========
! for close approaches
INTEGER verb_clo
! for propagator
INTEGER verb_pro
! for observartions
INTEGER verb_obs
! for differential corrections
INTEGER verb_dif
! for multiple solutions
INTEGER verb_mul
! for outlier rejection
INTEGER verb_rej
! for input/output
INTEGER verb_io
! for moid computation
INTEGER verb_moid
! for preliminary orbits
INTEGER verb_prelim
! for output of covariance
INTEGER verb_covariance
! for matrix inversion, conditioning
INTEGER verb_matrix
! for keplerien integrals
INTEGER verb_2body
! general rules: verb_x=10 in fitobs; verb_x=1 (default) in batch
! programs running many orbits; verb_x=5 to have more error messages
! on the consolle; verb_x>20 to produce huge output, for debugging only






END MODULE output_control

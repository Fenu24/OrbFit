MODULE multi_store

IMPLICIT NONE

PUBLIC
! to avoid duplications between multiple_sol and tp_trace (or circular dependency)
LOGICAL :: prob_sampl ! if true, sampling with step constatnt in probability

! deltasigma, max value of sigma - constant delta_sigma_lov case
DOUBLE PRECISION :: delta_sigma, sigma_max

! index of nominal solution 
INTEGER :: imi0, imul0

! contant IP steps
INTEGER, PARAMETER :: mulx = 99999
DOUBLE PRECISION   :: sigmavalv(mulx)
! end duplications

END MODULE multi_store

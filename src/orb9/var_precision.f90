! Precision for real variables
! Also long integers
MODULE var_precision
IMPLICIT NONE
INTEGER, PARAMETER        :: r_kind = SELECTED_REAL_KIND(15,50)
INTEGER, PARAMETER        :: long_int = SELECTED_INT_KIND(11)
END MODULE var_precision

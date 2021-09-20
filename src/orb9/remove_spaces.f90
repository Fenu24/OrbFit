! = = = = = = = = = = = = = = = = = = = = = = = = = = =
! Remove all the spaces from a character string

SUBROUTINE remove_spaces(string,length)
IMPLICIT NONE

INTEGER, PARAMETER :: n_spaces = 3 ! no. of chars considered as spaces
INTEGER, PARAMETER, DIMENSION(n_spaces) :: space_iachar = (/32,9,44/)! ASCII code of space chars (blank,TAB,comma)

CHARACTER(LEN=*),INTENT(INOUT) :: string ! character string
INTEGER, INTENT(OUT), OPTIONAL :: length ! string length on output
CHARACTER(LEN=LEN(string)) :: c1
INTEGER :: i,icc,k,l1,ll

ll = LEN(string)
l1 = 0
c1 = ' '
inp_char: DO i=1,ll
    icc = IACHAR(string(i:i))
    spaces: DO k=1,n_spaces
        IF(icc==space_iachar(k)) CYCLE inp_char
    END DO spaces
    l1 = l1 + 1
    c1(l1:l1) = string(i:i)
END DO inp_char
string = c1
IF(PRESENT(length)) length = l1

END SUBROUTINE remove_spaces

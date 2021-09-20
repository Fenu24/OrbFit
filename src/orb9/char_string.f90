! Elementary manipulation of character strings

! WARNING: throughout all the module, a "space" is any of the characters defined by PARAMETER space_iachar, a "quote" is any
!          of the characters defined by PARAMETER quote_iachar
!
! List of procedures:
!
! SUBROUTINE remove_spaces:		remove all the spaces from a character string
!
! FUNCTION lench:			length of a character string excluding trailing "spaces": it is similar to the
!					intrinsic LEN_TRIM, but the definition of "space" is different (see PARAMETER space_iachar)
!
! FUNCTION normal_string:		produce the "normal form" of a character string (leading spaces are removed, any succession
!					of one or more space characters is substituted with a single blank but "quoted" strings
!					are preserved literally)
!
! SUBROUTINE get_string_content:	get the "content" of a string (enclosed by quotes or not)
!
! SUBROUTINE cut_string:		cuts a character string into pieces, using the separator specified by the caller

MODULE char_string
IMPLICIT NONE
PRIVATE

INTEGER, PARAMETER :: max_char_len = 13 ! max length of character strings

INTEGER, PARAMETER :: n_spaces = 3 ! no. of chars considered as spaces
INTEGER, PARAMETER, DIMENSION(n_spaces) :: space_iachar = (/32,9,44/)! ASCII code of space chars (blank,TAB,comma)

INTEGER, PARAMETER :: n_quotes = 2 ! no. of chars considered as quotes
INTEGER, PARAMETER, DIMENSION(n_quotes) :: quote_iachar = (/39,34/)! ASCII code of quote chars (',")

! PUBLIC ENTITIES

PUBLIC :: max_char_len ! PARAMETERs
PUBLIC :: lench,normal_string ! FUNCTIONs
PUBLIC :: remove_spaces,get_string_content,cut_string ! SUBROUTINEs

CONTAINS

! = = = = = = = = = = = = = = = = = = = = = = = = = = =
! Remove all the spaces from a character string

SUBROUTINE remove_spaces(string,length)
IMPLICIT NONE
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

! = = = = = = = = = = = = = = = = = = = = = = = = = = =
! Length of a character string excluding trailing spaces

FUNCTION lench(string)
IMPLICIT NONE

CHARACTER(LEN=*),INTENT(IN):: string ! character string
INTEGER :: lench    ! string length

INTEGER :: i,icc,k

char: DO i=LEN(string),1,-1
    icc = IACHAR(string(i:i))
    spaces: DO k=1,n_spaces
        IF(icc==space_iachar(k)) CYCLE char
    END DO spaces
    lench = i
    RETURN
END DO char
lench = 0

END FUNCTION lench

! = = = = = = = = = = = = = = = = = = = = = = 
! Produce the "normal form" of a character string:
!	1) leading spaces are removed
!	2) any succession of one or more space characters is substituted with a single blank
!	3) "quoted" strings are preserved literally (spaces are not changed)

FUNCTION normal_string(string)
IMPLICIT NONE

CHARACTER(LEN=*),INTENT(IN):: string! character string
CHARACTER(LEN=LEN(string)):: normal_string
INTEGER :: i,icc,iq,k,l
LOGICAL :: prev_space,quoted

normal_string = ' '

l = 0
prev_space = .true.
quoted = .false.
ch: DO i=1,LEN(string)
    icc = IACHAR(string(i:i))
    IF(quoted) THEN
        l = l + 1
        normal_string(l:l) = string(i:i)
        IF(icc==iq) THEN
            quoted = .false.
            prev_space = .false.
        END IF
    ELSE
        quote: DO k=1,n_quotes
            IF(icc==quote_iachar(k)) THEN
                quoted = .true.
                iq = icc
                l = l + 1
                normal_string(l:l) = string(i:i)
                CYCLE ch
            END IF
        END DO quote
        spaces: DO k=1,n_spaces
            IF(icc==space_iachar(k)) THEN
                IF(.NOT.prev_space) THEN
                    l = l + 1
                    normal_string(l:l) = ' '
                    prev_space = .true.
                END IF
                CYCLE ch
            END IF
        END DO spaces
        l = l + 1
        normal_string(l:l) = string(i:i)
        prev_space = .false.
    END IF
END DO ch

END FUNCTION normal_string

! = = = = = = = = = = = = = = = = = = = = = = = = = = 
! Get the content of a quoted string, namely the characters between enclosing quotes (' or ", as defined by PARAMETER quote_iachar)
! or, if the first non-space characters of the string are not quotes, 
! For instance, an input of the form: >  "abc"  "cdf" <
!    gives on output:  >abc< (as CONTENT) and >  "cdf" < (as REST)
! an input of the form >  ABC  sdf<
!    gives on output:  >ABC< (as CONTENT) and >  sdf< (as REST)

SUBROUTINE get_string_content(string,content,rest,error)
IMPLICIT NONE
CHARACTER(LEN=*),INTENT(IN):: string! quoted string
CHARACTER(LEN=*),INTENT(OUT):: content ! content between quotes
CHARACTER(LEN=*),INTENT(OUT), OPTIONAL:: rest! rest of the input string (after closing quote)
CHARACTER(LEN=*),INTENT(OUT), OPTIONAL:: error! parsing error message
INTEGER :: i,i1,i2,icc,k,ls
CHARACTER(LEN=1):: qch
LOGICAL :: quoted

! Determine the beginning of the CONTENT string, by discarding leading spaces and finding the position of the opening quotes
! I1 is the position of the first character of CONTENT; QCH is the quote character found
ls = LEN(string)
open_quote: DO i=1,ls
    icc = IACHAR(string(i:i))
    spaces: DO k=1,n_spaces ! Skip leading space characters
        IF(icc==space_iachar(k)) CYCLE open_quote
    END DO spaces
    i1 = i + 1
    quote: DO k=1,n_quotes ! Look if the first non-space character is a quote
        IF(icc==quote_iachar(k)) THEN
            qch = string(i:i)
            quoted = .true.
            i1 = i + 1
            GOTO 1
        END IF
    END DO quote
! The first non-space character is not a quote
    quoted = .false.
    i1 = i
    GOTO 1
END DO open_quote

! The input string contains only space characters
content = ' '
IF(PRESENT(rest)) rest = ' '
RETURN

! Determine the end of the CONTENT string, by looking for a second occurrence of the quote character (QCH)
! I2 is the position of the last character of CONTENT
1 CONTINUE

IF(quoted) THEN
    close_quote: DO i=i1,ls
        IF(string(i:i)==qch) THEN
            i2 = i - 1
            GOTO 2
        END IF
    END DO close_quote
    IF(PRESENT(error)) THEN
        error = 'Missing closing quotes'
        content = ' '
        IF(PRESENT(rest)) rest = ' '
        RETURN
    END IF
    STOP ' **** char_string.get_string_content: missing closing quotes ****'
ELSE
    end_of_string: DO i=i1,ls
        icc = IACHAR(string(i:i))
        end_space: DO k=1,n_spaces
            IF(icc==space_iachar(k)) THEN
                i2 = i - 1
                GOTO 2
            END IF
        END DO end_space
    END DO end_of_string
    i2 = ls
END IF

! Extract CONTENT and REST
2 CONTINUE
IF(i2<i1) THEN
    content = ' '
ELSE
    IF(LEN(content)<i2-i1+1) STOP ' **** char_string.get_string_content: insufficient length for CONTENT ****'
    content = string(i1:i2)
END IF
IF(PRESENT(error)) error = ' '
IF(PRESENT(rest)) THEN
    IF(i2>ls-2) THEN
        rest = ' '
    ELSE
        IF(LEN(rest)<ls-i2-1) STOP ' **** char_string.get_string_content: insufficient length for REST ****'
        rest = string(i2+2:)
    END IF
END IF

END SUBROUTINE get_string_content

! = = = = = = = = = = = = = = = = = = = = = = = = = = = = 
! Cuts a character string into pieces, using the specified separator
! If the separator does not appear in the string, all the string is returned as "first_piece" and the "rest" is empty

SUBROUTINE cut_string(string,separator,first_piece,rest)
IMPLICIT NONE

CHARACTER(LEN=*),INTENT(IN):: string! input string to be cut
CHARACTER(LEN=*),INTENT(IN):: separator! separator
CHARACTER(LEN=*),INTENT(OUT):: first_piece! first piece of the string (preceding the separator)
CHARACTER(LEN=*),INTENT(OUT), OPTIONAL:: rest! rest of the string (following the separator)

INTEGER ::ipos,lensep

ipos = INDEX(string,separator)
IF(ipos==0) THEN
    IF(LEN(first_piece)<LEN(string)) STOP ' **** char_string.cut_string: insufficient length for FIRST_PIECE ****'
    first_piece = string
    IF(PRESENT(rest)) rest = ' '
    RETURN
END IF
lensep = LEN(separator)
IF(LEN(first_piece)<ipos-1) STOP ' **** char_string.cut_string: insufficient length for FIRST_PIECE ****'
first_piece = string(:ipos-1)
IF(PRESENT(rest)) THEN
    IF(LEN(rest)<LEN(string)-ipos-lensep+1) STOP ' **** char_string.cut_string: insufficient length for REST ****'
    rest = string(ipos+lensep:)
END IF

END SUBROUTINE cut_string

END MODULE char_string

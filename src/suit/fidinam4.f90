
! =======================================                               
! FILNAM                                                                
!                                                                       
! computes  file name in assigned directory                             
! ========================================                              
SUBROUTINE filnam(eledir,astnam,suffix,file,le) 
  IMPLICIT NONE 
  CHARACTER*(*) eledir 
  CHARACTER*(*) astnam 
  CHARACTER*(*) suffix 
  CHARACTER*(*) file 
  INTEGER le 
!==== end interface                                                     
  INTEGER ldir,l1,l2,l3 
  PARAMETER (ldir=60) 
  CHARACTER*200 tmpstr 
  INTEGER lench 
  INCLUDE 'sysdep.h90' 
  l1=lench(eledir) 
  l2=lench(astnam) 
  l3=lench(suffix) 
  IF(l1+l2+l3+2.gt.200)THEN 
     WRITE(*,*) 'filnam: Filename too long:',                       &
          &        eledir(1:l1),dircha,astnam(1:l2),'.',suffix(1:l3)         
     STOP 
  ENDIF
  tmpstr=eledir(1:l1)//dircha//astnam(1:l2)//'.'//suffix(1:l3) 
  CALL rmsp(tmpstr,le) 
  file=tmpstr(1:le) 
  IF(le.gt.ldir-4)THEN 
     WRITE(*,*)'filnam: possible file name truncation',astnam,file 
  ENDIF
END SUBROUTINE filnam
                                                                        
! ==============================================                        
! handling of name with multiple solutions                              
SUBROUTINE splinam(name0,name1,m) 
  IMPLICIT NONE 
  CHARACTER*(*) name0 
  CHARACTER*9 name1 
  INTEGER m 
  INTEGER i 
! find separator                                                        
  i=index(name0,'_') 
  IF(i.le.0)THEN 
!     WRITE(*,*)' name splitting problem for ', name0                
     name1=name0(1:9) 
     m=1 
     RETURN 
  ENDIF
  name1=' ' 
  name1=name0(1:i-1) 
  READ(name0(i+1:9),*)m 
  !     WRITE(*,*)name1,m                                                 
END SUBROUTINE splinam
! =======================================                               
! FIDINAM vers. 4
!                                                                       
! computes  file name in assigned directory, including ONS  
! ========================================                              
SUBROUTINE fidinam(eledir0,astnam0,suffix,file,le) 
  USE name_rules
  IMPLICIT NONE 
  INTEGER, PARAMETER :: ldir=100 
  CHARACTER*(60), INTENT(IN) :: eledir0 
  CHARACTER*(name_len), INTENT(IN) :: astnam0 
  CHARACTER*3, INTENT(IN) ::  suffix 
  CHARACTER*(*), INTENT(OUT) ::  file
  INTEGER, INTENT(OUT) :: le 
! END INTERFACE
  CHARACTER(name_len) nam0, astnam
  CHARACTER*60 fulldir
  CHARACTER*200 eledir
  INTEGER ld,ll,ld1 
  INCLUDE 'sysdep.h90' 
  INTEGER j 
  CHARACTER*10 dircom4,direct 
  eledir=eledir0
  CALL rmsp(eledir,le) 
! convert name in case it is of the form nam0=namp
  astnam=astnam0                      
  call rmsp(astnam,ld1) 
  ll=index(astnam,'=') 
  IF(ll.eq.0)THEN 
     nam0=astnam(1:9) 
  ELSE 
     nam0=astnam(1:ll-1) 
  ENDIF
  CALL rmsp(nam0,ld1) 
  direct=dircom4(nam0,ld)
  IF(eledir(le:le).eq.dircha)THEN 
     fulldir=eledir(1:le)//direct(1:ld) 
  ELSE 
     fulldir=eledir(1:le)//dircha//direct(1:ld) 
  ENDIF
  CALL rmsp(fulldir,le)
  IF(le.gt.ldir)THEN 
     WRITE(*,*)'fidinam: possible directory name truncation for ',astnam
     WRITE(*,*) fulldir 
  ENDIF
  file=fulldir(1:le)//dircha//astnam 
  CALL rmsp(file,le) 
  file=file(1:le)//'.'//suffix 
  CALL rmsp(file,le)
  IF(le.eq.60)THEN
    WRITE(*,*)'fidinam: possible file name truncation for ',astnam
     WRITE(*,*) file
  ENDIF
  DO j=le+1,60 
     file(j:j)=' ' 
  ENDDO
  RETURN 
END SUBROUTINE fidinam
! =============================                                         
! DIRCOM4                                                                
!                                                                       
! computes subdirectory name                                            
CHARACTER*10 FUNCTION dircom4(name0,le) 
  USE output_control 
  IMPLICIT NONE 
  CHARACTER*9 name0 
  INTEGER le 
  CHARACTER*9 name 
  CHARACTER*1 nm 
  INTEGER number,ll,ld,nyea,nodir
! new ONS part
  CHARACTER*3 dat
  DOUBLE PRECISION mjd
  INTEGER lun,lunation
! make sure it is left aligned                                          
  name=name0 
  CALL rmsp(name,ll) 
! find if numbered                                                      
  READ(name,FMT='(I9)',ERR=10)number 
! numbered: directories of 1000 each                                    
  nodir=number/1000 
  dircom4=' ' 
  IF(nodir.gt.99)THEN
     WRITE(dircom4,'(I3)')nodir 
     le=3
  ELSEIF(nodir.gt.9)THEN 
     WRITE(dircom4,'(I2)')nodir 
     le=2 
  ELSE 
     WRITE(dircom4,'(I1)')nodir 
     le=1 
  ENDIF
  RETURN 
! unnumbered                                                            
10 CONTINUE 
! select surveys                                                        
  ld=index(name,'-') 
  IF(ld.gt.0)THEN 
! surveys                                                               
     IF(index(name,'P-L').ne.0)THEN 
! Palomar-Leiden                                                        
        dircom4='P-L' 
        le=3 
     ELSEIF(index(name,'T-').ne.0)THEN 
! trojan surveys                                                        
        le=3 
        IF(index(name,'T-1').ne.0)THEN 
           dircom4='T-1' 
        ELSEIF(index(name,'T-2').ne.0)THEN 
           dircom4='T-2' 
        ELSEIF(index(name,'T-3').ne.0)THEN 
           dircom4='T-3' 
        ELSE 
           WRITE(*,*)'dircom4: name parsing error',name0 
           WRITE(ierrou,*)'dircom4: name parsing error',name0 
           numerr=numerr+1 
           dircom4='unknown' 
           le=7 
        ENDIF
     ELSE 
        WRITE(*,*)'dircom4: name parsing error',name0 
        WRITE(ierrou,*)'dircom4: name parsing error',name0 
        numerr=numerr+1 
        dircom4='unknown' 
        le=7 
     ENDIF
  ELSE 
! by year                                                               
     READ(name(1:4),FMT='(I4)',ERR=20) nyea 
     IF(nyea.lt.1920)THEN 
        dircom4='old' 
        le=3 
     ELSEIF(nyea.lt.1930)THEN 
        dircom4='1920' 
        le=4 
     ELSEIF(nyea.lt.1940)THEN 
        dircom4='1930' 
        le=4 
     ELSEIF(nyea.lt.1950)THEN 
        dircom4='1940' 
        le=4 
     ELSEIF(nyea.lt.1960)THEN 
        dircom4='1950' 
        le=4 
     ELSEIF(nyea.lt.1970)THEN 
        dircom4='1960' 
        le=4 
     ELSEIF(nyea.lt.1980)THEN 
        dircom4='1970' 
        le=4 
     ELSEIF(nyea.ge.1998)THEN 
        nm=name(5:5) 
        WRITE(dircom4,199)nyea,nm 
199     FORMAT(I4,A1) 
        le=5 
     ELSEIF(nyea.lt.0.or.nyea.gt.2100)THEN 
        WRITE(*,*)'dircom4: name parsing error',name0 
        WRITE(ierrou,*)'dircom4: name parsing error',name0 
        numerr=numerr+1 
        dircom4='unknown' 
        le=7 
     ELSE 
        WRITE(dircom4,FMT='(I4)') nyea 
        le=4 
     ENDIF
  ENDIF
  RETURN 
20 CONTINUE
! ONS?
  READ(name(5:7),'(A3)')dat
  CALL decode_mjd(dat,mjd)
  lun=lunation(mjd)
  IF(lun.lt.0.or.lun.gt.999)GOTO 40
  WRITE(dircom4,299)lun
299 FORMAT('ons_',I3)
  CALL rmsp(dircom4,le)
  RETURN
! error parsing year                                                    
40 CONTINUE 
  WRITE(*,*)'dircom4: name parsing error',name0 
  WRITE(ierrou,*)'dircom4: name parsing error',name0 
  numerr=numerr+1 
  dircom4='unknown' 
  le=7 
  RETURN 
END FUNCTION dircom4

INTEGER FUNCTION lunation(mjd) ! number of lunations, from some origin
DOUBLE PRECISION, INTENT(IN) :: mjd ! MJD date
DOUBLE PRECISION, PARAMETER :: per_syn=29.53059d0
DOUBLE PRECISION, PARAMETER :: new_moon0=5.2717d4 ! 19 march 2003
!
DOUBLE PRECISION cycles

cycles=(mjd-new_moon0)/per_syn
lunation=floor(cycles)+200
END FUNCTION lunation

SUBROUTINE lun_limits(lun,t1,t2)
IMPLICIT NONE
INTEGER, INTENT(IN) :: lun
DOUBLE PRECISION, INTENT(OUT) :: t1,t2
DOUBLE PRECISION, PARAMETER :: per_syn=29.53059d0
DOUBLE PRECISION, PARAMETER :: new_moon0=5.2717d4 ! 19 march 2003
t1=new_moon0+(lun-200)*per_syn
t2=t1+per_syn
END SUBROUTINE lun_limits

!!!!!!!!!!!!!!!!!!!!!
! Encoding routines !
!!!!!!!!!!!!!!!!!!!!!
!
! ons_code2 encodes designation
SUBROUTINE ons_code2(n, y, m , d, station, our_ons_code)
IMPLICIT NONE

INTEGER, INTENT(IN) :: n, y, m, d
CHARACTER (LEN=3), INTENT(IN)  :: station
CHARACTER (LEN=9), INTENT(OUT) :: our_ons_code

CHARACTER (LEN=4) :: string_number
CHARACTER (LEN=3) :: string_date
CHARACTER (LEN=2) :: string_station

CALL encode_number(n, string_number)
CALL encode_date(y,m,d, string_date)
CALL encode_station(station, string_station)

our_ons_code =  string_number// string_date // string_station

END SUBROUTINE ons_code2


!=========================================================
!              encode_number
!=========================================================
! INPUT  INTEGER :: 0< n < 54*64**3  
! OUTPUT CHARACTER (LEN=4) :: n_enc 
!
!=========================================================


SUBROUTINE encode_number(n, n_enc)
IMPLICIT NONE 

INTEGER, INTENT(IN) :: n
CHARACTER (LEN=4), INTENT(OUT) :: n_enc

CALL base_10_to_64(n+10*64**3,n_enc)

END SUBROUTINE encode_number

!=========================================================
!              encode_date
!=========================================================
! INPUT   INTEGER :: y, m, d
! OUTPUT  CHARACTER (LEN=3) :: date_code 
! 
!=========================================================

SUBROUTINE encode_date(y, m, d, date_code)
IMPLICIT NONE 

INTEGER, INTENT(IN) :: y, m, d
CHARACTER (LEN=3) :: date_code

DOUBLE PRECISION :: jd
CHARACTER (LEN=4) :: temp_date_code

CALL julian(y,m,d,0,0,0.d0, jd)
!WRITE(*,*) "julian day: ", jd-2400000.5   
CALL base_10_to_64(INT(jd-2400000.5), temp_date_code)

IF(LEN(TRIM(temp_date_code))==1) THEN
date_code = "00"//temp_date_code(1:1)
ELSEIF(LEN(TRIM(temp_date_code))==2) THEN
date_code = "0"//temp_date_code(1:2)
ELSEIF(LEN(TRIM(temp_date_code))==3) THEN
date_code = temp_date_code(1:3)
END IF

END SUBROUTINE encode_date


!=========================================================
!              encode_station
!=========================================================
! INPUT  CHARACTER (LEN=3) :: mpc_code 
! OUTPUT CHARACTER (LEN=2) :: our_code
!
!=========================================================

SUBROUTINE encode_station(mpc_code, our_code)
IMPLICIT NONE 
!INTERFACE
CHARACTER (LEN=3), INTENT(IN)  :: mpc_code 
CHARACTER (LEN=2), INTENT(OUT) :: our_code 
!END INTERFACE
INTEGER :: number, fb64_char_value
CHARACTER (LEN=4):: temp_code



  IF ( INDEX("0123456789", mpc_code(2:2))== 0 .OR. &
       INDEX("0123456789", mpc_code(3:3))== 0 ) THEN ! ERROR
     WRITE (*,*) "WARNING in encode_station: ", mpc_code , &
                   " can't be encoded  (last two characters are not digits!) "
     our_code = "**"
 ELSEIF ( INDEX("0123456789"// &
            "ABCDEFGHIJKLMNOPQRSTUV" // &
            "abcde" , mpc_code(1:1))== 0 ) THEN ! ERROR
    WRITE (*,*) "WARNING in encode_station: ", mpc_code , &
             " can't be encoded  (first digit not handled!) "
    our_code = "**"
 ELSE ! ALMOST RIGHT

    number = fb64_char_value(mpc_code(1:1)) * 100 + &
          fb64_char_value(mpc_code(2:2)) * 10  + &
          fb64_char_value(mpc_code(3:3))  

!WRITE (*,*) "DB:", number
   IF (number >= 4096) THEN ! ERROR
      WRITE (*,*) "WARNING in encode_station: ", mpc_code , " can't be encoded"
      our_code = "**"
   ELSE
      CALL base_10_to_64(number, temp_code)
      IF (LEN(TRIM(temp_code)) == 1) THEN
         our_code = "0"//temp_code(1:1)
      ELSEIF  (LEN(TRIM(temp_code)) == 2) THEN
         our_code = temp_code(1:2)
      END IF
   END IF
END IF

END SUBROUTINE encode_station



!=========================================================
!              decode_number
!=========================================================
! INPUT  CHARACTER (LEN=4) :: number_code 
! OUTPUT INTEGER :: number
!
!=========================================================

SUBROUTINE decode_number (number_code, number)
IMPLICIT NONE
!INTERFACE
CHARACTER (LEN=4), INTENT(IN) :: number_code
INTEGER, INTENT(OUT) :: number
!END INTERFACE

INTEGER :: fb64_char_value
CHARACTER (LEN=64):: basis_64="0123456789"// &
            "ABCDEFGHIJKLMNOPQRSTUVWXYZ"// &
            "abcdefghijklmnopqrstuvwxyz"// &
            "%&"

IF ((INDEX(basis_64, number_code(1:1)) == 0).OR. &
    (INDEX(basis_64, number_code(2:2)) == 0).OR. &
    (INDEX(basis_64, number_code(3:3)) == 0).OR. &
    (INDEX(basis_64, number_code(4:4)) == 0) ) THEN ! ERROR
   WRITE (*,*) "WARNING in decode_number: ->", number_code, "<- isn't a valid number code." 
   number= -1
ELSE ! RIGHT
   number = (fb64_char_value(number_code(1:1)) * 64**3 + &
             fb64_char_value(number_code(2:2)) * 64**2 + &
             fb64_char_value(number_code(3:3)) * 64    + &
             fb64_char_value(number_code(4:4)))        - &
             10*64**3  
END IF
END SUBROUTINE decode_number


!=========================================================
!              decode_date
!=========================================================
! INPUT  CHARACTER (LEN=3) :: date_code 
! OUTPUT INTEGER :: y, m, d
!
!=========================================================

SUBROUTINE decode_date(date_code, y, m, d)
IMPLICIT NONE
!INTERFACE
CHARACTER (LEN=3), INTENT(IN) :: date_code
INTEGER, INTENT(OUT) :: y, m, d
!END INTERFACE

DOUBLE PRECISION :: tjm, hour
!REAL (KIND=8) ::  tjm, hour
INTEGER :: fb64_char_value
CHARACTER (LEN=64):: basis_64="0123456789"// &
            "ABCDEFGHIJKLMNOPQRSTUVWXYZ"// &
            "abcdefghijklmnopqrstuvwxyz"// &
            "%&"

IF ((INDEX(basis_64, date_code(1:1)) == 0).OR. &
    (INDEX(basis_64, date_code(2:2)) == 0).OR. &
    (INDEX(basis_64, date_code(3:3)) == 0) ) THEN ! ERROR
   WRITE (*,*) "WARNING in decode_date: ->", date_code, "<- isn't a valid date code." 
   d= -1
   m= -1
   y= -1
   hour= -1
ELSE ! RIGHT
   tjm = (fb64_char_value(date_code(1:1)) * 64**2 + &
          fb64_char_value(date_code(2:2)) * 64    + &
          fb64_char_value(date_code(3:3))) *1.0 

   CALL mjddat(tjm, d, m, y, hour)
END IF

END SUBROUTINE decode_date

SUBROUTINE decode_date_mjd(date_code, tjm)
IMPLICIT NONE
!INTERFACE
CHARACTER (LEN=3), INTENT(IN) :: date_code
DOUBLE PRECISION, INTENT(OUT) :: tjm
!END INTERFACE

INTEGER :: fb64_char_value
CHARACTER (LEN=64):: basis_64="0123456789"// &
            "ABCDEFGHIJKLMNOPQRSTUVWXYZ"// &
            "abcdefghijklmnopqrstuvwxyz"// &
            "%&"

IF ((INDEX(basis_64, date_code(1:1)) == 0).OR. &
    (INDEX(basis_64, date_code(2:2)) == 0).OR. &
    (INDEX(basis_64, date_code(3:3)) == 0) ) THEN ! ERROR
   WRITE (*,*) "WARNING in decode_date_mjd: ->", date_code,    &
 &           "<- isn't a valid date code." 

ELSE ! RIGHT
   tjm = (fb64_char_value(date_code(1:1)) * 64**2 + &
          fb64_char_value(date_code(2:2)) * 64    + &
          fb64_char_value(date_code(3:3))) *1.0 
END IF

END SUBROUTINE decode_date_mjd

!=========================================================
!              decode_mjd
!=========================================================
! INPUT  CHARACTER (LEN=3) :: date_code 
! OUTPUT INTEGER :: mjd
!
!=========================================================

SUBROUTINE decode_mjd(date_code, tjm)
IMPLICIT NONE
!INTERFACE
CHARACTER (LEN=3), INTENT(IN) :: date_code
DOUBLE PRECISION, INTENT(OUT) :: tjm
!END INTERFACE

INTEGER :: fb64_char_value
CHARACTER (LEN=64):: basis_64="0123456789"// &
            "ABCDEFGHIJKLMNOPQRSTUVWXYZ"// &
            "abcdefghijklmnopqrstuvwxyz"// &
            "%&"

IF ((INDEX(basis_64, date_code(1:1)) == 0).OR. &
    (INDEX(basis_64, date_code(2:2)) == 0).OR. &
    (INDEX(basis_64, date_code(3:3)) == 0) ) THEN ! ERROR
   WRITE (*,*) "WARNING in decode_date: ->", date_code, "<- isn't a valid date code." 
ELSE ! RIGHT
   tjm = (fb64_char_value(date_code(1:1)) * 64.d0**2 + &
          fb64_char_value(date_code(2:2)) * 64.d0    + &
          fb64_char_value(date_code(3:3)))* 1.d0 
END IF

END SUBROUTINE decode_mjd

!=========================================================
!              decode_station
!=========================================================
! INPUT   CHARACTER (LEN=2) :: station_code
! OUTPUT  CHARACTER (LEN=3) :: station_mpc 
!
!=========================================================

SUBROUTINE decode_station(station_code, station_mpc)
IMPLICIT NONE
!INTERFACE
CHARACTER (LEN=2), INTENT(IN) :: station_code
CHARACTER (LEN=3), INTENT(OUT) :: station_mpc
!END INTERFACE 

INTEGER :: n
INTEGER :: fb64_char_value

INTEGER :: a, aa, aaa
CHARACTER (LEN=64):: basis_64="0123456789"// &
            "ABCDEFGHIJKLMNOPQRSTUVWXYZ"// &
            "abcdefghijklmnopqrstuvwxyz"// &
            "%&"

IF ((INDEX(basis_64, station_code(1:1)) == 0).OR. &
    (INDEX(basis_64, station_code(2:2)) == 0) ) THEN ! ERROR
   WRITE (*,*) "WARNING in decode_station: ->", station_code, "<- isn't a valid code"
   station_mpc = "***"
ELSE ! RIGHT
   n = fb64_char_value(station_code(1:1)) * 64 + &
        fb64_char_value(station_code(2:2))
!   WRITE(*,*) "number=", n   

   IF (n<10) THEN
      a = n
      station_mpc = "00"//basis_64(a+1:a+1)
   ELSEIF (n<100) THEN
      aa = INT((n*1.)/10.)
      a  = n - aa*10
      station_mpc = "0"//basis_64(aa+1:aa+1)//basis_64(a+1:a+1)
   ELSE
      aaa= INT((n*1.)/100.)
      aa = INT(((n - aaa *100)*1.)/10.)
      a  = n - aaa*100 - aa*10
      !WRITE(*,*) aaa, aa, a
      station_mpc=basis_64(aaa+1:aaa+1)//basis_64(aa+1:aa+1)//basis_64(a+1:a+1)
   END IF
END IF


END SUBROUTINE decode_station




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Auxiliar subroutines & functions !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!=========================================================
!              base_10_to_64    
!=========================================================
! INPUT            INTEGER :: b_10 < 64**4
! OUTPUT CHARACTER (LEN=4) :: b_64 
!
!=========================================================
SUBROUTINE base_10_to_64(b_10, b_64)
IMPLICIT NONE 
!INTERFACE
INTEGER, INTENT(IN):: b_10
CHARACTER (LEN=4), INTENT(OUT):: b_64
!END INTERFACE

INTEGER :: quotient, rest, temp
CHARACTER (LEN=64):: basis_64="0123456789"// &
     "ABCDEFGHIJKLMNOPQRSTUVWXYZ"// &
     "abcdefghijklmnopqrstuvwxyz"// &
     "%&"
IF ((b_10 >= 64**4).OR.(b_10 < 0)) THEN ! ERROR
   WRITE (*,*) "WARNING base_10_to_64: number ->", b_10, &
               "<- negative or too large!"
   b_64 = "****"
ELSE ! RIGHT
   b_64=""
   quotient = b_10
   DO
      IF (quotient < 64) EXIT 
      temp = INT((quotient*1.)/(64*1.))
      rest = quotient-temp *64
      quotient = temp
      
      b_64 = basis_64(rest+1:rest+1)//TRIM(b_64)
   END DO
   b_64 = basis_64(quotient+1:quotient+1)//TRIM(b_64)
END IF

END SUBROUTINE base_10_to_64

!=========================================================
!              base_64_to_10    
!=========================================================
! INPUT   CHARACTER (LEN=4) :: b_64 (char in 0..9A..Za..z%&)
! OUTPUT  INTEGER :: b_10 
!
!=========================================================

SUBROUTINE base_64_to_10(b_64, b_10)
IMPLICIT NONE 
!INTERFACE
CHARACTER (LEN=4), INTENT(IN):: b_64
INTEGER, INTENT(OUT):: b_10
!END INTERFACE

INTEGER :: n
INTEGER :: fb64_char_value
CHARACTER (LEN=64) :: char_set = "0123456789"// &
     "ABCDEFGHIJKLMNOPQRSTUVWXYZ"// &
     "abcdefghijklmnopqrstuvwxyz"// &
     "%&"

IF ( (LEN(b_64) /= 4) .OR. &
    (INDEX(char_set, b_64(1:1))==0) .OR. &
    (INDEX(char_set, b_64(2:2))==0) .OR. & 
    (INDEX(char_set, b_64(3:3))==0) .OR. &
    (INDEX(char_set, b_64(4:4))==0) ) THEN ! ERROR
!IF ( LEN(b_64) /= 4 ) THEN
WRITE (*,*) "base_64_to_10: Invalid string ->", b_64, "<-"
b_10 = -1
ELSE ! RIGHT 
b_10 = (fb64_char_value(b_64(1:1)) * 64**3 + &
         fb64_char_value(b_64(2:2)) * 64**2 + &
         fb64_char_value(b_64(3:3)) * 64    + &
         fb64_char_value(b_64(4:4)))       


!b_10 = 0
!DO n=0, LEN(TRIM(b_64))-1
!b_10 =  fb64_char_value(b_64(LEN(TRIM(b_64))-n:LEN(TRIM(b_64))-n)) * 64**n + b_10
!WRITE(*,*) b_64(LEN(TRIM(b_64))-n:LEN(TRIM(b_64))-n)
!END DO
END IF

END SUBROUTINE base_64_to_10




SUBROUTINE b64_char_value(char, value)
IMPLICIT NONE 

CHARACTER (LEN=1), INTENT(IN) :: char
INTEGER, INTENT(OUT) :: value

INTEGER :: temp_val

temp_val = IACHAR(char)


IF     ( (48<= temp_val).AND.( temp_val <=57) ) THEN ! 0-9 
value = temp_val - 48
ELSEIF ( (65<= temp_val).AND.( temp_val <=90 ) ) THEN ! A-Z
value = temp_val - 65 + 10
ELSEIF ( (97<= temp_val).AND.( temp_val <=122) ) THEN ! a-z
value = temp_val - 97 + 36
ELSEIF ( (37<= temp_val).AND.( temp_val <=38 ) ) THEN ! %-&
value = temp_val - 37 + 62
ELSE
value = -1
WRITE(*,*) "WARNING in b64_char_value: character ->", char, "<- not allowed"
END IF

END SUBROUTINE b64_char_value



!=========================================================
!                 fb64_char_value 
!=========================================================
! INPUT  CHARAACTER (LEN=1) :: char (in 0..9A..Za..z%&)
! OUTPUT CHARACTER (LEN=4) :: b_64 
!
!=========================================================

FUNCTION fb64_char_value(char)
IMPLICIT NONE 
!INTERFACE
CHARACTER (LEN=1), INTENT(IN) :: char
INTEGER :: fb64_char_value
!END INTERFACE

INTEGER :: temp_val



temp_val = IACHAR(char)

IF     ( (48<= temp_val).AND.( temp_val <=57) ) THEN ! 0-9 
fb64_char_value = temp_val - 48
ELSEIF ( (65<= temp_val).AND.( temp_val <=90 ) ) THEN ! A-Z
fb64_char_value = temp_val - 65 + 10
ELSEIF ( (97<= temp_val).AND.( temp_val <=122) ) THEN ! a-z
fb64_char_value = temp_val - 97 + 36
ELSEIF ( (37<= temp_val).AND.( temp_val <=38 ) ) THEN ! %-&
fb64_char_value = temp_val - 37 + 62
ELSE ! ERROR
fb64_char_value = -1
!WRITE(*,*) "WARNING in b64_char_value: character ->", char, "<- not allowed"
END IF

END FUNCTION fb64_char_value


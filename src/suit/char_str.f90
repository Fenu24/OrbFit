! MODULE CHAR_STR

MODULE char_str

! ---------------------------------------------------------------------
! former SPACES1.H and SPACES2.H
! List of characters to be considered as "spaces"
  INTEGER, PARAMETER :: nspac=2
! blank and horizontal TAB
  INTEGER, PARAMETER, DIMENSION(nspac) :: icspac=(/32 , 9/) 
  PUBLIC nspac,icspac
! ---------------------------------------------------------------------
! former PARCH.H
! Length of buffer character strings
  INTEGER, PARAMETER :: lchx=2048
  PUBLIC lchx
END MODULE char_str

! ==================LIBRARY char_str==========================          
!                                                                       
! note: many to be removed in fortran90/95                              
!                                                                       
! CONTAINS                                                              
! SUBROUTINES                                                           
! CHARACTER STRINGS:                                                    
!clench		length of a character string                                   
! rmsp		remove spaces from a a character string
! norstr	normal form of a character string                              
! strcnt	content of a character string                                  
! upcase	transforms a character string into uppercase                   
! locase	transforms a character string into lowercase                   
! stspli	split a character string using a given separator               
! isnum		tell whether a character string contains only digits           
! islett	tell whether a character string contains only letters          
! nitchs	number of items (separated by spaces) in a character string    
! spflds	split fields (separated by comma) in a character string        
! titast        strings with ast. names to use in graphics and file name
! sv2int         translation of string-valued keywords to integer       
! char_count    count how many times one character appears in a string

!      
! HEADERS:          
!      
! Copyright (C) 1996 by Mario Carpino (carpino@brera.mi.astro.it)       
! Version: August 9, 1996        
! --------------------------------------------------------------------- 
!      
!  ***************************************************************      
!  *                                                             *      
!  *                 L E N C H                                   *      
!  *                                                             *      
!  *       Computation of the length of a character string       *      
!  *     ignoring trailing spaces                                *      
!  *                                                             *      
!  ***************************************************************      
!      
INTEGER FUNCTION lench(c) 
  USE char_str
  IMPLICIT NONE      
  INTEGER icc,k 
  CHARACTER c*(*) 
  lench=LEN(c) 
1 CONTINUE 
  icc=ICHAR(c(lench:lench)) 
  DO k=1,nspac 
     IF(icc.eq.icspac(k)) GOTO 2 
  ENDDO
  RETURN      
2 lench=lench-1 
  IF(lench.le.0) RETURN 
  GOTO 1      
END FUNCTION lench
! Copyright (C) 1996 by Mario Carpino (carpino@brera.mi.astro.it)       
! Version: August 28, 1996       
! --------------------------------------------------------------------- 
!      
!  *****************************************************************    
!  *                                                               *    
!  *                           R M S P                             *    
!  *                                                               *    
!  *            Remove spaces from a character string              *    
!  *                                                               *    
!  *****************************************************************    
!                   
! INPUT:    C         -  Character string     
!      
! OUTPUT:   C         -  Character string without spaces   
!           L         -  Length of the output character string          
!      
SUBROUTINE rmsp(c,l) 
  USE char_str
  IMPLICIT NONE      
  INTEGER l,i,icc,k,l1 
  CHARACTER c*(*) 
  CHARACTER*(lchx) c1 
  l=LEN(c) 
  CALL chkpdf(l,lchx,'lchx') 
  IF(l.LE.0) STOP '**** rmsp: internal error (01) ****'     
  c1=' ' 
  l1=0   
  DO 1 i=1,l 
     icc=ICHAR(c(i:i)) 
     DO  k=1,nspac 
        IF(icc.EQ.icspac(k)) GOTO 1 
     ENDDO
     l1=l1+1 
     c1(l1:l1)=c(i:i) 
1 END DO     
  l=l1 
  IF(l.LE.0) THEN 
     c=' ' 
     RETURN 
  END IF     
  c=c1(1:l1) 
END SUBROUTINE rmsp
! Copyright (C) 1996 by Mario Carpino (carpino@brera.mi.astro.it)       
! Version: August 10, 1996       
! --------------------------------------------------------------------- 
!      
!  *****************************************************************    
!  *                                                               *    
!  *                         N O R S T R                           *    
!  *                                                               *    
!  *              Normal form of a character string                *    
!  *                                                               *    
!  *****************************************************************    
!      
! INPUT:    C         -  Character string     
!      
! OUTPUT:   C         -  Character string in normal form:  
!1) leading spaces are removed   
!2) tabs are substituted by blanks            
!3) multiple spaces are substituted by        
!   a single blank  
!4) quoted strings are preserved 
!           L         -  Length of the output character string          
!      
SUBROUTINE norstr(c,l)
  USE char_str 
  IMPLICIT NONE      
  CHARACTER*(*) c 
  CHARACTER*(lchx) c1 
  INTEGER l,i1,icc,k,l1,nsp,i 
  LOGICAL isspac,quoted  
  INTEGER lench 
  EXTERNAL lench 
  l=lench(c) 
  CALL chkpdf(l,lchx,'lchx') 
  IF(l.LE.0) THEN 
     c=' ' 
     RETURN 
  END IF     
  c1=' '    
  DO 1 i1=1,l 
     icc=ICHAR(c(i1:i1)) 
     DO k=1,nspac 
        IF(icc.EQ.icspac(k)) GOTO 1 
     ENDDO
     GOTO 2 
1 END DO
2 CONTINUE      
  l1=0 
  nsp=0 
  quoted=.false.      
  DO 3 i=i1,l    
     IF(c(i:i).EQ.'''') THEN 
        IF(quoted) THEN 
           quoted=.false. 
           nsp=0 
        ELSE 
           quoted=.true. 
        END IF
     END IF
     IF(quoted) GOTO 4 
     isspac=.false. 
     icc=ICHAR(c(i:i)) 
     DO 13 k=1,nspac 
        IF(icc.EQ.icspac(k)) THEN 
           isspac=.true. 
           GOTO 14 
        END IF
13   END DO
14   CONTINUE   
     IF(isspac) THEN 
        nsp=nsp+1 
        IF(nsp.EQ.1) THEN 
           l1=l1+1 
           c1(l1:l1)=' ' 
        END IF
        GOTO 3 
     END IF
            
4    l1=l1+1 
     c1(l1:l1)=c(i:i) 
     nsp=0 
            
3 END DO          
  l=l1 
  IF(l.LE.0) THEN 
     c=' ' 
     RETURN 
  END IF          
  c=c1(1:l1)         
END SUBROUTINE norstr
! Copyright (C) 1996 by Mario Carpino (carpino@brera.mi.astro.it)       
! Version: August 26, 1996                
! --------------------------------------------------------------------- 
!           
!  *****************************************************************    
!  *                                                               *    
!  *                         S T R C N T                           *    
!  *                                                               *    
!  *                Content of a character string                  *    
!  *               (enclosed within quotes or not)                 *    
!  *                                                               *    
!  *****************************************************************    
!           
! INPUT:    STRING    -  Input string     
!           
! OUTPUT:   CONT      -  String content   
!           REST      -  Rest of the input string (after the content)   
!           ERROR     -  Error flag       
!           
SUBROUTINE strcnt(string,cont,rest,error) 
  USE char_str
  IMPLICIT NONE  
  CHARACTER*(*) string,cont,rest 
  LOGICAL error 
  INTEGER ls,icc,i,k,i1,i2 
  LOGICAL quoted 
  CHARACTER*1 qch         
! Determine the beginning of the CONTENT string, by discarding leading  
! spaces and finding the position of the opening quotes.                
! I1 is the position of the first character of CONTENT;                 
! QCH is the quote character found        
  ls=LEN(string) 
  DO 1 i=1,ls 
     icc=ICHAR(string(i:i)) 
! Skip leading space characters           
     DO k=1,nspac 
        IF(icc.EQ.icspac(k)) GOTO 1 
     ENDDO
     i1=i 
! Look if the first non-space character is a quote                      
     IF(string(i1:i1).EQ.'''') THEN 
        quoted=.true. 
        qch='''' 
        i1=i1+1 
     ELSEIF(string(i1:i1).EQ.'"') THEN 
        quoted=.true. 
        qch='"' 
        i1=i1+1 
     ELSE 
        quoted=.false. 
     END IF
     GOTO 3 
1 END DO
! The input string contains only space characters                       
  error=.false. 
  cont=' ' 
  rest=' ' 
  RETURN 
3 CONTINUE 
! Determine the end of the CONTENT string, by looking for a             
! second occurrence of the quote character (QCH)                        
! I2 is the position of the last character of CONTENT                   
  IF(quoted) THEN 
     DO 4 i=i1,ls 
        IF(string(i:i).EQ.qch) THEN 
           i2=i-1 
           GOTO 7 
        END IF
4    ENDDO
! Missing closing quotes                  
     error=.true. 
     cont=' ' 
     rest=' ' 
     RETURN 
  ELSE 
     DO 5 i=i1,ls 
        icc=ICHAR(string(i:i)) 
        DO 6 k=1,nspac 
           IF(icc.EQ.icspac(k)) THEN 
              i2=i-1 
              GOTO 7 
           END IF
6       ENDDO 
5    ENDDO
     i2=ls 
  END IF
! Extract CONTENT and REST                
7 CONTINUE 
  IF(i2.LT.i1) THEN 
     cont=' ' 
  ELSE 
     IF(LEN(cont).LT.i2-i1+1)THEN 
        WRITE(*,*) LEN(cont), cont, i2,i1, string 
        STOP '**** strcnt: insufficient length for CONT ****' 
     ENDIF
     cont=string(i1:i2) 
  END IF
  error=.false. 
  IF(i2.GT.ls-2) THEN 
     rest=' ' 
  ELSE 
     IF(LEN(rest).LT.ls-i2-1)        &
     &        STOP '**** strcnt: insufficient length for REST ****'     
     rest=string(i2+2:) 
  END IF          
END SUBROUTINE strcnt
! Copyright (C) 1996 by Mario Carpino (carpino@brera.mi.astro.it)       
! Version: August 26, 1996                
! --------------------------------------------------------------------- 
!           
!  *****************************************************************    
!  *                                                               *    
!  *                         U P C A S E                           *    
!  *                                                               *    
!  *         Transforms a character string into uppercase          *    
!  *                                                               *    
!  *****************************************************************    
!           
! WARNING: assumes ASCII internal representation of characters          
!           
! INPUT:    STRING    -  Input string     
!           
! OUTPUT:   STRING    -  String transformed to uppercase                
!           
SUBROUTINE upcase(string) 
  IMPLICIT NONE 
  CHARACTER*(*) string 
  INTEGER i,icc 
  DO 1 i=1,LEN(string) 
     icc=ICHAR(string(i:i)) 
     IF(icc.GE.97 .AND. icc.LE.122) string(i:i)=CHAR(icc-32) 
1 END DO
END SUBROUTINE upcase
! Copyright (C) 1996 by Mario Carpino (carpino@brera.mi.astro.it)       
! Version: September 4, 1996              
! --------------------------------------------------------------------- 
!           
!  *****************************************************************    
!  *                                                               *    
!  *                         S T S P L I                           *    
!  *                                                               *    
!  *      Split a character string using a given separator         *    
!  *                                                               *    
!  *****************************************************************    
!           
! INPUT:    C         -  Character string 
!           SEP       -  Separator        
!           
! OUTPUT:   C         -  Character string without first item            
!           ITEM1     -  First item of the input string                 
!           NOSPLI    -  (LOG) if true, no splitting is possible        
!           
SUBROUTINE stspli(c,sep,item1,nospli) 
  IMPLICIT NONE           
  CHARACTER*(*) c,item1 
  CHARACTER sep*1,c1*200 
  LOGICAL nospli 
  INTEGER is           
  INTEGER lench 
  EXTERNAL lench        
  nospli=.false. 
  is=INDEX(c,sep) 
  IF(is.LE.0) THEN 
     IF(lench(c).LE.0) THEN 
        nospli=.true. 
        item1=' ' 
     ELSE 
        item1=c 
        c=' ' 
     END IF
  ELSE 
     item1=c(1:is-1) 
     c1=c(is+1:) 
     c=c1 
  END IF          
END SUBROUTINE stspli
! Copyright (C) 1996 by Mario Carpino (carpino@brera.mi.astro.it)       
! Version: August 26, 1996                
! --------------------------------------------------------------------- 
!           
!  *****************************************************************    
!  *                                                               *    
!  *                          I S N U M                            *    
!  *                                                               *    
!  *    Tells whether a character string contains only digits      *    
!  *                                                               *    
!  *****************************************************************    
!           
! WARNING: assumes ASCII internal representation of characters          
!           
! INPUT:    STRING    -  Input string     
!           
LOGICAL FUNCTION isnum(string) 
  IMPLICIT NONE           
  CHARACTER*(*) string 
  INTEGER i,icc 
  isnum=.false. 
  DO 1 i=1,LEN(string) 
     icc=ICHAR(string(i:i)) 
     IF(icc.LT.48 .OR. icc.GT.57) RETURN 
1 END DO
  isnum=.true.           
END FUNCTION isnum
! Copyright (C) 1996 by Mario Carpino (carpino@brera.mi.astro.it)       
! Version: August 26, 1996                
! --------------------------------------------------------------------- 
!           
!  *****************************************************************    
!  *                                                               *    
!  *                         I S L E T T                           *    
!  *                                                               *    
!  *    Tells whether a character string contains only letters     *    
!  *                          (a-z,A-Z)                            *    
!  *                                                               *    
!  *****************************************************************    
!           
! WARNING: assumes ASCII internal representation of characters          
!           
! INPUT:    STRING    -  Input string     
!           
LOGICAL FUNCTION islett(string) 
  IMPLICIT NONE     
  CHARACTER*(*) string 
  INTEGER i,icc 
  islett=.false. 
  DO 1 i=1,LEN(string) 
     icc=ICHAR(string(i:i)) 
     IF(icc.GE.65 .AND. icc.LE.90) GOTO 1 
     IF(icc.GE.97 .AND. icc.LE.122) GOTO 1 
     RETURN 
1 END DO
  islett=.true. 
END FUNCTION islett
! Copyright (C) 1996 by Mario Carpino (carpino@brera.mi.astro.it)       
! Version: August 10, 1996                
! --------------------------------------------------------------------- 
!           
!  *****************************************************************    
!  *                                                               *    
!  *                         N I T C H S                           *    
!  *                                                               *    
!  *  Number of items (separated by spaces) in a character string  *    
!  *                                                               *    
!  *****************************************************************    
!           
INTEGER FUNCTION nitchs(c)
  USE char_str 
  IMPLICIT NONE           
  INTEGER lc,i 
  CHARACTER*(*) c 
  CHARACTER*(lchx) c1 
  LOGICAL quoted 
  CALL chkpdf(LEN(c),lchx,'lchx') 
  c1=c 
  CALL norstr(c1,lc) 
  IF(lc.LE.0) THEN 
     nitchs=0 
     RETURN 
  END IF
  nitchs=1 
  quoted=.false. 
  DO 1 i=1,lc 
     IF(c1(i:i).EQ.'''') quoted=(.NOT.quoted) 
     IF(quoted) GOTO 1 
     IF(c1(i:i).EQ.' ') nitchs=nitchs+1 
1 END DO
END FUNCTION nitchs
! Copyright (C) 1996 by Mario Carpino (carpino@brera.mi.astro.it)       
! Version: August 26, 1996                
! --------------------------------------------------------------------- 
!           
!  *****************************************************************    
!  *                                                               *    
!  *                         L O C A S E                           *    
!  *                                                               *    
!  *         Transforms a character string into lowercase          *    
!  *                                                               *    
!  *****************************************************************    
!           
! WARNING: assumes ASCII internal representation of characters          
!           
! INPUT:    STRING    -  Input string     
!           
! OUTPUT:   STRING    -  String transformed to lowercase                
!           
SUBROUTINE locase(string) 
  IMPLICIT NONE           
  CHARACTER*(*) string 
  INTEGER i,icc 
  DO 1 i=1,LEN(string) 
     icc=ICHAR(string(i:i)) 
     IF(icc.GE.65 .AND. icc.LE.90) string(i:i)=CHAR(icc+32) 
1 END DO
END SUBROUTINE locase
! Copyright (C) 1997 by Mario Carpino (carpino@brera.mi.astro.it)       
! Version: October 13, 1997               
! --------------------------------------------------------------------- 
!           
!  ***************************************************************      
!  *                                                             *      
!  *                         S P F L D S                         *      
!  *                                                             *      
!  *       Split fields (separated by a comma) in a string       *      
!  *                                                             *      
!  ***************************************************************      
!           
! INPUT:    STRING    -  Input string     
!           NFX       -  Max dimension of FIELD vector                  
!           
! OUTPUT:   FIELD     -  Fields           
!           NF        -  Number of fields 
!           
SUBROUTINE spflds(string,field,nf,nfx) 
  IMPLICIT NONE 
  INTEGER nf,nfx 
  CHARACTER*(*) string,field(nfx) 
  INTEGER, PARAMETER :: nlx=200 
  INTEGER l 
  CHARACTER*(nlx) c,item1 
  LOGICAL nospli           
  INTEGER lench 
  EXTERNAL lench    
  l=lench(string) 
  IF(l.GT.nlx) STOP '**** spflds: l > nlx ****'    
  nf=0 
  c=string          
1 CONTINUE 
  CALL stspli(c,',',item1,nospli) 
  IF(nospli) RETURN 
  nf=nf+1 
  IF(nf.GT.nfx) STOP '**** spflds: nf > nfx ****' 
  field(nf)=item1 
  GOTO 1 
END SUBROUTINE spflds
! ====================================================                  
! TITAST handles asteroid names to form title string                    
! ====================================================                  
subroutine titast(iarc,astna0,astnap,titnam,filnam,let) 
  implicit none 
  integer iarc 
  character*(*) astna0,astnap 
  character*40 nam0,namp 
  character*(*) titnam 
  character*(*) filnam 
  integer le,le1,let,lench 
! ====================================================                  
!  shorten and remove blanks from name    
  nam0=astna0 
  CALL rmsp(nam0,le) 
  IF(lench(astnap).gt.0)THEN 
     namp=astnap 
     CALL rmsp(namp,le1) 
  ENDIF
  if(iarc.eq.1)then 
     titnam=astna0 
     filnam=nam0 
  elseif(iarc.eq.2)then 
     titnam=astnap 
     filnam=namp 
  elseif(iarc.eq.3)then 
     filnam=nam0(1:le)//'='//namp(1:le1) 
     let=le+le1+1 
     le=lench(astna0) 
     le1=lench(astnap) 
     titnam=astna0(1:le)//'='//astnap(1:le1) 
  else 
     write(*,*)'titast: this should not happen, iarc=',iarc 
  endif
  return 
END subroutine titast
! Copyright (C) 1996 by Mario Carpino (carpino@brera.mi.astro.it)       
! Version: November 15, 1999              
! --------------------------------------------------------------------- 
!           
!  *****************************************************************    
!  *                                                               *    
!  *                         I S B N U M                           *    
!  *                                                               *    
!  *    Tells whether a character string contains only digits      *    
!  *                 and/or blank characters                       *    
!  *                                                               *    
!  *****************************************************************    
!           
! WARNING: assumes ASCII internal representation of characters          
!           
! INPUT:    STRING    -  Input string     
!           
LOGICAL FUNCTION isbnum(string) 
  IMPLICIT NONE 
  CHARACTER*(*) string 
  INTEGER i,icc 
  isbnum=.false. 
  DO 1 i=1,LEN(string) 
     IF(string(i:i).EQ.' ') GOTO 1 
     icc=ICHAR(string(i:i)) 
     IF(icc.LT.48 .OR. icc.GT.57) RETURN 
1 END DO
  isbnum=.true. 
END FUNCTION isbnum
! Copyright (C) 2000 by Mario Carpino (carpino@brera.mi.astro.it)       
! Version: December 11, 2000              
! --------------------------------------------------------------------- 
!           
!  *****************************************************************    
!  *                                                               *    
!  *                         F I L S T R                           *    
!  *                                                               *    
!  *              Fill a string with another string                *    
!  *                                                               *    
!  *****************************************************************    
!           
! INPUT:    CIN       -  Input character string                         
!           LT        -  Total length of output string                  
!           NBI       -  Initial number of blanks (when possible)       
!           JUST      -  Justification (-1=left, 0=center, 1=right)     
!           
! OUTPUT:   COUT      -  Output character string                        
!           
SUBROUTINE filstr(cin,cout,lt,nbi,just) 
  IMPLICIT NONE 
  INTEGER lt,nbi,just 
  CHARACTER*(*) cin,cout 
  INTEGER li,nr1,nr2,nb1,nb2,nbt 
  CHARACTER*100 blank 
  INTEGER lench 
  EXTERNAL lench 
  blank=' ' 
  li=lench(cin) 
  ! NR1 = number of characters left in COUT after filling in CIN          
  nr1=lt-li 
  IF(nr1.LT.nbi) THEN 
     nb1=nr1 
     IF(nb1.GT.0) THEN 
        cout=blank(1:nb1)//cin 
     ELSE 
        cout=cin 
     END IF
     RETURN 
  END IF
! NR2 = number of characters left in COUT after filling in CIN and      
!       NBI initial blanks                
  nr2=lt-li-nbi 
  IF(nr2.LE.0) THEN 
     cout=blank(1:nbi)//cin 
     RETURN 
  END IF
  IF(just.EQ.0) THEN 
     nb2=nr2/2 
  ELSEIF(just.LT.0) THEN 
     nb2=0 
  ELSE 
     nb2=nr2 
  END IF
  nbt=nbi+nb2 
  IF(nbt.GT.0) THEN 
     cout=blank(1:nbt)//cin 
  ELSE 
     cout=cin 
  END IF
END SUBROUTINE filstr
! ==============================================================
!    REAL2STRING            
! Convert a real value to a string with given number of digits
! ==============================================================
SUBROUTINE real2string(value,ndigits_1,ndigits_2,string,error)
IMPLICIT NONE

DOUBLE PRECISION, INTENT(IN)  :: value      ! real value
INTEGER,          INTENT(IN)  :: ndigits_1  ! number of digits BEFORE decimal point
INTEGER,          INTENT(IN)  :: ndigits_2  ! number of digits AFTER  decimal point
CHARACTER(LEN=*), INTENT(OUT) :: string     ! character string
LOGICAL,          INTENT(OUT) :: error      ! conversion error

INTEGER :: required_length,lf,nd2
CHARACTER(LEN=20) :: format
!EXTERNAL :: rmsp

error=.true.
string=' '
nd2=MAX(ndigits_2,0)

required_length = ndigits_1 + nd2 + 1
IF(required_length > LEN(string)) RETURN

WRITE(format,100,ERR=10) required_length,nd2
100 FORMAT('(F',I3,'.',I3,')')
CALL rmsp(format,lf)
WRITE(string,FMT=format,ERR=10) value
IF(INDEX(string,'*') == 0) error=.false.
IF(string(1:1).eq.' ')THEN
   string(1:1)='0'
   IF(string(2:2).eq.' ')THEN
      string(2:2)='0'
   ENDIF
ENDIF

10 CONTINUE

END SUBROUTINE real2string

INTEGER FUNCTION char_count(string,le,cha)
  INTEGER, INTENT(IN) :: le ! length of string
  CHARACTER*(le), INTENT(IN) :: string
  CHARACTER*1, INTENT(IN) :: cha ! character to be counted in string
  INTEGER i
  char_count=0
  DO i=1,le
    IF(string(i:i).eq.cha) char_count=char_count+1
  ENDDO
END FUNCTION char_count

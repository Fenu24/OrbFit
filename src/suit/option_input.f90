! =======MODULE option_input========================================
MODULE option_input
! ---------------------------------------------------------------------
! former PARNAM.H
! Max number of namelist entries (main namelist)
  INTEGER, PARAMETER :: nnex=200
! Max number of namelist entries (simplified file-header namelist)
  INTEGER, PARAMETER :: nfnex=200
! Max number of entries in the list of keywords
  INTEGER, PARAMETER :: nklsx=500
! Max number of string to integer translations
  INTEGER, PARAMETER :: ns2itx=300
! Length of the "key" field
  INTEGER, PARAMETER :: lckx=30
! Length of the "value" field
  INTEGER, PARAMETER :: lcvx=2048
! Length of the "filename" field
  INTEGER, PARAMETER :: lcfx=300
  PUBLIC nnex,nfnex,nklsx,ns2itx,lckx,lcvx,lcfx

! ---------------------------------------------------------------------
! former COMNAM.H
! INPUT NAMELIST
! nne         -  Number of namelist items
! krecnm(i)   -  Record within input file (i=1,nne)
! kuord(i)    -  Order of use (i=1,nne)
! krtyp(i)    -  Keyword type:
!                      1 - integer
!                      2 - real
!                      3 - character string
!                      4 - logical
!                      5 - MJD
! kuorl       -  Last order of use
! iicnam      -  Initialization check
! keys(i)     -  Namelist keyword fields (i=1,nne)
! vals(i)     -  Namelist value fields (i=1,nne)
! namif(i)    -  Namelist input file (i=1,nne)
!
  INTEGER nne,krecnm(nnex),kuord(nnex),krtyp(nnex),kuorl,iicnam
  CHARACTER*(lckx), DIMENSION(nnex) :: keys
  CHARACTER*(lcvx), DIMENSION(nnex) :: vals
  CHARACTER*(lcfx), DIMENSION(nnex) :: namif
  PUBLIC nne,krecnm,kuord,krtyp,kuorl,iicnam,keys,vals,namif 
! ---------------------------------------------------------------------
! former COMFNM.H
! SIMPLIFIED FILE-HEADER NAMELIST
! nfne        -  Number of namelist items
! krcfnm(i)   -  Record within input file (i=1,nfne)
! kuorf(i)    -  Order of use (i=1,nfne)
! kuorlf      -  Last order of use
! iicfnm      -  Initialization check
! keysf(i)    -  Namelist keyword fields (i=1,nfne)
! valsf(i)    -  Namelist value fields (i=1,nfne)
! nmif        -  Namelist input file
! hnfuni      -  Namelist input unit
!
  INTEGER nfne,krcfnm(nfnex),kuorf(nfnex),kuorlf,hnfuni,iicfnm
  CHARACTER*(lckx), DIMENSION(nfnex) :: keysf
  CHARACTER*(lcvx), DIMENSION(nfnex) :: valsf
  CHARACTER nmif*(lcfx)
  PUBLIC nfne,krcfnm,kuorf,kuorlf,hnfuni,iicfnm,keysf,valsf,nmif

! ---------------------------------------------------------------------
! former COMKLS.H
! LIST OF RECOGNIZED KEYWORDS
!
! nkls        -  Number of entries in the list
! ns2it       -  Total number of string to integer translations 
! ns2i(i)     -  Number of string to integer translations 
!                for each keyword (i=1,nkls)
! ipos2i(i)   -  Position of the zero-th integer translation descriptors
!                for each keyword (i=1,nkls)
! intlst(k)   -  Integer value of the keyword (k=1,ns2it)
! keytyp      -  Keyword type:
!                      1 - integer
!                      2 - real
!                      3 - character string
!                      4 - logical
!                      5 - MJD
!                      6 - REF (reference system: see rdnref.f)
! iickls      -  Initialization check
! keylst      -  List of keywords
! vallst(k)   -  String value of the keyword (k=1,ns2it)
!
INTEGER nkls,ns2it,ns2i(nklsx),ipos2i(nklsx),intlst(ns2itx),keytyp(nklsx),iickls
CHARACTER*(lckx), DIMENSION(nklsx) :: keylst
CHARACTER*(lcvx), DIMENSION(ns2itx)::  vallst
PUBLIC nkls,ns2it,ns2i,ipos2i,intlst,keytyp,iickls,keylst,vallst

END MODULE option_input

! ======LIBRARY  option_input========================

!  MAIN OPTION NAMELIST  Carpino 1996-1998
! CONTAINS  
! ====================================================================  
!cnamini	initialization of namelist common block    
! rdnam		read a main option namelist      
! *rdnam1	read a main option namelist (INPUT level: 1)        
! *rdnam2	read a main option namelist (INPUT level: 2)        
! *rdnam3	read a main option namelist (INPUT level: 3)        
! rdklst	read the the list of valid keys  
! *rdkls1	read the the list of valid keys (INPUT level: 1)    
! *chkkey	check keyword validity
! rdncha	read a scalar character value    
! rdnint	read a scalar integer value      
! rdnlog	read a scalar logical value      
! rdnrea	read a scalar real value         
! rdntim	read a scalar time/date value    
! rdnref	read a reference system description        
! rdvint	read an integer vector of knowm dimension  
! rdvrea	read a real vector of knowm dimension      
! rdvcha	read a character vector of knowm dimension 
! rdmint	read an integer vector of unknowm dimension
! rdmrea	read a eeal vector of unknowm dimension    
! rdmcha	read a character vector of unknowm dimension         
! *getkv	locates a key/value pair in the input namelist       
! chkpdf        checks parameter definition         
! 
! HEADERS:  
! 
! 
! ================================================================      
! Copyright (C) 1996 by Mario Carpino (carpino@brera.mi.astro.it)       
! Version: August 10, 1996      
! --------------------------------------------------------------------- 
! 
!  *****************************************************************    
!  *   *    
!  *     N A M I N I       *    
!  *   *    
!  * Namelist common block initialization      *    
!  *   *    
!  *****************************************************************    
! 
SUBROUTINE namini 
  USE option_input
  IMPLICIT NONE 
! Common blocks to be initialized:        
  nne=0 
  kuorl=0 
  CALL chkpdf(1,nnex,'nnex') 
  iicnam=36 
END SUBROUTINE namini
! Copyright (C) 1997-1998 by Mario Carpino (carpino@brera.mi.astro.it)  
! Version: January 12, 1998     
! --------------------------------------------------------------------- 
! 
!  *****************************************************************    
!  *   *    
!  *      R D N A M        *    
!  *   *    
!  *     Reads a namelist from input file and stores in common     *    
!  *   *    
!  *****************************************************************    
! 
! INPUT:    IUN       -  Input FORTRAN unit         
! 
SUBROUTINE rdnam(iun) 
  USE option_input
  USE char_str
  IMPLICIT NONE 
  INTEGER iun 
  INCLUDE 'parcmc.h90' 
  CHARACTER*(lchx) rec,rec1,rec2,key1,key2,keyt,val1,infile,defcat 
  CHARACTER*(lchx) rest 
  CHARACTER*100 inpfl 
  LOGICAL opnd,error 
  INTEGER kr,ldc,lf,lr,ipc,lr2,iuna,ipu,lk,lv,ip1,ip2,lk2,i 
  INTEGER lench 
  EXTERNAL lench 
  IF(iicnam.NE.36) STOP '**** rdnam: internal error (01) ****' 
! Name of the input file        
  INQUIRE(iun,OPENED=opnd,NAME=infile) 
  IF(.NOT.opnd) STOP '**** rdnam: internal error (02) ****' 
  lf=lench(infile) 
  CALL chkpdf(lf,lcfx,'lcfx') 
  kr=0 
  defcat=' ' 
  ldc=0 
1 READ(iun,100,end=10) rec 
100 FORMAT(a) 
  kr=kr+1 
  rec1=rec 
  CALL rmsp(rec1,lr) 
  IF(lr.LE.0) GOTO 1 
! Compute length excluding comments       
  lr=lench(rec) 
!**   IF(lr.LT.1) GOTO 1        
  ipc=INDEX(rec(1:lr),comcha) 
  IF(ipc.EQ.1) GOTO 1 
  IF(ipc.EQ.0) THEN 
     rec1=rec(1:lr) 
  ELSE 
     rec1=rec(1:ipc-1) 
     lr=lench(rec1) 
     IF(lr.LT.1) GOTO 1 
  END IF
! Processing of "INPUT:" special keyword  
  rec2=rec1 
  CALL norstr(rec2,lr2) 
  IF(lr2.LT.6) GOTO 4 
  IF(rec2(1:6).EQ.'INPUT:') THEN 
     CALL strcnt(rec2(7:lr2),inpfl,rest,error) 
     IF(error) GOTO 21 
     CALL filopl(iuna,inpfl) 
     CALL rdnam1(iuna) 
     CALL filclo(iuna,' ') 
     GOTO 1 
  END IF
4 CONTINUE 
! Keyword field       
  ipu=INDEX(rec1(1:lr),'=') 
  IF(ipu.EQ.0) THEN 
     CALL rmsp(rec1,lk) 
     CALL chkfln(lk,lckx,'keyword',kr,infile) 
     key1=rec1(1:lk) 
     val1=' ' 
     GOTO 2 
  END IF
  key1=rec1(1:ipu-1) 
  CALL rmsp(key1,lk) 
  CALL chkfln(lk,lckx,'keyword',kr,infile) 
! Value field         
  val1=rec1(ipu+1:) 
  CALL norstr(val1,lv) 
  CALL chkfln(lv,lcvx,'value',kr,infile) 
2 CONTINUE 
! Handling of default category  
  IF(key1(lk:lk).EQ.'.') THEN 
     ldc=lk-1 
     defcat=key1(1:ldc) 
     GOTO 1 
  END IF
  IF(key1(1:1).EQ.'.') THEN 
     IF(ldc.LE.0) THEN 
        WRITE(*,101) key1(1:lk),infile(1:lf),kr 
        STOP '**** rdnam: abnormal end ****' 
     END IF
     keyt=defcat(1:ldc)//key1(1:lk) 
     key1=keyt 
     lk=lk+ldc 
  END IF
101 FORMAT(' ERROR: missing default category declaration'/  &
     &       '        ambiguous keyword "',a,'"'/   &
     &       '        (file "',a,'", record',i4,')')
! Check for parentheses         
  ip1=INDEX(key1,'(') 
  IF(ip1.EQ.0) THEN 
     ip2=INDEX(key1,')') 
     IF(ip2.NE.0) THEN 
        WRITE(*,102) key1(1:lk),infile(1:lf),kr 
        STOP '**** rdnam: abnormal end ****' 
     END IF
  ELSE 
     key2=key1(ip1+1:) 
     lk2=lench(key2) 
     IF(lk2.LE.0.OR.key2(lk2:lk2).NE.')') THEN 
        WRITE(*,102) key1(1:lk),infile(1:lf),kr 
        STOP '**** rdnam: abnormal end ****' 
     END IF
  END IF
102 FORMAT(' ERROR: illegal use of parentheses'/  &
     &       '        in keyword "',a,'"'/&
     &       '        (file "',a,'", record',i4,')')
! Look if the key is already present in the namelist
  DO 3 i=1,nne 
     IF(key1(1:lk).EQ.keys(i)) THEN 
        vals(i)=val1 
        namif(i)=infile 
        krecnm(i)=kr 
        kuord(i)=0 
        GOTO 1 
     END IF
3 END DO
  nne=nne+1 
  CALL chkpdf(nne,nnex,'nnex') 
  i=nne 
  keys(i)=key1 
  vals(i)=val1 
  namif(i)=infile 
  krecnm(i)=kr 
  kuord(i)=0 
  GOTO 1 
10 CONTINUE 
  RETURN 
! Error messages      
21 CONTINUE 
  WRITE(*,106) infile(1:lf),kr 
106 FORMAT(' ERROR: illegal "INPUT:" statement'/  &
         &       '        (file "',a,'", record',i4,')')
  STOP '**** rdnam: abnormal end ****' 
END SUBROUTINE rdnam
! Copyright (C) 1997-1998 by Mario Carpino (carpino@brera.mi.astro.it)  
! Version: January 12, 1998     
! --------------------------------------------------------------------- 
! 
!  *****************************************************************    
!  *   *    
!  *     R D N A M 1       *    
!  *   *    
!  *     Reads a namelist from input file and stores in common     *    
!  *   *    
!  *****************************************************************    
! 
! INPUT:    IUN       -  Input FORTRAN unit         
! 
SUBROUTINE rdnam1(iun) 
  USE option_input
  USE char_str
  IMPLICIT NONE 
  INTEGER iun 
  INCLUDE 'parcmc.h90' 
  CHARACTER*(lchx) rec,rec1,rec2,key1,key2,keyt,val1,infile,defcat 
  CHARACTER*(lchx) rest 
  CHARACTER*100 inpfl 
  LOGICAL opnd,error 
  INTEGER kr,ldc,lf,lr,ipc,lr2,iuna,ipu,lk,lv,ip1,ip2,lk2,i 
  INTEGER lench 
  EXTERNAL lench 
  IF(iicnam.NE.36) STOP '**** rdnam1: internal error (01) ****' 
! Name of the input file        
  INQUIRE(iun,OPENED=opnd,NAME=infile) 
  IF(.NOT.opnd) STOP '**** rdnam1: internal error (02) ****' 
  lf=lench(infile) 
  CALL chkpdf(lf,lcfx,'lcfx') 
  kr=0 
  defcat=' ' 
  ldc=0 
1 READ(iun,100,end=10) rec 
100 FORMAT(a) 
  kr=kr+1 
  rec1=rec 
  CALL rmsp(rec1,lr) 
  IF(lr.LE.0) GOTO 1 
! Compute length excluding comments       
  lr=lench(rec) 
!**   IF(lr.LT.1) GOTO 1        
  ipc=INDEX(rec(1:lr),comcha) 
  IF(ipc.EQ.1) GOTO 1 
  IF(ipc.EQ.0) THEN 
     rec1=rec(1:lr) 
  ELSE 
     rec1=rec(1:ipc-1) 
     lr=lench(rec1) 
     IF(lr.LT.1) GOTO 1 
  END IF
! Processing of "INPUT:" special keyword  
  rec2=rec1 
  CALL norstr(rec2,lr2) 
  IF(lr2.LT.6) GOTO 4 
  IF(rec2(1:6).EQ.'INPUT:') THEN 
     CALL strcnt(rec2(7:lr2),inpfl,rest,error) 
     IF(error) GOTO 21 
     CALL filopl(iuna,inpfl) 
     CALL rdnam2(iuna) 
     CALL filclo(iuna,' ') 
     GOTO 1 
  END IF
4 CONTINUE 
! Keyword field       
  ipu=INDEX(rec1(1:lr),'=') 
  IF(ipu.EQ.0) THEN 
     CALL rmsp(rec1,lk) 
     CALL chkfln(lk,lckx,'keyword',kr,infile) 
     key1=rec1(1:lk) 
     val1=' ' 
     GOTO 2 
  END IF
  key1=rec1(1:ipu-1) 
  CALL rmsp(key1,lk) 
  CALL chkfln(lk,lckx,'keyword',kr,infile) 
! Value field         
  val1=rec1(ipu+1:) 
  CALL norstr(val1,lv) 
  CALL chkfln(lv,lcvx,'value',kr,infile) 
2 CONTINUE 
! Handling of default category  
  IF(key1(lk:lk).EQ.'.') THEN 
     ldc=lk-1 
     defcat=key1(1:ldc) 
     GOTO 1 
  END IF
  IF(key1(1:1).EQ.'.') THEN 
     IF(ldc.LE.0) THEN 
        WRITE(*,101) key1(1:lk),infile(1:lf),kr 
        STOP '**** rdnam1: abnormal end ****' 
     END IF
     keyt=defcat(1:ldc)//key1(1:lk) 
     key1=keyt 
     lk=lk+ldc 
  END IF
101 FORMAT(' ERROR: missing default category declaration'/  &
     &       '        ambiguous keyword "',a,'"'/   &
     &       '        (file "',a,'", record',i4,')')
! Check for parentheses         
  ip1=INDEX(key1,'(') 
  IF(ip1.EQ.0) THEN 
     ip2=INDEX(key1,')') 
     IF(ip2.NE.0) THEN 
        WRITE(*,102) key1(1:lk),infile(1:lf),kr 
        STOP '**** rdnam1: abnormal end ****' 
     END IF
  ELSE 
     key2=key1(ip1+1:) 
     lk2=lench(key2) 
     IF(lk2.LE.0.OR.key2(lk2:lk2).NE.')') THEN 
        WRITE(*,102) key1(1:lk),infile(1:lf),kr 
        STOP '**** rdnam1: abnormal end ****' 
     END IF
  END IF
102 FORMAT(' ERROR: illegal use of parentheses'/  &
     &       '        in keyword "',a,'"'/&
     &       '        (file "',a,'", record',i4,')')
! Look if the key is already present in the namelist
  DO 3 i=1,nne 
     IF(key1(1:lk).EQ.keys(i)) THEN 
        vals(i)=val1 
        namif(i)=infile 
        krecnm(i)=kr 
        kuord(i)=0 
        GOTO 1 
     END IF
3 END DO
  nne=nne+1 
  CALL chkpdf(nne,nnex,'nnex') 
  i=nne 
  keys(i)=key1 
  vals(i)=val1 
  namif(i)=infile 
  krecnm(i)=kr 
  kuord(i)=0 
  GOTO 1 
10 CONTINUE 
  RETURN 
! Error messages      
21 CONTINUE 
  WRITE(*,106) infile(1:lf),kr 
106 FORMAT(' ERROR: illegal "INPUT:" statement'/  &
         &       '        (file "',a,'", record',i4,')')
  STOP '**** rdnam1: abnormal end ****' 
END SUBROUTINE rdnam1
! Copyright (C) 1997-1998 by Mario Carpino (carpino@brera.mi.astro.it)  
! Version: January 12, 1998     
! --------------------------------------------------------------------- 
! 
!  *****************************************************************    
!  *   *    
!  *     R D N A M 2       *    
!  *   *    
!  *     Read a namelist from input file and stores in common      *    
!  *   *    
!  *****************************************************************    
! 
! INPUT:    IUN       -  Input FORTRAN unit         
! 
SUBROUTINE rdnam2(iun) 
  USE option_input
  USE char_str
  IMPLICIT NONE 
  INTEGER iun 
  INCLUDE 'parcmc.h90' 
  CHARACTER*(lchx) rec,rec1,rec2,key1,key2,keyt,val1,infile,defcat 
  CHARACTER*(lchx) rest 
  CHARACTER*100 inpfl 
  LOGICAL opnd,error 
  INTEGER kr,ldc,lf,lr,ipc,lr2,iuna,ipu,lk,lv,ip1,ip2,lk2,i 
  INTEGER lench 
  EXTERNAL lench 
  IF(iicnam.NE.36) STOP '**** rdnam2: internal error (01) ****' 
! Name of the input file        
  INQUIRE(iun,OPENED=opnd,NAME=infile) 
  IF(.NOT.opnd) STOP '**** rdnam2: internal error (02) ****' 
  lf=lench(infile) 
  CALL chkpdf(lf,lcfx,'lcfx') 
  
  kr=0 
  defcat=' ' 
  ldc=0 
1 READ(iun,100,end=10) rec 
100 FORMAT(a) 
  kr=kr+1 
  
  rec1=rec 
  CALL rmsp(rec1,lr) 
  IF(lr.LE.0) GOTO 1 
  
! Compute length excluding comments       
  lr=lench(rec) 
!**   IF(lr.LT.1) GOTO 1        
  ipc=INDEX(rec(1:lr),comcha) 
  IF(ipc.EQ.1) GOTO 1 
  IF(ipc.EQ.0) THEN 
     rec1=rec(1:lr) 
  ELSE 
     rec1=rec(1:ipc-1) 
     lr=lench(rec1) 
     IF(lr.LT.1) GOTO 1 
  END IF
  
! Processing of "INPUT:" special keyword  
  rec2=rec1 
  CALL norstr(rec2,lr2) 
  IF(lr2.LT.6) GOTO 4 
  IF(rec2(1:6).EQ.'INPUT:') THEN 
     CALL strcnt(rec2(7:lr2),inpfl,rest,error) 
     IF(error) GOTO 21 
     CALL filopl(iuna,inpfl) 
     CALL rdnam3(iuna) 
     CALL filclo(iuna,' ') 
     GOTO 1 
  END IF
4 CONTINUE 
  
! Keyword field       
  ipu=INDEX(rec1(1:lr),'=') 
  IF(ipu.EQ.0) THEN 
     CALL rmsp(rec1,lk) 
     CALL chkfln(lk,lckx,'keyword',kr,infile) 
     key1=rec1(1:lk) 
     val1=' ' 
     GOTO 2 
  END IF
  key1=rec1(1:ipu-1) 
  CALL rmsp(key1,lk) 
  CALL chkfln(lk,lckx,'keyword',kr,infile) 
  
! Value field         
  val1=rec1(ipu+1:) 
  CALL norstr(val1,lv) 
  CALL chkfln(lv,lcvx,'value',kr,infile) 
2 CONTINUE 
  
! Handling of default category  
  IF(key1(lk:lk).EQ.'.') THEN 
     ldc=lk-1 
     defcat=key1(1:ldc) 
     GOTO 1 
  END IF
  IF(key1(1:1).EQ.'.') THEN 
     IF(ldc.LE.0) THEN 
        WRITE(*,101) key1(1:lk),infile(1:lf),kr 
        STOP '**** rdnam2: abnormal end ****' 
     END IF
     keyt=defcat(1:ldc)//key1(1:lk) 
     key1=keyt 
     lk=lk+ldc 
  END IF
101 FORMAT(' ERROR: missing default category declaration'/  &
     &       '        ambiguous keyword "',a,'"'/   &
     &       '        (file "',a,'", record',i4,')')
  
! Check for parentheses         
  ip1=INDEX(key1,'(') 
  IF(ip1.EQ.0) THEN 
     ip2=INDEX(key1,')') 
     IF(ip2.NE.0) THEN 
        WRITE(*,102) key1(1:lk),infile(1:lf),kr 
        STOP '**** rdnam2: abnormal end ****' 
     END IF
  ELSE 
     key2=key1(ip1+1:) 
     lk2=lench(key2) 
     IF(lk2.LE.0.OR.key2(lk2:lk2).NE.')') THEN 
        WRITE(*,102) key1(1:lk),infile(1:lf),kr 
        STOP '**** rdnam2: abnormal end ****' 
     END IF
  END IF
102 FORMAT(' ERROR: illegal use of parentheses'/  &
     &       '        in keyword "',a,'"'/&
     &       '        (file "',a,'", record',i4,')')
  
! Look if the key is already present in the namelist
  DO 3 i=1,nne 
     IF(key1(1:lk).EQ.keys(i)) THEN 
        vals(i)=val1 
        namif(i)=infile 
        krecnm(i)=kr 
        kuord(i)=0 
        GOTO 1 
     END IF
3 END DO
  
  nne=nne+1 
  CALL chkpdf(nne,nnex,'nnex') 
  i=nne 
  keys(i)=key1 
  vals(i)=val1 
  namif(i)=infile 
  krecnm(i)=kr 
  kuord(i)=0 
  GOTO 1 
  
10 CONTINUE 
  RETURN 
! Error messages      
21 CONTINUE 
  WRITE(*,106) infile(1:lf),kr 
106 FORMAT(' ERROR: illegal "INPUT:" statement'/  &
     &       '        (file "',a,'", record',i4,')')
  STOP '**** rdnam2: abnormal end ****'   
END SUBROUTINE rdnam2
! Copyright (C) 1996-1998 by Mario Carpino (carpino@brera.mi.astro.it)  
! Version: January 12, 1998     
! --------------------------------------------------------------------- 
! 
!  *****************************************************************    
!  *   *    
!  *     R D N A M 3       *    
!  *   *    
!  *     Read a namelist from input file and stores in common      *    
!  *   *    
!  *****************************************************************    
! 
! INPUT:    IUN       -  Input FORTRAN unit         
! 
SUBROUTINE rdnam3(iun) 
  USE option_input
  USE char_str
  IMPLICIT NONE 
  INTEGER iun 
  INCLUDE 'parcmc.h90' 
  CHARACTER*(lchx) rec,rec1,rec2,key1,key2,keyt,val1,infile,defcat 
  LOGICAL opnd 
  INTEGER kr,ldc,lf,lr,ipc,lr2,ipu,lk,lv,ip1,ip2,lk2,i 
  INTEGER lench 
  EXTERNAL lench 
  IF(iicnam.NE.36) STOP '**** rdnam3: internal error (01) ****' 
  
! Name of the input file        
  INQUIRE(iun,OPENED=opnd,NAME=infile) 
  IF(.NOT.opnd) STOP '**** rdnam3: internal error (02) ****' 
  lf=lench(infile) 
  CALL chkpdf(lf,lcfx,'lcfx') 
  
  kr=0 
  defcat=' ' 
  ldc=0 
1 READ(iun,100,end=10) rec 
100 FORMAT(a) 
  kr=kr+1 
  
  rec1=rec 
  CALL rmsp(rec1,lr) 
  IF(lr.LE.0) GOTO 1 
  
! Compute length excluding comments       
  lr=lench(rec) 
!**   IF(lr.LT.1) GOTO 1        
  ipc=INDEX(rec(1:lr),comcha) 
  IF(ipc.EQ.1) GOTO 1 
  IF(ipc.EQ.0) THEN 
     rec1=rec(1:lr) 
  ELSE 
     rec1=rec(1:ipc-1) 
     lr=lench(rec1) 
     IF(lr.LT.1) GOTO 1 
  END IF
  
! Processing of "INPUT:" special keyword  
  rec2=rec1 
  CALL norstr(rec2,lr2) 
  IF(lr2.LT.6) GOTO 4 
  IF(rec2(1:6).EQ.'INPUT:') &
       &  STOP '**** rdnam3: too many INPUT levels ****'        
4 CONTINUE 
  
! Keyword field       
  ipu=INDEX(rec1(1:lr),'=') 
  IF(ipu.EQ.0) THEN 
     CALL rmsp(rec1,lk) 
     CALL chkfln(lk,lckx,'keyword',kr,infile) 
     key1=rec1(1:lk) 
     val1=' ' 
     GOTO 2 
  END IF
  key1=rec1(1:ipu-1) 
  CALL rmsp(key1,lk) 
  CALL chkfln(lk,lckx,'keyword',kr,infile) 
  
! Value field         
  val1=rec1(ipu+1:) 
  CALL norstr(val1,lv) 
  CALL chkfln(lv,lcvx,'value',kr,infile) 
2 CONTINUE 
  
! Handling of default category  
  IF(key1(lk:lk).EQ.'.') THEN 
     ldc=lk-1 
     defcat=key1(1:ldc) 
     GOTO 1 
  END IF
  IF(key1(1:1).EQ.'.') THEN 
     IF(ldc.LE.0) THEN 
        WRITE(*,101) key1(1:lk),infile(1:lf),kr 
        STOP '**** rdnam3: abnormal end ****' 
     END IF
     keyt=defcat(1:ldc)//key1(1:lk) 
     key1=keyt 
     lk=lk+ldc 
  END IF
101 FORMAT(' ERROR: missing default category declaration'/  &
     &       '        ambiguous keyword "',a,'"'/   &
     &       '        (file "',a,'", record',i4,')')
! Check for parentheses         
  ip1=INDEX(key1,'(') 
  IF(ip1.EQ.0) THEN 
     ip2=INDEX(key1,')') 
     IF(ip2.NE.0) THEN 
        WRITE(*,102) key1(1:lk),infile(1:lf),kr 
        STOP '**** rdnam3: abnormal end ****' 
     END IF
  ELSE 
     key2=key1(ip1+1:) 
     lk2=lench(key2) 
     IF(lk2.LE.0.OR.key2(lk2:lk2).NE.')') THEN 
        WRITE(*,102) key1(1:lk),infile(1:lf),kr 
        STOP '**** rdnam3: abnormal end ****' 
     END IF
  END IF
102 FORMAT(' ERROR: illegal use of parentheses'/  &
     &       '        in keyword "',a,'"'/&
     &       '        (file "',a,'", record',i4,')')
! Look if the key is already present in the namelist
  DO 3 i=1,nne 
     IF(key1(1:lk).EQ.keys(i)) THEN 
        vals(i)=val1 
        namif(i)=infile 
        krecnm(i)=kr 
        kuord(i)=0 
        GOTO 1 
     END IF
3 END DO
  nne=nne+1 
  CALL chkpdf(nne,nnex,'nnex') 
  i=nne 
  keys(i)=key1 
  vals(i)=val1 
  namif(i)=infile 
  krecnm(i)=kr 
  kuord(i)=0 
  GOTO 1 
  
10 CONTINUE 
  RETURN 
! Error messages      
21 CONTINUE 
  WRITE(*,106) infile(1:lf),kr 
106 FORMAT(' ERROR: illegal "INPUT:" statement'/  &
     &       '        (file "',a,'", record',i4,')')
  STOP '**** rdnam3: abnormal end ****' 
END SUBROUTINE rdnam3
! Copyright (C) 1997-1998 by Mario Carpino (carpino@brera.mi.astro.it)  
! Version: January 13, 1998     
! --------------------------------------------------------------------- 
! 
!  *****************************************************************    
!  *   *    
!  *     R D K L S T       *    
!  *   *    
!  *     Read a the list of valid keys         *    
!  *   (with optional translation of string values to integers)    *    
!  *   *    
!  *****************************************************************    
! 
! INPUT:    FILE      -  Input file       
! 
SUBROUTINE rdklst(file) 
  USE option_input
  USE char_str
  IMPLICIT NONE 
  INCLUDE 'parcmc.h90' 
  CHARACTER*(*) file 
  CHARACTER*(lchx) rec,rest,rec1,key1,key2,stval,defcat,kt 
  CHARACTER*3 ktyp 
  CHARACTER*100 file1 
  LOGICAL nospli,error 
  INTEGER kr,ldc,iun,lv,lr,ipc,lk,ni,i,lkt 
  
  INTEGER lench,nitchs 
  EXTERNAL lench,nitchs 
  
  nkls=0 
  ns2it=0 
  kr=0 
  defcat=' ' 
  ldc=0 
  CALL filopl(iun,file) 
  
1 READ(iun,100,end=10) rec 
100 FORMAT(a) 
  kr=kr+1 
  
! Handling of integer translations of string values 
  IF(rec(1:2).EQ.'$I') THEN 
     IF(nkls.LE.0) STOP '**** rdklst: internal error (01) ****' 
     IF(keytyp(nkls).NE.3) THEN 
        WRITE(*,103) kr 
        STOP '**** rdklst: abnormal end ****' 
     END IF
     ns2it=ns2it+1 
     CALL chkpdf(ns2it,ns2itx,'ns2itx') 
     rest=rec(3:) 
     rec=rest 
     CALL strcnt(rec,stval,rest,error) 
     IF(error) STOP '**** rdklst: internal error (02) ****' 
     rec=rest 
     READ(rec,*,ERR=11) intlst(ns2it) 
     lv=lench(stval) 
     CALL chkpdf(lv,lcvx,'lcvx') 
     vallst(ns2it)=stval(1:lv) 
     ns2i(nkls)=ns2i(nkls)+1 
     GOTO 1 
  END IF
103 FORMAT(' **** rdklst: $I following a non-string keyword',         &
     &       'at line',i4,' ****')        
  
! Compute length excluding comments       
  lr=lench(rec) 
  IF(lr.LT.1) GOTO 1 
  ipc=INDEX(rec(1:lr),comcha) 
  IF(ipc.EQ.1) GOTO 1 
  IF(ipc.EQ.0) THEN 
     rec1=rec(1:lr) 
  ELSE 
     rec1=rec(1:ipc-1) 
     lr=lench(rec1) 
     IF(lr.LT.1) GOTO 1 
  END IF
  CALL norstr(rec1,lk) 
  
! Handling of "INPUT:" special keyword    
  IF(lk.LT.6) GOTO 2 
  IF(rec1(1:6).EQ.'INPUT:') THEN 
     CALL strcnt(rec1(7:),file1,rest,error) 
     IF(error) THEN 
        WRITE(*,105) kr 
        STOP '**** rdklst: abnormal end ****' 
     END IF
     CALL rdkls1(file1) 
     GOTO 1 
  END IF
105 FORMAT(' **** rdklst: input error at line',I4,' ****') 
2 CONTINUE 
  
! Keyword + type      
  CALL chkpdf(lk,lckx,'lckx') 
  kt=rec1(1:lk) 
  ni=nitchs(kt) 
  IF(ni.EQ.1) THEN 
     key1=kt(1:lk) 
     ktyp=' ' 
  ELSEIF(ni.EQ.2) THEN 
     CALL stspli(kt,' ',key1,nospli) 
     IF(nospli) STOP '**** rdklst: internal error (04) ****' 
     ktyp=kt 
  ELSE 
     WRITE(*,101) ni,kr 
     STOP '**** rdklst: abnormal end ****' 
  END IF
101 FORMAT(' **** rdklst: INTERNAL ERROR:',i3,' items at line',i5,    &
     &       ' ****') 
  
! Handling of default category  
  IF(ni.EQ.1) THEN 
     IF(key1(lk:lk).NE.'.') THEN 
        WRITE(*,102) kr 
        STOP '**** rdklst: abnormal end ****' 
     END IF
     ldc=lk-1 
     defcat=key1(1:ldc) 
     GOTO 1 
  END IF
  IF(key1(1:1).EQ.'.') THEN 
     IF(ldc.LE.0) STOP '**** rdklst: internal error (05) ****' 
     key2=defcat(1:ldc)//key1(1:lk) 
     key1=key2 
     lk=lk+ldc 
  END IF
102 FORMAT(' **** rdklst: ERROR: defcat without ending dot at line',  &
     &       i5,' ****')        
  
! Look if the key is already present in the namelist
  DO 3 i=1,nkls 
     IF(key1(1:lk).EQ.keylst(i)) THEN 
        WRITE(*,110) key1(1:lk),kr 
        STOP '**** rdklst: abnormal end ****' 
     END IF
3 END DO
110 FORMAT(' **** rdklst: duplicate key ****'/    &
     &       '      KEY=',a,' (line',i4,')')        
  
! Add the new keyword to the list         
  nkls=nkls+1 
  CALL chkpdf(nkls,nklsx,'nklsx') 
  keylst(nkls)=key1 
  ns2i(nkls)=0 
  ipos2i(nkls)=ns2it 
  IF(ktyp.EQ.'INT') THEN 
     keytyp(nkls)=1 
  ELSEIF(ktyp.EQ.'REA') THEN 
     keytyp(nkls)=2 
  ELSEIF(ktyp.EQ.'CHA') THEN 
     keytyp(nkls)=3 
  ELSEIF(ktyp.EQ.'LOG') THEN 
     keytyp(nkls)=4 
  ELSEIF(ktyp.EQ.'MJD' .OR. ktyp.EQ.'TIM') THEN 
     keytyp(nkls)=5 
  ELSEIF(ktyp.EQ.'REF') THEN 
     keytyp(nkls)=6 
  ELSE 
     lkt=lench(ktyp) 
     WRITE(*,111) ktyp(1:lkt),kr 
     STOP '**** rdklst: abnormal end ****' 
  END IF
111 FORMAT(' **** rdklst: unknown key type ****'/ &
     &       '      KTYPE="',a,'" (line',i4,')')    
  GOTO 1 
10 CONTINUE 
  iickls=36 
  CALL filclo(iun,' ') 
  RETURN 
11 CONTINUE 
  STOP '**** rdklst: internal error (03) ****' 
END SUBROUTINE rdklst
! Copyright (C) 1997-1998 by Mario Carpino (carpino@brera.mi.astro.it)  
! Version: January 13, 1998     
! --------------------------------------------------------------------- 
! 
!  *****************************************************************    
!  *   *    
!  *     R D K L S 1       *    
!  *   *    
!  *     Read a the list of valid keys         *    
!  *   (with optional translation of string values to integers)    *    
!  *   *    
!  *****************************************************************    
! 
! INPUT:    FILE      -  Input file       
! 
SUBROUTINE rdkls1(file)
  USE option_input
  USE char_str 
  IMPLICIT NONE 
  INCLUDE 'parcmc.h90' 
  CHARACTER*(*) file 
  CHARACTER*(lchx) rec,rest,rec1,key1,key2,stval,defcat,kt 
  CHARACTER*3 ktyp 
  LOGICAL nospli,error 
  INTEGER kr,ldc,iun,lv,lr,ipc,lk,ni,i,lkt,nkls1 
  
  INTEGER lench,nitchs 
  EXTERNAL lench,nitchs 
  
  kr=0 
  defcat=' ' 
  nkls1=0 
  ldc=0 
  CALL filopl(iun,file) 
  
1 READ(iun,100,end=10) rec 
100 FORMAT(a) 
  kr=kr+1 
  
! Handling of integer translations of string values 
  IF(rec(1:2).EQ.'$I') THEN 
     IF(nkls1.LE.0) STOP '**** rdkls1: internal error (01) ****' 
     IF(nkls.LE.0) STOP '**** rdkls1: internal error (01) ****' 
     IF(keytyp(nkls).NE.3) THEN 
        WRITE(*,103) kr 
        STOP '**** rdkls1: abnormal end ****' 
     END IF
     ns2it=ns2it+1 
     CALL chkpdf(ns2it,ns2itx,'ns2itx') 
     rest=rec(3:) 
     rec=rest 
     CALL strcnt(rec,stval,rest,error) 
     IF(error) STOP '**** rdkls1: internal error (02) ****' 
     rec=rest 
     READ(rec,*,ERR=11) intlst(ns2it) 
     lv=lench(stval) 
     CALL chkpdf(lv,lcvx,'lcvx') 
     vallst(ns2it)=stval(1:lv) 
     ns2i(nkls)=ns2i(nkls)+1 
     GOTO 1 
  END IF
103 FORMAT(' **** rdkls1: $I following a non-string keyword',         &
     &       'at line',i4,' ****')        
  
! Compute length excluding comments       
  lr=lench(rec) 
  IF(lr.LT.1) GOTO 1 
  ipc=INDEX(rec(1:lr),comcha) 
  IF(ipc.EQ.1) GOTO 1 
  IF(ipc.EQ.0) THEN 
     rec1=rec(1:lr) 
  ELSE 
     rec1=rec(1:ipc-1) 
     lr=lench(rec1) 
     IF(lr.LT.1) GOTO 1 
  END IF
  CALL norstr(rec1,lk) 
  
! Handling of "INPUT:" directive
  IF(lk.LT.6) GOTO 2 
  IF(rec1(1:6).EQ.'INPUT:') &
       &    STOP '**** rdkls1: too many nested INPUT levels ****'         
2 CONTINUE 
  
! Keyword + type      
  CALL chkpdf(lk,lckx,'lckx') 
  kt=rec1(1:lk) 
  ni=nitchs(kt) 
  IF(ni.EQ.1) THEN 
     key1=kt(1:lk) 
     ktyp=' ' 
  ELSEIF(ni.EQ.2) THEN 
     CALL stspli(kt,' ',key1,nospli) 
     IF(nospli) STOP '**** rdkls1: internal error (04) ****' 
     ktyp=kt 
  ELSE 
     WRITE(*,101) ni,kr 
     STOP '**** rdkls1: abnormal end ****' 
  END IF
101 FORMAT(' **** rdkls1: ERROR:',i3,' items at line',i5,' ****') 
  
! Handling of default category  
  IF(ni.EQ.1) THEN 
     IF(key1(lk:lk).NE.'.') THEN 
        WRITE(*,102) kr 
        STOP '**** rdkls1: abnormal end ****' 
     END IF
     ldc=lk-1 
     defcat=key1(1:ldc) 
     GOTO 1 
  END IF
  IF(key1(1:1).EQ.'.') THEN 
     IF(ldc.LE.0) STOP '**** rdkls1: missing defcat ****' 
     key2=defcat(1:ldc)//key1(1:lk) 
     key1=key2 
     lk=lk+ldc 
  END IF
102 FORMAT(' **** rdkls1: ERROR: defcat without ending dot at line',  &
     &       i5,' ****')        
  
! Look if the key is already present in the namelist
  DO 3 i=1,nkls 
     IF(key1(1:lk).EQ.keylst(i)) THEN 
        WRITE(*,110) key1(1:lk),kr 
        STOP '**** rdkls1: abnormal end ****' 
     END IF
3 END DO
110 FORMAT(' **** rdkls1: duplicate key ****'/    &
     &       '      KEY=',a,' (line',i4,')')        
  
! Add the new keyword to the list         
  nkls=nkls+1 
  nkls1=nkls1+1 
  CALL chkpdf(nkls,nklsx,'nklsx') 
  keylst(nkls)=key1 
  ns2i(nkls)=0 
  ipos2i(nkls)=ns2it 
  IF(ktyp.EQ.'INT') THEN 
     keytyp(nkls)=1 
  ELSEIF(ktyp.EQ.'REA') THEN 
     keytyp(nkls)=2 
  ELSEIF(ktyp.EQ.'CHA') THEN 
     keytyp(nkls)=3 
  ELSEIF(ktyp.EQ.'LOG') THEN 
     keytyp(nkls)=4 
  ELSEIF(ktyp.EQ.'MJD' .OR. ktyp.EQ.'TIM') THEN 
     keytyp(nkls)=5 
  ELSEIF(ktyp.EQ.'REF') THEN 
     keytyp(nkls)=6 
  ELSE 
     lkt=lench(ktyp) 
     WRITE(*,111) ktyp(1:lkt),kr 
     STOP '**** rdkls1: abnormal end ****' 
  END IF
111 FORMAT(' **** rdkls1: unknown key type ****'/ &
     &       '      KTYPE="',a,'" (line',i4,')')    
  GOTO 1 
  
10 CONTINUE 
  CALL filclo(iun,' ') 
  RETURN 
  
11 CONTINUE 
  STOP '**** rdkls1: internal error (03) ****' 
  
END SUBROUTINE rdkls1
! Copyright (C) 1996 by Mario Carpino (carpino@brera.mi.astro.it)       
! Version: September 19, 1996   
! --------------------------------------------------------------------- 
! 
!  *****************************************************************    
!  *   *    
!  *     C H K K E Y       *    
!  *   *    
!  *         Check keyword validity  *    
!  *    and transfer information on keyword type to namelist       *    
!  *   *    
!  *****************************************************************    
! 
SUBROUTINE chkkey 
  USE option_input
  IMPLICIT NONE 
  INTEGER lk2(nklsx),k,i,lk1,ipp,lf 
  CHARACTER keywp*(lckx) 
  LOGICAL fail,pare 
  
  INTEGER lench 
  EXTERNAL lench 
  
  fail=.false. 
  IF(iickls.NE.36) STOP '**** chkkey: internal error (01) ****' 
  IF(iicnam.NE.36) STOP '**** chkkey: internal error (02) ****' 
  
  DO  k=1,nkls 
    lk2(k)=lench(keylst(k)) 
  ENDDO  
  DO 3 i=1,nne 
     krtyp(i)=0 
     lk1=lench(keys(i)) 
     ipp=INDEX(keys(i)(1:lk1),'(') 
     IF(ipp.NE.0) THEN 
        pare=.true. 
        lk1=ipp-1 
     ELSE 
        pare=.false. 
     END IF
     DO 2 k=1,nkls 
        IF(pare) THEN 
           keywp=keylst(k) 
           ipp=INDEX(keywp,'(') 
           IF(ipp.EQ.0) GOTO 2 
           keywp=keywp(1:ipp-1) 
           IF(keys(i)(1:lk1).EQ.keywp) THEN 
              krtyp(i)=keytyp(k) 
              GOTO 3 
           END IF
        ELSE 
           IF(keys(i)(1:lk1).EQ.keylst(k)(1:lk2(k))) THEN 
              krtyp(i)=keytyp(k) 
              GOTO 3 
           END IF
        END IF
2    END DO
     lf=lench(namif(i)) 
     WRITE(*,100) keys(i)(1:lk1),krecnm(i),namif(i)(1:lf) 
100  FORMAT(' ERROR: unrecognized keyword "',a,'"'/&
     &        8x,'(record',i4,' in file "',a,'")')  
     fail=.true. 
3 END DO
  
  IF(fail) STOP '**** chkkey: abnormal end ***' 
END SUBROUTINE chkkey
! Copyright (C) 1997 by Mario Carpino (carpino@brera.mi.astro.it)       
! Version: March 5, 1997        
! --------------------------------------------------------------------- 
! 
!  *****************************************************************    
!  *   *    
!  *     R D N C H A       *    
!  *   *    
!  *Read a scalar quantity from the input namelist       *    
!  * (character string)    *    
!  *   *    
!  *****************************************************************    
! 
! INPUT:    CAT       -  Namelist keyword category  
! KEY       -  Namelist keyword name      
! REQRD     -  If true, when the keyword is not found in the  
!    namelist, an error message is issued and the   
!    FAIL flag is set 
! 
! OUTPUT:   V         -  Value associated with keyword        
! FOUND     -  True when the keyword has been found in the    
!    namelist         
! FAIL1     -  True if(.not.found.and.reqrd); false otherwise 
! FAIL      -  Set to true if(.not.found.and.reqrd);
!    otherwise, not changed     
! 
SUBROUTINE rdncha(cat,key,v,reqrd,found,fail1,fail) 
  USE option_input
  IMPLICIT NONE 
  CHARACTER*(*) cat,key 
  CHARACTER ck*(lckx),val*(lcvx),file*(lcfx),rest*(lcvx) 
  CHARACTER*50 vartyp 
  LOGICAL reqrd,found,fail,fail1,error 
  INTEGER ktyp0,ktyp1,lc,lk,lt,kr,lz,lf 
  
  CHARACTER*(*) v 
  
  INTEGER lench 
  EXTERNAL lench 
  
  vartyp='CHARACTER' 
  ktyp0=3 
  
  fail1=.false. 
! Find in namelist the value corresponding to the requested key         
  lc=lench(cat) 
  lk=lench(key) 
  lt=lc+lk 
  CALL chkpdf(lt,lckx,'lckx') 
  IF(lc.LE.0) THEN 
     ck=key(1:lk) 
  ELSE 
     ck=cat(1:lc)//key(1:lk) 
  END IF
  CALL getkv(ck(1:lt),val,ktyp1,file,kr,found) 
  lk=lench(key) 
  lz=lench(vartyp) 
  IF(.NOT.found) THEN 
     IF(reqrd) THEN 
        WRITE(*,100) ck(1:lt),vartyp(1:lz) 
        fail1=.true. 
        fail=.true. 
     END IF
     RETURN 
  END IF
  lf=lench(file) 
  IF(ktyp1.NE.ktyp0) THEN 
     WRITE(*,106) ck(1:lt),ktyp0,ktyp1 
     STOP '**** rdncha: abnormal end ****' 
  END IF
100 FORMAT(' ERROR: Missing definition of keyword "',a,'" (type: ',   &
     &       a,')')   
106 FORMAT(' ERROR: keyword type mismatch'/       &
     &       '        KEY = "',A,'"'/     &
     &       '        KTYP0 =',I2/        &
     &       '        KTYP1 =',I2)        
  
! Read value into output variable         
  CALL strcnt(val,v,rest,error) 
  IF(error) GOTO 1 
  
  RETURN 
  
! Error message       
1 CONTINUE 
  WRITE(*,101) ck(1:lt),vartyp(1:lz),kr,file(1:lf) 
101 FORMAT(' ERROR: Abnormal definition of keyword "',a,'" (type: ',  &
         &       a,')'/8x,'(record',i4,' in file "',a,'")')       
  fail1=.true. 
  fail=.true. 
  RETURN 
  
END SUBROUTINE rdncha
! Copyright (C) 1997 by Mario Carpino (carpino@brera.mi.astro.it)       
! Version: March 5, 1997        
! --------------------------------------------------------------------- 
! 
!  *****************************************************************    
!  *   *    
!  *     R D N I N T       *    
!  *   *    
!  *    Read an integer scalar quantity from the input namelist    *    
!  *   *    
!  *****************************************************************    
! 
! INPUT:    CAT       -  Namelist keyword category  
! KEY       -  Namelist keyword name      
! REQRD     -  If true, when the keyword is not found in the  
!    namelist, an error message is issued and the   
!    FAIL flag is set 
! 
! OUTPUT:   V         -  Value associated with keyword        
! FOUND     -  True when the keyword has been found in the    
!    namelist         
! FAIL1     -  True if(.not.found.and.reqrd); false otherwise 
! FAIL      -  Set to true if(.not.found.and.reqrd);
!    otherwise, not changed     
! 
SUBROUTINE rdnint(cat,key,v,reqrd,found,fail1,fail) 
  USE option_input
  IMPLICIT NONE 
  CHARACTER*(*) cat,key 
  CHARACTER ck*(lckx),val*(lcvx),file*(lcfx) 
  CHARACTER*50 vartyp 
  LOGICAL reqrd,found,fail,fail1 
  INTEGER ktyp0,ktyp1,lc,lk,lt,kr,lz,lf 
  INTEGER v 
  INTEGER lench 
  EXTERNAL lench 
  vartyp='INTEGER' 
  ktyp0=1 
  fail1=.false. 
! Find in namelist the value corresponding to the requested key         
  lc=lench(cat) 
  lk=lench(key) 
  lt=lc+lk 
  CALL chkpdf(lt,lckx,'lckx') 
  IF(lc.LE.0) THEN 
     ck=key(1:lk) 
  ELSE 
     ck=cat(1:lc)//key(1:lk) 
  END IF
  CALL getkv(ck(1:lt),val,ktyp1,file,kr,found) 
  lk=lench(key) 
  lz=lench(vartyp) 
  IF(.NOT.found) THEN 
     IF(reqrd) THEN 
        WRITE(*,100) ck(1:lt),vartyp(1:lz) 
        fail1=.true. 
        fail=.true. 
     END IF
     RETURN 
  END IF
  lf=lench(file) 
  IF(ktyp1.NE.ktyp0) THEN 
     WRITE(*,106) ck(1:lt),ktyp0,ktyp1 
     STOP '**** rdnxxx: abnormal end ****' 
  END IF
100 FORMAT(' ERROR: Missing definition of keyword "',a,'" (type: ',   &
     &       A,')')   
106 FORMAT(' ERROR: keyword type mismatch'/       &
     &       '        KEY = "',A,'"'/     &
     &       '        KTYP0 =',I2/        &
     &       '        KTYP1 =',I2)        
  ! Read value into output variable         
  READ(val,*,ERR=1) v 
  RETURN 
  ! Error message       
1 CONTINUE 
  WRITE(*,101) ck(1:lt),vartyp(1:lz),kr,file(1:lf) 
101 FORMAT(' ERROR: Abnormal definition of keyword "',a,'" (type: ',  &
     &       a,')'/8x,'(record',i4,' in file "',a,'")')       
  fail1=.true. 
  fail=.true. 
  
END SUBROUTINE rdnint
! Copyright (C) 1997 by Mario Carpino (carpino@brera.mi.astro.it)       
! Version: March 5, 1997        
! --------------------------------------------------------------------- 
! 
!  *****************************************************************    
!  *   *    
!  *     R D N L O G       *    
!  *   *    
!  *    Read a logical scalar quantity from the input namelist     *    
!  *   *    
!  *****************************************************************    
! 
! INPUT:    CAT       -  Namelist keyword category  
! KEY       -  Namelist keyword name      
! REQRD     -  If true, when the keyword is not found in the  
!    namelist, an error message is issued and the   
!    FAIL flag is set 
! 
! OUTPUT:   V         -  Value associated with keyword        
! FOUND     -  True when the keyword has been found in the    
!    namelist         
! FAIL1     -  True if(.not.found.and.reqrd); false otherwise 
! FAIL      -  Set to true if(.not.found.and.reqrd);
!    otherwise, not changed     
! 
SUBROUTINE rdnlog(cat,key,v,reqrd,found,fail1,fail) 
  USE option_input
  IMPLICIT NONE 
  CHARACTER*(*) cat,key 
  CHARACTER ck*(lckx),val*(lcvx),file*(lcfx) 
  CHARACTER*50 vartyp 
  LOGICAL reqrd,found,fail,fail1 
  INTEGER ktyp0,ktyp1,lc,lk,lt,kr,lz,lf 
  LOGICAL v 
  INTEGER lench 
  EXTERNAL lench 
  vartyp='LOGICAL' 
  ktyp0=4 
  fail1=.false. 
! Find in namelist the value corresponding to the requested key         
  lc=lench(cat) 
  lk=lench(key) 
  lt=lc+lk 
  CALL chkpdf(lt,lckx,'lckx') 
  IF(lc.LE.0) THEN 
     ck=key(1:lk) 
  ELSE 
     ck=cat(1:lc)//key(1:lk) 
  END IF
  CALL getkv(ck(1:lt),val,ktyp1,file,kr,found) 
  lk=lench(key) 
  lz=lench(vartyp) 
  IF(.NOT.found) THEN 
     IF(reqrd) THEN 
        WRITE(*,100) ck(1:lt),vartyp(1:lz) 
        fail1=.true. 
        fail=.true. 
     END IF
     RETURN 
  END IF
  lf=lench(file) 
  IF(ktyp1.NE.ktyp0) THEN 
     WRITE(*,106) ck(1:lt),ktyp0,ktyp1 
     STOP '**** rdnxxx: abnormal end ****' 
  END IF
100 FORMAT(' ERROR: Missing definition of keyword "',a,'" (type: ',   &
     &       a,')')   
106 FORMAT(' ERROR: keyword type mismatch'/       &
     &       '        KEY = "',A,'"'/     &
     &       '        KTYP0 =',I2/        &
     &       '        KTYP1 =',I2)        
  
! Read value into output variable   
  READ(val,*,ERR=1) v 
  
  RETURN 
! Error message       
1 CONTINUE 
  WRITE(*,101) ck(1:lt),vartyp(1:lz),kr,file(1:lf) 
101 FORMAT(' ERROR: Abnormal definition of keyword "',a,'" (type: ',  &
     &       a,')'/8x,'(record',i4,' in file "',a,'")')       
  fail1=.true. 
  fail=.true. 
  
END SUBROUTINE rdnlog
! Copyright (C) 1997 by Mario Carpino (carpino@brera.mi.astro.it)       
! Version: March 5, 1997        
! --------------------------------------------------------------------- 
! 
!  *****************************************************************    
!  *   *    
!  *     R D N R E A       *    
!  *   *    
!  *      Read a real scalar quantity from the input namelist      *    
!  *   *    
!  *****************************************************************    
! 
! INPUT:    CAT       -  Namelist keyword category  
! KEY       -  Namelist keyword name      
! REQRD     -  If true, when the keyword is not found in the  
!    namelist, an error message is issued and the   
!    FAIL flag is set 
! 
! OUTPUT:   V         -  Value associated with keyword        
! FOUND     -  True when the keyword has been found in the    
!    namelist         
! FAIL1     -  True if(.not.found.and.reqrd); false otherwise 
! FAIL      -  Set to true if(.not.found.and.reqrd);
!    otherwise, not changed     
! 
SUBROUTINE rdnrea(cat,key,v,reqrd,found,fail1,fail) 
  USE option_input
  IMPLICIT NONE 
  CHARACTER*(*) cat,key 
  CHARACTER ck*(lckx),val*(lcvx),file*(lcfx) 
  CHARACTER*50 vartyp 
  LOGICAL reqrd,found,fail,fail1 
  INTEGER ktyp0,ktyp1,lc,lk,lt,kr,lz,lf 
  DOUBLE PRECISION v 
  INTEGER lench 
  EXTERNAL lench 
  vartyp='REAL' 
  ktyp0=2 
  fail1=.false. 
! Find in namelist the value corresponding to the requested key         
  lc=lench(cat) 
  lk=lench(key) 
  lt=lc+lk 
  CALL chkpdf(lt,lckx,'lckx') 
  IF(lc.LE.0) THEN 
     ck=key(1:lk) 
  ELSE 
     ck=cat(1:lc)//key(1:lk) 
  END IF
  CALL getkv(ck(1:lt),val,ktyp1,file,kr,found) 
  lk=lench(key) 
  lz=lench(vartyp) 
  IF(.NOT.found) THEN 
     IF(reqrd) THEN 
        WRITE(*,100) ck(1:lt),vartyp(1:lz) 
        fail1=.true. 
        fail=.true. 
     END IF
     RETURN 
  END IF
  lf=lench(file) 
  IF(ktyp1.NE.ktyp0) THEN 
     WRITE(*,106) ck(1:lt),ktyp0,ktyp1 
     STOP '**** rdnxxx: abnormal end ****' 
  END IF
100 FORMAT(' ERROR: Missing definition of keyword "',a,'" (type: ',   &
     &       a,')')   
106 FORMAT(' ERROR: keyword type mismatch'/       &
     &       '        KEY = "',A,'"'/     &
     &       '        KTYP0 =',I2/        &
     &       '        KTYP1 =',I2)        
! Read value into output variable         
  READ(val,*,ERR=1) v 
  RETURN 
! Error message       
1 CONTINUE 
  WRITE(*,101) ck(1:lt),vartyp(1:lz),kr,file(1:lf) 
101 FORMAT(' ERROR: Abnormal definition of keyword "',a,'" (type: ',  &
     &       a,')'/8x,'(record',i4,' in file "',a,'")')       
  fail1=.true. 
  fail=.true. 
  
END SUBROUTINE rdnrea
! Copyright (C) 1997 by Mario Carpino (carpino@brera.mi.astro.it)       
! Version: March 5, 1997        
! --------------------------------------------------------------------- 
! 
!  *****************************************************************    
!  *   *    
!  *     R D N T I M       *    
!  *   *    
!  *   Read a time/date scalar quantity from the input namelist    *    
!  *   *    
!  *****************************************************************    
! 
! INPUT:    CAT       -  Namelist keyword category  
! KEY       -  Namelist keyword name      
! REQRD     -  If true, when the keyword is not found in the  
!    namelist, an error message is issued and the   
!    FAIL flag is set 
! 
! OUTPUT:   TIMSTR    -  Time value as character string       
! MJD       -  Modified Julian Date (integer part)  
! SEC       -  Seconds within the day     
! SCALE     -  Time scale (UTC/UT1/TAI/TDT/ET/GPS)  
! FOUND     -  True when the keyword has been found in the    
!    namelist         
! FAIL1     -  True if(.not.found.and.reqrd); false otherwise 
! FAIL      -  Set to true if(.not.found.and.reqrd);
!    otherwise, not changed     
! 
SUBROUTINE rdntim(cat,key,timstr,mjd,sec,scale,         &
     &        reqrd,found,fail1,fail)     
  USE option_input
  IMPLICIT NONE 
  CHARACTER*(*) cat,key,timstr,scale 
  CHARACTER ck*(lckx),val*(lcvx),file*(lcfx) 
  CHARACTER*50 vartyp 
  LOGICAL reqrd,found,fail,fail1,error 
  INTEGER ktyp0,ktyp1,lc,lk,lt,kr,lz,lf,mjd 
  DOUBLE PRECISION sec 
  INTEGER lench 
  EXTERNAL lench 
  vartyp='MJD' 
  ktyp0=5 
  fail1=.false. 
! Find in namelist the value corresponding to the requested key         
  lc=lench(cat) 
  lk=lench(key) 
  lt=lc+lk 
  CALL chkpdf(lt,lckx,'lckx') 
  IF(lc.LE.0) THEN 
     ck=key(1:lk) 
  ELSE 
     ck=cat(1:lc)//key(1:lk) 
  END IF
  CALL getkv(ck(1:lt),val,ktyp1,file,kr,found) 
  lk=lench(key) 
  lz=lench(vartyp) 
  IF(.NOT.found) THEN 
     IF(reqrd) THEN 
        WRITE(*,100) ck(1:lt),vartyp(1:lz) 
        fail1=.true. 
        fail=.true. 
     END IF
     RETURN 
  END IF
  IF(ktyp1.NE.ktyp0) THEN 
     WRITE(*,106) ck(1:lt),ktyp0,ktyp1 
     STOP '**** rdnxxx: abnormal end ****' 
  END IF
100 FORMAT(' ERROR: Missing definition of keyword "',a,'" (type: ',   &
     &       a,')')   
106 FORMAT(' ERROR: keyword type mismatch'/       &
     &       '        KEY = "',A,'"'/     &
     &       '        KTYP0 =',I2/        &
     &       '        KTYP1 =',I2)        
! Read value into output variables        
  CALL ch2tim(val,mjd,sec,scale,error) 
  IF(error) GOTO 1 
  timstr=val 
  RETURN 
! Error message       
1 CONTINUE 
  lf=lench(file) 
  WRITE(*,101) ck(1:lt),vartyp(1:lz),kr,file(1:lf) 
101 FORMAT(' ERROR: Abnormal definition of keyword "',a,'" (type: ',  &
     &       a,')'/8x,'(record',i4,' in file "',a,'")')       
  fail1=.true. 
  fail=.true. 
  
END SUBROUTINE rdntim
! Copyright (C) 1997 by Mario Carpino (carpino@brera.mi.astro.it)       
! Version: March 7, 1997        
! --------------------------------------------------------------------- 
! 
!  *****************************************************************    
!  *   *    
!  *     R D N R E F       *    
!  *   *    
!  *Read a scalar quantity from the input namelist       *    
!  *     (reference system description)        *    
!  *   *    
!  *****************************************************************    
! 
! INPUT:    CAT       -  Namelist keyword category  
! KEY       -  Namelist keyword name      
! REQRD     -  If true, when the keyword is not found in the  
!    namelist, an error message is issued and the   
!    FAIL flag is set 
! 
! OUTPUT:   RSYS      -  Reference system type (EQUM/EQUT/ECLM)         
! EPOCH     -  Epoch specification (J2000/OFDATE)   
! FOUND     -  True when the keyword has been found in the    
!    namelist         
! FAIL1     -  True if(.not.found.and.reqrd); false otherwise 
! FAIL      -  Set to true if(.not.found.and.reqrd);
!    otherwise, not changed     
! 
SUBROUTINE rdnref(cat,key,rsys,epoch,reqrd,found,fail1,fail) 
  USE reference_systems 
  USE option_input
  IMPLICIT NONE 
  CHARACTER*(*) cat,key 
  CHARACTER ck*(lckx),val*(lcvx),file*(lcfx) 
  CHARACTER*50 vartyp 
  LOGICAL reqrd,found,fail,fail1,error 
  INTEGER ktyp0,ktyp1,lc,lk,lt,kr,lz,lf 
  CHARACTER*(*) rsys,epoch 
  INTEGER lench 
  EXTERNAL lench 
  vartyp='REFERENCE SYSTEM' 
  ktyp0=6 
  fail1=.false. 
! Find in namelist the value corresponding to the requested key         
  lc=lench(cat) 
  lk=lench(key) 
  lt=lc+lk 
  CALL chkpdf(lt,lckx,'lckx') 
  IF(lc.LE.0) THEN 
     ck=key(1:lk) 
  ELSE 
     ck=cat(1:lc)//key(1:lk) 
  END IF
  CALL getkv(ck(1:lt),val,ktyp1,file,kr,found) 
  lk=lench(key) 
  lz=lench(vartyp) 
  IF(.NOT.found) THEN 
     IF(reqrd) THEN 
        WRITE(*,100) ck(1:lt),vartyp(1:lz) 
        fail1=.true. 
        fail=.true. 
     END IF
     rsys=' ' 
     epoch=' ' 
     RETURN 
  END IF
  IF(ktyp1.NE.ktyp0) THEN 
     WRITE(*,106) ck(1:lt),ktyp0,ktyp1 
     STOP '**** rdnref: abnormal end ****' 
  END IF
100 FORMAT(' ERROR: Missing definition of keyword "',a,'" (type: ',   &
     &       a,')')   
106 FORMAT(' ERROR: keyword type mismatch'/       &
     &       '        KEY = "',A,'"'/     &
     &       '        KTYP0 =',I2/        &
     &       '        KTYP1 =',I2)        
! Read value into output variable         
  CALL ch2ref(val,rsys,epoch,error) 
  IF(error) GOTO 1 
  RETURN 
! Error message       
1 CONTINUE 
  lf=lench(file) 
  WRITE(*,101) ck(1:lt),vartyp(1:lz),kr,file(1:lf) 
101 FORMAT(' ERROR: Abnormal definition of keyword "',a,'" (type: ',  &
     &       a,')'/8x,'(record',i4,' in file "',a,'")')       
  fail1=.true. 
  fail=.true. 
  RETURN   
END SUBROUTINE rdnref
! Copyright (C) 1997 by Mario Carpino (carpino@brera.mi.astro.it)       
! Version: March 5, 1997        
! --------------------------------------------------------------------- 
! 
!  *****************************************************************    
!  *   *    
!  *     R D V I N T       *    
!  *   *    
!  *    Read an integer vector quantity from the input namelist    *    
!  *   *    
!  *****************************************************************    
! 
! INPUT:    CAT       -  Namelist keyword category  
! KEY       -  Namelist keyword name      
! N         -  Dimension of the vector to be read   
! REQRD     -  If true, when the keyword is not found in the  
!    namelist, an error message is issued and the   
!    FAIL flag is set 
! 
! OUTPUT:   V         -  Value associated with keyword        
! FOUND     -  True when the keyword has been found in the    
!    namelist         
! FAIL1     -  True if(.not.found.and.reqrd); false otherwise 
! FAIL      -  Set to true if(.not.found.and.reqrd);
!    otherwise, not changed     
! 
      SUBROUTINE rdvint(cat,key,v,n,reqrd,found,fail1,fail) 
      USE option_input
      IMPLICIT NONE 
      CHARACTER*(*) cat,key 
      CHARACTER ck*(lckx),val*(lcvx),file*(lcfx) 
      CHARACTER*50 vartyp 
      LOGICAL reqrd,found,fail,fail1 
  
  
      INTEGER n,v(n) 
      INTEGER ktyp0,ktyp1,lc,lk,lt,kr,lz,lf,i 
  
      INTEGER lench 
      EXTERNAL lench 
  
      vartyp='INTEGER VECTOR' 
      ktyp0=1 
  
      fail1=.false. 
! Find in namelist the value corresponding to the requested key         
      lc=lench(cat) 
      lk=lench(key) 
      lt=lc+lk 
      CALL chkpdf(lt,lckx,'lckx') 
      IF(lc.LE.0) THEN 
ck=key(1:lk) 
      ELSE 
ck=cat(1:lc)//key(1:lk) 
      END IF 
      CALL getkv(ck(1:lt),val,ktyp1,file,kr,found) 
      lk=lench(key) 
      lz=lench(vartyp) 
      IF(.NOT.found) THEN 
IF(reqrd) THEN 
    WRITE(*,100) ck(1:lt),vartyp(1:lz) 
    fail1=.true. 
    fail=.true. 
END IF 
RETURN 
      END IF 
      IF(ktyp1.NE.ktyp0) THEN 
WRITE(*,106) ck(1:lt),ktyp0,ktyp1 
STOP '**** rdnint: abnormal end ****' 
      END IF 
  100 FORMAT(' ERROR: Missing definition of keyword "',a,'" (type: ',   &
     &       a,')')   
  106 FORMAT(' ERROR: keyword type mismatch'/       &
     &       '        KEY = "',A,'"'/     &
     &       '        KTYP0 =',I2/        &
     &       '        KTYP1 =',I2)        
  
! Read value into output variable         
  
      READ(val,*,ERR=1) (v(i),i=1,n) 
  
      RETURN 
  
! Error message       
    1 CONTINUE 
      lf=lench(file) 
      WRITE(*,101) ck(1:lt),vartyp(1:lz),kr,file(1:lf) 
  101 FORMAT(' ERROR: Abnormal definition of keyword "',a,'" (type: ',  &
     &       a,')'/8x,'(record',i4,' in file "',a,'")')       
      fail1=.true. 
      fail=.true. 
      END   
! Copyright (C) 1997 by Mario Carpino (carpino@brera.mi.astro.it)       
! Version: March 5, 1997        
! --------------------------------------------------------------------- 
! 
!  *****************************************************************    
!  *   *    
!  *     R D V R E A       *    
!  *   *    
!  *      Read a real vector quantity from the input namelist      *    
!  *   *    
!  *****************************************************************    
! 
! INPUT:    CAT       -  Namelist keyword category  
! KEY       -  Namelist keyword name      
! N         -  Dimension of the vector to be read   
! REQRD     -  If true, when the keyword is not found in the  
!    namelist, an error message is issued and the   
!    FAIL flag is set 
! 
! OUTPUT:   V         -  Value associated with keyword        
! FOUND     -  True when the keyword has been found in the    
!    namelist         
! FAIL1     -  True if(.not.found.and.reqrd); false otherwise 
! FAIL      -  Set to true if(.not.found.and.reqrd);
!    otherwise, not changed     
! 
      SUBROUTINE rdvrea(cat,key,v,n,reqrd,found,fail1,fail) 
      USE option_input
      IMPLICIT NONE 
      CHARACTER*(*) cat,key 
      CHARACTER ck*(lckx),val*(lcvx),file*(lcfx) 
      CHARACTER*50 vartyp 
      LOGICAL reqrd,found,fail,fail1 
  
  
      INTEGER n 
      DOUBLE PRECISION v(n) 
      INTEGER ktyp0,ktyp1,lc,lk,lt,kr,lz,lf,i 
  
      INTEGER lench 
      EXTERNAL lench 
  
      vartyp='REAL VECTOR' 
      ktyp0=2 
  
      fail1=.false. 
! Find in namelist the value corresponding to the requested key         
      lc=lench(cat) 
      lk=lench(key) 
      lt=lc+lk 
      CALL chkpdf(lt,lckx,'lckx') 
      IF(lc.LE.0) THEN 
ck=key(1:lk) 
      ELSE 
ck=cat(1:lc)//key(1:lk) 
      END IF 
      CALL getkv(ck(1:lt),val,ktyp1,file,kr,found) 
      lk=lench(key) 
      lz=lench(vartyp) 
      IF(.NOT.found) THEN 
IF(reqrd) THEN 
    WRITE(*,100) ck(1:lt),vartyp(1:lz) 
    fail1=.true. 
    fail=.true. 
END IF 
RETURN 
      END IF 
      IF(ktyp1.NE.ktyp0) THEN 
WRITE(*,106) ck(1:lt),ktyp0,ktyp1 
STOP '**** rdnrea: abnormal end ****' 
      END IF 
  100 FORMAT(' ERROR: Missing definition of keyword "',a,'" (type: ',   &
     &       a,')')   
  106 FORMAT(' ERROR: keyword type mismatch'/       &
     &       '        KEY = "',A,'"'/     &
     &       '        KTYP0 =',I2/        &
     &       '        KTYP1 =',I2)        
  
! Read value into output variable         
  
      READ(val,*,ERR=1) (v(i),i=1,n) 
  
      RETURN 
  
! Error message       
    1 CONTINUE 
      lf=lench(file) 
      WRITE(*,101) ck(1:lt),vartyp(1:lz),kr,file(1:lf) 
  101 FORMAT(' ERROR: Abnormal definition of keyword "',a,'" (type: ',  &
     &       a,')'/8x,'(record',i4,' in file "',a,'")')       
      fail1=.true. 
      fail=.true. 
      END   
! Copyright (C) 1997 by Mario Carpino (carpino@brera.mi.astro.it)       
! Version: March 5, 1997        
! --------------------------------------------------------------------- 
! 
!  *****************************************************************    
!  *   *    
!  *     R D V C H A       *    
!  *   *    
!  *   Read a character vector quantity from the input namelist    *    
!  *   *    
!  *****************************************************************    
! 
! INPUT:    CAT       -  Namelist keyword category  
! KEY       -  Namelist keyword name      
! N         -  Dimension of the vector to be read   
! REQRD     -  If true, when the keyword is not found in the  
!    namelist, an error message is issued and the   
!    FAIL flag is set 
! 
! OUTPUT:   V         -  Value associated with keyword        
! FOUND     -  True when the keyword has been found in the    
!    namelist         
! FAIL1     -  True if(.not.found.and.reqrd); false otherwise 
! FAIL      -  Set to true if(.not.found.and.reqrd);
!    otherwise, not changed     
! 
      SUBROUTINE rdvcha(cat,key,v,n,reqrd,found,fail1,fail) 
      USE option_input
      IMPLICIT NONE 
      CHARACTER*(*) cat,key 
      CHARACTER ck*(lckx),val*(lcvx),file*(lcfx) 
      CHARACTER*50 vartyp 
      LOGICAL reqrd,found,fail,fail1 
  
      INTEGER n 
      CHARACTER*(*) v(n) 
      CHARACTER*(lcvx) rest 
      INTEGER ktyp0,ktyp1,lc,lk,lt,kr,lz,lf,i 
      LOGICAL error 
  
      INTEGER lench 
      EXTERNAL lench 
  
      vartyp='CHARACTER VECTOR' 
      ktyp0=3 
  
      fail1=.false. 
! Find in namelist the value corresponding to the requested key         
      lc=lench(cat) 
      lk=lench(key) 
      lt=lc+lk 
      CALL chkpdf(lt,lckx,'lckx') 
      IF(lc.LE.0) THEN 
ck=key(1:lk) 
      ELSE 
ck=cat(1:lc)//key(1:lk) 
      END IF 
      CALL getkv(ck(1:lt),val,ktyp1,file,kr,found) 
      lk=lench(key) 
      lz=lench(vartyp) 
      IF(.NOT.found) THEN 
IF(reqrd) THEN 
    WRITE(*,100) ck(1:lt),vartyp(1:lz) 
    fail1=.true. 
    fail=.true. 
END IF 
RETURN 
      END IF 
      IF(ktyp1.NE.ktyp0) THEN 
WRITE(*,106) ck(1:lt),ktyp0,ktyp1 
STOP '**** rdncha: abnormal end ****' 
      END IF 
  100 FORMAT(' ERROR: Missing definition of keyword "',a,'" (type: ',   &
     &       a,')')   
  106 FORMAT(' ERROR: keyword type mismatch'/       &
     &       '        KEY = "',A,'"'/     &
     &       '        KTYP0 =',I2/        &
     &       '        KTYP1 =',I2)        
  
! Read value into output variable         
      DO 2 i=1,n 
      CALL strcnt(val,v(i),rest,error) 
      IF(error) GOTO 1 
      val=rest 
    2 END DO 
      RETURN 
  
! Error message       
    1 CONTINUE 
      lf=lench(file) 
      WRITE(*,101) ck(1:lt),vartyp(1:lz),kr,file(1:lf) 
  101 FORMAT(' ERROR: Abnormal definition of keyword "',a,'" (type: ',  &
     &       a,')'/8x,'(record',i4,' in file "',a,'")')       
      fail1=.true. 
      fail=.true. 
      END   
! Copyright (C) 1997 by Mario Carpino (carpino@brera.mi.astro.it)       
! Version: March 5, 1997        
! --------------------------------------------------------------------- 
! 
!  *****************************************************************    
!  *   *    
!  *     R D M I N T       *    
!  *   *    
!  *    Read an integer vector quantity from the input namelist    *    
!  *  inferring the dimension from the number of values supplied   *    
!  *   *    
!  *****************************************************************    
! 
! INPUT:    CAT       -  Namelist keyword category  
! KEY       -  Namelist keyword name      
! NX        -  Dimension of the vector to be read   
! NAMNX     -  Real name of the NX variable         
! REQRD     -  If true, when the keyword is not found in the  
!    namelist, an error message is issued and the   
!    FAIL flag is set 
! 
! OUTPUT:   V         -  Value associated with keyword        
! N         -  Actual number of elements read in the vector   
! FOUND     -  True when the keyword has been found in the    
!    namelist         
! FAIL1     -  True if(.not.found.and.reqrd); false otherwise 
! FAIL      -  Set to true if(.not.found.and.reqrd);
!    otherwise, not changed     
! 
      SUBROUTINE rdmint(cat,key,v,n,nx,namnx,reqrd,found,fail1,fail) 
      USE option_input
      IMPLICIT NONE 
      CHARACTER*(*) cat,key 
      CHARACTER ck*(lckx),val*(lcvx),file*(lcfx) 
      CHARACTER*50 vartyp 
      LOGICAL reqrd,found,fail,fail1 
      INTEGER n,nx,ktyp0,ktyp1,lc,lk,lt,kr,lz,lf,i 
  
  
      CHARACTER*(*) namnx 
      INTEGER v(nx) 
  
      INTEGER lench,nitchs 
      EXTERNAL lench,nitchs 
  
      vartyp='INTEGER VECTOR' 
      ktyp0=1 
      n=0 
  
      fail1=.false. 
! Find in namelist the value corresponding to the requested key         
      lc=lench(cat) 
      lk=lench(key) 
      lt=lc+lk 
      CALL chkpdf(lt,lckx,'lckx') 
      IF(lc.LE.0) THEN 
ck=key(1:lk) 
      ELSE 
ck=cat(1:lc)//key(1:lk) 
      END IF 
      CALL getkv(ck(1:lt),val,ktyp1,file,kr,found) 
      lk=lench(key) 
      lz=lench(vartyp) 
      IF(.NOT.found) THEN 
IF(reqrd) THEN 
    WRITE(*,100) ck(1:lt),vartyp(1:lz) 
    fail1=.true. 
    fail=.true. 
END IF 
RETURN 
      END IF 
      IF(ktyp1.NE.ktyp0) THEN 
WRITE(*,106) ck(1:lt),ktyp0,ktyp1 
STOP '**** rdnint: abnormal end ****' 
      END IF 
  100 FORMAT(' ERROR: Missing definition of keyword "',a,'" (type: ',   &
     &       a,')')   
  106 FORMAT(' ERROR: keyword type mismatch'/       &
     &       '        KEY = "',A,'"'/     &
     &       '        KTYP0 =',I2/        &
     &       '        KTYP1 =',I2)        
  
! Read value into output variable         
      n=nitchs(val) 
      CALL chkpdf(n,nx,namnx) 
      READ(val,*,ERR=1) (v(i),i=1,n) 
  
      RETURN 
  
! Error message       
    1 CONTINUE 
      lf=lench(file) 
      WRITE(*,101) ck(1:lt),vartyp(1:lz),kr,file(1:lf) 
  101 FORMAT(' ERROR: Abnormal definition of keyword "',a,'" (type: ',  &
     &       a,')'/8x,'(record',i4,' in file "',a,'")')       
      fail1=.true. 
      fail=.true. 
  
      END   
! Copyright (C) 1997 by Mario Carpino (carpino@brera.mi.astro.it)       
! Version: March 5, 1997        
! --------------------------------------------------------------------- 
! 
!  *****************************************************************    
!  *   *    
!  *     R D M R E A       *    
!  *   *    
!  *      Read a real vector quantity from the input namelist      *    
!  *  inferring the dimension from the number of values supplied   *    
!  *   *    
!  *****************************************************************    
! 
! INPUT:    CAT       -  Namelist keyword category  
! KEY       -  Namelist keyword name      
! NX        -  Dimension of the vector to be read   
! NAMNX     -  Real name of the NX variable         
! REQRD     -  If true, when the keyword is not found in the  
!    namelist, an error message is issued and the   
!    FAIL flag is set 
! 
! OUTPUT:   V         -  Value associated with keyword        
! N         -  Actual number of elements read in the vector   
! FOUND     -  True when the keyword has been found in the    
!    namelist         
! FAIL1     -  True if(.not.found.and.reqrd); false otherwise 
! FAIL      -  Set to true if(.not.found.and.reqrd);
!    otherwise, not changed     
! 
      SUBROUTINE rdmrea(cat,key,v,n,nx,namnx,reqrd,found,fail1,fail) 
      USE option_input
      IMPLICIT NONE 
      CHARACTER*(*) cat,key 
      CHARACTER ck*(lckx),val*(lcvx),file*(lcfx) 
      CHARACTER*50 vartyp 
      LOGICAL reqrd,found,fail,fail1 
      INTEGER n,nx,ktyp0,ktyp1,lc,lk,lt,kr,lz,lf,i 
  
  
      CHARACTER*(*) namnx 
      DOUBLE PRECISION v(nx) 
  
      INTEGER lench,nitchs 
      EXTERNAL lench,nitchs 
  
      vartyp='REAL VECTOR' 
      ktyp0=2 
      n=0 
  
      fail1=.false. 
! Find in namelist the value corresponding to the requested key         
      lc=lench(cat) 
      lk=lench(key) 
      lt=lc+lk 
      CALL chkpdf(lt,lckx,'lckx') 
      IF(lc.LE.0) THEN 
ck=key(1:lk) 
      ELSE 
ck=cat(1:lc)//key(1:lk) 
      END IF 
      CALL getkv(ck(1:lt),val,ktyp1,file,kr,found) 
      lk=lench(key) 
      lz=lench(vartyp) 
      IF(.NOT.found) THEN 
IF(reqrd) THEN 
    WRITE(*,100) ck(1:lt),vartyp(1:lz) 
    fail1=.true. 
    fail=.true. 
END IF 
RETURN 
      END IF 
      IF(ktyp1.NE.ktyp0) THEN 
WRITE(*,106) ck(1:lt),ktyp0,ktyp1 
STOP '**** rdnrea: abnormal end ****' 
      END IF 
  100 FORMAT(' ERROR: Missing definition of keyword "',a,'" (type: ',   &
     &       a,')')   
  106 FORMAT(' ERROR: keyword type mismatch'/       &
     &       '        KEY = "',A,'"'/     &
     &       '        KTYP0 =',I2/        &
     &       '        KTYP1 =',I2)        
  
! Read value into output variable         
      n=nitchs(val) 
      CALL chkpdf(n,nx,namnx) 
  
      READ(val,*,ERR=1) (v(i),i=1,n) 
  
      RETURN 
  
! Error message       
    1 CONTINUE 
      lf=lench(file) 
      WRITE(*,101) ck(1:lt),vartyp(1:lz),kr,file(1:lf) 
  101 FORMAT(' ERROR: Abnormal definition of keyword "',a,'" (type: ',  &
     &       a,')'/8x,'(record',i4,' in file "',a,'")')       
      fail1=.true. 
      fail=.true. 
  
      END   
! Copyright (C) 1997 by Mario Carpino (carpino@brera.mi.astro.it)       
! Version: March 5, 1997        
! --------------------------------------------------------------------- 
! 
!  *****************************************************************    
!  *   *    
!  *     R D M C H A       *    
!  *   *    
!  *    Read a character vector quantity from the input namelist   *    
!  *  inferring the dimension from the number of values supplied   *    
!  *   *    
!  *****************************************************************    
! 
! INPUT:    CAT       -  Namelist keyword category  
! KEY       -  Namelist keyword name      
! NX        -  Dimension of the vector to be read   
! NAMNX     -  Real name of the NX variable         
! REQRD     -  If true, when the keyword is not found in the  
!    namelist, an error message is issued and the   
!    FAIL flag is set 
! 
! OUTPUT:   V         -  Value associated with keyword        
! N         -  Actual number of elements read in the vector   
! FOUND     -  True when the keyword has been found in the    
!    namelist         
! FAIL1     -  True if(.not.found.and.reqrd); false otherwise 
! FAIL      -  Set to true if(.not.found.and.reqrd);
!    otherwise, not changed     
! 
      SUBROUTINE rdmcha(cat,key,v,n,nx,namnx,reqrd,found,fail1,fail) 
      USE option_input
      IMPLICIT NONE 
      CHARACTER*(*) cat,key 
      CHARACTER ck*(lckx),val*(lcvx),file*(lcfx),rest*(lcvx) 
      CHARACTER*50 vartyp 
      LOGICAL reqrd,found,fail,fail1,error 
      INTEGER n,nx,ktyp0,ktyp1,lc,lk,lt,kr,lz,lf,i 
  
      CHARACTER*(*) namnx 
      CHARACTER*(*) v(nx) 
  
      INTEGER lench,nitchs 
      EXTERNAL lench,nitchs 
  
      vartyp='CHARACTER VECTOR' 
      ktyp0=3 
      n=0 
  
      fail1=.false. 
! Find in namelist the value corresponding to the requested key         
      lc=lench(cat) 
      lk=lench(key) 
      lt=lc+lk 
      CALL chkpdf(lt,lckx,'lckx') 
      IF(lc.LE.0) THEN 
ck=key(1:lk) 
      ELSE 
ck=cat(1:lc)//key(1:lk) 
      END IF 
      CALL getkv(ck(1:lt),val,ktyp1,file,kr,found) 
      lk=lench(key) 
      lz=lench(vartyp) 
      IF(.NOT.found) THEN 
IF(reqrd) THEN 
    WRITE(*,100) ck(1:lt),vartyp(1:lz) 
    fail1=.true. 
    fail=.true. 
END IF 
RETURN 
      END IF 
      IF(ktyp1.NE.ktyp0) THEN 
WRITE(*,106) ck(1:lt),ktyp0,ktyp1 
STOP '**** rdncha: abnormal end ****' 
      END IF 
  100 FORMAT(' ERROR: Missing definition of keyword "',a,'" (type: ',   &
     &       a,')')   
  106 FORMAT(' ERROR: keyword type mismatch'/       &
     &       '        KEY = "',A,'"'/     &
     &       '        KTYP0 =',I2/        &
     &       '        KTYP1 =',I2)        
  
! Read value into output variable         
      n=nitchs(val) 
      CALL chkpdf(n,nx,namnx) 
      DO 2 i=1,n 
      CALL strcnt(val,v(i),rest,error) 
      IF(error) GOTO 1 
      val=rest 
    2 END DO 
      RETURN 
  
! Error message       
    1 CONTINUE 
      lf=lench(file) 
      WRITE(*,101) ck(1:lt),vartyp(1:lz),kr,file(1:lf) 
  101 FORMAT(' ERROR: Abnormal definition of keyword "',a,'" (type: ',  &
     &       a,')'/8x,'(record',i4,' in file "',a,'")')       
      fail1=.true. 
      fail=.true. 
  
      END   
! Copyright (C) 1996 by Mario Carpino (carpino@brera.mi.astro.it)       
! Version: September 19, 1996   
! --------------------------------------------------------------------- 
! 
!  *****************************************************************    
!  *   *    
!  *      G E T K V        *    
!  *   *    
!  *Locates a key/value pair in the input namelist       *    
!  *   *    
!  *****************************************************************    
! 
! INPUT:    KEY       -  Required key     
! 
! OUTPUT:   VAL       -  Value field associated with requested key      
! KTYP      -  Keyword type:    
!1 - integer
!2 - real   
!3 - character string 
!4 - logical
!5 - MJD    
! FILE      -  Name of the file where the key was read        
! KR        -  Record num. in the file where the key was read 
! FOUND     -  If .false., the key is not defined   
  
      SUBROUTINE getkv(key,val,ktyp,file,kr,found) 
      USE option_input
      IMPLICIT NONE 
      CHARACTER*(*) key,val,file 
      INTEGER ktyp,kr,lv,lf,lk,i 
      LOGICAL found 
  
      INTEGER lench 
      EXTERNAL lench 
  
      IF(iicnam.NE.36) STOP '**** getkv: internal error (01) ****' 
  
      lv=LEN(val) 
      lf=LEN(file) 
      CALL chkpdf(lv,lcvx,'lcvx') 
      CALL chkpdf(lf,lcfx,'lcfx') 
  
      lk=lench(key) 
      DO 1 i=1,nne 
      IF(key(1:lk).EQ.keys(i)) THEN 
val=vals(i) 
file=namif(i) 
kr=krecnm(i) 
kuorl=kuorl+1 
kuord(i)=kuorl 
ktyp=krtyp(i) 
IF(ktyp.LE.0.OR.ktyp.GT.5) THEN 
    WRITE(*,100) key(1:lk),ktyp 
    STOP '**** getkv: internal error (02) ****' 
END IF 
found=.true. 
RETURN 
      END IF 
    1 END DO 
  100 FORMAT(' ERROR: illegal keyword type'/        &
     &       '        KEY="',A,'"'/       &
     &       '        TYPE=',I2)
  
      found=.false. 
  
      END   
! Copyright (C) 1996-1998 by Mario Carpino (carpino@brera.mi.astro.it)  
! Version: January 12, 1998     
! --------------------------------------------------------------------- 
! 
!  *****************************************************************    
!  *   *    
!  *     C H K P D F       *    
!  *   *    
!  *       Check parameter definition*    
!  *   *    
!  *****************************************************************    
! 
! INPUT:    NA        -  Actual value required for the parameter        
! NX        -  Value declared for the parameter     
! NAME      -  Parameter name   
! 
SUBROUTINE chkpdf(na,nx,name) 
  IMPLICIT NONE 
  INTEGER na,nx,ll 
  CHARACTER*(*) name 
  INTEGER lench 
  EXTERNAL lench 
  IF(na.LE.nx) RETURN 
  ll=lench(name) 
  WRITE(*,100) name(1:ll),na 
100 FORMAT(' **** Insufficient PARAMETER definition ****'/  &
     &       ' **** Please set ',a,' >=',i7,' ****')
  STOP '**** chkpdf: abnormal end ****' 
END SUBROUTINE chkpdf

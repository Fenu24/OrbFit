! ===================MODULE header_input============================    
! HEADER NAMELIST: Carpino 1996-1998           
! CONTAINS:           
! SUBROUTINES         
! rdfnam		read and store header namelist       
! rdfcha		read a character value               
! rdfint		read an integer value                
! rdflog		read a logical value                 
! rdfrea		read a real value                    
! rdftim		read a date/time value               
! rdfref		read a reference system description  
! *splkvc		split a record into keyword+value+comment 
! *chkfln		check field length for input records
! *getrsc		get a data record skipping comment lines 
! sv2int                used by ephem_prop, iers-ser                    
! initopt               common part of all option input routines
! input_cha_opt etc     input character option etc
!                     
!                     
!                     
! Copyright (C) 1997 by Mario Carpino (carpino@brera.mi.astro.it)       
! Version: March 8, 1997                       
! --------------------------------------------------------------------- 
!                     
!  *****************************************************************    
!  *                                                               *    
!  *                         R D F N A M                           *    
!  *                                                               *    
!  *            Read simplified file-header namelist               *    
!  *                                                               *    
!  *****************************************************************    
!                     
! INPUT:    UNIT      -  Input unit            
!           FILE      -  Input filename (for error messages)            
!                     
! OUTPUT:   NR        -  Records read so far   
!                     
      SUBROUTINE rdfnam(unit,file,nr) 
      USE char_str
      USE option_input
      IMPLICIT NONE 
      INCLUDE 'parcmc.h90' 
      INTEGER unit,nr 
      INTEGER lf,lk,i 
      CHARACTER*(lchx) rec,key1,val1,comm 
      CHARACTER*(*) file 
      LOGICAL skip,end 
                      
      INTEGER lench 
      EXTERNAL lench 
                      
      iicfnm=0 
                      
      nfne=0 
      kuorlf=0 
      nmif=file 
      hnfuni=unit 
      lf=lench(file) 
      nr=0 
                      
! Read a record from the namelist              
    2 READ(unit,100,end=11) rec 
  100 FORMAT(a) 
      nr=nr+1 
      CALL splkvc(rec,key1,val1,comm,skip,end) 
                      
! Detection of the end of the namelist         
      IF(end) THEN 
          iicfnm=36 
          RETURN 
      END IF 
                      
! Skip comment lines  
      IF(skip) GOTO 2 
                      
! Look if the key is already present in the namelist                    
      lk=lench(key1) 
      DO 4 i=1,nfne 
      IF(key1(1:lk).EQ.keysf(i)) THEN 
          WRITE(*,210) key1(1:lk),file(1:lf),krcfnm(i),nr 
          STOP '**** rdfnam: abnormal end ****' 
      END IF 
    4 END DO 
  210 FORMAT(' rdfnam: duplicate definition of keyword "',a,            &
     &       '" in file "',a,'":'/             &
     &       '         first  occurence at line',i5/                    &
     &       '         second occurence at line',i5)                    
                      
      nfne=nfne+1 
      CALL chkpdf(nfne,nfnex,'nfnex') 
      i=nfne 
      keysf(i)=key1 
      valsf(i)=val1 
      krcfnm(i)=nr 
      kuorf(i)=0 
      GOTO 2 
                      
   11 CONTINUE 
                      
! Abort when no "END_OF_HEADER" record is found
!     WRITE(*,200)'end',file(1:lf),nr          
!     STOP '**** rdfnam: abnormal end ****'    
!200  FORMAT(' rdfnam: cannot find namelist ',a,':'/                    
!    +       '         unexpected END of file "',a,'" after record',i6) 
                      
      iicfnm=36 
                      
    END SUBROUTINE rdfnam
! Copyright (C) 1997 by Mario Carpino (carpino@brera.mi.astro.it)       
! Version: March 8, 1997                       
! --------------------------------------------------------------------- 
!                     
!  *****************************************************************    
!  *                                                               *    
!  *                         R D F C H A                           *    
!  *                                                               *    
!  *                  Read a character string                      *    
!  *           from simplified file-header namelist                *    
!  *                                                               *    
!  *****************************************************************    
!                     
! INPUT:    UNIT      -  Input unit            
!           KEY       -  Keyword               
!           REQRD     -  If true, when the keyword is not found in the  
!                        namelist, execution stops with an error message
!                     
! OUTPUT:   C         -  Character string      
!           FOUND     -  True when the keyword has been found in the    
!                        namelist              
!           KR        -  Record number in the input file                
!                     
      SUBROUTINE rdfcha(unit,key,reqrd,c,found,kr)
      USE option_input 
      IMPLICIT NONE 
      CHARACTER*(*) key,c 
      LOGICAL reqrd,found 
      INTEGER kr,unit 
                      
      INTEGER i,lk,lz,lf 
      CHARACTER vartyp*30,rest*(lcvx) 
      LOGICAL error 
                      
      INTEGER lench 
      EXTERNAL lench 
                      
      IF(iicfnm.NE.36) STOP '**** rdfcha: internal error (01) ****' 
      IF(unit.NE.hnfuni) STOP '**** rdfcha: internal error (02) ****' 
                      
      vartyp='CHARACTER' 
                      
      kr=0 
      DO 1 i=1,nfne 
      IF(key.EQ.keysf(i)) THEN 
          found=.true. 
          CALL strcnt(valsf(i),c,rest,error) 
          IF(error) GOTO 10 
          kr=krcfnm(i) 
          kuorlf=kuorlf+1 
          kuorf(i)=kuorlf 
          RETURN 
      END IF 
    1 END DO 
                      
      found=.false. 
      IF(reqrd) THEN 
          lk=lench(key) 
          lz=lench(vartyp) 
          lf=lench(nmif) 
          WRITE(*,100) key(1:lk),vartyp(1:lz),nmif(1:lf) 
          STOP '**** rdfcha: abnormal end ****' 
      END IF 
  100 FORMAT(' ERROR: missing definition of keyword "',a,'" (type: ',   &
     &       a,')'/8x,'in file "',a,'"')       
                      
      RETURN 
                      
! Error in reading namelist value              
   10 CONTINUE 
      lk=lench(key) 
      lz=lench(vartyp) 
      lf=lench(nmif) 
      WRITE(*,101) key(1:lk),vartyp(1:lz),krcfnm(i),nmif(1:lf) 
  101 FORMAT(' ERROR in reading value for keyword "',a,'" (type: ',     &
     &       a,')'/7x,'at record',i4,' in file "',a,'"')                
      STOP '**** rdfcha: abnormal end ****' 
    END SUBROUTINE rdfcha
! Copyright (C) 1997 by Mario Carpino (carpino@brera.mi.astro.it)       
! Version: March 8, 1997                       
! --------------------------------------------------------------------- 
!                     
!  *****************************************************************    
!  *                                                               *    
!  *                         R D F I N T                           *    
!  *                                                               *    
!  *                  Read an integer quantity                     *    
!  *           from simplified file-header namelist                *    
!  *                                                               *    
!  *****************************************************************    
!                     
! INPUT:    UNIT      -  Input unit            
!           KEY       -  Keyword               
!           REQRD     -  If true, when the keyword is not found in the  
!                        namelist, execution stops with an error message
!                     
! OUTPUT:   K         -  Integer value         
!           FOUND     -  True when the keyword has been found in the    
!                        namelist              
!           KR        -  Record number in the input file                
!                     
      SUBROUTINE rdfint(unit,key,reqrd,k,found,kr)
      USE option_input 
      IMPLICIT NONE 
      CHARACTER*(*) key 
      CHARACTER vartyp*30 
      LOGICAL reqrd,found 
      INTEGER k,kr,i,lk,lz,lf,unit 
                      
                      
      INTEGER lench 
      EXTERNAL lench 
                      
      IF(iicfnm.NE.36) STOP '**** rdfint: internal error (01) ****' 
      IF(unit.NE.hnfuni) STOP '**** rdfint: internal error (02) ****' 
                      
      vartyp='INTEGER' 
                      
      kr=0 
      DO 1 i=1,nfne 
      IF(key.EQ.keysf(i)) THEN 
          found=.true. 
          READ(valsf(i),*,ERR=10) k 
          kr=krcfnm(i) 
          kuorlf=kuorlf+1 
          kuorf(i)=kuorlf 
          RETURN 
      END IF 
    1 END DO 
                      
      found=.false. 
      IF(reqrd) THEN 
          lk=lench(key) 
          lz=lench(vartyp) 
          lf=lench(nmif) 
          WRITE(*,100) key(1:lk),vartyp(1:lz),nmif(1:lf) 
          STOP '**** rdfint: abnormal end ****' 
      END IF 
  100 FORMAT(' ERROR: missing definition of keyword "',a,'" (type: ',   &
     &       a,')'/8x,'in file "',a,'"')       
                      
      RETURN 
                      
! Error in reading namelist value              
   10 CONTINUE 
      lk=lench(key) 
      lz=lench(vartyp) 
      lf=lench(nmif) 
      WRITE(*,101) key(1:lk),vartyp(1:lz),krcfnm(i),nmif(1:lf) 
  101 FORMAT(' ERROR in reading value for keyword "',a,'" (type: ',     &
     &       a,')'/7x,'at record',i4,' in file "',a,'"')                
      STOP '**** rdfint: abnormal end ****' 
    END SUBROUTINE rdfint
! Copyright (C) 1997 by Mario Carpino (carpino@brera.mi.astro.it)       
! Version: March 8, 1997                       
! --------------------------------------------------------------------- 
!                     
!  *****************************************************************    
!  *                                                               *    
!  *                         R D F L O G                           *    
!  *                                                               *    
!  *                  Read a logical quantity                      *    
!  *           from simplified file-header namelist                *    
!  *                                                               *    
!  *****************************************************************    
!                     
! INPUT:    UNIT      -  Input unit            
!           KEY       -  Keyword               
!           REQRD     -  If true, when the keyword is not found in the  
!                        namelist, execution stops with an error message
!                     
! OUTPUT:   FLAG      -  Logical value         
!           FOUND     -  True when the keyword has been found in the    
!                        namelist              
!           KR        -  Record number in the input file                
!                     
      SUBROUTINE rdflog(unit,key,reqrd,flag,found,kr) 
      USE option_input
      IMPLICIT NONE 
      CHARACTER*(*) key 
      CHARACTER vartyp*30 
      LOGICAL reqrd,found,flag 
      INTEGER kr,i,lk,lz,lf,unit 
                      
                      
      INTEGER lench 
      EXTERNAL lench 
                      
      IF(iicfnm.NE.36) STOP '**** rdflog: internal error (01) ****' 
      IF(unit.NE.hnfuni) STOP '**** rdflog: internal error (02) ****' 
                      
      vartyp='LOGICAL' 
                      
      kr=0 
      DO 1 i=1,nfne 
      IF(key.EQ.keysf(i)) THEN 
          found=.true. 
          READ(valsf(i),*,ERR=10) flag 
          kr=krcfnm(i) 
          kuorlf=kuorlf+1 
          kuorf(i)=kuorlf 
          RETURN 
      END IF 
    1 END DO 
                      
      found=.false. 
      IF(reqrd) THEN 
          lk=lench(key) 
          lz=lench(vartyp) 
          lf=lench(nmif) 
          WRITE(*,100) key(1:lk),vartyp(1:lz),nmif(1:lf) 
          STOP '**** rdflog: abnormal end ****' 
      END IF 
  100 FORMAT(' ERROR: missing definition of keyword "',a,'" (type: ',   &
     &       a,')'/8x,'in file "',a,'"')       
                      
      RETURN 
                      
! Error in reading namelist value              
   10 CONTINUE 
      lk=lench(key) 
      lz=lench(vartyp) 
      lf=lench(nmif) 
      WRITE(*,101) key(1:lk),vartyp(1:lz),krcfnm(i),nmif(1:lf) 
  101 FORMAT(' ERROR in reading value for keyword "',a,'" (type: ',     &
     &       a,')'/7x,'at record',i4,' in file "',a,'"')                
      STOP '**** rdflog: abnormal end ****' 
    END SUBROUTINE rdflog
! Copyright (C) 1997 by Mario Carpino (carpino@brera.mi.astro.it)       
! Version: March 8, 1997                       
! --------------------------------------------------------------------- 
!                     
!  *****************************************************************    
!  *                                                               *    
!  *                         R D F R E A                           *    
!  *                                                               *    
!  *                    Read a real quantity                       *    
!  *           from simplified file-header namelist                *    
!  *                                                               *    
!  *****************************************************************    
!                     
! INPUT:    UNIT      -  Input unit            
!           KEY       -  Keyword               
!           REQRD     -  If true, when the keyword is not found in the  
!                        namelist, execution stops with an error message
!                     
! OUTPUT:   V         -  Real value            
!           FOUND     -  True when the keyword has been found in the    
!                        namelist              
!           KR        -  Record number in the input file                
!                     
      SUBROUTINE rdfrea(unit,key,reqrd,v,found,kr) 
      USE option_input
      IMPLICIT NONE 
      CHARACTER*(*) key 
      CHARACTER vartyp*30 
      LOGICAL reqrd,found 
      DOUBLE PRECISION v 
      INTEGER kr,i,lk,lz,lf,unit 
                      
                      
      INTEGER lench 
      EXTERNAL lench 
                      
      IF(iicfnm.NE.36) STOP '**** rdfrea: internal error (01) ****' 
      IF(unit.NE.hnfuni) STOP '**** rdfrea: internal error (02) ****' 
                      
      vartyp='REAL' 
                      
      kr=0 
      DO 1 i=1,nfne 
      IF(key.EQ.keysf(i)) THEN 
          found=.true. 
          READ(valsf(i),*,ERR=10) v 
          kr=krcfnm(i) 
          kuorlf=kuorlf+1 
          kuorf(i)=kuorlf 
          RETURN 
      END IF 
    1 END DO 
                      
      found=.false. 
      IF(reqrd) THEN 
          lk=lench(key) 
          lz=lench(vartyp) 
          lf=lench(nmif) 
          WRITE(*,100) key(1:lk),vartyp(1:lz),nmif(1:lf) 
          STOP '**** rdfrea: abnormal end ****' 
      END IF 
  100 FORMAT(' ERROR: missing definition of keyword "',a,'" (type: ',   &
     &       a,')'/8x,'in file "',a,'"')       
                      
      RETURN 
                      
! Error in reading namelist value              
   10 CONTINUE 
      lk=lench(key) 
      lz=lench(vartyp) 
      lf=lench(nmif) 
      WRITE(*,101) key(1:lk),vartyp(1:lz),krcfnm(i),nmif(1:lf) 
  101 FORMAT(' ERROR in reading value for keyword "',a,'" (type: ',     &
     &       a,')'/7x,'at record',i4,' in file "',a,'"')                
      STOP '**** rdfrea: abnormal end ****' 
    END SUBROUTINE rdfrea
! Copyright (C) 1997 by Mario Carpino (carpino@brera.mi.astro.it)       
! Version: March 8, 1997                       
! --------------------------------------------------------------------- 
!                     
!  *****************************************************************    
!  *                                                               *    
!  *                         R D F T I M                           *    
!  *                                                               *    
!  *                Read a time quantity (MJD+SEC)                 *    
!  *           from simplified file-header namelist                *    
!  *                                                               *    
!  *****************************************************************    
!                     
! INPUT:    UNIT      -  Input unit            
!           KEY       -  Keyword               
!           REQRD     -  If true, when the keyword is not found in the  
!                        namelist, execution stops with an error message
!                     
! OUTPUT:   TIMSTR    -  Time as a string      
!           MJD       -  Modified Julian Date (integer part)            
!           SEC       -  Seconds within the day
!           SCALE     -  Time scale            
!           FOUND     -  True when the keyword has been found in the    
!                        namelist              
!           KR        -  Record number in the input file                
!                     
      SUBROUTINE rdftim(unit,key,reqrd,timstr,mjd,sec,scale,found,kr) 
      USE option_input
      IMPLICIT NONE 
      CHARACTER*(*) key,timstr,scale 
      CHARACTER*30 vartyp 
      LOGICAL reqrd,found,error 
      DOUBLE PRECISION sec 
      INTEGER mjd,kr,i,lk,lz,lf,unit 
                      
      INTEGER lench 
      EXTERNAL lench 
                      
      IF(iicfnm.NE.36) STOP '**** rdftim: internal error (01) ****' 
      IF(unit.NE.hnfuni) STOP '**** rdftim: internal error (02) ****' 
                      
      vartyp='TIME' 
                      
      kr=0 
      DO 1 i=1,nfne 
      IF(key.EQ.keysf(i)) THEN 
          found=.true. 
          CALL ch2tim(valsf(i),mjd,sec,scale,error) 
          IF(error) GOTO 10 
          timstr=valsf(i) 
          kr=krcfnm(i) 
          kuorlf=kuorlf+1 
          kuorf(i)=kuorlf 
          RETURN 
      END IF 
    1 END DO 
                      
      found=.false. 
      IF(reqrd) THEN 
          lk=lench(key) 
          lz=lench(vartyp) 
          lf=lench(nmif) 
          WRITE(*,100) key(1:lk),vartyp(1:lz),nmif(1:lf) 
          STOP '**** rdftim: abnormal end ****' 
      END IF 
  100 FORMAT(' ERROR: missing definition of keyword "',a,'" (type: ',   &
     &       a,')'/8x,'in file "',a,'"')       
                      
      RETURN 
                      
! Error in reading namelist value              
   10 CONTINUE 
      lk=lench(key) 
      lz=lench(vartyp) 
      lf=lench(nmif) 
      WRITE(*,101) key(1:lk),vartyp(1:lz),krcfnm(i),nmif(1:lf) 
  101 FORMAT(' ERROR in reading value for keyword "',a,'" (type: ',     &
     &       a,')'/7x,'at record',i4,' in file "',a,'"')                
      STOP '**** rdftim: abnormal end ****' 
    END SUBROUTINE rdftim
! Copyright (C) 1997 by Mario Carpino (carpino@brera.mi.astro.it)       
! Version: March 8, 1997                       
! --------------------------------------------------------------------- 
!                     
!  *****************************************************************    
!  *                                      *    
!  *                         R D F R E F                           *    
!  *                                                               *    
!  *            Read a reference system description                *    
!  *           from simplified file-header namelist                *    
!  *                                                               *    
!  *****************************************************************    
!                                              
! INPUT:    UNIT      -  Input unit            
!           KEY       -  Keyword               
!           REQRD     -  If true, when the keyword is not found in the  
!                        namelist, execution stops with an error message
!                     
! OUTPUT:   RSYS      -  Reference system type (EQUM/EQUT/ECLM)         
!           EPOCH     -  Epoch specification (J2000/OFDATE)             
!           FOUND     -  True when the keyword has been found in the    
!                        namelist              
!           KR        -  Record number in the input file                
!                     
      SUBROUTINE rdfref(unit,key,reqrd,rsys,epoch,found,kr) 
      USE reference_systems 
      USE option_input
      IMPLICIT NONE 
                      
      CHARACTER*(*) key,rsys,epoch 
      LOGICAL reqrd,found 
      INTEGER kr,unit 
                      
      INTEGER i,lk,lz,lf 
      CHARACTER vartyp*30 
      LOGICAL error 
                      
      INTEGER lench 
      EXTERNAL lench 
                      
      IF(iicfnm.NE.36) STOP '**** rdfref: internal error (01) ****' 
      IF(unit.NE.hnfuni) STOP '**** rdfref: internal error (02) ****' 
                      
      vartyp='REFERENCE SYSTEM' 
                      
      kr=0 
      DO 1 i=1,nfne 
      IF(key.EQ.keysf(i)) THEN 
          found=.true. 
          CALL ch2ref(valsf(i),rsys,epoch,error) 
          IF(error) GOTO 10 
          kr=krcfnm(i) 
          kuorlf=kuorlf+1 
          kuorf(i)=kuorlf 
          RETURN 
      END IF 
    1 END DO 
                      
      found=.false. 
      IF(reqrd) THEN 
          lk=lench(key) 
          lz=lench(vartyp) 
          lf=lench(nmif) 
          WRITE(*,100) key(1:lk),vartyp(1:lz),nmif(1:lf) 
          STOP '**** rdfref: abnormal end ****' 
      END IF 
  100 FORMAT(' ERROR: missing definition of keyword "',a,'" (type: ',   &
     &       a,')'/8x,'in file "',a,'"')       
                      
      rsys=' ' 
      epoch=' ' 
                      
      RETURN 
                      
! Error in reading namelist value              
   10 CONTINUE 
      lk=lench(key) 
      lz=lench(vartyp) 
      lf=lench(nmif) 
      WRITE(*,101) key(1:lk),vartyp(1:lz),krcfnm(i),nmif(1:lf) 
  101 FORMAT(' ERROR in reading value for keyword "',a,'" (type: ',     &
     &       a,')'/7x,'at record',i4,' in file "',a,'"')                
      STOP '**** rdfref: abnormal end ****' 
    END SUBROUTINE rdfref
! Copyright (C) 1996 by Mario Carpino (carpino@brera.mi.astro.it)       
! Version: December 31, 1996                   
! --------------------------------------------------------------------- 
!                     
!  *****************************************************************    
!  *                                                               *    
!  *                         S P L K V C                           *    
!  *                                                               *    
!  *     Split a namelist record into keyword+value+comment        *    
!  *                                                               *    
!  *****************************************************************    
!                     
! INPUT:    REC       -  Record                
!                     
! OUTPUT:   KEY       -  Keyword field         
!           VAL       -  Value field           
!           COMM      -  Comment field         
!           SKIP      -  "Only comment" flag   
!           END       -  "End of namelist" flag
!                     
      SUBROUTINE splkvc(rec,key,val,comm,skip,end) 
      IMPLICIT NONE 
                      
      INCLUDE 'parcmc.h90' 
                      
      CHARACTER*(*) rec,key,val,comm 
      LOGICAL skip,end 
                      
      CHARACTER*200 rec1,tmp 
      INTEGER lr,lk,lv,ipc,ipu 
                      
      INTEGER lench 
      EXTERNAL lench 
                      
      key=' ' 
      val=' ' 
      comm=' ' 
                      
! Detection of the end of the namelist         
      rec1=rec 
      CALL rmsp(rec1,lr) 
      IF(rec1(1:lr).EQ.'END_OF_HEADER') THEN 
          end=.true. 
          skip=.false. 
          RETURN 
      ELSE 
          end=.false. 
      END IF 
                      
! Compute length excluding comments            
      lr=lench(rec) 
      IF(lr.LT.1) THEN 
          skip=.true. 
          RETURN 
      END IF 
      ipc=INDEX(rec(1:lr),comcha) 
      IF(ipc.EQ.1) THEN 
          skip=.true. 
          comm=rec(2:) 
          RETURN 
      END IF 
                      
! Separate comment from keyword+value          
      IF(ipc.EQ.0) THEN 
          rec1=rec(1:lr) 
      ELSE 
          rec1=rec(1:ipc-1) 
          comm=rec(ipc+1:) 
          lr=lench(rec1) 
          IF(lr.LT.1) THEN 
              skip=.true. 
              RETURN 
          END IF 
      END IF 
      skip=.false. 
                      
! Keyword field       
      ipu=INDEX(rec1(1:lr),'=') 
      IF(ipu.EQ.0) THEN 
          CALL rmsp(rec1,lk) 
          IF(lk.GT.LEN(key)) STOP '**** splkvc: lk > LEN(key) ****' 
          key=rec1(1:lk) 
          RETURN 
      END IF 
      tmp=rec1(1:ipu-1) 
      CALL rmsp(tmp,lk) 
      IF(lk.GT.LEN(key)) STOP '**** splkvc: lk > LEN(key) ****' 
      key=tmp(1:lk) 
                      
! Value field         
      tmp=rec1(ipu+1:) 
      CALL norstr(tmp,lv) 
      IF(lv.GT.LEN(val)) STOP '**** splkvc: lv > LEN(val) ****' 
      val=tmp(1:lv) 
                      
    END SUBROUTINE splkvc
! Copyright (C) 1996-1998 by Mario Carpino (carpino@brera.mi.astro.it)  
! Version: January 13, 1998                    
! --------------------------------------------------------------------- 
!                     
!  *****************************************************************    
!  *                                                               *    
!  *                         C H K F L N                           *    
!  *                                                               *    
!  *                     Check field length                        *    
!  *                                                               *    
!  *****************************************************************    
!                     
! INPUT:    NA        -  Actual value required for the parameter        
!           NX        -  Value declared for the parameter               
!           NAME      -  Parameter name        
!           KR        -  Record number in input file                    
!           FILE      -  Input file name       
!                     
      SUBROUTINE chkfln(na,nx,name,kr,file) 
      IMPLICIT NONE 
                      
      INTEGER na,nx,kr,ll,lf 
      CHARACTER*(*) name,file 
                      
      INTEGER lench 
      EXTERNAL lench 
                      
      IF(na.LE.nx) RETURN 
                      
      ll=lench(name) 
      lf=lench(file) 
      WRITE(*,100) name(1:ll),nx,kr,file(1:lf),na 
  100 FORMAT(' **** Insufficient PARAMETER definition ****'/            &
     & ' Sorry, the present version of the program does not',           &
     & ' allow'/      &
     & ' a length of the "',a,'" field longer than',i4,' characters'/   &
     & ' Please correct line',i4,' in file ',a/&
     & ' (length of the field =',i4,' characters)')                     
      STOP '**** chkfln: abnormal end ****' 
    END SUBROUTINE chkfln
                      
! Copyright (C) 1997 by Mario Carpino (carpino@brera.mi.astro.it)       
! Version: March 8, 1997                       
! --------------------------------------------------------------------- 
!                     
!  *****************************************************************    
!  *                                                               *    
!  *                         G E T R S C                           *    
!  *                                                               *    
!  *         Get a data record skipping comment lines              *    
!  *                                                               *    
!  *****************************************************************    
!                     
! INPUT:    UNIT      -  Input unit            
!           NR        -  Number of records read so far                  
!                     
! OUTPUT:   REC       -  Record                
!           NR        -  Number of records read so far (updated)        
!           END       -  End of file reached   
!                     
      SUBROUTINE getrsc(unit,rec,nr,end) 
      IMPLICIT NONE 
                      
      INCLUDE 'parcmc.h90' 
                      
      INTEGER unit,nr 
      CHARACTER*(*) rec 
      LOGICAL end 
                      
      end=.false. 
                      
    1 CONTINUE 
      READ(unit,100,END=2) rec 
  100 FORMAT(A) 
      nr=nr+1 
      IF(rec(1:1).EQ.comcha) GOTO 1 
      RETURN 
                      
    2 CONTINUE 
      end=.true. 
      rec=' ' 
                      
    END SUBROUTINE getrsc
! Copyright (C) 1997 by Mario Carpino (carpino@brera.mi.astro.it)       
! Version: June 23, 1997                       
! --------------------------------------------------------------------- 
!                     
!  *****************************************************************    
!  *                                                               *    
!  *                         S V 2 I N T                           *    
!  *                                                               *    
!  *        Translation of string-valued keywords to integer       *    
!  *                                                               *    
!  *****************************************************************    
!                     
! INPUT:    CAT       -  Namelist keyword category                      
!           KEY       -  Namelist keyword name 
!           REQRD     -  If true, when the keyword is not found in the  
!                        namelist, an error message is issued and the   
!                        FAIL flag is set      
!                     
! OUTPUT:   V         -  Value associated with keyword                  
!           IVK       -  Integer translation of the Value associated    
!                        with keyword          
!           FOUND     -  True when the keyword has been found in the    
!                        namelist              
!           FAIL1     -  True if(.NOT.found.AND.reqrd); false otherwise 
!           FAIL      -  Set to true if(.NOT.found.AND.reqrd);          
!                        otherwise, not changed
!                     
      SUBROUTINE sv2int(cat,key,v,ivk,reqrd,found,fail1,fail)
      USE option_input 
      IMPLICIT NONE 
      CHARACTER*(*) cat,key 
      CHARACTER ck*(lckx),val*(lcvx),file*(lcfx),rest*(lcvx) 
      CHARACTER*50 vartyp 
      LOGICAL reqrd,found,fail,fail1,error 
      CHARACTER v*(*) 
      INTEGER ivk,ktyp0,lc,lk,lt,ktyp1,kr,lz,lf,i,lk2,lv1,k,lv2,kk 
                      
      INTEGER lench 
      EXTERNAL lench 
                      
      vartyp='CHARACTER' 
      ktyp0=3 
                      
      IF(iickls.NE.36) STOP '**** sv2int: internal error (01) ****' 
                      
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
              ivk=0 
          END IF 
          RETURN 
      END IF 
      lf=lench(file) 
      IF(ktyp1.NE.ktyp0) THEN 
          WRITE(*,106) ck(1:lt),ktyp0,ktyp1 
          STOP '**** sv2int: abnormal end ****' 
      END IF 
  100 FORMAT(' ERROR: Missing definition of keyword "',a,'" (type: ',   &
     &       a,')')   
  106 FORMAT(' ERROR: keyword type mismatch'/  &
     &       '        KEY = "',A,'"'/          &
     &       '        KTYP0 =',I2/             &
     &       '        KTYP1 =',I2)             
      CALL strcnt(val,v,rest,error) 
      IF(error) GOTO 1 
                      
! Search keyword in list                       
      DO 10 i=1,nkls 
      lk2=lench(keylst(i)) 
      IF(ck(1:lt).EQ.keylst(i)(1:lk2)) GOTO 2 
   10 END DO 
      STOP '**** sv2int: internal error (02) ****' 
    2 CONTINUE 
                      
! Search string value 
      IF(ns2i(i).LE.0) STOP '**** sv2int: internal error (03) ****' 
      lv1=lench(v) 
      DO 3 k=1,ns2i(i) 
      kk=k+ipos2i(i) 
      lv2=lench(vallst(kk)) 
      IF(v(1:lv1).NE.vallst(kk)(1:lv2)) GOTO 3 
      ivk=intlst(kk) 
      RETURN 
    3 END DO 
                      
      WRITE(*,103) ck(1:lt),vartyp(1:lz),kr,file(1:lf),v(1:lv1) 
  103 FORMAT(' ERROR: unexpected value of keyword "',a,'" (type: ',     &
     &       a,')'/8x,'(record',i4,' in file "',a,'")'/                 &
     &       8x,'Present value: ''',a,''''/    &
     &       8x,'Possible values are:')        
      DO 4 k=1,ns2i(i) 
      kk=k+ipos2i(i) 
      lv2=lench(vallst(kk)) 
      WRITE(*,102) vallst(kk)(1:lv2) 
  102 FORMAT(13x,'''',a,'''') 
    4 END DO 
                      
      fail1=.true. 
      fail=.true. 
      RETURN 
                      
! Error message       
    1 CONTINUE 
      WRITE(*,101) ck(1:lt),vartyp(1:lz),kr,file(1:lf) 
  101 FORMAT(' ERROR: Abnormal definition of keyword "',a,'" (type: ',  &
     &       a,')'/8x,'(record',i4,' in file "',a,'")')                 
      fail1=.true. 
      fail=.true. 
                                               
    END SUBROUTINE sv2int
! ====================================================
! ROUTINES FOR OPTION INPUT
! ====================================================
! INITOPT
! initializes namelist inoput, force model, Gauss' method
! ====================================================
SUBROUTINE initopt(progna0,run0,suffix)  
  IMPLICIT NONE
  INCLUDE 'comlib.h90' 

  CHARACTER*6,INTENT(IN) :: progna0 
  CHARACTER*80,INTENT(IN) :: run0  
  CHARACTER*3, INTENT(IN) :: suffix
  CHARACTER*6 progna
  character*80 run
  character*120 file, nam1 
  integer iunit,lerun,le
  LOGICAL found
! read option for propagator                                            
  CALL libini 
  CALL namini 
! default options for propag                                   
  CALL filopl(iunit,'propag.def') 
  CALL rdnam(iunit) 
  CALL filclo(iunit,' ') 
! =============================                                         
  progna=progna0                                       
  call rmsp(progna,le) 
! default options for progna   
  file=progna(1:le)//'.def' 
  INQUIRE(FILE=file,EXIST=found)
  IF(.not.found)THEN
     nam1=libdir(1:lenld)//file 
     INQUIRE(FILE=nam1,EXIST=found)
  ENDIF
  IF(found)THEN
     CALL filopl(iunit,file(1:le+4)) 
     CALL rdnam(iunit) 
     CALL filclo(iunit,' ')
  ENDIF
! =============================                                         
! particular options for this run
  run=run0                                       
  CALL rmsp(run,lerun) 
  file=run(1:lerun)//'.'//suffix 
  INQUIRE(FILE=file,EXIST=found) 
  IF(found) THEN 
     CALL filopn(iunit,file,'OLD') 
     CALL rdnam(iunit) 
     CALL filclo(iunit,' ') 
  ELSE 
     write(*,*)'**** file not found: ',file 
     write(*,*)'******* ',progna0,' using defaults ****' 
  ENDIF
! possible options for progna
  file=progna(1:le)//'.key' 
  CALL rdklst(file) 
! check for non-existing options 
  CALL chkkey 
! ==============================  
END SUBROUTINE initopt

SUBROUTINE input_cha_opt(progna,nameopt,optval,ireq,found,comment,iunout)
  IMPLICIT NONE
  CHARACTER*6, INTENT(IN) :: progna ! name of the .mop file
  CHARACTER*(*), INTENT(IN) :: nameopt ! charcater option to be found
  CHARACTER*(*), INTENT(INOUT) :: optval ! value of the option; 

         ! inout because it could have a preassigned default value
  LOGICAL, INTENT(IN) :: ireq ! required or optional?
  LOGICAL, INTENT(OUT) :: found ! found in namelist input
  CHARACTER*(*), INTENT(IN) :: comment ! to be written in output
  INTEGER, INTENT(IN) :: iunout ! unit for log (if positive)
  CHARACTER*7 prognp
  INTEGER le,leco,lench,le1
  LOGICAL fail, fail1
! ----------------------------------------------
  fail=.false. 
  prognp=progna//'.' 
  CALL rmsp(prognp,le)
  leco=lench(comment)
  fail=.false.
  CALL rdncha(prognp,nameopt,optval,ireq,found,fail1,fail)
  IF(fail)THEN
     WRITE(*,*)prognp(1:le-1)//'.mop : ',comment(1:leco),' failure in reading'
     STOP
  ENDIF
  CALL rmsp(optval,le1)
  IF(.not.found)THEN
     IF(ireq)THEN 
        WRITE(*,*)prognp(1:le-1)//'.mop : ',comment(1:leco),' required' 
        STOP 
     ELSE
        WRITE(*,*)comment(1:leco),' default value ',optval(1:le1)
  IF(iunout.gt.0)WRITE(iunout,*)comment(1:leco),' default value= ',optval(1:le1)
     ENDIF
  ELSE
     WRITE(*,*) comment(1:leco),'= ', optval(1:le1) 
     IF(iunout.gt.0)WRITE(iunout,*) comment(1:leco),'= ', optval(1:le1)
  ENDIF 
END SUBROUTINE input_cha_opt

SUBROUTINE input_rea_opt(progna,nameopt,optval,ireq,found,comment,iunout)
  IMPLICIT NONE
  CHARACTER*6, INTENT(IN) :: progna ! name of the .mop file
  CHARACTER*(*), INTENT(IN) :: nameopt ! charcater option to be found
  DOUBLE PRECISION, INTENT(INOUT) :: optval ! value of the option; 
         ! inout because it could have a preassigned default value
  LOGICAL, INTENT(IN) :: ireq ! required or optional?
  LOGICAL, INTENT(OUT) :: found ! found in namelist input
  CHARACTER*(*), INTENT(IN) :: comment ! to be written in output
  INTEGER, INTENT(IN) :: iunout ! unit for log (if positive)
  CHARACTER*7 prognp
  INTEGER le,leco,lench
  LOGICAL fail, fail1
! -------------
  fail=.false. 
  prognp=progna//'.' 
  CALL rmsp(prognp,le)
  leco=lench(comment)
  fail=.false.
  CALL rdnrea(prognp,nameopt,optval,ireq,found,fail1,fail)
  IF(fail)THEN
     WRITE(*,*)prognp(1:le-1)//'.mop : ',comment(1:leco),' failure in reading'
     STOP
  ENDIF
  IF(.not.found)THEN
     IF(ireq)THEN 
        WRITE(*,*)prognp(1:le-1)//'.mop : ',comment(1:leco),' required' 
        STOP 
     ELSE
        WRITE(*,*)comment(1:leco),' default value ',optval
        IF(iunout.gt.0)WRITE(iunout,*)comment(1:leco),' default value ',optval
     ENDIF
  ELSE
     WRITE(*,*) comment(1:leco),' = ', optval 
     IF(iunout.gt.0)WRITE(iunout,*) comment(1:leco),' = ', optval
  ENDIF 
END SUBROUTINE input_rea_opt
!
SUBROUTINE input_int_opt(progna,nameopt,optval,ireq,found,comment,iunout)
  IMPLICIT NONE
  CHARACTER*6, INTENT(IN) :: progna ! name of the .mop file
  CHARACTER*(*), INTENT(IN) :: nameopt ! charcater option to be found
  INTEGER, INTENT(INOUT) :: optval ! value of the option; 
         ! inout because it could have a preassigned default value
  LOGICAL, INTENT(IN) :: ireq ! required or optional?
  LOGICAL, INTENT(OUT) :: found ! found in namelist input
  CHARACTER*(*), INTENT(IN) :: comment ! to be written in output
  INTEGER, INTENT(IN) :: iunout ! unit for log (if positive)
  CHARACTER*7 prognp
  INTEGER le,leco,lench
  LOGICAL fail, fail1
! -------------
  fail=.false. 
  prognp=progna//'.' 
  CALL rmsp(prognp,le)
  leco=lench(comment)
  fail=.false.
  CALL rdnint(prognp,nameopt,optval,ireq,found,fail1,fail)
  IF(fail)THEN
     WRITE(*,*)prognp(1:le-1)//'.mop : ',comment(1:leco),' failure in reading'
     STOP
  ENDIF
  IF(.not.found)THEN
     IF(ireq)THEN 
        WRITE(*,*)prognp(1:le-1)//'.mop : ',comment(1:leco),' required' 
        STOP 
     ELSE
        WRITE(*,*)comment(1:leco),' default value ',optval
        IF(iunout.gt.0)WRITE(iunout,*)comment(1:leco),' default value ',optval
     ENDIF
  ELSE
     WRITE(*,*) comment(1:leco),' = ', optval 
     IF(iunout.gt.0)WRITE(iunout,*) comment(1:leco),' = ', optval
  ENDIF 
END SUBROUTINE input_int_opt
!
SUBROUTINE input_log_opt(progna,nameopt,optval,ireq,found,comment,iunout)
  IMPLICIT NONE
  CHARACTER*6, INTENT(IN) :: progna ! name of the .mop file
  CHARACTER*(*), INTENT(IN) :: nameopt ! character option to be found
  LOGICAL, INTENT(INOUT) :: optval ! value of the option; 
         ! inout because it could have a preassigned default value
  LOGICAL, INTENT(IN) :: ireq ! required or optional?
  LOGICAL, INTENT(OUT) :: found ! found in namelist input
  CHARACTER*(*), INTENT(IN) :: comment ! to be written in output
  INTEGER, INTENT(IN) :: iunout ! unit for log (if positive)
  CHARACTER*7 prognp
  INTEGER le,leco,lench
  LOGICAL fail, fail1
! -------------
  fail=.false. 
  prognp=progna//'.' 
  CALL rmsp(prognp,le)
  leco=lench(comment)
  fail=.false.
  CALL rdnlog(prognp,nameopt,optval,ireq,found,fail1,fail)
  IF(fail)THEN
     WRITE(*,*)prognp(1:le-1)//'.mop : ',comment(1:leco),' failure in reading'
     STOP
  ENDIF
  IF(.not.found)THEN
     IF(ireq)THEN 
        WRITE(*,*)prognp(1:le-1)//'.mop : ',comment(1:leco),' required' 
        STOP 
     ELSE
        WRITE(*,*)comment(1:leco),' default value ',optval
        IF(iunout.gt.0)WRITE(iunout,*)comment(1:leco),' default value ',optval
     ENDIF
  ELSE
     WRITE(*,*) comment(1:leco),' is ', optval 
     IF(iunout.gt.0)WRITE(iunout,*) comment(1:leco),' is ', optval
  ENDIF 
END SUBROUTINE input_log_opt

! ===================================================================   
! TRIVOPT minimum option routine, without propagation!!!
! ===================================================================   
! input options, for the propagator and the specific main program       
! input: progna = program name (6 characters)                           
!        run    = run identifier (80 characters, up to 76 non-blank)    
      SUBROUTINE trivopt(progna,run,iun20) 
      implicit none 
      character*6 progna 
      character*80 run 
      integer iun20 
! ==========END INTERFACE============================================   
      integer le,iunit 
      character*100 file 
      logical found 
! =============================                                         
! read option for propagator                                            
      CALL libini 
      CALL namini 
! default options are only for propag                                   
      CALL filopl(iunit,'propag.def') 
      CALL rdnam(iunit) 
      CALL filclo(iunit,' ') 
! =============================                                         
! particular options for this run                                       
      file=run//'.mop' 
      CALL rmsp(file,le) 
      INQUIRE(FILE=file(1:le),EXIST=found) 
      IF(found) THEN 
        CALL filopn(iunit,file(1:le),'OLD') 
        CALL rdnam(iunit) 
        CALL filclo(iunit,' ') 
      ELSE 
        write(*,*)'**** file not found: ',file 
        write(*,*)'******* ',progna,' only default ****' 
      ENDIF 
! =============================                                         
! check for non-existing options                                        
      call rmsp(progna,le) 
      file=progna(1:le)//'.key' 
      CALL rdklst(file) 
      CALL chkkey                                                      
! ===================================================================== 
! Output files: for control and results, for covariance                 
      CALL rmsp(run,le) 
      file=run(1:le)//'.mou' 
      call filopn(iun20,file,'UNKNOWN') 
      END SUBROUTINE trivopt                                          

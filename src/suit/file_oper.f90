! ===========MODULE FILE_OPER=====================                      

MODULE file_oper
! former COMFIL.H
! Copyright (C) 1996 by Mario Carpino (carpino@brera.mi.astro.it)
! Version: August 10, 1996
! ---------------------------------------------------------------------
! File names and assigned units
!
! filnam      -  File names
! allunt      -  Allocation indicator
! iicfil      -  Initialization check
!
  INTEGER iicfil
  INTEGER, PARAMETER :: iunf1=50
  INTEGER, PARAMETER :: iunf2=99
  CHARACTER*80 filna_fo(iunf1:iunf2)
  LOGICAL allunt(iunf1:iunf2)
  PUBLIC iicfil,iunf1,iunf2,filnam,allunt

END MODULE file_oper

! out of module:
! FILE OPENING/CLOSE AND UNIT ASSIGNMENT:                               
! CONTAINS                                                              
! SUBROUTINES                                                           
!                                                                       
! filopn	file opening with unit assignment                              
! filclo	file closing with unit release                                 
! filass	unit assignment without file opening                           
! filopl	unit assignment and file opening (with library search)         
! filopf	unit assignment and file opening (with library search and      
!		no abort if the file does not exist) 
! dlifex	delete a file if esists                                        
! libini	inizialization of default library directory                    
! filna_fo        composition of file name from dir, name suffix          
! splinam       split asteroid name (identification, multiple solution) 
! fidinam       find path for file (moved to fidinam.f90) 

!                                                                       
!                                                                       
! Copyright (C) 1997 by Mario Carpino (carpino@brera.mi.astro.it)       
! Version: December 5, 1997                                             
! --------------------------------------------------------------------- 
!                                                                       
!  *****************************************************************    
!  *                                                               *    
!  *                         F I L O P N                           *    
!  *                                                               *    
!  *               Unit allocation and file opening                *    
!  *                                                               *    
!  *****************************************************************    
!                                                                       
! INPUT:    NAME      -  File name to be opened                         
!           STATUS    -  Open status                                    
!                                                                       
! OUTPUT:   IUN       -  Allocated unit                                 
!                                                                       
SUBROUTINE filopn(iun,name,status) 
  USE file_oper
  IMPLICIT NONE 
  INTEGER iun,ll,ls,i,i0,iunold 
  CHARACTER*(*) name,status 
  CHARACTER*100 name_rescue
  CHARACTER*40 name1
  LOGICAL opnd,exis, nmd 
                                                                        
  INTEGER lench,io_stat_no,ntry 
  EXTERNAL lench 
  LOGICAL first
  DATA first /.true./
  SAVE first,i0
  IF(first)THEN
     i0=iunf1-1 
     first=.false.
  ENDIF
  IF(iicfil.NE.36) THEN 
     DO  i=iunf1,iunf2 
        allunt(i)=.false.
     ENDDO 
     iicfil=36 
  END IF
  ntry=0
7 INQUIRE(FILE=name,OPENED=opnd,EXIST=exis,NUMBER=iunold) 
  IF(opnd) THEN 
     ll=lench(name) 
     WRITE(*,102) name(1:ll), opnd,exis,iunold,allunt(iunold)
102 FORMAT(' **** filopn: INQUIRE error ****'/                  &
     &       ' **** FILE: ',A,' opnd=',L1,' exis=',L1,' ****'/        &
     &       ' **** UNIT: ',i3,' allunt; ',i3)
     WRITE(*,*)'NUMBER=',iunold
     INQUIRE(UNIT=iunold,OPENED=opnd,NAMED=nmd,NAME=name1,EXIST=exis,ERR=99)
     WRITE(*,*)' UNIT= ',iunold,' NAMED= ',nmd,' NAME= ',name1
!     CLOSE(iunold,ERR=4,IOSTAT=io_stat_no)
!     GOTO 7
     GOTO 5
!!          STOP '**** filopn: abnormal end ****' 
  ELSE
     GOTO 5                       
  END IF
!4 WRITE(*,*)' **** filopn: CLOSE error, IOSTAT: ',io_stat_no
5 CONTINUE                                                          
  DO  iun=i0+1,iunf2 
     IF(allunt(iun)) CYCLE 
     OPEN(iun,FILE=name,STATUS=status,ERR=3,IOSTAT=io_stat_no) 
     filna_fo(iun)=name 
     allunt(iun)=.true. 
     i0=iun
     ll=lench(name)
!     WRITE(*,*)'filopn: ', iun, name(1:ll)
     RETURN 
  ENDDO
  DO  iun=iunf1,i0-1 
     IF(allunt(iun)) CYCLE 
     OPEN(iun,FILE=name,STATUS=status,ERR=3,IOSTAT=io_stat_no) 
     filna_fo(iun)=name 
     allunt(iun)=.true.
     i0=iun 
     ll=lench(name)
!     WRITE(*,*)'filopn: ', iun, name(1:ll)
     RETURN 
  ENDDO
  
  WRITE(*,*)'filopn: all units are already allocated'
  STOP
                                                                        
3 CONTINUE 
  ll=lench(name) 
  ls=lench(status) 
  WRITE(*,101) name(1:ll),status(1:ls),iun,exis,opnd,io_stat_no 
101 FORMAT(' **** filopn: cannot OPEN file "',A,'" (status=',A,       &
     &    '   as ',I4,' EXISTS=',L1,' OPENED=',L1, ' IOSTAT=',i8,') ****') 
  STOP '**** filopn: abnormal end ****' 
! try desperate move
!   name_rescue='rescue/'//name
!   CALL rmsp(name_rescue,ll)
!   OPEN(iun,FILE=name_rescue(1:ll),STATUS=status,ERR=33,IOSTAT=io_stat_no) 
!   WRITE(*,*)'filopn: rescue attempt ', iun, name_rescue(1:ll)
!   filna_fo(iun)=name 
!   allunt(iun)=.true.
!   i0=iun 
!    RETURN
!33  ls=lench(status) 
!    WRITE(*,103) name_rescue(1:ll),status(1:ls),iun,io_stat_no
!103 FORMAT(' **** filopn: cannot OPEN rescue file "',A,'" (status=',A,       &
!      &    '   as ',I4,' IOSTAT=',i8,') ****') 
!    ntry=ntry+1
!    IF(ntry.gt.10)  STOP '**** filopn: abnormal end ****'
!    CALL waste_time(100) 
!    GOTO 5
99 continue
  WRITE(*,*)' error in second INQUIRE, file=',name(1:ll)
  STOP '**** filopn: abnormal end ****' 
END SUBROUTINE filopn
SUBROUTINE waste_time(n)
  IMPLICIT NONE
  INTEGER,PARAMETER :: dim=100
  INTEGER,INTENT(IN) ::n
  DOUBLE PRECISION a(dim,dim),b(dim,dim),c(dim,dim)
  INTEGER i
  a=1.d0
  b=2.d0
  c=0.d0
  OPEN(22,file='waste.out',STATUS='UNKNOWN')
  DO i=1,n
     c=c+MATMUL(a,b)
     write(22,*)c
  ENDDO
  CLOSE(22)
END SUBROUTINE waste_time

! Copyright (C) 1996-1998 by Mario Carpino (carpino@brera.mi.astro.it)  
! Version: January 12, 1998                                             
! --------------------------------------------------------------------- 
!                                                                       
!  *****************************************************************    
!  *                                                               *    
!  *                         F I L C L O                           *    
!  *                                                               *    
!  *                File closing and unit release                  *    
!  *                                                               *    
!  *****************************************************************    
!                                                                       
! INPUT:    IUN       -  Unit to be closed                              
!           STATUS    -  Close status (if blank, it is not used)        
!                                                                       
SUBROUTINE filclo(iun,status) 
  USE file_oper
  IMPLICIT NONE 
  INTEGER iun,ls 
  CHARACTER*(*) status 
                                                                        
  INTEGER lench,le 
  EXTERNAL lench 
                                                                        
  IF(iicfil.NE.36) STOP '**** filclo: internal error (01) ****' 
  IF(iun.LT.iunf1.OR.iun.GT.iunf2)                                  &
       &          STOP '**** filclo: internal error (02) ****'            
  IF(.NOT.allunt(iun)) THEN 
     WRITE(*,200) iun 
     STOP 
  END IF
200 FORMAT(' **** filclo: unit',i4,' is not opened ****') 
  CALL flush(iun)  
  ls=lench(status) 
  IF(ls.LE.0) THEN 
     CLOSE(iun) 
  ELSE 
     CLOSE(iun,STATUS=status) 
  END IF
!  CALL rmsp(filna_fo(iun),le)
!  WRITE(*,*)'filclo: ',iun,filna_fo(iun)(1:le)
  allunt(iun)=.false. 
  filna_fo(iun)=' ' 
                                                                          
END SUBROUTINE filclo
! Copyright (C) 1996 by Mario Carpino (carpino@brera.mi.astro.it)       
! Version: September 19, 1996                                           
! --------------------------------------------------------------------- 
!                                                                       
!  *****************************************************************    
!  *                                                               *    
!  *                         F I L A S S                           *    
!  *                                                               *    
!  *               Unit allocation for file opening                *    
!  *                                                               *    
!  *****************************************************************    
!                                                                       
! INPUT:    NAME      -  File name to be opened                         
!                                                                       
! OUTPUT:   IUN       -  Allocated unit                                 
!                                                                       
SUBROUTINE filass(iun,name) 
  USE file_oper
  IMPLICIT NONE 
  INTEGER iun,i 
  CHARACTER*(*) name 
                                                                        
  IF(iicfil.NE.36) THEN 
     allunt(iunf1:iunf2)=.false.
     iicfil=36 
  END IF 
                                                                        
  DO 2 iun=iunf1,iunf2 
     IF(allunt(iun)) GOTO 2 
     filna_fo(iun)=name 
     allunt(iun)=.true. 
     RETURN 
2 END DO
  STOP '**** filass: all units are already allocated ****' 
END SUBROUTINE filass
! Copyright (C) 1996 by Mario Carpino (carpino@brera.mi.astro.it)       
! Version: September 19, 1996                                           
! --------------------------------------------------------------------- 
!                                                                       
!  *****************************************************************    
!  *                                                               *    
!  *                         F I L O P L                           *    
!  *                                                               *    
!  *        Unit allocation and file opening (STATUS='old')        *    
!  *         (searching also in default library directory)         *    
!  *                                                               *    
!  *****************************************************************    
!                                                                       
! INPUT:    NAME      -  File name to be opened                         
!                                                                       
! OUTPUT:   IUN       -  Allocated unit                                 
!                                                                       
SUBROUTINE filopl(iun,name) 
  USE file_oper
  IMPLICIT NONE 
! NEEDED common blocks:                                                 
  INCLUDE 'comlib.h90' 
  INTEGER iun,ll,i 
  CHARACTER*(*) name 
  CHARACTER*120 nam1 
  LOGICAL found 
  INTEGER lench 
  EXTERNAL lench 
  IF(iiclib.NE.36) STOP '**** filopl: internal error (01) ****' 
  IF(iicfil.NE.36) THEN 
    allunt(iunf1:iunf2)=.false. 
    iicfil=36 
 END IF
 nam1=name 
 INQUIRE(FILE=nam1,EXIST=found) 
 IF(.NOT.found) THEN 
    nam1=libdir(1:lenld)//name 
    INQUIRE(FILE=nam1,EXIST=found) 
    IF(.NOT.found) THEN 
       ll=lench(name) 
       WRITE(*,102) name(1:ll) 
       STOP '**** filopl: abnormal end ****' 
    END IF
 END IF
102 FORMAT(' **** filopl: cannot find file "',a,'" ****') 
                                                                        
 DO 2 iun=iunf1,iunf2 
    IF(allunt(iun)) GOTO 2 
    OPEN(iun,FILE=nam1,STATUS='old',ERR=3) 
    filna_fo(iun)=nam1 
    allunt(iun)=.true. 
    ll=lench(name)
!      WRITE(*,*)' filopl: ',iun,name(1:ll)
    RETURN 
2 END DO
 STOP '**** filopl: all units are already allocated ****' 
3 CONTINUE 
 ll=lench(nam1) 
 WRITE(*,101) nam1(1:ll) 
101 FORMAT(' **** filopl: cannot OPEN file "',a,'" (STATUS=old) ****') 
 STOP '**** filopl: abnormal end ****' 
END SUBROUTINE filopl
! Copyright (C) 1996-1998 by Mario Carpino (carpino@brera.mi.astro.it)  
! Version: January 12, 1998                                             
! --------------------------------------------------------------------- 
!                                                                       
!  *****************************************************************    
!  *                                                               *    
!  *                         L I B I N I                           *    
!  *                                                               *    
!  *         Inizialization of default library directory           *    
!  *                                                               *    
!  *****************************************************************    
!                                                                       
SUBROUTINE libini 
  IMPLICIT NONE 
  INCLUDE 'sysdep.h90' 
  INCLUDE 'parlib.h90' 
! Common blocks to be initialized:                                      
  INCLUDE 'comlib.h90' 
  LOGICAL found 
  INTEGER unit 
  INTEGER lench 
  EXTERNAL lench 
  CHARACTER(LEN=200) orbfit_home_env

! The file 'libdir.dat' can be used to modify the path of the           
! library directory (with respect to the built-in value contained       
! in parlib.h) without need of compiling again the software:            
! it is searched only in the current (working) directory                
  INQUIRE(FILE='libdir.dat',EXIST=found) 

! First try to get the home dir from environment, then
! file in current dir, then configured hard path.
!  CALL get_environment_variable("ORBFIT_HOME", dlibd, 99, env_status)
  CALL getenv("ORBFIT_DATA", orbfit_home_env)
!  WRITE(*,*)'ENV=',orbfit_home_env)
  IF(len(trim(orbfit_home_env)).gt.0) THEN
    libdir = orbfit_home_env
  ELSE
     libdir=dlibd
  END IF
  IF(found) THEN
     CALL filopn(unit,'libdir.dat','old')
     READ(unit,100,END=10,ERR=10) libdir
     CALL filclo(unit,' ')
  END IF
100 FORMAT(A)
!  WRITE(*,*)'libdir=',libdir
  lenld=lench(libdir)
  IF(libdir(lenld:lenld).NE.dircha) THEN
     lenld=lenld+1
     libdir(lenld:lenld)=dircha
  END IF
  iiclib=36
  RETURN
10 CONTINUE
  STOP '**** libini: error reading file "libdir.dat" ****'

END SUBROUTINE libini
! Copyright (C) 1996 by Mario Carpino (carpino@brera.mi.astro.it)       
! Version: September 19, 1996                                           
! --------------------------------------------------------------------- 
!                                                                       
!  *****************************************************************    
!  *                                                               *    
!  *                         F I L O P F                           *    
!  *                                                               *    
!  *        Unit allocation and file opening (STATUS='old')        *    
!  *         (searching also in default library directory)         *    
!  *                                                               *    
!  *****************************************************************    
!                                                                       
! INPUT:    NAME      -  File name to be opened                         
!                                                                       
! OUTPUT:   IUN       -  Allocated unit                                 
!           FOUND     -  Was the file found?                            
!                                                                       
SUBROUTINE filopf(iun,name,found) 
  USE file_oper
  IMPLICIT NONE 
! NEEDED common blocks:                                                 
  INCLUDE 'comlib.h90' 
  INTEGER iun,ll,i 
  CHARACTER*(*) name 
  CHARACTER*120 tname 
  LOGICAL found 
  INTEGER lench 
  EXTERNAL lench 
  IF(iiclib.NE.36) STOP '**** filopf: internal error (01) ****' 
  IF(iicfil.NE.36) THEN 
     allunt(iunf1:iunf2)=.false. 
     iicfil=36 
  END IF
  tname=name 
  INQUIRE(FILE=tname,EXIST=found) 
  IF(.NOT.found) THEN 
     tname=libdir(1:lenld)//name 
     INQUIRE(FILE=tname,EXIST=found) 
     IF(.NOT.found) THEN 
        iun=0 
        RETURN 
     END IF
  END IF
  DO 2 iun=iunf1,iunf2 
     IF(allunt(iun)) GOTO 2 
     OPEN(iun,FILE=tname,STATUS='old',ERR=3) 
     filna_fo(iun)=tname 
     allunt(iun)=.true. 
     RETURN 
2 END DO
  STOP '**** filopf: all units are already allocated ****' 
3 CONTINUE 
  ll=lench(tname) 
  WRITE(*,101) tname(1:ll) 
101 FORMAT(' **** filopf: cannot OPEN file "',a,'" (STATUS=old) ****') 
  STOP '**** filopf: abnormal end ****' 
END SUBROUTINE filopf
! Copyright (C) 1996 by Mario Carpino (carpino@brera.mi.astro.it)       
! Version: August 27, 1996                                              
! --------------------------------------------------------------------- 
!                                                                       
!  *****************************************************************    
!  *                                                               *    
!  *                         D L I F E X                           *    
!  *                                                               *    
!  *                Delete a file (if it exists)                   *    
!  *                                                               *    
!  *****************************************************************    
!                                                                       
! INPUT:    FILE      -  Name of the file to be deleted                 
!                                                                       
SUBROUTINE dlifex(file) 
  IMPLICIT NONE 
  CHARACTER*(*) file 
  INTEGER unit 
  LOGICAL found 
  INQUIRE(FILE=file,EXIST=found) 
  IF(found) THEN 
     CALL filopn(unit,file,'old') 
     CALL filclo(unit,'delete') 
  END IF
END SUBROUTINE dlifex

                           

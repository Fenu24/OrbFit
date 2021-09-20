MODULE RDOPT_OF

IMPLICIT NONE

PUBLIC :: rdopto, rdopte, rdopti, rdoptf, rdodin

CONTAINS

! Copyright (C) 1997-1999 by Mario Carpino (carpino@brera.mi.astro.it)
! Version: February 16, 1999
! ---------------------------------------------------------------------
!
!  *****************************************************************
!  *                                                               *
!  *                         R D O P T O                           *
!  *                                                               *
!  *                   Read options for ORBFIT                     *
!  *                                                               *
!  *****************************************************************
!
! INPUT:    RUN       -  Run name
!           NIFX      -  Physical dimension of arrays ELFT,ELF1
!           CNIFX     -  Real name of NIFX parameter
!
! OUTPUT:   OP        -  Selection of execution steps:
!                          1 = initial orbit determination
!                          2 = differential correction
!                          3 = identifications
!                          4 = ephemerides
!                          5 = magnitude determination
!           NAME      -  Object names
!           NAMEO     -  Names to be used for looking in orb. el. files
!           NAMOF     -  Names of observation files (w/o extension)
!           DIR       -  Observations directory
!           NOBJ      -  Number of objects
!           ELFT      -  Input orbital element files (for all objects)
!           NELFT     -  Number of input element files (for all objects)
!           ELF1      -  Input orbital element files (for each object)
!           NELF1     -  Number of input element files (for each object)
!           OEFILE    -  Output orbital element file (possibly blank)
!           OEPTIM    -  Epoch of output elements (MJD, TDT)
!           OEPSET    -  Flag stating that an output epoch is requested
!           OETYPE    -  Type of output elements (CAR/EQU/KEP/EQP)
!           ERROR_MODEL - Error model file name
!
  SUBROUTINE rdopto(run,op,name,nameo,namof,dir,nobj,elft,nelft,    &
       &                  elf1,nelf1,oefile,oeptim,oepset,oetype,         &
       &                  nifx,cnifx,error_model)
    
    INCLUDE 'sysdep.h90'
    CHARACTER(LEN=*), INTENT(IN)  :: run 
    INTEGER,          INTENT(IN)  :: nifx
    CHARACTER(LEN=*), INTENT(IN)  :: cnifx
    CHARACTER*(80),   INTENT(OUT) :: name(2),nameo(2),namof(2) 
    CHARACTER(LEN=*), INTENT(OUT) :: dir(2), elft(nifx), elf1(nifx,2)
    INTEGER,          INTENT(OUT) :: op(5),nobj,nelft,nelf1(2)
    CHARACTER(LEN=*), INTENT(OUT) :: oefile, oetype
    DOUBLE PRECISION, INTENT(OUT) :: oeptim
    LOGICAL,          INTENT(OUT) :: oepset
    CHARACTER(LEN=*), INTENT(OUT) :: error_model 
! ***********************************************************************************
    INTEGER          :: i,lr,mjd,mjde,ln,ld
    DOUBLE PRECISION :: sec,sece
    LOGICAL          :: found,fail1,fail
    CHARACTER        :: cc*100,scale*10
    
    INTEGER lench
    EXTERNAL lench

! Required execution steps (by default, all are selected but ephemeris)
    fail=.FALSE.
    CALL rdnint('operations.','init_orbdet',op(1),.false.,            &
         &            found,fail1,fail)
    CALL rdnint('operations.','diffcor',op(2),.false.,                &
         &            found,fail1,fail)
    CALL rdnint('operations.','ident',op(3),.false.,                  &
         &            found,fail1,fail)
    CALL rdnint('operations.','ephem',op(4),.false.,                  &
         &            found,fail1,fail)
    CALL rdnint('operations.','magfit',op(5),.false.,                 &
         &            found,fail1,fail)

! Input orbital element file (common to all objects)
    nelft=0
    CALL rdmcha('input_files.','incond',elft,nelft,nifx,              &
         &            cnifx,.false.,found,fail1,fail)
    
! Output orbital element file
    lr=lench(run)
    oefile=run(1:lr)//'.oel'
    CALL rdncha('output_files.','elem',oefile,.false.,                &
         &            found,fail1,fail)
    
    nobj=1
! First object
    name(1)=run
    CALL rdncha('object1.','name',name(1),.false.,                    &
         &            found,fail1,fail)
! Observations directory
    IF(.NOT.fail1) THEN
       CALL rdncha('object1.','obs_dir',dir(1),.false.,              &
            &                found,fail1,fail)
    END IF
! Initial condition file
    nelf1(1)=0
    CALL rdmcha('object1.','inc_files',elf1(1,1),nelf1(1),nifx,       &
         &            cnifx,.false.,found,fail1,fail)
    IF(.NOT.found) CALL rdodin(name(1),nelf1(1),elf1(1,1),nifx)
    CALL rdncha('object1.','inc_name',nameo(1),.false.,               &
         &            found,fail1,fail)
    IF(.NOT.found) nameo(1)=name(1)
    CALL rmsp(nameo(1),ln)
    CALL rdncha('object1.','obs_fname',namof(1),.false.,              &
         &            found,fail1,fail)
    IF(.NOT.found) THEN
       namof(1)=name(1)
       CALL rmsp(namof(1),ln)
    END IF
      
! Second object
    name(2)=' '
    nameo(2)=' '
    namof(2)=' '
    CALL rdncha('object2.','name',name(2),.false.,                    &
         &            found,fail1,fail)
    IF(found) THEN
       nobj=2
       CALL rdncha('object2.','obs_dir',dir(2),.false.,              &
            &                found,fail1,fail)
       nelf1(2)=0
       CALL rdmcha('object2.','inc_files',elf1(1,2),nelf1(2),nifx,   &
            &                cnifx,.false.,found,fail1,fail)
       IF(.NOT.found) CALL rdodin(name(2),nelf1(2),elf1(1,2),nifx)
       CALL rdncha('object2.','inc_name',nameo(2),.false.,           &
            &                found,fail1,fail)
       IF(.NOT.found) nameo(2)=name(2)
       CALL rmsp(nameo(2),ln)
       CALL rdncha('object2.','obs_fname',namof(2),.false.,          &
            &                found,fail1,fail)
       IF(.NOT.found) THEN
          namof(2)=name(2)
          CALL rmsp(namof(2),ln)
       END IF
    END IF
      
! Epoch and type of output elements
    oeptim=0.d0
    CALL rdntim('output.','epoch',cc,mjd,sec,scale,.false.,           &
         &            oepset,fail1,fail)
    IF(oepset) THEN
       CALL cnvtim(mjd,sec,scale,mjde,sece,'TDT')
       oeptim=mjde+sece/86400.d0
    END IF
    
    CALL rdncha('output.','elements',oetype,.false.,                  &
         &            found,fail1,fail)
    
    IF(fail) STOP '**** rdopto: abnormal end ****'
    
    DO 2 i=1,nobj
       ld=lench(dir(i))
       IF(ld.LE.0) GOTO 2
       IF(dir(i)(ld:ld).NE.dircha) dir(i)(ld+1:ld+1)=dircha
2   END DO
    
! Error model
    CALL rdncha('error_model.','name',error_model,.false.,found,fail1,fail)
    IF(fail) STOP '**** rdopto: abnormal end ****'
    
  END SUBROUTINE rdopto

! Copyright (C) 1997 by Mario Carpino (carpino@brera.mi.astro.it)
! Version: December 15, 1997
! ---------------------------------------------------------------------
!
!  *****************************************************************
!  *                                                               *
!  *                         R D O P T E                           *
!  *                                                               *
!  *         Read options for ephemeris generation (ORBFIT)        *
!  *                                                               *
!  *****************************************************************
!
  SUBROUTINE rdopte

! Common blocks to be initialized:
    INCLUDE 'comeph.h90'
    
    INTEGER           :: mjd,mjde
    DOUBLE PRECISION  :: sec,sece
    LOGICAL           :: found,fail1,fail
    CHARACTER         :: tmp*100,tsc*3
    
    fail=.false.

! List of objects
    CALL rdmint('ephem.','objects',kepobj,nepobj,3,'3',.false.,       &
         &            found,fail1,fail)
    IF(.NOT.found) THEN
       nepobj=3
       kepobj(1)=1
       kepobj(2)=2
       kepobj(3)=3
    END IF

! Limits and stepsize of ephemeris
    CALL rdntim('ephem.epoch.','start',tmp,mjd,sec,tsc,.true.,        &
         &            found,fail1,fail)
    IF(.NOT.fail1) THEN
       CALL cnvtim(mjd,sec,tsc,mjde,sece,'TDT')
       teph1=mjde+sece/86400.d0
    END IF
    CALL rdntim('ephem.epoch.','end',tmp,mjd,sec,tsc,.true.,          &
         &            found,fail1,fail)
    IF(.NOT.fail1) THEN
       CALL cnvtim(mjd,sec,tsc,mjde,sece,'TDT')
       teph2=mjde+sece/86400.d0
    END IF
    dteph=1
    CALL rdnrea('ephem.','step',dteph,.false.,found,fail1,fail)
    
! Observatory code
    idsta=500
    CALL rdnint('ephem.','obscode',idsta,.false.,found,fail1,fail)
    
! Timescale
    ephtsc='TDT'
    CALL rdncha('ephem.','timescale',ephtsc,.false.,found,fail1,fail)
    
! Output fields
!    ephfld='cal,coord,delta,r,elong,phase,mag'
    ephfld='cal,mjd,coord,mag,elev,airm,elsun,elong,mooel,glat,glon,r,delta,appmot,skyerr'
    CALL rdncha('ephem.','fields',ephfld,.false.,found,fail1,fail)
    
    IF(fail) STOP '**** rdopte: abnormal end ****'
    
    iiceph=36
    
  END SUBROUTINE rdopte

! Copyright (C) 1997-1998 by Mario Carpino (carpino@brera.mi.astro.it)
! Version: December 7, 1998
! ---------------------------------------------------------------------
!
!  *****************************************************************
!  *                                                               *
!  *                         R D O P T I                           *
!  *                                                               *
!  *        Read options for orbit identification (ORBFIT)         *
!  *                                                               *
!  *****************************************************************
!
  SUBROUTINE rdopti
    IMPLICIT NONE
    
! Common blocks to be initialized:
    INCLUDE 'comidn.h90'
    
    LOGICAL :: found,fail1,fail
    
    fail=.false.
    
    amfit=.true.
    CALL rdnlog('ident.','aM_fit',amfit,.false.,found,                &
         &            fail1,fail)

    IF(fail) STOP '**** rdopti: abnormal end ****'
    
    iicidn=36
    
  END SUBROUTINE rdopti

! Copyright (C) 1997-1998 by Mario Carpino (carpino@brera.mi.astro.it)
! Version: December 14, 1998
! ---------------------------------------------------------------------
!
!  *****************************************************************
!  *                                                               *
!  *                         R D O P T F                           *
!  *                                                               *
!  *      Read options for least squares orbital fit(ORBFIT)       *
!  *                                                               *
!  *****************************************************************
!
  SUBROUTINE rdoptf
    IMPLICIT NONE

! Common blocks to be initialized:
    INCLUDE 'comlsf.h90'
    
    LOGICAL :: found,fail1,fail
    
    fail=.false.

    IF(fail) STOP '**** rdoptf: abnormal end ****'
    
    iiclsf=36
    
  END SUBROUTINE rdoptf

! Copyright (C) 1997-1998 by Mario Carpino (carpino@brera.mi.astro.it)
! Version: December 14, 1998
! ---------------------------------------------------------------------
!
!  *****************************************************************
!  *                                                               *
!  *                         R D O D I N                           *
!  *                                                               *
!  *      Default names for input orbital element files (ORBFIT)   *
!  *                                                               *
!  *****************************************************************
!
! 
  SUBROUTINE rdodin(name,nelf,elf,nifx)
    IMPLICIT NONE
    
    INTEGER           :: nelf,nifx
    CHARACTER(LEN=*)  :: name,elf(nifx)
    
    INTEGER           :: ln
    CHARACTER(LEN=80) :: file
    LOGICAL           :: found
    
    INTEGER lench
    EXTERNAL lench
    
    IF(nelf.NE.0) RETURN
    
    ln=lench(name)
    file=name(1:ln)//'.ele'
    INQUIRE(FILE=file,EXIST=found)
    IF(found .AND. nelf.LT.nifx) THEN
       nelf=nelf+1
       elf(nelf)=file
    END IF
    
    file='astorb.dat'
    file=name(1:ln)//'.ele'
    INQUIRE(FILE=file,EXIST=found)
    IF(found .AND. nelf.LT.nifx) THEN
       nelf=nelf+1
       elf(nelf)=file
    END IF
    
  END SUBROUTINE rdodin
  
END MODULE RDOPT_OF

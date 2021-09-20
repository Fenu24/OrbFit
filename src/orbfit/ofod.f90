! Module OFOD
! Contains all the subroutines needed by orbfit:
! - ofinip
! - ofiobs
! - ofiorb
! - ofinod
! - ofofit
! - ofiden
! - ofephe

MODULE OFOD
  USE fund_const
  USE output_control
  USE astrometric_observations
  USE orbit_elements
  USE least_squares
  USE propag_state
  USE dyn_param
  IMPLICIT NONE
  PRIVATE
  
  PUBLIC :: ofinip, ofiobs, ofiorb, ofinod, ofofit, ofiden, ofprop, ofephe

! NEEDED common blocks:
  INCLUDE 'parobx.h90'
  INCLUDE 'parnob.h90'
  INCLUDE 'parcmc.h90'
  INCLUDE 'comlsf.h90'
  INCLUDE 'comidn.h90'
  INCLUDE 'comeph.h90'
CONTAINS

! Copyright (C) 1999 by Mario Carpino (carpino@brera.mi.astro.it)
! Version: February 8, 1999
! ---------------------------------------------------------------------
!
!  *****************************************************************
!  *                                                               *
!  *                         O F I N I P                           *
!  *                                                               *
!  *       Initializations required for orbit propagation          *
!  *                                                               *
!  *****************************************************************
!
! INPUT:    RUN       -  Run name
!
SUBROUTINE ofinip(run)
  
  CHARACTER(LEN=*), INTENT(IN) :: run
!*********************************************  
  INTEGER      :: lf
  CHARACTER*80 :: file
  
  INTEGER lench
  EXTERNAL lench
  
  lf=lench(run)
  
! Output file for errors
  file=run(1:lf)//'.err'
  CALL filopn(ierrou,file,'UNKNOWN')
  numerr=0

! Output file for close approaches
  file=run(1:lf)//'.clo'
  CALL filopn(iuncla,file,'UNKNOWN')
  numcla=0

! Output file for propagator parameters
  file=run(1:lf)//'.pro'
  CALL filopn(ipirip,file,'UNKNOWN')
  
END SUBROUTINE ofinip

! Copyright (C) 1997-1999 by Mario Carpino (carpino@brera.mi.astro.it)
! Version: February 11, 1999
! Version: November 4, 1999 (new inobs call, SRC)
! ---------------------------------------------------------------------
!
!  *****************************************************************
!  *                                                               *
!  *                         O F I O B S                           *
!  *                                                               *
!  *                Auxiliary routine for ORBFIT:                  *
!  *                    input of observations                      *
!  *                                                               *
!  *****************************************************************
!
! INPUT:    UNIREP    -  FORTRAN unit for report
!           NAME      -  Object file names
!           NAMOF     -  Names of observation files (w/o extension)
!           DIR       -  Directory containing observation/residual files
!           NOBJ      -  Number of objects
!           NOBX      -  Max number of observations
!
! OUTPUT:   N         -  Number of observations for each object
!           NT        -  Total number of observations
!           IP1       -  Pointer to first observation for each object
!           IP2       -  Pointer to first weight for each object
!           OBS,OBSW  -  Observations
!
SUBROUTINE ofiobs(unirep,name,namof,dir,nobj,n,nt,     &
     &            ip1,ip2,obs,obsw,error_model)
  

  
  INTEGER,                         INTENT(IN)  :: unirep,nobj
  CHARACTER*(*),                   INTENT(IN)  :: name(nobj),namof(nobj),dir(nobj)
  CHARACTER(LEN=*),                INTENT(IN)  :: error_model
  INTEGER,                         INTENT(OUT) :: nt,n(nobj), ip1(nobj), ip2(nobj)
  TYPE(ast_obs),  DIMENSION(nobx), INTENT(OUT) :: obs
  TYPE(ast_wbsr), DIMENSION(nobx), INTENT(OUT) :: obsw
!************************************************************************
  INTEGER :: obscod(nobx),sel(nobx),iobs(nobx)
  LOGICAL :: precob,ok, change
  INTEGER :: i,ln,ld
  
  INTEGER lench
  EXTERNAL lench
  
  WRITE(unirep,120)
120 FORMAT('Input of observations:')
  
  nt=0
  precob=.false.
! verbosity levels for an interactive program
  verb_io=10
  
  DO 1 i=1,nobj
     ip1(i)=nt+1
     ip2(i)=2*nt+1
! input data
     CALL input_obs(dir(i),namof(i),precob,error_model,ok,obs(ip1(i):),obsw(ip1(i):),n(i),unirep,change)
     IF(.NOT.ok) STOP '**** OFIOBS : abnormal end ****'
     nt=nt+n(i)
     ld=lench(dir(i))
     ln=lench(name(i))
     WRITE(unirep,200) n(i),name(i)(1:ln)
200  FORMAT(I9,' observations read for object ',A)
1 END DO

END SUBROUTINE ofiobs

! Copyright (C) 1997-2000 by Mario Carpino (carpino@brera.mi.astro.it)
! Version: June 7, 2000
! ---------------------------------------------------------------------
!
!  *****************************************************************
!  *                                                               *
!  *                         O F I O R B                           *
!  *                                                               *
!  *                Auxiliary routine for ORBFIT:                  *
!  *                  input of orbital elements                    *
!  *                                                               *
!  *****************************************************************
!
! INPUT:    UNIREP    -  FORTRAN unit for report
!           ELFT      -  Input file names (for all objects)
!           NELFT     -  Number of input files (for all objects)
!           ELF1      -  Input file names (for each object)
!           NELF1     -  Number of input files (for each object)
!           NAME      -  Object names (true)
!           NAMEO     -  Names to be searched in orbital element files
!           NOBJ      -  Number of objects
!           NFIX      -  First physical dimension of array ELF1
!
! OUTPUT:   ELEM      -  Orbital elements
!           ELEM_UNC  -  Orbital element uncertainty
!           DEFORB    -  Orbit definition flag
!           DEFCN     -  Tells whether covariance/normal matrices are defined
!           MASS      -  Mass
!           COMELE    -  Comment on orbital elements
!           ND        -  Dimension of the parameter space
!
SUBROUTINE ofiorb(unirep,elft,nelft,elf1,nelf1,name,nameo,nobj,elem,elem_unc,deforb,defcn,mass,comele,nfix,nd)
  

  
  INTEGER,                           INTENT(IN)  :: unirep,nobj,nelft,nelf1(nobj),nfix
  CHARACTER(LEN=*),                  INTENT(IN)  :: elft(nfix),elf1(nfix,nobj)
  CHARACTER(LEN=*),                  INTENT(IN)  :: name(nobj),nameo(nobj)
  TYPE(orbit_elem),DIMENSION(nobjx), INTENT(OUT) :: elem
  TYPE(orb_uncert),DIMENSION(nobjx), INTENT(OUT) :: elem_unc
  LOGICAL,                           INTENT(OUT) :: deforb(nobjx),defcn(nobjx)
  DOUBLE PRECISION,                  INTENT(OUT) :: mass(nobjx)
  CHARACTER(LEN=*),                  INTENT(OUT) :: comele(nobj)
  INTEGER,                           INTENT(OUT) :: nd 
!********************************************************************************************************     
  DOUBLE PRECISION :: telem(nobj1x),elemv(6,nobj1x)
  DOUBLE PRECISION :: h(nobj1x),g(nobj1x)
  DOUBLE PRECISION :: cove(ndimx,ndimx,nobj1x),nore(ndimx,ndimx,nobj1x)
  CHARACTER*3      :: eltype(nobj1x)
  INTEGER          :: i,j,ln,lc,lno
  
  INTEGER lench
  EXTERNAL lench
  
  WRITE(unirep,122)
122 FORMAT('Input of orbital elements:')
  
! No orbits defined yet
  DO 1 i=1,nobj
     deforb(i)=.false.
     defcn(i)=.false.
     mass(i)=0.d0
1 END DO

! Input for particular files for each object
  DO 2 i=1,nobj
     CALL rdelem(unirep,nameo(i),1,elf1(1,i),nelf1(i),deforb(i),       &
          &            defcn(i),eltype(i),telem(i),elemv(1,i),               &
          &            cove(1,1,i),nore(1,1,i),mass(i),h(i),g(i),comele(i))
     IF(.not.deforb(i))THEN 
        WRITE(*,*)'!!WARNING! WARNING! WARNING! WARNING!!' 
        WRITE(*,*)'Asteroid ',nameo(i),' not found in asteroid catalogs.'
        WRITE(*,*)nameo(i),'  not in catalogs.' 
        RETURN
     ENDIF
     IF(eltype(i).EQ.'KEP') THEN 
        WRITE(*,100) elemv(1,i),elemv(2,i),(elemv(j,i)*degrad,j=3,6)   
100     FORMAT(8X,'Semimajor axis     =',1P,E24.16,0P,' au'/          &
             &       8X,'Eccentricity       =',F20.15/                    &
             &       8X,'Inclination        =',F18.13,' deg'/             &
             &       8X,'Long. of node      =',F18.13,' deg'/             &
             &       8X,'Arg. of pericenter =',F18.13,' deg'/             &
             &       8X,'Mean anomaly       =',F18.13,' deg')  
     ELSEIF(eltype(i).EQ.'COM') THEN
        WRITE(*,104) elemv(1,i),elemv(2,i),(elemv(j,i)*degrad,j=3,5),elemv(6,i) 
104     FORMAT(8X,'Pericenter distance  =',1P,E23.14,0P,' au'/        &
             &       8X,'Eccentricity       =',F20.15/                    &
             &       8X,'Inclination        =',F18.13,' deg'/             &
             &       8X,'Long. of node      =',F18.13,' deg'/             &
             &       8X,'Arg. of pericenter =',F18.13,' deg'/             &
             &       8X,'Time of pericenter =',F18.13,' MJD')  
     ELSEIF(eltype(i).EQ.'EQU') THEN 
        WRITE(*,101) (elemv(j,i),j=1,5),elemv(6,i)*degrad
101     FORMAT(8X,'Semimajor axis     =',1P,E23.14,0P,' au'/     &
             &       8X,'h [e*sin(w)]       =',F20.15/                   &
             &       8X,'k [e*cos(w)]       =',F20.15/                   &
             &       8X,'P [tg(i/2)*sin(N)] =',F20.15/                   &
             &       8X,'Q [tg(i/2)*cos(N)] =',F20.15/                   &
             &       8X,'Mean longitude     =',F18.13,' deg')   
     ELSEIF(eltype(i).EQ.'CAR') THEN 
        WRITE(*,102) elemv(:,i) 
102     FORMAT(8X,'Position vector   =',1X,1P,3E22.14,' au'/              &
             &       8X,'Velocity vector   =',1X,3E22.14,' au/d') 
     END IF
     IF(dyn%nmod.GT.0.AND..NOT.ngr_opt)THEN
        WRITE(*,145) dyn%dp(1:dyn%ndp)
145     FORMAT('Non-gravitational perturbations: '/   &
             & 8X, 'A/M  = ',1X,1P,E22.14,' m^2/ton'/ &
             & 8X, ' A2  = ',1X,1P,E22.14,' 10^(-10) au/d^2')
     END IF
2 END DO

! Input from common files
  CALL rdelem(unirep,nameo,nobj,elft,nelft,deforb,defcn,eltype,telem,elemv,cove,nore,mass,h,g,comele)

  nd=6+nls
  DO i=1,nobj
     elem(i)=undefined_orbit_elem
     IF(deforb(i)) THEN
        elem(i)%coo=eltype(i)
        elem(i)%coord=elemv(1:6,i)
        elem(i)%t=telem(i)
        elem(i)%mag_set=(h(i) > 0.d0)
        elem(i)%h_mag=h(i)
        elem(i)%g_mag=g(i)
        ln=lench(name(i))
        IF(name(i).EQ.nameo(i)) THEN
           WRITE(unirep,123) name(i)(1:ln)
        ELSE
           lno=lench(nameo(i))
           WRITE(unirep,125) name(i)(1:ln),nameo(i)(1:lno)
        END IF
        lc=lench(comele(i))
        IF(lc.GT.0) WRITE(unirep,124) comele(i)(1:lc)
        IF(dyn%nmod.GT.0.AND..NOT.ngr_opt)THEN
           CALL outele_ngr(unirep,elemv(1,i),eltype(i),telem(i),' ',.true.,.false.,dyn%nmod,dyn%ndp,dyn%dp,ngr_opt)
        ELSE
           CALL outele(unirep,elemv(1,i),eltype(i),telem(i),' ',.true.,.false.)
        END IF
     END IF
     CALL undefined_orb_uncert(nd,elem_unc(i))
     IF(defcn(i)) THEN
        elem_unc(i)%c(1:nd,1:nd)=nore(1:nd,1:nd,i)
        elem_unc(i)%g(1:nd,1:nd)=cove(1:nd,1:nd,i)
        elem_unc(i)%succ=.true.
        elem_unc(i)%ndim=nd
     END IF
  END DO
123 FORMAT(5X,'Orbital elements for object ',A,':')
125 FORMAT(5X,'Orbital elements for object ',A,' (with name ',A,'):')
124 FORMAT(5X,'(origin: ',A,')')
  
END SUBROUTINE ofiorb

!
!  *****************************************************************
!  *                                                               *
!  *                         O F I N O D                           *
!  *                                                               *
!  *                Auxiliary routine for ORBFIT:                  *
!  *                 initial orbit determination                   *
!  *                                                               *
!  *****************************************************************
!
! INPUT:    UNIREP      -  FORTRAN unit for report
!           ELEOUT      -  Name of file for output of orbital elements
!           OPT         -  Initial orbit determination option
!           NAME        -  Object names
!           NAMOF       -  Names of observation files (w/o extension)
!           DEFORB      -  Orbit definition flag
!           DEFCN       -  Tells whether covariance/normal matrices
!                              are defined
!           DIR         -  Directory containing observation/residual file
!           NOBJ        -  Number of objects
!           OBS, OBSW   -  OBSERVATIONS
!           N           -  Number of observations for each object
!           NT          -  Total number of observations
!           IP1         -  Pointer to first observation for each object
!           ERROR_MODEL -  Error model
!           ND          -  Dimension of the parameter space (MUST BE 6)
!
! IN/OUT:   ELEM        -  Orbital elements
!           COMELE      -  Comment on orbital elements
!
! Input variables DEFORB, DEFCN, SEL are modified by the routine
!
SUBROUTINE ofinod(unirep,eleout,opt,name,namof,deforb,defcn,dir,nobj,obs,obsw,n,nt,ip1,elem,comele,error_model,nd)
  
  INTEGER                                        :: unirep,opt,nobj,ip1(nobj)
! OBSERVATIONS
  INTEGER                                        :: nt,n(nobj)
! Data types
  TYPE(ast_obs),DIMENSION(nt)                    :: obs
  TYPE(ast_wbsr),DIMENSION(nt)                   :: obsw
  CHARACTER(LEN=*)                               :: error_model ! weighing model
  
  TYPE(orbit_elem), DIMENSION(nobj), INTENT(OUT) :: elem
  CHARACTER(LEN=*),                  INTENT(OUT) :: comele(nobj)
  LOGICAL                                        :: deforb(nobj),defcn(nobj)
  CHARACTER(LEN=*)                               :: eleout,name(nobj),namof(nobj),dir(nobj)
  INTEGER,                           INTENT(IN)  :: nd
!*********************************************************************************************
  INTEGER                                        :: i,ln,uniele,lc,ld
  CHARACTER(LEN=150)                             :: rwofil
  LOGICAL                                        :: opnd,doit,fail
 
  INTEGER lench
  EXTERNAL lench

  IF(nd.NE.6) THEN
     WRITE(*,*) 'OFINOD: Let me explain why you are wrong: nd is', nd,', while it must be 6'
     STOP
  END IF
  
  opnd=.false.
  
  DO 1 i=1,nobj
     doit=.false.
     IF(opt.EQ.1) THEN
        doit=(.NOT.deforb(i))
     ELSEIF(opt.EQ.2) THEN
        doit=.true.
     END IF
     IF(.NOT.doit) GOTO 1
     ln=lench(namof(i))
     ld=lench(dir(i))
     IF(ld.GT.0) THEN
        rwofil=dir(i)(1:ld)//namof(i)(1:ln)//'.rwo'
     ELSE
        rwofil=namof(i)(1:ln)//'.rwo'
     END IF
     CALL io_det(unirep,rwofil,name(i),obs(ip1(i):ip1(i)+n(i)-1),  &
          & obsw(ip1(i):ip1(i)+n(i)-1),n(i),error_model,2,         &
          & elem(i)%coord,elem(i)%t,elem(i)%coo,comele(i),fail)
     deforb(i)=(.NOT.fail)
     defcn(i)=.false.
     IF(.NOT.fail) THEN
        IF(.NOT.opnd) THEN
           CALL filopn(uniele,eleout,'UNKNOWN')
           CALL wromlh(uniele,'ECLM','J2000')
           opnd=.true.
        END IF
        lc=lench(comele(i))
        WRITE(uniele,300) comcha,comele(i)(1:lc)
        CALL write_elems(elem(i),name(i),UNIT=uniele)
     END IF
300  FORMAT(A,1X,A)
     
1 END DO
  
  IF(opnd) CALL filclo(uniele,' ')
  
END SUBROUTINE ofinod


!
!  *****************************************************************
!  *                                                               *
!  *                         O F O F I T                           *
!  *                                                               *
!  *                Auxiliary routine for ORBFIT:                  *
!  *                  least squares orbital fit                    *
!  *                                                               *
!  *****************************************************************
!
! IN/OUT:   UNIREP    -  FORTRAN unit for report
!           UNIELE    -  FORTRAN unit for output of orbital elements
!           UNIDIF    -  FORTRAN unit for logfilr of DIFCOR
!           OPELE     -  Flag (need to open UNIELE)
!           ELEOUT    -  Name of file for output of orbital elements
!           ORBOPT    -  Orbit determination option
!           MAGOPT    -  Magnitude determination option
!           NAME      -  Object names
!           NAMOF     -  Names of observation files (w/o extension)
!           DEFORB    -  Orbit definition flag
!           DEFCN     -  Tells whether covariance/normal matrices
!                            are defined
!           ELEM      -  Orbital elements
!           ELEM_UNC  -  Orbital element uncertainty
!           COMELE    -  Comment on orbital elements
!           DIR       -  Directory containing observation/residual file
!           NOBJ      -  Number of objects
!           OBS, OBSW -  Observations
!           N         -  Number of observations for each object
!           NT        -  Total number of observations
!           IP1       -  Pointer to first observation for each object
!           IP2       -  Pointer to first weight for each object
!           ERROR_MODEL - Error model
!           OEPTIM    -  Epoch of output elements (MJD, TDT)
!           OEPSET    -  Flag stating that an output epoch is requested
!           OETYPE    -  Type of output elements (CAR/EQU/KEP/EQP)
!           ND        -  Dimension of the parameter space (INTENT IN)
!
SUBROUTINE ofofit(unirep,uniele,unidif,opele,eleout,orbopt,magopt,name,namof,deforb,defcn,elem,elem_unc, &
     &                  comele,dir,nobj,obs,obsw,n,nt,ip1,ip2,error_model,oeptim,oepset,oetype,nd)



! objects
  INTEGER                          :: unirep,uniele,unidif
  INTEGER                          :: orbopt,magopt,nobj,nt
  INTEGER                          :: ip1(nobj),ip2(nobj),sel(nt),iobs(nt),obscod(nt)
  DOUBLE PRECISION                 :: oeptim
  LOGICAL                          :: deforb(nobj),defcn(nobj),opele,oepset
  INTEGER                          :: fail
  CHARACTER(LEN=*)                 :: name(nobj),namof(nobj),eleout,dir(nobj)
  CHARACTER(LEN=*)                 :: comele(nobj),oetype
! OBSERVATIONS
  INTEGER                          :: n(nobj),nfer
! new data types
  TYPE(ast_obs),DIMENSION(nt)      :: obs
  TYPE(ast_wbsr),DIMENSION(nt)     :: obsw
  TYPE(orbit_elem),DIMENSION(nobj) :: elem
  TYPE(orb_uncert),DIMENSION(nobj) :: elem_unc
! error model 
  CHARACTER(LEN=*)                 :: error_model 
! Non-gravitational perturbations
  INTEGER, INTENT(IN)              :: nd
!********************************************************************************************  

! TUNING PARAMETERS
! Max RMS of a "good" observation (arcsec)
  DOUBLE PRECISION      ::maxrms
  PARAMETER (maxrms=5.d0)
! Small time (d)
  DOUBLE PRECISION      ::epst
  PARAMETER (epst=1.0D-9)
  
  INTEGER               :: i,ln,icor(6),j1,j2,ld
  DOUBLE PRECISION      :: gma1,elemt(6),enne,csinor,delnor,rms1,t1,t2
  DOUBLE PRECISION      :: dxde(6,6)
  DOUBLE PRECISION      :: rmsh
  DOUBLE PRECISION      :: hnew
  CHARACTER             :: file*150
  LOGICAL               :: ok,chep,error,defcov,defnor
  TYPE(orbit_elem)      :: elem0,elem1,elem_print
  TYPE(orb_uncert)      :: elem_unc_print
  
  INTEGER lench
  EXTERNAL lench
  
  IF(iiclsf.NE.36) STOP '**** ofofit: internal error (01) ****'
  rms1=maxrms*radsec
  
  DO 10 i=1,6
     icor(i)=1
10 END DO
  
  DO 1 i=1,nobj
     !IF(.NOT.deforb(i)) THEN
     !   WRITE(*,*) 'OFOFIT Warning: Orbit not defined for object ', i
     !   GOTO 1
     !END IF
! Set error model
     CALL errmod_set(error_model)
! Elements
     elem0=elem(i)
     ln=lench(name(i))
     WRITE(unirep,206) name(i)(1:ln)
206  FORMAT('Differential correction for object ',A,':')
     WRITE(unirep,123)
123  FORMAT(5X,'Starting orbital elements:')
     IF(dyn%nmod.GT.0.AND..NOT.ngr_opt)THEN
        CALL outele_ngr(unirep,elem0%coord,elem0%coo,elem0%t,' ',.true.,.false.,dyn%nmod,dyn%ndp,dyn%dp,ngr_opt)
     ELSE
        CALL outele(unirep,elem0%coord,elem0%coo,elem0%t,' ',.true.,.false.)
     END IF
     IF(magopt.NE.0 .AND. elem(i)%mag_set) THEN
        IF(elem0%mag_set .AND. (elem0%h_mag > -100.D0)) WRITE(unirep,140) elem0%h_mag
     END IF
140  FORMAT(10X,'Absolute mag. H    =',F9.2)
! Differential corrections     
     CALL diff_cor(n(i),obs(ip1(i):ip1(i)+n(i)-1),obsw(ip1(i):ip1(i)+n(i)-1), &
          & elem0,icor,unidif,elem1,elem_unc(i),csinor,delnor,ok,nd)
     IF(ok) THEN
        elem(i)=elem1
        defcn(i)=.true.
        comele(i)='computed with least squares fit'
! Magnitude estimation
        IF(magopt.NE.0) THEN
           CALL mag_est(n(i),obs(ip1(i):ip1(i)+n(i)-1),obsw(ip1(i):ip1(i)+n(i)-1),hnew,rmsh)
           IF(rmsh.GT.0.d0) elem(i)%h_mag=hnew
        ELSE
           rmsh=0.d0
        END IF
        IF(opele) THEN
           OPEN(uniele,FILE=eleout,STATUS='UNKNOWN')
           opele=.false.
           WRITE(uniele,301) comcha
           CALL wromlh(uniele,'ECLM','J2000')
        END IF
! Change coordinates (if needed)
        CALL coo_cha(elem(i),oetype,elem_print,fail,dxde)
        CALL convertunc(elem_unc(i),dxde,elem_unc_print)
        WRITE(unirep,124)
        IF(dyn%nmod.GT.0.AND..NOT.ngr_opt)THEN
           CALL outele_ngr(unirep,elem_print%coord,elem_print%coo,elem_print%t,' ',.true.,.false.,dyn%nmod,dyn%ndp,dyn%dp,ngr_opt)
        ELSE
           CALL outele(unirep,elem_print%coord,elem_print%coo,elem_print%t,' ',.true.,.false.)
        END IF
        IF(magopt.NE.0 .AND. elem(i)%mag_set .AND. (elem(i)%h_mag > -100.D0)) WRITE(unirep,140) elem(i)%h_mag
        WRITE(unirep,230) csinor,delnor
! Write elements
        IF(dyn%nmod.gt.0)THEN
           CALL write_elems(elem_print,name(i),DYNAM_PAR=dyn,NLSLOC=nls,LSLOC=ls,UNC=elem_unc_print,UNIT=uniele)
        ELSE
           CALL write_elems(elem_print,name(i),UNC=elem_unc_print,UNIT=uniele)
        END IF
        ld=lench(dir(i))
        ln=lench(namof(i))
        IF(ld.GT.0) THEN
           file=dir(i)(1:ld)//namof(i)(1:ln)//'.rwo'
        ELSE
           file=namof(i)(1:ln)//'.rwo'
        END IF
! Write rwo        
        CALL write_rwo(file,obs(ip1(i):ip1(i)+n(i)-1),obsw(ip1(i):ip1(i)+n(i)-1),n(i),error_model,csinor,rmsh)
     ELSE
        WRITE(unirep,202)
     END IF
     
1 END DO
124 FORMAT(5X,'Corrected orbital elements:')
301 FORMAT(A,' Orbits computed with differential corrections')
202 FORMAT(5X,'FAILED')
230 FORMAT(5X,'Residual norm   =',1P,E11.3/                           &
         &       5X,'Correction norm =',1P,E11.3)
  
END SUBROUTINE ofofit

!
!  *****************************************************************
!  *                                                               *
!  *                         O F I D E N                           *
!  *                                                               *
!  *                Auxiliary routine for ORBFIT:                  *
!  *                     orbit identification                      *
!  *                                                               *
!  *****************************************************************
!
! IN/OUT:   UNIREP    -  FORTRAN unit for report
!           UNIELE    -  FORTRAN unit for output of orbital elements
!           UNIDIF    -  FORTRAN unit for logfilr of DIFCOR
!           OPELE     -  Flag (need to open UNIELE)
!           ELEOUT    -  Name of file for output of orbital elements
!           MAGOPT    -  Magnitude determination option
!           NAME      -  Object names
!           NAMOF     -  Names of observation files (w/o extension)
!           DEFORB    -  Orbit definition flag
!           DEFCN     -  Tells whether covariance/normal matrices
!                            are defined
!           ELEM      -  Orbital elements
!           ELEM_UNC  -  Orbital element uncertainty
!           MASS      -  Object masses
!           COMELE    -  Comment on orbital elements
!           DIR       -  Directory containing observation/residual file
!           NOBJ      -  Number of objects
!           OBS, OBSW -  Observations
!           N         -  Number of observations for each object
!           NT        -  Total number of observations
!           IP1       -  Pointer to first observation for each object
!           IP2       -  Pointer to first weight for each object
!           GMSUN     -  G*Mass(Sun)
!           OETYPE    -  Type of output elements (CAR/EQU/KEP/EQP)
!           ND        -  Dimension of the parameter space
!
SUBROUTINE ofiden(unirep,uniele,unidif,opele,eleout,magopt,name,namof,deforb,defcn,elem,elem_unc,         &
     &                  mass,comele,dir,nobj,obs,obsw,n,nt,ip1,ip2,gmsun,oetype,error_model,nd)

  
  INTEGER                            :: unirep,uniele,unidif,nobj,magopt
  INTEGER                            :: ip1(nobjx),ip2(nobjx)
! OBSERVATIONS
  INTEGER                            :: nt,n(nobj)
! new data types
  TYPE(ast_obs),DIMENSION(nt)        :: obs
  TYPE(ast_wbsr),DIMENSION(nt)       :: obsw
  TYPE(orbit_elem),DIMENSION(nobj1x) :: elem
  TYPE(orb_uncert),DIMENSION(nobj1x) :: elem_unc
  
  CHARACTER(LEN=*)                   :: error_model ! weighing model
  
  DOUBLE PRECISION                   :: gmsun
  DOUBLE PRECISION                   :: mass(nobj1x)
  LOGICAL                            :: opele,deforb(nobj1x),defcn(nobj1x)
  CHARACTER(LEN=*)                   :: eleout,name(nobj1x),namof(nobj1x)
  CHARACTER(LEN=*)                   :: dir(nobj1x),comele(nobj1x),oetype
  INTEGER, INTENT(IN)                :: nd
!****************************************************************************  
  INTEGER                            :: i,k,ng,lnt,ln1,icor(6),nsm,fail_flag
  TYPE(orbit_elem),DIMENSION(nobjx)  :: elemp
  DOUBLE PRECISION                   :: hnew,rmsh
  DOUBLE PRECISION                   :: csinor,delnor,elemo(6),telemo
  DOUBLE PRECISION                   :: telemt,elemt(6),masst,ht,gt,ennet,gma1,enne
  DOUBLE PRECISION                   :: gmagc
  TYPE(orbit_elem)                   :: elem1
  CHARACTER                          :: namet*100,titnam*80,file*150
  LOGICAL                            :: ok,autrep

! TUNING PARAMETERS
! Small time (d)
  DOUBLE PRECISION epst
  PARAMETER (epst=1.0D-9)
  
  INTEGER lench
  EXTERNAL lench

  IF(nd.GT.6)THEN
     WRITE(*,*) 'OFIDEN: wrong dimension nd is ',nd,'while it should be 6'
     STOP
  END IF
  
  IF(nobj.LT.2) RETURN
!  DO i=1,nobj
!     IF(.NOT.(deforb(i).AND.defcn(i))) RETURN
!  END DO
  
! Set error model
  CALL errmod_set(error_model)
  nobj=2
  IF(nobj.GT.nobjx) STOP '**** ofiden: internal error (03) ****'
  IF(nobj1x.LT.3) STOP '**** ofiden: internal error (04) ****'
  
  WRITE(unirep,207)
207 FORMAT('Orbit identification:')
  
  DO i=1,nobj
     gma1=gmsun*(1+mass(i))
     CALL coo_cha(elem(i),'EQU',elemp(i),fail_flag)
     ln1=lench(name(i))
     WRITE(unirep,200) name(i)(1:ln1)
     CALL outele(unirep,elemp(i)%coord,'EQU',elemp(i)%t,' ',.true.,.false.)
     IF(magopt.NE.0 .AND. elem(i)%mag_set .AND. elem(i)%h_mag.GT.-100.D0) THEN
        WRITE(unirep,140) elem(i)%h_mag
     END IF
  END DO
200 FORMAT(5X,'Starting orbital elements for object ',A,':')
140 FORMAT(10X,'Absolute mag. H    =',F9.2)
  
! Computation of starting elements (average of single objects)
! First guess of a unique solution for both arcs, with averages
! and fit to longitudes
  telemt=0
  masst=0
  DO k=2,5
     elemt(k)=0
  END DO
  elemt(1)=elemp(1)%coord(1)
  elemt(6)=elemp(1)%coord(6)
  ht=0
  gt=0
  nsm=0
  ng=0
  ennet=0
! can be done                                                           
  WRITE(*,*) ' GUESS FROM LONGITUDE FIT='
  DO i=1,nobj
     telemt=(telemt+elemp(i)%t)
     masst=masst+mass(i)
     DO k=2,5
        elemt(k)=elemt(k)+elemp(i)%coord(k)
     END DO
     IF(elem(i)%mag_set) THEN
        IF(elem(i)%h_mag.GT.-100.d0) THEN
           ht=ht+elem(i)%h_mag
           gt=gt+elem(i)%g_mag
           nsm=nsm+1
        END IF
     END IF
  END DO
  telemt=telemt/nobj
  WRITE(*,*) telemt
  masst=masst/nobj
  DO k=2,5
     elemt(k)=elemt(k)/nobj
  END DO
  IF(nsm.GT.0) THEN
     ht=ht/nsm
     gt=gt/nsm
  ELSE
     ht=-1.D9
     gt=0.d0
  END IF
  gma1=gmsun*(1+masst)
! Estimate of the mean motion (hence semimajor axis)                    
! allowing to combine two arcs   
  CALL start(elemp(1),elemp(2),1,ng,ennet,elemt(1),elemt(6),ok)
  IF(ok) THEN
     WRITE(*,*)'OFIDEN: initial conditions available, now differential corrections'
  ELSE
     WRITE(*,*)'OFIDEN : problem in start (initial guess) '
  ENDIF
  CALL titast(3,namof(1),namof(2),titnam,namet,lnt)
  WRITE(unirep,200) namet(1:lnt)
  CALL outele(unirep,elemt,'EQU',telemt,' ',.true.,.false.)
  IF(magopt.NE.0 .AND. ht.GT.-100.D0) WRITE(unirep,140) ht
  WRITE(unirep,206) namet(1:lnt)
206 FORMAT('Differential correction for object ',A,':')
  dir(3)=' '
  
  DO i=1,6
     icor(i)=1
  END DO
  IF(magopt.NE.0 .AND. ht.GT.-100.D0) THEN
     gmagc=gt
  END IF
  elem1=undefined_orbit_elem
  elem1%t=telemt
  elem1%coo='EQU'
  elem1%coord=elemt
  elem1%h_mag=ht
  elem1%g_mag=gmagc
  elem1%mag_set=(elem1%h_mag > -100.D0)
  CALL diff_cor(nt,obs,obsw,elem1,icor,unidif,elem(3),elem_unc(3),csinor,delnor,ok,nd)
  IF(ok) THEN
     deforb(3)=.true.
     defcn(3)=.true.
     comele(3)='computed with least squares fit'
! Magnitude estimation
     IF(magopt.NE.0) THEN
        CALL mag_est(nt,obs,obsw,hnew,rmsh)
        IF(rmsh.GT.0.d0) ht=hnew
     END IF
     elem(3)%h_mag=ht
     elem(3)%g_mag=gt
     elem(3)%mag_set=(elem(3)%h_mag > -100.D0)
     name(3)=namet
     namof(3)=namet
     mass(3)=masst
     nobj=3
     IF(opele) THEN
        OPEN(uniele,FILE=eleout,STATUS='UNKNOWN')
        opele=.false.
        WRITE(uniele,301) comcha
        CALL wromlh(uniele,'ECLM','J2000')
     END IF
     CALL coocha(elem(3)%coord,elem(3)%coo,gma1,elemo,oetype,enne)
     telemo=elem(3)%t
     WRITE(unirep,124)
     CALL outele(unirep,elemo,oetype,telemo,' ',.true.,.false.)
     IF(magopt.NE.0 .AND. elem(3)%h_mag.GT.-100.D0) WRITE(unirep,140) elem(3)%h_mag
     WRITE(unirep,230) csinor,delnor
     CALL write_elems(elem(3),name(3),UNC=elem_unc(3),UNIT=uniele)
     file=namet(1:lnt)//'.rwo'
     
     CALL write_rwo(file,obs,obsw,nt,error_model,csinor,rmsh)
  ELSE
     WRITE(unirep,202)
  END IF
202 FORMAT(5X,'FAILED')
301 FORMAT(A,' Orbits computed with differential corrections')
124 FORMAT(5X,'Corrected orbital elements:')
125 FORMAT(5X,'Intermediate orbital elements (a-M correction):')
230 FORMAT(5X,'Residual norm   =',1P,E11.3/                           &
         &       5X,'Correction norm =',1P,E11.3)
  
END SUBROUTINE ofiden

!
!  *****************************************************************
!  *                                                               *
!  *                         O F P R O P                           *
!  *                                                               *
!  *                Auxiliary routine for ORBFIT:                  *
!  *               propagation of orbital elements                 *
!  *                                                               *
!  *****************************************************************
!
! INPUT:    UNIREP    -  FORTRAN unit for report
!           UNIELE    -  FORTRAN unit for output of orbital elements
!           OPELE     -  Flag (need to open UNIELE)
!           ELEOUT    -  Name of file for output of orbital elements
!           NAME      -  Object names
!           DEFORB    -  Orbit definition flag
!           DEFCN     -  Tells whether covariance/normal matrices
!                            are defined
!           ELEM      -  Orbital elements
!           ELEM_UNC  -  Orbital element uncertainty
!           MASS      -  Object masses
!           COMELE    -  Comment on orbital elements
!           NOBJ      -  Number of objects
!           GMSUN     -  G*Mass(Sun)
!           OEPTIM    -  Epoch of output elements (MJD, TDT)
!           OEPSET    -  Flag stating that an output epoch is requested
!           OETYPE    -  Type of output elements (CAR/EQU/KEP/EQP)
!           ND        -  Dimension of the parameter space
!
SUBROUTINE ofprop(unirep,uniele,opele,eleout,name,deforb,defcn,elem,elem_unc,mass,comele,nobj,gmsun,oeptim,oepset,oetype,nd)
  USE propag_state
  
  INTEGER                          :: unirep,uniele,nobj
  TYPE(orbit_elem),DIMENSION(nobj) :: elem
  TYPE(orb_uncert),DIMENSION(nobj) :: elem_unc
  DOUBLE PRECISION                 :: gmsun
  DOUBLE PRECISION                 :: mass(nobj),oeptim
  LOGICAL                          :: opele,deforb(nobj),defcn(nobj),oepset
  CHARACTER*(*)                    :: eleout,name(nobj),comele(nobj),oetype
  INTEGER, INTENT(IN)              :: nd
!********************************************************************************

  INTEGER          :: i,ln,lc,fail
  DOUBLE PRECISION :: gma1,enne,dxde(6,6)
  DOUBLE PRECISION :: elemt(6),covet(6,6),noret(6,6)
  DOUBLE PRECISION :: elemo(6),coveo(6,6),noreo(6,6)
  LOGICAL          :: error,defnro
  TYPE(orbit_elem) :: elem1,elem2
  TYPE(orb_uncert) :: unc1,unc2
  
  INTEGER lench
  EXTERNAL lench
  
  IF(.NOT.oepset) RETURN
  WRITE(unirep,100)
100 FORMAT('Propagation of orbital elements:')
  
  DO 1 i=1,nobj
     IF(.NOT.deforb(i)) GOTO 1
     gma1=gmsun*(1+mass(i))
     
     IF(defcn(i)) THEN
        CALL pro_ele(elem(i),oeptim,elem1,elem_unc(i),unc1)
        CALL coo_cha(elem1,oetype,elem2,fail,dxde)
        CALL convertunc(unc1,dxde,unc2)
     ELSE
        CALL pro_ele(elem(i),oeptim,elem1)
        CALL coo_cha(elem1,oetype,elem2,fail)
     END IF
     ln=lench(name(i))
     lc=lench(comele(i))
     IF(lc.LE.0) THEN
        WRITE(unirep,101) name(i)(1:ln)
     ELSE
        WRITE(unirep,102) name(i)(1:ln),comele(i)(1:lc)
     END IF
     IF(dyn%nmod.GT.0) THEN
        CALL outele_ngr(unirep,elem2%coord,elem2%coo,elem2%t,' ',.true.,.false.,dyn%nmod,dyn%ndp,dyn%dp,ngr_opt)
     ELSE
        CALL outele(unirep,elem2%coord,elem2%coo,elem2%t,' ',.true.,.false.)
     END IF
     IF(opele) THEN
        OPEN(uniele,FILE=eleout,STATUS='UNKNOWN')
        opele=.false.
        WRITE(uniele,301) comcha
        CALL wromlh(uniele,'ECLM','J2000')
     END IF
     IF(lc.LE.0) THEN
        WRITE(uniele,103) comcha,name(i)(1:ln)
     ELSE
        WRITE(uniele,104) comcha,name(i)(1:ln),comele(i)(1:lc)
     END IF
     IF(dyn%nmod.GT.0) THEN
        CALL write_elems(elem2,name(i),DYNAM_PAR=dyn,NLSLOC=nls,LSLOC=ls,UNC=unc2,UNIT=uniele)
     ELSE
        CALL write_elems(elem2,name(i),UNC=unc2,UNIT=uniele)
     END IF
1 END DO
101 FORMAT(5X,'Propagated orbital elements for object ',A,':')
102 FORMAT(5X,'Propagated orbital elements for object ',A,' (',A,'):')
103 FORMAT(A,' Orbital elements of object ',A)
104 FORMAT(A,' Orbital elements of object ',A,' (',A,')')
301 FORMAT(A,' Propagated orbits')
  
END SUBROUTINE ofprop

!
!  *****************************************************************
!  *                                                               *
!  *                         O F E P H E                           *
!  *                                                               *
!  *                Auxiliary routine for ORBFIT:                  *
!  *                          ephemerides                          *
!  *                                                               *
!  *****************************************************************
!
! INPUT:    UNIT      -  FORTRAN unit for output
!           NAME      -  Object names
!           DEFORB    -  Orbit definition flag
!           DEFCOV    -  Orbit covariance definition flag
!           ELEM      -  Orbital elements
!           ELEM_UNC  -  Orbital element uncertainty
!           MASS      -  Object masses
!           COMELE    -  Comment on orbital elements
!           NOBJ      -  Number of objects
!
 SUBROUTINE ofephe(unit,name,deforb,defcov,elem,elem_unc,mass,comele,nobj)
   USE ephem_prop

   INTEGER                          :: unit,nobj
   TYPE(orbit_elem),DIMENSION(nobj) :: elem
   TYPE(orb_uncert),DIMENSION(nobj) :: elem_unc
   DOUBLE PRECISION                 :: mass(nobj)
   LOGICAL                          :: deforb(nobj),defcov(nobj)
   CHARACTER(LEN=*)                 :: name(nobj),comele(nobj)
! NEEDED common blocks:

   
   INTEGER i,k,ln,lc
   CHARACTER*100 fields ! ephemerides output 
   CHARACTER*3 scale
   
   INTEGER lench
   EXTERNAL lench
   
   IF(iiceph.NE.36) STOP '**** ofephe: internal error (01) ****'
   
   DO 1 k=1,nepobj
      i=kepobj(k)
      IF(i.LT.1 .OR. i.GT.nobj) GOTO 1
      IF(.NOT.deforb(i)) GOTO 1
      ln=lench(name(i))
      WRITE(unit,100) name(i)(1:ln)
100   FORMAT(/'Ephemeris for object ',A)
      lc=lench(comele(i))
      IF(lc.GT.0) WRITE(unit,101) comele(i)(1:lc)
101   FORMAT(5X,'Origin of orbital elements: ',A)
      WRITE(unit,102)
102   FORMAT(/)
      scale='UTC' 
! old version: fields is set hard-coded
!      fields='cal,mjd,coord,mag,elev,airm,elsun,elong,mooel,glat,glon,r,delta,appmot,skyerr' 
!      CALL ephemc(unit,elem(i),elem_unc(i),defcov(i),teph1,teph2,dteph,idsta,scale,fields)
! new version: fields is read from option files .oop
      CALL ephemc(unit,elem(i),elem_unc(i),defcov(i),teph1,teph2,dteph,idsta,scale,ephfld)
      
1  END DO
   
 END SUBROUTINE ofephe

END MODULE OFOD

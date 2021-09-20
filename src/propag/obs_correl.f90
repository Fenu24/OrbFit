MODULE obs_correl
IMPLICIT NONE

PUBLIC :: obscor

! Models of autocorrelation functions for all the observatories
!
! Max number of parameters/functions (all observatories together)
INTEGER, PARAMETER :: nparx=500
! Min and max value for observatory code
INTEGER,PARAMETER :: obsc1=0,obsc2=2000
! Models of autocorrelation functions for all the observatories
! pto2f       -  pointer from observatory code to function
! pto2fm      -  pointer from observatory code to function for mixed class
! nfo         -  number of functions for the model of each station
! nfom        -  number of functions for the model of mixed class
! nfunt       -  total number of functions
! npart       -  total number of parameters
! nparf       -  number of parameters for each function
! nparo       -  number of parameters for each observatory
! nparom      -  number of parameters for mixed class
! kfun        -  function integer code
! ptf2p       -  pointer from function to parameters
! par         -  parameter values
! aprmx       -  max a-priori RMS (arcsec)
! minapw      -  min a-priori weight (rad**(-2))
! maxdst      -  max "decision step" (1=specific obs; 2=mixed class)
! iiccor      -  initialization check
!
INTEGER nfunt,npart,pto2f(obsc1:obsc2,2,2),nfo(obsc1:obsc2,2)
INTEGER nparf(nparx),kfun(nparx),ptf2p(2,nparx)
INTEGER nparo(obsc1:obsc2,2),maxdst,iiccor
INTEGER pto2fm(2,2),nparom(2),nfom(2)
DOUBLE PRECISION par(nparx),aprmx,minapw

CONTAINS

! Copyright (C) 2000 by Mario Carpino (carpino@brera.mi.astro.it)
! Version: November 28, 2000
! ---------------------------------------------------------------------
!
!  *****************************************************************
!  *                                                               *
!  *                         O B S C O R                           *
!  *                                                               *
!  * Correlation between two observations of the same observatory  *
!  *                                                               *
!  *****************************************************************
!
! INPUT:    ERRMOD    -  Error model
!           OBS1      -  Observation 1
!           W1        -  A-priori weight for RA/DEC (rad**(-2)) (obs 1)
!           OBS2      -  Observation 2
!           W2        -  A-priori weight for RA/DEC (rad**(-2)) (obs 2)
!
! OUTPUT:   COVRA     -  Off-diagonal covariance in RA (rad**2)
!           COVDEC    -  Off-diagonal covariance in DEC (rad**2)
!
SUBROUTINE obscor(errmod,obs1,w1,obs2,w2,covra,covdec)
USE astrometric_observations
IMPLICIT NONE

CHARACTER(LEN=*),               INTENT(IN)  :: errmod
TYPE(ast_obs),                  INTENT(IN)  :: obs1,obs2
DOUBLE PRECISION, DIMENSION(2), INTENT(IN)  :: w1,w2
DOUBLE PRECISION,               INTENT(OUT) :: covra,covdec

DOUBLE PRECISION :: rmsc1,rmsc2,bias1,bias2,dt
INTEGER, DIMENSION(2) :: idst1,idst2
LOGICAL,SAVE :: first=.true.,doit
CHARACTER*80 file
INTEGER :: lf,ic,ip1,if1
INTEGER, EXTERNAL :: lench
CHARACTER*20,SAVE :: errmodkeep
LOGICAL :: mixcl

! Initialization
IF(first) THEN
    first=.false.
    doit=(errmod.ne.' '.AND.errmod.ne.'cbm10'.AND.errmod.ne.'fcct14'.AND.errmod.ne.'vfcc17')
    errmodkeep=errmod
! Input of correlation model
    IF(doit) THEN
        file=errmod//'.cor'
        CALL rmsp(file,lf)
        WRITE(*,210) file(1:lf)
210 FORMAT('INFO(obscor): reading error correlation model from file "',A,'"')
        CALL rdcorm(file)
    ELSE
        WRITE(*,211)
    END IF
END IF

IF(errmod.ne.errmodkeep)THEN
   WRITE(*,*)' obscor: change in errmod from ',errmodkeep,' to ', errmod,' not allowed'
   STOP
ENDIF
211 FORMAT('INFO(obscor): NO error correlation model is used')

covra=0
covdec=0
IF(.NOT.doit) RETURN
IF(iiccor.NE.36) STOP '**** obscor: internal error (01) ****'
IF(obs1%obscod_i /= obs2%obscod_i) RETURN
IF(obs1%type /= obs2%type) RETURN
IF(obs1%type /= 'O') RETURN
IF(obs1%obscod_i.LT.obsc1.OR.obs1%obscod_i.GT.obsc2) STOP '**** obscor: internal error (02) ****'
dt=ABS(obs2%time_tdt-obs1%time_tdt)

! These calls are needed only for computing idst1 and idst2
CALL astrow(errmod,obs1%tech,obs1%time_tdt,obs1%obscod_i,obs1%acc_coord(1),obs1%acc_coord(2),  &
            rmsc1,rmsc2,bias1,bias2,idst1)
CALL astrow(errmod,obs2%tech,obs2%time_tdt,obs2%obscod_i,obs2%acc_coord(1),obs2%acc_coord(2),  &
            rmsc1,rmsc2,bias1,bias2,idst2)

! Right ascension
ic=1
IF(idst1(ic).EQ.0 .OR. idst2(ic).EQ.0) GOTO 1
IF(idst1(ic).GT.maxdst) GOTO 1
IF(idst2(ic).GT.maxdst) GOTO 1
IF(w1(1).LT.minapw) GOTO 1
IF(w2(1).LT.minapw) GOTO 1
if1=pto2f(obs1%obscod_i,ic,1)
mixcl=(if1.LE.0)
IF(mixcl) THEN
   if1=pto2fm(ic,1)
   IF(if1.LE.0) GOTO 1
   ip1=ptf2p(1,if1)
   CALL fcorob(kfun(if1),nfom(ic),par(ip1),nparf(if1),dt,covra)
ELSE
   ip1=ptf2p(1,if1)
   CALL fcorob(kfun(if1),nfo(obs1%obscod_i,ic),par(ip1),nparf(if1),dt,covra)
END IF
covra=covra/SQRT(w1(1)*w2(1))
1 CONTINUE

! Declination
ic=2
IF(idst1(ic).EQ.0 .OR. idst2(ic).EQ.0) GOTO 2
IF(idst1(ic).GT.maxdst) GOTO 2
IF(idst2(ic).GT.maxdst) GOTO 2
IF(w1(2).LT.minapw) GOTO 2
IF(w2(2).LT.minapw) GOTO 2
if1=pto2f(obs1%obscod_i,ic,1)
mixcl=(if1.LE.0)
IF(mixcl) THEN
   if1=pto2fm(ic,1)
   IF(if1.LE.0) GOTO 2
   ip1=ptf2p(1,if1)
   CALL fcorob(kfun(if1),nfom(ic),par(ip1),nparf(if1),dt,covdec)
ELSE
   ip1=ptf2p(1,if1)
   CALL fcorob(kfun(if1),nfo(obs1%obscod_i,ic),par(ip1),nparf(if1),dt,covdec)
END IF
covdec=covdec/SQRT(w1(2)*w2(2))
2 CONTINUE

END SUBROUTINE obscor

! Copyright (C) 2000 by Mario Carpino (carpino@brera.mi.astro.it)
! Version: November 16, 2000
! ---------------------------------------------------------------------
!
!  *****************************************************************
!  *                                                               *
!  *                         F C O R O B                           *
!  *                                                               *
!  *              Computation of correlation function              *
!  *                                                               *
!  *****************************************************************
!
! INPUT:    KFUN      -  Function integer codes
!           NFUN      -  Number of elementary functions
!           PAR       -  Parameter values
!           NPARF     -  Number of parameters per function
!           DT        -  Time lag
!
! OUTPUT:   FUN       -  Function value
!
      SUBROUTINE fcorob(kfun,nfun,par,nparf,dt,fun)
      IMPLICIT NONE

      INTEGER nfun
      INTEGER nparf(nfun),kfun(nfun)
      DOUBLE PRECISION par(*),dt,fun

      INTEGER if1,ip1,kf1
      DOUBLE PRECISION par1,ff1,fd1,term

      fun=0
      term=0

      ip1=1
      DO 1 if1=1,nfun
      kf1=kfun(if1)
      par1=par(ip1)

      INCLUDE 'fcfund.h90'

      ELSE
          STOP '**** fcorob: internal error (01) ****'
      END IF

      IF(kf1.EQ.1) THEN
          fun=fun+term
          term=ff1
      ELSE
          term=term*ff1
      END IF
      ip1=ip1+nparf(if1)

    1 CONTINUE
      fun=fun+term

      END SUBROUTINE fcorob
! Copyright (C) 2000 by Mario Carpino (carpino@brera.mi.astro.it)
! Version: November 28, 2000
! ---------------------------------------------------------------------
!
!  *****************************************************************
!  *                                                               *
!  *                         R D C O R M                           *
!  *                                                               *
!  *   Read the autocorrelation models for all the observatories   *
!  *                                                               *
!  *****************************************************************
!
! INPUT:    FILE      -  Input file
!
! OUTPUT:
!
      SUBROUTINE rdcorm(file)
      USE fund_const
      IMPLICIT NONE

      CHARACTER*(*) file

      INTEGER unit,io,line,ic,lf,ip,ivers,lr
      CHARACTER(LEN=80) :: rec,rec1
      LOGICAL :: mixcl

      INTEGER lench
      EXTERNAL lench

! Initialization
      DO 1 io=obsc1,obsc2
      pto2f(io,1,1)=0
      pto2f(io,2,1)=0
    1 END DO
      pto2fm(1,1)=0
      pto2fm(2,1)=0
      nfunt=0
      npart=0
      line=0

! Reading model
      lf=lench(file)
      CALL filopl(unit,file)
      READ(unit,*) ivers
      IF(ivers.NE.1) THEN
          WRITE(*,201) file(1:lf)
          STOP '**** Abnormal end ****'
      END IF
  201 FORMAT('ERROR(rdcorm): unsupported version of file "',A,'"')
      READ(unit,*) aprmx
      READ(unit,*) maxdst
      minapw=1.D0/((aprmx*radsec)**2)
      line=3

    2 CONTINUE
! Read a new observatory/coordinate
      READ(unit,100,END=10) rec
  100 FORMAT(A)
      line=line+1
! IC = coordinate (1=RA, 2=DEC)
      IF(rec(1:8).EQ.'TIME RA ') THEN
          ic=1
      ELSEIF(rec(1:8).EQ.'TIME DEC') THEN
          ic=2
      ELSE
          GOTO 20
      END IF
! IO = observatory code
      rec1=rec(9:)
      CALL norstr(rec1,lr)
      mixcl=(rec1.EQ.'MIX')
      IF(mixcl) THEN
          pto2fm(ic,1)=nfunt+1
          nfom(ic)=0
          nparom(ic)=0
      ELSE
          READ(rec1,*,ERR=20) io
          IF(io.LT.obsc1) STOP '**** rdcorm: obsc < obsc1 ****'
          IF(io.GT.obsc2) STOP '**** rdcorm: obsc > obsc2 ****'
          pto2f(io,ic,1)=nfunt+1
          nfo(io,ic)=0
          nparo(io,ic)=0
      END IF

    3 CONTINUE
! Read the model for one coordinate (IC) of one observatory (IO)
      READ(unit,100,ERR=20) rec
      line=line+1
      IF(rec.EQ.'END') THEN
          IF(mixcl) THEN
              IF(nfom(ic).LE.0) GOTO 20
          ELSE
              IF(nfo(io,ic).LE.0) GOTO 20
          END IF
          GOTO 2
      END IF
      nfunt=nfunt+1
      IF(nfunt.GT.nparx) STOP '**** rdcorm: nfunt > nparx ****'
      IF(mixcl) THEN
          pto2fm(ic,2)=nfunt
          nfom(ic)=nfom(ic)+1
      ELSE
          pto2f(io,ic,2)=nfunt
          nfo(io,ic)=nfo(io,ic)+1
      END IF
! Determine function integer code (KP1) and number of parameters (NP1)
      CALL fcsfun(rec,kfun(nfunt),nparf(nfunt))
      IF(kfun(nfunt).LE.0) GOTO 20
      IF(mixcl) THEN
          nparom(ic)=nparom(ic)+nparf(nfunt)
      ELSE
          nparo(io,ic)=nparo(io,ic)+nparf(nfunt)
      END IF
! Reading and storing parameter values and properties
      IF(npart+nparf(nfunt).GT.nparx) STOP '**** rdcorm: npart > nparx ****'
      ptf2p(1,nfunt)=npart+1
      DO 4 ip=1,nparf(nfunt)
      npart=npart+1
      READ(unit,*,ERR=20) par(npart)
    4 END DO
      ptf2p(2,nfunt)=npart
      GOTO 3

! Regular end
   10 CONTINUE
      CALL filclo(unit,' ')
      iiccor=36
      RETURN

! Error termination
   20 CONTINUE
      WRITE(*,200) file(1:lf),line
  200 FORMAT('rdcorm: INPUT ERROR from file "',A,'" at record',I5)
      STOP '**** rdcorm: abnormal end ****'

      END SUBROUTINE rdcorm
! Copyright (C) 2000 by Mario Carpino (carpino@brera.mi.astro.it)
! Version: November 3, 2000
! ---------------------------------------------------------------------
!
!  *****************************************************************
!  *                                                               *
!  *                         F C S F U N                           *
!  *                                                               *
!  *              LS fit of covariance functions:                  *
!  *       function integer codes and number of parameters         *
!  *                                                               *
!  *****************************************************************
!
! INPUT:    NAME      -  Function name
!
! OUTPUT:   KFUN      -  Function integer identificator (0=don't know)
!           NPAR      -  Number of parameters
!
      SUBROUTINE fcsfun(name,kfun,npar)
      IMPLICIT NONE

      INTEGER kfun,npar
      CHARACTER*(*) name

      kfun=0
      npar=0
! Multiplicative coefficient (constant)
      IF(name.EQ.'+COEF') THEN
          kfun=1
          npar=1
! Exponential function
      ELSEIF(name.EQ.'*EXP') THEN
          kfun=2
          npar=1
! Normal function
      ELSEIF(name.EQ.'*NORM') THEN
          kfun=3
          npar=1
! Parabola (to be multiplied by exponential or normal functions)
      ELSEIF(name.EQ.'*PARAB') THEN
          kfun=4
          npar=1
! ADD HERE NEW FUNCTIONS, using increasing integer identificator (KFUN)
! Remember to update accordingly also subroutines FCFUND and FCWPAR
      END IF

      END SUBROUTINE fcsfun
! Copyright (C) 1997-1999 by Mario Carpino (carpino@brera.mi.astro.it)
! Version: November 26, 1999
! ---------------------------------------------------------------------
!
!  *****************************************************************
!  *                                                               *
!  *                         A S T R O W                           *
!  *                                                               *
!  *     Computation of a-priori RMS of optical observations       *
!  *         from results of a statistical analysis or             *
!  *    from their accuracy (= number of digits in input file)     *
!  *                                                               *
!  *****************************************************************
!
! INPUT:    MPCTYP    -  Observation type (column 15 of MPC record)
!           TDT       -  Time of observation (MJD, TDT)
!           IDSTA     -  Observatory code 9numeric)
!           ACCA      -  Accuracy of right ascension (rad)
!           ACCD      -  Accuracy of declination (rad)
!
! OUTPUT:   RMSA      -  A-priori RMS of right ascension (rad)
!           RMSD      -  A-priori RMS of declination (rad)
!           BIASA     -  Bias in  right ascension (rad)
!           BIASD     -  Bias in  declination (rad)
!           DECSTEP   -  Decision steps of RMS assignment
!
      SUBROUTINE astrow(error_model,mpctyp,tdt,idsta,acca,accd,        &
     &        rmsa,rmsd,biasa,biasd,decstep)
      USE fund_const
      IMPLICIT NONE
      CHARACTER*(*),INTENT(IN)    :: error_model
      INTEGER, INTENT(IN):: idsta
      DOUBLE PRECISION, INTENT(IN) :: tdt,acca,accd
      CHARACTER*(*), INTENT(IN) :: mpctyp

      DOUBLE PRECISION, INTENT(OUT) :: rmsa,rmsd,biasa,biasd
      INTEGER, DIMENSION(2), INTENT(OUT) :: decstep

      INTEGER le
      DOUBLE PRECISION rmsmin
      CHARACTER*5 ads(2)
      CHARACTER*80 filea,filed,ermnam
      LOGICAL ermuse,error

! Additional info on how RMS are obtained
! ---------------------------------------------------------------------
! former ASTROW.H
! RMS classes of astrometric observations: additional information
!
! idcl        -  Class progressive number in index
! orstep      -  Decision step (RA/DEC) :
!                    0 = no information
!                    1 = single-station class
!                    2 = multi-station class
!                    3 = default
! tdtlim      -  Class limits (MJD, TDT)
!
      INTEGER idcl(2),orstep(2)
!      COMMON/cmaow1/idcl,orstep
      DOUBLE PRECISION tdtlim(2,2)
!      COMMON/cmaow2/tdtlim

      INTEGER lench
      EXTERNAL lench
      LOGICAL first
      DATA first/.true./
      SAVE first
! decide on use of error model file
      ermuse=(error_model.ne.' '.AND.error_model.ne.'cbm10')
! Input of RMS class definitions
      IF(first.and.ermuse) THEN
         le=lench(error_model)
         filea=error_model(1:le)//'.cla'
         filed=error_model(1:le)//'.cld'
         CALL rrmscl(filea,filed,.false.)
         first=.false.
      ENDIF
! Negative values means: not yet assigned
      rmsa=-1.D0
      rmsd=-1.D0
      biasa=0.d0
      biasd=0.d0
      orstep(1)=0
      orstep(2)=0
! STEPs 1/2: look for a specific RMS class for the observatory
      IF(ermuse) THEN
          CALL accstr(acca,accd,ads(1),ads(2),error)
          IF(error) GOTO 2
          CALL crmscl(idsta,ads,mpctyp,tdt,rmsa,rmsd,biasa,biasd,idcl,orstep,tdtlim)
          IF(orstep(1).GT.0)THEN
               rmsa=rmsa*radsec
               biasa=biasa*radsec
          ENDIF
          IF(orstep(2).GT.0)THEN
               rmsd=rmsd*radsec
               biasd=biasd*radsec
          ENDIF
      ENDIF

! STEP 3: default, rough rule-of-thumb error model
    2 CONTINUE
      IF(rmsa.LT.0.D0 .OR. rmsd.LT.0.D0) THEN
! Deafault value of RMS (time dependent)
! Before 1890
          IF(tdt.LT.11368.d0) THEN
              rmsmin=3.d0
! From 1890 to 1950
          ELSE IF(tdt.LT.33282.d0) THEN
              rmsmin=2.d0
! After 1950
          ELSE
              rmsmin=1.d0
          END IF
          rmsmin=rmsmin*radsec
          IF(rmsa.LT.0.D0) THEN
              rmsa=MAX(rmsmin,acca)
              orstep(1)=3
          END IF
          IF(rmsd.LT.0.D0) THEN
              rmsd=MAX(rmsmin,accd)
              orstep(2)=3
          END IF
      END IF
      decstep=orstep

      END SUBROUTINE astrow
! Copyright (C) 1999-2000 by Mario Carpino (carpino@brera.mi.astro.it)
! Version: November 28, 2000
! ---------------------------------------------------------------------
!
!  *****************************************************************
!  *                                                               *
!  *                         R R M S C L                           *
!  *                                                               *
!  *            Read from a file and store RMS classes             *
!  *                                                               *
!  *****************************************************************
!
! INPUT:    FILRA     -  File name (RA classes)
!           FILDEC    -  File name (DEC classes)
!           CLONLY    -  Load only class list, not accuracy description
!                        (this option is enabled when the routine is
!                        used by stat2.f)
!
! WARNING: the correct behaviour of the routine depends strictly on the
! sorting of input data. More precisely:
!   - "per-observatory" RMS classes must be sorted in increasing order
!      of observatory code (field 3);
!   -  "mixed" (all observatories together) classes must be sorted
!      according to the following sorting keys:
!         key 1: accuracy class (field 1)
!         key 2: type of observation (field 2)
!         key 3: starting time (field 4)
!
   SUBROUTINE rrmscl(filra,fildec,clonly)
      USE station_coordinates
      IMPLICIT NONE

      CHARACTER*(*) filra,fildec
      LOGICAL clonly

      INCLUDE 'parobc.h90'
      INCLUDE 'parrms.h90'

! Common blocks to be initialized:
      INCLUDE 'rmscl.h90'

      INTEGER ic,i,unit,n,obsc,obscp,nin,lf
      INTEGER n1,n2,n3
      DOUBLE PRECISION mjde1,mjde2,rms1,ave1,rms1a
      CHARACTER crmad1*5,crmod1*1,cobsc*3,file*100,rec*120
      LOGICAL mixobs,new1,new2

      INTEGER lench
      DOUBLE PRECISION tjm1
      EXTERNAL lench,tjm1

      DO 10 ic=1,2
      IF(ic.EQ.1) THEN
          file=filra
      ELSE
          file=fildec
      END IF
      lf=lench(file)
      CALL filopl(unit,file)

! Initializations (section 1)
      ncrm(ic)=0
      DO 1 i=obsc1,obsc2
      crmobp(1,i,ic)=0
      crmobp(2,i,ic)=0
    1 END DO
      nin=0
      n=1
      obscp=obsc1-999

! Initializations (section 2)
      n1=0
      n2=0
      n3=0

! Check file format (in order to avoid using files written according
! to old format, not containing bias information)
      IF(.NOT.clonly) THEN
          READ(unit,102,END=3) rec
          IF(lench(rec).LT.110) THEN
              WRITE(*,200) file(1:lf)
              STOP '**** rrmscl: abnormal end ****'
          END IF
          REWIND(unit)
      END IF
  102 FORMAT(A)
  200 FORMAT('**** ERROR: file "',A,'" is written in an obsolete ',     &
     &       'format: please use a more recente version ****')

! Start reading loop
    2 CONTINUE
      IF(clonly) THEN
          READ(unit,100,END=3) crmad1,crmod1,cobsc,mjde1,mjde2
          rms1=0
          ave1=0
          rms1a=0
      ELSE
          READ(unit,100,END=3) crmad1,crmod1,cobsc,mjde1,mjde2,         &
     &                         rms1,ave1,rms1a
      END IF
  100 FORMAT(A5,1X,A1,1X,A3,1X,F8.1,1X,F8.1,1X,F9.3,F11.5,F9.3,         &
     &       E11.3,F9.3,2I8,F7.2,I9)
      nin=nin+1
      IF(crmod1.EQ.' ') crmod1='P'

      mixobs=(cobsc.EQ.'ALL')

      IF(mixobs) THEN
          n3=n3+1
          IF(n3.GT.crx3nx) STOP '**** rrmscl: n3 > crx3nx ****'

! Check for change in accuracy descriptor
          IF(n1.LE.0) THEN
              new1=.true.
          ELSE
              new1=(crmad1.NE.crx1ad(n1,ic))
          END IF
          IF(new1) THEN
              n1=n1+1
              IF(n1.GT.crx1nx) STOP '**** rrmscl: n1 > crx1nx ****'
              n2=n2+1
              IF(n2.GT.crx2nx) STOP '**** rrmscl: n2 > crx2nx ****'
              crx1ad(n1,ic)=crmad1
              crx1pt(1,n1,ic)=n2
              crx1pt(2,n1,ic)=n2
              crx2ty(n2,ic)=crmod1
              crx2pt(1,n2,ic)=n3
          ELSE

! Check for change in observation type
              IF(n2.LE.0) THEN
                  new2=.true.
              ELSE
                  new2=(crmod1.NE.crx2ty(n2,ic))
              END IF
              IF(new2) THEN
                  n2=n2+1
                  IF(n2.GT.crx2nx) STOP '**** rrmscl: n2 > crx2nx ****'
                  crx1pt(2,n1,ic)=n2
                  crx2ty(n2,ic)=crmod1
                  crx2pt(1,n2,ic)=n3
              END IF
          END IF
          crx2pt(2,n2,ic)=n3
          crx3t(1,n3,ic)=mjde1
          crx3t(2,n3,ic)=mjde2
          crx3r(n3,ic)=rms1
          crx3a(n3,ic)=ave1
          crx3ra(n3,ic)=rms1a
      ELSE
          IF(n.GT.ncrmx) STOP '**** rrmscl: ncrm > ncrmx ****'
          CALL statcode(cobsc,obsc)
!         READ(cobsc,101) obsc
          IF(obsc.LT.obsc1 .OR. obsc.GT.obsc2)                          &
     &        STOP '**** rrmscl: input error (01) ****'
          crmad(n,ic)=crmad1
          crmotd(n,ic)=crmod1

! Integer codification of type descriptor
          crmoti(n,ic)=1000+ICHAR(crmotd(n,ic))
          crmot1(n,ic)=mjde1
          crmot2(n,ic)=mjde2
          crmrms(n,ic)=rms1
          crmave(n,ic)=ave1
          crmrma(n,ic)=rms1a

          IF(obsc.NE.obscp) THEN
              IF(crmobp(1,obsc,ic).NE.0)                                &
     &            STOP '**** rrmscl: input error (02) ****'
              crmobp(1,obsc,ic)=n
              obscp=obsc
          END IF
          crmobp(2,obsc,ic)=n
          n=n+1
      END IF
  101 FORMAT(I3)

      GOTO 2

    3 CONTINUE
      ncrm(ic)=n-1
      crx1n(ic)=n1
      crx2n(ic)=n2
      crx3n(ic)=n3
      CALL filclo(unit,' ')
   10 END DO

      iiccrm=36

    END SUBROUTINE rrmscl

! Copyright (C) 1999 by Mario Carpino (carpino@brera.mi.astro.it)
! Version: March 3, 1999
! ---------------------------------------------------------------------
!
!  *****************************************************************
!  *                                                               *
!  *                         A C C S T R                           *
!  *                                                               *
!  *               Accuracy description (string)                   *
!  *                                                               *
!  *****************************************************************
!
! INPUT:    ACCA      -  Accuracy of right ascension (rad)
!           ACCD      -  Accuracy of declination (rad)
!
! OUTPUT:   ADSA      -  Accuracy description string (RA)
!           ADSD      -  Accuracy description string (DEC)
!           ERROR     -  Error flag
!
   SUBROUTINE accstr(acca,accd,adsa,adsd,error)
      USE fund_const
      IMPLICIT NONE

      DOUBLE PRECISION acca,accd
      CHARACTER*(*) adsa,adsd
      LOGICAL error

      DOUBLE PRECISION accs
      INTEGER i,k,lads
      CHARACTER*10 ads

      error=.false.

      DO 1 i=1,2
      IF(i.EQ.1) THEN
          accs=acca*secrad/15.d0
      ELSE
          accs=accd*secrad
      END IF
      WRITE(ads,101) accs
  101 FORMAT(F10.4)
      CALL rmsp(ads,lads)
      DO 2 k=lads,1,-1
      IF(ads(k:k).EQ.'0') THEN
          ads(k:k)=' '
          lads=lads-1
      ELSE
          GOTO 3
      END IF
    2 END DO
    3 CONTINUE
      IF(ads(k:k).EQ.'.') THEN
          ads(k:k)=' '
          lads=lads-1
      END IF
      IF(i.EQ.1) THEN
          IF(lads.GT.LEN(adsa)) THEN
              error=.true.
              WRITE(*,200) 'ADSA',ads(1:lads)
          ELSE
              adsa=ads(1:lads)
          END IF
      ELSE
          IF(lads.GT.LEN(adsd)) THEN
              error=.true.
              WRITE(*,200) 'ADSD',ads(1:lads)
          ELSE
              adsd=ads(1:lads)
          END IF
      END IF
    1 END DO
  200 FORMAT('ERROR (accstr): ',A,' = "',A,'"')

 END SUBROUTINE accstr
! Copyright (C) 1999 by Mario Carpino (carpino@brera.mi.astro.it)
! Version: November 24, 1999
! ---------------------------------------------------------------------
!
!  *****************************************************************
!  *                                                               *
!  *                         C R M S C L                           *
!  *                                                               *
!  *   Computes a-priori observation RMS based on known classes    *
!  *                                                               *
!  *****************************************************************
!
! INPUT:    OBSCOD    -  Observatory code
!           ADS       -  Accuracy description string (RA/DEC)
!           MPCTYP    -  Observation type (column 15 of MPC record)
!           TDT       -  Time (MJD, TDT)
!
! OUTPUT:   RMSA      -  A-priori RMS of RA (arcsec)
!           RMSD      -  A-priori RMS of DEC (arcsec)
!           BIASA     -  Bias in RA (arcsec)
!           BIASD     -  Bias in DEC (arcsec)
!           IDCL      -  Class progressive number in index
!           STEP      -  RMS assignation steps:
!                           0 = not assigned
!                           1 = specific station class
!                           2 = mixed class
!           TDTLIM    -  Class limits (MJD, TDT)
!
! WARNING: if no valid class for the observation is found, the
!          corresponding output rms is set to a negative value
!
  SUBROUTINE crmscl(obscod,ads,mpctyp,tdt,rmsa,rmsd,biasa,biasd,idcl,step,tdtlim)
      IMPLICIT NONE

      INTEGER obscod,idcl(2),step(2)
      DOUBLE PRECISION rmsa,rmsd,biasa,biasd,tdt,tdtlim(2,2)
      CHARACTER*(*) mpctyp,ads(2)

      INCLUDE 'parobc.h90'
      INCLUDE 'parrms.h90'

! NEEDED common blocks:
      INCLUDE 'rmscl.h90'

      INTEGER ic,i,i1,i2,i3
      DOUBLE PRECISION rms1,bias1

      IF(iiccrm.NE.36) STOP '**** crmscl: internal error (01) ****'
      IF(obscod.LT.obsc1 .OR. obscod.GT.obsc2)                          &
     &    STOP '**** crmscl: input error (01) ****'

      DO 20 ic=1,2
      rms1=-1.d0
      bias1=0.d0
      step(ic)=0
      idcl(ic)=0
      tdtlim(1,ic)=-9.D9
      tdtlim(2,ic)=-9.D9
      IF(crmobp(1,obscod,ic).GT.0) THEN
          DO 1 i=crmobp(1,obscod,ic),crmobp(2,obscod,ic)
          IF(crmad(i,ic).NE.ads(ic)) GOTO 1
          IF(crmotd(i,ic).NE.mpctyp) GOTO 1
          IF(tdt.LT.crmot1(i,ic) .OR. tdt.GE.crmot2(i,ic)) GOTO 1
          rms1=crmrms(i,ic)
          bias1=crmave(i,ic)
          idcl(ic)=i
          step(ic)=1
          tdtlim(1,ic)=crmot1(i,ic)
          tdtlim(2,ic)=crmot2(i,ic)
          GOTO 2
    1     CONTINUE
    2     CONTINUE
      END IF
      IF(step(ic).GT.0) GOTO 10
      DO 3 i1=1,crx1n(ic)
      IF(ads(ic).EQ.crx1ad(i1,ic)) THEN
          DO 4 i2=crx1pt(1,i1,ic),crx1pt(2,i1,ic)
          IF(mpctyp.EQ.crx2ty(i2,ic)) THEN
              DO 5 i3=crx2pt(1,i2,ic),crx2pt(2,i2,ic)
              IF(tdt.GE.crx3t(1,i3,ic) .AND. tdt.LT.crx3t(2,i3,ic)) THEN
                  rms1=crx3r(i3,ic)
                  bias1=crx3a(i3,ic)
                  idcl(ic)=i3
                  step(ic)=2
                  tdtlim(1,ic)=crx3t(1,i3,ic)
                  tdtlim(2,ic)=crx3t(2,i3,ic)
                  GOTO 10
              END IF
    5         CONTINUE
          END IF
    4     CONTINUE
      END IF
    3 END DO

   10 CONTINUE
      IF(ic.EQ.1) THEN
          rmsa=rms1
          biasa=bias1
      ELSE
          rmsd=rms1
          biasd=bias1
      END IF
   20 END DO

 END SUBROUTINE crmscl


END MODULE obs_correl

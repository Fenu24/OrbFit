! Copyright (C) 2012 by OrbFit consortium       
! Version: January 4, 2012                                            
! --------------------------------------------------------------------- 
!                                                                       
! Generation of a binary (random access) ephemeris file                 
!                                                                       
PROGRAM bineph2 
  USE reference_systems 
  USE orbit_elements 
  USE fund_const 
  USE propag_state
  USE force_model
  USE astrometric_observations
  USE output_control
  USE least_squares
  IMPLICIT NONE 
                                                                        
  INCLUDE 'dim.h90' 
  INCLUDE 'parobx.h90' ! observation numbers: maximum
  CHARACTER*80 run,optfil,auxfil,binfil 
  INTEGER nbout,lr,nout,ln,i,k,recl,nhl,ll,nv,it
! no records, forward,backward
  INTEGER nrec,iptf1,iptb1
  INTEGER unit
  INTEGER :: fail_flag ! for coo_cha
  TYPE(orbit_elem), DIMENSION(nbmx) :: el0,el0aux ! initial conditions 
  DOUBLE PRECISION, DIMENSION(nbmx) :: astmas ! masses
  CHARACTER*30, DIMENSION(nbmx) :: astnam !names
  TYPE(orbit_elem), ALLOCATABLE ::  el1(:,:)
  DOUBLE PRECISION teph1,teph2,dteph 
  DOUBLE PRECISION t0,t 
  DOUBLE PRECISION tmp,xl,enne,elem(6,nbmx),gma,gma1
  LOGICAL found 
! for diffcor
  LOGICAL defobs, change, diffcor_bineph, ok, cov0, succ
  INTEGER le
  CHARACTER*200 rwofil
  CHARACTER*100 obsdir, file
  CHARACTER*20 error_model
  INTEGER m, usedobs ! number of observations
  TYPE(ast_obs),  DIMENSION(nobx) :: obs    ! observations
  TYPE(ast_wbsr), DIMENSION(nobx) :: obsw   ! weights and possibly residuals
  TYPE(orbit_elem) elc
  TYPE(orb_uncert) unc
  DOUBLE PRECISION csino,delno,rmsh
  INTEGER :: iundone ! done file
  !
  INTEGER lench              
  EXTERNAL lench 
                                                                        
  WRITE(*,110) 
110 FORMAT(' Run name =') 
  READ(*,100) run 
100 FORMAT(A) 
  lr=lench(run) 
  WRITE(*,111) run(1:lr) 
111 FORMAT(' Run name = ',A) 
                                                                        
! Input of general options    
  CALL initopt('bineph',run,'bop')
! options of propagator
  CALL rmodel(1)
! open log/err file
  file=run(1:lr)//'.log' 
  CALL filopn(iun_log,file,'UNKNOWN') 
  file=run(1:lr)//'.err' 
  CALL filopn(ierrou,file,'UNKNOWN')
  numerr=0
! open clo file
  file=run(1:lr)//'.clo' 
! close-approach file to be opened fresh                                
  CALL filopn(iuncla,file,'UNKNOWN') 
  numcla=0
! options of ephemerides
  CALL rdoptb2(nbout,el0aux,astmas,astnam,teph1,teph2,dteph,    &
&               diffcor_bineph,error_model,obsdir) 
  CALL errmod_set(error_model)
  DO k=1,nbout
     CALL coo_cha(el0aux(k),'KEP',el0(k),fail_flag)
     IF(fail_flag.ge.4)THEN
        WRITE(*,*)' bineph2: conversion to KEP for 2-body failed ',fail_flag,el0aux(k)
        STOP
     END IF
  END DO
  IF(nbout.LE.0) STOP '**** bineph: object list is empty ****' 
! Initializations of times
  dteph=ABS(dteph) 
  IF(teph1.GT.teph2) THEN 
     tmp=teph1 
     teph1=teph2 
     teph2=tmp 
  END IF
  nout=(teph2-teph1)/dteph+1 
  IF(teph2.GT.teph1+dteph*(nout-1)) nout=nout+1 
  teph2=teph1+dteph*(nout-1) 
! Deleting old copies of output files                                   
  auxfil=run(1:lr)//'.bai' 
  binfil=run(1:lr)//'.bep' 
  CALL dlifex(auxfil) 
  CALL dlifex(binfil) 
! Auxiliary file                                                        
  OPEN(1,FILE=auxfil,STATUS='UNKNOWN') 
  WRITE(1,200) nbout,teph1,teph2,dteph 
200 FORMAT(I5/3F13.5) 
  DO  i=1,nbout 
     ln=lench(astnam(i)) 
     WRITE(1,201) astmas(i),astnam(i)(1:ln) 
201  FORMAT(1P,E18.10,1X,A) 
  END DO
  CLOSE(1) 
! Binary ephemeris file                                                 
! Assumes RECL given in single precision words                          
  recl=8*6*nbout 
  OPEN(1,FILE=binfil,ACCESS='DIRECT',RECL=recl,STATUS='NEW') 
! First record: version number, number of bodies, number of records     
  WRITE(1,REC=1) 102,nbout,nout 
! Second record: first and last epoch (MJD, TDT), stepsize (d)          
  WRITE(1,REC=2) teph1,teph2,dteph 
! Records from 3 to 2+nbout: masses                                     
  DO i=1,nbout 
     gma=gms*astmas(i) 
     gma1=gma+gms 
     WRITE(1,REC=2+i) astmas(i),gma,gma1
  END DO
! Records from 3+nbout to 2+2*nbout: names                              
  DO i=1,nbout 
     WRITE(1,REC=2+nbout+i) astnam(i) 
  END DO
  nhl=2+2*nbout 
! Records from 3+nbout to end: keplerian elements (a,e,i,N,w,M);        
! units: AU, radians; refsys: ECLM J2000 
! allocation of temporary storage array
  nrec=CEILING((teph2-teph1)/dteph+1)
  ALLOCATE(el1(nbout,nrec))
! Loop on asteroids                                                 
  DO k=1,nbout
     t0=el0(k)%t
! First point to be computed in the forward integration 
     tmp=(t0-teph1)/dteph+1 
! records to be written backward from IC
     iptf1=tmp 
     IF(tmp-iptf1.GT.0.d0) iptf1=iptf1+1 
! backward integration to teph1
! avoiding self perturbations
     CALL selpert(astnam(k),found)
     IF(found)THEN 
        WRITE(*,*)astnam(k),' massive asteroid, no selfperturbation' 
     ENDIF
     CALL set_restart(.true.)
     IF(diffcor_bineph)THEN
! differential correction to make orbit consistent with propagation
!  model as used in the current run
        verb_dif=9
        verb_rej=19
! input observations
        CALL input_obs(obsdir,astnam(k),.false.,error_model,defobs             &
     &        ,obs,obsw,m,iun_log,change)
        IF(.not.defobs)THEN
           WRITE(*,*) "error in reading observations"
           STOP
        ENDIF
! name for output file
        CALL filnam(obsdir,astnam(k),'rwo',rwofil,le) 
! differential corrections  
        CALL coo_cha(el0(k),'EQU',el0(k),fail_flag)
        IF(fail_flag.ge.5)THEN
           WRITE(*,*) 'Impossible to convert to EQU'
           STOP
        ENDIF
        CALL tee(iun_log,' >>>>>>>DIFFERENTIAL CORRECTIONS=') 
        CALL fdiff_cor(.true.,0,defobs,.true.,ok,cov0,el0(k),m,obs,obsw, &
     &        usedobs,rwofil,elc,unc,csino,delno,rmsh,succ)
        IF(.not.succ)THEN
           WRITE(*,*) "error in differential corrections"
           STOP
        ENDIF
        el0(k)=elc
        CALL rmsp(astnam(k),ln)
        CALL write_elems(el0(k),astnam(k),'ML',UNC=unc,FILE=astnam(k)(1:ln)//".eq0")
        CALL coo_cha(el0(k),'KEP',el0(k),fail_flag)
     ENDIF
     IF(iptf1.LE.nout) THEN 
        t=teph1+dteph*(iptf1-1)
! propagation
        CALL pro_ele(el0(k),t,el1(k,iptf1))
        DO it=iptf1+1,nout 
           t=t+dteph
           CALL pro_ele(el1(k,iptf1),t,el1(k,it))
        ENDDO
        iptb1=iptf1-1 
     ELSE 
        iptb1=nout 
     END IF
! Backward integration
     CALL set_restart(.true.)
     IF(iptb1.GE.1) THEN 
        t=teph1+dteph*(iptb1-1)
        CALL pro_ele(el0(k),t,el1(k,iptb1))  
        DO it=iptb1-1,1,-1 
           t=t-dteph
           CALL pro_ele(el0(k),t,el1(k,it))
        ENDDO
     END IF
  ENDDO
! output on epehemerides files
!   WRITE(1,REC=nhl+iptf1) ((el1(k,iptf1)%coord(i),i=1,6),k=1,nbout)
!   WRITE(18,112) astnam(k),el1(k,iptf1)%t,el1(k,iptf1)%coord(1:6)
!   DO it=iptf1+1,nout
!      WRITE(1,REC=nhl+it) ((el1(k,it)%coord(i),i=1,6),k=1,nbout)
!      WRITE(18,112) astnam(k),el1(k,it)%t,el1(k,it)%coord(1:6)
!   ENDDO
!   WRITE(1,REC=nhl+iptb1) ((el1(k,iptb1)%coord(i),i=1,6),k=1,nbout)
!   WRITE(18,112) astnam(k),el1(k,iptb1)%t,el1(k,iptb1)%coord(1:6)
!   DO it=iptb1-1,1,-1
!     WRITE(1,REC=nhl+it) ((el1(k,it)%coord(i),i=1,6),k=1,nbout)
!     WRITE(18,112) astnam(k),el1(k,it)%t,el1(k,it)%coord(1:6)
!   ENDDO
112 FORMAT(A,7(1X,E18.10))
! alternative output
  DO it=1,nout
     WRITE(1,REC=nhl+it) ((el1(k,it)%coord(i),i=1,6),k=1,nbout)
!     DO k=1,nbout
!        WRITE(18,112) astnam(k),el1(k,it)%t,el1(k,it)%coord(1:6)
!     END DO
  ENDDO
! close close-approach file                                             
  IF(numcla.gt.0)THEN 
     CALL filclo(iuncla,' ') 
  ELSE 
     CALL filclo(iuncla,'DELETE') 
  ENDIF
! close err file (remove if empty)
  IF(numerr.eq.0)THEN
     CALL filclo(ierrou,'DELETE')
  ELSE 
     CALL filclo(ierrou,' ')
  ENDIF
  CALL filclo(iun_log,' ')
! close output file
  CLOSE(1) 

! done file
  CALL filopn(iundone,run(1:lr)//'.done','unknown')
  CALL filclo(iundone,' ')

  CONTAINS
       
! RDOPTB2
  SUBROUTINE rdoptb2(nbout,el0,astmas,astnam,t1,t2,dt, & 
&               diffcor_bineph,error_model,obsdir) 
      USE reference_systems 
      USE fund_const 
      USE orbit_elements
      USE dyn_param, ONLY:ndimx
      IMPLICIT NONE 
      INCLUDE 'dim.h90'
! interface
      INTEGER, INTENT(OUT) :: nbout ! number of asteroids to be used
      TYPE(orbit_elem), INTENT(OUT), DIMENSION(nbmx) :: el0 ! elements  
      DOUBLE PRECISION, INTENT(OUT) :: astmas(nbmx) ! masses
      CHARACTER*30, INTENT(OUT), DIMENSION(nbmx) :: astnam !names
      DOUBLE PRECISION, INTENT(OUT) :: t1,t2,dt ! beginning, end, step
      CHARACTER*100, INTENT(OUT) ::  obsdir
      CHARACTER*20, INTENT(OUT) :: error_model
      LOGICAL, INTENT(OUT) :: diffcor_bineph
! end interface      
! local variables
! Max number of initial condition files                                 
      INTEGER, PARAMETER :: niflx=10 
      CHARACTER*100 incfil(niflx),tmp 
      CHARACTER name1*30,eltype*3,rsys*10,epoch*10,scale*3 
      CHARACTER*30 astn1(nbmx),plnam(nbjx) 
      DOUBLE PRECISION elem(6),cove(ndimx,ndimx),nore(ndimx,ndimx),hmag,gmag,enne 
      DOUBLE PRECISION t0,t0t 
      DOUBLE PRECISION mass,sec,sece 
      INTEGER nifl,nbl,i,k,nbadd,kr,ln,lf,i1,i2,mjd,mjde, nbm,nba,nbt
      INTEGER fail_flag
      LOGICAL found,fail1,fail,defcov,defnor,end,loaded(nbmx),first 
      CHARACTER*60 comment
      LOGICAL ireq
!                                                         
      INTEGER lench 
      EXTERNAL lench 
! read options contained in jobname.bop
      fail=.false. 
      CALL rdmcha('input_files.','incond',incfil,nifl,niflx,'niflx',    &
     &            .true.,found,fail1,fail)                              
      CALL rdmcha('bineph.','eph_obj',astnam,nbout,nbmx,'nbmx',.true.,  &
     &            found,fail1,fail)                                     
      CALL rdmcha('bineph.','add_obj',astn1,nbadd,nbmx,'nbmx',.false.,  &
     &            found,fail1,fail)                                     
      IF(.NOT.found) nbadd=0 
      IF(fail) STOP '**** rdoptb: abnormal end 1 ****' 
! Limits and stepsize of binary ephemeris                               
      CALL rdntim('bineph.epoch.','start',tmp,mjd,sec,scale,.true.,     &
     &            found,fail1,fail)                                     
      IF(.NOT.fail1) THEN 
          CALL cnvtim(mjd,sec,scale,mjde,sece,'TDT') 
          t1=mjde+sece/86400.d0 
      END IF 
      CALL rdntim('bineph.epoch.','end',tmp,mjd,sec,scale,.true.,       &
     &            found,fail1,fail)                                     
      IF(.NOT.fail1) THEN 
          CALL cnvtim(mjd,sec,scale,mjde,sece,'TDT') 
          t2=mjde+sece/86400.d0 
      END IF 
      CALL rdnrea('bineph.','step',dt,.true.,found,fail1,fail) 
      IF(fail) STOP '**** rdoptb: abnormal end 2 ****' 
! total number of bodies 
      nbm=nbout+nbadd 
      CALL chkpdf(nbm,nbmx,'nbmx') 
      nba=0 
      nbt=nbm 
      DO 1 i=1,nbadd 
      astnam(nbout+i)=astn1(i) 
    1 END DO 
      DO 2 i=1,nbm 
      loaded(i)=.false. 
    2 END DO 
! Check for objects included more than once                             
      first=.true. 
      DO 10 i1=1,nbm-1 
         DO 11 i2=i1+1,nbm 
            IF(astnam(i1).EQ.astnam(i2)) THEN 
               IF(first) THEN 
                  WRITE(*,104) 
                  first=.false. 
               END IF
               ln=lench(astnam(i2)) 
               WRITE(*,103) astnam(i2)(1:ln) 
               fail=.true. 
            END IF
  104 FORMAT(' ERROR: the following objects are included more than',    &
     &       ' once'/                                                   &
     &       '        among initial conditions:')                       
  103 FORMAT(10X,A) 
11       END DO
10    END DO
      IF(fail) STOP '**** rdoptb2: abnormal end ****' 
! read initial conditions
      nbl=0 
      DO 4 i=1,nifl 
! open one of the files
         CALL oporbf(incfil(i),0) 
3        CONTINUE 
         CALL rd_orb(name1,elem,eltype,t0t,cove,defcov,nore,defnor,         &
     &          hmag,gmag,mass,rsys,epoch,kr,end)   
! what to do with rsys,epoch????
         IF(rsys.ne.'ECLM'.or.epoch.ne.'J2000')THEN
            WRITE(*,*)' rdoptb2: other reference systems not handled ', rsys,'  ', epoch
            STOP
         ENDIF
         IF(end) GOTO 6 
         DO 5 k=1,nbm 
            IF(loaded(k)) CYCLE 
            IF(name1.EQ.astnam(k)) THEN 
! store mass
               astmas(k)=mass  
               loaded(k)=.true. 
               nbl=nbl+1 
! store in array of data type elements
               el0(k)=undefined_orbit_elem
! copy into output elements
               el0(k)%coo=eltype
               el0(k)%coord=elem 
               el0(k)%t=t0t
               IF(eltype.eq.'ATT')THEN
! input obscod
                  WRITE(*,*) 'rdoptb2: no attributable elements'
                  STOP
               END IF
! absolute magnitude
               IF(hmag.lt.-10.d0)THEN
                  el0(k)%mag_set=.false.
               ELSE
                  el0(k)%mag_set=.true.
                  el0(k)%h_mag=hmag
                  el0(k)%g_mag=gmag
               ENDIF
! force to keplerian
               IF(el0(k)%coo.ne.'KEP')THEN
                  CALL coo_cha(el0(k),'KEP',el0(k),fail_flag)
               ENDIF
               IF(nbl.EQ.nbm) THEN 
                  CALL clorbf 
                  GOTO 7 
               END IF
            END IF
  100 FORMAT(' Epoch of initial conditions:',F13.5,' (MJD, TDT)') 
  101 FORMAT(' ERROR: object ',A,' (file ',A,', record',I6,')'/         &
     &       '        has orbital elements at a different epoch'/       &
     &       '        T0 =',F13.5,' (MJD, TDT)')                        
5        END DO
         GOTO 3 
6        CONTINUE 
         CALL clorbf 
4     END DO
      WRITE(*,102) 
  102 FORMAT(' ERROR: no initial conditions found for the following',   &
     &       ' objects:')                                               
      DO i=1,nbm 
         IF(.NOT.loaded(i)) THEN 
            ln=lench(astnam(i)) 
            WRITE(*,103) astnam(i)(1:ln) 
         END IF
      END DO
      STOP '**** rdoptb: abnormal end ****' 
7     CONTINUE 

! flag to decide whether or not perform differential corrections 
      ireq=.false.
      diffcor_bineph=.false.
      comment='differential corrections in bineph'
      CALL input_log_opt('bineph','diffcor_bineph',diffcor_bineph,ireq,&
           found,comment,iun_log)
! observations directory
      obsdir='.' !defaults as the input directory
      comment='observation directory'
      CALL input_cha_opt('bineph','obsdir',obsdir,ireq,found,comment,iun_log)
! error model
      error_model=' ' !defaults as the input directory
      comment='error model'
      CALL input_cha_opt('bineph','error_model',error_model,ireq,&
           found,comment,iun_log)

   END SUBROUTINE rdoptb2

                                                                        
END PROGRAM bineph2

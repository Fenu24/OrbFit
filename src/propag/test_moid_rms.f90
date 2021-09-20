! *********************************************************
! TEST PROGRAM for the COMPUTATION of the SIGNED 
! LOCAL MINIMA with their RMS
! written by  G.F. Gronchi  (January 2007)
! *********************************************************
PROGRAM test_moid_rms
  USE orbit_elements
  USE fund_const
  USE output_control
  USE critical_points
  USE dyn_param
  IMPLICIT NONE
  INTEGER, PARAMETER :: nastmax=10000
! for read_elems
  INTEGER :: h
  CHARACTER*80 :: catnam
  INTEGER :: iunin
  INTEGER :: astnum
! ==============OPTIONS============================
  CHARACTER(LEN=6) :: progna
  CHARACTER(LEN=80) :: run
  INTEGER :: iunrms
! ================PLANET===========================
  TYPE(orbit_elem) :: elem1
  TYPE(orb_uncert) :: unc1
  DOUBLE PRECISION :: pd_earth ! perihelion distance
! =================================================
  INTEGER :: j,m,enn,o,p,q
! ===============ASTEROID/COMET====================
  CHARACTER(LEN=9) :: name ! asteroid name 
  LOGICAL ::  eof, err
  CHARACTER(LEN=60) file,elefil,eledir ! file names 
  LOGICAL ::  chk_der
  INTEGER :: le
  LOGICAL :: neodys_elem
  CHARACTER(LEN=3) :: coord ! coord type
  LOGICAL, DIMENSION(3,nminx) :: comp_flag
  DOUBLE PRECISION, DIMENSION(6,6) :: dee 
  INTEGER :: nlsloc
  INTEGER :: lsloc(ndyx) ! solve-for dynamical
! -------------------------------------------------
  TYPE(orbit_elem) :: elem2
  TYPE(orb_uncert) :: unc2
! ================== MOID_RMS ===========================
  DOUBLE PRECISION,DIMENSION(nminx) :: dmintil
  DOUBLE PRECISION,DIMENSION(nminx) :: dmintrms
  DOUBLE PRECISION,DIMENSION(5,nminx) :: ddmintdel2
  DOUBLE PRECISION,DIMENSION(3,nminx) :: car1min,car2min
  DOUBLE PRECISION,DIMENSION(nminx) :: detH,detHrms
  DOUBLE PRECISION,DIMENSION(nminx) :: sint1t2,taurms
  DOUBLE PRECISION :: sinmutI,sinmutIrms
  DOUBLE PRECISION,DIMENSION(nminx) :: f1min,f2min
  DOUBLE PRECISION :: min_val,max_val
  DOUBLE PRECISION :: min_detH,max_detH
  DOUBLE PRECISION :: min_val_app_p,max_val_app_p
  DOUBLE PRECISION :: min_val_app_m,max_val_app_m
  INTEGER :: nummin,ncri ! number of relative minima and critical points
  DOUBLE PRECISION :: dnodp,dnodm,amoidp,amoidm
  DOUBLE PRECISION :: k1p,k2p,k1m,k2m
  DOUBLE PRECISION :: f1nodp,f2nodp,f1nodm,f2nodm!true anom. at mutual nodes
  DOUBLE PRECISION :: dnodprms,dnodmrms
  DOUBLE PRECISION :: amoidprms,amoidmrms
  INTEGER :: loop,ler
  LOGICAL :: ini2,cov2
  DOUBLE PRECISION :: atmp,etmp,itmp,Omtmp,omegtmp ! auxiliary
! =====================================================================
! options
  progna='compmo' 
  run='test_moid_rms'
  CALL compop
  CALL rmsp(run,ler)
  CALL filopn(ierrou,run(1:ler)//'.err','unknown')
  astnum = 0 
! =========== ENTER MAIN LOOP ===========================           
  DO  loop=1,nastmax 
     astnum = astnum + 1 
     ! WRITE(*,*) 'Asteroid/Comet Name?' 
5    READ(*,100,end=33) name 
!     write(*,*) name
100  FORMAT(a9) 
     ini2=.false. 
     cov2=.false. 
! =============INPUT ASTEROID ELEMENT SET================      
     CALL filnam(eledir,name,'eq1',elefil,le) 
     INQUIRE(file=elefil(1:le),exist=neodys_elem) 
       ! if data available
     IF(.not.neodys_elem)THEN
        WRITE(iun_log,*)'!!WARNING! WARNING! WARNING! WARNING!!' 
        WRITE(iun_log,*)'Asteroid ',name,' number ',astnum,&
             & ' initial conditions not found.'             
        STOP
     ENDIF

     CALL read_elems(elem2,name,eof,nlsloc,err,lsloc,UNC=unc2,FILE=elefil)
     write(*,*)'nlsloc',nlsloc

     If(.not.eof)THEN
        ini2=.true.
     ELSE
        WRITE(ierrou,*)'!!WARNING! WARNING! WARNING! WARNING!!' 
        WRITE(ierrou,*)'Asteroid ',name,' number ',astnum,&
             & ' initial conditions not found.' 
        numerr=numerr+1
        CYCLE
!        STOP
     ENDIF
! failure case: either .eq1 does not exist, or it is rotten             
     IF(.not.unc2%succ)THEN 
        WRITE(ierrou,*)'Asteroid ',name,' number ',astnum,&
             & ' covariance matrix not found.'                  
        WRITE(ierrou,*)'skipping moid computation'
        numerr=numerr+1
        CYCLE
!        STOP 
     ENDIF
! Earth orbital elements
     elem1=undefined_orbit_elem
     CALL undefined_orb_uncert(6,unc1)
     unc1%g=0.d0 ! to be sure Earth has negligible uncertainty
     CALL earth(elem2%t,elem1%coord)
     elem1%coo='EQU'
     elem1%t=elem2%t
     pd_earth=elem1%coord(1)*(1-SQRT(elem1%coord(2)**2+elem1%coord(3)**2))
! ------------------------------------------------------------------
     CALL dmintil_rms(elem1,elem2,nummin,dmintil,car1min,car2min,&
     & unc1,unc2,dmintrms,ddmintdel2,comp_flag,detH,detHrms,sint1t2,taurms,&
     & sinmutI,sinmutIrms,chk_der)
! test skipping optional variables
!     CALL dmintil_rms(elem1,elem2,nummin,dmintil)
! ------------------------------------------------------------------
111  FORMAT(5(a9,4x))
112  FORMAT(a9,4x,4(es11.4,1x))
! --- write output ---
     WRITE(*,111)name,'Nom.Val.','RMS','Min.Val.','Max.Val.'
     DO j=1,nummin
        write(*,*)'j=',j
        min_val=dmintil(j)-3.d0*dmintrms(j)
        max_val=dmintil(j)+3.d0*dmintrms(j)
        min_detH=detH(j)-3.d0*detHrms(j)
        max_detH=detH(j)+3.d0*detHrms(j)
        WRITE(*,112)'loc.min.',dmintil(j),dmintrms(j),min_val,max_val
!        WRITE(*,*)'loc.min.',dmintil(j)
        WRITE(*,112)'detH',detH(j),detHrms(j),min_detH,max_detH
        WRITE(*,*)'CAR Earth:',car1min(1:3,j)
        WRITE(*,*)'CAR asteroid:',car2min(1:3,j)
!        write(*,*)'comp_flag=',comp_flag(1:3,j)
        IF(.not.comp_flag(1,j))THEN
           WRITE(ierrou,*)'sinmutI may be zero: min=',sinmutI-3.d0*sinmutIrms,&
                & 'max=',sinmutI+3.d0*sinmutIrms 
           numerr=numerr+1
        ENDIF
        IF(.not.comp_flag(2,j))THEN
           WRITE(ierrou,*)'detH may be zero: min=',detH(j)-3.d0*detHrms(j),&
                & 'max=',detH(j)+3.d0*detHrms(j) 
           numerr=numerr+1
        ENDIF
        IF(.not.comp_flag(3,j))THEN
           WRITE(ierrou,*)'sint1t2 may be zero: min=', &
                & sint1t2(j)-3.d0*taurms(j),&
                & 'max=',sint1t2(j)+3.d0*taurms(j) 
           numerr=numerr+1
        ENDIF
     ENDDO     
  ENDDO
  WRITE(*,*)' number of orbits in input ', astnum !norb
33 CONTINUE  
  CALL filclo(ierrou,' ')
CONTAINS
! *****************************************************************  
  SUBROUTINE compop 
    IMPLICIT NONE
    INTEGER            :: le,iunout
    CHARACTER(LEN=100) :: file 
    LOGICAL            :: ireq,found
    CHARACTER(LEN=60)  :: comment 
! =============================  
    CALL initopt(progna,run,'mop')  
! read option for physical model and integration method                 
    CALL rmodel(1)    
    ! Output files: for control                 
    file=run//'.mou' 
    CALL rmsp(file,le) 
    call filopn(iunout,file(1:le),'UNKNOWN') 
! =============WHERE TO FIND ASTEROID ELEMENTS===================    
    ireq=.true. 
    comment='asteroid elements directory'
    CALL input_cha_opt(progna,'eledir',eledir,ireq,found,comment,iunout)
! =============CHECK DERIVATIVES===================
    ireq=.true.
    comment='check derivatives'
    CALL input_log_opt(progna,'chk_der',chk_der,ireq,found,comment,iunout)
  END SUBROUTINE compop
  
END PROGRAM test_moid_rms

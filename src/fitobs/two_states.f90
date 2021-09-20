MODULE two_states
  USE orbit_elements
  USE astrometric_observations
  IMPLICIT NONE
  PRIVATE
!============= asteroid names (18 CHARACTERs)========
  CHARACTER*18 astna0,astnap
  CHARACTER*25 astnaj,astnac 
! observation numbers: maximum and actual                               
  INCLUDE 'parobx.h90' 
  INTEGER m,mp,mall,mc 
! new data types
  TYPE(ast_obs),DIMENSION(nobx),SAVE :: obs,obsc
  TYPE(ast_wbsr),DIMENSION(nobx) :: obsw,obswc
  CHARACTER*20 ::  error_model ! weighing model
! successful input flags, for observations                              
  LOGICAL obs0,obsp,obstwo,obsflag
! =====state variables: first arc, second arc, joint, current========   
! equinoctal orbital elements, epoch times
  TYPE(orbit_elem) :: el0,elp,el,elc
 !  successful input flags, for elements                                 
  LOGICAL ini0,inip,initwo,inic 
! ===== differential corrections output ==============                  
! normal, covariance matrices, existence flags for covariance           
  TYPE(orb_uncert) :: unc0,uncp,unc,uncc
  LOGICAL cov0,covp,covtwo,covc ! ??? use unc0%succ
! norms of correction, of residuals                                     
  DOUBLE PRECISION delno0,delnop,delnor,delnoc 
  DOUBLE PRECISION csino0,csinop,csinor,csinoc
  DOUBLE PRECISION rmsh0,rmshp,rmsh,rmshc 
  INTEGER iob0,iobp,iobtwo,iobc ! number of obs used
! ============I/O =====================                        
  CHARACTER*80 elefi0,elefip,eletwo ! elements file names 
  INTEGER iunel0,iunelp,iunelt,iunelc ! elements logical units
  CHARACTER*60 rwofi0,rwofip,rwotwo,rwofil,rwofic ! weight file names   

PUBLIC el0,elp,el,elc,unc0,uncp,unc,uncc,ini0,inip,initwo,inic,cov0,covp,covc,covtwo
PUBLIC delno0,delnop,delnor,delnoc,csino0,csinop,csinor,csinoc,rmsh0,rmshp,rmsh,rmshc
PUBLIC m,mp,mall,mc,obs,obsw,obsc,obswc,obs0,obsp,obstwo,obsflag,iob0,iobp,iobtwo,iobc,error_model,nobx
PUBLIC astna0,astnap,astnac,astnaj
PUBLIC elefi0,elefip,eletwo,iunel0,iunelp,iunelt,iunelc,rwofi0,rwofip,rwotwo,rwofil,rwofic
! subroutines
PUBLIC sta_cop, obs_cop,set_state_def,orb_sel2,orb_sel


CONTAINS 
! =====================================================                 
! ORB_SEL                                                               
SUBROUTINE orb_sel2(batch,isel) 
  USE util_suit
  USE output_control
  LOGICAL, INTENT(IN) :: batch !preselected isel
! output: selection 1=arc1 2=arc2 3=joint solution                      
  INTEGER,INTENT(INOUT) :: isel 
! characters for menu                                                   
  CHARACTER*20 menunam
! two modes of operations:
  IF(batch)THEN
! controls????
     GOTO 2
  ENDIF 
! cases without questions                                               
  IF(ini0.and.(.not.inip.and..not.initwo))THEN 
     CALL tee(iun_log,' ARC 1=')
     isel=1 
     GOTO 2 
  ELSEIF(inip.and.(.not.ini0.and..not.initwo))THEN 
     CALL tee(iun_log,' ARC 2=')
     isel=2 
     GOTO 2 
  ELSEIF(.not.ini0.and..not.inip.and..not.initwo)THEN 
     CALL tee(iun_log,' No initial conditions available=')
     isel=0 
     RETURN
  ENDIF
! cases in which there is a choice                                      
  menunam='orbsel' 
3 CALL menu(isel,menunam,3,'which orbit?=',                     &
     &      'arc 1=','arc 2=',                                  &
     &      'joint computed orbit=')
  IF(isel.eq.0)RETURN
  IF(isel.eq.1.and..not.ini0)THEN
     WRITE(*,*)' intial conditions for arc 1 not available'
     GOTO 3
  ELSEIF(isel.eq.2.and..not.inip)THEN
     WRITE(*,*)' intial conditions for arc 2 not available'
     GOTO 3
  ELSEIF(isel.eq.3.and..not.initwo)THEN
     WRITE(*,*)' intial conditions for identification not available'
     GOTO 3
  ENDIF
2 CONTINUE
  CALL sta_cop(1,isel)
END SUBROUTINE orb_sel2
! ORB_SEL                                                               
SUBROUTINE orb_sel(batch,isel) 
  USE util_suit
  USE output_control
  LOGICAL, INTENT(IN) :: batch !preselected isel
! output: selection 1=arc1 2=arc2 3=joint solution                      
  INTEGER,INTENT(INOUT) :: isel 
! characters for menu                                                   
  CHARACTER*20 menunam
! cases in which there is a choice                                      
  menunam='orbsel' 
  CALL menu(isel,menunam,3,'which orbit?=',                     &
     &      'arc 1=','arc 2=',                                  &
     &      'joint computed orbit=')
END SUBROUTINE orb_sel
! =====================================================================
! SET_STATE_DEF
! sets defaults for state and observations, before any input
! =====================================================================
SUBROUTINE set_state_def
! setting of logical flags: nothing is available at the beginning       
  obs0=.false. 
  obsp=.false. 
  obstwo=.false. 
  ini0=.false. 
  inip=.false. 
  initwo=.false. 
  cov0=.false. 
  covp=.false. 
  covtwo=.false. 
! residual and correction norms also undefined
  csino0=-1.d0
  csinop=-1.d0
  csinor=-1.d0
  csinoc=-1.d0
  delno0=-1.d0
  delnop=-1.d0
  delnor=-1.d0
  delnoc=-1.d0
! setting to zero of counters for data not yet available                
  m=0 
  mp=0 
  mall=0 
  mc=0
  obs0=.false.
  obsp=.false.
  obstwo=.false.
  obsflag=.false.
  iob0=0 
  iobp=0 
  iobtwo=0 
  el0=undefined_orbit_elem
  el0%t=0.d0
  elp=undefined_orbit_elem
  elp%t=0.d0
  el=undefined_orbit_elem
  el%t=0.d0
  elc=undefined_orbit_elem
  elc%t=0.d0
END SUBROUTINE set_state_def
! ===================================================================== 
! STACOP                                                                
! ===================================================================== 
! copy state variables to/from current state                                 
! ================INTERFACE=========================================    
subroutine sta_cop(icop,iarc)
  integer, intent(IN) ::  icop,iarc
! =================END INTERFACE==================================== 
  CHARACTER*40 longname 
  INTEGER le  
! function: copy to/from
  IF(icop.eq.1)THEN ! copy from iarc to current
     IF(iarc.eq.1)THEN                                                
        elc=el0 
        inic=ini0
        covc=cov0
        IF(covc)uncc=unc0
        csinoc=csino0 
        delnoc=delno0 
        rmshc=rmsh0
        astnac=astna0
     ELSEIF(iarc.eq.2)THEN
        elc=elp
        inic=inip
        covc=covp
        IF(covc)uncc=uncp
        csinoc=csinop 
        delnoc=delnop 
        rmshc=rmshp
        astnac=astnap
     ELSEIF(iarc.eq.3)THEN
        elc=el
        inic=initwo
        covc=covtwo
        IF(covc)uncc=unc
        csinoc=csinor 
        delnoc=delnor 
        rmshc=rmsh
        longname=astna0//'='//astnap
        CALL rmsp(longname,le)
        astnac=longname(1:le)
     ELSE
        WRITE(*,*)' sta_cop: iarc not understood ',iarc
        STOP
     ENDIF
  ELSEIF(icop.eq.2)THEN! copy from current to iarc
      IF(iarc.eq.1)THEN                                                
        el0=elc 
        ini0=inic
        IF(covc)unc0=uncc
        cov0=covc
        csino0=csinoc 
        delno0=delnoc 
        rmsh0=rmshc
     ELSEIF(iarc.eq.2)THEN
        elp=elc
        inip=inic
        IF(covc)uncp=uncc
        covp=covc
        csinop=csinoc 
        delnop=delnoc 
        rmshp=rmshc
     ELSEIF(iarc.eq.3)THEN
        el=elc
        initwo=inic
        IF(covc)unc=uncc
        covtwo=covc
        csinor=csinoc 
        delnor=delnoc 
        rmsh=rmshc
     ELSE
        WRITE(*,*)' sta_cop: iarc not understood ',iarc
        STOP
     ENDIF
  ELSE 
     write(*,*)' sta_cop: copying not known, icop=',icop 
     STOP
  ENDIF
END SUBROUTINE sta_cop
! ===================================================================== 
! OBS_COP                                                                
! ===================================================================== 
! copy observations to/from current observations
! ================INTERFACE=========================================    
SUBROUTINE obs_cop(icop,iarc)
  INTEGER, INTENT(IN) ::  icop,iarc
! =================END INTERFACE====================================    
! function: copy to/from       
  IF(mall.ne.m+mp)THEN
     WRITE(*,*)' obs_cop: wrong obs. count ', m,mp,mall
     STOP
  ENDIF
  IF(icop.eq.1)THEN ! copy observations from iarc to current dir.
     IF(iarc.eq.1)THEN
        obsflag=obs0
        IF(obsflag)THEN
           mc=m
           obsc(1:mc)=obs(1:m) 
           obswc(1:mc)=obsw(1:m)
           iobc=iob0
        ENDIF
     ELSEIF(iarc.eq.2)THEN
        obsflag=obsp
        IF(obsflag)THEN
           mc=mp
           obsc(1:mc)=obs(m+1:m+mp) 
           obswc(1:mc)=obsw(m+1:m+mp)
           iobc=iobp
        ENDIF
     ELSEIF(iarc.eq.3)THEN
        obsflag=obstwo
        IF(obsflag)THEN
           mc=mall
           obsc(1:mc)=obs(1:mall) 
           obswc(1:mc)=obsw(1:mall)
           iobc=iobtwo
        ENDIF
     ELSE
        WRITE(*,*)'obs_cop: iarc not understood ', iarc
        STOP
     ENDIF
  ELSEIF(icop.eq.2)THEN ! copy residuals from current to iarc
     IF(iarc.eq.1)THEN
        IF(mc.ne.m)STOP 
        obsw(1:m)=obswc(1:mc)
        iob0=iobc
     ELSEIF(iarc.eq.2)THEN
        IF(mc.ne.mp)STOP 
        obsw(m+1:mall)=obswc(1:mc)
        iobp=iobc
     ELSEIF(iarc.eq.3)THEN
        IF(mc.ne.mall)STOP 
        obsw(1:mall)=obswc(1:mc)
        iobtwo=iobc
     ENDIF
  ELSE 
     write(*,*)' obs_cop: copying not known, icop=',icop 
     STOP
  ENDIF
 
END SUBROUTINE obs_cop

END MODULE two_states

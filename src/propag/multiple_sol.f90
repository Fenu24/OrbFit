! ========MODULE multiple_sol================================
!
MODULE multiple_sol  
  USE output_control
  USE fund_const
  USE astrometric_observations
  USE orbit_elements
  USE least_squares
  USE dyn_param
  USE multi_store
  IMPLICIT NONE
  PRIVATE
  
  ! LIST OF PUBLIC ENTITIES
  ! PUBLIC ROUTINEs
  PUBLIC :: f_multi, nmulti, nmulti_gen, mult_input, lin_multi, lin_multi_gen, multi_gen
  PUBLIC :: fmuobs, fmupro, fmuplo, outmul, prop_sig, step_fit2
  ! former lov_int.mod
  PUBLIC :: lovinit, lovobs, lovmagn, lovinterp, lovdensify, lovinterp3
  ! PUBLIC data
  ! maximum number of multiple solutions: moved to multi_store
  
  ! first, last, reference solution                            
  INTEGER, PUBLIC ::  imip,imim
  ! multiple solution arrays                                                
  TYPE(orbit_elem),                 PUBLIC :: elm(1:mulx), elm_loc(mulx)
  DOUBLE PRECISION,                 PUBLIC :: csinom(mulx),delnom(mulx)
  TYPE(orb_uncert),DIMENSION(mulx), PUBLIC :: unm
  DOUBLE PRECISION,                 PUBLIC :: dpm(ndyx,mulx), dpm_loc(ndyx,mulx)
  REAL(KIND=qkind)                        :: vvv(ndimx) ! TP scatter vector
  DOUBLE PRECISION,                 PUBLIC :: tdt_cat ! common epoch time
  DOUBLE PRECISION,                 PUBLIC :: v_inf(mulx)
  DOUBLE PRECISION,                 PUBLIC :: moid_m(mulx), dnp_m(mulx), dnm_m(mulx) ! moid and nodal dist
  INTEGER,DIMENSION(mulx),          PUBLIC :: imul ! imul(j) is the index of the j-th computed multiple solution 
  !                                                  (if the first VA in the ctc file has name va_N then imul(1)=N

  ! sigma_Q, next step in sigma_LOV
  DOUBLE PRECISION,DIMENSION(mulx), PUBLIC :: sigq, deltasigmav
  ! options for probability sampling: step in IP, max step in sigm_lov, 
  !                                   max extent in sigma_lov (moved to multi_store)
  DOUBLE PRECISION, PUBLIC :: prob_step,sigma_step_max
  ! ================options for use of different LOVs and samplings==========
  LOGICAL,          PUBLIC :: lin_lov      ! use of linear approximation to LOV
  LOGICAL,          PUBLIC :: no_outcov    ! output only orbits, no covariance, also in .ctc file
  LOGICAL,          PUBLIC :: scatterplane ! use of scatter TP metric for weak direction
  ! interval containing the scattering encounter
  ! WARNING : if there are many close app in the interval a mess occurs
  DOUBLE PRECISION, PUBLIC :: before_scatter,after_scatter
  
CONTAINS                                                                     
! =======================================
!  LIN_MULTI
! =======================================         
! minimum interface/computation version...problems with sigma_Q
!     CALL lin_multi(elc,uncert,sigma,imult,nd)
! maximum interface version
!     CALL lin_multi(batch,obsc,inic,ok,covc,elc,uncert,csinor,delnor,mc,obs,obsw,sigma,imult1,nd)
SUBROUTINE lin_multi(batch,el0,unc0,sigma,imult,nd,iwdir,m,obs,obsw,csinor)
  USE propag_state
  USE tp_trace
  USE close_app, ONLY: kill_propag
  USE semi_linear
! ================INPUT===========================                      
! initial conditions: epoch, elements, abs. magnitude, gmag  
  LOGICAL, INTENT(IN)          :: batch ! batch/interactive
  TYPE(orbit_elem), INTENT(IN) :: el0 ! 
! normal and covariance matrices, norms of residuals, of last correction
  TYPE(orb_uncert), INTENT(IN) :: unc0 ! 
  INTEGER, INTENT(IN)          :: m ! no. observations
  TYPE(ast_obs), DIMENSION(m) :: obs
  TYPE(ast_wbsr), DIMENSION(m) :: obsw 
  DOUBLE PRECISION, INTENT(IN) :: csinor 
! ===============INPUT============================
  DOUBLE PRECISION, INTENT(IN) :: sigma ! max of sigma along the line
  INTEGER, INTENT(IN)          :: imult ! number of alternate solutions
  INTEGER, INTENT(IN)          :: nd ! dimension of solve for parameters vector
  INTEGER, INTENT(IN)          :: iwdir ! unit for TP scatter vector output
! ===============OUTPUT=========================                
! ==============END INTERFACE=====================
  TYPE(tp_point)   tptr
  TYPE(orbit_elem) eltmp,elf          ! temporary orbit elements, final after scatter
  TYPE(orb_uncert) uncf               ! covariance after scatter           
  DOUBLE PRECISION :: wdir0(ndimx),sdir0 ! weak direction, rms along it
  DOUBLE PRECISION :: dn,delta_sigma     ! step in index, in sigma 
  DOUBLE PRECISION :: dp0(ndyx),dp(ndyx) ! non-grav parameters

! weak direction from scatterplane
  DOUBLE PRECISION :: gtp(2,2) ! covariance on target plane
  REAL(KIND=qkind) :: gtp_q(2,2) 
  REAL(KIND=qkind) :: ceicel(ndimx-2,2), v(ndimx,ndimx), b(2,2) ! for slinel
  REAL(KIND=qkind), DIMENSION(2) :: wtp, xy, wtp1, wtp2, wtp_ell
  REAL(KIND=qkind), DIMENSION(ndimx-2) :: hmhs
  REAL(KIND=qkind), DIMENSION(ndimx):: vvv1, vvv2
  REAL(KIND=qkind) :: dtpdetq(1:ndimx,1:2),dtpdeq(1:2,1:ndimx)
  REAL(KIND=qkind), DIMENSION(ndimx,ndimx) :: c,g
  REAL(KIND=qkind) :: chi1,chi2,chi12,chi,dchi
  REAL(KIND=qkind), DIMENSION(ndimx) :: cv1,cv2,cv
! for call to fdiff_cor
  INTEGER nused,itmaxold
  LOGICAL succ, ok, cov0
  DOUBLE PRECISION delno0,rmsh0
  CHARACTER(60) rwofi0
! eigenvalues, workspace, transposed                                    
  DOUBLE PRECISION :: eigval(2),fv1(2),fv2(2)
  INTEGER :: ierr
  DOUBLE PRECISION :: sig(2), axes(2,2), dynm(ndyx,mulx),elmcoo(6,mulx)
  LOGICAl bizarre                     ! control of bizarre orbits
  LOGICAL falsok
  INTEGER :: iun,iundump,iunint, iunrep, iuncat, iuncatngr, iunctc
  INTEGER :: ngr_case
  INTEGER ipla
  DOUBLE PRECISION, DIMENSION(ndimx) :: units ! for scaling
  DOUBLE PRECISION :: tafter,tbefore           ! epoch time
  DOUBLE PRECISION :: v_infty, v_inf0             ! velocity at the boundary of Earth's sphere of influence
  INTEGER :: j,i,imi, np, nm, jc
  DOUBLE PRECISION :: alpha, tprmin
  DOUBLE PRECISION prscal ! function
! Outputs
  DOUBLE PRECISION :: deltasigma
  DOUBLE PRECISION, DIMENSION(mulx) :: sigmalov
! hack for catalog                                                      
  INTEGER :: le 
  CHARACTER(LEN=20) :: astna0 
! file names
  CHARACTER(LEN=160) :: catname,repname,ctcname,catngrname

  iun=abs(iun_log) 
  IF(verb_mul.ge.30)THEN 
     iundump=abs(iun_log) 
  ELSE 
     iundump=-abs(iun_log) 
  ENDIF
  IF(verb_mul.ge.20)THEN 
     iunint=abs(iun_log) 
  ELSE 
     iunint=-abs(iun_log) 
  ENDIF
! ===================================================================== 
! nominal stepsize (in the index space)                                 
  dn=1.d0/float(imult) 
! nominal stepsize (in sigma-parameter space)
  delta_sigma=dn*sigma
! ===================================================================== 
! nominal solution at the center of the list                            
  imi0=imult+1 
  elm(imi0)=el0 
  unm(imi0)=unc0
  csinom(imi0)=csinor
  DO j=1,dyn%ndp
    dpm(j,imi0)=dyn%dp(j)
  ENDDO
  ! orbital distance                                                      
  CALL nomoid(el0%t,elm(imi0),moid_m(imi0),dnp_m(imi0),dnm_m(imi0))
! ===================================================================== 
! weak direction and weak axis
  IF(scatterplane)THEN
! select propagation time                                               
     ipla=3 
     v_inf0=v_infty(el0) 
     CALL aftclo2v(ipla,el0%t,before_scatter,after_scatter,v_inf0,tbefore,tafter) 
! availability of covariance for TP analysis online                     
     CALL cov_avai(unc0,el0%coo,el0%coord) 
! reset close approach storage                                          
     njc=0 
     iplam=0 
! propagate to find linear confidence ellipse in scattering TP
     CALL pro_ele(el0,tafter,elf,unc0,uncf) 
! check success in finding the required target plane
     IF(njc.eq.0)THEN 
        WRITE(*,*)' lin_lov failed; no close approach' 
        STOP
     ELSEIF(iplam.ne.3)THEN 
        WRITE(*,*)' lin_lov failed; close app. planet ', iplam, tafter
        STOP
     ELSEIF(tp_store(njc)%tcla.gt.after_scatter.or.tp_store(njc)%tcla.lt.before_scatter)THEN   
        IF(kill_propag)THEN
           WRITE(*,*)' lin_lov failed: previous collision at ',tcla(njc)
        ELSE
           WRITE(*,*)' lin_lov failed; tcla not in interval ',tcla(njc) 
        ENDIF
        STOP
     ELSE 
        IF(kill_propag)THEN
           WRITE(*,*)' lin_lov: collision found at nominal orbit'
        ENDIF
     ENDIF
! load into tp_point data type
! select lower minimum                                                  
     IF(njc.eq.1)THEN 
        jc=1 
     ELSE 
        tprmin=1.d10 
        DO j=1,njc
           IF(tp_store(j)%d.lt.tprmin)THEN 
              jc=j 
              tprmin=tp_store(j)%d 
           ENDIF
        ENDDO
     ENDIF
! ===quadruple precision portion=====================
! conversion to quadruple precision
     dtpdetq(1:nd,1:2)=tp_store(jc)%dtpdet(1:nd,1:2)
     dtpdeq(1:2,1:nd)=TRANSPOSE(dtpdetq(1:nd,1:2))
     g(1:nd,1:nd)=unc0%g(1:nd,1:nd)
     c(1:nd,1:nd)=unc0%c(1:nd,1:nd)
! covariance on target plane
     gtp_q=MATMUL(MATMUL(dtpdeq(1:2,1:nd),g(1:nd,1:nd)),dtpdetq(1:nd,1:2))
! WARNING: still in double precision
     gtp=gtp_q
! eigenvalues, major semiaxis of ellipse on TP     
     CALL rs(2,2,gtp,eigval,1,axes,fv1,fv2,ierr)
! END WARNING
! conversion to quadruple precision
     wtp_ell=axes(:,2)*sqrt(eigval(2)) 
! find counterimage in initial conditions
     CALL slinel_q(dtpdetq(1:nd,1:2),g(1:nd,1:nd),c(1:nd,1:nd),ceicel(1:nd-2,1:2),b,v(1:nd,1:nd),nd)
!     CALL slinel(dtpdetq(1:nd,1:2),g(1:nd,1:nd),c(1:nd,1:nd),ceicel(1:nd-2,1:2),b,v(1:nd,1:nd),nd)
     xy=MATMUL(b,wtp_ell) ! g-g^*
     hmhs(1:nd-2)=MATMUL(ceicel(1:nd-2,1:2),xy)! h-h^*
     vvv1(1:nd)=MATMUL(v(1:nd,1:2),xy) ! R^T [g-g^* 0]
     wtp1=MATMUL(dtpdeq(1:2,1:nd),vvv1(1:nd))
     vvv2(1:nd)=MATMUL(v(1:nd,3:nd),hmhs(1:nd-2)) ! R^T [0 (h-h^*)]   
     wtp2=MATMUL(dtpdeq(1:2,1:nd),vvv2(1:nd))
     wtp=wtp1-wtp2
     vvv(1:nd)=vvv1(1:nd)-vvv2(1:nd) ! R^T [0 h-h^*]
! axial vector continued as true vector
     IF(el0%coo.eq.'EQU'.or.el0%coo.eq.'KEP'.or.el0%coo.eq.'COM'.or.el0%coo.eq.'COT')THEN
! if coo='EQU', 'KEP' a growing with sigma, if 'COM' q growing with sigma 
        IF(vvv(1).lt.0.d0) vvv(1:nd)=-vvv(1:nd)
     ELSEIF(el0%coo.eq.'ATT')THEN
! if coo='ATT' r growing with sigma
        IF(vvv(5).lt.0.d0) vvv(1:nd)=-vvv(1:nd)
     ELSEIF(el0%coo.eq.'CAR')THEN
! if coo='CAR' r growing with sigma
        IF(prscal(vvv(1:3),el0%coord(1:3)).lt.0.d0) vvv(1:nd)=-vvv(1:nd)
     ENDIF
! output direction of LOV in initial conditions space for use by resret
     WRITE(iwdir,100) vvv(1:nd)
 100 FORMAT(1P,10E32.22)
  ELSE
! compute line of variations                                            
     CALL weak_dir(unc0%g(1:nd,1:nd),wdir0,sdir0,iunint,el0%coo,el0%coord,units,nd)
     vvv(1:nd)=units(1:nd)*wdir0(1:nd)*sdir0
  ENDIF
! prepare for iterations
  IF(nd.gt.6)THEN
     DO j=1,nls
        dp0(j)=dyn%dp(ls(j))
     ENDDO
     dp(1:nls)=dp0(1:nls)
  ENDIF
  sigq(imi0)=0.d0
  itmaxold=itmax
! ===================================================================== 
! loop on positive sigma-parameter
  DO i=1,imult
     imi=imi0+i
     WRITE(*,*)' alternate solution no. ',imi 
     WRITE(iun,*)' alternate solution no. ',imi 
     elm(imi)=elm(imi-1)
     elm(imi)%coord=elm(imi)%coord+vvv(1:6)*delta_sigma
     IF(nd.gt.6)THEN
        dpm(:,imi)=dpm(:,imi-1)
        DO j=1,nls
           dp(j)=dp(j)+vvv(6+j)*delta_sigma
           dpm(ls(j),imi)=dp(j)
! for possible computation of residuals
           dyn%dp(ls(j))=dpm(ls(j),imi)
        ENDDO
     ENDIF
     unm(imi)=unm(imi0)
! orbital distance                                                      
     CALL nomoid(el0%t,elm(imi),moid_m(imi),dnp_m(imi),dnm_m(imi))
! computation of residuals
     IF(MOD(imi,100).eq.1)THEN
        cov0=.false.
        itmax=0
        CALL fdiff_cor(.true.,-1,.true.,.true.,ok,cov0,elm(imi),m,obs,obsw,nused,   &
     &     rwofi0,elm(imi),unm(imi),csinom(imi),delno0,rmsh0,succ)
        sigq(imi)=sqrt(abs(csinom(imi)**2-csinom(imi0)**2)*nused)
        IF(sigq(imi-1).lt.0.d0) sigq(imi)=-sigq(imi)
     ELSEIF(MOD(imi,100).eq.2.and.imi.ne.imi0+1)THEN
         sigq(imi)=sigq(imi-2)+delta_sigma*2
     ELSE
        sigq(imi)=sigq(imi-1)+delta_sigma
     ENDIF
  ENDDO
  imip=imi0+imult
! ===================================================================== 
! loop on negative sigma-parameter
  dp(1:nls)=dp0(1:nls)
  DO i=1,imult
     imi=imi0-i
     WRITE(*,*)' alternate solution no. ',imi 
     WRITE(iun,*)' alternate solution no. ',imi 
     elm(imi)=elm(imi+1)
     elm(imi)%coord=elm(imi)%coord-vvv(1:6)*delta_sigma
     IF(nd.gt.6)THEN
        dpm(:,imi)=dpm(:,imi+1)
        DO j=1,nls
           dp(j)=dp(j)-vvv(6+j)*delta_sigma
           dpm(ls(j),imi)=dp(j)
! for possible computation of residuals
           dyn%dp(ls(j))=dpm(ls(j),imi)
        ENDDO
     ENDIF
     unm(imi)=unm(imi0)
! orbital distance                                                      
     CALL nomoid(el0%t,elm(imi),moid_m(imi),dnp_m(imi),dnm_m(imi))
! computation of residuals
     IF(MOD(imi,100).eq.1)THEN
        cov0=.false.
        itmax=0
        CALL fdiff_cor(.true.,-1,.true.,.true.,ok,cov0,elm(imi),m,obs,obsw,nused,   &
     &     rwofi0,elm(imi),unm(imi),csinom(imi),delno0,rmsh0,succ)
        sigq(imi)=sqrt(abs(csinom(imi)**2-csinom(imi0)**2)*nused)
        IF(sigq(imi+1).lt.0.d0) sigq(imi)=-sigq(imi)
     ELSEIF(MOD(imi,100).eq.0.and.imi.ne.imi0-1)THEN
         sigq(imi)=sigq(imi+2)-delta_sigma*2
     ELSE
        sigq(imi)=sigq(imi+1)-delta_sigma
     ENDIF
  ENDDO
  imim=imi0-imult
  itmax=itmaxold
! catalog reference time
  tdt_cat=el0%t
  IF(batch)RETURN
! count cases with non-grav parameters
  ngr_case=0
! Object name
  CALL rmsp(name_obj,le)    
! Catalog (name_obj.cat)
  catname=name_obj(1:le)//'.cat'
  CALL filopn(iuncat,catname,'unknown')
  CALL wro1lh2(iuncat,'ECLM','J2000',el0%coo,'OEF2.0')
! single line format for non-grav param                     
  catngrname=name_obj(1:le)//'.cat.ngr'
  CALL filopn(iuncatngr,catngrname,'unknown') 
  CALL wro1lh2_ngr(iuncatngr,'ECLM','J2000','OEF2.0',dyn%ndp)
! orbit file (name_obj.ctc)
  ctcname=name_obj(1:le)//'.ctc'
  CALL filopn(iunctc,ctcname,'unknown')
  CALL wromlh(iunctc,'ECLM','J2000') 
  DO i=imim,imip 
     WRITE(*,144)i,elm(i)%coord,dpm(1:dyn%ndp,i)
     WRITE(iun,144)i,elm(i)%coord,dpm(1:dyn%ndp,i)
144  FORMAT(I5,5F16.12,F16.9,1P,4(1X,D14.7)) 
  ENDDO
! report file (with MOID, RMS, etc.) 
  repname=name_obj(1:le)//'.mrep'                        
  CALL filopn(iunrep,repname,'unknown') 
! Header
  WRITE(iunrep,120)
120 FORMAT('MULTIPLE SOLUTIONS COMPUTATION, prob_sampl=.FALSE.')
  WRITE(iunrep,220)imult,sigma 
220 FORMAT('MULTIPLE SOLUTIONS SUMMARY, imult=',i5,' maxsigma=',F5.3) 
  WRITE(iunrep,320)scaling_lov,second_lov
320 FORMAT('Scaling LOV= ',L1/'Second  LOV= ',L1)
  WRITE(*,420)
  WRITE(iun,420)
  WRITE(iunrep,420)
420 FORMAT('   no      RMS        lastcor    mag    MOID       nod+       nod-       sigma_lov      delta_sigma        sigQ')
  deltasigma=sigma/imult
  DO i=imim,imip 
     sigmalov(i)=(i-(imult+1))*deltasigma
     WRITE(*,145)i,csinom(i),delnom(i),elm(i)%h_mag,moid_m(i),dnp_m(i),dnm_m(i),sigmalov(i),deltasigma,sigq(i)
     WRITE(iun,145)i,csinom(i),delnom(i),elm(i)%h_mag,moid_m(i),dnp_m(i),dnm_m(i),sigmalov(i),deltasigma,sigq(i)
     WRITE(iunrep,145)i,csinom(i),delnom(i),elm(i)%h_mag,moid_m(i),dnp_m(i),dnm_m(i),sigmalov(i),deltasigma,sigq(i)
145  FORMAT(I5,1X,1P,E13.5,E11.3,2X,0P,F5.2,1X,F8.5,1X,F10.5,1X,F10.5,3X,F12.8,4X,F12.8,3X,F12.8)
     astna0=' '
     IF(i.le.9)THEN              
        WRITE(astna0,101) i 
101     FORMAT('va_00000',I1)
     ELSEIF(i.le.99)THEN
        WRITE(astna0,102) i
102     FORMAT('va_0000',I2)
     ELSEIF(i.le.999)THEN
        WRITE(astna0,103) i
103     FORMAT('va_000',I3)
     ELSEIF(i.le.9999)THEN
        WRITE(astna0,104) i
104     FORMAT('va_00',I4)
     ELSEIF(i.le.99999)THEN
        WRITE(astna0,105) i
105     FORMAT('va_0',I5)
     ELSEIF(i.le.999999)THEN
        WRITE(astna0,106) i
106     FORMAT('va_',I6)
     ELSE
        WRITE(*,*)' fmulti: change format for imult>499999'
        STOP
     ENDIF
     CALL rmsp(astna0,le)
     IF(dyn%ndp.gt.0)THEN 
        IF(nd.GT.6)THEN
           dyn%dp(1:dyn%ndp)=dpm(1:dyn%ndp,i)
        END IF
        CALL write_elems(elm(i),astna0(1:le),'ML',dyn,nls,ls,unm(i),UNIT=iunctc)
     ELSE
        CALL write_elems(elm(i),astna0(1:le),'ML',UNC=unm(i),UNIT=iunctc)
     ENDIF
! Writing catalog (single line, with and without non-grav perturbations)
     IF(dyn%ndp.gt.0)THEN 
        IF(nd.GT.6)THEN
           dyn%dp(1:dyn%ndp)=dpm(1:dyn%ndp,i)
           ngr_case=ngr_case+1
        END IF
        CALL write_elems(elm(i),astna0(1:le),'1L',dyn,nls,ls,FILE=' ',UNIT=iuncat,UNIT2=iuncatngr)
     ELSE
        CALL write_elems(elm(i),astna0(1:le),'1L',FILE=' ',UNIT=iuncat)
     ENDIF
  ENDDO
  CALL filclo(iunctc,' ') 
  CALL filclo(iuncat,' ')
! Delete .cat.ngr if empty
  IF(ngr_case.EQ.0)THEN
     CALL filclo(iuncatngr,'DELETE') 
  ELSE
     CALL filclo(iuncatngr,' ')
  END IF
  CALL filclo(iunrep,' ') 
END SUBROUTINE lin_multi

! =======================================
!  LIN_MULTI_GEN
! =======================================         
SUBROUTINE lin_multi_gen(batch,el0,unc0,prob_step,sigma_step_max,sigma_max,nd,imult,iwdir,m,obs,obsw,csinor)
  USE propag_state
  USE tp_trace
  USE close_app, ONLY: kill_propag
  USE semi_linear
! ================INPUT====================================================================                      
  LOGICAL,          INTENT(IN) :: batch          ! batch/interactive
  TYPE(orbit_elem), INTENT(IN) :: el0            ! initial conditions: epoch, elements, 
                                                 ! abs. magnitude, gmag 
  TYPE(orb_uncert), INTENT(IN) :: unc0           ! normal and covariance matrices
  INTEGER,          INTENT(IN) :: m              ! no. observations
  TYPE(ast_obs), DIMENSION(m)  :: obs
  TYPE(ast_wbsr),DIMENSION(m)  :: obsw 
  DOUBLE PRECISION, INTENT(IN) :: csinor         ! norms of residuals, of last correction
  DOUBLE PRECISION, INTENT(IN) :: prob_step      ! completeness limit
  DOUBLE PRECISION, INTENT(IN) :: sigma_step_max ! max step in sigma
  DOUBLE PRECISION, INTENT(IN) :: sigma_max      ! max of sigma along the line
  INTEGER,          INTENT(IN) :: nd             ! dimension of solve for parameters vector
  INTEGER,          INTENT(IN) :: iwdir          ! unit for TP scatter vector output
! ===============OUTPUT=====================================================================  
  INTEGER, INTENT(OUT) :: imult ! number of alternate solutions
! ==============END INTERFACE===============================================================
  TYPE(tp_point)   :: tptr
  TYPE(orbit_elem) :: eltmp,elf               ! temporary orbit elements, final after scatter
  TYPE(orb_uncert) :: uncf                    ! covariance after scatter           
  DOUBLE PRECISION :: wdir0(ndimx),sdir0      ! weak direction, rms along it
  DOUBLE PRECISION :: dn,delta_sigma          ! step in index, in sigma 
  DOUBLE PRECISION :: dp0(ndyx),dp(ndyx)      ! non-grav parameters
  DOUBLE PRECISION :: eigval(2),fv1(2),fv2(2) ! eigenvalues, workspace, transposed     
  INTEGER          :: ierr
  DOUBLE PRECISION :: sig(2),axes(2,2)
  DOUBLE PRECISION :: dynm(ndyx,mulx)
  DOUBLE PRECISION :: elmcoo(6,mulx)
  LOGICAl          :: bizarre                 ! control of bizarre orbits
  LOGICAL          :: falsok
  INTEGER          :: ipla
  DOUBLE PRECISION :: tafter,tbefore          ! times for scatterplane
  DOUBLE PRECISION :: v_infty,v_inf0          ! velocity at the boundary of Earth's sphere of influence
  INTEGER          :: j,i,imi,np,nm,jc
  DOUBLE PRECISION :: alpha, tprmin
  DOUBLE PRECISION :: prscal                  ! function
! for scaling
  DOUBLE PRECISION,DIMENSION(ndimx) :: units
! weak direction from scatterplane
  DOUBLE PRECISION                        :: gtp(2,2) ! covariance on target plane
  REAL(KIND=qkind)                        :: gtp_q(2,2) 
  REAL(KIND=qkind)                        :: ceicel(ndimx-2,2),v(ndimx,ndimx),b(2,2) ! for slinel
  REAL(KIND=qkind),DIMENSION(2)           :: wtp, xy, wtp1, wtp2, wtp_ell
  REAL(KIND=qkind),DIMENSION(ndimx-2)     :: hmhs
  REAL(KIND=qkind),DIMENSION(ndimx)       :: vvv1, vvv2
  REAL(KIND=qkind)                        :: dtpdetq(1:ndimx,1:2),dtpdeq(1:2,1:ndimx)
  REAL(KIND=qkind),DIMENSION(ndimx,ndimx) :: c,g
  REAL(KIND=qkind)                        :: chi1,chi2,chi12,chi,dchi
  REAL(KIND=qkind),DIMENSION(ndimx)       :: cv1,cv2,cv
! for call to fdiff_cor
  INTEGER          :: nused,itmaxold
  LOGICAL          :: succ,ok,cov0
  DOUBLE PRECISION :: delno0,rmsh0
  CHARACTER(60)    :: rwofi0
! LoV sampling
  DOUBLE PRECISION :: disk_re, p1 ! constants in the generic completeness formula
  DOUBLE PRECISION, DIMENSION(0:(mulx-1)/2) :: pdf, sigmav, deltasigv ! temporary arrays
! file units
  INTEGER :: iun,iundump,iunint,iunctc,iuncat,iunrep,iuncatngr
! file names
  CHARACTER(LEN=160) :: catname,repname,ctcname,catngrname
! non-gravitational parameters
  INTEGER :: ngr_case
! hack for catalog                                                      
  INTEGER :: le 
  CHARACTER(LEN=20) :: astna0 
! =============================================================================================
  iun=abs(iun_log)
  IF(verb_mul.ge.20)THEN 
     iunint=abs(iun_log) 
  ELSE 
     iunint=-abs(iun_log) 
  ENDIF
! computation of sampling LOV points in terms of sigma_LOV
  disk_re=0.2*aukm/eradkm
  p1=0.5d0*disk_re*prob_step
  sigmav(0)=0.d0
  pdf(0)=1.d0/sqrt(dpig)
  deltasigv(0)=p1/pdf(0)
  DO i=1,mulx/2
    sigmav(i)=sigmav(i-1)+deltasigv(i-1)
    pdf(i)=pdf(0)*exp(-0.5d0*sigmav(i)**2)
    deltasigv(i)=MIN(p1/pdf(i),sigma_step_max)
    IF(sigmav(i).gt.sigma_max)GOTO 2
  ENDDO
  WRITE(*,*)' mult_gen: error in sampling:', i, mulx, pdf(i)
  STOP  
2 CONTINUE
  imult=i ! now the arrays sigmav,deltasigv are filled from 1 to imult
! ===================================================================== 
! nominal solution at the center of the list                            
  imi0=imult+1 
  elm(imi0)=el0 
  unm(imi0)=unc0
  sigmavalv(imi0)=sigmav(0)
  deltasigmav(imi0)=deltasigv(0)
  csinom(imi0)=csinor
  DO j=1,dyn%ndp
    dpm(j,imi0)=dyn%dp(j)
  ENDDO
  ! orbital distance                                                      
  CALL nomoid(el0%t,elm(imi0),moid_m(imi0),dnp_m(imi0),dnm_m(imi0))
! ===================================================================== 
! weak direction and weak axis
  IF(scatterplane)THEN
! select propagation time                                               
     ipla=3 
     v_inf0=v_infty(el0) 
     CALL aftclo2v(ipla,el0%t,before_scatter,after_scatter,v_inf0,tbefore,tafter) 
! availability of covariance for TP analysis online                     
     CALL cov_avai(unc0,el0%coo,el0%coord) 
! reset close approach storage                                          
     njc=0 
     iplam=0 
! propagate to find linear confidence ellipse in scattering TP
     CALL pro_ele(el0,tafter,elf,unc0,uncf) 
! check success in finding the required target plane
     IF(njc.eq.0)THEN 
        WRITE(*,*)' lin_lov failed; no close approach' 
        STOP
     ELSEIF(iplam.ne.3)THEN 
        WRITE(*,*)' lin_lov failed; close app. planet ', iplam, tafter
        STOP
     ELSEIF(tp_store(njc)%tcla.gt.after_scatter.or.tp_store(njc)%tcla.lt.before_scatter)THEN   
        IF(kill_propag)THEN
           WRITE(*,*)' lin_lov failed: previous collision at ',tcla(njc)
        ELSE
           WRITE(*,*)' lin_lov failed; tcla not in interval ',tcla(njc) 
        ENDIF
        STOP
     ELSE 
        IF(kill_propag)THEN
           WRITE(*,*)' lin_lov: collision found at nominal orbit'
        ENDIF
     ENDIF
! load into tp_point data type
! select lower minimum                                                  
     IF(njc.eq.1)THEN 
        jc=1 
     ELSE 
        tprmin=1.d10 
        DO j=1,njc
           IF(tp_store(j)%d.lt.tprmin)THEN 
              jc=j 
              tprmin=tp_store(j)%d 
           ENDIF
        ENDDO
     ENDIF
! ===quadruple precision portion=====================
! conversion to quadruple precision
     dtpdetq(1:nd,1:2)=tp_store(jc)%dtpdet(1:nd,1:2)
     dtpdeq(1:2,1:nd)=TRANSPOSE(dtpdetq(1:nd,1:2))
     g(1:nd,1:nd)=unc0%g(1:nd,1:nd)
     c(1:nd,1:nd)=unc0%c(1:nd,1:nd)
! covariance on target plane
     gtp_q=MATMUL(MATMUL(dtpdeq(1:2,1:nd),g(1:nd,1:nd)),dtpdetq(1:nd,1:2))
! WARNING: still in double precision
     gtp=gtp_q
! eigenvalues, major semiaxis of ellipse on TP     
     CALL rs(2,2,gtp,eigval,1,axes,fv1,fv2,ierr)
! END WARNING
! conversion to quadruple precision
     wtp_ell=axes(:,2)*sqrt(eigval(2)) 
! find counterimage in initial conditions
     CALL slinel_q(dtpdetq(1:nd,1:2),g(1:nd,1:nd),c(1:nd,1:nd),ceicel(1:nd-2,1:2),b,v(1:nd,1:nd),nd)
!     CALL slinel(dtpdetq(1:nd,1:2),g(1:nd,1:nd),c(1:nd,1:nd),ceicel(1:nd-2,1:2),b,v(1:nd,1:nd),nd)
     xy=MATMUL(b,wtp_ell) ! g-g^*
     hmhs(1:nd-2)=MATMUL(ceicel(1:nd-2,1:2),xy)! h-h^*
     vvv1(1:nd)=MATMUL(v(1:nd,1:2),xy) ! R^T [g-g^* 0]
     wtp1=MATMUL(dtpdeq(1:2,1:nd),vvv1(1:nd))
     vvv2(1:nd)=MATMUL(v(1:nd,3:nd),hmhs(1:nd-2)) ! R^T [0 (h-h^*)]   
     wtp2=MATMUL(dtpdeq(1:2,1:nd),vvv2(1:nd))
     wtp=wtp1-wtp2
     vvv(1:nd)=vvv1(1:nd)-vvv2(1:nd) ! R^T [0 h-h^*]
! axial vector continued as true vector
     IF(el0%coo.eq.'EQU'.or.el0%coo.eq.'KEP'.or.el0%coo.eq.'COM'.or.el0%coo.eq.'COT')THEN
! if coo='EQU', 'KEP' a growing with sigma, if 'COM' q growing with sigma 
        IF(vvv(1).lt.0.d0) vvv(1:nd)=-vvv(1:nd)
     ELSEIF(el0%coo.eq.'ATT')THEN
! if coo='ATT' r growing with sigma
        IF(vvv(5).lt.0.d0) vvv(1:nd)=-vvv(1:nd)
     ELSEIF(el0%coo.eq.'CAR')THEN
! if coo='CAR' r growing with sigma
        IF(prscal(vvv(1:3),el0%coord(1:3)).lt.0.d0) vvv(1:nd)=-vvv(1:nd)
     ENDIF
! output direction of LOV in initial conditions space for use by resret
     WRITE(iwdir,100) vvv(1:nd)
 100 FORMAT(1P,10E32.22)
  ELSE
! compute line of variations                                            
     CALL weak_dir(unc0%g(1:nd,1:nd),wdir0,sdir0,iunint,el0%coo,el0%coord,units,nd)
     vvv(1:nd)=units(1:nd)*wdir0(1:nd)*sdir0
  ENDIF
! prepare for iterations
  IF(nd.gt.6)THEN
     DO j=1,nls
        dp0(j)=dyn%dp(ls(j))
     ENDDO
     dp(1:nls)=dp0(1:nls)
  ENDIF
  sigq(imi0)=0.d0
  itmaxold=itmax
! ===================================================================== 
! loop on positive sigma-parameter
  DO i=1,imult
     imi=imi0+i
     WRITE(*,*)' alternate solution no. ',imi 
     WRITE(iun,*)' alternate solution no. ',imi 
     sigmavalv(imi)=sigmav(i)
     deltasigmav(imi)=deltasigv(i)
     elm(imi)=elm(imi-1)
     elm(imi)%coord=elm(imi)%coord+vvv(1:6)*deltasigmav(imi)
     IF(nd.gt.6)THEN
        dpm(:,imi)=dpm(:,imi-1)
        DO j=1,nls
           dp(j)=dp(j)+vvv(6+j)*deltasigmav(imi)
           dpm(ls(j),imi)=dp(j)
! for possible computation of residuals
           dyn%dp(ls(j))=dpm(ls(j),imi)
        ENDDO
     ENDIF
     unm(imi)=unm(imi0)
! orbital distance                                                      
     CALL nomoid(el0%t,elm(imi),moid_m(imi),dnp_m(imi),dnm_m(imi))
! computation of residuals
     IF(MOD(imi,100).eq.1)THEN
        cov0=.false.
        itmax=0
        CALL fdiff_cor(.true.,-1,.true.,.true.,ok,cov0,elm(imi),m,obs,obsw,nused,   &
     &     rwofi0,elm(imi),unm(imi),csinom(imi),delno0,rmsh0,succ)
        sigq(imi)=sqrt(abs(csinom(imi)**2-csinom(imi0)**2)*nused)
        IF(sigq(imi-1).lt.0.d0) sigq(imi)=-sigq(imi)
     ELSEIF(MOD(imi,100).eq.2.and.imi.ne.imi0+1)THEN
         sigq(imi)=sigq(imi-2)+deltasigmav(imi-2)*2
     ELSE
        sigq(imi)=sigq(imi-1)+deltasigmav(imi-1)
     ENDIF
  ENDDO
  imip=imi0+imult
! ===================================================================== 
! loop on negative sigma-parameter
  dp(1:nls)=dp0(1:nls)
  DO i=1,imult
     imi=imi0-i
     WRITE(*,*)' alternate solution no. ',imi 
     WRITE(iun,*)' alternate solution no. ',imi 
     elm(imi)=elm(imi+1)
     sigmavalv(imi)=-sigmav(i)
     deltasigmav(imi)=deltasigv(i)
     elm(imi)%coord=elm(imi)%coord-vvv(1:6)*deltasigmav(imi)
     IF(nd.gt.6)THEN
        dpm(:,imi)=dpm(:,imi+1)
        DO j=1,nls
           dp(j)=dp(j)-vvv(6+j)*deltasigmav(imi)
           dpm(ls(j),imi)=dp(j)
! for possible computation of residuals
           dyn%dp(ls(j))=dpm(ls(j),imi)
        ENDDO
     ENDIF
     unm(imi)=unm(imi0)
! orbital distance                                                      
     CALL nomoid(el0%t,elm(imi),moid_m(imi),dnp_m(imi),dnm_m(imi))
! computation of residuals
     IF(MOD(imi,100).eq.1)THEN
        cov0=.false.
        itmax=0
        CALL fdiff_cor(.true.,-1,.true.,.true.,ok,cov0,elm(imi),m,obs,obsw,nused,   &
     &     rwofi0,elm(imi),unm(imi),csinom(imi),delno0,rmsh0,succ)
        sigq(imi)=sqrt(abs(csinom(imi)**2-csinom(imi0)**2)*nused)
        IF(sigq(imi+1).lt.0.d0) sigq(imi)=-sigq(imi)
     ELSEIF(MOD(imi,100).eq.0.and.imi.ne.imi0-1)THEN
         sigq(imi)=sigq(imi+2)-deltasigmav(imi+2)*2
     ELSE
        sigq(imi)=sigq(imi+1)-deltasigmav(imi+1)
     ENDIF
  ENDDO
  imim=imi0-imult
  itmax=itmaxold
! catalog reference time
  tdt_cat=el0%t
! =============OUTPUT, IF ANY ============================
  IF(batch) RETURN
! count cases with non-grav parameters
  ngr_case=0
! Object name
  CALL rmsp(name_obj,le)    
! Catalog (name_obj.cat)
  catname=name_obj(1:le)//'.cat'
  CALL filopn(iuncat,catname,'unknown')
  CALL wro1lh2(iuncat,'ECLM','J2000',el0%coo,'OEF2.0')
! single line format for non-grav param                     
  catngrname=name_obj(1:le)//'.cat.ngr'
  CALL filopn(iuncatngr,catngrname,'unknown') 
  CALL wro1lh2_ngr(iuncatngr,'ECLM','J2000','OEF2.0',dyn%ndp)
! orbit file (name_obj.ctc)
  ctcname=name_obj(1:le)//'.ctc'
  CALL filopn(iunctc,ctcname,'unknown') 
  CALL wromlh(iunctc,'ECLM','J2000') 
  DO i=imim,imip 
     WRITE(*,144)i,elm(i)%coord,dpm(1:dyn%ndp,i) 
     WRITE(iun,144)i,elm(i)%coord,dpm(1:dyn%ndp,i)
144  FORMAT(I5,5F16.12,F16.9,1P,4(1X,D14.7)) 
  ENDDO
! report file (with MOID, RMS, etc.) 
  repname=name_obj(1:le)//'.mrep'                        
  CALL filopn(iunrep,repname,'unknown') 
! Header
  WRITE(iunrep,120)
120 FORMAT('MULTIPLE SOLUTIONS COMPUTATION, prob_sampl=.TRUE.')
  WRITE(iunrep,220)imult,sigma_max 
220 FORMAT('MULTIPLE SOLUTIONS SUMMARY, imult=',i5,' maxsigma=',F5.3) 
  WRITE(iunrep,320)scaling_lov,second_lov
320 FORMAT('Scaling LOV= ',L1/'Second  LOV= ',L1)
  WRITE(*,420)
  WRITE(iun,420)
  WRITE(iunrep,420)
420 FORMAT('   no      RMS        lastcor    mag    MOID       nod+       nod-       sigma_lov      delta_sigma        sigQ')
  DO i=imim,imip 
     WRITE(*,145)i,csinom(i),delnom(i),elm(i)%h_mag,moid_m(i),dnp_m(i),dnm_m(i),sigmavalv(i),deltasigmav(i),sigq(i)
     WRITE(iun,145)i,csinom(i),delnom(i),elm(i)%h_mag,moid_m(i),dnp_m(i),dnm_m(i),sigmavalv(i),deltasigmav(i),sigq(i)
     WRITE(iunrep,145)i,csinom(i),delnom(i),elm(i)%h_mag,moid_m(i),dnp_m(i),dnm_m(i),sigmavalv(i),deltasigmav(i),sigq(i)
145  FORMAT(I5,1X,1P,E13.5,E11.3,2X,0P,F5.2,1X,F8.5,1X,F10.5,1X,F10.5,3X,F12.8,4X,F12.8,3X,F12.8)
     astna0=' '
     IF(i.le.9)THEN              
        WRITE(astna0,101) i 
101     FORMAT('va_00000',I1)
     ELSEIF(i.le.99)THEN
        WRITE(astna0,102) i
102     FORMAT('va_0000',I2)
     ELSEIF(i.le.999)THEN
        WRITE(astna0,103) i
103     FORMAT('va_000',I3)
     ELSEIF(i.le.9999)THEN
        WRITE(astna0,104) i
104     FORMAT('va_00',I4)
     ELSEIF(i.le.99999)THEN
        WRITE(astna0,105) i
105     FORMAT('va_0',I5)
     ELSEIF(i.le.999999)THEN
        WRITE(astna0,106) i
106     FORMAT('va_',I6)
     ELSE
        WRITE(*,*)' fmulti: change format for imult>499999'
        STOP
     ENDIF
     CALL rmsp(astna0,le)
     IF(dyn%ndp.gt.0)THEN 
        IF(nd.GT.6)THEN
           dyn%dp(1:dyn%ndp)=dpm(1:dyn%ndp,i)
        END IF
        CALL write_elems(elm(i),astna0(1:le),'ML',dyn,nls,ls,unm(i),UNIT=iunctc)
     ELSE
        CALL write_elems(elm(i),astna0(1:le),'ML',UNC=unm(i),UNIT=iunctc)
     ENDIF
! Writing catalog (single line, with and without non-grav perturbations)
     IF(dyn%ndp.gt.0)THEN 
        IF(nd.GT.6)THEN
           dyn%dp(1:dyn%ndp)=dpm(1:dyn%ndp,i)
           ngr_case=ngr_case+1
        END IF
        CALL write_elems(elm(i),astna0(1:le),'1L',dyn,nls,ls,FILE=' ',UNIT=iuncat,UNIT2=iuncatngr)
     ELSE
        CALL write_elems(elm(i),astna0(1:le),'1L',FILE=' ',UNIT=iuncat)
     ENDIF
  ENDDO
  CALL filclo(iunctc,' ') 
  CALL filclo(iuncat,' ')
! Delete .cat.ngr if empty
  IF(ngr_case.EQ.0)THEN
     CALL filclo(iuncatngr,'DELETE') 
  ELSE
     CALL filclo(iuncatngr,' ')
  END IF
  CALL filclo(iunrep,' ') 
END SUBROUTINE lin_multi_gen

! =======================================                               
!  MULTI_GEN                                                               
! =======================================                               
! interface for multiple solutions               
! version with adaptive sigma-stepsize
! ===============INTERFACE========================
SUBROUTINE multi_gen(batch,el0,uncert,csinor,delnor,mc,obs,obsw,   &
     &     prob_step,sigma_step_max,sigma_max,        &
     &     nd,imult,csino0) 
! ================INPUT=================================================
  LOGICAL,          INTENT(IN)          :: batch ! batch/interactive
  INTEGER,          INTENT(IN)          :: mc    ! number of observations
! data types
  TYPE(ast_obs),DIMENSION(mc)           :: obs
  TYPE(ast_wbsr),DIMENSION(mc)          :: obsw
! initial conditions: epoch, elements, abs. magnitude, gmag             
  TYPE(orbit_elem), INTENT(IN)          :: el0  
! normal and covariance matrices, norms of residuals, of last correction
  TYPE(orb_uncert), INTENT(IN)          :: uncert 
  DOUBLE PRECISION, INTENT(IN)          :: csinor,delnor 
! step in maximum impact probability, max step in sigma_lov, max value of sigma_lov
  DOUBLE PRECISION, INTENT(IN)          :: prob_step,sigma_step_max,sigma_max
  INTEGER,          INTENT(IN)          :: nd    ! dimension of parameters space
! value of residual norm at the central solution, if different from csinor 
  DOUBLE PRECISION, INTENT(IN),OPTIONAL :: csino0
! ================OUTPUT===============================================
  INTEGER, INTENT(OUT),OPTIONAL         :: imult
! ==============END INTERFACE==========================================
! NOTE: el0, uncert, csinor, delnor, mc, obs, obsw have to be available
! might need to be protected by a paranoia-check interface (like in f_multi)
! =====================================================================
! local variables used in computation of sigmavalv, deltasigv vectors
  DOUBLE PRECISION,DIMENSION(0:(mulx-1)/2) :: pdf,sigmav,deltasigv ! temporary arrays
  INTEGER           :: i,j,imi,k
  DOUBLE PRECISION  :: disk_re,p1 ! constats in the generic completeness formula
  DOUBLE PRECISION  :: rescov  ! function to compute sigma_q
  DOUBLE PRECISION  :: csinoref ! reference value, should be minimum
  DOUBLE PRECISION  :: rmsh ! RMS photometry
  TYPE(orbit_elem)  :: eqf,elc ! temporary orbit elements
  INTEGER           :: ngr_case ! non-grav parameters
! weak direction, rms along it    
  DOUBLE PRECISION  :: wdir(ndimx),wdir0(ndimx),sdir,dp0(ndyx),dp(ndyx), ecc
  INTEGER           :: nused ! no obs, used
! control of bizarre orbits                                             
  LOGICAL           :: bizarre,succ,fail
! file names
  CHARACTER(LEN=160) :: catname,repname,ctcname,catngrname
! units for reading
  INTEGER           :: iunint,iun
  INTEGER           :: iuncat, iuncatngr, iunrep, iunctc
! hack for catalog                                                      
  INTEGER           :: le 
  CHARACTER(LEN=20) :: astna0 
! for scaling
  DOUBLE PRECISION,DIMENSION(ndimx) :: units
! ===============BEGIN EXECUTION=======================================
  iun=abs(iun_log)
  IF(verb_mul.ge.20)THEN 
     iunint=abs(iun_log) 
  ELSE 
     iunint=-abs(iun_log) 
  ENDIF
! computation of sampling LOV points in terms of sigma_LOV
  disk_re=0.2*aukm/eradkm
  p1=0.5d0*disk_re*prob_step
  sigmav(0)=0.d0
  pdf(0)=1.d0/sqrt(dpig)
  deltasigv(0)=p1/pdf(0)
  DO i=1,mulx/2
    sigmav(i)=sigmav(i-1)+deltasigv(i-1)
    pdf(i)=pdf(0)*exp(-0.5d0*sigmav(i)**2)
    deltasigv(i)=MIN(p1/pdf(i),sigma_step_max)
    IF(sigmav(i).gt.sigma_max)GOTO 2
  ENDDO
  WRITE(*,*)' mult_gen: error in sampling:', i, mulx, pdf(i)
  STOP  
2 CONTINUE
  imult=i ! now the arrays sigmav,deltasigv are filled from 1 to imult
! ===================================================================== 
! compute line of variations                                            
  CALL weak_dir(uncert%g(1:nd,1:nd),wdir,sdir,iunint,el0%coo,el0%coord,units,nd)
  wdir0=wdir 
! store in dyn_param the values of the non-grav parameters, if any
  DO k=1,nls
    dp0(k)=dyn%dp(ls(k))
  ENDDO
! ===================================================================== 
! attempt to compute 2*imult+1 LOV points
! ===================================================================== 
! nominal solution at the center of the list                            
  imi0=imult+1 
  elm(imi0)=el0 
  unm(imi0)=uncert
  csinom(imi0)=csinor 
  IF(PRESENT(csino0))THEN
     csinoref=csino0
  ELSE
     csinoref=csinor
  ENDIF 
  delnom(imi0)=delnor 
  sigq(imi0)=sqrt(abs(csinom(imi0)**2-csinoref**2)*nused)
  sigmavalv(imi0)=sigmav(0)
  deltasigmav(imi0)=deltasigv(0)
  IF(csinom(imi0).lt.csinoref)sigq(imi0)=-sigq(imi0)
  DO j=1,dyn%ndp
    dpm(j,imi0)=dyn%dp(j)
  ENDDO
! orbital distance                                                      
  CALL nomoid(el0%t,elm(imi0),moid_m(imi0),dnp_m(imi0),dnm_m(imi0))
! ===================================================================== 
! main loop on number of steps (positive side)                          
! ===================================================================== 
  DO 5 j=1,imult 
     imi=imult+j+1 
     sigmavalv(imi)=sigmav(j)
     deltasigmav(imi)=deltasigv(j)
     WRITE(*,*)' alternate solution no. ',imi 
     WRITE(iun,*)' alternate solution no. ',imi
     CALL prop_sig(.true.,elm(imi-1),elc,deltasigmav(imi),1.d0,   &
     &             mc,obs,obsw,wdir,sdir,units,fail,nd)
! check for hyperbolic                                                  
     IF(fail)THEN 
        WRITE(*,*)'prop_sig failed, imi= ',imi,sigmavalv(imi),deltasigmav(imi)
        imi=imi-1
        GOTO 6 
     ELSEIF(bizarre(elc,ecc))THEN 
        WRITE(*,*)'bizarre out of  prop_sig, imi= ',imi,ecc,sigmavalv(imi),deltasigmav(imi)
        imi=imi-1 
        GOTO 6 
     ELSE 
! constrained differential corrections:
         CALL constr_fit(mc,obs,obsw,elc,wdir,elm(imi),unm(imi),     &
     &       csinom(imi),delnom(imi),rmsh,nused,succ,nd)
! exit if not convergent                                                
        IF(.not.succ) THEN 
           WRITE(*,*)'fail in constr_fit, imi= ',imi,sigmavalv(imi),deltasigmav(imi)
           imi=imi-1 
           GOTO 6 
        ENDIF
! array of non-grav parameters
        DO k=1,dyn%ndp
           dpm(k,imi)=dyn%dp(k)
        ENDDO
! orbital distance                                                      
        CALL nomoid(elm(imi)%t,elm(imi),moid_m(imi),dnp_m(imi),dnm_m(imi))
! check for sigQ.le.sigma_max                                               
        sigq(imi)=sqrt(abs(csinom(imi)**2-csinoref**2)*nused)
! this has problems: if negative $\chi^2$ occurs, the values are scrambled.
        IF(csinom(imi).lt.csinoref)sigq(imi)=-sigq(imi)
        IF(sigq(imi).gt.sigma_max)THEN 
           WRITE(iun_log,*)'too large sigma_q ',imi,sigq(imi),sigmavalv(imi)
           GOTO 6 
        ENDIF
     ENDIF
5 END DO
6 imip=imi
  WRITE(*,*) ' end of forward LOV sampling, imip=',imip
! ===================================================================== 
! line of variations, negative sigma
  wdir=wdir0
  DO k=1,nls
    dyn%dp(ls(k))=dp0(k) ! restart from nominal values, if any non-grav
  ENDDO
  DO 7 j=1,imult 
! ===================================================================== 
! main loop on number of steps (negative side)                          
! ===================================================================== 
     imi=imult+1-j  
     sigmavalv(imi)=-sigmav(j)
     deltasigmav(imi)=deltasigv(j)
     WRITE(*,*)' alternate solution no. ',imi   
     WRITE(iun,*)' alternate solution no. ',imi
     CALL prop_sig(.true.,elm(imi+1),elc,deltasigmav(imi),-1.d0,mc,obs,obsw,wdir,sdir,units,fail,nd) 
     IF(fail)THEN ! check for divergence
        WRITE(*,*)'prop_sig fail, imi= ',imi,sigmavalv(imi),deltasigmav(imi)
        imi=imi+1 
        GOTO 8 
     ELSEIF(bizarre(elc,ecc))THEN ! check for hyperbolic
        WRITE(*,*)'step ',imi,' bizarre', ecc, sigmavalv(imi),deltasigmav(imi) 
        imi=imi+1 
        GOTO 8 
     ELSE 
! differential corrections:                                             
! constrained differential corrections:
        CALL constr_fit(mc,obs,obsw,elc,wdir,elm(imi),unm(imi),     &
     &       csinom(imi),delnom(imi),rmsh,nused,succ,nd)
        WRITE(*,*) imi           
! exit if not convergent 
        IF(.not.succ)THEN 
           WRITE(*,*)'constr_fit fail,imi= ',imi,sigmavalv(imi),deltasigmav(imi)
           imi=imi+1 
           GOTO 8 
        ENDIF
! array of non-grav parameters
        DO k=1,dyn%ndp
           dpm(k,imi)=dyn%dp(k)
        ENDDO
! orbital distance                                                      
        CALL nomoid(elm(imi)%t,elm(imi),moid_m(imi),dnp_m(imi),dnm_m(imi))
        sigq(imi)=sqrt(abs(csinom(imi)**2-csinoref**2)*nused)
        IF(csinom(imi).lt.csinoref)sigq(imi)=-sigq(imi)
        IF(sigq(imi).gt.sigma_max)THEN 
           WRITE(iun_log,*)'too large sigma_q ',imi,sigq(imi),sigmavalv(imi)
           GOTO 8 
        ENDIF
     ENDIF
7 ENDDO
! reset non-grav parameters to nominal
  DO k=1,nls
    dyn%dp(ls(k))=dp0(k)
  ENDDO
! =============OUTPUT, IF ANY ============================
! the mrep file needs to contain sigmavalv on top of sigq
8 imim=imi 
! catalog reference time
  tdt_cat=el0%t
  WRITE(*,*) ' MULTI_GEN: end of backward LOV sampling, imip=',imim
! Batch
  IF(batch) RETURN
! count cases with non-grav parameters
  ngr_case=0
! Catalog (name_obj.cat)
  CALL rmsp(name_obj,le)                       
  catname=name_obj(1:le)//'.cat'
  CALL filopn(iuncat,catname,'unknown')
  CALL wro1lh2(iuncat,'ECLM','J2000',el0%coo,'OEF2.0')
! single line format for non-grav param
  catngrname=name_obj(1:le)//'.cat.ngr'
  CALL filopn(iuncatngr,catngrname,'unknown') 
  CALL wro1lh2_ngr(iuncatngr,'ECLM','J2000','OEF2.0',dyn%ndp)
! orbit file (name_obj.ctc)
  ctcname=name_obj(1:le)//'.ctc'
  CALL filopn(iunctc,ctcname,'unknown') 
  CALL wromlh(iunctc,'ECLM','J2000') 
  DO i=imim,imip 
     WRITE(*,144)i,elm(i)%coord,dpm(1:dyn%ndp,i) 
     WRITE(iun,144)i,elm(i)%coord,dpm(1:dyn%ndp,i)
144  FORMAT(I5,5F16.12,F16.9,1P,4(1X,D14.7)) 
  ENDDO
! report file (with MOID, RMS, etc.)     
  repname=name_obj(1:le)//'.mrep'                        
  CALL filopn(iunrep,repname,'unknown') 
! Header
  WRITE(iunrep,120)
120 FORMAT('MULTIPLE SOLUTIONS COMPUTATION, prob_sampl=.TRUE.')
  WRITE(iunrep,199) prob_step,sigma_max,sigma_step_max, imult
199 FORMAT('MULTIPLE SOLUTIONS SUMMARY, prob_step=',1P,E7.1,0P,' sigma_max=',F5.3, ' sigma_step_max=',F5.3,' imult=',I5) 
  WRITE(iunrep,320)scaling_lov,second_lov
320 FORMAT('Scaling LOV= ',L1/'Second  LOV= ',L1)
  WRITE(*,420)
  WRITE(iun,420)
  WRITE(iunrep,420)
420 FORMAT('   no      RMS        lastcor    mag    MOID       nod+       nod-       sigma_lov      delta_sigma        sigQ')
  DO i=imim,imip 
     WRITE(*,145)i,csinom(i),delnom(i),elm(i)%h_mag,moid_m(i),dnp_m(i),dnm_m(i),sigmavalv(i),deltasigmav(i),sigq(i)
     WRITE(iun,145)i,csinom(i),delnom(i),elm(i)%h_mag,moid_m(i),dnp_m(i),dnm_m(i),sigmavalv(i),deltasigmav(i),sigq(i)
     WRITE(iunrep,145)i,csinom(i),delnom(i),elm(i)%h_mag,moid_m(i),dnp_m(i),dnm_m(i),sigmavalv(i),deltasigmav(i),sigq(i)
145  FORMAT(I5,1X,1P,E13.5,E11.3,2X,0P,F5.2,1X,F8.5,1X,F10.5,1X,F10.5,3X,F12.8,4X,F12.8,3X,F12.8)
     astna0=' '
     IF(i.le.9)THEN              
        WRITE(astna0,101) i 
101     FORMAT('va_00000',I1)
     ELSEIF(i.le.99)THEN
        WRITE(astna0,102) i
102     FORMAT('va_0000',I2)
     ELSEIF(i.le.999)THEN
        WRITE(astna0,103) i
103     FORMAT('va_000',I3)
     ELSEIF(i.le.9999)THEN
        WRITE(astna0,104) i
104     FORMAT('va_00',I4)
     ELSEIF(i.le.99999)THEN
        WRITE(astna0,105) i
105     FORMAT('va_0',I5)
     ELSEIF(i.le.999999)THEN
        WRITE(astna0,106) i
106     FORMAT('va_',I6)
     ELSE
        WRITE(*,*)' fmulti: change format for imult>499999'
        STOP
     ENDIF
     CALL rmsp(astna0,le)
     IF(dyn%ndp.gt.0)THEN 
        IF(nd.GT.6)THEN
           dyn%dp(1:dyn%ndp)=dpm(1:dyn%ndp,i)
        END IF
        CALL write_elems(elm(i),astna0(1:le),'ML',dyn,nls,ls,unm(i),UNIT=iunctc)
     ELSE
        CALL write_elems(elm(i),astna0(1:le),'ML',UNC=unm(i),UNIT=iunctc)
     ENDIF
! Writing catalog (single line, with and without non-grav perturbations)
     IF(dyn%ndp.gt.0)THEN 
        IF(nd.GT.6)THEN
           dyn%dp(1:dyn%ndp)=dpm(1:dyn%ndp,i)
           ngr_case=ngr_case+1
        END IF
        CALL write_elems(elm(i),astna0(1:le),'1L',dyn,nls,ls,FILE=' ',UNIT=iuncat,UNIT2=iuncatngr)
     ELSE
        CALL write_elems(elm(i),astna0(1:le),'1L',FILE=' ',UNIT=iuncat)
     ENDIF
  ENDDO
  CALL filclo(iunctc,' ') 
  CALL filclo(iuncat,' ')
! Delete .cat.ngr if empty
  IF(ngr_case.EQ.0)THEN
     CALL filclo(iuncatngr,'DELETE') 
  ELSE
     CALL filclo(iuncatngr,' ')
  END IF
  CALL filclo(iunrep,' ') 
END SUBROUTINE multi_gen
! =======================================                               
!  F_MULTI                                                               
! =======================================                               
! interface for multiple solutions               
! version with fixed sigma-stepsize                                        
! ===============INTERFACE=============================================================
SUBROUTINE f_multi(batch,obsc,inic,ok,covc,         &
     &     el0,uncert,csinor,delnor,                  &
     &     mc, obs, obsw,sigma,imult,nd,sig1,sig2,csino0) 
! ================INPUT================================================================
  LOGICAL, INTENT(IN)                   :: batch          ! batch/interactive
  INTEGER, INTENT(IN)                   :: mc             ! number of observations
! new data types
  TYPE(ast_obs),DIMENSION(mc)           :: obs 
  TYPE(ast_wbsr),DIMENSION(mc)          :: obsw
! initial conditions: epoch, elements, abs. magnitude, gmag             
  TYPE(orbit_elem), INTENT(IN)          :: el0 
! normal and covariance matrices, norms of residuals, of last correction
  TYPE(orb_uncert), INTENT(IN)          :: uncert  
  DOUBLE PRECISION, INTENT(IN)          :: csinor,delnor 
  LOGICAL, INTENT(IN)                   :: obsc,inic,covc ! logical flags  
! ========INPUT, but can also be assigned inside =======================================
  INTEGER, INTENT(IN)                   :: imult  ! number of alternate solutions
  INTEGER, INTENT(IN)                   :: nd     ! dimension of the parameters space
  DOUBLE PRECISION, INTENT(IN)          :: sigma  ! max of sigma along the line     
  DOUBLE PRECISION, INTENT(IN),OPTIONAL :: csino0 ! value of residual norm at the central solution,
                                                  ! if different from csinor 
! ===============OUTPUT=================================================================        
  LOGICAL, INTENT(OUT)                   :: ok        ! success flag     
  DOUBLE PRECISION, INTENT(OUT),OPTIONAL :: sig1,sig2 ! extreme values of sigma_Q
                                                       ! at imim, imip
! ==============END INTERFACE===========================================================
  DOUBLE PRECISION   :: rescov                  ! function to compute sigma_q
  DOUBLE PRECISION   :: csinoref                ! reference value, should be minimum
  TYPE(orbit_elem)   :: eltmp                   ! temporary orbit elements
  INTEGER            :: nused                   ! no obs, used 
! file names
  CHARACTER(LEN=160) :: catname,repname,ctcname,catngrname
! weak direction, rms along it
  DOUBLE PRECISION   :: wdir(ndimx),wdir0(ndimx),sdir,dp0(ndyx),dp(ndyx)
! ======== differential correction flags and controls ==================================
  DOUBLE PRECISION                   :: dn,sigmam    ! scalar temporaries       
  INTEGER                            :: i, j, imi,jj ! loop indexes 
  LOGICAL                            :: succ         ! success flag for diff. correction
  DOUBLE PRECISION                   :: rmsh         ! rms magnitudes  
  LOGICAL                            :: bizarre      ! control of bizarre orbits
  INTEGER                            :: ngr_case     ! non-grav parameters
! hack for catalog                                                      
  INTEGER                            :: iunctc,le 
  CHARACTER*20                       :: astna0 
! control of integration method                                         
  INTEGER                            :: imint,k 
  DOUBLE PRECISION                   :: hh,ratio,ecc
  TYPE(orbit_elem)                   :: eqf, elc 
  LOGICAL                            :: fail 
! file units
  INTEGER                            :: iun,iundump,iunint,iunrep,iuncat,iuncatngr
! for scaling
  DOUBLE PRECISION, DIMENSION(ndimx) :: units
! sigma lov and deltasigma (constant step in sigma)
  DOUBLE PRECISION                   :: deltasigma
  DOUBLE PRECISION, DIMENSION(mulx)  :: sigmalov

  iun=abs(iun_log) 
  IF(verb_mul.ge.30)THEN 
     iundump=abs(iun_log) 
  ELSE 
     iundump=-abs(iun_log) 
  ENDIF
  IF(verb_mul.ge.20)THEN 
     iunint=abs(iun_log) 
  ELSE 
     iunint=-abs(iun_log) 
  ENDIF
! ===================================================================== 
! check availability of observations and initial condition              
  IF(.not.obsc)THEN 
     WRITE(*,*)'f_multi: NO OBSERVATIONS' 
     ok=.false. 
     RETURN 
  ENDIF
  CALL chereq(2,inic,covc,el0%t,iunint,ok) 
  IF(.not.ok)THEN 
     WRITE(*,*)' f_multi: no initial data ',inic,covc,el0%t 
     RETURN 
  ENDIF
! ===================================================================== 
! check availability of JPL ephemerides and ET-UT table                 
  CALL chetim(obs(1)%time_tdt,obs(mc)%time_tdt,ok) 
  IF(.not.ok) THEN 
     WRITE(*,*)' f_multi: no JPL ephem/ET-UT table ',obs(1)%time_tdt,obs(mc)%time_tdt
     RETURN 
  ENDIF
! ===================================================================== 
! compute line of variations                                            
  CALL weak_dir(uncert%g(1:nd,1:nd),wdir,sdir,iunint,el0%coo,el0%coord,units,nd)
  wdir0=wdir 
  DO j=1,nls
    dp0(j)=dyn%dp(ls(j))
  ENDDO
! input parameters of segment on the variations line                    
  IF(batch)THEN 
! sigma and imult have to be passed in the call                         
     IF(2*imult+1.gt.mulx)THEN 
        WRITE(*,*)' too many; max is ',mulx 
        STOP ' fmulti: too many mutsol' 
     ENDIF
  ENDIF
! ===================================================================== 
! nominal stepsize (in the sigma space)                                 
  dn=1.d0/float(imult) 
! store information on stepsize
  delta_sigma=dn*sigma
! ===================================================================== 
! nominal solution at the center of the list                            
  imi0=imult+1 
  elm(imi0)=el0 
  unm(imi0)=uncert
  csinom(imi0)=csinor 
  IF(PRESENT(csino0))THEN
     csinoref=csino0
  ELSE
     csinoref=csinor
  ENDIF 
  delnom(imi0)=delnor 
  sigq(imi0)=sqrt(abs(csinom(imi0)**2-csinoref**2)*nused)
  IF(csinom(imi0).lt.csinoref)sigq(imi0)=-sigq(imi0)
  DO j=1,dyn%ndp
    dpm(j,imi0)=dyn%dp(j)
  ENDDO
  dp=0.d0
! orbital distance                                                      
  CALL nomoid(el0%t,elm(imi0),moid_m(imi0),dnp_m(imi0),dnm_m(imi0))
! ===================================================================== 
! main loop on number of steps (positive side)                          
! ===================================================================== 
  DO 5 i=1,imult 
     imi=imult+i+1 
     WRITE(*,*)' alternate solution no. ',imi 
     WRITE(iun,*)' alternate solution no. ',imi 
     CALL prop_sig(batch,elm(imi-1),elc,dn,sigma,mc,obs,obsw,wdir,sdir,units,fail,nd)
! check for hyperbolic                                                  
     IF(fail)THEN 
        IF(.not.batch)WRITE(*,*)'step ',imi,' failed' 
!        IF(.not.batch)WRITE(*,*)elc 
        WRITE(iun,*)'fail in prop_sig, imi= ',imi,dn
        imi=imi-1
        GOTO 6 
     ELSEIF(bizarre(elc,ecc))THEN 
        IF(.not.batch)WRITE(*,*)'step ',imi,' bizarre', 'ecc=', ecc 
        WRITE(iun,*)'bizarre out of  prop_sig, imi= ',imi,ecc,dn
        imi=imi-1 
        GOTO 6 
     ELSE 
! constrained differential corrections:
         CALL constr_fit(mc,obs,obsw,elc,wdir,elm(imi),unm(imi),     &
     &       csinom(imi),delnom(imi),rmsh,nused,succ,nd)
! exit if not convergent                                                
        IF(.not.succ) THEN 
           WRITE(iun,*)'fail in constr_fit, imi= ',imi,dn
           imi=imi-1 
           GOTO 6 
        ENDIF
! array of non-grav parameters
        DO j=1,dyn%ndp
           dpm(j,imi)=dyn%dp(j)
        ENDDO
! orbital distance                                                      
        CALL nomoid(elm(imi)%t,elm(imi),moid_m(imi),dnp_m(imi),dnm_m(imi))
! check for sigQ.le.sigma                                               
        sigq(imi)=sqrt(abs(csinom(imi)**2-csinoref**2)*nused)
        IF(csinom(imi).lt.csinoref)sigq(imi)=-sigq(imi)
        IF(sigq(imi).gt.sigma)THEN 
           WRITE(iun,*)'too large sigma_q ',imi,sigq(imi)
           GOTO 6 
        ENDIF
     ENDIF
5 END DO
! ===================================================================== 
! line of variations, negative sigma
6 imip=imi 
  wdir=wdir0
  DO j=1,nls
    dyn%dp(ls(j))=dp0(j)
  ENDDO
  DO 7 i=1,imult 
! ===================================================================== 
! main loop on number of steps (negative side)                          
! ===================================================================== 
     imi=imult-i+1 
     WRITE(*,*)' alternate solution no. ',imi 
     WRITE(iun,*)' alternate solution no. ',imi 
     sigmam=-sigma 
     CALL prop_sig(batch,elm(imi+1),elc,dn,sigmam,mc,obs,obsw,wdir,sdir,units,fail,nd) 
     IF(fail)THEN ! check for divergence
        IF(.not.batch)WRITE(*,*)'step ',imi,' failed' 
        WRITE(iun,*)'fail in prop_sig, imi= ',imi,dn
        imi=imi+1 
        GOTO 8 
     ELSEIF(bizarre(elc,ecc))THEN ! check for hyperbolic
        IF(.not.batch)WRITE(*,*)'step ',imi,' bizarre' 
        IF(.not.batch)WRITE(*,*)elc
        WRITE(iun,*)'bizarre out of prop_sig, imi= ',imi,ecc,dn
        imi=imi+1 
        GOTO 8 
     ELSE 
! differential corrections:                                             
! constrained differential corrections:
        CALL constr_fit(mc,obs,obsw,elc,wdir,elm(imi),unm(imi),     &
     &       csinom(imi),delnom(imi),rmsh,nused,succ,nd)           
! exit if not convergent                                                
        IF(.not.succ)THEN 
           WRITE(iun,*)'fail in constr_fit, imi= ',imi,dn
           imi=imi+1 
           GOTO 8 
        ENDIF
! array of non-grav parameters
        DO j=1,dyn%ndp
           dpm(j,imi)=dyn%dp(j)
        ENDDO
! orbital distance                                                      
        CALL nomoid(elm(imi)%t,elm(imi),moid_m(imi),dnp_m(imi),dnm_m(imi))
        sigq(imi)=sqrt(abs(csinom(imi)**2-csinoref**2)*nused)
        IF(csinom(imi).lt.csinoref)sigq(imi)=-sigq(imi)
        IF(sigq(imi).gt.sigma)THEN 
           WRITE(iun,*)'too large sigma_q ',imi,sigq(imi)
           GOTO 8 
        ENDIF
     ENDIF
7 ENDDO
! reset non-grav parameters to nominal
  DO j=1,nls
    dyn%dp(ls(j))=dp0(j)
  ENDDO
! ===================================================================== 
! summary table                                                         
! ===================================================================== 
8 imim=imi 
! extreme values of sigma_Q
  IF(PRESENT(sig1))sig1=sigq(imim)
  IF(PRESENT(sig2))sig2=sigq(imip)
! catalog reference time
  tdt_cat=el0%t
  IF(batch)RETURN
! count cases with non-grav parameters
  ngr_case=0
! Object name
  CALL rmsp(name_obj,le)    
! Catalog (name_obj.cat)
  catname=name_obj(1:le)//'.cat'
  CALL filopn(iuncat,catname,'unknown')
  CALL wro1lh2(iuncat,'ECLM','J2000',el0%coo,'OEF2.0')
! single line format for non-grav param                     
  catngrname=name_obj(1:le)//'.cat.ngr'
  CALL filopn(iuncatngr,catngrname,'unknown') 
  CALL wro1lh2_ngr(iuncatngr,'ECLM','J2000','OEF2.0',dyn%ndp)
! orbit file (name_obj.ctc)
  ctcname=name_obj(1:le)//'.ctc'
  CALL filopn(iunctc,ctcname,'unknown') 
  CALL wromlh(iunctc,'ECLM','J2000') 
  DO i=imim,imip 
     WRITE(*,144)i,elm(i)%coord,dpm(1:dyn%ndp,i) 
     WRITE(iun,144)i,elm(i)%coord,dpm(1:dyn%ndp,i)
144  FORMAT(I5,5F16.12,F16.9,1P,4(1X,D14.7)) 
  ENDDO
! report file (with MOID, RMS, etc.) 
  repname=name_obj(1:le)//'.mrep'                        
  CALL filopn(iunrep,repname,'unknown') 
! Header
  WRITE(iunrep,120)
120 FORMAT('MULTIPLE SOLUTIONS COMPUTATION, prob_sampl=.FALSE.')
  WRITE(iunrep,220)imult,sigma 
220 FORMAT('MULTIPLE SOLUTIONS SUMMARY, imult=',i5,' maxsigma=',F5.3) 
  WRITE(iunrep,320)scaling_lov,second_lov
320 FORMAT('Scaling LOV= ',L1/'Second  LOV= ',L1)
  WRITE(*,420)
  WRITE(iun,420)
  WRITE(iunrep,420)
420 FORMAT('   no      RMS        lastcor    mag    MOID       nod+       nod-       sigma_lov      delta_sigma        sigQ')
  deltasigma=sigma/imult
  DO i=imim,imip 
     sigmalov(i)=(i-(imult+1))*deltasigma
     WRITE(*,145)i,csinom(i),delnom(i),elm(i)%h_mag,moid_m(i),dnp_m(i),dnm_m(i),sigmalov(i),deltasigma,sigq(i)
     WRITE(iun,145)i,csinom(i),delnom(i),elm(i)%h_mag,moid_m(i),dnp_m(i),dnm_m(i),sigmalov(i),deltasigma,sigq(i)
     WRITE(iunrep,145)i,csinom(i),delnom(i),elm(i)%h_mag,moid_m(i),dnp_m(i),dnm_m(i),sigmalov(i),deltasigma,sigq(i)
145  FORMAT(I5,1X,1P,E13.5,E11.3,2X,0P,F5.2,1X,F8.5,1X,F10.5,1X,F10.5,3X,F12.8,4X,F12.8,3X,F12.8)
     astna0=' '
     IF(i.le.9)THEN              
        WRITE(astna0,101) i 
101     FORMAT('va_00000',I1)
     ELSEIF(i.le.99)THEN
        WRITE(astna0,102) i
102     FORMAT('va_0000',I2)
     ELSEIF(i.le.999)THEN
        WRITE(astna0,103) i
103     FORMAT('va_000',I3)
     ELSEIF(i.le.9999)THEN
        WRITE(astna0,104) i
104     FORMAT('va_00',I4)
     ELSEIF(i.le.99999)THEN
        WRITE(astna0,105) i
105     FORMAT('va_0',I5)
     ELSEIF(i.le.999999)THEN
        WRITE(astna0,106) i
106     FORMAT('va_',I6)
     ELSE
        WRITE(*,*)' fmulti: change format for imult>499999'
        STOP
     ENDIF
     CALL rmsp(astna0,le)
     IF(dyn%ndp.gt.0)THEN 
        IF(nd.GT.6)THEN
           dyn%dp(1:dyn%ndp)=dpm(1:dyn%ndp,i)
        END IF
        CALL write_elems(elm(i),astna0(1:le),'ML',dyn,nls,ls,unm(i),UNIT=iunctc)
     ELSE
        CALL write_elems(elm(i),astna0(1:le),'ML',UNC=unm(i),UNIT=iunctc)
     ENDIF
! Writing catalog (single line, with and without non-grav perturbations)
     IF(dyn%ndp.gt.0)THEN 
        IF(nd.GT.6)THEN
           dyn%dp(1:dyn%ndp)=dpm(1:dyn%ndp,i)
           ngr_case=ngr_case+1
        END IF
        CALL write_elems(elm(i),astna0(1:le),'1L',dyn,nls,ls,FILE=' ',UNIT=iuncat,UNIT2=iuncatngr)
     ELSE
        CALL write_elems(elm(i),astna0(1:le),'1L',FILE=' ',UNIT=iuncat)
     ENDIF
  ENDDO
  CALL filclo(iunctc,' ') 
  CALL filclo(iuncat,' ')
! Delete .cat.ngr if empty
  IF(ngr_case.EQ.0)THEN
     CALL filclo(iuncatngr,'DELETE') 
  ELSE
     CALL filclo(iuncatngr,' ')
  END IF
  CALL filclo(iunrep,' ') 
END  SUBROUTINE f_multi

     
!============================================================!                  
! NMULTI                                                     !
!============================================================!
! Multiple solutions for a differential correction problem   !
!============================================================!
SUBROUTINE nmulti(nam0,elc,uncert,csinor,delnor,            &
     &     mc,obs,obsw,imult,sigma,moid_min,nd) 
  USE name_rules, ONLY: name_len
  !=========================================================================================================
  INTEGER,                             INTENT(IN)    :: mc       ! Number of observations
  TYPE(ast_obs),        DIMENSION(mc), INTENT(IN)    :: obs      ! Observations
  TYPE(ast_wbsr),       DIMENSION(mc), INTENT(IN)    :: obsw     ! Weights/biases
  CHARACTER*(name_len),                INTENT(IN)    :: nam0     ! Asteroid name
  TYPE(orbit_elem),                    INTENT(IN)    :: elc      ! Orbital elements
  TYPE(orb_uncert),                    INTENT(IN)    :: uncert   ! Normal and covariance matrices
  DOUBLE PRECISION,                    INTENT(IN)    :: csinor   ! Norms of residuals
  DOUBLE PRECISION,                    INTENT(IN)    :: delnor   ! Norm of last correction
  INTEGER,                             INTENT(IN)    :: imult    ! The number of multiple sol. is 2*imult+1 
  DOUBLE PRECISION,                    INTENT(INOUT) :: sigma    ! Interval on LOV from -sigma to +sigma
  DOUBLE PRECISION,                    INTENT(OUT)   :: moid_min ! MOID
  INTEGER,                             INTENT(IN)    :: nd       ! Dimension of solve for parameters vector
  !=========================================================================================================
  INTEGER                           :: imult1                  ! To avoid intent problems             
  DOUBLE PRECISION                  :: rescov                  ! Renormalization factor       
  DOUBLE PRECISION                  :: dp(ndyx)                ! Local copy of dynamical parameters
  LOGICAL                           :: obsc                    ! Avaliability of observations
  LOGICAL                           :: inic                    ! Availability of initial conditions
  LOGICAL                           :: covc                    ! Availability of covariance
  LOGICAL                           :: ok                      ! Correct data
  INTEGER                           :: i                       ! Loop index
  LOGICAL                           :: batch                   ! Batch control
  LOGICAL                           :: linlov                  ! Computation of the linear LOV
  DOUBLE PRECISION                  :: deltasigma              ! Step in sigma
  DOUBLE PRECISION, DIMENSION(mulx) :: sigmalov                ! Sigma of the VAs
  !==== Output =============================================================================================
  CHARACTER(LEN=60)                 :: file                    ! Output file name
  INTEGER                           :: iunctc                  ! File .ctc (multiline catalog)
  INTEGER                           :: iunrep                  ! File .mrep (report file for multiple solutions)
  INTEGER                           :: iuncat                  ! File .cat (single line kep catalog)
  INTEGER                           :: iuncatngr               ! File .cat.ngr (single line catalog for non-grav par)
  INTEGER                           :: iwdir                   ! Multiple solution directory
  CHARACTER(LEN=20)                 :: astna0                  ! Catalog object name (with underscore)
  INTEGER                           :: le                      ! Length of astna0
  TYPE(orbit_elem)                  :: elk                     ! Keplerian elements for catalog (cat)
  INTEGER                           :: fail_flag               ! Flag for failure in coo_cha
  CHARACTER(LEN=6)                  :: vers                    ! OEF version
  INTEGER                           :: ngr_case                ! Counter for non-gravitational parameters cases
  !===== Magnitude correction ==============================================================================
  DOUBLE PRECISION                  :: vsize                   ! Norm of a vector (function)
  DOUBLE PRECISION                  :: d_ear_nom               ! Geocentric distance of the nominal
  DOUBLE PRECISION                  :: d_sun_nom               ! Heliocentric distance of the nominal
  DOUBLE PRECISION                  :: d_ear_VA                ! Geocentric distance of the VA
  DOUBLE PRECISION                  :: d_sun_VA                ! Heliocentric distance of the VA
  DOUBLE PRECISION                  :: phi1, phi2              ! Computation of the phase angle correction
  DOUBLE PRECISION                  :: f_alpha_nom, f_alpha_VA ! Phase angle correction
  DOUBLE PRECISION                  :: xear(6)                 ! State vector of the Earth
  DOUBLE PRECISION                  :: el_geo(3)               ! Geocentric elements of the nominal
  DOUBLE PRECISION                  :: el_geo_VA(3)            ! Geocentric elements of the VA
  DOUBLE PRECISION                  :: pha_nom, pha_VA         ! Phase angles (nominal and VAs)
  !=========================================================================================================
  !*******************!
  !  Inizializations  !
  !*******************!
  imult1=imult
  ! Set flags to true to avoid checks                                     
  obsc=.TRUE. 
  inic=.TRUE. 
  covc=.TRUE. 
  ! Batch mode                                                            
  batch=.TRUE. 
  !*******************************!
  !  Selection of LOV definition  !
  !*******************************!
  IF(lin_lov)THEN
     ! Linear LOV
     IF(scatterplane)THEN
        CALL filnam('multsol',nam0,'wdir',file,le) 
        CALL filopn(iwdir,file(1:le),'unknown') 
     ENDIF
     CALL lin_multi(batch,elc,uncert,sigma,imult,nd,iwdir,mc,obs,obsw,csinor)
     IF(scatterplane) CALL filclo(iwdir,' ')
  ELSE
     ! Constant step in sigma: call to f_multi (distribution version, but with batch mode)
     CALL f_multi(batch,obsc,inic,ok,covc,elc,uncert,csinor,delnor,mc,obs,obsw,sigma,imult1,nd)
  ENDIF
                     
  !############!
  !   OUTPUT   !
  !############!
  ! Count cases with non-grav parameters
  ngr_case=0                                           
  ! Multiple solutions catalog, multiline format                          
  CALL filnam('multsol',nam0,'ctc',file,le) 
  CALL filopn(iunctc,file(1:le),'unknown') 
  vers='OEF2.0'
  CALL wromlh2(iunctc,'ECLM','J2000',vers) 
  ! Single line, format KEP                                                    
  CALL filnam('multsol',nam0,'cat',file,le) 
  CALL filopn(iuncat,file(1:le),'unknown') 
  CALL wro1lh2(iuncat,'ECLM','J2000','KEP',vers) 
  ! Single line, format for non-grav param
  CALL filnam('multsol',nam0,'cat.ngr',file,le) 
  CALL filopn(iuncatngr,file(1:le),'unknown') 
  CALL wro1lh2_ngr(iuncatngr,'ECLM','J2000',vers,dyn%ndp) 
  ! Report file (with MOID, RMS, etc.)                                    
  CALL filnam('multsol',nam0,'mrep',file,le) 
  CALL filopn(iunrep,file(1:le),'unknown') 
  ! Information to be passed to resret2
  WRITE(iunrep,399)
399 FORMAT('MULTIPLE SOLUTIONS COMPUTATION, prob_sampl=.FALSE.')
  WRITE(iunrep,199)imult,sigma 
199 FORMAT('MULTIPLE SOLUTIONS SUMMARY, imult=',i5,' maxsigma=',F5.3) 
  WRITE(iunrep,299)scaling_lov,second_lov
299 FORMAT('Scaling LOV= ',L1/'Second  LOV= ',L1)
  WRITE(iunrep,499)
499 FORMAT('   no      RMS        lastcor    mag    MOID       nod+       nod-       sigma_lov      delta_sigma        sigQ')
  
  !#####################################################!
  !   LOOP ON MULTIPLE SOLUTIONS EFFECTIVELY COMPUTED   !
  !#####################################################!
  moid_min=1.d3          ! MOID inizialization
  deltasigma=sigma/imult ! Computation of the constant step
  !***************************************!
  !  Parameters for magnitude correction  !
  !***************************************!
  d_sun_nom=vsize(elm(imi0)%coord(1:3))  ! Heliocentric distance of the nominal
  CALL earcar(elm(imi0)%t,xear,1)
  el_geo=elm(imi0)%coord(1:3)-xear(1:3)  ! Geocentric elements of the nominal solution
  d_ear_nom=vsize(el_geo)                ! Geocentric distance of the nominal solution
  pha_nom=ACOS(DOT_PRODUCT(elm(imi0)%coord(1:3),el_geo)/(d_sun_nom*d_ear_nom)) ! Phase angle
  ! Correction depending upon the phase angle
  phi1=EXP(-3.33d0*TAN(pha_nom/2.d0)**0.63d0) 
  phi2=EXP(-1.87d0*TAN(pha_nom/2.d0)**1.22d0)
  f_alpha_nom=2.5d0*LOG10((1.d0-elm(imi0)%g_mag)*phi1+elm(imi0)%g_mag*phi2) 
  
  DO i=imim,imip 
     !**************************************!
     !  Correction for magnitude of the VA  !
     !**************************************!
     IF(elm(i)%h_mag.EQ.-9.99d0)THEN
        d_sun_VA=vsize(elm(i)%coord(1:3))      ! Heliocentric distance of the VA
        el_geo_VA=elm(i)%coord(1:3)-xear(1:3)  ! Geocentric elements of the VA
        d_ear_VA=vsize(el_geo_VA)              ! Geocentric distance of the VA
        ! Correction for the distances
        elm(i)%g_mag=elm(imi0)%g_mag           ! Slope parameter G
        elm(i)%h_mag=elm(imi0)%h_mag + 5.d0*LOG10((d_ear_nom*d_sun_nom)/(d_ear_VA*d_sun_VA))
        ! Correction depending upon the phase angle of the VA
        pha_VA=ACOS(DOT_PRODUCT(elm(i)%coord(1:3),el_geo_VA)/(d_sun_VA*d_ear_VA))
        phi1=EXP(-3.33d0*TAN(pha_VA/2.d0)**0.63d0) 
        phi2=EXP(-1.87d0*TAN(pha_VA/2.d0)**1.22d0)
        f_alpha_VA=2.5d0*LOG10((1.d0-elm(imi0)%g_mag)*phi1+elm(imi0)%g_mag*phi2) 
        elm(i)%h_mag=elm(i)%h_mag+f_alpha_nom-f_alpha_VA
     END IF
     !*****************!
     !  Output report  !
     !*****************!
     ! Orbital distance
     moid_min=MIN(moid_m(i),moid_min) 
     ! Sigma of the i-th multiple solution
     sigmalov(i)=(i-(imult+1))*deltasigma
     ! Report output
     WRITE(iunrep,145) i,csinom(i),delnom(i),elm(i)%h_mag,moid_m(i),dnp_m(i),dnm_m(i),sigmalov(i),deltasigma,sigq(i)
145  FORMAT(I5,1X,1P,E13.5,E11.3,2X,0P,F5.2,1X,F8.5,1X,F10.5,1X,F10.5,3X,F12.8,4X,F12.8,3X,F12.8)
     astna0=' '
     IF(i.LE.9)THEN              
        WRITE(astna0,111) i 
111     FORMAT('va_00000',I1)
     ELSEIF(i.LE.99)THEN
        WRITE(astna0,112) i
112     FORMAT('va_0000',I2)
     ELSEIF(i.LE.999)THEN
        WRITE(astna0,113) i
113     FORMAT('va_000',I3)
     ELSEIF(i.LE.9999)THEN
        WRITE(astna0,114) i
114     FORMAT('va_00',I4)
     ELSEIF(i.LE.99999)THEN
        WRITE(astna0,115) i
115     FORMAT('va_0',I5)
     ELSEIF(i.LE.999999)THEN
        WRITE(astna0,116) i
116     FORMAT('va_',I6)
     ELSE
        WRITE(*,*)' nmulti: change format for imult>499999'
        STOP
     ENDIF
     CALL rmsp(astna0,le)
     !****************************!
     !  Output multiline catalog  !
     !****************************!
     IF(no_outcov.AND.i.NE.imi0)THEN       
        IF(dyn%ndp.GT.0)THEN 
           IF(nd.GT.6)THEN
              dyn%dp(1:dyn%ndp)=dpm(1:dyn%ndp,i)
           END IF
           CALL write_elems(elm(i),astna0(1:le),'ML',dyn,nls,ls,UNIT=iunctc)
        ELSE
           CALL write_elems(elm(i),astna0(1:le),'ML',UNIT=iunctc)
        END IF
     ELSE
        IF(dyn%ndp.GT.0)THEN 
           IF(nd.GT.6)THEN
              dyn%dp(1:dyn%ndp)=dpm(1:dyn%ndp,i)
           END IF
           CALL write_elems(elm(i),astna0(1:le),'ML',dyn,nls,ls,unm(i),UNIT=iunctc)
        ELSE
           CALL write_elems(elm(i),astna0(1:le),'ML',UNC=unm(i),UNIT=iunctc)
        END IF
     ENDIF
     !****************************************!
     !  Output one line catalog, if possible  !
     !****************************************!
     CALL coo_cha(elm(i),'KEP',elk,fail_flag)
     IF(fail_flag.GE.4)THEN
        WRITE(iun_log,*) ' nmulti: hyperbolic orbit in output ', fail_flag, i     
     ELSE
        IF(dyn%ndp.GT.0)THEN 
           IF(nd.GT.6)THEN
              dyn%dp(1:dyn%ndp)=dpm(1:dyn%ndp,i)
              ngr_case=ngr_case+1
           END IF
           CALL write_elems(elk,astna0(1:le),'1L',dyn,nls,ls,FILE=' ',UNIT=iuncat,UNIT2=iuncatngr)
        ELSE
           CALL write_elems(elk,astna0(1:le),'1L',FILE=' ',UNIT=iuncat)
        ENDIF
     END IF
  ENDDO
  !*********************!
  !  Terminate cleanly  !
  !*********************!
  CALL filclo(iunctc,' ') 
  CALL filclo(iunrep,' ') 
  CALL filclo(iuncat,' ') 
  IF(ngr_case.EQ.0)THEN
     CALL filclo(iuncatngr,'DELETE') 
  ELSE
     CALL filclo(iuncatngr,' ')
  END IF
END SUBROUTINE nmulti


! =============================================================================          
!  NMULTI_GEN interface to (distribution) multi_gen
! =============================================================================           
! multiple solutions for a differential correction problem: adaptive step size              
! ===============INTERFACE=====================================================
SUBROUTINE nmulti_gen(nam0,elc,uncert,csinor,delnor,            &
     &     mc,obs,obsw,prob_step,sigma_step_max,sigma_max,moid_min,nd) 
! ================INPUT========================================================
  USE name_rules, ONLY: name_len
! begin interface
! number of observations
  INTEGER, INTENT(IN)                      ::  mc
! new data types
  TYPE(ast_obs),DIMENSION(mc), INTENT(IN)  :: obs
  TYPE(ast_wbsr),DIMENSION(mc), INTENT(IN) :: obsw                         
! asteroid name                                                         
  CHARACTER*(name_len), INTENT(IN)         ::  nam0 
! initial conditions: epoch, elements, abs. magnitude, gmag             
  TYPE(orbit_elem), INTENT(IN)             :: elc 
! normal and covariance matrices
  TYPE(orb_uncert), INTENT(IN)             :: uncert 
! norms of residuals, of last correction
  DOUBLE PRECISION, INTENT(IN)             :: csinor,delnor 
! options for generation of multiple solutions: constant step in IP
  DOUBLE PRECISION, INTENT(IN)             :: prob_step, sigma_step_max, sigma_max
! dimension of solve for parameters vector
  INTEGER, INTENT(IN)                      :: nd 
! end interface
!*******************************************************************************
! ==========OUTPUT==============                                        
  DOUBLE PRECISION, INTENT(OUT) ::  moid_min 
! renormalization factor                                                
  DOUBLE PRECISION rescov 
! local copy of dynamical parameters
  DOUBLE PRECISION dp(ndyx)
! logical variables used for call to fmulti                             
  LOGICAL obsc,inic,covc,ok 
! ==============END INTERFACE=====================                      
! loop indexes                                                          
  INTEGER i 
  LOGICAL batch ! batch control
  LOGICAL linlov,nooutcov
! catalog obj name (with underscore)                                    
  CHARACTER(LEN=20):: astna0 
  INTEGER le
! output files                                                          
  CHARACTER*60 file 
  INTEGER iunctc,iunrep,iuncat,iuncatngr,iwdir
  TYPE(orbit_elem) :: elk
  INTEGER fail_flag
  CHARACTER*6 vers
  INTEGER ngr_case
  INTEGER :: imult
!===== Magnitude correction ==============================================================================
  DOUBLE PRECISION                  :: vsize                   ! Norm of a vector (function)
  DOUBLE PRECISION                  :: d_ear_nom               ! Geocentric distance of the nominal
  DOUBLE PRECISION                  :: d_sun_nom               ! Heliocentric distance of the nominal
  DOUBLE PRECISION                  :: d_ear_VA                ! Geocentric distance of the VA
  DOUBLE PRECISION                  :: d_sun_VA                ! Heliocentric distance of the VA
  DOUBLE PRECISION                  :: phi1, phi2              ! Computation of the phase angle correction
  DOUBLE PRECISION                  :: f_alpha_nom, f_alpha_VA ! Phase angle correction
  DOUBLE PRECISION                  :: xear(6)                 ! State vector of the Earth
  DOUBLE PRECISION                  :: el_geo(3)               ! Geocentric elements of the nominal
  DOUBLE PRECISION                  :: el_geo_VA(3)            ! Geocentric elements of the VA
  DOUBLE PRECISION                  :: pha_nom, pha_VA         ! Phase angles (nominal and VAs)
  !=========================================================================================================
! set flags to true to avoid checks                                     
  obsc=.true. 
  inic=.true. 
  covc=.true. 
! batch mode                                                            
  batch=.true. 
! selection of LOV definition
  IF(lin_lov)THEN
     IF(scatterplane)THEN
        CALL filnam('multsol',nam0,'wdir',file,le) 
        CALL filopn(iwdir,file(1:le),'unknown') 
     ENDIF
! maximum interface version
!     CALL lin_multi(batch,obsc,inic,ok,covc,elc,uncert,csinor,delnor,mc,obs,obsw,sigma,imult1,nd)
! minimum interface/computation version
     CALL lin_multi_gen(batch,elc,uncert,prob_step,sigma_step_max,sigma_max,nd,imult,iwdir,mc,obs,obsw,csinor)
     IF(scatterplane) CALL filclo(iwdir,' ')
  ELSE
 ! Constant step in sigma: call to multi_gen (distribution version, but with batch mode)
     CALL multi_gen(batch,elc,uncert,csinor,delnor,mc,obs,obsw,   &
          &     prob_step,sigma_step_max,sigma_max,nd,imult)
  ENDIF
! ================================================                      
! output 
! count cases with non-grav parameters
  ngr_case=0                                                               
! multiple solutions catalog: multiline format                          
  CALL filnam('multsol',nam0,'ctc',file,le) 
  CALL filopn(iunctc,file(1:le),'unknown') 
  vers='OEF2.0'
  CALL wromlh2(iunctc,'ECLM','J2000',vers) 
! single line format KEP                                                    
  CALL filnam('multsol',nam0,'cat',file,le) 
  CALL filopn(iuncat,file(1:le),'unknown') 
  CALL wro1lh2(iuncat,'ECLM','J2000','KEP',vers) 
! single line format for non-grav param
  CALL filnam('multsol',nam0,'cat.ngr',file,le) 
  CALL filopn(iuncatngr,file(1:le),'unknown') 
  CALL wro1lh2_ngr(iuncatngr,'ECLM','J2000',vers,dyn%ndp) 
! report file (with MOID, RMS, etc.)                                    
  CALL filnam('multsol',nam0,'mrep',file,le) 
  CALL filopn(iunrep,file(1:le),'unknown') 
! information to be passed to resret2
  WRITE(iunrep,399)
399 FORMAT('MULTIPLE SOLUTIONS COMPUTATION, prob_sampl=.TRUE.')
  WRITE(iunrep,199) prob_step,sigma_max,sigma_step_max, imult
199 FORMAT('MULTIPLE SOLUTIONS SUMMARY, prob_step=',1P,E7.1,0P,' sigma_max=',F5.3, ' sigma_step_max=',F5.3,' imult=',I5) 
  WRITE(iunrep,299)scaling_lov,second_lov
299 FORMAT('Scaling LOV= ',L1/'Second  LOV= ',L1)
! header
  WRITE(iunrep,499)
499 FORMAT('   no      RMS        lastcor    mag    MOID       nod+       nod-       sigma_lov      delta_sigma        sigQ')
! loop on multiple solutions effectively computed                       
  moid_min=1.d3 
  
  !***************************************!
  !  Parameters for magnitude correction  !
  !***************************************!
  d_sun_nom=vsize(elm(imi0)%coord(1:3))  ! Heliocentric distance of the nominal
  CALL earcar(elm(imi0)%t,xear,1)
  el_geo=elm(imi0)%coord(1:3)-xear(1:3)  ! Geocentric elements of the nominal solution
  d_ear_nom=vsize(el_geo)                ! Geocentric distance of the nominal solution
  pha_nom=ACOS(DOT_PRODUCT(elm(imi0)%coord(1:3),el_geo)/(d_sun_nom*d_ear_nom)) ! Phase angle
  ! Correction depending upon the phase angle
  phi1=EXP(-3.33d0*TAN(pha_nom/2.d0)**0.63d0) 
  phi2=EXP(-1.87d0*TAN(pha_nom/2.d0)**1.22d0)
  f_alpha_nom=2.5d0*LOG10((1.d0-elm(imi0)%g_mag)*phi1+elm(imi0)%g_mag*phi2) 
  
  DO i=imim,imip 
     !**************************************!
     !  Correction for magnitude of the VA  !
     !**************************************!
     IF(elm(i)%h_mag.EQ.-9.99d0)THEN
        d_sun_VA=vsize(elm(i)%coord(1:3))      ! Heliocentric distance of the VA
        el_geo_VA=elm(i)%coord(1:3)-xear(1:3)  ! Geocentric elements of the VA
        d_ear_VA=vsize(el_geo_VA)              ! Geocentric distance of the VA
        ! Correction for the distances
        elm(i)%g_mag=elm(imi0)%g_mag           ! Slope parameter G
        elm(i)%h_mag=elm(imi0)%h_mag + 5.d0*LOG10((d_ear_nom*d_sun_nom)/(d_ear_VA*d_sun_VA))
        ! Correction depending upon the phase angle of the VA
        pha_VA=ACOS(DOT_PRODUCT(elm(i)%coord(1:3),el_geo_VA)/(d_sun_VA*d_ear_VA))
        phi1=EXP(-3.33d0*TAN(pha_VA/2.d0)**0.63d0) 
        phi2=EXP(-1.87d0*TAN(pha_VA/2.d0)**1.22d0)
        f_alpha_VA=2.5d0*LOG10((1.d0-elm(imi0)%g_mag)*phi1+elm(imi0)%g_mag*phi2) 
        elm(i)%h_mag=elm(i)%h_mag+f_alpha_nom-f_alpha_VA
     END IF
     
     ! orbital distance                                                      
     moid_min=min(moid_m(i),moid_min) 
     WRITE(iunrep,145)i,csinom(i),delnom(i),elm(i)%h_mag,moid_m(i),dnp_m(i),dnm_m(i),sigmavalv(i),deltasigmav(i),sigq(i)
145  FORMAT(I5,1X,1P,E13.5,E11.3,2x,0P,F5.2,1X,F8.5,1X,F10.5,1X,F10.5,3X,F12.8,4X,F12.8,3X,F12.8)
     astna0=' '
     IF(i.le.9)THEN              
        WRITE(astna0,111) i 
111     FORMAT('va_00000',I1)
     ELSEIF(i.le.99)THEN
        WRITE(astna0,112) i
112     FORMAT('va_0000',I2)
     ELSEIF(i.le.999)THEN
        WRITE(astna0,113) i
113     FORMAT('va_000',I3)
     ELSEIF(i.le.9999)THEN
        WRITE(astna0,114) i
114     FORMAT('va_00',I4)
     ELSEIF(i.le.99999)THEN
        WRITE(astna0,115) i
115     FORMAT('va_0',I5)
     ELSEIF(i.le.999999)THEN
        WRITE(astna0,116) i
116     FORMAT('va_',I6)
     ELSE
        WRITE(*,*)' nmulti: change format for imult>499999'
        STOP
     ENDIF
     CALL rmsp(astna0,le)
     ! output multiline catalog 
     IF(no_outcov.and.i.ne.imi0)THEN       
        IF(dyn%ndp.gt.0)THEN 
           IF(nd.GT.6)THEN
              dyn%dp(1:dyn%ndp)=dpm(1:dyn%ndp,i)
           END IF
           CALL write_elems(elm(i),astna0(1:le),'ML',dyn,nls,ls,UNIT=iunctc)
        ELSE
           CALL write_elems(elm(i),astna0(1:le),'ML',UNIT=iunctc)
        END IF
     ELSE
        IF(dyn%ndp.gt.0)THEN 
           IF(nd.GT.6)THEN
              dyn%dp(1:dyn%ndp)=dpm(1:dyn%ndp,i)
           END IF
           CALL write_elems(elm(i),astna0(1:le),'ML',dyn,nls,ls,unm(i),UNIT=iunctc)
        ELSE
           CALL write_elems(elm(i),astna0(1:le),'ML',UNC=unm(i),UNIT=iunctc)
        END IF
     ENDIF
! output one line catalog, if possible
     CALL coo_cha(elm(i),'KEP',elk,fail_flag)
     IF(fail_flag.ge.4)THEN
        WRITE(iun_log,*)' hyperbolic orbit in output ', fail_flag, i     
     ELSE
        IF(dyn%ndp.gt.0)THEN 
           IF(nd.GT.6)THEN
              dyn%dp(1:dyn%ndp)=dpm(1:dyn%ndp,i)
              ngr_case=ngr_case+1
           END IF
           CALL write_elems(elk,astna0(1:le),'1L',dyn,nls,ls,FILE=' ',UNIT=iuncat,UNIT2=iuncatngr)
        ELSE
           CALL write_elems(elk,astna0(1:le),'1L',FILE=' ',UNIT=iuncat)
        ENDIF
     END IF
  ENDDO
  CALL filclo(iunctc,' ') 
  CALL filclo(iunrep,' ') 
  CALL filclo(iuncat,' ') 
! Delete .cat.ngr if empty
  IF(ngr_case.EQ.0)THEN
     CALL filclo(iuncatngr,'DELETE') 
  ELSE
     CALL filclo(iuncatngr,' ')
  END IF
!=============================================                          
END SUBROUTINE nmulti_gen

! Copyright 2006, The Orbfit Consortium
!================================================================!                  
! STEP_FIT2                                                      !
!================================================================!
! Differential corrector slipping along the LOV by steps         !
!================================================================!
! Version 3.3.2, 10 Jan 2006, with prop_sig for smooth slipping. ! 
!                                                                !
! Remark. The weights are constant in this routine.              !
!================================================================!
SUBROUTINE step_fit2(m,obs,obsw,el0,unc0,csino0,delno0,nused,succ)
  USE pred_obs
  USE close_app, ONLY: kill_propag
  !=========================================================================================================
  INTEGER,                        INTENT(IN)    :: m      ! Number of observations
  TYPE(ast_obs),    DIMENSION(m), INTENT(IN)    :: obs    ! Observations
  TYPE(ast_wbsr),   DIMENSION(m), INTENT(INOUT) :: obsw   ! Weights/biases
  TYPE(orbit_elem),               INTENT(INOUT) :: el0    ! Orbital elements (assumed to be on the LOV)
  TYPE(orb_uncert),               INTENT(INOUT) :: unc0   ! Uncertainty matrices
  DOUBLE PRECISION,               INTENT(INOUT) :: delno0 ! Norm of last correction
  DOUBLE PRECISION,               INTENT(INOUT) :: csino0 ! Residuals norm
  INTEGER,                        INTENT(OUT)   :: nused  ! Number of obs. used
  LOGICAL,                        INTENT(OUT)   :: succ   ! Success flag (the solution is the nominal, or a
  !                                                         better one along LOV, if succ=.TRUE.)
  !===== Input data sorted =================================================================================
  INCLUDE 'parobx.h90' 
  TYPE(ast_obs),  DIMENSION(nobx) :: obs_s       ! Sorted observations
  TYPE(ast_wbsr), DIMENSION(nobx) :: obsw_s      ! Sorted weights/biases
  INTEGER                         :: iposs(nobx) ! Indexes for the sort
  INTEGER                         :: icor6(6)    ! Parameters to be solved
  TYPE(orbit_elem)                :: elc         ! Corrected elements
  !===== Loop indexes ======================================================================================
  INTEGER                         :: it          ! Counter of all iterations
  INTEGER                         :: itst        ! Counter of LOV steps only (limited by itstmax)
  INTEGER                         :: itc         ! Constrained correction iterations (limited by itcmax)
  INTEGER                         :: ir          ! Relaxation loop for correction
  INTEGER                         :: i,j         ! Indexes
  INTEGER, PARAMETER              :: itcmax=10   ! Max n. of constrained iterations
  INTEGER, PARAMETER              :: itstmax=20  ! Max n. of iterations along LOV
  !===== Corrections =======================================================================================
  DOUBLE PRECISION                :: deq6_lov(6)   ! Component along the LOV
  DOUBLE PRECISION                :: deq6_const(6) ! Component orthogonal to the LOV
  DOUBLE PRECISION                :: peq(6)        ! Weak direction
  DOUBLE PRECISION                :: deq6(6)       ! Correction of sin_cor
  DOUBLE PRECISION                :: deq(6)        ! Help variable for deq6
  DOUBLE PRECISION                :: sdir          ! Eigenvalues of the weak direction
  DOUBLE PRECISION                :: units(6)      ! Units for scaling
  TYPE(orb_uncert)                :: uncert        ! Uncertainty matrices (for sin_cor)
  DOUBLE PRECISION                :: delnor        ! Norm of last correction (for sin_cor)
  DOUBLE PRECISION                :: csinor        ! Residuals norm (for sin_cor)
  DOUBLE PRECISION                :: csino1        ! Storage of csino0
  DOUBLE PRECISION                :: delnor_lov    ! Norm of the correction along the LOV
  DOUBLE PRECISION                :: delnor_const  ! Norm of the correction orthogonal to the LOV
  !===== Functions =========================================================================================
  DOUBLE PRECISION                :: snorm         ! Function for the norm induced byy a matrix
  DOUBLE PRECISION                :: rescov        ! Function for scaling
  !===== Other variables ===================================================================================
  DOUBLE PRECISION                :: ecc           ! Eccentricity control for bizarre
  TYPE(orbit_elem)                :: elb           ! Temporary elements
  INTEGER                         :: nsolv         ! Number of parameters to be solved
  LOGICAL                         :: lov_step      ! Logical flag: if .TRUE., continue on the LOV
  LOGICAL                         :: bizarre       ! Bizarre control
  INTEGER                         :: nd            ! Solve-for parameters space dimension
  DOUBLE PRECISION                :: mu            ! Scaling factor
  DOUBLE PRECISION                :: dn            ! Sigma step (not used anymore)
  LOGICAL                         :: fail          ! Fail flag for prop_sig (not used anymore)
  !=========================================================================================================
  !******************!
  !  Initialization  !
  !******************!
  icor6=1 
  CALL blockset(m,obs,obsw)                          ! Definition of blocks if necessary
  CALL sort_obs(el0%t,obs,obsw,m,iposs,obs_s,obsw_s) ! Sort of observations time (w.r.t. el0%t)                 
  succ=.FALSE.
  ! Storage of the elements and of the residuals norm
  elc=el0 
  csino1=csino0 
  ! Verbosity control
  IF(verb_dif.GT.9) WRITE(iun_log,*)'starting values ', el0 
  IF(verb_dif.GT.19) WRITE(*,*)'starting values ',  el0 

  !=================================================================!
  !===== MAIN LOOP ON STEP ITERATIONS ==============================!
  !=================================================================!
  it=0   ! Counter of all iterations
  itst=0 ! Counter of LOV steps only (limited by itstmax)
  itc=0  ! Constrained correction iterations (limited by itcmax)
60 CONTINUE 
  ! Control on non-gravitational parameters
  nd=6+nls
  IF(nd.GT.6)THEN
     WRITE(*,*) ' step_fit2 not ready for non-grav, nd =',nd
     STOP
  ENDIF
  !*******************************************!
  !  Single step of differential corrections  !
  !*******************************************!
  CALL sin_cor(m,obs_s,obsw_s,elc,icor6,iun_log,.TRUE.,delnor,csinor,uncert,nd,deq6)
  IF(.NOT.uncert%succ) GOTO 9
  it=it+1 
  ! Compute weak direction (peq)
  CALL weak_dir(uncert%g(1:nd,1:nd),peq,sdir,-1,elc%coo,elc%coord,units,nd)
  deq6=deq6/units ! Scaling of the coordinates
  !********************************************!
  !  Compute the components of the correction  !
  !********************************************!
  IF(DOT_PRODUCT(peq,deq6).LT.0.d0) peq=-peq  ! Sign rule to decrease Q 
  deq6_lov=DOT_PRODUCT(peq,deq6)*peq          ! Component along the LOV
  deq6_const=deq6-deq6_lov                    ! Component orthogonal to the LOV
  deq6_lov=deq6_lov*units                     ! deq6_lov in unscaled coordinates
  deq6_const=deq6_const*units                 ! deq6_const in unscaled coordinates
  delnor_lov=snorm(deq6_lov,uncert%c,6,6)     ! Norm of the correction along the LOV
  delnor_const=snorm(deq6_const,uncert%c,6,6) ! Norm of the component normal to the LOV
  !**************************************************!
  !  Check if we need another constrained iteration  !
  !**************************************************!
  IF(delnor_const.GE.del_constr)THEN
     ! The norm of the orthogonal component is still too large
     ! Need to get back to the LOV
     itc=itc+1
     IF(verb_dif.GE.19) WRITE(*,*) ' constr.corr., delnor_const =',delnor_const
     IF(verb_dif.GE.9) WRITE(iun_log,*) ' constr.corr., delnor_const =',delnor_const
     lov_step=.FALSE.
     deq=deq6_const ! Save constrained corrections  
  ELSEIF(delnor_lov.LT.delcr)THEN
     ! NO, and also no need to move along the LOV -> Exit from the loop
     succ=.TRUE.
     IF(verb_dif.GE.19) WRITE(*,*) ' success, nominal solution, RMS, norms ',delnor_const,delnor_lov
     WRITE(iun_log,*) ' success, nominal solution, RMS, norms ',delnor_const,delnor_lov
     GOTO 70
  ELSE
     ! NO, compare the current result with the stored one
     IF(csinor.LT.csino1)THEN
        el0=elc
        unc0=uncert
        csino0=csinor
        delno0=delnor
        csino1=csinor
        succ=.TRUE.
        IF(verb_dif.GE.19) WRITE(*,*) ' success, LOV solution, RMS, norms ',csinor,delnor_const,delnor_lov
        WRITE(iun_log,*) ' success, LOV solution, RMS, norms ',csinor,delnor_const,delnor_lov
     ENDIF
     ! Try and move along the LOV again
     itc=0       ! Restart counter of constr. steps for the next time
     itst=itst+1 ! LOV step
     ! Use unconstrained corrections, of limited length in the weak direction
     IF(verb_dif.GE.19) WRITE(*,*) ' full corr., delnor_lov=',delnor_lov
     IF(verb_dif.GE.9) WRITE(iun_log,*) ' full corr., delnor_lov=',delnor_lov
     lov_step=.TRUE.
     IF(delnor_lov.GT.step_sig)THEN
        !dn=step_sig
        deq=deq6_const+deq6_lov*(step_sig/delnor_lov)
     ELSE
        !dn=delnor_lov
        deq=deq6_const+deq6_lov
     ENDIF
  ENDIF
  !******** ALTERNATIVE METHOD *********
  !     elb=elc
  !     CALL prop_sig(.true.,elb,elc,dn,sdir,m,obs,obsw,peq,sdir,units,fail)
  !***** END ALTERNATIVE METHOD ********
  !*************************!
  !  Appply the correction  !
  !*************************!
  elb=elc
  DO ir=1,4
     elb%coord=elc%coord+deq 
     IF(bizarre(elb,ecc))THEN
        IF(verb_dif.ge.19)WRITE(*,*) ' constr_fit: short step to avoid ecc =',ecc
        IF(verb_dif.ge.9)WRITE(iun_log,*) ' constr_fit: short step to avoid ecc =',ecc
        deq=deq/2.d0
        elb%coord=elc%coord+deq 
     ELSE
        EXIT
     ENDIF
  ENDDO
  elc=elb
  ! Full norm of the actual correction used
  delnor=snorm(deq,uncert%c,6,6) 
  IF(verb_dif.GE.9) WRITE(iun_log,200)it,itc,csinor,delnor,elc%coord
  IF(verb_dif.GE.19) WRITE(*,200)it,itc,csinor,delnor,elc%coord 
200 FORMAT(' *** iteration ',i3,2x,I3,' RMS residuals =',1p,d12.4,    &
         &        '   norm corr =',d12.4,'  new elem values:'/0p,6f13.7/)
  ! Control against hyperbolic and bizarre orbits                         
  IF(bizarre(elc,ecc))THEN 
     IF(verb_dif.GE.19) WRITE(*,*) ' step_fit: iter. ',it,' bizarre; e =',ecc
     IF(verb_dif.GE.9) WRITE(iun_log,*) ' step_fit: iter. ',it,' bizarre; e =',ecc
     RETURN
  ENDIF
  ! Check if iterations can go on
  IF(itc.GT.itcmax.OR.itst.GT.itstmax)THEN
     GOTO 70 ! Interrupt attempt
  ELSE
     GOTO 60 ! Keep trying
  ENDIF
70 CONTINUE
  !=================================================================!
  !===== END OF THE MAIN LOOP ======================================!
  !=================================================================!

  IF(.NOT.succ)THEN
     IF(verb_dif.GE.19)THEN
        WRITE(*,*) ' non convergent after ',it,itst,itc,' iterations ' 
     ELSEIF(verb_dif.GT.2)THEN
        WRITE(iun_log,*) ' non convergent after ',it,itst,itc,' iterations '
     ENDIF
  ENDIF
  !******************************************!
  !  Covariance and normal matrix rescaling  !
  !******************************************!
  nsolv=6 
  ! Count of the observations used
  nused=0 
  DO i=1,m 
     IF(obsw_s(i)%sel_coord.GT.0)THEN
        IF(obs_s(i)%type.EQ.'O'.OR.obs_s(i)%type.EQ.'S')THEN
           nused=nused+2
        ELSEIF(obs_s(i)%type.EQ.'R'.OR.obs_s(i)%type.EQ.'V')THEN
           nused=nused+1
        ENDIF
     ENDIF
  ENDDO
  ! Apply rescaling to both covariance and normal matrix 
  mu=rescov(nsolv,nused,csinor) 
  uncert%g=uncert%g*mu**2 
  uncert%c=uncert%c/mu**2 
  ! Output final result, with norm of the residuals                            
  IF(verb_dif.GE.9.AND.verb_dif.LT.19) WRITE(*,201) it,csinor,delnor 
201 FORMAT(' done constr. iter. ',i3,'  RMS =',1p,d12.4,' last corr. =',d12.4)
  ! Reordering the residuals for output                                   
  CALL unsort_obs(iposs,m,obsw_s,obsw)
  RETURN 
  ! Non convergent differential correction (sin_cor): impossible normal matrix inversion.
9 CONTINUE 
  IF(verb_mul.GE.9)THEN
     WRITE(*,*) ' step_fit2: inversion failed '
     WRITE(iun_log,*) ' step_fit2: inversion failed '
  ENDIF
  succ=.FALSE.
  csinor=csino0
END SUBROUTINE step_fit2


! ==============================================================================================                 
! mult_input initializes storage of multiple solution data  (used only in fitobs)            
! ==============================================================================================
SUBROUTINE mult_input(catname,ok,sigma_max,nd,sigma_step_max)  
  USE name_rules
! ------------INPUT------------------                                   
  CHARACTER(LEN=160), INTENT(IN) :: catname ! file with catalog of multiple solutions 
  INTEGER,            INTENT(IN) :: nd      ! number of solve for parameteres
! ------------OUTPUT-----------------                                   
  LOGICAL,                  INTENT(OUT) :: ok             ! Succesful input  
  DOUBLE PRECISION,         INTENT(OUT) :: sigma_max      ! Maximum value of sigma
  DOUBLE PRECISION,OPTIONAL,INTENT(OUT) :: sigma_step_max ! Maximum step in sigma
                                                          ! (If prob_sampl=false, sigma_step_max 
                                                          !  is equal to delta_sigma)
!***********************************************************************************************
  DOUBLE PRECISION :: v_infty 
  INTEGER          :: ndp,nmod,nlsloc,lsloc(1:ndyx)   ! LSP record
  INTEGER          :: ndloc   ! local dimension (paranoia check)
  CHARACTER*4      :: form ! 1L, ML
! index range of multiple solution
  INTEGER :: m1,m2,m0                                           
! names (for call to rdorb)                                                
  CHARACTER*(idname_len) :: name0 
  CHARACTER*(name_len)   :: name1 
  INTEGER                :: le
  CHARACTER*160          :: catname1, repname, ctcfile
! units for multiple solution catalog input 
  INTEGER :: iuncat,iunrep
! end of file, error in reading
  LOGICAL :: eof,err
! orbit elements, covariance, nongrav parameters
  TYPE(orbit_elem) :: eltmp
  TYPE(orb_uncert) :: untmp
  DOUBLE PRECISION :: dptmp(ndyx) 
! indexes of multiple solutions                                         
  INTEGER :: imul(mulx),norb,j, imult 
! velocity w.r.to Earth for each orbit                                  
  DOUBLE PRECISION :: v_max,v_min 
! temporary input arrays                                                
  DOUBLE PRECISION :: eq(6) 
  DOUBLE PRECISION :: g(nd,nd),c(nd,nd),h, zero
! loop indexes                                                          
  INTEGER :: i,ii 
! LOV computation
  LOGICAL :: second_lov,scaling_lov 
! system dependencies                                                   
  INCLUDE 'sysdep.h90' 
! -------------------------------------------------                     
  v_min=100.d0 
  v_max=0.d0 
! multiple solution report file    
  CALL rmsp(catname,le)                       
  repname=catname(1:le)//'.mrep' 
  CALL rmsp(repname,le) 
  CALL filopn(iunrep,repname(1:le),'old') 
!Reading mrep: prob_sampl
  READ(iunrep,5)prob_sampl
5 FORMAT(43X,L7)
  WRITE(*,55) prob_sampl
55 FORMAT(' prob_sampl  = ', L1)
  IF(prob_sampl)THEN
     IF(PRESENT(sigma_step_max))THEN
        READ(iunrep,399)prob_step, sigma_max, sigma_step_max, imult 
399     FORMAT(38X,E7.1,11X,F5.3,16X,F5.3,7X,I5)
        WRITE(*,400)prob_step,sigma_max, sigma_step_max, imult
400     FORMAT(' prob_step      = ',1P,E7.1,0P/,' sigma max      = ',F5.3/,' sigma step max = ',F5.3/,' imult         = ',I5) 
     ELSE
        WRITE(*,*) 'mult_input: missing argument sigma_step_max when prob_sampl is true'
        STOP
     ENDIF
  ELSE
     READ(iunrep,199)imult,sigma_max
199  FORMAT(34X,I5,10X,F5.3) 
     WRITE(*,200)imult, sigma_max
200  FORMAT(' imult       =',I5/,' sigma max   = ',F5.3) 
     sigma_step_max=sigma_max/imult
  END IF
! Reading scaling LOV and second LOV
  READ(iunrep,299)scaling_lov,second_lov
299 FORMAT(13X,L1/13X,L1)
  WRITE(*,499)scaling_lov,second_lov
499 FORMAT(' scaling_lov = ',L1/,' second_lov  = ',L1)
! Skipping header 
  READ(iunrep,*)
  imi0=imult+1 
! opening and reading multiple solution catalog                         
  catname1=catname
  CALL rmsp(catname1,le) 
  ctcfile=catname1(1:le)//'.ctc'
  CALL rmsp(ctcfile,le)
  INQUIRE(FILE=ctcfile(1:le),exist=ok) 
  IF(.not.ok)THEN 
     WRITE(*,*)'MULT_INPUT: ctc file requested. File ',ctcfile(1:le),' not found' 
     RETURN 
  ENDIF
  CALL oporbf(ctcfile(1:le),iuncat) 
  DO i=1,mulx 
     CALL read_elems(eltmp,name0,eof,nlsloc,err,lsloc,form,untmp,UNIT=iuncat)
     IF(err)THEN
        WRITE(*,*) 'mult_input: error in reading VA number i=',i
        STOP ' error in reading VA '
     ENDIF
!Save dyn%dp, nls, and ls
     IF(dyn%ndp.gt.0)THEN
        nls=nlsloc
        ls=lsloc
     ENDIF
! Dimension
     ndloc=6+nls
     IF(ndloc.NE.nd)THEN
        WRITE(*,*)'mult_input: ndloc is ', ndloc,' while nd is ', nd
        STOP
     END IF
     IF(eof)THEN 
        norb=i-1 
        GOTO 2 
     ENDIF
! control on time                                                       
     IF(i.eq.1)THEN 
        tdt_cat=eltmp%t 
     ELSE 
        IF(eltmp%t.ne.tdt_cat)THEN 
           WRITE(*,*)'mult_input: time discrepancy from tcat=',tdt_cat 
           WRITE(*,*)'mult_input: at record ',i,' t=',eltmp%t 
           ok=.false. 
           RETURN 
        ENDIF
     ENDIF
! availability of covariance                                            
     IF (.NOT.(untmp%succ)) THEN 
        WRITE(*,*)'mult_input',  name0, ' matrices not avalaible' 
        ok=.false. 
        RETURN 
     ENDIF
! handling of name; but note that the multiple solutions are assumed to 
! and with consecutive indexes!!! that is imul(i)=imul(1)+i-1           
     CALL splinam(name0,name1,imul(i)) 
     IF(imul(i).ne.imul(1)+i-1)THEN 
        WRITE(*,*)'mult_input: indexes not in order: ',             &
     &           imul(i),' at record ',i,' should be ',imul(1)+i-1      
!            RETURN                                                     
     ENDIF
     j=imul(i) 
! copy in output arrays
     elm(j)=eltmp
     unm(j)=untmp
     IF(dyn%ndp.gt.0) dpm(1:dyn%ndp,j)=dyn%dp(1:nd)
! Read mrep
     READ(iunrep,145)ii,csinom(j),delnom(j),h,                   &
          &  moid_m(j),dnp_m(j),dnm_m(j), sigmavalv(j), deltasigmav(j), sigq(j)
145  FORMAT(I5,1X,1P,E13.5,E11.3,2X,0P,F5.2,1X,                     &
          &  F8.5,1X,F10.5,1X,F10.5,3X,F12.8,4X,F12.8,3X,F12.8)
! v_infinity computation                                                
     v_inf(i)=v_infty(elm(j)) 
     v_min=min(v_min,v_inf(i)) 
     v_max=max(v_max,v_inf(i)) 
  ENDDO
! increase mulx                                                        
  WRITE(*,*)'mult_input: increase mulx, file is longer than',mulx 
33 CALL filclo(iunrep,' ')
2 CONTINUE 
  norb=i-1 
  ok=norb.gt.0 
  imim=imul(1) 
  imip=imul(norb) 
  WRITE(*,*)' mult_input: input of ',norb,' multiple solutions' 
  WRITE(*,*)' indexes between ', imim, ' and ', imip 
  WRITE(*,*)' max of v_infty ',v_max,' min ',v_min 
3 WRITE(*,*)' which one is the nominal?' 
  READ(*,*)imi0 
  IF(imi0.lt.imim.or.imi0.gt.imip)THEN 
     WRITE(*,*)' OUT OF RANGE ',imim,imip 
     GOTO 3 
  ELSE 
     WRITE(*,*)' nominal is ',imi0 
  ENDIF
END SUBROUTINE mult_input
!
! =======================================                               
!  FMUOBS                                                               
! =======================================                               
! output of multiple observations                                       
! ===============INTERFACE========================                      
SUBROUTINE fmuobs(type,ids,t1,tut1,sigma,         &
     &     aobs,dobs,iff,titnam,filnam,iun20) 
  USE pred_obs   
! sigma value                   
  DOUBLE PRECISION sigma 
! actual observations                                                   
  DOUBLE PRECISION aobs,dobs 
! strings with asteroid names, output unit                              
  CHARACTER*80 titnam 
  CHARACTER*60 filnam 
  INTEGER, INTENT(IN) ::  iun20 
! station code, flag for use of act.obs., observation type              
  INTEGER, INTENT(IN) ::  ids,iff
  CHARACTER*(1) type
! target time (TDT, UTC)                                                
  DOUBLE PRECISION, INTENT(IN) :: t1,tut1 
! =================END INTERFACE====================                    
  DOUBLE PRECISION alm(mulx),dem(mulx)
!       ,adotm(mulx),ddotm(mulx),disv(mulx) 
  INTEGER ng,i,npop, icov, inl
  DOUBLE PRECISION eqm1(6,mulx),amagn(mulx),eq1(6),alpha,delta 
!      DOUBLE PRECISION amagn1,alpha1,delta1 ! for test comparison
  DOUBLE PRECISION :: adot,ddot ! proper motion
  DOUBLE PRECISION :: pha,dis,dsun,elo,gallat ! phase, distance to Earth, distance to Sun
  DOUBLE PRECISION :: gamad(2,2),axes(2,2),sig(2) ! covariance on sky plane
! ====================================================                  
  type='O' ! and does not handle radar observations
  inl=1 !
! first compute nominal prediction                                      
  CALL  predic_obs(elm(imi0),ids,t1,type,       &
     &        alpha,delta,amagn(imi0),inl,        &
     &        GAMAD=gamad,SIG=sig,AXES=axes,    &
     &        ADOT0=adot,DDOT0=ddot,DIS0=dis,DSUN0=dsun,ELO0=elo,PHA0=pha)
  WRITE(*,*)' nominal solution ' 
  icov=1
  CALL outobc(iun20,type,ids,tut1,alpha,delta,amagn(imi0),adot,ddot,&
     &     elo,dis,icov,gamad,sig,axes,DSUN=dsun,PHA=pha)                 
  alm(imi0)=0.d0 
  dem(imi0)=0.d0 
! initialize revolution counter                                         
  ng=0 
! loop on existing multiple solutions: first forward, then backward     
  DO  i=imi0+1,imip 
     dyn%dp=dpm(:,i)
     CALL  predic_obs(elm(i),ids,t1,type,alm(i),dem(i),amagn(i),inl,   &
     &        ADOT0=adot,DDOT0=ddot,DIS0=dis,DSUN0=dsun,ELO0=elo,PHA0=pha)
     IF(alm(i).gt.pig)alm(i)=alm(i)-dpig 
     WRITE(*,*)' alternate obs.no. ',i
     icov=1 
     CALL outobc(iun20,type,ids,tut1,alm(i),dem(i),amagn(i),        &
     &     adot,ddot,elo,dis,icov,gamad,sig,axes,DSUN=dsun,PHA=pha)
     alm(i)=alm(i)-alpha 
     dem(i)=dem(i)-delta 
! update revolution counter                                             
     CALL angupd(alm(i),alm(i-1),ng) 
  ENDDO
  DO  i=imi0-1,imim,-1
     dyn%dp=dpm(:,i)
     CALL  predic_obs(elm(i),ids,t1,type,alm(i),dem(i),amagn(i),inl,  &
     &        ADOT0=adot,DDOT0=ddot,DIS0=dis,DSUN0=dsun,ELO0=elo,PHA0=pha)
     WRITE(*,*)' alternate obs.no. ',i 
     icov=1
     CALL outobc(iun20,type,ids,tut1,alm(i),dem(i),amagn(i),        &
     &     adot,ddot,elo,dis,icov,gamad,sig,axes,DSUN=dsun,PHA=pha)           
     alm(i)=alm(i)-alpha 
     dem(i)=dem(i)-delta 
! keep count of lost revolutions in alpha                               
     CALL angupd(alm(i),alm(i+1),ng) 
  END DO
  dyn%dp=dpm(:,imi0)
! ===============================================                       
! output multiple prediction of observations                            
  CALL outmul(titnam,filnam,tut1,sigma,alpha,delta,                 &
     &        alm,dem,amagn,imim,imip,imi0,iff,aobs,dobs,type)            
END SUBROUTINE fmuobs 

!                                                                       
! ====================================================                  
! FMUPRO multiple state propagation for FITOBS                          
! ====================================================                  
SUBROUTINE fmupro(iun20,tr) 
  USE propag_state
! =================INPUT========================================= 
  INTEGER, INTENT(IN) :: iun20! output units
  DOUBLE PRECISION, INTENT(IN) ::tr ! target time
! ================END INTERFACE==========================               
! ======== output moid =====================                            
  INTEGER j,i ! loop indexes  
! VA number
  CHARACTER*9 astna0
! ===================================================================== 
! main loop                                                             
! propagation to time tr                                                
  DO j=imim,imip 
     astna0=' '
     IF(j.le.9)THEN              
        WRITE(astna0,111) j 
111     FORMAT('va_00000',I1)
     ELSEIF(j.le.99)THEN
        WRITE(astna0,112) j
112     FORMAT('va_0000',I2)
     ELSEIF(j.le.999)THEN
        WRITE(astna0,113) j
113     FORMAT('va_000',I3)
     ELSEIF(j.le.9999)THEN
        WRITE(astna0,114) j
114     FORMAT('va_00',I4)
     ELSEIF(j.le.99999)THEN
        WRITE(astna0,115) j
115     FORMAT('va_0',I5)
     ELSEIF(j.le.999999)THEN
        WRITE(astna0,116) j
116     FORMAT('va_',I6)
     ELSE
        WRITE(*,*)' nmulti: change format for imult>499999'
        STOP
     ENDIF
     WRITE(iun20,*)astna0
     WRITE(*,*)astna0
     WRITE(iuncla,*)astna0
     dyn%dp=dpm(:,j)
     CALL pro_ele(elm(j),tr,elm(j),unm(j),unm(j)) 
! orbital distance                                                      
     CALL nomoid(tr,elm(j),moid_m(j),dnp_m(j),dnm_m(j))
  ENDDO
  dyn%dp=dpm(:,imi0)
! ===================================================================== 
! summary table                                                         
! ===================================================================== 
  CALL tee(iun20,'SUMMARY OF MULTIPLE SOLUTIONS=') 
  WRITE(iun20,223) tr
  WRITE(*,223) tr 
223 FORMAT(' elements at time ',f8.1,' (MJD):') 
  IF(elm(1)%coo.eq.'EQU')THEN 
     CALL tee(iun20,'no       a      h      k      p      q      lambda=')
  ELSEIF(elm(1)%coo.eq.'CAR')THEN 
     CALL tee(iun20,'no       x      y      z    xdot    ydot    zdot=')
  ELSEIF(elm(1)%coo.eq.'KEP')THEN 
     CALL tee(iun20,'no       a      e      I    Omeg    omeg    mean.an=')
  ELSEIF(elm(1)%coo.eq.'COM')THEN
    CALL tee(iun20,'no        q      e      I    Omeg    omeg    t.peri=')
  ELSEIF(elm(1)%coo.eq.'ATT')THEN
    CALL tee(iun20,'no      alpha  delta   adot  ddot     r      rdot=')
  ENDIF
  DO i=imim,imip 
     WRITE(*,144)i,elm(i)%coord 
144  FORMAT(i5,6f12.8) 
     WRITE(iun20,144)i,elm(i)%coord
  ENDDO

  CALL tee(iun20,'no.,  magn,  MOID ,  nod+  ,  nod-=') 
  DO i=imim,imip 
     WRITE(*,145)    i,elm(i)%h_mag,moid_m(i),dnp_m(i),dnm_m(i) 
     WRITE(iun20,145)i,elm(i)%h_mag,moid_m(i),dnp_m(i),dnm_m(i)
145  FORMAT(i4,2x,f5.2,1x,f8.5,1x,f8.5,1x,f8.5,1x,i2) 
  ENDDO
END SUBROUTINE fmupro
!                                                                       
! =======================================                               
!  FMUPLO                                                               
! =======================================                               
! graphic output of multiple orbits                                     
! ===============INTERFACE========================                      
SUBROUTINE fmuplo(titnam,sigma) 
  USE util_suit
  CHARACTER*80, INTENT(IN) :: titnam 
  DOUBLE PRECISION, INTENT(IN) :: sigma
! ============END INTERFACE==================                           
  DOUBLE PRECISION a(mulx),e(mulx),aa,ee ,q,qg,enne
  INTEGER i,j,k, fail_flag
  TYPE(orbit_elem) :: elcom
  CHARACTER*3 coo
  INTEGER icoo,numb
  CHARACTER*60 xlab,ylab
  CHARACTER*20 menunam ! characters for menu
! ===========================================
  menunam='plotcoo'  
  CALL menu(icoo,menunam,3,' which coordinates to plot?=',            &
     &      'keplerian a,e=',                                          &
     &      'cometary q,e=',                                           &
     &      'attributable r, rdot=')   
  IF(icoo.eq.0)RETURN
  numb=imip-imim+1
  coo=elm(1)%coo
  IF(icoo.eq.1)THEN      
     xlab='Semimajor axis (au)'
     ylab='Eccentricity'             
! a-e plot of multiple solutions
     IF(coo.eq.'EQU')THEN
        DO i=1,numb 
           j=i+imim-1
           a(i)=elm(j)%coord(1) 
           e(i)=sqrt(elm(j)%coord(2)**2+elm(j)%coord(3)**2) 
        ENDDO
        aa=elm(imi0)%coord(1)
        ee=sqrt(elm(imi0)%coord(2)**2+elm(imi0)%coord(3)**2) 
     ELSEIF(coo.eq.'KEP')THEN
        DO i=1,numb 
           j=i+imim-1
           a(i)=elm(j)%coord(1) 
           e(i)=elm(j)%coord(2)
        ENDDO
        aa=elm(imi0)%coord(1)
        ee=elm(imi0)%coord(2)
     ELSE
        DO i=1,numb 
           j=i+imim-1
           CALL coo_cha(elm(j),'COM',elcom,fail_flag)
           IF(fail_flag.le.10.and.fail_flag.ne.4)THEN
              e(i)=elcom%coord(2)
              a(i)=abs(elcom%coord(1)/(1.d0-e(i)))
           ELSE
              WRITE(*,*)' fmuplo: conversion failure ', elcom%coo 
              RETURN
           ENDIF
        ENDDO
        CALL coo_cha(elm(imi0),'COM',elcom,fail_flag)
        IF(fail_flag.le.10.and.fail_flag.ne.4)THEN
           ee=elcom%coord(2)
           aa=abs(elcom%coord(1)/(1.d0-ee))
        ELSE
           WRITE(*,*)' fmuplo: conversion failure ', elcom%coo 
           RETURN
        ENDIF
     ENDIF
  ELSEIF(icoo.eq.2)THEN
     xlab='Perihelion distance (au)'
     ylab='Eccentricity'     
     DO i=1,numb 
        j=i+imim-1
        CALL coo_cha(elm(j),'COM',elcom,fail_flag)
        IF(fail_flag.le.10.and.fail_flag.ne.4)THEN
           e(i)=elcom%coord(2)
           a(i)=elcom%coord(1)
        ELSE
           WRITE(*,*)' fmuplo: conversion failure ', fail_flag
           RETURN
        ENDIF
     ENDDO
     CALL coo_cha(elm(imi0),'COM',elcom,fail_flag)
     aa=elcom%coord(1)
     ee=elcom%coord(2)
  ELSEIF(icoo.eq.3)THEN
     xlab='Geocentric distance (au)'
     ylab='Geocentric range rate (au/d)'     
     DO k=1,numb 
        j=k+imim-1
        CALL coo_cha(elm(j),'ATT',elcom,fail_flag)
        IF(fail_flag.le.10.and.fail_flag.ne.4)THEN
           e(k)=elcom%coord(6)
           a(k)=elcom%coord(5)
        ELSE
           WRITE(*,*)' fmuplo: conversion failure ', fail_flag
           RETURN
        ENDIF
     ENDDO
     CALL coo_cha(elm(imi0),'ATT',elcom,fail_flag)
     ee=elcom%coord(6)
     aa=elcom%coord(5)
  ENDIF
  CALL ploae(elm(imi0)%t,a,e,aa,ee,sigma,numb,titnam,xlab,ylab)
   
END SUBROUTINE fmuplo 
! ===================================================================== 
!  OUTMUL                                                               
! ===================================================================== 
! output multiple observations                                          
! =============INTERFACE=============================================== 
SUBROUTINE outmul(titnam,filnam,t1,sigma,alpha,delta,             &
     &     alm,dem,hmagn,imloc,iploc,i0loc,iff,aobs,dobs,type)
! =============INPUT=================================================== 
! file name                                                             
  CHARACTER*80,INTENT(IN) :: titnam 
  CHARACTER*60,INTENT(IN) ::  filnam 
! first, last and central index control for closed curve,obs type             
  INTEGER,INTENT(IN) ::  imloc,iploc,i0loc, iff
  CHARACTER*(1), INTENT(OUT) :: type
! observation time MJD, sigma value, nominal prediction                 
  DOUBLE PRECISION,INTENT(IN) :: t1,sigma,alpha,delta 
! observation: predicted  value alpha, delta, magnitude, actual         
  DOUBLE PRECISION,INTENT(IN) ::  alm(iploc),dem(iploc),hmagn(iploc),aobs,dobs 
! =============END INTERFACE=========================================== 
! conversion to sessagesimal                                            
  DOUBLE PRECISION seca,secd 
  INTEGER inta,mina,intd,mind 
  CHARACTER*1 siga,sigd 
! conversion of time                                                    
  INTEGER iy,imo,iday 
  DOUBLE PRECISION hour 
! scalar temporaries, differences                                       
  DOUBLE PRECISION dee,daa,ado,ddo 
! file name                                                             
  INTEGER le 
  CHARACTER*90 file, filnam1 
! loop indexes, units                                                   
  INTEGER n,iun7 
! ======================================================================
! open output file
  filnam1=filnam                                                      
  CALL rmsp(filnam1,le) 
  file=filnam1(1:le)//'.cbd' 
  CALL filopn(iun7,file,'unknown') 
! date and sigma value                                                  
  CALL mjddat(t1,iday,imo,iy,hour) 
  WRITE(iun7,297)iday,imo,iy,hour,sigma 
297 FORMAT(i3,i3,i5,f8.4,f5.2) 
! line of variations                                                    
  DO n=imloc,iploc 
     daa=alpha+alm(n) 
     daa=mod(daa,dpig) 
     IF(daa.lt.0.d0)daa=daa+dpig 
     IF(daa.gt.dpig)daa=daa-dpig 
     daa=daa*degrad/15 
     IF(daa.lt.0.d0.or.daa.gt.24.d0)THEN 
        WRITE(*,*)' outmul: daa out of range ', daa 
     ENDIF
     CALL sessag(daa,siga,inta,mina,seca) 
     dee=(delta+dem(n))*degrad 
     CALL sessag(dee,sigd,intd,mind,secd) 
! proper motion in arcsec/hour                                          
!        ado=adotm(n)*secrad/24.d0 
!        ddo=ddotm(n)*secrad/24.d0 
! output                                                                
     IF(siga.eq.'+')siga=' ' 
     WRITE(iun7,396)n,siga,inta,mina,seca,sigd,intd,mind,secd,hmagn(n)
!     &       disv(n),ado,ddo,hmagn(n)                                   
396  FORMAT(i3,1x,a1,i2,1x,i2,1x,f4.1,2x,a1,i2,1x,i2,1x,f4.1,1x,f5.2)
!     &       1x,f8.5,1x,f8.2,1x,f8.2,                           
  ENDDO
  CALL filclo(iun7,' ') 
! graphics output                                                       
  IF(iff.eq.1)THEN 
     CALL plocbd(titnam,alpha,delta,sigma,t1,                       &
          &         alm(imloc),dem(imloc),iploc-imloc+1,type)
  ELSEIF(iff.eq.2)THEN 
     CALL ploobs(titnam,alpha,delta,sigma,t1,                       &
     &         alm(imloc),dem(imloc),iploc-imloc+1, aobs,dobs)
  ENDIF
END SUBROUTINE outmul
! =========================================                             
! PROP_SIG                                                               
! propagator in sigma space                                             
SUBROUTINE prop_sig(batch,el1,el2,dn,sigma,mc,obs,obsw,wdir,sdir,units,fail,nd)
! ==============input================= 
  LOGICAL batch ! batch control
! current elements and epoch; stepsize factor, target sigma 
  TYPE(orbit_elem), INTENT(IN) :: el1            
  DOUBLE PRECISION, INTENT(IN) :: dn,sigma 
! ======observations====                                                
  INTEGER, INTENT(IN) :: mc ! number of obs
! new data types
  TYPE(ast_obs),DIMENSION(mc), INTENT(IN) :: obs
  TYPE(ast_wbsr),DIMENSION(mc), INTENT(IN) :: obsw
! weak direction, rms along it                                          
  DOUBLE PRECISION, DIMENSION(nd), INTENT(INOUT) :: wdir, units
  DOUBLE PRECISION, INTENT(INOUT) :: sdir
! ==============output================   
  TYPE(orbit_elem), INTENT(OUT) :: el2 ! elements
  LOGICAl, INTENT(OUT) :: fail ! failure flag
  INTEGER, INTENT(IN) :: nd ! number of parameters to be solved
!  ===========end interface============                                  
  DOUBLE PRECISION hh,ch,h,ridmax 
  INTEGER iun,imint 
! intermediate point: state, covariance, normal matr., RMS, norm corr   
  TYPE(orbit_elem) :: elt,eltp
  TYPE(orb_uncert) :: uncert 
  DOUBLE PRECISION :: csinor,delnor, wdir0(ndimx) 
! obs,residuals control, flag (for difcor)         
  INTEGER inew, icor(ndimx),ncor,itmaxold,j 
! success control is dummy                                              
  LOGICAL succ 
! use weak direction to find initial conditions                         
  wdir0(1:nd)=wdir(1:nd)
  hh=dn*sigma 
  IF(.not.batch.and.verb_mul.ge.9)THEN 
     iun=abs(iun_log) 
  ELSE 
     iun=-abs(iun_log) 
  ENDIF
! 1=Euler 2= RK2                                                        
  imint=2 
! iteration on stepsize reductions                                      
!     ridmax=30.d0                                                      
  ridmax=130.d0 
  h=hh 
  ch=0.d0 
  elt=el1 
! try to do the step at once...                                         
1 CONTINUE 
!  write(*,*)'prop_sig 0: h,a, sdir',h,elt%coord(1),sdir 
!  WRITE(*,*)'mmulti: ch,h, eqt ', ch,h,(eqt(j),j=1,3) 
!  WRITE(*,*)' prop-sig: about to enter int_step',tc, obs(mc)%time_tdt
  CALL int_step(elt,el2,h,imint,mc,obs,obsw,iun,wdir,sdir,units,fail,nd) 
!  write(*,*)'prop_sig 1: a, sdir',elt%coord(1),sdir 
! if failed, half step                                                  
  IF(fail)THEN 
     IF(abs(h).ge.abs(hh)/ridmax)THEN 
        h=h/2.d0 
        IF(verb_mul.ge.20.and.iun.gt.0)THEN 
           WRITE(iun,*)'prop_sig: halving', h,hh 
        ENDIF
        GOTO 1 
     ELSE 
        IF(verb_mul.ge.5.and.iun.gt.0)THEN 
           WRITE(iun,*)'prop_sig: too many halvings ',h,hh,ridmax 
        ENDIF
        WRITE(*,*)'prop_sig: too many halvings ',h,hh,ridmax
        RETURN 
     ENDIF
  ELSE 
     ch=ch+h 
  ENDIF
  IF(abs(ch).lt.abs(hh))THEN 
     elt=el2 
! compute line of variations                                            
! compute vectorfield at intermediate point                             
     itmaxold=itmax 
     itmax=0 
     CALL whicor(0,nd,icor,ncor,inew) 
! compute covariance and residual norm     
!  WRITE(*,*)' prop-sig: about to enter diff-cor',tc, obs(mc)%time_tdt 
     CALL diff_cor(mc,obs,obsw,elt,icor,iun,eltp,uncert,csinor,delnor,succ,nd)
     itmax=itmaxold 
! compute weak direction and length                                     
     CALL weak_dir(uncert%g(1:nd,1:nd),wdir,sdir,iun,eltp%coo,eltp%coord,units,nd) 
!     write(*,*)'prop_sig 2: a, sdir',elt%coord(1),sdir 
     IF(DOT_PRODUCT(wdir,wdir0).lt.0.d0)wdir=-wdir
     GOTO 1 
  ENDIF
  IF(iun.gt.0.and.verb_mul.gt.20)WRITE(iun,*)' prop_sig: stepsize ',h 
END SUBROUTINE prop_sig
!=================================================                      
! INT_STEP integration step for propagation along LOV                    
SUBROUTINE int_step(el0,el1,hh,imint,mc,obs,obsw,iun,wdir,sdir,units,fail,nd)
! stepsize, integration method                                          
  DOUBLE PRECISION, INTENT(IN) ::  hh 
  INTEGER, INTENT(IN) :: imint 
  INTEGER,INTENT(IN) :: nd ! number of solve for params
! weak direction (in input computed at eq0)             
  DOUBLE PRECISION, INTENT(INOUT) :: wdir(nd),units(nd),sdir
! ======observations====                                                
! number of obs.
  INTEGER, INTENT(IN) ::  mc
 ! new data types
  TYPE(ast_obs),DIMENSION(mc), INTENT(IN) :: obs
  TYPE(ast_wbsr),DIMENSION(mc), INTENT(IN) :: obsw
! initial conditions: epoch, elements, new elements after step 
  TYPE(orbit_elem), INTENT(IN) :: el0
  TYPE(orbit_elem), INTENT(OUT) :: el1
! output unit                                                           
  INTEGER iun 
! failure flag                                                          
  LOGICAl fail
! =========end interface================                                
! intermediate point: state, covariance, normal matr., RMS, norm corr   
  TYPE(orbit_elem) ::  el12,el12p
  TYPE(orb_uncert) :: un12
  DOUBLE PRECISION :: csino12,delno12 
! weak direction (in input computed at eq0), length of axis             
  DOUBLE PRECISION wdir12(nd),sdir12,wdirst(nd),sdirst 
! obs,residuals control, flag (for difcor)         
  INTEGER inew, icor(nd),ncor,itmaxold 
! success control is dummy,bizarre is impor                             
  LOGICAL succ,bizarre 
! loop indexes                                                          
  INTEGER j 
! controls of fixed point iteration                                     
  INTEGER itx,it 
  DOUBLE PRECISION eps,cosa,ang,dsi,e1,e12,direc,q,qg,enne 
! intermediate and initial dynamical parametera
  DOUBLE PRECISION,DIMENSION(ndyx) :: dp12,dp0 
! ================================================                      
! initial dynamical parameters
  DO j=1,nls
     dp0(ls(j))=dyn%dp(ls(j))
  ENDDO
  
  IF(imint.eq.1)THEN 
! Euler method
     el1=el0                                                          
     el1%coord=el0%coord+units(1:6)*wdir(1:6)*sdir*hh
! update dynamical parameters
     DO j=1,nls
        dyn%dp(ls(j))=dp0(ls(j))+units(6+j)*wdir(6+j)*sdir*hh
     ENDDO
     fail=.false. 
  ELSEIF(imint.eq.2)THEN 
! Runge-Kutta-Gauss of order 2:                                         
     itx=5 
     eps=1.d-4 
! first store vectorfield at el0                            
     wdirst=wdir
     sdirst=sdir 
! compute intermediate point (first guess)                              
     el12=el0 
     el12%coord=el0%coord+units(1:6)*wdir(1:6)*sdir*hh*0.5d0
! intermediate dynamical parameters
     DO j=1,nls
        dp12(ls(j))=dp0(ls(j))+units(6+j)*wdir(6+j)*sdir*hh*0.5d0
     ENDDO
! before running differential corrections, tests that orbit is still good
     IF(bizarre(el12,e12))THEN 
        fail=.true. 
        IF(iun.gt.0)WRITE(iun,*)'int_step: bizarre ',e12 
        el1=el12
        DO j=1,nls
           dyn%dp(ls(j))=dp12(ls(j))
        ENDDO
        RETURN 
     ENDIF
! fixed point iteration                                                 
     DO it=1,itx 
! compute vectorfield at intermediate point                             
        itmaxold=itmax 
        itmax=0 
        CALL whicor(0,nd,icor,ncor,inew) 
! compute covariance and residual norm                                  
!       write(*,*)mc,iobs(1),ioco(1)      
! dynamical parameters are not passed, they must be set
        DO j=1,nls
           dyn%dp(ls(j))=dp12(ls(j))
        ENDDO
        CALL diff_cor(mc,obs,obsw,el12,icor,iun,el12p,un12,csino12,delno12,succ,nd)
        itmax=itmaxold 
! compute weak direction and length                                     
        CALL weak_dir(un12%g(1:nd,1:nd),wdir12,sdir12,iun,el12p%coo,el12p%coord,units,nd) 
        direc=DOT_PRODUCT(wdir12,wdir)  
        IF(direc.lt.0.d0) wdir12=-wdir12
! store vectorfield at current iteration                               
        wdirst=wdir12
        sdirst=sdir12 
! recompute intermediate point  
        el12%coord=el0%coord+units(1:6)*wdir12(1:6)*sdir12*hh*0.5d0
! non grav
        DO j=1,nls
           dp12(ls(j))=dp0(ls(j))+units(6+j)*wdir12(6+j)*sdir12*hh*0.5d0
        ENDDO
!    write(*,*)'weak sigma, iter, eq12',sdir12,it,(eq12(j),j=1,3) 
! before a new iteration, tests that orbit is still elliptic            
        IF(bizarre(el12,e12))THEN 
           IF(iun.gt.0)WRITE(iun,*)'int_step: bizarre ',e12 
           fail=.true. 
           el1=el12
           DO j=1,nls
              dyn%dp(ls(j))=dp12(ls(j))
           ENDDO
           RETURN 
        ENDIF
! convergence control                                                   
        cosa=DOT_PRODUCT(wdir12,wdirst) ! cosa=prscag(6,wdir12,wdirst) 
        IF(cosa.le.1.d0)THEN 
           ang=acos(cosa)
        ELSE 
           ang=0.d0 
        ENDIF
        dsi=abs(sdirst-sdir12)/sdirst 
!    WRITE(*,*)' int_step: it, ang, dsi, sdir12,sdirst ',          
!    +          it,ang,dsi,sdir12,sdirst                                
        IF(ang.lt.eps.and.dsi.lt.eps)THEN 
!             WRITE(*,*)' int_step: it, ang, dsir ',it,ang,dsi           
! convergence achieved                                                  
           GOTO 2 
        ENDIF
        wdirst=wdir12
        sdirst=sdir12 
     ENDDO
! convergence failed                                                    
! WARNING: very verbose
     WRITE(*,*)' int_step: failed convergence, it, ang, dsig'         
     WRITE(*,*)it-1,ang,dsi 
! WARNING                      
     fail=.true. 
     IF(iun.gt.0)WRITE(iun,*)'int_step: non convergent ',el12 
     el1=el12 
     DO j=1,nls
        dyn%dp(ls(j))=dp12(ls(j))
     ENDDO
     RETURN 
! convergence succeeded                                                 
2    CONTINUE 
! compute final point 
     el1=el0                                                          
     el1%coord=el0%coord+units(1:6)*wdir12(1:6)*sdir12*hh
! final non grav
     DO j=1,nls
        dyn%dp(ls(j))=dp0(ls(j))+units(6+j)*wdir12(6+j)*sdir12*hh
     ENDDO
! tests that orbit is still acceptable 
     IF(bizarre(el1,e1))THEN 
        fail=.true. 
        IF(iun.gt.0)WRITE(iun,*)'int_step: bizarre ',e1 
     ELSE 
        fail=.false. 
     ENDIF
! but we go on                                                          
  ELSE 
     WRITE(*,*)' int_step: not invented yet, imint=',imint 
     STOP 
  ENDIF
!     write(*,*)'int_step: at exit ',fail, (eq1(j),j=1,3)                
END SUBROUTINE int_step
                    
! ==================================================================
!  lovinit 
! initializes storage of multiple solution data                 
! ======================================================
SUBROUTINE lovinit(astname,mulsodir,obsdir,progna,iunout,      &
       &     succ,vel_inf,imult,sigma_max,sigma_step_max,no_outcov)
  USE obssto
  USE least_squares  
  USE dyn_param                     
! ==================INPUT=======================================================================  
  CHARACTER(LEN=9),  INTENT(IN)           :: astname          ! asteroid name 
  CHARACTER(LEN=80), INTENT(IN)           :: obsdir,mulsodir  ! obs and mult. sol. 
                                                              ! directory  name 
  INTEGER,           INTENT(IN)           :: iunout           ! output file 
  CHARACTER(LEN=*),  INTENT(IN)           :: progna           ! program name
  LOGICAL,           INTENT(IN), OPTIONAL :: no_outcov        ! accept orbits without covariance
!==================OUTPUT=======================================================================
  LOGICAL,          INTENT(OUT) :: succ             ! success flag 
  DOUBLE PRECISION, INTENT(OUT) :: vel_inf          ! velocity at infinite 
                                                    ! with respect to Earth (circular approx)
  INTEGER,          INTENT(OUT) :: imult            ! multiple solution specifications       
  DOUBLE PRECISION, INTENT(OUT) :: sigma_max        ! Maximum sigma value 
  DOUBLE PRECISION, INTENT(OUT) :: sigma_step_max   ! Maximum step in sigma
                                                    ! (If prob_sampl=false, sigma_step_max is equal to delta_sigma)
! ===============OUTPUT THROUGH MODULE MULTIPLE_SOL=============================================
! time of initial conditions 
! index range of multiple solution, of nominal solution  
! =============for call to read_elems=============================
  CHARACTER(LEN=160) :: catname,repname,file   ! file names 
  INTEGER            :: iuncat, iunrep         ! unity 
  CHARACTER(LEN=19) :: name0 
  CHARACTER(LEN=9) :: name1 
  INTEGER :: le, imtmp, norb
  DOUBLE PRECISION hmagn,tcat
  TYPE(orbit_elem) :: eltmp
  TYPE(orb_uncert) :: unctmp
  INTEGER ndp,nmod, nlsloc,lsloc(1:ndyx)   ! LSP record
  CHARACTER*4 form ! 1l, ML
  LOGICAL nooutcov
! non-grav parameters
  DOUBLE PRECISION :: dptmp(ndyx)
  INTEGER :: nd  
  INTEGER :: j                            ! VA number 
  LOGICAL :: eof, err                     ! end of file, error in reading  
  DOUBLE PRECISION :: v_infty,v_max,v_min ! velocity w.r.to Earth max and min
  CHARACTER(LEN=60) :: rwofi0             ! residuals and weights file name   
  LOGICAL :: obs0!,precob,change           ! logical flag et al
  CHARACTER(LEN=20) :: error_model ! error model file 
  INTEGER :: iun20                        ! unit 
  INTEGER :: i,ii                       ! loop indexes  
  INCLUDE 'sysdep.h90'                    ! system dependencies      
!===========================================================
  v_min=100.d0 
  v_max=0.d0 
  imul=0 ! flag to find if the elements have been read
! opening and reading multiple solution catalog                         
  succ=.false. 
  catname=mulsodir//dircha//astname//'.ctc' 
  CALL rmsp(catname,le)
  CALL filopn(iuncat,catname(1:le),'old')
  file=' '
  CALL oporbf(file,iuncat) 
! ================ INPUT MULTIPLE SOLUTIONS (REP FILE) ===========================
! multiple solution report file  
  repname=mulsodir//dircha//astname//'.mrep'                                   
  CALL rmsp(repname,le) 
  CALL filopn(iunrep,repname(1:le),'old')
!Reading mrep: prob_sampl
  READ(iunrep,5)prob_sampl
5 FORMAT(43X,L7)
  WRITE(*,55) prob_sampl
55 FORMAT(' prob_sampl  = ', L1)
  IF(prob_sampl)THEN
     READ(iunrep,399)prob_step, sigma_max, sigma_step_max, imult 
399 FORMAT(38X,E7.1,11X,F5.3,16X,F5.3,7X,I5)
     WRITE(iunout,400)prob_step,sigma_max, sigma_step_max, imult
400  FORMAT(' prob_step      = ',1P,E7.1,0P/,' sigma max      = ',F5.3/,' sigma step max = ',F5.3/,' imult         = ',I5) 
  ELSE
     READ(iunrep,199)imult,sigma_max
199  FORMAT(34X,I5,10X,F5.3) 
     WRITE(iunout,200)imult, sigma_max
200  FORMAT(' imult       =',I5/,' sigma max   = ',F5.3) 
     sigma_step_max=sigma_max/imult
  END IF
! Reading scaling LOV and second LOV
  READ(iunrep,299)scaling_lov,second_lov
299 FORMAT(13X,L1/13X,L1)
  WRITE(iunout,499)scaling_lov,second_lov
499 FORMAT(' scaling_lov = ',L1/,' second_lov  = ',L1)
  ! Skipping header 
  READ(iunrep,*)
  imi0=imult+1 
! option for not reading covariances
  IF(PRESENT(no_outcov))THEN
     nooutcov=no_outcov
  ELSE
     nooutcov=.false.
  ENDIF
! read orbit one by one                                                 
  DO i=1,mulx 
     CALL read_elems(eltmp,name0,eof,nlsloc,err,lsloc,form,unctmp,UNIT=iuncat)
     IF(eof) GOTO 2
! Check error
     IF(.not.err)THEN
! Save dyn%dp, nls, and ls
        IF(dyn%ndp.gt.0)THEN
           nls=nlsloc
           ls=lsloc
        ENDIF
! Dimension
        nd=6+nls
     ELSE
! failure case: either .eq0 does not exist, or it is rotten  
        WRITE(iun_log,*)'!!WARNING! WARNING! WARNING! WARNING!!' 
        WRITE(iun_log,*)'Asteroid ',astname,          &
             &           ' initial conditions not found.'                       
        STOP
     ENDIF
! check consistency of epoch times
     IF(i.eq.1)THEN 
        tcat=eltmp%t 
     ELSE 
        IF(eltmp%t.ne.tcat)THEN 
           WRITE(*,*)'lovinit: time discrepancy from tcat=',tcat 
           WRITE(*,*)'lovinit: at record ',i,' t=',eltmp%t 
           STOP 
        ENDIF
     ENDIF

! handling of name; but note that the multiple solutions are assumed to 
! and with consecutive indexes!!! that is imul(i)=imul(1)+i-1           
     CALL splinam(name0,name1,j) 
     IF(j.gt.0.and.j.le.mulx)THEN
        imul(i)=j
     ELSEIF(j.lt.0)THEN
        WRITE(*,*)'lovinit: VA with negative index', j
        STOP
     ELSE
        WRITE(*,*)'lovinit: VA with index too large=',j,' increase mulx=',mulx
        STOP
     ENDIF
! availability of covariance                                            
     IF(.NOT.(unctmp%succ))THEN 
        IF(nooutcov.and.j.ne.imi0)THEN
! ok not to have covariance
        ELSE         
           WRITE(*,*)'lovinit:',  name0, ' matrices not avalaible' 
           STOP
        ENDIF
     ENDIF

! note that the multiple solutions are assumed to be written
! with consecutive indexes!!! 
     IF(i.gt.1)THEN 
        IF(imul(i-1).eq.0)THEN
           WRITE(*,*)'lovinit: indexes not in order: VA ',j,' at record ',i,' missing ',j-1
        ENDIF
     ENDIF
! ok, copy in multiple_sol.mod
     elm(j)=eltmp
     unm(j)=unctmp
     IF(dyn%ndp.gt.0) dpm(1:ndyx,j)=dyn%dp
     READ(iunrep,145)ii,csinom(j),delnom(j),hmagn,                   &
          &  moid_m(j),dnp_m(j),dnm_m(j), sigmavalv(j), deltasigmav(j), sigq(j)
145  FORMAT(I5,1X,1P,E13.5,E11.3,2X,0P,F5.2,1X,                     &
          &  F8.5,1X,F10.5,1X,F10.5,3X,F12.8,4X,F12.8,3X,F12.8)
     IF(ii.ne.j)WRITE(*,*)'lovinit: discrepancy between mrep and ctc files'
! v_infinity computation                                                
     v_inf(j)=v_infty(elm(j)) 
     v_min=min(v_min,v_inf(j)) 
     v_max=max(v_max,v_inf(j)) 
! nominal solution                                                      
!     IF(j.eq.imi0)csinor=csinom(j) 
  ENDDO
  CALL filclo(iunrep,' ')
! increase mulx                                                        
  WRITE(*,*)'lovinit: increase nmax, file is longer than',mulx 
2   CONTINUE 
  CALL clorbf 
  norb=i-1 
  imim=imul(1) 
  imip=imul(norb) 
  WRITE(iunout,*)' lovinit: input of ',norb,' multiple solutions' 
  WRITE(iunout,*)' lovinit: max of v_infty ',v_max,' minimum ',v_min 
  vel_inf=v_min 
  tdt_cat=tcat
! input observation data                                                
!    precob=.false. 
  iun20=iunout 
! to get error_model from rwo file; does not work with .obs file
  CALL retinobs(obsdir,astname,obs0,error_model,rms_m,rmsmag_m)
  CALL errmod_set(error_model)
! find if there are radar data                                          
  CALL radar_ob(obs_m(1:m_m)%type,m_m) 
! forces the error_model, but can use .obs file
  IF(m_m.gt.0.and.obs0)succ=.true. 
END SUBROUTINE lovinit
! =====================================================                  
! LOVMAGN                                                               
! provides magnitude and v_inty for the VA nearest to the               
! real index x; also velocity at infinity                       
! =====================================================                  
SUBROUTINE lovmagn(x,v_i,h) 
! ================INPUT===============================
  DOUBLE PRECISION, INTENT(IN) :: x       ! fractional index 
! ================OUTPUT==============================           
  DOUBLE PRECISION, INTENT(OUT) :: v_i    ! velocity at 
                                            ! infinity (U in au/day)
  DOUBLE PRECISION, INTENT(OUT) :: h      ! H magnitude 
! =============END INTERFACE===========================                 
  INTEGER :: j                            ! loops index
! ===================================================
  DO j=imim,imip 
     IF(x.lt.j)THEN 
        v_i=v_inf(j) 
        h=elm(j)%h_mag 
        GOTO 2 
     ENDIF
  ENDDO
  h=elm(imip)%h_mag 
  v_i=v_inf(imip) 
2 CONTINUE 
END SUBROUTINE lovmagn

!====================================================
! LOVDENSIFY
! adds to elm_loc, dpm_loc the  intermediate VA
! ===================================================

SUBROUTINE lovdensify(ndiv,lre,ire,j1,j2,linear,dsig,newr,nadd)
  USE cla_store
  USE tp_trace, ONLY: vvv_tp
! INTERFACE
  INTEGER, INTENT(IN) :: ndiv                ! number of substeps for each step
  INTEGER, INTENT(IN) :: lre,ire                 ! length of array vas_traceloc
  INTEGER, INTENT(IN) :: j1,j2               ! range of indexes in array vas_traceloc  
  LOGICAL, INTENT(IN) :: linear              ! if true, use linear LOV
  DOUBLE PRECISION, INTENT(IN) :: dsig       ! current stepsize in sigma
  DOUBLE PRECISION, INTENT(OUT) :: newr(noxloc) ! new rindex
  INTEGER, INTENT(OUT) :: nadd               ! number of new vas_traceloc added
                                             ! lre(IN)+nadd=lre(OUT)
! END INTERFACE
  INTEGER j,i,k,jj,iii                       ! loop indexes
  INTEGER lrenew                             ! dimension for array elm
  DOUBLE PRECISION dsigma,dr                 ! current stepsize in sigma, in real index
 
!  TYPE(orbit_elem) eltmp                    ! temporary orbit elements
  DOUBLE PRECISION wdir0(ndimx),sdir0        ! weak direction, rms along it
  DOUBLE PRECISION vvv_used(ndimx)           ! increment to be actually used
  DOUBLE PRECISION,DIMENSION(ndimx) :: units ! for scaling
  DOUBLE PRECISION dnn,delta_sigma_new       ! new step in index, in sigma 
! for non-grav parameters
  INTEGER nd
  DOUBLE PRECISION dp0(ndyx),dp(ndyx)        ! non-grav parameters
  LOGICAl bizarre                            ! control of bizarre orbits
! paranoia checks
  IF(.not.lin_lov)THEN
     WRITE(*,*)' lovdensify: not ready for lin_lov=',lin_lov
     STOP 'lin_lov missing'
  ENDIF
  IF(j1.gt.j2.or.j1.lt.0.or.j2.gt.lre+1)THEN
     WRITE(*,*)' lovdensify: wrong indexes in call; lre=',lre,' j1=',j1,' j2=',j2
     STOP
  ENDIF 
  IF(.not.linear)THEN
     WRITE(*,*)' lovdensify: nonlinear LOV not ready'
     STOP
  ENDIF
! dimension of state vector being modified along the LOV
  nd=6+nls
! new stepsize (in sigma space)
  delta_sigma_new=dsig/float(ndiv)
  dr=delta_sigma_new/delta_sigma
! weak direction
  IF(scatterplane)THEN
     vvv_used(1:nd)=vvv_tp(1:nd)
  ELSE
     CALL weak_dir(unm(imi0)%g(1:nd,1:nd),wdir0,sdir0,0,elm(imi0)%coo,elm(imi0)%coord,units,nd)
     vvv_used(1:nd)=units(1:nd)*wdir0(1:nd)*sdir0
  ENDIF
  IF(nd.gt.6)THEN
     DO k=1,nls
        dp0(k)=dpm(ls(k),imi0)
     ENDDO
  ENDIF
! loop on intervals of the previous subdivision
  DO j = j1,j2-1
! check that the region to be densified has a constant step
     dsigma=vas_traceloc(j+1)%sigma-vas_traceloc(j)%sigma
     IF(ABS(dsig-dsigma).gt.10.d0*epsilon(1.d0))THEN
        WRITE(*,*)'lovdensify: impossible to densify if step is not constant'
        WRITE(*,*)' j1, j2=',j1,' ',j2,' dsig= ',dsig,' j=',j,' dsigma=',dsigma
        STOP
     ENDIF
  ENDDO
! setup the local array of elements, dynpar
!  iii=vas_traceloc(1)%rindex
!  elm_loc(1)=elm(iii)
!  elm_loc(1:lre)=elm(ire:ire+lre-1)
!  IF(nd.gt.6)THEN
!     dpm_loc(:,1:lre)=dpm(:,ire:ire+lre-1)
!  ENDIF
! counters
  nadd=0
  lrenew=lre
  newr(1:lre)=vas_traceloc(1:lre)%rindex
! loop on intervals of the previous subdivision
  DO j = j1,j2
! loop on new points
     IF(j.eq.j1)THEN 
        DO i=1,ndiv-1
! find the values of rindex for which a new point has to be added
           newr(lrenew+i)=vas_traceloc(j1)%rindex+dr*(i-ndiv)
!           elm_loc(lrenew+i)=elm_loc(j)
!           elm_loc(lrenew+i)%coord=elm_loc(lrenew+i)%coord-vvv_used(1:6)*dr*(ndiv-i)
!           IF(nd.gt.6)THEN
! find new non-grav parameters
!              dpm_loc(:,lrenew+i)=dpm_loc(:,j)
!              DO k=1,nls
!                 dp(k)=dp0(k)-vvv_used(6+k)*dr*(ndiv-i)
!                 dpm_loc(ls(k),lrenew+i)=dp(k)
!              ENDDO
!           ENDIF
! find new orbital elements
        ENDDO
        lrenew=lrenew+ndiv-1
        nadd=nadd+ndiv-1
     ENDIF
     DO i=1,ndiv-1
! find the values of rindex for which a new point has to be added
        newr(lrenew+i)=vas_traceloc(j)%rindex+dr*i
! find new orbital elements
!        elm_loc(lrenew+i)=elm_loc(j-1)
!        elm_loc(lrenew+i)%coord=elm_loc(lrenew+i)%coord+vvv_used(1:6)*dr*i
!        IF(nd.gt.6)THEN
! find new non-grav parameters
!           dpm_loc(:,lrenew+i)=dpm_loc(:,j-1)
!           DO k=1,nls
!              dp(k)=dp0(k)+vvv_used(6+k)*dr*i
!              dpm_loc(ls(k),lrenew+i)=dp(k)
!           ENDDO
!        ENDIF
     ENDDO
     lrenew=lrenew+ndiv-1
     nadd=nadd+ndiv-1
  ENDDO

END SUBROUTINE lovdensify

! ===================================================                    
! LOVINTERP                                                             
! provides an interpolated orbit along the LOV                          
! (with non-integer index rindex)                                       
! ===================================================   
                 
SUBROUTINE lovinterp(rindex,deltasig,el0,unc0,succ)
  USE obssto
  USE least_squares
  USE tp_trace, ONLY: vvv_tp
! ====================INPUT=============================      
  DOUBLE PRECISION, INTENT(IN) :: deltasig      ! interval in sigma 
                                                  ! between solutions
  DOUBLE PRECISION, INTENT(IN) ::  rindex       ! real index
! ====================OUTPUT============================ 
  TYPE(orbit_elem), INTENT(OUT) :: el0          ! elements     
  TYPE(orb_uncert), INTENT(OUT) :: unc0  ! normal and covariance matrices
                                                  ! corresponding to rindex
  LOGICAL, INTENT(OUT) :: succ                  ! success flag  
! ==================END INTERFACE=============================
  LOGICAL :: batch                              ! batch control
! ================LOV INTERPOLATION===========================
  DOUBLE PRECISION :: wdir(ndimx),units(ndimx),sdir,h0,diff(ndimx) 
  LOGICAL :: fail, bizarre 
!  INTEGER :: iun20 
  DOUBLE PRECISION :: csinew,delnew,rmshnew    ! for constr_fit   
  INTEGER :: nused 
  TYPE(orbit_elem) :: elc                ! corrected 
  INTEGER nd ! no. solve for parameters
! =============LOOP INDEXES, LOV INDEXES=======================
  INTEGER :: j,i,imu,ir 
  DOUBLE PRECISION :: s, ecc 
! =============================================================
! failure case ready                                                    
  succ=.false. 
! interpolate elements                                                  
  ir=NINT(rindex)

! extrapolation?                                                        
! WARNING: will result in failure if extrapolation is by a              
! long step beyond the end of the LOV string already computed           
  IF(ir.lt.imim)THEN 
     WRITE(*,*)'lovinterp: extrapolation rindex=',rindex,           &
            &        ' beyond ',imim                                        
     imu=imim 
  ELSEIF(ir.gt.imip)THEN 
     WRITE(*,*)'lovinterp: extrapolation rindex=',rindex,           &
            &        ' beyond ',imip                                     
     imu=imip 
  ELSE 
     imu=ir 
  ENDIF
! linear interpolation between the two consecutive ones                 
  s=rindex-imu 
!  iun20=-1  
  nd=6+nls
! non-grav parameters
! WARNING: use convex combination is easier
  IF(lin_lov.and..not.scatterplane)THEN  
     CALL weak_dir(unm(imu)%g(1:nd,1:nd),wdir,sdir,-1,elm(imu)%coo,elm(imu)%coord,units,nd) 
! linear shift in direction wdir
     el0=elm(imu)
     el0%coord=elm(imu)%coord+units(1:6)*wdir(1:6)*sdir*deltasig*s
     IF(nd.gt.6)THEN
        DO j=1,nls
           dyn%dp(ls(j))=dpm(ls(j),imu)+units(6+j)*wdir(6+j)*sdir*deltasig*s
        ENDDO
     ENDIF
     unc0=unm(imu)
     succ=.true.
  ELSEIF(lin_lov.and.scatterplane)THEN
 ! linear shift in direction vvv
     el0=elm(imu)
     el0%coord=elm(imu)%coord+vvv_tp(1:6)*delta_sigma*s
 !    el0%coord=(1-s)*elm(imu)%coord+s*elm(imu+1)%coord
     IF(nd.gt.6)THEN
        DO j=1,nls
           dyn%dp(ls(j))=dpm(ls(j),imu)+vvv_tp(6+j)*delta_sigma*s
 !          dyn%dp(ls(j))=(1-s)*dpm(ls(j),imu)+s*dpm(ls(j),imu+1)
        ENDDO
     ENDIF
     unc0=unm(imu)
     succ=.true.
  ELSEIF(.not.lin_lov)THEN
     CALL weak_dir(unm(imu)%g(1:nd,1:nd),wdir,sdir,-1,elm(imu)%coo,elm(imu)%coord,units,nd) 
! anti reversal check
     IF(imu.lt.imip)THEN
        diff(1:6)=elm(imu+1)%coord-elm(imu)%coord
        DO j=7,nd
           diff(j)=dpm(ls(j-6),imu+1)-dpm(ls(j-6),imu)
        ENDDO
     ELSE
        diff(1:6)=elm(imu)%coord-elm(imu-1)%coord
        DO j=7,nd
           diff(j)=dpm(ls(j-6),imu)-dpm(ls(j-6),imu-1)
        ENDDO
     ENDIF
     DO j=1,nls
        dyn%dp(ls(j))=dpm(ls(j),imu)
     END DO
     IF(DOT_PRODUCT(diff(1:nd),wdir(1:nd)).lt.0.d0)wdir(1:nd)=-wdir(1:nd)
! nonlinear interpolation   
     CALL prop_sig(batch,elm(imu),el0,s,deltasig,m_m,obs_m,obsw_m,wdir,sdir,units,fail,nd)
! check for hyperbolic                                                  
     IF(fail)THEN 
        WRITE(*,*)'step ',rindex,' hyperbolic' 
        WRITE(*,*)el0%coo,el0%coord 
        RETURN 
     ELSEIF(bizarre(el0,ecc))THEN 
        WRITE(*,*)'step ',rindex,' byzarre' 
        WRITE(*,*) el0%coo,el0%coord, ecc
        RETURN 
     ENDIF
! constrained corrections: 
     CALL constr_fit(m_m,obs_m,obsw_m,el0,wdir,elc,unc0,csinew,delnew,rmshnew,nused,succ,nd)
! exit if not convergent                                                
     IF(.not.succ) THEN 
        WRITE(*,*)'lovinterp: constr_fit failed for ',rindex 
        RETURN 
     ELSE
        el0=elc 
        succ=.true. 
     ENDIF
  ELSE
     WRITE(*,*)'lovinterp: case not yet invented'
     STOP
  ENDIF
END SUBROUTINE lovinterp

!====================================================!
! LOVINTERP3                                         !
!----------------------------------------------------!                   
! It provides an interpolated orbit along the LOV    !
! (with non-integer index rindex)                    !                  
!====================================================!
SUBROUTINE lovinterp3(rindex,sigmaout,el0,unc0,succ)
  USE obssto
  USE least_squares
  USE tp_trace, ONLY: vvv_tp
  DOUBLE PRECISION, INTENT(IN)  :: rindex   ! Real index
  DOUBLE PRECISION, INTENT(OUT) :: sigmaout ! Sigma and index values at end of step 
  TYPE(orbit_elem), INTENT(OUT) :: el0      ! Elements     
  TYPE(orb_uncert), INTENT(OUT) :: unc0     ! Normal and covariance matrices corresponding to rindex
  LOGICAL,          INTENT(OUT) :: succ     ! Success flag  
  ! Local variables
  LOGICAL          :: batch        ! Batch control
  DOUBLE PRECISION :: deltasig     ! Length in sigma_lov of the interval containing the index
  DOUBLE PRECISION :: wdir(ndimx)  ! Weak direction
  DOUBLE PRECISION :: units(ndimx) ! Units for scaling
  DOUBLE PRECISION :: sdir         ! Eigenvalue of the weak direction
  DOUBLE PRECISION :: diff(ndimx)  ! Element difference over the considered interval
  LOGICAL          :: fail         ! Failure in prop_sig
  LOGICAL          :: bizarre      ! Bizarre flag
  DOUBLE PRECISION :: csinew       ! Residuals norm after constr_fit
  DOUBLE PRECISION :: delnew       ! Last corr. norm after constr fit
  DOUBLE PRECISION :: rmshnew      ! RMS of H after constr_fit
  INTEGER          :: nused        ! No. of observations used by constr_fit
  TYPE(orbit_elem) :: elc          ! Corrected elements 
  INTEGER          :: nd           ! No. solve for parameters
  INTEGER          :: j,i          ! Indexes
  INTEGER          :: imu          ! Nearest LOV index to rindex
  INTEGER          :: ir           ! Nearest integer to rindex
  DOUBLE PRECISION :: s            ! Fractional part of rindex
  DOUBLE PRECISION :: ecc          ! Eccentricity for bizarre control

  IF(.NOT.prob_sampl)THEN
     WRITE(*,*) ' lovinterp3: it has to be called if prob_sampl is TRUE'
     STOP
  ENDIF
  ! Initialization
  succ = .FALSE. 
  nd   = 6 + nls
  ! Nearest integere to rindex
  ir = NINT(rindex)
  ! Compute imu, the nearest LOV index to ir
  IF(ir.LT.imim)THEN 
     imu = imim
  ELSEIF(ir.GT.imip)THEN 
     imu = imip
  ELSE  
     imu = ir 
  ENDIF
  ! Distance (with sign) from rindex to imu
  s = rindex-imu
  IF(ABS(s).GT.1.d0)THEN
     WRITE(*,*) ' lovinterp3: extrapolation index beyond computed LOV ', ir, imim, imip
     STOP
  ENDIF
  ! WARNING: use convex combination is easier
  IF(lin_lov .AND. .NOT.scatterplane)THEN  
     CALL weak_dir(unm(imu)%g(1:nd,1:nd),wdir,sdir,-1,elm(imu)%coo,elm(imu)%coord,units,nd)
     ! Find interpolation interval length
     IF(s.GT.0.d0)THEN
        deltasig = sigmavalv(imu+1)-sigmavalv(imu)
     ELSEIF(s.LT.0.d0)THEN
        deltasig = sigmavalv(imu)-sigmavalv(imu-1)
     ELSE 
        WRITE(*,*) ' lovinterp3: interpolation with integer rindex =', rindex
        STOP
     ENDIF
     ! Find interpolated orbit
     el0       = elm(imu)
     el0%coord = elm(imu)%coord + units(1:6)*wdir(1:6)*sdir*deltasig*s
     IF(nd.GT.6)THEN
        DO j=1,nls
           dyn%dp(ls(j)) = dpm(ls(j),imu)+units(6+j)*wdir(6+j)*sdir*deltasig*s
        ENDDO
     ENDIF
     unc0 = unm(imu)
     succ = .TRUE.
     IF(s.GT.0.d0)THEN
        sigmaout = (1.d0-s)*sigmavalv(imu)+s*sigmavalv(imu+1)
     ELSE
        sigmaout = (1.d0+s)*sigmavalv(imu)-s*sigmavalv(imu-1)
     ENDIF
  ELSEIF(lin_lov .AND. scatterplane)THEN
     ! Find interpolation interval length
     IF(s.GT.0.d0)THEN
        deltasig = sigmavalv(imu+1)-sigmavalv(imu)
     ELSEIF(s.LT.0.d0)THEN
        deltasig = sigmavalv(imu)-sigmavalv(imu-1)
     ELSE
        WRITE(*,*) ' lovinterp3: interpolation with integer rindex =', rindex
        STOP
     ENDIF
     ! Find interpolated orbit
     el0       = elm(imu)
     el0%coord = elm(imu)%coord+vvv_tp(1:6)*deltasig*s
     IF(nd.GT.6)THEN
        DO j=1,nls
           dyn%dp(ls(j)) = dpm(ls(j),imu)+vvv_tp(6+j)*deltasig*s
        ENDDO
     ENDIF
     unc0 = unm(imu)
     succ = .TRUE.
     IF(s.GT.0.d0)THEN
        sigmaout = (1.d0-s)*sigmavalv(imu)+s*sigmavalv(imu+1)
     ELSE
        sigmaout = (1.d0+s)*sigmavalv(imu)-s*sigmavalv(imu-1)
     ENDIF
  ELSEIF(.NOT.lin_lov)THEN
     CALL weak_dir(unm(imu)%g(1:nd,1:nd),wdir,sdir,-1,elm(imu)%coo,elm(imu)%coord,units,nd) 
     ! Anti-reversal check
     IF(imu.LT.imip)THEN
        diff(1:6) = elm(imu+1)%coord-elm(imu)%coord
        DO j=7,nd
           diff(j) = dpm(ls(j-6),imu+1)-dpm(ls(j-6),imu)
        ENDDO
     ELSE
        diff(1:6) = elm(imu)%coord-elm(imu-1)%coord
        DO j=7,nd
           diff(j) = dpm(ls(j-6),imu)-dpm(ls(j-6),imu-1)
        ENDDO
     ENDIF
     DO j=1,nls
        dyn%dp(ls(j)) = dpm(ls(j),imu)
     END DO
     IF(DOT_PRODUCT(diff(1:nd),wdir(1:nd)).LT.0.d0) wdir(1:nd) = -wdir(1:nd)
     ! Find interpolation interval
     IF(s.GT.0.d0)THEN
        deltasig = sigmavalv(imu+1)-sigmavalv(imu)
     ELSEIF(s.LT.0.d0)THEN
        deltasig = sigmavalv(imu)-sigmavalv(imu-1)
     ELSE
        WRITE(*,*) ' lovinterp3: interpolation with integer rindex =', rindex
        STOP
     ENDIF
     ! Non-linear interpolation   
     CALL prop_sig(batch,elm(imu),el0,s,deltasig,m_m,obs_m,obsw_m,wdir,sdir,units,fail,nd)
     IF(fail)THEN 
        WRITE(*,*) ' lovinterp3: called at index ',rindex,', and get a hyperbolic orbit' 
        WRITE(*,*) el0%coo,el0%coord 
        RETURN 
     ELSEIF(bizarre(el0,ecc))THEN 
        WRITE(*,*) ' lovinterp3: called at index ',rindex,' and get a bizarre orbit' 
        WRITE(*,*) el0%coo,el0%coord, ecc
        RETURN 
     ENDIF
     ! Constrained corrections
     CALL constr_fit(m_m,obs_m,obsw_m,el0,wdir,elc,unc0,csinew,delnew,rmshnew,nused,succ,nd)
     IF(.NOT.succ) THEN 
        WRITE(*,*) ' lovinterp3: constr_fit failed for index = ', rindex 
        RETURN 
     ELSE
        el0  = elc 
        succ = .TRUE. 
        IF(s.GT.0.d0)THEN
           sigmaout = (1.d0-s)*sigmavalv(imu)+s*sigmavalv(imu+1)
        ELSE
           sigmaout = (1.d0+s)*sigmavalv(imu)-s*sigmavalv(imu-1)
        ENDIF
     ENDIF
  ELSE
     WRITE(*,*) ' lovinterp3: case not yet invented'
     STOP
  ENDIF
END SUBROUTINE lovinterp3

! ====================================================                  
! LOVOBS                                                                
! ====================================================                  
! provides number and arc of observational data                         
SUBROUTINE lovobs(m0,nrej,calend1,calend2) 
  USE obssto
! ======================OUTPUT==========================
  INTEGER, INTENT(OUT) :: m0,nrej   ! number of obs., number rejected
  CHARACTER(LEN=14), INTENT(OUT) :: calend1,calend2 ! calendar dates of first 
                                                      ! and last
! =====================================================
!  INCLUDE 'parobx.h90'
  INTEGER :: j, ind(nobx)  ! loop index, sort index      
! =====================================================   
! calendar date
  CALL heapsort(obs_m(1:m_m)%time_tdt,m_m,ind)
  CALL calendwri(obs_m(ind(1))%time_tdt,calend1) 
  CALL calendwri(obs_m(ind(m_m))%time_tdt,calend2) 
  m0=m_m 
  nrej=0 
  DO j=1,m_m 
     IF(obsw_m(j)%sel_coord.eq.0)THEN
         nrej=nrej+1
     ENDIF
  ENDDO
END SUBROUTINE lovobs

END MODULE multiple_sol
! ====================================================================  
! Graham- Schmidt procedure to generate an orthonormal basis v          
! starting from 1  n-vector a                                           
! The new basis must be such that the first  vector is a                
SUBROUTINE graha_1(a,n,v) 
  implicit none 
  integer, intent(in) ::  n 
  double precision, intent(in) ::  a(n)
  double precision, intent(out) :: v(n,n)
! end interface 
  integer j,jok,jj,i 
  double precision cc,cc1,epsi,vl 
  integer,parameter :: nx=10  ! workspace
  double precision ws(nx) 
  logical ize 
! dimension check                                                       
  if(n.gt.nx)then 
     write(*,*)'n =',n,' larger than nx=',nx,' in graha' 
     stop 
  endif
! selection of the control for "zero" vectors                           
  cc1=sqrt(DOT_PRODUCT(a(1:n),a(1:n)))  ! cc1=sqrt(prscag(n,a,a)) 
  epsi=1.d-12*cc1 
  if(epsi.eq.0.d0)then 
     write(*,*)' a has rank zero' 
!        stop                                                           
  endif                                                                  
! V1 is the versor of A                                                 
  call versor(n,a,epsi,v(1,1),vl,ize) 
  if(ize)then 
     write(*,*)' first vector of a is too small' 
!        stop                                                           
  endif
! we now use the vectors of the canonic basis to supplement the span of 
  jok=0 
  do 1 j=1,n 
! remove the components along span(A), that is along V1                 
     cc1=-v(j,1) 
     do  i=1,n 
        ws(i)=cc1*v(i,1) 
     enddo
     ws(j)=ws(j)+1.d0 
     call versor(n,ws,epsi,v(1,2+jok),vl,ize) 
     if(.not.ize)then 
! now V(3+jok) is orthogonal to span(A); remove the components along    
! the previous ones (unless it is the first)                            
        if(jok.gt.0)then 
           do  jj=1,jok 
              cc=-DOT_PRODUCT(v(1:n,2+jok),v(1:n,1+jj)) 
                   ! cc=-prscag(n,v(1,2+jok),v(1,1+jj)) 
              v(1:n,2+jok)=v(1:n,2+jok)+cc*v(1:n,1+jj) 
                    ! call lincog(n,v(1,2+jok),1.d0,v(1,1+jj),cc,v(1,2+jok)) 
           enddo
           call versor(n,v(1,2+jok),epsi,v(1,2+jok),vl,ize) 
           if(ize)then 
              goto 1 
           endif
        endif
! the new versor is a good one                                          
        jok=jok+1 
        if(jok.eq.n-1)then 
           goto 2 
        endif
     endif
1 continue 
2 continue 
  if(jok.lt.n-1)then 
     write(*,*)' graha_1: something went wrong, jok=',jok 
  endif
END SUBROUTINE graha_1

! ================================================================      
! SUBROUTINE aftclo2v                                                  
! select time interval to get after the close approach                  
! taking into account the relative velocity w.r. to Earth
! ================================================================        
  SUBROUTINE aftclo2v(iplam,t0,t1,t2,v_inf0,tbefore,tafter) 
    USE planet_masses
! ======================INPUT===================================== 
    INTEGER, INTENT(IN) :: iplam             ! planet number
    DOUBLE PRECISION, INTENT(IN) :: t0       ! time of initial 
                                             ! conditions
    DOUBLE PRECISION, INTENT(IN) :: t1,t2    ! time of two known 
                                             ! close approaches
    DOUBLE PRECISION, INTENT(IN) :: v_inf0   ! velocity 
! =====================OUTPUT=====================================
    DOUBLE PRECISION, INTENT(OUT) :: tafter    ! time "after"
    DOUBLE PRECISION, INTENT(OUT) :: tbefore   ! time "before"
                                !  depending upon sense of propagation  
! ===================END INTERFACE=================================      
    DOUBLE PRECISION :: delt_tp       ! time interval to 
                                      ! exit from TP disk    
! ===================================================================   
! warning: really done only for Earth                                   
    IF(ordnam(iplam).ne.'EARTH') THEN 
       WRITE(*,*)' aftclo2v: not to be used for planet ',ordnam(iplam) 
    ENDIF
! time interval to exit from TP disk   
    IF(v_inf0.gt.0.d0)THEN                                 
       delt_tp=2*dmin(iplam)/v_inf0
!       delt_tp=4*dmin(iplam)/v_inf0
    ELSE
       delt_tp=365.25d0
    ENDIF
! forced to avoid infinte intervals for v-inf=0
    delt_tp=MIN(delt_tp,365.25d0)
! forced to avoid short intervals for fast encounters                   
    IF(delt_tp.lt.20.d0)delt_tp=20.d0 
!    IF(delt_tp.lt.50.d0)delt_tp=50.d0 
! time interval to be clear out of the TP disk is given as deltat       
    IF(t0.lt.t1)THEN 
! future close approaches                                               
       tafter=max(t1,t2)+delt_tp 
       tbefore=min(t1,t2)-delt_tp 
    ELSE 
! past close approaches                                                 
       tafter=min(t1,t2)-delt_tp 
       tbefore=max(t1,t2)+delt_tp 
    ENDIF
  END SUBROUTINE aftclo2v
! =========================================================================

! ================================================                      
! propagation of catalogs (multiline format)                            
program catpro 
  USE fund_const
  USE output_control
  USE orbit_elements
  USE propag_state
  USE dyn_param
  USE name_rules, ONLY: idname_len
! ================================================                      
  implicit none 
! ================================================                      
  CHARACTER*6 progna 
  CHARACTER*80 run
  INTEGER iuncat, iuncat0,iunout  ! input and output units
  CHARACTER*80 catnam0,catnam1,dummyfile !catalogs: input, output, dummyname
!asteroid name, short and long format, length
  CHARACTER*(idname_len) name 
  CHARACTER*(idname_len) namid
  INTEGER ld,le
  CHARACTER*10 rsys,epoch ! for catalog file header
  CHARACTER*3 eltype ! for rdorb call
  DOUBLE PRECISION t,eq(6),covar(ndimx,ndimx),normal(ndimx,ndimx),h,gmag,mass
  LOGICAL defcov,defnor
  INTEGER no
  DOUBLE PRECISION tref ! reference epoch for propagation 
  TYPE(orbit_elem) el0,el,ek !elements (in, out, keplerian)
  TYPE(orb_uncert) unc,unc0 ! uncertainties (in ,out)
  INTEGER fail_flag ! for coordinate changes
  LOGICAL eof,covpro,ok, err  ! end input file, propagate covariance, data available
  INTEGER i,ndone, nlsloc, lsloc(ndyx) ! counters
  INTEGER nd
! =================input options===============================
  progna='catpro' 
  run=progna                                        
  call cinopt(progna,run,iunout,tref,catnam0,catnam1,covpro) 
! check availability of JPL ephemrides                                  
  CALL chetim(tref,tref,ok) 
  IF(.not.ok)THEN 
     WRITE(*,*)' JPL ephemerides not available for t=',tref 
     STOP 
  ENDIF
! ================open generic output files ========================
  call filopn(ipirip,'catpro.pro','unknown') 
  call filopn(ierrou,'catpro.err','unknown') 
  numerr=0 
  call filopn(iuncla,'catpro.clo','unknown') 
  numcla=0 
! =============================================                         
! get reference system, etc. from old catalog to write into new catalog.
  CALL oporbf(catnam0,0) 
  CALL rdorb(namid,eq,eltype,t,covar,defcov,normal,defnor,h,gmag,   &
     &     mass,rsys,epoch,no,eof)                                      
  CALL clorbf 
! ===================================================================== 
! open new catalog and write header            
  call filopn(iuncat,catnam1,'unknown') 
  IF(covpro)THEN 
     CALL wromlh(iuncat,rsys,epoch) 
  ELSE 
     CALL wro1lh(iuncat,rsys,epoch,'KEP') 
  ENDIF
! ======================================================================
! Begin main loop  
  CALL filopn(iuncat0,catnam0,'old')
  CALL oporbf(catnam0,iuncat0) 
  ndone=0 
  i=0
10 CONTINUE 
! input orbit at common epoch, with covariance                          
  REWIND (ipirip) 
  CALL read_elems(el0,namid,eof,nlsloc,err,lsloc,UNC=unc0,FILE=dummyfile,UNIT=iuncat0)
  IF(dyn%nmod.GT.0)THEN
     nls=nlsloc
     ls=lsloc
  END IF
  nd=6+nls
  IF(eof) GOTO 5 
  i=i+1 ! counter of orbits read
! convert name in case it is of the form nam0=namp                      
  IF(covpro)THEN
     CALL rmsp(namid,ld)
     name=namid(1:ld)
  ELSE
     CALL rmsp(namid,le)
     ld=index(namid,'=')-1 
     IF(ld.lt.0)THEN 
        name=namid(1:le) 
     ELSE 
        name=namid(1:ld) 
     ENDIF
  ENDIF
! check availability of JPL ephemrides                                  
  CALL chetim(el0%t,tref,ok) 
  IF(.not.ok)THEN 
     WRITE(*,*)' JPL ephemerides not available for times ',el0%t,tref 
     GOTO 10 
  ENDIF
! ================================================================      
! propagation to time tref                                      
  IF(covpro.and.unc0%succ)THEN 
     CALL pro_ele(el0,tref,el,unc0,unc) 
  ELSEIF(covpro.and..not.unc0%succ)THEN
     WRITE(*,*)' matrix missing for ',name
     GOTO 10 
  ELSEIF(.not.covpro)THEN 
! state vector only                                                     
     CALL pro_ele(el0,tref,el)
  ENDIF 
! write output in multiline format (if covariance is required),         
! single line format otherwise.                                         
  IF(covpro)THEN 
     CALL write_elems(el,name,'ML',UNC=unc,FILE=dummyfile,UNIT=iuncat)
     !CALL write_elems(el,name,'ML',dummyfile,iuncat,unc) !INCFIT???
  ELSE 
! output single line catalog record                                     
     CALL coo_cha(el,'KEP',ek, fail_flag)
     IF(fail_flag.lt.4) THEN 
        CALL write_elems(ek,name,'1L',FILE=dummyfile,UNIT=iuncat)
     ELSE
        WRITE(*,*)' catpro: cannot handle comets in 1L format', el
     ENDIF
  ENDIF
  ndone=ndone+1 
  GOTO 10
!     end of input file                                                 
5 CALL clorbf 
  CALL filclo(iuncat,' ')
  CALL filclo(iuncla,' ')
  CALL filclo(ipirip,'DELETE')
  WRITE(*,*)' total number of orbits ',i,' propagated ',ndone 
  CALL rmsp(run,le)
  CALL final_rep(ndone,0,run,le,iunout,name)

  CONTAINS

! ===================================================================   
! CINOPT                                                                
! ===================================================================   
! input options, for the catalog propagator main program                
! input: progna = program name (6 characters)                           
!        run    = run identifier (80 characters, up to 76 non-blank)    
SUBROUTINE cinopt(progna,run,iunout,tref,catnam0,catnam1,covpro) 
  USE output_control
  implicit none   
  character*6,INTENT(IN) :: progna 
  character*80,INTENT(IN) :: run 
  integer, intent(out) :: iunout ! output units
  double precision, intent(out):: tref 
  CHARACTER*80, INTENT(OUT) :: catnam0,catnam1 
  logical, intent(out) :: covpro 
! ==========END INTERFACE============================================  
  integer le
  character*90 file 
  LOGICAL ireq,found,ex
  CHARACTER*60 comment 
! =============================
  CALL initopt(progna,run,'mop')                                         
! read option for physical model and integration method                 
  CALL rmodel(1) 
! ==============================                                        
! initialisations for Gauss method                                      
! CALL iodini                
! =============================                                         
! Output files: for control and results
  file=run//'.mou' 
  CALL rmsp(file,le) 
  call filopn(iunout,file,'UNKNOWN')
  iun_log=iunout ! later iunout to be abolished
! ================reference time============================            
  ireq=.true. 
  comment='reference time'
  CALL input_rea_opt(progna,'tref',tref,ireq,found,comment,iunout)
! ====================================================                  
! read the name of the elements file, inquire                           
  ireq=.true. 
  comment='input elements catalog'
  CALL input_cha_opt(progna,'catnam0',catnam0,ireq,found,comment,iunout)
  CALL rmsp(catnam0,le) 
  INQUIRE(file=catnam0(1:le),exist=ex) 
  IF(.not.ex)THEN 
     write(*,*) ex,' catalog not found: ',catnam0(1:le) 
     STOP 
  ENDIF
  ireq=.true. 
  comment='output elements catalog'
  CALL input_cha_opt(progna,'catnam1',catnam1,ireq,found,comment,iunout)
! ================generate catalog with covariances=============        
  ireq=.true.
  comment='generate covariance logical control'
  CALL input_log_opt(progna,'covpro',covpro,ireq,found,comment,iunout)
! ====================================================
END SUBROUTINE cinopt


END PROGRAM catpro

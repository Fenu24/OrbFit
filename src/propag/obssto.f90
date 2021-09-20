!=================MODULE obssto===================
! obssto.o: 
!	../suit/astrometric_observations.mod
!=================================================
MODULE obssto 
  USE astrometric_observations 
  USE orbit_elements
  IMPLICIT NONE 
  PUBLIC retinobs
!===============observation data needed by fmulti=======================       
  INCLUDE 'parobx.h90'
  INTEGER  :: m_m       ! actual number of observations
  TYPE(ast_obs), DIMENSION(nobx), SAVE  :: obs_m    ! observations
  TYPE(ast_wbsr), DIMENSION(nobx) :: obsw_m
  DOUBLE PRECISION :: rms_m,rmsmag_m 
CONTAINS

! ===========================================                           
! RETINOBS                                                                
! observation input control routine 
! version 3.0 fortran90 A. Milani June 2003                                    
! reads the .rwo,file with the errror_model given
  SUBROUTINE retinobs(obsdir,nam0,obs0,error_model,rms,rmsmag)
    USE astrometric_observations
    USE output_control
! ==============INPUT==================                                 
    CHARACTER*60, INTENT(IN) :: obsdir ! input directory (all files nam0.rwo must be there) 
    CHARACTER*9,INTENT(IN) :: nam0 ! asteroid name (9 characters)
! =============OUTPUT================================                   
    LOGICAL,INTENT(OUT) :: obs0 ! successfuj input flag
! ===== observational data ===========================                  
    CHARACTER*(20),INTENT(OUT) :: error_model ! error model file name
    DOUBLE PRECISION, INTENT(OUT) :: rms,rmsmag ! RMS of astrometric, photometric residuals
! ===========END INTERFACE========================= 
    CHARACTER*60 rwofi0 ! .rwo file name
    INTEGER lrwo,le 
    LOGICAL rwo 
    CHARACTER*60 :: obsdir1
    INCLUDE 'sysdep.h90' ! directory char
    INTEGER j ! loop index
! =============EXECUTION BEGINS======================
! check for existence of the corresponding .rwo file                    
! ========================================================              
!  compute rwo file name; no subdirectory structure
    obsdir1=obsdir
    CALL rmsp(obsdir1,le) 
    rwofi0=obsdir1(1:le)//'/'//nam0//'.rwo'
    CALL rmsp(rwofi0,lrwo)
! existence of .rwo file                                        
    INQUIRE(file=rwofi0(1:lrwo),exist=rwo) 
! ========================================================              
! select operations mode                                                
! ========================================================              
    IF(.not.rwo)THEN 
       WRITE(*,*) 'No .rwo file for ',nam0,' ',rwofi0(1:lrwo) 
       obs0=.false. 
       RETURN 
! ========================================================              
    ELSE 
! read .rwo, and store data in final array
       CALL read_rwo(rwofi0(1:lrwo),obs_m,obsw_m,m_m,error_model,rms,rmsmag)
       WRITE(iun_log,*)'rearwo: ',m_m,' obs from ',rwofi0(1:lrwo) 
       IF(m_m.eq.0)THEN 
          obs0=.false. 
       ELSE 
          obs0=.true. 
       ENDIF
    ENDIF
! get asteroid radius (if necessary) before returning
    CALL aster_radius(obs_m%objdes,obs_m%type,m_m)
  END SUBROUTINE retinobs

END MODULE obssto









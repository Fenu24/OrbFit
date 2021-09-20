MODULE force_sat
! module not in use in version 4.3
  USE fund_const
  USE planet_masses, ONLY: gmearth,gmoon,dmea
  IMPLICIT NONE 
  PRIVATE

! public routines
  PUBLIC :: eamoon_mass, forcesat
! shared data
  

CONTAINS

! ===========================================================           
! FORCESAT : accelerations acting on a satellite of Earth                   
! ===========================================================           
! version 3.6.1; 29 October 2008     
  SUBROUTINE forcesat(x,v,t0,f,nd,idc,xxpla,ips,imem, derf0)
!    USE perturbations
    USE spher_harm
    USE iers_ser
! ======INPUT===================                                        
! dimension of position vector                                          
    INTEGER, INTENT(IN) :: nd 
! Position, velocity,  time                                             
    REAL(KIND=dkind), INTENT(IN) :: x(nd), v(nd) ,t0 
! flag for recomputation, memory location                               
    INTEGER, INTENT(IN) :: ips,imem 
! WARNING: for now storage/partial recomputation not implemented
! ======OUTPUT===================                                       
! acceleration                                                          
    REAL(KIND=dkind) :: f(nd) 
! Positions and vel of the planet involved in  close-app                
! stored only if idc.ne.0                                               
    INTEGER :: idc 
    REAL(KIND=dkind) :: xxpla(6) 
    REAL(KIND=dkind), OPTIONAL :: derf0(3,3)
! ======END INTERFACE============== 
    STOP ' ****** force_sat not used *******'
  END SUBROUTINE forcesat

  SUBROUTINE eamoon_mass
    STOP ' ****** eamoon_mass not used ******'
  END SUBROUTINE eamoon_mass

END MODULE force_sat

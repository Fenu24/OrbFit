! ===================================================================== 
! STACOP                                                                
! ===================================================================== 
! copy state variables to current state                                 
! ================INTERFACE=========================================    
subroutine stacop(icop,el0,unc0,csino0,delno0,                   &
     &         el,unc,csinoc,delnoc)                                 
  USE orbit_elements
  implicit none 
  integer, intent(IN) ::  icop
  TYPE(orbit_elem) :: el0, el
  TYPE(orb_uncert) :: unc0, unc 
  double precision csino0,delno0,csinoc,delnoc 
! =================END INTERFACE====================================    
! function: copy to/from                                                
  if(icop.eq.1)then 
     el=el0 
     unc=unc0
     csinoc=csino0 
     delnoc=delno0 
  elseif(icop.eq.2)then 
     el0=el 
     unc0=unc
     csino0=csinoc 
     delno0=delnoc 
  else 
     write(*,*)' stacop: copying not known, icop=',icop 
     stop 
  endif
END SUBROUTINE stacop



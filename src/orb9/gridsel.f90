! ****************************************************************      
!                                                                       
!  proper elements assigned in a grid                                   
!                                                                       
SUBROUTINE gridsel(ngrid,ita,aa0,deltaa,ite,e0,delte,iti,aai0,delti,aw,ao)
  USE fund_const
  IMPLICIT NONE
  INTEGER, INTENT(OUT) :: ngrid, ita, ite,iti
  DOUBLE PRECISION, INTENT(OUT) :: aa0,deltaa,e0,delte,aai0,delti,aw,ao
! END INTERFACE
  INTEGER isel
  DOUBLE PRECISION amax, amin, emax, emin, aimax, aimin, aod, awd
!499  write(*,800)ngridx                                                
!800  format(' chose ngrid, max=',i5,' choose axes 1=ie 2=ai 3=ae')     
  read(1,701)ngrid 
  read(1,701)isel 
701 format(i8) 
  if(isel.eq.1)then 
     ita=1 
     ite=ngrid 
     iti=ngrid 
  else 
     if(isel.eq.2)then 
        ita=ngrid 
        ite=1 
        iti=ngrid 
     else 
        if(isel.eq.3)then 
           ita=ngrid 
           ite=ngrid 
           iti=1 
        else 
           write(*,501)isel 
501        format(' isel=',i4,'  option not supported') 
!          goto 499                                                   
           stop 
        endif
     endif
  endif
  if(ita.gt.1)then 
!550      write(*,602)                                                  
!602      format(' give range in a, amin, amax')                        
     read(1,702)amin,amax 
702  format(2f10.4) 
     write(*,663)amin,amax 
663  format(' range a',2x,2f10.4) 
!         if(amin.gt.amax)then                                          
!              write(*,662)amin,amax                                    
!662           format(' wrong a',2x,2f10.4)                             
!              stop                                                     
!         endif                                                         
     deltaa=(amax-amin)/(ita-1) 
     aa0=amin 
  else 
!         write(*,503)                                                  
!503      format(' give a')                                             
     read(1,702)aa0 
     deltaa=0.d0 
  endif
  if(ite.gt.1)then 
!560      write(*,512)                                                  
!512      format(' give range in e, emin, emax')                        
     read(1,702)emin,emax 
     write(*,563)emin,emax 
563  format(' range e',2x,2f10.4) 
!         if(emin.gt.emax)then                                          
!             write(*,562)emin,emax                                     
!562          format(' wrong e',2x,2f10.4)                              
!             goto 499                                                  
!             stop                                                      
!         endif                                                         
     delte=(emax-emin)/(ite-1) 
     e0=emin 
  else 
!        write(*,513)                                                   
!513      format(' give e')                                             
     read(1,702)e0 
     delte=0.d0 
  endif
  if(iti.gt.1)then 
!570      write(*,522)                                                  
!522      format(' give range in i, imin, imax')                        
     read(1,702)aimin,aimax 
     write(*,573)aimin,aimax 
573  format(' range i',2x,2f10.4) 
!         if(aimin.gt.aimax)then                                        
!              write(*,572)aimin,aimax                                  
!572           format(' wrong i',2x,2f10.4)                             
!              goto 499                                                 
!               stop                                                    
!          endif                                                        
     delti=(aimax-aimin)/(iti-1) 
     aai0=aimin 
  else 
!         write(*,523)                                                  
!523      format(' give i')                                             
     read(1,702)aai0 
     delti=0.d0 
  endif
!  angles are assigned at a fixed value (read in degrees)               
  read(1,710)awd,aod 
710 format(2f10.4) 
  aw=awd*radeg 
  ao=aod*radeg 
END SUBROUTINE gridsel

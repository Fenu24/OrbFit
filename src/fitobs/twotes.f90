!   twotes                                                              
!  test derivatives of element functions                                
SUBROUTINE twotes(m,tau,ioco,el,pos,vel,nd) 
  USE orbit_elements
  USE pred_obs
  USE dyn_param
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: m
  TYPE(orbit_elem), INTENT(IN) :: el
  INTEGER, INTENT(IN), DIMENSION(m) :: ioco
  DOUBLE PRECISION, INTENT(IN), DIMENSION(m) :: tau 
  DOUBLE PRECISION, INTENT(IN), DIMENSION(3,m) :: pos, vel
  INTEGER,INTENT(IN) :: nd ! number of parameters to be solved
  DOUBLE PRECISION eq(6),eqp(6)
  DOUBLE PRECISION dade(nd),ddde(nd), alj, dej, dal, dde, alp, dep, delta
  DOUBLE PRECISION dadep(nd),dddep(nd), erdera, erderd
  INTEGER jj, ider, k
  logical twobo 
  TYPE(orbit_elem) elp
  twobo=.false. 
  open(19,file='twotes.out',status='unknown') 
! size of increment                                                     
  write(*,*)' step size?' 
  read(*,*)delta 
  write(19,*)delta 
!                                                                       
  jj=1 
  ider=1
!  compute obs, first derivatives                     
99 CONTINUE
!  IF(iobs(jj)/1000.eq.1)THEN 
  call alph_del (el,tau(jj),ioco(jj),pos(:,jj),vel(:,jj),ider,twobo,nd,alj,dej,dade,ddde)
!      ELSEIF(iobs(jj)/1000.eq.2)THEN 
!         call rrdot(eq,iobs(jj),t0,tau(jj),ioco(jj),alj,dej,dade,ddde,  &
!     &        ider,twobo,ddade,dddde)                                   
!      ELSE 
!         write(*,*)' wrong iobs ', iobs(jj),jj 
!         RETURN 
!      ENDIF 
!  change elements and compute obs. and first derivatives               
!  change one element at the time
  DO 2 k=1,6
     elp=el
     elp%coord(k)=elp%coord(k)+delta
!  compute obs. and first deriv.                                        
     ider=1 
!  IF(iobs(jj)/1000.eq.1)THEN 
     call alph_del(elp,tau(jj),ioco(jj),pos(:,jj),vel(:,jj),ider,twobo,nd,alp,dep,dadep,dddep)
!        ELSEIF(iobs(jj)/1000.eq.2)THEN 
!           call rrdot(eqp,iobs(jj),t0,tau(jj),ioco(jj),alp,dep,dadep,   &
!     &          dddep,ider,twobo,ddade,dddde)                           
!        ENDIF 
!  rapporto incrememtale                                                
     dal=alp-alj 
     dde=dep-dej 
     erdera=dade(k)-dal/delta 
     erderd=ddde(k)-dde/delta 
     write(*,*)' partial w.r. to element',k
     write(*,*) ' values from var. eq. ',dade(k),ddde(k) 
     write(*,*) ' values from fin.dif. ',dal/delta,dde/delta 
     write(*,*) ' errors ',erdera,erderd 
     write(19,*)' partial w.r. to element',k 
     write(19,*) dade(k),ddde(k) 
     write(19,*) dal/delta,dde/delta 
     write(19,*) erdera,erderd 
  
2 ENDDO
  close(19) 
END SUBROUTINE twotes

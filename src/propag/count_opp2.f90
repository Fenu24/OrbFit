SUBROUTINE count_opp2(obs,obsw,m,el,nopp,napp,int_con,int_app,nopp_pot,napp_pot)
  USE astrometric_observations
  USE pred_obs
  USE fund_const
  USE orbit_elements
  USE propag_state, ONLY: pro_ele
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: m ! no. observations
  TYPE(ast_obs),DIMENSION(m),INTENT(IN) :: obs ! observations 
  TYPE(ast_wbsr),DIMENSION(m),INTENT(IN) :: obsw ! observation weights 
  TYPE(orbit_elem), INTENT(IN) :: el !  elements, epoch time
  INTEGER,INTENT(OUT) :: nopp,napp ! number of oppositons, of apparitions
  INTEGER,PARAMETER:: nstepx=10000
  DOUBLE PRECISION,DIMENSION(nstepx,2),INTENT(OUT) :: int_con,int_app
  INTEGER,INTENT(OUT) :: nopp_pot,napp_pot
! ==============end interface=================================
  INTEGER j,k,kk,inl,idt,nstep,nsup,np60,nm60,jopp,japp
  DOUBLE PRECISION h,g,alpha,delta,hmagn,dt
  DOUBLE PRECISION,DIMENSION(nstepx) :: pha,elo,dis,tt,tsup,tp60,tm60
  INTEGER,PARAMETER :: istep=10
  DOUBLE PRECISION,PARAMETER :: cos0=0.64278d0
  INCLUDE 'parobx.h90'
  INTEGER time_ord(nobx),sel(nobx)
  DOUBLE PRECISION tobs(nobx)
  DOUBLE PRECISION,DIMENSION(nstepx,2) :: int_con_pot,int_app_pot
  TYPE(orbit_elem) el1
! =============================================================
! sort times
  CALL heapsort(obs(1:m)%time_tdt,m,time_ord)
  DO j=1,m
    tobs(j)=obs(time_ord(j))%time_tdt
    sel(j)=obsw(time_ord(j))%sel_coord
  ENDDO
  h=0
  g=0.15
  inl=1
  tt(1)=tobs(1)
  CALL pro_ele(el,tt(1),el1)
  CALL set_restart(.true.) 
  CALL predic_obs(el1,500,tt(1),'O',alpha,delta,hmagn,inl,      &
     &    PHA0=pha(1),ELO0=elo(1),DIS0=dis(1))
  IF(elo(1).gt.0.d0)pha(1)=-pha(1)
  IF(cos(elo(1)).gt.cos0)THEN
     WRITE(*,*)'count_opp: first obs. close to the Sun', elo(1)*degrad
     STOP
  ENDIF   
  CALL set_restart(.false.) 
  dt=tobs(m)-tobs(1)
  idt=dt
  nstep=idt/istep+1
  DO j=1,nstep
     tt(j+1)=tobs(1)+j*istep
     CALL predic_obs(el1,500,tt(j+1),'O',alpha,delta,hmagn,inl,          &
     &    PHA0=pha(j+1),ELO0=elo(j+1),DIS0=dis(j+1))
     IF(elo(j+1).gt.0.d0)pha(j+1)=-pha(j+1)
     CALL set_restart(.false.)
  ENDDO
! count superior conjunctions and apparitions
  nsup=0
  np60=0
  nm60=0
  DO j=2,nstep+1
! count superior conjunctions
    IF(elo(j)*elo(j-1).le.0.d0.and.abs(elo(j)).le.pig/2.and.        &
  &    abs(elo(j-1)).le.pig/2.and.dis(j).gt.1.d0)THEN
! superior conjunction
       nsup=nsup+1
       tsup(nsup)=(tt(j)+tt(j-1))/2.d0
    ENDIF
! counting possible apparitions
    IF((cos(elo(j))-cos0).lt.0.d0.and.(cos(elo(j-1))-cos0).gt.0.d0)THEN
! crossing the elongation=+- 60 deg line becoming visible 
       np60=np60+1
       tp60(np60)=(tt(j)+tt(j-1))/2.d0
    ELSEIF((cos(elo(j))-cos0).gt.0.d0.and.(cos(elo(j-1))-cos0).lt.0.d0)THEN
! crossing the elongation=+- 60 deg line becoming invisible
       nm60=nm60+1
       tm60(nm60)=(tt(j)+tt(j-1))/2.d0
    ENDIF
  ENDDO
! create intervals for oppositions
  IF(nsup.gt.0)THEN
     int_con_pot(1,1)=tt(1)-100.d0
     int_con_pot(1,2)=tsup(1)
     DO k=2,nsup
        int_con_pot(k,1)=tsup(k-1)
        int_con_pot(k,2)=tsup(k)
     ENDDO
     int_con_pot(nsup+1,1)=tsup(nsup)
     int_con_pot(nsup+1,2)=tt(nstep)+100.d0
  ELSE
     int_con_pot(1,1)=tt(1)-100.d0 
     int_con_pot(1,2)=tt(nstep)+100.d0
  ENDIF
  nopp_pot=nsup+1
! count oppositions observed
  nopp=0
  jopp=0
  DO j=1,m
    DO k =jopp+1,nsup+1
       IF(tobs(j).gt.int_con_pot(k,1).and.tobs(j).lt.int_con_pot(k,2).and.sel(j).gt.0)THEN 
          nopp=nopp+1
          jopp=k
          int_con(k,1)=tobs(j)
          GOTO 1
       ELSEIF(tobs(j).lt.int_con_pot(jopp,2))THEN
          int_con(k-1,2)=tobs(j)
       ENDIF
    ENDDO
 1  CONTINUE
  ENDDO
! create intervals for apparitions
  IF(np60.gt.nm60+1)THEN
     WRITE(8,*)'count_opp: np60=',np60,' nm60=',nm60
  ENDIF
  napp_pot=np60+1
  IF(nm60.ge.1)THEN
     int_app_pot(1,1)=tt(1)-100.d0
     int_app_pot(1,2)=tm60(1)
     DO k=2,np60
        int_app_pot(k,1)=tp60(k-1)
        int_app_pot(k,2)=tm60(k)
     ENDDO
     int_app_pot(napp_pot,1)=tp60(np60)
     int_app_pot(napp_pot,2)=tt(nstep)+100.d0
  ELSE 
     int_app_pot(1,1)=tt(1)-100.d0
     int_app_pot(1,2)=tt(nstep)+100.d0
  ENDIF 
! count apparitions
  napp=0
  japp=0
  DO j=1,m
    DO k =japp+1,napp_pot
       IF(tobs(j).gt.int_app_pot(k,1).and.tobs(j).lt.int_app_pot(k,2).and.sel(j).gt.0)THEN 
          napp=napp+1
          japp=k
          int_app(k,1)=tobs(j)
          GOTO 2
       ELSEIF(tobs(j).lt.int_app_pot(jopp,2))THEN
          int_app(k,2)=tobs(j)
       ENDIF
    ENDDO
 2  CONTINUE
  ENDDO

END SUBROUTINE count_opp2


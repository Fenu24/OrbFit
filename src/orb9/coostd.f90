! ========================================                              
!   COOSTP case for planets                                             
!   reduction of input to standard coordinate system                    
!   also computation of initial number of revolutions                   
SUBROUTINE coostp(coox,sysx,refx,nbod,xp,ngp,                     &
     &        sysz,refz,ibar,barin,yp,ngip,bar,aap,eep,ennep)           
  USE massmod
  IMPLICIT NONE 
! number of bodies, baricenter flag                                     
  INTEGER, INTENT(IN) :: nbod,ibar 
!  identifiers of coordinate systems                                    
  character*3 :: coox,cooy,cooz,sysx,sysy,sysz 
  character*6 :: refx,refy,refz 
!  coordinates: input, output                                           
  DOUBLE PRECISION, INTENT(IN) :: xp(6,nbox),ennep(nbox),bar(6),yp(6,nbox) 
  DOUBLE PRECISION zp(6,nbox),zzp(6,nbox),zbp(6,nbox),wp(6,nbox) 
  DOUBLE PRECISION, INTENT(IN) :: barin(6) 
  DOUBLE PRECISION aap(nbox),eep(nbox) 
! revolution counters                                                   
  INTEGER ngip(nbox),ngp(nangx,nbox) 
! loop indexes                                                          
  INTEGER i,j 
  INTEGER npla 
! ==================================================                    
! number of planets                                                     
  npla=nbod-1 
! *******************************************************               
!  integration is in cartesian barycentric, output is always equinoctal 
  cooy='CAR' 
  sysy='BAR' 
  refy=refz 
  cooz='EQU' 
! *******************************************************               
!  coordinate change to integration system; also computation of         
!  initial number of revolutions, for the planets                       
  if(ibar.eq.0)then 
!  (6) computation of initial no. rev; zp is in the output system       
     DO j=1,npla 
        wp(1:6,j)=xp(1:6,j) 
     ENDDO
     call coord(wp,sysx,nbod,coox,refx,zp,sysz,cooz,refz,ennep,bar) 
     call numrev(zp,npla,wp,ngp,zp,ngp,coox,cooz) 
     DO j=1,npla 
        ngip(j)=ngp(1,j)
     ENDDO
!  (7) change to integration coordinate system yp                       
     call coord(xp,sysx,nbod,coox,refx,yp,sysy,cooy,refy,ennep,bar) 
  else 
!  (5) barycenter correction: zzp is cartesian,                         
!     zbp has the bar. correction added                                 
     call coord(xp,sysx,nbod,coox,refx,zzp,'HEL','CAR',refx,ennep,bar)
     DO j=1,nbod-1 
        zbp(1:6,j)=zzp(1:6,j)+barin(1:6) 
     ENDDO
!  (6) computation of initial number of revolutions; zp is in the output
     DO j=1,npla 
        wp(1:6,j)=zbp(1:6,j) 
     ENDDO
     call coord(wp,'HEL',nbod,'CAR',refx,zp,                         &
     &             sysz,cooz,refz,ennep,bar)                            
     call numrev(zp,npla,xp,ngp,zp,ngp,coox,cooz) 
     DO j=1,npla 
        ngip(j)=ngp(1,j)
     ENDDO
!  (7) change to integration coordinate system                          
     call coord(zbp,'HEL',nbod,'CAR',refx,yp,                        &
     &            sysy,cooy,refy,ennep,bar)                             
  endif
  DO j=1,npla 
     aap(j)=zp(1,j) 
     eep(j)=sqrt(zp(2,j)**2+zp(3,j)**2) 
  ENDDO
END SUBROUTINE coostp
! ========================================                              
!   COOSTA case for asteroids                                           
!   reduction of input to standard coordinate system                    
!   also computation of initial number of revolutions                   
SUBROUTINE coosta(coox,sysx,refx,nast,xa,nga,                     &
     &        sysz,refz,nbod,bar,ibar,barin,ya,ngia,aa,ee,ennea)        
  USE massmod
  IMPLICIT NONE 
! number of bodies, baricenter flag, number of massive bodies           
  INTEGER nast,ibar,nbod 
!  identifiers of coordinate systems                                    
  character*3 coox,cooy,cooz,sysx,sysy,sysz 
  character*6 refx,refy,refz 
  DOUBLE PRECISION xa(6,nastx),ennea(nastx),ya(6,nastx) 
  INTEGER nga(nangx,nastx),ngia(nastx) 
  DOUBLE PRECISION za(6,nastx),zza(6,nastx),zba(6,nastx),wa(6,nastx) 
  DOUBLE PRECISION barin(6),bar(6) 
  DOUBLE PRECISION aa(nastx),ee(nastx),sini,tgi2s 
! loop indexes                                                          
  INTEGER i,j 
! =================================================                     
!  integration is in cartesian barycentric, output is always equinoctal 
  cooy='CAR' 
  sysy='BAR' 
  refy=refz 
  cooz='EQU' 
! *******************************************************               
! change to the integration system;                                     
! also computation of the initial number of revolutions                 
  if(ibar.eq.0)then 
!  (9) computation of the initial number of revolutions; za is in the ou
     DO j=1,nast 
        wa(1:6,j)=xa(1:6,j) 
     ENDDO
     call cooast(wa,sysx,nast,coox,refx,bar,nbod,za,sysz,cooz,refz,ennea)
     call numrev(za,nast,wa,nga,za,nga,coox,cooz) 
     DO j=1,nast 
        ngia(j)=nga(1,j)
     ENDDO
!  (10) change to integration coordinate system ya                      
     call cooast(xa,sysx,nast,coox,refx,bar,nbod,ya,sysy,cooy,refy,ennea)
  else 
!  (8) adding the barycenter correction: zza is cartesian, zba is correc
     call cooast(xa,sysx,nast,coox,refx,bar,nbod,zza,'HEL','CAR',   &
     &              refx,ennea)                                         
     DO j=1,nast 
        zba(1:6,j)=zza(1:6,j)+barin(1:6) 
     ENDDO
!  (9) computation of the initial number of revolutions                 
     DO j=1,nast 
        wa(1:6,j)=zba(1:6,j) 
     ENDDO
     call cooast(wa,'HEL',nast,'CAR',refx,bar,nbod,za,sysz,cooz,refz,ennea)
     call numrev(za,nast,xa,nga,za,nga,coox,cooz) 
     DO j=1,nast 
        ngia(j)=nga(1,j)
     ENDDO
!  (10) change to integration coordinate system ya                      
     call cooast(zba,'HEL',nast,'CAR',refx,bar,nbod,ya,sysy,cooy,refy,ennea)
  endif
  DO j=1,nast 
     aa(j)=za(1,j) 
     tgi2s=za(4,j)**2+za(5,j)**2 
     sini=2.d0*sqrt(tgi2s)/(1.d0+tgi2s) 
     ee(j)=sqrt(za(2,j)**2+za(3,j)**2+sini**2) 
  ENDDO
END SUBROUTINE coosta

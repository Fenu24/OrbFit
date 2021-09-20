! ==========================================                            
!  compute size of short periodic perturbations                         
SUBROUTINE epsi(npla,el,enne) 
  USE massmod
  USE fund_const
  IMPLICIT NONE 
!  number of planets 
  INTEGER, PARAMETER ::  nplax=nbox-1
  INTEGER npla  ! actual number of planets                               
  DOUBLE PRECISION el(6,npla),enne(npla) 
! ============================================                          
!  Roy epsilon                                                          
  DOUBLE PRECISION eroy(nplax,nplax),eptot(nplax),epce(nplax) 
!  velocity perturbation                                                
  DOUBLE PRECISION eta(nplax,nplax),etatot(nplax),etce(nplax) 
!  position perturbation                                                
  DOUBLE PRECISION teta(nplax,nplax),tetatot(nplax),tece(nplax) 
! accumulators for Ceres                                                
  DOUBLE PRECISION epct,tecet,etcet,ennec,ac 
!  synodic period                                                       
  DOUBLE PRECISION sper(nplax,nplax),spece(nplax) 
! ============================================                          
!  loop indexes                                                         
  INTEGER i,j 
! output files                                                          
  OPEN(64,file='epsi.out',status='unknown') 
! ============================================                          
!  Roy epsilon               
  WRITE(64,*)'% Roy-Walker epsilon'
  WRITE(64,*)'% perturbed  epsilons tot position eff.'
  DO 1 i=1,npla 
     eptot(i)=0.d0 
!  interior perturbation                                                
     DO j=1,i-1 
        eroy(j,i)=pmu(j)*(1.d0-pmu(j))*(el(1,j)/el(1,i))**2 
        eptot(i)=eptot(i)+eroy(j,i)**2 
     ENDDO
     eroy(i,i)=0.d0 
! exterior perturbation                                                 
     DO j=i+1,npla 
        eroy(j,i)=(gm(j+1)/sm(i))*(el(1,i)/el(1,j))**3 
        eptot(i)=eptot(i)+eroy(j,i)**2 
     ENDDO
     eptot(i)=sqrt(eptot(i)) 
     WRITE(64,100)i,(eroy(j,i),j=1,npla),eptot(i)*el(1,i) 
100  FORMAT(i3,1p,9d10.2,10x,d11.3) 
1 ENDDO
! ============================================                          
!  velocity perturbation
  WRITE(64,*)'% heliocentric eta'
  WRITE(64,*)'% perturbed  eta tot position eff.'
  DO 2 i=1,npla 
     etatot(i)=0.d0 
!  interior perturbation                                                
     DO j=1,i-1 
        eta(j,i)=2.d0*(gm(j+1)/gm(1))*sqrt(el(1,i)/el(1,j)) 
        etatot(i)=etatot(i)+eta(j,i)**2 
     ENDDO
     eta(i,i)=0.d0 
! exterior perturbation                                                 
     DO j=i+1,npla 
        eta(j,i)=2.d0*(gm(j+1)/gm(1))*(el(1,i)/el(1,j))**2          
!        eta(j,i)=eroy(j,i) 
        etatot(i)=etatot(i)+eta(j,i)**2 
     ENDDO
     etatot(i)=sqrt(etatot(i)) 
     WRITE(64,101)i,(eta(j,i),j=1,npla),etatot(i)*el(1,i) 
101  FORMAT(i3,1p,9d10.2,10x,d11.3) 
2 ENDDO
! ============================================                          
!  position perturbation  
  WRITE(64,*)'% baricentric theta'
  WRITE(64,*)'% perturbed  theta tot position eff.'                         
  DO 3 i=1,npla 
     tetatot(i)=0.d0 
!  interior, exterior but not self perturbation                         
     DO j=1,npla 
        IF(j.ne.i)THEN 
           teta(j,i)=(gm(j+1)/gm(1))*(el(1,j)/el(1,i)) 
           tetatot(i)=tetatot(i)+teta(j,i)**2 
        ELSE 
           teta(j,i)=0.d0 
        ENDIF
     ENDDO
     tetatot(i)=sqrt(tetatot(i)) 
     WRITE(64,102)i,(teta(j,i),j=1,npla),tetatot(i)*el(1,i) 
102  FORMAT(i3,1p,9d10.2,10x,d11.3) 
3 ENDDO
! ============================================                          
!  synodic periods  
  WRITE(64,*)'% synodic periods'
  WRITE(64,*)'% perturbed  periods (years)'                                   
  DO 4 i=1,npla 
!  interior perturbation                                                
     DO j=1,i-1 
        sper(j,i)=dpig/(enne(j)-enne(i))/365.25 
     ENDDO
     sper(i,i)=0.d0 
! exterior perturbation                                                 
     DO j=i+1,npla 
        sper(j,i)=dpig/(enne(i)-enne(j))/365.25 
     ENDDO
     WRITE(64,104)i,(sper(j,i),j=1,npla) 
104  FORMAT(i3,9f8.3) 
4 ENDDO
! ===============================================                       
! total perturbation on Ceres                                           
  ac=2.767d0 
  ennec=sqrt(gm(1)/ac**3) 
! accumulators                                                          
  epct=0.d0 
  etcet=0.d0 
  tecet=0.d0 
! interior perturbations                                                
  DO j=1,4 
     epce(j)=pmu(j)*(1.d0-pmu(j))*(el(1,j)/ac)**2 
     epct=epct+epce(j)**2 
     etce(j)=2.d0*(gm(j+1)/gm(1))*sqrt(ac/el(1,j)) 
     etcet=etcet+etce(j)**2 
     tece(j)=(gm(j+1)/gm(1))*(el(1,j)/ac) 
     tecet=tecet+tece(j)**2 
     spece(j)=dpig/(enne(j)-ennec)/365.25 
  ENDDO
  DO j=5,9 
     epce(j)=(gm(j+1)/sm(4))*(ac/el(1,j))**3 
     epct=epct+epce(j)**2 
     etce(j)= 2.d0*(gm(j+1)/gm(1))*(ac/el(1,j))**2
     etcet=etcet+etce(j)**2 
     tece(j)=(gm(j+1)/gm(1))*(el(1,j)/ac) 
     tecet=tecet+tece(j)**2 
     spece(j)=dpig/(ennec-enne(j))/365.25 
  ENDDO
  epct=sqrt(epct) 
  etcet=sqrt(etcet) 
  tecet=sqrt(tecet)
  WRITE(64,*)'% Ceres Roy'
  WRITE(64,103)10,(epce(j),j=1,npla),epct*ac 
103 FORMAT(i3,1p,9d10.2,10x,d11.3) 
  WRITE(64,*)'% Ceres eta'
  WRITE(64,103)10,(etce(j),j=1,npla),etcet*ac
  WRITE(64,*)'% Ceres theta'
  WRITE(64,103)10,(tece(j),j=1,npla),tecet*ac 
  WRITE(64,*)'% Ceres synodic periods'
  WRITE(64,104)10,(spece(j),j=1,npla) 
END SUBROUTINE epsi

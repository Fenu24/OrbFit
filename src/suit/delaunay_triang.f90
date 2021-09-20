!npox,ntrix

! ================================================
!   written by MATTIA DE MICHIELI VITTURI 2002
! ================================================
!       last modified 19/11/2003
! ================================================
SUBROUTINE triangola(npox,ntrix,xRR,yRR,nodiest,noditot,triang,nt)
  USE fund_const
  IMPLICIT none
  INTEGER, INTENT(IN) :: npox,ntrix ! max number of points and triangles
  INTEGER, INTENT(IN) :: nodiest ! number of nodes on the boundary
  INTEGER, INTENT(OUT) :: noditot ! number of nodes
  DOUBLE PRECISION, DIMENSION(npox), INTENT(INOUT) :: xRR,yRR 
! input: nodiest boundary nodes
! output: all noditot nodes
  INTEGER,INTENT(OUT) :: nt ! number of triangles
  INTEGER, INTENT(OUT) :: triang(ntrix,4)
! END INTERFACE
! =========================================================================
  DOUBLE PRECISION xR(0:npox),yR(0:npox)
  INTEGER circmax,seg(ntrix,3,2),triangnew(4),nodo(4)
  INTEGER segnew3(3,2),segmento(npox,npox),nodoopposto(3)
  INTEGER escimax,h,l,i,m,j
  INTEGER intersez,badnodi,cambia(ntrix),vicini(3),quantivicini
  INTEGER quanti2,triangmax(4)
  DOUBLE PRECISION raggio,raggiont,xmax1,xmax2,xmax3,ymax1,ymax2,ymax3
  DOUBLE PRECISION angolo,raggiocirccentr(ntrix)
  DOUBLE PRECISION a11,a12,a21,a22
  DOUBLE PRECISION xcirccentr(ntrix),yverif1,yverif2,yverif3
  DOUBLE PRECISION ycirccentr(ntrix),xverif1,xverif2,xverif3
  DOUBLE PRECISION xbaric(ntrix),ybaric(ntrix)
  DOUBLE PRECISION density(npox),distmin
  DOUBLE PRECISION ITE, TA(2)
  DOUBLE PRECISION xcentro,ycentro,raggiomax,area2,etime
  INTEGER iun
  DOUBLE PRECISION cos1,cos2,sen1,sen2,angolo1,angolo2
  INTEGER exterior
  INTEGER :: count1,count2 ! counters
  DOUBLE PRECISION, PARAMETER :: raggiofix=0.58d0 
  DOUBLE PRECISION, PARAMETER :: epsilon0=1.D-3
  DOUBLE PRECISION, PARAMETER :: eps1=1.D-3
  INTEGER, PARAMETER :: ntmax=5000
  
  DO i = 1,nodiest
     xR(i)=xRR(i)
     yR(i)=yRR(i)
  ENDDO
  xR(nodiest+1)=xR(1)
  yR(nodiest+1)=yR(1)
  xR(0)=xR(nodiest)
  yR(0)=yR(nodiest)
  
  noditot=nodiest
      
  raggiocirccentr=0.d0      
            
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&      
! INITIAL TRIANGULATION
  nt=0
  DO i=1,noditot
     DO j=1,noditot
        segmento(i,j)=0
        IF ((j-i.EQ.1).OR.(j-i.EQ.nodiest-1)) segmento(i,j)=1
     ENDDO
  ENDDO
  
10 DO i=1,noditot-1
     DO j=i+1,noditot
        IF (segmento(i,j).EQ.1) THEN 
           GOTO 5
        ENDIF
     ENDDO
  ENDDO
  
  GOTO 39
  
! THEN we must add a triangle using the two nodes i,j     
5 nt=nt+1
  IF(nt.GT.ntrix)THEN
     write(*,*)'triangola: nt too large ',nt,ntrix
     STOP
  ENDIF
  
  escimax=noditot+1
  DO h=j+1,noditot
     exterior=0
     IF ((h.ne.j).AND.(h.ne.i)) THEN 
        
        nodo(1)=i
        nodo(2)=j
        nodo(3)=h
        nodo(4)=i
        
        DO m=1,3
           a11=yR(nodo(m)-1)-yR(nodo(m))
           a12=yR(nodo(m)+1)-yR(nodo(m))
           a21=xR(nodo(m)-1)-xR(nodo(m))
           a22=xR(nodo(m)+1)-xR(nodo(m))
           sen1=(a11*a22-a12*a21)/(sqrt(a11**2+a21**2)*sqrt(a12**2+a22**2))
           cos1=(a11*a12+a21*a22)/(sqrt(a11**2+a21**2)*sqrt(a12**2+a22**2))
           angolo1=mod(atan2(sen1,cos1)+2.d0*pig,2.d0*pig)
           a11=yR(nodo(m)-1)-yR(nodo(m))
           a12=yR(nodo(m+1))-yR(nodo(m))
           a21=xR(nodo(m)-1)-xR(nodo(m))
           a22=xR(nodo(m+1))-xR(nodo(m))
           sen2=(a11*a22-a12*a21)/(sqrt(a11**2+a21**2)*sqrt(a12**2+a22**2))
           cos2=(a11*a12+a21*a22)/(sqrt(a11**2+a21**2)*sqrt(a12**2+a22**2))
           angolo2=mod(atan2(sen2,cos2)+2.d0*pig,2.d0*pig)
           IF (angolo2.GT.angolo1+1.d-7) THEN
              exterior=1
           ENDIF
        ENDDO
        a11=yR(j)-yR(i)
        a12=-(yR(h)-yR(i))
        a21=-(xR(j)-xR(i))
        a22=xR(h)-xR(i)
        IF (abs(a11*a22-a12*a21)/(sqrt(a11**2+a21**2)&
             &  *sqrt(a12**2+a22**2)).le.epsilon0) exterior=2
        
        IF (exterior.EQ.0) THEN
           triangnew(1)=i
           triangnew(2)=j
           triangnew(3)=h
           triangnew(4)=i
           
           CALL nuovotriangolo2(triangnew,segnew3,&
                & xR(i),yR(i),xR(j),yR(j),xR(h),yR(h),&
                & xcentro,ycentro,raggiont)
           raggio=raggiont
           
           CALL verificaint(npox,ntrix,triang,nt,xR,yR,intersez,xcentro,&
                & ycentro,raggio,badnodi,noditot,triangnew,segmento)
           
           IF ((intersez.EQ.0).AND.((badnodi.LT.escimax)&
                & .OR.((badnodi.EQ.escimax).AND. &
                & (raggiont.lt.raggiocirccentr(nt)*(1.d0+ &
                & epsilon(1.d0)*100.d0))))) THEN
              escimax=badnodi
              
              CALL nuovotriangolo(triangnew,segnew3,&
                   & xR(i),yR(i),xR(j),yR(j),xR(h),yR(h),&
                   & xcentro,ycentro,raggiont,&
                   & xbaric(nt),ybaric(nt))
              raggio=raggiont
              
              triang(nt,1)=triangnew(1)
              triang(nt,2)=triangnew(2)
              triang(nt,3)=triangnew(3)
              xcirccentr(nt)=xcentro
              ycirccentr(nt)=ycentro
              raggiocirccentr(nt)=raggiont
              DO l=1,3
                 seg(nt,l,1)=segnew3(l,1)
                 seg(nt,l,2)=segnew3(l,2)
              ENDDO
           ENDIF
        ENDIF
     ENDIF
     
  ENDDO
!      write(*,*)'nt=',nt,triang(nt,1),triang(nt,2),triang(nt,3)   
  IF(triang(nt,1).EQ.0) THEN
     WRITE(*,*)'triangola: error for segment',i,j,&
          &' triang(nt,1)=',triang(nt,1)
     STOP
  ENDIF
  segmento(triang(nt,1),triang(nt,2))=segmento(triang(nt,1),triang(nt,2))+1
  segmento(triang(nt,2),triang(nt,3))=segmento(triang(nt,2),triang(nt,3))+1
  segmento(triang(nt,1),triang(nt,3))=segmento(triang(nt,1),triang(nt,3))+1
  
  GOTO 10
  !&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&      
  
  !&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  ! OPTIMIZATION OF THE INITIAL TRIANGULATION 
  ! (NEEDED BECAUSE DOMAIN NOT CONVEX)
39 count2 = 0
  
40 count2 = count2+1
  if(count2.GT.2**nt) THEN
     GOTO 36
  ENDIF
  
  DO i=1,nt
     
     CALL cercatriangolivicini(npox,ntrix,triang,nt,i,quantivicini,  &
          &   vicini,seg,nodoopposto)
     
     DO j=1,quantivicini
        
        IF (sqrt((xR(nodoopposto(j))-xcirccentr(i))           &
             &  **2+(yR(nodoopposto(j))-ycirccentr(i))**2)    &
             &  .LT.raggiocirccentr(i)*(1.d0-eps1)) THEN
           
           CALL listaswap(npox,ntrix,triang,xR,yR,xcirccentr,            &
                &   seg,ycirccentr,raggiocirccentr,           &
                &   i,quanti2,cambia,nt,eps1)
           
           CALL swappalista(npox,ntrix,triang,xR,yR,xcirccentr,seg,      &
                &   ycirccentr,raggiocirccentr,quanti2,cambia,&
                &   xbaric,ybaric,eps1)
           GOTO 40
           
        ENDIF
        
     ENDDO
     
  ENDDO

  !&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  ! COMPUTE DENSITY ON THE BOUNDARY
36 DO i=1,noditot
     
     density(i)=0.5d0*sqrt((xR(i)-xR(i-1))**2+(yR(i)-yR(i-1))**2)+ &
          &   0.5d0*sqrt((xR(i)-xR(i+1))**2+(yR(i)-yR(i+1))**2)
     
  ENDDO  
  
  !&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  ! STARTING TRIANGULATION REFINEMENT
  raggiomax=2.d0*raggiofix
  
  !      DO while (raggiomax.ge.raggiofix)
  DO while ((raggiomax.ge.raggiofix).AND.(nt.LT.ntmax-2))
     
     raggiomax=0
     DO i=1,nt
        xverif1=xR(triang(i,1))
        yverif1=yR(triang(i,1))
        xverif2=xR(triang(i,2))
        yverif2=yR(triang(i,2))
        xverif3=xR(triang(i,3))
        yverif3=yR(triang(i,3))
        
        
        distmin=1.d5
        DO j=1,3
           distmin = min(distmin,sqrt((xbaric(i)-xR(triang(i,j)))**2 &
                &  +(ybaric(i)-yR(triang(i,j)))**2))
        ENDDO
        
        area2 = 3.d0*distmin/( density(triang(i,1)) + &
             & density(triang(i,2))+ density(triang(i,3)) )
        
        IF (area2.ge.raggiomax) THEN
           raggiomax=area2
           circmax=i
        ENDIF
     ENDDO
     
     noditot=noditot+1
     
     xR(noditot) = xbaric(circmax)
     yR(noditot) = ybaric(circmax)
     
     density(noditot)=(density(triang(circmax,1))              &
          &  +density(triang(circmax,2))+density(triang(circmax,3)))/3
     
     xmax1=xR(triang(circmax,1))
     ymax1=yR(triang(circmax,1))
     xmax2=xR(triang(circmax,2))
     ymax2=yR(triang(circmax,2))
     xmax3=xR(triang(circmax,3))
     ymax3=yR(triang(circmax,3))
     triangmax(1)=triang(circmax,1)
     triangmax(2)=triang(circmax,2)
     triangmax(3)=triang(circmax,3)
     triangmax(4)=triang(circmax,1)
     
     DO i=circmax,nt-1
        triang(i,1)=triang(i+1,1)
        triang(i,2)=triang(i+1,2)
        triang(i,3)=triang(i+1,3)
        triang(i,4)=triang(i+1,3)
        seg(i,1,1)=seg(i+1,1,1)
        seg(i,1,2)=seg(i+1,1,2)
        seg(i,2,1)=seg(i+1,2,1)
        seg(i,2,2)=seg(i+1,2,2)
        seg(i,3,1)=seg(i+1,3,1)
        seg(i,3,2)=seg(i+1,3,2)
        raggiocirccentr(i)=raggiocirccentr(i+1)
        xcirccentr(i)=xcirccentr(i+1)
        ycirccentr(i)=ycirccentr(i+1)
        xbaric(i)=xbaric(i+1)
        ybaric(i)=ybaric(i+1)
        
     ENDDO
     nt=nt-1
     DO i=1,3
        nt=nt+1
        triangnew(1)=triangmax(i)
        triangnew(2)=triangmax(i+1)
        triangnew(3)=noditot
        triangnew(4)=triangmax(i)
        CALL nuovotriangolo(triangnew,segnew3,      &
             &         xR(triangnew(1)),yR(triangnew(1)),       &
             &         xR(triangnew(2)),yR(triangnew(2)),       &
             &         xR(triangnew(3)),yR(triangnew(3)),       &
             &         xcentro,ycentro,raggiont,                &
             &         xbaric(nt),ybaric(nt))
        
        DO l=1,3
           seg(nt,l,1)=segnew3(l,1)
           seg(nt,l,2)=segnew3(l,2)
        ENDDO
        triang(nt,1)=triangnew(1)
        triang(nt,2)=triangnew(2)
        triang(nt,3)=triangnew(3)
        xcirccentr(nt)=xcentro
        ycirccentr(nt)=ycentro
        raggiocirccentr(nt)=raggiont
        CALL cercatriangolivicini(npox,ntrix,triang,nt,nt,quantivicini, &
             &  vicini,seg,nodoopposto)
        DO j=1,quantivicini
           IF (sqrt((xR(nodoopposto(j))-xcirccentr(nt))**2  &
                &  +(yR(nodoopposto(j))-ycirccentr(nt))**2)   &
                &  .LT.raggiocirccentr(nt)*(1.d0-eps1)) THEN
              CALL listaswap(npox,ntrix,triang,xR,yR,xcirccentr, &
                   &   seg,ycirccentr,raggiocirccentr, &
                   &   nt,quanti2,cambia,nt,eps1)
              CALL swappalista(npox,ntrix,triang,xR,yR,xcirccentr,seg, &
                   &   ycirccentr,raggiocirccentr,quanti2,cambia, &
                   &   xbaric,ybaric,eps1)
           ENDIF
        ENDDO
     ENDDO
  ENDDO
  CALL smoothing(npox,ntrix,eps1,nodiest,noditot,nt,seg,xR,yR,density,triang)
  DO i = 1,noditot
     xRR(i)=xR(i)
     yRR(i)=yR(i)
  ENDDO
  
  RETURN
END SUBROUTINE triangola
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!                       !!!!!!!!!!!
!!!!!!!!!          NUOVOTRIANGOLO   !!!!!!!
!!!!!!!!!                           !!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE nuovotriangolo(triangnew1,segnew,xR1,yR1,xR2,yR2, &
     &          xR3,yR3,xcentronew,ycentronew,raggionew,xbaricx,ybaricx)
  IMPLICIT NONE
  INTEGER, INTENT(INOUT) :: triangnew1(4)
  INTEGER, INTENT(OUT) :: segnew(3,2)
  DOUBLE PRECISION, INTENT(IN) :: xR1,yR1,xR2,yR2,xR3,yR3
  DOUBLE PRECISION, INTENT(OUT)  :: xcentronew,ycentronew,raggionew
  DOUBLE PRECISION, INTENT(OUT)  :: xbaricx,ybaricx
  ! END INTERFACE
  DOUBLE PRECISION :: t,s
  DOUBLE PRECISION :: a11,a12,a21,a22,b1,b2
  DOUBLE PRECISION :: xcirccentr1,xcirccentr2,ycirccentr1,ycirccentr2
  ! =========================================================================
  
  IF (triangnew1(1).GT.triangnew1(2)) THEN
     triangnew1(4)=triangnew1(1)
     triangnew1(1)=triangnew1(2)
     triangnew1(2)=triangnew1(4)
  ENDIF
  IF (triangnew1(1).GT.triangnew1(3)) THEN
     triangnew1(4)=triangnew1(1)
     triangnew1(1)=triangnew1(3)
     triangnew1(3)=triangnew1(4)
  ENDIF
  IF (triangnew1(2).GT.triangnew1(3)) THEN
     triangnew1(4)=triangnew1(2)
     triangnew1(2)=triangnew1(3)
     triangnew1(3)=triangnew1(4)
  ENDIF
  triangnew1(4)=triangnew1(1)
  segnew(1,1)=triangnew1(1)
  segnew(1,2)=triangnew1(2)
  segnew(2,1)=triangnew1(2)
  segnew(2,2)=triangnew1(3) 
  segnew(3,1)=triangnew1(1)
  segnew(3,2)=triangnew1(3)
  a11=yR2-yR1
  a12=-(yR3-yR1)
  a21=-(xR2-xR1)
  a22=xR3-xR1
  b1=0.5*(xR3-xR2)
  b2=0.5*(yR3-yR2)
  
  IF (a11*a22-a12*a21.EQ.0.d0) THEN
     write(*,*) 'nuovotriangolo: zero denominator'
     stop
  ENDIF
  
  t=(b1*a22-a12*b2)/(a11*a22-a12*a21)
  s=(a11*b2-b1*a21)/(a11*a22-a12*a21)
  xcirccentr1=(xR1+xR2)/2.d0+t*(yR2-yR1)
  ycirccentr1=(yR1+yR2)/2.d0-t*(xR2-xR1)
  xcirccentr2=(xR1+xR3)/2.d0+s*(yR3-yR1)
  ycirccentr2=(yR1+yR3)/2.d0-s*(xR3-xR1)
  xcentronew=0.5d0*(xcirccentr1+xcirccentr2)
  ycentronew=0.5d0*(ycirccentr1+ycirccentr2)
  raggionew=sqrt((xcentronew-xR1)**2+(ycentronew-yR1)**2)
  xbaricx=(xR1+xR2+xR3)/3.d0
  ybaricx=(yR1+yR2+yR3)/3.d0
  
  RETURN
END SUBROUTINE nuovotriangolo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!                       !!!!!!!!!!!
!!!!!!!!!     NUOVOTRIANGOLO 2       !!!!!!!!!!!
!!!!!!!!!                         !!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE nuovotriangolo2(triangnew1,segnew,xR1,yR1,xR2,yR2, &
     &    xR3,yR3,xcentronew,ycentronew,raggionew)
  IMPLICIT NONE
  INTEGER, INTENT(INOUT) :: triangnew1(4)
  INTEGER, INTENT(OUT) :: segnew(3,2)
  DOUBLE PRECISION, INTENT(IN) :: xR1,yR1,xR2,yR2,xR3,yR3
  DOUBLE PRECISION, INTENT(OUT) :: xcentronew,ycentronew,raggionew
! END INTERFACE
  DOUBLE PRECISION :: t,s
  DOUBLE PRECISION :: a11,a12,a21,a22,b1,b2
  DOUBLE PRECISION :: xcirccentr1,xcirccentr2,ycirccentr1,ycirccentr2
  ! =========================================================================
  
  
  IF (triangnew1(1).GT.triangnew1(2)) THEN
     triangnew1(4)=triangnew1(1)
     triangnew1(1)=triangnew1(2)
     triangnew1(2)=triangnew1(4)
  ENDIF
  IF (triangnew1(1).GT.triangnew1(3)) THEN
     triangnew1(4)=triangnew1(1)
     triangnew1(1)=triangnew1(3)
     triangnew1(3)=triangnew1(4)
  ENDIF
  IF (triangnew1(2).GT.triangnew1(3)) THEN
     triangnew1(4)=triangnew1(2)
     triangnew1(2)=triangnew1(3)
     triangnew1(3)=triangnew1(4)
  ENDIF
  triangnew1(4)=triangnew1(1)
  segnew(1,1)=triangnew1(1)
  segnew(1,2)=triangnew1(2)
  segnew(2,1)=triangnew1(2)
  segnew(2,2)=triangnew1(3) 
  segnew(3,1)=triangnew1(1)
  segnew(3,2)=triangnew1(3)
  a11=yR2-yR1
  a12=-(yR3-yR1)
  a21=-(xR2-xR1)
  a22=xR3-xR1
  b1=0.5*(xR3-xR2)
  b2=0.5*(yR3-yR2)
  
  IF (a11*a22-a12*a21.EQ.0) THEN
     write(*,*) 'nuovotriangolo2: zero denominator'
     stop
  ENDIF
  
  t=(b1*a22-a12*b2)/(a11*a22-a12*a21)
  s=(a11*b2-b1*a21)/(a11*a22-a12*a21)
  xcirccentr1=(xR1+xR2)/2+t*(yR2-yR1)
  ycirccentr1=(yR1+yR2)/2-t*(xR2-xR1)
  xcirccentr2=(xR1+xR3)/2+s*(yR3-yR1)
  ycirccentr2=(yR1+yR3)/2-s*(xR3-xR1)
  xcentronew=0.5*(xcirccentr1+xcirccentr2)
  ycentronew=0.5*(ycirccentr1+ycirccentr2)
  raggionew=sqrt((xcentronew-xR1)**2+(ycentronew-yR1)**2)
  RETURN
END SUBROUTINE nuovotriangolo2


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!                       !!!!!!!!!!!
!!!!!!!!!      LISTASWAP            !!!!!!!
!!!!!!!!!                       !!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE listaswap(npox,ntrix,triangnew,xRnew,yRnew,xcirccentrnew,seg, &
     &    ycirccentrnew,raggiocirccentrnew,primo,quanti2,cambia,nt,eps1)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: npox,ntrix ! max number of points and triangles
  INTEGER, INTENT(IN) :: triangnew(ntrix,4),primo
  INTEGER, INTENT(OUT):: quanti2,cambia(ntrix)
  INTEGER, INTENT(IN) :: nt,seg(ntrix,3,2)
  DOUBLE PRECISION, INTENT(IN) :: xRnew(0:npox),yRnew(0:npox)
  DOUBLE PRECISION, INTENT(IN) :: xcirccentrnew(ntrix)
  DOUBLE PRECISION, INTENT(IN) :: raggiocirccentrnew(ntrix),&
       & ycirccentrnew(ntrix)
  DOUBLE PRECISION, INTENT(IN) :: eps1
  ! END INTERFACE
  INTEGER :: i,j ! loop indexes
  INTEGER :: esiste,quantivicini,vicini(3),nodoopposto(3)
! =========================================================================    
  quanti2=1
  cambia(1)=primo
900 CALL cercatriangolivicini(npox,ntrix,triangnew,nt,cambia(1),quantivicini, &
         &    vicini,seg,nodoopposto)
  cambia(quanti2+1)=cambia(1)
  DO j=1,quanti2
     cambia(j)=cambia(j+1)
  ENDDO
  DO i=1,quantivicini
     esiste=0
     DO j=1,quanti2
        IF (vicini(i).EQ.cambia(j)) esiste=1
     ENDDO
     IF ((sqrt((xRnew(nodoopposto(i))-xcirccentrnew(primo))**2 &
          & +(yRnew(nodoopposto(i))-ycirccentrnew(primo))**2) &
          & .LT.raggiocirccentrnew(primo)*(1.d0-eps1)).AND.(esiste.EQ.0))THEN
        quanti2=quanti2+1
        DO j=2,quanti2
           cambia(quanti2+2-j)=cambia(quanti2+1-j)
        ENDDO
        cambia(1)=vicini(i)
     ENDIF
  ENDDO
  IF (cambia(1).ne.primo) GOTO 900
  RETURN
END SUBROUTINE listaswap


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!                       !!!!!!!!!!!
!!!!!!!!!      SWAPPALISTA           !!!!!!
!!!!!!!!!                       !!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     

SUBROUTINE swappalista(npox,ntrix,triangnew,xRnew,yRnew,xcirccentrnew,seg, &
     &    ycirccentrnew,raggiocirccentrnew,quanti,cambiain,xbaric,ybaric, &
     &    eps1)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: npox,ntrix ! max number of points and triangles
  INTEGER, INTENT(INOUT) :: triangnew(ntrix,4)
  INTEGER, INTENT(IN) :: quanti,cambiain(ntrix)
  INTEGER, INTENT(INOUT) :: seg(ntrix,3,2)
  DOUBLE PRECISION, INTENT(IN) :: xRnew(0:npox),yRnew(0:npox)
  DOUBLE PRECISION,INTENT(INOUT) :: xcirccentrnew(ntrix)
  DOUBLE PRECISION,INTENT(INOUT):: raggiocirccentrnew(ntrix),&
       & ycirccentrnew(ntrix)
  DOUBLE PRECISION,INTENT(INOUT):: xbaric(ntrix),ybaric(ntrix)
  DOUBLE PRECISION,INTENT(IN) :: eps1
  ! END INTERFACE
  INTEGER :: cambia(ntrix)
  INTEGER :: h,i,j,l ! loop indexes
  INTEGER :: count,ntg ! counters
  INTEGER :: nodo2,nodo4
  INTEGER :: segnew4(3,2),segnew5(3,2),primo,triangnew4(4)
  DOUBLE PRECISION :: xcirco,ycirco,raggiocirco
! =========================================================================
  
  count = 0
  DO h = 1,quanti
     cambia(h) = cambiain(h)   
  ENDDO
  
  ntg=0
800 primo=cambia(1)
  
  count = count+1
  if(count.GT.2**quanti) THEN
     !       write(*,*) 'there is a loop'
     GOTO 313
  ENDIF
  
  IF (ntg.LT.quanti) THEN
     DO i=2,quanti
        DO j=1,3
           DO h=1,3
              IF ((seg(primo,j,1).EQ.seg(cambia(i),h,1)).AND.&
                   &  (seg(primo,j,2).EQ.seg(cambia(i),h,2))) THEN
                 IF (h.EQ.1) THEN
                    nodo4=triangnew(cambia(i),3)
                 ELSEIF (h.EQ.2) THEN
                    nodo4=triangnew(cambia(i),1)
                 ELSE
                    nodo4=triangnew(cambia(i),2)
                 ENDIF
                 IF (j.EQ.1) THEN
                    nodo2=triangnew(primo,3)
                 ELSEIF (j.EQ.2) THEN
                    nodo2=triangnew(primo,1)
                 ELSE
                    nodo2=triangnew(primo,2)
                 ENDIF
                 IF (sqrt((xRnew(nodo4)-xcirccentrnew(primo))**2     &
                      &  +(yRnew(nodo4)-ycirccentrnew(primo))**2).LT. &
                      &  raggiocirccentrnew(primo)*(1.d0-eps1)) THEN
                    triangnew4(1)=nodo2
                    triangnew4(2)=seg(primo,j,1)
                    triangnew4(3)=nodo4
                    CALL nuovotriangolo2(triangnew4,segnew4, &
                         &  xRnew(nodo2),yRnew(nodo2),          &
                         &  xRnew(seg(primo,j,1)),yRnew(seg(primo,j,1)),&
                         &  xRnew(nodo4),yRnew(nodo4),xcirco,ycirco, &
                         &  raggiocirco)
                    triangnew(primo,1)=triangnew4(1)
                    triangnew(primo,2)=triangnew4(2)
                    triangnew(primo,3)=triangnew4(3)
                    triangnew4(1)=nodo4
                    triangnew4(2)=nodo2
                    triangnew4(3)=seg(primo,j,2)
                    xcirccentrnew(primo)=xcirco
                    ycirccentrnew(primo)=ycirco
                    raggiocirccentrnew(primo)=raggiocirco
                    CALL nuovotriangolo2(triangnew4,segnew5, &
                         &  xRnew(nodo4),yRnew(nodo4),          &
                         &  xRnew(nodo2),yRnew(nodo2),          &
                         &  xRnew(seg(primo,j,2)),yRnew(seg(primo,j,2)),&
                         &  xcirco,ycirco,raggiocirco)
                    triangnew(cambia(i),1)=triangnew4(1)
                    triangnew(cambia(i),2)=triangnew4(2)
                    triangnew(cambia(i),3)=triangnew4(3)
                    xcirccentrnew(cambia(i))=xcirco
                    ycirccentrnew(cambia(i))=ycirco
                    raggiocirccentrnew(cambia(i))=raggiocirco
                    DO l=1,3
                       seg(primo,l,1)=segnew4(l,1)
                       seg(primo,l,2)=segnew4(l,2)
                    ENDDO
                    DO l=1,3
                       seg(cambia(i),l,1)=segnew5(l,1)
                       seg(cambia(i),l,2)=segnew5(l,2)
                    ENDDO
                    ntg=0
                    cambia(quanti+1)=cambia(1)
                    cambia(quanti+2)=cambia(i)
                    DO l=1,i-2
                       cambia(l)=cambia(l+1)
                    ENDDO
                    DO l=i-1,quanti
                       cambia(l)=cambia(l+2)
                    ENDDO
                    GOTO 800
                 ENDIF
              ENDIF
           ENDDO
        ENDDO
     ENDDO
     cambia(quanti+1)=cambia(1)
     DO l=1,quanti
        cambia(l)=cambia(l+1)
     ENDDO
     ntg=ntg+1
     GOTO 800
  ENDIF
  
313 CONTINUE
  
  DO i=1,quanti
     triangnew4(1)=triangnew(cambia(i),1)
     triangnew4(2)=triangnew(cambia(i),2)
     triangnew4(3)=triangnew(cambia(i),3)
     CALL nuovotriangolo(triangnew4,segnew5,xRnew(triangnew4(1)), &
          &  yRnew(triangnew4(1)),xRnew(triangnew4(2)),yRnew(triangnew4(2)) &
          &  ,xRnew(triangnew4(3)),yRnew(triangnew4(3)),   &
          &  xcirco,ycirco,raggiocirco,                    &
          &  xbaric(cambia(i)),ybaric(cambia(i)))
     
     IF ((triangnew4(1).EQ.2).AND.(triangnew4(2).EQ.3).AND. &
          & (triangnew4(3).EQ.20)) THEN
        
     ENDIF
     
  ENDDO
  
  RETURN
END SUBROUTINE swappalista

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!                       !!!!!!!!!!!
!!!!!!!!!   CERCATRIANGOLIVICINI     !!!!!!
!!!!!!!!!                       !!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE cercatriangolivicini(npox,ntrix,triang,nt,quale,quantivicini,&
     &  vicini,seg,nodoopposto)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: npox,ntrix ! max number of points and triangles
  INTEGER, INTENT(IN) :: triang(ntrix,4)
  INTEGER, INTENT(IN) :: nt,quale,seg(ntrix,3,2)
  INTEGER, INTENT(OUT) :: vicini(3),nodoopposto(3),quantivicini
! END INTERFACE
! =========================================================================
  INTEGER h,i,j
  
  quantivicini=0
  DO 1000 i=1,nt
     DO j=1,3
        DO h=1,3
           IF ((seg(quale,j,1).EQ.seg(i,h,1)).AND.               &
                &(seg(quale,j,2).EQ.seg(i,h,2)).AND.(i.ne.quale))THEN
              quantivicini=quantivicini+1
              IF (h.EQ.1) THEN
                 nodoopposto(quantivicini)=triang(i,3)
              ELSEIF (h.EQ.2) THEN
                 nodoopposto(quantivicini)=triang(i,1)
              ELSE
                 nodoopposto(quantivicini)=triang(i,2)
              ENDIF
              vicini(quantivicini)=i
              GOTO 1000
           ENDIF
        ENDDO
     ENDDO
1000 ENDDO
  RETURN
END SUBROUTINE cercatriangolivicini

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!                       !!!!!!!!!!!
!!!!!!!!!       VERIFICAINT         !!!!!!!
!!!!!!!!!                          !!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE verificaint(npox,ntrix,triang,nt,xR5,yR5,intersez2,&
     & xcentro,ycentro,raggio,nodibad,noditot,triangnew,segmento)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: npox,ntrix ! max number of points and triangles
  INTEGER, INTENT(IN) :: nt,noditot
  INTEGER, INTENT(OUT) :: intersez2,nodibad
  INTEGER, INTENT(IN) :: triangnew(4),triang(ntrix,4)
  DOUBLE PRECISION, INTENT(IN) :: xR5(0:npox),yR5(0:npox)
  INTEGER, INTENT(IN) :: segmento(npox,npox)
  ! END INTERFACE
  DOUBLE PRECISION :: x1,x2,x3,x4,ba11,ba12,ba21,ba22,bb1,bb2
  DOUBLE PRECISION :: y1,y2,y3,y4,lambda1,gamma1
  INTEGER :: i,j,h ! loop indexes
  DOUBLE PRECISION :: xcentro,ycentro,raggio
! =========================================================================
  
  intersez2=0
  DO i=1,noditot-1
     DO j=i,noditot
        IF ((segmento(i,j).ne.0).AND.&
             &(((i.ne.triangnew(1)).OR.(i.ne.triangnew(3))).OR.&
             &((j.ne.triangnew(3)).OR.(j.ne.triangnew(1))))) THEN
           x1=xR5(triangnew(1))
           x2=xR5(triangnew(3))
           x3=xR5(i)
           x4=xR5(j)
           y1=yR5(triangnew(1))
           y2=yR5(triangnew(3))
           y3=yR5(i)
           y4=yR5(j)
           ba11=x1-x2
           ba12=x4-x3
           ba21=y1-y2
           ba22=y4-y3
           IF (abs(ba11*ba22-ba12*ba21)/(sqrt(ba11**2+ba21**2)&
                & *sqrt(ba12**2+ba22**2)).ge.1.d-4) THEN
              bb1=x4-x2
              bb2=y4-y2
              lambda1=(bb1*ba22-bb2*ba12)/(ba11*ba22-ba12*ba21)
              gamma1=(ba11*bb2-ba21*bb1)/(ba11*ba22-ba12*ba21)
           ELSE
              lambda1=2
           ENDIF
           IF (((lambda1.GT.1.d-4).AND.(lambda1.LT.1-1.d-4)).AND.&
                & ((gamma1.GT.1.d-4).AND.(gamma1.LT.1-1.d-4))) THEN
              
              intersez2=1
              GOTO 2000
           ENDIF
        ENDIF
        IF ((segmento(i,j).ne.0).AND.&
             &(((i.ne.triangnew(2)).OR.(i.ne.triangnew(3))).OR.&
             &((j.ne.triangnew(3)).OR.(j.ne.triangnew(2))))) THEN
           x1=xR5(triangnew(2))
           x2=xR5(triangnew(3))
           x3=xR5(i)
           x4=xR5(j)
           y1=yR5(triangnew(2))
           y2=yR5(triangnew(3))
           y3=yR5(i)
           y4=yR5(j)
           ba11=x1-x2
           ba12=x4-x3
           ba21=y1-y2
           ba22=y4-y3
           IF (abs(ba11*ba22-ba12*ba21)/(sqrt(ba11**2+ba21**2) &
                & *sqrt(ba12**2+ba22**2)).ge.1.d-4) THEN
              bb1=x4-x2
              bb2=y4-y2
              lambda1=(bb1*ba22-bb2*ba12)/(ba11*ba22-ba12*ba21)
              gamma1=(ba11*bb2-ba21*bb1)/(ba11*ba22-ba12*ba21)
           ELSE
              lambda1=2
           ENDIF
           IF (((lambda1.GT.1.d-4).AND.(lambda1.LT.1-1.d-4)).AND.&
                &((gamma1.GT.1.d-4).AND.(gamma1.LT.1-1.d-4))) THEN
              intersez2=1
              
              GOTO 2000
           ENDIF
        ENDIF
     ENDDO
  ENDDO
  IF (nt.ge.2) THEN
     DO i=1,nt-1
        IF (((triang(i,1).EQ.triangnew(1)).AND.&
             &(triang(i,2).EQ.triangnew(2))).AND.&
             &(triang(i,3).EQ.triangnew(3))) THEN
           intersez2=1
           GOTO 2000
        ENDIF
     ENDDO
     
  ENDIF
  
  DO i =1,noditot
     IF((i.ne.triangnew(1)).AND.(i.ne.triangnew(2)).AND.&
          &(i.ne.triangnew(3))) THEN
        x1 = xR5(triangnew(1))
        x2 = xR5(triangnew(2))
        x3 = xR5(triangnew(3))
        x4 = xR5(i)
        y1 = yR5(triangnew(1))
        y2 = yR5(triangnew(2))
        y3 = yR5(triangnew(3))
        y4 = yR5(i)
        
        ba11=x2-x1
        ba12=x3-x1
        ba21=y2-y1
        ba22=y3-y1
        IF (abs(ba11*ba22-ba12*ba21)/(sqrt(ba11**2+ba21**2)&
             & *sqrt(ba12**2+ba22**2)).ge.1.d-4) THEN
           bb1=x4-x1
           bb2=y4-y1
           lambda1=(bb1*ba22-bb2*ba12)/(ba11*ba22-ba12*ba21)
           gamma1=(ba11*bb2-ba21*bb1)/(ba11*ba22-ba12*ba21)
        ELSE
           intersez2=1
           GOTO 2000
        ENDIF
        
        IF (((lambda1.GT.-1.d-4).AND.(lambda1.LT.1+1.d-4)).AND.&
             & ((gamma1.GT.-1.d-4).AND.(gamma1.LT.1+1.d-4)).AND.&
             &(lambda1+gamma1.le.1.E0+1.d-4)) THEN
           intersez2=1
           
           GOTO 2000
        ENDIF
     ENDIF
  ENDDO
  
  nodibad=0
  DO i=1,noditot
     IF (((i.EQ.triangnew(1)).OR.(i.EQ.triangnew(2)))&
          &.OR.(i.EQ.triangnew(3))) THEN
        nodibad=nodibad+1
     ELSE
        IF ((xcentro-xR5(i))**2+(ycentro-yR5(i))**2.le.raggio**2) &
             &THEN
           nodibad=nodibad+1
        ENDIF
     ENDIF
  ENDDO
2000 CONTINUE
  
  RETURN     
END SUBROUTINE verificaint

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!			     !!!!!!!!!!!
!!!!!!!!!   CERCANODIVICINI          !!!!!!!!!!!
!!!!!!!!!			     !!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE cercanodivicini(npox,ntrix,quale,seg,xR,yR,density,nt)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: npox,ntrix ! max number of points and triangles
  INTEGER,INTENT(IN) :: nt
  INTEGER,INTENT(IN) :: quale,seg(2*npox,3,2)
  DOUBLE PRECISION,INTENT(INOUT) :: xR(0:npox),yR(0:npox)
  DOUBLE PRECISION,INTENT(INOUT) :: density(npox)
  ! end interface
  DOUBLE PRECISION :: xnew,ynew,sum
  DOUBLE PRECISION :: densityNew
  INTEGER :: i,j,h
  INTEGER :: vertice(npox),newlist(npox),vertmax
  INTEGER :: indmax
! ==========================================================
  
  h=0
  DO i=1,nt
     DO j=1,3
        IF (seg(i,j,1).EQ.quale) THEN
           h=h+1
           vertice(h)=seg(i,j,2)
        ELSEIF (seg(i,j,2).EQ.quale) THEN
           h=h+1
           vertice(h)=seg(i,j,1)
        ENDIF
     ENDDO
  ENDDO
  
  DO i=1,h
     vertmax=0
     DO j=i,h
        IF (vertice(j).GT.vertmax) THEN
           vertmax=vertice(j)
           indmax=j
        ENDIF
     ENDDO
     vertice(indmax)=vertice(i)
     vertice(i)=vertmax
  ENDDO
  newlist(1)=vertice(1)
  j=1
  DO i=2,h
     IF (vertice(i).ne.newlist(j)) THEN
        j=j+1
        newlist(j)=vertice(i)
     ENDIF
  ENDDO
  xnew=0.
  ynew=0.
  sum=0.
  densityNew=0.
  DO i=1,j
     xnew=xnew+1/density(newlist(i))*xR(newlist(i))
     ynew=ynew+1/density(newlist(i))*yR(newlist(i))
     sum=sum+1/density(newlist(i))
     
  ENDDO
  xnew=xnew/sum
  ynew=ynew/sum
  densityNew=j/sum
  xR(quale)=xnew
  yR(quale)=ynew
  density(quale)=densityNew
END SUBROUTINE cercanodivicini

! *************************************************
! SMOOTHING and OPTIMIZATION of the TRIANGULATION
! *************************************************
SUBROUTINE smoothing(npox,ntrix,eps,nodiest,noditot,nt,seg, &
     & xR,yR,density,triang)
  IMPLICIT none
  INTEGER, INTENT(IN) :: npox,ntrix ! max number of points and triangles
  DOUBLE PRECISION,INTENT(IN) :: eps
  INTEGER,INTENT(IN) :: nodiest,noditot,nt
  INTEGER,INTENT(INOUT) :: seg(ntrix,3,2)
  DOUBLE PRECISION,INTENT(INOUT) :: xR(0:npox),yR(0:npox)
  DOUBLE PRECISION,INTENT(INOUT) :: density(npox)
  INTEGER,INTENT(INOUT) :: triang(ntrix,4)
! END INTERFACE
  INTEGER :: i,j,l,h ! loop indexes     
  INTEGER :: count ! counter
  INTEGER :: triangnew(4)
  INTEGER :: segnew(3,2) !only for CALL od nuovotriangolo, not needed here
  DOUBLE PRECISION :: xcirccentr(2*npox),ycirccentr(2*npox), &
       & raggiocirccentr(2*npox)
  DOUBLE PRECISION :: xcentro,ycentro,raggiont
  INTEGER :: nodoopposto(3),quantivicini,vicini(3),quanti2,cambia(2*npox)
  DOUBLE PRECISION :: xbaric(2*npox),ybaric(2*npox)
! ================================================================
  
! move the nodes in the barycenter of the neighbor polygon
  DO h=1,2
     DO 350 l=nodiest+1,noditot
        CALL cercanodivicini(npox,ntrix,l,seg,xR,yR,density,nt)
350  ENDDO
     DO 360 l=noditot,nodiest+1,-1
        CALL cercanodivicini(npox,ntrix,l,seg,xR,yR,density,nt)
360  ENDDO
  ENDDO
  
  DO i=1,nt
     
     triangnew(1)=triang(i,1)
     triangnew(2)=triang(i,2)
     triangnew(3)=triang(i,3)
     triangnew(4)=triang(i,1)
     CALL nuovotriangolo(triangnew,segnew, &
          &  xR(triangnew(1)),yR(triangnew(1)), &
          &  xR(triangnew(2)),yR(triangnew(2)), &
          &  xR(triangnew(3)),yR(triangnew(3)), &
          &  xcentro,ycentro,raggiont, &
          &  xbaric(i),ybaric(i))
     xcirccentr(i)=xcentro
     ycirccentr(i)=ycentro
     raggiocirccentr(i)=raggiont
     
  ENDDO
  
  count = 0
  
400 count = count+1
  !     if (count.GT.2**nt) THEN
  IF(count.gt.2*nt)THEN
     GOTO 133
  ENDIF
  
  DO i=1,nt
     CALL cercatriangolivicini(npox,ntrix,triang,nt,i,quantivicini, &
          &  vicini,seg,nodoopposto)
     
     DO j=1,quantivicini
        IF (sqrt((xR(nodoopposto(j))-xcirccentr(nt))**2  &
             &  +(yR(nodoopposto(j))-ycirccentr(nt))**2)   &
             &  .LT.raggiocirccentr(nt)*(1.d0-eps)) THEN
           CALL listaswap(npox,ntrix,triang,xR,yR,xcirccentr, &
                &   seg,ycirccentr,raggiocirccentr, &
                &   nt,quanti2,cambia,nt,eps)
           CALL swappalista(npox,ntrix,triang,xR,yR,xcirccentr,seg, &
                &   ycirccentr,raggiocirccentr,quanti2,cambia, &
                &   xbaric,ybaric,eps)
           
           GOTO 400
           
        ENDIF
     ENDDO
  ENDDO
  
133 CONTINUE 
  
END SUBROUTINE smoothing

! =============================
! READING TRIANGULATION FILES
! =========================================================================
SUBROUTINE outtriang(trifil,noditot,nt,xR,yR,triang)
  CHARACTER*(*),INTENT(IN) :: trifil
  INTEGER, INTENT(IN) :: noditot,nt,triang(nt,3)   !triang(nt,4)
  DOUBLE PRECISION, INTENT(IN) :: xR(noditot),yR(noditot)
  ! END INTERFACE
! =========================================================================
  DOUBLE PRECISION :: xn1,xn2,xn3,yn1,yn2,yn3
  INTEGER :: i
  INTEGER :: iun ! unit number      
  CALL filopn(iun,trifil,'unknown')  
  !       write total number of nodes and of triangles
  write(iun,100) noditot,nt, 0.d0 !,nfunc
100 FORMAT(i5,8x,i5,8x,i5)
  DO i=1,noditot
     write(iun,101) xR(i),yR(i), 0.d0
  ENDDO
101 FORMAT(f17.12,3x,f17.12,3x,f3.1)
  DO i=1,nt
     xn1=xR(triang(i,1))
     yn1=yR(triang(i,1))
     xn2=xR(triang(i,2))
     yn2=yR(triang(i,2))
     xn3=xR(triang(i,3))
     yn3=yR(triang(i,3))    
     IF ((xn2-xn1)*(yn3-yn1)-(xn3-xn1)*(yn2-yn1).LT.0) THEN
        write(iun,100) triang(i,1),triang(i,3),triang(i,2)
     ELSE
        write(iun,100) triang(i,1),triang(i,2),triang(i,3)
     ENDIF
  ENDDO
  CALL filclo(iun,' ')
END SUBROUTINE outtriang

SUBROUTINE intriang(npox,ntrix,trifil,noditot,nt,xR,yR,triang)
!  USE triangles
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: npox,ntrix ! max number of points and triangles
  CHARACTER*(*),INTENT(IN) :: trifil ! file name
  INTEGER, INTENT(OUT) :: noditot,nt ! total number of nodes and triangles
  INTEGER, INTENT(OUT) :: triang(3,ntrix) ! incidence relations between nodes
  DOUBLE PRECISION, INTENT(OUT) :: xR(npox),yR(npox) ! node coords
  ! END INTERFACE
! =========================================================================
  INTEGER :: i ! loop index
  INTEGER :: iun ! unit number      
  CALL filopn(iun,trifil,'old')  
  ! READ total number of nodes and of triangles
  READ(iun,100) noditot,nt
100 FORMAT(i5,8x,i5,8x,i5)
  ! READ  node coordinates
  DO i=1,noditot
     READ(iun,101) xR(i),yR(i)
  enddo
101 FORMAT(f17.12,3x,f17.12,3x,f3.1)
  ! READ incidence relations between nodes
  DO  i=1,nt
     READ(iun,100) triang(1,i),triang(2,i),triang(3,i)
  ENDDO
  CALL filclo(iun,' ')
END SUBROUTINE intriang

! =================================================================
! MERGE_TRIANG
! merges two triangulations into one, with list of nodes and list
! of triangles
! =================================================================
SUBROUTINE merge_triang(npox,ntrix,frv1,rdot1,npotot1,triang1,ntri1, &
     & frv2,rdot2,npotot2,triang2,ntri2, &
     & logrv,rdot,npotot,triang,ntri)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: npox,ntrix ! max number of points and triangles
! triangulations: node coordinates (log10(r), rdot)
  DOUBLE PRECISION,INTENT(IN),DIMENSION(npox) :: frv1,frv2,rdot1,rdot2
  DOUBLE PRECISION,INTENT(OUT),DIMENSION(npox) :: logrv,rdot
! triangle indexes
  INTEGER,INTENT(IN) :: triang1(ntrix,4),triang2(ntrix,4)
  INTEGER,INTENT(OUT) :: triang(ntrix,4)
! total number of points
  INTEGER,INTENT(IN) :: npotot1,npotot2 
  INTEGER,INTENT(OUT) :: npotot
! total number of triangles
  INTEGER,INTENT(IN) :: ntri1,ntri2
  INTEGER,INTENT(OUT) :: ntri
! =========== end interface =======================================
! to count the node repetitions
  INTEGER nrep,repnod2list(npotot2),repnod1list(npotot2)
! loop indexes
  INTEGER i,j,h,m,n
! =================================================================

! initialization
  nrep=0
  
! nodes and triangles of the first list remain unchanged
  DO i = 1,npotot1
     logrv(i) = frv1(i)
     rdot(i) = rdot1(i)
     
  ENDDO
  DO i=1,ntri1
     DO n = 1,3
        triang(i,n) = triang1(i,n)
     ENDDO
  ENDDO
  
  DO h = 1,npotot2
     DO i = 1,npotot1
        
        IF (abs((frv2(h) - frv1(i))/frv2(h)).lt.0.d0) THEN
           IF (abs((rdot2(h) - rdot1(i))).lt.0.d0) THEN
                  ! updating counters
              nrep = nrep + 1       ! count node repetitions
! number of the repeated node in 2nd list 
              repnod2list(nrep) = h
! number of the repeated node in 1st list
              repnod1list(nrep) = i
!                  WRITE(*,*)' doppio ', nrep,h,i                  
              GOTO 10
           ELSE
! add to the list as a new node
              logrv(h+npotot1-nrep) = frv2(h)
              rdot(h+npotot1-nrep) = rdot2(h)
           ENDIF
        ELSE
! add to the list as a new node
           logrv(h+npotot1-nrep) = frv2(h)
           rdot(h+npotot1-nrep) = rdot2(h)
        ENDIF
            
     ENDDO
         
10   CONTINUE
  ENDDO

! LOOP ON TRIANG      
  DO i = 1,ntri2
     DO n = 1,3
! no duplicates
        IF(nrep.eq.0)THEN
           triang(i+ntri1,n) = triang2(i,n)+npotot1
        ELSE
! there are duplicates            
           IF (triang2(i,n).lt.repnod2list(1)) THEN
              triang(i+ntri1,n) = triang2(i,n)+npotot1
           ENDIF
           
           DO m =1,nrep-1
              IF (triang2(i,n).eq.repnod2list(m)) THEN
                 triang(i+ntri1,n) = repnod1list(m)
              ELSEIF ((triang2(i,n).gt.repnod2list(m)).and. &
                   &             (triang2(i,n).lt.repnod2list(m+1))) THEN
                 triang(i+ntri1,n) = triang2(i,n)+npotot1-m
              ENDIF
           ENDDO
            
           IF (triang2(i,n).eq.repnod2list(nrep)) THEN
              triang(i+ntri1,n) = repnod1list(nrep)
           ELSEIF (triang2(i,n).gt.repnod2list(nrep)) THEN
              triang(i+ntri1,n) = triang2(i,n)+npotot1-nrep
           ENDIF
        ENDIF
     ENDDO
  ENDDO

! total number of triangles
  ntri = ntri1+ntri2
  npotot = npotot1 + npotot2 - nrep      

! writing on output file
  OPEN(unit=2,file='merge_tri.dat',status='unknown') 
  WRITE(2,100) npotot,ntri,0
  DO i=1,npotot 
     WRITE(2,101) logrv(i),rdot(i),0.d0
  ENDDO
  DO j=1,ntri
     WRITE(2,100) triang(j,1),triang(j,2),triang(j,3)
  ENDDO
  
100 FORMAT(i5,8x,i5,8x,i5)
101 FORMAT(f10.5,3x,f10.5,3x,f10.5)
  
  CLOSE(2)
  
  RETURN
END SUBROUTINE merge_triang

! ================================================
!   written by MATTIA DE MICHIELI VITTURI 2002
!       =======================================
!       =======================================
! last modified GFG 7/7/2003
SUBROUTINE uniform(xRR,yRR,nodiinit,nodifin,quali,npoxmul)
  implicit none
  INTEGER, INTENT(IN) :: nodiinit,nodifin,npoxmul
  INTEGER, INTENT(OUT) :: quali(npoxmul) 
  DOUBLE PRECISION ,INTENT(IN) :: xRR(nodiinit),yRR(nodiinit)      
!     END INTERFACE
  DOUBLE PRECISION :: lung(0:nodiinit),densita(nodiinit),lung2
  DOUBLE PRECISION :: lungtot,distmin,distanza(nodiinit)
  INTEGER :: nodoelim
  INTEGER :: quanti,i,j,h,dist
  INTEGER :: ninit,nfin
  DOUBLE PRECISION :: intervallo
  DOUBLE PRECISION, PARAMETER ::opt=1.d0
! ================================================      
  ninit = nodiinit
  nfin = nodifin
  quanti=ninit
  
  DO i=1,quanti
     quali(i)=i
  ENDDO
      
!      IF (ninit.eq.2) THEN
!         write(*,*)'xRR,yRR',xRR,yRR
!         write(*,*)'ninit,nfin',ninit,nfin
!      ENDIF
  IF(ninit.lt.nfin) THEN
     write(*,*)'ninit,nfin',ninit,nfin
     write(*,*)'ninit < nfin: forcing ninit = nfin'
     nfin = ninit
  ENDIF
  IF(nfin.lt.2) THEN
     IF(ninit.ge.2) THEN
        write(*,*)'nfin=',nfin,'  forcing nfin = 2'
        nfin = 2
     ELSE
        WRITE(*,*)'ERROR: ninit,nfin',ninit,nfin
     ENDIF
  ENDIF
! ======================================
! =======  CALCOLA ASC.CURV.  ==========
! ======================================
  lung(0)=0.d0
  lungtot=0.d0
  DO i = 1,ninit-1
     lung(i)=sqrt((xRR(i)-xRR(i+1))**2+(yRR(i)-yRR(i+1))**2)
     lungtot=lungtot+lung(i)
  ENDDO
  lung(ninit)=0.
  intervallo=(lungtot/(nfin-1))  
! ======================================
! ======  CALCOLA DIST E DENS  =========
! ======================================
  densita(1)=0.
  distanza(1)=0.
  DO i = 2,ninit-1
     densita(i)=opt*min(lung(i),lung(i-1))+(1-opt)*(lung(i)+lung(i-1))
     distanza(i)=lungtot
     lung2=0.
     DO j= 1,i
        lung2=lung2+lung(j-1)
     ENDDO
     DO j= 1,nfin-2
        distanza(i)=min(distanza(i),abs(lung2-lungtot*j/(nfin-1)))
     ENDDO
  ENDDO
  densita(ninit)=0.
  distanza(ninit)=0.  
! ======================================
! =========  ELIMINA I PUNTI  ==========
! ======================================
  DO h=1,ninit-nfin
     distmin=2*lungtot
     DO i=2,quanti-1
        IF (densita(i)/(1+sqrt(distanza(i))).lt.distmin) THEN
           distmin=densita(i)/(1+sqrt(distanza(i)))
           nodoelim=i   !indice punto da eliminare
        ENDIF
     ENDDO
     IF (nodoelim.eq.2) THEN
        lung(1)=lung(1)+lung(2)
        densita(3)=opt*min(lung(1),lung(3))+&
             &  (1-opt)*(lung(1)+lung(3))
        DO j=nodoelim,quanti-2
           quali(j)=quali(j+1)
           densita(j)=densita(j+1)
           lung(j)=lung(j+1)
           distanza(j)=distanza(j+1) 
        ENDDO
     ELSEIF (nodoelim.eq.quanti-1) THEN
        lung(quanti-2)=lung(quanti-2)+lung(quanti-1)
        densita(quanti-2)=opt*(min(lung(quanti-3),lung(quanti-2)))+&
             &  (1-opt)*(lung(quanti-3)+lung(quanti-2))
        quali(quanti-1)=quali(quanti)
     ELSE
        lung(nodoelim-1)=lung(nodoelim-1)+lung(nodoelim)
        densita(nodoelim-1)=opt*min(lung(nodoelim-2),&
             &   lung(nodoelim-1))+&
             &   (1-opt)*(lung(nodoelim-2)+lung(nodoelim-1))
        DO j=nodoelim,quanti-2
           quali(j)=quali(j+1)
           lung(j)=lung(j+1)
           densita(j)=opt*min(lung(j),lung(j-1))+&
                &   (1-opt)*(lung(j)+lung(j-1))
           distanza(j)=distanza(j+1) !!
        ENDDO
     ENDIF
     quali(quanti-1)=quali(quanti)
     quanti=quanti-1
  ENDDO
  
END SUBROUTINE uniform

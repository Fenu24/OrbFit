
! *************************************
!        S  O  L  V  P  O  L  Y        
! *************************************
! find the real roots of a polynomial with degree poldeg
! with real coefficients
 SUBROUTINE solvpoly(poldg,coef,roots,nroots,hzflag,multfl,radius)
   USE critical_points, ONLY: follow_roots
   USE output_control
   IMPLICIT NONE 
   INTEGER, INTENT(IN) :: poldg ! polynomial degree
   DOUBLE PRECISION, INTENT(IN) :: coef(0:poldg) ! z^i coeff. (i=0,poldg)   
   DOUBLE PRECISION, INTENT(OUT) :: roots(poldg) !real roots
   INTEGER, INTENT(OUT) :: nroots ! number of real roots
   LOGICAL,INTENT(INOUT) :: hzflag ! hzflag = .true.  OK!
                                   ! hzflag = .false. abs(root)>10^5
   LOGICAL,INTENT(INOUT) :: multfl ! multfl = .true.  OK!
                                   ! multfl = .false. 0 has multiplicity > 4
   DOUBLE PRECISION, INTENT(OUT)  :: radius(poldg) 
 ! =============== end interface =======================================
   REAL(KIND=8), PARAMETER :: epsil=10.d-30
   COMPLEX*16 :: zr1(poldg) 
   COMPLEX*16 :: poly(0:poldg)
   DOUBLE PRECISION :: epsm,big,small,rad(poldg),a1(poldg+1), &
        & a2(poldg+1),rim,rre 
   INTEGER :: nit,j,i,pdeg 
   DOUBLE PRECISION ,DIMENSION(poldg) :: zeros 
   INTEGER :: numzerosol,numzerocoe
   LOGICAL err(poldg+1) 
! =====================================================================
! initialization
   numzerosol = 0
   numzerocoe = 0
   roots(1:poldg) = 0.d0
   pdeg = poldg
   IF(verb_moid.ge.20)THEN
      WRITE(*,*)'poly coefficients:',coef(0:poldg)      
   ENDIF
! for polzeros
   epsm=2.D0**(-53) 
   big=2.D0**1023 
   small=2.D0**(-1022) 
   poly(0:pdeg)=coef(0:pdeg) 
   
! *************************
! look for zero solutions   
! *************************
   IF(abs(poly(0)).le.epsil) THEN
!      IF(verb_moid.ge.20) THEN
         WRITE(ierrou,*)'CONSTANT TERM is CLOSE to ZERO!',poly(0)
         numerr=numerr+1
!      ENDIF
      numzerosol = 1
      pdeg = pdeg-1
      IF(abs(poly(1)).le.epsil) THEN
!         IF(verb_moid.ge.20) THEN
            WRITE(ierrou,*)'FIRST DEGREE TERM is CLOSE to ZERO!',poly(1)
!         ENDIF
         numzerosol = 2
         pdeg = pdeg-1
         IF(abs(poly(2)).le.epsil) THEN
!            IF(verb_moid.ge.20) THEN
               WRITE(ierrou,*)'SECOND DEGREE TERM is CLOSE to ZERO!',poly(2)
!            ENDIF
            numzerosol = 3
            pdeg = pdeg-1
            IF(abs(poly(3)).le.epsil) THEN
!               IF(verb_moid.ge.20) THEN
                  WRITE(ierrou,*)'THIRD DEGREE TERM is CLOSE to ZERO!',poly(3)
!               ENDIF
               numzerosol = 4
               pdeg = pdeg-1
               IF(abs(poly(4)).le.epsil) THEN
!                  IF(verb_moid.ge.20) THEN
                     WRITE(ierrou,*)'solvpoly: ERROR! &
                          & zero solution has multiplicity > 4'
!                  ENDIF
                  multfl = .false.
               ENDIF
            ENDIF
         ENDIF
      ENDIF
   ENDIF
   
! ************************************************
! check the decrease of the degree up to poldg-4  
! ************************************************
   IF(abs(poly(poldg)).le.epsil) THEN
!      IF(verb_moid.ge.20) THEN
         WRITE(ierrou,*)'16th DEGREE TERM is CLOSE to ZERO!',poly(poldg)
         numerr=numerr+1
!      ENDIF
      pdeg = pdeg-1
      numzerocoe = 1
      IF(abs(poly(poldg-1)).le.epsil) THEN
!         IF(verb_moid.ge.20) THEN
            WRITE(ierrou,*)'15th DEGREE TERM is CLOSE to ZERO!',poly(poldg-1)
!         ENDIF
         pdeg = pdeg-1
         numzerocoe = 2
         IF(abs(poly(poldg-2)).le.epsil) THEN
!            IF(verb_moid.ge.20) THEN
               WRITE(ierrou,*)'14th DEGREE TERM is CLOSE to ZERO!', &
                    & poly(poldg-2)
!            ENDIF
            pdeg = pdeg-1
            numzerocoe = 3
            IF(abs(poly(poldg-3)).le.epsil) THEN
!               IF(verb_moid.ge.20) THEN
                  WRITE(ierrou,*)'13th DEGREE TERM is CLOSE to ZERO!', &
                       & poly(poldg-3)
!               ENDIF
               pdeg = pdeg-1
               numzerocoe = 4
            ENDIF
         ENDIF
      ENDIF
   ENDIF
   
! *********************************************************************   
   CALL polzeros(pdeg,poly(numzerosol:poldg-numzerocoe),epsm,big,small, &
        & 50,zr1,rad,err,nit,a1,a2) 
! *********************************************************************   
   
   nroots=0 
   
   zeros(1:poldg) = 0.d0
   IF(follow_roots) THEN
      WRITE(12,101) pdeg,zeros(2:poldg)
      WRITE(12,102) DBLE(zr1(1:pdeg)),zeros(pdeg+1:poldg)
      WRITE(12,102) DIMAG(zr1(1:pdeg)),zeros(pdeg+1:poldg)
   ENDIF
101 FORMAT(i20,15(2x,f20.8))
102 FORMAT(f20.8,15(2x,f20.8))

   DO 11 j = 1,pdeg

            if(verb_moid.ge.20) then
               WRITE(*,*)'ROOT(',j,')=',zr1(j)
            endif
!      rim=IMAG(zr1(j))                                                 
      rim=DIMAG(zr1(j)) 
!      IF(ABS(rim).LT.10.d-4) THEN                                      
!      IF(ABS(rim).LT.1.d-5) THEN                                      
      IF(ABS(rim).LT.rad(j)) THEN 
!      IF(ABS(rim).LT.2.d0*rad(j)) THEN 
!      IF(ABS(rim).LT.4.d0*rad(j)) THEN 
         rre=DBLE(zr1(j)) 
! uncomment if you want only positive real roots                    
!      IF(rre.LT.0.D0) GOTO 11                                      
         nroots=nroots+1 
         roots(nroots)=rre 
         radius(nroots)=rad(j) ! error estimate            

! ***************************
! check if some root is big
! ***************************
         IF(abs(roots(nroots)).gt.1.d8) THEN
!            if(verb_moid.ge.20) then
               write(ierrou,*)'large value for a root',roots(nroots)
               numerr=numerr+1
!            endif
            hzflag=.false.
!            write(*,*)'hzflag',hzflag
         ENDIF

      END IF
11 END DO
   
!   write(*,*)'roots:',roots(1:nroots)

   IF(numzerosol.gt.0) THEN
      DO i = 1,numzerosol
         nroots=nroots+1 
         roots(nroots)=0.d0 
         radius(nroots)=0.d0 ! dummy
!         write(*,*)'roots(',nroots,')=',roots(nroots)
      ENDDO
   ENDIF
   
 END SUBROUTINE solvpoly

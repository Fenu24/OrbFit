! =====================================================                 
!  INPTRO data input                                                    
SUBROUTINE inptwoway(iun1,iun2,t,el,nt,nbeg,naout,nax,    &
     &    number,elt,tt,jcon,ngfl,nt2) 
  USE fund_const
  USE name_rules, ONLY: name_len   
  USE synthcomp   
  IMPLICIT NONE 
  INTEGER, INTENT(IN) ::  iun1,iun2 
! max no asteroids in file
! pointer to the asteroid to be read first
  INTEGER, INTENT(IN) :: nax,nbeg
! number of asteroids output at once
  INTEGER, INTENT(OUT) :: naout 
! asteroid elements and time
  DOUBLE PRECISION, INTENT(OUT) :: el(7,nbbx,ntxx),t(ntxx)
  CHARACTER*(name_len), INTENT(OUT) ::  number(nax)
  INTEGER, INTENT(OUT) :: jcon(ntxx,nax)
! end interface
  DOUBLE PRECISION x(7),elt(7,ntxx),tt(ntxx),t0,t02  
  INTEGER ng(2)
  CHARACTER*(name_len) number2(nax) 
  CHARACTER coo*3,sys*3,ref*6,unit*3,comme*60,colhea*100 
  CHARACTER coo2*3,sys2*3,ref2*6,unit2*3,comme2*60,colhea2*100 
  LOGICAL consistent 
  INTEGER nt,na,ilce,na2,ilce2,nt2,ngfl 
  INTEGER i,j,n,no 
  CHARACTER*110 rec 
! ===========================================================           
!  read header of file on unit iun1                                     
  CALL reahea(iun1,number,comme,colhea,coo,sys,ref,unit,na,ilce,t0) 
! write the ones to be passed in this run                               
  IF(nbeg.gt.na)THEN 
     naout=0 
     RETURN 
  ELSE 
     naout=min(na-nbeg+1,nbb) 
  ENDIF
!     WRITE(*,100)(number(i),i=nbeg,nbeg+nbb-1)                         
  WRITE(*,100)number(nbeg),number(nbeg+naout-1) 
100 FORMAT(' 1st file, asteroids from :',1x,a9,'  to :',1x,a9) 
! ===========================================================           
!  read loop                                                            
  DO 1 j=1,ntx 
     READ(iun1,200,end=3,err=6)t(j) 
200  FORMAT(f12.2) 
     DO 10 n=1,na 
        jcon(j,n)=1 
        READ(iun1,101)rec 
101     FORMAT(a) 
! read in the record, with caution                                      
        READ(rec,201,end=4,err=4)(x(i),i=1,6),ng(1) 
201     FORMAT(f12.9,4f11.7,f11.8,1x,i9) 
        READ(rec,301,end=4,err=8)x(7) 
301     FORMAT(78x,e12.4) 
        GOTO 5 
!  second attempt                                                       
! 7        read(rec,202,end=4,err=4)(x(i),i=1,6),ng(1),x(7)             
! 202      FORMAT(f16.9,4f15.7,f11.8,1x,i9,1p,e12.4)                    
!          WRITE(*,*)'second attempt at record ',j,' ast ',number(n)    
!          GOTO 5                                                       
!  error cases                                                          
4       CONTINUE 
        IF(n.ge.nbeg.and.n.lt.nbeg+nbb)THEN 
!            WRITE(*,*)' error in elements, record ',j,' ast ',number(n)
           DO i=1,7 
              x(i)=0.d0 
           ENDDO
           ng(1)=0 
           jcon(j,n)=0 
           WRITE(22,*)' error in elements, record ',j,' ast ',number(n) 
           WRITE(22,101)rec 
        ENDIF
        GOTO 5 
8       IF(n.ge.nbeg.and.n.lt.nbeg+nbb)THEN 
!             WRITE(*,*)' error in gamma, record ',j,' ast ',number(n)  
           x(7)=0.d0 
           WRITE(22,*)' error in gamma, record ',j,' ast ',number(n) 
           WRITE(22,101)rec 
        ENDIF
! no error                                                              
5       CONTINUE 
        IF(n.ge.nbeg.and.n.lt.nbeg+nbb)THEN 
           no=n-nbeg+1 
           DO i=1,5 
              el(i,no,j)=x(i) 
           ENDDO
           el(6,no,j)=x(6)+dpig*ng(1) 
           el(7,no,j)=x(7) 
        ENDIF
10   ENDDO
1 ENDDO
! ending before end of file                                             
  WRITE(*,*)' too many records, max was',ntx 
! regular ending                                                        
3 nt=j-1 
  WRITE(*,*)' end of first file, number of records',nt 
! ===========================================================           
!  READ header of file on unit iun2                                     
  call reahea(iun2,number2,comme2,colhea2,coo2,sys2,ref2,unit2,na2, &
     & ilce2,t02)                                                       
! check consistency of the two files                                    
  IF(coo2.ne.coo.or.sys2.ne.sys.or.ref2.ne.ref.or.unit.ne.unit2     &
     &     .or.t0.ne.t02)THEN                                           
     WRITE(*,*)' inconsistency between files, coordinate system' 
     WRITE(*,*)coo,sys,ref,unit,t0,coo2,sys2,ref2,unit2,t02 
     STOP 
  ENDIF
  IF(na.ne.na2)THEN 
     WRITE(*,*)' inconsistency in content of files, no ast.',na,na2 
     STOP 
  ELSE 
     consistent=.true. 
     DO j=1,na2 
        IF(number(j).ne.number2(j))THEN 
           consistent=.false. 
           WRITE(*,*)' different asteroid: ',number(j),number2(j) 
        ENDIF
     ENDDO
     IF(.not.consistent)THEN 
        WRITE(*,*)' inconsistent list of asteroids ' 
        STOP 
     ENDIF
  ENDIF
  WRITE(*,111)number(nbeg),number(nbeg+naout-1) 
111 FORMAT('2nd file,  asteroids from :',1x,a9,'  to :',1x,a9) 
! ===========================================================           
!  read loop                                                            
  DO 11 j=1,ntx 
     READ(iun2,200,end=13,err=6)t(nt+j) 
     DO 110 n=1,na2 
        jcon(nt+j,n)=1 
        READ(iun2,101)rec 
! read in the record, with caution                                      
        READ(rec,201,end=14,err=14)(x(i),i=1,6),ng(1) 
        READ(rec,301,end=14,err=18)x(7) 
        GOTO 15 
!  second attempt                                                       
! 17       READ(rec,202,end=4,err=4)(x(i),i=1,6),ng(1),x(7)             
!          WRITE(*,*)'second attempt at record ',j,' ast ',number(n)    
!          GOTO 15                                                      
!  error cases                                                          
14      CONTINUE 
        IF(n.ge.nbeg.and.n.lt.nbeg+nbb)THEN 
           WRITE(*,*)' error in elements, record ',j,' ast ',number(n) 
           DO i=1,7 
              x(i)=0.d0 
           ENDDO
           ng(1)=0 
           jcon(nt+j,n)=0 
           WRITE(22,*)' error in elements, record ',j,' ast ',number(n) 
           WRITE(22,101)rec 
        ENDIF
        GOTO 15 
18      IF(n.ge.nbeg.and.n.lt.nbeg+nbb)THEN 
           WRITE(*,*)' error in gamma, record ',j,' ast ',number(n) 
           x(7)=0.d0 
           WRITE(22,*)' error in gamma, record ',j,' ast ',number(n) 
           WRITE(22,101)rec 
        ENDIF
! no error                                                              
15      CONTINUE 
        IF(n.ge.nbeg.and.n.lt.nbeg+nbb)THEN 
           no=n-nbeg+1 
           DO i=1,5 
              el(i,no,nt+j)=x(i) 
           ENDDO
           el(6,no,nt+j)=x(6)+dpig*ng(1) 
           el(7,no,nt+j)=x(7) 
        ENDIF
110  ENDDO
11 ENDDO
! ending before end of file                                             
  WRITE(*,*)' too many records, max was',ntx 
! regular ending                                                        
13 nt2=j-1 
  WRITE(*,*)' end of second file, number of records',nt2 
! =============================================================         
! now we assume that the first file is backward, the second is forward (
  DO 2 n=1,naout 
     DO 20 j=1,nt 
        DO i=1,7 
           elt(i,j)=el(i,n,j) 
        ENDDO
20   ENDDO
     DO 21 j=1,nt 
        DO i=1,7 
           el(i,n,nt+1-j)=elt(i,j) 
        ENDDO
21   ENDDO
2 ENDDO
  DO j=1,nt 
     tt(j)=t(j) 
  ENDDO
  DO j=1,nt 
     t(nt+1-j)=tt(j) 
  ENDDO
! =============================================                         
! reset the names                                                       
  DO no=1,naout 
     number(no)=number(nbeg+no-1) 
  ENDDO
! which file to compute LCE                                             
!      if(nt.gt.nt2)then                                                
  ngfl=1 
!      else                                                             
!         ngfl=2                                                        
!      endif                                                            
! total number of times                                                 
  nt=nt+nt2 
  RETURN 
! error return                                                          
6 WRITE(*,*)' error in time input, record ',j 
  STOP
END SUBROUTINE inptwoway

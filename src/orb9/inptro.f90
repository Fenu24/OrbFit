! =====================================================                 
!  INPTRO data input                                                    
SUBROUTINE inptro(iun,t,el,nt,nbeg,naout,nax,number,jcon)
  USE fund_const
  USE name_rules, ONLY: name_len
  USE synthcomp
  IMPLICIT NONE 
! input unit, opened
  INTEGER, INTENT(IN) :: iun
! max no ast. in file,
! pointer to the asteroid to be read first
  INTEGER, INTENT(IN) ::  nax,nbeg
! number of asteroids in output, number of times
  INTEGER, INTENT(OUT) :: naout,nt
! asteroid elements and time
  DOUBLE PRECISION, INTENT(OUT) ::  el(7,nbbx,ntxx),t(ntxx)
  CHARACTER*(name_len), INTENT(OUT) ::  number(nax) 
  INTEGER, INTENT(OUT) :: jcon(ntxx,nax)
! end interface
  DOUBLE PRECISION x(7), t0 
  INTEGER ng(2) 
  CHARACTER coo*3,sys*3,ref*6,unit*3,comme*60,colhea*100 
  INTEGER na,ilce 
  INTEGER i,j,n,no 
  CHARACTER*110 rec 
  INTEGER nerr 
! ===========================================================           
!  read header of file on unit iun                                      
  CALL reahea(iun,number,comme,colhea,coo,sys,ref,unit,na,ilce,t0) 
! WRITE the ones to be passed in this run                               
  IF(nbeg.gt.na)THEN 
     naout=0 
     RETURN 
  ELSE 
     naout=min(na-nbeg+1,nbb) 
  ENDIF
  WRITE(*,100)number(nbeg),number(nbeg+naout-1) 
100 FORMAT(' asteroids from :',1x,a9,'  to :',1x,a9) 
! ===========================================================           
  nerr=0 
!  read loop                                                            
  DO 1 j=1,ntx 
     READ(iun,200,end=3,err=6)t(j) 
200  FORMAT(f12.2) 
     DO 10 n=1,na 
        jcon(j,n)=1 
        READ(iun,101)rec 
101     FORMAT(a) 
! read in the record, with caution                                      
!          READ(rec,201,end=4,err=7)(x(i),i=1,6),ng(1)                  
        READ(rec,201,end=4,err=4)(x(i),i=1,6),ng(1) 
201     FORMAT(f12.9,4f11.7,f11.8,1x,i9) 
        READ(rec,301,end=4,err=8)x(7) 
301     FORMAT(78x,e12.4) 
        GOTO 5 
!  second attempt                                                       
! 7        READ(rec,202,end=4,err=4)(x(i),i=1,6),ng(1)                  
! 202      FORMAT(f16.9,4f15.7,f11.8,1x,i9)                             
!          READ(rec,302,end=4,err=8)x(7)                                
! 302      FORMAT(98x,e12.4)                                            
!          IF(n.ge.nbeg.and.n.lt.nbeg+nbb)THEN                          
!             WRITE(*,*)'second attempt at record ',j,' ast ',number(n) 
!          ENDIF                                                        
!          GOTO 5                                                       
!  error cases                                                          
4       CONTINUE 
        IF(n.ge.nbeg.and.n.lt.nbeg+nbb)THEN 
           nerr=nerr+1 
           IF(nerr.lt.100)THEN 
              WRITE(*,*)' error in elements, record ',j,' ast ',number(n) 
           ENDIF
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
           nerr=nerr+1 
           IF(nerr.lt.100)THEN 
              WRITE(*,*)' error in gamma, record ',j,' ast ',number(n) 
           ENDIF
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
  WRITE(*,*)' end of file, number of records',nt 
! reset the names                                                       
  DO no=1,naout 
     number(no)=number(nbeg+no-1) 
  ENDDO
  return 
! error return                                                          
6 WRITE(*,*)' error in time input, record ',j 
  STOP
END SUBROUTINE inptro

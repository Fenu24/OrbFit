! f90 quality control: done
! Copyright A. Milani 1997-2003                                   
MODULE util_suit
PUBLIC menu
CONTAINS
! ================================================================      
! MENU                                                                  
! ================================================================      
SUBROUTINE menu(ifl,menunam,nopt,s,                               &
     &             s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12) 
  IMPLICIT NONE 
  INTEGER, INTENT(IN) :: nopt
  CHARACTER*(*), INTENT(IN) ::  s,s1,s2
  CHARACTER*(*), INTENT(IN), OPTIONAL :: s3,s4,s5,s6,s7,s8,s9,s10,s11,s12 
  CHARACTER*20, INTENT(IN) :: menunam 
  INTEGER, INTENT(OUT) :: ifl
  CHARACTER*120 helpfi,ddocd1 
  INTEGER ll,iunit,lench
  INCLUDE 'doclib.h90' 
!                    
  IF(nopt.lt.2.or.nopt.gt.12)then 
     WRITE(*,*) ' this menu can only handle betw. 2 and 12 options' 
     ifl=0 
     return 
  ENDIF
3 continue 
  ll=index(s,'=') 
  WRITE(*,*) s(1:ll-1) 
  ll=index(s1,'=') 
  WRITE(*,*)' 1 = ', s1(1:ll-1) 
  ll=index(s2,'=') 
  WRITE(*,*)' 2 = ', s2(1:ll-1) 
  IF(nopt.lt.3) goto 2 
  ll=index(s3,'=') 
  WRITE(*,*)' 3 = ', s3(1:ll-1) 
  IF(nopt.lt.4) goto 2 
  ll=index(s4,'=') 
  WRITE(*,*)' 4 = ', s4(1:ll-1) 
  IF(nopt.lt.5) goto 2 
  ll=index(s5,'=') 
  WRITE(*,*)' 5 = ', s5(1:ll-1) 
  IF(nopt.lt.6) goto 2 
  ll=index(s6,'=') 
  WRITE(*,*)' 6 = ', s6(1:ll-1) 
  IF(nopt.lt.7) goto 2 
  ll=index(s7,'=') 
  WRITE(*,*)' 7 = ', s7(1:ll-1) 
  IF(nopt.lt.8) goto 2 
  ll=index(s8,'=') 
  WRITE(*,*)' 8 = ', s8(1:ll-1) 
  IF(nopt.lt.9) goto 2 
  ll=index(s9,'=') 
  WRITE(*,*)' 9 = ', s9(1:ll-1) 
  IF(nopt.lt.10) goto 2 
  ll=index(s10,'=') 
  WRITE(*,*)'10 = ', s10(1:ll-1)
 IF(nopt.lt.11) goto 2 
  ll=index(s11,'=') 
  WRITE(*,*)'11 = ', s11(1:ll-1) 
  IF(nopt.lt.12) goto 2 
  ll=index(s12,'=') 
  WRITE(*,*)'12 = ', s12(1:ll-1)
!                                                                       
! room to increase                                                      
!                                                                       
2 WRITE(*,103) 
103 format(' 0 = exit; -1=help') 
  WRITE(*,*)' selection?  ' 
  read(*,*,err=3)ifl 
! wrong flag and exit                                                   
4 IF(ifl.lt.-1.or.ifl.gt.nopt)THEN 
     WRITE(*,*)ifl,' option not understood' 
     goto 3 
  ELSEIF(ifl.eq.-1)THEN 
     ddocd1=ddocd 
     ll=lench(ddocd1) 
     helpfi=ddocd1(1:ll)//'/'//menunam 
     CALL rmsp(helpfi,ll) 
     helpfi=helpfi(1:ll)//'.help' 
     CALL filopn(iunit,helpfi,'OLD') 
     CALL filcat(iunit) 
     CALL filclo(iunit,' ') 
     WRITE(*,*)' selection?  ' 
     read(*,*,err=3)ifl 
     GOTO 4 
  ENDIF
  RETURN 
END SUBROUTINE menu
END module util_suit
! ===========================================                           
SUBROUTINE filcat(iunit) 
  IMPLICIT NONE 
  INTEGER iunit,i,imax,ll,lcom,lench 
  PARAMETER (imax=100) 
  CHARACTER*100 record 
! =========================================                             
  DO i=1,imax 
     READ(iunit,100,END=2)record 
100  FORMAT(a) 
     ll=lench(record) 
!        WRITE(*,*)ll                                                   
! comments begin with %, for TeX compatibility                          
     lcom=index(record,'%') 
     IF(lcom.gt.1)THEN 
        WRITE(*,100)record(1:lcom-1) 
     ELSEIF(lcom.eq.0)THEN 
        IF(ll.gt.0)THEN 
           WRITE(*,100)record(1:ll) 
        ELSE 
           WRITE(*,*) 
        ENDIF
     ENDIF
  ENDDO
2 RETURN 
END SUBROUTINE filcat
! ===================================================================== 
! TEE (as in unix shell)                                                
! ===================================================================== 
subroutine tee(iun,string) 
  implicit none 
  integer iun 
  character*(*) string 
  integer le 
  le=index(string,'=') 
  write(*,*)string(1:le-1) 
  write(iun,*)string(1:le-1) 
  return 
END subroutine tee

SUBROUTINE write_err(name0,iunout,string)
  USE output_control
  IMPLICIT NONE
  CHARACTER*(*),INTENT(IN):: name0
  CHARACTER*(*),INTENT(IN):: string
  INTEGER,INTENT(IN):: iunout
  WRITE(*,*)name0,string
  WRITE(iunout,*)name0,string
  WRITE(ierrou,*)name0,string
  numerr=numerr+1
END SUBROUTINE write_err

SUBROUTINE final_rep(i,ix,run,le,iunout,name0)
  USE output_control
  USE name_rules
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: i,ix,iunout,le
  CHARACTER*(name_len), INTENT(IN) :: name0
  CHARACTER*(*),INTENT(IN):: run
  INTEGER iundone,le1
  CHARACTER*20 run1
  run1=run
  CALL rmsp(run1,le1)
  CALL filopn(iundone,run1(1:le1)//'.done','unknown')
  IF((ix.gt.0.and.i.lt.ix).or.ix.eq.0)THEN
     WRITE(*,*) ' regular end; processed ',i,'  asteroids'
     WRITE(iundone,*) ' regular end; processed ',i,'  asteroids'
     CALL filclo(iundone,' ')
  ELSE
     CALL write_err(name0,iunout,' last asteroid processed, increase ix')
     CALL filclo(iundone,'DELETE')
  ENDIF
  IF(numerr.gt.0)THEN 
     CALL filclo(ierrou,' ') 
  ELSE 
     CALL filclo(ierrou,'DELETE') 
  ENDIF
END SUBROUTINE final_rep


MODULE name_rules

IMPLICIT NONE
PRIVATE
! length of official IAU designation, without spaces                    
  INTEGER, PUBLIC, PARAMETER :: name_len=9 
  INTEGER, PUBLIC, PARAMETER :: nmax_ids=20
  INTEGER, PUBLIC, PARAMETER :: idname_len=(name_len+1)*nmax_ids-1 
  INTEGER, PUBLIC, PARAMETER :: idnamvir_len=idname_len+6

PUBLIC fullnamparse, appears

CONTAINS
!=====================================================================
! APPEARS name appears in list?
! ====================================================================
 LOGICAL FUNCTION appears(secn,sec_nam,nsec)
   CHARACTER*(name_len), INTENT(IN) :: secn
   INTEGER, INTENT(IN) :: nsec
   CHARACTER*(name_len), DIMENSION(nsec), INTENT(IN) :: sec_nam
   INTEGER i
   appears=.false.
   DO i=1,nsec
     IF(secn.eq.sec_nam(i)) appears=.true.
   ENDDO
 END FUNCTION appears
! =======================================                               
! fullnamparse                                                              
! parse asteroid name to find identifications/multiple solutions        
! to be used by both vsa_attrib and vsa_attmore                         
! fullname is of the general form                                       
! name[=names(1)][=names(2)]=....=[_nvir]    [optional]
! if missing nvir is store as 0       
! nam0 has the _nvir removed                          
! =======================================                               
SUBROUTINE fullnamparse(fullname,srt,nam0,name,nid,nvir,ierror) 
! ------------input------------------                                   
  CHARACTER*(*), INTENT(IN) :: fullname
  LOGICAL, INTENT(IN) :: srt ! to sort or not to sort? 
! ------------output-----------------  
  CHARACTER*(name_len), INTENT(OUT) :: name(nmax_ids) 
  INTEGER, INTENT(OUT) :: nid 
  CHARACTER*(idname_len), INTENT(OUT) :: nam0 
  INTEGER nvir ! index of multiple solution, if applicable
  INTEGER ierror ! error flag 
! ------------end interface------------- 
  CHARACTER*(idname_len) longname
! length of fullname, of longname, position of _, position of =
  INTEGER le,le1,ll,ll1 
  INTEGER ieq,id
  CHARACTER*(name_len) nam
  INTEGER ind(nmax_ids), nvtr ! for sorting, which design. has triangles
  DOUBLE PRECISION times(nmax_ids) ! for sorting
  CHARACTER*(name_len) :: nameunsort(nmax_ids) 
! ==========================================                            
  ierror=0
! parse fullname to find index of virtual asteroid  
  longname=fullname                    
  CALL rmsp(longname,le) 
  IF(le.gt.idnamvir_len)THEN
     WRITE(*,*)'fullnamparse: too many ids. ',fullname(1:le),le 
     ierror=5
  ENDIF
  ll=index(longname,'_') 
  IF(ll.gt.0)THEN 
     IF(le-ll.eq.1)THEN 
        READ(longname(ll+1:ll+1),111)nvir 
111     format(i1) 
     ELSEIF(le-ll.eq.2)THEN 
        READ(longname(ll+1:ll+2),112)nvir 
112     format(i2) 
     ELSEIF(le-ll.eq.3)THEN 
        READ(longname(ll+1:ll+3),113)nvir 
113     format(i3) 
     ELSEIF(le-ll.eq.4)THEN 
        READ(longname(ll+1:ll+4),114)nvir 
114     format(i4) 
     ELSE 
        WRITE(*,*)'namparse: error in reading nvir ',longname,le,ll 
        ierror=1 
        RETURN 
     ENDIF
     longname(ll:)=' '
     le1=ll-1
  ELSEIF(ll.eq.0)THEN 
! single solution, give default number                                  
     nvir=0 
     le1=le
  ELSE 
     WRITE(*,*)'fullnamparse: error in _ sign ',longname(1:le),le,ll 
     ierror=2 
     STOP 
  ENDIF
! =============================================                         
! parse fullname to find identification marks
  DO nid=1,nmax_ids  
     ll1=index(longname,'=') 
     IF(ll1.eq.0)THEN 
! no more identifications
        IF(le1.gt.name_len)THEN
           WRITE(*,*)'fullnamparse: too long name truncated',longname(1:ll1-1) 
           ierror=3 
           nameunsort(nid)=longname(1:name_len)
           longname=' '
           EXIT
        ELSE
           nameunsort(nid)=longname(1:le1)
           longname=' '
           EXIT
        ENDIF        
     ELSEIF(ll1.gt.0)THEN
! check length of first name
        IF(ll1.gt.name_len+1)THEN
           WRITE(*,*)'fullnamparse: too long name truncated',longname(1:ll1-1) 
           ierror=3 
           nameunsort(nid)=longname(1:name_len)
           longname=longname(ll1+1:le1)
           le1=le1-ll1
        ELSEIF(ll1.gt.1)THEN
           nameunsort(nid)=longname(1:ll1-1)
           longname=longname(ll1+1:le1)
           le1=le1-ll1
        ELSE
           WRITE(*,*)'fullnamparse: error == ',fullname(1:le),le,ll,ll1
           ierror=4
           STOP
        ENDIF
     ELSE 
        WRITE(*,*)'fullnamparse: error in = sign ',fullname(1:le),le,ll,ll1 
        ierror=6 
        STOP 
     ENDIF
!     CALL rmsp(longname,le1)
!     IF(le1.eq.0)EXIT
  ENDDO
! nid=total number of designations contained in the string 
  IF(srt)THEN
!     WRITE(*,*) ' fullnamparse: sort function not ready'
!     STOP
! find day
     DO id=1,nid
        nam=nameunsort(id)
        CALL decode_date_mjd(nam(5:7),times(id))
     ENDDO
! sort by date
     CALL heapsortname(nameunsort,times,nid,ind)
! identification name, with ids sorted.
     nam0=nameunsort(ind(1))
     name(1)=nameunsort(ind(1))
     DO id=2,nid
        ieq=(name_len+1)*(id-1)
        nam0(ieq:ieq)='='
        nam0(ieq+1:ieq+name_len)=nameunsort(ind(id))
        name(id)=nameunsort(ind(id))
     ENDDO
! nvir encodes which triangulation was used in the first digit
     nvtr=nvir/1000
     nvir=nvir-nvtr*1000
     IF(nvtr.eq.0)THEN
        nvir=ind(1)*1000+nvir
     ELSE
        nvir=ind(nvtr)*1000+nvir
     ENDIF
  ELSE ! not to sort
! identification name, with ids not sorted.
     nam0=nameunsort(1)
     name(1)=nameunsort(1)
     DO id=2,nid
        ieq=(name_len+1)*(id-1)
        nam0(ieq:ieq)='='
        nam0(ieq+1:ieq+name_len)=nameunsort(id)
        name(id)=nameunsort(id)                
     ENDDO
  ENDIF
  DO id=nid+1,nmax_ids 
     name(id)=' '
  ENDDO
END SUBROUTINE fullnamparse

! ===================================                                   
! HEAPSORTNAME                                                             
! sorting by heapsort of an array a of names of lenght n         
! the input array a is not changed, in output                           
! ind is the indirect address of the sorted array                       
! ===================================                                   
SUBROUTINE heapsortname(a,times,n,ind) 
  INTEGER,INTENT(IN) ::  n 
  CHARACTER(name_len),INTENT(IN) :: a(n)
  DOUBLE PRECISION, INTENT(IN) :: times(n) 
  INTEGER, INTENT(OUT) :: ind(n) 
! =====end interface========                                            
  integer j,l,ir,indt,i 
  CHARACTER(name_len) q
  DOUBLE PRECISION tq
! initialise indexes to natural order                                   
  DO j=1,n 
     ind(j)=j 
  ENDDO
  IF(n.eq.1)RETURN 
! counters for recursion, length of heap                                
  l=n/2+1 
  ir=n 
! recursive loop                                                        
1 CONTINUE 
  IF(l.gt.1)THEN 
     l=l-1 
     indt=ind(l) 
     q=a(indt) 
     tq=times(indt)
  ELSE 
     indt=ind(ir) 
     q=a(indt)
     tq=times(indt)
     ind(ir)=ind(1) 
     ir=ir-1 
     IF(ir.eq.1)THEN 
        ind(1)=indt 
        RETURN 
     ENDIF
  ENDIF
  i=l 
  j=l+l 
2 IF(j.le.ir)THEN 
     IF(j.lt.ir)THEN 
        IF(namless(a(ind(j)),times(ind(j)),a(ind(j+1)),times(ind(j+1))))j=j+1 
     ENDIF
     IF(namless(q,tq,a(ind(j)),times(ind(j))))THEN 
        ind(i)=ind(j) 
        i=j 
        j=j+j 
     ELSE 
        j=ir+1 
     ENDIF
     GOTO 2 
  ENDIF
  ind(i)=indt 
  GOTO 1 
END SUBROUTINE heapsortname

LOGICAL FUNCTION namless(nama,tima,namb,timb)
  CHARACTER(name_len), INTENT(IN) :: nama, namb
  DOUBLE PRECISION, INTENT(IN) :: tima, timb 
  INTEGER n_a, n_b
  LOGICAL d_a, d_b 
  n_a=numberrr(nama) ! +infty if not numbered
  n_b=numberrr(namb)
  IF(n_a.lt.n_b)THEN
     namless=.true. ! a numbered always comes before 
     RETURN ! if both numbered the lower numbered is before
  ENDIF
  d_a=designatedd(nama)
  d_b=designatedd(namb)
  IF(d_a.and..not.d_b)THEN
     namless=.true. ! a designation comes before a ONS
     RETURN 
  ELSEIF(d_b.and..not.d_a)THEN
     namless=.false. ! a designation comes before a ONS
     RETURN
  ELSEIF(d_a.and.d_b)THEN
! temporary, prides to be used
     IF(tima.lt.timb)THEN
        namless=.true.
     ELSE
        namless=.false.
     ENDIF
  ELSE
     IF(tima.lt.timb)THEN
        namless=.true. ! between ONS, only time counts
     ELSE
        namless=.false.
     ENDIF
  ENDIF
END FUNCTION namless

INTEGER FUNCTION numberrr(nama)  
  CHARACTER(name_len), INTENT(IN) :: nama
  READ(nama,*, err=9) numberrr
  RETURN
9 numberrr=100000000
END FUNCTION numberrr

LOGICAL FUNCTION designatedd(nama)
  CHARACTER(name_len), INTENT(IN) :: nama
  INTEGER num
  READ(nama(1:4),*, err=9)num
  designatedd=.true. 
  RETURN
9 designatedd=.false.
END FUNCTION designatedd

END MODULE name_rules


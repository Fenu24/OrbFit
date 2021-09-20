SUBROUTINE bin_search_rea(a,ntot,i1,i2,alist,ordind,indnum)
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT(IN) :: a ! number to search for in the list
  INTEGER,INTENT(IN) :: ntot,i1,i2 ! number of elements in the list
! pointer to the first and last to be considered
  DOUBLE PRECISION,INTENT(IN) :: alist(ntot) ! unordered list of names
  INTEGER,INTENT(IN) :: ordind(ntot) ! ordered addresses of the names
  INTEGER, INTENT(OUT) :: indnum ! address of the max of weak minorants in the list
                                 ! sorted as indicated by ordind 
! end interface
  INTEGER nlist,fst,nmez
  INTEGER i ! loop indexes
  nlist=i2-i1+1
  fst=i1 ! firt extremum of the search interval
  DO i = 1,ntot
     IF(nlist.eq.1) THEN
        IF(fst.eq.i2)THEN
           IF(a.ge.alist(ordind(fst)))THEN
              indnum=fst
              EXIT
           ELSE
              indnum=fst-1 ! number is just before 
           ENDIF
        ELSEIF(a.ge.alist(ordind(fst)).and.a.lt.alist(ordind(fst+1)))THEN
           indnum=fst
        ELSE
           indnum=fst-1 ! number is just before 
        ENDIF   
        EXIT
     ENDIF
     nmez=nlist/2
     IF(a.lt.alist(ordind(fst+nmez))) THEN
        fst=fst
        nlist=nmez        
     ELSE
        fst=fst+nmez
        nlist=nlist-nmez        
     ENDIF
  ENDDO
  
END SUBROUTINE bin_search_rea

SUBROUTINE bin_search(name,ntot,namlist,ordind,indname) !,indordname)
  USE name_rules
  IMPLICIT NONE
  CHARACTER*(*),INTENT(IN) :: name ! name to search for in the list
  INTEGER,INTENT(IN) :: ntot ! number of elements in the list
  CHARACTER*(*),INTENT(IN) :: namlist(ntot) ! unordered list of names
  INTEGER,INTENT(IN) :: ordind(ntot) ! addresses of the names in the original
! list, sorted
  INTEGER, INTENT(OUT) :: indname ! address of the name we search for
                                  ! in the unordered list
!  INTEGER, INTENT(OUT) :: indordname ! address of the name we search for
                                  ! in the ordered list
! end interface
  INTEGER :: nlist,fst,nmez
  INTEGER i ! loop indexes

  nlist=ntot
  fst=1 ! firt extremum of the search interval
  DO i = 1,ntot
     IF(nlist.eq.1) THEN
        IF(name.eq.namlist(ordind(fst)))THEN
           indname=ordind(fst)
!           indordname=fst
        ELSE
!           WRITE(*,*) ' bin_search: ', name,' .ne. ', namlist(ordind(fst))
           indname=0 ! error message: name has not been found
!           indordname=0
        ENDIF   
        EXIT
     ENDIF
     nmez=nlist/2
     IF(LLT(name,namlist(ordind(fst+nmez)))) THEN
        fst=fst
        nlist=nmez        
     ELSE
        fst=fst+nmez
        nlist=nlist-nmez        
     ENDIF
  ENDDO

END SUBROUTINE bin_search

! ===================================                                   
! HEAPSORTNACH                                                            
! sorting by heapsort of a name array of lenght n 
! by using lexicographic order        
! the input array a is not changed, in output                           
! ind is the indirect address of the sorted array                       
! ===================================                                   
SUBROUTINE heapsortnach(a,n,ind) 
  USE name_rules
  IMPLICIT NONE
  INTEGER, INTENT(IN) ::  n 
  CHARACTER(name_len),INTENT(IN) :: a(n)
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
  ELSE 
     indt=ind(ir) 
     q=a(indt)
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
        IF(LLT(a(ind(j)),a(ind(j+1))))j=j+1 
     ENDIF
     IF(LLT(q,a(ind(j))))THEN 
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
END SUBROUTINE heapsortnach

! ===================================                                   
! HEAPSORTNACHLO                                                            
! sorting by heapsort of a CHARACTER*(le) array of lenght n 
! by using lexicographic order        
! the input array a is not changed, in output                           
! ind is the indirect address of the sorted array                       
! ===================================                                   
SUBROUTINE heapsortnachlo(a,n,le,q,ind) 
  IMPLICIT NONE
  INTEGER,INTENT(IN) ::  n, le 
  CHARACTER(le),INTENT(IN) :: a(n)
  CHARACTER(le), INTENT(OUT) :: q
  INTEGER, INTENT(OUT) :: ind(n) 
! =====end interface========                                            
  integer j,l,ir,indt,i 
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
  ELSE 
     indt=ind(ir) 
     q=a(indt)
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
        IF(LLT(a(ind(j)),a(ind(j+1))))j=j+1 
     ENDIF
     IF(LLT(q,a(ind(j))))THEN 
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
END SUBROUTINE heapsortnachlo

SUBROUTINE bin_search_int(iname,ntot,inamlist,ordind,indname)
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: iname ! number corresp. to name to search for in the list
  INTEGER,INTENT(IN) :: ntot ! number of elements in the list
  INTEGER,INTENT(IN) :: inamlist(ntot) ! unordered list of numbers(names)
  INTEGER,INTENT(IN) :: ordind(ntot) ! ordered addresses of the names
  INTEGER, INTENT(OUT) :: indname ! address of the name we search for
! end interface
  INTEGER :: nlist,fst,nmez
  INTEGER i ! loop indexes

  nlist=ntot
  fst=1 ! first extremum of the search interval
  DO i = 1,ntot
     IF(nlist.eq.1) THEN
        IF(iname.eq.inamlist(ordind(fst)))THEN
           indname=ordind(fst)
        ELSE
           WRITE(*,*) ' bin_search_int: ', iname,' .ne. ', inamlist(ordind(fst))
           indname=0 ! error message: name has not been found
        ENDIF   
        EXIT
     ENDIF
     nmez=nlist/2
     IF(iname.lt.inamlist(ordind(fst+nmez))) THEN
        fst=fst
        nlist=nmez        
     ELSE
        fst=fst+nmez
        nlist=nlist-nmez        
     ENDIF
   ENDDO
  
END SUBROUTINE bin_search_int

MODULE semi_linear
  USE fund_const, ONLY: qkind
  IMPLICIT NONE 

  PUBLIC linobs, ellips, ellipsoid, slinel
  PUBLIC slinel_q ! quadruple precision

! ========module semi_linear================== 
! routines shared by target_plane and by predict_obs                        
! CONTAINS                                                              
! PUBLIC  
!             linobs                                                
!             ellips  
!             ellipsoid                                              
!             elemov                                                
!             slinel                                                
!             graha 

CONTAINS

                          
! ===========================================================           
! LINOBS defines line of changes in orbital elements to be used for     
! confidence boundary on sky plane                                  
! ===========================================================           
  SUBROUTINE linobs(ibv,npo,el,axes,sig,b,v,sigma,ceicel,elm,dynm,npo1,nd) 
    USE fund_const
    USE orbit_elements
    USE dyn_param
! ====================INPUT==================================           
    INTEGER, INTENT(IN) :: ibv ! ibv=1 contour ibv=2 LOV
    INTEGER, INTENT(IN) :: npo ! requested points
    TYPE(orbit_elem), INTENT(IN) :: el  
    DOUBLE PRECISION, INTENT(IN) :: axes(2,2),sig(2),b(2,2),sigma 
! matrix defining the plane of the ellipse,new orthonormal reference    
    DOUBLE PRECISION, INTENT(IN) :: ceicel(nd-2,2),v(nd,nd) 
    INTEGER, INTENT(IN) :: nd ! dimension of covariance matrix
! ===================OUTPUT==================================           
    INTEGER, INTENT(OUT) :: npo1  ! obtained points
    DOUBLE PRECISION, INTENT(OUT) :: elm(6,npo),dynm(ndyx,npo) 
! ==================END INTERFACE============================           
    INTEGER nn,n,i,k
    DOUBLE PRECISION s,x,y,vad(2),xv,yv,dn,dth,theta,xa,yd 
    DOUBLE PRECISION eqnew(6), eq(6), del(ndimx)
    DOUBLE PRECISION alde(2),ecc 
! =====================================================================
    eq=el%coord
! line of maximum variation: in the alpha-delta plane                   
    vad(1:2)=axes(1:2,2)*sig(2) 
! in the elements space                                                 
    xv=(b(1,1)*vad(1)+b(1,2)*vad(2)) 
    yv=(b(2,1)*vad(1)+b(2,2)*vad(2)) 
! linear step for variation axis parametrisation                        
    dn=2.d0/float(npo-1) 
! angular step for ellipse parametrisation                              
    dth=dpig/float(npo) 
! ===========================================================           
! main loop on the number of output points                              
    nn=0 
    DO 7 n=1,npo 
! ===========================================================           
! choice between two output options                                     
       IF(ibv.eq.2)THEN 
! ===========================================================           
! line of maximum variation in the elements space                       
          s=(n-1)*dn-1.d0 
          x=sigma*s*xv 
          y=sigma*s*yv 
       ELSEIF(ibv.eq.1)THEN 
! ===================================================================== 
! parametrisation of the ellipse in the subspace of elements, based upon
! parametrisation of the ellipse in the alpha-delta plane               
! WARNING: npo must be divisible by 2, otherwise one tip of the         
! banana would be missed                                                
          theta=(n-1)*dth 
          xa=sig(1)*cos(theta)*sigma 
          yd=sig(2)*sin(theta)*sigma 
          alde=xa*axes(1:2,1)+yd*axes(1:2,2) 
! transfer of parametrisation in the V1,V2 plane                        
          x=(b(1,1)*alde(1)+b(1,2)*alde(2)) 
          y=(b(2,1)*alde(1)+b(2,2)*alde(2)) 
       ELSE 
          write(*,*)' linobs: this should not happen,ibv=',ibv 
       ENDIF
! compute displacement on the confidence ellipsoid corresponding to x,y 
       nn=nn+1
       CALL elemov(x,y,v,ceicel,del(1:nd),nd) 
       elm(1:6,nn)=del(1:6)
       IF(nd.gt.6)THEN
          dynm(1:ndyx,nn)=0.d0
          DO k=7,nd
             dynm(ls(k-6),nn)=del(k)          
          ENDDO
       ENDIF
! add to the original center of the ellipsoid of confidence             
       eqnew=eq+elm(1:6,nn) 
       IF(el%coo.eq.'EQU')THEN
          ecc=sqrt(eqnew(2)**2+eqnew(3)**2) 
          IF(ecc.ge.1.d0.or.eqnew(1).le.0.d0)THEN 
             write(*,*)'point ', nn,'  Hyperbolic, ecc=',ecc,' a=',eqnew(1) 
             nn=nn-1 
          ELSEIF(ecc.ge.0.99d0)THEN 
             write(*,*)'point ', nn,' Almost Hyperbolic, ecc=',ecc,' a=',eqnew(1)
             nn=nn-1 
          ENDIF
       ENDIF
7   ENDDO
! final count of non hyperbolic orbits                                  
    npo1=nn 
  END SUBROUTINE linobs
! ===================================================================== 
! ELEMOV                                                                
! compute displacement on the confidence ellipsoid corresponding to x,y 
! on the plane of the gradients of alpha-delta                          
! ===================================================================== 
  SUBROUTINE elemov(x,y,v,ceicel,del,nd) 
! inout/output                                                          
    DOUBLE PRECISION, INTENT(IN) :: x,y,v(nd,nd),ceicel(nd-2,2)
    INTEGER, INTENT(IN) :: nd
    DOUBLE PRECISION, INTENT(OUT) :: del(nd) 
! workspace                                                             
    DOUBLE PRECISION dee(nd-2),deel(nd) 
! ===================    
    del=x*v(1:nd,1)+y*v(1:nd,2) ! in the plane generated by rows of partials
    dee=-x*ceicel(1:nd-2,1)-y*ceicel(1:nd-2,2) ! 
    deel=MATMUL(v(1:nd,3:nd),dee) ! 
    del=del+deel ! 
  END SUBROUTINE elemov
! ===================================================================== 
! ELLIPS                                                                
! compute covariance ellipse of two observables                         
! ===================================================================== 
  SUBROUTINE ellips(daddet,gamm0,sig,axes,gamad,nd) 
! input covariance matrix                                               
    DOUBLE PRECISION, INTENT(IN) :: gamm0(nd,nd) 
! input partial derivatives of alpha, delta, w.r. to elements (by column
    DOUBLE PRECISION,INTENT(IN) :: daddet(nd,2) 
   INTEGER, INTENT(IN) :: nd ! dimension of covariance matrix
! output covariance                                                     
    DOUBLE PRECISION, INTENT(OUT) :: gamad(2,2),axes(2,2),sig(2) 
! ==============END INTERFACE========================================== 
! eigenvalues, workspace, transposed                                    
    DOUBLE PRECISION eigval(2),tmp26(2,nd),fv1(2),fv2(2),dadde(2,nd) 
! loop indexes                                                          
    INTEGER i 
! error flag                                                            
    INTEGER ierr 
! ===================================================================== 
    dadde=TRANSPOSE(daddet) 
    tmp26=MATMUL(dadde,gamm0)  
    gamad=MATMUL(tmp26,daddet) 
! ===================================================================== 
! compute ellipse of confidence                                         
! eigenvalues                                                           
    CALL rs(2,2,gamad,eigval,1,axes,fv1,fv2,ierr) 
    DO  i=1,2 
       IF(eigval(i).gt.0.d0)THEN 
          sig(i)=sqrt(eigval(i)) 
       ELSE 
          write(*,*) 'ellips: non positive eigenvalue' 
          sig(i)=0.d0 
       ENDIF
    ENDDO
  END SUBROUTINE ellips
! ===================================================================== 
! ELLIPSOID                                                                
! compute covariance ellipsoid of four observables                         
! ===================================================================== 
  SUBROUTINE ellipsoid(nd,daddet,gamm0,sig,axes,gamad) 
    USE output_control
    USE dyn_param, ONLY: ndimx
    INTEGER, INTENT(IN) :: nd ! dimension of covariance matrix
! input covariance matrix                                               
    DOUBLE PRECISION gamm0(nd,nd) 
! input partial derivatives of alpha, delta, w.r. to elements (by column
    DOUBLE PRECISION daddet(nd,4) 
! output covariance                                                     
    DOUBLE PRECISION gamad(4,4),axes(4,4),sig(4) 
! ==============END INTERFACE========================================== 
! eigenvalues, workspace, transposed                                    
    DOUBLE PRECISION eigval(4),tmp46(4,ndimx),fv1(4),fv2(4),dadde(4,ndimx) 
! loop indexes                                                          
    INTEGER i 
! error flag                                                            
    INTEGER ierr 
! ===================================================================== 
    dadde(1:4,1:nd)=TRANSPOSE(daddet) 
    tmp46(1:4,1:nd)=MATMUL(dadde(1:4,1:nd),gamm0)
    gamad=MATMUL(tmp46(1:4,1:nd),daddet)
! ===================================================================== 
! compute ellipsoid of confidence                                         
! eigenvalues                                                           
    CALL rs(4,4,gamad,eigval,1,axes,fv1,fv2,ierr) 
    DO  i=1,4
       IF(eigval(i).gt.0.d0)THEN 
          sig(i)=sqrt(eigval(i)) 
       ELSE 
!        write(ierrou,*) 'ellipsoid: non positive eigenvalue ',eigval(i)
!        numerr=numerr+1 
          sig(i)=0.d0 
       ENDIF
    ENDDO
  END SUBROUTINE ellipsoid
  !
  !=======================================================================!
  ! SLINEL                                                                !
  !=======================================================================!
  ! Semilinear boundary ellipse computation                               !
  !=======================================================================!
  SUBROUTINE slinel(dtpdet,gc,cc,ceicel,b,v,nd,bad_cond) 
! nd by 2 matrix with columns= gradients                                 
    DOUBLE PRECISION, INTENT(IN) :: dtpdet(nd,2) 
! normal and covariance matrices                                        
    DOUBLE PRECISION, INTENT(IN) :: gc(nd,nd),cc(nd,nd) 
! orthonormal basis                                                     
    DOUBLE PRECISION, INTENT(OUT):: v(nd,nd)
! inverted 2x2 matrix                                             
    DOUBLE PRECISION, INTENT(OUT):: b(2,2)
! regression formula
    DOUBLE PRECISION, INTENT(OUT):: ceicel(nd-2,2)
! number of variables to be determined
    INTEGER, INTENT(IN) :: nd
    INTEGER, INTENT(OUT), OPTIONAL :: bad_cond ! Flag for bad conditioning cases
! END INTERFACE
! orthonormal basis                                                     
    DOUBLE PRECISION vt(nd,nd),gamv(nd,nd),cv(nd,nd),tmp(nd,nd) 
! partial matrices                                                      
    DOUBLE PRECISION c4(nd-2,nd-2),cinv(nd-2,nd-2),c42(nd-2,2) 
! line of maximum variation                                             
    DOUBLE PRECISION a(2,2),deta 
! loop indexes ii=1,2, ij,ijj=1,nd-2                                       
    INTEGER ii, ij, ijj 
! for inversion with tcholevski: workspace, error flag                  
    DOUBLE PRECISION ws(nd-2) 
    INTEGER ierr 
    INTEGER :: cond_flag
! ===================================================================== 
    ! Initialization
    v=0.d0
    b=0.d0
    ceicel=0.d0
    IF(PRESENT(bad_cond)) bad_cond=0
    ! adapted orthonormal basis, covariance and normal matrix in the new bas
    CALL graha(dtpdet(1:nd,1:2),nd,v,cond_flag)
    ! Check against bad conditioning 
    IF(cond_flag.NE.0)THEN
       WRITE(*,*) ' slinel: bad conditioning case, return!'
       IF(PRESENT(bad_cond)) bad_cond=cond_flag
       RETURN
    END IF
    vt=TRANSPOSE(v) 
    tmp=MATMUL(vt,gc) 
    gamv=MATMUL(tmp,v) 
    tmp=MATMUL(vt,cc)  
    cv=MATMUL(tmp,v)   
! ===================================================================== 
! (nd-2)x(nd-2) and (nd-2)x2 submatrices of normal matrix 
    DO  ijj=1,nd-2
       DO ij=1,nd-2 
          c4(ijj,ij)=cv(ijj+2,ij+2) 
       ENDDO
       DO  ii=1,2 
          c42(ijj,ii)=cv(ijj+2,ii) 
       ENDDO
    ENDDO
! ===========================================================           
! Cholewski method for inversion                                        
    CALL tchinv(c4,nd-2,cinv,ws,ierr) 
    IF(ierr.NE.0)THEN 
       WRITE(*,*) ' slinel: degenerate case after tchinv, ierr=',ierr 
    ENDIF
! ===========================================================           
! matrix to be used for out of plane component                          
    ceicel=MATMUL(cinv,c42) ! CALL mulmat(cinv,4,4,c42,4,2,ceicel) 
! ===========================================================           
! linear map from the elements space (with base V) and the alpha-delta p
    a(1,1)=DOT_PRODUCT(dtpdet(1:nd,1),v(1:nd,1)) 
    a(1,2)=DOT_PRODUCT(dtpdet(1:nd,1),v(1:nd,2)) 
    a(2,1)=DOT_PRODUCT(dtpdet(1:nd,2),v(1:nd,1)) 
    a(2,2)=DOT_PRODUCT(dtpdet(1:nd,2),v(1:nd,2)) 
    CALL inv22(a,b,deta) 
  END SUBROUTINE slinel
  !
  !=============================================================================!
  ! GRAHA                                                                       !
  !=============================================================================!
  ! Graham-Schmidt procedure to generate an orthonormal basis v                 !
  ! starting from 2 n-vectors, given by means of a matrix a. The new            !
  ! basis must be such that the first 2 vectors are a basis for the             !
  ! space spanned by the 2 columns of a.                                        !
  !=============================================================================!
  SUBROUTINE graha(a,n,v,cond_flag)
    USE dyn_param, ONLY: ndimx 
    !=============================================================================================
    INTEGER,          INTENT(IN)            :: n         ! Number of vectors
    DOUBLE PRECISION, INTENT(IN)            :: a(n,2)    ! Vectors
    DOUBLE PRECISION, INTENT(OUT)           :: v(n,n)    ! Orthonormal basis
    INTEGER,          INTENT(OUT), OPTIONAL :: cond_flag ! Bad conditioning flag (0=OK)
    !=============================================================================================
    INTEGER :: j,jok,jj 
    DOUBLE PRECISION :: cc,cc1,cc2,epsi,vl 
    DOUBLE PRECISION :: ws(ndimx),vers(ndimx)
    LOGICAL :: ize 
    !=============================================================================================
    ! Inizialization
    v=0.d0 
    IF(PRESENT(cond_flag)) cond_flag=0
    ! Dimension check                                                       
    IF(n.GT.ndimx)then 
       WRITE(*,*)'graha: n =',n,' larger than ndimx=',ndimx,' in graha' 
       STOP 
    ENDIF
    ! Selection of the control for "zero" vectors                           
    cc1=SQRT(DOT_PRODUCT(a(1:n,1),a(1:n,1))) 
    cc2=SQRT(DOT_PRODUCT(a(1:n,2),a(1:n,2))) 
    epsi=1.d-12*MIN(cc1,cc2) 
    IF(epsi.EQ.0.d0)THEN 
       WRITE(*,*) ' a has rank zero' 
       IF(PRESENT(cond_flag)) cond_flag=1
       RETURN
    ENDIF
    ! start by orthonormalisation of the space spanned by the columns of a  
    !                                                                       
    ! V1 is the versor of A1                                                
    CALL versor(n,a(1:n,1),epsi,v(1:n,1),vl,ize) 
    IF(ize)THEN 
       WRITE(*,*) ' first vector of a is too small' 
       IF(PRESENT(cond_flag)) cond_flag=2
       RETURN
    ENDIF
    ! the following vectors are obtained                                    
    ! by removing the components along the previous ones, but first normalize
    CALL versor(n,a(1:n,2),epsi,vers(1:n),vl,ize)                    
    cc=-DOT_PRODUCT(v(1:n,1),vers(1:n)) 
    v(1:n,2)=vers(1:n)+cc*v(1:n,1) 
    CALL versor(n,v(1:n,2),epsi,v(1:n,2),vl,ize) 
    IF(ize)THEN 
       WRITE(*,*) ' a has practically rank one' 
       IF(PRESENT(cond_flag)) cond_flag=3
       RETURN
    ENDIF
    ! we now use the vectors of the canonic basis to supplement the span of 
    jok=0 
    DO 1 j=1,n 
       ! remove the components along span(A), that is along v1 and v2          
       ! ws=ej-<ej,v(1:n,1)>v(1:n,1)-<ej,v(1:n,2)>v(1:n,2)
       cc1=-v(j,1) 
       cc2=-v(j,2) 
       ws(1:n)=cc1*v(1:n,1)+cc2*v(1:n,2) 
       ws(j)=ws(j)+1.d0
       CALL versor(n,ws,epsi,v(1:n,3+jok),vl,ize) 
       IF(.NOT.ize)THEN 
          ! now v(3+jok) is orthogonal to span(A); remove the components along    
          ! the previous ones (unless it is the first)                            
          IF(jok.GT.0)THEN 
             DO  jj=1,jok 
                cc=-DOT_PRODUCT(v(1:n,3+jok),v(1:n,2+jj)) 
                v(1:n,3+jok)=v(1:n,3+jok)+ cc* v(1:n,2+jj)
             ENDDO
             CALL versor(n,v(1:n,3+jok),epsi,v(1:n,3+jok),vl,ize) 
             IF(ize)THEN 
                GOTO 1 
             ENDIF
          ENDIF
          ! the new versor is a good one                                          
          jok=jok+1 
          IF(jok.EQ.n-2)then 
             EXIT 
          ENDIF
       ENDIF
1   ENDDO
    IF(jok.LT.n-2)then 
       WRITE(*,*) ' graha: something went wrong, jok=',jok 
       IF(PRESENT(cond_flag)) cond_flag=4
       RETURN
    ENDIF
  END SUBROUTINE graha
! =======================================
! ===========================================================           
! SLINEL_Q                                                                
! semilinear boundary ellipse computation                               
! ===========================================================           
  SUBROUTINE slinel_q(dtpdet,gc,cc,ceicel,b,v,nd,bad_cond) 
! nd by 2 matrix with columns= gradients                                 
    REAL(KIND=qkind), INTENT(IN) :: dtpdet(nd,2) 
! normal and covariance matrices                                        
    REAL(KIND=qkind), INTENT(IN) :: gc(nd,nd),cc(nd,nd) 
! orthonormal basis                                                     
    REAL(KIND=qkind), INTENT(OUT):: v(nd,nd)
! inverted 2x2 matrix                                             
    REAL(KIND=qkind), INTENT(OUT):: b(2,2)
! regression formula
    REAL(KIND=qkind), INTENT(OUT):: ceicel(nd-2,2)
    INTEGER, INTENT(OUT), OPTIONAL :: bad_cond ! Flag for bad conditioning cases
! number of variables to be determined
    INTEGER, INTENT(IN) :: nd
! END INTERFACE
! orthonormal basis                                                     
    REAL(KIND=qkind) vt(nd,nd),gamv(nd,nd),cv(nd,nd),tmp(nd,nd) 
! partial matrices                                                      
    REAL(KIND=qkind) c4(nd-2,nd-2),cinv(nd-2,nd-2),c42(nd-2,2) 
! line of maximum variation                                             
    REAL(KIND=qkind) a(2,2),deta 
! loop indexes ii=1,2, ij,ijj=1,nd-2                                       
    INTEGER ii, ij, ijj 
! for inversion with tcholevski: workspace, error flag                  
    REAL(KIND=qkind) ws(nd-2) 
    INTEGER ierr 
    INTEGER cond_flag
! ===================================================================== 
! Initialization
    v=0.d0
    b=0.d0
    ceicel=0.d0
    IF(PRESENT(bad_cond)) bad_cond=0
! adapted orthonormal basis, covariance and normal matrix in the new bas
    CALL graha_q(dtpdet(1:nd,1:2),nd,v,cond_flag) 
    ! Check against bad conditioning 
    IF(cond_flag.NE.0)THEN
       WRITE(*,*) ' slinel_q: bad conditioning case, return!'
       IF(PRESENT(bad_cond)) bad_cond=cond_flag
       RETURN
    END IF
    vt=TRANSPOSE(v) 
    tmp=MATMUL(vt,gc) 
    gamv=MATMUL(tmp,v) 
    tmp=MATMUL(vt,cc)  
    cv=MATMUL(tmp,v)   
! ===================================================================== 
! (nd-2)x(nd-2) and (nd-2)x2 submatrices of normal matrix         
    DO  ijj=1,nd-2
       DO ij=1,nd-2 
          c4(ijj,ij)=cv(ijj+2,ij+2) ! C_hh
       ENDDO
       DO  ii=1,2 
          c42(ijj,ii)=cv(ijj+2,ii) ! C_hg
       ENDDO
    ENDDO
! ===========================================================           
! Cholewski method for inversion                                        
    CALL tchinv_q(c4,nd-2,cinv,ws,ierr) 
    IF(ierr.ne.0)THEN 
       WRITE(*,*) ' slinel_q: degenerate case after tchinv_q, ierr=',ierr 
    ENDIF
! ===========================================================           
! matrix to be used for out of plane component                          
    ceicel=MATMUL(cinv,c42) ! C_hh^-1 C_hg
! ===========================================================           
! linear map from the elements space (with base V) and the alpha-delta p
    a(1,1)=DOT_PRODUCT(dtpdet(1:nd,1),v(1:nd,1)) 
    a(1,2)=DOT_PRODUCT(dtpdet(1:nd,1),v(1:nd,2)) 
    a(2,1)=DOT_PRODUCT(dtpdet(1:nd,2),v(1:nd,1)) 
    a(2,2)=DOT_PRODUCT(dtpdet(1:nd,2),v(1:nd,2)) 
    CALL inv22_q(a,b,deta) 
  END SUBROUTINE slinel_q
! ====================================================================  
! Graham- Schmidt procedure to generate an orthonormal basis v          
! starting from 2  n-vectors a                                          
! The new basis must be such that the first 2 vectors are a basis       
! for the space spanned by the 2 columns of a                           
  SUBROUTINE graha_q(a,n,v,cond_flag)
    USE dyn_param, ONLY: ndimx 
    INTEGER, intent(in) :: n
    REAL(KIND=qkind), intent(in) :: a(n,2)
    REAL(KIND=qkind), intent(out) :: v(n,n) 
    INTEGER,          INTENT(OUT), OPTIONAL :: cond_flag ! Bad conditioning flag (0=OK)
! end interface
    integer j,jok,jj 
    REAL(KIND=qkind) cc,cc1,cc2,epsi,vl 
    REAL(KIND=qkind) ws(ndimx),vers(ndimx)
    logical ize 
! Inizialization
    v=0.d0 
    IF(PRESENT(cond_flag)) cond_flag=0
! dimension check                                                       
    if(n.gt.ndimx)then 
       write(*,*)'graha: n =',n,' larger than ndimx=',ndimx,' in graha' 
       stop 
    endif
! selection of the control for "zero" vectors                           
    cc1=sqrt(DOT_PRODUCT(a(1:n,1),a(1:n,1))) 
    cc2=sqrt(DOT_PRODUCT(a(1:n,2),a(1:n,2))) 
    epsi=1.d-20*min(cc1,cc2) 
    if(epsi.eq.0.d0)then 
       write(*,*)' a has rank zero' 
       IF(PRESENT(cond_flag)) cond_flag=1
       RETURN
!        stop                                                           
    endif
! start by orthonormalisation of the space spanned by the columns of a  
!                                                                       
! V1 is the versor of A1                                                
    call versor_q(n,a(1:n,1),epsi,v(1:n,1),vl,ize) 
    if(ize)then 
       write(*,*)' first vector of a is too small' 
       IF(PRESENT(cond_flag)) cond_flag=2
       RETURN
!        stop                                                           
    endif
! the following vectors are obtained                                    
! by removing the components along the previous ones, but first normalize
    CALL versor_q(n,a(1:n,2),epsi,vers(1:n),vl,ize)                    
    cc=-DOT_PRODUCT(v(1:n,1),vers(1:n)) 
    v(1:n,2)=vers(1:n)+cc*v(1:n,1) 
    call versor_q(n,v(1:n,2),epsi,v(1:n,2),vl,ize) 
    if(ize)then 
       write(*,*)' a has practically rank one' 
       IF(PRESENT(cond_flag)) cond_flag=3
       RETURN
!        stop                                                           
    endif
! we now use the vectors of the canonic basis to supplement the span of 
    jok=0 
    do 1 j=1,n 
! remove the components along span(A), that is along v1 and v2          
! ws=ej-<ej,v(1:n,1)>v(1:n,1)-<ej,v(1:n,2)>v(1:n,2)
       cc1=-v(j,1) 
       cc2=-v(j,2) 
       ws(1:n)=cc1*v(1:n,1)+cc2*v(1:n,2) 
       ws(j)=ws(j)+1.q0
       call versor_q(n,ws,epsi,v(1:n,3+jok),vl,ize) 
       if(.not.ize)then 
! now v(3+jok) is orthogonal to span(A); remove the components along    
! the previous ones (unless it is the first)                            
          if(jok.gt.0)then 
             do  jj=1,jok 
                cc=-DOT_PRODUCT(v(1:n,3+jok),v(1:n,2+jj)) 
                v(1:n,3+jok)=v(1:n,3+jok)+ cc* v(1:n,2+jj)
             enddo
             call versor_q(n,v(1:n,3+jok),epsi,v(1:n,3+jok),vl,ize) 
             if(ize)then 
                goto 1 
             endif
          endif
! the new versor is a good one                                          
          jok=jok+1 
          if(jok.eq.n-2)then 
             EXIT 
          endif
       endif
1   enddo
    if(jok.lt.n-2)then 
       write(*,*)'graha:  something went wrong, jok=',jok 
       IF(PRESENT(cond_flag)) cond_flag=4
       RETURN
    endif
  END SUBROUTINE graha_q
! =======================================

END MODULE semi_linear

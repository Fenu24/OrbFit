! ========LIBRARY semi_linear================== 
! routines shared by target_plane and by predict_obs                        
!             preob4 only for virimp                                    
!             linobs4  "                                                
!             ellip4   "                                                
!             elemo4   "                                                
!             slinel4  " 
! ===================================================================== 
! PREOB4- virtual impactor negative observation skyhole                 
! ============INTERFACE=================================================
SUBROUTINE preob4(el0,idsta,t1,unc0,gameq,                       &
     &    ceq,v6,sigma,npo,ibv,inl,al,de,hmagv,elm,                     &
     &    alpha,delta,hmagn,gamad,sig,axes,npo1)
  USE fund_const  
  USE pred_obs                      
  USE orbit_elements
  IMPLICIT NONE 
! ============= input ==================================================
! nominal initial conditions: epoch, elements  
  TYPE(orbit_elem), INTENT(IN) :: el0
! normal and covariance matrices                                        
  TYPE(orb_uncert), INTENT(IN) :: unc0
! sigmas for the boundary                                               
  DOUBLE PRECISION t1,gameq(6,6),ceq(6,6),sigma 
! number of points, flag for confidence bd/line of variation, nonlineari
  INTEGER npo,ibv,inl 
! magnitude                                                             
  DOUBLE PRECISION h,g,hmagn 
! station code                                                          
  INTEGER idsta 
! 6x6 rotation matrix giving the reference system adapted to the target 
  DOUBLE PRECISION v6(6,6) 
! ============= output =================================================
! points on the confidence boundary (difference w.r. to alpha,delta)    
! WARNING! the output number of points is npo1.le.npo;                  
! this beacuse hyperbolic points are discarded                          
  INTEGER npo1 
  DOUBLE PRECISION al(npo),de(npo),hmagv(npo) 
! line of elements                                                      
  DOUBLE PRECISION elm(6,npo) 
! best fit observations                                                 
  DOUBLE PRECISION alpha,delta 
! covariance                                                            
  DOUBLE PRECISION gamad(2,2),axes(2,2),sig(2) 
! ============END INTERFACE=============================================
! workspace                                                             
  DOUBLE PRECISION tmp(6,6),daddelt(6,2),v6t(6,6) 
  TYPE(orbit_elem) :: els
! partials in new reference system                                      
  DOUBLE PRECISION daddelt4(4,2),g4(4,4),gamv(6,6),c4(4,4),cv(6,6) 
  DOUBLE PRECISION ws(4) 
  INTEGER ii,jj,indp 
! partial derivatives of alpha, delta, w.r. to elements (by columns)    
  DOUBLE PRECISION daddet(6,2),dummy(6) 
! second derivatives of alpha, delta, w.r. to elements (not used)       
!  DOUBLE PRECISION ddade(6,6),dddde(6,6) 
! ===================================================================   
! orthonormal basis, matrix defining the plane of the ellipse           
  DOUBLE PRECISION v(4,4),ceicel(2,2) 
! transformation matrix between the two planes                          
  DOUBLE PRECISION b(2,2) 
! number of full revolutions around the sky                             
  INTEGER ng,nrev 
! functions                                                             
  DOUBLE PRECISION appmag
! elongation,distance to Earth, distance to Sun (to compute magnitude)  
  DOUBLE PRECISION pha,dis,rdot,dsun,elo,gallat
  DOUBLE PRECISION adot,ddot
  DOUBLE PRECISION adot0,ddot0
  DOUBLE PRECISION obs4(4), dobde(4,6) ! for alph_del2 
! ===================================================================   
! constant of gravitation, trigonometric constants from fund_const  
! temporaries, indexes                                                  
  DOUBLE PRECISION dal,ddl,maxsig,minsig 
  INTEGER n 
! flag for 2-body approximation; must be .false. for full n-body computa
  LOGICAL twobo 
  twobo=.false. 
!   static memory not required                                          
! ===================================================================== 
! compute observation; derivatives (of order 1) required                
  CALL  alph_del2 (el0,t1,idsta,obs4,1,6,dobde,   &
     &   TWOBO=twobo) 
!     &   pha,dis,rdot,dsun,elo,gallat,TWOBO=twobo) 
! store proper motion                                                   
  adot0=obs4(3) 
  ddot0=obs4(4)
! observations
  alpha=obs4(1)
  delta=obs4(2)      
! compute derivatives in the reference system adapted to the target plan
  v6t=TRANSPOSE(v6)
  daddet=TRANSPOSE(dobde(1:2,1:6)) 
  daddelt=MATMUL(v6t,daddet) 
! compute covariance matrix in the new reference system                 
  tmp=MATMUL(v6t,gameq) 
  gamv=MATMUL(tmp,v6) 
  tmp=MATMUL(v6t,ceq) 
  cv=MATMUL(tmp,v6) 
! reduce all to 4-d                                                     
  DO jj=1,4 
     DO ii=1,4 
!          g4(jj,ii)=gamv(jj+2,ii+2)                                    
        c4(jj,ii)=cv(jj+2,ii+2) 
     ENDDO
     DO ii=1,2 
        daddelt4(jj,ii)=daddelt(jj+2,ii) 
     ENDDO
  ENDDO
  CALL tchinv(c4,4,g4,ws,indp) 
!     WRITE(1,*)'daddelt4 ',daddelt4                                    
!     WRITE(1,*)'c4 ',c4                                                
!     WRITE(1,*)'g4 ',g4                                                
! compute ellipse of covariance of alpha,delta                          
  CALL ellip4(daddelt4,g4,sig,axes,gamad) 
!     write(1,*)' sig,axes,gamad ',sig,axes,gamad                       
! use nonlinear method                                                  
  IF(ibv.eq.0)THEN 
     maxsig=max(sig(1),sig(2)) 
     minsig=min(sig(1),sig(2)) 
     if(maxsig/minsig.le.200.d0)then 
        ibv=1 
     else 
        ibv=2 
     endif
  endif
! ===================================================================== 
! compute ellipse in the elements space                                 
  CALL slinel4(daddelt4,g4,c4,ceicel,b,v) 
!     WRITE(1,*)'ceicel,b,v4 ',ceicel,b,v                               
! ===========================================================           
! compute line of orbital elements                                      
  CALL linobs4(ibv,npo,el0%coord,axes,sig,b,v,sigma,ceicel,elm,v6,npo1) 
!     WRITE(1,*)' npo1',npo1                                            
!     DO jj=1,npo1                                                      
!        WRITE(1,*)'elm,jj=',jj,(elm(ii,jj),ii=1,6)                     
!     ENDDO                                                             
! ===========================================================           
  ng=0 
  DO 7 n=1,npo1 
! chose method to handle nonlinearity                                   
     IF(inl.eq.1)THEN 
! linear map from ellipse                                               
        dal=DOT_PRODUCT(elm(1:6,n),daddet(1:6,1)) 
        ddl=DOT_PRODUCT(elm(1:6,n),daddet(1:6,2)) 
        al(n)=dal 
        de(n)=ddl 
        elm(1:6,n)=el0%coord+elm(1:6,n) 
! apparent magnitude is the one of the nominal orbit                    
        hmagv(n)=hmagn 
     ELSEIF(inl.eq.2)THEN 
        write(*,*)' inl=2 is forbidden' 
        npo1=0 
        RETURN 
     ELSEIF(inl.eq.3)THEN 
! full n-body propagation from ellipse                                  
        elm(1:6,n)=el0%coord+elm(1:6,n) ! CALL vsumg(6,eq,elm(1,n),elm(1,n)) 
!           CALL alfdel (elm(1,n),t0,t1,idsta,al(n),de(n),               &
!     &          dummy,dummy,0,twobo,ddade,dddde)                        
        els=el0
        els%coord=elm(1:6,n)
        CALL  alph_del2 (els,t1,idsta,obs4,0,6,dobde,   &
     &         pha,dis,rdot,dsun,elo,gallat,TWOBO=twobo) 
        al(n)=obs4(1)-alpha 
        de(n)=obs4(2)-delta 
! compute apparent magnitude at time of observation                     
        hmagv(n)=appmag(h,g,dsun,dis,pha) 
     ELSE 
        WRITE(*,*)' preob4: this we have not invented yet ', inl 
        RETURN 
     ENDIF
! keep count of lost revolutions                                        
     IF(n.eq.1)THEN 
        IF(al(n).gt.pig)al(n)=al(n)-dpig 
     ELSE 
        CALL angupd(al(n),al(n-1),ng) 
     ENDIF
! temporary output                                                      
!       write(*,*)'Solution ',n,', RA/DEC (deg)',                       
!    +       al(n)*degrad,de(n)*degrad,ng                               
7 ENDDO
! ===================================================================== 
! ensure that LOV is consistent with nominal point                      
! first find midpoint of LOV, assume npo is even                        
  if(ibv.eq.2)then 
     nrev=nint((al(npo/2)+al(npo/2+1))/2.d0/dpig) 
!        write(*,*)'debug: nrev:',nrev                                  
     if(nrev.ne.0)then 
        do n=1,npo1 
           al(n)=al(n)-nrev*dpig 
        enddo
     endif
  endif
! restore original apparent motion                                      
  adot=adot0 
  ddot=ddot0 
END SUBROUTINE preob4
! ===========================================================           
! common subroutines for virmp                              
! A. Milani, 1999                              
! ===========================================================           
! LINOBS defines line of changes in orbital elements to be used for     
! confidence boundary/variations line                                   
! ===========================================================           
      SUBROUTINE linobs4(ibv,npo,eq,axes,sig,b,v,sigma,ceicel,elm,      &
     &  v6,npo1)
      USE fund_const                                                        
      IMPLICIT NONE 
! ====================INPUT==================================           
      INTEGER ibv,npo 
      DOUBLE PRECISION eq(6) 
      DOUBLE PRECISION axes(2,2),sig(2),b(2,2),sigma 
! matrix defining the plane of the ellipse,new orthonormal reference    
      DOUBLE PRECISION ceicel(2,2),v(4,4) 
! matrix defining the target plane adapted referemce system             
      DOUBLE PRECISION v6(6,6) 
! ===================OUTPUT==================================           
      INTEGER npo1 
      DOUBLE PRECISION elm(6,npo) 
! ==================END INTERFACE============================           
      INTEGER nn,n,i 
      DOUBLE PRECISION elm4(4) 
      DOUBLE PRECISION s,x,y,vad(2),xv,yv,dn,dth,theta,xa,yd 
      DOUBLE PRECISION eqnew(6) 
      DOUBLE PRECISION alde(2),ecc 
! ===================================================================== 
! line of maximum variation: in the alpha-delta plane                   
      DO i=1,2 
        vad(i)=axes(i,2)*sig(2) 
      ENDDO 
! in the elements space                                                 
      xv=(b(1,1)*vad(1)+b(1,2)*vad(2)) 
      yv=(b(2,1)*vad(1)+b(2,2)*vad(2)) 
!     WRITE(*,*)xv,yv                                                   
! direction not used any more                                           
!     theta0=atan2(yv,xv)                                               
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
           alde=xa*axes(1:2,1)+yd*axes(1:2,2) !
             ! CALL lincog(2,axes(1,1),xa,axes(1,2),yd,alde) 
! transfer of parametrisation in the V1,V2 plane                        
           x=(b(1,1)*alde(1)+b(1,2)*alde(2)) 
           y=(b(2,1)*alde(1)+b(2,2)*alde(2)) 
        ELSE 
           write(*,*)' linobs: this should not happen,ibv=',ibv 
        ENDIF 
! compute displacement on the confidence ellipsoid corresponding to x,y 
        nn=nn+1 
        CALL elemo4(x,y,v,ceicel,elm4) 
        elm(1:6,nn)=MATMUL(v6(1:6,3:6),elm4) 
            ! CALL mulmav(v6(1,3),6,4,elm4,4,elm(1,nn)) 
!       write(*,*)(elm(i,nn),i=1,6)                                     
! add to the original center of the ellipsoid of confidence             
        eqnew=eq+elm(1:6,nn) !CALL vsumg(6,eq,elm(1,nn),eqnew) 
        ecc=sqrt(eqnew(2)**2+eqnew(3)**2) 
        IF(ecc.ge.1.d0.or.eqnew(1).le.0.d0)THEN 
           write(*,*)' Hyperbolic, ecc=',ecc,' a=',eqnew(1) 
           nn=nn-1 
        ELSEIF(ecc.ge.0.99d0)THEN 
           write(*,*)' Almost Hyperbolic, ecc=',ecc,' a=',eqnew(1) 
           nn=nn-1 
        ENDIF 
    7 continue 
! final count of non hyperbolic orbits                                  
      npo1=nn 
      RETURN 
      END SUBROUTINE linobs4                                          
! ===================================================================== 
! ELLIPS                                                                
! compute covariance ellipse of two observables                         
! ===================================================================== 
      SUBROUTINE ellip4(daddet,gamm0,sig,axes,gamad) 
      IMPLICIT NONE 
! input covariance matrix                                               
      DOUBLE PRECISION gamm0(4,4) 
! input partial derivatives of alpha, delta, w.r. to elements (by column
      DOUBLE PRECISION daddet(4,2) 
! output covariance                                                     
      DOUBLE PRECISION gamad(2,2),axes(2,2),sig(2) 
! ==============END INTERFACE========================================== 
! eigenvalues, workspace, transposed                                    
      DOUBLE PRECISION eigval(2),tmp24(2,4),fv1(2),fv2(2),dadde(2,4) 
! loop indexes                                                          
      INTEGER i 
! error flag                                                            
      INTEGER ierr 
! ===================================================================== 
      dadde=TRANSPOSE(daddet) ! CALL transp(daddet,4,2,dadde) 
      tmp24=MATMUL(dadde,gamm0) ! CALL mulmat(dadde,2,4,gamm0,4,4,tmp24) 
      gamad=MATMUL(tmp24,daddet) ! CALL mulmat(tmp24,2,4,daddet,4,2,gamad) 
! ===================================================================== 
! compute ellipse of confidence                                         
! eigenvalues                                                           
      CALL rs(2,2,gamad,eigval,1,axes,fv1,fv2,ierr) 
      DO  i=1,2 
        IF(eigval(i).gt.0.d0)THEN 
           sig(i)=sqrt(eigval(i)) 
        ELSE 
           write(*,*) 'ellip4: non positive eigenvalue' 
           sig(i)=0.d0 
        ENDIF 
      ENDDO 
      RETURN 
      END SUBROUTINE ellip4                                          
! ===================================================================== 
! ELEMOV                                                                
! compute displacement on the confidence ellipsoid corresponding to x,y 
! on the plane of the gradients of alpha-delta                          
! ===================================================================== 
      SUBROUTINE elemo4(x,y,v,ceicel,del) 
      IMPLICIT NONE 
! inout/output                                                          
      DOUBLE PRECISION x,y,v(4,4),del(4) 
      DOUBLE PRECISION ceicel(2,2) 
! workspace                                                             
      DOUBLE PRECISION dee(2),deel(4) 
! ===================                                                   
      del=x*v(1:4,1)+y*v(1:4,2) ! CALL lincog(4,v(1,1),x,v(1,2),y,del) 
      dee=-x*ceicel(1:2,1)-y*ceicel(1:2,2) ! CALL lincog(2,ceicel(1,1),-x,ceicel(1,2),-y,dee) 
      deel=MATMUL(v(1:4,3:4),dee) ! CALL mulmav(v(1,3),4,2,dee,2,deel) 
      del=del+deel ! CALL vsumg(4,del,deel,del) 
      RETURN 
      END SUBROUTINE elemo4                                          
! ===========================================================           
! SLINEL                                                                
! semilinear boundary ellipse computation                               
! ===========================================================           
      SUBROUTINE slinel4(dtpdet,gc,cc,ceicel,b,v) 
        USE semi_linear, ONLY: graha
      IMPLICIT NONE 
      INTEGER ndim,ndimm2 
      PARAMETER (ndim=4,ndimm2=ndim-2) 
! 6 by 2 matrix with columns= gradients                                 
      DOUBLE PRECISION dtpdet(ndim,ndimm2) 
! normal and covariance matrices                                        
      DOUBLE PRECISION gc(ndim,ndim),cc(ndim,ndim) 
! orthonormal basis                                                     
      DOUBLE PRECISION v(ndim,ndim),vt(ndim,ndim) 
!       ,gamv(ndim,ndim)                                                
      DOUBLE PRECISION cv(ndim,ndim),tmp(ndim,ndim) 
! partial matrices                                                      
      DOUBLE PRECISION c4(ndimm2,ndimm2),cinv(ndimm2,ndimm2) 
      DOUBLE PRECISION c42(ndimm2,2),ceicel(ndimm2,2) 
! line of maximum variation                                             
      DOUBLE PRECISION a(2,2),b(2,2),deta 
! loop indexes ii=1,2, ij,ijj=1,ndim                                    
      INTEGER ii, ij, ijj 
! for inversion with tcholevski: workspace, error flag                  
      DOUBLE PRECISION ws(ndimm2) 
      INTEGER ierr 
! ===================================================================== 
! adapted orthonormal basis, covariance and normal matrix in the new bas
      CALL graha(dtpdet,ndim,v) 
      vt=TRANSPOSE(v) !CALL transp(v,ndim,ndim,vt) 
!     CALL mulmat(vt,ndim,ndim,gc,ndim,ndim,tmp)                        
!     CALL mulmat(tmp,ndim,ndim,v,ndim,ndim,gamv)                       
      tmp=MATMUL(vt,cc) ! CALL mulmat(vt,ndim,ndim,cc,ndim,ndim,tmp) 
      cv=MATMUL(tmp,v) ! CALL mulmat(tmp,ndim,ndim,v,ndim,ndim,cv) 
! ===================================================================== 
! 4x4 and 4x2 submatrices of normal matrix                              
      do 15 ijj=1,ndimm2 
        DO ij=1,ndimm2 
          c4(ijj,ij)=cv(ijj+2,ij+2) 
        ENDDO 
        DO  ii=1,2 
          c42(ijj,ii)=cv(ijj+2,ii) 
        ENDDO 
   15 continue 
! ===========================================================           
! Cholewski method for inversion                                        
      CALL tchinv(c4,ndimm2,cinv,ws,ierr) 
      IF(ierr.ne.0)THEN 
         WRITE(*,*) ' slinel4: degenerate case after tchinv, ierr=',ierr 
      ENDIF 
! ===========================================================           
! matrix to be used for out of plane component                          
      ceicel=MATMUL(cinv,c42) ! CALL mulmat(cinv,ndimm2,ndimm2,c42,ndimm2,2,ceicel) 
! ===========================================================           
! linear map from the elements space (with base V) and the alpha-delta p
      a(1,1)=DOT_PRODUCT(dtpdet(1:ndim,1),v(1:ndim,1)) 
      a(1,2)=DOT_PRODUCT(dtpdet(1:ndim,1),v(1:ndim,2)) 
      a(2,1)=DOT_PRODUCT(dtpdet(1:ndim,2),v(1:ndim,1)) 
      a(2,2)=DOT_PRODUCT(dtpdet(1:ndim,2),v(1:ndim,2)) 
      CALL inv22(a,b,deta) 
      RETURN 
    END SUBROUTINE slinel4                                          

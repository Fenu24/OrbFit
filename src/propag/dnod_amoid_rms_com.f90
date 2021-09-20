! *********************************************************
! COMPUTE the RMS of the NODAL DISTANCES and of the
! APPROXIMATED MOIDs (asc. and desc. node) using COM elements
! *********************************************************
! Written by GT Feb 2007
! *********************************************************
 SUBROUTINE dnod_amoid_rms_com(el1,el2,unc1,unc2,chk_der,dnodp,dnodm,&
      & dnodprms,dnodmrms,amoidp,amoidm,amoidprms,amoidmrms)
   USE orbit_elements
   USE output_control
   USE fund_const
   USE critical_points,ONLY : sign_dmin,mutual_ref
   IMPLICIT NONE
! =====================INTERFACE=======================
! Orbital elements with uncertainties
   TYPE(orbit_elem), INTENT(IN) :: el1,el2 ! orbital elements
   TYPE(orb_uncert), INTENT(IN) :: unc1,unc2
   LOGICAL, INTENT(IN) :: chk_der
! NODAL DISTANCES suffix p: ascending node suffix m: descending node
   DOUBLE PRECISION, INTENT(OUT) :: dnodp,dnodm
! Uncertainty of the nodal distances
   DOUBLE PRECISION, INTENT(OUT) :: dnodprms,dnodmrms
! Approximate MOID suffix p: ascending node suffix m: descending node 
   DOUBLE PRECISION, INTENT(OUT) :: amoidp,amoidm
! Uncertainty of the approximate MOID
   DOUBLE PRECISION, INTENT(OUT) :: amoidprms,amoidmrms
! =========================================================  
! ---------------------------------------------------------
! for the incremental ratio test
!    TYPE(orbit_elem) :: eq2var ! equinoctial elements
!    DOUBLE PRECISION :: dnodpvar,dnodmvar
!   DOUBLE PRECISION :: amoidpvar,amoidmvar
!    DOUBLE PRECISION :: eps_a2,eps_h2,eps_k2,eps_p2,eps_q2
! ---------------------------------------------------------
! 1: first mass (usually planet)  2:second mass (usually minor body)
   TYPE(orbit_elem) :: com1,com2  
   DOUBLE PRECISION :: q1,e1,i1,bigomega1,smallomega1  
   DOUBLE PRECISION :: q2,e2,i2,bigomega2,smallomega2  
! beta=SQRT(1-e**2)
   DOUBLE PRECISION :: beta1,beta2                     
! Mutual angles   
   DOUBLE PRECISION :: cosmutI,sinmutI                 
   DOUBLE PRECISION :: cmutom1,cmutom2                 
   DOUBLE PRECISION :: smutom1,smutom2
   DOUBLE PRECISION :: mutI,mutom1,mutom2
! Derivatives of sinmutI w.r.t. cosmutI
   DOUBLE PRECISION :: dsinmutIdcos
! Derivatives of cosmutI w.r.t. mutual angles
   DOUBLE PRECISION, DIMENSION(3,2) :: dcosmutI
! Derivatives of smutom w.r.t. cmutom
   DOUBLE PRECISION, DIMENSION(2) :: dsindcos
! Matrix of rotation from inertial to mutual reference frame
   DOUBLE PRECISION, DIMENSION(3,3):: rotinmut
! Derivatives of the matrix of rotation from inertial to 
! mutual reference frame
   DOUBLE PRECISION, DIMENSION(10,3,3):: drotinmut
! Derivatives of cosom1 and cosom2 (to compute the derivatives
! of the nodal distances)
   DOUBLE PRECISION, DIMENSION(2,3,2) :: dcosom 
   DOUBLE PRECISION, DIMENSION(10,10) :: covcom
! Derivatives of the nodal distances
   DOUBLE PRECISION :: ddnodpdcosom1,ddnodpdcosom2
   DOUBLE PRECISION :: ddnodpdq1,ddnodpde1,ddnodpdi1
   DOUBLE PRECISION :: ddnodpdbigom1,ddnodpdom1,ddnodpdom2
   DOUBLE PRECISION :: ddnodpdbigom2,ddnodpdq2,ddnodpde2,ddnodpdi2
   DOUBLE PRECISION :: ddnodmdcosom1,ddnodmdcosom2
   DOUBLE PRECISION :: ddnodmdq1,ddnodmde1,ddnodmdi1
   DOUBLE PRECISION :: ddnodmdbigom1,ddnodmdom1,ddnodmdom2
   DOUBLE PRECISION :: ddnodmdbigom2,ddnodmdq2,ddnodmde2,ddnodmdi2
! Matrices of transformation and covariance matrices
   DOUBLE PRECISION, DIMENSION(6,6) :: gel1,gel2
   DOUBLE PRECISION, DIMENSION(6,6) :: jaccomel1,jaccomel2
   DOUBLE PRECISION, DIMENSION(6,6) :: tjaccomel1,tjaccomel2
   DOUBLE PRECISION, DIMENSION(6,6) :: gamma1,gamma2
   DOUBLE PRECISION, DIMENSION(1,10) :: derdnodp
   DOUBLE PRECISION, DIMENSION(10,1) :: tderdnodp
   DOUBLE PRECISION, DIMENSION(1,1) :: dnodprmsvec
   DOUBLE PRECISION, DIMENSION(1,10) :: derdnodm
   DOUBLE PRECISION, DIMENSION(10,1) :: tderdnodm
   DOUBLE PRECISION, DIMENSION(1,1) :: dnodmrmsvec
! APPROXIMATE MOID (AMOID)
   DOUBLE PRECISION :: detAp,effe1p,effe2p,gi1p,gi2p
   DOUBLE PRECISION :: detAm,effe1m,effe2m,gi1m,gi2m
   DOUBLE PRECISION, DIMENSION(3) :: tau1p,tau2p,tau3p,tau3phat
   DOUBLE PRECISION, DIMENSION(3) :: newtau3p,newtau3m
   DOUBLE PRECISION, DIMENSION(10,3) :: dtau3m
   DOUBLE PRECISION :: modtau3p,modtau3m
   DOUBLE PRECISION, DIMENSION(3) :: tau1m,tau2m,tau3m,tau3mhat
   DOUBLE PRECISION :: x1barp,x2barp,x1barm,x2barm,k1p,k2p,k1m,k2m
   DOUBLE PRECISION, DIMENSION(3) :: xyz1p,xyz2p
   DOUBLE PRECISION, DIMENSION(3) :: xyz1m,xyz2m
   DOUBLE PRECISION :: dsignp,dsignm
   DOUBLE PRECISION, DIMENSION(3) :: deltaxyzp,deltaxyzm
   DOUBLE PRECISION :: ddeltaxdq1,ddeltaxdq2,ddeltaxde1,ddeltaxde2
   DOUBLE PRECISION :: dx1barpdcosmutompl,dx2barpdcosmutom
   DOUBLE PRECISION :: deffe1pdsmutompl,deffe2pdsmutom
   DOUBLE PRECISION :: dgi1pdcmutompl,dgi2pdcmutom
   DOUBLE PRECISION :: ddeltaxdi1,ddeltaxdi2
   DOUBLE PRECISION :: ddeltaxdbigom1,ddeltaxdbigom2
   DOUBLE PRECISION :: ddeltaxdom1,ddeltaxdom2
   DOUBLE PRECISION :: ddeltaydq1,ddeltaydq2,ddeltayde1,ddeltayde2
   DOUBLE PRECISION :: ddeltaydi1,ddeltaydi2
   DOUBLE PRECISION :: ddeltaydbigom1,ddeltaydbigom2
   DOUBLE PRECISION :: ddeltaydom1,ddeltaydom2
   DOUBLE PRECISION :: ddeltazdq1,ddeltazdq2,ddeltazde1,ddeltazde2
   DOUBLE PRECISION :: ddeltazdi1,ddeltazdi2
   DOUBLE PRECISION :: ddeltazdbigom1,ddeltazdbigom2
   DOUBLE PRECISION :: ddeltazdom1,ddeltazdom2
   DOUBLE PRECISION :: ddeltamxdq1,ddeltamxdq2,ddeltamxde1,ddeltamxde2
   DOUBLE PRECISION :: dx1barmdcosmutompl,dx2barmdcosmutom
   DOUBLE PRECISION :: deffe1mdsmutompl,deffe2mdsmutom
   DOUBLE PRECISION :: dgi1mdcmutompl,dgi2mdcmutom
   DOUBLE PRECISION :: ddeltamxdi1,ddeltamxdi2
   DOUBLE PRECISION :: ddeltamxdbigom1,ddeltamxdbigom2
   DOUBLE PRECISION :: ddeltamxdom1,ddeltamxdom2
   DOUBLE PRECISION :: ddeltamydq1,ddeltamydq2,ddeltamyde1,ddeltamyde2
   DOUBLE PRECISION :: ddeltamydi1,ddeltamydi2
   DOUBLE PRECISION :: ddeltamydbigom1,ddeltamydbigom2
   DOUBLE PRECISION :: ddeltamydom1,ddeltamydom2
   DOUBLE PRECISION :: ddeltamzdq1,ddeltamzdq2,ddeltamzde1,ddeltamzde2
   DOUBLE PRECISION :: ddeltamzdi1,ddeltamzdi2
   DOUBLE PRECISION :: ddeltamzdbigom1,ddeltamzdbigom2
   DOUBLE PRECISION :: ddeltamzdom1,ddeltamzdom2
   DOUBLE PRECISION, DIMENSION(10,3) :: ddelta,ddeltam
   DOUBLE PRECISION, DIMENSION(3) :: ddeltadkep,ddeltamdkep
   DOUBLE PRECISION, DIMENSION(3) :: deltarot,deltarotm
   DOUBLE PRECISION, DIMENSION(1,10) :: deramoidp,deramoidm
   DOUBLE PRECISION, DIMENSION(1,5) :: deramoidpdequ2,deramoidmdequ2
   DOUBLE PRECISION, DIMENSION(10,1) :: tderamoidp,tderamoidm
   DOUBLE PRECISION, DIMENSION(1,1) ::  amoidprmsvec,amoidmrmsvec
! auxiliary quantities
   INTEGER :: i,j,h,k,fail_flag     
   DOUBLE PRECISION :: vsize
   DOUBLE PRECISION :: s0,s1,s2,s3,s4,s5,s6,s7,s8,t0
! ============================================================
! Transformation to Cometary elements
   CALL coo_cha(el1,'COM',com1,fail_flag,jaccomel1)
   CALL coo_cha(el2,'COM',com2,fail_flag,jaccomel2)
! ------------------------------------------------------------
!    IF(chk_der) THEN
! to test the derivatives
!       OPEN(10,file='moidrms.opt',status='old')       
!       READ(10,*)eps_q2,eps_e2,eps_I2,eps_Om2,eps_omeg2
!       CLOSE(10)
!       com2var=com2
!       com2var%coord(1)=com2%coord(1)+eps_a2
!       com2var%coord(2)=com2%coord(2)+eps_h2
!       com2var%coord(3)=com2%coord(3)+eps_k2
!       com2var%coord(4)=com2%coord(4)+eps_p2
!       com2var%coord(5)=com2%coord(5)+eps_q2

!       CALL dnod_amoid(el1,el2var,dnodpvar,dnodmvar,amoidpvar,amoidmvar)
!       write(*,*)'amoidpvar,amoidmvar',amoidpvar,amoidmvar
!    ENDIF
! ------------------------------------------------------------    
   q1 = com1%coord(1) 
   e1 = com1%coord(2)
   i1 = com1%coord(3)
   bigomega1 = com1%coord(4)
   smallomega1 = com1%coord(5)
   q2 = com2%coord(1) 
   e2 = com2%coord(2) 
   i2 = com2%coord(3)
   bigomega2 = com2%coord(4)
   smallomega2 = com2%coord(5)
   beta1=sqrt(1-e1**2)
   beta2=sqrt(1-e2**2)
 ! Computation of mutual angles   
   CALL mutual_ref(com1,com2,mutI,mutom1,mutom2, &
        & dsindcos,dsinmutIdcos,dcosmutI,rotinmut,drotinmut,dcosom)
   cosmutI=cos(mutI)
   sinmutI=sin(mutI)
   cmutom1=cos(mutom1)
   smutom1=sin(mutom1)
   cmutom2=cos(mutom2)
   smutom2=sin(mutom2)

! ==================================================
! Computation of the ascending nodal distance dnodp
! ==================================================
   dnodp=(q2*(1+e2)/(1.d0+e2*cmutom2))- &
          & (q1*(1+e1)/(1.d0+e1*cmutom1))
! --------------------------------------------------------------------
! Computation of the derivatives of dnodp w.r.t. cosom1 and cosom2
! --------------------------------------------------------------------
   ddnodpdcosom1=q1*(1+e1)*e1/(1.d0+e1*cmutom1)**2
   ddnodpdcosom2=-q2*(1+e2)*e2/(1.d0+e2*cmutom2)**2
! --------------------------------------------------------------------
! Computation of the derivatives of dnodp w.r.t. the 10 elements
! --------------------------------------------------------------------
! Derivatives w.r.t. to planet elements
   ddnodpdq1=-(1+e1)/(1.d0+e1*cmutom1)
   ddnodpde1=(-q1*(1.d0+e1*cmutom1)+q1*(1+e1)*cmutom1)/ &
         & (1.d0+e1*cmutom1)**2
   ddnodpdi1= ddnodpdcosom1*dcosom(1,1,1)+ddnodpdcosom2*dcosom(2,1,1)
   ddnodpdbigom1=ddnodpdcosom1*dcosom(1,3,1)+ddnodpdcosom2*dcosom(2,3,1)
   ddnodpdom1=ddnodpdcosom1*dcosom(1,2,1)+ddnodpdcosom2*dcosom(2,2,1)
! Derivatives w.r.t. to asteroid/comet elements
   ddnodpdq2=(1+e2)/(1.d0+e2*cmutom2)
   ddnodpde2=(q2*(1.d0+e2*cmutom2)-q2*(1+e2)*cmutom2)/ &
        & (1.d0+e2*cmutom2)**2
   ddnodpdi2= ddnodpdcosom1*dcosom(1,1,2)+ddnodpdcosom2*dcosom(2,1,2)
   ddnodpdbigom2=ddnodpdcosom1*dcosom(1,3,2)+ddnodpdcosom2*dcosom(2,3,2)
   ddnodpdom2=ddnodpdcosom1*dcosom(1,2,2)+ddnodpdcosom2*dcosom(2,2,2)
! =====================================================================
! Construction of the covariance matrix in COM coordinates
   gel1=unc1%g(1:6,1:6)
   gel2=unc2%g(1:6,1:6)
   tjaccomel1=TRANSPOSE(jaccomel1)
   tjaccomel2=TRANSPOSE(jaccomel2)
! Covariance matrices
   gamma1=MATMUL(jaccomel1,MATMUL(gel1,tjaccomel1))
   gamma2=MATMUL(jaccomel2,MATMUL(gel2,tjaccomel2))
! Construction of matrix covcom
   covcom(1:10,1:10)=0.d0
   covcom(1:5,1:5)=gamma1(1:5,1:5)
   covcom(6:10,6:10)=gamma2(1:5,1:5)
! Derivatives of the ascending nodal distance
   derdnodp(1,1)=ddnodpdq1
   derdnodp(1,2)=ddnodpde1
   derdnodp(1,3)=ddnodpdi1
   derdnodp(1,4)=ddnodpdbigom1
   derdnodp(1,5)=ddnodpdom1
   derdnodp(1,6)=ddnodpdq2
   derdnodp(1,7)=ddnodpde2
   derdnodp(1,8)=ddnodpdi2
   derdnodp(1,9)=ddnodpdbigom2
   derdnodp(1,10)=ddnodpdom2
   tderdnodp=TRANSPOSE(derdnodp)
! RMS of the ascending nodal distance
   dnodprmsvec=MATMUL(derdnodp,MATMUL(covcom,tderdnodp)) 
   dnodprms=SQRT(dnodprmsvec(1,1))
! ===================================================
! Computation of the descending nodal distance dnodm
! ===================================================
   dnodm=(q2*(1+e2)/(1.d0-e2*cmutom2))- &
        & (q1*(1+e1)/(1.d0-e1*cmutom1))
! --------------------------------------------------------------------
! Computation of the derivatives of dnodp w.r.t. cosom1 and cosom2
! --------------------------------------------------------------------
   ddnodmdcosom1=-q1*(1+e1)*e1/(1.d0-e1*cmutom1)**2
   ddnodmdcosom2=q2*(1+e2)*e2/(1.d0-e2*cmutom2)**2
! --------------------------------------------------------------------
! Computation of the derivatives of dnodm w.r.t. the 10 elements
! --------------------------------------------------------------------
! Derivatives w.r.t. to planet elements
   ddnodmdq1=-(1+e1)/(1.d0-e1*cmutom1)
   ddnodmde1=-(q1*(1.d0-e1*cmutom1)+q1*(1+e1)*cmutom1)/ &
        & (1.d0-e1*cmutom1)**2
   ddnodmdi1= ddnodmdcosom1*dcosom(1,1,1)+ddnodmdcosom2*dcosom(2,1,1)
   ddnodmdbigom1=ddnodmdcosom1*dcosom(1,3,1)+ddnodmdcosom2*dcosom(2,3,1)
   ddnodmdom1=ddnodmdcosom1*dcosom(1,2,1)+ddnodmdcosom2*dcosom(2,2,1)
! Derivatives w.r.t. to asteroid/comet elements
   ddnodmdq2=(1+e2)/(1.d0-e2*cmutom2)
   ddnodmde2=(q2*(1.d0-e2*cmutom2)+q2*(1+e2)*cmutom2)/ &
        & (1.d0-e2*cmutom2)**2
   ddnodmdi2= ddnodmdcosom1*dcosom(1,1,2)+ddnodmdcosom2*dcosom(2,1,2)
   ddnodmdbigom2=ddnodmdcosom1*dcosom(1,3,2)+ddnodmdcosom2*dcosom(2,3,2)
   ddnodmdom2=ddnodmdcosom1*dcosom(1,2,2)+ddnodmdcosom2*dcosom(2,2,2)
! =====================================================================
! Derivatives of the descending nodal distance
   derdnodm(1,1)=ddnodmdq1
   derdnodm(1,2)=ddnodmde1
   derdnodm(1,3)=ddnodmdi1
   derdnodm(1,4)=ddnodmdbigom1
   derdnodm(1,5)=ddnodmdom1
   derdnodm(1,6)=ddnodmdq2
   derdnodm(1,7)=ddnodmde2
   derdnodm(1,8)=ddnodmdi2
   derdnodm(1,9)=ddnodmdbigom2
   derdnodm(1,10)=ddnodmdom2
   tderdnodm=TRANSPOSE(derdnodm)
! RMS of the descending nodal distance
   dnodmrmsvec=MATMUL(derdnodm,MATMUL(covcom,tderdnodm)) 
   dnodmrms=SQRT(dnodmrmsvec(1,1))
! =================================================================
! APPROXIMATE MOID
! =================================================================
! Approximate MOID at the ascending node
    effe1p=q1*e1*smutom1/(beta1*(1-e1))
    gi1p=q1*(1.d0+e1*cmutom1)/(beta1*(1-e1))
    effe2p=q2*e2*smutom2/(beta2*(1-e2))
    gi2p=q2*(1.d0+e2*cmutom2)/(beta2*(1-e2))
    detAp=effe2p**2*gi1p**2+effe1p**2*gi2p**2+gi1p**2*gi2p**2*sinmutI**2 &
         & -2.d0*effe1p*effe2p*gi1p*gi2p*cosmutI
    amoidp=ABS(dnodp)*SQRT(gi2p**2*gi1p**2*sinmutI**2/detAp)
! Computation of vectors tau1,tau2,tau3
    tau1p(1)=-effe1p
    tau1p(2)=gi1p
    tau1p(3)=0.d0
    tau2p(1)=-effe2p
    tau2p(2)=gi2p*cosmutI
    tau2p(3)=gi2p*sinmutI
    CALL prvec(tau1p,tau2p,tau3p)
    modtau3p=vsize(tau3p) 
    tau3phat(1:3)=tau3p(1:3)/modtau3p
! Computation of Delta x, Delta y and Delta z
    x1barp=q1*(1+e1)/(1.d0+e1*cmutom1)
    x2barp=q2*(1+e2)/(1.d0+e2*cmutom2)
    k1p=dnodp*gi2p*(effe2p*gi1p*cosmutI-effe1p*gi2p)/detAp
    k2p=dnodp*gi1p*(effe2p*gi1p-effe1p*gi2p*cosmutI)/detAp
    xyz1p(1)=x1barp-effe1p*k1p
    xyz1p(2)=gi1p*k1p
    xyz1p(3)=0.d0
    xyz2p(1)=x2barp-effe2p*k2p
    xyz2p(2)=gi2p*cosmutI*k2p
    xyz2p(3)=gi2p*sinmutI*k2p
    deltaxyzp(1)=xyz1p(1)-xyz2p(1)
    deltaxyzp(2)=xyz1p(2)-xyz2p(2)
    deltaxyzp(3)=xyz1p(3)-xyz2p(3)
    dx1barpdcosmutompl=-ddnodpdcosom1
    dx2barpdcosmutom=ddnodpdcosom2
    deffe1pdsmutompl=(q1/(1-e1))*e1/beta1
    deffe2pdsmutom=(q2/(1-e2))*e2/beta2
    dgi1pdcmutompl=deffe1pdsmutompl
    dgi2pdcmutom=deffe2pdsmutom
! give a sign to AMOIDp
     CALL sign_dmin(xyz1p,tau1p,xyz2p,tau2p,dsignp)
     amoidp=dsignp*amoidp  ! here amoidp is always >0
!------------------------------------------------------------------
! Derivatives of Delta x
!------------------------------------------------------------------
      t0 = (1.D0+e2*cmutom2)**2.D0*(1.D0+e1)*sinmutI**2*(1.D0+e1*cmutom1&
     &)/(e2**2*smutom2**2+2.D0*e2**2*smutom2**2*e1*cmutom1+e2**2*smutom2&
     &**2*e1**2*cmutom1**2+e1**2*smutom1**2+2.D0*e1**2*smutom1**2*e2*cmu&
     &tom2+e1**2*smutom1**2*e2**2*cmutom2**2+sinmutI**2+2.D0*sinmutI**2*&
     &e2*cmutom2+sinmutI**2*e2**2*cmutom2**2+2.D0*sinmutI**2*e1*cmutom1+&
     &4.D0*sinmutI**2*e1*cmutom1*e2*cmutom2+2.D0*sinmutI**2*e1*cmutom1*e&
     &2**2*cmutom2**2+sinmutI**2*e1**2*cmutom1**2+2.D0*sinmutI**2*e1**2*&
     &cmutom1**2*e2*cmutom2+sinmutI**2*e1**2*cmutom1**2*e2**2*cmutom2**2&
     &-2.D0*e1*smutom1*e2*smutom2*cosmutI-2.D0*e1*smutom1*e2**2*smutom2*&
     &cosmutI*cmutom2-2.D0*e1**2*smutom1*e2*smutom2*cosmutI*cmutom1-2.D0&
     &*e1**2*smutom1*e2**2*smutom2*cosmutI*cmutom1*cmutom2)             
      ddeltaxdq1=t0 
!                                                                        
      t0 = -(1.D0+e2)*sinmutI**2*(1.D0+e1*cmutom1)**2.D0*(1.D0+e2*cmutom&
     &2)/(e2**2*smutom2**2+2.D0*e2**2*smutom2**2*e1*cmutom1+e2**2*smutom&
     &2**2*e1**2*cmutom1**2+e1**2*smutom1**2+2.D0*e1**2*smutom1**2*e2*cm&
     &utom2+e1**2*smutom1**2*e2**2*cmutom2**2+sinmutI**2+2.D0*sinmutI**2&
     &*e2*cmutom2+sinmutI**2*e2**2*cmutom2**2+2.D0*sinmutI**2*e1*cmutom1&
     &+4.D0*sinmutI**2*e1*cmutom1*e2*cmutom2+2.D0*sinmutI**2*e1*cmutom1*&
     &e2**2*cmutom2**2+sinmutI**2*e1**2*cmutom1**2+2.D0*sinmutI**2*e1**2&
     &*cmutom1**2*e2*cmutom2+sinmutI**2*e1**2*cmutom1**2*e2**2*cmutom2**&
     &2-2.D0*e1*smutom1*e2*smutom2*cosmutI-2.D0*e1*smutom1*e2**2*smutom2&
     &*cosmutI*cmutom2-2.D0*e1**2*smutom1*e2*smutom2*cosmutI*cmutom1-2.D&
     &0*e1**2*smutom1*e2**2*smutom2*cosmutI*cmutom1*cmutom2)            
      ddeltaxdq2=t0 
 !                                                                        
      s1 = -(1.D0+e2*cmutom2)**2.D0 
      s5 = -q1*e2**2*smutom2**2*e1**2*cmutom1**2+2.D0*q1*e1**2*smutom1**&
     &2*e2*cmutom2+2.D0*q1*cmutom1**3*sinmutI**2*e1**2*e2*cmutom2+q1*cmu&
     &tom1**3*sinmutI**2*e1**2*e2**2*cmutom2**2+2.D0*q1*cmutom1**2*sinmu&
     &tI**2*e1*e2**2*cmutom2**2+q1*e1**2*smutom1**2*e2**2*cmutom2**2+q1*&
     &cmutom1*e1**2*smutom1**2*e2**2*cmutom2**2+2.D0*q1*cmutom1*e1**2*sm&
     &utom1**2*e2*cmutom2+q1*cmutom1**3*e2**2*smutom2**2*e1**2+2.D0*q1*c&
     &mutom1**2*e2**2*smutom2**2*e1+4.D0*q1*cmutom1**2*sinmutI**2*e1*e2*&
     &cmutom2-4.D0*q1*cmutom1*e1*smutom1*e2*smutom2*cosmutI+q1*cmutom1*s&
     &inmutI**2                                                         
      s4 = s5+2.D0*q1*cmutom1*sinmutI**2*e2*cmutom2+q1*cmutom1*sinmutI**&
     &2*e2**2*cmutom2**2-2.D0*q1*cmutom1**2*e1**2*smutom1*e2**2*smutom2*&
     &cosmutI*cmutom2-2.D0*q1*cmutom1**2*e1**2*smutom1*e2*smutom2*cosmut&
     &I-4.D0*q1*cmutom1*e1*smutom1*e2**2*smutom2*cosmutI*cmutom2-2.D0*q2&
     &*e1*smutom1**2-q1*e2**2*smutom2**2+2.D0*q1*e1*smutom1**2-2.D0*q2*e&
     &1**2*cmutom1*smutom1**2+2.D0*q1*cmutom1**2*sinmutI**2*e1-q1*sinmut&
     &I**2*e1**2*cmutom1**2+q1*cmutom1*e2**2*smutom2**2+q1*cmutom1**3*si&
     &nmutI**2*e1**2                                                    
      s5 = s4-q1*sinmutI**2-2.D0*q2*e1*smutom1**2*e2*cmutom2+4.D0*q1*e1*&
     &smutom1**2*e2*cmutom2-2.D0*q2*e2*e1**2*cmutom1*smutom1**2-2.D0*q2*&
     &e1**2*cmutom1*smutom1**2*e2*cmutom2-2.D0*q2*e2*e1*smutom1**2-2.D0*&
     &q2*e2**2*e1*smutom1**2*cmutom2+2.D0*smutom1*q2*e2*smutom2*cosmutI+&
     &2.D0*smutom1*q2*e1**2*cmutom1**2*e2*smutom2*cosmutI-2.D0*smutom1*q&
     &1*e2**2*smutom2*cosmutI*cmutom2+2.D0*smutom1*q2*e2**2*smutom2*cosm&
     &utI-2.D0*smutom1*q1*e2*smutom2*cosmutI+4.D0*smutom1*q2*e2*smutom2*&
     &cosmutI*e1*cmutom1                                                
      s3 = s5-2.D0*q2*e2**2*e1**2*cmutom1*smutom1**2*cmutom2+2.D0*q1*e2*&
     &*2*cmutom2**2*e1*smutom1**2-2.D0*q1*sinmutI**2*e2*cmutom2+2.D0*smu&
     &tom1*q2*e2**2*e1**2*cmutom1**2*smutom2*cosmutI-q1*sinmutI**2*e2**2&
     &*cmutom2**2-4.D0*q1*sinmutI**2*e1*cmutom1*e2*cmutom2+4.D0*smutom1*&
     &q2*e2**2*smutom2*cosmutI*e1*cmutom1-q1*sinmutI**2*e1**2*cmutom1**2&
     &*e2**2*cmutom2**2-2.D0*q1*sinmutI**2*e1*cmutom1+q1*e1**2*smutom1**&
     &2-2.D0*q1*e2**2*smutom2**2*e1*cmutom1+q1*cmutom1*e1**2*smutom1**2-&
     &2.D0*q1*sinmutI**2*e1**2*cmutom1**2*e2*cmutom2-2.D0*q1*sinmutI**2*&
     &e1*cmutom1*e2**2*cmutom2**2                                       
      s4 = sinmutI**2/(e2**2*smutom2**2+2.D0*e2**2*smutom2**2*e1*cmutom1&
     &+e2**2*smutom2**2*e1**2*cmutom1**2+e1**2*smutom1**2+2.D0*e1**2*smu&
     &tom1**2*e2*cmutom2+e1**2*smutom1**2*e2**2*cmutom2**2+sinmutI**2+2.&
     &D0*sinmutI**2*e2*cmutom2+sinmutI**2*e2**2*cmutom2**2+2.D0*sinmutI*&
     &*2*e1*cmutom1+4.D0*sinmutI**2*e1*cmutom1*e2*cmutom2+2.D0*sinmutI**&
     &2*e1*cmutom1*e2**2*cmutom2**2+sinmutI**2*e1**2*cmutom1**2+2.D0*sin&
     &mutI**2*e1**2*cmutom1**2*e2*cmutom2+sinmutI**2*e1**2*cmutom1**2*e2&
     &**2*cmutom2**2-2.D0*e1*smutom1*e2*smutom2*cosmutI-2.D0*e1*smutom1*&
     &e2**2*smutom2*cosmutI*cmutom2-2.D0*e1**2*smutom1*e2*smutom2*cosmut&
     &I*cmutom1-2.D0*e1**2*smutom1*e2**2*smutom2*cosmutI*cmutom1*cmutom2&
     &)**2.D0                                                           
      s2 = s3*s4 
      ddeltaxde1 = s1*s2 
  !                                                                        
      s2 = (1.D0+e1*cmutom1)**2.D0 
      s5 = -q2*sinmutI**2*e1**2*cmutom1**2-2.D0*q2*sinmutI**2*e2*cmutom2&
     &-2.D0*q2*sinmutI**2*e1*cmutom1-q2*sinmutI**2*e2**2*cmutom2**2+2.D0&
     &*q2*cmutom2**2*sinmutI**2*e2+q2*cmutom2**3*sinmutI**2*e2**2+e1**2*&
     &smutom1**2*q2*cmutom2-2.D0*q1*e2**2*smutom2**2*cmutom2+q2*cmutom2*&
     &e2**2*smutom2**2-2.D0*q1*e1*e2*smutom2**2+4.D0*q2*e2*smutom2**2*e1&
     &*cmutom1+2.D0*q2*e1**2*cmutom1**2*e2*smutom2**2-2.D0*q1*e1**2*e2*s&
     &mutom2**2*cmutom1                                                 
      s4 = s5+2.D0*e1**2*smutom1*q1*smutom2*cosmutI+4.D0*e1**2*smutom1*c&
     &mutom2*q1*e2*smutom2*cosmutI-2.D0*e1**2*smutom1*q2*smutom2*cosmutI&
     &*cmutom1+2.D0*e1*smutom1*q1*smutom2*cosmutI+2.D0*e1**2*smutom1*cmu&
     &tom2**2*q1*e2**2*smutom2*cosmutI+2.D0*e1*smutom1*cmutom2**2*q1*e2*&
     &*2*smutom2*cosmutI+4.D0*e1*smutom1*cmutom2*q1*e2*smutom2*cosmutI-2&
     &.D0*e1**2*smutom1**2*q2*e2*cmutom2-e1**2*smutom1**2*q2*e2**2*cmuto&
     &m2**2-4.D0*e1**2*smutom1*q2*cmutom2*e2*smutom2*cosmutI*cmutom1-2.D&
     &0*e1*smutom1*q2*cmutom2**2*e2**2*smutom2*cosmutI-2.D0*e1**2*smutom&
     &1*q2*cmutom2**2*e2**2*smutom2*cosmutI*cmutom1+2.D0*e1**2*smutom1**&
     &2*q2*e2*cmutom2**2                                                
      s5 = s4+e1**2*smutom1**2*q2*e2**2*cmutom2**3-q2*sinmutI**2-2.D0*q1&
     &*e1*e2**2*smutom2**2*cmutom2+q2*cmutom2**3*sinmutI**2*e1**2*cmutom&
     &1**2*e2**2+4.D0*q2*cmutom2**2*sinmutI**2*e1*cmutom1*e2+2.D0*q2*cmu&
     &tom2**3*sinmutI**2*e1*cmutom1*e2**2+q2*e2**2*smutom2**2*e1**2*cmut&
     &om1**2+2.D0*q2*e2**2*smutom2**2*e1*cmutom1-2.D0*q2*sinmutI**2*e1**&
     &2*cmutom1**2*e2*cmutom2-q2*sinmutI**2*e1**2*cmutom1**2*e2**2*cmuto&
     &m2**2+q2*cmutom2*e2**2*smutom2**2*e1**2*cmutom1**2+2.D0*q2*cmutom2&
     &*e2**2*smutom2**2*e1*cmutom1-2.D0*q2*sinmutI**2*e1*cmutom1*e2**2*c&
     &mutom2**2                                                         
      s3 = s5+2.D0*q2*cmutom2**2*sinmutI**2*e1**2*cmutom1**2*e2-2.D0*q1*&
     &e2*smutom2**2*e1*cmutom1-2.D0*q1*e2**2*smutom2**2*e1*cmutom1*cmuto&
     &m2+q2*cmutom2*sinmutI**2*e1**2*cmutom1**2+2.D0*q2*cmutom2*sinmutI*&
     &*2*e1*cmutom1-4.D0*q2*sinmutI**2*e1*cmutom1*e2*cmutom2-2.D0*e1*smu&
     &tom1*q2*smutom2*cosmutI-2.D0*q1*e1**2*e2**2*smutom2**2*cmutom1*cmu&
     &tom2-4.D0*e1*smutom1*q2*cmutom2*e2*smutom2*cosmutI-e1**2*smutom1**&
     &2*q2+q2*cmutom2*sinmutI**2-2.D0*q1*e2*smutom2**2+2.D0*q2*e2*smutom&
     &2**2+q2*e2**2*smutom2**2                                          
      s1 = s2*s3 
      s2 = sinmutI**2/(e2**2*smutom2**2+2.D0*e2**2*smutom2**2*e1*cmutom1&
     &+e2**2*smutom2**2*e1**2*cmutom1**2+e1**2*smutom1**2+2.D0*e1**2*smu&
     &tom1**2*e2*cmutom2+e1**2*smutom1**2*e2**2*cmutom2**2+sinmutI**2+2.&
     &D0*sinmutI**2*e2*cmutom2+sinmutI**2*e2**2*cmutom2**2+2.D0*sinmutI*&
     &*2*e1*cmutom1+4.D0*sinmutI**2*e1*cmutom1*e2*cmutom2+2.D0*sinmutI**&
     &2*e1*cmutom1*e2**2*cmutom2**2+sinmutI**2*e1**2*cmutom1**2+2.D0*sin&
     &mutI**2*e1**2*cmutom1**2*e2*cmutom2+sinmutI**2*e1**2*cmutom1**2*e2&
     &**2*cmutom2**2-2.D0*e1*smutom1*e2*smutom2*cosmutI-2.D0*e1*smutom1*&
     &e2**2*smutom2*cosmutI*cmutom2-2.D0*e1**2*smutom1*e2*smutom2*cosmut&
     &I*cmutom1-2.D0*e1**2*smutom1*e2**2*smutom2*cosmutI*cmutom1*cmutom2&
     &)**2.D0                                                           
      ddeltaxde2 = s1*s2 
!                                                                     
    ddeltaxdi1=dx1barpdcosmutompl*dcosom(1,1,1)-dx2barpdcosmutom*dcosom(2,1,1)+ &
         & k2p*deffe2pdsmutom*dsindcos(2)*dcosom(2,1,1)-k1p* &
         & deffe1pdsmutompl*dsindcos(1)*dcosom(1,1,1)
    ddeltaxdi2=dx1barpdcosmutompl*dcosom(1,1,2)-dx2barpdcosmutom*dcosom(2,1,2)+ &
         & k2p*deffe2pdsmutom*dsindcos(2)*dcosom(2,1,2)-k1p* &
         & deffe1pdsmutompl*dsindcos(1)*dcosom(1,1,2)
!
    ddeltaxdbigom1=dx1barpdcosmutompl*dcosom(1,3,1)-dx2barpdcosmutom* &
         & dcosom(2,3,1)+k2p*deffe2pdsmutom*dsindcos(2)*dcosom(2,3,1)-k1p* &
         & deffe1pdsmutompl*dsindcos(1)*dcosom(1,3,1)
    ddeltaxdbigom2=dx1barpdcosmutompl*dcosom(1,3,2)-dx2barpdcosmutom* &
         & dcosom(2,3,2)+k2p*deffe2pdsmutom*dsindcos(2)*dcosom(2,3,2)-k1p* &
         & deffe1pdsmutompl*dsindcos(1)*dcosom(1,3,2)
!
    ddeltaxdom1=dx1barpdcosmutompl*dcosom(1,2,1)-dx2barpdcosmutom* &
         & dcosom(2,2,1)+k2p*deffe2pdsmutom*dsindcos(2)*dcosom(2,2,1)-k1p* &
         & deffe1pdsmutompl*dsindcos(1)*dcosom(1,2,1)
    ddeltaxdom2=dx1barpdcosmutompl*dcosom(1,2,2)-dx2barpdcosmutom* &
         & dcosom(2,2,2)+k2p*deffe2pdsmutom*dsindcos(2)*dcosom(2,2,2)-k1p* &
         & deffe1pdsmutompl*dsindcos(1)*dcosom(1,2,2)
!------------------------------------------------------------------
! Derivatives of Delta y
!------------------------------------------------------------------
      t0 = -(1.D0+e2*cmutom2)**2.D0*(1.D0+e1)*e1*smutom1*(cosmutI-1.D0)*&
     &(cosmutI+1.D0)/(e2**2*smutom2**2+2.D0*e2**2*smutom2**2*e1*cmutom1+&
     &e2**2*smutom2**2*e1**2*cmutom1**2+e1**2*smutom1**2+2.D0*e1**2*smut&
     &om1**2*e2*cmutom2+e1**2*smutom1**2*e2**2*cmutom2**2+sinmutI**2+2.D&
     &0*sinmutI**2*e2*cmutom2+sinmutI**2*e2**2*cmutom2**2+2.D0*sinmutI**&
     &2*e1*cmutom1+4.D0*sinmutI**2*e1*cmutom1*e2*cmutom2+2.D0*sinmutI**2&
     &*e1*cmutom1*e2**2*cmutom2**2+sinmutI**2*e1**2*cmutom1**2+2.D0*sinm&
     &utI**2*e1**2*cmutom1**2*e2*cmutom2+sinmutI**2*e1**2*cmutom1**2*e2*&
     &*2*cmutom2**2-2.D0*e1*smutom1*e2*smutom2*cosmutI-2.D0*e1*smutom1*e&
     &2**2*smutom2*cosmutI*cmutom2-2.D0*e1**2*smutom1*e2*smutom2*cosmutI&
     &*cmutom1-2.D0*e1**2*smutom1*e2**2*smutom2*cosmutI*cmutom1*cmutom2)
      ddeltaydq1=t0 
!                                                                        
      t0 = (1.D0+e1*cmutom1)*(1.D0+e2)*e1*smutom1*(cosmutI-1.D0)*(cosmut&
     &I+1.D0)*(1.D0+e2*cmutom2)/(e2**2*smutom2**2+2.D0*e2**2*smutom2**2*&
     &e1*cmutom1+e2**2*smutom2**2*e1**2*cmutom1**2+e1**2*smutom1**2+2.D0&
     &*e1**2*smutom1**2*e2*cmutom2+e1**2*smutom1**2*e2**2*cmutom2**2+sin&
     &mutI**2+2.D0*sinmutI**2*e2*cmutom2+sinmutI**2*e2**2*cmutom2**2+2.D&
     &0*sinmutI**2*e1*cmutom1+4.D0*sinmutI**2*e1*cmutom1*e2*cmutom2+2.D0&
     &*sinmutI**2*e1*cmutom1*e2**2*cmutom2**2+sinmutI**2*e1**2*cmutom1**&
     &2+2.D0*sinmutI**2*e1**2*cmutom1**2*e2*cmutom2+sinmutI**2*e1**2*cmu&
     &tom1**2*e2**2*cmutom2**2-2.D0*e1*smutom1*e2*smutom2*cosmutI-2.D0*e&
     &1*smutom1*e2**2*smutom2*cosmutI*cmutom2-2.D0*e1**2*smutom1*e2*smut&
     &om2*cosmutI*cmutom1-2.D0*e1**2*smutom1*e2**2*smutom2*cosmutI*cmuto&
     &m1*cmutom2)                                                       
      ddeltaydq2=t0 
!
     s1 = smutom1*(cosmutI-1.D0)*(cosmutI+1.D0) 
      s3 = 1.D0+e2*cmutom2 
      s7 = q2*e2*e1**2*cmutom1**2*sinmutI**2+q2*sinmutI**2-q1*sinmutI**2&
     &+2.D0*cmutom1*q2*e2**3*e1*sinmutI**2*cmutom2**2+2.D0*q2*e2**2*e1**&
     &2*cmutom1**2*sinmutI**2*cmutom2+2.D0*cmutom1*q2*e2*e1*sinmutI**2+q&
     &2*e2**3*e1**2*cmutom1**2*sinmutI**2*cmutom2**2+q2*e2**3*e1**2*smut&
     &om2**2*cmutom1**2-2.D0*q1*e1*e2**3*cmutom2**3*sinmutI**2+4.D0*q1*e&
     &2**2*cmutom2*smutom2*cosmutI*e1**2*smutom1+2.D0*q1*e2**3*cmutom2**&
     &2*smutom2*cosmutI*e1**2*smutom1-2.D0*cmutom1*q1*e2*smutom2*cosmutI&
     &*e1**2*smutom1-q2*e2*smutom1**2*e1**2-q1*e2**3*cmutom2*smutom2**2+&
     &2.D0*q2*e2**2*sinmutI**2*cmutom2+q2*e2**3*sinmutI**2*cmutom2**2   
      s6 = s7-q1*e2**3*cmutom2**3*sinmutI**2-3.D0*q1*sinmutI**2*e2**2*cm&
     &utom2**2-3.D0*q1*sinmutI**2*e2*cmutom2+q1*sinmutI**2*e1**2*cmutom1&
     &**2+q1*e1**2*e2**3*cmutom2**3*sinmutI**2*cmutom1**2+q1*e1**2*e2**3&
     &*cmutom2*smutom2**2*cmutom1**2-6.D0*q1*e1**2*e2**2*cmutom2**2*sinm&
     &utI**2*cmutom1-6.D0*q1*e1**2*e2*cmutom2*sinmutI**2*cmutom1-2.D0*q1&
     &*e1**2*e2**3*cmutom2**3*sinmutI**2*cmutom1-6.D0*q1*e1*e2**2*cmutom&
     &2**2*sinmutI**2-2.D0*q1*e1*e2**3*cmutom2*smutom2**2-6.D0*q1*e1*e2*&
     &cmutom2*sinmutI**2-2.D0*q1*e1**2*e2**3*cmutom2*smutom2**2*cmutom1-&
     &2.D0*q1*e1**2*e2**2*smutom2**2*cmutom1+2.D0*q1*e2*smutom2*cosmutI*&
     &e1**2*smutom1+q2*e2**2*smutom2**2*e1**2*cmutom1**2                
      s7 = 2.D0*q2*e2**2*smutom2**2*e1*cmutom1+q2*sinmutI**2*e1**2*cmuto&
     &m1**2*e2**2*cmutom2**2-e1**2*smutom1**2*q2*e2**2*cmutom2**2+4.D0*q&
     &2*sinmutI**2*e1*cmutom1*e2*cmutom2+2.D0*q2*sinmutI**2*e1*cmutom1*e&
     &2**2*cmutom2**2+2.D0*q2*sinmutI**2*e1**2*cmutom1**2*e2*cmutom2-2.D&
     &0*e1**2*smutom1**2*q2*e2*cmutom2+q2*sinmutI**2*e1**2*cmutom1**2+2.&
     &D0*q2*sinmutI**2*e2*cmutom2+2.D0*q2*sinmutI**2*e1*cmutom1+q2*sinmu&
     &tI**2*e2**2*cmutom2**2-2.D0*q1*e1*e2**2*smutom2**2-2.D0*q1*e1**2*s&
     &inmutI**2*cmutom1+3.D0*q1*e1**2*smutom1**2*e2**2*cmutom2**2+3.D0*q&
     &1*e1**2*smutom1**2*e2*cmutom2+q1*e2**2*smutom2**2*e1**2*cmutom1**2
      s5 = s7+3.D0*q1*sinmutI**2*e1**2*cmutom1**2*e2*cmutom2+s6+3.D0*q1*&
     &sinmutI**2*e1**2*cmutom1**2*e2**2*cmutom2**2+q1*e2**3*cmutom2**3*s&
     &mutom1**2*e1**2-2.D0*q2*e2**2*smutom1**2*e1**2*cmutom2-q2*e2**3*sm&
     &utom1**2*e1**2*cmutom2**2+2.D0*cmutom1*q2*e2**3*e1*smutom2**2+4.D0&
     &*cmutom1*q2*e2**2*e1*sinmutI**2*cmutom2-4.D0*cmutom1*q1*e2**2*cmut&
     &om2*smutom2*cosmutI*e1**2*smutom1-2.D0*cmutom1*q1*e2**3*cmutom2**2&
     &*smutom2*cosmutI*e1**2*smutom1-e1**2*smutom1**2*q2-2.D0*q1*e1*sinm&
     &utI**2+q2*e2**2*smutom2**2+q2*e2**3*smutom2**2+q2*e2*sinmutI**2-q1&
     &*e2**2*smutom2**2+q1*e1**2*smutom1**2                             
      s6 = 1.D0/((e2**2*smutom2**2+2.D0*e2**2*smutom2**2*e1*cmutom1+e2**&
     &2*smutom2**2*e1**2*cmutom1**2+e1**2*smutom1**2+2.D0*e1**2*smutom1*&
     &*2*e2*cmutom2+e1**2*smutom1**2*e2**2*cmutom2**2+sinmutI**2+2.D0*si&
     &nmutI**2*e2*cmutom2+sinmutI**2*e2**2*cmutom2**2+2.D0*sinmutI**2*e1&
     &*cmutom1+4.D0*sinmutI**2*e1*cmutom1*e2*cmutom2+2.D0*sinmutI**2*e1*&
     &cmutom1*e2**2*cmutom2**2+sinmutI**2*e1**2*cmutom1**2+2.D0*sinmutI*&
     &*2*e1**2*cmutom1**2*e2*cmutom2+sinmutI**2*e1**2*cmutom1**2*e2**2*c&
     &mutom2**2-2.D0*e1*smutom1*e2*smutom2*cosmutI-2.D0*e1*smutom1*e2**2&
     &*smutom2*cosmutI*cmutom2-2.D0*e1**2*smutom1*e2*smutom2*cosmutI*cmu&
     &tom1-2.D0*e1**2*smutom1*e2**2*smutom2*cosmutI*cmutom1*cmutom2)**2.&
     &D0)                                                               
      s4 = s5*s6 
      s2 = s3*s4 
      ddeltayde1 = s1*s2 
 !                                                                       
      s1 = -smutom1*e1*(cosmutI-1.D0) 
      s3 = (cosmutI+1.D0)*(1.D0+e1*cmutom1) 
      s7 = -e1**2*smutom1**2*q2+q2*cmutom2*sinmutI**2-2.D0*q1*e2*smutom2&
     &**2+2.D0*q2*e2*smutom2**2+q2*e2**2*smutom2**2-q2*sinmutI**2-4.D0*e&
     &1*smutom1*q2*cmutom2*e2*smutom2*cosmutI+e1**2*smutom1**2*q2*cmutom&
     &2+q2*cmutom2*e2**2*smutom2**2+e1**2*smutom1**2*q2*e2**2*cmutom2**3&
     &+2.D0*e1**2*smutom1**2*q2*e2*cmutom2**2-e1**2*smutom1**2*q2*e2**2*&
     &cmutom2**2-2.D0*e1**2*smutom1**2*q2*e2*cmutom2                    
      s6 = s7-2.D0*e1*smutom1*q2*cmutom2**2*e2**2*smutom2*cosmutI-4.D0*e&
     &1**2*smutom1*q2*cmutom2*e2*smutom2*cosmutI*cmutom1-2.D0*e1**2*smut&
     &om1*q2*cmutom2**2*e2**2*smutom2*cosmutI*cmutom1+2.D0*e1*smutom1*cm&
     &utom2**2*q1*e2**2*smutom2*cosmutI-q2*sinmutI**2*e1**2*cmutom1**2+2&
     &.D0*e1**2*smutom1*cmutom2**2*q1*e2**2*smutom2*cosmutI+4.D0*e1**2*s&
     &mutom1*cmutom2*q1*e2*smutom2*cosmutI+2.D0*e1**2*smutom1*q1*smutom2&
     &*cosmutI-2.D0*q1*e2**2*smutom2**2*cmutom2+2.D0*e1*smutom1*q1*smuto&
     &m2*cosmutI-2.D0*e1**2*smutom1*q2*smutom2*cosmutI*cmutom1-2.D0*q1*e&
     &1**2*e2**2*smutom2**2*cmutom1*cmutom2-q2*sinmutI**2*e2**2*cmutom2*&
     &*2                                                                
      s7 = -2.D0*q2*sinmutI**2*e1*cmutom1+2.D0*q2*e1**2*cmutom1**2*e2*sm&
     &utom2**2-2.D0*q2*sinmutI**2*e2*cmutom2+4.D0*q2*e2*smutom2**2*e1*cm&
     &utom1-2.D0*q1*e1**2*e2*smutom2**2*cmutom1+4.D0*e1*smutom1*cmutom2*&
     &q1*e2*smutom2*cosmutI+q2*cmutom2*e2**2*smutom2**2*e1**2*cmutom1**2&
     &+2.D0*q2*cmutom2*sinmutI**2*e1*cmutom1+q2*cmutom2**3*sinmutI**2*e2&
     &**2+q2*cmutom2*sinmutI**2*e1**2*cmutom1**2-4.D0*q2*sinmutI**2*e1*c&
     &mutom1*e2*cmutom2+2.D0*q2*cmutom2**2*sinmutI**2*e2-2.D0*e1*smutom1&
     &*q2*smutom2*cosmutI+q2*e2**2*smutom2**2*e1**2*cmutom1**2          
      s5 = s7-2.D0*q1*e1*e2*smutom2**2+2.D0*q2*e2**2*smutom2**2*e1*cmuto&
     &m1-q2*sinmutI**2*e1**2*cmutom1**2*e2**2*cmutom2**2-2.D0*q2*sinmutI&
     &**2*e1*cmutom1*e2**2*cmutom2**2-2.D0*q2*sinmutI**2*e1**2*cmutom1**&
     &2*e2*cmutom2+2.D0*q2*cmutom2*e2**2*smutom2**2*e1*cmutom1-2.D0*q1*e&
     &2*smutom2**2*e1*cmutom1+2.D0*q2*cmutom2**2*sinmutI**2*e1**2*cmutom&
     &1**2*e2+2.D0*q2*cmutom2**3*sinmutI**2*e1*cmutom1*e2**2+4.D0*q2*cmu&
     &tom2**2*sinmutI**2*e1*cmutom1*e2-2.D0*q1*e1*e2**2*smutom2**2*cmuto&
     &m2-2.D0*q1*e2**2*smutom2**2*e1*cmutom1*cmutom2+q2*cmutom2**3*sinmu&
     &tI**2*e1**2*cmutom1**2*e2**2+s6                                   
      s6 = 1.D0/((e2**2*smutom2**2+2.D0*e2**2*smutom2**2*e1*cmutom1+e2**&
     &2*smutom2**2*e1**2*cmutom1**2+e1**2*smutom1**2+2.D0*e1**2*smutom1*&
     &*2*e2*cmutom2+e1**2*smutom1**2*e2**2*cmutom2**2+sinmutI**2+2.D0*si&
     &nmutI**2*e2*cmutom2+sinmutI**2*e2**2*cmutom2**2+2.D0*sinmutI**2*e1&
     &*cmutom1+4.D0*sinmutI**2*e1*cmutom1*e2*cmutom2+2.D0*sinmutI**2*e1*&
     &cmutom1*e2**2*cmutom2**2+sinmutI**2*e1**2*cmutom1**2+2.D0*sinmutI*&
     &*2*e1**2*cmutom1**2*e2*cmutom2+sinmutI**2*e1**2*cmutom1**2*e2**2*c&
     &mutom2**2-2.D0*e1*smutom1*e2*smutom2*cosmutI-2.D0*e1*smutom1*e2**2&
     &*smutom2*cosmutI*cmutom2-2.D0*e1**2*smutom1*e2*smutom2*cosmutI*cmu&
     &tom1-2.D0*e1**2*smutom1*e2**2*smutom2*cosmutI*cmutom1*cmutom2)**2.&
     &D0)                                                               
      s4 = s5*s6 
      s2 = s3*s4 
      ddeltayde2 = s1*s2 
!
    ddeltaydi1=k1p*dgi1pdcmutompl*dcosom(1,1,1)-k2p*(dcosmutI(1,1)*gi2p+ &
         & cosmutI*dgi2pdcmutom*dcosom(2,1,1))
    ddeltaydi2=k1p*dgi1pdcmutompl*dcosom(1,1,2)-k2p*(dcosmutI(1,2)*gi2p+ &
         & cosmutI*dgi2pdcmutom*dcosom(2,1,2))
!
    ddeltaydbigom1=k1p*dgi1pdcmutompl*dcosom(1,3,1)-k2p* &
         & (dcosmutI(3,1)*gi2p+cosmutI*dgi2pdcmutom*dcosom(2,3,1))
    ddeltaydbigom2=k1p*dgi1pdcmutompl*dcosom(1,3,2)-k2p* &
         & (dcosmutI(3,2)*gi2p+cosmutI*dgi2pdcmutom*dcosom(2,3,2))
!
    ddeltaydom1=k1p*dgi1pdcmutompl*dcosom(1,2,1)-k2p* &
         & (dcosmutI(2,1)*gi2p+cosmutI*dgi2pdcmutom*dcosom(2,2,1))
    ddeltaydom2=k1p*dgi1pdcmutompl*dcosom(1,2,2)-k2p* &
         & (dcosmutI(2,2)*gi2p+cosmutI*dgi2pdcmutom*dcosom(2,2,2))
!------------------------------------------------------------------
! Derivatives of Delta z
!------------------------------------------------------------------
      t0 = (1.D0+e2*cmutom2)*sinmutI*(1.D0+e1)*(e2*smutom2+e2*smutom2*e1&
     &*cmutom1-e1*smutom1*cosmutI-e1*smutom1*cosmutI*e2*cmutom2)/(e2**2*&
     &smutom2**2+2.D0*e2**2*smutom2**2*e1*cmutom1+e2**2*smutom2**2*e1**2&
     &*cmutom1**2+e1**2*smutom1**2+2.D0*e1**2*smutom1**2*e2*cmutom2+e1**&
     &2*smutom1**2*e2**2*cmutom2**2+sinmutI**2+2.D0*sinmutI**2*e2*cmutom&
     &2+sinmutI**2*e2**2*cmutom2**2+2.D0*sinmutI**2*e1*cmutom1+4.D0*sinm&
     &utI**2*e1*cmutom1*e2*cmutom2+2.D0*sinmutI**2*e1*cmutom1*e2**2*cmut&
     &om2**2+sinmutI**2*e1**2*cmutom1**2+2.D0*sinmutI**2*e1**2*cmutom1**&
     &2*e2*cmutom2+sinmutI**2*e1**2*cmutom1**2*e2**2*cmutom2**2-2.D0*e1*&
     &smutom1*e2*smutom2*cosmutI-2.D0*e1*smutom1*e2**2*smutom2*cosmutI*c&
     &mutom2-2.D0*e1**2*smutom1*e2*smutom2*cosmutI*cmutom1-2.D0*e1**2*sm&
     &utom1*e2**2*smutom2*cosmutI*cmutom1*cmutom2)                      
      ddeltazdq1=t0 
!                                                                        
       t0 = -sinmutI*(1.D0+e2)*(1.D0+e1*cmutom1)*(e2*smutom2+e2*smutom2*e&
     &1*cmutom1-e1*smutom1*cosmutI-e1*smutom1*cosmutI*e2*cmutom2)/(e2**2&
     &*smutom2**2+2.D0*e2**2*smutom2**2*e1*cmutom1+e2**2*smutom2**2*e1**&
     &2*cmutom1**2+e1**2*smutom1**2+2.D0*e1**2*smutom1**2*e2*cmutom2+e1*&
     &*2*smutom1**2*e2**2*cmutom2**2+sinmutI**2+2.D0*sinmutI**2*e2*cmuto&
     &m2+sinmutI**2*e2**2*cmutom2**2+2.D0*sinmutI**2*e1*cmutom1+4.D0*sin&
     &mutI**2*e1*cmutom1*e2*cmutom2+2.D0*sinmutI**2*e1*cmutom1*e2**2*cmu&
     &tom2**2+sinmutI**2*e1**2*cmutom1**2+2.D0*sinmutI**2*e1**2*cmutom1*&
     &*2*e2*cmutom2+sinmutI**2*e1**2*cmutom1**2*e2**2*cmutom2**2-2.D0*e1&
     &*smutom1*e2*smutom2*cosmutI-2.D0*e1*smutom1*e2**2*smutom2*cosmutI*&
     &cmutom2-2.D0*e1**2*smutom1*e2*smutom2*cosmutI*cmutom1-2.D0*e1**2*s&
     &mutom1*e2**2*smutom2*cosmutI*cmutom1*cmutom2)                     
      ddeltazdq2=t0 
!
      s1 = -sinmutI 
      s3 = 1.D0+e2*cmutom2 
      s8 = 2.D0*q1*e2**3*smutom2**2*e1**2*smutom1*cosmutI*cmutom1*cmutom&
     &2+2.D0*q1*e2**2*smutom2**2*e1**2*smutom1*cosmutI*cmutom1+q1*e2*smu&
     &tom2*cmutom1*e1**2*smutom1**2+2.D0*q1*e2**3*smutom2**2*e1*smutom1*&
     &cosmutI*cmutom2+2.D0*q1*cmutom1*e1**2*smutom1**2*cosmutI**2*e2*smu&
     &tom2+2.D0*q1*e1**2*smutom1*cosmutI*e2**3*cmutom2**3*sinmutI**2*cmu&
     &tom1-2.D0*q1*e1**2*smutom1**2*cosmutI**2*e2**3*cmutom2**2*smutom2-&
     &2.D0*q2*e2**2*smutom2*e1*smutom1**2*cmutom2-q1*e2*smutom2*sinmutI*&
     &*2-2.D0*e2**2*q2*smutom1*cosmutI*sinmutI**2*cmutom2+2.D0*q1*e2**2*&
     &smutom2**2*e1*smutom1*cosmutI+e2**3*q2*smutom1**3*cosmutI*e1**2*cm&
     &utom2**2-2.D0*q1*e1**2*smutom1**2*cosmutI**2*e2*smutom2           
      s7 = s8+2.D0*q1*e1*smutom1*cosmutI*e2**3*cmutom2**3*sinmutI**2-q1*&
     &e2*smutom2*sinmutI**2*e1**2*cmutom1**2+2.D0*q1*e2**2*smutom2*cmuto&
     &m1*sinmutI**2*cmutom2+3.D0*q1*smutom1*cosmutI*e2**2*cmutom2**2*sin&
     &mutI**2+3.D0*q1*smutom1*cosmutI*e2*cmutom2*sinmutI**2+q1*smutom1*c&
     &osmutI*sinmutI**2-q2*smutom1*cosmutI*sinmutI**2-q1*smutom1**3*cosm&
     &utI*e1**2+q2*smutom1**3*cosmutI*e1**2+2.D0*q1*e2**2*smutom2*e1**2*&
     &smutom1**2*cmutom2+q1*smutom1*cosmutI*e2**3*cmutom2**3*sinmutI**2+&
     &2.D0*q1*e1**2*smutom1*cosmutI*sinmutI**2*cmutom1-q1*smutom1*cosmut&
     &I*e2**3*cmutom2*smutom2**2                                        
      s8 = s7-q1*e2**3*smutom2**3+4.D0*q1*e2**2*smutom2*e1*smutom1**2*cm&
     &utom2-2.D0*e2**2*q2*smutom2*e1*smutom1**2-2.D0*cmutom1*q2*e1*smuto&
     &m1*cosmutI*sinmutI**2*e2**2*cmutom2**2-2.D0*e2**3*q2*smutom2*e1*sm&
     &utom1**2*cmutom2+2.D0*cmutom1*q2*e1*smutom1*cosmutI*e2**2*smutom2*&
     &*2-2.D0*cmutom1*e2**3*q2*smutom2*e1**2*smutom1**2*cmutom2+2.D0*cmu&
     &tom1*e2**3*q2*smutom2**2*e1*smutom1*cosmutI-2.D0*cmutom1*e2**2*q2*&
     &smutom2*e1**2*smutom1**2-4.D0*cmutom1*q2*e1*smutom1*cosmutI*sinmut&
     &I**2*e2*cmutom2-2.D0*cmutom1*q2*e1*smutom1*cosmutI*sinmutI**2-2.D0&
     &*cmutom1*e2**3*q2*e1*smutom1*cosmutI*sinmutI**2*cmutom2**2-2.D0*cm&
     &utom1*e2*q2*e1*smutom1*cosmutI*sinmutI**2                         
      s6 = s8-2.D0*cmutom1*q2*e2**2*smutom2*e1**2*smutom1**2*cmutom2-4.D&
     &0*cmutom1*e2**2*q2*e1*smutom1*cosmutI*sinmutI**2*cmutom2-2.D0*q2*e&
     &2*smutom2*e1*smutom1**2+e2**3*q2*smutom2**2*e1**2*smutom1*cosmutI*&
     &cmutom1**2-2.D0*cmutom1*q2*e2*smutom2*e1**2*smutom1**2-q2*e1**2*sm&
     &utom1*cosmutI*sinmutI**2*cmutom1**2*e2**2*cmutom2**2-2.D0*q2*e1**2&
     &*smutom1*cosmutI*sinmutI**2*cmutom1**2*e2*cmutom2+q2*e1**2*smutom1&
     &*cosmutI*e2**2*smutom2**2*cmutom1**2-2.D0*cmutom1**2*q2*e2**2*e1**&
     &2*smutom1*cosmutI*sinmutI**2*cmutom2-q2*e1**2*smutom1*cosmutI*sinm&
     &utI**2*cmutom1**2-cmutom1**2*q2*e2*e1**2*smutom1*cosmutI*sinmutI**&
     &2-cmutom1**2*q2*e2**3*e1**2*smutom1*cosmutI*sinmutI**2*cmutom2**2-&
     &2.D0*q1*e2**3*smutom2**3*e1*cmutom1-q1*e2**3*smutom2**3*e1**2*cmut&
     &om1**2                                                            
      s8 = s6+q1*e2*smutom2*e1**2*smutom1**2+6.D0*q1*e1**2*smutom1*cosmu&
     &tI*sinmutI**2*cmutom1*e2*cmutom2+2.D0*q1*cmutom1*e1**2*smutom1**2*&
     &cosmutI**2*e2**3*cmutom2**2*smutom2+6.D0*q1*e1*smutom1*cosmutI*sin&
     &mutI**2*e2*cmutom2-3.D0*q1*e2**3*smutom2**2*cmutom1**2*e1**2*smuto&
     &m1*cosmutI*cmutom2-q1*cmutom1**2*e1**2*smutom1*cosmutI*e2**3*cmuto&
     &m2**3*sinmutI**2+q1*e2**3*smutom2**3*cmutom1+6.D0*q1*e1*smutom1*co&
     &smutI*sinmutI**2*e2**2*cmutom2**2+6.D0*q1*e1**2*smutom1*cosmutI*si&
     &nmutI**2*cmutom1*e2**2*cmutom2**2+q1*e2**3*smutom2*cmutom1*sinmutI&
     &**2*cmutom2**2+2.D0*q1*e2**2*smutom2*cmutom1*e1**2*smutom1**2*cmut&
     &om2+q1*e2**3*smutom2*cmutom1*e1**2*smutom1**2*cmutom2**2          
      s7 = s8-q1*cmutom1**2*e1**2*smutom1*cosmutI*sinmutI**2-3.D0*q1*cmu&
     &tom1**2*e1**2*smutom1*cosmutI*sinmutI**2*e2*cmutom2+q2*smutom1*cos&
     &mutI*e2**2*smutom2**2+e2**3*q2*smutom1*cosmutI*smutom2**2+e2*q2*sm&
     &utom1**3*cosmutI*e1**2-q1*smutom1*cosmutI*e2**2*smutom2**2-3.D0*q1&
     &*cmutom1**2*e1**2*smutom1*cosmutI*sinmutI**2*e2**2*cmutom2**2+4.D0&
     &*q1*cmutom1*e1**2*smutom1**2*cosmutI**2*e2**2*smutom2*cmutom2+q2*s&
     &mutom1**3*cosmutI*e2**2*cmutom2**2*e1**2+2.D0*q2*smutom1**3*cosmut&
     &I*e2*cmutom2*e1**2-2.D0*q2*smutom1*cosmutI*e2*cmutom2*sinmutI**2-q&
     &2*smutom1*cosmutI*e2**2*cmutom2**2*sinmutI**2-3.D0*q1*smutom1**3*c&
     &osmutI*e2**2*cmutom2**2*e1**2-q1*smutom1**3*cosmutI*e2**3*cmutom2*&
     &*3*e1**2                                                          
      s8 = s7-3.D0*q1*smutom1**3*cosmutI*e2*cmutom2*e1**2+2.D0*q1*e2*smu&
     &tom2*e1*smutom1**2+2.D0*q1*e2**3*smutom2*e1*smutom1**2*cmutom2**2-&
     &3.D0*q1*e2**2*smutom2**2*cmutom1**2*e1**2*smutom1*cosmutI-4.D0*q1*&
     &e2**2*smutom2**2*cmutom1*e1*smutom1*cosmutI+q1*e2*smutom2*cmutom1*&
     &*3*sinmutI**2*e1**2+2.D0*q1*e2**2*smutom2*cmutom1**3*sinmutI**2*e1&
     &**2*cmutom2+q1*e2**3*smutom2*cmutom1**3*sinmutI**2*e1**2*cmutom2**&
     &2+4.D0*q1*e2**2*smutom2*cmutom1**2*sinmutI**2*e1*cmutom2+2.D0*q1*e&
     &2**3*smutom2*cmutom1**2*sinmutI**2*e1*cmutom2**2+2.D0*q1*e2*smutom&
     &2*cmutom1**2*sinmutI**2*e1-2.D0*q1*e2*smutom2*sinmutI**2*e1*cmutom&
     &1-4.D0*q1*e2**2*smutom2*sinmutI**2*e1*cmutom1*cmutom2             
      s5 = s8-2.D0*q1*e2**3*smutom2*sinmutI**2*e1*cmutom1*cmutom2**2+q1*&
     &e2**3*smutom2*e1**2*smutom1**2*cmutom2**2-4.D0*q1*e2**3*smutom2**2&
     &*cmutom1*e1*smutom1*cosmutI*cmutom2+2.D0*e2**2*q2*smutom1**3*cosmu&
     &tI*e1**2*cmutom2-e2*q2*smutom1*cosmutI*sinmutI**2-e2**3*q2*smutom1&
     &*cosmutI*sinmutI**2*cmutom2**2-4.D0*q1*e1**2*smutom1**2*cosmutI**2&
     &*e2**2*smutom2*cmutom2-2.D0*q1*e2**2*smutom2*sinmutI**2*e1**2*cmut&
     &om1**2*cmutom2+q1*e2*smutom2*cmutom1*sinmutI**2+2.D0*q1*e1*smutom1&
     &*cosmutI*sinmutI**2+2.D0*q1*e2**3*smutom2**3*cmutom1**2*e1+q1*e2**&
     &3*smutom2**3*cmutom1**3*e1**2-2.D0*q1*e2**2*smutom2*sinmutI**2*cmu&
     &tom2-q1*e2**3*smutom2*sinmutI**2*cmutom2**2-q1*e2**3*smutom2*sinmu&
     &tI**2*e1**2*cmutom1**2*cmutom2**2                                 
      s6 = 1.D0/((e2**2*smutom2**2+2.D0*e2**2*smutom2**2*e1*cmutom1+e2**&
     &2*smutom2**2*e1**2*cmutom1**2+e1**2*smutom1**2+2.D0*e1**2*smutom1*&
     &*2*e2*cmutom2+e1**2*smutom1**2*e2**2*cmutom2**2+sinmutI**2+2.D0*si&
     &nmutI**2*e2*cmutom2+sinmutI**2*e2**2*cmutom2**2+2.D0*sinmutI**2*e1&
     &*cmutom1+4.D0*sinmutI**2*e1*cmutom1*e2*cmutom2+2.D0*sinmutI**2*e1*&
     &cmutom1*e2**2*cmutom2**2+sinmutI**2*e1**2*cmutom1**2+2.D0*sinmutI*&
     &*2*e1**2*cmutom1**2*e2*cmutom2+sinmutI**2*e1**2*cmutom1**2*e2**2*c&
     &mutom2**2-2.D0*e1*smutom1*e2*smutom2*cosmutI-2.D0*e1*smutom1*e2**2&
     &*smutom2*cosmutI*cmutom2-2.D0*e1**2*smutom1*e2*smutom2*cosmutI*cmu&
     &tom1-2.D0*e1**2*smutom1*e2**2*smutom2*cosmutI*cmutom1*cmutom2)**2.&
     &D0)                                                               
      s4 = s5*s6 
      s2 = s3*s4 
      ddeltazde1 = s1*s2 
!
     s2 = 1.D0+e1*cmutom1 
      s6 = q1*smutom2*sinmutI**2-q1*smutom2**3*e2**2+q2*smutom2**3*e2**2&
     &+q1*e1**3*smutom2*smutom1**2+q1*e1*smutom2*sinmutI**2-q2*smutom2*e&
     &1**2*smutom1**2-q2*smutom2*sinmutI**2+q1*smutom2*e1**2*smutom1**2-&
     &q1*e1*smutom2**3*e2**2-2.D0*q2*e2*smutom2*sinmutI**2+2.D0*cmutom2*&
     &q1*e1**2*smutom1*cosmutI*e2**2*smutom2**2+2.D0*q2*e1**2*cmutom1*sm&
     &utom1*cosmutI*sinmutI**2-2.D0*q1*smutom2**3*e2**2*e1*cmutom1      
      s5 = s6-q1*smutom2**3*e2**2*e1**2*cmutom1**2-2.D0*q1*e1**2*smutom2&
     &**3*e2**2*cmutom1-q1*e1**3*smutom2**3*e2**2*cmutom1**2-cmutom2**3*&
     &q2*e1*smutom1*cosmutI*sinmutI**2*e2**2+2.D0*q1*smutom2**2*e1*smuto&
     &m1*e2*cosmutI-4.D0*cmutom2*q1*e1**2*smutom1**2*cosmutI**2*e2*smuto&
     &m2+3.D0*cmutom2**2*q2*e2**2*smutom2*sinmutI**2*e1**2*cmutom1**2+3.&
     &D0*q2*smutom2**3*e2**2*e1**2*cmutom1**2+3.D0*q2*smutom2**3*e2**2*e&
     &1*cmutom1+2.D0*cmutom2*q1*e2*smutom2*sinmutI**2+2.D0*cmutom2*q1*e1&
     &*e2*smutom2*sinmutI**2-2.D0*q2*smutom2**2*e1*smutom1*e2*cosmutI+q2&
     &*e1**3*cmutom1**3*smutom2**3*e2**2+2.D0*cmutom2*q1*e2*smutom2*e1**&
     &2*cmutom1**2*sinmutI**2                                           
      s6 = s5-4.D0*q2*smutom2**2*e1**2*smutom1*e2*cosmutI*cmutom1+2.D0*q&
     &1*smutom2**2*e1**2*smutom1*e2*cosmutI*cmutom1+2.D0*cmutom2**2*q2*e&
     &1**3*smutom1**2*cosmutI**2*e2**2*smutom2*cmutom1-4.D0*cmutom2**2*q&
     &2*e1**2*smutom1*cosmutI*sinmutI**2*cmutom1*e2+2.D0*cmutom2*q1*e1**&
     &3*e2*smutom2*smutom1**2+4.D0*cmutom2*q2*e1**3*smutom1**2*cosmutI**&
     &2*e2*smutom2*cmutom1+cmutom2**2*q1*e1**2*cmutom1**2*e2**2*smutom2*&
     &sinmutI**2+q2*e1**3*smutom1**3*cosmutI-2.D0*cmutom2**2*q1*e1**2*sm&
     &utom1**2*cosmutI**2*e2**2*smutom2-2.D0*cmutom2**3*q2*e1**2*smutom1&
     &*cosmutI*sinmutI**2*cmutom1*e2**2+cmutom2**2*q1*e1**3*e2**2*smutom&
     &2*smutom1**2-cmutom2**3*q2*e1**3*smutom1*cosmutI*sinmutI**2*cmutom&
     &1**2*e2**2-2.D0*cmutom2**2*q2*e1**3*smutom1*cosmutI*sinmutI**2*cmu&
     &tom1**2*e2                                                        
      s4 = s6+cmutom2**2*q1*e1*e2**2*smutom2*sinmutI**2-3.D0*cmutom2*q2*&
     &e1**3*smutom1*cosmutI*e2**2*smutom2**2*cmutom1**2-2.D0*cmutom2**2*&
     &q2*e1**3*smutom1**3*cosmutI*e2+q2*e1**3*cmutom1**2*smutom1*cosmutI&
     &*sinmutI**2+q1*smutom2*sinmutI**2*e1**2*cmutom1**2-cmutom2**3*q2*e&
     &1**3*smutom1**3*cosmutI*e2**2+2.D0*cmutom2*q1*e1**3*smutom1*cosmut&
     &I*e2**2*smutom2**2*cmutom1+2.D0*cmutom2*q1*e1*smutom1*cosmutI*e2**&
     &2*smutom2**2+4.D0*cmutom2*q1*e2*smutom2*e1*cmutom1*sinmutI**2+2.D0&
     &*cmutom2**2*q1*e2**2*smutom2*e1*cmutom1*sinmutI**2+2.D0*cmutom2*q1&
     &*e1**3*e2*smutom2*cmutom1**2*sinmutI**2+cmutom2**2*q2*e2**2*smutom&
     &2*e1**3*cmutom1**3*sinmutI**2-cmutom2*q2*e1*smutom1*cosmutI*sinmut&
     &I**2-cmutom2*q2*e1**3*smutom1**3*cosmutI+cmutom2**2*q2*e2**2*smuto&
     &m2*e1**2*smutom1**2                                               
      s6 = s4-6.D0*cmutom2*q2*e1**2*smutom1*cosmutI*e2**2*smutom2**2*cmu&
     &tom1+4.D0*q2*e1**2*smutom1*cosmutI*sinmutI**2*cmutom1*e2*cmutom2+q&
     &2*e1**3*smutom1*cosmutI*sinmutI**2*cmutom1**2*e2**2*cmutom2**2+2.D&
     &0*q2*e1**2*smutom1*cosmutI*sinmutI**2*cmutom1*e2**2*cmutom2**2-6.D&
     &0*q2*e2**2*smutom2*sinmutI**2*e1*cmutom1*cmutom2-6.D0*q2*e2*smutom&
     &2*sinmutI**2*e1*cmutom1-2.D0*q2*e2**2*smutom2*e1**3*cmutom1**3*sin&
     &mutI**2*cmutom2+2.D0*cmutom2**2*q1*e1**2*e2**2*smutom2*cmutom1*sin&
     &mutI**2-2.D0*q2*e2*smutom2*e1**3*cmutom1*smutom1**2+2.D0*q2*e1**3*&
     &smutom1*cosmutI*sinmutI**2*cmutom1**2*e2*cmutom2-2.D0*q2*e2*smutom&
     &2*e1**3*cmutom1**3*sinmutI**2-2.D0*q2*e2**2*smutom2*e1**3*cmutom1*&
     &smutom1**2*cmutom2-2.D0*q2*e2*smutom2*e1**2*smutom1**2            
      s5 = s6-2.D0*q2*e2**2*smutom2*sinmutI**2*cmutom2+q2*e1**3*smutom1*&
     &cosmutI*e2**2*smutom2**2*cmutom1**2+q2*e1*smutom1*cosmutI*e2**2*sm&
     &utom2**2+q2*e1**3*smutom1**3*cosmutI*e2**2*cmutom2**2+2.D0*q2*e1**&
     &3*smutom1**3*cosmutI*e2*cmutom2+q2*e1*smutom1*cosmutI*sinmutI**2*e&
     &2**2*cmutom2**2-2.D0*q2*e2**2*smutom2*e1**2*smutom1**2*cmutom2-6.D&
     &0*q2*e2*smutom2*sinmutI**2*e1**2*cmutom1**2+cmutom2**2*q1*e2**2*sm&
     &utom2*sinmutI**2+4.D0*cmutom2*q1*e1**2*e2*smutom2*cmutom1*sinmutI*&
     &*2-6.D0*q2*e2**2*smutom2*sinmutI**2*e1**2*cmutom1**2*cmutom2+2.D0*&
     &cmutom2*q1*e2*smutom2*e1**2*smutom1**2+cmutom2**2*q2*e2**2*smutom2&
     &*sinmutI**2+cmutom2**2*q1*e1**3*e2**2*smutom2*sinmutI**2*cmutom1**&
     &2                                                                 
      s6 = s5+cmutom2**2*q1*e2**2*smutom2*e1**2*smutom1**2+2.D0*q1*e1**2&
     &*smutom2**2*smutom1*e2*cosmutI+2.D0*cmutom2*q1*e1**2*smutom1*cosmu&
     &tI*e2**2*smutom2**2*cmutom1-4.D0*cmutom2*q1*e1**3*smutom1**2*cosmu&
     &tI**2*e2*smutom2-2.D0*q2*e1**3*cmutom1**2*smutom2**2*smutom1*e2*co&
     &smutI-2.D0*cmutom2**2*q1*e1**3*smutom1**2*cosmutI**2*e2**2*smutom2&
     &+3.D0*cmutom2**2*q2*e2**2*smutom2*sinmutI**2*e1*cmutom1+2.D0*q1*e1&
     &**3*smutom2**2*smutom1*e2*cosmutI*cmutom1+4.D0*cmutom2*q2*e1**2*sm&
     &utom1**2*cosmutI**2*e2*smutom2-cmutom2*q2*e1**3*smutom1*cosmutI*si&
     &nmutI**2*cmutom1**2+2.D0*cmutom2**2*q2*e1**2*smutom1**2*cosmutI**2&
     &*e2**2*smutom2-2.D0*cmutom2*q2*e1**2*smutom1*cosmutI*sinmutI**2*cm&
     &utom1+q2*e1*smutom1*cosmutI*sinmutI**2+q1*e1**3*smutom2*sinmutI**2&
     &*cmutom1**2                                                       
      s3 = s6-3.D0*cmutom2*q2*e1*smutom1*cosmutI*e2**2*smutom2**2-2.D0*c&
     &mutom2**2*q2*e1*smutom1*cosmutI*sinmutI**2*e2+cmutom2**2*q2*e2**2*&
     &smutom2*e1**3*cmutom1*smutom1**2+2.D0*q2*e1**2*smutom1*cosmutI*e2*&
     &*2*smutom2**2*cmutom1+2.D0*q2*e1*smutom1*cosmutI*sinmutI**2*e2*cmu&
     &tom2-q2*smutom2*e1**3*cmutom1*smutom1**2+2.D0*q2*e1**2*smutom1**2*&
     &cosmutI**2*smutom2+2.D0*q1*smutom2*e1*cmutom1*sinmutI**2+2.D0*q1*e&
     &1**2*smutom2*cmutom1*sinmutI**2-2.D0*q1*e1**2*smutom1**2*cosmutI**&
     &2*smutom2-2.D0*q1*e1**3*smutom1**2*cosmutI**2*smutom2-3.D0*q2*smut&
     &om2*sinmutI**2*e1*cmutom1-3.D0*q2*smutom2*sinmutI**2*e1**2*cmutom1&
     &**2-q2*smutom2*e1**3*cmutom1**3*sinmutI**2+2.D0*q2*e1**3*smutom1**&
     &2*cosmutI**2*smutom2*cmutom1                                      
      s1 = s2*s3 
      s2 = sinmutI/(e2**2*smutom2**2+2.D0*e2**2*smutom2**2*e1*cmutom1+e2&
     &**2*smutom2**2*e1**2*cmutom1**2+e1**2*smutom1**2+2.D0*e1**2*smutom&
     &1**2*e2*cmutom2+e1**2*smutom1**2*e2**2*cmutom2**2+sinmutI**2+2.D0*&
     &sinmutI**2*e2*cmutom2+sinmutI**2*e2**2*cmutom2**2+2.D0*sinmutI**2*&
     &e1*cmutom1+4.D0*sinmutI**2*e1*cmutom1*e2*cmutom2+2.D0*sinmutI**2*e&
     &1*cmutom1*e2**2*cmutom2**2+sinmutI**2*e1**2*cmutom1**2+2.D0*sinmut&
     &I**2*e1**2*cmutom1**2*e2*cmutom2+sinmutI**2*e1**2*cmutom1**2*e2**2&
     &*cmutom2**2-2.D0*e1*smutom1*e2*smutom2*cosmutI-2.D0*e1*smutom1*e2*&
     &*2*smutom2*cosmutI*cmutom2-2.D0*e1**2*smutom1*e2*smutom2*cosmutI*c&
     &mutom1-2.D0*e1**2*smutom1*e2**2*smutom2*cosmutI*cmutom1*cmutom2)**&
     &2.D0                                                              
      ddeltazde2 = s1*s2 
!
    ddeltazdi1=-k2p*(dgi2pdcmutom*dcosom(2,1,1)*sinmutI+ &
         & gi2p*dsinmutIdcos*dcosmutI(1,1))
    ddeltazdi2=-k2p*(dgi2pdcmutom*dcosom(2,1,2)*sinmutI+ &
         & gi2p*dsinmutIdcos*dcosmutI(1,2))
  
    ddeltazdbigom1=-k2p*(dgi2pdcmutom*dcosom(2,3,1)*sinmutI+ &
         & gi2p*dsinmutIdcos*dcosmutI(3,1))
    ddeltazdbigom2=-k2p*(dgi2pdcmutom*dcosom(2,3,2)*sinmutI+ &
         & gi2p*dsinmutIdcos*dcosmutI(3,2))

    ddeltazdom1=-k2p*(dgi2pdcmutom*dcosom(2,2,1)*sinmutI+ &
         & gi2p*dsinmutIdcos*dcosmutI(2,1))
    ddeltazdom2=-k2p*(dgi2pdcmutom*dcosom(2,2,2)*sinmutI+ &
         & gi2p*dsinmutIdcos*dcosmutI(2,2))
! -----------------------------------------------------------
! Derivatives of the difference of the Cartesian positions
! w.r.t. the 10 elements
! -----------------------------------------------------------
    ddelta(1,1)=ddeltaxdq1
    ddelta(1,2)=ddeltaydq1
    ddelta(1,3)=ddeltazdq1
    ddelta(2,1)=ddeltaxde1
    ddelta(2,2)=ddeltayde1
    ddelta(2,3)=ddeltazde1
    ddelta(3,1)=ddeltaxdi1
    ddelta(3,2)=ddeltaydi1
    ddelta(3,3)=ddeltazdi1
    ddelta(4,1)=ddeltaxdbigom1
    ddelta(4,2)=ddeltaydbigom1
    ddelta(4,3)=ddeltazdbigom1
    ddelta(5,1)=ddeltaxdom1
    ddelta(5,2)=ddeltaydom1
    ddelta(5,3)=ddeltazdom1
    ddelta(6,1)=ddeltaxdq2
    ddelta(6,2)=ddeltaydq2
    ddelta(6,3)=ddeltazdq2
    ddelta(7,1)=ddeltaxde2
    ddelta(7,2)=ddeltayde2
    ddelta(7,3)=ddeltazde2
    ddelta(8,1)=ddeltaxdi2
    ddelta(8,2)=ddeltaydi2
    ddelta(8,3)=ddeltazdi2 
    ddelta(9,1)=ddeltaxdbigom2
    ddelta(9,2)=ddeltaydbigom2
    ddelta(9,3)=ddeltazdbigom2
    ddelta(10,1)=ddeltaxdom2
    ddelta(10,2)=ddeltaydom2
    ddelta(10,3)=ddeltazdom2
! Computation of tau3p in the inertial reference frame
    DO j=1,3   
       newtau3p(j)=DOT_PRODUCT(rotinmut(j,1:3),tau3phat)
    END DO
! Computation of the derivatives of the AMOIDp
    DO h=1,10       
       DO k=1,3 
          ddeltadkep(k)=ddelta(h,k)
          deltarot(k)=DOT_PRODUCT(drotinmut(h,k,1:3),deltaxyzp)
       ENDDO
       deramoidp(1,h)=DOT_PRODUCT(tau3phat,ddeltadkep)+ &
       & DOT_PRODUCT(newtau3p,deltarot)
    ENDDO

!    IF(chk_der)THEN
!       deramoidpdequ2(1,1:5)= &
!            & MATMUL(deramoidp(1,6:10),jackepeq2(1:5,1:5)) 
!       write(*,113)'derivatives of AMOIDP:',deramoidpdequ2(1,1:5)
!113    FORMAT(A22,2x,5(es16.9,1x)) 
!    ENDIF

    tderamoidp=TRANSPOSE(deramoidp)
! Computation of the RMS of the AMOIDp
    amoidprmsvec=MATMUL(deramoidp,MATMUL(covcom,tderamoidp)) 
    amoidprms=SQRT(amoidprmsvec(1,1))
! -----------------------------------------
! Approximate MOID at the descending node
! -----------------------------------------
    effe1m=effe1p
    gi1m=-q1*(1.d0-e1*cmutom1)/(beta1*(1-e1))
    effe2m=effe2p
    gi2m=-q2*(1.d0-e2*cmutom2)/(beta2*(1-e2))
    detAm=effe2m**2*gi1m**2+effe1m**2*gi2m**2+gi1m**2*gi2m**2*sinmutI**2 &
         & -2.d0*effe1m*effe2m*gi1m*gi2m*cosmutI
    amoidm=ABS(dnodm)*SQRT(1.d0-(effe1m**2*gi2m**2-2.d0*effe1m*effe2m* &
         & gi1m*gi2m*cosmutI+effe2m**2*gi1m**2)/detAm)
! Computation of vectors tau1,tau2,tau3
    tau1m(1)=-effe1m
    tau1m(2)=gi1m
    tau1m(3)=0.d0
    tau2m(1)=-effe2m
    tau2m(2)=gi2m*cosmutI
    tau2m(3)=gi2m*sinmutI
    CALL prvec(tau1m,tau2m,tau3m)
    modtau3m=vsize(tau3m)
    tau3mhat(1:3)=tau3m(1:3)/modtau3m
! Computation of Delta x, Delta y, Delta z  
    x1barm=-q1*(1+e1)/(1.d0-e1*cmutom1)
    x2barm=-q2*(1+e2)/(1.d0-e2*cmutom2)
    k1m=dnodm*gi2m*(effe1m*gi2m-effe2m*gi1m*cosmutI)/detAm
    k2m=dnodm*gi1m*(effe1m*gi2m*cosmutI-effe2m*gi1m)/detAm
    xyz1m(1)=x1barm-effe1m*k1m
    xyz1m(2)=gi1m*k1m
    xyz1m(3)=0.d0
    xyz2m(1)=x2barm-effe2m*k2m
    xyz2m(2)=gi2m*cosmutI*k2m
    xyz2m(3)=gi2m*sinmutI*k2m
    deltaxyzm(1)=xyz1m(1)-xyz2m(1)
    deltaxyzm(2)=xyz1m(2)-xyz2m(2)
    deltaxyzm(3)=xyz1m(3)-xyz2m(3)
! give a sign to AMOIDm
     CALL sign_dmin(xyz1m,tau1m,xyz2m,tau2m,dsignm)
     amoidm=dsignm*amoidm! here amoidm is always >0
! Derivatives
    dx1barmdcosmutompl=ddnodmdcosom1
    dx2barmdcosmutom=-ddnodmdcosom2
    deffe1mdsmutompl=deffe1pdsmutompl
    deffe2mdsmutom=deffe2pdsmutom
    dgi1mdcmutompl=deffe1mdsmutompl
    dgi2mdcmutom=deffe2mdsmutom
!----------------------------------------------------------------
! Derivatives of Delta x
!----------------------------------------------------------------
      t0 = (1.D0+e1)*(-1.D0+e2*cmutom2)**2.D0*(-1.D0+e1*cmutom1)*sinmutI&
     &**2/(e2**2*smutom2**2-2.D0*e2**2*smutom2**2*e1*cmutom1+e2**2*smuto&
     &m2**2*e1**2*cmutom1**2+e1**2*smutom1**2-2.D0*e1**2*smutom1**2*e2*c&
     &mutom2+e1**2*smutom1**2*e2**2*cmutom2**2+sinmutI**2-2.D0*sinmutI**&
     &2*e2*cmutom2+sinmutI**2*e2**2*cmutom2**2-2.D0*sinmutI**2*e1*cmutom&
     &1+4.D0*sinmutI**2*e1*cmutom1*e2*cmutom2-2.D0*sinmutI**2*e1*cmutom1&
     &*e2**2*cmutom2**2+sinmutI**2*e1**2*cmutom1**2-2.D0*sinmutI**2*e1**&
     &2*cmutom1**2*e2*cmutom2+sinmutI**2*e1**2*cmutom1**2*e2**2*cmutom2*&
     &*2-2.D0*e1*smutom1*e2*smutom2*cosmutI+2.D0*e1*smutom1*e2**2*smutom&
     &2*cosmutI*cmutom2+2.D0*e1**2*smutom1*e2*smutom2*cosmutI*cmutom1-2.&
     &D0*e1**2*smutom1*e2**2*smutom2*cosmutI*cmutom1*cmutom2)           
      ddeltamxdq1=t0 
!                                                                        
      t0 = -(1.D0+e2)*(-1.D0+e1*cmutom1)**2.D0*(-1.D0+e2*cmutom2)*sinmut&
     &I**2/(e2**2*smutom2**2-2.D0*e2**2*smutom2**2*e1*cmutom1+e2**2*smut&
     &om2**2*e1**2*cmutom1**2+e1**2*smutom1**2-2.D0*e1**2*smutom1**2*e2*&
     &cmutom2+e1**2*smutom1**2*e2**2*cmutom2**2+sinmutI**2-2.D0*sinmutI*&
     &*2*e2*cmutom2+sinmutI**2*e2**2*cmutom2**2-2.D0*sinmutI**2*e1*cmuto&
     &m1+4.D0*sinmutI**2*e1*cmutom1*e2*cmutom2-2.D0*sinmutI**2*e1*cmutom&
     &1*e2**2*cmutom2**2+sinmutI**2*e1**2*cmutom1**2-2.D0*sinmutI**2*e1*&
     &*2*cmutom1**2*e2*cmutom2+sinmutI**2*e1**2*cmutom1**2*e2**2*cmutom2&
     &**2-2.D0*e1*smutom1*e2*smutom2*cosmutI+2.D0*e1*smutom1*e2**2*smuto&
     &m2*cosmutI*cmutom2+2.D0*e1**2*smutom1*e2*smutom2*cosmutI*cmutom1-2&
     &.D0*e1**2*smutom1*e2**2*smutom2*cosmutI*cmutom1*cmutom2)          
      ddeltamxdq2=t0 
 !                                                                       
      s1 = -sinmutI**2 
      s3 = (-1.D0+e2*cmutom2)**2.D0 
      s7 = -2.D0*q1*cmutom1**2*e2**2*smutom2**2*e1-2.D0*q1*e2**2*smutom2&
     &**2*e1*cmutom1+q1*cmutom1**3*e2**2*smutom2**2*e1**2+q1*sinmutI**2*&
     &e2**2*cmutom2**2+q1*sinmutI**2*e1**2*cmutom1**2-2.D0*q1*cmutom1**2&
     &*sinmutI**2*e1+q1*cmutom1**3*sinmutI**2*e1**2-2.D0*q1*sinmutI**2*e&
     &1*cmutom1-2.D0*q1*cmutom1**3*sinmutI**2*e1**2*e2*cmutom2+q1*cmutom&
     &1**3*sinmutI**2*e1**2*e2**2*cmutom2**2+e2**2*smutom2**2*q1*cmutom1&
     &-2.D0*q1*cmutom1**2*sinmutI**2*e1*e2**2*cmutom2**2+4.D0*q1*cmutom1&
     &**2*sinmutI**2*e1*e2*cmutom2                                      
      s6 = s7-2.D0*q1*cmutom1**2*e1**2*smutom1*e2**2*smutom2*cosmutI*cmu&
     &tom2+2.D0*q1*cmutom1**2*e1**2*smutom1*e2*smutom2*cosmutI-2.D0*q1*c&
     &mutom1*e1**2*smutom1**2*e2*cmutom2-2.D0*smutom1*q2*e2**2*smutom2*c&
     &osmutI-2.D0*smutom1*q2*e2*smutom2*cosmutI+4.D0*smutom1**2*q1*e1*e2&
     &*cmutom2+e2**2*smutom2**2*q1-4.D0*q1*cmutom1*e1*smutom1*e2*smutom2&
     &*cosmutI-2.D0*q1*cmutom1*sinmutI**2*e2*cmutom2+q1*cmutom1*sinmutI*&
     &*2*e2**2*cmutom2**2+q1*cmutom1*e1**2*smutom1**2-2.D0*smutom1**2*q2&
     &*e1**2*cmutom1+2.D0*smutom1**2*q2*e2*e1                           
      s7 = -2.D0*q1*sinmutI**2*e2*cmutom2-2.D0*smutom1**2*q1*e1*e2**2*cm&
     &utom2**2-2.D0*smutom1**2*q2*e2**2*e1*cmutom2-2.D0*smutom1*q2*e1**2&
     &*cmutom1**2*e2*smutom2*cosmutI+4.D0*smutom1*q2*e2**2*smutom2*cosmu&
     &tI*e1*cmutom1+2.D0*smutom1**2*q2*e1**2*cmutom1*e2*cmutom2+q1*sinmu&
     &tI**2-2.D0*smutom1*q1*e2**2*smutom2*cosmutI*cmutom2+2.D0*smutom1*q&
     &1*e2*smutom2*cosmutI+q1*cmutom1*e1**2*smutom1**2*e2**2*cmutom2**2+&
     &4.D0*q1*cmutom1*e1*smutom1*e2**2*smutom2*cosmutI*cmutom2-2.D0*smut&
     &om1*q2*e2**2*e1**2*cmutom1**2*smutom2*cosmutI+2.D0*smutom1**2*q2*e&
     &2**2*e1**2*cmutom1*cmutom2-2.D0*smutom1**2*q2*e2*e1**2*cmutom1    
      s5 = s7+4.D0*smutom1*q2*e2*smutom2*cosmutI*e1*cmutom1-q1*e1**2*smu&
     &tom1**2*e2**2*cmutom2**2+2.D0*q1*e1**2*smutom1**2*e2*cmutom2-2.D0*&
     &q1*sinmutI**2*e1**2*cmutom1**2*e2*cmutom2+q1*e2**2*smutom2**2*e1**&
     &2*cmutom1**2-2.D0*q1*sinmutI**2*e1*cmutom1*e2**2*cmutom2**2+4.D0*q&
     &1*sinmutI**2*e1*cmutom1*e2*cmutom2+q1*sinmutI**2*e1**2*cmutom1**2*&
     &e2**2*cmutom2**2+s6-2.D0*smutom1**2*q2*e1*e2*cmutom2-q1*e1**2*smut&
     &om1**2+2.D0*smutom1**2*q2*e1-2.D0*smutom1**2*q1*e1+q1*cmutom1*sinm&
     &utI**2                                                            
      s6 = 1.D0/((e2**2*smutom2**2-2.D0*e2**2*smutom2**2*e1*cmutom1+e2**&
     &2*smutom2**2*e1**2*cmutom1**2+e1**2*smutom1**2-2.D0*e1**2*smutom1*&
     &*2*e2*cmutom2+e1**2*smutom1**2*e2**2*cmutom2**2+sinmutI**2-2.D0*si&
     &nmutI**2*e2*cmutom2+sinmutI**2*e2**2*cmutom2**2-2.D0*sinmutI**2*e1&
     &*cmutom1+4.D0*sinmutI**2*e1*cmutom1*e2*cmutom2-2.D0*sinmutI**2*e1*&
     &cmutom1*e2**2*cmutom2**2+sinmutI**2*e1**2*cmutom1**2-2.D0*sinmutI*&
     &*2*e1**2*cmutom1**2*e2*cmutom2+sinmutI**2*e1**2*cmutom1**2*e2**2*c&
     &mutom2**2-2.D0*e1*smutom1*e2*smutom2*cosmutI+2.D0*e1*smutom1*e2**2&
     &*smutom2*cosmutI*cmutom2+2.D0*e1**2*smutom1*e2*smutom2*cosmutI*cmu&
     &tom1-2.D0*e1**2*smutom1*e2**2*smutom2*cosmutI*cmutom1*cmutom2)**2.&
     &D0)                                                               
      s4 = s5*s6 
      s2 = s3*s4 
      ddeltamxde1 = s1*s2 
 !                                                                       
      s1 = sinmutI**2*(-1.D0+e1*cmutom1)**2.D0 
      s5 = -e2**2*smutom2**2*q2*e1**2*cmutom1**2+4.D0*e1*smutom1*cmutom2&
     &*q1*e2*smutom2*cosmutI+e1**2*smutom1**2*q2-2.D0*smutom2**2*q1*e2**&
     &2*cmutom2+2.D0*smutom2**2*q1*e1*e2+smutom2**2*q2*e2**2*cmutom2+q2*&
     &sinmutI**2*e1**2*cmutom1**2+q2*cmutom2*sinmutI**2*e1**2*cmutom1**2&
     &+2.D0*e1*smutom1*q2*smutom2*cosmutI+q2*cmutom2**3*sinmutI**2*e1**2&
     &*cmutom1**2*e2**2+4.D0*q2*cmutom2**2*sinmutI**2*e1*cmutom1*e2+e1**&
     &2*smutom1**2*q2*e2**2*cmutom2**2+4.D0*smutom2**2*q2*e2*e1*cmutom1 
      s4 = s5-2.D0*q2*cmutom2**2*e1**2*smutom1*e2**2*smutom2*cosmutI*cmu&
     &tom1-2.D0*q2*cmutom2**2*sinmutI**2*e1**2*cmutom1**2*e2-2.D0*e1*smu&
     &tom1*q1*smutom2*cosmutI+4.D0*e1**2*smutom1*cmutom2*q1*e2*smutom2*c&
     &osmutI+q2*sinmutI**2*e1**2*cmutom1**2*e2**2*cmutom2**2-2.D0*e1**2*&
     &smutom1*q1*smutom2*cosmutI-2.D0*e1**2*smutom1*q2*smutom2*cosmutI*c&
     &mutom1-2.D0*smutom2**2*q1*e1*e2**2*cmutom2-2.D0*q2*sinmutI**2*e1**&
     &2*cmutom1**2*e2*cmutom2+smutom2**2*q2*e1**2*cmutom1**2*e2**2*cmuto&
     &m2+4.D0*q2*sinmutI**2*e1*cmutom1*e2*cmutom2-2.D0*q2*cmutom2*e2**2*&
     &smutom2**2*e1*cmutom1+2.D0*q2*cmutom2**2*e1*smutom1*e2**2*smutom2*&
     &cosmutI                                                           
      s5 = -2.D0*q2*sinmutI**2*e1*cmutom1*e2**2*cmutom2**2-4.D0*e1*smuto&
     &m1*q2*cmutom2*e2*smutom2*cosmutI-2.D0*smutom2**2*q1*e1**2*e2*cmuto&
     &m1-2.D0*e1*smutom1*cmutom2**2*q1*e2**2*smutom2*cosmutI-2.D0*e1**2*&
     &smutom1*cmutom2**2*q1*e2**2*smutom2*cosmutI+e1**2*smutom1**2*q2*cm&
     &utom2**3*e2**2+4.D0*e1**2*smutom1*q2*cmutom2*e2*smutom2*cosmutI*cm&
     &utom1-2.D0*q2*cmutom2**3*sinmutI**2*e1*cmutom1*e2**2-e2**2*smutom2&
     &**2*q2+q2*cmutom2*sinmutI**2-2.D0*smutom2**2*q2*e2+2.D0*smutom2**2&
     &*q1*e2-2.D0*q2*cmutom2*sinmutI**2*e1*cmutom1-2.D0*e1**2*smutom1**2&
     &*q2*cmutom2**2*e2                                                 
      s3 = s5+2.D0*smutom2**2*q1*e2**2*e1*cmutom1*cmutom2+q2*sinmutI**2*&
     &e2**2*cmutom2**2-2.D0*q2*sinmutI**2*e1*cmutom1+q2*cmutom2**3*sinmu&
     &tI**2*e2**2-2.D0*q2*cmutom2**2*sinmutI**2*e2+e1**2*smutom1**2*q2*c&
     &mutom2-2.D0*q2*sinmutI**2*e2*cmutom2+2.D0*e2**2*smutom2**2*q2*e1*c&
     &mutom1+2.D0*smutom2**2*q1*e1**2*e2**2*cmutom2*cmutom1-2.D0*smutom2&
     &**2*q2*e2*e1**2*cmutom1**2-2.D0*smutom2**2*q1*e2*e1*cmutom1+s4-2.D&
     &0*e1**2*smutom1**2*q2*e2*cmutom2+q2*sinmutI**2                    
      s4 = 1.D0/((e2**2*smutom2**2-2.D0*e2**2*smutom2**2*e1*cmutom1+e2**&
     &2*smutom2**2*e1**2*cmutom1**2+e1**2*smutom1**2-2.D0*e1**2*smutom1*&
     &*2*e2*cmutom2+e1**2*smutom1**2*e2**2*cmutom2**2+sinmutI**2-2.D0*si&
     &nmutI**2*e2*cmutom2+sinmutI**2*e2**2*cmutom2**2-2.D0*sinmutI**2*e1&
     &*cmutom1+4.D0*sinmutI**2*e1*cmutom1*e2*cmutom2-2.D0*sinmutI**2*e1*&
     &cmutom1*e2**2*cmutom2**2+sinmutI**2*e1**2*cmutom1**2-2.D0*sinmutI*&
     &*2*e1**2*cmutom1**2*e2*cmutom2+sinmutI**2*e1**2*cmutom1**2*e2**2*c&
     &mutom2**2-2.D0*e1*smutom1*e2*smutom2*cosmutI+2.D0*e1*smutom1*e2**2&
     &*smutom2*cosmutI*cmutom2+2.D0*e1**2*smutom1*e2*smutom2*cosmutI*cmu&
     &tom1-2.D0*e1**2*smutom1*e2**2*smutom2*cosmutI*cmutom1*cmutom2)**2.&
     &D0)                                                               
      s2 = s3*s4 
      ddeltamxde2 = s1*s2 
!                                                                        
    ddeltamxdi1=dx1barmdcosmutompl*dcosom(1,1,1)-dx2barmdcosmutom*dcosom(2,1,1)+ &
         & k2m*deffe2mdsmutom*dsindcos(2)*dcosom(2,1,1)-k1m* &
         & deffe1mdsmutompl*dsindcos(1)*dcosom(1,1,1)
    ddeltamxdi2=dx1barmdcosmutompl*dcosom(1,1,2)-dx2barmdcosmutom*dcosom(2,1,2)+ &
         & k2m*deffe2mdsmutom*dsindcos(2)*dcosom(2,1,2)-k1m* &
         & deffe1mdsmutompl*dsindcos(1)*dcosom(1,1,2)

    ddeltamxdbigom1=dx1barmdcosmutompl*dcosom(1,3,1)-dx2barmdcosmutom* &
         & dcosom(2,3,1)+k2m*deffe2mdsmutom*dsindcos(2)*dcosom(2,3,1)-k1m* &
         & deffe1mdsmutompl*dsindcos(1)*dcosom(1,3,1)
    ddeltamxdbigom2=dx1barmdcosmutompl*dcosom(1,3,2)-dx2barmdcosmutom* &
         & dcosom(2,3,2)+k2m*deffe2mdsmutom*dsindcos(2)*dcosom(2,3,2)-k1m* &
         & deffe1mdsmutompl*dsindcos(1)*dcosom(1,3,2)

    ddeltamxdom1=dx1barmdcosmutompl*dcosom(1,2,1)-dx2barmdcosmutom* &
         & dcosom(2,2,1)+k2m*deffe2mdsmutom*dsindcos(2)*dcosom(2,2,1)-k1m* &
         & deffe1mdsmutompl*dsindcos(1)*dcosom(1,2,1)
    ddeltamxdom2=dx1barmdcosmutompl*dcosom(1,2,2)-dx2barmdcosmutom* &
         & dcosom(2,2,2)+k2m*deffe2mdsmutom*dsindcos(2)*dcosom(2,2,2)-k1m* &
         & deffe1mdsmutompl*dsindcos(1)*dcosom(1,2,2)
!----------------------------------------------------------------
! Derivatives of Delta y
!----------------------------------------------------------------
      t0 = -(cosmutI-1.D0)*(cosmutI+1.D0)*(-1.D0+e2*cmutom2)**2.D0*smuto&
     &m1*e1*(1.D0+e1)/(e2**2*smutom2**2-2.D0*e2**2*smutom2**2*e1*cmutom1&
     &+e2**2*smutom2**2*e1**2*cmutom1**2+e1**2*smutom1**2-2.D0*e1**2*smu&
     &tom1**2*e2*cmutom2+e1**2*smutom1**2*e2**2*cmutom2**2+sinmutI**2-2.&
     &D0*sinmutI**2*e2*cmutom2+sinmutI**2*e2**2*cmutom2**2-2.D0*sinmutI*&
     &*2*e1*cmutom1+4.D0*sinmutI**2*e1*cmutom1*e2*cmutom2-2.D0*sinmutI**&
     &2*e1*cmutom1*e2**2*cmutom2**2+sinmutI**2*e1**2*cmutom1**2-2.D0*sin&
     &mutI**2*e1**2*cmutom1**2*e2*cmutom2+sinmutI**2*e1**2*cmutom1**2*e2&
     &**2*cmutom2**2-2.D0*e1*smutom1*e2*smutom2*cosmutI+2.D0*e1*smutom1*&
     &e2**2*smutom2*cosmutI*cmutom2+2.D0*e1**2*smutom1*e2*smutom2*cosmut&
     &I*cmutom1-2.D0*e1**2*smutom1*e2**2*smutom2*cosmutI*cmutom1*cmutom2&
     &)                                                                 
      ddeltamydq1=t0 
!
      t0 = (cosmutI-1.D0)*(cosmutI+1.D0)*(-1.D0+e2*cmutom2)*smutom1*e1*(&
     &1.D0+e2)*(-1.D0+e1*cmutom1)/(e2**2*smutom2**2-2.D0*e2**2*smutom2**&
     &2*e1*cmutom1+e2**2*smutom2**2*e1**2*cmutom1**2+e1**2*smutom1**2-2.&
     &D0*e1**2*smutom1**2*e2*cmutom2+e1**2*smutom1**2*e2**2*cmutom2**2+s&
     &inmutI**2-2.D0*sinmutI**2*e2*cmutom2+sinmutI**2*e2**2*cmutom2**2-2&
     &.D0*sinmutI**2*e1*cmutom1+4.D0*sinmutI**2*e1*cmutom1*e2*cmutom2-2.&
     &D0*sinmutI**2*e1*cmutom1*e2**2*cmutom2**2+sinmutI**2*e1**2*cmutom1&
     &**2-2.D0*sinmutI**2*e1**2*cmutom1**2*e2*cmutom2+sinmutI**2*e1**2*c&
     &mutom1**2*e2**2*cmutom2**2-2.D0*e1*smutom1*e2*smutom2*cosmutI+2.D0&
     &*e1*smutom1*e2**2*smutom2*cosmutI*cmutom2+2.D0*e1**2*smutom1*e2*sm&
     &utom2*cosmutI*cmutom1-2.D0*e1**2*smutom1*e2**2*smutom2*cosmutI*cmu&
     &tom1*cmutom2)                                                     
      ddeltamydq2=t0 
!
     s1 = smutom1*(cosmutI-1.D0)*(cosmutI+1.D0) 
      s3 = -1.D0+e2*cmutom2 
      s7 = -2.D0*q1*e1**2*e2**2*smutom2**2*cmutom1+q1*sinmutI**2+6.D0*q1&
     &*e1*e2**2*cmutom2**2*sinmutI**2-6.D0*e1**2*cmutom1*q1*e2**2*cmutom&
     &2**2*sinmutI**2-6.D0*q1*e1*e2*cmutom2*sinmutI**2+2.D0*e1**2*cmutom&
     &1*q1*e2**3*cmutom2**3*sinmutI**2+2.D0*q1*e1**2*e2**3*cmutom2*smuto&
     &m2**2*cmutom1-2.D0*q1*e2*smutom2*cosmutI*e1**2*smutom1-2.D0*q1*e2*&
     &*3*smutom2*cosmutI*e1**2*smutom1*cmutom2**2+2.D0*cmutom1*q2*e2**3*&
     &e1*smutom2**2-cmutom1**2*q2*e2*e1**2*sinmutI**2-cmutom1**2*q2*e2**&
     &3*e1**2*sinmutI**2*cmutom2**2+2.D0*cmutom1**2*q2*e2**2*e1**2*sinmu&
     &tI**2*cmutom2-2.D0*cmutom1*q1*e2**3*smutom2*cosmutI*e1**2*smutom1*&
     &cmutom2**2+3.D0*q1*e1**2*smutom1**2*e2*cmutom2-q2*sinmutI**2*e1**2&
     &*cmutom1**2*e2**2*cmutom2**2                                      
      s6 = s7+e1**2*smutom1**2*q2*e2**2*cmutom2**2-2.D0*e1**2*smutom1**2&
     &*q2*e2*cmutom2+2.D0*q2*sinmutI**2*e1**2*cmutom1**2*e2*cmutom2-e2**&
     &2*smutom2**2*q2*e1**2*cmutom1**2+6.D0*e1**2*cmutom1*q1*e2*cmutom2*&
     &sinmutI**2+q2*e2**3*smutom1**2*e1**2*cmutom2**2+q1*smutom1**2*e2**&
     &3*cmutom2**3*e1**2+cmutom1**2*q1*e1**2*e2**3*cmutom2*smutom2**2-2.&
     &D0*cmutom1*q1*e2*smutom2*cosmutI*e1**2*smutom1-2.D0*q2*e2**2*smuto&
     &m1**2*cmutom2*e1**2+2.D0*e2**2*smutom2**2*q2*e1*cmutom1+2.D0*q2*si&
     &nmutI**2*e1*cmutom1*e2**2*cmutom2**2-4.D0*q2*sinmutI**2*e1*cmutom1&
     &*e2*cmutom2+4.D0*cmutom1*q1*e2**2*smutom2*cosmutI*e1**2*smutom1*cm&
     &utom2+2.D0*cmutom1*q2*e2*e1*sinmutI**2+2.D0*cmutom1*q2*e2**3*e1*si&
     &nmutI**2*cmutom2**2                                               
      s7 = s6-4.D0*cmutom1*q2*e2**2*e1*sinmutI**2*cmutom2-q2*sinmutI**2-&
     &cmutom1**2*q2*e2**3*e1**2*smutom2**2+cmutom1**2*e1**2*q1*e2**3*cmu&
     &tom2**3*sinmutI**2-2.D0*q1*e1*e2**3*cmutom2**3*sinmutI**2-2.D0*q1*&
     &e1*e2**3*cmutom2*smutom2**2+4.D0*q1*e2**2*smutom2*cosmutI*e1**2*sm&
     &utom1*cmutom2-3.D0*q1*e1**2*smutom1**2*e2**2*cmutom2**2-q1*e2**2*s&
     &mutom2**2*e1**2*cmutom1**2+3.D0*q1*sinmutI**2*e1**2*cmutom1**2*e2*&
     &cmutom2-3.D0*q1*sinmutI**2*e1**2*cmutom1**2*e2**2*cmutom2**2-3.D0*&
     &q1*sinmutI**2*e2*cmutom2-q1*sinmutI**2*e1**2*cmutom1**2+3.D0*q1*si&
     &nmutI**2*e2**2*cmutom2**2+2.D0*q2*sinmutI**2*e1*cmutom1           
      s5 = s7+2.D0*q2*sinmutI**2*e2*cmutom2-q2*sinmutI**2*e2**2*cmutom2*&
     &*2-q2*sinmutI**2*e1**2*cmutom1**2-q2*e2**3*sinmutI**2*cmutom2**2+2&
     &.D0*q2*e2**2*sinmutI**2*cmutom2-q1*e2**3*cmutom2**3*sinmutI**2-q1*&
     &e2**3*cmutom2*smutom2**2+q2*e2*smutom1**2*e1**2+2.D0*q1*e1*e2**2*s&
     &mutom2**2-2.D0*e1**2*cmutom1*q1*sinmutI**2-q2*e2**3*smutom2**2+2.D&
     &0*q1*e1*sinmutI**2-q2*e2*sinmutI**2-e2**2*smutom2**2*q2+e1**2*smut&
     &om1**2*q2-q1*e1**2*smutom1**2+e2**2*smutom2**2*q1                 
      s6 = 1.D0/((e2**2*smutom2**2-2.D0*e2**2*smutom2**2*e1*cmutom1+e2**&
     &2*smutom2**2*e1**2*cmutom1**2+e1**2*smutom1**2-2.D0*e1**2*smutom1*&
     &*2*e2*cmutom2+e1**2*smutom1**2*e2**2*cmutom2**2+sinmutI**2-2.D0*si&
     &nmutI**2*e2*cmutom2+sinmutI**2*e2**2*cmutom2**2-2.D0*sinmutI**2*e1&
     &*cmutom1+4.D0*sinmutI**2*e1*cmutom1*e2*cmutom2-2.D0*sinmutI**2*e1*&
     &cmutom1*e2**2*cmutom2**2+sinmutI**2*e1**2*cmutom1**2-2.D0*sinmutI*&
     &*2*e1**2*cmutom1**2*e2*cmutom2+sinmutI**2*e1**2*cmutom1**2*e2**2*c&
     &mutom2**2-2.D0*e1*smutom1*e2*smutom2*cosmutI+2.D0*e1*smutom1*e2**2&
     &*smutom2*cosmutI*cmutom2+2.D0*e1**2*smutom1*e2*smutom2*cosmutI*cmu&
     &tom1-2.D0*e1**2*smutom1*e2**2*smutom2*cosmutI*cmutom1*cmutom2)**2.&
     &D0)                                                               
      s4 = s5*s6 
      s2 = s3*s4 
      ddeltamyde1 = s1*s2 
!
      s1 = -e1*smutom1*(cosmutI-1.D0) 
      s3 = (cosmutI+1.D0)*(-1.D0+e1*cmutom1) 
      s6 = -2.D0*smutom2**2*q2*e2-2.D0*q2*sinmutI**2*e1*cmutom1+q2*cmuto&
     &m2*sinmutI**2-2.D0*q2*sinmutI**2*e2*cmutom2-e2**2*smutom2**2*q2+q2&
     &*sinmutI**2*e2**2*cmutom2**2+q2*sinmutI**2*e1**2*cmutom1**2+smutom&
     &2**2*q2*e2**2*cmutom2+e1**2*smutom1**2*q2+2.D0*smutom2**2*q1*e1*e2&
     &-2.D0*smutom2**2*q1*e2**2*cmutom2-2.D0*q2*cmutom2**2*sinmutI**2*e2&
     &+2.D0*smutom2**2*q1*e2+q2*cmutom2**3*sinmutI**2*e2**2+e1**2*smutom&
     &1**2*q2*cmutom2+4.D0*e1*smutom1*cmutom2*q1*e2*smutom2*cosmutI+q2*s&
     &inmutI**2+4.D0*e1**2*smutom1*q2*cmutom2*e2*smutom2*cosmutI*cmutom1&
     &-2.D0*q2*cmutom2**3*sinmutI**2*e1*cmutom1*e2**2+e1**2*smutom1**2*q&
     &2*cmutom2**3*e2**2-2.D0*e1**2*smutom1**2*q2*cmutom2**2*e2-2.D0*q2*&
     &cmutom2*sinmutI**2*e1*cmutom1+2.D0*smutom2**2*q1*e2**2*e1*cmutom1*&
     &cmutom2+2.D0*smutom2**2*q1*e1**2*e2**2*cmutom2*cmutom1+2.D0*e2**2*&
     &smutom2**2*q2*e1*cmutom1-e2**2*smutom2**2*q2*e1**2*cmutom1**2     
      s7 = s6-2.D0*smutom2**2*q1*e2*e1*cmutom1-2.D0*smutom2**2*q2*e2*e1*&
     &*2*cmutom1**2-2.D0*e1*smutom1*q1*smutom2*cosmutI+2.D0*e1*smutom1*q&
     &2*smutom2*cosmutI+q2*cmutom2*sinmutI**2*e1**2*cmutom1**2+q2*cmutom&
     &2**3*sinmutI**2*e1**2*cmutom1**2*e2**2-2.D0*q2*cmutom2**2*sinmutI*&
     &*2*e1**2*cmutom1**2*e2-2.D0*q2*cmutom2**2*e1**2*smutom1*e2**2*smut&
     &om2*cosmutI*cmutom1-2.D0*e1**2*smutom1**2*q2*e2*cmutom2+e1**2*smut&
     &om1**2*q2*e2**2*cmutom2**2+4.D0*smutom2**2*q2*e2*e1*cmutom1+4.D0*q&
     &2*cmutom2**2*sinmutI**2*e1*cmutom1*e2-2.D0*e1**2*smutom1*cmutom2**&
     &2*q1*e2**2*smutom2*cosmutI                                        
      s5 = s7-2.D0*e1*smutom1*cmutom2**2*q1*e2**2*smutom2*cosmutI+4.D0*e&
     &1**2*smutom1*cmutom2*q1*e2*smutom2*cosmutI-2.D0*smutom2**2*q1*e1**&
     &2*e2*cmutom1-2.D0*e1**2*smutom1*q1*smutom2*cosmutI-2.D0*e1**2*smut&
     &om1*q2*smutom2*cosmutI*cmutom1+q2*sinmutI**2*e1**2*cmutom1**2*e2**&
     &2*cmutom2**2-2.D0*smutom2**2*q1*e1*e2**2*cmutom2-4.D0*e1*smutom1*q&
     &2*cmutom2*e2*smutom2*cosmutI+smutom2**2*q2*e1**2*cmutom1**2*e2**2*&
     &cmutom2+4.D0*q2*sinmutI**2*e1*cmutom1*e2*cmutom2-2.D0*q2*sinmutI**&
     &2*e1**2*cmutom1**2*e2*cmutom2+2.D0*q2*cmutom2**2*e1*smutom1*e2**2*&
     &smutom2*cosmutI-2.D0*q2*sinmutI**2*e1*cmutom1*e2**2*cmutom2**2-2.D&
     &0*q2*cmutom2*e2**2*smutom2**2*e1*cmutom1                          
      s6 = 1.D0/((e2**2*smutom2**2-2.D0*e2**2*smutom2**2*e1*cmutom1+e2**&
     &2*smutom2**2*e1**2*cmutom1**2+e1**2*smutom1**2-2.D0*e1**2*smutom1*&
     &*2*e2*cmutom2+e1**2*smutom1**2*e2**2*cmutom2**2+sinmutI**2-2.D0*si&
     &nmutI**2*e2*cmutom2+sinmutI**2*e2**2*cmutom2**2-2.D0*sinmutI**2*e1&
     &*cmutom1+4.D0*sinmutI**2*e1*cmutom1*e2*cmutom2-2.D0*sinmutI**2*e1*&
     &cmutom1*e2**2*cmutom2**2+sinmutI**2*e1**2*cmutom1**2-2.D0*sinmutI*&
     &*2*e1**2*cmutom1**2*e2*cmutom2+sinmutI**2*e1**2*cmutom1**2*e2**2*c&
     &mutom2**2-2.D0*e1*smutom1*e2*smutom2*cosmutI+2.D0*e1*smutom1*e2**2&
     &*smutom2*cosmutI*cmutom2+2.D0*e1**2*smutom1*e2*smutom2*cosmutI*cmu&
     &tom1-2.D0*e1**2*smutom1*e2**2*smutom2*cosmutI*cmutom1*cmutom2)**2.&
     &D0)                                                               
      s4 = s5*s6 
      s2 = s3*s4 
      ddeltamyde2 = s1*s2 
! 
    ddeltamydi1=k1m*dgi1mdcmutompl*dcosom(1,1,1)-k2m*(dcosmutI(1,1)*gi2m+ &
         & cosmutI*dgi2mdcmutom*dcosom(2,1,1))
    ddeltamydi2=k1m*dgi1mdcmutompl*dcosom(1,1,2)-k2m*(dcosmutI(1,2)*gi2m+ &
         & cosmutI*dgi2mdcmutom*dcosom(2,1,2))
!
    ddeltamydbigom1=k1m*dgi1mdcmutompl*dcosom(1,3,1)-k2m* &
         & (dcosmutI(3,1)*gi2m+cosmutI*dgi2mdcmutom*dcosom(2,3,1))
    ddeltamydbigom2=k1m*dgi1mdcmutompl*dcosom(1,3,2)-k2m* &
         & (dcosmutI(3,2)*gi2m+cosmutI*dgi2mdcmutom*dcosom(2,3,2))
!
    ddeltamydom1=k1m*dgi1mdcmutompl*dcosom(1,2,1)-k2m* &
         & (dcosmutI(2,1)*gi2m+cosmutI*dgi2mdcmutom*dcosom(2,2,1))
    ddeltamydom2=k1m*dgi1mdcmutompl*dcosom(1,2,2)-k2m* &
         & (dcosmutI(2,2)*gi2m+cosmutI*dgi2mdcmutom*dcosom(2,2,2))
!----------------------------------------------------------------
! Derivatives of Delta z
!----------------------------------------------------------------
      t0 = -(-e1*smutom1*cosmutI+e1*smutom1*cosmutI*e2*cmutom2+e2*smutom&
     &2-e2*smutom2*e1*cmutom1)*(1.D0+e1)*sinmutI*(-1.D0+e2*cmutom2)/(e2*&
     &*2*smutom2**2-2.D0*e2**2*smutom2**2*e1*cmutom1+e2**2*smutom2**2*e1&
     &**2*cmutom1**2+e1**2*smutom1**2-2.D0*e1**2*smutom1**2*e2*cmutom2+e&
     &1**2*smutom1**2*e2**2*cmutom2**2+sinmutI**2-2.D0*sinmutI**2*e2*cmu&
     &tom2+sinmutI**2*e2**2*cmutom2**2-2.D0*sinmutI**2*e1*cmutom1+4.D0*s&
     &inmutI**2*e1*cmutom1*e2*cmutom2-2.D0*sinmutI**2*e1*cmutom1*e2**2*c&
     &mutom2**2+sinmutI**2*e1**2*cmutom1**2-2.D0*sinmutI**2*e1**2*cmutom&
     &1**2*e2*cmutom2+sinmutI**2*e1**2*cmutom1**2*e2**2*cmutom2**2-2.D0*&
     &e1*smutom1*e2*smutom2*cosmutI+2.D0*e1*smutom1*e2**2*smutom2*cosmut&
     &I*cmutom2+2.D0*e1**2*smutom1*e2*smutom2*cosmutI*cmutom1-2.D0*e1**2&
     &*smutom1*e2**2*smutom2*cosmutI*cmutom1*cmutom2)                   
      ddeltamzdq1=t0 
!                                                                        
      t0 = (-e1*smutom1*cosmutI+e1*smutom1*cosmutI*e2*cmutom2+e2*smutom2&
     &-e2*smutom2*e1*cmutom1)*(-1.D0+e1*cmutom1)*(1.D0+e2)*sinmutI/(e2**&
     &2*smutom2**2-2.D0*e2**2*smutom2**2*e1*cmutom1+e2**2*smutom2**2*e1*&
     &*2*cmutom1**2+e1**2*smutom1**2-2.D0*e1**2*smutom1**2*e2*cmutom2+e1&
     &**2*smutom1**2*e2**2*cmutom2**2+sinmutI**2-2.D0*sinmutI**2*e2*cmut&
     &om2+sinmutI**2*e2**2*cmutom2**2-2.D0*sinmutI**2*e1*cmutom1+4.D0*si&
     &nmutI**2*e1*cmutom1*e2*cmutom2-2.D0*sinmutI**2*e1*cmutom1*e2**2*cm&
     &utom2**2+sinmutI**2*e1**2*cmutom1**2-2.D0*sinmutI**2*e1**2*cmutom1&
     &**2*e2*cmutom2+sinmutI**2*e1**2*cmutom1**2*e2**2*cmutom2**2-2.D0*e&
     &1*smutom1*e2*smutom2*cosmutI+2.D0*e1*smutom1*e2**2*smutom2*cosmutI&
     &*cmutom2+2.D0*e1**2*smutom1*e2*smutom2*cosmutI*cmutom1-2.D0*e1**2*&
     &smutom1*e2**2*smutom2*cosmutI*cmutom1*cmutom2)                    
      ddeltamzdq2=t0    
!
      s1 = sinmutI*(-1.D0+e2*cmutom2) 
      s6 = q2*e2**3*smutom1**3*cosmutI*e1**2*cmutom2**2+2.D0*q1*e2**3*sm&
     &utom2**3*e1*cmutom1**2+2.D0*q2*e2**2*smutom2*e1*smutom1**2*cmutom2&
     &+q2*e2**3*smutom1*cosmutI*smutom2**2+2.D0*q2*e2**3*smutom2*e1*smut&
     &om1**2*cmutom2-2.D0*q2*e2**2*smutom1**3*cosmutI*cmutom2*e1**2+q2*e&
     &2*smutom1**3*cosmutI*e1**2+q2*smutom1**3*cosmutI*e2**2*cmutom2**2*&
     &e1**2-2.D0*q2*smutom1**3*cosmutI*e2*cmutom2*e1**2-2.D0*q2*e2**2*sm&
     &utom2*e1*smutom1**2-q1*e2*smutom2*sinmutI**2-2.D0*q2*e2*smutom2*e1&
     &*smutom1**2+2.D0*q2*e2**2*smutom1*cosmutI*sinmutI**2*cmutom2      
      s5 = s6-2.D0*q1*e2**3*cmutom2**3*e1*smutom1*cosmutI*sinmutI**2+2.D&
     &0*q1*e1*smutom1*cosmutI*e2**2*smutom2**2+4.D0*q1*e1**2*smutom1**2*&
     &cosmutI**2*e2**2*smutom2*cmutom2+2.D0*q1*e1**2*cmutom1**2*e2**2*cm&
     &utom2*smutom2*sinmutI**2+q1*e2**3*smutom2*e1**2*smutom1**2*cmutom2&
     &**2-3.D0*q1*cmutom1**2*e1**2*smutom1*cosmutI*e2**2*smutom2**2-2.D0&
     &*q1*e1**2*cmutom1*smutom1*cosmutI*e2**2*smutom2**2+2.D0*q1*cmutom1&
     &**2*e2**3*cmutom2**2*smutom2*e1*sinmutI**2-q1*cmutom1*e2**3*smutom&
     &2**3+2.D0*q1*e1**2*cmutom1*smutom1*cosmutI*e2**3*cmutom2**3*sinmut&
     &I**2-4.D0*q1*cmutom1**2*e2**2*cmutom2*smutom2*e1*sinmutI**2+4.D0*q&
     &1*cmutom1*e1*smutom1*cosmutI*e2**2*smutom2**2-q1*cmutom1*e2**3*cmu&
     &tom2**2*smutom2*e1**2*smutom1**2                                  
      s6 = 2.D0*q1*cmutom1*e2**2*cmutom2*smutom2*e1**2*smutom1**2-q1*e1*&
     &*2*cmutom1**2*e2**3*cmutom2**2*smutom2*sinmutI**2-4.D0*q1*cmutom1*&
     &e1*smutom1*cosmutI*e2**3*cmutom2*smutom2**2+6.D0*q1*e1**2*cmutom1*&
     &smutom1*cosmutI*e2*cmutom2*sinmutI**2-q1*cmutom1*e2*smutom2*e1**2*&
     &smutom1**2+q1*smutom1*cosmutI*sinmutI**2-6.D0*q1*e1**2*cmutom1*smu&
     &tom1*cosmutI*e2**2*cmutom2**2*sinmutI**2-q2*e2*smutom1*cosmutI*sin&
     &mutI**2-4.D0*q1*e2**2*smutom2*e1*smutom1**2*cmutom2+3.D0*q1*cmutom&
     &1**2*e1**2*smutom1*cosmutI*e2**3*cmutom2*smutom2**2-q1*cmutom1*e2*&
     &smutom2*sinmutI**2+2.D0*q1*cmutom1*e2**2*cmutom2*smutom2*sinmutI**&
     &2-q1*e2*smutom2*e1**2*cmutom1**3*sinmutI**2-q1*cmutom1*e2**3*cmuto&
     &m2**2*smutom2*sinmutI**2                                          
      s4 = s6+2.D0*cmutom1*q2*e2**2*smutom2*e1**2*smutom1**2-q1*e2**3*sm&
     &utom2**3-4.D0*cmutom1*q2*e1*smutom1*cosmutI*sinmutI**2*e2*cmutom2-&
     &2.D0*cmutom1*q2*e1*smutom1*cosmutI*e2**2*smutom2**2-2.D0*q1*cmutom&
     &1*e1**2*smutom1**2*cosmutI**2*e2*smutom2-4.D0*q1*e2**2*cmutom2*smu&
     &tom2*e1*cmutom1*sinmutI**2-q1*cmutom1**2*e1**2*smutom1*cosmutI*sin&
     &mutI**2-2.D0*q1*e1**2*smutom1**2*cosmutI**2*e2*smutom2-2.D0*q1*e2*&
     &*2*smutom2*e1**2*smutom1**2*cmutom2+2.D0*q1*e2**3*cmutom2**2*smuto&
     &m2*e1*cmutom1*sinmutI**2+6.D0*q1*e2**2*cmutom2**2*e1*smutom1*cosmu&
     &tI*sinmutI**2-cmutom1**2*q2*e1**2*smutom1*cosmutI*sinmutI**2+2.D0*&
     &q1*e1**2*cmutom1*smutom1*cosmutI*e2**3*cmutom2*smutom2**2+s5      
      s6 = q1*e2*smutom2*e1**2*smutom1**2+2.D0*q1*e2**3*smutom2**3*e1*cm&
     &utom1+2.D0*q1*e2**2*smutom2*sinmutI**2*cmutom2+2.D0*q1*e1*smutom1*&
     &cosmutI*sinmutI**2-3.D0*q1*e2**2*cmutom2**2*smutom1**3*cosmutI*e1*&
     &*2-q2*e2**3*smutom1*cosmutI*sinmutI**2*cmutom2**2-3.D0*q1*cmutom1*&
     &*2*e2**2*cmutom2**2*e1**2*smutom1*cosmutI*sinmutI**2-2.D0*cmutom1*&
     &q2*e2**3*smutom2*e1**2*smutom1**2*cmutom2-2.D0*q1*e2**3*cmutom2**2&
     &*e1**2*smutom1**2*cosmutI**2*smutom2-cmutom1**2*q2*e1**2*smutom1*c&
     &osmutI*sinmutI**2*e2**2*cmutom2**2-2.D0*q1*cmutom1*e2**3*cmutom2**&
     &2*e1**2*smutom1**2*cosmutI**2*smutom2+3.D0*q1*cmutom1**2*e1**2*smu&
     &tom1*cosmutI*e2*cmutom2*sinmutI**2-q1*cmutom1**3*e2**3*cmutom2**2*&
     &smutom2*e1**2*sinmutI**2                                          
      s5 = s6-q1*e2**3*smutom2*sinmutI**2*cmutom2**2+4.D0*q1*cmutom1*e1*&
     &*2*smutom1**2*cosmutI**2*e2**2*cmutom2*smutom2+q1*cmutom1**2*e2**3&
     &*cmutom2**3*e1**2*smutom1*cosmutI*sinmutI**2+2.D0*cmutom1*q2*e2*sm&
     &utom2*e1**2*smutom1**2+2.D0*q1*cmutom1**3*e2**2*cmutom2*smutom2*e1&
     &**2*sinmutI**2+2.D0*cmutom1*q2*e1*smutom1*cosmutI*sinmutI**2*e2**2&
     &*cmutom2**2-6.D0*q1*e1*smutom1*cosmutI*sinmutI**2*e2*cmutom2-2.D0*&
     &q1*e1*smutom1*cosmutI*e2**3*cmutom2*smutom2**2+2.D0*q1*e2*smutom2*&
     &e1*cmutom1*sinmutI**2+2.D0*cmutom1*q2*e1*smutom1*cosmutI*sinmutI**&
     &2+2.D0*cmutom1**2*q2*e1**2*smutom1*cosmutI*sinmutI**2*e2*cmutom2+c&
     &mutom1**2*q2*e1**2*smutom1*cosmutI*e2**2*smutom2**2+s4-q1*e2**3*sm&
     &utom2**3*e1**2*cmutom1**3                                         
      s6 = s5-2.D0*q1*e1**2*cmutom1*smutom1*cosmutI*sinmutI**2-cmutom1**&
     &2*q2*e2*e1**2*smutom1*cosmutI*sinmutI**2-q1*e1**2*cmutom1**2*e2*sm&
     &utom2*sinmutI**2+cmutom1**2*q2*e2**3*e1**2*smutom1*cosmutI*smutom2&
     &**2+2.D0*cmutom1*q2*e2**3*e1*smutom1*cosmutI*sinmutI**2*cmutom2**2&
     &-4.D0*cmutom1*q2*e2**2*e1*smutom1*cosmutI*sinmutI**2*cmutom2+q2*sm&
     &utom1**3*cosmutI*e1**2-cmutom1**2*q2*e2**3*e1**2*smutom1*cosmutI*s&
     &inmutI**2*cmutom2**2-q1*smutom1**3*cosmutI*e1**2+2.D0*cmutom1*q2*e&
     &2*e1*smutom1*cosmutI*sinmutI**2-2.D0*cmutom1*q2*e2**3*e1*smutom1*c&
     &osmutI*smutom2**2+2.D0*cmutom1**2*q2*e2**2*e1**2*smutom1*cosmutI*s&
     &inmutI**2*cmutom2-2.D0*cmutom1*q2*e2**2*smutom2*e1**2*smutom1**2*c&
     &mutom2                                                            
      s3 = s6+q2*smutom1*cosmutI*e2**2*smutom2**2+2.D0*q1*e2*smutom2*e1*&
     &cmutom1**2*sinmutI**2-q2*smutom1*cosmutI*sinmutI**2-q1*smutom1*cos&
     &mutI*e2**2*smutom2**2+3.D0*q1*smutom1**3*cosmutI*e2*cmutom2*e1**2+&
     &q1*smutom1*cosmutI*e2**3*cmutom2*smutom2**2-q1*e1**2*cmutom1**2*e2&
     &**3*smutom2**3+2.D0*q1*e2*smutom2*e1*smutom1**2+2.D0*q1*e2**3*smut&
     &om2*e1*smutom1**2*cmutom2**2-q2*smutom1*cosmutI*e2**2*cmutom2**2*s&
     &inmutI**2+2.D0*q2*smutom1*cosmutI*e2*cmutom2*sinmutI**2-3.D0*q1*sm&
     &utom1*cosmutI*e2*cmutom2*sinmutI**2-q1*e2**3*cmutom2**3*smutom1*co&
     &smutI*sinmutI**2+3.D0*q1*e2**2*cmutom2**2*smutom1*cosmutI*sinmutI*&
     &*2+q1*e2**3*cmutom2**3*smutom1**3*cosmutI*e1**2                   
      s4 = 1.D0/((e2**2*smutom2**2-2.D0*e2**2*smutom2**2*e1*cmutom1+e2**&
     &2*smutom2**2*e1**2*cmutom1**2+e1**2*smutom1**2-2.D0*e1**2*smutom1*&
     &*2*e2*cmutom2+e1**2*smutom1**2*e2**2*cmutom2**2+sinmutI**2-2.D0*si&
     &nmutI**2*e2*cmutom2+sinmutI**2*e2**2*cmutom2**2-2.D0*sinmutI**2*e1&
     &*cmutom1+4.D0*sinmutI**2*e1*cmutom1*e2*cmutom2-2.D0*sinmutI**2*e1*&
     &cmutom1*e2**2*cmutom2**2+sinmutI**2*e1**2*cmutom1**2-2.D0*sinmutI*&
     &*2*e1**2*cmutom1**2*e2*cmutom2+sinmutI**2*e1**2*cmutom1**2*e2**2*c&
     &mutom2**2-2.D0*e1*smutom1*e2*smutom2*cosmutI+2.D0*e1*smutom1*e2**2&
     &*smutom2*cosmutI*cmutom2+2.D0*e1**2*smutom1*e2*smutom2*cosmutI*cmu&
     &tom1-2.D0*e1**2*smutom1*e2**2*smutom2*cosmutI*cmutom1*cmutom2)**2.&
     &D0)                                                               
      s2 = s3*s4 
      ddeltamzde1 = s1*s2 
!
      s1 = -sinmutI 
      s3 = -1.D0+e1*cmutom1 
      s8 = q1*e1*smutom2*sinmutI**2-q2*smutom2*e1**2*smutom1**2-q1*e1*sm&
     &utom2**3*e2**2+q1*smutom2*e1**2*smutom1**2-2.D0*q2*e2*smutom2*sinm&
     &utI**2+q1*e1**3*smutom2*smutom1**2-2.D0*q2*e1**3*smutom1**3*cosmut&
     &I*e2*cmutom2+2.D0*q1*e1**2*smutom2**3*cmutom1*e2**2-2.D0*q1*smutom&
     &2*e1*cmutom1*sinmutI**2+q1*e1**3*smutom2*cmutom1**2*sinmutI**2+q1*&
     &smutom2*sinmutI**2-q2*smutom2*sinmutI**2-2.D0*q1*smutom2**2*e1**2*&
     &cmutom1*smutom1*e2*cosmutI                                        
      s7 = s8-2.D0*cmutom2*q1*e2*smutom2*sinmutI**2-2.D0*cmutom2*q1*e1**&
     &2*smutom1*cosmutI*e2**2*smutom2**2+4.D0*cmutom2*q1*e1**2*smutom1**&
     &2*cosmutI**2*e2*smutom2+q2*smutom2**3*e2**2-2.D0*cmutom2**3*e2**2*&
     &cmutom1*q2*e1**2*smutom1*cosmutI*sinmutI**2-2.D0*q1*e1**2*smutom2*&
     &cmutom1*sinmutI**2-2.D0*q1*e1**3*smutom1**2*cosmutI**2*smutom2-2.D&
     &0*q1*e1**2*smutom1**2*cosmutI**2*smutom2-2.D0*cmutom1*q2*e1**2*smu&
     &tom1*cosmutI*e2**2*smutom2**2-2.D0*q2*smutom2**2*e1*smutom1*e2*cos&
     &mutI-2.D0*q1*e1**3*smutom2**2*cmutom1*smutom1*e2*cosmutI-2.D0*cmut&
     &om2**2*q1*e1**3*smutom1**2*cosmutI**2*e2**2*smutom2+cmutom2**3*e2*&
     &*2*q2*e1**3*smutom1**3*cosmutI+4.D0*cmutom2*q1*e1**3*smutom1**2*co&
     &smutI**2*e2*smutom2                                               
      s8 = cmutom2**2*q1*e2**2*smutom2*sinmutI**2-2.D0*cmutom2**2*cmutom&
     &1**2*q2*e1**3*smutom1*cosmutI*e2*sinmutI**2-2.D0*cmutom2**2*cmutom&
     &1*q2*e1**3*smutom1**2*cosmutI**2*e2**2*smutom2-q1*smutom2**3*e2**2&
     &-2.D0*cmutom2**2*q2*e1**3*smutom1**3*cosmutI*e2-2.D0*cmutom2**2*q2&
     &*e1*smutom1*cosmutI*e2*sinmutI**2+q2*e1*smutom1*cosmutI*sinmutI**2&
     &*cmutom2+q1*smutom2*sinmutI**2*e1**2*cmutom1**2-2.D0*cmutom1*q2*e1&
     &**3*smutom1**2*cosmutI**2*smutom2+cmutom1**2*q2*e1**3*smutom1*cosm&
     &utI*sinmutI**2+cmutom1**2*q2*e1**3*smutom1*cosmutI*cmutom2*sinmutI&
     &**2-2.D0*cmutom1*q2*e1**2*smutom1*cosmutI*sinmutI**2-2.D0*cmutom1*&
     &q2*e1**2*smutom1*cosmutI*cmutom2*sinmutI**2+3.D0*cmutom2**2*e2**2*&
     &cmutom1**2*q2*e1**2*smutom2*sinmutI**2                            
      s6 = s8-cmutom2**2*e2**2*q2*smutom2*e1**3*cmutom1*smutom1**2-2.D0*&
     &cmutom2*q1*e1**2*cmutom1**2*e2*smutom2*sinmutI**2-4.D0*cmutom2*q2*&
     &e1**2*smutom1**2*cosmutI**2*e2*smutom2+cmutom2**2*e2**2*q2*smutom2&
     &*sinmutI**2+cmutom1**2*q2*e1**3*smutom1*cosmutI*e2**2*smutom2**2+4&
     &.D0*q2*smutom2**2*e1**2*cmutom1*smutom1*e2*cosmutI+q2*e1**3*smutom&
     &1**3*cosmutI-2.D0*q2*smutom2**2*e1**3*cmutom1**2*smutom1*e2*cosmut&
     &I-2.D0*e2**2*cmutom2*q2*smutom2*e1**3*cmutom1*smutom1**2+3.D0*cmut&
     &om2*q2*e1*smutom1*cosmutI*e2**2*smutom2**2+6.D0*e2**2*cmutom2*cmut&
     &om1**2*q2*e1**2*smutom2*sinmutI**2+2.D0*cmutom2**2*q2*e1**2*smutom&
     &1**2*cosmutI**2*e2**2*smutom2+q2*e1*smutom1*cosmutI*sinmutI**2-3.D&
     &0*cmutom1**2*q2*e1**2*smutom2*sinmutI**2+s7  
      s8 = s6-2.D0*e2**2*cmutom2**2*cmutom1*q2*e1**2*smutom1*cosmutI*sin&
     &mutI**2-2.D0*cmutom2*q1*e2*smutom2*e1**2*smutom1**2+3.D0*cmutom2*c&
     &mutom1**2*q2*e1**3*smutom1*cosmutI*e2**2*smutom2**2-6.D0*cmutom2*c&
     &mutom1*q2*e1**2*smutom1*cosmutI*e2**2*smutom2**2-2.D0*cmutom2**2*q&
     &1*e1**2*e2**2*smutom2*cmutom1*sinmutI**2+4.D0*cmutom2**2*cmutom1*q&
     &2*e1**2*smutom1*cosmutI*e2*sinmutI**2+cmutom2**2*e2**2*q2*smutom2*&
     &e1**2*smutom1**2-3.D0*cmutom2**2*e2**2*q2*smutom2*e1*cmutom1*sinmu&
     &tI**2+2.D0*cmutom2*q1*e1**2*cmutom1*smutom1*cosmutI*e2**2*smutom2*&
     &*2-2.D0*cmutom2*q1*e1*e2*smutom2*sinmutI**2+cmutom2**3*e2**2*q2*e1&
     &*smutom1*cosmutI*sinmutI**2+cmutom2**2*q1*e1**2*cmutom1**2*e2**2*s&
     &mutom2*sinmutI**2-2.D0*cmutom2**2*q1*e1**2*smutom1**2*cosmutI**2*e&
     &2**2*smutom2                                                      
      s7 = s8-2.D0*cmutom2*q1*e1*smutom1*cosmutI*e2**2*smutom2**2-2.D0*c&
     &mutom2**2*q1*e2**2*smutom2*e1*cmutom1*sinmutI**2+cmutom2**2*q1*e2*&
     &*2*smutom2*e1**2*smutom1**2+cmutom2**2*q1*e1**3*e2**2*smutom2*smut&
     &om1**2+4.D0*cmutom2*q1*e1**2*e2*smutom2*cmutom1*sinmutI**2+cmutom2&
     &**2*q1*e1*e2**2*smutom2*sinmutI**2+2.D0*cmutom2*q1*e1**3*smutom1*c&
     &osmutI*e2**2*smutom2**2*cmutom1+4.D0*cmutom2*cmutom1*q2*e1**3*smut&
     &om1**2*cosmutI**2*e2*smutom2-2.D0*e2**2*cmutom2*cmutom1**3*q2*e1**&
     &3*smutom2*sinmutI**2+cmutom1**3*q2*e1**3*smutom2*sinmutI**2+q2*e1*&
     &smutom1*cosmutI*e2**2*smutom2**2+2.D0*e2**2*cmutom2*q2*smutom2*e1*&
     &*2*smutom1**2+e2**2*cmutom2**2*cmutom1**2*q2*e1**3*smutom1*cosmutI&
     &*sinmutI**2-6.D0*e2**2*cmutom2*q2*smutom2*e1*cmutom1*sinmutI**2   
      s8 = 6.D0*q2*e2*smutom2*e1*cmutom1*sinmutI**2+2.D0*q2*e2*smutom2*e&
     &1**3*cmutom1*smutom1**2-2.D0*q2*e1*smutom1*cosmutI*e2*cmutom2*sinm&
     &utI**2+q2*e1**3*smutom1**3*cosmutI*cmutom2+cmutom2**3*e2**2*cmutom&
     &1**2*q2*e1**3*smutom1*cosmutI*sinmutI**2+2.D0*e2**2*cmutom2*q2*smu&
     &tom2*sinmutI**2-2.D0*cmutom2*q1*e1**3*e2*smutom2*cmutom1**2*sinmut&
     &I**2+e2**2*cmutom2**2*q2*e1**3*smutom1**3*cosmutI-cmutom2**2*e2**2&
     &*cmutom1**3*q2*e1**3*smutom2*sinmutI**2+4.D0*cmutom2*q1*e2*smutom2&
     &*e1*cmutom1*sinmutI**2+cmutom2**2*q1*e1**3*e2**2*smutom2*cmutom1**&
     &2*sinmutI**2+3.D0*q2*smutom2*e1*cmutom1*sinmutI**2-2.D0*cmutom2*q1&
     &*e1**3*e2*smutom2*smutom1**2+s7+e2**2*cmutom2**2*q2*e1*smutom1*cos&
     &mutI*sinmutI**2  
      s5 = s8+2.D0*q1*e1**2*smutom2**2*smutom1*e2*cosmutI+4.D0*cmutom1*q&
     &2*e1**2*smutom1*cosmutI*e2*cmutom2*sinmutI**2+2.D0*q1*smutom2**2*e&
     &1*smutom1*e2*cosmutI-2.D0*q2*e2*smutom2*e1**2*smutom1**2-q1*e1**3*&
     &smutom2**3*cmutom1**2*e2**2-2.D0*cmutom1**2*q2*e1**3*smutom1*cosmu&
     &tI*e2*cmutom2*sinmutI**2-q1*smutom2**3*e1**2*cmutom1**2*e2**2-3.D0&
     &*q2*smutom2**3*e1*cmutom1*e2**2+2.D0*q1*smutom2**3*e1*cmutom1*e2**&
     &2-q2*smutom2**3*e1**3*cmutom1**3*e2**2+3.D0*q2*smutom2**3*e1**2*cm&
     &utom1**2*e2**2-6.D0*cmutom1**2*q2*e1**2*e2*smutom2*sinmutI**2+2.D0&
     &*q2*e1**2*smutom1**2*cosmutI**2*smutom2+q2*smutom2*e1**3*cmutom1*s&
     &mutom1**2+2.D0*cmutom1**3*q2*e1**3*e2*smutom2*sinmutI**2          
      s6 = 1.D0/((e2**2*smutom2**2-2.D0*e2**2*smutom2**2*e1*cmutom1+e2**&
     &2*smutom2**2*e1**2*cmutom1**2+e1**2*smutom1**2-2.D0*e1**2*smutom1*&
     &*2*e2*cmutom2+e1**2*smutom1**2*e2**2*cmutom2**2+sinmutI**2-2.D0*si&
     &nmutI**2*e2*cmutom2+sinmutI**2*e2**2*cmutom2**2-2.D0*sinmutI**2*e1&
     &*cmutom1+4.D0*sinmutI**2*e1*cmutom1*e2*cmutom2-2.D0*sinmutI**2*e1*&
     &cmutom1*e2**2*cmutom2**2+sinmutI**2*e1**2*cmutom1**2-2.D0*sinmutI*&
     &*2*e1**2*cmutom1**2*e2*cmutom2+sinmutI**2*e1**2*cmutom1**2*e2**2*c&
     &mutom2**2-2.D0*e1*smutom1*e2*smutom2*cosmutI+2.D0*e1*smutom1*e2**2&
     &*smutom2*cosmutI*cmutom2+2.D0*e1**2*smutom1*e2*smutom2*cosmutI*cmu&
     &tom1-2.D0*e1**2*smutom1*e2**2*smutom2*cosmutI*cmutom1*cmutom2)**2.&
     &D0)                                                               
      s4 = s5*s6 
      s2 = s3*s4 
      ddeltamzde2 = s1*s2 
!
    ddeltamzdi1=-k2m*(dgi2mdcmutom*dcosom(2,1,1)*sinmutI+ &
         & gi2m*dsinmutIdcos*dcosmutI(1,1))
    ddeltamzdi2=-k2m*(dgi2mdcmutom*dcosom(2,1,2)*sinmutI+ &
         & gi2m*dsinmutIdcos*dcosmutI(1,2))
  
    ddeltamzdbigom1=-k2m*(dgi2mdcmutom*dcosom(2,3,1)*sinmutI+ &
         & gi2m*dsinmutIdcos*dcosmutI(3,1))
    ddeltamzdbigom2=-k2m*(dgi2mdcmutom*dcosom(2,3,2)*sinmutI+ &
         & gi2m*dsinmutIdcos*dcosmutI(3,2))

    ddeltamzdom1=-k2m*(dgi2mdcmutom*dcosom(2,2,1)*sinmutI+ &
         & gi2m*dsinmutIdcos*dcosmutI(2,1))
    ddeltamzdom2=-k2m*(dgi2mdcmutom*dcosom(2,2,2)*sinmutI+ &
         & gi2m*dsinmutIdcos*dcosmutI(2,2))
! ----------------------------------------------------------------
! Derivatives of the difference of the Cartesian positions
! w.r.t. the 10 elements
! ----------------------------------------------------------------
    ddeltam(1,1)=ddeltamxdq1
    ddeltam(1,2)=ddeltamydq1
    ddeltam(1,3)=ddeltamzdq1
    ddeltam(2,1)=ddeltamxde1
    ddeltam(2,2)=ddeltamyde1
    ddeltam(2,3)=ddeltamzde1
    ddeltam(3,1)=ddeltamxdi1
    ddeltam(3,2)=ddeltamydi1
    ddeltam(3,3)=ddeltamzdi1
    ddeltam(4,1)=ddeltamxdbigom1
    ddeltam(4,2)=ddeltamydbigom1
    ddeltam(4,3)=ddeltamzdbigom1
    ddeltam(5,1)=ddeltamxdom1
    ddeltam(5,2)=ddeltamydom1
    ddeltam(5,3)=ddeltamzdom1
    ddeltam(6,1)=ddeltamxdq2
    ddeltam(6,2)=ddeltamydq2
    ddeltam(6,3)=ddeltamzdq2
    ddeltam(7,1)=ddeltamxde2
    ddeltam(7,2)=ddeltamyde2
    ddeltam(7,3)=ddeltamzde2
    ddeltam(8,1)=ddeltamxdi2
    ddeltam(8,2)=ddeltamydi2
    ddeltam(8,3)=ddeltamzdi2
    ddeltam(9,1)=ddeltamxdbigom2
    ddeltam(9,2)=ddeltamydbigom2
    ddeltam(9,3)=ddeltamzdbigom2 
    ddeltam(10,1)=ddeltamxdom2
    ddeltam(10,2)=ddeltamydom2 
    ddeltam(10,3)=ddeltamzdom2
! Computation of tau3m in the inertial reference frame
    DO j=1,3   
       newtau3m(j)=DOT_PRODUCT(rotinmut(j,1:3),tau3mhat)
    END DO
! Computation of the derivatives of the AMOIDm 
    DO h=1,10       
       DO k=1,3 
          ddeltamdkep(k)=ddeltam(h,k)
          deltarotm(k)=DOT_PRODUCT(drotinmut(h,k,1:3),deltaxyzm)
       ENDDO
        deramoidm(1,h)=DOT_PRODUCT(tau3mhat,ddeltamdkep) +  &
             & DOT_PRODUCT(newtau3m,deltarotm)
    ENDDO
! ----------------------------------------------------------------
!    IF(chk_der)THEN
! note that here amoidpvar,amoidmvar are >0, while amoidp, amoidm 
! are given a sign; to compute the incremental ratio we decide to give 
! them the same (positive) sign
!       deramoidmdequ2(1,1:5)= &
!            & MATMUL(deramoidm(1,6:10),jackepeq2(1:5,1:5)) 
!       write(*,113)'derivatives of AMOIDM:',deramoidmdequ2(1,1:5)
!       write(*,*)'derivative with incremental ratio: der(amoidp)', &
!            & (amoidpvar - amoidp)/ &
!            & (eps_a2+eps_h2+eps_k2+eps_p2+eps_q2),'  der(amoidm)', &
!            & (amoidmvar - amoidm)/ &
!            & (eps_a2+eps_h2+eps_k2+eps_p2+eps_q2)
!    ENDIF
! ----------------------------------------------------------------
    tderamoidm=TRANSPOSE(deramoidm)
! Computation of the RMS of the AMOIDm 
    amoidmrmsvec=MATMUL(deramoidm,MATMUL(covcom,tderamoidm)) 
    amoidmrms=SQRT(amoidmrmsvec(1,1))

  END SUBROUTINE dnod_amoid_rms_com

    

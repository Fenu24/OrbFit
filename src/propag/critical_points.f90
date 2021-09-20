! ******************************************************************
! ********* M O D U L E   C R I T I C A L _ P O I N T S ************
! *********** for the COMPUTATION of the CRITICAL POINTS *********** 
! ************* of the KEPLERIAN DISTANCE ******** *****************
! **** written by G.F. GRONCHI, G. TOMMEI, A.MILANI (Mar. 2007) ****
! ==================================================================
! WITH THE FOLLOWING ROUTINES:
! ROUTINE moid_rms,Habs_rms,crit_pts,orbitcoe,matrixdat,compmodsylv16
! ROUTINE solvesystem,hessian,int_eval,d2eval
! ROUTINE sign_dmin,mutual_ref,CP_newton_raphson,dmintil_rms,check_rms
! ROUTINE comp_rms_com,q_Q_T_rms,tau1_tau2
! ROUTINE derdmintest,sin_mutI,det_H,sintau1tau2
! ===================================================================
MODULE critical_points
  USE fund_const
  USE orbit_elements
  USE output_control
  IMPLICIT NONE
  PRIVATE

  INTEGER, PARAMETER, PUBLIC :: nminx = 7 ! max num of minima
  INTEGER, PARAMETER, PUBLIC :: poldeg = 16 ! polynomial degree
! to write all the complex roots of 1-D poly r(w)
  LOGICAL, PUBLIC :: follow_roots ! follow_roots = .false. do not write
                  ! follow_roots = .true. write the roots of r(w)
! elements in the first 2 rows of the cracovians (see Sitarski 1968)
  REAL(KIND=8) :: PPx,PPy,PPz,QQx,QQy,QQz
  REAL(KIND=8) :: px,py,pz,qx,qy,qz
  REAL(KIND=8) :: A1x,A1y,A1z,B1x,B1y,B1z,a2x,a2y,a2z,b2x,b2y,b2z 
!                
  REAL(KIND=8) :: KK,LL,MM,NN
  REAL(KIND=8) :: Ppl,pcom,Epl,ecom
  REAL(KIND=8) :: a2Q1,b2Q1,B1a2,B1b2,A1a2,A1b2
  REAL(KIND=8) :: A1q2,B1q2,b2A1,b2B1,a2A1,a2B1
! ---------------- Sylvester matrix elements ----------
  INTEGER, PARAMETER :: Nev = 16 ! number of evaluations of the 15th deg poly
  INTEGER, PARAMETER :: expo = 4 ! 2^expo = Nev

  REAL(KIND=8),DIMENSION(Nev) :: p0,p1,p2
  REAL(KIND=8),DIMENSION(Nev) :: q0,q1,q2,q3,q4
  REAL(KIND=8),DIMENSION(Nev) :: r31,r32,r33,r34,r35
  REAL(KIND=8),DIMENSION(Nev) :: r36,r41,r45,r46
! -----------------------------------------------------

! public routines
  PUBLIC :: crit_pts,d2eval 
  PUBLIC :: sign_dmin,mutual_ref,dmintil_rms,q_Q_T_rms, moid_rms,Habs_rms
CONTAINS

! =================================================================
! MOID_RMS
! main interface to compute the MOID, optionally its covariance
! also with the absolute magnitude covariance
! =================================================================
  SUBROUTINE moid_rms(elem,moid,uncel,covmh,rmsh)
   TYPE(orbit_elem), INTENT(IN) :: elem ! asteroid elements
   DOUBLE PRECISION, INTENT(OUT) :: moid
   TYPE(orb_uncert), INTENT(IN), OPTIONAL :: uncel ! covariance
   DOUBLE PRECISION, INTENT(OUT), OPTIONAL :: covmh(2,2) ! covariance
   DOUBLE PRECISION, INTENT(IN), OPTIONAL :: rmsh ! photometry fit RMS
      ! matrix of (MOID, H), containing also the photometry uncertainty
      ! if rmsh is supplied in input WARNING: errors in photmetry and
      ! in astrometry are assumed to be independent, thus their effects
      ! are summed in a RMS sense
! end interface
   TYPE(orbit_elem) :: elea
   TYPE(orb_uncert) :: uncea
   INTEGER nummin
   DOUBLE PRECISION,DIMENSION(5,nminx) :: ddmintdel2
   DOUBLE PRECISION,DIMENSION(6) :: dHdel2
   DOUBLE PRECISION dmintil(nminx),dmintrms(nminx), Hrms
   DOUBLE PRECISION,DIMENSION(5) :: vectmp
! Earth's elements at the same time                                    
   elea=undefined_orbit_elem
   elea%t=elem%t
   elea%coo='CAR'
   CALL earcar(elea%t,elea%coord,1)
   IF(PRESENT(uncel).and.PRESENT(covmh))THEN
      covmh=0.d0
      IF(PRESENT(rmsh))THEN
         covmh(2,2)=rmsh**2
      ENDIF
      CALL undefined_orb_uncert(6,uncea)
      uncea%g=0.d0
! compute MOID with sign and its variance
      CALL dmintil_rms(elea,elem,nummin,dmintil,UNC1=uncea,UNC2=uncel,&
&          DMINTRMS=dmintrms,DDMINTDEL2=ddmintdel2)
      moid=abs(dmintil(1))
      covmh(1,1)= dmintrms(1)**2
      CALL Habs_rms(elea,elem,uncel,Hrms,dHdel2)
      covmh(2,2)= covmh(2,2)+Hrms**2
      vectmp=MATMUL(uncel%g(1:5,1:5),ddmintdel2(1:5,1))
! WARNING: we are using ONLY the first 5 elements of dHdel2
      covmh(1,2)= DOT_PRODUCT(dHdel2(1:5),vectmp)
      covmh(2,1)= covmh(1,2)
   ELSE
! compute MOID with sign
      CALL dmintil_rms(elea,elem,nummin,dmintil)
      moid=abs(dmintil(1))
   ENDIF
 END SUBROUTINE moid_rms

! **************************************************************
! RMS of the absolute magnitude
! **************************************************************
! written by 
! **************************************************************
  SUBROUTINE Habs_rms(el1,el2,unc2,Hrms,dHdel2)
! =====================INTERFACE=======================
    TYPE(orbit_elem), INTENT(IN) :: el1,el2 ! elements of the Earth 
                                            ! and of the asteroid
    TYPE(orb_uncert), INTENT(IN) :: unc2    ! covariance of the asteroid
    DOUBLE PRECISION, INTENT(OUT) :: Hrms   ! RMS of H due to astrometry
    DOUBLE PRECISION, INTENT(OUT), OPTIONAL :: dHdel2(6)
! ------ end interface ---------------------------------------
    TYPE(orbit_elem) :: elatt2 !ATT elements of the asteroid
    TYPE(orbit_elem) :: elcar1,elcar2 !CAR elements
    DOUBLE PRECISION, DIMENSION(6,6) :: dcardel2, dattdel2,gel
    DOUBLE PRECISION, DIMENSION(6) :: drhodel2,drdel2,dbetadel2,dHtmp 
    DOUBLE PRECISION :: rho,dHdrho,r,dHdr,cbeta,dHdbeta,Hvar
    DOUBLE PRECISION :: dcbetadr,tb2,phi1,phi2,gmag
    DOUBLE PRECISION :: dphi1db,dphi2db
    DOUBLE PRECISION, DIMENSION(3) :: dbetadx
    DOUBLE PRECISION :: vsize,rea_sun
    DOUBLE PRECISION,PARAMETER :: a1=3.33d0, a2=1.87d0 
    DOUBLE PRECISION,PARAMETER :: b1=0.63d0, b2=1.22d0 
    INTEGER  fail_flag
! ============================================================
    CALL coo_cha(el2,'ATT',elatt2,fail_flag,dattdel2) ! to compute rho
    CALL coo_cha(el2,'CAR',elcar2,fail_flag,dcardel2) ! to compute r
! geocentric distance and partials
    rho = elatt2%coord(5)
    drhodel2=dattdel2(5,1:6)
! heliocentric distance and partials
    r=vsize(elcar2%coord(1:3))
    drdel2=MATMUL(elcar2%coord(1:3)*(1.d0/r),dcardel2(1:3,1:6))
! derivatives w.r. to rho,r
    dHdrho=-5.d0/(log(1.d1)*rho)
    dHdr=-5.d0/(log(1.d1)*r)
! computation of beta (phate)
    CALL coo_cha(el1,'CAR',elcar1,fail_flag) ! used only for beta (phase)
    rea_sun=vsize(elcar1%coord(1:3))
    cbeta=(r**2+rho**2-rea_sun**2)/(2.d0*r*rho)
!    beta=acos(cbeta) ! beta in [0, pi]
    dcbetadr=(r**2-rho**2+rea_sun**2)/(2*rho*r**2)
    dbetadx=(-1.d0/sqrt(1.d0-cbeta**2)*dcbetadr/r)*elcar2%coord(1:3)
! WARNING: missing dependency on phase
!    tb2=TAN(beta/2.d0) 
    tb2=sqrt((1.d0-cbeta)/(1.d0+cbeta))
    phi1=EXP(-a1*tb2**b1) 
    phi2=EXP(-a2*tb2**b2) 
    gmag=elatt2%g_mag
! ***************************************************************
!    appmag=5.d0*LOG10(ds*dt)+h-2.5d0*LOG10((1.d0-g)*phi1+g*phi2) 
! *** add reference to Bowell paper for this formula ***
    dphi1db=-a1*b1*phi1*tb2**(b1-1.d0)/(1.d0+cbeta)
    dphi2db=-a2*b2*phi2*tb2**(b2-1.d0)/(1.d0+cbeta)
    dHdbeta=2.5d0/(log(10.d0)*((1.d0-gmag)*phi1+gmag*phi2))* &
         & ((1.d0-gmag)*dphi1db+gmag*dphi2db)  
    dbetadel2=MATMUL(dbetadx,dcardel2(1:3,1:6))
! total derivatives of H w.r. to original elements
    dHdel2=dHdrho*drhodel2+dHdr*drdel2 +dHdbeta*dbetadel2
! ------ covariance computation ------
    gel = unc2%g(1:6,1:6)
    dHtmp = MATMUL(gel,dHdel2)
    Hvar = DOT_PRODUCT(dHdel2,dHtmp)
    Hrms = sqrt(Hvar)
  END SUBROUTINE Habs_rms

! ******************************************************************
! ** MAIN ROUTINE for the COMPUTATION of the STATIONARY POINTS **
! *********** written by GIOVANNI F. GRONCHI (Oct.2004) ************
! last modified 31/01/2007 GFG
! ==================================================================
  SUBROUTINE crit_pts(elpl,elcom,fpl,fcom,nroots, &
       & nummin,nummax,answer,warnflag,sflag,morse,weier,hzflag,hwflag, & 
       & multfl,hevalflag)
! cometary elements of the planet and of the comet
    REAL(KIND=8),INTENT(IN),DIMENSION(5) :: elpl
    REAL(KIND=8),INTENT(IN),DIMENSION(5) :: elcom 
    REAL(KIND=8),INTENT(OUT) :: fpl(poldeg),fcom(poldeg)! true anomalies
                                                        ! in [-180,180]  
    INTEGER,INTENT(OUT) :: nroots,nummin,nummax ! number of roots, 
                                                ! minima and maxima 
                                                !(ERROR: nroots=-1)
! ------- error flags ----------------------------------------------
    INTEGER,INTENT(OUT) :: answer(poldeg) ! type of critical point:
!                                         answer =  1      maximum
!                                         answer =  0      saddle
!                                         answer = -1      minimum     
!                                         answer = -2      cannot decide
    LOGICAL,INTENT(INOUT) :: morse ! check with Morse theory 
!                        morse = true      OK   
!                        morse = false     ERROR  
    LOGICAL,INTENT(INOUT) :: weier ! check with Weierstrass theory:
!                        weier = true      OK 
!                        weier = false     ERROR
    LOGICAL,INTENT(INOUT) :: warnflag(3) ! program warning flags:
! warnflag(1) = true    OK
! warnflag(1) = false   leading coefficient of \tilde{r}(t) is very small
! warnflag(2) = true    OK
! warnflag(2) = false   higher degree terms \tilde{r}(t) are not small
! warnflag(3) = true    OK
! warnflag(3) = false   low precision for \tilde{r}(t) coefficients     
    LOGICAL,INTENT(INOUT) :: sflag(6) ! solving system messages:
! sflag(1) = true    OK
! sflag(1) = false   there are two good solutions             
! sflag(2) = true    OK
! sflag(2) = false   neither of the two evaluation is close to 0
! sflag(3) = true    OK
! sflag(3) = false   the s-component of the solution is complex
! sflag(4) = true    OK
! sflag(4) = false   leading coeff. of the 2nd degree pol. is small
! sflag(5) = true    OK
! sflag(5) = false   1st and 2nd coeffs of 2nd degree poly are small
! sflag(6) = true    OK
! sflag(6) = false   2nd degree poly have all coefficients small
    LOGICAL,INTENT(INOUT) :: hzflag ! hzflag = .true.  OK!
                                    ! hzflag = .false. abs(root)>10^5
    LOGICAL,INTENT(INOUT) :: hwflag ! hwflag = .true.  OK!
                                    ! hwflag = .false. abs(root)>10^5
    LOGICAL,INTENT(INOUT) :: multfl ! multfl = .true.  OK!
                                    ! multfl = .false. 0 has multiplicity > 4
    LOGICAL,INTENT(INOUT) :: hevalflag ! hevalflag = .true. OK!
                                       ! hevalflag = .false. unsuccessful 
                                       !                     hessian evaluation
! =========== end interface ========================================
    DOUBLE PRECISION :: fpltmp(poldeg),fcomtmp(poldeg),radtmp(poldeg) 
    INTEGER :: anstmp(poldeg),nstattmp
!
    REAL(KIND=8),DIMENSION(5) :: elc1,elc2 
    INTEGER :: nsol
! resultant evaluations (Nev) and coefficients
    REAL(KIND=8) :: polycoe(Nev),polycoe1(Nev+1),evalpoly(Nev),evalp0(Nev) 
!
    COMPLEX(KIND=8) :: evpol15(Nev)
! for a test of the coefficient found                               
    REAL(KIND=8) :: evalrt(Nev),testevalpoly(Nev) 
    REAL(KIND=8) :: zzero(poldeg) ! roots of the resultant
    REAL(KIND=8) :: radius(poldeg)! error estimate  
    REAL(KIND=8) :: wzero(poldeg) ! corresponding value of s 
! for the eccentric anomalies                                       
    REAL(KIND=8) :: sfpl(poldeg),cfpl(poldeg) 
    REAL(KIND=8) :: sfcom(poldeg),cfcom(poldeg) 
! for the angle translations
    REAL(KIND=8) :: alpha,sa,ca,beta,sb,cb
!    REAL(KIND=8) :: fpltil,sfpltil,cfpltil,fcomtil,sfcomtil,cfcomtil
    REAL(KIND=8) :: sx,cx,sy,cy
! Sylvester matrix elements                                         
    REAL(KIND=8) :: evp0(Nev),evp1(Nev),evp2(Nev)
    REAL(KIND=8) :: evq0(Nev),evq1(Nev),evq2(Nev),evq3(Nev),evq4(Nev)
    REAL(KIND=8) :: evr31(Nev),evr32(Nev),evr33(Nev),evr34(Nev),evr35(Nev)
    REAL(KIND=8) :: evr36(Nev),evr41(Nev),evr45(Nev),evr46(Nev)
    COMPLEX(KIND=8) :: cp0,cp1,cp2 
    COMPLEX(KIND=8) :: cq0,cq1,cq2,cq3,cq4
    COMPLEX(KIND=8) :: cr31,cr32,cr33,cr34,cr35
    COMPLEX(KIND=8) :: cr36,cr41,cr45,cr46

! constant term of det(S_1)
    COMPLEX(KIND=8) :: S10(Nev,Nev),detS10
! complex variables                                                  
    COMPLEX(KIND=8) :: complp0(Nev),complp1(Nev),complp2(Nev)
    COMPLEX(KIND=8) :: complq0(Nev),complq1(Nev),complq2(Nev)
    COMPLEX(KIND=8) :: complq3(Nev),complq4(Nev)
    COMPLEX(KIND=8) :: complr31(Nev),complr32(Nev),complr33(Nev)
    COMPLEX(KIND=8) :: complr34(Nev),complr35(Nev),complr36(Nev)
    COMPLEX(KIND=8) :: complr41(Nev),complr45(Nev),complr46(Nev)
! Sylvester matrix 
    REAL(KIND=8) :: SYLV(6,6,Nev)
! complex variables                                                 
    COMPLEX(KIND=8) :: evalSYLV(6,6,Nev),evalSYLVj(6,6) 
    COMPLEX(KIND=8) :: complSYLVline(6) 
    COMPLEX(KIND=8) :: detevalSYLV,deteval(Nev) ! for the determinant 
! N-th roots of unity
    REAL(KIND=8) :: puroot(Nev),puroottil(Nev)
    COMPLEX(KIND=8) :: compluroot(Nev)

    COMPLEX(KIND=8) :: z_1,z_2,z_3,z_4,z_5,z_6
    COMPLEX(KIND=8) :: w_1,w_2,w_3,w_4,w_5,w_6

    REAL(KIND=8) :: s,t ! polynomial system variables 
    INTEGER :: count ! counter
    INTEGER :: j,l,nr,jj ! loop indexes  
    INTEGER :: h,k ! matrix indexes 
    INTEGER :: ncoe,leftcoe ! indexes for the powers of the resultant
    INTEGER :: ans !  ans = 1: maximum 
!                     ans = 0: saddle 
!                     ans= -1: minimum 
!                     ans= -2: cannot decide
!   ==================================================================
    elc1(1:5)=elpl(1:5)
    elc2(1:5)=elcom(1:5)
! angular shift initialization
    IF(elc1(2).ge.1.d0.or.elc2(2).ge.1.d0) THEN
       ! no angular shift if one orbit is unbounded
       alpha = 0.d0
       beta = 0.d0
    ELSE 
       ! better begin with an angular shift if both orbits are bounded
       alpha = 0.5d0
       beta = 0.5d0
    ENDIF
    count = 0 ! counter initialization
!    count = 9 ! to skip the shift

!    goto 10 ! jump writing the following roots
!    if(verb_moid.ge.20) then
!       write(*,*)'1.d0-elc1(2)*cos(alpha)=',1.d0-elc1(2)*cos(alpha)
!       z_1= (elc1(2)*sin(alpha)+sqrt(elc1(2)**2-1.d0))/ &
!            & (1.d0-elc1(2)*cos(alpha))
!       write(*,*)'z_1=',z_1
!       write(*,*)'z_2=',(elc1(2)*sin(alpha)+sqrt(elc1(2)**2-1.d0))/ &
!            & (1.d0-elc1(2)*cos(alpha))
!       write(*,*)'z_3=',(cos(alpha)-1.d0)/sin(alpha)
!       write(*,*)'z_4=',(cos(alpha)+1.d0)/sin(alpha)
!       write(*,*)'z_5=',(sin(alpha)+sqrt(1.d0-elc1(2)**2))/(elc1(2)+cos(alpha))
!       write(*,*)'z_6=',(sin(alpha)-sqrt(1.d0-elc1(2)**2))/(elc1(2)+cos(alpha))
!       
!       write(*,*)'w_1=',(elc2(2)*sin(beta)+sqrt(elc2(2)**2-1.d0)) &
!            & /(1.d0-elc2(2)*cos(beta))
!       write(*,*)'w_2=',(elc2(2)*sin(beta)+sqrt(elc2(2)**2-1.d0)) &
!            & /(1.d0-elc2(2)*cos(beta))
!       write(*,*)'w_3=',(cos(beta)-1.d0)/sin(beta)
!       write(*,*)'w_4=',(cos(beta)+1.d0)/sin(beta)
!       write(*,*)'w_5=',(sin(beta)+sqrt(1.d0-elc2(2)**2))/(elc2(2)+cos(beta))
!       write(*,*)'w_6=',(sin(beta)-sqrt(1.d0-elc2(2)**2))/(elc2(2)+cos(beta))
!       
!       write(*,*)'z_+=',(cos(alpha)+1.d0)/sin(alpha)
!       write(*,*)'w_+=',(cos(beta)+1.d0)/sin(beta)
!       write(*,*)'z_-=',(cos(alpha)-1.d0)/sin(alpha)
!       write(*,*)'w_-=',(cos(beta)-1.d0)/sin(beta)
!    endif

10  CONTINUE
    CALL orbitcoe(alpha,beta,elc1,elc2)
! initialization
    evp0(1:Nev) = 0.d0
    evp1(1:Nev) = 0.d0
    evp2(1:Nev) = 0.d0
    evq0(1:Nev) = 0.d0
    evq1(1:Nev) = 0.d0
    evq2(1:Nev) = 0.d0
    evq3(1:Nev) = 0.d0
    evq4(1:Nev) = 0.d0
    evr31(1:Nev) = 0.d0
    evr32(1:Nev) = 0.d0
    evr33(1:Nev) = 0.d0
    evr34(1:Nev) = 0.d0
    evr35(1:Nev) = 0.d0
    evr36(1:Nev) = 0.d0
    evr41(1:Nev) = 0.d0
    evr45(1:Nev) = 0.d0
    evr46(1:Nev) = 0.d0
! initialization for polynomial r(t)
    polycoe(1:Nev)=0.d0 
    polycoe1(1:Nev+1) = 0.d0
    evalrt(1:Nev)=0.d0 
    testevalpoly(1:Nev)=0.d0 
!
    nummin = 0 
    nummax = 0 

! READ MATRIX DATA
    CALL matrixdat(alpha,beta)
! remember p0,p1,p2,q0,q1,q2,q3,q4 and use
! evp0,evp1,evp2,evq0,evq1,evq2,evq3,evq4 as
! in/out variables
    evp0(1:5) = p0(1:5)
    evp1(1:5) = p1(1:5)
    evp2(1:5) = p2(1:5)
    evq0(1:5) = q0(1:5)
    evq1(1:5) = q1(1:5)
    evq2(1:5) = q2(1:5)
    evq3(1:5) = q3(1:5)
    evq4(1:5) = q4(1:5)    
    evr31(1:5) = r31(1:5)
    evr32(1:5) = r32(1:5)
    evr33(1:5) = r33(1:5)
    evr34(1:5) = r34(1:5)
    evr35(1:5) = r35(1:5)
    evr36(1:5) = r36(1:5)
    evr41(1:5) = r41(1:5)
    evr45(1:5) = r45(1:5)
    evr46(1:5) = r46(1:5)

! COMPUTE CONSTANT TERM of r(t)
    cp0 = evp0(1)
    cp1 = evp1(1)
    cp2 = evp2(1)
    cq0 = evq0(1)
    cq1 = evq1(1)
    cq2 = evq2(1)
    cq3 = evq3(1)
    cq4 = evq4(1)
    cr31 = evr31(1) 
    cr32 = evr32(1) 
    cr33 = evr33(1) 
    cr34 = evr34(1) 
    cr35 = evr35(1) 
    cr36 = evr36(1) 
    cr41 = evr41(1) 
    cr45 = evr45(1) 
    cr46 = evr46(1) 
    CALL compmodsylv16(cp0,cp1,cp2,cq0,cq1,cq2,cq3,cq4,cr31,cr32,cr33, &
         & cr34,cr35,cr36,cr41,cr45,cr46,S10)
    CALL cdetcomp(S10,detS10)

! EVALUATE the COEFFICIENTS of the MATRIX \tilde{S}_1
    CALL rvfft(evp0,Nev,expo) 
    CALL rvfft(evp1,Nev,expo) 
    CALL rvfft(evp2,Nev,expo) 
    CALL rvfft(evq0,Nev,expo) 
    CALL rvfft(evq1,Nev,expo) 
    CALL rvfft(evq3,Nev,expo) 
    CALL rvfft(evq4,Nev,expo) 
    CALL rvfft(evr31,Nev,expo) 
    CALL rvfft(evr32,Nev,expo)
    CALL rvfft(evr33,Nev,expo)
    CALL rvfft(evr34,Nev,expo)
    CALL rvfft(evr35,Nev,expo) ! is constant: avoid evaluation
    CALL rvfft(evr36,Nev,expo) ! is constant: avoid evaluation
    CALL rvfft(evr41,Nev,expo) 
    CALL rvfft(evr45,Nev,expo) ! is constant: avoid evaluation
    CALL rvfft(evr46,Nev,expo) ! is constant: avoid evaluation

    complp0(1) = DCMPLX(evp0(1),0.d0)
    complp1(1) = DCMPLX(evp1(1),0.d0)
    complp2(1) = DCMPLX(evp2(1),0.d0)
    complq0(1) = DCMPLX(evq0(1),0.d0)
    complq1(1) = DCMPLX(evq1(1),0.d0)
    complq2(1) = DCMPLX(evq2(1),0.d0)
    complq3(1) = DCMPLX(evq3(1),0.d0)
    complq4(1) = DCMPLX(evq4(1),0.d0) 
    complr31(1) = DCMPLX(evr31(1),0.d0) 
    complr32(1) = DCMPLX(evr32(1),0.d0) 
    complr33(1) = DCMPLX(evr33(1),0.d0) 
    complr34(1) = DCMPLX(evr34(1),0.d0) 
    complr35(1) = DCMPLX(evr35(1),0.d0) 
    complr36(1) = DCMPLX(evr36(1),0.d0) 
    complr41(1) = DCMPLX(evr41(1),0.d0) 
    complr45(1) = DCMPLX(evr45(1),0.d0) 
    complr46(1) = DCMPLX(evr46(1),0.d0) 
!
    complp0(Nev/2+1) = DCMPLX(evp0(Nev/2+1),0.d0)
    complp1(Nev/2+1) = DCMPLX(evp1(Nev/2+1),0.d0)
    complp2(Nev/2+1) = DCMPLX(evp2(Nev/2+1),0.d0)
    complq0(Nev/2+1) = DCMPLX(evq0(Nev/2+1),0.d0)
    complq1(Nev/2+1) = DCMPLX(evq1(Nev/2+1),0.d0)
    complq2(Nev/2+1) = DCMPLX(evq2(Nev/2+1),0.d0)
    complq3(Nev/2+1) = DCMPLX(evq3(Nev/2+1),0.d0)
    complq4(Nev/2+1) = DCMPLX(evq4(Nev/2+1),0.d0)
    complr31(Nev/2+1) = DCMPLX(evr31(Nev/2+1),0.d0) 
    complr32(Nev/2+1) = DCMPLX(evr32(Nev/2+1),0.d0) 
    complr33(Nev/2+1) = DCMPLX(evr33(Nev/2+1),0.d0) 
    complr34(Nev/2+1) = DCMPLX(evr34(Nev/2+1),0.d0) 
    complr35(Nev/2+1) = DCMPLX(evr35(Nev/2+1),0.d0) 
    complr36(Nev/2+1) = DCMPLX(evr36(Nev/2+1),0.d0) 
    complr41(Nev/2+1) = DCMPLX(evr41(Nev/2+1),0.d0) 
    complr45(Nev/2+1) = DCMPLX(evr45(Nev/2+1),0.d0) 
    complr46(Nev/2+1) = DCMPLX(evr46(Nev/2+1),0.d0) 
!
    DO j = 1,Nev/2-1
       complp0(j+1) = DCMPLX(evp0(j+1),evp0(Nev-j+1))
       complp0(Nev/2+j+1) = DCMPLX(evp0(Nev/2-j+1),-evp0(Nev/2+1+j))        

       complp1(j+1) = DCMPLX(evp1(j+1),evp1(Nev-j+1))
       complp1(Nev/2+j+1) = DCMPLX(evp1(Nev/2-j+1),-evp1(Nev/2+1+j))        

       complp2(j+1) = DCMPLX(evp2(j+1),evp2(Nev-j+1))
       complp2(Nev/2+j+1) = DCMPLX(evp2(Nev/2-j+1),-evp2(Nev/2+1+j))        

       complq0(j+1) = DCMPLX(evq0(j+1),evq0(Nev-j+1))
       complq0(Nev/2+j+1) = DCMPLX(evq0(Nev/2-j+1),-evq0(Nev/2+1+j))        

       complq1(j+1) = DCMPLX(evq1(j+1),evq1(Nev-j+1))
       complq1(Nev/2+j+1) = DCMPLX(evq1(Nev/2-j+1),-evq1(Nev/2+1+j))        

       complq2(j+1) = DCMPLX(evq2(j+1),evq2(Nev-j+1))
       complq2(Nev/2+j+1) = DCMPLX(evq2(Nev/2-j+1),-evq2(Nev/2+1+j))        

       complq3(j+1) = DCMPLX(evq3(j+1),evq3(Nev-j+1))
       complq3(Nev/2+j+1) = DCMPLX(evq3(Nev/2-j+1),-evq3(Nev/2+1+j))        

       complq4(j+1) = DCMPLX(evq4(j+1),evq4(Nev-j+1))
       complq4(Nev/2+j+1) = DCMPLX(evq4(Nev/2-j+1),-evq4(Nev/2+1+j))        

       complr31(j+1) = DCMPLX(evr31(j+1),evr31(Nev-j+1))
       complr31(Nev/2+j+1) = DCMPLX(evr31(Nev/2-j+1),-evr31(Nev/2+1+j))        

       complr32(j+1) = DCMPLX(evr32(j+1),evr32(Nev-j+1))
       complr32(Nev/2+j+1) = DCMPLX(evr32(Nev/2-j+1),-evr32(Nev/2+1+j))
        
       complr33(j+1) = DCMPLX(evr33(j+1),evr33(Nev-j+1))
       complr33(Nev/2+j+1) = DCMPLX(evr33(Nev/2-j+1),-evr33(Nev/2+1+j))        

       complr34(j+1) = DCMPLX(evr34(j+1),evr34(Nev-j+1))
       complr34(Nev/2+j+1) = DCMPLX(evr34(Nev/2-j+1),-evr34(Nev/2+1+j))        

       complr35(j+1) = DCMPLX(evr35(j+1),evr35(Nev-j+1))
       complr35(Nev/2+j+1) = DCMPLX(evr35(Nev/2-j+1),-evr35(Nev/2+1+j))        

       complr36(j+1) = DCMPLX(evr36(j+1),evr36(Nev-j+1))
       complr36(Nev/2+j+1) = DCMPLX(evr36(Nev/2-j+1),-evr36(Nev/2+1+j))        

       complr41(j+1) = DCMPLX(evr41(j+1),evr41(Nev-j+1))
       complr41(Nev/2+j+1) = DCMPLX(evr41(Nev/2-j+1),-evr41(Nev/2+1+j))        

       complr45(j+1) = DCMPLX(evr45(j+1),evr45(Nev-j+1))
       complr45(Nev/2+j+1) = DCMPLX(evr45(Nev/2-j+1),-evr45(Nev/2+1+j))        

       complr46(j+1) = DCMPLX(evr46(j+1),evr46(Nev-j+1))
       complr46(Nev/2+j+1) = DCMPLX(evr46(Nev/2-j+1),-evr46(Nev/2+1+j))        
    ENDDO
! p(x)=x POLYNOMIAL (for roots of unity)
    puroot(1) = 0.d0
    puroot(2) = 1.d0
    DO j = 3,Nev
       puroot(j) = 0.d0
    ENDDO
! === computing N-th roots of unity ===
    CALL rvfft(puroot,Nev,expo) 
    compluroot(1) = DCMPLX(puroot(1),0.d0)
    compluroot(Nev/2+1) = DCMPLX(puroot(Nev/2+1),0.d0)
    DO j = 1,Nev/2-1
       compluroot(j+1) = DCMPLX(puroot(j+1),puroot(Nev-j+1))
       compluroot(Nev/2+j+1) = DCMPLX(puroot(Nev/2-j+1),-puroot(Nev/2+1+j))        
    ENDDO
! EVALUATED MATRIX evalSYLV(h,k,j)
    DO j = 1,Nev 
       CALL compmodsylv16(complp0(j),complp1(j),complp2(j), &
            & complq0(j),complq1(j),complq2(j),complq3(j),complq4(j), &
            & complr31(j),complr32(j),complr33(j),complr34(j),complr35(j), &
            & complr36(j),complr41(j),complr45(j),complr46(j),evalSYLVj)
       DO h = 1,6 
          DO k = 1,6 
!     defining evaluated matrixes                                       
             evalSYLV(h,k,j) = evalSYLVj(h,k) 
          ENDDO
       ENDDO
    ENDDO
! evaluation of r(t) at the N-th roots o unity
    DO j = 1,Nev 
! defining matrixes whose determinant is to compute
       DO h = 1,6 
          DO k = 1,6 
             evalSYLVj(h,k)=evalSYLV(h,k,j) 
          ENDDO
       ENDDO
! for a given value of j (1<j<Nev) compute det(evalSYLV(h,k,j))
       CALL cdetcomp(evalSYLVj,detevalSYLV)
! evaluations of the resultant in the Nev-th roots of unity
       deteval(j)=detevalSYLV
    ENDDO
! compute evaluation \tilde{r}: (deteval(j)-r0)/uroot(j)
    DO j = 1,Nev
       evpol15(j) = (deteval(j)-detS10)/compluroot(j)
    ENDDO
! compute coefficients of \tilde{r}(t)
    CALL code_input(Nev,evpol15,evalpoly) 
!   remember the evaluations of the resultant for test 2
    DO j = 1,Nev 
       evalrt(j)=evalpoly(j) 
    ENDDO
    CALL irvfft(evalpoly,Nev,expo) 
!   check leading coefficient (poldeg=16)
    IF (abs(evalpoly(poldeg)).le.1d-15) THEN 
       warnflag(1) = .false. 
       if(verb_moid.ge.20) then
          WRITE(ierrou,*)'WARNING! very small leading coefficient '
          WRITE(ierrou,*)'of tilde{r} polynomial (< 1D-15)',evalpoly(poldeg)
          numerr=numerr+1
       endif
       GOTO 13
    ENDIF
! compute coefficients of r(t)
    polycoe1(1) = DBLE(detS10)
    DO ncoe = 1,poldeg
       polycoe1(ncoe+1) = evalpoly(ncoe)
    ENDDO
101 FORMAT (32(es16.6,1x))
!   get a monic polynomial (polycoe1(poldeg+1) = 1)     
    IF(abs(polycoe1(poldeg+1)).ge.1.d-5) THEN
       DO ncoe = 1,poldeg+1 
          polycoe1(ncoe) = polycoe1(ncoe)/polycoe1(poldeg+1) 
       ENDDO
    ELSE
! leave the coeficients unnormalized
       if(verb_moid.ge.20)then
          write(ierrou,*)'leave the coefficients unnormalized'
          write(ierrou,*)'leading coefficient:',polycoe1(poldeg+1)
          numerr=numerr+1
       endif
    ENDIF

! skip test 
    GOTO 111
! TEST:  EVALUATING THE OBTAINED POLYNOMIAL OBTAINED             
    CALL rvfft(testevalpoly,Nev,expo) 
    DO j = 1,Nev 
       IF ( abs(testevalpoly(j)-evalrt(j)).gt.1.d-5) THEN 
!   HINT: polycoe(j) are coded, so we have to test them with evalrt(j)
          warnflag(3) = .false. 
          if (verb_moid.ge.20)then
             WRITE(ierrou,*)'WARNING: TEST FAILED!'
             WRITE(ierrou,*)'DIFFERENCE IN THE EVALUATIONS OF &
                  & THE RESULTANT IS NOT SMALL (> 1.D-5):'                 
             WRITE(ierrou,*)'diff=',abs(testevalpoly(j)-evalrt(j))          
             numerr=numerr+1
          endif
       ENDIF
    ENDDO
111 CONTINUE
                                                                        
! FIND REAL ROOTS OF THE POLYNOMIAL res[ P1(t,s),P2(t,s),t ]        
    CALL solvpoly(poldeg,polycoe1(1:poldeg+1),wzero,nsol,hzflag, &
         & multfl,radius)
    if(.not.hzflag) then
       if(verb_moid.ge.20) then
          write(ierrou,*)'warning: hzflag=',hzflag
          numerr=numerr+1
       endif
       goto 13
    endif

! COMPUTE FOR EACH w THE CORRESPONDING VALUE of z                   
    CALL solvesystem(alpha,beta,nsol,wzero(1:nsol), &
         & zzero(1:nsol),sflag,hwflag) 
    if(.not.hwflag) then
       goto 13
    endif
    sa = sin(alpha)
    ca = cos(alpha)
    sb = sin(beta)
    cb = cos(beta)
    
    DO nr = 1,nsol 
! conversion from (z,w) to the true anomalies (Vpl,vcom)
       sx =  2.d0*zzero(nr)/(1.d0+(zzero(nr))**2) 
       cx = (1.d0-(zzero(nr))**2)/(1.d0+(zzero(nr))**2) 
       sfpl(nr) = sx*ca + cx*sa
       cfpl(nr) = cx*ca - sx*sa
!
       sy = 2.d0*wzero(nr)/(1.d0+(wzero(nr))**2) 
       cy = (1.d0-(wzero(nr))**2)/(1.d0+(wzero(nr))**2) 
       sfcom(nr) = sy*cb + cy*sb
       cfcom(nr) = cy*cb - sy*sb
! critical points (true anomalies) in [-pi,pi]                           
       fpl(nr) = datan2(sfpl(nr),cfpl(nr)) 
       fcom(nr) = datan2(sfcom(nr),cfcom(nr)) 

       DO j=1,2
          CALL CP_newton_raphson(elc1(1:5),elc2(1:5),fpl(nr),fcom(nr))
       ENDDO
    ENDDO

! *********************************************************
! in case of one or two hyperbolas select only one branch
! *********************************************************
    IF((elc1(2).gt.1.d0).or.(elc2(2).gt.1.d0)) THEN
       nstattmp=0
       DO jj = 1,nsol
          IF((1.d0+elc1(2)*cos(radeg*fpl(jj)).gt.0.d0).and. &
               & (1.d0+elc2(2)*cos(radeg*fcom(jj)).gt.0.d0)) THEN
! accept critical point
             nstattmp = nstattmp+1
             fpltmp(nstattmp)=fpl(jj)
             fcomtmp(nstattmp)=fcom(jj)
             radtmp(nstattmp)=radius(jj)
          ENDIF
       ENDDO
!             write(*,*)'nstat,nstattmp',nstat,nstattmp
       nsol=nstattmp
       DO jj=1,nsol
          fpl(jj)=fpltmp(jj)
          fcom(jj)=fcomtmp(jj)
          radius(jj)=radtmp(jj)
       ENDDO
    ENDIF
    nroots=nsol
    IF(nroots.eq.0)THEN
       write(ierrou,*)'crit_pts: ERROR! no critical point found! nroots=',nroots
       numerr=numerr+1
       nummin=0
       nummax=0
       RETURN
    ENDIF

    nummin=0
    nummax=0
    DO nr = 1,nroots 
! compute type of critical points
       CALL hessian(fpl(nr),fcom(nr),ans) 
! compute type of critical points (the determinant is 4 times the previous one)
!       CALL hess_ta_new(fpl(nr),fcom(nr),ans) 
!       ans = -2 !dummy to force int_eval
       IF(ans.eq.-2) THEN
          hevalflag = .false.
!          WRITE(ierrou,*)'cannot decide type of critical point:  &
!               & calling int_eval'
!          numerr=numerr+1
          CALL int_eval(radius(nr),fpl(nr),fcom(nr),ans)
       ENDIF
       IF (ans.eq.-1) THEN 
          nummin = nummin + 1 
       ELSEIF (ans.eq.1)THEN 
          nummax = nummax + 1 
       ENDIF
       answer(nr) = ans 
! conversion into degrees (in [-180,180])                      
       fpl(nr) = fpl(nr)*degrad 
       fcom(nr) =  fcom(nr)*degrad 
    ENDDO

! ***********************
! CHECK the COMPUTATION 
! ***********************
! if both orbits are cometary  ! STILL TO BE WRITTEN (but not interesting...)
    IF (elc1(2).ge.1.d0.and.elc2(2).ge.1.d0) THEN
       WRITE(ierrou,*)'both orbits are cometary! skipping check...'
       numerr=numerr+1
       GOTO 1300 ! skip check
! if at least one orbit is bounded
    ELSE !IF ((elc1(2).lt.1.d0.or.elc2(2).lt.1.d0))THEN
! ============= CHECK WITH MORSE THEORY ===============             
       IF (2*(nummax+nummin).ne.nroots) THEN 
          morse = .false.
       ENDIF
! ============ CHECK WITH WEIERSTRASS THEOREM =========             
       IF(nummin.eq.0)THEN !the absolute minimum should exist 
          weier = .false.
       ELSEIF(nummax.eq.0.and.(elc2(2).lt.1.d0.and.elc1(2).lt.1.d0))THEN
          weier=.false. ! in compatc case the maximum also should exist
       ENDIF
! ============ WRITE ERROR MESSAGE=====================
       IF(.not.morse.or..not.weier)THEN
          IF(verb_moid.gt.20)THEN
             IF(.not.morse)WRITE(ierrou,*)'Morse control failed!'
             IF(.not.weier)WRITE(ierrou,*)'Weierstrass control failed!'
             WRITE(ierrou,*)'nroots,nummin,nummax',nroots,nummin,nummax
!          WRITE(ierrou,*)'ec1 elements:',elc1(1:2),elc1(3:5)*degrad
!          WRITE(ierrou,*)'ec2 elements:',elc2(1:2),elc2(3:5)*degrad
             numerr=numerr+1
             IF(2*(nummin+1).ge.nroots.and.nummin.ge.1)THEN
                WRITE(ierrou,*)'consistent no of minimum points'
             ELSE
                WRITE(ierrou,*)'not consistent no of minimum points'
             ENDIF
             IF(nummin.ge.1)THEN
                WRITE(ierrou,*)'at least one minimum point exists: acceptable solutions'
             ELSE
                WRITE(ierrou,*)'ERROR: there is no minimum point!'
             ENDIF
          ENDIF
       ENDIF
    ENDIF
!    morse=.false. ! dummy: to force the shift
! ======================================================================
! ******************  A N G U L A R   S H I F T  ***********************
! ************** V = \Xi + \alpha  ;  v = \xi + \beta ******************
13  IF((.not.morse).or. &
         & (.not.weier).or. &
         & (.not.warnflag(1)).or.&
!         & (.not.sflag(2)).or. &
!         & (.not.sflag(3)).or. &
         & (.not.sflag(4)).or. &
         & (.not.sflag(5)).or. &
         & (.not.sflag(6)).or. &
         & (.not.hzflag).or. &
         & (.not.hwflag) ) THEN
!         & (.not.multfl)
       count = count + 1
       IF(count.eq.10) THEN
          if(verb_moid.ge.20) then
             WRITE(ierrou,*)'COMPUTATION FAILED:',count,'ANGULAR SHIFT TRIED!'
             WRITE(ierrou,*)'ec1 elements:',elc1(1:2),elc1(3:5)*degrad
             WRITE(ierrou,*)'ec2 elements:',elc2(1:2),elc2(3:5)*degrad
             numerr=numerr+1
          endif
          GOTO 1300
       ENDIF
!          if(verb_moid.ge.20) then
!             WRITE(*,*)'morse,weier,warnflag(1)',morse,weier,warnflag(1)
!             WRITE(*,*)'sflag(2:6)',sflag(2:6)
!             WRITE(*,*)'hzflag,hwflag,multfl',hzflag,hwflag,multfl
!          endif
       alpha = alpha + dpig/60.d0
       beta = beta + dpig/50.d0
       if(verb_moid.ge.20) then
          WRITE(ierrou,*)'APPLY ANGULAR SHIFT: count=',count
          WRITE(ierrou,*)'               alpha,beta:',alpha,beta
          numerr=numerr+1
       endif
! restoring flags
       morse = .true.
       weier = .true.
       warnflag(1) = .true.
       warnflag(2) = .true.
       warnflag(3) = .true.
       sflag(1) = .true.
       sflag(2) = .true.
       sflag(3) = .true.
       sflag(4) = .true.
       sflag(5) = .true.
       sflag(6) = .true.
       hzflag = .true.
       multfl = .true.
       hevalflag = .true.
       GOTO 10 
    ELSE
       IF(verb_moid.ge.20) THEN
          WRITE(ierrou,*)'computation OK! count = ',count
          WRITE(ierrou,*)'morse,weier',morse,weier
          WRITE(ierrou,*)'warnflag(1)',warnflag(1)
          WRITE(ierrou,*)'sflag(2:6)',sflag(2:6)
          WRITE(ierrou,*)'hzflag',hzflag
          numerr=numerr+1
          IF(count.gt.0)THEN
             WRITE(ierrou,*)'Morse and Weierstrass control succeded!'
             WRITE(ierrou,*)'nroots,nummin,nummax',nroots,nummin,nummax
             WRITE(ierrou,*)'ec1 elements:',elc1(1:2),elc1(3:5)*degrad
             WRITE(ierrou,*)'ec2 elements:',elc2(1:2),elc2(3:5)*degrad
          ENDIF
       ENDIF
! continue the computation
    ENDIF

1300 CONTINUE ! exit point
 
  END SUBROUTINE crit_pts

! =====================================================
  SUBROUTINE orbitcoe(alpha,beta,elpl,elcom) 
! HINT: angles must be passed in radians 
    IMPLICIT NONE 
! shift angles
    REAL(KIND=8),INTENT(IN) :: alpha,beta
! cometary elements of the planet and of the comet
    REAL(KIND=8),INTENT(IN),DIMENSION(5) :: elpl
    REAL(KIND=8),INTENT(IN),DIMENSION(5) :: elcom 
! ---------------- end interface ---------------------
    REAL(KIND=8) :: Qpl,i1,Omnod1,omega1 ! planet cometary elements
    REAL(KIND=8) :: q,i2,Omnod2,omega2 ! comet cometary elements
    REAL(KIND=8) :: deltaOm ! diff of node longitudes
! elements in the first 2 rows of the cracovians (see Sitarski 1968)
    REAL(KIND=8) :: sa,ca,sb,cb
! ==================================================================

    Qpl    = elpl(1)
    Epl    = elpl(2)
    i1     = elpl(3)
    Omnod1 = elpl(4)
    omega1 = elpl(5)
    q      = elcom(1)
    ecom   = elcom(2)
    i2     = elcom(3)
    Omnod2 = elcom(4)
    omega2 = elcom(5)
    deltaOm = Omnod2-Omnod1

    PPx = cos(omega1) 
    PPy = sin(omega1)*cos(i1)
    PPz = sin(omega1)*sin(i1)

    QQx = -sin(omega1) 
    QQy = cos(omega1)*cos(i1)
    QQz = cos(omega1)*sin(i1)

    px = cos(omega2)*cos(deltaOm) - &
         & sin(omega2)*cos(i2)*sin(deltaOm)
    py = cos(omega2)*sin(deltaOm) + &
         & sin(omega2)*cos(i2)*cos(deltaOm)
    pz = sin(omega2)*sin(i2)

    qx = -sin(omega2)*cos(deltaOm) - &
         & cos(omega2)*cos(i2)*sin(deltaOm)
    qy = -sin(omega2)*sin(deltaOm) + &
         & cos(omega2)*cos(i2)*cos(deltaOm)
    qz = cos(omega2)*sin(i2)

    KK = PPx*px + PPy*py + PPz*pz
    LL = QQx*px + QQy*py + QQz*pz
    MM = PPx*qx + PPy*qy + PPz*qz
    NN = QQx*qx + QQy*qy + QQz*qz
    Ppl = Qpl*(1.d0+Epl)
    pcom = q*(1.d0+ecom)

! -------- SHIFT COEFFICIENTS ---------
    ca = cos(alpha)
    sa = sin(alpha)
    cb = cos(beta)
    sb = sin(beta)
    A1x = ca*PPx + sa*QQx
    A1y = ca*PPy + sa*QQy
    A1z = ca*PPz + sa*QQz

    B1x = -sa*PPx + ca*QQx
    B1y = -sa*PPy + ca*QQy
    B1z = -sa*PPz + ca*QQz

    a2x = cb*px + sb*qx
    a2y = cb*py + sb*qy
    a2z = cb*pz + sb*qz

    b2x = -sb*px + cb*qx
    b2y = -sb*py + cb*qy
    b2z = -sb*pz + cb*qz

!    A1Q1 = A1x*QQx + A1y*QQy + A1z*QQz
!    B1Q1 = B1x*QQx + B1y*QQy + B1z*QQz
    a2Q1 = a2x*QQx + a2y*QQy + a2z*QQz
    b2Q1 = b2x*QQx + b2y*QQy + b2z*QQz
    B1a2 = B1x*a2x + B1y*a2y + B1z*a2z
    B1b2 = B1x*b2x + B1y*b2y + B1z*b2z
    A1a2 = A1x*a2x + A1y*a2y + A1z*a2z 
    A1b2 = A1x*b2x + A1y*b2y + A1z*b2z 
!    a2q2 = a2x*qx + a2y*qy + a2z*qz
!    b2q2 = b2x*qx + b2y*qy + b2z*qz
    A1q2 = A1x*qx + A1y*qy + A1z*qz
    B1q2 = B1x*qx + B1y*qy + B1z*qz
    b2A1 = A1b2
    b2B1 = B1b2
    a2A1 = A1a2
    a2B1 = B1a2
  END SUBROUTINE orbitcoe
  
! ******************************************************************
! ***** Given the orbcoe's, get the coeff's of the polynomials *****
! ***** of the modified Sylvester matrix ***************************
! *********** written by GIOVANNI F. GRONCHI (Oct.2004) ************
! ==================================================================
  SUBROUTINE matrixdat(alpha,beta)
    IMPLICIT NONE 
    REAL(KIND=8),INTENT(IN) :: alpha,beta ! shift angles
! ========== end interface =========================================
    REAL(KIND=8) :: AE,BE,CE,DE,EE,FE ! coeff of the linear combination 
                                      ! of the rows of Sylvester's matrix
    INTEGER :: i,j,l ! loop indexes   
    REAL(KIND=8) :: sa,ca,sb,cb ! auxiliary
! ==================================================================
 
    sa=sin(alpha)
    ca=cos(alpha)
    sb=sin(beta)
    cb=cos(beta)

    AE = Epl*sa/(1.d0+Epl*ca)
    BE = Epl*sa/(1.d0-Epl*ca)
    CE = (Epl**2-1.d0+Epl**2*sa**2)/((1.d0+Epl*ca)**2)
    DE = (Epl**2-1.d0+Epl**2*sa**2)/((1.d0-Epl*ca)**2)
    EE = (Epl*sa/((1.d0+Epl*ca)**3))*(3.d0*(Epl**2-1.d0) + Epl**2*sa**2)
    FE = (Epl*sa/((1.d0-Epl*ca)**3))*(3.d0*(Epl**2-1.d0) + Epl**2*sa**2)

! p0(t) = p0(5)*t^4 + p0(4)*t^3 + p0(3)*t^2 + p0(2)*t + p0(1)
    p0(5) = -Ppl*A1b2+ecom*pcom*sb-Ppl*ecom**2*cb*A1q2+ &
         &ecom*pcom*sb*Epl*ca+Ppl*ecom*cb*A1b2+Ppl*ecom*A1q2
    p0(4) = -2*Ppl*A1a2+2*Ppl*ecom*cb*A1a2-2*Ppl*ecom**2*sb*A1q2- &
         & 2*ecom*pcom*cb*Epl*ca-2*ecom*pcom*cb+2*Ppl*ecom*sb*A1b2
    p0(3) = 2*Ppl*ecom*(-cb*A1b2+2*sb*A1a2+A1q2)
    p0(2) = -2*Ppl*A1a2-2*Ppl*ecom*cb*A1a2-2*Ppl*ecom**2*sb*A1q2- &
         & 2*ecom*pcom*cb-2*ecom*pcom*cb*Epl*ca-2*Ppl*ecom*sb*A1b2
    p0(1) = -ecom*pcom*sb+Ppl*A1b2+Ppl*ecom*A1q2+Ppl*ecom**2*cb*A1q2+ &
         & Ppl*ecom*cb*A1b2-ecom*pcom*sb*Epl*ca

! p1(t) = p1(5)*t^4 + p1(4)*t^3 + p1(3)*t^2 + p1(2)*t + p1(1)
    p1(5) = -2*ecom*pcom*sb*Epl*sa-2*Ppl*ecom**2*cb*B1q2+ &
         & 2*Ppl*ecom*cb*B1b2+2*Ppl*ecom*B1q2-2*Ppl*B1b2
    p1(4) = 4*ecom*pcom*cb*Epl*sa- &
         & 4*Ppl*ecom**2*sb*B1q2+4*Ppl*ecom*sb*B1b2+ &
         & 4*Ppl*ecom*cb*a2B1-4*Ppl*a2B1
    p1(3) = 4*Ppl*ecom*(2*sb*a2B1-cb*B1b2+B1q2)
    p1(2) = -4*Ppl*a2B1+4*ecom*pcom*cb*Epl*sa- &
         & 4*Ppl*ecom**2*sb*B1q2-4*Ppl*ecom*sb*B1b2-4*Ppl*ecom*cb*a2B1
    p1(1) = 2*Ppl*B1b2+2*Ppl*ecom*B1q2+2*Ppl*ecom**2*cb*B1q2+ &
         & 2*Ppl*ecom*cb*B1b2+2*ecom*pcom*sb*Epl*sa

! p2(t) = p2(5)*t^4 + p2(4)*t^3 + p2(3)*t^2 + p2(2)*t + p2(1)
    p2(5) = Ppl*A1b2+ecom*pcom*sb+Ppl*ecom**2*cb*A1q2- &
         & ecom*pcom*sb*Epl*ca-Ppl*ecom*cb*A1b2-Ppl*ecom*A1q2
    p2(4) = 2*Ppl*A1a2-2*Ppl*ecom*cb*A1a2+2*Ppl*ecom**2*sb*A1q2 &
         & +2*ecom*pcom*cb*Epl*ca-2*ecom*pcom*cb-2*Ppl*ecom*sb*A1b2
    p2(3) = -2*Ppl*ecom*(-cb*A1b2+2*sb*A1a2+A1q2)
    p2(2) = 2*Ppl*A1a2+2*Ppl*ecom*cb*A1a2+2*Ppl*ecom**2*sb*A1q2 &
         & -2*ecom*pcom*cb+2*ecom*pcom*cb*Epl*ca+2*Ppl*ecom*sb*A1b2
    p2(1) = -ecom*pcom*sb-Ppl*A1b2-Ppl*ecom*A1q2-Ppl*ecom**2*cb*A1q2 &
         & -Ppl*ecom*cb*A1b2+ecom*pcom*sb*Epl*ca

! q0(t) = q0(3)*t^2 + q0(2)*t + q0(1)
    q0(3) = pcom*Epl*a2Q1+Epl*Ppl*sa+pcom*Epl*ca*B1a2+ &
         & pcom*Epl**2*ca*a2Q1-Epl*Ppl*sa*ecom*cb+pcom*B1a2
    q0(2) = -2*pcom*Epl*b2Q1-2*pcom*Epl*ca*B1b2-2*pcom*Epl**2*ca*b2Q1 &
         & -2*Epl*Ppl*sa*ecom*sb-2*pcom*B1b2
    q0(1) = Epl*Ppl*sa-pcom*B1a2-pcom*Epl*ca*B1a2-pcom*Epl**2*ca*a2Q1 &
           & +Epl*Ppl*sa*ecom*cb-pcom*Epl*a2Q1
            
! q1(t) = q1(3)*t^2 + q1(2)*t + q1(1)
   q1(3) = 2*Epl*Ppl*ca-2*pcom*Epl**2*sa*a2Q1-2*Epl*Ppl*ca*ecom*cb- &
        &2*pcom*Epl*ca*A1a2-2*pcom*Epl*sa*B1a2-2*pcom*A1a2 
   q1(2) = 4*pcom*A1b2-4*Epl*Ppl*ca*ecom*sb+4*pcom*Epl**2*sa*b2Q1+ &
        & 4*pcom*Epl*sa*B1b2+4*pcom*Epl*ca*A1b2
   q1(1) = 2*Epl*Ppl*ca+2*pcom*Epl**2*sa*a2Q1+2*pcom*Epl*ca*A1a2+ &
        & 2*pcom*Epl*sa*B1a2+2*Epl*Ppl*ca*ecom*cb+2*pcom*A1a2
                     
! q2(t) = q2(3)*t^2 + q2(2)*t + q2(1)
   q2(3) = 2*pcom*Epl*(a2Q1-ca*B1a2+2*sa*A1a2)      
   q2(2) = -4*pcom*Epl*(b2Q1+2*sa*A1b2-ca*B1b2)
   q2(1) = -2*pcom*Epl*(a2Q1-ca*B1a2+2*sa*A1a2)
                      
! q3(t) = q3(3)*t^2 + q3(2)*t + q3(1)
   q3(3) = -2*pcom*A1a2-2*Epl*Ppl*ca*ecom*cb+2*Epl*Ppl*ca+2*pcom*Epl*sa*B1a2+ &
        & 2*pcom*Epl*ca*A1a2-2*pcom*Epl**2*sa*a2Q1 
   q3(2) = 4*pcom*A1b2+4*pcom*Epl**2*sa*b2Q1-4*Epl*Ppl*ca*ecom*sb- &
        & 4*pcom*Epl*ca*A1b2-4*pcom*Epl*sa*B1b2
   q3(1) = 2*Epl*Ppl*ca+2*pcom*A1a2+2*Epl*Ppl*ca*ecom*cb-2*pcom*Epl*sa*B1a2 &
        & -2*pcom*Epl*ca*A1a2+2*pcom*Epl**2*sa*a2Q1

! q4(t) = q4(3)*t^2 + q4(2)*t + q4(1)
   q4(3) = pcom*Epl*a2Q1-pcom*B1a2-Epl*Ppl*sa+Epl*Ppl*sa*ecom*cb- &
        & pcom*Epl**2*ca*a2Q1+pcom*Epl*ca*B1a2
   q4(2) = 2*pcom*B1b2-2*pcom*Epl*b2Q1+2*pcom*Epl**2*ca*b2Q1+ &
        & 2*Epl*Ppl*sa*ecom*sb-2*pcom*Epl*ca*B1b2
   q4(1) = pcom*B1a2-pcom*Epl*a2Q1-Epl*Ppl*sa+pcom*Epl**2*ca*a2Q1 &
        & -Epl*Ppl*sa*ecom*cb-pcom*Epl*ca*B1a2

! ===============================================================
   
   r31(3) = -2*Ppl*(ecom*A1q2-Epl**2*ecom*A1q2+ecom*A1q2*Epl**2*ca**2- &
    & ecom*A1q2*Epl*ca-Epl**2*sa*ecom*B1q2*ca+Epl*sa*ecom*B1q2+A1b2*Epl*ca- &
    & Epl*sa*B1b2-A1b2*Epl**2*ca**2+Epl**2*sa*B1b2*ca-A1b2+Epl**2*A1b2)/(-1+Epl*ca)**2

   r31(2) = -4*Ppl*(-A1a2+Epl**2*sa*a2B1*ca-A1a2*Epl**2*ca**2- &
     & Epl*sa*a2B1+A1a2*Epl*ca+Epl**2*A1a2)/(-1+Epl*ca)**2

   r31(1) = -2*Ppl*(A1b2*Epl**2*ca**2+ecom*A1q2-Epl**2*A1b2-A1b2*Epl*ca- &
     & Epl**2*sa*ecom*B1q2*ca+ecom*A1q2*Epl**2*ca**2+Epl*sa*ecom*B1q2+Epl*sa*B1b2- &
     & Epl**2*sa*B1b2*ca-ecom*A1q2*Epl*ca+A1b2-Epl**2*ecom*A1q2)/(-1+Epl*ca)**2

   r32(3) = -2*Ppl*(-ecom*B1q2+ecom*B1q2*Epl**2*ca**2+Epl**2*sa*ecom*A1q2*ca- &
     & Epl**2*sa*A1b2*ca-B1b2*Epl**2*ca**2+B1b2)/(-1+Epl*ca)/(Epl*ca+1)

   r32(2) = 4*Ppl*(Epl**2*sa*A1a2*ca-a2B1+a2B1*Epl**2*ca**2)/(-1+Epl*ca)/(Epl*ca+1)
   
   r32(1) = -2*Ppl*(-ecom*B1q2+ecom*B1q2*Epl**2*ca**2+Epl**2*sa*ecom*A1q2*ca+ &
     & Epl**2*sa*A1b2*ca-B1b2+B1b2*Epl**2*ca**2)/(-1+Epl*ca)/(Epl*ca+1)
   
   r33(3) = 2*Ppl*(ecom*A1q2-Epl**2*ecom*A1q2+ecom*A1q2*Epl**2*ca**2+ &
     & ecom*A1q2*Epl*ca-Epl**2*sa*ecom*B1q2*ca-Epl*sa*ecom*B1q2-A1b2*Epl*ca+ &
     & Epl*sa*B1b2-A1b2*Epl**2*ca**2+Epl**2*sa*B1b2*ca-A1b2+Epl**2*A1b2)/(Epl*ca+1)**2

   r33(2) = 4*Ppl*(-A1a2+Epl**2*sa*a2B1*ca-A1a2*Epl**2*ca**2+Epl*sa*a2B1- &
     & A1a2*Epl*ca+ Epl**2*A1a2)/(Epl*ca+1)**2
   
   r33(1) = 2*Ppl*(A1b2*Epl**2*ca**2+ecom*A1q2-Epl**2*A1b2+A1b2*Epl*ca- &
     & Epl**2*sa*ecom*B1q2*ca+ecom*A1q2*Epl**2*ca**2-Epl*sa*ecom*B1q2-Epl*sa*B1b2- &
     & Epl**2*sa*B1b2*ca+ecom*A1q2*Epl*ca+A1b2-Epl**2*ecom*A1q2)/(Epl*ca+1)**2

   r34(3) = 2*Ppl*(-B1b2*Epl**2*ca**2-Epl**2*sa*A1b2*ca+Epl**3*sa*ecom*A1q2*ca**2- &
     & 2*Epl**2*ecom*B1q2+Epl**3*ecom*B1q2*ca**3-2*Epl**3*sa*ecom*A1q2+ &
     & 2*Epl*sa*ecom*A1q2- 2*Epl**3*ecom*B1q2*ca+ecom*B1q2*Epl*ca+ &
     & 2*Epl**2*B1b2+ecom*B1q2-B1b2- Epl**3*sa*A1b2*ca**2-2*Epl*sa*A1b2+ &
     & 2*Epl**3*B1b2*ca-B1b2*Epl*ca-Epl**3*B1b2*ca**3+ &
     & 2*Epl**3*sa*A1b2+Epl**2*sa*ecom*A1q2*ca+ecom*B1q2*Epl**2*ca**2)/(Epl*ca+1)**3

   r34(2) = -4*Ppl*(-2*Epl**2*a2B1-2*Epl**3*a2B1*ca+a2B1*Epl*ca+2*Epl*sa*A1a2- &
     & 2*Epl**3*sa*A1a2+Epl**3*a2B1*ca**3+Epl**3*sa*A1a2*ca**2+a2B1+ &
     & a2B1*Epl**2*ca**2+Epl**2*sa*A1a2*ca)/(Epl*ca+1)**3

   r34(1) = 2*Ppl*(-2*Epl**2*B1b2-2*Epl**2*ecom*B1q2+Epl**3*sa*A1b2*ca**2- &
     & 2*Epl**3*ecom*B1q2*ca+ecom*B1q2*Epl*ca+2*Epl*sa*ecom*A1q2- &
     & 2*Epl**3*sa*ecom*A1q2+B1b2*Epl*ca+Epl**3*B1b2*ca**3+Epl**3*ecom*B1q2*ca**3- &
     & 2*Epl**3*B1b2*ca+2*Epl*sa*A1b2-2*Epl**3*sa*A1b2+Epl**3*sa*ecom*A1q2*ca**2+ &
     & ecom*B1q2+B1b2*Epl**2*ca**2+Epl**2*sa*A1b2*ca+ecom*B1q2*Epl**2*ca**2+ &
     & Epl**2*sa*ecom*A1q2*ca+B1b2)/(Epl*ca+1)**3
   
   r35(1) = 8*Epl**2*sa*Ppl*ca*(Epl-1)*(Epl+1)/(-1+Epl*ca)**2/(Epl*ca+1)**2
   
   r36(1) = -4*Epl**2*Ppl*(Epl-1)*(Epl+1)*(Epl*ca-1+2*ca**2)/(-1+Epl*ca)/(Epl*ca+1)**3
   
   r41(3) = 2*Ppl*(B1b2*Epl**2*ca**2+Epl**2*sa*A1b2*ca+Epl**3*sa*ecom*A1q2*ca**2+ &
     & 2*Epl**2*ecom*B1q2+Epl**3*ecom*B1q2*ca**3-2*Epl**3*sa*ecom*A1q2+ &
     & 2*Epl*sa*ecom*A1q2-2*Epl**3*ecom*B1q2*ca+ecom*B1q2*Epl*ca-2*Epl**2*B1b2- &
     & ecom*B1q2+B1b2-Epl**3*sa*A1b2*ca**2-2*Epl*sa*A1b2+2*Epl**3*B1b2*ca- &
     & B1b2*Epl*ca- Epl**3*B1b2*ca**3+2*Epl**3*sa*A1b2-Epl**2*sa*ecom*A1q2*ca- &
     & ecom*B1q2*Epl**2*ca**2)/(-1+Epl*ca)**3

   r41(2) = -4*Ppl*(2*Epl**2*a2B1-2*Epl**3*a2B1*ca+a2B1*Epl*ca+2*Epl*sa*A1a2- &
     & 2*Epl**3*sa*A1a2+Epl**3*a2B1*ca**3+Epl**3*sa*A1a2*ca**2-a2B1- &
     & a2B1*Epl**2*ca**2-Epl**2*sa*A1a2*ca)/(-1+Epl*ca)**3

   r41(1) = 2*Ppl*(2*Epl**2*B1b2+2*Epl**2*ecom*B1q2+Epl**3*sa*A1b2*ca**2- &
     & 2*Epl**3*ecom*B1q2*ca+ecom*B1q2*Epl*ca+2*Epl*sa*ecom*A1q2- &
     & 2*Epl**3*sa*ecom*A1q2 +B1b2*Epl*ca+Epl**3*B1b2*ca**3+ &
     & Epl**3*ecom*B1q2*ca**3-2*Epl**3*B1b2*ca+2*Epl*sa*A1b2 &
     & -2*Epl**3*sa*A1b2+Epl**3*sa*ecom*A1q2*ca**2-ecom*B1q2-B1b2*Epl**2*ca**2- &
     & Epl**2*sa*A1b2*ca-ecom*B1q2*Epl**2*ca**2-Epl**2*sa*ecom*A1q2*ca- &
     & B1b2)/(-1+Epl*ca)**3

   r45(1) = 4*Epl**2*Ppl*(Epl-1)*(Epl+1)*(-Epl*ca+2*ca**2-1)/(-1+Epl*ca)**3/(Epl*ca+1)
     
!   r46(1) = 8*Epl**2*sa*Ppl*ca*(Epl-1)*(Epl+1)/(-1+Epl*ca)**2/(Epl*ca+1)**2
   r46(1) = r35(1)

! set to zero the coefficients of higher degree terms
    p0(6:Nev) = 0.d0 
    p1(6:Nev) = 0.d0 
    p2(6:Nev) = 0.d0
 
    q0(4:Nev) = 0.d0 
    q1(4:Nev) = 0.d0 
    q2(4:Nev) = 0.d0 
    q3(4:Nev) = 0.d0 
    q4(4:Nev) = 0.d0 
    r31(4:Nev) = 0.d0
    r32(4:Nev) = 0.d0
    r33(4:Nev) = 0.d0
    r34(4:Nev) = 0.d0
    r41(4:Nev) = 0.d0

    r35(2:Nev) = 0.d0
    r36(2:Nev) = 0.d0
    r45(2:Nev) = 0.d0
    r46(2:Nev) = 0.d0

  END SUBROUTINE matrixdat

! ******************************************************************
! ***** COMPUTE COEFFICIENTS OF MODIFIED SYLVESTER MATRIX **********
! ******************************************************************
! *********** written by GIOVANNI F. GRONCHI (Oct.2004) ************
! ==================================================================
  SUBROUTINE compmodsylv16(pp0,pp1,pp2,qq0,qq1,qq2,qq3,qq4,rr31,rr32,rr33, &
         & rr34,rr35,rr36,rr41,rr45,rr46,SSYLV) 
    IMPLICIT NONE 
! Sylvester matrix elements                                         
    COMPLEX*16,INTENT(IN) :: pp0,pp1,pp2
    COMPLEX*16,INTENT(IN) :: qq0,qq1,qq2,qq3,qq4
    COMPLEX*16,INTENT(IN) :: rr31,rr32,rr33,rr34,rr35
    COMPLEX*16,INTENT(IN) :: rr36,rr41,rr45,rr46
    COMPLEX*16,INTENT(OUT) :: SSYLV(6,6) 
! ==================================================================

    SSYLV(1,1) = pp2 
    SSYLV(1,2) = (0.d0,0.d0) 
    SSYLV(1,3) = (0.d0,0.d0) 
    SSYLV(1,4) = (0.d0,0.d0) 
    SSYLV(1,5) = qq4 
    SSYLV(1,6) = (0.d0,0.d0) 
!                                                                       
    SSYLV(2,1) = pp1 
    SSYLV(2,2) = pp2 
    SSYLV(2,3) = (0.d0,0.d0) 
    SSYLV(2,4) = (0.d0,0.d0) 
    SSYLV(2,5) = qq3
    SSYLV(2,6) = qq4
!                                                                       
    SSYLV(3,1) = rr31
    SSYLV(3,2) = rr32
    SSYLV(3,3) = rr33 
    SSYLV(3,4) = rr34
    SSYLV(3,5) = rr35
    SSYLV(3,6) = rr36
!                                                                       
    SSYLV(4,1) = rr41
    SSYLV(4,2) = rr31
    SSYLV(4,3) = rr32
    SSYLV(4,4) = rr33
    SSYLV(4,5) = rr45
    SSYLV(4,6) = rr46
!                                                                       
    SSYLV(5,1) = (0.d0,0.d0) 
    SSYLV(5,2) = (0.d0,0.d0) 
    SSYLV(5,3) = pp0 
    SSYLV(5,4) = pp1 
    SSYLV(5,5) = qq0 
    SSYLV(5,6) = qq1 
!                                                                       
    SSYLV(6,1) = (0.d0,0.d0) 
    SSYLV(6,2) = (0.d0,0.d0) 
    SSYLV(6,3) = (0.d0,0.d0) 
    SSYLV(6,4) = pp0 
    SSYLV(6,5) = (0.d0,0.d0) 
    SSYLV(6,6) = qq0 
    
  END SUBROUTINE compmodsylv16

! ******************************************************************
! GIVEN t COMPUTE s SUCH THAT (s,t) SATISFIES THE SYSTEM p = q = 0
! ************ written by GIOVANNI F. GRONCHI (Sept.2004) **********
! ==================================================================
  SUBROUTINE solvesystem(Vpltil,vcomtil,nroots,wzero,zzero,sflag,hwflag) 
    IMPLICIT NONE 
    REAL(KIND=8),INTENT(IN) :: Vpltil,vcomtil
    INTEGER,INTENT(INOUT) :: nroots 
    REAL(KIND=8),INTENT(INOUT), DIMENSION(nroots) :: wzero 
! is inout to eliminate components if the corresponding zzero is complex
    REAL(KIND=8),INTENT(OUT), DIMENSION(nroots) :: zzero 
    LOGICAL,INTENT(INOUT) :: sflag(6) ! solving system messages:
!         sflag(1) = true    OK
!         sflag(1) = false   there are two good solutions             
!         sflag(2) = true    OK                                               
!         sflag(2) = false   neither of the two evaluation is close to 0
!         sflag(3) = true    OK(the discriminant of the 2nd-deg poly is non-negative)
!         sflag(3) = false   the z-component of the solution is complex
!         sflag(4) = true    OK    
!         sflag(4) = false   leading coeff. of the 2nd degree pol. is small
!         sflag(5) = true    OK    
!         sflag(5) = false   1st and 2nd coeffs except of 2nd degree poly
!                            are small                       
!         sflag(6) = true    OK    
!         sflag(6) = false   2nd degree poly have all coefficients small
    LOGICAL,INTENT(INOUT) :: hwflag ! hwflag = .true.  OK!
                                    ! hwflag = .false. abs(root)>10^5
!   2nd degree poly smallness parameter
    REAL(KIND=8),PARAMETER :: eps_2deg=1.d-25
!   2nd degree poly bigness parameter
    REAL(KIND=8),PARAMETER :: big_root=1.d25
!   evaluation smallness parameter
!    REAL(KIND=8),PARAMETER :: eps_eval=1.d-3  
    REAL(KIND=8),PARAMETER :: eps_eval=1.d-5 
!   --------- end interface ------------------------------------------
    REAL(KIND=8) :: svcomtil,cvcomtil,sVpltil,cVpltil
    REAL(KIND=8) :: sx,cx,sy,cy
    REAL(KIND=8) :: stmp(2)
    REAL(KIND=8) :: sf,cf,vcom,sfpl(2),cfpl(2),Vpl(2)
    REAL(KIND=8) :: sfplchk(Nev),cfplchk(Nev)
    REAL(KIND=8) :: sfplchktmp(Nev,Nev),cfplchktmp(Nev,Nev)
    REAL(KIND=8) :: sfchk(Nev),cfchk(Nev)
    REAL(KIND=8) :: vvcom(Nev),VVpl
    REAL(KIND=8) :: deg2z,deg1z,deg0z
!    REAL(KIND=8) :: deg2zold,deg1zold,deg0wold
    REAL(KIND=8) :: deg2ztmp,deg0wtmp
!                        
    REAL(KIND=8) :: z,w,evalst(2),evalstdeg2w(2)
!   choice flag ( choice = 0 ------> init value
!               ( choice = 1 ------> stmp(1)
!               ( choice = 2 ------> stmp(2)
    INTEGER :: choice
!   auxiliary
    REAL(KIND=8), DIMENSION(20) :: wzerotmp(Nev),zzerotmp(Nev) 
    LOGICAL :: discheck(Nev) ! check the discriminant of the 2nd-deg poly for each root
                           ! discheck(i) = .true. OK
                           ! discheck(i) = .false. negative discriminant 
    INTEGER :: nrealroots
!   loop indexes                                                      
    INTEGER :: i,j,k 
!   ==================================================================

! initialization
    discheck(1:Nev)=.true.
    sfplchktmp(1:Nev,1:Nev)=0.d0
    cfplchktmp(1:Nev,1:Nev)=0.d0
! initialization
    sfplchk(1:Nev)=0.d0
    cfplchk(1:Nev)=0.d0
    sfchk(1:Nev)=0.d0
    cfchk(1:Nev)=0.d0                        
    vvcom(1:Nev)=0.d0
    
    svcomtil = sin(vcomtil)
    cvcomtil = cos(vcomtil)
    sVpltil = sin(Vpltil)
    cVpltil = cos(Vpltil)
    
    DO i = 1,nroots     
       w=wzero(i) 
!       write(*,*)'root number',i,' w=',w
! ==============================================                    
! deg2z(w)*z**2 + deg1z(w)*z + deg0z(w) = 0 
! (p2(w)*z**2 + p1(w)*z + p0(w) = 0)                     
       deg2z = p2(5)*w**4 + p2(4)*w**3 + p2(3)*w**2 + p2(2)*w + p2(1)
       deg1z = p1(5)*w**4 + p1(4)*w**3 + p1(3)*w**2 + p1(2)*w + p1(1)
       deg0z = p0(5)*w**4 + p0(4)*w**3 + p0(3)*w**2 + p0(2)*w + p0(1)
!       write(*,*)'deg2z,deg1z,deg0z',deg2z,deg1z,deg0z

! CHECK SMALLNESS of the COEFFICIENTS
       IF(abs(deg2z).lt.eps_2deg)THEN 
          ! linear equation
          sflag(4) = .false.
          if(verb_moid.ge.20) then
             WRITE(ierrou,*)'small leading coeff of 2nd degree poly',deg2z 
             numerr=numerr+1
          endif
          IF(abs(deg1z).ge.eps_2deg) THEN
             zzero(i) = -deg0z/deg1z
             discheck(i) = .true. ! dummy value
             GOTO 133
          ELSEIF(abs(deg1z).lt.eps_2deg)THEN 
             sflag(5) = .false.
             if(verb_moid.ge.20) then
                WRITE(ierrou,*)'small also the coeff of the linear term',deg1z 
                numerr=numerr+1
             endif
             IF(abs(deg0z).lt.eps_2deg)THEN 
                sflag(6) = .false.
                if(verb_moid.ge.20) then
                   WRITE(ierrou,*)'almost completely degenerate &
                        & 2nd degree poly',deg0z,deg1z,deg2z
                   numerr=numerr+1
                endif
             ENDIF
             ! to skip the computation of this root
             GOTO 124 !actually should compute the roots of the 4th-deg poly...
          ENDIF
       ENDIF

! NORMALIZING COEFFICIENTS (the order is IMPORTANT! deg2z must be the last)
       deg0z = deg0z/deg2z
       deg1z = deg1z/deg2z
       deg2z = deg2z/deg2z
!       write(*,*)'deg2z,deg1z,deg0z after normalization:',deg2z,deg1z,deg0z

! check the discriminant
       IF(deg1z**2 - 4.d0*deg0z*deg2z.lt.0.d0) THEN
          discheck(i) = .false.
          sflag(3) = .false.
          if (verb_moid.ge.20) then
             write(ierrou,*)'negative discriminant',deg1z**2 - 4.d0*deg0z*deg2z
             numerr=numerr+1
          endif
! *****************************************
          GOTO 133 ! skip this root
! *****************************************
       ELSE
          discheck(i) = .true.
          sflag(3) = .true.
       ENDIF
! solutions                                                         
       stmp(1) = (-deg1z + dsqrt(deg1z**2 - 4.d0*deg0z*deg2z))/       &
            &        (2.d0*deg2z)                                              
       stmp(2) = (-deg1z - dsqrt(deg1z**2 - 4.d0*deg0z*deg2z))/       &
            &        (2.d0*deg2z)                                              
       IF((abs(stmp(1)).gt.big_root).or.(abs(stmp(2)).gt.big_root)) THEN
          hwflag = .false.
          if(verb_moid.ge.20) then
             write(ierrou,*)'solvesystem: large roots of 2nd deg poly!',&
                  & abs(stmp(1:2))
             numerr=numerr+1
          endif
          GOTO 124
       ENDIF

!     ======================================================            
!     deg4z(w)*z**4 + deg3z(w)*z**3 + deg2z(w)*z**2 + deg1z(w)*z + deg0z(w) = 0
!     (q4(w)*z**4 + q3(w)*z**3 + q2(w)*z**2 + q1(w)*z + q0(w) = 0)

!     conversion into true anomalies                               
       sx = 2.d0*wzero(i)/(1.d0+(wzero(i))**2) 
       cx = (1.d0-(wzero(i))**2)/(1.d0+(wzero(i))**2) 
       sf = sx*cvcomtil + cx*svcomtil
       cf = cx*cvcomtil - sx*svcomtil
       vcom = datan2(sf,cf) 
       vvcom(i) = vcom ! unshifted true anomaly of the comet   
       sfchk(i) = sx 
       cfchk(i) = cx 
!                                                                       
       DO j = 1,2 
          sy =  2.d0*stmp(j)/(1.d0+(stmp(j))**2) 
          cy = (1.d0-(stmp(j))**2)/(1.d0+(stmp(j))**2) 
          sfpl(j) = sy*cVpltil + cy*sVpltil
          cfpl(j) = cy*cVpltil - sy*sVpltil
          Vpl(j) = datan2(sfpl(j),cfpl(j)) 
          VVpl = Vpl(j) ! unshifted true anomaly of the planet
          sfplchktmp(i,j) =  sy 
          cfplchktmp(i,j) =  cy 

          evalstdeg2w(j) = -Ppl*(1.d0+ecom*cos(vvcom(i)))* &
               & ( -ecom*MM*cos(VVpl) - ecom*NN*sin(VVpl) + &
               & sin(vvcom(i))*KK*cos(VVpl) + sin(vvcom(i))*LL*sin(VVpl) - &
               & cos(vvcom(i))*MM*cos(VVpl) - cos(vvcom(i))*NN*sin(VVpl)) -&
               & ecom*pcom*sin(vvcom(i))*(1.d0+Epl*cos(VVpl))
!
          evalst(j) = -pcom*(1.d0+Epl*cos(VVpl))* &
               & ( Epl*LL*cos(vvcom(i)) + Epl*NN*sin(vvcom(i)) + &
               & cos(VVpl)*LL*cos(vvcom(i)) + cos(VVpl)*NN*sin(vvcom(i)) -&
               & sin(VVpl)*KK*cos(vvcom(i)) - sin(VVpl)*MM*sin(vvcom(i))) +&
               & Epl*Ppl*sin(VVpl)*(1.d0+ecom*cos(vvcom(i)))
       ENDDO
!       write(*,*)'************* root(',i,')=',w
!       WRITE(*,*)'EVALSTdeg2w(1)=',evalstdeg2w(1), &
!            &'  EVALSTdeg2w(2)=',evalstdeg2w(2)
!       WRITE(*,*)'EVALST(1)=',evalst(1),'  EVALST(2)=',evalst(2)
!       write(*,*)'--------------------------------------------'

! -------------------------------------------
!     SELECTING the CORRESPONDING SOLUTION                              
       choice = 0 ! initialization 
!     ************ checking smallness of the evaluations *************  
       IF((abs(evalst(1)).lt.eps_eval).or.(abs(evalst(2)).lt.eps_eval))THEN 
!     choosing the one that gives the minimum value                     
          IF (abs(evalst(1)).le.abs(evalst(2))) THEN 
             choice = 1 
             zzero(i) = stmp(1) 
             sfplchk(i) = sfplchktmp(i,1) 
             cfplchk(i) = cfplchktmp(i,1)  
          ELSEIF (abs(evalst(1)).gt.abs(evalst(2))) THEN 
             choice = 2 
             zzero(i) = stmp(2) 
             sfplchk(i) = sfplchktmp(i,2) 
             cfplchk(i) = cfplchktmp(i,2)  
          ENDIF
!     if they are BOTH CLOSE TO ZERO                                    
          IF((abs(evalst(1)).lt.eps_eval).and.(abs(evalst(2)).lt.eps_eval))THEN
!     **************************                                         
             sflag(1) = .false.
!     **************************                                           
             if(verb_moid.ge.20) then
!     writing on screen                                                 
               WRITE(ierrou,*)'SOLVING SYSTEM WARNING: TWO GOOD SOLUTIONS'
               WRITE(ierrou,*)'EVALST(1)=',evalst(1),'  EVALST(2)=',evalst(2)
               numerr=numerr+1
            endif
         ENDIF
!     ================== DOUBLE CORRESPONDENCE CHECK ===================
!     if sflag(1)-.false., then do this check for all the following values of i
!     -------------------------------------------------------------     
!     CASE 1 (i=1)                
          IF((.not.sflag(1)).and.(i.eq.1)) THEN 
!     definitively choosing previously chosen value of z
!     -------------------------------------------------------------     
!     CASE 2 (i>1)                    
          ELSEIF((.not.sflag(1)).and.(i.gt.1))THEN 
!     at first choosing previously chosen value of z
             DO k = 1,i-1 
                IF( (abs(sfchk(i)-sfchk(k)).lt.1.d-1).and.     &
                     &  (abs(cfchk(i)-cfchk(k)).lt.1.d-1).and.         &
                     &  (abs(sfplchk(i)-sfplchk(k)).lt.1.d-1).and.     &
                     &  (abs(cfplchk(i)-cfplchk(k)).lt.1.d-1) )THEN      
!                   write(*,*)'OK, MODIFYING zzero to',stmp(2)
                   IF(choice.eq.1) THEN 
                      zzero(i) = stmp(2) 
                      sfplchk(i) = sfpl(2) 
                      cfplchk(i) = cfpl(2) 
                   ELSEIF(choice.eq.2) THEN 
                      zzero(i) = stmp(1) 
                      sfplchk(i) = sfpl(1) 
                      cfplchk(i) = cfpl(1) 
                   ELSE
                      WRITE(ierrou,*)'ERROR: choice=',choice
                      numerr=numerr+1
                      STOP
                   ENDIF
                ELSE 
!     select previous value of z                                              
                ENDIF
             ENDDO
          ENDIF
!     neither of the evaluation is close to zero                        
       ELSE 
!     ***********************                                           
          sflag(2) = .false. 
!     ***********************                                           
          if(verb_moid.ge.20) then
             WRITE(ierrou,*)' NEITHER OF THE EVALUATIONS IS CLOSE TO ZERO'   
             WRITE(ierrou,*)' SELECTING THE ONE WITH THE SMALLEST VALUE '    
             WRITE(ierrou,*)'EVALST(1)=',evalst(1),'EVALST(2)=',evalst(2)    
             numerr=numerr+1
!            WRITE(*,*)'THE SMALLEST IS',min(abs(evalst(1)),            
!     *           abs(evalst(2)))                                       
          endif
!     choosing the one that gives the minimum value                     
          IF (abs(evalst(1)).le.abs(evalst(2))) THEN 
             zzero(i) = stmp(1) 
          ELSEIF (abs(evalst(1)).gt.abs(evalst(2))) THEN 
             zzero(i) = stmp(2) 
          ENDIF
       ENDIF
133    CONTINUE
    ENDDO
!     selecting roots with both real components
    nrealroots = 0
    DO i=1,nroots
       IF(discheck(i)) THEN
          nrealroots=nrealroots+1
          wzerotmp(nrealroots)=wzero(i)
          zzerotmp(nrealroots)=zzero(i)
       ENDIF
    ENDDO
! renumbering wzero,zzero,nroots
    DO i = 1,nrealroots
       wzero(i) = wzerotmp(i)
       zzero(i) = zzerotmp(i)
    ENDDO
    nroots = nrealroots

124 CONTINUE ! to skip computations in case of degeneration of 2nd deg poly
             ! or in case of large roots
  END SUBROUTINE solvesystem

! ******************************************************************
! ********** SELECT AMONG MINIMA, MAXIMA and SADLE POINTS **********
! ******************************************************************
! *********** written by GIOVANNI F. GRONCHI (Sept. 2004) **********
! last modified 21/10/2005 GFG
! ==================================================================
  SUBROUTINE hessian(Vpl,vcom,ans) 
!    USE fund_const
    IMPLICIT NONE 
!   input angles Vpl,vcom are in radians
    REAL(KIND=8),INTENT(IN) :: Vpl,vcom
    INTEGER,INTENT(OUT) :: ans !  ans = 1: maximum 
!                                 ans = 0: saddle 
!                                 ans= -1: minimum 
!                                 ans= -2: cannot decide
!   -------- end interface -----------------------------------------
    REAL(KIND=8) :: H11,H12,H21,H22,trH,detH 
    REAL(KIND=8) :: Rpl,rcom
    REAL(KIND=8) :: Xpl,Ypl,xcom,ycom
!   eigenvalues                                                       
    REAL(KIND=8) :: discrim
    REAL(KIND=8) :: lam1,lam2 
    REAL(KIND=8), PARAMETER :: eps=epsilon(1.d0)*100 !tolerance parameter
!   ==================================================================

!   HESSIAN MATRIX (of d^2/2, not of d^2, see Sitarski 1968 Acta Astron.) 
!(we have eliminated some terms that are zero at critical points)
      Rpl = Ppl/(1.d0 + Epl*cos(Vpl))
      rcom = pcom/(1.d0 + ecom*cos(vcom))
      Xpl = Rpl*cos(Vpl)
      Ypl = Rpl*sin(Vpl)
      xcom = rcom*cos(vcom)
      ycom = rcom*sin(vcom)
      H11 = (Rpl/Ppl)*( Epl*Rpl*(Rpl/Ppl)*(Epl*Rpl+Xpl) + &
           & Xpl*(KK*xcom+MM*ycom) + Ypl*(LL*xcom+NN*ycom) )
      H12 = (Rpl*rcom/(Ppl*pcom))*( (ecom*rcom+xcom)*(NN*(Epl*Rpl+Xpl) - &
           & MM*Ypl) - ycom*(LL*(Epl*Rpl+Xpl) - KK*Ypl))
      H21 = H12
      H22 = (rcom/pcom)*( ecom*rcom*(rcom/pcom)*(ecom*rcom+xcom) + &
           & xcom*(KK*Xpl+LL*Ypl) + ycom*(MM*Xpl+NN*Ypl) )
!   TRACE OF THE HESSIAN                                              
    trH = H11 + H22 
!   DETERMINANT OF THE HESSIAN                                        
    detH = H11*H22 - H12*H21 
!   check                                                             
    IF (abs(detH).le.eps) THEN 
       WRITE(ierrou,*)'hess: det(H) is very small',detH
       numerr=numerr+1
       ans = -2
       GOTO 10 
    ENDIF
!   compute eigenvalues                                               
    discrim=trH**2-4.d0*detH
    IF (discrim.ge.0.d0) THEN 
       lam1 = (trH + sqrt(discrim))/2.d0 
       lam2 = (trH - sqrt(discrim))/2.d0 
       IF ((lam1.gt.eps).and.(lam2.gt.eps)) THEN 
          ans = -1 
       ELSEIF ((lam1.lt.-eps).and.(lam2.lt.-eps)) THEN 
          ans = 1 
       ELSEIF ((lam1.gt.eps).and.(lam2.lt.-eps)) THEN 
          ans = 0 
       ELSEIF ((lam1.lt.-eps).and.(lam2.gt.eps)) THEN 
          ans = 0 
       ELSE 
          ans = -2
       ENDIF
    ELSE 
       WRITE(ierrou,*)'negative discriminant: complex eigenvalues!' 
       numerr=numerr+1
       ans = -2
    ENDIF
10  CONTINUE 
  END SUBROUTINE hessian

! ******************************************************************
! ********** SELECT AMONG MINIMA, MAXIMA and SADLE POINTS **********
! **** by EVALUATING D2 in a NEIGHBORHOOD of the CRITICAL POiNT ****
! ******************************************************************
! ************* written by GIOVANNI F. GRONCHI (2003) **************
! ********** Department of Mathematics, UNIVERSITY of PISA *********
! ==================================================================
  SUBROUTINE int_eval(rad,Vpl,vcom,ans)
!    USE fund_const
    IMPLICIT NONE 
!   u,upl are passed in radians                                       
    REAL(KIND=8),INTENT(IN) :: rad,Vpl,vcom 
    INTEGER,INTENT(OUT) :: ans !  ans = 1: maximum 
!                                 ans = 0: saddle 
!                                 ans= -1: minimum 
!                                 ans= -2: cannot decide
! ----------- end interface ---------------------------------------
    INTEGER,PARAMETER :: n = 30 
    REAL(KIND=8) :: fpl,fcom
    REAL(KIND=8) :: vD2
    REAL(KIND=8) :: fD2, rr,eps
    INTEGER :: count ! counter
    INTEGER :: j
! =================================================================
    count = 0
    eps=epsilon(1.d0)
    rr=MAX(eps*10000,rad)
    CALL d2eval(Vpl,vcom,vD2) 
    DO j = 0,n-1
       fpl = Vpl + rr*sin(j*dpig/n)
       fcom = vcom + rr*cos(j*dpig/n)
       CALL d2eval(fpl,fcom,fD2) 
       IF(vD2.lt.fD2) THEN
          count = count - 1
       ELSEIF(vD2.gt.fD2) THEN
          count = count + 1
       ELSE
!          WRITE(ierrou,*)'equal values of D2'
!          numerr=numerr+1
       ENDIF
    ENDDO
    IF(count.eq.-n) THEN
       ans = -1
    ELSEIF(count.eq.n) THEN
       ans = 1
    ELSE
       ans = 0
    ENDIF
  END SUBROUTINE int_eval

! ******************************************************************
! ****** EVALUATE THE SQUARED DISTANCE IN THE POINTS fpl,fcom ******
! ******************************************************************
! *********** written by GIOVANNI F. GRONCHI (Sept.2004) ***********
! ==================================================================
  SUBROUTINE d2eval(Vpl,vcom,D2) 
!    USE fund_const
    IMPLICIT NONE                            
    DOUBLE PRECISION,INTENT(IN) :: Vpl,vcom ! true anomalies (in rad.)
    DOUBLE PRECISION,INTENT(OUT):: D2 ! SQUARED DISTANCE function 
!   ------------- end interface --------------------------------------
    DOUBLE PRECISION :: Rpl,rcom
    DOUBLE PRECISION :: Xpl,Ypl,xcom,ycom
    DOUBLE PRECISION :: Xpl1,Xpl2,Xpl3
    DOUBLE PRECISION :: xcom1,xcom2,xcom3
    DOUBLE PRECISION :: XX,YY,ZZ
!   ==================================================================
    Rpl = Ppl/(1.d0 + Epl*cos(Vpl))
    rcom = pcom/(1.d0 + ecom*cos(vcom))
    Xpl = Rpl*cos(Vpl)
    Ypl = Rpl*sin(Vpl)
    xcom = rcom*cos(vcom)
    ycom = rcom*sin(vcom)
! space coords of the planet orbit
    Xpl1 = Xpl*PPx + Ypl*QQx
    Xpl2 = Xpl*PPy + Ypl*QQy
    Xpl3 = Xpl*PPz + Ypl*QQz
! space coords of the comet orbit
    xcom1 = xcom*px + ycom*qx
    xcom2 = xcom*py + ycom*qy
    xcom3 = xcom*pz + ycom*qz
! difference vector
    XX = Xpl1 - xcom1
    YY = Xpl2 - xcom2
    ZZ = Xpl3 - xcom3
    D2 = XX**2 + YY**2 + ZZ**2
!   check if D2 i positive; otherwise set it to zero                  
    IF (D2.ge.0) THEN 
!   do nothing                                                        
    ELSEIF (D2.lt.0) THEN 
       D2 = 0.d0 
    ENDIF
  END SUBROUTINE d2eval

! ********************************************************
! ************* R M S  of  the  M O I D ******************
! ********************************************************
! ================================================
! to select a sign for d_min
SUBROUTINE sign_dmin(pos1,tau1,pos2,tau2,dsign)
  IMPLICIT NONE
! Cartesian coordinates at minimum and tangent vectors
! to the orbits
  DOUBLE PRECISION,INTENT(IN),DIMENSION(3) :: pos1,tau1
  DOUBLE PRECISION,INTENT(IN),DIMENSION(3) :: pos2,tau2
  DOUBLE PRECISION,INTENT(OUT) :: dsign ! dmin sign
! end interface
  DOUBLE PRECISION,DIMENSION(3) :: deltapos,tau3
  DOUBLE PRECISION :: dotprod
  deltapos=pos1-pos2
  CALL prvec(tau1,tau2,tau3)
  dotprod=DOT_PRODUCT(deltapos,tau3)
!  write(ierrou,*)dotprod
  IF(dotprod.gt.0.d0)THEN
     dsign = 1.d0
  ELSEIF(dotprod.lt.0.d0)THEN
     dsign = -1.d0
  ELSE
     write(*,*)'sign_dmin error!!'
     dsign = 2.d0
!     STOP
  ENDIF
END SUBROUTINE sign_dmin

! ==============================================================
! to select a sign for d_min (old version, it shows someproblems)
SUBROUTINE sign_dmin2(pos1,tau1,pos2,tau2,dsign)
  IMPLICIT NONE
  DOUBLE PRECISION,PARAMETER :: eps_S=100.d0*epsilon(1.d0)
! Cartesian coordinates at minimum and tangent vectors
! to the orbits
  DOUBLE PRECISION,INTENT(IN),DIMENSION(3) :: pos1,tau1
  DOUBLE PRECISION,INTENT(IN),DIMENSION(3) :: pos2,tau2
  DOUBLE PRECISION,INTENT(OUT) :: dsign ! dmin sign
!
  DOUBLE PRECISION :: detT1, detT2, detT3
  DOUBLE PRECISION :: S1,S2,S3
  DOUBLE PRECISION :: Smax

  detT1 = tau1(2)*tau2(3)-tau2(2)*tau1(3)
  detT2 = tau1(3)*tau2(1)-tau2(3)*tau1(1)
  detT3 = tau1(1)*tau2(2)-tau2(1)*tau1(2)

  S1=(pos1(1)-pos2(1))*detT1
  S2=(pos1(2)-pos2(2))*detT2
  S3=(pos1(3)-pos2(3))*detT3

! added 5/2/2012
  Smax = max(abs(S1),abs(S2),abs(S3))
  IF(Smax.gt.eps_S)THEN
     IF(Smax.eq.abs(S1))THEN
        dsign = SIGN(1.d0,S1)
     ELSEIF(Smax.eq.abs(S2))THEN
        dsign = SIGN(1.d0,S2)
     ELSEIF(Smax.eq.abs(S3))THEN
        dsign = SIGN(1.d0,S3)
     ENDIF
  ELSE
     WRITE(ierrou,*)'sign_dmin: WARNING! ALMOST DEGENERATE CASE: &
          & S1,S2,S3',S1,S2,S3
  ENDIF
  RETURN
! -----------------------------------

  IF(abs(S1).gt.eps_S)THEN
     dsign=SIGN(1.d0,S1)
  ELSEIF(abs(S2).gt.eps_S)THEN
     dsign=SIGN(1.d0,S2)
  ELSEIF(abs(S3).gt.eps_S)THEN
     dsign=SIGN(1.d0,S3)
  ELSE
     WRITE(ierrou,*)'sign_dmin: WARNING! ALMOST DEGENERATE CASE: &
          & S1,S2,S3',S1,S2,S3
     numerr=numerr+1
     IF(S1.ne.0.d0)THEN
        dsign=SIGN(1.d0,S1)
     ELSEIF(S2.ne.0.d0)THEN
        dsign=SIGN(1.d0,S2)
     ELSEIF(S3.ne.0.d0)THEN
        dsign=SIGN(1.d0,S3)
     ELSE
        WRITE(ierrou,*)'sign_dmin: ERROR! DEGENERATE CASE: S1,S2,S3',S1,S2,S3
        numerr=numerr+1
     ENDIF
  ENDIF
END SUBROUTINE sign_dmin2

!==========================================================
! COMPUTATION of the MUTUAL ELEMENTS mutI,mutom1,mutom2
! and of the ROTATION MATRIX (with derivatives) from
! inertial to mutual reference frame
! written by G. Tommei 01/02/2007
! (see Gronchi, Ph.D. thesis, pp. 113-114)
!==========================================================
SUBROUTINE mutual_ref(el1,el2,mutI,mutom1,mutom2, &
     & dsindcos,dsinmutIdcos,dcosmutI,rotinmut,drotinmut,dcosom)
!  USE orbit_elements
!  USE fund_const
  IMPLICIT NONE
!=====================INTERFACE=======================
! Orbital elements   
! 1: first mass (usually planet)  2:second mass (usually minor body)
  TYPE(orbit_elem), INTENT(IN) :: el1,el2
! Mutual angles and derivatives
  DOUBLE PRECISION, INTENT(OUT) :: mutI,mutom1,mutom2
  DOUBLE PRECISION,OPTIONAL,DIMENSION(2),INTENT(OUT) :: dsindcos
  DOUBLE PRECISION,OPTIONAL,INTENT(OUT) :: dsinmutIdcos
  DOUBLE PRECISION,OPTIONAL,DIMENSION(3,2),INTENT(OUT) :: dcosmutI
! Matrix of rotation from inertial to mutual reference frame
  DOUBLE PRECISION,OPTIONAL,DIMENSION(3,3),INTENT(OUT) :: rotinmut
! Derivatives of the matrix of rotation from inertial to 
! mutual reference frame
  DOUBLE PRECISION,OPTIONAL,DIMENSION(10,3,3),INTENT(OUT) :: drotinmut
! Derivatives of cosom1 and cosom2 (to compute the derivatives
! of the nodal distances)
  DOUBLE PRECISION,OPTIONAL,DIMENSION(2,3,2),INTENT(OUT) :: dcosom
!====== end interface ============================
  DOUBLE PRECISION :: cosmutI,sinmutI
  DOUBLE PRECISION :: cmutom1,smutom1
  DOUBLE PRECISION :: cmutom2,smutom2
! Cometary orbital elements
  TYPE(orbit_elem) :: com1,com2  
  DOUBLE PRECISION :: i1,bigomega1,smallomega1  
  DOUBLE PRECISION :: i2,bigomega2,smallomega2  
! Direction of the angular momentum 
  DOUBLE PRECISION, DIMENSION(3) :: Enne1,Enne2
! Derivatives
  DOUBLE PRECISION, DIMENSION(3) :: dEnne1di1,dEnne1dbigom1
  DOUBLE PRECISION, DIMENSION(3) :: dEnne2di2,dEnne2dbigom2
! Direction of the ascending mutual node
  DOUBLE PRECISION, DIMENSION(3) :: Anod,Anodver
! Direction of the perihelion    
  DOUBLE PRECISION, DIMENSION(3) :: chi1,chi2
! Derivatives
  DOUBLE PRECISION, DIMENSION(3) :: dchi1di1,dchi2di2
  DOUBLE PRECISION, DIMENSION(3) :: dchi1dom1,dchi1dbigom1 
  DOUBLE PRECISION, DIMENSION(3) :: dchi1dbigom2 
  DOUBLE PRECISION, DIMENSION(3) :: dchi2dom2,dchi2dbigom2
  DOUBLE PRECISION, DIMENSION(3) :: dAnoddi1,dAnoddi2
  DOUBLE PRECISION, DIMENSION(3) :: dAnoddom1,dAnoddom2 
  DOUBLE PRECISION, DIMENSION(3) :: dAnoddbigom1,dAnoddbigom2
  DOUBLE PRECISION, DIMENSION(3) :: dAnodverdi1,dAnodverdi2
  DOUBLE PRECISION, DIMENSION(3) :: dAnodverdom1,dAnodverdom2
  DOUBLE PRECISION, DIMENSION(3) :: dAnodverdbigom1,dAnodverdbigom2
  DOUBLE PRECISION :: dcosom1di1, dcosom1dom1,dcosom1dbigom1
  DOUBLE PRECISION :: dcosom1di2, dcosom1dom2,dcosom1dbigom2
  DOUBLE PRECISION :: dcosom2di1, dcosom2dom1,dcosom2dbigom1
  DOUBLE PRECISION :: dcosom2di2, dcosom2dom2,dcosom2dbigom2
! Cross product of Enne1 and Anod 
  DOUBLE PRECISION, DIMENSION(3) :: Enne1Anod,Enne1Anodver
  DOUBLE PRECISION, DIMENSION(3) :: Anodverchi1,Anodverchi2
  DOUBLE PRECISION :: vsize,modAnod,modEnne1Anod
! Derivatives
  DOUBLE PRECISION, DIMENSION(3) :: dEnne1Anodverdi1,dEnne1Anodverdi2
  DOUBLE PRECISION, DIMENSION(3) :: dEnne1Anodverdbigom1,dEnne1Anodverdbigom2
  DOUBLE PRECISION, DIMENSION(3) :: dEnne1Anodverdom1,dEnne1Anodverdom2
! Derivatives of mutual inclination
  DOUBLE PRECISION :: dcosmutIdi1, dcosmutIdi2
  DOUBLE PRECISION :: dcosmutIdbigom1,dcosmutIdbigom2
  DOUBLE PRECISION :: dcosmutIdom1,dcosmutIdom2
! Derivatives of rotation matrix
  DOUBLE PRECISION, DIMENSION(3,3) :: drotinmutdi1,drotinmutdi2
  DOUBLE PRECISION, DIMENSION(3,3) :: drotinmutdbigom1,drotinmutdbigom2
! Auxiliary quantities
  DOUBLE PRECISION :: cOmn1,sOmn1,cOmn2,sOmn2,comeg1,someg1
  DOUBLE PRECISION :: comeg2,someg2,ci1,si1,ci2,si2
  INTEGER :: fail_flag
  DOUBLE PRECISION, DIMENSION(3) :: ii1,iii1
  DOUBLE PRECISION, DIMENSION(3) :: bo1,boo1 
!==========================================================
! Tranformation to Cometary elements
  CALL coo_cha(el1,'COT',com1,fail_flag)
  CALL coo_cha(el2,'COT',com2,fail_flag)
  i1 = com1%coord(3)
  bigomega1 = com1%coord(4)
  smallomega1 = com1%coord(5)
  i2 = com2%coord(3)
  bigomega2 = com2%coord(4)
  smallomega2 = com2%coord(5)
! Angular variables
  cOmn1=cos(bigomega1)
  sOmn1=sin(bigomega1)
  cOmn2=cos(bigomega2)
  sOmn2=sin(bigomega2)
  comeg1=cos(smallomega1)
  someg1=sin(smallomega1)
  comeg2=cos(smallomega2)
  someg2=sin(smallomega2)
  ci1=cos(i1)
  si1=sin(i1)
  ci2=cos(i2)
  si2=sin(i2)
! N1 and N2
  Enne1(1)=sOmn1*si1
  Enne1(2)=-cOmn1*si1
  Enne1(3)=ci1
  Enne2(1)=sOmn2*si2
  Enne2(2)=-cOmn2*si2
  Enne2(3)=ci2
! Mutual inclination
  cosmutI=DOT_PRODUCT(Enne1,Enne2)
  sinmutI=SQRT(1-cosmutI**2)
  mutI=atan2(sinmutI,cosmutI)
! Computation od Anod
  CALL prvec(Enne1,Enne2,Anod)
  modAnod=vsize(Anod)
  Anodver(1:3)=Anod(1:3)/modAnod
! Computation of N1 x Anod
  CALL prvec(Enne1,Anodver,Enne1Anodver)
! Computation of Chi1 and Chi2
  chi1(1)=cOmn1*comeg1-sOmn1*someg1*ci1
  chi1(2)=sOmn1*comeg1+cOmn1*someg1*ci1
  chi1(3)=someg1*si1
  chi2(1)=cOmn2*comeg2-sOmn2*someg2*ci2
  chi2(2)=sOmn2*comeg2+cOmn2*someg2*ci2
  chi2(3)=someg2*si2
! Mutual perihelion arguments
  cmutom1=DOT_PRODUCT(Anodver,chi1)
  cmutom2=DOT_PRODUCT(Anodver,chi2)
  CALL prvec(Anodver,chi1,Anodverchi1)
  CALL prvec(Anodver,chi2,Anodverchi2)
  smutom1=DOT_PRODUCT(Anodverchi1,Enne1)
  smutom2=DOT_PRODUCT(Anodverchi2,Enne2)
  mutom1 = atan2(smutom1,cmutom1)
  mutom2 = atan2(smutom2,cmutom2)

  IF(PRESENT(dsindcos).and.PRESENT(dsinmutIdcos).and.PRESENT(dcosmutI) &
       & .and.PRESENT(rotinmut).and.PRESENT(drotinmut) &
       & .and.PRESENT(dcosom)) THEN
!----------------------------------------------------------------
! Derivatives
!----------------------------------------------------------------
    dsindcos(1)=-cmutom1/smutom1
    dsindcos(2)=-cmutom1/smutom2
    dEnne1di1(1)=sOmn1*ci1
    dEnne1di1(2)=-cOmn1*ci1
    dEnne1di1(3)=-si1
    dEnne1dbigom1(1)=cOmn1*si1
    dEnne1dbigom1(2)=sOmn1*si1
    dEnne1dbigom1(3)=0.d0
    dEnne2di2(1)=sOmn2*ci2
    dEnne2di2(2)=-cOmn2*ci2
    dEnne2di2(3)=-si2
    dEnne2dbigom2(1)=cOmn2*si2
    dEnne2dbigom2(2)=sOmn2*si2
    dEnne2dbigom2(3)=0.d0
    dsinmutIdcos=-cosmutI/sinmutI
! Derivatives of cosmutI=<N1,N2>
    dcosmutIdi1=DOT_PRODUCT(dEnne1di1,Enne2)
    dcosmutIdi2=DOT_PRODUCT(Enne1,dEnne2di2)
    dcosmutIdom1=0.d0
    dcosmutIdom2=0.d0
    dcosmutIdbigom1=DOT_PRODUCT(dEnne1dbigom1,Enne2)
    dcosmutIdbigom2=DOT_PRODUCT(Enne1,dEnne2dbigom2)
! Output array
    dcosmutI(1,1)=dcosmutIdi1
    dcosmutI(1,2)=dcosmutIdi2
    dcosmutI(2,1)=dcosmutIdom1
    dcosmutI(2,2)=dcosmutIdom2
    dcosmutI(3,1)=dcosmutIdbigom1
    dcosmutI(3,2)=dcosmutIdbigom2
! ---------------------------------------------------------------
! Computation of the derivatives of cos(mutompl)=cosom1=
! =<Anodver,chi1> w.r.t. I1    
! ---------------------------------------------------------------
! Derivative of Chi1
    dchi1di1(1)=sOmn1*someg1*si1
    dchi1di1(2)=-cOmn1*someg1*si1
    dchi1di1(3)=someg1*ci1
! Derivative of Anod
    dAnoddi1(1)=-si1*cOmn2*si2-ci2*cOmn1*ci1
    dAnoddi1(2)=-si1*sOmn2*si2-ci2*sOmn1*ci1
    dAnoddi1(3)=cOmn1*ci1*sOmn2*si2-sOmn1*ci1*cOmn2*si2
! Derivative of Anodver
    dAnodverdi1(1)=(1.d0/modAnod)*dAnoddi1(1)- &
         & Anod(1)*DOT_PRODUCT(Anod,dAnoddi1)/modAnod**3
    dAnodverdi1(2)=(1.d0/modAnod)*dAnoddi1(2)- &
         & Anod(2)*DOT_PRODUCT(Anod,dAnoddi1)/modAnod**3
    dAnodverdi1(3)=(1.d0/modAnod)*dAnoddi1(3)- &
         & Anod(3)*DOT_PRODUCT(Anod,dAnoddi1)/modAnod**3
! Derivative of cosom1
    dcosom1di1=DOT_PRODUCT(dAnodverdi1,chi1)+DOT_PRODUCT(Anodver,dchi1di1)
    dcosom(1,1,1)=dcosom1di1
! --------------------------------------------------------------------
! Computation of the derivatives of cos(mutompl)=cosom1 w.r.t. omega1
! --------------------------------------------------------------------
! Derivative of Chi1
    dchi1dom1(1)=-someg1*cOmn1-comeg1*sOmn1*ci1
    dchi1dom1(2)=-someg1*sOmn1+comeg1*cOmn1*ci1
    dchi1dom1(3)=comeg1*si1
! Derivative of cosom1
    dcosom1dom1=DOT_PRODUCT(Anodver,dchi1dom1)
    dcosom(1,2,1)=dcosom1dom1  
! --------------------------------------------------------------------
! Computation of the derivatives of cos(mutompl)=cosom1 w.r.t. Omega1
! --------------------------------------------------------------------
! Derivative of Chi1
    dchi1dbigom1(1)=-sOmn1*comeg1-cOmn1*someg1*ci1
    dchi1dbigom1(2)=cOmn1*comeg1-sOmn1*someg1*ci1
    dchi1dbigom1(3)=0.d0
! Derivative of Anod    
    dAnoddbigom1(1)=sOmn1*ci2*si1
    dAnoddbigom1(2)=-cOmn1*ci2*si1
    dAnoddbigom1(3)=-sOmn1*si1*sOmn2*si2-cOmn1*si1*cOmn2*si2
! Derivative of Anodver
    dAnodverdbigom1(1)=(1.d0/modAnod)*dAnoddbigom1(1)- &
         & Anod(1)*DOT_PRODUCT(Anod,dAnoddbigom1)/modAnod**3
    dAnodverdbigom1(2)=(1.d0/modAnod)*dAnoddbigom1(2)- &
         & Anod(2)*DOT_PRODUCT(Anod,dAnoddbigom1)/modAnod**3
    dAnodverdbigom1(3)=(1.d0/modAnod)*dAnoddbigom1(3)- &
         & Anod(3)*DOT_PRODUCT(Anod,dAnoddbigom1)/modAnod**3
! Derivative of cosom1
    dcosom1dbigom1=DOT_PRODUCT(dAnodverdbigom1,chi1)+ &
         & DOT_PRODUCT(Anodver,dchi1dbigom1)
     dcosom(1,3,1)=dcosom1dom1
! --------------------------------------------------------------------
! Computation of the derivatives of cos(mutompl)=cosom1 w.r.t. I2
! --------------------------------------------------------------------
! Derivative of Anod 
    dAnoddi2(1)=ci1*cOmn2*ci2+si2*cOmn1*si1
    dAnoddi2(2)=ci1*sOmn2*ci2+si2*sOmn1*si1
    dAnoddi2(3)=cOmn1*si1*sOmn2*ci2-sOmn1*si1*cOmn2*ci2
! Derivative of Anodver
    dAnodverdi2(1)=(1.d0/modAnod)*dAnoddi2(1)- &
         & Anod(1)*DOT_PRODUCT(Anod,dAnoddi2)/modAnod**3
    dAnodverdi2(2)=(1.d0/modAnod)*dAnoddi2(2)- &
         & Anod(2)*DOT_PRODUCT(Anod,dAnoddi2)/modAnod**3
    dAnodverdi2(3)=(1.d0/modAnod)*dAnoddi2(3)- &
         & Anod(3)*DOT_PRODUCT(Anod,dAnoddi2)/modAnod**3
! Derivative of cosom1
    dcosom1di2=DOT_PRODUCT(dAnodverdi2,chi1)
    dcosom(1,1,2)=dcosom1di2
! --------------------------------------------------------------------
! Computation of the derivatives of cos(mutompl)=cosom1 w.r.t. omega2
! --------------------------------------------------------------------
    dcosom1dom2=0.d0
    dcosom(1,2,2)=dcosom1dom2
! --------------------------------------------------------------------
! Computation of the derivatives of cos(mutompl)=cosom1 w.r.t. Omega2
! --------------------------------------------------------------------
! Derivative of Anod
    dAnoddbigom2(1)=-ci1*sOmn2*si2
    dAnoddbigom2(2)=ci1*cOmn2*si2
    dAnoddbigom2(3)=cOmn1*si1*cOmn2*si2+sOmn1*si1*sOmn2*si2
! Derivative of Anodver
    dAnodverdbigom2(1)=(1.d0/modAnod)*dAnoddbigom2(1)- &
         & Anod(1)*DOT_PRODUCT(Anod,dAnoddbigom2)/modAnod**3
    dAnodverdbigom2(2)=(1.d0/modAnod)*dAnoddbigom2(2)- &
         & Anod(2)*DOT_PRODUCT(Anod,dAnoddbigom2)/modAnod**3
    dAnodverdbigom2(3)=(1.d0/modAnod)*dAnoddbigom2(3)- &
         & Anod(3)*DOT_PRODUCT(Anod,dAnoddbigom2)/modAnod**3
! Derivative of cosom1
    dcosom1dbigom2=DOT_PRODUCT(dAnodverdbigom2,chi1)
    dcosom(1,3,2)=dcosom1dbigom2
! --------------------------------------------------------------------
! Computation of the derivatives of cos(mutom)=cosom2=
! <Anod,chi2> w.r.t. I2    
! --------------------------------------------------------------------
! Derivative of Chi2
    dchi2di2(1)=someg2*sOmn2*si2
    dchi2di2(2)=-someg2*cOmn2*si2
    dchi2di2(3)=someg2*ci2
! Derivative of cosom2
    dcosom2di2=DOT_PRODUCT(dAnodverdi2,chi2)+DOT_PRODUCT(Anodver,dchi2di2)
    dcosom(2,1,2)=dcosom2di2 
! --------------------------------------------------------------------
! Computation of the derivatives of cos(mutom)=cosom2 w.r.t. omega2
! --------------------------------------------------------------------
! Derivative of Chi2
    dchi2dom2(1)=-someg2*cOmn2-comeg2*sOmn2*ci2
    dchi2dom2(2)=-someg2*sOmn2+comeg2*cOmn2*ci2
    dchi2dom2(3)=comeg2*si2
! Derivative of cosom2
    dcosom2dom2=DOT_PRODUCT(Anodver,dchi2dom2)
    dcosom(2,2,2)=dcosom2dom2
! --------------------------------------------------------------------
! Computation of the derivatives of cos(mutom)=cosom2 w.r.t. Omega2
! --------------------------------------------------------------------
! Derivative of Chi2
    dchi2dbigom2(1)=-sOmn2*comeg2-cOmn2*someg2*ci2
    dchi2dbigom2(2)=cOmn2*comeg2-sOmn2*someg2*ci2
    dchi2dbigom2(3)=0.d0
! Derivative of cosom2
    dcosom2dbigom2=DOT_PRODUCT(dAnodverdbigom2,chi2)+ &
         & DOT_PRODUCT(Anodver,dchi2dbigom2)
    dcosom(2,3,2)=dcosom2dbigom2
! --------------------------------------------------------------------
! Computation of the derivatives of cos(mutom)=cosom2 w.r.t. I1
! --------------------------------------------------------------------
! Derivative of cosom2
    dcosom2di1=DOT_PRODUCT(dAnodverdi1,chi2)
    dcosom(2,1,1)=dcosom2di1
! --------------------------------------------------------------------
! Computation of the derivatives of cos(mutom)=cosom2 w.r.t. omega1
! --------------------------------------------------------------------
! Derivative of cosom2
    dcosom2dom1=0.d0
    dcosom(2,2,1)=dcosom2dom1
! --------------------------------------------------------------------
! Computation of the derivatives of cos(mutom)=cosom2 w.r.t. Omega1
! --------------------------------------------------------------------
! Derivative of cosom2
    dcosom2dbigom1=DOT_PRODUCT(dAnodverdbigom1,chi2)
    dcosom(2,3,1)=dcosom2dbigom1
! ----------------------------------------------------------------
!------------------------------------------------------------------
! Rotation matrix from mutual to inertial coordinates: rotinmut
! (from inertial to mutual reference frame)
! It does not depend on a1,e1,a2,e2,om1,om2
!------------------------------------------------------------------
    rotinmut(1:3,1)=Anodver
    rotinmut(1:3,2)=Enne1Anodver
    rotinmut(1:3,3)=Enne1
! Derivatives of the second column of matrix rotinmut w.r.t. I1,I2
    CALL prvec(dEnne1di1,Anodver,ii1)
    CALL prvec(Enne1,dAnodverdi1,iii1)
    dEnne1Anodverdi1=ii1+iii1
    CALL prvec(Enne1,dAnodverdi2,dEnne1Anodverdi2)
! Derivatives of the second column of matrix rotinmut w.r.t. bigom1,bigom2
    CALL prvec(dEnne1dbigom1,Anodver,bo1)
    CALL prvec(Enne1,dAnodverdbigom1,boo1)
    dEnne1Anodverdbigom1=bo1+boo1
    CALL prvec(Enne1,dAnodverdbigom2,dEnne1Anodverdbigom2)
! Derivatives of the second column of matrix rotinmut w.r.t. om1,om2
    dEnne1Anodverdom1(1:3)=0.d0
    dEnne1Anodverdom2(1:3)=0.d0
! Derivatives of the matrix rotinmut w.r.t. I1,I2
    drotinmutdi1(1:3,1)=dAnodverdi1
    drotinmutdi1(1:3,2)=dEnne1Anodverdi1
    drotinmutdi1(1:3,3)=dEnne1di1
    drotinmutdi2(1:3,1)=dAnodverdi2
    drotinmutdi2(1:3,2)=dEnne1Anodverdi2
    drotinmutdi2(1:3,3)=0.d0
! Derivatives of the matrix rotinmut w.r.t. bigom1,bigom2
    drotinmutdbigom1(1:3,1)=dAnodverdbigom1
    drotinmutdbigom1(1:3,2)=dEnne1Anodverdbigom1
    drotinmutdbigom1(1:3,3)=dEnne1dbigom1
    drotinmutdbigom2(1:3,1)=dAnodverdbigom2
    drotinmutdbigom2(1:3,2)=dEnne1Anodverdbigom2
    drotinmutdbigom2(1:3,3)=0.d0
! 3-index array of derivatives od rotinmut
    drotinmut(1:10,1:3,1:3)=0.d0      ! Initialization
    drotinmut(3,1:3,1:3)=drotinmutdi1
    drotinmut(4,1:3,1:3)=drotinmutdbigom1
    drotinmut(8,1:3,1:3)=drotinmutdi2
    drotinmut(9,1:3,1:3)=drotinmutdbigom2
  ENDIF
END SUBROUTINE mutual_ref

! =================================================================
SUBROUTINE CP_newton_raphson(com1c,com2c,iof1,iof2)
  DOUBLE PRECISION,INTENT(IN),DIMENSION(5) :: com1c,com2c
  DOUBLE PRECISION,INTENT(INOUT) :: iof1,iof2
! ----- end interface ----
  DOUBLE PRECISION :: gradd2f1,gradd2f2,detH
  DOUBLE PRECISION,DIMENSION(2,2) :: Hess,invH
  DOUBLE PRECISION :: q1,e1,i1,Om1,omeg1,q2,e2,i2,Om2,omeg2
  DOUBLE PRECISION :: f1,f2,diff1,diff2 ! auxiliary

! elements
  q1 = com1c(1) 
  e1 = com1c(2)
  i1 = com1c(3)
  Om1 = com1c(4)
  omeg1 = com1c(5)
  q2 = com2c(1) 
  e2 = com2c(2) 
  i2 = com2c(3)
  Om2 = com2c(4)
  omeg2 = com2c(5)

  f1=iof1
  f2=iof2
  ! GRADIENT of d^2
  gradd2f1                                                          &
     &= 2.D0/(1.D0+e1*cos(f1))*(e1*q1**2*(1.D0+e1)**2.D0/(1.D0+e1*cos   &
     &(f1))**2.D0*sin(f1)+q1*(1.D0+e1)/(1.D0+e1*cos(f1))*sin(f1)*(KK*q2*&
     &(1.D0+e2)/(1.D0+e2*cos(f2))*cos(f2)+MM*q2*(1.D0+e2)/(1.D0+e2*cos(f&
     &2))*sin(f2))-(e1*q1*(1.D0+e1)/(1.D0+e1*cos(f1))+q1*(1.D0+e1)/(1.D0&
     &+e1*cos(f1))*cos(f1))*(LL*q2*(1.D0+e2)/(1.D0+e2*cos(f2))*cos(f2)+N&
     &N*q2*(1.D0+e2)/(1.D0+e2*cos(f2))*sin(f2)))                        
                                                                        
  gradd2f2                                                          &
     &= 2.D0/(1.D0+e2*cos(f2))*(e2*q2**2*(1.D0+e2)**2.D0/(1.D0+e2*cos   &
     &(f2))**2.D0*sin(f2)+q2*(1.D0+e2)/(1.D0+e2*cos(f2))*sin(f2)*(KK*q1*&
     &(1.D0+e1)/(1.D0+e1*cos(f1))*cos(f1)+LL*q1*(1.D0+e1)/(1.D0+e1*cos(f&
     &1))*sin(f1))-(e2*q2*(1.D0+e2)/(1.D0+e2*cos(f2))+q2*(1.D0+e2)/(1.D0&
     &+e2*cos(f2))*cos(f2))*(MM*q1*(1.D0+e1)/(1.D0+e1*cos(f1))*cos(f1)+N&
     &N*q1*(1.D0+e1)/(1.D0+e1*cos(f1))*sin(f1)))                        
! force to zero (critical point)
!      gradd2f1=0.D0
!      gradd2f2=0.D0
! Hessian Matrix H of d^2

  Hess(1,1)&
     &= e1/(1.D0+e1*cos(f1))*sin(f1)*gradd2f1+2.D0/(1.D0+e1*cos(f1))*&
     &(e1**2*q1**2*(1.D0+e1)**2.D0/(1.D0+e1*cos(f1))**3.D0*sin(f1)**2.D0&
     &+(e1*q1*(1.D0+e1)/(1.D0+e1*cos(f1))**2.D0*sin(f1)**2.D0+q1*(1.D0+e&
     &1)/(1.D0+e1*cos(f1))*cos(f1))*(e1*q1*(1.D0+e1)/(1.D0+e1*cos(f1))+K&
     &K*q2*(1.D0+e2)/(1.D0+e2*cos(f2))*cos(f2)+MM*q2*(1.D0+e2)/(1.D0+e2*&
     &cos(f2))*sin(f2))-(e1**2*q1*(1.D0+e1)/(1.D0+e1*cos(f1))**2.D0*sin(&
     &f1)+e1*q1*(1.D0+e1)/(1.D0+e1*cos(f1))**2.D0*cos(f1)*sin(f1)-q1*(1.&
     &D0+e1)/(1.D0+e1*cos(f1))*sin(f1))*(LL*q2*(1.D0+e2)/(1.D0+e2*cos(f2&
     &))*cos(f2)+NN*q2*(1.D0+e2)/(1.D0+e2*cos(f2))*sin(f2)))            
                                                                        
  Hess(2,2)&
     &= e2/(1.D0+e2*cos(f2))*sin(f2)*gradd2f2+2.D0/(1.D0+e2*cos(f2))*&
     &(e2**2*q2**2*(1.D0+e2)**2.D0/(1.D0+e2*cos(f2))**3.D0*sin(f2)**2.D0&
     &+(e2*q2*(1.D0+e2)/(1.D0+e2*cos(f2))**2.D0*sin(f2)**2.D0+q2*(1.D0+e&
     &2)/(1.D0+e2*cos(f2))*cos(f2))*(e2*q2*(1.D0+e2)/(1.D0+e2*cos(f2))+K&
     &K*q1*(1.D0+e1)/(1.D0+e1*cos(f1))*cos(f1)+LL*q1*(1.D0+e1)/(1.D0+e1*&
     &cos(f1))*sin(f1))-(e2**2*q2*(1.D0+e2)/(1.D0+e2*cos(f2))**2.D0*sin(&
     &f2)+e2*q2*(1.D0+e2)/(1.D0+e2*cos(f2))**2.D0*cos(f2)*sin(f2)-q2*(1.&
     &D0+e2)/(1.D0+e2*cos(f2))*sin(f2))*(MM*q1*(1.D0+e1)/(1.D0+e1*cos(f1&
     &))*cos(f1)+NN*q1*(1.D0+e1)/(1.D0+e1*cos(f1))*sin(f1)))            
                                                                        
  Hess(1,2)&
     &= 2.D0/(1.D0+e1*cos(f1))*(q1*(1.D0+e1)/(1.D0+e1*cos(f1))*sin(f1&
     &)*(KK*(e2*q2*(1.D0+e2)/(1.D0+e2*cos(f2))**2.D0*cos(f2)*sin(f2)-q2*&
     &(1.D0+e2)/(1.D0+e2*cos(f2))*sin(f2))+MM*(e2*q2*(1.D0+e2)/(1.D0+e2*&
     &cos(f2))**2.D0*sin(f2)**2.D0+q2*(1.D0+e2)/(1.D0+e2*cos(f2))*cos(f2&
     &)))-(e1*q1*(1.D0+e1)/(1.D0+e1*cos(f1))+q1*(1.D0+e1)/(1.D0+e1*cos(f&
     &1))*cos(f1))*(LL*(e2*q2*(1.D0+e2)/(1.D0+e2*cos(f2))**2.D0*cos(f2)*&
     &sin(f2)-q2*(1.D0+e2)/(1.D0+e2*cos(f2))*sin(f2))+NN*(e2*q2*(1.D0+e2&
     &)/(1.D0+e2*cos(f2))**2.D0*sin(f2)**2.D0+q2*(1.D0+e2)/(1.D0+e2*cos(&
     &f2))*cos(f2))))                                                   
                                                                        
  Hess(2,1)=Hess(1,2)
  detH=Hess(1,1)*Hess(2,2)-Hess(1,2)**2

  invH(1,1)=Hess(2,2)/detH
  invH(2,1)=-Hess(2,1)/detH
  invH(1,2)=-Hess(1,2)/detH
  invH(2,2)=Hess(1,1)/detH
  
  diff1 = invH(1,1)*gradd2f1 + invH(1,2)*gradd2f2
  diff2 = invH(2,1)*gradd2f1 + invH(2,2)*gradd2f2     
  iof1= f1 - diff1
  iof2= f2 - diff2

END SUBROUTINE CP_newton_raphson
! =================================================================
SUBROUTINE dmintil_rms(el1,el2,nummin,dmintil,c1min,c2min,&
     & unc1,unc2,dmintrms,ddmintdel2,comp_flag,detH,detHrms, &
     & sint1t2,taurms,sinmutI,sinmutIrms,chk_der,v1min,v2min,tau_1,tau_2,dminflag)
  IMPLICIT NONE
! ===================== INTERFACE =======================
  TYPE(orbit_elem),INTENT(IN) :: el1,el2 ! orbital elements
                                         ! el1: planet  el2: asteroid
  INTEGER,INTENT(OUT) :: nummin ! number of relative minima found 
  DOUBLE PRECISION,DIMENSION(nminx),INTENT(OUT) :: dmintil ! local minima
                                                            ! with sign
! ================ OPTIONAL VARIABLES ===================
  DOUBLE PRECISION,OPTIONAL,DIMENSION(3,nminx),INTENT(OUT) :: c1min,c2min
  TYPE(orb_uncert),OPTIONAL,INTENT(IN) :: unc1,unc2
  DOUBLE PRECISION,OPTIONAL,DIMENSION(nminx),INTENT(OUT) :: dmintrms
  DOUBLE PRECISION,OPTIONAL,DIMENSION(5,nminx),INTENT(OUT) :: ddmintdel2
! comp_flag(1) --> sinmutI; comp_flag(2) --> detH; comp_flag(3) --> st1t2
  LOGICAL,INTENT(OUT),OPTIONAL,DIMENSION(3,nminx) :: comp_flag ! check 
! detH, detHrms, sint1t2, taurms at local minima
  DOUBLE PRECISION,OPTIONAL,DIMENSION(nminx),INTENT(OUT) :: detH,detHrms 
  DOUBLE PRECISION,OPTIONAL,DIMENSION(nminx),INTENT(OUT) :: sint1t2
  DOUBLE PRECISION,OPTIONAL,DIMENSION(nminx),INTENT(OUT) :: taurms
!
  DOUBLE PRECISION,OPTIONAL,INTENT(OUT) :: sinmutI
  DOUBLE PRECISION,OPTIONAL,INTENT(OUT) :: sinmutIrms
  LOGICAL,OPTIONAL,INTENT(IN) :: chk_der
  DOUBLE PRECISION,OPTIONAL,DIMENSION(3,nminx),INTENT(OUT) :: tau_1,tau_2
  DOUBLE PRECISION,OPTIONAL,DIMENSION(nminx),INTENT(OUT):: v1min,v2min
! failure in computing moid with sign
! dminflag = 1 could not compute dmin
! dminflag = 2 could not give a sign to dmin
  INTEGER,INTENT(OUT),OPTIONAL :: dminflag(nminx) 
! =================== END INTERFACE =====================
  LOGICAL :: chk_der_aux   
  LOGICAL,DIMENSION(3) :: comp_flag_aux
  DOUBLE PRECISION :: detH_aux,detHrms_aux,taurms_aux
  DOUBLE PRECISION :: sinmutI_aux,sinmutIrms_aux
  DOUBLE PRECISION,DIMENSION(nminx) :: sint1t2_aux
! ------------------------------------------------------------
  DOUBLE PRECISION :: dsign ! selected sign
  DOUBLE PRECISION, DIMENSION(1,1) ::dmint_cov
! ------Elements and Jacobian Matrices------------------------   
  DOUBLE PRECISION, DIMENSION(6,6) :: gcom1,gcom2
  INTEGER :: fail_flag
! -----Covariance matrix of equinoctial and cometary elements-----
  DOUBLE PRECISION,DIMENSION(10,10) :: covcom
  DOUBLE PRECISION, DIMENSION(6,6) :: jaccomel1,jaccomel2
! -------------------------------------------
  TYPE(orbit_elem) :: com1,com2 ! auxiliary
  TYPE(orb_uncert) :: unccom1,unccom2
! ----- for crit_pts -----
  INTEGER, PARAMETER :: iplam=3
  REAL(KIND=8):: fpl(poldeg),fcom(poldeg)! true anomalies 
  INTEGER :: nstat,nummax ! number of statpts,minima and maxim
! error/warning flags for crit_pts
  INTEGER :: answer(poldeg) ! type of critical point
  LOGICAL :: morse       ! check with Morse theory 
  LOGICAL :: weier       ! check with Weierstrass theory:
  LOGICAL :: warnflag(3) ! program warning flags:
  LOGICAL :: sflag(6)    ! solving system messages:
  LOGICAL :: hzflag
  LOGICAL :: hwflag
  LOGICAL :: multfl
  LOGICAL :: hevalflag
! -------------------------------------------------------------------
  DOUBLE PRECISION, DIMENSION(poldeg):: D2 ! squared distance function 
  DOUBLE PRECISION, DIMENSION(nminx):: D2min ! squared distance function 
  DOUBLE PRECISION :: DD2
  DOUBLE PRECISION, DIMENSION(poldeg) :: f1,f2
  DOUBLE PRECISION, DIMENSION(nminx) :: f1min,f2min
! for the Hessian matrix
  DOUBLE PRECISION,DIMENSION(2,2) :: Hess
  DOUBLE PRECISION,DIMENSION(poldeg) :: trH
  DOUBLE PRECISION :: discrim 
  DOUBLE PRECISION :: lam1,lam2 ! eigenvalues                      
! ordered 3-uples,for sorting 
  INTEGER, DIMENSION(poldeg) :: srtnum
!
  DOUBLE PRECISION, DIMENSION(poldeg) ::lambda1,lambda2
! to compute derivatives of dmintil
  TYPE(orbit_elem), DIMENSION(nminx) :: el1min,el2min
  TYPE(orbit_elem), DIMENSION(nminx) :: car1min,car2min
  TYPE(orbit_elem) :: com1min,com2min
  DOUBLE PRECISION,DIMENSION(1,10)  :: derdmintil
  DOUBLE PRECISION,DIMENSION(10,1)  :: tderdmintil
  DOUBLE PRECISION,DIMENSION(3)  :: ddeltadcom
  DOUBLE PRECISION :: vsize
! ----------------------------------------------------------------
  DOUBLE PRECISION, DIMENSION(6,6) :: dcardel1,dcardel2
  DOUBLE PRECISION, DIMENSION(6,10) :: dcardeq
  DOUBLE PRECISION, DIMENSION(6,6) :: dcardcom1,dcardcom2
  DOUBLE PRECISION, DIMENSION(6,10) :: dcardcom
  DOUBLE PRECISION, DIMENSION(3) :: tau1,tau2,tau3,tau3hat
!
  INTEGER ::i,j,k,h,kmin   ! loop index 
  DOUBLE PRECISION, DIMENSION(poldeg) :: u1,l1,u2,l2
! ----------------------------------------------------------------------

!  chk_der = .false. ! no control of derivatives

  IF(PRESENT(dminflag))THEN
     dminflag=0 !initialization
  ENDIF

!  IF(verb_moid.ge.20) THEN
!     WRITE(*,*)'Cometary elements:'
!     WRITE(*,*)'COM1:',com1%coord(1:2),com1%coord(3:5)*degrad
!     WRITE(*,*)'COM2:',com2%coord(1:2),com2%coord(3:5)*degrad
!  ENDIF
     
  IF(PRESENT(unc1).and.PRESENT(unc2))THEN
! conversion into cometary elements
     CALL coo_cha(el1,'COT',com1,fail_flag,jaccomel1)
     CALL coo_cha(el2,'COT',com2,fail_flag,jaccomel2)
! convert into uncertainty relative to cometary elems
     CALL convertunc(unc1,jaccomel1,unccom1)
     CALL convertunc(unc2,jaccomel2,unccom2)
! covariance matrix w.r.t. COM elems
     gcom1=unccom1%g(1:6,1:6)
     gcom2=unccom2%g(1:6,1:6)
     covcom(1:10,1:10)=0.d0
     covcom(1:5,1:5)=gcom1(1:5,1:5)
     covcom(6:10,6:10)=gcom2(1:5,1:5)
  ELSE
! conversion into cometary elements 
! ***********************************************************
! HINT: COM(1:5)=COT(1:5); 
! COM(6) and COT(6), which are different, are NOT used here!
! ***********************************************************
!     WRITE(*,*) el2%coo, el1%coo
     IF(el2%coo.eq.'COM')THEN
        com2=el2
     ELSE
        CALL coo_cha(el2,'COT',com2,fail_flag)
     ENDIF
     IF(el1%coo.eq.'COM')THEN
        com1=el1
     ELSE
        CALL coo_cha(el1,'COT',com1,fail_flag)
     ENDIF
  ENDIF

! check for circular coplanar orbits
  IF((com1%coord(2).eq.0.d0).and.(com2%coord(2).eq.0.d0)) THEN
     IF((com1%coord(3).eq.0.d0.and.com2%coord(3).eq.0.d0).or. &
          &(com1%coord(3).eq.com2%coord(3).and. &
          & com1%coord(4).eq.com2%coord(4))) THEN
        write(ierrou,*)'circular coplanar orbits: skip computation'
113     FORMAT(a3,1x,5(f10.5,1x))
        write(ierrou,113)'com1',com1%coord(1:2),com1%coord(3:5)*degrad
        write(ierrou,113)'com2',com2%coord(1:2),com2%coord(3:5)*degrad
        numerr=numerr+1
        RETURN
     ENDIF
  ENDIF
! check for overlapping ellipses
  IF((com1%coord(1).eq.com2%coord(1)).and. &
       & (com1%coord(2).eq.com2%coord(2)).and. &
       & (com1%coord(3).eq.com2%coord(3)).and. &
       & (com1%coord(4).eq.com2%coord(4)).and. &
       & (com1%coord(5).eq.com2%coord(5))) THEN
     write(ierrou,*)'overlapping orbits: skip computation'
     write(ierrou,113)'com1',com1%coord(1:2),com1%coord(3:5)*degrad
     write(ierrou,113)'com2',com2%coord(1:2),com2%coord(3:5)*degrad
     numerr = numerr+1
     RETURN
  ENDIF

! initialization 
  srtnum(1:poldeg)=0.d0 
! flags initialization  for crit_pts
  weier = .true.
  morse = .true.
  sflag(1) = .true. 
  sflag(2) = .true. 
  sflag(3) = .true.
  sflag(4) = .true.
  sflag(5) = .true.
  sflag(6) = .true.
  warnflag(1) = .true.
  warnflag(2) = .true.
  warnflag(3) = .true.
  hzflag = .true.
  hwflag = .true.
  multfl = .true.
  hevalflag = .true.
! -----------------------------------------------------------------------
  CALL crit_pts(com1%coord(1:5),com2%coord(1:5), &
       & f1,f2,nstat,nummin,nummax, &
       & answer,warnflag,sflag,morse,weier,hzflag,hwflag,multfl,hevalflag)
! -----------------------------------------------------------------------  
! f1 and f2 are given in degrees (in [-180,180])
  IF(.not.hevalflag) THEN
     write(ierrou,*)'dmintil_rms: hessian evaluation failed'
     numerr=numerr+1
  ENDIF
  IF(nstat.le.0.or.nummin.le.0) THEN
     WRITE(ierrou,*)'dmintil_rms: error! nstat=',nstat,'nummin=',nummin
     numerr=numerr+1
!     nummin=0 !non serve
     IF(PRESENT(dminflag)) dminflag(1)=1
     RETURN
  ENDIF
  
  DO j = 1,nstat 
! compute d^2 at critical (fpl,fcom)
     f1(j)=f1(j)*radeg 
     f2(j)=f2(j)*radeg 
     CALL d2eval(f1(j),f2(j),DD2)
     D2(j)=DD2 
  ENDDO
! sorting 3-uples (fcom,fpl,D2) according to value of D2               
  CALL heapsort(D2,nstat,srtnum)
  kmin = 0 !counting minimum points
  DO k = 1,nstat 
     IF(answer(srtnum(k)).eq.-1) THEN
        kmin=kmin+1
! added 2/6/2014, GFG ---------------------------------
        IF(kmin.gt.nminx)THEN
           WRITE(*,*)'dmintil_rms: nminx exceeded!',kmin
           WRITE(ierrou,*)'dmintil_rms: nminx exceeded!',kmin
           numerr=numerr+1
           nummin = kmin-1 !setting nummin to nminx
           CYCLE
        ENDIF
! -----------------------------------------------------
        f1min(kmin) = f1(srtnum(k))
        f2min(kmin) = f2(srtnum(k))
        D2min(kmin) = D2(srtnum(k))
        IF(kmin.gt.4)THEN
           WRITE(ierrou,*)' dmintil_rms: too many minima ', kmin,nummin
           numerr=numerr+1
        ENDIF
     ENDIF
  ENDDO
  IF(kmin.ne.nummin) THEN
     WRITE(ierrou,*)'dmintil_rms: kmin .ne. nummin',kmin,nummin
     numerr=numerr+1
  ENDIF
  IF(kmin.eq.0)THEN
     f1min(1) = f1(srtnum(1))
     f2min(1) = f2(srtnum(1))
     D2min(1) = D2(srtnum(1))
  ENDIF

  IF(PRESENT(unc1).and.(.not.PRESENT(unc2)))THEN
     WRITE(ierrou,*)'unc2 not present'
     numerr=numerr+1
     STOP
  ELSEIF(PRESENT(unc2).and.(.not.PRESENT(unc1)))THEN
     WRITE(ierrou,*)'unc1 not present'
     numerr=numerr+1
     STOP
  ELSEIF(PRESENT(unc1).and.PRESENT(unc2).and.PRESENT(dmintrms))THEN
     IF(PRESENT(chk_der))THEN
        chk_der_aux=chk_der
     ELSE
        chk_der_aux=.false.
     ENDIF
! -------------------------------------------------------------
! eigenvalue monitoring of the Hessian matrix
!(we have eliminated some terms that are zero at critical points)
     DO i=1,MAX(nummin,1)
!       CALL comp_rms_com(com1,com2,unccom1,unccom2,chk_der,f1min(i),f2min(i),&
!             & sinmutI,sinmutIrms,Hess,detH(i),trH(i),detHrms(i),car1min(i),&
!             & car2min(i),dcardcom1,dcardcom2,tau1,tau2,tau3hat,sint1t2(i),&
!             & taurms(i))
        CALL comp_rms_com(com1,com2,unccom1,unccom2,chk_der_aux, &
             & f1min(i),f2min(i),sinmutI_aux,sinmutIrms_aux,Hess, &
             & detH_aux,trH(i),detHrms_aux,car1min(i),car2min(i), &
             & dcardcom1,dcardcom2,tau1,tau2,tau3hat,sint1t2_aux(i),&
             & taurms_aux)
        IF(PRESENT(tau_1).AND.PRESENT(tau_2))THEN
           tau_1(1:3,i)=tau1(1:3)
           tau_2(1:3,i)=tau2(1:3)
        ENDIF
        IF(PRESENT(detH))THEN
           detH(i)=detH_aux
        ENDIF
        IF(PRESENT(detHrms))THEN
           detHrms(i)=detHrms_aux
        ENDIF
        IF(PRESENT(sint1t2))THEN
           sint1t2(i)=sint1t2_aux(i)
        ENDIF
        IF(PRESENT(taurms))THEN
           taurms(i)=taurms_aux
        ENDIF
        IF(PRESENT(sinmutI))THEN
           sinmutI=sinmutI_aux
        ENDIF
        IF(PRESENT(sinmutIrms))THEN
           sinmutIrms=sinmutIrms_aux
        ENDIF
        CALL check_rms(sinmutI_aux,sinmutIrms_aux,detH_aux,detHrms_aux, &
             & sint1t2_aux(i),taurms_aux,comp_flag_aux(1:3))
        IF(PRESENT(comp_flag))THEN
           comp_flag(1:3,i)=comp_flag_aux(1:3)
        ENDIF
        dcardcom(1:6,1:10)=0.d0! initialization
        dcardcom(1:3,1:5)=dcardcom1(1:3,1:5)
        dcardcom(4:6,6:10)=dcardcom2(1:3,1:5)
        DO h=1,10       
           DO k=1,3 
              ddeltadcom(k)=dcardcom(k,h)-dcardcom(k+3,h)
           ENDDO
           derdmintil(1,h)=DOT_PRODUCT(tau3hat,ddeltadcom)!derivs w.r.t.COM
        ENDDO
! give a sign to the minimal distance
        CALL sign_dmin(car1min(i)%coord(1:3),tau1, &
             & car2min(i)%coord(1:3),tau2,dsign)
        IF(PRESENT(dminflag))THEN
           IF(dsign.eq.2.d0)THEN
              dminflag(i)=2
           ENDIF
        ENDIF
        dmintil(i)=dsign*SQRT(D2min(i))
        tderdmintil=TRANSPOSE(derdmintil)
        dmint_cov=MATMUL(derdmintil,MATMUL(covcom,tderdmintil)) 
        dmintrms(i)=SQRT(dmint_cov(1,1)) ! scalar value
        IF(PRESENT(c1min).and.PRESENT(c2min))THEN
           c1min(1:3,i) = car1min(i)%coord(1:3)
           c2min(1:3,i) = car2min(i)%coord(1:3)
        ENDIF
        IF(PRESENT(ddmintdel2))THEN
           ddmintdel2(1:5,i)=derdmintil(1,6:10)
        ENDIF
        IF(chk_der_aux)THEN
           write(ierrou,*)'derivatives of dmintilde w.r.t COM:'
           write(ierrou,100) derdmintil(1,6:10)
           numerr=numerr+1
        ENDIF
     ENDDO


  ELSEIF(PRESENT(ddmintdel2))THEN
! -------------------------------------------------------------
! eigenvalue monitoring of the Hessian matrix
!(we have eliminated some terms that are zero at critical points)
     DO i=1,MAX(nummin,1)
        com1min=com1
        com2min=com2
! substitute last element time of perihelion passage
        com1min%coord(6)=f1min(i)
        com2min%coord(6)=f2min(i)
        CALL coo_cha(com1min,'CAR',car1min(i),fail_flag,dcardcom1)
        CALL coo_cha(com2min,'CAR',car2min(i),fail_flag,dcardcom2)       

        CALL tau1_tau2(com1,com2,f1min(i),f2min(i),car1min(i),car2min(i),&
             & tau1,tau2,tau3hat,sint1t2_aux(i))
        IF(PRESENT(tau_1).AND.PRESENT(tau_2))THEN
           tau_1(1:3,i)=tau1(1:3)
           tau_2(1:3,i)=tau2(1:3)
        ENDIF
!        CALL comp_rms_com(com1,com2,unccom1,unccom2,chk_der_aux, &
!             & f1min(i),f2min(i),sinmutI_aux,sinmutIrms_aux,Hess, &
!             & detH_aux,trH(i),detHrms_aux,car1min(i),car2min(i), &
!             & dcardcom1,dcardcom2,tau1,tau2,tau3hat,sint1t2_aux(i),&
!             & taurms_aux)

        dcardcom(1:6,1:10)=0.d0! initialization
        dcardcom(1:3,1:5)=dcardcom1(1:3,1:5)
        dcardcom(4:6,6:10)=dcardcom2(1:3,1:5)
        DO h=1,10       
           DO k=1,3 
              ddeltadcom(k)=dcardcom(k,h)-dcardcom(k+3,h)
           ENDDO
           derdmintil(1,h)=DOT_PRODUCT(tau3hat,ddeltadcom)!derivs w.r.t.COM
        ENDDO
! give a sign to the minimal distance
        CALL sign_dmin(car1min(i)%coord(1:3),tau1, &
             & car2min(i)%coord(1:3),tau2,dsign)
        IF(PRESENT(dminflag))THEN
           IF(dsign.eq.2.d0)THEN
              dminflag(i)=2
           ENDIF
        ENDIF
        dmintil(i)=dsign*SQRT(D2min(i))
        tderdmintil=TRANSPOSE(derdmintil)
        IF(PRESENT(c1min).and.PRESENT(c2min))THEN
           c1min(1:3,i) = car1min(i)%coord(1:3)
           c2min(1:3,i) = car2min(i)%coord(1:3)
        ENDIF

        ddmintdel2(1:5,i)=derdmintil(1,6:10)
     ENDDO

  ELSE
     DO i=1,nummin
        CALL tau1_tau2(com1,com2,f1min(i),f2min(i),car1min(i),car2min(i),&
             & tau1,tau2,tau3hat,sint1t2_aux(i))
        IF(PRESENT(tau_1).AND.PRESENT(tau_2))THEN
           tau_1(1:3,i)=tau1(1:3)
           tau_2(1:3,i)=tau2(1:3)
        ENDIF
        CALL sign_dmin(car1min(i)%coord(1:3),tau1, &
             & car2min(i)%coord(1:3),tau2,dsign)

! ======================================
! ADDED 12/2/2012
        IF(PRESENT(dminflag))THEN
           IF(dsign.eq.2.d0)THEN
              dminflag(i) = 2
              CYCLE
           ENDIF
        ENDIF
! =======================================
        dmintil(i)=dsign*SQRT(D2min(i)) 
        IF(PRESENT(c1min).and.PRESENT(c2min))THEN
           c1min(1:3,i) = car1min(i)%coord(1:3)
           c2min(1:3,i) = car2min(i)%coord(1:3)
        ENDIF
     ENDDO
  ENDIF
  IF(PRESENT(sint1t2))THEN
     sint1t2(1:nummin)=sint1t2_aux(1:nummin)
  ENDIF

  IF(PRESENT(chk_der).and.PRESENT(sinmutI).and.PRESENT(detH) &
       & .and.PRESENT(sint1t2))THEN
     IF(chk_der)THEN
        CALL derdmintest(com1,com2,nstat,nummin,f1min,f2min, &
             & sinmutI,D2min,detH,sint1t2)
!     write(*,*)'derivatives of dmintilde w.r.t COM:'
!     write(*,100) derdmintil(1,6:10)
 100 FORMAT(5(f12.8,1x))
     ENDIF
  ENDIF

  IF(PRESENT(v1min).AND.PRESENT(v2min)) THEN
     v1min(1:nummin)=f1min(1:nummin)
     v2min(1:nummin)=f2min(1:nummin)
  ENDIF

END SUBROUTINE dmintil_rms
! ============================================================
SUBROUTINE check_rms(sinmutI,sinmutIrms,detH,detHrms,st1t2,&
          & st1t2rms,comp_flag)
  IMPLICIT NONE
  DOUBLE PRECISION,INTENT(IN) :: sinmutI,sinmutIrms,detH,detHrms,st1t2,st1t2rms
  LOGICAL,DIMENSION(3),INTENT(OUT) :: comp_flag
! -------------- end interface -----------------
  comp_flag(1:3)=.true.! initialization
  IF(SIGN(1.d0,(sinmutI-3.d0*sinmutIrms)*(sinmutI+3.d0*sinmutIrms)).le.0.d0)THEN
     comp_flag(1)=.false.
  ENDIF
  IF(SIGN(1.d0,(detH-3.d0*detHrms)*(detH+3.d0*detHrms)).le.0.d0)THEN
     comp_flag(2)=.false.
  ENDIF
  IF(SIGN(1.d0,(st1t2-3.d0*st1t2rms)*(st1t2+3.d0*st1t2rms)).le.0.d0)THEN
     comp_flag(3)=.false.
  ENDIF

END SUBROUTINE check_rms

! =====================================================================
! compute rms of sinmutI,detH,sint1t2 by covariance propagation
! written by G.F. Gronchi and G. Tommei (10/2/2006)
! last modified 31/01/2007
SUBROUTINE comp_rms_com(com1,com2,unc1,unc2,chk_der,f1,f2,sinmutI,&
     &sinmutIrms,Hess,detH,trH,detHrms,car1min,car2min,dcardcom1,dcardcom2,&
     &tau1,tau2,tau3hat,sint1t2,taurms)
  USE ever_pitkin
  IMPLICIT NONE
!===================================================

  TYPE(orbit_elem), INTENT(IN) :: com1,com2 ! cometary elements
                                          ! 1: planet 
                                          ! 2: asteroid
  TYPE(orb_uncert), INTENT(IN) :: unc1,unc2
  LOGICAL,INTENT(IN) :: chk_der
  DOUBLE PRECISION, INTENT(IN) :: f1,f2   ! anomalies at the critical point
                                          ! in [0,2*pi]
  DOUBLE PRECISION, INTENT(OUT) :: sinmutI
  DOUBLE PRECISION,INTENT(OUT) :: sinmutIrms
  DOUBLE PRECISION, DIMENSION(2,2), INTENT(OUT) :: Hess 
  DOUBLE PRECISION, INTENT(OUT) :: detH,trH
  DOUBLE PRECISION, INTENT(OUT) :: detHrms
  TYPE(orbit_elem),INTENT(OUT) :: car1min,car2min
  DOUBLE PRECISION,DIMENSION(6,6),INTENT(OUT) :: dcardcom1,dcardcom2
  DOUBLE PRECISION,DIMENSION(3),INTENT(OUT) :: tau1,tau2,tau3hat
! absolute val of the sinus of the angle between tau1 and tau2
  DOUBLE PRECISION,INTENT(OUT) :: sint1t2
  DOUBLE PRECISION,INTENT(OUT) ::taurms
!===================================================
  DOUBLE PRECISION :: cosmutI
  DOUBLE PRECISION :: dcosmutIdi1,dcosmutIdOm1
  DOUBLE PRECISION :: dcosmutIdi2,dcosmutIdOm2
  DOUBLE PRECISION, DIMENSION(1,10) :: dersinmutIdcom
  DOUBLE PRECISION, DIMENSION(1,5) :: dersinmutIdcom2
  DOUBLE PRECISION, DIMENSION(10,1) :: tdersinmutIdcom
  DOUBLE PRECISION, DIMENSION(1,1) :: sinmutIrmsvec
!===============================================================
!  TYPE(orbit_elem) :: elcpl,elc
  INTEGER :: fail_flag
  DOUBLE PRECISION, DIMENSION(6,6) :: gcom1,gcom2
!  DOUBLE PRECISION, DIMENSION(6,6) :: jackepcom1,jackepcom2
!  DOUBLE PRECISION, DIMENSION(6,6) :: tjackepcom1,tjackepcom2
  DOUBLE PRECISION, DIMENSION(6,6) :: gamma1,gamma2
  DOUBLE PRECISION, DIMENSION(10,10) :: covcom
  DOUBLE PRECISION :: q1,e1,i1,Om1,omeg1  ! Planet
  DOUBLE PRECISION :: q2,e2,i2,Om2,omeg2  ! Asteroid/Comet
  DOUBLE PRECISION :: s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12,s13,s14,s15
  DOUBLE PRECISION, DIMENSION(2,2) :: invH 
!---------------------------------------------------
  DOUBLE PRECISION :: dKKdq1,dKKde1,dKKdi1,dKKdOm1,dKKdomeg1
  DOUBLE PRECISION :: dKKdq2,dKKde2,dKKdi2,dKKdOm2,dKKdomeg2
  DOUBLE PRECISION :: dLLdq1,dLLde1,dLLdi1,dLLdOm1,dLLdomeg1
  DOUBLE PRECISION :: dLLdq2,dLLde2,dLLdi2,dLLdOm2,dLLdomeg2
  DOUBLE PRECISION :: dMMdq1,dMMde1,dMMdi1,dMMdOm1,dMMdomeg1
  DOUBLE PRECISION :: dMMdq2,dMMde2,dMMdi2,dMMdOm2,dMMdomeg2
  DOUBLE PRECISION :: dNNdq1,dNNde1,dNNdi1,dNNdOm1,dNNdomeg1
  DOUBLE PRECISION :: dNNdq2,dNNde2,dNNdi2,dNNdOm2,dNNdomeg2
!----------------------------------------------------
  DOUBLE PRECISION :: gradd2f1,gradd2f2
  DOUBLE PRECISION :: dgradd2f1_dKK,dgradd2f1_dLL,dgradd2f1_dMM,dgradd2f1_dNN
  DOUBLE PRECISION :: dgradd2f2_dKK,dgradd2f2_dLL,dgradd2f2_dMM,dgradd2f2_dNN
  DOUBLE PRECISION :: dgradd2f1_dq1,dgradd2f1_dq2
  DOUBLE PRECISION :: dgradd2f1_de1,dgradd2f1_de2
  DOUBLE PRECISION :: dgradd2f1_di1,dgradd2f1_di2
  DOUBLE PRECISION :: dgradd2f1_dOm1,dgradd2f1_dOm2
  DOUBLE PRECISION :: dgradd2f1_domeg1,dgradd2f1_domeg2
  DOUBLE PRECISION :: dgradd2f2_dq1,dgradd2f2_dq2
  DOUBLE PRECISION :: dgradd2f2_de1,dgradd2f2_de2
  DOUBLE PRECISION :: dgradd2f2_di1,dgradd2f2_di2
  DOUBLE PRECISION :: dgradd2f2_dOm1,dgradd2f2_dOm2
  DOUBLE PRECISION :: dgradd2f2_domeg1,dgradd2f2_domeg2
!----------------------------------------------------
  DOUBLE PRECISION :: df1_dq1,df2_dq1,df1_de1,df2_de1,df1_di1,df2_di1
  DOUBLE PRECISION :: df1_dOm1,df2_dOm1,df1_domeg1,df2_domeg1
  DOUBLE PRECISION :: df1_dq2,df2_dq2,df1_de2,df2_de2,df1_di2,df2_di2
  DOUBLE PRECISION :: df1_dOm2,df2_dOm2,df1_domeg2,df2_domeg2
!----------------------------------------------------
  DOUBLE PRECISION :: dH11_df1,dH11_df2,dH22_df1,dH22_df2
  DOUBLE PRECISION :: dH12_df1,dH12_df2
!----------------------------------------------------
  DOUBLE PRECISION :: dH11_dKK,dH11_dLL,dH11_dMM,dH11_dNN
  DOUBLE PRECISION :: dH22_dKK,dH22_dLL,dH22_dMM,dH22_dNN
  DOUBLE PRECISION :: dH12_dKK,dH12_dLL,dH12_dMM,dH12_dNN
!----------------------------------------------------
  DOUBLE PRECISION :: dH11_dq1,dH11_de1,dH11_di1,dH11_dOm1,dH11_domeg1
  DOUBLE PRECISION :: dH11_dq2,dH11_de2,dH11_di2,dH11_dOm2,dH11_domeg2
  DOUBLE PRECISION :: dH22_dq1,dH22_de1,dH22_di1,dH22_dOm1,dH22_domeg1
  DOUBLE PRECISION :: dH22_dq2,dH22_de2,dH22_di2,dH22_dOm2,dH22_domeg2
  DOUBLE PRECISION :: dH12_dq1,dH12_de1,dH12_di1,dH12_dOm1,dH12_domeg1
  DOUBLE PRECISION :: dH12_dq2,dH12_de2,dH12_di2,dH12_dOm2,dH12_domeg2
!----------------------------------------------------
  DOUBLE PRECISION :: derH11_dq1,derH11_de1,derH11_di1
  DOUBLE PRECISION :: derH11_dOm1,derH11_domeg1
  DOUBLE PRECISION :: derH11_dq2,derH11_de2,derH11_di2
  DOUBLE PRECISION :: derH11_dOm2,derH11_domeg2
  DOUBLE PRECISION :: derH22_dq1,derH22_de1,derH22_di1
  DOUBLE PRECISION :: derH22_dOm1,derH22_domeg1
  DOUBLE PRECISION :: derH22_dq2,derH22_de2,derH22_di2
  DOUBLE PRECISION :: derH22_dOm2,derH22_domeg2
  DOUBLE PRECISION :: derH12_dq1,derH12_de1,derH12_di1
  DOUBLE PRECISION :: derH12_dOm1,derH12_domeg1
  DOUBLE PRECISION :: derH12_dq2,derH12_de2,derH12_di2
  DOUBLE PRECISION :: derH12_dOm2,derH12_domeg2
!----------------------------------------------------
  DOUBLE PRECISION :: ddetH_dq1,ddetH_de1,ddetH_di1
  DOUBLE PRECISION :: ddetH_dOm1,ddetH_domeg1
  DOUBLE PRECISION :: ddetH_dq2,ddetH_de2,ddetH_di2
  DOUBLE PRECISION :: ddetH_dOm2,ddetH_domeg2
!----------------------------------------------------
  DOUBLE PRECISION, DIMENSION(1,10) :: ddetH
!  DOUBLE PRECISION, DIMENSION(1,5) :: ddetHdcom2
  DOUBLE PRECISION, DIMENSION(10,1) :: tddetH
  DOUBLE PRECISION, DIMENSION(1,1) :: detHrmsvec
!====================================================
 DOUBLE PRECISION ::l1,l2,lambda1,lambda2,u1,u2
! to compute derivatives of dmintil
  TYPE(orbit_elem) :: com1min,com2min
  DOUBLE PRECISION :: vsize,ltau1,ltau2 
  DOUBLE PRECISION,DIMENSION(3) :: tau1vers,tau2vers,tau3vec,tau3
  DOUBLE PRECISION,DIMENSION(3) :: dtau1versdq1,dtau1versde1,dtau1versdi1
  DOUBLE PRECISION,DIMENSION(3) :: dtau1versdOm1,dtau1versdomeg1,dtau1versdf1
  DOUBLE PRECISION,DIMENSION(3) :: dtau2versdq2,dtau2versde2,dtau2versdi2
  DOUBLE PRECISION,DIMENSION(3) :: dtau2versdOm2,dtau2versdomeg2,dtau2versdf2
  DOUBLE PRECISION,DIMENSION(3) :: dtau3vecdq1,dtau3vecde1,dtau3vecdi1
  DOUBLE PRECISION,DIMENSION(3) :: dtau3vecdOm1,dtau3vecdomeg1,dtau3vecdf1
  DOUBLE PRECISION,DIMENSION(3) :: dtau3vecdq2,dtau3vecde2,dtau3vecdi2
  DOUBLE PRECISION,DIMENSION(3) :: dtau3vecdOm2,dtau3vecdomeg2,dtau3vecdf2
!
  DOUBLE PRECISION :: dsint1t2dq1,dsint1t2de1,dsint1t2di1
  DOUBLE PRECISION :: dsint1t2dOm1,dsint1t2domeg1,dsint1t2df1
  DOUBLE PRECISION :: dsint1t2dq2,dsint1t2de2,dsint1t2di2
  DOUBLE PRECISION :: dsint1t2dOm2,dsint1t2domeg2,dsint1t2df2
  DOUBLE PRECISION,DIMENSION(3) :: dtau1dq1,dtau1de1,dtau1di1
  DOUBLE PRECISION,DIMENSION(3) :: dtau1dOm1,dtau1domeg1,dtau1df1
  DOUBLE PRECISION,DIMENSION(3) :: dtau2dq2,dtau2de2,dtau2di2
  DOUBLE PRECISION,DIMENSION(3) :: dtau2dOm2,dtau2domeg2,dtau2df2
  DOUBLE PRECISION, DIMENSION(1,10) :: dersint1t2dcom
!  DOUBLE PRECISION, DIMENSION(1,5) :: dersint1t2dcom2
  DOUBLE PRECISION, DIMENSION(10,1) :: tdersint1t2dcom
  DOUBLE PRECISION, DIMENSION(1,1) :: taurmsvec
! to compute time of passage at perihelion
  DOUBLE PRECISION ::cosn1,sinn1,coso1,sino1,cosi1,sini1
  DOUBLE PRECISION ::cosn2,sinn2,coso2,sino2,cosi2,sini2
  DOUBLE PRECISION,DIMENSION(3) :: x01,y01,x02,y02 
  DOUBLE PRECISION :: r01,v01,r1,mu1,alpha1,psi1,dt1,ptime1
  DOUBLE PRECISION :: r02,v02,r2,mu2,alpha2,psi2,dt2,ptime2
  DOUBLE PRECISION :: psi1_ini,psi2_ini,u1x,u2x
! ==========================================================

! covariance matrix
    gcom1=unc1%g(1:6,1:6)
    gcom2=unc2%g(1:6,1:6)
    covcom(1:10,1:10)=0.d0
    covcom(1:5,1:5)=gcom1(1:5,1:5)
    covcom(6:10,6:10)=gcom2(1:5,1:5)

! elements
    q1 = com1%coord(1) 
    e1 = com1%coord(2)
    i1 = com1%coord(3)
    Om1 = com1%coord(4)
    omeg1 = com1%coord(5)
    q2 = com2%coord(1) 
    e2 = com2%coord(2) 
    i2 = com2%coord(3)
    Om2 = com2%coord(4)
    omeg2 = com2%coord(5)

! ================================
! RMS of the mutual inclination
! ================================
    cosmutI= sin(Om1)*sin(i1)*sin(Om2)*sin(i2)+cos(Om1)* &
         & sin(i1)*cos(Om2)*sin(i2)+cos(i1)*cos(i2)
    sinmutI=sqrt(1-cosmutI**2)
    dcosmutIdi1 = sin(Om1)*cos(i1)*sin(Om2)*sin(i2)+ &
         & cos(Om1)*cos(i1)*cos(Om2)*sin(i2)-sin(i1)*cos(i2)
    dcosmutIdOm1 = cos(Om1)*sin(i1)*sin(Om2)*sin(i2)- &
         & sin(Om1)*sin(i1)*cos(Om2)*sin(i2)
    dcosmutIdi2= sin(Om1)*sin(i1)*sin(Om2)*cos(i2)+ &
         & cos(Om1)*sin(i1)*cos(Om2)*cos(i2)-cos(i1)*sin(i2)
    dcosmutIdOm2= sin(Om1)*sin(i1)*cos(Om2)*sin(i2)- &
         & cos(Om1)*sin(i1)*sin(Om2)*sin(i2)
    dersinmutIdcom(1,1)=0.d0
    dersinmutIdcom(1,2)=0.d0
    dersinmutIdcom(1,3)=(-cosmutI/sinmutI)*dcosmutIdi1
    dersinmutIdcom(1,4)=(-cosmutI/sinmutI)*dcosmutIdOm1
    dersinmutIdcom(1,5)=0.d0
    dersinmutIdcom(1,6)=0.d0
    dersinmutIdcom(1,7)=0.d0
    dersinmutIdcom(1,8)=(-cosmutI/sinmutI)*dcosmutIdi2
    dersinmutIdcom(1,9)=(-cosmutI/sinmutI)*dcosmutIdOm2
    dersinmutIdcom(1,10)=0.d0
    tdersinmutIdcom=TRANSPOSE(dersinmutIdcom)
    sinmutIrmsvec=MATMUL(dersinmutIdcom,MATMUL(covcom,tdersinmutIdcom)) 
    sinmutIrms=SQRT(sinmutIrmsvec(1,1))

    IF(chk_der) THEN
!       dersinmutIdcom2(1:1,1:5)= dersinmutIdcom(1:1,6:10)
       write(ierrou,*)'-----------------------------------------&
            &-------------------------'
       write(ierrou,*)'derivatives of sinmutI:'
       write(ierrou,100) dersinmutIdcom(1:1,6:10)
       numerr=numerr+1
    ENDIF

!===================================================================
! RMS of the determinant of the Hessian matrix
!====================================================================
! GRADIENT of d^2
!====================================================================
    gradd2f1                                                          &
     &= 2.D0/(1.D0+e1*cos(f1))*(e1*q1**2*(1.D0+e1)**2.D0/(1.D0+e1*cos   &
     &(f1))**2.D0*sin(f1)+q1*(1.D0+e1)/(1.D0+e1*cos(f1))*sin(f1)*(KK*q2*&
     &(1.D0+e2)/(1.D0+e2*cos(f2))*cos(f2)+MM*q2*(1.D0+e2)/(1.D0+e2*cos(f&
     &2))*sin(f2))-(e1*q1*(1.D0+e1)/(1.D0+e1*cos(f1))+q1*(1.D0+e1)/(1.D0&
     &+e1*cos(f1))*cos(f1))*(LL*q2*(1.D0+e2)/(1.D0+e2*cos(f2))*cos(f2)+N&
     &N*q2*(1.D0+e2)/(1.D0+e2*cos(f2))*sin(f2)))                        
                                                                        
    gradd2f2                                                          &
     &= 2.D0/(1.D0+e2*cos(f2))*(e2*q2**2*(1.D0+e2)**2.D0/(1.D0+e2*cos   &
     &(f2))**2.D0*sin(f2)+q2*(1.D0+e2)/(1.D0+e2*cos(f2))*sin(f2)*(KK*q1*&
     &(1.D0+e1)/(1.D0+e1*cos(f1))*cos(f1)+LL*q1*(1.D0+e1)/(1.D0+e1*cos(f&
     &1))*sin(f1))-(e2*q2*(1.D0+e2)/(1.D0+e2*cos(f2))+q2*(1.D0+e2)/(1.D0&
     &+e2*cos(f2))*cos(f2))*(MM*q1*(1.D0+e1)/(1.D0+e1*cos(f1))*cos(f1)+N&
     &N*q1*(1.D0+e1)/(1.D0+e1*cos(f1))*sin(f1)))                        

! force to zero (critical point)
!      gradd2f1=0.D0
!      WRITE(*,*)'gradd2f1:',gradd2f1
!      gradd2f2=0.D0
!      WRITE(*,*)'gradd2f2:',gradd2f2

!====================================================================
! Hessian Matrix H of d^2
!====================================================================

    Hess(1,1)&
     &= e1/(1.D0+e1*cos(f1))*sin(f1)*gradd2f1+2.D0/(1.D0+e1*cos(f1))*&
     &(e1**2*q1**2*(1.D0+e1)**2.D0/(1.D0+e1*cos(f1))**3.D0*sin(f1)**2.D0&
     &+(e1*q1*(1.D0+e1)/(1.D0+e1*cos(f1))**2.D0*sin(f1)**2.D0+q1*(1.D0+e&
     &1)/(1.D0+e1*cos(f1))*cos(f1))*(e1*q1*(1.D0+e1)/(1.D0+e1*cos(f1))+K&
     &K*q2*(1.D0+e2)/(1.D0+e2*cos(f2))*cos(f2)+MM*q2*(1.D0+e2)/(1.D0+e2*&
     &cos(f2))*sin(f2))-(e1**2*q1*(1.D0+e1)/(1.D0+e1*cos(f1))**2.D0*sin(&
     &f1)+e1*q1*(1.D0+e1)/(1.D0+e1*cos(f1))**2.D0*cos(f1)*sin(f1)-q1*(1.&
     &D0+e1)/(1.D0+e1*cos(f1))*sin(f1))*(LL*q2*(1.D0+e2)/(1.D0+e2*cos(f2&
     &))*cos(f2)+NN*q2*(1.D0+e2)/(1.D0+e2*cos(f2))*sin(f2)))            
                                                                        
    Hess(2,2)&
     &= e2/(1.D0+e2*cos(f2))*sin(f2)*gradd2f2+2.D0/(1.D0+e2*cos(f2))*&
     &(e2**2*q2**2*(1.D0+e2)**2.D0/(1.D0+e2*cos(f2))**3.D0*sin(f2)**2.D0&
     &+(e2*q2*(1.D0+e2)/(1.D0+e2*cos(f2))**2.D0*sin(f2)**2.D0+q2*(1.D0+e&
     &2)/(1.D0+e2*cos(f2))*cos(f2))*(e2*q2*(1.D0+e2)/(1.D0+e2*cos(f2))+K&
     &K*q1*(1.D0+e1)/(1.D0+e1*cos(f1))*cos(f1)+LL*q1*(1.D0+e1)/(1.D0+e1*&
     &cos(f1))*sin(f1))-(e2**2*q2*(1.D0+e2)/(1.D0+e2*cos(f2))**2.D0*sin(&
     &f2)+e2*q2*(1.D0+e2)/(1.D0+e2*cos(f2))**2.D0*cos(f2)*sin(f2)-q2*(1.&
     &D0+e2)/(1.D0+e2*cos(f2))*sin(f2))*(MM*q1*(1.D0+e1)/(1.D0+e1*cos(f1&
     &))*cos(f1)+NN*q1*(1.D0+e1)/(1.D0+e1*cos(f1))*sin(f1)))            
                                                                        
    Hess(1,2)&
     &= 2.D0/(1.D0+e1*cos(f1))*(q1*(1.D0+e1)/(1.D0+e1*cos(f1))*sin(f1&
     &)*(KK*(e2*q2*(1.D0+e2)/(1.D0+e2*cos(f2))**2.D0*cos(f2)*sin(f2)-q2*&
     &(1.D0+e2)/(1.D0+e2*cos(f2))*sin(f2))+MM*(e2*q2*(1.D0+e2)/(1.D0+e2*&
     &cos(f2))**2.D0*sin(f2)**2.D0+q2*(1.D0+e2)/(1.D0+e2*cos(f2))*cos(f2&
     &)))-(e1*q1*(1.D0+e1)/(1.D0+e1*cos(f1))+q1*(1.D0+e1)/(1.D0+e1*cos(f&
     &1))*cos(f1))*(LL*(e2*q2*(1.D0+e2)/(1.D0+e2*cos(f2))**2.D0*cos(f2)*&
     &sin(f2)-q2*(1.D0+e2)/(1.D0+e2*cos(f2))*sin(f2))+NN*(e2*q2*(1.D0+e2&
     &)/(1.D0+e2*cos(f2))**2.D0*sin(f2)**2.D0+q2*(1.D0+e2)/(1.D0+e2*cos(&
     &f2))*cos(f2))))                                                   
                                                                        
    Hess(2,1)=Hess(1,2)

! Determinant of Hess
    detH=Hess(1,1)*Hess(2,2)-Hess(1,2)**2
! Trace of Hess
    trH=Hess(1,1)+Hess(2,2)
! Inverse matrix of H
    IF(detH.ne.0)THEN
       invH(1,1)=Hess(2,2)/detH
       invH(2,2)=Hess(1,1)/detH
       invH(1,2)=-Hess(1,2)/detH
       invH(2,1)=-Hess(2,1)/detH
    ELSE
       WRITE(ierrou,*)'detH_rms: detH= ',detH
       numerr=numerr+1
    ENDIF

!=========================================================
! Derivatives of KK,LL,MM,NN w.r.t. the elements
!=========================================================
    dKKdq1=0.d0
    dKKde1=0.d0
    dKKdq2=0.d0
    dKKde2=0.d0
!------------------------------------------------------------
    dLLdq1=0.d0
    dLLde1=0.d0
    dLLdq2=0.d0
    dLLde2=0.d0
!------------------------------------------------------------
    dMMdq1=0.d0
    dMMde1=0.d0
    dMMdq2=0.d0
    dMMde2=0.d0
!-------------------------------------------------------------
    dNNdq1=0.d0
    dNNde1=0.d0
    dNNdq2=0.d0
    dNNde2=0.d0

    dKKdi1                                                            &
     &=-sin(omeg1)*sin(i1)*(-cos(omeg2)*sin(-Om2+Om1)+sin(omeg2)*cos    &
     &(i2)*cos(-Om2+Om1))+sin(omeg1)*cos(i1)*sin(omeg2)*sin(i2)         
                                                                        
    dKKdi2                                                            &
     &= -cos(omeg1)*sin(omeg2)*sin(i2)*sin(-Om2+Om1)-sin(omeg1)*cos(i&
     &1)*sin(omeg2)*sin(i2)*cos(-Om2+Om1)+sin(omeg1)*sin(i1)*sin(omeg2)*&
     &cos(i2)                                                           
                                                                        
    dKKdOm1                                                           &
     &= cos(omeg1)*(-cos(omeg2)*sin(-Om2+Om1)+sin(omeg2)*cos(i2)*cos(   &
     &-Om2+Om1))+sin(omeg1)*cos(i1)*(-cos(omeg2)*cos(-Om2+Om1)-sin(omeg2&
     &)*cos(i2)*sin(-Om2+Om1))                                          
                                                                        
    dKKdOm2                                                           &
     &= cos(omeg1)*(cos(omeg2)*sin(-Om2+Om1)-sin(omeg2)*cos(i2)*cos(-   &
     &Om2+Om1))+sin(omeg1)*cos(i1)*(cos(omeg2)*cos(-Om2+Om1)+sin(omeg2)*&
     &cos(i2)*sin(-Om2+Om1))                                            
                                                                        
    dKKdomeg1                                                         &
     &= -sin(omeg1)*(cos(omeg2)*cos(-Om2+Om1)+sin(omeg2)*cos(i2)*sin(   &
     &-Om2+Om1))+cos(omeg1)*cos(i1)*(-cos(omeg2)*sin(-Om2+Om1)+sin(omeg2&
     &)*cos(i2)*cos(-Om2+Om1))+cos(omeg1)*sin(i1)*sin(omeg2)*sin(i2)    
                                                                        
    dKKdomeg2                                                         &
     &= cos(omeg1)*(-sin(omeg2)*cos(-Om2+Om1)+cos(omeg2)*cos(i2)*sin(   &
     &-Om2+Om1))+sin(omeg1)*cos(i1)*(sin(omeg2)*sin(-Om2+Om1)+cos(omeg2)&
     &*cos(i2)*cos(-Om2+Om1))+sin(omeg1)*sin(i1)*cos(omeg2)*sin(i2)     
                                                                        
    dLLdi1                                                            &
     &= -cos(omeg1)*sin(i1)*(-cos(omeg2)*sin(-Om2+Om1)+sin(omeg2)*cos   &
     &(i2)*cos(-Om2+Om1))+cos(omeg1)*cos(i1)*sin(omeg2)*sin(i2)         
                                                                        
    dLLdi2                                                            &
     &= sin(omeg1)*sin(omeg2)*sin(i2)*sin(-Om2+Om1)-cos(omeg1)*cos(i1   &
     &)*sin(omeg2)*sin(i2)*cos(-Om2+Om1)+cos(omeg1)*sin(i1)*sin(omeg2)*c&
     &os(i2)                                                            
                                                                        
    dLLdOm1                                                           &
     &= -sin(omeg1)*(-cos(omeg2)*sin(-Om2+Om1)+sin(omeg2)*cos(i2)*cos&
     &(-Om2+Om1))+cos(omeg1)*cos(i1)*(-cos(omeg2)*cos(-Om2+Om1)-sin(omeg&
     &2)*cos(i2)*sin(-Om2+Om1))                                         
                                                                        
    dLLdOm2                                                           &
     &= -sin(omeg1)*(cos(omeg2)*sin(-Om2+Om1)-sin(omeg2)*cos(i2)*cos(   &
     &-Om2+Om1))+cos(omeg1)*cos(i1)*(cos(omeg2)*cos(-Om2+Om1)+sin(omeg2)&
     &*cos(i2)*sin(-Om2+Om1))                                           
                                                                        
    dLLdomeg1                                                         &
     &= -cos(omeg1)*(cos(omeg2)*cos(-Om2+Om1)+sin(omeg2)*cos(i2)*sin(   &
     &-Om2+Om1))-sin(omeg1)*cos(i1)*(-cos(omeg2)*sin(-Om2+Om1)+sin(omeg2&
     &)*cos(i2)*cos(-Om2+Om1))-sin(omeg1)*sin(i1)*sin(omeg2)*sin(i2)    
    
    dLLdomeg2                                                         &
     &= -sin(omeg1)*(-sin(omeg2)*cos(-Om2+Om1)+cos(omeg2)*cos(i2)*sin   &
     &(-Om2+Om1))+cos(omeg1)*cos(i1)*(sin(omeg2)*sin(-Om2+Om1)+cos(omeg2&
     &)*cos(i2)*cos(-Om2+Om1))+cos(omeg1)*sin(i1)*cos(omeg2)*sin(i2)    
                                                                        
    dMMdi1                                                            &
     &= -sin(omeg1)*sin(i1)*(sin(omeg2)*sin(-Om2+Om1)+cos(omeg2)*cos(   &
     &i2)*cos(-Om2+Om1))+sin(omeg1)*cos(i1)*cos(omeg2)*sin(i2)          
                                                                        
    dMMdi2                                                            &
     &= -cos(omeg1)*cos(omeg2)*sin(i2)*sin(-Om2+Om1)-sin(omeg1)*cos(i&
     &1)*cos(omeg2)*sin(i2)*cos(-Om2+Om1)+sin(omeg1)*sin(i1)*cos(omeg2)*&
     &cos(i2)                                                           
                                                                        
    dMMdOm1                                                           &
     &= cos(omeg1)*(sin(omeg2)*sin(-Om2+Om1)+cos(omeg2)*cos(i2)*cos(-   &
     &Om2+Om1))+sin(omeg1)*cos(i1)*(sin(omeg2)*cos(-Om2+Om1)-cos(omeg2)*&
     &cos(i2)*sin(-Om2+Om1))                                            
                                                                        
    dMMdOm2                                                           &
     &= cos(omeg1)*(-sin(omeg2)*sin(-Om2+Om1)-cos(omeg2)*cos(i2)*cos(   &
     &-Om2+Om1))+sin(omeg1)*cos(i1)*(-sin(omeg2)*cos(-Om2+Om1)+cos(omeg2&
     &)*cos(i2)*sin(-Om2+Om1))                                          
                                                                        
    dMMdomeg1                                                         &
     &= -sin(omeg1)*(-sin(omeg2)*cos(-Om2+Om1)+cos(omeg2)*cos(i2)*sin   &
     &(-Om2+Om1))+cos(omeg1)*cos(i1)*(sin(omeg2)*sin(-Om2+Om1)+cos(omeg2&
     &)*cos(i2)*cos(-Om2+Om1))+cos(omeg1)*sin(i1)*cos(omeg2)*sin(i2)    
                                                                        
    dMMdomeg2                                                         &
     &= cos(omeg1)*(-cos(omeg2)*cos(-Om2+Om1)-sin(omeg2)*cos(i2)*sin(   &
     &-Om2+Om1))+sin(omeg1)*cos(i1)*(cos(omeg2)*sin(-Om2+Om1)-sin(omeg2)&
     &*cos(i2)*cos(-Om2+Om1))-sin(omeg1)*sin(i1)*sin(omeg2)*sin(i2)     
                                                                        
    dNNdi1                                                            &
     &= -cos(omeg1)*sin(i1)*(sin(omeg2)*sin(-Om2+Om1)+cos(omeg2)*cos(   &
     &i2)*cos(-Om2+Om1))+cos(omeg1)*cos(i1)*cos(omeg2)*sin(i2)          
                                                                        
    dNNdi2                                                            &
     &= sin(omeg1)*cos(omeg2)*sin(i2)*sin(-Om2+Om1)-cos(omeg1)*cos(i1   &
     &)*cos(omeg2)*sin(i2)*cos(-Om2+Om1)+cos(omeg1)*sin(i1)*cos(omeg2)*c&
     &os(i2)                                                            
                                                                        
    dNNdOm1                                                           &
     &= -sin(omeg1)*(sin(omeg2)*sin(-Om2+Om1)+cos(omeg2)*cos(i2)*cos(   &
     &-Om2+Om1))+cos(omeg1)*cos(i1)*(sin(omeg2)*cos(-Om2+Om1)-cos(omeg2)&
     &*cos(i2)*sin(-Om2+Om1))                                           
                                                                        
    dNNdOm2                                                           &
     &= -sin(omeg1)*(-sin(omeg2)*sin(-Om2+Om1)-cos(omeg2)*cos(i2)*cos   &
     &(-Om2+Om1))+cos(omeg1)*cos(i1)*(-sin(omeg2)*cos(-Om2+Om1)+cos(omeg&
     &2)*cos(i2)*sin(-Om2+Om1))                                         
                                                                        
    dNNdomeg1                                                         &
     &= -cos(omeg1)*(-sin(omeg2)*cos(-Om2+Om1)+cos(omeg2)*cos(i2)*sin   &
     &(-Om2+Om1))-sin(omeg1)*cos(i1)*(sin(omeg2)*sin(-Om2+Om1)+cos(omeg2&
     &)*cos(i2)*cos(-Om2+Om1))-sin(omeg1)*sin(i1)*cos(omeg2)*sin(i2)    
                                                                        
    dNNdomeg2                                                         &
     &= -sin(omeg1)*(-cos(omeg2)*cos(-Om2+Om1)-sin(omeg2)*cos(i2)*sin   &
     &(-Om2+Om1))+cos(omeg1)*cos(i1)*(cos(omeg2)*sin(-Om2+Om1)-sin(omeg2&
     &)*cos(i2)*cos(-Om2+Om1))-cos(omeg1)*sin(i1)*sin(omeg2)*sin(i2)    

!======================================================================
! Derivatives of gradd2f1 e gradd2f2 w.r.t. the elements
! (by formuledetH2com.f90)
!======================================================================

      dgradd2f1_dKK = 2.D0/(1.D0+e1*cos(f1))**2.D0*q1*(1.D0+e1)*sin(f1)*&
     &q2*(1.D0+e2)/(1.D0+e2*cos(f2))*cos(f2)                            
                                                                        
      dgradd2f1_dLL = -2.D0/(1.D0+e1*cos(f1))*(e1*q1*(1.D0+e1)/(1.D0+e1*&
     &cos(f1))+q1*(1.D0+e1)/(1.D0+e1*cos(f1))*cos(f1))*q2*(1.D0+e2)/(1.D&
     &0+e2*cos(f2))*cos(f2)                                             
                                                                        
      dgradd2f1_dMM = 2.D0/(1.D0+e1*cos(f1))**2.D0*q1*(1.D0+e1)*sin(f1)*&
     &q2*(1.D0+e2)/(1.D0+e2*cos(f2))*sin(f2)                            
                                                                        
      dgradd2f1_dNN = -2.D0/(1.D0+e1*cos(f1))*(e1*q1*(1.D0+e1)/(1.D0+e1*&
     &cos(f1))+q1*(1.D0+e1)/(1.D0+e1*cos(f1))*cos(f1))*q2*(1.D0+e2)/(1.D&
     &0+e2*cos(f2))*sin(f2)                                             
                                                                        
      dgradd2f2_dKK = 2.D0/(1.D0+e2*cos(f2))**2.D0*q2*(1.D0+e2)*sin(f2)*&
     &q1*(1.D0+e1)/(1.D0+e1*cos(f1))*cos(f1)                            
                                                                        
      dgradd2f2_dLL = 2.D0/(1.D0+e2*cos(f2))**2.D0*q1*(1.D0+e1)/(1.D0+e1&
     &*cos(f1))*sin(f1)*q2*(1.D0+e2)*sin(f2)                            
                                                                        
      dgradd2f2_dMM = -2.D0/(1.D0+e2*cos(f2))*(e2*q2*(1.D0+e2)/(1.D0+e2*&
     &cos(f2))+q2*(1.D0+e2)/(1.D0+e2*cos(f2))*cos(f2))*q1*(1.D0+e1)/(1.D&
     &0+e1*cos(f1))*cos(f1)                                             
                                                                        
      dgradd2f2_dNN = -2.D0/(1.D0+e2*cos(f2))*(e2*q2*(1.D0+e2)/(1.D0+e2*&
     &cos(f2))+q2*(1.D0+e2)/(1.D0+e2*cos(f2))*cos(f2))*q1*(1.D0+e1)/(1.D&
     &0+e1*cos(f1))*sin(f1)                                             
                                                                        
      dgradd2f1_dq1 = 2.D0/(1.D0+e1*cos(f1))*(2.D0*e1*q1*(1.D0+e1)**2.D0&
     &/(1.D0+e1*cos(f1))**2.D0*sin(f1)+(1.D0+e1)/(1.D0+e1*cos(f1))*sin(f&
     &1)*(KK*q2*(1.D0+e2)/(1.D0+e2*cos(f2))*cos(f2)+MM*q2*(1.D0+e2)/(1.D&
     &0+e2*cos(f2))*sin(f2))-(e1*(1.D0+e1)/(1.D0+e1*cos(f1))+(1.D0+e1)/(&
     &1.D0+e1*cos(f1))*cos(f1))*(LL*q2*(1.D0+e2)/(1.D0+e2*cos(f2))*cos(f&
     &2)+NN*q2*(1.D0+e2)/(1.D0+e2*cos(f2))*sin(f2)))                    
                                                                        
      s1 = -2.D0/(1.D0+e1*cos(f1))**2.D0*(e1*q1**2*(1.D0+e1)**2.D0/(1.D0&
     &+e1*cos(f1))**2.D0*sin(f1)+q1*(1.D0+e1)/(1.D0+e1*cos(f1))*sin(f1)*&
     &(KK*q2*(1.D0+e2)/(1.D0+e2*cos(f2))*cos(f2)+MM*q2*(1.D0+e2)/(1.D0+e&
     &2*cos(f2))*sin(f2))-(e1*q1*(1.D0+e1)/(1.D0+e1*cos(f1))+q1*(1.D0+e1&
     &)/(1.D0+e1*cos(f1))*cos(f1))*(LL*q2*(1.D0+e2)/(1.D0+e2*cos(f2))*co&
     &s(f2)+NN*q2*(1.D0+e2)/(1.D0+e2*cos(f2))*sin(f2)))*cos(f1)         
      s3 = 2.D0 
      s5 = 1.D0/(1.D0+e1*cos(f1)) 
      s7 = q1**2*(1.D0+e1)**2.D0/(1.D0+e1*cos(f1))**2.D0*sin(f1)+2.D0*e1&
     &*q1**2*(1.D0+e1)/(1.D0+e1*cos(f1))**2.D0*sin(f1)-2.D0*e1*q1**2*(1.&
     &D0+e1)**2.D0/(1.D0+e1*cos(f1))**3.D0*sin(f1)*cos(f1)              
      s6 = s7+q1/(1.D0+e1*cos(f1))*sin(f1)*(KK*q2*(1.D0+e2)/(1.D0+e2*cos&
     &(f2))*cos(f2)+MM*q2*(1.D0+e2)/(1.D0+e2*cos(f2))*sin(f2))-q1*(1.D0+&
     &e1)/(1.D0+e1*cos(f1))**2.D0*sin(f1)*(KK*q2*(1.D0+e2)/(1.D0+e2*cos(&
     &f2))*cos(f2)+MM*q2*(1.D0+e2)/(1.D0+e2*cos(f2))*sin(f2))*cos(f1)-(q&
     &1*(1.D0+e1)/(1.D0+e1*cos(f1))+e1*q1/(1.D0+e1*cos(f1))-e1*q1*(1.D0+&
     &e1)/(1.D0+e1*cos(f1))**2.D0*cos(f1)+q1/(1.D0+e1*cos(f1))*cos(f1)-q&
     &1*(1.D0+e1)/(1.D0+e1*cos(f1))**2.D0*cos(f1)**2.D0)*(LL*q2*(1.D0+e2&
     &)/(1.D0+e2*cos(f2))*cos(f2)+NN*q2*(1.D0+e2)/(1.D0+e2*cos(f2))*sin(&
     &f2))                                                              
      s4 = s5*s6 
      s2 = s3*s4 
      dgradd2f1_de1 = s1+s2 
                                                                        
      dgradd2f1_di1 = 0.D0 
                                                                        
      dgradd2f1_dOm1 = 0.D0 
                                                                        
      dgradd2f1_domeg1 = 0.D0 
                                                                        
      dgradd2f1_dq2 = 2.D0/(1.D0+e1*cos(f1))*(q1*(1.D0+e1)/(1.D0+e1*cos(&
     &f1))*sin(f1)*(KK*(1.D0+e2)/(1.D0+e2*cos(f2))*cos(f2)+MM*(1.D0+e2)/&
     &(1.D0+e2*cos(f2))*sin(f2))-(e1*q1*(1.D0+e1)/(1.D0+e1*cos(f1))+q1*(&
     &1.D0+e1)/(1.D0+e1*cos(f1))*cos(f1))*(LL*(1.D0+e2)/(1.D0+e2*cos(f2)&
     &)*cos(f2)+NN*(1.D0+e2)/(1.D0+e2*cos(f2))*sin(f2)))                
                                                                        
      dgradd2f1_de2 = 2.D0/(1.D0+e1*cos(f1))*(q1*(1.D0+e1)/(1.D0+e1*cos(&
     &f1))*sin(f1)*(KK*q2/(1.D0+e2*cos(f2))*cos(f2)-KK*q2*(1.D0+e2)/(1.D&
     &0+e2*cos(f2))**2.D0*cos(f2)**2.D0+MM*q2/(1.D0+e2*cos(f2))*sin(f2)-&
     &MM*q2*(1.D0+e2)/(1.D0+e2*cos(f2))**2.D0*sin(f2)*cos(f2))-(e1*q1*(1&
     &.D0+e1)/(1.D0+e1*cos(f1))+q1*(1.D0+e1)/(1.D0+e1*cos(f1))*cos(f1))*&
     &(LL*q2/(1.D0+e2*cos(f2))*cos(f2)-LL*q2*(1.D0+e2)/(1.D0+e2*cos(f2))&
     &**2.D0*cos(f2)**2.D0+NN*q2/(1.D0+e2*cos(f2))*sin(f2)-NN*q2*(1.D0+e&
     &2)/(1.D0+e2*cos(f2))**2.D0*sin(f2)*cos(f2)))                      
                                                                        
      dgradd2f1_di2 = 0.D0 
                                                                        
      dgradd2f1_dOm2 = 0.D0 
                                                                        
      dgradd2f1_domeg2 = 0.D0 
                                                                        
      dgradd2f2_dq1 = 2.D0/(1.D0+e2*cos(f2))*(q2*(1.D0+e2)/(1.D0+e2*cos(&
     &f2))*sin(f2)*(KK*(1.D0+e1)/(1.D0+e1*cos(f1))*cos(f1)+LL*(1.D0+e1)/&
     &(1.D0+e1*cos(f1))*sin(f1))-(e2*q2*(1.D0+e2)/(1.D0+e2*cos(f2))+q2*(&
     &1.D0+e2)/(1.D0+e2*cos(f2))*cos(f2))*(MM*(1.D0+e1)/(1.D0+e1*cos(f1)&
     &)*cos(f1)+NN*(1.D0+e1)/(1.D0+e1*cos(f1))*sin(f1)))                
                                                                        
      dgradd2f2_de1 = 2.D0/(1.D0+e2*cos(f2))*(q2*(1.D0+e2)/(1.D0+e2*cos(&
     &f2))*sin(f2)*(KK*q1/(1.D0+e1*cos(f1))*cos(f1)-KK*q1*(1.D0+e1)/(1.D&
     &0+e1*cos(f1))**2.D0*cos(f1)**2.D0+LL*q1/(1.D0+e1*cos(f1))*sin(f1)-&
     &LL*q1*(1.D0+e1)/(1.D0+e1*cos(f1))**2.D0*sin(f1)*cos(f1))-(e2*q2*(1&
     &.D0+e2)/(1.D0+e2*cos(f2))+q2*(1.D0+e2)/(1.D0+e2*cos(f2))*cos(f2))*&
     &(MM*q1/(1.D0+e1*cos(f1))*cos(f1)-MM*q1*(1.D0+e1)/(1.D0+e1*cos(f1))&
     &**2.D0*cos(f1)**2.D0+NN*q1/(1.D0+e1*cos(f1))*sin(f1)-NN*q1*(1.D0+e&
     &1)/(1.D0+e1*cos(f1))**2.D0*sin(f1)*cos(f1)))                      
                                                                        
      dgradd2f2_di1 = 0.D0 
                                                                        
      dgradd2f2_dOm1 = 0.D0 
                                                                        
      dgradd2f2_domeg1 = 0.D0 
                                                                        
      dgradd2f2_dq2 = 2.D0/(1.D0+e2*cos(f2))*(2.D0*e2*q2*(1.D0+e2)**2.D0&
     &/(1.D0+e2*cos(f2))**2.D0*sin(f2)+(1.D0+e2)/(1.D0+e2*cos(f2))*sin(f&
     &2)*(KK*q1*(1.D0+e1)/(1.D0+e1*cos(f1))*cos(f1)+LL*q1*(1.D0+e1)/(1.D&
     &0+e1*cos(f1))*sin(f1))-(e2*(1.D0+e2)/(1.D0+e2*cos(f2))+(1.D0+e2)/(&
     &1.D0+e2*cos(f2))*cos(f2))*(MM*q1*(1.D0+e1)/(1.D0+e1*cos(f1))*cos(f&
     &1)+NN*q1*(1.D0+e1)/(1.D0+e1*cos(f1))*sin(f1)))                    
                                                                        
      s1 = -2.D0/(1.D0+e2*cos(f2))**2.D0*(e2*q2**2*(1.D0+e2)**2.D0/(1.D0&
     &+e2*cos(f2))**2.D0*sin(f2)+q2*(1.D0+e2)/(1.D0+e2*cos(f2))*sin(f2)*&
     &(KK*q1*(1.D0+e1)/(1.D0+e1*cos(f1))*cos(f1)+LL*q1*(1.D0+e1)/(1.D0+e&
     &1*cos(f1))*sin(f1))-(e2*q2*(1.D0+e2)/(1.D0+e2*cos(f2))+q2*(1.D0+e2&
     &)/(1.D0+e2*cos(f2))*cos(f2))*(MM*q1*(1.D0+e1)/(1.D0+e1*cos(f1))*co&
     &s(f1)+NN*q1*(1.D0+e1)/(1.D0+e1*cos(f1))*sin(f1)))*cos(f2)         
      s3 = 2.D0 
      s5 = 1.D0/(1.D0+e2*cos(f2)) 
      s7 = q2**2*(1.D0+e2)**2.D0/(1.D0+e2*cos(f2))**2.D0*sin(f2)+2.D0*e2&
     &*q2**2*(1.D0+e2)/(1.D0+e2*cos(f2))**2.D0*sin(f2)-2.D0*e2*q2**2*(1.&
     &D0+e2)**2.D0/(1.D0+e2*cos(f2))**3.D0*sin(f2)*cos(f2)              
      s6 = s7+q2/(1.D0+e2*cos(f2))*sin(f2)*(KK*q1*(1.D0+e1)/(1.D0+e1*cos&
     &(f1))*cos(f1)+LL*q1*(1.D0+e1)/(1.D0+e1*cos(f1))*sin(f1))-q2*(1.D0+&
     &e2)/(1.D0+e2*cos(f2))**2.D0*sin(f2)*(KK*q1*(1.D0+e1)/(1.D0+e1*cos(&
     &f1))*cos(f1)+LL*q1*(1.D0+e1)/(1.D0+e1*cos(f1))*sin(f1))*cos(f2)-(q&
     &2*(1.D0+e2)/(1.D0+e2*cos(f2))+e2*q2/(1.D0+e2*cos(f2))-e2*q2*(1.D0+&
     &e2)/(1.D0+e2*cos(f2))**2.D0*cos(f2)+q2/(1.D0+e2*cos(f2))*cos(f2)-q&
     &2*(1.D0+e2)/(1.D0+e2*cos(f2))**2.D0*cos(f2)**2.D0)*(MM*q1*(1.D0+e1&
     &)/(1.D0+e1*cos(f1))*cos(f1)+NN*q1*(1.D0+e1)/(1.D0+e1*cos(f1))*sin(&
     &f1))                                                              
      s4 = s5*s6 
      s2 = s3*s4 
      dgradd2f2_de2 = s1+s2 
                                                                        
      dgradd2f2_di2 = 0.D0 
                                                                        
      dgradd2f2_dOm2 = 0.D0 
                                                                        
      dgradd2f2_domeg2 = 0.D0 
                                                                        
! ====== chain rule ======
      dgradd2f1_di1 = dgradd2f1_dKK*dKKdi1+dgradd2f1_dLL*dLLdi1+ &
           & dgradd2f1_dMM*dMMdi1+dgradd2f1_dNN*dNNdi1

      dgradd2f1_dOm1 = dgradd2f1_dKK*dKKdOm1+dgradd2f1_dLL*dLLdOm1+ &
           & dgradd2f1_dMM*dMMdOm1+dgradd2f1_dNN*dNNdOm1

      dgradd2f1_domeg1 = dgradd2f1_dKK*dKKdomeg1+dgradd2f1_dLL*dLLdomeg1+ &
           & dgradd2f1_dMM*dMMdomeg1+dgradd2f1_dNN*dNNdomeg1
                    
      dgradd2f1_di2 = dgradd2f1_dKK*dKKdi2+dgradd2f1_dLL*dLLdi2+ &
           & dgradd2f1_dMM*dMMdi2+dgradd2f1_dNN*dNNdi2

      dgradd2f1_dOm2 = dgradd2f1_dKK*dKKdOm2+dgradd2f1_dLL*dLLdOm2+ &
           & dgradd2f1_dMM*dMMdOm2+dgradd2f1_dNN*dNNdOm2

      dgradd2f1_domeg2 = dgradd2f1_dKK*dKKdomeg2+dgradd2f1_dLL*dLLdomeg2+ &
           & dgradd2f1_dMM*dMMdomeg2+dgradd2f1_dNN*dNNdomeg2

      dgradd2f2_di1 = dgradd2f2_dKK*dKKdi1+dgradd2f2_dLL*dLLdi1+ &
           & dgradd2f2_dMM*dMMdi1+dgradd2f2_dNN*dNNdi1

      dgradd2f2_dOm1 = dgradd2f2_dKK*dKKdOm1+dgradd2f2_dLL*dLLdOm1+ &
           & dgradd2f2_dMM*dMMdOm1+dgradd2f2_dNN*dNNdOm1

      dgradd2f2_domeg1 = dgradd2f2_dKK*dKKdomeg1+dgradd2f2_dLL*dLLdomeg1+ &
           & dgradd2f2_dMM*dMMdomeg1+dgradd2f2_dNN*dNNdomeg1
                    
      dgradd2f2_di2 = dgradd2f2_dKK*dKKdi2+dgradd2f2_dLL*dLLdi2+ &
           & dgradd2f2_dMM*dMMdi2+dgradd2f2_dNN*dNNdi2

      dgradd2f2_dOm2 = dgradd2f2_dKK*dKKdOm2+dgradd2f2_dLL*dLLdOm2+ &
           & dgradd2f2_dMM*dMMdOm2+dgradd2f2_dNN*dNNdOm2

      dgradd2f2_domeg2 = dgradd2f2_dKK*dKKdomeg2+dgradd2f2_dLL*dLLdomeg2+ &
           & dgradd2f2_dMM*dMMdomeg2+dgradd2f2_dNN*dNNdomeg2

!======================================================================
! Derivatives of f1 e f2 w.r.t. the elements
!======================================================================
    df1_dq1=-(invH(1,1)*dgradd2f1_dq1+invH(1,2)*dgradd2f2_dq1)
    df2_dq1=-(invH(2,1)*dgradd2f1_dq1+invH(2,2)*dgradd2f2_dq1)

    df1_de1=-(invH(1,1)*dgradd2f1_de1+invH(1,2)*dgradd2f2_de1)
    df2_de1=-(invH(2,1)*dgradd2f1_de1+invH(2,2)*dgradd2f2_de1)

    df1_di1=-(invH(1,1)*dgradd2f1_di1+invH(1,2)*dgradd2f2_di1)
    df2_di1=-(invH(2,1)*dgradd2f1_di1+invH(2,2)*dgradd2f2_di1)

    df1_dOm1=-(invH(1,1)*dgradd2f1_dOm1+invH(1,2)*dgradd2f2_dOm1)
    df2_dOm1=-(invH(2,1)*dgradd2f1_dOm1+invH(2,2)*dgradd2f2_dOm1)

    df1_domeg1=-(invH(1,1)*dgradd2f1_domeg1+invH(1,2)*dgradd2f2_domeg1)
    df2_domeg1=-(invH(2,1)*dgradd2f1_domeg1+invH(2,2)*dgradd2f2_domeg1)

    df1_dq2=-(invH(1,1)*dgradd2f1_dq2+invH(1,2)*dgradd2f2_dq2)
    df2_dq2=-(invH(2,1)*dgradd2f1_dq2+invH(2,2)*dgradd2f2_dq2)

    df1_de2=-(invH(1,1)*dgradd2f1_de2+invH(1,2)*dgradd2f2_de2)
    df2_de2=-(invH(2,1)*dgradd2f1_de2+invH(2,2)*dgradd2f2_de2)

    df1_di2=-(invH(1,1)*dgradd2f1_di2+invH(1,2)*dgradd2f2_di2)
    df2_di2=-(invH(2,1)*dgradd2f1_di2+invH(2,2)*dgradd2f2_di2)

    df1_dOm2=-(invH(1,1)*dgradd2f1_dOm2+invH(1,2)*dgradd2f2_dOm2)
    df2_dOm2=-(invH(2,1)*dgradd2f1_dOm2+invH(2,2)*dgradd2f2_dOm2)

    df1_domeg2=-(invH(1,1)*dgradd2f1_domeg2+invH(1,2)*dgradd2f2_domeg2)
    df2_domeg2=-(invH(2,1)*dgradd2f1_domeg2+invH(2,2)*dgradd2f2_domeg2)

!=========================================================================
! Derivatives of elements of matrix H (by formuledetHcom.f90)
!=========================================================================
! Derivatives w.r.t. f1 and f2
      s1 = e1**2/(1.D0+e1*cos(f1))**2.D0*sin(f1)**2.D0*gradd2f1+e1/(1.D0&
     &+e1*cos(f1))*cos(f1)*gradd2f1                                     
      s2 = s1 
      s4 = 2.D0/(1.D0+e1*cos(f1))**2.D0*(e1**2*q1**2*(1.D0+e1)**2.D0/(1.&
     &D0+e1*cos(f1))**3.D0*sin(f1)**2.D0+(e1*q1*(1.D0+e1)/(1.D0+e1*cos(f&
     &1))**2.D0*sin(f1)**2.D0+q1*(1.D0+e1)/(1.D0+e1*cos(f1))*cos(f1))*(e&
     &1*q1*(1.D0+e1)/(1.D0+e1*cos(f1))+KK*q2*(1.D0+e2)/(1.D0+e2*cos(f2))&
     &*cos(f2)+MM*q2*(1.D0+e2)/(1.D0+e2*cos(f2))*sin(f2))-(e1**2*q1*(1.D&
     &0+e1)/(1.D0+e1*cos(f1))**2.D0*sin(f1)+e1*q1*(1.D0+e1)/(1.D0+e1*cos&
     &(f1))**2.D0*cos(f1)*sin(f1)-q1*(1.D0+e1)/(1.D0+e1*cos(f1))*sin(f1)&
     &)*(LL*q2*(1.D0+e2)/(1.D0+e2*cos(f2))*cos(f2)+NN*q2*(1.D0+e2)/(1.D0&
     &+e2*cos(f2))*sin(f2)))*e1*sin(f1)                                 
      s6 = 2.D0 
      s8 = 1.D0/(1.D0+e1*cos(f1)) 
      s10 = 3.D0*e1**3*q1**2*(1.D0+e1)**2.D0/(1.D0+e1*cos(f1))**4.D0*sin&
     &(f1)**3.D0+2.D0*e1**2*q1**2*(1.D0+e1)**2.D0/(1.D0+e1*cos(f1))**3.D&
     &0*sin(f1)*cos(f1)                                                 
      s11 = s10+(2.D0*e1**2*q1*(1.D0+e1)/(1.D0+e1*cos(f1))**3.D0*sin(f1)&
     &**3.D0+3.D0*e1*q1*(1.D0+e1)/(1.D0+e1*cos(f1))**2.D0*cos(f1)*sin(f1&
     &)-q1*(1.D0+e1)/(1.D0+e1*cos(f1))*sin(f1))*(e1*q1*(1.D0+e1)/(1.D0+e&
     &1*cos(f1))+KK*q2*(1.D0+e2)/(1.D0+e2*cos(f2))*cos(f2)+MM*q2*(1.D0+e&
     &2)/(1.D0+e2*cos(f2))*sin(f2))                                     
      s9 = s11+(e1*q1*(1.D0+e1)/(1.D0+e1*cos(f1))**2.D0*sin(f1)**2.D0+q1&
     &*(1.D0+e1)/(1.D0+e1*cos(f1))*cos(f1))*e1**2*q1*(1.D0+e1)/(1.D0+e1*&
     &cos(f1))**2.D0*sin(f1)-(2.D0*e1**3*q1*(1.D0+e1)/(1.D0+e1*cos(f1))*&
     &*3.D0*sin(f1)**2.D0+e1**2*q1*(1.D0+e1)/(1.D0+e1*cos(f1))**2.D0*cos&
     &(f1)+2.D0*e1**2*q1*(1.D0+e1)/(1.D0+e1*cos(f1))**3.D0*cos(f1)*sin(f&
     &1)**2.D0-2.D0*e1*q1*(1.D0+e1)/(1.D0+e1*cos(f1))**2.D0*sin(f1)**2.D&
     &0+e1*q1*(1.D0+e1)/(1.D0+e1*cos(f1))**2.D0*cos(f1)**2.D0-q1*(1.D0+e&
     &1)/(1.D0+e1*cos(f1))*cos(f1))*(LL*q2*(1.D0+e2)/(1.D0+e2*cos(f2))*c&
     &os(f2)+NN*q2*(1.D0+e2)/(1.D0+e2*cos(f2))*sin(f2))                 
      s7 = s8*s9 
      s5 = s6*s7 
      s3 = s4+s5 
      dH11_df1 = s2+s3 

      s1 = 2.D0 
      s3 = 1.D0/(1.D0+e1*cos(f1)) 
      s4 = (e1*q1*(1.D0+e1)/(1.D0+e1*cos(f1))**2.D0*sin(f1)**2.D0+q1*(1.&
     &D0+e1)/(1.D0+e1*cos(f1))*cos(f1))*(KK*q2*(1.D0+e2)/(1.D0+e2*cos(f2&
     &))**2.D0*cos(f2)*e2*sin(f2)-KK*q2*(1.D0+e2)/(1.D0+e2*cos(f2))*sin(&
     &f2)+MM*q2*(1.D0+e2)/(1.D0+e2*cos(f2))**2.D0*sin(f2)**2.D0*e2+MM*q2&
     &*(1.D0+e2)/(1.D0+e2*cos(f2))*cos(f2))-(e1**2*q1*(1.D0+e1)/(1.D0+e1&
     &*cos(f1))**2.D0*sin(f1)+e1*q1*(1.D0+e1)/(1.D0+e1*cos(f1))**2.D0*co&
     &s(f1)*sin(f1)-q1*(1.D0+e1)/(1.D0+e1*cos(f1))*sin(f1))*(LL*q2*(1.D0&
     &+e2)/(1.D0+e2*cos(f2))**2.D0*cos(f2)*e2*sin(f2)-LL*q2*(1.D0+e2)/(1&
     &.D0+e2*cos(f2))*sin(f2)+NN*q2*(1.D0+e2)/(1.D0+e2*cos(f2))**2.D0*si&
     &n(f2)**2.D0*e2+NN*q2*(1.D0+e2)/(1.D0+e2*cos(f2))*cos(f2))         
      s2 = s3*s4 
      dH11_df2 = s1*s2 
                                                                        
      s1 = 2.D0 
      s3 = 1.D0/(1.D0+e2*cos(f2)) 
      s4 = (e2*q2*(1.D0+e2)/(1.D0+e2*cos(f2))**2.D0*sin(f2)**2.D0+q2*(1.&
     &D0+e2)/(1.D0+e2*cos(f2))*cos(f2))*(KK*q1*(1.D0+e1)/(1.D0+e1*cos(f1&
     &))**2.D0*cos(f1)*e1*sin(f1)-KK*q1*(1.D0+e1)/(1.D0+e1*cos(f1))*sin(&
     &f1)+LL*q1*(1.D0+e1)/(1.D0+e1*cos(f1))**2.D0*sin(f1)**2.D0*e1+LL*q1&
     &*(1.D0+e1)/(1.D0+e1*cos(f1))*cos(f1))-(e2**2*q2*(1.D0+e2)/(1.D0+e2&
     &*cos(f2))**2.D0*sin(f2)+e2*q2*(1.D0+e2)/(1.D0+e2*cos(f2))**2.D0*co&
     &s(f2)*sin(f2)-q2*(1.D0+e2)/(1.D0+e2*cos(f2))*sin(f2))*(MM*q1*(1.D0&
     &+e1)/(1.D0+e1*cos(f1))**2.D0*cos(f1)*e1*sin(f1)-MM*q1*(1.D0+e1)/(1&
     &.D0+e1*cos(f1))*sin(f1)+NN*q1*(1.D0+e1)/(1.D0+e1*cos(f1))**2.D0*si&
     &n(f1)**2.D0*e1+NN*q1*(1.D0+e1)/(1.D0+e1*cos(f1))*cos(f1))         
      s2 = s3*s4 
      dH22_df1 = s1*s2 
                                                                        
      s1 = e2**2/(1.D0+e2*cos(f2))**2.D0*sin(f2)**2.D0*gradd2f2+e2/(1.D0&
     &+e2*cos(f2))*cos(f2)*gradd2f2                                     
      s2 = s1 
      s4 = 2.D0/(1.D0+e2*cos(f2))**2.D0*(e2**2*q2**2*(1.D0+e2)**2.D0/(1.&
     &D0+e2*cos(f2))**3.D0*sin(f2)**2.D0+(e2*q2*(1.D0+e2)/(1.D0+e2*cos(f&
     &2))**2.D0*sin(f2)**2.D0+q2*(1.D0+e2)/(1.D0+e2*cos(f2))*cos(f2))*(e&
     &2*q2*(1.D0+e2)/(1.D0+e2*cos(f2))+KK*q1*(1.D0+e1)/(1.D0+e1*cos(f1))&
     &*cos(f1)+LL*q1*(1.D0+e1)/(1.D0+e1*cos(f1))*sin(f1))-(e2**2*q2*(1.D&
     &0+e2)/(1.D0+e2*cos(f2))**2.D0*sin(f2)+e2*q2*(1.D0+e2)/(1.D0+e2*cos&
     &(f2))**2.D0*cos(f2)*sin(f2)-q2*(1.D0+e2)/(1.D0+e2*cos(f2))*sin(f2)&
     &)*(MM*q1*(1.D0+e1)/(1.D0+e1*cos(f1))*cos(f1)+NN*q1*(1.D0+e1)/(1.D0&
     &+e1*cos(f1))*sin(f1)))*e2*sin(f2)                                 
      s6 = 2.D0 
      s8 = 1.D0/(1.D0+e2*cos(f2)) 
      s10 = 3.D0*e2**3*q2**2*(1.D0+e2)**2.D0/(1.D0+e2*cos(f2))**4.D0*sin&
     &(f2)**3.D0+2.D0*e2**2*q2**2*(1.D0+e2)**2.D0/(1.D0+e2*cos(f2))**3.D&
     &0*sin(f2)*cos(f2)                                                 
      s11 = s10+(2.D0*e2**2*q2*(1.D0+e2)/(1.D0+e2*cos(f2))**3.D0*sin(f2)&
     &**3.D0+3.D0*e2*q2*(1.D0+e2)/(1.D0+e2*cos(f2))**2.D0*cos(f2)*sin(f2&
     &)-q2*(1.D0+e2)/(1.D0+e2*cos(f2))*sin(f2))*(e2*q2*(1.D0+e2)/(1.D0+e&
     &2*cos(f2))+KK*q1*(1.D0+e1)/(1.D0+e1*cos(f1))*cos(f1)+LL*q1*(1.D0+e&
     &1)/(1.D0+e1*cos(f1))*sin(f1))                                     
      s9 = s11+(e2*q2*(1.D0+e2)/(1.D0+e2*cos(f2))**2.D0*sin(f2)**2.D0+q2&
     &*(1.D0+e2)/(1.D0+e2*cos(f2))*cos(f2))*e2**2*q2*(1.D0+e2)/(1.D0+e2*&
     &cos(f2))**2.D0*sin(f2)-(2.D0*e2**3*q2*(1.D0+e2)/(1.D0+e2*cos(f2))*&
     &*3.D0*sin(f2)**2.D0+e2**2*q2*(1.D0+e2)/(1.D0+e2*cos(f2))**2.D0*cos&
     &(f2)+2.D0*e2**2*q2*(1.D0+e2)/(1.D0+e2*cos(f2))**3.D0*cos(f2)*sin(f&
     &2)**2.D0-2.D0*e2*q2*(1.D0+e2)/(1.D0+e2*cos(f2))**2.D0*sin(f2)**2.D&
     &0+e2*q2*(1.D0+e2)/(1.D0+e2*cos(f2))**2.D0*cos(f2)**2.D0-q2*(1.D0+e&
     &2)/(1.D0+e2*cos(f2))*cos(f2))*(MM*q1*(1.D0+e1)/(1.D0+e1*cos(f1))*c&
     &os(f1)+NN*q1*(1.D0+e1)/(1.D0+e1*cos(f1))*sin(f1))                 
      s7 = s8*s9 
      s5 = s6*s7 
      s3 = s4+s5 
      dH22_df2 = s2+s3 
                                                                        
      s1 = 2.D0/(1.D0+e1*cos(f1))**2.D0*(q1*(1.D0+e1)/(1.D0+e1*cos(f1))*&
     &sin(f1)*(KK*(e2*q2*(1.D0+e2)/(1.D0+e2*cos(f2))**2.D0*cos(f2)*sin(f&
     &2)-q2*(1.D0+e2)/(1.D0+e2*cos(f2))*sin(f2))+MM*(e2*q2*(1.D0+e2)/(1.&
     &D0+e2*cos(f2))**2.D0*sin(f2)**2.D0+q2*(1.D0+e2)/(1.D0+e2*cos(f2))*&
     &cos(f2)))-(e1*q1*(1.D0+e1)/(1.D0+e1*cos(f1))+q1*(1.D0+e1)/(1.D0+e1&
     &*cos(f1))*cos(f1))*(LL*(e2*q2*(1.D0+e2)/(1.D0+e2*cos(f2))**2.D0*co&
     &s(f2)*sin(f2)-q2*(1.D0+e2)/(1.D0+e2*cos(f2))*sin(f2))+NN*(e2*q2*(1&
     &.D0+e2)/(1.D0+e2*cos(f2))**2.D0*sin(f2)**2.D0+q2*(1.D0+e2)/(1.D0+e&
     &2*cos(f2))*cos(f2))))*e1*sin(f1)                                  
      s3 = 2.D0 
      s5 = 1.D0/(1.D0+e1*cos(f1)) 
      s7 = q1*(1.D0+e1)/(1.D0+e1*cos(f1))**2.D0*sin(f1)**2.D0*(KK*(e2*q2&
     &*(1.D0+e2)/(1.D0+e2*cos(f2))**2.D0*cos(f2)*sin(f2)-q2*(1.D0+e2)/(1&
     &.D0+e2*cos(f2))*sin(f2))+MM*(e2*q2*(1.D0+e2)/(1.D0+e2*cos(f2))**2.&
     &D0*sin(f2)**2.D0+q2*(1.D0+e2)/(1.D0+e2*cos(f2))*cos(f2)))*e1      
      s8 = q1*(1.D0+e1)/(1.D0+e1*cos(f1))*cos(f1)*(KK*(e2*q2*(1.D0+e2)/(&
     &1.D0+e2*cos(f2))**2.D0*cos(f2)*sin(f2)-q2*(1.D0+e2)/(1.D0+e2*cos(f&
     &2))*sin(f2))+MM*(e2*q2*(1.D0+e2)/(1.D0+e2*cos(f2))**2.D0*sin(f2)**&
     &2.D0+q2*(1.D0+e2)/(1.D0+e2*cos(f2))*cos(f2)))-(e1**2*q1*(1.D0+e1)/&
     &(1.D0+e1*cos(f1))**2.D0*sin(f1)+e1*q1*(1.D0+e1)/(1.D0+e1*cos(f1))*&
     &*2.D0*cos(f1)*sin(f1)-q1*(1.D0+e1)/(1.D0+e1*cos(f1))*sin(f1))*(LL*&
     &(e2*q2*(1.D0+e2)/(1.D0+e2*cos(f2))**2.D0*cos(f2)*sin(f2)-q2*(1.D0+&
     &e2)/(1.D0+e2*cos(f2))*sin(f2))+NN*(e2*q2*(1.D0+e2)/(1.D0+e2*cos(f2&
     &))**2.D0*sin(f2)**2.D0+q2*(1.D0+e2)/(1.D0+e2*cos(f2))*cos(f2)))   
      s6 = s7+s8 
      s4 = s5*s6 
      s2 = s3*s4 
      dH12_df1 = s1+s2 
                                                                        
      s1 = 2.D0 
      s3 = 1.D0/(1.D0+e1*cos(f1)) 
      s5 = q1*(1.D0+e1)/(1.D0+e1*cos(f1))*sin(f1)*(KK*(2.D0*e2**2*q2*(1.&
     &D0+e2)/(1.D0+e2*cos(f2))**3.D0*cos(f2)*sin(f2)**2.D0-2.D0*e2*q2*(1&
     &.D0+e2)/(1.D0+e2*cos(f2))**2.D0*sin(f2)**2.D0+e2*q2*(1.D0+e2)/(1.D&
     &0+e2*cos(f2))**2.D0*cos(f2)**2.D0-q2*(1.D0+e2)/(1.D0+e2*cos(f2))*c&
     &os(f2))+MM*(2.D0*e2**2*q2*(1.D0+e2)/(1.D0+e2*cos(f2))**3.D0*sin(f2&
     &)**3.D0+3.D0*e2*q2*(1.D0+e2)/(1.D0+e2*cos(f2))**2.D0*cos(f2)*sin(f&
     &2)-q2*(1.D0+e2)/(1.D0+e2*cos(f2))*sin(f2)))                       
      s6 = -(e1*q1*(1.D0+e1)/(1.D0+e1*cos(f1))+q1*(1.D0+e1)/(1.D0+e1*cos&
     &(f1))*cos(f1))*(LL*(2.D0*e2**2*q2*(1.D0+e2)/(1.D0+e2*cos(f2))**3.D&
     &0*cos(f2)*sin(f2)**2.D0-2.D0*e2*q2*(1.D0+e2)/(1.D0+e2*cos(f2))**2.&
     &D0*sin(f2)**2.D0+e2*q2*(1.D0+e2)/(1.D0+e2*cos(f2))**2.D0*cos(f2)**&
     &2.D0-q2*(1.D0+e2)/(1.D0+e2*cos(f2))*cos(f2))+NN*(2.D0*e2**2*q2*(1.&
     &D0+e2)/(1.D0+e2*cos(f2))**3.D0*sin(f2)**3.D0+3.D0*e2*q2*(1.D0+e2)/&
     &(1.D0+e2*cos(f2))**2.D0*cos(f2)*sin(f2)-q2*(1.D0+e2)/(1.D0+e2*cos(&
     &f2))*sin(f2)))                                                    
      s4 = s5+s6 
      s2 = s3*s4 
      dH12_df2 = s1*s2 
                                                                        
! Derivatives w.r.t. KK,LL,MM,NN
      dH11_dKK = 2.D0/(1.D0+e1*cos(f1))*(e1*q1*(1.D0+e1)/(1.D0+e1*cos(f1&
     &))**2.D0*sin(f1)**2.D0+q1*(1.D0+e1)/(1.D0+e1*cos(f1))*cos(f1))*q2*&
     &(1.D0+e2)/(1.D0+e2*cos(f2))*cos(f2)                               
                                                                        
      dH11_dLL = -2.D0/(1.D0+e1*cos(f1))*(e1**2*q1*(1.D0+e1)/(1.D0+e1*co&
     &s(f1))**2.D0*sin(f1)+e1*q1*(1.D0+e1)/(1.D0+e1*cos(f1))**2.D0*cos(f&
     &1)*sin(f1)-q1*(1.D0+e1)/(1.D0+e1*cos(f1))*sin(f1))*q2*(1.D0+e2)/(1&
     &.D0+e2*cos(f2))*cos(f2)                                           
                                                                        
      dH11_dMM = 2.D0/(1.D0+e1*cos(f1))*(e1*q1*(1.D0+e1)/(1.D0+e1*cos(f1&
     &))**2.D0*sin(f1)**2.D0+q1*(1.D0+e1)/(1.D0+e1*cos(f1))*cos(f1))*q2*&
     &(1.D0+e2)/(1.D0+e2*cos(f2))*sin(f2)                               
                                                                        
      dH11_dNN = -2.D0/(1.D0+e1*cos(f1))*(e1**2*q1*(1.D0+e1)/(1.D0+e1*co&
     &s(f1))**2.D0*sin(f1)+e1*q1*(1.D0+e1)/(1.D0+e1*cos(f1))**2.D0*cos(f&
     &1)*sin(f1)-q1*(1.D0+e1)/(1.D0+e1*cos(f1))*sin(f1))*q2*(1.D0+e2)/(1&
     &.D0+e2*cos(f2))*sin(f2)                                           
                                                                        
      dH22_dKK = 2.D0/(1.D0+e2*cos(f2))*(e2*q2*(1.D0+e2)/(1.D0+e2*cos(f2&
     &))**2.D0*sin(f2)**2.D0+q2*(1.D0+e2)/(1.D0+e2*cos(f2))*cos(f2))*q1*&
     &(1.D0+e1)/(1.D0+e1*cos(f1))*cos(f1)                               
                                                                        
      dH22_dLL = 2.D0/(1.D0+e2*cos(f2))*(e2*q2*(1.D0+e2)/(1.D0+e2*cos(f2&
     &))**2.D0*sin(f2)**2.D0+q2*(1.D0+e2)/(1.D0+e2*cos(f2))*cos(f2))*q1*&
     &(1.D0+e1)/(1.D0+e1*cos(f1))*sin(f1)                               
                                                                        
      dH22_dMM = -2.D0/(1.D0+e2*cos(f2))*(e2**2*q2*(1.D0+e2)/(1.D0+e2*co&
     &s(f2))**2.D0*sin(f2)+e2*q2*(1.D0+e2)/(1.D0+e2*cos(f2))**2.D0*cos(f&
     &2)*sin(f2)-q2*(1.D0+e2)/(1.D0+e2*cos(f2))*sin(f2))*q1*(1.D0+e1)/(1&
     &.D0+e1*cos(f1))*cos(f1)                                           
                                                                        
      dH22_dNN = -2.D0/(1.D0+e2*cos(f2))*(e2**2*q2*(1.D0+e2)/(1.D0+e2*co&
     &s(f2))**2.D0*sin(f2)+e2*q2*(1.D0+e2)/(1.D0+e2*cos(f2))**2.D0*cos(f&
     &2)*sin(f2)-q2*(1.D0+e2)/(1.D0+e2*cos(f2))*sin(f2))*q1*(1.D0+e1)/(1&
     &.D0+e1*cos(f1))*sin(f1)                                           
                                                                        
      dH12_dKK = 2.D0/(1.D0+e1*cos(f1))**2.D0*q1*(1.D0+e1)*sin(f1)*(e2*q&
     &2*(1.D0+e2)/(1.D0+e2*cos(f2))**2.D0*cos(f2)*sin(f2)-q2*(1.D0+e2)/(&
     &1.D0+e2*cos(f2))*sin(f2))                                         
                                                                        
      dH12_dLL = -2.D0/(1.D0+e1*cos(f1))*(e1*q1*(1.D0+e1)/(1.D0+e1*cos(f&
     &1))+q1*(1.D0+e1)/(1.D0+e1*cos(f1))*cos(f1))*(e2*q2*(1.D0+e2)/(1.D0&
     &+e2*cos(f2))**2.D0*cos(f2)*sin(f2)-q2*(1.D0+e2)/(1.D0+e2*cos(f2))*&
     &sin(f2))                                                          
                                                                        
      dH12_dMM = 2.D0/(1.D0+e1*cos(f1))**2.D0*(e2*q2*(1.D0+e2)/(1.D0+e2*&
     &cos(f2))**2.D0*sin(f2)**2.D0+q2*(1.D0+e2)/(1.D0+e2*cos(f2))*cos(f2&
     &))*q1*(1.D0+e1)*sin(f1)                                           
                                                                        
      dH12_dNN = -2.D0/(1.D0+e1*cos(f1))*(e1*q1*(1.D0+e1)/(1.D0+e1*cos(f&
     &1))+q1*(1.D0+e1)/(1.D0+e1*cos(f1))*cos(f1))*(e2*q2*(1.D0+e2)/(1.D0&
     &+e2*cos(f2))**2.D0*sin(f2)**2.D0+q2*(1.D0+e2)/(1.D0+e2*cos(f2))*co&
     &s(f2))                                                            
                                                                        
! Derivatives w.r.t. the elements
      s1 = 2.D0 
      s2 = 1.D0/(1.D0+e1*cos(f1))*(2.D0*e1**2*q1*(1.D0+e1)**2.D0/(1.D0+e&
     &1*cos(f1))**3.D0*sin(f1)**2.D0+(e1*(1.D0+e1)/(1.D0+e1*cos(f1))**2.&
     &D0*sin(f1)**2.D0+(1.D0+e1)/(1.D0+e1*cos(f1))*cos(f1))*(e1*q1*(1.D0&
     &+e1)/(1.D0+e1*cos(f1))+KK*q2*(1.D0+e2)/(1.D0+e2*cos(f2))*cos(f2)+M&
     &M*q2*(1.D0+e2)/(1.D0+e2*cos(f2))*sin(f2))+(e1*q1*(1.D0+e1)/(1.D0+e&
     &1*cos(f1))**2.D0*sin(f1)**2.D0+q1*(1.D0+e1)/(1.D0+e1*cos(f1))*cos(&
     &f1))*e1*(1.D0+e1)/(1.D0+e1*cos(f1))-(e1**2*(1.D0+e1)/(1.D0+e1*cos(&
     &f1))**2.D0*sin(f1)+e1*(1.D0+e1)/(1.D0+e1*cos(f1))**2.D0*cos(f1)*si&
     &n(f1)-(1.D0+e1)/(1.D0+e1*cos(f1))*sin(f1))*(LL*q2*(1.D0+e2)/(1.D0+&
     &e2*cos(f2))*cos(f2)+NN*q2*(1.D0+e2)/(1.D0+e2*cos(f2))*sin(f2)))   
      dH11_dq1 = s1*s2 
                                                                        
      s1 = 1.D0/(1.D0+e1*cos(f1))*sin(f1)*gradd2f1-e1/(1.D0+e1*cos(f1))*&
     &*2.D0*sin(f1)*gradd2f1*cos(f1)                                    
      s2 = s1 
      s4 = -2.D0/(1.D0+e1*cos(f1))**2.D0*(e1**2*q1**2*(1.D0+e1)**2.D0/(1&
     &.D0+e1*cos(f1))**3.D0*sin(f1)**2.D0+(e1*q1*(1.D0+e1)/(1.D0+e1*cos(&
     &f1))**2.D0*sin(f1)**2.D0+q1*(1.D0+e1)/(1.D0+e1*cos(f1))*cos(f1))*(&
     &e1*q1*(1.D0+e1)/(1.D0+e1*cos(f1))+KK*q2*(1.D0+e2)/(1.D0+e2*cos(f2)&
     &)*cos(f2)+MM*q2*(1.D0+e2)/(1.D0+e2*cos(f2))*sin(f2))-(e1**2*q1*(1.&
     &D0+e1)/(1.D0+e1*cos(f1))**2.D0*sin(f1)+e1*q1*(1.D0+e1)/(1.D0+e1*co&
     &s(f1))**2.D0*cos(f1)*sin(f1)-q1*(1.D0+e1)/(1.D0+e1*cos(f1))*sin(f1&
     &))*(LL*q2*(1.D0+e2)/(1.D0+e2*cos(f2))*cos(f2)+NN*q2*(1.D0+e2)/(1.D&
     &0+e2*cos(f2))*sin(f2)))*cos(f1)                                   
      s6 = 2.D0 
      s8 = 1.D0/(1.D0+e1*cos(f1)) 
      s10 = 2.D0*e1*q1**2*(1.D0+e1)**2.D0/(1.D0+e1*cos(f1))**3.D0*sin(f1&
     &)**2.D0+2.D0*e1**2*q1**2*(1.D0+e1)/(1.D0+e1*cos(f1))**3.D0*sin(f1)&
     &**2.D0-3.D0*e1**2*q1**2*(1.D0+e1)**2.D0/(1.D0+e1*cos(f1))**4.D0*si&
     &n(f1)**2.D0*cos(f1)                                               
      s11 = s10+(q1*(1.D0+e1)/(1.D0+e1*cos(f1))**2.D0*sin(f1)**2.D0+e1*q&
     &1/(1.D0+e1*cos(f1))**2.D0*sin(f1)**2.D0-2.D0*e1*q1*(1.D0+e1)/(1.D0&
     &+e1*cos(f1))**3.D0*sin(f1)**2.D0*cos(f1)+q1/(1.D0+e1*cos(f1))*cos(&
     &f1)-q1*(1.D0+e1)/(1.D0+e1*cos(f1))**2.D0*cos(f1)**2.D0)*(e1*q1*(1.&
     &D0+e1)/(1.D0+e1*cos(f1))+KK*q2*(1.D0+e2)/(1.D0+e2*cos(f2))*cos(f2)&
     &+MM*q2*(1.D0+e2)/(1.D0+e2*cos(f2))*sin(f2))                       
      s9 = s11+(e1*q1*(1.D0+e1)/(1.D0+e1*cos(f1))**2.D0*sin(f1)**2.D0+q1&
     &*(1.D0+e1)/(1.D0+e1*cos(f1))*cos(f1))*(q1*(1.D0+e1)/(1.D0+e1*cos(f&
     &1))+e1*q1/(1.D0+e1*cos(f1))-e1*q1*(1.D0+e1)/(1.D0+e1*cos(f1))**2.D&
     &0*cos(f1))-(2.D0*e1*q1*(1.D0+e1)/(1.D0+e1*cos(f1))**2.D0*sin(f1)+e&
     &1**2*q1/(1.D0+e1*cos(f1))**2.D0*sin(f1)-2.D0*e1**2*q1*(1.D0+e1)/(1&
     &.D0+e1*cos(f1))**3.D0*sin(f1)*cos(f1)+2.D0*q1*(1.D0+e1)/(1.D0+e1*c&
     &os(f1))**2.D0*cos(f1)*sin(f1)+e1*q1/(1.D0+e1*cos(f1))**2.D0*cos(f1&
     &)*sin(f1)-2.D0*e1*q1*(1.D0+e1)/(1.D0+e1*cos(f1))**3.D0*cos(f1)**2.&
     &D0*sin(f1)-q1/(1.D0+e1*cos(f1))*sin(f1))*(LL*q2*(1.D0+e2)/(1.D0+e2&
     &*cos(f2))*cos(f2)+NN*q2*(1.D0+e2)/(1.D0+e2*cos(f2))*sin(f2))      
      s7 = s8*s9 
      s5 = s6*s7 
      s3 = s4+s5 
      dH11_de1 = s2+s3 
                                                                        
      dH11_di1 = 0.D0 
                                                                        
      dH11_dOm1 = 0.D0 
                                                                        
      dH11_domeg1 = 0.D0 
                                                                        
      dH11_dq2 = 2.D0/(1.D0+e1*cos(f1))*((e1*q1*(1.D0+e1)/(1.D0+e1*cos(f&
     &1))**2.D0*sin(f1)**2.D0+q1*(1.D0+e1)/(1.D0+e1*cos(f1))*cos(f1))*(K&
     &K*(1.D0+e2)/(1.D0+e2*cos(f2))*cos(f2)+MM*(1.D0+e2)/(1.D0+e2*cos(f2&
     &))*sin(f2))-(e1**2*q1*(1.D0+e1)/(1.D0+e1*cos(f1))**2.D0*sin(f1)+e1&
     &*q1*(1.D0+e1)/(1.D0+e1*cos(f1))**2.D0*cos(f1)*sin(f1)-q1*(1.D0+e1)&
     &/(1.D0+e1*cos(f1))*sin(f1))*(LL*(1.D0+e2)/(1.D0+e2*cos(f2))*cos(f2&
     &)+NN*(1.D0+e2)/(1.D0+e2*cos(f2))*sin(f2)))                        
                                                                        
      dH11_de2 = 2.D0/(1.D0+e1*cos(f1))*((e1*q1*(1.D0+e1)/(1.D0+e1*cos(f&
     &1))**2.D0*sin(f1)**2.D0+q1*(1.D0+e1)/(1.D0+e1*cos(f1))*cos(f1))*(K&
     &K*q2/(1.D0+e2*cos(f2))*cos(f2)-KK*q2*(1.D0+e2)/(1.D0+e2*cos(f2))**&
     &2.D0*cos(f2)**2.D0+MM*q2/(1.D0+e2*cos(f2))*sin(f2)-MM*q2*(1.D0+e2)&
     &/(1.D0+e2*cos(f2))**2.D0*sin(f2)*cos(f2))-(e1**2*q1*(1.D0+e1)/(1.D&
     &0+e1*cos(f1))**2.D0*sin(f1)+e1*q1*(1.D0+e1)/(1.D0+e1*cos(f1))**2.D&
     &0*cos(f1)*sin(f1)-q1*(1.D0+e1)/(1.D0+e1*cos(f1))*sin(f1))*(LL*q2/(&
     &1.D0+e2*cos(f2))*cos(f2)-LL*q2*(1.D0+e2)/(1.D0+e2*cos(f2))**2.D0*c&
     &os(f2)**2.D0+NN*q2/(1.D0+e2*cos(f2))*sin(f2)-NN*q2*(1.D0+e2)/(1.D0&
     &+e2*cos(f2))**2.D0*sin(f2)*cos(f2)))                              
                                                                        
      dH11_di2 = 0.D0 
                                                                        
      dH11_dOm2 = 0.D0 
                                                                        
      dH11_domeg2 = 0.D0 
                                                                        
      dH22_dq1 = 2.D0/(1.D0+e2*cos(f2))*((e2*q2*(1.D0+e2)/(1.D0+e2*cos(f&
     &2))**2.D0*sin(f2)**2.D0+q2*(1.D0+e2)/(1.D0+e2*cos(f2))*cos(f2))*(K&
     &K*(1.D0+e1)/(1.D0+e1*cos(f1))*cos(f1)+LL*(1.D0+e1)/(1.D0+e1*cos(f1&
     &))*sin(f1))-(e2**2*q2*(1.D0+e2)/(1.D0+e2*cos(f2))**2.D0*sin(f2)+e2&
     &*q2*(1.D0+e2)/(1.D0+e2*cos(f2))**2.D0*cos(f2)*sin(f2)-q2*(1.D0+e2)&
     &/(1.D0+e2*cos(f2))*sin(f2))*(MM*(1.D0+e1)/(1.D0+e1*cos(f1))*cos(f1&
     &)+NN*(1.D0+e1)/(1.D0+e1*cos(f1))*sin(f1)))                        
                                                                        
      dH22_de1 = 2.D0/(1.D0+e2*cos(f2))*((e2*q2*(1.D0+e2)/(1.D0+e2*cos(f&
     &2))**2.D0*sin(f2)**2.D0+q2*(1.D0+e2)/(1.D0+e2*cos(f2))*cos(f2))*(K&
     &K*q1/(1.D0+e1*cos(f1))*cos(f1)-KK*q1*(1.D0+e1)/(1.D0+e1*cos(f1))**&
     &2.D0*cos(f1)**2.D0+LL*q1/(1.D0+e1*cos(f1))*sin(f1)-LL*q1*(1.D0+e1)&
     &/(1.D0+e1*cos(f1))**2.D0*sin(f1)*cos(f1))-(e2**2*q2*(1.D0+e2)/(1.D&
     &0+e2*cos(f2))**2.D0*sin(f2)+e2*q2*(1.D0+e2)/(1.D0+e2*cos(f2))**2.D&
     &0*cos(f2)*sin(f2)-q2*(1.D0+e2)/(1.D0+e2*cos(f2))*sin(f2))*(MM*q1/(&
     &1.D0+e1*cos(f1))*cos(f1)-MM*q1*(1.D0+e1)/(1.D0+e1*cos(f1))**2.D0*c&
     &os(f1)**2.D0+NN*q1/(1.D0+e1*cos(f1))*sin(f1)-NN*q1*(1.D0+e1)/(1.D0&
     &+e1*cos(f1))**2.D0*sin(f1)*cos(f1)))                              
                                                                        
      dH22_di1 = 0.D0 
                                                                        
      dH22_dOm1 = 0.D0 
                                                                        
      dH22_domeg1 = 0.D0 
                                                                        
      s1 = 2.D0 
      s2 = 1.D0/(1.D0+e2*cos(f2))*(2.D0*e2**2*q2*(1.D0+e2)**2.D0/(1.D0+e&
     &2*cos(f2))**3.D0*sin(f2)**2.D0+(e2*(1.D0+e2)/(1.D0+e2*cos(f2))**2.&
     &D0*sin(f2)**2.D0+(1.D0+e2)/(1.D0+e2*cos(f2))*cos(f2))*(e2*q2*(1.D0&
     &+e2)/(1.D0+e2*cos(f2))+KK*q1*(1.D0+e1)/(1.D0+e1*cos(f1))*cos(f1)+L&
     &L*q1*(1.D0+e1)/(1.D0+e1*cos(f1))*sin(f1))+(e2*q2*(1.D0+e2)/(1.D0+e&
     &2*cos(f2))**2.D0*sin(f2)**2.D0+q2*(1.D0+e2)/(1.D0+e2*cos(f2))*cos(&
     &f2))*e2*(1.D0+e2)/(1.D0+e2*cos(f2))-(e2**2*(1.D0+e2)/(1.D0+e2*cos(&
     &f2))**2.D0*sin(f2)+e2*(1.D0+e2)/(1.D0+e2*cos(f2))**2.D0*cos(f2)*si&
     &n(f2)-(1.D0+e2)/(1.D0+e2*cos(f2))*sin(f2))*(MM*q1*(1.D0+e1)/(1.D0+&
     &e1*cos(f1))*cos(f1)+NN*q1*(1.D0+e1)/(1.D0+e1*cos(f1))*sin(f1)))   
      dH22_dq2 = s1*s2 
                                                                        
      s1 = 1.D0/(1.D0+e2*cos(f2))*sin(f2)*gradd2f2-e2/(1.D0+e2*cos(f2))*&
     &*2.D0*sin(f2)*gradd2f2*cos(f2)                                    
      s2 = s1 
      s4 = -2.D0/(1.D0+e2*cos(f2))**2.D0*(e2**2*q2**2*(1.D0+e2)**2.D0/(1&
     &.D0+e2*cos(f2))**3.D0*sin(f2)**2.D0+(e2*q2*(1.D0+e2)/(1.D0+e2*cos(&
     &f2))**2.D0*sin(f2)**2.D0+q2*(1.D0+e2)/(1.D0+e2*cos(f2))*cos(f2))*(&
     &e2*q2*(1.D0+e2)/(1.D0+e2*cos(f2))+KK*q1*(1.D0+e1)/(1.D0+e1*cos(f1)&
     &)*cos(f1)+LL*q1*(1.D0+e1)/(1.D0+e1*cos(f1))*sin(f1))-(e2**2*q2*(1.&
     &D0+e2)/(1.D0+e2*cos(f2))**2.D0*sin(f2)+e2*q2*(1.D0+e2)/(1.D0+e2*co&
     &s(f2))**2.D0*cos(f2)*sin(f2)-q2*(1.D0+e2)/(1.D0+e2*cos(f2))*sin(f2&
     &))*(MM*q1*(1.D0+e1)/(1.D0+e1*cos(f1))*cos(f1)+NN*q1*(1.D0+e1)/(1.D&
     &0+e1*cos(f1))*sin(f1)))*cos(f2)                                   
      s6 = 2.D0 
      s8 = 1.D0/(1.D0+e2*cos(f2)) 
      s10 = 2.D0*e2*q2**2*(1.D0+e2)**2.D0/(1.D0+e2*cos(f2))**3.D0*sin(f2&
     &)**2.D0+2.D0*e2**2*q2**2*(1.D0+e2)/(1.D0+e2*cos(f2))**3.D0*sin(f2)&
     &**2.D0-3.D0*e2**2*q2**2*(1.D0+e2)**2.D0/(1.D0+e2*cos(f2))**4.D0*si&
     &n(f2)**2.D0*cos(f2)                                               
      s11 = s10+(q2*(1.D0+e2)/(1.D0+e2*cos(f2))**2.D0*sin(f2)**2.D0+e2*q&
     &2/(1.D0+e2*cos(f2))**2.D0*sin(f2)**2.D0-2.D0*e2*q2*(1.D0+e2)/(1.D0&
     &+e2*cos(f2))**3.D0*sin(f2)**2.D0*cos(f2)+q2/(1.D0+e2*cos(f2))*cos(&
     &f2)-q2*(1.D0+e2)/(1.D0+e2*cos(f2))**2.D0*cos(f2)**2.D0)*(e2*q2*(1.&
     &D0+e2)/(1.D0+e2*cos(f2))+KK*q1*(1.D0+e1)/(1.D0+e1*cos(f1))*cos(f1)&
     &+LL*q1*(1.D0+e1)/(1.D0+e1*cos(f1))*sin(f1))                       
      s9 = s11+(e2*q2*(1.D0+e2)/(1.D0+e2*cos(f2))**2.D0*sin(f2)**2.D0+q2&
     &*(1.D0+e2)/(1.D0+e2*cos(f2))*cos(f2))*(q2*(1.D0+e2)/(1.D0+e2*cos(f&
     &2))+e2*q2/(1.D0+e2*cos(f2))-e2*q2*(1.D0+e2)/(1.D0+e2*cos(f2))**2.D&
     &0*cos(f2))-(2.D0*e2*q2*(1.D0+e2)/(1.D0+e2*cos(f2))**2.D0*sin(f2)+e&
     &2**2*q2/(1.D0+e2*cos(f2))**2.D0*sin(f2)-2.D0*e2**2*q2*(1.D0+e2)/(1&
     &.D0+e2*cos(f2))**3.D0*sin(f2)*cos(f2)+2.D0*q2*(1.D0+e2)/(1.D0+e2*c&
     &os(f2))**2.D0*cos(f2)*sin(f2)+e2*q2/(1.D0+e2*cos(f2))**2.D0*cos(f2&
     &)*sin(f2)-2.D0*e2*q2*(1.D0+e2)/(1.D0+e2*cos(f2))**3.D0*cos(f2)**2.&
     &D0*sin(f2)-q2/(1.D0+e2*cos(f2))*sin(f2))*(MM*q1*(1.D0+e1)/(1.D0+e1&
     &*cos(f1))*cos(f1)+NN*q1*(1.D0+e1)/(1.D0+e1*cos(f1))*sin(f1))      
      s7 = s8*s9 
      s5 = s6*s7 
      s3 = s4+s5 
      dH22_de2 = s2+s3 
                                                                        
      dH22_di2 = 0.D0 
                                                                        
      dH22_dOm2 = 0.D0 
                                                                        
      dH22_domeg2 = 0.D0 
                                                                        
      dH12_dq1 = 2.D0/(1.D0+e1*cos(f1))*((1.D0+e1)/(1.D0+e1*cos(f1))*sin&
     &(f1)*(KK*(e2*q2*(1.D0+e2)/(1.D0+e2*cos(f2))**2.D0*cos(f2)*sin(f2)-&
     &q2*(1.D0+e2)/(1.D0+e2*cos(f2))*sin(f2))+MM*(e2*q2*(1.D0+e2)/(1.D0+&
     &e2*cos(f2))**2.D0*sin(f2)**2.D0+q2*(1.D0+e2)/(1.D0+e2*cos(f2))*cos&
     &(f2)))-(e1*(1.D0+e1)/(1.D0+e1*cos(f1))+(1.D0+e1)/(1.D0+e1*cos(f1))&
     &*cos(f1))*(LL*(e2*q2*(1.D0+e2)/(1.D0+e2*cos(f2))**2.D0*cos(f2)*sin&
     &(f2)-q2*(1.D0+e2)/(1.D0+e2*cos(f2))*sin(f2))+NN*(e2*q2*(1.D0+e2)/(&
     &1.D0+e2*cos(f2))**2.D0*sin(f2)**2.D0+q2*(1.D0+e2)/(1.D0+e2*cos(f2)&
     &)*cos(f2))))                                                      
                                                                        
      s1 = -2.D0/(1.D0+e1*cos(f1))**2.D0*(q1*(1.D0+e1)/(1.D0+e1*cos(f1))&
     &*sin(f1)*(KK*(e2*q2*(1.D0+e2)/(1.D0+e2*cos(f2))**2.D0*cos(f2)*sin(&
     &f2)-q2*(1.D0+e2)/(1.D0+e2*cos(f2))*sin(f2))+MM*(e2*q2*(1.D0+e2)/(1&
     &.D0+e2*cos(f2))**2.D0*sin(f2)**2.D0+q2*(1.D0+e2)/(1.D0+e2*cos(f2))&
     &*cos(f2)))-(e1*q1*(1.D0+e1)/(1.D0+e1*cos(f1))+q1*(1.D0+e1)/(1.D0+e&
     &1*cos(f1))*cos(f1))*(LL*(e2*q2*(1.D0+e2)/(1.D0+e2*cos(f2))**2.D0*c&
     &os(f2)*sin(f2)-q2*(1.D0+e2)/(1.D0+e2*cos(f2))*sin(f2))+NN*(e2*q2*(&
     &1.D0+e2)/(1.D0+e2*cos(f2))**2.D0*sin(f2)**2.D0+q2*(1.D0+e2)/(1.D0+&
     &e2*cos(f2))*cos(f2))))*cos(f1)                                    
      s3 = 2.D0 
      s5 = 1.D0/(1.D0+e1*cos(f1)) 
      s7 = q1/(1.D0+e1*cos(f1))*sin(f1)*(KK*(e2*q2*(1.D0+e2)/(1.D0+e2*co&
     &s(f2))**2.D0*cos(f2)*sin(f2)-q2*(1.D0+e2)/(1.D0+e2*cos(f2))*sin(f2&
     &))+MM*(e2*q2*(1.D0+e2)/(1.D0+e2*cos(f2))**2.D0*sin(f2)**2.D0+q2*(1&
     &.D0+e2)/(1.D0+e2*cos(f2))*cos(f2)))                               
      s9 = -q1*(1.D0+e1)/(1.D0+e1*cos(f1))**2.D0*sin(f1)*(KK*(e2*q2*(1.D&
     &0+e2)/(1.D0+e2*cos(f2))**2.D0*cos(f2)*sin(f2)-q2*(1.D0+e2)/(1.D0+e&
     &2*cos(f2))*sin(f2))+MM*(e2*q2*(1.D0+e2)/(1.D0+e2*cos(f2))**2.D0*si&
     &n(f2)**2.D0+q2*(1.D0+e2)/(1.D0+e2*cos(f2))*cos(f2)))*cos(f1)      
      s10 = -(q1*(1.D0+e1)/(1.D0+e1*cos(f1))+e1*q1/(1.D0+e1*cos(f1))-e1*&
     &q1*(1.D0+e1)/(1.D0+e1*cos(f1))**2.D0*cos(f1)+q1/(1.D0+e1*cos(f1))*&
     &cos(f1)-q1*(1.D0+e1)/(1.D0+e1*cos(f1))**2.D0*cos(f1)**2.D0)*(LL*(e&
     &2*q2*(1.D0+e2)/(1.D0+e2*cos(f2))**2.D0*cos(f2)*sin(f2)-q2*(1.D0+e2&
     &)/(1.D0+e2*cos(f2))*sin(f2))+NN*(e2*q2*(1.D0+e2)/(1.D0+e2*cos(f2))&
     &**2.D0*sin(f2)**2.D0+q2*(1.D0+e2)/(1.D0+e2*cos(f2))*cos(f2)))     
      s8 = s9+s10 
      s6 = s7+s8 
      s4 = s5*s6 
      s2 = s3*s4 
      dH12_de1 = s1+s2 
                                                                        
      dH12_di1 = 0.D0 
                                                                        
      dH12_dOm1 = 0.D0 
                                                                        
      dH12_domeg1 = 0.D0 
                                                                        
      dH12_dq2 = 2.D0/(1.D0+e1*cos(f1))*(q1*(1.D0+e1)/(1.D0+e1*cos(f1))*&
     &sin(f1)*(KK*(e2*(1.D0+e2)/(1.D0+e2*cos(f2))**2.D0*cos(f2)*sin(f2)-&
     &(1.D0+e2)/(1.D0+e2*cos(f2))*sin(f2))+MM*(e2*(1.D0+e2)/(1.D0+e2*cos&
     &(f2))**2.D0*sin(f2)**2.D0+(1.D0+e2)/(1.D0+e2*cos(f2))*cos(f2)))-(e&
     &1*q1*(1.D0+e1)/(1.D0+e1*cos(f1))+q1*(1.D0+e1)/(1.D0+e1*cos(f1))*co&
     &s(f1))*(LL*(e2*(1.D0+e2)/(1.D0+e2*cos(f2))**2.D0*cos(f2)*sin(f2)-(&
     &1.D0+e2)/(1.D0+e2*cos(f2))*sin(f2))+NN*(e2*(1.D0+e2)/(1.D0+e2*cos(&
     &f2))**2.D0*sin(f2)**2.D0+(1.D0+e2)/(1.D0+e2*cos(f2))*cos(f2))))   
                                                                        
      s1 = 2.D0 
      s3 = 1.D0/(1.D0+e1*cos(f1)) 
      s5 = q1*(1.D0+e1)/(1.D0+e1*cos(f1))*sin(f1)*(KK*(2.D0*q2*(1.D0+e2)&
     &/(1.D0+e2*cos(f2))**2.D0*cos(f2)*sin(f2)+e2*q2/(1.D0+e2*cos(f2))**&
     &2.D0*cos(f2)*sin(f2)-2.D0*e2*q2*(1.D0+e2)/(1.D0+e2*cos(f2))**3.D0*&
     &cos(f2)**2.D0*sin(f2)-q2/(1.D0+e2*cos(f2))*sin(f2))+MM*(q2*(1.D0+e&
     &2)/(1.D0+e2*cos(f2))**2.D0*sin(f2)**2.D0+e2*q2/(1.D0+e2*cos(f2))**&
     &2.D0*sin(f2)**2.D0-2.D0*e2*q2*(1.D0+e2)/(1.D0+e2*cos(f2))**3.D0*si&
     &n(f2)**2.D0*cos(f2)+q2/(1.D0+e2*cos(f2))*cos(f2)-q2*(1.D0+e2)/(1.D&
     &0+e2*cos(f2))**2.D0*cos(f2)**2.D0))                               
      s6 = -(e1*q1*(1.D0+e1)/(1.D0+e1*cos(f1))+q1*(1.D0+e1)/(1.D0+e1*cos&
     &(f1))*cos(f1))*(LL*(2.D0*q2*(1.D0+e2)/(1.D0+e2*cos(f2))**2.D0*cos(&
     &f2)*sin(f2)+e2*q2/(1.D0+e2*cos(f2))**2.D0*cos(f2)*sin(f2)-2.D0*e2*&
     &q2*(1.D0+e2)/(1.D0+e2*cos(f2))**3.D0*cos(f2)**2.D0*sin(f2)-q2/(1.D&
     &0+e2*cos(f2))*sin(f2))+NN*(q2*(1.D0+e2)/(1.D0+e2*cos(f2))**2.D0*si&
     &n(f2)**2.D0+e2*q2/(1.D0+e2*cos(f2))**2.D0*sin(f2)**2.D0-2.D0*e2*q2&
     &*(1.D0+e2)/(1.D0+e2*cos(f2))**3.D0*sin(f2)**2.D0*cos(f2)+q2/(1.D0+&
     &e2*cos(f2))*cos(f2)-q2*(1.D0+e2)/(1.D0+e2*cos(f2))**2.D0*cos(f2)**&
     &2.D0))                                                            
      s4 = s5+s6 
      s2 = s3*s4 
      dH12_de2 = s1*s2 
                                                                        
      dH12_di2 = 0.D0 
                                                                        
      dH12_dOm2 = 0.D0 
                                                                        
      dH12_domeg2 = 0.D0 
                                                                        
!=====================================================================
! "Total" partial derivatives of elements of matrix H w.r.t. the elements
! computed by the chain rule
! NOTE: we neglect the part of the chain rule containing the derivatives 
! of the gradient of d^2
!=====================================================================
   derH11_dq1 = dH11_dq1 + dH11_df1*df1_dq1+dH11_df2*df2_dq1

   derH11_de1 = dH11_de1 + dH11_df1*df1_de1+dH11_df2*df2_de1

   derH11_di1 = dH11_dKK*dKKdi1+dH11_dLL*dLLdi1+ &
        & dH11_dMM*dMMdi1+dH11_dNN*dNNdi1+dH11_df1*df1_di1+dH11_df2*df2_di1

   derH11_dOm1 = dH11_dKK*dKKdOm1+dH11_dLL*dLLdOm1+ &
        & dH11_dMM*dMMdOm1+dH11_dNN*dNNdOm1+dH11_df1*df1_dOm1+&
        & dH11_df2*df2_dOm1

   derH11_domeg1 = dH11_dKK*dKKdomeg1+dH11_dLL*dLLdomeg1+ &
        & dH11_dMM*dMMdomeg1+dH11_dNN*dNNdomeg1+dH11_df1*df1_domeg1+ &
        & dH11_df2*df2_domeg1 

   derH11_dq2 = dH11_dq2 + dH11_df1*df1_dq2 + dH11_df2*df2_dq2 

   derH11_de2 = dH11_de2 + dH11_df1*df1_de2+dH11_df2*df2_de2

   derH11_di2 = dH11_dKK*dKKdi2+dH11_dLL*dLLdi2+ &
        & dH11_dMM*dMMdi2+dH11_dNN*dNNdi2+dH11_df1*df1_di2+dH11_df2*df2_di2

   derH11_dOm2 = dH11_dKK*dKKdOm2+dH11_dLL*dLLdOm2+ &
        & dH11_dMM*dMMdOm2+dH11_dNN*dNNdOm2+dH11_df1*df1_dOm2+&
        & dH11_df2*df2_dOm2 

   derH11_domeg2 = dH11_dKK*dKKdomeg2+dH11_dLL*dLLdomeg2+ &
        & dH11_dMM*dMMdomeg2+dH11_dNN*dNNdomeg2+ dH11_df1*df1_domeg2+ &
        & dH11_df2*df2_domeg2 

   derH22_dq1 = dH22_dq1 + dH22_df1*df1_dq1+dH22_df2*df2_dq1 

   derH22_de1 = dH22_de1 + dH22_df1*df1_de1+dH22_df2*df2_de1

   derH22_di1 = dH22_dKK*dKKdi1+dH22_dLL*dLLdi1+ &
        & dH22_dMM*dMMdi1+dH22_dNN*dNNdi1+dH22_df1*df1_di1+dH22_df2*df2_di1

   derH22_dOm1 = dH22_dKK*dKKdOm1+dH22_dLL*dLLdOm1+ &
        & dH22_dMM*dMMdOm1+dH22_dNN*dNNdOm1+dH22_df1*df1_dOm1+&
        & dH22_df2*df2_dOm1 

   derH22_domeg1 = dH22_dKK*dKKdomeg1+dH22_dLL*dLLdomeg1+ &
        & dH22_dMM*dMMdomeg1+dH22_dNN*dNNdomeg1+dH22_df1*df1_domeg1+ &
        & dH22_df2*df2_domeg1

   derH22_dq2 = dH22_dq2 + dH22_df1*df1_dq2 + dH22_df2*df2_dq2 

   derH22_de2 = dH22_de2 + dH22_df1*df1_de2+dH22_df2*df2_de2

   derH22_di2 = dH22_dKK*dKKdi2+dH22_dLL*dLLdi2+ &
        & dH22_dMM*dMMdi2+dH22_dNN*dNNdi2+dH22_df1*df1_di2+dH22_df2*df2_di2

   derH22_dOm2 = dH22_dKK*dKKdOm2+dH22_dLL*dLLdOm2+ &
        & dH22_dMM*dMMdOm2+dH22_dNN*dNNdOm2+dH22_df1*df1_dOm2+&
        & dH22_df2*df2_dOm2

   derH22_domeg2 = dH22_dKK*dKKdomeg2+dH22_dLL*dLLdomeg2+ &
        & dH22_dMM*dMMdomeg2+dH22_dNN*dNNdomeg2+dH22_df1*df1_domeg2+ &
        & dH22_df2*df2_domeg2

! Note that H12 does NOT depend on gradd2f1 or gradd2f2
   derH12_dq1 = dH12_dq1 + dH12_df1*df1_dq1+dH12_df2*df2_dq1

   derH12_de1 = dH12_de1 + dH12_df1*df1_de1+dH12_df2*df2_de1

   derH12_di1 = dH12_dKK*dKKdi1+dH12_dLL*dLLdi1+ &
        & dH12_dMM*dMMdi1+dH12_dNN*dNNdi1+dH12_df1*df1_di1+dH12_df2*df2_di1

   derH12_dOm1 = dH12_dKK*dKKdOm1+dH12_dLL*dLLdOm1+ &
        & dH12_dMM*dMMdOm1+dH12_dNN*dNNdOm1+dH12_df1*df1_dOm1+ &
        & dH12_df2*df2_dOm1

   derH12_domeg1 = dH12_dKK*dKKdomeg1+dH12_dLL*dLLdomeg1+ &
        & dH12_dMM*dMMdomeg1+dH12_dNN*dNNdomeg1+dH12_df1*df1_domeg1+ &
        & dH12_df2*df2_domeg1

   derH12_dq2 = dH12_dq2 + dH12_df1*df1_dq2+dH12_df2*df2_dq2

   derH12_de2 = dH12_de2 + dH12_df1*df1_de2+dH12_df2*df2_de2

   derH12_di2 = dH12_dKK*dKKdi2+dH12_dLL*dLLdi2+ &
        & dH12_dMM*dMMdi2+dH12_dNN*dNNdi2+dH12_df1*df1_di2+dH12_df2*df2_di2

   derH12_dOm2 = dH12_dKK*dKKdOm2+dH12_dLL*dLLdOm2+ &
        & dH12_dMM*dMMdOm2+dH12_dNN*dNNdOm2+dH12_df1*df1_dOm2+ &
        & dH12_df2*df2_dOm2

   derH12_domeg2 = dH12_dKK*dKKdomeg2+dH12_dLL*dLLdomeg2+ &
        & dH12_dMM*dMMdomeg2+dH12_dNN*dNNdomeg2+dH12_df1*df1_domeg2+ &
        & dH12_df2*df2_domeg2

!=====================================================================
! Derivatives of the determinant
!=====================================================================
  ddetH_dq1=derH11_dq1*Hess(2,2)+Hess(1,1)*derH22_dq1- &
       & 2.d0*Hess(1,2)*derH12_dq1
  ddetH_de1=derH11_de1*Hess(2,2)+Hess(1,1)*derH22_de1- &
       & 2.d0*Hess(1,2)*derH12_de1
  ddetH_di1=derH11_di1*Hess(2,2)+Hess(1,1)*derH22_di1- &
       & 2.d0*Hess(1,2)*derH12_di1
  ddetH_dOm1=derH11_dOm1*Hess(2,2)+Hess(1,1)*derH22_dOm1-2.d0*Hess(1,2)* &
       & derH12_dOm1
  ddetH_domeg1=derH11_domeg1*Hess(2,2)+Hess(1,1)*derH22_domeg1- &
       & 2.d0*Hess(1,2)*derH12_domeg1
  ddetH_dq2=derH11_dq2*Hess(2,2)+Hess(1,1)*derH22_dq2- &
       & 2.d0*Hess(1,2)*derH12_dq2
  ddetH_de2=derH11_de2*Hess(2,2)+Hess(1,1)*derH22_de2- &
       & 2.d0*Hess(1,2)*derH12_de2
  ddetH_di2=derH11_di2*Hess(2,2)+Hess(1,1)*derH22_di2- &
       & 2.d0*Hess(1,2)*derH12_di2
  ddetH_dOm2=derH11_dOm2*Hess(2,2)+Hess(1,1)*derH22_dOm2-2.d0*Hess(1,2)* &
       & derH12_dOm2
  ddetH_domeg2=derH11_domeg2*Hess(2,2)+Hess(1,1)*derH22_domeg2- &
       & 2.d0*Hess(1,2)*derH12_domeg2

  ddetH(1,1)=ddetH_dq1
  ddetH(1,2)=ddetH_de1
  ddetH(1,3)=ddetH_di1
  ddetH(1,4)=ddetH_dOm1
  ddetH(1,5)=ddetH_domeg1
  ddetH(1,6)=ddetH_dq2
  ddetH(1,7)=ddetH_de2
  ddetH(1,8)=ddetH_di2
  ddetH(1,9)=ddetH_dOm2
  ddetH(1,10)=ddetH_domeg2
  tddetH=TRANSPOSE(ddetH)
  detHrmsvec=MATMUL(ddetH,MATMUL(covcom,tddetH)) 
  detHrms=SQRT(detHrmsvec(1,1))

  IF(chk_der) THEN
!     ddetHdcom2(1:1,1:5)= ddetH(1:1,6:10)
     write(ierrou,*)'derivatives of detH:'
     write(ierrou,100) ddetH(1:1,6:10)
     numerr=numerr+1
  ENDIF

!=================================================================
! RMS of tau1 || tau2
!=================================================================
! pericenter pos. and vel. of 1st orbit  
  r01=q1
  mu1=centerbody_mass(com1%center)
  v01=sqrt(mu1*(1.d0+e1)/q1)
  cosn1=cos(com1%coord(4))
  sinn1=sin(com1%coord(4))
  coso1=cos(com1%coord(5))
  sino1=sin(com1%coord(5))
  cosi1=cos(com1%coord(3))
  sini1=sin(com1%coord(3))
  x01=(/coso1*cosn1-sino1*sinn1*cosi1, coso1*sinn1+sino1*cosn1*cosi1,&
       & sino1*sini1/)*r01
  y01=(/-sino1*cosn1-coso1*sinn1*cosi1, -sino1*sinn1+coso1*cosn1*cosi1,&
       & coso1*sini1/)*v01
  alpha1 = vsize(y01)**2-2.d0*mu1/r01
  r1 = q1*(1.d0+e1)/(1.d0+e1*cos(f1))
  IF(r1.lt.0.d0.and.verb_moid.ge.20)THEN
     WRITE(ierrou,*)'negative r1',r1
     numerr=numerr+1
  ENDIF
  
! pericenter pos. and vel. of 2nd orbit  
  r02=q2
  mu2=centerbody_mass(com2%center)
!  mu2=gms
  v02=sqrt(mu2*(1.d0+e2)/q2)
  cosn2=cos(com2%coord(4))
  sinn2=sin(com2%coord(4))
  coso2=cos(com2%coord(5))
  sino2=sin(com2%coord(5))
  cosi2=cos(com2%coord(3))
  sini2=sin(com2%coord(3))
  x02=(/coso2*cosn2-sino2*sinn2*cosi2, coso2*sinn2+sino2*cosn2*cosi2,&
       & sino2*sini2/)*r02
  y02=(/-sino2*cosn2-coso2*sinn2*cosi2, -sino2*sinn2+coso2*cosn2*cosi2,&
       & coso2*sini2/)*v02
  alpha2 = vsize(y02)**2-2.d0*mu2/r02
  r2 = q2*(1.d0+e2)/(1.d0+e2*cos(f2))
  IF(r2.lt.0.d0.and.verb_moid.ge.20)THEN
     WRITE(ierrou,*)'negative r2',r2
     numerr=numerr+1
  ENDIF

! starting guess for the Earth
!  IF(e1.lt.1.d0) THEN
!     u1=2*ATAN(TAN(f1/2.d0)*SQRT((1.d0-e1)/(1.d0+e1)))
!     write(*,*)'f1=',f1*degrad,'u1=',u1*degrad
!     psi1_ini = u1/sqrt(-alpha1)
!     write(*,*)'psi1_ini=',psi1_ini,'u1=',u1,'alpha1=',alpha1
!  ELSEIF(e1.eq.1.d0) THEN
!     u1 = tan(f1/2.d0)
!     psi1_ini = u1*sqrt(2.d0*q1/mu1)
!  ELSEIF(e1.gt.1.d0) THEN
!     u1x=TAN(f1/2.d0)*SQRT((e1-1.d0)/(1.d0+e1))
!     u1=log((1.d0+u1x)/(1.d0-u1x))
!     psi1_ini = u1/sqrt(alpha1)
!  ENDIF
! starting guess for the asteroid
!  IF(e2.lt.1.d0) THEN
!     u2=2*ATAN(TAN(f2/2.d0)*SQRT((1.d0-e2)/(1.d0+e2)))
!     psi2_ini = u2/sqrt(-alpha2)
!  ELSEIF(e2.eq.1.d0) THEN
!     u2 = tan(f2/2.d0)
!     psi2_ini = u2*sqrt(2.d0*q2/mu2)
!  ELSEIF(e2.gt.1.d0) THEN
!     u2x=TAN(f2/2.d0)*SQRT((e2-1.d0)/(1.d0+e2))
!     u2=log((1.d0+u2x)/(1.d0-u2x))
!     psi2_ini = u2/sqrt(alpha2)
!  ENDIF
! time of passage at perihelion
!  CALL pdtime(q1,r1,mu1,alpha1,psi1_ini,psi1,dt1,1.d-10)
!  ptime1 = -dt1+com1%t
!  CALL pdtime(q2,r2,mu2,alpha2,psi2_ini,psi2,dt2,1.d-10)
!  ptime2 = -dt2+com2%t

  com1min=com1
  com2min=com2
! substitute last element time of perihelion passage
  com1min%coord(6)=f1
  com2min%coord(6)=f2
  CALL coo_cha(com1min,'CAR',car1min,fail_flag,dcardcom1)
  CALL coo_cha(com2min,'CAR',car2min,fail_flag,dcardcom2)       

! tangent vectors at minimum points
  tau1(1) = (cos(Om1)*cos(omeg1)-sin(Om1)*cos(i1)*sin(omeg1))*q1*(1.D&
     &0+e1)/(1.D0+e1*cos(f1))**2.D0*cos(f1)*e1*sin(f1)-(cos(Om1)*cos(ome&
     &g1)-sin(Om1)*cos(i1)*sin(omeg1))*q1*(1.D0+e1)/(1.D0+e1*cos(f1))*si&
     &n(f1)+(-cos(Om1)*sin(omeg1)-sin(Om1)*cos(i1)*cos(omeg1))*q1*(1.D0+&
     &e1)/(1.D0+e1*cos(f1))**2.D0*sin(f1)**2.D0*e1+(-cos(Om1)*sin(omeg1)&
     &-sin(Om1)*cos(i1)*cos(omeg1))*q1*(1.D0+e1)/(1.D0+e1*cos(f1))*cos(f1)
                                                                        
  tau1(2) = (sin(Om1)*cos(omeg1)+cos(Om1)*cos(i1)*sin(omeg1))*q1*(1.D&
     &0+e1)/(1.D0+e1*cos(f1))**2.D0*cos(f1)*e1*sin(f1)-(sin(Om1)*cos(ome&
     &g1)+cos(Om1)*cos(i1)*sin(omeg1))*q1*(1.D0+e1)/(1.D0+e1*cos(f1))*si&
     &n(f1)+(-sin(Om1)*sin(omeg1)+cos(Om1)*cos(i1)*cos(omeg1))*q1*(1.D0+&
     &e1)/(1.D0+e1*cos(f1))**2.D0*sin(f1)**2.D0*e1+(-sin(Om1)*sin(omeg1)&
     &+cos(Om1)*cos(i1)*cos(omeg1))*q1*(1.D0+e1)/(1.D0+e1*cos(f1))*cos(f1)
                                                                        
  tau1(3) = sin(i1)*sin(omeg1)*q1*(1.D0+e1)/(1.D0+e1*cos(f1))**2.D0*c&
     &os(f1)*e1*sin(f1)-sin(i1)*sin(omeg1)*q1*(1.D0+e1)/(1.D0+e1*cos(f1)&
     &)*sin(f1)+sin(i1)*cos(omeg1)*q1*(1.D0+e1)/(1.D0+e1*cos(f1))**2.D0*&
     &sin(f1)**2.D0*e1+sin(i1)*cos(omeg1)*q1*(1.D0+e1)/(1.D0+e1*cos(f1))&
     &*cos(f1)                                                          
                                                                        
  tau2(1) = (cos(Om2)*cos(omeg2)-sin(Om2)*cos(i2)*sin(omeg2))*q2*(1.D&
     &0+e2)/(1.D0+e2*cos(f2))**2.D0*cos(f2)*e2*sin(f2)-(cos(Om2)*cos(ome&
     &g2)-sin(Om2)*cos(i2)*sin(omeg2))*q2*(1.D0+e2)/(1.D0+e2*cos(f2))*si&
     &n(f2)+(-cos(Om2)*sin(omeg2)-sin(Om2)*cos(i2)*cos(omeg2))*q2*(1.D0+&
     &e2)/(1.D0+e2*cos(f2))**2.D0*sin(f2)**2.D0*e2+(-cos(Om2)*sin(omeg2)&
     &-sin(Om2)*cos(i2)*cos(omeg2))*q2*(1.D0+e2)/(1.D0+e2*cos(f2))*cos(f2)
                                                                        
  tau2(2) = (sin(Om2)*cos(omeg2)+cos(Om2)*cos(i2)*sin(omeg2))*q2*(1.D&
     &0+e2)/(1.D0+e2*cos(f2))**2.D0*cos(f2)*e2*sin(f2)-(sin(Om2)*cos(ome&
     &g2)+cos(Om2)*cos(i2)*sin(omeg2))*q2*(1.D0+e2)/(1.D0+e2*cos(f2))*si&
     &n(f2)+(-sin(Om2)*sin(omeg2)+cos(Om2)*cos(i2)*cos(omeg2))*q2*(1.D0+&
     &e2)/(1.D0+e2*cos(f2))**2.D0*sin(f2)**2.D0*e2+(-sin(Om2)*sin(omeg2)&
     &+cos(Om2)*cos(i2)*cos(omeg2))*q2*(1.D0+e2)/(1.D0+e2*cos(f2))*cos(f2)
                                                                        
  tau2(3) = sin(i2)*sin(omeg2)*q2*(1.D0+e2)/(1.D0+e2*cos(f2))**2.D0*c&
     &os(f2)*e2*sin(f2)-sin(i2)*sin(omeg2)*q2*(1.D0+e2)/(1.D0+e2*cos(f2)&
     &)*sin(f2)+sin(i2)*cos(omeg2)*q2*(1.D0+e2)/(1.D0+e2*cos(f2))**2.D0*&
     &sin(f2)**2.D0*e2+sin(i2)*cos(omeg2)*q2*(1.D0+e2)/(1.D0+e2*cos(f2))&
     &*cos(f2)                                                          
                                                                        
! compute joining unit vector tau3hat
  CALL prvec(tau1,tau2,tau3)
  tau3hat=tau3/vsize(tau3)  

! -------------------------------------
  ltau1=vsize(tau1)
  tau1vers(1:3)=tau1(1:3)/ltau1

  ltau2=vsize(tau2)
  tau2vers(1:3)=tau2(1:3)/ltau2
!  write(*,*)'tau2vers',tau2vers(1:3)

  CALL prvec(tau1vers,tau2vers,tau3vec)  ! tau3vec has the same direction
                                         ! as tau3 but it is not the same
  sint1t2=vsize(tau3vec)

! ======================================================================

      dtau1dq1(1) = (cos(Om1)*cos(omeg1)-sin(Om1)*cos(i1)*sin(omeg1))*(1.&
     &D0+e1)/(1.D0+e1*cos(f1))**2.D0*cos(f1)*e1*sin(f1)-(cos(Om1)*cos(om&
     &eg1)-sin(Om1)*cos(i1)*sin(omeg1))*(1.D0+e1)/(1.D0+e1*cos(f1))*sin(&
     &f1)+(-cos(Om1)*sin(omeg1)-sin(Om1)*cos(i1)*cos(omeg1))*(1.D0+e1)/(&
     &1.D0+e1*cos(f1))**2.D0*sin(f1)**2.D0*e1+(-cos(Om1)*sin(omeg1)-sin(&
     &Om1)*cos(i1)*cos(omeg1))*(1.D0+e1)/(1.D0+e1*cos(f1))*cos(f1)      
                                                                        
      s1 = (cos(Om1)*cos(omeg1)-sin(Om1)*cos(i1)*sin(omeg1))*q1/(1.D0+e1&
     &*cos(f1))**2.D0*cos(f1)*e1*sin(f1)-2.D0*(cos(Om1)*cos(omeg1)-sin(O&
     &m1)*cos(i1)*sin(omeg1))*q1*(1.D0+e1)/(1.D0+e1*cos(f1))**3.D0*cos(f&
     &1)**2.D0*e1*sin(f1)+2.D0*(cos(Om1)*cos(omeg1)-sin(Om1)*cos(i1)*sin&
     &(omeg1))*q1*(1.D0+e1)/(1.D0+e1*cos(f1))**2.D0*cos(f1)*sin(f1)-(cos&
     &(Om1)*cos(omeg1)-sin(Om1)*cos(i1)*sin(omeg1))*q1/(1.D0+e1*cos(f1))&
     &*sin(f1)                                                          
      dtau1de1(1) = s1+(-cos(Om1)*sin(omeg1)-sin(Om1)*cos(i1)*cos(omeg1))&
     &*q1/(1.D0+e1*cos(f1))**2.D0*sin(f1)**2.D0*e1-2.D0*(-cos(Om1)*sin(o&
     &meg1)-sin(Om1)*cos(i1)*cos(omeg1))*q1*(1.D0+e1)/(1.D0+e1*cos(f1))*&
     &*3.D0*sin(f1)**2.D0*e1*cos(f1)+(-cos(Om1)*sin(omeg1)-sin(Om1)*cos(&
     &i1)*cos(omeg1))*q1*(1.D0+e1)/(1.D0+e1*cos(f1))**2.D0*sin(f1)**2.D0&
     &+(-cos(Om1)*sin(omeg1)-sin(Om1)*cos(i1)*cos(omeg1))*q1/(1.D0+e1*co&
     &s(f1))*cos(f1)-(-cos(Om1)*sin(omeg1)-sin(Om1)*cos(i1)*cos(omeg1))*&
     &q1*(1.D0+e1)/(1.D0+e1*cos(f1))**2.D0*cos(f1)**2.D0                
                                                                        
      dtau1di1(1) = sin(Om1)*sin(i1)*sin(omeg1)*q1*(1.D0+e1)/(1.D0+e1*cos&
     &(f1))**2.D0*cos(f1)*e1*sin(f1)-sin(Om1)*sin(i1)*sin(omeg1)*q1*(1.D&
     &0+e1)/(1.D0+e1*cos(f1))*sin(f1)+sin(Om1)*sin(i1)*cos(omeg1)*q1*(1.&
     &D0+e1)/(1.D0+e1*cos(f1))**2.D0*sin(f1)**2.D0*e1+sin(Om1)*sin(i1)*c&
     &os(omeg1)*q1*(1.D0+e1)/(1.D0+e1*cos(f1))*cos(f1)                  
                                                                        
      dtau1dOm1(1) = (-sin(Om1)*cos(omeg1)-cos(Om1)*cos(i1)*sin(omeg1))*q&
     &1*(1.D0+e1)/(1.D0+e1*cos(f1))**2.D0*cos(f1)*e1*sin(f1)-(-sin(Om1)*&
     &cos(omeg1)-cos(Om1)*cos(i1)*sin(omeg1))*q1*(1.D0+e1)/(1.D0+e1*cos(&
     &f1))*sin(f1)+(sin(Om1)*sin(omeg1)-cos(Om1)*cos(i1)*cos(omeg1))*q1*&
     &(1.D0+e1)/(1.D0+e1*cos(f1))**2.D0*sin(f1)**2.D0*e1+(sin(Om1)*sin(o&
     &meg1)-cos(Om1)*cos(i1)*cos(omeg1))*q1*(1.D0+e1)/(1.D0+e1*cos(f1))*&
     &cos(f1)                                                           
                                                                        
      dtau1domeg1(1) = (-cos(Om1)*sin(omeg1)-sin(Om1)*cos(i1)*cos(omeg1))&
     &*q1*(1.D0+e1)/(1.D0+e1*cos(f1))**2.D0*cos(f1)*e1*sin(f1)-(-cos(Om1&
     &)*sin(omeg1)-sin(Om1)*cos(i1)*cos(omeg1))*q1*(1.D0+e1)/(1.D0+e1*co&
     &s(f1))*sin(f1)+(-cos(Om1)*cos(omeg1)+sin(Om1)*cos(i1)*sin(omeg1))*&
     &q1*(1.D0+e1)/(1.D0+e1*cos(f1))**2.D0*sin(f1)**2.D0*e1+(-cos(Om1)*c&
     &os(omeg1)+sin(Om1)*cos(i1)*sin(omeg1))*q1*(1.D0+e1)/(1.D0+e1*cos(f&
     &1))*cos(f1)                                                       
                                                                        
      s1 = 2.D0*(cos(Om1)*cos(omeg1)-sin(Om1)*cos(i1)*sin(omeg1))*q1*(1.&
     &D0+e1)/(1.D0+e1*cos(f1))**3.D0*cos(f1)*e1**2*sin(f1)**2.D0-2.D0*(c&
     &os(Om1)*cos(omeg1)-sin(Om1)*cos(i1)*sin(omeg1))*q1*(1.D0+e1)/(1.D0&
     &+e1*cos(f1))**2.D0*sin(f1)**2.D0*e1+(cos(Om1)*cos(omeg1)-sin(Om1)*&
     &cos(i1)*sin(omeg1))*q1*(1.D0+e1)/(1.D0+e1*cos(f1))**2.D0*cos(f1)**&
     &2.D0*e1                                                           
      dtau1df1(1) = s1-(cos(Om1)*cos(omeg1)-sin(Om1)*cos(i1)*sin(omeg1))*&
     &q1*(1.D0+e1)/(1.D0+e1*cos(f1))*cos(f1)+2.D0*(-cos(Om1)*sin(omeg1)-&
     &sin(Om1)*cos(i1)*cos(omeg1))*q1*(1.D0+e1)/(1.D0+e1*cos(f1))**3.D0*&
     &sin(f1)**3.D0*e1**2+3.D0*(-cos(Om1)*sin(omeg1)-sin(Om1)*cos(i1)*co&
     &s(omeg1))*q1*(1.D0+e1)/(1.D0+e1*cos(f1))**2.D0*cos(f1)*e1*sin(f1)-&
     &(-cos(Om1)*sin(omeg1)-sin(Om1)*cos(i1)*cos(omeg1))*q1*(1.D0+e1)/(1&
     &.D0+e1*cos(f1))*sin(f1)                                           
                                                                        
      dtau1dq1(2) = (sin(Om1)*cos(omeg1)+cos(Om1)*cos(i1)*sin(omeg1))*(1.&
     &D0+e1)/(1.D0+e1*cos(f1))**2.D0*cos(f1)*e1*sin(f1)-(sin(Om1)*cos(om&
     &eg1)+cos(Om1)*cos(i1)*sin(omeg1))*(1.D0+e1)/(1.D0+e1*cos(f1))*sin(&
     &f1)+(-sin(Om1)*sin(omeg1)+cos(Om1)*cos(i1)*cos(omeg1))*(1.D0+e1)/(&
     &1.D0+e1*cos(f1))**2.D0*sin(f1)**2.D0*e1+(-sin(Om1)*sin(omeg1)+cos(&
     &Om1)*cos(i1)*cos(omeg1))*(1.D0+e1)/(1.D0+e1*cos(f1))*cos(f1)      
                                                                        
      s1 = (sin(Om1)*cos(omeg1)+cos(Om1)*cos(i1)*sin(omeg1))*q1/(1.D0+e1&
     &*cos(f1))**2.D0*cos(f1)*e1*sin(f1)-2.D0*(sin(Om1)*cos(omeg1)+cos(O&
     &m1)*cos(i1)*sin(omeg1))*q1*(1.D0+e1)/(1.D0+e1*cos(f1))**3.D0*cos(f&
     &1)**2.D0*e1*sin(f1)+2.D0*(sin(Om1)*cos(omeg1)+cos(Om1)*cos(i1)*sin&
     &(omeg1))*q1*(1.D0+e1)/(1.D0+e1*cos(f1))**2.D0*cos(f1)*sin(f1)-(sin&
     &(Om1)*cos(omeg1)+cos(Om1)*cos(i1)*sin(omeg1))*q1/(1.D0+e1*cos(f1))&
     &*sin(f1)                                                          
      dtau1de1(2) = s1+(-sin(Om1)*sin(omeg1)+cos(Om1)*cos(i1)*cos(omeg1))&
     &*q1/(1.D0+e1*cos(f1))**2.D0*sin(f1)**2.D0*e1-2.D0*(-sin(Om1)*sin(o&
     &meg1)+cos(Om1)*cos(i1)*cos(omeg1))*q1*(1.D0+e1)/(1.D0+e1*cos(f1))*&
     &*3.D0*sin(f1)**2.D0*e1*cos(f1)+(-sin(Om1)*sin(omeg1)+cos(Om1)*cos(&
     &i1)*cos(omeg1))*q1*(1.D0+e1)/(1.D0+e1*cos(f1))**2.D0*sin(f1)**2.D0&
     &+(-sin(Om1)*sin(omeg1)+cos(Om1)*cos(i1)*cos(omeg1))*q1/(1.D0+e1*co&
     &s(f1))*cos(f1)-(-sin(Om1)*sin(omeg1)+cos(Om1)*cos(i1)*cos(omeg1))*&
     &q1*(1.D0+e1)/(1.D0+e1*cos(f1))**2.D0*cos(f1)**2.D0                
                                                                        
      dtau1di1(2) = -cos(Om1)*sin(i1)*sin(omeg1)*q1*(1.D0+e1)/(1.D0+e1*co&
     &s(f1))**2.D0*cos(f1)*e1*sin(f1)+cos(Om1)*sin(i1)*sin(omeg1)*q1*(1.&
     &D0+e1)/(1.D0+e1*cos(f1))*sin(f1)-cos(Om1)*sin(i1)*cos(omeg1)*q1*(1&
     &.D0+e1)/(1.D0+e1*cos(f1))**2.D0*sin(f1)**2.D0*e1-cos(Om1)*sin(i1)*&
     &cos(omeg1)*q1*(1.D0+e1)/(1.D0+e1*cos(f1))*cos(f1)                 
                                                                        
      dtau1dOm1(2) = (cos(Om1)*cos(omeg1)-sin(Om1)*cos(i1)*sin(omeg1))*q1&
     &*(1.D0+e1)/(1.D0+e1*cos(f1))**2.D0*cos(f1)*e1*sin(f1)-(cos(Om1)*co&
     &s(omeg1)-sin(Om1)*cos(i1)*sin(omeg1))*q1*(1.D0+e1)/(1.D0+e1*cos(f1&
     &))*sin(f1)+(-cos(Om1)*sin(omeg1)-sin(Om1)*cos(i1)*cos(omeg1))*q1*(&
     &1.D0+e1)/(1.D0+e1*cos(f1))**2.D0*sin(f1)**2.D0*e1+(-cos(Om1)*sin(o&
     &meg1)-sin(Om1)*cos(i1)*cos(omeg1))*q1*(1.D0+e1)/(1.D0+e1*cos(f1))*&
     &cos(f1)                                                           
                                                                        
      dtau1domeg1(2) = (-sin(Om1)*sin(omeg1)+cos(Om1)*cos(i1)*cos(omeg1))&
     &*q1*(1.D0+e1)/(1.D0+e1*cos(f1))**2.D0*cos(f1)*e1*sin(f1)-(-sin(Om1&
     &)*sin(omeg1)+cos(Om1)*cos(i1)*cos(omeg1))*q1*(1.D0+e1)/(1.D0+e1*co&
     &s(f1))*sin(f1)+(-sin(Om1)*cos(omeg1)-cos(Om1)*cos(i1)*sin(omeg1))*&
     &q1*(1.D0+e1)/(1.D0+e1*cos(f1))**2.D0*sin(f1)**2.D0*e1+(-sin(Om1)*c&
     &os(omeg1)-cos(Om1)*cos(i1)*sin(omeg1))*q1*(1.D0+e1)/(1.D0+e1*cos(f&
     &1))*cos(f1)                                                       
                                                                        
      s1 = 2.D0*(sin(Om1)*cos(omeg1)+cos(Om1)*cos(i1)*sin(omeg1))*q1*(1.&
     &D0+e1)/(1.D0+e1*cos(f1))**3.D0*cos(f1)*e1**2*sin(f1)**2.D0-2.D0*(s&
     &in(Om1)*cos(omeg1)+cos(Om1)*cos(i1)*sin(omeg1))*q1*(1.D0+e1)/(1.D0&
     &+e1*cos(f1))**2.D0*sin(f1)**2.D0*e1+(sin(Om1)*cos(omeg1)+cos(Om1)*&
     &cos(i1)*sin(omeg1))*q1*(1.D0+e1)/(1.D0+e1*cos(f1))**2.D0*cos(f1)**&
     &2.D0*e1                                                           
      dtau1df1(2) = s1-(sin(Om1)*cos(omeg1)+cos(Om1)*cos(i1)*sin(omeg1))*&
     &q1*(1.D0+e1)/(1.D0+e1*cos(f1))*cos(f1)+2.D0*(-sin(Om1)*sin(omeg1)+&
     &cos(Om1)*cos(i1)*cos(omeg1))*q1*(1.D0+e1)/(1.D0+e1*cos(f1))**3.D0*&
     &sin(f1)**3.D0*e1**2+3.D0*(-sin(Om1)*sin(omeg1)+cos(Om1)*cos(i1)*co&
     &s(omeg1))*q1*(1.D0+e1)/(1.D0+e1*cos(f1))**2.D0*cos(f1)*e1*sin(f1)-&
     &(-sin(Om1)*sin(omeg1)+cos(Om1)*cos(i1)*cos(omeg1))*q1*(1.D0+e1)/(1&
     &.D0+e1*cos(f1))*sin(f1)                                           
                                                                        
      dtau1dq1(3) = sin(i1)*sin(omeg1)*(1.D0+e1)/(1.D0+e1*cos(f1))**2.D0*&
     &cos(f1)*e1*sin(f1)-sin(i1)*sin(omeg1)*(1.D0+e1)/(1.D0+e1*cos(f1))*&
     &sin(f1)+sin(i1)*cos(omeg1)*(1.D0+e1)/(1.D0+e1*cos(f1))**2.D0*sin(f&
     &1)**2.D0*e1+sin(i1)*cos(omeg1)*(1.D0+e1)/(1.D0+e1*cos(f1))*cos(f1)
                                                                        
      dtau1de1(3) = sin(i1)*sin(omeg1)*q1/(1.D0+e1*cos(f1))**2.D0*cos(f1)&
     &*e1*sin(f1)-2.D0*sin(i1)*sin(omeg1)*q1*(1.D0+e1)/(1.D0+e1*cos(f1))&
     &**3.D0*cos(f1)**2.D0*e1*sin(f1)+2.D0*sin(i1)*sin(omeg1)*q1*(1.D0+e&
     &1)/(1.D0+e1*cos(f1))**2.D0*cos(f1)*sin(f1)-sin(i1)*sin(omeg1)*q1/(&
     &1.D0+e1*cos(f1))*sin(f1)+sin(i1)*cos(omeg1)*q1/(1.D0+e1*cos(f1))**&
     &2.D0*sin(f1)**2.D0*e1-2.D0*sin(i1)*cos(omeg1)*q1*(1.D0+e1)/(1.D0+e&
     &1*cos(f1))**3.D0*sin(f1)**2.D0*e1*cos(f1)+sin(i1)*cos(omeg1)*q1*(1&
     &.D0+e1)/(1.D0+e1*cos(f1))**2.D0*sin(f1)**2.D0+sin(i1)*cos(omeg1)*q&
     &1/(1.D0+e1*cos(f1))*cos(f1)-sin(i1)*cos(omeg1)*q1*(1.D0+e1)/(1.D0+&
     &e1*cos(f1))**2.D0*cos(f1)**2.D0                                   
                                                                        
      dtau1di1(3) = cos(i1)*sin(omeg1)*q1*(1.D0+e1)/(1.D0+e1*cos(f1))**2.&
     &D0*cos(f1)*e1*sin(f1)-cos(i1)*sin(omeg1)*q1*(1.D0+e1)/(1.D0+e1*cos&
     &(f1))*sin(f1)+cos(i1)*cos(omeg1)*q1*(1.D0+e1)/(1.D0+e1*cos(f1))**2&
     &.D0*sin(f1)**2.D0*e1+cos(i1)*cos(omeg1)*q1*(1.D0+e1)/(1.D0+e1*cos(&
     &f1))*cos(f1)                                                      
                                                                        
      dtau1dOm1(3) = 0.D0 
                                                                        
      dtau1domeg1(3) = sin(i1)*cos(omeg1)*q1*(1.D0+e1)/(1.D0+e1*cos(f1))*&
     &*2.D0*cos(f1)*e1*sin(f1)-sin(i1)*cos(omeg1)*q1*(1.D0+e1)/(1.D0+e1*&
     &cos(f1))*sin(f1)-sin(i1)*sin(omeg1)*q1*(1.D0+e1)/(1.D0+e1*cos(f1))&
     &**2.D0*sin(f1)**2.D0*e1-sin(i1)*sin(omeg1)*q1*(1.D0+e1)/(1.D0+e1*c&
     &os(f1))*cos(f1)                                                   
                                                                        
      dtau1df1(3) = 2.D0*sin(i1)*sin(omeg1)*q1*(1.D0+e1)/(1.D0+e1*cos(f1)&
     &)**3.D0*cos(f1)*e1**2*sin(f1)**2.D0-2.D0*sin(i1)*sin(omeg1)*q1*(1.&
     &D0+e1)/(1.D0+e1*cos(f1))**2.D0*sin(f1)**2.D0*e1+sin(i1)*sin(omeg1)&
     &*q1*(1.D0+e1)/(1.D0+e1*cos(f1))**2.D0*cos(f1)**2.D0*e1-sin(i1)*sin&
     &(omeg1)*q1*(1.D0+e1)/(1.D0+e1*cos(f1))*cos(f1)+2.D0*sin(i1)*cos(om&
     &eg1)*q1*(1.D0+e1)/(1.D0+e1*cos(f1))**3.D0*sin(f1)**3.D0*e1**2+3.D0&
     &*sin(i1)*cos(omeg1)*q1*(1.D0+e1)/(1.D0+e1*cos(f1))**2.D0*cos(f1)*e&
     &1*sin(f1)-sin(i1)*cos(omeg1)*q1*(1.D0+e1)/(1.D0+e1*cos(f1))*sin(f1)

 ! --------------------------------------------------------------------
      dtau2dq2(1) = (cos(Om2)*cos(omeg2)-sin(Om2)*cos(i2)*sin(omeg2))*(1.&
     &D0+e2)/(1.D0+e2*cos(f2))**2.D0*cos(f2)*e2*sin(f2)-(cos(Om2)*cos(om&
     &eg2)-sin(Om2)*cos(i2)*sin(omeg2))*(1.D0+e2)/(1.D0+e2*cos(f2))*sin(&
     &f2)+(-cos(Om2)*sin(omeg2)-sin(Om2)*cos(i2)*cos(omeg2))*(1.D0+e2)/(&
     &1.D0+e2*cos(f2))**2.D0*sin(f2)**2.D0*e2+(-cos(Om2)*sin(omeg2)-sin(&
     &Om2)*cos(i2)*cos(omeg2))*(1.D0+e2)/(1.D0+e2*cos(f2))*cos(f2)      
                                                                        
      s1 = (cos(Om2)*cos(omeg2)-sin(Om2)*cos(i2)*sin(omeg2))*q2/(1.D0+e2&
     &*cos(f2))**2.D0*cos(f2)*e2*sin(f2)-2.D0*(cos(Om2)*cos(omeg2)-sin(O&
     &m2)*cos(i2)*sin(omeg2))*q2*(1.D0+e2)/(1.D0+e2*cos(f2))**3.D0*cos(f&
     &2)**2.D0*e2*sin(f2)+2.D0*(cos(Om2)*cos(omeg2)-sin(Om2)*cos(i2)*sin&
     &(omeg2))*q2*(1.D0+e2)/(1.D0+e2*cos(f2))**2.D0*cos(f2)*sin(f2)-(cos&
     &(Om2)*cos(omeg2)-sin(Om2)*cos(i2)*sin(omeg2))*q2/(1.D0+e2*cos(f2))&
     &*sin(f2)                                                          
      dtau2de2(1) = s1+(-cos(Om2)*sin(omeg2)-sin(Om2)*cos(i2)*cos(omeg2))&
     &*q2/(1.D0+e2*cos(f2))**2.D0*sin(f2)**2.D0*e2-2.D0*(-cos(Om2)*sin(o&
     &meg2)-sin(Om2)*cos(i2)*cos(omeg2))*q2*(1.D0+e2)/(1.D0+e2*cos(f2))*&
     &*3.D0*sin(f2)**2.D0*e2*cos(f2)+(-cos(Om2)*sin(omeg2)-sin(Om2)*cos(&
     &i2)*cos(omeg2))*q2*(1.D0+e2)/(1.D0+e2*cos(f2))**2.D0*sin(f2)**2.D0&
     &+(-cos(Om2)*sin(omeg2)-sin(Om2)*cos(i2)*cos(omeg2))*q2/(1.D0+e2*co&
     &s(f2))*cos(f2)-(-cos(Om2)*sin(omeg2)-sin(Om2)*cos(i2)*cos(omeg2))*&
     &q2*(1.D0+e2)/(1.D0+e2*cos(f2))**2.D0*cos(f2)**2.D0                
                                                                        
      dtau2di2(1) = sin(Om2)*sin(i2)*sin(omeg2)*q2*(1.D0+e2)/(1.D0+e2*cos&
     &(f2))**2.D0*cos(f2)*e2*sin(f2)-sin(Om2)*sin(i2)*sin(omeg2)*q2*(1.D&
     &0+e2)/(1.D0+e2*cos(f2))*sin(f2)+sin(Om2)*sin(i2)*cos(omeg2)*q2*(1.&
     &D0+e2)/(1.D0+e2*cos(f2))**2.D0*sin(f2)**2.D0*e2+sin(Om2)*sin(i2)*c&
     &os(omeg2)*q2*(1.D0+e2)/(1.D0+e2*cos(f2))*cos(f2)                  
                                                                        
      dtau2dOm2(1) = (-sin(Om2)*cos(omeg2)-cos(Om2)*cos(i2)*sin(omeg2))*q&
     &2*(1.D0+e2)/(1.D0+e2*cos(f2))**2.D0*cos(f2)*e2*sin(f2)-(-sin(Om2)*&
     &cos(omeg2)-cos(Om2)*cos(i2)*sin(omeg2))*q2*(1.D0+e2)/(1.D0+e2*cos(&
     &f2))*sin(f2)+(sin(Om2)*sin(omeg2)-cos(Om2)*cos(i2)*cos(omeg2))*q2*&
     &(1.D0+e2)/(1.D0+e2*cos(f2))**2.D0*sin(f2)**2.D0*e2+(sin(Om2)*sin(o&
     &meg2)-cos(Om2)*cos(i2)*cos(omeg2))*q2*(1.D0+e2)/(1.D0+e2*cos(f2))*&
     &cos(f2)                                                           
                                                                        
      dtau2domeg2(1) = (-cos(Om2)*sin(omeg2)-sin(Om2)*cos(i2)*cos(omeg2))&
     &*q2*(1.D0+e2)/(1.D0+e2*cos(f2))**2.D0*cos(f2)*e2*sin(f2)-(-cos(Om2&
     &)*sin(omeg2)-sin(Om2)*cos(i2)*cos(omeg2))*q2*(1.D0+e2)/(1.D0+e2*co&
     &s(f2))*sin(f2)+(-cos(Om2)*cos(omeg2)+sin(Om2)*cos(i2)*sin(omeg2))*&
     &q2*(1.D0+e2)/(1.D0+e2*cos(f2))**2.D0*sin(f2)**2.D0*e2+(-cos(Om2)*c&
     &os(omeg2)+sin(Om2)*cos(i2)*sin(omeg2))*q2*(1.D0+e2)/(1.D0+e2*cos(f&
     &2))*cos(f2)                                                       
                                                                        
      s1 = 2.D0*(cos(Om2)*cos(omeg2)-sin(Om2)*cos(i2)*sin(omeg2))*q2*(1.&
     &D0+e2)/(1.D0+e2*cos(f2))**3.D0*cos(f2)*e2**2*sin(f2)**2.D0-2.D0*(c&
     &os(Om2)*cos(omeg2)-sin(Om2)*cos(i2)*sin(omeg2))*q2*(1.D0+e2)/(1.D0&
     &+e2*cos(f2))**2.D0*sin(f2)**2.D0*e2+(cos(Om2)*cos(omeg2)-sin(Om2)*&
     &cos(i2)*sin(omeg2))*q2*(1.D0+e2)/(1.D0+e2*cos(f2))**2.D0*cos(f2)**&
     &2.D0*e2                                                           
      dtau2df2(1) = s1-(cos(Om2)*cos(omeg2)-sin(Om2)*cos(i2)*sin(omeg2))*&
     &q2*(1.D0+e2)/(1.D0+e2*cos(f2))*cos(f2)+2.D0*(-cos(Om2)*sin(omeg2)-&
     &sin(Om2)*cos(i2)*cos(omeg2))*q2*(1.D0+e2)/(1.D0+e2*cos(f2))**3.D0*&
     &sin(f2)**3.D0*e2**2+3.D0*(-cos(Om2)*sin(omeg2)-sin(Om2)*cos(i2)*co&
     &s(omeg2))*q2*(1.D0+e2)/(1.D0+e2*cos(f2))**2.D0*cos(f2)*e2*sin(f2)-&
     &(-cos(Om2)*sin(omeg2)-sin(Om2)*cos(i2)*cos(omeg2))*q2*(1.D0+e2)/(1&
     &.D0+e2*cos(f2))*sin(f2)                                           
                                                                        
      dtau2dq2(2) = (sin(Om2)*cos(omeg2)+cos(Om2)*cos(i2)*sin(omeg2))*(1.&
     &D0+e2)/(1.D0+e2*cos(f2))**2.D0*cos(f2)*e2*sin(f2)-(sin(Om2)*cos(om&
     &eg2)+cos(Om2)*cos(i2)*sin(omeg2))*(1.D0+e2)/(1.D0+e2*cos(f2))*sin(&
     &f2)+(-sin(Om2)*sin(omeg2)+cos(Om2)*cos(i2)*cos(omeg2))*(1.D0+e2)/(&
     &1.D0+e2*cos(f2))**2.D0*sin(f2)**2.D0*e2+(-sin(Om2)*sin(omeg2)+cos(&
     &Om2)*cos(i2)*cos(omeg2))*(1.D0+e2)/(1.D0+e2*cos(f2))*cos(f2)      
                                                                        
      s1 = (sin(Om2)*cos(omeg2)+cos(Om2)*cos(i2)*sin(omeg2))*q2/(1.D0+e2&
     &*cos(f2))**2.D0*cos(f2)*e2*sin(f2)-2.D0*(sin(Om2)*cos(omeg2)+cos(O&
     &m2)*cos(i2)*sin(omeg2))*q2*(1.D0+e2)/(1.D0+e2*cos(f2))**3.D0*cos(f&
     &2)**2.D0*e2*sin(f2)+2.D0*(sin(Om2)*cos(omeg2)+cos(Om2)*cos(i2)*sin&
     &(omeg2))*q2*(1.D0+e2)/(1.D0+e2*cos(f2))**2.D0*cos(f2)*sin(f2)-(sin&
     &(Om2)*cos(omeg2)+cos(Om2)*cos(i2)*sin(omeg2))*q2/(1.D0+e2*cos(f2))&
     &*sin(f2)                                                          
      dtau2de2(2) = s1+(-sin(Om2)*sin(omeg2)+cos(Om2)*cos(i2)*cos(omeg2))&
     &*q2/(1.D0+e2*cos(f2))**2.D0*sin(f2)**2.D0*e2-2.D0*(-sin(Om2)*sin(o&
     &meg2)+cos(Om2)*cos(i2)*cos(omeg2))*q2*(1.D0+e2)/(1.D0+e2*cos(f2))*&
     &*3.D0*sin(f2)**2.D0*e2*cos(f2)+(-sin(Om2)*sin(omeg2)+cos(Om2)*cos(&
     &i2)*cos(omeg2))*q2*(1.D0+e2)/(1.D0+e2*cos(f2))**2.D0*sin(f2)**2.D0&
     &+(-sin(Om2)*sin(omeg2)+cos(Om2)*cos(i2)*cos(omeg2))*q2/(1.D0+e2*co&
     &s(f2))*cos(f2)-(-sin(Om2)*sin(omeg2)+cos(Om2)*cos(i2)*cos(omeg2))*&
     &q2*(1.D0+e2)/(1.D0+e2*cos(f2))**2.D0*cos(f2)**2.D0                
                                                                        
      dtau2di2(2) = -cos(Om2)*sin(i2)*sin(omeg2)*q2*(1.D0+e2)/(1.D0+e2*co&
     &s(f2))**2.D0*cos(f2)*e2*sin(f2)+cos(Om2)*sin(i2)*sin(omeg2)*q2*(1.&
     &D0+e2)/(1.D0+e2*cos(f2))*sin(f2)-cos(Om2)*sin(i2)*cos(omeg2)*q2*(1&
     &.D0+e2)/(1.D0+e2*cos(f2))**2.D0*sin(f2)**2.D0*e2-cos(Om2)*sin(i2)*&
     &cos(omeg2)*q2*(1.D0+e2)/(1.D0+e2*cos(f2))*cos(f2)                 
                                                                        
      dtau2dOm2(2) = (cos(Om2)*cos(omeg2)-sin(Om2)*cos(i2)*sin(omeg2))*q2&
     &*(1.D0+e2)/(1.D0+e2*cos(f2))**2.D0*cos(f2)*e2*sin(f2)-(cos(Om2)*co&
     &s(omeg2)-sin(Om2)*cos(i2)*sin(omeg2))*q2*(1.D0+e2)/(1.D0+e2*cos(f2&
     &))*sin(f2)+(-cos(Om2)*sin(omeg2)-sin(Om2)*cos(i2)*cos(omeg2))*q2*(&
     &1.D0+e2)/(1.D0+e2*cos(f2))**2.D0*sin(f2)**2.D0*e2+(-cos(Om2)*sin(o&
     &meg2)-sin(Om2)*cos(i2)*cos(omeg2))*q2*(1.D0+e2)/(1.D0+e2*cos(f2))*&
     &cos(f2)                                                           
                                                                        
      dtau2domeg2(2) = (-sin(Om2)*sin(omeg2)+cos(Om2)*cos(i2)*cos(omeg2))&
     &*q2*(1.D0+e2)/(1.D0+e2*cos(f2))**2.D0*cos(f2)*e2*sin(f2)-(-sin(Om2&
     &)*sin(omeg2)+cos(Om2)*cos(i2)*cos(omeg2))*q2*(1.D0+e2)/(1.D0+e2*co&
     &s(f2))*sin(f2)+(-sin(Om2)*cos(omeg2)-cos(Om2)*cos(i2)*sin(omeg2))*&
     &q2*(1.D0+e2)/(1.D0+e2*cos(f2))**2.D0*sin(f2)**2.D0*e2+(-sin(Om2)*c&
     &os(omeg2)-cos(Om2)*cos(i2)*sin(omeg2))*q2*(1.D0+e2)/(1.D0+e2*cos(f&
     &2))*cos(f2)                                                       
                                                                        
      s1 = 2.D0*(sin(Om2)*cos(omeg2)+cos(Om2)*cos(i2)*sin(omeg2))*q2*(1.&
     &D0+e2)/(1.D0+e2*cos(f2))**3.D0*cos(f2)*e2**2*sin(f2)**2.D0-2.D0*(s&
     &in(Om2)*cos(omeg2)+cos(Om2)*cos(i2)*sin(omeg2))*q2*(1.D0+e2)/(1.D0&
     &+e2*cos(f2))**2.D0*sin(f2)**2.D0*e2+(sin(Om2)*cos(omeg2)+cos(Om2)*&
     &cos(i2)*sin(omeg2))*q2*(1.D0+e2)/(1.D0+e2*cos(f2))**2.D0*cos(f2)**&
     &2.D0*e2                                                           
      dtau2df2(2) = s1-(sin(Om2)*cos(omeg2)+cos(Om2)*cos(i2)*sin(omeg2))*q&
     &2*(1.D0+e2)/(1.D0+e2*cos(f2))*cos(f2)+2.D0*(-sin(Om2)*sin(omeg2)+c&
     &os(Om2)*cos(i2)*cos(omeg2))*q2*(1.D0+e2)/(1.D0+e2*cos(f2))**3.D0*s&
     &in(f2)**3.D0*e2**2+3.D0*(-sin(Om2)*sin(omeg2)+cos(Om2)*cos(i2)*cos&
     &(omeg2))*q2*(1.D0+e2)/(1.D0+e2*cos(f2))**2.D0*cos(f2)*e2*sin(f2)-(&
     &-sin(Om2)*sin(omeg2)+cos(Om2)*cos(i2)*cos(omeg2))*q2*(1.D0+e2)/(1.&
     &D0+e2*cos(f2))*sin(f2)                                            
                                                                        
      dtau2dq2(3) = sin(i2)*sin(omeg2)*(1.D0+e2)/(1.D0+e2*cos(f2))**2.D0*&
     &cos(f2)*e2*sin(f2)-sin(i2)*sin(omeg2)*(1.D0+e2)/(1.D0+e2*cos(f2))*&
     &sin(f2)+sin(i2)*cos(omeg2)*(1.D0+e2)/(1.D0+e2*cos(f2))**2.D0*sin(f&
     &2)**2.D0*e2+sin(i2)*cos(omeg2)*(1.D0+e2)/(1.D0+e2*cos(f2))*cos(f2)
                                                                        
      dtau2de2(3) = sin(i2)*sin(omeg2)*q2/(1.D0+e2*cos(f2))**2.D0*cos(f2)&
     &*e2*sin(f2)-2.D0*sin(i2)*sin(omeg2)*q2*(1.D0+e2)/(1.D0+e2*cos(f2))&
     &**3.D0*cos(f2)**2.D0*e2*sin(f2)+2.D0*sin(i2)*sin(omeg2)*q2*(1.D0+e&
     &2)/(1.D0+e2*cos(f2))**2.D0*cos(f2)*sin(f2)-sin(i2)*sin(omeg2)*q2/(&
     &1.D0+e2*cos(f2))*sin(f2)+sin(i2)*cos(omeg2)*q2/(1.D0+e2*cos(f2))**&
     &2.D0*sin(f2)**2.D0*e2-2.D0*sin(i2)*cos(omeg2)*q2*(1.D0+e2)/(1.D0+e&
     &2*cos(f2))**3.D0*sin(f2)**2.D0*e2*cos(f2)+sin(i2)*cos(omeg2)*q2*(1&
     &.D0+e2)/(1.D0+e2*cos(f2))**2.D0*sin(f2)**2.D0+sin(i2)*cos(omeg2)*q&
     &2/(1.D0+e2*cos(f2))*cos(f2)-sin(i2)*cos(omeg2)*q2*(1.D0+e2)/(1.D0+&
     &e2*cos(f2))**2.D0*cos(f2)**2.D0                                   
                                                                        
      dtau2di2(3) = cos(i2)*sin(omeg2)*q2*(1.D0+e2)/(1.D0+e2*cos(f2))**2.&
     &D0*cos(f2)*e2*sin(f2)-cos(i2)*sin(omeg2)*q2*(1.D0+e2)/(1.D0+e2*cos&
     &(f2))*sin(f2)+cos(i2)*cos(omeg2)*q2*(1.D0+e2)/(1.D0+e2*cos(f2))**2&
     &.D0*sin(f2)**2.D0*e2+cos(i2)*cos(omeg2)*q2*(1.D0+e2)/(1.D0+e2*cos(&
     &f2))*cos(f2)                                                      
                                                                        
      dtau2dOm2(3) = 0.D0 
                                                                        
      dtau2domeg2(3) = sin(i2)*cos(omeg2)*q2*(1.D0+e2)/(1.D0+e2*cos(f2))*&
     &*2.D0*cos(f2)*e2*sin(f2)-sin(i2)*cos(omeg2)*q2*(1.D0+e2)/(1.D0+e2*&
     &cos(f2))*sin(f2)-sin(i2)*sin(omeg2)*q2*(1.D0+e2)/(1.D0+e2*cos(f2))&
     &**2.D0*sin(f2)**2.D0*e2-sin(i2)*sin(omeg2)*q2*(1.D0+e2)/(1.D0+e2*c&
     &os(f2))*cos(f2)                                                   
                                                                        
      dtau2df2(3) = 2.D0*sin(i2)*sin(omeg2)*q2*(1.D0+e2)/(1.D0+e2*cos(f2)&
     &)**3.D0*cos(f2)*e2**2*sin(f2)**2.D0-2.D0*sin(i2)*sin(omeg2)*q2*(1.&
     &D0+e2)/(1.D0+e2*cos(f2))**2.D0*sin(f2)**2.D0*e2+sin(i2)*sin(omeg2)&
     &*q2*(1.D0+e2)/(1.D0+e2*cos(f2))**2.D0*cos(f2)**2.D0*e2-sin(i2)*sin&
     &(omeg2)*q2*(1.D0+e2)/(1.D0+e2*cos(f2))*cos(f2)+2.D0*sin(i2)*cos(om&
     &eg2)*q2*(1.D0+e2)/(1.D0+e2*cos(f2))**3.D0*sin(f2)**3.D0*e2**2+3.D0&
     &*sin(i2)*cos(omeg2)*q2*(1.D0+e2)/(1.D0+e2*cos(f2))**2.D0*cos(f2)*e&
     &2*sin(f2)-sin(i2)*cos(omeg2)*q2*(1.D0+e2)/(1.D0+e2*cos(f2))*sin(f2)

! ====================================================================
  dtau1versdq1(1)=1.d0/ltau1*dtau1dq1(1) - 1.d0/(ltau1**3)* &
       & DOT_PRODUCT(tau1,dtau1dq1)*tau1(1)
  dtau1versdq1(2)=1.d0/ltau1*dtau1dq1(2) - 1.d0/(ltau1**3)* &
       & DOT_PRODUCT(tau1,dtau1dq1)*tau1(2)
  dtau1versdq1(3)=1.d0/ltau1*dtau1dq1(3) - 1.d0/(ltau1**3)* &
       & DOT_PRODUCT(tau1,dtau1dq1)*tau1(3)

  dtau1versde1(1)=1.d0/ltau1*dtau1de1(1) - 1.d0/(ltau1**3)* &
       & DOT_PRODUCT(tau1,dtau1de1)*tau1(1)
  dtau1versde1(2)=1.d0/ltau1*dtau1de1(2) - 1.d0/(ltau1**3)* &
       & DOT_PRODUCT(tau1,dtau1de1)*tau1(2)
  dtau1versde1(3)=1.d0/ltau1*dtau1de1(3) - 1.d0/(ltau1**3)* &
       & DOT_PRODUCT(tau1,dtau1de1)*tau1(3)

  dtau1versdi1(1)=1.d0/ltau1*dtau1di1(1) - 1.d0/(ltau1**3)* &
       & DOT_PRODUCT(tau1,dtau1di1)*tau1(1)
  dtau1versdi1(2)=1.d0/ltau1*dtau1di1(2) - 1.d0/(ltau1**3)* &
       & DOT_PRODUCT(tau1,dtau1di1)*tau1(2)
  dtau1versdi1(3)=1.d0/ltau1*dtau1di1(3) - 1.d0/(ltau1**3)* &
       & DOT_PRODUCT(tau1,dtau1di1)*tau1(3)

  dtau1versdOm1(1)=1.d0/ltau1*dtau1dOm1(1) - 1.d0/(ltau1**3)* &
       & DOT_PRODUCT(tau1,dtau1dOm1)*tau1(1)
  dtau1versdOm1(2)=1.d0/ltau1*dtau1dOm1(2) - 1.d0/(ltau1**3)* &
       & DOT_PRODUCT(tau1,dtau1dOm1)*tau1(2)
  dtau1versdOm1(3)=1.d0/ltau1*dtau1dOm1(3) - 1.d0/(ltau1**3)* &
       & DOT_PRODUCT(tau1,dtau1dOm1)*tau1(3)

  dtau1versdomeg1(1)=1.d0/ltau1*dtau1domeg1(1) - 1.d0/(ltau1**3)* &
       & DOT_PRODUCT(tau1,dtau1domeg1)*tau1(1)
  dtau1versdomeg1(2)=1.d0/ltau1*dtau1domeg1(2) - 1.d0/(ltau1**3)* &
       & DOT_PRODUCT(tau1,dtau1domeg1)*tau1(2)
  dtau1versdomeg1(3)=1.d0/ltau1*dtau1domeg1(3) - 1.d0/(ltau1**3)* &
       & DOT_PRODUCT(tau1,dtau1domeg1)*tau1(3)

  dtau1versdf1(1)=1.d0/ltau1*dtau1df1(1) - 1.d0/(ltau1**3)* &
       & DOT_PRODUCT(tau1,dtau1df1)*tau1(1)
  dtau1versdf1(2)=1.d0/ltau1*dtau1df1(2) - 1.d0/(ltau1**3)* &
       & DOT_PRODUCT(tau1,dtau1df1)*tau1(2)
  dtau1versdf1(3)=1.d0/ltau1*dtau1df1(3) - 1.d0/(ltau1**3)* &
       & DOT_PRODUCT(tau1,dtau1df1)*tau1(3)

! ------------------------------------------------------------
  dtau2versdq2(1)=1.d0/ltau2*dtau2dq2(1) - 1.d0/(ltau2**3)* &
       & DOT_PRODUCT(tau2,dtau2dq2)*tau2(1)
  dtau2versdq2(2)=1.d0/ltau2*dtau2dq2(2) - 1.d0/(ltau2**3)* &
       & DOT_PRODUCT(tau2,dtau2dq2)*tau2(2)
  dtau2versdq2(3)=1.d0/ltau2*dtau2dq2(3) - 1.d0/(ltau2**3)* &
       & DOT_PRODUCT(tau2,dtau2dq2)*tau2(3)

  dtau2versde2(1)=1.d0/ltau2*dtau2de2(1) - 1.d0/(ltau2**3)* &
       & DOT_PRODUCT(tau2,dtau2de2)*tau2(1)
  dtau2versde2(2)=1.d0/ltau2*dtau2de2(2) - 1.d0/(ltau2**3)* &
       & DOT_PRODUCT(tau2,dtau2de2)*tau2(2)
  dtau2versde2(3)=1.d0/ltau2*dtau2de2(3) - 1.d0/(ltau2**3)* &
       & DOT_PRODUCT(tau2,dtau2de2)*tau2(3)

  dtau2versdi2(1)=1.d0/ltau2*dtau2di2(1) - 1.d0/(ltau2**3)* &
       & DOT_PRODUCT(tau2,dtau2di2)*tau2(1)
  dtau2versdi2(2)=1.d0/ltau2*dtau2di2(2) - 1.d0/(ltau2**3)* &
       & DOT_PRODUCT(tau2,dtau2di2)*tau2(2)
  dtau2versdi2(3)=1.d0/ltau2*dtau2di2(3) - 1.d0/(ltau2**3)* &
       & DOT_PRODUCT(tau2,dtau2di2)*tau2(3)

  dtau2versdOm2(1)=1.d0/ltau2*dtau2dOm2(1) - 1.d0/(ltau2**3)* &
       & DOT_PRODUCT(tau2,dtau2dOm2)*tau2(1)
  dtau2versdOm2(2)=1.d0/ltau2*dtau2dOm2(2) - 1.d0/(ltau2**3)* &
       & DOT_PRODUCT(tau2,dtau2dOm2)*tau2(2)
  dtau2versdOm2(3)=1.d0/ltau2*dtau2dOm2(3) - 1.d0/(ltau2**3)* &
       & DOT_PRODUCT(tau2,dtau2dOm2)*tau2(3)

  dtau2versdomeg2(1)=1.d0/ltau2*dtau2domeg2(1) - 1.d0/(ltau2**3)* &
       & DOT_PRODUCT(tau2,dtau2domeg2)*tau2(1)
  dtau2versdomeg2(2)=1.d0/ltau2*dtau2domeg2(2) - 1.d0/(ltau2**3)* &
       & DOT_PRODUCT(tau2,dtau2domeg2)*tau2(2)
  dtau2versdomeg2(3)=1.d0/ltau2*dtau2domeg2(3) - 1.d0/(ltau2**3)* &
       & DOT_PRODUCT(tau2,dtau2domeg2)*tau2(3)

  dtau2versdf2(1)=1.d0/ltau2*dtau2df2(1) - 1.d0/(ltau2**3)* &
       & DOT_PRODUCT(tau2,dtau2df2)*tau2(1)
  dtau2versdf2(2)=1.d0/ltau2*dtau2df2(2) - 1.d0/(ltau2**3)* &
       & DOT_PRODUCT(tau2,dtau2df2)*tau2(2)
  dtau2versdf2(3)=1.d0/ltau2*dtau2df2(3) - 1.d0/(ltau2**3)* &
       & DOT_PRODUCT(tau2,dtau2df2)*tau2(3)

! ===========================================================
  CALL prvec(dtau1versdq1,tau2vers,dtau3vecdq1)
  CALL prvec(dtau1versde1,tau2vers,dtau3vecde1)
  CALL prvec(dtau1versdi1,tau2vers,dtau3vecdi1)
  CALL prvec(dtau1versdOm1,tau2vers,dtau3vecdOm1)
  CALL prvec(dtau1versdomeg1,tau2vers,dtau3vecdomeg1)
  CALL prvec(dtau1versdf1,tau2vers,dtau3vecdf1)

  CALL prvec(tau1vers,dtau2versdq2,dtau3vecdq2)
  CALL prvec(tau1vers,dtau2versde2,dtau3vecde2)
  CALL prvec(tau1vers,dtau2versdi2,dtau3vecdi2)
  CALL prvec(tau1vers,dtau2versdOm2,dtau3vecdOm2)
  CALL prvec(tau1vers,dtau2versdomeg2,dtau3vecdomeg2)
  CALL prvec(tau1vers,dtau2versdf2,dtau3vecdf2)

! ============================================================
  dsint1t2dq1= 1.d0/sint1t2*DOT_PRODUCT(tau3vec,dtau3vecdq1)
  dsint1t2de1= 1.d0/sint1t2*DOT_PRODUCT(tau3vec,dtau3vecde1)
  dsint1t2di1= 1.d0/sint1t2*DOT_PRODUCT(tau3vec,dtau3vecdi1)
  dsint1t2dOm1= 1.d0/sint1t2*DOT_PRODUCT(tau3vec,dtau3vecdOm1)
  dsint1t2domeg1= 1.d0/sint1t2*DOT_PRODUCT(tau3vec,dtau3vecdomeg1)
  dsint1t2df1= 1.d0/sint1t2*DOT_PRODUCT(tau3vec,dtau3vecdf1)

  dsint1t2dq2= 1.d0/sint1t2*DOT_PRODUCT(tau3vec,dtau3vecdq2)
  dsint1t2de2= 1.d0/sint1t2*DOT_PRODUCT(tau3vec,dtau3vecde2)
  dsint1t2di2= 1.d0/sint1t2*DOT_PRODUCT(tau3vec,dtau3vecdi2)
  dsint1t2dOm2= 1.d0/sint1t2*DOT_PRODUCT(tau3vec,dtau3vecdOm2)
  dsint1t2domeg2= 1.d0/sint1t2*DOT_PRODUCT(tau3vec,dtau3vecdomeg2)
  dsint1t2df2= 1.d0/sint1t2*DOT_PRODUCT(tau3vec,dtau3vecdf2)

! ===========  CHAIN RULE ===================================
  dersint1t2dcom(1,1)=dsint1t2dq1+dsint1t2df1*df1_dq1+dsint1t2df2*df2_dq1
  dersint1t2dcom(1,2)=dsint1t2de1+dsint1t2df1*df1_de1+dsint1t2df2*df2_de1
  dersint1t2dcom(1,3)=dsint1t2di1+dsint1t2df1*df1_di1+dsint1t2df2*df2_di1
  dersint1t2dcom(1,4)=dsint1t2dOm1+dsint1t2df1*df1_dOm1+dsint1t2df2*df2_dOm1
  dersint1t2dcom(1,5)=dsint1t2domeg1+dsint1t2df1*df1_domeg1+ &
       & dsint1t2df2*df2_domeg1
  dersint1t2dcom(1,6)=dsint1t2dq2 +dsint1t2df1*df1_dq2+dsint1t2df2*df2_dq2
  dersint1t2dcom(1,7)=dsint1t2de2+dsint1t2df1*df1_de2+dsint1t2df2*df2_de2
  dersint1t2dcom(1,8)=dsint1t2di2+dsint1t2df1*df1_di2+dsint1t2df2*df2_di2
  dersint1t2dcom(1,9)=dsint1t2dOm2+dsint1t2df1*df1_dOm2+dsint1t2df2*df2_dOm2
  dersint1t2dcom(1,10)=dsint1t2domeg2+dsint1t2df1*df1_domeg2+ &
       & dsint1t2df2*df2_domeg2
  tdersint1t2dcom=TRANSPOSE(dersint1t2dcom)
  taurmsvec=MATMUL(dersint1t2dcom,MATMUL(covcom,tdersint1t2dcom)) 
  taurms=SQRT(taurmsvec(1,1))
  
  IF(chk_der) THEN
!     dersint1t2dequ2(1:1,1:5)= &
!         & MATMUL(dersint1t2dcom(1:1,6:10),jaccomeq2(1:5,1:5))
     write(ierrou,*)'derivatives of sint1t2:'
     write(ierrou,100) dersint1t2dcom(1,6:10)
     write(ierrou,*)'-----------------------------------------&
          &-------------------------'
     numerr=numerr+1
  ENDIF
  
100 FORMAT(5(f13.8,1x))

END SUBROUTINE comp_rms_com

! **************************************************************
! RMS of the perihelion distance q, of the afhelion distance Q 
! and of the orbital period T
! **************************************************************
! written by G.F. Gronchi
! **************************************************************
  SUBROUTINE q_Q_T_rms(eq,unc,perih,apoh,period, &
       & perihrms,apohrms,periodrms)
    IMPLICIT NONE
! =====================INTERFACE=======================
    TYPE(orbit_elem), INTENT(IN) :: eq !EQU elements of the asteroid
    TYPE(orb_uncert), INTENT(IN) :: unc
    DOUBLE PRECISION, INTENT(OUT) :: perih,apoh,period
    DOUBLE PRECISION, INTENT(OUT) :: perihrms,apohrms,periodrms
! =========================================================  
    DOUBLE PRECISION, DIMENSION(6,6) :: geq
    DOUBLE PRECISION :: vsize
    DOUBLE PRECISION :: aa,ee,hh,kk
    DOUBLE PRECISION, DIMENSION(1,3) :: derperih,derapoh
    DOUBLE PRECISION :: derperiod
    DOUBLE PRECISION, DIMENSION(3,1) :: tderperih,tderapoh,tderperiod
! auxiliary
    DOUBLE PRECISION, DIMENSION(1,1) :: perihrmsvec,apohrmsvec
    INTEGER :: i,j,h,k,fail_flag     ! Loop indexes and flags
! ============================================================

    aa=eq%coord(1)
    hh=eq%coord(2)
    kk=eq%coord(3)
    ee=sqrt(hh**2+kk**2)

    perih = aa*(1.d0-ee)
    apoh = aa*(1.d0+ee)
    derperih(1,1) = 1.d0-ee
    derapoh(1,1) = 1.d0+ee
    IF(ee.gt.0.d0)THEN
       derperih(1,2) = -aa*hh/ee
       derperih(1,3) = -aa*kk/ee
       derapoh(1,2) = aa*hh/ee
       derapoh(1,3) = aa*kk/ee
    ELSE
       WRITE(ierrou,*)'q_Q_T_rms:ERROR! zero eccentricity',ee
       numerr=numerr+1
    ENDIF
    tderperih=TRANSPOSE(derperih)
    tderapoh=TRANSPOSE(derapoh)

    period=(2.d0*pig/gk)*aa**(3.d0/2.d0)
    derperiod = (3.d0*pig/gk)*sqrt(aa)

! ------ covariance computation ------
    geq = unc%g(1:6,1:6)
    perihrmsvec = MATMUL(derperih,MATMUL(geq(1:3,1:3),tderperih))
    apohrmsvec = MATMUL(derapoh,MATMUL(geq(1:3,1:3),tderapoh))
    perihrms = sqrt(perihrmsvec(1,1))
    apohrms = sqrt(apohrmsvec(1,1))
!
    periodrms = geq(1,1)*derperiod**2

  END SUBROUTINE q_Q_T_rms

! =====================================================================
! written by G.F. Gronchi and G. Tommei (10/2/2006)
! last modified 6/2/2007
SUBROUTINE tau1_tau2(com1,com2,f1,f2,car1min,car2min,tau1,tau2,tau3hat,sint1t2)
  USE ever_pitkin
  TYPE(orbit_elem), INTENT(IN) :: com1,com2 ! cometary elements
                                          ! 1: planet 
                                          ! 2: asteroid
  DOUBLE PRECISION, INTENT(IN) :: f1,f2   ! anomalies at the critical point
                                          ! in [0,2*pi]
  TYPE(orbit_elem),INTENT(OUT) :: car1min,car2min
  DOUBLE PRECISION,DIMENSION(3),INTENT(OUT) :: tau1,tau2,tau3hat
  DOUBLE PRECISION,INTENT(OUT) :: sint1t2
! ----- end interface -----
  INTEGER :: fail_flag
  DOUBLE PRECISION :: q1,e1,i1,Om1,omeg1  ! Planet
  DOUBLE PRECISION :: q2,e2,i2,Om2,omeg2  ! Asteroid/Comet
  DOUBLE PRECISION :: s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12,s13,s14,s15
  DOUBLE PRECISION, DIMENSION(2,2) :: invH 
!====================================================
 DOUBLE PRECISION ::l1,l2,lambda1,lambda2,u1,u2
! to compute derivatives of dmintil
  TYPE(orbit_elem) :: com1min,com2min
  DOUBLE PRECISION :: vsize,ltau1,ltau2 
  DOUBLE PRECISION,DIMENSION(3) :: tau1vers,tau2vers,tau3vec,tau3
! to compute time of passage at perihelion
  DOUBLE PRECISION ::cosn1,sinn1,coso1,sino1,cosi1,sini1
  DOUBLE PRECISION ::cosn2,sinn2,coso2,sino2,cosi2,sini2
  DOUBLE PRECISION,DIMENSION(3) :: x01,y01,x02,y02 
  DOUBLE PRECISION :: r01,v01,r1,mu1,alpha1,psi1,dt1,ptime1
  DOUBLE PRECISION :: r02,v02,r2,mu2,alpha2,psi2,dt2,ptime2
  DOUBLE PRECISION :: psi1_ini,psi2_ini,u1x,u2x
! ==========================================================

! elements
  q1 = com1%coord(1) 
  e1 = com1%coord(2)
  i1 = com1%coord(3)
  Om1 = com1%coord(4)
  omeg1 = com1%coord(5)
  q2 = com2%coord(1) 
  e2 = com2%coord(2) 
  i2 = com2%coord(3)
  Om2 = com2%coord(4)
  omeg2 = com2%coord(5)

!=================================================================
! RMS of tau1 || tau2
!=================================================================
! pericenter pos. and vel. of 1st orbit  
  r01=q1
  mu1=centerbody_mass(com1%center)
  v01=sqrt(mu1*(1.d0+e1)/q1)
  cosn1=cos(com1%coord(4))
  sinn1=sin(com1%coord(4))
  coso1=cos(com1%coord(5))
  sino1=sin(com1%coord(5))
  cosi1=cos(com1%coord(3))
  sini1=sin(com1%coord(3))
  x01=(/coso1*cosn1-sino1*sinn1*cosi1, coso1*sinn1+sino1*cosn1*cosi1,&
       & sino1*sini1/)*r01
  y01=(/-sino1*cosn1-coso1*sinn1*cosi1, -sino1*sinn1+coso1*cosn1*cosi1,&
       & coso1*sini1/)*v01
  alpha1 = vsize(y01)**2-2.d0*mu1/r01
  r1 = q1*(1.d0+e1)/(1.d0+e1*cos(f1))
  IF(r1.lt.0.d0.and.verb_moid.ge.20)THEN
     WRITE(ierrou,*)'negative r1',r1
     numerr=numerr+1
  ENDIF
  
! pericenter pos. and vel. of 2nd orbit  
  r02=q2
  mu2=centerbody_mass(com2%center)
!  mu2=gms
  v02=sqrt(mu2*(1.d0+e2)/q2)
  cosn2=cos(com2%coord(4))
  sinn2=sin(com2%coord(4))
  coso2=cos(com2%coord(5))
  sino2=sin(com2%coord(5))
  cosi2=cos(com2%coord(3))
  sini2=sin(com2%coord(3))
  x02=(/coso2*cosn2-sino2*sinn2*cosi2, coso2*sinn2+sino2*cosn2*cosi2,&
       & sino2*sini2/)*r02
  y02=(/-sino2*cosn2-coso2*sinn2*cosi2, -sino2*sinn2+coso2*cosn2*cosi2,&
       & coso2*sini2/)*v02
  alpha2 = vsize(y02)**2-2.d0*mu2/r02
  r2 = q2*(1.d0+e2)/(1.d0+e2*cos(f2))
  IF(r2.lt.0.d0.and.verb_moid.ge.20)THEN
     WRITE(ierrou,*)'negative r2',r2
     numerr=numerr+1
  ENDIF

! starting guess for the Earth
!  IF(e1.lt.1.d0) THEN
!     u1=2*ATAN(TAN(f1/2.d0)*SQRT((1.d0-e1)/(1.d0+e1)))
!     write(*,*)'f1=',f1*degrad,'u1=',u1*degrad
!     psi1_ini = u1/sqrt(-alpha1)
!     write(*,*)'psi1_ini=',psi1_ini,'u1=',u1,'alpha1=',alpha1
!  ELSEIF(e1.eq.1.d0) THEN
!     u1 = tan(f1/2.d0)
!     psi1_ini = u1*sqrt(2.d0*q1/mu1)
!  ELSEIF(e1.gt.1.d0) THEN
!     u1x=TAN(f1/2.d0)*SQRT((e1-1.d0)/(1.d0+e1))
!     u1=log((1.d0+u1x)/(1.d0-u1x))
!     psi1_ini = u1/sqrt(alpha1)
!  ENDIF
! starting guess for the asteroid
!  IF(e2.lt.1.d0) THEN
!     u2=2*ATAN(TAN(f2/2.d0)*SQRT((1.d0-e2)/(1.d0+e2)))
!     psi2_ini = u2/sqrt(-alpha2)
!  ELSEIF(e2.eq.1.d0) THEN
!     u2 = tan(f2/2.d0)
!     psi2_ini = u2*sqrt(2.d0*q2/mu2)
!  ELSEIF(e2.gt.1.d0) THEN
!     u2x=TAN(f2/2.d0)*SQRT((e2-1.d0)/(1.d0+e2))
!     u2=log((1.d0+u2x)/(1.d0-u2x))
!     psi2_ini = u2/sqrt(alpha2)
!  ENDIF
! time of passage at perihelion
!  CALL pdtime(q1,r1,mu1,alpha1,psi1_ini,psi1,dt1,1.d-10)
!  ptime1 = -dt1+com1%t
!  CALL pdtime(q2,r2,mu2,alpha2,psi2_ini,psi2,dt2,1.d-10)
!  ptime2 = -dt2+com2%t

  com1min=com1
  com2min=com2
! substitute last element time of perihelion passage
  com1min%coord(6)=f1
  com2min%coord(6)=f2
  CALL coo_cha(com1min,'CAR',car1min,fail_flag)
  CALL coo_cha(com2min,'CAR',car2min,fail_flag)       

! tangent vectors at minimum points
  tau1(1) = (cos(Om1)*cos(omeg1)-sin(Om1)*cos(i1)*sin(omeg1))*q1*(1.D&
     &0+e1)/(1.D0+e1*cos(f1))**2.D0*cos(f1)*e1*sin(f1)-(cos(Om1)*cos(ome&
     &g1)-sin(Om1)*cos(i1)*sin(omeg1))*q1*(1.D0+e1)/(1.D0+e1*cos(f1))*si&
     &n(f1)+(-cos(Om1)*sin(omeg1)-sin(Om1)*cos(i1)*cos(omeg1))*q1*(1.D0+&
     &e1)/(1.D0+e1*cos(f1))**2.D0*sin(f1)**2.D0*e1+(-cos(Om1)*sin(omeg1)&
     &-sin(Om1)*cos(i1)*cos(omeg1))*q1*(1.D0+e1)/(1.D0+e1*cos(f1))*cos(f1)
                                                                        
  tau1(2) = (sin(Om1)*cos(omeg1)+cos(Om1)*cos(i1)*sin(omeg1))*q1*(1.D&
     &0+e1)/(1.D0+e1*cos(f1))**2.D0*cos(f1)*e1*sin(f1)-(sin(Om1)*cos(ome&
     &g1)+cos(Om1)*cos(i1)*sin(omeg1))*q1*(1.D0+e1)/(1.D0+e1*cos(f1))*si&
     &n(f1)+(-sin(Om1)*sin(omeg1)+cos(Om1)*cos(i1)*cos(omeg1))*q1*(1.D0+&
     &e1)/(1.D0+e1*cos(f1))**2.D0*sin(f1)**2.D0*e1+(-sin(Om1)*sin(omeg1)&
     &+cos(Om1)*cos(i1)*cos(omeg1))*q1*(1.D0+e1)/(1.D0+e1*cos(f1))*cos(f1)
                                                                        
  tau1(3) = sin(i1)*sin(omeg1)*q1*(1.D0+e1)/(1.D0+e1*cos(f1))**2.D0*c&
     &os(f1)*e1*sin(f1)-sin(i1)*sin(omeg1)*q1*(1.D0+e1)/(1.D0+e1*cos(f1)&
     &)*sin(f1)+sin(i1)*cos(omeg1)*q1*(1.D0+e1)/(1.D0+e1*cos(f1))**2.D0*&
     &sin(f1)**2.D0*e1+sin(i1)*cos(omeg1)*q1*(1.D0+e1)/(1.D0+e1*cos(f1))&
     &*cos(f1)                                                          
                                                                        
  tau2(1) = (cos(Om2)*cos(omeg2)-sin(Om2)*cos(i2)*sin(omeg2))*q2*(1.D&
     &0+e2)/(1.D0+e2*cos(f2))**2.D0*cos(f2)*e2*sin(f2)-(cos(Om2)*cos(ome&
     &g2)-sin(Om2)*cos(i2)*sin(omeg2))*q2*(1.D0+e2)/(1.D0+e2*cos(f2))*si&
     &n(f2)+(-cos(Om2)*sin(omeg2)-sin(Om2)*cos(i2)*cos(omeg2))*q2*(1.D0+&
     &e2)/(1.D0+e2*cos(f2))**2.D0*sin(f2)**2.D0*e2+(-cos(Om2)*sin(omeg2)&
     &-sin(Om2)*cos(i2)*cos(omeg2))*q2*(1.D0+e2)/(1.D0+e2*cos(f2))*cos(f2)
                                                                        
  tau2(2) = (sin(Om2)*cos(omeg2)+cos(Om2)*cos(i2)*sin(omeg2))*q2*(1.D&
     &0+e2)/(1.D0+e2*cos(f2))**2.D0*cos(f2)*e2*sin(f2)-(sin(Om2)*cos(ome&
     &g2)+cos(Om2)*cos(i2)*sin(omeg2))*q2*(1.D0+e2)/(1.D0+e2*cos(f2))*si&
     &n(f2)+(-sin(Om2)*sin(omeg2)+cos(Om2)*cos(i2)*cos(omeg2))*q2*(1.D0+&
     &e2)/(1.D0+e2*cos(f2))**2.D0*sin(f2)**2.D0*e2+(-sin(Om2)*sin(omeg2)&
     &+cos(Om2)*cos(i2)*cos(omeg2))*q2*(1.D0+e2)/(1.D0+e2*cos(f2))*cos(f2)
                                                                        
  tau2(3) = sin(i2)*sin(omeg2)*q2*(1.D0+e2)/(1.D0+e2*cos(f2))**2.D0*c&
     &os(f2)*e2*sin(f2)-sin(i2)*sin(omeg2)*q2*(1.D0+e2)/(1.D0+e2*cos(f2)&
     &)*sin(f2)+sin(i2)*cos(omeg2)*q2*(1.D0+e2)/(1.D0+e2*cos(f2))**2.D0*&
     &sin(f2)**2.D0*e2+sin(i2)*cos(omeg2)*q2*(1.D0+e2)/(1.D0+e2*cos(f2))&
     &*cos(f2)                                                          
                                                                        
! compute joining unit vector tau3hat
  CALL prvec(tau1,tau2,tau3)
  tau3hat=tau3/vsize(tau3)  

! -------------------------------------
  ltau1=vsize(tau1)
  tau1vers(1:3)=tau1(1:3)/ltau1

  ltau2=vsize(tau2)
  tau2vers(1:3)=tau2(1:3)/ltau2
!  write(*,*)'tau2vers',tau2vers(1:3)

  CALL prvec(tau1vers,tau2vers,tau3vec)  ! tau3vec has the same direction
                                         ! as tau3 but it is not the same
  sint1t2=vsize(tau3vec)

END SUBROUTINE tau1_tau2

! =================================================================
  SUBROUTINE derdmintest(com1,com2,nstat,nummin,f1min,f2min,&
       & sinmutI,D2min,detH,sint1t2)
    TYPE(orbit_elem),INTENT(IN) :: com1,com2 ! COM elements
                                         ! 1: planet  2: asteroid
    INTEGER,INTENT(IN) :: nstat,nummin
    DOUBLE PRECISION,DIMENSION(nminx),INTENT(IN) :: f1min,f2min
    DOUBLE PRECISION,INTENT(IN) :: sinmutI
    DOUBLE PRECISION,DIMENSION(nminx),INTENT(IN) :: D2min 
    DOUBLE PRECISION,DIMENSION(nminx),INTENT(IN) :: detH
    DOUBLE PRECISION,DIMENSION(nminx),INTENT(IN) :: sint1t2
! ----- end interface ------
    TYPE(orbit_elem) :: com2var ! equinoctial elements
    DOUBLE PRECISION :: a2,h2,k2,p2,q2
    DOUBLE PRECISION :: a1,e1
    DOUBLE PRECISION :: pd2var,e2var,i2var,bigomega2var,smallomega2var,a2var
! -----test derivatives of detH--------------
    DOUBLE PRECISION :: eps_q2,eps_e2,eps_I2,eps_Om2,eps_omeg2
    TYPE(orbit_elem) :: com1var1,com2var1 ! auxiliary COM elements
    TYPE(orbit_elem) :: ekplvar1,elkepvar1 ! auxiliary KEP elements
    DOUBLE PRECISION,DIMENSION(nminx) :: detHvar,trHvar
    DOUBLE PRECISION,DIMENSION(2,2) :: Hessvar
    DOUBLE PRECISION :: sinmutIvar
! number of relative minima found 
    INTEGER :: numminvar
! number of statpts,minima and maxima
    INTEGER :: nstatvar,nummaxvar
! error flags for crit_pts
    INTEGER :: answervar(poldeg)
    LOGICAL :: morsevar
    LOGICAL :: weiervar
    LOGICAL :: warnflagvar(3)
    LOGICAL :: sflagvar(6)
    LOGICAL :: hzflagvar
    LOGICAL :: hwflagvar
    LOGICAL :: multflvar
    LOGICAL :: hevalflagvar
! -------------------------------------------------------------------
    DOUBLE PRECISION :: ff1var,ff2var,chf1var,chf2var 
    DOUBLE PRECISION, DIMENSION(poldeg):: D2var
    DOUBLE PRECISION, DIMENSION(nminx) :: D2minv
    DOUBLE PRECISION :: DD2var
    DOUBLE PRECISION, DIMENSION(poldeg) :: f1var,f2var
    DOUBLE PRECISION, DIMENSION(nminx) :: f1minvar,f2minvar
!ordered 3-uples,for sorting 
    INTEGER, DIMENSION(poldeg) :: srtnumvar,srtnumvar1 
    TYPE(orbit_elem) :: car1min,car2min
    DOUBLE PRECISION,DIMENSION(3) :: tau1,tau2,tau3hat
    DOUBLE PRECISION, DIMENSION(nminx) :: sint1t2var
    INTEGER :: kminvar ! loop index 
    INTEGER :: fail_flag
!    DOUBLE PRECISION :: pd2
    INTEGER :: i,j,k
! ====================================================

! to test the derivatives
    OPEN(10,file='moidrms.opt',status='old')
    READ(10,*)eps_q2,eps_e2,eps_I2,eps_Om2,eps_omeg2
    CLOSE(10)

    com2var=com2
    com2var%coord(1)=com2%coord(1)+eps_q2
    com2var%coord(2)=com2%coord(2)+eps_e2
    com2var%coord(3)=com2%coord(3)+eps_I2
    com2var%coord(4)=com2%coord(4)+eps_Om2
    com2var%coord(5)=com2%coord(5)+eps_omeg2
! -----------------------------------------------------------------------
! test derivatives of mutual inclination
    CALL sin_mutI(com1%coord(1:5),com2var%coord(1:5),sinmutIvar)
!    write(*,*)'sinmutIvar,sinmutI',sinmutIvar,sinmutI
    write(ierrou,*)'derivative of sinmutI computed by the incremental ratio:',&
         & (sinmutIvar-sinmutI)/(eps_q2+eps_e2+eps_I2+eps_Om2+eps_omeg2)
    numerr=numerr+1
  !   initialization
    srtnumvar(1:poldeg)=0.d0 
  !   flags initialization                                              
    weiervar = .true.
    morsevar = .true.
    sflagvar(1) = .true. 
    sflagvar(2) = .true. 
    sflagvar(3) = .true.
    sflagvar(4) = .true.
    sflagvar(5) = .true.
    sflagvar(6) = .true.
    warnflagvar(1) = .true.
    warnflagvar(2) = .true.
    warnflagvar(3) = .true.
    hzflagvar = .true.
    hwflagvar = .true.
    multflvar = .true.
    hevalflagvar = .true.
! -----------------------------------------------------------------------
    CALL crit_pts(com1%coord(1:5),com2var%coord(1:5),f1var,f2var,nstatvar, &
         & numminvar,nummaxvar,answervar,warnflagvar,sflagvar, &
         & morsevar,weiervar,hzflagvar,hwflagvar,multflvar,hevalflagvar)
! -----------------------------------------------------------------------  
    IF(.not.hevalflagvar) THEN
       write(ierrou,*)'derdmintest: hessian evaluation failed'
       numerr=numerr+1
    ENDIF
    DO j = 1,nstatvar
       f1var(j)=f1var(j)*radeg 
       f2var(j)=f2var(j)*radeg 
       CALL d2eval(f1var(j),f2var(j),DD2var)
       D2var(j)=DD2var 
    ENDDO
    CALL heapsort(D2var,nstatvar,srtnumvar)
    kminvar = 0 !counting minimum points
    DO k = 1,nstatvar 
       IF(answervar(srtnumvar(k)).eq.-1) THEN
          kminvar=kminvar+1
          f1minvar(kminvar) = f1var(srtnumvar(k))
          f2minvar(kminvar) = f2var(srtnumvar(k))
          D2minv(kminvar) = D2var(srtnumvar(k))
       ENDIF
    ENDDO
    IF((nstat.ne.nstatvar).or.(nummin.ne.numminvar)) THEN
       WRITE(*,*)'derdmintest: ERROR! nstat',nstat,'nstatvar',nstatvar
       WRITE(*,*)'                nummin',nummin,'numminvar',numminvar
       STOP
    ENDIF
!    write(*,*)'nummin,D2min',nummin,D2min(1:nummin)
!    write(*,*)'numminvar,D2minv',numminvar,D2minv(1:nummin)
    DO j = 1,numminvar 
!       write(*,*)'variated distance, distance:',sqrt(D2minv(j)),sqrt(D2min(j))
       write(ierrou,*)'derivative of dmin computed by the incremental  &
            & ratio:',(sqrt(D2minv(j))-sqrt(D2min(j)))/ &
            & (eps_q2+eps_e2+eps_I2+eps_Om2+eps_omeg2)
       numerr=numerr+1
    ENDDO
! loop on minimum points
    DO i=1,numminvar
! test derivatives of detH
       CALL det_H(com1,com2var,f1minvar(i),f2minvar(i), &
            & Hessvar,detHvar(i),trHvar(i)) 

!       write(*,*)'detHvar',detHvar(i)
       write(ierrou,*)'derivative of detH computed by the &
            & incremental ratio:', &
            & (detHvar(i)-detH(i))/(eps_q2+eps_e2+eps_I2+eps_Om2+eps_omeg2)
       numerr=numerr+1
! test derivatives of st1t2
       CALL tau1_tau2(com1,com2var,f1minvar(i),f2minvar(i),car1min,car2min,&
            & tau1,tau2,tau3hat,sint1t2var(i))

        write(ierrou,*)'derivative of sint1t2 computed by the  &
        & incremental ratio:', &
             & (sint1t2var(i)-sint1t2(i))/ &
             & (eps_q2+eps_e2+eps_I2+eps_Om2+eps_omeg2)
        numerr=numerr+1
     ENDDO

   END SUBROUTINE derdmintest

! ==============================================================
! COMPUTE the RMS of the mutual inclination
! written by G.F. Gronchi and G. Tommei 28/02/2006
SUBROUTINE sin_mutI(com1,com2,sinmutI)
  DOUBLE PRECISION,DIMENSION(5),INTENT(IN):: com1,com2
! cometary elements ---  1: planet   2: asteroid
  DOUBLE PRECISION,INTENT(OUT):: sinmutI
! ------ end interface ------
  DOUBLE PRECISION :: i1,Om1 ! Planet
  DOUBLE PRECISION :: i2,Om2  ! Asteroid/Comet
  DOUBLE PRECISION :: cosmutI
!===============================================================
  i1=com1(3)
  Om1=com1(4)
  i2=com2(3)
  Om2=com2(4)
  cosmutI= sin(Om1)*sin(i1)*sin(Om2)*sin(i2)+cos(Om1)* &
       & sin(i1)*cos(Om2)*sin(i2)+cos(i1)*cos(i2)
  sinmutI=sqrt(1.d0-cosmutI**2)
  END SUBROUTINE sin_mutI

! ==================================================
! COMPUTE the DETERMINANT of HESSIAN MATRIX
! written by Giovanni Federico Gronchi & Giacomo Tommei, 7/2/2007
! ==================================================
SUBROUTINE det_H(com1,com2,f1,f2,Hess,detH,trH)
  IMPLICIT NONE
!===================================================
  TYPE(orbit_elem), INTENT(IN) :: com1,com2 ! equinoctial elements
                                          ! 1: planet 
                                          ! 2: asteroid
  DOUBLE PRECISION, INTENT(IN) :: f1,f2   ! anomalies at the critical point
                                          ! in [0,2*pi]
  DOUBLE PRECISION, DIMENSION(2,2), INTENT(OUT) :: Hess 
  DOUBLE PRECISION, INTENT(OUT) :: detH,trH
!===================================================
  INTEGER :: fail_flag
  DOUBLE PRECISION :: q1,e1,i1,Om1,omeg1  ! Planet
  DOUBLE PRECISION :: q2,e2,i2,Om2,omeg2  ! Asteroid/Comet
!  DOUBLE PRECISION :: deltaOm,KK,LL,MM,NN
!  DOUBLE PRECISION :: PPx,PPy,PPz,QQx,QQy,QQz,px,py,pz,qx,qy,qz
!---------------------------------------------------
  DOUBLE PRECISION :: gradd2f1,gradd2f2
!====================================================

    q1 = com1%coord(1) 
    e1 = com1%coord(2)
    i1 = com1%coord(3)
    Om1 = com1%coord(4)
    omeg1 = com1%coord(5)
    q2 = com2%coord(1) 
    e2 = com2%coord(2) 
    i2 = com2%coord(3)
    Om2 = com2%coord(4)
    omeg2 = com2%coord(5)
!    deltaOm = Om2-Om1

!    PPx = cos(omeg1) 
!    PPy = sin(omeg1)*cos(i1)
!    PPz = sin(omeg1)*sin(i1)

!    QQx = -sin(omeg1) 
!    QQy = cos(omeg1)*cos(i1)
!    QQz = cos(omeg1)*sin(i1)

!    px = cos(omeg2)*cos(deltaOm) - &
!         & sin(omeg2)*cos(i2)*sin(deltaOm)
!    py = cos(omeg2)*sin(deltaOm) + &
!         & sin(omeg2)*cos(i2)*cos(deltaOm)
!    pz = sin(omeg2)*sin(i2)

!    qx = -sin(omeg2)*cos(deltaOm) - &
!         & cos(omeg2)*cos(i2)*sin(deltaOm)
!    qy = -sin(omeg2)*sin(deltaOm) + &
!         & cos(omeg2)*cos(i2)*cos(deltaOm)
!    qz = cos(omeg2)*sin(i2)

!    KK = PPx*px + PPy*py + PPz*pz
!    LL = QQx*px + QQy*py + QQz*pz
!    MM = PPx*qx + PPy*qy + PPz*qz
!    NN = QQx*qx + QQy*qy + QQz*qz

!====================================================================
! GRADIENT of d^2
!====================================================================
    gradd2f1                                                          &
     &= 2.D0/(1.D0+e1*cos(f1))*(e1*q1**2*(1.D0+e1)**2.D0/(1.D0+e1*cos   &
     &(f1))**2.D0*sin(f1)+q1*(1.D0+e1)/(1.D0+e1*cos(f1))*sin(f1)*(KK*q2*&
     &(1.D0+e2)/(1.D0+e2*cos(f2))*cos(f2)+MM*q2*(1.D0+e2)/(1.D0+e2*cos(f&
     &2))*sin(f2))-(e1*q1*(1.D0+e1)/(1.D0+e1*cos(f1))+q1*(1.D0+e1)/(1.D0&
     &+e1*cos(f1))*cos(f1))*(LL*q2*(1.D0+e2)/(1.D0+e2*cos(f2))*cos(f2)+N&
     &N*q2*(1.D0+e2)/(1.D0+e2*cos(f2))*sin(f2)))                        
                                                                        
    gradd2f2                                                          &
     &= 2.D0/(1.D0+e2*cos(f2))*(e2*q2**2*(1.D0+e2)**2.D0/(1.D0+e2*cos   &
     &(f2))**2.D0*sin(f2)+q2*(1.D0+e2)/(1.D0+e2*cos(f2))*sin(f2)*(KK*q1*&
     &(1.D0+e1)/(1.D0+e1*cos(f1))*cos(f1)+LL*q1*(1.D0+e1)/(1.D0+e1*cos(f&
     &1))*sin(f1))-(e2*q2*(1.D0+e2)/(1.D0+e2*cos(f2))+q2*(1.D0+e2)/(1.D0&
     &+e2*cos(f2))*cos(f2))*(MM*q1*(1.D0+e1)/(1.D0+e1*cos(f1))*cos(f1)+N&
     &N*q1*(1.D0+e1)/(1.D0+e1*cos(f1))*sin(f1)))                        

! force to zero (critical point)
!      gradd2f1=0.D0
!      WRITE(*,*)'gradd2f1:',gradd2f1
!      gradd2f2=0.D0
!      WRITE(*,*)'gradd2f2:',gradd2f2

!====================================================================
! Hessian Matrix H of d^2
!====================================================================
    Hess(1,1)&
     &= e1/(1.D0+e1*cos(f1))*sin(f1)*gradd2f1+2.D0/(1.D0+e1*cos(f1))*&
     &(e1**2*q1**2*(1.D0+e1)**2.D0/(1.D0+e1*cos(f1))**3.D0*sin(f1)**2.D0&
     &+(e1*q1*(1.D0+e1)/(1.D0+e1*cos(f1))**2.D0*sin(f1)**2.D0+q1*(1.D0+e&
     &1)/(1.D0+e1*cos(f1))*cos(f1))*(e1*q1*(1.D0+e1)/(1.D0+e1*cos(f1))+K&
     &K*q2*(1.D0+e2)/(1.D0+e2*cos(f2))*cos(f2)+MM*q2*(1.D0+e2)/(1.D0+e2*&
     &cos(f2))*sin(f2))-(e1**2*q1*(1.D0+e1)/(1.D0+e1*cos(f1))**2.D0*sin(&
     &f1)+e1*q1*(1.D0+e1)/(1.D0+e1*cos(f1))**2.D0*cos(f1)*sin(f1)-q1*(1.&
     &D0+e1)/(1.D0+e1*cos(f1))*sin(f1))*(LL*q2*(1.D0+e2)/(1.D0+e2*cos(f2&
     &))*cos(f2)+NN*q2*(1.D0+e2)/(1.D0+e2*cos(f2))*sin(f2)))            
    Hess(2,2)&
     &= e2/(1.D0+e2*cos(f2))*sin(f2)*gradd2f2+2.D0/(1.D0+e2*cos(f2))*&
     &(e2**2*q2**2*(1.D0+e2)**2.D0/(1.D0+e2*cos(f2))**3.D0*sin(f2)**2.D0&
     &+(e2*q2*(1.D0+e2)/(1.D0+e2*cos(f2))**2.D0*sin(f2)**2.D0+q2*(1.D0+e&
     &2)/(1.D0+e2*cos(f2))*cos(f2))*(e2*q2*(1.D0+e2)/(1.D0+e2*cos(f2))+K&
     &K*q1*(1.D0+e1)/(1.D0+e1*cos(f1))*cos(f1)+LL*q1*(1.D0+e1)/(1.D0+e1*&
     &cos(f1))*sin(f1))-(e2**2*q2*(1.D0+e2)/(1.D0+e2*cos(f2))**2.D0*sin(&
     &f2)+e2*q2*(1.D0+e2)/(1.D0+e2*cos(f2))**2.D0*cos(f2)*sin(f2)-q2*(1.&
     &D0+e2)/(1.D0+e2*cos(f2))*sin(f2))*(MM*q1*(1.D0+e1)/(1.D0+e1*cos(f1&
     &))*cos(f1)+NN*q1*(1.D0+e1)/(1.D0+e1*cos(f1))*sin(f1)))            
    Hess(1,2)&
     &= 2.D0/(1.D0+e1*cos(f1))*(q1*(1.D0+e1)/(1.D0+e1*cos(f1))*sin(f1&
     &)*(KK*(e2*q2*(1.D0+e2)/(1.D0+e2*cos(f2))**2.D0*cos(f2)*sin(f2)-q2*&
     &(1.D0+e2)/(1.D0+e2*cos(f2))*sin(f2))+MM*(e2*q2*(1.D0+e2)/(1.D0+e2*&
     &cos(f2))**2.D0*sin(f2)**2.D0+q2*(1.D0+e2)/(1.D0+e2*cos(f2))*cos(f2&
     &)))-(e1*q1*(1.D0+e1)/(1.D0+e1*cos(f1))+q1*(1.D0+e1)/(1.D0+e1*cos(f&
     &1))*cos(f1))*(LL*(e2*q2*(1.D0+e2)/(1.D0+e2*cos(f2))**2.D0*cos(f2)*&
     &sin(f2)-q2*(1.D0+e2)/(1.D0+e2*cos(f2))*sin(f2))+NN*(e2*q2*(1.D0+e2&
     &)/(1.D0+e2*cos(f2))**2.D0*sin(f2)**2.D0+q2*(1.D0+e2)/(1.D0+e2*cos(&
     &f2))*cos(f2))))                                                   
    Hess(2,1)=Hess(1,2)
! Determinant of Hess
    detH=Hess(1,1)*Hess(2,2)-Hess(1,2)**2
! Trace of Hess
    trH=Hess(1,1)+Hess(2,2)
! Inverse matrix of H
!    IF(detH.ne.0)THEN
!       invH(1,1)=Hess(2,2)/detH
!       invH(2,2)=Hess(1,1)/detH
!       invH(1,2)=-Hess(1,2)/detH
!       invH(2,1)=-Hess(2,1)/detH
!    ELSE
!       WRITE(*,*)'det_H: error! determinant =',detH
!    ENDIF
END SUBROUTINE det_H

END MODULE critical_points

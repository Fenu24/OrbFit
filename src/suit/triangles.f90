MODULE triangles

PRIVATE

! public data
INTEGER, PARAMETER :: npox=200, ntrix=2*npox ! ??????
PUBLIC npox, ntrix

INTEGER nfunc ! choice of metric function for triangulation
PUBLIC nfunc ! 1= adaptive 2= log10(r) 3=r

! public routines
PUBLIC triangulate,outtriang,shstar_lim,triang_outer,sample_bound_ne1!, sorting

CONTAINS

SUBROUTINE triang_outer(npo,pmult,nroots,roots,c,refs,xo,vo, &
     & E_bound,att0,nrootsd,rootsd,rrdot,triang,npt,ntri,name0)
  USE attributable
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: npo,pmult !controls on no of pts, densification
  INTEGER,INTENT(IN) :: nroots ! number of positive roots; 
  DOUBLE PRECISION, DIMENSION(3), INTENT(IN) :: roots
  DOUBLE PRECISION, INTENT(IN) ::  c(0:5)! polynomial coefficients
! matrix of the r eps theta ref. system, topocentric equatorial
  DOUBLE PRECISION, INTENT(IN) :: refs(3,3)
! observer pos/vel, asteroid pos/vel, equatorial heliocentric
  DOUBLE PRECISION, INTENT(IN) :: xo(3),vo(3)
  DOUBLE PRECISION, INTENT(IN) :: E_bound ! bound for the energy
  TYPE(attrib), INTENT(IN) :: att0
  INTEGER, INTENT(IN) :: nrootsd ! number of roots of the poly derivative 
  DOUBLE PRECISION,DIMENSION(8),INTENT(IN) :: rootsd
! node coordinates (r,rdot)
  DOUBLE PRECISION,DIMENSION(2,npox),INTENT(OUT) :: rrdot 
  INTEGER,DIMENSION(ntrix,3),INTENT(OUT) :: triang! triangle indexes
  INTEGER,INTENT(INOUT) :: npt,ntri ! number of triangles, of points  
  CHARACTER*(*), INTENT(IN) :: name0 ! name 
! ============== end interface ====================
  DOUBLE PRECISION :: apm ! apparent magnitude
! ----- for sample_bound --------------------------
  DOUBLE PRECISION :: r1,r2,r2out
  DOUBLE PRECISION :: eta2 ! squared proper motion
! output of sample_bound: npoout,frv,rdot,rdw
  INTEGER :: nposb
  DOUBLE PRECISION,DIMENSION(npox) :: frv,rdot
  DOUBLE PRECISION :: rdw
  DOUBLE PRECISION,DIMENSION(npox) :: fr,rd
! ----- for triangola -----------------------------
! output of triangola: nptout,triangout,ntriout
  INTEGER :: nptouter,ntriouter ! no of points and of triangles
  INTEGER,DIMENSION(ntrix,4) :: triangouter! triangle indexes
!
  DOUBLE PRECISION, DIMENSION(npox) :: rr
  INTEGER :: i ! loop index
! =================================================
  rrdot=0.d0
  triang=0
  apm = att0%apm
  eta2 = att0%eta**2
  r1=roots(2)
  r2=roots(3)
  r2out = r2

  CALL sample_bound(E_bound,eta2,r1,r2out,2*npo,pmult,nposb,c,&
       &frv,rdot,npox,rdw)
  IF(nposb.lt.1)THEN
     npt=0
     ntri=0
     RETURN
  ENDIF
  fr(1:nposb)=frv(1:nposb)
  rd(1:nposb)=rdot(1:nposb)*rdw ! scaling before calling triangola
  
  CALL triangola(npox,ntrix,fr,rd,nposb,nptouter,triangouter,ntriouter)

  rd(1:nptouter)=rd(1:nptouter)*(1.d0/rdw)
  DO i = 1,nptouter
     rr(i) = f_inv(fr(i),r2out)
  ENDDO
  
  rrdot(1,1:nptouter) = rr(1:nptouter) 
  rrdot(2,1:nptouter) = rd(1:nptouter)  
  triang(1:ntriouter,3) = triangouter(1:ntriouter,3)
  npt = nptouter
  ntri = ntriouter 

END SUBROUTINE triang_outer

SUBROUTINE triangulate(att,npo,pmult,h_max,roots3,nroots,rshoot &
&              ,ncomp,rrdot,triango,npt,ntri,name0,iunmat,iunpol)
  USE attributable
  USE output_control
  IMPLICIT NONE
  TYPE(attrib), INTENT(IN) :: att ! attributable
! options for triangulation
  INTEGER, INTENT(IN) :: npo,pmult ! number of points lower boundary
                                   ! multiplier for uniformization
  DOUBLE PRECISION, INTENT(IN) :: h_max ! absolute magnitude limit
  CHARACTER*(*), INTENT(IN) :: name0 ! name 
  INTEGER, INTENT(IN) :: iunmat,iunpol ! for matlab processing
! TRIANGULATION
  DOUBLE PRECISION, DIMENSION(2,npox), INTENT(OUT) :: rrdot ! nodes (r,rdot)
  INTEGER, INTENT(OUT) ::  triango(ntrix,3)! triangle indexes
  INTEGER, INTENT(OUT) :: nroots ! number of positive roots; 
!                       0 is used as error message
  INTEGER, INTENT(OUT) ::  ncomp ! no. connected components (can be 1,2) 
  DOUBLE PRECISION, INTENT(OUT) ::  roots3(3),rshoot ! range limits
! ============== end interface ====================
  INTEGER ntri,npt ! current number of triangles, of points
!  INTEGER iunm, iunp
!  CHARACTER*(19) name

! admissible region: interface with admis_reg
! ---------------------------------------------
  INTEGER :: nrootsd ! number of positive roots of the derivative 
  DOUBLE PRECISION, DIMENSION(3) :: roots
                                        ! should be not more than 3.
  INTEGER, DIMENSION(6) :: indsrt 
  DOUBLE PRECISION, DIMENSION(8) :: rootsd
                                        ! roots of the derivative
  INTEGER, DIMENSION(8) :: inddsrt 
  DOUBLE PRECISION ::  c(0:5)! polynomial coefficients

  INTEGER :: no,no1,no2,no3
  DOUBLE PRECISION :: r1,r2,r2inn,r2out,rder0
! triangulations: node coordinates (f(r),rdot)
  DOUBLE PRECISION, DIMENSION(npox) :: frv1,frv2,frv3,frv4
  DOUBLE PRECISION, DIMENSION(npox) :: rdot,rdot1,rdot2,rdot3,rdot4
! same for triangulation
  DOUBLE PRECISION, DIMENSION(npox) :: fr,fr1,fr2,fr3,fr4
  DOUBLE PRECISION, DIMENSION(npox) :: rd,rd1,rd2,rd3,rd4
  DOUBLE PRECISION, DIMENSION(npox) :: rr,rr1,rr2,rr3,rr4
!triangle indexes
  INTEGER :: triang(ntrix,4),triang1(ntrix,4),triang2(ntrix,4), &
     & triang3(ntrix,4),triang4(ntrix,4)
! partial number of points
  INTEGER :: npt1,npt2,npt3,npt4
! partial number of triangles
  INTEGER :: ntri1,ntri2,ntri3,ntri4 
! 1/weight for rdot scaling
  DOUBLE PRECISION :: rdw
! proper motion squared
  DOUBLE PRECISION :: eta2
! matrix of the r eps theta ref. system, topocentric equatorial
  DOUBLE PRECISION :: refs(3,3)
! observer pos/vel, equatorial heliocentric
  DOUBLE PRECISION :: xo(3),vo(3)
! computation of magnitude
  DOUBLE PRECISION :: phase,rsun,apm,apmf
! -----------------------------------------
! bound for the energy
  DOUBLE PRECISION :: a_orb,E_bound
! loop indexes
  INTEGER i,ind
  ! functions
  DOUBLE PRECISION :: vsize,prscal,appmag!,f_inv,f_func,shstar_lim

! initialization
  roots3(1:3) = 0.d0

! &&&&&&&&&&&&&&&&&&&&&&&&&
! find admissible region
! &&&&&&&&&&&&&&&&&&&&&&&&&
  a_orb = 100.d0 ! comet boundary
  eta2=(att%eta)**2

  CALL admis_reg(name0,iunmat,att%tdtobs,att%angles,att%obscod,nroots,&
           & roots,c,refs,xo,vo,a_orb,E_bound,nrootsd,rootsd)
! roots and rootsd have been sorted in admis_reg

  IF(nroots.gt.0.and.nroots.le.3)THEN 
     roots3(1:nroots)=roots(1:nroots) ! output of triangulate
     ncomp=1
     IF(nroots.eq.3)ncomp=2
  ELSE
     WRITE(*,*)' wrong case: ',name0,' nroots ',nroots
  ENDIF
  apm=att%apm
  IF(apm.le.0.d0.or.apm.gt.30.d0)THEN
     apmf=12
  ELSE
     apmf=apm
  ENDIF

  eta2=(att%eta)**2 ! seems redundant

! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! TWO connected components with non-empty interior
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  IF(nroots.eq.3)THEN
!     write(*,*)'two c.c.'

! &&&&&&&&&&&&&&
! OUTER REGION
! &&&&&&&&&&&&&&
     r1=roots(2)
     r2=roots(3)
     r2out = r2
     CALL sample_bound(E_bound,eta2,r1,r2out,2*npo,pmult,no3,c,&
              &frv3,rdot3,npox,rdw)
! if second component too small (no3=0) then skip to one c.c.
     IF(no3.le.1) THEN
        nroots=1
        WRITE(ierrou,*)'  too small second component!:'
        WRITE(ierrou,*)name0,nroots,no3
        numerr=numerr+1
! WARNING: compiler message: A jump into a block from outside the block 
! check if it may happen!
        GOTO 16        
     ENDIF

     fr3(1:no3)=frv3(1:no3)
     rd3(1:no3)=rdot3(1:no3)*rdw       
     CALL triangola(npox,ntrix,fr3,rd3,no3,npt3,triang3,ntri3)
     rd3(1:npt3)=rd3(1:npt3)*(1.d0/rdw)
     DO i = 1,npt3
        rr3(i) = f_inv(fr3(i),r2out)
     ENDDO
! &&&&&&&&&&&&&&&&&&&
! NEAR EARTH REGION
! &&&&&&&&&&&&&&&&&&&
! approximate phase computed as pi-elongation (valid for small r)
     phase=acos(prscal(refs(1:3,1),xo)/vsize(xo))
     rshoot=shstar_lim(apmf,h_max,phase)
     r1=rshoot
     IF(r1.lt.3.d-7)THEN
        r1=3.d-7
     ENDIF
     r2=roots(1)
     r2inn = r2
     IF(r1.gt.r2)THEN
        WRITE(ierrou,199)name0, r1,r2
199     FORMAT(a9,' triangulate fails r1>r2 ',2f10.6)
        numerr=numerr+1
! *********** temporary ******************
        nroots=0
! copy rr3,rd3,npt3,ntri3  into output arrays 
        RETURN
     ENDIF
! search for a zero of the derivative between r1,r2 
     rder0 = -1.d0
     ind = 0
     IF(nrootsd.ge.1) THEN
        IF(nrootsd.eq.1) THEN
               ! do nothing
        ELSEIF(nrootsd.eq.2) THEN
               ! do nothing
        ELSEIF(nrootsd.eq.3) THEN
           DO i=1,nrootsd
              IF((rootsd(i).gt.r1).and.(rootsd(i).lt.r2)) THEN
                 ind=i
              ENDIF
           ENDDO
           IF(ind.gt.2)THEN
              IF(rootsd(ind-2).gt.r1) THEN
                 rder0=rootsd(ind-2)
                 write(*,*)'r1,rder0,r2',r1,rder0,r2
              ENDIF
           ENDIF
           IF(rder0.ge.r2) THEN
              WRITE(*,*)'ERROR!!: rder0 > r2'
           ENDIF
               !            write(*,*)'r1,rder0,r2',r1,rder0,r2
        ELSE
           WRITE(*,*)'not included case: r1,r2 nrootsd=',nrootsd
        ENDIF
     ENDIF
     CALL sample_bound_ne1(E_bound,eta2,r1,rder0,r2inn,npo,pmult,no1,c,&
              & frv1,rdot1,npox,rdw)
     IF(no1.le.1)THEN
        WRITE(ierrou,*)'too small first component!:'
        WRITE(ierrou,*)name0,nroots,no1
        numerr=numerr+1
        nroots=0
        RETURN
     ENDIF
     fr1(1:no1)=frv1(1:no1)
     rd1(1:no1)=rdot1(1:no1)*rdw
     CALL triangola(npox,ntrix,fr1,rd1,no1,npt1,triang1,ntri1)
     rd1(1:npt1)=rd1(1:npt1)*(1.d0/rdw)
     DO i = 1,npt1
        rr1(i) = f_inv(fr1(i),r2inn)
     ENDDO
     npt=npt1+npt3   
     CALL merge_triang(npox,ntrix,rr3,rd3,npt3,triang3,ntri3, &
              &  rr1,rd1,npt1,triang1,ntri1, &
              &  rr,rd,npt,triang,ntri)
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! ONE connected component, near Earth
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  ELSEIF(nroots.eq.1)THEN
16   CONTINUE
! approximate phase computed as pi-elongation (valid for small r)
     phase=acos(prscal(refs(1:3,1),xo)/vsize(xo))
     rshoot=shstar_lim(apmf,h_max,phase)
     r1=rshoot
     IF(r1.lt.3.d-7)THEN
        r1=3.d-7
     ENDIF
     r2=roots(1)
     r2inn = r2
     IF(r1.gt.r2)THEN
        nroots=0
        WRITE(ierrou,199)name0,r1,r2
        numerr=numerr+1
        RETURN
     ENDIF
     r2out=0.d0
     rder0 = -1.d0
     ind = 0
     IF(nrootsd.ge.1) THEN
        DO i=1,nrootsd
           IF((rootsd(i).gt.r1).and.(rootsd(i).lt.r2)) THEN
              ind=i
           ENDIF
        ENDDO
        IF(ind.gt.1)THEN
           IF(rootsd(ind-1).gt.r1) THEN
              rder0=rootsd(ind-1)
           ENDIF
        ENDIF
        IF(rder0.ge.r2) THEN
           WRITE(*,*)'ERROR!!: rder0>r2'
        ENDIF
     ENDIF
     CALL sample_bound_ne1(E_bound,eta2,r1,rder0,r2inn,npo,pmult,no1,c,&
              & frv1,rdot1,npox,rdw)
     IF(no1.le.1)THEN
        WRITE(ierrou,*)'too small unique component!:'
        WRITE(ierrou,*)name0,nroots,no1
        numerr=numerr+1
        nroots=0
        RETURN
     ENDIF
     fr(1:no1)=frv1(1:no1)
     rd(1:no1)=rdot1(1:no1)*rdw
     CALL triangola(npox,ntrix,fr,rd,no1,npt,triang,ntri)
     rd(1:npt)=rd(1:npt)*(1.d0/rdw)
     DO i = 1,npt
        rr(i) = f_inv(fr(i),r2)
     ENDDO   
  ELSEIF(nroots.eq.2) THEN ! some double root
     WRITE(*,*)' this should be rare, nroots= ',nroots,roots(1:nroots)
  ELSE ! karakiri
     WRITE(*,*)' this should not happen, nroots= ',nroots,roots(1:nroots)
     STOP
  ENDIF 
! write on vsaeph.pol
  IF(iunpol.gt.0)THEN
     WRITE(iunpol,177)name0,r1,r2inn,r2out,c
177  FORMAT(a9,9(1x,d14.6)) 
  ENDIF
! copy on output (to be simplified later)
  triango(1:ntri,1:3)=triang(1:ntri,1:3)
  rrdot(1,1:npt)=rr(1:npt)
  rrdot(2,1:npt)=rd(1:npt)
END SUBROUTINE triangulate

! =============================================
!  shstar_lim
! provides the minimum range consistent with the maximum allowed 
! absolute magnitude (shooting star limit), as a function
! of apparent magnitude, max allowed abs. magnitude, and phase
! uses the low range limit
! ============================================
DOUBLE PRECISION FUNCTION shstar_lim(apm,h_max,phase)
  USE fund_const, ONLY: reau
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT(IN) :: apm, h_max, phase ! apparent mag, 
   ! maximum allowed absolute magnitude, phase
  DOUBLE PRECISION appmag, x0_corr
  DOUBLE PRECISION, PARAMETER :: g_mag=0.15d0  
  x0_corr = appmag(0.d0,g_mag,1.d0,1.d0,phase)
  shstar_lim = 10.d0**((apm-h_max-x0_corr)/5.d0)
  shstar_lim=MAX(shstar_lim,reau)
END FUNCTION shstar_lim






! =======================================
!  B O U N D A R Y    S A M P L I N G 
!                of the 
!  A D M I S S I B L E   R E G I O N 
! =======================================
! written by G.F.Gronchi
! February 2003

! ===================================
! SAMPLING of the OUTER REGION
! ===================================
SUBROUTINE sample_bound(E_bound,eta2,r1,r2,npotmp,pmult,npoout,c0,&
           & frv,rdot,npox,rdw)
  USE fund_const
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT(IN) :: E_bound
  DOUBLE PRECISION, INTENT(IN) :: eta2,r1,r2
  DOUBLE PRECISION,INTENT(IN) ::  c0(0:5) ! polynomial coefficients
  INTEGER,INTENT(IN):: npotmp, pmult
  INTEGER,INTENT(OUT):: npoout
  INTEGER,INTENT(IN) :: npox
  DOUBLE PRECISION,INTENT(OUT),DIMENSION(npox) :: frv,rdot
  DOUBLE PRECISION, INTENT(OUT) :: rdw
! END INTERFACE
! ---------------------------------------------------------------------------
  DOUBLE PRECISION,DIMENSION(npox*pmult) :: frv1,rdot1,frv2,rdot2
  DOUBLE PRECISION,DIMENSION(npox*pmult) :: frv3,rdot3
  INTEGER :: npoout3
  INTEGER :: npo
  DOUBLE PRECISION :: c(0:5)
! functions
!  DOUBLE PRECISION :: rdot_comp,f_func
!  DOUBLE PRECISION :: r(npox)
! mass of the Earth+Moon system
  DOUBLE PRECISION, PARAMETER :: meinv=3.2890056140000012D+05
  DOUBLE PRECISION :: me
! radius of the influence sphere of the Earth
  DOUBLE PRECISION :: radSI
! counters
  INTEGER :: npts
! loop indexes
  INTEGER :: h,i
! steps
  DOUBLE PRECISION :: bstep
  DOUBLE PRECISION :: hstep
!
! ***********************************                      
  INTEGER :: extrac(pmult*npox)
! =====================================================================
  IF(r1.ge.r2) THEN
     WRITE(*,199)r1,r2
199  FORMAT(' triangulate fails r1>r2 ',2f10.6)
     STOP
  ENDIF
! for uniformization of the distance on the boundary
  npo=npotmp*pmult
! initialization
  npts=0
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! modified c(4) for strictly negative energy
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  DO i = 0,5
     c(i) = c0(i)
     IF(i.eq.4) THEN
        c(i) = c0(i) - 2.d0*E_bound
     ENDIF
  ENDDO
! &&&&&&&&&&&&
! FIRST STEP
! &&&&&&&&&&&&
  frv1(1) = f_func(r1,r2)
  rdot1(1) = -c(1)/2.d0
  frv2(1) = f_func(r1,r2)
  rdot2(1) = -c(1)/2.d0
  npts = npts+1
  bstep = (f_func(r2,r2)-f_func(r1,r2))/(npo-1)
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! &&&&&&&&&&&&&&
! MAIN DO LOOP 
! &&&&&&&&&&&&&&
  DO i = 1,npox-2  
! &&&&&&&&&&&&&&&&&&&&&&
! COMPUTE LOGRV,RDOT
! &&&&&&&&&&&&&&&&&&&&&&
! sampling in logr
     hstep = bstep
     frv1(i+1) = frv1(i) + hstep
     frv2(i+1) = frv2(i) + hstep
     IF((frv1(i+1).gt.f_func(r2,r2)-bstep/2.d0).or.&
          &(frv2(i+1).gt.f_func(r2,r2)-bstep/2.d0)) THEN 
        GOTO 10
     ENDIF
     rdot1(i+1) = rdot_comp(c,frv1(i+1),-1,r2)   
     rdot2(i+1) = rdot_comp(c,frv2(i+1),1,r2)
     IF((rdot1(i+1).ge.-c(1)/2).or.(rdot2(i+1).le.-c(1)/2)) THEN
        GOTO 10
     ENDIF
     npts = npts+1
     IF((npts.eq.npox-2)) THEN
        write(*,*)'maximum of iterations reached: npts=',npts
     ENDIF   
  ENDDO
10 CONTINUE  
! &&&&&&&&&&&&&&&&&&&       
! ADDING LAST POINT
! &&&&&&&&&&&&&&&&&&& 
  rdot1(npts+1) = -c(1)/2.d0
  rdot2(npts+1) = -c(1)/2.d0
  frv1(npts+1) = f_func(r2,r2)
  frv2(npts+1) = f_func(r2,r2)    
  npts = npts + 1
! &&&&&&&&&&&&&&&&
! scaling factor
! &&&&&&&&&&&&&&&&
 IF ( abs( (MAXVAL(rdot2(1:npts))-MINVAL(rdot1(1:npts))) )/ &
      & MAX(abs(MAXVAL(rdot2(1:npts))),1.d-20) .gt. 1.d-10) THEN
    rdw=(MAXVAL(frv1(1:npts))-MINVAL(frv1(1:npts)))/         &
         &       (MAXVAL(rdot2(1:npts))-MINVAL(rdot1(1:npts))) 
 ELSE
    npoout3 = 0
    write(*,*)'empty second component'
    GOTO 15
 ENDIF
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! joining the upper/lower region 
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  frv3(1:npts)=frv1(1:npts)
  rdot3(1:npts)=rdot1(1:npts)
  DO i =1,npts-1
     frv3(npts+i)=frv2(npts-i)
     rdot3(npts+i)=rdot2(npts-i)
  ENDDO
! &&&&&&&&&&&&&&&&&&&&&&&
! total number of nodes
! &&&&&&&&&&&&&&&&&&&&&&&
  npoout3 = 2*npts-1
! to skip uniformization
  IF(pmult.eq.1) THEN
     GOTO 15
  ENDIF
! &&&&&&&&&&&&&&&&
! uniformization
! &&&&&&&&&&&&&&&&
! --------------------------------------------------------------------
  CALL uniform(frv3(1:npoout3),rdw*rdot3(1:npoout3),npoout3,&
           &npts/pmult+1,extrac,npox*pmult)
  npts = npts/pmult
  DO i = 1,npts
!     frv3(i) = frv3(extrac(i))   
!     rdot3(i) = rdot3(extrac(i))
     frv(i) = frv3(extrac(i))   
     rdot(i) = rdot3(extrac(i))
  ENDDO
  npoout = npts
! --------------------------------------------------------------------
15 CONTINUE
END SUBROUTINE sample_bound


! =======================================
!  B O U N D A R Y    S A M P L I N G 
!                of the 
!  A D M I S S I B L E   R E G I O N 
! =======================================
! written by G.F.Gronchi, October 2003 
! E-mail:gronchi@dm.unipi.it

! ===================================
! SAMPLING of the NEAR EARTH REGION
! ===================================
  SUBROUTINE sample_bound_ne1(E_bound,eta2,r1,rder0,r2,npotmp,pmult,&
       &ntot,c0,frv,rdot,npox,rdw)
    USE fund_const
    USE output_control
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN) :: E_bound
    DOUBLE PRECISION, INTENT(IN) :: eta2,r1,r2,rder0
    DOUBLE PRECISION,INTENT(IN) ::  c0(0:5) ! polynomial coefficients
    INTEGER,INTENT(IN):: npotmp,pmult
    INTEGER,INTENT(OUT):: ntot ! total no of points
    INTEGER,INTENT(IN) :: npox
    DOUBLE PRECISION,INTENT(OUT),DIMENSION(npox) :: frv,rdot
    DOUBLE PRECISION, INTENT(OUT) :: rdw
! END INTERFACE
! ---------------------------------------------------------------------------
    DOUBLE PRECISION :: rdwtmp
    DOUBLE PRECISION,DIMENSION(npox*pmult) :: frv1,rdot1,frv2,rdot2
    DOUBLE PRECISION,DIMENSION(npox*pmult) :: frv3,rdot3,frv4,rdot4
! auxiliary for last call of uniforma
    DOUBLE PRECISION,DIMENSION(npox) :: frvtmp,rdottmp
    DOUBLE PRECISION :: smallpar,dist 
    INTEGER :: ntmp
!
    DOUBLE PRECISION :: c(0:5)
    DOUBLE PRECISION :: delta,bigC,bigS
    DOUBLE PRECISION :: r0,rmin
    DOUBLE PRECISION :: rdotfr1p,rdotfr1m,rdotfrder0
! functions
!    DOUBLE PRECISION :: f_r_rdot0,rdot_comp,f_func,f_inv
    DOUBLE PRECISION :: fr0,fradSI,frmin,fr1,fr2,frder0,frvr2
!    DOUBLE PRECISION :: r(npox)
! mass of the Earth+Moon system
    DOUBLE PRECISION, PARAMETER :: meinv=3.2890056140000012D+05
    DOUBLE PRECISION :: me
! radius of the influence sphere of the Earth
    DOUBLE PRECISION :: radSI
! rdot for the symmetry line of E_sun = -GM/2a
    DOUBLE PRECISION :: rdotsym
! counters
    INTEGER :: npo,npo1,npo2
    INTEGER :: npts1 ! no of points from r1 to r2
    INTEGER :: nout1,nout2 ! no of output of uniforma
    INTEGER :: nout ! no of points from r1 to r2 after uniformization
    INTEGER :: noutsun ! no of points on the E_oplus curve
    INTEGER :: npoearth
    INTEGER :: noutearth ! 1/2 no of points on the E_odot curve
    INTEGER :: npo4,npo5
    INTEGER :: npts4,npts5
    INTEGER :: npts
! loop indexes
    INTEGER :: h,i,j
! steps
    DOUBLE PRECISION :: bstep,bstep1,bstep2,bstep3,bstep4,bstep5
! for the inside curve
    DOUBLE PRECISION :: radic
! flags
    INTEGER :: sampflag1,sampflag2
    INTEGER :: c1mez_flag != 1 if over  {rdot = -c1/2}
                          !=-1 if below {rdot = -c1/2}
! for uniforma
    INTEGER :: extrac(pmult*npox)
!
    DOUBLE PRECISION :: diff,diff4,diff5,newstep1,newstep2
    INTEGER :: numstep1,numstep2
    LOGICAL :: r0flag,crossflag,skipflag
! =====================================================================

! densifying the sampling of the boundary
    npo = npotmp*pmult
    IF(npo.lt.3) THEN
       write(*,*)'not enough points: npo =',npo
    ENDIF

! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! modified c(4) for strictly negative energy
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    DO i = 0,5
       c(i) = c0(i)
       !IF(i.eq.4) THEN
       !   c(i) = c0(i) - 0.5d0*E_bound
       !ENDIF
    ENDDO

! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! symmetry data of the region
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    rdotsym = -c(1)/2.d0  

    fr1 = f_func(r1,r2)
    fr2 = f_func(r2,r2)
    rdotfr1m = rdot_comp(c,fr1,-1,r2,E_bound)
    rdotfr1p = rdot_comp(c,fr1,+1,r2,E_bound)

    IF(fr2.le.fr1) THEN
       WRITE(*,*)'sample_bound: ERROR r1>r2',r1,r2    
       STOP
    ENDIF

    IF(rder0.le.r1) THEN
!       write(*,102)'fr1,fr2:',fr1,fr2
       npo1 = npo
       npo2 = 0
       bstep1 = (fr2-fr1)/(npo-1)
       bstep2 = 0.d0
    ELSEIF(rder0.gt.r1) THEN
!       WRITE(ierrou,*) name_obj(1:lobjnam),' min rdot at ',rder0
!       WRITE(*,*) name_obj(1:lobjnam),' min rdot at ',rder0
!       numerr=numerr+1
       frder0 = f_func(rder0,r2)
       rdotfrder0 = rdot_comp(c,frder0,-1,r2,E_bound)

       rdwtmp=(fr2-fr1)/(rdotfr1p-rdotfr1m) 

       npo1 = MAX(MIN(NINT(npo*(abs(frder0-fr1) + &
            & rdwtmp*abs(rdotfrder0 - rdotfr1m) ) &
            & /(abs(fr2-fr1) + rdwtmp*abs(rdotsym - rdotfr1m)) ),npo-3),2)

       if(npo-3.lt.2) then
          write(*,*) 'very poor sampling!'
       endif
       npo2 = npo - npo1
       bstep1 = (frder0-fr1)/(npo1-1)
!       bstep2 = (fr2-frder0)/(npo2-1)
       bstep2 = (fr2-frder0)/npo2
    ENDIF

102 FORMAT(a8,2x,f15.8,2x,f15.8)
103 FORMAT(a15,2x,f13.5,2x,f13.5,2x,f13.5)

!    write(*,*)'npo,npo1,npo2',npo,npo1,npo2
!    write(*,*)'bstep1=',bstep1,'  bstep2=',bstep2
!    write(*,*)'fr1,fr2:',fr1,fr2
!    write(*,*)'r1,r2:',r1,r2

! initialization
    sampflag1 = 1
    sampflag2 = 1
    r0flag = .false.
    crossflag =.false.

! ************************************************************         
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! nodes on the boundary E_{sun}= -gk**2/2*(a_orb) 
! (we consider elliptic orbits up to a semimajor axis a_orb)
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

! &&&&&&&&&&&&
! FIRST STEP
! &&&&&&&&&&&&
    frv1(1) = f_func(r1,r2)
    c1mez_flag = -1
    rdot1(1) = rdot_comp(c,frv1(1),c1mez_flag,r2,E_bound)

    frv2(1) = f_func(r1,r2)
    c1mez_flag = 1
    rdot2(1) = rdot_comp(c,frv2(1),c1mez_flag,r2,E_bound)

! &&&&&&&&&&&&&&
! MAIN DO LOOP
! &&&&&&&&&&&&&&
    DO i = 1,npo1-1
       frv1(i+1) = frv1(i) + bstep1
       c1mez_flag=-1
       rdot1(i+1) = rdot_comp(c,frv1(i+1),c1mez_flag,r2,E_bound)
       ! write(*,*)'frv1,rdot1',frv1(i+1),rdot1(i+1)         
         
       frv2(i+1) = frv2(i) + bstep1
       c1mez_flag=1
       rdot2(i+1) = rdot_comp(c,frv2(i+1),c1mez_flag,r2,E_bound)
       ! write(*,*)'frv2,rdot2',frv2(i+1),rdot2(i+1)         
    ENDDO

    IF(rder0.le.r1) THEN
       npts1 = npo-1 
    ELSEIF(rder0.gt.r1) THEN
!       npts1 = npo-2
       npts1 = npo-1
    ENDIF

    DO i = npo1,npts1 ! note that if (rder0.le.r1) this is skipped
       frv1(i+1) = frv1(i) + bstep2
       c1mez_flag=-1
       rdot1(i+1) = rdot_comp(c,frv1(i+1),c1mez_flag,r2,E_bound)
       ! write(*,*)'frv1,rdot1',frv1(i+1),rdot1(i+1)         
       
       frv2(i+1) = frv2(i) + bstep2
       c1mez_flag=1
       rdot2(i+1) = rdot_comp(c,frv2(i+1),c1mez_flag,r2,E_bound)
       ! write(*,*)'frv2,rdot2',frv2(i+1),rdot2(i+1)         
    ENDDO
    frvr2=f_func(r2,r2)

! symmetry point
    frv1(npo) = fr2 
    rdot1(npo) = rdotsym 
    frv2(npo) = frv1(npo)
    rdot2(npo) = rdotsym

! &&&&&&&&&&&&&&&&
! scaling factor
! &&&&&&&&&&&&&&&&
    rdw=(MAXVAL(frv2(1:npts1))-MINVAL(frv1(1:npts1)))/    &
         & (MAXVAL(rdot2(1:npts1))-MINVAL(rdot1(1:npts1))) 

! &&&&&&&&&&&&&&&&
! uniformization
! &&&&&&&&&&&&&&&&
    nout1 = min(max(npo1/pmult,3),npo1)
! first call to uniform,ato be done in both cases
    CALL uniform(frv1(1:npo1),rdw*rdot1(1:npo1),npo1,&
         & nout1,extrac,npox)

    DO i = 1,nout1
       frv1(i) = frv1(extrac(i))   
       rdot1(i) = rdot1(extrac(i))
       frv2(i) = frv2(extrac(i))   
       rdot2(i) = rdot2(extrac(i))
    ENDDO

! second call to uniforma, to be done if (rder0.gt.r1)
    IF(rder0.gt.r1) THEN 
       IF(npo2.eq.1) THEN
       ! skip uniformization: note that necessarily nout2 = npo2
          frv1(nout1+1) = frv1(npo)   
          rdot1(nout1+1) = rdot1(npo)
          frv2(nout1+1) = frv2(npo)   
          rdot2(nout1+1) = rdot2(npo)
          nout2=2
       ELSE
!          CALL uniform(frv1(npo1:npo-1),rdw*rdot1(npo1:npo-1),npo2,&
!               & nout2,extrac,npox)
          nout2 = min(max(npo2/pmult+1,3),npo2) !steep slope: one point more
          CALL uniform(frv1(npo1:npo),rdw*rdot1(npo1:npo),npo2+1,&
               & nout2,extrac,npox)
          DO i = 2,nout2
             frv1(nout1-1+i) = frv1(npo1-1+extrac(i))   
             rdot1(nout1-1+i) = rdot1(npo1-1+extrac(i))
             frv2(nout1-1+i) = frv2(npo1-1+extrac(i))   
             rdot2(nout1-1+i) = rdot2(npo1-1+extrac(i))
          ENDDO
       ENDIF
    ENDIF
    
    IF(rder0.le.r1) THEN
       nout = nout1
    ELSEIF(rder0.gt.r1) THEN 
       nout = nout1 + nout2 - 1
    ENDIF
    
! &&&&&&&&&&&&&&&&&&&&&&&&         
! re-ordering the points
! &&&&&&&&&&&&&&&&&&&&&&&&         
    DO j = 1,nout-1 
       frv1(nout+j) = frv2(nout-j) 
       rdot1(nout+j) = rdot2(nout-j)
    ENDDO
    noutsun = 2*nout - 1
       
!    write(*,*)'ordered points'
!    write(*,*)'frv, rdot'
    DO i = 1,noutsun
       frv(i) = frv1(i)
       rdot(i) = rdot1(i)
!       write(*,104) frv(i),rdot(i)
    ENDDO
!    write(*,*)''
104 FORMAT(f15.8,2x,f15.8)

! ****************************************************************
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! nodes on the boundary E_{earth}=0
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    me=1.d0/meinv
    r0 = (2.d0*gk**2*me/c0(2))**(1.d0/3.d0)
    fr0 = f_func(r0,r2)
! radius of the influence sphere of the Earth (assume a_earth = 1)      
    radSI = (me/3.d0)**(1.d0/3.d0)
    fradSI = f_func(radSI,r2)
    rmin = min(r0,radSI)
    frmin = min(fr0,fradSI)
    npoearth = max(NINT(npo*(frmin-fr1)/(fr2-fr1)),2)
    bstep3 = abs((frmin-fr1)/(npoearth-1))

!      write(*,*)'fr0=',fr0,'  r0=',r0
!      write(*,*)'fradSI=',fradSI
!      write(*,*)'rmin=',rmin,'frmin=',frmin
!      write(*,*)'npo',npo
!      write(*,*)'npoearth=',npoearth
!      write(*,*)'bstep3=',bstep3

    IF(rmin.gt.r1) THEN
! check if adding very close points
       IF(bstep3.gt.bstep1/5.d0) THEN
          ! initial step
          frv3(1) = frmin
          frv4(1) = frv3(1)
          radic=2.d0*gk**2*me/rmin - eta2*rmin**2
          IF(rmin.eq.r0)THEN
             rdot3(1)=0.d0
             rdot4(1)=0.d0
          ELSE
             rdot3(1) = -sqrt(radic)   
             rdot4(1) = -rdot3(1)
          ENDIF
          ! main do loop (from rmin to r1  <--- )
          DO h = 1,npoearth-1
             frv3(h+1) = frv3(h) - bstep3
             frv4(h+1) = frv3(h+1)
             radic=2.d0*gk**2*me/f_inv(frv3(h+1),r2) - &
                  & eta2*f_inv(frv3(h+1),r2)**2
             rdot3(h+1) = -sqrt(radic)   
             rdot4(h+1) = -rdot3(h+1)
          ENDDO
!           write(*,*)'frv3,rdot3,rdot4'
!           DO i = 1,npoearth
!              write(*,*) frv3(i),rdot3(i),rdot4(i)
!           ENDDO
!           write(*,*)''
          
          noutearth = max(npoearth/pmult,2)
          
          CALL uniform(frv3(1:npoearth),rdw*rdot3(1:npoearth),npoearth,&
               & noutearth,extrac,npox)
          DO i = 1,noutearth
             frv3(i) = frv3(extrac(i))   
             rdot3(i) = rdot3(extrac(i))
             frv4(i) = frv4(extrac(i))   
             rdot4(i) = rdot4(extrac(i))
          ENDDO

       ELSEIF(bstep3.le.bstep1/5.d0) THEN
          ! basestep too small: draw a straight line 
          ! from rdot2(1) to rdot1(1)
          rdot3(1) = rdot2(1)  
          frv3(1) = fr1
       ENDIF

! this should not happen
    ELSEIF(rmin.le.r1) THEN
       ! write(*,*)'sample_bound: min(radSI,r0) < r1 !!!'
       ! draw a straight line from rdot2(1) to rdot1(1)
          rdot3(1) = rdot2(1)  
          frv3(1) = fr1

    ENDIF

! ***************************************************************
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! adding nodes on the vertical line r = r1
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    diff = abs(rdot2(1) - rdot1(1))

    IF((bstep3.gt.bstep1/5.d0).and.(rmin.gt.r1)) THEN
       diff4 = abs(rdot3(noutearth) - rdot1(1))
       diff5 = abs(rdot4(noutearth) - rdot2(1))
       npo4 = max(NINT(nout*diff4/diff),2)
       npo5 = max(NINT(nout*diff5/diff),2)
       bstep4 = diff4/(npo4-1)
       bstep5 = diff5/(npo5-1)
       
       DO h = noutearth,noutearth+npo4-1
          rdot3(h+1) = rdot3(h) - bstep4
          frv3(h+1) = frv3(h)
       ENDDO
       DO h = noutearth,noutearth+npo5-1
          rdot4(h+1) = rdot4(h) + bstep5
          frv4(h+1) = frv4(h)
       ENDDO
       
       npts4 = noutearth + npo4 -1
       npts5 = noutearth + npo5 -1
      
! &&&&&&&&&&&&&&&&&&&&&&&&
! re-ordering the points
! &&&&&&&&&&&&&&&&&&&&&&&&
       DO j = 1,npts5-1
          frv(noutsun+j) = frv4(npts5-j) 
          rdot(noutsun+j) = rdot4(npts5-j)
       ENDDO

       IF((r0.lt.radSI).or.(radSI.le.r1)) THEN
          npts = noutsun+npts5-2
       ELSE
          npts = noutsun+npts5-1
       ENDIF
       
       DO j = 1,npts4-1
          frv(npts+j) = frv3(j) 
          rdot(npts+j) = rdot3(j)
       ENDDO
       
! &&&&&&&&&&&&&&&&&&&&&&&
! total number of nodes
! &&&&&&&&&&&&&&&&&&&&&&&
       ntot = npts + npts4 - 1

! in order to avoid adding of very close points 
! draw a straight line ignoring the E_{Earth} curve 
    ELSEIF((bstep3.le.bstep1/5.d0).or.(rmin.le.r1)) THEN   
!       write(*,*)'*** draw a vertical line ***'
       bstep4 = diff/(nout-1) 
            
       DO h = 1,nout-1 
          rdot3(h+1) = rdot3(h) - bstep4
          frv3(h+1) = fr1
       ENDDO
       
! &&&&&&&&&&&&&&&&&&&&&&&&
! re-ordering the points
! &&&&&&&&&&&&&&&&&&&&&&&&
       DO j = 1,nout-2 
          frv(noutsun+j) = frv3(j+1) 
          rdot(noutsun+j) = rdot3(j+1)
       ENDDO
       
! &&&&&&&&&&&&&&&&&&&&&&&
! total number of nodes
! &&&&&&&&&&&&&&&&&&&&&&&
       ntot = noutsun + nout - 2 

    ENDIF

! to skip last uniformization
!    GOTO 111

! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! check for number of very close points
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    ntmp = ntot
    frvtmp(1) = frv(1) 
    rdottmp(1) = rdot(1)
    DO j =1,ntot-1 
       frvtmp(j+1) = frv(j+1) 
       rdottmp(j+1) = rdot(j+1)
       smallpar = sqrt((rdw*diff)**2 + (abs(fr1-fr2))**2)/(ntot*1.d2)
       dist = sqrt((frv(j+1)-frv(j))**2 + (rdot(j+1)-rdot(j))**2)
       IF(dist.le.smallpar) THEN
          ntmp = ntmp-1           
       ENDIF
    ENDDO
!    write(*,*)rdot(1:ntot) 
!    write(*,*)'ntot,ntmp',ntot,ntmp

! &&&&&&&&&&&&&&&&&&&&&&&
! last call to uniforma 
! &&&&&&&&&&&&&&&&&&&&&&&
    CALL uniform(frvtmp(1:ntot),rdw*rdottmp(1:ntot), &
         & ntot,ntmp,extrac,npox)
    DO i = 1,ntmp
       frv(i) = frvtmp(extrac(i))   
       rdot(i) = rdottmp(extrac(i))
    ENDDO
    ntot = ntmp
!    write(*,*)rdot(1:ntot) 
 
111 CONTINUE
    
  END SUBROUTINE sample_bound_ne1

! ==================================
! computation of rdot, given f(r)
! ==================================
    DOUBLE PRECISION FUNCTION rdot_comp(c,fr,c1mez_flag,r2,E_bound)
      USE fund_const
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) ::  c(0:5) ! polynomial coefficients
      DOUBLE PRECISION, INTENT(IN) :: fr 
      INTEGER, INTENT(IN) :: c1mez_flag !=1 if over -c1/2
                                        !=-1 if below -c1/2
      DOUBLE PRECISION, INTENT(IN) :: r2 
      DOUBLE PRECISION, INTENT(IN), OPTIONAL :: E_bound
! END INTERFACE
! =====================================================================
!      DOUBLE PRECISION :: f_inv
      DOUBLE PRECISION :: r
      DOUBLE PRECISION :: bigC,bigS,delta
      DOUBLE PRECISION :: E_bound_loc
! =====================================================================

      r = f_inv(fr,r2)

      IF(PRESENT(E_bound)) THEN
         E_bound_loc = E_bound
      ELSE
         E_bound_loc = 0.d0
      END IF
      bigC = c(2)*r**2 + c(3)*r + c(4)
      bigS = r**2 + c(5)*r + c(0)
      IF(E_bound_loc.GT.0) THEN
         delta = c(1)**2 - 4.d0*bigC + 8.d0*gk**2/sqrt(bigS)+2.d0*E_bound
      ELSE
         delta = c(1)**2 - 4.d0*bigC + 8.d0*gk**2/sqrt(bigS)
      END IF

      IF(delta.lt.0.d0)THEN
         IF(abs(delta).lt.1.d-8) THEN
            delta = 0.d0
         ELSEIF((abs(delta).ge.1.d-8).and.(abs(delta).lt.1.d-3)) THEN
            WRITE(*,*)'rdot_comp: WARNING negative discriminant=',delta
            delta = 0.d0
         ELSE
            WRITE(*,*)'rdot-comp: ERROR negative discriminant=',delta
            delta=0.d0
         ENDIF
      ENDIF

      IF(c1mez_flag.eq.1) THEN
         rdot_comp = -(c(1) - dsqrt(delta))/2.d0
      ELSEIF(c1mez_flag.eq.-1) THEN
         rdot_comp = -(c(1) + dsqrt(delta))/2.d0
      ELSE
         WRITE(*,*)'ERROR in rdot_comp'
         STOP
      ENDIF

    END FUNCTION rdot_comp


! =====================================================================
! ====================================
! function used to define the metric
! ====================================
    DOUBLE PRECISION FUNCTION f_func(r,sigma)
      USE fund_const
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) :: r 
      DOUBLE PRECISION, INTENT(IN) :: sigma 
! ============= end interface =========================================
      IF(nfunc.eq.1)THEN      
! adaptive metric
         f_func  = 1.d0 - exp(-r**2/(2.0*(sigma)**2))
      ELSEIF(nfunc.eq.2)THEN
! near earth
         f_func  = log10(r)
      ELSEIF(nfunc.eq.3)THEN
! simplest metric
         f_func = r
      ELSE
         WRITE(*,*)' function used to define triangulation metric unknown, nfunc= ', nfunc
         STOP
      ENDIF    
    END FUNCTION f_func
    
! =====================================================================
! ====================================
! inverse of the function used 
! to define the metric
! ====================================
    DOUBLE PRECISION FUNCTION f_inv(y,sigma)
      USE fund_const
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) :: y 
      DOUBLE PRECISION, INTENT(IN) :: sigma 
! ============= end interface =========================================
      IF(nfunc.eq.1)THEN 
! adaptive metric
         IF(y.le.0.d0) THEN 
            WRITE(*,*)'f_inv: negative argument for sqrt! y:',y
            f_inv=0.d0
         ELSEIF(y.ge.1.d0)THEN
            WRITE(*,*)'f_inv: infinite argument for sqrt! y:',y
            f_inv=1.d10
         ELSE
            f_inv  = sqrt(2.0*(sigma)**2*log(1.d0/(1.d0-y)))
         ENDIF
      ELSEIF(nfunc.eq.2)THEN
! near earth
         f_inv  = 10.d0**y
      ELSEIF(nfunc.eq.3)THEN
! simplest metric
         f_inv = y 
      ELSE
         WRITE(*,*)' function used to define triangulation metric unknown, nfunc= ', nfunc
           STOP
        ENDIF
      END FUNCTION f_inv

! =====================================================================
! =======================================
! value of f(r) corresponding to rdot=0 
! =======================================
     DOUBLE PRECISION FUNCTION f_r_rdot0(c,r2)
      USE fund_const
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) ::  c(0:5) ! polynomial coefficients
      DOUBLE PRECISION, INTENT(IN) ::  r2
! =========== end interface ===========================================
!      DOUBLE PRECISION :: f_func
      DOUBLE PRECISION :: rdot !current rdot
      DOUBLE PRECISION :: d(0:6),alpha
!     to call find_pos_roots
!      DOUBLE PRECISION :: radius(6)
      LOGICAL :: hzflag ! hzflag = .true.  OK!
                        ! hzflag = .false. abs(root)>10^5
      LOGICAL :: multfl ! multfl = .true.  OK!
                        ! multfl = .false. 0 has multiplicity > 4
      DOUBLE PRECISION :: roots(6)
      DOUBLE PRECISION :: frtmp(6)
      INTEGER nroots
      INTEGER :: h ! loop index
! =====================================================================

      nroots = 0 
      rdot = 0.d0
      alpha = rdot**2 + c(1)*rdot
      
      d(0) = 2*c(0)*alpha*c(4)-4*gk**4+c(0)*c(4)**2+c(0)*alpha**2

      d(1) = 2*c(5)*alpha*c(4)+2*c(0)*c(3)*c(4)+2*c(0)*alpha*c(3)+ &
           & c(5)*alpha**2+c(5)*c(4)**2
      
      d(2) = 2*c(5)*c(3)*c(4)+2*c(5)*alpha*c(3)+2*c(0)*c(2)*c(4)+ &
           & 2*c(0)*alpha*c(2)+2*alpha*c(4)+c(0)*c(3)**2+c(4)**2+alpha**2
      
      d(3) = 2*c(5)*c(2)*c(4)+2*c(5)*alpha*c(2)+2*c(0)*c(2)*c(3)+ &
           & c(5)*c(3)**2+2*alpha*c(3)+2*c(3)*c(4)
      
      d(4) = 2*alpha*c(2)+2*c(5)*c(2)*c(3)+2*c(2)*c(4)+ &
           & c(0)*c(2)**2+c(3)**2
      
      d(5) = c(5)*c(2)**2+2*c(2)*c(3)
      
      d(6) = c(2)**2
      
! compute positive real roots
      CALL find_pos_roots(6,d,roots,nroots,hzflag,multfl)
!      CALL solv6(d,roots,nroots,radius)
!      write(*,*)'roots',roots(1:nroots)
      IF(nroots.eq.0) THEN
         WRITE(*,*)'ERROR in samplebound: nroots=',nroots
         STOP
      ELSEIF(nroots.eq.1) THEN
         f_r_rdot0  = f_func(roots(1),r2)
      ELSEIF((nroots.gt.1).and.(nroots.le.3)) THEN
         f_r_rdot0  =  f_func(MINVAL(roots(1:3)),r2)
      ELSE
         WRITE(*,*)'ERROR in sample_bound: nroots=',nroots
         STOP
      ENDIF
      
    END FUNCTION f_r_rdot0
END MODULE triangles

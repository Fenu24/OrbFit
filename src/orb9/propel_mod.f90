MODULE propel_mod
  USE name_rules, ONLY: name_len
  IMPLICIT NONE
  PRIVATE

! proper elements, synthetic version
TYPE proper_els
  CHARACTER*(name_len) :: name
  DOUBLE PRECISION               :: hmag  ! absol. mag
  DOUBLE PRECISION, DIMENSION(3) :: pr_el ! a(au), e, sin I
  DOUBLE PRECISION, DIMENSION(3) :: pr_fr ! n, g, s in rad/y
  DOUBLE PRECISION               :: lce   ! max. Lyap exp. (1/y)
  DOUBLE PRECISION               :: t_int ! integration time span (y) 
END TYPE proper_els

TYPE close_couple
  CHARACTER*(name_len), DIMENSION(2) :: names
  INTEGER, DIMENSION(2)              :: addr
  DOUBLE PRECISION, DIMENSION(3)     :: delta_el ! da/a, de, dsinI
  DOUBLE PRECISION                   :: sigma
!  DOUBLE PRECISION, DIMENSION(3,2)   :: pr_el ! a(au), e, sin I
!  DOUBLE PRECISION, DIMENSION(3,2)   :: pr_fr ! n, g, s (rad/y)
  DOUBLE PRECISION, DIMENSION(3,2)   :: angles ! omegatilda, Omega, elle (rad)
END TYPE close_couple

TYPE familyatt
  CHARACTER*(name_len) name
  DOUBLE PRECISION hmag
  INTEGER status ! 0 no known fam 1 second step family 2 attributed to core 3 core 4 multiple
  CHARACTER*(name_len), DIMENSION(2) :: family
  INTEGER, DIMENSION(2) :: famno
  CHARACTER*(name_len), DIMENSION(2) ::  near
  DOUBLE PRECISION, DIMENSION(2) :: dv
  CHARACTER*6 ::  rescode
END type familyatt

PUBLIC proper_els, close_couple, familyatt

! public routines
PUBLIC input_propels, taxostep, input_families, input_propeltro  
!
! shared data
!
! proper elements
!
INTEGER, PARAMETER :: nprox=800000 ! max size of array
TYPE(proper_els), DIMENSION(nprox) :: propel ! storage array
INTEGER    :: npro ! actual number
! pointers for sort by sine of inclination, by name
INTEGER, DIMENSION(nprox) :: isrti, isrtnam 
!
! couples
!
INTEGER, PARAMETER :: ncoupx=1000000 ! max size of array
TYPE(close_couple), DIMENSION(ncoupx) :: clos 
INTEGER, DIMENSION(ncoupx) :: isrts ! sort by sigma
!
! families
!
! list of all members of all families
CHARACTER*(name_len), DIMENSION(nprox) :: memberlist 
! max number of families, of members in one family
INTEGER, PARAMETER :: nfamx=1000,nmembx=10000 
! pointer in memberlist to first of family, number of members in family
INTEGER, DIMENSION(nfamx) :: fampoint, ninfam 
! actual number of families, of members for all families
INTEGER  :: nfam, nmemb 
! fix for numeric code of multiopposition used in families lists 
INTEGER, PARAMETER :: ncodx=100000
INTEGER icod(ncodx), nconv  ! numeric codes, total number of them 
CHARACTER*(name_len) namul(ncodx) ! corresponding alphanumeric names

!
! classification in families
!
TYPE(familyatt), DIMENSION(nprox) :: famrec 

PUBLIC ncodx, icod, nconv, namul
PUBLIC nprox, propel, npro, isrti, isrtnam
PUBLIC nfamx,nmembx,nfam,nmemb,fampoint, ninfam, memberlist
PUBLIC ncoupx, clos,isrts, famrec
 
CONTAINS

! input of family list
 SUBROUTINE input_families(iun,err_line)
  INTEGER, INTENT(IN)  :: iun ! input unit
  INTEGER, INTENT(OUT) :: err_line ! error flag, if >0 is error location
! end interface 
  INTEGER j,k,l,n,nline, le, lench ! loop indexes, record counter, character counter
  CHARACTER*6 intmemb(13)
  CHARACTER*256 record
  INTEGER icode,icode6
  EXTERNAL lench
! main loop on families
  err_line=0
  nline=0   
  DO j=nfam+1,nfamx
    IF(nfam.eq.0)THEN
       fampoint(j)=1
    ELSEIF(nfam.gt.0)THEN
       fampoint(j)=fampoint(j-1)+ninfam(j-1)
    ELSE
       STOP '*** nfam wrong*****'
    ENDIF
    ninfam(j)=0
    DO k=1,nmembx,13
       READ(iun,101,END=3)record
101    FORMAT(A256)
       le=lench(record)
       IF(le.eq.0)THEN
! blank line means family ended
          EXIT
       ENDIF
       READ(record,100)intmemb
100    FORMAT(A6,12(1X,A6)) 
       nline=nline+1
       DO l=13,0,-1
         IF(intmemb(l).ne.'     0')THEN
            DO n=1,l
! fix for multiopposition added as numbers above 700000
              IF(intmemb(n).eq.'703259')THEN
                 memberlist(nmemb+n)='2001XS146'
              ELSEIF(intmemb(n).eq.'761869')THEN
                 memberlist(nmemb+n)='2010OX76'
              ELSEIF(intmemb(n).eq.'731055')THEN
                 memberlist(nmemb+n)='2017AA'
              ELSEIF(intmemb(n).eq.'747159')THEN
                 memberlist(nmemb+n)='2017BB'
              ELSE
! read a 6-digit integer from the string 
                 READ(intmemb(n),*)icode6
                 IF(icode6.gt.700000)THEN
! this must be a multiopposition given with a numeric code
                    icode=icode6-700000
                    IF(icode.gt.nconv)THEN
                       WRITE(*,*)' input_families: too many multiopp in conversion list?', icode,' ',nconv
                       STOP
                    ENDIF
                    memberlist(nmemb+n)=namul(icode)
                 ELSE
! this is a numbered asteroid (WARNING: numbers up to 699999)
                    memberlist(nmemb+n)='   '//intmemb(n)
                 ENDIF
              ENDIF
              CALL rmsp(memberlist(nmemb+n),le)
            ENDDO
            nmemb=nmemb+l
            EXIT
         ENDIF
       ENDDO
       IF(l.gt.0.and.l.lt.13)THEN
          ninfam(j)=ninfam(j)+l
! zeros in the line means family ended, but skip one blank line
          READ(iun,101,END=3)record
          le=lench(record)
          IF(le.eq.0)THEN
! should be blank line
             EXIT
          ELSE
! error in file
             err_line=nline
             RETURN
          ENDIF
       ELSEIF(l.eq.13)THEN
          ninfam(j)=ninfam(j)+13
! read another record
       ELSEIF(l.eq.0)THEN
          WRITE(*,*)' line of zeros at line ', nline
          err_line=nline
       ENDIF
    ENDDO
! end of family, skipped one line already
    nfam=nfam+1
 ENDDO
! can get here only if too many families?
 WRITE(*,*)' too many families, nfam=',nfam,' nfamx=',nfamx
 STOP
3 WRITE(*,*)' end k=',k, ' nfam=', nfam
 IF(k.gt.1)THEN
    nfam=nfam+1
 ENDIF
END SUBROUTINE input_families
! input of synthetic proper element
 SUBROUTINE input_propels(iun,npr,err_line,secre)
  INTEGER, INTENT(IN)  :: iun ! input unit
  INTEGER, INTENT(INOUT)  :: npr ! previously assigned in input
  INTEGER, INTENT(OUT) :: err_line ! error flag, if >0 is error location
  LOGICAL, INTENT(IN), OPTIONAL :: secre ! T=secular resonant F=sordinary syntetic
! 
  CHARACTER*256 record
  INTEGER j
  CHARACTER*(name_len) nam
  DOUBLE PRECISION hmag, pr_el(3),pr_fr(3),lce,t_int
  LOGICAL secrenow
! skip header
  READ(iun,*)
  READ(iun,*)
! initialization
  npro=npr
  err_line=0
! secular resonant case?
  IF(PRESENT(secre))THEN
     secrenow=secre
  ELSE
     secrenow=.false.
  ENDIF
! main loop
  DO j=npr+1,nprox
    READ(iun,'(A)',END=3) record
    READ(record,*,ERR=2)nam, hmag,pr_el, pr_fr, lce, t_int
    npro=npro+1
    propel(npro)%name=nam
    propel(npro)%hmag=hmag
    If(secrenow)pr_el(2)=pr_el(2)+1.d0
    propel(npro)%pr_el=pr_el
    propel(npro)%pr_fr=pr_fr
    propel(npro)%lce=lce
    propel(npro)%t_int=t_int
  ENDDO 
! too low nprox
  WRITE(*,*)' increase nprox, was ',nprox
  npr=npro 
  RETURN
2 WRITE(*,*)' read error in proper elements file, line ',j
  err_line=j
  npr=npro  ! stored anyway
  RETURN
! regular ending at eof
3 WRITE(*,*)' read ',npro,' proper elements'
  npr=npro
 END SUBROUTINE input_propels
 
! input of synthetic proper element
 SUBROUTINE input_propeltro(iun,iuns,npr,err_line)
  USE fund_const, ONLY: gk
  INTEGER, INTENT(IN)  :: iun, iuns ! input units (proper elements, sigmas)
  INTEGER, INTENT(INOUT)  :: npr ! previously assigned in input
  INTEGER, INTENT(OUT) :: err_line ! error flag, if >0 is error location
! 
  CHARACTER*256 record, records
  INTEGER j
  CHARACTER*(name_len) nam, nams
  DOUBLE PRECISION hmag, pr_el(3),pr_fr(3),lce,t_int
! specific for Trojans
  INTEGER l45, my, outsyn, outsig
  DOUBLE PRECISION da, cr, f, pre, g, prsi, s
  DOUBLE PRECISION sda, scr, sf, spre,sg,sprsi, ss
  DOUBLE PRECISION :: nj
  DOUBLE PRECISION, PARAMETER :: aj=5.2025696d0
! inizialization of nj  
  nj=gk*sqrt(1.d0/aj**3)
! remove  headers
  READ(iun,*)
  READ(iun,*)
  READ(iuns,*)
  READ(iuns,*)
  READ(iuns,*)
! initialization
  npro=npr
  err_line=0
! open output files 
  CALL filopn(outsyn,'tromod.syn','unknown')
  CALL filopn(outsig,'tromod.sig','unknown')
! main loop
  DO j=npr+1,nprox
    READ(iun,'(A)',END=3) record
    READ(record,*,ERR=2)nam,hmag,da,cr,f,pre,g,prsi,s,l45,my
    t_int=my
    READ(iuns,'(A)',END=5) records
    READ(records,*,ERR=4)nams,sda,scr,sf,spre,sg,sprsi,ss,lce
    IF(nam.ne.nams)THEN
       WRITE(*,*)' input_propeltro: inconsistent inoput in .syn and .sig: ', nam,' ',nams
       STOP
    ENDIF 
    npro=npro+1
! insert in array propel
    propel(npro)%name=nam
    propel(npro)%hmag=hmag
    IF(l45.eq.4)THEN
       propel(npro)%pr_el(1)=aj-da
       propel(npro)%pr_fr(1)=nj+f
    ELSEIF(l45.eq.5)THEN
       propel(npro)%pr_el(1)=aj+da
       propel(npro)%pr_fr(1)=nj-f
    ELSE
       WRITE(*,*)' input_propeltro: l45=', l45
       STOP ' ***** error in Trojan proper elements input *********'
    ENDIF
    propel(npro)%pr_el(2)=pre
    propel(npro)%pr_fr(2)=g
    propel(npro)%pr_el(3)=prsi
    propel(npro)%pr_fr(3)=s
    propel(npro)%lce=lce
    propel(npro)%t_int=t_int
! output in tromod.syn, tromod.sig
    WRITE(outsyn,110)propel(npro)%name,propel(npro)%hmag,propel(npro)%pr_el, propel(npro)%pr_fr, lce, my
110 FORMAT(a9,1x,F5.2,1X,f11.7,1x,f10.7,1x,f9.7,1x,f11.6,1x,f12.6,1x,f12.6,1x,f7.2,1x,I3)
!  rmsang0(1) is the RMS of residuals w.r. to linear fit of mean longitude, not available for Trojans
    WRITE(outsig,112)propel(npro)%name,sda,spre,sprsi,sf,sg,ss, 0.d0, my
112 FORMAT(a9,1x,f11.8,1x,f11.7,1x,f9.7,1x,f9.6,1x,f10.6,1x,f9.6,1x,f12.4,1x,I3)
  ENDDO 
! too low nprox
  WRITE(*,*)' increase nprox, was ',nprox
  npr=npro 
  GOTO 10
2 WRITE(*,*)' read error in proper elements file, line ',j
  err_line=j
  npr=npro  ! stored anyway
  GOTO 10
! errors in tro.sig
4 WRITE(*,*)' read error in proper sigmas file, line ',j
  err_line=j
  npr=npro  ! stored anyway
  GOTO 10
! errors in tro.sig
5 WRITE(*,*)' early termination in proper sigmas file, line ',j
  err_line=j
  npr=npro  ! stored anyway
  GOTO 10
! regular ending at eof
3 WRITE(*,*)' read ',npro,' proper elements'
  npr=npro
! close files
10 CONTINUE
  CALL filclo(outsyn,' ')
  CALL filclo(outsig,' ')
 END SUBROUTINE input_propeltro

  SUBROUTINE taxostep(ii,pro,nam0,hmag,d_max,n1)
    USE fund_const, ONLY: gk
    INTEGER, INTENT(IN)          :: ii     ! address in array of pro
    DOUBLE PRECISION, INTENT(IN) :: pro(3), hmag  ! current proper elements, magnitude
    CHARACTER*(*), INTENT(IN)    :: nam0    ! orbit name
    DOUBLE PRECISION, INTENT(IN) :: d_max   ! control in m/s
    INTEGER, INTENT(INOUT)       :: n1      ! number of filter 1 ids
! ==============end interface============================
    DOUBLE PRECISION, PARAMETER :: wa=1.25d0, we=2.d0, wi=2.d0 !metrics
    INTEGER indal
    DOUBLE PRECISION d2_maxna, vmin, vmax, dmaxna, ennea
!
    vmax=pro(3)+d_max/wi
    vmin=pro(3)-d_max/wi
! mean orbital velocity in m/s
    ennea=gk/sqrt(pro(1))*149597870.7d3/86400.d0
! controls divided by mean velocity
    dmaxna=d_max/ennea
    d2_maxna=dmaxna*dmaxna
    CALL bin_search_rea(pro(3),npro,1,npro,propel(1:npro)%pr_el(3),isrti,indal)
    CALL scan_list(.true.)
    CALL scan_list(.false.)
  CONTAINS
! =========================================                             
!  SCAN_LIST
! ========================================= 
  SUBROUTINE scan_list(forward)
    LOGICAl, INTENT(IN) :: forward
! end interface
    INTEGER lmin, lmax, lstep, na, ni,i,nn1 ! loop indexes and control
    DOUBLE PRECISION di,de,da, sigma2,sigma
    CHARACTER*(name_len) nam
    INTEGER le, ll ! character manipulations
! find extremes of loop
    IF(forward)THEN
       lmin=indal+1
       lmax=npro
       lstep=1
    ELSE
       lmin=indal
       lmax=1
       lstep=-1
    ENDIF
! loop on list sorted by alpha
    DO 22 na=lmin, lmax, lstep ! na is the index of the list sorted by alpha
       ni=isrti(na) ! ni is the index of the list sorted by sinI
! remove duplication, keep with the first earlier in the original file
       IF(ni.le.ii) CYCLE
! test for distance in sinI
          di=propel(ni)%pr_el(3)
       IF(forward)THEN
          IF(di.gt.vmax) RETURN
       ELSE
          IF(di.lt.vmin) RETURN
       ENDIF
       di=di-pro(3)
       de=propel(ni)%pr_el(2)-pro(2)
       IF(ABS(de*we).gt.dmaxna) CYCLE
       da=(propel(ni)%pr_el(1)-pro(1))/pro(1)
       IF(ABS(da*wa).gt.dmaxna) CYCLE
       sigma2=(da*wa)**2+(de*we)**2+(di*wi)**2
       IF(sigma2.gt.d2_maxna) CYCLE
       sigma=sqrt(sigma2)*ennea
! check for self identification
       IF(ii.eq.ni) CYCLE
! store for possible relationship
       IF(n1.eq.ncoupx)THEN
          WRITE(*,*)' increase ncoupx, was ',ncoupx
       ELSE
          n1=n1+1
          clos(n1)%names(1)=nam0
          clos(n1)%names(2)=propel(ni)%name
          clos(n1)%addr(1)=ii
          clos(n1)%addr(2)=ni
          clos(n1)%sigma=sigma
          clos(n1)%delta_el(1)=da
          clos(n1)%delta_el(2)=de
          clos(n1)%delta_el(3)=di
      ENDIF
22  ENDDO
  END SUBROUTINE scan_list
END SUBROUTINE taxostep



END MODULE propel_mod

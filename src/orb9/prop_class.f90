!-------------------------------------------------------------------
! PROP_CLASS
!
! Classification in existing families
!  version 4.3, 1 December 2012
! with only 1 step of nearest neighbour
! but including (from vers. of Feb. 27, 2013) step 2 families
!  as separate families (later to be merged)
! ------------------------------------------------------------------
PROGRAM prop_class
 USE propel_mod
 USE name_rules, ONLY: name_len
 IMPLICIT NONE
 CHARACTER*6 progna ! option control: progname
 CHARACTER*80 run ! option control: current run 
 INTEGER ler,le ! character counters
 INTEGER in_pro, in_fam, in_res, in_sig ! input units
 INTEGER iun_log,iun_fam,iun_tab, iun_rec ! output units
 INTEGER err_line, npr
 LOGICAL isthere, istheres
! loop indexes and counters
 INTEGER n1, n1f, j, i, i1, i2, jj, jjj
! counters for zones
 INTEGER nfam2, nfam3,nfam4, nfam2i, nfam3i,nfam4i,nfam5,nfam1,nfam6
 INTEGER nfamh2,nfamh3,nfamh4,nfam7,nfam8,nfam9
 INTEGER kz, nf, nf1, nadd, ndouble
 DOUBLE PRECISION box(3,2,nfamx)
 INTEGER clopoint(nfamx),naddfam(nfamx),nnewfam(nfamx),nsmall(nfamx),ncorefam(nfamx)
 DOUBLE PRECISION d_max ! control
 DOUBLE PRECISION pro(3), hmag, ap
 CHARACTER*(name_len) name, famname, famname1, merg(2,nfamx),merg_from, merge_to
 INTEGER in_merg, nmerg, nm, ist1,ist2
 INTEGER istat_mer(2,nfamx),famno_merg(2,nfamx), status
 INTEGER merge_flag(nfamx) ! 0=no merge, else nf of parent family 
 INTEGER r1,r2,r3,r4,pl,tp,pa,l45,pl1,pl2,nr2,nr3 ! resonances 2body-3body
 INTEGER in_conver,icode
 CHARACTER*(name_len) nam
 CHARACTER*1 pllab
 TYPE(familyatt) f
 LOGICAL secre ! T when reading secular resonant proper elements
! **********BEGIN EXEXCUTION******************************** 
! progna='??'
 WRITE(*,*)' run name?'
 READ(*,*) run 
 CALL rmsp(run,ler)
! open log file
 CALL filopn(iun_log,run(1:ler)//'.clarep','unknown')
! -------------------------------------------------------------
! read proper elements
! -------------------------------------------------------------
 npr=0
! open input file, numbered 
 INQUIRE(file='numb_res.syn', exist=isthere)
 IF(isthere)THEN
    CALL filopn(in_pro,'numb_res.syn','old')
    CALL  input_propels(in_pro,npr,err_line)
    WRITE(iun_log,*)' input ',npr,' prop. els of numb. main belt, including resonant'
! error return?
    IF(err_line.gt.0)THEN
       WRITE(*,*)' error return in proper elements numbered from record ',err_line
       STOP
    ENDIF
    CALL filclo(in_pro,' ')
 ELSE
    WRITE(iun_log,*)'no proper elements for numbered with secres'
    STOP
 ENDIF
!
! open input file, Trojans
 INQUIRE(file='tro.syn', exist=isthere)
 INQUIRE(file='tro.sig', exist=istheres)
 IF(isthere.and.istheres)THEN
    CALL filopn(in_pro,'tro.syn','old')
    CALL filopn(in_sig,'tro.sig','old')
    CALL  input_propeltro(in_pro,in_sig,npr,err_line)
    WRITE(iun_log,*)' input ',npr,' prop. els of main belt and Trojans'
!  error return?
    IF(err_line.gt.0)THEN
       WRITE(*,*)' error return in trojan proper elements from record ',err_line
       STOP
    ENDIF
    CALL filclo(in_pro,' ')
    CALL filclo(in_sig,' ')
 ELSE
    WRITE(iun_log,*)'no proper elements for Trojans'
    STOP
 ENDIF
! open input file, multiopposition
 INQUIRE(file='mult_res.syn', exist=isthere)
 IF(isthere)THEN
    CALL filopn(in_pro,'mult_res.syn','old')
    CALL  input_propels(in_pro,npr,err_line)
    WRITE(iun_log,*)' input ',npr,' prop. els including multioppostion'
! error return?
    IF(err_line.gt.0)THEN
       WRITE(*,*)' error return in proper elements multiopp from record ',err_line
       STOP
    ENDIF
    CALL filclo(in_pro,' ')
 ELSE
    WRITE(iun_log,*)'no proper elements for multiopposition'
 ENDIF

! -------------------------------------------------------------
! initialize family records for all the asteroids with proper elements 
! -------------------------------------------------------------
 famrec(1:npr)%name=propel(1:npr)%name
 famrec(1:npr)%hmag=propel(1:npr)%hmag
 famrec%status=0
 DO i=1,2
    famrec%family(i)='0         '
    famrec%famno(i)=0
    famrec%near(i)='0        '
    famrec%dv(i)=0.d0
 ENDDO
 famrec%rescode='0     '
! -------------------------------------------------------------
! sorting
!
! sort by sinI all asteroids
!
! -------------------------------------------------------------
 CALL heapsort(propel(1:npro)%pr_el(3),npro,isrti)
 WRITE(iun_log,*)' end sorting by sinI'
! 
! sort by name; used for resonances
!
 CALL heapsortnach(propel(1:npro)%name,npro,isrtnam)
 WRITE(iun_log,*)' end sorting by name'
! -------------------------------------------------------------
! add resonance flag
! -------------------------------------------------------------
 CALL filopn(in_res,'famdat/twobody.fla','old')
 READ(in_res,*)
 nr2=0
 DO j=1,npro ! resonant are less than proper elements
    READ(in_res,*,END=19)jj,r1,r2,r3,r4,pl,tp,pa,l45
    WRITE(name,'(I9)')jj
    CALL rmsp(name,le)
    CALL bin_search(name,npro,propel(1:npro)%name,isrtnam,jj)
    IF(jj.le.0)THEN
!       WRITE(*,*)name,' resonant but has no proper elements'
    ELSE
       nr2=nr2+1
       IF(pl.eq.5)THEN
          pllab='J'
       ELSEIF(pl.eq.4)THEN
          pllab='M'
       ELSE
          pllab=' '
          WRITE(*,*)name, pl,' resonance with unknown planet'
       ENDIF
       WRITE(famrec(jj)%rescode,199)r1,abs(r2),pllab
 199   FORMAT(I2,'/',I1,' ',A1)
       CALL rmsp(famrec(jj)%rescode,le)
    ENDIF
 ENDDO
19 CALL filclo(in_res,' ')
 WRITE(*,*)' read twobody resonances, with proper elements are ',nr2
 CALL filopn(in_res,'famdat/threebody.fla','old')
 READ(in_res,*)
 nr3=0
 DO j=1,npro ! resonant are less than proper elements
    READ(in_res,*,END=29)jj,r1,r2,r3,pl1,pl2,tp,pa
    WRITE(name,'(I9)')jj
    CALL rmsp(name,le)
    CALL bin_search(name,npro,propel(1:npro)%name,isrtnam,jj)
    IF(jj.le.0)THEN
!       WRITE(*,*)name,' resonant but has no proper elements'
    ELSE
       nr3=nr3+1
       IF(pl1.eq.5.and.pl2.eq.6)THEN
!   ok default
       ELSE
          WRITE(*,*)name, pl1, pl2,' 3-body resonance with unknown planets'
       ENDIF
       WRITE(famrec(jj)%rescode,299)r1,r2,r3
 299   FORMAT(I2,I2,I2)
       IF(r2.gt.0)THEN
          WRITE(famrec(jj)%rescode(3:3),298)
 298      FORMAT('+')
       ENDIF
       IF(r3.gt.0)THEN
          WRITE(famrec(jj)%rescode(5:5),297)
 297      FORMAT('+')
       ENDIF
       CALL rmsp(famrec(jj)%rescode,le)
    ENDIF
 ENDDO
29 CALL filclo(in_res,' ')
 WRITE(*,*)' read threebody resonances, with proper elements are ',nr3
! -------------------------------------------------------------
! add resonance code for Trojans, Griquas and Hildas
! -------------------------------------------------------------
 DO jj=1,npro
   ap=propel(jj)%pr_el(1)
   IF(ap.gt.4.8d0.and.ap.lt.5.6d0)THEN
     famrec(jj)%rescode='1/1J  '
   ELSEIF(ap.gt.3.945d0.and.ap.lt.3.97d0)THEN
      famrec(jj)%rescode='3/2J  '
   ELSEIF(ap.gt.3.26d0.and.ap.lt.3.29d0)THEN
      famrec(jj)%rescode='2/1J  '
   ENDIF
 ENDDO
! -------------------------------------------------------------
! input list of mergers
! -------------------------------------------------------------
 CALL filopn(in_merg,'famdat/mergefam.out','old')
 READ(in_merg,*) 
 nmerg=0
 DO nf=1,nfamx
    READ(in_merg,*,END=2)merg_from, ist1,merge_to, ist2 
!
    nmerg=nmerg+1 
    merg(1,nmerg)=merg_from
    merg(2,nmerg)=merge_to
    istat_mer(1,nmerg)=ist1
    istat_mer(2,nmerg)=ist2
 ENDDO
2 CONTINUE
 CALL filclo(in_merg,' ')
 WRITE(iun_log,*)' merged families no. ',nmerg 
! -------------------------------------------------------------
!
! input list of family members
!
! -------------------------------------------------------------
 nfam=0
 nmemb=0
! zone 2
 CALL filopn(in_fam,'famdat/core_zona2_QRL.out','old')
 READ(in_fam,*)
 CALL input_families(in_fam,err_line)
 IF(err_line.ne.0)THEN
    WRITE(*,*)'core_zona2_QRL.out: err_line=',err_line
 ENDIF
 CALL filclo(in_fam,' ')
 nfam2=nfam
! zone 2 high I
 CALL filopn(in_fam,'famdat/highI_zone2_QRL.out','old')
 READ(in_fam,*)
 CALL input_families(in_fam,err_line)
 IF(err_line.ne.0)THEN
    WRITE(*,*)'highI_zone2_QRL.out: err_line=',err_line
 ENDIF
 CALL filclo(in_fam,' ')
 nfam2i=nfam
! zone 3
 CALL filopn(in_fam,'famdat/core_zona3_QRL.out','old')
 READ(in_fam,*)
 CALL input_families(in_fam,err_line)
 IF(err_line.ne.0)THEN
    WRITE(*,*)'core_zona3_QRL.out: err_line=',err_line
 ENDIF
 CALL filclo(in_fam,' ')
 nfam3=nfam
! zone 3 high I
 CALL filopn(in_fam,'famdat/highI_zone3_QRL.out','old')
 READ(in_fam,*)
 CALL input_families(in_fam,err_line)
 IF(err_line.ne.0)THEN
    WRITE(*,*)'highI_zone3_QRL.out: err_line=',err_line
 ENDIF
 CALL filclo(in_fam,' ')
 nfam3i=nfam
! zone 4
 CALL filopn(in_fam,'famdat/core_zona4_QRL.out','old')
 READ(in_fam,*)
 CALL input_families(in_fam,err_line)
 IF(err_line.ne.0)THEN
    WRITE(*,*)'core_zona4_QRL.out: err_line=',err_line
 ENDIF
 CALL filclo(in_fam,' ')
 nfam4=nfam
! zone 4 high I
 CALL filopn(in_fam,'famdat/highI_zone4_QRL.out','old')
 READ(in_fam,*)
 CALL input_families(in_fam,err_line)
 IF(err_line.ne.0)THEN
    WRITE(*,*)'highI_zone4_QRL.out: err_line=',err_line
 ENDIF
 CALL filclo(in_fam,' ')
 nfam4i=nfam
! zone 5, all inclinations
 CALL filopn(in_fam,'famdat/zona5_QRL.out','old')
 READ(in_fam,*)
 CALL input_families(in_fam,err_line)
 IF(err_line.ne.0)THEN
    WRITE(*,*)'zona5_QRL.out: err_line=',err_line
 ENDIF
 CALL filclo(in_fam,' ')
 nfam5=nfam
! zone Hungaria
 CALL filopn(in_fam,'famdat/zona1_QRL.out','old')
 READ(in_fam,*)
 CALL input_families(in_fam,err_line)
 IF(err_line.ne.0)THEN
    WRITE(*,*)'famdat/zona1_QRL.out: err_line=',err_line
 ENDIF
 CALL filclo(in_fam,' ')
 nfam1=nfam
! zone Hilda, new version 2015 
 CALL filopn(in_fam,'famdat/zona6_QRL_060.out','old')
  READ(in_fam,*)
  CALL input_families(in_fam,err_line)
 IF(err_line.ne.0)THEN
    WRITE(*,*)'famdat/zona6_QRL_060.out: err_line=',err_line
 ENDIF
 CALL filclo(in_fam,' ')
 nfam6=nfam
!---------------------------------------------------------
! modifications December 2015
! Trojans zone 
!-----------------------------------------------------------
 ! interpretation of numeric codes for multiopposition
 INQUIRE(file='famdat/conversiontroj.tab', exist=isthere)
 IF(isthere)THEN
! load conversion table
    CALL filopn(in_conver,'famdat/conversiontroj.tab','old')
    nconv=0
    DO i=1,ncodx
      READ(in_conver, *,END=17, ERR=16) nam, icode
      nconv=nconv+1
      namul(nconv)=nam
      icod(nconv)=icode
 16   CYCLE 
 17   EXIT
   ENDDO
   CALL filclo(in_conver,' ')  
 ELSE
    WRITE(*,*)' missing file of conversion for numeric multiopp codes'
    STOP 
 ENDIF
! actually read the family list; conversion is done in input_families
 CALL filopn(in_fam,'famdat/trojansL4.out','old')
 READ(in_fam,*)
 CALL input_families(in_fam,err_line)
 IF(err_line.ne.0)THEN
    WRITE(*,*)'trojansL4.out: err_line=',err_line
 ENDIF
 nfam7=nfam
 CALL filopn(in_fam,'famdat/trojansL5.out','old')
 READ(in_fam,*)
 CALL input_families(in_fam,err_line)
 IF(err_line.ne.0)THEN
    WRITE(*,*)'trojansL5.out: err_line=',err_line
 ENDIF
 nfam8=nfam
 CALL filopn(in_fam,'famdat/Griquas.out','old')
 READ(in_fam,*)
 CALL input_families(in_fam,err_line)
 IF(err_line.ne.0)THEN
    WRITE(*,*)'Griquas.out: err_line=',err_line
 ENDIF
 nfam9=nfam
! core families stop here; from now on, halo (or small) families with status=1
! -----------end modifications dec. 2015-----------------
! halo zone 2 
 CALL filopn(in_fam,'famdat/small_fam_zone2_040.out','old')
 READ(in_fam,*)
 CALL input_families(in_fam,err_line)
 IF(err_line.ne.0)THEN
    WRITE(*,*)'small_fam_zone2_040.out: err_line=',err_line
 ENDIF
 CALL filclo(in_fam,' ')
 nfamh2=nfam
 ! halo zone 3 
 CALL filopn(in_fam,'famdat/small_fam_zone3_040.out','old')
 READ(in_fam,*)
 CALL input_families(in_fam,err_line)
 IF(err_line.ne.0)THEN
    WRITE(*,*)'small_fam_zone3_040.out: err_line=',err_line
 ENDIF
 CALL filclo(in_fam,' ')
 nfamh3=nfam
! halo zone 4 
 CALL filopn(in_fam,'famdat/small_fam_zone4_040.out','old')
 READ(in_fam,*)
 CALL input_families(in_fam,err_line)
 IF(err_line.ne.0)THEN
    WRITE(*,*)'small_fam_zone4_040.out: err_line=',err_line
 ENDIF
 CALL filclo(in_fam,' ')
 nfamh4=nfam

 WRITE(iun_log,*)' input ', nfam, 'families with a total of ',nmemb,' members'
 WRITE(iun_log,*) ' core families ',nfam2,nfam2i,nfam3,nfam3i,nfam4,nfam4i,nfam5,nfam1,nfam6
 WRITE(iun_log,*) ' Trojans and Griquas families', nfam9
 WRITE(iun_log,*) ' small families ', nfamh2,nfamh3,nfamh4
 
! -------------------------------------------------------------------
! preparation to change family name to merged families
! -------------------------------------------------------------------
 DO nm=1,nmerg
    DO nf=1,nfam
       famname=memberlist(fampoint(nf))
       IF(famname.eq.merg(1,nm)) famno_merg(1,nm)=nf
       IF(famname.eq.merg(2,nm)) famno_merg(2,nm)=nf
    ENDDO
 ENDDO
! flag merged families
 DO nf=1,nfam
    famname=memberlist(fampoint(nf))
    merge_flag(nf)=0
    DO nm=1,nmerg
       IF(famname.eq.merg(1,nm))THEN
          IF(merge_flag(nf).ne.0)THEN
             WRITE(*,*)' attempt to merge family twice ', merg(:,nm)
             STOP
          ENDIF
          merge_flag(nf)=famno_merg(2,nm)
       ENDIF
    ENDDO
 ENDDO

! -------------------------------------------------------------
!
! loop on core families
!
! -------------------------------------------------------------
! -------------------------------------------------------------
! loop on known families; select QRL from Cellino et al. (2012) and from
! Novakovic et al. (2011), update 2013
 n1=0
 DO nf=1,nfam
    CALL zone_qrl(nf,kz,d_max)
    famname=memberlist(fampoint(nf))
    WRITE(*,*)famname, kz, d_max 
    n1f=n1
! -------------------------------------------------------------
!
! loop on established family members
!
! -------------------------------------------------------------
    DO j=1,ninfam(nf)
! loop searching for new members
       jjj=fampoint(nf)+j-1
       name=memberlist(jjj)
       CALL rmsp(name,le)
       CALL bin_search(name,npro,propel(1:npro)%name,isrtnam,jj)
       IF(jj.le.0.or.jj.gt.npro)THEN
          WRITE(iun_log,*)' name ', name, ' not found proper elements, jj=',jj
          CYCLE
       ENDIF
       pro=propel(jj)%pr_el
       name=propel(jj)%name
       hmag=propel(jj)%hmag
       CALL taxostep(jj,pro,name,hmag,d_max,n1)
! summary added members for family
       naddfam(nf)=n1-n1f 
! warning: list of additions could be empty
       clopoint(nf)=n1f+1  
  ! take care of merged families
       IF(merge_flag(nf).eq.0)THEN
          nf1=nf
          famname1=famname
       ELSE
          nf1=merge_flag(nf)
          famname1=memberlist(fampoint(nf1))
       ENDIF
! update family record for core/small members
       famrec(jj)%family(1)=famname1
       famrec(jj)%famno=nf1
       IF(nf.le.nfam9)THEN
! core families  keep their status even after merge, including trojans
          famrec(jj)%status=3
       ELSE
! halo/small families keep their status even after merge
          famrec(jj)%status=1
       ENDIF
    ENDDO
    WRITE(iun_log,*)' couples for family ',memberlist(fampoint(nf)), ' are ',n1-n1f
! find added members
    DO i=clopoint(nf),clopoint(nf)+naddfam(nf)-1
       i1=clos(i)%addr(1)
       i2=clos(i)%addr(2)
! update family record for additional members
       IF(famrec(i2)%status.eq.3.or.famrec(i2)%status.eq.1)THEN
! core members are kept as such
       ELSEIF(famrec(i2)%status.eq.2)THEN
! if already attributed to some family, is it the same? 
! taking into account of mergers, that is use famname1 and nf1
          IF(famrec(i2)%family(1).eq.famname1)THEN
! same family, different neighbour
             IF(famrec(i2)%dv(1).gt.clos(i)%sigma)THEN
! if the new is nearer, replace the old one
                famrec(i2)%dv(1)=clos(i)%sigma
                famrec(i2)%near(1)=clos(i)%names(1)
             ELSE
! if the old one is still nearer, leave it
             ENDIF
          ELSE
! this asteroid belongs to 2 families
             famrec(i2)%status=4
             IF(clos(i)%sigma.lt.famrec(i2)%dv(1))THEN
! the previous classification becomes secondary
                famrec(i2)%family(2)=famrec(i2)%family(1)
                famrec(i2)%dv(2)=famrec(i2)%dv(1)
                famrec(i2)%near(2)=famrec(i2)%near(1)
                famrec(i2)%famno(2)=famrec(i2)%famno(1)
! this becomes the "best" family classification
                famrec(i2)%family(1)=famname1
                famrec(i2)%dv(1)=clos(i)%sigma
                famrec(i2)%near(1)=clos(i)%names(1)
                famrec(i2)%famno(1)=nf1
             ELSE
! this becomes the secondary classification
                famrec(i2)%family(2)=famname
                famrec(i2)%famno(2)=nf1
                famrec(i2)%dv(2)=clos(i)%sigma
                famrec(i2)%near(2)=clos(i)%names(1)
             ENDIF
          ENDIF
       ELSEIF(famrec(i2)%status.eq.4)THEN
! this asteroid belongs to >1 families
! but this was known already INCOMPLETE
          IF(famrec(i2)%family(2).eq.famname1)THEN
! no new family
          ELSEIF(famrec(i2)%family(1).eq.famname1)THEN
! no new family
          ELSE
             famrec(i2)%status= famrec(i2)%status+1
          ENDIF
       ELSEIF(famrec(i2)%status.eq.0)THEN
! insert new status and new data
          famrec(i2)%status=2
          famrec(i2)%family(1)=famname1
          famrec(i2)%famno=nf1
          famrec(i2)%dv(1)=clos(i)%sigma
          famrec(i2)%near(1)=clos(i)%names(1)
       ENDIF         
    ENDDO
 ENDDO
 WRITE(iun_log,*)' end selection of close couples, total=',n1
 WRITE(*,*)' end selection of close couples, total=',n1
! ------------------------------------------------------------
!
! output the family records for all, including the zeros
!
! ------------------------------------------------------------
! initialize boxes before loop
 box(:,1,:)=100.d0
 box(:,2,:)=0.d0
! open output file
 CALL filopn(iun_rec,run(1:ler)//'.famrec','unknown')
! loop on all the asteroids with proper elements
 nadd=0
 ndouble=0
 nnewfam(:)=0
 nsmall(:)=0
 ncorefam(:)=0
 WRITE(iun_rec,500)
500 FORMAT('%ast.name',2X,'Hmag',3X,'status',1X,'family1 ',5X,'dv_fam1',&
&    1X,'near1',2X,' family2 ',4X,'dv_fam2', 1X,'near2',4X,'rescod') 
! loop on family classification records
 DO j=1,npr
    f=famrec(j)
    WRITE(iun_rec,400)f%name,f%hmag,f%status,f%family(1),f%dv(1),f%near(1), &
 &                               f%family(2),f%dv(2),f%near(2),f%rescode
400 FORMAT(A9,2X,F5.2,3X,I3,2X,2(2X,A9,2X,F7.2,2X,A9),2X,A6)
! ------------------------------------------------------------
! extended family summary
! -----------------------------------------------------------
    nf=f%famno(1)
    IF(nf.le.0) CYCLE
    CALL bin_search(f%name,npro,propel(1:npro)%name,isrtnam,jj)
    IF(jj.le.0.or.jj.gt.npro)THEN
       WRITE(*,*)' name ', name, ' not found in proper elements, jj=',jj
       CYCLE
    ENDIF
    pro=propel(jj)%pr_el
    DO i=1,3
       box(i,1,nf)=MIN(box(i,1,nf),pro(i))
       box(i,2,nf)=MAX(box(i,2,nf),pro(i))
    ENDDO
    IF(f%status.eq.2.or.f%status.eq.4)THEN 
          nnewfam(nf)=nnewfam(nf)+1
! ------------------------------------------------------------
       IF(f%status.eq.2)THEN
          nadd=nadd+1
       ELSEIF(f%status.eq.4)THEN
          ndouble=ndouble+1
       ENDIF
    ELSEIF(f%status.eq.1)THEN
       nsmall(nf)=nsmall(nf)+1
    ELSEIF(f%status.eq.3)THEN
       ncorefam(nf)=ncorefam(nf)+1
    ENDIF
 ENDDO
 CALL filclo(iun_rec,' ')
! -------------------------------------------------------------
!
! open output files: tables
!
! table with numbers in family, boxes, atc
!
! ------------------------------------------------------------
 CALL filopn(iun_tab,run(1:ler)//'.famtab','unknown')
 WRITE(iun_tab,302)'%core_fam zone d_max no_core no_small no_add no_tot  a_min  a_max  e_min  e_max sI_min sI_max'
302 FORMAT(A92)
!300 FORMAT(A9,2X,I2,2X,F5.0,2X,I5,2X,I5,2X,I5,2(2X,F5.3),2(2X,F5.3),2(2X,F5.3))
 WRITE(iun_tab,301)
301 FORMAT('%---------------------------------------------------------------------------------------------')
 DO nf=1,nfam
   IF(merge_flag(nf).ne.0) CYCLE
   CALL zone_qrl(nf,kz,d_max)
   famname=memberlist(fampoint(nf))
! write core family summary table
    WRITE(iun_tab,600)famname, kz, d_max, ncorefam(nf), nsmall(nf), nnewfam(nf), &
  &    ncorefam(nf)+nsmall(nf)+nnewfam(nf), (box(i,1,nf), box(i,2,nf), i=1,3)
600 FORMAT(A9,2X,I2,2X,F5.0,2X,I5,2X,I5,2X,I5,2X,I5,2(2X,F5.3),2(2X,F5.3),2(2X,F5.3))
 ENDDO
 CALL filclo(iun_tab,' ')
! ---------------------------------------------------------------------
! end message
! ---------------------------------------------------------------------
 WRITE(iun_log,*)' end selection of additional members, total=',nadd
 WRITE(iun_log,*)' multiple family members, total=',ndouble
 CALL filclo(iun_log,' ')
! ---------------------------------------------------------------------
CONTAINS
  SUBROUTINE zone_qrl(n,kz,d_max)
  INTEGER, INTENT(IN) :: n
  INTEGER, INTENT(OUT) :: kz
  DOUBLE PRECISION, INTENT(OUT) :: d_max
  IF(n.le.nfam2i)THEN
       kz=2
       IF(n.le.nfam2)THEN
! Cellino 2012
          d_max=70.d0
       ELSE
! Novakovic 2013
          d_max=130.d0
       ENDIF
    ELSEIF(n.le.nfam3i)THEN
       kz=3
       famname=memberlist(fampoint(n))
       IF(famname.eq.'5        ')THEN
          d_max=60.d0
       ELSEIF(n.le.nfam3)THEN
! Cellino 2012
          d_max=90.d0
       ELSE
! Novakovic 2013
          d_max=140.d0
       ENDIF
    ELSEIF(n.le.nfam4i)THEN
       kz=4
       IF(n.le.nfam4)THEN
! Cellino 2012
          d_max=100.d0
       ELSE
! Novakovic 2013
          d_max=80.d0
       ENDIF
    ELSEIF(n.le.nfam5)THEN
       kz=5
! Cellino 2012
       d_max=120.d0
    ELSEIF(n.le.nfam1)THEN
       kz=1
! Hungaria: Cellino 2013
       d_max=70.d0
    ELSEIF(n.le.nfam6)THEN
       kz=6
! Hilda: Cellino 2013-2015 mixed, use 2013
       d_max=60.d0
    ELSEIF(n.le.nfam7)THEN
       kz=7
! TrojansL4: Cellino 2015
       d_max=40.d0
    ELSEIF(n.le.nfam8)THEN
       kz=8 
! TrojansL5: Cellino 2015
       d_max=60.d0
    ELSEIF(n.le.nfam9)THEN
       kz=9
! Griqua: Cellino 2015
       d_max=150.d0
    ELSEIF(n.le.nfamh2)THEN
       kz=2
! Milani 2013
       d_max=40.d0
    ELSEIF(n.le.nfamh3)THEN
       kz=3
! Milani 2013
       d_max=40.d0
    ELSEIF(n.le.nfamh4)THEN
       kz=4
! Milani 2013
       d_max=40.d0
    ELSE
       WRITE(*,*)' inconsistent list of families nf=', n,nfam
       STOP
    ENDIF
  END SUBROUTINE zone_qrl

END PROGRAM prop_class

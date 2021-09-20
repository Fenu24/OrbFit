PROGRAM prop_tax
 USE propel_mod
 USE name_rules, ONLY: name_len
 USE fund_const
 USE orbit_elements
 USE dyn_param
 IMPLICIT NONE
 CHARACTER*6 progna ! option control: progname
 CHARACTER*80 run ! option control: current run 
 INTEGER ler ! character counters
 INTEGER in_pro ! input units
 INTEGER iun_log,nearout,neardat ! output units
 INTEGER err_line, npr
 INTEGER n1, j, i, i1, i2
 DOUBLE PRECISION d_max ! control
 DOUBLE PRECISION pro(3), hmag
 CHARACTER*(name_len) name
! arrays of osculating elements 
 INTEGER nosc, iii  ! total length, current position 
 INTEGER, PARAMETER :: noscx=1000000
 TYPE(orbit_elem), DIMENSION(noscx) :: oscels
 CHARACTER*(name_len), DIMENSION(noscx) :: namosc
 INTEGER, DIMENSION(noscx) :: indname
 CHARACTER(name_len) qnam
! output auxiliary data
 INTEGER kk
 DOUBLE PRECISION, DIMENSION(3):: dangles, dfreq, meaprop
 DOUBLE PRECISION princ, pridif
! **********BEGIN EXEXCUTION******************************** 
! progna='protax'
 WRITE(*,*)' run name?'
 READ(*,*) run 
 CALL rmsp(run,ler)
 WRITE(*,*)' control d_max?'
 READ(*,*) d_max 
 WRITE(*,*)' control ', d_max 
! open log file
 CALL filopn(iun_log,run(1:ler)//'.prta','unknown')
 WRITE(iun_log,*)' control ', d_max 
! open input file 
 CALL filopn(in_pro,'numb_res.syn','old')
! read proper elements
 npr=0
 CALL  input_propels(in_pro,npr,err_line)
 CALL filclo(in_pro,' ')
 WRITE(iun_log,*)' input ',npr,' proper elements of numbered'
! error return?
 IF(err_line.gt.0)THEN
    STOP
 ENDIF
! open input file NO: multiopposition do not have such accurate oscels
!  CALL filopn(in_pro,'mult.syn','old')
! read proper elements
!  CALL  input_propels(in_pro,npr,err_line)
!  WRITE(iun_log,*)' input ',npr,' prop. els of numbered and multiopp'
!  CALL filclo(in_pro,' ')
! read osculating elements
 CALL filopn(in_pro,'allnum.cat','old')
 CALL oporbf('allnum.cat',in_pro)
 CALL input_oscels(in_pro,nosc)
! CALL filclo(in_pro,' ')
 CALL clorbf
 WRITE(iun_log,*)' input ',nosc,' osculating elements of numbered'
! open output files 
 CALL filopn(nearout,run(1:ler)//'.near','unknown')
 WRITE(nearout,200)
200 FORMAT('% name   ',1X,'  H  ',1X,'  name   ',1X,'  H  ',1X,'  sigma   ',1X,&
&      '   da/a   ',1X,'    de    ',1X,'   dsinI   ')
!  FORMAT(A9,1X,F5.2,2X,A9,1X,F5.2,1X,F10.7,3(1X,F10.7))
 CALL filopn(neardat,run(1:ler)//'.ndat','unknown')
 WRITE(neardat,201)
 201 FORMAT('% name   ',1X,'  name   ',1X,'   pr.a(au) ',1X,' pr.e  ',1X,' pr. sinI'&
&      '   dn(deg/y) ',1X,' dg("/y) ',1X,' ds("/y) ',&
& 1X,'dlper(deg)',1X,'dnod(deg),'1X ' dm.an.(deg)')

! sort by sinI
 CALL heapsort(propel(1:npro)%pr_el(3),npro,isrti)
 WRITE(iun_log,*)' end sorting by sinI'
! main loop searching for neighbours
 n1=0
! select nearby cases
 DO j=1,npro-1
    pro=propel(j)%pr_el
    name=propel(j)%name
    hmag=propel(j)%hmag
    CALL taxostep(j,pro,name,hmag,d_max,n1)
    IF(mod(j,1000).eq.0) WRITE(iun_log,*)' tried ', j, ' objects, couples ',n1 
 ENDDO
 WRITE(iun_log,*)' end selection of close couples, n1=',n1
! sort by sigma
 CALL heapsort(clos(1:n1)%sigma,n1,isrts)
! write .near file
 DO j=1,n1
    i=isrts(j)
    i1=clos(i)%addr(1)
    i2=clos(i)%addr(2)
! write close couple with magnitudes, distance in proper elements
    WRITE(nearout,100)clos(i)%names(1),propel(i1)%hmag,clos(i)%names(2),propel(i2)%hmag,&
&           clos(i)%sigma,clos(i)%delta_el(1),clos(i)%delta_el(2),clos(i)%delta_el(3)
100 FORMAT(A9,1X,F5.2,2X,A9,1X,F5.2,1X,F10.7,3(1X,F10.7))
! auxiliary data include osculating elements data
    CALL  bin_search(clos(i)%names(1),nosc,namosc,indname,iii)
    IF(iii.le.0)THEN
       WRITE(*,*)' missing asteroid ',clos(i)%names(1),' in osculating elements catalog' 
    ELSE
       clos(i)%angles(2:3,1)=oscels(iii)%coord(5:6)
! use longitude of pericenter
       clos(i)%angles(1,1)=princ(oscels(iii)%coord(4)+oscels(iii)%coord(5))
    ENDIF
    CALL  bin_search(clos(i)%names(2),nosc,namosc,indname,iii)
    IF(iii.le.0)THEN
       WRITE(*,*)' missing asteroid ',clos(i)%names(2),' in osculating elements catalog' 
    ELSE
       clos(i)%angles(2:3,2)=oscels(iii)%coord(5:6)
! use longitude of pericenter
       clos(i)%angles(1,2)=princ(oscels(iii)%coord(4)+oscels(iii)%coord(5))
    ENDIF
! differences in delaunay angle variables in degrees, mean of proper elements, diff in freq (same units as in file)
    dangles=0.d0;meaprop=0.d0;dfreq=0.d0
    DO kk=1,3
! differences in delaunay angle variables in degrees
      dangles(kk)=degrad*pridif(clos(i)%angles(kk,1),clos(i)%angles(kk,2))
! mean of proper elements
      meaprop(kk)=0.5d0*(propel(i1)%pr_el(kk)+propel(i2)%pr_el(kk))
! diff in freq (same units as in file)
      dfreq(kk)=propel(i1)%pr_fr(kk)-propel(i2)%pr_fr(kk)
    ENDDO
! write names, mean of proper elements, differences in proper frequencies, differences in Delaunay angles
    WRITE(neardat,101)clos(i)%names(1),clos(i)%names(2),meaprop,dfreq,dangles
101 FORMAT(A9,2X,A9,1X,3(1X,F9.7),3(1X,F12.6),3(1X,F11.6))
 ENDDO
 CALL filclo(nearout,' ')
 CALL filclo(neardat,' ')
! summary
 WRITE(*,*)' found ', n1,' couples closer than ', d_max
 WRITE(iun_log,*)' found ', n1,' couples closer than ', d_max

CONTAINS

! read catalog of osculating elements, assumed all at same epoch, and store
! in array, with indexing by name
   SUBROUTINE input_oscels(iun,n)
    USE dyn_param, ONLY: ndyx
    INTEGER, INTENT(IN) :: iun   ! input unit
    INTEGER, INTENT(OUT) :: n    ! number actually read
    INTEGER i, nlsloc, lsloc(ndyx)
    LOGICAL eof, err
    DO i=1,noscx
       CALL read_elems(oscels(i),namosc(i),eof,nlsloc,err,lsloc,FILE='    ',UNIT=iun)
       IF(eof) EXIT
    ENDDO
    n=i-1
    IF(n.eq.noscx)THEN
       WRITE(*,*)' input:oscels: name array too long, max is  ', noscx
    ENDIF
    CALL heapsortnach(namosc,n,indname)
    WRITE(*,*)' read ',n,' osculating elements'
   END SUBROUTINE input_oscels

END PROGRAM prop_tax

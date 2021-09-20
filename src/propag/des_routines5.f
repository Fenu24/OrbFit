c ================================================================
c
c   LIBRARY DES_ROUTINES
c
c FORTRAN 77 subroutines for I/O of data according
c to the Pan-STARRS data exchange standard
c see PSDC-xxx for definitions.
c
c ================================================================
c
c EPHEMERIS
c
c ================================================================
c
c WRITE_EPHEMERIS
c
c ================================================================
      SUBROUTINE write_ephemeris(iun,id_oid,index,time,att4,g4,appma,
     + rmsmag)
      IMPLICIT NONE
c INPUT: unit (already opened)
      INTEGER iun
c INPUT: identification/orbit OID
      CHARACTER*(*) id_oid
c INPUT: index of solutions
      INTEGER index   
c INPUT: epoch time (MJD), apparent magnitude
      DOUBLE PRECISION time, appma
c INPUT: attributable, (assumed all in DEG, DEG/day)
      DOUBLE PRECISION att4(4)
c covariance matrix (assumed in arcsec, arcsec/day)
      DOUBLE PRECISION  g4(4,4) 
c INPUT: photometric uncertainty 
      DOUBLE PRECISION rmsmag
c =================================================================
      DOUBLE PRECISION rms4(4), cr(4,4)
      INTEGER i,j
      DO i=1,4
        rms4(i)=sqrt(g4(i,i))
      ENDDO
      DO i=1,4
        DO j=1,i-1
          cr(i,j)=g4(i,j)/(rms4(i)*rms4(j))          
        ENDDO
      ENDDO
      WRITE(iun,100)id_oid, index, time, att4, appma, rms4, 
     + cr(2,1),cr(3,1),cr(4,1),cr(3,2),cr(4,2),cr(4,3),rmsmag 
 100  FORMAT(A,1X,I3,1X,F14.7,1X,2(F16.12,1X),2(F15.10,1X),F13.10,1X,
     + 1P,2(D11.4,1X),2(D11.4,1X),0P,6(F9.5,1X),F7.4)
      RETURN
      END
c ================================================================
c
c EPHEM_HEADER
c
c ================================================================
      SUBROUTINE ephem_header(iun)
      IMPLICIT NONE
c INPUT: unit (already opened)
      INTEGER iun
c column headings
      WRITE(iun,100)
 100  FORMAT('!!ident_oid index time      ',
     + '     RA(deg)       ','    DEC(deg)      ','   RAdot(deg/d)   ',
     + '   DECdot(deg/d)     ','   APP.MAG   ','    RMS(RA)   ',
     + ' RMS(DEC)    ',' RMS(RAdot)   ',' RMS(DECdot)   ',
     + '  CORRELS  ',' RMS(App.MAG) ')
      RETURN
      END
c ================================================================
c
c METRICS
c
c ================================================================
c
c READ_METRICS
c
c ================================================================
      SUBROUTINE read_metrics(iun,id_oid,index,npar,nrej,n_met,met,eof)
      IMPLICIT NONE
c INPUT: unit (already opened)
      INTEGER iun
c OUPUT: identification/orbit OID
      CHARACTER*(*) id_oid
c OUTPUT: index of solutions, of fit parameter
      INTEGER index,npar,nrej
c OUTPUT: metrics
      INTEGER n_met
      DOUBLE PRECISION met(n_met)
c OUTPUT: end of file 
      LOGICAL eof
c =================================================================
      CHARACTER*512 record
      eof=.false.
 1    CONTINUE
      READ(iun,'(A)',END=2) record
      IF(record(1:1).eq.'!')THEN
         IF(record(2:2).eq.'!')THEN
!            WRITE(*,*)'read_metrics: COLUMN HEADINGS'
c len_trim or rmsp  would be needed...
!            WRITE(*,*) record
         ENDIF
         GOTO 1
      ENDIF
      READ(record,*,ERR=21,END=21) id_oid, index, npar, nrej, met
      RETURN
 2    eof=.true.
c note that on eof the elements are undefined in exit
      RETURN
 21   WRITE(*,*)' read_metrics: error in reading file at unit ', iun
      WRITE(*,*) id_oid, index, npar, nrej, met
      END
c ================================================================
c
c WRITE_METRICS
c
c ================================================================
      SUBROUTINE write_metrics(iun,id_oid,index,npar,nrej,n_met,met)
      IMPLICIT NONE
c INPUT: unit (already opened)
      INTEGER iun
c INPUT: identification/orbit OID
      CHARACTER*(*) id_oid
c INPUT: index of solutions, of fit parameter
      INTEGER index,npar,nrej
c INPUT: metrics
      INTEGER n_met
      DOUBLE PRECISION met(n_met)
c =================================================================
      INTEGER le
      CALL remspace(id_oid,le)
      WRITE(iun,100,ERR=20) id_oid(1:le), index, npar, nrej, met
 100  FORMAT(A,1X,I3,1X,I2,1X,I4,(2X,3(1X,F9.5)),2(2X,4(1x,F9.5))) 
      RETURN
 20   WRITE(*,*)' write_metrics: error in output format on unit ',iun
      WRITE(*,*)  id_oid, index, npar, nrej, met
c      STOP
      END
c ================================================================
c
c METRICS_HEADER
c
c ================================================================
      SUBROUTINE metrics_header(iun)
      IMPLICIT NONE
c INPUT: unit (already opened)
      INTEGER iun
c column headings
      WRITE(iun,100)
 100  FORMAT('!!ident_oid index npar nrej ',
     + '   RMS   ',' RMSmag  ','Max% rej ',
     + 'RA Res.Zsi.',' Res.Curv',' Res.trend',' Res.bias',
     + ' DECRes.Zsi.',' Res.Curv',' Res.trend',' Res.bias')
      RETURN
      END
c ================================================================
c
c COVARIANCE
c
c ================================================================
c
c READ_COVAR(iance and normal matrices)
c
c ================================================================
      SUBROUTINE read_covar(iun, cove, nore, orb_oid, ind, eof)
      IMPLICIT NONE
c INPUT: unit (already opened)
      INTEGER iun
c OUTPUT: covariance matrix, normal matrix, symmetric, for COM elems
c        with angles I, Omega, argperi in degrees
      DOUBLE PRECISION cove(6,6), nore(6,6)
c OUTPUT: unique identifier, can have many characters, and index 
      CHARACTER*(*) orb_oid
      INTEGER ind
c OUTPUT: end of file
      LOGICAL eof
c ================================================================
      INTEGER le, i, k, jrec,ii
c ================================================================
      CHARACTER*512 record
      CHARACTER* 512 oid
      eof=.false.
 1    CONTINUE
      READ(iun,'(A)',END=2) record
      IF(record(1:1).eq.'!')THEN
         IF(record(2:2).eq.'!')THEN
!            WRITE(*,*)'read_covar: COLUMN HEADINGS'
c len_trim or rmsp  would be needed...
!            WRITE(*,*) record
         ENDIF
         GOTO 1
      ENDIF
      READ(record,'(A)') oid
      ii=index(oid,' ')
      orb_oid=' '
      orb_oid=oid(1:ii-1)
      READ(record(ii+1:),*) ind
      IF(ind.eq.0) ind=1
      jrec=0
c limitation: no comments after the orb_oid line
      READ(iun,108, ERR=4) ((cove(i,k),k=i,6),i=1,6) 
 108  FORMAT(5X,3E23.15)
      READ(iun,109, END=4) ((nore(i,k),k=i,6),i=1,6) 
 109  FORMAT(5X,3E23.15)
c should be all read      
      DO i=1,6
        DO k=1,i-1
          cove(i,k)=cove(k,i)
          nore(i,k)=nore(k,i)
        ENDDO
      ENDDO
      RETURN
 2    CONTINUE
      eof=.true.
      RETURN
 4    CONTINUE
      WRITE(*,*)' read_covar: error in reading format'
      WRITE(*,*) record
      STOP 
      END
c ================================================================
c
c WRITE_COVAR(iance and normal matrices)
c
c ================================================================
      SUBROUTINE write_covar(iun, cove, nore, orb_oid, index)
      IMPLICIT NONE
c INPUT: unit (already opened)
      INTEGER iun
c INPUT: covariance matrix, normal matrix, symmetric, for COM elems
c        with angles I, Omega, argperi in degrees
      DOUBLE PRECISION cove(6,6), nore(6,6)
c INPUT: unique identifier, can have many characters
      CHARACTER*(*) orb_oid
c INPUT: index of multiple solution; defaults to 1
      INTEGER index
c ================================================================
      INTEGER le, i, k
      CALL remspace(orb_oid,le)
      WRITE(iun,'(A,2X,I4)') orb_oid(1:le),index
      WRITE(iun,108) ((cove(i,k),k=i,6),i=1,6) 
108   FORMAT(' COV ',1P,3E23.15)
      WRITE(iun,109) ((nore(i,k),k=i,6),i=1,6) 
109   FORMAT(' NOR ',1P,3E23.15)
      RETURN
      END
c ================================================================
c
c COVAR HEADER
c
c ================================================================
      SUBROUTINE covar_header(iun)
      IMPLICIT NONE
c INPUT: unit (already opened)
      INTEGER iun
c column headings
      WRITE(iun,100)
 100  FORMAT('!! orbit_oid and covariance/normal matrices, 6x6, '
     + ,'symmetric, lower triangular,')
      WRITE(iun,101)
 101  FORMAT('!! with column index growing first; elements '
     +  ,'cometary, angles in degrees')
      RETURN
      END
c ================================================================
c
c ORBIT
c
c ================================================================
c  WARNING: what about G (opposition effect)? now not stored,
c           assumed at default value 0.15
c ================================================================
c
c READ_ORBIT
c
c ================================================================
      SUBROUTINE read_orbit(iun, coo, ele, mag, epoch, orb_oid, eof,
     +    index,n_par,moid, orbcc)
      IMPLICIT NONE
c INPUT: unit (already opened)
      INTEGER iun
c OUTPUT: coordinates found in input
      CHARACTER*3 coo
c OUTPUT: elements (6 coords), absolute magnitude, epoch time (MJD TDT)
      DOUBLE PRECISION ele(6), mag, epoch
c OUTPUT: unique identifier, can have many characters
      CHARACTER*(*) orb_oid
c OUTPUT: end of file
      LOGICAL eof
c OUTPUT: : index of multiple solution, no. parameters (4,5,6)
      INTEGER index,n_par
c OUTPUT: MOID 
      DOUBLE PRECISION moid 
c OUTPUT: orbit computer code
      CHARACTER*(*) orbcc
c ================================================================
      CHARACTER*512 record
      eof=.false.
 1    CONTINUE
      READ(iun,'(A)',END=2) record
      IF(record(1:1).eq.'!')THEN
         IF(record(2:2).eq.'!')THEN
!            WRITE(*,*)'read_orbit: COLUMN HEADINGS'
c len_trim or rmsp  would be needed...
!            WRITE(*,*) record
         ENDIF
         GOTO 1
      ENDIF
      READ(record,*,ERR=21)orb_oid,coo,ele,mag,epoch, 
     +            index,n_par,moid,orbcc
      RETURN
 2    eof=.true.
c note that on eof the elements are undefined in exit
      RETURN
 21   WRITE(*,*)' read_orbit: error in reading file at unit ', iun
      WRITE(*,*) coo, ele, mag, epoch, orb_oid, n_par, moid, orbcc
      WRITE(*,*) record
      END
c ================================================================
c
c WRITE_ORBIT
c
c ================================================================
      SUBROUTINE write_orbit(iun,coo,ele,mag,epoch,orb_oid,   
     +       index,n_par,moid, orbcc)
      IMPLICIT NONE
c INPUT: unit (already opened)
      INTEGER iun
c INPUT: type of elements, elements (6 coords), absolute magnitude, epoch time (MJD TDT)
      CHARACTER*3 coo
      DOUBLE PRECISION ele(6), mag, epoch
c INPUT: unique identifier, can have many characters
      CHARACTER*(*) orb_oid
c INPUT: index of multiple solution, no. parameters (4,5,6)
      INTEGER index,n_par
c INPUT: MOID 
      DOUBLE PRECISION moid 
c INPUT: orbit computer code
      CHARACTER*(*) orbcc
c ================================================================
      INTEGER le
      DOUBLE PRECISION moidsafe
      CALL remspace(orb_oid,le)
      moidsafe=MIN(moid,1.d6)
      IF(coo.eq.'COM')THEN
         WRITE(iun,100,ERR=20) orb_oid(1:le),coo,ele, mag, epoch,
     +        index,n_par,moidsafe,orbcc
 100     FORMAT(A,1X,A3,1X,F16.10,1X,F14.10,1X,F11.7,1X,F11.7,1X,F11.7,
     +       1X,F16.7,1X,F6.3,1X,F14.7,1X,I4,1X,I2,1X,F14.5,1X,A)
c     ele(1)=q (perihelion distance, in AU)
c     ele(2)=e (eccentricity)
c     ele(3)=I (inclination on ecliptic, in deg)
c     ele(4)=Omega (longitude of ascending node, in deg)
c     ele(5)=argperi (angle between node and direction of pericenter, in deg)
c     ele(6)=t_peri (time of nearest perihelion, MJD TDT)
      ELSEIF(coo.eq.'COT')THEN
         WRITE(iun,102,ERR=20) orb_oid(1:le), coo, ele, mag, epoch, 
     +        index,n_par,moidsafe,orbcc
 102     FORMAT(A,1X,A3,1X,F16.10,1X,F14.10,4(1X,F12.8),
     +       1X,F6.3,1X,F14.7,1X,I4,1X,I2,1X,F14.5,1X,A)
c     as above, but
c     ele(6)=true_an (true anomaly, in deg)
      ELSEIF(coo.eq.'KEP')THEN
         WRITE(iun,101,ERR=20) orb_oid(1:le), coo, ele, mag, epoch,
     +        index,n_par,moidsafe,orbcc
 101     FORMAT(A,1X,A3,1X,F16.10,1X,F12.10,4(1X,F12.8),
     +       1X,F6.3,1X,F14.7,1X,I4,1X,I2,1X,F14.5,1X,A)
c     ele(1)=a (semimajor axis, in AU)
c     ele(2)=e (eccentricity)
c     ele(3)=I (inclination on ecliptic, in deg)
c     ele(4)=Omega (longitude of ascending node, in deg)
c     ele(5)=argperi (angle between node and direction of pericenter, in deg)
c     ele(6)=mean_an (mean anomaly, in deg)
       ELSEIF(coo.eq.'EQU')THEN
         WRITE(iun,105,ERR=20) orb_oid(1:le), coo, ele, mag, epoch,
     +        index,n_par,moidsafe,orbcc
 105     FORMAT(A,1X,A3,1X,F16.10,2(1X,F12.10),2(1X,F14.10),1X,F12.8,
     +       1X,F6.3,1X,F14.7,1X,I4,1X,I2,1X,F14.5,1X,A)
c     ele(1)=a (semimajor axis, in AU)
c     ele(2)= h = e sin(varpi) ; varpi = Omega + omega
c     ele(3)= k = e cos(varpi)
c     ele(4)= p = tg(I/2) sin(Omega)
c     ele(5)= q = tg(I/2) cos(Omega)
c     ele(6)= mean longitude, measured from fb vector (mean an + varpi)
      ELSEIF(coo.eq.'CAR')THEN
         WRITE(iun,103,ERR=20) orb_oid(1:le), coo, ele, mag, epoch,
     +        index,n_par,moidsafe,orbcc
 103     FORMAT(A,1X,A3,1X,6(F16.10,1X),
     +       F6.3,1X,F14.7,1X,I4,1X,I2,1X,F14.5,1X,A)
c     ele(1:3)= cartesian positions in AU, heliocentric
c     ele(4:6)= cartesian velocity in AU/day, heliocentric
      ELSEIF(coo.eq.'ATT')THEN
         WRITE(iun,104,ERR=20) orb_oid(1:le), coo, ele, mag, epoch,
     +        index,n_par,moidsafe,orbcc
 104     FORMAT(A,1X,A3,1X,6(F16.10,1X),
     +       F6.3,1X,F14.7,1X,I4,1X,I2,1X,F14.5,1X,A)
      ELSE
         WRITE(*,*)' write_orbit: error in coordinate type ', coo
         STOP
      ENDIF
      RETURN
 20   WRITE(*,*)' write_orbit_opt: error in output format on unit ',iun
      WRITE(*,*) ele, mag, epoch, orb_oid,orbcc
      STOP
      END
c ================================================================
c
c ORBIT HEADER
c
c ================================================================
      SUBROUTINE orbit_header(iun,cooy)
      IMPLICIT NONE
c INPUT: unit (already opened)
      INTEGER iun
c INPUT: coordinates
      CHARACTER*3 cooy
c column headings
      IF(cooy.eq.'COM')THEN
         WRITE(iun,100)
 100     FORMAT('!!   oid     q(AU)       e     I(deg)     Omega(deg)'
     + ,' argperi(deg) t_peri(MJD)  absmag  '
     + ,'epoch(MJD)    index n_par moid')
      ELSEIF(cooy.eq.'COT')THEN
         WRITE(iun,102)
 102     FORMAT('!!  oid      q(AU)       e     I(deg)     Omega(deg)'
     + ,' argperi(deg) true_an(deg)  absmag  '
     + ,'epoch(MJD)    oid index n_par moid')
      ELSEIF(cooy.eq.'KEP')THEN
         WRITE(iun,101)
 101     FORMAT('!!  oid      a(AU)       e      I(deg)     Omega(deg)'
     + ,' argperi(deg) mean_an(deg)  absmag  '
     + ,'epoch(MJD)    oid index n_par moid')
      ELSEIF(cooy.eq.'CAR')THEN
         WRITE(iun,103)
 103     FORMAT('!!  oid      x(AU)     y(AU)     z(AU)   xdot(AU/d) '
     + ,' ydot(AU/d)   zdot(AU/d)  absmag  '
     + ,'epoch(MJD)    index n_par moid')
      ELSEIF(cooy.eq.'ATT')THEN
         WRITE(iun,104)
 104     FORMAT('!!  oid    RA(deg)    DEC(deg)    RAdot(deg/d)   ',
     + '  DECdot(deg/d)   rho(AU)   rhodot(AU/d)  absmag  '
     + ,'epoch(MJD)    index n_par moid')

      ELSE
         WRITE(*,*)' orbit_header: error in coordinate type ', cooy
         STOP
      ENDIF
      RETURN
      END
c ================================================================
c
c IDENT_HEADER
c
c ================================================================
c
c READ_IDHEAD
c
c ================================================================
      SUBROUTINE read_idhead(iunids,eof,oid,nidx,nid,tranam,op_code, 
     +  counters, param)
      IMPLICIT NONE
c INPUT: unit number (must be opened already) 
      INTEGER iunids
c OUTPUT: end-of-file, unique identifier
      LOGICAL eof 
      CHARACTER*(*) oid
c INPUT: max no of ids
      INTEGER nidx
c OUTPUT: number of tracklets, list of their oid
      INTEGER nid
      CHARACTER*(*) tranam(nidx)
c OUTPUT: op code (what to do)
      CHARACTER*(*) op_code
c OUTPUT: counters: N_OBS, N_SOLUTIONS,  N_NIGHTS,  ARC_TYPE
      INTEGER counters(5)
c OUTPUT: parameters (for requests): CHI_TRACK, S2N_TRACK, CHI_SINGLE, S2N_SINGLE
      DOUBLE PRECISION param(4)
c ===================================================================
c END INTERFACE
      CHARACTER*512 record
      INTEGER j
c ==================================================================
      eof=.false.
c read one record
 1    CONTINUE
      READ(iunids,'(A)',END=2) record
      IF(record(1:1).eq.'!')THEN
         IF(record(2:2).eq.'!')THEN
!            WRITE(*,*)'read_idhead: COLUMN HEADINGS'
c len_trim or rmsp  would be needed...
!            WRITE(*,*) record
         ENDIF
         GOTO 1
      ENDIF
      READ(record,*,ERR=3)oid,nid
      IF(nid.gt.nidx)THEN
         WRITE(*,*) ' read_idhead: ', nid, ' tracklets, max was ', nidx
         STOP '****** error in reading ident_header*****'
      ENDIF
      READ(record,*,ERR=3)oid,nid,(tranam(j),j=1,nid),op_code,counters
      DO j=nid+1,nidx
         tranam(j)=' '
      ENDDO
c WARNING: this tolerance for missing fields in the record
c is wrong, and should be removed. Actually, the first param has
c meaning (rms_prelim) also for REQUEST_PRELIM.
      param(1)=0.d0
      param(2)=0.d0
      param(3)=0.d0
      param(4)=0.d0
c      IF(op_code.ne.'REQUEST_PRELIM'.and.
c     +      op_code.ne.'REQUEST_ORBIT')THEN
         READ(record,*)oid,nid,(tranam(j),j=1,nid),op_code, 
     +      counters,param
c      ENDIF
      RETURN
 3    WRITE(*,*)' read_idhead: error in reading ', record 
      STOP '****** error in reading ident_header*****'
 2    CONTINUE
      eof=.true.
      RETURN
      END
c ================================================================
c
c WRITE_IDHEAD
c
c ================================================================
      SUBROUTINE write_idhead(iunids,oid,nid,tranam,op_code, 
     +  counters, parameters)
      IMPLICIT NONE
c INPUT: unit number (must be opened already) 
      INTEGER iunids
c INPUT: unique identifier
      CHARACTER*(*) oid
c INPUT: number of tracklets, list of their oid
      INTEGER nid
      CHARACTER*(9) tranam(nid)
c INPUT: op code (what to do)
      CHARACTER*(*) op_code
c INPUT: counters: N_OBS, N_SOLUTIONS,  N_NIGHTS,  ARC_TYPE
      INTEGER counters(5)
c INPUT: parameters (for requests): CHI_TRACK, S2N_TRACK, CHI_SINGLE, S2N_SINGLE
c if the op_code is not REQUEST_PRECOVERY, not written
      DOUBLE PRECISION parameters(4)
c ===================================================================
c END INTERFACE
      CHARACTER*512 record, record1, record2
      CHARACTER*512 recout
      INTEGER j,k,k1,k2,kk
      CHARACTER*200 form
c ==================================================================
c write one record
c      IF(nid.gt.16)THEN
c         WRITE(*,*) ' write_idhead: ', nid, ' tracklets, max was ', 16
c         STOP '****** error in writinging ident_header*****'
c      ENDIF
      WRITE(record,'(A)')oid
      CALL remspace(record,k)
      WRITE(form,120)nid
 120  FORMAT('(1X,I3,1X,',I2,'(A9,1X))')
      WRITE(record1,form,ERR=3)nid,(tranam(j),j=1,nid)
 100  FORMAT(1X,I3,1X,18(A9,1X))
      k1=5+10*nid
      IF(parameters(1).eq.0.d0.and.parameters(2).eq.0.d0.and.
     +     parameters(3).eq.0.d0.and.parameters(4).eq.0.d0)THEN
         WRITE(record2,201,ERR=3)op_code,counters,0,0,0,0
 201     FORMAT(1X,A,1X,I4,1X,I3,1X,I4,2(1X,I3),4(1X,I1))
         k2=21+16+20
      ELSE
         WRITE(record2,101,ERR=3)op_code,counters,parameters
 101     FORMAT(1X,A,1X,I4,1X,I3,1X,I4,2(1X,I3),4(1X,F10.4))
         k2=21+20+44
      ENDIF
      recout=record(1:k)//record1(1:k1)//record2(1:k2)
      kk=k+k1+k2
      WRITE(iunids,'(A)') recout(1:kk)
      RETURN
 3    WRITE(*,*)' write_idhead: error in writing ', record 
      STOP '****** error in writing ident_header*****'
      END
c ================================================================
c
c IDHEAD_HEADER
c
c ================================================================
      SUBROUTINE idhead_header(iunids)
      IMPLICIT NONE
c INPUT: unit (already opened)
      INTEGER iunids
c column headings
      WRITE(iunids,100)
 100  FORMAT('!!  id_oid   nid  tracklet_oids ',
     + ' OP_CODE  N_OBS N_SOL N_NIGHTS ARC_TYPE NO_RAD PARAM')
      RETURN
      END

c ================================================================
c
c DETECTIONS, TRACKLETS AND RESIDUALS
c
c ================================================================
c 
c READ_RESIDUAL
c
c ================================================================
      SUBROUTINE read_residual(iunobs,eof,oid,tutc, resra,resdec,
     +    resmag,observatory,wei_ra,wei_dec,wei_mag,
     +    chi,sel_astr,sel_phot,secnam)
      IMPLICIT NONE
c INPUT: unit number (must be opened already) 
      INTEGER iunobs
c OUTPUT: end-of-file, unique identifier
      LOGICAL eof 
      CHARACTER*(*) oid
c OUTPUT: obscode
      CHARACTER*3 observatory
c OUTPUT: time, UTC MJD
      DOUBLE PRECISION tutc
c OUTPUT: residual in R.A.(times cos(DEC)), in DEC (arcsec), in apparent magnitude (mag)
      DOUBLE PRECISION resra,resdec,resmag
c OUTPUT: WEIGHT of RA*cos(DEC), of DEC, both in arcsec, WEIGHT of APPMAG (mag)
      DOUBLE PRECISION wei_ra, wei_dec, wei_mag
c OUTPUT: chi value from outlier rejection, astrometry used in fit (if >0) , 
c           photometry used in fit (if >0)
      DOUBLE PRECISION chi
      INTEGER sel_astr,sel_phot
c INPUT: secret name in simulations/true name for known objects, NA otherwise
      CHARACTER*(*) secnam
c END INTERFACE
      CHARACTER*512 record
c ==================================================================
      eof=.false.
c read one record
 1    CONTINUE
      READ(iunobs,'(A)',END=2) record
      IF(record(1:1).eq.'!')THEN
         IF(record(2:2).eq.'!')THEN
!            WRITE(*,*)'read_residual: COLUMN HEADINGS'
c len_trim or rmsp  would be needed...
!            WRITE(*,*) record
         ENDIF
         GOTO 1
      ENDIF
      READ(record,*,ERR=3) oid,tutc,observatory,resra,resdec,resmag,
     +     wei_ra,wei_dec,wei_mag,chi,sel_astr,sel_phot,secnam
      RETURN
 2    CONTINUE
      eof=.true.
      RETURN
c WARNING: the data are undefined in output if eof
 3    CONTINUE
      WRITE(*,*) ' read-residual: error in reading the record'
      WRITE(*,*) record
      STOP ' ***** error in reading detections****'
      END
c ================================================================
c 
c WRITE_RESIDUAL
c
c ================================================================
      SUBROUTINE write_residual(iunobs,oid,tutc, resra,resdec,
     +    resmag,observatory,wei_ra,wei_dec,wei_mag,
     +    chi,sel_astr,sel_phot,secnam)
      IMPLICIT NONE
c INPUT: unit number (must be opened already) 
      INTEGER iunobs
c     INPUT: unique identifier of the trackletdata
      CHARACTER*(*) oid
c INPUT: obscode
      CHARACTER*3 observatory
c INPUT: time, UTC MJD
      DOUBLE PRECISION tutc
c INPUT: residual in R.A.(times cos(DEC)), in DEC (arcsec), 
c       in apparent magnitude (mag)
      DOUBLE PRECISION resra,resdec,resmag
c INPUT: WEIGHT of RA*cos(DEC), of DEC, in arcsec, WEIGHT APPMAG (mag)
      DOUBLE PRECISION wei_ra, wei_dec, wei_mag
c INPUT: chi value from outlier rejection, astrometry used in fit (if >0) , 
c           photometry used in fit (if >0)
      DOUBLE PRECISION chi
      INTEGER sel_astr,sel_phot
c INPUT: secret name in simulations/true name for known objects, NA otherwise
      CHARACTER*(*) secnam
c ==================================================================
c write one record
      IF(abs(resra).lt.9.d3.and.abs(resdec).lt.9.d3.and.abs(chi).
     +        lt.9.d2)THEN
         WRITE(iunobs,100) oid,tutc,observatory,resra,resdec,resmag,
     +     wei_ra,wei_dec,wei_mag,chi,sel_astr,sel_phot,secnam
 100  FORMAT(A,1x,f17.10,1x,A3,1x,f10.4,1x,f10.4,1x,f7.3,1x,
     +       f9.4,1x,f9.4,1x,f7.3,1x,f9.4,1x,I2,1x,I2,1x,A)
      ELSE
         WRITE(iunobs,101) oid,tutc,observatory,resra,resdec,resmag,
     +     wei_ra,wei_dec,wei_mag,chi,sel_astr,sel_phot,secnam
 101     FORMAT(A,1x,f17.10,1x,A3,1x,1P,D10.3,1x,D10.3,0P,1x,f7.3,1x,
     +       f9.4,1x,f9.4,1x,f7.3,1x,1P,D9.2,0P,1x,I2,1x,I2,1x,A)
      ENDIF
c later we could implement a smart format, with the appropriate number 
c of digits computed on the basis of the weight.
      RETURN
      END
c ================================================================
c
c RESIDUAL HEADER
c
c ================================================================
      SUBROUTINE residual_header(iunres)
      IMPLICIT NONE
c INPUT: unit (already opened)
      INTEGER iunres
c column headings
      WRITE(iunres,100)
 100  FORMAT('!! oid    time (UTC,MJD)  obscod  Res. RA    Res. DEC'
     +    ,' Res.Mag  Wei. RA  Wei. DEC  Wei.Mag '
     +    ,' chi sel.A. sel.P. name')
      RETURN
      END
c ================================================================
c 
c WRITE_RADAR_RESIDUAL
c
c ================================================================
      SUBROUTINE write_radar_residual(iunobs,oid,tutc,observ, type, 
     +    res,weight,chi,selflag)
      IMPLICIT NONE
c INPUT: unit number (must be opened already) 
      INTEGER iunobs
c     INPUT: unique identifier of the trackletdata
      CHARACTER*(*) oid
c INPUT: obscode
      CHARACTER*3 observ
c INPUT: time, UTC MJD
      DOUBLE PRECISION tutc
c INPUT: R for range, V for range-rate
      CHARACTER*1 type
c INPUT: residual in range/range-rate
      DOUBLE PRECISION res
c INPUT: WEIGHT 
      DOUBLE PRECISION weight
c INPUT: chi value from outlier rejection, radar astrometry used in fit (if >0)
      DOUBLE PRECISION chi
      INTEGER selflag
c ==================================================================
c write one record
      WRITE(iunobs,100) oid,tutc,observ,type,res,weight,chi,selflag
 100  FORMAT(A,1x,f17.10,1x,A3,1x,A1,1x,f12.6,1x,f10.6,1x,f9.4,1x,I2)
c later we could implement a smart format, with the appropriate number 
c of digits computed on the basis of the weight.
      RETURN
      END
c ================================================================
c
c RESIDUAL HEADER
c
c ================================================================
      SUBROUTINE radres_header(iunres)
      IMPLICIT NONE
c INPUT: unit (already opened)
      INTEGER iunres
c column headings
      WRITE(iunres,100)
 100  FORMAT('!! oid    time (UTC,MJD)  obscod  typ   Residual  '
     +    ,' Weight   chi   sel ')
      RETURN
      END
c ================================================================
c
c RADAR ASTROMETRY HEADER
c
c ================================================================
      SUBROUTINE radarast_header(iunrada)
      IMPLICIT NONE
c INPUT: unit (already opened)
      INTEGER iunrada
c column  headings
c number  name            date       hour       value         rms unit surf freq transm    receiv
c 100085 1992 UY4         2005-08-04 07:20:00   47331721.00   2.000 us COM  8560 DSS 14    DSS 14 
      WRITE(iunrada,100)
 100  FORMAT('!!number  name            date       hour       value  ',
     + '       rms unit surf freq transm    receiv') 
      RETURN
      END
c ================================================================
c
c READ_DETECTION
c
c ================================================================
      SUBROUTINE read_detection(iunobs,eof,oid,tutc,obs_type,ra,dec,
     +   appmag,filter,observatory,rms_ra,rms_dec,rms_mag,s2n,secnam)
      IMPLICIT NONE
c INPUT: unit number (must be opened already) 
      INTEGER iunobs
c OUTPUT: end-of-file, unique identifier
      LOGICAL eof 
      CHARACTER*(*) oid
c OUTPUT: observation type, filter, obscode
      CHARACTER*1 obs_type
      CHARACTER*1 filter
      CHARACTER*3 observatory
c OUTPUT: time, UTC MJD
      DOUBLE PRECISION tutc
c OUTPUT: R.A., DEC (degrees), apparent magnitude, filter, signal/noise
      DOUBLE PRECISION ra,dec,appmag, s2n
c OUTPUT: RMS of RA*cos(DEC), RMS of DEC, both in deg, RMS of APPMAG
c WARNING: rms_mag must not be used for weighting in a fit to a constant 
c    absolute magnitude, because it only
c    contains the photometric uncertainty, not the lightcurve!
      DOUBLE PRECISION rms_ra, rms_dec, rms_mag
c POSSIBLE OUPUT (NORMALLY IGNORED) object name (known only in simulations)
      CHARACTER*(*) secnam
c END INTERFACE
      CHARACTER*512 record
c ==================================================================
      eof=.false.
c read one record
 1    CONTINUE
      READ(iunobs,'(A)',END=2) record
      IF(record(1:1).eq.'!')THEN
         IF(record(2:2).eq.'!')THEN
!            WRITE(*,*)'read_detection: COLUMN HEADINGS'
c len_trim or rmsp  would be needed...
!            WRITE(*,*) record
         ENDIF
         GOTO 1
      ENDIF
      READ(record,*,ERR=3) oid,tutc,obs_type,ra,dec,appmag,filter,
     +     observatory, rms_ra, rms_dec, rms_mag, s2n, secnam
!      IF(secnam.eq.'NA') secnam=' '
      RETURN
 2    CONTINUE
      eof=.true.
      RETURN
c WARNING: the data are undefined in output if eof
 3    CONTINUE
      WRITE(*,*) ' read-detection: error in reading the record'
      WRITE(*,*) record
      STOP ' ***** error in reading detections****'
      END
c ================================================================
c
c READ_TRACKLET
c
c WARNING: when reading multiple tracklets, use always on the same file
c after eof, a new file can be opened and the routine reused
c ================================================================
      SUBROUTINE read_tracklet(iunobs, nx, eof, oid, nobs,
     +     tutc, obs_type, ra, dec, appmag,
     +     filter, observatory, rms_ra, rms_dec, rms_mag, s2n, secnam)
      IMPLICIT NONE
c INPUT: unit number (must be opened already), max no detections in tracklet 
      INTEGER iunobs,nx
c OUTPUT: end-of-file, unique identifier
      LOGICAL eof 
      CHARACTER*(*) oid
c OUTPUT: no. detections in tracklet
      INTEGER nobs
c OUTPUT: observation type, filter, obscode (these should be the same, but...)
      CHARACTER*1 obs_type(nx)
      CHARACTER*1 filter(nx)
      CHARACTER*3 observatory(nx)
c OUTPUT: time, UTC MJD
      DOUBLE PRECISION tutc(nx)
c OUTPUT: R.A., DEC (degrees), apparent magnitude, filter, signal/noise
      DOUBLE PRECISION ra(nx),dec(nx),appmag(nx), s2n(nx)
c OUTPUT: RMS of RA*cos(DEC), RMS of DEC, both in deg, RMS of APPMAG
c WARNING: rms_mag must not be used for weighting in a fit to a constant 
c    absolute magnitude, because it only
c    contains the photometric uncertainty, not the lightcurve!
      DOUBLE PRECISION rms_ra(nx), rms_dec(nx), rms_mag(nx)
c POSSIBLE OUPUT (NORMALLY IGNORED) object name (known only in simulations)
      CHARACTER*(*) secnam(nx)
c END INTERFACE
      INTEGER j,i
      DATA j/0/
      CHARACTER*9 oid_det
c LOCAL:observation type, filter, obscode (these should be the same, but...)
      CHARACTER*1 obs_t
      CHARACTER*2 filt
      CHARACTER*3 observ
c LOCAL: time, UTC MJD
      DOUBLE PRECISION tut
c LOCAL: R.A., DEC (degrees), apparent magnitude, filter, signal/noise
      DOUBLE PRECISION ras,decs,appmags, s2ns
c LOCAL: RMS of RA*cos(DEC), RMS of DEC, both in deg, RMS of APPMAG
      DOUBLE PRECISION rms_ras, rms_decs, rms_mags
c LOCAL: secret name
      CHARACTER*(9) secna
c save local variables to keep the stored record
      SAVE
c ================================================================== 
      IF(j.eq.1)THEN
         oid=oid_det
c use record from last read
         tutc(j)=tut
         obs_type(j)=obs_t
         filter(j)=filt
         observatory(j)=observ
         ra(j)=ras
         dec(j)=decs
         appmag(j)=appmags
         s2n(j)=s2ns
         rms_ra(j)=rms_ras
         rms_dec(j)=rms_decs
         rms_mag(j)=rms_mags
         secnam(j)=secna
      ENDIF     
 1    CONTINUE
c read one new record
      CALL read_detection(iunobs,eof,oid_det,tut,obs_t,
     +     ras,decs,appmags,filt,observ,    
     +     rms_ras, rms_decs, rms_mags, s2ns, secna)
      IF(eof)THEN
c end of file: terminate current tracklet
         nobs=j
         j=0
         RETURN
      ELSE
c new detection read
         j=j+1
      ENDIF 
      IF(j.eq.1)THEN
c tracklet oid
         oid=oid_det
      ELSEIF(oid.ne.oid_det)THEN
c new tracklet: output the one finished 
         nobs=j-1
c for next time, save the record already read
         j=1
         RETURN
      ENDIF
      IF(j.gt.nx)THEN
        WRITE(*,*)' tracklet too long, max was ', nx, j 
      ELSE
c if same tracklet, use record just read
         tutc(j)=tut
         obs_type(j)=obs_t
         filter(j)=filt
         observatory(j)=observ
         ra(j)=ras
         dec(j)=decs
         appmag(j)=appmags
         s2n(j)=s2ns
         rms_ra(j)=rms_ras
         rms_dec(j)=rms_decs
         rms_mag(j)=rms_mags
         secnam(j)=secna
      ENDIF
C read next, if required
      IF(j.le.nx)THEN  
         GOTO 1
      ELSE
         WRITE(*,*)' read_tracklet: too many detections'
         WRITE(*,*)' read_tracklet: oid=',oid,' nx was ',nx
         STOP '***** error in reading tracklet ****' 
      ENDIF
      END
c ================================================================
c
c WRITE_DETECTION
c
c ================================================================
      SUBROUTINE write_detection(iunobs,oid,tutc,obs_type,ra,dec,
     +    appmag,filter,observatory,rms_ra,rms_dec,rms_mag,s2n,secnam)
      IMPLICIT NONE
c INPUT: unit number (must be opened already) 
      INTEGER iunobs
c INPUT: unique identifier
      CHARACTER*(*) oid
c INPUT: observation type, filter, obscode
      CHARACTER*1 obs_type
      CHARACTER*1 filter
      CHARACTER*3 observatory
c INPUT: time, UTC MJD
      DOUBLE PRECISION tutc
c INPUT: R.A., DEC (degrees), apparent magnitude, filter, signal/noise
      DOUBLE PRECISION ra,dec,appmag, s2n
c INPUT: RMS of RA*cos(DEC), RMS of DEC, both in deg, RMS of APPMAG
      DOUBLE PRECISION rms_ra, rms_dec, rms_mag
c POSSIBLE OUPUT (NORMALLY NOT AVAILABLE) object name (only in simulations)
      CHARACTER*(*) secnam
c ==================================================================
      INTEGER le
      CALL remspace(secnam,le)
      IF(le.eq.0)THEN
         IF(obs_type.EQ.'R'.OR.obs_type.EQ.'V')THEN
            WRITE(iunobs,110) oid,tutc,obs_type,ra,dec,appmag,filter,
     +           observatory,rms_ra,rms_dec,rms_mag,s2n,'NA'
         ELSE
            WRITE(iunobs,100) oid,tutc,obs_type,ra,dec,appmag,filter,
     +           observatory,rms_ra,rms_dec,rms_mag,s2n,'NA'
         ENDIF
      ELSE
c write one record
         IF(obs_type.EQ.'R'.OR.obs_type.EQ.'V')THEN
            WRITE(iunobs,110) oid,tutc,obs_type,ra,dec,appmag,filter,
     +           observatory,rms_ra,rms_dec,rms_mag,s2n,secnam(1:le)
         ELSE
            WRITE(iunobs,100) oid,tutc,obs_type,ra,dec,appmag,filter,
     +           observatory,rms_ra,rms_dec,rms_mag,s2n,secnam(1:le)
         ENDIF
      ENDIF
 100  FORMAT(A,1x,f17.10,1x,A1,1x,f16.12,1x,f16.12,1x,f13.10,1x,A1,1x
     +       ,A3,1x,f9.4,1x,f9.4,1x,f6.3,1x,f11.3,1X,A)
 110  FORMAT(A,1x,f17.10,1x,A1,1x,f16.4,1x,f16.4,1x,f13.10,1x,A1,1x
     +       ,A3,1x,f9.4,1x,f9.4,1x,f6.3,1x,f11.3,1X,A)
c later we could implement a smart format, with the appropriate number 
c of digits computed on the basis of the rms.
      RETURN
      END
c ================================================================
c
c WRITE_OBSPOS
c
c ================================================================
      SUBROUTINE write_obspos(iunobs,obspos)
      IMPLICIT NONE
      INTEGER iunobs
      DOUBLE PRECISION obspos(3)
      WRITE(iunobs,100) obspos
 100  FORMAT(3(2X,F12.5)) 
      END

c ================================================================
c
c READ_OBSPOS
c
c ================================================================
      SUBROUTINE read_obspos(iunobs,obspos,eof)
      IMPLICIT NONE
      INTEGER iunobs
      DOUBLE PRECISION obspos(3)
      LOGICAL eof
      eof=.false.
      READ(iunobs,100,END=3) obspos
 100  FORMAT(3(2X,F12.5)) 
      RETURN
 3    CONTINUE
      eof=.true.
      END

c ================================================================
c
c WRITE_TRACKLET
c
c ================================================================
      SUBROUTINE write_tracklet(iunobs, oid, nx,
     +     tutc, obs_type, ra, dec, appmag,
     +     filter, observatory, rms_ra, rms_dec, rms_mag, s2n,secnam)
      IMPLICIT NONE
c INPUT: unit number (must be opened already)
      INTEGER iunobs
c INPUT: unique identifier
      CHARACTER*(*) oid
c INPUT: no. detections in tracklet
      INTEGER nx
c INPUT: observation type, filter, obscode (these should be the same, but...)
      CHARACTER*1 obs_type(nx)
      CHARACTER*1 filter(nx)
      CHARACTER*3 observatory(nx)
c INPUT: time, UTC MJD
      DOUBLE PRECISION tutc(nx)
c INPUT: R.A., DEC (degrees), apparent magnitude, filter, signal/noise
      DOUBLE PRECISION ra(nx),dec(nx),appmag(nx), s2n(nx)
c INPUT: RMS of RA*cos(DEC), RMS of DEC, both in deg, RMS of APPMAG
      DOUBLE PRECISION rms_ra(nx), rms_dec(nx), rms_mag(nx)
c POSSIBLE INPUT (NORMALLY NOT AVAILABLE) object name (known only in simulations)
      CHARACTER*(*) secnam(nx)
c END INTERFACE
      INTEGER j
c ==================================================================
      DO j=1,nx
         CALL write_detection(iunobs,oid,tutc(j),obs_type(j),
     +    ra(j),dec(j),appmag(j),filter(j),observatory(j), 
     +    rms_ra(j), rms_dec(j), rms_mag(j), s2n(j), secnam(j))
      ENDDO
      RETURN
      END
c ================================================================
c
c DETECTION HEADER
c
c ================================================================
      SUBROUTINE detection_header(iunobs)
      IMPLICIT NONE
c INPUT: unit (already opened)
      INTEGER iunobs
c column headings
      WRITE(iunobs,100)
 100  FORMAT('!! oid  time(UTC,MJD)  obs_type      ra             dec  '
     +    ,'        app_mag  filt obscod  rms_ra   rms_dec '
     +    ,' rms_mag     s2n  ')
      RETURN
      END
c ---------------------------------------------------------------------
c
c  *****************************************************************
c  *                                                               *
c  *                      REMSPACE                                 *
c  *                                                               *
c  *            Remove spaces from a character string              *
c  *                                                               *
c  *****************************************************************
c
c INPUT:    C         -  Character string
c
c OUTPUT:   C         -  Character string without spaces
c           L         -  Length of the output character string
c
      SUBROUTINE remspace(c,l)
      IMPLICIT NONE
      INTEGER l,i,l1
      CHARACTER c*(*)
      CHARACTER*(512) c1
      l=LEN(c)
      IF(l.LE.0.or.l.gt.512) STOP '**** remspace: problems ****'
      c1=' '
      l1=0
      DO i=1,l
         IF(c(i:i).ne.' ')THEN
            l1=l1+1
            c1(l1:l1)=c(i:i)
         ENDIF
      ENDDO
      l=l1
      IF(l.LE.0) THEN
         c=' '
      ELSE
         c=c1(1:l1)
      ENDIF
      END

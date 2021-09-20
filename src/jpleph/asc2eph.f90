      PROGRAM ASC2EPH 
!                                                                       
!      ASC2EPH creates a binary format JPL Planetary Ephemeris file from
!      one or more ascii text files.                                    
!                                                                       
!$ Disclaimer                                                           
!                                                                       
!     THIS SOFTWARE AND ANY RELATED MATERIALS WERE CREATED BY THE       
!     CALIFORNIA INSTITUTE OF TECHNOLOGY (CALTECH) UNDER A U.S.         
!     GOVERNMENT CONTRACT WITH THE NATIONAL AERONAUTICS AND SPACE       
!     ADMINISTRATION (NASA). THE SOFTWARE IS TECHNOLOGY AND SOFTWARE    
!     PUBLICLY AVAILABLE UNDER U.S. EXPORT LAWS AND IS PROVIDED "AS-IS" 
!     TO THE RECIPIENT WITHOUT WARRANTY OF ANY KIND, INCLUDING ANY      
!     WARRANTIES OF PERFORMANCE OR MERCHANTABILITY OR FITNESS FOR A     
!     PARTICULAR USE OR PURPOSE (AS SET FORTH IN UNITED STATES UCC      
!     SECTIONS 2312-2313) OR FOR ANY PURPOSE WHATSOEVER, FOR THE        
!     SOFTWARE AND RELATED MATERIALS, HOWEVER USED.                     
!                                                                       
!     IN NO EVENT SHALL CALTECH, ITS JET PROPULSION LABORATORY, OR NASA 
!     BE LIABLE FOR ANY DAMAGES AND/OR COSTS, INCLUDING, BUT NOT        
!     LIMITED TO, INCIDENTAL OR CONSEQUENTIAL DAMAGES OF ANY KIND,      
!     INCLUDING ECONOMIC DAMAGE OR INJURY TO PROPERTY AND LOST PROFITS, 
!     REGARDLESS OF WHETHER CALTECH, JPL, OR NASA BE ADVISED, HAVE      
!     REASON TO KNOW, OR, IN FACT, SHALL KNOW OF THE POSSIBILITY.       
!                                                                       
!     RECIPIENT BEARS ALL RISK RELATING TO QUALITY AND PERFORMANCE OF   
!     THE SOFTWARE AND ANY RELATED MATERIALS, AND AGREES TO INDEMNIFY   
!     CALTECH AND NASA FOR ALL THIRD-PARTY CLAIMS RESULTING FROM THE    
!     ACTIONS OF RECIPIENT IN THE USE OF THE SOFTWARE.                  
!                                                                       
!                                                                       
!      This program, 'asc2eph', requires (via standard input) an ascii  
!      header file ('header.XXX'), followed by one or more ascii ephemer
!      data files ('ascSYYYY.XXX').  All files must have the same epheme
!      number, XXX.  Further, the data files must be consecutive in time
!      with no gaps between them.                                       
!                                                                       
!      By default, the output ephemeris will span the same interval as t
!      text file(s).  If you are interested in only a portion of data, s
!      below T1 and T2 to the begin and end times of the span you desire
!      and T2 must be specified in  Julian Ephemeris Days (ET).         
!                                                                       
!      A sample sequence of files might be:                             
!                                                                       
!        header.405  asc+1920.405 asc+1940.405 asc+1960.405 asc+1980.405
!                                                                       
!      This program is written in standard Fortran-77.                  
!                                                                       
! **********************************************************************
!                                                                       
!                                    *** NOTE ***                       
!                                                                       
!      However, the units in which the length of a direct access record 
!      are PROCESSOR DEPENDENT.  The parameter NRECL, the number of unit
!      controls the length of a record in the direct access ephemeris.  
!      The user MUST select the correct value of NRECL by editing one of
!      comemnted lines defining PARAMETER (NRECL) below.                
!                                                                       
! **********************************************************************
!                                                                       
!     Updated 02 March  2013 to accommodate more than 400 dynamical para
!     Updated 02 August 2013 to accommodate reading of TT-TDB           
!     Updated 15 August 2013 to accommodate negative Julian dates       
!                                                                       
! **********************************************************************
                                                                        
      IMPLICIT NONE 
                                                                        
      INTEGER NRECL 
                                                                        
! *****  The user must choose one of the following statements  *****    
!            ( Usually NRECL = 4 is used on Unix platforms)             
                                                                        
      PARAMETER ( NRECL = 4 ) 
!      PARAMETER ( NRECL = 1 )                                          
                                                                        
      INTEGER OLDMAX 
      PARAMETER ( OLDMAX = 400) 
      INTEGER NMAX 
      PARAMETER ( NMAX = 1000) 
                                                                        
! **********************************************************************
                                                                        
       CHARACTER*6  CNAM (NMAX) 
       CHARACTER*6  TTL  (14,3) 
       CHARACTER*12 HEADER 
                                                                        
       DOUBLE PRECISION  AU 
       DOUBLE PRECISION  CVAL (NMAX) 
       DOUBLE PRECISION  DB(3000) 
       DOUBLE PRECISION  DB2Z 
       DOUBLE PRECISION  EMRAT 
       DOUBLE PRECISION  SS   (3) 
       DOUBLE PRECISION  T1, T2 
                                                                        
       INTEGER  I, IRECSZ, IN 
       INTEGER  J 
       INTEGER  K,KK,KP2,KSIZE 
                                                   ! Pointers to number 
       INTEGER  IPT (3,12) 
                                                   ! Pointer to number o
       INTEGER  LPT(3) 
                                                   ! Pointer to number o
       INTEGER  RPT(3) 
                                                   ! Pointer to number o
       INTEGER  TPT(3) 
       INTEGER  N, NCON, NROUT, NCOEFF, NRW, NUMDE 
       INTEGER  OUT 
                                                                        
       LOGICAL  FIRST 
       DATA     FIRST / .TRUE. / 
                                                                        
! **********************************************************************
!                                                                       
!     By default, the output ephemeris will span the same interval as th
!     input ascii data file(s).  The user may reset these to other JED's
!                                                                       
       DB2Z = -99999999.d0 
       T1   = -99999999.d0 
       T2   =  99999999.d0 
                                                                        
       IF(NRECL.NE.1.AND.NRECL.NE.4) THEN 
         WRITE(*,*)'*** ERROR: User did not set NRECL ***' 
         STOP 
       ENDIF 
                                                                        
!      Write a fingerprint to the screen.                               
                                                                        
       WRITE(*,*) ' JPL ASCII-TO-DIRECT-I/O program. ' //               &
     &            ' Last modified 15-Aug-2013.'                         
                                                                        
!      Read the size and number of main ephemeris records.              
!        (from header.xxx)                                              
                                                                        
       READ  (*,'(6X,I6)')  KSIZE 
       WRITE (*,100)              KSIZE 
  100  FORMAT(/'KSIZE =',I6) 
                                                                        
       IRECSZ = NRECL * KSIZE 
                                                                        
!      Now for the alphameric heading records (GROUP 1010)              
                                                                        
       CALL  NXTGRP ( HEADER ) 
                                                                        
       IF (HEADER.NE.'GROUP   1010')  THEN 
          CALL  ERRPRT ( 1010, 'NOT HEADER' ) 
       ENDIF 
                                                                        
       READ (*,'(14A6)')    TTL 
       WRITE(*,'(/(14A6))') TTL 
                                                                        
!      Read start, end and record span  (GROUP 1030)                    
                                                                        
       CALL  NXTGRP ( HEADER ) 
                                                                        
       IF ( HEADER.NE.'GROUP   1030' ) THEN 
          CALL ERRPRT ( 1030, 'NOT HEADER' ) 
       ENDIF 
                                                                        
       READ (*,'(3D12.0)')  SS 
                                                                        
!      Read number of constants and names of constants (GROUP 1040/4).  
                                                                        
       CALL  NXTGRP ( HEADER ) 
                                                                        
       IF ( HEADER.NE.'GROUP   1040 ') THEN 
          CALL ERRPRT ( 1040, 'NOT HEADER' ) 
       ENDIF 
                                                                        
       READ (*,'(I6)')    N 
       READ (*,'(10A8)') (CNAM(I),I=1,N) 
                                                                        
       NCON = N 
                                                                        
!      Read number of values and values (GROUP 1041/4)                  
                                                                        
       CALL  NXTGRP ( HEADER ) 
                                                                        
       IF ( HEADER.NE.'GROUP   1041' ) THEN 
          CALL ERRPRT ( 1041, 'NOT HEADER' ) 
       ENDIF 
                                                                        
       READ (*,'(I6)')       N 
       READ (*,*)  (CVAL(I),I=1,N) 
                                                                        
       DO  I = 1, N 
           IF ( CNAM(I).EQ.'AU    ' )  AU    = CVAL(I) 
           IF ( CNAM(I).EQ.'EMRAT ' )  EMRAT = CVAL(I) 
           IF ( CNAM(I).EQ.'DENUM ' )  NUMDE = CVAL(I) 
       END DO 
                                                                        
       WRITE (*,'(500(/2(A8,D24.16)))') (CNAM(I),CVAL(I),I=1,N) 
                                                                        
!      Zero out pointer arrays                                          
                                                                        
       DO I = 1,3 
         DO J = 1,12 
           IPT(I,J) = 0 
         ENDDO 
         LPT(I) = 0 
         RPT(I) = 0 
         TPT(I) = 0 
       ENDDO 
                                                                        
!      Read pointers needed by INTERP (GROUP 1050)                      
                                                                        
       CALL  NXTGRP ( HEADER ) 
                                                                        
       IF ( HEADER.NE.'GROUP   1050' ) THEN 
          CALL ERRPRT ( 1050, 'NOT HEADER' ) 
       ENDIF 
                                                                        
       WRITE(*,'(/)') 
                                                                        
       DO I=1,3 
         READ (*,'(15I6)') (IPT(I,J),J=1,12),LPT(I),RPT(I),TPT(I) 
         WRITE(*,'(15I5)')(IPT(I,J),J=1,12),LPT(I),RPT(I),TPT(I) 
       ENDDO 
                                                                        
!      Open direct-access output file ('JPLEPH')                        
                                                                        
       OPEN ( UNIT   = 12,                                              &
     &        FILE   = 'JPLEPH',                                        &
     &        ACCESS = 'DIRECT',                                        &
     &        FORM   = 'UNFORMATTED',                                   &
     &        RECL   = IRECSZ,                                          &
     &        STATUS = 'NEW' )                                          
                                                                        
!     Read and write the ephemeris data records (GROUP 1070).           
                                                                        
      CALL  NXTGRP ( HEADER ) 
                                                                        
      IF ( HEADER.NE.'GROUP   1070' ) CALL ERRPRT(1070,'NOT HEADER') 
                                                                        
      NROUT  = 0 
      IN     = 0 
      OUT    = 0 
                                                                        
    1 continue 
                                                                        
      READ(*,'(2I6)')NRW,NCOEFF 
      IF(NRW.EQ.0) GO TO 1 
                                                                        
      DO K=1,NCOEFF,3 
        KP2 = MIN(K+2,NCOEFF) 
        READ(*,*,IOSTAT = IN)(DB(KK),KK=K,KP2) 
        IF(IN.NE.0)STOP ' Error reading 1st set of coeffs' 
      ENDDO 
                                                                        
      DO WHILE (       ( IN.EQ.0 ).AND.( DB(2).LT.T2) )                               
                                                                        
          IF ( 2*NCOEFF.NE.KSIZE ) THEN 
             CALL ERRPRT(NCOEFF,' 2*NCOEFF not equal to KSIZE') 
          ENDIF 
                                                                        
!         Skip this data block if the end of the interval is less       
!         than the specified start time or if the it does not begin     
!         where the previous block ended.                               
                                                                        
          IF  ( (DB(2).GE.T1) .AND. (DB(1).GE.DB2Z) ) THEN 
                                                                        
             IF ( FIRST ) THEN 
                                                                        
!               Don't worry about the intervals overlapping             
!               or abutting if this is the first applicable             
!               interval.                                               
                                                                        
                DB2Z  = DB(1) 
                FIRST = .FALSE. 
             ENDIF 
                                                                        
             IF (DB(1).NE.DB2Z ) THEN 
                                                                        
!               Beginning of current interval is past the end           
!               of the previous one.                                    
                                                                        
                CALL ERRPRT (NRW, 'Records do not overlap or abut') 
             ENDIF 
                                                                        
             DB2Z  = DB(2) 
             NROUT = NROUT + 1 
                                                                        
             WRITE (12,REC=NROUT+2,IOSTAT=OUT) (DB(K),K=1,NCOEFF) 
                                                                        
             IF ( OUT.NE.0 ) THEN 
                CALL ERRPRT (NROUT,                                     &
     &                     'th record not written because of error')    
             ENDIF 
                                                                        
!            Save this block's starting date, its interval span, and its
!            date.                                                      
                                                                        
             IF (NROUT.EQ.1) THEN 
                SS(1) = DB(1) 
                SS(3) = DB(2) - DB(1) 
             ENDIF 
             SS(2) = DB(2) 
                                                                        
!            Update the user as to our progress every 100th block.      
                                                                        
             IF ( MOD(NROUT,100).EQ.1 ) THEN 
                IF ( DB(1).GE.T1 ) THEN 
                   WRITE (*,271) NROUT, DB(2) 
                ELSE 
                   WRITE (*,*)                                          &
     &           ' Searching for first requested record...'             
                ENDIF 
             ENDIF 
  271        FORMAT (I6,                                                &
     &              ' EPHEMERIS RECORDS WRITTEN.  LAST JED = ',         &
     &              F12.2)                                              
                                                                        
          ENDIF 
                                                                        
          READ (*,'(2I6)',IOSTAT =IN) NRW, NCOEFF 
                                                                        
          IF(IN.EQ.0)then 
            DO K=1,NCOEFF,3 
              KP2 = MIN(K+2,NCOEFF) 
              READ(*,*,IOSTAT = IN)(DB(KK),KK=K,KP2) 
              IF(IN.NE.0)STOP ' Error reading nth set of coeffs' 
            ENDDO 
          ENDIF 
                                                                        
      END DO 
                                                                        
      WRITE (*,275) NROUT, DB(2) 
  275 FORMAT(I6, ' EPHEMERIS RECORDS WRITTEN.  LAST JED = ',F12.2) 
                                                                        
!     Write header records onto output file.                            
                                                                        
      NROUT = 1 
                                                                        
      IF(NCON .LE. OLDMAX)THEN 
         WRITE(12,REC=1,IOSTAT=OUT) TTL,(CNAM(I),I=1,OLDMAX),SS,NCON,AU,&
     &              EMRAT,IPT,NUMDE,LPT,RPT,TPT                         
      ELSE 
         K = OLDMAX+1 
         WRITE(12,REC=1,IOSTAT=OUT) TTL,(CNAM(I),I=1,OLDMAX),SS,NCON,AU,&
     &              EMRAT,IPT,NUMDE,LPT,(CNAM(J),J=K,NCON),RPT,TPT      
      ENDIF 
                                                                        
      IF ( OUT .NE. 0 ) THEN 
          CALL ERRPRT ( NROUT, 'st record not written because of error') 
      ENDIF 
                                                                        
      NROUT = 2 
                                                                        
      IF(NCON .LE. OLDMAX) THEN 
         WRITE(12,REC=2,IOSTAT=OUT)(CVAL(I),I=1,OLDMAX) 
      ELSE 
         WRITE(12,REC=2,IOSTAT=OUT)(CVAL(I),I=1,NCON) 
      ENDIF 
                                                                        
      IF ( OUT .NE. 0 ) THEN 
          CALL ERRPRT ( NROUT, 'nd record not written because of error') 
      ENDIF 
                                                                        
!     We're through.  Wrap it up.                                       
                                                                        
      CLOSE (12) 
      STOP ' OK' 
                                                                        
      END                                           
                                                                        
                                                                        
      SUBROUTINE  ERRPRT (I, MSG) 
                                                                        
      CHARACTER*(*)  MSG 
      INTEGER        I 
                                                                        
      WRITE (*,200)  I, MSG 
  200 FORMAT('ERROR #',I8,2X,A50) 
                                                                        
      STOP ' ERROR ' 
      END                                           
                                                                        
                                                                        
      SUBROUTINE  NXTGRP ( HEADER ) 
                                                                        
      CHARACTER*(*)  HEADER 
      CHARACTER*12   BLANK 
                                                                        
!     Start with nothing.                                               
                                                                        
      HEADER = ' ' 
                                                                        
!     The next non-blank line we encounter is a header record.          
!     The group header and data are seperated by a blank line.          
                                                                        
      DO WHILE ( HEADER .EQ. ' ' ) 
          READ (*,'(A)') HEADER 
      ENDDO 
                                                                        
!     Found the header.  Read the blank line so we can get at the data. 
                                                                        
      IF ( HEADER .NE. 'GROUP   1070' )READ (*, '(A)') BLANK 
                                                                        
      RETURN 
      END                                           

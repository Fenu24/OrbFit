! MODULE io_elems

MODULE io_elems

! ---------------------------------------------------------------------
! former COMORB.H
! Information on opened orbital element file
!
! orbunt      -  FORTRAN unit
! orbnr       -  Number of records read so far
! orbfn       -  File name
! dfrsty      -  default reference system type
! dfrsep      -  default reference system epoch
! rectyp      -  record type
! deltyp      -  default orbital element type
! depstr      -  default epoch (character string)
! dept0       -  default epoch (MJD, TDT)
! deft0       -  is dept0 defined?
! nxtend      -  end encountered at previous call
! iicorb      -  initialization check
!
  INTEGER orbunt,orbnr,iicorb
  CHARACTER*100 orbfn,depstr
  DOUBLE PRECISION dept0
  CHARACTER*4 dfrsty,rectyp,deltyp
  CHARACTER*10 dfrsep
  LOGICAL deft0,nxtend
  CHARACTER*20 oef_vers
  PUBLIC orbunt,orbnr,iicorb,orbfn,depstr,dept0,dfrsty,rectyp,deltyp, &
  &      dfrsep,deft0,nxtend,oef_vers
! dirty fix for passing the obscode (for ATT elelemnts) from rdorb without changing interface 
  CHARACTER*3, PUBLIC ::  obscod

  !PUBLIC :: oefdet, oporbf, clorbf, fixcnm

END MODULE io_elems

!  LIBRARY io_elems ELEMENTS  INPUT
!       TOP LEVEL
! rdelem	read orbital els for a list of objects from a list of files    
! oefdet	auto-detects format of orbital element files
!        READING DIFFERENT FORMATS  
! rdast2	read orbital elements for list of obj (Bowell's astorb.dat format) 
! rdmpca	read orbital elements for list of obj (MPC format for asteroids)  
! rdmpca3       read all orbital elements (MPC format for asteroids) 
! rdoef		read orbital elements (OEF format)   
! oporbf	open an orbital element file (OEF format)                     
! rdorb	        read orbital element record (OEF format)                 
! clorbf	close an orbital element file (OEF format)                    
!  ELEMENTS OUTPUT
! outele	verbose output of a set orbital elements to a report file      
! wro1lh	writes the header of an orbital element file (1L format)       
! wro1lr	writes an orbital element record in an orbital element file (1L
! wro1lr_ngr    writes a record of dyn parameters in 1L catalog
! wromlh	writes the header of an orbital element file (ML format)       
!
!  CONVERSIONS                                                                     
! mpcdat	computes MJD from MPC-style packed dates   
! fixcnm	generates normal matrix from covariance matrix or viceversa   


! Copyright (C) 1997-1999 by Mario Carpino (carpino@brera.mi.astro.it)  
!                            Andrea Milani (milani@dm.unipi.it)         
!                                                                       
! Version: February 12, 1999                                            
! Revised by Genny on 5 November 2005 to deal with the new format of 
! astorb.dat to include numbered objects > 99999
! --------------------------------------------------------------------- 
!                                                                       
!  *****************************************************************    
!  *                                                               *    
!  *                         O E F D E T                           *    
!  *                                                               *    
!  *         Auto-detects format of orbital element files          *    
!  *                                                               *    
!  *****************************************************************    
!                                                                       
! INPUT:    UNIT      -  Input FORTRAN unit (must be already opened)    
!           FILNAM    -  Input file name (for error messages)           
!                                                                       
! OUTPUT:   FORM      -  Format:                                        
!                          1) OEF:   ORBFIT orbital element file        
!                          2) BA1:   Bowell's astorb.dat (pre-1999)     
!                          3) BA2:   Bowell's astorb.dat (post-1999)    
!                          4) MPC-A: MPC (asteroids)                    
!                          5) BAC:   Bowell private format with C       
! in the future:                                                        
!                          6) MPC-C: MPC (comets)                       
!                                                                       
SUBROUTINE oefdet(unit,filnam,form) 
  IMPLICIT NONE                                                                       
  INTEGER unit 
  CHARACTER*(*) filnam,form                                                                 
  INCLUDE 'parcmc.h90'                                                  
! Number of supported formats                                           
  INTEGER, PARAMETER :: nfx=5                                            
  INTEGER lf,lr,i,k,npf,ipf,lrwb 
  LOGICAL poss(nfx),neoh,error 
  CHARACTER rec*300,recwb*200,b4*4,p4*4 
  CHARACTER*50 tmp,tmp1,tmp2                                                           
  INTEGER, EXTERNAL :: lench 
! Names of supported formats                                            
  CHARACTER*10 fname(nfx) 
  DATA fname/'OEF','BA1','BA2','MPC-A','BAC'/ 
  form=' '                                                               
! METHOD: at the beginning, we flag as possible all the supported       
! formats. Then we scan the first records of the file and discard       
! formats as we find lines which are not compatible with them.          
! At the end, detection is successful if one and only one format        
! remains.                                                              
  poss(1:nfx)=.true.                                                                       
  neoh=.true. 
! Scan only first 100 records of the file                               
  DO 2 k=1,100 
      READ(unit,100,END=3) rec 
100   FORMAT(A) 
      lr=lench(rec)                                                              
! Format 1 (OEF)                                                        
      IF(poss(1).AND.neoh) THEN 
         recwb=rec 
         CALL rmsp(recwb,lrwb) 
         IF(lrwb.GT.0) THEN 
            IF(recwb.EQ.'END_OF_HEADER') THEN 
               neoh=.false. 
            ELSEIF(recwb(1:1).NE.comcha) THEN 
               IF(index(recwb,'=').EQ.0) THEN 
                  poss(1)=.false. 
               ELSEIF(recwb(1:7).EQ.'format=') THEN 
                  tmp=recwb(8:) 
                  i=index(tmp,comcha) 
                  IF(i.GT.0) THEN 
                     tmp1=tmp(1:i-1) 
                     tmp=tmp1 
                  END IF
                  CALL strcnt(tmp,tmp1,tmp2,error) 
                  IF(error) THEN 
                     poss(1)=.false. 
                  ELSEIF(tmp1.NE.'OEF1.1'.AND.tmp1.NE.'OEF2.0') THEN 
                     poss(1)=.false. 
                  END IF
               END IF
            END IF
         END IF
      END IF                                                            
! Format 2 (BA1)                                                        
      IF(poss(2)) THEN 
         IF(lr.NE.265) poss(2)=.false. 
         b4(1:1)=rec(5:5) 
         b4(2:2)=rec(24:24) 
         b4(3:3)=rec(40:40) 
         b4(4:4)=rec(46:46) 
         IF(b4.NE.'    ') poss(2)=.false. 
         p4(1:1)=rec(117:117) 
         p4(2:2)=rec(128:128) 
         p4(3:3)=rec(139:139) 
         p4(4:4)=rec(149:149) 
         IF(p4.NE.'....') poss(2)=.false. 
      END IF                                                               
! Format 3 (BA2)                                                        
      IF(poss(3)) THEN 
         IF(lr.NE.267) poss(3)=.false. 
         b4(1:1)=rec(7:7) 
         b4(2:2)=rec(26:26) 
         b4(3:3)=rec(42:42) 
         b4(4:4)=rec(48:48) 
         IF(b4.NE.'    ') poss(3)=.false. 
         p4(1:1)=rec(119:119) 
         p4(2:2)=rec(130:130) 
         p4(3:3)=rec(141:141) 
         p4(4:4)=rec(151:151) 
         IF(p4.NE.'....') poss(3)=.false. 
      END IF                                                               
! Format 4 (MPC-A)                                                      
      IF(poss(4)) THEN 
         b4(1:1)=rec(8:8) 
         b4(2:2)=rec(26:26) 
         b4(3:3)=rec(58:58) 
         b4(4:4)=rec(80:80) 
         IF(b4.NE.'    ') poss(4)=.false. 
         p4(1:1)=rec(30:30) 
         p4(2:2)=rec(41:41) 
         p4(3:3)=rec(52:52) 
         p4(4:4)=rec(63:63) 
         IF(p4.NE.'....') poss(4)=.false. 
      END IF                                                               
! Format 5 (BAC)                                                        
      IF(poss(5)) THEN 
          IF(mod(k,9).eq.2) THEN 
              IF(rec(1:1).eq.'B')THEN 
                  poss(5)=.true. 
              ELSEIF(rec(1:1).eq.'M')THEN 
                  poss(5)=.true. 
              ELSE 
                  poss(5)=.false. 
              END IF 
          ELSE 
              poss(5)=.false. 
          END IF 
      END IF 
                                                                        
2  END DO
3  CONTINUE 
   IF(neoh) poss(1)=.false. 
                                                                        
! Final check                                                           
   ipf=0 
   npf=0 
   DO 10 i=1,nfx 
      IF(poss(i)) THEN 
         npf=npf+1 
         ipf=i 
      END IF
10 END DO
   IF(npf.EQ.1) form=fname(ipf) 
   IF(form.EQ.' ') THEN 
      lf=lench(filnam) 
      WRITE(*,200) filnam(1:lf) 
   END IF
200 FORMAT('ERROR: format auto-detection failed for file "',A,'":'/   &
     &       '       please specify format explicitly')  
 END SUBROUTINE oefdet

! Copyright (C) 1997-1999 by Mario Carpino (carpino@brera.mi.astro.it)  
! Version: February 12, 1999                                            
! Revised by Genny on 5 November 2005 to deal with the new format of 
! astorb.dat to include numbered objects > 99999
! --------------------------------------------------------------------- 
!                                                                       
!  *****************************************************************    
!  *                                                               *    
!  *                         R D A S T 2                           *    
!  *                                                               *    
!  *          Read orbital elements for a list of objects          *    
!  *       from a file written in Bowell's astorb.dat format       *    
!  *          (version reading NEW, post-1999 format)              *    
!  *                                                               *    
!  *****************************************************************    
!                                                                       
! INPUT:    UNIT      -  Input FORTRAN unit (must be already opened)    
!           FILNAM    -  Input file name (for error messages)           
!           OBJNAM    -  Object names                                   
!           NOBJ      -  Number of objects                              
!                                                                       
! OUTPUT:   DEFORB    -  Tells whether orbital elements are defined     
!           DEFCN     -  Tells whether covariance/normal matrices       
!                            are defined                                
!           ELTYPE    -  Type of orbital elements (EQU/KEP/CAR)         
!           TELEM     -  Epoch of orbital elements (MJD, TDT)           
!           ELEM      -  Orbital elements (ECLM J2000)                  
!           COVE      -  Covariance matrix of orbital elements          
!           NORE      -  Normal matrix of orbital elements              
!           MASS      -  Mass (solar masses)                            
!           H         -  H absolute magnitude (if <-100, missing)       
!           G         -  G slope parameter                              
!           COMELE    -  Comment on orbital elements                    
!                                                                       
! WARNING: the routine assumes that objects having DEFORB=.true.        
!          have already orbital elements defined (possibly from another 
!          input file) and does not overwrite them                      
!                                                                       
! OBJECT NAME TRANSLATION: all names of objects appearing in the        
! file are modified (removing all blanks) before comparison with        
! the name requested by the calling module                              
!                                                                       
SUBROUTINE rdast2(unit,filnam,objnam,nobj,deforb,defcn,           &
     &                  eltype,telem,elem,cove,nore,mass,h,g,comele) 
  USE fund_const   
  IMPLICIT NONE                                                                       
  INTEGER unit,nobj 
  DOUBLE PRECISION telem(nobj),elem(6,nobj),cove(6,6,nobj) 
  DOUBLE PRECISION nore(6,6,nobj),mass(nobj),h(nobj),g(nobj) 
  CHARACTER*(*) filnam,objnam(nobj),eltype(nobj),comele(nobj) 
  LOGICAL deforb(nobj),defcn(nobj)                                                                       
  INTEGER ln,nr,lf,nrem,k,flags(6),year,month,day,lc 
  DOUBLE PRECISION el1(6) 
  CHARACTER n1*6,n2*18,name*18,hc*5,gc*5,krc*10                                                   
  INTEGER lench 
  DOUBLE PRECISION tjm1 
  EXTERNAL lench,tjm1                                                                 
! Number of remaining object (orbit not yet found)                      
  nrem=0 
  DO 10 k=1,nobj 
     IF(deforb(k)) GOTO 10 
     nrem=nrem+1 
10 END DO
  IF(nrem.LE.0) RETURN 
  lf=lench(filnam)                                                                       
  nr=0 
1 CONTINUE 
  READ(unit,100,END=2) n1,n2 
  nr=nr+1 
  IF(n1.EQ.'     ') THEN 
     name=n2 
  ELSE 
     name=n1 
  END IF
  CALL rmsp(name,ln) 
  IF(ln.LE.0) THEN 
     WRITE(*,200) filnam(1:lf),nr 
  200 FORMAT('ERROR in reading file "',A,'": no object name at record',I6)
     GOTO 1 
  END IF
  DO 3 k=1,nobj 
     IF(deforb(k)) GOTO 3 
     IF(name.EQ.objnam(k)) THEN 
        BACKSPACE(unit) 
        READ(unit,100) n1,n2,hc,gc,flags,year,month,day,el1 
100     FORMAT(A6,1X,A18,17X,A5,1X,A5,17X,6I4,12X,I4,2I2,1X,              &
     &       3(F10.6,1X),F9.6,1X,F10.8,1X,F12.8)    
        deforb(k)=.true. 
        defcn(k)=.false. 
        eltype(k)='KEP' 
        telem(k)=tjm1(day,month,year,0.D0) 
        elem(1,k)=el1(6) 
        elem(2,k)=el1(5) 
        elem(3,k)=el1(4)*radeg 
        elem(4,k)=el1(3)*radeg 
        elem(5,k)=el1(2)*radeg 
        elem(6,k)=el1(1)*radeg 
        mass(k)=0.d0 
        IF(hc.EQ.'     ') THEN 
           h(k)=-1.D9 
        ELSE 
           READ(hc,101) h(k) 
        END IF
        IF(gc.EQ.'     ') THEN 
           g(k)=0.15D0 
        ELSE 
           READ(gc,101) g(k) 
101        FORMAT(F5.2) 
        END IF
        WRITE(krc,107) nr 
107     FORMAT(I6) 
        CALL rmsp(krc,lc) 
        comele(k)='read from file "'//filnam(1:lf)//                  &
     &              '" at record '//krc(1:lc)  
        nrem=nrem-1 
        IF(nrem.LE.0) RETURN 
     END IF
3 END DO                                                                      
  GOTO 1 
2 CONTINUE 
END SUBROUTINE rdast2                                       
! Copyright (C) 1997-1998 by Mario Carpino (carpino@brera.mi.astro.it)  
! Version: June 16, 1998                                                
! --------------------------------------------------------------------- 
!                                                                       
!  *****************************************************************    
!  *                                                               *    
!  *                         R D M P C A                           *    
!  *                                                               *    
!  *          Read orbital elements for a list of objects          *    
!  *        from a file written in MPC format for asteroids        *    
!  *                                                               *    
!  *****************************************************************    
!                                                                       
! INPUT:    UNIT      -  Input FORTRAN unit (must be already opened)    
!           FILNAM    -  Input file name (for error messages)           
!           OBJNAM    -  Object names                                   
!           NOBJ      -  Number of objects                              
!                                                                       
! OUTPUT:   DEFORB    -  Tells whether orbital elements are defined     
!           DEFCN     -  Tells whether covariance/normal matrices       
!                            are defined                                
!           ELTYPE    -  Type of orbital elements (EQU/KEP/CAR)         
!           TELEM     -  Epoch of orbital elements (MJD, TDT)           
!           ELEM      -  Orbital elements (ECLM J2000)                  
!           COVE      -  Covariance matrix of orbital elements          
!           NORE      -  Normal matrix of orbital elements              
!           MASS      -  Mass (solar masses)                            
!           H         -  H absolute magnitude (if <-100, missing)       
!           G         -  G slope parameter                              
!           COMELE    -  Comment on orbital elements                    
!                                                                       
! WARNING: the routine assumes that objects having DEFORB=.true.        
!          have already orbital elements defined (possibly from another 
!          input file) and does not overwrite them                      
!                                                                       
SUBROUTINE rdmpca(unit,filnam,objnam,nobj,deforb,defcn,           &
     &                  eltype,telem,elem,cove,nore,mass,h,g,comele) 
  USE fund_const   
  IMPLICIT NONE 
  INTEGER unit,nobj 
  DOUBLE PRECISION telem(nobj),elem(6,nobj),cove(6,6,nobj) 
  DOUBLE PRECISION nore(6,6,nobj),mass(nobj),h(nobj),g(nobj) 
  CHARACTER*(*) filnam,objnam(nobj),eltype(nobj),comele(nobj) 
  LOGICAL deforb(nobj),defcn(nobj)
  INTEGER,PARAMETER :: nobjx=10                                          
  INTEGER nrem,k,lf,nr,ln,lc 
  DOUBLE PRECISION el1(6) 
  CHARACTER*7 nmpc(nobjx),nmpc1 
  CHARACTER hc*5,ep5*5,krc*10,gc*5 
  LOGICAL error                                                             
  INTEGER lench 
  EXTERNAL lench 
  IF(nobj.GT.nobjx) STOP '**** rdmpca: nobj > nobjx ****' 
! Number of remaining object (orbit not yet found)                      
  nrem=0 
  DO 10 k=1,nobj 
     IF(deforb(k)) GOTO 10 
     nrem=nrem+1 
     CALL mpcpds(objnam(k),nmpc(k),error) 
     IF(error) THEN 
        ln=lench(objnam(k)) 
        WRITE(*,110) objnam(k)(1:ln) 
110     FORMAT('rdmpca: cannot understand asteroid code "',A,'"') 
     END IF
10 ENDDO
  IF(nrem.LE.0) RETURN 
  lf=lench(filnam)
  nr=0 
1 CONTINUE 
  READ(unit,100,END=2) nmpc1 
100 FORMAT(A7,A5,2X,A5,1X,A5,1X,F9.5,2X,F9.5,2X,F9.5,2X,F9.5,         &
     &       2X,F9.7,13X,F11.7)
! handling alphanumeric asteroid number
  nr=nr+1                                                               
  DO 3 k=1,nobj 
     IF(deforb(k)) GOTO 3 
     IF(nmpc1.EQ.nmpc(k)) THEN 
        BACKSPACE(unit) 
        READ(unit,100) nmpc1,hc,gc,ep5,el1 
        deforb(k)=.true. 
        defcn(k)=.false. 
        eltype(k)='KEP' 
        CALL mpcdat(ep5,telem(k),error) 
        IF(error) THEN 
           WRITE(*,111) nr,filnam(1:lf) 
111        FORMAT('rdmpca: illegal date code at record',I6,' of file "',A,'"')
           STOP '**** rdmpca: illegal date code ****' 
        END IF
        elem(1,k)=el1(6) 
        elem(2,k)=el1(5) 
        elem(3,k)=el1(4)*radeg 
        elem(4,k)=el1(3)*radeg 
        elem(5,k)=el1(2)*radeg 
        elem(6,k)=el1(1)*radeg 
        mass(k)=0.d0 
        IF(hc.EQ.'     ') THEN 
           h(k)=-1.D9 
        ELSE 
           READ(hc,101) h(k) 
        END IF
        IF(gc.EQ.'     ') THEN 
           g(k)=-1.D9 
        ELSE 
           READ(gc,101) g(k) 
101        FORMAT(F5.2) 
        END IF
        WRITE(krc,107) nr 
107     FORMAT(I6) 
        CALL rmsp(krc,lc) 
        comele(k)='read from file "'//filnam(1:lf)//'" at record '//krc(1:lc) 
        nrem=nrem-1 
        IF(nrem.LE.0) RETURN 
     END IF
3 END DO                                                                 
  GOTO 1 
2 CONTINUE 
END SUBROUTINE rdmpca
! adapted from rdmpca2 by A. Milani (to be used by mfitobs/ufitobs)
! with the purpose of saving memory space (by not having covariance
! and normal matrices, not supplied by MPC anyway) 
!                                                                       
!  *****************************************************************    
!  *                                                               *    
!  *                         R D M P C A 3                         *    
!  *                                                               *    
!  *          Read orbital elements of all objects                 *    
!  *        as given in a file written in MPC format for asteroids *    
!  *                                                               *    
!  *****************************************************************    
!                                                                       
! INPUT:    UNIT      -  Input FORTRAN unit (must be already opened)    
!           FILNAM    -  Input file name (for error messages)           
!           NOBJX      -  Number of objects (maximum)                   
!                                                                       
! OUTPUT:   OBJNAM    -  Object names                                   
!           DEFCN     -  Tells whether covariance/normal matrices       
!                            are defined                                
!           ELTYPE    -  Type of orbital elements (EQU/KEP/CAR)         
!           TELEM     -  Epoch of orbital elements (MJD, TDT)           
!           ELEM      -  Orbital elements (ECLM J2000)                  
!           COVE      -  Covariance matrix of orbital elements          
!           NORE      -  Normal matrix of orbital elements              
!           MASS      -  Mass (solar masses)                            
!           H         -  H absolute magnitude (if <-100, missing)       
!           G         -  G slope parameter                              
!           COMELE    -  Comment on orbital elements                    
!           NOBJ      -  Number of objects, actual number found    

SUBROUTINE rdmpca3(unit,filnam,objnam,nobjx,nobj,defcn,           &
     &                  eltype,telem,elem,h,g,comele)
  USE fund_const                   
  IMPLICIT NONE
  INTEGER unit,nobj,nobjx 
  DOUBLE PRECISION telem(nobjx),elem(6,nobjx) 
!  DOUBLE PRECISION nore(6,6,nobjx),mass(nobjx),cove(6,6,nobjx) 
  DOUBLE PRECISION h(nobjx),g(nobjx) 
  CHARACTER*(*) filnam,objnam(nobjx),eltype,comele 
  LOGICAL defcn
  INTEGER nrem,k,lf,nr,ln,lc 
  DOUBLE PRECISION el1(6) 
  CHARACTER*7 nmpc1 
  CHARACTER hc*5,ep5*5,krc*10,gc*5 
  LOGICAL error
  INTEGER lench 
  EXTERNAL lench 
  CALL rmsp(filnam,lf)
  nr=0 
  DO 3 k=1,nobjx 
     READ(unit,100,END=2) nmpc1,hc,gc,ep5,el1 
100  FORMAT(A7,A5,2X,A5,1X,A5,1X,F9.5,2X,F9.5,2X,F9.5,2X,F9.5,      &
     &        2X,F9.7,13X,F11.7)                                        
     nr=nr+1 
     nobj=nr 
! name conversion                                                       
     CALL iaucod2(nmpc1,objnam(k),error) 
     IF(error) THEN 
        WRITE(*,112) nr,filnam(1:lf) 
112     FORMAT('INPUT ERROR: illegal name code at record',I6,       &
     &       ' of file "',A,'"')                                        
        STOP '**** rdmpca3: name conversion error ****' 
     END IF
     CALL mpcdat(ep5,telem(k),error) 
     IF(error) THEN 
        WRITE(*,111) nr,filnam(1:lf) 
111     FORMAT('INPUT ERROR: illegal date code at record',I6,       &
     &       ' of file "',A,'"')                                        
        STOP '**** rdmpca3: date conversion error ****' 
     END IF
     elem(1,k)=el1(6) 
     elem(2,k)=el1(5) 
     elem(3,k)=el1(4)*radeg 
     elem(4,k)=el1(3)*radeg 
     elem(5,k)=el1(2)*radeg 
     elem(6,k)=el1(1)*radeg 
!     mass(k)=0.d0 
     IF(hc.EQ.'     ') THEN 
        h(k)=-1.D9 
     ELSE 
        READ(hc,101) h(k) 
101     FORMAT(F5.2) 
     END IF
     IF(gc.EQ.'     ') THEN 
        g(k)=-1.D9 
     ELSE 
        READ(gc,101,ERR=99) g(k) 
     END IF
     GOTO 98
99   WRITE(*,*)objnam(k),h(k), g(k)
     g(k)=-1.D9
98   CONTINUE
     WRITE(krc,107) nr 
107  FORMAT(I6) 
     CALL rmsp(krc,lc)
3 ENDDO
  WRITE(*,*)' file not completely read, record ',nr 
2 CONTINUE 
  comele='read from file "'//filnam(1:lf) 
  defcn=.false. 
  eltype='KEP' 
END SUBROUTINE rdmpca3                           

                                 
! Copyright (C) 1997-1998 by Mario Carpino (carpino@brera.mi.astro.it)  
! Version: June 21, 1998                                                
! --------------------------------------------------------------------- 
!                                                                       
!  *****************************************************************    
!  *                                                               *    
!  *                         O P O R B F                           *    
!  *                                                               *    
!  *                Open an orbital element file                   *    
!  *                                                               *    
!  *****************************************************************    
!                                                                       
! INPUT:    FILE      -  File name                                      
!           INUNIT    -  if it is .le.0, the file must be opened        
!                        if not, it is the unit                         
!                                                                       
 SUBROUTINE oporbf(file,inunit)
   USE io_elems 
   IMPLICIT NONE  
   CHARACTER*(*) file 
   INTEGER, INTENT(IN) :: inunit 
   INTEGER kr,lf,mjd,mjde,nn 
   LOGICAL first,found,end 
   CHARACTER scale*3,rec*200 
   DOUBLE PRECISION sec,sece 
   SAVE first    
   INTEGER lench,nitchs 
   EXTERNAL lench,nitchs                                                       
   DATA first/.true./ 
   IF(first) THEN 
      orbunt=0 
      orbfn=' ' 
      first=.false. 
   END IF      
   IF(orbunt.NE.0) STOP '**** oporbf: error (01), second OEF file ****' 
   orbfn=file 
   IF(inunit.gt.0)THEN 
! no need to open, the file is opened                                   
      orbunt=inunit 
   ELSE 
      CALL filopn(orbunt,orbfn,'OLD') 
   ENDIF
   CALL rdfnam(orbunt,orbfn,orbnr) 
   lf=lench(orbfn)
! Format                                                                
   CALL rdfcha(orbunt,'format',.true.,oef_vers,found,kr) 
   IF(oef_vers.NE.'OEF1.1'.and.oef_vers.NE.'OEF2.0') THEN 
      WRITE(*,100) orbfn(1:lf) 
100   FORMAT(' ERROR: unsupported format in file ',A)
      STOP '**** oporbf: wrong file format ****' 
   END IF
! Record type and default orbital element type                          
   CALL rdfcha(orbunt,'rectype',.false.,rectyp,found,kr) 
   IF(.NOT.found) THEN 
1     CONTINUE 
      CALL getrsc(orbunt,rec,orbnr,end) 
      IF(end) THEN 
         WRITE(*,104) orbfn(1:lf) 
104      FORMAT(' ERROR: file ',A,' is empty') 
         STOP '**** oporbf: abnormal end ****' 
      END IF
      orbnr=orbnr-1 
      BACKSPACE(orbunt) 
      nn=nitchs(rec) 
      IF(nn.EQ.1) THEN 
         rectyp='ML' 
      ELSEIF(nn.GE.7) THEN 
         rectyp='1L' 
      ELSE 
         orbnr=orbnr+1 
         GOTO 10 
      END IF
   END IF
   IF(rectyp.EQ.'1L') THEN 
      CALL rdfcha(orbunt,'elem',.true.,deltyp,found,kr) 
   ELSEIF(rectyp.EQ.'ML') THEN 
      deltyp=' ' 
   ELSE 
      WRITE(*,101) orbfn(1:lf) 
101   FORMAT('ERROR: unsupported record type in file ',A) 
      STOP '**** oporbf: abnormal end ****' 
   END IF                                                             
! Default reference system                                              
   CALL rdfref(orbunt,'refsys',.false.,dfrsty,dfrsep,found,kr) 
   IF(.NOT.found) THEN 
      IF(rectyp.EQ.'1L') THEN 
         WRITE(*,105) orbfn(1:lf) 
         STOP '**** oporbf: abnormal end ****' 
      END IF
      dfrsty=' ' 
      dfrsep=' ' 
   END IF
105 FORMAT(' ERROR: missing keyword "refsys" in file ',A) 
                                                                        
! Default epoch for orbital elements                                    
   CALL rdftim(orbunt,'epoch',.false.,depstr,mjd,sec,scale,deft0,kr) 
   IF(deft0) THEN 
      CALL cnvtim(mjd,sec,scale,mjde,sece,'TDT') 
      dept0=mjde+sece/86400.d0 
   END IF                                                                     
   nxtend=.false. 
   iicorb=36 
   RETURN 
                                                                        
10 CONTINUE 
   WRITE(*,102) orbfn(1:lf),orbnr 
102 FORMAT(' FORMAT ERROR in file ',A,' at line',I5) 
   STOP '**** oporbf: abnormal end ****' 
 END SUBROUTINE oporbf
! Copyright (C) 1997 by Mario Carpino (carpino@brera.mi.astro.it)       
! Version: March 7, 1997                                                
! --------------------------------------------------------------------- 
!                                                                       
!  *****************************************************************    
!  *                                                               *    
!  *                         C L O R B F                           *    
!  *                                                               *    
!  *                Close an orbital element file                  *    
!  *                                                               *    
!  *****************************************************************    
!                                                                       
! INPUT:    FILE      -  File name                                      
!                                                                       
 SUBROUTINE clorbf
   USE io_elems 
   IMPLICIT NONE
   IF(orbunt.LE.0) STOP '**** clorbf: internal error (01) ****' 
   CALL filclo(orbunt,' ') 
   orbunt=0 
   orbfn=' ' 
   iicorb=0
 END SUBROUTINE clorbf
! Copyright (C) 1997-1998 by Mario Carpino (carpino@brera.mi.astro.it)  
! Version: June 21, 1998                                                
! modified to handle COM/ATT A. Milani June 2005                                                
! --------------------------------------------------------------------- 
!                                                                       
!  *****************************************************************    
!  *                                                               *    
!  *                          R D O R B                            *    
!  *                                                               *    
!  *     Read orbital elements from a file opened with OPORBF      *    
!  *                                                               *    
!  *****************************************************************    
!                                                                       
! The routine operates in a sequential way, returning each time         
! the orbital elements of the next object contained in the file         
!                                                                       
! OUTPUT:   NAME      -  Name of planet/asteroid/comet                  
!           ELEM(6)   -  Orbital element vector                         
!           ELTYPE    -  Type of orbital elements (KEP/EQU/CAR/COM/ATT)
!           T0        -  Epoch of orbital elements (MJD, TDT)           
!           COVE      -  Covariance matrix of orbital elements          
!           DEFCOV    -  Tells whether the covariance matrix is defined 
!           NORE      -  Normal matrix of orbital elements              
!           DEFNOR    -  Tells whether the normal matrix is defined     
!           H         -  H absolute magnitude (if <-100, missing)       
!           G         -  G slope parameter                              
!           MASS      -  Mass (solar masses)                            
!           RSYS      -  Reference system type (EQUM/EQUT/ECLM)         
!           EPOCH     -  Epoch specification (J2000/OFDATE)             
!           KR        -  Record number at which object is found         
!           EOF       -  End-of-file flag                               
!                                                                       
 SUBROUTINE rdorb(name,elem,eltype,t0,cove,defcov,nore,defnor,     &
     &                 h,g,mass,rsys,epoch,kr,eof) 
   USE fund_const
   USE io_elems
   IMPLICIT NONE
   DOUBLE PRECISION elem(6),t0,h,g,cove(6,6),nore(6,6),mass 
   CHARACTER*(*) name,eltype,rsys,epoch 
   LOGICAL defcov,defnor,eof 
   INTEGER kr
   CHARACTER rec*200,rest*200,scale*3 
   INTEGER lf,nit,i,k,mjd,mjde,ik 
   DOUBLE PRECISION sec,sece,tmp(21),cnv(6) 
   LOGICAL error,end1,noep
   INTEGER lench,nitchs 
   EXTERNAL lench,nitchs                                                       
   IF(iicorb.NE.36) STOP '**** rdorb: internal error (01) ****' 
   IF(orbunt.LE.0) STOP '**** rdorb: internal error (02) ****' 
   rsys=dfrsty 
   epoch=dfrsep 
   mass=0.d0 
   DO 1 i=1,6 
      DO 2 k=1,6 
         cove(i,k)=0 
         nore(i,k)=0 
2     END DO
1  END DO
   defcov=.false. 
   defnor=.false. 
   h=-1.d9
   g=-1.d9
   IF(rectyp.EQ.'1L') THEN 
      CALL getrsc(orbunt,rec,orbnr,eof) 
      IF(eof) RETURN 
      CALL strcnt(rec,name,rest,error) 
      IF(error) GOTO 20 
      kr=orbnr 
      nit=nitchs(rest) 
      IF(deft0) THEN 
         t0=dept0 
         IF(nit.LT.8) THEN 
            READ(rest,*,ERR=20) elem 
         ELSE 
            READ(rest,*,ERR=20) elem,h,g 
         END IF
      ELSE 
         IF(nit.LT.9) THEN 
            READ(rest,*,ERR=20) t0,elem 
         ELSE 
            READ(rest,*,ERR=20) t0,elem,h,g 
         END IF
      END IF
      eltype=deltyp 
   ELSEIF(rectyp.EQ.'ML') THEN 
      IF(nxtend) THEN 
         eof=.true. 
         RETURN 
      END IF
      noep=.true. 
! Name                                                                  
      CALL getrsc(orbunt,name,orbnr,eof) 
      IF(eof) RETURN 
      kr=orbnr 
! Orbital elements (mandatory, immediately after the name)              
      CALL getrsc(orbunt,rec,orbnr,end1) 
      IF(end1) GOTO 20
      IF(rec(1:1).ne.' ') GOTO 20 
      eltype=rec(2:4)
      READ(rec(5:),*,ERR=20) elem
! Other keywords                                                        
3     CONTINUE 
      CALL getrsc(orbunt,rec,orbnr,end1) 
      IF(end1) THEN 
         nxtend=.true. 
         GOTO 4 
      END IF
      IF(rec(1:1).NE.' ') THEN 
         BACKSPACE(orbunt) 
         orbnr=orbnr-1 
         GOTO 4 
      END IF
! Epoch of elements                                                     
      IF(rec(1:4).EQ.' MJD' .OR. rec(1:4).EQ.' JD ' .OR.            &
           &       rec(1:4).EQ.' CAL') THEN                                   
         CALL ch2tim(rec,mjd,sec,scale,error) 
         IF(error) GOTO 20 
         CALL cnvtim(mjd,sec,scale,mjde,sece,'TDT') 
         t0=mjde+sece/86400.d0 
         noep=.false. 
      ELSEIF(rec(1:4).EQ.' MAG') THEN 
         READ(rec(5:),*,ERR=20) h,g 
      ELSEIF(rec(1:4).EQ.' STA')THEN ! station code for ATT only
         READ(rec(5:),*,ERR=20) obscod
      ELSEIF(rec(1:4).EQ.' MAS') THEN 
         READ(rec(5:),*,ERR=20) mass 
      ELSEIF(rec(1:4).EQ.' COV') THEN 
         READ(rec(5:),*,ERR=40) (tmp(i),i=1,3) 
         DO 17 k=1,6 
            CALL getrsc(orbunt,rec,orbnr,end1) 
            IF(end1) GOTO 40 
            IF(rec(1:4).NE.' COV') GOTO 40 
            READ(rec(5:),*,ERR=20) (tmp(i),i=3*k+1,3*k+3) 
17       CONTINUE 
         ik=0 
         DO 8 i=1,6 
            DO 7 k=i,6 
               ik=ik+1 
               cove(i,k)=tmp(ik) 
7           ENDDO 
8        ENDDO 
         IF(ik.NE.21) STOP '**** rdorb: internal error (03) ****' 
         defcov=.true. 
      ELSEIF(rec(1:4).EQ.' NOR') THEN 
         READ(rec(5:),*,ERR=40) (tmp(i),i=1,3) 
         DO 27 k=1,6 
            CALL getrsc(orbunt,rec,orbnr,end1) 
            IF(end1) GOTO 40 
            IF(rec(1:4).NE.' NOR') GOTO 40 
            READ(rec(5:),*,ERR=40) (tmp(i),i=3*k+1,3*k+3) 
   27    ENDDO 
         ik=0 
         DO 38 i=1,6 
            DO 37 k=i,6 
               ik=ik+1 
               nore(i,k)=tmp(ik) 
37          ENDDO
38       ENDDO
         IF(ik.NE.21) STOP '**** rdorb: internal error (04) ****' 
         defnor=.true. 
      ELSE 
         GOTO 40 
      END IF 
      GOTO 3 

40    CONTINUE
      lf=lench(orbfn) 
      WRITE(*,240) orbfn(1:lf),orbnr 
240   FORMAT(' ERROR in covariance, file ',A,' at line',I6) 

4     CONTINUE 
      IF(noep) THEN 
         IF(deft0) THEN 
            t0=dept0 
         ELSE 
            GOTO 20 
         END IF
      END IF
   ELSE 
      STOP '**** rdorb: internal error (05) ****' 
   END IF                                                                     
! Transformation of angles in orbital elements                          
! and covariance matrix                                                 
   cnv(1:6)=1 
   IF(eltype.EQ.'KEP') THEN 
      cnv(3:6)=radeg 
   ELSEIF(eltype.EQ.'COM') THEN 
      cnv(3:5)=radeg 
   ELSEIF(eltype.EQ.'COT')THEN
      cnv(3:6)=radeg
   ELSEIF(eltype.EQ.'EQU') THEN 
      cnv(6)=radeg 
   ELSEIF(eltype.EQ.'CAR') THEN 
      CONTINUE 
   ELSEIF(eltype.EQ.'ATT') THEN 
      cnv(1:4)=radeg 
   ELSE 
      STOP '**** rdorb: internal error (06) ****' 
   END IF
   DO  i=1,6 
      elem(i)=elem(i)*cnv(i) 
   END DO 
   IF(defcov) THEN 
      DO 13 i=1,6 
         cove(i,i)=cove(i,i)*(cnv(i)**2) 
         DO 12 k=i+1,6 
            cove(i,k)=cove(i,k)*cnv(i)*cnv(k) 
            cove(k,i)=cove(i,k) 
12       ENDDO
13    ENDDO
   END IF
   IF(defnor) THEN 
      DO 15 i=1,6 
         nore(i,i)=nore(i,i)/(cnv(i)**2) 
         DO 14 k=i+1,6 
            nore(i,k)=nore(i,k)/(cnv(i)*cnv(k)) 
            nore(k,i)=nore(i,k) 
14       ENDDO
15    ENDDO
   END IF                                                                
   eof=.false. 
   RETURN                                                                      
20 CONTINUE 
   lf=lench(orbfn) 
   WRITE(*,200) orbfn(1:lf),orbnr 
200 FORMAT(' ERROR in file ',A,' at line',I6) 
   STOP '**** rdorb: abnormal end ****' 
END SUBROUTINE rdorb
! ===============================================================
! OUTPUT ELEMENTS
! ===============================================================
! Copyright (C) 1997-1998 by Mario Carpino (carpino@brera.mi.astro.it)  
! Version: June 21, 1998                                                
! --------------------------------------------------------------------- 
!                                                                       
!  *****************************************************************    
!  *                                                               *    
!  *                         W R O 1 L H                           *    
!  *                                                               *    
!  *       Writes the header of an orbital element file            *    
!  *   (single-line format, different epochs, keplerian elements)  *    
!  *                                                               *    
!  *****************************************************************    
!                                                                       
! INPUT:   UNIT      -  Output FORTRAN unit                            
!           RSYS      -  Reference system type (EQUM/EQUT/ECLM)         
!           EPOCH     -  Reference system epoch (J2000/OFDATE)          
!           ELTYPE    -  Type of orbital elements (EQU/KEP/CAR)         
!                                                                       
SUBROUTINE wro1lh(unit,rsys,epoch,eltype) 
  IMPLICIT NONE
  INTEGER unit 
  CHARACTER*(*) rsys,epoch,eltype
  INCLUDE 'parcmc.h90'                                                  
  INTEGER l1,l2,l3,nb 
  CHARACTER*100 bl
  INTEGER lench 
  EXTERNAL lench                                                          
  l1=lench(rsys) 
  l2=lench(epoch) 
  l3=lench(eltype) 
  nb=MAX(14-l1-l2,1) 
  bl=' '                                                                
  WRITE(unit,100) comcha,comcha,eltype(1:l3),comcha,                &
     &                rsys(1:l1),epoch(1:l2),                           &
     &                bl(1:nb),comcha                                   
100 FORMAT('format  = ''OEF2.0''       ',A,' file format'/            &
     &       'rectype = ''1L''           ',A,' record type (1L/ML)'/    &
     &       'elem    = ''',A,'''          ',A,                         &
     &                   ' type of orbital elements'/                   &
     &       'refsys  = ',A,1X,A,A,A,' default reference system'/       &
     &       'END_OF_HEADER')                                           
                                                                        
  IF(eltype.EQ.'KEP') THEN 
     WRITE(unit,201) comcha 
201  FORMAT(A,' Name, Epoch(MJD), a, e, i, long. node,',' arg. peric., mean anomaly') 
  ELSEIF(eltype.EQ.'CAR') THEN 
     WRITE(unit,202) comcha 
202  FORMAT(A,' Name, Epoch(MJD), cartesian position and velocity',' vectors')
  ELSEIF(eltype.EQ.'EQU') THEN 
     WRITE(unit,203) comcha 
203  FORMAT(A,' Name, Epoch(MJD), a, e*sin(LP), e*cos(LP),',' tan(i/2)*sin(LN), tan(i/2)*cos(LN), mean long.')
  ELSEIF(eltype.EQ.'COM') THEN
     WRITE(unit,204) comcha 
204  FORMAT(A,' Name, Epoch(MJD), q, e, i, long. node,',' arg. peric., perihelion time')
  ELSEIF(eltype.EQ.'COT') THEN
     WRITE(unit,205) comcha 
205  FORMAT(A,' Name, Epoch(MJD), q, e, i, long. node,',' arg. peric., true anomaly')
  ELSEIF(eltype.eq.'ATT')THEN
     WRITE(unit,206) comcha 
206  FORMAT(A,' Name, Epoch(MJD), R.A., DEC, R.A.dot, DECdot, r, rdot,')
  END IF
END SUBROUTINE wro1lh

!  *****************************************************************    
!  *                                                               *    
!  *                         W R O 1 L H _ M A T L A B             *    
!  *                                                               *    
!  *       Writes the header of an orbital element file            *    
!  *   (single-line format, different epochs, keplerian elements)  *    
!  *                                                               *    
!  *****************************************************************    
!                                                                       
! INPUT:   UNIT      -  Output FORTRAN unit                            
!           RSYS      -  Reference system type (EQUM/EQUT/ECLM)         
!           EPOCH     -  Reference system epoch (J2000/OFDATE)          
!           ELTYPE    -  Type of orbital elements (EQU/KEP/CAR)         
!                                                                       
SUBROUTINE wro1lh_matlab(unit,rsys,epoch,eltype) 
  IMPLICIT NONE
  INTEGER unit 
  CHARACTER*(*) rsys,epoch,eltype
  INCLUDE 'parcmc.h90'                                                  
  INTEGER l1,l2,l3,nb 
  CHARACTER*100 bl
  INTEGER lench 
  EXTERNAL lench                                                          
  l1=lench(rsys) 
  l2=lench(epoch) 
  l3=lench(eltype) 
  nb=MAX(14-l1-l2,1) 
  bl=' '                                                                
  WRITE(unit,100) comcha,comcha,eltype(1:l3),comcha,                &
     &                rsys(1:l1),epoch(1:l2),                           &
     &                bl(1:nb),comcha                                   
100 FORMAT('%format  = ''OEF2.0''       ',A,' file format'/            &
     &       '%rectype = ''1L''           ',A,' record type (1L/ML)'/    &
     &       '%elem    = ''',A,'''          ',A,                         &
     &                   ' type of orbital elements'/                   &
     &       '%refsys  = ',A,1X,A,A,A,' default reference system'/       &
     &       '%END_OF_HEADER')                                           
                                                                        
  IF(eltype.EQ.'KEP') THEN 
     WRITE(unit,201) comcha 
201  FORMAT('%',A,' Name, Epoch(MJD), a, e, i, long. node,',' arg. peric., mean anomaly, H mag, NGR, chi, RMS, &
          &connected component') 
  ELSEIF(eltype.EQ.'CAR') THEN 
     WRITE(unit,202) comcha 
202  FORMAT('%',A,' Name, Epoch(MJD), cartesian position and velocity',' vectors')
  ELSEIF(eltype.EQ.'EQU') THEN 
     WRITE(unit,203) comcha 
203  FORMAT('%',A,' Name, Epoch(MJD), a, e*sin(LP), e*cos(LP),',' tan(i/2)*sin(LN), tan(i/2)*cos(LN), mean long.')
  ELSEIF(eltype.EQ.'COM') THEN
     WRITE(unit,204) comcha 
204  FORMAT('%',A,'Name       Epoch(MJD)     a (AU)                   q (AU)                   e                        &
          &i (deg)                  long. node (deg)         arg. peric. (deg)       peric. time  &
          &H      NGR  chi          RMS         conn. comp.')
  ELSEIF(eltype.EQ.'COT') THEN
     WRITE(unit,205) comcha 
205  FORMAT('%',A,' Name, Epoch(MJD), q, e, i, long. node,',' arg. peric., true anomaly')
  ELSEIF(eltype.eq.'ATT')THEN
     WRITE(unit,206) comcha 
206  FORMAT('%',A,' Name, Epoch(MJD), R.A., DEC, R.A.dot, DECdot, r, rdot,')
  END IF
END SUBROUTINE wro1lh_matlab

! Version May 20, 2015 Federica Spoto (spoto@mail.dm.unipi.it)
! 1 line catalogs for non-grav parameters
! The catalog contains:
! - the name of the asteroid (or VA)
! - the dynamical model used
! - the number of parameters in use
! - the values of the parameters
! --------------------------------------------------------------------- 
!                                                                       
!  *****************************************************************    
!  *                                                               *    
!  *                 W R O 1 L R _ N G R                           *    
!  *                                                               *    
!  *  Writes an orbital element record in an orbital element file  *    
!  *   (single-line format, different epochs, keplerian elements)  *    
!  *                     (multi-line format)                       *    
!  *                                                               *    
!  *****************************************************************    
!                                                                       
! WARNING: the routine does not write the header of the file: this      
!          must be generated by calling subroutine wro1lh2_ngr               
!                                                                       

SUBROUTINE wro1lr_ngr(unit,name,n_mod,n_dp,ngr_par)
  USE fund_const
  IMPLICIT NONE
! Begin interface
  INTEGER,          INTENT(IN) :: unit 
  CHARACTER*(*),    INTENT(IN) :: name
  INTEGER,          INTENT(IN) :: n_mod, n_dp
  DOUBLE PRECISION, INTENT(IN) :: ngr_par(n_dp)
! End interface
! Expected max length of name                                           
  INTEGER, PARAMETER :: namtl=12 
  INTEGER            ::ln,nb,i 
  CHARACTER*(namtl)  :: blanks
  INTEGER            ::lench 
  EXTERNAL lench 
!============================================================
! Name                                                                  
  ln=MAX(1,lench(name)) 
  nb=MAX(1,namtl-ln) 
  blanks=' ' 
  IF(n_dp.EQ.2) THEN
     WRITE(unit,600) name(1:ln),blanks(1:nb),n_mod,n_dp,ngr_par
600  FORMAT('''',A,'''',A,2X,I2,3X,I2,1P,2E22.14)
  ELSE IF(n_dp.EQ.3) THEN
     WRITE(unit,601) name(1:ln),blanks(1:nb),n_mod,n_dp,ngr_par
601  FORMAT('''',A,'''',A,2X,I2,3X,I2,1P,3E22.14)
  ELSE IF(n_dp.EQ.4) THEN
     WRITE(unit,602) name(1:ln),blanks(1:nb),n_mod,n_dp,ngr_par
602  FORMAT('''',A,'''',A,2X,I2,3X,I2,1P,4E22.14)
  ELSE IF(n_dp.EQ.1) THEN
     WRITE(unit,603) name(1:ln),blanks(1:nb),n_mod,n_dp,ngr_par
603  FORMAT('''',A,'''',A,2X,I2,3X,I2,1P,E22.14)
  ELSE
     WRITE(*,*) 'wro1lr_ngr: n_dp, wrong value ', n_dp
     STOP
  END IF
END SUBROUTINE wro1lr_ngr


! Copyright (C) 1997-2000 by Mario Carpino (carpino@brera.mi.astro.it)  
! Version: June 7, 2000                                                 
! --------------------------------------------------------------------- 
!                                                                       
!  *****************************************************************    
!  *                                                               *    
!  *                         O U T E L E                           *    
!  *                                                               *    
!  *           Verbose output of a set orbital elements            *    
!  *                       to a report file                        *    
!  *                                                               *    
!  *****************************************************************    
!                                                                       
! INPUT:    UNIT      -  Output FORTRAN unit (report file)              
!           ELEM      -  Orbital elements (ECLM J2000)                  
!           ELTYPE    -  Type of orbital elements (EQU/KEP/CAR)         
!           T0        -  Epoch of orbital elements (MJD, TDT)           
!           LABEL     -  Label                                          
!           MULTI     -  Multi-line output                              
!           STDOUT    -  Standard output                                
!                                                                       
SUBROUTINE outele(unit,elem,eltype,t0,label,multi,stdout) 
  USE fund_const
  IMPLICIT NONE                                                             

  INTEGER unit 
  DOUBLE PRECISION elem(6),t0 
  CHARACTER*(*) eltype,label 
  LOGICAL multi,stdout
  INTEGER lt,i,day,month,year,ll 
  CHARACTER cm*3 
  DOUBLE PRECISION hour
  INTEGER lench 
  CHARACTER*3 chmon 

  EXTERNAL lench,chmon
  ll=lench(label)                                                      
  IF(multi .AND. (ll.GT.0)) THEN 
     IF(unit.GT.0) WRITE(unit,133) label(1:ll) 
     IF(stdout) WRITE(*,133) label(1:ll) 
133  FORMAT(8X,'Orbital elements for ',A,':') 
  END IF
  IF(eltype.EQ.'KEP') THEN 
     IF(multi) THEN 
        IF(unit.GT.0) WRITE(unit,100) elem(1),elem(2),(elem(i)*degrad,i=3,6)
        IF(stdout) WRITE(*,100) elem(1),elem(2),(elem(i)*degrad,i=3,6)   
100     FORMAT(8X,'Semimajor axis     =',1P,E24.16,0P,' au'/          &
             &       8X,'Eccentricity       =',F20.15/                    &
             &       8X,'Inclination        =',F18.13,' deg'/             &
             &       8X,'Long. of node      =',F18.13,' deg'/             &
             &       8X,'Arg. of pericenter =',F18.13,' deg'/             &
             &       8X,'Mean anomaly       =',F18.13,' deg')  
     ELSE 
        IF(ll.GT.0) THEN 
           IF(unit.GT.0) WRITE(unit,120) label(1:ll),elem(1),elem(2),(elem(i)*degrad,i=3,6),t0
           IF(stdout) WRITE(*,120) label(1:ll),elem(1),elem(2),(elem(i)*degrad,i=3,6),t0     
        ELSE 
           IF(unit.GT.0) WRITE(unit,130) elem(1),elem(2),(elem(i)*degrad,i=3,6),t0
           IF(stdout) WRITE(unit,130) elem(1),elem(2),(elem(i)*degrad,i=3,6),t0
120        FORMAT(8X,'KepElem(',A,'):',1P,E15.7,0P,F13.8,4F10.5,' (T=',F10.3,')')
130        FORMAT(8X,'KepElem:',1P,E15.7,0P,F13.8,4F10.5,' (T=',F10.3,')')   
        END IF
     END IF
  ELSEIF(eltype.EQ.'COM') THEN 
     IF(multi) THEN 
        IF(unit.GT.0) WRITE(unit,104) elem(1),elem(2),(elem(i)*degrad,i=3,5),elem(6)      
        IF(stdout) WRITE(*,104) elem(1),elem(2),(elem(i)*degrad,i=3,5),elem(6) 
104     FORMAT(8X,'Pericenter distance  =',1P,E23.14,0P,' au'/        &
             &       8X,'Eccentricity       =',F20.15/                    &
             &       8X,'Inclination        =',F18.13,' deg'/             &
             &       8X,'Long. of node      =',F18.13,' deg'/             &
             &       8X,'Arg. of pericenter =',F18.13,' deg'/             &
             &       8X,'Time of pericenter =',F18.13,' MJD')  
     ELSE 
        IF(ll.GT.0) THEN 
           IF(unit.GT.0) WRITE(unit,124) label(1:ll),elem(1),elem(2),(elem(i)*degrad,i=3,6),t0
           IF(stdout) WRITE(*,124) label(1:ll),elem(1),elem(2),(elem(i)*degrad,i=3,6),t0     
        ELSE 
           IF(unit.GT.0) WRITE(unit,134) elem(1),elem(2),(elem(i)*degrad,i=3,6),t0
           IF(stdout) WRITE(unit,134) elem(1),elem(2),(elem(i)*degrad,i=3,6),t0 
124        FORMAT(8X,'ComElem(',A,'):',1P,E15.7,0P,F13.8,3F10.5,F13.6,' (T=',F10.3,')')
134        FORMAT(8X,'ComElem:',1P,E15.7,0P,F13.8,3F10.5,F13.6,' (T=',F10.3,')')   
        END IF
     END IF
  ELSEIF(eltype.EQ.'EQU') THEN 
     IF(multi) THEN 
        IF(unit.GT.0) WRITE(unit,101) (elem(i),i=1,5),elem(6)*degrad   
        IF(stdout) WRITE(*,101) (elem(i),i=1,5),elem(6)*degrad
101     FORMAT(8X,'Semimajor axis     =',1P,E23.14,0P,' au'/     &
             &       8X,'h [e*sin(w)]       =',F20.15/                   &
             &       8X,'k [e*cos(w)]       =',F20.15/                   &
             &       8X,'P [tg(i/2)*sin(N)] =',F20.15/                   &
             &       8X,'Q [tg(i/2)*cos(N)] =',F20.15/                   &
             &       8X,'Mean longitude     =',F18.13,' deg')   
     ELSE 
        IF(ll.GT.0) THEN 
           IF(unit.GT.0) WRITE(unit,121) label(1:ll),(elem(i),i=1,5),elem(6)*degrad,t0       
           IF(stdout) WRITE(*,121) label(1:ll),(elem(i),i=1,5),elem(6)*degrad,t0             
        ELSE 
           IF(unit.GT.0) WRITE(unit,131) (elem(i),i=1,5),elem(6)*degrad,t0  
           IF(stdout) WRITE(*,131) (elem(i),i=1,5),elem(6)*degrad,t0
121        FORMAT(8X,'EQUElem(',A,'):',1P,E15.7,0P,4F13.8,F10.5,' (T=',F10.3,')')
131        FORMAT(8X,'EQUElem:',1P,E15.7,0P,4F13.8,F10.5,' (T=',F10.3,')')      
        END IF
     END IF
  ELSEIF(eltype.EQ.'CAR') THEN 
     IF(multi) THEN 
        IF(unit.GT.0) WRITE(unit,102) elem 
        IF(stdout) WRITE(*,102) elem 
102     FORMAT(8X,'Position vector   =',1X,1P,3E22.14,' au'/              &
             &       8X,'Velocity vector   =',1X,3E22.14,' au/d') 
     ELSE 
        IF(ll.GT.0) THEN 
           IF(unit.GT.0) WRITE(unit,122) label(1:ll),elem,t0 
           IF(stdout) WRITE(*,122) label(1:ll),elem,t0 
        ELSE 
           IF(unit.GT.0) WRITE(unit,132) elem,t0 
           IF(stdout) WRITE(*,132) elem,t0 
122        FORMAT(8X,'PosVel(',A,'):',1P,6E15.7,' (T=',F10.3,')') 
132        FORMAT(8X,'PosVel:',1P,6E15.7,' (T=',F10.3,')')  
        END IF
     END IF
  ELSE 
     lt=lench(eltype) 
     WRITE(*,200) eltype(1:lt) 
200  FORMAT('ERROR: unknown type "',A,'" of orbital elements') 
     STOP '**** outele: unknown type of orbital elements ****' 
  END IF
! Write non-grav parameters (if any)
!  IF(n_mod.GT.0.AND.multi.AND.PRESENT(ngropt).AND..NOT.ngropt)THEN
!  IF(n_mod.GT.0.AND.multi)THEN
!     IF(unit.GT.0) THEN
!        WRITE(unit,140) ngr_par(1:n_dp)
!     END IF
!     IF(stdout) THEN
!        WRITE(unit,140) ngr_par(1:n_dp)
!     END IF
!140  FORMAT('Non-gravitational perturbations: '/   &
!          & 8X, 'A/M  = ',1X,1P,E22.14,' m^2/ton'/ &
!          & 8X, ' A2  = ',1X,1P,E22.14,' 10^(-10) au/d^2')
!  END IF
  IF(multi) THEN 
     CALL mjddat(t0,day,month,year,hour) 
     cm=chmon(month) 
     WRITE(unit,110) t0,cm,day,year,hour 
     IF(stdout) WRITE(*,110) t0,cm,day,year,hour 
110  FORMAT(8X,'Epoch of elements  : MJD',F17.8,' TDT (',A,I3,',',I5,',',F10.6,' h)')  
  END IF
END SUBROUTINE outele

! Copyright (C) 1997-2000 by Mario Carpino (carpino@brera.mi.astro.it)  
! Version: June 7, 2000                                                 
! --------------------------------------------------------------------- 
!                                                                       
!  *****************************************************************    
!  *                                                               *    
!  *                         O U T E L E                           *    
!  *                                                               *    
!  *           Verbose output of a set orbital elements            *    
!  *                       to a report file                        *    
!  *                                                               *    
!  *****************************************************************    
!                                                                       
! INPUT:    UNIT      -  Output FORTRAN unit (report file)              
!           ELEM      -  Orbital elements (ECLM J2000)                  
!           ELTYPE    -  Type of orbital elements (EQU/KEP/CAR)         
!           T0        -  Epoch of orbital elements (MJD, TDT)           
!           LABEL     -  Label                                          
!           MULTI     -  Multi-line output                              
!           STDOUT    -  Standard output                                
!                                                                       
SUBROUTINE outele_ngr(unit,elem,eltype,t0,label,multi,stdout,n_mod,n_dp,ngr_par,ngropt) 
  USE fund_const
  IMPLICIT NONE                                                             

  INTEGER unit 
  DOUBLE PRECISION elem(6),t0 
  CHARACTER*(*) eltype,label 
  LOGICAL multi,stdout
  INTEGER lt,i,day,month,year,ll 
  CHARACTER cm*3 
  DOUBLE PRECISION hour
  INTEGER lench 
  CHARACTER*3 chmon 
! Non-gravitational parameters
  INTEGER,          INTENT(IN) :: n_mod, n_dp
  DOUBLE PRECISION, INTENT(IN) :: ngr_par(n_dp)  
  LOGICAL,          INTENT(IN) :: ngropt
  EXTERNAL lench,chmon
  ll=lench(label)                                                      
  IF(multi .AND. (ll.GT.0)) THEN 
     IF(unit.GT.0) WRITE(unit,133) label(1:ll) 
     IF(stdout) WRITE(*,133) label(1:ll) 
133  FORMAT(8X,'Orbital elements for ',A,':') 
  END IF
  IF(eltype.EQ.'KEP') THEN 
     IF(multi) THEN 
        IF(unit.GT.0) WRITE(unit,100) elem(1),elem(2),(elem(i)*degrad,i=3,6)
        IF(stdout) WRITE(*,100) elem(1),elem(2),(elem(i)*degrad,i=3,6)   
100     FORMAT(8X,'Semimajor axis     =',1P,E24.16,0P,' au'/          &
             &       8X,'Eccentricity       =',F20.15/                    &
             &       8X,'Inclination        =',F18.13,' deg'/             &
             &       8X,'Long. of node      =',F18.13,' deg'/             &
             &       8X,'Arg. of pericenter =',F18.13,' deg'/             &
             &       8X,'Mean anomaly       =',F18.13,' deg')  
     ELSE 
        IF(ll.GT.0) THEN 
           IF(unit.GT.0) WRITE(unit,120) label(1:ll),elem(1),elem(2),(elem(i)*degrad,i=3,6),t0
           IF(stdout) WRITE(*,120) label(1:ll),elem(1),elem(2),(elem(i)*degrad,i=3,6),t0     
        ELSE 
           IF(unit.GT.0) WRITE(unit,130) elem(1),elem(2),(elem(i)*degrad,i=3,6),t0
           IF(stdout) WRITE(unit,130) elem(1),elem(2),(elem(i)*degrad,i=3,6),t0
120        FORMAT(8X,'KepElem(',A,'):',1P,E15.7,0P,F13.8,4F10.5,' (T=',F10.3,')')
130        FORMAT(8X,'KepElem:',1P,E15.7,0P,F13.8,4F10.5,' (T=',F10.3,')')   
        END IF
     END IF
  ELSEIF(eltype.EQ.'COM') THEN 
     IF(multi) THEN 
        IF(unit.GT.0) WRITE(unit,104) elem(1),elem(2),(elem(i)*degrad,i=3,5),elem(6)      
        IF(stdout) WRITE(*,104) elem(1),elem(2),(elem(i)*degrad,i=3,5),elem(6) 
104     FORMAT(8X,'Pericenter distance  =',1P,E23.14,0P,' au'/        &
             &       8X,'Eccentricity       =',F20.15/                    &
             &       8X,'Inclination        =',F18.13,' deg'/             &
             &       8X,'Long. of node      =',F18.13,' deg'/             &
             &       8X,'Arg. of pericenter =',F18.13,' deg'/             &
             &       8X,'Time of pericenter =',F18.13,' MJD')  
     ELSE 
        IF(ll.GT.0) THEN 
           IF(unit.GT.0) WRITE(unit,124) label(1:ll),elem(1),elem(2),(elem(i)*degrad,i=3,6),t0
           IF(stdout) WRITE(*,124) label(1:ll),elem(1),elem(2),(elem(i)*degrad,i=3,6),t0     
        ELSE 
           IF(unit.GT.0) WRITE(unit,134) elem(1),elem(2),(elem(i)*degrad,i=3,6),t0
           IF(stdout) WRITE(unit,134) elem(1),elem(2),(elem(i)*degrad,i=3,6),t0 
124        FORMAT(8X,'ComElem(',A,'):',1P,E15.7,0P,F13.8,3F10.5,F13.6,' (T=',F10.3,')')
134        FORMAT(8X,'ComElem:',1P,E15.7,0P,F13.8,3F10.5,F13.6,' (T=',F10.3,')')   
        END IF
     END IF
  ELSEIF(eltype.EQ.'EQU') THEN 
     IF(multi) THEN 
        IF(unit.GT.0) WRITE(unit,101) (elem(i),i=1,5),elem(6)*degrad   
        IF(stdout) WRITE(*,101) (elem(i),i=1,5),elem(6)*degrad
101     FORMAT(8X,'Semimajor axis     =',1P,E23.14,0P,' au'/     &
             &       8X,'h [e*sin(w)]       =',F20.15/                   &
             &       8X,'k [e*cos(w)]       =',F20.15/                   &
             &       8X,'P [tg(i/2)*sin(N)] =',F20.15/                   &
             &       8X,'Q [tg(i/2)*cos(N)] =',F20.15/                   &
             &       8X,'Mean longitude     =',F18.13,' deg')   
     ELSE 
        IF(ll.GT.0) THEN 
           IF(unit.GT.0) WRITE(unit,121) label(1:ll),(elem(i),i=1,5),elem(6)*degrad,t0       
           IF(stdout) WRITE(*,121) label(1:ll),(elem(i),i=1,5),elem(6)*degrad,t0             
        ELSE 
           IF(unit.GT.0) WRITE(unit,131) (elem(i),i=1,5),elem(6)*degrad,t0  
           IF(stdout) WRITE(*,131) (elem(i),i=1,5),elem(6)*degrad,t0
121        FORMAT(8X,'EQUElem(',A,'):',1P,E15.7,0P,4F13.8,F10.5,' (T=',F10.3,')')
131        FORMAT(8X,'EQUElem:',1P,E15.7,0P,4F13.8,F10.5,' (T=',F10.3,')')      
        END IF
     END IF
  ELSEIF(eltype.EQ.'CAR') THEN 
     IF(multi) THEN 
        IF(unit.GT.0) WRITE(unit,102) elem 
        IF(stdout) WRITE(*,102) elem 
102     FORMAT(8X,'Position vector   =',1X,1P,3E22.14,' au'/              &
             &       8X,'Velocity vector   =',1X,3E22.14,' au/d') 
     ELSE 
        IF(ll.GT.0) THEN 
           IF(unit.GT.0) WRITE(unit,122) label(1:ll),elem,t0 
           IF(stdout) WRITE(*,122) label(1:ll),elem,t0 
        ELSE 
           IF(unit.GT.0) WRITE(unit,132) elem,t0 
           IF(stdout) WRITE(*,132) elem,t0 
122        FORMAT(8X,'PosVel(',A,'):',1P,6E15.7,' (T=',F10.3,')') 
132        FORMAT(8X,'PosVel:',1P,6E15.7,' (T=',F10.3,')')  
        END IF
     END IF
  ELSE 
     lt=lench(eltype) 
     WRITE(*,200) eltype(1:lt) 
200  FORMAT('ERROR: unknown type "',A,'" of orbital elements') 
     STOP '**** outele: unknown type of orbital elements ****' 
  END IF
! Write non-grav parameters (if any)
  IF(n_mod.GT.0.AND..NOT.ngropt)THEN
     IF(unit.GT.0) THEN
        WRITE(unit,140) ngr_par(1:n_dp)
     END IF
     IF(stdout) THEN
        WRITE(unit,140) ngr_par(1:n_dp)
     END IF
140  FORMAT('Non-gravitational perturbations: '/   &
          & 8X, 'A/M  = ',1X,1P,E22.14,' m^2/ton'/ &
          & 8X, ' A2  = ',1X,1P,E22.14,' 10^(-10) au/d^2')
  END IF
  IF(multi) THEN 
     CALL mjddat(t0,day,month,year,hour) 
     cm=chmon(month) 
     WRITE(unit,110) t0,cm,day,year,hour 
     IF(stdout) WRITE(*,110) t0,cm,day,year,hour 
110  FORMAT(8X,'Epoch of elements  : MJD',F17.8,' TDT (',A,I3,',',I5,',',F10.6,' h)')  
  END IF
END SUBROUTINE outele_ngr

! Copyright (C) 1997-1998 by Mario Carpino (carpino@brera.mi.astro.it)  
! Version: June 21, 1998                                                
! --------------------------------------------------------------------- 
!                                                                       
!  *****************************************************************    
!  *                                                               *    
!  *                         W R O M L H                           *    
!  *                                                               *    
!  *       Writes the header of an orbital element file            *    
!  *                     (multi-line format)                       *    
!  *                                                               *    
!  *****************************************************************    
!                                                                       
! OUTPUT:   UNIT      -  Output FORTRAN unit                            
!           RSYS      -  Reference system type (EQUM/EQUT/ECLM)         
!           EPOCH     -  Reference system epoch (J2000/OFDATE)          
!                                                                       
SUBROUTINE wromlh(unit,rsys,epoch) 
  IMPLICIT NONE
  INTEGER unit 
  CHARACTER*(*) rsys,epoch
  INCLUDE 'parcmc.h90'                                                  
  INTEGER l1,l2,nb 
  CHARACTER*100 bl                                                          
  INTEGER lench 
  EXTERNAL lench                                                              
  l1=lench(rsys) 
  l2=lench(epoch) 
  nb=MAX(14-l1-l2,1) 
  bl=' '                                                                       
  WRITE(unit,100) comcha,comcha,rsys(1:l1),epoch(1:l2),             &
     &                bl(1:nb),comcha                                   
100 FORMAT('format  = ''OEF2.0''       ',A,' file format'/            &
     &       'rectype = ''ML''           ',A,' record type (1L/ML)'/    &
     &       'refsys  = ',A,1X,A,A,A,' default reference system'/       &
     &       'END_OF_HEADER')   
END SUBROUTINE wromlh

! ===============================================================
! CONVERSIONS
! ===============================================================
! Copyright (C) 1997 by Mario Carpino (carpino@brera.mi.astro.it)       

! Copyright (C) 1997 by Mario Carpino (carpino@brera.mi.astro.it)       
! Version: December 10, 1997                                            
! --------------------------------------------------------------------- 
!                                                                       
!  *****************************************************************    
!  *                                                               *    
!  *                         M P C D A T                           *    
!  *                                                               *    
!  *         Computes MJD from MPC-style packed dates              *    
!  *                                                               *    
!  *****************************************************************    
!                                                                       
! INPUT:    PDATE     -  MPC-style packed date                          
!                                                                       
! OUTPUT:   TJM       -  Modified Julian Date (TDT)                     
!           ERROR     -  Error flag (cannot understand input)           
!                                                                       
SUBROUTINE mpcdat(pdate,tjm,error) 
  IMPLICIT NONE                                                                       
  CHARACTER*(*) pdate 
  DOUBLE PRECISION tjm 
  LOGICAL error                                                                     
  INTEGER year,yy,month,day                                                             
  INTEGER lench 
  LOGICAL isnum 
  DOUBLE PRECISION tjm1 
  EXTERNAL lench,isnum,tjm1                                                                     
  error=.true. 
  tjm=0.d0 
  IF(lench(pdate).NE.5) RETURN                                                                  
! Year                                                                  
  IF(pdate(1:1).EQ.'I') THEN 
     year=1800 
  ELSEIF(pdate(1:1).EQ.'J') THEN 
     year=1900 
  ELSEIF(pdate(1:1).EQ.'K') THEN 
     year=2000 
  ELSE 
     RETURN 
  END IF
  READ(pdate(2:3),100,ERR=10) yy 
100 FORMAT(I2) 
  year=year+yy                                                                       
! Month                                                                 
  IF(isnum(pdate(4:4))) THEN 
     READ(pdate(4:4),101,ERR=10) month 
  ELSEIF(pdate(4:4).EQ.'A') THEN 
     month=10 
  ELSEIF(pdate(4:4).EQ.'B') THEN 
     month=11 
  ELSEIF(pdate(4:4).EQ.'C') THEN 
     month=12 
  ELSE 
     RETURN 
  END IF
101 FORMAT(I1)                                                                     
! Day                                                                   
  IF(isnum(pdate(5:5))) THEN 
     READ(pdate(5:5),101,ERR=10) day 
  ELSE 
     day=ichar(pdate(5:5))-55 
     IF(day.LT.10 .OR. day.GT.31) GOTO 10 
  END IF
  tjm=tjm1(day,month,year,0.d0) 
  error=.false.                                                                       
10 CONTINUE                                                             
END SUBROUTINE mpcdat
! Copyright (C) 1998 by Mario Carpino (carpino@brera.mi.astro.it),      
!                                                                       
! Version: June 19, 1998                                                
! --------------------------------------------------------------------- 
!                                                                       
!  *****************************************************************    
!  *                                                               *    
!  *                         F I X C N M                           *    
!  *                                                               *    
!  *                Fix covariance/normal matrices                 *    
!  *                                                               *    
!  *****************************************************************    
!                                                                       
! IN/OUT:   DEFCOV    -  Tells whether the covariance matrix is defined 
!           DEFNOR    -  Tells whether the normal matrix is defined     
!           COVE      -  Covariance matrix of orbital elements          
!           NORE      -  Normal matrix of orbital elements              
!                                                                       
! OUTPUT:   DEFCN     -  Tells whether covariance/normal matrices       
!                            are defined                                
!                                                                       
! The purpose of this routine is to compute the covariance or           
! normal matrix of orbital elements when only one of the two is         
! available                                                             
!                                                                       
SUBROUTINE fixcnm(defcov,defnor,defcn,cove,nore) 
   IMPLICIT NONE 
   DOUBLE PRECISION cove(6,6),nore(6,6) 
   LOGICAL defcov,defnor,defcn 
   DOUBLE PRECISION err
   DOUBLE PRECISION tmp(6) 
   INTEGER i,k,indp 
   err=epsilon(1.d0)*100 
   defcn=(defcov.AND.defnor) 
   IF(defcn) RETURN 
                                                                        
   IF(defcov) THEN 
      nore=cove
      CALL tchol(nore,6,6,indp,err) 
      defnor=(indp.EQ.0) 
      IF(defnor) CALL inver(nore,tmp,6,6) 
   END IF
   IF(defnor) THEN 
      cove=nore
      CALL tchol(cove,6,6,indp,err) 
      defcov=(indp.EQ.0) 
      IF(defcov) CALL inver(cove,tmp,6,6) 
   END IF
   defcn=(defcov.AND.defnor) 
   IF(defcn) RETURN
! failure case 
   cove=0.d0
   nore=0.d0
   defcov=.false. 
   defnor=.false. 
END SUBROUTINE fixcnm



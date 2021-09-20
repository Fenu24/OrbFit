! Copyright (C) 1997 by Mario Carpino (carpino@brera.mi.astro.it)       
! Version: October 16, 1997                                             
! --------------------------------------------------------------------- 
!                                                                       
!  *****************************************************************    
!  *                                                               *    
!  *                         R D O P T B                           *    
!  *                                                               *    
!  *                   Read options for BINEPH                     *    
!  *                                                               *    
!  *****************************************************************    
!                                                                       
! INPUT:    NIFLX     -  Dimension of vector INCFIL                     
!                                                                       
! OUTPUT:   NBOUT     -  Number of objects in binary ephemeris file     
!           T0        -  Epoch of initial conditions (MJD, TDT)         
!           X0,V0     -  Initial conditions (pos/vel vectors)           
!           T1        -  Starting time for ephemeris (MJD, TDT)         
!           T2        -  Ending time for ephemeris (MJD, TDT)           
!           DT        -  Ephemeris stepsize (d)                         
!                                                                       
      SUBROUTINE rdoptb(nbout,t0,x0,v0,t1,t2,dt) 
      USE reference_systems 
      USE fund_const 
      IMPLICIT NONE 
                                                                        
      INCLUDE 'dim.h90' 
                                                                        
! Max number of initial condition files                                 
      INTEGER niflx 
      PARAMETER (niflx=10) 
                                                                        
                                                                        
      INCLUDE 'comast.h90' 
      INCLUDE 'comjpl.h90' 
                                                                        
      CHARACTER*100 incfil(niflx),tmp 
      CHARACTER name1*30,eltype*3,rsys*10,epoch*10,scale*3 
      CHARACTER*30 astn1(nbmx),plnam(nbjx) 
      DOUBLE PRECISION xv(6),cove(6,6),nore(6,6),hmag,gmag,enne 
      DOUBLE PRECISION t0,x0(3,nbmx),v0(3,nbmx),t1,t2,dt,t0t 
      DOUBLE PRECISION elem(6),rot(3,3),mass,sec,sece 
      INTEGER nifl,nbl,i,k,nbout,nbadd,kr,ln,lf,i1,i2,mjd,mjde 
      LOGICAL found,fail1,fail,defcov,defnor,end,loaded(nbmx),first 
                                                                        
      INTEGER lench 
      EXTERNAL lench 
                                                                        
      fail=.false. 
      CALL rdmcha('input_files.','incond',incfil,nifl,niflx,'niflx',    &
     &            .true.,found,fail1,fail)                              
      CALL rdmcha('bineph.','eph_obj',astnam,nbout,nbmx,'nbmx',.true.,  &
     &            found,fail1,fail)                                     
      CALL rdmcha('bineph.','add_obj',astn1,nbadd,nbmx,'nbmx',.false.,  &
     &            found,fail1,fail)                                     
      IF(.NOT.found) nbadd=0 
                                                                        
      IF(fail) STOP '**** rdoptb: abnormal end ****' 
                                                                        
      nbm=nbout+nbadd 
      CALL chkpdf(nbm,nbmx,'nbmx') 
      nba=0 
      nbt=nbm 
      DO 1 i=1,nbadd 
      astnam(nbout+i)=astn1(i) 
    1 END DO 
      DO 2 i=1,nbm 
      loaded(i)=.false. 
    2 END DO 
! Check for objects included more than once                             
      first=.true. 
      DO 10 i1=1,nbm-1 
      DO 11 i2=i1+1,nbm 
      IF(astnam(i1).EQ.astnam(i2)) THEN 
          IF(first) THEN 
              WRITE(*,104) 
              first=.false. 
          END IF 
          ln=lench(astnam(i2)) 
          WRITE(*,103) astnam(i2)(1:ln) 
          fail=.true. 
      END IF 
  104 FORMAT(' ERROR: the following objects are included more than',    &
     &       ' once'/                                                   &
     &       '        among initial conditions:')                       
  103 FORMAT(10X,A) 
   11 END DO 
   10 END DO 
      IF(fail) STOP '**** rdoptb: abnormal end ****' 
      gmsun=0.01720209895d0**2 
                                                                        
      nbl=0 
      DO 4 i=1,nifl 
      CALL oporbf(incfil(i),0) 
    3 CONTINUE 
      CALL rdorb(name1,elem,eltype,t0t,cove,defcov,nore,defnor,         &
     &     hmag,gmag,mass,rsys,epoch,kr,end)                            
      IF(end) GOTO 6 
      DO 5 k=1,nbm 
      IF(loaded(k)) GOTO 5 
      IF(name1.EQ.astnam(k)) THEN 
          astmas(k)=mass 
          gma(k)=gmsun*astmas(k) 
          gma1(k)=gma(k)+gmsun 
          CALL coocha(elem,eltype,gma1(k),xv,'CAR',enne) 
          CALL rotpn(rot,rsys,epoch,t0t,'EQUM','J2000',0.d0) 
          CALL prodmv(x0(1,k),rot,xv(1)) 
          CALL prodmv(v0(1,k),rot,xv(4)) 
          loaded(k)=.true. 
          nbl=nbl+1 
          IF(nbl.EQ.1) THEN 
              t0=t0t 
              WRITE(*,100) t0 
          ELSEIF(ABS(t0-t0t).GT.1.d-6) THEN 
              ln=lench(name1) 
              lf=lench(incfil(i)) 
              WRITE(*,101) name1(1:ln),incfil(i)(1:lf),kr,t0t 
              STOP '**** rdoptb: abnormal end ****' 
          END IF 
          IF(nbl.EQ.nbm) THEN 
              CALL clorbf 
              GOTO 7 
          END IF 
      END IF 
  100 FORMAT(' Epoch of initial conditions:',F13.5,' (MJD, TDT)') 
  101 FORMAT(' ERROR: object ',A,' (file ',A,', record',I6,')'/         &
     &       '        has orbital elements at a different epoch'/       &
     &       '        T0 =',F13.5,' (MJD, TDT)')                        
    5 END DO 
      GOTO 3 
    6 CONTINUE 
      CALL clorbf 
    4 END DO 
      WRITE(*,102) 
  102 FORMAT(' ERROR: no initial conditions found for the following',   &
     &       ' objects:')                                               
      DO 9 i=1,nbm 
      IF(.NOT.loaded(i)) THEN 
          ln=lench(astnam(i)) 
          WRITE(*,103) astnam(i)(1:ln) 
      END IF 
    9 END DO 
      STOP '**** rdoptb: abnormal end ****' 
    7 CONTINUE 
                                                                        
! Limits and stepsize of binary ephemeris                               
      CALL rdntim('bineph.epoch.','start',tmp,mjd,sec,scale,.true.,     &
     &            found,fail1,fail)                                     
      IF(.NOT.fail1) THEN 
          CALL cnvtim(mjd,sec,scale,mjde,sece,'TDT') 
          t1=mjde+sece/86400.d0 
      END IF 
      CALL rdntim('bineph.epoch.','end',tmp,mjd,sec,scale,.true.,       &
     &            found,fail1,fail)                                     
      IF(.NOT.fail1) THEN 
          CALL cnvtim(mjd,sec,scale,mjde,sece,'TDT') 
          t2=mjde+sece/86400.d0 
      END IF 
      CALL rdnrea('bineph.','step',dt,.true.,found,fail1,fail) 
                                                                        
! Planets from JPL Digital Ephemeris                                    
      CALL rdmcha('force.JPLDE.','planets',plnam,nbj,nbjx,'nbjx',.true.,&
     &            found,fail1,fail)                                     
      CALL jpllis(plnam,nbj,gmj,jplid,fail) 
                                                                        
      IF(fail) STOP '**** rdoptb: abnormal end ****' 
                                                                        
      END SUBROUTINE rdoptb                                          

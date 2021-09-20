MODULE iodet

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: io_det_ades

CONTAINS                                           


! Copyright (C) 1998-2000 by Mario Carpino (carpino@brera.mi.astro.it)  
! Version: June 9, 2000
! Updated: July 2020   AB, GFG
! --------------------------------------------------------------------- 
!                                                                       
!  *****************************************************************    
!  *                                                               *    
!  *                     I O D E T _ A D E S                       *    
!  *                                                               *    
!  *      Initial orbit determination (Gauss, Vaisala, etc.)       *    
!  *                                                               *    
!  *****************************************************************    
!                                                                       
! INPUT:    UNIREP    -  FORTRAN unit for report (0 for no report)      
!           RWOFIL    -  File for residual output (or blank)            
!           NAME      -  Object name   
!           OBS, OBSW - observations and weights         
!           OBSW%SEL_COORD   -  Selection index (0=don't use; 1=use for fit;   
!                                         2=use for fit & Gauss method) 
!           N         -  Number of observations                         
!           IFO       -  Vaisala interactive level:                     
!                           1=interactive (fitobs)                      
!                           2=noninteractive (orbfit)                   
!                                                                       
! OUTPUT:   ELEM      -  Orbital elements                               
!           TELEM     -  Epoch of orbital elements (MJD, TDT)           
!           RESA      -  Residuals in RA (rad)                          
!           RESD      -  Residuals in DEC (rad)                         
!           ELTYPE    -  Type of orbital elements (EQU/KEP/CAR)         
!           COMELE    -  Comment on orbital elements                    
!           FAIL      -  Error flag                                     
!                                                                       
! CONCEPT: the routines tries to compute a preliminary orbit using      
! different triplets of observations and different method, possibly     
! applying also small random increments to the astrometric measurements 
! in order to counteract observation errors, but returns only one       
! solution (the first elliptic solution encountered)                    
!                    

  SUBROUTINE io_det_ades(unirep,rwofil,name,obs,obsw,n,error_model,ifo,&
       &                elem,telem,eltype,comele,fail,ades_type) 
    USE fund_const
    USE astrometric_observations
    USE ades
    IMPLICIT NONE 
    ! INPUT observations: new data types
    INTEGER, INTENT(IN) :: n  !number of observations
    TYPE(ast_obs),DIMENSION(n),INTENT(IN) :: obs
    TYPE(ast_wbsr),DIMENSION(n),INTENT(INOUT) :: obsw
    CHARACTER*20, INTENT(IN) ::  error_model ! weighing model
    ! other input 
    CHARACTER*(*), INTENT(IN) :: name
    CHARACTER*(*), INTENT(IN) :: rwofil
    INTEGER,INTENT(IN) :: unirep ! unit for output
    INTEGER, INTENT(IN) :: ifo ! Vaisala interactive level
    ! OUTPUT 
    DOUBLE PRECISION, INTENT(OUT) :: elem(6),telem       
    CHARACTER*(*), INTENT(OUT) :: eltype,comele
    LOGICAL, INTENT(OUT) :: fail
    INTEGER, INTENT(IN),OPTIONAL :: ades_type   ! flag for type of ADES output file (.psv/.xml)
    INTEGER ades_type_loc
    ! END INTERFACE
    ! NEEDED parameters, common blocks
    INCLUDE 'pariod.h90' 
    INCLUDE 'comiod.h90' 

    INTEGER, PARAMETER :: ncx=7 

    DOUBLE PRECISION elemv(6,8),t0v(8),tdt3(3),alpha3(3),delta3(3),rtop(3) 
    DOUBLE PRECISION rmsb,rms1,noisea(3),noised(3),dt1,dt2,dtc 
    INTEGER ln,lm,selipt(3),itry,nsol,i,lmsg,ngr,is1,is2,is3,is
    INTEGER ir(2),irb(2),obsc3(3),nkep,in,selb(3) 
    INTEGER nc(ncx) 
    CHARACTER msg*100,methb*20,label*50 
    CHARACTER*4 ieltyv(8),ieltyb 
    LOGICAL eot,existb
    INTEGER, EXTERNAL :: lench 
    DOUBLE PRECISION, EXTERNAL :: rvnorm 
    DOUBLE PRECISION rmsmag ! aggiunto io

    IF(.NOT.PRESENT(ades_type)) THEN
       ades_type_loc = 0   ! if flag is not present, we set it to print only .rwo file
    ELSE
       ades_type_loc = ades_type
    ENDIF

    IF(iiciod.NE.36) STOP '**** iodet: internal error (01) ****' 

    fail=.true. 
    existb=.false. 
    rmsb=0 
    methb=' '

    ln=lench(name) 
    WRITE(*,202) name(1:ln) 
    IF(unirep.GT.0) WRITE(unirep,202) name(1:ln) 
202 FORMAT('Preliminary orbit for object ',A,':') 

    ! Computation of best triplets of observations                          
    CALL sel3_mc(obs%time_tdt,obs%type,    &
         &     obsw%rms_coord(1),obsw%rms_coord(2),obsw%sel_coord,n,iodntr) 

    ! Statistics:                                                           
    ! NC(1) = total number of triplets tried                                
    ! NC(2) = total number of noise sample added (including no noise)       
    ! NC(3) = total number of roots (Gauss' method)                         
    ! NC(4) = total number of acceptable solutions (Gauss' method)          
    ! NC(5) = total number of elliptic solutions (Gauss' method)            
    ! NC(6) = total number of trials (Vaisala method)                       
    ! NC(7) = total number of acceptable solutions (Vaisala method)         
    DO i=1,ncx 
       nc(i)=0 
    ENDDO

    ! LOOP 1 (on triplets of observations)                                  
10  CONTINUE 
    CALL sel3mg(selipt,eot) 
    IF(eot) GOTO 20 
    ! Total number of triplets tried                                        
    nc(1)=nc(1)+1 
    ! Selected points                                                       
    DO  i=1,3 
       is=selipt(i) 
       tdt3(i)=obs(is)%time_tdt 
       obsc3(i)=obs(is)%obscod_i 
    ENDDO
    is1=selipt(1) 
    is2=selipt(2) 
    is3=selipt(3)
    IF(tdt3(3).eq.tdt3(2).or.tdt3(2).eq.tdt3(1))THEN
       WRITE(*,203) 'Failure',selipt,tdt3(2)-tdt3(1),tdt3(3)-tdt3(2)
       IF(unirep.GT.0) WRITE(unirep,203) 'Failure',selipt, &
            &        tdt3(2)-tdt3(1),tdt3(3)-tdt3(2)
       fail=.TRUE.
       RETURN
    ENDIF
    IF(iodvrb.GE.2) THEN 
       WRITE(*,203) 'Trying',selipt,                                 &
            &                 tdt3(2)-tdt3(1),tdt3(3)-tdt3(2)              
       IF(unirep.GT.0) WRITE(unirep,203) 'Trying',selipt,            &
            &                    tdt3(2)-tdt3(1),tdt3(3)-tdt3(2)           
    END IF
203 FORMAT(4X,A,' observations:',3I6,' (DT=',2F10.2,' d)') 
    CALL iodsdt(selipt,obs%time_tdt,n,iodexp,ioddtm,ir) 
    IF(iodvrb.GE.2) THEN 
       WRITE(*,204) ir,obs(ir(2))%time_tdt-obs(ir(1))%time_tdt 
       IF(unirep.GT.0) WRITE(unirep,204) ir,                         &
            &    obs(ir(2))%time_tdt-obs(ir(1))%time_tdt
    END IF
204 FORMAT(4X,'RMS check observations:', I6,' -',I6,                  &
         &          ' (DT=',F10.2,' d)')                                    

    ! LOOP 2 (on different realizations of noise)                           
    DO 11 in=0,iodnit 
       IF(in.EQ.0) THEN 
          DO 14 i=1,3 
             is=selipt(i) 
             noisea(i)=0 
             noised(i)=0 
             alpha3(i)=obs(is)%coord(1) 
             delta3(i)=obs(is)%coord(2) 
14        ENDDO
       ELSE 
          DO 15 i=1,3 
             is=selipt(i) 
             noisea(i)=obsw(is)%rms_coord(1)*iodksi*rvnorm() 
             noised(i)=obsw(is)%rms_coord(2)*iodksi*rvnorm() 
             alpha3(i)=obs(is)%coord(1)+noisea(i) 
             delta3(i)=obs(is)%coord(2)+noised(i) 
15        ENDDO
       END IF
       nc(2)=nc(2)+1 

       ! LOOP 3 (on different methods)                                         
       DO 1 itry=1,iodnm 
          lm=lench(iodmen(itry)) 
          msg=' ' 
          IF((iodvrb.GE.2) .AND. iodmul) THEN 
             WRITE(*,210) iodmen(itry)(1:lm) 
             IF(unirep.GT.0) WRITE(unirep,210) iodmen(itry)(1:lm) 
          END IF
210       FORMAT(8X,'Trying ',A,' method') 
          IF(iodvrb.GE.3) THEN 
             WRITE(*,240) iodmen(itry)(1:lm),name(1:ln),in,                &
                  &                 (noisea(i),noised(i),i=1,3)        
             IF(unirep.GT.0) WRITE(unirep,240) iodmen(itry)(1:lm),         &
                  &                    name(1:ln),in,                                &
                  &                    (noisea(i),noised(i),i=1,3)     
          END IF
          IF(iodmet(itry).EQ.1) THEN 
             CALL gaussn(tdt3,alpha3,delta3,obsc3,elemv,ieltyv,t0v,        &
                  &                ngr,nsol,rtop,fail,msg,(iodvrb.GE.3),iodmul)           
             nc(3)=nc(3)+ngr 
             nc(4)=nc(4)+nsol 
             nkep=0 
             DO 2 is=1,nsol 
                IF(ieltyv(is).EQ.'KEP') nkep=nkep+1 
2            ENDDO
             nc(5)=nc(5)+nkep 

             IF(iodvrb.GE.2) THEN 
                lmsg=lench(msg) 
                IF(lmsg.GT.0) THEN 
                   WRITE(*,220) 'Gauss',name(1:ln),ngr,nsol,nkep,        &
                        &                         msg(1:lmsg)                              
                   IF(unirep.GT.0) WRITE(unirep,220) 'Gauss',name(1:ln), &
                        &                            ngr,nsol,nkep,msg(1:lmsg)             
                ELSE 
                   WRITE(*,221) 'Gauss',name(1:ln),ngr,nsol,nkep 
                   IF(unirep.GT.0) WRITE(unirep,221)'Gauss', name(1:ln), &
                        &                            ngr,nsol,nkep                         
                END IF
             END IF
             DO 3 is=1,nsol
                IF(elemv(2,is).gt.1.d0.and.ieltyv(is).eq.'KEP')THEN
                   WRITE(*,*)' io_det: inconsistent elements type ', &
                        &                   elemv(1:6,is),ieltyv(is)
                   CYCLE
                ELSEIF(elemv(2,is).gt.0.99)THEN
                   WRITE(*,*)' io_det: bizarre orbit ',elemv(1:6,is)
                   CYCLE
                ENDIF
                CALL iod_rms(elemv(1,is),ieltyv(is),t0v(is),obs,obsw,         &
                     &              ir(1),ir(2),rms1)       
                IF(iodvrb.GE.2) THEN 
                   IF(ieltyv(is).EQ.'KEP') THEN 
                      WRITE(*,222) 'Gauss',name(1:ln),is,'a',elemv(1,is),   &
                           &                          elemv(2,is),rms1                 
                      IF(unirep.GT.0) WRITE(unirep,222)'Gauss', name(1:ln), &
                           &                    is,'a',elemv(1,is),elemv(2,is),rms1   
                   ELSE 
                      WRITE(*,222) 'Gauss',name(1:ln),is,'q',elemv(1,is),   &
                           &                          elemv(2,is),rms1                 
                      IF(unirep.GT.0) WRITE(unirep,222) 'Gauss',name(1:ln), &
                           &                     is,'q',elemv(1,is),elemv(2,is),rms1  
                   END IF
                END IF
                IF((ieltyv(is).EQ.'KEP') .AND. (rms1.LE.iodrmx))THEN
                   CALL iodsbs(elem,ieltyb,telem,rmsb,methb,selb,irb, &
                        &                    elemv(1,is),ieltyv(is),t0v(is),rms1,'Gauss',  &
                        &                    selipt,ir,existb)
                ENDIF
3            ENDDO
             IF(existb .AND. (rmsb.LE.iodrok)) GOTO 20 

          ELSEIF(iodmet(itry).EQ.2) THEN 
             nc(6)=nc(6)+1 
             CALL vaisala(tdt3,alpha3,delta3,obsc3,ifo,elemv(1,1),t0v(1),fail)
             nsol=1 
             IF(fail) nsol=0 
             nc(7)=nc(7)+nsol 
             IF(iodvrb.GE.2) THEN 
                WRITE(*,221) 'Vaisala',name(1:ln),nsol,nsol,nsol 
                IF(unirep.GT.0) WRITE(unirep,221) 'Vaisala',name(1:ln),   &
                     &                                           nsol,nsol,nsol         
             END IF
             ieltyv(1)='KEP' 
             DO 5 is=1,nsol 
                IF(elemv(2,is).gt.1.d0.and.ieltyv(is).eq.'KEP')THEN
                   WRITE(*,*)' iodet_vaisala: inconsistent elements type ',&
                        &            elemv(1:6,is),ieltyv(is)
                   CYCLE
                ELSEIF(elemv(2,is).gt.0.99)THEN
                   WRITE(*,*)' iodet_vaisala: bizarre orbit ',elemv(1:6,is)
                   CYCLE
                ENDIF
                CALL iod_rms(elemv(1,is),ieltyv(is),t0v(is),obs,obsw,   &
                     &                ir(1),ir(2),rms1)       
                IF(iodvrb.GE.2) THEN 
                   IF(ieltyv(is).EQ.'KEP') THEN 
                      WRITE(*,222) 'Vaisala',name(1:ln),is,'a',elemv(1,is), &
                           &                         elemv(2,is),rms1                  
                      IF(unirep.GT.0) WRITE(unirep,222) 'Vaisala',          &
                           &                            name(1:ln),is,'a',elemv(1,is),        &
                           &                            elemv(2,is),rms1               
                   ELSE 
                      WRITE(*,222) 'Vaisala',name(1:ln),is,'q',elemv(1,is), &
                           &                          elemv(2,is),rms1                 
                      IF(unirep.GT.0) WRITE(unirep,222) 'Vaisala',          &
                           &                            name(1:ln),is,'q',elemv(1,is),        &
                           &                            elemv(2,is),rms1               
                   END IF
                END IF
                IF((ieltyv(is).EQ.'KEP') .AND. (rms1.LE.iodrmx))     &
                     &        CALL iodsbs(elem,ieltyb,telem,rmsb,methb,selb,irb,        &
                     &                    elemv(1,is),ieltyv(is),t0v(is),rms1,'Vaisala',&
                     &                    selipt,ir,existb)                             
5            ENDDO
             IF(existb .AND. (rmsb.LE.iodrok)) GOTO 20 
6            CONTINUE 
          END IF
220       FORMAT(8X,A,'(',A,'): NR=',I1,' NC=',I1,' NK=',I1,' (',A,')') 
221       FORMAT(8X,A,'(',A,'): NR=',I1,' NC=',I1,' NK=',I1) 
222       FORMAT(12X,A,'(',A,'#',I1,'): ',A,'=',F9.4,' e=',F9.4,' RMS=',    &
               &       1P,E10.2,' arcsec')                                        
240       FORMAT(8X,A,'(',A,'): Iter =',I4,'; Noise =',6F9.3) 

          ! END OF LOOP 3                                                         
1      END DO

       ! END OF LOOP 1                                                         
11  END DO

    ! END OF LOOP 1                                                         
    GOTO 10 
20  CONTINUE 
    IF(existb) THEN 
       dt1=obs(selb(2))%time_tdt-obs(selb(1))%time_tdt 
       dt2=obs(selb(3))%time_tdt-obs(selb(2))%time_tdt 
       dtc=obs(irb(2))%time_tdt-obs(irb(1))%time_tdt 
       WRITE(*,230) name(1:ln),nc,selb,irb,dt1,dt2,dtc 
       IF(unirep.GT.0) WRITE(unirep,230)                             &
            &                        name(1:ln),nc,selb,ir,dt1,dt2,dtc         
    ELSE 
       WRITE(*,235) name(1:ln),nc 
       IF(unirep.GT.0) WRITE(unirep,235) name(1:ln),nc 
    END IF
230 FORMAT(4X,'IOD(',A,') STAT:',7I6,' SEL:',3I6,' IR:',2I6,' DT:',   &
         &       3F10.2)                                                    
235 FORMAT(4X,'IOD(',A,') STAT:',7I6) 

    fail=(.NOT.existb) 
    IF(existb) THEN 
       lm=lench(methb) 
       IF(in.LE.0) THEN 
          WRITE(*,232) name(1:ln),methb(1:lm),rmsb
          IF(unirep.GT.0) WRITE(unirep,232) name(1:ln),methb(1:lm),rmsb
          comele=methb(1:lm) 
       ELSE 
          WRITE(*,233) name(1:ln),methb(1:lm),rmsb 
          IF(unirep.GT.0) WRITE(unirep,233) name(1:ln),methb(1:lm),rmsb
          comele=methb(1:lm)//'(wN)' 
       END IF
       label=name(1:ln)//'/'//methb(1:lm) 
       CALL outele(unirep,elem,ieltyb,telem,label,iodmul,.true.) 
       IF(ieltyb.NE.'KEP')                                           &
            &        STOP '**** iodet: internal error (02) ****'               
       eltype='KEP' 
       IF(elem(2).gt.1.d0.and.ieltyb.eq.'KEP')THEN
          WRITE(*,*)' io_det: inconsistent elements type ',elem(1:6),ieltyb
          STOP '**** iodet: error new *********'
       ENDIF
       CALL iod_rms(elem,ieltyb,telem,obs,obsw,                &
            &                1,n,rms1)               
       obsw(selb(1))%sel_coord=2 
       obsw(selb(2))%sel_coord=2 
       obsw(selb(3))%sel_coord=2 
!!!!!!??????????
       rmsmag=0.d0
!!!!!!??????????
       IF(rwofil.NE.' ' .AND. ades_type_loc .GE. 0) THEN
          CALL write_rwo(rwofil,obs,obsw,n,error_model,rmsb,rmsmag)
       ENDIF
       IF (ades_type_loc .NE. 0) THEN
          CALL convert_fitobs_ades(obs,obsw,.TRUE.,.TRUE.)
          CALL write_ades(psv_name(1:len_ades_name),error_model,ades_type_loc,rmsb,rmsmag)
       END IF
    ELSE 
       WRITE(*,231) name(1:ln) 
       IF(unirep.GT.0) WRITE(unirep,231) name(1:ln) 
       comele='FAIL' 
    END IF
231 FORMAT(4X,'IOD(',A,'): FAILED') 
232 FORMAT(4X,'IOD(',A,'): solution with ',A,' method (no noise): ',  &
         &           'RMS=',1P,E10.2,' arcsec')                             
233 FORMAT(4X,'IOD(',A,'): solution with ',A,' method (with noise): ',&
         &           'RMS=',1P,E10.2,' arcsec')                             

  END SUBROUTINE io_det_ades

END MODULE iodet

  
! Copyright (C) 1998-2000 by Mario Carpino (carpino@brera.mi.astro.it)  
! Version: June 9, 2000                                                 
! --------------------------------------------------------------------- 
!                                                                       
!  *****************************************************************    
!  *                                                               *    
!  *                          I O D E T                            *    
!  *                                                               *    
!  *      Initial orbit determination (Gauss, Vaisala, etc.)       *    
!  *                                                               *    
!  *****************************************************************    
!                                                                       
! INPUT:    UNIREP    -  FORTRAN unit for report (0 for no report)      
!           RWOFIL    -  File for residual output (or blank)            
!           NAME      -  Object name   
!           OBS, OBSW - observations and weights         
!           OBSW%SEL_COORD   -  Selection index (0=don't use; 1=use for fit;   
!                                         2=use for fit & Gauss method) 
!           N         -  Number of observations                         
!           IFO       -  Vaisala interactive level:                     
!                           1=interactive (fitobs)                      
!                           2=noninteractive (orbfit)                   
!                                                                       
! OUTPUT:   ELEM      -  Orbital elements                               
!           TELEM     -  Epoch of orbital elements (MJD, TDT)           
!           RESA      -  Residuals in RA (rad)                          
!           RESD      -  Residuals in DEC (rad)                         
!           ELTYPE    -  Type of orbital elements (EQU/KEP/CAR)         
!           COMELE    -  Comment on orbital elements                    
!           FAIL      -  Error flag                                     
!                                                                       
! CONCEPT: the routines tries to compute a preliminary orbit using      
! different triplets of observations and different method, possibly     
! applying also small random increments to the astrometric measurements 
! in order to counteract observation errors, but returns only one       
! solution (the first elliptic solution encountered)                    
!                                                                       
SUBROUTINE io_det(unirep,rwofil,name,obs,obsw,n,error_model,ifo,&
     &                elem,telem,eltype,comele,fail) 
  USE fund_const
  USE astrometric_observations
  IMPLICIT NONE 
 ! INPUT observations: new data types
  INTEGER, INTENT(IN) :: n  !number of observations
  TYPE(ast_obs),DIMENSION(n),INTENT(IN) :: obs
  TYPE(ast_wbsr),DIMENSION(n),INTENT(INOUT) :: obsw
  CHARACTER*20, INTENT(IN) ::  error_model ! weighing model
! other input 
  CHARACTER*(*), INTENT(IN) :: name
  CHARACTER*(*), INTENT(IN) :: rwofil
  INTEGER,INTENT(IN) :: unirep ! unit for output
  INTEGER, INTENT(IN) :: ifo ! Vaisala interactive level
! OUTPUT 
  DOUBLE PRECISION, INTENT(OUT) :: elem(6),telem       
  CHARACTER*(*), INTENT(OUT) :: eltype,comele
  LOGICAL, INTENT(OUT) :: fail 
! END INTERFACE
! NEEDED parameters, common blocks
  INCLUDE 'pariod.h90' 
  INCLUDE 'comiod.h90' 
                                                                        
  INTEGER, PARAMETER :: ncx=7 
                                                                 
  DOUBLE PRECISION elemv(6,8),t0v(8),tdt3(3),alpha3(3),delta3(3),rtop(3) 
  DOUBLE PRECISION rmsb,rms1,noisea(3),noised(3),dt1,dt2,dtc 
  INTEGER ln,lm,selipt(3),itry,nsol,i,lmsg,ngr,is1,is2,is3,is 
  INTEGER ir(2),irb(2),obsc3(3),nkep,in,selb(3) 
  INTEGER nc(ncx) 
  CHARACTER msg*100,methb*20,label*50 
  CHARACTER*4 ieltyv(8),ieltyb 
  LOGICAL eot,existb 
                                                                       
  INTEGER, EXTERNAL :: lench 
  DOUBLE PRECISION, EXTERNAL :: rvnorm 
  DOUBLE PRECISION rmsmag ! aggiunto io

  IF(iiciod.NE.36) STOP '**** iodet: internal error (01) ****' 
  
  fail=.true. 
  existb=.false. 
  rmsb=0 
  methb=' ' 
                                                                        
  ln=lench(name) 
  WRITE(*,202) name(1:ln) 
  IF(unirep.GT.0) WRITE(unirep,202) name(1:ln) 
202 FORMAT('Preliminary orbit for object ',A,':') 
                                                                        
! Computation of best triplets of observations                          
  CALL sel3_mc(obs%time_tdt,obs%type,    &
     &     obsw%rms_coord(1),obsw%rms_coord(2),obsw%sel_coord,n,iodntr) 
                                                                        
! Statistics:                                                           
! NC(1) = total number of triplets tried                                
! NC(2) = total number of noise sample added (including no noise)       
! NC(3) = total number of roots (Gauss' method)                         
! NC(4) = total number of acceptable solutions (Gauss' method)          
! NC(5) = total number of elliptic solutions (Gauss' method)            
! NC(6) = total number of trials (Vaisala method)                       
! NC(7) = total number of acceptable solutions (Vaisala method)         
  DO i=1,ncx 
     nc(i)=0 
  ENDDO
                                                                        
! LOOP 1 (on triplets of observations)                                  
10 CONTINUE 
  CALL sel3mg(selipt,eot) 
  IF(eot) GOTO 20 
! Total number of triplets tried                                        
  nc(1)=nc(1)+1 
! Selected points                                                       
  DO  i=1,3 
     is=selipt(i) 
     tdt3(i)=obs(is)%time_tdt 
     obsc3(i)=obs(is)%obscod_i 
  ENDDO
  is1=selipt(1) 
  is2=selipt(2) 
  is3=selipt(3)
  IF(tdt3(3).eq.tdt3(2).or.tdt3(2).eq.tdt3(1))THEN
     WRITE(*,203) 'Failure',selipt,tdt3(2)-tdt3(1),tdt3(3)-tdt3(2)
     IF(unirep.GT.0) WRITE(unirep,203) 'Failure',selipt, &
     &        tdt3(2)-tdt3(1),tdt3(3)-tdt3(2)
     fail=.TRUE.
     RETURN
  ENDIF
  IF(iodvrb.GE.2) THEN 
     WRITE(*,203) 'Trying',selipt,                                 &
     &                 tdt3(2)-tdt3(1),tdt3(3)-tdt3(2)              
     IF(unirep.GT.0) WRITE(unirep,203) 'Trying',selipt,            &
     &                    tdt3(2)-tdt3(1),tdt3(3)-tdt3(2)           
  END IF
203 FORMAT(4X,A,' observations:',3I6,' (DT=',2F10.2,' d)') 
  CALL iodsdt(selipt,obs%time_tdt,n,iodexp,ioddtm,ir) 
  IF(iodvrb.GE.2) THEN 
     WRITE(*,204) ir,obs(ir(2))%time_tdt-obs(ir(1))%time_tdt 
     IF(unirep.GT.0) WRITE(unirep,204) ir,                         &
     &    obs(ir(2))%time_tdt-obs(ir(1))%time_tdt
  END IF
204 FORMAT(4X,'RMS check observations:', I6,' -',I6,                  &
     &          ' (DT=',F10.2,' d)')                                    
                                                                        
! LOOP 2 (on different realizations of noise)                           
  DO 11 in=0,iodnit 
     IF(in.EQ.0) THEN 
        DO 14 i=1,3 
           is=selipt(i) 
           noisea(i)=0 
           noised(i)=0 
           alpha3(i)=obs(is)%coord(1) 
           delta3(i)=obs(is)%coord(2) 
14      ENDDO 
     ELSE 
        DO 15 i=1,3 
           is=selipt(i) 
           noisea(i)=obsw(is)%rms_coord(1)*iodksi*rvnorm() 
           noised(i)=obsw(is)%rms_coord(2)*iodksi*rvnorm() 
           alpha3(i)=obs(is)%coord(1)+noisea(i) 
           delta3(i)=obs(is)%coord(2)+noised(i) 
   15   ENDDO 
     END IF
     nc(2)=nc(2)+1 
                                                                        
! LOOP 3 (on different methods)                                         
     DO 1 itry=1,iodnm 
        lm=lench(iodmen(itry)) 
        msg=' ' 
        IF((iodvrb.GE.2) .AND. iodmul) THEN 
           WRITE(*,210) iodmen(itry)(1:lm) 
           IF(unirep.GT.0) WRITE(unirep,210) iodmen(itry)(1:lm) 
        END IF
210     FORMAT(8X,'Trying ',A,' method') 
        IF(iodvrb.GE.3) THEN 
           WRITE(*,240) iodmen(itry)(1:lm),name(1:ln),in,                &
     &                 (noisea(i),noised(i),i=1,3)        
           IF(unirep.GT.0) WRITE(unirep,240) iodmen(itry)(1:lm),         &
     &                    name(1:ln),in,                                &
     &                    (noisea(i),noised(i),i=1,3)     
        END IF
        IF(iodmet(itry).EQ.1) THEN 
           CALL gaussn(tdt3,alpha3,delta3,obsc3,elemv,ieltyv,t0v,        &
     &                ngr,nsol,rtop,fail,msg,(iodvrb.GE.3),iodmul)           
           nc(3)=nc(3)+ngr 
           nc(4)=nc(4)+nsol 
           nkep=0 
           DO 2 is=1,nsol 
              IF(ieltyv(is).EQ.'KEP') nkep=nkep+1 
2          ENDDO
           nc(5)=nc(5)+nkep 
                                                                        
           IF(iodvrb.GE.2) THEN 
              lmsg=lench(msg) 
              IF(lmsg.GT.0) THEN 
                 WRITE(*,220) 'Gauss',name(1:ln),ngr,nsol,nkep,        &
     &                         msg(1:lmsg)                              
                 IF(unirep.GT.0) WRITE(unirep,220) 'Gauss',name(1:ln), &
     &                            ngr,nsol,nkep,msg(1:lmsg)             
              ELSE 
                 WRITE(*,221) 'Gauss',name(1:ln),ngr,nsol,nkep 
                 IF(unirep.GT.0) WRITE(unirep,221)'Gauss', name(1:ln), &
     &                            ngr,nsol,nkep                         
              END IF
           END IF
           DO 3 is=1,nsol
              IF(elemv(2,is).gt.1.d0.and.ieltyv(is).eq.'KEP')THEN
                 WRITE(*,*)' io_det: inconsistent elements type ', &
     &                   elemv(1:6,is),ieltyv(is)
                 CYCLE
              ELSEIF(elemv(2,is).gt.0.99)THEN
                 WRITE(*,*)' io_det: bizarre orbit ',elemv(1:6,is)
                 CYCLE
              ENDIF
              CALL iod_rms(elemv(1,is),ieltyv(is),t0v(is),obs,obsw,         &
     &              ir(1),ir(2),rms1)       
              IF(iodvrb.GE.2) THEN 
                 IF(ieltyv(is).EQ.'KEP') THEN 
                    WRITE(*,222) 'Gauss',name(1:ln),is,'a',elemv(1,is),   &
     &                          elemv(2,is),rms1                 
                    IF(unirep.GT.0) WRITE(unirep,222)'Gauss', name(1:ln), &
     &                    is,'a',elemv(1,is),elemv(2,is),rms1   
                 ELSE 
                    WRITE(*,222) 'Gauss',name(1:ln),is,'q',elemv(1,is),   &
     &                          elemv(2,is),rms1                 
                    IF(unirep.GT.0) WRITE(unirep,222) 'Gauss',name(1:ln), &
     &                     is,'q',elemv(1,is),elemv(2,is),rms1  
                 END IF
              END IF
              IF((ieltyv(is).EQ.'KEP') .AND. (rms1.LE.iodrmx))THEN
                 CALL iodsbs(elem,ieltyb,telem,rmsb,methb,selb,irb, &
     &                    elemv(1,is),ieltyv(is),t0v(is),rms1,'Gauss',  &
     &                    selipt,ir,existb)
              ENDIF                             
3          ENDDO
           IF(existb .AND. (rmsb.LE.iodrok)) GOTO 20 
                                                                        
        ELSEIF(iodmet(itry).EQ.2) THEN 
           nc(6)=nc(6)+1 
           CALL vaisala(tdt3,alpha3,delta3,obsc3,ifo,elemv(1,1),t0v(1),fail)
           nsol=1 
           IF(fail) nsol=0 
           nc(7)=nc(7)+nsol 
           IF(iodvrb.GE.2) THEN 
              WRITE(*,221) 'Vaisala',name(1:ln),nsol,nsol,nsol 
              IF(unirep.GT.0) WRITE(unirep,221) 'Vaisala',name(1:ln),   &
     &                                           nsol,nsol,nsol         
           END IF
           ieltyv(1)='KEP' 
           DO 5 is=1,nsol 
              IF(elemv(2,is).gt.1.d0.and.ieltyv(is).eq.'KEP')THEN
                 WRITE(*,*)' iodet_vaisala: inconsistent elements type ',&
     &            elemv(1:6,is),ieltyv(is)
                 CYCLE
              ELSEIF(elemv(2,is).gt.0.99)THEN
                 WRITE(*,*)' iodet_vaisala: bizarre orbit ',elemv(1:6,is)
                 CYCLE
              ENDIF
              CALL iod_rms(elemv(1,is),ieltyv(is),t0v(is),obs,obsw,   &
     &                ir(1),ir(2),rms1)       
              IF(iodvrb.GE.2) THEN 
                 IF(ieltyv(is).EQ.'KEP') THEN 
                    WRITE(*,222) 'Vaisala',name(1:ln),is,'a',elemv(1,is), &
     &                         elemv(2,is),rms1                  
                    IF(unirep.GT.0) WRITE(unirep,222) 'Vaisala',          &
     &                            name(1:ln),is,'a',elemv(1,is),        &
     &                            elemv(2,is),rms1               
                 ELSE 
                    WRITE(*,222) 'Vaisala',name(1:ln),is,'q',elemv(1,is), &
     &                          elemv(2,is),rms1                 
                    IF(unirep.GT.0) WRITE(unirep,222) 'Vaisala',          &
     &                            name(1:ln),is,'q',elemv(1,is),        &
     &                            elemv(2,is),rms1               
                 END IF
              END IF
              IF((ieltyv(is).EQ.'KEP') .AND. (rms1.LE.iodrmx))     &
     &        CALL iodsbs(elem,ieltyb,telem,rmsb,methb,selb,irb,        &
     &                    elemv(1,is),ieltyv(is),t0v(is),rms1,'Vaisala',&
     &                    selipt,ir,existb)                             
5          ENDDO
           IF(existb .AND. (rmsb.LE.iodrok)) GOTO 20 
6          CONTINUE 
        END IF
  220 FORMAT(8X,A,'(',A,'): NR=',I1,' NC=',I1,' NK=',I1,' (',A,')') 
  221 FORMAT(8X,A,'(',A,'): NR=',I1,' NC=',I1,' NK=',I1) 
  222 FORMAT(12X,A,'(',A,'#',I1,'): ',A,'=',F9.4,' e=',F9.4,' RMS=',    &
     &       1P,E10.2,' arcsec')                                        
  240 FORMAT(8X,A,'(',A,'): Iter =',I4,'; Noise =',6F9.3) 
                                                                        
! END OF LOOP 3                                                         
1    END DO
                                                                        
! END OF LOOP 1                                                         
11 END DO
                                                                        
! END OF LOOP 1                                                         
  GOTO 10 
20 CONTINUE 
  IF(existb) THEN 
     dt1=obs(selb(2))%time_tdt-obs(selb(1))%time_tdt 
     dt2=obs(selb(3))%time_tdt-obs(selb(2))%time_tdt 
     dtc=obs(irb(2))%time_tdt-obs(irb(1))%time_tdt 
     WRITE(*,230) name(1:ln),nc,selb,irb,dt1,dt2,dtc 
     IF(unirep.GT.0) WRITE(unirep,230)                             &
     &                        name(1:ln),nc,selb,ir,dt1,dt2,dtc         
  ELSE 
     WRITE(*,235) name(1:ln),nc 
     IF(unirep.GT.0) WRITE(unirep,235) name(1:ln),nc 
  END IF
230 FORMAT(4X,'IOD(',A,') STAT:',7I6,' SEL:',3I6,' IR:',2I6,' DT:',   &
     &       3F10.2)                                                    
235 FORMAT(4X,'IOD(',A,') STAT:',7I6) 
                                                                        
  fail=(.NOT.existb) 
  IF(existb) THEN 
     lm=lench(methb) 
     IF(in.LE.0) THEN 
        WRITE(*,232) name(1:ln),methb(1:lm),rmsb
        IF(unirep.GT.0) WRITE(unirep,232) name(1:ln),methb(1:lm),rmsb
        comele=methb(1:lm) 
     ELSE 
        WRITE(*,233) name(1:ln),methb(1:lm),rmsb 
        IF(unirep.GT.0) WRITE(unirep,233) name(1:ln),methb(1:lm),rmsb
        comele=methb(1:lm)//'(wN)' 
     END IF
     label=name(1:ln)//'/'//methb(1:lm) 
     CALL outele(unirep,elem,ieltyb,telem,label,iodmul,.true.) 
     IF(ieltyb.NE.'KEP')                                           &
     &        STOP '**** iodet: internal error (02) ****'               
     eltype='KEP' 
     IF(elem(2).gt.1.d0.and.ieltyb.eq.'KEP')THEN
        WRITE(*,*)' io_det: inconsistent elements type ',elem(1:6),ieltyb
        STOP '**** iodet: error new *********'
     ENDIF
     CALL iod_rms(elem,ieltyb,telem,obs,obsw,                &
     &                1,n,rms1)               
     obsw(selb(1))%sel_coord=2 
     obsw(selb(2))%sel_coord=2 
     obsw(selb(3))%sel_coord=2 
!!!!!!??????????
     rmsmag=0.d0
!!!!!!??????????
     IF(rwofil.NE.' ')                                             &
     &        CALL write_rwo(rwofil,obs,obsw,n,error_model,rmsb,rmsmag)  
  ELSE 
     WRITE(*,231) name(1:ln) 
     IF(unirep.GT.0) WRITE(unirep,231) name(1:ln) 
     comele='FAIL' 
  END IF
231 FORMAT(4X,'IOD(',A,'): FAILED') 
232 FORMAT(4X,'IOD(',A,'): solution with ',A,' method (no noise): ',  &
         &           'RMS=',1P,E10.2,' arcsec')                             
233 FORMAT(4X,'IOD(',A,'): solution with ',A,' method (with noise): ',&
     &           'RMS=',1P,E10.2,' arcsec')                             
                                                                        
END SUBROUTINE io_det



!!$ -*- f90 -*-
MODULE fir_filters
USE var_precision
IMPLICIT NONE
SAVE

  PUBLIC filter, filope ! public routines

! ===============================================
! {\bf parfil.h} ORBIT8V
! max filter length
      INTEGER, PARAMETER :: nfilma=1500
! ===============================================
! filter coefficients and size
      REAL(KIND=r_kind), DIMENSION(nfilma) :: filter_h, filter_h2 ! filter coefficients
      INTEGER :: nfil, nfil2 ! filter length
      INTEGER :: nsamp, nsamp2 ! sampling ratio
!      REAL(KIND=r_kind), PUBLIC :: hlen, hlen2 ! filter half length 
! maximum number of different data to be filtered
! this should be controlled by some module like dimension_masses
! anyway there is a control which stops if the dimensions are not enough 
      INTEGER, PARAMETER :: ndx=5501 
! maximum number of parallel filtering stages (it must be greater than nfil/nsamp)
      INTEGER, PARAMETER :: npx=20
CONTAINS
!                                                                       
!                                                                       
!  ***************************************************************      
!  *                                                             *      
!  *                       f i l t e r                           *      
!  *                                                             *      
!  *      general purpose implementation of a fir filter         *      
!  *                     with decimation                         *      
!  *                                                             *      
!  ***************************************************************      
!                                                                       
!   input:   timin  -  time of input data                               
!            datin  -  input data (vector)                              
!            ndat   -  number of input data                             
!            fstep   -  filter input step                               
!            nskip  -  number of data to be skipped the first time      
!                            the filter is called                       
!                warning: must be a variable (its value is changed)     
!                  must be >=0 the first time the routine is called     
!                                                                       
!                                                                       
!   output:  timout -  time of output data                              
!            datout -  output data (vector)                             
!            outdat -  (logical) output data available when .true.      
!                                                                       
!                                                                       
  SUBROUTINE filter(timin,datin,ndat,nskip,fstep,timout,datout,outdat)
    SAVE
    INTEGER, INTENT(in)            :: ndat
    REAL(KIND=r_kind), INTENT(in)  :: datin(ndat),timin,fstep
    REAL(KIND=r_kind), INTENT(out) :: datout(ndat),timout
    INTEGER, INTENT(inout)         :: nskip
    LOGICAL, INTENT(out)           :: outdat
! END INTERFACE
    INTEGER num_int(npx),ncont,naff
    INTEGER npar, j, k ! indexes for do loops 
    REAL(KIND=r_kind) s(ndx,npx) ! accumulator of partial sums
    REAL(KIND=r_kind) :: hlen !filter half length
    REAL(KIND=r_kind) :: hac 
    IF(ndat.gt.ndx) stop ' **** filter: ndat > ndx ****' 
! initialisation and skip                                               
    outdat=.false.
    IF(nskip.ge.0)THEN 
       s(1:ndat,1:npx)=0._r_kind  ! is this initialisation at zero
       num_int(1:npx)=0              ! really necessary?         
       ncont=nsamp-nskip 
       nskip=-1 
       npar=0 
       hlen=fstep*(nfil-1)/2
    ENDIF
! every nsamp input data, introduces a new filtering stage              
    IF(ncont.ge.nsamp)THEN 
       ncont=1 
       npar=npar+1 
       if(npar.gt.npx)stop ' **** filter: npar > npx ****' 
       s(1:ndat,npar)=0._r_kind
       num_int(npar)=0 
    ELSE 
       ncont=ncont+1 
    ENDIF
    IF(npar.le.0)RETURN 
! introduces the new data into the filters:                             
!        s(j,k) is the accumulator for the k-th filtering stage of the  
!        j-th datum                                                     
    DO  k=1,npar 
       naff=num_int(k)+1 
       num_int(k)=naff 
       hac=filter_h(naff) 
       DO  j=1,ndat 
          s(j,k)=s(j,k)+hac*datin(j)
       END DO
    END DO
! if the oldest filtering line is completed, passes the result          
!        to the following filtering stages and updates the buffer       
    if(num_int(1).ge.nfil)then 
       timout=timin-hlen       
       datout(1:ndat)=s(1:ndat,1) 
       npar=npar-1 
! rotation of filters                                                   
       DO k=1,npar 
          naff=k+1 
          num_int(k)=num_int(naff) 
          s(1:ndat,k)=s(1:ndat,naff)
       END DO
       outdat=.true. 
    else 
       outdat=.false. 
    end if
  END SUBROUTINE filter
!                                                                       
!  ***************************************************************      
!  *                                                             *      
!  *                       f i l t e r 2                         *      
!  *                                                             *      
!  *      general purpose implementation of a fir filter         *      
!  *                     with decimation                         *      
!  *                                                             *      
!  ***************************************************************      
!                                                                       
!   input:   timin  -  time of input data                               
!            datin  -  input data (vector)                              
!            ndat   -  number of input data                             
!            nskip2 -  number of data to be skipped the first time      
!                            the filter is called                       
!                      warning: must be a variable (its value is changed
!                               must be >=0 the first time the routine i
!            fstep   -  filter input step

!                                                                       
!   output:  timout -  time of output data                              
!            datout -  output data (vector)                             
!            outdat -  (logical) output data available when .true.      
!                                                                       
!                                                                       
  SUBROUTINE filter2(timin,datin,ndat,nskip2,fstep,timout,datout,outdat)
    SAVE
    INTEGER, INTENT(in)            :: ndat
    REAL(KIND=r_kind), INTENT(in)  :: datin(ndat),timin,fstep
    REAL(KIND=r_kind), INTENT(out) :: datout(ndat),timout
    INTEGER, INTENT(inout)         :: nskip2
    LOGICAL, INTENT(out)           :: outdat 
! END INTERFACE
    REAL(KIND=r_kind) s(ndx,npx)
    INTEGER  num_int(npx), ncont, naff
    INTEGER npar, j, k ! indexes for do loops 
    REAL(KIND=r_kind) :: hac 
    REAL(KIND=r_kind) :: hlen2 !filter half length
    if(ndat.gt.ndx)stop ' **** filter2: ndat > ndx ****' 
! initialisation and skip                                               
    outdat=.false. 
    if(nskip2.ge.0)then 
       s(1:ndat,1:npx)=0._r_kind
       num_int(1:npx)=0 
       ncont=nsamp2-nskip2 
       nskip2=-1 
       npar=0 
       hlen2=fstep*(nfil2-1)/2
    end if
! every nsamp input data, introduces a new filtering stage              
    if(ncont.ge.nsamp2)then 
       ncont=1 
       npar=npar+1 
       if(npar.gt.npx)stop ' **** filter2: npar > npx ****' 
       s(1:ndat,npar)=0.d0 
       num_int(npar)=0 
    else 
       ncont=ncont+1 
    end if
    if(npar.le.0)return 
! introduces the new data into the filters:                             
!        s(j,k) is the accumulator for the k-th filtering stage of the  
!        j-th datum                                                     
    DO k=1,npar 
       naff=num_int(k)+1 
       num_int(k)=naff 
       hac=filter_h2(naff) 
       DO  j=1,ndat 
          s(j,k)=s(j,k)+hac*datin(j)
       ENDDO
    ENDDO
! if the oldest filtering line is completed, passes the result          
!        to the following filtering stages and updates the buffer       
    if(num_int(1).ge.nfil2)then 
       timout=timin-hlen2 
       DO j=1,ndat 
          datout(j)=s(j,1) 
       ENDDO
       npar=npar-1 
! rotation of filters                                                   
       DO k=1,npar 
          naff=k+1 
          num_int(k)=num_int(naff) 
          DO j=1,ndat 
             s(j,k)=s(j,naff) 
          ENDDO
       ENDDO
       outdat=.true. 
    else 
       outdat=.false. 
    end if
  END SUBROUTINE filter2
!
!  ***************************************************************      
!  *                                                             *      
!  *                       f i l o p e                           *      
!  *                                                             *      
!  *                choose  filter input file                    *      
!  *                                                             *      
!  ***************************************************************
!
  SUBROUTINE filope(nsamp1,ifilt)
    INTEGER, INTENT(IN) :: nsamp1, ifilt ! sampling ratio, number of filter
    CHARACTER (LEN=3) :: sampl
    CHARACTER (LEN=72) :: filfil
    CHARACTER (LEN=72), PARAMETER :: filter_home ='./'
    LOGICAL isther
    IF (nsamp1.lt.10) then 
       WRITE (sampl, 196) nsamp1 
196    FORMAT   (i1) 
       filfil = trim(filter_home)//'filter.d'//sampl 
    ELSEIF (nsamp1.lt.100) then 
       WRITE (sampl, 197) nsamp1 
197    FORMAT   (i2) 
       filfil = trim(filter_home)//'filter.d'//sampl 
    ELSEIF (nsamp1.ge.100) then 
       WRITE (sampl, 198) nsamp1 
198    FORMAT   (i3) 
       filfil = trim(filter_home)//'filter.'//sampl 
    ENDIF
    INQUIRE (file = filfil, exist = isther) 
    IF (.not.isther) then 
       WRITE (*, * ) ' filter for sampling ', nsamp1, ' not available' 
       STOP 
    ENDIF
    OPEN (8, file = filfil, status = 'old') 
    IF(ifilt.eq.1)THEN
       CALL legfil (8)
       IF (nsamp1.ne.nsamp) THEN 
          WRITE (*,*)' filter coeff.file  ',filfil,'  contains  ',nsamp
          STOP 
       ENDIF
    ELSE
       CALL legfil2(8)
       IF (nsamp1.ne.nsamp2) THEN 
          WRITE (*,*)' filter coeff.file  ',filfil,'  contains  ',nsamp2
          STOP 
       ENDIF
    ENDIF
    CLOSE (8) 
  END SUBROUTINE filope
!  ***************************************************************      
!  *                                                             *      
!  *                  filter_skip                                *      
!  *                                                             *      
!  *             skip to keep in phase the filtered              *
!  *             and the sampled file                            *      
!  *                                                             *      
!  ***************************************************************      
!    
  SUBROUTINE filter_skip(nskip,nskipa)
!       skip needed for the two filter routines
    INTEGER, INTENT(out)  :: nskip, nskipa 
    INTEGER nskk
    nskk = mod ( (nfil + 1) / 2, nsamp) 
    IF (nskk.eq.0) then 
       nskip = 0 
    ELSE 
       nskip = nsamp - nskk 
    ENDIF
    nskipa = nskip 
  END SUBROUTINE filter_skip
!  ***************************************************************      
!  *                                                             *      
!  *                       l e g f i l                           *      
!  *                                                             *      
!  *               input of filter parameters                    *      
!  *                                                             *      
!  ***************************************************************      
!                                                                       
!                                                                       
!   input:   iun    -  input unit                                       
!                                                                       
!   output:  h      -  filter impulse response                          
!            nfil   -  filter length                                    
!            nsamp  -  sampling periodicity                             
!                                                                       
!                                                                       
  SUBROUTINE legfil(iun)
    INTEGER, INTENT(in) :: iun ! input unit, already opened  
    INTEGER :: nr ! half length (including center)
    INTEGER j
    read(iun,*)nsamp 
    read(iun,*)nfil 
    read(iun,*)nr 
    read(iun,*)(filter_h(j),j=1,nr) 
    DO  j=nr+1,nfil 
       filter_h(j)=filter_h(nfil+1-j)
    END DO
  END SUBROUTINE legfil
! second copy for second filter, if needed (for Trojans)
  SUBROUTINE legfil2(iun)
    INTEGER, INTENT(in) :: iun ! input unit, already opened  
    INTEGER :: nr2 ! half length (including center)
    INTEGER j
    read(iun,*)nsamp2 
    read(iun,*)nfil2 
    read(iun,*)nr2 
    read(iun,*)(filter_h2(j),j=1,nr2) 
    DO  j=nr2+1,nfil2 
       filter_h2(j)=filter_h2(nfil2+1-j)
    END DO
  END SUBROUTINE legfil2

!                                                                       
!  ***************************************************************      
!  *                                                             *      
!  *                       f i l t a n                           *      
!  *                                                             *      
!  *      general purpose implementation of a fir filter         *      
!  *                     with decimation                         *      
!  *                                                             *      
!  *   (for angular variables represented as integer number of   *      
!  *               revolution + angle in radians)                *      
!  *  TO BE REWRITTEN USING THE DERIVED TYPE ANGLE
!  ***************************************************************      
!                                                                       
!   input:   timin  -  time of input data                               
!            radin  -  input data vector (radians angle)                
!            nrin   -  input data vector (integer number of revolutions)
!            ndat   -  number of input data                             
!            nskip  -  number of data to be skipped the first time      
!                            the filter is called                       
!               warning: must be a variable (its value is changed)      
!                     must be >=0 the first time the routine is called
!            fstep   -  filter input step
!                                                                       
!                                                                       
!   output:  timout -  time of output data                              
!            radout -  output data vector (radians angle)               
!            nrout  -  output data vector (integer no of revolutions)   
!            outdat -  (logical) output data available when .true.      
!                                                                       
!                                                                       
  SUBROUTINE filtan(timin,radin,nrin,ndat,nskip,fstep,  &
     &                  timout,radout,nrout,outdat)  
    USE fund_const
    SAVE
    INTEGER, INTENT(in)            :: ndat
    REAL(KIND=r_kind), INTENT(in)  :: radin(ndat),timin,fstep
    REAL(KIND=r_kind), INTENT(out) :: radout(ndat),timout
    INTEGER, INTENT(in)            :: nrin(ndat)
    INTEGER, INTENT(out)           :: nrout(ndat)
    INTEGER, INTENT(inout)         :: nskip
    LOGICAL, INTENT(out)           :: outdat 
    INTEGER num_int(npx),ncont,naff,nrin1(ndx,npx)
    INTEGER npar, k ! indexes for do loops 
    REAL(KIND=r_kind) s(ndx,npx) ! accumulator of partial sums
    REAL(KIND=r_kind) :: hlen !filter half length
    REAL(KIND=r_kind) :: hac 
! execution
    IF(ndat.gt.ndx) stop ' **** filter: ndat > ndx ****' 
    outdat=.false.
    IF(nskip.ge.0)THEN 
       ncont=nsamp-nskip 
       nskip=-1 
       npar=0 
       hlen=fstep*(nfil-1)/2
    END IF
! every nsamp input data, introduces a new filtering stage              
    if(ncont.ge.nsamp)then 
       ncont=1 
       npar=npar+1 
       if(npar.gt.npx)stop ' **** filtan: npar > npx ****' 
       nrin1(1:ndat,npar)=nrin(1:ndat) 
       s(1:ndat,npar)=0.d0 
       num_int(npar)=0 
    else 
       ncont=ncont+1 
    end if
    if(npar.le.0)return 
! introduces the new data into the filters:                             
!        s(j,k) is the accumulator for the k-th filtering stage of the  
!        j-th datum                                                     
    DO  k=1,npar 
       naff=num_int(k)+1 
       num_int(k)=naff 
       hac=filter_h(naff) 
       s(1:ndat,k)=s(1:ndat,k)+hac*(radin(1:ndat)  &
&            +dpig*(nrin(1:ndat)-nrin1(1:ndat,k))) 
    END DO
! if the oldest filtering line is completed, passes the result          
!        to the following filtering stages and updates the buffer       
    if(num_int(1).ge.nfil)then 
       timout=timin-hlen 
       radout(1:ndat)=s(1:ndat,1)-floor(s(1:ndat,1)/dpig)*dpig 
       nrout(1:ndat)=floor(s(1:ndat,1)/dpig)+nrin1(1:ndat,1) 
       npar=npar-1 
! rotation of filters                                                   
       DO k=1,npar 
          naff=k+1 
          num_int(k)=num_int(naff) 
          nrin1(1:ndat,k)=nrin1(1:ndat,naff) 
          s(1:ndat,k)=s(1:ndat,naff) 
       END DO
       outdat=.true. 
    else 
       outdat=.false. 
    end if
  END SUBROUTINE filtan
!                                                                       
!  ***************************************************************      
!  *                                                             *      
!  *                       f i l t a n 2                         *      
!  *                                                             *      
!  *      general purpose implementation of a fir filter         *      
!  *                     with decimation                         *      
!  *                                                             *      
!  *   (for angular variables represented as integer number of   *      
!  *               revolution + angle in radians)                *      
!  *                                                             *      
!  ***************************************************************      
!                                                                       
!   input:   timin  -  time of input data                               
!            radin  -  input data vector (radians angle)                
!            nrin   -  input data vector (integer number of revolutions)
!            ndat   -  number of input data                             
!            nskip2 -  number of data to be skipped the first time      
!                            the filter is called                       
!                      warning: must be a variable (its value is changed
!                               must be >=0 the first time the routine  
!                               is called
!            fstep   -  filter input step       

!                                                                       
!   output:  timout -  time of output data                              
!            radout -  output data vector (radians angle)               
!            nrout  -  output data vector (integer number of revolutions
!            outdat -  (logical) output data available when .true.      
!                                                                       
!                                                                       
  SUBROUTINE filtan2(timin,radin,nrin,ndat,nskip2,fstep, &
     &                  timout,radout,nrout,outdat)                     
    USE fund_const
    SAVE
    INTEGER, INTENT(in)            :: ndat
    REAL(KIND=r_kind), INTENT(in)  :: radin(ndat),timin,fstep
    REAL(KIND=r_kind), INTENT(out) :: radout(ndat),timout
    INTEGER, INTENT(in)            :: nrin(ndat)
    INTEGER, INTENT(out)           :: nrout(ndat)
    INTEGER, INTENT(inout)         :: nskip2
    LOGICAL, INTENT(out)           :: outdat 
! END INTERFACE
    INTEGER num_int(npx),ncont,naff,nrin1(ndx,npx)
    INTEGER npar, k ! indexes for do loops 
    REAL(KIND=r_kind) s(ndx,npx) ! accumulator of partial sums
    REAL(KIND=r_kind) :: hlen2 !filter half length
    REAL(KIND=r_kind) :: hac
! execution 
    if(ndat.gt.ndx)stop ' **** filtan2: ndat > ndx ****' 
    outdat=.false. 
   if(nskip2.ge.0)then 
       s(1:ndat,1:npx)=0._r_kind
       num_int(1:npx)=0 
       ncont=nsamp2-nskip2 
       nskip2=-1 
       npar=0 
       hlen2=fstep*(nfil2-1)/2
    end if
! every nsamp input data, introduces a new filtering stage              
    if(ncont.ge.nsamp2)then 
       ncont=1 
       npar=npar+1 
       if(npar.gt.npx)stop ' **** filtan2: npar > npx ****' 
       nrin1(1:ndat,npar)=nrin(1:ndat) 
       s(1:ndat,npar)=0.d0 
       num_int(npar)=0 
    else 
       ncont=ncont+1 
    end if
    if(npar.le.0)return 
! introduces the new data into the filters:                             
!        s(j,k) is the accumulator for the k-th filtering stage of the  
!        j-th datum                                                     
    DO k=1,npar 
       naff=num_int(k)+1 
       num_int(k)=naff 
       hac=filter_h2(naff) 
       s(1:ndat,k)=s(1:ndat,k)+hac*(radin(1:ndat)  &
&            +dpig*(nrin(1:ndat)-nrin1(1:ndat,k))) 
    ENDDO
! if the oldest filtering line is completed, passes the result          
!        to the following filtering stages and updates the buffer       
    if(num_int(1).ge.nfil2)then 
       timout=timin-hlen2 
       radout(1:ndat)=s(1:ndat,1)-floor(s(1:ndat,1)/dpig)*dpig 
       nrout(1:ndat)=floor(s(1:ndat,1)/dpig)+nrin1(1:ndat,1)
       npar=npar-1 
! rotation of filters                                                   
       DO k=1,npar 
          naff=k+1 
          num_int(k)=num_int(naff) 
          nrin1(1:ndat,k)=nrin1(1:ndat,naff) 
          s(1:ndat,k)=s(1:ndat,naff) 
       ENDDO
       outdat=.true. 
    else 
       outdat=.false. 
    end if
  END SUBROUTINE filtan2

END MODULE fir_filters













! ========== FORCE_MODEL===============                          
! OUT OF MODULE
!               radar_ob  to adjust model if there is radar 
!               selpert    avoiding self perturb. of massive asteroids  
!               selpert2        "
!               vel_required  (if there is relativity)
! PUBLIC ROUTINES                                                          
!               masjpl    to adjust planetary masses etc.          
!               force     right hand side                               
!                 planast  planetary and asteroid positions             
!                 j2sun    sun oblateness perturbation                  
!                 genrel   general relativity, sun's contribution only  
!                 eihrel   general relativity, n-body  
!                                                                       
! HEADERS    
! force_model.o: \
!	../include/jplhdr.h90 \ only in masjpl, planast
!	../suit/FUND_CONST.mod \
!	../suit/PLANET_MASSES.mod \
!	../suit/REFERENCE_SYSTEMS.mod \
!	force_model.o \ out of module routines require the module
!	yark_pert.o   only in force
                                         
MODULE force_model
  USE planet_masses
  USE fund_const 
  USE dyn_param
  IMPLICIT NONE
  PRIVATE
  
  PUBLIC masjpl, force
  
! PUBLIC DATA
  LOGICAL :: velo_req ! velocity required from jpl_ephem, in force
  PUBLIC icrel,velo_req  !for vel_required (used from outside), control of istate

! control of close approach software
  INTEGER, PUBLIC :: iclap ! close approach control (0=off)

! radar quality of the orbit propagation
  LOGICAL :: radar
  INTEGER :: norad_obs
  PUBLIC radar, norad_obs

CONTAINS     
! =================================================                                
! MASJPL                                                                
!                                                                       
! version 1.4, A. Milani, Nov. 10, 1997                                 
! 1.5.2 corrected January 3, 1997                                       
!    input of planetary masses from jpl header                            
!    it assumes that dpleph has already been called                       
! =================================================                                      
  SUBROUTINE masjpl 
! for gms 
    USE fund_const
    USE planet_masses, ONLY:iatrue,iast
! planetary masses, constants, etc.                                     
    double precision gmerc,gm1,gmsun 
! ======== JPL EPHEM ============================                       
! standard JPL header                                                   
    INCLUDE 'jplhdr.h90' 
! ===============================================                       
!  loop indexes, lengths, output units                                  
    integer m,mm,ia,i,j,ln,iabe 
! functions                                                             
    integer lench 
!  string with asteroid name                                         
    character*30 string 
    integer mlun,mea,mjup,mnep 
!  strings with planet names                                            
    character*30 nomi(11) 
    data nomi/'MERCURY','VENUS','EARTH_MOON','MARS',                  &
   & 'JUPITER','SATURN','URANUS','NEPTUNE','PLUTO','MOON','EARTH'/    
    LOGICAL first
! call to jpllis2 
    DOUBLE PRECISION gmp(13)
    LOGICAL fail
    DATA first/.true./
    SAVE
    fail=.FALSE.
! Call to jpllis2: masses from JPL ephemerides
    CALL jpllis2(gmp,fail)
    IF(fail) THEN
       STOP ' **** masjpl: failed call to jpllis2 ********'
    ENDIF
    gmearth=gmp(11)
!  number of planets in our integration                                 
!  also sets the distance for close approaches                          
    npla=7 
    if(imerc.gt.0)then 
       npla=npla+1 
       m=0 
       mlun=9 
       mea=3 
       mjup=5 
       mnep=8 
    else 
       m=1 
       mlun=8 
       mea=2 
       mjup=4 
       mnep=7 
    endif
    if(iplut.gt.0)then 
       npla=npla+1 
       mlun=mlun+1 
    endif
!  distance defining a close approach: choise in model.opt              
    DO j=1,npla 
! let dist for mercury, venus and mars be determined by 'dter' (default 
! while for the earth it is determined by 'dmea' (default 0.1 AU)       
! if non-default values are desired, set the propag-options .dter and .d
! in the option file (e.g., .dter=0.05d0 and .dmea=0.2d0 for NEAs)      
       IF(imerc.gt.0)THEN 
          IF(j.eq.1.or.j.eq.2.or.j.eq.4)dmin(j)=dter 
       ELSE 
          IF(j.eq.1.or.j.eq.3)dmin(j)=dter 
       ENDIF
       IF(j.eq.mea)dmin(j)=dmea 
       IF(j.ge.mjup.and.j.le.mnep)dmin(j)=dmjup 
    ENDDO 
!                                                                       
!  gms in solar masses                                                  
    gmsun=gmp(12)
!  flags to require data                                                
    do  mm=1,12 
       listpl(mm)=0 
    enddo
!  setup of the planetary masses                                        
    DO mm=1,npla 
       m=m+1 
       gm(mm)=gmp(m) 
       ordnam(mm)=nomi(m) 
       if(m.eq.3)then 
          itarg(mm)=13 
       else 
          itarg(mm)=m 
       endif
       listpl(m)=2 
    ENDDO
! store mass of Earth-Moon plus Sun (in solar masses)                   
    gmse=gmsun+gmp(3)
! store mercury/sun mass ratio; Mercury position is needed anyway,     
! to use Sun-Mercury center of mass                                     
    if(imerc.eq.0)then 
       gmerc=gmp(1)
       listpl(1)=2 
    endif
!  the moon as possible additional massive body (then replace E+M center
!  of mass with E).                                                     
    if(iclun.gt.0)then 
       npla=npla+1 
!  close approaches to the Moon are not reported separately             
!  from close approaches to Earth                                       
       dmin(mlun)=dmoon 
!  mass of the Moon is computed from Earth/Moon mass ratio              
       emrat=gmp(13)
       gm(mlun)=gmp(10)
       ordnam(mlun)=nomi(10) 
       itarg(mlun)=10 
       listpl(10)=2 
!  if the Moon is included, the Earth is included as a single body      
!  and not as E-M center of mass                                        
       gm(mea)=gmp(11)
       ordnam(mea)=nomi(11) 
       itarg(mea)=3 
       listpl(3)=2 
    endif
!                                                                       
    if(iast.ne.0)then 
! part relative to massives asteroids                                   
       IF(first)THEN
! read names of asteroids and asteroid masses from ascii file 'CPV.abe' 
          call filopn(iabe,filbec,'old') 
          read(iabe,*) 
          read(iabe,*) 
          do  ia=1,iast 
             read(iabe,201,err=202,end=202)masbep(ia),string 
201          FORMAT(1P,E18.10,1X,A)
             ln=lench(string) 
             astnam(ia)=string(1:ln) 
          enddo
          goto 203 
202       WRITE(*,*)'masjpl: too many asteroids requested, iast=',iast 
          iast=ia-1
          iatrue=iast 
          WRITE(*,*)'masjpl: asteroids available ',iast 
203       call filclo(iabe,' ')
          first=.false.
       ENDIF 
       DO ia=1,iatrue 
!  asteroid close approach distance 
          dmin(npla+ia)=dmast
!  asteroid names                                                       
          ordnam(npla+ia)=astnam(astid(ia)) 
!  asteroid masses (unita` masse solari)                                
          gm(npla+ia)=gmsun*masbep(astid(ia))
       ENDDO
    endif
!  gmsun in the chosen units (au, day)                                  
    gm0=gmsun 
    nmass=npla+iatrue 
!  if mercury is not included, its mass is included in the mass of the s
    if(imerc.eq.0) then 
       gm1=gm0+gmerc 
       gmu=gmerc 
       gm0=gm1 
    endif
    
    return 
  END SUBROUTINE masjpl
!                                                                       
! =============================================================           
! FORCE : accelerations acting on a massless asteroid                   
! =============================================================           
! version 2.2.9; A. Milani, January 9, 2002                             
! with no recomputation of JPL ephemerides                              
! if ips=0 reconsult JPL ephem (and asteroid .bep file)                 
!          and store in location imem                                      
! if ips>0 use JPL ephem in location imem                               
!          and recompute all                                           
! if ips<0 use JPL ephem in location imem                               
!          partial recomputation (not implemented yet)               
! Accelerations on the asteroid, computed by using                      
! JPL ephemerides as source for the position of the planets             
! (variational equation second part)                                    
! close approach control                                                
! asteroid perturbations with trick to avoid self-perturbations
! ============================================================= 
! A. Del Vigna (March 2016) - Models for outgassing forces and
!                             its variational equations
! =============================================================           
  SUBROUTINE force(x,v,t0,f,ndf,idc,xxpla,ips,imem) 
    USE yark_pert
    USE non_grav
    USE eval_risk, ONLY: givenmass, massgiven
    USE astrometric_observations, ONLY: radius
    USE iers_ser
    USE spher_harm
! ======INPUT======================================================================
    INTEGER,          INTENT(IN) :: ndf              ! dimension of position vector  
    DOUBLE PRECISION, INTENT(IN) :: x(ndf),v(ndf),t0 ! position, velocity, time     
    INTEGER,          INTENT(IN) :: ips              ! flag for recomputation
    INTEGER,          INTENT(IN) :: imem             ! memory location          
! ======OUTPUT======================================================================
    DOUBLE PRECISION, INTENT(OUT) :: f(ndf)   ! acceleration     
    INTEGER,          INTENT(OUT) :: idc      ! integer flag
    DOUBLE PRECISION, INTENT(OUT) :: xxpla(6) ! positions and vel of the planet involved in close-app 
                                              ! stored only if idc.ne.0
! ======END INTERFACE===============================================================
! control of derivatives                                                
!    integer ide 
!    common/deriv/ide 
! ==================================================================================
! J2 of the Earth
    DOUBLE PRECISION, PARAMETER :: pj2=-1.0826267d-3
! gravitational self energy of the Sun
    DOUBLE PRECISION, PARAMETER :: omega0=-3.52d-6
! for spherical harmonics 
    DOUBLE PRECISION, DIMENSION(3) :: y                            ! body-fixed position
    DOUBLE PRECISION               :: rot(3,3),rot1(3,3),rot2(3,3) ! rotation matrices to BF
    DOUBLE PRECISION               :: rott(3,3)                    ! rotation matrix from BF
    LOGICAL                        :: partials                     ! to compute derivatives for spher harm
    INTEGER                        :: lmax,lmi                     ! max number of spherical harmonics, index
    DOUBLE PRECISION               :: plmi                         ! coefficient of spherical harmonics
! temporaries to store harmonics                                       
    DOUBLE PRECISION :: army(3,4,narmx), accx(3),accy(3),dxeq(3),dxecl(3),armaccy(3,4)
! conversion of times 
    DOUBLE PRECISION               :: sec1
    INTEGER                        :: mjd1
    DOUBLE PRECISION, DIMENSION(3) :: temp ! to produce time series
! flag for outgassing model
    INTEGER :: outgas_sym !(def=0, 1=symm model, 2=asymm model)
! other variables
    DOUBLE PRECISION :: r(nmassx),d(nmassx),rast  ! distances 
    DOUBLE PRECISION :: xpla(6,nmassx)            ! positions of the planets also vel.            
    DOUBLE PRECISION :: drgr(3,7),frel(3)         ! relativistic perturbations
    DOUBLE PRECISION :: yarkv(21),yarkvs(21)      ! Yarkovsky force, diurnal and seasonal
    DOUBLE PRECISION :: accrad(3)                 ! solar radiation pressure
    DOUBLE PRECISION :: secacc(3)                 ! secular perturbation along the velocity
    DOUBLE PRECISION :: secacc_der(3,6)           ! derivatives of secular perturbation along the velocity
    DOUBLE PRECISION :: nongrav(3,6)              ! non-gravitational forces for outgassing (Marsden model)
                                                  ! columns 4-6 are the partials of col. 1,2,3 w.r.t.
                                                  ! the delay parameter
    INTEGER          :: i,j,k,ir,ic,kk            ! loop indexes i=1,npla j=1,3
    DOUBLE PRECISION :: sum,var1                  ! scalar temporaries
    DOUBLE PRECISION :: enne,semiaxis,drift,vsize 
    DOUBLE PRECISION :: derf(3,6),dfb(3,3),dfc(3,3),dfd(3,3) 
! ==================================================================================
    outgas_sym=0
! ===========================================================
! JPL ephemerides, and asteroid ephemerides if required                 
    CALL planast(t0,ips,imem,velo_req,xpla) 
! ===========================================================           
! Computation of planet vector lengths                                  
    DO i=1,nmass
       r(i)=sqrt(xpla(1,i)**2+xpla(2,i)**2+xpla(3,i)**2) 
   ENDDO
   rast=sqrt(x(1)**2+x(2)**2+x(3)**2)
! ===========================================================           
! Computation of planets-asteroid distances                             
   idc=0 
   DO 30 i=1,nmass 
      d(i)=sqrt((x(1)-xpla(1,i))**2+(x(2)-xpla(2,i))**2+              &
     &             (x(3)-xpla(3,i))**2)                                 
      if(iclap.gt.0)then 
         if(d(i).lt.dmin(i))then 
            if(idc.eq.0)then 
               idc=i 
               xxpla(1:3)=xpla(1:3,i) 
               IF(velo_req)xxpla(4:6)=xpla(4:6,i) 
            else 
               write(*,*)' force: this should not happen',t0, idc,i 
               write(*,*)nmass,(d(j),dmin(j),j=1,nmass) 
               stop 
            endif
         endif
      endif
30 ENDDO
! ===========================================================           
! initialize force and derivatives                                        
   f=0.d0 
   derf=0.d0
! ===========================================================           
! general relativistic correction                                       
   if(icrel.eq.1)then
      call genrel(x,v,drgr) 
      temp=drgr(1:3,1)
      f(1:3)=f(1:3)+drgr(1:3,1) 
   elseif(icrel.eq.2)then
      call eihrel2(x,v,xpla,d,r,rast,frel)!,t0) 
      temp=frel(1:3)
      f(1:3)=f(1:3)+frel(1:3) 
! if we need the refined relativity then also include J_2 for the Sun   
      call j2sun(x,frel) 
      temp=temp+frel(1:3)
      f(1:3)=f(1:3)+frel(1:3) 
   endif
! ===========================================================
! direct radiation pressure, if A/M is available (or being estimated)
   IF(dyn%nmod.EQ.1.AND.dyn%active(1))THEN
      CALL radpast(x(1:3),accrad)
      f(1:3)=f(1:3)+accrad*dyn%dp(1)
   END IF
! Yarkovsky effect modeled by seiempirical parameter
! 22.11.2016 LDimare: added secacc_der
   IF(dyn%nmod.EQ.1.AND.dyn%active(2))THEN
      CALL secular_nongrav(x(1:3),v(1:3),secacc,secacc_der)
      temp=secacc*dyn%dp(2)
      f(1:3)=f(1:3)+secacc*dyn%dp(2)
      derf=derf+secacc_der(1:3,1:6)*dyn%dp(2)
   ENDIF
! ===========================================================
! Non gravitational contribution for comet like effects
! Symmetric outgassing model from Marsden formalism (1973)
! Asymmetric outgassing model from Yeomans and Chodas (1989)
! ===========================================================
! computation of force components
   nongrav=0.d0
   IF(dyn%nmod.EQ.2)THEN
      IF(dyn%active(1).OR.dyn%active(2).OR.dyn%active(3))THEN
         outgas_sym=1
      END IF
      IF(dyn%active(4))THEN
         outgas_sym=2
      END IF
      IF(outgas_sym.EQ.1)THEN
         CALL comet_non_grav_symm(x,v,nongrav) ! Call to symmetric model
         DO kk=1,3
            f(1:3)=f(1:3)+nongrav(1:3,kk)*dyn%dp(kk)
         ENDDO
      ELSEIF(outgas_sym.EQ.2)THEN
         CALL comet_non_grav_asymm(x,v,nongrav) ! Call to asymmetric model
         DO kk=1,3
            f(1:3)=f(1:3)+nongrav(1:3,kk)*dyn%dp(kk)
         ENDDO
      ENDIF
   END IF
! ===========================================================           
! Yarkovsky effect, if required and data are available (not being estimated)
   IF(iyark.GE.1.AND.iyark.LE.2.AND.yarfil)THEN
      IF(.not.yarini)THEN 
         WRITE(*,*)' contradiction in non gravitational parameters' 
         WRITE(*,*)' iyark=',iyark,' yarfil=',yarfil,' yarini=',yarini 
         STOP 
      ENDIF
! diurnal                                                               
      CALL yarkdi(x(1:3),yarkv,iyarpt)
      temp=yarkv(1:3)
      f(1:3)=f(1:3)+yarkv(1:3) 
! seasonal                                                              
      IF(iyark.EQ.2)THEN 
         CALL yarkse(x(1:3),v(1:3),yarkvs,iyarpt) 
         temp=temp+yarkvs(1:3)
         f(1:3)=f(1:3)+yarkvs(1:3) 
      ENDIF  
   END IF

! ===========================================================           
! Computation of indirect force FI
   temp=0.d0
   DO i=1,nmass
      IF(dyn%nmod.EQ.3)THEN
         f(1:3)=f(1:3)-gm(i)/r(i)**3*xpla(1:3,i)*(1.d0+dyn%dp(1)*omega0)
         temp=temp-gm(i)/r(i)**3*xpla(1:3,i)*(1.d0+dyn%dp(1)*omega0)
      ELSE
         f(1:3)=f(1:3)-gm(i)/r(i)**3*xpla(1:3,i) 
         temp=temp-gm(i)/r(i)**3*xpla(1:3,i)
      ENDIF
   ENDDO      
! ===========================================================           
! Adding planets-asteroid attractions                                   
   temp=0.d0
   DO i=1,nmass
      f(1:3)=f(1:3)+gm(i)/d(i)**3*(xpla(1:3,i)-x(1:3)) 
      temp=temp+gm(i)/d(i)**3*(xpla(1:3,i)-x(1:3)) 
! only for Earth, and for asteroid orbits
      IF(imerc.eq.1.and.i.eq.3)THEN
! if distance less than 0.1 AU, that is 23.5 Earth radii
! Special case: (29075) 1950 DA
!         IF(d(i).lt.1d-2)THEN 
! Original version
         IF(d(i).lt.1d-1)THEN ! new limit
! place radius, mass of planet in spher\_harm.mod
! should be consistent with gravity model
            partials=.TRUE.
            r_planet=reau
            gm_planet=gmearth
            mjd1=t0
            sec1=(t0-mjd1)*86400.d0
            CALL rotsys('ECLM',mj2000,s2000,'BF  ',mjd1,sec1,rot,rot1,rot2,1)
            rott=TRANSPOSE(rot)
            dxecl=x(1:3)-xpla(1:3,i)
            y=MATMUL(rot,dxecl)
! spherical harmonics
            CALL parm10(y,2,army,partials)
! adding up spherical harmonics            
            armaccy=0.d0
            armaccy(1:3,1)=armaccy(1:3,1)+pj2*army(1:3,1,5)
            IF(partials)THEN
               armaccy(1:3,2:4)=armaccy(1:3,2:4)+pj2*army(1:3,2:4,5)
            ENDIF
            f(1:3)=f(1:3)+MATMUL(rott,armaccy(1:3,1))
            !WRITE(37,170) t0,MATMUL(rott,armaccy(1:3,1)),dxecl
            ! derivatives
            IF(partials)THEN
               derf(1:3,1:3)=derf(1:3,1:3)+ &
                    MATMUL(rott,MATMUL(armaccy(1:3,2:4),rot))
            ENDIF
         ENDIF
      ENDIF
   ENDDO
! ===========================================================           
! Adding solar attraction                                               
   f(1:3)=f(1:3)-gm0*x(1:3)/rast**3
   temp=-gm0*x(1:3)/rast**3
! ===========================================================
! VARIATIONAL EQUATION
   IF(ndf.eq.3)RETURN 
! Computation of partial derivatives matrix                             
   do 1 ir=1,3 
      do 2 ic=ir,3 
         var1=3.d0*gm0*x(ir)*x(ic)/rast**5 
         if(ir.eq.ic)var1=var1-(gm0/(rast**3)) 
         sum=0.d0 
         DO i=1,nmass 
            sum=sum+3.d0*gm(i)*                                        &
     &             (xpla(ir,i)-x(ir))*(xpla(ic,i)-x(ic))/(d(i)**5)      
            if(ir.eq.ic)sum=sum-gm(i)/(d(i)**3) 
         ENDDO 
         derf(ir,ic)=derf(ir,ic)+var1+sum 
         if(ir.ne.ic) then 
            derf(ic,ir)=derf(ir,ic) 
         endif
2     enddo
1  enddo
! ===========================================================           
! Computation of variational equations second part                      
   DO 24 j=1,3 
      DO 25 k=1,3 
         dfb(j,k)=0.d0 
         dfc(j,k)=0.d0 
         DO i=1,3 
            ! 22.11.2016 LDimare: added terms wrt velocity (for nongrav. pert.)
            dfb(j,k)=dfb(j,k)+derf(j,i)*x(i+3*k)+derf(j,i+3)*v(i+3*k) 
            dfc(j,k)=dfc(j,k)+derf(j,i)*x(i+3*k+9)+derf(j,i+3)*v(i+3*k+9) 
         ENDDO
25    ENDDO 
24 ENDDO 
   DO j=1,3 
      DO k=1,3 
         f(j+3*k)=dfb(j,k) 
         f(j+3*k+9)=dfc(j,k) 
      ENDDO 
   ENDDO 
   IF(ndf.le.21)RETURN
! variational equation for nongrav
   IF(ndf.ne.21+3*nls)THEN
      WRITE(*,*) ' force: wrong dimensions, ndf,nls ',ndf,nls
      STOP
   ENDIF
   IF(sep_viol)THEN
      WRITE(*,*) ' force: variational eq. not ready for sep_viol ',sep_viol
      STOP
   ENDIF   
! right hand side of variational eq. for determination of dynpar
   DO kk=1,nls   
! Direct radiation pressure and Yarkovsky
! 22.11.2016 Added missing terms, LDimare
      IF(dyn%nmod.EQ.1)THEN
         IF(ls(kk).eq.1)THEN
            CALL prodmv(dfd(1:3,kk),derf(1:3,1:3),x(21+3*(kk-1)+1:21+3*kk))
            dfd(1:3,kk)=dfd(1:3,kk)+MATMUL(derf(1:3,4:6),v(21+3*(kk-1)+1:21+3*kk))
            dfd(1:3,kk)=dfd(1:3,kk)+accrad
            f(21+3*(kk-1)+1:21+3*kk)=dfd(1:3,kk)
         ELSEIF(ls(kk).eq.2)THEN
            CALL prodmv(dfd(1:3,kk),derf(1:3,1:3),x(21+3*(kk-1)+1:21+3*kk))
            dfd(1:3,kk)=dfd(1:3,kk)+MATMUL(derf(1:3,4:6),v(21+3*(kk-1)+1:21+3*kk))
            dfd(1:3,kk)=dfd(1:3,kk)+secacc
            f(21+3*(kk-1)+1:21+3*kk)=dfd(1:3,kk)
         ELSE
            WRITE(*,*)' error in controls, nls, ls, kk ', nls, ls(1:nls),kk
            STOP 
         ENDIF
! Outgassing model, cometary style
      ELSEIF(dyn%nmod.EQ.2)THEN
         IF(ls(kk).GT.0.AND.ls(kk).LE.3)THEN
            CALL prodmv(dfd(1:3,kk),derf(1:3,1:3),x(21+3*(kk-1)+1:21+3*kk))
            dfd(1:3,kk)=dfd(1:3,kk)+nongrav(1:3,ls(kk))
            f(21+3*(kk-1)+1:21+3*kk)=dfd(1:3,kk)
         ELSEIF(ls(kk).EQ.4.AND.outgas_sym.EQ.2)THEN
            CALL prodmv(dfd(1:3,kk),derf(1:3,1:3),x(21+3*(kk-1)+1:21+3*kk))
            dfd(1:3,kk)=dfd(1:3,kk)+nongrav(1:3,4)*dyn%dp(1)+&
                 & nongrav(1:3,5)*dyn%dp(2)+nongrav(1:3,6)*dyn%dp(3)
            f(21+3*(kk-1)+1:21+3*kk)=dfd(1:3,kk)
         ENDIF
      ENDIF
   ENDDO
 CONTAINS

! ======================================                                
! PLANAST                                                               
! subroutine providing planets and asteroids                            
! in ecliptic coordinates at given time t0                                                      
! ======================================                                
   SUBROUTINE planast(t0,ips,imem,velo_intrp,xpla) 
     USE fund_const
! input                                                                 
     DOUBLE PRECISION, INTENT(IN) :: t0 !   time MJD  
! flag for recomputation, memory location, flag for velocities, no. bo
     INTEGER, INTENT(IN) :: ips,imem
     INTEGER istate
     LOGICAL, INTENT(IN) :: velo_intrp 
! output                                                                
     DOUBLE PRECISION,INTENT(OUT) :: xpla(6,nmassx) 
! hidden input:  JPL header
     include 'jplhdr.h90' 
! end interface
! Workspace for state call
     DOUBLE PRECISION et0(2),rot(3,3),pv(6,12),pnut(4),xx(3),xxp(3) 
     INTEGER, PARAMETER:: memx=30 ! stored planets array 
     DOUBLE PRECISION xp(6,nmassx,memx) 
     DOUBLE PRECISION xast(3,nbepx),vast(3,nbepx) ! positions of the massive asteroids
     INTEGER i,j,k,ia ! loop indexes; ia=1,nast 
!****************                                                       
     SAVE ! static memory required 
! ===========================================================           
! reference system rotation matrix; from reference_systems.mod  
     IF(velo_intrp)THEN
        istate=2
     ELSE
        istate=1
     ENDIF
! computations to be done only if ips>0                                 
     IF(ips.eq.0)THEN                                                  
! Read planetary positions from JPL files                               
        et0(1)=2400000.5d0                                             
        et0(2)=t0                                                      
        bary=.false.                                                   
        CALL state(et0,listpl,pv,pnut,istate)                          
! reorder data (see pleph)                                              
        DO 9 i=1,npla                                    
           IF(itarg(i).eq.13)THEN                                  
              if(istate.eq.2)xxp(1:3)=pv(4:6,3)                       
              xx(1:3)=pv(1:3,3)
           ELSEIF(itarg(i).eq.10)THEN   
              if(istate.eq.2) xxp(1:3)=pv(4:6,10)*emrat/(1.d0+emrat)+pv(4:6,3)
              xx(1:3)=pv(1:3,10)*emrat/(1.d0+emrat)+pv(1:3,3)
           ELSEIF(itarg(i).eq.3)THEN  
              if(istate.eq.2)xxp(1:3)=pv(4:6,3)-pv(4:6,10)/(1.d0+emrat)
              xx(1:3)=pv(1:3,3)-pv(1:3,10)/(1.d0+emrat)  
           ELSE
              xx(1:3)=pv(1:3,itarg(i))
              if(istate.eq.2)xxp(1:3)=pv(4:6,itarg(i))
           ENDIF
! Change of reference system EQUM00 ---> ECLM00
           xp(1:3,i,imem)=MATMUL(roteqec,xx)
           IF(istate.eq.2)xp(4:6,i,imem)=MATMUL(roteqec,xxp)
9       ENDDO
! ===========================================================           
        IF(iatrue.gt.0)THEN                                            
! read asteroid positions from binary ephemerides                       
           CALL rdbep(t0,iatrue,astid,xast,vast)                       
! stacking selected asteroids in the planets array                      
           xp(1:3,npla+1:npla+iatrue,imem)=xast(1:3,1:iatrue)
           IF(istate.eq.2)xp(4:6,npla+1:npla+iatrue,imem)=vast(1:3,1:iatrue)
        ENDIF
     ENDIF
! ========================================================              
! copy into output array                                                
     xpla(1:3,1:nmass)=xp(1:3,1:nmass,imem)
     if(istate.eq.2)xpla(4:6,1:nmass)=xp(4:6,1:nmass,imem)
   END SUBROUTINE planast
                                                                        
! ***************************************************************       
! J2SUN                                                                 
! ***************************************************************       
! This subroutine computes the J_2 acceleration on an asteroid          
! due to the Sun; J_2(Sun) = 2 x 10^-7 according to JPL DE405.          
! Neglects tilt of the solar spin axis to the ecliptic.                 
! ***************************************************************       
   SUBROUTINE j2sun(x,accj2) 
! planetary masses and other model parameters                           
     DOUBLE PRECISION solj2,radsun 
     PARAMETER (solj2=2.d-7,radsun=4.6527174d-3) 
! scalars                                                               
     DOUBLE PRECISION xsun2,xsun5,ratio,brac1,brac2 
! vectors                                                               
     DOUBLE PRECISION x(3),accj2(3) 
! function                                                              
     DOUBLE PRECISION prscal 
! ---------------------------------------------------------------       
     xsun2=prscal(x,x) 
     xsun5=xsun2*xsun2*dsqrt(xsun2) 
     ratio=1.5d0*gm0*solj2*(radsun*radsun/xsun2)/xsun5 
     brac1=5.d0*x(3)*x(3)-xsun2 
     brac2=5.d0*x(3)*x(3)-3.d0*xsun2 
     accj2(1)=ratio*x(1)*brac1 
     accj2(2)=ratio*x(2)*brac1 
     accj2(3)=ratio*x(3)*brac2
!     accj2=0.d0
!      write(*,*)accj2(1),accj2(2)                                      
     RETURN 
   END SUBROUTINE j2sun
                                                                
! *******************************************************************  
! GENREL
! *******************************************************************                                                   
! This routine treates PN terms in the kinematically                   
! nonrotating frame according to the scheme given by Damour,           
! Soffel and Xu (DSX IV, Phys Rev D, Jan 1994).                        
!                                                                      
! only the Sun monopole term is accepted                                
! output format: 1) drgr(i,1) ... components of the acceleration       
!                2) drgr(i,j+1) ... derivatives with respect to the    
!                                 initial conditions                   
!                                                                      
!                                    D. Vokrouhlick\'y, 1/3/94         
!  Modified for ORBFIT, Milani & Baccili 1997                          
!  computation of derivatives disabled                                 
! *******************************************************************  
   SUBROUTINE genrel(x,vs,drgr) 
! vlight is required  
     USE fund_const
! position, velocity, relativistic effects                              
     DOUBLE PRECISION x(3),vs(3),drgr(3,7) 
! intermediate for partials                                             
!     double precision dmgrx(3,6)                                       
! scalar temporaries                                                    
     DOUBLE PRECISION rsate2,rsate,xv,v2,c2,eafac,fac,brac 
! loop indexes                                                          
     INTEGER i                                                        
! --------------------------------------------------------------------- 
     c2=vlight*vlight 
!                                                                       
     rsate2=x(1)*x(1)+x(2)*x(2)+x(3)*x(3) 
     rsate=dsqrt(rsate2) 
     xv=x(1)*vs(1)+x(2)*vs(2)+x(3)*vs(3) 
     v2=vs(1)*vs(1)+vs(2)*vs(2)+vs(3)*vs(3) 
!  this refer to the Sun, not the earth                                 
     eafac=gm0/rsate 
!                                                                       
!  Internal terms:                                                      
!     --  monopole (DSX IV, 3.13)                                       
!         (rem. also accepted as the IERS Standard)                     
     brac=4.d0*eafac-v2 
     fac=eafac/rsate2/c2 
!  acceleration and partial derivatives with respect to the X and V     
     drgr(1:3,1)=fac*(brac*x(1:3)+4.d0*xv*vs(1:3)) 
! --------------------------------------------------------------------- 
!  partial derivatives available for later use                          
!     do 21 i=1,3                                                       
!       do 22 j=1,3                                                     
!         dmgrx(i,j)=fac*(4.d0*vs(i)*(vs(j)-3.d0*xv*x(j)/rsate2)-       
!    .       (3.d0*brac+4.d0*eafac)*x(i)*x(j)/rsate2)                   
!         dmgrx(i,j+3)=2.d0*fac*(2.d0*vs(i)*x(j)-x(i)*vs(j))            
!         if (i.eq.j) then                                              
!            dmgrx(i,j)=dmgrx(i,j)+brac*fac                             
!            dmgrx(i,j+3)=dmgrx(i,j+3)+4.d0*xv*fac                      
!         endif                                                         
!22     enddo                                                           
!21   enddo                                                             
!                                                                       
!     do 31 i=1,3                                                       
!       drgr(i,2)=dmgrx(i,1)                                            
!       drgr(i,3)=dmgrx(i,2)                                            
!       drgr(i,4)=dmgrx(i,3)                                            
!       do 30 j=4,6                                                     
!         drgr(i,j+1)=dmgrx(i,j)                                        
!30     enddo                                                           
!31   enddo                                                             
! --------------------------------------------------------------------- 
!                                                                       
     RETURN 
   END SUBROUTINE genrel

! ***************************************************************       
! EIHREL                                                                
! ***************************************************************       
! This subroutine computes heliocentric acceleration of the             
! asteroid on the (1/c^2)(= PN) level. Both the direct solar term       
! (~Schwarzschild contribution) and the planetary contributions         
! are included. The monopole terms considered only (= EIH               
! approximation).                                                       
!                                                                       
! Written by D Vokrouhlicky and S Chesley, Nov 3, 1999                  
! ***************************************************************       
   SUBROUTINE eihrel(x,vs,xpla,d,rpla,xsun,drgr) 
! vlight is required  
      USE fund_const
! position, velocity, relativistic effects                              
      double precision x(3),vs(3),xpla(6,nmassx),d(nmassx),drgr(3) 
      double precision rpla(nmassx) 
! scalar temporaries                                                    
      double precision temp,rplapla,strqua,potall,potalls,gmall,tempp 
      double precision xsun,xsun3,xdsun2,vast2,vsun2,rplapla2,scals 
      double precision vlight2 
! loop indexes                                                          
      integer i,ii,j,ib 
      double precision g(3),vsun(3),vast(3),vpla(3,nmassx),tbtvec(3) 
      double precision xastpla(3),nastpla(3),vastpla(3) 
      double precision rplavec(3),tmpvec(3),rpla3(nmassx) 
! --------------------------------------------------------------------- 
! initialization                                                        
      drgr=0.d0 
! compute potential of all bodies on asteroid, planets on Sun           
! and heliocentric distances of planets and the total mass of the       
! system                                                                
      potall=gm0/xsun 
      gmall=gm0 
      potalls=0.d0 
      do i=1,npla 
         potall=potall+gm(i)/d(i) 
         rpla3(i)=rpla(i)**3 
         potalls=potalls+gm(i)/rpla(i) 
         gmall=gmall+gm(i) 
      enddo 
! compute barycentric velocities of sun, asteroid, and planets          
      do j=1,3 
         vsun(j)=0.d0 
         do i=1,npla 
            vsun(j)=vsun(j)-gm(i)*xpla(j+3,i)/gmall
         enddo 
         vast(j)=vsun(j)+vs(j) 
      enddo 
      vsun2=vsun(1)*vsun(1)+vsun(2)*vsun(2)+vsun(3)*vsun(3) 
      vast2=vast(1)*vast(1)+vast(2)*vast(2)+vast(3)*vast(3) 
      do i=1,npla 
         do j=1,3 
            vpla(j,i)=vsun(j)+xpla(3+j,i) 
         enddo 
      enddo 
! ................................................................      
! compute indirect (1/c^2) planetary accelerations:                     
      do ib=1,npla 
         g(1:3)=-gm(ib)*xpla(1:3,ib)/rpla3(ib)
!        strange quantity                                               
         strqua=1.5d0*gm0/rpla(ib) 
         do ii=1,npla 
            if(ii.ne.ib) then 
               rplavec(1:3)=xpla(1:3,ii)-xpla(1:3,ib) 
               rplapla2=rplavec(1)**2+rplavec(2)**2+rplavec(3)**2 
               rplapla=sqrt(rplapla2) 
               scals=xpla(1,ib)*rplavec(1)+xpla(2,ib)*rplavec(2)+       &
     &               xpla(3,ib)*rplavec(3)                              
               strqua=strqua+gm(ii)*(1.d0-0.5d0*scals/rplapla2)/rplapla 
            endif 
         enddo 
         tempp=(xpla(1,ib)*vpla(1,ib)+xpla(2,ib)*vpla(2,ib)+            &
     &          xpla(3,ib)*vpla(3,ib))/rpla(ib)                         
         temp=2.d0*(xpla(4,ib)*xpla(4,ib)+xpla(5,ib)*xpla(5,ib)+        &
     &        xpla(6,ib)*xpla(6,ib))-vsun2-1.5d0*tempp*tempp-           &
     &        4.d0*potalls-strqua                                       
! add "g*temp" terms                                                    
         do i=1,3 
            drgr(i)=drgr(i)+temp*g(i) 
         enddo 
! add "velocity-dependent" terms                                        
         do i=1,3 
            tmpvec(i)=-4.d0*xpla(3+i,ib)+vpla(i,ib) 
         enddo 
         temp=g(1)*tmpvec(1)+g(2)*tmpvec(2)+g(3)*tmpvec(3) 
         do i=1,3 
            drgr(i)=drgr(i)+temp*xpla(3+i,ib) 
         enddo 
! "third-body" term (tbt)                                               
         do i=1,3 
            tbtvec(i)=gm0*xpla(i,ib)/rpla3(ib) 
         enddo 
         do ii=1,npla 
            if(ii.ne.ib) then 
               do j=1,3 
                  rplavec(j)=xpla(j,ib)-xpla(j,ii) 
               enddo 
               rplapla=dsqrt(rplavec(1)**2+rplavec(2)**2+rplavec(3)**2) 
               temp=gm(ii)/(rplapla**3) 
               do i=1,3 
                  tbtvec(i)=tbtvec(i)+temp*rplavec(i) 
               enddo 
            endif 
         enddo 
         temp=3.5d0*gm(ib)/rpla(ib) 
         do i=1,3 
            drgr(i)=drgr(i)+temp*tbtvec(i) 
         enddo 
      enddo 
! compute direct (1/c^2) accelerations:                                 
! -- planetary terms                                                    
      do ib=1,npla 
         do j=1,3 
            g(j)=gm(ib)*(x(j)-xpla(j,ib))/(d(ib)**3) 
         enddo 
!        pos rel to planet (unit vector)                                
         do j=1,3 
            xastpla(j)=x(j)-xpla(j,ib) 
            nastpla(j)=xastpla(j)/d(ib) 
         enddo 
!        vast rel to planet                                             
         do j=1,3 
            vastpla(j)=vast(j)-vpla(j,ib) 
         enddo 
!        strange quantity                                               
         scals=xastpla(1)*xpla(1,ib)+xastpla(2)*xpla(2,ib)+             &
     &         xastpla(3)*xpla(3,ib)                                    
         strqua=gm0*(1.d0-0.5d0*scals/(rpla(ib)**2))/rpla(ib) 
         do ii=1,npla 
            if(ii.ne.ib)then 
               do j=1,3 
                  rplavec(j)=xpla(j,ii)-xpla(j,ib) 
               enddo 
               rplapla2=rplavec(1)**2+rplavec(2)**2+rplavec(3)**2 
               rplapla=dsqrt(rplapla2) 
               scals=xastpla(1)*rplavec(1)+xastpla(2)*rplavec(2)+       &
     &               xastpla(3)*rplavec(3)                              
               strqua=strqua+gm(ii)/rplapla*(1.d0+0.5d0*scals/rplapla2) 
            endif 
         enddo 
         tempp=nastpla(1)*vpla(1,ib)+nastpla(2)*vpla(2,ib)+             &
     &         nastpla(3)*vpla(3,ib)                                    
         temp=2.d0*(vastpla(1)**2+vastpla(2)**2+vastpla(3)**2)-         &
     &        vast2-1.5d0*tempp*tempp-4.d0*potall-strqua                
!        add g*temp to force                                            
         do i=1,3 
            drgr(i)=drgr(i)-g(i)*temp 
         enddo 
! add "velocity-dependent" terms                                        
         do i=1,3 
            tmpvec(i)=4.d0*vastpla(i)+vpla(i,ib) 
         enddo 
         temp=g(1)*tmpvec(1)+g(2)*tmpvec(2)+g(3)*tmpvec(3) 
         do i=1,3 
            drgr(i)=drgr(i)+vastpla(i)*temp 
         enddo 
! "third-body" term (tbt)                                               
         do i=1,3 
            tbtvec(i)=gm0*xpla(i,ib)/rpla3(ib) 
         enddo 
         do ii=1,npla 
            if(ii.ne.ib)then 
               do j=1,3 
                  rplavec(j)=xpla(j,ib)-xpla(j,ii) 
               enddo 
               rplapla=dsqrt(rplavec(1)**2+rplavec(2)**2+rplavec(3)**2) 
               temp=gm(ii)/(rplapla**3) 
               do i=1,3 
                  tbtvec(i)=tbtvec(i)+temp*rplavec(i) 
               enddo 
            endif 
         enddo 
         temp=-3.5d0*gm(ib)/d(ib) 
         do i=1,3 
            drgr(i)=drgr(i)+temp*tbtvec(i) 
         enddo 
      enddo 
! compute direct (1/c^2) accelerations:                                 
! -- solar (~ Schwarzschild) term                                       
      xsun3=xsun**3 
      do i=1,3 
         g(i)=gm0*x(i)/xsun3 
      enddo 
      xdsun2=vs(1)*vs(1)+vs(2)*vs(2)+vs(3)*vs(3) 
!     strange quantity                                                  
      strqua=0.d0 
      do i=1,npla 
         scals=x(1)*xpla(1,i)+x(2)*xpla(2,i)+x(3)*xpla(3,i) 
         strqua=strqua+gm(i)*(1.d0+0.5d0*scals/(rpla(i)**2))/rpla(i) 
      enddo 
      tempp=(x(1)*vsun(1)+x(2)*vsun(2)+x(3)*vsun(3))/xsun 
      temp=2.d0*xdsun2-vast2-1.5d0*tempp*tempp-4.d0*potall-strqua 
! add g*temp to force                                                   
      drgr(1:3)=drgr(1:3)-g(1:3)*temp 
! add velocity terms ...                                                
      tmpvec(1:3)=4.d0*vs(1:3)+vsun(1:3) 
      temp=DOT_PRODUCT(g(1:3),tmpvec(1:3))
      drgr(1:3)=drgr(1:3)+temp*vs(1:3) 
! "third-body" term (tbt)                                               
      do i=1,3 
         tbtvec(i)=0.d0 
      enddo 
      do ii=1,npla 
         temp=gm(ii)/rpla3(ii) 
         do i=1,3 
            tbtvec(i)=tbtvec(i)+temp*xpla(i,ii) 
         enddo 
      enddo 
      temp=3.5d0*gm0/xsun 
      do i=1,3 
         drgr(i)=drgr(i)+temp*tbtvec(i) 
      enddo 
! ..................................................................    
! scale by 1/c^2                                                        
      vlight2=vlight**2 
      do i=1,3 
         drgr(i)=drgr(i)/vlight2 
      enddo 
      return 
    END SUBROUTINE eihrel

! ***************************************************************       
! EIHREL2                                                                
! ***************************************************************       
! This subroutine computes heliocentric acceleration of the             
! asteroid on the (1/c^2)(= PN) level. Both the direct solar term       
! (~Schwarzschild contribution) and the planetary contributions         
! are included. The monopole terms considered only (= EIH               
! approximation).                                                       
! Equation (4-26) from 
! http://descanso.jpl.nasa.gov/Monograph/series2/Descanso2_all.pdf
! Written by. D Farnocchia, Feb 17, 2012                  
! ***************************************************************       
    SUBROUTINE eihrel2(x,v,xpla,d,rpla,rast,drgr)!,t) 
      USE fund_const
!!! BEGIN INTERFACE !!!
! heliocentric position and velocity of the asteorid 
      DOUBLE PRECISION,INTENT(IN),DIMENSION(3) :: x,v
! heliocentric position and velocity of the planets
      DOUBLE PRECISION,INTENT(IN),DIMENSION(6,nmassx) :: xpla
! distance between asteroid and planets
      DOUBLE PRECISION,INTENT(IN),DIMENSION(nmassx) :: d
! planet heliocentric distances
      DOUBLE PRECISION,INTENT(IN),DIMENSION(nmassx) :: rpla
! asteroid heliocentric distance
      DOUBLE PRECISION,INTENT(IN) :: rast
! relativistic acceleration
      DOUBLE PRECISION,INTENT(OUT),DIMENSION(3) :: drgr
! current time (only for time series)
!      DOUBLE PRECISION,INTENT(IN) :: t
!!! END INTERFACE !!!
! sun + planet masses 
      DOUBLE PRECISION :: gmall
! Newtonian acceleration on planets
      DOUBLE PRECISION,DIMENSION(3,nmassx) :: g
! Newtonian acceleration on the Sun
      DOUBLE PRECISION,DIMENSION(3) :: g0
! Baricentric velocity of Sun and asteroid, and their norm^2
      DOUBLE PRECISION,DIMENSION(3) :: vsun,vast
      DOUBLE PRECISION :: vsun2,vast2
! Baricentric velocity of planets
      DOUBLE PRECISION,DIMENSION(3,nmassx) :: vpla
! Auxiliary variable
      DOUBLE PRECISION :: term
! Function for norm of 3-dim vectors
      DOUBLE PRECISION :: vsize
! Loop indexes
      INTEGER :: j,k,l
! for time series
!      DOUBLE PRECISION :: temp(0:nmassx,3),drgr_ind(3)
! --------------------------------------------------------------------- 
! initialization                                                        
      drgr=0.d0
!      temp=0.d0
! compute barycentric velocities of sun, asteroid, and planets          
      vsun(1:3)=0.d0
      gmall=gm0
      DO k=1,npla
         gmall=gmall+gm(k)
      ENDDO
      DO k=1,npla 
         vsun(1:3)=vsun(1:3)-gm(k)*xpla(4:6,k)/gmall
      ENDDO
      vast(1:3)=vsun(1:3)+v(1:3) 
      vsun2=DOT_PRODUCT(vsun,vsun)
      vast2=DOT_PRODUCT(vast,vast)
      DO k=1,npla 
         vpla(1:3,k)=vsun(1:3)+xpla(4:6,k) 
      ENDDO 

! PLANETS       
! first term (first 3 rows)
      DO j=1,npla
         term=-4.d0*gm0/rast-gm0/rpla(j)
! first term, first row
         DO l=1,npla
            term=term-4.d0*gm(l)/d(l)
         ENDDO
! second term first row
         DO k=1,npla
            IF(k.EQ.j) CYCLE
            term=term-gm(k)/vsize(xpla(1:3,k)-xpla(1:3,j))
         ENDDO
! second row
         term=term+vast2+2.d0*DOT_PRODUCT(vpla(1:3,j),vpla(1:3,j))&
              -4.d0*DOT_PRODUCT(vast,vpla(1:3,j))
! third row
         g(1:3,j)=-gm0*xpla(1:3,j)/rpla(j)**3
         DO k=1,npla  ! acceleration on planet ib
            IF(k.EQ.j) CYCLE
            g(1:3,j)=g(1:3,j)+ &
                 gm(k)/vsize(xpla(1:3,k)-xpla(1:3,j))**3&
                 *(xpla(1:3,k)-xpla(1:3,j))
         ENDDO
         term=term-1.5d0*(DOT_PRODUCT(x(1:3)-&
              xpla(1:3,j),vpla(1:3,j))/d(j))**2 &
              +0.5d0*DOT_PRODUCT(xpla(1:3,j)-x(1:3),g(1:3,j))
         drgr=drgr+term*gm(j)*(xpla(1:3,j)-x(1:3))/d(j)**3
!         temp(j,:)=temp(j,:)+term*gm(j)*(xpla(1:3,j)-x(1:3))/d(j)**3
      END DO
! secon term (4th row)
      DO j=1,npla
         drgr=drgr+gm(j)/d(j)**3*&
              DOT_PRODUCT(x(1:3)-xpla(1:3,j),4.d0*vast(1:3)-&
              3.d0*vpla(1:3,j))*(vast(1:3)-vpla(1:3,j))
!         temp(j,:)=temp(j,:)+gm(j)/d(j)**3*&
!              DOT_PRODUCT(x(1:3)-xpla(1:3,j),4.d0*vast(1:3)-&
!              3.d0*vpla(1:3,j))*(vast(1:3)-vpla(1:3,j))
      ENDDO
! third term (5th row)
      DO j=1,npla
         drgr=drgr+3.5d0*gm(j)/d(j)*g(1:3,j)
!         temp(j,:)=temp(j,:)+3.5d0*gm(j)/d(j)*g(1:3,j)
      ENDDO

! SUN
! first term (first 3 rows)
      term=-4.d0*gm0/rast
! first term, first row
      DO l=1,npla
         term=term-4.d0*gm(l)/d(l)
      ENDDO
! second term first row
      DO k=1,npla
         term=term-gm(k)/rpla(k)
      ENDDO
! second row
      term=term+vast2+2.d0*vsun2&
              -4.d0*DOT_PRODUCT(vast,vsun)
! third row
      g0=0.d0
      DO k=1,npla  ! acceleration on sun
         g0=g0+gm(k)/rpla(k)**3*(xpla(1:3,k))
      ENDDO
      term=term-1.5d0*(DOT_PRODUCT(x(1:3),vsun)/rast)**2 &
              +0.5d0*DOT_PRODUCT(-x(1:3),g0)
      drgr=drgr-term*gm0*x(1:3)/rast**3.d0
!      temp(0,:)=temp(0,:)-term*gm0*x(1:3)/rast**3.d0
! secon term (4th row)
      drgr=drgr+gm0/rast**3*DOT_PRODUCT(x(1:3),4.d0*vast&
           -3.d0*vsun)*(vast-vsun)
!      temp(0,:)=temp(0,:)+gm0/rast**3*DOT_PRODUCT(x(1:3),4.d0*vast&
!           -3.d0*vsun)*(vast-vsun)
! third term (5th row)
      drgr=drgr+3.5d0*gm0/rast*g0
!      temp(0,:)=temp(0,:)+3.5d0*gm0/rast*g0
! dividing by vlight**2
      drgr=drgr/vlight**2.d0
!      temp=temp/vlight**2.d0
!! print time series
!!      168 FORMAT(I3,1X,F16.10,3(1X,E24.15))
!!      DO j=0,npla
!!         WRITE(34,168) j,t,temp(j,:)
!!      ENDDO


!! INDIRECT TERM ON THE SUN      
!! first term (first 3 rows)
!!      DO j=1,npla
!!         term=-gm0/rpla(j)
!! first term, first row
!!         DO l=1,npla
!!            term=term-4.d0*gm(l)/rpla(l)
!!         ENDDO
!! second term first row
!!         DO k=1,npla
!!            IF(k.EQ.j) CYCLE
!!            term=term-gm(k)/vsize(xpla(1:3,k)-xpla(1:3,j))
!!         ENDDO
!! second row
!!         term=term+vsun2+2.d0*DOT_PRODUCT(vpla(1:3,j),vpla(1:3,j))&
!!              -4.d0*DOT_PRODUCT(vsun,vpla(1:3,j))
!! third row
!!         g(1:3,j)=-gm0*xpla(1:3,j)/rpla(j)**3
!!         DO k=1,npla  ! acceleration on planet ib
!!            IF(k.EQ.j) CYCLE
!!            g(1:3,j)=g(1:3,j)+ &
!!                 gm(k)/vsize(xpla(1:3,k)-xpla(1:3,j))**3&
!!                 *(xpla(1:3,k)-xpla(1:3,j))
!!         ENDDO
!!         term=term-1.5d0*(DOT_PRODUCT(-&
!!              xpla(1:3,j),vpla(1:3,j))/rpla(j))**2 &
!!              +0.5d0*DOT_PRODUCT(xpla(1:3,j),g(1:3,j))
!!         drgr_ind=drgr_ind+term*gm(j)*(xpla(1:3,j))/rpla(j)**3
!!         temp(j,:)=temp(j,:)+term*gm(j)*(xpla(1:3,j)-x(1:3))/d(j)**3
!!      END DO
!! secon term (4th row)
!!      DO j=1,npla
!!         drgr_ind=drgr_ind+gm(j)/rpla(j)**3*&
!!              DOT_PRODUCT(-xpla(1:3,j),4.d0*vsun(1:3)-&
!!              3.d0*vpla(1:3,j))*(vsun(1:3)-vpla(1:3,j))
!!         temp(j,:)=temp(j,:)+gm(j)/d(j)**3*&
!!              DOT_PRODUCT(x(1:3)-xpla(1:3,j),4.d0*vast(1:3)-&
!!              3.d0*vpla(1:3,j))*(vast(1:3)-vpla(1:3,j))
!!      ENDDO
!! third term (5th row)
!!      DO j=1,npla
!!         drgr_ind=drgr_ind+3.5d0*gm(j)/rpla(j)*g(1:3,j)
!!         temp(j,:)=temp(j,:)+3.5d0*gm(j)/d(j)*g(1:3,j)
!!      ENDDO
!!      drgr_ind=drgr_ind/vlight**2.d0
!!      167 FORMAT(F16.10,3(1X,E24.15))
!!      WRITE(44,167) t,drgr_ind

    END SUBROUTINE eihrel2

  END SUBROUTINE force
  
END MODULE force_model
! ====================================================== 
! vel_required
! logical informing the propagator that velocity is required
! because relativistic perturbations are being computed
! ======================================================
LOGICAL FUNCTION velocity_req()
  USE force_model
  IMPLICIT NONE
  IF(icrel.gt.0)THEN
     velocity_req=.true.
  ELSE
     velocity_req=.false.
  ENDIF
  RETURN
END FUNCTION velocity_req
! ======================================================                
!  RADAR_OB                                                             
!  sets the radar flag if there is radar data, to be used               
!  by selmet to adjust the force model                                  
! ======================================================                
SUBROUTINE radar_ob(type,m) 
  USE force_model, ONLY: radar, norad_obs
  IMPLICIT NONE 
! INPUT                                                                 
  INTEGER m 
  CHARACTER*(1) type(m) 
! HIDDEN OUTPUT : in PUBLIC DATA  radar, norad_obs 
! END INTERFACE                                                         
  INTEGER j 
  radar=.false. 
  norad_obs=0
  DO j=1,m 
     IF(type(j).eq.'R'.or.type(j).eq.'V')THEN
        radar=.true. 
        norad_obs=norad_obs+1
     ENDIF
  ENDDO
  RETURN 
END SUBROUTINE radar_ob

!  self perturbation control 
SUBROUTINE selpert(name,found) 
  USE planet_masses
  USE force_model
  USE name_rules, ONLY: name_len
  IMPLICIT NONE 
  LOGICAL found 
  CHARACTER*(*) name 
  CHARACTER*(name_len) nam1 
  CHARACTER*30 string 
  INTEGER iabe,ln,ls,ia,iat
  LOGICAL first
  DATA first /.true./
  SAVE first
  IF(first)THEN
     CALL masjpl
     first=.false.
  ENDIF
!                                                                       
  found=.false. 
  IF(iast.eq.0)RETURN 
!                                                                       
  nam1=name 
  CALL rmsp(nam1,ln) 
  iat=0 
  do  ia=1,iast 
     IF(nam1.eq.astnam(ia))THEN
        WRITE(*,*)' self perturbation of ',nam1(1:ln),' avoided' 
        found=.true. 
     ELSE 
        iat=iat+1 
        astid(iat)=ia 
     ENDIF
  enddo
  iatrue=iat 
  RETURN 
201 FORMAT(1P,E18.10,1X,A) 
202 WRITE(*,*)'selpert: too many asteroids requested, iast=',iast 
  iast=ia-1 
  WRITE(*,*)'selpert: asteroids available ',iast 
  RETURN 
END  SUBROUTINE selpert

SUBROUTINE selpert2(nam0,namp,nfound) 
  USE force_model
  USE planet_masses, ONLY: iast
  IMPLICIT NONE 
  CHARACTER*(*) nam0,namp 
  LOGICAL found0,foundp 
  INTEGER nfound 
  CALL selpert(nam0,found0) 
  CALL selpert(namp,foundp) 
  IF(found0.and.foundp)THEN 
     WRITE(*,*)' please do not try to identify ',nam0,' with ', namp 
     WRITE(*,*)' All perturbations by massive asteroids disabled' 
     nfound=2 
     iast=0 
  ELSEIF(found0)THEN 
     WRITE(*,*)' you should not do identification with an' 
     WRITE(*,*)' asteroid with mass, such as ',nam0 
     WRITE(*,*)' perturbations by ',nam0,' disabled' 
     nfound=1 
  ELSEIF(foundp)THEN 
     WRITE(*,*)' you should not do identification with an' 
     WRITE(*,*)' asteroid with mass, such as ',namp 
     WRITE(*,*)' perturbations by ',namp,' disabled' 
     nfound=1 
  ELSE 
     nfound=0 
  ENDIF
  RETURN 
END SUBROUTINE selpert2

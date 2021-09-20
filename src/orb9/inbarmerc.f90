!  =========================================                            
!  inbarmerc                                                            
!                                                                       
!  computes the barycenter of sun plus mercury                          
!                                                                       
!  input: masses and elements of all the planets from JPL ephemerides   
!  output: barycenter of the inner planets, in EQUM00                   
!  =========================================                            
PROGRAM inbarmerc 
  USE fund_const
  USE massmod
  IMPLICIT NONE 
  INTEGER, PARAMETER ::  nplax=nbox-1
  INTEGER npla  ! actual number of planets
! arrays                                                                
  DOUBLE PRECISION x(6,nbox),y(6,nbox) 
  DOUBLE PRECISION bx(6),by(6),bz(6),xb(6,nbox) 
  DOUBLE PRECISION enne(nbox),gmi(4) 
! J2 and relativistic corrections                                       
  INTEGER nbodin,nplain 
!  DOUBLE PRECISION vlight 
!  names and numbers of input planets and asteroids, comments           
  character*9 nombar(1) 
  character*60 compla 
  character*100 colhea 
!  identifiers of coordinate systems                                    
  character*3 coox,cooy,sysx,sysy,unitx 
  character*3 cooxb,cooyb,sysxb,sysyb 
  character*6 refx,refy 
  character*6 refxb,refyb 
! times                                                                 
  INTEGER iday,month,iyear,ifla 
  DOUBLE PRECISION t0,dmjd,hour,dt 
! constants                                                             
  DOUBLE PRECISION gjyr,gjyr2 
! functions                                                             
  DOUBLE PRECISION tjm1 
! integer scalars                                                       
  INTEGER ncen,nbo,noutp,nbod 
! loop indexes                                                          
  INTEGER i,m,j,len 
! fail flag in reading JPL ephemerides
  LOGICAL fail
!  =========================================                            
!  arrays for consulting jpl software                                   
  INCLUDE 'jplhdr.h90' 
  DOUBLE PRECISION rout(6,nbox),et(2),elem(6,nbox) 
! planet names                                                          
  character*9 nomi(nplax) 
  data nomi/'MERCURY  ','VENUS    ','EARTHMOON','MARS     ',        &
     &    'JUPITER  ','SATURN   ','URANUS   ','NEPTUNE  ','PLUTO    '/  
!  ========================================                             
  character cmont*3,cmv(12)*3 
  data cmv/'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep',   &
     &         'Oct','Non','Dec'/                                       
!  ========================================                             
! library directory                                                     
  call libini 
!  (1) choice of date                                                   
  open(1,file='inbaric.opt',status='old') 
  read(1,*)ifla 
  if(ifla.eq.0)then 
     read(1,*)iday,month,iyear,hour 
     dmjd=tjm1(iday,month,iyear,hour) 
     cmont=cmv(month) 
  else 
     read(1,*)dmjd 
     call datjm1(dmjd,iday,month,iyear,hour,cmont) 
  endif
!   ========================================                            
!   (2)  convert from MJD to JD                                         
  et(1)=2400000.5d0 
  et(2)=dmjd 
  t0=et(1)+et(2) 
!   ========================================                            
!   (3) consult jpl ephemerides                                         
!     heliocentric positions and velocities                             
  ncen=11 
  nbod=10 
  npla=nbod-1 
  DO 10 i=1,npla 
!  we want earth-moon barycenter and not the earth                      
     if(i.eq.3)then 
        nbo=13 
     else 
        nbo=i 
     endif
     call dpleph(et,nbo,ncen,rout(1,i),2) 
10 ENDDO
!   ========================================                            
!   (4) copy masses, in au**3/day**2 in JPL ephem                       
!  gm(1)=cval(18) 
!  DO  m=1,npla 
!     gm(m+1)=cval(8+m) 
!  ENDDO
  gm=0.d0
  DO i=1,30 ! maybe more???
     IF(cnam(i).eq.'GMS')THEN
        gm(1)=cval(i)
     ELSEIF(cnam(i).eq.'GM1')THEN
        gm(2)=cval(i)
     ELSEIF(cnam(i).eq.'GM2')THEN
        gm(3)=cval(i)
     ELSEIF(cnam(i).eq.'GMB')THEN
        gm(4)=cval(i)
     ELSEIF(cnam(i).eq.'GM4')THEN
        gm(5)=cval(i)
     ELSEIF(cnam(i).eq.'GM5')THEN
        gm(6)=cval(i)
     ELSEIF(cnam(i).eq.'GM6')THEN
        gm(7)=cval(i)
     ELSEIF(cnam(i).eq.'GM7')THEN
        gm(8)=cval(i)
     ELSEIF(cnam(i).eq.'GM8')THEN
        gm(9)=cval(i)
     ELSEIF(cnam(i).eq.'GM9')THEN
        gm(10)=cval(i) ! actually Pluto's mass is not anymore used
     ENDIF
  ENDDO
  DO i=1,10
    IF(gm(i).eq.0.d0) fail=.true.
  ENDDO
  call masses(gm,nbox,pmu,rm,rmt,sm,smp) 
!   ========================================                            
!  (5) conversion to cartesian coordinates, heliocentric and ecliptic   
!     jpl data are cartesian, heliocentric, equatorial                  
  coox='CAR' 
  sysx='HEL' 
  refx='EQUM00' 
  cooy='CAR' 
  sysy='HEL' 
  refy='ECLM00' 
  call coord(rout,sysx,nbox,coox,refx,y,sysy,cooy,refy,enne,by) 
!   ========================================                            
!  (6) barycenter: computation                                          
  cooxb='CAR' 
  sysxb='HEL' 
  refxb='ECLM00' 
  cooyb='CAR' 
  sysyb='JAC' 
  refyb='ECLM00' 
! this is for excluding from the integration Mercury only               
  nbodin=2 
  call coord(y,sysxb,nbodin,cooxb,refxb,xb,sysyb,cooyb,refyb,enne,bx)
!   ========================================                            
!   (7) conversion to output system for all planets                     
  cooy='CAR' 
  sysy='HEL' 
  refy='ECLM00' 
  coox='KEP' 
  sysx='HEL' 
  refx='ECLM00' 
  call coord(y,sysy,nbox,cooy,refy,x,sysx,coox,refx,enne,bz) 
! ==========================================                            
!  compute size of short periodic perturbations                         
  call epsi(npla,x,enne) 
! ==============================================                        
!  (6d) inner planets J2 and relativistic corrections                   
! see Nobili et al., A&A 1989 vol. 210 pag313-336; page 316             
!    conversion factor to au**3/year**2                                 
  gjyr=365.25d0 
  gjyr2=gjyr*gjyr 
! astronomical unit (Km) taken from JPL ephemerides                     
!      au = 1.4959787069100000d13                                       
!     write(*,*)au                                                      
! speed of light (Km/s)                                                 
  vlight = 2.99792458d5 
! speed of light (au/y)                                                 
  vlight=vlight/au*8.64d5*365.25 
! small parameter of Schwarzschild potential correction                 
  scw=3.d0*(gm(1)*gjyr2/vlight)**2 
! fictitious J2 of the Sun (for standard radius of 1 AU)                
  cj2=0.d0 
  nplain=nbodin-1 
  DO m=1,nplain 
!       write(*,*)gm(m+1),x(1,m)                                        
     cj2=cj2+gm(m+1)*gjyr2*x(1,m)**2/2.d0 
  ENDDO
! ===============================================                       
!  (6)  output file barsunxx                                            
  open(2,file='barmerxx.inc',status='unknown') 
  compla=                                                           &
     &' JPL DE 431; vector from baryc.inner sol.system to the sun '     
  colhea=                                                           &
     &'   x    y    z     dx/dt   dy/dt  dz/dt'                         
  nombar(1)='BARI     ' 
  call wrihea(2,nombar,compla,colhea,cooxb,sysxb,refxb,'RAD',1,0,t0)
! masses of inner planets                                               
  write(2,107) 
  107 format('ninp=1; no. inner planets used, with inv.masses:')
  DO j=1,nbodin-1 
     gmi(j)=1.d0/rm(j) 
  ENDDO
  write(2,116)(gmi(j),j=1,nbodin-1) 
116 format(4f17.8) 
! output coefficients for J2 and  relativistic correction               
  write(2,117)cj2,scw 
117 format('; J2 and relativistic correction constants'/2d15.7/       &
     & '; Sun - barycenter of inner solar system , pos and vel')        
! convert velocities to au/yr from au/day                               
  DO i=4,6 
     bx(i)=bx(i)*3.6525d2 
  ENDDO
  write(2,106)bx 
106 format(3d18.10) 
  close(2) 
!   ========================================                            
!   (8) output all planets: header                                      
  open(2,file='allplmxx.inc',status='unknown') 
  DO j=1,npla 
     nompla(j)=nomi(j) 
  ENDDO
  write(compla,128)iday,cmont,iyear,hour 
128 format(' JPL DE 431 ',i3,1x,a3,1x,i4,1x,f8.4) 
  colhea= '  a(au)    e    I    Node  nrev   arg.peri nrev  an.med.  nrev '
  call wrihea(2,nompla,compla,colhea,coox,sysx,refx,'DEG',npla,0,t0)
!   ========================================                            
!   (9) output all planets: masses                                      
!    convert to au**3/year**2                                           
  gjyr=365.25d0 
  gjyr2=gjyr*gjyr 
  DO j=1,nbod 
     gm(j)=gm(j)*gjyr2 
  ENDDO
  call wrimas(2,nbox) 
  dt=0.d0 
  write(2,118)dt 
118 format(f12.8) 
!   ========================================                            
!   (10) output all planets: elements                                   
  call unrad(x,npla,coox,'DEG',dpig) 
  DO j=1,nplax 
     write(2,108)(x(i,j),i=1,4),0,x(5,j),0,x(6,j),0 
  ENDDO
108 format(f14.10,f13.10,f13.8,2(f13.7,1x,i2),f14.8,1x,i3) 
  close(2) 
!   ========================================                            
!   (11) output outer planets: header                                   
  open(2,file='pvenexx.inc',status='unknown') 
  noutp=7 
  DO j=1,noutp 
     nompla(j)=nomi(j+1) 
  ENDDO
  write(compla,128)iday,cmont,iyear,hour 
  colhea='  a(au)    e    I    Node  nrev   arg.peri nrev  an.med.  nrev '
  call wrihea(2,nompla,compla,colhea,coox,sysx,refx,'DEG',noutp,0,t0)
!   ========================================                            
!   (12) output outer planets: masses                                   
  nbod=8 
  DO j=2,nbod 
     gm(j)=gm(j+1) 
  ENDDO
  call wrimas(2,nbod) 
  dt=0.d0 
  write(2,118)dt 
!   ========================================                            
!   (13) output outer planets: elements                                 
!      call unrad(x,noutp,coox,'DEG',dpig)                              
  DO j=1,noutp 
     write(2,108)(x(i,j+1),i=1,4),0,x(5,j+1),0,x(6,j+1),0 
  ENDDO
  close(2) 
END PROGRAM inbarmerc

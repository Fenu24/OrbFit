!============================================== 
!   ANGLES:  f90 NORMALIZED
!      could become a module with
!      type(angle_var) and overloading of +,-,= 
! princ		principal value of an angle
! angupd        update of angle (princ+rev)     
! primea	mean of two angle variables
! pridif        difference of two angles  
! angvcf        conversion factor for angular velocities
! sessag        angle into sessagesimal
! polar         transformation from cartesian to polar coordinates 
!============================================== 
!           
! Computes the principal value of an angle      
!           
double precision function princ (a) 
  USE fund_const
  double precision, intent(in) :: a 
  princ=dmod(a,dpig) 
  if(princ.lt.0.d0)princ=princ+dpig 
END FUNCTION princ
!
!========================
!  primea
!  arithmetic mean of two angles,
!  expressed as principal values
double precision function primea(a,b)
  USE fund_const 
  implicit none 
  double precision a,b,princ 
  a=princ(a) 
  b=princ(b) 
  if(abs(b-a).gt.pig)then 
     primea=(a+b)/2.d0+pig 
     primea=princ(primea) 
  else 
     primea=(a+b)/2.d0 
  endif
END function primea
!========================
!  pridif
!  difference of two angles,
!  expressed as principal values
double precision function pridif(a,b) 
  USE fund_const 
  implicit none  
  double precision a,b,princ           
  a=princ(a) 
  b=princ(b) 
  pridif=a-b 
  if(pridif.gt.pig)then 
     pridif=pridif-dpig 
  elseif(pridif.lt.-pig)then 
     pridif=pridif+dpig 
  endif
END function pridif   
! **********************************************************
!  ANGUPD
!   given the principal value ang of an angle, and an angle vang        
!   with ng revolutions, the routine
!   finds a new new ng in the assumption that   
!   less than half a revolution has occurred; ang is modified to        
!   include the multiples of dpig   
!   in case ang is not between -pig and pig, e.g. because it is         
!   between 0 and dpig, it is first reduced to principal value          
SUBROUTINE angupd(ang,vang,ng)
  USE fund_const 
  IMPLICIT NONE 
  DOUBLE PRECISION ang,vang,d 
  INTEGER ng,ig 
! first reduce to principal value   
  if(ang.gt.pig)ang=ang-dpig 
  if(ang.lt.-pig)ang=ang+dpig 
! count revolutions     
  ig=nint((vang-ang)/dpig) 
  ang=ang+dpig*ig 
  d=ang-vang 
! correct by one revolution if necessary        
  if(d.gt.pig)then 
     ng=ng-1 
     ang=ang-dpig 
  elseif(d.lt.-pig)then 
     ng=ng+1 
     ang=ang+dpig 
  endif
END SUBROUTINE angupd
! Copyright (C) 2000 by Mario Carpino (carpino@brera.mi.astro.it)       
! Version: November 30, 2000        
! --------------------------------------------------------------------- 
!           
!  *****************************************************************    
!  *                                                               *    
!  *                         A N G V C F                           *    
!  *                                                               *    
!  *           Conversion factor for angular velocities            *    
!  *                                                               *    
!  *****************************************************************    
!           
! INPUT:    UNIT      -  Description of unit    
!           
! OUTPUT:   CF        -  Conversion factor from rad/d to "unit"         
!           ERROR     -  Error flag 
!
SUBROUTINE angvcf(unit,cf,error) 
  USE fund_const 
  IMPLICIT NONE           
  CHARACTER*(*) unit 
  DOUBLE PRECISION cf 
  LOGICAL error         
  DOUBLE PRECISION ct 
  CHARACTER*20 unang,untim 
  LOGICAL nospli         
  error=.true.
! Input string must be of the form 'unang/untim',           
! where: "unang" is the angular unit (rad/deg/arcmin/arcsec)
!        "untim" is the time unit (d/h/min/s);  
! for instance: 'deg/d' or 'arcsec/min'         
  untim=unit 
  CALL stspli(untim,'/',unang,nospli) 
  IF(nospli) RETURN 
  CALL locase(unang) 
  CALL locase(untim)           
! List of supported angular units   
  IF(unang.EQ.'rad') THEN 
     cf=1 
  ELSEIF(unang.EQ.'deg') THEN 
     cf=degrad 
  ELSEIF(unang.EQ.'arcmin') THEN 
     cf=degrad*60 
  ELSEIF(unang.EQ.'''') THEN 
     cf=degrad*60 
  ELSEIF(unang.EQ.'arcsec') THEN 
     cf=degrad*3600 
  ELSEIF(unang.EQ.'"') THEN 
     cf=degrad*3600 
  ELSE 
     RETURN 
  END IF
! List of supported time units      
  IF(untim.EQ.'d') THEN 
     ct=1 
  ELSEIF(untim.EQ.'day') THEN 
     ct=1 
  ELSEIF(untim.EQ.'h') THEN 
     ct=24 
  ELSEIF(untim.EQ.'hour') THEN 
     ct=24 
  ELSEIF(untim.EQ.'min') THEN 
     ct=24*60 
  ELSEIF(untim.EQ.'s') THEN 
     ct=24*3600 
  ELSE 
     RETURN 
  END IF
  cf=cf/ct 
  error=.false.           
END SUBROUTINE angvcf
! Copyright (C) 1996 by Mario Carpino (carpino@brera.mi.astro.it)       
! Version: August 27, 1996          
! --------------------------------------------------------------------- 
!           
!  *****************************************************************    
!  *                                                               *    
!  *                         S E S S A G                           *    
!  *                                                               *    
!  *       Transform an angle into sessagesimal notation           *    
!  *                                                               *    
!  *****************************************************************    
!           
! INPUT:    ANG       -  Angle      
!           
! OUTPUT:   SIGA      -  Sign       
!           INTA      -  Integer part           
!           MINA      -  Minutes    
!           SECA      -  Seconds    
!
SUBROUTINE sessag(ang,siga,inta,mina,seca) 
  IMPLICIT NONE           
  DOUBLE PRECISION ang,seca 
  INTEGER inta,mina,truncat 
  CHARACTER*1 siga         
  DOUBLE PRECISION anga,u 
  IF(ang.GE.0.d0) THEN 
     anga=ang 
     siga='+' 
  ELSE 
     anga=ABS(ang) 
     siga='-' 
  END IF          
  inta=truncat(anga,1d-10) 
  u=(anga-inta)*60.d0 
  u=ABS(u) 
  mina=truncat(u,1d-8) 
  seca=(u-mina)*60.d0           
END SUBROUTINE sessag 
! Copyright (C) 1998 by Mario Carpino (carpino@brera.mi.astro.it)       
! Version: December 11, 1998        
! --------------------------------------------------------------------- 
!           
!  *****************************************************************    
!  *                                                               *    
!  *                          P O L A R                            *    
!  *                                                               *    
!  *    Transformation from cartesian to polar coordinates         *    
!  *                                                               *    
!  *****************************************************************    
!           
! INPUT:    X         -  Cartesian coordinates  
!           
! OUTPUT:   ALPHA     -  Right ascension (rad)  
!           DELTA     -  Declination (rad)      
!           R         -  Distance   
!           
SUBROUTINE polar(x,alpha,delta,r) 
  USE fund_const
  IMPLICIT NONE           
  DOUBLE PRECISION x(3),alpha,delta,r 
  DOUBLE PRECISION send,cosd,sina,cosa     
  r=SQRT(x(1)**2+x(2)**2+x(3)**2) 
  IF(r.EQ.0.d0) THEN 
     alpha=0.d0 
     delta=0.d0 
     RETURN 
  END IF          
  send=x(3)/r 
  delta=ASIN(send) 
  cosd=COS(delta) 
  IF(cosd.EQ.0.d0) THEN 
     alpha=0.d0 
     RETURN 
  END IF          
  cosa=x(1)/(r*cosd) 
  sina=x(2)/(r*cosd) 
  alpha=ATAN2(sina,cosa) 
  IF(alpha.LT.0.d0) alpha=alpha+dpig 
END SUBROUTINE polar

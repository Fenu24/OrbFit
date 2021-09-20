! ==================================================                    
!  GENERAL PURPOSE COORDINATE CHANGE PACKAGE                            
!                                                                       
!  Vers. 2.0, Meudon, Septemmber 1991                                   
!                                                                       
!  Coordinate systems supported:\par                                    
!          'BAR'=barycentric j=1,nbod-1\par                             
!          'HEL'=heliocentric; since the Sun is assumed to              
!                be at the origin, j=1,nbod-1 and x(i,nbod) is dummy\par
!          'JAC'=jacobian; since the center of mass is assumed to       
!                be at the origin, j=1,nbod-1 and x(i,nbod) is dummy\par
!          'HEC'=heliocentric canonical: positions are heliocentric,    
!                velocitys are barycentric\par                          
!  Coordinate types supported: \par                                     
!          'CAR'= cartesian positions and velocities\par                
!          'KEP'= keplerian elements, a,e,I,Omega,omega,mean an.\par    
!          'EQU'= equinoctal elements a,h,k,p,q,lambda\par              
!                  all the angles in radians\par                        
!  Masses required: gm(i), i=1,nbod; the routine masses has to          
!           be called before the first call to coord/cooast to          
!           initialise the mass parameters \par                             
!  Reference systems supported:\par                                     
!      'EQUM00' mean equator and gamma point of J2000.0  \par           
!      'EQUM50' mean equator and gamma point of B1950.0  \par           
!      'ECLM50' mean ecliptic and gamma point of B1950.0 \par           
!      'ECLM00' mean ecliptic and gamma point of J2000.0 \par           
!      'INVL1B' invariable plane of  the outer planets                  
!                from LONGSTOP 1B\par                                   
! ================================================                      
!  {\bf  masses}                                                        
!  computes mass ratios for use of coord, cooast, coosys                
!  input:                                                               
!        gm(j), j=1,nbod (list of masses; first one is the sun)\par     
!  output:                                                              
!        pmu(j), j=1,nbod-1 (mass ratio of the jacobian binary j-1)\par 
!        rm(j), j=1,nbod-1 (mass ratio of planet j to mass of the sun)\p
!        rmt(j), j=1,nbod-1 (mass ratio of planet j to total mass)\par  
!        sm(j), j=1,nbod-1 (sum of all the masses up to j+1,to change to
!             elements in jacobian coordinates)\par                     
!        smp(j), j=1,nbod-1 (mass of j+1 plus mass of the sun)\par      
! all dimensioned nbod, the last component unused in pmu,rm,smp,rmt,sm  
! ================================================                      
SUBROUTINE masses(gm,nbod,pmu,rm,rmt,sm,smp) 
   IMPLICIT NONE
   INTEGER, INTENT(IN) :: nbod
   DOUBLE PRECISION, INTENT(IN)  ::  gm(nbod)
   DOUBLE PRECISION, INTENT(OUT)  ::  pmu(nbod),rm(nbod),sm(nbod),smp(nbod),rmt(nbod)
! end interface 
   INTEGER i
   DOUBLE PRECISION smtot
   sm(1)=gm(1)+gm(2) 
   DO i=2,nbod-1 
      sm(i)=sm(i-1)+gm(i+1)
   ENDDO 
   DO  i=2,nbod 
      smp(i-1)=gm(i)+gm(1) 
      rm(i-1)=gm(i)/gm(1) 
      pmu(i-1)=gm(i)/sm(i-1)
   ENDDO 
   smtot=sm(nbod-1) 
   DO i=1,nbod-1 
      rmt(i)=gm(i+1)/smtot 
   ENDDO
 END SUBROUTINE masses
! \vfill\eject                                                          
! ================================================                      
!  {\bf coord} vers. 2.0, Meudon, September 1991                        
!                                                                       
!   general purpose coordinate change for planets                       
!                                                                       
!  INPUT: x(i,j), i=1,6, j=1,nbod\par                                   
!         sysx,sysy= coordinate system\par                              
!         nbod= total number of bodies, Sun included; the number of     
!               vectors given both in input and output is nbod-1\par    
!         coox,cooy= coordinate types\par                               
!  OUTPUT: y(i,j), i=1,6, j=1,nbod-1 as specified by cooy, sysy \par    
!          bar(i), i=1,6 barycentric coordinates of the Sun\par         
!          enne(j), j=1,nbod-1 mean motions (not computed               
!                when coox=cooy='CAR'; take care that it                
!                refers to the output coordinate system                 
!                only if cooy.neq.'CAR')\par                            
! ===============================================                       
 SUBROUTINE coord(x,sysx,nbod,coox,refx,y,sysy,cooy,refy,enne,bar) 
   USE massmod
   IMPLICIT NONE
   INTEGER, INTENT(IN) :: nbod
! ==============================================                        
   DOUBLE PRECISION ::  x(6,nbod),y(6,nbod),enne(nbod) 
   DOUBLE PRECISION ::  bar(6),z(6,nbox) 
   character*3 coox,cooy,sysx,sysy 
   character*6 refx,refy 
! end interface
   INTEGER j
!  mass ratios etc.; this routine has been called already               
!c    call masses(gm,nbod,pmu,rm,rmt,sm,smp)                            
! ===============================================                       
!  transformation to cartesian coordinates of input x:                  
   if(coox.eq.'CAR')then 
      y(1:6,1:nbod-1)=x(1:6,1:nbod-1)
   else 
      if(sysx.eq.'BAR')then 
!  barycentric elements have sun mass as the mass of the center         
! WARNING: do we still agree?
         DO j=1,nbod-1 
            call coocha(x(1,j),coox,gm(1),y(1,j),'CAR',enne(j)) 
         ENDDO
      elseif(sysx.eq.'HEL'.or.sysx.eq.'HEC')then 
!  heliocentric elements use mass of the sun + mass of current planet   
         DO j=1,nbod-1 
            call coocha(x(1,j),coox,smp(j),y(1,j),'CAR',enne(j))
         ENDDO 
      elseif(sysx.eq.'JAC')then 
!  jacobian elements use the sum of the masses up to (and included) the 
!  current planet as the mass of the center                             
         DO j=1,nbod-1 
            call coocha(x(1,j),coox,sm(j),y(1,j),'CAR',enne(j)) 
         ENDDO
      endif
   endif
! ===============================================                       
!  change of coordinate system                                          
   call coosys(y,sysx,refx,nbod,pmu,rm,rmt,sysy,refy,z,bar) 
! ===============================================                       
!  transformation from cartesian coordinates of output y:               
   if(cooy.eq.'CAR')then 
      y(1:6,1:nbod-1)=z(1:6,1:nbod-1)
   else 
      if(sysy.eq.'BAR')then 
!  barycentric elements have Sun mass as the mass of the center         
         DO j=1,nbod-1 
            call coocha(z(1,j),'CAR',gm(1),y(1,j),cooy,enne(j)) 
         ENDDO
      elseif(sysy.eq.'HEL'.or.sysy.eq.'HEC')then 
!  heliocentric elements use mass of the sun + mass of current planet   
         DO j=1,nbod-1 
            call coocha(z(1,j),'CAR',smp(j),y(1,j),cooy,enne(j)) 
         ENDDO
      elseif(sysy.eq.'JAC')then 
!  jacobian elements use the sum of the masses up to (and included) the 
!  current planet as the mass of the center                             
         DO j=1,nbod-1 
            call coocha(z(1,j),'CAR',sm(j),y(1,j),cooy,enne(j)) 
         ENDDO
      endif
   endif
  END SUBROUTINE coord
! \vfill\eject                                                          
! ================================================                      
!  {\bf cooast} general purpose coordinate change for asteroids         
!                                                                       
!  INPUT: x(i,j), i=1,6, j=1,nast\par                                   
!         bar(i),i=1,6 vector from the center of mass to the sun,       
!           always in the input ref. system refx\par                    
!         sysx,sysy= coordinate system\par                              
!         nast= total number of asteroids\par                           
!         coox,cooy= coordinate types\par                               
!  OUTPUT: y(i,j), i=1,6, j=1,nast as specified by cooy, sysy \par      
!          enne(j), j=1,nast mean motions (not computed                 
!                when coox=cooy='CAR'; take care that it                
!                refers to the output coordinate system                 
!                only if cooy.neq.'CAR')\par                            
! ===============================================                       
  SUBROUTINE cooast(x,sysx,nast,coox,refx,bar,nbod,y,sysy,cooy,     &
     &          refy,enne)                                              
    USE massmod
    IMPLICIT NONE
! ==============================================                        
    INTEGER, INTENT(IN) :: nast, nbod
    DOUBLE PRECISION ::  x(6,nast),y(6,nast),enne(nast) 
    DOUBLE PRECISION ::  bar(6)
    character*3 coox,cooy,sysx,sysy 
    character*6 refx,refy 
! end interface
    DOUBLE PRECISION :: z(6,nastx) 
    INTEGER j
    DOUBLE PRECISION smtot
! ===============================================                       
!  case of trivial tranformation                                        
    if(sysx.eq.sysy.and.coox.eq.cooy.and.refx.eq.refy)then 
       y(1:6,1:nast)=x(1:6,1:nast)
       return 
    endif
!  transformation to cartesian coordinates of input x:                  
    if(coox.eq.'CAR')then 
       y(1:6,1:nast)=x(1:6,1:nast)
    else 
       if(sysx.eq.'BAR')then 
!  barycentric elements have total mass as the mass of the center       
          smtot=sm(nbod-1) 
          DO j=1,nast 
             call coocha(x(1,j),coox,smtot,y(1,j),'CAR',enne(j)) 
          ENDDO
       elseif(sysx.eq.'HEL'.or.sysx.eq.'JAC'.or.sysy.eq.'HEC')then 
!  heliocentric elements use mass of the sun                            
          DO j=1,nast 
             call coocha(x(1,j),coox,gm(1),y(1,j),'CAR',enne(j)) 
          ENDDO
       else 
          write(9,*)' coordinate system ',sysx,' not supported' 
          stop 
       endif
    endif
! ===============================================                       
!  change of coordinate system                                          
    if(sysx.eq.'HEL'.or.sysx.eq.'JAC')then 
!  input heliocentric (or jacobian, by def. the same for                
!  massless bodies):                                                    
       if(sysy.eq.'HEL'.or.sysy.eq.'JAC')then 
!  no transformation 
          z(1:6,1:nast)=y(1:6,1:nast)
       elseif(sysy.eq.'BAR')then 
!  from heliocentric to barycentric: add barycentric                    
!  coordinates of the sun
          DO j=1,nast
             z(1:6,j)=y(1:6,j)+bar(1:6)
          ENDDO
       elseif(sysy.eq.'HEC')then 
!  from heliocentric to hel.canonical: add barycentric                  
!  velocity of the sun to the velocity                                  
          DO j=1,nast 
             z(1:3,j)=y(1:3,j) 
             z(4:6,j)=y(4:6,j)+bar(4:6) 
          ENDDO
       endif
    elseif(sysx.eq.'BAR')then 
!  barycentric input                                                    
       if(sysy.eq.'BAR')then 
!  no transformation
          z(1:6,1:nast)=y(1:6,1:nast)
       elseif(sysy.eq.'HEL'.or.sysy.eq.'JAC')then 
!  from barycentric to heliocentric: subtract barycentric               
!  coordinates of the sun                                               
          DO j=1,nast 
             z(1:6,j)=y(1:6,j)-bar(1:6)
          ENDDO
       elseif(sysy.eq.'HEC')then 
!  from barycentric to hel.canonical: subtract barycentric              
!  position of the sun from the position                                
          DO j=1,nast 
             z(1:3,j)=y(1:3,j)-bar(1:3) 
             z(4:6,j)=y(4:6,j)
          ENDDO
       endif
    elseif(sysx.eq.'HEC')then 
!  heliocentric canonical input                                         
       if(sysy.eq.'HEC')then 
!  no transformation    
          z(1:6,1:nast)=y(1:6,1:nast)                        
       elseif(sysy.eq.'HEL'.or.sysy.eq.'JAC')then 
!  from hel.canonical to heliocentric: subtract barycentric             
!  velocity of the sun from velocity                                    
          DO j=1,nast 
             z(1:3,j)=y(1:3,j) 
             z(4:6,j)=y(4:6,j)-bar(1:3)
          ENDDO
       elseif(sysy.eq.'BAR')then 
!  from hel.canonical to barycentric: add barycentric                   
!  position of the sun to the position                                  
          DO j=1,nast 
             z(4:6,j)=y(4:6,j) 
             z(1:3,j)=y(1:3,j)+bar(1:3)
          ENDDO 
       endif
    endif
! ===============================================                       
!  change of reference system                                           
    call refer(z,nast,refx,z,refy) 
! ===============================================                       
!  transformation from cartesian coordinates of output y:               
    if(cooy.eq.'CAR')then 
       y(1:6,1:nast)=z(1:6,1:nast)     
    else 
       if(sysy.eq.'BAR')then 
!  barycentric elements have total mass as the mass of the center       
          smtot=sm(nbod-1) 
          DO j=1,nast 
             call coocha(z(1,j),'CAR',smtot,y(1,j),cooy,enne(j)) 
          ENDDO
       elseif(sysy.eq.'HEL'.or.sysy.eq.'JAC'.or.sysy.eq.'HEC')then 
!  heliocentric elements use mass of the sun                            
          DO j=1,nast 
             call coocha(z(1,j),'CAR',gm(1),y(1,j),cooy,enne(j)) 
          ENDDO
       else 
          write(9,*)' coordinate system ',sysy,' not supported' 
          stop 
       endif
    endif
    return 
  END SUBROUTINE cooast
! ====================================================                  
!  {\bf coosys}                                                         
!                                                                       
!  general purpose coordinate system change routine \par                
!  given a set of nbod vectors x in a coordinate system coox:\par       
!     'BAR': barycentric \par                                           
!     'HEL': heliocentric \par                                          
!     'HEC': heliocentric canonical \par                                
!     'JAC': jacobian \par                                              
!  and the set of the masses, outputs y  according to cooy.\par         
!  There are always nbod-1 vectors, nbod always includes the sun.       
!                                                                       
!  The mass vectors (see routine 'masses') are:\par                     
!  pmu (mass ratio of the jacobian binaries)\par                        
!  rm (mass ratio of planets to the sun)\par                            
!  rmt (mass ratio of planets to total mass)\par                        
!  all dimensioned nbod with the last component unused.\par             
!  the arrays have to be dimensioned (6,nbod); however, the last        
!  vector is not used, therefore dimensions (6,nbod-1) should not       
!  create problems.\par                                                 
!  b in output it is always the barycentic coordinates of the Sun       
!  in the input reference system                                        
! ==============================================                        
  SUBROUTINE coosys(x,coox,refx,nbod,pmu,rm,rmt,cooy,refy,y,b) 
    IMPLICIT NONE  
    INTEGER, INTENt(IN) :: nbod
    DOUBLE PRECISION ::  y(6,nbod),x(6,nbod),b(6),pmu(nbod) 
    DOUBLE PRECISION ::  rmt(nbod),rm(nbod) 
    character*3 coox,cooy 
    character*6 refx,refy
! end interface
    INTEGER :: i,n
! =============================================                         
!   jacobian input:                                                     
    if(coox.eq.'JAC')then 
!  jacobian coordinates are untangled into heliocentric\par             
!  after step n, b is the barycenter of the bodies up to planet n;      
!  at the end, the  vector in b goes from the Sun to the center of mass 
!  and has to be changed of sign to output the Sun's baryc. vector\par  
!  pmu(n) is the mass ratio of the jacobian vector ending in            
!  planet number n (the index are one less than Milani and Nobili 1983);
! $$x_j=y_j+b_{j-1} \eqspa b_j=b_{j-1}+\mu_j y_j $$                     
! $$ \mu_j=\frac{m_j}{N_j} \eqspa N_j=\sum_{i=0}^{j}m_i $$             
       b(1:6)=0.d0 
       DO n=1,nbod-1 
          y(1:6,n)=x(1:6,n)+b(1:6) 
          b(1:6)=pmu(n)*x(1:6,n)+b(1:6)
       ENDDO 
       b=-b
! =============================================                         
!   heliocentric input:                                                 
    elseif(coox.eq.'HEL')then 
!   the first transformation is the identity, the center of mass has to 
!   be found:                                                           
! $$ b=-b_{N-1}=-\sum_{j=1}^{N-1} \frac{m_{j}}{M_{tot}}x_j $$           
       DO n=1,nbod-1 
          y(1:6,n)=x(1:6,n)
       ENDDO
       DO 14 i=1,6 
          b(i)=0.d0 
          DO n=1,nbod-1 
             b(i)=b(i)-rmt(n)*x(i,n)
          ENDDO
14     ENDDO
! ===============================================                       
!   barycentric input:                                                  
    elseif(coox.eq.'BAR')then 
!  b = barycentric coordinates of the Sun; if x are barycentric         
! $$ b=-\sum_{j=1}^{N-1} \frac{m_{j}}{m_1}x_j $$                        
       DO 26 i=1,6 
          b(i)=0.d0 
          DO n=1,nbod-1 
             b(i)=b(i)-rm(n)*x(i,n)
          ENDDO
26     ENDDO
!   the vectors are reduced to heliocentric, by subtracting b           
       DO 17 n=1,nbod-1 
          DO i=1,6 
             y(i,n)=x(i,n)-b(i)
          ENDDO
17     ENDDO
    elseif(coox.eq.'HEC')then 
!  the barycenter has to be computed with heliocentric formula          
!  for position, barycentric formula for velocities                     
       DO 27 i=1,3 
          b(i)=0.d0 
          b(i+3)=0.d0 
          DO n=1,nbod-1 
             b(i)=b(i)-rmt(n)*x(i,n) 
             b(i+3)=b(i+3)-rm(n)*x(i+3,n)
          ENDDO
27     ENDDO
!  velocities are reduced to heliocentric                               
       DO 28 i=1,3 
          DO n=1,nbod-1 
             y(i,n)=x(i,n) 
             y(i+3,n)=x(i+3,n)-b(i+3)
          ENDDO
28     ENDDO
    else 
       write(9,*)'******** coosys: input coord.sys ',coox,' unknown' 
       stop 
    endif
! ================================================                      
!   now y is anyway heliocentric, b is the barycentric  vector          
!   of the Sun                                                          
! ================================================                      
!   barycentric output:                                                 
    if(cooy.eq.'BAR')then 
!   if output is barycentric, b has to be added to all                  
       DO n=1,nbod-1 
          y(1:6,n)=y(1:6,n)+b(1:6)
       ENDDO
! ================================================                      
!   jacobian output:                                                    
    elseif(cooy.eq.'JAC')then 
       b=0.d0 
       do 22 n=1,nbod-1 
           DO i=1,6 
             y(i,n)=y(i,n)-b(i) 
             b(i)=pmu(n)*y(i,n)+b(i)
          ENDDO
22     ENDDO
       b=-b
! ================================================                      
!   heliocentric output: nothing to do                                  
    elseif(cooy.eq.'HEL')then 
!   heliocentric canonical output:                                      
    elseif(cooy.eq.'HEC')then 
       DO n=1,nbod-1 
          y(4:6,n)=y(4:6,n)+b(4:6)
       ENDDO 
    else 
       write(9,*)'******** coosys: out coord.sys ',cooy,' unknown' 
       stop 
    endif
!  ===============================================                      
!   change of reference system from refx to refy                        
    CALL refer(y,nbod-1,refx,y,refy) 
  END SUBROUTINE coosys
! \vfill\eject                                                          
! *************************************************                     
!  read/write inverse masses and unit conversion                        
! ==========================================================            
!  {\bf wrimas} writes planetary (inverse) masses                       
!  as used                                                              
SUBROUTINE wrimas(iun,nbod) 
  USE massmod
  USE fund_const
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: iun, nbod
! =================================================                     
  DOUBLE PRECISION ::  gmout(nbox), gjyr, gjyr2, g
  INTEGER i
!  Gauss constant from fund_const.mod      gk=0.01720209895d0 
!  conversion to internal units: 1AU, 1JYR=365.25 d(Julian year)        
  gjyr=365.25d0 
  gjyr2=gjyr*gjyr 
  g=gms*gjyr2 
  DO i=1,nbod 
     gmout(i)=gm(i)/g 
     gmout(i)=1.d0/gmout(i) 
  ENDDO
  write(iun,100) 
100 format('; planetary inverse masses as used') 
  write(iun,101)(gmout(i),i=1,nbod) 
101 format(1p,3d24.16) 
END SUBROUTINE wrimas
! ==========================================================            
! {\bf reamas} reads planetary (inverse) masses and                     
! computes auxiliary ratios, reduced masses, etc.                       
SUBROUTINE reamas(iun,nbod,ibar,pm,npin) 
  USE massmod
  USE fund_const
  IMPLICIT NONE        
  INTEGER, INTENT(IN) :: iun, nbod, ibar, npin
! end interface
  character*1 cc 
  DOUBLE PRECISION ::  pm(npin),gjyr,gjyr2,g
  INTEGER j,i  
! ====================================                                  
  read(iun,100)cc 
100 format(a1) 
  read(iun,*,end=3,err=3)(gm(j),j=1,nbod) 
  gm(1)=1.d0/gm(1) 
  if(ibar.ne.0)then 
     DO i=1,npin 
        gm(1)=gm(1)+1/pm(i)
     ENDDO
  endif
  DO j=2,nbod 
     gm(j)=1.d0/gm(j)
  ENDDO
!  Gauss constant: from fund_const.mod gk, gms=gk*gk
!  conversion to internal units: 1AU, 1JYR=365.25 d(Julian year)        
  gjyr=365.25d0 
  gjyr2=gjyr*gjyr 
  g=gms*gjyr2 
  DO i=1,nbod 
     gm(i)=gm(i)*g
  ENDDO
  call masses(gm,nbod,pmu,rm,rmt,sm,smp) 
  return 
3 write(9,*)' error in mass input, gm=',gm 
  stop 
END SUBROUTINE reamas
!  \vfill\eject                                                         
! =======================================================               
!  {\bf refer} change of reference system\par                           
!  input: vectors x(1:6,j), j=1:nn, in refx\par                         
!  output: vectors y(1:6,j), j=1:nn, in refy\par                        
!  x and y in the same memory position is allowed \par                  
!  reference systems supported:\par                                     
!  EQUM00 mean equator and gamma point of J2000.0  \par                 
!  EQUM50 mean equator and gamma point of B1950.0  \par                 
!  ECLM50 mean ecliptic and gamma point of B1950.0 \par                 
!  ECLM00 mean ecliptic and gamma point of J2000.0 \par                 
!  INVL1B invariable plane of  the outer planets from LONGSTOP 1B\par   
!  the transformations are performed along the tree: \par               
!  ECLM00 -- EQUM00 $<$-a- EQUM50 -- ECLM50 \par                        
!  INVL1B $<$-b- EQUM00 \par                                            
!  WARNING: this routine is optimized to perform many coordinate        
!  changes within the same run; static memory is essential              
!  (remember the -static flag in compilers with default                 
!  memory mode -automatic, e.g. MIPS)                                   
 SUBROUTINE refer(x,nn,refx,y,refy)
   IMPLICIT NONE 
   INTEGER, INTENT(IN) :: nn
   CHARACTER*6, INTENT(IN) :: refx,refy 
   DOUBLE PRECISION, INTENT(IN) ::  x(6,nn)
   DOUBLE PRECISION, INTENT(INOUT) :: y(6,nn)
! end interface 
   DOUBLE PRECISION ::  a(3,3),b(3,3),z(6),w(6)
   DOUBLE PRECISION ::  eps50, eps00, ra, dec, cra, sra,cde,sde, tm, tmv
   INTEGER lflag, i, j, n
   SAVE a,b,eps50,eps00,lflag 
! ===========================================                           
!   transformation matrix a is given:                                   
!   precession matrix from 1950.0 to 2000.0 from  Murray,C.A. (1989),   
!   Astron. Astrophys. 218, 325-329; this is X(0), equation (29),       
!   page 328 (given by column);                                         
   data a/                                                           &
     & 0.9999256794956877d0, 0.0111814832391717d0, 0.0048590037723143d0,&
     &-0.0111814832204662d0, 0.9999374848933135d0,-0.0000271702937440d0,&
     &-0.0048590038153592d0,-0.0000271625947142d0, 0.9999881946023742d0/
!  flag for first call                                                  
   data lflag/0/ 
   IF(lflag.eq.0)THEN 
      lflag=1 
!   transformation matrix b is computed at the first call:              
!   rotation from equator J2000 to invariable plane LONGSTOP 1B         
      ra=273.852730d0*(4.d0*datan(1.d0))/1.8d2 
      dec=66.991493d0*(4.d0*datan(1.d0))/1.8d2 
      cra=dcos(ra) 
      sra=dsin(ra) 
      cde=dcos(dec) 
      sde=dsin(dec) 
      b(1,1)=-sra 
      b(1,2)=cra 
      b(1,3)=0.d0 
      b(2,1)=-cra*sde 
      b(2,2)=-sra*sde 
      b(2,3)=cde 
      b(3,1)=cra*cde 
      b(3,2)=sra*cde 
      b(3,3)=sde 
!   obliquity of the mean equator B1950.0 (from recommendations         
!   from D.K.Yeomans, December 5, 1990, Appendix E; also Explanatory    
!   Supplement to the Ephemeris p. 98 and 1983 Astronomical Almanac     
!   p. B6)                                                              
      eps50=23.44578787d0*atan(1.d0)/4.5d1 
!   obliquity of the mean equator J2000.0 (from Yeomans cit.; also      
!   1990 Astronomical Almanac (p. K6, L2)                               
      eps00=23.43929111d0*atan(1.d0)/4.5d1 
   ENDIF 
! =====================================================                 
!  computations to be done at each call, for each object:               
   if(refx.eq.refy)then 
!  no change    
      y(1:6,1:nn)=x(1:6,1:nn)
      return 
   endif
   do 1 n=1,nn 
! ============================================                          
!  first step: everything is transformed to EQUM00                      
      if(refx.eq.'EQUM00')then 
!  first step is identity
         w(1:6)=x(1:6,n) 
      elseif(refx.eq.'ECLM00')then 
!  from ecliptic to equator of 2000                                     
         call rotyz(x(1:6,n),w,eps00) 
      elseif(refx.eq.'INVL1B')then 
!  from invariable plane to equator 2000                                
         do 4 i=1,3 
            tm=0.d0 
            tmv=0.d0 
            DO j=1,3 
               tm=tm+b(j,i)*x(j,n) 
               tmv=tmv+b(j,i)*x(j+3,n) 
            ENDDO
            w(i)=tm 
            w(i+3)=tmv
4        ENDDO
      elseif(refx.eq.'ECLM50')then 
!  from ecliptic to equator 1950                                        
         call rotyz(x(1:6,n),z,eps50) 
!   rotation to mean equator of 2000.0                                  
         DO 11 i=1,3 
            tm=0.d0 
            tmv=0.d0 
            DO j=1,3 
               tm=tm+a(i,j)*z(j) 
               tmv=tmv+a(i,j)*z(j+3)
            ENDDO
            w(i)=tm 
            w(i+3)=tmv
11       ENDDO
      elseif(refx.eq.'EQUM50')then 
!   rotation to mean equator of 2000.0                                  
         DO 13 i=1,3 
            tm=0.d0 
            tmv=0.d0 
            DO j=1,3 
               tm=tm+a(i,j)*x(j,n) 
               tmv=tmv+a(i,j)*x(j+3,n)
            ENDDO
            w(i)=tm 
            w(i+3)=tmv
13       ENDDO
      else 
         write(*,*)' reference system ',refx,' not supported' 
         stop 
      endif
! ==================================================                    
!  second step: now w is EQUM00, transformation to refy                 
      if(refy.eq.'EQUM00')then 
!  second step is identity
         y(1:6,n)=w(1:6) 
      elseif(refy.eq.'ECLM00')then 
!  from equator to ecliptic of 2000                                     
         call rotyz(w,y(1:6,n),-eps00) 
      elseif(refy.eq.'INVL1B')then 
!  from equator 2000 to invariable plane                                
         DO 20 i=1,3 
            tm=0.d0 
            tmv=0.d0 
            DO j=1,3 
               tm=tm+b(i,j)*w(j) 
               tmv=tmv+b(i,j)*w(j+3)
            ENDDO
            y(i,n)=tm 
            y(i+3,n)=tmv
20       ENDDO
      elseif(refy.eq.'ECLM50')then 
!   rotation to mean equator of 1950.0                                  
         DO 31 i=1,3 
            tm=0.d0 
            tmv=0.d0 
            DO j=1,3 
               tm=tm+a(j,i)*w(j) 
               tmv=tmv+a(j,i)*w(j+3) 
            ENDDO
            z(i)=tm 
            z(i+3)=tmv
31       ENDDO
!  from  equator to ecliptic 1950                                       
         call rotyz(z,y(1:6,n),-eps50) 
      elseif(refy.eq.'EQUM50')then 
!   rotation to mean equator of 1950.0                                  
         DO 33 i=1,3 
            tm=0.d0 
            tmv=0.d0 
            DO j=1,3 
               tm=tm+a(j,i)*w(j) 
               tmv=tmv+a(j,i)*w(j+3)
            ENDDO
            y(i,n)=tm 
            y(i+3,n)=tmv 
33       ENDDO
      else 
         write(*,*)' reference system ',refy,' not supported' 
         stop 
      endif
1  ENDDO
   return 
 END SUBROUTINE refer
! ========================================================              
!  {\bf rotyz}                                                          
! rotation around x axis by theta                                       
! not to be used with x,y in same address                               
SUBROUTINE rotyz(x,y,theta) 
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT(IN) :: x(6), theta
  DOUBLE PRECISION, INTENT(OUT) :: y(6) 
! end interface
  DOUBLE PRECISION :: c,s 
  c=dcos(theta) 
  s=dsin(theta) 
  y(2)=c*x(2)-s*x(3) 
  y(3)=s*x(2)+c*x(3) 
  y(5)=c*x(5)-s*x(6) 
  y(6)=s*x(5)+c*x(6) 
  y(1)=x(1) 
  y(4)=x(4) 
END SUBROUTINE rotyz

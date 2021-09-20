! ========================================================              
!  secthe vers. 8 (May 16, 1997)                                        
!   secular synthetic theory for the outer planets                      
!   from LONGSTOP 1B (Nobili, Milani & Carpino, A&A 1989)               
!   however for numerical integration tests, the frequencies            
!   are rederived from numerical integrations (p1, 20My)                
!   for semimajor axes: use of mean values from numerical               
!   integration (p1, 20 My)                                             
! ========================================================              
SUBROUTINE secthe(t,t0,pm,pa,pni,pg,pl,pmi,ps,pd,nnl,fnl,pninl,phanl)
  USE fund_const
  USE controlmod
  IMPLICIT NONE
! INPUT
  DOUBLE PRECISION, INTENT(IN) :: t, t0 ! current time shift, reference time (JD)
! OUTPUT
  INTEGER, INTENT(OUT) :: nnl ! number of nonlinear terms
  DOUBLE PRECISION, INTENT(OUT) :: pni(8,8),pg(8),pl(8),pmi(8,8),ps(8),pd(8),pa(8),pm(8)
  DOUBLE PRECISION, INTENT(OUT) :: fnl(nnlx),phanl(nnlx),pninl(2,nnlx)
! END INTERFACE
!  frequencies (arcsec/yr), amplitudes (adimensional), phases (deg)     
!  fundamental frequencies in arcsec/yr (from LONGSTOP)                 
  DOUBLE PRECISION g(8), s(8)
  data g(5),g(6),g(7),g(8),s(6),s(7),s(8)/                          &
     &   4.25749319d0,28.24552984,3.08675577,0.67255084,                &
     &   -26.34496354d0,-2.99266093,-0.69251386/                        
!  phases in degrees
  DOUBLE PRECISION  ag(8,8),fag(8),as(8,8),fas(8)
  data fag(5),fag(6),fag(7),fag(8),fas(6),fas(7),fas(8)/            &
     &   27.0005d0,124.1994d0,117.0516d0,70.7508d0,                     &
     &   303.9699d0,316.7946d0,200.3717d0/                              
!  amplitudes for jupiter                                               
  data ag(5,5),ag(5,6),ag(5,7),ag(5,8),as(5,6),as(5,7),as(5,8)/     &
     &   4.41872d-2,-1.57002d-2,1.8139d-3,5.8d-5,                       &
     &   3.15329d-3,-4.8478d-4,-5.8448d-4/                              
!  amplitudes for saturn                                                
  data ag(6,5),ag(6,6),ag(6,7),ag(6,8)                              &
     &   /3.29583d-2,4.82093d-2,1.5113d-3,5.75d-5/                      
  data as(6,6),as(6,7),as(6,8)/-7.85795d-3,-3.9446d-4,-5.6391d-4/ 
!   amplitudes for uranus                                               
  data ag(7,5),ag(7,6),ag(7,7),ag(7,8)/-3.75866d-2,-1.5471d-3,      &
     &          2.9033d-2,1.6665d-3/                                    
  data as(7,6),as(7,7),as(7,8)/3.5339d-4,8.88675d-3,5.4325d-4/ 
!   amplitudes for neptune                                              
  data ag(8,5),ag(8,6),ag(8,7),ag(8,8)/1.88143d-3,-1.0309d-4,       &
     &        -3.69711d-3,9.11787d-3/                                   
  data as(8,6),as(8,7),as(8,8)/3.782d-5,-1.06159d-3,5.78995d-3/ 
!  semimajor axes  - mean LONGSTOP                                      
  DOUBLE PRECISION aj0,as0,au0,an0
  data aj0,as0,au0,an0/5.2025696d0,9.5456831d0,19.194220d0,30.070998d0/
!  the s(5) term has by definition zero amplitude with respect to the   
!  invariable plane                                                     
  data s(5),as(5,5),as(6,5),as(7,5),as(8,5),fas(5)/6*0.d0/ 
!  amplitudes of the main nonlinear terms in the synthetic theory       
!  for Jupiter and Saturn 
  DOUBLE PRECISION        anl(2,nnlx)                                       
  data anl(1,1),anl(2,1)/-5.735d-4,1.9194d-3/ 
  data anl(1,2),anl(2,2)/-4.96d-5,1.666d-4/ 
  data anl(1,3),anl(2,3)/-2.95d-5,9.47d-5/ 
  data anl(1,4),anl(2,4)/1.982d-4,-6.052d-4/ 
  data anl(1,5),anl(2,5)/-1.936d-4,5.982d-4/ 
  data anl(1,6),anl(2,6)/-1.226d-4,3.792d-4/ 
  data anl(1,7),anl(2,7)/1.104d-4,-3.368d-4/ 
  data anl(1,8),anl(2,8)/-6.56d-5,-4.44d-5/ 
! LONGSTOP synthetic theory reference epoch, shifts
  DOUBLE PRECISION tlstp, deltjd, deltyr, tc
  data tlstp/2.4404005d6/
! scalar temporaries
  DOUBLE PRECISION summ, corr
! loop indexes
  INTEGER j,k,n
! startup
  INTEGER iflag 
  data iflag/0/ 
!  if(iflag.eq.0)then 
! time shift to align with LONGSTOP initial time
     deltjd=t0-tlstp 
     deltyr=deltjd/365.25d0 
!  masses of the outer planets ; mass of sun=1                          
     if(isec.eq.1)then 
!  standish & campbell values                                           
        pm(5)=1./1047.349e0 
        pm(6)=1./3497.915e0 
        pm(7)=1./22941.0e0 
        pm(8)=1./19432.0e0 
!  semimajor axes in AU                                                 
!  LONGSTOP 1B values                                                   
        pa(5)=aj0 
        pa(6)=as0 
        pa(7)=au0 
        pa(8)=an0 
!  correction from jacobian to heliocentric                             
        summ=1.d0 
        DO j=5,8 
           summ=summ+pm(j) 
           corr=summ/(1.d0+pm(j)) 
           pa(j)=pa(j)*corr
        ENDDO
     elseif(isec.eq.2)then 
!  De403 values                                                         
        pm(5)=1.d0/1047.349d0 
        pm(6)=1.d0/3497.898d0 
        pm(7)=1.d0/22902.98d0 
        pm(8)=1.d0/19412.24d0 
!  semimajor axes in AU                                                 
!  numerical integration value (p1, 20My)                               
        pa(5)=5.20259867695D+00 
        pa(6)=9.55493438259D+00 
        pa(7)=1.92184342840D+01 
        pa(8)=3.01104468747D+01 
!  fundamental frequencies in arcsec/yr                                 
!  numerical integration value (p1, 20My)                               
        g(5)=4.245016d0 
        g(6)=28.237937d0 
        g(7)=3.087924d0 
        g(8)=0.672952d0 
        s(6)=-26.339342d0 
        s(7)=-2.993419d0 
        s(8)=-0.692000d0 
     endif
!  frequencies, converted to rad/yr
     pg=0.d0
     ps=0.d0
     DO j=nplin,nplou 
        pg(j)=g(j)*radsec 
        ps(j)=s(j)*radsec 
     ENDDO 
!  nonlinear frequencies; the 8 biggest (non-chaotic) terms are used    
     nnl=8 
     fnl(1)=2.d0*pg(6)-pg(5) 
     fnl(2)=2.d0*pg(6)-pg(7) 
     fnl(3)=-2.d0*pg(5)+3.d0*pg(6) 
     fnl(4)=-pg(5)+pg(6)+pg(7) 
     fnl(5)=pg(5)+pg(6)-pg(7) 
     fnl(6)=-pg(5)+2.d0*pg(6)+ps(6)-ps(7) 
     fnl(7)=pg(5)-ps(6)+ps(7) 
!  mistake found by M.Moons, March 1994.                                
     fnl(8)=2.d0*pg(5)-pg(7) 
     DO n=1,nnl 
        pninl(1,n)=anl(1,n) 
        pninl(2,n)=anl(2,n)
     ENDDO 
!  adimensional amplitudes                                              
     DO j=nplin,nplou 
       DO k=nplin,nplou 
          pni(j,k)=ag(j,k) 
!  warning! longstop uses sin(i/2), yuasa uses sin(i)!!!!!              
!  this is an approximation because the series for                      
!  one variable cannot be converted into a series for                   
!  another variable term by term!                                       
          pmi(j,k)=2.d0*as(j,k)
       ENDDO
     ENDDO 
! end of startup
     iflag=1 
!  endif
! computations to be done at each call of this routine                  
  tc=t+deltyr 
! phases of the linear terms, converted to rad
! warning: the LONGSTOP phases are used anyway; maybe this should be corrected
! but the difference is small                                           
  DO j=nplin,nplou 
     pl(j)=fag(j)*radeg+tc*pg(j) 
     pd(j)=fas(j)*radeg+tc*ps(j)
  ENDDO
! phases of the nonlinear terms                                         
  phanl(1)=2.d0*pl(6)-pl(5) 
  phanl(2)=2.d0*pl(6)-pl(7) 
  phanl(3)=-2.d0*pl(5)+3.d0*pl(6) 
  phanl(4)=-pl(5)+pl(6)+pl(7) 
  phanl(5)=pl(5)+pl(6)-pl(7) 
  phanl(6)=-pl(5)+2.d0*pl(6)+pd(6)-pd(7) 
  phanl(7)=pl(5)-pd(6)+pd(7) 
  phanl(8)=2.d0*pl(5)-pl(6) 
END SUBROUTINE secthe

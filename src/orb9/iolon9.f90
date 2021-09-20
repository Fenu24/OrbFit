! ==================================================================    
!  input/output routines for longit8x; version 8                        
! ==================================================================    
!  outpro                                                               
!  print proper elements, frequencies                                   
! ==================================================================    
! in this version, proper elements are printed anyway;                  
! however, the quality code contains warnings and some are in a         
! separate file; also added quality code for mean elements              
SUBROUTINE outpro(inp,t,nam0,aa,pre,prsini,gsy,ssy,iqce,          &
     &    api1,ate1,g0sy,s0sy,iqcf,npas,ngper,apiv,ngnod,atev           &
     &    ,ipro,iang,iqcm,iqco,hm,irfl,irfl2)                           
  USE fund_const
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: inp ! input mode
! string for alphanumeric asteroid name                                 
  CHARACTER*9, INTENT(IN) :: nam0 
  DOUBLE PRECISION, INTENT(INOUT ) :: aa,pre,prsini ! proper elems.
  DOUBLE PRECISION, INTENT(INOUT) :: api1,ate1,apiv,atev,g0sy,s0sy,gsy,ssy ! proper angles, frequencies
  DOUBLE PRECISION, INTENT(IN) :: hm,t ! magnitude, epoch
  INTEGER, INTENT(IN) :: iqcm,iqco,irfl,irfl2 ! quality codes, flags
  INTEGER, INTENT(INOUT) :: iqce,iqcf
  INTEGER, INTENT(IN) :: npas,ngper,ngnod ! no. iter, no. revolutions
  INTEGER, INTENT(IN) :: ipro,iang  ! output units
! END INTERFACES
  INTEGER ir
  DOUBLE PRECISION api, ate
! safeties against format errors (absurd cases)                         
  if(pre.gt.0.9999d0)pre=0.9999d0 
  if(prsini.gt.0.9999d0)prsini=0.9999d0 
  if(gsy.gt.9999.d0)gsy=9999.d0 
  if(gsy.lt.-9999.d0)gsy=-9999.d0 
  if(ssy.gt.9999.d0)ssy=9999.d0 
  if(ssy.lt.-9999.d0)ssy=-9999.d0 
  if(g0sy.gt.9999.d0)g0sy=9999.d0 
  if(g0sy.lt.-9999.d0)g0sy=-9999.d0 
  if(s0sy.gt.9999.d0)s0sy=9999.d0 
  if(s0sy.lt.-9999.d0)s0sy=-9999.d0 
  if(inp.eq.1.or.inp.eq.4)then 
! ****************************************************                  
!  catalogue of proper elements                                         
     if(iqce.lt.0)iqce=0 
     if(iqcf.lt.0)iqcf=0 
     if(irfl.eq.0)then 
        if(iqce.gt.10)then 
           ir=90+(iqce-10)/3 
        else 
           ir=0 
        endif
     else 
        ir=irfl 
     endif
     api=api1*360/dpig 
     if(api.gt.360)api=api-360.d0 
     ate=ate1*360/dpig 
     if(ate.gt.360)ate=ate-360.d0 
     write (ipro,200) nam0,aa,pre,prsini,gsy,ssy,                &
     &      ir,iqcm,iqco,hm                                             
200  format (a9,f8.5,1x,f5.4,1x,f5.4,2f10.3,i4,i3,i3,1x,f5.2) 
     write(iang,400)nam0,api,ate                                 &
     &       ,g0sy,s0sy,npas,iqcf,iqce,irfl2                            
400  format(a9,2f9.3,2f10.3,4i4) 
  elseif(inp.eq.2.or.inp.eq.3)then 
! ****************************************************                  
!   output of numerical integration                                     
!   the proper angles are computed for current time, with number of     
!   revolutions accounted for; output in degrees                        
     call prngup(ate1,atev,ngnod) 
     ate1=ate1*360/dpig 
     ate1=ate1+ngnod*360 
!   for reduction to t=t0:   -(t-t0)*ssy/3600                           
     call prngup(api1,apiv,ngper) 
     api1=api1*360/dpig 
     api1=api1+ngper*360 
     if(irfl.eq.0)then 
        if(iqce.gt.10)then 
           ir=90+(iqce-10)/3 
        else 
           ir=0 
        endif
     else 
        ir=irfl 
     endif
!   for reduction to t=t0:   -(t-t0)*gsy/3600                           
     write (ipro,201) t,aa,pre,prsini,gsy,ssy,ir,iqcm 
201  format (f12.2,f12.8,2f9.6,2f12.4,i4,i3) 
     write (iang,202) t,api1,ate1,g0sy,s0sy,npas,iqcf,              &
     &        iqce,irfl2                                                
202  format (f12.2,2f15.5,2f12.4,4i4) 
  endif
END SUBROUTINE outpro
! ---------------------------------------------------------------       
! input of filtered output and conversion to keplerian                  
SUBROUTINE inpfil(iun,tt,em,eof) 
  USE fund_const
  IMPLICIT NONE 
  INTEGER iun 
  DOUBLE PRECISION em(6),tt 
  LOGICAL eof 
  DOUBLE PRECISION eh,ek,eq,ep,awt 
  read(iun,*,end=99)tt,em(1),eh,ek,ep,eq 
  eof=.false. 
  if(eh.eq.0.d0.and.ek.eq.0.d0)then 
     awt=0.d0 
  else 
     awt=atan2(eh,ek)*degrad 
  endif
  if(ep.eq.0.d0.and.eq.eq.0.d0)then 
     em(4)=0.d0 
  else 
     em(4)=atan2(ep,eq)*degrad 
  endif
  em(5)=awt-em(4) 
  em(2)=sqrt(eh**2+ek**2) 
  em(3)=atan(sqrt(ep**2+eq**2))*2.d0*degrad 
  RETURN 
99 eof=.true. 
END SUBROUTINE inpfil
! ==================================================================    
!   openfi                                                              
!   open and intest routine for longit7                                 
! ********************************************************************* 
SUBROUTINE openfi(inp,t0,filast,filfil,filpro,filpla,filbar       &
     &    ,iun2,iun3,iun8,iun9,iun10,iun11,iun19,iun21,           &
!iun18,
     &     iun7,iun17)     
  IMPLICIT NONE 
! input type, file names                                                
  INTEGER inp 
  CHARACTER*60 filpro,filast,filfil,filpla,filbar 
! output                                                                
! reference time, only for inp=3                                        
  DOUBLE PRECISION t0 
! unit numbers                                                          
  INTEGER iun2,iun3,iun8,iun9,iun10,iun11,iun19,iun21,iun18,iun7,iun17
! end interface                                                         
  INTEGER len,nast,ilce 
  CHARACTER*60 filang,filres,filmea,filnam,fildis,filadi,filrem,filprp  
  CHARACTER*3 coox,unitx,sysx 
  CHARACTER*6 refx 
  CHARACTER*4 nomast(20) 
  CHARACTER*60 comast 
  CHARACTER*100 colhea 
  INTEGER ibar 
! open input files                                                      
  IF(inp.eq.1.or.inp.eq.2)THEN 
! input either catalog or time series, but of osculating elements anyway
! open files for input of osculating elements                           
     call rmbl(filast,len) 
     IF(inp.eq.1)THEN
        CALL oporbf(filast(1:len), -1)
     ELSEIF(inp.eq.2)THEN 
        call filopn(iun8,filast(1:len),'old')                       
     ENDIF
! open input file for planetary data                                    
     call rmbl(filpla,len) 
     call filopn(iun7,filpla,'old') 
! open input file for barycenter of inner solar system                  
     if(inp.eq.2)then 
        ibar=0 
        write(*,*)' no barycentric correction in time series' 
     elseif(inp.eq.1)then 
        ibar=1 
        write(*,*)' baricenter from file ',filbar 
        call rmbl(filbar,len) 
        call filopn(iun17,filbar,'old') 
     endif
!  ELSEIF(inp.eq.3)THEN 
!  input directly from orbit9 (digitally filtered)                      
!     WRITE(*,*)' inp=3 NOT SUPPORTED'
!     STOP 
!     call rmbl(filfil,len) 
!     call filopn(iun18,filfil(1:len),'old') 
!     call reahea(iun18,nomast,comast,colhea,coox,sysx,refx,      &
!     &           unitx,nast,ilce,t0)                                    
!     if(nast.gt.1)then 
!        write(*,*)nast,' orbts at once not allowed; use conv8v' 
!        stop 
!     endif
!     if(coox.ne.'EQU'.or.sysx.ne.'HEL'.or.refx.ne.'INVL1B')then 
!        write(*,*)' wrong input coordinates/ref.system, use conv8v' 
!        stop 
!     endif
  ELSE 
     WRITE(*,*)' inp =',inp,' not allowed' 
     STOP 
  ENDIF
! open output files                                                     
  IF(inp.eq.1.or.inp.eq.2)THEN 
     ibar=22 ! arbitrary number to avoid undefined
! open files for output of mean elements                                
     call rmbl(filpro,len) 
     filmea=filpro(1:len)//'.mea' 
     call filopn(iun19,filmea(1:len+4),'unknown') 
     filrem=filpro(1:len)//'.rem' 
     call filopn(iun21,filrem(1:len+4),'unknown') 
  ENDIF
!   open files with names from prop.opt                                 
!   these always exist                                                  
  call rmbl(filpro,len) 
  filprp=filpro(1:len)//'.pro' 
  call filopn(iun2,filprp(1:len+4),'unknown') 
  filang=filpro(1:len)//'.ang' 
  call filopn(iun9,filang(1:len+4),'unknown') 
  filres=filpro(1:len)//'.res' 
  call filopn(iun10,filres(1:len+4),'unknown') 
  IF(inp.eq.1)THEN 
!   discarded cases in a separate file                                  
     fildis=filpro(1:len)//'.dis' 
     call filopn(iun3,fildis(1:len+4),'unknown') 
     filadi=filpro(1:len)//'.adi' 
     call filopn(iun11,filadi(1:len+4),'unknown') 
  ENDIF
! *********************************************                         
!   output headers for all input types                                  
  write(iun2,592) 
592 format(' PROPER ELEMENTS MILANI and KNEZEVIC') 
  write(iun9,592) 
  write(iun10,592) 
  if(inp.eq.1)then 
! original input from catalog of osculating elements                    
! no input to long periodic part because  discarded,
!  output in a separate file                                  
     write(iun3,592) 
     write(iun11,592) 
  endif
END SUBROUTINE openfi
! ==================================================================    
!  selout                                                               
!  conversions, output selection, output                                
SUBROUTINE selout(inp,t,nam0,                                     &
     &    gv,sv,aa,pre,prsini,iqce,iqcm,iqco,                           &
     &    hm,api,ate,g0,s0,ipas,iqcf,irfl,irfl2,                        &
     &        ngper,apiv,ngnod,atev,conv,iun2,iun3,iun9,iun11)          
  IMPLICIT NONE
  DOUBLE PRECISION gv,sv,aa,pre,prsini,hm,api,ate, g0,s0,apiv,atev,conv,t
  INTEGER iqce,iqcm,iqco,iun2,iun3,iun9,iun11,ipas,iqcf,irfl,irfl2, ngper,inp,ngnod
  LOGICAL discard 
  CHARACTER*9 nam0 
! END INTERFACE
  DOUBLE PRECISION gsy, ssy, g0sy, s0sy
  INTEGER ipro, iang
!  conversions to arcsec/yr                                             
  gsy=gv*conv 
  ssy=sv*conv 
  g0sy=g0*conv 
  s0sy=s0*conv 
  call select(inp,aa,pre,prsini,gsy,ssy,iqce,iqcf,iqcm,discard)
  IF(discard)THEN 
     ipro=iun3 
     iang=iun11 
  ELSE 
     ipro=iun2 
     iang=iun9 
  ENDIF
  call outpro(inp,t,nam0,aa,pre,prsini,gsy,ssy,iqce,                &
     &    api,ate,g0sy,s0sy,iqcf,ipas,ngper,apiv,ngnod,atev             &
     &    ,ipro,iang,iqcm,iqco,hm,irfl,irfl2)                           
END SUBROUTINE selout
! ==================================================================    
!  elems8                                                               
!  if it is meaningful, compute proper e and sin I                      
!  if blow up occurs (e>1 and/or sin I>1) write strong warning          
!    input: rho(4) proper Poincare' elements                            
!    output: ami,ani proper actions*2; api,ate proper angles            
!            pre,prsini proper keplerian el                             
!            idisc discarded elements code (see below)                  
SUBROUTINE elems8(rho,ani,ami,api,ate,pre,prsini,idisc) 
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT(IN) :: rho(4)
  DOUBLE PRECISION, INTENT(OUT) :: ani,ami,api,ate,pre,prsini
  INTEGER, INTENT(OUT) :: idisc 
! END INTERFACE
  DOUBLE PRECISION cani, cami
!   proper actions (doubled)                                            
  ani=rho(1)**2+rho(2)**2 
  ami=rho(3)**2+rho(4)**2 
! proper angles                                                         
  api=datan2(rho(2),rho(1)) 
  ate=datan2(rho(4),rho(3)) 
!  proper e                                                             
  if(ani.lt.2.d0)then 
     cani=ani-ani**2/4.d0 
     pre=sqrt(cani) 
     idisc=0 
  else 
!  e>1: idisc=1                                                         
     idisc=1 
  endif
!  proper sin I                                                         
  if(ami/sqrt(1-cani).lt.2.d0)then 
     cami=ami/sqrt(1-cani)-ami**2/4.d0/(1-cani) 
     prsini=sqrt(cami) 
  else 
!  sin I >1: idisc=2; sin I >1.and.e>1: idisc=3                         
     idisc=idisc+2 
  endif
END SUBROUTINE elems8

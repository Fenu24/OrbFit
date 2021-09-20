! =============MODULE runge_kutta_gauss=====
! CONTAINS
! public routines
!            rkimp
!            kintrp
!            legnum
!            rkstep
!            rkg    used by falsi
! private routines
!            fct    reduction to ord.1 for rkimp
!            fctcl    "            "   for rkg
!HEADERS, MUDULES used
! runge_kutta_gauss.o: \
!	../include/nvarx.h90 \
!	../suit/OUTPUT_CONTROL.mod \
!	../suit/PLANET_MASSES.mod \
!	force_model.o 
!
! to be done: rkimp and rkg are partially duplicates...
! 
MODULE runge_kutta_gauss
USE output_control
USE fund_const
IMPLICIT NONE

PRIVATE

! public routines

PUBLIC rkg, rkimp, kintrp, legnum, rkstep

! public data (options set by inipro)
! former comint.h
! controls for numerical integration: rkimp 
DOUBLE PRECISION eprk,eprk_c
INTEGER isrk,lit1,lit2,isrk_c,lit1_c
PUBLIC isrk,lit1,lit2,eprk,isrk_c,lit1_c,eprk_c
! for dimensioning in propin
INTEGER,PARAMETER:: ismax=20
INTEGER, PARAMETER :: itmax=50 ! also for close app.
PUBLIC ismax, itmax 

! common private data
! rkcoef.h
! ======================================
!  coefficients RKG for rkimp
      double precision ark(ismax,ismax),brk(ismax),crk(ismax)
!  coefficients interpolator for k array, for kintrp
      double precision a1(ismax,ismax)

CONTAINS

! ***************************************************************       
!  RKIMP ORBFIT                                                         
!                                                                       
!  scopo : propagatore che compie un passo con                          
!       il metodo di runge kutta implicito (formula di gauss)   
!       ad is passi. il numero di iterazioni non deve superare          
!       lit .  fct e' il secondo membro, che riduce al primo            
!       ordine il vero primo membro fct2                                
!  input:                                                               
!      t1 : tempo iniziale                                              
!      h : lunghezza del passo da eseguire                              
!      y1(nvar) : vettore di stato al tempo t1                          
!      ck(ismax,nvar) : valori iniziali per l'equazione implicita       
!      epsi : controllo di convergenza per l'equazione implicita        
!      ndim : no. variabili da usare per la norma da cfr. con epsi      
!  output :                                                             
!      t1: invariato                                                    
!      y3(nvar) : stato al tempo t1+h (non puo' avere lo stesso         
!               indirizzo di y1)                                        
!      ep(i) : controlli numerici;per i=1,lit converg.in                
!              gauss-seidel                                             
!      dery(nvar) : spazio di lavoro per il calcolo del secondo membro  
!      lf : flag che e' >0 se c'e' stata soddisfacente convergenza      
!              altrimenti segnala problemi                              
!  codice indici: i=1,nvar; j=1,is; id=1,ndim; it=indice di iterazione  
      SUBROUTINE rkimp(t1,h,y1,dery,ck,is,y3,lit,nvar,epsi,ep,lfleps,ndim)  
!INPUT                                                       
      double precision, intent(in) ::  t1,h !  tempo, passo 
      integer,intent(in):: nvar !  dimensioni variabili 
      double precision,intent(in):: y1(nvar) ! current state
      double precision,intent(in) ::  epsi !  controlli 
! OUTPUT
      double precision, intent(out) :: dery(nvar),y3(nvar) ! derivatives, new state
      double precision, intent(out) :: ep(itmax) ! controls on convergence
      double precision, intent(inout) :: ck(ismax,nvar) ! stack
! end interface
!  indici di iterazione, di sottopasso, di dimensione, flag             
      integer it,lit,j,i,is,jj,ndim,id,lfleps 
      INTEGER ips,imem ! controllo memoria
      double precision de,t(ismax) !  temporanei 
!****************                                                       
!   static memory not required                                          
!****************                                                       
! ===============================================                       
      do 8 j=1,is 
    8   t(j)=t1+h*crk(j) 
      do 7 it=1,itmax 
    7   ep(it)=0.d0 
!                                                                       
!  gauss-seidel per i ck                                                
!  inizio iterazioni per i ck                                           
      ips=0 
      it=1 
!  main loop                                                            
    1 do 11 j=1,is 
        do 12 i=1,nvar 
          de=0.d0 
          do 13 jj=1,is 
   13       de=de+ark(j,jj)*ck(jj,i) 
   12     y3(i)=de*h+y1(i) 
        imem=j 
        call fct(t(j),y3,dery,nvar,ips,imem) 
        do 14 i=1,ndim 
   14      ep(it)=ep(it)+dabs(dery(i)-ck(j,i)) 
        do 15 id=1,nvar 
   15      ck(j,id)=dery(id) 
   11 continue 
!  controllo se le iterazioni g-s sono finite                           
      ep(it)=ep(it)/is 
      lfleps=it 
      if(ep(it).gt.epsi)then 
         if(it.ge.lit)then 
!  troppe iterazioni in gauss-seidel                                    
!  il nuovo valore non viene calcolato                                  
            lfleps=-it 
            return 
          else 
            it=it+1 
            ips=-1 
            goto 1 
          endif 
      endif 
!                                                                       
!  calcolo nuovo punto y3                                               
      do 41 i=1,nvar 
        de=0.d0 
        do 41 j=1,is 
          de=de+brk(j)*ck(j,i) 
   41     y3(i)=y1(i)+h*de 
      return 
      END SUBROUTINE rkimp

! ===========================================================           
!   FCT                                                                 
!   right hand side routine                                             
!   reduces equations to order 1                                        
!   starting from accelerations computed by force                       
SUBROUTINE fct(t,y,dery,nvar,ips,imem)
  USE planet_masses 
  USE force_model
  USE force_sat
  USE force9d
  INTEGER, INTENT(IN) :: nvar
  DOUBLE PRECISION, INTENT(IN) :: y(nvar)
  DOUBLE PRECISION, INTENT(IN) :: t 
  INTEGER, INTENT(IN) :: ips,imem
  DOUBLE PRECISION, INTENT(OUT) :: dery(nvar) 
! END INTERFACE 
  DOUBLE PRECISION xxpla(6) 
  INTEGER nvar2,idc,i 
! close app. control options from force_model.mod
!****************                                                       
!   static memory not required                                          
!****************                                                       
  nvar2=nvar/2 
  IF(rhs.eq.1)THEN
     call force(y,y(nvar2+1),t,dery(nvar2+1),nvar2,idc,xxpla,ips,imem)
  ELSEIF(rhs.eq.2)THEN
     call forcesat(y,y(nvar2+1),t,dery(nvar2+1),nvar2,idc,xxpla,ips,imem)
  ELSEIF(rhs.eq.3)THEN
     CALL force9(y,y(nvar2+1),t,dery(nvar2+1),nvar2,idc,xxpla,ips,imem)
  ENDIF
  IF(rhs.eq.1)THEN 
     if(iclap.ne.0.and.idc.ne.0)then 
! to be improved with a real close approach subroutine like closapp/fals
!            write(*,*)'t =',t,' close approach to planet=',             &
!     &           ordnam(idc)                                            
!            write(iuncla,*)'t =',t,' close approach to planet=',        &
!     &           ordnam(idc)                                            
     endif
  ELSEIF(rhs.eq.3)THEN 
!         if(idc.ne.0)then 
!            write(*,*)'t =',t,' close approach code=',idc 
!            write(iuncla,*)'t =',t,' close approach code =',idc 
!         endif 
  ELSEIF(rhs.eq.2)THEN
! handle satellite case (ballistic orbit)
  ENDIF
  do i=1,nvar2 
     dery(i)=y(nvar2+i) 
  enddo
END SUBROUTINE fct
                                       
! ==========================================================            
! RKSTEP                                                                
! automated RKG stepsize change                                         
! ==========================================================            
SUBROUTINE rkstep(ep,npas,nrk,lf,h) 
  DOUBLE PRECISION ep(itmax),h 
  INTEGER npas,nrk,lf,l,i 
!  printout data on aborted step                                        
  l=iabs(lf) 
  write (*,*) 'Non-convergence in rk-gauss. See .pro file.' 
  write (ipirip,100)npas,nrk,lf,ep(l),eprk,h,isrk,(ep(i),i=1,l) 
100 format(' non convergence in rk at step ',i4,' nrk ',i3,' lf= ',i5 &
     &/' last control ',d14.5,' convergence required ',d12.3/           &
     &' stepsize ',d14.6,'  order 2*',i2/                               &
     &' controls ',5d12.3/(5d12.3/))                                    
! change stepsize and retry                                             
  h=0.8d0*h 
! with h changed, kintrp can not be used;                               
! also the catatst prepared for the mulktistep needs to be redone       
  nrk=0 
  lf=2 
  npas=0 
!
END SUBROUTINE rkstep

! ***************************************************************       
!  {\bf  kintrp} ORB8V                                                  
!                                                                       
!  scopo : predizione dei valori di ck per interpolazione da quelli     
!          del passo precedente (con polinomio di grado is-1)           
!          serve per avere condizioni iniziali vicine al punto unito    
!          nel procedimento iterativo per risolvere le equazioni        
!          implicite la cui soluzione e' la matrice ck.                 
!  input :                                                              
!       ck1(ismax,nvar) : ck al passo precedente                        
!       is : numero di stadi del metodo rk                              
!       nvar : numero di variabili nell'equazione di moto               
!  output:                                                              
!       ck(ismax,nvar) : valori interpolati                             
!  osservazione : non va usato nel passo iniziale e se rispetto al      
!                 passo precedente e' cambiato il passo, oppure is.     
!****************                                                       
!   static memory not required                                          
!****************                                                       
      SUBROUTINE kintrp(ck1,ck,is,nvar) 
!  dimensioni dipendenti da ismax (qui=ismax)                           
      integer nvar 
      double precision ck1(ismax,nvar),ck(ismax,nvar) 
!  numero di step intermedi                                             
      integer is 
!  indici di loop                                                       
      integer j,n,jj 
!  temporanei                                                           
      double precision de 
      do 1 j=1,is 
        do 1 n=1,nvar 
        de=0.d0 
        do 2 jj=1,is 
    2     de=de+a1(j,jj)*ck1(jj,n) 
    1   ck(j,n)=de 
      return 
      END SUBROUTINE kintrp                                          
! ***************************************************************       
!  LEGNUM                                                               
! ***************************************************************       
! reads  Runge--Kutta--gauss coefficients                               
! to be read in file ./lib/rk.coef                                      
!  is=required number of substeps; isfl=0 if found, otherwise           
!  uses closest available                                               
! ==============INTERFACE===========================                    
      SUBROUTINE legnum(is,isfl) 
! =========INPUT============                                            
      integer is 
! ========OUTPUT============                                            
      integer isfl 
! input unit, current is                                                
      integer iun,is1 
! loop indexes                                                          
      integer i,j,jj 
! skip trick                                                            
      character*1 cc(ismax) 
!****************                                                       
!   static memory not required                                          
!****************                                                       
! reads RKG coefficients rk                                             
      isfl=-1 
      call filopl(iun,'rk.coe') 
  198 read(iun,100,end=199)is1 
  100 format(6x,i4) 
!  control on is compatible wtih parameter ismax                        
      if(is1.gt.ismax.or.is1.le.0)goto 199 
      if(is1.eq.is)then 
         read(iun,101)(crk(j),j=1,is1) 
  101    format(7x,5d24.16) 
         read(iun,102)(brk(j),j=1,is1) 
  102    format(7x,5d24.16) 
         do 103 j=1,is1 
           read(iun,104)(ark(i,j),i=1,is1) 
  103    continue 
  104    format(3x,5d24.16) 
         do 105 j=1,is1 
           read(iun,106)(a1(i,j),i=1,is1) 
  106      format(4x,5d24.16) 
  105    continue 
         isfl=0 
         goto 199 
       else 
         read(iun,201)(cc(j),j=1,is1) 
  201    format(7x,5(23x,a1)) 
         read(iun,201)(cc(j),j=1,is1) 
         do 203 j=1,is1 
           read(iun,204)(cc(jj),jj=1,is1) 
  204    format(3x,5(23x,a1)) 
  203    continue 
         do 205 j=1,is1 
           read(iun,206)(cc(jj),jj=1,is1) 
  206    format(4x,5(23x,a1)) 
  205    continue 
         isfl=is1 
         goto 198 
      endif 
! end read                                                              
  199 call filclo(iun,' ') 
      return 
      END SUBROUTINE legnum  
! =================================================                     
! RKG                                                                   
!  Runge-Kutta-Gauss, used as a pure single step                        
! =================================================                     
      SUBROUTINE rkg(t1,xa,va,nv,h,xat,vat,xplat) 
! INPUT 
      DOUBLE PRECISION,INTENT(IN) ::  t1,h ! time, stepsize 
      INTEGER,INTENT(IN) :: nv ! state: position, velocities   
      DOUBLE PRECISION, INTENT(IN)::  xa(nv),va(nv) 
! OUTPUT                                                                
! state: position, velocities, approached planet position               
      DOUBLE PRECISION,INTENT(OUT) :: xat(nv),vat(nv),xplat(6) 
 !     END INTERFACE                                                         
      INCLUDE 'nvarx.h90' 
      DOUBLE PRECISION y1(nvarx),dery(nvarx),yat(nvarx),de 
      INTEGER nvar 
      DOUBLE PRECISION ep(itmax),ck(ismax,nvarx),t(ismax) 
      INTEGER i,j,it,jj 
      INTEGER imem,ips 
! control of derivatives: left as it is                                 
!     INTEGER ide,ideold                                                
!     COMMON/deriv/ide                                                  
!     ideold=ide                                                        
!     ide=0                                           
      nvar=2*nv ! state vector length
! state vector before the step                                          
      y1(1:nv)=xa(1:nv) 
      y1(nv+1:2*nv)=va(1:nv)    
      DO j=1,isrk_c 
        t(j)=t1+h*crk(j) ! set intermediate times 
      ENDDO 
      ck(1:isrk_c,1:nvar)=0.d0 ! array ck initialized to zero                      
      ep(1:itmax)=0.d0 ! controls on convergence inititalized to zero
!  gauss-seidel iterations for array ck                                 
      it=1 
      ips=0 
!  main loop on intermediate points                                     
    1 DO 11 j=1,isrk_c 
        DO 12 i=1,nvar 
          de=0.d0 
          DO  jj=1,isrk_c 
            de=de+ark(j,jj)*ck(jj,i) 
          ENDDO 
          yat(i)=de*h+y1(i) 
   12   continue 
! memory location not used by ra15v                                     
        imem=j+10 
        CALL fctcl(t(j),yat,dery,nvar,xplat,ips,imem) 
        DO i=1,nvar 
          ep(it)=ep(it)+dabs(dery(i)-ck(j,i)) 
        ENDDO 
        ck(j,1:nvar)=dery(1:nvar) 
   11 continue 
!  control on end of gauss-seidel iterations                            
      ep(it)=ep(it)/isrk_c 
      IF(it.gt.1)THEN
         IF(ep(it).gt.ep(it-1)*1.1d0)THEN 
            WRITE(iun_log,*)' rkg: stop at iteration ',it, ' before too late' 
            WRITE(iun_log,*)t1,(ep(jj),jj=it-1,it) 
            GOTO 77 
         ENDIF
      ENDIF 
!  new state vector yat                                                  
      DO i=1,nvar 
        de=0.d0 
        DO  j=1,isrk_c 
          de=de+brk(j)*ck(j,i) 
          yat(i)=y1(i)+h*de 
        ENDDO 
      ENDDO 
      IF(ep(it).gt.eprk)THEN 
         IF(it.ge.lit1_c)THEN 
!  too many gauss-seidel iterations                                     
            WRITE(iun_log,*)' rkg: non convergent after ',it,' iterations' 
            WRITE(iun_log,*)t1,ep(it) 
            GOTO 77 ! ADDED 11/6/2014
         ENDIF 
         it=it+1 
         ips=-1 
         GOTO 1 
      ENDIF 
!                                                                       
! right hand side at new state                                          
   77 CALL fctcl(t1+h,yat,dery,nvar,xplat,0,1) 
! copy pos. vel.                                                        
      DO i=1,nv 
        xat(i)=yat(i) 
        vat(i)=yat(i+nv) 
      ENDDO 
! control of derivatives reset: NO                                      
!     ide=ideold                                                        
      RETURN 
      END SUBROUTINE rkg                                          
! ==========================================================            
!   FCTCL                                                               
! reduction to first order, to use with RKG                             
! ==========================================================            
!   subroutine secondo membro                                           
!   equazioni ridotte all'ordine 1                                      
!   a partire dall'eq del secondo ordine                                
!   calcolata da force                                                  
    SUBROUTINE fctcl(t,y,dery,nvar,xxpla,ips,imem)
      USE force_model
      USE force_sat 
      USE force9d
      integer nvar 
      double precision y(nvar),dery(nvar) 
      double precision xxpla(6) 
      double precision t 
      INTEGER ips,imem 
! end interface                                                         
      integer nvar2,idc,i 
!****************                                                       
!   static memory not required                                          
!****************                                                       
      nvar2=nvar/2 
      IF(rhs.eq.1)THEN
         call force(y,y(nvar2+1),t,dery(nvar2+1),nvar2,idc,xxpla,ips,imem) 
      ELSEIF(rhs.eq.2)THEN
         call forcesat(y,y(nvar2+1),t,dery(nvar2+1),nvar2,idc,xxpla,ips,imem) 
      ELSEIF(rhs.eq.3)THEN
          CALL force9(y,y(nvar2+1),t,dery(nvar2+1),nvar2,idc,xxpla,ips,imem)
      ENDIF
      do  i=1,nvar2 
        dery(i)=y(nvar2+i) 
      enddo 
      return 
    END SUBROUTINE  fctcl   
                                  
END MODULE runge_kutta_gauss

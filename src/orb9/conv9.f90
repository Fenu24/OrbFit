! **********************************************************            
!                                                                       
!           c o n v 9                                                   
!                                                                       
! ***********************************************************           
PROGRAM conv9 
  USE fund_const
  USE massmod
  IMPLICIT NONE
  INTEGER, PARAMETER :: ndatx=100001 
!      parameter (nastx=100)                                            
!      parameter (nbodx=10)                                             
  INTEGER, PARAMETER :: nelex=6 
!      parameter (nangx=3)                                              
! integers                                                              
  INTEGER ifil,ifilp,iorfl(nastx) 
  INTEGER iia1,iia2 
! arrays for  data                                                      
  DOUBLE PRECISION x(nelex,nastx),y(nelex,nastx)
  INTEGER ng(nangx,nastx),ngp(nangx,nbox),ngold(nangx,nastx),ngs(nangx,nastx)
  DOUBLE PRECISION xp(nelex,nbox),yp(nelex,nbox) 
  DOUBLE PRECISION yold(nelex,nastx),supp(nangx,nastx), gamma(nastx)
  INTEGER iast(nastx) 
  DOUBLE PRECISION enne(nastx),bar(6) 
!  character arrays                                                     
  CHARACTER*3 coox,cooy,sysx,sysy,unitx,uniy 
  CHARACTER*6 refx,refy 
  CHARACTER*9 number(nastx) 
  CHARACTER*60 formx,formy,formx1,formy1 
! scalar variables
  INTEGER ilce, nast, imas, nbod, ifo, nia, iall, nangy, nacty, inter
  INTEGER isup, i, j, nerr, nn, nact1, nang1, jj, iia, iun, n, ndat
  DOUBLE PRECISION t, tp
! ***************************************************************       
!  input job specification                                              
  CALL jobinp(inter,nast,imas,nbod,ifo,number,iast,nia,coox,sysx,   &
     &      refx,unitx,ilce,cooy,sysy,refy,uniy,iall,ifil,iorfl)        
  write(*,*)'ilce=',ilce 
!c      write(*,*)' no. supplementary angles?'                          
!c      read(*,*)isup                                                   
  isup=0 
!  computation of formats (not really used in input, only number        
!  of variables is used)                                                
  call forma(isup,1,coox,formx,formx1,nact1,nang1) 
!  computation of the output format                                     
  call forma(isup,iall,cooy,formy,formy1,nacty,nangy) 
!  initialise revolution counters                                       
   ng=0
   ngold=0
   yold=0.d0
!  read-write loop                                                      
   nerr=0 
   DO 1 nn=1,ndatx 
      if(ifo.ge.1)then 
         read(1,*,end=3,err=4)t 
         DO 2 i=1,nast 
            if(ilce.lt.i)then 
!               read(1,formx,end=4,err=70)(x(n,i),n=1,nact1),           
               read(1,*,end=4,err=70)(x(n,i),n=1,nact1),               &
     &            (x(n+nact1,i),ng(n,i),n=1,nang1),                     &
     &            (supp(n,i),ngs(n,i),n=1,isup)                         
            elseif(ilce.ge.i)then 
!               read(1,formx1,end=4,err=70)(x(n,i),n=1,nact1),          
               read(1,*,end=4,err=70)(x(n,i),n=1,nact1),               &
     &            (x(n+nact1,i),ng(n,i),n=1,nang1),                     &
     &            (supp(n,i),ngs(n,i),n=1,isup)                         &
     &           ,gamma(i)                                              
            endif
            CYCLE 
70          nerr=nerr+1 
            if(nerr.lt.100) write(*,*) ' error, t=',t,' ast no ',i 
!             do 77 n=1,nact1                                           
! 77            x(n,i)=0.d0                                             
!             do 78 n=1,nang1                                           
!               x(n+nact1,i)=0.d0                                       
! 78            ng(n,i)=ngold(n,i)                                      
!             do 79 n=1,isup                                            
!              ngs(n,i)=ngs(n,i)                                        
! 79            supp(n,i)=0.d0                                          
            if(ilce.ge.i)then 
               gamma(i)=0 
            endif
2        ENDDO
      elseif(ifo.eq.0)then 
         if(ilce.le.0)then 
            read(1,*,end=3,err=4)t,(x(n,1),n=1,nact1),                &
     &         (x(n+nact1,1),ng(n,1),n=1,nang1),                        &
     &         (supp(n,1),ngs(n,1),n=1,isup)                            
         elseif(ilce.ge.1)then 
            read(1,*,end=4)t,(x(n,1),n=1,nact1),                      &
     &         (x(n+nact1,1),ng(n,1),n=1,nang1),                        &
     &         (supp(n,1),ngs(n,1),n=1,isup)                            &
     &         ,gamma(1)                                                
         endif
      endif
      DO j=ilce+1,nast 
         gamma(j)=0.d0
      ENDDO
!  conversion: angles in radiants                                       
      CALL radang(x,nast,coox,unitx,dpig) 
!  conversion of coordinates and coord. system                          
      if(imas.eq.1)then 
         CALL coord(x,sysx,nbod,coox,refx,y,sysy,cooy,refy,enne,bar)
      else 
         if(sysx.ne.sysy)then 
            read(2,*,end=33)tp 
            if(tp.ne.t)then 
               write(*,*)' planet time tp=',tp,' ast. time t=',t 
               stop 
            endif
            do 9 j=1,nbod-1 
               read(2,*,end=33)(xp(i,j),i=1,nact1),                    &
     &              (xp(i+nact1,j),ngp(i,j),i=1,nang1)                  
9           ENDDO
            CALL coord(xp,sysx,nbod,coox,refx,yp,sysy,cooy,refy,enne,bar)
         endif
         CALL cooast(x,sysx,nast,coox,refx,bar,nbod,y,sysy,cooy,refy,enne)
      endif
!  accounting of the number of revolutions                              
      call numrev(y,nast,x,ng,yold,ngold,coox,cooy) 
      DO 155 jj=1,nast 
         if(abs(ng(nang1,jj)-ngold(nang1,jj)).ge.10000000)then 
            write(*,*)' too many revolutions ',ng(nang1,jj),ngold(nang1,jj),jj
         endif
  155   ENDDO 
!  conversion of angles to selected units                               
        call unrad(y,nast,cooy,uniy,dpig) 
!  size control (to avoid format problems with wildly hyperbolic orbits)
        DO 55 iia=1,nia 
           DO n=1,nacty 
              if(y(n,iast(iia)).gt.99.9d00)then 
                 y(n,iast(iia))=99.9d0 
              elseif(y(n,iast(iia)).lt.-9.99d0)then 
                 y(n,iast(iia))=9.99d0 
              endif
           ENDDO
55      ENDDO
!  write                                                                
        if(iall.eq.0)then 
           DO 7 iia=1,nia 
              iun=10+iia 
              if(ilce.eq.0)then 
                 write(iun,formy)t,(y(n,iast(iia)),n=1,nacty),           &
     &          (y(n+nacty,iast(iia)),ng(n,iast(iia)),n=1,nangy),       &
     &          (supp(n,iast(iia)),ngs(n,iast(iia)),n=1,isup)           
              else 
                write(iun,formy1)t,(y(n,iast(iia)),n=1,nacty),          &
     &          (y(n+nacty,iast(iia)),ng(n,iast(iia)),n=1,nangy),       &
     &          (supp(n,iast(iia)),ngs(n,iast(iia)),n=1,isup)           &
     &          ,gamma(iast(iia))                                       
             endif
7         ENDDO
       elseif(iall.eq.1)then 
          ifilp=1 
          iia2=0 
   57     iun=10+ifilp 
          write(iun,*)t 
          iia1=iia2+1 
          iia2=iia2+iorfl(ifilp) 
          DO 58 iia=iia1,iia2 
             if(ilce.eq.0)then 
                write(iun,formy)(y(n,iast(iia)),n=1,nacty),             &
     &          (y(n+nacty,iast(iia)),ng(n,iast(iia)),n=1,nangy),       &
     &          (supp(n,iast(iia)),ngs(n,iast(iia)),n=1,isup)           
             else 
                write(iun,formy1)(y(n,iast(iia)),n=1,nacty),            &
     &          (y(n+nacty,iast(iia)),ng(n,iast(iia)),n=1,nangy),       &
     &          (supp(n,iast(iia)),ngs(n,iast(iia)),n=1,isup)           &
     &          ,gamma(iast(iia))                                       
             endif
   58     ENDDO 
          ifilp=ifilp+1 
          if(ifilp.gt.ifil-1) goto 1 
          goto 57 
       else 
          iun=11 
          write(iun,*)t 
          DO 8 iia=1,nia 
             if(ilce.eq.0)then 
                write(iun,formy)(y(n,iast(iia)),n=1,nacty),             &
     &          (y(n+nacty,iast(iia)),ng(n,iast(iia)),n=1,nangy),       &
     &          (supp(n,iast(iia)),ngs(n,iast(iia)),n=1,isup)           
             else 
                write(iun,formy1)(y(n,iast(iia)),n=1,nacty),            &
     &          (y(n+nacty,iast(iia)),ng(n,iast(iia)),n=1,nangy),       &
     &          (supp(n,iast(iia)),ngs(n,iast(iia)),n=1,isup)           &
     &          ,gamma(iast(iia))                                       
             endif
8         ENDDO
       endif
1   ENDDO
    write(*,*)' too many records, max was ',ndatx 
    stop 
3   ndat=nn-1 
    write(*,*)' regular end after ',ndat,' records' 
    stop 
4   write(*,*)' unreadable record no. ',nn 
    stop 
!  error in input of planets for system conversion                      
33  write(*,*)' insufficient planet data,t=',t 
  END PROGRAM conv9
! =========================================================             
!  {\bf forma}                                                          
!  computation of formats                                               
  SUBROUTINE forma(isup,iall,coo,form,form1,nact,nang)
    INTEGER, INTENT(IN) :: isup, iall
    INTEGER, INTENT(OUT) :: nact, nang
    CHARACTER*3, INTENT(IN) :: coo 
    CHARACTER*60, INTENT(OUT) :: form,form1 
    nact=numact(coo) 
    nang=numang(coo) 
    if(isup.eq.0.and.iall.eq.0)then 
!  no supplementary angle,one orbit per file                            
       if(coo.eq.'KEP')then 
         form='(f16.4,f13.9,f11.7,f11.6,2(f12.6,1x,i5),f12.6,1x,i9)' 
         form1='(f16.4,f13.9,f11.7,f11.6,2(f12.6,1x,i5),f12.6,1x,i9,e12.4)' 
      elseif(coo.eq.'EQU')then 
         form='(f16.4,f13.9,4f11.7,f12.7,1x,i9)' 
         form1='(f16.4,f13.9,4f11.7,f12.7,1x,i9,e12.4)' 
      elseif(coo.eq.'CAR')then 
         form='(f16.4,6f11.7)' 
         form1='(f16.4,6f11.7,e12.4)' 
      else 
         write(*,*)'  coordinates ',coo,' not supported' 
         stop 
      endif 
   elseif(isup.eq.1.and.iall.eq.0)then 
!  1 supplementary angle, 1 orbit per file                              
      if(coo.eq.'KEP')then 
         form=                                                          &
     & '(f16.4,f13.9,f11.7,f11.6,2(f12.6,1x,i5),2(f12.6,1x,i6))'        
         form1=                                                         &
     & '(f16.4,f13.9,f11.7,f11.6,2(f12.6,i6),2(f12.6,1x,i6),e12.4)'     
      elseif(coo.eq.'EQU')then 
         form='(f16.4,f13.9,4f11.7,f12.7,1x,i9,f12.7,i9)' 
         form1='(f16.4,f13.9,4f11.7,f12.7,1x,i9,f12.7,i9,e12.4)' 
      elseif(coo.eq.'CAR')then 
         form='(f16.4,6f11.7,f12.6,i9)' 
         form1='(f16.4,6f11.7,f12.6,i9,e12.4)' 
      else 
         write(*,*)'  coordinates ',coo,' not supported' 
         stop 
      endif 
   elseif(isup.eq.0.and.iall.ge.1)then 
!  no supplementary angle, many orbits per file                         
      if(coo.eq.'KEP')then 
         form='(f12.9,f11.7,f11.6,2(f12.6,1x,i5),f12.6,1x,i9)' 
         form1='(f12.9,f11.7,f11.6,2(f12.6,1x,i5),f12.6,1x,i9,1p,e12.4)' 
      elseif(coo.eq.'EQU')then 
         form='(f12.9,4f11.7,f11.8,1x,i9)' 
         form1='(f12.9,4f11.7,f11.8,1x,i9,1p,e12.4)' 
      elseif(coo.eq.'CAR')then 
         form='(6f11.7)' 
         form1='(6f11.7,e12.4)' 
      else 
         write(*,*)'  coordinates ',coo,' not supported' 
         stop 
      endif 
   elseif(isup.eq.1.and.iall.ge.1)then 
!  1 supplementary angle, many orbits per file                          
      if(coo.eq.'KEP')then 
         form=                                                          &
     & '(f12.9,f11.7,f11.6,2(f12.6,1x,i5),2(f12.6,1x,i6))'              
         form1=                                                         &
     & '(f12.9,f11.7,f11.6,2(f12.6,i6),2(f12.6,1x,i6),e12.4)'           
      elseif(coo.eq.'EQU')then 
         form='(f12.9,4f11.7,f11.8,1x,i9,f12.7,i9)' 
         form1='(f12.9,4f11.7,f11.8,1x,i9,f12.7,i9,e12.4)' 
      elseif(coo.eq.'CAR')then 
         form='(6f11.7,f12.6,i9)' 
         form1='(6f11.7,f12.6,i9,e12.4)' 
      else 
         write(*,*)'  coordinates ',coo,' not supported' 
         stop 
      endif 
   else 
      write(*,*)' forma, option not supported, isup=',isup,          &
     &      ' iall=',iall                                               
      stop 
   endif
 END SUBROUTINE forma

! *********************************************************             
!   {\bf outsub}  (vers. v.2.0)                                         
!  correction 14/5/93 to avoid r.o. errors in successive reading        
!  should be general purpose with output                                
!  always equinoctal, angles in radians                                 
!  if(fform)then formatted, otherwise unformatted for dump              
!                                                                       
! *********************************************************             
SUBROUTINE outsub(ip,ia,tf,eqf,angf,ngf,ncof,                     &
     &      norb,na,nvz,fform)                                          
  USE fund_const
      implicit double precision(a-h,o-z) 
      character*60 form,form1,form2 
      dimension eqf(ncof),angf(norb) 
      dimension ngf(norb) 
      logical fform 
      npla=norb-na 
!  planets                                                              
      if(fform)then 
         form='(f12.9,4f11.7,f11.8,1x,i9)' 
         form1='(f12.9,4f11.7,f11.8,1x,i9,1p,e12.4)' 
         form2='(f16.9,4f15.7,f11.8,1x,i9,1p,e12.4)' 
         write(ip,101)tf 
      else 
         form='(6d24.16,1x,i9)' 
         form1='(6d24.16,1x,i9,d12.4)' 
         form2=form1 
         write(ip,102)tf 
      endif 
  101 format(f16.4) 
  102 format(d24.16) 
      jq=0 
      DO j=1,npla 
        write(ip,form)(eqf(i+jq),i=1,5),angf(j),ngf(j) 
       jq=jq+5 
      ENDDO
!  asteroids and var.eq.                                                
      if(fform)then 
         write(ia,101)tf 
      else 
         write(ia,102)tf 
      endif 
      DO 2 j=1,na 
        if(j.gt.nvz)then 
           ga=0.d0 
        else 
           ga=eqf(ncof-nvz+j) 
        endif 
        if(eqf(jq+1).lt.-9.999d0.or.eqf(jq+1).gt.99.9d0.or.             &
     &     abs(eqf(jq+2)).gt.99.99d0)then                               
           write(ia,form2)(eqf(i+jq),i=1,5),angf(j+npla)                &
     &        ,ngf(j+npla),ga                                           
        else 
           write(ia,form1)(eqf(i+jq),i=1,5),angf(j+npla)                &
     &        ,ngf(j+npla),ga                                           
        endif 
       jq=jq+5 
    2 ENDDO
  END SUBROUTINE outsub

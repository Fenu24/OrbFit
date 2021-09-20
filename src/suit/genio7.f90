!=================================                                      
!   I/O library; longit7 version (Mar. 1995)                            
!   modified for orbit9 (October 1998)                                  
! ==========================                                            
!   skip                                                                
! to skip nlin lines in the input file (unit iun)                       
      subroutine skip(iun,nlin) 
      character*1 cc 
      do 1 n=1,nlin 
        read(iun,100)cc 
  100   format(a1) 
    1 continue 
      return 
      END                                           
! ================================                                      
!   named constants:                                                    
! ================================                                      
!   write named floating                                                
      subroutine wriflo(iun,name,c,lcom) 
      double precision c 
      character*(*) name,lcom 
      ln=lench(name) 
      lc=lench(lcom) 
      write(iun,100)name(1:ln),c,lcom(1:lc) 
  100 format(a,'= ',1p,d24.16,' ; ',a) 
      return 
      END                                           
! ================================                                      
!   write named integer                                                 
      subroutine wriint(iun,name,n,lcom) 
      character*(*) name,lcom 
      ln=lench(name) 
      lc=lench(lcom) 
      write(iun,100)name(1:ln),n,lcom(1:lc) 
  100 format(a,'= ',i8,' ; ',a) 
      return 
      END                                           
! ================================                                      
!  {\bf reaint}                                                         
!   read named integer                                                  
      subroutine reaint(iun,name,n) 
      character*(*) name 
      character*80 line,lnam,lval,lcom 
      call rmbl(name,nnam) 
      read(iun,100)line 
  100 format(a80) 
      call splili(line,ll,lnam,inam,lval,ival,lcom,icom) 
!   check variable identity                                             
      if(name(1:nnam).ne.lnam(1:inam))then 
         write(*,*)' input variable does not match ',name,' unit ',iun 
         write(*,*)' lnam=',lnam,' inam=',inam 
         write(*,*)' lval=',lval,'ival=',ival 
         write(*,*)' lcom=',lcom,' icom=',icom 
         stop 
      else 
         read(lval,*)n 
         return 
      endif 
      END                                           
! ================================                                      
! {\bf reaflo}                                                          
!   read named floating                                                 
      subroutine reaflo(iun,name,c) 
      double precision c 
      character*(*) name 
      character*80 line,lnam,lval,lcom 
      call rmbl(name,nnam) 
      read(iun,100)line 
  100 format(a80) 
      call splili(line,ll,lnam,inam,lval,ival,lcom,icom) 
!   check variable identity                                             
      if(name(1:nnam).ne.lnam(1:inam))then 
         write(*,*)' input variable does not match ',name,' unit ',iun 
         write(*,*)' lnam=',lnam,' inam=',inam 
         write(*,*)' lval=',lval,'ival=',ival 
         write(*,*)' lcom=',lcom,' icom=',icom 
         stop 
      else 
         read(lval,*)c 
         return 
      endif 
      END                                           
! ================================                                      
!  {\bf reaflc}                                                         
!   read named floating, with comment (max 60 characters)               
      subroutine reaflc(iun,name,c,commen) 
      double precision c 
      character*(*) name 
      character*80 line,lnam,lval,lcom 
      character*60 commen 
      call rmbl(name,nnam) 
      read(iun,100)line 
  100 format(a80) 
      call splili(line,ll,lnam,inam,lval,ival,lcom,icom) 
!   check variable identity                                             
      if(name(1:nnam).ne.lnam(1:inam))then 
         write(*,*)' input variable does not match ',name,' unit ',iun 
         write(*,*)' lnam=',lnam,' inam=',inam 
         write(*,*)' lval=',lval,'ival=',ival 
         write(*,*)' lcom=',lcom,' icom=',icom 
         stop 
      else 
         read(lval,*)c 
         commen=' ' 
         commen=lcom(1:icom) 
         return 
      endif 
      END                                           
! ================================                                      
!   read named string                                                   
      subroutine reastr(iun,name,str) 
      character*(*) name,str 
      character*80 line,lnam,lval,lcom,str1 
      call rmbl(name,nnam) 
      read(iun,100)line 
  100 format(a) 
      call splili(line,ll,lnam,inam,lval,ival,lcom,icom) 
!   check variable identity                                             
      if(name(1:nnam).ne.lnam(1:inam))then 
         write(*,*)' input variable does not match ',name,' unit ',iun 
         write(*,*)' lnam=',lnam,' inam=',inam 
         write(*,*)' lval=',lval,'ival=',ival 
         write(*,*)' lcom=',lcom,' icom=',icom 
         stop 
      else 
         if(ival.gt.2)then 
            str1=lval(2:ival-1) 
!           read(lval(1:ival),*)str1                                    
         else 
!  control against empty strings; they need at least one blank          
            str1=' ' 
         endif 
         l1=lench(str1) 
         l2=len(str) 
         if(l2.lt.l1)then 
             write(*,*)' Stringa letta troppo lunga' 
             write(*,*)' L(letta)=',l1,' L(max)=',l2 
         end if 
         str=str1(1:l1) 
         return 
      endif 
      END                                           
! ================================                                      
!   read named string, with comment                                     
      subroutine reastc(iun,name,str,lcom) 
      character*(*) name,str 
      character*80 line,lnam,lval,lcom,str1 
      call rmbl(name,nnam) 
      read(iun,100)line 
  100 format(a) 
      call splili(line,ll,lnam,inam,lval,ival,lcom,icom) 
!   check variable identity                                             
      if(name(1:nnam).ne.lnam(1:inam))then 
         write(*,*)' input variable does not match ',name,' unit ',iun 
         write(*,*)' lnam=',lnam,' inam=',inam 
         write(*,*)' lval=',lval,'ival=',ival 
         write(*,*)' lcom=',lcom,' icom=',icom 
         stop 
      else 
         if(ival.gt.2)then 
            read(lval(1:ival),*)str1 
         else 
!  control against empty strings; they need at least one blank          
            str1=' ' 
         endif 
         l1=lench(str1) 
         l2=len(str) 
         if(l2.lt.l1)then 
             write(*,*)' Stringa letta troppo lunga' 
             write(*,*)' L(letta)=',l1,' L(max)=',l2 
         end if 
         if(l1.ge.1)then 
            str=str1(1:l1) 
         else 
            str=' ' 
         endif 
         return 
      endif 
      END                                           
! ================================                                      
!   write named string; the strings must be dimensioned and of length 60
      subroutine wristr(iun,name,str,lcom) 
      character*60 name,str,lcom 
      ln=lench(name) 
      ls=lench(str) 
      lc=lench(lcom) 
      if(lc.eq.0)then 
         write(iun,101)name(1:ln),str(1:ls) 
  101    format(a,'= ''',a,''' ; ') 
      else 
         write(iun,100)name(1:ln),str(1:ls),lcom(1:lc) 
  100     format(a,'= ''',a,''' ; ',a) 
      endif 
      return 
      END                                           
! ================================                                      
!   read named logical                                                  
      subroutine realog(iun,name,flag) 
      character*(*) name 
      character*80 line,lnam,lval,lcom 
      logical flag 
      call rmbl(name,nnam) 
      read(iun,100)line 
  100 format(a) 
      call splili(line,ll,lnam,inam,lval,ival,lcom,icom) 
!   check variable identity                                             
      if(name(1:nnam).ne.lnam(1:inam))then 
         write(*,*)' input variable does not match ',name,' unit ',iun 
         write(*,*)' lnam=',lnam,' inam=',inam 
         write(*,*)' lval=',lval,'ival=',ival 
         write(*,*)' lcom=',lcom,' icom=',icom 
         stop 
      else 
         read(lval,*)flag 
         return 
      endif 
      END                                           
!=================================                                      
!  {\bf wrilog}                                                         
!   write named logical                                                 
      subroutine wrilog(iun,name,flag,lcom) 
      character*8 name*(*),lcom*(*) 
      logical flag 
      ln=lench(name) 
      lc=lench(lcom) 
      write(iun,100)name(1:ln),flag,lcom(1:lc) 
  100 format(a,'= ',l1,' ; ',a) 
      return 
      END                                           
!  ***************************************************************      
!                                                                       
!                         r m b l                                       
!                                                                       
!         remove blanks in a string                                     
!                                                                       
!  ***************************************************************      
!                                                                       
! input:    c         -  string                                         
!                                                                       
! output:   c         -  string without blanks                          
!           len       -  length of the compressed string                
!                                                                       
      subroutine rmbl(c,len) 
      character c*(*) 
    5 len=lench(c) 
      ip1=index(c,' ') 
      if(ip1.eq.0.or.ip1.gt.len)return 
      do 1 j=ip1,len-1 
    1 c(j:j)=c(j+1:j+1) 
      c(len:len)=' ' 
      goto 5 
      END                                           
!                                                                       
!  ***************************************************************      
!                                                                       
!                         c m p s t r                                   
!                                                                       
!         compression of multiple blanks in a string                    
!                                                                       
!  ***************************************************************      
!                                                                       
!                                                                       
! input:    c         -  string                                         
!                                                                       
! output:   c         -  compressed string                              
!           len       -  length of the compressed string                
!                                                                       
      subroutine cmpstr(c,len) 
      character c*(*) 
    5 len=lench(c) 
      ip1=index(c,'  ') 
      if(ip1.eq.0.or.ip1.gt.len)return 
      do 1 j=ip1,len 
      if(c(j:j).ne.' ')goto 2 
    1 continue 
    2 ip2=j 
      idif=ip2-ip1-1 
      do 3 j=ip1+1,len-idif 
      j1=j+idif 
    3 c(j:j)=c(j1:j1) 
      do 4 j=len-idif+1,len 
    4 c(j:j)=' ' 
      goto 5 
      END                                           
! =============================================                         
!   parsing of an input line of the form: "name=value;comment"          
      subroutine splili(line,ll,lnam,inam,lval,ival,lcom,icom) 
      character*80 line,lnam,lval,lcom 
      ll=lench(line) 
      call cmpstr(line,ll) 
      isc=index(line,';') 
      if(isc.eq.0)then 
          write(*,*)' missing ; in line ' 
          write(*,*)line 
          stop 
      endif 
      ieq=index(line,'=') 
      if(ieq.eq.0)then 
          write(*,*)' missing = in line ' 
          write(*,*)line 
          stop 
      endif 
      inam=ieq-1 
      lnam=line(1:inam) 
      call rmbl(lnam,inam) 
      ival=isc-1 
      lval=line(ieq+1:ival) 
      ival=ival-ieq 
      call rmbl(lval,ival) 
      icom=ll-isc 
      if(icom.gt.0)then 
         lcom=line(isc+1:ll) 
      else 
         lcom=' ' 
         icom=1 
      endif 
      return 
      END                                           

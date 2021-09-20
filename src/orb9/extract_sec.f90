program extract_sec
  implicit none 
  character*9 name,name1 
  character*100 a1,b1,c1 
  double precision a(8),b(8),c(7),d(7) 
  integer i,l,k,k1,m,m1 
                                                                        
! open files with proper elements and standard deviations 
! for INPUT              
  open(6,file='numb.syn',status='old') ! new numbered, all together (de>1)
  open(10,file='numb.sig',status='old') ! with header
  open(16,file='mult.syn',status='old') ! nuw multiopposition, all together
  open(20,file='mult.sig',status='old') ! without header
  open(7,file='secres0.syn',status='old')
  open(27,file='secres0.sig',status='old')
! for OUTPUT for create and users
  open(8,file='numb_withoutres.syn',status='unknown') ! non resonant 
  open(12,file='numb_withoutres.sig',status='unknown') 
  open(18,file='mult_withoutres.syn',status='unknown') 
  open(22,file='mult_withoutres.sig',status='unknown') 
  open(9,file='numb_secres.syn',status='unknown') ! resonant only de<0.1
  open(13,file='numb_secres.sig',status='unknown')
  open(19,file='mult_secres.syn',status='unknown') 
  open(23,file='mult_secres.sig',status='unknown')
! for OUTPUT for visualizer and matlab
  open(14,file='numb_res.syn',status='unknown') ! all together, de>0.7
  open(15,file='numb_res.sig',status='unknown')
  open(24,file='mult_res.syn',status='unknown')
  open(25,file='mult_res.sig',status='unknown')
!----------------------------------------------------------
! handle numbered only
!----------------------------------------------------------
! MODIFICATION 16/8/2017: 
! copy the header (2 lines) of numb.syn and numb.sig
  DO i=1,2 
     read(6,400)a1 
     write(8,400)a1 
     read(7,400)c1
     write(9,400)c1 
     write(14,400)a1
     read(10,400)b1 
     write(12,400)b1
     read(27,400)c1
     write(13,400)c1
     write(15,400)b1
  ENDDO
400 format(a100) 
! beginning of read loop
2 read(6,100,end=5) name,b,k1 ! read in numb.syn
  call rmsp(name,l) 
5 read(10,101,end=6) name1,d,m1 ! read in numb.sig
  call rmsp(name1,l) 
! decide the output                                                     
  if(name.eq.name1)then   
! find the secular resonant elements (e>1)
     if (b(3).ge.1.d0)then
        b(3)=b(3)-1.d0
        write(9,100)name,b,k1  ! write in numb_secres.syn
        write(13,101)name,d,m1 ! write in numb_secres.sig
        b(3)=b(3)+0.7d0
        write(14,100)name,b,k1 ! write in numb_res.syn
        write(15,101)name,d,m1 ! write in numb_res.sig
     else
        write(8,100)name,b,k1  ! write in numb_withoutres.syn
        write(12,101)name,d,m1 ! write in numb_withoutres.sig
        write(14,100)name,b,k1 ! write in numb_res.syn
        write(15,101)name,d,m1 ! write in numb_res.sig
     end if
100 format(a9,f6.2,f12.7,f11.7,f10.7,f12.6,f13.6,f13.6,f8.2,i4) 
101 format(a9,f12.8,f11.7,f10.7,f10.6,f11.6,f10.6,f13.4,i4) 
  else
     WRITE(*,*) 'misalignememt'
     WRITE(*,*) name,' not equal to ',name1
     STOP 
  endif
  goto 2
! normal file termination
6 close (6) 
  close (8) 
  close (9) 
  close (10)
  close (11)
  close (12) 
  close (13)
  close (14)
  close (15)
!----------------------------------------------------------
! handle multiopposition only
!----------------------------------------------------------
! no header; but fictitious header is needed for mult_res.syn, because
! it si assumed by prop_class
   write(24,'(A1/A1)')'%','%'
! beginning of read loop
12 read(16,100,end=15) name,b,k1 ! read in mult.syn
  call rmsp(name,l) 
15 read(20,101,end=16) name1,d,m1 ! read in mult.sig
  call rmsp(name1,l) 
! decide the output                                                     
  if(name.eq.name1)then   
! find the secular resonant elements (e>1)
     if (b(3).ge.1.d0)then
        b(3)=b(3)-1.d0
        write(19,100)name,b,k1  ! write in mult_secres.syn
        write(23,101)name,d,m1 ! write in mult_secres.sig
        b(3)=b(3)+0.7d0
        write(24,100)name,b,k1 ! write in mult_res.syn
        write(25,101)name,d,m1 ! write in mult_res.sig
     else
        write(18,100)name,b,k1  ! write in mult_withoutres.syn
        write(22,101)name,d,m1 ! write in mult_withoutres.sig
        write(24,100)name,b,k1 ! write in mult_res.syn
        write(25,101)name,d,m1 ! write in mult_res.sig
     end if
  else
     WRITE(*,*) 'misalignememt'
     WRITE(*,*) name,' not equal to ',name1
     STOP 
  endif
  goto 12
! normal file termination
16 close (16) 
  close (18) 
  close (19) 
  close (20)
  close (22) 
  close (23)
  close (24)
  close (25)
END program extract_sec

! program to merge outputs from synthetic and resonant proper elements  
! files and to                                                          
! replace e-1 into the resulting file, also to purge 2.5<a<2.7|sinI<0.3 
! from multi props                                                      
        program splitsec 
        implicit none 
        character*9 name,name1 
        character*100 a1,b1 
        double precision a(8),b(8),c(7),d(7) 
        integer i,l,k,k1,m,m1 
                                                                        
! open files with proper elements and standard deviations               
        open(6,file='numb.syn',status='old') 
        open(10,file='numb.sig',status='old') 
        open(7,file='secres0.syn',status='old') 
        open(11,file='secres0.sig',status='old') 
        open(8,file='numb_withoutres.syn',status='unknown') 
        open(12,file='numb_withoutres.sig',status='unknown') 
        open(9,file='secres.syn',status='unknown') 
        open(13,file='secres.sig',status='unknown') 
! copy the header of numb.syn and numb.sig                              
        do i=1,2 
           read(6,400)a1 
           write(8,400)a1 
           read(10,400)b1 
           write(12,400)b1 
        enddo 
  400   format(a100) 
                                                                        
! beginning of the read loop                                            
    1   read(7,100,end=4) name,a,k ! numb to be purged
        call rmsp(name,l) 
    4   read(11,101,end=2) name,c,m  ! 
        call rmsp(name,l) 
    2   read(6,100,end=5) name1,b,k1 
        call rmsp(name1,l) 
    5   read(10,101,end=6) name1,d,m1 
        call rmsp(name1,l) 
! decide the output                                                     
        if(name.eq.name1)then 
          if(a(3).ge.1.d0)then ! Delta e in place of e
              a(3)=a(3)-1.d0
              write(9,100)name,a,k ! Write in secres.syn/.sig
              write(13,101)name,c,m 
          else  ! non-resonant
              write(8,100)name1,b,k1  ! Write in numb_withoutres.syn/.sig
              write(12,101)name1,d,m1 
          endif 
          goto 1 ! get another from numb.syn/.sig
        else 
           write(8,100)name1,b,k1 ! Write in numb_withoutres.syn/.sig
           write(12,101)name1,d,m1 
           goto 2 ! get another from secres0.syn/.sig
        endif 
                                                                        
  100   format(a9,f6.2,f12.7,f11.7,f10.7,f12.6,f13.6,f13.6,f8.2,i4) 
  101   format(a9,f12.8,f11.7,f10.7,f10.6,f11.6,f10.6,f13.4,i4) 
                                                                        
    6   close (6) 
        close (7) 
        close (8) 
        close (9) 
        close (10) 
                                                                        
      END program splitsec

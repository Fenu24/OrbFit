! program to merge outputs from synthetic and resonant proper elements f
! replace e-1 into the resulting file, also to purge 2.5<a<2.7|sinI<0.3 
        program mergsec 
        implicit none 
        character*9 name,name1 
        character*100 a1 
        double precision a(9),b(9),c(8),d(8),x(9),y(8), g5, g6
        integer i,l 
! open files with proper elements                                       
        open(6,file='numb.syn',status='old') 
        open(7,file='secres.syn',status='old') 
        open(8,file='numb_res.syn',status='unknown') 
! copy the header of numb.syn                                           
        do i=1,2 
           read(6,400)a1 
           write(8,400)a1 
        enddo 
  400   format(a100) 
! beginning of the read loop                                            
    1   read(7,100,end=2) name,a 
        call rmsp(name,l) 
    2   read(6,100,end=3) name1,b 
        call rmsp(name1,l) 
!	write(*,*) name, name1                                                
!	pause                                                                 
        if(name.eq.name1)then 
           if(a(3).ge.1.d0)then 
              a(3)=a(3)-0.3d0 
           endif 
           write(8,100)name,a 
           goto 1 
        else 
           write(8,100)name1,b 
           goto 2 
        endif 
    3   continue 
  100   format(a9,f6.2,f12.7,f11.7,f10.7,f12.6,f13.6,f13.6,f8.2,i4) 
        close (6) 
        close (7) 
        close (8) 
                                                                        
! open files with standard deviations of proper elements                
        open(6,file='numb.sig',status='old') 
        open(7,file='secres.sig',status='old') 
        open(8,file='numb_res.sig',status='unknown') 
! copy the header of numb.sig                                           
        do i=1,2 
           read(6,400)a1 
           write(8,400)a1 
        enddo 
! beginning of the read loop                                            
    4   read(7,101,end=5) name,c 
        call rmsp(name,l) 
    5   read(6,101,end=6) name1,d 
        call rmsp(name1,l) 
!	write(*,*) name, name1                                                
!	pause                                                                 
        if(name.eq.name1)then 
           write(8,101)name,c 
           goto 4 
        else 
           write(8,101)name1,d 
           goto 5 
        endif 
    6   continue 
  101   format(a9,f12.8,f11.7,f10.7,f10.6,f11.6,f10.6,f13.4,i4) 
        close (6) 
        close (8) 
                                                                        
! open multiopposition files with proper elements                       
        open(6,file='mult.syn',status='old') 
        open(8,file='mult_purged.syn',status='unknown') 
        open(7,file='mult.sig',status='old') 
        open(9,file='mult_purged.sig',status='unknown') 
! copy the header of mult.syn/sig                                       
        do i=1,2 
           read(6,400)a1 
           write(8,400)a1 
           read(7,400)a1 
           write(9,400)a1 
        enddo 
        g5=4.2574d0 
        g6=28.2456d0 
    7   read(6,100,end=8) name,x 
        call rmsp(name,l)
    8   read(7,101,end=9) name1,y 
        call rmsp(name1,l)
        IF(name.ne.name1)THEN
           WRITE(*,*)' syn= ',name,' sig= ',name1
           GOTO 7
        ENDIF
        if(x(2).gt.2.5d0.and.x(2).lt.2.7d0.and.x(4).lt.0.3d0)then 
            if(abs(x(6)+g5-2*g6).lt.0.4d0)then 
              goto 7 
            else 
              write(8,100)name,x ! write in mult_purged.syn/.sig
              write(9,101)name1,y 
            endif 
         else 
            write(8,100)name,x ! write in mult_purged.syn/.sig
            write(9,101)name1,y 
        endif 
        goto 7 
    9   close (6) 
        close (7) 
        close (8) 
        close (9) 
                                                                        
      END program mergsec

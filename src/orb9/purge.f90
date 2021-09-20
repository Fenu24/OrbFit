! program to clean the multiopposition proper element files             
! from the objects numbered meanwhile or anyhow missing in the          
! current multiopposition object list                                   
        program purge 
        implicit none 
! initializations                                                       
        character*1 comcha,a1 
        character*10 name 
        character*11 name1(50000),rec 
        character*87 record 
!        character*109 record1                                          
        integer i,ipc,lr,imax 
! ==========================================================            
! files opening                                                         
        open(1,file='multall.pro', status='OLD') 
!	open(1,file='multall.sig', status='OLD')                              
!	open(1,file='multall.del', status='OLD')                              
        open(2,file='ufitobs.cat', status='OLD') 
        open(3,file='mult.pro', status='UNKNOWN') 
!	open(3,file='mult.sig', status='UNKNOWN')                             
!	open(3,file='mult.del', status='UNKNOWN')                             
! skip the headers; write the header in the output file                 
        do i=1,2 
           read(1,100)name,record 
           write(3,100)name,record 
        enddo 
        do i=1,6 
           read(2,101)a1,name1(i) 
        enddo 
! ========================================================              
! reading in osculating elements and                                    
! parse name to get rid of the trailing "'"                             
        comcha="'" 
        imax=1 
    5   read(2,101,END=10,ERR=999)a1,name1(imax) 
        rec=name1(imax) 
! compute length excluding comments                                     
        lr=len(rec) 
        ipc=INDEX(rec(1:lr),comcha) 
        IF(ipc.EQ.0) THEN 
          name1(imax)=rec(1:lr) 
        ELSE 
          name1(imax)=rec(1:ipc-1) 
          lr=len(name1(imax)) 
          IF(lr.LT.1) GOTO 998 
        END IF 
          imax=imax+1 
        goto 5 
! ========================================================              
! reading in proper elements                                            
   10   read (1,100,END=999,ERR=998)name,record 
        rec=name 
! compute length excluding comments                                     
        lr=len(rec) 
        ipc=INDEX(rec(1:lr),' ') 
        IF(ipc.EQ.0) THEN 
            name=rec(1:lr) 
        ELSE 
            name=rec(1:ipc-1) 
            lr=len(name) 
            IF(lr.LT.1) GOTO 998 
        END IF 
        write(*,*)name 
! ==========================================================            
! start comparison                                                      
        do i=1,imax 
        if(name.eq.name1(i))then 
           write(3,100)name,record 
           goto 10 
        endif 
        enddo 
        goto 10 
! formats                                                               
  100   format(a11,a81) 
  101   format(a1,a11) 
  998   write(*,*) 'error' 
999     CONTINUE
      END program purge

      subroutine selcon(nams,nns,vals) 
                                                                        
!  Input                                                                
                                                                        
!           nams :  a list of names (character*6)                       
                                                                        
!           nns  :  number of names in the list                         
                                                                        
!  Output                                                               
                                                                        
!           vals :  values corresponding to the input names             
                                                                        
                                                                        
                                                                        
      double precision vlc(400),sss(3),vals(1) 
                                                                        
      character*6 nmc(400),nams(1) 
                                                                        
      call const(nmc,vlc,sss,nnc) 
                                                                        
      do 2 i=1,nns 
                                                                        
      do j=1,nnc 
      if(nams(i) .eq. nmc(j)) go to 1 
      enddo 
                                                                        
      vals(i)=0.0 
      go to 2 
                                                                        
    1 vals(i)=vlc(j) 
                                                                        
    2 continue 
                                                                        
      return 
      END                                           

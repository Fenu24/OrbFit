! **************************************************************        
!  {\bf plotcr}                                                         
! plotting routine using GNUPLOT                                        
!                                                                       
      SUBROUTINE plotcr (xd, yd, n, xlab, ylab, title, idev, istyle)
      USE var_precision 
      USE char_string
      REAL x (n), y (n) 
      REAL (KIND=r_kind), INTENT(in) ::  xd (n), yd (n) 
      CHARACTER(80) ylab 
      CHARACTER(80) xlab, title 
      CHARACTER(82) label 
      CHARACTER(10) style, tempfi 
! write temporary file with data                                        
      tempfi = 'giffv.tmp' 
      DO 2 i = 1, n 
         x (i) = xd (i) 
         y (i) = yd (i)
    2 END DO 
      OPEN (22, file = tempfi, status = 'unknown') 
      DO 1 i = 1, n 
         WRITE (22, * ) x (i), y (i)
    1 END DO 
      CLOSE (22) 
! plot style                                                            
      IF (istyle.eq.1) then 
         style = 'lines' 
      ELSEIF (istyle.eq.0) then 
         style = 'dots' 
      ELSEIF (istyle.lt.0) then 
         style = 'points' 
      ELSE 
         WRITE ( * , * ) ' style code ', istyle, ' not known' 
         RETURN 
      ENDIF 
!  preparation of command file                                          
      OPEN (22, file = 'giffv.gnu') 
!  xlable,ylabel, title                                                 
      WRITE (22, * ) 'set nokey' 
      nc = lench (xlab) 
      label = '"'//xlab (1:nc) //'"' 
      WRITE (22, * ) 'set xlabel ', label (1:nc + 2) 
      nc = lench (ylab) 
      label = '"'//ylab (1:nc) //'"' 
      WRITE (22, 122) label (1:nc + 2) 
  122 FORMAT('set ylabel ',a82) 
      nc = lench (title) 
      label = '"'//title (1:nc) //'"' 
      WRITE (22, * ) 'set title ', label (1:nc + 2) 
!  set graphics terminal type and output device                         
      IF (idev.eq.3) then 
         WRITE (22, * ) 'set terminal X11' 
         WRITE (22, * ) 'set output' 
      ELSEIF (idev.eq.2) then 
         WRITE (22, * ) 'set terminal vgalib' 
         WRITE (22, * ) 'set output' 
      ELSEIF (idev.eq.4) then 
         WRITE (22, * ) 'set terminal tek40xx' 
         WRITE (22, * ) 'set output' 
      ELSEIF (idev.eq. - 5) then 
         WRITE (22, * ) 'set terminal postscript monochrome' 
         WRITE (22, * ) 'set output "giffvbw.eps"' 
      ELSEIF (idev.eq. - 6) then 
         WRITE (22, * ) 'set terminal postscript monochrome' 
         WRITE (22, * ) 'set output "|lpr -h"' 
      ELSEIF (idev.eq. - 7) then 
         WRITE (22, * ) 'set terminal postscript color' 
         WRITE (22, * ) 'set output "giffvc.eps"' 
      ELSE 
         WRITE ( * , * ) ' this device flag ', idev, ' not known' 
         RETURN 
      ENDIF 
      WRITE (22,  * ) 'plot ''giffv.tmp''  with ', style 
!  pause for screen images                                              
      IF (idev.gt.0) then 
         WRITE (22, * ) 'pause -1' 
      ENDIF 
      CLOSE (22) 
!  on SONY NEWS system is a function                                    
!     ii=system('gnuplot giffv.gnu')                                    
!     write(*,*)ii                                                      
!  on IBM RISC system is a subroutine                                   
      IF (idev.eq.4) then 
         CALL system ('xterm -t -e gnuplot giffv.gnu') 
      ELSE 
         CALL system ('gnuplot giffv.gnu') 
      ENDIF 
!  hard copy (if not done already)                                      
      IF (idev.eq. - 7) then 
!        ii=system('lpr giffv.eps')                                      
!        write(*,*)ii                                                   
         CALL system ('lpr giffvc.eps') 
      ELSEIF (idev.eq. - 5) then
!  option moved in front of filename 
         CALL system ('lpr -h giffvbw.eps') 
      ENDIF 
!  ADDED removal of giffv.tmp and giffv.gnu to prevent program crush
!      CALL system ('rm -f giffv.tmp')
!      CALL system ('rm -f giffv.gnu') 
      RETURN 
      END SUBROUTINE plotcr                         
! ***************************************************************       
! dummy subroutines (with GNUPLOT each page is a separate run           
! of GNUPLOT, including automatic initialisation)                       
!  {\bf ploflu} flush graphics device                                   
      SUBROUTINE ploflu (idev) 
      RETURN 
      END SUBROUTINE ploflu                         
! ***************************************************************       
!  {\bf newpag} new page/erase; GNUPLOT version                         
      SUBROUTINE newpag (idev, ipag) 
      RETURN 
      END SUBROUTINE newpag                         
! *************************************************************         
!  {\bf grafcl} close graphic device; GNUPLOT version                   
      SUBROUTINE grafcl (idev) 
      RETURN 
      END SUBROUTINE grafcl                         
! **************************************************************        
!  {\bf grafin} graphic device initialisation; GNUPLOT version          
      SUBROUTINE grafin (idev) 
      RETURN 
      END SUBROUTINE grafin                         


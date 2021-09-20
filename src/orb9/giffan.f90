!!$ -*- f90 -*-
!  *******************************************                          
!  {\sl GIFFAN, data analysis subroutines for GIFF}                     
!  *******************************************                          
!  {\bf peri2}                                                          
!                                                                       
!  peri1 - scritto da christian veillet 1981                            
!        ref.ferraz-mello a.j. 1981 vol.86 p.619-624                    
!  aggiornamento peri2- milani 1985                                     
!  con algoritmo basato su formule di addizione                         
!  e senza pesi- molto piu' rapido se tutti gli intervalli              
!  sono uguali (tolleranza=cambiamento relativo di 1.e-4)               
!  inoltre non si assume media zero e non si                            
!  alterano i dati se itest.eq.0                                        
!      input:t(n) tempi yy(n) dati                                      
!            per periodo                                                
!      workspace:vs(n),vc(n),y(n)                                       
!      output: d0 termine costante d1 coseno d2 seno                    
!             sp  densita' spettrale alla freq. data                    
!             yy(n) residui se itest.ne.0 ,altrim.immutati              
!  ************************************************                     
      SUBROUTINE peri2 (t, yy, n, per, sp, itest, d0, d1, d2)
      USE var_precision
      USE fund_const 
      IMPLICIT REAL(KIND=r_kind) (a - h, o - z) 
      REAL (KIND=r_kind), INTENT(in) :: t (n)
      REAL (KIND=r_kind), INTENT(inout) ::  yy (n) 
      REAL (KIND=r_kind) y (n), vs (n), vc (n) 
      xk = dpig  / per 
!    inizializzazione delle somme                                       
      sw = 0. 
      swc = 0. 
      swc2 = 0. 
      sws = 0. 
      sws2 = 0. 
      swcs = 0. 
      swcy = 0. 
      swsy = 0. 
      swy2 = 0. 
!    normalizzazione e riduzione a media pesata zero     
      swy= SUM (yy)               
      sw = n 
      ym = swy / sw 
      y(1:n)=yy(1:n) - ym 
!     fast fourier transform-type computation                           
      tcur = t (1) 
      x = xk * t (1) 
      vc (1) = cos (x) 
      vs (1) = sin (x) 
      dt = t (2) - t (1) 
      tol = abs (dt * 1.e-4) 
      xd = xk * dt 
      cd = cos (xd) 
      sd = sin (xd) 
!    sommatorie                                                         
      DO  i = 1, n 
         IF (abs (t (i) - tcur) .GT.tol) THEN 
            tcur = t (i) 
            x = tcur * xk 
            vc (i) = cos (x) 
            vs (i) = sin (x) 
         ENDIF
         c = vc (i) 
         s = vs (i) 
         swc = swc + c 
         swc2 = swc2 + c**2 
         sws = sws + s 
         sws2 = sws2 + s**2 
         swcs = swcs + c * s 
         swcy = swcy + c * y (i) 
         swsy = swsy + s * y (i) 
         swy2 = swy2 + y (i) **2 
         IF (i.eq.n) EXIT 
         vc (i + 1) = c * cd - s * sd 
         vs (i + 1) = s * cd + c * sd 
         tcur = tcur + dt 
      END DO 
!   calcolo densita' spettrale                                          
      a0 = 1 / sqrt (sw) 
      a1 = 1 / sqrt (swc2 - a0**2 * swc**2) 
      a2 = sws2 - a0**2 * sws**2 - a1**2 * swcs**2 - a0**4 * a1**2 *    &
      (swc * sws) **2                                                   
      a2 = a2 + 2 * a0**2 * a1**2 * swc * sws * swcs 
      a2 = 1. / sqrt (a2) 
      c1 = a1 * swcy 
      c2 = a2 * swsy - a1 * a2 * c1 * (swcs - a0**2 * swc * sws) 
      sp = c1**2 + c2**2 
      sp = sp / swy2 
!    calcolo elemento periodico                                         
      d2 = a2 * c2 
      d1 = a1 * c1 + a2 * a1**2 * c2 * (a0**2 * swc * sws - swcs) 
      d0 = - a0**2 * (d1 * swc + d2 * sws) + ym 
!    filtraggio se itest.ne.0                                           
      IF (itest.eq.0) return 
      yy (1:n) = y (1:n) - (d0 + d1 * vc (1:n) + d2 * vs (1:n) ) 
      RETURN 
      END SUBROUTINE peri2                          
!  *******************************************                          
!  {\bf peri1}                                                          
!                                                                       
!  peri1 - scritto da christian veillet 1981                            
!        ref.ferraz-mello a.j. 1981 vol.86 p.619-624                    
!  versione che non presuppone i dati equispaziati                      
!  inoltre non si assume media zero e non si                            
!  alterano i dati se itest.eq.0                                        
!      input:t(n) tempi yy(n) dati                                      
!            per periodo                                                
!      workspace:vs(n),vc(n),y(n)                                       
!      output: d0 termine costante d1 coseno d2 seno                    
!             sp  densita' spettrale alla freq. data                    
!             yy(n) residui se itest.ne.0 ,altrim.immutati              
!  ************************************************                     
      SUBROUTINE peri1 (t, yy, n, per, sp, itest, d0, d1, d2) 
      USE var_precision
      IMPLICIT REAl(KIND=r_kind) (a - h, o - z) 
      DIMENSION t (n), yy (n), y (n), vs (n), vc (n) 
      xk = 2.d0 * 3.1415926535d0 / per 
!    inizializzazione delle somme                                       
      sw = 0. 
      swc = 0. 
      swc2 = 0. 
      sws = 0. 
      sws2 = 0. 
      swcs = 0. 
      swcy = 0. 
      swsy = 0. 
      swy2 = 0. 
!    normalizzazione e riduzione a media pesata zero 
      swy = SUM (yy)                   
           sw = n 
      ym = swy / sw 
      y (1:n) = yy (1:n) - ym 
!    sommatorie                                                         
      DO  i = 1, n 
         tcur = t (i) 
         x = tcur * xk 
         vc (i) = cos (x) 
         vs (i) = sin (x) 
         c = vc (i) 
         s = vs (i) 
         swc = swc + c 
         swc2 = swc2 + c**2 
         sws = sws + s 
         sws2 = sws2 + s**2 
         swcs = swcs + c * s 
         swcy = swcy + c * y (i) 
         swsy = swsy + s * y (i) 
         swy2 = swy2 + y (i) **2 
      END DO 
!   calcolo densita' spettrale                                          
      a0 = 1 / sqrt (sw) 
      a1 = 1 / sqrt (swc2 - a0**2 * swc**2) 
      a2 = sws2 - a0**2 * sws**2 - a1**2 * swcs**2 - a0**4 * a1**2 *    &
      (swc * sws) **2                                                   
      a2 = a2 + 2 * a0**2 * a1**2 * swc * sws * swcs 
      a2 = 1. / sqrt (a2) 
      c1 = a1 * swcy 
      c2 = a2 * swsy - a1 * a2 * c1 * (swcs - a0**2 * swc * sws) 
      sp = c1**2 + c2**2 
      sp = sp / swy2 
!    calcolo elemento periodico                                         
      d2 = a2 * c2 
      d1 = a1 * c1 + a2 * a1**2 * c2 * (a0**2 * swc * sws - swcs) 
      d0 = - a0**2 * (d1 * swc + d2 * sws) + ym 
!    filtraggio se itest.ne.0                                           
      IF (itest.eq.0) return 
      yy (1:n) = y (1:n) - (d0 + d1 * vc (1:n) + d2 * vs (1:n) ) 
      RETURN 
      END SUBROUTINE peri1                          
                     











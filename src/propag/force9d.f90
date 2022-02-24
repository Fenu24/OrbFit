MODULE force9d
  USE massmod
  IMPLICIT NONE 
  PRIVATE

! shared data

INTEGER, PUBLIC :: irelj2 ! flag for J2 and relativistic perturbations 

! public routines
PUBLIC force9, reord, stack

CONTAINS

! **********************************************************            
!   FORCE9D  ORBIT9D                                                    
!                                                                       
!       general purpose right hand side with J2 and relativistic corr.  
!                                                                       
!   (vers. N+M bodies, 3d, no regular., nvz var.eq., closapp control)   
! **********************************************************            
SUBROUTINE force9(x,v,t,f,nd,idc,xxpla,ips,imem) 
  USE planet_masses, ONLY: dmin
  USE yark_pert, ONLY:  sec_nong9
  USE dyn_param, ONLY: iyark
  USE fund_const
  USE yorp_module
  INCLUDE 'comnbo.h90' 
!  INCLUDE 'comdis.h90' 
! **********************************************************            
! nvar must be 6*(nbod-1+na+nvz); nvar2=nvar/2                          
! the velocity is available to allow for                                
! a velocity dependent right hand side; not used here                   
! warning: x begins with first planet; the coordinates of the Sun       
! are not among the dynamical variables. To better undertand the        
! allocation of the coordinates in the stack x, see stack, reord below  
! INPUT
  INTEGER, INTENT(IN) :: nd 
  DOUBLE PRECISION, INTENT(IN) :: x(nd),v(nd),t
  INTEGER, INTENT(IN) :: ips, imem ! for storage, so far dummy
! OUTPUT
  DOUBLE PRECISION, INTENT(OUT) :: f(nd)
! close approach data                                                   
  INTEGER, INTENT(OUT) :: idc 
  DOUBLE PRECISION, INTENT(OUT) :: xxpla(6) 
! END INTERFACE
  DOUBLE PRECISION dy(3),dyy(3,norbx),dd3(norbx),d2(norbx)
  DOUBLE PRECISION, DIMENSION(3):: s, sv, xb, vb, secacc
! variables added for J2 and relativistic computation                   
  DOUBLE PRECISION dr3(norbx),z2,x2y2,dr4,cj2r2,x2,y2,cj2r9,sc2r6,r6 
  DOUBLE PRECISION d2udx(3,3) 
! indexes for addresses                                                 
  INTEGER nvar,nvar2,norb,nvz1,nvz2,ia,iv,jad,ja,j1,j,jj 
! scalar temporaries                                                    
  DOUBLE PRECISION sum,st,d5,st1,d3 
  ! =========================
  ! VARIABLES FOR YORP EFFECT
  ! =========================
  ! Parameters for YORP vector field
  real(kind=dkind) :: rpar(5)
  ! Variables for integration and output
  real(kind=dkind) :: y_yorpnew(2)
  real(kind=dkind) :: y_out(2), t_aux
  real(kind=dkind) :: gamnew, Pnew
  real(kind=dkind) :: gamma_new, omega_new
  logical out_flag
  ! Variables for computation of the Yarkovsky drift
  real(kind=dkind) :: pos(3), vel(3), energy, sma
  real(kind=dkind) :: alpha, epsi
  real(kind=dkind) :: dadt
  ! Collision reorientation flag
  logical reorient_flag
  ! Mass fraction for mass shedding event and critical period for
  ! mass shedding
  real(kind=dkind) :: q_fact
  real(kind=dkind) :: crit_period
  ! Variables for the stochastic YORP
  real(kind=dkind) :: tau_yorp
  real(kind=dkind) :: elap_time
  ! Asteroid name
  character*9 astname

! address precomputations                                               
  nvar=6*(nbod-1+na+nvz) 
  nvar2=nd 
  IF(nvar2*2.ne.nvar)THEN 
     WRITE(*,*)' force9d: nvar,nvar2,nbod,na,nvz=',nvar,nvar2,nbod,na,nvz
     STOP 
  ENDIF
  norb=nbod-1+na 
! close approach monitor                                                
  idc=0 
! first and last variational equation                                   
  nvz1=nbod 
  nvz2=nbod-1+nvz 
  ia=3*nbod-3 
  iv=ia+3*na 
! position of the sun in the barycentric system                         
!  $$s=-\sum_{j=2}^{nbod} \frac{m_j}{m_1} x_j $$                        
! (nbod-1)* (3 M + 3 S) unrolled inner loop                             
  s=0.d0 
  DO j=1,nbod-1 
     s(1:3)=s(1:3)-rm(j)*x(3*j-2:3*j) 
  ENDDO 
  IF(iyark.eq.3)THEN
     sv=0.d0 
     DO j=1,nbod-1 
        sv(1:3)=sv(1:3)-rm(j)*v(3*j-2:3*j) 
     ENDDO
  ENDIF 
! acceleration from the sun                                             
!  $$f^\circ_j=\frac{G m_1}{|x_j-x_1|^3}(x_1-x_j) $$                    
! (nbod-1+na)(8 M + 5 S + 1 SQRT)  unrolled inner loop                  
!$DIR NO_RECURRENCE                                                     
  DO 33 j=1,norb 
     st=s(1)-x(3*j-2) 
     sum=st*st 
     dyy(1,j)=st 
     st=s(2)-x(3*j-1) 
     sum=st*st+sum 
     dyy(2,j)=st 
     st=s(3)-x(3*j) 
     sum=st*st+sum 
     dyy(3,j)=st 
     d2(j)=sum 
     dr3(j)=1.d0/(sum*sqrt(sum)) 
     dd3(j)=gm(1)*dr3(j) 
     f(3*j-2)=dd3(j)*dyy(1,j) 
     f(3*j-1)=dd3(j)*dyy(2,j) 
     f(3*j)=dd3(j)*dyy(3,j) 
33 ENDDO
! Yarkovsky effect
!  IF(iyark.eq.3.or.iyark.eq.4)THEN
!     DO jj=1,na
!        j=nbod+jj-1 ! loop on asteroid only
!        xb=x(3*j-2:3*j)
!        vb=v(3*j-2:3*j)
!        CALL sec_nong9(xb,vb,s,sv,secacc,iyark,dadt_My(jj))
!        f(3*j-2:3*j)=f(3*j-2:3*j)+secacc
!     ENDDO
!  ENDIF

  !====================================
  !     YARKOVSKY/YORP EFFECT          
  !====================================
  ! Author: Marco Fenucci
  !   Date: September 2021

  !=====================================================
  !            SPIN AXIS INTEGRATION STEP
  !=====================================================
  ! NOTE: do this only if the integration with Yarkovsky 
  ! force is required AND if spin-axis is required
  !
  ! Check if we need to make a step in (\omega, \gamma)
  out_flag = .false.
  if(iyark.eq.3.and.yorp_flag.eq.1.and.abs(t*y2d-t_yorp).gt.h_yorp)then
     ! Update the time first
     t_yorp = t_yorp + h_yorp
     ! If we have to make a step, update the states of (\omega,
     ! \gamma) for all the bodies, and write the output on a file.
     do j=1, na
        ! =========================================
        !    CHECK FOR COLLISION RE-ORIENTATION
        ! =========================================
        ! Check if collision reorientation is occurring for P > 1000h
        if(1.d0/y_yorp(j, 1)*d2h .ge. 1000.d0)then
           call reor_probability(y_yorp(j,1), D_ast(j), h_yorp, reorient_flag)
           ! If a collisional reorientation is occurring, generate new
           ! initial conditions for \gamma and \omega and take new
           ! functions f and g
           if(reorient_flag)then
              call spin_state_new(gamma_new, omega_new)
              y_yorp(j, 1) = omega_new
              y_yorp(j, 2) = gamma_new
              ! Choose new f and g
              call shape_gen(K_ast(j), coeff_ast(j,1), coeff_ast(j,2))
           endif
        endif
        ! ===========================================
        !      CHECK FOR MASS SHEDDING EVENTS       
        ! ===========================================
        ! Assume mass shedding if the rotation period is too small.
        ! For this case, we divide two cases: if the density is low, 
        ! we use a low cohesion P = 10 Pa, while if it is large we assume 
        ! a large cohesion of P = 1000 Pa. The critical spin limit
        ! at small size is computed using the analytical formulas by
        ! Hu et al. 2021.
        call crit_rot_period(rho_ast(j), D_ast(j), crit_period)
        if(1.d0/y_yorp(j,1)*d2h.lt.crit_period)then
           ! \gamma is unchanged, while \omega is computed anew
           ! Generate a random ratio for the mass shed
           call random_ratio(q_fact)
           omega_new = sqrt(abs(y_yorp(j,2)**2.d0 - K_fact*q_fact))
           y_yorp(j, 1) = omega_new
           ! Choose new f and g
           call shape_gen(K_ast(j), coeff_ast(j,1), coeff_ast(j,2))
        endif
        ! ===========================================
        !  CHECK FOR STOCHASTIC EVENTS, IF REQUIRED  
        ! ===========================================
        if(stoc_yorp_flag.eq.1)then
           ! Compute tau_yorp (in years) as a function of the size
           ! tau_yorp = 0.25 My * D(km)^(4/3)*3/4
           tau_yorp  = 0.25d6*(D_ast(j)/1000.d0)**(4.d0/3.d0)*c_STOC
           ! Compute the time elapsed from the last event, in year
           elap_time = abs(t_yorp*d2y-last_stoc_event(j)) 
           ! Check if a time of tau_yorp has elapsed from the last
           ! event
           if(elap_time.gt.tau_yorp)then
              ! If yes, choose new f and g
              call shape_gen(K_ast(j), coeff_ast(j,1), coeff_ast(j,2))
              ! Update the time of the last stochastic event, saved
              ! in years
              last_stoc_event(j) = t_yorp*d2y
           endif
        endif
        ! ===========================================
        !    INTEGRATION STEP OF SPIN DYNAMICS       
        ! ===========================================
        ! Take the parameters for YORP vector field
        rpar(1) = rho_ast(j)
        rpar(2) = sma_ast(j)
        rpar(3) =   D_ast(j)
        rpar(4) = coeff_ast(j, 1)
        rpar(5) = coeff_ast(j, 2)
        ! Update the state with rk4
        call rk4(y_yorp(j, 1:2), h_yorp, y_yorpnew, rpar)  
        ! ===========================================
        !        CHECK FOR REQUIRED OUTPUT           
        ! ===========================================
        if(  time_out .le. (t_yorp-tstart_copy)*d2y        .and. &
         &   time_out .ge. (t_yorp-h_yorp-tstart_copy)*d2y .and. &
         & enable_out .eq. 1)then
           ! If that is the case, write the current integration
           ! in the file. First of all, compute the output value
           ! by linear interpolation between the timestep at time t
           ! and the timestep at time t+h
           t_aux = (time_out-(t_yorp-h_yorp-tstart_copy)*d2y)/(h_yorp*d2y)
           y_out(1:2) = y_yorp(j, 1:2)*(1.d0-t_aux) + y_yorpnew(1:2)*t_aux
           ! Open the output file
           astname = id_copy(j)
           open(unit=127, file=astname(1:len_trim(astname))//'.yorp',status='old', access='append')
           ! Write period in h, and gamma in deg
           write(127,1001) time_out, (1.d0/y_out(1))*d2h, y_out(2)*rad2deg, dadt_MY(j)
           close(127)
           if(.not.out_flag)then
              out_flag = .true.
           endif
1001 format(4(e12.6,1x))
        endif
        ! Update the variable
        y_yorp(j, 1:2) = y_yorpnew(1:2)
        ! ===========================================
        !         UPDATE THE YARKOVSKY EFFECT        
        ! ===========================================
        alpha = alpha_ast(j)
        epsi  =  epsi_ast(j)
        ! Compute the new semimajor axis
        jj  = nbod + j-1 
        pos(1:3) = x(3*jj-2:3*jj) - s(1:3)
        vel(1:3) = v(3*jj-2:3*jj) - sv(1:3)
        energy   = 0.5d0*norm2(vel)**2 - gms*365.25d0**2.d0/norm2(pos)
        sma      = -gms*365.25d0**2.d0/(2.d0*energy)
        ! Update also the global variable with the new semimajor axis
        sma_ast(j) = sma
        ! Compute new gamma and new P
        Pnew   = 1.d0/y_yorp(j, 1)*d2h
        gamnew = y_yorp(j, 2)*rad2deg
        ! Update Yarkovsky drift
        if(energy.lt.0.d0)then
           call dadt_comp(rho_ast(j), K_ast(j), C_ast(j), &
                          0.5d0*D_ast(j), sma, gamnew,    &
                          Pnew, alpha, epsi, dadt)
        else
           dadt = 0.d0
        endif
        dadt_My(j) = dadt
     enddo
     ! If output was performed, increase the variable keeping track
     ! of the output time
     if(out_flag)then
        time_out = time_out + dt_out
     endif
  endif

  ! Add the Yarkovsky effect to the gravitational vector field
  if(iyark.eq.3.or.iyark.eq.4)then
     do jj=1,na
        j  = nbod + jj-1 ! loop on asteroid only
        xb = x(3*j-2:3*j)
        vb = v(3*j-2:3*j)
        CALL sec_nong9(xb,vb,s,sv,secacc,iyark,dadt_My(jj))
        f(3*j-2:3*j) = f(3*j-2:3*j)+secacc
     enddo
  endif

!====================================
  IF(irelj2.gt.0)THEN
! J2 and relativistic corrections                                       
! from Nobili et al., A&A 1989 vol. 210 pag313-336; page 316            
     DO 34 j=1,norb 
        dr4=1.d0/d2(j)**2 
        x2y2=dyy(1,j)**2+dyy(2,j)**2 
        z2=dyy(3,j)**2 
        cj2r2=cj2*dr3(j)*3.d0 
        f(3*j-2)=f(3*j-2)+dyy(1,j)*(                                    &
     &    cj2r2*(2.d0*z2-0.5d0*x2y2)-                                   &
     &    2.d0*scw)*dr4                                                 
        f(3*j-1)=f(3*j-1)+dyy(2,j)*(                                    &
     &    cj2r2*(2.d0*z2-0.5d0*x2y2)-                                   &
     &    2.d0*scw)*dr4                                                 
        f(3*j)=f(3*j)+dyy(3,j)*(                                        &
     &    cj2r2*(z2-1.5d0*x2y2)-                                        &
     &    2.d0*scw)*dr4                                                 
34   ENDDO
  ENDIF
! variational equations for all asteroids 
   IF(irelj2.eq.0)THEN
! two body formula
! $$(\derparz{f^\circ_j}{x_j})_{kh}=-\frac{Gm_1}{|x_j-x_1|^3}\delta_{kh}
!   +\frac{3Gm_1}{|x_j-x_1|^5}(x_1-x_j)_k(x_1-x_j)_h$$                  
! nvz*(14 M + 5 S)  
      DO 3 j=nvz1,nvz2 
         jad=iv+3*(j-nvz1) 
         d5=3.d0*dd3(j)/d2(j) 
         st=dyy(1,j)*x(jad+1)+dyy(2,j)*x(jad+2)+dyy(3,j)*x(jad+3) 
         f(1+jad)=d5*st*dyy(1,j)-dd3(j)*x(1+jad) 
         f(2+jad)=d5*st*dyy(2,j)-dd3(j)*x(2+jad) 
         f(3+jad)=d5*st*dyy(3,j)-dd3(j)*x(3+jad) 
3     ENDDO
   ELSE            
! relativistic and J2 correction included
     DO 30 j=nvz1,nvz2 
!ccc       if(j.ge.nvz1.and.j.le.nvz2)then                              
! used in var. eq. for unperturbed orbit                                
        jad=iv+3*(j-nvz1) 
        d5=3.d0*dd3(j)/d2(j) 
        st=dyy(1,j)*x(jad+1)+dyy(2,j)*x(jad+2)+dyy(3,j)*x(jad+3) 
! used in the variational equations for J2, rel                         
        r6=dr3(j)*dr3(j) 
        x2=dyy(1,j)*dyy(1,j) 
        y2=dyy(2,j)*dyy(2,j) 
        z2=dyy(3,j)*dyy(3,j) 
        cj2r9=-3.d0*cj2*dr3(j)*r6 
        sc2r6=2.d0*scw*r6 
        d2udx(1,1)=x2*(cj2r9*(2.d0*x2+1.5d0*y2-13.5d0*z2)-3.d0*sc2r6)+  &
     &             y2*(cj2r9*(1.5d0*z2-0.5d0*y2)+sc2r6)+                &
     &             z2*(cj2r9*2.d0*z2+sc2r6)                             
        d2udx(2,2)=x2*(cj2r9*(1.5d0*z2-0.5d0*x2+1.5d0*y2)+sc2r6)+       &
     &             y2*(cj2r9*(2.d0*y2-13.5d0*z2)-3.d0*sc2r6)+           &
     &             z2*(cj2r9*2.d0*z2+sc2r6)                             
        d2udx(3,3)=x2*(cj2r9*(12.d0*z2-1.5d0*x2-3.d0*y2)+sc2r6)+        &
     &             y2*(cj2r9*(12.d0*z2-1.5d0*y2)+sc2r6)+                &
     &             z2*(-cj2r9*4.d0*z2-3.d0*sc2r6)                       
        d2udx(1,2)=dyy(1,j)*dyy(2,j)*(cj2r9*(2.5d0*x2+2.5*y2-15.d0*z2)  &
     &             -4.d0*sc2r6)                                         
        d2udx(2,1)=d2udx(1,2) 
        d2udx(1,3)=dyy(1,j)*dyy(3,j)*(cj2r9*(7.5d0*x2+7.5*y2-10.d0*z2)  &
     &             -4.d0*sc2r6)                                         
        d2udx(3,1)=d2udx(1,3) 
        d2udx(2,3)=dyy(2,j)*dyy(3,j)*(cj2r9*(7.5d0*x2+7.5*y2-10.d0*z2)  &
     &             -4.d0*sc2r6)                                         
        d2udx(3,2)=d2udx(2,3) 
! right hand side of variational equations                              
        f(1+jad)=d5*st*dyy(1,j)-dd3(j)*x(1+jad)+                        &
     &    d2udx(1,1)*x(1+jad)+d2udx(1,2)*x(2+jad)+                      &
     &    d2udx(1,3)*x(3+jad)                                           
        f(2+jad)=d5*st*dyy(2,j)-dd3(j)*x(2+jad)+                        &
     &    d2udx(2,1)*x(1+jad)+d2udx(2,2)*x(2+jad)+                      &
     &    d2udx(2,3)*x(3+jad)                                           
        f(3+jad)=d5*st*dyy(3,j)-dd3(j)*x(3+jad)+                        &
     &    d2udx(3,1)*x(1+jad)+d2udx(3,2)*x(2+jad)+                      &
     &    d2udx(3,3)*x(3+jad)                                           
!ccc       endif                                                        
30   ENDDO
  ENDIF
! iteration on the number of planets                                    
  DO 6 j=1,nbod-1 
! planet-planet interaction                                             
!  $$ f_{ji}=\frac{x_j-x_i}{|x_j-x_i|^3}\ \ \ for\ i<j$$                
!  $$f_j=f^\circ_j+\sum_{i<j}Gm_i f_{ji} + \sum_{i>j}Gm_i f_{ij}$$      
! (nbod-1)*(nbod-2)/2 * (13 M +11 S + 1 SQRT) unrolled inner loop       
!ccccC$DIR NO_RECURRENCE                                                
     DO 7 j1=1,j-1 
! vector differences and distances                                      
        st=x(3*j-2)-x(3*j1-2) 
        sum=st*st 
        dy(1)=st 
        st=x(3*j-1)-x(3*j1-1) 
        sum=sum+st*st 
        dy(2)=st 
        st=x(3*j)-x(3*j1) 
        sum=sum+st*st 
        dy(3)=st 
        d3=1.d0/(sum*sqrt(sum)) 
! mutual perturbations                                                  
        st=d3*gm(j+1) 
        st1=d3*gm(j1+1) 
        f(3*j1-2)=f(3*j1-2)+st*dy(1) 
        f(3*j-2)=f(3*j-2)-st1*dy(1) 
        f(3*j1-1)=f(3*j1-1)+st*dy(2) 
        f(3*j-1)=f(3*j-1)-st1*dy(2) 
        f(3*j1)=f(3*j1)+st*dy(3) 
        f(3*j)=f(3*j)-st1*dy(3) 
7    ENDDO
! planet-asteroid interaction                                           
! na*(nbod-1)*(8 M + 8 S + 1 SQRT) unrolled inner loop                  
!$DIR NO_RECURRENCE                                                     
     DO 8 ja=1,na 
! vector differences and distances                                      
        st=x(3*j-2)-x(3*ja-2+ia) 
        sum=st*st 
        dy(1)=st 
        st=x(3*j-1)-x(3*ja-1+ia) 
        sum=sum+st*st 
        dy(2)=st 
        st=x(3*j)-x(3*ja+ia) 
        sum=sum+st*st 
        dy(3)=st 
        d3=gm(j+1)/(sum*sqrt(sum)) 
! perturbations                                                         
        f(3*ja-2+ia)=f(3*ja-2+ia)+dy(1)*d3 
        f(3*ja-1+ia)=f(3*ja-1+ia)+dy(2)*d3 
        f(3*ja+ia)=f(3*ja+ia)+dy(3)*d3 
!  distance control                                                     
        if(sum.le.dmin(j))then 
           if(idc.eq.0)then 
              idc=j+(nbod-1)*ja 
              xxpla(1:3)=x(3*j-2:3*j)
              xxpla(4:6)=v(3*j-2:3*j)
           else 
              write(*,*)'force: double encounter ',t,idc,j,ia 
           endif
        endif
! variational equations for some asteroids                              
! nvz*(nbod-1)*(14 M + 8 S)                                             
!c         if(ja.le.nvz)then                                            
        jad=iv+3*(ja-1) 
        d5=3.d0*d3/sum 
        st=dy(1)*x(jad+1)+dy(2)*x(jad+2)+dy(3)*x(jad+3) 
        f(1+jad)=d5*st*dy(1)-d3*x(1+jad)+f(1+jad) 
        f(2+jad)=d5*st*dy(2)-d3*x(2+jad)+f(2+jad) 
        f(3+jad)=d5*st*dy(3)-d3*x(3+jad)+f(3+jad) 
!c         endif                                                        
    8 ENDDO 
! end loop on number of planets                                         
6 ENDDO
! =================================================                     
!  flop summary                                                         
!                                                                       
! do 2: (nbod-1)* 8                                                     
!                                                                       
! do 3 :(nbod-1+na)(13 + SQRT) + nvz*19                                 
!                                                                       
! do 7 :(nbod-1)*(nbod-2)/2 * (24 + SQRT)                               
!                                                                       
! do 8: na*(nbod-1)*(16 + SQRT) + nvz*(nbod-1)*22                       
!                                                                       
! tot: nvz*(nbod*22 -3) + na*(nbod*(13+SQRT)-3) +                       
! (nbod-1)(nbod-2)/2 * (24 + SQRT) + (nbod-1)*8                         
END SUBROUTINE force9
! *********************************************************             
!  {\bf reord} ORB8V output: reordering of the arrays                   
SUBROUTINE reord(y2,nvar,yp,npla,ya,na,yv,nvz,ndim,ndim2) 
  IMPLICIT NONE 
  INTEGER, INTENT(IN) :: nvar,ndim,ndim2,na,nvz,npla
  DOUBLE PRECISION, INTENT(IN) :: y2(nvar)
  DOUBLE PRECISION, INTENT(OUT) :: yp(ndim2,npla),ya(ndim2,na),yv(ndim2,nvz) 
! END INTERFACE  
  INTEGER nbod,nvar2                
  integer jj,j,i 
  nvar2=nvar/2 
  nbod=npla+1 
!  planets: reordering of the arrays                                    
  DO 11 j=1,npla 
     jj=3*j-3 
     DO i=1,ndim 
        yp(i+ndim,j)=y2(nvar2+jj+i) 
        yp(i,j)=y2(i+jj)
     ENDDO
11 ENDDO
!  asteroids: reordering of the arrays                                  
  DO 12 j=1,na 
     jj=3*(nbod-1+j)-3 
     DO i=1,ndim 
        ya(i+ndim,j)=y2(nvar2+jj+i) 
        ya(i,j)=y2(i+jj)
     ENDDO
12 ENDDO
!  variational equations                                                
  DO 20 j=1,nvz 
     jj=3*(nbod-1+na+j)-3 
     DO i=1,ndim 
        yv(i+ndim,j)=y2(nvar2+jj+i) 
        yv(i,j)=y2(i+jj)
     ENDDO
20 ENDDO
END SUBROUTINE reord
! *********************************************************             
!  {\bf stack} ORB8V initialisation: reordering of the arrays           
SUBROUTINE stack(y1,nvar,yp,npla,ya,na,yv,nvz,ndim,ndim2) 
  IMPLICIT NONE 
  INTEGER, INTENT(IN) ::  nvar,npla,na,nvz,ndim,ndim2
  DOUBLE PRECISION, INTENT(IN) :: yp(ndim2,npla),ya(ndim2,na),yv(ndim2,nvz)
  DOUBLE PRECISION, INTENT(OUT) :: y1(nvar)
  ! END INTERFACE                                               
  integer jj,j,i,nvar2 
  nvar2=nvar/2 
!  planets                                                              
  DO 1 j=1,npla 
     jj=3*j-3 
     DO i=1,ndim 
        y1(nvar2+jj+i)=yp(i+ndim,j) 
        y1(i+jj)=yp(i,j)
     ENDDO
1 ENDDO
! asteroids                                                             
  DO 2 j=1,na 
     jj=3*(npla+j)-3 
     DO i=1,ndim 
        y1(nvar2+jj+i)=ya(i+ndim,j) 
        y1(i+jj)=ya(i,j)
     ENDDO
2 ENDDO
! variational equations                                                 
  DO 3 j=1,nvz 
     jj=3*(npla+na+j)-3 
     DO i=1,ndim 
        y1(nvar2+jj+i)=yv(i+ndim,j) 
        y1(i+jj)=yv(i,j)
     ENDDO
3 ENDDO
END SUBROUTINE stack


END MODULE force9d

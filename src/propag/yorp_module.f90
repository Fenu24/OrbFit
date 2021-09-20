! =============================================================================================
!       Author : Marco Fenucci 
!         Date : September 2021
!  Institution : Department of Astronomy, Faculty of Mathematics, University of Belgrade
! =============================================================================================
module yorp_module
   implicit none
   private
   ! For FORTRAN programming
   integer, parameter :: dkind=kind(1.d0)
   ! Mathematical Constants 
   real(kind=dkind), parameter :: pi = 4.d0*atan(1.d0) 
   real(kind=dkind), parameter :: duepi = 2.d0*pi
   real(kind=dkind), parameter :: deg2rad = pi/180.d0
   real(kind=dkind), parameter :: rad2deg = 180.d0/pi
   ! ===========================================
   ! Physical parameters: standard units [kg, m, s, K]
   ! ===========================================
   ! Universal gravitational constant [m^3 kg^-1 s^-2]
   real(kind=dkind), parameter :: uGc = 6.67408d-11 
   ! Gravitational parameters G*m [m3 s-2]
   real(kind=dkind), parameter :: gmSun = 1.32712440018d20
   ! Stefan-Boltzmann constant [W m^-2 K^-4] = [kg s^-3 K^-4]
   real(kind=dkind), parameter :: sBoltz = 5.670374419d-8
   ! Solar luminosity [W] = [kg m^2 s^-3]
   real(kind=dkind), parameter :: lumSun = 3.828d26
   ! Astronomical Unit [m]
   real(kind=dkind), parameter :: AU   = 149597870700.d0
   real(kind=dkind), parameter :: au2m = AU
   real(kind=dkind), parameter :: m2au = 1.d0/AU
   ! Speed of light [m s^-1]
   real(kind=dkind), parameter :: clight = 299792458.d0
   ! Conversions
   real(kind=dkind), parameter :: h2s = 3600.d0
   real(kind=dkind), parameter :: s2h = 1.d0/h2s
   real(kind=dkind), parameter :: d2s = 24.d0*h2s
   real(kind=dkind), parameter :: s2d = 1.d0/d2s
   real(kind=dkind), parameter :: d2h = 24.d0
   real(kind=dkind), parameter :: h2d = 1.d0/d2h
   real(kind=dkind), parameter :: y2d = 365.25d0
   real(kind=dkind), parameter :: d2y = 1.d0/365.25d0
   real(kind=dkind), parameter :: y2s = 365.25d0*d2s
   real(kind=dkind), parameter :: s2y = 1.d0/y2s
   real(kind=dkind), parameter :: yr2my = 1.d-6
   real(kind=dkind), parameter :: my2yr = 1/yr2my
   real(kind=dkind), parameter :: s2my  = s2y*yr2my
   real(kind=dkind), parameter :: my2s  = 1.d0/s2my
   real(kind=dkind), parameter :: my2d  = my2yr*y2d 
   real(kind=dkind), parameter :: d2my  = 1.d0/my2d 
   ! ============================================
   !              FOR YORP EFFECT                
   ! ============================================
   ! Variables for the spline interpolation
   integer   :: nf, ng
   real(kind=dkind), allocatable, dimension(:) :: gamma_f, gamma_g
   real(kind=dkind), allocatable, dimension(:) :: f_fun,   g_fun
   real(kind=dkind), allocatable, dimension(:) :: f_b, f_c, f_d
   real(kind=dkind), allocatable, dimension(:) :: g_b, g_c, g_d
   ! ============================================
   ! Global variables for the coupled evolution  
   ! of Yarkovsky drift and spin axis            
   ! ============================================
   ! This should be equal to nastx
   integer, parameter :: nmax = 1100
   real(kind=dkind)   ::         dadt_MY(nmax)
   real(kind=dkind)   ::         rho_ast(nmax)
   real(kind=dkind)   ::           K_ast(nmax)
   real(kind=dkind)   ::           C_ast(nmax)
   real(kind=dkind)   ::           D_ast(nmax)
   real(kind=dkind)   ::       gamma_ast(nmax)
   real(kind=dkind)   ::       omega_ast(nmax)
   real(kind=dkind)   ::           P_ast(nmax)
   real(kind=dkind)   ::       alpha_ast(nmax)
   real(kind=dkind)   ::        epsi_ast(nmax)
   real(kind=dkind)   ::         sma_ast(nmax)
   real(kind=dkind)   ::       coeff_ast(nmax, 2)
   real(kind=dkind)   :: last_stoc_event(nmax)
   ! For YORP model
   real(kind=dkind)   :: c_YORP, c_REOR, c_STOC
   ! Variables for the integration:
   !     t_yorp: current time
   !     h_yorp: time step
   integer          :: yorp_flag
   integer          :: stoc_yorp_flag
   real(kind=dkind) :: dt_out, time_out
   real(kind=dkind) :: tstart_copy
   real(kind=dkind) :: t_yorp
   real(kind=dkind) :: h_yorp
   real(kind=dkind) :: y_yorp(nmax, 2)
   integer          :: enable_out
   ! Parameter for spin-up mass shedding
   real(kind=dkind), parameter :: K_fact=20.d0*pi*uGc*2000.d0/6.d0/s2d**2.d0
   ! Copy of the names of the bodies
   character*9     :: id_copy(nmax)
   ! ============================================
   ! SET VISIBILITY OF SUBROUTINES AND VARIABLES
   ! ============================================
   ! Public subroutines
   public :: read_f_g_spline, rk4, dadt_comp, yorp_vf, rad0pi, init_random_seed, &
          &  reor_probability, spin_state_new, random_ratio, shape_gen, crit_rot_period
   ! Public variables
   public :: dadt_MY, rho_ast, K_ast, C_ast, D_ast, gamma_ast, omega_ast, P_ast, sma_ast, &
          &  alpha_ast, epsi_ast, yorp_flag, dt_out, time_out, stoc_yorp_flag, last_stoc_event, &
          &  K_fact, enable_out, c_YORP, c_REOR, c_STOC
   public :: t_yorp, h_yorp, y_yorp, id_copy, tstart_copy, coeff_ast
   public :: h2s, s2h, d2s, s2d, y2d, d2y, y2s, s2y, yr2my, s2my, my2s, h2d, d2h, deg2rad, rad2deg
   contains

   ! ==============================
   !  INTEGRATION OF SPIN DYNAMICS 
   ! ==============================

   ! PURPOSE: integration step using an explicit Runge-Kutta integrator of order 4 (RK4)
   ! INPUT
   !  - y    : the current state vector y = (\omega, \gamma)
   !  - h    : the timestep
   !  - rpar : parameters for the YORP vector field. 
   !              rpar(1) = rho
   !              rpar(2) = sma
   !              rpar(3) = diameter (m)
   ! OUTPUT
   !  - ytph : the timestep at time t+h 
   subroutine rk4(y, h, ytph, rpar)
      real(kind=dkind), intent(in)  :: h
      real(kind=dkind), intent(in)  :: y(2)
      real(kind=dkind), intent(out) :: ytph(2)
      real(kind=dkind), intent(in)  :: rpar(5)
      ! end interface
      real(kind=dkind) :: f1(2), f2(2), f3(2), f4(2)
      real(kind=dkind) :: d1, d2, d3, d0
      real(kind=dkind) :: aux
      ! NOTE: the YORP vector field is independent from the time
      call yorp_vf(y,             f1, rpar)
      call yorp_vf(y + h*f1/2.d0, f2, rpar)
      call yorp_vf(y + h*f2/2.d0, f3, rpar)
      call yorp_vf(y + h*f3,      f4, rpar)
      ytph = y + h*(f1 + 2.d0*f2 + 2.d0*f3 + f4)/6.d0
      ! NOTE: this is done to avoid having 1.d-300 in \gamma
      if(abs(ytph(2)).lt.1.d-10)then
         ytph(2) = 0.d0
      endif
      ! If the period is too long, put it to 1000
      if(1.d0/ytph(1)*d2h.gt.1000.d0)then
         ytph(1) = 1.d0/(1000.d0*h2d)
      endif
      ! NOTE: this is done to handle the cases that goes asymptotically too fast!
      !       Indeed, in the step of 50 years is too long and they are decelerating very fast,
      !       we simply set a rotation period of 1000 h and put \gamma to the closest asymptotic
      !       value. 
      ! NOTE: the failure of the step is recognized because \omega becomes negative, which does not 
      !       make any physical sense!
      if(ytph(1).lt.0.d0)then
         ytph(1) = 1.d0/(1000.d0*h2d)
         d1 = abs(y(2))
         d2 = abs(y(2) - pi/2.d0)
         d3 = abs(y(2) - pi)
         d0 = min(d1, d2, d3)
         if(d0.eq.d1)then
            ytph(2) = 0.d0
         elseif(d0.eq.d1)then
            ytph(2) = pi/2.d0
         else
            ytph(2) = pi
         endif
      endif
   end subroutine rk4

   ! PURPOSE: compute the vector field due to YORP effect for the evolution of (\gamma, \omega), 
   !          using equations (7), (8) of Morbidelli & Vokrouhlicky 2003.
   !          Recall that here
   !              [\gamma] = rad
   !              [\omega] = 1/d
   !              [f, g]   = 1/s/My
   ! INPUT
   !  - x    : the current state vector x = (\omega, \gamma)
   !  - rpar : parameters of the YORP vector field, as in rk4 routine
   ! OUTPUT
   !  - vf   : the vector field
   subroutine yorp_vf(x, vf, rpar)
      real(kind=dkind), intent(in)  :: x(2)
      real(kind=dkind), intent(out) :: vf(2)
      real(kind=dkind), intent(in)  :: rpar(5)
      ! end interface
      real(kind=dkind) :: omega, gam
      real(kind=dkind) :: f, g
      real(kind=dkind) :: aux
      real(kind=dkind) :: rho, sma, diam, fact 
      real(kind=dkind) :: coe1, coe2
      ! Take the variables, [omega] = 1/d, [gam] = rad
      omega = x(1)
      gam   = x(2)
      ! Put gamma in the interval [0, pi]
      call rad0pi(gam, aux)
      gam = aux
      ! Compute f and g
      call comp_f_g(gam, f, g)
      ! Compute the rescaling factor for f and g
      rho  = rpar(1)
      sma  = rpar(2)
      diam = rpar(3)
      coe1 = rpar(4)
      coe2 = rpar(5)
      ! Rescaling factor with semimajor axis, density and diameter
      fact = c_YORP*(2500.d0/rho)*(2.5d0/sma)**2.d0*(2000.d0/diam)**2.d0
      ! Compute the vector field. Note that f, g are in units 1/s/My, it needs to be
      ! converted into 1/d^2
      vf(1) = f/s2d/my2d*fact*coe1
      vf(2) = g/omega/my2d/s2d*fact*coe2
      ! If the rotation period is too long (> 1000 h) or too small (< 1 m), 
      ! then freeze the evolution of \omega
      if(1.d0/(omega/d2h).gt.1.d3 .or. 1.d0/(omega/d2h) .lt. 30.d0*s2h)then
         vf(1) = 0.d0 
         vf(2) = 0.d0 
      endif
   end subroutine

   ! PURPOSE: read the files yorp_f.txt, yorp_g needed to compute the vector field, 
   !          and compute the coefficients of the cubic spline
   subroutine read_f_g_spline
      integer          :: j
      real(kind=dkind) :: x,y
      ! Compute the length of the files
      open(unit=200, file='input/yorp_f.txt', action='read')
      do j=1,1000
         read(200, *, end=33) x, y 
      enddo
      33 continue
      close(200)
      nf = j-1
      open(unit=200, file='input/yorp_g.txt', action='read')
      do j=1,1000
         read(200, *, end=34) x, y 
      enddo
      34 continue
      close(200)
      ng = j-1
      ! Allocate the memory
      allocate(gamma_f(nf), gamma_g(ng))
      allocate(f_fun(nf),   g_fun(ng))
      allocate(f_b(nf), f_c(nf), f_d(nf))
      allocate(g_b(ng), g_c(ng), g_d(ng))
      ! Read and store the content of the files
      open(unit=200, file='input/yorp_f.txt', action='read')
      do j=1, nf
         read(200, *) gamma_f(j), f_fun(j)
      enddo
      close(200)
      gamma_f = gamma_f*deg2rad
      ! Compute the spline coefficients
      call spline(nf, gamma_f, f_fun, f_b, f_c, f_d)
      open(unit=200, file='input/yorp_g.txt', action='read')
      do j=1, ng
         read(200, *) gamma_g(j), g_fun(j)
      enddo
      close(200)
      gamma_g = gamma_g*deg2rad
      ! Compute the spline coefficients
      call spline(ng, gamma_g, g_fun, g_b, g_c, g_d)
   end subroutine

   ! PURPOSE: Compute the functions f and g at the point gam (in degrees) using
   !          the cubic spline.
   ! INPUT
   !   - gam  : the value of the obliquity, in radians
   ! OUTPUT
   !   - f, g : the value of the torques
   subroutine comp_f_g(gam, f, g)
      real(kind=dkind), intent(in)  :: gam
      real(kind=dkind), intent(out) :: f, g
      ! end interface
      real(kind=dkind) :: gammap
      real(kind=dkind) :: dg
      integer          :: j
      ! Set f and g to zero
      f = 0.d0
      g = 0.d0
      ! Translate the value of gamma in [0, pi]
      call rad0pi(gam, gammap)
      ! Interpolate the value of f
      do j=1, nf-1
         if(gammap.ge.gamma_f(j) .and. gammap.le.gamma_f(j+1))then
            dg = gammap - gamma_f(j)
            f = f_fun(j) + f_b(j)*dg + f_c(j)*dg**2 + f_d(j)*dg**3
         endif
      enddo
      ! Interpolate the value of g
      do j=1, ng-1
         if(gammap.ge.gamma_g(j) .and. gammap.le.gamma_g(j+1))then
            dg = gammap - gamma_g(j)
            g = g_fun(j) + g_b(j)*dg + g_c(j)*dg**2 + g_d(j)*dg**3
         endif
      enddo
   end subroutine

   ! =======================================
   !   FUNCTIONS TO HANDLE THE YORP CYCLES  
   ! =======================================

   ! PURPOSE: generate random coefficients to determine the functions f,g. To
   !          this purpose, we use the statistics presented in Capek & Vokrouhlicky 2004
   ! INPUT
   !   - K : the value of thermal conductivity of the asteroid, in W/m/K
   ! OUTPUT
   !   - coeff1, coeff2 : coefficients for f,g respectively
   subroutine shape_gen(K, coeff1, coeff2)
      real(kind=dkind), intent(in)  :: K
      real(kind=dkind), intent(out) :: coeff1, coeff2
      ! end interface
      real(kind=dkind), parameter :: K_t = 0.005d0
      real(kind=dkind), parameter :: max_g = 1.8d0/1.1d0
      real(kind=dkind), parameter :: min_g = 0.4d0/1.1d0
      real(kind=dkind), parameter :: std_g = abs(1.1d0 - 1.8d0/1.1d0)/3.d0
      real(kind=dkind), parameter :: max_f =  3.d0/2.d0 
      real(kind=dkind), parameter :: min_f = -3.d0/2.d0 
      real(kind=dkind), parameter :: std_f = 0.5d0**2.0
      real(kind=dkind) :: rand_num, rand_g(1), sgn
      ! If K <= K_t, we assume:
      !    - 80% probability to reach 0/180 (g, coeff2)
      !    - 40% probability to accelerate  (f, coeff1)
      ! If K > K_t,  we assume
      !    - 100% probability to reach 0/180 (g, coeff2)
      !    - 50% probability to accelerate   (f, coeff1)
      if(K.le.K_t)then
         ! Draw a Gaussian random number with average equal to 1
         ! and std equal to std_g
         call gauss_random(1, rand_g)
         coeff2 = rand_g(1)*std_g + 1.d0
         if(coeff2.gt.max_g)then
            coeff2 = max_g
         elseif(coeff2.lt.min_g)then
            coeff2 = min_g
         endif 
         ! Draw a random number to decide whether to reach 
         ! 0/180 or 90 deg. If rand_num < 0.8, then we reach
         ! 0/180 (hence coeff2 > 0 and we do nothing),
         ! otherwise we reach 90  (hence coeff1 < 0, i.e. we change sign)
         call random_number(rand_num)
         if(rand_num.gt.0.8d0)then
            coeff2 = -coeff2
         endif
         ! If the asymptotic state is 90 (i.e. coeff2 < 0), then we always decelerate
         ! the rotation rate! See Capek & Vokrouhlicky 2004, Fig 7.
         if(coeff2.lt.0.d0)then
            ! Draw a Gaussian with mean equal to 1 and std equal to 0.5^2
            call gauss_random(1, rand_g)
            sgn    = -1.d0
            coeff1 = sgn*(rand_g(1)*std_f + 1.d0)
            if(coeff1.gt.max_f)then
               coeff1 = max_f
            elseif(coeff1.lt.min_f)then
               coeff1 = min_f
            endif
         else
            ! If the asymtotitc state is 0/180, use the 60/40 statistics.
            ! Draw a number to decide wheter to asymptotically accelerate
            ! or decelerate. 
            !  - if rand_num <= 0.6 we decelerate (hence we generate coeff2 between 0 and max_f)
            !  - if rand_num  > 0.6 we accelerate (hence we generate coeff2 between min_f and 0)
            call random_number(rand_num)
            if(rand_num.le.0.6d0)then
               sgn =  1.d0
            else
               sgn = -1.d0
            endif
            ! Draw a Gaussian with mean equal to 1 and std equal to 0.5^2
            call gauss_random(1, rand_g)
            coeff1 = sgn*(rand_g(1)*std_f + 1.d0)
            if(coeff1.gt.max_f)then
               coeff1 = max_f
            elseif(coeff1.lt.min_f)then
               coeff1 = min_f
            endif
         endif
      else
         ! Draw a random number between min_f and max_f, with Gaussian probability
         ! Choose the sign first
         call random_number(rand_num)
         if(rand_num.gt.0.5)then
            sgn = 1.d0
         else
            sgn = -1.d0
         endif
         ! Draw a Gaussian with mean equal to 1 and std equal to 0.5^2
         call gauss_random(1, rand_g)
         coeff1 = sgn*(rand_g(1)*std_f + 1.d0)
         if(coeff1.gt.max_f)then
            coeff1 = max_f
         elseif(coeff1.lt.min_f)then
            coeff1 = min_f
         endif
         ! Draw a Gaussian random number with average equal to 1
         ! and std equal to std_g
         call gauss_random(1, rand_g)
         coeff2 = rand_g(1)*std_g + 1.d0
         if(coeff2.gt.max_g)then
            coeff2 = max_g
         elseif(coeff2.lt.min_g)then
            coeff2 = min_g
         endif 
      endif
   end subroutine

   ! PURPOSE: compute a probability for a collisional reorientation to occur. The units are:
   !          [omega]  = 1/days
   !          [D]      = meters
   !          [h_yorp] = days
   ! INPUT
   !   - omega  : the rotation rate
   !   - D      : asteroid's diameter
   !   - h_yorp : timestep of the spin-axis integrator
   ! OUTPUT
   !   - reorient_flag : logical variable saying if a collisional reorientation happened or not
   subroutine reor_probability(omega, D, h_yorp, reorient_flag)
      real(kind=dkind), intent(in)  :: omega
      real(kind=dkind), intent(in)  :: D
      real(kind=dkind), intent(in)  :: h_yorp
      logical, intent(out)          :: reorient_flag
      ! end interface
      real(kind=dkind), parameter   :: B = 84.5d3*y2d
      real(kind=dkind), parameter   :: beta1  = 5.0d0/6.0d0
      real(kind=dkind), parameter   :: beta2  = 4.0d0/3.0d0
      real(kind=dkind), parameter   :: omega0 = 1.d0/(5.d0*h2d)
      real(kind=dkind), parameter   :: D0 = 2.0d0
      real(kind=dkind)              :: t_reor
      real(kind=dkind)              :: rand_num
      ! Compute the characteristic timescale for collisional reorientation.
      ! For details, see Broz et al 2011, Eq. (9)
      t_reor = c_reor*B*(omega/omega0)**beta1*(D/D0)**beta2
      ! Take a random number between 0 and 1
      call random_number(rand_num)
!      write(*,*) "t_reor ", t_reor*d2y/1.d6
!      write(*,*) " fact  ", 1.d0-exp(-h_yorp/t_reor)
!      write(*,*) "random ", rand_num
      ! Decide if reorientation is happening or not. For the strategy, 
      ! see Morbidelli & Vokrouhlicky 2003
      if(rand_num .lt. 1.d0-exp(-h_yorp/t_reor))then
         reorient_flag = .true.
      else
         reorient_flag = .false.
      endif
   end subroutine 

   ! PURPOSE: generate a new initial condition for the spin state.
   !          The obliquity is such that cos\gamma is equidistributed,
   !          while the rotation period is Gaussian with peak at 6 h and
   !          standard deviation of 4 h. Values smaller than 4 h and larger than 24 h 
   !          are discarded
   ! OUTPUT:
   !   - gamma_new : new wobliquity, in radians
   !   - omega_new : new rotation rate, in 1/days
   subroutine spin_state_new(gamma_new, omega_new)
      real(kind=dkind), intent(out) :: gamma_new
      real(kind=dkind), intent(out) :: omega_new
      ! end interface
      real(kind=dkind)              :: rand_num
      real(kind=dkind)              :: rand_gauss(3)
      real(kind=dkind)              :: P
      real(kind=dkind), parameter   :: P_peak = 10.d0
      ! Generate a new value of gamma
      call random_number(rand_num)
      gamma_new = acos(2.d0*rand_num - 1.d0)
      ! Generate the rotation period. The Maxwellian disribution is obtained 
      ! as a sum of three Gaussian distributions.
      P = 0.0
      ! Generate three numbers with standard Gaussian distribution
      call gauss_random(3, rand_gauss)
      ! Transform to Maxwellian
      P = P_peak*sqrt(rand_gauss(1)**2.d0 + rand_gauss(2)**2.d0 + rand_gauss(3)**2.d0)
      ! Transform to rotation rate
      omega_new = 1.d0/(P*h2d)
   end subroutine

   ! PURPOSE: compute the critical rotation period as given by Hu et al. 2021. 
   !          When the period returned by the formula is too large, we set a 
   !          hard limit of 2.44 h.
   ! INPUT
   !   - rho : density of the asteroid, in kg/m^2
   !   - D   : diameter of the asteroid, in meters
   ! OUTPUT
   !   - crit_period : critical rotation period, in hours
   subroutine crit_rot_period(rho, D, crit_period)
      real(kind=dkind), intent(in)  :: rho, D
      real(kind=dkind), intent(out) :: crit_period
      ! end interface
      real(kind=dkind), parameter :: phi = 32.5d0*deg2rad
      real(kind=dkind), parameter :: k = (1.d0+sin(phi))/(2.d0*cos(phi))
      real(kind=dkind)            :: T_c
      real(kind=dkind)            :: C
      ! Assume a constant cohesion of 100 Pa
      C = 100.d0
      ! Compute the critical rotation period with the formula by Hu et al. 2021
      ! Here D is in m and C in Pa, hence we have to convert it in hours
      T_c = duepi*sqrt(rho*k/C/5.d0)*D*s2h
      if(T_c.gt.2.44d0)then
         crit_period = 2.44d0
      else
         crit_period = T_c
      endif
   end subroutine 

   ! PURPOSE: generate a random mass ration between 0.002 and 0.02
   ! OUTPUT 
   !   - q : mass ratio
   subroutine random_ratio(q)
      real(kind=dkind), intent(out) :: q
      ! end interface
      real(kind=dkind)            :: rand_num, ee
      real(kind=dkind), parameter :: a = log10(0.002d0)
      real(kind=dkind), parameter :: b = log10(0.2d0)
      ! Generate a random number
      call random_number(rand_num)
      ee = (b-a)*rand_num + a
      q  = 10.d0**ee
   end subroutine

   ! =======================================
   !   COMPUTATION OF THE YARKOVSKY DRIFT   
   !            VOKROUHLICKY 1999           
   ! =======================================

   ! PURPOSE: compute the semimajor axis drift da/dt in au/my. The units of the quantities in
   !          input are as follows
   !           [rho]     = km/m^3
   !           [K]       = W/m/K
   !           [C]       = J/kg/K
   !           [radius]  = m
   !           [semiaxm] = au
   !           [gam]     = deg
   !           [rotPer]  = h
   ! INPUT
   !   - rho     : density of the asteroid
   !   - K       : thermal conductivity
   !   - C       : heat capacity
   !   - radius  : radius of the asteroid
   !   - semiaxm : semimajor axis of the asteroid's orbit
   !   - gam     : obliquity of the asteroid
   !   - rotPer  : rotation period of the asteroid
   !   - alpha   : absorption coefficient
   !   - epsi    : emissivity coefficient
   ! OUTPUT
   !   - yarko   : semimajor axis drift due to Yarkovsky effect
   subroutine dadt_comp(rho, K, C, radius, semiaxm, gam, rotPer, alpha, epsi, yarko)
      real(kind=dkind), intent(in)  :: rho, K, C, radius, semiaxm, rotPer, gam, alpha, epsi
      real(kind=dkind), intent(out) :: yarko
      ! end interface
      real(kind=dkind) :: dad, das
      call yarko_seasonal_vokrouhlicky(rho, K, C, radius, semiaxm, gam, rotPer, alpha, epsi, das)
      call yarko_diurnal_vokrouhlicky(rho, K, C, radius, semiaxm, gam, rotPer, alpha, epsi, dad)
      yarko = (das+dad)*m2au/s2my
   end subroutine dadt_comp

   ! PURPOSE: compute the seasonal component of the Yarkovsky effect. 
   !          Input and output variables are the same as in dadt_comp
   subroutine yarko_seasonal_vokrouhlicky(rho, K, C, R, a0, gam, rotPer, alpha, epsi, yarko_s)
      real(kind=dkind), intent(in)  :: rho, K, C, R, a0, gam, rotPer, alpha, epsi
      real(kind=dkind), intent(out) :: yarko_s
      ! end interface
      real(kind=dkind) :: mAst, mu
      real(kind=dkind) :: gam_rad, omega_rev
      real(kind=dkind) :: ls 
      real(kind=dkind) :: Theta, phi
      real(kind=dkind) :: E_star, T_star
      real(kind=dkind) :: dist
      real(kind=dkind) :: F_omega_rev
      real(kind=dkind) :: Rps
      ! Convert gamma in radians
      gam_rad = gam*deg2rad
      ! Compute mass of the asteroid
      mAst = 4.d0*pi*rho*R**3.d0/3.d0
      ! Compute omega_rev 
      mu = gmsun+uGc*mAst
      omega_rev = sqrt(mu/(a0*au2m)**3.d0)
      ls = sqrt(K/(rho*C*omega_rev))
      ! Compute R'
      Rps = R/ls
      ! Convert distance in meters
      dist   = a0*au2m
      ! Compute the solar radiation flux at distance dist
      E_star = lumSun/(4.d0*pi*dist**2.d0)
      ! Compute the radiation force
      phi    = pi*R**2.d0*E_star/(mAst*cLight)
      ! Compute the subsolar temperature
      T_star = (alpha*E_star/(epsi*sBoltz))**0.25d0
      ! Compute the thermal parameter
      Theta = sqrt(rho*K*C*omega_rev)/(epsi*sBoltz*T_star**3.d0)
      ! Compute F_omega_rev
      call Fnu_eval(Rps, Theta, F_omega_rev)
      ! Compute the Yarkovsky acceleration
      yarko_s = 4.d0*alpha*phi*F_omega_rev*sin(gam_rad)**2.d0/(9.d0*omega_rev)
   end subroutine yarko_seasonal_vokrouhlicky

   ! PURPOSE: compute the diurnal component of the Yarkovsky effect.
   !          Input and output variables are the same as in dadt_comp
   subroutine yarko_diurnal_vokrouhlicky(rho, K, C, R, a0, gam, rotPer, alpha, epsi, yarko_d)
      real(kind=dkind), intent(in)  :: rho, K, C, R, a0, gam, rotPer, alpha, epsi
      real(kind=dkind), intent(out) :: yarko_d
      ! end interface
      real(kind=dkind) :: mAst, mu
      real(kind=dkind) :: gam_rad, omega_rot, omega_rev
      real(kind=dkind) :: ls, ld
      real(kind=dkind) :: phi
      real(kind=dkind) :: Rpd, Rsmall, Rbig
      real(kind=dkind) :: Theta
      real(kind=dkind) :: E_star, T_star
      real(kind=dkind) :: dist
      real(kind=dkind) :: F_omega_rot
      ! Convert gamma in radians
      gam_rad = gam*deg2rad
      ! Compute mass of the asteroid
      mAst = 4.d0*pi*rho*R**3.d0/3.d0
      ! Compute omega_rev 
      mu = gmsun+uGc*mAst
      omega_rev = sqrt(mu/(a0*au2m)**3.d0)
      omega_rot = duepi/(rotPer*h2s)
      ld = sqrt(K/(rho*C*omega_rot))
      ! Compute R'
      Rpd = R/ld
      ! Convert distance in meters
      dist   = a0*au2m
      ! Compute the solar radiation flux at distance dist
      E_star = lumSun/(4.d0*pi*dist**2.d0)
      ! Compute the radiation force
      phi    = pi*R**2.d0*E_star/(mAst*cLight)
      ! Compute the subsolar temperature
      T_star = (alpha*E_star/(epsi*sBoltz))**0.25d0
      ! Compute the thermal parameter
      Theta = sqrt(rho*K*C*omega_rot)/(epsi*sBoltz*T_star**3.d0)
      ! Compute F_omega_rev
      call Fnu_eval(Rpd, Theta, F_omega_rot)
      ! Compute the Yarkovsky acceleration
      yarko_d = -8.d0*alpha*phi*F_omega_rot*cos(gam_rad)/(9.d0*omega_rev)
   end subroutine yarko_diurnal_vokrouhlicky

   ! PURPOSE: compute the value of the function F_\nu used in Vokrouhlicky 1999
   ! INPUT
   !   - R     : radius of the asteroid, in meters
   !   - Theta : thermal parameter
   ! OUTPUT
   !   - Fnu : the value of the function
   subroutine Fnu_eval(R, Theta, Fnu)
      real(kind=dkind), intent(in)  :: R, Theta
      real(kind=dkind), intent(out) :: Fnu
      ! end interface
      real(kind=dkind) :: k1, k2, k3
      real(kind=dkind) :: Rb, Rs
      real(kind=dkind) :: x, chi
      real(kind=dkind) :: A, B, C, D, U, V, den
      real(kind=dkind) :: xx1, xx2
      ! Formulation by Vokrouhlicky 1999 
      call k1k2k3_eval(R, Theta, k1, k2, k3)
      Fnu = -k1*Theta/(1.d0 + 2.d0*k2*Theta + k3*Theta**2.d0)  
   end subroutine Fnu_eval

   ! PURPOSE: compute the coefficients needed for the evaluation of F_\nu
   ! INPUT
   !   - R     : radius of the asteroid, in meters
   !   - Theta : thermal parameter
   ! OUTPUT
   !   - k1, k2, k3 : the three coefficients
   subroutine k1k2k3_eval(R, Theta, k1, k2, k3)
      real(kind=dkind), intent(in)  :: R, Theta
      real(kind=dkind), intent(out) :: k1, k2, k3
      ! end interface
      real(kind=dkind) :: Rb, Rs
      real(kind=dkind) :: x, chi
      real(kind=dkind) :: A, B, C, D, U, V, den
      real(kind=dkind) :: xx1, xx2
      ! Formulation by Vokrouhlicky 1999 
      if(R.gt.30d0)then
         k1  = 0.5d0
         k2  = 0.5d0
         k3  = 0.5d0
      else
         x   = sqrt(2.d0)*R
         A = -(x+2.d0) - exp(x)*((x-2.d0)*cos(x) - x*sin(x))
         B = -x - exp(x)*(x*cos(x) + (x-2.d0)*sin(x))
         U = 3.d0*(x+2.d0) + exp(x)*(3.d0*(x-2.d0)*cos(x)+x*(x-3.d0)*sin(x))
         V = x*(x+3.d0) - exp(x)*(x*(x-3.d0)*cos(x) - 3.d0*(x-2.d0)*sin(x))
         den = x*(A**2.d0+B**2.d0);
         k1 = (A*V-B*U)/den
         k2 = (A*(A+U) + B*(B+V))/den
         k3 = ((A+U)**2.d0 + (B+V)**2.d0)/(den*x)
      endif
   end subroutine k1k2k3_eval

   ! =======================================
   !         GENERAL SUBROUTINES            
   ! =======================================

   ! PURPOSE: Generate a random vector of real number with Gaussian distribution, 
   !          with mean zero and unitary standard deviation, using the Box-Muller 
   !          algorithm
   ! INPUT
   !   - n : dimension of the vector in output
   !   - u : vector of randomly gaussian numbers
   subroutine gauss_random(n, u)
      integer, intent(in)           :: n
      real(kind=dkind), intent(out) :: u(n)
      ! end interface
      real(kind=dkind) :: u1, u2
      integer :: i
      do i=1, n
         call random_number(u1)
         call random_number(u2)
         u(i) = sqrt(-2.d0*log(u1))*cos(duepi*u2) 
      enddo
   end subroutine

   ! PURPOSE: Subroutine to get angles (passed in radians) in [0,pi] interval  
   subroutine rad0pi(lambda,chl) 
      real(kind=dkind), intent(in)  :: lambda 
      real(kind=dkind), intent(out) :: chl 
      ! end interface 
      real(kind=dkind) :: x,y 
      x = sin(lambda) 
      y = cos(lambda) 
      !     chl in [-pi,pi]
      chl = atan2(x,y) 
      !     chl in [0,360]
      if((chl.ge.0.d0).and.(chl.le.pi)) then 
         chl = chl 
      elseif((chl.le.0.d0).and.(chl.ge.-pi)) then 
         chl = duepi + chl 
      else 
         write(*,*)'error!!!!, chl ', chl 
         write(*,*) "lambda = ", lambda
      endif
      !     control   
      if(chl.gt.duepi) then 
         write(*,*)'error error error !!!' 
      endif

      if(chl.gt.pi)then
         chl = chl - pi
      endif
   end subroutine rad0pi

   ! PURPOSE: cubic spline interpolation
   subroutine spline (n, x, y, b, c, d) 
     integer n 
     double precision x(n), y(n), b(n), c(n), d(n) 
   !                                                                      
   !  the coefficients b(i), c(i), and d(i), i=1,2,...,n are computed     
   !  for a cubic interpolating spline                                    
   !                                                                      
   !    s(x) = y(i) + b(i)*(x-x(i)) + c(i)*(x-x(i))**2 + d(i)*(x-x(i))**3 
   !                                                                      
   !    for  x(i) .le. x .le. x(i+1)                                      
   !                                                                      
   !  input..                                                             
   !                                                                      
   !    n = the number of data points or knots (n.ge.2)                   
   !    x = the abscissas of the knots in strictly increasing order       
   !    y = the ordinates of the knots                                    
   !                                                                      
   !  output..                                                            
   !                                                                      
   !    b, c, d  = arrays of spline coefficients as defined above.        
   !                                                                      
   !  using  p  to denote differentiation,                                
   !                                                                      
   !    y(i) = s(x(i))                                                    
   !    b(i) = sp(x(i))                                                   
   !    c(i) = spp(x(i))/2                                                
   !    d(i) = sppp(x(i))/6  (derivative from the right)                  
   !                                                                      
   !  the accompanying function subprogram  seval  can be used            
   !  to evaluate the spline.                                             
         integer nm1, ib, i 
         double precision t 
   !                                                                      
         nm1 = n-1 
         if ( n .lt. 2 ) return 
         if ( n .lt. 3 ) go to 50 
   !                                                                      
   !  set up tridiagonal system                                           
   !                                                                      
   !  b = diagonal, d = offdiagonal, c = right hand side.                 
   !                                                                      
         d(1) = x(2) - x(1) 
         c(2) = (y(2) - y(1))/d(1) 
         do 10 i = 2, nm1 
            d(i) = x(i+1) - x(i) 
            b(i) = 2.*(d(i-1) + d(i)) 
            c(i+1) = (y(i+1) - y(i))/d(i) 
            c(i) = c(i+1) - c(i) 
      10 continue 
   !                                                                      
   !  end conditions.  third derivatives at  x(1)  and  x(n)              
   !  obtained from divided differences                                   
   !                                                                      
         b(1) = -d(1) 
         b(n) = -d(n-1) 
         c(1) = 0. 
         c(n) = 0. 
         if ( n .eq. 3 ) go to 15 
         c(1) = c(3)/(x(4)-x(2)) - c(2)/(x(3)-x(1)) 
         c(n) = c(n-1)/(x(n)-x(n-2)) - c(n-2)/(x(n-1)-x(n-3)) 
         c(1) = c(1)*d(1)**2/(x(4)-x(1)) 
         c(n) = -c(n)*d(n-1)**2/(x(n)-x(n-3)) 
   !                                                                      
   !  forward elimination                                                 
   !                                                                      
      15 do 20 i = 2, n 
            t = d(i-1)/b(i-1) 
            b(i) = b(i) - t*d(i-1) 
            c(i) = c(i) - t*c(i-1) 
      20 continue 
   !                                                                      
   !  back substitution                                                   
   !                                                                      
         c(n) = c(n)/b(n) 
         do 30 ib = 1, nm1 
            i = n-ib 
            c(i) = (c(i) - d(i)*c(i+1))/b(i) 
      30 continue 
   !                                                                      
   !  c(i) is now the sigma(i) of the text                                
   !                                                                      
   !  compute polynomial coefficients                                     
   !                                                                      
         b(n) = (y(n) - y(nm1))/d(nm1) + d(nm1)*(c(nm1) + 2.*c(n)) 
         do 40 i = 1, nm1 
            b(i) = (y(i+1) - y(i))/d(i) - d(i)*(c(i+1) + 2.*c(i)) 
            d(i) = (c(i+1) - c(i))/d(i) 
            c(i) = 3.*c(i) 
      40 continue 
         c(n) = 3.*c(n) 
         d(n) = d(n-1) 
         return 
   !                                                                      
      50 b(1) = (y(2)-y(1))/(x(2)-x(1)) 
         c(1) = 0. 
         d(1) = 0. 
         b(2) = b(1) 
         c(2) = 0. 
         d(2) = 0. 
         return 
   end subroutine                                           

   ! PURPOSE: initial seed for the random number generator
   subroutine init_random_seed()
         integer :: i, n, clock
         integer, dimension(:), allocatable :: seed
         call random_seed(size = n)
         allocate(seed(n))
         call system_clock(count=clock)
         seed = clock + 37 * (/ (i - 1, i = 1, n) /)
         call random_seed(put = seed)
         deallocate(seed)
   end subroutine init_random_seed

end module yorp_module

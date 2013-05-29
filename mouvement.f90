subroutine mouvement(time, alpha, alpha_t, alpha_tt, LeadingEdge)
  use share_vars
  implicit none
  real (kind=pr), intent(out) :: alpha, alpha_t, alpha_tt 
  real (kind=pr), intent(in) :: time
  real (kind=pr), dimension(1:6), intent(out) :: LeadingEdge !LeadingEdge: x, y, vx, vy, ax, ay (Array)
  real (kind=pr) :: angle_max, t_max, tau1, u0, f, A0, beta, sigma_t, sigma_r, phi, G_translation, G_rotation, C_startup
  real (kind=pr) :: y_max, tau2, k, kt, ktt,a,b,c,d, y,yt,ytt
  integer :: n

  select case (iMotion)
    case (0) ! fixed beam (no leading edge mouvement)    
      !--------------------------------------------------------------------------------------------------------
      ! fixed
      !--------------------------------------------------------------------------------------------------------    
      LeadingEdge = 0.0
      LeadingEdge(1) = x0
      LeadingEdge(2) = y0
      alpha    = AngleBeam*pi/180.0
      alpha_t  = 0.0
      alpha_tt = 0.0
    case (1)
      !--------------------------------------------------------------------------------------------------------
      ! heaving: soft-startup version
      !--------------------------------------------------------------------------------------------------------
      y_max = 0.5
      t_max = 1.4 
     
      if (time <= 1.0) then
	a = -20.0; b= 70.0; c=-84.0; d=35.0;
	k    = a*time**7 + b*time**6 + c*time**5 + d*time**4
	kt  = 7.0*a*time**6 + 6.0*b*time**5 + 5.0*c*time**4 + 4.0*d*time**3
	ktt = 42.*a*time**5 + 30.*b*time**4 + 20.*c*time**3 + 12.*d*time**2	
      else
	k   = 1.0
	kt  = 0.0
	ktt = 0.0
      endif
      
      y   = y_max*sin(2.0*pi*time/t_max+pi/2.)
      yt  = y_max*(2.0*pi/t_max)*cos(2.0*pi*time/t_max+pi/2.)
      ytt =-y_max*(2.0*pi/t_max)**2*sin(2.0*pi*time/t_max+pi/2.)       

      LeadingEdge = 0.0
      LeadingEdge(1) = x0 
      LeadingEdge(2) = y0 + k*y
      LeadingEdge(3) = 0.0
      LeadingEdge(4) = kt*y + yt*k
      LeadingEdge(5) = 0.0
      LeadingEdge(6) = ktt*y + ytt*k + 2.*kt*yt
      
      alpha=AngleBeam*pi/180.0
      alpha_t=0.0
      alpha_tt=0.0
    case (2)
      !--------------------------------------------------------------------------------------------------------
      ! heaving
      !--------------------------------------------------------------------------------------------------------      
      y_max = 0.5
      t_max = 2.5
      LeadingEdge = 0.0
      LeadingEdge(1) = x0 
      LeadingEdge(2) = y0+y_max*sin(2.0*pi*time/t_max+pi/2)
      LeadingEdge(3) = 0.0
      LeadingEdge(4) = y_max*(2.0*pi/t_max)*cos(2.0*pi*time/t_max+pi/2)
      LeadingEdge(5) = 0.0
      LeadingEdge(6) = -y_max*(2.0*pi/t_max)**2*sin(2.0*pi*time/t_max+pi/2)
      alpha=AngleBeam*pi/180.0
      alpha_t=0.0
      alpha_tt=0.0
    case (3)
      !--------------------------------------------------------------------------------------------------------
      ! flapping
      !--------------------------------------------------------------------------------------------------------      
      angle_max = 80.0*pi/180.d0
      f = 0.5d0 ! frequency
      
      LeadingEdge = 0.d0
      LeadingEdge(1) = x0
      LeadingEdge(2) = y0
      
      alpha    = angle_max * sin(2.d0*pi*f*time)
      alpha_t  = angle_max * cos(2.d0*pi*f*time) * (2.d0*pi*f)
      alpha_tt = -1.d0 * angle_max * sin(2.d0*pi*f*time) * (2.d0*pi*f)**2
    case (4) 
      !--------------------------
      ! implusive translation
      !--------------------------
      u0 = -1.0
      LeadingEdge = 0.0
      LeadingEdge(1) = x0 + u0*time
      LeadingEdge(3) = u0 ! leading edge velocity
      LeadingEdge(2) = y0
      LeadingEdge(4) = 0.0
      alpha=AngleBeam*pi/180.0
      alpha_t=0.0
      alpha_tt=0.0

!     case (2,10) ! flapping motion (10: with active FSI coupling)
!       !--------------------------------------------------------------------------------------------------------
!       ! flapping
!       !--------------------------------------------------------------------------------------------------------          
!       angle_max = AngleBeam*pi/180.d0
!       f = 1.d0 ! frequency
!       
!       LeadingEdge = 0.d0
!       LeadingEdge(1) = x0
!       LeadingEdge(2) = y0
!       
!       alpha    = angle_max * sin(2.d0*pi*f*time)
!       alpha_t  = angle_max * cos(2.d0*pi*f*time) * (2.d0*pi*f)
!       alpha_tt = -1.d0 * angle_max * sin(2.d0*pi*f*time) * (2.d0*pi*f)**2
      !--------------------------------------------------------------------------------------------------------
      ! flapping, soft startup version
      !--------------------------------------------------------------------------------------------------------      
! ! !       if (time <= 1.0) then
! ! ! 	a = -20.0; b= 70.0; c=-84.0; d=35.0;
! ! ! 	k    = a*time**7 + b*time**6 + c*time**5 + d*time**4
! ! ! 	kt  = 7.0*a*time**6 + 6.0*b*time**5 + 5.0*c*time**4 + 4.0*d*time**3
! ! ! 	ktt = 42.*a*time**5 + 30.*b*time**4 + 20.*c*time**3 + 12.*d*time**2	
! ! !       else
! ! ! 	k   = 1.0
! ! ! 	kt  = 0.0
! ! ! 	ktt = 0.0
! ! !       endif     
! ! !       
! ! !       y    = angle_max * sin(2.d0*pi*f*time)
! ! !       yt   = angle_max * cos(2.d0*pi*f*time) * (2.d0*pi*f)
! ! !       ytt  = -1.d0 * angle_max * sin(2.d0*pi*f*time) * (2.d0*pi*f)**2
! ! !       
! ! !       alpha = k*y
! ! !       alpha_t = kt*y + yt*k
! ! !       alpha_tt = ktt*y + ytt*k + 2.*kt*yt      
      
      
      
!     case (3) ! impulsive translation
!       
!     case (101:107) !eldrege kinematics, case 1
!     write(*,*) "edrege won't work as mouvement is called sveral times per time step..."
!     stop
!       sigma_r=0.628
!       sigma_t=0.628
!       phi=0.0
!       call eldrege(time, alpha, alpha_t, alpha_tt, LeadingEdge, sigma_r, sigma_t, phi)
!     case (102) !eldrege kinematics, case 2
!       sigma_r=1.885
!       sigma_t=1.885
!       phi=0.0
!       call eldrege(time, alpha, alpha_t, alpha_tt, LeadingEdge, sigma_r, sigma_t, phi)
!     case (103) !eldrege kinematics, case 3
!       sigma_r=1.885
!       sigma_t=1.885
!       phi=pi/4.
!       call eldrege(time, alpha, alpha_t, alpha_tt, LeadingEdge, sigma_r, sigma_t, phi)
!     case (104) !eldrege kinematics, case 4
!       sigma_r=3.770
!       sigma_t=3.770
!       phi=0.0
!       call eldrege(time, alpha, alpha_t, alpha_tt, LeadingEdge, sigma_r, sigma_t, phi)
!     case (105) !eldrege kinematics, case 5
!       sigma_r=3.770
!       sigma_t=3.770
!       phi=pi/4.
!       call eldrege(time, alpha, alpha_t, alpha_tt, LeadingEdge, sigma_r, sigma_t, phi)
!     case (106) !eldrege kinematics, case 6
!       sigma_r=0.628
!       sigma_t=3.770
!       phi=0.0
!       call eldrege(time, alpha, alpha_t, alpha_tt, LeadingEdge, sigma_r, sigma_t, phi)  
!     case (107) !eldrege kinematics, case 6
!       sigma_r=3.770
!       sigma_t=0.628
!       phi=0.0
!       call eldrege(time, alpha, alpha_t, alpha_tt, LeadingEdge, sigma_r, sigma_t, phi) 
  end select

end subroutine

!===================================================================== 

subroutine eldrege(time, alpha, alpha_t, alpha_tt, LeadingEdge,sigma_r,sigma_t,phi)
  use share_vars
  implicit none
  real (kind=pr), intent(out) :: alpha, alpha_t, alpha_tt
  real (kind=pr), intent(in) :: time,sigma_r,sigma_t,phi
  real (kind=pr), dimension(1:6), intent(out) :: LeadingEdge !LeadingEdge: x, y, vx, vy, ax, ay (Array)
  real (kind=pr) :: G_translation, G_rotation, C_startup, X_cg, Y_cg, L, G_translation_max, dt_derivative=1.0e-5, G,C
  real (kind=pr) :: X_cg_plus, X_cg_minus, G_plus, G_minus, alpha_minus, alpha_plus,vx,vy,ax,ay,x,y,x_plus,y_plus,x_minus,y_minus
  real (kind=pr) :: C_plus,C_minus, beta, A0
  beta=pi/4.
  
  A0 = 0.7 ! note our normalization is different. our beam is L=2c
  L = 0.25 ! L_half = 0.5*c = 0.25*L = 0.25
  
!   this function computes eldreges kinematics. note that the original formulation deals with the center of gravity, while
!   we actually need the leading edge. thanks god that's simple. however, we also need the velocity and acceleration at the
!   leading edge, and we compute it with centered finite differences. 
!   note you cannot call G_translation 3 times (it has an internal counter)
!   the following relations hold:
!        xup=xcg+L*sin(alpha); leading edge
!        yup=L*cos(alpha);
!        xdo=xcg-L*sin(alpha); trailing edge
!        ydo=-L*cos(alpha);
  
  
  ! determine which case we are. note that ONLY the following cases are allowed
  ! these values are found numerically. the maximum of the function is very close to t=0.25
  if (abs(sigma_t-0.628)<1.e-2) G_translation_max= 9.216961e-02 ! found with help of MATLAB
  if (abs(sigma_t-1.885)<1.e-2) G_translation_max= 1.855609e-01
  if (abs(sigma_t-3.770)<1.e-2) G_translation_max= 2.199116e-01
  if (abs(G_translation_max)<1.0e-3) then
    write(*,*) "!!! subroutine mouvement->eldrege:  case for sigma_t not known"
    write(*,*) "for normalization, we need to know à priori the maximum value of the integral"
    stop
  endif
    
  
  G       = G_translation(time, sigma_r, sigma_t, phi, beta) * C_startup(time, sigma_r, sigma_t, phi, beta)
  G_plus  = G + dt_derivative * tanh(sigma_r*cos(2.0*pi*time)) !advance a little bit to compute numerically derivative
  G_minus = G - dt_derivative * tanh(sigma_r*cos(2.0*pi*time)) !reverse a little bit to compute numerically derivative
  
  !------------------------ first lets deal with the angle
  alpha       = -beta * G_rotation(time, sigma_r, sigma_t, phi, beta)
  alpha_minus = -beta * G_rotation(time-dt_derivative, sigma_r, sigma_t, phi, beta)
  alpha_plus  = -beta * G_rotation(time+dt_derivative, sigma_r, sigma_t, phi, beta)
  
  alpha_t  = -beta *( alpha_plus-alpha_minus)/2.0/dt_derivative  
  alpha_tt = -beta *( alpha_minus -2.0*alpha + alpha_plus )/dt_derivative**2
  
  !------------------------------- then with the position  
  C       = C_startup(time, sigma_r, sigma_t, phi, beta)
  C_plus  = C_startup(time+dt_derivative, sigma_r, sigma_t, phi, beta)
  C_minus = C_startup(time-dt_derivative, sigma_r, sigma_t, phi, beta)
  
  ! variables defined by eldrege (center of gravity...)
  Y_cg = y0 ! does not change in time
  X_cg       = C       * 0.5*A0 * G /  G_translation_max
  X_cg_plus  = C_plus  * 0.5*A0 * G_plus /  G_translation_max
  X_cg_minus = C_minus * 0.5*A0 * G_minus /  G_translation_max
  
  x       = X_cg       + L*sin(alpha)
  x_plus  = X_cg_plus  + L*sin(alpha_plus);
  x_minus = X_cg_minus + L*sin(alpha_minus);
  
  y       = Y_cg + L*cos(alpha) ! note Y_cg does not depend on time
  y_plus  = Y_cg + L*cos(alpha_plus);
  y_minus = Y_cg + L*cos(alpha_minus);
  
  
  vx = (x_plus - x_minus)/2.0/dt_derivative 
  vy = (y_plus - y_minus)/2.0/dt_derivative
  
  ax = (x_plus - 2.0*x + x_minus)/dt_derivative**2
  ay = (y_plus - 2.0*y + y_minus)/dt_derivative**2
  
  LeadingEdge(1)=x;
  LeadingEdge(2)=y;
  LeadingEdge(3)=vx;
  LeadingEdge(4)=vy;
  LeadingEdge(5)=ax;
  LeadingEdge(6)=ay;
  

end subroutine

!===================================================================== 

real(kind=pr) function G_translation(time, sigma_r, sigma_t, phi, beta)
  use share_vars
  implicit none
  real (kind=pr), intent(in) :: time,sigma_r,sigma_t,phi,beta
  real (kind=pr) :: dt_real, dt, t, G_translation_max
  real (kind=pr), save :: G_translation_old = 0.0, time_old=0.0
  integer :: it,nt

  !----------------------------------
  !   How does this function work? 
  !   You need to numerically integrate the analytical function, but you already computed a value in the last time step
  !   so the idea is that you shouldn't start at t=0, but you can start at the last time step. then you can choose a very small time step
  !   so that your solution is quite precise, even though you use a silly algorithm. it is much cheaper this way.
  !----------------------------------


  ! the integral is computed between one time level and the other, by the simplest available integration rule. as the time step of the solver may
  ! be relatively large, we subcycle the integration of the G_t function, so that its local time step is about dt=1.0e-5
  dt_real = time-time_old  ! don't worry about the first time step, it's fine
  nt = nint(dt_real/1.0e-5)
  dt = dt_real/real(nt) ! round it, so you reach new time in nt steps
  
 
  G_translation = G_translation_old ! start value is value of previous call
  t=time_old
  
  do it=1,nt    
    G_translation = G_translation + dt * tanh(sigma_r*cos(2.0*pi*t))
    t = t + dt
  enddo

  ! iterate save variables
  time_old = time ! new time level
  G_translation_old = G_translation
  
end function

!=====================================================================

real(kind=pr) function G_rotation(time, sigma_r, sigma_t, phi, beta)
  use share_vars
  implicit none
  real (kind=pr), intent(in) :: time,sigma_r,sigma_t,phi,beta
  real (kind=pr) :: G_rotation_max
  
  G_rotation_max = (exp(2.0*sigma_r)-1.0) / (exp(2.0*sigma_r)+1.0) ! determined analytically
  
  G_rotation = tanh(sigma_r*cos(2.0*pi*time+phi)) / G_rotation_max
  
end function

!=====================================================================

real(kind=pr) function C_startup(time, sigma_r, sigma_t, phi, beta)
  use share_vars
  implicit none
  real (kind=pr), intent(in) :: time,sigma_r,sigma_t,phi,beta
  
  C_startup = ( tanh(8.0*time-2.0)+tanh(2.0) )/(1.0+tanh(2.0))
  
end function
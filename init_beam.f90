subroutine init_beam (beam, pressure_beam)
  use share_vars
  implicit none
  integer :: n
  real (kind=pr), dimension (0:ns-1, 1:6), intent (out) :: beam
  ! Beam indices: 1=beam_x 2=beam_y 3=beam_vx 4=beam_vy 5=theta 6=theta_dot
  real (kind=pr), dimension (0:ns-1), intent (out) :: pressure_beam
  real (kind=pr) :: alpha, alpha_t, alpha_tt,dx,dy
  real (kind=pr), dimension(1:6) :: LeadingEdge !LeadingEdge: x, y, vx, vy, ax, ay (Array)  
  call mouvement(0.0, alpha, alpha_t, alpha_tt, LeadingEdge )

  beam = 0.0
  pressure_beam = 0.00
  dx=xl/real(nx)
  dy=yl/real(ny)



  if ((iMultires.ne.666).and.(iMultires.ne.777)) then 
  write (*,'(" --- Beam width (2*t_beam) covers ",(f4.1)," points")') 2.0*t_beam/max(dy,dx)
  endif

!	n=0
!	beam(n,1) = LeadingEdge(1) 
!	beam(n,2) = LeadingEdge(2) 

   
  do n = 0, ns-1 
! 	beam(n,5) = 0.0*pi*(real(n) / real(ns-1))/180.0
!	beam(n,1) = ds*cos(beam(n,5)) + beam(n-1,1)
!	beam(n,2) = ds*sin(beam(n,5)) + beam(n-1,2)
! 	
! 

 	beam(n,1) = LeadingEdge(1) + real(n) *ds
 	beam(n,2) = LeadingEdge(2) !+ 0.5*(real(n)/real(ns))**2
! beam(n,4) = 1.9*(real(n)/real(ns))**2
 
  enddo
  call integrate_position (0.0, beam) ! to take static angle into account

  
  if (inicond==99) beam = beam_init ! for batch resolution runs
  !----------------------------------------------------------------
  ! set pressure for the beam to the last value of the previous simulation, if required
  !---------------------------------------------------------------- 
  if (inicond==99) pressure_beam=pressure_beam_init


end subroutine init_beam

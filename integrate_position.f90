subroutine integrate_position (time, beam) !attention! parameter time added (for mouvement)
  use share_vars
! this function integrates the incoming beam(:,5) and beam(:,6) (angle and angle_dot) and returns the new positions
! in the same array
  implicit none
  integer :: i,s
  real (kind=pr) :: xrel, yrel, vxrel, vyrel
  real (kind=pr), dimension (0:ns-1, 1:6), intent (inout) ::  beam
  real (kind=pr) :: alpha, alpha_t, alpha_tt 
  real (kind=pr), intent(in) :: time
  real (kind=pr), dimension(1:6) :: LeadingEdge !LeadingEdge: x, y, vx, vy, ax, ay (Array)  
  
  call mouvement(time, alpha, alpha_t, alpha_tt, LeadingEdge )
  
  xrel = 0.0
  yrel = 0.0
  vxrel = 0.0
  vyrel = 0.0
  beam(:,1:4) = 0.0
  
  do s=0,ns-1
    do i=0,s-1
      xrel = xrel + cos(beam(i,5))*ds
      yrel = yrel + sin(beam(i,5))*ds
      vxrel = vxrel - beam(i,6)*sin(beam(i,5))*ds 
      vyrel = vyrel + beam(i,6)*cos(beam(i,5))*ds 
    enddo 
    ! these sums now contain relative pos / vel 
    ! convert them to absolute values:
    beam(s,1) = LeadingEdge(1) + xrel*cos(alpha) - yrel*sin(alpha)
    beam(s,2) = LeadingEdge(2) + yrel*cos(alpha) + xrel*sin(alpha)
    beam(s,3) = LeadingEdge(3) + cos(alpha)*(vxrel - alpha_t*yrel) - sin(alpha)*(vyrel + alpha_t*xrel)
    beam(s,4) = LeadingEdge(4) + cos(alpha)*(vyrel + alpha_t*xrel) + sin(alpha)*(vxrel - alpha_t*yrel)
    xrel = 0.0
    yrel = 0.0
    vxrel = 0.0
    vyrel = 0.0    
  enddo
  

end subroutine integrate_position

module BeamForces
  implicit none
  contains

subroutine GetForces (time, beam, pressure_beam, press, force_pressure )
  use share_vars
  use FieldExport
  use Interpolation
  implicit none
  integer :: n,nx_new,ny_new, upsampling, smoothing
  integer :: i,j,iymin,iymax,ixmin,ixmax,beam_cg_y,beam_cg_x
  real (kind=pr), intent (in)   :: time
  integer                       :: SPL
  real (kind=pr), intent (in)   :: beam (0:ns-1, 1:6)
  real (kind=pr), intent (out)  :: pressure_beam (0:ns-1)
  real (kind=pr), intent (out)  :: force_pressure (1:2)
  real (kind=pr), intent (in)   :: press (0:nx-1, 0:ny-1)
  
  real (kind=pr)                :: press_dx (0:up*nxs-1, 0:up*nys-1)  ,press_dy (0:up*nxs-1, 0:up*nys-1)
  real (kind=pr)                :: press_dxdy (0:up*nxs-1, 0:up*nys-1),press_k (0:up*nxs-1, 0:up*nys-1)
  real (kind=pr)                :: temp (0:up*nxs-1, 0:up*nys-1)
  
  real (kind=pr)                :: pressure_interp (0:up*nxs-1,0:up*nys-1), filter(0:nxs-1,0:nys-1)
  real (kind=pr)                :: xu,yu,xb,yb,x1,y1,x2,y2,delta_smoothing
  real (kind=pr)                :: alpha, alpha_t, alpha_tt, fx, fy,soft_startup,dx,dy
  real (kind=pr), dimension(1:6) :: LeadingEdge !LeadingEdge: x, y, vx, vy, ax, ay (Array)  
        
  call mouvement(time, alpha, alpha_t, alpha_tt, LeadingEdge )

  !------------------------------------------------------------------------------------------------------------------------------
  ! soft_startup: pressure can be slowly activated
  !------------------------------------------------------------------------------------------------------------------------------
  if (time <= T_release) then
      soft_startup = 0.0
  elseif ( ( time >T_release ).and.(time<(T_release + tau)) ) then
      soft_startup =  ((time-T_release)**3)/(-0.5*tau**3)   + 3.*((time-T_release)**2)/tau**2
  else
      soft_startup = 1.0
  endif

  dx=xl/real(nx)
  dy=yl/real(ny)
  
!  smoothing = 20 ! in gridpoints
!  delta_smoothing = real(smoothing)/2.0 !we use this layer for smoothing
  
!  ! roughly estimate center of gravity. is only first order approach, but sufficient
!  beam_cg_x = nint(ds*sum(beam(0:ns-2,1))/length/dx)
!  beam_cg_y = nint(ds*sum(beam(0:ns-2,2))/length/dy)
  
!  ! first lets choose the box without taking care of borders
!  iymin = beam_cg_y - nys/2
!  ixmin = beam_cg_x - nxs/2
!  iymax = iymin + nys - 1 
!  ixmax = ixmin + nxs - 1 
  
!  ! then lets correct possible mistakes
!  if (iymax>ny-1) then
!    iymax = ny-1
!    iymin = iymax - nys + 1
!  endif
!  if (ixmax>nx-1) then
!    ixmax = nx-1
!    ixmin = ixmax - nxs + 1
!  endif
  
!  if (iymin<0) then
!    iymin = 0
!    iymax = nys - 1
!  endif
!  if (ixmin<0) then
!    ixmin = 0
!    ixmax = nxs - 1
!  endif  
  
!  x1=real(ixmin)*dx
!  y1=real(iymin)*dy
!  x2=real(ixmax+1)*dx ! remember in peridoc grid, the last point is ny-1 but the dx still dx=xl/nx
!  y2=real(iymax+1)*dy
  
!  nx_new = ixmax-ixmin 
!  ny_new = iymax-iymin
  
!  write (*,'("time:",es11.4," dx=",es11.4," dy=",es11.4," box: ",i4,":",i4," ",i4,":",i4," BOX: ",4(es11.4,1x), " CG ",i4,i4)'   ) &
!  time,dx,dy,ixmin,ixmax,iymin,iymax, x1,y1,x2,y2,beam_cg_x,beam_cg_y
 
!  call PeriodizeFilter ( filter,delta_smoothing )  
!  filter=1.0
!  call SpectralInterpolation ( filter*press(ixmin:ixmax,iymin:iymax), pressure_interp, press_k )    
  !------------------------------------------------------------------------------------------------------------------
  !globales upsampling
  call SpectralInterpolation ( press, pressure_interp, press_k )    
    x1=0.
    y1=0.
    x2=xl
    y2=yl
  !------------------------------------------------------------------------------------------------------------------
!  pressure_interp=press(ixmin:ixmax,iymin:iymax)


!  call SaveGif(press,'source',1)
!  call SaveGif(pressure_interp,'target',1)
!  stop
 
  !------------------------------------------------------------------------------------------------------------------------------
  !--		Compute interpolated values of the pressure at the beampoints (pressure jump accross the beam)
  !------------------------------------------------------------------------------------------------------------------------------
SPL=iRPM

  if (SPL==1) then
!     call coftxy_big(press,press_k)
     
     call cofdx_big(press_k,temp)
     call cofitxy_big(temp,press_dx)
!     call SaveGif(press_dx,"dx")     
     
     call cofdy_big(press_k,temp)
     call cofitxy_big(temp,press_dy)
!     call SaveGif(press_dy,"dy")
     
     call cofdxdy_big(press_k,temp)
     call cofitxy_big(temp,press_dxdy)
!     call SaveGif(press_dxdy,"dxdy")
  endif
!  stop

  do n=0, ns-1
     !top point
    xu = beam(n,1)  - t_beam * sin(beam(n,5)+alpha)
    yu = beam(n,2)  + t_beam * cos(beam(n,5)+alpha)
    !bottom point
    xb = beam(n,1)  + t_beam * sin(beam(n,5)+alpha)
    yb = beam(n,2)  - t_beam * cos(beam(n,5)+alpha)
    
    !pressure jump across the beam
    if (SPL==1) then
    pressure_beam(n) = BicubicInterpolation (xu, yu, pressure_interp,press_dx,press_dy,press_dxdy, x1,y1,x2,y2) &
                     - BicubicInterpolation (xb, yb, pressure_interp,press_dx,press_dy,press_dxdy, x1,y1,x2,y2)
    elseif (SPL==0) then
    pressure_beam(n) = LinearInterpolation (xu, yu, pressure_interp, x1,y1,x2,y2) &
                     - LinearInterpolation (xb, yb, pressure_interp, x1,y1,x2,y2)                     
    elseif (SPL==2) then
    pressure_beam(n) = SplineInterpolation (xu, yu, pressure_interp, x1,y1,x2,y2) &
                     - SplineInterpolation (xb, yb, pressure_interp, x1,y1,x2,y2)                                          
    endif
  enddo 
  
  !------------------------------------------------------------------------------------------------------------------------------
  !--		startup (slow coupling or values from another simulation)
  !------------------------------------------------------------------------------------------------------------------------------ 
  if (inicond==99) then
    pressure_beam = pressure_beam_init*(1.0-soft_startup)+pressure_beam*soft_startup
  else
    pressure_beam = pressure_beam*soft_startup ! normal, without starting pressure
  endif
  
  !------------------------------------------------------------------------------------------------------------------------------
  !--		Compute Drag and Lift resulting from the pressure on the beam (w/o viscous tensions)
  !------------------------------------------------------------------------------------------------------------------------------

  fx=0.0
  fy=0.0
  !beam segments
  do i=0, ns-2
    fx=fx - sin(beam(i,5)+alpha)*ds*0.5*(pressure_beam(i)+pressure_beam(i+1))
    fy=fy + cos(beam(i,5)+alpha)*ds*0.5*(pressure_beam(i)+pressure_beam(i+1))
  enddo
  !endpoints
  
  fx=fx-pressure_beam(0)*2.0*t_beam*cos(beam(0,5)+alpha)
  fy=fy-pressure_beam(0)*2.0*t_beam*sin(beam(0,5)+alpha)
  
  fx=fx+pressure_beam(ns-1)*2.0*t_beam*cos(beam(ns-1,5)+alpha)
  fy=fy+pressure_beam(ns-1)*2.0*t_beam*sin(beam(ns-1,5)+alpha)
  
  force_pressure(1)=fx
  force_pressure(2)=fy

end subroutine GetForces

!=================================================================================================================================

subroutine PeriodizeFilter(filter,delta_filter)
  use share_vars
  implicit none
  real (kind=pr), intent (inout) :: filter(:,:)
  real (kind=pr), intent (in) :: delta_filter
  integer :: nx_small,ny_small,ix,iy
  real (kind=pr) :: tmp1,tmp2,tmp3,tmp4
  
  nx_small = size(filter,1)
  ny_small = size(filter,2)
  
  filter = 0.0    
  do ix = 1,nx_small
    do iy = 1,ny_small
      call FilterStep (tmp1,real(ix-1)        ,delta_filter,delta_filter) !left
      call FilterStep (tmp2,real(nx_small-ix) ,delta_filter,delta_filter) !right
      call FilterStep (tmp3,real(iy-1)        ,delta_filter,delta_filter) !left
      call FilterStep (tmp4,real(ny_small-iy) ,delta_filter,delta_filter) !right 
          
      filter(ix,iy) = max(tmp1,tmp2,tmp3,tmp4  )
    enddo
  enddo
  
  filter=1.0-filter
end subroutine PeriodizeFilter

!=================================================================================================================================

subroutine FilterStep (f,x,t,h)
  use share_vars
  implicit none
  !-----------------------------------------------------------------
  !-- This subroutine returns the value f of a smooth step function
  !-- The sharp step function would be 1 if x<=t and 0 if x>t
  !-- h is the semi-size of the smoothing area, so
  !-- f is 1 if x<=t-h
  !-- f is 0 if x>t+h
  !-- f is variable (smooth) in between
  !-----------------------------------------------------------------
  real (kind=pr), intent (out) :: f
  real (kind=pr), intent (in)  :: x,t,h
  real (kind=pr) :: a,b,c,d, delta, GradientERF

  if (x<=t-h) then
    f = 1.0
  elseif (((t-h)<x).and.(x<(t+h))) then
    f = 0.5*(1.+cos((x-t+h)*pi/(2.0*h)) )
  else
    f = 0.0
  endif

end subroutine FilterStep


!===================================================================================================================


end module beamforces

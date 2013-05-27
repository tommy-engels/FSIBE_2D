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
  real (kind=pr)                :: press_dx (0:nx-1, 0:ny-1),press_dy (0:nx-1, 0:ny-1)
  real (kind=pr)                :: press_dxdy (0:nx-1, 0:ny-1),press_k (0:nx-1, 0:ny-1)
  real (kind=pr)                :: temp (0:nx-1, 0:ny-1)
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
  
!!  smoothing = 20 ! in gridpoints
!!  delta_smoothing = real(smoothing)/2.0 !we use this layer for smoothing
  
!!  ! roughly estimate center of gravity. is only first order approach, but sufficient
!!  beam_cg_x = nint(ds*sum(beam(0:ns-2,1))/length/dx)
!!  beam_cg_y = nint(ds*sum(beam(0:ns-2,2))/length/dy)
  
!!  ! first lets choose the box without taking care of borders
!!  iymin = beam_cg_y - nys/2
!!  ixmin = beam_cg_x - nxs/2
!!  iymax = iymin + nys - 1 
!!  ixmax = ixmin + nxs - 1 
  
!!  ! then lets correct possible mistakes
!!  if (iymax>ny-1) then
!!    iymax = ny-1
!!    iymin = iymax - nys + 1
!!  endif
!!  if (ixmax>nx-1) then
!!    ixmax = nx-1
!!    ixmin = ixmax - nxs + 1
!!  endif
  
!!  if (iymin<0) then
!!    iymin = 0
!!    iymax = nys - 1
!!  endif
!!  if (ixmin<0) then
!!    ixmin = 0
!!    ixmax = nxs - 1
!!  endif  
  
!!  x1=real(ixmin)*dx
!!  y1=real(iymin)*dy
!!  x2=real(ixmax)*dx
!!  y2=real(iymax)*dy
  
!!  nx_new = ixmax-ixmin 
!!  ny_new = iymax-iymin
  
!!  write (*,'("time:",es11.4," dx=",es11.4," dy=",es11.4," box: ",i4,":",i4," ",i4,":",i4," BOX: ",4(es11.4,1x), " CG ",i4,i4)'   ) &
!!  time,dx,dy,ixmin,ixmax,iymin,iymax, x1,y1,x2,y2,beam_cg_x,beam_cg_y
 
!!!  call PeriodizeFilter ( filter,delta_smoothing )  
!!!  call SpectralInterpolation ( filter*press(ixmin:ixmax,iymin:iymax), pressure_interp )    
!!  pressure_interp=press(ixmin:ixmax,iymin:iymax)


!!  call SaveGif(filter*press(ixmin:ixmax,iymin:iymax),'source')
!!  call SaveGif(pressure_interp,'target')
 
!  !------------------------------------------------------------------------------------------------------------------------------
  !--		Compute interpolated values of the pressure at the beampoints (pressure jump accross the beam)
  !------------------------------------------------------------------------------------------------------------------------------
SPL=iRPM

  if (SPL==1) then
     call coftxy(press,press_k)
     
     call cofdx(press_k,temp)
     call cofitxy(temp,press_dx)
     
     call cofdy(press_k,temp)
     call cofitxy(temp,press_dy)
     
     call cofdxdy(press_k,temp)
     call cofitxy(temp,press_dxdy)
  endif

  do n=0, ns-1
  
    x1=0.
    y1=0.
    x2=xl
    y2=yl

    !top point
    xu = beam(n,1)  - t_beam * sin(beam(n,5)+alpha)
    yu = beam(n,2)  + t_beam * cos(beam(n,5)+alpha)
    !bottom point
    xb = beam(n,1)  + t_beam * sin(beam(n,5)+alpha)
    yb = beam(n,2)  - t_beam * cos(beam(n,5)+alpha)
    
    !pressure jump across the beam
    if (SPL==1) then
    pressure_beam(n) = BicubicInterpolation (xu, yu, press,press_dx,press_dy,press_dxdy, x1,y1,x2,y2) &
                     - BicubicInterpolation (xb, yb, press,press_dx,press_dy,press_dxdy, x1,y1,x2,y2)
    elseif (SPL==0) then
    pressure_beam(n) = LinearInterpolation (xu, yu, press, x1,y1,x2,y2) &
                     - LinearInterpolation (xb, yb, press, x1,y1,x2,y2)                     
    elseif (SPL==2) then
    pressure_beam(n) = SplineInterpolation (xu, yu, press, x1,y1,x2,y2) &
                     - SplineInterpolation (xb, yb, press, x1,y1,x2,y2)                                          
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


end module beamforces

module BeamForces
  implicit none
  contains
  
  
subroutine GetForces (time, beam, pressure_beam, press, force_pressure )
  use share_vars
  use FieldExport
  implicit none
  real (kind=pr), intent (in)   :: time
  integer                       :: SPL,nxl,nyl,a,b,c
  real (kind=pr), intent (in)   :: beam (0:ns-1, 1:6)
  real (kind=pr), intent (out)  :: pressure_beam (0:ns-1)
  real (kind=pr), intent (out)  :: force_pressure (1:2)
  real (kind=pr), intent (in)   :: press (0:nx-1, 0:ny-1)
  
  real (kind=pr), dimension(0:ns-1)  :: p1,p2,p3,p4,p5,p6,p7,p66
  
!  ! zuerst mit box und upsampling
!  call GetForcesBox (time, beam, p1, press, force_pressure,0  )! linear
!  call GetForcesBox (time, beam, p2, press, force_pressure,1  )! bicubic
  
!  !dann mit box, aber ohne upsamling
!  call GetForcesBoxNoUp (time, beam, p3, press, force_pressure,0  )! linear
!  call GetForcesBoxNoUp (time, beam, p4, press, force_pressure,1  )! bicubic  
  
!  !dann ohne box, ohne alles
!  call GetForcesOhneAlles (time, beam, p5, press, force_pressure,0  )! linear
!  call GetForcesOhneAlles (time, beam, p6, press, force_pressure,1  )! bicubic  
  call GetForcesOhneAlles (time, beam, pressure_beam, press, force_pressure, 2  )! bicubic, exact derivs
  
  ! dann noch uralte versin
!  call GetForcesAlteVersion (time, beam, p7, press, force_pressure)
  
!  a=ns-1
!  b=3*ns/4
!  c=ns/2
!  open  (10, file=trim(simulation_name)//"pressures.aaa", status = 'unknown', access = 'append')
!  write (10, '(9(es15.8,1x))') time,p1(a), p2(a), p3(a), p4(a), p5(a), p6(a),p66(a), p7(a)
!  close (10)
  
!  open  (10, file=trim(simulation_name)//"pressures.bbb", status = 'unknown', access = 'append')
!  write (10, '(9(es15.8,1x))') time,p1(b), p2(b), p3(b), p4(b), p5(b), p6(b),p66(b), p7(b)
!  close (10)
  
!  open  (10, file=trim(simulation_name)//"pressures.ccc", status = 'unknown', access = 'append')
!  write (10, '(9(es15.8,1x))') time,p1(c), p2(c), p3(c), p4(c), p5(c), p6(c),p66(c), p7(c)
!  close (10)
  
end subroutine GetForces



subroutine VelocityResiduum (time,u, beam)
! This subroutine computes the surface integral of the velocity for a turek test (cylinder with beam). The integral should actually be zero
! but due to the volume penalization method it is not exactly (see discussion with romain dirichilet vs navier slip boundary conditions

  use share_vars
  use FieldExport
  use Interpolation
  implicit none
  integer, parameter 		 :: nc=512 ! points on the cylinder surface
  integer			 :: i,n
  real (kind=pr), intent (in)   :: time
  real (kind=pr)	        :: dr, e_cyl_norm, e_cyl_tang, e_bea_norm, e_bea_tang, alpha=0.0, alpha_t,alpha_tt,dx,dy, res
  real (kind=pr), intent (in)   :: beam (0:ns-1, 1:6) 
  
  real (kind=pr), dimension(0:nc)   :: phi, ux_c, uy_c, xc, yc !cylinder
  real (kind=pr), dimension(0:ns-1) :: yt,yb,xt,xb, ux_bb, uy_bb, ux_bt, uy_bt ! beam  
  real (kind=pr), dimension (0:nx-1, 0:ny-1, 1:2), intent (inout) :: u
  real (kind=pr), dimension (0:nx-1, 0:ny-1) :: work1, work2, work3, u_dx, u_dy, u_dxdy
  
!   u(:,:,1) = 1.0
!   u(:,:,2) = 1.0
  
  ! the cylinder
  do i = 0, nc
    phi(i) = 2.d0*pi*real(i) / real(nc)
    xc (i) = R_cylinder*sin(phi(i)) + x0 - R_cylinder ! note the cylinder is not at (x0|y0)
    yc (i) = R_cylinder*cos(phi(i)) + y0
  enddo
  
  ! the beam (top and bottom side)
!  call mouvement(time, alpha, alpha_t, alpha_tt, LeadingEdge )
  do n=0, ns-1
    !top point
    xt(n) = beam(n,1) 
    yt(n) = beam(n,2)  + t_beam 
    !bottom point
    xb(n) = beam(n,1)  
    yb(n) = beam(n,2)  - t_beam
  enddo  
  
  dx = xl/real(nx)
  dy = yl/real(ny)
  dr = 2.0*pi*R_cylinder/real(nc)
  
  !---------------------------------------------------------------
  ! first we interpolate the ux values on the surface
  !---------------------------------------------------------------
  call coftxy	(u(:,:,1),work1)  
  call cofdx	(work1,work2)
  call cofitxy	(work2,u_dx)  
  call cofdy	(work1,work2)
  call cofitxy	(work2,u_dy)  
  call cofdy	(work1,work2)
  call cofitxy	(work2,u_dxdy)

  do i = 0, nc
    ux_c (i) = BicubicInterpolation (xc(i), yc(i), u(:,:,1) ,u_dx,u_dy,u_dxdy, 0.0,0.0,xl-dx,yl-dy)
  enddo
  
  do i = 0, ns-1
    ux_bt (i) = BicubicInterpolation (xt(i), yt(i), u(:,:,1) ,u_dx,u_dy,u_dxdy, 0.0,0.0,xl-dx,yl-dy)
    ux_bb (i) = BicubicInterpolation (xb(i), yb(i), u(:,:,1) ,u_dx,u_dy,u_dxdy, 0.0,0.0,xl-dx,yl-dy)
  enddo
  
  !---------------------------------------------------------------
  ! then UY
  !---------------------------------------------------------------
  
  call coftxy	(u(:,:,2),work1)  
  call cofdx	(work1,work2)
  call cofitxy	(work2,u_dx)  
  call cofdy	(work1,work2)
  call cofitxy	(work2,u_dy)  
  call cofdy	(work1,work2)
  call cofitxy	(work2,u_dxdy)

  do i = 0, nc
    uy_c (i) = BicubicInterpolation (xc(i), yc(i), u(:,:,2),u_dx,u_dy,u_dxdy, 0.0,0.0,xl-dx,yl-dy)
  enddo
  
  do i = 0, ns-1
    uy_bt (i) = BicubicInterpolation (xt(i), yt(i), u(:,:,2) ,u_dx,u_dy,u_dxdy, 0.0,0.0,xl-dx,yl-dy)
    uy_bb (i) = BicubicInterpolation (xb(i), yb(i), u(:,:,2) ,u_dx,u_dy,u_dxdy, 0.0,0.0,xl-dx,yl-dy)
  enddo

  !---------------------------------------------------------------
  ! now lets integrate (on the surface)
  !---------------------------------------------------------------
  
  e_cyl_norm = dr*sum( ( (cos(phi)*ux_c) + (sin(phi)*uy_c) )**2 )
  e_cyl_tang = dr*sum( ( (sin(phi)*ux_c) - (cos(phi)*uy_c) )**2 )
  
  e_bea_norm = ds*sum(uy_bt**2) + ds*sum(uy_bb**2 )  !for the moment, only straight beams
  e_bea_tang = ds*sum(ux_bt**2) + ds*sum(ux_bb**2 )
  
  res = dr*sum( sqrt (ux_c**2 + uy_c**2) ) + ds* sum( sqrt( ux_bt**2 + uy_bt**2)) + ds* sum( sqrt( ux_bb**2 + uy_bb**2))
  
  open (14, file = trim(dir_name)//'/'//trim(simulation_name)//'residuum', status = 'unknown', access = 'append') ! Append output data file
  write (14, '(6(es15.8,1x))') time, e_cyl_norm, e_cyl_tang, e_bea_norm, e_bea_tang, res
!   write (14,*) '--- cylinder'
!   write (14,'(513(es11.4,1x))') phi
!   write (14,'(513(es11.4,1x))') xc
!   write (14,'(513(es11.4,1x))') yc
!   write (14,'(513(es11.4,1x))') ux_c
!   write (14,'(513(es11.4,1x))') uy_c
!   write (14,'(513(es11.4,1x))') ( (cos(phi)*ux_c) + (sin(phi)*uy_c) )**2
!   write (14,'(513(es11.4,1x))') ( (sin(phi)*ux_c) - (cos(phi)*uy_c) )**2  
!   write (14,*) '---beam'
!   write (14,'("xt",256(es11.4,1x))') xt
!   write (14,'("yt",256(es11.4,1x))') yt
!   write (14,'(256(es11.4,1x))') ( uy_bt**2 + uy_bb**2 ) 
!   write (14,'(256(es11.4,1x))') ( ux_bt**2 + ux_bb**2 )
!   write (14,'(256(es11.4,1x))') ux_bt
!   write (14,'(256(es11.4,1x))') uy_bt
!   write (14,'(256(es11.4,1x))') ux_bb
!   write (14,'(256(es11.4,1x))') uy_bb 
  
  close (14) ! Close the file to protect data
  
!   stop

end subroutine VelocityResiduum

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

subroutine GetForcesBoxNoUp (time, beam, pressure_beam, press, force_pressure,SPL )
  use share_vars
  use FieldExport
  use Interpolation
  implicit none
  integer :: n,nx_new,ny_new, upsampling, smoothing,safety,nxs,nys
  integer :: i,j,iymin,iymax,ixmin,ixmax,beam_cg_y,beam_cg_x
  real (kind=pr), intent (in)   :: time
  integer       , intent (in)   :: SPL
  integer                       :: nxl,nyl
  real (kind=pr), intent (in)   :: beam (0:ns-1, 1:6)
  real (kind=pr), intent (out)  :: pressure_beam (0:ns-1)
  real (kind=pr), intent (out)  :: force_pressure (1:2)
  real (kind=pr), intent (in)   :: press (0:nx-1, 0:ny-1)
  
  real (kind=pr), allocatable   :: press_smaller(:,:)
  
  real (kind=pr), allocatable   :: pressure_interp (:,:)
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
  !------------------------------------------------------------------------------------------------------------------------------
  ! get a smaller box containing the beam region
  !------------------------------------------------------------------------------------------------------------------------------
  !this is the same everywhere in the code:
  dx=xl/real(nx)
  dy=yl/real(ny)

  safety=20 !in gridpoints
  iymin = max (0   , int((minval(beam(:,2))-2.0*t_beam)/dy)-safety )
  ixmin = max (0   , int((minval(beam(:,1))-2.0*t_beam)/dx)-safety )
  iymax = min (ny-1, int((maxval(beam(:,2))+2.0*t_beam)/dy)+safety )
  ixmax = min (nx-1, int((maxval(beam(:,1))+2.0*t_beam)/dx)+safety )

  nxs = ixmax-ixmin + 1
  nys = iymax-iymin + 1 
  
  x1=real(ixmin)*dx
  y1=real(iymin)*dy
  x2=real(ixmax)*dx 
  y2=real(iymax)*dy 
  
  allocate ( press_smaller(0:nxs-1,0:nys-1) )
  press_smaller = press(ixmin:ixmax,iymin:iymax)
 
  !------------------------------------------------------------------------------------------------------------------------------
  !--		Compute interpolated values of the pressure at the beampoints (pressure jump accross the beam)
  !------------------------------------------------------------------------------------------------------------------------------
!  SPL=iRPM

  do n=0, ns-1
    !top point
    xu = beam(n,1)  - t_beam * sin(beam(n,5)+alpha)
    yu = beam(n,2)  + t_beam * cos(beam(n,5)+alpha)
    !bottom point
    xb = beam(n,1)  + t_beam * sin(beam(n,5)+alpha)
    yb = beam(n,2)  - t_beam * cos(beam(n,5)+alpha)

    if (SPL==0) then
    pressure_beam(n) = LinearInterpolation (xu, yu, press_smaller, x1,y1,x2,y2) &
                     - LinearInterpolation (xb, yb, press_smaller, x1,y1,x2,y2)  
    elseif(SPL==1) then
    pressure_beam(n) = BicubicInterpolationFD (xu, yu, press_smaller, x1,y1,x2,y2) &
                     - BicubicInterpolationFD (xb, yb, press_smaller, x1,y1,x2,y2)                           
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

end subroutine GetForcesBoxNoUp
  
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

subroutine GetForcesOhneAlles(time, beam, pressure_beam, press, force_pressure, SPL )
  use share_vars
  use FieldExport
  use Interpolation
  implicit none
  integer :: n,nx_new,ny_new, upsampling, smoothing,safety,nxs,nys
  integer :: i,j,iymin,iymax,ixmin,ixmax,beam_cg_y,beam_cg_x
  real (kind=pr), intent (in)   :: time
  integer       , intent (in)   :: SPL
  integer                       :: nxl,nyl
  real (kind=pr), intent (in)   :: beam (0:ns-1, 1:6)
  real (kind=pr), intent (out)  :: pressure_beam (0:ns-1)
  real (kind=pr), intent (out)  :: force_pressure (1:2)
  real (kind=pr), intent (in)   :: press (0:nx-1, 0:ny-1)
  real (kind=pr)                :: xu,yu,xb,yb,x1,y1,x2,y2,delta_smoothing
  real (kind=pr)                :: alpha, alpha_t, alpha_tt, fx, fy,soft_startup,dx,dy
  real (kind=pr), dimension(1:6) :: LeadingEdge !LeadingEdge: x, y, vx, vy, ax, ay (Array)  
  real (kind=pr)                :: press_dx (0:nx-1, 0:ny-1),press_dy (0:nx-1, 0:ny-1)
  real (kind=pr)                :: press_dxdy (0:nx-1, 0:ny-1),press_k (0:nx-1, 0:ny-1)
  real (kind=pr)                :: temp (0:nx-1, 0:ny-1)
        
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
  !------------------------------------------------------------------------------------------------------------------------------
  !--		Compute interpolated values of the pressure at the beampoints (pressure jump accross the beam)
  !------------------------------------------------------------------------------------------------------------------------------

  if (SPL==2)then
     call coftxy(press,press_k)
 
     call cofdx(press_k,temp)
     call cofitxy(temp,press_dx)
     
     call cofdy(press_k,temp)
     call cofitxy(temp,press_dy)
     
     call cofdxdy(press_k,temp)
     call cofitxy(temp,press_dxdy)
  endif

  
  do n=0, ns-1
    !top point
    xu = beam(n,1)  - t_beam * sin(beam(n,5)+alpha)
    yu = beam(n,2)  + t_beam * cos(beam(n,5)+alpha)
    !bottom point
    xb = beam(n,1)  + t_beam * sin(beam(n,5)+alpha)
    yb = beam(n,2)  - t_beam * cos(beam(n,5)+alpha)

    if (SPL==0) then
    pressure_beam(n) = LinearInterpolation (xu, yu, press, 0.0,0.0,xl-dx,yl-dy) &
                     - LinearInterpolation (xb, yb, press, 0.0,0.0,xl-dx,yl-dy)  
    elseif (SPL==1) then
    pressure_beam(n) = BicubicInterpolationFD (xu, yu, press, 0.0,0.0,xl-dx,yl-dy) &
                     - BicubicInterpolationFD (xb, yb, press, 0.0,0.0,xl-dx,yl-dy)                           
    elseif (SPL==2) then       
    pressure_beam(n) = BicubicInterpolation (xu, yu, press,press_dx,press_dy,press_dxdy, 0.0, 0.0, xl-dx,yl-dy) &
                     - BicubicInterpolation (xb, yb, press,press_dx,press_dy,press_dxdy, 0.0, 0.0, xl-dx,yl-dy)                     
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
  
end subroutine GetForcesOhneAlles
  
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  
subroutine GetForcesAlteVersion (time, beam, pressure_beam, press, force_pressure )

!  Routine to get pressure forces on the beam. solves poisson for pressure.
!  "pressure_beam" is the pressure on the beam "beam" in field "vortk"

  use share_vars
  use FieldExport
  implicit none
  integer :: n
  integer :: i,j
  real (kind=pr) :: x, y,x_1,y_1,x_2,y_2, p_up,p_lo,dx, dy, R1,R2, soft_startup
  real (kind=pr), dimension (0:ns-1, 1:6), intent (in) :: beam
  real (kind=pr), dimension (0:ns-1), intent (out) :: pressure_beam
  real (kind=pr), dimension (1:2), intent (out) :: force_pressure
  real (kind=pr), dimension (0:nx-1, 0:ny-1), intent (in) :: press
  real (kind=pr), intent (in) :: time
  character(len=4) :: name
  character(len=11) :: name_time
  integer, save :: n_snapshot = 1
  real (kind=pr) :: alpha, alpha_t, alpha_tt, fx, fy
  real (kind=pr), dimension(1:6) :: LeadingEdge !LeadingEdge: x, y, vx, vy, ax, ay (Array)  
        
  call mouvement(time, alpha, alpha_t, alpha_tt, LeadingEdge )

 
  !------------------------------------------------------------------------------------------------------------------------------
  !--		Compute interpolated values of the pressure at the beampoints (pressure jumpo accross the beam)
  !------------------------------------------------------------------------------------------------------------------------------
  dx=xl/real(nx)
  dy=yl/real(ny)

  do n=0, ns-1
    !upper point
      x = beam(n,1)  - t_beam * sin(beam(n,5)+alpha)
      y = beam(n,2)  + t_beam * cos(beam(n,5)+alpha)
      i=int(x/dx)
      j=int(y/dy)
      x_1= real(i)*dx
      y_1= real(j)*dy
      x_2= dx*real(i+1)
      y_2= dy*real(j+1)
      R1 = (x_2-x)*press(i,j)/dx   + (x-x_1)*press(i+1,j)/dx
      R2 = (x_2-x)*press(i,j+1)/dx + (x-x_1)*press(i+1,j+1)/dx
      p_up = (y_2-y)*R1/dy + (y-y_1)*R2/dy
    !lower point
      x = beam(n,1)  + t_beam * sin(beam(n,5)+alpha)
      y = beam(n,2)  - t_beam * cos(beam(n,5)+alpha)
      i=int(x/dx)
      j=int(y/dy)
      x_1= real(i)*dx
      y_1= real(j)*dy
      x_2= dx*real(i+1)
      y_2= dy*real(j+1)
      R1 = (x_2-x)*press(i,j)/dx   + (x-x_1)*press(i+1,j)/dx
      R2 = (x_2-x)*press(i,j+1)/dx + (x-x_1)*press(i+1,j+1)/dx
      p_lo = (y_2-y)*R1/dy + (y-y_1)*R2/dy
  !pressure jump across the beam
    pressure_beam(n) = p_up - p_lo
 
   enddo !end of loop over beampoints
end subroutine GetForcesAlteVersion


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
subroutine GetForcesBox (time, beam, pressure_beam, press, force_pressure,SPL )
  use share_vars
  use FieldExport
  use Interpolation
  implicit none
  integer :: n,nx_new,ny_new, upsampling, smoothing,safety,nxs,nys
  integer :: i,j,iymin,iymax,ixmin,ixmax,beam_cg_y,beam_cg_x
  real (kind=pr), intent (in)   :: time
  integer       , intent (in)   :: SPL
  integer                       :: nxl,nyl
  real (kind=pr), intent (in)   :: beam (0:ns-1, 1:6)
  real (kind=pr), intent (out)  :: pressure_beam (0:ns-1)
  real (kind=pr), intent (out)  :: force_pressure (1:2)
  real (kind=pr), intent (in)   :: press (0:nx-1, 0:ny-1)  
  real (kind=pr), allocatable   :: level1(:,:),level2(:,:),level3(:,:),level4(:,:),level5(:,:),press_smaller(:,:)  
  real (kind=pr), allocatable   :: pressure_interp (:,:)
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
  !------------------------------------------------------------------------------------------------------------------------------
  ! get a smaller box containing the beam region
  !------------------------------------------------------------------------------------------------------------------------------
  !this is the same everywhere in the code:
  dx=xl/real(nx)
  dy=yl/real(ny)

  safety=20 !in gridpoints
  iymin = max (0   , int((minval(beam(:,2))-2.0*t_beam)/dy)-safety )
  ixmin = max (0   , int((minval(beam(:,1))-2.0*t_beam)/dx)-safety )
  iymax = min (ny-1, int((maxval(beam(:,2))+2.0*t_beam)/dy)+safety )
  ixmax = min (nx-1, int((maxval(beam(:,1))+2.0*t_beam)/dx)+safety )

  nxs = ixmax-ixmin + 1
  nys = iymax-iymin + 1 
  
  x1=real(ixmin)*dx
  y1=real(iymin)*dy
  x2=real(ixmax)*dx 
  y2=real(iymax)*dy 
  
  allocate ( press_smaller(0:nxs-1,0:nys-1) )
  press_smaller = press(ixmin:ixmax,iymin:iymax)
 
  !------------------------------------------------------------------------------------------------------------------------------
  ! sample it up using Deslauriers-Dubuc interpolation
  !------------------------------------------------------------------------------------------------------------------------------
  nxl = 2*nxs-1; nyl = 2*nys-1 
  allocate ( level1(0:nxl-1,0:nyl-1) )
  call PolynomUpsampling(press_smaller,level1)  
  ! runde 2
  nxl=2*nxl-1; nyl=2*nyl-1
  allocate ( level2(0:nxl-1,0:nyl-1) )
  call PolynomUpsampling(level1,level2)
  ! runde 3
  nxl=2*nxl-1; nyl=2*nyl-1
  allocate ( level3(0:nxl-1,0:nyl-1) )
  call PolynomUpsampling(level2,level3)
  ! runde 4
  nxl=2*nxl-1; nyl=2*nyl-1
  allocate ( level4(0:nxl-1,0:nyl-1) )
  call PolynomUpsampling(level3,level4)
    
  !------------------------------------------------------------------------------------------------------------------------------
  !--		Compute interpolated values of the pressure at the beampoints (pressure jump accross the beam)
  !------------------------------------------------------------------------------------------------------------------------------

  do n=0, ns-1
    !top point
    xu = beam(n,1)  - t_beam * sin(beam(n,5)+alpha)
    yu = beam(n,2)  + t_beam * cos(beam(n,5)+alpha)
    !bottom point
    xb = beam(n,1)  + t_beam * sin(beam(n,5)+alpha)
    yb = beam(n,2)  - t_beam * cos(beam(n,5)+alpha)

    if (SPL==0) then
    pressure_beam(n) = LinearInterpolation (xu, yu, level4, x1,y1,x2,y2) &
                     - LinearInterpolation (xb, yb, level4, x1,y1,x2,y2)  
    elseif(SPL==1)                 then
    pressure_beam(n) = BicubicInterpolationFD (xu, yu, level4, x1,y1,x2,y2) &
                     - BicubicInterpolationFD (xb, yb, level4, x1,y1,x2,y2)                           
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

end subroutine GetForcesBox

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

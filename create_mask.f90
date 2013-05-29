!==============================================================================================================================================
subroutine create_mask ( time, beam )
  use share_vars
  use FieldExport
  implicit none
  integer :: i, j, iymin, iymax,ix,iy
  real (kind=pr) :: mask_cutoff, dx, dy, a, b, volume, ds_inter, N, h_star, temp, alpha1,alpha2,gamma,ConvertAngle, sponge_size, sponge_size_star, tmp,x
  real (kind=pr) :: y_chan, H_effective, Mean_ux, Mean_uy
  integer :: iixmax, iixmin, iiymax, iiymin, n_inter, n_subpoints
  real (kind=pr), intent (in) :: time
! Beam indices: 1=beam_x 2=beam_y 3=beam_vx 4=beam_vy 5=theta 6=theta_dot
  real (kind=pr), dimension (0:ns-1, 1:6), intent (in) :: beam
  real (kind=pr) :: alpha, alpha_t, alpha_tt 
  real (kind=pr), dimension(1:6) :: LeadingEdge !LeadingEdge: x, y, vx, vy, ax, ay (Array)  
  call mouvement(time, alpha, alpha_t, alpha_tt, LeadingEdge )
  
  !--------------------------------------------------------------------------------
  !		Initialization
  !--------------------------------------------------------------------------------
  N=N_smooth
  dx=xl/real(nx)
  dy=yl/real(ny)

  ! check if everthing seems to be ok
  if ((maxval(beam(:,1))>xl).or.(minval(beam(:,1))<0.0).or.(maxval(beam(:,2))>yl).or.(minval(beam(:,2))<0.0)) then
    write (*,*) "fatal in mask: out of domain"
    stop
  endif

  h_star=h_channel+N*dy ! channel walls are always in y-direction

  mask = 0.0 !initialize masks as 0.0
  maskvx = 0.0
  maskvy = 0.0
  
 
  !--------------------------------------------------------------------------------
  !		Draw the walls
  !--------------------------------------------------------------------------------
  if (iWalls>0) then
      !-- bottom wall
      iymin = 0
      iymax = nint( h_star/ dy )
      !$omp parallel do private(j,temp)
      do j=iymin, iymax
        call SmoothStep(temp, real(j)*dy, h_channel, N*dy )
        mask(:,j)= temp
      enddo
      !$omp end parallel do
      
      !-- top wall
      iymin = nint( (yl-h_star)/dy)
      iymax = ny-1
      !$omp parallel do private(j,temp)
      do j=iymax,iymin,-1
        call SmoothStep(temp,abs( real(j)*dy-yl), h_channel, N*dy )
        mask(:,j)= temp
      enddo
      !$omp end parallel do

      !------------------------------------
      ! modified version (sharp walls!)
      !------------------------------------
!       iymin = nint ( h_channel / dy ) + 1 ! note we set the actual BC to one
!       iymax = nint ( (yl-h_channel) / dy ) - 1 ! note we set the actual BC to one
!       
!       mask = 1.0      
!       mask(:,iymin:iymax) = 0.0
      
      
  endif

  !--------------------------------------------------------------------------------
  !		Draw the cylinder
  !--------------------------------------------------------------------------------
  if ((iCylinder==1).and.(iFSI.ne.8)) then  !usually the case for a cylinder at the leading edge
    call DrawADot( x0-R_cylinder, y0, 0.0, 0.0, R_cylinder, N)
  endif


  !--------------------------------------------------------------------------------
  !		Draw the beam
  !--------------------------------------------------------------------------------
  if (iBeam>0) then
    !---- Draw beam segments
    do i=0, ns-2
      call DrawBeamSegment( beam(i,1),  beam(i,2),&
			    beam(i+1,1),beam(i+1,2),&
			    beam(i,3)  ,beam(i,4),&
			    beam(i+1,3),beam(i+1,4),&
			    t_beam, N )
    enddo
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (Mask_end==1) then
      call DrawSharpEnd( beam(ns-1,1:2), beam(ns-1,5)+alpha, beam(ns-1,3:4),t_beam, N )
    endif
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
    
        
    !---- Fill the gaps with tiny circle segments
    do i=1,ns-2
      alpha1 = (atan2( (beam(i,2)-beam(i-1,2)),(beam(i,1)-beam(i-1,1)) ))
      alpha2 = (atan2( (beam(i+1,2)-beam(i,2)),(beam(i+1,1)-beam(i,1)) ))
      call DrawHinge( beam(i,1),beam(i,2),alpha1,alpha2, beam(i,3), beam(i,4), t_beam, N )
    enddo
    !---- Leading edge endpoint
    call DrawADot( beam(0,1),beam(0,2), beam(0,3), beam(0,4), t_beam, N )
    !---- Trailing edge endpoint
    if (Mask_end==0) then
      call DrawADot( beam(ns-1,1),beam(ns-1,2), beam(ns-1,3), beam(ns-1,4), t_beam, N)
    endif

  endif

  !-------------------------------------------------------------------------------------------
  !		Velocity sponge
  !-------------------------------------------------------------------------------------------
  ! in the current version, the "left" side of the sponge, the actual outflow, is not interesting.
  ! we therefore use a very smooth transient function, 25% of the sponge thickness is used for the
  ! transition. this helps reducing Gibbs perturbations owe to the sponge. Astonishingly, the 25% solution
  ! appears to be the best; a linear function is much worse, a sin over half the sponge thickness also.
  ! let's keep it like this.

    if (iSponge==4) then
      sponge_size_star = SpongeSize ! what you set is the effective sponge size
      sponge_size = SpongeSize - 2.0*N*dx   !"core" size of the sponge
      H_effective = yl-2.0*h_channel

      !-------------------------
      !--left half of the sponge
      !-------------------------
      iixmax=nint ((0.5*sponge_size_star)/dx)
      iixmin=0
      !$omp parallel do private(ix,iy,y_chan,tmp,x)
      do ix=iixmin, iixmax
      do iy=0,ny-1
        y_chan = (real(iy)*dy-h_channel)
        x = 0.75*sponge_size_star - real(ix)*dx
!        call SmoothStep (tmp, x, sponge_size*0.5, N*real(dx) ) !for a transition based on the smoothing layer thickness
        call SmoothStep (tmp, x, sponge_size*0.5, 0.25*sponge_size )
        if (mask(ix,iy)<tmp) mask(ix,iy)=tmp
        ! set parabolic velocity profile
        if ((y_chan>0.0).and.(y_chan<H_effective)) then
          ! attention: the maximum velocity here is always 1, but the mean velocity is smaller!
          maskvx(ix,iy) =  1.5*Mean_ux(time)*y_chan*(H_effective-y_chan)/((0.5*H_effective)**2)
        else
          maskvx(ix,iy) =  0.0
        endif
      enddo
      enddo
      !$omp end parallel do
      
      !-------------------------
      !--right half of the sponge
      !-------------------------
      iixmin=nint ((0.5*sponge_size_star)/dx)
      iixmax=nint ((1.0*sponge_size_star)/dx)
      !$omp parallel do private(ix,iy,y_chan,tmp,x)
      do ix=iixmin, iixmax
      do iy=0,ny-1
        y_chan = (real(iy)*dy-h_channel)
        x = real(ix)*dx - 0.5*sponge_size_star
        call SmoothStep (tmp, x, sponge_size*0.5, N*real(dx) )
        if (mask(ix,iy)<tmp) mask(ix,iy)=tmp
        ! set parabolic velocity profile
        if ((y_chan>0.0).and.(y_chan<H_effective)) then
          ! attention: the velocity here is always 1, but the mean velocity is smaller!
	  maskvx(ix,iy) =  1.5*Mean_ux(time)*y_chan*(H_effective-y_chan)/((0.5*H_effective)**2)
        else
          maskvx(ix,iy) =  0.0
        endif
      enddo
      enddo
      !$omp end parallel do
      
      !----------------------------
      ! modified version, sharp mask function     
      !----------------------------
!       iixmin = 0
!       iixmax = nint (SpongeSize/dx)
!       H_effective = yl-2.0*h_channel
!       
!       !$omp parallel do private (iy,y_chan)
!       do iy=0,ny-1
!         y_chan = (real(iy)*dy-h_channel)        
!         ! set parabolic velocity profile
!         if ((y_chan>0.0).and.(y_chan<H_effective)) then
! 	  maskvx( iixmin:iixmax, iy ) =  1.5*Mean_ux(time)*y_chan*(H_effective-y_chan)/((0.5*H_effective)**2)
!         else
!           maskvx( iixmin:iixmax, iy ) =  0.0
!         endif        
! 	mask ( iixmin:iixmax, iy ) = 1.0
!       enddo
!       !$omp end parallel do
      
   endif
   
  !-------------------------------------------------------------------------------------------
  !		Divide by penalization parameter
  !-------------------------------------------------------------------------------------------
  
  !$omp parallel do private(iy)
  do iy=0,ny-1
     mask(:,iy) = mask(:,iy) / eps
  enddo
  !$omp end parallel do
  


end subroutine create_mask


!==============================================================================================================

subroutine create_sponge_mask()
  use share_vars
  use FieldExport
  implicit none
  real (kind=pr) :: dy,dx,sponge_begin, epsilon_sponge
  real (kind=pr) :: sponge_size, tmp,x, N, sponge_size_star,tmp1,tmp2,tmp3,tmp4
  integer :: ix, iixmin, iixmax, iy

  epsilon_sponge=eps_sponge
  N=N_smooth !smoothing layer thickness
  dx=xl/real(nx)
  dy=yl/real(ny)
  mask_sponge=0.0

  if (iSponge==1) then        !sponge in outflow
      sponge_size_star = SpongeSize ! what you set is the effective sponge size
      sponge_size = sponge_size_star - 2.0*N*dx
      sponge_begin = xl - sponge_size_star
      ! left half of the sponge
      iixmax=nint ((xl-0.5*sponge_size_star)/dx)
      iixmin=nint (sponge_begin/dx)
      do ix=iixmin, iixmax
	x = xl-0.5*sponge_size_star - real(ix)*dx
	call SmoothStep (tmp, x, sponge_size*0.5, N*real(dx) )
	mask_sponge(ix,:)=tmp
      enddo

      !right half of the sponge
      iixmin=nint ((xl-0.5*sponge_size_star)/dx)
      iixmax=nx-1
      do ix=iixmin, iixmax
	x = real(ix)*dx - xl +0.5*sponge_size_star-1.0*dx
	call SmoothStep (tmp, x, sponge_size*0.4, N*real(dx) )
	mask_sponge(ix,:)=tmp
      enddo
  elseif (iSponge==2) then ! sponge in inflow
      sponge_size_star = SpongeSize ! what you set is the effective sponge size
      sponge_size = sponge_size_star - 2.0*N*dx
      sponge_begin = xl - sponge_size_star

      ! left half of the sponge
      iixmax=nint ((0.5*sponge_size_star)/dx)
      iixmin=0
      do ix=iixmin, iixmax
        x = 0.5*sponge_size_star - real(ix)*dx
        call SmoothStep (tmp, x, sponge_size*0.5, N*real(dx) )
        mask_sponge(ix,:)=tmp
      enddo

      !right half of the sponge
      iixmin=nint ((0.5*sponge_size_star)/dx)
      iixmax=nint ((1.0*sponge_size_star)/dx)
      do ix=iixmin, iixmax
        x = real(ix)*dx - 0.5*sponge_size_star
        call SmoothStep (tmp, x, sponge_size*0.5, N*real(dx) )
        mask_sponge(ix,:)=tmp
      enddo

      do ix=0,iixmax
      do iy=0,ny-1
        if (mask(ix,iy)*eps>0.0) mask_sponge(ix,iy)= mask_sponge(ix,iy)*(1.0-mask(ix,iy)*eps)
      enddo
      enddo
  elseif (iSponge==3) then
      ! sponge around the domain
      sponge_size_star = SpongeSize ! what you set is the effective sponge size
      sponge_size = sponge_size_star - N*dx
      do ix=0, nx-1
      do iy=0, ny-1
          call SmoothStep (tmp1, real(ix)*dx   , sponge_size, N*real(dx) )
          call SmoothStep (tmp2, xl-real(ix)*dx, sponge_size, N*real(dx) )
          call SmoothStep (tmp3, real(iy)*dy   , sponge_size, N*real(dy) )
          call SmoothStep (tmp4, yl-real(iy)*dy, sponge_size, N*real(dy) )
          mask_sponge(ix,iy) = tmp1+tmp2+tmp3+tmp4
          if (mask_sponge(ix,iy)>1.0) mask_sponge(ix,iy)=1.0
      enddo
      enddo
  else
      mask_sponge=0.0
  endif

  mask_sponge=mask_sponge/epsilon_sponge
  call SaveField(trim(dir_name)//'/fields/'//trim(simulation_name)//'masksponge', mask_sponge, 1, xl,yl, "mask")


  
  
end subroutine create_sponge_mask




!==============================================================
! AUXILARY FUNCTIONS FOR THE MASK GENERATION
!===============================================================


subroutine SmoothStep (f,x,t,h)
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
  !--polynomial coefficients:
!   a =  1.0 / (4.0*(h**3))
!   b = -3.0*t / (4.0*(h**3))
!   c =  3.0*(t+h)*(t-h)/(4.0*(h**3))
!   d =  ((t+h)**2)*(2.0*h-t)/(4.0*(h**3))

!  if (x<=t-h) then
!    f = 1.0
!  elseif (((t-h)<x).and.(x<(t+h))) then
! !     f = a*(x**3) + b*(x**2) + c*x + d
! !     f = 1.0 - (x-t+h)/(2.0*h)
!    f = 0.5*(1.+cos((x-t+h)*pi/(2.0*h)) )
!  else
!    f = 0.0
!  endif

  !-----------------------------------
  ! version 14 - error function as non-oscilatory shape
  !-----------------------------------
  ! h - delta (gradient thickness)
  ! t - thickness (radius)
  if (sharp==.false.) then
      GradientERF = abs( ( exp(-(2.0*1.0)**2)  - 1.0 )/sqrt(pi) )
      delta = h*GradientERF
      f = 0.5*( erf( (t-x)/delta ) + erf( (x+t)/delta )  )
  else
      if (x<=t) f=1.d0
      if (x>t)  f=0.d0
  endif


end subroutine SmoothStep

!==============================================================================================================================================

real (kind=pr) function ConvertAngle(angle1)
  use share_vars
  implicit none

  real(kind=pr), intent (in) :: angle1
    !converts an angle from fortran standard (-pi<angle<pi) to real (0<angle<2pi)
  ConvertAngle = angle1
  if (angle1<0.0) then
    ConvertAngle = 2.0*pi+angle1
  endif

  return
end function ConvertAngle

!==============================================================================================================================================

subroutine DrawHinge (pointx, pointy, alpha1, alpha2, vx, vy, t, N)
! DrawHinge draws just a circle sector, not an entire circle.
  use share_vars
  implicit none
  real (kind=pr), intent (in)  :: pointx, pointy, vx, vy, t, N, alpha1, alpha2
  real (kind=pr) :: t_star, dx, dy, R, temp, gamma, ConvertAngle
  integer :: ixmax, ixmin, iymax, iymin, i, j

  dx=xl/real(nx)
  dy=yl/real(ny)
  t_star = t + 10.*N*max(dx,dy) !effective beam thickness including the smoothing zone
  ixmin = nint( (pointx-t_star)/dx)
  ixmax = nint( (pointx+t_star)/dx)
  iymin = nint( (pointy-t_star)/dy)
  iymax = nint( (pointy+t_star)/dy)

  ixmin = max( ixmin, 0)  
  ixmax = min( ixmax,nx-1)
  iymin = max( iymin, 0)  
  iymax = min( iymax,ny-1)  
  

  
  do i=ixmin, ixmax
    do j=iymin, iymax
      R = sqrt ( (real(i)*dx-pointx)**2 + (real(j)*dy-pointy)**2 )
      gamma = (atan2( (real(j)*dy-pointy),(real(i)*dx-pointx) ))

      if ((R <= t_star).and. &  !inside the beam
      ( ((0.5*pi+alpha2<=gamma).and.(gamma<=pi)).or.((-pi<=gamma).and.(gamma<=alpha2-pi*0.5)) ) & !half circle 2: backward facing
      .and.((gamma>=alpha1-0.5*pi).and.(gamma<=alpha1+0.5*pi)) )   then !half circle 1: foreward facing

      call SmoothStep(temp, R, t, N*max(dx,dy) )
      if (temp>mask(i,j)) then
        mask(i,j) = temp !existing mask points are overwritten if they are smaller than the new one
        maskvx(i,j)=vx
        maskvy(i,j)=vy
      endif

      endif
    enddo
  enddo
  

end subroutine DrawHinge

!==============================================================================================================================================

subroutine DrawEndpointLeft (pointx, pointy, alpha1, alpha2, vx, vy, t, N)
! DrawHinge draws just a circle sector, not an entire circle, between alpha1 and alpha2. fills NOT only the tiny gap.
  use share_vars
  implicit none
  real (kind=pr), intent (in)  :: pointx, pointy, vx, vy, t, N, alpha1, alpha2
  real (kind=pr) :: t_star, dx, dy, R, temp, gamma, ConvertAngle
  integer :: ixmax, ixmin, iymax, iymin, i, j

  dx=xl/real(nx)
  dy=yl/real(ny)
  t_star = t + N*max(dx,dy) !effective beam thickness including the smoothing zone
  ixmin = nint( (pointx-t_star)/dx)
  ixmax = nint( (pointx+t_star)/dx)
  iymin = nint( (pointy-t_star)/dy)
  iymax = nint( (pointy+t_star)/dy)


  if ((ixmin<0).or.(ixmax>nx).or.(iymin<0).or.(iymax>ny)) then
    write (*,*) "Mask:: Beam Hinge. Index out of region"
    stop
  endif


  do i=ixmin, ixmax
    do j=iymin, iymax
      R = sqrt ( (real(i)*dx-pointx)**2 + (real(j)*dy-pointy)**2 )
      gamma = ConvertAngle(atan2( (real(j)*dy-pointy),(real(i)*dx-pointx) ))
      if (((R <= t_star).and.((gamma>=alpha1).and.(gamma<=alpha2))).or.(R<=max(dx,dy)))   then

	call SmoothStep(temp, R, t, N*max(dx,dy) )
	if (temp>mask(i,j)) then
	  mask(i,j) = temp !existing mask points are overwritten if they are smaller than the new one
	  if (mask(i,j)>0.0) then ! if the mask is set there, give it a velocity
	    maskvx(i,j)=vx!*mask(i,j)*eps
	    maskvy(i,j)=vy!*mask(i,j)*eps
	  endif
	endif

      endif
    enddo
  enddo


end subroutine DrawEndpointLeft

!==============================================================================================================================================

subroutine DrawEndpointRight (pointx, pointy, alpha1, alpha2, vx, vy, t, N)
! DrawHinge draws just a circle sector, not an entire circle, between alpha1 and alpha2. fills NOT only the tiny gap.
  use share_vars
  implicit none
  real (kind=pr), intent (in)  :: pointx, pointy, vx, vy, t, N, alpha1, alpha2
  real (kind=pr) :: t_star, dx, dy, R, temp, gamma, ConvertAngle
  integer :: ixmax, ixmin, iymax, iymin, i, j

  dx=xl/real(nx)
  dy=yl/real(ny)
  t_star = t + N*max(dx,dy) !effective beam thickness including the smoothing zone
  ixmin = nint( (pointx-t_star)/dx)
  ixmax = nint( (pointx+t_star)/dx)
  iymin = nint( (pointy-t_star)/dy)
  iymax = nint( (pointy+t_star)/dy)


  if ((ixmin<0).or.(ixmax>nx).or.(iymin<0).or.(iymax>ny)) then
    write (*,*) "Mask:: Beam Hinge. Index out of region"
    stop
  endif

  
  do i=ixmin, ixmax
    do j=iymin, iymax
      R = sqrt ( (real(i)*dx-pointx)**2 + (real(j)*dy-pointy)**2 )
      gamma = (atan2( (real(j)*dy-pointy),(real(i)*dx-pointx) ))
      if (((R <= t_star).and.((gamma>=alpha1).and.(gamma<=alpha2))).or.(R<=max(dx,dy)))   then
	call SmoothStep(temp, R, t, N*max(dx,dy) )
	if (temp>mask(i,j)) then
	  mask(i,j) = temp !existing mask points are overwritten if they are smaller than the new one
	  if (mask(i,j)>0.0) then ! if the mask is set there, give it a velocity
	    maskvx(i,j)=vx
	    maskvy(i,j)=vy
	  endif
	endif
      endif
    enddo
  enddo
  
end subroutine DrawEndpointRight

!==============================================================================================================================================

subroutine DrawADot (pointx, pointy, vx, vy, t, N)
  use share_vars
  implicit none
  !-----------------------------------------------------------------
  ! This subroutine draws a smoothed point on the mask and two sharp ones on
  ! maskvx and maskvy.
  ! It is important that DrawADot is called before DrawBeamSegment
  ! because DrawBeamSegment overrides the points which are on the segment AND on the point.
  ! like this, the hinge fills only the Gap with the local velocity.
  !--------------------------------------------------------------------
  real (kind=pr), intent (in)  :: pointx, pointy, vx, vy, t, N
  real (kind=pr) :: t_star, dx, dy, R, temp
  integer :: ixmax, ixmin, iymax, iymin, i, j

  dx=xl/real(nx)
  dy=yl/real(ny)
  t_star = t + 3.0*N*max(dx,dy) !effective beam thickness including the smoothing zone
  ixmin = 0!max( nint( (pointx-t_star)/dx), 0)
  ixmax = nx-1!min( nint( (pointx+t_star)/dx), nx)
  iymin = 0!max( nint( (pointy-t_star)/dy), 0)
  iymax = ny-1!min( nint( (pointy+t_star)/dy), ny)

  if ((ixmin<0).or.(ixmax>nx-1).or.(iymin<0).or.(iymax>ny-1).or.(pointy<0.0)) then
    write (*,*) "Mask:: Beam Hinge. Index out of region"
    stop
  endif

  !$omp parallel do private(i,j,R,temp)
  do i=ixmin, ixmax
    do j=iymin, iymax
      R = sqrt ( (real(i)*dx-pointx)**2 + (real(j)*dy-pointy)**2 )
      call SmoothStep(temp, R, t, N*max(dx,dy) )
      if ( temp>mask(i,j) ) then!existing mask points are overwritten if they are smaller than the new one
        mask(i,j) = temp
        maskvx(i,j)=vx
        maskvy(i,j)=vy
      endif
    enddo
  enddo
  !$omp end parallel do
  
end subroutine DrawADot

!==============================================================================================================================================

subroutine DrawBeamSegment ( p1x,p1y,p2x,p2y, v1x, v1y, v2x,v2y, t, N )
  use share_vars
  implicit none
  real (kind=pr), intent (in)  :: t, N
  real (kind=pr), intent (in) :: p1x,p1y,p2x,p2y, v1x, v1y, v2x,v2y
  real (kind=pr), dimension(1:2) :: point1, point2, v1, v2
  real (kind=pr), dimension(1:2) :: point1_star, point2_star
  real (kind=pr) :: t_star, dx, dy, R, x_star, y_star, temp, alpha
  integer :: ixmax, ixmin, iymax, iymin, i, j

  point1(1) = p1x
  point1(2) = p1y
  point2(1) = p2x
  point2(2) = p2y
  v1(1) = v1x
  v1(2) = v1y
  v2(1) = v2x
  v2(2) = v2y
  
  !--------------------------------------------------------------
  ! This routine draws a smooth rectangle under an arbitrary angle
  ! defined by the two points.
  ! This is done by turning the coordinates by "angle", in this the box is
  ! calculated analytically.
  ! For the velocity, the rectangle is filled with linear
  ! interpolated velocitys of the 2 beampoints
  !--------------------------------------------------------------
  dx=xl/real(nx)
  dy=yl/real(ny)
  t_star = t + 3.0*N*max(dx,dy) !effective half beam thickness including the smoothing zone
  !real angle between the two points
  
  alpha = -atan2( (point2(2)-point1(2)),(point2(1)-point1(1)) )

  !domain for calculation (to reduce comp. costs)
  ixmin = nint( (min(point1(1),point2(1))-t_star)/dx)
  ixmax = nint( (max(point1(1),point2(1))+t_star)/dx)

  iymin = nint( (min(point1(2),point2(2))-t_star)/dy)
  iymax = nint( (max(point1(2),point2(2))+t_star)/dy)
  
  ixmin = max( ixmin, 0)  
  ixmax = min( ixmax,nx-1)
  iymin = max( iymin, 0)  
  iymax = min( iymax,ny-1)  

  ! star points are in rotated coordinates 
  !!!!!!!!!!
  ! note there is a sign error. cf DrawSharpEnd
  !!!!!!!!!!
  point1_star(1) = cos(alpha)*point1(1) - sin(alpha)*point1(2)
  point1_star(2) = sin(alpha)*point1(1) + cos(alpha)*point1(2)
  point2_star(1) = cos(alpha)*point2(1) - sin(alpha)*point2(2)
  point2_star(2) = sin(alpha)*point2(1) + cos(alpha)*point2(2)

  
  do i=ixmin, ixmax
    do j=iymin, iymax
      x_star = real(i)*dx*cos(alpha) -real(j)*dy*sin(alpha)
      y_star = real(i)*dx*sin(alpha) +real(j)*dy*cos(alpha)
      if ((x_star>=point1_star(1)).and.(x_star<=point2_star(1)).and.(y_star<=(point1_star(2)+t_star)).and.(y_star>=(point2_star(2)-t_star))) then
      call SmoothStep(temp, abs(y_star - point1_star(2)), t, N*max(dx,dy) )
      if (temp>mask(i,j)) then
        mask(i,j) = temp !existing mask points are overwritten if they are smaller than the new one
        ! linear interpolated velocitys
        maskvx(i,j) = (v1(1) + (x_star-point1_star(1)) * (v2(1)-v1(1)) / ds)
        maskvy(i,j) = (v1(2) + (x_star-point1_star(1)) * (v2(2)-v1(2)) / ds)
      endif
      endif
    enddo
  enddo
  
end subroutine DrawBeamSegment



subroutine DrawSharpEnd(point_end, alpha, v,t,N)
  use share_vars
  implicit none
  real (kind=pr), intent (in)  :: N,t, alpha
  real (kind=pr), dimension(1:2), intent (in)  :: point_end, v
  real (kind=pr), dimension(1:2) :: point1_star, point2_star, point2, point1
  real (kind=pr) :: t_star, dx, dy, R, x_star, y_star, temp,tmp2 
  integer :: ixmax, ixmin, iymax, iymin, i, j

  dx=xl/real(nx)
  dy=yl/real(ny)
  t_star = t + 3.0*N*max(dx,dy) !effective half beam thickness including the smoothing zone
  
  ! the first point lies INSIDE the beam
  point1(1) = point_end(1) - 5.0*max(dx,dy) * cos(alpha)
  point1(2) = point_end(2) - 5.0*max(dx,dy) * sin(alpha)
   
  ! the second point lies OUTSIDE the beam 
  point2(1) = point_end(1) + 5.0*max(dx,dy) * cos(alpha) 
  point2(2) = point_end(2) + 5.0*max(dx,dy) * sin(alpha)

  !domain for calculation (to reduce comp. costs)
  ixmin = nint( (min(point1(1),point2(1))-t_star)/dx)
  ixmax = nint( (max(point1(1),point2(1))+t_star)/dx)
  iymin = nint( (min(point1(2),point2(2))-t_star)/dy)
  iymax = nint( (max(point1(2),point2(2))+t_star)/dy)  
  
  ixmin = max( ixmin, 0)  ! ensure valid coordinates
  ixmax = min( ixmax,nx-1)
  iymin = max( iymin, 0)  
  iymax = min( iymax,ny-1)  

  ! star points are in rotated coordinates (note: I think I corrected the sign error made in 2009 here)
  point1_star(1) = cos(alpha)*point1(1) + sin(alpha)*point1(2)
  point1_star(2) =-sin(alpha)*point1(1) + cos(alpha)*point1(2)
  point2_star(1) = cos(alpha)*point2(1) + sin(alpha)*point2(2)
  point2_star(2) =-sin(alpha)*point2(1) + cos(alpha)*point2(2)

    
  do i=ixmin, ixmax
    do j=iymin, iymax
      x_star =  real(i)*dx*cos(alpha) + real(j)*dy*sin(alpha)
      y_star = -real(i)*dx*sin(alpha) + real(j)*dy*cos(alpha)
      
      if ((x_star>=point1_star(1)).and.(x_star<=point2_star(1)).and.(y_star<=(point1_star(2)+t_star)).and.(y_star>=(point2_star(2)-t_star))) then
	! note we start earlier in the beam, therefore 'zero' is actually inside the beam, thus we HAVE a thickness
	call SmoothStep (tmp2, x_star-point1_star(1), 2.0*max(dx,dy) , N*max(dx,dy) )
	! this is smoothing in normal direction
	call SmoothStep (temp, abs(y_star - point1_star(2)), t, N*max(dx,dy) )	
	
	mask(i,j) = temp*tmp2
	maskvx(i,j) = v(1)
	maskvy(i,j) = v(2)
      endif
    enddo
  enddo
  
    
  
end subroutine DrawSharpEnd


!===================================================================================================================================================
!===================================================================================================================================================
!===================================================================================================================================================
!  In the following you can find the subroutines used for the mask transport equation. they are currently not used, but maye useful one day
!===================================================================================================================================================
!===================================================================================================================================================
!===================================================================================================================================================

! subroutine RK4 ( dt, cylinder )
!   use share_vars
!   implicit none
!   real (kind=pr), dimension (0:nx-1,0:ny-1) :: mask1, mask2, mask3, mask4
!   real (kind=pr), intent (in) :: dt
!   real (kind=pr), dimension (1:6),         intent (in) :: cylinder
! 
! 
!   mask1 = mask
!   call MaskRHS(mask1,cylinder)
!   mask2     = mask + 0.5*dt*mask1
!   call MaskRHS(mask2,cylinder) 
!   mask3     = mask + 0.5*dt*mask2
!   call MaskRHS(mask3,cylinder)  
!   mask4     = mask + dt*mask3
!   call MaskRHS(mask4,cylinder) 
!   mask       = mask + dt*(mask1 + 2.*mask2 + 2.*mask3 + mask4 )/6.0
!   
! end subroutine RK4
! 
! subroutine MaskRHS(mask1,cylinder)
!   use share_vars
!   implicit none  
!   real (kind=pr), dimension (0:nx-1,0:ny-1), intent (inout) :: mask1
!   real (kind=pr), dimension (1:6),         intent (in) :: cylinder
!   real (kind=pr), dimension (0:nx-1,0:ny-1) :: mask_k, mask_dy
!   real (kind=pr) :: u0=-1.0
! 
!   call coftxy(mask1,mask_k)
!   call cofdy(mask_k, mask_dy)
!   call cofitxy(mask_dy,mask1)
!   mask1=-1.0*mask1*u0
! 
! end subroutine MaskRHS
! 
! subroutine EvolveMaskAB2(dt0,dt1,n0,n1,mask_RHS,cylinder)
!   use share_vars
!   use FieldExport
!   implicit none
!   real (kind=pr), intent (in) :: dt0,dt1
!   integer, intent (in):: n0,n1
!   real (kind=pr), dimension(0:nx-1,0:ny-1) :: mask_k!, mask_dy
!   real (kind=pr), dimension(0:nx-1,0:ny-1,0:1),intent(inout) :: mask_RHS
!   ! Beam indices: 1=beam_x 2=beam_y 3=beam_vx 4=beam_vy 5=theta 6=theta_dot
!   real (kind=pr), dimension (1:6),         intent (in) :: cylinder
!   real (kind=pr) :: b10,b11
!   integer ::i
! 
!   mask_RHS(:,:,n0)=mask
!   call MaskRHS(mask_RHS(:,:,n0),cylinder)
! 
!   b10 = dt1/dt0 * (0.5*dt1 + dt0)
!   b11 = -0.5 * dt1**2 / dt0
! 
!   mask = mask + b10*mask_RHS(:,:,n0) + b11*mask_RHS(:,:,n1)
!   
!   maskvx = cylinder(3)
!   maskvy = cylinder(4)
! 
! 
! end subroutine EvolveMaskAB2
! 
! 
! subroutine EvolveMask(dt1,time,n0,n1,mask_RHS,cylinder)
!   use share_vars
!   use FieldExport
!   implicit none
!   real (kind=pr), intent (in) :: dt1,time
!   integer, intent (in):: n0,n1
!   real (kind=pr), dimension(0:nx-1,0:ny-1) :: mask_k, mask_dy
!   complex (kind=pr), dimension(0:ny-1) :: expx0
!   real (kind=pr), dimension(0:nx-1,0:ny-1,0:1),intent(inout) :: mask_RHS
!   ! Beam indices: 1=beam_x 2=beam_y 3=beam_vx 4=beam_vy 5=theta 6=theta_dot
!   real (kind=pr), dimension (1:6),         intent (inout) :: cylinder
!   real (kind=pr) :: coefx, scale1,a,b,c,d,k
! integer ::i, Nsub=10,iy,ix,kx,ky
! 
! 
! !!!!!!!!!!!!!!!!
! cylinder(4)=-1.0
! !!!!!!!!!!!!!!!!!!
! 
! 
! ! when used as startup for AB2:
! mask_RHS(:,:,n0)=mask
! call MaskRHS(mask_RHS(:,:,n0),cylinder)
! mask=mask+dt1*mask_RHS(:,:,n0)
! 
! !do i=1,Nsub
! !call RK4(dt1/real(Nsub),cylinder)
! !call RK4(dt1,cylinder)
! !enddo 
! 
! !----------------------------------------
! ! exact version
! 
! 
! 
! 
! call coftxy(mask_sponge,mask_k) !here, mask_sponge is the mask at time t=0, in order not to accumulate an error
! 
! !     ck = ak + i bk with ak= fk(k,2l) and bk= fk(k,2l+1)
! !                         for l=0, ky-1 and all k=0,nx-1
!   scale1 = 2.0*pi/yl
!   !$omp parallel do private(kx,ky)
!   do kx = 0, nx-1
!      do ky = 0, ny-2, 2 ! loop over REAL parts of the FFT
!        ! is a mutliplication by a complex exponential, nothing but a multiplication of two complex numbers (a+ib)*(c+id)
!         k = real(ky/2)*scale1
!         a = cos(-k*cylinder(4)*time)
!         b = sin(-k*cylinder(4)*time)
!         c = mask_k(kx, ky)
!         d = mask_k(kx, ky+1) !imag         
!         mask_dy(kx, ky) = a*c - b*d
!      end do
!      do ky = 1, ny-1, 2 ! loop over imaginary parts of the FFT
!       ! actually, both (real/imag) have of course the same wavenumber. 
!         k = real((ky-1)/2)*scale1
!         a = cos(-k*cylinder(4)*time)
!         b = sin(-k*cylinder(4)*time)
!         c = mask_k(kx, ky-1) !real
!         d = mask_k(kx, ky) !imag
!         mask_dy(kx, ky) = a*d + b*c
!      end do
!   end do
!   !$omp end parallel do
! 
!   call cofitxy(mask_dy, mask)
! 
!   if (maxval(mask)<0.9) then
!   write (*,*) "no value"
!   stop
!   endif
! !----------------------------------------
! 
! 
! maskvx = cylinder(3)
! maskvy = cylinder(4)
! 
! 
! end subroutine EvolveMask


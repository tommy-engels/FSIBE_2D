module SpectralSolver

  implicit none
  contains

! ==========================================================================================================================

subroutine cal_nlk (time, dt1, nlk, vortk, beam, u, f_mean)
  !---------------------------------------------------------------
  !     determine fourier coeffs of the r.h.s. that is due to
  !     the NON-linear parts and forcing (ie penalisation)
  !---------------------------------------------------------------
  use share_vars
  use FieldExport
  implicit none
  ! input 
  real (kind=pr), intent (in) :: 					time
  real (kind=pr), dimension (0:nx-1, 0:ny-1), intent (in) :: 		vortk
  real (kind=pr), dimension (0:nx-1, 0:ny-1,1:2), intent (in) :: 	u
  real (kind=pr), dimension (0:nx-1, 0:ny-1,1:2) :: 			penal
  real (kind=pr), dimension (0:ns-1, 1:6), intent (in) :: 		beam
  ! output
  real (kind=pr), dimension (0:nx-1, 0:ny-1), intent (out) :: 		nlk
  real (kind=pr), dimension (0:nx-1, 0:ny-1) :: 			vort_init
  real (kind=pr), intent (out) :: 					dt1
  real (kind=pr), intent (out), dimension(1:2) :: 			f_mean
  ! local variables 
  real (kind=pr), dimension (0:nx-1, 0:ny-1) :: vort_dx, vort_dy, stream, work2, penvortk
  real (kind=pr) :: 				u_max, Mean_ux, Mean_uy,  y_chan, H_effective,dx, dy
  integer :: 					iy,ix, GetIndex
  
  dx = xl / real (nx)
  dy = yl / real (ny)

  !-----------------------------------------------------------------------------------------
  ! Calculate nonlinear term in physical space      (-u * nabla(vort)   
  !-----------------------------------------------------------------------------------------
  
  !--Calculate dx and dy of the vorticity
  call cofdx (vortk, work2)
  call cofitxy (work2, vort_dx)

  call cofdy (vortk, work2)
  call cofitxy (work2, vort_dy)

  !$omp parallel do private(iy)
  do iy=0,ny-1
     work2(:,iy) = - u(:,iy,1)*vort_dx(:,iy) - u(:,iy,2)*vort_dy(:,iy) !sign changed
  enddo
  !$omp end parallel do

  call coftxy (work2, nlk) !--Calculate fourier coefficient of nonlinear term
  
  !-----------------------------------------------------------------------------------------
  ! Calculate time step
  !-----------------------------------------------------------------------------------------
  !-- CFL condition for the fluid
  
  !$omp parallel do private(iy)
  do iy=0,ny-1
     work2(:,iy) = u(:,iy,1)**2 + u(:,iy,2)**2
  enddo
  !$omp end parallel do

  u_max = sqrt(maxval(work2))
  dt1 = cfl * min (dx, dy) / u_max

  !-- u_max is very very small
  if (u_max < 1.e-10) dt1 = 1.e-2 
  !-- Time stepping control for volume penalization
  if ((iBeam>0).or.(iWalls>0).or.(iCylinder>0)) dt1 = min(0.9*eps,dt1)
  !-- Fixed time step
  if (dt_fixed>0.0) dt1 = min(dt_fixed, dt1)  
  !-- perfectly reach target time
  if ((Time_end - time) < dt1) then
    dt1 = Time_end - time
    write (*,*) "last time step:", dt1
  endif
    
  open (14, file = trim(dir_name)//'/'//trim(simulation_name)//'u_max', status = 'unknown', access = 'append') ! Append output data file
  write (14,'(1x, 5(es12.4,1x))') time, dt1, u_max, sum(u(:,:,1))/real(nx*ny), sum(u(:,:,2))/real(nx*ny)
  close (14)
  
  !-----------------------------------------------------------------------------------------
  ! Calculate boundary penalization term in Fourier space
  !-----------------------------------------------------------------------------------------

  if ((iBeam>0).or.(iWalls>0).or.(iCylinder>0)) then
    !$omp parallel do private(iy)
    do iy=0,ny-1
      penal(:,iy,1) = -mask(:,iy) * (u(:,iy,1) - maskvx(:,iy))  ! ux: x-component of the relative velocity
      penal(:,iy,2) = -mask(:,iy) * (u(:,iy,2) - maskvy(:,iy))  ! uy: y-component of the relative velocity
    enddo
    !$omp end parallel do
    
    !--Compute mean value of penalization term, if imposing zero mean pressure
    if ( iMeanVelocity == 11 ) then
      f_mean(1) = sum(penal(:,:,1)) / real(nx*ny)
      f_mean(2) = sum(penal(:,:,2)) / real(nx*ny)
    endif

    !--Calculate x-derivative in Fourier space
    call coftxy (penal(:,:,2), work2) ! fourier transformation (work2::out)
    call cofdx (work2, penal(:,:,2) )  ! derivative dPenal_y / dx

    !--Calculate y-derivative in Fourier space
    call coftxy (penal(:,:,1), work2)
    call cofdy (work2, penal(:,:,1) ) !derivative dPenal_x / dy
  else
    penal = 0.0
  endif
    
  
  !-----------------------------------------------------------------------------------------
  ! Calculate vorticity sponge term
  !-----------------------------------------------------------------------------------------

  if ((iSponge>0).and.(iSponge<4)) then
  
      !------------------------------------------------------------------------------------
      !--Sponge modification:
      !------------------------------------------------------------------------------------
      ! mask_sponge is now global and only set once at startup to save a bit of time
      penvortk=vortk !stupid: you cannot modify an input argument, therefore duplicate
      call cofitxy(penvortk, work2) !transformation into physical space      

      if (iSpongeType ==1) then !poisseuille sponge - forces linear vorticity profile
      
          !set the vorticity profile
          H_effective=yl-2.0*h_channel
          
          !$omp parallel do private(iy,y_chan)
          do iy=0,ny-1
            ! the Poisseuille flow requires a linear vorticity
            y_chan = real(iy)*dy-h_channel
            if ((y_chan>0.0).and.(y_chan<H_effective)) then
            vort_init(:,iy) = (-6.0*Mean_ux(time)/H_effective)*(1.0-2.0*y_chan/H_effective)
            else
            vort_init(:,iy) = 0.0
            endif
          enddo
          !$omp end parallel do
          
          !$omp parallel do private(iy)
          do iy=0,ny-1
            work2(:,iy) = mask_sponge(:,iy) * (work2(:,iy)-vort_init(:,iy)) ! perform penalization = multiply with the spongemask
          enddo
          !$omp end parallel do
          
      elseif (iSpongeType ==2) then !No vorticity sponge
      
          !$omp parallel do private(iy)
          do iy=0,ny-1
            work2(:,iy) = mask_sponge(:,iy) * work2(:,iy) ! perform penalization = multiply with the spongemask
          enddo
          !$omp end parallel do
          
      endif     

      call coftxy(work2,penvortk) ! and back to fourier space. now penvortk is penalized with a sponge
      
      ! add the actual term (in Fourier space)
      !$omp parallel do private(iy)
      do iy=0,ny-1
	!nlk already contains the nonlinear term (u*NABLA(omega))
        nlk(:,iy) = nlk(:,iy) + (penal(:,iy,2) - penal(:,iy,1)) - penvortk(:,iy) !now the penalized vorticity goes on the RHS (sponge)
      enddo
      !$omp end parallel do      
      
  else !without vorticity sponge  
  
      !$omp parallel do private(iy)
      do iy=0,ny-1
        !nlk already contains the nonlinear term (u*NABLA(omega))
        nlk(:,iy) = nlk(:,iy) + (penal(:,iy,2) - penal(:,iy,1)) !sign changed
      enddo
      !$omp end parallel do
      
  endif

end subroutine cal_nlk

!===================================================================================================

subroutine cal_velocity(time, vortk, u, stream )
  ! this routine computes the streamfunction and the velocity field given a vorticity in fourier space
  use share_vars
  implicit none
  real (kind=pr), dimension (0:nx-1, 0:ny-1,1:2), intent (out) :: u
  real (kind=pr), dimension (0:nx-1, 0:ny-1), intent (in)  :: vortk
  real (kind=pr), dimension (0:nx-1, 0:ny-1), intent (out) :: stream
  real (kind=pr), dimension (0:nx-1, 0:ny-1) :: work2
  real (kind=pr), intent(in) :: time
  real (kind=pr) :: Mean_ux, Mean_uy
  integer :: iy

  !--Solve poisson equation for stream function
  call poisson (vortk, stream)

  !--Calculate x-derivative of stream function
  call cofdx (stream, work2)
  call cofitxy (work2, u(:,:,2)) ! u(:,:,2) here contains -uy = dpsi/dx

  !--Calculate y-derivative of stream function
  call cofdy (stream, work2)
  call cofitxy (work2, u(:,:,1)) ! u(:,:,1) here contains ux = dpsi/dy

  !--Add a mean velocity
  !$omp parallel do private(iy)
  do iy=0,ny-1
      u(:,iy,1) =   u(:,iy,1) + U_mean_true*Mean_ux(time)  ! +ux
      u(:,iy,2) = - u(:,iy,2) + U_mean_true*Mean_uy(time)  ! +uy, U_mean_true is used for channel flows
  enddo
  !$omp end parallel do

end subroutine cal_velocity

!===================================================================================================

subroutine pressure ( time, dt1, vortk, press, u, beam, Forces, FluidIntegralQuantities )
  !-----------------------------------------------------------------------------------------------
  ! --> Calculates the pressure field given vorticity in Fourier space <--
  !
  ! This version makes use of the Navier-Stokes eqn with the non-linear term in its
  ! rotational formulation. That means, we are left with vort x velocity (omega x u)
  ! and the total pressure q = p + 0,5*(u^T \cdot u). This formulation is faster to compute
  ! than the standard convevtive form. 
  ! Moreover, computational time can be saved noticing that the new velocity field, u^(n+1),
  ! can be re-used in the next step when computing the RHS (in cal_nlk). Therefore this 
  ! routine also returns the velocity field at the new time level.
  !-----------------------------------------------------------------------------------------------
  use share_vars
  use PerformanceMeasurement
  use FieldExport
  implicit none
  integer :: k, l, iy,ix,ixmin,ixmax,iymin,iymax
  real (kind=pr) :: quot
  real (kind=pr), intent(in) :: time,dt1
  real (kind=pr), dimension (0:nx-1, 0:ny-1), intent (out) :: press
  real (kind=pr), dimension (0:nx-1, 0:ny-1, 1:2), intent (out) :: u
  real (kind=pr), dimension (0:nx-1, 0:ny-1, 1:2) :: penal
  real (kind=pr), dimension (0:nx-1, 0:ny-1), intent (in)  :: vortk
  real (kind=pr), dimension (0:ns-1, 1:6), intent (in) :: beam
  real (kind=pr), dimension (1:4), intent(out) :: Forces !Forces: drag lift drag_unst lift_unst
  real (kind=pr), dimension (1:7), intent(out) :: FluidIntegralQuantities ! 1= vor_rms 2=vor_rms_dot 3=Fluid kinetic Enegry 4=enstrophy
  real (kind=pr) :: Mean_ux, Mean_uy
  real (kind=pr), save :: vor_rms_old=0.0
!--Working files
  real (kind=pr), dimension (0:nx-1, 0:ny-1) :: work1, work_dy, work_dx, work2
  real (kind=pr) :: dx, dy
  dx = xl / real (nx)
  dy = yl / real (ny)

  !------------------------------------------------
  ! compute velocity
  !------------------------------------------------
  call cal_velocity( time+dt1, vortk, u, work1) ! velocity u at time n+1, will be used for RHS at t_n in the next time step


  !------------------------------------------------
  ! compute penalization term
  !------------------------------------------------
  !$omp parallel do private(iy)
  do iy=0,ny-1
    !penalty term is positive in poisson-pressure eqn
    penal(:,iy,1) = mask(:,iy) * (u(:,iy,1) - maskvx(:,iy))  ! ux: x-component of the relative velocity
    penal(:,iy,2) = mask(:,iy) * (u(:,iy,2) - maskvy(:,iy))  ! uy: y-component of the relative velocity    
  enddo
  !$omp end parallel do


  !------------------------------------------------
  ! compute hydrodynamic forces
  ! as this is just the penalty term, we do it here. 
  !------------------------------------------------
  call cal_drag (time, dt1, penal, beam, Forces)
  
  !------------------------------------------------
  ! vor to physical space
  ! you need it for the non-linear term in the 
  ! pressure eqn, and you compute vor_rms
  !------------------------------------------------
  call cofitxy (vortk, work1) !-- vort in physical space
  
  !$omp sections
  !$omp section
  ! at this occasion, compute rms value and its derivative
  FluidIntegralQuantities(1) = dx*dy*sqrt ( sum(work1**2) )
  FluidIntegralQuantities(2) = (FluidIntegralQuantities(1)-vor_rms_old)/dt1      
  
  if ((vor_rms_old==0.0).and.(FluidIntegralQuantities(1).ne.0.0)) then
    FluidIntegralQuantities(2) =0.0 !first time step, we can't compute vor_rms_dot
  endif
    
  vor_rms_old = FluidIntegralQuantities(1)
  !$omp section
  ! fluid kinetic energy
  FluidIntegralQuantities(3) = 0.5d0*dx*dy*sum( ((u(:,:,1)-Mean_ux(time+dt1))**2 + (u(:,:,2)-Mean_uy(time+dt1))**2) )
  !$omp section
  ! enstrophy (energy dissipation term)
  FluidIntegralQuantities(4) = nu*dx*dy*sum( work1**2 )
  !$omp section
  ! dissipation in penalty term
  FluidIntegralQuantities(5) = dx*dy*sum( mask*((u(:,:,1)-maskvx)**2+(u(:,:,2)-maskvy)**2) )
  !$omp section
  ! energy inflow due to penalization
  FluidIntegralQuantities(6) = dx*dy*sum( mask*( maskvx*(u(:,:,1)-maskvx) + maskvy*(u(:,:,2)-maskvy)  ) )
  !$omp section
  ! energy inflow (note we use the full mask, not only the actual object)
  FluidIntegralQuantities(7) = dx*dy*u_mean_true*sum( Mean_ux(time+dt1)*penal(:,:,1)+Mean_uy(time+dt1)*penal(:,:,2) )
  !$omp end sections

  
  !------------------------------------------------
  ! for fixed beams (like the CFD test) it is not nessesairy to compute the pressure, so we skip what follows.
  !------------------------------------------------
  
  if ((iFLUSI==1).or.(iSaveBeam>0)) then     
      !------------------------------------------------
      ! DIVERGENCE of penalty term
      !------------------------------------------------
      !--Calculate x-derivative in Fourier space
      call coftxy (penal(:,:,1), work2) !-- d Px / dx
      call cofdx  (work2, penal(:,:,1))
      !--Calculate y-derivative in Fourier space
      call coftxy (penal(:,:,2), work2) !-- d Py / dy
      call cofdy  (work2, penal(:,:,2))

      !------------------------------------------------
      ! DIVERGENCE of non-linear term
      !------------------------------------------------
      !$omp parallel do private(iy)
      do iy=0,ny-1
	work_dx(:,iy) = -work1(:,iy)*u(:,iy,2) ! -vor*uy
	work_dy(:,iy) = +work1(:,iy)*u(:,iy,1) ! +vor*ux
      enddo
      !$omp end parallel do

      !compute divergence
      call coftxy(work_dx, work1)     !-- to fourier space
      call cofdx (work1  , work_dx)   !-- d(-vor*uy)/dx
      call coftxy(work_dy, work1)     !-- to fourier space
      call cofdy (work1  , work_dy)   !-- d(+vor*ux)/dy

      !--right-hand side of poisson eqn 
      !$omp parallel do private(iy)
      do iy=0,ny-1
	work1(:,iy) = penal (:,iy,1) + penal (:,iy,2)  & !-- divergence of penalty term (fourier-space)
		     +work_dy(:,iy)  + work_dx(:,iy)     !-- divergence of non-linear term (fourier-space)
      enddo
      !$omp end parallel do            
      
      !------------------------------------------------
      ! solve poisson eqn
      !------------------------------------------------
      call poisson (work1, press) !--Solve Poisson equation to get the pressure q. q is static pressure p + 0.5*v_abs**2  
      call cofitxy (press, work1) !--Transform pressure to physical space
      
      
      !$omp parallel do private(iy)
      do iy=0,ny-1
	press(:,iy) = work1(:,iy) - 0.5*(u(:,iy,1)**2 + u(:,iy,2)**2) !--substract dynamical pressure to get static pressure
      enddo
      !$omp end parallel do      
   else !------------fixed run, no pressure required. 
      !$omp parallel do private(iy)
      do iy=0,ny-1
	press(:,iy) = 0.0
      enddo
      !$omp end parallel do      
   endif

end subroutine pressure

  

!===================================================================================================

integer function GetIndex(ix,nx)
  implicit none
  integer, intent (in) ::ix,nx
  integer :: tmp
  tmp=ix
  if (tmp<0) tmp = tmp+nx
  if (tmp>nx-1) tmp = tmp-nx
  GetIndex=tmp
  return
end function GetIndex

!===================================================================================================

subroutine EvolveFluidExplicit(time, dt1, n0, n1, beam, vortk, workvis, nlk, press, u, Forces, FluidIntegralQuantities)
  use share_vars
  use FieldExport
  implicit none
  real (kind=pr), dimension (0:nx-1, 0:ny-1, 0:1), intent (inout) :: vortk, workvis, nlk
  real (kind=pr), dimension (0:nx-1, 0:ny-1), intent (inout) :: press
  real (kind=pr), dimension (0:nx-1, 0:ny-1) :: stream
  real (kind=pr), dimension (0:nx-1, 0:ny-1, 1:2), intent (inout) :: u
  real (kind=pr), dimension (0:ns-1, 1:6), intent (in) :: beam
  real (kind=pr), dimension (1:4), intent(out) :: Forces !Forces: drag lift drag_unst lift_unst
  ! 1= vor_rms 2=vor_rms_dot 3=Fluid kinetic Enegry 4=Dissipation 5=Dissipation(Penalty) 6=EnergyInflow
  real (kind=pr), dimension (1:7), intent(out) :: FluidIntegralQuantities
  real (kind=pr), intent (inout) :: dt1, time
  integer, intent (in) :: n0, n1
  integer :: iy
  !---------------------------------------------------------
  !     n0 = n   n1 = n-1
  !	n0 = n   n1 = n+1 (for vortk and workvis)
  !---------------------------------------------------------

  !--dealiase
  vortk(:,:,n0) = dealiase * vortk(:,:,n0)

  !--for the first time step, compute velocity and penalization term explicitly
  call cal_velocity(time, vortk(:,:,n0), u, stream)

  !--calculate RHS
  call cal_nlk ( time, dt1, nlk (:,:,n0), vortk(:,:,n0), beam , u, fmean(:,n0) )

  call vis ( dt1, workvis(:,:,n1) )
  
  !$omp parallel do private(iy)
  do iy=0,ny-1
    vortk(:,iy,n1) = (vortk(:,iy,n0) + dt1*nlk(:,iy,n0)) * workvis(:,iy,n1) * dealiase(:,iy)
  enddo
  !$omp end parallel do

  !--if imposing zero mean pressure, solve the corresponding eqn
  if (iMeanVelocity == 11 ) then
    u_mean(1) = u_mean(1) + fmean(1,n0)*dt1
    u_mean(2) = u_mean(2) + fmean(2,n0)*dt1
  endif
  
  
  call pressure(time,dt1, vortk(:,:,n1), press, u, beam, Forces, FluidIntegralQuantities)

  forces(3:4)=0.0 ! first time step: unsteady correction fails
end subroutine EvolveFluidExplicit

!==============================================================================

subroutine EvolveFluidAB2(time, dt0,dt1, n0, n1, beam,  vortk, workvis, nlk, u, press, Forces, FluidIntegralQuantities)
  use share_vars
  use FieldExport
  implicit none
! adams bashforth step
  real (kind=pr), dimension (0:nx-1, 0:ny-1, 0:1), intent (inout) :: vortk, workvis, nlk
  real (kind=pr), dimension (0:nx-1, 0:ny-1, 1:2), intent (inout) :: u
  real (kind=pr), dimension (0:nx-1, 0:ny-1), intent (out) :: press
  real (kind=pr), dimension (0:ns-1, 1:6), intent (in) :: beam
  real (kind=pr), dimension (1:4), intent(out) :: Forces !Forces: drag lift drag_unst lift_unst
  ! 1= vor_rms 2=vor_rms_dot 3=Fluid kinetic Enegry 4=Dissipation 5=Dissipation(Penalty) 6=EnergyInflow
  real (kind=pr), dimension (1:7), intent(out) :: FluidIntegralQuantities 
  real (kind=pr), intent (inout) :: dt0, dt1, time
  real (kind=pr) ::b10,b11
  integer, intent (in) :: n0, n1
  integer :: iy,ix


  !$omp parallel do private(iy)
  do iy=0,ny-1
      nlk(:,iy,n1) = nlk(:,iy,n1)*workvis(:,iy,n0) !OLD workvis, OLD NLK
  enddo
  !$omp end parallel do

  !--Calculate fourier coeffs of nonlinear rhs and forcing (is RHS_n !! )
  call cal_nlk ( time, dt1, nlk (:,:,n0), vortk(:,:,n0), beam, u, fmean(:,n0) )

  
  b10 = dt1/dt0 * (0.5*dt1 + dt0)
  b11 = -0.5 * dt1**2 / dt0
  
  
  call vis( dt1, workvis(:,:,n1) )
  
  !$omp parallel do private(iy)
  do iy=0,ny-1
      vortk(:,iy,n1) = (vortk(:,iy,n0) + b10*nlk(:,iy,n0) + b11*nlk(:,iy,n1)) * workvis(:,iy,n1)*dealiase(:,iy)
  enddo
  !$omp end parallel do

  !-- if imposing zero mean pressure, solve corresponding eqn
  if ( iMeanVelocity == 11 ) then
    u_mean = u_mean + b10*fmean(:,n0) + b11*fmean(:,n1) ! im not sure about the ordering b10 and b11. it seems to be consistent with fluid.
  endif 
  
  ! in pressure, we now compute the u, penal terms for the next time step.
  call pressure(time, dt1, vortk(:,:,n1), press, u, beam, Forces, FluidIntegralQuantities)
  
end subroutine EvolveFluidAB2

!==========================================================================================================================

subroutine cal_drag (time, dt1, penal, beam, Forces)
  use share_vars
  use FieldExport
  implicit none
  real (kind=pr), intent (in) 					:: time, dt1
  real (kind=pr), dimension (0:nx-1, 0:ny-1,1:2), intent (in) 	:: penal
  real (kind=pr), dimension (1:4), intent(out) :: Forces !Forces: drag lift drag_unst lift_unst
  real (kind=pr), dimension (1:6) :: LeadingEdge !LeadingEdge: x0, y0, vx, vy, ax, ay (Array)
  real (kind=pr), dimension (0:ns-1, 1:6), intent (in) :: beam
  real (kind=pr), dimension (0:ns-1, 1:2) :: accel
  real (kind=pr), save :: drag_unst_old, lift_unst_old, drag, lift, lift_unst, drag_unst, drag_unst_new, lift_unst_new
  real (kind=pr) 	:: dx, dy, norm, alpha, alpha_t, alpha_tt, h_star
  integer 		    :: ix, iy,ixmin,ixmax,iymin,iymax,i

  call mouvement(time, alpha, alpha_t, alpha_tt, LeadingEdge )
  

  !--Calculate lift and drag (relative to direction of mean flow)
  dx = xl/real(nx)
  dy = yl/real(ny)
  norm = dx*dy

  !------------------------------------------------------------------------------------------------------------------------------
  !--			Compute Drag and Lift forces & unsteady correction terms
  !------------------------------------------------------------------------------------------------------------------------------
  ! optimized 10/11/2011 now using penal from pressure subroutine and the smallest (best) integration windows possible  
  if (iBeam>0) then ! case of a beam with or without cylinder at the leading edge
    iymax = int( (max( maxval(beam(:,2))+4.0*t_beam, LeadingEdge(2)+2.0*r_cylinder))/dy)
    iymin = int( (min( minval(beam(:,2))-4.0*t_beam, LeadingEdge(2)-2.0*r_cylinder))/dy)
    ! absolute upper limit for xmin is the leading edge - 4 times thickness
    ixmin = int( min( (LeadingEdge(1)-max(2.5*r_cylinder, 4.0*t_beam)), (minval(beam(:,1))-4.0*t_beam) )  /dx)
    ! even if the beam is hanging, its minimum x value is the leading edge
    ixmax = int( (maxval(beam(:,1))+4.0*t_beam)  /dx)
  else
    iymax=ny-1
    iymin=0
    ixmax=nx-1
    ixmin=0
  endif

  h_star=h_channel+N_smooth*dy ! channel walls are always in y-direction
  
  if (iWalls>0) then
    if (ixmin<0)    			ixmin = 0	! for safety reasons, should never actually take place
    if (ixmax>nx-1) 			ixmax = nx-1
    if (iymin<nint(     h_star/dy)+2)   iymin = nint(h_star/dy)+2
    if (iymax>nint((yl-h_star)/dy)-2)   iymax = nint((yl-h_star)/dy)-2
  endif

  if (ixmin<0)    ixmin = 0
  if (ixmax>nx-1) ixmax = nx-1
  if (iymin<0)    iymin = 0
  if (iymax>ny-1) iymax = ny-1


  drag = sum(penal(ixmin:ixmax,iymin:iymax,1)) * norm
  lift = sum(penal(ixmin:ixmax,iymin:iymax,2)) * norm

  !------------------------------------------------------------------------------------------------------------------------------
  !--			compute unsteady lift/drag corrections
  !------------------------------------------------------------------------------------------------------------------------------
  ! numercial correction, works, no significant disadvantage w.r.t beamßelement method 28-3-2012
  drag_unst_new = sum (maskvx(ixmin:ixmax,iymin:iymax)*mask(ixmin:ixmax,iymin:iymax)) * norm * eps
  lift_unst_new = sum (maskvy(ixmin:ixmax,iymin:iymax)*mask(ixmin:ixmax,iymin:iymax)) * norm * eps  
  
  drag_unst   = (drag_unst_new - drag_unst_old)/dt1 !compute unsteady correction
  lift_unst   = (lift_unst_new - lift_unst_old)/dt1  
  drag_unst_old = drag_unst_new  !iterate
  lift_unst_old = lift_unst_new  !iterate
 

  Forces(1) = drag
  Forces(2) = lift
  Forces(3) = drag_unst
  Forces(4) = lift_unst  
  
! !   if (iBeam>0) then
! !       ! compute unsteady corrections via the beam 
! !       accel = (beam(:,3:4)-beam_tmp(:,3:4))/dt1 ! first order derivative
! !       Forces(3) = sum(accel(:,1))*2.0*ds*t_beam
! !       Forces(4) = sum(accel(:,2))*2.0*ds*t_beam
! !       beam_tmp=beam !beam_tmp is global. why? because an automatic object must not have a save attribute. sometimes, I hate fortran.
! !   endif
  
end subroutine cal_drag

!==========================================================================================================================================================


subroutine ComputePressureSnapshot ( time, vortk, press )
  !---------------------------------------------------------------
  !     Calculates the pressure given vorticity in Fourier space
  !     and mean velocity
  !     Corrected, see a note on 28 may 2008 (Dmitry)
  !
  !   Actually, this subroutine is just to compute the pressure to make a video snapshot; normally, it is computed in "pressure" (when using it in every time step)
  !
  !---------------------------------------------------------------

  use share_vars
  use PerformanceMeasurement
  use FieldExport
  implicit none
  integer :: k, l, iy,ix,ixmin,ixmax,iymin,iymax
  real (kind=pr), intent(in) :: time
  real (kind=pr), dimension (0:nx-1, 0:ny-1), intent (out) :: press
  real (kind=pr), dimension (0:nx-1, 0:ny-1, 1:2) :: penal, u
  real (kind=pr), dimension (0:nx-1, 0:ny-1), intent (in)  :: vortk  
  real (kind=pr) :: Mean_ux, Mean_uy
!--Working files
  real (kind=pr), dimension (0:nx-1, 0:ny-1) :: work1, work_dy, work_dx, work2
  real (kind=pr) :: dx, dy

  dx = xl / real (nx)
  dy = yl / real (ny)

  call cal_velocity( time, vortk, u, work1) ! velocity u at time n+1, will be used for RHS at t_n in the next time step
                                            ! work1 is streamfunction for now
                                            
  !penalty term is positive in poisson-pressure eqn
  penal(:,:,1) = mask * (u(:,:,1) - maskvx)  ! ux: x-component of the relative velocity
  penal(:,:,2) = mask * (u(:,:,2) - maskvy)  ! uy: y-component of the relative velocity  
  !--------
  ! compute vor_rms and its time derivative
  call cofitxy (vortk, work1) !-- vort in physical space
  
  !========================================================================================================================
  ! DIVERGENCE of penalty term
  !========================================================================================================================
  !--Calculate x-derivative in Fourier space
  call coftxy (penal(:,:,1), work2) !-- d Px / dx
  call cofdx  (work2, penal(:,:,1))
  !--Calculate y-derivative in Fourier space
  call coftxy (penal(:,:,2), work2) !-- d Py / dy
  call cofdy  (work2, penal(:,:,2))

  !========================================================================================================================
  ! DIVERGENCE of non-linear term
  !========================================================================================================================
  work_dx = -work1*u(:,:,2) ! -vor*uy
  work_dy = +work1*u(:,:,1) ! +vor*ux
  
  !compute divergence
  call coftxy(work_dx, work1)     !-- to fourier space
  call cofdx (work1  , work_dx)   !-- d(-vor*uy)/dx
  call coftxy(work_dy, work1)     !-- to fourier space
  call cofdy (work1  , work_dy)   !-- d(+vor*ux)/dy

  !--right-hand side of poisson eqn  
  work1 = penal (:,:,1) + penal (:,:,2)  & !-- divergence of penalty term (fourier-space)
	+ work_dy  + work_dx     !-- divergence of non-linear term (fourier-space)
  !========================================================================================================================
  ! solve poisson eqn
  !========================================================================================================================  
  call poisson (work1, press) !--Solve Poisson equation to get the pressure q. q is static pressure p + 0.5*v_abs**2  
  call cofitxy (press, work1) !--Transform pressure to physical space
  press = work1 - 0.5*(u(:,:,1)**2 + u(:,:,2)**2) !--substract dynamical pressure to get static pressure
   
end subroutine ComputePressureSnapshot


end module SpectralSolver

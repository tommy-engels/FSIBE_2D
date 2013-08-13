module SpectralSolver

  implicit none
  contains

! ==========================================================================================================================

subroutine cal_nlk (time, dt1, nlk, vortk, u, f_mean)
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
  ! output
  real (kind=pr), dimension (0:nx-1, 0:ny-1), intent (out) :: 		nlk
  real (kind=pr), dimension (0:nx-1, 0:ny-1) :: 			vort_init
  real (kind=pr), intent (out) :: 					dt1
  real (kind=pr), intent (out), dimension(1:2) :: 			f_mean
  ! local variables 
  real (kind=pr), dimension (0:nx-1, 0:ny-1) :: vort_dx, vort_dy, work2, penvortk
  real (kind=pr) :: 				u_max,Mean_ux, y_chan, H_effective,dx, dy
  integer :: 					iy
  
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




subroutine cal_energy ( time, dt1, vor, u, penal )
  use share_vars
  implicit none
  real (kind=pr), intent (in) :: time, dt1
  real (kind=pr), intent (in) :: vor(0:nx-1,0:ny-1)
  real (kind=pr), intent (in) :: u(0:nx-1,0:ny-1,1:2)
  real (kind=pr), intent (in) :: penal(0:nx-1,0:ny-1,1:2)
  real (kind=pr) :: dx,dy, vor_rms, dissipation, porous_diss, inflow_solid, inflow_mean, kinetic, Mean_ux, Mean_uy
  dx = xl / real (nx)
  dy = yl / real (ny)
    
  !$omp parallel sections
    !$omp section
    ! at this occasion, compute rms value and its derivative
    vor_rms = dx*dy*sqrt ( sum(vor**2) )

    !$omp section
    ! fluid kinetic energy
    kinetic = 0.5d0*dx*dy*sum( ((u(:,:,1)-Mean_ux(time+dt1))**2 + (u(:,:,2)-Mean_uy(time+dt1))**2) )
    
    !$omp section
    ! enstrophy (energy dissipation term)
    dissipation = nu*dx*dy*sum( vor**2 )
    
    !$omp section  
    ! dissipation in penalty term
    porous_diss = dx*dy*sum( mask*((u(:,:,1)-maskvx)**2+(u(:,:,2)-maskvy)**2) )
    
    !$omp section  
    ! energy inflow due to penalization
    inflow_solid = dx*dy*sum( mask*( maskvx*(u(:,:,1)-maskvx) + maskvy*(u(:,:,2)-maskvy)  ) )
    
    !$omp section  
    ! energy inflow (note we use the full mask, not only the actual object)
    inflow_mean = dx*dy*u_mean_true*sum( Mean_ux(time+dt1)*penal(:,:,1)+Mean_uy(time+dt1)*penal(:,:,2) )
  !$omp end parallel sections  
  
  open (14, file = trim(dir_name)//'/'//trim(simulation_name)//'energies', status = 'unknown', access = 'append')
  write (14, '(8(es15.8,1x))') time, dt1, vor_rms, kinetic, dissipation, porous_diss, inflow_solid, inflow_mean
  close (14)
  
end subroutine




!===================================================================================================




subroutine pressure ( time, dt1, vortk, press, u, beams )
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
  integer :: iy
  real (kind=pr), intent(in) :: time,dt1
  real (kind=pr), dimension (0:nx-1, 0:ny-1), intent (out) :: press
  real (kind=pr), dimension (0:nx-1, 0:ny-1, 1:2), intent (out) :: u
  real (kind=pr), dimension (0:nx-1, 0:ny-1, 1:2) :: penal
  real (kind=pr), dimension (0:nx-1, 0:ny-1), intent (in)  :: vortk
  type (solid), dimension(1:iBeam), intent(inout) :: beams
  !--Work arrays
  real (kind=pr), dimension (0:nx-1, 0:ny-1) :: work1, work2
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
  call cal_drag ( time, dt1, penal, beams )
  
  !------------------------------------------------
  ! vor to physical space
  ! you need it for the non-linear term in the 
  ! pressure eqn, and you compute vor_rms
  !------------------------------------------------
  call cofitxy (vortk, work1) !-- vort in physical space
  
  !------------------------------------------------
  ! compute energy contributions
  !------------------------------------------------  
  call cal_energy ( time, dt1, work1, u, penal)
  
  
  !------------------------------------------------
  ! for fixed beams (like the CFD test) it is not nessesairy to compute the pressure, so we skip what follows.
  !-----------------------------------------------  
  if ((iFLUSI==1).or.(iSaveBeam>0)) then     
      !------------------------------------------------
      ! right hand side of the pressure poisson eqn. 
      ! note we can compute everything in physical space
      ! before taking the divergence - this is more efficient.
      !------------------------------------------------
      
      !------------------------------------------------
      ! we already have the penalty term, add the NL term:
      !------------------------------------------------      
      !$omp parallel do private(iy)
      do iy=0,ny-1
        !penalty term is positive in poisson-pressure eqn
        penal(:,iy,1) = penal(:,iy,1) - work1(:,iy)*u(:,iy,2) ! -vor*uy
        penal(:,iy,2) = penal(:,iy,2) + work1(:,iy)*u(:,iy,1) ! +vor*ux
      enddo
      !$omp end parallel do      
      
      !------------------------------------------------
      ! compute divergence
      !------------------------------------------------      
      !--Calculate x-derivative in Fourier space
      call coftxy (penal(:,:,1), work2) !-- d Px / dx
      call cofdx  (work2, penal(:,:,1))
      !--Calculate y-derivative in Fourier space
      call coftxy (penal(:,:,2), work2) !-- d Py / dy
      call cofdy  (work2, penal(:,:,2))
      
      !--right-hand side of poisson eqn 
      !$omp parallel do private(iy)
      do iy=0,ny-1
        work1(:,iy) = penal(:,iy,1) + penal(:,iy,2)  !-- divergence of penalty+NL term (fourier-space)
      enddo
      !$omp end parallel do      
      
      !------------------------------------------------
      ! solve poisson eqn
      !------------------------------------------------
      call poisson (work1, press) !--Solve Poisson equation to get the pressure q. q is static pressure p + 0.5*v_abs**2  
      call cofitxy (press, work1) !--Transform pressure to physical space
      
      !------------------------------------------------
      ! substract dynamical pressure to get static pressure
      !------------------------------------------------
      !$omp parallel do private(iy)
      do iy=0,ny-1
        press(:,iy) = work1(:,iy) - 0.5*(u(:,iy,1)**2 + u(:,iy,2)**2)
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

subroutine EvolveFluidExplicit(time, dt1, n0, n1, beams, vortk, workvis, nlk, u, press )
  use share_vars
  use FieldExport
  implicit none
  real (kind=pr), dimension (0:nx-1, 0:ny-1, 0:1), intent (inout) :: vortk, workvis, nlk
  real (kind=pr), dimension (0:nx-1, 0:ny-1), intent (inout) :: press
  real (kind=pr), dimension (0:nx-1, 0:ny-1) :: stream
  real (kind=pr), dimension (0:nx-1, 0:ny-1, 1:2), intent (inout) :: u
  type (solid), dimension(1:iBeam), intent(inout) :: beams
  real (kind=pr), intent (inout) :: dt1, time
  integer, intent (in) :: n0, n1
  integer :: iy
  !---------------------------------------------------------
  !     n0 = n   n1 = n-1
  !	n0 = n   n1 = n+1 (for vortk and workvis)
  !---------------------------------------------------------

  !--for the first time step, compute velocity and penalization term explicitly
  call cal_velocity ( time, vortk(:,:,n0), u, stream )
  !--calculate RHS
  call cal_nlk ( time, dt1, nlk (:,:,n0), vortk(:,:,n0), u, fmean(:,n0) )
  !--compute exp term for viscosity
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
    
  call pressure ( time,dt1, vortk(:,:,n1), press, u, beams )  
end subroutine EvolveFluidExplicit

!==============================================================================

subroutine EvolveFluidAB2(time, dt0, dt1, n0, n1, beams, vortk, workvis, nlk, u, press)
  use share_vars
  use FieldExport
  implicit none
! adams bashforth step
  real (kind=pr), dimension (0:nx-1, 0:ny-1, 0:1), intent (inout) :: vortk, workvis, nlk
  real (kind=pr), dimension (0:nx-1, 0:ny-1, 1:2), intent (inout) :: u
  real (kind=pr), dimension (0:nx-1, 0:ny-1), intent (out) :: press
  type (solid), dimension(1:iBeam), intent(inout) :: beams
  real (kind=pr), intent (inout) :: dt0, dt1, time
  real (kind=pr) ::b10,b11
  integer, intent (in) :: n0, n1
  integer :: iy


  !$omp parallel do private(iy)
  do iy=0,ny-1
      nlk(:,iy,n1) = nlk(:,iy,n1)*workvis(:,iy,n0) !OLD workvis, OLD NLK
  enddo
  !$omp end parallel do

  !--Calculate fourier coeffs of nonlinear rhs and forcing (is RHS_n !! )
  call cal_nlk ( time, dt1, nlk (:,:,n0), vortk(:,:,n0), u, fmean(:,n0) )

  
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
  call pressure(time, dt1, vortk(:,:,n1), press, u, beams )
  
end subroutine EvolveFluidAB2

!==========================================================================================================================

subroutine cal_drag ( time, dt1, penal, beams )
  use share_vars
  use FieldExport
  implicit none
  real (kind=pr), intent (in)                                   :: time, dt1
  real (kind=pr), dimension (0:nx-1, 0:ny-1,1:2), intent (in)   :: penal
  real (kind=pr), dimension (1:6)                               :: LeadingEdge
  type (solid), dimension(1:iBeam), intent(inout)               :: beams
  real (kind=pr)                                                :: dx, dy, norm, alpha, alpha_t, alpha_tt, h_star
  integer                                                       :: ixmin,ixmax,iymin,iymax,i 
  
  dx = xl / real (nx)
  dy = yl / real (ny)
  norm = dx*dy
  
  do i=1, iBeam  
  call mouvement(time, alpha, alpha_t, alpha_tt, LeadingEdge, beams(i) )

  !------------------------------------------------------------
  !-- Compute Drag and Lift forces
  !------------------------------------------------------------
  ! optimized 10/11/2011 now using penal from pressure subroutine and the smallest (best) integration windows possible  
  if (iBeam>0) then ! case of a beam with or without cylinder at the leading edge
    iymax = int( (max( maxval(beams(i)%y)+4.0*t_beam, LeadingEdge(2)+2.0*r_cylinder))/dy)
    iymin = int( (min( minval(beams(i)%y)-4.0*t_beam, LeadingEdge(2)-2.0*r_cylinder))/dy)
    ! absolute upper limit for xmin is the leading edge - 4 times thickness
    ixmin = int( min( (LeadingEdge(1)-max(2.5*r_cylinder, 4.0*t_beam)), (minval(beams(i)%x)-4.0*t_beam) )  /dx)
    ! even if the beam is hanging, its minimum x value is the leading edge
    ixmax = int( (maxval(beams(i)%x)+4.0*t_beam)  /dx)
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

  ! lift / drag:
  beams(i)%Force(1) = sum(penal(ixmin:ixmax,iymin:iymax,1)) * norm
  beams(i)%Force(2) = sum(penal(ixmin:ixmax,iymin:iymax,2)) * norm

  !------------------------------------------------------------
  !-- Compute unsteady lift/drag corrections
  !------------------------------------------------------------
  ! all save variables go to disk; there shouldn't be a problem with restarting
  beams(i)%drag_unst_new = sum (maskvx(ixmin:ixmax,iymin:iymax)*mask(ixmin:ixmax,iymin:iymax)) * norm * eps
  beams(i)%lift_unst_new = sum (maskvy(ixmin:ixmax,iymin:iymax)*mask(ixmin:ixmax,iymin:iymax)) * norm * eps  
  
  !compute unsteady correction
  if (beams(i)%UnsteadyCorrectionsReady) then
  beams(i)%Force_unst(1)   = (beams(i)%drag_unst_new - beams(i)%drag_unst_old)/dt1 
  beams(i)%Force_unst(2)   = (beams(i)%lift_unst_new - beams(i)%lift_unst_old)/dt1  
  else
  beams(i)%Force_unst = 0.0
  endif
  
  beams(i)%drag_unst_old = beams(i)%drag_unst_new  !iterate
  beams(i)%lift_unst_old = beams(i)%lift_unst_new  !iterate
  
  ! now, in any case, we're ready!
  beams(i)%UnsteadyCorrectionsReady = .true.


  enddo ! end loop over beams
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
  real (kind=pr), intent(in) :: time
  real (kind=pr), dimension (0:nx-1, 0:ny-1), intent (out) :: press
  real (kind=pr), dimension (0:nx-1, 0:ny-1, 1:2) :: penal, u
  real (kind=pr), dimension (0:nx-1, 0:ny-1), intent (in)  :: vortk  
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

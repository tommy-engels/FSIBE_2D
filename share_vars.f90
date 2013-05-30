module share_vars
  implicit none
  integer, parameter :: pr = kind (0.d0) 
  

  character (len=40) :: simulation_name, dir_name, simulation_name_org, inicond_file
  integer, save :: nx, ny,ns, iMeanVelocity, iSpongeType, nPalettes=14, iMultiRes
  integer, save :: idealis, ihypvisc, inicond, iobst, iSponge, iBeam, iWalls, iCylinder, iImpulse
  !!!!!!!!!!!
  integer, save :: iFLUSI, iMotion, iViscous
  !!!!!!!!!!!
  ! save switches
  integer, save :: iSaveBeam, iSaveMask, iSaveMaskVel, iSaveSTR, iSavePress, iSaveVel, iSaveVort, iSaveStress
  integer, dimension (10), save :: ifaxx, ifaxy
  real (kind=pr), save :: cfl, nu, eps, pi, scalex, scaley, length, xl, yl, eta, mue, t_beam, N_smooth, h_channel, AngleBeam
  real (kind=pr), save :: U_mean_true, theta_inf, ds, T_release, tau, dt_fixed, sigma, R_cylinder, Time_end, SpongeSize, T_fullspeed
  real (kind=pr), dimension (:), allocatable, save :: trigsx, trigsy, T_beam_save
  real (kind=pr), dimension (:,:), allocatable, save :: dealiase, mask, maskvx, maskvy, mask_sponge, vor_init, vor_init_high
  real (kind=pr), save :: x0, y0, grav, eps_sponge, tdrag, tsave, Time_init  
  
  
  real (kind=pr), dimension(1:2), save :: u_mean
  real (kind=pr), dimension(1:2,0:1), save :: fmean
  
  real (kind=pr), save :: colorscale = 0.0 ! scaling for farge palette (vorticity)
 
  ! ----------------------------------------------
  ! mask parameter
  ! ----------------------------------------------
  integer :: iSharpTrailing 
  logical, save::  sharp= .false. ! if the smoothing is 0, then we switch to a sharp mask.
  
  
  real (kind=pr), save :: time_pardiso, time_A, time_B, time_C
  
  
  real (kind=pr), dimension (:), allocatable, save :: pressure_beam_init
  real (kind=pr), dimension (:,:), allocatable, save :: beam_init, beam_tmp
  
  character (len=256), dimension (:), allocatable, save :: Params_Header
  logical, save:: continue_timestep = .true.

  ! ----------------------------------------------
  ! Solid Solver variables
  ! ----------------------------------------------
  
  integer, parameter :: EulerImplicit = 1, CrankNicholson = 2, RungeKutta4 = 3, EulerExplicit=4, BDF2 = 5 ! integers for different solid solvers. do not modify!!!!!
  integer, save :: TimeMethodSolid != BDF2  ! set values from list above; NOT VALUE IS OVERWRITTEN IN PARAMS
  logical, save :: StartupStep=.true.
  integer, parameter :: N_cpu_solid = 4 ! number of CPU's for PARDISO solver
  
  ! ----------------------------------------------
  ! Solid datatype
  ! ----------------------------------------------
  type Solid
    real(kind=pr), dimension (:), allocatable :: x,y,vx,vy, theta, theta_dot
    real(kind=pr), dimension (:), allocatable :: pressure_old, pressure_new, tau_old, tau_new
    real(kind=pr), dimension(1:2) :: Force, Force_unst, Force_press
    real(kind=pr) :: E_kinetic, E_pot, E_elastic
    real(kind=pr) :: x0, y0
    !MOUVEMENT MUST GO HERE ALSO
  end type

  
end module share_vars

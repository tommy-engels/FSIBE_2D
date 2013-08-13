program dns
  use share_vars
  use FieldExport
  use SolidSolver
  implicit none
  integer :: N_runs
  real(kind=pr) :: dx,dy  
  integer :: i
  
  
  write (*,'(A)') "O--o--|--o--|--o--|--o--|--o--|--o--|--o--|--o--|--o--|--o--|--o--|--o--|--o--|--o--|--o---O" 
  write (*,'(A)')" ______ _____ _____ ____  ______    _____       _                                   __ _____ "
  write (*,'(A)')"|  ____/ ____|_   _|  _ \|  ____|  / ____|     | |  Engels Kolomenskiy Schneider   /_ | ____|"
  write (*,'(A)')"| |__ | (___   | | | |_) | |__    | (___   ___ | |_   _____ _ __              __   _| | |__  "
  write (*,'(A)')"|  __| \___ \  | | |  _ <|  __|    \___ \ / _ \| \ \ / / _ \ '__|             \ \ / / |___ \ "
  write (*,'(A)')"| |    ____) |_| |_| |_) | |____   ____) | (_) | |\ V /  __/ |                 \ V /| |___) |"
  write (*,'(A)')"|_|   |_____/|_____|____/|______| |_____/ \___/|_| \_/ \___|_|                  \_/ |_|____/ "
  write (*,'(A)')"                                                                                             "
  write (*,'(A)') "O--o--|--o--|--o--|--o--|--o--|--o--|--o--|--o--|--o--|--o--|--o--|--o--|--o--|--o--|--o---O"
  
  call params()
  write (*,*) "*** information: passed params routine"
  
  dx = xl / real (nx)
  dy = yl / real (ny)
  

  if (iMultiRes==0) then !just a normal single run
      dir_name='.'
      write (*,*) "*** Information: single run"
      call StartSimulation() 
  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  !@@@     run only solid solver   
  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  elseif (iMultiRes==666) then
      dir_name='.'
      write (*,*) "*** Information: running solid solver only"
      call OnlySolidSimulation()
!   elseif (iMultiRes==777) then
!       dir_name='.'
!       write (*,*) "*** Information: running solid solver only"
!       call SolidTimeConvergence()      
  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  !@@@     galilean change of reference frame    
  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  elseif (iMultiRes==5) then
!       write (*,*) "*** Information: iMultiRes==5 (Gallilean run)"
!       !----------------------------
!       ! fixed beam
!       !------------------------------
!       dir_name='./mvt'
!       simulation_name=trim(simulation_name_org)//'_mvt.'
!       iFSI = 3
!       iMeanVelocity = 4
!       
!       call StartSimulation() 
!       !----------------------------
!       ! moving beam
!       !------------------------------
!       dir_name='./fix'
!       simulation_name=trim(simulation_name_org)//'_fix.'
!       iFSI=0
!       iMeanVelocity=3
!       
!       call StartSimulation()    
write (*,*) "to do, fix me here..:"
  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  !@@@     batch runs, always starting at T=0, not with coarse-to-fine
  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@      
  elseif (iMultiRes==1) then
      write (*,*) "*** Information: iMultiRes==1 (batch run)"
      ! read the different epsilons from the file
      open (78, file = 'batch_resolutions.in', status = 'old')
      read (78,*) N_runs ! first line in file should contain number of elements
      
      do i=1, N_runs ! already did the first one
	  ! read in parameters for the run
	  read (78,*) nx, ny, eps, Time_end, dt_fixed

	  write (*,'(A)') "::::::::::::::::::::::::::::::::::::::::::::::::"
	  write (*,'(" Batch run  ------> ", 2(i4, 1x), 2(es11.5,1x)," dt=",es11.5 )') nx, ny, eps, Time_end, dt_fixed
	  write (*,'(A)') "::::::::::::::::::::::::::::::::::::::::::::::::"

	  ! directory
! 	  write (dir_name       ,'("nx",i5.5,"_ny",i5.5,"_eps",es7.1)') nx, ny, eps
! 	  write (simulation_name,'("nx",i5.5,"_ny",i5.5,"_eps",es7.1)') nx, ny, eps
	  
! 	  write (dir_name       ,'(i2.2,"_dt",es7.1)') i, dt_fixed
! 	  write (simulation_name,'(i2.2,"_dt",es7.1)') i, dt_fixed
	  
	  write (dir_name       ,'(i2.2,"_eps",es7.1)') i, eps
	  write (simulation_name,'(i2.2,"_eps",es7.1)') i, eps
	  
	  
	  simulation_name = trim(simulation_name_org) // trim(simulation_name) // '.'
	  write (*,*) dir_name
	  write (*,*) simulation_name
	  
	  call StartSimulation() 	  
      enddo
      close (78)   

  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  !@@@     turek tests (batch runs, higher resolution starts where the lower one finnished)
  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@      
  elseif (iMultiRes==2) then
      write (*,*) "*** Information: iMultiRes==2 CURRENTLY UNAVAILABLE!!!!!!"
!       ! read the different epsilons from the file
!       open (78, file = 'batch_resolutions.in', status = 'old')
!       read (78,*) N_runs ! first line in file should contain number of elements
!       
!       ! --------------------------------------------------------------------------------------
!       ! first run on coarsest grid
!       read (78,*) nx, ny, eps, Time_end
!       write (*,'(A)') "::::::::::::::::::::::::::::::::::::::::::::::::"
!       write (*,'(" First Batch run  ------> ", 2(i4, 1x), 2(es11.5,1x) )') nx, ny, eps, Time_end
!       write (*,'(A)') "::::::::::::::::::::::::::::::::::::::::::::::::"
!       allocate ( vor_init(0:nx-1, 0:ny-1) )  
! 
!       ! directory
!       write (dir_name       ,'("nx",i5.5,"_ny",i5.5,"_eps",es7.1)') nx, ny, eps
!       write (simulation_name,'("nx",i5.5,"_ny",i5.5,"_eps",es7.1)') nx, ny, eps
!       simulation_name = trim(simulation_name_org) // trim(simulation_name) // '.'
!       write (*,*) dir_name
!       write (*,*) simulation_name
! 
!       call StartSimulation() 
!       ! --------------------------------------------------------------------------------------
! 
!       inicond = 99 ! right initial condition
!       T_release = 0.0
! 
!       do i=2, N_runs ! already did the first one
! 	  ! read in parameters for the run
! 	  read (78,*) nx, ny, eps, Time_end
! 	  Time_end = Time_end + Time_init
! 	  write (*,'(A)') "::::::::::::::::::::::::::::::::::::::::::::::::"
! 	  write (*,'(" Batch run  ------> ", 2(i4, 1x), 2(es11.5,1x) )') nx, ny, eps, Time_end
! 	  write (*,'(A)') "::::::::::::::::::::::::::::::::::::::::::::::::"
! 
! 	  ! resize the final field on the coarser grid
! 	  allocate   ( vor_init_high(0:nx-1,0:ny-1) )
! 	  call ResizeField (vor_init, vor_init_high,xl,yl)
! 	  call SaveGif (vor_init_high, trim(dir_name)//"/initial_condition_next_resolution")
! 	  deallocate ( vor_init )
!           allocate   ( vor_init (0:nx-1,0:ny-1) )
! 	  vor_init = vor_init_high
! 	  deallocate ( vor_init_high )
! 
! 	  ! directory
! 	  write (dir_name       ,'("nx",i5.5,"_ny",i5.5,"_eps",es7.1)') nx, ny, eps
! 	  write (simulation_name,'("nx",i5.5,"_ny",i5.5,"_eps",es7.1)') nx, ny, eps
! 	  simulation_name = trim(simulation_name_org) // trim(simulation_name) // '.'
! 	  write (*,*) dir_name
! 	  write (*,*) simulation_name
! 
! 	  call StartSimulation() 	  
!       enddo
!       close (78)   
	  
  endif  


  call Goodbye()

end program dns



subroutine StartSimulation()
  use share_vars
  implicit none
  write (*,*) "*** information: entering StartSimulation"

  call system('mkdir -p '//trim(dir_name) )
  call system('mkdir -p '//trim(dir_name)//'/fields' )
  call system('mkdir -p '//trim(dir_name)//'/vor' )
  call system('mkdir -p '//trim(dir_name)//'/mask' )
  call system('mkdir -p '//trim(dir_name)//'/press' )

  write (*,*) "*** information: created directories"

  allocate ( dealiase(0:nx-1,0:ny-1) )
  allocate ( mask(0:nx-1,0:ny-1) )    
  allocate ( maskvx(0:nx-1,0:ny-1) )
  allocate ( maskvy(0:nx-1,0:ny-1) )
  allocate ( mask_sponge(0:nx-1,0:ny-1) )
  write (*,*) "*** information: allocated memory"

! Initialize fft
  call fft_initialize
  write (*,*) "*** information: did fft_initialize"
! Set up mask for dealiasing
  dealiase = 1.0d0
  if (idealis == 1) call dealiase_mask
! Step forward in time
  write (*,*) "*** information: entering time_step"
  call time_step  !after time_step, the last vort field is saved in mask_sponge
  deallocate (dealiase, mask_sponge, maskvx, maskvy, Params_Header, mask) ! its important to free all these fields so the solver won't be confused
  call fft_free        

end subroutine StartSimulation


! ! ! ! subroutine hinge_debugging()
! ! ! ! ! I used this 29.05.2013 at 23:40 to fin a bug in create_mask. its found, all smooth.
! ! ! !   use share_vars
! ! ! !   use FieldExport
! ! ! !   use SolidSolver
! ! ! !   implicit none
! ! ! !   real (kind=pr) :: alpha, beta, time,dt
! ! ! !   integer :: ivideo
! ! ! !   real (kind=pr), dimension (0:ns-1, 1:6) :: beam
! ! ! !   character(len=17) :: namestring
! ! ! !   
! ! ! !   ivideo = 0
! ! ! !   allocate ( mask(0:nx-1,0:ny-1)   )      
! ! ! !   allocate ( maskvx(0:nx-1,0:ny-1) )
! ! ! !   allocate ( maskvy(0:nx-1,0:ny-1) )
! ! ! !   
! ! ! !   dt = 5.e-2
! ! ! !   time = 0.0
! ! ! !   
! ! ! !   do while (time<1.)
! ! ! !   alpha = 2.*pi*sin(time)
! ! ! !   beta  = 2.*pi*cos(1.5*time)
! ! ! !   beam = 0
! ! ! !   beam(0,5) = alpha
! ! ! !   beam(1:2,5) = beta
! ! ! !   write (namestring,'(i4.4)')    ivideo
! ! ! !   call integrate_position(0.0, beam)
! ! ! !   call create_mask(0.0, beam)
! ! ! !   call SaveGIF(mask, trim(namestring), 13, 0.0, 1.0/eps)
! ! ! !   ivideo = ivideo + 1
! ! ! !   time = time +dt
! ! ! !   enddo
! ! ! !   
! ! ! !   
! ! ! !   deallocate (maskvx, maskvy, mask)
! ! ! ! 
! ! ! ! 
! ! ! ! end subroutine 


! ********************************************************************************************************************************************
! **
! **	A single run for the solid model, without surrounding fluid
! **
! ********************************************************************************************************************************************

subroutine OnlySolidSimulation()
  use share_vars
  use SolidSolver
  use PerformanceMeasurement
  implicit none
  type(solid), dimension(1:iBeam) :: beams
  real (kind=pr) :: time
  integer :: it,ndrag,nsave

  write (*,*) "*** information: starting OnlySolidSimulation"
  call system('mkdir '//trim(dir_name) )
  time = GetRuntime('start') !for runtime measurement
  time = 0.0
  it = 0
  
  continue_timestep = .true.
  
  if ((dt_fixed)==0.0) dt_fixed = 1.0e-4 ! in case you forgot to set it

  ! initialization
  call InitBeamFiles()	!init files for data, override
  call init_beams ( beams )
  
  ndrag = int(tdrag/dt_fixed)
  nsave = max(int(1.0e-4/dt_fixed),1)


  do while ((time<Time_end).and.(continue_timestep==.true.))  
    !---------------------
    ! actual time step
    !---------------------
    call SolidSolverWrapper( time, dt_fixed , beams )
    
    
    it   = it+1
    time = real(it)*dt_fixed


    if (mod(it,nsave)==0 ) call SaveBeamData( time, beams, dt_fixed )
    
    if ( (mod(it,ndrag)==0).or.(it==22) ) then
      call SaveDeflectionLine (time, beams )
    endif
  enddo

end subroutine



! ********************************************************************************************************************************************
! **
! **	A convergence test for the solid model (dt varies)
! **
! ********************************************************************************************************************************************



! subroutine SolidTimeConvergence()
!   use share_vars
!   use PerformanceMeasurement
!   use SolidSolver
!   implicit none
!   real (kind=pr), dimension (0:ns-1, 1:6) :: beam, beam_ref
!   real (kind=pr), dimension (0:ns-1) :: pressure, tau_beam
!   real (kind=pr), dimension (1:9) :: dts
!   real (kind=pr) :: time, time_dt, T_lastdrag
!   integer (kind=8) :: it=0,i,nt
! 
!   write (*,*) "*** information: starting SolidTimeConvergence"
!   x0=0.0
!   y0=0.0
!   
!   time = GetRuntime('start') !for runtime measurement
!   
!   dts(1) = 1.0e-3
!   dts(2) = 5.0e-4
!   dts(3) = 2.5e-4
!   dts(4) = 2.0e-4
!   dts(5) = 1.0e-4
!   dts(6) = 5.0e-5
!   dts(7) = 2.5e-5
!   dts(8) = 2.0e-5  
!   dts(9) = 1.0e-5  
! 
!   ! initialization
!   call InitBeamFiles()	!init files for data, override
!   
!   call InitializeSolidSolver() ! to tell the BDF2 solver to use CN2 in the first step
!   
!   tau_beam = 0.0
!   pressure = 0.0
!   
!   ! --------------------------
!   !	reference solution
!   ! --------------------------
!   write (*,*) "*** information: starting Computing reference solution"
!   dt_fixed = 5.0e-6
!   time = 0.0
!   nt=nint(Time_end/dt_fixed)
!   call init_beam (beam_ref, pressure)
!   do it=1,nt
!     time_dt=performance("start",10)    
! 
!     call SolidSolverWrapper( time, dt_fixed , beam_ref, pressure, pressure, tau_beam, tau_beam)
!     time=real(it)*dt_fixed
!     
!     time_dt=performance("stop",10)    
!     if ( (time-T_lastdrag>tdrag).or.(it==22) ) then
!       T_lastdrag=time
!       call SavePerformance(time, time_dt, 0.0, 0.0, 0.0, 0.0, nint((Time_end-time)/dt_fixed), dt_fixed )
!     endif
!   enddo
! 
!   write(*,*) "reference solution!!!!", dt_fixed, time
!   do it=1,6
!   write(*,'(1024(es15.8,1x))') beam_ref(:,it)
!   enddo
!   write(*,*) "---"
! 
!   
!   ! --------------------------
!   !	loop over dt's
!   ! --------------------------
!   write (*,*) "*** information: beginning loop over dt"
!   do i=1,9
!       dt_fixed = dts(i); time = 0.0; T_lastdrag = 0.0
!       nt=nint(Time_end/dt_fixed)
!       
!       call init_beam (beam, pressure)
!       call InitializeSolidSolver() ! to tell the BDF2 solver to use CN2 in the first step
!       
!       do it=1,nt
! 	if (continue_timestep) then
! 	    time_dt=performance("start",10)    
! 	    call SolidSolverWrapper( time, dt_fixed , beam, pressure, pressure, tau_beam, tau_beam)
! 	    time=real(it)*dt_fixed
!     
! 	    time_dt=performance("stop",10)    
! 	    if ( (time-T_lastdrag>tdrag).or.(it==22) ) then
! 	      T_lastdrag=time
! 	      call SavePerformance(time, time_dt, 0.0, 0.0, 0.0, 0.0, nint((Time_end-time)/dt_fixed), dt_fixed )
! 	    endif
! 	endif
!       enddo
! 
!       write(*,*) dt_fixed, time
!       do it=1,6
!       write(*,'(1024(es15.8,1x))') beam(:,it)
!       enddo
!       write(*,*) "---"
!       if (continue_timestep) then
!       open  (10, file = "convergence", status = 'unknown', access = 'append')
!       write (10,'(12(1x,es15.8,1x))') time, dt_fixed, &
! 	  sqrt ( sum( beam_ref(:,5)-beam(:,5) )**2 ), &
! 	  maxval (abs(beam_ref(:,5)-beam(:,5) ) ),&
! 	  sqrt ( sum( beam_ref(:,6)-beam(:,6) )**2 ), &
! 	  maxval (abs(beam_ref(:,6)-beam(:,6) ) )
!       close (10) 
!       endif
!       continue_timestep = .true.
!   enddo
!   
!   
! 
! end subroutine
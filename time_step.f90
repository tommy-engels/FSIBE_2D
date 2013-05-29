subroutine time_step 
  use share_vars
  use FieldExport
  use SpectralSolver
  use BeamForces
  use PerformanceMeasurement

  implicit none
  integer :: inter, it,iStop, ivideo=1
  integer :: n0 = 0, n1 = 1
  integer :: q
  integer :: nbackup = 0 ! 0 - backup to file runtime1_backup0, 1 - to runtime1_backup1, 2 - no backup
  real (kind=pr), dimension (0:nx-1, 0:ny-1) :: vort, press
  real (kind=pr), dimension (0:nx-1, 0:ny-1, 0:1) :: nlk, vortk, workvis
  real (kind=pr), dimension (0:nx-1, 0:ny-1, 1:2) :: u
  real (kind=pr), dimension (1:2) :: force_pressure=0.0
  !Forces: 1=drag 2=lift 3=drag_unst 4=lift_unst
  real (kind=pr), dimension (1:4) :: Forces_old=0.0, Forces_new=0.0 
  ! 1= vor_rms 2=vor_rms_dot 3=Fluid kinetic Enegry 4=Dissipation 5=Dissipation(Penalty) 6=EnergyInflow
  real (kind=pr), dimension (1:7) :: FluidIntegralQuantities 
  ! Beam indices: 1=beam_x 2=beam_y 3=beam_vx 4=beam_vy 5=theta 6=theta_dot
  real (kind=pr), dimension (0:ns-1, 1:6) :: beam
  real (kind=pr), dimension (0:ns-1) :: bpressure_old, bpressure_new
  real (kind=pr), dimension (0:ns-1) :: tau_beam_new, tau_beam_old
  integer*8 , dimension (1:2)     ::  total_time,dummy2 ! for performance measurement
  integer*8			  ::  systemtime_lastbackup
  real (kind=pr) :: dt0 = 0, dt1, ROC, time_start
  real (kind=pr) :: time, pmin, pmax
  real (kind=pr) :: dx, dy
  real (kind=pr) :: T_lastsave, T_lastdrag
  real (kind=pr) :: time_dt, time_nst, time_solid, time_mask, time_pressure ! times for performance measurement
  character(len=17) :: namestring, timestring
  logical :: file_existent, Colorscale_done = .false.


  write(*,'(A)') "=================================================<( time_step.f90 )>==================================================="
 
  time = GetRuntime('start') !for runtime1 measurement
  time = 0.0

    
  if (inicond.ne.2) then
    T_lastsave = 0.0
    T_lastdrag = 0.0
    !-- if imposing zero mean pressure, initialize u_mean
    u_mean = 0.0
  endif
  
  
  q=0;dt1=0.;ROC=0.
  dx = xl/real(nx)
  dy = yl/real(ny)



  call CreateHeader() !stores the PARAMS file in a array in order to write it to the output files as a header
  write(*,*) "*** information: allocated and created header for output files"
  
  !----------------------------------------------------------------
  ! Initialize output files and beam
  !----------------------------------------------------------------  
  if (inicond.ne.2) then 	!not retaking a backup
    call InitBeamFiles()	!init files for data, override
    if (iBeam>0) call init_beam (beam, bpressure_old)
  endif
  
  !mark in performance file that you restarted or started now
  if (inicond==2) then
  open  (10, file = trim(dir_name)//'/'//trim(simulation_name)//"performance_details", status = 'unknown', access = 'append')
  write (10,'(A)' ) "----------restarting!"
  close (10)
  endif
  
  !----------------------------------------------------------------
  ! Initialize the solid solver
  !----------------------------------------------------------------  
  call InitializeSolidSolver() ! to tell the BDF2 solver to use CN2 in the first step
  
  !----------------------------------------------------------------  
  ! Initialize vorticity or read values from a backup file
  ! init must be called BEFORE create_sponge_mask when running multiple-resolutions!!!!
  !----------------------------------------------------------------
  call init_fields (n1, time, dt1, vortk, nlk, workvis, beam, bpressure_old, tau_beam_old, ivideo, u)
  n0 = 1 - n1
  
  !----------------------------------------------------------------
  ! create startup mask
  !---------------------------------------------------------------- 
  call create_mask (time, beam)
  call SaveGIF(mask, trim(dir_name)//'/'//trim(simulation_name)//"startup_mask", 13, 0.0, 1.0/eps)
  
  !----------------------------------------------------------------
  ! vorticity sponge
  !---------------------------------------------------------------- 
  call create_sponge_mask() !contains saving the field  
  if ((iSponge>0).and.(iSponge<4)) call SaveGIF(mask_sponge, trim(dir_name)//'/'//trim(simulation_name)//"mask_sponge", 13, minval(mask_sponge), maxval(mask_sponge))
  
  !----------------------------------------------------------------
  ! compute true mean speed for channels
  !----------------------------------------------------------------
  U_mean_true = 1.d0
  if ((iWalls==1)) then
    U_mean_true=1.0-2.0*h_channel/yl
    write (*,*) "--- True mean speed:", U_mean_true
  endif
  
  !---------------------------------------------------------
  ! Save initial values (fields)
  !---------------------------------------------------------
  nbackup = 2 ! no backup at the beginning
  if (inicond .ne. 2 ) then
  call save_fields (n1, time, dt1, vortk, nlk, workvis, nbackup, beam, bpressure_old, ivideo, u, press, tau_beam_old)
  write (*,*) "*** information: saved startup fields"
  endif
  
  
  call SaveBeamData(time,beam,bpressure_new*0.0, dt1,q,ROC,Forces_new*0.0,force_pressure*0.0, FluidIntegralQuantities*0.0,it)

  !------------------------------------------------------------------------------------------------------------------
  !     STARTUP - EULER EXPLICIT
  !------------------------------------------------------------------------------------------------------------------
  write (*,*) "*** information: Initialization done. Ready for startup."
  
  if (inicond .ne. 2) then ! inicond = 2 means to retake a backup, no startup nessesaire
      !------------------------------------------------------------------------------------------------------------------
      ! step one: create_mask. nessesairy if the beam moves and if there is a velocity sponge with time-dependend mean velocity
      !------------------------------------------------------------------------------------------------------------------
      if (((iSponge==4).and.(time<T_fullspeed)).or.(iFLUSI==1).or.(iMotion>0)) call create_mask( time, beam )
      !--------------------------------------------------------------------------------------------------------------------
      ! step two: solve navier-stokes (always)
      !--------------------------------------------------------------------------------------------------------------------
      call EvolveFluidExplicit(time, dt1, n0, n1, beam, vortk, workvis, nlk, press, u, Forces_new, FluidIntegralQuantities)
      !--------------------------------------------------------------------------------------------------------------------
      ! step three: get forces. required for FSI runs or if you wish to save the forces
      !--------------------------------------------------------------------------------------------------------------------
      if ((iFLUSI==1).or.(iSaveBeam>0)) call GetForces  ( time, beam, bpressure_new, tau_beam_new, press, force_pressure,u )
      !--------------------------------------------------------------------------------------------------------------------
      ! step four: get the new beam. either by imposed motion of by solving the beam eqn
      !--------------------------------------------------------------------------------------------------------------------
      if (iFLUSI==1) then
	  !-----------------------------
	  ! active FSI coupling
	  !-----------------------------
	  call GravityImpulse(time)
	  call SolidSolverWrapper( time, dt1 , beam, bpressure_old, bpressure_new, tau_beam_old, tau_beam_new) 	
      elseif ((iFLUSI==0).and.(iMotion>0)) then
	  !-----------------------------
	  ! imposed motion only
	  !-----------------------------
	  call integrate_position (time, beam)
      endif
     

      inter 		= n1 ; n1 = n0 ; n0 = inter 			! Switch time levels
      time 		= time + dt1  					! Advance in time
      bpressure_old 	= bpressure_new			!--Advance beam forces
      tau_beam_old  	= tau_beam_new
      call SaveBeamData(time,beam,bpressure_new, dt1,q,ROC,Forces_new,force_pressure, FluidIntegralQuantities,it) 	!--bpressure new is already at the new time step, ok to be here
  endif

  Forces_old = Forces_new

  write(*,*) "*** information: Startup done, beginning loop over time steps"
  !------------------------------------------------------------------------------------------------------------------
  !     BEGIN LOOP OVER TIME STEPS
  !     n0 = n   n1 = n-1
  !	    n0 = n   n1 = n+1 (for vortk and workvis)
  !------------------------------------------------------------------------------------------------------------------
  call system_clock(total_time(1),dummy2(1),dummy2(2))  !start date of the loop
  systemtime_lastbackup = total_time(1)			!when did I do the last backup? now!

  it = 0
  time_start = time
  continue_timestep=.true.
  
  do while ((time<Time_end).and.(continue_timestep == .true.))
	time_dt = performance("start",10)
	!--Advance the time step
	dt0 = dt1
	
	!------------------------------------------------------------------------------------------------------------------
	! step one: create_mask. nessesairy if the beam moves and if there is a velocity sponge with time-dependend mean velocity
	!------------------------------------------------------------------------------------------------------------------	    
	time_mask = performance("start",2)
	if (((iSponge==4).and.(time<T_fullspeed)).or.(iFLUSI==1).or.(iMotion>0)) call create_mask( time, beam )
	time_mask = performance("stop",2)
	!--------------------------------------------------------------------------------------------------------------------
	! step two: solve navier-stokes (always)
	!--------------------------------------------------------------------------------------------------------------------
	time_nst = performance("start",1)
	call EvolveFluidAB2(time, dt0, dt1, n0, n1, beam, vortk, workvis, nlk, u, press, Forces_new, FluidIntegralQuantities) !u: in/out
	time_nst = performance("stop", 1)
	!--------------------------------------------------------------------------------------------------------------------
	! step three: get forces. required for FSI runs or if you wish to save the forces
	!--------------------------------------------------------------------------------------------------------------------
	time_pressure=performance("start",3)
	if ((iFLUSI==1).or.(iSaveBeam>0)) call GetForces  ( time, beam, bpressure_new, tau_beam_new,  press, force_pressure,u)
	time_pressure=performance("stop",3)	    
	!--------------------------------------------------------------------------------------------------------------------
	! step four: get the new beam. either by imposed motion or by solving the beam eqn
	!--------------------------------------------------------------------------------------------------------------------
	time_solid=performance("start",4)
	if (iFLUSI==1) then
	    !-----------------------------
	    ! active FSI coupling
	    !-----------------------------
	    call GravityImpulse(time)
	    call SolidSolverWrapper( time, dt1 , beam, bpressure_old, bpressure_new, tau_beam_old, tau_beam_new) 	
	elseif ((iFLUSI==0).and.(iMotion>0)) then
	    !-----------------------------
	    ! imposed motion only
	    !-----------------------------
	    call integrate_position (time, beam)
	endif
	time_solid=performance("stop",4)   

  
	inter 		= n1 ; n1 = n0 ; n0 = inter 		!--Switch time levels
	time 		= time + dt1 				!--Advance in time
	bpressure_old 	= bpressure_new			!--Advance beam forces
	tau_beam_old 	= tau_beam_new
	time_dt 	= performance("stop",10)
	
	!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
	!--output
	!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  
	call SaveBeamData(time,beam,bpressure_new,dt1,q,ROC,Forces_new,force_pressure, FluidIntegralQuantities,it) 	!--bpressure new is already at the new time step, ok to be here
	if (it==22) then !this gives the first performance estimation
	  call SavePerformance(time, time_dt, time_nst, time_pressure, time_mask, time_solid, nint((Time_end-time)/dt1), dt1 )
	endif

 
	if ((time-T_lastdrag>tdrag).or.(continue_timestep==.false.)) then ! this block is also executed when finnishing
	    !--------------------------------------------------------------------------------------------
	    !--		Backuping after 4 hours
	    !--------------------------------------------------------------------------------------------
	    call system_clock(total_time(2),dummy2(1),dummy2(2))
 	    if ((real(total_time(2)-systemtime_lastbackup)/real(dummy2(1)))>4.0*60.0*60.0) then
	      write(*,'(" +++ Its time for a backup! it=",i7," time=",es11.4," elapsed time since last backup= ",es11.4)') it, time, real(total_time(2)-systemtime_lastbackup)/real(dummy2(1))
	      call MakeRuntimeBackup (n1, time, dt1, vortk, nlk, workvis, nbackup, beam, bpressure_old, ivideo, u, tau_beam_old)
	      systemtime_lastbackup = total_time(2)
	    end if
	    !--------------------------------------------------------------------------------------------
	    !--		video snapshots
	    !--------------------------------------------------------------------------------------------
	    write (namestring,'(i4.4)')    ivideo
	    write (timestring, '(es10.4)') time
	    call cofitxy (vortk(:,:,n1), vort)
	    
	    ! remove vorticity from inside the solid, making the video nicer to watch
	    where (mask*eps>0.5)
	      vort = 0.0
	    end where


	    if ((iFLUSI==1).or.(iSaveBeam>0)) then 
	    ! in these runs, the pressure is computed in every time step, so we can just save it
	    else ! in other runs, we don't need the pressure, so we have to compute it here to make a snapshot
	      call ComputePressureSnapshot(time, vortk(:,:,n1), press)
	    endif

	    if (((time<T_fullspeed+6.0*tdrag).or.(Colorscale_done==.false.)).or.(inicond==22)) then ! the color scaling shouldnt alter after T_fullspeed+6*0.05 (sixth snapshot after fullspeed)
	      colorscale=0.25*max(maxval(vort),abs(minval(vort)))
	      pmin = minval(press)
	      pmax = maxval(press)
	      Colorscale_done = .true.
	    endif
	    
	    colorscale = max(colorscale, 0.10*max(maxval(vort),abs(minval(vort))) )
	    write (*,'("time=",es12.4," color=",es12.4)') time, colorscale 
	    ! with fixed scaling
	    call SaveGIF(vort, trim(dir_name)//"/vor/"//trim(namestring)//"_T="//trim(timestring)//".vor", 1, -colorscale, colorscale)
	    call SaveGIF(press, trim(dir_name)//"/press/"//trim(namestring)//"_T="//trim(timestring)//".press", 14, pmin, pmax)
! 	    call SaveGIF(press, trim(dir_name)//"/press2/"//trim(namestring)//"_T="//trim(timestring)//".press")	    
! 	    call SaveGIF(vort, trim(dir_name)//"/vor2/"//trim(namestring)//"_T="//trim(timestring)//".vor", 13, -colorscale, colorscale)	    
	    call SaveGIF(mask, trim(dir_name)//"/press2/"//trim(namestring)//"_T="//trim(timestring)//".mask", 13, 0.0, 1.0/eps)
! 	    call SaveGIF(abs(maskvx), trim(dir_name)//"/press2/vx"//trim(namestring)//"_T="//trim(timestring)//".mask", 13, 0.0, maxval(abs(maskvx)))
! 	    call SaveGIF(abs(maskvy), trim(dir_name)//"/press2/vy"//trim(namestring)//"_T="//trim(timestring)//".mask", 13, 0.0, maxval(abs(maskvy)))
	    ! without fixed scaling
! 	    call SaveGIF(vort, trim(dir_name)//"/vor/"//trim(namestring)//"_T="//trim(timestring)//".vor", 1, minval(vort), maxval(vort))
! 	    call SaveGIF(press, trim(dir_name)//"/press/"//trim(namestring)//"_T="//trim(timestring)//".press", 14, minval(press), maxval(press))

      
	    ivideo=ivideo+1
	    T_lastdrag=time

	    !--------------------------------------------------------------------------------------------
	    !--		save complete deflection line
	    !--------------------------------------------------------------------------------------------
	    if ((iBeam>0)) call SaveBeamPositions(time, beam)
	    if (it>22) then
	      call SavePerformance(time, time_dt, time_nst, time_pressure, time_mask, time_solid, nint((Time_end-time)/dt1), dt1 )
	    endif
	    endif
	    !--------------------------------------------------------------------------------------------
	    !--   save fields
	    !--------------------------------------------------------------------------------------------
	    if (time-T_lastsave>tsave) then
		if ((iFLUSI==1).or.(iSaveBeam>0)) then 
		! in these runs, the pressure is computed in every time step, so we can just save it
		else ! in other runs, we don't need the pressure, so we have to compute it here to make a snapshot
		  call ComputePressureSnapshot(time, vortk(:,:,n1), press)
		endif
		call save_fields (n1, time, dt1, vortk, nlk, workvis, nbackup, beam, bpressure_old, ivideo, u, press, tau_beam_old)
		! if saving a field, we save also the deflection line for postprocessing.
		if ((iBeam>0)) call SaveBeamPositions(time, beam)
		T_lastsave=time
	    endif  
	    !-------------------------------------------------------------------------------------------------------------
	    ! 	runtime remote control
	    !-------------------------------------------------------------------------------------------------------------
	    if (modulo(it,25)==0) then ! check every 25 time steps if you should do something
	    inquire(file = trim(dir_name)//'/STOP_SIMULATION',exist=file_existent)
	    if (file_existent==.true.) then
	      open  (10, file = trim(dir_name)//'/STOP_SIMULATION', status = 'old')
	      read (10, *) iStop
	      close (10)
	      if (iStop == 99) then
		call save_fields (n1, time, dt1, vortk, nlk, workvis, nbackup, beam, bpressure_old, ivideo, u, press, tau_beam_old)
		T_lastsave=time
		write (*,*) "Stopp command recieved. Time=", time, "it=", it
		call ReWriteStopFile
		stop
	      elseif (iStop == 33) then
		call ReLoadParams
		call ReWriteStopFile
	      elseif (iStop == 666) then
		call ReWriteStopFile
		! abort this run (when doing multires) and got to the next
		write (*,*) " !!! interuption command, stopping this run now. whatever."
		continue_timestep=.false.
	      elseif (iStop == 55) then
		call save_fields (n1, time, dt1, vortk, nlk, workvis, nbackup, beam, bpressure_old, ivideo, u, press, tau_beam_old)
		call ReWriteStopFile
	      endif
	    else
	      call ReWriteStopFile
	    endif
  endif
  
  it=it+1

  enddo

  !------------------save before exiting
  call save_fields (n1, time, dt1, vortk, nlk, workvis, nbackup, beam, bpressure_old, ivideo, u, press, tau_beam_old)
  write (*,*) "*** information: Simulation ended without errors. it=", it, " time=", time," . Did backup before exiting"


  !----- save some things you need if you restart this simulation with a better resolution
  if (iMultiRes==2) then
    call cofitxy (vortk(:,:,n1), vor_init)
    beam_init = beam ! for multires FSI
    pressure_beam_init = bpressure_new
    time_init = time
  endif

 
  !---- save a the final vorticity field
  call cofitxy(vortk,press)
  call SaveField( trim(simulation_name) //"vor", press, 1, xl, yl, "precision")

  open (10, file=trim(simulation_name)//"inicond", form='unformatted', status='replace')
  write (10) nx, ny, press
  close (10)
  
end subroutine time_step






subroutine ReWriteStopFile
use share_vars
          open  (10, file = trim(dir_name)//'/STOP_SIMULATION', status = 'replace')
          write (10, *) "00 = iStop. Set 99 to interupt the simulation at the next drag interval."
          write (10, *) "Set 33 to partially reload the PARAMS.m file or 55 to save now or 666 to skip this run and go to the next (if this is a batch run) !"
          close (10)
end subroutine 












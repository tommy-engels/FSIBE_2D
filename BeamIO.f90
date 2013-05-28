! ================================================================================================================================================
! 							Beam Input/Output subroutines
! ================================================================================================================================================



subroutine InitBeamFiles()
  use share_vars
  implicit none
  !-------------------------------
  if ((iSaveBeam==2).or.(iSaveBeam==3).or.(iSaveBeam==1)) then
      open (14, file = trim(dir_name)//'/'//trim(simulation_name)//'beam_data', status = 'replace')
      call WriteHeaderToFile(14)
      write (14,'(A)') "%  beam data file, all quantities for the end point "
      write (14,'(A)') "%         time   pressure(ns-1)  pressure(3ns/4) pressure(ns/2)  x-coordinate    y-coordinate    x-velocity      y-velocity      theta           theta_dot       time step       ROC            iter  drag            lift            drag_unst       lift_unst       F_press_x       F_press_y       vor_rms         vor_rms_dot     E_kin_fluid     Dissipation     DissipPenalty   Energ_solid_in  Energ_meanf_in  E_kin_solid     E_elastic_sol   E_potential     Mask Volume"
      close (14)
      open  (14, file = trim(dir_name)//'/'//trim(simulation_name)//'beam.complete', status = 'replace')
      call WriteHeaderToFile(14)
      close (14)
  elseif (iSaveBeam==1) then
      open (14, file = trim(dir_name)//'/'//trim(simulation_name)//'beam_x', status = 'replace')
      call WriteHeaderToFile(14)
      close (14)
      open (14, file = trim(dir_name)//'/'//trim(simulation_name)//'beam_y', status = 'replace')
      call WriteHeaderToFile(14)
      close (14)
      open (14, file = trim(dir_name)//'/'//trim(simulation_name)//'iter', status = 'replace')
      call WriteHeaderToFile(14)
      close (14)
      open (14, file = trim(dir_name)//'/'//trim(simulation_name)//'beam_vx', status = 'replace')
      call WriteHeaderToFile(14)
      close (14)
      open (14, file = trim(dir_name)//'/'//trim(simulation_name)//'beam_vy', status = 'replace')
      call WriteHeaderToFile(14)
      close (14)
      open (14, file = trim(dir_name)//'/'//trim(simulation_name)//'beam_pressure', status = 'replace')
      call WriteHeaderToFile(14)
      close (14)
      open (14, file = trim(dir_name)//'/'//trim(simulation_name)//'beam_pressure_s', status = 'replace')
      call WriteHeaderToFile(14)
      close (14)
      open (14, file = trim(dir_name)//'/'//trim(simulation_name)//'beam_theta', status = 'replace')
      call WriteHeaderToFile(14)
      close (14)
      open (14, file = trim(dir_name)//'/'//trim(simulation_name)//'beam_theta_dot', status = 'replace')
      call WriteHeaderToFile(14)
      close (14)
  elseif (iSaveBeam==0) then
      open (14, file = trim(dir_name)//'/'//trim(simulation_name)//'beam_data', status = 'replace')
      call WriteHeaderToFile(14)
      write (14,'(A)') "%  beam data file, all quantities for the end point "
      write (14,'(A)') " %         time  dt              drag            lift            drag_unst       lift_unst       F_press_x       F_press_y       vor_rms         vor_rms_dot     E_kin_fluid     Dissipation     DissipPenalty   Energ_solid_in  Energ_meanf_in"
      close (14)
  endif
  
  !----------------------------------------------------------------------------------------------------------------------------
  if (iFSI==8) then
      open (14, file = trim(dir_name)//'/'//trim(simulation_name)//'cylinder', status = 'replace')
      call WriteHeaderToFile(14)
      write (14,'(A)') " %         time     timestep       position x      position y      velocity x      velocity y    acceleration x  acceleration y     force x        force y        force_unst x    force_unst y"
      close (14)
  endif
  !----------------------------------------------------------------------------------------------------------------------------
  
!   if (iFSI==2) then
      open (14, file = trim(dir_name)//'/'//trim(simulation_name)//'mouvement', status = 'replace')
      call WriteHeaderToFile(14)
      write (14,'(A)') " %         time  alpha           alpha_t         alpha_tt        x-pos           y-pos           x-vel          y-vel            x-accel        y-accel"
      close (14)    
!   endif 
    !----------------------------------------------------------------------------------------------------------------------------
    open  (10, file = trim(dir_name)//'/STOP_SIMULATION', status = 'replace')
    write (10, *) "00 = iStop. Set 99 to interupt the simulation at the next drag interval."
    write (10, *) "Set 33 to partially reload the PARAMS.m file or 55 to save now!"
    close (10)
    !----------------------------------------------------------------------------------------------------------------------------
    open  (10, file = trim(dir_name)//'/'//trim(simulation_name)//"performance_details", status = 'replace')
    write (10,'(A)') "%    time_left     runtime         time        dt             sec / dt    sum         time_NST         time_pressure    time_mask        time_solid"
    close (10)
end subroutine InitBeamFiles

! =================================================================================================================================================================================

subroutine SaveCylinderData(time,dt1,cylinder,forces,volume)
  use share_vars
  implicit none
  real (kind=pr), intent (in) :: time,dt1,volume
  real (kind=pr), dimension (1:4), intent (in) :: Forces
  real (kind=pr), dimension (1:6), intent (in) :: cylinder
  real (kind=pr) :: dx,dy
  integer :: n

  
  dx=xl/real(nx)
  dy=yl/real(ny)

  open (14, file = trim(dir_name)//'/'//trim(simulation_name)//'cylinder', status = 'unknown', access = 'append') ! Append output data file
  ! the funny (1.-mue) term comes from the buoyancy
  write (14, '(15(es15.8,1x))') time, dt1, &
  (cylinder(n), n=1,6), (Forces(n), n=1,4),&
  volume,&
  grav*(y0-cylinder(2))*(1.0-mue)*pi,&
  0.5*mue*pi*(cylinder(3)**2 + cylinder(4)**2)
  close (14)
  
end subroutine


! =================================================================================================================================================================================


subroutine SaveBeamData( time, beam, bpressure, dt1, q, ROC, Forces, force_pressure, FluidIntegralQuantities,it)
  use share_vars
  implicit none
  real (kind=pr), intent (in) :: time,dt1,ROC
  integer, intent (in) :: q,it
  real (kind=pr), dimension (0:ns-1, 1:6), intent (in) :: beam ! Beam indices: 1=beam_x 2=beam_y 3=beam_vx 4=beam_vy 5=theta 6=theta_dot
  real (kind=pr), dimension (0:ns-1),      intent (in) :: bpressure
  real (kind=pr), dimension (1:4), intent (in) :: Forces
  real (kind=pr), dimension (1:2), intent (in) :: force_pressure
  ! 1= vor_rms 2=vor_rms_dot 3=Fluid kinetic Enegry 4=Dissipation 5=Dissipation(Penalty) 6=EnergyInflow
  real (kind=pr), dimension (1:7), intent (in) :: FluidIntegralQuantities 
  real (kind=pr), dimension (0:ns-1) :: p_s, theta_s
  real (kind=pr), dimension (0:ns-1, 0:ns-1) :: D
  character(len=16) :: format_ns1, format_ns, format_ns1_HP
  character(len=3)  :: ns_string, ns1_string
  real (kind=pr) :: alpha, alpha_t, alpha_tt , dx, dy
  real (kind=pr), dimension(1:6) :: LeadingEdge !LeadingEdge: x, y, vx, vy, ax, ay (Array)  
  integer :: n,i, step
  
  dx=xl/real(nx)
  dy=yl/real(ny)  
  
  
  call mouvement(time, alpha, alpha_t, alpha_tt, LeadingEdge )

  step=1 ! save the entire beam? or only every 2..3 points?
  ! set up formats
  write(ns_string, '(I3)') ns
  write(ns1_string, '(I3)') ns+1
  format_ns  = '('//ns_string//'(es12.5,1x))'
  format_ns1 = '('//ns1_string//'(es12.5,1x))'
  format_ns1_HP = '('//ns1_string//'(es15.8,1x))'

  if (iSaveBeam==1) then !save the entire beam
      ! Save the pressure-jumps at the beampoint, after multiplying with startup conditioner
      open  (14, file = trim(dir_name)//'/'//trim(simulation_name)//'beam_pressure', status = 'unknown', access = 'append') ! Append output data file
      write (14, format_ns1_HP) time, (bpressure(n), n=0,ns-1,step)
      close (14) ! Close the file to protect data
      open (14, file = trim(dir_name)//'/'//trim(simulation_name)//'beam_theta', status = 'unknown', access = 'append') ! Append output data file
      write (14, format_ns1_HP) time, (beam(n,5), n=0,ns-1,step)
      close (14) ! Close the file to protect data
      open (14, file = trim(dir_name)//'/'//trim(simulation_name)//'beam_theta_dot', status = 'unknown', access = 'append') ! Append output data file
      write (14, format_ns1_HP) time, (beam(n,6), n=0,ns-1,step)
      close (14) ! Close the file to protect data
      open (14, file = trim(dir_name)//'/'//trim(simulation_name)//'beam_x', status = 'unknown', access = 'append') ! Append output data file
      write (14, format_ns1) time, (beam(n,1), n=0,ns-1,step)
      close (14) ! Close the file to protect data
      open (14, file = trim(dir_name)//'/'//trim(simulation_name)//'beam_y', status = 'unknown', access = 'append') ! Append output data file
      write (14, format_ns1) time, (beam(n,2), n=0,ns-1,step)
      close (14) ! Close the file to protect data
      open (14, file = trim(dir_name)//'/'//trim(simulation_name)//'beam_vx', status = 'unknown', access = 'append') ! Append output data file
      write (14, format_ns1) time, (beam(n,3), n=0,ns-1,step)
      close (14) ! Close the file to protect data
      open (14, file = trim(dir_name)//'/'//trim(simulation_name)//'beam_vy', status = 'unknown', access = 'append') ! Append output data file
      write (14, format_ns1) time, (beam(n,4), n=0,ns-1,step)
      close (14) ! Close the file to protect data
      !---------------------------------------------------------
      call Differentiate1D ( bpressure, p_s, ns, ds, 1)
      open  (14, file = trim(dir_name)//'/'//trim(simulation_name)//'beam_pressure_s', status = 'unknown', access = 'append') ! Append output data file
      write (14, format_ns1_HP) time, (p_s(n), n=0,ns-1,step)
      close (14) 
      !-------------------------------------------------------------
      open (14, file = trim(dir_name)//'/'//trim(simulation_name)//'iter', status = 'unknown', access = 'append') ! Append output data file
      write (14, '(3(es15.8,1x),i6)') time, dt1, ROC, q
      close (14) ! Close the file to protect data
      !------------------------------------------------------- drag/lift and so on
      ! for beam elastic energy, we need theta_s over the beam  
      call Differentiate1D ( beam(:,5), theta_s, ns, ds, 1)      
      open (14, file = trim(dir_name)//'/'//trim(simulation_name)//'beam_data', status = 'unknown', access = 'append') ! Append output data file
      write (14, '(12(es15.8,1x),i3,2x,17(es15.8,1x))') &
             time, bpressure(ns-1),bpressure(3*ns/4),bpressure(ns/2),&
             (beam(ns-1,n), n=1,6),&
             dt1, ROC, q,&
             (Forces(n), n=1,4), (force_pressure(n), n=1,2), &
             FluidIntegralQuantities(1), & 				! 1= vor_rms
             FluidIntegralQuantities(2), &				! 2= vor_rms_dot
             FluidIntegralQuantities(3), &				! 3= Fluid kinetic Enegry
             FluidIntegralQuantities(4), &				! 4= Dissipation
             FluidIntegralQuantities(5), &				! 5= Dissipation(Penalty)
             FluidIntegralQuantities(6), &				! 6= Energy Inflow from solid
             FluidIntegralQuantities(7), &				! 7= Energy Inflow through mean flow forcing
             mue*0.5d0*ds*sum( beam(:,3)**2+beam(:,4)**2 ),& 		! E_kin_solid
             eta*0.5d0*ds*sum(theta_s**2)     ,&			! E_elastic_solid
             mue*grav*ds*sum(beam(:,2)-y0),&				! E_potential_solid
             eps*dx*dy*sum(mask)					! mask volume
      close (14)

  
  elseif ((iSaveBeam==2).or.(iSaveBeam==3)) then !save only the last point
  
      ! for beam elastic energy, we need theta_s over the beam  
      call Differentiate1D ( beam(:,5), theta_s, ns, ds, 1)
      
      open (14, file = trim(dir_name)//'/'//trim(simulation_name)//'beam_data', status = 'unknown', access = 'append') ! Append output data file
      write (14, '(12(es15.8,1x),i3,2x,17(es15.8,1x))') &
             time, bpressure(ns-1),bpressure(3*ns/4),bpressure(ns/2),&
             (beam(ns-1,n), n=1,6),&
             dt1, ROC, q,&
             (Forces(n), n=1,4), (force_pressure(n), n=1,2), &
             FluidIntegralQuantities(1), & 				! 1= vor_rms
             FluidIntegralQuantities(2), &				! 2= vor_rms_dot
             FluidIntegralQuantities(3), &				! 3= Fluid kinetic Enegry
             FluidIntegralQuantities(4), &				! 4= Dissipation
             FluidIntegralQuantities(5), &				! 5= Dissipation(Penalty)
             FluidIntegralQuantities(6), &				! 6= Energy Inflow from solid
             FluidIntegralQuantities(7), &				! 7= Energy Inflow through mean flow forcing
             mue*0.5d0*ds*sum( beam(:,3)**2+beam(:,4)**2 ),& 		! E_kin_solid
             0.5d0*ds*eta*sum(theta_s**2)     ,&			! E_elastic_solid
             mue*grav*ds*sum(beam(:,2)-y0),&				! E_potential_solid
             sum(mask)*eps*(yl/real(ny))*(xl/real(nx))			! mask volume
      close (14)

  elseif ((iSaveBeam==0).or.(iBeam==0)) then !======================================== if you don't want to save the beam, then at least save lift/drag and so on

      open (14, file = trim(dir_name)//'/'//trim(simulation_name)//'beam_data', status = 'unknown', access = 'append') ! Append output data file
      write (14, '(15(es15.8,1x))') time, dt1, (Forces(n), n=1,4), (force_pressure(n), n=1,2), (FluidIntegralQuantities(n),n=1,7)
      close (14)
  endif
  
  ! for beam elastic energy, we need theta_s over the beam  
  call Differentiate1D ( beam(:,5), theta_s, ns, ds, 1)
  
    open (14, file = trim(dir_name)//'/'//trim(simulation_name)//'energies', status = 'unknown', access = 'append') ! Append output data file
    write (14, '(9(es15.8,1x))') &
	    time, &
	    FluidIntegralQuantities(3), &				! 2= Fluid kinetic Enegry
	    FluidIntegralQuantities(4), &				! 3= Dissipation
	    FluidIntegralQuantities(5), &				! 4= Dissipation(Penalty)
	    FluidIntegralQuantities(6), &				! 5= Energy Inflow from solid
	    FluidIntegralQuantities(7), &				! 6= Energy Inflow through mean flow forcing
	    mue*0.5d0*ds*sum( beam(:,3)**2+beam(:,4)**2 ),& 		! 7=E_kin_solid
	    eta*0.5d0*ds*sum( theta_s**2 )     ,&			! 8=E_elastic_solid
	    mue*grav*ds*sum(beam(:,2)-y0)				! 9=E_potential_solid
    close (14)
  
  
  
!   if (iFSI>2) then 
    open (14, file = trim(dir_name)//'/'//trim(simulation_name)//'mouvement', status = 'unknown', access = 'append')
    write (14, '(10(es15.8,1x))') time, alpha, alpha_t, alpha_tt, (LeadingEdge(n), n=1,6)
    close (14)
!   end if 

end subroutine SaveBeamData




subroutine SaveDeflectionLine(time,beam)
  ! --------------------------------------------------------------------------------------------------------
  ! just save deflection line (complete beam, full resolution)
  ! --------------------------------------------------------------------------------------------------------
  use share_vars
  implicit none
  real (kind=pr), intent (in) :: time
  real (kind=pr), dimension (0:ns-1, 1:6), intent (in) :: beam
  ! Beam indices: 1=beam_x 2=beam_y 3=beam_vx 4=beam_vy 5=theta 6=theta_dot
  character(len=22) :: format_ns
  integer :: n,i, nsave=100, step, points

    
  write (format_ns, '("(",(i4.4),"(es12.5,1x))")') 2*ns+1   
  open  (14, file = trim(dir_name)//'/'//trim(simulation_name)//'beam.complete', status = 'unknown', access = 'append')
  write (14, format_ns) time, (beam(n,1), n=0,ns-1), (beam(n,2), n=0,ns-1)
  close (14)

end subroutine SaveDeflectionLine




subroutine SaveBeamPositions(time,beam)
  ! --------------------------------------------------------------------------------------------------------
  ! this subroutine saves a complete deflection line, but it uses only 100 points (if there are enough)
  ! higher resolutions are not really required (and there may be buffer overflows in that case)
  !
  ! corrected 05 apr 2012 now working correctly
  !
  ! --------------------------------------------------------------------------------------------------------
  use share_vars
  implicit none
  real (kind=pr), intent (in) :: time
  real (kind=pr), dimension (0:ns-1, 1:6), intent (in) :: beam
  ! Beam indices: 1=beam_x 2=beam_y 3=beam_vx 4=beam_vy 5=theta 6=theta_dot
  character(len=16) :: format_ns1, format_ns, format_ns1_HP
  character(len=3)  :: ns_string, ns1_string
  integer :: n,i, nsave=100, step, points

  if ((iSaveBeam==3).and.(ns>=100)) then
  !more than 100 points on the beam
      step=ns/nsave ! integer division
      points = ns/step !integer division
      
      write(ns_string, '(I3)') 2*points+1+4 !take a bit more place in the file (time+safety)
      format_ns  = '('//ns_string//'(es12.5,1x))'  
      open  (14, file = trim(dir_name)//'/'//trim(simulation_name)//'beam.complete', status = 'unknown', access = 'append')
      write (14, format_ns) time, (beam(n,1), n=0,ns-1,step), (beam(n,2), n=0,ns-1,step)
      close (14)
  elseif ((iSaveBeam==3).and.(ns<100)) then
  !less than 100 points
      open  (14, file = trim(dir_name)//'/'//trim(simulation_name)//'beam.complete', status = 'unknown', access = 'append')
      write (14, '(201(es12.5,1x))') time, (beam(n,1), n=0, ns-1), (beam(n,2), n=0, ns-1)
      close (14)
  endif

end subroutine SaveBeamPositions








subroutine CreateHeader()
  use share_vars
  !----------------------------------------------------------------
  ! This routine reads in the PARAMS file and stores it in a
  ! global array. This array can now be written in every
  ! output file, for disambiguation.
  !----------------------------------------------------------------
  implicit none
  integer, parameter :: Nchar = 256 ! remember to change in share_vars, too, if required.
  integer :: io_error=0, linenumber=0,i=0
  character (len=Nchar) :: dummy
  !-------------------------------------------------------------------------------------------------------------
  ! first, lets determine how many entrys the PARAMS file has
      open (14, file = 'PARAMS.m', status='old', action='read')
      do while (io_error==0)
        read (14,'(A)',iostat=io_error) dummy
        linenumber=linenumber+1
      enddo
      linenumber=linenumber-1 !counted one too far
      close (14)
  !-------------------------------------------------------------------------------------------------------------
  ! then allocate the table for the header and fill the array
      open (14, file = 'PARAMS.m', status='old', action='read')
      allocate (Params_Header(1:linenumber+2))
      Params_Header = ""
      Params_Header(1)            ="% ======================================PARAMS begin======================================================="
      Params_Header(linenumber+2) ="% ======================================PARAMS END========================================================="
      do i=2,linenumber+1
        read (14,'(A)',iostat=io_error) dummy
        Params_Header(i) = "% "//dummy
      enddo
      close (14)
  !-------------------------------------------------------------------------------------------------------------
end subroutine CreateHeader




!================================================================================================================================================




subroutine WriteHeaderToFile(file_identifier)
  !----------------------------------------------------------------
  ! This routine writes the PARAMS file as a backup in the specified file_identifier
  !----------------------------------------------------------------
  use share_vars
  implicit none
  integer, intent (in) :: file_identifier
  integer :: i
  !------------------------------------------------------------------
  do i=1,size(Params_Header)
    write(file_identifier,'(A)') Params_Header(i)
  enddo
end subroutine WriteHeaderToFile
subroutine save_fields (n1, time, dt1, vortk, nlk, workvis, nbackup, ivideo, u, press, beams)
  use share_vars
  use FieldExport
  use SpectralSolver
  implicit none
  integer :: ix, iy
  integer, intent (inout) :: nbackup
  integer, intent (in) :: n1,ivideo
  real (kind=pr), intent (in) :: time, dt1
  real (kind=pr), dimension (0:nx-1, 0:ny-1, 0:1), intent (in)  :: nlk, vortk, workvis
  real (kind=pr), dimension (0:nx-1, 0:ny-1, 1:2), intent (in) :: u
  real (kind=pr), dimension (0:nx-1,0:ny-1) ::  vort, work1, work2, work3
  real (kind=pr), dimension (0:nx-1, 0:ny-1), intent (in) :: press
  type (solid), dimension(1:iBeam), intent(in) :: beams
  character (len=11) :: name
  character (len=1) :: name1
  real (kind=pr) :: Mean_ux, Mean_uy, delx, dely
  delx = xl / real (nx)
  dely = yl / real (ny)
!--Set up file name base
  write (name, '(es10.4)') time


  !=================================================================================
  !--Save pressure
  if (iSavePress > 0) then
  call SaveField( trim(dir_name)//'/fields/'//trim(simulation_name)//'p_'//name, press, iSavePress, xl,yl, "pressure")
  endif
  !=================================================================================
  !--Save velocity and streamfunction
  if (iSaveVel>0) then
  call SaveField( trim(dir_name)//'/fields/'//trim(simulation_name)//'ux_'//name, u(:,:,1), iSaveVel, xl,yl, "x-velocity")
  call SaveField( trim(dir_name)//'/fields/'//trim(simulation_name)//'uy_'//name, u(:,:,2), iSaveVel, xl,yl, "y-velocity")
  endif
  
  if (iSaveSTR>0) then
  call poisson (vortk, vort) ! vort = streamfunction for the moment
  call SaveField( trim(dir_name)//'/fields/'//trim(simulation_name)//'str_'//name, vort, iSaveSTR, xl,yl, "streamfunction")  
  endif
  !=================================================================================
  !--Save vorticity
  if (iSaveVort > 0) then
  call cofitxy (vortk, vort)
  call SaveField( trim(dir_name)//'/fields/'//trim(simulation_name)//'vor_'//name, vort, iSaveVort, xl,yl,"vorticity")
  call SaveField( trim(dir_name)//'/fields/'//trim(simulation_name)//'vor_mean_'//name, vor_mean, iSaveVort, xl,yl,"vorticity")
  endif
  !=================================================================================
  !save mask
  if (iSaveMask > 0) call SaveField( trim(dir_name)//'/fields/'//trim(simulation_name)//'mask_'//name, mask, iSaveMask, xl,yl,"mask")

  !=================================================================================
  if (iSaveMaskVel >0) then
  call SaveField( trim(dir_name)//'/fields/'//trim(simulation_name)//'mask_vy_'//name, maskvy, iSaveMaskVel, xl,yl, "mask_vx")
  call SaveField( trim(dir_name)//'/fields/'//trim(simulation_name)//'mask_vx_'//name, maskvx, iSaveMaskVel, xl,yl, "mask_vy")
  endif

!=================================================================================
   if (iSaveStress>0) then
       call coftxy  ( u(:,:,1), vort  )
       call cofdx   ( vort    , work1 ) 
       call cofitxy ( work1   , work3 ) 
       
       work3 = 2.0 * nu * work3
       call SaveField( './fields/'//trim(simulation_name)//'stress_a_'//name, work3, iSaveStress, xl,yl,"viscous-stress-a")
       !------------------------
       call coftxy  ( u(:,:,1), vort  )
       call cofdy   ( vort    , work1 ) 
       call cofitxy ( work1   , work2 ) 
       
       call coftxy  ( u(:,:,2), vort  )
       call cofdx   ( vort    , work1 ) 
       call cofitxy ( work1   , work3 ) 
       
       work3 = nu * ( work2 + work3 )
       call SaveField( './fields/'//trim(simulation_name)//'stress_b_'//name, work3, iSaveStress, xl,yl,"viscous-stress-b")
       !----------finite differences
       do ix=1,nx-2
       do iy=1,ny-2
	work3(ix,iy)=(u(ix+1,iy,1)-u(ix-1,iy,1))/2.0/delx
       enddo
       enddo
       work3 = 2.0 * nu * work3
       call SaveField( './fields/'//trim(simulation_name)//'stress_a_FD_'//name, work3, iSaveStress, xl,yl,"viscous-stress-a")
       
       do ix=1,nx-2
       do iy=1,ny-2
	work2(ix,iy)=(u(ix,iy+1,1)-u(ix,iy-1,1))/2.0/dely
	work3(ix,iy)=(u(ix+1,iy,2)-u(ix-1,iy,2))/2.0/delx
       enddo
       enddo
       work3 = nu * ( work2 + work3 )
       call SaveField( './fields/'//trim(simulation_name)//'stress_b_FD_'//name, work3, iSaveStress, xl,yl,"viscous-stress-b")

   endif


!=================================================================================
!--Backup data
  if (nbackup == 2) then
     nbackup = 0 ! do not backup this time
  else
    call MakeRuntimeBackup(n1, time, dt1, vortk, nlk, workvis, nbackup, ivideo, u, beams)
  endif

end subroutine save_fields

!=====================================================================================================================================

subroutine MakeRuntimeBackup(n1, time, dt1, vortk, nlk, workvis, nbackup, ivideo, u, beams)
  use share_vars
  implicit none
  integer, 					   intent (inout) :: nbackup
  integer, 					   intent (in) :: n1,ivideo
  real (kind=pr),                                  intent (in) :: time, dt1
  real (kind=pr), dimension (0:nx-1, 0:ny-1, 0:1), intent (in) :: nlk, vortk, workvis, u
  type (solid), dimension(1:iBeam),                intent (in) :: beams
  integer :: i
  character (len=1) :: name1

  write (*,*) '*** Making a backup, time=',time
  
  write (name1, '(I1)') nbackup
  open (15, file = trim(dir_name)//'/runtime_backup'//name1//'.in', form='unformatted', status='replace')
  write (15) time
  write (15) 123.0 ! checksum
  write (15) n1, dt1, vortk, nlk, workvis, mask, maskvx, maskvy, ivideo, u
  
  write (15) ns  
  write (15) colorscale, colorscale_done
  
  ! dump all beams to disk
  do i = 1, iBeam
    write (15) beams(i)%x, beams(i)%y, beams(i)%vx, beams(i)%vy, beams(i)%theta, beams(i)%theta_dot
    write (15) beams(i)%pressure_new, beams(i)%pressure_old, beams(i)%tau_new, beams(i)%tau_old
    write (15) beams(i)%Force, beams(i)%Force_unst, beams(i)%Force_press, beams(i)%E_kinetic, beams(i)%E_pot, beams(i)%E_elastic , beams(i)%x0, beams(i)%y0
    write (15) beams(i)%AngleBeam, beams(i)%iMouvement, beams(i)%drag_unst_new, beams(i)%drag_unst_old, beams(i)%lift_unst_new, beams(i)%lift_unst_old
    write (15) beams(i)%UnsteadyCorrectionsReady, beams(i)%dt_old, beams(i)%beam_oldold
    write (15) beams(i)%ax, beams(i)%ay, beams(i)%Inertial_Force
  enddo
  
  write (15) vor_mean
  
  close (15)
  nbackup = 1 - nbackup
  
  write (*,*) '*** excecuting: tar czf backup_'//name1//'.tar.gz '//trim(simulation_name)//'*'
  call system('tar czf backup_'//name1//'.tar.gz '//trim(simulation_name)//'*')
   
  
  write (*,*) '*** information: backup done'

end subroutine MakeRuntimeBackup


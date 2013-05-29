subroutine init_fields (n1, time, dt1, vortk, nlk, workvis, beam, bpressure, tau_beam_old, ivideo, u)
  use share_vars
  use FieldExport
  use SpectralSolver
  implicit none
!=======================================================
! This subroutine set up the inital flow which is normally
! a mean flow.
!======================================================
  integer :: ibackup, ix, iy, nx_file, ny_file
  integer :: ierr = 0
  integer, intent (out) :: ivideo
  real (kind=pr) :: time1
  real (kind=pr) :: x, y
  real (kind=pr) :: x1v, x2v, x3v,r0,we,d,r1,r2
  real (kind=pr) :: y1v, y2v, y3v
  real (kind=pr) :: a1, a2, a3
  real (kind=pr) :: s1, s2, s3,dx,dy
  real (kind=pr) :: w1, w2, w3,  y_chan, H_effective, Mean_Ux
  ! Beam indices: 1=beam_x 2=beam_y 3=beam_vx 4=beam_vy 5=theta 6=theta_dot
  real (kind=pr), dimension (0:ns-1, 1:6), intent (inout) :: beam
  real (kind=pr), dimension (0:nx-1, 0:ny-1) :: vort, work1
  real (kind=pr), dimension (0:ns-1), intent (inout) :: bpressure, tau_beam_old
  integer, intent (inout) :: n1
  real (kind=pr), intent (inout) :: time, dt1
  real (kind=pr), dimension (0:nx-1, 0:ny-1, 0:1), intent (out)  :: nlk, vortk, workvis, u
  character (len=1) :: name1

! Assign zero values
  time 		= 0.
  dt1 		= 0.
  vortk 	= 0.
  nlk 		= 0.
  dx 		= xl/real(nx)
  dy 		= yl/real(ny)
  bpressure 	= 0.0
  tau_beam_old 	= 0.0
  ivideo 	= 1 ! counter for the snapshots
  
  write (*,'(" --- inicond = ",i1)') inicond

  if (inicond == 1) then !==================================================================================================
      ! normal startup
      vort 		= 0.0
      vortk 		= 0.0
      workvis 		= 0.0
      u 		= 0.0
  elseif (inicond == 2) then !==================================================================================================
      ! Read from backup file
      do ibackup = 0, 1
      write (name1, '(I1)') ibackup
      open (15, file = trim(dir_name)//'/runtime_backup'//name1//'.in', form='unformatted', status='old', iostat=ierr, err=100)
        read (15) time1
        if (time1 > time) then
          time = time1
          write (*,*)'last record in runtime_backup'//name1//'.in', ' is at time = ', time
          read (15) n1, dt1, vortk, nlk, workvis, mask, maskvx, maskvy, beam, bpressure, ivideo, u, tau_beam_old
        endif
      close (15)
      100 enddo
      
      if (time1 == 0) then
      write (*,*) 'Unable to resume'
      open (14, file = 'Unable to resume- change inicond', status = 'replace')
      close (14)
      stop
      endif
  elseif (inicond == 3) then !==================================================================================================
      ! channel initial condition, linear vorticity profile
      vort = 0.0
      !set the vorticity profile
      H_effective=yl-2.0*h_channel
      
      do iy=0,ny-1
      ! the Poisseuille flow requires a linear vorticity
      y_chan = real(iy)*dy-h_channel
      if ((y_chan>0.0).and.(y_chan<H_effective)) then
      vort(:,iy) = (-6.0*Mean_ux(time)/H_effective)*(1.0-2.0*y_chan/H_effective)
      else
      vort(:,iy) = 0.0
      endif
      enddo
      vort = vort * (1.0-mask*eps) ! no initial vorticity inside solid region
      call coftxy(vort,vortk)
  elseif (inicond == 99) then !==================================================================================================
      call coftxy( vor_init, vortk(:,:,1) ) ! in this scenario, mask_sponge contains the last vort field in a previous simulation, so just transform it
      call coftxy( vor_init, vortk(:,:,0) ) ! in this scenario, mask_sponge contains the last vort field in a previous simulation, so just transform it
      workvis=0.0
      call SaveGif (vor_init, trim(dir_name)//"/initial_condition")
      u=0.0
      time=time_init
  elseif (inicond == 111) then
      ! ------------------------------------
      ! reads in a binary vorticity file.
      ! ------------------------------------
      if ( len_trim(inicond_file)==0 ) then
        write (*,*) 'inicond_file not set.. suicide!'
        stop
      endif
      
      open (10, file=trim(inicond_file), form='unformatted', status='old')
      read (10) nx_file, ny_file, vort
      close (10)
      
      if ((nx_file==nx).and.(ny_file==ny)) then
	  call coftxy(vort,vortk)
	  call cal_velocity(0.d0, vortk, u, vort )
      else
	  write (*,*) 'inicond field resolution does not match present resolution. suicide!'
	  stop
      endif
      
  elseif (inicond==4) then
      vort = 0.d0
      do ix=0,nx-1
      do iy=0,ny-1
        vort(ix,iy)= 100.d0*exp(- ( (real(ix)*dx-0.5*xl)**2 + (real(iy)*dy-0.75*yl)**2 )/ (10.0d0*dx)**2 )
      enddo
      enddo
      vort=vort*(1.d0-eps*mask)
      call SaveGif(vort,"vort_init")
      call coftxy(vort,vortk)
      workvis = 0.0
      call cal_velocity(0.d0, vortk, u, vort )
      ivideo=1 ! counter for the snapshots
  elseif (inicond==7) then
      ! ----------------------------------------------------
      ! vortex "wake" arangement 06.05.2012
      ! ----------------------------------------------------
      vort = 0.d0
      a1 = 5.d0 ! single vortex width
      w1 =35.d0 ! vortex strength
      do ix=0,nx-1
      do iy=0,ny-1
        vort(ix,iy) =+w1*1.0001d0*exp(- ( (real(ix)*dx-1.0)**2 + (real(iy)*dy-2.5)**2 )/ (a1*dx)**2 ) &
                     -w1*1.d0*exp(- ( (real(ix)*dx-2.5)**2 + (real(iy)*dy-1.5)**2 )/ (a1*dx)**2 ) &
                     +w1*1.d0*exp(- ( (real(ix)*dx-4.0)**2 + (real(iy)*dy-2.5)**2 )/ (a1*dx)**2 ) &
                     -w1*1.d0*exp(- ( (real(ix)*dx-5.5)**2 + (real(iy)*dy-1.5)**2 )/ (a1*dx)**2 ) &
                     +w1*1.d0*exp(- ( (real(ix)*dx-7.0)**2 + (real(iy)*dy-2.5)**2 )/ (a1*dx)**2 ) &
                     -w1*1.d0*exp(- ( (real(ix)*dx-8.5)**2 + (real(iy)*dy-1.5)**2 )/ (a1*dx)**2 ) &
                     +w1*1.d0*exp(- ( (real(ix)*dx-10.0)**2 + (real(iy)*dy-2.5)**2 )/ (a1*dx)**2 ) &
                     -w1*1.d0*exp(- ( (real(ix)*dx-11.4)**2 + (real(iy)*dy-1.5)**2 )/ (a1*dx)**2 )             
      enddo
      enddo

      call SaveGif(vort,"vort_init")
      call coftxy(vort,vortk)
      workvis = 0.0
      call cal_velocity(0.d0, vortk, u, vort )
      ivideo=1 ! counter for the snapshots    
  
  elseif (inicond==9) then
      ! ----------------------------------------------------
      ! vortex "street" arangement 06.05.2012
      ! ----------------------------------------------------
      vort = 0.d0
      a1 = 5.d0 ! single vortex width
      w1 =35.d0 ! vortex strength
      do ix=0,nx-1
      do iy=0,ny-1
        vort(ix,iy) =+w1*1.0001d0*exp(- ( (real(ix)*dx-1.0)**2 + (real(iy)*dy-2.0)**2 )/ (a1*dx)**2 ) &
                     +w1*1.d0*exp(- ( (real(ix)*dx-2.5)**2 + (real(iy)*dy-2.0)**2 )/ (a1*dx)**2 ) &
                     +w1*1.d0*exp(- ( (real(ix)*dx-4.0)**2 + (real(iy)*dy-2.0)**2 )/ (a1*dx)**2 ) &
                     +w1*1.d0*exp(- ( (real(ix)*dx-5.5)**2 + (real(iy)*dy-2.0)**2 )/ (a1*dx)**2 ) &
                     +w1*1.d0*exp(- ( (real(ix)*dx-7.0)**2 + (real(iy)*dy-2.0)**2 )/ (a1*dx)**2 ) &
                     +w1*1.d0*exp(- ( (real(ix)*dx-8.5)**2 + (real(iy)*dy-2.0)**2 )/ (a1*dx)**2 ) &
                     +w1*1.d0*exp(- ( (real(ix)*dx-10.0)**2 + (real(iy)*dy-2.0)**2 )/ (a1*dx)**2 ) &
                     +w1*1.d0*exp(- ( (real(ix)*dx-11.4)**2 + (real(iy)*dy-2.0)**2 )/ (a1*dx)**2 )                     
      enddo
      enddo

      call SaveGif(vort,"vort_init")
      call coftxy(vort,vortk)
      workvis = 0.0
      call cal_velocity(0.d0, vortk, u, vort )
      ivideo=1 ! counter for the snapshots     
elseif (inicond==8) then
      !------------------------------------------------
      ! vortex pair 06.05.2012
      !------------------------------------------------
      vort = 0.d0
      a1 = 5.d0 ! single vortex width
      w1 =35.d0 ! vortex strength
      do ix=0,nx-1
      do iy=0,ny-1
        vort(ix,iy) =+w1*exp(- ( (real(ix)*dx-1.0)**2 + (real(iy)*dy-2.5)**2 )/ (a1*dx)**2 ) &
                     +w1*exp(- ( (real(ix)*dx-2.5)**2 + (real(iy)*dy-1.5)**2 )/ (a1*dx)**2 )       
      enddo
      enddo

      call SaveGif(vort,"vort_init")
      call coftxy(vort,vortk)
      workvis = 0.0
      call cal_velocity(0.d0, vortk, u, vort )
      call SaveGif(u(:,:,1),"ux_init")
      ivideo=1 ! counter for the snapshots   
elseif (inicond==22) then
      !------------------------------------------------
      ! random inital condition 27.09.2012
      !------------------------------------------------
      do ix=0,nx-1
      do iy=0,ny-1
        call random_number(w1)
        vort(ix,iy) = 10.0*w1
      enddo
      enddo
      call SaveGif(vort,"vort_init")
      call coftxy(vort,vortk)
      
      workvis = 0.0
      call cal_velocity(0.d0, vortk, u, vort )
      call SaveGif(u(:,:,1),"ux_init")
      ivideo=1 ! counter for the snapshots         
elseif (inicond == 55 ) then
    x0 = 0.5*xl
    y0 = 0.5*yl
    r0 = 0.1
    we = 299.528385375226
    d = 0.1
    
    do ix=0,nx-1
    do iy=0,ny-1
       r1 = sqrt( (real(ix)*dx-x0)**2 + (real(iy)*dy-y0-d)**2 ) / r0
       r2 = sqrt( (real(ix)*dx-x0)**2 + (real(iy)*dy-y0+d)**2 ) / r0
      vort(ix,iy) = we * (1.0-r1**2)*exp(-r1**2) - we * (1.0-r2**2)*exp(-r2**2)
    enddo
    enddo
    
    call SaveGIF(vort,'inicond.vor')
    call coftxy (vort, vortk)
    call cal_velocity(0.0, vortk, u, work1 )
    call SaveGIF(u(:,:,0),'inicond.ux')
    call SaveGIF(u(:,:,1),'inicond.uy')
    ivideo = 1
  else
      write (*,'(A)') '!!! Error: Invalid initial condition (inicond)'
  endif

end subroutine init_fields













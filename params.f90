subroutine params ()
  use share_vars
  implicit none
  real (kind=pr)                       :: checksum
  real (kind=pr)                       :: A, dummy

  real (kind=pr) :: dx, dy

  read *, dummy !----------------
  read *, nx
  read *, ny
  read *, ns
  write (*,'(" --- Resolution: nx=",i5," ny=",i5," ns=",i4)') nx,ny,ns
  read *, Time_end  
  read *, dt_fixed
  read *, cfl
  read *, dummy !----------------
  read *, eta  
  read *, mue
  read *, T_release 	!time to begin the release-process for the beam  
  read *, tau 	! duration for the soft-startup. during TAU, couppling is set from 0 to 1 by a second degree polynome.
  read *, iFSI
  read *, AngleBeam
  read *, grav
  read *, iImpulse
  read *, sigma
  read *, TimeMethodSolid
  read *, dummy !----------------
  read *, iCylinder
  read *, iWalls
  read *, iBeam
  read *, dummy !----------------
  read *, tsave
  read *, tdrag
  read *, iSavePress
  read *, iSaveVel
  read *, iSaveVort
  read *, iSaveSTR
  read *, iSaveStress
  read *, iSaveMask
  read *, iSaveMaskVel
  read *, iSaveBeam
  read *, dummy !----------------
  read *, nu
  nu=1.d0/nu
  if (nu>1.0) then
   write (*,'(A)') "Oh Oh. nu>1.0. Probably you work with an old PARAMS.m file, in this version, please set the REYNOLDS number here"
   write (*,'(A)') "so instead of nu=1e-2 set Re=100. its easier like this."
   stop
  endif 
  read *, eps
  read *, iMeanVelocity
  read *, T_fullspeed
  read *, theta_inf
  read *, idealis
  read *, iSponge
  read *, iSpongeType
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  eps_sponge=eps
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  read *, dummy !----------------
  read *, inicond
  read *, iIteration
  read *, iMultiRes
  read *, dummy !----------------
  read *, xl
  read *, yl
  read *, R_cylinder
  read *, x0
  read *, y0
  read *, length
  read *, t_beam
  read *, h_channel
  read *, N_smooth
  read *, SpongeSize  
  read *, dummy !----------------
  read *, checksum
  read *, simulation_name
  if (inicond == 111) then
  read *, inicond_file
  endif
  
  if (mod(nx,2)>0) then
    call DisplayError
    write(*,*) "!!! crucial: use even number for grid points"
  stop
  endif
  if (mod(ny,2)>0) then
    call DisplayError
    write(*,*) "!!! crucial: use even number for grid points"
  stop
  endif

  if ( (ns<32).and.((iFSI==1).or.(iFSI==10)) ) then
    call DisplayError
    write(*,*) "??? you set ns<32 and you want to do FSI. this is probably not enough, please correct this"
  stop
  endif 
  
  if ((iSaveBeam==0).and.(iFSI==1)) then
    call DisplayError
    write (*,*) " ??? You want to compute FSI without storing the beam data? Really?"
    write (*,*) " I won't cry and abandon, but that's strange you admit."
  endif
  
  if (iIteration==0) write(*,*) "*** Sequential Staggered Algorithm"
  if (iIteration==1) write(*,*) "*** Iteration Algorithm"
  if ((iIteration.ne.0).and.(iIteration.ne.1)) then 
    call DisplayError
    write(*,*) "!!! Error: Wrong choice of algorithm (iIteration)"
    stop
  endif

  ! fixed values; unused
  ihypvisc = 0
  iobst = 1
  
  if ((iBeam==0).and.(iWalls==0).and.(iCylinder==0)) then
    call DisplayError
    write (*,*) "!!! PARAMS: you didn't set any flag for the mask. "
  endif

  if ( (iMeanVelocity==4).and.(iFSI==0).and.(inicond==1) ) then
    call DisplayError
    write (*,*) "!!! PARAMS: fluid at rest and a fixed obstacle? That is completely senseless... "
    stop
  endif

  if     (TimeMethodSolid == EulerImplicit) then
  write(*,'(A)') " $$$ The solid solver is using EulerImplicit "
  elseif (TimeMethodSolid == CrankNicholson) then
  write(*,'(A)') " $$$ The solid solver is using CrankNicholson "
  elseif (TimeMethodSolid == BDF2) then
  write(*,'(A)') " $$$ The solid solver is using BDF2 "
  elseif (TimeMethodSolid == RungeKutta4) then
  write(*,'(A)') " $$$ The solid solver is using RK4 "
  else
  write(*,'(A)') " !!! Solid Solver: invalid choice of time stepping algorithm "
  stop
  endif

! Set up various parameters
  pi     = 4.0 * atan (1.0)
  scalex = (2.0*pi / xl)**2
  scaley = (2.0*pi / yl)**2
  

  if (iBeam>0) then
    if ((x0>xl).or.(y0>yl).or.(x0+length>xl)) then
      write (*,*) "!!! stupid error: beam out of domain"
      stop
    endif
  endif
  
  if (iMultiRes==5) then !galilean run
    if (T_release>0.0) then 
      write (*,*) "??? you want to do a galilean comparison, but T_release is not 0.0"
      stop
    endif
    if (tau>0.0) then 
      write (*,*) "??? you want to do a galilean comparison, but tau is not 0.0"
      stop
    endif
  endif

  
  if (checksum .ne. 555.) then
   call DisplayError
   write(*,*) "!!! bad params-file! the checksum number should be 555"
   stop
  endif
  theta_inf = theta_inf * pi/180.0 
    
  if (eps>1.0E-1) then
    write (*,*) "??? eps too large", eps
!    stop
  endif 
  
  if (((iFSI==1).or.(iFSI>9)).and.(iSaveBeam==0)) then
  write (*,*) "??? you simulate FSI without saving the beam data. That seems likely to be not what you intended"
  iSaveBeam = 3
  endif

  if (eps<1.0E-10) then
    write (*,*) "??? eps too small", eps
!     stop
  endif 
  
  if ( (iSponge>0).and.(iSponge.ne.4).and.(eps_sponge<1.0e-9) ) then
    write (*,'(A)') "???  You specified the usage of a vorticity sponge without giving me an epsilon for it. I leave here and let you correct your mistake."
    stop
  endif

  ds = length/(ns-1)
  dx = xl / real (nx)
  dy = yl / real (ny)

  
  simulation_name=trim(simulation_name)
  simulation_name_org=simulation_name !copy the name (for multi-res batch runs!)

  if ((iFSI.ne.1).and.(iFSI.ne.10).and.(iIteration==1)) then
    call DisplayError
    write (*,*) "!!! I'm really sorry but the combination of iFSI and iIteration does not make sense. I'm sure you err and I abort now."
    stop
  endif
  

  if ((iWalls==1).and.(theta_inf.ne.0.0)) then
    call DisplayError
    write (*,*) "??? Really? You impose a channel but also an angle for the mean flow? I'm quite sure you didn't intent that. I abort here."
    stop
  endif
  
  if (abs(N_smooth)<1e-4) then
    write (*,*) "*** information: very small smoothing, so I switch to sharp masks!"
    sharp=.true.
  endif
  
  if ( (iSaveBeam.ne.0).and.( (iFSI==0).or.(iFSI==2).or.(iFSI==3).or.(iFSI==8) ) ) then
      open  (10, file = 'WARNINGS', status = 'replace')
      write (10,*)   " you save the beam; therefore at every time step will compute the pressure, although this is not nessesairy"
      write (10,*)   " you can save time by setting iSaveBeam=0"
      close (10)
      
  endif
  
end subroutine params




subroutine ReLoadParams
  use share_vars
  implicit none
  real (kind=pr)                       :: checksum
  real (kind=pr)                       :: A, dummy
  real (kind=pr) :: dx, dy
  
  open  (10, file = 'PARAMS.m', status = 'old')
  write (*,*)   "$$$  ---------- Re-Load PARAMS ----------"
  read (10,*) dummy !----------------
  read (10,*) dummy !nx
  read (10,*) dummy !ny
  read (10,*) dummy !ns
  read (10,*) Time_end  
  read (10,*) dt_fixed
  read (10,*) cfl
  read (10,*) dummy !----------------
  read (10,*) eta  
  read (10,*) mue
  read (10,*) T_release 	!time to begin the release-process for the beam  
  read (10,*) tau 	! duration for the soft-startup. during TAU, couppling is set from 0 to 1 by a second degree polynome.
  read (10,*) dummy !iFSI
  read (10,*) dummy !AngleBeam
  read (10,*) grav
  read (10,*) iImpulse
  read (10,*) sigma
  read (10,*) dummy
  read (10,*) dummy !----------------
  read (10,*) dummy !iCylinder
  read (10,*) dummy !iWalls
  read (10,*) dummy !iBeam
  read (10,*) dummy !----------------
  read (10,*) tsave
  read (10,*) tdrag
  read (10,*) iSavePress
  read (10,*) iSaveVel
  read (10,*) iSaveVort
  read (10,*) iSaveSTR
  read (10,*) iSaveStress
  read (10,*) iSaveMask
  read (10,*) iSaveMaskVel
  read (10,*) iSaveBeam
  read (10,*) dummy !----------------
  read (10,*) nu
  nu=1.d0/nu
  if (nu>1.0) then
   write (*,'(A)') "Oh Oh. nu>1.0. Probably you work with an old PARAMS.m file, in this version, please set the REYNOLDS number here"
   write (*,'(A)') "so instead of nu=1e-2 set nu=100. its easier like this."
   stop
  endif 
  
  read (10,*) eps
  read (10,*) dummy !iMeanVelocity
  read (10,*) T_fullspeed
  read (10,*) dummy !theta_inf
  read (10,*) idealis
  read (10,*) dummy !iSponge
  read (10,*) dummy !iSpongeType
  eps_sponge=eps
  read (10,*) dummy !----------------
  read (10,*) dummy !inicond
  read (10,*) dummy !iIteration
  read (10,*) dummy !iMultiRes
  read (10,*) dummy !----------------
  read (10,*) dummy !xl
  read (10,*) dummy !yl
  read (10,*) dummy !R_cylinder
  read (10,*) dummy !x0
  read (10,*) dummy !y0
  read (10,*) dummy !length
  read (10,*) dummy !t_beam
  read (10,*) dummy !h_channel
  read (10,*) N_smooth
  read (10,*) dummy !SpongeSize  
  read (10,*) dummy !----------------
  read (10,*) checksum

  close(10)
 
end subroutine ReLoadParams







subroutine GoodBye()
  write(*,'(A)')"  ____                 _ _                "
  write(*,'(A)')" / ___| ___   ___   __| | |__  _   _  ___ "
  write(*,'(A)')"| |  _ / _ \ / _ \ / _` | '_ \| | | |/ _ \"
  write(*,'(A)')"| |_| | (_) | (_) | (_| | |_) | |_| |  __/"
  write(*,'(A)')" \____|\___/ \___/ \__,_|_.__/ \__, |\___|"
  write(*,'(A)')"                               |___/      "
end subroutine 

subroutine DisplayError()
    write (*,'(A)') " _____ ____  ____   ___  ____  _ "
    write (*,'(A)') "| ____|  _ \|  _ \ / _ \|  _ \| |"
    write (*,'(A)') "|  _| | |_) | |_) | | | | |_) | |"
    write (*,'(A)') "| |___|  _ <|  _ <| |_| |  _ <|_|"
    write (*,'(A)') "|_____|_| \_\_| \_\\___/|_| \_(_)"
end subroutine

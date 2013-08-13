! returns the time and position dependend mean velocity profile
! ------------------------------------------------------------

real (kind=pr) function Mean_Ux(time)
  use share_vars
  implicit none
  real (kind=pr), intent(in) :: time
  real (kind=pr) :: timefactor, tau1
  
  if (iMeanVelocity==1) then
        tau1 = T_fullspeed !full speed reached at this time
        !-----------------------------------------------------------------------------------------------------
        !in order to overcome the startup-difficulty couppling is slowly activatd
        if (time<tau1) then
            timefactor =  (time**3)/(-0.5*tau1**3)   + 3*(time**2)/tau1**2
        else
            timefactor = 1.0
        endif
        Mean_Ux = cos (theta_inf) * timefactor
  elseif (iMeanVelocity==2) then
        tau1 = T_fullspeed !full speed reached at this time
        if (time<tau1) then
          timefactor = 0.5*(1.0-cos(0.5*pi*time/tau1))
        else 
          timefactor = 1.0
        endif
        Mean_Ux=cos (theta_inf) *timefactor
  elseif (iMeanVelocity==3) then
        Mean_Ux = cos (theta_inf)
  elseif (iMeanVelocity==4) then !fluid does not move at all
        Mean_Ux = 0.0
  elseif (iMeanVelocity==6) then
	Mean_ux = time 
  elseif (iMeanVelocity==11) then ! zero mean pressure forcing
	Mean_ux = u_mean(1)
  elseif (iMeanVelocity==22) then 
        Mean_ux = .1
  endif


  return  
end function Mean_Ux


!-------------------------------------------------------
real (kind=pr) function Mean_Uy(time)
  use share_vars
  implicit none
  real (kind=pr), intent(in) :: time
  real (kind=pr) :: timefactor, tau1
  
  if (iMeanVelocity==1) then
        tau1 = T_fullspeed !full speed reached at this time
        !-----------------------------------------------------------------------------------------------------
        !in order to overcome the startup-difficulty couppling is slowly activatd
        if (time<tau1) then
            timefactor =  (time**3)/(-0.5*tau1**3)   + 3*(time**2)/tau1**2
        else
            timefactor = 1.0
        endif
        Mean_Uy = sin (theta_inf) * timefactor
  elseif (iMeanVelocity==2) then
        tau1 = T_fullspeed !full speed reached at this time
        if (time<tau1) then
          timefactor = 0.5*(1.0-cos(0.5*pi*time/tau1))
        else 
          timefactor = 1.0
        endif
        Mean_Uy = sin (theta_inf) *timefactor
  elseif (iMeanVelocity==3) then
        Mean_Uy = sin (theta_inf)
  elseif (iMeanVelocity==4) then !fluid does not move at all
        Mean_Uy = 0.0!uy_mean
  elseif (iMeanVelocity==6) then
	Mean_uy = 0.0
  elseif (iMeanVelocity==11) then
	Mean_uy = u_mean(2)	
  else
        Mean_uy = 0.d0
  endif


  return
end function Mean_Uy

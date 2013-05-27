subroutine poisson (f, ans)
! Calculate solution to the Poisson equation f = -grad^2 ans in Fourier space
  use share_vars
  implicit none
  real (kind=pr), dimension (0:nx-1, 0:ny-1), intent (in) :: f
  real (kind=pr), dimension (0:nx-1, 0:ny-1), intent (out) :: ans
  integer :: kx, kx_max, ky, ky_max
  real (kind=pr) :: quot

  kx_max = (nx/2-1) 
  ky_max = (ny/2-1)
  !$omp parallel do private(kx,ky,quot) 
  do ky = 0, ky_max
     do kx = 0, kx_max
        if ( (kx == 0) .and. (ky == 0) ) then
           ans (2*kx, 2*ky)     = 0.0
           ans (2*kx+1, 2*ky)   = 0.0
           ans (2*kx, 2*ky+1)   = 0.0
           ans (2*kx+1, 2*ky+1) = 0.0
        else
           quot = (real(kx**2)*scalex + real(ky**2)*scaley)
           ans (2*kx, 2*ky)     =  f (2*kx, 2*ky) / quot
           ans (2*kx+1, 2*ky)   =  f (2*kx+1, 2*ky) / quot
           ans (2*kx, 2*ky+1)   =  f (2*kx, 2*ky+1) / quot
           ans (2*kx+1, 2*ky+1) =  f (2*kx+1, 2*ky+1) / quot
        end if
     end do
  end do
  !$omp end parallel do
end subroutine poisson

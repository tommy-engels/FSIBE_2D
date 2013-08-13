subroutine vis ( dt, work)
!---------------------------------------------------------------
!  Calculate viscous term for time advancement
!  exp (-nu*k^2*dt)
!  Optimized for vectorization (Dmitry, Feb 1, 2008)
!---------------------------------------------------------------
  use share_vars
  implicit none
  integer :: kx, kx_max, ky, ky_max
  real (kind=pr), dimension (0:nx-1, 0:ny-1), intent(out) :: work
  real (kind=pr), intent (in) :: dt
  real (kind=pr) :: coefx, coefy

  kx_max = (nx/2-1) 
  ky_max = (ny/2-1)
  coefx = - dt * nu * scalex
  coefy = - dt * nu * scaley

  !$omp parallel do private(kx,ky)
  do ky = 0, ky_max
    do kx = 0, kx_max
        work (2*kx, 2*ky) = exp( coefx * real(kx**2) + coefy * real(ky**2) )
    enddo
  enddo
  !$omp end parallel do 


  !$omp parallel do private(kx,ky)
  do ky=0,ny/2-1
     work (1:nx-1:2,2*ky) = work (0:nx-2:2,2*ky)
     work (0:nx-2:2,2*ky+1) = work (0:nx-2:2,2*ky)
     work (1:nx-1:2,2*ky+1) = work (0:nx-2:2,2*ky)
  end do
  !$omp end parallel do 

end subroutine vis

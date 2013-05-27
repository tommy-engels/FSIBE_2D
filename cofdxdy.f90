subroutine cofdxdy (fk, fk_dxdy)
!---------------------------------------------------------------
!     Calculation of d/dx in the Fourier-space ==> -kx^2
!     for the first index
!     Scaling included
!     Ck = Ak + i Bk with Ak= Fk(2k,l) and Bk= Fk(2k+1,l)
!                         for k=0, KX-1 and all l=0,NY-1
!---------------------------------------------------------------
  use share_vars
  implicit none
  integer kx, kx_max, ky, ky_max
  real (kind=pr), dimension (0:nx-1, 0:ny-1), intent (in) :: fk
  real (kind=pr), dimension (0:nx-1, 0:ny-1), intent (out) :: fk_dxdy 
  real (kind=pr) :: fac, scale

  scale = (2.0*pi/xl) * (2.0*pi/yl)
  kx_max = (nx/2-1) 
  ky_max = (ny/2-1)

  !$omp parallel do private(kx,ky)
  do ky = 0, ky_max
     do kx = 0, kx_max
        fac = - real(kx) * real (ky) * scale
        fk_dxdy (2*kx, 2*ky)     = fac * fk (2*kx+1, 2*ky+1)
        fk_dxdy (2*kx+1, 2*ky)   = - fac * fk (2*kx, 2*ky+1)
        fk_dxdy (2*kx, 2*ky+1)   = - fac * fk (2*kx+1, 2*ky)
        fk_dxdy (2*kx+1, 2*ky+1) = fac * fk (2*kx, 2*ky)
     end do
  end do
  !$omp end parallel do

end subroutine



!subroutine cofdxdy_big (fk, fk_dxdy)
!!---------------------------------------------------------------
!!     Calculation of d/dx in the Fourier-space ==> -kx^2
!!     for the first index
!!     Scaling included
!!     Ck = Ak + i Bk with Ak= Fk(2k,l) and Bk= Fk(2k+1,l)
!!                         for k=0, KX-1 and all l=0,NY-1
!!---------------------------------------------------------------
!  use share_vars
!  implicit none
!  integer kx, kx_max, ky, ky_max,nx_tmp,ny_tmp
!  real (kind=pr), dimension (0:up*nxs-1, 0:up*nys-1), intent (in) :: fk
!  real (kind=pr), dimension (0:up*nxs-1, 0:up*nys-1), intent (out) :: fk_dxdy 
!  real (kind=pr) :: fac, scale

!  scale = (2.0*pi/xl) * (2.0*pi/yl)
  
!      ! LAZY!
!  nx_tmp = nx
!  ny_tmp = ny
!  nx=up*nxs
!  ny=up*nys
  
!  kx_max = (nx/2-1) 
!  ky_max = (ny/2-1)

!  !$omp parallel do private(kx,ky)
!  do ky = 0, ky_max
!     do kx = 0, kx_max
!        fac = - real(kx) * real (ky) * scale
!        fk_dxdy (2*kx, 2*ky)     = fac * fk (2*kx+1, 2*ky+1)
!        fk_dxdy (2*kx+1, 2*ky)   = - fac * fk (2*kx, 2*ky+1)
!        fk_dxdy (2*kx, 2*ky+1)   = - fac * fk (2*kx+1, 2*ky)
!        fk_dxdy (2*kx+1, 2*ky+1) = fac * fk (2*kx, 2*ky)
!     end do
!  end do
!  !$omp end parallel do
!  nx=nx_tmp
!  ny=ny_tmp
!end subroutine

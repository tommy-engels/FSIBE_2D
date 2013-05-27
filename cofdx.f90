subroutine cofdx (fk, fk_dx)
!---------------------------------------------------------------
!     calculation of d/dx in the fourier-space ==> *ik
!     for the first index
!     scaling included
!     ck = ak + i bk with ak= fk(2k,l) and bk= fk(2k+1,l)
!                         for k=0, kx-1 and all l=0,ny-1
!---------------------------------------------------------------
  use share_vars
  implicit none
  integer kx, ky
  real (kind=pr), dimension (0:nx-1, 0:ny-1), intent (in) :: fk   
  real (kind=pr), dimension (0:nx-1, 0:ny-1), intent (out) :: fk_dx   
  real (kind=pr) :: scale1

  scale1 = 2.0*pi/xl

  !$omp parallel do private(kx,ky)
  do ky = 0, ny-1
     do kx = 0, nx-2, 2
        fk_dx (kx, ky) = -fk (kx+1, ky) * real (kx/2) * scale1
     end do  
     do kx = 1, nx-1, 2
        fk_dx (kx, ky) = fk (kx-1, ky) * real ((kx-1)/2) * scale1
     end do
  end do
  !$omp end parallel do

end subroutine cofdx



!subroutine cofdx_big (fk, fk_dx)
!!---------------------------------------------------------------
!!     calculation of d/dx in the fourier-space ==> *ik
!!     for the first index
!!     scaling included
!!     ck = ak + i bk with ak= fk(2k,l) and bk= fk(2k+1,l)
!!                         for k=0, kx-1 and all l=0,ny-1
!!---------------------------------------------------------------
!  use share_vars
!  implicit none
!  integer kx, ky,nx_tmp,ny_tmp
!  real (kind=pr), dimension (0:up*nxs-1, 0:up*nys-1), intent (in) :: fk   
!  real (kind=pr), dimension (0:up*nxs-1, 0:up*nys-1), intent (out) :: fk_dx   
!  real (kind=pr) :: scale1

!  scale1 = 2.0*pi/xl
!      ! LAZY!
!  nx_tmp = nx
!  ny_tmp = ny
!  nx=up*nxs
!  ny=up*nys

!  !$omp parallel do private(kx,ky)
!  do ky = 0, ny-1
!     do kx = 0, nx-2, 2
!        fk_dx (kx, ky) = -fk (kx+1, ky) * real (kx/2) * scale1
!     end do  
!     do kx = 1, nx-1, 2
!        fk_dx (kx, ky) = fk (kx-1, ky) * real ((kx-1)/2) * scale1
!     end do
!  end do
!  !$omp end parallel do
!  nx=nx_tmp
!  ny=ny_tmp
!end subroutine cofdx_big

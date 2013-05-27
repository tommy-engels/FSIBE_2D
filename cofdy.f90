subroutine cofdy (fk, fk_dy)
!---------------------------------------------------------------
!     calculation of d/dy in the fourier-space ==> *ik
!     for the second index
!     scaling included
!     ck = ak + i bk with ak= fk(k,2l) and bk= fk(k,2l+1)
!                         for l=0, ky-1 and all k=0,nx-1
!---------------------------------------------------------------
  use share_vars
  implicit none
  integer :: kx, ky
  real (kind=pr), dimension (0:nx-1, 0:ny-1), intent (in) :: fk   
  real (kind=pr), dimension (0:nx-1, 0:ny-1), intent (out) :: fk_dy   
  real (kind=pr) :: scale1

  scale1 = 2.0*pi/yl

  !$omp parallel do private(kx,ky)
  do kx = 0, nx-1
     do ky = 0, ny-2, 2 ! loop over real parts of the FFT signal
        fk_dy (kx, ky) = -fk (kx, ky+1) * real (ky/2) * scale1
     end do
     do ky = 1, ny-1, 2 ! loop over imaginary parts of the FFT
      ! actually, both (real/imag) have of course the same wavenumber. 
        fk_dy (kx, ky) = fk (kx, ky-1) * real ((ky-1)/2) * scale1
     end do
  end do
  !$omp end parallel do

end subroutine


!subroutine cofdy_big (fk, fk_dy)
!!---------------------------------------------------------------
!!     calculation of d/dy in the fourier-space ==> *ik
!!     for the second index
!!     scaling included
!!     ck = ak + i bk with ak= fk(k,2l) and bk= fk(k,2l+1)
!!                         for l=0, ky-1 and all k=0,nx-1
!!---------------------------------------------------------------
!  use share_vars
!  implicit none
!  integer :: kx, ky, nx_tmp,ny_tmp
!  real (kind=pr), dimension (0:up*nxs-1, 0:up*nys-1), intent (in) :: fk   
!  real (kind=pr), dimension (0:up*nxs-1, 0:up*nys-1), intent (out) :: fk_dy   
!  real (kind=pr) :: scale1

!  scale1 = 2.0*pi/yl
  
!    ! LAZY!
!  nx_tmp = nx
!  ny_tmp = ny
!  nx=up*nxs
!  ny=up*nys

!  !$omp parallel do private(kx,ky)
!  do kx = 0, nx-1
!     do ky = 0, ny-2, 2 ! loop over real parts of the FFT signal
!        fk_dy (kx, ky) = -fk (kx, ky+1) * real (ky/2) * scale1
!     end do
!     do ky = 1, ny-1, 2 ! loop over imaginary parts of the FFT
!      ! actually, both (real/imag) have of course the same wavenumber. 
!        fk_dy (kx, ky) = fk (kx, ky-1) * real ((ky-1)/2) * scale1
!     end do
!  end do
!  !$omp end parallel do
!  nx=nx_tmp
!  ny=ny_tmp

!end subroutine

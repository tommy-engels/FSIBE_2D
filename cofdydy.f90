subroutine cofdydy (fk, fk_dydy)
!---------------------------------------------------------------
!     Calculation of d/dx in the Fourier-space ==> -kx^2
!     for the first index
!     Scaling included
!     Ck = Ak + i Bk with Ak= Fk(2k,l) and Bk= Fk(2k+1,l)
!                         for k=0, KX-1 and all l=0,NY-1
!---------------------------------------------------------------
  use share_vars
  implicit none
  integer kx, ky
  real (kind=pr), dimension (0:nx-1, 0:ny-1), intent (in)  :: fk
  real (kind=pr), dimension (0:nx-1, 0:ny-1), intent (out) :: fk_dydy 
  real (kind=pr) :: scale

  scale = (2.0*pi/yl)**2

  !$omp parallel do private(kx,ky)
  do kx = 0, nx-1
     do ky = 0, ny-2, 2
        fk_dydy (kx, ky)   = - fk (kx, ky)   * real (ky/2)**2 * scale
        fk_dydy (kx, ky+1) = - fk (kx, ky+1) * real (ky/2)**2 * scale
     end do
  end do
  !$omp end parallel do

end subroutine

subroutine cofdxdx (fk, fk_dxdx)
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
  real (kind=pr), dimension (0:nx-1, 0:ny-1), intent (in) :: fk
  real (kind=pr), dimension (0:nx-1, 0:ny-1), intent (out) :: fk_dxdx 
  real (kind=pr) :: scale

  scale = (2.0*pi/xl)**2

  !$omp parallel do private(kx,ky)
  do ky = 0, ny-1
     do kx = 0, nx-2, 2
        fk_dxdx (kx, ky)   = - fk (kx, ky)   * real (kx/2)**2 * scale
        fk_dxdx (kx+1, ky) = - fk (kx+1, ky) * real (kx/2)**2 * scale
     end do
  end do
  !$omp end parallel do

end subroutine

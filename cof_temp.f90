!====================================================================
!====================================================================
!
!        Default Fourier transform subroutines
!
!====================================================================
!====================================================================


subroutine fft_initialize
  use share_vars
  implicit none

  allocate ( trigsx (3*nx/2+1) )
  allocate ( trigsy (3*ny/2+1) )

  call set99 (trigsx, ifaxx, nx)
  call set99 (trigsy, ifaxy, ny)

end subroutine fft_initialize



subroutine fft_free
  use share_vars
  implicit none

  deallocate (trigsx, trigsy)

end subroutine fft_free



subroutine coftx (f, fk)
!====================================================================
!     Calculation of the Fourier-coefficients of a real function
!     with the FFT de Temperton / TRANSFORMATION at first INDEX !!!
!     FILTERING
!====================================================================
  use share_vars
  implicit none
  integer, parameter :: isign = -1, jump = 1
  integer :: kx, ky
  real (kind=pr), dimension (0:ny-1, 0:nx+1) :: ft
  real (kind=pr), dimension (2*(nx+2)*ny)  :: work
  real (kind=pr), dimension (0:nx-1, 0:ny-1), intent (in) ::  f  
  real (kind=pr), dimension (0:nx-1, 0:ny-1), intent (out) ::  fk  

!   FFT Plan physique ---> plan spectral
  do ky = 0, ny-1
     do kx = 0, nx-1
        ft (ky, kx) = f (kx, ky)
     end do
  end do
  do ky = 0, ny-1
     ft(ky,nx)   = 0.0
     ft(ky,nx+1) = 0.0
  end do
  call fft991 (ft, work, trigsx, ifaxx, ny, jump, nx, ny, isign)
  do ky = 0, ny-1
     do kx = 0, nx-1
        fk (kx, ky) = ft (ky, kx)
     end do
  end do
!      last mode M=NX, NX+1; mode KF left unconsidered => filtering
end subroutine coftx



subroutine cofty (f, fk)
!====================================================================
!     calculation of the fourier-coefficients of a real function
!     with the fft de temperton / transformation at second index
!     filtering
!====================================================================
  use share_vars
  implicit none
  integer, parameter :: isign = -1, jump = 1
  integer :: kx, ky
  real (kind=pr), dimension (0:nx-1, 0:ny+1) :: ft
  real (kind=pr), dimension (2*(ny+2)*nx)  :: work
  real (kind=pr), dimension (0:nx-1, 0:ny-1), intent (in) ::  f  
  real (kind=pr), dimension (0:nx-1, 0:ny-1), intent (out) ::  fk  

!   fft plan physique ---> plan spectral
  do ky = 0, ny-1
     do kx = 0, nx-1
        ft (kx, ky) = f (kx, ky)
     end do
  end do
  do kx = 0, nx-1
     ft(kx,ny)   = 0.0
     ft(kx,ny+1) = 0.0
  end do
  call fft991 (ft, work, trigsy, ifaxy, nx, jump, ny, nx, isign)
  do ky = 0, ny-1
     do kx = 0, nx-1
        fk (kx, ky) = ft (kx, ky)
     end do
  end do
!     last mode m=ny, ny+1; mode kf left unconsidered => filtering
end subroutine



subroutine cofitx (fk, f)
!====================================================================
!     calculation of a real function from its fourier-coefficients
!     with the fft de temperton / transformation sur premiere index !!!
!     filtering
!====================================================================
  use share_vars
  implicit none
  integer, parameter :: isign = 1, jump = 1
  integer :: kx, ky
  real (kind=pr), dimension (0:ny-1, 0:nx+1) :: ft
  real (kind=pr), dimension (2*(nx+2)*ny) :: work
  real (kind=pr), dimension (0:nx-1, 0:ny-1), intent (out) :: f  
  real (kind=pr), dimension (0:nx-1, 0:ny-1), intent (in)  :: fk 

!   fft plan spectral ---> plan physique
  do ky = 0, ny-1
     do kx = 0, nx-1
        ft (ky, kx) = fk (kx, ky)
     end do
  end do
!     filtering
  do ky = 0, ny-1
     do kx = nx, nx+1
        ft (ky, kx) = 0.0
     end do
  end do
  call fft991 (ft, work, trigsx, ifaxx, ny, jump, nx, ny, isign)
  do ky = 0, ny-1
     do kx = 0, nx-1
        f (kx, ky) = ft (ky, kx)
     end do
  end do
end subroutine cofitx



subroutine cofity (fk,f)
!====================================================================
!     calculation of a real function from its fourier-coefficients
!     with the fft de temperton / transformation sur deuxieme index !!!
!     filtering
!====================================================================
  use share_vars
  implicit none
  integer, parameter :: isign = 1, jump = 1
  integer kx, ky
  real (kind=pr), dimension (0:nx-1, 0:ny+1) :: ft
  real (kind=pr), dimension (2*(ny+2)*nx)   :: work
  real (kind=pr), dimension (0:nx-1, 0:ny-1), intent (out) :: f  
  real (kind=pr), dimension (0:nx-1, 0:ny-1), intent (in)  :: fk 

!   fft plan spectral ---> plan physique
  do ky = 0, ny-1
     do kx = 0, nx-1
        ft (kx, ky) = fk (kx, ky)
     end do
  end do
!     filtering
  do ky = ny, ny+1
     do kx = 0, nx-1
        ft (kx, ky) = 0.0
     end do
  end do
  call fft991 (ft,work,trigsy,ifaxy,nx,jump,ny,nx,isign)
  do ky = 0, ny-1
     do kx = 0, nx-1
        f (kx, ky) = ft (kx, ky)
     end do
  end do
end subroutine cofity


subroutine coftxy (f, fk)
!====================================================================
!     Calculation of the Fourier-coefficients of a real function
!     with the Temperton FFT
!     In x and y directions
!     FILTERING
!====================================================================
  use share_vars
  implicit none
  real (kind=pr), dimension (0:nx-1, 0:ny-1) :: ft
  real (kind=pr), dimension (0:nx-1, 0:ny-1), intent (in) ::  f  
  real (kind=pr), dimension (0:nx-1, 0:ny-1), intent (out) ::  fk  

  call coftx (f, ft)
  call cofty (ft, fk)

end subroutine coftxy



subroutine cofitxy (fk, f)
!====================================================================
!     calculation of a real function from its fourier-coefficients
!     with the temperton fft
!     In x and y directions
!     filtering
!====================================================================
  use share_vars
  implicit none
  real (kind=pr), dimension (0:nx-1, 0:ny-1) :: ft
  real (kind=pr), dimension (0:nx-1, 0:ny-1), intent (out) :: f  
  real (kind=pr), dimension (0:nx-1, 0:ny-1), intent (in)  :: fk 

  call cofity (fk, ft)
  call cofitx (ft, f)

end subroutine cofitxy



subroutine cofts (f, fk, L, n)
!====================================================================
!     Calculation of the Fourier-coefficients of a real function
!     with the FFT de Temperton / TRANSFORMATION at first INDEX !!!
!     FILTERING
!====================================================================
  use share_vars
  implicit none
  integer, parameter :: isign = -1, jump = 1
  integer, intent (in) :: L, n 
  integer :: kL, kn
  integer, dimension (10) :: ifaxs
  real (kind=pr), dimension (0:n-1, 0:L+1) :: ft
  real (kind=pr), dimension (2*(L+2)*n)  :: work
  real (kind=pr), dimension (0:L-1, 0:n-1), intent (in) ::  f  
  real (kind=pr), dimension (0:L-1, 0:n-1), intent (out) ::  fk
  real (kind=pr), dimension (:), allocatable :: trigss  

  ! Initialize 1D FFT
  allocate ( trigss (3*L/2+1) )
  call set99 (trigss, ifaxs, L)

  do kn = 0, n-1
     do kL = 0, L-1
        ft (kn, kL) = f (kL, kn)
     end do
  end do
  do kn = 0, n-1
     ft(kn,L)   = 0.0
     ft(kn,L+1) = 0.0
  end do

  !   FFT Plan physique ---> plan spectral
  call fft991 (ft, work, trigss, ifaxs, n, jump, L, n, isign)

  do kn = 0, n-1
     do kL = 0, L-1
        fk (kL, kn) = ft (kn, kL)
     end do
  end do

  ! Deallocate memory
  deallocate (trigss)

!      last mode M=NX, NX+1; mode KF left unconsidered => filtering
end subroutine cofts



subroutine cofits (fk, f, L, n)
!====================================================================
!     calculation of a real function from its fourier-coefficients
!     with the fft de temperton / transformation sur premiere index !!!
!     filtering
!====================================================================
  use share_vars
  implicit none
  integer, parameter :: isign = 1, jump = 1
  integer, intent (in) :: L, n
  integer :: kn, kL
  integer, dimension (10) :: ifaxs
  real (kind=pr), dimension (0:n-1, 0:L+1) :: ft
  real (kind=pr), dimension (2*(L+2)*n) :: work
  real (kind=pr), dimension (0:L-1, 0:n-1), intent (out) :: f  
  real (kind=pr), dimension (0:L-1, 0:n-1), intent (in)  :: fk 
  real (kind=pr), dimension (:), allocatable :: trigss

  ! Initialize 1D FFT
  allocate ( trigss (3*L/2+1) )
  call set99 (trigss, ifaxs, L)

  do kn = 0, n-1
     do kL = 0, L-1
        ft (kn, kL) = fk (kL, kn)
     end do
  end do
!     filtering
  do kn = 0, n-1
     do kL = L, L+1
        ft (kn, kL) = 0.0
     end do
  end do

  !   fft plan spectral ---> plan physique
  call fft991 (ft, work, trigss, ifaxs, n, jump, L, n, isign)

  do kn = 0, n-1
     do kL = 0, L-1
        f (kL, kn) = ft (kn, kL)
     end do
  end do

  ! Deallocate memory
  deallocate (trigss)

end subroutine cofits

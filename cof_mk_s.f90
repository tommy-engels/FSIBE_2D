!====================================================================
!====================================================================
!
! Fourier transform subroutines using MathKeisan library on IDRIS NEC
!
!====================================================================
!====================================================================


subroutine fft_initialize
! Not required for MathKeisan library
end subroutine fft_initialize



subroutine fft_free
! Not required for MathKeisan library
end subroutine fft_free



subroutine coftx (f, fk)
!====================================================================
!     Calculation of the Fourier-coefficients of a real function
!     along x (1st index)
!     FILTERING
!====================================================================
  use share_vars
  implicit none
  integer, parameter :: iopt = 1, jump = 1
  integer :: ierr
  real (kind=pr), dimension (0:nx+1, 0:ny-1) :: ft
  real (kind=pr), dimension (0:nx-1, 0:ny-1), intent (in) ::  f  
  real (kind=pr), dimension (0:nx-1, 0:ny-1), intent (out) ::  fk  
  real (kind=pr) :: norm

!   FFT Plan physique ---> plan spectral

  !$omp parallel workshare
  ft(0:nx-1,:) = f
  ft(nx:nx+1,:) = 0.0
  !$omp end parallel workshare
 
  call srcfts (ft, nx, jump, ny, nx+2, iopt, ierr)

  norm = 1.0 / real(nx)

  !$omp parallel workshare
  fk = ft(0:nx-1,:) * norm
  !$omp end parallel workshare

!      last mode M=NX, NX+1; mode KF left unconsidered => filtering
end subroutine coftx



subroutine cofitx (fk, f)
!====================================================================
!     Calculation of a real function from its Fourier coefficients
!     along x (1st index)
!     FILTERING
!====================================================================
  use share_vars
  implicit none
  integer, parameter :: iopt = -1, jump = 1
  integer :: ierr
  real (kind=pr), dimension (0:nx+1, 0:ny-1) :: ft
  real (kind=pr), dimension (0:nx-1, 0:ny-1), intent (in) ::  fk  
  real (kind=pr), dimension (0:nx-1, 0:ny-1), intent (out) ::  f  
  real (kind=pr) :: norm

!   Inverse FFT 

  norm = real(nx)

  !$omp parallel workshare
  ft(0:nx-1,:) = fk * norm
  ft(nx:nx+1,:) = 0.0
  !$omp end parallel workshare
 
  call srcfts (ft, nx, jump, ny, nx+2, iopt, ierr)

  !$omp parallel workshare
  f = ft(0:nx-1,:)
  !$omp end parallel workshare

!      last mode M=NX, NX+1; mode KF left unconsidered => filtering
end subroutine cofitx



subroutine cofty (f, fk)
!====================================================================
!     Calculation of the Fourier-coefficients of a real function
!     along y (2nd index)
!     FILTERING
!====================================================================
  use share_vars
  implicit none
  integer, parameter :: iopt = 1, jump = 1
  integer :: ierr
  real (kind=pr), dimension (0:nx-1, 0:ny+1) :: ft
  real (kind=pr), dimension (0:nx-1, 0:ny-1), intent (in) ::  f  
  real (kind=pr), dimension (0:nx-1, 0:ny-1), intent (out) ::  fk  
  real (kind=pr) :: norm

!   FFT Plan physique ---> plan spectral

  !$omp parallel workshare
  ft(:,0:ny-1) = f
  ft(:,ny:ny+1) = 0.0
  !$omp end parallel workshare
 
  call srcfts (ft, ny, nx, nx, jump, iopt, ierr)

  norm = 1.0 / real(ny)

  !$omp parallel workshare
  fk = ft(:,0:ny-1) * norm
  !$omp end parallel workshare

!      last mode M=NY, NY+1; mode KF left unconsidered => filtering
end subroutine cofty



subroutine cofity (fk, f)
!====================================================================
!     Calculation of a real function from its Fourier coefficients
!     along y (2nd index)
!     FILTERING
!====================================================================
  use share_vars
  implicit none
  integer, parameter :: iopt = -1, jump = 1
  integer :: ierr
  real (kind=pr), dimension (0:nx-1, 0:ny+1) :: ft
  real (kind=pr), dimension (0:nx-1, 0:ny-1), intent (in) ::  fk  
  real (kind=pr), dimension (0:nx-1, 0:ny-1), intent (out) ::  f  
  real (kind=pr) :: norm

!   Inverse FFT 

  norm = real(ny)

  !$omp parallel workshare
  ft(:,0:ny-1) = fk * norm
  ft(:,ny:ny+1) = 0.0
  !$omp end parallel workshare
 
  call srcfts (ft, ny, nx, nx, jump, iopt, ierr)

  !$omp parallel workshare
  f = ft(:,0:ny-1)
  !$omp end parallel workshare

!      last mode M=NY, NY+1; mode KF left unconsidered => filtering
end subroutine cofity



subroutine coftxy (f, fk)
!====================================================================
!     Calculation of the Fourier-coefficients of a real function
!     along x (1st index) and y (2nd index)
!     FILTERING
!     MK library is used
!====================================================================
  use share_vars
  implicit none
  integer, parameter :: iopt = 1
  integer :: ierr
  complex (kind=pr), dimension (0:nx, 0:ny-1) :: ft
  real (kind=pr), dimension (0:nx-1, 0:ny-1), intent (in) ::  f  
  real (kind=pr), dimension (0:nx-1, 0:ny-1), intent (out) ::  fk  
  real (kind=pr) :: norm

!   FFT Plan physique ---> plan spectral

  !$omp parallel workshare
  ft(0:nx-1,:) = cmplx( f, 0.0, pr )
  ft(nx,:) = 0.0
  !$omp end parallel workshare
 
  ! call real-to-complex 2d fft
  ! complex array is used for storage (redundant but faster)
  call crc2ft (ft, nx, ny, nx+1, iopt, ierr)

  ! normalization
  norm = 1.0 / ( real(nx) * real(ny) )

  ! complex 2d Fourier coefficients are used to 
  ! compute coefficients before sin and cos (Temperton style)

  !$omp parallel workshare
  fk(0:nx-2:2,0) = real( ft(0:nx/2-1,0) )
  fk(1:nx-1:2,0) = imag( ft(0:nx/2-1,0) )

  fk(:,1) = 0.0
  
  fk(0:nx-2:2,2:ny-2:2) = 0.5 * ( real( ft(0:nx/2-1,1:ny/2-1)) &
                         + real( ft(0:nx/2-1,ny-1:ny/2+1:-1)) )
  fk(0:nx-2:2,3:ny-1:2) = 0.5 * ( imag( ft(0:nx/2-1,1:ny/2-1)) &
                         - imag( ft(0:nx/2-1,ny-1:ny/2+1:-1)) )
  fk(1:nx-1:2,2:ny-2:2) = 0.5 * ( imag( ft(0:nx/2-1,1:ny/2-1)) &
                         + imag( ft(0:nx/2-1,ny-1:ny/2+1:-1)) )
  fk(1:nx-1:2,3:ny-1:2) = 0.5 * ( - real( ft(0:nx/2-1,1:ny/2-1)) &
                         + real( ft(0:nx/2-1,ny-1:ny/2+1:-1)) )

  fk = fk * norm
  !$omp end parallel workshare

!      mode KF left unconsidered => filtering
end subroutine coftxy



subroutine cofitxy (fk, f)
!====================================================================
!     Calculation of a real function from its Fourier coefficients
!     along x (1st index) and y (2nd index)
!     FILTERING
!     MK library is used
!====================================================================
  use share_vars
  implicit none
  integer, parameter :: iopt = -1
  integer :: ierr
  complex (kind=pr), dimension (0:nx, 0:ny-1) :: ft
  real (kind=pr), dimension (0:nx-1, 0:ny-1), intent (in) ::  fk  
  real (kind=pr), dimension (0:nx-1, 0:ny-1), intent (out) ::  f  
  real (kind=pr) :: norm

!   Inverse FFT 

  norm = real(nx) * real(ny)

  !$omp parallel workshare
  ft = 0.0

  ft(0:nx/2-1,0) = cmplx( fk(0:nx-2:2,0), fk(1:nx-1:2,0), pr )

  ft(0:nx/2-1,1:ny/2-1) = cmplx( fk(0:nx-2:2,2:ny-2:2) - fk(1:nx-1:2,3:ny-1:2), &
                       fk(0:nx-2:2,3:ny-1:2) + fk(1:nx-1:2,2:ny-2:2) )

  ft(0:nx/2-1,ny-1:ny/2+1:-1) = cmplx( fk(0:nx-2:2,2:ny-2:2) + fk(1:nx-1:2,3:ny-1:2), &
                      -fk(0:nx-2:2,3:ny-1:2) + fk(1:nx-1:2,2:ny-2:2) )

  ft = ft * norm
  !$omp end parallel workshare
 
  call crc2ft (ft, nx, ny, nx+1, iopt, ierr)

  !$omp parallel workshare
  f = real( ft(0:nx-1,:) )
  !$omp end parallel workshare

!      mode KF left unconsidered => filtering
end subroutine cofitxy



subroutine cofts (f, fk, L, n)
!====================================================================
!     Calculation of the Fourier coefficients of n real functions
!     aligned in columns
!     FILTERING
!     MK library is used
!====================================================================
  use share_vars
  implicit none
  integer, parameter :: iopt = 1, jump = 1
  integer, intent (in) :: l, n
  integer :: ierr
  real (kind=pr), dimension (0:L+1, 0:n-1) :: ft
  real (kind=pr), dimension (0:L-1, 0:n-1), intent (in) ::  f  
  real (kind=pr), dimension (0:L-1, 0:n-1), intent (out) ::  fk  
  real (kind=pr) :: norm

!   FFT Plan physique ---> plan spectral

  !$omp parallel workshare
  ft(0:L-1,:) = f
  ft(L:L+1,:) = 0.0
  !$omp end parallel workshare
 
  call srcfts (ft, L, jump, n, L+2, iopt, ierr)

  norm = 1.0 / real(L)

  !$omp parallel workshare
  fk = ft(0:L-1,:) * norm
  !$omp end parallel workshare

!      last mode M=L, L+1; mode KF left unconsidered => filtering
end subroutine cofts



subroutine cofits (fk, f, L, n)
!====================================================================
!     Calculation of n real functions from their Fourier coefficients
!     aligned in columns
!     FILTERING
!     MK library is used
!====================================================================
  use share_vars
  implicit none
  integer, parameter :: iopt = -1, jump = 1
  integer, intent (in) :: L, n ! length of the array and number of arrays
  integer :: ierr
  real (kind=pr), dimension (0:L+1, 0:n-1) :: ft
  real (kind=pr), dimension (0:L-1, 0:n-1), intent (in) ::  fk  
  real (kind=pr), dimension (0:L-1, 0:n-1), intent (out) ::  f  
  real (kind=pr) :: norm

!   Inverse FFT 

  norm = real(L)

  !$omp parallel workshare
  ft(0:L-1,:) = fk * norm
  ft(L:L+1,:) = 0.0
  !$omp end parallel workshare
 
  call srcfts (ft, L, jump, n, L+2, iopt, ierr)

  !$omp parallel workshare
  f = ft(0:L-1,:)
  !$omp end parallel workshare

!      last mode M=L, L+1; mode KF left unconsidered => filtering
end subroutine cofits

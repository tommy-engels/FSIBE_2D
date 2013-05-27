module Interpolation
  implicit none
  contains
  
  
real (kind=pr) function DeltaInterpolation (x_target, y_target, field2, x1_box, y1_box, x2_box, y2_box,typ )
  use share_vars
  implicit none
!  interpolation with delta kernels. see also other interpolation routines for comments on spacing,  domain size and so on..
  integer :: i,j,ix,iy, N_support
  real (kind=pr) :: x,y,x_local,y_local,x_2,y_2,dx, dy, R1,R2
  real (kind=pr), intent (in) :: field2(0:,0:), x_target, y_target, x1_box, y1_box, x2_box, y2_box
  character*(*), intent (in) :: typ

  dx = (x2_box-x1_box)/real(size(field2,1)-1 )
  dy = (y2_box-y1_box)/real(size(field2,2)-1 )


  if ( (x_target > x2_box).or.(x_target < x1_box).or.(y_target > y2_box).or.(y_target < y1_box) ) then
    write(*,'("target: (",es11.4,"|",es11.4,") but box: (",es11.4,"|",es11.4,") x (",es11.4,"|",es11.4,")")') &
    x_target, y_target, x1_box, y1_box, x2_box, y2_box
    write (*,*) "!!! DeltaInterpolation: target coordinates not in the field."
    stop
  endif
  
  ! i,j are the coordinates of the lower left point of the grid containing the target point  
  i = int((x_target-x1_box)/dx)  ! attention on index shift because of automatic array
  j = int((y_target-y1_box)/dy)
  
  R1=0.0
  
  N_support = 4
  if (typ=="gaussian_3") N_support=15
  if (typ=="gaussian_4") N_support=15 !very expensive support!
  
  do ix=i-N_support,i+N_support ! the box size around the point
  do iy=j-N_support,j+N_support
    x_local= real(ix)*dx + x1_box !coordinates of the point we're looking at
    y_local= real(iy)*dy + y1_box
    R1 = R1 + delta(abs(x_local - x_target),dx,typ)*delta(abs(y_local - y_target),dy,typ)*field2(ix,iy)
  enddo
  enddo

  DeltaInterpolation = R1

  return

end function DeltaInterpolation


real (kind=pr) function delta(x,dx,typ)
  use share_vars
  real(kind=pr), intent(in) :: x,dx
  real(kind=pr) :: r,s
  character*(*), intent (in) :: typ
  ! ----------------------------------
  !  This function returns different delta kernels
  ! ----------------------------------
  
  r=abs(x/dx)  
       
  select case (typ)
      case ("gaussian_1")   
	  delta=exp(-0.5*r**2)/sqrt(2.0*pi)
      case ("gaussian_2")   
          s=0.75
	  delta=exp(-0.5*(r/s)**2)/sqrt(2.0*pi)/s
      case ("gaussian_3")   
          s=2.0
	  delta=exp(-0.5*(r/s)**2)/sqrt(2.0*pi)/s
      case ("gaussian_4")   
          s=3.0
	  delta=exp(-0.5*(r/s)**2)/sqrt(2.0*pi)/s
      case ("phi_4norm")  ! phi 4 not smoothed (see JCP 228)
	  if (r<1.0) then
	      delta = (1./8.)*(3.-2.*r+sqrt(1.+4.*r-4.*r**2))
	  elseif ( (r>=1.0) .and. (r<=2.0)     ) then
	      delta = (1./8.)*(5.-2.*r-sqrt(-7.+12.*r-4.*r**2))
	  elseif ( (r>=2.0)    ) then
	      delta = 0.0
	  endif	  
      case ("phi_4star") ! see Yang, Zhang, Li: A smoothing technique for discrete delta functions [...] JCP 228 (2009)
	  if (r<0.5) then
	      delta = (3./8.)+(pi/32.)-0.25*r**2
	  elseif ( (r>=0.5) .and. (r<=1.5)   ) then
	      delta = 0.25 + (1-r)/8.  *sqrt(-2.+8.*r-4.*r**2) -asin(sqrt(2.)*(r-1.))/8.
	  elseif ( (r>=1.5) .and. (r<=2.5)   ) then
	      delta = (17./16.) - (pi/64.) - (3.*r/4.) + ((r**2)/8.) + (r-2.)*sqrt(-14.+16.*r-4.*r**2)/16. +asin(sqrt(2.)*(r-2.))/16.
	  elseif ( (r>=2.5)    ) then
	      delta = 0.0;
	  endif  
      case ("peskin") ! Peskin's original cosine function
	  if (r<2.0) then
	      delta = 0.25*(1.0+cos(0.5*pi*r))
	  else
	      delta = 0.0
	  endif
  end select
  
  return
end function
! ----------------------------------------------------------------------------------------


!=================================================================================================================================

real (kind=pr) function LinearInterpolation (x_target, y_target, field2, x1_box, y1_box, x2_box, y2_box )
!  LINEAR Interpolation in a field. The field is of automatic size, indices starting with 0 both. The domain is 
!  defined by x1_box,y1_box and x2_box,y2_box. The target coordinates should lie within that box.
!  NOTE: attention on the upper point of the box. In the rest of the code, which is periodic, the grid is 0:nx-1
!        but the lattice spacing is yl/nx. This means that the point (nx-1) has NOT the coordinate yl but yl-dx
!        (otherwise this point would exist two times!)
!  NOTE3: Coordinates in the box are a constant source for errors. be careful and note that x1_box is NOT ZERO
  use share_vars
  implicit none
  integer :: i,j
  real (kind=pr) :: x,y,x_1,y_1,x_2,y_2,dx, dy, R1,R2
  real (kind=pr), intent (in) :: field2(0:,0:), x_target, y_target, x1_box, y1_box, x2_box, y2_box

  dx = (x2_box-x1_box)/real(size(field2,1)-1 )
  dy = (y2_box-y1_box)/real(size(field2,2)-1 )


  if ( (x_target > x2_box).or.(x_target < x1_box).or.(y_target > y2_box).or.(y_target < y1_box) ) then
    write(*,'("target: (",es11.4,"|",es11.4,") but box: (",es11.4,"|",es11.4,") x (",es11.4,"|",es11.4,")")') &
    x_target, y_target, x1_box, y1_box, x2_box, y2_box
    write (*,*) "!!! LinearInterpolation: target coordinates not in the field."
    stop
  endif

  i=int((x_target-x1_box)/dx)  ! attention on index shift because of automatic array
  j=int((y_target-y1_box)/dy)

  
  x_1= real(i)*dx + x1_box
  y_1= real(j)*dy + y1_box
  x_2= dx*real(i+1) + x1_box
  y_2= dy*real(j+1) + y1_box
  R1 = (x_2-x_target)*field2(i,j)/dx   + (x_target-x_1)*field2(i+1,j)/dx
  R2 = (x_2-x_target)*field2(i,j+1)/dx + (x_target-x_1)*field2(i+1,j+1)/dx

  LinearInterpolation = (y_2-y_target)*R1/dy + (y_target-y_1)*R2/dy

  return

end function LinearInterpolation


!=================================================================================================================================

real (kind=pr) function BicubicInterpolation (x_target, y_target, field,field_dx, field_dy, field_dxdy, x1_box, y1_box, x2_box, y2_box )
!  Bicubic Interpolation in a field. The field is of automatic size, indices starting with 0 both. The domain is 
!  defined by x1_box,y1_box and x2_box,y2_box. The target coordinates should lie within that box.
!  NOTE: attention on the upper point of the box. In the rest of the code, which is periodic, the grid is 0:nx-1
!        but the lattice spacing is yl/nx. This means that the point (nx-1) has NOT the coordinate yl but yl-dx
!        (otherwise this point would exist two times!)
!  NOTE2: This version works with finite Differences for the derivatives. If you want to, you can specify other derivatives 
!         when using BicubicInterpolation.     
!  NOTE3: Coordinates in the box are a constant source for errors. be careful and note that x1_box is NOT ZERO
  use share_vars
  implicit none
  integer :: ix,iy
  real (kind=pr)                  :: f_interp,f_dx_interp,f_dy_interp,dx,dy
  real (kind=pr), intent (in)     :: field(0:,0:),field_dx(0:,0:),field_dy(0:,0:),field_dxdy(0:,0:)
  real (kind=pr), intent (in)     :: x_target, y_target, x1_box, y1_box, x2_box, y2_box
  real (kind=pr), dimension(1:4)  :: y,y1,y2,y12
  real (kind=pr)                  :: x1l,x1u,x2l,x2u,x1,x2

  if ( (x_target > x2_box).or.(x_target < x1_box).or.(y_target > y2_box).or.(y_target < y1_box) ) then
    write(*,'("target: (",es11.4,"|",es11.4,") but box: (",es11.4,"|",es11.4,") x (",es11.4,"|",es11.4,")")') &
    x_target, y_target, x1_box, y1_box, x2_box, y2_box 
    write (*,*) "!!! BicubicInterpolation: target coordinates not in the field."
    stop
  endif

  dx = (x2_box-x1_box)/real(size(field,1)-1)
  dy = (y2_box-y1_box)/real(size(field,2)-1)

  ix=int((x_target-x1_box)/dx)
  iy=int((y_target-y1_box)/dy)

  y(1) = field(ix  ,iy)
  y(2) = field(ix+1,iy)
  y(3) = field(ix+1,iy+1)
  y(4) = field(ix  ,iy+1)

  y1(1) = field_dx(ix  ,iy)
  y1(2) = field_dx(ix+1,iy)
  y1(3) = field_dx(ix+1,iy+1)
  y1(4) = field_dx(ix  ,iy+1)

  y2(1) = field_dy(ix  ,iy)
  y2(2) = field_dy(ix+1,iy)
  y2(3) = field_dy(ix+1,iy+1)
  y2(4) = field_dy(ix  ,iy+1)

  y12(1) = field_dxdy(ix  ,iy)
  y12(2) = field_dxdy(ix+1,iy)
  y12(3) = field_dxdy(ix+1,iy+1)
  y12(4) = field_dxdy(ix  ,iy+1)

  x1l = real(ix)*dx   + x1_box
  x1u = real(ix+1)*dx + x1_box

  x2l = real(iy)*dy   + y1_box
  x2u = real(iy+1)*dy + y1_box

  call bcuint(y,y1,y2,y12,x1l,x1u,x2l,x2u,x_target,y_target, f_interp,f_dx_interp,f_dy_interp)

  BicubicInterpolation = f_interp
  
!  write(*,*) "-------------------------"
!  write (*,'("y  :",4(es11.4,1x))') y
!  write (*,'("y1 :",4(es11.4,1x))') y1
!  write (*,'("y2 :",4(es11.4,1x))') y2
!  write (*,'("y12:",4(es11.4,1x))') y12
!  write (*,'(4(es11.4,3x))') x1l,x1u,x2l,x2u
!  write (*,'(5(es11.4,1x))') x_target,y_target, f_interp,f_dx_interp,f_dy_interp
!  write(*,*) "------------------------------------------"
  
  return
end function BicubicInterpolation

!=================================================================================================================================



real (kind=pr) function BicubicInterpolationFD (x_target, y_target, field, x1_box, y1_box, x2_box, y2_box )
!  Bicubic Interpolation in a field. The field is of automatic size, indices starting with 0 both. The domain is 
!  defined by x1_box,y1_box and x2_box,y2_box. The target coordinates should lie within that box.
!  NOTE: attention on the upper point of the box. In the rest of the code, which is periodic, the grid is 0:nx-1
!        but the lattice spacing is yl/nx. This means that the point (nx-1) has NOT the coordinate yl but yl-dx
!        (otherwise this point would exist two times!)
!  NOTE2: This version works with finite Differences for the derivatives. If you want to, you can specify other derivatives 
!         when using BicubicInterpolation.     
!  NOTE3: Coordinates in the box are a constant source for errors. be careful and note that x1_box is NOT ZERO
  use share_vars
  implicit none
  integer :: ix,iy
  real (kind=pr)                  :: f_interp,f_dx_interp,f_dy_interp,dx,dy
  real (kind=pr), intent (in)     :: field(0:,0:)
  real (kind=pr), intent (in)     :: x_target, y_target, x1_box, y1_box, x2_box, y2_box
  real (kind=pr), dimension(1:4)  :: y,y1,y2,y12
  real (kind=pr)                  :: x1l,x1u,x2l,x2u,x1,x2

  if ( (x_target > x2_box).or.(x_target < x1_box).or.(y_target > y2_box).or.(y_target < y1_box) ) then
    write(*,'("target: (",es11.4,"|",es11.4,") but box: (",es11.4,"|",es11.4,") x (",es11.4,"|",es11.4,")")') &
    x_target, y_target, x1_box, y1_box, x2_box, y2_box 
    write (*,*) "!!! BicubicInterpolation: target coordinates not in the field."
    stop
  endif

  dx = (x2_box-x1_box)/real(size(field,1)-1)
  dy = (y2_box-y1_box)/real(size(field,2)-1)

  ix=int((x_target-x1_box)/dx)
  iy=int((y_target-y1_box)/dy)


  y(1) = field(ix  ,iy)
  y(2) = field(ix+1,iy)
  y(3) = field(ix+1,iy+1)
  y(4) = field(ix  ,iy+1)

  y1(1) = (field(ix+1,iy)-field(ix-1,iy)     )/dx/2.d0    !d/dx
  y1(2) = (field(ix+2,iy)-field(ix,iy)       )/dx/2.d0
  y1(3) = (field(ix+2,iy+1)-field(ix,iy+1)   )/dx/2.d0
  y1(4) = (field(ix+1,iy+1)-field(ix-1,iy+1) )/dx/2.d0

  y2(1) = (field(ix,iy+1)-field(ix,iy-1)    )/dy/2.d0   !d/dy
  y2(2) = (field(ix+1,iy+1)-field(ix+1,iy-1))/dy/2.d0
  y2(3) = (field(ix+1,iy+2)-field(ix+1,iy)  )/dy/2.d0
  y2(4) = (field(ix,iy+2)-field(ix,iy)      )/dy/2.d0

  y12(1) = (field(ix+1,iy+1)-field(ix-1,iy+1)-field(ix+1,iy-1)+field(ix-1,iy-1) )/(4.d0*dx*dy)  !d/dxdy
  y12(2) = (field(ix+2,iy+1)-field(ix,iy+1)  -field(ix+2,iy-1)+field(ix,iy-1)   )/(4.d0*dx*dy)
  y12(3) = (field(ix+2,iy+2)-field(ix,iy+2)  -field(ix+2,iy)  +field(ix,iy)     )/(4.d0*dx*dy)
  y12(4) = (field(ix+1,iy+2)-field(ix-1,iy+2)-field(ix+1,iy)  +field(ix-1,iy)   )/(4.d0*dx*dy)

  x1l = real(ix)*dx   + x1_box
  x1u = real(ix+1)*dx + x1_box

  x2l = real(iy)*dy   + y1_box
  x2u = real(iy+1)*dy + y1_box

  call bcuint(y,y1,y2,y12,x1l,x1u,x2l,x2u, x_target,y_target, f_interp,f_dx_interp,f_dy_interp)

  BicubicInterpolationFD = f_interp
  
!  write(*,*) "-------------------------"
!  write (*,'("y  :",4(es11.4,1x))') y
!  write (*,'("y1 :",4(es11.4,1x))') y1
!  write (*,'("y2 :",4(es11.4,1x))') y2
!  write (*,'("y12:",4(es11.4,1x))') y12
!  write (*,'(4(es11.4,3x))') x1l,x1u,x2l,x2u
!  write (*,'(5(es11.4,1x))') x_target,y_target, f_interp,f_dx_interp,f_dy_interp
!  write(*,*) "------------------------------------------"
!  stop
  return
end function BicubicInterpolationFD

!========================================================================================================================



real (kind=pr) function SplineInterpolation (x_target, y_target, field2, x1, y1, x2, y2 )
  use share_vars
  implicit none
  integer :: i,j,nx1,ny1
  real (kind=pr) :: dx, dy
  real (kind=pr), intent (in) :: field2(0:,0:), x_target, y_target, x1, y1, x2, y2
  real (kind=pr), allocatable :: x(:), y2a(:,:), x1a(:,:), y(:)

  dx = (x2-x1)/real(size(field2,1)-1)
  dy = (y2-y1)/real(size(field2,2)-1)

  nx1=size(field2,1)
  ny1=size(field2,2)

  if ( (x_target > x2).or.(x_target < x1).or.(y_target > y2).or.(y_target < y1) ) then
    write(*,'("target: (",es11.4,"|",es11.4,") but box: (",es11.4,"|",es11.4,") x (",es11.4,"|",es11.4,")")') &
    x_target, y_target, x1, y1, x2, y2
    write (*,*) "!!! SplineInterpolation: target coordinates not in the field."
    stop
  endif

  allocate ( y(0:ny1-1),x(0:nx1-1), y2a(0:nx1-1,0:ny1-1), x1a(0:nx1-1,0:ny1-1) )

  do i=0,nx1-1
    x(i)=real(i)*dx + x1 ! the grid
  enddo
  do i=0,ny1-1
    y(i)=real(i)*dy + y1 ! the grid
  enddo
  
  call splie2(x,y,field2,y2a)

 
  SplineInterpolation=splin2(x,y,field2,y2a,x_target,y_target)

  return

end function SplineInterpolation

!=================================================================================================================================

	subroutine bcucof(y,y1,y2,y12,d1,d2,c)
	use nrtype
	implicit none
	real(sp), intent(in) :: d1,d2
	real(sp), dimension(4), intent(in) :: y,y1,y2,y12
	real(sp), dimension(4,4), intent(out) :: c
	real(sp), dimension(16) :: x
	real(sp), dimension(16,16) :: wt
	data wt /1,0,-3,2,4*0,-3,0,9,-6,2,0,-6,4,&
		8*0,3,0,-9,6,-2,0,6,-4,10*0,9,-6,2*0,-6,4,2*0,3,-2,6*0,-9,6,&
		2*0,6,-4,4*0,1,0,-3,2,-2,0,6,-4,1,0,-3,2,8*0,-1,0,3,-2,1,0,-3,&
		2,10*0,-3,2,2*0,3,-2,6*0,3,-2,2*0,-6,4,2*0,3,-2,0,1,-2,1,5*0,&
		-3,6,-3,0,2,-4,2,9*0,3,-6,3,0,-2,4,-2,10*0,-3,3,2*0,2,-2,2*0,&
		-1,1,6*0,3,-3,2*0,-2,2,5*0,1,-2,1,0,-2,4,-2,0,1,-2,1,9*0,-1,2,&
		-1,0,1,-2,1,10*0,1,-1,2*0,-1,1,6*0,-1,1,2*0,2,-2,2*0,-1,1/
	x(1:4)=y
	x(5:8)=y1*d1
	x(9:12)=y2*d2
	x(13:16)=y12*d1*d2
	x=matmul(wt,x)
	c=reshape(x,(/4,4/),order=(/2,1/))
	end subroutine bcucof

!=================================================================================================================================
	subroutine bcuint(y,y1,y2,y12,x1l,x1u,x2l,x2u,x1,x2,ansy,ansy1,ansy2)
! 	use nrtype; use nrutil, only : nrerror
use nrtype 
use nrutil
	implicit none
	real(sp), dimension(4), intent(in) :: y,y1,y2,y12
	real(sp), intent(in) :: x1l,x1u,x2l,x2u,x1,x2
	real(sp), intent(out) :: ansy,ansy1,ansy2
	integer(i4b) :: i
	real(sp) :: t,u
	real(sp), dimension(4,4) :: c
	call bcucof(y,y1,y2,y12,x1u-x1l,x2u-x2l,c)
	if (x1u == x1l .or. x2u == x2l) call &
		nrerror('bcuint: problem with input values - boundary pair equal?')
	t=(x1-x1l)/(x1u-x1l)
	u=(x2-x2l)/(x2u-x2l)
	ansy=0.0
	ansy2=0.0
	ansy1=0.0
	do i=4,1,-1
		ansy=t*ansy+((c(i,4)*u+c(i,3))*u+c(i,2))*u+c(i,1)
		ansy2=t*ansy2+(3.0_sp*c(i,4)*u+2.0_sp*c(i,3))*u+c(i,2)
		ansy1=u*ansy1+(3.0_sp*c(4,i)*t+2.0_sp*c(3,i))*t+c(2,i)
	end do
	ansy1=ansy1/(x1u-x1l)
	ansy2=ansy2/(x2u-x2l)
	end subroutine bcuint


!=================================================================================================================================
	function locate(xx,x)
	use nrtype
	implicit none
	real(sp), dimension(:), intent(in) :: xx
	real(sp), intent(in) :: x
	integer(i4b) :: locate
	integer(i4b) :: n,jl,jm,ju
	logical :: ascnd
	n=size(xx)
	ascnd = (xx(n) >= xx(1))
	jl=0
	ju=n+1
	do
		if (ju-jl <= 1) exit
		jm=(ju+jl)/2
		if (ascnd .eqv. (x >= xx(jm))) then
			jl=jm
		else
			ju=jm
		end if
	end do
	if (x == xx(1)) then
		locate=1
	else if (x == xx(n)) then
		locate=n-1
	else
		locate=jl
	end if
	end function locate
  
  !=================================================================================================================================
  
  	subroutine tridag_ser(a,b,c,r,u)
	use nrtype; use nrutil, only : assert_eq,nrerror
	implicit none
	real(sp), dimension(:), intent(in) :: a,b,c,r
	real(sp), dimension(:), intent(out) :: u
	real(sp), dimension(size(b)) :: gam
	integer(i4b) :: n,j
	real(sp) :: bet
	n=assert_eq((/size(a)+1,size(b),size(c)+1,size(r),size(u)/),'tridag_ser')
	bet=b(1)
	if (bet == 0.0) call nrerror('tridag_ser: error at code stage 1')
	u(1)=r(1)/bet
	do j=2,n
		gam(j)=c(j-1)/bet
		bet=b(j)-a(j-1)*gam(j)
		if (bet == 0.0) &
			call nrerror('tridag_ser: error at code stage 2')
		u(j)=(r(j)-a(j-1)*u(j-1))/bet
	end do
	do j=n-1,1,-1
		u(j)=u(j)-gam(j+1)*u(j+1)
	end do
	end subroutine tridag_ser
  
  !=================================================================================================================================

	recursive subroutine tridag_par(a,b,c,r,u)
	use nrtype; use nrutil, only : assert_eq,nrerror
!	use nr, only : tridag_ser
	implicit none
	real(sp), dimension(:), intent(in) :: a,b,c,r
	real(sp), dimension(:), intent(out) :: u
	integer(i4b), parameter :: npar_tridag=4
	integer(i4b) :: n,n2,nm,nx1
	real(sp), dimension(size(b)/2) :: y,q,piva
	real(sp), dimension(size(b)/2-1) :: x,z
	real(sp), dimension(size(a)/2) :: pivc
	n=assert_eq((/size(a)+1,size(b),size(c)+1,size(r),size(u)/),'tridag_par')
	if (n < npar_tridag) then
		call tridag_ser(a,b,c,r,u)
	else
		if (maxval(abs(b(1:n))) == 0.0) &
			call nrerror('tridag_par: possible singular matrix')
		n2=size(y)
		nm=size(pivc)
		nx1=size(x)
		piva = a(1:n-1:2)/b(1:n-1:2)
		pivc = c(2:n-1:2)/b(3:n:2)
		y(1:nm) = b(2:n-1:2)-piva(1:nm)*c(1:n-2:2)-pivc*a(2:n-1:2)
		q(1:nm) = r(2:n-1:2)-piva(1:nm)*r(1:n-2:2)-pivc*r(3:n:2)
		if (nm < n2) then
			y(n2) = b(n)-piva(n2)*c(n-1)
			q(n2) = r(n)-piva(n2)*r(n-1)
		end if
		x = -piva(2:n2)*a(2:n-2:2)
		z = -pivc(1:nx1)*c(3:n-1:2)
		call tridag_par(x,y,z,q,u(2:n:2))
		u(1) = (r(1)-c(1)*u(2))/b(1)
		u(3:n-1:2) = (r(3:n-1:2)-a(2:n-2:2)*u(2:n-2:2) &
			-c(3:n-1:2)*u(4:n:2))/b(3:n-1:2)
		if (nm == n2) u(n)=(r(n)-a(n-1)*u(n-1))/b(n)
	end if
	end subroutine tridag_par
  
  !=================================================================================================================================
  
  function splint(xa,ya,y2a,x)
	use nrtype; use nrutil, only : assert_eq,nrerror
!	use nr, only: locate
	implicit none
	real(sp), dimension(:), intent(in) :: xa,ya,y2a
	real(sp), intent(in) :: x
	real(sp) :: splint
	integer(i4b) :: khi,klo,n
	real(sp) :: a,b,h
	n=assert_eq(size(xa),size(ya),size(y2a),'splint')
	klo=max(min(locate(xa,x),n-1),1)
	khi=klo+1
	h=xa(khi)-xa(klo)
	if (h == 0.0) call nrerror('bad xa input in splint')
	a=(xa(khi)-x)/h
	b=(x-xa(klo))/h
	splint=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.0_sp
	end function splint

!=================================================================================================================================

	subroutine spline(x,y,yp1,ypn,y2)
	use nrtype; use nrutil, only : assert_eq
!	use nr, only : tridag
	implicit none
	real(sp), dimension(:), intent(in) :: x,y
	real(sp), intent(in) :: yp1,ypn
	real(sp), dimension(:), intent(out) :: y2
	integer(i4b) :: n
	real(sp), dimension(size(x)) :: a,b,c,r
	n=assert_eq(size(x),size(y),size(y2),'spline')
	c(1:n-1)=x(2:n)-x(1:n-1)
	r(1:n-1)=6.0_sp*((y(2:n)-y(1:n-1))/c(1:n-1))
	r(2:n-1)=r(2:n-1)-r(1:n-2)
	a(2:n-1)=c(1:n-2)
	b(2:n-1)=2.0_sp*(c(2:n-1)+a(2:n-1))
	b(1)=1.0
	b(n)=1.0
	if (yp1 > 0.99e30_sp) then
		r(1)=0.0
		c(1)=0.0
	else
		r(1)=(3.0_sp/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
		c(1)=0.5
	end if
	if (ypn > 0.99e30_sp) then
		r(n)=0.0
		a(n)=0.0
	else
		r(n)=(-3.0_sp/(x(n)-x(n-1)))*((y(n)-y(n-1))/(x(n)-x(n-1))-ypn)
		a(n)=0.5
	end if
	call tridag_par(a(2:n),b(1:n),c(1:n-1),r(1:n),y2(1:n))
	end subroutine spline
  
  !=================================================================================================================================
  
  	function splin2(x1a,x2a,ya,y2a,x1,x2)
	use nrtype; use nrutil, only : assert_eq
!	use nr, only : spline,splint
	implicit none
	real(sp), dimension(:), intent(in) :: x1a,x2a
	real(sp), dimension(:,:), intent(in) :: ya,y2a
	real(sp), intent(in) :: x1,x2
	real(sp) :: splin2
	integer(i4b) :: j,m,ndum
	real(sp), dimension(size(x1a)) :: yytmp,y2tmp2
	m=assert_eq(size(x1a),size(ya,1),size(y2a,1),'splin2: m')
	ndum=assert_eq(size(x2a),size(ya,2),size(y2a,2),'splin2: ndum')
	do j=1,m
		yytmp(j)=splint(x2a,ya(j,:),y2a(j,:),x2)
	end do
	call spline(x1a,yytmp,1.0e30_sp,1.0e30_sp,y2tmp2)
	splin2=splint(x1a,yytmp,y2tmp2,x1)
	end function splin2
  
  !=================================================================================================================================

	subroutine splie2(x1a,x2a,ya,y2a)
	use nrtype; use nrutil, only : assert_eq
!	use nr, only : spline
	implicit none
	real(sp), dimension(:), intent(in) :: x1a,x2a
	real(sp), dimension(:,:), intent(in) :: ya
	real(sp), dimension(:,:), intent(out) :: y2a
	integer(i4b) :: j,m,ndum
	m=assert_eq(size(x1a),size(ya,1),size(y2a,1),'splie2: m')
	ndum=assert_eq(size(x2a),size(ya,2),size(y2a,2),'splie2: ndum')
	do j=1,m
		call spline(x2a,ya(j,:),1.0e30_sp,1.0e30_sp,y2a(j,:))
	end do
	end subroutine splie2




end module interpolation

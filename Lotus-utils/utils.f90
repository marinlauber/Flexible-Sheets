module uMod
  use fieldMod, only: field
  use vectorMod, only: vfield
  use fluidMod, only: fluid
!
! -- average field type
  type :: avrgField
    type(vfield) :: velocity
    type(field)  :: pressure
    integer      :: iter=0
   contains
    procedure :: init, add, write
  end type avrgField

  public :: avrgField,profile

  contains
!
! -- initializes an average field
  subroutine init(a,flow)
    implicit none
    class(avrgField),intent(inout) :: a
    type(fluid),intent(inout)     ::  flow
    call a%velocity%init(flow%velocity%e(1)%size())
    call a%pressure%init(flow%pressure%size())
    call a%add(flow)
  end subroutine init
!
! -- add realization
  subroutine add(a,flow)
    implicit none
    class(avrgField),intent(inout) :: a
    type(fluid),intent(inout)     ::  flow
    integer :: d
    ! average values
    do d=1,a%velocity%ndims
      a%velocity%e(d)%p = a%velocity%e(d)%p+flow%velocity%e(d)%p
    end do
    a%pressure%p = a%pressure%p+flow%pressure%p
    a%iter = a%iter+1
  end subroutine add
!
! -- write average fields
  subroutine write(a)
    use ioMod, only: write_vtk
    implicit none
    class(avrgField),intent(inout) :: a
    integer :: d
    ! average values
    do d=1,a%velocity%ndims
      a%velocity%e(d)%p = a%velocity%e(d)%p/a%iter
    end do
    a%pressure%p = a%pressure%p/a%iter
    call write_vtk(a%velocity,a%pressure,name_in='avrgF.')
  end subroutine write
!
! -- write profile along a given segment
  subroutine profile(sfield,x0,x1,n,fn)
    use fieldMod, only: field
    implicit none
    type(field),intent(in) :: sfield
    real,intent(in) :: x0(3),x1(3)
    integer,intent(in) :: n
    integer,intent(in) :: fn ! fort.fn file to write to
    real :: p(n),dx(3),x(3)
    integer :: i
    dx = (x1-x0)/(n-1)
    do i=1,n
    x = x0+dx*(i-1)
    p(i) = sfield%at(x)
    end do
    write(fn,'(*(e16.8))') p(:)
    flush(fn)
  end subroutine profile
end module uMod
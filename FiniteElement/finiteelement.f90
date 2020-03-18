module FiniteElementMod

  type,public :: FEAsolver
    real :: zeta=1,omega=1
    real :: xi(3),dt,dxi(3)
    real,allocatable :: IEN(:,:) ! element node array, how the elements are linked to the nodes.
    real,allocatable :: ID(:)    ! destination array, links the global node number to the global equation number in the final linear system
  contains
    procedure,public :: init,update
  end type

 contains

  subroutine init(a,xi0,zeta,omega)
    class(FEAsolver),intent(inout) :: a
    real,intent(in) :: xi0,zeta,omega
    a%zeta = zeta; a%omega = omega
    a%xi = (/xi0,xi0,xi0/)
  end subroutine init

  subroutine update(a,q,dt)
    class(FEAsolver),intent(inout) :: a
    real,intent(in) :: q,dt    
    real :: xi,xi_1,xi_2
    xi=a%xi(1); xi_1=a%xi(2); xi_2=a%xi(3)
! -- update in time using second order forward difference
    xi_n = (q + (5+4*a%zeta*a%omega*dt)*xi-(4+a%zeta*a%omega*dt)*xi_1+xi_2) &
            /(2+3*a%zeta*a%omega*dt+a%omega**2*dt**2)
! -- prepare next time step
    a%xi = (/xi_n, xi, xi_1/)
  end subroutine update
!
! -- test on a oscillator, this writes to a file
  subroutine FEAsolver_test
    implicit none
    type(FEAsolver) :: FEsolver
    real,allocatable :: vals(:),ts(:)
    real :: dt,time=0,stop=10
    integer :: i = 1
    call FEsolver%init(xi0=1.,zeta=sqrt(10.)/40,omega=sqrt(10.))
    dt = 0.001
    allocate(vals(int(stop/dt)),ts(int(stop/dt)))
    do while(time<stop)
      call FEsolver%update(q=0.,dt=dt)
      time = time+dt
      vals(i) = FEsolver%xi(1)
      ts(i) = time
      i=i+1
    end do
    open(14,file='out.csv')
    write(14,*) vals
    write(14,*) ts
    close(14)
  end subroutine FEAsolver_test

end module FiniteElementMod
program test
  use FiniteElementMod
  call FEAsolver_test
end program test

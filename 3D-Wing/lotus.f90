program ellipsoidalwing
  use fluidMod,   only: fluid
  use fieldMod,   only: field
  use vectorMod,  only: vfield
  use bodyMod,    only: body
  use mympiMod
  use gridMod
  use imageMod,   only: display
  use geom_shape

  implicit none
!
! -- numerical parameters
  real,parameter    :: C=32               ! resolution (pnts per chord)
  integer,parameter :: ndims=3
  real,parameter    :: dr=2.
  integer           :: n(3)               ! number of points
  logical           :: p(3) = (/.FALSE.,.FALSE.,.FALSE./) 

! -- motion params
  real,parameter  :: A_theta=0., A_phi=0.35*pi, A_alpha=pi/4.
  real,parameter  :: fr=.5, f_theta=.5, xi=0.
  real,parameter  :: T = 1./fr*C
  real,parameter  :: Rg = 1.581
  real,parameter  :: c_bar = pi/4.
  real,parameter  :: visc = 0.0273
  real,parameter  :: Re=4*A_phi*fr*Rg*c_bar/visc
  real            :: nu=C/Re

! -- MPI utils
  integer :: b(3) = (/1,1,1/)                                   ! 1 blocks
  integer :: m(3)                                               ! number of points per block
  logical :: root                                               ! root processor
  logical :: there=.false.                                      ! exits the time loop
! -- utils                      
  real,parameter :: Np = 1.25                                    ! number of period
  real  :: finish =  2*Np*T                                      ! number of timesteps
  real,parameter :: dtPrint = 0.1*T                              ! print rate
! -- variables
  type(fluid) :: flow
  type(body)  :: geom
  real :: dt,t1,pforce(ndims),vforce(ndims),Area

! -- Initialize MPI (if MPI is OFF, b is set to 1)
  call init_mympi(ndims,set_blocks=b(:ndims),set_periodic=p(:ndims))
  root = mympi_rank()==0
  eps2 = 0.5
!
! -- Initialize the 3D ellipse
  n = composite(C*(/4.,3.,3./), prnt=root)
  call xg(1)%stretch(n(1), -3.*C,-1.5*C,  1.5*C,  3*C,h_min=1.,prnt=root)
  call xg(2)%stretch(n(2),-3.5*C,-1.5*C,  0.5*C,2.5*C,h_min=1.,prnt=root)
  call xg(3)%stretch(n(3),-2.25*C,-0.25*C,1.75*C,3.75*C,h_min=1.,prnt=root)
  ! geom = sphere(radius=0.5*C,center=0)&
  !         .map.init_scale(1,lengthx,rate)&
  !         .map.init_scale(2,lengthy,rate)&
  !         .map.init_rigid(3,z)&
  !         .map.init_rigid(6,alpha)&
  !         .map.init_rigid(5,phi)&
  !         .map.init_rigid(4,theta)
  geom =  cylinder(axis=2,radius=0.5*C,center=0).and. &
          plane(norm=(/0,-1,0/),center=(/0.,-dr/2.,0./)).and. &
          plane(norm=(/0, 1,0/),center=(/0., dr/2.,0./))&
          .map.init_scale(1,lengthx,rate)&
          .map.init_rigid(3,z)&
          .map.init_rigid(6,alpha)&
          .map.init_rigid(5,phi)&
          .map.init_rigid(4,theta)
          
  if(root) then
    print*, 'Simulation parameters:'
    print*, 'DOF:', product(n)
    print*, 'finish', finish
    print*, 'A_alpha',A_alpha
    print*, 'A_theta',A_theta
    print*, 'A_phi',A_phi
    print*, 'fr',fr
    print*, 'f_theta',f_theta
    print*, 'Period',T
    print*, 'Reynolds number',Re
  end if

  Area = pi/2.*C
!
! -- Initialize fluid
  m = n/b
  call flow%init(m, geom, V=[0.,0.,0.], nu=nu)
  call flow%write(geom,lambda=.true.)

  if(root) print *, '-- init complete --'
! -- Time update loop
  finish = finish+flow%time
  flow%dt = 1.
  do while (flow%time<finish .and. .not.there)
    dt = flow%dt                    ! time step
    t1 = flow%time+dt               ! time at the end of this step
    call geom%update(flow%time)
    call flow%update(geom)          ! update N-S
    flow%dt = min(1., flow%dt)
    pforce = -2.*geom%pforce(flow%pressure)/Area
    vforce =  2.*nu*geom%vforce(flow%velocity)/Area
    if(root) write(9,'(f10.4,f8.4,6e16.8,3f10.4)') t1/T,dt,pforce,vforce,&
                   phi(real(flow%time,8)),theta(real(flow%time,8)),alpha(real(flow%time,8))
    if(root) flush(9)
    if(mod(t1,dtPrint)<dt) then
      if(root) print "('Time:',f15.3,'. Time remaining:',f15.3)",t1/C,finish/C-t1/C
      call display(flow%velocity%vorticity_Z(), 'out_vort', lim = 0.5)
      call flow%write(geom,lambda=.true.)
    end if
!
! -- check for kill flag
    inquire(file='.kill', exist=there)
  end do

  call flow%write(geom,lambda=.true.)
  call mympi_end()

  contains
  real(8) pure function lengthx(time)
    real(8), intent(in) :: time
    lengthx = .5
  end function lengthx
  real(8) pure function lengthy(time)
    real(8), intent(in) :: time
    lengthy = .05
  end function lengthy
  real(8) pure function z(time)
    real(8),intent(in) :: time
    z = 0.75*C
  end function z
  real(8) pure function rate(time)
    real(8), intent(in) :: time
    rate = 0.0
  end function rate
  real(8) pure function theta(time)
    real(8),intent(in) :: time
    theta = A_theta*cos(2*pi*f_theta*time/T)
  end function theta
  real(8) pure function phi(time)
    real(8),intent(in) :: time
    phi = A_phi*sin(2*pi*fr*time/T)
  end function phi
  real(8) pure function alpha(time)
    real(8),intent(in) :: time
    alpha = pi/2. - A_alpha*cos(2*pi*fr*time/T + xi)
  end function alpha
end program ellipsoidalwing

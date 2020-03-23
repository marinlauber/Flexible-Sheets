program sphere_stl
  use flexMod,    only: flexBody
  use fluidMod,   only: fluid
  use mympiMod,   only: init_mympi,mympi_end,mympi_rank
  use gridMod,    only: xg
  use imageMod,   only: display
  use geom_shape
  use geom_global,only: pi
  
  implicit none
!
! -- Define parameter, declare variables
  real,parameter :: D = 2*32, Re = 30e3 ! diameter and Reynolds number
  integer :: n(3) = 32*(/4,2,2/)  ! numer of points
  integer :: b(3) = (/1,1,1/)     ! MPI blocks (product must equal n_procs)
  logical :: root                 ! root processor
  type(fluid) :: flow
  type(flexBody) :: geom
!
! -- Initialize MPI (if MPI is OFF, b is set to 1)
  call init_mympi(3,set_blocks=b)
  root = mympi_rank()==0
!
  if(root) print *,'       Setting up the grid, body and fluid      '
  if(root) print *,'------------------------------------------------'
  call xg(1)%stretch(n(1), -2*D, -0.5*D, 0.5*D, 5*D, h_min=2., h_max=12.,prnt=root)
  call xg(2)%stretch(n(2), -2*D, -0.5*D, 0.5*D, 2*D, h_min=2.,prnt=root)
  call xg(3)%stretch(n(3), -2*D, -0.5*D, 0.5*D, 2*D, h_min=2.,prnt=root)

! build geometry
  call geom%init(fname='../sphere.stl',s=D*(/.005,.005,.005/),xi0=(/0.,-16.,0./))
  call flow%init(n/b,geom,V=(/1.,0.,0./),nu=D/Re)
!
  if(root) print *,'            Starting time update loop           '
  if(root) print *,'------------------------------------------------'
  if(root) print *,' -t- , -dt- '

  do while(flow%time<10*D)
    call geom%updateFlex(flow)   ! update position
    call flow%update(geom)
    write(9,'(f10.4,f8.4,3e16.8)') flow%time/D,flow%dt, 2.*geom%pforce(flow%pressure)/(pi*D**2/4.)
    if(mod(flow%time,0.1*D)<flow%dt) then
      if(root) print '(f6.1,",",f6.3)',flow%time/D,flow%dt
      call geom%writePoints(flow%pressure,flow%time)
      call flow%write(geom)
      call display(flow%velocity%vorticity_Z(),'out_vort',lim=20./D)
    end if
  end do
  
  call geom%writeSTL('moved')

  if(root) print *,'Loop complete: writing restart files and exiting'
  if(root) print *,'------------------------------------------------'
  call mympi_end

end program sphere_stl

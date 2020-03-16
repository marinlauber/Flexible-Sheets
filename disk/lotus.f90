program Disk
    use fluidMod,   only: fluid
    use fieldMod,   only: field
    use bodyMod,    only: body
    use geom_global,only: eps
    use mympiMod
    use gridMod
    use imageMod,   only: display
    use geom_shape
  
    implicit none
  !
  ! -- numerical parameters
    real,parameter  :: C=32           ! resolution (pnts per chord)
    real,parameter  :: Re=100
    real,parameter  :: U=1
    real            :: dr=1+sqrt(3.)
    integer         :: n(3)            ! number of points
  ! -- MPI utils
    integer :: b(3) = (/4,2,1/)        ! 32 blocks
    integer :: d                       ! number of points per block
    logical :: root                    ! root processor
    logical :: kill=.false.            ! exits the time loop
  ! -- variables
    type(fluid) :: flow
    type(body)  :: geom
    real :: pforce(3), vforce(3), h(3)
  !
  ! -- change kernel width
    eps=0.5

  ! -- Initialize MPI (if MPI is OFF, b is set to 1)
    call init_mympi(3,set_blocks=b)
    root = mympi_rank()==0
  !
    geom = disc(C)

  ! -- Initialize grid (uniform gird)
    if(root) print *,'Setting up the grid, body and fluid'
    if(root) print *,'-----------------------------------'
    n = composite(C*(/4,4,4/), prnt=root)
    call xg(1)%stretch(n(1), -6*C, -2*C, 2*C, 6*C, h_min=2., h_max=12.,prnt=root)
    call xg(2)%stretch(n(2), -6*C, -2*C, 2*C, 6*C, h_min=2.,prnt=root)
    call xg(3)%stretch(n(3), -6*C, -2*C, 2*C, 6*C, h_min=2.,prnt=root)
  
  ! -- Initialize fluid
    call flow%init(n/b,geom,V=(/0.,0.,0./),nu=C/Re)
    call flow%write(geom)
  
    if(root) print *,'Starting time update loop'
    if(root) print *,'-----------------------------------'
    if(root) print *,' -t- , -dt- '
  
    do while (flow%time<2*U*C)
      call geom%update(flow%time+flow%dt)
      call flow%update(geom) ! update N-S
      pforce = -2.*geom%pforce(flow%pressure)/C
      vforce = 2.*Re*geom%vforce(flow%velocity)
      if(root) write(9,'(f10.4,f8.4,6e16.8)') flow%time/C,flow%dt,pforce,vforce
      if(mod(flow%time,0.1*C)<flow%dt) then
        if(root) print '(f6.1,",",f6.3)',flow%time/C,flow%dt
        call display(flow%velocity%vorticity_Z(), 'out_vort', lim = 20./C)
        call flow%write(geom)
      end if
    end do

    if(root) print *,'Loop complete: writing restart files and exiting'
    if(root) print *,'-----------------------------------'
    call mympi_end()
  
    contains

    real(8) pure function y(t)
      real(8),intent(in) :: t
      y = C-U*t
    end function y
    type(body) function disc(D)
      real,intent(in) :: D
        disc = (cylinder(axis=1,radius=D,center=(/0,0,0/)) &
               -plane(norm=(/-1,0,0/), center=dr/2.) - plane(norm=(/1,0,0/), center=-dr/2.)).map.init_rigid(1,y)
    end function disc
  end program Disk
  
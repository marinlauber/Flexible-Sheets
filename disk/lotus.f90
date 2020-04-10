!-------------------------------------------------------!
!------------------- Impulsive Disk --------------------!
!-------------------------------------------------------!
program disk_impulse
  use fieldMod,   only: field
  use fluidMod,   only: fluid
  use bodyMod,    only: body
  use mympiMod,   only: init_mympi,mympi_end,mympi_rank
  use gridMod,    only: xg,composite
  use imageMod,   only: display
  use geom_shape  ! to define geom (set,eps,plane, etc)
  implicit none
  real,parameter     :: R = 128            ! length scale
  real,parameter     :: Re = 1.25e5        ! Reynolds number
  real,parameter     :: a0 = 0.5           ! acceleration
  real,parameter     :: nu = 2*R/Re        ! viscosity to achieve Re
  integer            :: b(3) = (/1,1,1/)   ! blocks
  integer            :: n(3) = (/2,1,1/)*R ! number of points
  logical            :: root,there=.false. ! flags
  real               :: dt,t,force(3)
  real               :: s=0,U=0,U0,a,dr
!
  type(fluid)        :: flow
  type(body)         :: disk
  type(field)        :: f
  real,pointer       :: p(:,:,:)
!
! -- reduce kernel
  eps2=0.5
  dr=1+sqrt(2.)
!
! -- Initialize geometry
  call init_mympi(3,set_blocks=b)
  root = mympi_rank()==0
  call xg(1)%stretch(n(1),-10*R,-0.25*R,R,10*R,prnt=root)
  call xg(2)%stretch(n(2),0.,0.,1.25*R,6.66*R,h_min=2.,prnt=root)
  call xg(3)%stretch(n(3),0.,0.,1.25*R,6.66*R,h_min=2.,prnt=root)
  disk = cylinder(axis=1,radius=R,center=0.).and. &
         plane(norm=(/-1,0,0/),center=(/-dr/2.,0.,0./)).and. &
         plane(norm=(/ 1,0,0/),center=(/ dr/2.,0.,0./))
!
! -- Initialize fluid
  call flow%init(n/b,disk,nu=nu,exit=.false.,V=(/U,0.,0./))
  if(root) print *, '-- init complete --'
  call disk%writeFields(disk%kernel)
!
! -- Time update loop
  do while (s<3.01 .and..not.there)
    flow%dt = min(1.,flow%dt)
    dt = flow%dt
    t = (flow%time+dt)/R
    U0 = U
    U = a0*t+(1+tanh(31.4*(t-1./a0)))/2.*(1-a0*t)
    s = s+0.5*(U+U0)*dt/R; a = (U-U0)/dt*R
    call flow%update(V=(/U,0.,0./))
!
! -- measure & print force
    force = -4.0*disk%pforce(flow%pressure)/R**2
    write(9,1) t,dt,force(1),s,U,a
1   format(f10.4,f8.4,e16.8,3f10.4)
    flush(9)
    call vort_props(flow%velocity)
!
! -- print and write image
    if(mod(s,0.1)<U*dt/R) then
      f = flow%velocity%vorticity_Z()
      p => f%point()
      call display(p(:,:,1:1),'out_vort',box=int((/-R,0.,3*R,2*R/)),lim=0.2)
      if(root) print 1,t,dt,force(1),s,U,a
    end if
!
! -- check for kill flag
    inquire(file='.kill', exist=there)
  end do

! -- write to file and finalize
  call flow%write(lambda=.true.)
  if(root) write(6,*) '--- complete --- '
  call mympi_end
contains
  subroutine vort_props(v)
    use gridMod, only: dv
    use geom_global, only: pi
    use mympiMod, only: mympi_sum
    use vectorMod, only: vfield,vorticity
    type(vfield),intent(in) :: v
    type(vfield) :: S(3),omega
    real :: x(3),Rhat(3),Rmag,omega_a,omega_z,RdRdz
    real :: Pp,Pn,Lp,Ln,Gp,Gn,Xp,Xn
    integer :: i,j,k,is,ie,js,je,ks,ke
    Pp=0;Pn=0;Lp=0;Ln=0;Gp=0;Gn=0;Xp=0;Xn=0
!
! -- get vorticity and loop over array
    call v%gradTensor(S)
    call vorticity(S,omega)
    call omega%e(1)%limits(is,ie,js,je,ks,ke)
    do i=is,ie; do j=js,je; do k=ks,ke
!
! -- get cylindrical coords and vector components
      x = omega%e(1)%pos(i,j,k)
      Rhat = (/0.,x(2),x(3)/); Rmag = sqrt(sum(Rhat**2))
      Rhat = Rhat/Rmag; Rmag = Rmag/R
      omega_a = -omega%e(2)%p(i,j,k)*Rhat(3)+omega%e(3)%p(i,j,k)*Rhat(2)
      omega_z = omega%e(1)%p(i,j,k)
      RdRdz = 4.0*dv(i-2,j-2,k-2)/R**2
!
! -- get integrals for impulse and circulation
      if(x(1)>4) then
        Xp = Xp+omega_a*Rmag*RdRdz*x(1)/R
        Pp = Pp+omega_a*Rmag*RdRdz
        Lp = Lp+omega_z*Rmag**2*RdRdz
        Gp = Gp+omega_a/Rmag*RdRdz
      else
        Xn = Xn+omega_a*Rmag*RdRdz*x(1)/R
        Pn = Pn+omega_a*Rmag*RdRdz
        Ln = Ln+omega_z*Rmag**2*RdRdz
        Gn = Gn+omega_a/Rmag*RdRdz
      end if
    end do; end do; end do
!
! -- MPI sum
    call mympi_sum(Xp)
    call mympi_sum(Pp)
    call mympi_sum(Lp)
    call mympi_sum(Gp)
    call mympi_sum(Xn)
    call mympi_sum(Pn)
    call mympi_sum(Ln)
    call mympi_sum(Gn)
!
! -- write
    write(14,2) t,0.5*Pp,-0.5*Lp,Gp/(2.*pi),sqrt(Pp/Gp),Xp/Pp, &
                  0.5*Pn,-0.5*Ln,Gn/(2.*pi),sqrt(Pn/Gn),Xn/Pn
2   format(f10.4,10e16.8)
    flush(14)
  end subroutine vort_props
end program disk_impulse

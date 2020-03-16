# Lotus Flexible Module

#### (See Lotus/solver/oop/flexible.f90)

This module extends the standard `type(body)` for flexible bodies. This is also the module you want to use to probe data points on the surface of an `.stl`  geometry. To use this module you need to import it instead of `bodyMod` using

```fortran
use flexMod, only: flexBody
```
and then create a geometry based on this new body type

```fortran
type(flexBody) :: geom
```

## The standard way

The main difference between this module and the standard `body` is the wat you initialize a geometry based on an `.stl` file. The standard way of doing this is to have a function of with the following structure
```fortran
type(set) function sphere(R)
    use geom_shape
    real,intent(in) :: R
    type(model_info) :: mod_info

    mod_info%file = '../sphere.stl'     ! file to use
    mod_info%s = R*(/1.,1.,1./)         ! scale it
    mod_info%x = (/0.,0.,0./)           ! define origin
    mod_info%r = (/0.,0.,0./)           ! rotation

    sphere = model_init(mod_info)

end function sphere
```
This function return a `type(set)` that is then assigned to the geometry onto which we can add a mapping
```fortran
type(body) :: geom

geom = sphere(0.5*D).map.init_rigid(2,y)
```

## The other way
The way you have to initialize the geometry if you want to be able to probe nodal values of the flow on the body is as follows.
```fortran
subroutine sphere(geom)
    use geom_shape
    class(flexBody),intent(inout) :: geom
    class(model),pointer :: sphere_model
    type(stlth),pointer :: ptr
    type(set) :: sets
    type(model_info) :: mod_info

    mmod_info%file = '../sphere.stl'     ! file to use
    mod_info%s = R*(/1.,1.,1./)          ! scale it
    mod_info%x = (/0.,0.,0./)            ! define origin
    mod_info%r = (/0.,0.,0./)            ! rotation

    call model_init_ptr(mod_info, sphere_model, ptr)
    sets = sphere_model
    geom = sets.map.init_rigid(2,y)      ! apply mapping
    geom%srf => ptr                      ! store nodes

end subroutine sphere
```
We now have a __subroutine__ to initialize the geometry for reasons that will not be detailed. The geom is initialzed in the main part of the code as
```fortran
type(flexBody) :: geom

call sphere(geom)
```

The rest of the fluid `lotus.f90` file is the same, except that now you can output nodal values of, say the pressure, every N time step
```fortran
do while(flow%time<finish)

    call geom%update(flow%time+flow%dt)   ! update position
    call flow%update(geom)

    write(9,'(f10.4,f8.4,3e16.8)') flow%time/D,flow%dt, 2.*geom%pforce(flow%pressure)/(pi*D**2/4.)
    
    if(mod(flow%time,wtime)<flow%dt) then
      if(root) print '(f6.1,",",f6.3)',flow%time/D,flow%dt

      call geom%writePoints(flow%pressure,flow%time)        ! write pressure point to "surf.vtp"
      call flow%write(geom)
      call display(flow%velocity%vorticity_Z(),'out_vort',lim=20./D)
    end if

end do
```
The main loop using the standard body type would be exactly the same, without the `call geom%writePoints(...)` call.

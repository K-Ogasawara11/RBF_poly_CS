module md_params
#include <petsc/finclude/petscksp.h>
    use petscksp
    use petscisdef
    implicit none 
    
    !node
    character(len=20) :: node='md'
    !character(len=20) :: node='bauer'
    !character(len=20) :: node='icos_nomodify'

    integer,parameter :: n=10000 ! number of node
    integer,parameter ::  ss=55 !Stencil Size  
    integer,parameter :: ie= n*ss
    
    !fillter 
    double precision,parameter :: kk=4d0 ! oder of hyperviscosity
    double precision,parameter :: gg=0d0 ! constants of hyperviscosity
    double precision,parameter :: gamma =-gg/(dble(n)**kk) !
    
    !time 
    double precision,parameter :: dt =600d0 !time step
    double precision,parameter :: day=5d0 ! calc day
    double precision,parameter :: nd=day * (60d0*60d0*24d0) !  calc day in sec 
    integer,parameter ::  step_in_hour=int(3600/dt) ! steps in a hour
    integer,parameter ::  steps=int(nd/dt) ! number of steps
    !integer,parameter ::  steps=step_in_hour*6 ! number of steps by hour
    integer,parameter ::  nt=steps ! number of steps
    integer,parameter ::  rec=step_in_hour*24! 
    integer,parameter ::  output_time=step_in_hour*1!


    ! test case
    character(len=20) :: notc='2' !nonlinear zonal geostropic flow
    !character(len=20) :: notc='6' !Rossby-Haurwitz wave
    double precision,parameter ::pi=dble(dacos(-1d0))
    double precision,parameter ::alpha=pi/2d0 !angle of rotation

    !constants
    double precision,parameter ::g=9.80616d0 !grav
    double precision,parameter ::ea=6.37122d6 !earth radius
    double precision,parameter ::omega=7.292d-5 !earth angular velocity
    double precision,parameter ::ca=1d0/dsqrt(3d0)
    
    PetscErrorCode :: ierr
    
endmodule
program main   
#include <petsc/finclude/petscksp.h>
                
    use petscksp
    ! use mpi
    use petscisdef
    use md_init
    use md_params
    use md_vector, only :allocate_base,allocate_invec
    use md_dmatrix
    use md_integration
    use md_tendency
    use md_io
    use initial_condition

    implicit none
    PetscLogDouble :: stime,etime
    PetscInt :: zs  
    PetscBool :: flg
    
    call PetscInitialize(PETSC_NULL_CHARACTER, ierr)
    call PetscOptionsGetInt(PETSC_NULL_OPTIONS, PETSC_NULL_CHARACTER, "-n", zs, flg, ierr)          

    print*,'steps',steps,'rec',rec,'output_time',output_time
    call init()

    !set intial condition
    if(notc=='2')then
        call case2()
    elseif(notc=='6')then
        call case6()
    endif

    call output_init()
    call init_cube()


    !integration
    call PetscTime(stime,ierr)
    print*,'start integration,dt=',dt 
    call RK4() 
    call PetscTime(etime,ierr)
    print*,etime-stime,'sec'
    call output()
    call PetscFinalize(ierr)
endprogram
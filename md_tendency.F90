module md_tendency
    implicit none

contains
    subroutine discretization_2d(Dx,x,dxa,A,B,size)
#include <petsc/finclude/petscksp.h>
        use petscmat
        
        implicit none
        integer size,i
        double precision A(size),B(size)
        PetscScalar, pointer :: xl(:),dxal(:)
        PetscErrorCode :: ierr
        Mat:: Dx
        Vec:: x,dxa,ix
        
        PetscCallA(VecGetArrayF90(x,xl,ierr))
        do i=1,size
            xl(i)=A(i)
        enddo
        

        PetscCallA(VecRestoreArrayF90(x,xl,ierr))
        PetscCallA(Matmult(Dx,x,dxa,ierr))
        PetscCallA(VecGetArrayF90(dxa,dxal,ierr))
        do i=1,size
            B(i)=dxal(i)
        enddo
        call VecRestoreArrayF90(dxa,dxal,ierr)
        
    endsubroutine



    subroutine tendency(size,Dx,Dy,hiv,x,dxa,cui,cvi,ui,vi,hi,col,du,dv,dh,gi)
        use md_vector
        use md_dmatrix
        use md_params
        implicit none
        integer size,pn
        Mat:: Dx,Dy,hiv
        Vec:: x,dxa
        double precision,dimension(size)::vor,E,dudy,dvdx,dhdx,dhdy,gi,&
        cui,cvi,ui,vi,hi,du,dv,dh,col,dedx,dedy,hivu,hivv,divx,divy,hivh


        call discretization_2d(Dy,x,dxa,cui,dudy,size)
        call discretization_2d(Dx,x,dxa,cvi,dvdx,size)

        vor=(dvdx-dudy)/gi
        E=(((cui*ui)+(cvi*vi))/2d0)+(g*hi)

        call discretization_2d(Dx,x,dxa,E,dedx,size)
        call discretization_2d(hiv,x,dxa,cui,hivu,size)
        call discretization_2d(Dy,x,dxa,E,dedy,size)
        call discretization_2d(hiv,x,dxa,cvi,hivv,size)
        du=(-dedx+gi*vi*(col+vor)+hivu)*dt
        dv=(-dedy-gi*ui*(col+vor)+hivv)*dt
        

        call discretization_2d(Dx,x,dxa,gi*ui*hi,divx,size)
        call discretization_2d(Dy,x,dxa,gi*vi*hi,divy,size)
        call discretization_2d(hiv,x,dxa,hi,hivh,size)

        dh=(-(divx+divy)/gi+hivh)*dt

    endsubroutine

endmodule
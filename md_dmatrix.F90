module md_dmatrix
#include <petsc/finclude/petscksp.h>
    use md_vector,only :pvec_size1,pvec_size2,pvec_size3,pvec_size4,&
    pvec_size5,pvec_size6,inda,indb,indc,dxb,dyb,hivb,plane1
    use md_params
    use petscksp
    ! use mpi
    use petscisdef

    implicit none
    Mat ::Dx1,Dx2,Dx3,Dx4,Dx5,Dx6
    Mat ::Dy1,Dy2,Dy3,Dy4,Dy5,Dy6
    Mat ::hiv1,hiv2,hiv3,hiv4,hiv5,hiv6
    Vec ::invec1,invec2,invec3,invec4,invec5,invec6,testvec
    Vec ::dxa1,dxa2,dxa3,dxa4,dxa5,dxa6
    
    contains

    subroutine make_inputvector()
#include <petsc/finclude/petscksp.h>
        implicit none
        integer i
        real(kind=8)func
       

        call VecCreateMPI(PETSC_COMM_WORLD,PETSC_DECIDE,pvec_size1,invec1,ierr)
        call VecCreateMPI(PETSC_COMM_WORLD,PETSC_DECIDE,pvec_size2,invec2,ierr)
        call VecCreateMPI(PETSC_COMM_WORLD,PETSC_DECIDE,pvec_size3,invec3,ierr)
        call VecCreateMPI(PETSC_COMM_WORLD,PETSC_DECIDE,pvec_size4,invec4,ierr)
        call VecCreateMPI(PETSC_COMM_WORLD,PETSC_DECIDE,pvec_size5,invec5,ierr)
        call VecCreateMPI(PETSC_COMM_WORLD,PETSC_DECIDE,pvec_size6,invec6,ierr)
        call VecCreateMPI(PETSC_COMM_WORLD,PETSC_DECIDE,pvec_size1,testvec,ierr)
        call VecSetFromOptions(invec1,ierr)
        call VecSetFromOptions(invec2,ierr)
        call VecSetFromOptions(invec3,ierr)
        call VecSetFromOptions(invec4,ierr)
        call VecSetFromOptions(invec5,ierr)
        call VecSetFromOptions(invec6,ierr)
        call VecSetFromOptions(testvec,ierr)
        call VecAssemblyBegin(invec1,ierr)
        call VecAssemblyEnd(  invec1,ierr)
        call VecAssemblyBegin(invec2,ierr)
        call VecAssemblyEnd(  invec2,ierr)
        call VecAssemblyBegin(invec3,ierr)
        call VecAssemblyEnd(  invec3,ierr)
        call VecAssemblyBegin(invec4,ierr)
        call VecAssemblyEnd(  invec4,ierr)
        call VecAssemblyBegin(invec5,ierr)
        call VecAssemblyEnd(  invec5,ierr)
        call VecAssemblyBegin(invec6,ierr)
        call VecAssemblyEnd(  invec6,ierr)
        
        !func=1d0
        !do i=1,pvec_size1
        !    func=1d0+func
        !    call VecSetValue(testvec,i-1,func,INSERT_VALUES,ierr) 
        !enddo
        call VecAssemblyBegin(testvec,ierr)
        call VecAssemblyEnd(testvec,ierr)

        call VecCreateMPI(PETSC_COMM_WORLD,PETSC_DECIDE,pvec_size1,dxa1,ierr)
        call VecCreateMPI(PETSC_COMM_WORLD,PETSC_DECIDE,pvec_size2,dxa2,ierr)
        call VecCreateMPI(PETSC_COMM_WORLD,PETSC_DECIDE,pvec_size3,dxa3,ierr)
        call VecCreateMPI(PETSC_COMM_WORLD,PETSC_DECIDE,pvec_size4,dxa4,ierr)
        call VecCreateMPI(PETSC_COMM_WORLD,PETSC_DECIDE,pvec_size5,dxa5,ierr)
        call VecCreateMPI(PETSC_COMM_WORLD,PETSC_DECIDE,pvec_size6,dxa6,ierr)
        call VecSetFromOptions(dxa1,ierr)
        call VecSetFromOptions(dxa2,ierr)
        call VecSetFromOptions(dxa3,ierr)
        call VecSetFromOptions(dxa4,ierr)
        call VecSetFromOptions(dxa5,ierr)
        call VecSetFromOptions(dxa6,ierr)
        call VecAssemblyBegin(dxa1,ierr)
        call VecAssemblyEnd(  dxa1,ierr)
        call VecAssemblyBegin(dxa2,ierr)
        call VecAssemblyEnd(  dxa2,ierr)
        call VecAssemblyBegin(dxa3,ierr)
        call VecAssemblyEnd(  dxa3,ierr)
        call VecAssemblyBegin(dxa4,ierr)
        call VecAssemblyEnd(  dxa4,ierr)
        call VecAssemblyBegin(dxa5,ierr)
        call VecAssemblyEnd(  dxa5,ierr)
        call VecAssemblyBegin(dxa6,ierr)
        call VecAssemblyEnd(  dxa6,ierr)

    endsubroutine




    subroutine load_dmatrix()
#include <petsc/finclude/petscksp.h>
        use md_vector, only:edx1,edy1,row1,coll1
        implicit none
        integer i
        double precision c,d,e,f
        PetscInt a,b
        PetscScalar pc,pd,pe
        
        open(11,file='dmatrixP1.txt',status='old')
        open(12,file='dmatrixP2.txt',status='old')
        open(13,file='dmatrixP3.txt',status='old')
        open(14,file='dmatrixP4.txt',status='old')
        open(15,file='dmatrixP5.txt',status='old')
        open(16,file='dmatrixP6.txt',status='old')
        
        !P1の設定
        call MatCreate(PETSC_COMM_WORLD,Dx1 ,ierr)
        call MatCreate(PETSC_COMM_WORLD,Dy1 ,ierr)
        call MatCreate(PETSC_COMM_WORLD,hiv1 ,ierr)
        call MatSetSizes(Dx1, PETSC_DECIDE,PETSC_DECIDE, pvec_size1, pvec_size1,ierr)
        call MatSetSizes(Dy1,PETSC_DECIDE,PETSC_DECIDE, pvec_size1, pvec_size1,ierr)
        call MatSetSizes(hiv1,PETSC_DECIDE,PETSC_DECIDE, pvec_size1, pvec_size1,ierr)
        call MatSetUp(Dx1,ierr);CHKERRA(ierr)
        call MatSetUp(Dy1,ierr);CHKERRA(ierr)
        call MatSetUp(hiv1,ierr);CHKERRA(ierr)
        
        do i=1,plane1*ss
            !print*,i
            read(11,*)a,b,c,d,e
            row1(i)=a;coll1(i)=b;edx1(i)=c;edy1(i)=d
            f=e*gamma
            call Matsetvalue(Dx1, a, b, c, INSERT_VALUES, ierr)
            call Matsetvalue(Dy1, a, b, d, INSERT_VALUES, ierr)
            call Matsetvalue(hiv1,a, b, e*gamma, INSERT_VALUES, ierr)
        enddo
        
        !201 continue
        call MatAssemblyBegin(Dx1, MAT_FINAL_ASSEMBLY, ierr)
        call MatAssemblyEND(Dx1, MAT_FINAL_ASSEMBLY, ierr)
        call MatAssemblyBegin(Dy1, MAT_FINAL_ASSEMBLY, ierr)
        call MatAssemblyEND(Dy1, MAT_FINAL_ASSEMBLY, ierr)
        call MatAssemblyBegin(hiv1, MAT_FINAL_ASSEMBLY, ierr)
        call MatAssemblyEND(hiv1, MAT_FINAL_ASSEMBLY, ierr)
        
        !P2の設定
        call MatCreate(PETSC_COMM_WORLD,Dx2 ,ierr)
        call MatCreate(PETSC_COMM_WORLD,Dy2 ,ierr)
        call MatCreate(PETSC_COMM_WORLD,hiv2 ,ierr)
        call MatSetSizes(Dx2, PETSC_DECIDE,PETSC_DECIDE, pvec_size2, pvec_size2,ierr)
        call MatSetSizes(Dy2,PETSC_DECIDE,PETSC_DECIDE, pvec_size2, pvec_size2,ierr)
        call MatSetSizes(hiv2,PETSC_DECIDE,PETSC_DECIDE, pvec_size2, pvec_size2,ierr)
        call MatSetUp(Dx2,ierr);CHKERRA(ierr)
        call MatSetUp(Dy2,ierr);CHKERRA(ierr)
        call MatSetUp(hiv2,ierr);CHKERRA(ierr)
        do i=1,ie
            read(12,*,end=202)a,b,c,d,e
            call Matsetvalue(Dx2, a, b, c, INSERT_VALUES, ierr)
            call Matsetvalue(Dy2, a, b, d, INSERT_VALUES, ierr)
            call Matsetvalue(hiv2,a, b, e*gamma, INSERT_VALUES, ierr)
        enddo
        202 continue
        call MatAssemblyBegin(Dx2, MAT_FINAL_ASSEMBLY, ierr)
        call MatAssemblyEND(Dx2, MAT_FINAL_ASSEMBLY, ierr)
        call MatAssemblyBegin(Dy2, MAT_FINAL_ASSEMBLY, ierr)
        call MatAssemblyEND(Dy2, MAT_FINAL_ASSEMBLY, ierr)
        call MatAssemblyBegin(hiv2, MAT_FINAL_ASSEMBLY, ierr)
        call MatAssemblyEND(hiv2, MAT_FINAL_ASSEMBLY, ierr)

        !P3の設定
        call MatCreate(PETSC_COMM_WORLD,Dx3 ,ierr)
        call MatCreate(PETSC_COMM_WORLD,Dy3 ,ierr)
        call MatCreate(PETSC_COMM_WORLD,hiv3 ,ierr)
        call MatSetSizes(Dx3, PETSC_DECIDE,PETSC_DECIDE, pvec_size3, pvec_size3,ierr)
        call MatSetSizes(Dy3,PETSC_DECIDE,PETSC_DECIDE,  pvec_size3, pvec_size3,ierr)
        call MatSetSizes(hiv3,PETSC_DECIDE,PETSC_DECIDE, pvec_size3, pvec_size3,ierr)
        call MatSetUp( Dx3,ierr);CHKERRA(ierr)
        call MatSetUp( Dy3,ierr);CHKERRA(ierr)
        call MatSetUp(hiv3,ierr);CHKERRA(ierr)
        do i=1,ie
            read(13,*,end=203)a,b,c,d,e
            call Matsetvalue( Dx3, a, b, c, INSERT_VALUES, ierr)
            call Matsetvalue( Dy3, a, b, d, INSERT_VALUES, ierr)
            call Matsetvalue(hiv3,a, b, e*gamma, INSERT_VALUES, ierr)
        enddo
        203 continue
        call MatAssemblyBegin(Dx3, MAT_FINAL_ASSEMBLY, ierr)
        call MatAssemblyEND(  Dx3, MAT_FINAL_ASSEMBLY, ierr)
        call MatAssemblyBegin(Dy3, MAT_FINAL_ASSEMBLY, ierr)
        call MatAssemblyEND(  Dy3, MAT_FINAL_ASSEMBLY, ierr)
        call MatAssemblyBegin(hiv3, MAT_FINAL_ASSEMBLY, ierr)
        call MatAssemblyEND(  hiv3, MAT_FINAL_ASSEMBLY, ierr)

        !P4の設定
        call MatCreate(PETSC_COMM_WORLD,Dx4 ,ierr)
        call MatCreate(PETSC_COMM_WORLD,Dy4 ,ierr)
        call MatCreate(PETSC_COMM_WORLD,hiv4 ,ierr)
        call MatSetSizes(Dx4, PETSC_DECIDE,PETSC_DECIDE, pvec_size4, pvec_size4,ierr)
        call MatSetSizes(Dy4,PETSC_DECIDE,PETSC_DECIDE,  pvec_size4, pvec_size4,ierr)
        call MatSetSizes(hiv4,PETSC_DECIDE,PETSC_DECIDE, pvec_size4, pvec_size4,ierr)
        call MatSetUp( Dx4,ierr);CHKERRA(ierr)
        call MatSetUp( Dy4,ierr);CHKERRA(ierr)
        call MatSetUp(hiv4,ierr);CHKERRA(ierr)
        do i=1,ie
            read(14,*,end=204)a,b,c,d,e
            call Matsetvalue( Dx4, a, b, c, INSERT_VALUES, ierr)
            call Matsetvalue( Dy4, a, b, d, INSERT_VALUES, ierr)
            call Matsetvalue(hiv4,a, b, e*gamma, INSERT_VALUES, ierr)
        enddo
        204 continue
        call MatAssemblyBegin(Dx4, MAT_FINAL_ASSEMBLY, ierr)
        call MatAssemblyEND(  Dx4, MAT_FINAL_ASSEMBLY, ierr)
        call MatAssemblyBegin(Dy4, MAT_FINAL_ASSEMBLY, ierr)
        call MatAssemblyEND(  Dy4, MAT_FINAL_ASSEMBLY, ierr)
        call MatAssemblyBegin(hiv4, MAT_FINAL_ASSEMBLY, ierr)
        call MatAssemblyEND(  hiv4, MAT_FINAL_ASSEMBLY, ierr)

        !P5の設定
        call MatCreate(PETSC_COMM_WORLD,Dx5 ,ierr)
        call MatCreate(PETSC_COMM_WORLD,Dy5 ,ierr)
        call MatCreate(PETSC_COMM_WORLD,hiv5 ,ierr)
        call MatSetSizes(Dx5, PETSC_DECIDE,PETSC_DECIDE, pvec_size5, pvec_size5,ierr)
        call MatSetSizes(Dy5,PETSC_DECIDE,PETSC_DECIDE,  pvec_size5, pvec_size5,ierr)
        call MatSetSizes(hiv5,PETSC_DECIDE,PETSC_DECIDE, pvec_size5, pvec_size5,ierr)
        call MatSetUp( Dx5,ierr);CHKERRA(ierr)
        call MatSetUp( Dy5,ierr);CHKERRA(ierr)
        call MatSetUp(hiv5,ierr);CHKERRA(ierr)
        do i=1,ie
            read(15,*,end=205)a,b,c,d,e
            call Matsetvalue( Dx5, a, b, c, INSERT_VALUES, ierr)
            call Matsetvalue( Dy5, a, b, d, INSERT_VALUES, ierr)
            call Matsetvalue(hiv5,a, b, e*gamma, INSERT_VALUES, ierr)
        enddo
        205 continue
        call MatAssemblyBegin(Dx5, MAT_FINAL_ASSEMBLY, ierr)
        call MatAssemblyEND(  Dx5, MAT_FINAL_ASSEMBLY, ierr)
        call MatAssemblyBegin(Dy5, MAT_FINAL_ASSEMBLY, ierr)
        call MatAssemblyEND(  Dy5, MAT_FINAL_ASSEMBLY, ierr)
        call MatAssemblyBegin(hiv5, MAT_FINAL_ASSEMBLY, ierr)
        call MatAssemblyEND(  hiv5, MAT_FINAL_ASSEMBLY, ierr)

        !P6の設定
        call MatCreate(PETSC_COMM_WORLD,Dx6 ,ierr)
        call MatCreate(PETSC_COMM_WORLD,Dy6 ,ierr)
        call MatCreate(PETSC_COMM_WORLD,hiv6 ,ierr)
        call MatSetSizes(Dx6, PETSC_DECIDE,PETSC_DECIDE, pvec_size6, pvec_size6,ierr)
        call MatSetSizes(Dy6,PETSC_DECIDE,PETSC_DECIDE,  pvec_size6, pvec_size6,ierr)
        call MatSetSizes(hiv6,PETSC_DECIDE,PETSC_DECIDE, pvec_size6, pvec_size6,ierr)
        call MatSetUp( Dx6,ierr);CHKERRA(ierr)
        call MatSetUp( Dy6,ierr);CHKERRA(ierr)
        call MatSetUp(hiv6,ierr);CHKERRA(ierr)
        do i=1,ie
            read(16,*,end=206)a,b,c,d,e
            call Matsetvalue( Dx6, a, b, c, INSERT_VALUES, ierr)
            call Matsetvalue( Dy6, a, b, d, INSERT_VALUES, ierr)
            call Matsetvalue(hiv6,a, b, e*gamma, INSERT_VALUES, ierr)
        enddo
        206 continue

        call MatAssemblyBegin(Dx6, MAT_FINAL_ASSEMBLY, ierr)
        call MatAssemblyEND(  Dx6, MAT_FINAL_ASSEMBLY, ierr)
        call MatAssemblyBegin(Dy6, MAT_FINAL_ASSEMBLY, ierr)
        call MatAssemblyEND(  Dy6, MAT_FINAL_ASSEMBLY, ierr)
        call MatAssemblyBegin(hiv6, MAT_FINAL_ASSEMBLY, ierr)
        call MatAssemblyEND(  hiv6, MAT_FINAL_ASSEMBLY, ierr)

        


    endsubroutine
endmodule
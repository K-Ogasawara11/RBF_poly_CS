module md_io 
    implicit none
    

contains


    subroutine load_vec()
#include <petsc/finclude/petscksp.h>
        use md_vector ,only: fl,lon,lat,nx,ny,nz,wgt
        use md_params
        use petscksp
 !       ! use mpi
        use petscisdef
        implicit none
        integer i,count
        real(kind=8) a,b,c,d,e,f
        character(len=20)nnn  
        !print*,'level*n=',n*levels
        write(nnn,*)n
        open(40,file="../nodes/knn_"//trim(adjustl(node))//"_N"//trim(adjustl(nnn))//".txt",status='old')
        
       
        
        do i=1,n
            read (40,*,end=300)a,b,c,d,e,f
            lon(i)=a
            lat(i)=b
            nx(i)=c
            ny(i)=d
            nz(i)=e
            wgt(i)=f

        end do
        300 continue
    
        close(40)
    
        fl=2d0*omega*dsin(lat)
    endsubroutine

    subroutine load_cindex()
        use md_vector, only: cindexb,cindex1,cindex2,cindex3&
        ,cindex4,cindex5,cindex6,pvec_size1,pvec_size2,pvec_size3&
        ,pvec_size4,pvec_size5,pvec_size6
        use md_params
        implicit none
        integer i,count
        integer(kind=4) a,b,c,d,e,f
        open(11,file='cindex1.txt',status='old')
        open(12,file='cindex2.txt',status='old')
        open(13,file='cindex3.txt',status='old')
        open(14,file='cindex4.txt',status='old')
        open(15,file='cindex5.txt',status='old')
        open(16,file='cindex6.txt',status='old')

        count=0
        do i=1,n
            read (11,*,end=200)a
            cindexb(i)=a+1
            count=count+1
        end do
        200 continue

        allocate(cindex1(count))
        do i=1,count
            cindex1(i)=cindexb(i)
        end do
        pvec_size1=count
        !print*,cindex1(1:10)

        count=0
        do i=1,n
            read (12,*,end=201)a
            cindexb(i)=a+1
            count=count+1
        end do
        201 continue
        
        
        allocate(cindex2(count))
        do i=1,count
            cindex2(i)=cindexb(i)
        end do
        pvec_size2=count
        count=0
        do i=1,n
            read (13,*,end=202)a
            cindexb(i)=a+1
            count=count+1
        end do
        202 continue

        allocate(cindex3(count))
        do i=1,count
            cindex3(i)=cindexb(i)
        end do
        pvec_size3=count

        count=0
        do i=1,n
            read (14,*,end=203)a
            cindexb(i)=a+1
            count=count+1
        end do
        203 continue

        allocate(cindex4(count))
        do i=1,count
            cindex4(i)=cindexb(i)
        end do
        pvec_size4=count

        count=0
        do i=1,n
            read (15,*,end=204)a
            cindexb(i)=a+1
            count=count+1
        end do
        204 continue

        allocate(cindex5(count))
        do i=1,count
            cindex5(i)=cindexb(i)
        end do
        pvec_size5=count

        count=0
        do i=1,n
            read (16,*,end=205)a
            cindexb(i)=a+1
            count=count+1
        end do
        205 continue

        allocate(cindex6(count))
        do i=1,count
            cindex6(i)=cindexb(i)
        end do 
        deallocate(cindexb)
        pvec_size6=count

        close(11)
        close(12)
        close(13)
        close(14)
        close(15)
        close(16)
    endsubroutine

    subroutine load_overlaparea()
        use md_vector, only:overlapareab,overlaparea,olon,olat,us,vs,ho&
        ,oversize
        use md_params, only:n
        implicit none 
        integer i,count
        integer(kind=4) a,b
        open(11,file='overlaparea.txt',status='old')
        count=0
        do i=1,n
            read(11,*,end=201)a
            overlapareab(i)=a
            count=count+1
        enddo
        201 continue
        allocate(overlaparea(count),olon(count),olat(count),&
        us(count),vs(count),ho(count))
        oversize=count
        do i=1,count
            overlaparea(i)=overlapareab(i)+1
        enddo
        close(11)
        deallocate(overlapareab)

    endsubroutine

    subroutine load_b2overlap()
        use md_vector, only:b2overlap1,b2overlap2,b2overlap3&
        ,b2overlap4,b2overlap5,b2overlap6,b2overlapb,&
        c2o_size1,c2o_size2,c2o_size3,c2o_size4,c2o_size5,c2o_size6
        use md_params, only:n
        implicit none 
        integer i,count
        integer(kind=4) a,b
        open(11,file='b2overlap1.txt',status='old')
        open(12,file='b2overlap2.txt',status='old')
        open(13,file='b2overlap3.txt',status='old')
        open(14,file='b2overlap4.txt',status='old')
        open(15,file='b2overlap5.txt',status='old')
        open(16,file='b2overlap6.txt',status='old')

        count=0
        do i=1,n
            read(11,*,end=201)a
            b2overlapb(i)=a
            count=count+1
        enddo
        201 continue
        c2o_size1=count
        allocate(b2overlap1(count))
        do i=1,count
            b2overlap1(i)=b2overlapb(i)+1
        enddo

        count=0
        do i=1,n
            read(12,*,end=202)a
            b2overlapb(i)=a
            count=count+1
        enddo
        202 continue
        c2o_size2=count
        allocate(b2overlap2(count))
        do i=1,count
            b2overlap2(i)=b2overlapb(i)+1
        enddo

        count=0
        do i=1,n
            read(13,*,end=203)a
            b2overlapb(i)=a
            count=count+1
        enddo
        203 continue
        c2o_size3=count
        allocate(b2overlap3(count))
        do i=1,count
            b2overlap3(i)=b2overlapb(i)+1
        enddo

        count=0
        do i=1,n
            read(14,*,end=204)a
            b2overlapb(i)=a
            count=count+1
        enddo
        204 continue
        c2o_size4=count
        allocate(b2overlap4(count))
        do i=1,count
            b2overlap4(i)=b2overlapb(i)+1
        enddo

        count=0
        do i=1,n
            read(15,*,end=205)a
            b2overlapb(i)=a
            count=count+1
        enddo
        205 continue
        c2o_size5=count
        allocate(b2overlap5(count))
        do i=1,count
            b2overlap5(i)=b2overlapb(i)+1
        enddo

        count=0
        do i=1,n
            read(16,*,end=206)a
            b2overlapb(i)=a
            count=count+1
        enddo
        206 continue
        c2o_size6=count
        allocate(b2overlap6(count))
        do i=1,count
            b2overlap6(i)=b2overlapb(i)+1
        enddo
        deallocate(b2overlapb)
        close(11)
        close(12)
        close(13)
        close(14)
        close(15)
        close(16)

    endsubroutine


    subroutine load_cidim_and_c2in()
        use md_vector, only:cidim1,cidim2,cidim3,cidim4,cidim5,cidim6,cidimb&
        ,c2innerL1,c2innerL2,c2innerL3,c2innerL4,c2innerL5,c2innerL6,c2innerLb&
        ,plane1,plane2,plane3,plane4,plane5,plane6
        use md_params, only:n
        implicit none 
        integer i,count
        integer(kind=4) a,b
        open(11,file='cidim+c2innerLP1.txt',status='old')
        open(12,file='cidim+c2innerLP2.txt',status='old')
        open(13,file='cidim+c2innerLP3.txt',status='old')
        open(14,file='cidim+c2innerLP4.txt',status='old')
        open(15,file='cidim+c2innerLP5.txt',status='old')
        open(16,file='cidim+c2innerLP6.txt',status='old')

        count=0
        do i=1,n
            read(11,*,end=201)a,b
            cidimb(i)=a
            c2innerLb(i)=b
            count=count+1
        enddo
        201 continue
        plane1=count
        allocate(cidim1(count),c2innerL1(count))
        do i=1,count
            cidim1(i)=cidimb(i)+1
            c2innerL1(i)=c2innerLb(i)+1
        enddo
        
        count=0
        do i=1,n
            read(12,*,end=202)a,b
            cidimb(i)=a
            c2innerLb(i)=b
            count=count+1
        enddo
        202 continue
        plane2=count
        allocate(cidim2(count),c2innerL2(count))
        do i=1,count
            cidim2(i)=cidimb(i)+1
            c2innerL2(i)=c2innerLb(i)+1
        enddo

        count=0
        do i=1,n
            read(13,*,end=203)a,b
            cidimb(i)=a
            c2innerLb(i)=b
            count=count+1
        enddo
        203 continue
        plane3=count
        allocate(cidim3(count),c2innerL3(count))
        do i=1,count
            cidim3(i)=cidimb(i)+1
            c2innerL3(i)=c2innerLb(i)+1
        enddo

        count=0
        do i=1,n
            read(14,*,end=204)a,b
            cidimb(i)=a
            c2innerLb(i)=b
            count=count+1
        enddo
        204 continue
        plane4=count
        allocate(cidim4(count),c2innerL4(count))
        do i=1,count
            cidim4(i)=cidimb(i)+1
            c2innerL4(i)=c2innerLb(i)+1
        enddo

        count=0
        do i=1,n
            read(15,*,end=205)a,b
            cidimb(i)=a
            c2innerLb(i)=b
            count=count+1
        enddo
        205 continue
        plane5=count
        allocate(cidim5(count),c2innerL5(count))
        do i=1,count
            cidim5(i)=cidimb(i)+1
            c2innerL5(i)=c2innerLb(i)+1
        enddo

        count=0
        do i=1,n
            read(16,*,end=206)a,b
            cidimb(i)=a
            c2innerLb(i)=b
            count=count+1
        enddo
        206 continue
        plane6=count
        allocate(cidim6(count),c2innerL6(count))
        do i=1,count
            cidim6(i)=cidimb(i)+1
            c2innerL6(i)=c2innerLb(i)+1
        enddo
        deallocate(cidimb,c2innerLb)
        close(11)
        close(12)
        close(13)
        close(14)
        close(15)
        close(16)


    endsubroutine

    subroutine load_over2in()
        use md_vector, only:over2inL1,over2inL2,over2inL3,over2inL4,over2inL5,over2inL6&
        ,over2inLb,over2inR1,over2inR2,over2inR3,over2inR4,over2inR5,over2inR6,over2inRb&
        ,o2in_size1,o2in_size2,o2in_size3,o2in_size4,o2in_size5,o2in_size6
        use md_params, only:n
        implicit none 
        integer i,count
        integer(kind=4) a,b
        open(11,file='over2inP1.txt',status='old')
        open(12,file='over2inP2.txt',status='old')
        open(13,file='over2inP3.txt',status='old')
        open(14,file='over2inP4.txt',status='old')
        open(15,file='over2inP5.txt',status='old')
        open(16,file='over2inP6.txt',status='old')

        count=0
        do i=1,n
            read(11,*,end=201)a,b
            over2inRb(i)=b
            over2inLb(i)=a
            count=count+1
        enddo
        201 continue
        o2in_size1=count
        allocate(over2inR1(count),over2inL1(count))
        do i=1,count
            over2inR1(i)=over2inRb(i)+1
            over2inL1(i)=over2inLb(i)+1
        enddo
        
        count=0
        do i=1,n
            read(12,*,end=202)a,b
            over2inRb(i)=b
            over2inLb(i)=a
            count=count+1
        enddo
        202 continue
        o2in_size2=count
        allocate(over2inR2(count),over2inL2(count))
        do i=1,count
            over2inR2(i)=over2inRb(i)+1
            over2inL2(i)=over2inLb(i)+1
        enddo

        count=0
        do i=1,n
            read(13,*,end=203)a,b
            over2inRb(i)=b
            over2inLb(i)=a
            count=count+1
        enddo
        203 continue
        o2in_size3=count
        allocate(over2inR3(count),over2inL3(count))
        do i=1,count
            over2inR3(i)=over2inRb(i)+1
            over2inL3(i)=over2inLb(i)+1
        enddo

        count=0
        do i=1,n
            read(14,*,end=204)a,b
            over2inRb(i)=b
            over2inLb(i)=a
            count=count+1
        enddo
        204 continue
        o2in_size4=count
        allocate(over2inR4(count),over2inL4(count))
        do i=1,count
            over2inR4(i)=over2inRb(i)+1
            over2inL4(i)=over2inLb(i)+1
        enddo

        count=0
        do i=1,n
            read(15,*,end=205)a,b
            over2inRb(i)=b
            over2inLb(i)=a
            count=count+1
        enddo
        205 continue
        o2in_size5=count
        allocate(over2inR5(count),over2inL5(count))
        do i=1,count
            over2inR5(i)=over2inRb(i)+1
            over2inL5(i)=over2inLb(i)+1
        enddo

        count=0
        do i=1,n
            read(16,*,end=206)a,b
            over2inRb(i)=b
            over2inLb(i)=a
            count=count+1
        enddo
        206 continue
        o2in_size6=count
        allocate(over2inR6(count),over2inL6(count))
        do i=1,count
            over2inR6(i)=over2inRb(i)+1
            over2inL6(i)=over2inLb(i)+1
        enddo

        close(11)
        close(12)
        close(13)
        close(14)
        close(15)
        close(16)
        deallocate(over2inLb,over2inRb)


    endsubroutine

    subroutine load_b2R()
        use md_vector, only:c2bR1,c2bR2,c2bR3,c2bR4,c2bR5,c2bR6&
        ,c2bRb,in2bR1,in2bR2,in2bR3,in2bR4,in2bR5,in2bR6,in2bRb
        use md_params, only:n
        implicit none 
        integer i,count
        integer(kind=4) a,b
        open(11,file='2bR1.txt',status='old')
        open(12,file='2bR2.txt',status='old')
        open(13,file='2bR3.txt',status='old')
        open(14,file='2bR4.txt',status='old')
        open(15,file='2bR5.txt',status='old')
        open(16,file='2bR6.txt',status='old')

        count=0
        do i=1,n
            read(11,*,end=201)a,b
            in2bRb(i)=b
            c2bRb(i)=a
            count=count+1
        enddo
        201 continue
        allocate(in2bR1(count),c2bR1(count))
        do i=1,count
            in2bR1(i)=in2bRb(i)+1
            c2bR1(i)=c2bRb(i)+1
        enddo
        
        count=0
        do i=1,n
            read(12,*,end=202)a,b
            in2bRb(i)=b
            c2bRb(i)=a
            count=count+1
        enddo
        202 continue
        allocate(in2bR2(count),c2bR2(count))
        do i=1,count
            in2bR2(i)=in2bRb(i)+1
            c2bR2(i)=c2bRb(i)+1
        enddo

        count=0
        do i=1,n
            read(13,*,end=203)a,b
            in2bRb(i)=b
            c2bRb(i)=a
            count=count+1
        enddo
        203 continue
        allocate(in2bR3(count),c2bR3(count))
        do i=1,count
            in2bR3(i)=in2bRb(i)+1
            c2bR3(i)=c2bRb(i)+1
        enddo

        count=0
        do i=1,n
            read(14,*,end=204)a,b
            in2bRb(i)=b
            c2bRb(i)=a
            count=count+1
        enddo
        204 continue
        allocate(in2bR4(count),c2bR4(count))
        do i=1,count
            in2bR4(i)=in2bRb(i)+1
            c2bR4(i)=c2bRb(i)+1
        enddo

        count=0
        do i=1,n
            read(15,*,end=205)a,b
            in2bRb(i)=b
            c2bRb(i)=a
            count=count+1
        enddo
        205 continue
        allocate(in2bR5(count),c2bR5(count))
        do i=1,count
            in2bR5(i)=in2bRb(i)+1
            c2bR5(i)=c2bRb(i)+1
        enddo

        count=0
        do i=1,n
            read(16,*,end=206)a,b
            in2bRb(i)=b
            c2bRb(i)=a
            count=count+1
        enddo
        206 continue
        allocate(in2bR6(count),c2bR6(count))
        do i=1,count
            in2bR6(i)=in2bRb(i)+1
            c2bR6(i)=c2bRb(i)+1
        enddo

        deallocate(in2bRb,c2bRb)

        close(11)
        close(12)
        close(13)
        close(14)
        close(15)
        close(16)


    endsubroutine

    subroutine load_surf(hb,Fhb,n,nproc,myid,ie,pnn)
        implicit none
        integer i,j,n,pnn,nproc,myid,geta,ss,ie,inda(ie),indb(ie)
        real(kind=8) pi,is,a,b,c,d,e
        real(kind = 8)Fhb(n),hb(n)
        
        character(len=20) latfile,lonfile,nn,sss,ep
        write(nn,*)n
        open(20+myid,file="hs_"//trim(adjustl(nn))//".txt",status='old')
        
       
        do i=1,n
            read (20+myid,*,end=200)a
            Fhb(i)=a
            hb(i)=hb(i)-a
        end do
        
        200 continue
        
        
        close(20+myid)
    endsubroutine

    

    subroutine output()
        use md_vector
        use md_params
        implicit none
        double precision, dimension(n)::acu,acv,ah 
        character(len=20)wn
        integer i
        write(wn,*)n
        acu(cidim1)=cui1(c2innerL1);acv(cidim1)=cvi1(c2innerL1)
        acu(cidim2)=cui2(c2innerL2);acv(cidim2)=cvi2(c2innerL2)
        acu(cidim3)=cui3(c2innerL3);acv(cidim3)=cvi3(c2innerL3)
        acu(cidim4)=cui4(c2innerL4);acv(cidim4)=cvi4(c2innerL4)
        acu(cidim5)=cui5(c2innerL5);acv(cidim5)=cvi5(c2innerL5)
        acu(cidim6)=cui6(c2innerL6);acv(cidim6)=cvi6(c2innerL6)
        
        ah(cidim1)=hi1(c2innerL1);ah(cidim2)=hi2(c2innerL2);ah(cidim3)=hi3(c2innerL3)
        ah(cidim4)=hi4(c2innerL4);ah(cidim5)=hi5(c2innerL5);ah(cidim6)=hi6(c2innerL6)
        
        open(11,file="case"//trim(adjustl(notc))//"_h_N"//trim(adjustl(wn))//".txt",status='replace')
        open(12,file="case"//trim(adjustl(notc))//"_cu_N"//trim(adjustl(wn))//".txt",status='replace')
        open(13,file="case"//trim(adjustl(notc))//"_cv_N"//trim(adjustl(wn))//".txt",status='replace')

        do i=1,n
            write(11,*)ah(i)
            write(12,*)acu(i)
            write(13,*)acv(i)
        end do
    endsubroutine

    subroutine output_in_integration(step)
        use md_vector
        use md_params
        implicit none
        double precision, dimension(n)::acu,acv,ah 
        character(len=20)wn,wtime,wss
        integer k,step
        write(wn,*)n
        write(wtime,*)int(step*dt/(60d0*60d0))
        write(wss,*)ss
        acu(cidim1)=cui1(c2innerL1);acv(cidim1)=cvi1(c2innerL1)
        acu(cidim2)=cui2(c2innerL2);acv(cidim2)=cvi2(c2innerL2)
        acu(cidim3)=cui3(c2innerL3);acv(cidim3)=cvi3(c2innerL3)
        acu(cidim4)=cui4(c2innerL4);acv(cidim4)=cvi4(c2innerL4)
        acu(cidim5)=cui5(c2innerL5);acv(cidim5)=cvi5(c2innerL5)
        acu(cidim6)=cui6(c2innerL6);acv(cidim6)=cvi6(c2innerL6)

        ah(cidim1)=hi1(c2innerL1);ah(cidim2)=hi2(c2innerL2);ah(cidim3)=hi3(c2innerL3)
        ah(cidim4)=hi4(c2innerL4);ah(cidim5)=hi5(c2innerL5);ah(cidim6)=hi6(c2innerL6)
        
        open(11,file="case"//trim(adjustl(notc))//"_h_N"//trim(adjustl(wn))//"n"//trim(adjustl(wss))//"_"//trim(adjustl(wtime))//&
        "hour.txt",status='replace')
        open(12,file="case"//trim(adjustl(notc))//"_cu_N"//trim(adjustl(wn))//"n"//trim(adjustl(wss))//"_"//trim(adjustl(wtime))//&
        "hour.txt",status='replace')
        open(13,file="case"//trim(adjustl(notc))//"_cv_N"//trim(adjustl(wn))//"n"//trim(adjustl(wss))//"_"//trim(adjustl(wtime))//&
        "hour.txt",status='replace')
        !print*,'n',n
        do k=1,n
            write(11,*)ah(k)
        end do

        do k=1,n
            write(13,*)acv(k)
        end do

        do k=1,n
            write(12,*)acu(k)
        end do
        close(11)
        close(12)
        close(13)
    endsubroutine


    subroutine output_init()
        use md_vector
        use md_params
        implicit none
        double precision, dimension(n)::acu,acv,ah 
        integer i
        character(len=20)nnn  
        write(nnn,*)n
        
            
            
        open(11,file="case"//trim(adjustl(notc))//"_init_h_N"//trim(adjustl(nnn))//".txt",status='replace')
        open(12,file="case"//trim(adjustl(notc))//"_init_u_N"//trim(adjustl(nnn))//".txt",status='replace')
        open(13,file="case"//trim(adjustl(notc))//"_init_v_N"//trim(adjustl(nnn))//".txt",status='replace')
    
        do i=1,n
            !write(11,*)ah(i)
            write(11,*)phib(i)
            write(12,*)ub(i)
            write(13,*)vb(i)
        end do

        close(11)
        close(12)
        close(13)


        endsubroutine
    
        


    
endmodule md_io

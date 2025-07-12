module md_vector
    use md_params ,only :n,ie
    implicit none
    double precision,allocatable:: nx(:),ny(:),nz(:),lat(:),lon(:)&
    ,fl(:),dudx(:),dxpsl(:),f(:)
    double precision,allocatable::ub(:),vb(:),phib(:),dxb(:),dyb(:),hivb(:)
    double precision,allocatable::visb(:),wgt(:)
    double precision,allocatable::du1(:),du2(:),du3(:),du4(:),du5(:),du6(:)
    double precision,allocatable::dv1(:),dv2(:),dv3(:),dv4(:),dv5(:),dv6(:)
    double precision,allocatable::ui1(:),ui2(:),ui3(:),ui4(:),ui5(:),ui6(:)&
    ,vi1(:),vi2(:),vi3(:),vi4(:),vi5(:),vi6(:)
    double precision,allocatable::us1(:),us2(:),us3(:),us4(:),us5(:),us6(:)&
    ,vs1(:),vs2(:),vs3(:),vs4(:),vs5(:),vs6(:),us(:),vs(:),ho(:)
    double precision,allocatable::cui1(:),cui2(:),cui3(:),cui4(:),cui5(:),cui6(:)&
    ,cvi1(:),cvi2(:),cvi3(:),cvi4(:),cvi5(:),cvi6(:),olat(:),olon(:)
    double precision,allocatable::hi1(:),hi2(:),hi3(:),hi4(:),hi5(:),hi6(:)
    double precision,allocatable::slon1(:),slon2(:),slon3(:),slon4(:),slon5(:),slon6(:)
    double precision,allocatable::slat1(:),slat2(:),slat3(:),slat4(:),slat5(:),slat6(:)
    double precision,allocatable::cslat1(:),cslat2(:),cslat3(:),cslat4(:),cslat5(:),cslat6(:)
    double precision,allocatable::cslon1(:),cslon2(:),cslon3(:),cslon4(:),cslon5(:),cslon6(:)
    double precision,allocatable::sslat1(:),sslat2(:),sslat3(:),sslat4(:),sslat5(:),sslat6(:)
    double precision,allocatable::sslon1(:),sslon2(:),sslon3(:),sslon4(:),sslon5(:),sslon6(:)
    double precision,allocatable::tslat1(:),tslat2(:),tslat3(:),tslat4(:)
    double precision,allocatable::tslon1(:),tslon2(:),tslon3(:),tslon4(:)
    
    
    double precision,allocatable:: edx1(:),edy1(:),edx2(:),edy2(:),edx3(:),edy3(:),edx4(:),edy4(:),edx5(:),edy5(:),edx6(:),edy6(:)

    double precision,allocatable::o2inlat1(:),o2inlat2(:),o2inlat3(:),o2inlat4(:),o2inlat5(:),o2inlat6(:)
    double precision,allocatable::o2inlon1(:),o2inlon2(:),o2inlon3(:),o2inlon4(:),o2inlon5(:),o2inlon6(:)
    double precision,allocatable::cos_o2inlat1(:),cos_o2inlat2(:),cos_o2inlat3(:)&
    ,cos_o2inlat4(:),cos_o2inlat5(:),cos_o2inlat6(:)
    double precision,allocatable::cos_o2inlon1(:),cos_o2inlon2(:),cos_o2inlon3(:)&
    ,cos_o2inlon4(:),cos_o2inlon5(:),cos_o2inlon6(:)
    double precision,allocatable::tan_o2inlat1(:),tan_o2inlat2(:),tan_o2inlat3(:),tan_o2inlat4(:)
    double precision,allocatable::tan_o2inlon1(:),tan_o2inlon2(:),tan_o2inlon3(:),tan_o2inlon4(:)
    double precision,allocatable::sin_o2inlat1(:),sin_o2inlat2(:),sin_o2inlat3(:),&
    sin_o2inlat4(:),sin_o2inlat5(:),sin_o2inlat6(:)
    double precision,allocatable::sin_o2inlon1(:),sin_o2inlon2(:),sin_o2inlon3(:)&
    ,sin_o2inlon4(:),sin_o2inlon5(:),sin_o2inlon6(:)

    double precision,allocatable::c2overlat1(:),c2overlat2(:),c2overlat3(:),c2overlat4(:),c2overlat5(:),c2overlat6(:)
    double precision,allocatable::c2overlon1(:),c2overlon2(:),c2overlon3(:),c2overlon4(:),c2overlon5(:),c2overlon6(:)
    double precision,allocatable::sin_c2overlat1(:),sin_c2overlat2(:),sin_c2overlat3(:),&
    sin_c2overlat4(:),sin_c2overlat5(:),sin_c2overlat6(:)
    double precision,allocatable::sin_c2overlon1(:),sin_c2overlon2(:),sin_c2overlon3(:),&
    sin_c2overlon4(:),sin_c2overlon5(:),sin_c2overlon6(:)
    double precision,allocatable::cos_c2overlat1(:),cos_c2overlat2(:),cos_c2overlat3(:),&
    cos_c2overlat4(:),cos_c2overlat5(:),cos_c2overlat6(:)
    double precision,allocatable::cos_c2overlon1(:),cos_c2overlon2(:),cos_c2overlon3(:),&
    cos_c2overlon4(:),cos_c2overlon5(:),cos_c2overlon6(:)
    double precision,allocatable::g1(:),g2(:),g3(:),g4(:),g5(:),g6(:)
    double precision,allocatable::col1(:),col2(:),col3(:),col4(:),col5(:),col6(:)
    double precision,allocatable::dh1(:),dh2(:),dh3(:),dh4(:),dh5(:),dh6(:)
    integer(kind=4),allocatable:: inda(:),indb(:),indc(:)
    integer(kind=4),allocatable:: overlaparea(:),overlapareab(:)
    integer(kind=4),allocatable:: b2overlap1(:),b2overlap2(:),b2overlap3(:),b2overlap4(:)&
    ,b2overlap5(:),b2overlap6(:),b2overlapb(:)
    integer(kind=4),allocatable:: cindex1(:),cindex2(:),cindex3(:),&
    cindex4(:),cindex5(:),cindex6(:),cindexb(:)
    integer(kind=4),allocatable:: cidim1(:),cidim2(:),cidim3(:),&
    cidim4(:),cidim5(:),cidim6(:),cidimb(:)
    integer(kind=4),allocatable:: c2innerL1(:),c2innerL2(:),c2innerL3(:),&
    c2innerL4(:),c2innerL5(:),c2innerL6(:),c2innerLb(:)
    integer(kind=4),allocatable:: over2inR1(:),over2inR2(:),over2inR3(:),&
    over2inR4(:),over2inR5(:),over2inR6(:),over2inRb(:)
    integer(kind=4),allocatable:: over2inL1(:),over2inL2(:),over2inL3(:),&
    over2inL4(:),over2inL5(:),over2inL6(:),over2inLb(:)
    integer(kind=4),allocatable:: c2bR1(:),c2bR2(:),c2bR3(:),&
    c2bR4(:),c2bR5(:),c2bR6(:),c2bRb(:)
    integer(kind=4),allocatable:: in2bR1(:),in2bR2(:),in2bR3(:),&
    in2bR4(:),in2bR5(:),in2bR6(:),in2bRb(:)

    integer(kind=4),allocatable::row1(:),coll1(:),row2(:),coll2(:),row3(:),coll3(:)
    integer(kind=4),allocatable::row4(:),coll4(:),row5(:),coll5(:),row6(:),coll6(:)

    integer pvec_size1,pvec_size2,pvec_size3,pvec_size4,pvec_size5,pvec_size6
    integer o2in_size1,o2in_size2,o2in_size3,o2in_size4,o2in_size5,o2in_size6
    integer c2o_size1,c2o_size2,c2o_size3,c2o_size4,c2o_size5,c2o_size6
    integer plane1,plane2,plane3,plane4,plane5,plane6
    integer oversize
    
   


    contains

    subroutine allocate_base()
        
        allocate(nx(n),ny(n),nz(n))
        allocate(ub(n),vb(n),phib(n),fl(n),lon(n),lat(n),wgt(n))
        allocate(dxb(ie),dyb(ie),hivb(ie),visb(ie),inda(ie),indb(ie),indc(ie))

        allocate(b2overlapb(n/2),cindexb(n/2),overlapareab(n/2),cidimb(n/2))
        allocate(in2bRb(n/2),c2bRb(n/2),over2inLb(n/2),over2inRb(n/2),c2innerLb(n/2))
    endsubroutine

    subroutine allocate_invec()
        use md_params
        allocate(us1(pvec_size1),us2(pvec_size2),us3(pvec_size3),us4(pvec_size4))
        allocate(us5(pvec_size5),us6(pvec_size6),vs1(pvec_size1),vs2(pvec_size2),vs3(pvec_size3)&
        ,vs4(pvec_size4),vs5(pvec_size5),vs6(pvec_size6))
        allocate(ui1(pvec_size1),ui2(pvec_size2),ui3(pvec_size3),ui4(pvec_size4))
        allocate(ui5(pvec_size5),ui6(pvec_size6),vi1(pvec_size1),vi2(pvec_size2),vi3(pvec_size3)&
        ,vi4(pvec_size4),vi5(pvec_size5),vi6(pvec_size6))
        allocate(cui1(pvec_size1),cui2(pvec_size2),cui3(pvec_size3),cui4(pvec_size4))
        allocate(cui5(pvec_size5),cui6(pvec_size6),cvi1(pvec_size1),cvi2(pvec_size2),cvi3(pvec_size3)&
        ,cvi4(pvec_size4),cvi5(pvec_size5),cvi6(pvec_size6))
        allocate(hi1(pvec_size1),hi2(pvec_size2),hi3(pvec_size3),hi4(pvec_size4))
        allocate(hi5(pvec_size5),hi6(pvec_size6))
        allocate(dh1(pvec_size1),dh2(pvec_size2),dh3(pvec_size3),dh4(pvec_size4))
        allocate(dh5(pvec_size5),dh6(pvec_size6))
        allocate(g1(pvec_size1),g2(pvec_size2),g3(pvec_size3),g4(pvec_size4))
        allocate(g5(pvec_size5),g6(pvec_size6))
        allocate(col1(pvec_size1),col2(pvec_size2),col3(pvec_size3),col4(pvec_size4))
        allocate(col5(pvec_size5),col6(pvec_size6))
        allocate(du1(pvec_size1),du2(pvec_size2),du3(pvec_size3),du4(pvec_size4))
        allocate(du5(pvec_size5),du6(pvec_size6),dv1(pvec_size1),dv2(pvec_size2),dv3(pvec_size3)&
        ,dv4(pvec_size4),dv5(pvec_size5),dv6(pvec_size6))
        
        
        
        allocate(slon1(pvec_size1),slon2(pvec_size2),slon3(pvec_size3),slon4(pvec_size4))
        allocate(slon5(pvec_size5),slon6(pvec_size6))
        allocate(slat1(pvec_size1),slat2(pvec_size2),slat3(pvec_size3),slat4(pvec_size4))
        allocate(slat5(pvec_size5),slat6(pvec_size6))
        allocate(cslon1(pvec_size1),cslon2(pvec_size2),cslon3(pvec_size3),cslon4(pvec_size4))
        allocate(cslon5(pvec_size5),cslon6(pvec_size6))
        allocate(cslat1(pvec_size1),cslat2(pvec_size2),cslat3(pvec_size3),cslat4(pvec_size4))
        allocate(cslat5(pvec_size5),cslat6(pvec_size6))
        allocate(sslon1(pvec_size1),sslon2(pvec_size2),sslon3(pvec_size3),sslon4(pvec_size4))
        allocate(sslon5(pvec_size5),sslon6(pvec_size6))
        allocate(sslat1(pvec_size1),sslat2(pvec_size2),sslat3(pvec_size3),sslat4(pvec_size4))
        allocate(sslat5(pvec_size5),sslat6(pvec_size6))
        allocate(tslat1(pvec_size1),tslat2(pvec_size2),tslat3(pvec_size3),tslat4(pvec_size4))
        allocate(tslon1(pvec_size1),tslon2(pvec_size2),tslon3(pvec_size3),tslon4(pvec_size4))
        
        allocate(edx1(plane1*ss),edy1(plane1*ss),row1(plane1*ss),coll1(plane1*ss))
        allocate(edx2(plane2*ss),edy2(plane2*ss),row2(plane2*ss),coll2(plane2*ss))
        allocate(edx3(plane3*ss),edy3(plane3*ss),row3(plane3*ss),coll3(plane3*ss))
        allocate(edx4(plane4*ss),edy4(plane4*ss),row4(plane4*ss),coll4(plane4*ss))
        allocate(edx5(plane5*ss),edy5(plane5*ss),row5(plane5*ss),coll5(plane5*ss))
        allocate(edx6(plane6*ss),edy6(plane6*ss),row6(plane6*ss),coll6(plane6*ss))
    endsubroutine

    subroutine prepare_sincos_stensil()
        use md_params
        implicit none
        integer i
        double precision gb
        double precision,allocatable:: slonb(:)
        double precision,dimension(pvec_size2):: col_cslon2,col_cslat2,col_sslat2
        double precision,dimension(pvec_size3):: col_cslon3,col_cslat3,col_sslat3
        double precision,dimension(pvec_size4):: col_cslon4,col_cslat4,col_sslat4
        
        slon1=lon(cindex1);slon2=lon(cindex2)
        slon3=lon(cindex3);slon4=lon(cindex4)
        slon5=lon(cindex5);slon6=lon(cindex6)
        slat1=lat(cindex1);slat2=lat(cindex2)
        slat3=lat(cindex3);slat4=lat(cindex4)
        slat5=lat(cindex5);slat6=lat(cindex6)
        sslat1=dsin(slat1);sslat2=dsin(slat2)
        sslat3=dsin(slat3);sslat4=dsin(slat4)
        sslat5=dsin(slat5);sslat6=dsin(slat6)
        cslat1=dcos(slat1);cslat2=dcos(slat2)
        cslat3=dcos(slat3);cslat4=dcos(slat4)
        cslat5=dcos(slat5);cslat6=dcos(slat6)
        tslat1=dtan(slat1);tslat2=dtan(slat2)
        tslat3=dtan(slat3);tslat4=dtan(slat4)

        col_cslat2=dcos(slat2);col_cslon2=dcos(slon2);col_sslat2=dsin(slat2)
        col_cslat3=dcos(slat3);col_cslon3=dcos(slon3);col_sslat3=dsin(slat3)
        col_cslat4=dcos(slat4);col_cslon4=dcos(slon4);col_sslat4=dsin(slat4)


        sslon1=dsin(slon1);cslon1=dcos(slon1);tslon1=dtan(slon1)
        g1=dcos(slon1)**3d0*dcos(slat1)**3d0/ca**2d0
        allocate(slonb(pvec_size2))
        slonb=slon2-pi/2d0
        sslon2=dsin(slonb);cslon2=dcos(slonb);tslon2=dtan(slonb)
        g2=dcos(slonb)**3d0*dcos(slat2)**3d0/ca**2d0
        deallocate(slonb)
        allocate(slonb(pvec_size3))
        slonb=slon3
        do i=1,pvec_size3
            if(slonb(i)<0d0)then
                slonb(i)=slonb(i)+pi*2d0
            endif
        enddo
        slonb=slonb-pi
        sslon3=dsin(slonb);cslon3=dcos(slonb);tslon3=dtan(slonb)
        g3=dcos(slonb)**3d0*dcos(slat3)**3d0/ca**2d0
        deallocate(slonb)
        allocate(slonb(pvec_size4))
        slonb=slon4+pi/2d0
        g4=dcos(slonb)**3d0*dcos(slat4)**3d0/ca**2d0
        sslon4=dsin(slonb);cslon4=dcos(slonb);tslon4=dtan(slonb)
        sslon5=dsin(slon5);cslon5=dcos(slon5)
        g5=dsin(slat5)**3d0/ca**2d0
        sslon6=dsin(slon6);cslon6=dcos(slon6)
        g6=dsqrt(dsin(slat6)**6d0/ca**4d0)

        col1=(2d0*omega)*((-cslon1*cslat1*dsin(alpha))+(sslat1*dcos(alpha)))
        col2=(2d0*omega)*((-col_cslon2*col_cslat2*dsin(alpha))+(col_sslat2*dcos(alpha)))
        col3=(2d0*omega)*((-col_cslon3*col_cslat3*dsin(alpha))+(col_sslat3*dcos(alpha)))
        col4=(2d0*omega)*((-col_cslon4*col_cslat4*dsin(alpha))+(col_sslat4*dcos(alpha)))
        col5=(2d0*omega)*((-cslon5*cslat5*dsin(alpha))+(sslat5*dcos(alpha)))
        col6=(2d0*omega)*((-cslon6*cslat6*dsin(alpha))+(sslat6*dcos(alpha)))
        

    endsubroutine

    subroutine prepare_sincos_overlap()
        use md_params
        implicit none
        integer i
        double precision,allocatable:: slonb(:),slatb(:),olonb(:),clonb(:)
        allocate(c2overlat1(c2o_size1),c2overlon1(c2o_size1)&
        ,cos_c2overlat1(c2o_size1)&
        ,cos_c2overlon1(c2o_size1),sin_c2overlon1(c2o_size1))
       
        allocate(c2overlat2(c2o_size2),c2overlon2(c2o_size2),&
        cos_c2overlat2(c2o_size2),&
        cos_c2overlon2(c2o_size2),sin_c2overlon2(c2o_size2))
        allocate(c2overlat3(c2o_size3),c2overlon3(c2o_size3),&
        cos_c2overlat3(c2o_size3),&
        cos_c2overlon3(c2o_size3),sin_c2overlon3(c2o_size3))
        allocate(c2overlat4(c2o_size4),c2overlon4(c2o_size4),&
        cos_c2overlat4(c2o_size4),&
        cos_c2overlon4(c2o_size4),sin_c2overlon4(c2o_size4))
        allocate(c2overlat5(c2o_size5),c2overlon5(c2o_size5),&
        cos_c2overlat5(c2o_size5),&
        cos_c2overlon5(c2o_size5),sin_c2overlon5(c2o_size5))
        allocate(c2overlat6(c2o_size6),c2overlon6(c2o_size6)&
        ,cos_c2overlat6(c2o_size6),&
        cos_c2overlon6(c2o_size6),sin_c2overlon6(c2o_size6))

        
        allocate(slatb(plane1),slonb(plane1))
        slatb=lat(cidim1);slonb=lon(cidim1)
        c2overlat1=slatb(c2bR1);c2overlon1=slonb(c2bR1)
        deallocate(slatb,slonb)

        allocate(slatb(plane2),slonb(plane2))
        slatb=lat(cidim2);slonb=lon(cidim2)-pi/2d0
        c2overlat2=slatb(c2bR2);c2overlon2=slonb(c2bR2)
        deallocate(slatb,slonb)

        allocate(slatb(plane3),slonb(plane3))
        slatb=lat(cidim3);slonb=lon(cidim3)
        do i=1,plane3
            if(slonb(i)<0)then
                slonb(i)=slonb(i)+pi*2d0
            endif
        enddo
        slonb=slonb-pi
        c2overlat3=slatb(c2bR3);c2overlon3=slonb(c2bR3)
        deallocate(slatb,slonb)

        allocate(slatb(plane4),slonb(plane4))
        slatb=lat(cidim4);slonb=lon(cidim4)+pi/2d0
        c2overlat4=slatb(c2bR4);c2overlon4=slonb(c2bR4)
        deallocate(slatb,slonb)

        allocate(slatb(plane5),slonb(plane5))
        slatb=lat(cidim5);slonb=lon(cidim5)
        c2overlat5=slatb(c2bR5);c2overlon5=slonb(c2bR5)
        deallocate(slatb,slonb)
        allocate(slatb(plane6),slonb(plane6))
        slatb=lat(cidim6);slonb=lon(cidim6)
        c2overlat6=slatb(c2bR6);c2overlon6=slonb(c2bR6)
        deallocate(slatb,slonb)

        cos_c2overlat1=dcos(c2overlat1);cos_c2overlat2=dcos(c2overlat2)
        cos_c2overlat3=dcos(c2overlat3);cos_c2overlat4=dcos(c2overlat4)
        cos_c2overlat5=dcos(c2overlat5);cos_c2overlat6=dcos(c2overlat6)
        sin_c2overlat1=dsin(c2overlat1);sin_c2overlat2=dsin(c2overlat2)
        sin_c2overlat3=dsin(c2overlat3);sin_c2overlat4=dsin(c2overlat4)
        sin_c2overlat5=dsin(c2overlat5);sin_c2overlat6=dsin(c2overlat6)

        cos_c2overlon1=dcos(c2overlon1);cos_c2overlon2=dcos(c2overlon2)
        cos_c2overlon3=dcos(c2overlon3);cos_c2overlon4=dcos(c2overlon4)
        cos_c2overlon5=dcos(c2overlon5);cos_c2overlon6=dcos(c2overlon6)
        sin_c2overlon1=dsin(c2overlon1);sin_c2overlon2=dsin(c2overlon2)
        sin_c2overlon3=dsin(c2overlon3);sin_c2overlon4=dsin(c2overlon4)
        sin_c2overlon5=dsin(c2overlon5);sin_c2overlon6=dsin(c2overlon6)


        
        allocate(o2inlat1(o2in_size1),o2inlat2(o2in_size2),o2inlat3(o2in_size3),&
        o2inlat4(o2in_size4),o2inlat5(o2in_size5),o2inlat6(o2in_size6))

        allocate(cos_o2inlat1(o2in_size1),cos_o2inlat2(o2in_size2),cos_o2inlat3(o2in_size3),&
        cos_o2inlat4(o2in_size4),cos_o2inlat5(o2in_size5),cos_o2inlat6(o2in_size6))

        allocate(sin_o2inlat1(o2in_size1),sin_o2inlat2(o2in_size2),sin_o2inlat3(o2in_size3),&
        sin_o2inlat4(o2in_size4),sin_o2inlat5(o2in_size5),sin_o2inlat6(o2in_size6))

        allocate(tan_o2inlat1(o2in_size1),tan_o2inlat2(o2in_size2),tan_o2inlat3(o2in_size3),&
        tan_o2inlat4(o2in_size4))

        allocate(o2inlon1(o2in_size1),o2inlon2(o2in_size2),o2inlon3(o2in_size3),&
        o2inlon4(o2in_size4),o2inlon5(o2in_size5),o2inlon6(o2in_size6))

        allocate(sin_o2inlon1(o2in_size1),sin_o2inlon2(o2in_size2),sin_o2inlon3(o2in_size3),&
        sin_o2inlon4(o2in_size4),sin_o2inlon5(o2in_size5),sin_o2inlon6(o2in_size6))

        allocate(cos_o2inlon1(o2in_size1),cos_o2inlon2(o2in_size2),cos_o2inlon3(o2in_size3),&
        cos_o2inlon4(o2in_size4),cos_o2inlon5(o2in_size5),cos_o2inlon6(o2in_size6))

        allocate(tan_o2inlon1(o2in_size1),tan_o2inlon2(o2in_size2),tan_o2inlon3(o2in_size3),&
        tan_o2inlon4(o2in_size4))

        olon=lon(overlaparea)
        olat=lat(overlaparea)
        !print*,'over2inR1',size(over2inR1)
        o2inlat1=olat(over2inR1);o2inlat2=olat(over2inR2)
        o2inlat3=olat(over2inR3);o2inlat4=olat(over2inR4)
        o2inlat5=olat(over2inR5);o2inlat6=olat(over2inR6)

        o2inlon1=olon(over2inR1);o2inlon2=olon(over2inR2)
        o2inlon3=olon(over2inR3);o2inlon4=olon(over2inR4)
        o2inlon5=olon(over2inR5);o2inlon6=olon(over2inR6)
        
        sin_o2inlat1=dsin(o2inlat1);cos_o2inlat1=dcos(o2inlat1);tan_o2inlat1=dtan(o2inlat1)
        sin_o2inlat2=dsin(o2inlat2);cos_o2inlat2=dcos(o2inlat2);tan_o2inlat2=dtan(o2inlat2)
        sin_o2inlat3=dsin(o2inlat3);cos_o2inlat3=dcos(o2inlat3);tan_o2inlat3=dtan(o2inlat3)
        sin_o2inlat4=dsin(o2inlat4);cos_o2inlat4=dcos(o2inlat4);tan_o2inlat4=dtan(o2inlat4)
        sin_o2inlat5=dsin(o2inlat5);cos_o2inlat5=dcos(o2inlat5)
        sin_o2inlat6=dsin(o2inlat6);cos_o2inlat6=dcos(o2inlat6)

        !P1
        sin_o2inlon1=dsin(o2inlon1);cos_o2inlon1=dcos(o2inlon1);tan_o2inlon1=dtan(o2inlon1)
        !P2
        allocate(olonb(o2in_size2))
        olonb=o2inlon2-pi/2d0
        sin_o2inlon2=dsin(olonb);cos_o2inlon2=dcos(olonb);tan_o2inlon2=dtan(olonb)
        deallocate(olonb)

        !P3
        allocate(olonb(o2in_size3))
        olonb=o2inlon3
        do i=1,o2in_size3
            if(olonb(i)<0)then
                olonb(i)=olonb(i)+pi*2d0
            endif
        enddo
        olonb=olonb-pi
        sin_o2inlon3=dsin(olonb);cos_o2inlon3=dcos(olonb);tan_o2inlon3=dtan(olonb)
        deallocate(olonb)

        !P4
        allocate(olonb(o2in_size4))
        olonb=o2inlon4+pi/2d0
        sin_o2inlon4=dsin(olonb);cos_o2inlon4=dcos(olonb);tan_o2inlon4=dtan(olonb)
        deallocate(olonb)
        !P5
        sin_o2inlon5=dsin(o2inlon5);cos_o2inlon5=dcos(o2inlon5)
        !P6
        sin_o2inlon6=dsin(o2inlon6);cos_o2inlon6=dcos(o2inlon6)

        



        
        

    endsubroutine

endmodule
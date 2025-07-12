module md_init
    use md_params
    implicit none
    
    contains

    subroutine init()
        use md_io, only: load_vec,load_cindex,load_b2overlap,load_overlaparea&
        ,load_cidim_and_c2in,load_over2in,load_b2R
        use md_dmatrix
        use md_vector,only:allocate_base,allocate_invec,prepare_sincos_stensil
        implicit none
        
        call allocate_base()
        call load_vec()
        call load_cindex()
        call load_overlaparea()
        call load_b2overlap()
        call load_cidim_and_c2in()
        call load_over2in()
        call load_b2R()
        call make_inputvector()
        call allocate_invec()
        call prepare_sincos_stensil()
        call load_dmatrix()


    endsubroutine

    subroutine init_cube()
        use md_dmatrix
        use md_vector
        use md_cube
        implicit none
        integer k
        double precision l
        double precision,dimension(o2in_size1)::o2in_h1,usb1,vsb1
        double precision,dimension(o2in_size2)::o2in_h2,usb2,vsb2
        double precision,dimension(o2in_size3)::o2in_h3,usb3,vsb3
        double precision,dimension(o2in_size4)::o2in_h4,usb4,vsb4
        double precision,dimension(o2in_size5)::o2in_h5,usb5,vsb5
        double precision,dimension(o2in_size6)::o2in_h6,usb6,vsb6

        hi1=phib(cindex1);hi2=phib(cindex2)
        hi3=phib(cindex3);hi4=phib(cindex4)
        hi5=phib(cindex5);hi6=phib(cindex6)

        us1=ub(cindex1);us2=ub(cindex2)
        us3=ub(cindex3);us4=ub(cindex4)
        us5=ub(cindex5);us6=ub(cindex6)
        
        vs1=vb(cindex1);vs2=vb(cindex2)
        vs3=vb(cindex3);vs4=vb(cindex4)
        vs5=vb(cindex5);vs6=vb(cindex6)

        us=ub(overlaparea)
        vs=vb(overlaparea)
        ho=phib(overlaparea)

        o2in_h1=ho(over2inR1)
        o2in_h2=ho(over2inR2)
        o2in_h3=ho(over2inR3)
        o2in_h4=ho(over2inR4)
        o2in_h5=ho(over2inR5)
        o2in_h6=ho(over2inR6)

        hi1(over2inL1)=ho(over2inR1)
        hi2(over2inL2)=ho(over2inR2)
        hi3(over2inL3)=ho(over2inR3)
        hi4(over2inL4)=ho(over2inR4)
        hi5(over2inL5)=ho(over2inR5)
        hi6(over2inL6)=ho(over2inR6)

        usb1=us(over2inR1)
        usb2=us(over2inR2)
        usb3=us(over2inR3)
        usb4=us(over2inR4)
        usb5=us(over2inR5)
        usb6=us(over2inR6)

        vsb1=vs(over2inR1)
        vsb2=vs(over2inR2)
        vsb3=vs(over2inR3)
        vsb4=vs(over2inR4)
        vsb5=vs(over2inR5)
        vsb6=vs(over2inR6)

        
        call prepare_sincos_overlap()

        call ca2cs_eqv2(pvec_size1,ui1,vi1,us1,vs1,cslat1,cslon1,tslat1,tslon1)
        call ca2cs_eqv2(pvec_size2,ui2,vi2,us2,vs2,cslat2,cslon2,tslat2,tslon2)
        call ca2cs_eqv2(pvec_size3,ui3,vi3,us3,vs3,cslat3,cslon3,tslat3,tslon3)
        call ca2cs_eqv2(pvec_size4,ui4,vi4,us4,vs4,cslat4,cslon4,tslat4,tslon4)
        call ca2cs_np(pvec_size5,ui5,vi5,us5,vs5,cslon5,sslat5,sslon5)
        call ca2cs_sp(pvec_size6,ui6,vi6,us6,vs6,cslon6,sslat6,sslon6)

        call calc_covariant_eqv2(pvec_size1,cui1,cvi1,ui1,vi1,cslat1,cslon1,sslat1,sslon1)
        call calc_covariant_eqv2(pvec_size2,cui2,cvi2,ui2,vi2,cslat2,cslon2,sslat2,sslon2)
        call calc_covariant_eqv2(pvec_size3,cui3,cvi3,ui3,vi3,cslat3,cslon3,sslat3,sslon3)
        call calc_covariant_eqv2(pvec_size4,cui4,cvi4,ui4,vi4,cslat4,cslon4,sslat4,sslon4)
        call calc_covariant_np(pvec_size5,ui5,vi5,cui5,cvi5,cslat5,cslon5,sslat5,sslon5)
        call calc_covariant_sp(pvec_size6,ui6,vi6,cui6,cvi6,cslat6,cslon6,sslat6,sslon6)
        
        
    endsubroutine


endmodule
module md_integration
    use md_vector
    use md_params
    implicit none
   
    
    contains

    subroutine RK4()
        use md_vector
        use md_params
        use md_cube
        use md_dmatrix
        use md_tendency
        use md_io,only: output_in_integration
        implicit none
        
        double precision,dimension(pvec_size1)::h1c,cu1c,cv1c,cut1,cvt1,ht1
        double precision,dimension(pvec_size2)::h2c,cu2c,cv2c,cut2,cvt2,ht2
        double precision,dimension(pvec_size3)::h3c,cu3c,cv3c,cut3,cvt3,ht3
        double precision,dimension(pvec_size4)::h4c,cu4c,cv4c,cut4,cvt4,ht4
        double precision,dimension(pvec_size5)::h5c,cu5c,cv5c,cut5,cvt5,ht5
        double precision,dimension(pvec_size6)::h6c,cu6c,cv6c,cut6,cvt6,ht6

        double precision,dimension(o2in_size1)::cu1b,cv1b,usb1,vsb1,ui1b,vi1b
        double precision,dimension(o2in_size2)::cu2b,cv2b,usb2,vsb2,ui2b,vi2b
        double precision,dimension(o2in_size3)::cu3b,cv3b,usb3,vsb3,ui3b,vi3b
        double precision,dimension(o2in_size4)::cu4b,cv4b,usb4,vsb4,ui4b,vi4b
        double precision,dimension(o2in_size5)::cu5b,cv5b,usb5,vsb5,ui5b,vi5b
        double precision,dimension(o2in_size6)::cu6b,cv6b,usb6,vsb6,ui6b,vi6b
        
        double precision,dimension(c2o_size1)::u1bc,v1bc,u1b,v1b
        double precision,dimension(c2o_size2)::u2bc,v2bc,u2b,v2b
        double precision,dimension(c2o_size3)::u3bc,v3bc,u3b,v3b
        double precision,dimension(c2o_size4)::u4bc,v4bc,u4b,v4b
        double precision,dimension(c2o_size5)::u5bc,v5bc,u5b,v5b
        double precision,dimension(c2o_size6)::u6bc,v6bc,u6b,v6b
        double precision, dimension(n)::acu,acv,ah,au,av 
        
        integer i,srec0,srec

        srec0=rec
        srec=rec

        print*,'start integration dt=',dt
        
        
        do i=1,steps
            h1c=hi1;h2c=hi2;h3c=hi3;h4c=hi4;h5c=hi5;h6c=hi6
            cu1c=cui1;cu2c=cui2;cu3c=cui3;cu4c=cui4;cu5c=cui5;cu6c=cui6
            cv1c=cvi1;cv2c=cvi2;cv3c=cvi3;cv4c=cvi4;cv5c=cvi5;cv6c=cvi6
        
        !K1
            call co2contra_eqv2(pvec_size1,ui1,vi1,cui1,cvi1,cslat1,cslon1,sslat1,sslon1)
            call co2contra_eqv2(pvec_size2,ui2,vi2,cui2,cvi2,cslat2,cslon2,sslat2,sslon2)
            call co2contra_eqv2(pvec_size3,ui3,vi3,cui3,cvi3,cslat3,cslon3,sslat3,sslon3)
            call co2contra_eqv2(pvec_size4,ui4,vi4,cui4,cvi4,cslat4,cslon4,sslat4,sslon4)
            call co2contra_np(pvec_size5,ui5,vi5,cui5,cvi5,cslat5,cslon5,sslat5,sslon5)
            call co2contra_sp(pvec_size6,ui6,vi6,cui6,cvi6,cslat6,cslon6,sslat6,sslon6)
            

            hi1(over2inL1)=ho(over2inR1);hi2(over2inL2)=ho(over2inR2);hi3(over2inL3)=ho(over2inR3)
            hi4(over2inL4)=ho(over2inR4);hi5(over2inL5)=ho(over2inR5);hi6(over2inL6)=ho(over2inR6)

            usb1=us(over2inR1);usb2=us(over2inR2);usb3=us(over2inR3);usb4=us(over2inR4);usb5=us(over2inR5);usb6=us(over2inR6)
            vsb1=vs(over2inR1);vsb2=vs(over2inR2);vsb3=vs(over2inR3);vsb4=vs(over2inR4);vsb5=vs(over2inR5);vsb6=vs(over2inR6)
        

            call ca2cs_eqv2(o2in_size1,ui1b,vi1b,usb1,vsb1,cos_o2inlat1,cos_o2inlon1,tan_o2inlat1,tan_o2inlon1)
            call ca2cs_eqv2(o2in_size2,ui2b,vi2b,usb2,vsb2,cos_o2inlat2,cos_o2inlon2,tan_o2inlat2,tan_o2inlon2)
            call ca2cs_eqv2(o2in_size3,ui3b,vi3b,usb3,vsb3,cos_o2inlat3,cos_o2inlon3,tan_o2inlat3,tan_o2inlon3)
            call ca2cs_eqv2(o2in_size4,ui4b,vi4b,usb4,vsb4,cos_o2inlat4,cos_o2inlon4,tan_o2inlat4,tan_o2inlon4)
            call ca2cs_np(o2in_size5,ui5b,vi5b,usb5,vsb5,cos_o2inlon5,sin_o2inlat5,sin_o2inlon5)
            call ca2cs_sp(o2in_size6,ui6b,vi6b,usb6,vsb6,cos_o2inlon6,sin_o2inlat6,sin_o2inlon6)
            
            
            
            ui1(over2inL1)=ui1b;vi1(over2inL1)=vi1b;ui2(over2inL2)=ui2b;vi2(over2inL2)=vi2b
            ui3(over2inL3)=ui3b;vi3(over2inL3)=vi3b;ui4(over2inL4)=ui4b;vi4(over2inL4)=vi4b
            ui5(over2inL5)=ui5b;vi5(over2inL5)=vi5b;ui6(over2inL6)=ui6b;vi6(over2inL6)=vi6b

            call calc_covariant_eqv2(o2in_size1,cu1b,cv1b,ui1b,vi1b,cos_o2inlat1,cos_o2inlon1,sin_o2inlat1,sin_o2inlon1)
            call calc_covariant_eqv2(o2in_size2,cu2b,cv2b,ui2b,vi2b,cos_o2inlat2,cos_o2inlon2,sin_o2inlat2,sin_o2inlon2)
            call calc_covariant_eqv2(o2in_size3,cu3b,cv3b,ui3b,vi3b,cos_o2inlat3,cos_o2inlon3,sin_o2inlat3,sin_o2inlon3)
            call calc_covariant_eqv2(o2in_size4,cu4b,cv4b,ui4b,vi4b,cos_o2inlat4,cos_o2inlon4,sin_o2inlat4,sin_o2inlon4)
            call calc_covariant_np(o2in_size5,ui5b,vi5b,cu5b,cv5b,cos_o2inlat5,cos_o2inlon5,sin_o2inlat5,sin_o2inlon5)
            call calc_covariant_sp(o2in_size6,ui6b,vi6b,cu6b,cv6b,cos_o2inlat6,cos_o2inlon6,sin_o2inlat6,sin_o2inlon6)
            

            cui1(over2inL1)=cu1b;cvi1(over2inL1)=cv1b
            cui2(over2inL2)=cu2b;cvi2(over2inL2)=cv2b
            cui3(over2inL3)=cu3b;cvi3(over2inL3)=cv3b
            cui4(over2inL4)=cu4b;cvi4(over2inL4)=cv4b
            cui5(over2inL5)=cu5b;cvi5(over2inL5)=cv5b
            cui6(over2inL6)=cu6b;cvi6(over2inL6)=cv6b  
        
            call tendency(pvec_size1,Dx1,Dy1,hiv1,invec1,dxa1,cui1,cvi1,ui1,vi1,hi1,col1,du1,dv1,dh1,g1)
            call tendency(pvec_size2,Dx2,Dy2,hiv2,invec2,dxa2,cui2,cvi2,ui2,vi2,hi2,col2,du2,dv2,dh2,g2)
            call tendency(pvec_size3,Dx3,Dy3,hiv3,invec3,dxa3,cui3,cvi3,ui3,vi3,hi3,col3,du3,dv3,dh3,g3)
            call tendency(pvec_size4,Dx4,Dy4,hiv4,invec4,dxa4,cui4,cvi4,ui4,vi4,hi4,col4,du4,dv4,dh4,g4)
            call tendency(pvec_size5,Dx5,Dy5,hiv5,invec5,dxa5,cui5,cvi5,ui5,vi5,hi5,col5,du5,dv5,dh5,g5)
            call tendency(pvec_size6,Dx6,Dy6,hiv6,invec6,dxa6,cui6,cvi6,ui6,vi6,hi6,col6,du6,dv6,dh6,g6)

            cui1=cui1+(du1/6d0);cvi1=cvi1+(dv1/6d0);hi1=hi1+(dh1/6d0)
            cui2=cui2+(du2/6d0);cvi2=cvi2+(dv2/6d0);hi2=hi2+(dh2/6d0)
            cui3=cui3+(du3/6d0);cvi3=cvi3+(dv3/6d0);hi3=hi3+(dh3/6d0)
            cui4=cui4+(du4/6d0);cvi4=cvi4+(dv4/6d0);hi4=hi4+(dh4/6d0)
            cui5=cui5+(du5/6d0);cvi5=cvi5+(dv5/6d0);hi5=hi5+(dh5/6d0)
            cui6=cui6+(du6/6d0);cvi6=cvi6+(dv6/6d0);hi6=hi6+(dh6/6d0)

        !K2
            cut1=cu1c+(du1/2d0);cvt1=cv1c+(dv1/2d0);ht1=h1c+(dh1/2d0)
            cut2=cu2c+(du2/2d0);cvt2=cv2c+(dv2/2d0);ht2=h2c+(dh2/2d0)
            cut3=cu3c+(du3/2d0);cvt3=cv3c+(dv3/2d0);ht3=h3c+(dh3/2d0)
            cut4=cu4c+(du4/2d0);cvt4=cv4c+(dv4/2d0);ht4=h4c+(dh4/2d0)
            cut5=cu5c+(du5/2d0);cvt5=cv5c+(dv5/2d0);ht5=h5c+(dh5/2d0)
            cut6=cu6c+(du6/2d0);cvt6=cv6c+(dv6/2d0);ht6=h6c+(dh6/2d0)

            ho(b2overlap1)=ht1(in2bR1);ho(b2overlap2)=ht2(in2bR2);ho(b2overlap3)=ht3(in2bR3)
            ho(b2overlap4)=ht4(in2bR4);ho(b2overlap5)=ht5(in2bR5);ho(b2overlap6)=ht6(in2bR6)

            call co2contra_eqv2(pvec_size1,ui1,vi1,cut1,cvt1,cslat1,cslon1,sslat1,sslon1)
            call co2contra_eqv2(pvec_size2,ui2,vi2,cut2,cvt2,cslat2,cslon2,sslat2,sslon2)
            call co2contra_eqv2(pvec_size3,ui3,vi3,cut3,cvt3,cslat3,cslon3,sslat3,sslon3)
            call co2contra_eqv2(pvec_size4,ui4,vi4,cut4,cvt4,cslat4,cslon4,sslat4,sslon4)
            call co2contra_np(pvec_size5,ui5,vi5,cut5,cvt5,cslat5,cslon5,sslat5,sslon5)
            call co2contra_sp(pvec_size6,ui6,vi6,cut6,cvt6,cslat6,cslon6,sslat6,sslon6)
            
        
            u1b=ui1(in2bR1);u2b=ui2(in2bR2);u3b=ui3(in2bR3);u4b=ui4(in2bR4);u5b=ui5(in2bR5);u6b=ui6(in2bR6)
            v1b=vi1(in2bR1);v2b=vi2(in2bR2);v3b=vi3(in2bR3);v4b=vi4(in2bR4);v5b=vi5(in2bR5);v6b=vi6(in2bR6)

            call cs2ca_eqv2(c2o_size1,u1b,v1b,u1bc,v1bc,cos_c2overlat1,cos_c2overlon1,sin_c2overlat1,sin_c2overlon1)
            call cs2ca_eqv2(c2o_size2,u2b,v2b,u2bc,v2bc,cos_c2overlat2,cos_c2overlon2,sin_c2overlat2,sin_c2overlon2)
            call cs2ca_eqv2(c2o_size3,u3b,v3b,u3bc,v3bc,cos_c2overlat3,cos_c2overlon3,sin_c2overlat3,sin_c2overlon3)
            call cs2ca_eqv2(c2o_size4,u4b,v4b,u4bc,v4bc,cos_c2overlat4,cos_c2overlon4,sin_c2overlat4,sin_c2overlon4)
            call cs2ca_np(c2o_size5,u5b,v5b,u5bc,v5bc,cos_c2overlat5,cos_c2overlon5,sin_c2overlat5,sin_c2overlon5)
            call cs2ca_sp(c2o_size6,u6b,v6b,u6bc,v6bc,cos_c2overlat6,cos_c2overlon6,sin_c2overlat6,sin_c2overlon6)

            us(b2overlap1)=u1bc;us(b2overlap2)=u2bc;us(b2overlap3)=u3bc
            us(b2overlap4)=u4bc;us(b2overlap5)=u5bc;us(b2overlap6)=u6bc
            vs(b2overlap1)=v1bc;vs(b2overlap2)=v2bc;vs(b2overlap3)=v3bc
            vs(b2overlap4)=v4bc;vs(b2overlap5)=v5bc;vs(b2overlap6)=v6bc

            ht1(over2inL1)=ho(over2inR1);ht2(over2inL2)=ho(over2inR2);ht3(over2inL3)=ho(over2inR3)
            ht4(over2inL4)=ho(over2inR4);ht5(over2inL5)=ho(over2inR5);ht6(over2inL6)=ho(over2inR6)

            usb1=us(over2inR1);usb2=us(over2inR2);usb3=us(over2inR3);usb4=us(over2inR4);usb5=us(over2inR5);usb6=us(over2inR6)
            vsb1=vs(over2inR1);vsb2=vs(over2inR2);vsb3=vs(over2inR3);vsb4=vs(over2inR4);vsb5=vs(over2inR5);vsb6=vs(over2inR6)
            call ca2cs_eqv2(o2in_size1,ui1b,vi1b,usb1,vsb1,cos_o2inlat1,cos_o2inlon1,tan_o2inlat1,tan_o2inlon1)
            call ca2cs_eqv2(o2in_size2,ui2b,vi2b,usb2,vsb2,cos_o2inlat2,cos_o2inlon2,tan_o2inlat2,tan_o2inlon2)
            call ca2cs_eqv2(o2in_size3,ui3b,vi3b,usb3,vsb3,cos_o2inlat3,cos_o2inlon3,tan_o2inlat3,tan_o2inlon3)
            call ca2cs_eqv2(o2in_size4,ui4b,vi4b,usb4,vsb4,cos_o2inlat4,cos_o2inlon4,tan_o2inlat4,tan_o2inlon4)
            call ca2cs_np(o2in_size5,ui5b,vi5b,usb5,vsb5,cos_o2inlon5,sin_o2inlat5,sin_o2inlon5)
            call ca2cs_sp(o2in_size6,ui6b,vi6b,usb6,vsb6,cos_o2inlon6,sin_o2inlat6,sin_o2inlon6)

            ui1(over2inL1)=ui1b;vi1(over2inL1)=vi1b;ui2(over2inL2)=ui2b;vi2(over2inL2)=vi2b
            ui3(over2inL3)=ui3b;vi3(over2inL3)=vi3b;ui4(over2inL4)=ui4b;vi4(over2inL4)=vi4b
            ui5(over2inL5)=ui5b;vi5(over2inL5)=vi5b;ui6(over2inL6)=ui6b;vi6(over2inL6)=vi6b

            call calc_covariant_eqv2(o2in_size1,cu1b,cv1b,ui1b,vi1b,cos_o2inlat1,cos_o2inlon1,sin_o2inlat1,sin_o2inlon1)
            call calc_covariant_eqv2(o2in_size2,cu2b,cv2b,ui2b,vi2b,cos_o2inlat2,cos_o2inlon2,sin_o2inlat2,sin_o2inlon2)
            call calc_covariant_eqv2(o2in_size3,cu3b,cv3b,ui3b,vi3b,cos_o2inlat3,cos_o2inlon3,sin_o2inlat3,sin_o2inlon3)
            call calc_covariant_eqv2(o2in_size4,cu4b,cv4b,ui4b,vi4b,cos_o2inlat4,cos_o2inlon4,sin_o2inlat4,sin_o2inlon4)
            call calc_covariant_np(o2in_size5,ui5b,vi5b,cu5b,cv5b,cos_o2inlat5,cos_o2inlon5,sin_o2inlat5,sin_o2inlon5)
            call calc_covariant_sp(o2in_size6,ui6b,vi6b,cu6b,cv6b,cos_o2inlat6,cos_o2inlon6,sin_o2inlat6,sin_o2inlon6)

            cut1(over2inL1)=cu1b;cvt1(over2inL1)=cv1b
            cut2(over2inL2)=cu2b;cvt2(over2inL2)=cv2b
            cut3(over2inL3)=cu3b;cvt3(over2inL3)=cv3b
            cut4(over2inL4)=cu4b;cvt4(over2inL4)=cv4b
            cut5(over2inL5)=cu5b;cvt5(over2inL5)=cv5b
            cut6(over2inL6)=cu6b;cvt6(over2inL6)=cv6b
        
            call tendency(pvec_size1,Dx1,Dy1,hiv1,invec1,dxa1,cut1,cvt1,ui1,vi1,ht1,col1,du1,dv1,dh1,g1)
            call tendency(pvec_size2,Dx2,Dy2,hiv2,invec2,dxa2,cut2,cvt2,ui2,vi2,ht2,col2,du2,dv2,dh2,g2)
            call tendency(pvec_size3,Dx3,Dy3,hiv3,invec3,dxa3,cut3,cvt3,ui3,vi3,ht3,col3,du3,dv3,dh3,g3)
            call tendency(pvec_size4,Dx4,Dy4,hiv4,invec4,dxa4,cut4,cvt4,ui4,vi4,ht4,col4,du4,dv4,dh4,g4)
            call tendency(pvec_size5,Dx5,Dy5,hiv5,invec5,dxa5,cut5,cvt5,ui5,vi5,ht5,col5,du5,dv5,dh5,g5)
            call tendency(pvec_size6,Dx6,Dy6,hiv6,invec6,dxa6,cut6,cvt6,ui6,vi6,ht6,col6,du6,dv6,dh6,g6)

        
            cui1=cui1+(du1/3d0);cvi1=cvi1+(dv1/3d0);hi1=hi1+(dh1/3d0)
            cui2=cui2+(du2/3d0);cvi2=cvi2+(dv2/3d0);hi2=hi2+(dh2/3d0)
            cui3=cui3+(du3/3d0);cvi3=cvi3+(dv3/3d0);hi3=hi3+(dh3/3d0)
            cui4=cui4+(du4/3d0);cvi4=cvi4+(dv4/3d0);hi4=hi4+(dh4/3d0)
            cui5=cui5+(du5/3d0);cvi5=cvi5+(dv5/3d0);hi5=hi5+(dh5/3d0)
            cui6=cui6+(du6/3d0);cvi6=cvi6+(dv6/3d0);hi6=hi6+(dh6/3d0)

        !k3

            cut1=cu1c+(du1/2d0);cvt1=cv1c+(dv1/2d0);ht1=h1c+(dh1/2d0)
            cut2=cu2c+(du2/2d0);cvt2=cv2c+(dv2/2d0);ht2=h2c+(dh2/2d0)
            cut3=cu3c+(du3/2d0);cvt3=cv3c+(dv3/2d0);ht3=h3c+(dh3/2d0)
            cut4=cu4c+(du4/2d0);cvt4=cv4c+(dv4/2d0);ht4=h4c+(dh4/2d0)
            cut5=cu5c+(du5/2d0);cvt5=cv5c+(dv5/2d0);ht5=h5c+(dh5/2d0)
            cut6=cu6c+(du6/2d0);cvt6=cv6c+(dv6/2d0);ht6=h6c+(dh6/2d0)

   
            ho(b2overlap1)=ht1(in2bR1);ho(b2overlap2)=ht2(in2bR2);ho(b2overlap3)=ht3(in2bR3)
            ho(b2overlap4)=ht4(in2bR4);ho(b2overlap5)=ht5(in2bR5);ho(b2overlap6)=ht6(in2bR6)

            call co2contra_eqv2(pvec_size1,ui1,vi1,cut1,cvt1,cslat1,cslon1,sslat1,sslon1)
            call co2contra_eqv2(pvec_size2,ui2,vi2,cut2,cvt2,cslat2,cslon2,sslat2,sslon2)
            call co2contra_eqv2(pvec_size3,ui3,vi3,cut3,cvt3,cslat3,cslon3,sslat3,sslon3)
            call co2contra_eqv2(pvec_size4,ui4,vi4,cut4,cvt4,cslat4,cslon4,sslat4,sslon4)
            call co2contra_np(pvec_size5,ui5,vi5,cut5,cvt5,cslat5,cslon5,sslat5,sslon5)
            call co2contra_sp(pvec_size6,ui6,vi6,cut6,cvt6,cslat6,cslon6,sslat6,sslon6)

            
            u1b=ui1(in2bR1);u2b=ui2(in2bR2);u3b=ui3(in2bR3);u4b=ui4(in2bR4);u5b=ui5(in2bR5);u6b=ui6(in2bR6)
            v1b=vi1(in2bR1);v2b=vi2(in2bR2);v3b=vi3(in2bR3);v4b=vi4(in2bR4);v5b=vi5(in2bR5);v6b=vi6(in2bR6)

            call cs2ca_eqv2(c2o_size1,u1b,v1b,u1bc,v1bc,cos_c2overlat1,cos_c2overlon1,sin_c2overlat1,sin_c2overlon1)
            call cs2ca_eqv2(c2o_size2,u2b,v2b,u2bc,v2bc,cos_c2overlat2,cos_c2overlon2,sin_c2overlat2,sin_c2overlon2)
            call cs2ca_eqv2(c2o_size3,u3b,v3b,u3bc,v3bc,cos_c2overlat3,cos_c2overlon3,sin_c2overlat3,sin_c2overlon3)
            call cs2ca_eqv2(c2o_size4,u4b,v4b,u4bc,v4bc,cos_c2overlat4,cos_c2overlon4,sin_c2overlat4,sin_c2overlon4)
            call cs2ca_np(c2o_size5,u5b,v5b,u5bc,v5bc,cos_c2overlat5,cos_c2overlon5,sin_c2overlat5,sin_c2overlon5)
            call cs2ca_sp(c2o_size6,u6b,v6b,u6bc,v6bc,cos_c2overlat6,cos_c2overlon6,sin_c2overlat6,sin_c2overlon6)


            us(b2overlap1)=u1bc;us(b2overlap2)=u2bc;us(b2overlap3)=u3bc
            us(b2overlap4)=u4bc;us(b2overlap5)=u5bc;us(b2overlap6)=u6bc
            vs(b2overlap1)=v1bc;vs(b2overlap2)=v2bc;vs(b2overlap3)=v3bc
            vs(b2overlap4)=v4bc;vs(b2overlap5)=v5bc;vs(b2overlap6)=v6bc


            ht1(over2inL1)=ho(over2inR1);ht2(over2inL2)=ho(over2inR2);ht3(over2inL3)=ho(over2inR3)
            ht4(over2inL4)=ho(over2inR4);ht5(over2inL5)=ho(over2inR5);ht6(over2inL6)=ho(over2inR6)

            usb1=us(over2inR1);usb2=us(over2inR2);usb3=us(over2inR3);usb4=us(over2inR4);usb5=us(over2inR5);usb6=us(over2inR6)
            vsb1=vs(over2inR1);vsb2=vs(over2inR2);vsb3=vs(over2inR3);vsb4=vs(over2inR4);vsb5=vs(over2inR5);vsb6=vs(over2inR6)
            call ca2cs_eqv2(o2in_size1,ui1b,vi1b,usb1,vsb1,cos_o2inlat1,cos_o2inlon1,tan_o2inlat1,tan_o2inlon1)
            call ca2cs_eqv2(o2in_size2,ui2b,vi2b,usb2,vsb2,cos_o2inlat2,cos_o2inlon2,tan_o2inlat2,tan_o2inlon2)
            call ca2cs_eqv2(o2in_size3,ui3b,vi3b,usb3,vsb3,cos_o2inlat3,cos_o2inlon3,tan_o2inlat3,tan_o2inlon3)
            call ca2cs_eqv2(o2in_size4,ui4b,vi4b,usb4,vsb4,cos_o2inlat4,cos_o2inlon4,tan_o2inlat4,tan_o2inlon4)
            call ca2cs_np(o2in_size5,ui5b,vi5b,usb5,vsb5,cos_o2inlon5,sin_o2inlat5,sin_o2inlon5)
            call ca2cs_sp(o2in_size6,ui6b,vi6b,usb6,vsb6,cos_o2inlon6,sin_o2inlat6,sin_o2inlon6)

            ui1(over2inL1)=ui1b;vi1(over2inL1)=vi1b;ui2(over2inL2)=ui2b;vi2(over2inL2)=vi2b
            ui3(over2inL3)=ui3b;vi3(over2inL3)=vi3b;ui4(over2inL4)=ui4b;vi4(over2inL4)=vi4b
            ui5(over2inL5)=ui5b;vi5(over2inL5)=vi5b;ui6(over2inL6)=ui6b;vi6(over2inL6)=vi6b

            call calc_covariant_eqv2(o2in_size1,cu1b,cv1b,ui1b,vi1b,cos_o2inlat1,cos_o2inlon1,sin_o2inlat1,sin_o2inlon1)
            call calc_covariant_eqv2(o2in_size2,cu2b,cv2b,ui2b,vi2b,cos_o2inlat2,cos_o2inlon2,sin_o2inlat2,sin_o2inlon2)
            call calc_covariant_eqv2(o2in_size3,cu3b,cv3b,ui3b,vi3b,cos_o2inlat3,cos_o2inlon3,sin_o2inlat3,sin_o2inlon3)
            call calc_covariant_eqv2(o2in_size4,cu4b,cv4b,ui4b,vi4b,cos_o2inlat4,cos_o2inlon4,sin_o2inlat4,sin_o2inlon4)
            call calc_covariant_np(o2in_size5,ui5b,vi5b,cu5b,cv5b,cos_o2inlat5,cos_o2inlon5,sin_o2inlat5,sin_o2inlon5)
            call calc_covariant_sp(o2in_size6,ui6b,vi6b,cu6b,cv6b,cos_o2inlat6,cos_o2inlon6,sin_o2inlat6,sin_o2inlon6)

            cut1(over2inL1)=cu1b;cvt1(over2inL1)=cv1b
            cut2(over2inL2)=cu2b;cvt2(over2inL2)=cv2b
            cut3(over2inL3)=cu3b;cvt3(over2inL3)=cv3b
            cut4(over2inL4)=cu4b;cvt4(over2inL4)=cv4b
            cut5(over2inL5)=cu5b;cvt5(over2inL5)=cv5b
            cut6(over2inL6)=cu6b;cvt6(over2inL6)=cv6b
        
            call tendency(pvec_size1,Dx1,Dy1,hiv1,invec1,dxa1,cut1,cvt1,ui1,vi1,ht1,col1,du1,dv1,dh1,g1)
            call tendency(pvec_size2,Dx2,Dy2,hiv2,invec2,dxa2,cut2,cvt2,ui2,vi2,ht2,col2,du2,dv2,dh2,g2)
            call tendency(pvec_size3,Dx3,Dy3,hiv3,invec3,dxa3,cut3,cvt3,ui3,vi3,ht3,col3,du3,dv3,dh3,g3)
            call tendency(pvec_size4,Dx4,Dy4,hiv4,invec4,dxa4,cut4,cvt4,ui4,vi4,ht4,col4,du4,dv4,dh4,g4)
            call tendency(pvec_size5,Dx5,Dy5,hiv5,invec5,dxa5,cut5,cvt5,ui5,vi5,ht5,col5,du5,dv5,dh5,g5)
            call tendency(pvec_size6,Dx6,Dy6,hiv6,invec6,dxa6,cut6,cvt6,ui6,vi6,ht6,col6,du6,dv6,dh6,g6)
        
            cui1=cui1+(du1/3d0);cvi1=cvi1+(dv1/3d0);hi1=hi1+(dh1/3d0)
            cui2=cui2+(du2/3d0);cvi2=cvi2+(dv2/3d0);hi2=hi2+(dh2/3d0)
            cui3=cui3+(du3/3d0);cvi3=cvi3+(dv3/3d0);hi3=hi3+(dh3/3d0)
            cui4=cui4+(du4/3d0);cvi4=cvi4+(dv4/3d0);hi4=hi4+(dh4/3d0)
            cui5=cui5+(du5/3d0);cvi5=cvi5+(dv5/3d0);hi5=hi5+(dh5/3d0)
            cui6=cui6+(du6/3d0);cvi6=cvi6+(dv6/3d0);hi6=hi6+(dh6/3d0)

        !k4

            cut1=cu1c+du1;cvt1=cv1c+dv1;ht1=h1c+dh1
            cut2=cu2c+du2;cvt2=cv2c+dv2;ht2=h2c+dh2
            cut3=cu3c+du3;cvt3=cv3c+dv3;ht3=h3c+dh3
            cut4=cu4c+du4;cvt4=cv4c+dv4;ht4=h4c+dh4
            cut5=cu5c+du5;cvt5=cv5c+dv5;ht5=h5c+dh5
            cut6=cu6c+du6;cvt6=cv6c+dv6;ht6=h6c+dh6


            ho(b2overlap1)=ht1(in2bR1);ho(b2overlap2)=ht2(in2bR2);ho(b2overlap3)=ht3(in2bR3)
            ho(b2overlap4)=ht4(in2bR4);ho(b2overlap5)=ht5(in2bR5);ho(b2overlap6)=ht6(in2bR6)

            call co2contra_eqv2(pvec_size1,ui1,vi1,cut1,cvt1,cslat1,cslon1,sslat1,sslon1)
            call co2contra_eqv2(pvec_size2,ui2,vi2,cut2,cvt2,cslat2,cslon2,sslat2,sslon2)
            call co2contra_eqv2(pvec_size3,ui3,vi3,cut3,cvt3,cslat3,cslon3,sslat3,sslon3)
            call co2contra_eqv2(pvec_size4,ui4,vi4,cut4,cvt4,cslat4,cslon4,sslat4,sslon4)
            call co2contra_np(pvec_size5,ui5,vi5,cut5,cvt5,cslat5,cslon5,sslat5,sslon5)
            call co2contra_sp(pvec_size6,ui6,vi6,cut6,cvt6,cslat6,cslon6,sslat6,sslon6)

            u1b=ui1(in2bR1);u2b=ui2(in2bR2);u3b=ui3(in2bR3);u4b=ui4(in2bR4);u5b=ui5(in2bR5);u6b=ui6(in2bR6)
            v1b=vi1(in2bR1);v2b=vi2(in2bR2);v3b=vi3(in2bR3);v4b=vi4(in2bR4);v5b=vi5(in2bR5);v6b=vi6(in2bR6)

            call cs2ca_eqv2(c2o_size1,u1b,v1b,u1bc,v1bc,cos_c2overlat1,cos_c2overlon1,sin_c2overlat1,sin_c2overlon1)
            call cs2ca_eqv2(c2o_size2,u2b,v2b,u2bc,v2bc,cos_c2overlat2,cos_c2overlon2,sin_c2overlat2,sin_c2overlon2)
            call cs2ca_eqv2(c2o_size3,u3b,v3b,u3bc,v3bc,cos_c2overlat3,cos_c2overlon3,sin_c2overlat3,sin_c2overlon3)
            call cs2ca_eqv2(c2o_size4,u4b,v4b,u4bc,v4bc,cos_c2overlat4,cos_c2overlon4,sin_c2overlat4,sin_c2overlon4)
            call cs2ca_np(c2o_size5,u5b,v5b,u5bc,v5bc,cos_c2overlat5,cos_c2overlon5,sin_c2overlat5,sin_c2overlon5)
            call cs2ca_sp(c2o_size6,u6b,v6b,u6bc,v6bc,cos_c2overlat6,cos_c2overlon6,sin_c2overlat6,sin_c2overlon6)


            us(b2overlap1)=u1bc;us(b2overlap2)=u2bc;us(b2overlap3)=u3bc
            us(b2overlap4)=u4bc;us(b2overlap5)=u5bc;us(b2overlap6)=u6bc
            vs(b2overlap1)=v1bc;vs(b2overlap2)=v2bc;vs(b2overlap3)=v3bc
            vs(b2overlap4)=v4bc;vs(b2overlap5)=v5bc;vs(b2overlap6)=v6bc

   
            ht1(over2inL1)=ho(over2inR1);ht2(over2inL2)=ho(over2inR2);ht3(over2inL3)=ho(over2inR3)
            ht4(over2inL4)=ho(over2inR4);ht5(over2inL5)=ho(over2inR5);ht6(over2inL6)=ho(over2inR6)

            usb1=us(over2inR1);usb2=us(over2inR2);usb3=us(over2inR3);usb4=us(over2inR4);usb5=us(over2inR5);usb6=us(over2inR6)
            vsb1=vs(over2inR1);vsb2=vs(over2inR2);vsb3=vs(over2inR3);vsb4=vs(over2inR4);vsb5=vs(over2inR5);vsb6=vs(over2inR6)
            call ca2cs_eqv2(o2in_size1,ui1b,vi1b,usb1,vsb1,cos_o2inlat1,cos_o2inlon1,tan_o2inlat1,tan_o2inlon1)
            call ca2cs_eqv2(o2in_size2,ui2b,vi2b,usb2,vsb2,cos_o2inlat2,cos_o2inlon2,tan_o2inlat2,tan_o2inlon2)
            call ca2cs_eqv2(o2in_size3,ui3b,vi3b,usb3,vsb3,cos_o2inlat3,cos_o2inlon3,tan_o2inlat3,tan_o2inlon3)
            call ca2cs_eqv2(o2in_size4,ui4b,vi4b,usb4,vsb4,cos_o2inlat4,cos_o2inlon4,tan_o2inlat4,tan_o2inlon4)
            call ca2cs_np(o2in_size5,ui5b,vi5b,usb5,vsb5,cos_o2inlon5,sin_o2inlat5,sin_o2inlon5)
            call ca2cs_sp(o2in_size6,ui6b,vi6b,usb6,vsb6,cos_o2inlon6,sin_o2inlat6,sin_o2inlon6)

            ui1(over2inL1)=ui1b;vi1(over2inL1)=vi1b;ui2(over2inL2)=ui2b;vi2(over2inL2)=vi2b
            ui3(over2inL3)=ui3b;vi3(over2inL3)=vi3b;ui4(over2inL4)=ui4b;vi4(over2inL4)=vi4b
            ui5(over2inL5)=ui5b;vi5(over2inL5)=vi5b;ui6(over2inL6)=ui6b;vi6(over2inL6)=vi6b

            call calc_covariant_eqv2(o2in_size1,cu1b,cv1b,ui1b,vi1b,cos_o2inlat1,cos_o2inlon1,sin_o2inlat1,sin_o2inlon1)
            call calc_covariant_eqv2(o2in_size2,cu2b,cv2b,ui2b,vi2b,cos_o2inlat2,cos_o2inlon2,sin_o2inlat2,sin_o2inlon2)
            call calc_covariant_eqv2(o2in_size3,cu3b,cv3b,ui3b,vi3b,cos_o2inlat3,cos_o2inlon3,sin_o2inlat3,sin_o2inlon3)
            call calc_covariant_eqv2(o2in_size4,cu4b,cv4b,ui4b,vi4b,cos_o2inlat4,cos_o2inlon4,sin_o2inlat4,sin_o2inlon4)
            call calc_covariant_np(o2in_size5,ui5b,vi5b,cu5b,cv5b,cos_o2inlat5,cos_o2inlon5,sin_o2inlat5,sin_o2inlon5)
            call calc_covariant_sp(o2in_size6,ui6b,vi6b,cu6b,cv6b,cos_o2inlat6,cos_o2inlon6,sin_o2inlat6,sin_o2inlon6)


            cut1(over2inL1)=cu1b;cvt1(over2inL1)=cv1b
            cut2(over2inL2)=cu2b;cvt2(over2inL2)=cv2b
            cut3(over2inL3)=cu3b;cvt3(over2inL3)=cv3b
            cut4(over2inL4)=cu4b;cvt4(over2inL4)=cv4b
            cut5(over2inL5)=cu5b;cvt5(over2inL5)=cv5b
            cut6(over2inL6)=cu6b;cvt6(over2inL6)=cv6b
        
            call tendency(pvec_size1,Dx1,Dy1,hiv1,invec1,dxa1,cut1,cvt1,ui1,vi1,ht1,col1,du1,dv1,dh1,g1)
            call tendency(pvec_size2,Dx2,Dy2,hiv2,invec2,dxa2,cut2,cvt2,ui2,vi2,ht2,col2,du2,dv2,dh2,g2)
            call tendency(pvec_size3,Dx3,Dy3,hiv3,invec3,dxa3,cut3,cvt3,ui3,vi3,ht3,col3,du3,dv3,dh3,g3)
            call tendency(pvec_size4,Dx4,Dy4,hiv4,invec4,dxa4,cut4,cvt4,ui4,vi4,ht4,col4,du4,dv4,dh4,g4)
            call tendency(pvec_size5,Dx5,Dy5,hiv5,invec5,dxa5,cut5,cvt5,ui5,vi5,ht5,col5,du5,dv5,dh5,g5)
            call tendency(pvec_size6,Dx6,Dy6,hiv6,invec6,dxa6,cut6,cvt6,ui6,vi6,ht6,col6,du6,dv6,dh6,g6)
        
            cui1=cui1+(du1/6d0);cvi1=cvi1+(dv1/6d0);hi1=hi1+(dh1/6d0)
            cui2=cui2+(du2/6d0);cvi2=cvi2+(dv2/6d0);hi2=hi2+(dh2/6d0)
            cui3=cui3+(du3/6d0);cvi3=cvi3+(dv3/6d0);hi3=hi3+(dh3/6d0)
            cui4=cui4+(du4/6d0);cvi4=cvi4+(dv4/6d0);hi4=hi4+(dh4/6d0)
            cui5=cui5+(du5/6d0);cvi5=cvi5+(dv5/6d0);hi5=hi5+(dh5/6d0)
            cui6=cui6+(du6/6d0);cvi6=cvi6+(dv6/6d0);hi6=hi6+(dh6/6d0)


         
            ho(b2overlap1)=hi1(in2bR1);ho(b2overlap2)=hi2(in2bR2);ho(b2overlap3)=hi3(in2bR3)
            ho(b2overlap4)=hi4(in2bR4);ho(b2overlap5)=hi5(in2bR5);ho(b2overlap6)=hi6(in2bR6)

            call co2contra_eqv2(pvec_size1,ui1,vi1,cui1,cvi1,cslat1,cslon1,sslat1,sslon1)
            call co2contra_eqv2(pvec_size2,ui2,vi2,cui2,cvi2,cslat2,cslon2,sslat2,sslon2)
            call co2contra_eqv2(pvec_size3,ui3,vi3,cui3,cvi3,cslat3,cslon3,sslat3,sslon3)
            call co2contra_eqv2(pvec_size4,ui4,vi4,cui4,cvi4,cslat4,cslon4,sslat4,sslon4)
            call co2contra_np(pvec_size5,ui5,vi5,cui5,cvi5,cslat5,cslon5,sslat5,sslon5)
            call co2contra_sp(pvec_size6,ui6,vi6,cui6,cvi6,cslat6,cslon6,sslat6,sslon6)


            u1b=ui1(in2bR1);u2b=ui2(in2bR2);u3b=ui3(in2bR3);u4b=ui4(in2bR4);u5b=ui5(in2bR5);u6b=ui6(in2bR6)
            v1b=vi1(in2bR1);v2b=vi2(in2bR2);v3b=vi3(in2bR3);v4b=vi4(in2bR4);v5b=vi5(in2bR5);v6b=vi6(in2bR6)

            call cs2ca_eqv2(c2o_size1,u1b,v1b,u1bc,v1bc,cos_c2overlat1,cos_c2overlon1,sin_c2overlat1,sin_c2overlon1)
            call cs2ca_eqv2(c2o_size2,u2b,v2b,u2bc,v2bc,cos_c2overlat2,cos_c2overlon2,sin_c2overlat2,sin_c2overlon2)
            call cs2ca_eqv2(c2o_size3,u3b,v3b,u3bc,v3bc,cos_c2overlat3,cos_c2overlon3,sin_c2overlat3,sin_c2overlon3)
            call cs2ca_eqv2(c2o_size4,u4b,v4b,u4bc,v4bc,cos_c2overlat4,cos_c2overlon4,sin_c2overlat4,sin_c2overlon4)
            call cs2ca_np(c2o_size5,u5b,v5b,u5bc,v5bc,cos_c2overlat5,cos_c2overlon5,sin_c2overlat5,sin_c2overlon5)
            call cs2ca_sp(c2o_size6,u6b,v6b,u6bc,v6bc,cos_c2overlat6,cos_c2overlon6,sin_c2overlat6,sin_c2overlon6)

            us(b2overlap1)=u1bc;us(b2overlap2)=u2bc;us(b2overlap3)=u3bc
            us(b2overlap4)=u4bc;us(b2overlap5)=u5bc;us(b2overlap6)=u6bc
            vs(b2overlap1)=v1bc;vs(b2overlap2)=v2bc;vs(b2overlap3)=v3bc
            vs(b2overlap4)=v4bc;vs(b2overlap5)=v5bc;vs(b2overlap6)=v6bc


            if(i==srec)then
                print*,srec/rec*1,cui6(1)
                call output_in_integration(i)
                srec=srec+rec
            endif
            if(i<500)then
                print*,cui1(1)
            endif
        enddo

    endsubroutine

    

    
    

endmodule
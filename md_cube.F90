module md_cube
    use md_vector
    use md_params
    implicit none
    
    contains


    subroutine co2contra_eqv2(size,u1,v1,cu1,cv1,clat,clon,slat,slon)
        implicit none
        integer size
        real(kind=8) u1(size),v1(size),clon(size),clat(size),slat(size),slon(size),cu1(size),cv1(size) 
        u1=ca**2d0*(((clon**2d0*clat**4d0)*cu1)+((clat**3d0*clon**2d0*slat*slon)*cv1))/(clat**6d0*clon**6d0)
        v1=ca**2d0*(((clat**3d0*clon**2d0*slat*slon)*cu1)+((clat**2d0*clon**2d0*((clon**2d0)+(slat**2d0*slon**2d0)))*cv1))&
        /(clat**6d0*clon**6d0)
    endsubroutine


    subroutine calc_covariant_eqv2(size,cu1,cv1,u1,v1,clat,clon,slat,slon)
        implicit none
        integer size
        real(kind=8) u1(size),v1(size),clon(size),clat(size),slat(size),slon(size),cu1(size),cv1(size)
        cu1=(1d0/ca**2d0)*((((clat**2d0*clon**2d0)*(clon**2d0+(slat**2d0*slon**2d0)))*u1)-((clat**3d0*clon**2d0*slat*slon)*v1))
        cv1=(1d0/ca**2d0)*((-clat**3d0*clon**2d0*slat*slon*u1)+(clon**2d0*clat**4d0*v1))
        
    endsubroutine

    subroutine ca2cs_eqv2(size,u1,v1,us1,vs1,clat,clon,tlat,tlon)
        implicit none
        integer size
        real(kind=8) u1(size),v1(size),clon(size),clat(size),tlat(size),tlon(size),us1(size),vs1(size) 
        u1=ca*(us1/(clat*clon**2d0))
        v1=ca*((us1*tlon*tlat)+(vs1/clat))/(clat*clon)
    endsubroutine

    subroutine cs2ca_eqv2(size,u1,v1,us1,vs1,clat,clon,slat,slon)
        implicit none
        integer size
        real(kind=8) u1(size),v1(size),clon(size),clat(size),slat(size),slon(size),us1(size),vs1(size) 
        us1=(u1*clat*(clon**2d0))/ca
        vs1=((-u1*slat*slon*clat*clon)+(v1*clon*clat**2d0))/ca
    endsubroutine

    subroutine ca2cs_np(size,u1,v1,us1,vs1,coslon,sinlat,sinlon)
        implicit none
        integer size
        real(kind=8) u1(size),v1(size),coslon(size),sinlat(size),sinlon(size),us1(size),vs1(size) 
        
        u1=((us1*sinlat*coslon)-(vs1*sinlon))/(sinlat**2d0)*ca
        v1=((us1*sinlat*sinlon)+(vs1*coslon))/(sinlat**2d0)*ca
    endsubroutine

    subroutine cs2ca_np(size,u1,v1,us1,vs1,clat,clon,slat,slon)
        implicit none
        integer size
        real(kind=8) u1(size),v1(size),clon(size),clat(size),slat(size),slon(size),us1(size),vs1(size) 
        us1=((slat)*((u1*clon)+(v1*slon)))/ca
        vs1=((slat)*((-u1*slat*slon)+(v1*slat*clon)))/ca
    endsubroutine

    subroutine calc_covariant_np(size,u1,v1,cu1,cv1,clat,clon,slat,slon)
        implicit none
        integer size
        real(kind=8) u1(size),v1(size),clon(size),clat(size),slat(size),slon(size),cu1(size),cv1(size)  
        cu1=((slat**2d0*((clon**2d0)+(slat**2d0*slon**2d0))*u1)+((clat**2d0*slat**2d0*slon*clon)*v1))/ca**2d0
        cv1=(((clat**2d0*slat**2d0*slon*clon)*u1)+((slat**2d0*((slon**2d0)+(clon**2d0*slat**2d0)))*v1))/ca**2d0
    endsubroutine

    subroutine co2contra_np(size,u1,v1,cu1,cv1,clat,clon,slat,slon)
        implicit none
        integer size
        real(kind=8) u1(size),v1(size),clon(size),clat(size),slat(size),slon(size),cu1(size),cv1(size)
        u1=ca**2d0*(((slat**2d0*((slon**2d0)+(clon**2d0*slat**2d0)))*cu1)&
        -((clat**2d0*slat**2d0*slon*clon)*cv1))/(slat**6d0)
        v1=ca**2d0*(((-clat**2d0*slat**2d0*slon*clon)*cu1)+&
        ((slat**2d0*((clon**2d0)+(slat**2d0*slon**2d0)))*cv1))/(slat**6d0)

    endsubroutine

    subroutine ca2cs_sp(size,u1,v1,us1,vs1,coslon,sinlat,sinlon)
        implicit none
        integer size
        real(kind=8) u1(size),v1(size),coslon(size),sinlat(size),sinlon(size),us1(size),vs1(size) 
        
        u1=((-us1*sinlat*coslon)+(vs1*sinlon))/(sinlat**2d0)*ca
        v1=((us1*sinlat*sinlon)+(vs1*coslon))/(sinlat**2d0)*ca
    endsubroutine

    subroutine cs2ca_sp(size,u1,v1,us1,vs1,clat,clon,slat,slon)
        implicit none
        integer size
        real(kind=8) u1(size),v1(size),clon(size),clat(size),slat(size),slon(size),us1(size),vs1(size) 
        us1=((slat)*((-u1*clon)+(v1*slon)))/ca
        vs1=((slat)*((u1*slat*slon)+(v1*slat*clon)))/ca
    endsubroutine

    subroutine calc_covariant_sp(size,u1,v1,cu1,cv1,clat,clon,slat,slon)
        implicit none
        integer size
        real(kind=8) u1(size),v1(size),clon(size),clat(size),slat(size),slon(size),cu1(size),cv1(size)  
        cu1=((slat**2d0*((clon**2d0)+(slat**2d0*slon**2d0))*u1)-(clat**2d0*slat**2d0*slon*clon*v1))/ca**2d0
        cv1=((-clat**2d0*slat**2d0*slon*clon*u1)+(slat**2d0*((slon**2d0)+(clon**2d0*slat**2d0))*v1))/ca**2d0
    endsubroutine

    subroutine co2contra_sp(size,u1,v1,cu1,cv1,clat,clon,slat,slon)
        implicit none
        integer size
        real(kind=8) u1(size),v1(size),clon(size),clat(size),slat(size),slon(size),cu1(size),cv1(size)
        u1=ca**2d0*(((slat**2d0*((slon**2d0)+(clon**2d0*slat**2d0)))*cu1)&
        +((clat**2d0*slat**2d0*slon*clon)*cv1))/(slat**6d0)
        v1=ca**2d0*(((clat**2d0*slat**2d0*slon*clon)*cu1)&
        +((slat**2d0*((clon**2d0)+(slat**2d0*slon**2d0)))*cv1))/(slat**6d0)

    endsubroutine

endmodule
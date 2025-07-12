module initial_condition
    use md_vector, only:lat,lon,ub,vb,phib,fl,us,vs
    use md_params
    implicit none


    contains

    subroutine case2()
        use md_params
        implicit none
        real(kind=8)u0,gh0
        print*,'alpha=',alpha/pi*180d0
        gh0=2.94d4
        u0=2d0*pi*ea/(12d0*24d0*60d0*60d0)
        ub=u0*(dcos(lat)*dcos(alpha)+dcos(lon)*dsin(lat)*dsin(alpha))
        vb=-u0*dsin(lon)*dsin(alpha)
        phib=gh0-(ea*omega*u0+(u0**2d0)/2d0)*(-dcos(lon)*dcos(lat)*dsin(alpha)+dsin(lat)*dcos(alpha))**2d0
        phib=phib/g
        
    endsubroutine

    subroutine case6()
        use md_params
        implicit none
        real(kind=8) omega6,K,R,h0
        real(kind=8),dimension(n)::Agh,Bgh,Cgh
        omega6=7.848d-6
        K=omega6
        R=4d0
        h0=8d3
        ub=ea*omega6*dcos(lat)+ea*K*dcos(lat)**(R-1d0)*(R*dsin(lat)**2d0-dcos(lat)**2d0)*dcos(R*lon)
        vb=-ea*K*R*dcos(lat)**(R-1d0)*dsin(lat)*dsin(R*lon)

        Agh=omega6*0.5d0*(2d0*omega+omega6)*dcos(lat)**2d0&
            +0.25d0*K**2d0*dcos(lat)**(2d0*R)&
                *((R+1d0)*dcos(lat)**2d0+((2d0*R**2d0)-R-2d0)-2d0*R**2d0/dcos(lat)**(2d0))
        Bgh=2d0*(omega+omega6)*K/((R+1d0)*(R+2d0))&
            *dcos(lat)**(R)&
                *((R**2d0+2d0*R+2d0)-(R+1d0)**2d0*dcos(lat)**2d0)
        Cgh=0.25d0*K**2d0*dcos(lat)**(2d0*R)*&
            ((R+1d0)*dcos(lat)**2d0-(R+2d0))

        phib=(g*h0+ea**2d0*Agh+ea**2d0*Bgh*dcos(R*lon)+ea**2d0*Cgh*dcos(2d0*R*lon))/g

    endsubroutine


endmodule
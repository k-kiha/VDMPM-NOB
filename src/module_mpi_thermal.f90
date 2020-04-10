module mpi_thermal
    ! use debug
    use global
    implicit none
    integer :: m1sb,m2sb,m3sb
    double precision, allocatable, dimension(:,:,:) :: T    
    double precision, allocatable, dimension(:,:) :: TBCbt_sub,TBCup_sub
    double precision, allocatable, dimension(:,:,:) :: KPP,dKPP,invRhoCp,dinvRhoCp
contains    
    subroutine mpi_thermal_allocation(n1in,n2in,n3in)
        implicit none
        integer :: n1in,n2in,n3in

        m1sb=n1in
        m2sb=n2in
        m3sb=n3in

        allocate(T(0:m1sb,0:m2sb,0:m3sb))
        allocate(TBCup_sub(0:m1sb,0:m3sb),TBCbt_sub(0:m1sb,0:m3sb))

        T(:,:,:) = 0.d0
        TBCup_sub(:,:) = 0.d0
        TBCbt_sub(:,:) = 0.d0

    end subroutine mpi_thermal_allocation

    subroutine mpi_thermal_clean()
        implicit none
        
        deallocate(T)
        deallocate(TBCup_sub,TBCbt_sub)
        
    end subroutine mpi_thermal_clean

    subroutine mpi_thermal_initial(x1,x2,x3)
        implicit none
        double precision:: x1(0:m1sb),x2(0:m2sb),x3(0:m3sb)
        integer :: i,j,k
        
        T(0:m1sb,0:m2sb,0:m3sb)=0.d0

        ! do k=1,m3sb-1
        !     do j=1,m2sb-1
        !         ! do i=1,m1sb-1
        !             T(:,j,k)=(Tup-Tbt)/L2*(x2(j+1)+x2(j))*0.5d0+Tbt
        !         ! enddo
        !     enddo
        ! enddo
    end subroutine mpi_thermal_initial

    subroutine mpi_thermal_boundary(x1,x2,x3, comm_1d_x1, comm_1d_x2, comm_1d_x3)
        use mpi_topology, only : cart_comm_1d
        implicit none
        type(cart_comm_1d), intent(in)  :: comm_1d_x1, comm_1d_x2, comm_1d_x3
        double precision, dimension(0:m1sb, 0:m2sb, 0:m3sb) :: x1,x2,x3
        integer :: k
        
        TBCup_sub(:,:)=Tup
        TBCbt_sub(:,:)=Tbt
        
        if(comm_1d_x2%myrank==0) then   !No i,j Dirichlet 
            do k=0, m3sb
                T(0:m1sb,0,k)=TBCbt_sub(:,k)
            enddo
        else if(comm_1d_x2%myrank==comm_1d_x2%nprocs-1) then
            do k=0, m3sb
                T(0:m1sb,m2sb,k)=TBCup_sub(:,k)
            enddo
        endif    
    end subroutine mpi_thermal_boundary

    subroutine mpi_thermal_coeffi()
        implicit none
        double precision, allocatable, dimension(:) :: tmp
        integer :: i,j,k
        allocate(KPP(0:m1sb,0:m2sb,0:m3sb),dKPP(0:m1sb,0:m2sb,0:m3sb))
        allocate(invRhoCp(0:m1sb,0:m2sb,0:m3sb),dinvRhoCp(0:m1sb,0:m2sb,0:m3sb))
        allocate(tmp(0:m1sb))
        
        do k=0,m3sb
            do j=0,m2sb
                tmp(:)=DeltaT*T(:,j,k)
                KPP(:,j,k)  = c10*(1.d0 + c11*tmp(:) +      c12*tmp(:)**2.d0 +      c13*tmp(:)**3.d0 +      c14*tmp(:)**4.d0 +      c15*tmp(:)**5.d0 )/KPP0            
                dKPP(:,j,k) = c10*(       c11        + 2.d0*c12*tmp(:)       + 3.d0*c13*tmp(:)**2.d0 + 4.d0*c14*tmp(:)**3.d0 + 5.d0*c15*tmp(:)**4.d0 )/KPP0
                    
                invRhoCp(:,j,k) = Rho0*Cp0  &
                                 /(a10 + a10*(a11*tmp(:) + a12*tmp(:)**2.d0 + a13*tmp(:)**3.d0 + a14*tmp(:)**4.d0 + a15*tmp(:)**5.d0))   &
                                 /(b10 + b10*(b11*tmp(:) + b12*tmp(:)**2.d0 + b13*tmp(:)**3.d0 + b14*tmp(:)**4.d0 + b15*tmp(:)**5.d0))
                
                dinvRhoCp(:,j,k) = -( b11*(1.d0 + 3.d0*a12*tmp(:)**2.d0 + 4.d0*a13*tmp(:)**3.d0)                    &
                                     +a11*(1.d0 + 2.d0*b11*tmp(:) + 3.d0*b12*tmp(:)**2.d0 + 4.d0*b13*tmp(:)**3.d0)  &
                                     +tmp(:)*( b12*(2.d0 + 5.d0*a13*tmp(:)**3.d0)                                   &
                                              +a12*(2.d0 + 4.d0*b12*tmp(:)**2.d0 + 5.d0*b13*tmp(:)**3.d0)           &
                                              +3.d0*tmp(:)*(a13 + b13 + 2.d0*a13*b13*tmp(:)**3.d0)          )       &
                                    )/(1.d0 + a11*tmp(:) + a12*tmp(:)**2.d0 + a13*tmp(:)**3.d0)**2.d0               &
                                     /(1.d0 + b11*tmp(:) + b12*tmp(:)**2.d0 + b13*tmp(:)**3.d0)**2.d0
                
                                    
            enddo
        enddo
        deallocate(tmp)
    end subroutine mpi_thermal_coeffi

    subroutine mpi_thermal_coeffi_clean()
        implicit none
        
        deallocate(KPP,dKPP,invRhoCp,dinvRhoCp)
        
    end subroutine mpi_thermal_coeffi_clean
    
    subroutine mpi_thermal_solver(U,V,W,dx1,dx2,dx3,dmx1,dmx2,dmx3, comm_1d_x1, comm_1d_x2, comm_1d_x3)
        use mpi_topology, only : cart_comm_1d
        use PaScaL_TDMA
        implicit none
        double precision, dimension(0:m1sb,0:m2sb,0:m3sb) ::  U,V,W
        double precision :: dx1(0:m1sb),dx2(0:m2sb),dx3(0:m3sb)
        double precision :: dmx1(0:m1sb),dmx2(0:m2sb),dmx3(0:m3sb)
        type(cart_comm_1d), intent(in)  :: comm_1d_x1, comm_1d_x2, comm_1d_x3
        type(ptdma_plan_many)     :: ptdma_plan
        
        integer :: i,j,k
        integer :: im,jm,km
        integer :: ip,jp,kp
        integer :: jep,jem

        integer :: Sc,Ec
        integer :: Sm,Em
        integer :: Sp,Ep

        double precision, allocatable, dimension(:,:,:) :: dTc,RHS
        double precision, allocatable, dimension(:,:,:) :: am,ac,ap

        double precision, allocatable, dimension(:) :: e1,e2,e3,e4,e5,e6
        double precision, allocatable, dimension(:) :: dedx1,dedx2,dedy3,dedy4,dedz5,dedz6
        double precision, allocatable, dimension(:) :: convect_e1,convect_e2,convect_e3
        double precision, allocatable, dimension(:) :: KPP1,KPP2,KPP3,KPP4,KPP5,KPP6
        double precision, allocatable, dimension(:) :: dKPP1,dKPP2,dKPP3,dKPP4,dKPP5,dKPP6
        double precision, allocatable, dimension(:) :: viscous_e1,viscous_e2,viscous_e3
        double precision, allocatable, dimension(:) :: ebc_up,ebc_down
        double precision, allocatable, dimension(:) :: eACK,eAPK,eAMK,eACI,eAPI,eAMI,eACJ,eAPJ,eAMJ    ! 1D-32

        Sc=i_indexC  ;Ec=m1sb-1
        Sm=i_indexC-1;Em=m1sb-1-1
        Sp=i_indexC+1;Ep=m1sb-1+1
        
        allocate(e1(Sc:Ec),e2(Sc:Ec),e3(Sc:Ec))
        allocate(e4(Sc:Ec),e5(Sc:Ec),e6(Sc:Ec))
        allocate(dedx1(Sc:Ec),dedx2(Sc:Ec),dedy3(Sc:Ec))
        allocate(dedy4(Sc:Ec),dedz5(Sc:Ec),dedz6(Sc:Ec))
        allocate(convect_e1(Sc:Ec),convect_e2(Sc:Ec),convect_e3(Sc:Ec))
        allocate(KPP1(Sc:Ec),KPP2(Sc:Ec),KPP3(Sc:Ec))
        allocate(KPP4(Sc:Ec),KPP5(Sc:Ec),KPP6(Sc:Ec))
        allocate(dKPP1(Sc:Ec),dKPP2(Sc:Ec),dKPP3(Sc:Ec))
        allocate(dKPP4(Sc:Ec),dKPP5(Sc:Ec),dKPP6(Sc:Ec))
        allocate(viscous_e1(Sc:Ec),viscous_e2(Sc:Ec),viscous_e3(Sc:Ec))
        allocate(ebc_up(Sc:Ec),ebc_down(Sc:Ec))
        allocate(eACK(Sc:Ec),eAPK(Sc:Ec),eAMK(Sc:Ec))
        allocate(eACI(Sc:Ec),eAPI(Sc:Ec),eAMI(Sc:Ec))
        allocate(eACJ(Sc:Ec),eAPJ(Sc:Ec),eAMJ(Sc:Ec))

        allocate(RHS(1:m1sb-1,1:m2sb-1,1:m3sb-1))        

        ! i-j-k:RHS-dTc
        ! i-j-k
        allocate(am(1:m1sb-1,1:m2sb-1,1:m3sb-1))
        allocate(ac(1:m1sb-1,1:m2sb-1,1:m3sb-1))
        allocate(ap(1:m1sb-1,1:m2sb-1,1:m3sb-1))
        do k = k_indexC, m3sb-1
            km=k-1
            kp=k+1
            do j = j_indexC, m2sb-1
                jm=j-1
                jp=j+1
                jem = jC_BC(jm)
                jep = jC_BC(jp)
                    
                ! CONVECTION TERM                
                e1(Sc:Ec) = (dx1(Sc:Ec)*T(Sm:Em,j,k) + dx1(Sm:Em)*T(Sc:Ec,j,k))*(0.5d0/dmx1(Sc:Ec))
                e2(Sc:Ec) = (dx1(Sp:Ep)*T(Sc:Ec,j,k) + dx1(Sc:Ec)*T(Sp:Ep,j,k))*(0.5d0/dmx1(Sp:Ep))
                e3(Sc:Ec) = (dx2(j )*T(Sc:Ec,jm,k) + dx2(jm)*T(Sc:Ec,j,k ))*(0.5d0/dmx2(j ))
                e4(Sc:Ec) = (dx2(jp)*T(Sc:Ec,j,k ) + dx2(j )*T(Sc:Ec,jp,k))*(0.5d0/dmx2(jp))  
                e5(Sc:Ec) = (dx3(k)*T(Sc:Ec,j,km) + dx3(km)*T(Sc:Ec,j,k))*(0.5d0/dmx3(k ))
                e6(Sc:Ec) = (dx3(k)*T(Sc:Ec,j,kp) + dx3(kp)*T(Sc:Ec,j,k))*(0.5d0/dmx3(kp))
                dedx1(Sc:Ec) = (T(Sc:Ec,j,k) - T(Sm:Em,j,k))/dmx1(Sc:Ec)
                dedx2(Sc:Ec) = (T(Sp:Ep,j,k) - T(Sc:Ec,j,k))/dmx1(Sp:Ep)                      
                dedy3(Sc:Ec) = (T(Sc:Ec,j ,k) - T(Sc:Ec,jm,k))/dmx2(j )
                dedy4(Sc:Ec) = (T(Sc:Ec,jp,k) - T(Sc:Ec,j ,k))/dmx2(jp)
                dedz5(Sc:Ec) = (T(Sc:Ec,j,k ) - T(Sc:Ec,j,km))/dmx3(k )
                dedz6(Sc:Ec) = (T(Sc:Ec,j,kp) - T(Sc:Ec,j,k ))/dmx3(kp)

                ! 1/2 ( U(Tp)x + V(Tp)y + W(Tp)z )
                convect_e1(Sc:Ec) = (U(Sc:Ec,j,k)*dedx1(Sc:Ec) + U(Sp:Ep,j,k)*dedx2(Sc:Ec))*0.5d0
                convect_e2(Sc:Ec) = (V(Sc:Ec,j,k)*dedy3(Sc:Ec) + V(Sc:Ec,jp,k)*dedy4(Sc:Ec))*0.5d0 
                convect_e3(Sc:Ec) = (W(Sc:Ec,j,k)*dedz5(Sc:Ec) + W(Sc:Ec,j,kp)*dedz6(Sc:Ec))*0.5d0
                ! RHS(Sc:Ec,j,k) = -convect + viscous + ebc -RHS_e
                RHS(Sc:Ec,j,k) = -0.5d0*(convect_e1(Sc:Ec) + convect_e2(Sc:Ec) + convect_e3(Sc:Ec))
                
                ! DIFFUSION TERM
                KPP1(Sc:Ec) = 0.5d0/dmx1(Sc:Ec)*(dx1(Sc:Ec)*KPP(Sm:Em,j ,k ) + dx1(Sm:Em)*KPP(Sc:Ec,j ,k ))
                KPP2(Sc:Ec) = 0.5d0/dmx1(Sp:Ep)*(dx1(Sp:Ep)*KPP(Sc:Ec,j ,k ) + dx1(Sc:Ec)*KPP(Sp:Ep,j ,k )) 
                KPP3(Sc:Ec) = 0.5d0/dmx2(j )*(dx2(j )*KPP(Sc:Ec,jm,k ) + dx2(jm)*KPP(Sc:Ec,j ,k ))
                KPP4(Sc:Ec) = 0.5d0/dmx2(jp)*(dx2(jp)*KPP(Sc:Ec,j ,k ) + dx2(j )*KPP(Sc:Ec,jp,k )) 
                KPP5(Sc:Ec) = 0.5d0/dmx3(k )*(dx3(k )*KPP(Sc:Ec,j ,km) + dx3(km)*KPP(Sc:Ec,j ,k ))
                KPP6(Sc:Ec) = 0.5d0/dmx3(kp)*(dx3(kp)*KPP(Sc:Ec,j ,k ) + dx3(k )*KPP(Sc:Ec,j ,kp))

                dKPP1(Sc:Ec) = 0.5d0/dmx1(Sc:Ec)*(dx1(Sc:Ec)*dKPP(Sm:Em,j ,k ) + dx1(Sm:Em)*dKPP(Sc:Ec,j ,k ))
                dKPP2(Sc:Ec) = 0.5d0/dmx1(Sp:Ep)*(dx1(Sp:Ep)*dKPP(Sc:Ec,j ,k ) + dx1(Sc:Ec)*dKPP(Sp:Ep,j ,k )) 
                dKPP3(Sc:Ec) = 0.5d0/dmx2(j )*(dx2(j )*dKPP(Sc:Ec,jm,k ) + dx2(jm)*dKPP(Sc:Ec,j ,k ))
                dKPP4(Sc:Ec) = 0.5d0/dmx2(jp)*(dx2(jp)*dKPP(Sc:Ec,j ,k ) + dx2(j )*dKPP(Sc:Ec,jp,k )) 
                dKPP5(Sc:Ec) = 0.5d0/dmx3(k )*(dx3(k )*dKPP(Sc:Ec,j ,km) + dx3(km)*dKPP(Sc:Ec,j ,k ))
                dKPP6(Sc:Ec) = 0.5d0/dmx3(kp)*(dx3(kp)*dKPP(Sc:Ec,j ,k ) + dx3(k )*dKPP(Sc:Ec,j ,kp))

                viscous_e1(Sc:Ec) = (KPP2(Sc:Ec)*dedx2(Sc:Ec) - KPP1(Sc:Ec)*dedx1(Sc:Ec) - e2(Sc:Ec)*dKPP2(Sc:Ec)*dedx2(Sc:Ec) + e1(Sc:Ec)*dKPP1(Sc:Ec)*dedx1(Sc:Ec))/dx1(Sc:Ec)
                viscous_e2(Sc:Ec) = (KPP4(Sc:Ec)*dedy4(Sc:Ec) - KPP3(Sc:Ec)*dedy3(Sc:Ec) - e4(Sc:Ec)*dKPP4(Sc:Ec)*dedy4(Sc:Ec) + e3(Sc:Ec)*dKPP3(Sc:Ec)*dedy3(Sc:Ec))/dx2(j)
                viscous_e3(Sc:Ec) = (KPP6(Sc:Ec)*dedz6(Sc:Ec) - KPP5(Sc:Ec)*dedz5(Sc:Ec) - e6(Sc:Ec)*dKPP6(Sc:Ec)*dedz6(Sc:Ec) + e5(Sc:Ec)*dKPP5(Sc:Ec)*dedz5(Sc:Ec))/dx3(k)
                
                ! RHS(Sc:Ec,j,k) = -convect + viscous + ebc -RHS_e
                RHS(Sc:Ec,j,k) = RHS(Sc:Ec,j,k)     &
                                + 0.5d0*Ct* invRhoCp(Sc:Ec,j,k)*(viscous_e1(Sc:Ec) + viscous_e2(Sc:Ec) + viscous_e3(Sc:Ec))                      &
                                - 0.5d0*Ct*dinvRhoCp(Sc:Ec,j,k)*T(Sc:Ec,j,k)*( (KPP2(Sc:Ec)*dedx2(Sc:Ec) - KPP1(Sc:Ec)*dedx1(Sc:Ec))/dx1(Sc:Ec)  &
                                                                              +(KPP4(Sc:Ec)*dedy4(Sc:Ec) - KPP3(Sc:Ec)*dedy3(Sc:Ec))/dx2(j)      &
                                                                              +(KPP6(Sc:Ec)*dedz6(Sc:Ec) - KPP5(Sc:Ec)*dedz5(Sc:Ec))/dx3(k)      )  

                ! ebc For Y-direction- From Convection Terms
                ebc_up(Sc:Ec)   =-0.25d0*V(Sc:Ec,jp,k)/dmx2(jp)*TBCup_sub(Sc:Ec,k)
                ebc_down(Sc:Ec) = 0.25d0*V(Sc:Ec,j ,k)/dmx2(j )*TBCbt_sub(Sc:Ec,k)
                ! ebc For Y-direction- From Diffusion Terms
                ebc_up(Sc:Ec) = ebc_up (Sc:Ec)                                                                  &
                               + 0.5d0*Ct*invRhoCp(Sc:Ec,j,k)/dx2(j)* KPP4(Sc:Ec)/dmx2(jp)*TBCup_sub(Sc:Ec,k)    & 
                               + 0.5d0*Ct*invRhoCp(Sc:Ec,j,k)/dx2(j)*dKPP4(Sc:Ec)*dedy4(Sc:Ec)*TBCup_sub(Sc:Ec,k)   
                ebc_down(Sc:Ec) = ebc_down(Sc:Ec)                                                                   &
                                + 0.5d0*Ct*invRhoCp(Sc:Ec,j,k)/dx2(j)*KPP3(Sc:Ec)/dmx2(j)*TBCbt_sub(Sc:Ec,k)        &
                                - 0.5d0*Ct*invRhoCp(Sc:Ec,j,k)/dx2(j)*dKPP3(Sc:Ec)*dedy3(Sc:Ec)*TBCbt_sub(Sc:Ec,k)

                ! ebc For Y-direction
                ! RHS(Sc:Ec,j,k) = -convect + viscous + ebc -RHS_e
                RHS(Sc:Ec,j,k) = RHS(Sc:Ec,j,k)     &
                                + dble(1-jem)*ebc_down(Sc:Ec)   &
                                + dble(1-jep)*ebc_up(Sc:Ec)
                
                ! Z-DIRECTION
                eACK(Sc:Ec) = 0.5d0*Ct*invRhoCp(Sc:Ec,j,k)/dx3(k)*( KPP6(Sc:Ec)/dmx3(kp) + KPP5(Sc:Ec)/dmx3(k)             &
                                                                   -dx3(kp)*(0.5d0/dmx3(kp))*dKPP6(Sc:Ec)*dedz6(Sc:Ec)     &
                                                                   +dx3(km)*(0.5d0/dmx3(k ))*dKPP5(Sc:Ec)*dedz5(Sc:Ec) )   &
                            + (-0.25d0*W(Sc:Ec,j,kp)/dmx3(kp) + 0.25d0*W(Sc:Ec,j,k)/dmx3(k ))
                eAPK(Sc:Ec) =-0.5d0*Ct*invRhoCp(Sc:Ec,j,k)/dx3(k)*(KPP6(Sc:Ec)/dmx3(kp) + dx3(k)*(0.5d0/dmx3(kp))*dKPP6(Sc:Ec)*dedz6(Sc:Ec))    &
                            + 0.25d0*W(Sc:Ec,j,kp)/dmx3(kp)
                eAMK(Sc:Ec) =-0.5d0*Ct*invRhoCp(Sc:Ec,j,k)/dx3(k)*(KPP5(Sc:Ec)/dmx3(k ) - dx3(k)*(0.5d0/dmx3(k ))*dKPP5(Sc:Ec)*dedz5(Sc:Ec))    &
                            - 0.25d0*W(Sc:Ec,j,k )/dmx3(k )

                ! X-DIRECTION 
                eACI(Sc:Ec) = 0.5d0*Ct*invRhoCp(Sc:Ec,j,k)/dx1(Sc:Ec)*( KPP2(Sc:Ec)/dmx1(Sp:Ep) + KPP1(Sc:Ec)/dmx1(Sc:Ec)          &
                                                                       -dx1(Sp:Ep)*(0.5d0/dmx1(Sp:Ep))*dKPP2(Sc:Ec)*dedx2(Sc:Ec)   &
                                                                       +dx1(Sm:Em)*(0.5d0/dmx1(Sc:Ec))*dKPP1(Sc:Ec)*dedx1(Sc:Ec) ) &
                            + (0.25d0*U(Sc:Ec,j,k)/dmx1(Sc:Ec) - 0.25d0*U(Sp:Ep,j,k)/dmx1(Sp:Ep))   
                eAPI(Sc:Ec) =-0.5d0*Ct*invRhoCp(Sc:Ec,j,k)/dx1(Sc:Ec)*(KPP2(Sc:Ec)/dmx1(Sp:Ep) + dx1(Sc:Ec)*(0.5d0/dmx1(Sp:Ep))*dKPP2(Sc:Ec)*dedx2(Sc:Ec))    &
                            + 0.25d0*U(Sp:Ep,j,k)/dmx1(Sp:Ep)  
                eAMI(Sc:Ec) =-0.5d0*Ct*invRhoCp(Sc:Ec,j,k)/dx1(Sc:Ec)*(KPP1(Sc:Ec)/dmx1(Sc:Ec) - dx1(Sc:Ec)*(0.5d0/dmx1(Sc:Ec))*dKPP1(Sc:Ec)*dedx1(Sc:Ec))    &
                            - 0.25d0*U(Sc:Ec,j,k)/dmx1(Sc:Ec)
                        
                ! Y-DIRECTION   
                eACJ(Sc:Ec) = 0.5d0*Ct*invRhoCp(Sc:Ec,j,k)/dx2(j)*( KPP4(Sc:Ec)/dmx2(jp) + KPP3(Sc:Ec)/dmx2(j)                      &
                                                                   -dx2(jp)*(0.5d0/dmx2(jp))*dKPP4(Sc:Ec)*dedy4(Sc:Ec)*dble(jep)    &
                                                                   +dx2(jm)*(0.5d0/dmx2(j ))*dKPP3(Sc:Ec)*dedy3(Sc:Ec)*dble(jem) )  &
                            + (0.25d0*V(Sc:Ec,j,k)/dmx2(j ) - 0.25d0*V(Sc:Ec,jp,k)/dmx2(jp))   
                eAPJ(Sc:Ec) =-0.5d0*Ct*invRhoCp(Sc:Ec,j,k)/dx2(j)*(KPP4(Sc:Ec)/dmx2(jp) + dx2(j)*(0.5d0/dmx2(jp))*dKPP4(Sc:Ec)*dedy4(Sc:Ec))    &
                            + 0.25d0*V(Sc:Ec,jp,k)/dmx2(jp)
                eAMJ(Sc:Ec) =-0.5d0*Ct*invRhoCp(Sc:Ec,j,k)/dx2(j)*(KPP3(Sc:Ec)/dmx2(j ) - dx2(j)*(0.5d0/dmx2(j ))*dKPP3(Sc:Ec)*dedy3(Sc:Ec))    &
                            - 0.25d0*V(Sc:Ec,j ,k)/dmx2(j )                           
                eAPJ(Sc:Ec) = eAPJ(Sc:Ec)*dble(jep)    
                eAMJ(Sc:Ec) = eAMJ(Sc:Ec)*dble(jem)

                ! RHS(Sc:Ec,j,k) = -convect + viscous + ebc -RHS_e
                RHS(Sc:Ec,j,k) = RHS(Sc:Ec,j,k)     &
                                -( eAPK(Sc:Ec)*T(Sc:Ec,j ,kp) + eACK(Sc:Ec)*T(Sc:Ec,j ,k ) + eAMK(Sc:Ec)*T(Sc:Ec,j ,km)                                  &
                                  +eAPJ(Sc:Ec)*T(Sc:Ec,jp,k ) + eACJ(Sc:Ec)*T(Sc:Ec,j ,k ) + eAMJ(Sc:Ec)*T(Sc:Ec,jm,k )                                  &
                                  +eAPI(Sc:Ec)*T(Sp:Ep,j ,k ) + eACI(Sc:Ec)*T(Sc:Ec,j ,k ) + eAMI(Sc:Ec)*T(Sm:Em,j ,k )                                  &
                                  -0.5d0*Ct*T(Sc:Ec,j ,k )*dinvRhoCp(Sc:Ec,j ,k )*( (KPP2(Sc:Ec)*dedx2(Sc:Ec) - KPP1(Sc:Ec)*dedx1(Sc:Ec))/dx1(Sc:Ec)     &
                                                                                   +(KPP4(Sc:Ec)*dedy4(Sc:Ec) - KPP3(Sc:Ec)*dedy3(Sc:Ec))/dx2(j)         &
                                                                                   +(KPP6(Sc:Ec)*dedz6(Sc:Ec) - KPP5(Sc:Ec)*dedz5(Sc:Ec))/dx3(k)    )    )

                ! TDMA A- K Direction
                ac(Sc:Ec,j,k)= eACK(Sc:Ec)*dt   &
                             - 0.25d0*Ct*dinvRhoCp(Sc:Ec,j ,k )*( (KPP2(Sc:Ec)*dedx2(Sc:Ec) - KPP1(Sc:Ec)*dedx1(Sc:Ec))/dx1(Sc:Ec)      & 
                                                                 +(KPP4(Sc:Ec)*dedy4(Sc:Ec) - KPP3(Sc:Ec)*dedy3(Sc:Ec))/dx2(j)          &
                                                                 +(KPP6(Sc:Ec)*dedz6(Sc:Ec) - KPP5(Sc:Ec)*dedz5(Sc:Ec))/dx3(k)   )*dt   &
                             + 1.d0
                ap(Sc:Ec,j,k)=eAPK(Sc:Ec)*dt
                am(Sc:Ec,j,k)=eAMK(Sc:Ec)*dt
                RHS(Sc:Ec,j,k)=RHS(Sc:Ec,j,k)*dt

            enddo
        enddo
        ! TDMA-i*j-K
        call PaScaL_TDMA_plan_many_create(ptdma_plan, (m1sb-1)*(m2sb-1), comm_1d_x3%myrank, comm_1d_x3%nprocs, comm_1d_x3%mpi_comm)
        call PaScaL_TDMA_many_solve_cycle(ptdma_plan, am, ac, ap, RHS,(m1sb-1)*(m2sb-1),(m3sb-1))
        call PaScaL_TDMA_plan_many_destroy(ptdma_plan,comm_1d_x3%nprocs)
        deallocate(am,ac,ap)
        
        ! Trans i-j-K->j-k-i
        allocate(dTc(1:m2sb-1,1:m1sb-1,1:m3sb-1))
        do k = k_indexC, m3sb-1
            dTc(:,:,k)=TRANSPOSE(RHS(:,:,k))
        enddo
        deallocate(RHS) 
        allocate(RHS(1:m2sb-1,1:m3sb-1,1:m1sb-1))
        do k = k_indexC, m3sb-1
            do i = i_indexC, m1sb-1
                RHS(1:m2sb-1,k,i)=dTc(1:m2sb-1,i,k)
            enddo
        enddo
        deallocate(dTc)

        ! j-k-i
        allocate(am(1:m2sb-1,1:m3sb-1,1:m1sb-1),ac(1:m2sb-1,1:m3sb-1,1:m1sb-1),ap(1:m2sb-1,1:m3sb-1,1:m1sb-1))
        do k = k_indexC, m3sb-1
            do j = j_indexC, m2sb-1
                KPP1(Sc:Ec) = 0.5d0/dmx1(Sc:Ec)*(dx1(Sc:Ec)*KPP(Sm:Em,j ,k ) + dx1(Sm:Em)*KPP(Sc:Ec,j ,k ))
                KPP2(Sc:Ec) = 0.5d0/dmx1(Sp:Ep)*(dx1(Sp:Ep)*KPP(Sc:Ec,j ,k ) + dx1(Sc:Ec)*KPP(Sp:Ep,j ,k )) 
                dKPP1(Sc:Ec) = 0.5d0/dmx1(Sc:Ec)*(dx1(Sc:Ec)*dKPP(Sm:Em,j ,k ) + dx1(Sm:Em)*dKPP(Sc:Ec,j ,k ))
                dKPP2(Sc:Ec) = 0.5d0/dmx1(Sp:Ep)*(dx1(Sp:Ep)*dKPP(Sc:Ec,j ,k ) + dx1(Sc:Ec)*dKPP(Sp:Ep,j ,k )) 

                dedx1(Sc:Ec) = (T(Sc:Ec,j,k) - T(Sm:Em,j,k))/dmx1(Sc:Ec)
                dedx2(Sc:Ec) = (T(Sp:Ep,j,k) - T(Sc:Ec,j,k))/dmx1(Sp:Ep)                      
                dedy3(Sc:Ec) = (T(Sc:Ec,j ,k) - T(Sc:Ec,jm,k))/dmx2(j )
                dedy4(Sc:Ec) = (T(Sc:Ec,jp,k) - T(Sc:Ec,j ,k))/dmx2(jp)
                dedz5(Sc:Ec) = (T(Sc:Ec,j,k ) - T(Sc:Ec,j,km))/dmx3(k )
                dedz6(Sc:Ec) = (T(Sc:Ec,j,kp) - T(Sc:Ec,j,k ))/dmx3(kp)
                ! X-DIRECTION 
                eACI(Sc:Ec) = 0.5d0*Ct*invRhoCp(Sc:Ec,j,k)/dx1(Sc:Ec)*( KPP2(Sc:Ec)/dmx1(Sp:Ep) + KPP1(Sc:Ec)/dmx1(Sc:Ec)           &
                                                                       -dx1(Sp:Ep)*(0.5d0/dmx1(Sp:Ep))*dKPP2(Sc:Ec)*dedx2(Sc:Ec)    &
                                                                       +dx1(Sm:Em)*(0.5d0/dmx1(Sc:Ec))*dKPP1(Sc:Ec)*dedx1(Sc:Ec) )  &
                            + (0.25d0*U(Sc:Ec,j,k)/dmx1(Sc:Ec) - 0.25d0*U(Sp:Ep,j,k)/dmx1(Sp:Ep))                                   &                            
                            - 0.25d0*Ct*dinvRhoCp(Sc:Ec,j ,k )*( (KPP2(Sc:Ec)*dedx2(Sc:Ec) - KPP1(Sc:Ec)*dedx1(Sc:Ec))/dx1(i)       & 
                                                                +(KPP4(Sc:Ec)*dedy4(Sc:Ec) - KPP3(Sc:Ec)*dedy3(Sc:Ec))/dx2(j)       &
                                                                +(KPP6(Sc:Ec)*dedz6(Sc:Ec) - KPP5(Sc:Ec)*dedz5(Sc:Ec))/dx3(k)   )
                eAPI(Sc:Ec) =-0.5d0*Ct*invRhoCp(Sc:Ec,j,k)/dx1(Sc:Ec)*(KPP2(Sc:Ec)/dmx1(Sp:Ep) + dx1(Sc:Ec)*(0.5d0/dmx1(Sp:Ep))*dKPP2(Sc:Ec)*dedx2(Sc:Ec))    &
                            + 0.25d0*U(Sp:Ep,j,k)/dmx1(Sp:Ep)
                eAMI(Sc:Ec) =-0.5d0*Ct*invRhoCp(Sc:Ec,j,k)/dx1(Sc:Ec)*(KPP1(Sc:Ec)/dmx1(Sc:Ec) - dx1(Sc:Ec)*(0.5d0/dmx1(Sc:Ec))*dKPP1(Sc:Ec)*dedx1(Sc:Ec))    &
                            - 0.25d0*U(Sc:Ec,j,k)/dmx1(Sc:Ec)
                            
                do i = i_indexC, m1sb-1
                    ac(j,k,i)=eACI(i)*dt + 1.d0
                    ap(j,k,i)=eAPI(i)*dt
                    am(j,k,i)=eAMI(i)*dt
                enddo
            enddo
        enddo
        ! TDMA-j*K-i
        call PaScaL_TDMA_plan_many_create(ptdma_plan, (m2sb-1)*(m3sb-1), comm_1d_x1%myrank, comm_1d_x1%nprocs, comm_1d_x1%mpi_comm)
        call PaScaL_TDMA_many_solve_cycle(ptdma_plan, am, ac, ap, RHS,(m2sb-1)*(m3sb-1),(m1sb-1))
        call PaScaL_TDMA_plan_many_destroy(ptdma_plan,comm_1d_x1%nprocs)
        deallocate(am,ac,ap)
        
        ! Trans j-k-i->k-i-j 
        allocate(dTc(1:m3sb-1,1:m2sb-1,1:m1sb-1))
        do i = i_indexC, m1sb-1
            dTc(:,:,i)=TRANSPOSE(RHS(:,:,i))
        enddo
        deallocate(RHS) 
        allocate(RHS(1:m3sb-1,1:m1sb-1,1:m2sb-1))
        do j = j_indexC, m2sb-1
            do i = i_indexC, m1sb-1
                RHS(1:m3sb-1,i,j)=dTc(1:m3sb-1,j,i)
            enddo
        enddo
        deallocate(dTc)

        ! k-i-j
        allocate(am(1:m3sb-1,1:m1sb-1,1:m2sb-1),ac(1:m3sb-1,1:m1sb-1,1:m2sb-1),ap(1:m3sb-1,1:m1sb-1,1:m2sb-1))
        do k = k_indexC, m3sb-1
            do j = j_indexC, m2sb-1
                jm=j-1
                jp=j+1
                jem = jC_BC(jm)
                jep = jC_BC(jp)
                
                KPP3(Sc:Ec) = 0.5d0/dmx2(j )*(dx2(j )*KPP(Sc:Ec,jm,k ) + dx2(jm)*KPP(Sc:Ec,j ,k ))
                KPP4(Sc:Ec) = 0.5d0/dmx2(jp)*(dx2(jp)*KPP(Sc:Ec,j ,k ) + dx2(j )*KPP(Sc:Ec,jp,k )) 
                dKPP3(Sc:Ec) = 0.5d0/dmx2(j )*(dx2(j )*dKPP(Sc:Ec,jm,k ) + dx2(jm)*dKPP(Sc:Ec,j ,k ))
                dKPP4(Sc:Ec) = 0.5d0/dmx2(jp)*(dx2(jp)*dKPP(Sc:Ec,j ,k ) + dx2(j )*dKPP(Sc:Ec,jp,k )) 

                dedx1(Sc:Ec) = (T(Sc:Ec,j,k) - T(Sm:Em,j,k))/dmx1(Sc:Ec)
                dedx2(Sc:Ec) = (T(Sp:Ep,j,k) - T(Sc:Ec,j,k))/dmx1(Sp:Ep)                      
                dedy3(Sc:Ec) = (T(Sc:Ec,j ,k) - T(Sc:Ec,jm,k))/dmx2(j )
                dedy4(Sc:Ec) = (T(Sc:Ec,jp,k) - T(Sc:Ec,j ,k))/dmx2(jp)
                dedz5(Sc:Ec) = (T(Sc:Ec,j,k ) - T(Sc:Ec,j,km))/dmx3(k )
                dedz6(Sc:Ec) = (T(Sc:Ec,j,kp) - T(Sc:Ec,j,k ))/dmx3(kp)
                
                ! Y-DIRECTION   
                eACJ(Sc:Ec) = 0.5d0*Ct*invRhoCp(Sc:Ec,j,k)/dx2(j)*( KPP4(Sc:Ec)/dmx2(jp) + KPP3(Sc:Ec)/dmx2(j)                      &
                                                                   -dx2(jp)*(0.5d0/dmx2(jp))*dKPP4(Sc:Ec)*dedy4(Sc:Ec)*dble(jep)    &
                                                                   +dx2(jm)*(0.5d0/dmx2(j ))*dKPP3(Sc:Ec)*dedy3(Sc:Ec)*dble(jem) )  &
                            + (0.25d0*V(Sc:Ec,j,k)/dmx2(j ) - 0.25d0*V(Sc:Ec,jp,k)/dmx2(jp))   
                eAPJ(Sc:Ec) =-0.5d0*Ct*invRhoCp(Sc:Ec,j,k)/dx2(j)*(KPP4(Sc:Ec)/dmx2(jp) + dx2(j)*(0.5d0/dmx2(jp))*dKPP4(Sc:Ec)*dedy4(Sc:Ec))    &
                            + 0.25d0*V(Sc:Ec,jp,k)/dmx2(jp)
                eAMJ(Sc:Ec) =-0.5d0*Ct*invRhoCp(Sc:Ec,j,k)/dx2(j)*(KPP3(Sc:Ec)/dmx2(j ) - dx2(j)*(0.5d0/dmx2(j ))*dKPP3(Sc:Ec)*dedy3(Sc:Ec))    &
                            - 0.25d0*V(Sc:Ec,j ,k)/dmx2(j )                           
                eAPJ(Sc:Ec) = eAPJ(Sc:Ec)*dble(jep)    
                eAMJ(Sc:Ec) = eAMJ(Sc:Ec)*dble(jem)

                do i = i_indexC, m1sb-1
                    ac(k,i,j)=eACJ(i)*dt + 1.d0
                    ap(k,i,j)=eAPJ(i)*dt
                    am(k,i,j)=eAMJ(i)*dt
                enddo
            enddo
        enddo
        ! TDMA-k*i-j
        call PaScaL_TDMA_plan_many_create(ptdma_plan, (m3sb-1)*(m1sb-1), comm_1d_x2%myrank, comm_1d_x2%nprocs, comm_1d_x2%mpi_comm)
        call PaScaL_TDMA_many_solve(ptdma_plan, am, ac, ap, RHS,(m3sb-1)*(m1sb-1),(m2sb-1))
        call PaScaL_TDMA_plan_many_destroy(ptdma_plan,comm_1d_x2%nprocs)
        deallocate(am,ac,ap)

        ! Trans k-i-j->i-j-k
        allocate(dTc(1:m1sb-1,1:m3sb-1,1:m2sb-1))
        do j = j_indexC, m2sb-1
            dTc(:,:,j)=TRANSPOSE(RHS(:,:,j))
        enddo
        deallocate(RHS)

        ! Update
        do k = k_indexC, m3sb-1
            do j = j_indexC, m2sb-1
                T(1:m1sb-1,j,k) = T(1:m1sb-1,j,k) + dTc(1:m1sb-1,k,j)
            enddo
        enddo
        deallocate(dTc)

        deallocate(e1,e2,e3,e4,e5,e6)
        deallocate(dedx1,dedx2,dedy3,dedy4,dedz5,dedz6)
        deallocate(convect_e1,convect_e2,convect_e3)
        deallocate(KPP1,KPP2,KPP3,KPP4,KPP5,KPP6)
        deallocate(dKPP1,dKPP2,dKPP3,dKPP4,dKPP5,dKPP6)
        deallocate(viscous_e1,viscous_e2,viscous_e3)
        deallocate(ebc_up,ebc_down)
        deallocate(eACK,eAPK,eAMK,eACI,eAPI,eAMI,eACJ,eAPJ,eAMJ)

    end subroutine mpi_thermal_solver

end module mpi_thermal

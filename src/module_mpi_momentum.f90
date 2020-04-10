module mpi_momentum
    ! use debug
    use global
    implicit none
    integer :: n1sb,n2sb,n3sb
    double precision, allocatable, dimension(:,:,:) :: U,V,W,P
    double precision, allocatable, dimension(:,:,:) :: dU,dV,dW
    double precision, allocatable, dimension(:,:) :: UBCup_sub,VBCup_sub,WBCup_sub
    double precision, allocatable, dimension(:,:) :: UBCbt_sub,VBCbt_sub,WBCbt_sub
    
    double precision, allocatable, dimension(:,:,:) :: invRho,Mu
contains    
    subroutine mpi_momentum_allocation(n1in,n2in,n3in)
        implicit none
        integer :: n1in,n2in,n3in

        n1sb=n1in
        n2sb=n2in
        n3sb=n3in

        allocate(U(0:n1sb,0:n2sb,0:n3sb),V(0:n1sb,0:n2sb,0:n3sb),W(0:n1sb,0:n2sb,0:n3sb))
        allocate(P(0:n1sb,0:n2sb,0:n3sb))

        allocate( UBCup_sub(0:n1sb,0:n3sb),VBCup_sub(0:n1sb,0:n3sb),WBCup_sub(0:n1sb,0:n3sb))
        allocate( UBCbt_sub(0:n1sb,0:n3sb),VBCbt_sub(0:n1sb,0:n3sb),WBCbt_sub(0:n1sb,0:n3sb))

    end subroutine mpi_momentum_allocation

    subroutine mpi_momentum_clean()
        implicit none
        
        deallocate(U,V,W)
        deallocate(P)

        deallocate( UBCup_sub,VBCup_sub,WBCup_sub)
        deallocate( UBCbt_sub,VBCbt_sub,WBCbt_sub)
        
    end subroutine mpi_momentum_clean

    subroutine mpi_momentum_initial(x1,x2,x3)
        implicit none
        double precision:: x1(0:n1sb),x2(0:n2sb),x3(0:n3sb)
        integer :: i,j,k

            U(0:n1sb,0:n2sb,0:n3sb)=0.d0
            V(0:n1sb,0:n2sb,0:n3sb)=0.d0
            W(0:n1sb,0:n2sb,0:n3sb)=0.d0
            P(0:n1sb,0:n2sb,0:n3sb)=0.d0

            ! call random_seed() 
            ! do k=1,n3sb-1
            !     do j=1,n2sb-1
            !         do i=1,n1sb-1                    
            !             call random_number(U(i,j,k))
            !             call random_number(V(i,j,k))
            !             call random_number(W(i,j,k))
            !         enddo
            !     enddo
            ! enddo

    end subroutine mpi_momentum_initial

    subroutine mpi_momentum_boundary(comm_1d_x1, comm_1d_x2, comm_1d_x3)
        use mpi_topology, only : cart_comm_1d
        implicit none
        type(cart_comm_1d), intent(in)  :: comm_1d_x1, comm_1d_x2, comm_1d_x3
        integer :: k
        
        UBCup_sub(:,:)=Uup
        VBCup_sub(:,:)=Vup
        WBCup_sub(:,:)=Wup
        
        UBCbt_sub(:,:)=Ubt
        VBCbt_sub(:,:)=Vbt
        WBCbt_sub(:,:)=Wbt

        if(comm_1d_x2%myrank==0) then   !No i,j Dirichlet 
            do k=0, n3sb
                U(0:n1sb,0,k)=UBCbt_sub(:,k)
                W(0:n1sb,0,k)=WBCbt_sub(:,k)
                V(0:n1sb,0,k)=VBCbt_sub(:,k)
                V(0:n1sb,1,k)=VBCbt_sub(:,k)
            enddo
        else if(comm_1d_x2%myrank==comm_1d_x2%nprocs-1) then
            do k=0, n3sb
                U(0:n1sb,n2sb,k)=UBCup_sub(:,k)
                W(0:n1sb,n2sb,k)=WBCup_sub(:,k)
                V(0:n1sb,n2sb,k)=VBCup_sub(:,k)
            enddo
        endif    
    end subroutine mpi_momentum_boundary

    subroutine mpi_momentum_coeffi(T)
        implicit none
        double precision :: T(0:n1sb,0:n2sb,0:n3sb)
        double precision, allocatable, dimension(:) :: tmp
        integer :: i,j,k

        allocate(invRho(0:n1sb,0:n2sb,0:n3sb),Mu(0:n1sb,0:n2sb,0:n3sb))
        allocate(tmp(0:n1sb))
        do k=0,n3sb
            do j=0,n2sb
                tmp(:)=DeltaT*T(:,j,k)

                Mu(:,j,k) = (a10 + a10*(a11*tmp(:) + a12*tmp(:)**2. + a13*tmp(:)**3. + a14*tmp(:)**4. + a15*tmp(:)**5.))   &
                           *(d10 + d10*(d11*tmp(:) + d12*tmp(:)**2. + d13*tmp(:)**3. + d14*tmp(:)**4. + d15*tmp(:)**5.))/Mu0

                invRho(:,j,k) = Rho0/(a10 + a10*(a11*tmp(:) + a12*tmp(:)**2.d0 + a13*tmp(:)**3.d0 + a14*tmp(:)**4.d0 + a15*tmp(:)**5.d0))


            enddo
        enddo
        deallocate(tmp)

    end subroutine mpi_momentum_coeffi

    subroutine mpi_momentum_pseudoupdateUVW()
        implicit none
        integer :: i,j,k

        do k = 1, n3sb-1
            do j = 2, n2sb-1
                U(1:n1sb-1,j,k)=U(1:n1sb-1,j,k)+dU(1:n1sb-1,j,k)
                V(1:n1sb-1,j,k)=V(1:n1sb-1,j,k)+dV(1:n1sb-1,j,k)
                W(1:n1sb-1,j,k)=W(1:n1sb-1,j,k)+dW(1:n1sb-1,j,k)
            end do
            j=1
            U(1:n1sb-1,j,k)=U(1:n1sb-1,j,k)+dU(1:n1sb-1,j,k)
            V(1:n1sb-1,j,k)=V(1:n1sb-1,j,k)+dV(1:n1sb-1,j,k)*dble(2-j_indexS)
            W(1:n1sb-1,j,k)=W(1:n1sb-1,j,k)+dW(1:n1sb-1,j,k)
        end do
        deallocate(dU)
        deallocate(dV)
        deallocate(dW)
        
    end subroutine mpi_momentum_pseudoupdateUVW

    subroutine mpi_momentum_coeffi_clean()
        implicit none
        deallocate(invRho,Mu)
        
    end subroutine mpi_momentum_coeffi_clean

    subroutine mpi_momentum_solvedU(T,dx1,dx2,dx3,dmx1,dmx2,dmx3, comm_1d_x1, comm_1d_x2, comm_1d_x3)
        use mpi_topology, only : cart_comm_1d
        use PaScaL_TDMA
        implicit none
        double precision, dimension(0:n1sb,0:n2sb,0:n3sb) ::  T
        double precision :: dx1(0:n1sb),dx2(0:n2sb),dx3(0:n3sb)
        double precision :: dmx1(0:n1sb),dmx2(0:n2sb),dmx3(0:n3sb)
        type(cart_comm_1d), intent(in)  :: comm_1d_x1, comm_1d_x2, comm_1d_x3
        type(ptdma_plan_many)     :: ptdma_plan
            
        integer :: i,j,k
        integer :: im,jm,km
        integer :: ip,jp,kp
        integer :: jup,jum

        integer :: Sc,Ec
        integer :: Sm,Em
        integer :: Sp,Ep

        double precision, allocatable, dimension(:,:,:) :: ddU,RHS
        double precision, allocatable, dimension(:,:,:) :: am,ac,ap

        double precision, allocatable, dimension(:) :: u1,u2,u5,u6,v3,v4,w5,w6,Tc
        double precision, allocatable, dimension(:) :: dudx1,dudx2,dudy3,dudy4,dudz5,dudz6,dvdx3,dvdx4,dwdx5,dwdx6
        double precision, allocatable, dimension(:) :: mua,mub,muc,mu3,mu4,mu5,mu6
        double precision, allocatable, dimension(:) :: invRhoc,viscous_u1,viscous_u2,viscous_u3,viscous_u12,viscous_u13
        double precision, allocatable, dimension(:) :: ubc_up,ubc_down
        double precision, allocatable, dimension(:) :: mACI,mAPI,mAMI,mACJ,mAPJ,mAMJ,mACK,mAPK,mAMK

        Sc=i_indexS  ;Ec=n1sb-1
        Sm=i_indexS-1;Em=n1sb-1-1
        Sp=i_indexS+1;Ep=n1sb-1+1
        
        allocate(u1(Sc:Ec),u2(Sc:Ec),u5(Sc:Ec),u6(Sc:Ec),v3(Sc:Ec),v4(Sc:Ec),w5(Sc:Ec),w6(Sc:Ec),Tc(Sc:Ec))
        allocate(dudx1(Sc:Ec),dudx2(Sc:Ec),dudy3(Sc:Ec),dudy4(Sc:Ec),dudz5(Sc:Ec),dudz6(Sc:Ec),dvdx3(Sc:Ec),dvdx4(Sc:Ec),dwdx5(Sc:Ec),dwdx6(Sc:Ec))
        allocate(mua(Sc:Ec),mub(Sc:Ec),muc(Sc:Ec),mu3(Sc:Ec),mu4(Sc:Ec),mu5(Sc:Ec),mu6(Sc:Ec))
        allocate(invRhoc(Sc:Ec),viscous_u1(Sc:Ec),viscous_u2(Sc:Ec),viscous_u3(Sc:Ec),viscous_u12(Sc:Ec),viscous_u13(Sc:Ec))
        allocate(ubc_up(Sc:Ec),ubc_down(Sc:Ec))
        allocate(mACI(Sc:Ec),mAPI(Sc:Ec),mAMI(Sc:Ec),mACJ(Sc:Ec),mAPJ(Sc:Ec),mAMJ(Sc:Ec),mACK(Sc:Ec),mAPK(Sc:Ec),mAMK(Sc:Ec))

        allocate(RHS(1:n1sb-1,1:n2sb-1,1:n3sb-1))        

        ! i-j-k:RHS-ddU
        ! i-j-k
        allocate(am(1:n1sb-1,1:n2sb-1,1:n3sb-1))
        allocate(ac(1:n1sb-1,1:n2sb-1,1:n3sb-1))
        allocate(ap(1:n1sb-1,1:n2sb-1,1:n3sb-1))
        do k = k_indexC, n3sb-1
            km=k-1
            kp=k+1
            do j = j_indexC, n2sb-1
                jm = j-1
                jp = j+1
                jum = jC_BC(jm)
                jup = jC_BC(jp)
                            
                u1(Sc:Ec) = 0.5d0*(U(Sm:Em,j,k) + U(Sc:Ec,j,k))            
                u2(Sc:Ec) = 0.5d0*(U(Sp:Ep,j,k) + U(Sc:Ec,j,k))             

                v3(Sc:Ec) = 0.5d0*(dx1(Sm:Em)*V(Sc:Ec,j ,k) + dx1(Sc:Ec)*V(Sm:Em,j ,k))/dmx1(Sc:Ec)            
                v4(Sc:Ec) = 0.5d0*(dx1(Sm:Em)*V(Sc:Ec,jp,k) + dx1(Sc:Ec)*V(Sm:Em,jp,k))/dmx1(Sc:Ec)
                
                w5(Sc:Ec) = 0.5d0*(dx1(Sm:Em)*W(Sc:Ec,j,k ) + dx1(Sc:Ec)*W(Sm:Em,j,k ))/dmx1(Sc:Ec)
                w6(Sc:Ec) = 0.5d0*(dx1(Sm:Em)*W(Sc:Ec,j,kp) + dx1(Sc:Ec)*W(Sm:Em,j,kp))/dmx1(Sc:Ec)

                dudx1(Sc:Ec) = (U(Sc:Ec,j ,k ) - U(Sm:Em,j ,k ))/dx1(Sm:Em)
                dudx2(Sc:Ec) = (U(Sp:Ep,j ,k ) - U(Sc:Ec,j ,k ))/dx1(Sc:Ec)
                dudy3(Sc:Ec) = (U(Sc:Ec,j ,k ) - U(Sc:Ec,jm,k ))/dmx2(j )
                dudy4(Sc:Ec) = (U(Sc:Ec,jp,k ) - U(Sc:Ec,j ,k ))/dmx2(jp)
                dudz5(Sc:Ec) = (U(Sc:Ec,j ,k ) - U(Sc:Ec,j ,km))/dmx3(k )
                dudz6(Sc:Ec) = (U(Sc:Ec,j ,kp) - U(Sc:Ec,j ,k ))/dmx3(kp)

                dvdx3(Sc:Ec) = (V(Sc:Ec,j ,k) - V(Sm:Em,j ,k))/dmx1(Sc:Ec)
                dvdx4(Sc:Ec) = (V(Sc:Ec,jp,k) - V(Sm:Em,jp,k))/dmx1(Sc:Ec)

                dwdx5(Sc:Ec) = (W(Sc:Ec,j,k ) - W(Sm:Em,j,k ))/dmx1(Sc:Ec)
                dwdx6(Sc:Ec) = (W(Sc:Ec,j,kp) - W(Sm:Em,j,kp))/dmx1(Sc:Ec)
                                    
                mua(Sc:Ec) = 0.5d0*(dx1(Sm:Em)*Mu(Sc:Ec,jm,k ) + dx1(Sc:Ec)*Mu(Sm:Em,jm,k ))/dmx1(Sc:Ec)
                muc(Sc:Ec) = 0.5d0*(dx1(Sm:Em)*Mu(Sc:Ec,j ,k ) + dx1(Sc:Ec)*Mu(Sm:Em,j ,k ))/dmx1(Sc:Ec)
                mub(Sc:Ec) = 0.5d0*(dx1(Sm:Em)*Mu(Sc:Ec,jp,k ) + dx1(Sc:Ec)*Mu(Sm:Em,jp,k ))/dmx1(Sc:Ec)
                
                mu3(Sc:Ec) = 0.5d0*(dx2(jm)*muc + dx2(j)*mua)/dmx2(j )
                mu4(Sc:Ec) = 0.5d0*(dx2(jp)*muc + dx2(j)*mub)/dmx2(jp)
                
                mua(Sc:Ec) = 0.5d0*(dx1(Sm:Em)*Mu(Sc:Ec,j ,km) + dx1(Sc:Ec)*Mu(Sm:Em,j ,km))/dmx1(Sc:Ec)
                muc(Sc:Ec) = 0.5d0*(dx1(Sm:Em)*Mu(Sc:Ec,j ,k ) + dx1(Sc:Ec)*Mu(Sm:Em,j ,k ))/dmx1(Sc:Ec)
                mub(Sc:Ec) = 0.5d0*(dx1(Sm:Em)*Mu(Sc:Ec,j ,kp) + dx1(Sc:Ec)*Mu(Sm:Em,j ,kp))/dmx1(Sc:Ec)
                
                mu5(Sc:Ec) = 0.5d0*(dx3(km)*muc + dx3(k)*mua)/dmx3(k )
                mu6(Sc:Ec) = 0.5d0*(dx3(kp)*muc + dx3(k)*mub)/dmx3(kp)
                
                invRhoc(Sc:Ec) = 0.5d0*(dx1(Sm:Em)*invRho(Sc:Ec,j ,k )+dx1(Sc:Ec)*invRho(Sm:Em,j ,k ) )/dmx1(Sc:Ec)

                !------viscous
                viscous_u1(Sc:Ec) = 1.d0*(Mu(Sc:Ec,j ,k )*dudx2(Sc:Ec) - Mu(Sm:Em,j ,k )*dudx1(Sc:Ec))/dmx1(Sc:Ec)
                viscous_u2(Sc:Ec) = 1.d0*(mu4(Sc:Ec)*dudy4(Sc:Ec) - mu3(Sc:Ec)*dudy3(Sc:Ec))/dx2(j )
                viscous_u3(Sc:Ec) = 1.d0*(mu6(Sc:Ec)*dudz6(Sc:Ec) - mu5(Sc:Ec)*dudz5(Sc:Ec))/dx3(k )            
                viscous_u12(Sc:Ec) = 1.d0*(mu4(Sc:Ec)*dvdx4(Sc:Ec) - mu3(Sc:Ec)*dvdx3(Sc:Ec))/dx2(j)
                viscous_u13(Sc:Ec) = 1.d0*(mu6(Sc:Ec)*dwdx6(Sc:Ec) - mu5(Sc:Ec)*dwdx5(Sc:Ec))/dx3(k)

                RHS(Sc:Ec,j,k) = 0.5d0*Cmu*invRhoc(Sc:Ec)*(2.*viscous_u1(Sc:Ec) + viscous_u2(Sc:Ec) + viscous_u3(Sc:Ec) + viscous_u12(Sc:Ec) + viscous_u13(Sc:Ec))

                !------PRESSURE
                RHS(Sc:Ec,j,k) = RHS(Sc:Ec,j,k) &
                            - Cmp*invRhoc(Sc:Ec)*(P(Sc:Ec,j,k) - P(Sm:Em,j,k))/dmx1(Sc:Ec)

                !---!BUOYANCY TERM (we use Cmt*theta for the buoyancy term)
                Tc(Sc:Ec) = 0.5d0*(dx1(Sm:Em)*T(Sc:Ec,j,k) + dx1(Sc:Ec)*T(Sm:Em,j,k))/dmx1(Sc:Ec)
                !THETAy = Cmt*((Tc(Sc:Ec) - Thetam))**b
                !THETAy = Cmt*(Tc(Sc:Ec))
                RHS(Sc:Ec,j,k) = RHS(Sc:Ec,j,k)     &
                            + Cmt*(Tc(Sc:Ec) + a12pera11*Tc(Sc:Ec)**2.*DeltaT)*invRhoc(Sc:Ec)

                !------UBC
                !FROM CONVECTION TERM
                ubc_down(Sc:Ec) = (0.25d0*v3(Sc:Ec)/dmx2(j)*UBCbt_sub(Sc:Ec,k)                                             &
                    &     - 0.25d0*dudy3(Sc:Ec)*(dx1(Sc:Ec)*VBCbt_sub(Sm:Em,k) + dx1(Sm:Em)*VBCbt_sub(Sc:Ec,k))/dmx1(Sc:Ec)*0.5d0)
                ubc_up(Sc:Ec) = (-0.25d0*v4(Sc:Ec)/dmx2(jp)*UBCup_sub(Sc:Ec,k)                                             &
                    &     -0.25d0*dudy4(Sc:Ec)*(dx1(Sc:Ec)*VBCup_sub(Sm:Em,k) + dx1(Sm:Em)*VBCup_sub(Sc:Ec,k))/dmx1(Sc:Ec)*0.5d0)

                !FROM DIFFUSION TERM              
                ubc_down(Sc:Ec) = ubc_down(Sc:Ec) + 0.5d0*Cmu*invRhoc(Sc:Ec)/dx2(j)*mu3(Sc:Ec)/dmx2(j )*UBCbt_sub(Sc:Ec,k)
                ubc_down(Sc:Ec) = ubc_down(Sc:Ec) - 0.5d0*Cmu*invRhoc(Sc:Ec)/dx2(j)*mu3(Sc:Ec)/dmx1(Sc:Ec)*(VBCbt_sub(Sc:Ec,k) - VBCbt_sub(Sm:Em,k))

                ubc_up(Sc:Ec)   = ubc_up(Sc:Ec)   + 0.5d0*Cmu*invRhoc(Sc:Ec)/dx2(j)*mu4(Sc:Ec)/dmx2(jp)*UBCup_sub(Sc:Ec,k)
                ubc_up(Sc:Ec)   = ubc_up(Sc:Ec)   + 0.5d0*Cmu*invRhoc(Sc:Ec)/dx2(j)*mu4(Sc:Ec)/dmx1(Sc:Ec)*(VBCup_sub(Sc:Ec,k) - VBCup_sub(Sm:Em,k))

                RHS(Sc:Ec,j,k) = RHS(Sc:Ec,j,k) &
                            + dble(1.d0 - jum)*ubc_down(Sc:Ec)    &
                            + dble(1.d0 - jup)*ubc_up(Sc:Ec)

                !M11Un
                !X-direction
                mACI(Sc:Ec) = 0.5d0*Cmu*invRhoc(Sc:Ec)/dmx1(Sc:Ec)*(Mu(Sc:Ec,j ,k )/dx1(Sc:Ec) + Mu(Sc:Ec,j ,k )/dx1(Sm:Em))* 2.   &
                            + 0.25d0*(dudx2(Sc:Ec)*0.5d0 + dudx1(Sc:Ec)*0.5d0 - u2(Sc:Ec)/dx1(Sc:Ec) + u1(Sc:Ec)/dx1(Sm:Em) )
                mAPI(Sc:Ec) = -0.5d0*Cmu*invRhoc(Sc:Ec)/dmx1(Sc:Ec)*Mu(Sc:Ec,j ,k )/dx1(Sc:Ec)* 2.  &
                            + 0.25d0*( u2(Sc:Ec)/dx1(Sc:Ec) + dudx2(Sc:Ec)*0.5d0) 
                mAMI(Sc:Ec) = -0.5d0*Cmu*invRhoc(Sc:Ec)/dmx1(Sc:Ec)*Mu(Sc:Ec,j ,k )/dx1(Sm:Em)* 2.  &
                            + 0.25d0*(-u1(Sc:Ec)/dx1(Sm:Em) + dudx1(Sc:Ec)*0.5d0)

                !Y-direction
                mACJ(Sc:Ec) = 0.5d0*Cmu*invRhoc(Sc:Ec)/dx2(j)*(mu4(Sc:Ec)/dmx2(jp) + mu3(Sc:Ec)/dmx2(j))  &
                            + 0.25d0*(-v4(Sc:Ec)/dmx2(jp) + v3(Sc:Ec)/dmx2(j ))
                mAPJ(Sc:Ec) = -0.5d0*Cmu*invRhoc(Sc:Ec)/dx2(j)*mu4(Sc:Ec)/dmx2(jp)                   &
                            + 0.25d0*( v4(Sc:Ec)/dmx2(jp))
                mAMJ(Sc:Ec) = -0.5d0*Cmu*invRhoc(Sc:Ec)/dx2(j)*mu3(Sc:Ec)/dmx2(j)                        &
                            + 0.25d0*(-v3(Sc:Ec)/dmx2(j))
                mAPJ(Sc:Ec) = mAPJ(Sc:Ec)*dble(jup)
                mAMJ(Sc:Ec) = mAMJ(Sc:Ec)*dble(jum)
                
                !Z-direction
                mACK(Sc:Ec) = 0.5d0*Cmu*invRhoc(Sc:Ec)/dx3(k)*(mu6(Sc:Ec)/dmx3(kp) + mu5(Sc:Ec)/dmx3(k))    &
                            + 0.25d0*(-w6(Sc:Ec)/dmx3(kp) + w5(Sc:Ec)/dmx3(k))
                mAPK(Sc:Ec) = -0.5d0*Cmu*invRhoc(Sc:Ec)/dx3(k)*mu6(Sc:Ec)/dmx3(kp)  &
                            + 0.25d0*( w6(Sc:Ec)/dmx3(kp))
                mAMK(Sc:Ec) = -0.5d0*Cmu*invRhoc(Sc:Ec)/dx3(k)*mu5(Sc:Ec)/dmx3(k )  &
                            + 0.25d0*(-w5(Sc:Ec)/dmx3(k))

                RHS(Sc:Ec,j,k) = RHS(Sc:Ec,j,k)                     &
                            - ( mAPI(Sc:Ec)*U(Sp:Ep,j ,k ) + mACI(Sc:Ec)*U(Sc:Ec,j ,k) + mAMI(Sc:Ec)*U(Sm:Em,j ,k )      &
                                +mAPJ(Sc:Ec)*U(Sc:Ec,jp,k ) + mACJ(Sc:Ec)*U(Sc:Ec,j ,k) + mAMJ(Sc:Ec)*U(Sc:Ec,jm,k )      &
                                +mAPK(Sc:Ec)*U(Sc:Ec,j ,kp) + mACK(Sc:Ec)*U(Sc:Ec,j ,k) + mAMK(Sc:Ec)*U(Sc:Ec,j ,km) )

                RHS(Sc:Ec,j,k) = RHS(Sc:Ec,j,k)                     &
                            - ( 0.25d0*(v4(Sc:Ec)*dble(jup)*dudy4 + v3(Sc:Ec)*dble(jum)*dudy3)   &
                                -0.5d0*Cmu*invRhoc(Sc:Ec)/dx2(j)*(mu4(Sc:Ec)*dvdx4*dble(jup) - mu3(Sc:Ec)*dvdx3*dble(jum))    )

                RHS(Sc:Ec,j,k) = RHS(Sc:Ec,j,k)                     &
                            - ( 0.25d0*(w6(Sc:Ec)*dudz6(Sc:Ec) + w5(Sc:Ec)*dudz5(Sc:Ec)) &
                                -0.5d0*Cmu*invRhoc(Sc:Ec)/dx3(k)*(mu6(Sc:Ec)*dwdx6(Sc:Ec) - mu5(Sc:Ec)*dwdx5(Sc:Ec))  )
                    
                    ! RHS(Sc:Ec,j,k) =   {viscous} {- Cmp*dpdx}  + {ubc} {- M11Un} {- M12Vn} [- M13Wn]
            
                ! M11Un: Z-direction
                ac(Sc:Ec,j,k) =mACK(Sc:Ec)*dt + 1.d0
                ap(Sc:Ec,j,k) =mAPK(Sc:Ec)*dt
                am(Sc:Ec,j,k) =mAMK(Sc:Ec)*dt
                RHS(Sc:Ec,j,k)=RHS(Sc:Ec,j,k)*dt
            enddo
        enddo
        
        ! TDMA-i*j-K
        call PaScaL_TDMA_plan_many_create(ptdma_plan, (n1sb-1)*(n2sb-1), comm_1d_x3%myrank, comm_1d_x3%nprocs, comm_1d_x3%mpi_comm)
        call PaScaL_TDMA_many_solve_cycle(ptdma_plan, am, ac, ap, RHS,(n1sb-1)*(n2sb-1),(n3sb-1))
        call PaScaL_TDMA_plan_many_destroy(ptdma_plan,comm_1d_x3%nprocs)
        deallocate(am,ac,ap)

            
        ! Trans i-j-K->j-k-i
        allocate(ddU(1:n2sb-1,1:n1sb-1,1:n3sb-1))
        do k = k_indexC, n3sb-1
            ddU(:,:,k)=TRANSPOSE(RHS(:,:,k))
        enddo
        deallocate(RHS) 
        allocate(RHS(1:n2sb-1,1:n3sb-1,1:n1sb-1))
        do k = k_indexC, n3sb-1
            do i = i_indexS, n1sb-1
                RHS(1:n2sb-1,k,i)=ddU(1:n2sb-1,i,k)
            enddo
        enddo
        deallocate(ddU)
        

        ! j-k-i
        
        allocate(am(1:n2sb-1,1:n3sb-1,1:n1sb-1),ac(1:n2sb-1,1:n3sb-1,1:n1sb-1),ap(1:n2sb-1,1:n3sb-1,1:n1sb-1))
        do k = k_indexC, n3sb-1
            do j = j_indexC, n2sb-1
                
                u1(Sc:Ec) = 0.5d0*(U(Sm:Em,j,k) + U(Sc:Ec,j,k))            
                u2(Sc:Ec) = 0.5d0*(U(Sp:Ep,j,k) + U(Sc:Ec,j,k))             

                dudx1(Sc:Ec) = (U(Sc:Ec,j ,k ) - U(Sm:Em,j ,k ))/dx1(Sm:Em)
                dudx2(Sc:Ec) = (U(Sp:Ep,j ,k ) - U(Sc:Ec,j ,k ))/dx1(Sc:Ec)
                
                invRhoc(Sc:Ec) = 0.5d0*(dx1(Sm:Em)*invRho(Sc:Ec,j ,k )+dx1(Sc:Ec)*invRho(Sm:Em,j ,k ) )/dmx1(Sc:Ec)
            
                !X-direction
                mACI(Sc:Ec) = 0.5d0*Cmu*invRhoc(Sc:Ec)/dmx1(Sc:Ec)*(Mu(Sc:Ec,j ,k )/dx1(Sc:Ec) + Mu(Sc:Ec,j ,k )/dx1(Sm:Em))* 2.   &
                            + 0.25d0*(dudx2(Sc:Ec)*0.5d0 + dudx1(Sc:Ec)*0.5d0 - u2(Sc:Ec)/dx1(Sc:Ec) + u1(Sc:Ec)/dx1(Sm:Em) )
                mAPI(Sc:Ec) = -0.5d0*Cmu*invRhoc(Sc:Ec)/dmx1(Sc:Ec)*Mu(Sc:Ec,j ,k )/dx1(Sc:Ec)* 2.  &
                            + 0.25d0*( u2(Sc:Ec)/dx1(Sc:Ec) + dudx2(Sc:Ec)*0.5d0) 
                mAMI(Sc:Ec) = -0.5d0*Cmu*invRhoc(Sc:Ec)/dmx1(Sc:Ec)*Mu(Sc:Ec,j ,k )/dx1(Sm:Em)* 2.  &
                            + 0.25d0*(-u1(Sc:Ec)/dx1(Sm:Em) + dudx1(Sc:Ec)*0.5d0)
                            
                do i = i_indexS, n1sb-1
                    ac(j,k,i)=mACI(i)*dt + 1.d0
                    ap(j,k,i)=mAPI(i)*dt
                    am(j,k,i)=mAMI(i)*dt
                enddo
            enddo
        enddo
        
        ! TDMA-j*K-i        
        call PaScaL_TDMA_plan_many_create(ptdma_plan, (n2sb-1)*(n3sb-1), comm_1d_x1%myrank, comm_1d_x1%nprocs, comm_1d_x1%mpi_comm)
        call PaScaL_TDMA_many_solve_cycle(ptdma_plan, am, ac, ap, RHS,(n2sb-1)*(n3sb-1),(n1sb-1))
        call PaScaL_TDMA_plan_many_destroy(ptdma_plan,comm_1d_x1%nprocs)
        
        deallocate(am,ac,ap)
        
        ! Trans j-k-i->k-i-j 
        allocate(ddU(1:n3sb-1,1:n2sb-1,1:n1sb-1))
        do i = i_indexS, n1sb-1
            ddU(:,:,i)=TRANSPOSE(RHS(:,:,i))
        enddo
        deallocate(RHS) 
        allocate(RHS(1:n3sb-1,1:n1sb-1,1:n2sb-1))
        do j = j_indexC, n2sb-1
            do i = i_indexS, n1sb-1
                RHS(1:n3sb-1,i,j)=ddU(1:n3sb-1,j,i)
            enddo
        enddo
        deallocate(ddU)
        

        ! k-i-j
        
        allocate(am(1:n3sb-1,1:n1sb-1,1:n2sb-1),ac(1:n3sb-1,1:n1sb-1,1:n2sb-1),ap(1:n3sb-1,1:n1sb-1,1:n2sb-1))
        do k = k_indexC, n3sb-1
            do j = j_indexC, n2sb-1
                jm=j-1
                jp=j+1
                jum = jC_BC(jm)
                jup = jC_BC(jp)
                
                v3(Sc:Ec) = 0.5d0*(dx1(Sm:Em)*V(Sc:Ec,j ,k) + dx1(Sc:Ec)*V(Sm:Em,j ,k))/dmx1(Sc:Ec)
                v4(Sc:Ec) = 0.5d0*(dx1(Sm:Em)*V(Sc:Ec,jp,k) + dx1(Sc:Ec)*V(Sm:Em,jp,k))/dmx1(Sc:Ec)
            
                mua(Sc:Ec) = 0.5d0*(dx1(Sm:Em)*Mu(Sc:Ec,jm,k ) + dx1(Sc:Ec)*Mu(Sm:Em,jm,k ))/dmx1(Sc:Ec)
                muc(Sc:Ec) = 0.5d0*(dx1(Sm:Em)*Mu(Sc:Ec,j ,k ) + dx1(Sc:Ec)*Mu(Sm:Em,j ,k ))/dmx1(Sc:Ec)
                mub(Sc:Ec) = 0.5d0*(dx1(Sm:Em)*Mu(Sc:Ec,jp,k ) + dx1(Sc:Ec)*Mu(Sm:Em,jp,k ))/dmx1(Sc:Ec)
                mu3(Sc:Ec) = 0.5d0*(dx2(jm)*muc(Sc:Ec) + dx2(j)*mua(Sc:Ec))/dmx2(j )
                mu4(Sc:Ec) = 0.5d0*(dx2(jp)*muc(Sc:Ec) + dx2(j)*mub(Sc:Ec))/dmx2(jp)
                
                invRhoc(Sc:Ec) = 0.5d0*(dx1(Sm:Em)*invRho(Sc:Ec,j ,k )+dx1(Sc:Ec)*invRho(Sm:Em,j ,k ) )/dmx1(Sc:Ec)
                !Y-direction
                mACJ(Sc:Ec) = 0.5d0*Cmu*invRhoc(Sc:Ec)/dx2(j)*(mu4(Sc:Ec)/dmx2(jp) + mu3(Sc:Ec)/dmx2(j))  &
                            + 0.25d0*(-v4(Sc:Ec)/dmx2(jp) + v3(Sc:Ec)/dmx2(j ))
                mAPJ(Sc:Ec) = -0.5d0*Cmu*invRhoc(Sc:Ec)/dx2(j)*mu4(Sc:Ec)/dmx2(jp)                   &
                            + 0.25d0*( v4(Sc:Ec)/dmx2(jp))
                mAMJ(Sc:Ec) = -0.5d0*Cmu*invRhoc(Sc:Ec)/dx2(j)*mu3(Sc:Ec)/dmx2(j)                        &
                            + 0.25d0*(-v3(Sc:Ec)/dmx2(j))
                mAPJ(Sc:Ec) = mAPJ(Sc:Ec)*dble(jup)
                mAMJ(Sc:Ec) = mAMJ(Sc:Ec)*dble(jum)

                do i = i_indexS, n1sb-1
                    ac(k,i,j)=mACJ(i)*dt + 1.d0
                    ap(k,i,j)=mAPJ(i)*dt
                    am(k,i,j)=mAMJ(i)*dt
                enddo
            enddo
        enddo
        
        ! TDMA-k*i-j
        
        call PaScaL_TDMA_plan_many_create(ptdma_plan, (n3sb-1)*(n1sb-1), comm_1d_x2%myrank, comm_1d_x2%nprocs, comm_1d_x2%mpi_comm)
        call PaScaL_TDMA_many_solve(ptdma_plan, am, ac, ap, RHS,(n3sb-1)*(n1sb-1),(n2sb-1))
        call PaScaL_TDMA_plan_many_destroy(ptdma_plan,comm_1d_x2%nprocs)        
        deallocate(am,ac,ap)

        ! Trans k-i-j->i-j-k
        
        allocate(ddU(1:n1sb-1,1:n3sb-1,1:n2sb-1))
        do j = j_indexC, n2sb-1
            ddU(:,:,j)=TRANSPOSE(RHS(:,:,j))
        enddo
        deallocate(RHS)
        allocate(dU(0:n1sb,0:n2sb,0:n3sb))
        do k = k_indexC, n3sb-1
            do j = j_indexC, n2sb-1
                dU(1:n1sb-1,j,k)=ddU(1:n1sb-1,k,j)
            enddo
        enddo
        deallocate(ddU)
        

        deallocate(u1,u2,u5,u6,v3,v4,w5,w6,Tc)
        deallocate(dudx1,dudx2,dudy3,dudy4,dudz5,dudz6,dvdx3,dvdx4,dwdx5,dwdx6)
        deallocate(mua,mub,muc,mu3,mu4,mu5,mu6)
        deallocate(invRhoc,viscous_u1,viscous_u2,viscous_u3,viscous_u12,viscous_u13)
        deallocate(ubc_up,ubc_down)
        deallocate(mACI,mAPI,mAMI,mACJ,mAPJ,mAMJ,mACK,mAPK,mAMK)

        ! deallocate(RHS)
    end subroutine mpi_momentum_solvedU

    subroutine mpi_momentum_solvedV(T,dx1,dx2,dx3,dmx1,dmx2,dmx3, comm_1d_x1, comm_1d_x2, comm_1d_x3)
        use mpi_topology, only : cart_comm_1d
        use PaScaL_TDMA
        implicit none
        double precision, dimension(0:n1sb,0:n2sb,0:n3sb) ::  T
        double precision :: dx1(0:n1sb),dx2(0:n2sb),dx3(0:n3sb)
        double precision :: dmx1(0:n1sb),dmx2(0:n2sb),dmx3(0:n3sb)
        type(cart_comm_1d), intent(in)  :: comm_1d_x1, comm_1d_x2, comm_1d_x3
        type(ptdma_plan_many)     :: ptdma_plan
            
        integer :: i,j,k
        integer :: im,jm,km
        integer :: ip,jp,kp
        integer :: jvp,jvm

        integer :: Sc,Ec
        integer :: Sm,Em
        integer :: Sp,Ep

        double precision, allocatable, dimension(:,:,:) :: ddV,RHS
        double precision, allocatable, dimension(:,:,:) :: am,ac,ap

        double precision, allocatable, dimension(:) :: u1,u2,v3,v4,w5,w6,Tc
        double precision, allocatable, dimension(:) :: dudy1,dudy2,dvdx1,dvdx2,dvdy3,dvdy4,dvdz5,dvdz6,dwdy5,dwdy6
        double precision, allocatable, dimension(:) :: ddu1,ddu2,ddudy1,ddudy2
        double precision, allocatable, dimension(:) :: mua,muc,mub,mu1,mu2,mu5,mu6
        double precision, allocatable, dimension(:) :: invRhoc,viscous_v1,viscous_v2,viscous_v3,viscous_v21,viscous_v23
        double precision, allocatable, dimension(:) :: vbc_up,vbc_down
        double precision, allocatable, dimension(:) :: mACI,mAPI,mAMI,mACJ,mAPJ,mAMJ,mACK,mAPK,mAMK

        Sc=i_indexC  ;Ec=n1sb-1
        Sm=i_indexC-1;Em=n1sb-1-1
        Sp=i_indexC+1;Ep=n1sb-1+1
        
        allocate(u1(Sc:Ec),u2(Sc:Ec),v3(Sc:Ec),v4(Sc:Ec),w5(Sc:Ec),w6(Sc:Ec),Tc(Sc:Ec))
        allocate(dudy1(Sc:Ec),dudy2(Sc:Ec),dvdx1(Sc:Ec),dvdx2(Sc:Ec),dvdy3(Sc:Ec),dvdy4(Sc:Ec),dvdz5(Sc:Ec),dvdz6(Sc:Ec),dwdy5(Sc:Ec),dwdy6(Sc:Ec))
        allocate(ddu1(Sc:Ec),ddu2(Sc:Ec),ddudy1(Sc:Ec),ddudy2(Sc:Ec))
        allocate(mua(Sc:Ec),muc(Sc:Ec),mub(Sc:Ec),mu1(Sc:Ec),mu2(Sc:Ec),mu5(Sc:Ec),mu6(Sc:Ec))
        
        allocate(invRhoc(Sc:Ec),viscous_v1(Sc:Ec),viscous_v2(Sc:Ec),viscous_v3(Sc:Ec),viscous_v21(Sc:Ec),viscous_v23(Sc:Ec))
        allocate(vbc_up(Sc:Ec),vbc_down(Sc:Ec))
        allocate(mACI(Sc:Ec),mAPI(Sc:Ec),mAMI(Sc:Ec),mACJ(Sc:Ec),mAPJ(Sc:Ec),mAMJ(Sc:Ec),mACK(Sc:Ec),mAPK(Sc:Ec),mAMK(Sc:Ec))

        ! i-j-k:RHS-ddV
        ! i-j-k
        
        allocate(RHS(1:n1sb-1,1:n2sb-1,1:n3sb-1))
        allocate(am(1:n1sb-1,1:n2sb-1,1:n3sb-1))
        allocate(ac(1:n1sb-1,1:n2sb-1,1:n3sb-1))
        allocate(ap(1:n1sb-1,1:n2sb-1,1:n3sb-1))
        do k = 1, n3sb-1
            kp = k+1
            km = k-1
            do j =  j_indexS, n2sb-1
            jp = j + 1
            jm = j - 1
            jvp = jS_BC(jm)
            jvm = jS_BC(jp)
                
                v3(Sc:Ec) = 0.5d0*(V(Sc:Ec,jm,k) + V(Sc:Ec,j ,k))
                v4(Sc:Ec) = 0.5d0*(V(Sc:Ec,j ,k) + V(Sc:Ec,jp,k))

                u1(Sc:Ec) = (dx2(jm)*U(Sc:Ec, j,k) + dx2(j)*U(Sc:Ec, jm,k))/dmx2(j)*0.5d0
                u2(Sc:Ec) = (dx2(jm)*U(Sp:Ep,j,k) + dx2(j)*U(Sp:Ep,jm,k))/dmx2(j)*0.5d0

                w5(Sc:Ec) = (dx2(jm)*W(Sc:Ec,j,k ) + dx2(j)*W(Sc:Ec,jm,k ))/dmx2(j)*0.5d0
                w6(Sc:Ec) = (dx2(jm)*W(Sc:Ec,j,kp) + dx2(j)*W(Sc:Ec,jm,kp))/dmx2(j)*0.5d0
                
                ddu1(Sc:Ec) = (dx2(jm)*dU(Sc:Ec,j,k) + dx2(j)*dU(Sc:Ec,jm,k))/dmx2(j)*0.5d0
                ddu2(Sc:Ec) = (dx2(jm)*dU(Sp:Ep,j,k) + dx2(j)*dU(Sp:Ep,jm,k))/dmx2(j)*0.5d0

                
                dvdx1(Sc:Ec) = (V(Sc:Ec, j,k) - V(Sm:Em,j,k))/dmx1(Sc:Ec)
                dvdx2(Sc:Ec) = (V(Sp:Ep,j,k) - V(Sc:Ec, j,k))/dmx1(Sp:Ep)
                dvdy3(Sc:Ec) = (V(Sc:Ec,j ,k) - V(Sc:Ec,jm,k))/dx2(jm)
                dvdy4(Sc:Ec) = (V(Sc:Ec,jp,k) - V(Sc:Ec,j ,k))/dx2(j )
                dvdz5(Sc:Ec) = (V(Sc:Ec,j,k ) - V(Sc:Ec,j,km))/dmx3(k )
                dvdz6(Sc:Ec) = (V(Sc:Ec,j,kp) - V(Sc:Ec,j,k ))/dmx3(kp)
                

                dudy1(Sc:Ec) = (U(Sc:Ec, j,k) - U(Sc:Ec, jm,k))/dmx2(j)
                dudy2(Sc:Ec) = (U(Sp:Ep,j,k) - U(Sp:Ep,jm,k))/dmx2(j)

                dwdy5(Sc:Ec) = (W(Sc:Ec,j,k ) - W(Sc:Ec,jm,k ))/dmx2(j)
                dwdy6(Sc:Ec) = (W(Sc:Ec,j,kp) - W(Sc:Ec,jm,kp))/dmx2(j)

                dvdx1(Sc:Ec) = (V(Sc:Ec, j,k) - V(Sm:Em,j,k))/dmx1(Sc:Ec)
                dvdx2(Sc:Ec) = (V(Sp:Ep,j,k) - V(Sc:Ec, j,k))/dmx1(Sp:Ep)
                ddudy1(Sc:Ec) = (dU(Sc:Ec, j,k) - dU(Sc:Ec, jm,k))/dmx2(j)
                ddudy2(Sc:Ec) = (dU(Sp:Ep,j,k) - dU(Sp:Ep,jm,k))/dmx2(j) 

                mua(Sc:Ec) = 0.5d0*(dx2(jm)*Mu(Sm:Em,j ,k ) + dx2(j)*Mu(Sm:Em,jm,k ))/dmx2(j)
                muc(Sc:Ec) = 0.5d0*(dx2(jm)*Mu(Sc:Ec,j ,k ) + dx2(j)*Mu(Sc:Ec,jm,k ))/dmx2(j)
                mub(Sc:Ec) = 0.5d0*(dx2(jm)*Mu(Sp:Ep,j ,k ) + dx2(j)*Mu(Sp:Ep,jm,k ))/dmx2(j)
                mu1(Sc:Ec) = 0.5d0*(dx1(Sm:Em)*muc(Sc:Ec) + dx1(Sc:Ec)*mua(Sc:Ec))/dmx1(Sc:Ec)
                mu2(Sc:Ec) = 0.5d0*(dx1(Sp:Ep)*muc(Sc:Ec) + dx1(Sc:Ec)*mub(Sc:Ec))/dmx1(Sp:Ep)

                mua(Sc:Ec) = 0.5d0*(dx2(jm)*Mu(Sc:Ec,j ,km) + dx2(j)*Mu(Sc:Ec,jm,km))/dmx2(j)
                muc(Sc:Ec) = 0.5d0*(dx2(jm)*Mu(Sc:Ec,j ,k ) + dx2(j)*Mu(Sc:Ec,jm,k ))/dmx2(j)
                mub(Sc:Ec) = 0.5d0*(dx2(jm)*Mu(Sc:Ec,j ,kp) + dx2(j)*Mu(Sc:Ec,jm,kp))/dmx2(j)
                mu5(Sc:Ec) = 0.5d0*(dx3(km)*muc(Sc:Ec) + dx3(k)*mua(Sc:Ec))/dmx3(k )
                mu6(Sc:Ec) = 0.5d0*(dx3(kp)*muc(Sc:Ec) + dx3(k)*mub(Sc:Ec))/dmx3(kp)

                invRhoc(Sc:Ec) = 0.5d0*(dx2(jm)*invRho(Sc:Ec,j ,k ) + dx2(j )*invRho(Sc:Ec,jm,k ))/dmx2(j)

                !---DIFFUSION TERM
                viscous_v1(Sc:Ec) = 1.d0*(mu2(Sc:Ec)*dvdx2(Sc:Ec) - mu1(Sc:Ec)*dvdx1(Sc:Ec))/dx1(Sc:Ec)
                viscous_v2(Sc:Ec) = 1.d0*(Mu(Sc:Ec,j ,k )*dvdy4(Sc:Ec) - Mu(Sc:Ec,jm,k )*dvdy3(Sc:Ec))/dmx2(j)
                viscous_v3(Sc:Ec) = 1.d0*(mu6(Sc:Ec)*dvdz6(Sc:Ec) - mu5(Sc:Ec)*dvdz5(Sc:Ec))/dx3(k )
                viscous_v21(Sc:Ec) = 1.d0*(mu2(Sc:Ec)*dudy2(Sc:Ec) - mu1(Sc:Ec)*dudy1(Sc:Ec))/dx1(Sc:Ec)
                viscous_v23(Sc:Ec) = 1.d0*(mu6(Sc:Ec)*dwdy6(Sc:Ec) - mu5(Sc:Ec)*dwdy5(Sc:Ec))/dx3(k)

                RHS(Sc:Ec,j,k) = 0.5d0*Cmu*invRhoc(Sc:Ec)*(viscous_v1(Sc:Ec) + 2.*viscous_v2(Sc:Ec) + viscous_v3(Sc:Ec) + viscous_v21(Sc:Ec) + viscous_v23(Sc:Ec))          ! GE

                !---PRESSURE TERM (invRhoc(Sc:Ec) is included for the NOB case)
                RHS(Sc:Ec,j,k) = RHS(Sc:Ec,j,k)     &
                            - Cmp*invRhoc(Sc:Ec)*(P(Sc:Ec,j,k) - P(Sc:Ec,jm,k))/dmx2(j)

                !---!BUOYANCY TERM (we use Cmt*theta for the buoyancy term)
                ! Tc(Sc:Ec) = 0.5d0*(dx2(jm)*T(Sc:Ec,j,k) + dx2(j)*T(Sc:Ec,jm,k))/dmx2(j)
                ! !THETAy = Cmt*((Tc(Sc:Ec) - Thetam))**b
                ! !THETAy = Cmt*(Tc(Sc:Ec))
                ! RHS(Sc:Ec,j,k) = RHS(Sc:Ec,j,k)     &
                !             + Cmt*(Tc(Sc:Ec) + a12pera11*Tc(Sc:Ec)**2.*DeltaT)*invRhoc(Sc:Ec)
        
                !---vbc
                vbc_down(Sc:Ec) = ( 0.25d0/dx2(jm)*v3(Sc:Ec)*VBCbt_sub(Sc:Ec,k)    &
                                - 0.25d0*0.5d0*dvdy3(Sc:Ec)*VBCbt_sub(Sc:Ec,k)   )

                vbc_up(Sc:Ec) = (-0.25d0/dx2(j )*v4(Sc:Ec)*VBCup_sub(Sc:Ec,k)  &
                            - 0.25d0*0.5d0*dvdy4(Sc:Ec)*VBCup_sub(Sc:Ec,k) )

                vbc_down(Sc:Ec) = vbc_down(Sc:Ec) + 0.5d0*Cmu*invRhoc(Sc:Ec)/dmx2(j)*Mu(Sc:Ec,jm,k )/dx2(jm)*VBCbt_sub(Sc:Ec,k)*2.           !'*2.' is needed since the modified GE
                vbc_up(Sc:Ec)   = vbc_up(Sc:Ec)   + 0.5d0*Cmu*invRhoc(Sc:Ec)/dmx2(j)*Mu(Sc:Ec,j ,k )/dx2(j )*VBCup_sub(Sc:Ec,k)*2.

                RHS(Sc:Ec,j,k) = RHS(Sc:Ec,j,k)     &
                            + dble(1.d0 - jvm)*vbc_down(Sc:Ec)    &
                            + dble(1.d0 - jvp)*vbc_up(Sc:Ec)

                !----M22Vn
                mACI(Sc:Ec) =  0.5d0*Cmu*invRhoc(Sc:Ec)/dx1(Sc:Ec)*(mu2(Sc:Ec)/dmx1(Sp:Ep) + mu1(Sc:Ec)/dmx1(Sc:Ec))       &
                            + 0.25d0*(-u2(Sc:Ec)/dmx1(Sp:Ep)+ u1(Sc:Ec)/dmx1(Sc:Ec))
                mAPI(Sc:Ec) = -0.5d0*Cmu*invRhoc(Sc:Ec)/dx1(Sc:Ec)*mu2(Sc:Ec)/dmx1(Sp:Ep)   &
                            + 0.25d0*( u2(Sc:Ec)/dmx1(Sp:Ep))
                mAMI(Sc:Ec) = -0.5d0*Cmu*invRhoc(Sc:Ec)/dx1(Sc:Ec)*mu1(Sc:Ec)/dmx1(Sc:Ec)   &
                            + 0.25d0*(-u1(Sc:Ec)/dmx1(Sc:Ec))

                mACJ(Sc:Ec) = 0.5d0*Cmu*invRhoc(Sc:Ec)/dmx2(j)*(Mu(Sc:Ec,j ,k )/dx2(j) + Mu(Sc:Ec,jm,k )/dx2(jm))*2.   &
                            + 0.25d0*(dvdy4(Sc:Ec)*0.5d0 + dvdy3(Sc:Ec)*0.5d0-v4(Sc:Ec)/dx2(j) + v3(Sc:Ec)/dx2(jm) )
                mAPJ(Sc:Ec) = -0.5d0*Cmu*invRhoc(Sc:Ec)/dmx2(j)*Mu(Sc:Ec,j ,k )/dx2(j )*2.   &
                            + 0.25d0*( v4(Sc:Ec)/dx2(j ) + dvdy4(Sc:Ec)*0.5d0)
                mAPJ(Sc:Ec) = mAPJ(Sc:Ec)*dble(jvp)
                mAMJ(Sc:Ec) = -0.5d0*Cmu*invRhoc(Sc:Ec)/dmx2(j)*Mu(Sc:Ec,jm,k )/dx2(jm)*2.   &
                            + 0.25d0*(-v3(Sc:Ec)/dx2(jm) + dvdy3(Sc:Ec)*0.5d0)
                mAMJ(Sc:Ec) = mAMJ(Sc:Ec)*dble(jvm)

                mACK(Sc:Ec) = 0.5d0*Cmu*invRhoc(Sc:Ec)/dx3(k)*(mu6(Sc:Ec)/dmx3(kp) + mu5(Sc:Ec)/dmx3(k))    &
                            + 0.25d0*(-w6(Sc:Ec)/dmx3(kp) + w5(Sc:Ec)/dmx3(k))
                mAPK(Sc:Ec) = -0.5d0*Cmu*invRhoc(Sc:Ec)/dx3(k)*mu6(Sc:Ec)/dmx3(kp)  &
                            + 0.25d0*( w6(Sc:Ec)/dmx3(kp))
                mAMK(Sc:Ec) = -0.5d0*Cmu*invRhoc(Sc:Ec)/dx3(k)*mu5(Sc:Ec)/dmx3(k )  &
                            + 0.25d0*(-w5(Sc:Ec)/dmx3(k ))

                RHS(Sc:Ec,j,k) = RHS(Sc:Ec,j,k)     &
                            -(mAPI(Sc:Ec)*V(Sp:Ep,j, k ) + mACI(Sc:Ec)*V(Sc:Ec,j,k) + mAMI(Sc:Ec)*V(Sm:Em,j, k )      &
                            + mAPJ(Sc:Ec)*V(Sc:Ec, jp,k ) + mACJ(Sc:Ec)*V(Sc:Ec,j,k) + mAMJ(Sc:Ec)*V(Sc:Ec, jm,k )      &
                            + mAPK(Sc:Ec)*V(Sc:Ec, j, kp) + mACK(Sc:Ec)*V(Sc:Ec,j,k) + mAMK(Sc:Ec)*V(Sc:Ec, j, km)   )

                !----M21Un 
                RHS(Sc:Ec,j,k) = RHS(Sc:Ec,j,k)     &
                            -(0.25d0*(u2(Sc:Ec)*dvdx2(Sc:Ec) + u1(Sc:Ec)*dvdx1(Sc:Ec))    &
                            - 0.5d0*Cmu*invRhoc(Sc:Ec)/dx1(Sc:Ec)*(mu2(Sc:Ec)*dudy2(Sc:Ec) - mu1(Sc:Ec)*dudy1(Sc:Ec)) )
                    
                !----M23Wn
                RHS(Sc:Ec,j,k) = RHS(Sc:Ec,j,k)     &
                            -(0.25d0*(w6(Sc:Ec)*dvdz6(Sc:Ec) + w5(Sc:Ec)*dvdz5(Sc:Ec))       &
                            - 0.5d0*Cmu*invRhoc(Sc:Ec)/dx3(k)*(mu6(Sc:Ec)*dwdy6(Sc:Ec) - mu5(Sc:Ec)*dwdy5(Sc:Ec)))
            
                !----M21 ddU
                RHS(Sc:Ec,j,k) = RHS(Sc:Ec,j,k)     &
                            -(0.25d0*(ddu2(Sc:Ec)*dvdx2(Sc:Ec)+ddu1(Sc:Ec)*dvdx1(Sc:Ec))     &
                            - 0.5d0*Cmu*invRhoc(Sc:Ec)/dx1(Sc:Ec)*(mu2(Sc:Ec)*ddudy2(Sc:Ec) - mu1(Sc:Ec)*ddudy1(Sc:Ec)) )

                ! M22Vn: Z-direction
                ac(Sc:Ec,j,k) =mACK(Sc:Ec)*dt + 1.d0
                ap(Sc:Ec,j,k) =mAPK(Sc:Ec)*dt
                am(Sc:Ec,j,k) =mAMK(Sc:Ec)*dt
                RHS(Sc:Ec,j,k)=RHS(Sc:Ec,j,k)*dt

            end do
        end do
        
        ! TDMA-i*j-K
        
        call PaScaL_TDMA_plan_many_create(ptdma_plan, (n1sb-1)*(n2sb-1), comm_1d_x3%myrank, comm_1d_x3%nprocs, comm_1d_x3%mpi_comm)
        call PaScaL_TDMA_many_solve_cycle(ptdma_plan, am, ac, ap, RHS,(n1sb-1)*(n2sb-1),(n3sb-1))
        call PaScaL_TDMA_plan_many_destroy(ptdma_plan,comm_1d_x3%nprocs)
        
        deallocate(am,ac,ap)

        ! Trans i-j-K->j-k-i
        
        allocate(ddV(1:n2sb-1,1:n1sb-1,1:n3sb-1))
        do k = k_indexC, n3sb-1
            ddV(:,:,k)=TRANSPOSE(RHS(:,:,k))
        enddo
        deallocate(RHS) 
        allocate(RHS(1:n2sb-1,1:n3sb-1,1:n1sb-1))
        do k = k_indexC, n3sb-1
            do i = i_indexC, n1sb-1
                RHS(1:n2sb-1,k,i)=ddV(1:n2sb-1,i,k)
            enddo
        enddo
        deallocate(ddV)
        

        ! j-k-i
        
        allocate(am(1:n2sb-1,1:n3sb-1,1:n1sb-1),ac(1:n2sb-1,1:n3sb-1,1:n1sb-1),ap(1:n2sb-1,1:n3sb-1,1:n1sb-1))
        do k = k_indexC, n3sb-1
            do j = j_indexS, n2sb-1
                jp = j + 1
                jm = j - 1
                jvp = jS_BC(jm)
                jvm = jS_BC(jp)
                    
                u1(Sc:Ec) = (dx2(jm)*U(Sc:Ec, j,k) + dx2(j)*U(Sc:Ec, jm,k))/dmx2(j)*0.5d0
                u2(Sc:Ec) = (dx2(jm)*U(Sp:Ep,j,k) + dx2(j)*U(Sp:Ep,jm,k))/dmx2(j)*0.5d0
            
                invRhoc(Sc:Ec) = 0.5d0*(dx2(jm)*invRho(Sc:Ec,j ,k ) + dx2(j )*invRho(Sc:Ec,jm,k ))/dmx2(j)

                ! M22Vn:X-Direction
                mACI(Sc:Ec) =  0.5d0*Cmu*invRhoc(Sc:Ec)/dx1(Sc:Ec)*(mu2(Sc:Ec)/dmx1(Sp:Ep) + mu1(Sc:Ec)/dmx1(Sc:Ec))       &
                            + 0.25d0*(-u2(Sc:Ec)/dmx1(Sp:Ep)+ u1(Sc:Ec)/dmx1(Sc:Ec))
                mAPI(Sc:Ec) = -0.5d0*Cmu*invRhoc(Sc:Ec)/dx1(Sc:Ec)*mu2(Sc:Ec)/dmx1(Sp:Ep)   &
                            + 0.25d0*( u2(Sc:Ec)/dmx1(Sp:Ep))
                mAMI(Sc:Ec) = -0.5d0*Cmu*invRhoc(Sc:Ec)/dx1(Sc:Ec)*mu1(Sc:Ec)/dmx1(Sc:Ec)   &
                            + 0.25d0*(-u1(Sc:Ec)/dmx1(Sc:Ec))
                            
                do i = i_indexC, n1sb-1
                    ac(j,k,i)=mACI(i)*dt + 1.d0
                    ap(j,k,i)=mAPI(i)*dt
                    am(j,k,i)=mAMI(i)*dt
                enddo
            enddo
        enddo
        
        
        ! TDMA-j*K-i
        
        call PaScaL_TDMA_plan_many_create(ptdma_plan, (n2sb-1)*(n3sb-1), comm_1d_x1%myrank, comm_1d_x1%nprocs, comm_1d_x1%mpi_comm)
        call PaScaL_TDMA_many_solve_cycle(ptdma_plan, am, ac, ap, RHS,(n2sb-1)*(n3sb-1),(n1sb-1))
        call PaScaL_TDMA_plan_many_destroy(ptdma_plan,comm_1d_x1%nprocs)
        
        deallocate(am,ac,ap)
        
        ! Trans j-k-i->k-i-j 
        
        allocate(ddV(1:n3sb-1,1:n2sb-1,1:n1sb-1))
        do i = i_indexC, n1sb-1
            ddV(:,:,i)=TRANSPOSE(RHS(:,:,i))
        enddo
        deallocate(RHS) 
        allocate(RHS(1:n3sb-1,1:n1sb-1,1:n2sb-j_indexS))
        do j = j_indexS, n2sb-1
            do i = i_indexC, n1sb-1
                RHS(1:n3sb-1,i,j-j_indexS+1)=ddV(1:n3sb-1,j,i)
            enddo
        enddo
        deallocate(ddV)
        

        ! k-i-j
        
        allocate(am(1:n3sb-1,1:n1sb-1,1:n2sb-j_indexS),ac(1:n3sb-1,1:n1sb-1,1:n2sb-j_indexS),ap(1:n3sb-1,1:n1sb-1,1:n2sb-j_indexS))
        do k = k_indexC, n3sb-1
            do j = j_indexS, n2sb-1
                jm=j-1
                jp=j+1
                jvp = jS_BC(jm)
                jvm = jS_BC(jp)

                v3(Sc:Ec) = 0.5d0*(V(Sc:Ec,jm,k) + V(Sc:Ec,j ,k))
                v4(Sc:Ec) = 0.5d0*(V(Sc:Ec,j ,k) + V(Sc:Ec,jp,k))

                dvdy3(Sc:Ec) = (V(Sc:Ec,j ,k) - V(Sc:Ec,jm,k))/dx2(jm)
                dvdy4(Sc:Ec) = (V(Sc:Ec,jp,k) - V(Sc:Ec,j ,k))/dx2(j )

                invRhoc(Sc:Ec) = 0.5d0*(dx2(jm)*invRho(Sc:Ec,j ,k ) + dx2(j )*invRho(Sc:Ec,jm,k ))/dmx2(j)

                ! M22Vn:Y-DIRECTION
                mACJ(Sc:Ec) = 0.5d0*Cmu*invRhoc(Sc:Ec)/dmx2(j)*(Mu(Sc:Ec,j ,k )/dx2(j) + Mu(Sc:Ec,jm,k )/dx2(jm))*2.   &
                            + 0.25d0*(dvdy4(Sc:Ec)*0.5d0 + dvdy3(Sc:Ec)*0.5d0-v4(Sc:Ec)/dx2(j) + v3(Sc:Ec)/dx2(jm) )
                mAPJ(Sc:Ec) = -0.5d0*Cmu*invRhoc(Sc:Ec)/dmx2(j)*Mu(Sc:Ec,j ,k )/dx2(j )*2.   &
                            + 0.25d0*( v4(Sc:Ec)/dx2(j ) + dvdy4(Sc:Ec)*0.5d0)
                mAPJ(Sc:Ec) = mAPJ(Sc:Ec)*dble(jvp)
                mAMJ(Sc:Ec) = -0.5d0*Cmu*invRhoc(Sc:Ec)/dmx2(j)*Mu(Sc:Ec,jm,k )/dx2(jm)*2.   &
                            + 0.25d0*(-v3(Sc:Ec)/dx2(jm) + dvdy3(Sc:Ec)*0.5d0)
                mAMJ(Sc:Ec) = mAMJ(Sc:Ec)*dble(jvm)

                do i = i_indexC, n1sb-1
                    ac(k,i,j-j_indexS+1)=mACJ(i)*dt + 1.d0
                    ap(k,i,j-j_indexS+1)=mAPJ(i)*dt
                    am(k,i,j-j_indexS+1)=mAMJ(i)*dt
                enddo
            enddo
        enddo
        

        ! TDMA-k*i-j
        
        call PaScaL_TDMA_plan_many_create(ptdma_plan, (n3sb-1)*(n1sb-1), comm_1d_x2%myrank, comm_1d_x2%nprocs, comm_1d_x2%mpi_comm)
        call PaScaL_TDMA_many_solve(ptdma_plan, am, ac, ap, RHS,(n3sb-1)*(n1sb-1),(n2sb-j_indexS))
        call PaScaL_TDMA_plan_many_destroy(ptdma_plan,comm_1d_x2%nprocs)
        
        deallocate(am,ac,ap)

        ! Trans k-i-j->i-j-k
        
        allocate(ddV(1:n1sb-1,1:n3sb-1,1:n2sb-1))
        do j = j_indexS, n2sb-1
            ddV(:,:,j)=TRANSPOSE(RHS(:,:,j-j_indexS+1))
        enddo
        deallocate(RHS)
        allocate(dV(0:n1sb,0:n2sb,0:n3sb))
        dV(0:n1sb,0:1,0:n3sb)=0.d0;dV(0:n1sb,n2sb,0:n3sb)=0.d0
        do k = k_indexC, n3sb-1
            do j = j_indexS, n2sb-1
                dV(1:n1sb-1,j,k)=ddV(1:n1sb-1,k,j)
            enddo
        enddo
        deallocate(ddV)
        

        deallocate(u1,u2,v3,v4,w5,w6,Tc)
        deallocate(dudy1,dudy2,dvdx1,dvdx2,dvdy3,dvdy4,dvdz5,dvdz6,dwdy5,dwdy6)
        deallocate(ddu1,ddu2,ddudy1,ddudy2)
        deallocate(mua,muc,mub,mu1,mu2,mu5,mu6)
        deallocate(invRhoc,viscous_v1,viscous_v2,viscous_v3,viscous_v21,viscous_v23)
        deallocate(vbc_up,vbc_down)
        deallocate(mACI,mAPI,mAMI,mACJ,mAPJ,mAMJ,mACK,mAPK,mAMK)

    end subroutine mpi_momentum_solvedV
    
    subroutine mpi_momentum_solvedW(T,dx1,dx2,dx3,dmx1,dmx2,dmx3, comm_1d_x1, comm_1d_x2, comm_1d_x3)
        use mpi_topology, only : cart_comm_1d
        use PaScaL_TDMA
        implicit none
        double precision, dimension(0:n1sb,0:n2sb,0:n3sb) ::  T
        double precision :: dx1(0:n1sb),dx2(0:n2sb),dx3(0:n3sb)
        double precision :: dmx1(0:n1sb),dmx2(0:n2sb),dmx3(0:n3sb)
        type(cart_comm_1d), intent(in)  :: comm_1d_x1, comm_1d_x2, comm_1d_x3
        type(ptdma_plan_many)     :: ptdma_plan
            
        integer :: i,j,k
        integer :: im,jm,km
        integer :: ip,jp,kp
        integer :: jwp,jwm

        integer :: Sc,Ec
        integer :: Sm,Em
        integer :: Sp,Ep

        double precision, allocatable, dimension(:,:,:) :: ddW,RHS
        double precision, allocatable, dimension(:,:,:) :: am,ac,ap

        double precision, allocatable, dimension(:) :: u1,u2,v3,v4,w6,w5
        double precision, allocatable, dimension(:) :: dudz1,dudz2,dvdz3,dvdz4,dwdx1,dwdx2,dwdy3,dwdy4,dwdz5,dwdz6
        double precision, allocatable, dimension(:) :: ddu1,ddu2,ddudz1,ddudz2
        double precision, allocatable, dimension(:) :: ddv3,ddv4,ddvdz3,ddvdz4
        double precision, allocatable, dimension(:) :: mua,muc,mub,mu1,mu2,mu3,mu4,invRhoc
        double precision, allocatable, dimension(:) :: viscous_w1,viscous_w2,viscous_w3,viscous_w31,viscous_w32
        double precision, allocatable, dimension(:) :: zbc_up,zbc_down
        double precision, allocatable, dimension(:) :: mACI,mAPI,mAMI,mACJ,mAPJ,mAMJ,mACK,mAPK,mAMK

        Sc=i_indexC  ;Ec=n1sb-1
        Sm=i_indexC-1;Em=n1sb-1-1
        Sp=i_indexC+1;Ep=n1sb-1+1

        allocate(u1(Sc:Ec),u2(Sc:Ec),v3(Sc:Ec),v4(Sc:Ec),w6(Sc:Ec),w5(Sc:Ec))
        allocate(dudz1(Sc:Ec),dudz2(Sc:Ec),dvdz3(Sc:Ec),dvdz4(Sc:Ec))
        allocate(dwdx1(Sc:Ec),dwdx2(Sc:Ec),dwdy3(Sc:Ec),dwdy4(Sc:Ec),dwdz5(Sc:Ec),dwdz6(Sc:Ec))
        allocate(ddu1(Sc:Ec),ddu2(Sc:Ec),ddudz1(Sc:Ec),ddudz2(Sc:Ec))
        allocate(ddv3(Sc:Ec),ddv4(Sc:Ec),ddvdz3(Sc:Ec),ddvdz4(Sc:Ec))
        allocate(mua(Sc:Ec),muc(Sc:Ec),mub(Sc:Ec),mu1(Sc:Ec),mu2(Sc:Ec),mu3(Sc:Ec),mu4(Sc:Ec),invRhoc(Sc:Ec))
        allocate(viscous_w1(Sc:Ec),viscous_w2(Sc:Ec),viscous_w3(Sc:Ec),viscous_w31(Sc:Ec),viscous_w32(Sc:Ec))
        allocate(zbc_up(Sc:Ec),zbc_down(Sc:Ec))
        allocate(mACI(Sc:Ec),mAPI(Sc:Ec),mAMI(Sc:Ec),mACJ(Sc:Ec),mAPJ(Sc:Ec),mAMJ(Sc:Ec),mACK(Sc:Ec),mAPK(Sc:Ec),mAMK(Sc:Ec))
        
        ! i-j-k:RHS-ddV
        ! i-j-k
        
        allocate(RHS(1:n1sb-1,1:n2sb-1,1:n3sb-1))
        allocate(am(1:n1sb-1,1:n2sb-1,1:n3sb-1))
        allocate(ac(1:n1sb-1,1:n2sb-1,1:n3sb-1))
        allocate(ap(1:n1sb-1,1:n2sb-1,1:n3sb-1))
        do k = k_indexS, n3sb-1
            kp = k+1
            km = k-1
            do j = j_indexC, n2sb-1
                jp = j + 1
                jm = j - 1
                jwp = jC_BC(jp)
                jwm = jC_BC(jm)

                w5(Sc:Ec) = 0.5d0*(W(Sc:Ec,j,k) + W(Sc:Ec,j,km))
                w6(Sc:Ec) = 0.5d0*(W(Sc:Ec,j,k) + W(Sc:Ec,j,kp))
                
                u1(Sc:Ec) = (dx3(km)*U(Sc:Ec,j,k) + dx3(k)*U(Sc:Ec, j,km))/dmx3(k)*0.5d0
                u2(Sc:Ec) = (dx3(km)*U(Sp:Ep,j,k) + dx3(k)*U(Sp:Ep,j,km))/dmx3(k)*0.5d0
                v3(Sc:Ec) = (dx3(km)*V(Sc:Ec,j, k) + dx3(k)*V(Sc:Ec,j, km))/dmx3(k)*0.5d0
                v4(Sc:Ec) = (dx3(km)*V(Sc:Ec,jp,k) + dx3(k)*V(Sc:Ec,jp,km))/dmx3(k)*0.5d0
                
                dwdx1(Sc:Ec) = (W(Sc:Ec, j,k) - W(Sm:Em,j,k))/dmx1(Sc:Ec)
                dwdx2(Sc:Ec) = (W(Sp:Ep,j,k) - W(Sc:Ec, j,k))/dmx1(Sp:Ep)
                dwdy3(Sc:Ec) = (W(Sc:Ec,j, k) - W(Sc:Ec,jm,k))/dmx2(j )
                dwdy4(Sc:Ec) = (W(Sc:Ec,jp,k) - W(Sc:Ec,j, k))/dmx2(jp)
                dwdz5(Sc:Ec) = (W(Sc:Ec,j,k ) - W(Sc:Ec,j,km))/dx3(km)
                dwdz6(Sc:Ec) = (W(Sc:Ec,j,kp) - W(Sc:Ec,j,k ))/dx3(k )
                
                dudz1(Sc:Ec) = (U(Sc:Ec,j,k) - U(Sc:Ec,j,km))/dmx3(k)
                dudz2(Sc:Ec) = (U(Sp:Ep,j,k) - U(Sp:Ep,j,km))/dmx3(k)
                                    
                dvdz3(Sc:Ec) = (V(Sc:Ec,j ,k) - V(Sc:Ec,j ,km))/dmx3(k)
                dvdz4(Sc:Ec) = (V(Sc:Ec,jp,k) - V(Sc:Ec,jp,km))/dmx3(k)
                
                ddudz1(Sc:Ec) = (dU(Sc:Ec,j,k) - dU(Sc:Ec,j,km))/dmx3(k)
                ddudz2(Sc:Ec) = (dU(Sp:Ep,j,k) - dU(Sp:Ep,j,km))/dmx3(k) 
                ddvdz3(Sc:Ec) = (dV(Sc:Ec,j ,k) - dV(Sc:Ec,j ,km))/dmx3(k)
                ddvdz4(Sc:Ec) = (dV(Sc:Ec,jp,k) - dV(Sc:Ec,jp,km))/dmx3(k)
                
                ddu1(Sc:Ec) = (dx3(km)*dU(Sc:Ec,j,k) + dx3(k)*dU(Sc:Ec, j,km))/dmx3(k)*0.5d0
                ddu2(Sc:Ec) = (dx3(km)*dU(Sp:Ep,j,k) + dx3(k)*dU(Sp:Ep,j,km))/dmx3(k)*0.5d0
                ddv3(Sc:Ec) = (dx3(km)*dV(Sc:Ec,j, k) + dx3(k)*dV(Sc:Ec,j, km))/dmx3(k)*0.5d0*dble(jwm)
                ddv4(Sc:Ec) = (dx3(km)*dV(Sc:Ec,jp,k) + dx3(k)*dV(Sc:Ec,jp,km))/dmx3(k)*0.5d0*dble(jwp)
                
                mua(Sc:Ec) = 0.5d0*(dx3(km)*Mu(Sm:Em,j ,k ) + dx3(k)*Mu(Sm:Em,j ,km))/dmx3(k)
                muc(Sc:Ec) = 0.5d0*(dx3(km)*Mu(Sc:Ec,j ,k ) + dx3(k)*Mu(Sc:Ec,j ,km))/dmx3(k)
                mub(Sc:Ec) = 0.5d0*(dx3(km)*Mu(Sp:Ep,j ,k ) + dx3(k)*Mu(Sp:Ep,j ,km))/dmx3(k)
                mu1(Sc:Ec) = 0.5d0*(dx1(Sm:Em)*muc(Sc:Ec) + dx1(Sc:Ec)*mua(Sc:Ec))/dmx1(Sc:Ec)
                mu2(Sc:Ec) = 0.5d0*(dx1(Sp:Ep)*muc(Sc:Ec) + dx1(Sc:Ec)*mub(Sc:Ec))/dmx1(Sp:Ep)
                
                mua(Sc:Ec) = 0.5d0*(dx3(km)*Mu(Sc:Ec,jm,k ) + dx3(k)*Mu(Sc:Ec,jm,km))/dmx3(k)
                muc(Sc:Ec) = 0.5d0*(dx3(km)*Mu(Sc:Ec,j ,k ) + dx3(k)*Mu(Sc:Ec,j ,km))/dmx3(k)
                mub(Sc:Ec) = 0.5d0*(dx3(km)*Mu(Sc:Ec,jp,k ) + dx3(k)*Mu(Sc:Ec,jp,km))/dmx3(k)
                mu3(Sc:Ec) = 0.5d0*(dx2(jm)*muc(Sc:Ec) + dx2(j)*mua(Sc:Ec))/dmx2(j )
                mu4(Sc:Ec) = 0.5d0*(dx2(jp)*muc(Sc:Ec) + dx2(j)*mub(Sc:Ec))/dmx2(jp)
                
                invRhoc(Sc:Ec) = 0.5d0*(dx3(km)*invRho(Sc:Ec,j ,k ) + dx3(k )*invRho(Sc:Ec,j ,km) )/dmx3(k)
                                
                !---Diffusion TERM
                viscous_w1(Sc:Ec) = 1.d0*(mu2(Sc:Ec)*dwdx2(Sc:Ec) - mu1(Sc:Ec)*dwdx1(Sc:Ec))/dx1(Sc:Ec)
                viscous_w2(Sc:Ec) = 1.d0*(mu4(Sc:Ec)*dwdy4(Sc:Ec) - mu3(Sc:Ec)*dwdy3(Sc:Ec))/dx2(j )
                viscous_w3(Sc:Ec) = 1.d0*(Mu(Sc:Ec,j ,k )*dwdz6(Sc:Ec) - Mu(Sc:Ec,j ,km)*dwdz5(Sc:Ec))/dmx3(k)
                viscous_w31(Sc:Ec) = 1.d0*(mu2(Sc:Ec)*dudz2(Sc:Ec) - mu1(Sc:Ec)*dudz1(Sc:Ec))/dx1(Sc:Ec)
                viscous_w32(Sc:Ec) = 1.d0*(mu4(Sc:Ec)*dvdz4(Sc:Ec) - mu3(Sc:Ec)*dvdz3(Sc:Ec))/dx2(j)

                RHS(Sc:Ec,j,k) =  0.5d0*Cmu*invRhoc(Sc:Ec)*(viscous_w1(Sc:Ec) + viscous_w2(Sc:Ec) + 2.*viscous_w3(Sc:Ec) + viscous_w31(Sc:Ec) + viscous_w32(Sc:Ec))

                !---PRESSURE TERM
                RHS(Sc:Ec,j,k) = RHS(Sc:Ec,j,k)   &
                               - Cmp*invRhoc(Sc:Ec)*(P(Sc:Ec,j,k) - P(Sc:Ec,j,km))/dmx3(k)

                !---wbc
                zbc_down(Sc:Ec) = 0.25d0*v3(Sc:Ec)/dmx2(j)*WBCbt_sub(Sc:Ec,k)    &
                                - 0.25d0*dwdy3(Sc:Ec)*(dx3(km)*VBCbt_sub(Sc:Ec,k) + dx3(k)*VBCbt_sub(Sc:Ec,km))/dmx3(k)*0.5d0     

                zbc_up(Sc:Ec) = - 0.25d0*v4(Sc:Ec)/dmx2(jp)*WBCup_sub(Sc:Ec,k)   &
                                - 0.25d0*dwdy4(Sc:Ec)*(dx3(km)*VBCup_sub(Sc:Ec,k) + dx3(k)*VBCup_sub(Sc:Ec,km))/dmx3(k)*0.5d0

                zbc_down(Sc:Ec) = zbc_down(Sc:Ec) + 0.5d0*Cmu*invRhoc(Sc:Ec)/dx2(j)*mu3(Sc:Ec)/dmx2(jp)*WBCbt_sub(Sc:Ec,k)
                zbc_up(Sc:Ec)   = zbc_up(Sc:Ec)   + 0.5d0*Cmu*invRhoc(Sc:Ec)/dx2(j)*mu4(Sc:Ec)/dmx2(j )*WBCup_sub(Sc:Ec,k)

                zbc_down(Sc:Ec) = zbc_down(Sc:Ec) - 0.5d0*Cmu*invRhoc(Sc:Ec)/dx2(j)*mu3(Sc:Ec)/dmx3(k)*(VBCbt_sub(Sc:Ec,k) - VBCbt_sub(Sc:Ec,km))
                zbc_up(Sc:Ec)   = zbc_up(Sc:Ec)   + 0.5d0*Cmu*invRhoc(Sc:Ec)/dx2(j)*mu4(Sc:Ec)/dmx3(k)*(VBCup_sub(Sc:Ec,k) - VBCup_sub(Sc:Ec,km))

                RHS(Sc:Ec,j,k) = RHS(Sc:Ec,j,k)                     &
                               + dble(1.d0 - jwm)*zbc_down(Sc:Ec)   &
                               + dble(1.d0 - jwp)*zbc_up(Sc:Ec)

                !---M33Wn
                mACI(Sc:Ec) = 0.5d0*Cmu*invRhoc(Sc:Ec)/dx1(Sc:Ec)*(mu2(Sc:Ec)/dmx1(Sp:Ep) + mu1(Sc:Ec)/dmx1(Sc:Ec))     &
                            + 0.25d0*(-u2(Sc:Ec)/dmx1(Sp:Ep) + u1(Sc:Ec)/dmx1(Sc:Ec))
                mAPI(Sc:Ec) = -0.5d0*Cmu*invRhoc(Sc:Ec)/dx1(Sc:Ec)*mu2(Sc:Ec)/dmx1(Sp:Ep)   &
                            + 0.25d0*( u2(Sc:Ec)/dmx1(Sp:Ep))
                mAMI(Sc:Ec) = -0.5d0*Cmu*invRhoc(Sc:Ec)/dx1(Sc:Ec)*mu1(Sc:Ec)/dmx1(Sc:Ec)   &
                            + 0.25d0*(-u1(Sc:Ec)/dmx1(Sc:Ec))

                mACJ(Sc:Ec) = 0.5d0*Cmu*invRhoc(Sc:Ec)/dx2(j)*(mu4(Sc:Ec)/dmx2(jp) + mu3(Sc:Ec)/dmx2(j))    &
                            + 0.25d0*(-v4(Sc:Ec)/dmx2(jp) + v3(Sc:Ec)/dmx2(j ))
                mAPJ(Sc:Ec) = -0.5d0*Cmu*invRhoc(Sc:Ec)/dx2(j)*mu4(Sc:Ec)/dmx2(jp)  &
                            + 0.25d0*( v4(Sc:Ec)/dmx2(jp))
                mAPJ(Sc:Ec) = mAPJ(Sc:Ec)*dble(jwp)
                mAMJ(Sc:Ec) = -0.5d0*Cmu*invRhoc(Sc:Ec)/dx2(j)*mu3(Sc:Ec)/dmx2(j)   &
                            + 0.25d0*(-v3(Sc:Ec)/dmx2(j))
                mAMJ(Sc:Ec) = mAMJ(Sc:Ec)*dble(jwm)
                
                mACK(Sc:Ec) =  0.5d0*Cmu*invRhoc(Sc:Ec)/dmx3(k)*(Mu(Sc:Ec,j ,k )/dx3(k) + Mu(Sc:Ec,j ,km)/dx3(km))*2.    &
                            + 0.25d0*(-w6(Sc:Ec)/dx3(k) + w5(Sc:Ec)/dx3(km) + 0.5d0*dwdz6(Sc:Ec) + 0.5d0*dwdz5(Sc:Ec))
                mAPK(Sc:Ec) = -0.5d0*Cmu*invRhoc(Sc:Ec)/dmx3(k)*Mu(Sc:Ec,j ,k )/dx3(k)*2.             &
                            + 0.25d0*( w6(Sc:Ec)/dx3(k) + 0.5d0*dwdz6(Sc:Ec))
                mAMK(Sc:Ec) = -0.5d0*Cmu*invRhoc(Sc:Ec)/dmx3(k)*Mu(Sc:Ec,j ,km)/dx3(km)*2.              &
                            + 0.25d0*(-w5(Sc:Ec)/dx3(km) + 0.5d0*dwdz5(Sc:Ec))

                RHS(Sc:Ec,j,k) = RHS(Sc:Ec,j,k)                     &
                               - (mAPI(Sc:Ec)*W(Sp:Ep,j, k ) + mACI(Sc:Ec)*W(Sc:Ec,j,k) + mAMI(Sc:Ec)*W(Sm:Em,j, k )    &
                                 +mAPJ(Sc:Ec)*W(Sc:Ec,jp,k ) + mACJ(Sc:Ec)*W(Sc:Ec,j,k) + mAMJ(Sc:Ec)*W(Sc:Ec,jm,k )    &
                                 +mAPK(Sc:Ec)*W(Sc:Ec,j, kp) + mACK(Sc:Ec)*W(Sc:Ec,j,k) + mAMK(Sc:Ec)*W(Sc:Ec,j, km)    )
                        
                !---M31Un    
                RHS(Sc:Ec,j,k) = RHS(Sc:Ec,j,k)                     &
                               - (0.25d0*(u2(Sc:Ec)*dwdx2(Sc:Ec) + u1(Sc:Ec)*dwdx1(Sc:Ec))                                  & 
                                 -0.5d0*Cmu*invRhoc(Sc:Ec)/dx1(Sc:Ec)*(mu2(Sc:Ec)*dudz2(Sc:Ec) - mu1(Sc:Ec)*dudz1(Sc:Ec))   )

                !---M32Vn
                RHS(Sc:Ec,j,k) = RHS(Sc:Ec,j,k)                     &
                               -(0.25d0*(v3(Sc:Ec)*dble(jwm)*dwdy3(Sc:Ec) + v4(Sc:Ec)*dble(jwp)*dwdy4(Sc:Ec))                                                   &
                                -0.5d0*Cmu*invRhoc(Sc:Ec)/dx2(j)*(mu4(Sc:Ec)*dvdz4(Sc:Ec)*dble(jwp) - mu3(Sc:Ec)*dvdz3(Sc:Ec)*dble(jwm))    )

                !---M31 ddU
                RHS(Sc:Ec,j,k) = RHS(Sc:Ec,j,k)                     &
                               -(0.25d0*(ddu2(Sc:Ec)*dwdx2(Sc:Ec) + ddu1(Sc:Ec)*dwdx1(Sc:Ec))                               &
                                -0.5d0*Cmu*invRhoc(Sc:Ec)/dx1(Sc:Ec)*(mu2(Sc:Ec)*ddudz2(Sc:Ec) - mu1(Sc:Ec)*ddudz1(Sc:Ec))  )

                !---M32 ddV
                RHS(Sc:Ec,j,k) = RHS(Sc:Ec,j,k)                     &
                               -(0.25d0*(ddv3(Sc:Ec)*dble(jwm)*dwdy3(Sc:Ec) + ddv4(Sc:Ec)*dble(jwp)*dwdy4(Sc:Ec))                                               &
                                -0.5d0*Cmu*invRhoc(Sc:Ec)/dx2(j)*(mu4(Sc:Ec)*ddvdz4(Sc:Ec)*dble(jwp) - mu3(Sc:Ec)*ddvdz3(Sc:Ec)*dble(jwm))  )

                ! M33WVn: Z-direction
                ac(Sc:Ec,j,k) =mACK(Sc:Ec)*dt + 1.d0
                ap(Sc:Ec,j,k) =mAPK(Sc:Ec)*dt
                am(Sc:Ec,j,k) =mAMK(Sc:Ec)*dt
                RHS(Sc:Ec,j,k)=RHS(Sc:Ec,j,k)*dt

            end do  
        end do
        

        ! TDMA-i*j-K
        
        call PaScaL_TDMA_plan_many_create(ptdma_plan, (n1sb-1)*(n2sb-1), comm_1d_x3%myrank, comm_1d_x3%nprocs, comm_1d_x3%mpi_comm)
        call PaScaL_TDMA_many_solve_cycle(ptdma_plan, am, ac, ap, RHS,(n1sb-1)*(n2sb-1),(n3sb-1))
        call PaScaL_TDMA_plan_many_destroy(ptdma_plan,comm_1d_x3%nprocs)
        
        deallocate(am,ac,ap)

        ! Trans i-j-K->j-k-i
        
        allocate(ddW(1:n2sb-1,1:n1sb-1,1:n3sb-1))
        do k = k_indexS, n3sb-1
            ddW(:,:,k)=TRANSPOSE(RHS(:,:,k))
        enddo
        deallocate(RHS) 
        allocate(RHS(1:n2sb-1,1:n3sb-1,1:n1sb-1))
        do k = k_indexS, n3sb-1
            do i = i_indexC, n1sb-1
                RHS(1:n2sb-1,k,i)=ddW(1:n2sb-1,i,k)
            enddo
        enddo
        deallocate(ddW)
        

        ! j-k-i
        
        allocate(am(1:n2sb-1,1:n3sb-1,1:n1sb-1),ac(1:n2sb-1,1:n3sb-1,1:n1sb-1),ap(1:n2sb-1,1:n3sb-1,1:n1sb-1))
        do k = k_indexS, n3sb-1
            kp = k+1
            km = k-1
            do j = j_indexC, n2sb-1
                jp = j + 1
                jm = j - 1
                jwp = jC_BC(jp)
                jwm = jC_BC(jm)
                    
                u1(Sc:Ec) = (dx3(km)*U(Sc:Ec,j,k) + dx3(k)*U(Sc:Ec, j,km))/dmx3(k)*0.5d0
                u2(Sc:Ec) = (dx3(km)*U(Sp:Ep,j,k) + dx3(k)*U(Sp:Ep,j,km))/dmx3(k)*0.5d0                
            
                invRhoc(Sc:Ec) = 0.5d0*(dx3(km)*invRho(Sc:Ec,j ,k ) + dx3(k )*invRho(Sc:Ec,j ,km) )/dmx3(k)
                                
                ! M33Wn: X-direction
                mACI(Sc:Ec) = 0.5d0*Cmu*invRhoc(Sc:Ec)/dx1(Sc:Ec)*(mu2(Sc:Ec)/dmx1(Sp:Ep) + mu1(Sc:Ec)/dmx1(Sc:Ec))     &
                            + 0.25d0*(-u2(Sc:Ec)/dmx1(Sp:Ep) + u1(Sc:Ec)/dmx1(Sc:Ec))
                mAPI(Sc:Ec) = -0.5d0*Cmu*invRhoc(Sc:Ec)/dx1(Sc:Ec)*mu2(Sc:Ec)/dmx1(Sp:Ep)   &
                            + 0.25d0*( u2(Sc:Ec)/dmx1(Sp:Ep))
                mAMI(Sc:Ec) = -0.5d0*Cmu*invRhoc(Sc:Ec)/dx1(Sc:Ec)*mu1(Sc:Ec)/dmx1(Sc:Ec)   &
                            + 0.25d0*(-u1(Sc:Ec)/dmx1(Sc:Ec))
                            
                do i = i_indexC, n1sb-1
                    ac(j,k,i)=mACI(i)*dt + 1.d0
                    ap(j,k,i)=mAPI(i)*dt
                    am(j,k,i)=mAMI(i)*dt
                enddo
            enddo
        enddo
        
        
        ! TDMA-j*K-i
        
        call PaScaL_TDMA_plan_many_create(ptdma_plan, (n2sb-1)*(n3sb-1), comm_1d_x1%myrank, comm_1d_x1%nprocs, comm_1d_x1%mpi_comm)
        call PaScaL_TDMA_many_solve_cycle(ptdma_plan, am, ac, ap, RHS,(n2sb-1)*(n3sb-1),(n1sb-1))
        call PaScaL_TDMA_plan_many_destroy(ptdma_plan,comm_1d_x1%nprocs)
        
        deallocate(am,ac,ap)
        
        ! Trans j-k-i->k-i-j 
        
        allocate(ddW(1:n3sb-1,1:n2sb-1,1:n1sb-1))
        do i = i_indexC, n1sb-1
            ddW(:,:,i)=TRANSPOSE(RHS(:,:,i))
        enddo
        deallocate(RHS) 
        allocate(RHS(1:n3sb-1,1:n1sb-1,1:n2sb-1))
        do j = j_indexC, n2sb-1
            do i = i_indexC, n1sb-1
                RHS(1:n3sb-1,i,j)=ddW(1:n3sb-1,j,i)
            enddo
        enddo
        deallocate(ddW)
        

        ! k-i-j
        
        allocate(am(1:n3sb-1,1:n1sb-1,1:n2sb-1),ac(1:n3sb-1,1:n1sb-1,1:n2sb-1),ap(1:n3sb-1,1:n1sb-1,1:n2sb-1))
        do k = k_indexS, n3sb-1
            do j = j_indexC, n2sb-1
                jm=j-1
                jp=j+1
                jwp = jC_BC(jp)
                jwm = jC_BC(jm)

                v3(Sc:Ec) = (dx3(km)*V(Sc:Ec,j, k) + dx3(k)*V(Sc:Ec,j, km))/dmx3(k)*0.5d0
                v4(Sc:Ec) = (dx3(km)*V(Sc:Ec,jp,k) + dx3(k)*V(Sc:Ec,jp,km))/dmx3(k)*0.5d0                
            
                invRhoc(Sc:Ec) = 0.5d0*(dx3(km)*invRho(Sc:Ec,j ,k ) + dx3(k )*invRho(Sc:Ec,j ,km) )/dmx3(k)

                ! M33Wn: Y-direction
                    mACJ(Sc:Ec) = 0.5d0*Cmu*invRhoc(Sc:Ec)/dx2(j)*(mu4(Sc:Ec)/dmx2(jp) + mu3(Sc:Ec)/dmx2(j))    &
                                + 0.25d0*(-v4(Sc:Ec)/dmx2(jp) + v3(Sc:Ec)/dmx2(j ))
                    mAPJ(Sc:Ec) = -0.5d0*Cmu*invRhoc(Sc:Ec)/dx2(j)*mu4(Sc:Ec)/dmx2(jp)  &
                                + 0.25d0*( v4(Sc:Ec)/dmx2(jp))
                    mAPJ(Sc:Ec) = mAPJ(Sc:Ec)*dble(jwp)
                    mAMJ(Sc:Ec) = -0.5d0*Cmu*invRhoc(Sc:Ec)/dx2(j)*mu3(Sc:Ec)/dmx2(j)   &
                                + 0.25d0*(-v3(Sc:Ec)/dmx2(j))
                    mAMJ(Sc:Ec) = mAMJ(Sc:Ec)*dble(jwm)

                do i = i_indexC, n1sb-1
                    ac(k,i,j)=mACJ(i)*dt + 1.d0
                    ap(k,i,j)=mAPJ(i)*dt
                    am(k,i,j)=mAMJ(i)*dt
                enddo
            enddo
        enddo
        

        ! TDMA-k*i-j
        
        call PaScaL_TDMA_plan_many_create(ptdma_plan, (n3sb-1)*(n1sb-1), comm_1d_x2%myrank, comm_1d_x2%nprocs, comm_1d_x2%mpi_comm)
        call PaScaL_TDMA_many_solve(ptdma_plan, am, ac, ap, RHS,(n3sb-1)*(n1sb-1),(n2sb-1))
        call PaScaL_TDMA_plan_many_destroy(ptdma_plan,comm_1d_x2%nprocs)
        
        deallocate(am,ac,ap)

        ! Trans k-i-j->i-j-k
        
        allocate(ddW(1:n1sb-1,1:n3sb-1,1:n2sb-1))
        do j = j_indexC, n2sb-1
            ddW(:,:,j)=TRANSPOSE(RHS(:,:,j))
        enddo
        deallocate(RHS)
        allocate(dW(0:n1sb,0:n2sb,0:n3sb))
        do k = k_indexS, n3sb-1
            do j = j_indexC, n2sb-1
                dW(1:n1sb-1,j,k)=ddW(1:n1sb-1,k,j)
            enddo
        enddo
        deallocate(ddW)
        

        deallocate(u1,u2,v3,v4,w6,w5)
        deallocate(dudz1,dudz2,dvdz3,dvdz4,dwdx1,dwdx2,dwdy3,dwdy4,dwdz5,dwdz6)
        deallocate(ddu1,ddu2,ddudz1,ddudz2)
        deallocate(ddv3,ddv4,ddvdz3,ddvdz4)
        deallocate(mua,muc,mub,mu1,mu2,mu3,mu4,invRhoc)
        deallocate(viscous_w1,viscous_w2,viscous_w3,viscous_w31,viscous_w32)
        deallocate(zbc_up,zbc_down)
        deallocate(mACI,mAPI,mAMI,mACJ,mAPJ,mAMJ,mACK,mAPK,mAMK)

    end subroutine mpi_momentum_solvedW
    
    subroutine mpi_momentum_blockLdV(T,dx1,dx2,dx3,dmx1,dmx2,dmx3)
        implicit none
        double precision, dimension(0:n1sb,0:n2sb,0:n3sb) ::  T
        double precision :: dx1(0:n1sb),dx2(0:n2sb),dx3(0:n3sb)
        double precision :: dmx1(0:n1sb),dmx2(0:n2sb),dmx3(0:n3sb)
            
        integer :: i,j,k
        integer :: im,jm,km
        integer :: ip,jp,kp
        integer :: jvp,jvm

        integer :: Sc,Ec
        integer :: Sm,Em
        integer :: Sp,Ep
        double precision, allocatable, dimension(:) :: dwm5,dwm6
        double precision, allocatable, dimension(:) :: dvdz5,dvdz6,ddwdy5,ddwdy6
        double precision, allocatable, dimension(:) :: mua,muc,mub,mu5,mu6
        double precision, allocatable, dimension(:) :: invRhoc

        Sc=i_indexC  ;Ec=n1sb-1
        Sm=i_indexC-1;Em=n1sb-1-1
        Sp=i_indexC+1;Ep=n1sb-1+1
        
        allocate(dwm5(Sc:Ec),dwm6(Sc:Ec))
        allocate(dvdz5(Sc:Ec),dvdz6(Sc:Ec),ddwdy5(Sc:Ec),ddwdy6(Sc:Ec))
        allocate(mua(Sc:Ec),muc(Sc:Ec),mub(Sc:Ec),mu5(Sc:Ec),mu6(Sc:Ec))
        allocate(invRhoc(Sc:Ec))

        ! i-j-k
        do k = 1, n3sb-1
            kp = k+1
            km = k-1
            do j =  j_indexS, n2sb-1
            jp = j + 1
            jm = j - 1
            jvp = jS_BC(jm)
            jvm = jS_BC(jp)
        
            dwm5(Sc:Ec) = (dx2(jm)*dW(Sc:Ec,j,k ) + dx2(j)*dW(Sc:Ec,jm,k ))/dmx2(j)*0.5d0
            dwm6(Sc:Ec) = (dx2(jm)*dW(Sc:Ec,j,kp) + dx2(j)*dW(Sc:Ec,jm,kp))/dmx2(j)*0.5d0

            ddwdy5(Sc:Ec) = (dW(Sc:Ec,j,k ) - dW(Sc:Ec,jm,k ))/dmx2(j)
            ddwdy6(Sc:Ec) = (dW(Sc:Ec,j,kp) - dW(Sc:Ec,jm,kp))/dmx2(j)  
            
            dvdz5(Sc:Ec) = (V(Sc:Ec,j,k ) - V(Sc:Ec,j,km))/dmx3(k )
            dvdz6(Sc:Ec) = (V(Sc:Ec,j,kp) - V(Sc:Ec,j,k ))/dmx3(kp)

            invRhoc(Sc:Ec) = 0.5d0*(dx2(jm)*invRho(Sc:Ec,j ,k ) + dx2(j )*invRho(Sc:Ec,jm,k ))/dmx2(j)
            
            mua(Sc:Ec) = 0.5d0*(dx2(jm)*Mu(Sc:Ec,j ,km) + dx2(j)*Mu(Sc:Ec,jm,km))/dmx2(j)
            muc(Sc:Ec) = 0.5d0*(dx2(jm)*Mu(Sc:Ec,j ,k ) + dx2(j)*Mu(Sc:Ec,jm,k ))/dmx2(j)
            mub(Sc:Ec) = 0.5d0*(dx2(jm)*Mu(Sc:Ec,j ,kp) + dx2(j)*Mu(Sc:Ec,jm,kp))/dmx2(j)

            mu5(Sc:Ec) = 0.5d0*(dx3(km)*muc(Sc:Ec) + dx3(k)*mua(Sc:Ec))/dmx3(k )
            mu6(Sc:Ec) = 0.5d0*(dx3(kp)*muc(Sc:Ec) + dx3(k)*mub(Sc:Ec))/dmx3(kp)

            invRhoc(Sc:Ec) = 0.5d0*(dx2(jm)*invRho(Sc:Ec,j ,k ) + dx2(j )*invRho(Sc:Ec,jm,k ))/dmx2(j)
            !>  dV(i,j,k) = dV(i,j,k) - dt*M23dWm        
            dV(Sc:Ec,j,k) = dV(Sc:Ec,j,k)     &
                          - dt*(0.25d0*(dwm5(Sc:Ec)*dvdz5(Sc:Ec) + dwm6(Sc:Ec)*dvdz6(Sc:Ec))       &
                               -0.5d0*Cmu*invRhoc(Sc:Ec)/dx3(k)*(mu6(Sc:Ec)*ddwdy6(Sc:Ec) - mu5(Sc:Ec)*ddwdy5(Sc:Ec))   )

            end do
        end do

        deallocate(dwm5,dwm6)
        deallocate(dvdz5,dvdz6,ddwdy5,ddwdy6)
        deallocate(mua,muc,mub,mu5,mu6)
        deallocate(invRhoc)
    end subroutine mpi_momentum_blockLdV

    subroutine mpi_momentum_blockLdU(T,dx1,dx2,dx3,dmx1,dmx2,dmx3)
        implicit none
        double precision, dimension(0:n1sb,0:n2sb,0:n3sb) ::  T
        double precision :: dx1(0:n1sb),dx2(0:n2sb),dx3(0:n3sb)
        double precision :: dmx1(0:n1sb),dmx2(0:n2sb),dmx3(0:n3sb)
            
        integer :: i,j,k
        integer :: im,jm,km
        integer :: ip,jp,kp
        integer :: jup,jum

        integer :: Sc,Ec
        integer :: Sm,Em
        integer :: Sp,Ep
        double precision, allocatable, dimension(:) :: dvm3,dvm4
        double precision, allocatable, dimension(:) :: dudy3,dudy4,dvmdx3,dvmdx4
        double precision, allocatable, dimension(:) :: mua,mub,muc,mu3,mu4
        double precision, allocatable, dimension(:) :: invRhoc,viscous_u1

        double precision, allocatable, dimension(:) :: dwm5,dwm6
        double precision, allocatable, dimension(:) :: dudz5,dudz6,dwmdx5,dwmdx6
        double precision, allocatable, dimension(:) :: mu5,mu6

        Sc=i_indexS  ;Ec=n1sb-1
        Sm=i_indexS-1;Em=n1sb-1-1
        Sp=i_indexS+1;Ep=n1sb-1+1
        
        allocate(dvm3(Sc:Ec),dvm4(Sc:Ec))
        allocate(dudy3(Sc:Ec),dudy4(Sc:Ec),dvmdx3(Sc:Ec),dvmdx4(Sc:Ec))
        allocate(mua(Sc:Ec),mub(Sc:Ec),muc(Sc:Ec),mu3(Sc:Ec),mu4(Sc:Ec))
        allocate(dwm5(Sc:Ec),dwm6(Sc:Ec))
        allocate(dudz5(Sc:Ec),dudz6(Sc:Ec),dwmdx5(Sc:Ec),dwmdx6(Sc:Ec))
        allocate(mu5(Sc:Ec),mu6(Sc:Ec))
        allocate(invRhoc(Sc:Ec),viscous_u1(Sc:Ec))

        ! i-j-k:RHS-ddU
        ! i-j-k
        do k = k_indexC, n3sb-1
            km=k-1
            kp=k+1
            do j = j_indexC, n2sb-1
                jm = j-1
                jp = j+1
                jum = jC_BC(jm)
                jup = jC_BC(jp)
               
                invRhoc(Sc:Ec) = 0.5d0*(dx1(Sm:Em)*invRho(Sc:Ec,j ,k ) + dx1(Sc:Ec)*invRho(Sm:Em,j ,k ) )/dmx1(Sc:Ec)

                dvm3(Sc:Ec) = (dx1(Sc:Ec)*dV(Sm:Em,j, k) + dx1(Sm:Em)*dV(Sc:Ec,j, k))/dmx1(Sc:Ec)*0.5d0*dble(jum)
                dvm4(Sc:Ec) = (dx1(Sc:Ec)*dV(Sm:Em,jp,k) + dx1(Sm:Em)*dV(Sc:Ec,jp,k))/dmx1(Sc:Ec)*0.5d0*dble(jup)
            
                dudy3(Sc:Ec) = (U(Sc:Ec,j, k) - U(Sc:Ec,jm,k))/dmx2(j )
                dudy4(Sc:Ec) = (U(Sc:Ec,jp,k) - U(Sc:Ec,j, k))/dmx2(jp)
            
                dvmdx3 = (dV(Sc:Ec,j ,k) - dV(Sm:Em,j ,k))/dmx1(Sc:Ec)*dble(jum)
                dvmdx4 = (dV(Sc:Ec,jp,k) - dV(Sm:Em,jp,k))/dmx1(Sc:Ec)*dble(jup)
            
                mua(Sc:Ec) = 0.5d0*(dx1(Sm:Em)*Mu(Sc:Ec,jm,k ) + dx1(Sc:Ec)*Mu(Sm:Em,jm,k ))/dmx1(Sc:Ec)
                muc(Sc:Ec) = 0.5d0*(dx1(Sm:Em)*Mu(Sc:Ec,j ,k ) + dx1(Sc:Ec)*Mu(Sm:Em,j ,k ))/dmx1(Sc:Ec)
                mub(Sc:Ec) = 0.5d0*(dx1(Sm:Em)*Mu(Sc:Ec,jp,k ) + dx1(Sc:Ec)*Mu(Sm:Em,jp,k ))/dmx1(Sc:Ec)
                mu3(Sc:Ec) = 0.5d0*(dx2(jm)*muc(Sc:Ec) + dx2(j)*mua(Sc:Ec))/dmx2(j )
                mu4(Sc:Ec) = 0.5d0*(dx2(jp)*muc(Sc:Ec) + dx2(j)*mub(Sc:Ec))/dmx2(jp)
            
                dwm5(Sc:Ec) = (dx1(Sc:Ec)*dW(Sm:Em,j,k ) + dx1(Sm:Em)*dW(Sc:Ec,j,k ))/dmx1(Sc:Ec)*0.5d0 
                dwm6(Sc:Ec) = (dx1(Sc:Ec)*dW(Sm:Em,j,kp) + dx1(Sm:Em)*dW(Sc:Ec,j,kp))/dmx1(Sc:Ec)*0.5d0 
            
                dudz5(Sc:Ec) = (U(Sc:Ec,j,k ) - U(Sc:Ec,j,km))/dmx3(k )
                dudz6(Sc:Ec) = (U(Sc:Ec,j,kp) - U(Sc:Ec,j,k ))/dmx3(kp)
            
                dwmdx5(Sc:Ec) = (dW(Sc:Ec,j,k ) - dW(Sm:Em,j,k ))/dmx1(Sc:Ec)
                dwmdx6(Sc:Ec) = (dW(Sc:Ec,j,kp) - dW(Sm:Em,j,kp))/dmx1(Sc:Ec)
            
                mua(Sc:Ec) = 0.5d0*(dx1(Sm:Em)*Mu(Sc:Ec,j ,km) + dx1(Sc:Ec)*Mu(Sm:Em,j ,km))/dmx1(Sc:Ec)
                muc(Sc:Ec) = 0.5d0*(dx1(Sm:Em)*Mu(Sc:Ec,j ,k ) + dx1(Sc:Ec)*Mu(Sm:Em,j ,k ))/dmx1(Sc:Ec)
                mub(Sc:Ec) = 0.5d0*(dx1(Sm:Em)*Mu(Sc:Ec,j ,kp) + dx1(Sc:Ec)*Mu(Sm:Em,j ,kp))/dmx1(Sc:Ec)            
                mu5(Sc:Ec) = 0.5d0*(dx3(km)*muc(Sc:Ec) + dx3(k)*mua(Sc:Ec))/dmx3(k )
                mu6(Sc:Ec) = 0.5d0*(dx3(kp)*muc(Sc:Ec) + dx3(k)*mub(Sc:Ec))/dmx3(kp)
                
                !>  dU(i,j,k) = dU(i,j,k) - dt*M12dVm - dt*M13dWm
                dU(Sc:Ec,j,k) = dU(Sc:Ec,j,k)     &
                            -dt*( 0.25d0*(dvm4(Sc:Ec)*dudy4(Sc:Ec) + dvm3(Sc:Ec)*dudy3(Sc:Ec))                      &
                                -0.5d0*Cmu*invRhoc(Sc:Ec)/dx2(j)*(mu4(Sc:Ec)*dvmdx4 - mu3(Sc:Ec)*dvmdx3)    )     &
                            -dt*( 0.25d0*(dwm6(Sc:Ec)*dudz6(Sc:Ec) + dwm5(Sc:Ec)*dudz5(Sc:Ec))                      &
                                -0.5d0*Cmu*invRhoc(Sc:Ec)/dx3(k)*(mu6(Sc:Ec)*dwmdx6 - mu5(Sc:Ec)*dwmdx5)   )

            enddo
        enddo



        deallocate(dvm3,dvm4)
        deallocate(dudy3,dudy4,dvmdx3,dvmdx4)
        deallocate(mua,mub,muc,mu3,mu4)
        deallocate(dwm5,dwm6)
        deallocate(dudz5,dudz6,dwmdx5,dwmdx6)
        deallocate(mu5,mu6)
        deallocate(invRhoc,viscous_u1)
    end subroutine mpi_momentum_blockLdU

end module mpi_momentum
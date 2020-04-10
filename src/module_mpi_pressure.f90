module mpi_pressure
    ! use debug
    use global
    implicit none
    integer :: nn1sb,nn2sb,nn3sb
    integer :: n1msub,n3msub,n3mssub,wv1,wv1sub,wv1ssub,wv1ssubA,wv1ssubB
    integer :: nn1mI,wv1I,nn3mIssb
    integer :: wv1sb,nn3msb
    integer :: wv1Kssb,nn3mK,wv3K
    integer, allocatable, dimension(:) :: ddtype_C2I_dble,ddtype_I2C_dble
    integer, allocatable, dimension(:) :: ddtype_C2I_cplx,ddtype_I2C_cplx
    integer, allocatable, dimension(:) :: ddtype_C2K_cplx,ddtype_K2C_cplx

    integer, allocatable, dimension(:) :: countsendI, countdistI
    integer, allocatable, dimension(:) :: countsendK, countdistK

    double precision, allocatable, dimension(:,:,:) :: dPhat,dP0,PRHS,dP
contains    
    subroutine mpi_pressure_allocation(n1in,n2in,n3in, comm_1d_x1, comm_1d_x2, comm_1d_x3)
        use mpi_topology, only : cart_comm_1d
        implicit none
        integer :: n1in,n2in,n3in
        type(cart_comm_1d), intent(in)  :: comm_1d_x1, comm_1d_x2, comm_1d_x3

        nn1sb=n1in
        nn2sb=n2in
        nn3sb=n3in
        
        allocate(ddtype_C2I_dble(0:comm_1d_x1%nprocs-1),ddtype_I2C_dble(0:comm_1d_x1%nprocs-1))
        allocate(ddtype_C2I_cplx(0:comm_1d_x1%nprocs-1),ddtype_I2C_cplx(0:comm_1d_x1%nprocs-1))
        allocate(ddtype_C2K_cplx(0:comm_1d_x3%nprocs-1),ddtype_K2C_cplx(0:comm_1d_x3%nprocs-1))
        allocate(countsendI(0:comm_1d_x1%nprocs-1), countdistI(0:comm_1d_x1%nprocs-1))
        allocate(countsendK(0:comm_1d_x3%nprocs-1), countdistK(0:comm_1d_x3%nprocs-1))

        allocate(dPhat(0:nn1sb,0:nn2sb,0:nn3sb),dP0(0:nn1sb,0:nn2sb,0:nn3sb))

        countsendI(:)=1
        countdistI(:)=0
        countsendK(:)=1
        countdistK(:)=0

        dPhat(:,:,:) = 0.d0
        dP0(:,:,:) = 0.d0

    end subroutine mpi_pressure_allocation

    subroutine mpi_pressure_clean()
        implicit none
        
        deallocate(ddtype_C2I_dble,ddtype_I2C_dble)
        deallocate(ddtype_C2I_cplx,ddtype_I2C_cplx)
        deallocate(ddtype_C2K_cplx,ddtype_K2C_cplx)
        deallocate(countsendI, countdistI)
        deallocate(countsendK, countdistK)
        deallocate(dPhat,dP0)
        
    end subroutine mpi_pressure_clean

    subroutine mpi_pressure_initial(x1,x2,x3)
        implicit none
        double precision:: x1(0:nn1sb),x2(0:nn2sb),x3(0:nn3sb)
        integer :: i,j,k
        
        dPhat(0:nn1sb,0:nn2sb,0:nn3sb)=0.d0
        dP0(0:nn1sb,0:nn2sb,0:nn3sb) = 0.d0
    end subroutine mpi_pressure_initial

    subroutine mpi_pressure_RHS(invRho,U,V,W,VBCup_sub,VBCbt_sub,dx1,dx2,dx3,dmx1,dmx2,dmx3)
        implicit none
        double precision, dimension(0:nn1sb,0:nn2sb,0:nn3sb) :: invRho,U,V,W
        double precision, dimension(0:nn1sb,0:nn3sb) :: VBCup_sub,VBCbt_sub
        double precision :: dx1(0:nn1sb),dx2(0:nn2sb),dx3(0:nn3sb)
        double precision :: dmx1(0:nn1sb),dmx2(0:nn2sb),dmx3(0:nn3sb)
        integer :: i,j,k
        integer :: kp,km  
        integer :: jp,jm,jvm,jvp
        integer :: Sc,Ec
        integer :: Sm,Em
        integer :: Sp,Ep

        double precision, allocatable, dimension(:) :: DivUm,cbc
        double precision, allocatable, dimension(:) :: ddpdx1,ddpdx2,ddpdy3,ddpdy4,ddpdz5,ddpdz6
        double precision, allocatable, dimension(:) :: invRho1,invRho2,invRho3,invRho4,invRho5,invRho6
        double precision, allocatable, dimension(:) :: ExtraTerm

        Sc=i_indexS  ;Ec=nn1sb-1
        Sm=i_indexS-1;Em=nn1sb-1-1
        Sp=i_indexS+1;Ep=nn1sb-1+1

        allocate(PRHS(1:nn1sb-1,1:nn2sb-1,1:nn3sb-1))
        
        allocate(DivUm(Sc:Ec),cbc(Sc:Ec))
        allocate(ddpdx1(Sc:Ec),ddpdx2(Sc:Ec),ddpdy3(Sc:Ec),ddpdy4(Sc:Ec),ddpdz5(Sc:Ec),ddpdz6(Sc:Ec))
        allocate(invRho1(Sc:Ec),invRho2(Sc:Ec),invRho3(Sc:Ec),invRho4(Sc:Ec),invRho5(Sc:Ec),invRho6(Sc:Ec))
        allocate(ExtraTerm(Sc:Ec))

        do k = 1, nn3sb-1
            kp = k+1
            km = k-1
            do j = 1, nn2sb-1
                jp = j + 1
                jm = j - 1
                jvm = jC_BC(jm)
                jvp = jC_BC(jp)
            
                DivUm(Sc:Ec) = (            U(Sp:Ep,j,k) -            U(Sc:Ec,j,k)  )/dx1(Sc:Ec)    &
                             + ( dble(jvp) *V(Sc:Ec,jp,k) - dble(jvm)*V(Sc:Ec,j,k)  )/dx2(j)        &
                             + (            W(Sc:Ec,j,kp) -           W(Sc:Ec,j,k)  )/dx3(k)

                cbc(Sc:Ec) = dble(1. - jvm)*VBCbt_sub(Sc:Ec,k)/dx2(j)    &
                           - dble(1. - jvp)*VBCup_sub(Sc:Ec,k)/dx2(j)

                ddpdx1(Sc:Ec) = (dPhat(Sc:Ec, j, k) - dPhat(Sm:Em,j, k))/dmx1(Sc:Ec)
                ddpdx2(Sc:Ec) = (dPhat(Sp:Ep,j, k) - dPhat(Sc:Ec,j, k))/dmx1(Sp:Ep)

                ddpdy3(Sc:Ec) = (dPhat(Sc:Ec,j , k) - dPhat(Sc:Ec,jm, k))/dmx2(j )
                ddpdy4(Sc:Ec) = (dPhat(Sc:Ec,jp, k) - dPhat(Sc:Ec,j , k))/dmx2(jp)
                
                ddpdz5(Sc:Ec) = (dPhat(Sc:Ec,j , k) - dPhat(Sc:Ec,j, km))/dmx3(k )
                ddpdz6(Sc:Ec) = (dPhat(Sc:Ec,j, kp) - dPhat(Sc:Ec,j , k))/dmx3(kp)

                invRho1(Sc:Ec) = 0.5/dmx1(Sc:Ec) *(dx1(Sm:Em)*invRho(Sc:Ec,j ,k ) + dx1(Sc:Ec)*invRho(Sm:Em,j ,k ))
                invRho2(Sc:Ec) = 0.5/dmx1(Sp:Ep)*(dx1(Sp:Ep)*invRho(Sc:Ec,j ,k ) + dx1(Sc:Ec)*invRho(Sp:Ep,j ,k ))
                invRho3(Sc:Ec) = 0.5/dmx2(j) *(dx2(jm)*invRho(Sc:Ec,j ,k ) + dx2(j)*invRho(Sc:Ec,jm,k ))
                invRho4(Sc:Ec) = 0.5/dmx2(jp)*(dx2(jp)*invRho(Sc:Ec,j ,k ) + dx2(j)*invRho(Sc:Ec,jp,k ))
                invRho5(Sc:Ec) = 0.5/dmx3(k) *(dx3(km)*invRho(Sc:Ec,j ,k ) + dx3(K)*invRho(Sc:Ec,j ,km))
                invRho6(Sc:Ec) = 0.5/dmx3(kp)*(dx3(kp)*invRho(Sc:Ec,j ,k ) + dx3(K)*invRho(Sc:Ec,j ,kp))

                ExtraTerm(Sc:Ec) = ((1.d0-invRho2)*ddpdx2 - (1.d0-invRho1)*ddpdx1)/dx1(Sc:Ec)    &
                                 + ((1.d0-invRho4)*ddpdy4 - (1.d0-invRho3)*ddpdy3)/dx2(j)        &
                                 + ((1.d0-invRho6)*ddpdz6 - (1.d0-invRho5)*ddpdz5)/dx3(k) 

                PRHS(Sc:Ec,j,k) = (DivUm(Sc:Ec) - cbc(Sc:Ec))/dt/Cmp + ExtraTerm(Sc:Ec)
                
            end do
        end do
        
        deallocate(DivUm,cbc)
        deallocate(ddpdx1,ddpdx2,ddpdy3,ddpdy4,ddpdz5,ddpdz6)
        deallocate(invRho1,invRho2,invRho3,invRho4,invRho5,invRho6)
        deallocate(ExtraTerm)
    end subroutine mpi_pressure_RHS

    subroutine mpi_pressure_Poisson_FFT(dx2,dmx2,comm_1d_x1, comm_1d_x2, comm_1d_x3)
        use mpi_topology, only : cart_comm_1d
        use PaScaL_TDMA
        use MPI
        implicit none
        include 'fftw3.f'
        double precision :: dx2(0:nn2sb),dmx2(0:nn2sb)
        type(cart_comm_1d), intent(in)  :: comm_1d_x1, comm_1d_x2, comm_1d_x3
        type(ptdma_plan_many)     :: ptdma_plan

        double precision, allocatable, dimension(:,:,:) :: RHS_Iline,tmp
        complex(8), allocatable, dimension(:,:,:) :: RHSIhat_Iline,RHSIhat,RHSIhat_Kline,RHSIKhat_Kline

        double precision, allocatable, dimension(:,:) :: FFT_f_r2c_in,FFT_f_c2r_out
        complex(8), allocatable, dimension(:,:) :: FFT_f_r2c_out,FFT_f_c2r_in
        complex(8), allocatable, dimension(:,:) :: FFT_f_c2c_in,FFT_f_c2c_out

        double precision, allocatable, dimension(:,:,:) :: Am_r,Ac_r,Ap_r,Be_r
        double precision, allocatable, dimension(:,:,:) :: Am_c,Ac_c,Ap_c,Be_c
        double precision, allocatable, dimension(:) :: dzk2
        integer :: nh3m
        double precision :: ddx,ddz,Lz,Lx,dxk2
        double precision :: fft_amj,fft_apj,fft_acj

        complex(8), parameter :: II = DCMPLX(0.d0, 1.d0)
        integer(8) :: plan1,plan2
        
        integer :: i,j,k,im,jm,km,jp,isub,ierr        
        double precision :: Pbc_a,Pbc_b

        double precision :: AVERsub, AVERmpi_I,AVERmpi_J,AVERmpi_K,AVERmpi
        
        !== Alltoall C to I
        allocate(RHS_Iline(1:n1-1,1:nn2sb-1,1:n3mssub))
        call MPI_ALLTOALLW(PRHS     , countsendI, countdistI, ddtype_C2I_dble   &
                          ,RHS_Iline, countsendI, countdistI, ddtype_I2C_dble   &
                          ,comm_1d_x1%mpi_comm, ierr         )
        deallocate(PRHS)
            

            !== FFT Forward I
            allocate(FFT_f_r2c_in(1:n1-1,1:nn2sb-1),FFT_f_r2c_out(1:wv1,1:nn2sb-1))
            allocate(RHSIhat_Iline(1:wv1,1:nn2sb-1,1:n3mssub))
            ! FFT I: RHS_Iline>FFT_f_r2c_in --> FFT_f_r2c_out>RHSIhat_Iline
            call dfftw_plan_dft_r2c_1d(plan1,n1-1,FFT_f_r2c_in(:,1),FFT_f_r2c_out(:,1),FFTW_ESTIMATE)
            do k=1,n3mssub
                FFT_f_r2c_in(1:n1-1,1:nn2sb-1)=RHS_Iline(1:n1-1,1:nn2sb-1,k)
                    do j=1,nn2sb-1
                        call dfftw_execute_dft_r2c(plan1,FFT_f_r2c_in(:,j),FFT_f_r2c_out(:,j))
                    enddo
                RHSIhat_Iline(1:wv1,1:nn2sb-1,k)=FFT_f_r2c_out(1:wv1,1:nn2sb-1)
            enddo
            call dfftw_destroy_plan(plan1)
            deallocate(RHS_Iline)
            deallocate(FFT_f_r2c_in,FFT_f_r2c_out)
                
            
                !== Alltoall I to C
                allocate(RHSIhat(1:wv1sub,1:nn2sb-1,1:n3msub))
                call MPI_ALLTOALLW(RHSIhat_Iline, countsendI, countdistI, ddtype_I2C_cplx   &
                                ,RHSIhat      , countsendI, countdistI, ddtype_C2I_cplx   &
                                ,comm_1d_x1%mpi_comm, ierr         )
                deallocate(RHSIhat_Iline)
                
                    
                    !== Alltoall C to K
                    allocate(RHSIhat_Kline(1:wv1ssub,1:nn2sb-1,1:n3-1))
                    call MPI_ALLTOALLW(RHSIhat      , countsendK, countdistK, ddtype_C2K_cplx   &
                                    ,RHSIhat_Kline, countsendK, countdistK, ddtype_K2C_cplx   &
                                    ,comm_1d_x3%mpi_comm, ierr         )
                    deallocate(RHSIhat)


                        !== FFT Forward K
                        allocate(FFT_f_c2c_in(1:n3-1,1:wv1ssub),FFT_f_c2c_out(1:n3-1,1:wv1ssub))
                        allocate(RHSIKhat_Kline(1:n3-1,1:wv1ssub,1:nn2sb-1))
                        ! FFT K: RHSIhat_Kline>FFT_f_c2c_in --> FFT_f_c2c_out>RHSIKhat_Kline
                        call dfftw_plan_dft_1d(plan1,n3-1,FFT_f_c2c_in(:,1),FFT_f_c2c_out(:,1),fftw_forward,FFTW_ESTIMATE)
                        do j=1,nn2sb-1
                            do k=1,n3-1
                            do i=1,wv1ssub
                                FFT_f_c2c_in(k,i)=RHSIhat_Kline(i,j,k)
                            enddo
                            enddo
                                do i=1,wv1ssub
                                    call dfftw_execute_dft(plan1,  FFT_f_c2c_in(:,i),  FFT_f_c2c_out(:,i))
                                enddo
                            RHSIKhat_Kline(1:n3-1,1:wv1ssub,j)=FFT_f_c2c_out(1:n3-1,1:wv1ssub)
                        enddo
                        call dfftw_destroy_plan(plan1)
                        deallocate(RHSIhat_Kline)
                        deallocate(FFT_f_c2c_in,FFT_f_c2c_out)

                            !== TDMA RHS=RHSIKhat_Kline
                            Lx=L1;  ddx=Lx/dble(n1-1)
                            Lz=L3;  ddz=Lz/dble(n3-1)
                            nh3m=int((n3-1)/2)  
                            allocate(Am_r(1:n3-1,1:wv1ssub,nn2sb-1),Ac_r(1:n3-1,1:wv1ssub,nn2sb-1))
                            allocate(Ap_r(1:n3-1,1:wv1ssub,nn2sb-1),Be_r(1:n3-1,1:wv1ssub,nn2sb-1))
                            allocate(Am_c(1:n3-1,1:wv1ssub,nn2sb-1),Ac_c(1:n3-1,1:wv1ssub,nn2sb-1))
                            allocate(Ap_c(1:n3-1,1:wv1ssub,nn2sb-1),Be_c(1:n3-1,1:wv1ssub,nn2sb-1))
                            allocate(dzk2(1:n3-1))                          
                            do k = 1, n3-1
                                if(k <= nh3m+1)then
                                    km = k - 1
                                else
                                    km = N3m - k + 1
                                end if
                                dzk2(k) = 2.*(1.-cos(2.*PI*dble(km)*ddz/Lz))/(ddz*ddz) 
                            enddo

                            do j = 1, nn2sb-1     
                                jp = j + 1   
                                fft_amj = 1./dx2(j)/dmx2(j )
                                fft_apj = 1./dx2(j)/dmx2(jp)    
                                if(comm_1d_x2%myrank==0.and.j==1  )  fft_amj = 0.d0
                                if(comm_1d_x2%myrank==comm_1d_x2%nprocs-1.and.j==nn2sb-1 )  fft_apj = 0.d0

                                fft_acj =  - fft_amj - fft_apj
                        
                                do isub = 1, wv1ssub
                                    i = wv1ssubA+isub-1
                                    im = i - 1
                                    dxk2 = 2.*(1.-cos(2.*PI*dble(im)*ddx/Lx))/(ddx*ddx)
                                        
                                    ! Define the RHS for both 'real' and 'complex' 
                                    Be_r(1:n3-1,isub,j) = dble(RHSIKhat_Kline(1:n3-1,isub,j))
                                    Be_c(1:n3-1,isub,j) = aimag(RHSIKhat_Kline(1:n3-1,isub,j))
                        
                                    Am_r(1:n3-1,isub,j) = fft_amj
                                    Ac_r(1:n3-1,isub,j) = fft_acj - dxk2 - dzk2(1:n3-1)
                                    Ap_r(1:n3-1,isub,j) = fft_apj
                        
                                    Am_c(1:n3-1,isub,j) = fft_amj
                                    Ac_c(1:n3-1,isub,j) = fft_acj - dxk2 - dzk2(1:n3-1)
                                    Ap_c(1:n3-1,isub,j) = fft_apj
                                    
                                end do
                            end do

                            if(comm_1d_x1%myrank==0.and.comm_1d_x2%myrank==0.and.comm_1d_x3%myrank==0) then
                                Am_r(1,1,1) = 0.d0
                                Ac_r(1,1,1) = 1.d0
                                Ap_r(1,1,1) = 0.d0
                                Be_r(1,1,1) = 0.d0
                            endif
                        
                            call PaScaL_TDMA_plan_many_create(ptdma_plan, (n3-1)*(wv1ssub), comm_1d_x2%myrank, comm_1d_x2%nprocs, comm_1d_x2%mpi_comm)
                            call PaScaL_TDMA_many_solve(ptdma_plan,Am_r,Ac_r,Ap_r,Be_r,(n3-1)*(wv1ssub),nn2sb-1)
                            call PaScaL_TDMA_many_solve(ptdma_plan,Am_c,Ac_c,Ap_c,Be_c,(n3-1)*(wv1ssub),nn2sb-1)
                            call PaScaL_TDMA_plan_many_destroy(ptdma_plan,comm_1d_x2%nprocs)
                        
                            do j = 1, nn2sb-1   
                                if(wv1ssubA==1) then
                                    Be_c(1,1,j)= 0.d0;
                                    Be_c(nh3m+1,1,j)= 0.d0;
                                else if(wv1ssubB==wv1) then
                                    Be_c(1,wv1ssub,j)= 0.d0;
                                    Be_c(nh3m+1,wv1ssub,j)= 0.d0;
                                endif
                                RHSIKhat_Kline(1:n3-1,1:wv1ssub,j) = DCMPLX(Be_r(1:n3-1,1:wv1ssub,j), Be_c(1:n3-1,1:wv1ssub,j))
                            end do
                            deallocate(dzk2)
                            deallocate(Am_r,Ac_r,Ap_r,Be_r)
                            deallocate(Am_c,Ac_c,Ap_c,Be_c)

                        !== FFT backward K
                        allocate(FFT_f_c2c_out(1:n3-1,1:wv1ssub),FFT_f_c2c_in(1:n3-1,1:wv1ssub))
                        allocate(RHSIhat_Kline(1:wv1ssub,1:nn2sb-1,1:n3-1))
                        ! FFT back K: RHSIKhat_Kline>FFT_f_c2c_in --> FFT_f_c2c_out>RHSIhat_Kline
                        call dfftw_plan_dft_1d(plan2,n3-1,FFT_f_c2c_in(:,1),FFT_f_c2c_out(:,1),fftw_backward,fftw_estimate)
                        do j=1,nn2sb-1
                            FFT_f_c2c_in(1:n3-1,1:wv1ssub)=RHSIKhat_Kline(1:n3-1,1:wv1ssub,j)
                                do i=1,wv1ssub
                                    call dfftw_execute_dft(plan2,  FFT_f_c2c_in(:,i),  FFT_f_c2c_out(:,i))
                                enddo
                            do k=1,n3-1
                            do i=1,wv1ssub
                                RHSIhat_Kline(i,j,k)=FFT_f_c2c_out(k,i)/dble(n3-1)
                            enddo
                            enddo
                        enddo
                        call dfftw_destroy_plan(plan2)
                        deallocate(RHSIKhat_Kline)
                        deallocate(FFT_f_c2c_out,FFT_f_c2c_in)

                        
                    !== Alltoall K to C
                    allocate(RHSIhat(1:wv1sub,1:nn2sb-1,1:n3msub));
                    call MPI_ALLTOALLW(RHSIhat_Kline, countsendK, countdistK, ddtype_K2C_cplx   &
                                    ,RHSIhat      , countsendK, countdistK, ddtype_C2K_cplx   &
                                    ,comm_1d_x3%mpi_comm, ierr         )
                    deallocate(RHSIhat_Kline)
                    

                !== Alltoall C to I
                allocate(RHSIhat_Iline(1:wv1,1:nn2sb-1,1:n3mssub))
                call MPI_ALLTOALLW(RHSIhat      , countsendI, countdistI, ddtype_C2I_cplx   &
                                ,RHSIhat_Iline, countsendI, countdistI, ddtype_I2C_cplx   &
                                ,comm_1d_x1%mpi_comm, ierr         )
                deallocate(RHSIhat)

            !== FFT backward I
            allocate(FFT_f_c2r_out(1:n1-1,1:nn2sb-1),FFT_f_c2r_in(1:wv1,1:nn2sb-1))
            allocate(RHS_Iline(1:n1-1,1:nn2sb-1,1:n3mssub))
            ! FFT I: RHSIhat_Iline>FFT_f_c2r_in --> FFT_f_c2r_out>RHS_Iline
            call dfftw_plan_dft_c2r_1d(plan2,n1-1,FFT_f_c2r_in(:,1),FFT_f_c2r_out(:,1),FFTW_ESTIMATE)
            do k=1,n3mssub
                FFT_f_c2r_in(1:wv1,1:nn2sb-1)=RHSIhat_Iline(1:wv1,1:nn2sb-1,k)
                    do j=1,nn2sb-1
                        call dfftw_execute_dft_c2r(plan2,FFT_f_c2r_in(:,j),FFT_f_c2r_out(:,j))
                    enddo
                RHS_Iline(1:n1-1,1:nn2sb-1,k)=FFT_f_c2r_out(1:n1-1,1:nn2sb-1)/dble(n1-1)
            enddo
            call dfftw_destroy_plan(plan2)
            deallocate(RHSIhat_Iline)
            deallocate(FFT_f_c2r_in,FFT_f_c2r_out)

        !== Alltoall I to C
        allocate(tmp(1:nn1sb-1,1:nn2sb-1,1:nn3sb-1))
        call MPI_ALLTOALLW(RHS_Iline, countsendI, countdistI, ddtype_I2C_dble   &
                          ,tmp      , countsendI, countdistI, ddtype_C2I_dble   &
                          ,comm_1d_x1%mpi_comm, ierr         )
        deallocate(RHS_Iline)

        AVERsub=0.d0
        do k=1,nn3sb-1
            AVERsub = AVERsub+sum(tmp(1:nn1sb-1,1:nn2sb-1,k))
        enddo
        
        AVERmpi_I=0.d0;AVERmpi_J=0.d0;AVERmpi_K=0.d0;AVERmpi=0.d0;
        call MPI_ALLREDUCE(AVERsub  , AVERmpi_I, 1, MPI_DOUBLE, MPI_SUM, comm_1d_x1%mpi_comm, ierr)
        call MPI_ALLREDUCE(AVERmpi_I, AVERmpi_K, 1, MPI_DOUBLE, MPI_SUM, comm_1d_x2%mpi_comm, ierr)
        call MPI_ALLREDUCE(AVERmpi_K, AVERmpi  , 1, MPI_DOUBLE, MPI_SUM, comm_1d_x3%mpi_comm, ierr)    
        AVERmpi=AVERmpi/dble(n1-1)/dble(n2-1)/dble(n3-1)

        allocate(dP(0:nn1sb,0:nn2sb,0:nn3sb))
        dP(1:nn1sb-1,1:nn2sb-1,1:nn3sb-1)=tmp(1:nn1sb-1,1:nn2sb-1,1:nn3sb-1)- AVERmpi
        deallocate(tmp)
        
        if(comm_1d_x2%myrank == 0) then
            Pbc_a = (dmx2(1 )+dmx2(2  ))**2./((dmx2(1 )+dmx2(2  ))**2.-dmx2(1 )**2.)
            Pbc_b = dmx2(1 )**2.           /((dmx2(1 )+dmx2(2  ))**2.-dmx2(1 )**2.)
            
            do k=1,nn3sb
                dP(1:nn1sb-1,0, k) = Pbc_a*dP(1:nn1sb-1,1,   k) - Pbc_b*dP(1:nn1sb-1,2,   k)
            enddo
            
        endif
        if(comm_1d_x2%myrank==comm_1d_x2%nprocs-1) then
            Pbc_a = (dmx2(nn2sb)+dmx2(nn2sb-1))**2./((dmx2(nn2sb)+dmx2(nn2sb-1))**2.-dmx2(nn2sb)**2.)
            Pbc_b = dmx2(nn2sb)**2.           /((dmx2(nn2sb)+dmx2(nn2sb-1))**2.-dmx2(nn2sb)**2.)
            
            do k=1,nn3sb
                dP(1:nn1sb-1,nn2sb,k) = Pbc_a*dP(1:nn1sb-1,nn2sb-1,k) - Pbc_b*dP(1:nn1sb-1,nn2sb-2,k)
            enddo
            
        endif
        
    end subroutine mpi_pressure_Poisson_FFT

    subroutine mpi_pressure_Projection(invRho,U,V,W,P,dx1,dx2,dx3,dmx1,dmx2,dmx3,Jrank,Jnprocs)
        implicit none
        double precision, dimension(0:nn1sb,0:nn2sb,0:nn3sb) :: invRho,U,V,W,P
        double precision :: dx1(0:nn1sb),dx2(0:nn2sb),dx3(0:nn3sb)
        double precision :: dmx1(0:nn1sb),dmx2(0:nn2sb),dmx3(0:nn3sb)
        integer :: Jrank,Jnprocs
        integer :: i,j,k
        integer :: kp,km  
        integer :: jp,jm,jvm,jvp
        integer :: Sc,Ec
        integer :: Sm,Em
        integer :: Sp,Ep

        double precision, allocatable, dimension(:) :: invRhoc        
        double precision :: Pbc_a,Pbc_b

        Sc=1  ;Ec=nn1sb-1
        Sm=1-1;Em=nn1sb-1-1
        Sp=1+1;Ep=nn1sb-1+1

        allocate(invRhoc(Sc:Ec))
        
        do k = 1, nn3sb-1
            kp = k+1
            km = k-1
            do j = 2, nn2sb-1
                jp = j + 1
                jm = j - 1

                invRhoc(Sc:Ec) = 0.5d0/dmx1(Sc:Ec)*(dx1(Sc:Ec)*invRho(Sm:Em,j ,k ) + dx1(Sm:Em)*invRho(Sc:Ec,j ,k ))
                U(Sc:Ec,j ,k ) = U(Sc:Ec,j ,k )     &
                               - dt*Cmp*( (dP(Sc:Ec,j ,k ) - dP(Sm:Em,j , k))/dmx1(Sc:Ec)                              &
                                         +(invRhoc(Sc:Ec) - 1.d0)*(dPhat(Sc:Ec,j ,k )- dPhat(Sm:Em,j ,k ))/dmx1(Sc:Ec) )
                
                invRhoc(Sc:Ec) = 0.5d0/dmx2(j)*(dx2(j)*invRho(Sc:Ec,jm,k ) + dx2(jm)*invRho(Sc:Ec,j ,k ))
                V(Sc:Ec,j ,k ) = V(Sc:Ec,j ,k )     &
                               - dt*Cmp*( (dP(Sc:Ec,j ,k ) - dP(Sc:Ec,jm,k ))/dmx2(j)                                  &
                                         +(invRhoc(Sc:Ec) - 1.d0)*(dPhat(Sc:Ec,j ,k ) - dPhat(Sc:Ec,jm,k ))/dmx2(j)    )
                
                invRhoc(Sc:Ec) = 0.5d0/dmx3(k)*(dx3(k)*invRho(Sc:Ec,j ,km) + dx3(km)*invRho(Sc:Ec,j ,k ))
                W(Sc:Ec,j ,k ) = W(Sc:Ec,j ,k )     &
                               - dt*Cmp*( (dP(Sc:Ec,j ,k ) - dP(Sc:Ec, j, km))/dmx3(k)                                 &
                                         +(invRhoc(Sc:Ec) - 1.d0)*(dPhat(Sc:Ec,j ,k )- dPhat(Sc:Ec,j ,km))/dmx3(k)     )
                P(Sc:Ec,j,k) = P(Sc:Ec,j,k) + dP(Sc:Ec,j,k)

                dP0(Sc:Ec,j,k)=(dPhat(Sc:Ec,j,k)+dP0(Sc:Ec,j,k))/2.d0
                dPhat(Sc:Ec,j,k)=2.d0*dP(Sc:Ec,j,k)-dP0(Sc:Ec,j,k)
            enddo

            j=1
                jp = j + 1
                jm = j - 1

                invRhoc(Sc:Ec) = 0.5d0/dmx1(Sc:Ec)*(dx1(Sc:Ec)*invRho(Sm:Em,j ,k ) + dx1(Sm:Em)*invRho(Sc:Ec,j ,k ))
                U(Sc:Ec,j ,k ) = U(Sc:Ec,j ,k )     &
                               - dt*Cmp*( (dP(Sc:Ec,j ,k ) - dP(Sm:Em,j , k))/dmx1(Sc:Ec)                              &
                                         +(invRhoc(Sc:Ec) - 1.d0)*(dPhat(Sc:Ec,j ,k )- dPhat(Sm:Em,j ,k ))/dmx1(Sc:Ec) )

                invRhoc(Sc:Ec) = 0.5d0/dmx2(j)*(dx2(j)*invRho(Sc:Ec,jm,k ) + dx2(jm)*invRho(Sc:Ec,j ,k ))
                V(Sc:Ec,j ,k ) = V(Sc:Ec,j ,k )     &
                               - dt*Cmp*( (dP(Sc:Ec,j ,k ) - dP(Sc:Ec,jm,k ))/dmx2(j)                                  &
                                         +(invRhoc(Sc:Ec) - 1.d0)*(dPhat(Sc:Ec,j ,k ) - dPhat(Sc:Ec,jm,k ))/dmx2(j)    )*dble(2-j_indexS)
                
                invRhoc(Sc:Ec) = 0.5d0/dmx3(k)*(dx3(k)*invRho(Sc:Ec,j ,km) + dx3(km)*invRho(Sc:Ec,j ,k ))
                W(Sc:Ec,j ,k ) = W(Sc:Ec,j ,k )     &
                               - dt*Cmp*( (dP(Sc:Ec,j ,k ) - dP(Sc:Ec, j, km))/dmx3(k)                                 &
                                         +(invRhoc(Sc:Ec) - 1.d0)*(dPhat(Sc:Ec,j ,k )- dPhat(Sc:Ec,j ,km))/dmx3(k)     )
                P(Sc:Ec,j,k) = P(Sc:Ec,j,k) + dP(Sc:Ec,j,k)

                dP0(Sc:Ec,j,k)=(dPhat(Sc:Ec,j,k)+dP0(Sc:Ec,j,k))/2.d0
                dPhat(Sc:Ec,j,k)=2.d0*dP(Sc:Ec,j,k)-dP0(Sc:Ec,j,k)
                 
            if(Jrank == 0) then
                Pbc_a = (dmx2(1 )+dmx2(2  ))**2./((dmx2(1 )+dmx2(2  ))**2.-dmx2(1 )**2.)
                Pbc_b = dmx2(1 )**2.           /((dmx2(1 )+dmx2(2  ))**2.-dmx2(1 )**2.)
                
                P(Sc:Ec,0, k) = Pbc_a*P(Sc:Ec,1,   k) - Pbc_b*P(Sc:Ec,2,   k)
                
            endif
            if(Jrank==Jnprocs-1) then
                Pbc_a = (dmx2(nn2sb)+dmx2(nn2sb-1))**2./((dmx2(nn2sb)+dmx2(nn2sb-1))**2.-dmx2(nn2sb)**2.)
                Pbc_b = dmx2(nn2sb)**2.           /((dmx2(nn2sb)+dmx2(nn2sb-1))**2.-dmx2(nn2sb)**2.)
                
                P(Sc:Ec,nn2sb,k) = Pbc_a*P(Sc:Ec,nn2sb-1,k) - Pbc_b*P(Sc:Ec,nn2sb-2,k)
                
            endif
        enddo

        deallocate(invRhoc)
        deallocate(dP)
    end subroutine mpi_pressure_Projection

end module mpi_pressure

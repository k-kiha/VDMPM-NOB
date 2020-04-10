module mpi_subdomain
    use MPI
    use global
    use mpi_topology

    implicit none
    
    !> @{ Grid numbers in the subdomain
    integer, public :: n1sub,n2sub,n3sub
    !> @}
    !> @{ Grid indices of the assigned range
    integer, public :: ista, iend, jsta, jend, ksta, kend
    !> @}

    !> @{ Coordinates of grid points in the subdomain
    double precision, allocatable, dimension(:), public     :: x1_sub, x2_sub, x3_sub
    !> @}
    !> @{ Grid lengths in the subdomain
    double precision, allocatable, dimension(:), public     ::  dx1_sub,  dx2_sub,  dx3_sub
    double precision, allocatable, dimension(:), public     :: dmx1_sub, dmx2_sub, dmx3_sub
    !> @}

    !> @{ Derived datatype for communication between x-neighbor subdomains
    integer :: ddtype_sendto_E, ddtype_recvfrom_W, ddtype_sendto_W, ddtype_recvfrom_E
    !> @}
    !> @{ Derived datatype for communication between y-neighbor subdomains
    integer :: ddtype_sendto_N, ddtype_recvfrom_S, ddtype_sendto_S, ddtype_recvfrom_N
    !> @}
    !> @{ Derived datatype for communication between z-neighbor subdomains
    integer :: ddtype_sendto_F, ddtype_recvfrom_B, ddtype_sendto_B, ddtype_recvfrom_F
    !> @}
    
    private :: subdomain_para_range

    public  :: mpi_subdomain_make
    public  :: mpi_subdomain_clean
    public  :: mpi_subdomain_mesh
    public  :: mpi_subdomain_DDT_ghostcell
    public  :: mpi_subdomain_DDT_FFT
    public  :: mpi_subdomain_ghostcell_update

    contains

    subroutine mpi_subdomain_make()
        implicit none

        ! Assigning grid numbers and grid indices of my subdomain.
        call subdomain_para_range(1, n1-1, comm_1d_x1%nprocs, comm_1d_x1%myrank, ista, iend)
        n1sub = iend - ista + 2
        call subdomain_para_range(1, n2-1, comm_1d_x2%nprocs, comm_1d_x2%myrank, jsta, jend)
        n2sub = jend - jsta + 2
        call subdomain_para_range(1, n3-1, comm_1d_x3%nprocs, comm_1d_x3%myrank, ksta, kend)
        n3sub = kend - ksta + 2
        
        ! Allocate subdomain variables.
        allocate( x1_sub(0:n1sub), dmx1_sub(0:n1sub), dx1_sub(0:n1sub))
        allocate( x2_sub(0:n2sub), dmx2_sub(0:n2sub), dx2_sub(0:n2sub))
        allocate( x3_sub(0:n3sub), dmx3_sub(0:n3sub), dx3_sub(0:n3sub))
        
    end subroutine mpi_subdomain_make

    subroutine mpi_subdomain_clean
        implicit none
        
        deallocate(x1_sub, dmx1_sub, dx1_sub)
        deallocate(x2_sub, dmx2_sub, dx2_sub)
        deallocate(x3_sub, dmx3_sub, dx3_sub)

    end subroutine mpi_subdomain_clean

    subroutine mpi_subdomain_mesh()
        implicit none

        ! integer :: i,j,k,level
        integer :: i,j,k
        double precision :: temp
        integer :: request_S2E, request_S2W, ierr
        integer :: STATUS(MPI_STATUS_SIZE)
        
        x1_sub(0:n1sub)=0.d0; dmx1_sub(0:n1sub)=0.d0; dx1_sub(0:n1sub)=0.d0;
        x2_sub(0:n2sub)=0.d0; dmx2_sub(0:n2sub)=0.d0; dx2_sub(0:n2sub)=0.d0;
        x3_sub(0:n3sub)=0.d0; dmx3_sub(0:n3sub)=0.d0; dx3_sub(0:n3sub)=0.d0;
        ! MGDX(0:N1,MaxLevel)=0.d0
        ! MGDY(0:N2,MaxLevel)=0.d0
        ! MGDZ(0:N3,MaxLevel)=0.d0
        
        if(comm_1d_x1%myrank==0) x1_sub(1)=x1_start;x1_sub(0)=x1_sub(1)
        if(comm_1d_x2%myrank==0) x2_sub(1)=x2_start;x2_sub(0)=x2_sub(1)
        if(comm_1d_x3%myrank==0) x3_sub(1)=x3_start;x3_sub(0)=x3_sub(1)

        if(period(0)==.true.) then
            if(comm_1d_x1%myrank==0) x1_sub(0)= -(x1_end-x1_start)/dble(N1m)
            ! if(comm_1d_x1%myrank==comm_1d_x1%nprocs-1) x1_sub(0)=
        endif
        if(period(1)==.true.) then
            if(comm_1d_x2%myrank==0) x2_sub(0)= -(x2_end-x2_start)/dble(N2m)
            ! if(comm_1d_x2%myrank==comm_1d_x2%nprocs-1) x2_sub(0)=
        endif
        if(period(2)==.true.) then
            if(comm_1d_x3%myrank==0) x3_sub(0)= -(x3_end-x3_start)/dble(N3m)
            ! if(comm_1d_x3%myrank==comm_1d_x3%nprocs-1) x3_sub(0)=
        endif
        
        ! X-DIRECTION
        if(UNIFORM1 == 1) then 
            temp = (x1_end-x1_start)/dble(N1m)
            do i = ista-1, iend+1
                x1_sub(i-ista+1) = dble(i-1)*temp + x1_start
            end do

        else
            do i = ista-1, iend+1
                x1_sub(i-ista+1) = (x1_end-x1_start)*0.5*(1. + tanh(0.5*GAMMA1*(2.*dble(i-1)/dble(N1m)-1.0))/tanh(GAMMA1*0.5) ) + x1_start
            end do

        end if
        if(period(0)==.false. .and. comm_1d_x1%myrank==0) x1_sub(0)=x1_sub(1)

        do i = 1, n1sub-1
            dx1_sub(i) = x1_sub(i+1) - x1_sub(i)
        end do
        call MPI_ISEND(dx1_sub(n1sub-1),1, MPI_DOUBLE, comm_1d_x1%east_rank, 111, comm_1d_x1%mpi_comm, request_S2E, ierr)
        call MPI_IRECV(dx1_sub(0)      ,1, MPI_DOUBLE, comm_1d_x1%west_rank, 111, comm_1d_x1%mpi_comm, request_S2E, ierr)
        call MPI_WAIT(request_S2E,STATUS,ierr)
        call MPI_ISEND(dx1_sub(1)     ,1, MPI_DOUBLE, comm_1d_x1%west_rank, 111, comm_1d_x1%mpi_comm, request_S2W, ierr)
        call MPI_IRECV(dx1_sub(n1sub),1, MPI_DOUBLE, comm_1d_x1%east_rank, 111, comm_1d_x1%mpi_comm, request_S2W, ierr)
        call MPI_WAIT(request_S2W,STATUS,ierr)

        if(period(0)==.false. .and. comm_1d_x1%myrank==0) dx1_sub(0)=0.d0
        if(period(0)==.false. .and. comm_1d_x1%myrank==comm_1d_x1%nprocs-1) dx1_sub(n1sub)=0.d0

        do i = 1, n1sub-1
            dmx1_sub(i) = 0.5*(dx1_sub(i-1)+dx1_sub(i))
        end do
        if(period(0)==.false. .and. comm_1d_x1%myrank==comm_1d_x1%nprocs-1) dmx1_sub(n1sub)=0.5*(dx1_sub(n1sub-1)+dx1_sub(n1sub))

        call MPI_ISEND(dmx1_sub(n1sub-1),1, MPI_DOUBLE, comm_1d_x1%east_rank, 111, comm_1d_x1%mpi_comm, request_S2E, ierr)
        call MPI_IRECV(dmx1_sub(0)       ,1, MPI_DOUBLE, comm_1d_x1%west_rank, 111, comm_1d_x1%mpi_comm, request_S2E, ierr)
        call MPI_WAIT(request_S2E,STATUS,ierr)
        call MPI_ISEND(dmx1_sub(1)     ,1, MPI_DOUBLE, comm_1d_x1%west_rank, 111, comm_1d_x1%mpi_comm, request_S2W, ierr)
        call MPI_IRECV(dmx1_sub(n1sub),1, MPI_DOUBLE, comm_1d_x1%east_rank, 111, comm_1d_x1%mpi_comm, request_S2W, ierr)
        call MPI_WAIT(request_S2W,STATUS,ierr)

        if(period(0)==.true. .and. comm_1d_x1%myrank==0) dmx1_sub(0)=(x1_end-x1_start)/dble(N1m)
        if(period(0)==.true. .and. comm_1d_x1%myrank==comm_1d_x1%nprocs-1) dmx1_sub(n1sub)=(x1_end-x1_start)/dble(N1m)

        ! Y-DIRECTION
        if(UNIFORM2 == 1) then 
            temp = (x2_end-x2_start)/dble(N2m)
            do j = jsta-1, jend+1
                x2_sub(j-jsta+1) = dble(j-1)*temp + x2_start
            end do

        else
            do j = jsta-1, jend+1
                x2_sub(j-jsta+1) = (x2_end-x2_start)*0.5*(1. + tanh(0.5*GAMMA2*(2.*dble(j-1)/dble(N2m)-1.0))/tanh(GAMMA2*0.5) ) + x2_start
            end do
        end if
        if(period(1)==.false. .and. comm_1d_x2%myrank==0) x2_sub(0)=x2_sub(1)

        do j = 1, n2sub-1
            dx2_sub(j) = x2_sub(j+1) - x2_sub(j)
        end do
        call MPI_ISEND(dx2_sub(n2sub-1),1, MPI_DOUBLE, comm_1d_x2%east_rank, 111, comm_1d_x2%mpi_comm, request_S2E, ierr)
        call MPI_IRECV(dx2_sub(0)      ,1, MPI_DOUBLE, comm_1d_x2%west_rank, 111, comm_1d_x2%mpi_comm, request_S2E, ierr)
        call MPI_WAIT(request_S2E,STATUS,ierr)
        call MPI_ISEND(dx2_sub(1)     ,1, MPI_DOUBLE, comm_1d_x2%west_rank, 111, comm_1d_x2%mpi_comm, request_S2W, ierr)
        call MPI_IRECV(dx2_sub(n2sub),1, MPI_DOUBLE, comm_1d_x2%east_rank, 111, comm_1d_x2%mpi_comm, request_S2W, ierr)
        call MPI_WAIT(request_S2W,STATUS,ierr)

        if(period(1)==.false. .and. comm_1d_x2%myrank==0) dx2_sub(0)=0.d0
        if(period(1)==.false. .and. comm_1d_x2%myrank==comm_1d_x2%nprocs-1) dx2_sub(n2sub)=0.d0

        do j = 1, n2sub-1
            dmx2_sub(j) = 0.5*(dx2_sub(j-1)+dx2_sub(j))
        end do
        if(period(1)==.false. .and. comm_1d_x2%myrank==comm_1d_x2%nprocs-1) dmx2_sub(n2sub)=0.5*(dx2_sub(n2sub-1)+dx2_sub(n2sub))

        call MPI_ISEND(dmx2_sub(n2sub-1),1, MPI_DOUBLE, comm_1d_x2%east_rank, 111, comm_1d_x2%mpi_comm, request_S2E, ierr)
        call MPI_IRECV(dmx2_sub(0)      ,1, MPI_DOUBLE, comm_1d_x2%west_rank, 111, comm_1d_x2%mpi_comm, request_S2E, ierr)
        call MPI_WAIT(request_S2E,STATUS,ierr)
        call MPI_ISEND(dmx2_sub(1)     ,1, MPI_DOUBLE, comm_1d_x2%west_rank, 111, comm_1d_x2%mpi_comm, request_S2W, ierr)
        call MPI_IRECV(dmx2_sub(n2sub),1, MPI_DOUBLE, comm_1d_x2%east_rank, 111, comm_1d_x2%mpi_comm, request_S2W, ierr)
        call MPI_WAIT(request_S2W,STATUS,ierr)

        if(period(1)==.true. .and. comm_1d_x2%myrank==0) dmx2_sub(0)=(x2_end-x2_start)/dble(N2m)
        if(period(1)==.true. .and. comm_1d_x2%myrank==comm_1d_x2%nprocs-1) dmx2_sub(n2sub)=(x2_end-x2_start)/dble(N2m)


        ! Z-DIRECTION
        if(UNIFORM3 == 1) then 
            temp = (x3_end-x3_start)/dble(N3m)
            do k = ksta-1, kend+1
                x3_sub(k-ksta+1) = dble(k-1)*temp + x3_start
            end do

        else
            do k = ksta-1, kend+1
                x3_sub(k-ksta+1) = (x3_end-x3_start)*0.5*(1. + tanh(0.5*GAMMA3*(2.*dble(k-1)/dble(N3m)-1.0))/tanh(GAMMA3*0.5) ) + x3_start
            end do
        end if
        if(period(2)==.false. .and. comm_1d_x3%myrank==0) x3_sub(0)=x3_sub(1)

        do k = 1, n3sub-1
            dx3_sub(k) = x3_sub(k+1) - x3_sub(k)
        end do
        call MPI_ISEND(dx3_sub(n3sub-1),1, MPI_DOUBLE, comm_1d_x3%east_rank, 111, comm_1d_x3%mpi_comm, request_S2E, ierr)
        call MPI_IRECV(dx3_sub(0)      ,1, MPI_DOUBLE, comm_1d_x3%west_rank, 111, comm_1d_x3%mpi_comm, request_S2E, ierr)
        call MPI_WAIT(request_S2E,STATUS,ierr)
        call MPI_ISEND(dx3_sub(1)     ,1, MPI_DOUBLE, comm_1d_x3%west_rank, 111, comm_1d_x3%mpi_comm, request_S2W, ierr)
        call MPI_IRECV(dx3_sub(n3sub),1, MPI_DOUBLE, comm_1d_x3%east_rank, 111, comm_1d_x3%mpi_comm, request_S2W, ierr)
        call MPI_WAIT(request_S2W,STATUS,ierr)

        if(period(2)==.false. .and. comm_1d_x3%myrank==0) dx3_sub(0)=0.d0
        if(period(2)==.false. .and. comm_1d_x3%myrank==comm_1d_x3%nprocs-1) dx3_sub(n3sub)=0.d0

        do k = 1, n3sub-1
            dmx3_sub(k) = 0.5*(dx3_sub(k-1)+dx3_sub(k))
        end do
        if(period(2)==.false. .and. comm_1d_x3%myrank==comm_1d_x3%nprocs-1) dmx3_sub(n3sub)=0.5*(dx3_sub(n3sub-1)+dx3_sub(n3sub))
        call MPI_ISEND(dmx3_sub(n3sub-1),1, MPI_DOUBLE, comm_1d_x3%east_rank, 111, comm_1d_x3%mpi_comm, request_S2E, ierr)
        call MPI_IRECV(dmx3_sub(0)      ,1, MPI_DOUBLE, comm_1d_x3%west_rank, 111, comm_1d_x3%mpi_comm, request_S2E, ierr)
        call MPI_WAIT(request_S2E,STATUS,ierr)
        call MPI_ISEND(dmx3_sub(1)     ,1, MPI_DOUBLE, comm_1d_x3%west_rank, 111, comm_1d_x3%mpi_comm, request_S2W, ierr)
        call MPI_IRECV(dmx3_sub(n3sub),1, MPI_DOUBLE, comm_1d_x3%east_rank, 111, comm_1d_x3%mpi_comm, request_S2W, ierr)
        call MPI_WAIT(request_S2W,STATUS,ierr)

        if(period(2)==.true. .and. comm_1d_x3%myrank==0) dmx3_sub(0)=(x3_end-x3_start)/dble(N3m)
        if(period(2)==.true. .and. comm_1d_x3%myrank==comm_1d_x3%nprocs-1) dmx3_sub(n3sub)=(x3_end-x3_start)/dble(N3m)


        ! ! FOR THE MULTIGRID PPOISSON SOLVER
        ! MGDX(0:N1,1:MaxLevel) = 0.
        ! MGDY(0:N2,1:MaxLevel) = 0.
        ! MGDZ(0:N3,1:MaxLevel) = 0.

        ! MGDX(0:N1,1) = DX(0:N1)
        ! MGDY(0:N2,1) = DY(0:N2)
        ! MGDZ(0:N3,1) = DZ(0:N3)

        ! !  GENERATE MGDX, MGDY, MGDZ
        ! do level = 2, MaxLevel
            
        !     do i = 1, N1m/(2**(level-1))
        !         MGDX(i, level) = MGDX(2*i-1, level-1) + MGDX(2*i, level-1)
        !     end do

        !     do j = 1, N2m/(2**(level-1))
        !         MGDY(j, level) = MGDY(2*j-1, level-1) + MGDY(2*j, level-1)
        !     end do    

        !     do k = 1, N3m/(2**(level-1))
        !         MGDZ(k, level) = MGDZ(2*k-1, level-1) + MGDZ(2*k, level-1)
        !     end do   

        ! end do
    
    end subroutine mpi_subdomain_mesh


    subroutine mpi_subdomain_DDT_ghostcell

        implicit none
        integer :: sizes(0:2), subsizes(0:2), starts(0:2), ierr     ! Local variables for MPI_Type_create_subarray

        ! ddtype sending data to east MPI process (x+ neighbor)
        sizes    = (/n1sub+1,n2sub+1,n3sub+1/)
        subsizes = (/      1,n2sub+1,n3sub+1/)
        starts   = (/n1sub-1,      0,      0/)
        call MPI_Type_create_subarray(3, sizes, subsizes, starts, MPI_ORDER_FORTRAN, &
                                    MPI_DOUBLE_PRECISION, ddtype_sendto_E, ierr)
        call MPI_Type_commit(ddtype_sendto_E,ierr)

        ! ddtype receiving data from west MPI process (x- neighbor)
        starts   = (/      0,      0,      0/)
        call MPI_Type_create_subarray(3, sizes, subsizes, starts, MPI_ORDER_FORTRAN, &
                                    MPI_DOUBLE_PRECISION, ddtype_recvfrom_W, ierr)
        call MPI_Type_commit(ddtype_recvfrom_W,ierr)

        ! ddtype sending data to west MPI process (x- neighbor)
        starts   = (/      1,      0,      0/)
        call MPI_Type_create_subarray(3, sizes, subsizes, starts, MPI_ORDER_FORTRAN, &
                                    MPI_DOUBLE_PRECISION, ddtype_sendto_W, ierr)
        call MPI_Type_commit(ddtype_sendto_W,ierr)

        ! ddtype receiving data from east MPI process (x+ neighbor)
        starts   = (/  n1sub,      0,      0/)
        call MPI_Type_create_subarray(3, sizes, subsizes, starts, MPI_ORDER_FORTRAN, &
                                    MPI_DOUBLE_PRECISION, ddtype_recvfrom_E, ierr)
        call MPI_Type_commit(ddtype_recvfrom_E,ierr)

        ! ddtype sending data to north MPI process (y+ neighbor)
        sizes    = (/n1sub+1,n2sub+1,n3sub+1/)
        subsizes = (/n1sub+1,      1,n3sub+1/)
        starts   = (/      0,n2sub-1,      0/)
        call MPI_Type_create_subarray(3, sizes, subsizes, starts, MPI_ORDER_FORTRAN, &
                                    MPI_DOUBLE_PRECISION, ddtype_sendto_N, ierr)
        call MPI_Type_commit(ddtype_sendto_N,ierr)

        ! ddtype receiving data from south MPI process (y- neighbor)
        starts   = (/      0,      0,      0/)
        call MPI_Type_create_subarray(3, sizes, subsizes, starts, MPI_ORDER_FORTRAN, &
                                    MPI_DOUBLE_PRECISION, ddtype_recvfrom_S, ierr)
        call MPI_Type_commit(ddtype_recvfrom_S,ierr)

        ! ddtype sending data to south MPI process (y- neighbor)
        starts   = (/      0,      1,      0/)
        call MPI_Type_create_subarray(3, sizes, subsizes, starts, MPI_ORDER_FORTRAN, &
                                    MPI_DOUBLE_PRECISION, ddtype_sendto_S, ierr)
        call MPI_Type_commit(ddtype_sendto_S,ierr)

        ! ddtype receiving data from north MPI process (y+ neighbor)
        starts   = (/      0,  n2sub,      0/)
        call MPI_Type_create_subarray(3, sizes, subsizes, starts, MPI_ORDER_FORTRAN, &
                                    MPI_DOUBLE_PRECISION, ddtype_recvfrom_N, ierr)
        call MPI_Type_commit(ddtype_recvfrom_N,ierr)

        ! ddtype sending data to forth MPI process (z+ neighbor)
        sizes    = (/n1sub+1,n2sub+1,n3sub+1/)
        subsizes = (/n1sub+1,n2sub+1,      1/)
        starts   = (/      0,      0,n3sub-1/)
        call MPI_Type_create_subarray(3, sizes, subsizes, starts, MPI_ORDER_FORTRAN, &
                                    MPI_DOUBLE_PRECISION, ddtype_sendto_F, ierr)
        call MPI_Type_commit(ddtype_sendto_F,ierr)

        ! ddtype receiving data from back MPI process (z- neighbor)
        starts   = (/      0,      0,      0/)
        call MPI_Type_create_subarray(3, sizes, subsizes, starts, MPI_ORDER_FORTRAN, &
                                    MPI_DOUBLE_PRECISION, ddtype_recvfrom_B, ierr)
        call MPI_Type_commit(ddtype_recvfrom_B,ierr)

        ! ddtype sending data to back MPI process (z- neighbor)
        starts   = (/      0,      0,      1/)
        call MPI_Type_create_subarray(3, sizes, subsizes, starts, MPI_ORDER_FORTRAN, &
                                    MPI_DOUBLE_PRECISION, ddtype_sendto_B, ierr)
        call MPI_Type_commit(  ddtype_sendto_B,ierr)

        ! ddtype receiving data from forth MPI process (z+ neighbor)
        starts   = (/      0,      0,  n3sub/)
        call MPI_Type_create_subarray(3, sizes, subsizes, starts, MPI_ORDER_FORTRAN, &
                                    MPI_DOUBLE_PRECISION, ddtype_recvfrom_F, ierr)
        call MPI_Type_commit(ddtype_recvfrom_F,ierr)

    end subroutine mpi_subdomain_DDT_ghostcell

    subroutine mpi_subdomain_DDT_FFT( n1msub,n3msub,n3mssub,wv1,wv1sub,wv1ssub,wv1ssubA,wv1ssubB  &
                                    , ddtype_C2I_dble,ddtype_C2I_cplx,ddtype_C2K_cplx &
                                    , ddtype_I2C_dble,ddtype_I2C_cplx,ddtype_K2C_cplx )
        implicit none
        integer :: n1msub,n3msub,n3mssub,wv1,wv1sub,wv1ssub,wv1ssubA,wv1ssubB
        integer :: ddtype_C2I_dble(0:comm_1d_x1%nprocs-1),ddtype_I2C_dble(0:comm_1d_x1%nprocs-1)
        integer :: ddtype_C2I_cplx(0:comm_1d_x1%nprocs-1),ddtype_I2C_cplx(0:comm_1d_x1%nprocs-1)
        integer :: ddtype_C2K_cplx(0:comm_1d_x3%nprocs-1),ddtype_K2C_cplx(0:comm_1d_x3%nprocs-1)
        
        integer :: i
        integer :: bigsize(3), subsize(3), start(3), ierr 
        integer :: indexA,indexB
        integer, allocatable, dimension(:) :: n3mssubAll,n1msubAll,wv1subAll,wv1ssubAll,n3msubAll

        allocate(n3mssubAll(0:comm_1d_x1%nprocs-1),n1msubAll(0:comm_1d_x1%nprocs-1),wv1subAll(0:comm_1d_x1%nprocs-1))
        allocate(wv1ssubAll(0:comm_1d_x3%nprocs-1),n3msubAll(0:comm_1d_x3%nprocs-1))
        
        n1msub=n1sub-1
        n3msub=n3sub-1
        call subdomain_para_range(1, n3msub, comm_1d_x1%nprocs, comm_1d_x1%myrank, indexA, indexB)
        n3mssub= indexB - indexA + 1     
        
        wv1=int((n1-1)/2)+1
        call subdomain_para_range(1, wv1, comm_1d_x1%nprocs, comm_1d_x1%myrank, indexA, indexB)
        wv1sub= indexB - indexA + 1

        call subdomain_para_range(1, wv1sub, comm_1d_x3%nprocs, comm_1d_x3%myrank, indexA, indexB)
        wv1ssub= indexB - indexA + 1
        
        do i=0,comm_1d_x1%nprocs-1
            call subdomain_para_range(1, n3msub, comm_1d_x1%nprocs, i, indexA, indexB)
            n3mssubAll(i)= indexB - indexA + 1        
            call subdomain_para_range(1, n1-1, comm_1d_x1%nprocs, i, indexA, indexB)
            n1msubAll(i)= indexB - indexA + 1        
            call subdomain_para_range(1, wv1, comm_1d_x1%nprocs, i, indexA, indexB)
            wv1subAll(i)= indexB - indexA + 1
        enddo
        
        do i=0,comm_1d_x3%nprocs-1             
            call subdomain_para_range(1, wv1sub, comm_1d_x3%nprocs, i, indexA, indexB)
            wv1ssubAll(i) =indexB - indexA + 1 
            call subdomain_para_range(1, n3-1, comm_1d_x3%nprocs, i, indexA, indexB)
            n3msubAll(i) =indexB - indexA + 1 
        enddo

        wv1ssubA=sum(wv1subAll(0:comm_1d_x1%myrank))-wv1subAll(comm_1d_x1%myrank)+1 &
                +sum(wv1ssubAll(0:comm_1d_x3%myrank))-wv1ssubAll(comm_1d_x3%myrank)
        wv1ssubB=wv1ssubA+wv1ssubAll(comm_1d_x3%myrank)-1
        
        do i=0,comm_1d_x1%nprocs-1
            bigsize(1) = n1msub
            bigsize(2) = n2sub-1
            bigsize(3) = n3msub
            subsize(1) = n1msub
            subsize(2) = n2sub-1
            subsize(3) = n3mssubAll(i)
            start(1) = 0
            start(2) = 0
            start(3) = sum(n3mssubAll(0:i)) - n3mssubAll(i)
                        
            call MPI_TYPE_CREATE_SUBARRAY( 3, bigsize, subsize, start, MPI_ORDER_FORTRAN    &
                                         , MPI_DOUBLE, ddtype_C2I_dble(i), ierr             )
            call MPI_TYPE_COMMIT(ddtype_C2I_dble(i),ierr)
                        
            bigsize(1) = n1-1
            bigsize(2) = n2sub-1
            bigsize(3) = n3mssub
            subsize(1) = n1msubAll(i)
            subsize(2) = n2sub-1
            subsize(3) = n3mssub
            start(1) = sum(n1msubAll(0:i)) - n1msubAll(i)
            start(2) = 0
            start(3) = 0
                        
            call MPI_TYPE_CREATE_SUBARRAY( 3, bigsize, subsize, start, MPI_ORDER_FORTRAN    &
                                         , MPI_DOUBLE, ddtype_I2C_dble(i), ierr     )
            call MPI_TYPE_COMMIT(ddtype_I2C_dble(i),ierr)
        enddo

        do i=0,comm_1d_x1%nprocs-1
            bigsize(1) = wv1sub
            bigsize(2) = n2sub-1
            bigsize(3) = n3msub
            subsize(1) = wv1sub
            subsize(2) = n2sub-1
            subsize(3) = n3mssubAll(i)
            start(1) = 0
            start(2) = 0
            start(3) = sum(n3mssubAll(0:i)) - n3mssubAll(i)
            call MPI_TYPE_CREATE_SUBARRAY( 3, bigsize, subsize, start, MPI_ORDER_FORTRAN    &
                                         , MPI_DOUBLE_COMPLEX, ddtype_C2I_cplx(i), ierr             )
            call MPI_TYPE_COMMIT(ddtype_C2I_cplx(i),ierr)
                                                
            bigsize(1) = wv1
            bigsize(2) = n2sub-1
            bigsize(3) = n3mssub
            subsize(1) = wv1subAll(i)
            subsize(2) = n2sub-1
            subsize(3) = n3mssub
            start(1) = sum(wv1subAll(0:i)) - wv1subAll(i)
            start(2) = 0
            start(3) = 0
            call MPI_TYPE_CREATE_SUBARRAY( 3, bigsize, subsize, start, MPI_ORDER_FORTRAN    &
                                         , MPI_DOUBLE_COMPLEX, ddtype_I2C_cplx(i), ierr     )
            call MPI_TYPE_COMMIT(ddtype_I2C_cplx(i),ierr)
        enddo

        do i=0,comm_1d_x3%nprocs-1
            bigsize(1) = wv1sub
            bigsize(2) = n2sub-1
            bigsize(3) = n3msub
            subsize(1) = wv1ssubAll(i)
            subsize(2) = n2sub-1
            subsize(3) = n3msub
            start(1) = sum(wv1ssubAll(0:i)) - wv1ssubAll(i)
            start(2) = 0
            start(3) = 0
            call MPI_TYPE_CREATE_SUBARRAY( 3, bigsize, subsize, start, MPI_ORDER_FORTRAN    &
                                         , MPI_DOUBLE_COMPLEX, ddtype_C2K_cplx(i), ierr     )
            call MPI_TYPE_COMMIT(ddtype_C2K_cplx(i),ierr)
                                                
            bigsize(1) = wv1ssub
            bigsize(2) = n2sub-1
            bigsize(3) = n3-1
            subsize(1) = wv1ssub
            subsize(2) = n2sub-1
            subsize(3) = n3msubAll(i)
            start(1) = 0
            start(2) = 0
            start(3) = sum(n3msubAll(0:i)) - n3msubAll(i)
            call MPI_TYPE_CREATE_SUBARRAY( 3, bigsize, subsize, start, MPI_ORDER_FORTRAN    &
                                         , MPI_DOUBLE_COMPLEX, ddtype_K2C_cplx(i), ierr     )
            call MPI_TYPE_COMMIT(ddtype_K2C_cplx(i),ierr)
        enddo
        deallocate(n3mssubAll,n1msubAll,wv1subAll,wv1ssubAll,n3msubAll)
    end subroutine mpi_subdomain_DDT_FFT

    subroutine mpi_subdomain_ghostcell_update(Value_sub)
        implicit none
        double precision, dimension(0:n1sub, 0:n2sub, 0:n3sub), intent(inout) :: Value_sub

        integer :: ierr
        integer :: request(4)

        ! Update the ghostcells in the x-direction using derived datatypes and subcommunicator.
        call MPI_Isend(Value_sub,1, ddtype_sendto_E  , comm_1d_x1%east_rank, 111, comm_1d_x1%mpi_comm, request(1), ierr)
        call MPI_Irecv(Value_sub,1, ddtype_recvfrom_W, comm_1d_x1%west_rank, 111, comm_1d_x1%mpi_comm, request(2), ierr)
        call MPI_Isend(Value_sub,1, ddtype_sendto_W  , comm_1d_x1%west_rank, 222, comm_1d_x1%mpi_comm, request(3), ierr)
        call MPI_Irecv(Value_sub,1, ddtype_recvfrom_E, comm_1d_x1%east_rank, 222, comm_1d_x1%mpi_comm, request(4), ierr)
        call MPI_Waitall(4, request, MPI_STATUSES_IGNORE, ierr)
        
        ! Update the ghostcells in the y-direction using derived datatypes and subcommunicator.
        call MPI_Isend(Value_sub,1, ddtype_sendto_N  , comm_1d_x2%east_rank, 111, comm_1d_x2%mpi_comm, request(1), ierr)
        call MPI_Irecv(Value_sub,1, ddtype_recvfrom_S, comm_1d_x2%west_rank, 111, comm_1d_x2%mpi_comm, request(2), ierr)
        call MPI_Isend(Value_sub,1, ddtype_sendto_S  , comm_1d_x2%west_rank, 222, comm_1d_x2%mpi_comm, request(3), ierr)
        call MPI_Irecv(Value_sub,1, ddtype_recvfrom_N, comm_1d_x2%east_rank, 222, comm_1d_x2%mpi_comm, request(4), ierr)
        call MPI_Waitall(4, request, MPI_STATUSES_IGNORE, ierr)

        ! Update the ghostcells in the z-direction using derived datatypes and subcommunicator.
        call MPI_Isend(Value_sub,1, ddtype_sendto_F  , comm_1d_x3%east_rank, 111, comm_1d_x3%mpi_comm, request(1), ierr)
        call MPI_Irecv(Value_sub,1, ddtype_recvfrom_B, comm_1d_x3%west_rank, 111, comm_1d_x3%mpi_comm, request(2), ierr)
        call MPI_Isend(Value_sub,1, ddtype_sendto_B  , comm_1d_x3%west_rank, 222, comm_1d_x3%mpi_comm, request(3), ierr)
        call MPI_Irecv(Value_sub,1, ddtype_recvfrom_F, comm_1d_x3%east_rank, 222, comm_1d_x3%mpi_comm, request(4), ierr)
        call MPI_Waitall(4, request, MPI_STATUSES_IGNORE, ierr)
        
        
    end subroutine mpi_subdomain_ghostcell_update

    subroutine subdomain_para_range(n1, n2, nprocs, myrank, ista, iend)

        implicit none

        integer, intent(in)     :: n1, n2, nprocs, myrank
        integer, intent(out)    :: ista, iend
        integer :: iwork1, iwork2

        iwork1 = int((n2 - n1 + 1) / nprocs)
        iwork2 = mod(n2 - n1 + 1, nprocs)
        ista = myrank * iwork1 + n1 + min(myrank, iwork2)
        iend = ista + iwork1 - 1
        if (iwork2 > myrank) iend = iend + 1

    end subroutine subdomain_para_range

end module mpi_subdomain
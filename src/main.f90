program main
    ! use debug
    use mpi
    use mpi_topology
    use global
    use mpi_subdomain !<- depend 'global' and 'mpi_topology'
    
    use mpi_thermal
    use mpi_momentum
    use mpi_pressure

    use mpi_Post !<- depend 'global' and 'mpi_topology' and 'mpi_subdomain'

    implicit none

    integer :: i,j,k,TimeStep
    integer :: nprocs, myrank, ierr

    double precision :: maxDivU,timer_a,timer_b
    
    call MPI_INIT(ierr)
    call MPI_COMM_SIZE( MPI_COMM_WORLD, nprocs, ierr)
    call MPI_COMM_RANK( MPI_COMM_WORLD, myrank, ierr)
    if(myrank==0) write(*,*) '[Main] The simulation starts!'

    if(myrank==0) call system('mkdir ./data')
    if(myrank==0) call system('mkdir ./data/1_continu')
    if(myrank==0) call system('mkdir ./data/2_instanfield')
    call MPI_Barrier(MPI_COMM_WORLD,ierr)
    
    call global_inputpara(np_dim) !: np_dim(:) <- 'mpi_topology'

    period(0)=.true.; period(1)=.false.; period(2)=.true. !: period(:) <- 'mpi_topology'
    call mpi_topology_make()
    call mpi_subdomain_make()
    call mpi_subdomain_mesh()

    call global_mpi_indices(n1sub,n2sub,n3sub) !<- depend 'mpi_topology'

    call mpi_subdomain_DDT_ghostcell()

    call mpi_thermal_allocation(n1sub,n2sub,n3sub)
    call mpi_thermal_initial(x1_sub, x2_sub, x3_sub)
    call mpi_thermal_boundary(x1_sub, x2_sub, x3_sub, comm_1d_x1, comm_1d_x2, comm_1d_x3)
    
    call mpi_momentum_allocation(n1sub,n2sub,n3sub)
    call mpi_momentum_initial(x1_sub, x2_sub, x3_sub)
    call mpi_momentum_boundary(comm_1d_x1, comm_1d_x2, comm_1d_x3)

    call mpi_pressure_allocation(n1sub,n2sub,n3sub, comm_1d_x1, comm_1d_x2, comm_1d_x3)
    call mpi_pressure_initial(x1_sub, x2_sub, x3_sub)
    call mpi_subdomain_DDT_FFT( n1msub,n3msub,n3mssub,wv1,wv1sub,wv1ssub,wv1ssubA,wv1ssubB    &
                              , ddtype_C2I_dble,ddtype_C2I_cplx,ddtype_C2K_cplx &
                              , ddtype_I2C_dble,ddtype_I2C_cplx,ddtype_K2C_cplx )

    CFL = MaxCFL
    time = tStart 
    dt = dtStart

    call mpi_subdomain_ghostcell_update(T)
    call mpi_subdomain_ghostcell_update(U)
    call mpi_subdomain_ghostcell_update(V)
    call mpi_subdomain_ghostcell_update(W)
    call mpi_subdomain_ghostcell_update(P)

    call MPI_Barrier(MPI_COMM_WORLD,ierr)

    call mpi_Post_allocation()

    if ( ConinueQ==1 ) then
        call mpi_Post_FileIn_Continue(myrank,tStart,dt,U,V,W,P,T)
        call mpi_Post_CFL(U,V,W,CFL,dt)
        time = tStart
    end if

            ! call debug_infor(n1sub,n2sub,n3sub,period,ista,iend,jsta,jend,ksta,kend    &
            !                 ,x1_sub, x2_sub, x3_sub                                    &
            !                 ,comm_1d_x1%nprocs,comm_1d_x2%nprocs,comm_1d_x3%nprocs     &
            !                 ,comm_1d_x1%myrank,comm_1d_x2%myrank,comm_1d_x3%myrank,myrank)

            ! Dprint_range(1:2,1)=(/0,comm_1d_x1%nprocs-1/)
            ! Dprint_range(1:2,2)=(/0,comm_1d_x2%nprocs-1/)
            ! Dprint_range(1:2,3)=(/0,comm_1d_x3%nprocs-1/)
            
    !================= All Setup Finish ===============]]]

    do TimeStep = 1, Timestepmax
        timer_a=MPI_WTIME()
        
        call mpi_thermal_coeffi() ! alloc :: KPP, dKPP, invRhoCp, dinvRhoCp
        call mpi_thermal_solver(U,V,W,dx1_sub,dx2_sub,dx3_sub,dmx1_sub,dmx2_sub,dmx3_sub, comm_1d_x1, comm_1d_x2, comm_1d_x3)
        call mpi_thermal_coeffi_clean() ! dealloc :: KPP, dKPP, invRhoCp, dinvRhoCp
        call mpi_subdomain_ghostcell_update(T)
        
        call mpi_momentum_coeffi(T) ! alloc :: Mu,invRho
        
        call mpi_momentum_solvedU(T,dx1_sub,dx2_sub,dx3_sub,dmx1_sub,dmx2_sub,dmx3_sub, comm_1d_x1, comm_1d_x2, comm_1d_x3) ! alloc :: dU
        call mpi_subdomain_ghostcell_update(dU)
        
        call mpi_momentum_solvedV(T,dx1_sub,dx2_sub,dx3_sub,dmx1_sub,dmx2_sub,dmx3_sub, comm_1d_x1, comm_1d_x2, comm_1d_x3) ! alloc :: dV
        call mpi_subdomain_ghostcell_update(dV)
        
        call mpi_momentum_solvedW(T,dx1_sub,dx2_sub,dx3_sub,dmx1_sub,dmx2_sub,dmx3_sub, comm_1d_x1, comm_1d_x2, comm_1d_x3) ! alloc :: dW
        call mpi_subdomain_ghostcell_update(dW)

        call mpi_momentum_blockLdV(T,dx1_sub,dx2_sub,dx3_sub,dmx1_sub,dmx2_sub,dmx3_sub)
        call mpi_subdomain_ghostcell_update(dV)        
        call mpi_momentum_blockLdU(T,dx1_sub,dx2_sub,dx3_sub,dmx1_sub,dmx2_sub,dmx3_sub)

        call mpi_momentum_pseudoupdateUVW() ! dealloc :: dU,dV,dW
        call mpi_subdomain_ghostcell_update(U)
        call mpi_subdomain_ghostcell_update(V)
        call mpi_subdomain_ghostcell_update(W)

        call mpi_pressure_RHS(invRho,U,V,W,VBCup_sub,VBCbt_sub,dx1_sub,dx2_sub,dx3_sub,dmx1_sub,dmx2_sub,dmx3_sub) ! alloc :: PRHS
        call mpi_pressure_Poisson_FFT(dx2_sub,dmx2_sub, comm_1d_x1, comm_1d_x2, comm_1d_x3) ! dealloc :: PRHS, alloc :: dP
        call mpi_subdomain_ghostcell_update(dP)
        call mpi_pressure_Projection(invRho,U,V,W,P,dx1_sub,dx2_sub,dx3_sub,dmx1_sub,dmx2_sub,dmx3_sub,comm_1d_x2%myrank,comm_1d_x2%nprocs-1) ! dealloc :: dP
        call mpi_subdomain_ghostcell_update(U)
        call mpi_subdomain_ghostcell_update(V)
        call mpi_subdomain_ghostcell_update(W)
        call mpi_subdomain_ghostcell_update(P)
        call mpi_subdomain_ghostcell_update(dPhat)
        
        call mpi_momentum_coeffi_clean() ! dealloc :: Mu,invRho
        
        call mpi_Post_Div(U,V,W,maxDivU)
        
        if(maxDivU>=1.d-3) then
            call mpi_Post_FileOut_InstanField(myrank,TimeStep,time,U,V,W,P,T)
            if(myrank==0) write(*,*) 'BLOW UP',TimeStep
            exit
        endif

        call mpi_Post_FileOut_InstanField(myrank,TimeStep,time,U,V,W,P,T)
        
        timer_b=MPI_WTIME()
        call mpi_Post_MonitorOut(myrank,TimeStep,time,dt,CFL,maxDivU,timer_b-timer_a)
        call mpi_Post_CFL(U,V,W,CFL,dt)

        time=time+dt

    enddo

    call mpi_Post_FileOut_Continue(myrank,TimeStep,time,dt,U,V,W,P,T)

    ! call debug_clean()
    call mpi_pressure_clean()
    call mpi_momentum_clean()    
    call mpi_thermal_clean()

    call global_mpi_indices_clean()

    call mpi_subdomain_clean()
    call mpi_topology_clean()
    call MPI_FINALIZE(ierr)
    if(myrank==0) write(*,*) '[Main] The main simulation complete! '
end
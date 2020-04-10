module global
    implicit none
    double precision, parameter :: PI = acos(-1.d0)

    ! PHYSICAL PARAMETERS
    double precision :: Pr, Ra, Gr
    double precision :: Cmu, Cmt, Ct, Cmp

    ! Mesh
    integer :: n1,n2,n3
    integer :: n1m,n2m,n3m
    integer :: n1p,n2p,n3p

    integer :: UNIFORM1,UNIFORM2,UNIFORM3
    double precision :: GAMMA1, GAMMA2, GAMMA3

    integer :: Timestepmax
    double precision :: time,dt
    double precision :: CFL,MaxCFL, dtStart, tStart, dtMax

    double precision :: H,Aspect1,Aspect3
    double precision :: L1,L2,L3
    double precision :: x1_start,x2_start,x3_start
    double precision :: x1_end,x2_end,x3_end

    double precision :: Uup,Vup,Wup,Tup
    double precision :: Ubt,Vbt,Wbt,Tbt

    integer, public :: i_indexC,i_indexS,j_indexC,j_indexS,k_indexC,k_indexS
    integer, allocatable, dimension(:), public :: iC_BC,iS_BC,jC_BC,jS_BC,kC_BC,kS_BC
    
    double precision :: DeltaT
    double precision :: a10, a11, a12, a13, a14, a15, a12pera11
    double precision :: b10, b11, b12, b13, b14, b15
    double precision :: c10, c11, c12, c13, c14, c15
    double precision :: d10, d11, d12, d13, d14, d15
    double precision :: KPP0, Rho0, Cp0, Mu0

    character(len=17) dir_continu
    character(len=21) dir_instantfield
    character(len=512) filedirectory
    integer :: print_start_step,print_interval_step,print_j_index_wall,print_j_index_bulk,print_k_index
    integer :: ConinueQ,steamL
contains
    
    subroutine global_inputpara(np_dim)
        implicit none
        integer, intent(out) :: np_dim(0:2)
        double precision :: tttmp(1:3)
        character*512 :: temp_char
        integer :: i

            !Using file input
            open(unit = 1, file = "PARA_INPUT.dat")
            read(1,*) temp_char,ConinueQ,filedirectory
            read(1,*) temp_char                    !!--
            read(1,*) temp_char,np_dim(0)          !Xprocs
            read(1,*) temp_char,np_dim(1)          !Yprocs
            read(1,*) temp_char,np_dim(2)          !Zprocs
            read(1,*) temp_char                    !!--
            read(1,*) temp_char,n1m                !N1m
            read(1,*) temp_char,n2m                !N2m
            read(1,*) temp_char,n3m                !N3m
            read(1,*) temp_char,UNIFORM1,temp_char,GAMMA1    !gamma
            read(1,*) temp_char,UNIFORM2,temp_char,GAMMA2    !gamma
            read(1,*) temp_char,UNIFORM3,temp_char,GAMMA3    !gamma
            read(1,*) temp_char,Aspect1            !AspectRatio1
            read(1,*) temp_char,H                  !AspectRatio2
            read(1,*) temp_char,Aspect3            !AspectRatio3
            read(1,*) temp_char                    !!--
            read(1,*) temp_char,Ra                 !Ra
            read(1,*) temp_char,Pr                 !Pr
            read(1,*) temp_char,DeltaT             !DeltaT
            read(1,*) temp_char                    !!--
            read(1,*) temp_char,MaxCFL             !MaxCFL
            read(1,*) temp_char                    !!--
            read(1,*) temp_char,Timestepmax        !TotalTimeStep
            read(1,*) temp_char,print_start_step   !print_start_step
            read(1,*) temp_char,print_interval_step!PrintSteplength
            close(1)
                
            steamL=0
            do i=1,512
                if(filedirectory(i:i)=='%') EXIT
                steamL=steamL+1
            enddo

                dtStart = 5.0D-2; tStart = 0.d0; 
                    
                Uup = 0.d0;Ubt = 0.d0
                Vup = 0.d0;Vbt = 0.d0
                Wup = 0.d0;Wbt = 0.d0
                Tup = -0.5d0;Tbt =0.5d0
                
                ! Coefficients for Water
                ! a10 = 0.9922d+3; a11 = -3.736d-4;  a12 = -3.98d-6;  a13 = 0.;         a14 = 0.;  a15 = 0.
                ! b10 = 4.1690d+3; b11 =  0.084d-4;  b12 =  4.60d-6;  b13 = 0.;         b14 = 0.;  b15 = 0.
                ! c10 = 0.6297;    c11 = 21.99d-4;   c12 = -17.8d-6;  c13 = 0.;         c14 = 0.;  c15 = 0. 
                ! d10 = 0.6690d-6; d11 = -175.9d-4;  d12 = 295.8d-6;  d13 = -460.d-8;   d14 = 0.;  d15 = 0.
                a10 = 1.d0; a11 =0.d0; a12 = 0.d0;  a13 = 0.d0; a14 = 0.;  a15 = 0.
                b10 = 1.d0; b11 =0.d0; b12 = 0.d0;  b13 = 0.d0; b14 = 0.;  b15 = 0.
                c10 = 1.d0; c11 =0.d0; c12 = 0.d0;  c13 = 0.d0; c14 = 0.;  c15 = 0. 
                d10 = 1.d0; d11 =0.d0; d12 = 0.d0;  d13 = 0.d0; d14 = 0.;  d15 = 0.
                
                if(a12<1.0D-14) then
                    a12pera11=0.d0
                else
                    a12pera11=a12/a11
                endif
                dir_continu='./data/1_continu/'
                dir_instantfield='./data/2_instanfield/'
                
                ! Coefficients for Glycerol
                ! a10 = 1.2477d+3; a11 = -4.789d-4;  a12 = -0.3795d-6;  a13 = 0.;         a14 = 0.;  a15 = 0.
                ! b10 = 2.5108d+3; b11 = 22.511d-4;  b12 =  0.       ;  b13 = 0.;         b14 = 0.;  b15 = 0.
                ! c10 = 2.9351d-3; c11 = 3.863d-4;   c12 =  0.       ;  c13 = 0.;         c14 = 0.;  c15 = 0. 
                ! d10 = 238.71d-6; d11 =-702.83d-4; d12 = 2393.1d-6;    d13 = -6923.0d-8; d14 = 33131.3d-10; d15 = -71517.5d-12

            !----
            
            Rho0  = a10
            Cp0   = b10
            KPP0  = c10
            Mu0   = a10*d10            

            Cmu = sqrt(Pr/Ra); Cmt = 1.d0; Ct  = 1.d0/sqrt(Ra*Pr)
            Gr  = Ra/Pr
            Cmp = 1.d0

            n1=n1m+1;n1p=n1+1;
            n2=n2m+1;n2p=n2+1;
            n3=n3m+1;n3p=n3+1;
            
            L1=H*Aspect1
            L2=H
            L3=H*Aspect3

            x1_end=x1_start+L1
            x2_end=x2_start+L2
            x3_end=x3_start+L3
            
            tttmp(1)=L1/dble(n1-1)
            tttmp(2)=L2/dble(n2-1)
            tttmp(3)=L3/dble(n3-1)

            dtMax=minval(tttmp)*100
            ! dtMax=5.d-2
            
            !--------

    end subroutine global_inputpara

    subroutine global_mpi_indices(n1sub,n2sub,n3sub)
        use mpi_topology
        implicit none
        integer :: n1sub,n2sub,n3sub

        allocate(iC_BC(0:n1sub),iS_BC(0:n1sub))
        allocate(jC_BC(0:n2sub),jS_BC(0:n2sub))
        allocate(kC_BC(0:n3sub),kS_BC(0:n3sub))

        i_indexC=1; i_indexS=1
        j_indexC=1; j_indexS=1
        k_indexC=1; k_indexS=1
        
        if(period(0)==.false. .and. comm_1d_x1%myrank==0) i_indexS=2
        if(period(1)==.false. .and. comm_1d_x2%myrank==0) j_indexS=2
        if(period(2)==.false. .and. comm_1d_x3%myrank==0) k_indexS=2

        iC_BC(0:n1sub)=1;iS_BC(0:n1sub)=1;
        jC_BC(0:n2sub)=1;jS_BC(0:n2sub)=1;
        kC_BC(0:n3sub)=1;kS_BC(0:n3sub)=1;

        if(period(0)==.false. .and. comm_1d_x1%myrank==0) iC_BC(0:0)=0!iC_BC(0:1)=0
        if(period(0)==.false. .and. comm_1d_x1%myrank==0) iS_BC(0:1)=0!iS_BC(0:2)=0
        if(period(1)==.false. .and. comm_1d_x2%myrank==0) jC_BC(0:0)=0!jC_BC(0:1)=0
        if(period(1)==.false. .and. comm_1d_x2%myrank==0) jS_BC(0:1)=0!jS_BC(0:2)=0
        if(period(2)==.false. .and. comm_1d_x3%myrank==0) kC_BC(0:0)=0!kC_BC(0:1)=0
        if(period(2)==.false. .and. comm_1d_x3%myrank==0) kS_BC(0:1)=0!kS_BC(0:2)=0

        if(period(0)==.false. .and. comm_1d_x1%myrank==comm_1d_x1%nprocs-1) iC_BC(n1sub:n1sub)=0!iC_BC(n1sub-1:n1sub)=0
        if(period(0)==.false. .and. comm_1d_x1%myrank==comm_1d_x1%nprocs-1) iS_BC(n1sub:n1sub)=0!iS_BC(n1sub-1:n1sub)=0
        if(period(1)==.false. .and. comm_1d_x2%myrank==comm_1d_x2%nprocs-1) jC_BC(n2sub:n2sub)=0!jC_BC(n2sub-1:n2sub)=0
        if(period(1)==.false. .and. comm_1d_x2%myrank==comm_1d_x2%nprocs-1) jS_BC(n2sub:n2sub)=0!jS_BC(n2sub-1:n2sub)=0
        if(period(2)==.false. .and. comm_1d_x3%myrank==comm_1d_x3%nprocs-1) kC_BC(n3sub:n3sub)=0!kC_BC(n3sub-1:n3sub)=0
        if(period(2)==.false. .and. comm_1d_x3%myrank==comm_1d_x3%nprocs-1) kS_BC(n3sub:n3sub)=0!kS_BC(n3sub-1:n3sub)=0
        
    end subroutine global_mpi_indices

    subroutine global_mpi_indices_clean()
        implicit none

        deallocate(iC_BC,iS_BC)
        deallocate(jC_BC,jS_BC)
        deallocate(kC_BC,kS_BC)
        
    end subroutine global_mpi_indices_clean
end module global
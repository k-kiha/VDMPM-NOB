module mpi_topology
    use mpi
    implicit none

    integer, public :: mpi_world_cart
    integer, public :: np_dim(0:2)
    logical, public :: period(0:2)

    !> @brief   Type variable for the information of 1D communicator
    type, public :: cart_comm_1d
        integer :: myrank
        integer :: nprocs
        integer :: west_rank
        integer :: east_rank
        integer :: mpi_comm
    end type cart_comm_1d

    type(cart_comm_1d), public :: comm_1d_x1
    type(cart_comm_1d), public :: comm_1d_x2
    type(cart_comm_1d), public :: comm_1d_x3

    private

    public  :: mpi_topology_make
    public  :: mpi_topology_clean

    contains

    !>
    !> @brief       Destroy the communicator for cartesian topology.
    !>
    subroutine mpi_topology_clean()

        implicit none
        integer :: ierr

        call MPI_Comm_free(mpi_world_cart, ierr)

    end subroutine mpi_topology_clean

    !>
    !> @brief       Create the cartesian topology for the MPI processes and subcommunicators.
    !>
    subroutine mpi_topology_make()
        implicit none
        logical :: remain(0:2)
        integer :: ierr

        ! Create the cartesian topology.
        call MPI_Cart_create( MPI_COMM_WORLD    &!  input  | integer      | Input communicator (handle).
                            , 3                 &!  input  | integer      | Number of dimensions of Cartesian grid (integer).
                            , np_dim            &!  input  | integer(1:3) | Integer array of size ndims specifying the number of processes in each dimension.
                            , period            &!  input  | logical(1:3) | Logical array of size ndims specifying whether the grid is periodic (true=1) or not (false=0) in each dimension.
                            , .false.           &!  input  | logical      | Ranking may be reordered (true=1) or not (false=0) (logical).
                            , mpi_world_cart    &! *output | integer      | Communicator with new Cartesian topology (handle).
                            , ierr              &!  output | integer      | Fortran only: Error status
                            )

        ! Create subcommunicators and assign two neighboring processes in the x-direction.
        remain(0) = .true.
        remain(1) = .false.
        remain(2) = .false.
        call MPI_Cart_sub( mpi_world_cart, remain, comm_1d_x1%mpi_comm, ierr)
        call MPI_Comm_rank(comm_1d_x1%mpi_comm, comm_1d_x1%myrank, ierr)
        call MPI_Comm_size(comm_1d_x1%mpi_comm, comm_1d_x1%nprocs, ierr)
        call MPI_Cart_shift(comm_1d_x1%mpi_comm, 0, 1, comm_1d_x1%west_rank, comm_1d_x1%east_rank, ierr)

        ! Create subcommunicators and assign two neighboring processes in the y-direction
        remain(0) = .false.
        remain(1) = .true.
        remain(2) = .false.
        call MPI_Cart_sub( mpi_world_cart, remain, comm_1d_x2%mpi_comm, ierr)
        call MPI_Comm_rank(comm_1d_x2%mpi_comm, comm_1d_x2%myrank, ierr)
        call MPI_Comm_size(comm_1d_x2%mpi_comm, comm_1d_x2%nprocs, ierr)
        call MPI_Cart_shift(comm_1d_x2%mpi_comm, 0, 1, comm_1d_x2%west_rank, comm_1d_x2%east_rank, ierr)

        ! Create subcommunicators and assign two neighboring processes in the z-direction
        remain(0) = .false.
        remain(1) = .false.
        remain(2) = .true.
        call MPI_Cart_sub( mpi_world_cart, remain, comm_1d_x3%mpi_comm, ierr)
        call MPI_Comm_rank(comm_1d_x3%mpi_comm, comm_1d_x3%myrank, ierr)
        call MPI_Comm_size(comm_1d_x3%mpi_comm, comm_1d_x3%nprocs, ierr)
        call MPI_Cart_shift(comm_1d_x3%mpi_comm, 0, 1, comm_1d_x3%west_rank, comm_1d_x3%east_rank, ierr)

    end subroutine mpi_topology_make

end module mpi_topology

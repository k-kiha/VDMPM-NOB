module mpi_Post
    ! use debug
    use mpi_topology
    use global
    use mpi_subdomain
    implicit none
    logical :: FirstPrint1=.true.,FirstPrint2=.true.,FirstPrint3=.true.
    character(len=17) :: xyzrank

contains
    subroutine mpi_Post_allocation()
        implicit none
        
        integer :: r1_0,r2_0,r3_0
        integer :: r1_1,r2_1,r3_1
        integer :: r1_2,r2_2,r3_2
        integer :: r1_3,r2_3,r3_3
        integer :: r1_4,r2_4,r3_4
        
        r1_4=int((comm_1d_x1%myrank)/10000)
        r1_3=int((comm_1d_x1%myrank-r1_4*10000)/1000 )
        r1_2=int((comm_1d_x1%myrank-r1_4*10000-r1_3*1000)/100 )
        r1_1=int((comm_1d_x1%myrank-r1_4*10000-r1_3*1000-r1_2*100)/10 )
        r1_0=int((comm_1d_x1%myrank-r1_4*10000-r1_3*1000-r1_2*100-r1_1*10)/1)

        r2_4=int((comm_1d_x2%myrank)/10000)
        r2_3=int((comm_1d_x2%myrank-r2_4*10000)/1000 )
        r2_2=int((comm_1d_x2%myrank-r2_4*10000-r2_3*1000)/100 )
        r2_1=int((comm_1d_x2%myrank-r2_4*10000-r2_3*1000-r2_2*100)/10 )
        r2_0=int((comm_1d_x2%myrank-r2_4*10000-r2_3*1000-r2_2*100-r2_1*10)/1)
        
        r3_4=int((comm_1d_x3%myrank)/10000)
        r3_3=int((comm_1d_x3%myrank-r3_4*10000)/1000 )
        r3_2=int((comm_1d_x3%myrank-r3_4*10000-r3_3*1000)/100 )
        r3_1=int((comm_1d_x3%myrank-r3_4*10000-r3_3*1000-r3_2*100)/10 )
        r3_0=int((comm_1d_x3%myrank-r3_4*10000-r3_3*1000-r3_2*100-r3_1*10)/1)

        
        write(xyzrank, '(5I1,1A1,5I1,1A1,5I1)' ) r1_4,r1_3,r1_2,r1_1,r1_0,'_'   &
                                                 ,r2_4,r2_3,r2_2,r2_1,r2_0,'_'   &
                                                 ,r3_4,r3_3,r3_2,r3_1,r3_0
    end subroutine mpi_Post_allocation

    subroutine mpi_Post_Div(U,V,W,maxDivU)
        use MPI
        implicit none
        double precision, dimension(0:n1sub,0:n2sub,0:n3sub) :: U,V,W
        double precision :: maxDivU

        integer :: i,j,k
        integer :: im,jm,km
        integer :: ip,jp,kp
        integer :: Sc,Ec
        integer :: Sm,Em
        integer :: Sp,Ep
        double precision, allocatable, dimension(:,:,:) :: DivU
        double precision :: maxI,maxIJ,maxIJK
        double precision, allocatable, dimension(:) :: findmaxk
        double precision, allocatable, dimension(:,:) :: findmaxIJ
        integer :: ierr
    
        allocate(findmaxIJ(1:n1sub-1,1:n2sub-1),findmaxk(1:n3sub-1))
        
        Sc=1  ;Ec=n1sub-1
        Sm=1-1;Em=n1sub-1-1
        Sp=1+1;Ep=n1sub-1+1
        
        do k = 1, n3sub-1
            kp = k+1
            km = k-1
            do j = 1, n2sub-1
                jp = j + 1
                jm = j - 1
            
                findmaxIJ(Sc:Ec,j) = ( U(Sp:Ep,j ,k ) - U(Sc:Ec,j,k)  )/dx1_sub(Sc:Ec)    &
                                   + ( V(Sc:Ec,jp,k ) - V(Sc:Ec,j,k)  )/dx2_sub(j)        &
                                   + ( W(Sc:Ec,j ,kp) - W(Sc:Ec,j,k)  )/dx3_sub(k)

            end do
            findmaxk(k)=maxval(findmaxIJ)
        end do

        maxDivU=maxval(findmaxk)
        call MPI_ALLREDUCE(maxDivU, maxI  , 1, MPI_DOUBLE, MPI_MAX, comm_1d_x1%mpi_comm, ierr)
        call MPI_ALLREDUCE(maxI   , maxIJ , 1, MPI_DOUBLE, MPI_MAX, comm_1d_x2%mpi_comm, ierr)
        call MPI_ALLREDUCE(maxIJ  , maxIJK, 1, MPI_DOUBLE, MPI_MAX, comm_1d_x3%mpi_comm, ierr)
        maxDivU=maxIJK

        deallocate(findmaxIJ,findmaxk)

    end subroutine mpi_Post_Div

    subroutine mpi_Post_CFL(U,V,W,newCFL,newdt)
        use MPI
        implicit none
        double precision, dimension(0:n1sub,0:n2sub,0:n3sub) :: U,V,W
        double precision :: newCFL,newdt

        integer :: i,j,k
        integer :: im,jm,km
        integer :: ip,jp,kp
        integer :: Sc,Ec
        integer :: Sm,Em
        integer :: Sp,Ep
        double precision, allocatable, dimension(:,:,:) :: DivU
        double precision :: maxI,maxIJ,maxIJK,maxUovX
        double precision, allocatable, dimension(:) :: findmaxk
        double precision, allocatable, dimension(:,:) :: findmaxIJ
        integer :: ierr
    
        allocate(findmaxIJ(1:n1sub-1,1:n2sub-1),findmaxk(1:n3sub-1))
        
        Sc=1  ;Ec=n1sub-1
        Sm=1-1;Em=n1sub-1-1
        Sp=1+1;Ep=n1sub-1+1
        
        do k = 1, n3sub-1
            kp = k+1
            km = k-1
            do j = 1, n2sub-1
                jp = j + 1
                jm = j - 1
            
                findmaxIJ(Sc:Ec,j) = (ABS( U(Sp:Ep,j ,k ) + U(Sc:Ec,j,k)  )/2.d0)/dx1_sub(Sc:Ec)    &
                                   + (ABS( V(Sc:Ec,jp,k ) + V(Sc:Ec,j,k)  )/2.d0)/dx2_sub(j)        &
                                   + (ABS( W(Sc:Ec,j ,kp) + W(Sc:Ec,j,k)  )/2.d0)/dx3_sub(k)

            end do
            findmaxk(k)=maxval(findmaxIJ)
        end do

        maxUovX=maxval(findmaxk)
        call MPI_ALLREDUCE(maxUovX, maxI  , 1, MPI_DOUBLE, MPI_MAX, comm_1d_x1%mpi_comm, ierr)
        call MPI_ALLREDUCE(maxI   , maxIJ , 1, MPI_DOUBLE, MPI_MAX, comm_1d_x2%mpi_comm, ierr)
        call MPI_ALLREDUCE(maxIJ  , maxIJK, 1, MPI_DOUBLE, MPI_MAX, comm_1d_x3%mpi_comm, ierr)
        maxUovX=maxIJK
        newCFL=newdt*maxUovX
        newdt=MaxCFL/maxUovX
        if(newdt>dtMax) newdt=dtMax
        newCFL=newdt*maxUovX

        deallocate(findmaxIJ,findmaxk)

    end subroutine mpi_Post_CFL

    subroutine mpi_Post_MonitorOut(myrank,outTimeStep,outtime,outdt,outCFL,outmaxDivU,outtimer)
        implicit none
        integer :: outTimeStep,myrank
        double precision :: outtime,outdt,outCFL,outmaxDivU,outtimer

        if(mod(outTimeStep,20)==1) then
            if(myrank==0) write(*,*)
            if(myrank==0) write(*,'(6A15)') 'Timestep', 'Time', 'dt', 'CFL', 'max_DivU', 'WTime/step'
        endif
        if(myrank==0) write(*,'(1I15,5E15.5)') outTimeStep,outtime,outdt,outCFL,outmaxDivU,outtimer
    end subroutine mpi_Post_MonitorOut

    subroutine mpi_Post_FileOut_InstanField(myrank,tstepin,timein,Uin,Vin,Win,Pin,Tin)
        use mpi
        implicit none
        integer :: myrank,tstepin
        double precision :: timein
        double precision, dimension(0:n1sub, 0:n2sub, 0:n3sub) :: Uin,Vin,Win,Pin,Tin
        character(len=22) :: filename_instantfieldXY
        character(len=27) :: filename_instantfieldXZ_wall,filename_instantfieldXZ_bulk

        integer :: i,j,k,onebyone,ierr
        integer :: print_rankj_wall,print_rankj_bulk,print_rankk

        filename_instantfieldXY='Output_instantfield_XY'
        filename_instantfieldXZ_wall='Output_instantfield_XZ_wall'
        filename_instantfieldXZ_bulk='Output_instantfield_XZ_bulk'

        if(tstepin==0.or.(tstepin>=print_start_step.and.mod(tstepin-print_start_step,print_interval_step)==0)) then

            if(comm_1d_x3%myrank==comm_1d_x3%nprocs-1) then
                open(unit=myrank,file=dir_instantfield//filename_instantfieldXY//xyzrank//'.plt', position='append')
                    if(FirstPrint1) then
                        write(myrank,*) 'VARIABLES="X","Y","Z","U","V","W","P","T"'  !--
                        write(myrank,*) 'zone t="',tstepin,'"','i=',n1sub+1,'j=',n2sub+1,'k=',1
                        write(myrank,*) 'STRANDID=', tstepin
                        write(myrank,*) 'SOLUTIONTIME=', timein
                        k= (n3sub-1)/2
                        do j=0,n2sub
                        do i=0,n1sub
                            write(myrank,'(3D20.10,5D30.20)') x1_sub(i),x2_sub(j),x3_sub(k)    &
                                                        &, Uin(i,j,k),Vin(i,j,k),Win(i,j,k),Pin(i,j,k),Tin(i,j,k)
                        enddo
                        enddo

                        FirstPrint1=.false.
                    else
                        write(myrank,*) 'zone t="',tstepin,'"','i=',n1sub+1,'j=',n2sub+1,'k=',1
                        write(myrank,*) 'VARSHARELIST= ([1-3]=1)'
                        write(myrank,*) 'STRANDID=', tstepin
                        write(myrank,*) 'SOLUTIONTIME=', timein
                        k= (n3sub-1)/2
                        do j=0,n2sub
                        do i=0,n1sub
                            write(myrank,'(5D30.20)') Uin(i,j,k),Vin(i,j,k),Win(i,j,k),Pin(i,j,k),Tin(i,j,k)
                        enddo
                        enddo

                    endif
                close(myrank)
            endif

            if(comm_1d_x2%myrank==0) then
                open(unit=myrank,file=dir_instantfield//filename_instantfieldXZ_wall//xyzrank//'.plt', position='append')
                    if(FirstPrint2) then
                        write(myrank,*) 'VARIABLES="X","Y","Z","U","V","W","P","T"'  !--
                        write(myrank,*) 'zone t="',tstepin,'"','i=',n1sub+1,'j=',1,'k=',n3sub+1
                        write(myrank,*) 'STRANDID=', tstepin
                        write(myrank,*) 'SOLUTIONTIME=', timein
                        j=3
                        do k=0,n3sub
                        do i=0,n1sub
                            write(myrank,'(3D20.10,5D30.20)') x1_sub(i),x2_sub(j),x3_sub(k)    &
                                                        &, Uin(i,j,k),Vin(i,j,k),Win(i,j,k),Pin(i,j,k),Tin(i,j,k)
                        enddo
                        enddo

                        FirstPrint2=.false.
                    else
                        write(myrank,*) 'zone t="',tstepin,'"','i=',n1sub+1,'j=',1,'k=',n3sub+1
                        write(myrank,*) 'VARSHARELIST= ([1-3]=1)'
                        write(myrank,*) 'STRANDID=', tstepin
                        write(myrank,*) 'SOLUTIONTIME=', timein
                        j=3
                        do k=0,n3sub
                        do i=0,n1sub
                            write(myrank,'(5D30.20)') Uin(i,j,k),Vin(i,j,k),Win(i,j,k),Pin(i,j,k),Tin(i,j,k)
                        enddo
                        enddo

                    endif
                close(myrank)
            endif
        endif

    end subroutine mpi_Post_FileOut_InstanField

    subroutine mpi_Post_FileIn_Continue(myrank,timein,dt,Uin,Vin,Win,Pin,Tin)
        use mpi
        implicit none
        integer :: myrank
        double precision :: timein,dt
        double precision, dimension(0:n1sub, 0:n2sub, 0:n3sub) :: Uin,Vin,Win,Pin,Tin

        integer :: i,j,k,onebyone,ierr
        
        do onebyone=0,comm_1d_x2%nprocs-1
            if(onebyone==comm_1d_x2%myrank) then
                OPEN(myrank, FILE=filedirectory(1:steamL)//'binarary_time_'//xyzrank//'.data'   &
                        &, FORM='UNFORMATTED', ACCESS='DIRECT', RECL=sizeof(timein)*2        )
                    read(myrank, REC=1) timein,dt
                CLOSE(myrank)
                OPEN(myrank, FILE=filedirectory(1:steamL)//'binarary_U_'//xyzrank//'.data'   &
                        &, FORM='UNFORMATTED', ACCESS='DIRECT', RECL=sizeof(Uin)        )
                    read(myrank, REC=1) Uin
                CLOSE(myrank)
                OPEN(myrank, FILE=filedirectory(1:steamL)//'binarary_V_'//xyzrank//'.data'   &
                        &, FORM='UNFORMATTED', ACCESS='DIRECT', RECL=sizeof(Vin)        )
                    read(myrank, REC=1) Vin
                CLOSE(myrank)
                OPEN(myrank, FILE=filedirectory(1:steamL)//'binarary_W_'//xyzrank//'.data'   &
                        &, FORM='UNFORMATTED', ACCESS='DIRECT', RECL=sizeof(Win)        )
                    read(myrank, REC=1) Win
                CLOSE(myrank)
                OPEN(myrank, FILE=filedirectory(1:steamL)//'binarary_P_'//xyzrank//'.data'   &
                        &, FORM='UNFORMATTED', ACCESS='DIRECT', RECL=sizeof(Pin)        )
                    read(myrank, REC=1) Pin
                CLOSE(myrank)
                OPEN(myrank, FILE=filedirectory(1:steamL)//'binarary_THETA_'//xyzrank//'.data'   &
                        &, FORM='UNFORMATTED', ACCESS='DIRECT', RECL=sizeof(Tin)        )
                    read(myrank, REC=1) Tin
                CLOSE(myrank)
            endif
            call MPI_Barrier(MPI_COMM_WORLD,ierr)
        enddo

    end subroutine mpi_Post_FileIn_Continue

    subroutine mpi_Post_FileOut_Continue(myrank,tstepin,timein,dt,Uin,Vin,Win,Pin,Tin)
        use mpi
        implicit none
        integer :: myrank,tstepin
        double precision :: timein,dt
        double precision, dimension(0:n1sub, 0:n2sub, 0:n3sub) :: Uin,Vin,Win,Pin,Tin

        integer :: i,j,k,onebyone,ierr
        
        do onebyone=0,comm_1d_x2%nprocs-1
            if(onebyone==comm_1d_x2%myrank) then
                OPEN(myrank, FILE=dir_continu//'binarary_time_'//xyzrank//'.data'    &
                        &, FORM='UNFORMATTED', ACCESS='DIRECT', RECL=sizeof(timein)*2              )
                    WRITE(myrank, REC=1) timein,dt
                CLOSE(myrank)
                OPEN(myrank, FILE=dir_continu//'binarary_U_'//xyzrank//'.data'   &
                        &, FORM='UNFORMATTED', ACCESS='DIRECT', RECL=sizeof(Uin)        )
                    WRITE(myrank, REC=1) Uin
                CLOSE(myrank)
                OPEN(myrank, FILE=dir_continu//'binarary_V_'//xyzrank//'.data'   &
                        &, FORM='UNFORMATTED', ACCESS='DIRECT', RECL=sizeof(Vin)        )
                    WRITE(myrank, REC=1) Vin
                CLOSE(myrank)
                OPEN(myrank, FILE=dir_continu//'binarary_W_'//xyzrank//'.data'   &
                        &, FORM='UNFORMATTED', ACCESS='DIRECT', RECL=sizeof(Win)        )
                    WRITE(myrank, REC=1) Win
                CLOSE(myrank)
                OPEN(myrank, FILE=dir_continu//'binarary_P_'//xyzrank//'.data'   &
                        &, FORM='UNFORMATTED', ACCESS='DIRECT', RECL=sizeof(Pin)        )
                    WRITE(myrank, REC=1) Pin
                CLOSE(myrank)
                OPEN(myrank, FILE=dir_continu//'binarary_THETA_'//xyzrank//'.data'   &
                        &, FORM='UNFORMATTED', ACCESS='DIRECT', RECL=sizeof(Tin)        )
                    WRITE(myrank, REC=1) Tin
                CLOSE(myrank)
            endif
            call MPI_Barrier(MPI_COMM_WORLD,ierr)
        enddo

    end subroutine mpi_Post_FileOut_Continue


end module mpi_Post
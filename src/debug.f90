module debug
    implicit none
    integer :: Dn1sub,Dn2sub,Dn3sub
    logical :: Dperiod(0:2)
    integer :: Dista,Diend,Djsta,Djend,Dksta,Dkend
    integer :: Dnprocs_1,Dnprocs_2,Dnprocs_3
    integer :: Dmyrank_1,Dmyrank_2,Dmyrank_3,Dmyrank
    integer :: Dprint_range(1:2,1:3)

    double precision, allocatable, dimension(:), public :: y1_sub, y2_sub, y3_sub
    double precision, allocatable, dimension(:,:,:), public :: debug_tmp1,debug_tmp2,debug_tmp3,debug_tmp4,debug_tmp5
    complex(8), allocatable, dimension(:,:,:), public :: debug_tmp_cplx1,debug_tmp_cplx2,debug_tmp_cplx3,debug_tmp_cplx4,debug_tmp_cplx5
    
    character(len=17) :: filerank

    double precision, dimension(0:2) :: timer_1Toeffi=0.d0
    double precision, dimension(0:2) :: timer_2Tsolver=0.d0,timer_2Tsolver_1RHSnAk=0.d0
    double precision, dimension(0:2) :: timer_2Tsolver_2TDMAk=0.d0
    double precision, dimension(0:2) :: timer_2Tsolver_3Trans1=0.d0
    double precision, dimension(0:2) :: timer_2Tsolver_4Ai=0.d0
    double precision, dimension(0:2) :: timer_2Tsolver_5TDMAi=0.d0
    double precision, dimension(0:2) :: timer_2Tsolver_6Trans2=0.d0
    double precision, dimension(0:2) :: timer_2Tsolver_7Aj=0.d0
    double precision, dimension(0:2) :: timer_2Tsolver_8TDMAj=0.d0
    double precision, dimension(0:2) :: timer_2Tsolver_9Update=0.d0

    double precision, dimension(0:2) :: timer_3Ucoeffi=0.d0
    double precision, dimension(0:2) :: timer_4dUsolver=0.d0,timer_4dUsolver_1RHSnAk=0.d0
    double precision, dimension(0:2) :: timer_4dUsolver_2TDMAk=0.d0
    double precision, dimension(0:2) :: timer_4dUsolver_3Trans1=0.d0
    double precision, dimension(0:2) :: timer_4dUsolver_4Ai=0.d0
    double precision, dimension(0:2) :: timer_4dUsolver_5TDMAi=0.d0
    double precision, dimension(0:2) :: timer_4dUsolver_6Trans2=0.d0
    double precision, dimension(0:2) :: timer_4dUsolver_7Aj=0.d0
    double precision, dimension(0:2) :: timer_4dUsolver_8TDMAj=0.d0
    double precision, dimension(0:2) :: timer_4dUsolver_9Update=0.d0

    double precision, dimension(0:2) :: timer_5dVsolver=0.d0,timer_5dVsolver_1RHSnAk=0.d0
    double precision, dimension(0:2) :: timer_5dVsolver_2TDMAk=0.d0
    double precision, dimension(0:2) :: timer_5dVsolver_3Trans1=0.d0
    double precision, dimension(0:2) :: timer_5dVsolver_4Ai=0.d0
    double precision, dimension(0:2) :: timer_5dVsolver_5TDMAi=0.d0
    double precision, dimension(0:2) :: timer_5dVsolver_6Trans2=0.d0
    double precision, dimension(0:2) :: timer_5dVsolver_7Aj=0.d0
    double precision, dimension(0:2) :: timer_5dVsolver_8TDMAj=0.d0
    double precision, dimension(0:2) :: timer_5dVsolver_9Update=0.d0

    double precision, dimension(0:2) :: timer_6dWsolver=0.d0,timer_6dWsolver_1RHSnAk=0.d0
    double precision, dimension(0:2) :: timer_6dWsolver_2TDMAk=0.d0
    double precision, dimension(0:2) :: timer_6dWsolver_3Trans1=0.d0
    double precision, dimension(0:2) :: timer_6dWsolver_4Ai=0.d0
    double precision, dimension(0:2) :: timer_6dWsolver_5TDMAi=0.d0
    double precision, dimension(0:2) :: timer_6dWsolver_6Trans2=0.d0
    double precision, dimension(0:2) :: timer_6dWsolver_7Aj=0.d0
    double precision, dimension(0:2) :: timer_6dWsolver_8TDMAj=0.d0
    double precision, dimension(0:2) :: timer_6dWsolver_9Update=0.d0

contains
    
    subroutine debug_infor(n1sub_in,n2sub_in,n3sub_in,period_in             &
                          ,ista_in,iend_in,jsta_in,jend_in,ksta_in,kend_in  &
                          ,x1_sub_in, x2_sub_in, x3_sub_in                  &
                          ,nprocs_1_in,nprocs_2_in,nprocs_3_in              &
                          ,myrank_1_in,myrank_2_in,myrank_3_in,myrank_in    )
        implicit none
        integer :: n1sub_in,n2sub_in,n3sub_in
        logical :: period_in(0:2)
        integer :: ista_in,iend_in,jsta_in,jend_in,ksta_in,kend_in
        double precision :: x1_sub_in(0:n1sub_in), x2_sub_in(0:n2sub_in), x3_sub_in(0:n3sub_in)
        integer :: nprocs_1_in,nprocs_2_in,nprocs_3_in
        integer :: myrank_1_in,myrank_2_in,myrank_3_in,myrank_in
        
        integer :: r1_0,r2_0,r3_0
        integer :: r1_1,r2_1,r3_1
        integer :: r1_2,r2_2,r3_2
        integer :: r1_3,r2_3,r3_3
        integer :: r1_4,r2_4,r3_4
        
        allocate(y1_sub(0:n1sub_in), y2_sub(0:n2sub_in), y3_sub(0:n3sub_in))
        allocate(debug_tmp1(0:n1sub_in,0:n2sub_in,0:n3sub_in))
        allocate(debug_tmp2(0:n1sub_in,0:n2sub_in,0:n3sub_in))
        allocate(debug_tmp3(0:n1sub_in,0:n2sub_in,0:n3sub_in))
        allocate(debug_tmp4(0:n1sub_in,0:n2sub_in,0:n3sub_in))
        allocate(debug_tmp5(0:n1sub_in,0:n2sub_in,0:n3sub_in))

        debug_tmp1=0.d0
        debug_tmp2=0.d0
        debug_tmp3=0.d0
        debug_tmp4=0.d0
        debug_tmp5=0.d0
        
        Dn1sub=n1sub_in
        Dn2sub=n2sub_in
        Dn3sub=n3sub_in
        Dperiod(0:2)=period_in(0:2)
        
        Dista=ista_in; Diend=iend_in
        Djsta=jsta_in; Djend=jend_in
        Dksta=ksta_in; Dkend=kend_in
        
        y1_sub(0:n1sub_in)=x1_sub_in(0:n1sub_in)
        y2_sub(0:n2sub_in)=x2_sub_in(0:n2sub_in)
        y3_sub(0:n3sub_in)=x3_sub_in(0:n3sub_in)

        Dnprocs_1=nprocs_1_in
        Dnprocs_2=nprocs_2_in
        Dnprocs_3=nprocs_3_in
        Dmyrank_1=myrank_1_in
        Dmyrank_2=myrank_2_in
        Dmyrank_3=myrank_3_in
        Dmyrank=myrank_in

        r1_4=int((Dmyrank_1)/10000)
        r1_3=int((Dmyrank_1-r1_4*10000)/1000 )
        r1_2=int((Dmyrank_1-r1_4*10000-r1_3*1000)/100 )
        r1_1=int((Dmyrank_1-r1_4*10000-r1_3*1000-r1_2*100)/10 )
        r1_0=int((Dmyrank_1-r1_4*10000-r1_3*1000-r1_2*100-r1_1*10)/1)

        r2_4=int((Dmyrank_2)/10000)
        r2_3=int((Dmyrank_2-r2_4*10000)/1000 )
        r2_2=int((Dmyrank_2-r2_4*10000-r2_3*1000)/100 )
        r2_1=int((Dmyrank_2-r2_4*10000-r2_3*1000-r2_2*100)/10 )
        r2_0=int((Dmyrank_2-r2_4*10000-r2_3*1000-r2_2*100-r2_1*10)/1)
        
        r3_4=int((Dmyrank_3)/10000)
        r3_3=int((Dmyrank_3-r3_4*10000)/1000 )
        r3_2=int((Dmyrank_3-r3_4*10000-r3_3*1000)/100 )
        r3_1=int((Dmyrank_3-r3_4*10000-r3_3*1000-r3_2*100)/10 )
        r3_0=int((Dmyrank_3-r3_4*10000-r3_3*1000-r3_2*100-r3_1*10)/1)

        
        write(filerank, '(5I1,1A1,5I1,1A1,5I1)' ) r1_4,r1_3,r1_2,r1_1,r1_0,'_'   &
                                                 ,r2_4,r2_3,r2_2,r2_1,r2_0,'_'   &
                                                 ,r3_4,r3_3,r3_2,r3_1,r3_0       


    end subroutine debug_infor

    subroutine debug_fileout(cha,Dtimestep,Dprint_range_in,V1,V2,V3,V4,V5)
        implicit none
        character(len=6) :: cha
        integer :: Dtimestep,Dprint_range_in(1:2,1:3)
        double precision :: V1(0:Dn1sub,0:Dn2sub,0:Dn3sub),V2(0:Dn1sub,0:Dn2sub,0:Dn3sub),V3(0:Dn1sub,0:Dn2sub,0:Dn3sub)
        double precision :: V4(0:Dn1sub,0:Dn2sub,0:Dn3sub),V5(0:Dn1sub,0:Dn2sub,0:Dn3sub)
                
        integer :: i,j,k
        
        if ( Dmyrank_1>=Dprint_range_in(1,1).and.Dmyrank_1<=Dprint_range_in(2,1)  &
            .and. Dmyrank_2>=Dprint_range_in(1,2).and.Dmyrank_2<=Dprint_range_in(2,2)  &
            .and. Dmyrank_3>=Dprint_range_in(1,3).and.Dmyrank_3<=Dprint_range_in(2,3)  ) then
            !=====================================================
                
                open(unit=Dmyrank,file=cha//filerank//'.plt', position='append')
                        write(Dmyrank,*) 'VARIABLES="X","Y","Z","U","V","W","P","T"'  !--
                        write(Dmyrank,*) 'zone t="',filerank,'"','i=',Dn1sub+1,'j=',Dn2sub+1,'k=',Dn3sub+1
                        write(Dmyrank,*) 'STRANDID=', Dtimestep
                        write(Dmyrank,*) 'SOLUTIONTIME=', Dtimestep
                        do k=0,Dn3sub
                        do j=0,Dn2sub
                        do i=0,Dn1sub
                            write(Dmyrank,'(3D20.10,5D30.20)') y1_sub(i),y2_sub(j),y3_sub(k)    &
                                                        &, V1(i,j,k),V2(i,j,k),V3(i,j,k),V4(i,j,k),V5(i,j,k)
                        enddo
                        enddo
                        enddo
                close(Dmyrank)  
                
            !=====================================================
        end if
    end subroutine debug_fileout

    subroutine debug_fileout2(cha,Dtimestep,Dprint_range_in,V1,V2,V3,V4,V5)
        implicit none
        character(len=6) :: cha
        integer :: Dtimestep,Dprint_range_in(1:2,1:3)
        double precision :: V1(0:Dn1sub,0:Dn2sub,0:Dn3sub),V2(0:Dn1sub,0:Dn2sub,0:Dn3sub),V3(0:Dn1sub,0:Dn2sub,0:Dn3sub)
        double precision :: V4(0:Dn1sub,0:Dn2sub,0:Dn3sub),V5(0:Dn1sub,0:Dn2sub,0:Dn3sub)
                
        integer :: i,j,k
        
        if ( Dmyrank_1>=Dprint_range_in(1,1).and.Dmyrank_1<=Dprint_range_in(2,1)  &
            .and. Dmyrank_2>=Dprint_range_in(1,2).and.Dmyrank_2<=Dprint_range_in(2,2)  &
            .and. Dmyrank_3>=Dprint_range_in(1,3).and.Dmyrank_3<=Dprint_range_in(2,3)  ) then
            !=====================================================
                
                open(unit=Dmyrank,file=cha//filerank//'.plt', position='append')
                        write(Dmyrank,*) 'VARIABLES="X","Y","Z","U","V","W","P","T"'  !--
                        write(Dmyrank,*) 'zone t="',filerank,'"','i=',Dn1sub-1,'j=',Dn2sub-1,'k=',Dn3sub-1
                        write(Dmyrank,*) 'STRANDID=', Dtimestep
                        write(Dmyrank,*) 'SOLUTIONTIME=', Dtimestep
                        do k=1,Dn3sub-1
                        do j=1,Dn2sub-1
                        do i=1,Dn1sub-1
                            ! write(Dmyrank,'(3I4,5D30.20)') i,j,k    &
                            write(Dmyrank,'(3D20.10,5D30.20)') y1_sub(i),y2_sub(j),y3_sub(k)    &
                                                        &, V1(i,j,k),V2(i,j,k),V3(i,j,k),V4(i,j,k),V5(i,j,k)
                            ! if(abs(V3(i,j,k))>1.0d-10) write(Dmyrank,'(3I4,5D30.20)') i,j,k, V1(i,j,k),V2(i,j,k),V3(i,j,k),V4(i,j,k),V5(i,j,k)
                        enddo
                        enddo
                        enddo
                close(Dmyrank)  
                
            !=====================================================
        end if
    end subroutine debug_fileout2
    subroutine debug_clean()
        implicit none
        deallocate(y1_sub, y2_sub, y3_sub)
        deallocate(debug_tmp1)
        deallocate(debug_tmp2)
        deallocate(debug_tmp3)
        deallocate(debug_tmp4)
        deallocate(debug_tmp5)
    end subroutine debug_clean
    
    subroutine debug_mpitime(times,Checker)
        use mpi
        implicit none
        double precision :: times(0:2)
        integer :: Checker
        double precision :: timer

        times(Checker)=MPI_WTIME()
        times(2)=times(2)+times(Checker)-times(0)

    end subroutine debug_mpitime

end module debug
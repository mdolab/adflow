!
!      ******************************************************************
!      *                                                                *
!      * File:          wOverset.f90                                    *
!      * Author:        Steve Repsher                                   *
!      * Starting date: 02-04-2005                                      *
!      * Last modified: 09-16-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine wOverset(level, start, end, commPressure,       &
                           commVarGamma, commLamVis, commEddyVis, &
                           commPattern, internal, nlev)
!
!      ******************************************************************
!      *                                                                *
!      * wOverset controls the communication between overset halos      *
!      * for the cell-centered variables by interpolating the solution  *
!      * from other blocks consistent with the chimera approach. A tri- *
!      * linear interpolation is used as per the input from cgns. It    *
!      * is possible to send a range of variables and not the entire    *
!      * set, e.g. only the flow variables or only the turbulent        *
!      * variables. This is controlled by the arguments start, end,     *
!      * commPressure and commViscous. The exchange takes place for     *
!      * the given grid level.                                          *
!      *                                                                *
!      ******************************************************************
!
       use block
       use communication
       use inputTimeSpectral
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType), intent(in) :: level, start, end, nlev
       logical, intent(in) :: commPressure, commVarGamma
       logical, intent(in) :: commLamVis, commEddyVis

       type(commType), dimension(nlev,*), intent(in) :: commPattern
       type(internalCommType), dimension(nlev,*), intent(in) :: internal
!
!      Local variables.
!
       integer :: size, procId, ierr, index
       integer, dimension(mpi_status_size) :: status

       integer(kind=intType) :: nVar, mm
       integer(kind=intType) :: i, j, k, ii, jj
       integer(kind=intType) :: d1, i1, j1, k1, d2, i2, j2, k2

       real(kind=realType), dimension(:), pointer :: weight
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Determine the number of variables per cell to be sent.

       nVar = max(0_intType,(end - start + 1))
       if( commPressure ) nVar = nVar + 1
       if( commVarGamma ) nVar = nVar + 1
       if( commLamVis )   nVar = nVar + 1
       if( commEddyVis )  nVar = nVar + 1

       if(nVar == 0) return

       ! Loop over the number of spectral solutions.

       spectralModes: do mm=1,nTimeIntervalsSpectral

         ! Send the variables. The data is first copied into
         ! the send buffer after which the buffer is sent asap.

         ii = 1
         sends: do i=1,commPattern(level,mm)%nProcSend

           ! Store the processor id and the size of the message
           ! a bit easier.

           procID = commPattern(level,mm)%sendProc(i)
           size    = nVar*commPattern(level,mm)%nsend(i)

           ! Copy the data in the correct part of the send buffer.

           jj = ii
           do j=1,commPattern(level,mm)%nsend(i)

             ! Store the block id and the indices of the donor
             ! a bit easier.

             d1 = commPattern(level,mm)%sendList(i)%block(j)
             i1 = commPattern(level,mm)%sendList(i)%indices(j,1)
             j1 = commPattern(level,mm)%sendList(i)%indices(j,2)
             k1 = commPattern(level,mm)%sendList(i)%indices(j,3)

             weight => commPattern(level,mm)%sendList(i)%interp(j,:)

             ! Copy the given range of the working variables for
             ! this cell in the buffer. Update the counter jj.

             do k=start,end
               sendBuffer(jj) = &
                  weight(1)*flowDoms(d1,level,mm)%w(i1  ,j1  ,k1  ,k) + &
                  weight(2)*flowDoms(d1,level,mm)%w(i1+1,j1  ,k1  ,k) + &
                  weight(3)*flowDoms(d1,level,mm)%w(i1  ,j1+1,k1  ,k) + &
                  weight(4)*flowDoms(d1,level,mm)%w(i1+1,j1+1,k1  ,k) + &
                  weight(5)*flowDoms(d1,level,mm)%w(i1  ,j1  ,k1+1,k) + &
                  weight(6)*flowDoms(d1,level,mm)%w(i1+1,j1  ,k1+1,k) + &
                  weight(7)*flowDoms(d1,level,mm)%w(i1  ,j1+1,k1+1,k) + &
                  weight(8)*flowDoms(d1,level,mm)%w(i1+1,j1+1,k1+1,k)
               jj = jj + 1
             enddo

             ! The pressure, if needed.

             if( commPressure ) then
               sendBuffer(jj) = &
                  weight(1)*flowDoms(d1,level,mm)%p(i1  ,j1  ,k1  ) + &
                  weight(2)*flowDoms(d1,level,mm)%p(i1+1,j1  ,k1  ) + &
                  weight(3)*flowDoms(d1,level,mm)%p(i1  ,j1+1,k1  ) + &
                  weight(4)*flowDoms(d1,level,mm)%p(i1+1,j1+1,k1  ) + &
                  weight(5)*flowDoms(d1,level,mm)%p(i1  ,j1  ,k1+1) + &
                  weight(6)*flowDoms(d1,level,mm)%p(i1+1,j1  ,k1+1) + &
                  weight(7)*flowDoms(d1,level,mm)%p(i1  ,j1+1,k1+1) + &
                  weight(8)*flowDoms(d1,level,mm)%p(i1+1,j1+1,k1+1)
               jj = jj + 1
             endif

             ! The specific heat ratio, if needed. Note that level == 1.

             if( commVarGamma ) then
               sendBuffer(jj) = &
                  weight(1)*flowDoms(d1,1,mm)%gamma(i1  ,j1  ,k1  ) + &
                  weight(2)*flowDoms(d1,1,mm)%gamma(i1+1,j1  ,k1  ) + &
                  weight(3)*flowDoms(d1,1,mm)%gamma(i1  ,j1+1,k1  ) + &
                  weight(4)*flowDoms(d1,1,mm)%gamma(i1+1,j1+1,k1  ) + &
                  weight(5)*flowDoms(d1,1,mm)%gamma(i1  ,j1  ,k1+1) + &
                  weight(6)*flowDoms(d1,1,mm)%gamma(i1+1,j1  ,k1+1) + &
                  weight(7)*flowDoms(d1,1,mm)%gamma(i1  ,j1+1,k1+1) + &
                  weight(8)*flowDoms(d1,1,mm)%gamma(i1+1,j1+1,k1+1)
               jj = jj + 1
             endif

             ! The laminar viscosity for a viscous computation.
             ! Again level == 1.

             if( commLamVis ) then
               sendBuffer(jj) = &
                  weight(1)*flowDoms(d1,1,mm)%rlv(i1  ,j1  ,k1  ) + &
                  weight(2)*flowDoms(d1,1,mm)%rlv(i1+1,j1  ,k1  ) + &
                  weight(3)*flowDoms(d1,1,mm)%rlv(i1  ,j1+1,k1  ) + &
                  weight(4)*flowDoms(d1,1,mm)%rlv(i1+1,j1+1,k1  ) + &
                  weight(5)*flowDoms(d1,1,mm)%rlv(i1  ,j1  ,k1+1) + &
                  weight(6)*flowDoms(d1,1,mm)%rlv(i1+1,j1  ,k1+1) + &
                  weight(7)*flowDoms(d1,1,mm)%rlv(i1  ,j1+1,k1+1) + &
                  weight(8)*flowDoms(d1,1,mm)%rlv(i1+1,j1+1,k1+1)
               jj = jj + 1
             endif

             ! The eddy viscosity for eddy viscosity models.
             ! Level is the true multigrid level, because the eddy
             ! viscosity is allocated on all grid levels.

             if( commEddyVis ) then
               sendBuffer(jj) = &
                  weight(1)*flowDoms(d1,level,mm)%rev(i1  ,j1  ,k1  ) + &
                  weight(2)*flowDoms(d1,level,mm)%rev(i1+1,j1  ,k1  ) + &
                  weight(3)*flowDoms(d1,level,mm)%rev(i1  ,j1+1,k1  ) + &
                  weight(4)*flowDoms(d1,level,mm)%rev(i1+1,j1+1,k1  ) + &
                  weight(5)*flowDoms(d1,level,mm)%rev(i1  ,j1  ,k1+1) + &
                  weight(6)*flowDoms(d1,level,mm)%rev(i1+1,j1  ,k1+1) + &
                  weight(7)*flowDoms(d1,level,mm)%rev(i1  ,j1+1,k1+1) + &
                  weight(8)*flowDoms(d1,level,mm)%rev(i1+1,j1+1,k1+1)
               jj = jj + 1
             endif

           enddo

           ! Send the data.

           call mpi_isend(sendBuffer(ii), size, sumb_real, procId,  &
                          procId, SUmb_comm_world, sendRequests(i), &
                          ierr)

           ! Set ii to jj for the next processor.

           ii = jj

         enddo sends

         ! Post the nonblocking receives.

         ii = 1
         receives: do i=1,commPattern(level,mm)%nProcRecv

           ! Store the processor id and the size of the message
           ! a bit easier.

           procID = commPattern(level,mm)%recvProc(i)
           size    = nVar*commPattern(level,mm)%nrecv(i)

           ! Post the receive.

           call mpi_irecv(recvBuffer(ii), size, sumb_real, procId, &
                          myId, SUmb_comm_world, recvRequests(i), ierr)

           ! And update ii.

           ii = ii + size

         enddo receives

         ! Do the local interpolation.

         localInterp: do i=1,internal(level,mm)%ncopy

           ! Store the block and the indices of the donor a bit easier.

           d1 = internal(level,mm)%donorBlock(i)
           i1 = internal(level,mm)%donorIndices(i,1)
           j1 = internal(level,mm)%donorIndices(i,2)
           k1 = internal(level,mm)%donorIndices(i,3)

           weight => internal(level,mm)%donorInterp(i,:)

           ! Idem for the halo's.

           d2 = internal(level,mm)%haloBlock(i)
           i2 = internal(level,mm)%haloIndices(i,1)
           j2 = internal(level,mm)%haloIndices(i,2)
           k2 = internal(level,mm)%haloIndices(i,3)

           ! Copy the given range of working variables.

           do k=start,end
             flowDoms(d2,level,mm)%w(i2,j2,k2,k) = &
                  weight(1)*flowDoms(d1,level,mm)%w(i1  ,j1  ,k1  ,k) + &
                  weight(2)*flowDoms(d1,level,mm)%w(i1+1,j1  ,k1  ,k) + &
                  weight(3)*flowDoms(d1,level,mm)%w(i1  ,j1+1,k1  ,k) + &
                  weight(4)*flowDoms(d1,level,mm)%w(i1+1,j1+1,k1  ,k) + &
                  weight(5)*flowDoms(d1,level,mm)%w(i1  ,j1  ,k1+1,k) + &
                  weight(6)*flowDoms(d1,level,mm)%w(i1+1,j1  ,k1+1,k) + &
                  weight(7)*flowDoms(d1,level,mm)%w(i1  ,j1+1,k1+1,k) + &
                  weight(8)*flowDoms(d1,level,mm)%w(i1+1,j1+1,k1+1,k)
           enddo

           ! The pressure, if needed.

           if( commPressure ) then
             flowDoms(d2,level,mm)%p(i2,j2,k2) = &
                  weight(1)*flowDoms(d1,level,mm)%p(i1  ,j1  ,k1  ) + &
                  weight(2)*flowDoms(d1,level,mm)%p(i1+1,j1  ,k1  ) + &
                  weight(3)*flowDoms(d1,level,mm)%p(i1  ,j1+1,k1  ) + &
                  weight(4)*flowDoms(d1,level,mm)%p(i1+1,j1+1,k1  ) + &
                  weight(5)*flowDoms(d1,level,mm)%p(i1  ,j1  ,k1+1) + &
                  weight(6)*flowDoms(d1,level,mm)%p(i1+1,j1  ,k1+1) + &
                  weight(7)*flowDoms(d1,level,mm)%p(i1  ,j1+1,k1+1) + &
                  weight(8)*flowDoms(d1,level,mm)%p(i1+1,j1+1,k1+1)
           end if

           ! The specific heat ratio, if needed. Note that level == 1.

           if( commVarGamma ) then
             flowDoms(d2,1,mm)%gamma(i2,j2,k2) = &
                  weight(1)*flowDoms(d1,1,mm)%gamma(i1  ,j1  ,k1  ) + &
                  weight(2)*flowDoms(d1,1,mm)%gamma(i1+1,j1  ,k1  ) + &
                  weight(3)*flowDoms(d1,1,mm)%gamma(i1  ,j1+1,k1  ) + &
                  weight(4)*flowDoms(d1,1,mm)%gamma(i1+1,j1+1,k1  ) + &
                  weight(5)*flowDoms(d1,1,mm)%gamma(i1  ,j1  ,k1+1) + &
                  weight(6)*flowDoms(d1,1,mm)%gamma(i1+1,j1  ,k1+1) + &
                  weight(7)*flowDoms(d1,1,mm)%gamma(i1  ,j1+1,k1+1) + &
                  weight(8)*flowDoms(d1,1,mm)%gamma(i1+1,j1+1,k1+1)
           end if

           ! The laminar viscosity for viscous computations.
           ! Again level == 1.

           if( commLamVis ) then
             flowDoms(d2,1,mm)%rlv(i2,j2,k2) = &
                  weight(1)*flowDoms(d1,1,mm)%rlv(i1  ,j1  ,k1  ) + &
                  weight(2)*flowDoms(d1,1,mm)%rlv(i1+1,j1  ,k1  ) + &
                  weight(3)*flowDoms(d1,1,mm)%rlv(i1  ,j1+1,k1  ) + &
                  weight(4)*flowDoms(d1,1,mm)%rlv(i1+1,j1+1,k1  ) + &
                  weight(5)*flowDoms(d1,1,mm)%rlv(i1  ,j1  ,k1+1) + &
                  weight(6)*flowDoms(d1,1,mm)%rlv(i1+1,j1  ,k1+1) + &
                  weight(7)*flowDoms(d1,1,mm)%rlv(i1  ,j1+1,k1+1) + &
                  weight(8)*flowDoms(d1,1,mm)%rlv(i1+1,j1+1,k1+1)
           end if

           ! The eddy viscosity for eddy viscosity models.
           ! Level is the true multigrid level, because the eddy
           ! viscosity is allocated on all grid levels.

           if( commEddyVis ) then
             flowDoms(d2,level,mm)%rev(i2,j2,k2) = &
                  weight(1)*flowDoms(d1,level,mm)%rev(i1  ,j1  ,k1  ) + &
                  weight(2)*flowDoms(d1,level,mm)%rev(i1+1,j1  ,k1  ) + &
                  weight(3)*flowDoms(d1,level,mm)%rev(i1  ,j1+1,k1  ) + &
                  weight(4)*flowDoms(d1,level,mm)%rev(i1+1,j1+1,k1  ) + &
                  weight(5)*flowDoms(d1,level,mm)%rev(i1  ,j1  ,k1+1) + &
                  weight(6)*flowDoms(d1,level,mm)%rev(i1+1,j1  ,k1+1) + &
                  weight(7)*flowDoms(d1,level,mm)%rev(i1  ,j1+1,k1+1) + &
                  weight(8)*flowDoms(d1,level,mm)%rev(i1+1,j1+1,k1+1)
           end if

         enddo localInterp

         ! Complete the nonblocking receives in an arbitrary sequence and
         ! copy the variables from the buffer into the halo's.

         size = commPattern(level,mm)%nProcRecv
         completeRecvs: do i=1,commPattern(level,mm)%nProcRecv

           ! Complete any of the requests.

           call mpi_waitany(size, recvRequests, index, status, ierr)

           ! Copy the data just arrived in the halo's.

           ii = index
           jj = nVar*commPattern(level,mm)%nrecvCum(ii-1)
           do j=1,commPattern(level,mm)%nrecv(ii)

             ! Store the block and the indices of the halo a bit easier.

             d2 = commPattern(level,mm)%recvList(ii)%block(j)
             i2 = commPattern(level,mm)%recvList(ii)%indices(j,1)
             j2 = commPattern(level,mm)%recvList(ii)%indices(j,2)
             k2 = commPattern(level,mm)%recvList(ii)%indices(j,3)

             ! Copy the conservative variables.

             do k=start,end
               jj = jj + 1
               flowDoms(d2,level,mm)%w(i2,j2,k2,k) = recvBuffer(jj)
             enddo

             ! The pressure, if needed.

             if( commPressure ) then
               jj = jj + 1
               flowDoms(d2,level,mm)%p(i2,j2,k2) = recvBuffer(jj)
             endif

             ! The specific heat ratio, if needed. Note that level == 1.

             if( commVarGamma ) then
               jj = jj + 1
               flowDoms(d2,1,mm)%gamma(i2,j2,k2) = recvBuffer(jj)
             endif

             ! The laminar viscosity for viscous computations.
             ! Again level == 1.

             if( commLamVis ) then
               jj = jj + 1
               flowDoms(d2,1,mm)%rlv(i2,j2,k2) = recvBuffer(jj)
             endif

             ! The eddy viscosity ratio for eddy viscosity models.
             ! Level is the true multigrid level, because the eddy
             ! viscosity is allocated on all grid levels.

             if( commEddyVis ) then
               jj = jj + 1
               flowDoms(d2,level,mm)%rev(i2,j2,k2) = recvBuffer(jj)
             endif

            enddo

         enddo completeRecvs

         ! Complete the nonblocking sends.

         size = commPattern(level,mm)%nProcSend
         do i=1,commPattern(level,mm)%nProcSend
           call mpi_waitany(size, sendRequests, index, status, ierr)
         enddo

       enddo spectralModes

       end subroutine wOverset

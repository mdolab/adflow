!
!      ******************************************************************
!      *                                                                *
!      * File:          indexHaloNode1to1.f90                               *
!      * Author:        C.A.(Sandy) Mader                               *
!      * Starting date: 07-07-2008                                      *
!      * Last modified: 07-18-2008                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine indexHaloNode1to1(level)
!
!      ******************************************************************
!      *                                                                *
!      * indexHalo1to1 exchanges the global node indicies for the       *
!      * 1 to 1 matching subfaces. based on exchangecoor.f90            *
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
       integer(kind=intType), intent(in) :: level

!
!      Local variables.
!
       integer :: size, procID, ierr, index
       integer, dimension(mpi_status_size) :: status

       integer(kind=intType) :: mm
       integer(kind=intType) :: i, j, ii, jj
       integer(kind=intType) :: d1, i1, j1, k1, d2, i2, j2, k2


!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!

       ! Loop over the number of spectral solutions.

       spectralModes: do mm=1,nTimeIntervalsSpectral

         ! Send the variables. The data is first copied into
         ! the send buffer after which the buffer is sent asap.

         ii = 1
         sends: do i=1,commPatternNode_1st(level)%nProcSend

           ! Store the processor id and the size of the message
           ! a bit easier.

           procID = commPatternNode_1st(level)%sendProc(i)
           size    = 1*commPatternNode_1st(level)%nsend(i)

           ! Copy the data in the correct part of the send buffer.

           jj = ii
           do j=1,commPatternNode_1st(level)%nsend(i)

             ! Store the block id and the indices of the donor
             ! a bit easier.

             d1 = commPatternNode_1st(level)%sendList(i)%block(j)
             i1 = commPatternNode_1st(level)%sendList(i)%indices(j,1)
             j1 = commPatternNode_1st(level)%sendList(i)%indices(j,2)
             k1 = commPatternNode_1st(level)%sendList(i)%indices(j,3)

             ! Copy the index for
             ! this cell in the buffer. Update the counter jj.

             sendBuffer(jj) = flowDoms(d1,level,mm)%globalNode(i1,j1,k1)
             jj = jj + 1
             

           enddo

           ! Send the data.

           call mpi_isend(sendBuffer(ii), size, sumb_real, procID,  &
                          procID, SUmb_comm_world, sendRequests(i), &
                          ierr)

           ! Set ii to jj for the next processor.

           ii = jj

         enddo sends

         ! Post the nonblocking receives.

         ii = 1
         receives: do i=1,commPatternNode_1st(level)%nProcRecv

           ! Store the processor id and the size of the message
           ! a bit easier.

           procID = commPatternNode_1st(level)%recvProc(i)
           size    = 1*commPatternNode_1st(level)%nrecv(i)

           ! Post the receive.

           call mpi_irecv(recvBuffer(ii), size, sumb_real, procID, &
                          myID, SUmb_comm_world, recvRequests(i), ierr)

           ! And update ii.

           ii = ii + size

         enddo receives

         ! Copy the local data.

         localCopy: do i=1,internalNode_1st(level)%ncopy

           ! Store the block and the indices of the donor a bit easier.

           d1 = internalNode_1st(level)%donorBlock(i)
           i1 = internalNode_1st(level)%donorIndices(i,1)
           j1 = internalNode_1st(level)%donorIndices(i,2)
           k1 = internalNode_1st(level)%donorIndices(i,3)

           ! Idem for the halo's.

           d2 = internalNode_1st(level)%haloBlock(i)
           i2 = internalNode_1st(level)%haloIndices(i,1)
           j2 = internalNode_1st(level)%haloIndices(i,2)
           k2 = internalNode_1st(level)%haloIndices(i,3)

           ! Copy the given range of indicies.

           flowDoms(d2,level,mm)%globalNode(i2,j2,k2) = &
                flowDoms(d1,level,mm)%globalNode(i1,j1,k1)

         enddo localCopy
!!$
!!$         ! Correct the periodic halo's of the internal communication
!!$         ! pattern, if needed.
!!$
!!$            call correctPeriodicCoor(level, mm,                 &
!!$                                        internalNode_1st(level)%nPeriodic, &
!!$                                        internalNode_1st(level)%periodicData)

         ! Complete the nonblocking receives in an arbitrary sequence and
         ! copy the variables from the buffer into the halo's.

         size = commPatternNode_1st(level)%nProcRecv
         completeRecvs: do i=1,commPatternNode_1st(level)%nProcRecv

           ! Complete any of the requests.

           call mpi_waitany(size, recvRequests, index, status, ierr)

           ! Copy the data just arrived in the halo's.

           ii = index
           jj = 1*commPatternNode_1st(level)%nrecvCum(ii-1)
           do j=1,commPatternNode_1st(level)%nrecv(ii)

             ! Store the block and the indices of the halo a bit easier.

             d2 = commPatternNode_1st(level)%recvList(ii)%block(j)
             i2 = commPatternNode_1st(level)%recvList(ii)%indices(j,1)
             j2 = commPatternNode_1st(level)%recvList(ii)%indices(j,2)
             k2 = commPatternNode_1st(level)%recvList(ii)%indices(j,3)

             ! Copy the indicies.

             jj = jj + 1
             flowDoms(d2,level,mm)%globalNode(i2,j2,k2) = recvBuffer(jj)
            

            enddo

         enddo completeRecvs
!!$
!!$         ! Correct the periodic halo's of the external communication
!!$         ! pattern, if needed.
!!$
!!$           call correctPeriodicCoor(level, mm,                    &
!!$                                        commPatternNode_1st(level)%nPeriodic, &
!!$                                        commPatternNode_1st(level)%periodicData)

         ! Complete the nonblocking sends.

         size = commPatternNode_1st(level)%nProcSend
         do i=1,commPatternNode_1st(level)%nProcSend
           call mpi_waitany(size, sendRequests, index, status, ierr)
         enddo

       enddo spectralModes

     end subroutine indexHaloNode1to1

!      ==================================================================


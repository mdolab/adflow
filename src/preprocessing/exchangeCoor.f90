!
!      ******************************************************************
!      *                                                                *
!      * File:          exchangeCoor.F90                                *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 02-24-2003                                      *
!      * Last modified: 06-12-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine exchangeCoor(level)
!
!      ******************************************************************
!      *                                                                *
!      * ExchangeCoor exchanges the coordinates of the given grid       *
!      * level.                                                         *
!      *                                                                *
!      ******************************************************************
!
       use block
       use communication
       use inputTimeSpectral
!
!      Subroutine arguments.
!
       integer(kind=intType), intent(in) :: level
!
!      Local variables.
!
       integer :: size, procID, ierr, index
       integer, dimension(mpi_status_size) :: status

       integer(kind=intType) :: i, j, ii, jj, mm
       integer(kind=intType) :: d1, i1, j1, k1, d2, i2, j2, k2
!
!      Interfaces
!
       interface
         subroutine correctPeriodicCoor(level, sp, nPeriodic, &
                                        periodicData)
         use block
         use communication
         implicit none

         integer(kind=intType), intent(in) :: level, sp, nPeriodic
         type(periodicDataType), dimension(:), pointer :: periodicData

         end subroutine correctPeriodicCoor
       end interface
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Loop over the number of spectral solutions.

       spectralLoop: do mm=1,nTimeIntervalsSpectral

         ! Send the coordinates i have to send. The data is first copied
         ! into the send buffer and this buffer is sent.

         ii = 1
         sends: do i=1,commPatternNode_1st(level)%nProcSend

           ! Store the processor id and the size of the message
           ! a bit easier.

           procID = commPatternNode_1st(level)%sendProc(i)
           size   = 3*commPatternNode_1st(level)%nSend(i)

           ! Copy the data in the correct part of the send buffer.

           jj = ii
           do j=1,commPatternNode_1st(level)%nSend(i)

             ! Store the block id and the indices of the donor
             ! a bit easier.

             d1 = commPatternNode_1st(level)%sendList(i)%block(j)
             i1 = commPatternNode_1st(level)%sendList(i)%indices(j,1)
             j1 = commPatternNode_1st(level)%sendList(i)%indices(j,2)
             k1 = commPatternNode_1st(level)%sendList(i)%indices(j,3)

             ! Copy the coordinates of this point in the buffer.
             ! Update the counter jj accordingly.

             sendBuffer(jj)   = flowDoms(d1,level,mm)%x(i1,j1,k1,1)
             sendBuffer(jj+1) = flowDoms(d1,level,mm)%x(i1,j1,k1,2)
             sendBuffer(jj+2) = flowDoms(d1,level,mm)%x(i1,j1,k1,3)
             jj = jj + 3

           enddo

           ! Send the data.

           call mpi_isend(sendBuffer(ii), size, sumb_real, procID,    &
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
           size   = 3*commPatternNode_1st(level)%nRecv(i)

           ! Post the receive.

           call mpi_irecv(recvBuffer(ii), size, sumb_real, procID, &
                           myID, SUmb_comm_world, recvRequests(i), ierr)

           ! And update ii.

           ii = ii + size

         enddo receives

         ! Copy the local data.

         localCopy: do i=1,internalNode_1st(level)%nCopy

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

           ! Copy the coordinates.

           flowDoms(d2,level,mm)%x(i2,j2,k2,1) = &
                                flowDoms(d1,level,mm)%x(i1,j1,k1,1)
           flowDoms(d2,level,mm)%x(i2,j2,k2,2) = &
                                flowDoms(d1,level,mm)%x(i1,j1,k1,2)
           flowDoms(d2,level,mm)%x(i2,j2,k2,3) = &
                                flowDoms(d1,level,mm)%x(i1,j1,k1,3)

         enddo localCopy

         ! Correct the periodic halos of the internal communication
         ! pattern

         call correctPeriodicCoor(level, mm,                          &
                                  internalNode_1st(level)%nPeriodic,  &
                                  internalNode_1st(level)%periodicData)

         ! Complete the nonblocking receives in an arbitrary sequence and
         ! copy the coordinates from the buffer into the halo's.

         size = commPatternNode_1st(level)%nProcRecv
         completeRecvs: do i=1,commPatternNode_1st(level)%nProcRecv

           ! Complete any of the requests.

           call mpi_waitany(size, recvRequests, index, status, ierr)

           ! Copy the data just arrived in the halo's.

           ii = index
           jj = 3*commPatternNode_1st(level)%nRecvCum(ii-1) +1
           do j=1,commPatternNode_1st(level)%nRecv(ii)

             ! Store the block and the indices of the halo a bit easier.

             d2 = commPatternNode_1st(level)%recvList(ii)%block(j)
             i2 = commPatternNode_1st(level)%recvList(ii)%indices(j,1)
             j2 = commPatternNode_1st(level)%recvList(ii)%indices(j,2)
             k2 = commPatternNode_1st(level)%recvList(ii)%indices(j,3)

             ! Copy the data.

             flowDoms(d2,level,mm)%x(i2,j2,k2,1) = recvBuffer(jj)
             flowDoms(d2,level,mm)%x(i2,j2,k2,2) = recvBuffer(jj+1)
             flowDoms(d2,level,mm)%x(i2,j2,k2,3) = recvBuffer(jj+2)
             jj = jj + 3

           enddo

         enddo completeRecvs

         ! Correct the periodic halos of the external communication
         ! pattern.

         call correctPeriodicCoor(level, mm,                            &
                                  commPatternNode_1st(level)%nPeriodic, &
                                  commPatternNode_1st(level)%periodicData)

         ! Complete the nonblocking sends.

         size = commPatternNode_1st(level)%nProcSend
         do i=1,commPatternNode_1st(level)%nProcSend
           call mpi_waitany(size, sendRequests, index, status, ierr)
         enddo

       enddo spectralLoop

       end subroutine exchangeCoor

!      ==================================================================

       subroutine correctPeriodicCoor(level, sp, nPeriodic, periodicData)
!
!      ******************************************************************
!      *                                                                *
!      * correctPeriodicCoor applies the periodic transformation to     *
!      * the coordinates of the nodal halo's in periodicData.           *
!      *                                                                *
!      ******************************************************************
!
       use block
       use communication
       implicit none
!
!      Subroutine arguments
!
       integer(kind=intType), intent(in) :: level, sp, nPeriodic
       type(periodicDataType), dimension(:), pointer :: periodicData
!
!      Local variables.
!
       integer(kind=intType) :: nn, mm, ii, i, j, k
       real(kind=realType)   :: dx, dy, dz

       real(kind=realType), dimension(3,3) :: rotMatrix
       real(kind=realType), dimension(3)   :: rotCenter, translation
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Loop over the number of periodic transformations.

       do nn=1,nPeriodic

         ! Store the rotation matrix, rotation center and translation
         ! vector a bit easier.

         rotMatrix   = periodicData(nn)%rotMatrix
         rotCenter   = periodicData(nn)%rotCenter
         translation = periodicData(nn)%translation + rotCenter

         ! Loop over the number of halo nodes for this transformation.

         do ii=1,periodicData(nn)%nHalos

           ! Store the block and the indices a bit easier.

           mm = periodicData(nn)%block(ii)
           i  = periodicData(nn)%indices(ii,1)
           j  = periodicData(nn)%indices(ii,2)
           k  = periodicData(nn)%indices(ii,3)

           ! Determine the vector from the center of rotation to the
           ! uncorrected halo value.

           dx = flowDoms(mm,level,sp)%x(i,j,k,1) - rotCenter(1)
           dy = flowDoms(mm,level,sp)%x(i,j,k,2) - rotCenter(2)
           dz = flowDoms(mm,level,sp)%x(i,j,k,3) - rotCenter(3)

           ! Compute the corrected coordinates.

           flowDoms(mm,level,sp)%x(i,j,k,1) = rotMatrix(1,1)*dx &
                                            + rotMatrix(1,2)*dy &
                                            + rotMatrix(1,3)*dz &
                                            + translation(1)
           flowDoms(mm,level,sp)%x(i,j,k,2) = rotMatrix(2,1)*dx &
                                            + rotMatrix(2,2)*dy &
                                            + rotMatrix(2,3)*dz &
                                            + translation(2)
           flowDoms(mm,level,sp)%x(i,j,k,3) = rotMatrix(3,1)*dx &
                                            + rotMatrix(3,2)*dy &
                                            + rotMatrix(3,3)*dz &
                                            + translation(3)
         enddo
       enddo

       end subroutine correctPeriodicCoor

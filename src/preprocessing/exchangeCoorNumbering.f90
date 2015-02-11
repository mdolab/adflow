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
subroutine exchangeCoorNumbering(level)
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
  integer :: ssize, procID, ierr, index
  integer, dimension(mpi_status_size) :: status

  integer(kind=intType) :: i, j, ii, jj, mm
  integer(kind=intType) :: d1, i1, j1, k1, d2, i2, j2, k2
  integer(kind=intType), allocatable, dimension(:) :: intRecvBuffer, intSendBuffer
  !
  !      ******************************************************************
  !      *                                                                *
  !      * Begin execution                                                *
  !      *                                                                *
  !      ******************************************************************
  !

  ! Allocate a new integer buffer. 
  allocate(intRecvBuffer(size(recvBuffer)), &
       intSendBuffer(size(sendBuffer)))


  ! Loop over the number of spectral solutions.
  spectralLoop: do mm=1,nTimeIntervalsSpectral

     ! Send the coordinates i have to send. The data is first copied
     ! into the send buffer and this buffer is sent.

     ii = 1
     sends: do i=1,commPatternNode_1st(level)%nProcSend

        ! Store the processor id and the size of the message
        ! a bit easier.

        procID = commPatternNode_1st(level)%sendProc(i)
        ssize   = commPatternNode_1st(level)%nSend(i)

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

           intSendBuffer(jj)   = flowDoms(d1,level,mm)%globalNode(i1,j1,k1)
           jj = jj + 1

        enddo

        ! Send the data.

        call mpi_isend(intSendBuffer(ii), ssize, sumb_integer, procID,    &
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
        ssize   = commPatternNode_1st(level)%nRecv(i)

        ! Post the receive.

        call mpi_irecv(intRecvBuffer(ii), ssize, sumb_integer, procID, &
             myID, SUmb_comm_world, recvRequests(i), ierr)

        ! And update ii.

        ii = ii + ssize

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
        flowDoms(d2,level,mm)%globalNode(i2,j2,k2) = &
             flowDoms(d1,level,mm)%globalNode(i1,j1,k1)

     enddo localCopy

     ! Complete the nonblocking receives in an arbitrary sequence and
     ! copy the coordinates from the buffer into the halo's.

     ssize = commPatternNode_1st(level)%nProcRecv
     completeRecvs: do i=1,commPatternNode_1st(level)%nProcRecv

        ! Complete any of the requests.

        call mpi_waitany(ssize, recvRequests, index, status, ierr)

        ! Copy the data just arrived in the halo's.

        ii = index
        jj = commPatternNode_1st(level)%nRecvCum(ii-1) +1
        do j=1,commPatternNode_1st(level)%nRecv(ii)

           ! Store the block and the indices of the halo a bit easier.

           d2 = commPatternNode_1st(level)%recvList(ii)%block(j)
           i2 = commPatternNode_1st(level)%recvList(ii)%indices(j,1)
           j2 = commPatternNode_1st(level)%recvList(ii)%indices(j,2)
           k2 = commPatternNode_1st(level)%recvList(ii)%indices(j,3)

           ! Copy the data.

           flowDoms(d2,level,mm)%globalNode(i2,j2,k2) = intRecvBuffer(jj)
           jj = jj + 1

        enddo

     enddo completeRecvs

     ! Complete the nonblocking sends.

     ssize = commPatternNode_1st(level)%nProcSend
     do i=1,commPatternNode_1st(level)%nProcSend
        call mpi_waitany(ssize, sendRequests, index, status, ierr)
     enddo

  enddo spectralLoop
  deallocate(intSendBuffer, intRecvBuffer)
end subroutine exchangeCoorNumbering

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
subroutine exchangeCellNumbering(level)
  !
  !      ******************************************************************
  !      *                                                                *
  !      * ExchangeCellNumbering exchanges the cell indices of the given  *
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
  spectralLoop: do mm=1, nTimeIntervalsSpectral

     ! Send the coordinates i have to send. The data is first copied
     ! into the send buffer and this buffer is sent.

     ii = 1
     sends: do i=1,commPatternCell_2nd(level)%nProcSend

        ! Store the processor id and the size of the message
        ! a bit easier.

        procID = commPatternCell_2nd(level)%sendProc(i)
        ssize   = commPatternCell_2nd(level)%nSend(i)

        ! Copy the data in the correct part of the send buffer.

        jj = ii
        do j=1,commPatternCell_2nd(level)%nSend(i)

           ! Store the block id and the indices of the donor
           ! a bit easier.

           d1 = commPatternCell_2nd(level)%sendList(i)%block(j)
           i1 = commPatternCell_2nd(level)%sendList(i)%indices(j,1)
           j1 = commPatternCell_2nd(level)%sendList(i)%indices(j,2)
           k1 = commPatternCell_2nd(level)%sendList(i)%indices(j,3)

           ! Copy the coordinates of this point in the buffer.
           ! Update the counter jj accordingly.

           intSendBuffer(jj)   = flowDoms(d1,level,mm)%globalCell(i1,j1,k1)
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
     receives: do i=1,commPatternCell_2nd(level)%nProcRecv

        ! Store the processor id and the size of the message
        ! a bit easier.

        procID = commPatternCell_2nd(level)%recvProc(i)
        ssize   = commPatternCell_2nd(level)%nRecv(i)

        ! Post the receive.

        call mpi_irecv(intRecvBuffer(ii), ssize, sumb_integer, procID, &
             myID, SUmb_comm_world, recvRequests(i), ierr)

        ! And update ii.

        ii = ii + ssize

     enddo receives

     ! Copy the local data.

     localCopy: do i=1,internalCell_2nd(level)%nCopy

        ! Store the block and the indices of the donor a bit easier.

        d1 = internalCell_2nd(level)%donorBlock(i)
        i1 = internalCell_2nd(level)%donorIndices(i,1)
        j1 = internalCell_2nd(level)%donorIndices(i,2)
        k1 = internalCell_2nd(level)%donorIndices(i,3)
        ! Idem for the halo's.

        d2 = internalCell_2nd(level)%haloBlock(i)
        i2 = internalCell_2nd(level)%haloIndices(i,1)
        j2 = internalCell_2nd(level)%haloIndices(i,2)
        k2 = internalCell_2nd(level)%haloIndices(i,3)
        ! Copy the coordinates.
        flowDoms(d2,level,mm)%globalCell(i2,j2,k2) = &
             flowDoms(d1,level,mm)%globalCell(i1,j1,k1)

     enddo localCopy

     ! Complete the nonblocking receives in an arbitrary sequence and
     ! copy the coordinates from the buffer into the halo's.

     ssize = commPatternCell_2nd(level)%nProcRecv
     completeRecvs: do i=1,commPatternCell_2nd(level)%nProcRecv

        ! Complete any of the requests.

        call mpi_waitany(ssize, recvRequests, index, status, ierr)

        ! Copy the data just arrived in the halo's.

        ii = index
        jj = commPatternCell_2nd(level)%nRecvCum(ii-1) +1
        do j=1,commPatternCell_2nd(level)%nRecv(ii)

           ! Store the block and the indices of the halo a bit easier.

           d2 = commPatternCell_2nd(level)%recvList(ii)%block(j)
           i2 = commPatternCell_2nd(level)%recvList(ii)%indices(j,1)
           j2 = commPatternCell_2nd(level)%recvList(ii)%indices(j,2)
           k2 = commPatternCell_2nd(level)%recvList(ii)%indices(j,3)

           ! Copy the data.

           flowDoms(d2,level,mm)%globalCell(i2,j2,k2) = intRecvBuffer(jj)
           jj = jj + 1

        enddo

     enddo completeRecvs

     ! Complete the nonblocking sends.

     ssize = commPatternCell_2nd(level)%nProcSend
     do i=1,commPatternCell_2nd(level)%nProcSend
        call mpi_waitany(ssize, sendRequests, index, status, ierr)
     enddo

  enddo spectralLoop
  deallocate(intSendBuffer, intRecvBuffer)
end subroutine exchangeCellNumbering

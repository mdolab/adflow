subroutine exchangeFringes(level, sps, commPattern, internal)
  !
  !       ExchangeFringes exchanges the donorInformation of the fringes: 
  !       donorProc, donorBlock, dI, dJ, dK and donorFrac. It does this  
  !       the 1:1 halos for the given level and spectral instance. Since 
  !       we have real values and integer values we will do all the ints 
  !       first and then the reals.                                      
  !
  use constants
  use blockPointers, only : flowDoms
  use communication, only : commType, internalCommType, recvBuffer, sendBuffer, myid, &
       sumb_comm_world, sendRequests, recvRequests
  !use overset
  implicit none
  !
  !      Subroutine arguments.
  !
  integer(kind=intType), intent(in) :: level, sps

  type(commType),          dimension(*), intent(in) :: commPattern
  type(internalCommType), dimension(*), intent(in) :: internal
  !
  !      Local variables.
  !
  integer :: size, procId, ierr, index
  integer, dimension(mpi_status_size) :: status

  integer(kind=intType) :: i, j, ii, jj, nVar
  integer(kind=intType) :: d1, i1, j1, k1, d2, i2, j2, k2

  integer(kind=intType), dimension(:), allocatable :: sendBufInt
  integer(kind=intType), dimension(:), allocatable :: recvBufInt

  ! Allocate the memory for the sending and receiving buffers.
  nVar = 5
  ii = commPattern(level)%nProcSend
  ii = commPattern(level)%nsendCum(ii)
  jj = commPattern(level)%nProcRecv
  jj = commPattern(level)%nrecvCum(jj)

  allocate(sendBufInt(ii*nVar), recvBufInt(jj*nVar), stat=ierr)

  ! Send the variables. The data is first copied into
  ! the send buffer after which the buffer is sent asap.

  ii = 1
  intSends: do i=1,commPattern(level)%nProcSend

     ! Store the processor id and the size of the message
     ! a bit easier.

     procID = commPattern(level)%sendProc(i)
     size   = nVar*commPattern(level)%nsend(i)

     ! Copy the data in the correct part of the send buffer.

     jj = ii
     do j=1,commPattern(level)%nsend(i)

        ! Store the block id and the indices of the donor
        ! a bit easier.

        d1 = commPattern(level)%sendList(i)%block(j)
        i1 = commPattern(level)%sendList(i)%indices(j,1)
        j1 = commPattern(level)%sendList(i)%indices(j,2)
        k1 = commPattern(level)%sendList(i)%indices(j,3)

        ! Copy integer values to buffer

        sendBufInt(jj  ) = flowDoms(d1,level,sps)%fringes(i1,j1,k1)%donorProc
        sendBufInt(jj+1) = flowDoms(d1,level,sps)%fringes(i1,j1,k1)%donorBlock
        sendBufInt(jj+2) = flowDoms(d1,level,sps)%fringes(i1,j1,k1)%dI
        sendBufInt(jj+3) = flowDoms(d1,level,sps)%fringes(i1,j1,k1)%dJ
        sendBufInt(jj+4) = flowDoms(d1,level,sps)%fringes(i1,j1,k1)%dK

        jj = jj + 5

     enddo

     ! Send the data.

     call mpi_isend(sendBufInt(ii), size, sumb_integer, procId, &
          procId, SUmb_comm_world, sendRequests(i),   &
          ierr)

     ! Set ii to jj for the next processor.

     ii = jj

  enddo intSends

  ! Post the nonblocking receives.

  ii = 1
  intReceives: do i=1,commPattern(level)%nProcRecv

     ! Store the processor id and the size of the message
     ! a bit easier.

     procID = commPattern(level)%recvProc(i)
     size    = nVar*commPattern(level)%nrecv(i)

     ! Post the receive.

     call mpi_irecv(recvBufInt(ii), size, sumb_integer, procId, &
          myId, SUmb_comm_world, recvRequests(i), ierr)

     ! And update ii.

     ii = ii + size

  enddo intReceives

  ! Copy the local data.

  intLocalCopy: do i=1,internal(level)%ncopy

     ! Store the block and the indices of the donor a bit easier.

     d1 = internal(level)%donorBlock(i)
     i1 = internal(level)%donorIndices(i,1)
     j1 = internal(level)%donorIndices(i,2)
     k1 = internal(level)%donorIndices(i,3)

     ! Idem for the halo's.

     d2 = internal(level)%haloBlock(i)
     i2 = internal(level)%haloIndices(i,1)
     j2 = internal(level)%haloIndices(i,2)
     k2 = internal(level)%haloIndices(i,3)

     flowDoms(d2,level,sps)%fringes(i2,j2,k2)%donorProc = &
          flowDoms(d1,level,sps)%fringes(i1,j1,k1)%donorProc

     flowDoms(d2,level,sps)%fringes(i2,j2,k2)%donorBlock = &
          flowDoms(d1,level,sps)%fringes(i1,j1,k1)%donorBlock

     flowDoms(d2,level,sps)%fringes(i2,j2,k2)%dI = &
          flowDoms(d1,level,sps)%fringes(i1,j1,k1)%dI

     flowDoms(d2,level,sps)%fringes(i2,j2,k2)%dJ = &
          flowDoms(d1,level,sps)%fringes(i1,j1,k1)%dJ

     flowDoms(d2,level,sps)%fringes(i2,j2,k2)%dK = &
          flowDoms(d1,level,sps)%fringes(i1,j1,k1)%dK
  enddo intLocalCopy

  ! Complete the nonblocking receives in an arbitrary sequence and
  ! copy the variables from the buffer into the halo's.

  size = commPattern(level)%nProcRecv
  intCompleteRecvs: do i=1,commPattern(level)%nProcRecv
     
     ! Complete any of the requests.

     call mpi_waitany(size, recvRequests, index, status, ierr)

     ! Copy the data just arrived in the halo's.

     ii = index
     jj = nVar*commPattern(level)%nrecvCum(ii-1)
     do j=1,commPattern(level)%nrecv(ii)

        ! Store the block and the indices of the halo a bit easier.

        d2 = commPattern(level)%recvList(ii)%block(j)
        i2 = commPattern(level)%recvList(ii)%indices(j,1)
        j2 = commPattern(level)%recvList(ii)%indices(j,2)
        k2 = commPattern(level)%recvList(ii)%indices(j,3)

        flowDoms(d2,level,sps)%fringes(i2,j2,k2)%donorProc  = recvBufInt(jj+1)
        flowDoms(d2,level,sps)%fringes(i2,j2,k2)%donorBlock = recvBufInt(jj+2)
        flowDoms(d2,level,sps)%fringes(i2,j2,k2)%dI = recvBufInt(jj+3)
        flowDoms(d2,level,sps)%fringes(i2,j2,k2)%dJ = recvBufInt(jj+4)
        flowDoms(d2,level,sps)%fringes(i2,j2,k2)%dK = recvBufInt(jj+5)

        jj = jj + 5
     enddo

  enddo intCompleteRecvs

  ! Complete the nonblocking sends.

  size = commPattern(level)%nProcSend
  do i=1,commPattern(level)%nProcSend
     call mpi_waitany(size, sendRequests, index, status, ierr)
  enddo

  ! Done with the integer memory. 

  deallocate(sendBufInt, recvBufInt)

  ! Now do the real exchange. We can use the regular real buffers here
  ! since they are large enough


  ! ================================================================================


 ! Allocate the memory for the sending and receiving buffers.
  nVar = 3

  ! Send the variables. The data is first copied into
  ! the send buffer after which the buffer is sent asap.

  ii = 1
  sends: do i=1,commPattern(level)%nProcSend

     ! Store the processor id and the size of the message
     ! a bit easier.

     procID = commPattern(level)%sendProc(i)
     size   = nVar*commPattern(level)%nsend(i)

     ! Copy the data in the correct part of the send buffer.

     jj = ii
     do j=1,commPattern(level)%nsend(i)

        ! Store the block id and the indices of the donor
        ! a bit easier.

        d1 = commPattern(level)%sendList(i)%block(j)
        i1 = commPattern(level)%sendList(i)%indices(j,1)
        j1 = commPattern(level)%sendList(i)%indices(j,2)
        k1 = commPattern(level)%sendList(i)%indices(j,3)

        ! Copy integer values to buffer

        sendBuffer(jj:jj+2) = flowDoms(d1,level,sps)%fringes(i1,j1,k1)%donorFrac

        jj = jj + 3

     enddo

     ! Send the data.

     call mpi_isend(sendBuffer(ii), size, sumb_real, procId, &
          procId, SUmb_comm_world, sendRequests(i),   &
          ierr)

     ! Set ii to jj for the next processor.

     ii = jj

  enddo sends

  ! Post the nonblocking receives.

  ii = 1
  receives: do i=1,commPattern(level)%nProcRecv

     ! Store the processor id and the size of the message
     ! a bit easier.

     procID = commPattern(level)%recvProc(i)
     size    = nVar*commPattern(level)%nrecv(i)

     ! Post the receive.

     call mpi_irecv(recvBuffer(ii), size, sumb_real, procId, &
          myId, SUmb_comm_world, recvRequests(i), ierr)

     ! And update ii.

     ii = ii + size

  enddo receives

  ! Copy the local data.

  localCopy: do i=1,internal(level)%ncopy

     ! Store the block and the indices of the donor a bit easier.

     d1 = internal(level)%donorBlock(i)
     i1 = internal(level)%donorIndices(i,1)
     j1 = internal(level)%donorIndices(i,2)
     k1 = internal(level)%donorIndices(i,3)

     ! Idem for the halo's.

     d2 = internal(level)%haloBlock(i)
     i2 = internal(level)%haloIndices(i,1)
     j2 = internal(level)%haloIndices(i,2)
     k2 = internal(level)%haloIndices(i,3)

     flowDoms(d2,level,sps)%fringes(i2,j2,k2)%donorFrac = &
          flowDoms(d1,level,sps)%fringes(i1,j1,k1)%donorFrac

  enddo localCopy

  ! Complete the nonblocking receives in an arbitrary sequence and
  ! copy the variables from the buffer into the halo's.

  size = commPattern(level)%nProcRecv
  completeRecvs: do i=1,commPattern(level)%nProcRecv

     ! Complete any of the requests.

     call mpi_waitany(size, recvRequests, index, status, ierr)

     ! Copy the data just arrived in the halo's.

     ii = index
     jj = nVar*commPattern(level)%nrecvCum(ii-1)
     do j=1,commPattern(level)%nrecv(ii)

        ! Store the block and the indices of the halo a bit easier.

        d2 = commPattern(level)%recvList(ii)%block(j)
        i2 = commPattern(level)%recvList(ii)%indices(j,1)
        j2 = commPattern(level)%recvList(ii)%indices(j,2)
        k2 = commPattern(level)%recvList(ii)%indices(j,3)

        flowDoms(d2,level,sps)%fringes(i2,j2,k2)%donorFrac = recvBuffer(jj+1:jj+3)

        jj = jj + 3
     enddo

  enddo completeRecvs

  ! Complete the nonblocking sends.

  size = commPattern(level)%nProcSend
  do i=1,commPattern(level)%nProcSend
     call mpi_waitany(size, sendRequests, index, status, ierr)
  enddo

end subroutine exchangeFringes

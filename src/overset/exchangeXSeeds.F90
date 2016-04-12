subroutine exchangeXSeeds(level, sps, commPattern, internal)

  use block
  use communication
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

  nVar = 3
  ii = commPattern(level)%nProcSend
  ii = commPattern(level)%nsendCum(ii)
  jj = commPattern(level)%nProcRecv
  jj = commPattern(level)%nrecvCum(jj)

  ! Send the variables. The data is first copied into
  ! the send buffer after which the buffer is sent asap.

  ii = 1
  sends: do i=1,commPattern(level)%nProcSend

     ! Store the processor id and the size of the message
     ! a bit easier.

     procID = commPattern(level)%sendProc(i)
     size    = nVar*commPattern(level)%nsend(i)

     ! Copy the data in the correct part of the send buffer.

     jj = ii
     do j=1,commPattern(level)%nsend(i)

        ! Store the block id and the indices of the donor
        ! a bit easier.

        d1 = commPattern(level)%sendList(i)%block(j)
        i1 = commPattern(level)%sendList(i)%indices(j,1)
        j1 = commPattern(level)%sendList(i)%indices(j,2)
        k1 = commPattern(level)%sendList(i)%indices(j,3)

        ! Copy xSeed values to buffer

        sendBuffer(jj    ) = flowDoms(d1,level,sps)%xSeed(i1, j1, k1, 1)
        sendBuffer(jj + 1) = flowDoms(d1,level,sps)%xSeed(i1, j1, k1, 2)
        sendBuffer(jj + 2) = flowDoms(d1,level,sps)%xSeed(i1, j1, k1, 3)
        jj = jj + 3

     enddo

     ! Send the data.

     call mpi_isend(sendBuffer(ii), size, sumb_integer, procId, &
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

     call mpi_irecv(recvBuffer(ii), size, sumb_integer, procId, &
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

     ! Copy the 3 values directly
     flowDoms(d2,level,sps)%xSeed(i2,j2,k2,:) = flowDoms(d1,level,sps)%xSeed(i1,j1,k1,:)

  enddo localCopy

  ! Complete the nonblocking receives in an arbitrary sequence and
  ! copy the variables from the buffer into the halo's.

  size = commPattern(level)%nProcRecv
  completeRecvs: do i=1,commPattern(level)%nProcRecv

     ! Complete any of the requests.

     call mpi_waitany(size, recvRequests, index, status, ierr)

     ! Copy the data just arrived in the halo's.

     ii = index
     jj = commPattern(level)%nrecvCum(ii-1)
     do j=1,commPattern(level)%nrecv(ii)

        ! Store the block and the indices of the halo a bit easier.

        d2 = commPattern(level)%recvList(ii)%block(j)
        i2 = commPattern(level)%recvList(ii)%indices(j,1)
        j2 = commPattern(level)%recvList(ii)%indices(j,2)
        k2 = commPattern(level)%recvList(ii)%indices(j,3)

        ! Copy the minimum of the iblank value and 9 (see above).

        jj = jj + 1
        flowDoms(d2,level,sps)%xSeed(i2,j2,k2,1) = recvBuffer(jj)
        jj = jj + 1
        flowDoms(d2,level,sps)%xSeed(i2,j2,k2,2) = recvBuffer(jj)
        jj = jj + 1
        flowDoms(d2,level,sps)%xSeed(i2,j2,k2,3) = recvBuffer(jj)

     enddo

  enddo completeRecvs

  ! Complete the nonblocking sends.

  size = commPattern(level)%nProcSend
  do i=1,commPattern(level)%nProcSend
     call mpi_waitany(size, sendRequests, index, status, ierr)
  enddo

end subroutine exchangeXSeeds

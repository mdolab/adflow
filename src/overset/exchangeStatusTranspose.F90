subroutine exchangeStatusTranspose(level, sps, commPattern, internal)

  ! exchangeStatusTranspose performs the *TRANSPOSE* of the normal
  ! halo exchange. That means it takes information *in the halo cells*
  ! and accumulate it into the *owned cells*. In this particular case,
  ! we are transmitting the isDonor and isWallDonor information from
  ! the halos to the owned cells. The "accumulate" operation will be
  ! an MPI_LOR.  Note that this actually hast he same comm structure
  ! as 'whalo1to1_b'.

  use block
  use communication
  use inputTimeSpectral
  implicit none
  !
  !      Subroutine arguments.
  !
  integer(kind=intType), intent(in) :: level, sps
  type(commType), dimension(*), intent(in)         :: commPattern
  type(internalCommType), dimension(*), intent(in) :: internal
  !
  !      Local variables.
  !
  integer :: size, procID, ierr, index
  integer, dimension(mpi_status_size) :: status

  integer(kind=intType) :: mm
  integer(kind=intType) :: i, j, k, ii, jj
  integer(kind=intType) :: d1, i1, j1, k1, d2, i2, j2, k2
  integer(kind=intType), dimension(:), allocatable :: sendBuf, recvBuf
  logical :: CisDonor, CisHole, CisCompute, CisFloodSeed, CisFlooded, CisWall, CisWallDonor
  logical :: DisDonor, DisHole, DisCompute, DisFloodSeed, DisFlooded, DisWall, DisWallDonor
integer(kind=intType) :: cellStatus, donorStatus

  !
  !      ******************************************************************
  !      *                                                                *
  !      * Begin execution                                                *
  !      *                                                                *
  !      ******************************************************************
  !

  ii = commPattern(level)%nProcSend
  ii = commPattern(level)%nsendCum(ii)
  jj = commPattern(level)%nProcRecv
  jj = commPattern(level)%nrecvCum(jj)

  ! We are exchanging 1 piece of information
  allocate(sendBuf(ii), recvBuf(jj), stat=ierr)

  ! Gather up the seeds into the *recv* buffer. Note we loop
  ! over nProcRECV here! After the buffer is assembled it is
  ! sent off.

  jj = 1
  ii = 1
  recvs: do i=1,commPattern(level)%nProcRecv

     ! Store the processor id and the size of the message
     ! a bit easier.

     procID = commPattern(level)%recvProc(i)
     size    = commPattern(level)%nrecv(i)

     ! Copy the data into the buffer

     do j=1,commPattern(level)%nrecv(i)

        ! Store the block and the indices of the halo a bit easier.

        d2 = commPattern(level)%recvList(i)%block(j)
        i2 = commPattern(level)%recvList(i)%indices(j,1)
        j2 = commPattern(level)%recvList(i)%indices(j,2)
        k2 = commPattern(level)%recvList(i)%indices(j,3)

        recvBuf(jj) = flowDoms(d2, level, sps)%fringes(i2, j2, k2)%status
        jj = jj + 1

     enddo

     ! Send the data.
     call mpi_isend(recvBuf(ii), size, sumb_integer, procID,  &
          procID, SUmb_comm_world, sendRequests(i), &
          ierr)

     ! Set ii to jj for the next processor.

     ii = jj

  enddo recvs

  ! Post the nonblocking receives.

  ii = 1
  sends: do i=1,commPattern(level)%nProcSend

     ! Store the processor id and the size of the message
     ! a bit easier.

     procID = commPattern(level)%sendProc(i)
     size    = commPattern(level)%nsend(i)

     ! Post the receive.

     call mpi_irecv(sendBuf(ii), size, sumb_integer, procID, &
          myID, SUmb_comm_world, recvRequests(i), ierr)

     ! And update ii.

     ii = ii + size

  enddo sends

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

     ! OR operation. Note we modify the '1' values ie. the 'donors'
     ! which are now receivers because of the transpose operation.
     cellStatus = flowDoms(d1, level, sps)%fringes(i1, j1, k1)%status
     call getStatus(cellStatus, CisDonor, CisHole, CisCompute, &
          CisFloodSeed, CisFlooded, CisWall, CisWallDonor)

     donorStatus = flowDoms(d2, level, sps)%fringes(i2, j2, k2)%status
     call getStatus(donorStatus, DisDonor, DisHole, DisCompute, &
          DisFloodSeed, DisFlooded, DisWall, DisWallDonor)
     
     call setIsDonor(flowDoms(d1, level, sps)%fringes(i1, j1, k1)%status, &
          CIsDonor .or. DisDonor)

     call setIsWallDonor(flowDoms(d1, level, sps)%fringes(i1, j1, k1)%status, &
          CIsWallDonor .or. DisWallDonor)

  enddo localCopy

  ! Complete the nonblocking receives in an arbitrary sequence and
  ! copy the variables from the buffer into the halo's.

  size = commPattern(level)%nProcSend
  completeSends: do i=1,commPattern(level)%nProcSend

     ! Complete any of the requests.

     call mpi_waitany(size, recvRequests, index, status, ierr)

     ! Copy the data just arrived in the halo's.

     ii = index

     jj = commPattern(level)%nsendCum(ii-1)

     do j=1,commPattern(level)%nsend(ii)

        ! Store the block and the indices of the halo a bit easier.

        d1 = commPattern(level)%sendList(ii)%block(j)
        i1 = commPattern(level)%sendList(ii)%indices(j,1)
        j1 = commPattern(level)%sendList(ii)%indices(j,2)
        k1 = commPattern(level)%sendList(ii)%indices(j,3)


        cellStatus = flowDoms(d1, level, sps)%fringes(i1, j1, k1)%status
        call getStatus(cellStatus, CisDonor, CisHole, CisCompute, &
             CisFloodSeed, CisFlooded, CisWall, CisWallDonor)
        jj = jj + 1
        donorStatus = sendBuf(jj)
        call getStatus(donorStatus, DisDonor, DisHole, DisCompute, &
             DisFloodSeed, DisFlooded, DisWall, DisWallDonor)
        
        call setIsDonor(flowDoms(d1, level, sps)%fringes(i1, j1, k1)%status, &
             CIsDonor .or. DisDonor)
        
        call setIsWallDonor(flowDoms(d1, level, sps)%fringes(i1, j1, k1)%status, &
             CIsWallDonor .or. DisWallDonor)
     enddo
  enddo completeSends

  ! Complete the nonblocking sends.

  size = commPattern(level)%nProcRecv
  do i=1,commPattern(level)%nProcRecv
     call mpi_waitany(size, sendRequests, index, status, ierr)
  enddo

  deallocate(recvBuf, sendBuf)

end subroutine exchangeStatusTranspose

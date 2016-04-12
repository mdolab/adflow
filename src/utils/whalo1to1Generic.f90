subroutine whalo1to1RealGeneric(nVar, level, sps, commPattern, internal)
  !
  !      ******************************************************************
  !      *                                                                *
  !      * whalo1to1 exchanges the 1 to 1 internal halo's for the cell    *
  !      * centered variables for the given communication pattern.        *
  !      * Pointers must be set for var1, var2...varN                     *
  !      *                                                                *
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
  integer(kind=intType), intent(in) :: level, sps

  type(commType), dimension(*), intent(in)         :: commPattern
  type(internalCommType), dimension(*), intent(in) :: internal
  !
  !      Local variables.
  !

  integer :: size, procID, ierr, index
  integer, dimension(mpi_status_size) :: status

  integer(kind=intType) :: nVar, mm
  integer(kind=intType) :: i, j, k, ii, jj
  integer(kind=intType) :: d1, i1, j1, k1, d2, i2, j2, k2
  !
  !      ******************************************************************
  !      *                                                                *
  !      * Begin execution                                                *
  !      *                                                                *
  !      ******************************************************************
  !

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

        ! Copy the given range of the working variables for
        ! this cell in the buffer. Update the counter jj.
        if (nVar == 1) then 
           sendBuffer(jj) = flowDoms(d1,level,sps)%rVar1(i1,j1,k1)
           jj = jj + 1

        else if (nVar == 2) then 
           sendBuffer(jj  ) = flowDoms(d1,level,sps)%rVar1(i1,j1,k1)
           sendBuffer(jj+1) = flowDoms(d1,level,sps)%rVar2(i1,j1,k1)
           jj = jj + 2

        else if (nVar == 3) then 
           sendBuffer(jj  ) = flowDoms(d1,level,sps)%rVar1(i1,j1,k1)
           sendBuffer(jj+1) = flowDoms(d1,level,sps)%rVar2(i1,j1,k1)
           sendBuffer(jj+2) = flowDoms(d1,level,sps)%rVar3(i1,j1,k1)
           jj = jj + 3

        else if (nVar == 4) then 
           sendBuffer(jj  ) = flowDoms(d1,level,sps)%rVar1(i1,j1,k1)
           sendBuffer(jj+1) = flowDoms(d1,level,sps)%rVar2(i1,j1,k1)
           sendBuffer(jj+2) = flowDoms(d1,level,sps)%rVar3(i1,j1,k1)
           sendBuffer(jj+3) = flowDoms(d1,level,sps)%rVar4(i1,j1,k1)
           jj = jj + 4

        else if (nVar == 5) then 
           sendBuffer(jj  ) = flowDoms(d1,level,sps)%rVar1(i1,j1,k1)
           sendBuffer(jj+1) = flowDoms(d1,level,sps)%rVar2(i1,j1,k1)
           sendBuffer(jj+2) = flowDoms(d1,level,sps)%rVar3(i1,j1,k1)
           sendBuffer(jj+3) = flowDoms(d1,level,sps)%rVar4(i1,j1,k1)
           sendBuffer(jj+4) = flowDoms(d1,level,sps)%rVar5(i1,j1,k1)
           jj = jj + 5

        else if (nVar == 6) then 
           sendBuffer(jj  ) = flowDoms(d1,level,sps)%rVar1(i1,j1,k1)
           sendBuffer(jj+1) = flowDoms(d1,level,sps)%rVar2(i1,j1,k1)
           sendBuffer(jj+2) = flowDoms(d1,level,sps)%rVar3(i1,j1,k1)
           sendBuffer(jj+3) = flowDoms(d1,level,sps)%rVar4(i1,j1,k1)
           sendBuffer(jj+4) = flowDoms(d1,level,sps)%rVar5(i1,j1,k1)
           sendBuffer(jj+5) = flowDoms(d1,level,sps)%rVar6(i1,j1,k1)
           jj = jj + 6

        else if (nVar == 7) then 
           sendBuffer(jj  ) = flowDoms(d1,level,sps)%rVar1(i1,j1,k1)
           sendBuffer(jj+1) = flowDoms(d1,level,sps)%rVar2(i1,j1,k1)
           sendBuffer(jj+2) = flowDoms(d1,level,sps)%rVar3(i1,j1,k1)
           sendBuffer(jj+3) = flowDoms(d1,level,sps)%rVar4(i1,j1,k1)
           sendBuffer(jj+4) = flowDoms(d1,level,sps)%rVar5(i1,j1,k1)
           sendBuffer(jj+5) = flowDoms(d1,level,sps)%rVar6(i1,j1,k1)
           sendBuffer(jj+6) = flowDoms(d1,level,sps)%rVar7(i1,j1,k1)
           jj = jj + 7

        else if (nVar == 8) then 
           sendBuffer(jj  ) = flowDoms(d1,level,sps)%rVar1(i1,j1,k1)
           sendBuffer(jj+1) = flowDoms(d1,level,sps)%rVar2(i1,j1,k1)
           sendBuffer(jj+2) = flowDoms(d1,level,sps)%rVar3(i1,j1,k1)
           sendBuffer(jj+3) = flowDoms(d1,level,sps)%rVar4(i1,j1,k1)
           sendBuffer(jj+4) = flowDoms(d1,level,sps)%rVar5(i1,j1,k1)
           sendBuffer(jj+5) = flowDoms(d1,level,sps)%rVar6(i1,j1,k1)
           sendBuffer(jj+6) = flowDoms(d1,level,sps)%rVar7(i1,j1,k1)
           sendBuffer(jj+7) = flowDoms(d1,level,sps)%rVar8(i1,j1,k1)
           jj = jj + 8

        else if (nVar == 9) then 
           sendBuffer(jj  ) = flowDoms(d1,level,sps)%rVar1(i1,j1,k1)
           sendBuffer(jj+1) = flowDoms(d1,level,sps)%rVar2(i1,j1,k1)
           sendBuffer(jj+2) = flowDoms(d1,level,sps)%rVar3(i1,j1,k1)
           sendBuffer(jj+3) = flowDoms(d1,level,sps)%rVar4(i1,j1,k1)
           sendBuffer(jj+4) = flowDoms(d1,level,sps)%rVar5(i1,j1,k1)
           sendBuffer(jj+5) = flowDoms(d1,level,sps)%rVar6(i1,j1,k1)
           sendBuffer(jj+6) = flowDoms(d1,level,sps)%rVar7(i1,j1,k1)
           sendBuffer(jj+7) = flowDoms(d1,level,sps)%rVar8(i1,j1,k1)
           sendBuffer(jj+8) = flowDoms(d1,level,sps)%rVar9(i1,j1,k1)
           jj = jj + 9

        end if
     end do

     ! Send the data.

     call mpi_isend(sendBuffer(ii), size, sumb_real, procID,  &
          procID, SUmb_comm_world, sendRequests(i), &
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

     call mpi_irecv(recvBuffer(ii), size, sumb_real, procID, &
          myID, SUmb_comm_world, recvRequests(i), ierr)

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

     if (nVar == 1) then 
        flowDoms(d2,level,sps)%rVar1(i2, j2, k2) = flowDoms(d1,level,sps)%rVar1(i1,j1,k1)
        
     else if (nVar == 2) then 
        flowDoms(d2,level,sps)%rVar1(i2, j2, k2) = flowDoms(d1,level,sps)%rVar1(i1,j1,k1)
        flowDoms(d2,level,sps)%rVar2(i2, j2, k2) = flowDoms(d1,level,sps)%rVar2(i1,j1,k1)

     else if (nVar == 3) then 
        flowDoms(d2,level,sps)%rVar1(i2, j2, k2) = flowDoms(d1,level,sps)%rVar1(i1,j1,k1)
        flowDoms(d2,level,sps)%rVar2(i2, j2, k2) = flowDoms(d1,level,sps)%rVar2(i1,j1,k1)
        flowDoms(d2,level,sps)%rVar3(i2, j2, k2) = flowDoms(d1,level,sps)%rVar3(i1,j1,k1)

     else if (nVar == 4) then 
        flowDoms(d2,level,sps)%rVar1(i2, j2, k2) = flowDoms(d1,level,sps)%rVar1(i1,j1,k1)
        flowDoms(d2,level,sps)%rVar2(i2, j2, k2) = flowDoms(d1,level,sps)%rVar2(i1,j1,k1)
        flowDoms(d2,level,sps)%rVar3(i2, j2, k2) = flowDoms(d1,level,sps)%rVar3(i1,j1,k1)
        flowDoms(d2,level,sps)%rVar4(i2, j2, k2) = flowDoms(d1,level,sps)%rVar4(i1,j1,k1)

     else if (nVar == 5) then 
        flowDoms(d2,level,sps)%rVar1(i2, j2, k2) = flowDoms(d1,level,sps)%rVar1(i1,j1,k1)
        flowDoms(d2,level,sps)%rVar2(i2, j2, k2) = flowDoms(d1,level,sps)%rVar2(i1,j1,k1)
        flowDoms(d2,level,sps)%rVar3(i2, j2, k2) = flowDoms(d1,level,sps)%rVar3(i1,j1,k1)
        flowDoms(d2,level,sps)%rVar4(i2, j2, k2) = flowDoms(d1,level,sps)%rVar4(i1,j1,k1)
        flowDoms(d2,level,sps)%rVar5(i2, j2, k2) = flowDoms(d1,level,sps)%rVar5(i1,j1,k1)

     else if (nVar == 6) then 
        flowDoms(d2,level,sps)%rVar1(i2, j2, k2) = flowDoms(d1,level,sps)%rVar1(i1,j1,k1)
        flowDoms(d2,level,sps)%rVar2(i2, j2, k2) = flowDoms(d1,level,sps)%rVar2(i1,j1,k1)
        flowDoms(d2,level,sps)%rVar3(i2, j2, k2) = flowDoms(d1,level,sps)%rVar3(i1,j1,k1)
        flowDoms(d2,level,sps)%rVar4(i2, j2, k2) = flowDoms(d1,level,sps)%rVar4(i1,j1,k1)
        flowDoms(d2,level,sps)%rVar5(i2, j2, k2) = flowDoms(d1,level,sps)%rVar5(i1,j1,k1)
        flowDoms(d2,level,sps)%rVar6(i2, j2, k2) = flowDoms(d1,level,sps)%rVar6(i1,j1,k1)

     else if (nVar == 7) then 
        flowDoms(d2,level,sps)%rVar1(i2, j2, k2) = flowDoms(d1,level,sps)%rVar1(i1,j1,k1)
        flowDoms(d2,level,sps)%rVar2(i2, j2, k2) = flowDoms(d1,level,sps)%rVar2(i1,j1,k1)
        flowDoms(d2,level,sps)%rVar3(i2, j2, k2) = flowDoms(d1,level,sps)%rVar3(i1,j1,k1)
        flowDoms(d2,level,sps)%rVar4(i2, j2, k2) = flowDoms(d1,level,sps)%rVar4(i1,j1,k1)
        flowDoms(d2,level,sps)%rVar5(i2, j2, k2) = flowDoms(d1,level,sps)%rVar5(i1,j1,k1)
        flowDoms(d2,level,sps)%rVar6(i2, j2, k2) = flowDoms(d1,level,sps)%rVar6(i1,j1,k1)
        flowDoms(d2,level,sps)%rVar7(i2, j2, k2) = flowDoms(d1,level,sps)%rVar7(i1,j1,k1)

     else if (nVar == 8) then 
        flowDoms(d2,level,sps)%rVar1(i2, j2, k2) = flowDoms(d1,level,sps)%rVar1(i1,j1,k1)
        flowDoms(d2,level,sps)%rVar2(i2, j2, k2) = flowDoms(d1,level,sps)%rVar2(i1,j1,k1)
        flowDoms(d2,level,sps)%rVar3(i2, j2, k2) = flowDoms(d1,level,sps)%rVar3(i1,j1,k1)
        flowDoms(d2,level,sps)%rVar4(i2, j2, k2) = flowDoms(d1,level,sps)%rVar4(i1,j1,k1)
        flowDoms(d2,level,sps)%rVar5(i2, j2, k2) = flowDoms(d1,level,sps)%rVar5(i1,j1,k1)
        flowDoms(d2,level,sps)%rVar6(i2, j2, k2) = flowDoms(d1,level,sps)%rVar6(i1,j1,k1)
        flowDoms(d2,level,sps)%rVar7(i2, j2, k2) = flowDoms(d1,level,sps)%rVar7(i1,j1,k1)
        flowDoms(d2,level,sps)%rVar8(i2, j2, k2) = flowDoms(d1,level,sps)%rVar8(i1,j1,k1)

     else if (nVar == 9) then 
        flowDoms(d2,level,sps)%rVar1(i2, j2, k2) = flowDoms(d1,level,sps)%rVar1(i1,j1,k1)
        flowDoms(d2,level,sps)%rVar2(i2, j2, k2) = flowDoms(d1,level,sps)%rVar2(i1,j1,k1)
        flowDoms(d2,level,sps)%rVar3(i2, j2, k2) = flowDoms(d1,level,sps)%rVar3(i1,j1,k1)
        flowDoms(d2,level,sps)%rVar4(i2, j2, k2) = flowDoms(d1,level,sps)%rVar4(i1,j1,k1)
        flowDoms(d2,level,sps)%rVar5(i2, j2, k2) = flowDoms(d1,level,sps)%rVar5(i1,j1,k1)
        flowDoms(d2,level,sps)%rVar6(i2, j2, k2) = flowDoms(d1,level,sps)%rVar6(i1,j1,k1)
        flowDoms(d2,level,sps)%rVar7(i2, j2, k2) = flowDoms(d1,level,sps)%rVar7(i1,j1,k1)
        flowDoms(d2,level,sps)%rVar8(i2, j2, k2) = flowDoms(d1,level,sps)%rVar8(i1,j1,k1)
        flowDoms(d2,level,sps)%rVar9(i2, j2, k2) = flowDoms(d1,level,sps)%rVar9(i1,j1,k1)
     end if

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

        if (nVar == 1) then 
           jj = jj + 1
           flowDoms(d2,level,sps)%rVar1(i2, j2, k2) = recvBuffer(jj)

        else if (nVar == 2) then 
           jj = jj + 1
           flowDoms(d2,level,sps)%rVar1(i2, j2, k2) = recvBuffer(jj)
           jj = jj + 1
           flowDoms(d2,level,sps)%rVar2(i2, j2, k2) = recvBuffer(jj)

        else if (nVar == 3) then 
           jj = jj + 1
           flowDoms(d2,level,sps)%rVar1(i2, j2, k2) = recvBuffer(jj)
           jj = jj + 1
           flowDoms(d2,level,sps)%rVar2(i2, j2, k2) = recvBuffer(jj)
           jj = jj + 1
           flowDoms(d2,level,sps)%rVar3(i2, j2, k2) = recvBuffer(jj)

        else if (nVar == 4) then 
           jj = jj + 1
           flowDoms(d2,level,sps)%rVar1(i2, j2, k2) = recvBuffer(jj)
           jj = jj + 1
           flowDoms(d2,level,sps)%rVar2(i2, j2, k2) = recvBuffer(jj)
           jj = jj + 1
           flowDoms(d2,level,sps)%rVar3(i2, j2, k2) = recvBuffer(jj)
           jj = jj + 1
           flowDoms(d2,level,sps)%rVar4(i2, j2, k2) = recvBuffer(jj)

        else if (nVar == 5) then 
           jj = jj + 1
           flowDoms(d2,level,sps)%rVar1(i2, j2, k2) = recvBuffer(jj)
           jj = jj + 1
           flowDoms(d2,level,sps)%rVar2(i2, j2, k2) = recvBuffer(jj)
           jj = jj + 1
           flowDoms(d2,level,sps)%rVar3(i2, j2, k2) = recvBuffer(jj)
           jj = jj + 1
           flowDoms(d2,level,sps)%rVar4(i2, j2, k2) = recvBuffer(jj)
           jj = jj + 1
           flowDoms(d2,level,sps)%rVar5(i2, j2, k2) = recvBuffer(jj)

        else if (nVar == 6) then 
           jj = jj + 1
           flowDoms(d2,level,sps)%rVar1(i2, j2, k2) = recvBuffer(jj)
           jj = jj + 1
           flowDoms(d2,level,sps)%rVar2(i2, j2, k2) = recvBuffer(jj)
           jj = jj + 1
           flowDoms(d2,level,sps)%rVar3(i2, j2, k2) = recvBuffer(jj)
           jj = jj + 1
           flowDoms(d2,level,sps)%rVar4(i2, j2, k2) = recvBuffer(jj)
           jj = jj + 1
           flowDoms(d2,level,sps)%rVar5(i2, j2, k2) = recvBuffer(jj)
           jj = jj + 1
           flowDoms(d2,level,sps)%rVar6(i2, j2, k2) = recvBuffer(jj)

        else if (nVar == 7) then 
           jj = jj + 1
           flowDoms(d2,level,sps)%rVar1(i2, j2, k2) = recvBuffer(jj)
           jj = jj + 1
           flowDoms(d2,level,sps)%rVar2(i2, j2, k2) = recvBuffer(jj)
           jj = jj + 1
           flowDoms(d2,level,sps)%rVar3(i2, j2, k2) = recvBuffer(jj)
           jj = jj + 1
           flowDoms(d2,level,sps)%rVar4(i2, j2, k2) = recvBuffer(jj)
           jj = jj + 1
           flowDoms(d2,level,sps)%rVar5(i2, j2, k2) = recvBuffer(jj)
           jj = jj + 1
           flowDoms(d2,level,sps)%rVar6(i2, j2, k2) = recvBuffer(jj)
           jj = jj + 1
           flowDoms(d2,level,sps)%rVar7(i2, j2, k2) = recvBuffer(jj)

        else if (nVar == 8) then 
           jj = jj + 1
           flowDoms(d2,level,sps)%rVar1(i2, j2, k2) = recvBuffer(jj)
           jj = jj + 1
           flowDoms(d2,level,sps)%rVar2(i2, j2, k2) = recvBuffer(jj)
           jj = jj + 1
           flowDoms(d2,level,sps)%rVar3(i2, j2, k2) = recvBuffer(jj)
           jj = jj + 1
           flowDoms(d2,level,sps)%rVar4(i2, j2, k2) = recvBuffer(jj)
           jj = jj + 1
           flowDoms(d2,level,sps)%rVar5(i2, j2, k2) = recvBuffer(jj)
           jj = jj + 1
           flowDoms(d2,level,sps)%rVar6(i2, j2, k2) = recvBuffer(jj)
           jj = jj + 1
           flowDoms(d2,level,sps)%rVar7(i2, j2, k2) = recvBuffer(jj)
           jj = jj + 1
           flowDoms(d2,level,sps)%rVar8(i2, j2, k2) = recvBuffer(jj)

        else if (nVar == 9) then 
           jj = jj + 1
           flowDoms(d2,level,sps)%rVar1(i2, j2, k2) = recvBuffer(jj)
           jj = jj + 1
           flowDoms(d2,level,sps)%rVar2(i2, j2, k2) = recvBuffer(jj)
           jj = jj + 1
           flowDoms(d2,level,sps)%rVar3(i2, j2, k2) = recvBuffer(jj)
           jj = jj + 1
           flowDoms(d2,level,sps)%rVar4(i2, j2, k2) = recvBuffer(jj)
           jj = jj + 1
           flowDoms(d2,level,sps)%rVar5(i2, j2, k2) = recvBuffer(jj)
           jj = jj + 1
           flowDoms(d2,level,sps)%rVar6(i2, j2, k2) = recvBuffer(jj)
           jj = jj + 1
           flowDoms(d2,level,sps)%rVar7(i2, j2, k2) = recvBuffer(jj)
           jj = jj + 1
           flowDoms(d2,level,sps)%rVar8(i2, j2, k2) = recvBuffer(jj)
           jj = jj + 1
           flowDoms(d2,level,sps)%rVar9(i2, j2, k2) = recvBuffer(jj)
        end if
     enddo

  enddo completeRecvs


  ! Complete the nonblocking sends.

  size = commPattern(level)%nProcSend
  do i=1,commPattern(level)%nProcSend
     call mpi_waitany(size, sendRequests, index, status, ierr)
  enddo


end subroutine whalo1to1RealGeneric

subroutine whalo1to1IntGeneric(nVar, level, sps, commPattern, internal)
  !
  !      ******************************************************************
  !      *                                                                *
  !      * whalo1to1 exchanges the 1 to 1 internal halo's for the cell    *
  !      * centered variables for the given communication pattern.        *
  !      * Pointers must be set for var1, var2...varN                     *
  !      *                                                                *
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
  integer(kind=intType), intent(in) :: level, sps

  type(commType), dimension(*), intent(in)         :: commPattern
  type(internalCommType), dimension(*), intent(in) :: internal

  !
  !      Local variables.
  !
  integer :: size, procID, ierr, index
  integer, dimension(mpi_status_size) :: status

  integer(kind=intType) :: nVar, mm
  integer(kind=intType) :: i, j, k, ii, jj
  integer(kind=intType) :: d1, i1, j1, k1, d2, i2, j2, k2

  integer(kind=intType), dimension(:), allocatable :: sendBufInt
  integer(kind=intType), dimension(:), allocatable :: recvBufInt

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

  allocate(sendBufInt(ii), recvBufInt(jj), stat=ierr)

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

        ! Copy the given range of the working variables for
        ! this cell in the buffer. Update the counter jj.
        if (nVar == 1) then 
           sendBufInt(jj) = flowDoms(d1,level,sps)%iVar1(i1,j1,k1)
           jj = jj + 1

        else if (nVar == 2) then 
           sendBufInt(jj  ) = flowDoms(d1,level,sps)%iVar1(i1,j1,k1)
           sendBufInt(jj+1) = flowDoms(d1,level,sps)%iVar2(i1,j1,k1)
           jj = jj + 2

        else if (nVar == 3) then 
           sendBufInt(jj  ) = flowDoms(d1,level,sps)%iVar1(i1,j1,k1)
           sendBufInt(jj+1) = flowDoms(d1,level,sps)%iVar2(i1,j1,k1)
           sendBufInt(jj+2) = flowDoms(d1,level,sps)%iVar3(i1,j1,k1)
           jj = jj + 3

        else if (nVar == 4) then 
           sendBufInt(jj  ) = flowDoms(d1,level,sps)%iVar1(i1,j1,k1)
           sendBufInt(jj+1) = flowDoms(d1,level,sps)%iVar2(i1,j1,k1)
           sendBufInt(jj+2) = flowDoms(d1,level,sps)%iVar3(i1,j1,k1)
           sendBufInt(jj+3) = flowDoms(d1,level,sps)%iVar4(i1,j1,k1)
           jj = jj + 4

        else if (nVar == 5) then 
           sendBufInt(jj  ) = flowDoms(d1,level,sps)%iVar1(i1,j1,k1)
           sendBufInt(jj+1) = flowDoms(d1,level,sps)%iVar2(i1,j1,k1)
           sendBufInt(jj+2) = flowDoms(d1,level,sps)%iVar3(i1,j1,k1)
           sendBufInt(jj+3) = flowDoms(d1,level,sps)%iVar4(i1,j1,k1)
           sendBufInt(jj+4) = flowDoms(d1,level,sps)%iVar5(i1,j1,k1)
           jj = jj + 5

        else if (nVar == 6) then 
           sendBufInt(jj  ) = flowDoms(d1,level,sps)%iVar1(i1,j1,k1)
           sendBufInt(jj+1) = flowDoms(d1,level,sps)%iVar2(i1,j1,k1)
           sendBufInt(jj+2) = flowDoms(d1,level,sps)%iVar3(i1,j1,k1)
           sendBufInt(jj+3) = flowDoms(d1,level,sps)%iVar4(i1,j1,k1)
           sendBufInt(jj+4) = flowDoms(d1,level,sps)%iVar5(i1,j1,k1)
           sendBufInt(jj+5) = flowDoms(d1,level,sps)%iVar6(i1,j1,k1)
           jj = jj + 6

        else if (nVar == 7) then 
           sendBufInt(jj  ) = flowDoms(d1,level,sps)%iVar1(i1,j1,k1)
           sendBufInt(jj+1) = flowDoms(d1,level,sps)%iVar2(i1,j1,k1)
           sendBufInt(jj+2) = flowDoms(d1,level,sps)%iVar3(i1,j1,k1)
           sendBufInt(jj+3) = flowDoms(d1,level,sps)%iVar4(i1,j1,k1)
           sendBufInt(jj+4) = flowDoms(d1,level,sps)%iVar5(i1,j1,k1)
           sendBufInt(jj+5) = flowDoms(d1,level,sps)%iVar6(i1,j1,k1)
           sendBufInt(jj+6) = flowDoms(d1,level,sps)%iVar7(i1,j1,k1)
           jj = jj + 7

        else if (nVar == 8) then 
           sendBufInt(jj  ) = flowDoms(d1,level,sps)%iVar1(i1,j1,k1)
           sendBufInt(jj+1) = flowDoms(d1,level,sps)%iVar2(i1,j1,k1)
           sendBufInt(jj+2) = flowDoms(d1,level,sps)%iVar3(i1,j1,k1)
           sendBufInt(jj+3) = flowDoms(d1,level,sps)%iVar4(i1,j1,k1)
           sendBufInt(jj+4) = flowDoms(d1,level,sps)%iVar5(i1,j1,k1)
           sendBufInt(jj+5) = flowDoms(d1,level,sps)%iVar6(i1,j1,k1)
           sendBufInt(jj+6) = flowDoms(d1,level,sps)%iVar7(i1,j1,k1)
           sendBufInt(jj+7) = flowDoms(d1,level,sps)%iVar8(i1,j1,k1)
           jj = jj + 8

        else if (nVar == 9) then 
           sendBufInt(jj  ) = flowDoms(d1,level,sps)%iVar1(i1,j1,k1)
           sendBufInt(jj+1) = flowDoms(d1,level,sps)%iVar2(i1,j1,k1)
           sendBufInt(jj+2) = flowDoms(d1,level,sps)%iVar3(i1,j1,k1)
           sendBufInt(jj+3) = flowDoms(d1,level,sps)%iVar4(i1,j1,k1)
           sendBufInt(jj+4) = flowDoms(d1,level,sps)%iVar5(i1,j1,k1)
           sendBufInt(jj+5) = flowDoms(d1,level,sps)%iVar6(i1,j1,k1)
           sendBufInt(jj+6) = flowDoms(d1,level,sps)%iVar7(i1,j1,k1)
           sendBufInt(jj+7) = flowDoms(d1,level,sps)%iVar8(i1,j1,k1)
           sendBufInt(jj+8) = flowDoms(d1,level,sps)%iVar9(i1,j1,k1)
           jj = jj + 9

        end if
     end do

     ! Send the data.

     call mpi_isend(sendBufInt(ii), size, sumb_real, procID,  &
          procID, SUmb_comm_world, sendRequests(i), &
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

     call mpi_irecv(recvBufInt(ii), size, sumb_real, procID, &
          myID, SUmb_comm_world, recvRequests(i), ierr)

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

     if (nVar == 1) then 
        flowDoms(d2,level,sps)%iVar1(i2, j2, k2) = flowDoms(d1,level,sps)%iVar1(i1,j1,k1)

     else if (nVar == 2) then 
        flowDoms(d2,level,sps)%iVar1(i2, j2, k2) = flowDoms(d1,level,sps)%iVar1(i1,j1,k1)
        flowDoms(d2,level,sps)%iVar2(i2, j2, k2) = flowDoms(d1,level,sps)%iVar2(i1,j1,k1)

     else if (nVar == 3) then 
        flowDoms(d2,level,sps)%iVar1(i2, j2, k2) = flowDoms(d1,level,sps)%iVar1(i1,j1,k1)
        flowDoms(d2,level,sps)%iVar2(i2, j2, k2) = flowDoms(d1,level,sps)%iVar2(i1,j1,k1)
        flowDoms(d2,level,sps)%iVar3(i2, j2, k2) = flowDoms(d1,level,sps)%iVar3(i1,j1,k1)

     else if (nVar == 4) then 
        flowDoms(d2,level,sps)%iVar1(i2, j2, k2) = flowDoms(d1,level,sps)%iVar1(i1,j1,k1)
        flowDoms(d2,level,sps)%iVar2(i2, j2, k2) = flowDoms(d1,level,sps)%iVar2(i1,j1,k1)
        flowDoms(d2,level,sps)%iVar3(i2, j2, k2) = flowDoms(d1,level,sps)%iVar3(i1,j1,k1)
        flowDoms(d2,level,sps)%iVar4(i2, j2, k2) = flowDoms(d1,level,sps)%iVar4(i1,j1,k1)

     else if (nVar == 5) then 
        flowDoms(d2,level,sps)%iVar1(i2, j2, k2) = flowDoms(d1,level,sps)%iVar1(i1,j1,k1)
        flowDoms(d2,level,sps)%iVar2(i2, j2, k2) = flowDoms(d1,level,sps)%iVar2(i1,j1,k1)
        flowDoms(d2,level,sps)%iVar3(i2, j2, k2) = flowDoms(d1,level,sps)%iVar3(i1,j1,k1)
        flowDoms(d2,level,sps)%iVar4(i2, j2, k2) = flowDoms(d1,level,sps)%iVar4(i1,j1,k1)
        flowDoms(d2,level,sps)%iVar5(i2, j2, k2) = flowDoms(d1,level,sps)%iVar5(i1,j1,k1)

     else if (nVar == 6) then 
        flowDoms(d2,level,sps)%iVar1(i2, j2, k2) = flowDoms(d1,level,sps)%iVar1(i1,j1,k1)
        flowDoms(d2,level,sps)%iVar2(i2, j2, k2) = flowDoms(d1,level,sps)%iVar2(i1,j1,k1)
        flowDoms(d2,level,sps)%iVar3(i2, j2, k2) = flowDoms(d1,level,sps)%iVar3(i1,j1,k1)
        flowDoms(d2,level,sps)%iVar4(i2, j2, k2) = flowDoms(d1,level,sps)%iVar4(i1,j1,k1)
        flowDoms(d2,level,sps)%iVar5(i2, j2, k2) = flowDoms(d1,level,sps)%iVar5(i1,j1,k1)
        flowDoms(d2,level,sps)%iVar6(i2, j2, k2) = flowDoms(d1,level,sps)%iVar6(i1,j1,k1)

     else if (nVar == 7) then 
        flowDoms(d2,level,sps)%iVar1(i2, j2, k2) = flowDoms(d1,level,sps)%iVar1(i1,j1,k1)
        flowDoms(d2,level,sps)%iVar2(i2, j2, k2) = flowDoms(d1,level,sps)%iVar2(i1,j1,k1)
        flowDoms(d2,level,sps)%iVar3(i2, j2, k2) = flowDoms(d1,level,sps)%iVar3(i1,j1,k1)
        flowDoms(d2,level,sps)%iVar4(i2, j2, k2) = flowDoms(d1,level,sps)%iVar4(i1,j1,k1)
        flowDoms(d2,level,sps)%iVar5(i2, j2, k2) = flowDoms(d1,level,sps)%iVar5(i1,j1,k1)
        flowDoms(d2,level,sps)%iVar6(i2, j2, k2) = flowDoms(d1,level,sps)%iVar6(i1,j1,k1)
        flowDoms(d2,level,sps)%iVar7(i2, j2, k2) = flowDoms(d1,level,sps)%iVar7(i1,j1,k1)

     else if (nVar == 8) then 
        flowDoms(d2,level,sps)%iVar1(i2, j2, k2) = flowDoms(d1,level,sps)%iVar1(i1,j1,k1)
        flowDoms(d2,level,sps)%iVar2(i2, j2, k2) = flowDoms(d1,level,sps)%iVar2(i1,j1,k1)
        flowDoms(d2,level,sps)%iVar3(i2, j2, k2) = flowDoms(d1,level,sps)%iVar3(i1,j1,k1)
        flowDoms(d2,level,sps)%iVar4(i2, j2, k2) = flowDoms(d1,level,sps)%iVar4(i1,j1,k1)
        flowDoms(d2,level,sps)%iVar5(i2, j2, k2) = flowDoms(d1,level,sps)%iVar5(i1,j1,k1)
        flowDoms(d2,level,sps)%iVar6(i2, j2, k2) = flowDoms(d1,level,sps)%iVar6(i1,j1,k1)
        flowDoms(d2,level,sps)%iVar7(i2, j2, k2) = flowDoms(d1,level,sps)%iVar7(i1,j1,k1)
        flowDoms(d2,level,sps)%iVar8(i2, j2, k2) = flowDoms(d1,level,sps)%iVar8(i1,j1,k1)

     else if (nVar == 9) then 
        flowDoms(d2,level,sps)%iVar1(i2, j2, k2) = flowDoms(d1,level,sps)%iVar1(i1,j1,k1)
        flowDoms(d2,level,sps)%iVar2(i2, j2, k2) = flowDoms(d1,level,sps)%iVar2(i1,j1,k1)
        flowDoms(d2,level,sps)%iVar3(i2, j2, k2) = flowDoms(d1,level,sps)%iVar3(i1,j1,k1)
        flowDoms(d2,level,sps)%iVar4(i2, j2, k2) = flowDoms(d1,level,sps)%iVar4(i1,j1,k1)
        flowDoms(d2,level,sps)%iVar5(i2, j2, k2) = flowDoms(d1,level,sps)%iVar5(i1,j1,k1)
        flowDoms(d2,level,sps)%iVar6(i2, j2, k2) = flowDoms(d1,level,sps)%iVar6(i1,j1,k1)
        flowDoms(d2,level,sps)%iVar7(i2, j2, k2) = flowDoms(d1,level,sps)%iVar7(i1,j1,k1)
        flowDoms(d2,level,sps)%iVar8(i2, j2, k2) = flowDoms(d1,level,sps)%iVar8(i1,j1,k1)
        flowDoms(d2,level,sps)%iVar9(i2, j2, k2) = flowDoms(d1,level,sps)%iVar9(i1,j1,k1)
     end if

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

        if (nVar == 1) then 
           jj = jj + 1
           flowDoms(d2,level,sps)%iVar1(i2, j2, k2) = recvBufInt(jj)

        else if (nVar == 2) then 
           jj = jj + 1
           flowDoms(d2,level,sps)%iVar1(i2, j2, k2) = recvBufInt(jj)
           jj = jj + 1
           flowDoms(d2,level,sps)%iVar2(i2, j2, k2) = recvBufInt(jj)

        else if (nVar == 3) then 
           jj = jj + 1
           flowDoms(d2,level,sps)%iVar1(i2, j2, k2) = recvBufInt(jj)
           jj = jj + 1
           flowDoms(d2,level,sps)%iVar2(i2, j2, k2) = recvBufInt(jj)
           jj = jj + 1
           flowDoms(d2,level,sps)%iVar3(i2, j2, k2) = recvBufInt(jj)

        else if (nVar == 4) then 
           jj = jj + 1
           flowDoms(d2,level,sps)%iVar1(i2, j2, k2) = recvBufInt(jj)
           jj = jj + 1
           flowDoms(d2,level,sps)%iVar2(i2, j2, k2) = recvBufInt(jj)
           jj = jj + 1
           flowDoms(d2,level,sps)%iVar3(i2, j2, k2) = recvBufInt(jj)
           jj = jj + 1
           flowDoms(d2,level,sps)%iVar4(i2, j2, k2) = recvBufInt(jj)

        else if (nVar == 5) then 
           jj = jj + 1
           flowDoms(d2,level,sps)%iVar1(i2, j2, k2) = recvBufInt(jj)
           jj = jj + 1
           flowDoms(d2,level,sps)%iVar2(i2, j2, k2) = recvBufInt(jj)
           jj = jj + 1
           flowDoms(d2,level,sps)%iVar3(i2, j2, k2) = recvBufInt(jj)
           jj = jj + 1
           flowDoms(d2,level,sps)%iVar4(i2, j2, k2) = recvBufInt(jj)
           jj = jj + 1
           flowDoms(d2,level,sps)%iVar5(i2, j2, k2) = recvBufInt(jj)

        else if (nVar == 6) then 
           jj = jj + 1
           flowDoms(d2,level,sps)%iVar1(i2, j2, k2) = recvBufInt(jj)
           jj = jj + 1
           flowDoms(d2,level,sps)%iVar2(i2, j2, k2) = recvBufInt(jj)
           jj = jj + 1
           flowDoms(d2,level,sps)%iVar3(i2, j2, k2) = recvBufInt(jj)
           jj = jj + 1
           flowDoms(d2,level,sps)%iVar4(i2, j2, k2) = recvBufInt(jj)
           jj = jj + 1
           flowDoms(d2,level,sps)%iVar5(i2, j2, k2) = recvBufInt(jj)
           jj = jj + 1
           flowDoms(d2,level,sps)%iVar6(i2, j2, k2) = recvBufInt(jj)

        else if (nVar == 7) then 
           jj = jj + 1
           flowDoms(d2,level,sps)%iVar1(i2, j2, k2) = recvBufInt(jj)
           jj = jj + 1
           flowDoms(d2,level,sps)%iVar2(i2, j2, k2) = recvBufInt(jj)
           jj = jj + 1
           flowDoms(d2,level,sps)%iVar3(i2, j2, k2) = recvBufInt(jj)
           jj = jj + 1
           flowDoms(d2,level,sps)%iVar4(i2, j2, k2) = recvBufInt(jj)
           jj = jj + 1
           flowDoms(d2,level,sps)%iVar5(i2, j2, k2) = recvBufInt(jj)
           jj = jj + 1
           flowDoms(d2,level,sps)%iVar6(i2, j2, k2) = recvBufInt(jj)
           jj = jj + 1
           flowDoms(d2,level,sps)%iVar7(i2, j2, k2) = recvBufInt(jj)

        else if (nVar == 8) then 
           jj = jj + 1
           flowDoms(d2,level,sps)%iVar1(i2, j2, k2) = recvBufInt(jj)
           jj = jj + 1
           flowDoms(d2,level,sps)%iVar2(i2, j2, k2) = recvBufInt(jj)
           jj = jj + 1
           flowDoms(d2,level,sps)%iVar3(i2, j2, k2) = recvBufInt(jj)
           jj = jj + 1
           flowDoms(d2,level,sps)%iVar4(i2, j2, k2) = recvBufInt(jj)
           jj = jj + 1
           flowDoms(d2,level,sps)%iVar5(i2, j2, k2) = recvBufInt(jj)
           jj = jj + 1
           flowDoms(d2,level,sps)%iVar6(i2, j2, k2) = recvBufInt(jj)
           jj = jj + 1
           flowDoms(d2,level,sps)%iVar7(i2, j2, k2) = recvBufInt(jj)
           jj = jj + 1
           flowDoms(d2,level,sps)%iVar8(i2, j2, k2) = recvBufInt(jj)

        else if (nVar == 9) then 
           jj = jj + 1
           flowDoms(d2,level,sps)%iVar1(i2, j2, k2) = recvBufInt(jj)
           jj = jj + 1
           flowDoms(d2,level,sps)%iVar2(i2, j2, k2) = recvBufInt(jj)
           jj = jj + 1
           flowDoms(d2,level,sps)%iVar3(i2, j2, k2) = recvBufInt(jj)
           jj = jj + 1
           flowDoms(d2,level,sps)%iVar4(i2, j2, k2) = recvBufInt(jj)
           jj = jj + 1
           flowDoms(d2,level,sps)%iVar5(i2, j2, k2) = recvBufInt(jj)
           jj = jj + 1
           flowDoms(d2,level,sps)%iVar6(i2, j2, k2) = recvBufInt(jj)
           jj = jj + 1
           flowDoms(d2,level,sps)%iVar7(i2, j2, k2) = recvBufInt(jj)
           jj = jj + 1
           flowDoms(d2,level,sps)%iVar8(i2, j2, k2) = recvBufInt(jj)
           jj = jj + 1
           flowDoms(d2,level,sps)%iVar9(i2, j2, k2) = recvBufInt(jj)
        end if
     enddo

  enddo completeRecvs

  ! Complete the nonblocking sends.

  size = commPattern(level)%nProcSend
  do i=1,commPattern(level)%nProcSend
     call mpi_waitany(size, sendRequests, index, status, ierr)
  enddo


end subroutine whalo1to1IntGeneric

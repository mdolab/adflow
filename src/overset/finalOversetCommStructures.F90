subroutine finalOversetCommStructures(level, sps, MAGIC)


  ! We need to fill in the following information in the comm patterns:
  ! sendProc, nProcSend, nSend, and sendList for each proc
  ! recvProc, nProcRecv, nRecv, and recvList for each proc

  use blockPointers
  use communication
  use overset

  implicit none

  ! Input Parameters
  integer(kind=intType), intent(in) :: level, sps, MAGIC

  ! Working Parameters
  integer(kind=intType) :: i, j, k, ii, jj, kk, nn, tag, ierr, nLocalFringe
  integer(kind=intType) :: n, nFringeProc, sendCount, iProc
  integer(kind=intType) :: nCopy,  iSize, iStart, iEnd, iProcRecv, iSendProc, iRecvProc
  type(fringeType), dimension(:), allocatable :: localFringes
  integer(kind=intType), dimension(:), allocatable :: nProcSend, nProcSendLocal
  integer(kind=intType), dimension(:), allocatable :: fringeProc, cumFringeProc

  ! MPI Stuff
  integer status(MPI_STATUS_SIZE) 
  integer, allocatable, dimension(:) :: sendRequests1
  integer(kind=intTYpe) :: nMsgToSend, iRecv, nRecv
  logical :: flag

  ! We need to fill in the following information in the comm patterns:
  ! sendProc, nProcSend, nSend, and sendList for each proc
  ! recvProc, nProcRecv, nRecv, and recvList for each proc

  ! Count ALL fringes including the halos
  nLocalFringe = 0
  do nn=1, nDom
     call setPointers(nn, level, sps)
     do k=0, kb
        do j=0, jb
           do i=0, ib
              if (fringes(i, j, k)%donorProc /= -1) then 
                 nLocalFringe = nLocalFringe + 1
              end if
           end do
        end do
     end do
  end do

  allocate(localFringes(nLocalFringe))
  ! Now add in all of our local fringes
  nLocalFringe = 0
  do nn=1, nDom
     call setPointers(nn, level, sps)
     do k=0, kb
        do j=0, jb
           do i=0, ib
              if (fringes(i, j, k)%donorProc /= -1) then 
                 nLocalFringe = nLocalFringe + 1
                 localFringes(nLocalFringe) = fringes(i, j, k)
              end if
           end do
        end do
     end do
  end do

  ! Sort the actual fringes
  call qsortFringeType(localFringes, nLocalFringe)

  ! On the first pass with our local data we can determine the
  ! internal nCopy as well as the nProcRecv. We also increment
  ! nProcSend, and then run an allreduce to determine the number of
  ! processors I need to send stuff to. This will then allow us to
  ! allocate the sendData as well.
  allocate(fringeProc(nProc), cumFringeProc(1:nProc+1))

  call computeFringeProcArray(localFringes, nLocalFringe, &
       fringeProc, cumFringeProc, nFringeProc)

  nCopy = 0
  nProcRecv = 0
  allocate(nProcSend(nProc), nProcSendLocal(nProc))
  nProcSendLocal = 0

  do i=1, nFringeProc ! The numer of processors i'm dealing with
     if (fringeProc(i) == myid) then 
        nCopy = cumFringeProc(i+1) - cumFringeProc(i)
     else
        nProcRecv = nProcRecv + 1
        nProcSendLocal(fringeProc(i)+1) = 1
     end if
  end do

  ! This will determine the number of procs each processor has to
  ! send donors to. The value we care about is nSend(myid). Also,
  ! this counts as a barrier so we know that all comms have
  ! finished. That means we can restart our MAGIC tag back to 1.
  call mpi_allreduce(nProcSendLocal, nProcSend, nProc, sumb_integer, MPI_SUM, &
       sumb_comm_world, ierr)
  call ECHK(ierr, __FILE__, __LINE__)

  ! ---------------------------------------------------------------------
  ! Since we are setting up...ie doing the transpose of the actual
  ! final comm structure, the number of sends we make in this routine,
  ! is the same is the number of receives we actually want to do in
  ! the final comm structure. We need to allocate the sendRequests1
  ! and to know the number of recevies and since this information is
  ! already computed we assign it here. I know it looks confusing. 
  allocate(sendRequests1(nProcRecv))
  nRecv = nProcSend(myid+1)
  ! ---------------------------------------------------------------------

  ! We can allocate all necessary space for the send and receive information
  commPatternOverset(level, sps)%nProcRecv = nProcRecv
  allocate(commPatternOverset(level, sps)%recvProc(nProcRecv))
  allocate(commPatternOverset(level, sps)%nRecv(nProcRecv))
  allocate(commPatternOverset(level, sps)%recvList(nProcRecv))

  commPatternOverset(level, sps)%nProcSend = nProcSend(myid+1)
  allocate(commPatternOverset(level, sps)%sendProc(nProcSend(myid+1)))
  allocate(commPatternOverset(level, sps)%nSend(nProcSend(myid+1)))
  allocate(commPatternOverset(level, sps)%sendList(nProcSend(myid+1)))

  ! As well as the copy information
  internalOverset(level, sps)%nCopy = nCopy
  allocate(internalOverset(level, sps)%donorBlock(nCopy))
  allocate(internalOverset(level, sps)%donorIndices(nCopy, 3))
  allocate(internalOverset(level, sps)%donorInterp(nCopy, 8))
  allocate(internalOverset(level, sps)%haloBlock(nCopy))
  allocate(internalOverset(level, sps)%haloIndices(nCopy, 3))

  sendCount = 0
  iRecvProc = 0
  do i=1, nFringeProc ! The numer of processors i'm dealing with
     iSize = cumFringeProc(i+1) - cumFringeProc(i)
     iStart = cumFringeProc(i)
     iEnd   = cumFringeProc(i+1)-1

     if (fringeProc(i) == myid) then 
        ii =0
        do j=iStart,iEnd
           ii = ii + 1
           ! This is the donor information
           internalOverset(level, sps)%donorBlock(ii) = localFringes(j)%donorBlock
           internalOverset(level, sps)%donorIndices(ii, 1) = localFringes(j)%dI
           internalOverset(level, sps)%donorIndices(ii, 2) = localFringes(j)%dJ
           internalOverset(level, sps)%donorIndices(ii, 3) = localFringes(j)%dK
           call fracToWeights(localFringes(j)%donorFrac, internalOverset(level, sps)%donorInterp(ii, :))

           ! And the receiver (halo) information
           internalOverset(level, sps)%haloBlock(ii) = localFringes(j)%myBlock
           internalOverset(level, sps)%haloIndices(ii, 1) = localFringes(j)%myI
           internalOverset(level, sps)%haloIndices(ii, 2) = localFringes(j)%myJ
           internalOverset(level, sps)%haloIndices(ii, 3) = localFringes(j)%myK

        end do

     else
        ! Set the receiver info and send the donor info to the donor proc
        iRecvProc = iRecvProc + 1
        commPatternOverset(level, sps)%recvProc(iRecvProc) = fringeProc(i)
        commPatternOverset(level, sps)%nRecv(iRecvProc) = iSize
        allocate(commPatternOverset(level, sps)%recvList(iRecvProc)%block(iSize))
        allocate(commPatternOverset(level, sps)%recvList(iRecvProc)%indices(iSize, 3))

        ! Fill up the local information from myBlock and my{I,j,k}
        ii = 0
        do j=iStart, iEnd
           ii = ii + 1
           commPatternOverset(level, sps)%recvList(iRecvProc)%block(ii) = localFringes(j)%myBlock
           commPatternOverset(level, sps)%recvList(iRecvProc)%indices(ii,1) = localFringes(j)%myI
           commPatternOverset(level, sps)%recvList(iRecvProc)%indices(ii,2) = localFringes(j)%myJ
           commPatternOverset(level, sps)%recvList(iRecvProc)%indices(ii,3) = localFringes(j)%myK
        end do

        ! Now iSend these fringes to where they need to go: Start
        ! the tags at MAGIC not 1 since we will call
        ! exchangeIBlanks below which have tags from 0:nProc
        tag = MAGIC + myID + 1
        sendCount = sendCount + 1
        call mpi_isend(localFringes(iStart:iend), iSize, oversetMPIFringe, &
             fringeProc(i), tag, SUmb_comm_world, sendRequests1(sendCount), ierr)
     end if
  end do

  iSendProc = 0
  iRecv = 0
  do while (iRecv < nRecv)

     call MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, SUmb_comm_world, flag, status, ierr)
     call ECHK(ierr, __FILE__, __LINE__)

     ! Check if a message is ready
     tag = status(MPI_TAG)
     if (flag .and. tag >= MAGIC+1 .and. tag <= 2*MAGIC) then
        iRecv = iRecv + 1

        call MPI_Get_count(status, oversetMPIFringe, n, ierr)
        call ECHK(ierr, __FILE__, __LINE__)

        ! Allocate space for temporary fringes if not big enough
        if (n > size(tmpFringes)) then 
           deallocate(tmpFringes)
           allocate(tmpFringes(n))
        end if

        ! Now actually receive the fringes
        call mpi_recv(tmpFringes, n, oversetMPIFringe, status(MPI_SOURCE), &
             status(MPI_TAG), SUmb_comm_world, status, ierr)
        call ECHK(ierr, __FILE__, __LINE__)

        ! We don't actually know the order that the fringes will be
        ! received. this is ok. We just put them in the order we get
        ! them. Increment the iSendProc. 
        iSendProc = iSendProc + 1

        commPatternOverset(level, sps)%sendProc(iSendProc) = status(MPI_TAG) - MAGIC -1 
        commPatternOverset(level, sps)%nSend(iSendProc) = n
        allocate(& 
             commPatternOverset(level, sps)%sendList(iSendProc)%block(n), &
             commPatternOverset(level, sps)%sendList(iSendProc)%indices(n, 3), &
             commPatternOverset(level, sps)%sendList(iSendProc)%interp(n, 8))

        ! Now set the data
        do i=1, n
           commPatternOverset(level, sps)%sendList(iSendProc)%block(i) =  tmpFringes(i)%donorBlock
           commPatternOverset(level, sps)%sendList(iSendProc)%indices(i, :) = &
                (/tmpFringes(i)%dI, tmpFringes(i)%dJ, tmpFringes(i)%dK/)
           call fracToWeights(tmpFringes(i)%donorFrac, &
                commPatternOverset(level, sps)%sendList(iSendProc)%interp(i, :))
        end do
     end if
  end do

  ! One last thing to do is to create the cumulative forms of nSend
  ! and nRecv (nSendCum and nRecvCum)
  allocate(commPatternOverset(level, sps)%nSendCum(0:nProcSend(myid+1)), &
       commPatternOverset(level, sps)%nRecvCum(0:nProcRecv))

  call getCumulativeForm(commPatternOverset(level, sps)%nSend, commPatternOverset(level, sps)%nProcSend, &
       commPatternOverset(level, sps)%nSendCum)

  call getCumulativeForm(commPatternOverset(level, sps)%nRecv, commPatternOverset(level, sps)%nProcRecv, &
       commPatternOverset(level, sps)%nRecvCum)

  call mpi_waitall(sendCount, sendRequests1, MPI_STATUSES_IGNORE, ierr)
  call ECHK(ierr, __FILE__, __LINE__)

  deallocate(localFringes, fringeProc, cumFringeProc, nProcSend)

  ! Put a barrier here to make sure that the iprobe only get the
  ! messages from this routine
  call mpi_barrier(sumb_comm_world, ierr)
  call ECHK(ierr, __FILE__, __LINE__)

end subroutine finalOversetCommStructures

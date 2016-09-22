subroutine finalOversetCommStructures(level, sps)


  ! We need to fill in the following information in the comm patterns:
  ! sendProc, nProcSend, nSend, and sendList for each proc
  ! recvProc, nProcRecv, nRecv, and recvList for each proc
  use constants
  use blockPointers, only : nDom, fringes, ib, jb, kb
  use overset, only: fringeType
  use communication, only : adflow_comm_world, myid, nProc, sendRequests, recvRequests, &
       commPatternOverset, internalOverset
  use utils, only : setPointers, terminate, EChk

  implicit none

  ! Input Parameters
  integer(kind=intType), intent(in) :: level, sps

  ! Working Parameters
  integer(kind=intType) :: i, j, k, ii, jj, kk, nn, tag, ierr, nLocalFringe
  integer(kind=intType) :: n, nFringeProc, sendCount, recvCount, iProc
  integer(kind=intType) :: iSize, iStart, iEnd, iProcRecv, iSendProc, iRecvProc
  integer(kind=intType) :: nProcSend, nProcRecv, nCopy, totalRecvSize
  type(fringeType), dimension(:), allocatable :: localFringes
  integer(kind=intType), dimension(:), allocatable :: tmpInt
  integer(kind=intType), dimension(:), allocatable :: recvSizes
  integer(kind=intType), dimension(:), allocatable :: nProcSendLocal
  integer(kind=intType), dimension(:), allocatable :: nProcSendLocaltmp
  integer(kind=intType), dimension(:), allocatable :: fringeProc, cumFringeProc
  integer(kind=intType), dimension(:), allocatable :: intSendBuf, intRecvBuf
  real(kind=realType), dimension(:), allocatable :: realSendBuf, realRecvBuf
  integer status(MPI_STATUS_SIZE) 

  ! We need to fill in the following information in the comm patterns:
  ! sendProc, nProcSend, nSend, and sendList for each proc
  ! recvProc, nProcRecv, nRecv, and recvList for each proc

  ! Count all actual fringes (including second level halos)
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

  allocate(tmpInt(0:nProc-1), recvSizes(0:nProc-1))
  tmpInt = 0
  do j=1, nFringeProc
     iProc = fringeProc(j)
     if (iProc /= myid) then 
        tmpInt(iProc) = (cumFringeProc(j+1) - cumFringeProc(j))
     end if
  end do

  ! Sum how much data we must receive from each processor. 
  call mpi_alltoall(tmpInt, 1, adflow_integer, recvSizes, 1, adflow_integer, &
       adflow_comm_world, ierr)
  call ECHK(ierr, __FILE__, __LINE__)
  deallocate(tmpInt)

  ! Allocate space for the sending and receiving buffers
  totalRecvSize = sum(recvSizes)
  allocate(intSendBuf(4*nLocalFringe), intRecvBuf(totalRecvSize*4), &
       realSendBuf(3*nLocalFringe), realRecvBuf(totalRecvSize*3))
  
  ! Pack the real and integer buffers with donorBlock, dI, dJ, dK and
  ! donorFrac. We are putting everything in here, including our
  ! own. That's ok.
  do j=1, nLocalFringe
     intSendBuf(4*j-3) = localFringes(j)%donorBlock
     intSendbuf(4*j-2) = localFringes(j)%dI
     intSendBuf(4*j-1) = localFringes(j)%dJ
     intSendBuf(4*j  ) = localFringes(j)%dK
     realSendBuf(3*j-2:3*j) = localFringes(j)%donorFrac
  end do

  nCopy = 0
  nProcRecv = 0
  allocate(nProcSendLocal(0:nProc-1))
  nProcSendLocal = 0

  do i=1, nFringeProc ! The numer of processors i'm dealing with
     if (fringeProc(i) == myid) then 
        nCopy = cumFringeProc(i+1) - cumFringeProc(i)
     else
        nProcRecv = nProcRecv + 1  ! I will receive something from this proc
        nProcSendLocal(fringeProc(i)) = 1 ! That proc will send it to me
     end if
  end do

  ! This will sum up the nProcSendLocal array and then send out the
  ! number of sends I have to do. 
  ! call mpi_reduce_scatter_block(nProcSendLocal, nProcSend, 1, &
  !      adflow_integer, MPI_SUM, adflow_comm_world, ierr)

 !
  ! The following is done for MPI 2.0 compatibility. 
  ! ------------------------------------------------
  ! Break mpi_reduce_scatter_block to two steps: mpi_reduce and mpi_scatter
  allocate(nProcSendLocaltmp(0:nProc-1))
  nProcSendLocaltmp = 0

  ! Step 1: Reduce at root proc
  call mpi_reduce(nProcSendLocal, nProcSendLocaltmp, nProc, &
           adflow_integer, MPI_SUM, 0, adflow_comm_world, ierr)

  ! Step 2: Scatter from root proc
  ! sendbuf = nProcSendLocaltmp, sendcount = 1, 
  ! recvbuf = nProcSend, recvcount = 1, source = 0
  call mpi_scatter(nProcSendLocaltmp, 1, adflow_integer, nProcSend, 1, &
           adflow_integer, 0, adflow_comm_world, ierr)

  deallocate(nProcSendLocaltmp, nProcSendLocal)

  ! We can allocate all necessary space for the send and receive information
  commPatternOverset(level, sps)%nProcRecv = nProcRecv
  allocate(commPatternOverset(level, sps)%recvProc(nProcRecv))
  allocate(commPatternOverset(level, sps)%nRecv(nProcRecv))
  allocate(commPatternOverset(level, sps)%recvList(nProcRecv))

  commPatternOverset(level, sps)%nProcSend = nProcSend
  allocate(commPatternOverset(level, sps)%sendProc(nProcSend))
  allocate(commPatternOverset(level, sps)%nSend(nProcSend))
  allocate(commPatternOverset(level, sps)%sendList(nProcSend))

  ! As well as the copy information
  internalOverset(level, sps)%nCopy = nCopy
  allocate(internalOverset(level, sps)%donorBlock(nCopy))
  allocate(internalOverset(level, sps)%donorIndices(nCopy, 3))
  allocate(internalOverset(level, sps)%donorInterp(nCopy, 8))
  allocate(internalOverset(level, sps)%haloBlock(nCopy))
  allocate(internalOverset(level, sps)%haloIndices(nCopy, 3))

  ! Send the donors back to their own processors.
  sendCount = 0
  do j=1, nFringeProc
     
     iProc = fringeProc(j)
     iStart = cumFringeProc(j)-1
     iSize = cumFringeProc(j+1) - cumFringeProc(j)
     
     if (iProc /= myid) then 
        sendCount = sendCount + 1
        call mpi_isend(intSendBuf(iStart*4+1), 4*iSize, adflow_integer, iProc, myid, &
             adflow_comm_world, sendRequests(sendCount), ierr)
        call ECHK(ierr, __FILE__, __LINE__)

        sendCount = sendCount + 1
        call mpi_isend(realSendBuf(iStart*3+1), 3*iSize, adflow_real, iProc, myid, &
             adflow_comm_world, sendRequests(sendCount), ierr)
        call ECHK(ierr, __FILE__, __LINE__)
     end if
  end do
  
  ! Non-blocking receives
  recvCount = 0
  ii = 1
  jj = 1
  do iProc=0, nProc-1
     
     if (recvSizes(iProc) > 0) then
        recvCount = recvCount + 1
        call mpi_irecv(intRecvBuf(ii), 4*recvSizes(iProc), adflow_integer, &
             iProc, iProc, adflow_comm_world, recvRequests(recvCount), ierr) 
        call ECHK(ierr, __FILE__, __LINE__) 
       
        ii = ii + recvSizes(iProc)*4

        recvCount = recvCount + 1
        call mpi_irecv(realRecvBuf(jj), 3*recvSizes(iProc), adflow_real, &
             iProc, iProc, adflow_comm_world, recvRequests(recvCount), ierr) 
        call ECHK(ierr, __FILE__, __LINE__) 
        jj = jj + recvSizes(iProc)*3
     end if
  end do

  ! Do a little local work while we wait for the data to send/recv
  iRecvProc = 0
  do i=1, nFringeProc
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

        ! Set the receiver info. The info is already sent and in flight
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
     end if
  end do

  ! Complete all the sends/receives. We could do overlapping here
  ! like the frist comm for the fringes/blocks. 
  call mpi_waitall(recvCount, recvRequests, MPI_STATUSES_IGNORE, ierr)
  call ECHK(ierr, __FILE__, __LINE__)
  
  call mpi_waitall(sendCount, sendRequests, MPI_STATUSES_IGNORE, ierr)
  call ECHK(ierr, __FILE__, __LINE__)

  ! All of our data has now arrived we can now finish completing the send information.
  ii = 0
  jj = 0
  iSendProc = 0 ! running counter of the ith processor
  do iProc=0, nProc-1
     if (recvSizes(iProc)> 0) then
        ! We should have received something from this processor
        iSendProc = iSendProc + 1
        commPatternOverset(level, sps)%sendProc(iSendProc) = iProc
        n = recvSizes(iProc)
        commPatternOverset(level, sps)%nSend(iSendProc) = n
        allocate(& 
             commPatternOverset(level, sps)%sendList(iSendProc)%block(n), &
             commPatternOverset(level, sps)%sendList(iSendProc)%indices(n, 3), &
             commPatternOverset(level, sps)%sendList(iSendProc)%interp(n, 8))

        ! Now set the data
        do i=1, n
           commPatternOverset(level, sps)%sendList(iSendProc)%block(i) = intRecvBuf(ii+1)
           commPatternOverset(level, sps)%sendList(iSendProc)%indices(i, :) = &
                intRecvBuf(ii+2:ii+4)
           ii = ii + 4

           call fracToWeights(realRecvBuf(jj+1:jj+3), &
                commPatternOverset(level, sps)%sendList(iSendProc)%interp(i, :))
           jj = jj + 3
        end do
     end if
  end do

  ! One last thing to do is to create the cumulative forms of nSend
  ! and nRecv (nSendCum and nRecvCum)
  allocate(commPatternOverset(level, sps)%nSendCum(0:nProcSend), &
       commPatternOverset(level, sps)%nRecvCum(0:nProcRecv))

  call getCumulativeForm(commPatternOverset(level, sps)%nSend, commPatternOverset(level, sps)%nProcSend, &
       commPatternOverset(level, sps)%nSendCum)

  call getCumulativeForm(commPatternOverset(level, sps)%nRecv, commPatternOverset(level, sps)%nProcRecv, &
       commPatternOverset(level, sps)%nRecvCum)

  deallocate(localFringes, fringeProc, cumFringeProc, &
       intRecvBuf, intSendBuf, realRecvBuf, realSendBuf)

end subroutine finalOversetCommStructures

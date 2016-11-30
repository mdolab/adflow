subroutine finalOversetCommStructures(level, sps)


  ! We need to fill in the following information in the comm patterns:
  ! sendProc, nProcSend, nSend, and sendList for each proc
  ! recvProc, nProcRecv, nRecv, and recvList for each proc
  use constants
  use block, only : flowDoms
  use blockPointers, only : nDom, fringes, ib, jb, kb, status, fringePtr
  use overset, only: fringeType
  use communication, only : adflow_comm_world, myid, nProc, sendRequests, recvRequests, &
       commPatternOverset, internalOverset
  use utils, only : setPointers, terminate, EChk
  use oversetUtilities, only : fracToWeights, qsortFringeType, getCumulativeForm, computeFringeProcArray, &
       isReceiver, unwindIndex
  implicit none

  ! Input Parameters
  integer(kind=intType), intent(in) :: level, sps

  ! Working Parameters
  integer(kind=intType) :: i, j, k, ii, jj, kk, nn, tag, ierr, nLocalFringe
  integer(kind=intType) :: myI, myJ, myK, dI, dJ, dK, il, jl, kl, myBlock, dBlock
  integer(kind=intType) :: n, nFringeProc, sendCount, recvCount, iProc
  integer(kind=intType) :: iSize, iStart, iEnd, iProcRecv, iSendProc, iRecvProc
  integer(kind=intType) :: nProcSend, nProcRecv, nCopy, totalRecvSize, index
  type(fringeType), dimension(:), allocatable :: localFringes
  integer(kind=intType), dimension(:), allocatable :: tmpInt
  integer(kind=intType), dimension(:), allocatable :: recvSizes
  integer(kind=intType), dimension(:), allocatable :: nProcSendLocal
  integer(kind=intType), dimension(:), allocatable :: nProcSendLocaltmp
  integer(kind=intType), dimension(:), allocatable :: fringeProc, cumFringeProc
  integer(kind=intType), dimension(:), allocatable :: intSendBuf, intRecvBuf
  real(kind=realType), dimension(:), allocatable :: realSendBuf, realRecvBuf
  integer mpiStatus(MPI_STATUS_SIZE) 

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
              if (isReceiver(status(i, j, k))) then 
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
              if (isReceiver(status(i, j, k))) then 
                 nLocalFringe = nLocalFringe + 1
                 ii = fringePtr(1, i, j, k)
                 localFringes(nLocalFringe) = fringes(ii)
              end if
           end do
        end do
     end do
  end do

  ! Sort the actual fringes
  call qsortFringeType(localFringes, nLocalFringe, sortByDonor)

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
  allocate(intSendBuf(2*nLocalFringe), intRecvBuf(totalRecvSize*2), &
       realSendBuf(3*nLocalFringe), realRecvBuf(totalRecvSize*3))
  
  ! Pack the real and integer buffers with donorBlock, dIndex and
  ! donorFrac. We are putting everything in here, including our
  ! own. That's ok.
  do j=1, nLocalFringe
     intSendBuf(2*j-1) = localFringes(j)%donorBlock
     intSendbuf(2*j  ) = localFringes(j)%dindex
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
  allocate(internalOverset(level, sps)%donorInterpd(nCopy, 8))
  allocate(internalOverset(level, sps)%XCen(nCopy, 3))
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
        call mpi_isend(intSendBuf(iStart*2+1), 2*iSize, adflow_integer, iProc, myid, &
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
        call mpi_irecv(intRecvBuf(ii), 2*recvSizes(iProc), adflow_integer, &
             iProc, iProc, adflow_comm_world, recvRequests(recvCount), ierr) 
        call ECHK(ierr, __FILE__, __LINE__) 
       
        ii = ii + recvSizes(iProc)*2

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
           dBlock = localFringes(j)%donorBlock
           internalOverset(level, sps)%donorBlock(ii) = dBlock
           il = flowDoms(dBlock, level, sps)%il
           jl = flowDoms(dBlock, level, sps)%jl
           kl = flowDoms(dBlock, level, sps)%kl
           call unwindIndex(localFringes(j)%dIndex, il, jl, kl, dI, dJ, dK)

           internalOverset(level, sps)%donorIndices(ii, :) = (/dI, dJ, dk/)
           call fracToWeights(localFringes(j)%donorFrac, internalOverset(level, sps)%donorInterp(ii, :))
           internalOverset(level, sps)%donorInterpd(ii, :) = zero
           internalOverset(level, sps)%xCen(ii, :) = zero

           ! And the receiver (halo) information
           myBlock = localFringes(j)%myBlock
           internalOverset(level, sps)%haloBlock(ii) = myBlock

           il = flowDoms(myBlock, level, sps)%il
           jl = flowDoms(myBlock, level, sps)%jl
           kl = flowDoms(myBlock, level, sps)%kl
           call unwindIndex(localFringes(j)%myIndex, il, jl, kl, myI, myJ, myK)
           internalOverset(level, sps)%haloIndices(ii, :) = (/myI, myJ, myK/)
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
           myBlock = localFringes(j)%myBlock
           commPatternOverset(level, sps)%recvList(iRecvProc)%block(ii) = myBlock

           il = flowDoms(myBlock, level, sps)%il
           jl = flowDoms(myBlock, level, sps)%jl
           kl = flowDoms(myBlock, level, sps)%kl
           call unwindIndex(localFringes(j)%myIndex, il, jl, kl, myI, myJ, myK)
           commPatternOverset(level, sps)%recvList(iRecvProc)%indices(ii,:) = (/myI, myJ, myK/)
        end do
     end if
  end do

  ! Complete all the sends/receives. We could do overlapping here
  ! like the frist comm for the fringes/blocks. 

  do i=1, recvCount
     call mpi_waitany(recvCount, recvRequests, index, mpiStatus, ierr)
     call EChk(ierr,__FILE__,__LINE__)
  enddo

  do i=1, sendCount
     call mpi_waitany(sendCount, sendRequests, index, mpiStatus, ierr)
     call EChk(ierr,__FILE__,__LINE__)
  enddo

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
             commPatternOverset(level, sps)%sendList(iSendProc)%interp(n, 8), &
             commPatternOverset(level, sps)%sendList(iSendProc)%interpd(n, 8))

        ! Now set the data
        do i=1, n
           myBlock = intRecvBuf(ii+1)
           commPatternOverset(level, sps)%sendList(iSendProc)%block(i) = myBlock

           il = flowDoms(myBlock, level, sps)%il
           jl = flowDoms(myBlock, level, sps)%jl
           kl = flowDoms(myBlock, level, sps)%kl

           index = intRecvBuf(ii+2)
           call unWindIndex(index, il, jl, kl, dI, dJ, dK)
           commPatternOverset(level, sps)%sendList(iSendProc)%indices(i, :) = &
                (/dI, dJ, dK/)
           ii = ii + 2

           call fracToWeights(realRecvBuf(jj+1:jj+3), &
                commPatternOverset(level, sps)%sendList(iSendProc)%interp(i, :))
           commPatternOverset(level, sps)%sendList(iSendProc)%interpd(i, :) = zero
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

!      ******************************************************************
!      *                                                                *
!      * computeOversetInterpolation is the top level routine that      *
!      * implements the implicit hole cutting method for determing      *
!      * overset grid connectivitiies. It operates on a given multigrid *
!      * level and spectral instance                                    *
!      *                                                                *
!      ******************************************************************

subroutine oversetComm(level, firstTime, coarseLevel)
  use constants
  use communication
  use blockPointers
  use overset
  use stencils
  use BCTypes
  use inputTimeSpectral
  use ADTapi
  use inputOverset
  implicit none

  ! Input Parameters
  integer(kind=intType), intent(in) :: level
  logical :: firstTime, coarseLevel
  integer(kind=intType) :: sps

  ! Local Variables
  integer(kind=intType) :: i, ii, j, jj, k, kk, i_stencil
  integer(kind=intType) :: m, iSize, iStart, iEnd, index, rSize, nClusters
  integer(kind=intType) :: iDom, jDom, iDim, iPtr, rPtr
  integer(kind=intType) :: nn, mm, n, ierr, iProc, myIndex
  integer(kind=intType) :: iWork, nWork, nFringeProc, nLocalFringe
  real(kind=realType) :: startTime, endTime, quality
  logical :: computeCellFound, oversetPresent, isCompute

  type(CSRMatrix), pointer :: overlap
  type(CSRMatrix) :: overlapTranspose

  integer(kind=intType), dimension(:), allocatable :: clusters
  integer(kind=intType), dimension(:), allocatable :: cumFringeRecv, fringeRecvSizes
  integer(kind=intType), dimension(:, :), allocatable :: work, tmpInt2D

  real(kind=realType), dimension(:), allocatable :: tmpReal
  real(kind=realType), dimension(:, :), allocatable :: xMin, xMax

  logical, dimension(:), allocatable :: oBlockReady, oFringeReady, oWallReady

  type(oversetBlock), dimension(:), allocatable :: oBlocks
  type(oversetFringe), dimension(:), allocatable :: oFringes
  type(oversetWall), dimension(:), allocatable :: oWalls
  type(fringeType), dimension(:), allocatable :: localFringes

  ! MPI/Communication related
  integer status(MPI_STATUS_SIZE) 
  integer(kind=intType) :: MAGIC, source, tag, sendCount, recvCount
  integer(kind=intType) :: nOFringeSend, nOFringeRecv
  integer(kind=intType) :: nOBlockSend, nOBlockRecv
  integer(kind=intType) :: nOWallSend, nOWallRecv
  logical :: flag

  integer(kind=intType), dimension(:, :), allocatable :: oBlockSendList, oBlockRecvList
  integer(kind=intType), dimension(:, :), allocatable :: oFringeSendList, oFringeRecvList
  integer(kind=intType), dimension(:, :), allocatable :: oWallSendList, oWallRecvList
  integer(kind=intType), dimension(:, :), allocatable :: bufSizes, recvInfo
  integer(kind=intType), dimension(:), allocatable :: intRecvBuf
  real(kind=realType), dimension(:), allocatable :: realRecvBuf


  ! If there is not overset meshes present, just make an empty comm
  ! structure and call it a day. 
  if (.not. oversetPresent()) then 
     do sps=1,nTimeIntervalsSpectral
        call emptyOversetComm(level, sps)

        do nn=1, nDom
           flowDoms(nn,level,sps)%nOrphans = 0
           if (.not. associated(flowDoms(nn,level,sps)%iblank)) then
              i = flowDoms(nn,level,sps)%ib
              j = flowDoms(nn,level,sps)%jb
              k = flowDoms(nn,level,sps)%kb

              allocate(flowDoms(nn,level,sps)%iblank(0:i,0:j,0:k))
              flowDoms(nn, level, sps)%iblank = 1
           end if
        end do
     end do
     return
  end if

  ! -----------------------------------------------------------------
  ! Step 1: Initializaion: Make sure the stencils are initialized. 
  ! -----------------------------------------------------------------

  call initialize_stencils()

  ! -----------------------------------------------------------------
  ! Step 2: Communicate the block size info to everyone. Also generate
  ! cumDomProc which is the cumulative form. This will make our lives
  ! easier a little later on with indexing. (Routine below)
  ! -----------------------------------------------------------------

  allocate(nDomProc(0:nProc-1), cumDomProc(0:nProc))
  call blockSizeInfo

  ! Determine the magic number which is actually the same as
  ! nDomTotal. nDomTotal is guaranteed to be greater than or equal to
  ! nProc. This will be used to space out tags on communications to
  ! make sure they do not overlap. 
  MAGIC = nDomTotal

  ! Master SPS loop
  spectralLoop: do sps=1,nTimeIntervalsSpectral

     ! Set a pointer to make the code easier to read
     overlap => overlapMatrix(level, sps)

     ! -----------------------------------------------------------------
     ! Step 3: Next we determine the "clusters" that are present in the
     ! grid. We will use the original CGNS information since this is
     ! slightly more readily available. Futhermore, there will generally
     ! be fewer blocks in the CGNS file that compute blocks due to
     ! splitting. Once we have labelled the CGNS blocks, we can simply
     ! label the compute blocks by looking at what CGNS block they came
     ! from. Since all procs have this information we can simply run on
     ! all procs. 
     ! -----------------------------------------------------------------
     allocate(clusters(nDomTotal))
     call determineClusters(clusters, nDomTotal, cumDomProc, nClusters)

     ! -----------------------------------------------------------------
     ! Step 4: Compute the 3D axis oriented bounding boxes for each block
     ! and communicate then to everyone. Also determine communicate the
     ! minimum volume for each block to everyone.  (Routine below)
     ! -----------------------------------------------------------------

     allocate(xMin(3, nDomTotal), xMax(3, nDomTotal))
     call computeDomainBoundingBoxes
     ! -----------------------------------------------------------------
     ! Step 8: Build a global sparse matrix representation of the overlap
     ! matrix.  Every processor will have the same sparse matrix
     ! representation when it is finished. (Routine below)
     ! -----------------------------------------------------------------
     if (firstTime) then 
        call deallocateCSRMatrix(overlap)
        call buildGlobalSparseOverlap(overlap)
     end if

     ! -----------------------------------------------------------------
     ! Step 8: This is going to put the number of searches (coordinates)
     ! in for the costs.  Eventually we need to use the previous time,
     ! but that is going to be tricky since the sparsity structure of the
     ! overlap matrix could change....:-(
     ! -----------------------------------------------------------------

     ! Loop over by blocks and the owned cells on my blocks. This will
     ! determine which of my coordinates need to be searched on a
     ! given block.

     ! Also, since this is the last of the collective blocking
     ! communication for a while do some of the other collective comm
     ! we need to do as well. Specifically we want to tell all the
     ! processors the size of the int Buffer, real Buffer and the
     ! number of fringes we can expect to receive.

     allocate(bufSizes(nDomTotal, 6), tmpInt2D(nDomTotal, 6), &
          tmpReal(size(overlap%data)))

     ! Initialization 
     tmpReal = zero
     tmpInt2D = 0
     do nn=1, nDom
        call setPointers(nn, level, sps)
        iDom = cumDomProc(myid) + nn
        if (firstTime) then 
           do jj=overlap%rowPtr(iDom), overlap%rowPtr(iDom+1)-1
              tmpReal(jj) = real(nx*ny*nz)
           end do
        end if
        ! Sizes
        call getOBlockBufferSizes (il, jl, kl, tmpInt2D(iDom, 1), tmpInt2D(iDom, 2))
        call getOFringeBufferSizes(il, jl, kl, tmpInt2D(iDom, 3), tmpInt2D(iDom, 4))
        call getOWallBufferSizes  (il, jl, kl, tmpInt2D(iDom, 5), tmpInt2D(iDom, 6), .True.)
     end do

     if (.not. firstTime) then 
        tmpReal = overlap%data
     end if

     ! Determine the total search costs for each proc and all the bufferSizes
     call mpi_allreduce(tmpReal, overlap%data, overlap%nnz, sumb_real, MPI_SUM, &
          sumb_comm_world, ierr)
     call ECHK(ierr, __FILE__, __LINE__)

     call mpi_allreduce(tmpInt2D, bufSizes, 6*nDomTotal, sumb_integer, MPI_SUM, &
          sumb_comm_world, ierr)
     call ECHK(ierr, __FILE__, __LINE__)

     ! Done with the tmp arrays. This should be the last of the
     ! blocking collectives for a while. 
     deallocate(tmpReal, tmpInt2D)

     ! -----------------------------------------------------------------
     ! Step 8: We are now ready to partiaion and loadbalance the work
     ! based on the costs stored in the overlap matrix. These costs
     ! may be the search estimates from initializeOverlapCosts OR they
     ! may be actual timings from a previous assembly. Also create a
     ! transpose of the matrix which is useful to use for the fringe
     ! sending operation (it is transpose of the block sending)
     ! -----------------------------------------------------------------
     if (.not. lowOversetMemory) then 
        call oversetLoadBalance(overlap)
     end if
     call transposeOverlap(overlap, overlapTranspose)

     ! -----------------------------------------------------------------
     !  Step 8: Section out just the intersections we have to
     !  do. Essentially this is just the entries in the matrix that we
     !  have been assigned according to assignedProc. This will keep
     !  track of the work we have to do and if it is yet completed or
     !  not. We want this small so we can constantlly loop over it quickly. 
     !  -----------------------------------------------------------------

     nWork = 0
     do jj=1,overlap%nnz
        if (overlap%assignedProc(jj) == myid) then 
           nWork = nWork + 1
        end if
     end do
     allocate(work(4, nWork))

     nWork = 0
     do iDom=1, nDomTotal
        do jj=overlap%rowPtr(iDom), overlap%rowPtr(iDom+1)-1
           jDom = overlap%colInd(jj)
           if (overlap%assignedProc(jj) == myID) then 
              nWork = nWork + 1
              work(1, nWork) = iDom
              work(2, nWork) = jDom
              work(3, nWork) = jj
              work(4, nWork) = 0
           end if
        end do
     end do

     ! -----------------------------------------------------------------
     ! Step 8: Allocation of temporary data structures: oBlocks and fringeList
     !
     ! oBlocks: These contain the AD trees. We allocate the array to
     ! size of the total number of domains (nDomTotal). Firstly, we
     ! just add the range of domains we own (cumDomProc(myid)+1 :
     ! cumDomProc(myid)+nDom). If we need to receive a domain from
     ! another processor, we just put it in the it's global ordering
     ! slot.
     !
     ! fringeList: Similar logic for the fringe list. Allocated to
     ! size nDomTotal, add the fringes we own, and the allocate
     ! additional space for any that we need to receive from other
     ! processors. Note that we have to be really careful with the
     ! fringes: We 'isend' our owned fringes and may modify them
     ! locally as well. The MPI standard says that you cannot modify
     ! the send buffer until the isend completes. This is why we use
     ! the overMPISearchCoord type that sends only 'x', 'origQuality'
     ! and 'isWall'. These are guaranteed not to be changed by the
     ! local process while the send is on-going so we should be
     ! ok. Using the same send buffer multiple times should
     ! technically be ok according to:
     ! http://stackoverflow.com/questions/17074884/ok-to-call-mpi-isend-multiple-times-on-one-buffer
     ! the MPI 2.1 standard allows this to be ok.
     !
     ! Note that oBlocks and fringeList are allocated to size
     ! nDomTotal...which is not technically scalable, but there since
     ! there are only a few scattered variables and no large arrays it
     ! should be ok.
     ! -----------------------------------------------------------------

     allocate(oBlocks(nDomTotal), oFringes(nDomTotal), oWalls(nDomTotal))

     ! Thse variables keep track of if the block/fringes are
     ! ready. Initialized to false and only flipped when we are sure
     ! they are ready to be used. 

     allocate(oBlockReady(nDomTotal), oFringeReady(nDomTotal), oWallReady(nDomTotal))
     oBlockReady = .False.
     oFringeReady = .False.
     ! Flag the oWalls that do not actually have a wall, as being
     ! ready.
     do iDom=1,nDomtotal
        if (bufSizes(iDom, 6) == 0) then 
           oWallReady(iDom) = .True.
        else
           oWallReady(iDom) = .False.
        end if
     end do

     ! Allocate space for the localWallFringes. localWallFringes keeps
     ! track of donors for cells that are next to a wall. These must
     ! be recorded independently of the actual donors since we don't
     ! actually care what the interpolation stencil is, rather just
     ! who the donor is such that we can use that information for the
     ! flooding process. We arbitrarily set a size here and it will be
     ! automatically expanded as necessary in the fringeSearch
     ! routine.
     allocate(localWallFringes(1000))
     nLocalWallFringe = 0

     ! Determine the cells that are near wall. We have a special routine for this. 
     call flagNearWallCells(level, sps)

     do nn=1, nDom
        call setPointers(nn, level, sps)
        iDom = cumDomProc(myid) + nn

        call initializeOBlock(oBlocks(iDom), nn, level, sps)
        oBlockReady(iDom) = .True.

        call initializeOFringes(oFringes(iDom), nn)
        oFringeReady(iDom) = .True. 

        call initializeOWall(oWalls(iDom), .True., clusters(iDom))
        call packOWall(oWalls(iDom), .true.)
        oWallReady(iDom) = .True. 
     end do

     ! Call the generic routines to determine the send/receive pattern
     ! for oBlock comm and the fringe comm. These are transpose of
     ! each other. Just overestimate the sizes of the lists. 

     ! For sending, the worse case is sending all my blocks/fringes/walls to
     ! everyone but myself:
     ii = nDom*(nProc-1)
     allocate(oBlockSendList(2, ii), oFringeSendList(2, ii), oWallSendList(2, ii))

     ! For receiving, the worse receive is all the blocks/fringes/wall I
     ! don't already have:
     ii = nDomTotal - nDom
     allocate(oBlockRecvList(2, ii), oFringeRecvList(2, ii), oWallRecvList(2, ii))

     call getCommPattern(overlap, &
          oblockSendList, size(oBlockSendList, 2),  nOblockSend, &
          oBlockRecvList, size(oBlockRecvList, 2), nOblockRecv)

     call getCommPattern(overlapTranspose, &
          oFringeSendList, size(oFringeSendList, 2), nOFringeSend, &
          oFringeRecvList, size(oFringeRecvList, 2), nOFringeRecv)

     ! The wall send/recv list is essentially the merging of the
     ! oBlock and oFringe send/recv lists. Essentially if we have an
     ! oBlock OR an oFringe we need to have the oWall for it as well. 
     call getOWallCommPattern(overlap, overlapTranspose, &
          oWallSendList, size(oWallSendList, 2), nOWallSend, &
          oWallRecvList, size(oWallRecvList, 2), nOWallRecv, bufSizes(:, 6))

     ! Done with the transposed matrix
     call deallocateCSRMatrix(overlapTranspose)

     ! Zero out the overlap data since we will be doing new timings in
     ! doMyWork()
     overlap%data = zero

     ! Allocate the exact space for our send and recv requests. Note
     ! that for the oBlocks, two values are set, real and integer. 
     nn = max(nProc, &
          2*nOBlockSend + 2*nOFringeSend + 2*nOWallSend, &
          2*nOBlockRecv + 2*nOfringeRecv + 2*nOWallRecv)
     if (allocated(sendRequests)) then 
        deallocate(sendRequests, recvRequests)
     end if

     allocate(sendRequests(nn), recvRequests(nn))
     allocate(recvInfo(2, nn))

     ! Post all the oBlock/oFringe/oWall iSends
     sendCount = 0
     do jj=1, nOblockSend
        iProc = oBlockSendList(1, jj)
        iDom = oBlockSendList(2, jj)
        call packOBlock(oBlocks(iDom))
        call sendOBlock(oBlocks(iDom), iDom, iProc, 0, sendCount)
     end do

     do jj=1, nOFringeSend
        iProc = oFringeSendList(1, jj)
        iDom = oFringeSendList(2, jj)
        call packOFringe(oFringes(iDom))
        call sendOFringe(oFringes(iDom), iDom, iProc, MAGIC, sendCount)
     end do

     do jj=1, nOWallSend
        iProc = oWallSendList(1, jj)
        iDom = oWallSendList(2, jj)
        call sendOWall(oWalls(iDom), iDom, iProc, 2*MAGIC, sendCount)
     end do

     ! Post all the oBlock/oFringe/oWall receives. Before posting the actual
     ! receive, allocate the receiving buffer. 
     recvCount = 0
     do jj=1, nOBlockRecv
        iProc = oBlockRecvList(1, jj)
        iDom = oBlockRecvList(2, jj)
        call recvOBlock(oBlocks(iDom), iDom, iProc, 0, &
             bufSizes(iDom, 1), bufSizes(iDom, 2), recvCount, recvInfo)
     end do

     do jj=1, nOFringeRecv
        iProc = oFringeRecvList(1, jj)
        iDom = oFringeRecvList(2, jj)
        call recvOFringe(oFringes(iDom), iDom, iProc, MAGIC, &
             bufSizes(iDom, 3), bufSizes(iDom, 4), recvCount, recvInfo)
     end do

     do jj=1, nOWallRecv
        iProc = oWallRecvList(1, jj)
        iDom = oWallRecvList(2, jj)
        call recvOWall(oWalls(iDom), iDom, iProc, 2*MAGIC, &
             bufSizes(iDom, 5), bufSizes(iDom, 6), recvCount, recvInfo)
     end do


     ! Before we start waiting for the receives to finish, we can see
     ! if we can do any searches with the blocks/fringes we already
     ! have. Call the internal routine for this.
     call doMyWork(flag)

     ! Complete all the recives
     do i=1, recvCount

        ! Complete any one of the recv requests
        call mpi_waitany(recvCount, recvRequests, index, status, ierr)
        call ECHK(ierr, __FILE__, __LINE__)

        ! Global domain index of the recv that finished
        iDom = recvInfo(1, index)

        ! Check which type of receive just finished and flag them as
        ! being complete.
        if     (recvInfo(2, index) == 1) then 
           oBlocks(iDom)%realBufferReady = .True. 
        else if (recvInfo(2, index) == 2) then 
           oBlocks(iDom)%intBufferReady = .True. 
        else if (recvInfo(2, index) == 3) then 
           oFringes(iDom)%realBufferReady = .True. 
        else if (recvInfo(2, index) == 4) then 
           oFringes(iDOm)%intBufferReady = .True. 
        else if (recvInfo(2, index) == 5) then 
           oWalls(iDom)%realBufferReady = .True. 
        else if (recvInfo(2, index) == 6) then 
           oWalls(iDom)%intBufferReady = .True. 
        end if

        ! If both int and real buffers are received, we can unpack the
        ! oblock and flag it as ready.
        if (oBlocks(iDom)%realBufferReady .and. oBlocks(iDom)%intBufferReady .and. &
             .not.oBlocks(iDom)%allocated) then 
           call unpackOBlock(oBlocks(iDom))           
           oBlockReady(iDom) = .True.
        end if

        ! If both int and real buffers are received, we can unpack the
        ! oFringe and flag it as ready.
        if (oFringes(iDom)%realBufferReady .and. oFringes(iDom)%intBufferReady .and. &
             .not.oFringes(iDom)%allocated) then 
           call unpackOFringe(oFringes(iDom))
           oFringeReady(iDom) = .True.
        end if

        ! If both int and real buffers are received, we can unpack the
        ! oWall and flag it as ready.
        if (oWalls(iDom)%realBufferReady .and. oWalls(iDom)%intBufferReady .and. &
             .not.oWalls(iDom)%allocated) then 
           call unpackOWall(oWalls(iDom))
           oWallReady(iDom) = .True.
        end if

        ! Now see if we can do any more of the work, ie the searches. 
        call doMyWork(flag)

        ! Sanity check. flag better be true when i=recvCount
        if (i==recvCount .and. .not. flag) then 
           call terminate("computeInterpolationParallel", "Inconsistent Comm pattern detected.")
        end if
     end do

     ! Last thing to do wait for all the sends to finish 
     do i=1,sendCount
        call mpi_waitany(sendCount, sendRequests, index, status, ierr)
        call ECHK(ierr, __FILE__, __LINE__)
     end do

     ! We are now completely finished with oBlocks and owalls so
     ! delete before we allocate space for all the fringes
     call deallocateOBlocks(oBlocks, size(oBlocks))
     call deallocateOWalls(oWalls, size(oWalls))
     deallocate(oBlocks, oWalls)

     ! Make sure all oFringe buffers are delloacted before we allocate
     ! space for the large fringe arraay
     do iDom=1, nDomTotal
        if (allocated(oFringes(iDom)%iBuffer)) then 
           deallocate(oFringes(iDom)%iBuffer, &
                oFringes(iDom)%rBuffer)
        end if
     end do
 
     ! Now create the fringes
     do nn=1, nDom
        call setPointers(nn, level, sps)
        call initializeFringes(nn, level, sps)
     end do

     ! -----------------------------------------------------------------
     ! Step 9: Well, all the searches are done, so now we can now send
     ! the fringes back to where they came from. However, since we
     ! have a *lot* of data to send back, we won't send it all
     ! back. Instead, we will go through the oFringes we we've had to
     ! deal with and "prune" them...that is make a compacted list with
     ! just the ones that are not useless. This way we have to send
     ! less data and the receiving process has to deal with less data
     ! as well. We will be a little sneaky here since we will put the
     ! values into the i/r Buffers at the same time so don't need to
     ! make a separate pass. Essentially this *is* the packing routine.
     ! -----------------------------------------------------------------

     do iDom=1, nDomTotal
        if (oFringes(iDom)%allocated) then 
           ! Fringe is allocated so check it
           iPtr = 0
           rPtr = 0
           ! First pass count up the sizes
           do i=1, size(oFringes(iDom)%donorProc)
              if (oFringes(iDom)%donorProc(i) /= -1) then 
                 iPtr = iPtr + 14
                 rPtr = rPtr + 4
              end if
           end do

           allocate(oFringes(iDom)%iBuffer(iPtr), oFringes(iDom)%rBuffer(rPtr))
           iPtr = 0
           rPtr = 0
           ! Second pass add the values
           do i=1, size(oFringes(iDom)%donorProc)
              if (oFringes(iDom)%donorProc(i) /= -1) then 
                 
                 ! Copy the values for this fringe into the fringe's
                 ! buffer
                 oFringes(iDom)%iBuffer(iPtr+1) = oFringes(iDom)%donorProc(i)
                 oFringes(iDom)%iBuffer(iPtr+2) = oFringes(iDom)%donorBlock(i)
                 oFringes(iDom)%iBuffer(iPtr+3) = oFringes(iDom)%dI(i)
                 oFringes(iDom)%iBuffer(iPtr+4) = oFringes(iDom)%dJ(i)
                 oFringes(iDom)%iBuffer(iPtr+5) = oFringes(iDom)%dK(i)
                 oFringes(iDom)%iBuffer(iPtr+6) = oFringes(iDom)%gInd(1, i)
                 oFringes(iDom)%iBuffer(iPtr+7) = oFringes(iDom)%gInd(2, i)
                 oFringes(iDom)%iBuffer(iPtr+8) = oFringes(iDom)%gInd(3, i)
                 oFringes(iDom)%iBuffer(iPtr+9) = oFringes(iDom)%gInd(4, i)
                 oFringes(iDom)%iBuffer(iPtr+10) = oFringes(iDom)%gInd(5, i)
                 oFringes(iDom)%iBuffer(iPtr+11) = oFringes(iDom)%gInd(6, i)
                 oFringes(iDom)%iBuffer(iPtr+12) = oFringes(iDom)%gInd(7, i)
                 oFringes(iDom)%iBuffer(iPtr+13) = oFringes(iDom)%gInd(8, i)
                 oFringes(iDom)%iBuffer(iPtr+14) = oFringes(iDom)%myIndex(i)
                 iPtr = iPtr + 14

                 oFringes(iDom)%rBuffer(rPtr+1) = oFringes(iDom)%donorFrac(1, i)
                 oFringes(iDom)%rBuffer(rPtr+2) = oFringes(iDom)%donorFrac(2, i)
                 oFringes(iDom)%rBuffer(rPtr+3) = oFringes(iDom)%donorFrac(3, i)
                 oFringes(iDom)%rBuffer(rPtr+4) = oFringes(iDom)%quality(i)

                 rPtr = rPtr + 4
              end if
           end do
           oFringes(iDom)%fringeReturnSize = rPtr/4

        end if
     end do

     ! -----------------------------------------------------------------
     ! For this data exchange we use the exact *reverse* of fringe
     ! communication pattern from the previous data exchange. We
     ! actually do two exchanges: The purpose of the first exchange is
     ! to just communicate the sizes. Then the receiving processors
     ! can allocated sufficent buffer space for the incoming fringe
     ! information. This is necessary since we want to use a
     ! non-blocking receive and we don't what to do a collective comm
     ! here to get the sizes. So the only way is do another
     ! point-to-point.
     ! -----------------------------------------------------------------

     ! Post all the fringe iSends
     sendCount = 0
     do jj=1, nOFringeRecv

        iProc = oFringeRecvList(1, jj)
        iDom = oFringeRecvList(2, jj)

        sendCount = sendCount + 1
        call mpi_isend(oFringes(iDom)%fringeReturnSize, 1, sumb_integer, &
             iproc, iDom, sumb_comm_world, sendRequests(sendCount), ierr)
        call ECHK(ierr, __FILE__, __LINE__)
     end do

     allocate(fringeRecvSizes(nOfringeSend))

     ! Non-blocking receives
     recvCount = 0
     do jj=1, nOFringeSend

        iProc = oFringeSendList(1, jj)
        iDom = oFringeSendList(2, jj)
        recvCount = recvCount + 1

        call mpi_irecv(fringeRecvSizes(jj), 1, sumb_integer, &
             iProc, iDom, sumb_comm_world, recvRequests(recvCount), ierr)
        call ECHK(ierr, __FILE__, __LINE__)
     end do

     ! Last thing to do wait for all the sends and receives to finish 
     ! Last thing to do wait for all the sends to finish 
     do i=1,sendCount
        call mpi_waitany(sendCount, sendRequests, index, status, ierr)
        call ECHK(ierr, __FILE__, __LINE__)
     end do

     do i=1,recvCount
        call mpi_waitany(recvCount, recvRequests, index, status, ierr)
        call ECHK(ierr, __FILE__, __LINE__)
     end do

     ! Now before we do the actual receives, before we need to
     ! allocate space intRecvBuff and realRecvBuff for the receive. We
     ! also compute the cumulative offsets so that we know what to
     ! check when a particular receive completes.

     allocate(cumFringeRecv(1:nOFringeSend+1))
     cumFringeRecv(1) = 1
     do jj=1, nOFringeSend ! These are the fringes we *sent*
        ! originally, now are going to receive them
        ! back
        cumFringeRecv(jj+1) = cumFringeRecv(jj) + fringeRecvSizes(jj)
     end do

     ! Now alocate the integer space. Note we are receiving 4 real
     ! values and 14 int values:
     ii = cumFringeRecv(nOfringeSend+1)-1
     allocate(intRecvBuf(ii*14), realRecvBuf(ii*4))   

     ! We are now ready to actually receive our fringes
     sendCount = 0
     do jj=1, nOFringeRecv

        iProc = oFringeRecvList(1, jj)
        iDom = oFringeRecvList(2, jj)
        iSize = oFringes(iDom)%fringeReturnSize 
        if (iSize > 0) then 
           tag = iDom + MAGIC
           sendCount = sendCount + 1
           call mpi_isend(oFringes(iDom)%rBuffer, iSize*4, sumb_real, &
                iproc, tag, sumb_comm_world, sendRequests(sendCount), ierr)
           call ECHK(ierr, __FILE__, __LINE__)

           tag = iDom + 2*MAGIC
           sendCount = sendCount + 1
           call mpi_isend(oFringes(iDom)%iBuffer, iSize*14, sumb_integer, &
                iproc, tag, sumb_comm_world, sendRequests(sendCount), ierr)
           call ECHK(ierr, __FILE__, __LINE__)
        end if
     end do

     ! Non-blocking receives
     recvCount = 0
     do jj=1, nOfringeSend

        iProc = oFringeSendList(1, jj)
        iDom = oFringeSendList(2, jj)
        iSize = cumFringeRecv(jj+1) - cumFringeRecv(jj) 
        if (iSize > 0) then 

           iStart = (cumFringeRecv(jj  )-1)*4 + 1
           tag = iDom + MAGIC
           recvCount = recvCount + 1       
           call mpi_irecv(realRecvBuf(iStart), iSize*4, sumb_real, &
                iProc, tag, sumb_comm_world, recvRequests(recvCount), ierr)
           call ECHK(ierr, __FILE__, __LINE__)
           recvInfo(:, recvCount) = (/iDom, 1/) ! 1 for real recv

           iStart = (cumFringeRecv(jj  )-1)*14 + 1
           tag = iDom + 2*MAGIC
           recvCount = recvCount + 1                
           call mpi_irecv(intRecvBuf(iStart), iSize*14, sumb_integer, &
                iProc, tag, sumb_comm_world, recvRequests(recvCount), ierr)
           call ECHK(ierr, __FILE__, __LINE__)
           recvInfo(:, recvCount) = (/iDom, 2/) ! 2 for int recv
        end if
     end do

     ! We can do some useful work while the fringes are
     ! communicating. Specifically we can process the local
     ! fringes. This is essentially the operation as below we perform
     ! with the incoming data from other processors. Keep a running
     ! count of the number of *acutal* local fringes this processor
     ! owns. 
     nLocalFringe = 0

     do nn=1,nDom
        call setPointers(nn, level, sps)

        iDom = cumDomProc(myid) + nn

        ! We can cheat here and just do a nice triple loop, this is
        ! because these fringes are local and we know we still have
        ! all of them and are still in the "right order"
        ii = 0
        do k=2, kl
           do j=2, jl
              do i=2, il
                 ii =ii + 1

                 if (oFringes(iDom)%donorProc(ii) /= -1) then 

                    nLocalFringe = nLocalFringe + 1

                    ! We have a donor to use:
                    fringes(i, j, k)%donorProc  = oFringes(iDom)%donorProc(ii)
                    fringes(i, j, k)%donorBlock = oFringes(iDom)%donorBlock(ii)
                    fringes(i, j, k)%donorFrac = oFringes(iDom)%donorFrac(:, ii)
                    fringes(i, j, k)%dI = oFringes(iDom)%dI(ii)
                    fringes(i, j, k)%dJ = oFringes(iDom)%dJ(ii)
                    fringes(i, j, k)%dK = oFringes(iDom)%dK(ii)
                    fringes(i, j, k)%gInd(1) = oFringes(iDom)%gInd(1, ii)
                    fringes(i, j, k)%gInd(2) = oFringes(iDom)%gInd(2, ii)
                    fringes(i, j, k)%gInd(3) = oFringes(iDom)%gInd(3, ii)
                    fringes(i, j, k)%gInd(4) = oFringes(iDom)%gInd(4, ii)
                    fringes(i, j, k)%gInd(5) = oFringes(iDom)%gInd(5, ii)
                    fringes(i, j, k)%gInd(6) = oFringes(iDom)%gInd(6, ii)
                    fringes(i, j, k)%gInd(7) = oFringes(iDom)%gInd(7, ii)
                    fringes(i, j, k)%gInd(8) = oFringes(iDom)%gInd(8, ii)

                    ! Now unwind the index of *donor*

                    ! Remove the compute status of this cell
                    call setIsCompute(fringes(i, j, k)%status, .False.)

                    ! Very Important --- also store the quality of the
                    ! donor...we need to compare this with the quality
                    ! of the more potential donors comming in on the comm.
                    fringes(i, j, k)%quality = oFringes(iDom)%quality(ii)

                 end if
              end do
           end do
        end do
     end do

     ! Now wait for the sends and receives to finish
     do i=1,sendCount
        call mpi_waitany(sendCount, sendRequests, index, status, ierr)
        call ECHK(ierr, __FILE__, __LINE__)
     end do

     do i=1,recvCount
        call mpi_waitany(recvCount, recvRequests, index, status, ierr)
        call ECHK(ierr, __FILE__, __LINE__)
     end do

     ! Process the data we just received. 
     do kk=1, nOfringeSend

        ! Local block index of the fringes
        iDom = oFringeSendList(2, kk)
        nn = iDom - cumDomProc(myid)

        ! Set the block pointers for the local block we are dealing
        ! with:
        call setPointers(nn, level, sps)

        ! This is the range of fringes that are now ready. 
        do jj=cumFringeRecv(kk), cumFringeRecv(kk+1)-1

           ! We need to unwind *my* index. The reason why we sent
           ! the index in the first place is we are not getting the
           ! same number of fringes back as we sent so the myIndex
           ! lets of know which ones are actuall coming back. 

           ! myindex is 1 based so we need the -1 at the end 
           myIndex = intRecvBuf(14*(jj-1) + 14) - 1
           i = mod(myIndex, nx) + 2
           j = mod(myIndex/nx, ny) + 2
           k = myIndex/(nx*ny) + 2

           ! Extract the quality value from the buffer
           quality = realRecvBuf(4*(jj-1) + 4)

           ! This is the acutal implict hole cutting "less than"
           ! operation. 
           if (quality < fringes(i, j, k)%quality) then 

              ! Only count this a new local fringe if it doesn't
              ! already have one
              if (fringes(i, j, k)%donorProc == -1) then 
                 nLocalFringe = nLocalFringe + 1
              end if

              ! Accept the incoming fringe. 
              iStart = 14*(jj-1)
              fringes(i, j, k)%donorProc  = intRecvBuf(iStart + 1)
              fringes(i, j, k)%donorBlock = intRecvBuf(iStart + 2)
              fringes(i, j, k)%dI =         intRecvBuf(iStart + 3)
              fringes(i, j, k)%dJ =         intRecvBuf(iStart + 4)
              fringes(i, j, k)%dK =         intRecvBuf(iStart + 5)

              fringes(i, j, k)%gInd(1)    = intRecvBuf(iStart + 6)
              fringes(i, j, k)%gInd(2)    = intRecvBuf(iStart + 7)
              fringes(i, j, k)%gInd(3)    = intRecvBuf(iStart + 8)
              fringes(i, j, k)%gInd(4)    = intRecvBuf(iStart + 9)
              fringes(i, j, k)%gInd(5)    = intRecvBuf(iStart + 10)
              fringes(i, j, k)%gInd(6)    = intRecvBuf(iStart + 11)
              fringes(i, j, k)%gInd(7)    = intRecvBuf(iStart + 12)
              fringes(i, j, k)%gInd(8)    = intRecvBuf(iStart + 13)

              fringes(i, j, k)%donorFrac = realRecvBuf(4*jj-3:4*jj-1)

              ! Set this new quality
              fringes(i, j, k)%quality = quality

              ! Remove the compute status of this cell
              call setIsCompute(fringes(i, j, k)%status, .False.)

           end if
        end do
     end do

     ! ------------------------------------------------------------------
     ! We are now completely finished with oFringes and the buffers
     call deallocateOFringes(oFringes, size(oFringes))
     deallocate(oFringes, intRecvBuf, realRecvBuf)

     ! -----------------------------------------------------------------
     ! Step 9: We now have computed all the fringes that we can. Some of
     ! them may be 'irrigular' as described in "A highly automated
     ! parallel Chimera method for overset grids based on the implicit
     ! hole cutting technique". If a cell is both a fringe and a donor,
     ! we need to force it to be a compute cell and remove it's receiver
     ! status. Do this as long as it isn't a forced receiver...we can't
     ! touch the forced receivers. 
     !
     ! How we do this is as follows: On each processor we gather up the
     ! fringes on each block on each processor that are actual fringes,
     ! ie, they have donorProc /= -1. That is, there is a
     ! potential donor that is better than its own cell. Once we get a
     ! list of these, we sort them, by processor, then block, then
     ! index. Then we can send this fringe information back to the
     ! processors where the donors came from.
     ! -----------------------------------------------------------------

     ! Allocate the a new local 1D fringe list that just has our local
     ! fringes (from all local blocks) that are *actually* fringes. Do
     ! not include the halos.
     allocate(localFringes(nLocalFringe))


     ! Fill up these fringes
     nLocalFringe = 0
     do nn=1, nDom
        call setPointers(nn, level, sps)
        do k=2, kl
           do j=2, jl
              do i=2, il

                 ! Check if this cell is a fringe:
                 if (fringes(i, j, k)%donorProc /= -1) then 

                    ! Now check if this cell *really* needs a
                    ! donor...it all its neighbours are also
                    ! interpolated it will get blanked so we can just
                    ! forget about it. 
                    computeCellFound = .False.
                    stencilLoop2: do i_stencil=1, N_visc_drdw
                       ii = visc_drdw_stencil(i_stencil, 1) + i
                       jj = visc_drdw_stencil(i_stencil, 2) + j
                       kk = visc_drdw_stencil(i_stencil, 3) + k

                       if (isCompute(fringes(ii, jj, kk)%status)) then 
                          ! This is a compute cell
                          computeCellFound = .True.
                       end if
                    end do stencilLoop2

                    if (computeCellFound) then 
                       nLocalFringe = nLocalFringe + 1
                       localFringes(nLocalFringe) = fringes(i, j, k)
                    end if
                 end if
              end do
           end do
        end do
     end do

     call determineDonors(level, sps, localFringes, nLocalFringe, .False.)

     !==================================================================================

     ! -----------------------------------------------------------------
     ! Step 10: We also have to send the wall fringes to their
     ! respective donor procs. We use a very similar communication
     ! structure as we used for the regular fringes. The difference
     ! here is that we unlike the regular case we keep *all* fringes
     ! we receive and just put them in a big list:
     ! wallFringes. wallFringes is the list of fringes ON THE DONOR
     ! processor whos receiver is a wall point. This is precisely the
     ! information we need to start the flooding process later on. 
     ! -----------------------------------------------------------------

     call determineDonors(level, sps, localWallFringes, nLocalWallFringe, .True.)

     ! !=================================================================================
     ! ! -----------------------------------------------------------------
     ! ! Step 10: We can now locally perform the irregular cell correction
     ! ! by looping over the fringes on my proc and just checking if
     ! ! donorProc is not -1 and isDonor is True.  If so, we force it back
     ! ! to be compute, by cancelling the donor information. Update the
     ! ! fringes when we're done so everyone has up to date information.
     ! ! -----------------------------------------------------------------

     call exchangeStatusTranspose(level, sps, commPatternCell_2nd, internalCell_2nd)
     call irregularCellCorrection(level, sps)

     ! Next we have to perfrom the interior cell flooding. We already
     ! have the information we need: we have isWallFringe defined in
     ! the fringes as well as knowing if a cell is a compute. We
     ! should probably only flood compute cells that are not also
     ! donors, since that would get a little complicated. 

     call floodInteriorCells(level, sps)


     ! The fringeReduction just needs to be isCompute flag so exchange
     ! the status this as these may have been changed by the flooding

     call exchangeStatus(level, sps, commPatternCell_2nd, internalCell_2nd)
     !-----------------------------------------------------------------
     ! Step 15: Reduction of the number of fringes. What we do is look at
     ! all the fringes and see if all the cells in its stencil are also
     ! fringes or holes. If so we can flag that particular cell as a
     ! hole.
     ! -----------------------------------------------------------------

     call fringeReduction(level, sps)

     ! Before we can do the final comm structures, we need to make
     ! sure that every processor's halo have any donor information
     ! necessary to build its own comm pattern. For this will need to
     ! send donorProc, donorBlock, dI, dJ, dK and donorFrac. 

     call exchangeFringes(level, sps, commPatternCell_2nd, internalCell_2nd)

     ! -----------------------------------------------------------------
     ! Step 17: We can now create the final required comm structures
     ! based on our interpolation. This is relatively straight forward:
     ! The fringes we have local are the cells the "receiver" cells. We
     ! sort those by donor processor (as we did before) and then send
     ! them to donor proc. This then form the sending information. The
     ! internal copy is formed from the part that is on-processor. 
     ! -----------------------------------------------------------------

     call finalOversetCommStructures(level, sps)

     ! VERY last thing is to update iBlank based on the status of our local fringes. 
     call setIblankArray(level, sps)

     ! -----------------------------------------------------------------
     ! Step 16: The algorithm is now complete. Run the checkOverset
     ! algorithm to verify that we actually have a valid interpolation
     ! -----------------------------------------------------------------
     call checkOverset(level, sps)
     
     ! -----------------------------------------------------------------
     ! Step 18: Create the zipper mesh. We pass in a few arrays
     ! dealing with wall exchange since there is no need to recompute them. 
     ! -----------------------------------------------------------------
     call createZipperMesh(level, sps, oWallSendList, oWallRecvList, &
          nOwallSend, nOwallRecv, size(oWallSendList, 2), &
          size(oWallRecvList, 2), work, nWork)

     ! Setup the buffer sizes
     call setBufferSizes(level, sps, .false., .false., .true.)

     ! Deallocate some data we no longer need
     deallocate(Xmin, Xmax, oBlockReady, oFringeReady, oWallReady,  work)

     ! Done with the clusters
     deallocate(clusters)

  end do spectralLoop



  ! Free the buffer and make a new one that includes necessary sizes
  ! for the overset comm
  deallocate(sendBuffer, recvBuffer)
  allocate(sendBuffer(sendBufferSize), recvBuffer(recvBufferSize))
  deallocate(cumdomproc, ndomproc)
contains

  ! Simple utility-type routines that make the main subroutine
  ! easier to read

  subroutine blockSizeInfo

    implicit none

    ! Gather the dimensions of all blocks to everyone
    call mpi_allreduce(nDom, nDomTotal, 1, sumb_integer, MPI_SUM, &
         sumb_comm_world, ierr)
    call ECHK(ierr, __FILE__, __LINE__)

    ! Receive the number of domains from each proc using an allgather.
    call mpi_allgather(nDom, 1, sumb_integer, nDomProc, 1, sumb_integer, &
         sumb_comm_world, ierr)
    call ECHK(ierr, __FILE__, __LINE__)

    ! Compute the cumulative format:
    cumDomProc(0) = 0
    do iProc=1, nProc
       cumDomProc(iProc) = cumDomProc(iProc-1) + nDomProc(iProc-1)
    end do

  end subroutine blockSizeInfo

  subroutine computeDomainBoundingBoxes

    implicit none

    ! Working Variables
    real(kind=realType), dimension(3, nDom) :: xMinLocal, xMaxLocal

    do nn=1,nDom
       call setPointers(nn, level, sps)

       xMinLocal(1, nn) = minval(x(:, :, :, 1))
       xMinLocal(2, nn) = minval(x(:, :, :, 2))
       xMinLocal(3, nn) = minval(x(:, :, :, 3))

       xMaxLocal(1, nn) = maxval(x(:, :, :, 1))
       xMaxLocal(2, nn) = maxval(x(:, :, :, 2))
       xMaxLocal(3, nn) = maxval(x(:, :, :, 3))

    end do

    ! Now we can allgather the xMin and xMax  from each
    ! processor to everyone
    call mpi_allgatherV(xMinLocal, nDom*3, sumb_real, xMin, 3*nDomProc, &
         3*cumDomProc, sumb_real, sumb_comm_world, ierr)
    call ECHK(ierr, __FILE__, __LINE__)

    call mpi_allgatherV(xMaxLocal, nDom*3, sumb_real, xMax, 3*nDomProc, &
         3*cumDomProc, sumb_real, sumb_comm_world, ierr)
    call ECHK(ierr, __FILE__, __LINE__)

  end subroutine computeDomainBoundingBoxes

  subroutine buildGLobalSparseOverlap(overlap)

    use overset
    implicit none

    ! Input/Output
    type(CSRMatrix), intent(inout) :: overlap

    ! Working
    integer(kind=intType), dimension(:), allocatable :: colIndLocal, rowPtrLocal
    logical , dimension(:, :), allocatable :: localOverlap
    integer(kind=intType) :: nnzLocal
    integer(kind=intType), dimension(:), allocatable :: nnzProc, cumNNzProc

    ! Allocate the space for my processor's compoent of iOverlap. The
    ! number of rows is the number of blocks I own (in my own numbering)
    ! and the coluns are the total number of blocks. It is initialized
    ! to True, and remains that way unless we can conclusively prove the
    ! blocks don't overlap. 
    allocate(localOverlap(nDom, nDomTotal))
    localOverlap = .True.

    ! Assume all rows are dense, and decrement when we eliminate an
    ! entry
    nnzLocal = nDom*nDomTotal

    ! Bounding box check: We will check our *own* blocks against the
    ! bounding boxes for for each domain.  This will fill up the
    ! localOverlapArray and determine nnzLocal -- the number of
    ! non-zero local entries in the global sparse matrix.

    do nn=1, nDom
       iDom = cumDomProc(myid) + nn

       ! Now Loop over *all* of the other blocks
       do jDom=1, nDomTotal 

          ! We can eliminate some of pairs using the cluser
          ! analysis. Note that the cluster analysis will take care
          ! of making sure that a block doesn't interset itself. ie,
          ! the diagonal of the full matrix must be .False.

          if (clusters(iDom) == clusters(jDom)) then 
             localOverlap(nn, jDom) = .False.
             nnzLocal = nnzLocal - 1
          end if

          ! Only do the spatial check if we haven't elminated the
          ! connection through the cluster check
          if (localOverlap(nn, jDom)) then 

             ! Now do the box overlap check
             if ( &
                  xMin(1, iDom) >= xMax(1, jDom) .or. &
                  xMax(1, iDom) <= xMin(1, jDom) .or. &
                  xMin(2, iDom) >= xMax(2, jDom) .or. &
                  xMax(2, iDom) <= xMin(2, jDom) .or. &     
                  xMin(3, iDom) >= xMax(3, jDom) .or. &
                  xMax(3, iDom) <= xMin(3, jDom)) then 

                ! These bounding boxes do not intersect. 
                localOverlap(nn, jDom) = .False.
                nnzLocal = nnzLocal -1 
             end if
          end if
       end do
    end do

    ! Now, create a sparse matrix representation of the local part

    allocate(colIndLocal(nnzLocal), rowPtrLocal(nDom + 1))
    rowPtrLocal(1) = 1
    i = 0
    do nn=1, nDom
       do jDom=1,nDomTotal
          if (localOverlap(nn, jDom)) then 
             i = i + 1
             colIndLocal(i) = jDom
          end if
       end do
       rowPtrLocal(nn+1)  = i + 1 
    end do

    ! Now we want to assemble the global (sparse!) connectivity
    ! matrix for all processors. This is going to require a little
    ! communication of sizes first followed by the actual data.

    ! Determine distribution of non-zero locations
    allocate(nnzProc(nProc), cumNNZProc(0:nProc))
    call mpi_allgather(nnzLocal, 1, sumb_integer, nnzProc, 1, sumb_integer, &
         sumb_comm_world, ierr)

    overlap%nnz = sum(nnzProc)
    overlap%nRow = nDomTotal
    overlap%nCol = nDomTotal
    overlap%nnzLocal = nnzLocal
    ! We can now know how big the global data will be
    allocate(&
         overlap%data(overlap%nnz), &
         overlap%colInd(overlap%nnz), &
         overlap%rowPtr(overlap%nRow + 1), &
         overlap%assignedProc(overlap%nnz))

    overlap%allocated = .True.
    cumNNZProc(0) = 0
    do iProc=1,nProc
       cumNNZProc(iProc) = cumNNZProc(iProc-1) + nnzProc(iProc)
    end do

    ! Correct for rowPtrLocal for the cumulative number of nnz up to
    ! this proc:
    rowPtrLocal = rowPtrLocal + cumNNZProc(myid)

    ! Gather the column indicies
    call mpi_allgatherV(colIndLocal, nnzLocal, sumb_integer, overlap%colInd, &
         nnzProc, cumNNZProc, sumb_integer, sumb_comm_world, ierr)

    ! Now we gather the rowPtr to everone
    overlap%rowPtr(1) = 1
    call mpi_allgatherV(rowPtrLocal(2:nDom), nDom, sumb_integer, &
         overlap%rowPtr(2:nDomTotal+1), &
         nDomProc, cumDomProc, sumb_integer, sumb_comm_world, ierr)

    ! Initialize the assignedProc to the owned rows of a processor
    do iProc=0,nProc-1
        overlap%assignedProc(cumNNZProc(iProc)+1:cumNNZProc(iProc+1)) = iProc
     end do

    deallocate(colIndLocal, rowPtrLocal, localOverlap, cumNNzProc)

  end subroutine buildGlobalSparseOverlap

  subroutine doMyWork(flag)

    ! This internal subroutine which may be called repeadly, performs
    ! as many of the searches my processor is responsible for. It
    ! returns true when all the work has been completed. 

    implicit none
    logical, intent(out) :: flag
    integer(kind=intType) :: iDom, jDom, jj

    flag = .True.
    do iWork=1, nWork

       iDom = work(1, iWork)
       jDom = work(2, iWork)
       jj   = work(3, iWork) ! Index from the overlap matrix

       ! Check if I have the oBlock and fringes i need to do this
       ! intersection and I haven't already done it.
       if (oBlockReady(iDom) .and. oFringeReady(jDom) .and. &
            oWallReady(iDom) .and. oWallReady(jDom) .and. &
            work(4, iWork) == 0) then 

          startTime = mpi_wtime()
          call fringeSearch(oBlocks(iDom), oFringes(jDom), oWalls(iDom), oWalls(jDom))
          endTime = mpi_wtime()
          overlap%data(jj) = endTime - startTime

          ! Flag this work as being done:
          work(4, iWork) = 1
       end if

       ! If any one is not done, flag is flipped to False.
       if (work(4, iWork) /= 1) then 
          flag = .False. 
       end if
    end do

  end subroutine doMyWork

end subroutine oversetComm

! ----------------------------------------------------------------
! Other subroutines
! ----------------------------------------------------------------



subroutine test
  use overset
  use blockPointers
  implicit none

  integer(kind=intType) :: nn, i, j, k, mm
  real(kind=realType) :: val1, val2, coor(3)

  do nn=1,1

     call setPointers(nn, 1, 1)
     do k=1, ke
        do j=1, je
           do i=1, ie

              ! Get cell center coordinate
              coor = eighth*(&
                   x(i-1, j-1, k-1, :) + &
                   x(i  , j-1, k-1, :) + &
                   x(i-1, j  , k-1, :) + &
                   x(i  , j  , k-1, :) + &
                   x(i-1, j-1, k  , :) + &
                   x(i  , j-1, k  , :) + &
                   x(i-1, j  , k  , :) + &
                   x(i  , j  , k  , :))
              w(i, j, k, iVx) = coor(1) + 2*coor(2) + 3*coor(3)
           end do
        end do
     end do
  end do

  call setPointers(2, 1, 1)
  w(:, :, :, iVx) = 0.0

  call whalo2(1, 1, 5, .False., .False. , .False. )

end subroutine test


subroutine writePartionedMesh(fileName)

  ! This is a debugging routine for writing out meshes *as they are
  ! partioned*. This can be useful for debugging overset issues.

  use communication
  use blockPointers

  implicit none

  character(len=*), intent(in) :: fileName
  integer(kind=intType) :: nDomTotal, iProc, nn, i, j, k, iDim, iDom, ierr, ii
  integer(kind=intType) :: bufSize, maxSize, ibufSize, imaxSize
  integer(kind=intType), dimension(3, nDom) :: localDim
  integer(kind=intType), dimension(:), allocatable :: nDomProc, cumDomProc
  integer(kind=intType), dimension(:, :), allocatable :: dims
  real(kind=realType), dimension(:), allocatable :: buffer
  integer(kind=intType), dimension(:), allocatable :: ibuffer !iblank
  character*40 :: tmpStr

  integer status(MPI_STATUS_SIZE) 
  ! Gather the dimensions of all blocks to everyone
  call mpi_allreduce(nDom, nDomTotal, 1, sumb_integer, MPI_SUM, &
       sumb_comm_world, ierr)
  call ECHK(ierr, __FILE__, __LINE__)

  ! Store the sizes of the local blocks
  do nn=1,nDom
     call setPointers(nn, 1, 1)

     ! Store the 'l' sizes for transferring
     localDim(1, nn) = il
     localDim(2, nn) = jl
     localDim(3, nn) = kl
  end do

  ! Allocate the space we need for the numbers and cumulative form
  allocate(nDomProc(0:nProc-1), cumDomProc(0:nProc), dims(3, nDomTotal))

  ! Receive the number of domains from each proc using an allgather.
  call mpi_allgather(nDom, 1, sumb_integer, nDomProc, 1, sumb_integer, &
       sumb_comm_world, ierr)
  call ECHK(ierr, __FILE__, __LINE__)

  ! Compute the cumulative format:
  cumDomProc(0) = 0
  do iProc=1, nProc
     cumDomProc(iProc) = cumDomProc(iProc-1) + nDomProc(iProc-1)
  end do

  ! We will also allgather all of the block sizes which will make
  ! things a little easier since everyone will know the proper sizes
  ! for the sends
  call mpi_allgatherV(localDim, nDom*3, sumb_integer, dims, 3*nDomProc, &
       3*cumDomProc, sumb_integer, sumb_comm_world, ierr)
  call ECHK(ierr, __FILE__, __LINE__)

  maxSize = 0
  imaxSize = 0
  do i=1,nDomTotal
     maxSize = max(maxSize, dims(1, i)*dims(2,i)*dims(3,i)*3)
     imaxSize = max(imaxSize, dims(1, i)*dims(2,i)*dims(3,i))
  end do

  allocate(buffer(maxSize))
  allocate(ibuffer(imaxSize))

  if (myid == 0) then 
     print *,'writing mesh...'
     ! Root proc does all the writing. Just dump to ascii tecplot
     ! file---really slow.

     open(unit=1, file=fileName, form='formatted', status='unknown')
     !write(1, *) "Variables = X Y Z"
     write(1, *) "Variables = X Y Z IBLANK NEARWALL"

     ! Write my own blocks first
     do nn=1,nDom
        call setPointers(nn, 1, 1)
        !write(tmpStr, *) "Proc ", 0, " Local ID", nn
        write(tmpStr, "(a,I3.3,a,I3.3,a)"), """Proc ", 0, " Local ID ", nn , """"
        write(1,*) "ZONE I=", il, " J=",jl, "K=", kl, "T=", trim(tmpStr)
        write(1, *) "DATAPACKING=BLOCK"
        write(1, *) "VARLOCATION=([1,2,3,5]=NODAL, [4]=CELLCENTERED)"

        do iDim=1, 3
           do k=1, kl
              do j=1, jl
                 do i=1, il
                    write(1, *) x(i, j, k, idim)
                 end do
              end do
           end do
        end do

        ! Iblanks are cell center values,  imposed on primal nodes, for plotting
        ! purpose only
        do k=2, kl
           do j=2, jl
              do i=2, il
                 write(1, *) iBlank(i, j, k) 
              end do
           end do
        end do

        do k=1, kl
           do j=1, jl
              do i=1, il
                 write(1, *) flowDoms(nn, 1,1)%nearWall(i, j, k)
              end do
           end do
        end do

     end do

     ! Now loop over the remaining blocks...receiving each and writing:

     do iProc=1, nProc-1
        do nn=1, nDomProc(iProc)
           iDom = cumDomProc(iProc) + nn
           bufSize = dims(1, iDom)*dims(2, iDom)*dims(3,iDom)*3
           ibufSize = (dims(1, iDom)-1)*(dims(2, iDom)-1)*(dims(3,iDom)-1)*2

           call MPI_Recv(buffer, bufSize, sumb_real, iProc, iProc, &
                sumb_comm_world, status, ierr)

           call MPI_Recv(ibuffer, ibufSize, sumb_integer, iProc, 10*iProc, &
                sumb_comm_world, status, ierr)

           write(tmpStr, "(a,I3.3,a,I3.3,a)"), """Proc ", iProc, " Local ID ", nn ,""""
           write(1,*) "ZONE I=", dims(1, iDom), " J=", dims(2, iDom), "K=", dims(3, iDom), "T=", trim(tmpStr)
           write(1, *) "DATAPACKING=BLOCK"
           write(1, *) "VARLOCATION=([1,2,3]=NODAL, [4,5]=CELLCENTERED)"

           ! Dump directly...already in the right order
           do i=1, bufSize
              write(1, *), buffer(i)
           end do
           do i=1, ibufSize
              write(1, *), ibuffer(i)
           end do
        end do
     end do
     close(1)
     print *,'Done writing mesh...'
  else 

     ! Pand and send my stuff:
     do nn=1, nDom
        call setPointers(nn, 1, 1)
        ii = 0
        do iDim=1, 3
           do k=1, kl
              do j=1, jl
                 do i=1, il
                    ii = ii + 1
                    buffer(ii) = x(i, j, k, iDim)
                 end do
              end do
           end do
        end do

        call mpi_send(buffer, ii, sumb_real, 0, myid, &
             sumb_comm_world, ierr)

        ! Iblanks are cell center values,  imposed on primal nodes, for plotting
        ! purpose only
        ii = 0
        do k=2, kl
           do j=2, jl
              do i=2, il
                 ii = ii + 1
                 ibuffer(ii) = iBlank(i, j, k)
              end do
           end do
        end do

        do k=2, kl
           do j=2, jl
              do i=2, il
                 ii = ii + 1
                 ibuffer(ii) = flowDoms(nn, 1,1)%nearWall(i,j,k)
              end do
           end do
        end do

        call mpi_send(ibuffer, ii, sumb_integer, 0, 10*myid, &
             sumb_comm_world, ierr)
     end do
  end if

  deallocate(buffer, ibuffer, nDomProc, cumDomProc, dims)

end subroutine writePartionedMesh

subroutine writeWalls

  use communication
  use overset
  use constants
  use blockPointers
  use BCTypes
  implicit none

  character(80) :: fileName, zoneName
  integer(kind=intType) :: i, j, nn, iDom, iBeg, iEnd, jBeg, jEnd, mm, iDim
  real(kind=realType), dimension(:, :, :), pointer :: xx

  write (fileName,"(a,I2.2,a)") "wall_", myid, ".dat"

  open(unit=101,file=trim(fileName),form='formatted')
  write(101,*) 'TITLE = "mywalls"'
  write(101,*) 'Variables = "X", "Y", "Z", "CellIBlank"'

  do nn=1,nDom
     iDom = nn + cumDomProc(myid)
     call setPointers(nn, 1, 1)
     if (nBocos > 0) then 
        do mm=1, nBocos
           jBeg = BCData(mm)%jnBeg ; jEnd = BCData(mm)%jnEnd
           iBeg = BCData(mm)%inBeg ; iEnd = BCData(mm)%inEnd
           if (BCType(mm) == EulerWall .or. &
                BCType(mm) == NSWallAdiabatic .or. &
                BCType(mm) == NSWallIsothermal) then
              select case (BCFaceID(mm))
              case (iMin)
                 xx => x(1,:,:,:)
              case (iMax)
                 xx => x(il,:,:,:)
              case (jMin)
                 xx => x(:,1,:,:)
              case (jMax)
                 xx => x(:,jl,:,:)
              case (kMin)
                 xx => x(:,:,1,:)
              case (kMax)
                 xx => x(:,:,kl,:)
              end select

              write(zoneName, "(a,I2.2,a,I2.2)") "Zone", iDom, "_Proc_", myid
110           format('ZONE T=',a, " I=", i5, " J=", i5)
              write(101, 110), trim(zoneName), iEnd-iBeg+1, jEnd-jBeg+1
              write (101,*) "DATAPACKING=BLOCK, VARLOCATION=([1,2,3]=NODAL, [4]=CELLCENTERED)"

13            format (E14.6)
              do iDim=1,3
                 do j=jBeg, jEnd
                    do i=iBeg, iEnd
                       write(101, *) xx(i+1, j+1, iDim)
                    end do
                 end do
              end do

              do j=jBeg+1, jEnd
                 do i=iBeg+1, iEnd
                    write(101, *) BCData(mm)%iBlank(i, j)
                 end do
              end do
           end if
        end do
     else
        ! Write dummy zone
        write(zoneName, "(a,I1,a,I1)") "Zone", 0, "_Proc_", myid
        write(101, 110), trim(zoneName), 1, 1
        write (101,*) "DATAPACKING=POINT"
        write(101, *) zero, zero, zero, one, one
     end if
  end do
  close(101)
end subroutine writeWalls

subroutine writeOversetString(string, fileID)

  use communication
  use overset
  implicit none

  type(oversetString), intent(inout) :: string
  integer(kind=intType), intent(in) :: fileID
  integer(kind=intType) :: i, j
  character(80) :: zoneName  

 
  write (zoneName,"(a,I2.2)") "Zone T=gap_", string%myID
  write (fileID, *) trim(zoneName)
  
  write (fileID,*) "Nodes = ", string%nNodes, " Elements= ", string%nElems, " ZONETYPE=FELINESEG"
  write (fileID,*) "DATAPACKING=POINT"
13 format (E20.12)
  
  do i=1, string%nNodes
     ! Write the coordinates
     do j=1, 3
        write(fileID,13, advance='no') string%x(j, i)
     end do
     if (associated(string%otherX)) then 
        do j=1, 3
           write(fileID,13, advance='no') string%otherX(j, i) - string%x(j, i) !string%norm(j, i)
        end do
     else
        do j=1, 3
           write(fileID,13, advance='no') string%norm(j, i)
        end do
     end if
     write(fileID,13, advance='no') real(string%ind(i))
     write(fileID,13, advance='no') real(string%myID)
     write(fileID,13, advance='no') real(i)
     if (associated(string%otherID)) then 
        write(fileID,13, advance='no') real(string%otherID(1, i))
        write(fileID,13, advance='no') real(string%otherID(2, i))
     else
        write(fileID,13, advance='no') real(zero)
        write(fileID,13, advance='no') real(zero)
     end if

     write(fileID,"(1x)")
  end do
  
15 format(I5, I5)
  do i=1, string%nElems
     write(fileID, 15) string%conn(1, i), string%conn(2, i)
  end do
  
end subroutine writeOversetString



subroutine writeOversetTriangles(string, fileName)

  use communication
  use overset
  implicit none

  type(oversetString), intent(inout) :: string
  character(*) :: fileName
  integer(kind=intType) :: i, j
  character(80) :: zoneName  

  open(unit=101, file=trim(fileName), form='formatted')
  write(101,*) 'TITLE = "Triangles"'
  write(101,*) 'Variables = "X", "Y", "Z"'

  write (zoneName,"(a,I2.2)") "Zone T=triangles_", string%myID
  write (101, *) trim(zoneName)
  
  write (101,*) "Nodes = ", string%nNodes, " Elements= ", string%nTris, " ZONETYPE=FETRIANGLE"
  write (101,*) "DATAPACKING=POINT"
13 format (E20.12)

  ! Write all the coordinates  
  do i=1, string%nNodes
     do j=1, 3
        write(101,13, advance='no') string%x(j, i)
     end do
     write(101,"(1x)")
  end do

  
15 format(I5, I5)
  do i=1, string%nTris
     write(101, 15) string%tris(1, i), string%tris(2, i), string%tris(3, i)
  end do
  close(101)
end subroutine writeOversetTriangles

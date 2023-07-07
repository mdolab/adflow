module oversetAPI
    !       computeOversetInterpolation is the top level routine that
    !       implements the implicit hole cutting method for determing
    !       overset grid connectivitiies. It operates on a given multigrid
    !       level and spectral instance
contains
    subroutine oversetComm(level, firstTime, coarseLevel, closedFamList)

        use constants
        use communication, only: adflow_comm_world, sendRequests, &
                                 recvRequests, sendBuffer, recvBuffer, commPatternCell_2nd, &
                                 internalCell_2nd, sendBufferSize, recvBufferSize, myid, &
                                 nProc
        use blockPointers, only: flowDoms, nDom, fringeType, fringes, &
                                 il, jl, kl, ie, je, ke, x, nx, ny, nz, iBlank, globalCell, ib, jb, kb, nDonors, &
                                 vol, fringePtr, forcedRecv, status, nbkglobal, si, sj, sk
        use oversetData, only: CSRMatrix, oversetBlock, oversetFringe, &
                               oversetWall, nClusters, cumDomProc, localWallFringes, nDomTotal, &
                               nLocalWallFringe, clusterWalls, oversetPresent, nDomProc, &
                               overlapMatrix, tmpFringePtr, oversetTimes
        use cgnsGrid, only: cgnsDoms
        use stencils, only: N_visc_drdw, visc_drdw_stencil
        use inputTimeSpectral, only: nTimeIntervalsSpectral
        use adtBuild, only: destroySerialQuad
        use inputOverset, onlY: useoversetLoadBalance, overlapFactor, nRefine, backgroundVolScale, &
                                useOversetWallScaling
        use utils, only: EChk, setPointers, setBufferSizes, terminate, returnFail, mynorm2
        use surfaceFamilies, only: BCFamGroups
        use kdtree2_module, onlY: kdtree2_create, kdtree2destroy
        use oversetInitialization, only: initializeOBlock, initializeOFringes, initializeStatus, &
                                         reInitializeStatus
        use oversetCommUtilities, only: recvOBlock, recvOFringe, getCommPattern, getOSurfCommPattern, &
                                        emptyOversetComm, exchangeStatusTranspose, exchangeStatus, oversetLoadBalance, &
                                        exchangeFringes, sendOFringe, sendOBlock, setupFringeGlobalInd, &
                                        getFringeReturnSizes
        use oversetUtilities, only: isCompute, checkOverset, irregularCellCorrection, &
                                    fringeReduction, transposeOverlap, setIBlankArray, &
                                    deallocateOFringes, deallocateoBlocks, &
                                    deallocateOSurfs, deallocateCSRMatrix, &
                                    setIsCompute, getWorkArray, flagForcedRecv, &
                                    qsortFringeType, isReceiver, setIsReceiver, &
                                    addToFringeList, printOverlapMatrix, &
                                    tic, toc, unwindIndex, windIndex, isFloodSeed, &
                                    isFlooded, wallsOnBlock
        use oversetPackingRoutines, only: packOFringe, packOBlock, unpackOFringe, unpackOBlock, &
                                          getOFringeBufferSizes, getOBlockBufferSizes, getOSurfBufferSizes
        implicit none

        ! Input Parameters
        integer(kind=intType), intent(in) :: level
        logical :: firstTime, coarseLevel
        integer(kind=intType) :: sps
        integer(kind=intType), intent(in) :: closedFamList(:)

        ! Local Variables
        integer(kind=intType) :: i, ii, j, jj, k, kk, i_stencil, curI, curJ, curK, mInt
        integer(kind=intType) :: m, iSize, iStart, iEnd, index, rSize, iFringe
        integer(kind=intType) :: iDom, jDom, iDim, nReal, iCol
        integer(kind=intType) :: nn, mm, n, ierr, iProc, iRefine
        integer(kind=intType) :: iWork, nWork, nLocalFringe, totalOrphans
        integer(kind=intType) :: myBlock, myINdex, dIndex, donorBlock, donorProc, absDBlock
        real(kind=realType) :: startTime, endTime, curQuality, aspect(3), fact
        logical :: localChanged, globalChanged, wallsPresent

        type(CSRMatrix), pointer :: overlap
        type(CSRMatrix) :: overlapTranspose

        integer(kind=intType), dimension(:), allocatable :: cumFringeRecv, fringeRecvSizes
        integer(kind=intType), dimension(:, :), allocatable :: work, tmpInt2D
        integer(kind=intType), dimension(:), pointer :: wallFamList
        real(kind=realType), dimension(:), allocatable :: tmpReal
        real(kind=realType), dimension(:, :), allocatable :: xMin, xMax

        logical, dimension(:), allocatable :: oBlockReady, oFringeReady

        type(oversetBlock), dimension(:), allocatable :: oBlocks
        type(oversetFringe), dimension(:), allocatable :: oFringes
        type(oversetWall), pointer :: wall
        type(fringeType), dimension(:), allocatable :: localFringes
        type(fringeType) :: fringe

        ! MPI/Communication related
        integer mpiStatus(MPI_STATUS_SIZE)
        integer(kind=intType) :: MAGIC, source, tag, sendCount, recvCount
        integer(kind=intType) :: nOFringeSend, nOFringeRecv
        integer(kind=intType) :: nOBlockSend, nOBlockRecv
        integer(kind=intType) :: nOSurfSend, nOSurfRecv
        logical :: flag, invalid

        integer(kind=intType), dimension(:, :), allocatable :: oBlockSendList, oBlockRecvList
        integer(kind=intType), dimension(:, :), allocatable :: oFringeSendList, oFringeRecvList
        integer(kind=intType), dimension(:, :), allocatable :: bufSizes, recvInfo
        integer(kind=intType), dimension(:), allocatable :: intRecvBuf
        real(kind=realType), dimension(:), allocatable :: realRecvBuf

        ! -----------------------------------------------------------------
        ! Step 1: Initializaion: Make sure the stencils are initialized.
        ! -----------------------------------------------------------------

        call initialize_stencils()

        ! -----------------------------------------------------------------
        ! Step 2: Communicate the block size info to everyone. Also generate
        ! cumDomProc which is the cumulative form. This will make our lives
        ! easier a little later on with indexing. (Routine below)
        ! -----------------------------------------------------------------

        ! If there is not overset meshes present, just make an empty comm
        ! structure and call it a day.
        if (.not. oversetPresent) then
            do sps = 1, nTimeIntervalsSpectral
                call emptyOversetComm(level, sps)

                do nn = 1, nDom
                    flowDoms(nn, level, sps)%nOrphans = 0
                    if (.not. associated(flowDoms(nn, level, sps)%iblank)) then
                        i = flowDoms(nn, level, sps)%ib
                        j = flowDoms(nn, level, sps)%jb
                        k = flowDoms(nn, level, sps)%kb

                        allocate (flowDoms(nn, level, sps)%iblank(0:i, 0:j, 0:k))
                        flowDoms(nn, level, sps)%iblank = 1
                    end if
                    do mm = 1, flowDoms(nn, level, sps)%nBocos
                        flowDoms(nn, level, sps)%BCData(mm)%iblank = 1
                    end do
                end do
            end do
            return
        end if

        ! Determine the magic number which is actually the same as
        ! nDomTotal. nDomTotal is guaranteed to be greater than or equal to
        ! nProc. This will be used to space out tags on communications to
        ! make sure they do not overlap.
        MAGIC = nDomTotal

        ! Zero out all the timers
        oversetTimes = zero

        call MPI_barrier(adflow_comm_world, ierr)
        call tic(iTotal)
        ! Master SPS loop
        spectralLoop: do sps = 1, nTimeIntervalsSpectral

            ! Set a pointer to make the code easier to read
            overlap => overlapMatrix(level, sps)

            ! -----------------------------------------------------------------
            ! Step 4: Compute the 3D axis oriented bounding boxes for each block
            ! and communicate then to everyone. Also determine communicate the
            ! minimum volume for each block to everyone.  (Routine below)
            ! -----------------------------------------------------------------

            call tic(iBoundingBox)
            allocate (xMin(3, nDomTotal), xMax(3, nDomTotal))
            call computeDomainBoundingBoxes
            call toc(iBoundingBox)

            ! -----------------------------------------------------------------
            ! Step 8: Build a global sparse matrix representation of the overlap
            ! matrix.  Every processor will have the same sparse matrix
            ! representation when it is finished. (Routine below)
            ! -----------------------------------------------------------------
            call tic(iBuildOverlap)
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

            allocate (bufSizes(nDomTotal, 6), tmpInt2D(nDomTotal, 6), &
                      tmpReal(size(overlap%data)))

            wallFamList => BCFamGroups(iBCGroupWalls)%famList
            ! Initialization
            tmpReal = zero
            tmpInt2D = 0
            do nn = 1, nDom
                call setPointers(nn, level, sps)
                iDom = cumDomProc(myid) + nn
                ! ------------------------------
                ! Old Code for setting the data
                ! if (firstTime) then
                !    do jj=overlap%rowPtr(iDom), overlap%rowPtr(iDom+1)-1
                !       tmpReal(jj) = real(nx*ny*nz)
                !    end do
                ! end if
                ! ------------------------------

                ! Sizes
                call getOBlockBufferSizes(il, jl, kl, tmpInt2D(iDom, 1), tmpInt2D(iDom, 2))
                call getOFringeBufferSizes(il, jl, kl, tmpInt2D(iDom, 3), tmpInt2D(iDom, 4))
                call getOSurfBufferSizes(wallFamList, il, jl, kl, tmpInt2D(iDom, 5), tmpInt2D(iDom, 6), .True.)

            end do

            ! Set the tmpReal variable to to be the number of fringes we
            ! need to search for. This is inefficient becuase we loop over
            ! the all the matrix rows:
            if (firstTime) then
                do i = 1, overlap%nrow
                    do jj = overlap%rowPtr(i), overlap%rowPtr(i + 1) - 1
                        iCol = overlap%colInd(jj)
                        if (iCol > cumDomProc(myid) .and. iCol <= cumDomProc(myid + 1)) then
                            nn = iCol - cumDomProc(myid)
                            tmpReal(jj) = &
                                flowDoms(nn, level, sps)%nx * &
                                flowDoms(nn, level, sps)%ny * &
                                flowDoms(nn, level, sps)%nz
                        end if
                    end do
                end do
            end if

            if (.not. firstTime) then
                tmpReal = overlap%data
            end if

            ! Determine the total search costs for each proc and all the bufferSizes
            call mpi_allreduce(tmpReal, overlap%data, overlap%nnz, adflow_real, MPI_SUM, &
                               adflow_comm_world, ierr)
            call ECHK(ierr, __FILE__, __LINE__)

            call mpi_allreduce(tmpInt2D, bufSizes, 6 * nDomTotal, adflow_integer, MPI_SUM, &
                               adflow_comm_world, ierr)
            call ECHK(ierr, __FILE__, __LINE__)

            ! Done with the tmp arrays. This should be the last of the
            ! blocking collectives for a while.
            deallocate (tmpReal, tmpInt2D)

            ! -----------------------------------------------------------------
            ! Step 8: We are now ready to partiaion and loadbalance the work
            ! based on the costs stored in the overlap matrix. These costs
            ! may be the search estimates from initializeOverlapCosts OR they
            ! may be actual timings from a previous assembly. Also create a
            ! transpose of the matrix which is useful to use for the fringe
            ! sending operation (it is transpose of the block sending)
            ! -----------------------------------------------------------------
            useOversetLoadBalance = .False.
            if (useOversetLoadBalance) then
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

            call getWorkArray(overlap, work)
            nWork = size(work, 2)

            ! Call the generic routines to determine the send/receive pattern
            ! for oBlock comm and the fringe comm. These are transpose of
            ! each other. Just overestimate the sizes of the lists.

            ! For sending, the worse case is sending all my blocks/fringes/walls to
            ! everyone but myself:
            ii = nDom * (nProc - 1)
            allocate (oBlockSendList(2, ii), oFringeSendList(2, ii))

            ! For receiving, the worse receive is all the blocks/fringes/wall I
            ! don't already have:
            ii = nDomTotal - nDom
            allocate (oBlockRecvList(2, ii), oFringeRecvList(2, ii))

            call getCommPattern(overlap, oblockSendList, nOblockSend, &
                                oBlockRecvList, nOblockRecv)

            call getCommPattern(overlapTranspose, oFringeSendList, nOFringeSend, &
                                oFringeRecvList, nOFringeRecv)

            ! Done with the transposed matrix
            call deallocateCSRMatrix(overlapTranspose)

            ! Zero out the overlap data since we will be doing new timings in
            ! doMyWork()
            overlap%data = zero

            ! Allocate the exact space for our send and recv requests. Note
            ! that for the oBlocks, two values are set, real and integer.
            nn = max(2 * nProc, &
                     2 * nOBlockSend + 2 * nOFringeSend, &
                     2 * nOBlockRecv + 2 * nOfringeRecv)
            if (allocated(sendRequests)) then
                deallocate (sendRequests, recvRequests)
            end if

            allocate (sendRequests(nn), recvRequests(nn))
            allocate (recvInfo(2, nn))

            ! On the first pass we need to get an estimate of what is
            ! inside the body and what isn't. This method isn't
            ! perfect; some cells that are actually inside the true
            ! surface won't be flagged, but that's ok.
            work(4, :) = 0

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

            allocate (oBlocks(nDomTotal), oFringes(nDomTotal))

            ! Thse variables keep track of if the block/fringes are
            ! ready. Initialized to false and only flipped when we are sure
            ! they are ready to be used.

            allocate (oBlockReady(nDomTotal), oFringeReady(nDomTotal))
            oBlockReady = .False.
            oFringeReady = .False.

            ! Allocate space for the localWallFringes. localWallFringes keeps
            ! track of donors for cells that are next to a wall. These must
            ! be recorded independently of the actual donors since we don't
            ! actually care what the interpolation stencil is, rather just
            ! who the donor is such that we can use that information for the
            ! flooding process. We arbitrarily set a size here and it will be
            ! automatically expanded as necessary in the fringeSearch
            ! routine.
            allocate (localWallFringes(1000))
            nLocalWallFringe = 0

            call toc(iBuildOverlap)

            call tic(iBuildClusterWalls)
            allocate (clusterWalls(nClusters))
            call buildClusterWalls(level, sps, .True., clusterWalls, wallFamList, size(wallFamList))
            call toc(iBuildClusterWalls)

            ! Determine the cells that are near wall. We have a special
            ! routine for this.
            call tic(iComputeCellWallPoint)
            call computeCellWallPoint(level, sps, clusterWalls)

            ! We need a couple of extra things that buildCluster wall
            ! doesn't do:
            do ii = 1, nClusters
                wall => clusterWalls(ii)
                if (wall%nNodes > 0) then
                    wall%tree => kdtree2_create(wall%x(:, 1:wall%nNodes))
                end if

                ! Build the inverse of the connectivity, the nodeToElem array.
                allocate (wall%nte(4, wall%nNodes))
                wall%nte = 0
                do i = 1, wall%nCells
                    do j = 1, 4
                        n = wall%conn(j, i)
                        inner: do k = 1, 4
                            if (wall%nte(k, n) == 0) then
                                wall%nte(k, n) = i
                                exit inner
                            end if
                        end do inner
                    end do
                end do
            end do
            call toc(iComputeCellWallPoint)

            ! Flag all the cells that we know are forced receivers.
            call flagForcedRecv()

            ! Initialize the overset-specific data structures for
            ! performing searches.
            do nn = 1, nDom
                call setPointers(nn, level, sps)
                iDom = cumDomProc(myid) + nn

                call tic(iBuildADT)
                call initializeOBlock(oBlocks(iDom), nn, level, sps)
                oBlockReady(iDom) = .True.
                call toc(iBuildADT)

                call tic(iBuildSearchPoints)
                call initializeOFringes(oFringes(iDom), nn, closedFamList)
                oFringeReady(iDom) = .True.
                call toc(iBuildSearchPoints)
            end do

            ! Post all the oBlock/oFringe iSends
            sendCount = 0
            do jj = 1, nOblockSend
                iProc = oBlockSendList(1, jj)
                iDom = oBlockSendList(2, jj)
                call packOBlock(oBlocks(iDom))
                call sendOBlock(oBlocks(iDom), iDom, iProc, 0, sendCount)
            end do

            do jj = 1, nOFringeSend
                iProc = oFringeSendList(1, jj)
                iDom = oFringeSendList(2, jj)
                call packOFringe(oFringes(iDom))
                call sendOFringe(oFringes(iDom), iDom, iProc, MAGIC, sendCount)
            end do

            ! Post all the oBlock/oFringe receives. Before posting the actual
            ! receive, allocate the receiving buffer.
            recvCount = 0
            do jj = 1, nOBlockRecv
                iProc = oBlockRecvList(1, jj)
                iDom = oBlockRecvList(2, jj)
                call recvOBlock(oBlocks(iDom), iDom, iProc, 0, &
                                bufSizes(iDom, 1), bufSizes(iDom, 2), recvCount, recvInfo)
            end do

            do jj = 1, nOFringeRecv
                iProc = oFringeRecvList(1, jj)
                iDom = oFringeRecvList(2, jj)
                call recvOFringe(oFringes(iDom), iDom, iProc, MAGIC, &
                                 bufSizes(iDom, 3), bufSizes(iDom, 4), recvCount, recvInfo)
            end do

            ! Before we start waiting for the receives to finish, we can see
            ! if we can do any searches with the blocks/fringes we already
            ! have. Call the internal routine for this.
            call doMyWork(flag)

            ! Complete all the recives
            do i = 1, recvCount

                ! Complete any one of the recv requests
                call mpi_waitany(recvCount, recvRequests, index, mpiStatus, ierr)
                call ECHK(ierr, __FILE__, __LINE__)

                ! Global domain index of the recv that finished
                iDom = recvInfo(1, index)

                ! Check which type of receive just finished and flag them as
                ! being complete.
                if (recvInfo(2, index) == 1) then
                    oBlocks(iDom)%realBufferReady = .True.
                else if (recvInfo(2, index) == 2) then
                    oBlocks(iDom)%intBufferReady = .True.
                else if (recvInfo(2, index) == 3) then
                    oFringes(iDom)%realBufferReady = .True.
                else if (recvInfo(2, index) == 4) then
                    oFringes(iDOm)%intBufferReady = .True.
                end if

                ! If both int and real buffers are received, we can unpack the
                ! oblock and flag it as ready.
                if (oBlocks(iDom)%realBufferReady .and. oBlocks(iDom)%intBufferReady .and. &
                    .not. oBlocks(iDom)%allocated) then
                    call unpackOBlock(oBlocks(iDom))
                    oBlockReady(iDom) = .True.
                end if

                ! If both int and real buffers are received, we can unpack the
                ! oFringe and flag it as ready.
                if (oFringes(iDom)%realBufferReady .and. oFringes(iDom)%intBufferReady .and. &
                    .not. oFringes(iDom)%allocated) then
                    call unpackOFringe(oFringes(iDom))
                    oFringeReady(iDom) = .True.
                end if

                ! Now see if we can do any more of the work, ie the searches.
                call doMyWork(flag)

                ! Sanity check. flag better be true when i=recvCount
                if (i == recvCount .and. .not. flag) then
                    call terminate("computeInterpolationParallel", "Inconsistent Comm pattern detected.")
                end if
            end do

            ! Last thing to do wait for all the sends to finish
            do i = 1, sendCount
                call mpi_waitany(sendCount, sendRequests, index, mpiStatus, ierr)
                call ECHK(ierr, __FILE__, __LINE__)
            end do

            ! We are now completely finished with oBlocks so
            ! delete before we allocate space for all the fringes
            call deallocateOBlocks(oBlocks, size(oBlocks))
            deallocate (oBlocks)
            deallocate (oBlockReady, oFringeReady)

            ! Destroy the walls we had made for the surface-overlap
            ! searching.
            do i = 1, nClusters
                wall => clusterWalls(i)
                call destroySerialQuad(wall%ADT)
                if (wall%nNodes > 0) then
                    call kdtree2destroy(wall%tree)
                end if
                deallocate (wall%x, wall%conn, wall%ind, wall%nte)
            end do
            deallocate (clusterWalls)

            ! -----------------------------------------------------------------
            ! Step 9: Well, all the searches are done, so now we can now send
            ! the fringes back to where they came from.
            ! -----------------------------------------------------------------

            ! We are done with most of the data of the oFringe. Nuke that
            ! right now to save memory. The fringeIntBuffer and
            ! fringeRealBuffer have to stick around a little longer.
            do iDim = 1, nDomTotal
                if (allocated(oFringes(iDom)%x)) &
                    deallocate (oFringes(iDom)%x)

                if (allocated(oFringes(iDom)%isWall)) &
                    deallocate (oFringes(iDom)%isWall)

                if (allocated(oFringes(iDom)%xSeed)) &
                    deallocate (oFringes(iDom)%xSeed)

                if (allocated(oFringes(iDom)%wallind)) &
                    deallocate (oFringes(iDom)%wallInd)
            end do

            call tic(iFringeProcessing)
            nReal = 4
            mInt = 5
            do iDom = 1, nDomTotal
                if (oFringes(iDom)%allocated) then
                    ! Fringe is allocated so check it
                    oFringes(iDom)%fringeReturnSize = oFringes(iDom)%nDonor

                    ! Check if this domain is one I own. If so, we just copy
                    ! the the fringes into the block-based list:

                    nn = iDom - cumDomProc(myid)
                    if (nn >= 1 .and. nn <= nDom) then
                        call setPointers(nn, level, sps)
                        if (associated(flowDoms(nn, level, sps)%fringes)) then
                            deallocate (flowDoms(nn, level, sps)%fringes)
                            deallocate (flowDoms(nn, level, sps)%fringePtr)
                            deallocate (flowDoms(nn, level, sps)%nDonors)
                        end if
                        if (associated(flowDoms(nn, level, sps)%gInd)) then
                            deallocate (flowDoms(nn, level, sps)%gInd)
                        end if

                        ! Estimate about 1 donors for every cell.
                        mm = nx * ny * nz
                        allocate (flowDoms(nn, level, sps)%fringes(mm), &
                                  flowDoms(nn, level, sps)%fringePtr(3, 0:ib, 0:jb, 0:kb), &
                                  flowDoms(nn, level, sps)%nDonors)
                        flowDoms(nn, level, sps)%fringePtr = 0
                        flowDoms(nn, level, sps)%nDonors = 0

                        ! Reset the pointers due to the allocation
                        call setPointers(nn, level, sps)

                        ! Copy over the data from oFringes and then nuke the oFringe
                        do ii = 1, oFringes(iDom)%nDonor
                            ! Reconsitute the fringe itself.
                            fringe%donorProc = oFringes(iDom)%fringeIntBuffer(1, ii)
                            fringe%donorBlock = oFringes(iDom)%fringeIntBuffer(2, ii)
                            fringe%dIndex = oFringes(iDom)%fringeIntBuffer(3, ii)
                            fringe%myBlock = oFringes(iDom)%fringeIntBuffer(4, ii)
                            fringe%myIndex = oFringes(iDom)%fringeIntBuffer(5, ii)
                            fringe%donorFrac = oFringes(iDom)%fringeRealBuffer(1:3, ii)
                            fringe%quality = oFringes(iDom)%fringeRealBuffer(4, ii)

                            call addToFringeList(flowDoms(nn, level, sps)%fringes, nDonors, &
                                                 fringe)
                        end do
                    end if
                end if
            end do
            call toc(iFringeProcessing)

            ! This will compute the amount to data we expect to receive
            ! back for each OFringe we set out.
            call getFringeReturnSizes(oFringeSendList, oFringeRecvList, &
                                      nOFringeSend, nOFringeRecv, oFringes, fringeRecvSizes, cumFringeRecv)

            ! Now alocate the integer and real space. Note we are receiving nReal real
            ! values and mInt int values:
            ii = cumFringeRecv(nOfringeSend + 1) - 1
            allocate (intRecvBuf(ii * mInt), realRecvBuf(ii * nReal))

            ! We are now ready to actually send and receive our fringes
            sendCount = 0
            do jj = 1, nOFringeRecv

                iProc = oFringeRecvList(1, jj)
                iDom = oFringeRecvList(2, jj)
                iSize = oFringes(iDom)%fringeReturnSize
                if (iSize > 0) then
                    tag = iDom + MAGIC
                    sendCount = sendCount + 1
                    call mpi_isend(oFringes(iDom)%fringeRealBuffer, iSize * nReal, adflow_real, &
                                   iproc, tag, adflow_comm_world, sendRequests(sendCount), ierr)
                    call ECHK(ierr, __FILE__, __LINE__)

                    tag = iDom + 2 * MAGIC
                    sendCount = sendCount + 1
                    call mpi_isend(oFringes(iDom)%fringeIntBuffer, iSize * mInt, adflow_integer, &
                                   iproc, tag, adflow_comm_world, sendRequests(sendCount), ierr)
                    call ECHK(ierr, __FILE__, __LINE__)
                end if
            end do

            ! Non-blocking receives
            recvCount = 0
            do jj = 1, nOfringeSend

                iProc = oFringeSendList(1, jj)
                iDom = oFringeSendList(2, jj)
                iSize = cumFringeRecv(jj + 1) - cumFringeRecv(jj)
                if (iSize > 0) then

                    iStart = (cumFringeRecv(jj) - 1) * nReal + 1
                    tag = iDom + MAGIC
                    recvCount = recvCount + 1
                    call mpi_irecv(realRecvBuf(iStart), iSize * nReal, adflow_real, &
                                   iProc, tag, adflow_comm_world, recvRequests(recvCount), ierr)
                    call ECHK(ierr, __FILE__, __LINE__)
                    recvInfo(:, recvCount) = (/iDom, 1/) ! 1 for real recv

                    iStart = (cumFringeRecv(jj) - 1) * mInt + 1
                    tag = iDom + 2 * MAGIC
                    recvCount = recvCount + 1
                    call mpi_irecv(intRecvBuf(iStart), iSize * mInt, adflow_integer, &
                                   iProc, tag, adflow_comm_world, recvRequests(recvCount), ierr)
                    call ECHK(ierr, __FILE__, __LINE__)
                    recvInfo(:, recvCount) = (/iDom, 2/) ! 2 for int recv
                end if
            end do

            ! Now wait for the sends and receives to finish
            do i = 1, sendCount
                call mpi_waitany(sendCount, sendRequests, index, mpiStatus, ierr)
                call ECHK(ierr, __FILE__, __LINE__)
            end do

            ! With the send finished, we can ditch the oFringe real buffer
            ! that contains the frac and quality. That is no longer needed.
            do iDom = 1, nDom
                if (associated(oFringes(iDom)%fringeRealBuffer)) then
                    deallocate (oFringes(iDom)%fringeRealBuffer)
                end if
            end do

            do i = 1, recvCount
                call mpi_waitany(recvCount, recvRequests, index, mpiStatus, ierr)
                call ECHK(ierr, __FILE__, __LINE__)
            end do

            call tic(iFringeProcessing)

            ! Process the data we just received.
            do kk = 1, nOfringeSend

                ! Local block index of the fringes
                iDom = oFringeSendList(2, kk)
                nn = iDom - cumDomProc(myid)

                ! Set the block pointers for the local block we are dealing
                ! with:
                call setPointers(nn, level, sps)

                ! This is the range of fringes that are now ready.
                do jj = cumFringeRecv(kk), cumFringeRecv(kk + 1) - 1

                    ! Recreate the fringe type
                    iStart = mInt * (jj - 1)
                    fringe%donorProc = intRecvBuf(iStart + 1)
                    fringe%donorBlock = intRecvBuf(iStart + 2)
                    fringe%dIndex = intRecvBuf(iStart + 3)
                    fringe%myBlock = intRecvBuf(iStart + 4)
                    fringe%myIndex = intRecvBuf(iStart + 5)
                    iStart = nReal * (jj - 1)
                    fringe%donorFrac = realRecvBuf(iStart + 1:iStart + 3)
                    fringe%quality = realRecvBuf(iStart + 4)

                    ! Add directly to the list.
                    call addToFringeList(flowDoms(nn, level, sps)%fringes, nDonors, fringe)
                end do
            end do

            ! Ditch the temporary receiving memory for reals. We still the
            ! integer receiving buffer for inside the loop.
            deallocate (realRecvBuf)

            ! Before we start the loop, for each block, we record the total
            ! number of donors, that is the current value of nDonors. This
            ! is necessary because on sebsequent iterations, the
            ! exchangeFringes will have added the extra ones that are
            ! necessary for the local halo cells. So we save the current
            ! value of nDonors so that when we start the loop we reset the
            ! value back this. We also sort the fringes as this only needs
            ! to be completed once.

            do nn = 1, nDom
                call setPointers(nn, level, sps)

                ! First we need to sort the fringes by RECEIVER. That will
                ! put all the cells in order.
                call qsortFringeType(fringes, nDonors, sortByReceiver)

                flowDoms(nn, level, sps)%nDonorsOnOwnedCells = flowDoms(nn, level, sps)%nDonors

                ! In this loop we assign the start and end indices for the
                ! fringePtr. This also is invarient of the subsequent
                ! interations. Note that the '1' index is not set, that is
                ! we're not saying which of hte potential frings it could
                ! use, only the range of the potential fringes.

                curI = 0
                curJ = 0
                curK = 0

                do ii = 1, nDonors

                    ! fringePtr has the following 3 integers:

                    ! 1: The actual fringe this cell is using. 0 if it
                    ! doesn't have a fringe

                    ! 2: The start index into fringes of all of this cell's
                    ! fringes. 0 if cell has no fringes.

                    ! 3: The end index into fringes of all of this cell's
                    ! fringes. 0 if cell has no fringes.

                    call unwindIndex(fringes(ii)%myIndex, il, jl, kl, i, j, k)
                    ! Set the start and end to the current index
                    if (curI /= i .or. curJ /= j .or. curK /= k) then

                        fringePtr(2, i, j, k) = ii
                        fringePtr(3, i, j, k) = ii
                        curI = i
                        curJ = j
                        curK = k
                    else
                        ! We have the same i,j,k. Set the end index to ii
                        fringePtr(3, i, j, k) = ii
                    end if
                end do
            end do

            ! Allocate space for iblankLast which is used to keep track of
            ! the previous iblank values to determine when we can sto pcd the
            ! loop. It is initialized to 1 such that at least 2 iteations
            ! will always be done.
            do nn = 1, nDom
                call setPointers(nn, level, sps)
                ! Note that we only allocate the owned cells since it
                ! sufficient to check just the owned cells for if any one of
                ! them has changed.
                allocate (flowDoms(nn, level, sps)%iBlankLast(2:il, 2:jl, 2:kl))
                flowDoms(nn, level, sps)%iBlankLast = 1
            end do

            call toc(iFringeProcessing)

            ! Wall donors are now set. So we can do them once before the
            ! loop starts and then ditch the memory. Status must be
            ! initialized first.
            call initializeStatus(level, sps)
            call exchangeStatus(level, sps, commPatternCell_2nd, internalCell_2nd)
            call determineDonors(level, sps, localWallFringes, nLocalWallFringe, .True.)

            call toc(iDetermineDonors)
            call exchangeStatusTranspose(level, sps, commPatternCell_2nd, internalCell_2nd)
            call exchangeStatus(level, sps, commPatternCell_2nd, internalCell_2nd)
            deallocate (localWallFringes)

            ! ---------------------------
            ! Start the refinement loop:
            ! ---------------------------
            refineLoop: do iRefine = 1, nRefine

                ! re-Initialize the status array.
                call reInitializeStatus(level, sps)

                ! Exchange the status so that the halos get the right
                call exchangeStatus(level, sps, commPatternCell_2nd, internalCell_2nd)

                call tic(iCheckDonors)
                nLocalFringe = 0
                do nn = 1, nDom
                    call setPointers(nn, level, sps)

                    ! Reset the nDonors to the total number that were found for
                    ! *just* the owned cells.

                    nDonors = flowDoms(nn, level, sps)%nDonorsOnOwnedCells

                    ! Now we *actually* perform the implicit hole cut. We
                    ! compute our own quality, and then loop over the
                    ! perspective donors and see if we can get one
                    ! better. The purpose of this code is to set the "1"
                    ! index in fringePtr.
                    call wallsOnBlock(wallsPresent)
                    do k = 2, kl
                        do j = 2, jl
                            do i = 2, il
                                ! Reset the pointer to the donor cell.
                                fringePtr(1, i, j, k) = 0

                                if (iblank(i, j, k) <= -2) then
                                    ! This is a flooded cell, don't try to get a donor
                                    cycle
                                end if

                                ! This is my original quality, accounting for
                                ! the overlap factor.
                                if (wallsPresent) then
                                    aspect = one
                                    if (useOversetWallScaling) then
                                        if (CGNSDoms(nbkGlobal)%viscousDir(1)) &
                                            aspect(1) = (half * (mynorm2(si(i - 1, j, k, :)) + &
                                                                 mynorm2(si(i, j, k, :)))) / vol(i, j, k)
                                        if (CGNSDoms(nbkGlobal)%viscousDir(2)) &
                                            aspect(2) = (half * (mynorm2(sj(i, j - 1, k, :)) + &
                                                                 mynorm2(sj(i, j, k, :)))) / vol(i, j, k)
                                        if (CGNSDoms(nbkGlobal)%viscousDir(3)) &
                                            aspect(3) = (half * (mynorm2(sk(i, j, k - 1, :)) + &
                                                                 mynorm2(sk(i, j, k, :)))) / vol(i, j, k)
                                    end if
                                    fact = min(aspect(1) * aspect(2) * aspect(3), 100.0_realType)

                                    curQuality = (vol(i, j, k)**third / fact) * overlapFactor
                                else
                                    curQuality = (backGroundVolScale * vol(i, j, k))**third * overlapFactor
                                end if
                                ! Account for explict priority
                                curQuality = curQuality * cgnsDoms(nbkGlobal)%priority

                                ! For forced receivers, they get a large quality
                                ! such that ANY donor is selected.
                                if (forcedRecv(i, j, k) > 0) then
                                    curQuality = large
                                end if

                                ! Loop over the potential donors
                                iStart = fringePtr(2, i, j, k)
                                iEnd = fringePtr(3, i, j, k)

                                ! iStart of 0 means no donor:
                                if (iStart > 0) then
                                    do ii = iStart, iEnd
                                        if (fringes(ii)%quality < curQuality) then

                                            ! Update the pointer to the fringe we want
                                            fringePtr(1, i, j, k) = ii

                                            ! Update the current quality to the one
                                            ! we just found, in case we have
                                            ! multiple ones to choose from.
                                            curQuality = fringes(ii)%quality

                                            ! Flag this cell as a receiver and remove
                                            ! it's compute status.
                                            call setIsReceiver(status(i, j, k), .True.)
                                            call setIsCompute(status(i, j, k), .False.)
                                        end if
                                    end do
                                end if
                            end do
                        end do
                    end do

                    ! Now count up the number of actual fringes we had.
                    do k = 2, kl
                        do j = 2, jl
                            do i = 2, il
                                if (isReceiver(status(i, j, k))) then
                                    nLocalFringe = nLocalFringe + 1
                                end if
                            end do
                        end do
                    end do
                end do
                call toc(iCheckDonors)

                ! -------------------------------
                !            Flooding
                ! -------------------------------
                call tic(iFlooding)
                call floodInteriorCells(level, sps)
                call toc(iFlooding)

                ! ---------------------------------
                ! Determine which cells are donors
                ! ---------------------------------
                allocate (localFringes(nLocalFringe))
                ! Load up the fringes
                nLocalFringe = 0
                do nn = 1, nDom
                    call setPointers(nn, level, sps)
                    do k = 2, kl
                        do j = 2, jl
                            do i = 2, il
                                if (isReceiver(status(i, j, k))) then
                                    nLocalFringe = nLocalFringe + 1
                                    ii = fringePtr(1, i, j, k)
                                    localFringes(nLocalFringe) = fringes(ii)
                                end if
                            end do
                        end do
                    end do
                end do
                call tic(iDetermineDonors)
                call determineDonors(level, sps, localFringes, nLocalFringe, .False.)
                call toc(iDetermineDonors)
                call exchangeStatusTranspose(level, sps, commPatternCell_2nd, internalCell_2nd)
                call exchangeStatus(level, sps, commPatternCell_2nd, internalCell_2nd)
                deallocate (localFringes)

                ! Determine which cells need to be forced receivers. This
                ! may have changed due to flooding.
                call flagForcedRecv()

                ! Flag cells that are no longer valid donors. This is
                ! slightly tricky: What we do here is

                ! 1. Firstly we go through the oFringes array and mark all
                ! the fringes that are now invalid.

                do iDom = 1, nDomTotal
                    if (oFringes(iDom)%allocated) then
                        do iFringe = 1, oFringes(iDom)%nDonor
                            donorBlock = oFringes(iDom)%fringeIntBuffer(2, iFringe)
                            dIndex = oFringes(iDom)%fringeIntBuffer(3, iFringe)
                            nn = abs(donorBlock)

                            if (nn < 1) then
                                print *, 'donor:', myid, donorBlock, cumDomProc(myID), nn
                                stop
                            end if
                            il = flowDoms(nn, level, sps)%il
                            jl = flowDoms(nn, level, sps)%jl
                            kl = flowDoms(nn, level, sps)%kl

                            call unwindIndex(dIndex, il, jl, kl, i, j, k)

                            ! Check if this is invalid:
                            invalid = .False.
                            do kk = k, k + 1
                                do jj = j, j + 1
                                    do ii = i, i + 1

                                        if (isFlooded(flowDoms(nn, level, sps)%status(ii, jj, kk)) .or. &
                                            isFloodSeed(flowDoms(nn, level, sps)%status(ii, jj, kk)) .or. &
                                            flowDoms(nn, level, sps)%iblank(ii, jj, kk) == -4 .or. &
                                            flowDoms(nn, level, sps)%forcedRecv(ii, jj, kk) > 0) then
                                            invalid = .True.
                                        end if
                                    end do
                                end do
                            end do

                            if (invalid) then
                                ! Change the donorBlock to be negative if it is
                                ! not already so.  The negative block index
                                ! indicates that this donor is no longer valid.
                                if (donorBlock > 0) then
                                    oFringes(iDom)%fringeIntBuffer(2, iFringe) = -donorBlock
                                end if
                            end if
                        end do
                    end if
                end do

                ! 2. We send the fringes back to where they came from like
                ! during the search. The only difference is we don't add
                ! them to the list again, but rather modify what is already there.

                sendCount = 0
                do jj = 1, nOFringeRecv

                    iProc = oFringeRecvList(1, jj)
                    iDom = oFringeRecvList(2, jj)
                    iSize = oFringes(iDom)%fringeReturnSize
                    if (iSize > 0) then
                        tag = iDom + 2 * MAGIC
                        sendCount = sendCount + 1
                        call mpi_isend(oFringes(iDom)%fringeIntBuffer, iSize * mInt, adflow_integer, &
                                       iproc, tag, adflow_comm_world, sendRequests(sendCount), ierr)
                        call ECHK(ierr, __FILE__, __LINE__)
                    end if
                end do

                ! Non-blocking receives
                recvCount = 0
                do jj = 1, nOfringeSend

                    iProc = oFringeSendList(1, jj)
                    iDom = oFringeSendList(2, jj)
                    iSize = cumFringeRecv(jj + 1) - cumFringeRecv(jj)
                    if (iSize > 0) then

                        iStart = (cumFringeRecv(jj) - 1) * mInt + 1
                        tag = iDom + 2 * MAGIC
                        recvCount = recvCount + 1
                        call mpi_irecv(intRecvBuf(iStart), iSize * mInt, adflow_integer, &
                                       iProc, tag, adflow_comm_world, recvRequests(recvCount), ierr)
                        call ECHK(ierr, __FILE__, __LINE__)
                        recvInfo(:, recvCount) = (/iDom, 2/) ! 2 for int recv
                    end if
                end do

                ! Process any data we can locally while we wait.
                do iDom = 1, nDomTotal
                    if (oFringes(iDom)%allocated) then

                        nn = iDom - cumDomProc(myid)
                        if (nn >= 1 .and. nn <= nDom) then
                            call setPointers(nn, level, sps)

                            do jj = 1, oFringes(iDom)%nDonor

                                donorProc = oFringes(iDom)%fringeIntBuffer(1, jj)
                                donorBlock = oFringes(iDom)%fringeIntBuffer(2, jj)
                                dIndex = oFringes(iDom)%fringeIntBuffer(3, jj)
                                myBlock = oFringes(iDom)%fringeIntBuffer(4, jj)
                                myIndex = oFringes(iDom)%fringeIntBuffer(5, jj)

                                if (donorBlock < 0) then
                                    ! This fringe is no longer valid! Modify this
                                    ! fringe in our list.
                                    call unwindIndex(myIndex, il, jl, kl, i, j, k)

                                    ! Now loop over the donors this cells has to match
                                    ! up the one that we are dealing with here.

                                    ! Loop over the potential donors.
                                    iStart = fringePtr(2, i, j, k)
                                    iEnd = fringePtr(3, i, j, k)

                                    do ii = iStart, iEnd
                                        if (fringes(ii)%donorBlock == abs(donorBlock) .and. &
                                            fringes(ii)%dIndex == dIndex) then

                                            ! This is the fringe. Flag as bad.
                                            fringes(ii)%quality = 2 * Large
                                        end if
                                    end do
                                end if
                            end do
                        end if
                    end if
                end do

                ! Now wait for the sends and receives to finish
                do i = 1, sendCount
                    call mpi_waitany(sendCount, sendRequests, index, mpiStatus, ierr)
                    call ECHK(ierr, __FILE__, __LINE__)
                end do

                do i = 1, recvCount
                    call mpi_waitany(recvCount, recvRequests, index, mpiStatus, ierr)
                    call ECHK(ierr, __FILE__, __LINE__)
                end do

                ! Process the data we just received.
                do kk = 1, nOfringeSend

                    ! Local block index of the fringes
                    iDom = oFringeSendList(2, kk)
                    nn = iDom - cumDomProc(myid)

                    ! Set the block pointers for the local block we are dealing
                    ! with:
                    call setPointers(nn, level, sps)

                    ! This is the range of fringes that are now ready.
                    do jj = cumFringeRecv(kk), cumFringeRecv(kk + 1) - 1

                        iStart = mInt * (jj - 1)
                        donorProc = intRecvBuf(iStart + 1)
                        donorBlock = intRecvBuf(iStart + 2)
                        dIndex = intRecvBuf(iStart + 3)
                        myBlock = intRecvBuf(iStart + 4)
                        myIndex = intRecvBuf(iStart + 5)

                        if (donorBlock < 0) then
                            ! This fringe is no longer valid! Modify this
                            ! fringe in our list.
                            call unwindIndex(myIndex, il, jl, kl, i, j, k)

                            ! Now loop over the donors this cells has to match
                            ! up the one that we are dealing with here.

                            ! Loop over the potential donors.
                            iStart = fringePtr(2, i, j, k)
                            iEnd = fringePtr(3, i, j, k)

                            if (iStart == 0) then
                                print *, 'Something bad happened. The block does not think'&
                                     &'It has any donors but there is an incoming one.'
                                stop
                            end if

                            do ii = iStart, iEnd
                                if (fringes(ii)%donorProc == donorProc .and. &
                                    fringes(ii)%donorBlock == abs(donorBlock) .and. &
                                    fringes(ii)%dIndex == dIndex) then

                                    ! This is the fringe. Flag as bad.
                                    fringes(ii)%quality = 2 * Large
                                end if
                            end do
                        end if
                    end do
                end do

                ! -------------------------------------------------------------

                ! Correct irregular cells
                call irregularCellCorrection(level, sps)

                ! Set UPdate the iblank array.
                call setIblankArray(level, sps)

                ! -----------------------------------------------------------------
                ! Step 16: The algorithm is now complete. Run the checkOverset
                ! algorithm to verify that we actually have a valid interpolation
                ! -----------------------------------------------------------------
                call checkOverset(level, sps, totalOrphans, .false.)

                ! Determine if we can exit the loop. To do this we need to
                ! check if *any* of the iblank values has changed on *any*
                ! processor since the last iteration.
                localChanged = .False.
                do nn = 1, nDom
                    call setPointers(nn, level, sps)
                    do k = 2, kl
                        do j = 2, jl
                            do i = 2, il
                                if (flowDoms(nn, level, sps)%iblankLast(i, j, k) /= iblank(i, j, k)) then
                                    localChanged = .True.
                                end if
                                ! Now save the value
                                flowDoms(nn, level, sps)%iblankLast(i, j, k) = iblank(i, j, k)
                            end do
                        end do
                    end do
                end do
                globalChanged = .False.
                call mpi_allreduce(localChanged, globalChanged, 1, mpi_logical, MPI_LOR, &
                                   adflow_comm_world, ierr)
                call ECHK(ierr, __FILE__, __LINE__)

                if (.not. globalChanged) then
                    ! We're done!
                    exit refineLoop
                end if

            end do refineLoop

            if (globalChanged) then
                if (myID == 0) then
                    print "(2(A), I4, A)", "Warning: The overset connectivity loop exited ", &
                        "before the connectivity was complete after running: ", nRefine, " iterations."
                    print "(a)", "         Increase the number of iterations by setting nRefine."
                end if
            end if
            ! Final operations after the interpolations have
            ! stabilized. Status must be exchanged due to the last
            ! irregular cell correction.
            call tic(iFringeReduction)
            call exchangeStatus(level, sps, commPatternCell_2nd, internalCell_2nd)
            call fringeReduction(level, sps)
            call toc(iFringeReduction)

            call exchangeStatus(level, sps, commPatternCell_2nd, internalCell_2nd)
            call exchangeFringes(level, sps, commPatternCell_2nd, internalCell_2nd)

            call tic(iFinalCommStructures)
            call finalOversetCommStructures(level, sps)
            call setIblankArray(level, sps)
            call checkOverset(level, sps, totalOrphans, .True.)

            do nn = 1, nDom
                call setPointers(nn, level, sps)
                allocate (flowDoms(nn, level, sps)%gInd(8, 0:ib, 0:jb, 0:kb))
                flowDoms(nn, level, sps)%gInd = -1
            end do

            ! Setup the buffer sizes
            call setBufferSizes(level, sps, .false., .True.)

            ! Set up the gInd using the final overset comm structure.
            call setupFringeGlobalInd(level, sps)

            ! Deallocate some data we no longer need
            deallocate (Xmin, Xmax, work)
            do nn = 1, nDom
                deallocate (flowDoms(nn, level, sps)%iBlankLast)
            end do
            call toc(iFinalCommStructures)

            ! Done with the oFringes
            call deallocateOFringes(oFringes, size(oFringes))

            ! Ditch the temporary receiving memory and the oFringe array itself.
            deallocate (oFringes, fringeRecvSizes, cumFringeRecv, intRecvBuf)

        end do spectralLoop

        ! Add fail flag if we get orphans
        if (totalOrphans > 0) then
            if (myID == 0) then
                print *, 'Orphans present in the grid. Setting fail flags to True'
            end if
            call returnFail("oversetComm", "Orphans present in grid.")
            call mpi_barrier(ADflow_comm_world, ierr)
        end if

        ! Free the buffer and make a new one that includes necessary sizes
        ! for the overset comm
        deallocate (sendBuffer, recvBuffer)
        allocate (sendBuffer(sendBufferSize), recvBuffer(recvBufferSize))

        call MPI_barrier(adflow_comm_world, ierr)
        call ECHK(ierr, __FILE__, __LINE__)
        call toc(iTotal)

    contains

        subroutine computeDomainBoundingBoxes

            implicit none

            ! Working Variables
            real(kind=realType), dimension(3, nDom) :: xMinLocal, xMaxLocal

            xMinLocal = huge(1.0d0)
            xMaxLocal = -huge(1.0d0)
            do nn = 1, nDom
                call setPointers(nn, level, sps)
                do k = 1, kl
                    do j = 1, jl
                        do i = 1, il
                            do iDim = 1, 3
                                xMinLocal(iDim, nn) = &
                                    min(xMinLocal(iDim, nn), x(i, j, k, iDim))
                                xMaxLocal(iDim, nn) = &
                                    max(xMaxLocal(iDim, nn), x(i, j, k, iDim))
                            end do
                        end do
                    end do
                end do
            end do

            ! Now we can allgather the xMin and xMax  from each
            ! processor to everyone
            call mpi_allgatherV(xMinLocal, nDom * 3, adflow_real, xMin, 3 * nDomProc, &
                                3 * cumDomProc, adflow_real, adflow_comm_world, ierr)
            call ECHK(ierr, __FILE__, __LINE__)

            call mpi_allgatherV(xMaxLocal, nDom * 3, adflow_real, xMax, 3 * nDomProc, &
                                3 * cumDomProc, adflow_real, adflow_comm_world, ierr)
            call ECHK(ierr, __FILE__, __LINE__)

        end subroutine computeDomainBoundingBoxes

        subroutine buildGLobalSparseOverlap(overlap)

            use oversetData, only: clusters
            implicit none

            ! Input/Output
            type(CSRMatrix), intent(inout) :: overlap

            ! Working
            integer(kind=intType), dimension(:), allocatable :: colIndLocal, rowPtrLocal
            logical, dimension(:, :), allocatable :: localOverlap
            integer(kind=intType) :: nnzLocal
            integer(kind=intType), dimension(:), allocatable :: nnzProc, cumNNzProc

            ! Allocate the space for my processor's compoent of iOverlap. The
            ! number of rows is the number of blocks I own (in my own numbering)
            ! and the coluns are the total number of blocks. It is initialized
            ! to True, and remains that way unless we can conclusively prove the
            ! blocks don't overlap.
            allocate (localOverlap(nDom, nDomTotal))
            localOverlap = .True.

            ! Assume all rows are dense, and decrement when we eliminate an
            ! entry
            nnzLocal = nDom * nDomTotal

            ! Bounding box check: We will check our *own* blocks against the
            ! bounding boxes for for each domain.  This will fill up the
            ! localOverlapArray and determine nnzLocal -- the number of
            ! non-zero local entries in the global sparse matrix.

            do nn = 1, nDom
                iDom = cumDomProc(myid) + nn

                ! Now Loop over *all* of the other blocks
                do jDom = 1, nDomTotal

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
                            nnzLocal = nnzLocal - 1
                        end if
                    end if
                end do
            end do

            ! Now, create a sparse matrix representation of the local part

            allocate (colIndLocal(nnzLocal), rowPtrLocal(nDom + 1))
            rowPtrLocal(1) = 1
            i = 0
            do nn = 1, nDom
                do jDom = 1, nDomTotal
                    if (localOverlap(nn, jDom)) then
                        i = i + 1
                        colIndLocal(i) = jDom
                    end if
                end do
                rowPtrLocal(nn + 1) = i + 1
            end do

            ! Now we want to assemble the global (sparse!) connectivity
            ! matrix for all processors. This is going to require a little
            ! communication of sizes first followed by the actual data.

            ! Determine distribution of non-zero locations
            allocate (nnzProc(nProc), cumNNZProc(0:nProc))
            call mpi_allgather(nnzLocal, 1, adflow_integer, nnzProc, 1, adflow_integer, &
                               adflow_comm_world, ierr)

            overlap%nnz = sum(nnzProc)
            overlap%nRow = nDomTotal
            overlap%nCol = nDomTotal
            overlap%nnzLocal = nnzLocal
            ! We can now know how big the global data will be
            allocate ( &
                overlap%data(overlap%nnz), &
                overlap%colInd(overlap%nnz), &
                overlap%rowPtr(overlap%nRow + 1), &
                overlap%assignedProc(overlap%nnz))

            overlap%allocated = .True.
            cumNNZProc(0) = 0
            do iProc = 1, nProc
                cumNNZProc(iProc) = cumNNZProc(iProc - 1) + nnzProc(iProc)
            end do

            ! Correct for rowPtrLocal for the cumulative number of nnz up to
            ! this proc:
            rowPtrLocal = rowPtrLocal + cumNNZProc(myid)

            ! Gather the column indicies
            call mpi_allgatherV(colIndLocal, nnzLocal, adflow_integer, overlap%colInd, &
                                nnzProc, cumNNZProc, adflow_integer, adflow_comm_world, ierr)

            ! Now we gather the rowPtr to everone
            overlap%rowPtr(1) = 1
            call mpi_allgatherV(rowPtrLocal(2:nDom), nDom, adflow_integer, &
                                overlap%rowPtr(2:nDomTotal + 1), &
                                nDomProc, cumDomProc, adflow_integer, adflow_comm_world, ierr)

            ! Initialize the assignedProc to the owned rows of a processor
            do iProc = 0, nProc - 1
                overlap%assignedProc(cumNNZProc(iProc) + 1:cumNNZProc(iProc + 1)) = iProc
            end do

            deallocate (colIndLocal, rowPtrLocal, localOverlap, cumNNzProc)

        end subroutine buildGlobalSparseOverlap

        subroutine doMyWork(flag)

            ! This internal subroutine which may be called repeadly, performs
            ! as many of the searches my processor is responsible for. It
            ! returns true when all the work has been completed.

            implicit none
            logical, intent(out) :: flag
            integer(kind=intType) :: iDom, jDom, jj

            flag = .True.
            do iWork = 1, nWork

                iDom = work(1, iWork)
                jDom = work(2, iWork)
                jj = work(3, iWork) ! Index from the overlap matrix

                ! Check if I have the oBlock and fringes i need to do this
                ! intersection and I haven't already done it.
                if (oBlockReady(iDom) .and. oFringeReady(jDom) .and. &
                    work(4, iWork) == 0) then

                    startTime = mpi_wtime()

                    call fringeSearch(oBlocks(iDom), oFringes(jDom))
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

    subroutine writePartitionedMesh(fileName)

        ! This is a debugging routine for writing out meshes *as they are
        ! partioned*. This can be useful for debugging overset issues. Only
        ! the grid coordinates are written...these will have to be post
        ! processed to get connectivity information if the grid is to be
        ! used as input again.

        use constants
        use communication, only: adflow_comm_world, myID, nProc
        use blockPointers, only: il, jl, kl, nx, ny, nz, x, nDom
        use utils, only: EChk, setPointers
        use su_cgns
        implicit none

        character(len=*), intent(in) :: fileName
        integer(kind=intType) :: nDomTotal, iProc, nn, i, j, k, iDim, iDom, ierr, ii
        integer(kind=intType) :: iii, jjj, kkk
        integer(kind=intType) :: bufSize, maxSize, ibufSize, imaxSize
        integer(kind=intType), dimension(3, nDom) :: localDim
        integer(kind=intType), dimension(:), allocatable :: nDomProc, cumDomProc
        integer(kind=intType), dimension(:, :), allocatable :: dims
        real(kind=realType), dimension(:), allocatable :: buffer
        real(kind=realType), dimension(:, :, :, :), allocatable :: xtmp
        integer(kind=intType) :: ier, zoneCOunter, base, zoneID, coordID, cg, zone
        integer(kind=cgsize_t) :: sizes(9)
        integer(kind=intType) :: ifield, iSol
        character(len=40) :: tmpStr, zoneName
        character(len=32) :: coorNames(3)
        integer mpiStatus(MPI_STATUS_SIZE)
        character(len=maxStringLen) :: zoneProcFormat = "(A, I5.5, A, I3.3)"

        coorNames(1) = "CoordinateX"
        coorNames(2) = "CoordinateY"
        coorNames(3) = "CoordinateZ"

        call MPI_BARRIER(adflow_comm_world, ierr)

        ! Gather the dimensions of all blocks to everyone
        call mpi_allreduce(nDom, nDomTotal, 1, adflow_integer, MPI_SUM, &
                           adflow_comm_world, ierr)
        call ECHK(ierr, __FILE__, __LINE__)

        ! Store the sizes of the local blocks
        do nn = 1, nDom
            call setPointers(nn, 1, 1)

            ! Store the 'l' sizes for transferring
            localDim(1, nn) = il
            localDim(2, nn) = jl
            localDim(3, nn) = kl
        end do

        ! Allocate the space we need for the numbers and cumulative form
        allocate (nDomProc(0:nProc - 1), cumDomProc(0:nProc), dims(3, nDomTotal))

        ! Receive the number of domains from each proc using an allgather.
        call mpi_allgather(nDom, 1, adflow_integer, nDomProc, 1, adflow_integer, &
                           adflow_comm_world, ierr)
        call ECHK(ierr, __FILE__, __LINE__)

        ! Compute the cumulative format:
        cumDomProc(0) = 0
        do iProc = 1, nProc
            cumDomProc(iProc) = cumDomProc(iProc - 1) + nDomProc(iProc - 1)
        end do

        ! We will also allgather all of the block sizes which will make
        ! things a little easier since everyone will know the proper sizes
        ! for the sends
        call mpi_allgatherV(localDim, nDom * 3, adflow_integer, dims, 3 * nDomProc, &
                            3 * cumDomProc, adflow_integer, adflow_comm_world, ierr)
        call ECHK(ierr, __FILE__, __LINE__)

        maxSize = 0
        do i = 1, nDomTotal
            maxSize = max(maxSize, dims(1, i) * dims(2, i) * dims(3, i) * 3)
        end do

        allocate (buffer(maxSize))

        if (myid == 0) then

            ! Open the CGNS File
            call cg_open_f(fileName, mode_write, cg, ier)
            base = 1
            call cg_base_write_f(cg, "Base#1", 3, 3, base, ier)

            zoneCounter = 0
            ! Write my own blocks first
            do nn = 1, nDom
                call setPointers(nn, 1, 1)

                sizes(1) = il
                sizes(2) = jl
                sizes(3) = kl
                sizes(4) = nx
                sizes(5) = ny
                sizes(6) = nz
                sizes(7) = 0
                sizes(8) = 0
                sizes(9) = 0

                zoneCounter = zoneCounter + 1
                write (zonename, zoneProcFormat) 'domain.', zoneCounter, 'proc.', myid

                call cg_zone_write_f(cg, base, zonename, sizes, Structured, zoneID, ier)

                allocate (xtmp(sizes(1), sizes(2), sizes(3), 3))

                do k = 1, kl
                    do j = 1, jl
                        do i = 1, il
                            xtmp(i, j, k, 1:3) = x(i, j, k, 1:3)
                        end do
                    end do
                end do

                do idim = 1, 3
                    call cg_coord_write_f(cg, base, zoneID, realDouble, coorNames(idim), &
                                          xtmp(:, :, :, idim), coordID, ier)
                end do
                deallocate (xtmp)
            end do

            ! Now loop over the remaining blocks...receiving each and writing:

            do iProc = 1, nProc - 1
                do nn = 1, nDomProc(iProc)
                    iDom = cumDomProc(iProc) + nn
                    bufSize = dims(1, iDom) * dims(2, iDom) * dims(3, iDom) * 3

                    call MPI_Recv(buffer, bufSize, adflow_real, iProc, iProc, &
                                  adflow_comm_world, mpiStatus, ierr)

                    zoneCounter = zoneCounter + 1
                    write (zonename, zoneProcFormat) 'domain.', zoneCounter, 'proc.', iProc
                    sizes(1) = dims(1, iDom)
                    sizes(2) = dims(2, iDom)
                    sizes(3) = dims(3, iDom)
                    sizes(4) = dims(1, iDom) - 1
                    sizes(5) = dims(2, iDom) - 1
                    sizes(6) = dims(3, iDom) - 1
                    sizes(7) = 0
                    sizes(8) = 0
                    sizes(9) = 0
                    call cg_zone_write_f(cg, base, zonename, sizes, Structured, zoneID, ier)
                    ii = 0
                    allocate (xtmp(sizes(1), sizes(2), sizes(3), 3))
                    do k = 1, sizes(3)
                        do j = 1, sizes(2)
                            do i = 1, sizes(1)
                                xtmp(i, j, k, 1) = buffer(ii + 1)
                                xtmp(i, j, k, 2) = buffer(ii + 2)
                                xtmp(i, j, k, 3) = buffer(ii + 3)
                                ii = ii + 3
                            end do
                        end do
                    end do

                    do idim = 1, 3
                        call cg_coord_write_f(cg, base, zoneID, realDouble, coorNames(idim), &
                                              xtmp(:, :, :, idim), coordID, ier)
                    end do
                    deallocate (xtmp)
                end do
            end do
        else
            ! Pack and send my stuff:
            do nn = 1, nDom
                call setPointers(nn, 1, 1)
                ii = 0
                do k = 1, kl
                    do j = 1, jl
                        do i = 1, il
                            do iDim = 1, 3
                                buffer(ii + idim) = x(i, j, k, iDim)
                            end do
                            ii = ii + 3
                        end do
                    end do
                end do

                call mpi_send(buffer, ii, adflow_real, 0, myid, &
                              adflow_comm_world, ierr)
            end do
        end if

        deallocate (buffer, nDomProc, cumDomProc, dims)

    end subroutine writePartitionedMesh

    ! Debugging routine for writing the dual grids along with the volume
    ! in CGNS format to help debug.

    subroutine writeDualMesh(fileName)

        ! This is a debugging routine for writing out meshes *as they are
        ! partioned*. This can be useful for debugging overset issues.
        use constants
        use communication, only: adflow_comm_world, myid, nProc
        use blockPointers, only: ie, je, ke, il, jl, kl, x, globalCell, vol, &
                                 nDom, iblank
        use utils, only: setPointers, EChk
        use su_cgns
        implicit none

        character(len=*), intent(in) :: fileName
        integer(kind=intType) :: nDomTotal, iProc, nn, i, j, k, iDim, iDom, ierr, ii
        integer(kind=intType) :: iii, jjj, kkk
        integer(kind=intType) :: bufSize, maxSize, ibufSize, imaxSize
        integer(kind=intType), dimension(3, nDom) :: localDim
        integer(kind=intType), dimension(:), allocatable :: nDomProc, cumDomProc
        integer(kind=intType), dimension(:, :), allocatable :: dims
        real(kind=realType), dimension(:), allocatable :: buffer
        real(kind=realType), dimension(:, :, :, :), allocatable :: xtmp
        integer(kind=intType) :: ier, zoneCOunter, base, zoneID, coordID, cg, zone
        integer(kind=cgsize_t) :: sizes(9)
        integer(kind=intType) :: ifield, iSol
        character(len=40) :: tmpStr, zoneName
        character(len=32) :: coorNames(3)
        integer mpiStatus(MPI_STATUS_SIZE)
        character(len=maxStringLen) :: zoneFormat = "(A, I5.5)"

        coorNames(1) = "CoordinateX"
        coorNames(2) = "CoordinateY"
        coorNames(3) = "CoordinateZ"

        ! Gather the dimensions of all blocks to everyone
        call mpi_allreduce(nDom, nDomTotal, 1, adflow_integer, MPI_SUM, &
                           adflow_comm_world, ierr)
        call ECHK(ierr, __FILE__, __LINE__)

        ! Store the sizes of the local blocks
        do nn = 1, nDom
            call setPointers(nn, 1, 1)
            localDim(1, nn) = ie
            localDim(2, nn) = je
            localDim(3, nn) = ke
        end do

        ! Allocate the space we need for the numbers and cumulative form
        allocate (nDomProc(0:nProc - 1), cumDomProc(0:nProc), dims(3, nDomTotal))

        ! Receive the number of domains from each proc using an allgather.
        call mpi_allgather(nDom, 1, adflow_integer, nDomProc, 1, adflow_integer, &
                           adflow_comm_world, ierr)
        call ECHK(ierr, __FILE__, __LINE__)

        ! Compute the cumulative format:
        cumDomProc(0) = 0
        do iProc = 1, nProc
            cumDomProc(iProc) = cumDomProc(iProc - 1) + nDomProc(iProc - 1)
        end do

        ! We will also allgather all of the block sizes which will make
        ! things a little easier since everyone will know the proper sizes
        ! for the sends
        call mpi_allgatherV(localDim, nDom * 3, adflow_integer, dims, 3 * nDomProc, &
                            3 * cumDomProc, adflow_integer, adflow_comm_world, ierr)
        call ECHK(ierr, __FILE__, __LINE__)

        maxSize = 0
        do i = 1, nDomTotal
            maxSize = max(maxSize, dims(1, i) * dims(2, i) * dims(3, i) * 5)
        end do

        allocate (buffer(maxSize))

        if (myid == 0) then
            call cg_open_f(fileName, mode_write, cg, ier)
            base = 1
            call cg_base_write_f(cg, "Base#1", 3, 3, base, ier)

            zoneCounter = 0
            ! Write my own blocks first
            do nn = 1, nDom
                call setPointers(nn, 1, 1)

                sizes(1) = ie
                sizes(2) = je
                sizes(3) = ke
                sizes(4) = il
                sizes(5) = jl
                sizes(6) = kl
                sizes(7) = 0
                sizes(8) = 0
                sizes(9) = 0

                zoneCounter = zoneCounter + 1
                write (zonename, zoneFormat) 'domain.', zoneCounter

                call cg_zone_write_f(cg, base, zonename, sizes, Structured, zoneID, ier)

                allocate (xtmp(sizes(1), sizes(2), sizes(3), 5))

                do k = 1, ke
                    do j = 1, je
                        do i = 1, ie
                            xtmp(i, j, k, 1:3) = eighth * ( &
                                                 x(i - 1, j - 1, k - 1, :) + &
                                                 x(i, j - 1, k - 1, :) + &
                                                 x(i - 1, j, k - 1, :) + &
                                                 x(i, j, k - 1, :) + &
                                                 x(i - 1, j - 1, k, :) + &
                                                 x(i, j - 1, k, :) + &
                                                 x(i - 1, j, k, :) + &
                                                 x(i, j, k, :))
                            xtmp(i, j, k, 4) = vol(i, j, k)
                            if (globalCell(i, j, k) >= 0) then
                                xtmp(i, j, k, 5) = dble(iblank(i, j, k))
                            else
                                xtmp(i, j, k, 5) = zero
                            end if
                        end do
                    end do
                end do

                do idim = 1, 3
                    call cg_coord_write_f(cg, base, zoneID, realDouble, coorNames(idim), &
                                          xtmp(:, :, :, idim), coordID, ier)
                end do

                call cg_sol_write_f(cg, base, zoneID, "flowSolution", Vertex, iSol, ier)

                call cg_field_write_f(cg, base, zoneID, iSol, realDouble, "volume", &
                                      xtmp(:, :, :, 4), iField, ier)

                call cg_field_write_f(cg, base, zoneID, iSol, realDouble, "iBlank", &
                                      xtmp(:, :, :, 5), iField, ier)

                deallocate (xtmp)
            end do

            ! Now loop over the remaining blocks...receiving each and writing:

            do iProc = 1, nProc - 1
                do nn = 1, nDomProc(iProc)
                    iDom = cumDomProc(iProc) + nn
                    bufSize = dims(1, iDom) * dims(2, iDom) * dims(3, iDom) * 5

                    call MPI_Recv(buffer, bufSize, adflow_real, iProc, iProc, &
                                  adflow_comm_world, mpiStatus, ierr)

                    zoneCounter = zoneCounter + 1
                    write (zonename, zoneFormat) 'domain.', zoneCounter
                    sizes(1) = dims(1, iDom)
                    sizes(2) = dims(2, iDom)
                    sizes(3) = dims(3, iDom)
                    sizes(4) = dims(1, iDom) - 1
                    sizes(5) = dims(2, iDom) - 1
                    sizes(6) = dims(3, iDom) - 1
                    sizes(7) = 0
                    sizes(8) = 0
                    sizes(9) = 0
                    call cg_zone_write_f(cg, base, zonename, sizes, Structured, zoneID, ier)
                    ii = 0
                    allocate (xtmp(sizes(1), sizes(2), sizes(3), 5))
                    do k = 1, sizes(3)
                        do j = 1, sizes(2)
                            do i = 1, sizes(1)
                                xtmp(i, j, k, 1) = buffer(ii + 1)
                                xtmp(i, j, k, 2) = buffer(ii + 2)
                                xtmp(i, j, k, 3) = buffer(ii + 3)
                                xtmp(i, j, k, 4) = buffer(ii + 4)
                                xtmp(i, j, k, 5) = buffer(ii + 5)
                                ii = ii + 5
                            end do
                        end do
                    end do

                    do idim = 1, 3
                        call cg_coord_write_f(cg, base, zoneID, realDouble, coorNames(idim), &
                                              xtmp(:, :, :, idim), coordID, ier)
                    end do

                    call cg_sol_write_f(cg, base, zoneID, "flowSolution", Vertex, iSol, ier)
                    call cg_field_write_f(cg, base, zoneID, iSol, realDouble, "volume", &
                                          xtmp(:, :, :, 4), iField, ier)
                    call cg_field_write_f(cg, base, zoneID, iSol, realDouble, "iBlank", &
                                          xtmp(:, :, :, 5), iField, ier)

                    deallocate (xtmp)
                end do
            end do

        else

            ! Pack and send my stuff:
            do nn = 1, nDom
                call setPointers(nn, 1, 1)
                ii = 0
                do k = 1, ke
                    do j = 1, je
                        do i = 1, ie
                            do iDim = 1, 3
                                buffer(ii + idim) = eighth * ( &
                                                    x(i - 1, j - 1, k - 1, idim) + &
                                                    x(i, j - 1, k - 1, idim) + &
                                                    x(i - 1, j, k - 1, idim) + &
                                                    x(i, j, k - 1, idim) + &
                                                    x(i - 1, j - 1, k, idim) + &
                                                    x(i, j - 1, k, idim) + &
                                                    x(i - 1, j, k, idim) + &
                                                    x(i, j, k, idim))
                            end do
                            buffer(ii + 4) = vol(i, j, k)
                            if (globalCell(i, j, k) > 0) then
                                buffer(ii + 5) = dble(iblank(i, j, k))
                            else
                                buffer(ii + 5) = zero
                            end if

                            ii = ii + 5
                        end do
                    end do
                end do

                call mpi_send(buffer, ii, adflow_real, 0, myid, &
                              adflow_comm_world, ierr)
            end do
        end if

        deallocate (buffer, nDomProc, cumDomProc, dims)

    end subroutine writeDualMesh

    !
    !       determineClusters determines which blocks are connected with
    !       1to1 cgns connections. Essentially what we are doing is
    !       identifying the consitutive multiblock meshes that make up
    !       an overset mesh. There should be precisely 1 cluster for a
    !       face mached mesh and 2 or more for a overset mesh

    subroutine determineClusters()

        use constants
        use blockPointers, only: nDom, flowDoms
        use cgnsGrid, only: CGNSDoms, cgnsNDom
        use communication, only: adflow_comm_world, myID
        use oversetData, only: clusters, nDomTotal, nClusters, cumDomProc
        implicit none

        ! Working variables
        integer(kind=intType) :: numBlocks, blockID, cgnsBlk, ierr, clusterID
        integer(kind=intType) :: i, nn
        integer(kind=intType), dimension(nDomTotal) :: clustersLocal
        logical :: blocksAvailable

        ! Initialize the cluster of each of the CGNSDoms to 0
        do i = 1, cgnsNDom
            cgnsDoms(i)%cluster = 0
        end do

        ! Allocate clusters (defined in overset)
        allocate (clusters(nDomTotal))

        ! Initialize cluster counter
        clusterID = 0

        ! Initialize counter of classified blocks
        blockID = 0

        ! Initialize variable to state that we have unclassified blocks
        blocksAvailable = .True.

        ! Loop until all blocks are checked
        do while (blocksAvailable)

            ! Find position of the available block
            blocksAvailable = .false.

            do while ((.not. blocksAvailable) .and. (blockID .lt. cgnsnDom))
                blockID = blockID + 1 ! Increment counter
                if (cgnsDoms(blockID)%cluster == 0) then
                    blocksAvailable = .true.
                end if
            end do

            ! If we have blocks available, we start the search
            if (blocksAvailable) then
                clusterID = clusterID + 1 ! Increment the running cluser counter
                cgnsDoms(blockID)%cluster = clusterID
                call clusterSearch(blockID)
            end if

        end do

        ! Set the clusters to 0 so we can just all reduce
        clustersLocal = 0

        ! Set the cluster ID for all my blocks:
        do nn = 1, nDom
            cgnsBlk = flowDoms(nn, 1, 1)%cgnsBlockID
            clustersLocal(cumDomProc(myid) + nn) = cgnsDoms(cgnsBlk)%cluster
        end do
        call MPI_Allreduce(clustersLocal, clusters, nDomTotal, adflow_integer, MPI_SUM, &
                           adflow_comm_world, ierr)

        ! Finally, set the total number of clusters
        nClusters = clusterID

    contains

        recursive subroutine clusterSearch(blockID)

            ! This is the recursive part of cluster search
            implicit none

            ! Subroutine inputs
            integer(kind=intType), intent(in) :: blockID

            ! Working variables
            integer(kind=intTYpe) :: clusterID, connID, connBlock

            ! Get the cluster ID from the reference block
            clusterID = cgnsDoms(blockID)%cluster

            ! Loop over all connections of this block
            do connID = 1, cgnsDoms(blockID)%n1to1

                connBlock = cgnsDoms(blockID)%conn1to1(connID)%donorBlock

                ! Check if connected block is already classified
                if (cgnsDoms(connBlock)%cluster == 0) then
                    cgnsDoms(connBlock)%cluster = clusterID ! Assign block to the same cluster
                    call clusterSearch(connBlock) ! Start search on the new block

                else if (cgnsDoms(connBlock)%cluster .ne. clusterID) then ! Check symmetry
                    print *, 'Non-symmetric connection between CGNS blocks:', blockID, ' and', connBlock
                    stop
                end if
            end do

        end subroutine clusterSearch
    end subroutine determineClusters

    subroutine determineViscousDirs()

        ! Set the viscousDir flags in the CGNS grid based on the CGNS grid
        ! boundary conditions
        use constants
        use cgnsGrid, only: CGNSDoms, cgnsNDom
        use communication, only: myid
        implicit none

        ! Working variables
        integer(kind=intType) :: i, j, bc

        do i = 1, cgnsNDom
            do j = 1, cgnsDoms(i)%nBocos
                bc = cgnsDoms(i)%bocoInfo(j)%BCType

                if (bc == NSWallAdiabatic .or. bc == NSWallIsoThermal .or. bc == EulerWall) then

                    if (cgnsDoms(i)%bocoInfo(j)%iBeg == cgnsDoms(i)%bocoInfo(j)%iEnd) &
                        cgnsDoms(i)%viscousDir(1) = .True.

                    if (cgnsDoms(i)%bocoInfo(j)%jBeg == cgnsDoms(i)%bocoInfo(j)%jEnd) &
                        cgnsDoms(i)%viscousDir(2) = .True.

                    if (cgnsDoms(i)%bocoInfo(j)%kBeg == cgnsDoms(i)%bocoInfo(j)%kEnd) &
                        cgnsDoms(i)%viscousDir(3) = .True.
                end if
            end do
        end do
    end subroutine determineViscousDirs

    subroutine setExplicitHoleCut(flag)

        ! This is meant to be the gateway for doing any explict hole
        ! cutting. Right now we have a call-back approach

        use constants
        use adjointVars, only: nCellsGlobal
        use inputTimeSpectral, only: nTimeIntervalsSpectral
        use blockPointers, only: nDom, il, jl, kl, x, iBlank, flowDoms
        use utils, only: setPointers
        use communication, only: commPatternCell_2nd, internalCell_2nd
        use haloExchange, only: whalo1to1IntGeneric
        use adjointvars, only: nCellsLocal
        implicit none

        ! Input/output
        integer(kind=intType), dimension(:), intent(in) :: flag

        ! Working
        integer(kind=intType) :: i, j, k, ii, nn, sps, level

        level = 1
        ! Set iblank to -4 if the flag is true:
        ii = 0
        do nn = 1, nDom
            do sps = 1, nTimeIntervalsSpectral
                call setPointers(nn, 1, sps)

                do k = 2, kl
                    do j = 2, jl
                        do i = 2, il
                            ii = ii + 1
                            if (flag(ii) /= 0) then
                                iblank(i, j, k) = -4
                            end if
                        end do
                    end do
                end do
            end do
        end do

        do sps = 1, nTimeIntervalsSpectral
            ! Exchange iblanks
            domainLoop: do nn = 1, nDom
                flowDoms(nn, level, sps)%intCommVars(1)%var => &
                    flowDoms(nn, level, sps)%iblank(:, :, :)
            end do domainLoop

            ! Run the generic integer exchange
            call wHalo1to1IntGeneric(1, level, sps, commPatternCell_2nd, internalCell_2nd)
        end do

    end subroutine setExplicitHoleCut

    subroutine flagCellsInSurface(pts, npts, conn, nconn, flag, ncell, blockids, nblocks, k_min)

        use constants
        use communication, only: myID
        use adtBuild, only: buildSerialQuad, destroySerialQuad
        use adtLocalSearch, only: minDistanceTreeSearchSinglePoint
        use ADTUtils, only: stack
        use ADTData
        use blockPointers, only: x, il, jl, kl, nDom, iBlank, vol, nbkGlobal, kBegOr
        use adjointVars, only: nCellsLocal
        use utils, only: setPointers, EChk
        use sorting, only: famInList
        implicit none

        ! Input/Output
        integer(kind=intType), intent(in) :: npts, nconn, ncell, nblocks
        real(kind=realType), intent(in), dimension(3, npts) :: pts
        integer(kind=intType), intent(in), dimension(4, nconn) :: conn
        integer(kind=intType), intent(inout), dimension(ncell) :: flag
        integer(kind=intType), intent(in), dimension(nblocks) :: blockids
        integer(kind=intType), intent(in) :: k_min

        ! Working variables
        integer(kind=intType) :: i, j, k, l, nn, iDim, cellID, intInfo(3), sps, level, iii, ierr
        integer(kind=intType) :: iblock, cell_counter, k_cgns
        real(kind=realType) :: dStar, frac, volLocal
        real(kind=realType), dimension(3) :: minX, maxX, v1, v2, v3, axisVec
        real(kind=realType), dimension(3) :: diag1, diag2, diag3, diag4
        real(kind=realType) :: dd1, dd2, dd3, dd4, diag_max
        type(adtType) :: ADT
        real(kind=realType), dimension(:, :), allocatable :: norm
        integer(kind=intType), dimension(:), allocatable :: normCount
        integer(kind=intType), dimension(:, :), pointer :: tmp

        ! ADT Type required data
        integer(kind=intType), dimension(:), pointer :: frontLeaves, frontLeavesNew
        type(adtBBoxTargetType), dimension(:), pointer :: BB
        real(kind=realType) :: coor(4), uvw(5)
        real(kind=realType) :: dummy(3, 2)

        ! we use the same approach as the actuator zone addition; we project the cell centers to the
        ! surface provided. This way, we can fairly quickly determine the coordinates
        ! that are inside the closed volumes. The surface mesh is duplicated on all procs, whereas
        ! the nodes projected are distributed based on the grid paritioning.

        ! Since this is effectively a wall-distance calc it gets super
        ! costly for the points far away. Luckly, we can do a fairly
        ! simple shortcut: Just compute the bounding box of the region and
        ! use that as the "already found" distance in the closest point
        ! search. This will eliminate all the points further away
        ! immediately and this should be sufficiently fast.

        ! So...compute that bounding box:
        do iDim = 1, 3
            minX(iDim) = minval(pts(iDim, :))
            maxX(iDim) = maxval(pts(iDim, :))
        end do

        ! Get the max distance. This should be quite conservative.
        dStar = (maxX(1) - minx(1))**2 + (maxX(2) - minX(2))**2 + (maxX(3) - minX(3))**2

        ! Now build the tree.
        call buildSerialQuad(size(conn, 2), size(pts, 2), pts, conn, ADT)

        ! Compute the (averaged) unique nodal vectors:
        allocate (norm(3, size(pts, 2)), normCount(size(pts, 2)))

        norm = zero
        normCount = 0

        do i = 1, size(conn, 2)

            ! Compute cross product normal and normalize
            v1 = pts(:, conn(3, i)) - pts(:, conn(1, i))
            v2 = pts(:, conn(4, i)) - pts(:, conn(2, i))

            v3(1) = (v1(2) * v2(3) - v1(3) * v2(2))
            v3(2) = (v1(3) * v2(1) - v1(1) * v2(3))
            v3(3) = (v1(1) * v2(2) - v1(2) * v2(1))
            v3 = v3 / sqrt(v3(1)**2 + v3(2)**2 + v3(3)**2)

            ! Add to each of the four pts and increment the number added
            do j = 1, 4
                norm(:, conn(j, i)) = norm(:, conn(j, i)) + v3
                normCount(conn(j, i)) = normCount(conn(j, i)) + 1
            end do
        end do

        ! Now just divide by the norm count
        do i = 1, size(pts, 2)
            norm(:, i) = norm(:, i) / normCount(i)
        end do

        ! Norm count is no longer needed
        deallocate (normCount)

        ! Allocate the extra data the tree search requires.
        allocate (stack(100), BB(20), frontLeaves(25), frontLeavesNew(25))

        ! Now search for all the coordinate. Note that We have explictly
        ! set sps to 1 becuase it is only implemented for single grid.
        sps = 1
        level = 1
        cell_counter = 1

        do nn = 1, nDom
            call setPointers(nn, level, sps)

            ! only check this cell if it is within one of the block IDs we are asked to look at
            ! here, we reuse the famInList function from sorting.F90. The reason is just having
            ! if (any(blockids == nbkGlobal)) does not work with complexify; complexify changes the
            ! == (or .eq.) to .ceq., which does not work with integers. The famInList is originally
            ! written for surface integrations, but the same idea applies here.
            if (famInList(nbkGlobal, blockids)) then
                do k = 2, kl
                    do j = 2, jl
                        do i = 2, il

                            ! get the k index of this cell in the cgns grid
                            k_cgns = k + kBegOr - 2

                            ! check if we are above the kmin range needed. if no kmin is provided, then the value
                            ! defaults to -1 from python, so all points satisfy the check
                            if (k_cgns .gt. k_min) then

                                ! calculate the 4 diagonals of the cell
                                diag1 = x(i - 1, j - 1, k - 1, :) - x(i, j, k, :)
                                diag2 = x(i, j - 1, k - 1, :) - x(i - 1, j, k, :)
                                diag3 = x(i, j, k - 1, :) - x(i - 1, j - 1, k, :)
                                diag4 = x(i - 1, j, k - 1, :) - x(i, j - 1, k, :)

                                dd1 = diag1(1) * diag1(1) + diag1(2) * diag1(2) + diag1(3) * diag1(3)
                                dd2 = diag2(1) * diag2(1) + diag2(2) * diag2(2) + diag2(3) * diag2(3)
                                dd3 = diag3(1) * diag3(1) + diag3(2) * diag3(2) + diag3(3) * diag3(3)
                                dd4 = diag4(1) * diag4(1) + diag4(2) * diag4(2) + diag4(3) * diag4(3)

                                ! get the max
                                diag_max = max(dd1, dd2, dd3, dd4)
                                ! if the projection is greater than this for any node, we stop testing the current point.
                                ! we can work with squared distances for the sake of efficiency. the projection routine
                                ! will also return the squared distance for the same reason.

                                ! actually test each node
                                do l = 1, 8
                                    select case (l)
                                    case (1)
                                        coor(1:3) = x(i - 1, j - 1, k - 1, :)
                                    case (2)
                                        coor(1:3) = x(i, j - 1, k - 1, :)
                                    case (3)
                                        coor(1:3) = x(i, j, k - 1, :)
                                    case (4)
                                        coor(1:3) = x(i - 1, j, k - 1, :)
                                    case (5)
                                        coor(1:3) = x(i - 1, j - 1, k, :)
                                    case (6)
                                        coor(1:3) = x(i, j - 1, k, :)
                                    case (7)
                                        coor(1:3) = x(i, j, k, :)
                                    case (8)
                                        coor(1:3) = x(i - 1, j, k, :)
                                    end select

                                    ! reset the "closest point already found" variable.
                                    coor(4) = dStar
                                    intInfo(3) = 0
                                    call minDistancetreeSearchSinglePoint(ADT, coor, intInfo, &
                                                                          uvw, dummy, 0, BB, &
                                                                          frontLeaves, frontLeavesNew)
                                    cellID = intInfo(3)
                                    if (cellID > 0) then
                                        ! Now check if this was successful or not:
                                        if (checkInside()) then
                                            ! Whoohoo! We are inside the region. Flag this cell
                                            flag(cell_counter) = 1
                                            exit
                                        else
                                            ! we are outside. now check if the projection distance is larger than
                                            ! the max diagonal. if so, we can quit early here.
                                            if (uvw(4) .gt. diag_max) then
                                                ! projection is larger than our biggest diagonal.
                                                ! other nodes wont be in the surface, so we can exit the cell early here
                                                exit
                                            end if
                                        end if
                                    end if
                                end do
                            end if

                            cell_counter = cell_counter + 1
                        end do
                    end do
                end do
            else
                ! we dont want to consider block, but we need to increment the counter
                cell_counter = cell_counter + (kl - 1) * (jl - 1) * (il - 1)
            end if
        end do

        ! Final memory cleanup
        deallocate (stack, norm, frontLeaves, frontLeavesNew, BB)
        call destroySerialQuad(ADT)

    contains

        function checkInside()

            implicit none
            logical :: checkInside
            integer(kind=intType) :: jj
            real(kind=realType) :: shp(4), xp(3), normal(3), v1(3), dp

            ! bi-linear shape functions (CCW ordering)
            shp(1) = (one - uvw(1)) * (one - uvw(2))
            shp(2) = (uvw(1)) * (one - uvw(2))
            shp(3) = (uvw(1)) * (uvw(2))
            shp(4) = (one - uvw(1)) * (uvw(2))

            xp = zero
            normal = zero
            do jj = 1, 4
                xp = xp + shp(jj) * pts(:, conn(jj, cellID))
                normal = normal + shp(jj) * norm(:, conn(jj, cellID))
            end do

            ! Compute the dot product of normal with cell center
            ! (stored in coor) with the point on the surface.
            v1 = coor(1:3) - xp
            dp = normal(1) * v1(1) + normal(2) * v1(2) + normal(3) * v1(3)

            if (dp < zero) then
                checkInside = .True.
            else
                checkInside = .False.
            end if
        end function checkInside

    end subroutine flagCellsInSurface

    subroutine updateOverset(flag, n, closedFamList, nFam)

        ! This is the main gateway routine for updating the overset
        ! connectivity during a solution. It is wrapped and intended to be
        ! called from Python. What this routine does depends on the value
        ! of oversetUpdateMode:

        ! updateFrozen: Nothing happens. The initial
        ! connectivity computed during initialization is kept.

        ! updateFast: Update just the weight, but leave the donors
        ! unchanged. This is only applicable when the entire mesh is
        ! warped at the same time like with USMesh in IDWarp.

        ! updateFull: Complete from scratch update. Run the full
        ! oversetComm routine.

        use constants
        use inputOverset, only: oversetUpdateMode
        use block, only: flowDOms
        use inputTimeSpectral, only: nTimeIntervalsSpectral
        use oversetCommUtilities, onlY: updateOversetConnectivity
        use oversetData, only: oversetPresent
        implicit none

        ! Input/Output
        integer(kind=intType), intent(in) :: n, nFam
        integer(kind=intType), dimension(n) :: flag
        integer(kind=intType), dimension(nFam) :: closedFamList

        ! Working
        integer(kind=intType) :: nLevels, level, sps, nn, nDoms, mm

        ! If no overset, don't do anything.
        if (.not. oversetPresent) then
            return
        end if

        select case (oversetUpdateMode)
        case (updateFrozen)
            return

        case (updateFast)
            ! Run our fast version.
            nLevels = ubound(flowDoms, 2)
            do level = 1, nLevels
                do sps = 1, nTimeIntervalsSpectral
                    call updateOversetConnectivity(level, sps)
                end do
            end do

        case (updateFull)
            ! Do the full update like we would have done during
            ! initialization.
            nLevels = ubound(flowDoms, 2)

            ! Reinitialize domain iblanks
            nDoms = ubound(flowDoms, 1)
            do nn = 1, nDoms
                do level = 1, nLevels
                    do sps = 1, nTimeIntervalsSpectral
                        flowDoms(nn, level, sps)%iBlank = 1
                        flowDoms(nn, level, sps)%forcedRecv = 0
                        flowDoms(nn, level, sps)%status = 0
                        do mm = 1, flowDoms(nn, level, sps)%nBocos
                            flowDoms(nn, level, sps)%BCData(mm)%iblank = 1
                        end do
                    end do
                end do
            end do

            do level = 1, nLevels
                if (level == 1) then
                    call setExplicitHoleCut(flag)
                    call oversetComm(level, .False., .false., closedFamList)
                else
                    call oversetComm(level, .False., .True., closedFamList)
                end if
            end do

        end select

    end subroutine updateOverset

    subroutine setBlockPriority(blkName, value, setValue)

        ! Set the CGNSblock with blkName to have a priority of "value"

        use constants
        use cgnsGrid, only: cgnsDoms, cgnsNDom
        use communication
        implicit none

        ! Input Parameters
        character(len=*), intent(in) :: blkName
        real(kind=realType), intent(in) :: value
        logical, intent(out) :: setValue

        ! Working
        integer(kind=intType) :: iDom

        ! We should do a binary search for the block names, but for now
        ! just loop over the number of CGNS Doms.
        if (cgnsNDom > 0) then
            setValue = .False.
            do iDom = 1, CGNSnDom
                if (trim(cgnsDoms(iDom)%zoneName) == trim(blkName)) then
                    setValue = .True.
                    cgnsDoms(iDom)%priority = value
                    exit
                end if
            end do
        else
            ! If there are no domains yet, just say we did it so we don't
            ! raise an error in python.
            setValue = .True.
        end if
    end subroutine setBlockPriority

end module oversetAPI

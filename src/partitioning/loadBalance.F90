module loadBalance

contains

  subroutine loadBalanceGrid
    !
    !       loadBalance determines the mapping of the blocks onto the
    !       processors. If the user allows so blocks my be split to obtain
    !       a better load balance.
    !
    use constants
    use block, only : flowDoms, nDom
    use cgnsGrid, only : cgnsDoms, cgnsNDom
    use communication, only : adflow_comm_world, nProc, myID
    use inputMotion, only : gridMotionSpecified
    use inputParallel, only : partitionLikeNProc
    use inputPhysics, only : equationMode
    use inputTimeSpectral, only : nTimeIntervalsSpectral
    use iteration, only : deforming_Grid
    use partitionMod, only : subBlocksofCGNSType, blocks, part, nBlocks, &
         sortRangesSplitInfo, qsortSubblocksOfCGNSType
    use utils, only : terminate
    implicit none
    !
    !      Local variables.
    !
    integer :: ierr

    integer(kind=intType) :: i, j, k, nn, mm, ii, jj, kk
    integer(kind=intType) :: nViscBocos

    integer(kind=intType), dimension(0:nProc-1) :: nBlockPerProc

    integer(kind=intType), dimension(:), allocatable :: oldSubfaceID

    type(subblocksOfCGNSType), dimension(:), allocatable :: &
         subblocksOfCGNS


    ! Determine the block distribution over the processors.


    if (partitionLikeNProc > nProc) then
       nProc = partitionLikenProc
    end if

    call blockDistribution

    ! Restore the size of the comm
    call mpi_Comm_Size(adflow_comm_world, nProc, ierr)

    ! We ned to modify what comes back if we are using the
    ! partitionLikenProc option:
    if (partitionLikenProc > nProc) then
       do i=1, nBlocks
          part(i) = mod(part(i), nProc)
       end do
    end if

    !
    !       Determine the local block info.
    !
    ! Initialize nBlockPerProc to 0.

    nBlockPerProc = 0

    ! Determine the number of blocks the current processor will store
    ! and the local block number for every block.

    nDom = 0
    do i=1,nBlocks
       if(part(i) == myID) nDom = nDom +1

       nBlockPerProc(part(i)) = nBlockPerProc(part(i)) + 1
       blocks(i)%blockID = nBlockPerProc(part(i))
    enddo

    ! Allocate the memory for flowDoms and initialize its pointers
    ! to null pointers.

    call initFlowDoms

    ! Repeat the loop, but now store the info of the blocks
    ! in flowDoms. Store the number of time intervals for the spectral
    ! method a bit easier in mm. Note that this number is 1 for the
    ! steady and unsteady modes.

    nn = 0
    mm = nTimeIntervalsSpectral
    domains: do i=1,nBlocks
       myBlock: if(part(i) == myID) then

          ! Update the counter nn.

          nn = nn + 1

          ! Copy the dimensions of the block.

          flowDoms(nn,1,1:mm)%nx = blocks(i)%nx
          flowDoms(nn,1,1:mm)%ny = blocks(i)%ny
          flowDoms(nn,1,1:mm)%nz = blocks(i)%nz

          flowDoms(nn,1,1:mm)%il = blocks(i)%il
          flowDoms(nn,1,1:mm)%jl = blocks(i)%jl
          flowDoms(nn,1,1:mm)%kl = blocks(i)%kl

          ! The number of single halo quantities.

          flowDoms(nn,1,1:mm)%ie = blocks(i)%il + 1
          flowDoms(nn,1,1:mm)%je = blocks(i)%jl + 1
          flowDoms(nn,1,1:mm)%ke = blocks(i)%kl + 1

          ! The number of double halo quantities.

          flowDoms(nn,1,1:mm)%ib = blocks(i)%il + 2
          flowDoms(nn,1,1:mm)%jb = blocks(i)%jl + 2
          flowDoms(nn,1,1:mm)%kb = blocks(i)%kl + 2

          ! Relation to the original cgns grid.

          flowDoms(nn,1,1:mm)%cgnsBlockID = blocks(i)%cgnsBlockID

          flowDoms(nn,1,1:mm)%iBegOr = blocks(i)%iBegOr
          flowDoms(nn,1,1:mm)%jBegOr = blocks(i)%jBegOr
          flowDoms(nn,1,1:mm)%kBegOr = blocks(i)%kBegOr

          flowDoms(nn,1,1:mm)%iEndOr = blocks(i)%iEndOr
          flowDoms(nn,1,1:mm)%jEndOr = blocks(i)%jEndOr
          flowDoms(nn,1,1:mm)%kEndOr = blocks(i)%kEndOr

          ! Determine whether or not the block is moving.
          ! First initialize it to gridMotionSpecified. This is
          ! .true. if a rigid body motion was specified for the
          ! entire grid; otherwise it is .false.

          flowDoms(nn,1,1:mm)%blockIsMoving = gridMotionSpecified

          ! Check whether the corresponding cgns block is moving.
          ! Although it is possible that boundaries of a block rotate
          ! differently than the block itself, this should not be
          ! taken into account here; that's a matter of BC's.
          ! Here only the internal block structure is looked at.

          k = flowDoms(nn,1,1)%cgnsBlockID
          if( cgnsDoms(k)%rotatingFrameSpecified ) &
               flowDoms(nn,1,1:mm)%blockIsMoving = .true.

          ! For an unsteady computation on a deforming mesh
          ! blockIsMoving is always .true. Note that the time spectral
          ! method is also an unsteady computation.

          if(deforming_Grid .and.            &
               (equationMode == unsteady .or. &
               equationMode == timeSpectral)) &
               flowDoms(nn,1,1:mm)%blockIsMoving = .true.

          ! Set addGridVelocities to blockIsMoving. This could be
          ! overwritten later when the code is running in python mode.

          flowDoms(nn,1,1:mm)%addGridVelocities = &
               flowDoms(nn,1,1:mm)%blockIsMoving

          ! Set the number of subfaces and allocate the memory for the
          ! subface info. Note that this memory is only allocated for
          ! the first spectral time value; the other ones are identical.

          flowDoms(nn,1,1:mm)%nBocos   = blocks(i)%nBocos
          flowDoms(nn,1,1:mm)%n1to1    = blocks(i)%n1to1
          flowDoms(nn,1,1:mm)%nSubface = blocks(i)%nSubface
          j                           = blocks(i)%nSubface

          allocate(flowDoms(nn,1,1)%BCType(j),      &
               flowDoms(nn,1,1)%BCFaceID(j),    &
               flowDoms(nn,1,1)%cgnsSubface(j), &
               flowDoms(nn,1,1)%inBeg(j),       &
               flowDoms(nn,1,1)%jnBeg(j),       &
               flowDoms(nn,1,1)%knBeg(j),       &
               flowDoms(nn,1,1)%inEnd(j),       &
               flowDoms(nn,1,1)%jnEnd(j),       &
               flowDoms(nn,1,1)%knEnd(j),       &
               flowDoms(nn,1,1)%dinBeg(j),      &
               flowDoms(nn,1,1)%djnBeg(j),      &
               flowDoms(nn,1,1)%dknBeg(j),      &
               flowDoms(nn,1,1)%dinEnd(j),      &
               flowDoms(nn,1,1)%djnEnd(j),      &
               flowDoms(nn,1,1)%dknEnd(j),      &
               flowDoms(nn,1,1)%neighProc(j),   &
               flowDoms(nn,1,1)%neighBlock(j),  &
               flowDoms(nn,1,1)%l1(j),          &
               flowDoms(nn,1,1)%l2(j),          &
               flowDoms(nn,1,1)%l3(j),          &
               flowDoms(nn,1,1)%groupNum(j),    &
               stat=ierr)


          if(ierr /= 0)                   &
               call terminate("loadBalance", &
               "Memory allocation failure for subface info")

          ! Determine the new numbering of the boundary subfaces, such
          ! that the viscous subfaces are numbered first, followed by
          ! the inViscid subfaces, etc.

          allocate(oldSubfaceID(blocks(i)%nBocos), stat=ierr)
          if(ierr /= 0)                   &
               call terminate("loadBalance", &
               "Memory allocation failure for oldSubfaceID")

          call sortSubfaces(oldSubfaceID, blocks(i))

          ! Initialize the number of viscous boundary subfaces to 0.

          nViscBocos = 0

          ! Copy the info. Set the neighboring proc and block id to -1
          ! and 0 repectively in this loop. This is okay for boundary
          ! faces, but must be corrected for the internal block
          ! boundaries.

          do j=1,blocks(i)%nSubface

             ! Store the old subface id in k. For boundary faces the
             ! sorting is taken into account; for 1 to 1 subfaces the
             ! number is identical to the subface id in block.

             k = j
             if(j <= blocks(i)%nBocos) k = oldSubfaceID(j)

             ! Copy the info.

             flowDoms(nn,1,1)%BCType(j)      = blocks(i)%BCType(k)
             flowDoms(nn,1,1)%BCFaceID(j)    = blocks(i)%BCFaceID(k)
             flowDoms(nn,1,1)%cgnsSubface(j) = blocks(i)%cgnsSubface(k)

             flowDoms(nn,1,1)%inBeg(j) = blocks(i)%inBeg(k)
             flowDoms(nn,1,1)%jnBeg(j) = blocks(i)%jnBeg(k)
             flowDoms(nn,1,1)%knBeg(j) = blocks(i)%knBeg(k)
             flowDoms(nn,1,1)%inEnd(j) = blocks(i)%inEnd(k)
             flowDoms(nn,1,1)%jnEnd(j) = blocks(i)%jnEnd(k)
             flowDoms(nn,1,1)%knEnd(j) = blocks(i)%knEnd(k)

             flowDoms(nn,1,1)%dinBeg(j) = blocks(i)%dinBeg(k)
             flowDoms(nn,1,1)%djnBeg(j) = blocks(i)%djnBeg(k)
             flowDoms(nn,1,1)%dknBeg(j) = blocks(i)%dknBeg(k)
             flowDoms(nn,1,1)%dinEnd(j) = blocks(i)%dinEnd(k)
             flowDoms(nn,1,1)%djnEnd(j) = blocks(i)%djnEnd(k)
             flowDoms(nn,1,1)%dknEnd(j) = blocks(i)%dknEnd(k)

             flowDoms(nn,1,1)%neighProc(j)  = -1
             flowDoms(nn,1,1)%neighBlock(j) =  0

             flowDoms(nn,1,1)%l1(j) = blocks(i)%l1(k)
             flowDoms(nn,1,1)%l2(j) = blocks(i)%l2(k)
             flowDoms(nn,1,1)%l3(j) = blocks(i)%l3(k)

             flowDoms(nn,1,1)%groupNum(j) = blocks(i)%groupNum(k)

             ! Update the number of viscous boundaries if this
             ! is a viscous subface.

             if(flowDoms(nn,1,1)%BCType(j) == NSWallAdiabatic .or. &
                  flowDoms(nn,1,1)%BCType(j) == NSWallIsothermal)    &
                  nViscBocos = nViscBocos + 1
          enddo

          flowDoms(nn,1,1:mm)%nViscBocos = nViscBocos

          ! Correct the neighboring block and proc ID for internal
          ! block boundaries.

          do k=1,blocks(i)%n1to1
             j = blocks(i)%nBocos + k

             flowDoms(nn,1,1)%neighProc(j)  = part(blocks(i)%neighBlock(j))
             flowDoms(nn,1,1)%neighBlock(j) = &
                  blocks(blocks(i)%neighBlock(j))%blockID
          enddo

          ! Release the memory of oldSubfaceID.

          deallocate(oldSubfaceID, stat=ierr)
          if(ierr /= 0)                   &
               call terminate("loadBalance", &
               "Deallocation error for oldSubfaceID")

       endif myBlock

       ! Release the memory of the subface on this block.

       deallocate(blocks(i)%bcType,      blocks(i)%bcFaceid,    &
            blocks(i)%cgnsSubface, blocks(i)%inBeg,       &
            blocks(i)%jnBeg,       blocks(i)%knBeg,       &
            blocks(i)%inEnd,       blocks(i)%jnEnd,       &
            blocks(i)%knEnd,       blocks(i)%dinBeg,      &
            blocks(i)%djnBeg,      blocks(i)%dknBeg,      &
            blocks(i)%dinEnd,      blocks(i)%djnEnd,      &
            blocks(i)%dknEnd,      blocks(i)%neighBlock,  &
            blocks(i)%l1,          blocks(i)%l2,          &
            blocks(i)%l3,          blocks(i)%groupNum,    &
            stat=ierr)
       if(ierr /= 0)                   &
            call terminate("loadBalance", &
            "Deallocation error for boundary info")
    enddo domains

    !       Determine the number of processors, the processor ID's on
    !       which the original cgns blocks are stored, the local
    !       block ID's and the nodal ranges of the subblocks. As blocks
    !       can be split during run-time, multiple processors can store a
    !       part of original block.
    !
    ! Allocate the memory for subblocksOfCGNS.

    allocate(subblocksOfCGNS(nBlocks), stat=ierr)
    if(ierr /= 0)                   &
         call terminate("loadBalance", &
         "Memory allocation failure for subblocksOfCGNS")

    ! Copy the data into subblocksOfCGNS.

    do nn=1,nBlocks
       subblocksOfCGNS(nn)%cgnsBlockID = blocks(nn)%cgnsBlockID
       subblocksOfCGNS(nn)%procID      = part(nn)
       subblocksOfCGNS(nn)%blockID     = blocks(nn)%blockID

       subblocksOfCGNS(nn)%iBegOr = blocks(nn)%iBegOr
       subblocksOfCGNS(nn)%iEndOr = blocks(nn)%iEndOr
       subblocksOfCGNS(nn)%jBegOr = blocks(nn)%jBegOr
       subblocksOfCGNS(nn)%jEndOr = blocks(nn)%jEndOr
       subblocksOfCGNS(nn)%kBegOr = blocks(nn)%kBegOr
       subblocksOfCGNS(nn)%kEndOr = blocks(nn)%kEndOr
    enddo

    ! Sort subblocksOfCGNS in increasing order.

    call qsortSubblocksOfCGNSType(subblocksOfCGNS, nBlocks)

    ! Loop over the number of cgns blocks and find out the number of
    ! subblocks it contains.

    ii = 1
    subBlockLoop: do nn=1,cgnsNDom

       ! Determine the ending index jj in subblocksOfCGNS for this
       ! CGNS block. The starting index is ii.

       if(nn == cgnsNDom) then
          jj = nBlocks
       else
          jj = ii
          do
             if(subblocksOfCGNS(jj+1)%cgnsBlockID > nn) exit
             jj = jj + 1
          enddo
       endif

       ! Set nSubBlocks and allocate the memory for procStored,
       ! localBlockID, iBegOr, iEndOr, etc.

       cgnsDoms(nn)%nSubBlocks = jj - ii + 1
       k = cgnsDoms(nn)%nSubBlocks
       allocate(cgnsDoms(nn)%iBegOr(k),       cgnsDoms(nn)%iEndOr(k), &
            cgnsDoms(nn)%jBegOr(k),       cgnsDoms(nn)%jEndOr(k), &
            cgnsDoms(nn)%kBegOr(k),       cgnsDoms(nn)%kEndOr(k), &
            cgnsDoms(nn)%procStored(k),                           &
            cgnsDoms(nn)%localBlockID(k), stat=ierr)
       if(ierr /= 0)                   &
            call terminate("loadBalance", &
            "Memory allocation failure for procStored, &
            &localBlockID, iBegOr, iEndOr, etc.")

       ! Copy the processor ID's, the local block ID's
       ! and the subranges.

       do i=1,cgnsDoms(nn)%nSubBlocks
          j = i + ii - 1
          cgnsDoms(nn)%procStored(i)   = subblocksOfCGNS(j)%procID
          cgnsDoms(nn)%localBlockID(i) = subblocksOfCGNS(j)%blockID

          cgnsDoms(nn)%iBegOr(i) = subblocksOfCGNS(j)%iBegOr
          cgnsDoms(nn)%iEndOr(i) = subblocksOfCGNS(j)%iEndOr
          cgnsDoms(nn)%jBegOr(i) = subblocksOfCGNS(j)%jBegOr
          cgnsDoms(nn)%jEndOr(i) = subblocksOfCGNS(j)%jEndOr
          cgnsDoms(nn)%kBegOr(i) = subblocksOfCGNS(j)%kBegOr
          cgnsDoms(nn)%kEndOr(i) = subblocksOfCGNS(j)%kEndOr
       enddo

       ! Set ii for the next CGNS block.

       ii = jj + 1

    enddo subBlockLoop

    ! Release the memory of blocks, part and subblocksOfCGNS.

    deallocate(blocks, part, subblocksOfCGNS, stat=ierr)
    if(ierr /= 0)                   &
         call terminate("loadBalance", &
         "Deallocation error for blocks, part and &
         &subblocksOfCGNS")

    !j = 20+myID
    !do nn=1,ndom
    !  write(j,"(8I4)") nn, flowDoms(nn,1,1)%cgnsBlockID, &
    !                   flowDoms(nn,1,1)%iBegOr, flowDoms(nn,1,1)%iEndOr, &
    !                   flowDoms(nn,1,1)%jBegOr, flowDoms(nn,1,1)%jEndOr, &
    !                   flowDoms(nn,1,1)%kBegOr, flowDoms(nn,1,1)%kEndOr
    !enddo

  end subroutine loadBalanceGrid

  subroutine blockDistribution
    !
    !       blockDistribution determines the distribution of the blocks
    !       over the processors. If blocks must be split to obtain a good
    !       load balance an iterative algorithm is used to determine the
    !       best way to split them.
    !
    use constants
    use cgnsGrid, only : cgnsDoms, cgnsNDom
    use communication, only : myID, adflow_comm_world, nProc
    use inputParallel, only : loadBalanceIter, loadImbalance, splitBlocks
    use partitionMod, onlY : splitCGNSType, distributionBlockType, ubVec, &
         blocks, part, nBlocks, sortRangesSplitInfo
    use utils, only : terminate
    implicit none
    !
    !      Local variables.
    !
    integer :: ierr

    integer(kind=intType) :: nn, nx, ny, nz
    integer(kind=intType) :: nCellsTot, nCellsEven, nCellsUpper
    integer(kind=intType) :: nCellsPerProcMax
    integer(kind=intType) :: iter, iterMax

    logical :: cellsBalanced,   facesBalanced
    logical :: emptyPartitions, commNeglected

    type(splitCGNSType), dimension(cgnsNDom) :: splitInfo

    ! If it is not allowed to split the blocks, check that the
    ! number of blocks is equal or larger than the number of
    ! processors. If not, print an error message and exit.

    if(.not. splitBlocks .and. nProc > cgnsNDom) then

       ! Processor 0 prints the error message, while the others wait
       ! to get killed.

       if(myID == 0)                         &
            call terminate("blockDistribution", &
            "Number of processors is larger than number &
            &of blocks, but it is not allowed to split &
            &blocks")
       call mpi_barrier(ADflow_comm_world, ierr)
    endif
    !
    !       Determine how the blocks must be split (if allowed) for
    !       load balancing reasons.
    !
    ! Initialize splitInfo to the original cgns blocks and determine
    ! the total number of cells.

    nCellsTot = 0
    nBlocks   = cgnsNDom

    do nn=1,cgnsNDom
       splitInfo(nn)%nSubBlocks = 1

       allocate(splitInfo(nn)%ranges(1,3,2), stat=ierr)
       if(ierr /= 0)                         &
            call terminate("blockDistribution", &
            "Memory allocation failure for ranges")

       splitInfo(nn)%ranges(1,1,1) = 1
       splitInfo(nn)%ranges(1,2,1) = 1
       splitInfo(nn)%ranges(1,3,1) = 1

       splitInfo(nn)%ranges(1,1,2) = cgnsDoms(nn)%il
       splitInfo(nn)%ranges(1,2,2) = cgnsDoms(nn)%jl
       splitInfo(nn)%ranges(1,3,2) = cgnsDoms(nn)%kl

       nx = cgnsDoms(nn)%nx
       ny = cgnsDoms(nn)%ny
       nz = cgnsDoms(nn)%nz

       nCellsTot = nCellsTot + nx*ny*nz
    enddo

    ! Determine the desirable number of cells per processor and
    ! determine the upper bound based on the imbalance tolerance.
    ! Initialize nCellsPerProcMax, which is used to decide whether
    ! or not a block should be split, to this upper limit.

    nCellsEven       = nCellsTot/nProc
    nCellsUpper      = nCellsEven + loadImbalance*nCellsEven
    nCellsPerProcMax = nCellsUpper

    ! Start the loop to determine the block topology to be used
    ! in the computation.

    initSplitLoop: do iter=1,2

       ! Loop over the number of original cgns blocks and determine
       ! whether or not they must be split further.

       do nn=1,cgnsNDom
          if(.not. splittingIsOkay(nn) ) &
               call splitBlockInitialization(nn)
       enddo

       ! Exit the loop if the number of processors is smaller than or
       ! equal to the number of blocks.

       if(nBlocks >= nProc) exit

       ! The number of blocks is smaller than the number of processors.
       ! Set the tolerance for splitting to the desired number of
       ! cells on a processor.

       nCellsPerProcMax  = nCellsEven

    enddo initSplitLoop

    ! Set the number of iterations to determine the distribution over
    ! the processors. If blocks cannot be split this is set to 1;
    ! otherwise it is set to 2.

    iterMax = 1
    if( splitBlocks ) iterMax = loadBalanceIter

    ! Loop to determine a good load balance.

    distributionLoop: do iter=1,iterMax

       ! Determine the computational blocks from the splitting info of
       ! the original blocks.

       call determineComputeBlocks(splitInfo)
       ! Apply the graph partitioning to the computational blocks.

       call graphPartitioning(emptyPartitions, commNeglected)

       ! Determine whether the load balance is okay. If empty
       ! partitions are present the load balance is per definition
       ! not okay and there is no need to call checkLoadBalance.

       if( emptyPartitions ) then
          cellsBalanced = .false.
          facesBalanced = .false.
       else
          call checkLoadBalance(cellsBalanced, facesBalanced)
       endif

       ! Exit the loop if the cells or the faces are load balanced
       ! or if the maximum number of iterations have been reached.

       if(cellsBalanced .or. facesBalanced .or. &
            iter == iterMax) exit

       ! exit

       ! Split some blocks on the processors with too many cells/faces.
       call splitBlocksLoadBalance

    enddo distributionLoop

    ! Deallocate the memory for ranges.

    do nn=1,cgnsNDom
       deallocate(splitInfo(nn)%ranges, stat=ierr)
       if(ierr /= 0)                         &
            call terminate("blockDistribution", &
            "Deallocation failure for ranges")
    enddo

    ! If empty partitions are present print an error message and
    ! exit. Only processor 0 prints the message to avoid a mess.

    if( emptyPartitions ) then
       if(myID == 0)                         &
            call terminate("blockDistribution", &
            "Empty partitions present")
       call mpi_barrier(ADflow_comm_world, ierr)
    endif

    ! Processor 0 prints a warning if the communication was
    ! neglected to obtain a valid partitioning.

    if(commNeglected .and. myID == 0) then
       print "(a)", "#"
       print "(a)", "#                    Warning"
       print "(a)", "# Communication costs neglected to obtain&
            & a valid partitioning."
    endif

    ! If the load imbalance tolerance was not met, print a warning
    ! message if I am processor 0.

    if (myid == 0) then
       if(.not.(cellsBalanced .and. facesBalanced)) then
          print "(a)", "#"
          print "(a)", "#                    Warning"
          print 100, loadImbalance
          print 101, ubvec(1), ubvec(2)
          print "(a)", "#"
100       format("# Specified load imbalance tolerance",1X,F6.3,1X,"not &
               &achieved.")
101       format("# I continue with",1X,F6.3,1X,"load imbalance for the &
               &cells and",1X,F6.3,1X,"for the faces")
       else
          print "(a)","#"

          print 102, loadImbalance
          print 103, ubvec(1), ubvec(2)
102       format("# Specified load imbalance tolerance",1X,F6.3,1X,"acheived")
103       format("# Continuing with",1X,F6.3,1X,"load imbalance for the &
               &cells and",1X,F6.3,1X,"for the faces")
          print "(a)", "#"
       end if

    endif

    !=================================================================

  contains

    !===============================================================

    logical function splittingIsOkay(cgnsID)
      !
      !         splittingIsOkay determines whether or not the splitting of
      !         the given cgns block is okay in the sense that all subblocks
      !         are smaller than the allowed number of cells and faces.
      !
      implicit none
      !
      !        Function argument
      !
      integer(kind=intType), intent(in) :: cgnsID
      !
      !        Local variables.
      !
      integer(kind=intType) :: i, nx, ny, nz, nCells

      ! Initialize splittingIsOkay to .true.

      splittingIsOkay = .true.

      ! Loop over the subblocks.

      do i=1,splitInfo(cgnsID)%nSubBlocks

         ! Determine the number of cells in the three directions.

         nx = splitInfo(cgnsID)%ranges(i,1,2) &
              - splitInfo(cgnsID)%ranges(i,1,1)
         ny = splitInfo(cgnsID)%ranges(i,2,2) &
              - splitInfo(cgnsID)%ranges(i,2,1)
         nz = splitInfo(cgnsID)%ranges(i,3,2) &
              - splitInfo(cgnsID)%ranges(i,3,1)

         ! Determine the number of cells for this subblock.

         nCells = nx*ny*nz

         ! Check whether this number is smaller or equal to the
         ! maximum allowed number. If not, set splittingIsOkay to
         ! .false. and exit the loop.

         if(nCells > nCellsPerProcMax) then
            splittingIsOkay = .false.
            exit
         endif

      enddo

    end function splittingIsOkay

    !===============================================================

    subroutine splitBlockInitialization(cgnsID)
      !
      !         splitBlockInitialization splits the given cgns block ID
      !         into a number of subbocks during the initialization phase.
      !
      implicit none
      !
      !        Subroutine argument.
      !
      integer(kind=intType), intent(in) :: cgnsID
      !
      !        Local variables.
      !
      integer :: ierr

      integer(kind=intType) :: nn, mm
      integer(kind=intType) :: nx, ny, nz, nCells, nSub

      integer(kind=intType), dimension(:,:,:), allocatable :: tmpRange

      type(distributionBlockType) :: tmpBlock

      ! Check whether it is allowed to split blocks.

      if(.not. splitBlocks)                        &
           call terminate("splitBlockInitialization", &
           "Block must be split for load balance, &
           &but I am not allowed to do so")

      ! Store the number of cells of this block a bit easier.

      nx     = cgnsDoms(cgnsID)%nx
      ny     = cgnsDoms(cgnsID)%ny
      nz     = cgnsDoms(cgnsID)%nz
      nCells = nx*ny*nz

      ! Copy the information from the given cgns block into
      ! tmpBlock, such that the routine splitBlock can be used.

      ! First the scalar info.

      tmpBlock%nx = nx; tmpBlock%il = nx + 1
      tmpBlock%ny = ny; tmpBlock%jl = ny + 1
      tmpBlock%nz = nz; tmpBlock%kl = nz + 1

      tmpBlock%nCell = nCells
      tmpBlock%nface = (nx+1)*ny*nz + (ny+1)*nx*nz &
           + (nz+1)*nx*ny

      tmpBlock%cgnsBlockID = cgnsID

      tmpBlock%iBegor = 1; tmpBlock%iEndor = tmpBlock%il
      tmpBlock%jBegor = 1; tmpBlock%jEndor = tmpBlock%jl
      tmpBlock%kBegor = 1; tmpBlock%kEndor = tmpBlock%kl

      tmpBlock%nBocos   = cgnsDoms(cgnsID)%nBocos
      tmpBlock%n1to1    = cgnsDoms(cgnsID)%n1to1
      tmpBlock%nSubface = tmpBlock%nBocos + tmpBlock%n1to1

      ! Allocate the memory for the nodal ranges of the subfaces
      ! and nullify the other pointers.

      nSub = tmpBlock%nSubface

      allocate(tmpBlock%BCType(nSub), tmpBlock%BCFaceID(nSub), &
           tmpBlock%inBeg(nSub),  tmpBlock%inEnd(nSub),    &
           tmpBlock%jnBeg(nSub),  tmpBlock%jnEnd(nSub),    &
           tmpBlock%knBeg(nSub),  tmpBlock%knEnd(nSub),    &
           stat=ierr)
      if(ierr /= 0) &
           call terminate("splitBlockInitialization", &
           "Deallocation failure for the subface &
           &info in tmpBlock")

      nullify(tmpBlock%cgnsSubface, tmpBlock%neighBlock, &
           tmpBlock%dinBeg,      tmpBlock%dinEnd,     &
           tmpBlock%djnBeg,      tmpBlock%djnEnd,     &
           tmpBlock%dknBeg,      tmpBlock%dknEnd,     &
           tmpBlock%l1,          tmpBlock%l2,         &
           tmpBlock%l3,          tmpBlock%groupNum)

      ! Copy the range of the subfaces into tmpBlock and set the
      ! corresponding boundary condition.
      ! First the boundary faces.

      do nn=1,cgnsDoms(cgnsID)%nBocos
         tmpBlock%inBeg(nn) = cgnsDoms(cgnsID)%bocoInfo(nn)%iBeg
         tmpBlock%jnBeg(nn) = cgnsDoms(cgnsID)%bocoInfo(nn)%jBeg
         tmpBlock%knBeg(nn) = cgnsDoms(cgnsID)%bocoInfo(nn)%kBeg

         tmpBlock%inEnd(nn) = cgnsDoms(cgnsID)%bocoInfo(nn)%iEnd
         tmpBlock%jnEnd(nn) = cgnsDoms(cgnsID)%bocoInfo(nn)%jEnd
         tmpBlock%knEnd(nn) = cgnsDoms(cgnsID)%bocoInfo(nn)%kEnd

         tmpBlock%BCType(nn) = cgnsDoms(cgnsID)%bocoInfo(nn)%BCType
      enddo

      ! And the internal block faces; set the boundary condition to
      ! B2BMatch.

      do mm=1,cgnsDoms(cgnsID)%n1to1
         nn = mm + cgnsDoms(cgnsID)%nBocos

         tmpBlock%inBeg(nn) = cgnsDoms(cgnsID)%conn1to1(mm)%iBeg
         tmpBlock%jnBeg(nn) = cgnsDoms(cgnsID)%conn1to1(mm)%jBeg
         tmpBlock%knBeg(nn) = cgnsDoms(cgnsID)%conn1to1(mm)%kBeg

         tmpBlock%inEnd(nn) = cgnsDoms(cgnsID)%conn1to1(mm)%iEnd
         tmpBlock%jnEnd(nn) = cgnsDoms(cgnsID)%conn1to1(mm)%jEnd
         tmpBlock%knEnd(nn) = cgnsDoms(cgnsID)%conn1to1(mm)%kEnd

         tmpBlock%BCType(nn) = B2BMatch
      enddo

      ! Determine for the block face on which the subface is located.
      ! Assume that the subface connectivity is correct. This will be
      ! tested later on in determineComputeBlocks.

      do nn=1,tmpBlock%nSubface
         if(tmpBlock%inBeg(nn) == tmpBlock%inEnd(nn)) then
            tmpBlock%BCFaceID(nn) = iMax
            if(tmpBlock%inBeg(nn) == 1) tmpBlock%BCFaceID(nn) = iMin
         else if(tmpBlock%jnBeg(nn) == tmpBlock%jnEnd(nn)) then
            tmpBlock%BCFaceID(nn) = jMax
            if(tmpBlock%jnBeg(nn) == 1) tmpBlock%BCFaceID(nn) = jMin
         else
            tmpBlock%BCFaceID(nn) = kMax
            if(tmpBlock%knBeg(nn) == 1) tmpBlock%BCFaceID(nn) = kMin
         endif
      enddo

      ! Determine the number of subblocks into which this block is
      ! to be split.

      nSub = nint(real(nCells,realType)/real(nCellsEven,realType))
      nSub = max(nSub,2_intType)

      ! Deallocate the memory of the splitting info and allocate
      ! the memory of tmpRange.

      deallocate(splitInfo(cgnsID)%ranges, stat=ierr)
      if(ierr /= 0) &
           call terminate("splitBlockInitialization", &
           "Deallocation failure for ranges")

      allocate(tmpRange(nSub,3,2), stat=ierr)
      if(ierr /= 0) &
           call terminate("splitBlockInitialization", &
           "Memory allocation failure for tmpRange")

      ! Split the block into the subblocks. It is possible that
      ! the block is split into less subblocks than the desired
      ! number. The actual number of subblocks is returned in nSub.

      call splitBlock(tmpBlock, nSub, nCellsEven, tmpRange)

      ! Allocate the memory for ranges.

      allocate(splitInfo(cgnsID)%ranges(nSub,3,2), stat=ierr)
      if(ierr /= 0) &
           call terminate("splitBlockInitialization", &
           "Memory allocation failure for ranges")

      ! Determine the new number of computational blocks. This is
      ! the old number plus the number or extra blocks created
      ! by the splitting.

      nBlocks = nBlocks + nSub - splitInfo(cgnsID)%nSubBlocks

      ! Copy tmpRange into ranges of splitInfo.

      splitInfo(cgnsID)%nSubBlocks = nSub

      do nn=1,nSub
         splitInfo(cgnsID)%ranges(nn,1,1) = tmpRange(nn,1,1)
         splitInfo(cgnsID)%ranges(nn,2,1) = tmpRange(nn,2,1)
         splitInfo(cgnsID)%ranges(nn,3,1) = tmpRange(nn,3,1)

         splitInfo(cgnsID)%ranges(nn,1,2) = tmpRange(nn,1,2)
         splitInfo(cgnsID)%ranges(nn,2,2) = tmpRange(nn,2,2)
         splitInfo(cgnsID)%ranges(nn,3,2) = tmpRange(nn,3,2)
      enddo

      ! Deallocate the memory allocated for tmpBlock and tmpRange.

      deallocate(tmpBlock%BCType, tmpBlock%BCFaceID, &
           tmpBlock%inBeg,  tmpBlock%inEnd,    &
           tmpBlock%jnBeg,  tmpBlock%jnEnd,    &
           tmpBlock%knBeg,  tmpBlock%knEnd,    &
           tmpRange, stat=ierr)
      if(ierr /= 0) &
           call terminate("splitBlockInitialization", &
           "Deallocation failure for tmpBlock &
           &and tmpRange")

      ! Sort the subranges of this block in increasing order.

      call sortRangesSplitInfo(splitInfo(cgnsID))

    end subroutine splitBlockInitialization

    !===============================================================

    subroutine splitBlocksLoadBalance
      !
      !         splitBlocksLoadBalance splits some (sub)blocks even
      !         further to obtain a better load balance.
      !
      use sorting, only : bsearchIntegers, qsortIntegers
      implicit none
      !
      !        Local variables.
      !
      integer :: ierr

      integer(kind=intType) :: i, j, k, kk, jj
      integer(kind=intType) :: nCellDiff, proc, nSub, cgnsID
      integer(kind=intType) :: ncgnsSort

      integer(kind=intType), dimension(0:cgnsNDom) :: nSubBlocks
      integer(kind=intType), dimension(cgnsNDom)   :: cgnsSort

      integer(kind=intType), dimension(nProc)   :: nCell, nCellOr
      integer(kind=intType), dimension(nProc)   :: tmp, procNCell
      integer(kind=intType), dimension(0:nProc) :: nBlockProc
      integer(kind=intType), dimension(0:nProc) :: multNCell
      integer(kind=intType), dimension(nBlocks) :: blockProc

      integer(kind=intType), dimension(:,:,:), allocatable :: tmpRange
      integer(kind=intType), dimension(:,:,:), pointer :: oldRanges

      logical, dimension(cgnsNDom) :: cgnsBlockFlagged

      ! Determine the number of blocks per cgns block in cumulative
      ! storage format. These values will serve as an offset to
      ! determine the local subblock ID.
      ! Initialize in the same loop cgnsBlockFlagged to .false.

      nSubBlocks(0) = 0
      do i=1,cgnsNDom
         nSubBlocks(i) = nSubBlocks(i-1) + splitInfo(i)%nSubBlocks
         cgnsBlockFlagged(i) = .false.
      enddo

      ! Determine the number of cells and blocks per processor.
      ! Note that part(i) contains the processor ID, which start at 0.

      nCell      = 0
      nBlockProc = 0

      do i=1,nBlocks
         j = part(i) + 1
         nCell(j)      = nCell(j) + blocks(i)%nCell
         nBlockProc(j) = nBlockProc(j) + 1
      enddo

      ! Put nBlockProc in cumulative storage format and store the
      ! starting entries in tmp, which will be used as a counter to
      ! put blockProc in the correct place.

      do i=1,nProc
         tmp(i)        = nBlockProc(i-1)
         nBlockProc(i) = nBlockProc(i) + tmp(i)
      enddo

      ! Determine the block ID's for every processor. Again note
      ! that part(i) starts at 0; tmp is used as a counter variable.

      do i=1,nBlocks
         j                 = part(i) + 1
         tmp(j)            = tmp(j) + 1
         blockProc(tmp(j)) = i
      enddo

      ! Sort nCell in increasing order. Store the original values.

      nCellOr   = nCell
      nCellDiff = nProc
      call qsortIntegers(nCell, nCellDiff)

      ! Determine the different number of cells per proc and store
      ! the multiplicity.

      multNCell(0) = 0
      multNCell(1) = 1
      nCellDiff    = 1

      do i=2,nProc
         if(nCell(i) == nCell(nCellDiff)) then
            multNCell(nCellDiff) = multNCell(nCellDiff) + 1
         else
            nCellDiff = nCellDiff + 1
            multNCell(nCellDiff) = multNCell(nCellDiff-1) + 1
            nCell(nCellDiff)     = nCell(i)
         endif
      enddo

      ! Set tmp to multNCell; it will serve as a counter to determine
      ! the corresponding processor ID's of the sorted cell numbers

      do i=1,nCellDiff
         tmp(i) = multNCell(i-1)
      enddo

      ! Determine the processor ID's corresponding to the sorted
      ! number of cells. Note that the processor ID's start at 0.

      do i=1,nProc
         j = bsearchIntegers(nCellOr(i), nCell(1:nCellDiff))
         tmp(j) = tmp(j) + 1
         procNCell(tmp(j)) = i - 1
      enddo

      ! Loop to subdivide blocks. As nCell is sorted in increasing
      ! order, the largest number of cells are at the end of this
      ! array. Therefore this loop starts at the back.

      i = nCellDiff
      ncgnsSort = 0
      divisionLoop: do

         ! Condition to exit the loop. The number of cells per
         ! processor is allowed in the load balance.

         if(nCell(i) <= nCellsUpper) exit

         ! Loop over the processors with this number of cells.

         processorLoop: do j=(multNCell(i-1)+1), multNCell(i)

            ! Store the processor ID a bit easier.

            proc = procNCell(j)

            ! Determine the largest block of this processor, which will
            ! be stored in index jj. As only the largest is needed a
            ! linear search algorithm is okay.

            k  = nBlockProc(proc) + 1
            jj = blockProc(k)

            do k=(nBlockProc(proc) + 2),nBlockProc(proc+1)
               kk = blockProc(k)
               if(blocks(kk)%nCell > blocks(jj)%nCell) jj = kk
            enddo

            ! Determine the local subblock number of computational block
            ! jj in the original cgns block. This will be stored in kk.

            cgnsID = blocks(jj)%cgnsBlockID
            kk     = jj - nSubBlocks(cgnsID-1)

            ! Determine the number of cells, which should be split from
            ! block jj, such that the load balance of processor proc
            ! would be optimal (in terms of cells).

            nCellDiff = nCell(i) - nCellsEven

            ! A very theoretical case would be that nCellDiff is larger
            ! than the number of cells of block jj. This could occur,
            ! because the partitioning is a multi-constraint
            ! partitioning; however, it is not likely. In this case
            ! nothing is done and the next processor is treated.

            if(nCellDiff >= blocks(jj)%nCell) cycle

            ! Determine the situation we are having here. If nCellDiff
            ! is less or equal to the number of cells to obtain a good
            ! load balance for the least loaded processor the block will
            ! be split into 2. Test this.

            if(nCellDiff <= (nCellsUpper - nCell(1))) then

               nSub = 2

            else

               ! nCellDiff is larger than the number of cells to create
               ! a good load balance for the least loaded processor.
               ! This probably means that the block should be split in
               ! more than 2 subblocks.
               ! Store the number of cells which should be kept in block
               ! jj for an optimal load balance in nCellDiff.

               nCellDiff = blocks(jj)%nCell - nCellDiff

               ! Determine the number of subblocks.

               nSub = nint(real(blocks(jj)%nCell,realType) &
                    /      real(nCellDiff,realType))
               nSub = max(nSub,2_intType)

            endif

            ! Allocate the memory for tmpRange, which will store the
            ! block splitting of block jj.

            allocate(tmpRange(nSub,3,2), stat=ierr)
            if(ierr /= 0)                              &
                 call terminate("splitBlocksLoadBalance", &
                 "Memory allocation failure for tmpRange")

            ! Split computational block into the desired number of
            ! subblocks, if possible. The actual number of splittings
            ! is returned in nSub.

            call splitBlock(blocks(jj), nSub, nCellDiff, tmpRange)

            ! Store the old ranges, store the old number of subblocks
            ! in jj, determine the new number of subblocks and allocate
            ! the memory for ranges.

            jj         = splitInfo(cgnsID)%nSubBlocks
            oldRanges => splitInfo(cgnsID)%ranges

            splitInfo(cgnsID)%nSubBlocks = jj + nSub - 1

            k = splitInfo(cgnsID)%nSubBlocks
            allocate(splitInfo(cgnsID)%ranges(k,3,2), stat=ierr)
            if(ierr /= 0)                              &
                 call terminate("splitBlocksLoadBalance", &
                 "Memory allocation failure for ranges")

            ! Determine the new number of computational blocks.
            ! This is the old number plus the number or extra blocks
            ! created by the splitting.

            nBlocks = nBlocks + splitInfo(cgnsID)%nSubBlocks - jj

            ! Copy the original info back into ranges.

            do k=1,jj
               splitInfo(cgnsID)%ranges(k,1,1) = oldRanges(k,1,1)
               splitInfo(cgnsID)%ranges(k,2,1) = oldRanges(k,2,1)
               splitInfo(cgnsID)%ranges(k,3,1) = oldRanges(k,3,1)

               splitInfo(cgnsID)%ranges(k,1,2) = oldRanges(k,1,2)
               splitInfo(cgnsID)%ranges(k,2,2) = oldRanges(k,2,2)
               splitInfo(cgnsID)%ranges(k,3,2) = oldRanges(k,3,2)
            enddo

            ! The splitting kk has disappeared. Replace it by the
            ! first splitting in tmpRange.

            splitInfo(cgnsID)%ranges(kk,1,1) = tmpRange(1,1,1)
            splitInfo(cgnsID)%ranges(kk,2,1) = tmpRange(1,2,1)
            splitInfo(cgnsID)%ranges(kk,3,1) = tmpRange(1,3,1)

            splitInfo(cgnsID)%ranges(kk,1,2) = tmpRange(1,1,2)
            splitInfo(cgnsID)%ranges(kk,2,2) = tmpRange(1,2,2)
            splitInfo(cgnsID)%ranges(kk,3,2) = tmpRange(1,3,2)

            ! Add the other splittings to the end of ranges.

            do k=2,nSub
               jj = jj + 1

               splitInfo(cgnsID)%ranges(jj,1,1) = tmpRange(k,1,1)
               splitInfo(cgnsID)%ranges(jj,2,1) = tmpRange(k,2,1)
               splitInfo(cgnsID)%ranges(jj,3,1) = tmpRange(k,3,1)

               splitInfo(cgnsID)%ranges(jj,1,2) = tmpRange(k,1,2)
               splitInfo(cgnsID)%ranges(jj,2,2) = tmpRange(k,2,2)
               splitInfo(cgnsID)%ranges(jj,3,2) = tmpRange(k,3,2)
            enddo

            ! Deallocate the memory of oldRanges and tmpRange.

            deallocate(oldRanges, tmpRange, stat=ierr)
            if(ierr /= 0)                                 &
                 call terminate("splitBlocksLoadBalance", &
                 "Deallocation failure for oldRanges and &
                 &tmpRange")

            ! Store the cgns ID, if not already stored. In this way
            ! it is known for which cgns blocks the subranges must be
            ! sorted.

            if(.not. cgnsBlockFlagged(cgnsID)) then
               cgnsBlockFlagged(cgnsID) = .true.
               ncgnsSort = ncgnsSort + 1
               cgnsSort(ncgnsSort) = cgnsID
            endif

         enddo processorLoop

         ! Decrease the counter i for the next number of cells.

         i = i - 1

      enddo divisionLoop

      ! Sort the subranges of the CGNS blocks that were split in this
      ! routine in increasing order.

      do i=1,ncgnsSort
         call sortRangesSplitInfo(splitInfo(cgnsSort(i)))
      enddo

    end subroutine splitBlocksLoadBalance

  end subroutine blockDistribution

  subroutine initFlowDoms
    !
    !       initFlowDoms allocates the memory for flowDoms and initializes
    !       its pointers to null pointers, such that they do not have
    !       random targets.
    !
    use constants
    use block, only : flowDoms, nDom
    use inputIteration, only : nMGLevels, mgStartLevel
    use inputTimeSpectral, only : nTimeIntervalsSpectral
    use utils, only : terminate, nullifyFlowDomPointers
    implicit none
    !
    !      Local variables.
    !
    integer :: ierr

    integer(kind=intType) :: i, j, k, nn

    ! Allocate the memory for flowDoms. Set nn to the maximum of the
    ! number of mg levels needed in the cycle and mg start level.
    ! This is namely the amount of grid levels the solver needs.

    nn = max(nMGLevels, mgStartlevel)
    allocate(flowDoms(nDom, nn, nTimeIntervalsSpectral), stat=ierr)
    if(ierr /= 0)                    &
         call terminate("initFlowDoms", &
         "Memory allocation failure for flowDoms")

    ! Loop over all the blocks and initialize its pointers to the
    ! null-pointer.

    do k=1,nTimeIntervalsSpectral
       do j=1,nn
          do i=1,nDom
             call nullifyFlowDomPointers(i,j,k)
          enddo
       enddo
    enddo

  end subroutine initFlowDoms
  subroutine sortSubfaces(oldSubfaceID, blockID)
    !
    !       sortSubfaces sorts the boundary subfaces of the given block
    !       such that viscous subfaces are numbered first, followed by
    !       inviscid, etc.
    !
    use constants
    use partitionMod, only : distributionBlockType
    use sorting, only : qsortIntegers, bsearchIntegers
    implicit none
    !
    !      Subroutine arguments
    !
    integer(kind=intType), dimension(*), intent(out) :: oldSubfaceID
    type(distributionBlockType), intent(in)          :: blockID
    !
    !      Local variables.
    !
    integer(kind=intType) :: i, ii, nDiff

    integer(kind=intType), dimension(blockID%nBocos) :: bcPrior
    integer(kind=intType), dimension(blockID%nBocos) :: bcPriorSort

    integer(kind=intType), dimension(0:blockID%nBocos) :: mult

    ! Loop over the boundary subfaces and determine the priorities.
    ! Store the priorities in both bcPrior and bcPriorSort.

    do i=1,blockID%nBocos

       select case (blockID%BCType(i))

       case (NSWallAdiabatic)
          bcPrior(i) = 1

       case (NSWallIsothermal)
          bcPrior(i) = 2

       case (EulerWall)
          bcPrior(i) = 3

       case (Symm)
          bcPrior(i) = 4

       case (SymmPolar)
          bcPrior(i) = 5

       case (FarField)
          bcPrior(i) = 6

       case (SupersonicInflow)
          bcPrior(i) = 7

       case (SubsonicInflow)
          bcPrior(i) = 8

       case (SupersonicOutflow)
          bcPrior(i) = 9

       case (SubsonicOutflow)
          bcPrior(i) = 10

       case (MassBleedInflow)
          bcPrior(i) = 11

       case (MassBleedOutflow)
          bcPrior(i) = 12

       case (mDot)
          bcPrior(i) = 13

       case (bcThrust)
          bcPrior(i) = 14

       case (Extrap)
          bcPrior(i) = 15

       case (SlidingInterface)
          bcPrior(i) = 19

       case (OversetOuterBound)
          bcPrior(i) = 20

       case (DomainInterfaceAll)
          bcPrior(i) = 21

       case (DomainInterfaceRhoUVW)
          bcPrior(i) = 22

       case (DomainInterfaceP)
          bcPrior(i) = 23

       case (DomainInterfaceRho)
          bcPrior(i) = 24

       case (DomainInterfaceTotal)
          bcPrior(i) = 25

       end select

       bcPriorSort(i) = bcPrior(i)

    enddo

    ! Sort bcPriorSort in increasing order.

    call qsortIntegers(bcPriorSort, blockID%nBocos)

    ! Get rid of the multiple entries and store the multiplicity in
    ! cumulative storage format. nDiff contains the number of
    ! different boundary conditions for this block. The initialization
    ! with the min function is necessary to be able to treat blocks
    ! without a boundary condition correctly.

    nDiff       = min(1_intType, blockID%nBocos)
    mult(0)     = 0
    mult(nDiff) = 1

    do i=2,blockID%nBocos
       if(bcPriorSort(i) == bcPriorSort(nDiff)) then
          mult(nDiff) = mult(nDiff) + 1
       else
          nDiff              = nDiff + 1
          mult(nDiff)        = mult(nDiff-1) + 1
          bcPriorSort(nDiff) = bcPriorSort(i)
       endif
    enddo

    ! Determine the old subface ID by searching in the sorted
    ! priorities.

    do i=1,blockID%nBocos

       ii = bsearchIntegers(bcPrior(i), bcPriorSort(1:nDiff))

       ! Update the lower boundary in the multiplicity for this
       ! boundary condition and store the entry a bit easier.

       mult(ii-1) = mult(ii-1) + 1
       ii         = mult(ii-1)

       ! Store the mapping of the new to old subface id.

       oldSubfaceID(ii) = i

    enddo

  end subroutine sortSubfaces

  subroutine determineComputeBlocks(splitInfo)
    !
    !       determineComputeBlocks determines the computational blocks
    !       from the original grid and the given information how to split
    !       these blocks.
    !
    use constants
    use cgnsGrid, only :cgnsnDom, cgnsDoms
    use partitionMod, only : blocks, nBlocks, splitCGNSType
    use utils, only : terminate
    implicit none
    !
    !      Subroutine arguments.
    !
    type(splitCGNSType), dimension(cgnsNDom), intent(in) :: splitInfo
    !
    !      Local variables.
    !
    integer :: ierr

    integer(kind=intType) :: i, ii, jj, mm
    integer(kind=intType) :: nx, ny, nz, nAlloc

    integer(kind=intType), dimension(0:cgnsNDom) :: nSubPerCGNS

    ! Determine the number of subblocks per cgns block in
    ! cumulative storage format.

    nSubPerCGNS(0) = 0
    do i=1,cgnsNDom
       nSubPerCGNS(i) = nSubPerCGNS(i-1) + splitInfo(i)%nSubblocks
    enddo

    ! Check whether blocks are already allocated. If so, release
    ! the memory.
    if( allocated(blocks) ) then

       ! Loop over the number of old blocks and release the memory.

       mm = ubound(blocks,1)
       do i = 1,mm
          deallocate(blocks(i)%BCType,      blocks(i)%BCFaceID,   &
               blocks(i)%cgnsSubface, blocks(i)%inBeg,      &
               blocks(i)%jnBeg,       blocks(i)%knBeg,      &
               blocks(i)%inEnd,       blocks(i)%jnEnd,      &
               blocks(i)%knEnd,       blocks(i)%dinBeg,     &
               blocks(i)%djnBeg,      blocks(i)%dknBeg,     &
               blocks(i)%dinEnd,      blocks(i)%djnEnd,     &
               blocks(i)%dknEnd,      blocks(i)%neighBlock, &
               blocks(i)%l1,          blocks(i)%l2,         &
               blocks(i)%l3,          blocks(i)%groupNum,   &
               stat=ierr)
          if(ierr /= 0) &
               call terminate("determineComputeBlocks", &
               "Deallocation error for subface info")
       enddo

       ! Release the memory of blocks itself.

       deallocate(blocks, stat=ierr)
       if(ierr /= 0) &
            call terminate("determineComputeBlocks", &
            "Deallocation error for blocks")
    endif

    ! Allocate the memory for blocks.

    allocate(blocks(nBlocks), stat=ierr)
    if(ierr /= 0) &
         call terminate("determineComputeBlocks", &
         "Memory allocation failure for blocks")

    ! Set the counter ii for the global computational block number
    ! and loop over the cgns blocks.

    ii = 0
    cgnsLoop: do i=1,cgnsNDom

       ! Loop over the number of subblocks of this cgns block.

       subblockLoop: do mm=1,splitInfo(i)%nSubblocks

          ! Update the counter ii and store the number of cells in
          ! the three directions in nx, ny and nz.

          ii = ii + 1
          nx = splitInfo(i)%ranges(mm,1,2) &
               - splitInfo(i)%ranges(mm,1,1)
          ny = splitInfo(i)%ranges(mm,2,2) &
               - splitInfo(i)%ranges(mm,2,1)
          nz = splitInfo(i)%ranges(mm,3,2) &
               - splitInfo(i)%ranges(mm,3,1)

          ! Initialize the scalar variables of blocks(ii).

          blocks(ii)%nx = nx; blocks(ii)%il = nx + 1
          blocks(ii)%ny = ny; blocks(ii)%jl = ny + 1
          blocks(ii)%nz = nz; blocks(ii)%kl = nz + 1

          blocks(ii)%ncell = nx*ny*nz
          blocks(ii)%nface = (nx+1)*ny*nz + (ny+1)*nx*nz &
               + (nz+1)*nx*ny

          blocks(ii)%cgnsBlockID = i

          blocks(ii)%iBegor = splitInfo(i)%ranges(mm,1,1)
          blocks(ii)%jBegor = splitInfo(i)%ranges(mm,2,1)
          blocks(ii)%kBegor = splitInfo(i)%ranges(mm,3,1)

          blocks(ii)%iEndor = splitInfo(i)%ranges(mm,1,2)
          blocks(ii)%jEndor = splitInfo(i)%ranges(mm,2,2)
          blocks(ii)%kEndor = splitInfo(i)%ranges(mm,3,2)

          blocks(ii)%nBocos   = 0
          blocks(ii)%nSubface = 0
          blocks(ii)%n1to1    = 0

          ! Do an allocation for the subface info. NAlloc is such
          ! that no reallocation is needed for the boundary info.

          nAlloc = cgnsDoms(i)%nBocos + cgnsDoms(i)%n1to1

          allocate(blocks(ii)%BCType(nAlloc),      &
               blocks(ii)%BCFaceID(nAlloc),    &
               blocks(ii)%cgnsSubface(nAlloc), &
               blocks(ii)%inBeg(nAlloc),       &
               blocks(ii)%jnBeg(nAlloc),       &
               blocks(ii)%knBeg(nAlloc),       &
               blocks(ii)%inEnd(nAlloc),       &
               blocks(ii)%jnEnd(nAlloc),       &
               blocks(ii)%knEnd(nAlloc),       &
               blocks(ii)%dinBeg(nAlloc),      &
               blocks(ii)%djnBeg(nAlloc),      &
               blocks(ii)%dknBeg(nAlloc),      &
               blocks(ii)%dinEnd(nAlloc),      &
               blocks(ii)%djnEnd(nAlloc),      &
               blocks(ii)%dknEnd(nAlloc),      &
               blocks(ii)%neighBlock(nAlloc),  &
               blocks(ii)%l1(nAlloc),          &
               blocks(ii)%l2(nAlloc),          &
               blocks(ii)%l3(nAlloc),          &
               blocks(ii)%groupNum(nAlloc),    &
               stat=ierr)
          if(ierr /= 0) &
               call terminate("determineComputeBlocks", &
               "Memory allocation failure for &
               &subface info")

          ! Determine the boundary condition subfaces and the subfaces
          ! of the subblock, which are located on the block boundaries
          ! of the original cgns block.

          jj = 0
          call BCFacesSubblock(i, ii, jj)
          blocks(ii)%nBocos = jj

          call externalFacesSubblock(i, ii, jj, nSubPerCGNS, &
               nAlloc, splitInfo)

          ! Determine the subfaces of the subblock created by the
          ! splitting of the original block.

          call internalFacesSubblock(i, ii, jj, nSubPerCGNS, &
               nAlloc, splitInfo(i))
          blocks(ii)%nSubface = jj
          blocks(ii)%n1to1    = jj - blocks(ii)%nBocos

       enddo subblockLoop
    enddo cgnsLoop

  end subroutine determineComputeBlocks

  !========================================================================

  subroutine BCFacesSubblock(cgnsID, ii, jj)
    !
    !       BCFacesSubblock determines the boundary subfaces of compute
    !       block ii, which is a subblock of the given cgns block.
    !       Jj is the counter for the number of subfaces.
    !
    use cgnsGrid
    use communication
    use inputPhysics
    use partitionMod
    use utils, only : terminate
    implicit none
    !
    !      Subroutine arguments.
    !
    integer(kind=intType), intent(in)    :: cgnsID, ii
    integer(kind=intType), intent(inout) :: jj
    !
    !      Local variables.
    !
    integer :: ierr

    integer(kind=intType) :: j, mm
    integer(kind=intType) :: iBeg, iEnd, jBeg, jEnd, kBeg, kEnd

    character(len=maxCGNSNameLen) :: zoneName, subName
    character(len=2*maxStringLen) :: errorMessage

    ! Loop over the physical boundaries of the original block.

    bocoLoop: do j=1,cgnsDoms(cgnsID)%nBocos

       ! Continue with the next boundary subface if this is a
       ! degenerated subface.

       if(.not. cgnsDoms(cgnsID)%bocoInfo(j)%actualFace) cycle

       ! Store the subface range a bit easier. Make sure that the
       ! indices run from low to high.

       iBeg = min(cgnsDoms(cgnsID)%bocoInfo(j)%iBeg, &
            cgnsDoms(cgnsID)%bocoInfo(j)%iEnd)
       iEnd = max(cgnsDoms(cgnsID)%bocoInfo(j)%iBeg, &
            cgnsDoms(cgnsID)%bocoInfo(j)%iEnd)

       jBeg = min(cgnsDoms(cgnsID)%bocoInfo(j)%jBeg, &
            cgnsDoms(cgnsID)%bocoInfo(j)%jEnd)
       jEnd = max(cgnsDoms(cgnsID)%bocoInfo(j)%jBeg, &
            cgnsDoms(cgnsID)%bocoInfo(j)%jEnd)

       kBeg = min(cgnsDoms(cgnsID)%bocoInfo(j)%kBeg, &
            cgnsDoms(cgnsID)%bocoInfo(j)%kEnd)
       kEnd = max(cgnsDoms(cgnsID)%bocoInfo(j)%kBeg, &
            cgnsDoms(cgnsID)%bocoInfo(j)%kEnd)

       ! Check for a possible overlap between the current boundary
       ! subface and subblock ii.

       overlap: if(iBeg <= blocks(ii)%iEndor .and. &
            iEnd >= blocks(ii)%iBegor .and. &
            jBeg <= blocks(ii)%jEndor .and. &
            jEnd >= blocks(ii)%jBegor .and. &
            kBeg <= blocks(ii)%kEndor .and. &
            kEnd >= blocks(ii)%kBegor) then

          ! Determine the overlap region between the current boundary
          ! face and subblock ii.

          iBeg = max(blocks(ii)%iBegor, iBeg)
          iEnd = min(blocks(ii)%iEndor, iEnd)

          jBeg = max(blocks(ii)%jBegor, jBeg)
          jEnd = min(blocks(ii)%jEndor, jEnd)

          kBeg = max(blocks(ii)%kBegor, kBeg)
          kEnd = min(blocks(ii)%kEndor, kEnd)

          ! Check the number of equal indices, which is stored in mm.

          mm = 0
          if(iBeg == iEnd) mm = mm + 1
          if(jBeg == jEnd) mm = mm + 1
          if(kBeg == kEnd) mm = mm + 1

          ! If no constant index is found something is wrong with the
          ! grid. Processor 0 prints an error message, while the
          ! others wait until they are killed.

          if(mm == 0) then
             if(myID == 0) then
                zoneName = cgnsDoms(cgnsID)%zoneName
                subName  = cgnsDoms(cgnsID)%bocoInfo(j)%bocoName
                write(errorMessage,100) trim(zoneName), trim(subName)
100             format("Zone",1X,A,", boundary subface",1X,A, &
                     ": No constant index found for subface")
                call terminate("BCFacesSubblock", errorMessage)
             endif

             call mpi_barrier(ADflow_comm_world, ierr)
          endif

          ! Continue with the next subface if there is more than
          ! one constant index. This means that there is no overlap,
          ! but just an adjacency.

          if(mm > 1) cycle

          ! Update the counter jj and determine the range of the
          ! subface in the subblock.

          jj = jj + 1

          blocks(ii)%inBeg(jj) = iBeg - blocks(ii)%iBegor + 1
          blocks(ii)%inEnd(jj) = iEnd - blocks(ii)%iBegor + 1

          blocks(ii)%jnBeg(jj) = jBeg - blocks(ii)%jBegor + 1
          blocks(ii)%jnEnd(jj) = jEnd - blocks(ii)%jBegor + 1

          blocks(ii)%knBeg(jj) = kBeg - blocks(ii)%kBegor + 1
          blocks(ii)%knEnd(jj) = kEnd - blocks(ii)%kBegor + 1

          ! Determine the block face id on which this subface
          ! is located.

          if(iBeg == iEnd) then

             blocks(ii)%BCFaceID(jj) = iMax
             if(iBeg == blocks(ii)%iBegor) blocks(ii)%BCFaceID(jj) = iMin

          else if(jBeg == jEnd) then

             blocks(ii)%BCFaceID(jj) = jMax
             if(jBeg == blocks(ii)%jBegor) blocks(ii)%BCFaceID(jj) = jMin

          else

             blocks(ii)%BCFaceID(jj) = kMax
             if(kBeg == blocks(ii)%kBegor) blocks(ii)%BCFaceID(jj) = kMin

          endif

          ! Set some variables to 0, which are not relevant
          ! for boundary subfaces.

          blocks(ii)%dinBeg(jj) = 0; blocks(ii)%dinEnd(jj) = 0
          blocks(ii)%djnBeg(jj) = 0; blocks(ii)%djnEnd(jj) = 0
          blocks(ii)%dknBeg(jj) = 0; blocks(ii)%dknEnd(jj) = 0

          blocks(ii)%neighBlock(jj) = 0

          blocks(ii)%l1(jj) = 0
          blocks(ii)%l2(jj) = 0
          blocks(ii)%l3(jj) = 0

          ! Set the boundary condition and store to which original cgns
          ! subface this subface belongs.

          blocks(ii)%BCType(jj) = cgnsDoms(cgnsID)%bocoInfo(j)%BCType

          blocks(ii)%cgnsSubface(jj) = j

          ! Check whether this is a valid boundary condition for
          ! the current simulation.

          if(blocks(ii)%BCType(jj) == BCNotValid) then

             ! To avoid a messy output only processor 0 calls
             ! terminate. The other processors will wait until
             ! they are killed.

             if(myID == 0) then
                zoneName = cgnsDoms(cgnsID)%zoneName
                subName  = cgnsDoms(cgnsID)%bocoInfo(j)%bocoName

                ! Check whether this is an internal or an external
                ! flow problem and create the error message
                ! accordingly.

                if(flowType == internalFlow) then
                   write(errorMessage,120) trim(zoneName), trim(subName)
120                format("Zone",1X,A,", boundary subface",1X,A, &
                        ": Not a valid boundary condition for &
                        &internal flow")
                else
                   write(errorMessage,130) trim(zoneName), trim(subName)
130                format("Zone",1X,A,", boundary subface",1X,A, &
                        ": Not a valid boundary condition for &
                        &external flow")
                endif

                call terminate("BCFacesSubblock", errorMessage)
             endif

             call mpi_barrier(ADflow_comm_world, ierr)
          endif

          ! Store the corresponding family a bit easier.

          mm = cgnsDoms(cgnsID)%bocoInfo(j)%familyID

          ! Check if this is either a sliding mesh interface or a
          ! bleed flow region. If so the group nummer is set to the
          ! sliding interface ID or bleed flow region ID respectivily.
          ! Otherwise the group nummer is the family nummer, which
          ! is 0 if the subface does not belong to a family.

          select case (blocks(ii)%BCType(jj))
          case (SlidingInterface)
             blocks(ii)%groupNum(jj) = &
                  cgnsDoms(cgnsID)%bocoInfo(j)%slidingID

          case (MassBleedInflow, MassBleedOutflow)
             blocks(ii)%groupNum(jj) = cgnsFamilies(mm)%bleedRegionID

          case default
             blocks(ii)%groupNum(jj) = mm
          end select

       endif overlap

    enddo bocoLoop

  end subroutine BCFacesSubblock

  !========================================================================

  subroutine externalFacesSubblock(cgnsID, ii, jj, nSubPerCGNS, &
       nAlloc, splitInfo)
    !
    !       externalFacesSubblock determines the block boundaries of
    !       the compute block ii which are located on the boundaries of
    !       the given original cgns block. As it is possible that due to
    !       a splitting of a neighboring block the number of block
    !       boundaries is larger than the original number, it must be
    !       checked whether enough memory has been allocated.
    !       jj is the counter for the number of subfaces.
    !
    use cgnsGrid
    use communication
    use partitionMod
    use utils, only : delta, terminate
    implicit none
    !
    !      Subroutine arguments
    !
    integer(kind=intType), intent(in)    :: cgnsID, ii
    integer(kind=intType), intent(inout) :: jj, nAlloc

    integer(kind=intType), dimension(0:cgnsNDom), intent(in) :: &
         nSubPerCGNS
    type(splitCGNSType), dimension(cgnsNDom), intent(in) :: splitInfo
    !
    !      Local variables.
    !
    integer :: ierr

    integer(kind=intType) :: j, k, kk, mm
    integer(kind=intType) :: l1, L2, l3
    integer(kind=intType) :: iBeg, iEnd, jBeg, jEnd, kBeg, kEnd
    integer(kind=intType) :: diBeg, diEnd, djBeg, djEnd
    integer(kind=intType) :: dkBeg, dkEnd

    integer(kind=intType), dimension(3,3) :: trMat

    integer(kind=intType), dimension(:,:,:), pointer :: ranges

    character(len=maxCGNSNameLen) :: zoneName, subName
    character(len=2*maxStringLen) :: errorMessage

    logical :: diSwap, djSwap, dkSwap
    !
    ! Loop over the 1 to 1 block boundaries of the original block.

    n1to1Loop: do j=1,cgnsDoms(cgnsID)%n1to1

       ! Store the subface range a bit easier. Make sure that the
       ! indices run from low to high.

       iBeg = min(cgnsDoms(cgnsID)%conn1to1(j)%iBeg, &
            cgnsDoms(cgnsID)%conn1to1(j)%iEnd)
       iEnd = max(cgnsDoms(cgnsID)%conn1to1(j)%iBeg, &
            cgnsDoms(cgnsID)%conn1to1(j)%iEnd)

       jBeg = min(cgnsDoms(cgnsID)%conn1to1(j)%jBeg, &
            cgnsDoms(cgnsID)%conn1to1(j)%jEnd)
       jEnd = max(cgnsDoms(cgnsID)%conn1to1(j)%jBeg, &
            cgnsDoms(cgnsID)%conn1to1(j)%jEnd)

       kBeg = min(cgnsDoms(cgnsID)%conn1to1(j)%kBeg, &
            cgnsDoms(cgnsID)%conn1to1(j)%kEnd)
       kEnd = max(cgnsDoms(cgnsID)%conn1to1(j)%kBeg, &
            cgnsDoms(cgnsID)%conn1to1(j)%kEnd)

       ! Check for a possible overlap between the current boundary
       ! subface and subblock ii.

       overlap: if(iBeg <= blocks(ii)%iEndor .and. &
            iEnd >= blocks(ii)%iBegor .and. &
            jBeg <= blocks(ii)%jEndor .and. &
            jEnd >= blocks(ii)%jBegor .and. &
            kBeg <= blocks(ii)%kEndor .and. &
            kEnd >= blocks(ii)%kBegor) then

          ! Determine the overlap region between the current boundary
          ! face and subblock ii.

          iBeg = max(blocks(ii)%iBegor, iBeg)
          iEnd = min(blocks(ii)%iEndor, iEnd)

          jBeg = max(blocks(ii)%jBegor, jBeg)
          jEnd = min(blocks(ii)%jEndor, jEnd)

          kBeg = max(blocks(ii)%kBegor, kBeg)
          kEnd = min(blocks(ii)%kEndor, kEnd)

          ! Check the number of equal indices, which is stored in kk.

          kk = 0
          if(iBeg == iEnd) kk = kk + 1
          if(jBeg == jEnd) kk = kk + 1
          if(kBeg == kEnd) kk = kk + 1

          ! If no constant index is found something is wrong with the
          ! grid. Processor 0 prints an error message, while the
          ! others wait until they are killed.

          if(kk == 0) then
             if(myID == 0) then
                zoneName = cgnsDoms(cgnsID)%zoneName
                subName  = cgnsDoms(cgnsID)%bocoInfo(j)%bocoName
                write(errorMessage,140) trim(zoneName), trim(subName)
140             format("Zone",1X,A,", 1 to 1 block connectivity",1X,A, &
                     ": No constant index found for subface")
                call terminate("externalFacesSubblock", errorMessage)
             endif

             call mpi_barrier(ADflow_comm_world, ierr)
          endif

          ! Continue with the next subface if there is more than
          ! one constant index. This means that there is no overlap,
          ! but just an adjacency.

          if(kk > 1) cycle

          ! Preserve negative running indices of the subface.

          if(cgnsDoms(cgnsID)%conn1to1(j)%iEnd < &
               cgnsDoms(cgnsID)%conn1to1(j)%iBeg) then
             mm   = iBeg
             iBeg = iEnd
             iEnd = mm
          endif

          if(cgnsDoms(cgnsID)%conn1to1(j)%jEnd < &
               cgnsDoms(cgnsID)%conn1to1(j)%jBeg) then
             mm   = jBeg
             jBeg = jEnd
             jEnd = mm
          endif

          if(cgnsDoms(cgnsID)%conn1to1(j)%kEnd < &
               cgnsDoms(cgnsID)%conn1to1(j)%kBeg) then
             mm   = kBeg
             kBeg = kEnd
             kEnd = mm
          endif

          ! Determine the transformation matrix between the
          ! current subface and the donor subface.

          l1 = cgnsDoms(cgnsID)%conn1to1(j)%l1
          L2 = cgnsDoms(cgnsID)%conn1to1(j)%l2
          l3 = cgnsDoms(cgnsID)%conn1to1(j)%l3

          trMat(1,1) = sign(1_intType,l1) * delta(l1,1_intType)
          trMat(2,1) = sign(1_intType,l1) * delta(l1,2_intType)
          trMat(3,1) = sign(1_intType,l1) * delta(l1,3_intType)

          trMat(1,2) = sign(1_intType,l2) * delta(l2,1_intType)
          trMat(2,2) = sign(1_intType,l2) * delta(l2,2_intType)
          trMat(3,2) = sign(1_intType,l2) * delta(l2,3_intType)

          trMat(1,3) = sign(1_intType,l3) * delta(l3,1_intType)
          trMat(2,3) = sign(1_intType,l3) * delta(l3,2_intType)
          trMat(3,3) = sign(1_intType,l3) * delta(l3,3_intType)

          ! Determine the corresponding donor range of the subface
          ! iBeg, iEnd; jBeg, jEnd; kBeg, kEnd.

          l1 = iBeg - cgnsDoms(cgnsID)%conn1to1(j)%iBeg
          L2 = jBeg - cgnsDoms(cgnsID)%conn1to1(j)%jBeg
          l3 = kBeg - cgnsDoms(cgnsID)%conn1to1(j)%kBeg

          diBeg = cgnsDoms(cgnsID)%conn1to1(j)%diBeg &
               + trMat(1,1)*l1 + trMat(1,2)*l2 + trMat(1,3)*l3
          djBeg = cgnsDoms(cgnsID)%conn1to1(j)%djBeg &
               + trMat(2,1)*l1 + trMat(2,2)*l2 + trMat(2,3)*l3
          dkBeg = cgnsDoms(cgnsID)%conn1to1(j)%dkBeg &
               + trMat(3,1)*l1 + trMat(3,2)*l2 + trMat(3,3)*l3

          l1 = iEnd - cgnsDoms(cgnsID)%conn1to1(j)%iBeg
          L2 = jEnd - cgnsDoms(cgnsID)%conn1to1(j)%jBeg
          l3 = kEnd - cgnsDoms(cgnsID)%conn1to1(j)%kBeg

          diEnd = cgnsDoms(cgnsID)%conn1to1(j)%diBeg &
               + trMat(1,1)*l1 + trMat(1,2)*l2 + trMat(1,3)*l3
          djEnd = cgnsDoms(cgnsID)%conn1to1(j)%djBeg &
               + trMat(2,1)*l1 + trMat(2,2)*l2 + trMat(2,3)*l3
          dkEnd = cgnsDoms(cgnsID)%conn1to1(j)%dkBeg &
               + trMat(3,1)*l1 + trMat(3,2)*l2 + trMat(3,3)*l3

          ! Make sure that the donor indices are positive running
          ! indices. If they must be swapped, the corresponding
          ! logical is set to .true.

          diSwap = .false.
          if(diBeg > diEnd) then
             mm = diBeg; diBeg = diEnd; diEnd = mm; diSwap = .true.
          endif

          djSwap = .false.
          if(djBeg > djEnd) then
             mm = djBeg; djBeg = djEnd; djEnd = mm; djSwap = .true.
          endif

          dkSwap = .false.
          if(dkBeg > dkEnd) then
             mm = dkBeg; dkBeg = dkEnd; dkEnd = mm; dkSwap = .true.
          endif

          ! Store the index of the donor block a bit easier and loop
          ! over its subblocks to find the donor range.

          mm = cgnsDoms(cgnsID)%conn1to1(j)%donorBlock
          ranges => splitInfo(mm)%ranges

          donorLoop: do k=1,splitInfo(mm)%nSubblocks

             ! Check whether this subblock and the given donor range
             ! overlap.

             donorOverlap: if(diBeg <= ranges(k,1,2) .and. &
                  diEnd >= ranges(k,1,1) .and. &
                  djBeg <= ranges(k,2,2) .and. &
                  djEnd >= ranges(k,2,1) .and. &
                  dkBeg <= ranges(k,3,2) .and. &
                  dkEnd >= ranges(k,3,1)) then

                ! Determine the range of the donor face, which is
                ! stored in iBeg, iEnd, etc.

                iBeg = max(diBeg,ranges(k,1,1))
                iEnd = min(diEnd,ranges(k,1,2))

                jBeg = max(djBeg,ranges(k,2,1))
                jEnd = min(djEnd,ranges(k,2,2))

                kBeg = max(dkBeg,ranges(k,3,1))
                kEnd = min(dkEnd,ranges(k,3,2))

                ! Check whether the subfaces are truely overlapping
                ! or just adjacent.

                kk = 0
                if(iBeg == iEnd) kk = kk + 1
                if(jBeg == jEnd) kk = kk + 1
                if(kBeg == kEnd) kk = kk + 1

                if(kk > 1) cycle

                ! Update the counter jj and check whether enough memory
                ! has been allocated. If not, reallocate.

                jj = jj + 1
                if(jj > nAlloc) call reallocSubfaceMemory(ii,nAlloc)

                ! Set some info for this subface, which can be
                ! determined relatively easily.

                blocks(ii)%BCType(jj)      = B2BMatch
                blocks(ii)%cgnsSubface(jj) = j
                blocks(ii)%groupNum(jj)    = 0

                blocks(ii)%l1(jj) = cgnsDoms(cgnsID)%conn1to1(j)%l1
                blocks(ii)%l2(jj) = cgnsDoms(cgnsID)%conn1to1(j)%l2
                blocks(ii)%l3(jj) = cgnsDoms(cgnsID)%conn1to1(j)%l3

                ! Determine the neighboring block id.

                blocks(ii)%neighBlock(jj) = nSubPerCGNS(mm-1) + k

                ! Determine the range of the donor. First switch the
                ! indices if the original indices were swapped.

                if( diSwap ) then
                   kk = iBeg; iBeg = iEnd; iEnd = kk
                endif

                if( djSwap ) then
                   kk = jBeg; jBeg = jEnd; jEnd = kk
                endif

                if( dkSwap ) then
                   kk = kBeg; kBeg = kEnd; kEnd = kk
                endif

                ! Determine the local range of the donor, i.e. the
                ! offset in the original block must be substracted.

                blocks(ii)%dinBeg(jj) = iBeg - ranges(k,1,1) + 1
                blocks(ii)%djnBeg(jj) = jBeg - ranges(k,2,1) + 1
                blocks(ii)%dknBeg(jj) = kBeg - ranges(k,3,1) + 1

                blocks(ii)%dinEnd(jj) = iEnd - ranges(k,1,1) + 1
                blocks(ii)%djnEnd(jj) = jEnd - ranges(k,2,1) + 1
                blocks(ii)%dknEnd(jj) = kEnd - ranges(k,3,1) + 1

                ! Transform the donor range in the original donor block
                ! back the a subface range in the original cgns block.
                ! The inverse of the transformation matrix trMat is
                ! the transpose.

                l1 = iBeg - cgnsDoms(cgnsID)%conn1to1(j)%diBeg
                L2 = jBeg - cgnsDoms(cgnsID)%conn1to1(j)%djBeg
                l3 = kBeg - cgnsDoms(cgnsID)%conn1to1(j)%dkBeg

                iBeg = cgnsDoms(cgnsID)%conn1to1(j)%iBeg &
                     + trMat(1,1)*l1 + trMat(2,1)*l2 + trMat(3,1)*l3
                jBeg = cgnsDoms(cgnsID)%conn1to1(j)%jBeg &
                     + trMat(1,2)*l1 + trMat(2,2)*l2 + trMat(3,2)*l3
                kBeg = cgnsDoms(cgnsID)%conn1to1(j)%kBeg &
                     + trMat(1,3)*l1 + trMat(2,3)*l2 + trMat(3,3)*l3

                l1 = iEnd - cgnsDoms(cgnsID)%conn1to1(j)%diBeg
                L2 = jEnd - cgnsDoms(cgnsID)%conn1to1(j)%djBeg
                l3 = kEnd - cgnsDoms(cgnsID)%conn1to1(j)%dkBeg

                iEnd = cgnsDoms(cgnsID)%conn1to1(j)%iBeg &
                     + trMat(1,1)*l1 + trMat(2,1)*l2 + trMat(3,1)*l3
                jEnd = cgnsDoms(cgnsID)%conn1to1(j)%jBeg &
                     + trMat(1,2)*l1 + trMat(2,2)*l2 + trMat(3,2)*l3
                kEnd = cgnsDoms(cgnsID)%conn1to1(j)%kBeg &
                     + trMat(1,3)*l1 + trMat(2,3)*l2 + trMat(3,3)*l3

                ! Store the subface range of the new block, i.e.
                ! An offset must be subtracted.

                blocks(ii)%inBeg(jj) = iBeg - blocks(ii)%iBegor + 1
                blocks(ii)%jnBeg(jj) = jBeg - blocks(ii)%jBegor + 1
                blocks(ii)%knBeg(jj) = kBeg - blocks(ii)%kBegor + 1

                blocks(ii)%inEnd(jj) = iEnd - blocks(ii)%iBegor + 1
                blocks(ii)%jnEnd(jj) = jEnd - blocks(ii)%jBegor + 1
                blocks(ii)%knEnd(jj) = kEnd - blocks(ii)%kBegor + 1

                ! Determine the block face id on which this subface
                ! is located.

                if(iBeg == iEnd) then

                   blocks(ii)%BCFaceID(jj) = iMax
                   if(iBeg == blocks(ii)%iBegor) &
                        blocks(ii)%BCFaceID(jj) = iMin

                else if(jBeg == jEnd) then

                   blocks(ii)%BCFaceID(jj) = jMax
                   if(jBeg == blocks(ii)%jBegor) &
                        blocks(ii)%BCFaceID(jj) = jMin

                else

                   blocks(ii)%BCFaceID(jj) = kMax
                   if(kBeg == blocks(ii)%kBegor) &
                        blocks(ii)%BCFaceID(jj) = kMin

                endif

             endif donorOverlap
          enddo donorLoop
       endif overlap
    enddo n1to1Loop

  end subroutine externalFacesSubblock

  !========================================================================

  subroutine internalFacesSubblock(cgnsID, ii, jj, nSubPerCGNS, &
       nAlloc, splitInfo)
    !
    !       internalFacesSubblock determines the block boundaries of
    !       the compute block ii which are created due to the splitting of
    !       the original block into subblock. As the number of these
    !       internal boundaries is not known, it must be checked whether
    !       enough memory has been allocated. jj is the counter for the
    !       number of subfaces.
    !
    use cgnsGrid
    use communication
    use partitionMod
    implicit none
    !
    !      Subroutine arguments
    !
    integer(kind=intType), intent(in)    :: cgnsID, ii
    integer(kind=intType), intent(inout) :: jj, nAlloc

    integer(kind=intType), dimension(0:cgnsNDom), intent(in) :: &
         nSubPerCGNS
    type(splitCGNSType), intent(in) :: splitInfo
    !
    !      Local variables.
    !
    integer(kind=intType) :: indFace, jBeg, jEnd, kBeg, kEnd
    integer(kind=intType) :: i, i2, j, k, faceID

    ! iMin face.

    if(blocks(ii)%iBegor > 1) then

       ! Imin face is created through splitting. Set some variables
       ! for the general treatment.

       indFace = blocks(ii)%iBegor
       jBeg    = blocks(ii)%jBegor
       jEnd    = blocks(ii)%jEndor
       kBeg    = blocks(ii)%kBegor
       kEnd    = blocks(ii)%kEndor

       i = 1; j = 2; k = 3
       i2 = 2; faceID = iMin

       ! Search for neighbors in the subblocks of the given
       ! cgns blocks.

       call searchInternalNeighbors

    endif

    ! iMax face.

    if(blocks(ii)%iEndor < cgnsDoms(cgnsID)%il) then

       ! Imax face is created through splitting. Set some variables
       ! for the general treatment.

       indFace = blocks(ii)%iEndor
       jBeg    = blocks(ii)%jBegor
       jEnd    = blocks(ii)%jEndor
       kBeg    = blocks(ii)%kBegor
       kEnd    = blocks(ii)%kEndor

       i = 1; j = 2; k = 3
       i2 = 1; faceID = iMax

       ! Search for neighbors in the subblocks of the given
       ! cgns blocks.

       call searchInternalNeighbors

    endif

    ! jMin face.

    if(blocks(ii)%jBegor > 1) then

       ! Jmin face is created through splitting. Set some variables
       ! for the general treatment.

       indFace = blocks(ii)%jBegor
       jBeg    = blocks(ii)%iBegor
       jEnd    = blocks(ii)%iEndor
       kBeg    = blocks(ii)%kBegor
       kEnd    = blocks(ii)%kEndor

       i = 2; j = 1; k = 3
       i2 = 2; faceID = jMin

       ! Search for neighbors in the subblocks of the given
       ! cgns blocks.

       call searchInternalNeighbors

    endif

    ! jMax face.

    if(blocks(ii)%jEndor < cgnsDoms(cgnsID)%jl) then

       ! Jmax face is created through splitting. Set some variables
       ! for the general treatment.

       indFace = blocks(ii)%jEndor
       jBeg    = blocks(ii)%iBegor
       jEnd    = blocks(ii)%iEndor
       kBeg    = blocks(ii)%kBegor
       kEnd    = blocks(ii)%kEndor

       i = 2; j = 1; k = 3
       i2 = 1; faceID = jMax

       ! Search for neighbors in the subblocks of the given
       ! cgns blocks.

       call searchInternalNeighbors

    endif

    ! kMin face.

    if(blocks(ii)%kBegor > 1) then

       ! Kmin face is created through splitting. Set some variables
       ! for the general treatment.

       indFace = blocks(ii)%kBegor
       jBeg    = blocks(ii)%iBegor
       jEnd    = blocks(ii)%iEndor
       kBeg    = blocks(ii)%jBegor
       kEnd    = blocks(ii)%jEndor

       i = 3; j = 1; k = 2
       i2 = 2; faceID = kMin

       ! Search for neighbors in the subblocks of the given
       ! cgns blocks.

       call searchInternalNeighbors

    endif

    ! kMax face.

    if(blocks(ii)%kEndor < cgnsDoms(cgnsID)%kl) then

       ! Kmax face is created through splitting. Set some variables
       ! for the general treatment.

       indFace = blocks(ii)%kEndor
       jBeg    = blocks(ii)%iBegor
       jEnd    = blocks(ii)%iEndor
       kBeg    = blocks(ii)%jBegor
       kEnd    = blocks(ii)%jEndor

       i = 3; j = 1; k = 2
       i2 = 1; faceID = kMax

       ! Search for neighbors in the subblocks of the given
       ! cgns blocks.

       call searchInternalNeighbors

    endif

    !=================================================================

  contains

    !===============================================================

    subroutine searchInternalNeighbors
      !
      !         searchInternalNeighbors determines block faces created by
      !         by the splitting of the original block. The variables set in
      !         internalFacesSubblock are used such that a general
      !         treatment is possible.
      !
      implicit none
      !
      !        Local variables
      !
      integer(kind=intType) :: mm, jnBeg, jnEnd, knBeg, knEnd

      integer(kind=intType), dimension(3,2) :: subRange

      integer(kind=intType), dimension(:,:,:), pointer :: ranges

      ! Set the pointer for ranges to make the code more readable.

      ranges => splitInfo%ranges

      ! Loop over the number of blocks into the original block
      ! is split.

      do mm=1,splitInfo%nSubblocks

         ! Check whether the constant index of the face matches
         ! the given index of subblock.

         if(ranges(mm,i,i2) == indFace) then

            ! Check whether the faces overlap.

            if(jBeg <= ranges(mm,j,2) .and. &
                 jEnd >= ranges(mm,j,1) .and. &
                 kBeg <= ranges(mm,k,2) .and. &
                 kEnd >= ranges(mm,k,1) ) then

               ! There is a possible overlap. Determine the nodal
               ! range of the subface.

               jnBeg = max(jBeg,ranges(mm,j,1))
               jnEnd = min(jEnd,ranges(mm,j,2))

               knBeg = max(kBeg,ranges(mm,k,1))
               knEnd = min(kEnd,ranges(mm,k,2))

               ! Check whether this is a true subface.

               if(jnEnd > jnBeg .and. knEnd > knBeg) then

                  ! An overlap occurs. Update the counter jj and check
                  ! whether enough memory has been allocated.
                  ! If not, reallocate.

                  jj = jj + 1
                  if(jj > nAlloc) call reallocSubfaceMemory(ii,nAlloc)

                  ! Set the information of the BCType, BCFaceID,
                  ! cgnsSubface and the group number. As this face is
                  ! created internally the latter two variables are
                  ! set to 0.

                  blocks(ii)%BCType(jj)      = B2BMatch
                  blocks(ii)%BCFaceID(jj)    = faceID
                  blocks(ii)%cgnsSubface(jj) = 0
                  blocks(ii)%groupNum(jj)    = 0

                  ! Determine the subRange of the subface in the
                  ! original block.

                  subRange(i,1) = indFace; subRange(i,2) = indFace
                  subRange(j,1) = jnBeg;   subRange(j,2) = jnEnd
                  subRange(k,1) = knBeg;   subRange(k,2) = knEnd

                  ! Determine the nodal range in the current subblock.

                  blocks(ii)%inBeg(jj) = subRange(1,1) &
                       - blocks(ii)%iBegor + 1
                  blocks(ii)%inEnd(jj) = subRange(1,2) &
                       - blocks(ii)%iBegor + 1

                  blocks(ii)%jnBeg(jj) = subRange(2,1) &
                       - blocks(ii)%jBegor + 1
                  blocks(ii)%jnEnd(jj) = subRange(2,2) &
                       - blocks(ii)%jBegor + 1

                  blocks(ii)%knBeg(jj) = subRange(3,1) &
                       - blocks(ii)%kBegor + 1
                  blocks(ii)%knEnd(jj) = subRange(3,2) &
                       - blocks(ii)%kBegor + 1

                  ! Determine the nodal range in the donor block.

                  blocks(ii)%dinBeg(jj) = subRange(1,1) &
                       - ranges(mm,1,1) + 1
                  blocks(ii)%dinEnd(jj) = subRange(1,2) &
                       - ranges(mm,1,1) + 1

                  blocks(ii)%djnBeg(jj) = subRange(2,1) &
                       - ranges(mm,2,1) + 1
                  blocks(ii)%djnEnd(jj) = subRange(2,2) &
                       - ranges(mm,2,1) + 1

                  blocks(ii)%dknBeg(jj) = subRange(3,1) &
                       - ranges(mm,3,1) + 1
                  blocks(ii)%dknEnd(jj) = subRange(3,2) &
                       - ranges(mm,3,1) + 1

                  ! Set the neighboring block to mm plus the offset
                  ! for the current cgns block and set the transformation
                  ! matrix. The latter is simply 1-2-3, because the
                  ! orientation of the subblocks is identical to the
                  ! original block.

                  blocks(ii)%neighBlock(jj) = mm + nSubPerCGNS(cgnsID-1)

                  blocks(ii)%l1(jj) = 1
                  blocks(ii)%l2(jj) = 2
                  blocks(ii)%l3(jj) = 3

               endif
            endif
         endif
      enddo

    end subroutine searchInternalNeighbors

  end subroutine internalFacesSubblock
  subroutine graphPartitioning(emptyPartitions, commNeglected)
    !
    !       graphPartitioning partitions the corresponding graph of the
    !       computational blocks such that both the number of cells and
    !       number of faces is about equal on all processors.
    !
    use constants
    use communication, only : myID, adflow_comm_world, nProc
    use partitionMod, only : part, blocks, ubvec, nBlocks
    use inputParallel, only : loadImbalance
    use utils, only : terminate
    use sorting, only : qsortIntegers, bsearchIntegers
    implicit none
    !
    !      Subroutine arguments.
    !
    logical, intent(out) :: emptyPartitions, commNeglected
    !
    !      Variables to store the graph in metis format.
    !
    ! nVertex        : Number of vertices in the graph, equals nBlocks.
    ! nCon           : Number of contraints, 2.
    ! xadj(0:nVertex): Number of edges per vertex, cumulative storage
    !                  format.
    ! adjncy(:)      : End vertex of the edge; the size of adjncy is
    !                  xadj(nVertex).
    ! vwgt(:,:)      : Vertex weights, size equals nCon,nVertex. The
    !                  vertex weights are stored contiguously.
    ! adjwgt(:)      : Edge weights, size equals xadj(nVertex). Note
    !                  that the edge weights of edge i-j can be
    !                  different from the weight of edge j-i.
    ! wgtflag        : Whether or not to use weights on edges. Here
    !                  wgtflag should always be 1 to indicate that
    !                  edge weights are used.
    ! numflag        : Flag to indicate the numbering convention,
    !                  starting from 0 or 1. Here we start from 0.
    ! nParts         : Number of parts to split the graph. This is
    !                  nProc.
    ! ubvec(2)       : Tolerance for the constraints. Stored in the
    !                  module dpartitionMod.
    ! options(5)     : Option array; normally the default is used
    !                  indicated by options(1) = 0.
    ! edgecut        : On return it contains the edge cut of the
    !                  distributed graph.
    ! part(nVertex)  : On return the processor ID for each block.
    !                  It will be returned in fortran numbering,
    !                  i.e. starting at 1.  Stored in the module
    !                  distributionMod.

    integer :: nVertex, nCon, wgtflag, numflag, nParts, edgecut
    integer, dimension(5) :: options

    integer(kind=intType), dimension(:),   allocatable :: xadj, adjncy
    integer(kind=intType), dimension(:),   allocatable :: adjwgt
    integer(kind=intType), dimension(:,:), allocatable :: vwgt
    !
    !      Local variables.
    !
    integer :: ierr

    integer(kind=intType) :: i, j
    integer(kind=intType) :: nEdges, nEdgesMax, ii, jj, kk

    integer(kind=intType), dimension(0:nProc-1) :: nBlockPerProc

    integer(kind=intType), dimension(:), allocatable :: tmp

    integer(kind=8) :: nCellsTotal    ! 8 byte integers to avoid
    integer(kind=8) :: nFacesTotal    ! overflow.

    real(kind=4) :: ubvec_temp(2) ! Explict float for ubvec
    !
    ! Check whether part is allocated from a previous call. If so,
    ! release the memory.

    if( allocated(part) ) then
       deallocate(part, stat=ierr)
       if(ierr /= 0)                       &
            call terminate("graphPartitioning", &
            "Deallocation failure for part")
    endif

    ! Determine the number of edges in the graph and the maximum
    ! number for a vertex in the graph.

    nEdges = 0
    nEdgesMax = 0
    do i=1,nBlocks
       ii        = blocks(i)%n1to1
       nEdges    = nedges + ii
       nEdgesMax = max(nedgesMax, ii)
    enddo

    ! Initialize some values for the graph.

    nVertex  = nBlocks
    nCon     = 2
    wgtflag  = 1
    numflag  = 0
    nParts   = nProc
    ubvec(1) = one + loadImbalance
    ubvec(2) = one + loadImbalance
    options  = 0

    !  Allocate the memory to store/build the graph.

    allocate(xadj(0:nVertex), vwgt(nCon,nVertex), adjncy(nEdges), &
         adjwgt(nEdges), part(nVertex), tmp(nEdgesMax), stat=ierr)
    if(ierr /= 0)                         &
         call terminate("graphPartitioning", &
         "Memory allocation failure for graph variables")

    ! Initialize xadj(0) to 0.
    ! Furthermore initialize adjwgt to 0, as these values are
    ! accumulated due to multiple subfaces between blocks.

    xadj(0) = 0
    adjwgt  = 0

    ! Loop over the number of blocks to build the graph.

    graphVertex: do i=1,nBlocks

       ! Store the both vertex weights.

       vwgt(1,i) = blocks(i)%nCell
       vwgt(2,i) = blocks(i)%nFace

       ! Sort the neighbors in increasing order and neglect the
       ! communication to myself, i.e. do not allow an edge to myself.
       ! The sorting is necessary, because a block might have several
       ! subfaces with another block

       nEdges = 0
       do j=1,blocks(i)%n1to1
          ii = blocks(i)%nBocos + j
          if(blocks(i)%neighBlock(ii) /= i) then
             nEdges = nEdges +1
             tmp(nEdges) = blocks(i)%neighBlock(ii)
          endif
       enddo

       ! Sort tmp in increasing order and get rid of the possible
       ! multiple entries.

       call qsortIntegers(tmp, nEdges)

       ii = min(nEdges,1_intType)  ! Be aware of nEdges == 0
       do j=2,nEdges
          if(tmp(j) /= tmp(ii)) then
             ii = ii + 1
             tmp(ii) = tmp(j)
          endif
       enddo

       ! Set nEdges to ii and update xadj(i).

       nEdges  = ii
       xadj(i) = xadj(i-1) + nEdges

       ! Repeat the loop over the subfaces, but now store
       ! the edge info.

       Edges1to1: do j=1,blocks(i)%n1to1

          ii = blocks(i)%nBocos + j
          if(blocks(i)%neighBlock(ii) /= i) then

             ! Search for the block ID and add the offset of xadj(i-1)
             ! to obtain the correct index to store the edge info.
             ! The -1 to the adjncy is present because C-numbering
             ! is used when calling Metis.

             jj = xadj(i-1) &
                  + bsearchIntegers(blocks(i)%neighBlock(ii), tmp(1:nEdges))
             adjncy(jj) = blocks(i)%neighBlock(ii) - 1

             ! The weight equals the number of 1st and 2nd level halo
             ! cells to be communicated between the blocks. The weights
             ! are accumulated, as multiple subfaces between blocks are
             ! possible.

             kk = 1
             adjwgt(jj) = adjwgt(jj) + 2 * &
                  ( max(abs(blocks(i)%inEnd(ii) - blocks(i)%inBeg(ii)), kk) &
                  * max(abs(blocks(i)%jnEnd(ii) - blocks(i)%jnBeg(ii)), kk) &
                  * max(abs(blocks(i)%knEnd(ii) - blocks(i)%knBeg(ii)), kk) )
          endif

       enddo Edges1to1

    enddo graphVertex

    ! Metis has problems when the total number of cells or faces
    ! used in the weights exceeds 2Gb. Therefore the sum of these
    ! values is determined and an appropriate weight factor is
    ! determined. Note that the type of nCellsTotal and nFacesTotal
    ! is integer*8.

    nCellsTotal = 0
    nFacesTotal = 0
    do i=1,nBlocks
       nCellsTotal = nCellsTotal + vwgt(1,i)
       nFacesTotal = nFacesTotal + vwgt(2,i)
    enddo

    if(nCellsTotal > 2147483647 .or. nFacesTotal > 2147483647) then
       nCellsTotal = nCellsTotal/2147483647 + 1
       nFacesTotal = nFacesTotal/2147483647 + 1

       do i=1,nBlocks
          vwgt(1,i) = vwgt(1,i)/nCellsTotal
          vwgt(2,i) = vwgt(2,i)/nFacesTotal
       enddo
    endif

    ! Loop over the number of attempts to partition the graph.
    ! In the first attempt the communication is taken into account.
    ! If not successful, i.e. empty partitions present, the metis
    ! routine is called once more, but now with zero adjwgt. This
    ! means that the communication cost is neglected and metis
    ! normally gives a valid partitioning.
    ! Initialize commNeglected to .false. This will change if in
    ! the loop below the first call to metis is not successful.

    commNeglected = .false.
    attemptLoop: do ii=1,2

       ! Copy ubvec to a float since if you use -r8 flag it gets
       ! converted to real
       ubvec_temp(1) = real(ubvec(1))
       ubvec_temp(2) = real(ubvec(2))

       call metisInterface(nVertex, nCon, xadj, adjncy, vwgt, &
            adjwgt, wgtflag, numflag, nParts,  &
            ubvec_temp, options, edgecut, part)
       ! Determine the number of blocks per processor.

       nBlockPerProc = 0
       do i=1,nBlocks
          nBlockPerProc(part(i)) = nBlockPerProc(part(i)) + 1
       enddo

       ! Check for empty partitions.

       emptyPartitions = .false.
       do i=0,nProc-1
          if(nBlockPerProc(i) == 0) emptyPartitions = .true.
       enddo

       ! Exit the loop if no empty partitions are present or if
       ! this is the second time this loop is executed.

       if(ii == 2 .or. (.not. emptyPartitions)) exit attemptLoop

       ! The first call to metis resulted in empty partitions.
       ! Ignore the communication, i.e. set the number of
       ! neighbors to 0, and try again.

       commNeglected = .true.
       xadj          = 0

    enddo attemptLoop

    ! Deallocate the memory for the graph except part.

    deallocate(xadj, vwgt, adjncy, adjwgt, tmp, stat=ierr)
    if(ierr /= 0)                         &
         call terminate("graphPartitioning", &
         "Deallocation failure for graph variables")

  end subroutine graphPartitioning

  subroutine checkLoadBalance(cellsBalanced, facesBalanced)
    !
    !       checkLoadBalance determines whether or not the load balance
    !       for the cells and faces is met.
    !
    use constants
    use communication, only : adflow_comm_world, myid, nProc
    use partitionMod, only : blocks, ubvec, part, nBlocks
    implicit none
    !
    !      Subroutine arguments.
    !
    logical, intent(out) :: cellsBalanced, facesBalanced
    !
    !      Local variables
    !
    integer(kind=intType) :: i, j
    integer(kind=intType) :: nCellMax, nCellTol
    integer(kind=intType) :: nFaceMax, nFaceTol

    integer(kind=intType), dimension(nProc) :: nCell, nFace

    integer(kind=8) :: nCellsEven, nFacesEven   ! 8 byte integers to
    ! avoid overflow.

    ! Initialize nCell and nFace to 0. These variables will contain
    ! the number of cells and faces per partition (== processor)
    ! respectively.

    nCell = 0
    nFace = 0

    ! Determine the number of cells and faces per partition.
    ! Note that part(i) is the processor id, which starts at 0.

    do i=1,nblocks
       j = part(i) + 1
       nCell(j) = nCell(j) + blocks(i)%nCell
       nFace(j) = nFace(j) + blocks(i)%nFace
    enddo

    ! Determine the desirable number of cells and faces per processor.

    nCellsEven = nCell(1)
    nFacesEven = nFace(1)

    do i=2,nProc
       nCellsEven = nCellsEven + nCell(i)
       nFacesEven = nFacesEven + nFace(i)
    enddo

    nCellsEven = nCellsEven/nProc
    nFacesEven = nFacesEven/nProc

    ! Determine the maximum value of nCell and nFace
    ! and substract the optimal value.

    nCellMax = abs(maxval(nCell) - nCellsEven)
    nFaceMax = abs(maxval(nFace) - nFacesEven)

    ! Determine the tolerance for the cells and faces.

    nCellTol = (ubvec(1) - one)*nCellsEven
    nFaceTol = (ubvec(2) - one)*nFacesEven

    ! Check whether the load balance values for the cells and faces
    ! are met.

    cellsBalanced = .true.
    facesBalanced = .true.

    if(nCellMax > nCellTol) cellsBalanced = .false.
    if(nFaceMax > nFaceTol) facesBalanced = .false.

    ! Determine the load imbalances for the cells and faces
    ! and store it in ubvec.

    ubvec(1) = real(nCellMax,realType) &
         / real(nCellsEven,realType)
    ubvec(2) = real(nFaceMax,realType) &
         / real(nFacesEven,realType)

  end subroutine checkLoadBalance

  subroutine splitBlock(compBlock, nSub, nCells, ranges)
    !
    !       splitBlock tries to split the given computational block into
    !       the desired number of subblocks nSub. However it can happen
    !       that nSub is a strange number and a different splitting is
    !       performed. On return, nSub contains the actual number into the
    !       block is split. This number is smaller or equal to nSub on
    !       entry. As it is possible that the computational block itself
    !       is a subblock of an original cgns block, on return ranges will
    !       contain the nodal ranges of the subblocks in the original cgns
    !       block. The splitting attempts to keep the needed multigrid
    !       capabilities as much as possible.
    !       A recursive bisection algorithm is used.
    !
    use constants
    use inputIteration, only : smoother
    use partitionMod, onlY : distributionBlockType
    implicit none
    !
    !      Subroutine arguments.
    !
    type(distributionBlockType), intent(in)    :: compBlock
    integer(kind=intType),       intent(in)    :: nCells
    integer(kind=intType),       intent(inout) :: nSub

    integer(kind=intType), dimension(nSub,3,2), intent(out) :: ranges
    !
    !      Local variables.
    !
    integer(kind=intType) :: nLevels, level, nn, mm, nTarget
    integer(kind=intType) :: nSplit, nSplitNew

    integer(kind=intType), dimension(nSub) :: nSubblocks

    logical, dimension(3) :: viscousDir

    ! Determine the viscous directions of the block.

    viscousDir = .false.
    do nn=1,compBlock%nBocos

       ! Check for viscous boundary conditions and if found set the
       ! corresponding direction to viscous.

       if(compBlock%BCType(nn) == NSWallAdiabatic .or. &
            compBlock%BCType(nn) == NSWallIsothermal) then

          select case (compBlock%BCFaceID(nn))
          case (iMin,iMax)
             viscousDir(1) = .true.

          case (jMin,jMax)
             viscousDir(2) = .true.

          case (kMin,kMax)
             viscousDir(3) = .true.
          end select
       endif
    enddo

    ! Set viscousDir to .false. if an explicit smoother is used.

    if(smoother == RungeKutta) viscousDir = .false.

    ! Determine the number of levels in the bisection.

    nLevels = log(real(nSub,realType))/log(two)
    if(2**nLevels < nSub) nLevels = nLevels + 1

    ! Initialize the range of the first splittable block to the
    ! entire block.

    ranges(1,1,1) = 1; ranges(1,1,2) = compBlock%il
    ranges(1,2,1) = 1; ranges(1,2,2) = compBlock%jl
    ranges(1,3,1) = 1; ranges(1,3,2) = compBlock%kl

    ! Initialize the number of blocks to be split to 1 and the
    ! number of blocks it should be split into to nSub.

    nSplit = 1
    nSubblocks(1) = nSub

    ! Loop over the number of levels

    levelLoop: do level=1,nLevels

       ! Initialize nSplitNew to nSplit; it will be used to store the
       ! position of the new subblocks. At the end of this loop this
       ! value will be set to nSplit for the new level.

       nSplitNew = nSplit

       ! Loop over the number of blocks to be split.

       blockSplitLoop: do nn=1,nSplit

          ! Determine the situation we are having it.

          select case (nSubblocks(nn))

          case (2_intType,3_intType)

             ! Subblock must be split in either two or three. Store
             ! the number into the new subblocks must be in nSubblocks.
             ! As the routine split2block puts the subblock with the
             ! target number of cells in position nSplitNew, this
             ! block should not be split any further and thus gets a
             ! value of 1. Position nn gets the old value -1, which
             ! will be either 1 or 2

             nSplitNew             = nSplitNew + 1
             nSubblocks(nSplitNew) = 1
             nSubblocks(nn)        = nSubblocks(nn) - 1

             ! Split the block into two, where a block with nCells
             ! should be split off.

             call split2block(nSub, nn, nSplitNew, nCells, ranges, &
                  viscousDir)

             !===========================================================

         case (4_intType:)


             ! Subblock must be split in 4 or more. First determine the
             ! number of subblocks into the two new subblocks must be
             ! split in the next round. Save the old value of
             ! nSubblocks(nn) in mm.

             mm = nSubblocks(nn)

             nSplitNew             = nSplitNew + 1
             nSubblocks(nSplitNew) = nSubblocks(nn)/2
             nSubblocks(nn)        = nSubblocks(nn) &
                  - nSubblocks(nSplitNew)

             ! Determine the number of cells for the new subblock
             ! nSplitNew.

             nTarget = (ranges(nn,1,2) - ranges(nn,1,1)) &
                  * (ranges(nn,2,2) - ranges(nn,2,1)) &
                  * (ranges(nn,3,2) - ranges(nn,3,1))
             nTarget = nTarget*(real(nSubblocks(nSplitNew),realType) &
                  /          real(mm,realType))

             ! Split the block into two, where one block should get
             ! approximately nTarget cells.

             call split2block(nSub, nn, nSplitNew, nTarget, ranges, &
                  viscousDir)

          end select

       enddo blockSplitLoop

       ! Set nSplit to nSplitNew for the next level of bisection.

       nSplit = nSplitNew

    enddo levelLoop

    ! Get rid of the possible zero cell ranges. In that case
    ! nSub will change as well. Add the offset of the computational
    ! block, as this might be a subblock itself.

    mm   = nSub
    nSub = 0
    do nn=1,mm

       ! Test whether this is a non-zero cell subblock.

       if(ranges(nn,1,2) > ranges(nn,1,1) .and. &
            ranges(nn,2,2) > ranges(nn,2,1) .and. &
            ranges(nn,3,2) > ranges(nn,3,1)) then

          ! This is a valid subblock. Update nSub and copy the range
          ! including the offset of the computational block.

          nSub = nSub + 1

          ranges(nSub,1,1) = ranges(nn,1,1) + compBlock%iBegor - 1
          ranges(nSub,2,1) = ranges(nn,2,1) + compBlock%jBegor - 1
          ranges(nSub,3,1) = ranges(nn,3,1) + compBlock%kBegor - 1

          ranges(nSub,1,2) = ranges(nn,1,2) + compBlock%iBegor - 1
          ranges(nSub,2,2) = ranges(nn,2,2) + compBlock%jBegor - 1
          ranges(nSub,3,2) = ranges(nn,3,2) + compBlock%kBegor - 1

       endif
    enddo

  end subroutine splitBlock

  !========================================================================

  subroutine split2block(nSub, n1, n2, nTarget, ranges, &
       viscousDir)
    !
    !       split2block splits the block stored in ranges(n1,:,:) into
    !       two. The new blocks are stored in ranges(n1,:,:) and
    !       ranges(n2,:,:), where the n2 will store the block with the
    !      * number of cells closest to nTarget.
    !
    use inputIteration
    implicit none
    !
    !      Subroutine arguments.
    !
    integer(kind=intType), intent(in) :: nSub, n1, n2, nTarget

    integer(kind=intType), dimension(nSub,3,2), intent(inout) :: ranges

    logical, dimension(3), intent(in) :: viscousDir
    !
    !      Local variables.
    !
    integer(kind=intType) :: nMG, mm, ii, jj, i, iiOpt, iOpt
    integer(kind=intType) :: nCellOpt

    integer(kind=intType), dimension(3) :: nc, mc, c, nf, prefDir

    ! Determine the maximum of the number of mg levels needed in the
    ! cycle and the mg start level. Store this value in nMG.

    nMG = max(nMGLevels, mgStartlevel)

    ! Store the number of cells of the original subblock in nc.

    nc(1) = ranges(n1,1,2) - ranges(n1,1,1)
    nc(2) = ranges(n1,2,2) - ranges(n1,2,1)
    nc(3) = ranges(n1,3,2) - ranges(n1,3,1)

    ! Determine its multigrid capabilities in the three directions.

    loopDir: do mm=1,3

       ! Initialize nf(mm), number of mg levels, to 1 and ii to 2.
       ! nf is used as a temporary buffer.

       nf(mm) = 1
       ii     = 2

       ! Loop to determine the number of mg levels.

       do
          ! Conditions to exit the loop. The first condition is the
          ! truncation to the maximum number of mg levels needed by the
          ! iterative algorithm of the solver. The second condition is
          ! simply that the number of cells do not allow more standard
          ! mg levels.

          if(nf(mm) == nMG) exit
          if(mod(nc(mm),ii) /= 0) exit

          ! Update nf(mm) and multiply ii by 2 to test for the
          ! next level.

          nf(mm) = nf(mm) +1
          ii     = 2*ii
       enddo

       ! Determine the constants c(mm) and mc(mm), such that
       ! nc(mm) = c(mm)*mc(mm). Mc(mm) is the number of cells in every
       ! direction per supercell. A supercell is the smallest possible
       ! block able to do multigridding, i.e. mc(mm) = 2**(nf(mm)-1)
       ! if the block allows for nMG multigrid levels.

       c(mm)  = 2*nc(mm)/ii
       mc(mm) = nc(mm)/c(mm)

    enddo loopDir

    ! Test whether the multigrid requirements are not too restrictive
    ! for this subblock. If they are, relax them, if possible.

    if(nMG > 1 .and. max(c(1),c(2),c(3)) == 1) then
       c(1) =  c(1)*2;  c(2) =  c(2)*2;  c(3) =  c(3)*2
       mc(1) = mc(1)/2; mc(2) = mc(2)/2; mc(3) = mc(3)/2
    endif

    ! Determine the number of faces on a slice in the
    ! three directions.

    nf(1) = nc(2)*nc(3)
    nf(2) = nc(1)*nc(3)
    nf(3) = nc(1)*nc(2)

    ! Determine the preferred split direction; this the direction
    ! which results in the least amount of additional faces. This is
    ! important, because the number of flux evaluations is
    ! proportional to the number of faces.

    if(nf(1) < nf(2)) then
       if(nf(1) < nf(3)) then
          prefDir(1) = 1
          if(nf(2) < nf(3)) then
             prefDir(2) = 2
             prefDir(3) = 3
          else
             prefDir(2) = 3
             prefDir(3) = 2
          endif
       else
          prefDir(1) = 3
          prefDir(2) = 1
          prefDir(3) = 2
       endif
    else if(nf(2) < nf(3)) then
       prefDir(1) = 2
       if(nf(1) < nf(3)) then
          prefDir(2) = 1
          prefDir(3) = 3
       else
          prefDir(2) = 3
          prefDir(3) = 1
       endif
    else
       prefDir(1) = 3
       prefDir(2) = 2
       prefDir(3) = 1
    endif

    ! Take viscousDir into account to determine the preference
    ! direction. For implicit methods it may not be a good idea to
    ! split the block parallel to a viscous wall.

    ! Check the 1st preference direction. If it is viscous swap
    ! it with the best inviscid direction, if available.

    if( viscousDir(prefDir(1)) ) then
       if(.not. viscousDir(prefDir(2)) ) then
          mm         = prefDir(1)
          prefDir(1) = prefDir(2)
          prefDir(2) = mm
       else if(.not. viscousDir(prefDir(3)) ) then
          mm         = prefDir(1)
          prefDir(1) = prefDir(3)
          prefDir(3) = mm
       endif
    endif

    ! Check the second preference direction. If it is viscous
    ! try to swap it whith the third direction.

    if(viscousDir(prefDir(2)) .and. &
         .not. viscousDir(prefDir(3)) ) then
       mm         = prefDir(2)
       prefDir(2) = prefDir(3)
       prefDir(3) = mm
    endif

    ! Determine the splitting which is best. Initialize nCellOpt
    ! to a ridiculously high number.

    nCellOpt = 10*nc(1)*nc(2)*nc(3)

    ! Loop over the three directions.

    do mm=1,3

       ! Store the current direction considered in ii and determine
       ! the index i, which gives a number of cells as close as
       ! possible to the desired number.

       ii = prefDir(mm)
       i  = nint(real(nTarget,realType) &
            /      real(mc(ii)*nf(ii),realType),intType)
       i  = max(i,1_intType)

       ! Determine whether the corresponding number of cells is
       ! closer to the target number than the currently stored value.
       ! If so, store the settings for this splitting.

       jj = i*mc(ii)*nf(ii)
       if(abs(jj - nTarget) < abs(nCellOpt - nTarget)) then
          nCellOpt = jj
          iiOpt    = ii
          iOpt     = i
       endif

    enddo

    ! Copy the range of subblock n1 into n2 as initialization.

    ranges(n2,1,1) = ranges(n1,1,1); ranges(n2,1,2) = ranges(n1,1,2)
    ranges(n2,2,1) = ranges(n1,2,1); ranges(n2,2,2) = ranges(n1,2,2)
    ranges(n2,3,1) = ranges(n1,3,1); ranges(n2,3,2) = ranges(n1,3,2)

    ! Determine the nodal index where the splitting takes place.

    jj = iOpt*mc(iiOpt) + ranges(n1,iiOpt,1)

    ! Adapt the corresponding indices of the new subblocks n1 and
    ! n2, such that it corresponds to the new situation; n2 should
    ! contain the subblock, which contains a number of cells as
    ! close as possible to nTarget.

    ranges(n2,iiOpt,2) = jj
    ranges(n1,iiOpt,1) = jj

  end subroutine split2block

  subroutine reallocSubfaceMemory(ii,nAlloc)
    !
    !       reallocSubfaceMemory reallocates the memory to store the
    !       subface information for the given block ii. On entry nAlloc
    !       contains the current number of allocated subfaces, on exit
    !       this is updated to the new number.
    !
    use constants
    use partitionmod, only : blocks
    use utils, only : reallocateInteger
    implicit none
    !
    !      Subroutine arguments.
    !
    integer(kind=intType), intent(in)    :: ii
    integer(kind=intType), intent(inout) :: nAlloc
    !
    !      Local arguments.
    !
    integer(kind=intType) :: nOld, nNew
    !
    ! Store the old value of the allocated array and determine the
    ! new one. Store this new value in nAlloc.

    nOld   = nAlloc
    nNew   = nOld + 6
    nAlloc = nNew

    ! Reallocate the memory.

    call reallocateInteger(blocks(ii)%BCType,      nNew, nOld, .true.)
    call reallocateInteger(blocks(ii)%BCFaceID,    nNew, nOld, .true.)
    call reallocateInteger(blocks(ii)%cgnsSubface, nNew, nOld, .true.)

    call reallocateInteger(blocks(ii)%inBeg, nNew, nOld, .true.)
    call reallocateInteger(blocks(ii)%jnBeg, nNew, nOld, .true.)
    call reallocateInteger(blocks(ii)%knBeg, nNew, nOld, .true.)

    call reallocateInteger(blocks(ii)%inEnd, nNew, nOld, .true.)
    call reallocateInteger(blocks(ii)%jnEnd, nNew, nOld, .true.)
    call reallocateInteger(blocks(ii)%knEnd, nNew, nOld, .true.)

    call reallocateInteger(blocks(ii)%dinBeg, nNew, nOld, .true.)
    call reallocateInteger(blocks(ii)%djnBeg, nNew, nOld, .true.)
    call reallocateInteger(blocks(ii)%dknBeg, nNew, nOld, .true.)

    call reallocateInteger(blocks(ii)%dinEnd, nNew, nOld, .true.)
    call reallocateInteger(blocks(ii)%djnEnd, nNew, nOld, .true.)
    call reallocateInteger(blocks(ii)%dknEnd, nNew, nOld, .true.)

    call reallocateInteger(blocks(ii)%neighBlock, nNew, nOld, .true.)
    call reallocateInteger(blocks(ii)%groupNum,   nNew, nOld, .true.)

    call reallocateInteger(blocks(ii)%l1, nNew, nOld, .true.)
    call reallocateInteger(blocks(ii)%l2, nNew, nOld, .true.)
    call reallocateInteger(blocks(ii)%l3, nNew, nOld, .true.)

  end subroutine reallocSubfaceMemory

end module loadBalance

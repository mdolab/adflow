!
!      ******************************************************************
!      *                                                                *
!      * File:          blockDistribution.f90                           *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 03-13-2004                                      *
!      * Last modified: 10-10-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine blockDistribution
!
!      ******************************************************************
!      *                                                                *
!      * blockDistribution determines the distribution of the blocks    *
!      * over the processors. If blocks must be split to obtain a good  *
!      * load balance an iterative algorithm is used to determine the   *
!      * best way to split them.                                        *
!      *                                                                *
!      ******************************************************************
!
       use cgnsGrid
       use communication
       use inputParallel
       use partitionMod
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
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
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
         call mpi_barrier(SUmb_comm_world, ierr)
       endif
!
!      ******************************************************************
!      *                                                                *
!      * Determine how the blocks must be split (if allowed) for        *
!      * load balancing reasons.                                        *
!      *                                                                *
!      ******************************************************************
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
       if( splitBlocks ) iterMax = 2

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
         call mpi_barrier(SUmb_comm_world, ierr)
       endif

       ! Processor 0 prints a warning if the communication was
       ! neglected to obtain a valid partitioning.

       if(commNeglected .and. myID == 0) then
         print "(a)", "#"
         print "(a)", "#                    Warning"
         print "(a)", "# Communication costs neglected to obtain&
                       & a valid partitioning."
         print "(a)", "#"
       endif

       ! If the load imbalance tolerance was not met, print a warning
       ! message if I am processor 0.

       if(.not.(cellsBalanced .and. facesBalanced) .and. &
          myID == 0) then
         print "(a)", "#"
         print "(a)", "#                    Warning"
         print 100, loadImbalance
         print 101, ubvec(1), ubvec(2)
         print "(a)", "#"
 100     format("# Specified load imbalance tolerance",1X,F6.3,1X,"not &
                &achieved.")
 101     format("# I continue with",1X,F6.3,1X,"load imbalance for the &
                &cells and",1X,F6.3,1X,"for the faces")
       endif

       !=================================================================

       contains

         !===============================================================

         logical function splittingIsOkay(cgnsID)
!
!        ****************************************************************
!        *                                                              *
!        * splittingIsOkay determines whether or not the splitting of   *
!        * the given cgns block is okay in the sense that all subblocks *
!        * are smaller than the allowed number of cells and faces.      *
!        *                                                              *
!        ****************************************************************
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
!
!        ****************************************************************
!        *                                                              *
!        * Begin execution                                              *
!        *                                                              *
!        ****************************************************************
!
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
!        ****************************************************************
!        *                                                              *
!        * splitBlockInitialization splits the given cgns block ID      *
!        * into a number of subbocks during the initialization phase.   *
!        *                                                              *
!        ****************************************************************
!
         use BCTypes
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
!
!        ****************************************************************
!        *                                                              *
!        * Begin execution                                              *
!        *                                                              *
!        ****************************************************************
!
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
!        ****************************************************************
!        *                                                              *
!        * splitBlocksLoadBalance splits some (sub)blocks even          *
!        * further to obtain a better load balance.                     *
!        *                                                              *
!        ****************************************************************
!
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
!
!        Function definition
!
         integer(kind=intType) :: bsearchIntegers
!
!        ****************************************************************
!        *                                                              *
!        * Begin execution                                              *
!        *                                                              *
!        ****************************************************************
!
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
           j = bsearchIntegers(nCellOr(i), nCell, nCellDiff)
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

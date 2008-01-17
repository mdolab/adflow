!
!      ******************************************************************
!      *                                                                *
!      * File:          determineComputeBlocks.f90                      *
!      * Author:        Edwin van der Weide, Steve Repsher              *
!      * Starting date: 03-19-2004                                      *
!      * Last modified: 10-10-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine determineComputeBlocks(splitInfo)
!
!      ******************************************************************
!      *                                                                *
!      * determineComputeBlocks determines the computational blocks     *
!      * from the original grid and the given information how to split  *
!      * these blocks.                                                  *
!      *                                                                *
!      ******************************************************************
!
       use cgnsGrid
       use partitionMod
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
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
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
                      blocks(i)%cgnsOver,    blocks(i)%ipntOver,   &
                      blocks(i)%neighOver,   blocks(i)%overComm,   &
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
                    blocks(ii)%groupNum(nAlloc),    stat=ierr)
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

       ! Loop over the cgns blocks one more time to distribute all of
       ! its overset cells to the sublocks.

       do i=1,cgnsNDom
         call distributeOversetCells(i, nsubPerCGNS, splitInfo)
       end do

       end subroutine determineComputeBlocks

!========================================================================

       subroutine BCFacesSubblock(cgnsID, ii, jj)
!
!      ******************************************************************
!      *                                                                *
!      * BCFacesSubblock determines the boundary subfaces of compute    *
!      * block ii, which is a subblock of the given cgns block.         *
!      * Jj is the counter for the number of subfaces.                  *
!      *                                                                *
!      ******************************************************************
!
       use BCTypes
       use cgnsGrid
       use communication
       use inputPhysics
       use partitionMod
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
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
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
 100           format("Zone",1X,A,", boundary subface",1X,A, &
                      ": No constant index found for subface")
               call terminate("BCFacesSubblock", errorMessage)
             endif

             call mpi_barrier(SUmb_comm_world, ierr)
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
 120             format("Zone",1X,A,", boundary subface",1X,A, &
                        ": Not a valid boundary condition for &
                        &internal flow")
               else
                 write(errorMessage,130) trim(zoneName), trim(subName)
 130             format("Zone",1X,A,", boundary subface",1X,A, &
                        ": Not a valid boundary condition for &
                        &external flow")
               endif

               call terminate("BCFacesSubblock", errorMessage)
             endif

             call mpi_barrier(SUmb_comm_world, ierr)
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
!      ******************************************************************
!      *                                                                *
!      * externalFacesSubblock determines the block boundaries of       *
!      * the compute block ii which are located on the boundaries of    *
!      * the given original cgns block. As it is possible that due to   *
!      * a splitting of a neighboring block the number of block         *
!      * boundaries is larger than the original number, it must be      *
!      * checked whether enough memory has been allocated.              *
!      * jj is the counter for the number of subfaces.                  *
!      *                                                                *
!      ******************************************************************
!
       use BCTypes
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

       integer(kind=intType), dimension(3,3) :: tMat

       integer(kind=intType), dimension(:,:,:), pointer :: ranges

       character(len=maxCGNSNameLen) :: zoneName, subName
       character(len=2*maxStringLen) :: errorMessage

       logical :: diSwap, djSwap, dkSwap
!
!      Function definitions.
!
       integer(kind=intType) :: delta
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
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
 140           format("Zone",1X,A,", 1 to 1 block connectivity",1X,A, &
                      ": No constant index found for subface")
               call terminate("externalFacesSubblock", errorMessage)
             endif

             call mpi_barrier(SUmb_comm_world, ierr)
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

           tMat(1,1) = sign(1_intType,l1) * delta(l1,1_intType)
           tMat(2,1) = sign(1_intType,l1) * delta(l1,2_intType)
           tMat(3,1) = sign(1_intType,l1) * delta(l1,3_intType)

           tMat(1,2) = sign(1_intType,l2) * delta(l2,1_intType)
           tMat(2,2) = sign(1_intType,l2) * delta(l2,2_intType)
           tMat(3,2) = sign(1_intType,l2) * delta(l2,3_intType)

           tMat(1,3) = sign(1_intType,l3) * delta(l3,1_intType)
           tMat(2,3) = sign(1_intType,l3) * delta(l3,2_intType)
           tMat(3,3) = sign(1_intType,l3) * delta(l3,3_intType)

           ! Determine the corresponding donor range of the subface
           ! iBeg, iEnd; jBeg, jEnd; kBeg, kEnd.

           l1 = iBeg - cgnsDoms(cgnsID)%conn1to1(j)%iBeg
           L2 = jBeg - cgnsDoms(cgnsID)%conn1to1(j)%jBeg
           l3 = kBeg - cgnsDoms(cgnsID)%conn1to1(j)%kBeg

           diBeg = cgnsDoms(cgnsID)%conn1to1(j)%diBeg &
                 + tMat(1,1)*l1 + tMat(1,2)*l2 + tMat(1,3)*l3
           djBeg = cgnsDoms(cgnsID)%conn1to1(j)%djBeg &
                 + tMat(2,1)*l1 + tMat(2,2)*l2 + tMat(2,3)*l3
           dkBeg = cgnsDoms(cgnsID)%conn1to1(j)%dkBeg &
                 + tMat(3,1)*l1 + tMat(3,2)*l2 + tMat(3,3)*l3

           l1 = iEnd - cgnsDoms(cgnsID)%conn1to1(j)%iBeg
           L2 = jEnd - cgnsDoms(cgnsID)%conn1to1(j)%jBeg
           l3 = kEnd - cgnsDoms(cgnsID)%conn1to1(j)%kBeg

           diEnd = cgnsDoms(cgnsID)%conn1to1(j)%diBeg &
                 + tMat(1,1)*l1 + tMat(1,2)*l2 + tMat(1,3)*l3
           djEnd = cgnsDoms(cgnsID)%conn1to1(j)%djBeg &
                 + tMat(2,1)*l1 + tMat(2,2)*l2 + tMat(2,3)*l3
           dkEnd = cgnsDoms(cgnsID)%conn1to1(j)%dkBeg &
                 + tMat(3,1)*l1 + tMat(3,2)*l2 + tMat(3,3)*l3

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
               ! The inverse of the transformation matrix tMat is
               ! the transpose.

               l1 = iBeg - cgnsDoms(cgnsID)%conn1to1(j)%diBeg
               L2 = jBeg - cgnsDoms(cgnsID)%conn1to1(j)%djBeg
               l3 = kBeg - cgnsDoms(cgnsID)%conn1to1(j)%dkBeg

               iBeg = cgnsDoms(cgnsID)%conn1to1(j)%iBeg &
                    + tMat(1,1)*l1 + tMat(2,1)*l2 + tMat(3,1)*l3
               jBeg = cgnsDoms(cgnsID)%conn1to1(j)%jBeg &
                    + tMat(1,2)*l1 + tMat(2,2)*l2 + tMat(3,2)*l3
               kBeg = cgnsDoms(cgnsID)%conn1to1(j)%kBeg &
                    + tMat(1,3)*l1 + tMat(2,3)*l2 + tMat(3,3)*l3

               l1 = iEnd - cgnsDoms(cgnsID)%conn1to1(j)%diBeg
               L2 = jEnd - cgnsDoms(cgnsID)%conn1to1(j)%djBeg
               l3 = kEnd - cgnsDoms(cgnsID)%conn1to1(j)%dkBeg

               iEnd = cgnsDoms(cgnsID)%conn1to1(j)%iBeg &
                    + tMat(1,1)*l1 + tMat(2,1)*l2 + tMat(3,1)*l3
               jEnd = cgnsDoms(cgnsID)%conn1to1(j)%jBeg &
                    + tMat(1,2)*l1 + tMat(2,2)*l2 + tMat(3,2)*l3
               kEnd = cgnsDoms(cgnsID)%conn1to1(j)%kBeg &
                    + tMat(1,3)*l1 + tMat(2,3)*l2 + tMat(3,3)*l3

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
!      ******************************************************************
!      *                                                                *
!      * internalFacesSubblock determines the block boundaries of       *
!      * the compute block ii which are created due to the splitting of *
!      * the original block into subblock. As the number of these       *
!      * internal boundaries is not known, it must be checked whether   *
!      * enough memory has been allocated. jj is the counter for the    *
!      * number of subfaces.                                            *
!      *                                                                *
!      ******************************************************************
!
       use BCTypes
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
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
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
!        ****************************************************************
!        *                                                              *
!        * searchInternalNeighbors determines block faces created by    *
!        * by the splitting of the original block. The variables set in *
!        * internalFacesSubblock are used such that a general           *
!        * treatment is possible.                                       *
!        *                                                              *
!        ****************************************************************
!
         implicit none
!
!        Local variables
!
         integer(kind=intType) :: mm, jnBeg, jnEnd, knBeg, knEnd

         integer(kind=intType), dimension(3,2) :: subRange

         integer(kind=intType), dimension(:,:,:), pointer :: ranges
!
!        ****************************************************************
!        *                                                              *
!        * Begin execution                                              *
!        *                                                              *
!        ****************************************************************
!
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

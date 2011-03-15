!
!      ******************************************************************
!      *                                                                *
!      * File:          determineVarLocalIO.f90                         *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 10-27-2005                                      *
!      * Last modified: 09-13-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine determineVarLocalIO(nGlobalBlocks,  nChunks,      &
                                      sizeP3D_Real,   myChunk,      &
                                      globalBlockIDs, startBlock,   &
                                      blockDimCGNS,   cgnsOffset,   &
                                      chunks,         dataIsRead,   &
                                      includeHalos)
!
!      ******************************************************************
!      *                                                                *
!      * determineVarLocalIO determines the variables needed to copy    *
!      * data between the IO buffer and the local IOVar. These          *
!      * variables are the same for both reading and writing.           *
!      *                                                                *
!      ******************************************************************
!
       use cgnsGrid
       use communication
       use IOModule
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType), intent(in) :: nGlobalBlocks,  nChunks
       integer(kind=intType), intent(in) :: sizeP3D_Real,   myChunk

       integer(kind=intType), dimension(nGlobalBlocks), intent(in) :: &
                                                          globalBlockIDs
       integer(kind=intType), dimension(3,nGlobalBlocks), &
                                              intent(in) :: blockDimCGNS

       integer(kind=intType), dimension(nChunks), intent(in) :: &
                                                              startBlock

       integer(kind=mpi_offset_kind), dimension(0:nGlobalBlocks), &
                                               intent(in) :: cgnsOffset
       integer(kind=mpi_offset_kind), dimension(0:nChunks), &
                                               intent(in) :: chunks

       logical, intent(in) :: dataIsRead, includeHalos
!
!      Local variables.
!
       integer :: ierr

       integer(kind=intType) :: ii, jj, kk, ll, mm, nn, i
       integer(kind=intType) :: offset, iBeg, jBeg, kBeg, lBeg
       integer(kind=intType) :: iGMin, jGMin, kGMin, iGMax, jGMax, kGMax
       integer(kind=intType) :: iMin, jMin, kMin, iMax, jMax, kMax
       integer(kind=intType) :: iOffset, jOffset, kOffset
       integer(kind=intType) :: ic, jc, kc, lc

       integer(kind=mpi_offset_kind) :: lowerBound
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
!      ******************************************************************
!      *                                                                *
!      * Determine the variables used for the actual local IO. These    *
!      * are stored in the module IOModule.                             *
!      *                                                                *
!      ******************************************************************
!
       ! Loop over the number of global (sub)blocks in my chunk to
       ! determine the local data to be copied.

       globalBlockLoop: do nn=1,P3D_nIOParts

         ! Store in mm the entry globalBlockIDs and in jj the
         ! corresponding CGNS block ID..

         mm = startBlock(myChunk) + nn - 1
         jj = globalBlockIDs(mm)

         ! Set the index range for this global block.

         iGMin = 1; iGMax = cgnsDoms(jj)%il
         jGMin = 1; jGMax = cgnsDoms(jj)%jl
         kGMin = 1; kGMax = cgnsDoms(jj)%kl

         if(P3D_DataStorage == cellDataNoHalo) then
           iGMin = 2; jGMin = 2; kGMin = 2;
         else if(P3D_DataStorage == cellDataPlusHalo) then
           iGMax = iGMax + 1; jGMax = jGMax + 1; kGMax = kGMax + 1
         endif

         ! Determine the lower bound for this global block part inside
         ! this chunk. The record integer is taken into account,
         ! because it does contain information of interest.

         lowerBound = cgnsOffset(mm-1) + nBytesPerRecordIntPLOT3D
         lowerBound = max(lowerBound, chunks(myChunk-1))

         ! Determine the indices iBeg, jBeg, kBeg and lBeg to which
         ! the lower bound corresponds.

         offset = lowerBound - cgnsOffset(mm-1) &
                - nBytesPerRecordIntPLOT3D

         call getBegIndices(offset, sizeP3D_Real, blockDimCGNS(1,mm), &
                            iBeg, jBeg, kBeg, lBeg)

         ! Allocate the memory for the data of the local copy, i.e. 
         ! blockID, indices, indexW and posLocal.

         ii = P3D_IOParts(nn)%nItemsLocal
         allocate(P3D_IOParts(nn)%blockID(ii),  &
                  P3D_IOParts(nn)%indexW(ii),   &
                  P3D_IOParts(nn)%posLocal(ii), &
                  P3D_IOParts(nn)%indices(ii,3), stat=ierr)
         if(ierr /= 0)                           &
           call terminate("determineVarLocalIO", &
                          "Memory allocation failure for blockID, etc.")
!
!        ****************************************************************
!        *                                                              *
!        * Loop over the local subblocks of global block jj and         *
!        * determine the mapping between the IO buffer and the indices  *
!        * of IOVar.                                                    *
!        *                                                              *
!        ****************************************************************
!
         ! Loop over the local subblocks of the global block.
         ! kk is the counter for the local entries.

         kk = 0

         localSubBlockLoop: do ll=1,cgnsDoms(jj)%nSubBlocks

           ! Check if this subblock is stored locally.

           testLocal: if(cgnsDoms(jj)%procStored(ll) == myID) then

             ! Store the local block ID a bit easier in ii.

             ii = cgnsDoms(jj)%localBlockID(ll)

             ! Determine the range of this subblock in the original block.
             ! This depends on what type of data is read/written.

             call getMinMaxIndices(cgnsDoms(jj)%iBegOr(ll), &
                                   cgnsDoms(jj)%jBegOr(ll), &
                                   cgnsDoms(jj)%kBegOr(ll), &
                                   cgnsDoms(jj)%iEndOr(ll), &
                                   cgnsDoms(jj)%jEndOr(ll), &
                                   cgnsDoms(jj)%kEndOr(ll), &
                                   cgnsDoms(jj)%il,         &
                                   cgnsDoms(jj)%jl,         &
                                   cgnsDoms(jj)%kl,         &
                                   dataIsRead,              &
                                   includeHalos,            &
                                   iMin, jMin, kMin, iMax, jMax, kMax)

             ! Determine the offset between global and local indices.

             iOffset = cgnsDoms(jj)%iBegOr(ll) - 1
             jOffset = cgnsDoms(jj)%jBegOr(ll) - 1
             kOffset = cgnsDoms(jj)%kBegOr(ll) - 1

             ! Loop over the total number of entities in this IO part
             ! and store the indices which coincide.

             ic = iBeg; jc = jBeg; kc = kBeg; lc = lBeg

             itemsLocalLoop: do i=1,P3D_IOParts(nn)%nItemsTotal

               ! Check if the current i, j, k indices coincide with the
               ! subblock considered.

               if(ic >= iMin .and. jc >= jMin .and. kc >= kMin .and. &
                  ic <= iMax .and. jc <= jMax .and. kc <= kMax) then

                 ! Entity is part of the subblock. Store it in
                 ! P3D_IOParts(nn).

                 kk = kk + 1

                 P3D_IOParts(nn)%blockID(kk) = cgnsDoms(jj)%localBlockID(ll)
                 P3D_IOParts(nn)%indexW(kk)  = lc - 1

                 P3D_IOParts(nn)%indices(kk,1) = ic - iOffset
                 P3D_IOParts(nn)%indices(kk,2) = jc - jOffset
                 P3D_IOParts(nn)%indices(kk,3) = kc - kOffset

                 P3D_IOParts(nn)%posLocal(kk) = i
               endif

               ! Update the working indices ic, jc, kc and lc.

               ic = ic + 1
               if(ic > iGMax) then
                 ic = iGMin; jc = jc + 1
               endif

               if(jc > jGMax) then
                 jc = jGMin; kc = kc + 1
               endif

               if(kc > kGMax) then
                 kc = kGMin; lc = lc + 1
               endif

             enddo itemsLocalLoop
           endif testLocal
         enddo localSubBlockLoop

         ! Test in debug mode that the desired number of data is
         ! stored in P3D_IOParts(nn).

         if( debug ) then
           if(kk /= P3D_IOParts(nn)%nItemsLocal)   &
             call terminate("determineVarLocalIO", &
                            "kk differs from nItemsLocal")
         endif

       enddo globalBlockLoop

       end subroutine determineVarLocalIO

!
!      ******************************************************************
!      *                                                                *
!      * File:          determineBytesCommChunk.f90                     *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 10-26-2005                                      *
!      * Last modified: 09-13-2007                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine determineBytesCommChunk(nGlobalBlocks,    &
                                          sizeP3D_Real,     &
                                          nChunks,          &
                                          globalBlockIDs,   &
                                          nPartsPerChunk,   &
                                          startBlock,       &
                                          nBytesCommChunk,  &
                                          blockDimCGNS,     &
                                          cgnsOffset,       &
                                          chunks,           &
                                          dataIsRead,       &
                                          includeHalos,     &
                                          IParticipateInIO)
!
!      ******************************************************************
!      *                                                                *
!      * determineBytesCommChunk determines the number of bytes of each *
!      * chunk to be read/written that would be communicated to other   *
!      * processors.                                                    *
!      *                                                                *
!      ******************************************************************
!
       use block
       use cgnsGrid
       use communication
       use IOModule
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType), intent(in)  :: nGlobalBlocks
       integer(kind=intType), intent(in)  :: sizeP3D_Real
       integer(kind=intType), intent(out) :: nChunks

       integer(kind=intType), dimension(nGlobalBlocks), intent(in) :: &
                                                           globalBlockIDs
       integer(kind=intType), dimension(3,nGlobalBlocks), &
                                              intent(out) :: blockDimCGNS

       integer(kind=intType), dimension(nProc), intent(out) :: nPartsPerChunk
       integer(kind=intType), dimension(nProc), intent(out) :: startBlock
       integer(kind=intType), dimension(nProc), intent(out) :: nBytesCommChunk

       integer(kind=mpi_offset_kind), dimension(0:nGlobalBlocks), &
                                               intent(out) :: cgnsOffset
       integer(kind=mpi_offset_kind), dimension(0:nProc), &
                                               intent(out) :: chunks

       logical, intent(in)  :: dataIsRead, includeHalos
       logical, intent(out) :: IParticipateInIO
!
!      Local variables.
!
       integer :: ierr

       integer(kind=intType) :: ii, jj, ll, mm, nn
       integer(kind=intType) :: sizePerEntity, sizeComm, offset
       integer(kind=intType) :: sizeRel, sizeRealPart
       integer(kind=intType) :: iBeg, jBeg, kBeg, lBeg
       integer(kind=intType) :: iEnd, jEnd, kEnd, lEnd
       integer(kind=intType) :: iiBeg, jjBeg, kkBeg, llBeg
       integer(kind=intType) :: iiEnd, jjEnd, kkEnd, llEnd
       integer(kind=intType) :: iMin, jMin, kMin, iMax, jMax, kMax

       integer(kind=mpi_offset_kind) :: sizeAvg, sizeTarget
       integer(kind=mpi_offset_kind) :: lowerBound, upperBound

       logical, dimension(cgnsNDom) :: blockIsWritten
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Determine the number of chunks and whether or not I participate
       ! at all in the IO.

       testAllBlocks: if(nGlobalBlocks == cgnsNDom) then

         ! Everything must be read/written. All processors participate
         ! and the number of chunks equals the number of processors.

         IParticipateInIO = .true.
         nChunks          = nProc

       else testAllBlocks

         ! Only a subset of the blocks must be read/written.
         ! Determine if I participate at all.

         blockIsWritten = .false.
         do nn=1,nGlobalBlocks
           blockIsWritten(globalBlockIDs(nn)) = .true.
         enddo

         IParticipateInIO = .false.
         do nn=1,nDom
           ii = flowDoms(nn,1,1)%cgnsBlockID
           if( blockIsWritten(ii) ) IParticipateInIO = .true.
         enddo

         ! Determine the number of processors which participate in the
         ! IO. This will be the number of chunks.

         ii = 0
         if( IParticipateInIO ) ii = 1

         call mpi_allreduce(ii, nChunks, 1, sumb_integer, mpi_sum, &
                            SUmb_comm_world, ierr)

       endif testAllBlocks
!
!      ******************************************************************
!      *                                                                *
!      * Determine the offset of the blocks to be written relative to   *
!      * the start of the first block.                                  *
!      *                                                                *
!      ******************************************************************
!
       ! Determine the total size in bytes per global block to be
       ! read/written. For historic reasons this is called a CGNS block
       ! and therefore the variable name is cgnsOffset. Cumulative
       ! storage format is used. The actual size depends on quite a few
       ! things, like the size of a floating point variable, what type
       ! of data is read and whether or not an iblanking variable is
       ! present.

       select case (P3D_DataStorage)
         case (cellDataNoHalo)
           ii = -1
         case (nodeData)
           ii =  0
         case (cellDataPlusHalo)
           ii =  1
       end select

       sizePerEntity = P3D_nVar*sizeP3D_Real
       if( P3D_iblank ) &
        sizePerEntity = sizePerEntity + nBytesPerIntPLOT3D

       cgnsOffset(0) = 0

       do nn=1,nGlobalBlocks
         jj = globalBlockIDs(nn)
         blockDimCGNS(1,nn)   = cgnsDoms(jj)%il + ii
         blockDimCGNS(2,nn)   = cgnsDoms(jj)%jl + ii
         blockDimCGNS(3,nn)   = cgnsDoms(jj)%kl + ii

         cgnsOffset(nn) = cgnsOffset(nn-1) + 2*nBytesPerRecordIntPLOT3D &
                        + sizePerEntity      * blockDimCGNS(1,nn)       &
                        * blockDimCGNS(2,nn) * blockDimCGNS(3,nn)
       enddo

       ! Determine the average amount of bytes per chunk and check if
       ! the amount is acceptable. Actually the check is made on an upper
       ! bound, because it is just to check if the 2 Gbyte per processor
       ! is not exceeded. In practice this is a ridiculous limit, but
       ! just to be sure.

       sizeAvg = cgnsOffset(cgnsNDom)/nChunks

       if((sizeAvg + 10*sizeP3D_Real) >= maxSizeIO) then
         if(myID == 0)                               &
           call terminate("determineBytesCommChunk", &
                          "Amount of data to read per processor is &
                          &too much. Increase the number of processors")
         call mpi_barrier(SUmb_comm_world, ierr)
       endif
!
!      ******************************************************************
!      *                                                                *
!      * Determine the byte offset for the chunks in cumulative storage *
!      * format.                                                        *
!      *                                                                *
!      ******************************************************************
!
       chunks(0)  = 0
       sizeTarget = 0
       mm = 1

       ! Loop over the nChunks - 1 inner boundaries that need to
       ! be determined.

       chunkLoop: do nn=1,nChunks-1

         ! Determine the ideal end for this chunk.

         sizeTarget = sizeTarget + sizeAvg
         chunks(nn) = sizeTarget

         ! Find out in which global block this end is situated.

         ii = mm
         do
           if(chunks(nn) <= cgnsOffset(mm)) exit
           mm = mm + 1
         enddo

         ! Make sure that chunks(nn) corresponds to the end of an entity.
         ! Several situations should be distinguished.

         if((cgnsOffset(mm) - chunks(nn)) <= &
            (nBytesPerRecordIntPLOT3D + 10*sizeP3D_Real)) then

           ! End of the chunk is pretty close to the end of the record
           ! for global block mm. Simply set it to the end of the block.

           chunks(nn) = cgnsOffset(mm)

         else if((chunks(nn) - cgnsOffset(mm-1)) <= &
                 (nBytesPerRecordIntPLOT3D + 10*sizeP3D_Real)) then

           ! End of the chunk is pretty close to the beginning of the
           ! record for global block mm. Simply set it to the end of
           ! the previous block. This means that mm must be decremented.

           mm = mm - 1
           chunks(nn) = cgnsOffset(mm)

         else

           ! The current chunk ends somewhere in the middle of the record
           ! of global block mm. Make sure that it ends appropriately,
           ! i.e. at the end of an entity.

           ! Determine the number of bytes relative to start of global
           ! block mm. Remove the integer of the record as well.
           ! Set chunks(nn) to the beginning of the data of the block.

           sizeRel    = chunks(nn) - cgnsOffset(mm-1) &
                      - nBytesPerRecordIntPLOT3D
           chunks(nn) = cgnsOffset(mm-1) + nBytesPerRecordIntPLOT3D

           ! Check the situation we are dealing with.

           if( P3D_iblank ) then

             ! I-blanking is present. Check if the end of the chunk is in
             ! the floating point or integer part.

             sizeRealPart = P3D_nVar*sizeP3D_Real * blockDimCGNS(1,mm) &
                          * blockDimCGNS(2,mm) * blockDimCGNS(3,mm)

             if(sizeRel > sizeRealPart) then

               ! End of chunk is in the integer part. Truncate it to
               ! the nearest integer.

               sizeRel    = (sizeRel - sizeRealPart) &
                          / nBytesPerRecordIntPLOT3D
               chunks(nn) = chunks(nn) + sizeRealPart &
                          + sizeRel*nBytesPerIntPLOT3D
             else

               ! End of chunk is in the floating point part. Truncate
               ! it to the nearest real.

               sizeRel    = sizeRel/sizeP3D_Real
               chunks(nn) = chunks(nn) + sizeRel*sizeP3D_Real

             endif

           else

             ! No i-blanking used. Truncate it to the nearest real.

             sizeRel    = sizeRel/sizeP3D_Real
             chunks(nn) = chunks(nn) + sizeRel*sizeP3D_Real

           endif

         endif

         ! Find out how many block parts are stored in this chunk.
         ! If the end of this chunk corresponds to the end of the record
         ! update mm, because the next chunk does not contain anything
         ! from that global block. Also store the starting index in the
         ! global blockID array for this chunk.

         nPartsPerChunk(nn) = mm - ii + 1
         if(chunks(nn) == cgnsOffset(mm)) mm = mm + 1

         startBlock(nn) = ii

       enddo chunkLoop

       ! Fix the end of the last chunk, determine the number of block
       ! parts and the starting block ID of the last chunk.

       chunks(nChunks)         = cgnsOffset(nGlobalBlocks)
       nPartsPerChunk(nChunks) = nGlobalBlocks - mm + 1
       startBlock(nChunks)     = mm
!
!      ******************************************************************
!      *                                                                *
!      * Determine the number of bytes per chunk that would be          *
!      * communicated, i.e. for which the IO is not local.              *
!      *                                                                *
!      ******************************************************************
!
       ! Loop over the number of chunks.

       chunkCommLoop: do nn=1,nChunks

         nBytesCommChunk(nn) = 0

         ! Loop over the number of global block parts present in
         ! this chunk.

         ii = startBlock(nn) + nPartsPerChunk(nn) - 1
         globalCGNSLoop: do mm=startBlock(nn),ii

           ! Determine the beginning and ending byte distance for this
           ! global block in this chunk.

           lowerBound = max(chunks(nn-1), cgnsOffset(mm-1))
           upperBound = min(chunks(nn),   cgnsOffset(mm))

           ! Determine the indices iBeg, jBeg, kBeg and lBeg to which
           ! the lower bound corresponds.

           offset = lowerBound - cgnsOffset(mm-1) &
                  - nBytesPerRecordIntPLOT3D

           call getBegIndices(offset, sizeP3D_Real, blockDimCGNS(1,mm), &
                              iBeg, jBeg, kBeg, lBeg)

           ! Determine the indices iEnd, jEnd, kEnd and lEnd to which
           ! the upper bound corresponds. Correct the offset if
           ! i-blanking is present. I-blanking arrays are simply ignored
           ! during the reading and are not written.

           offset = cgnsOffset(mm) - upperBound &
                  - nBytesPerRecordIntPLOT3D

           if( P3D_iblank ) &
             offset = offset - nBytesPerIntPLOT3D*blockDimCGNS(1,mm) &
                    * blockDimCGNS(2,mm)*blockDimCGNS(3,mm)

           call getEndIndices(offset, sizeP3D_Real, blockDimCGNS(1,mm), &
                              iEnd, jEnd, kEnd, lEnd)

           ! Store the CGNS block ID corresponding to the currently
           ! active block a bit easier in jj.

           jj = globalBlockIDs(mm)

           ! Loop over the subblocks of global block jj. If the
           ! subblock is stored on a different processor, update
           ! nBytesCommChunk for this chunk.

           subBlockLoop: do ll=1,cgnsDoms(jj)%nSubBlocks

             ! Check if the subblock needs to be communicated.

             commSubBlock: if(cgnsDoms(jj)%procStored(ll) /= myID) then

               ! Determine the range of this subblock in the original
               ! block. This depends on what type of data is
               ! read/written.

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

               ! Adapt the indices iBeg, jBeg, kBeg and lBeg such that
               ! they correspond to the first entity that coincides with
               ! the subblock. This makes the evaluation of the number of
               ! entities a lot easier.
               ! The result are stored in iiBeg, jjBeg, kkBeg and llBeg.

               iiBeg = iBeg; jjBeg = jBeg; kkBeg = kBeg; llBeg = lBeg

               if(iiBeg < iMin) iiBeg = iMin

               if(jjBeg < jMin) then
                 jjBeg = jMin; iiBeg = iMin
               endif

               if(kkBeg < kMin) then
                 kkBeg = kMin; jjBeg = jMin; iiBeg = iMin
               endif

               if(iiBeg > iMax) then
                 iiBeg = iMin; jjBeg = jjBeg + 1
               endif

               if(jjBeg > jMax) then
                 iiBeg = iMin; jjBeg = jMin; kkBeg = kkBeg + 1
               endif

               if(kkBeg > kMax) then
                 iiBeg = iMin; jjBeg = jMin; kkBeg = kMin;
                 llBeg = llBeg + 1
               endif

               ! Adapt the indices iEnd, jEnd, kEnd and lEnd such that
               ! they correspond to the last entity that coincides with
               ! the subblock. This makes the evaluation of the number of
               ! entities a lot easier.
               ! The result are stored in iiEnd, jjEnd, kkEnd and llEnd.

               iiEnd = iEnd; jjEnd = jEnd; kkEnd = kEnd; llEnd = lEnd

               if(iiEnd > iMax) iiEnd = iMax

               if(jjEnd > jMax) then
                 jjEnd = jMax; iiEnd = iMax
               endif

               if(kkEnd > kMax) then
                 kkEnd = kMax; jjEnd = jMax; iiEnd = iMax
               endif

               if(iiEnd < iMin) then
                 iiEnd = iMax; jjEnd = jjEnd - 1
               endif

               if(jjEnd < jMin) then
                 iiEnd = iMax; jjEnd = jMax; kkEnd = kkEnd - 1
               endif

               if(kkEnd < kMin) then
                 iiEnd = iMax; jjEnd = jMax; kkEnd = kMax;
                 llEnd = llEnd - 1
               endif

               ! Determine the number of entities of the subblock that
               ! overlap with the subrange of the global block inside
               ! this chunk.

               sizeComm = iMax - iiBeg + 1 + iiEnd - iMin + 1   &
                        + (iMax - iMin + 1)                     &
                        * (jMax - jjBeg + jjEnd - jMin)         &
                        + (iMax - iMin + 1) * (jMax - jMin + 1) &
                        * (kMax - kkBeg + kkEnd - kMin)         &
                        + (iMax - iMin + 1) * (jMax - jMin + 1) &
                        * (kMax - kMin + 1) * (llEnd - llBeg - 1)

               ! Add the corresponding number of bytes to
               ! nBytesCommChunk(nn) if there is an overlap, i.e. if
               ! sizeComm is positive.

               if(sizeComm > 0) &
                 nBytesCommChunk(nn) = nBytesCommChunk(nn) &
                                     + sizeComm*sizeP3D_Real

             endif commSubBlock
           enddo subBlockLoop
         enddo globalCGNSLoop
       enddo chunkCommLoop

       end subroutine determineBytesCommChunk

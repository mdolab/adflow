!
!      ******************************************************************
!      *                                                                *
!      * File:          determineCommSizesIO.f90                        *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 10-27-2005                                      *
!      * Last modified: 09-13-2007                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine determineCommSizesIO(nGlobalBlocks,  nChunks,        &
                                       sizeP3D_Real,   myChunk,        &
                                       globalBlockIDs, nPartsPerChunk, &
                                       startBlock,     blockDimCGNS,   &
                                       cgnsOffset,     chunks,         &
                                       dataIsRead,     includeHalos)
!
!      ******************************************************************
!      *                                                                *
!      * determineCommSizesIO determines the number of variables to be  *
!      * copied to/from the IO buffer from/to local variables and the   *
!      * number of variables to be communicated to/from other           *
!      * processors. Both these numbers are the same for reading and    *
!      * writing, but the data to be communicated is stored differently.*
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
       integer(kind=intType), intent(in) :: nGlobalBlocks, nChunks
       integer(kind=intType), intent(in) :: sizeP3D_Real,  myChunk

       integer(kind=intType), dimension(nGlobalBlocks), intent(in) :: &
                                                          globalBlockIDs
       integer(kind=intType), dimension(3,nGlobalBlocks), &
                                              intent(in) :: blockDimCGNS

       integer(kind=intType), dimension(nChunks), intent(inout) :: &
                                               nPartsPerChunk, startBlock

       integer(kind=mpi_offset_kind), dimension(0:nGlobalBlocks), &
                                               intent(in) :: cgnsOffset
       integer(kind=mpi_offset_kind), dimension(0:nChunks), &
                                               intent(inout) :: chunks

       logical, intent(in) :: dataIsRead, includeHalos
!
!      Local variables.
!
       integer :: ierr

       integer(kind=intType) :: ii, jj, ll, mm, nn
       integer(kind=intType) :: sizeRel, sizeRealPart
       integer(kind=intType) :: offset, iBeg, jBeg, kBeg, lBeg
       integer(kind=intType) :: iEnd, jEnd, kEnd, lEnd
       integer(kind=intType) :: iiBeg, jjBeg, kkBeg, llBeg
       integer(kind=intType) :: iiEnd, jjEnd, kkEnd, llEnd
       integer(kind=intType) :: iMin, jMin, kMin, iMax, jMax, kMax

       integer(kind=intType), dimension(nProc) :: nDataComm

       integer(kind=mpi_offset_kind) :: lowerBound, upperBound
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! If I do not participate in the IO initialize some variables
       ! and return. The rest of this subroutine is only suited for
       ! myChunk > 0, i.e. processors that participate in the IO.

       if(myChunk == 0) then
         P3D_mySizeIO             = 0
         P3D_myOffset             = 0
         P3D_nIOParts             = 0
         P3D_nProcSend            = 0
         P3D_nProcRecv            = 0
         P3D_nRecordIntegersWrite = 0

         P3D_commPart%offsetIO       = 0
         P3D_commPart%nItemsTotal    = 0
         P3D_commPart%nItemsLocal    = 0
         P3D_commPart%nItemsNonLocal = 0

         nullify(P3D_commPart%blockID,  P3D_commPart%indexW,      &
                 P3D_commPart%posLocal, P3D_commPart%indices,     &
                 P3D_commPart%posComm,  P3D_commPart%posNonLocal)

         ! Dummy allocation for P3D_IOParts to avoid problems during
         ! releasing the memory.

         allocate(P3D_IOParts(P3D_nIOParts), stat=ierr)
         if(ierr /= 0)                            &
           call terminate("determineCommSizesIO", &
                          "Memory allocation failure for P3D_IOParts")

         ! Dummy allocations for the processor ID's and message sizes
         ! to avoid problems later on. This depends whether data is read
         ! or written.

         if( dataIsRead ) then

           allocate(P3D_procSend(P3D_nProcSend), &
                    P3D_sendSize(0:P3D_nProcSend), stat=ierr)
           if(ierr /= 0)                            &
             call terminate("determineCommSizesIO", &
                            "Memory allocation failure for P3D_procSend &
                            &and P3D_sendSize")
           P3D_sendSize(0) = 0

         else

           allocate(P3D_procRecv(P3D_nProcRecv), &
                  P3D_recvSize(0:P3D_nProcRecv), stat=ierr)
           if(ierr /= 0)                            &
             call terminate("determineCommSizesIO", &
                            "Memory allocation failure for P3D_procRecv &
                            &and P3D_recvSize")
           P3D_recvSize(0) = 0

         endif

         ! Make the return.

         return
       endif

       ! Check if my starting position is not somewhere in an iblank
       ! array.

       if( P3D_iblank ) then

         ! Store the starting global block ID a bit easier in nn and
         ! determine the starting position of my chunk relative to
         ! the starting position of the data of that global block
         ! as well as the size of the data we are interested in.

         nn           = startBlock(myChunk)
         sizeRel      = chunks(myChunk-1) - nBytesPerRecordIntPLOT3D &
                      - cgnsOffset(nn-1)
         sizeRealPart = P3D_nVar*sizeP3D_Real*blockDimCGNS(1,nn) &
                      * blockDimCGNS(2,nn)   *blockDimCGNS(3,nn)

         ! Check if the starting position is in the iblank part.

         if(sizeRel >= sizeRealPart) then

           ! Starting position is in the iblanking. Modify the starting
           ! position of my chunk.

           chunks(myChunk-1)       = cgnsOffset(nn)
           nPartsPerChunk(myChunk) = nPartsPerChunk(myChunk) - 1
           startBlock(myChunk)     = startBlock(myChunk) + 1
         endif
       endif
!
!      ******************************************************************
!      *                                                                *
!      * Determine the number of variables to be copied locally and to  *
!      * be communicated.                                               *
!      *                                                                *
!      ******************************************************************
!
       ! Determine the amount of data I have to read and the offset
       ! from the beginning. The max function for P3D_mySizeIO is there
       ! because of possible shift done for iblanking, see above. It is
       ! very unlikely that this clipping is ever necessary.

       P3D_myOffset = chunks(myChunk-1)
       P3D_mySizeIO = chunks(myChunk) - P3D_myOffset
       P3D_mySizeIO = max(P3D_mySizeIO, 0)

       ! Determine the number of global (sub)blocks in this chunk and
       ! allocate the memory for P3D_IOParts.

       P3D_nIOParts = nPartsPerChunk(myChunk)

       allocate(P3D_IOParts(P3D_nIOParts), stat=ierr)
       if(ierr /= 0)                            &
         call terminate("determineCommSizesIO", &
                        "Memory allocation failure for P3D_IOParts")

       ! Nullify the pointers of P3D_IOParts.

       do nn=1,P3D_nIOParts
         nullify(P3D_IOParts(nn)%blockID,  P3D_IOParts(nn)%indexW,      &
                 P3D_IOParts(nn)%posLocal, P3D_IOParts(nn)%indices,     &
                 P3D_IOParts(nn)%posComm,  P3D_IOParts(nn)%posNonLocal)
       enddo

       ! Also nullify the pointers of P3D_commPart.

       nullify(P3D_commPart%blockID,  P3D_commPart%indexW,      &
               P3D_commPart%posLocal, P3D_commPart%indices,     &
               P3D_commPart%posComm,  P3D_commPart%posNonLocal)

       ! Loop over the number of global (sub)blocks in my chunk to
       ! determine the actual numbers. In nDataComm the number of items
       ! which I have to send to/receive from other processors is stored.

       nDataComm = 0

       globalBlockLoop: do nn=1,P3D_nIOParts

         ! Store in mm the entry globalBlockIDs and in jj the
         ! corresponding CGNS block ID.

         mm = startBlock(myChunk) + nn - 1
         jj = globalBlockIDs(mm)

         ! Determine the lower and upper bound for this global block
         ! part inside this chunk. In the lower bound the record integer
         ! is taken into account, because it does contain information
         ! of interest. The formulation used for the upper bound is also
         ! valid if iblank data is present.

         lowerBound = cgnsOffset(mm-1) + nBytesPerRecordIntPLOT3D
         lowerBound = max(lowerBound, chunks(myChunk-1))

         upperBound = cgnsOffset(mm-1) + nBytesPerRecordIntPLOT3D &
                    + P3D_nVar*sizeP3D_Real*blockDimCGNS(1,mm)    &
                    * blockDimCGNS(2,mm)   *blockDimCGNS(3,mm)
         upperBound = min(upperBound, chunks(myChunk))

         ! Determine the offset of this part relative to the starting
         ! position of this chunk and the total number of items for
         ! this global block in this chunk.

         P3D_IOParts(nn)%offsetIO    = lowerBound - chunks(myChunk-1)
         P3D_IOParts(nn)%nItemsTotal = (upperBound - lowerBound) &
                                     / sizeP3D_Real
         P3D_IOParts(nn)%nItemsTotal = max(P3D_IOParts(nn)%nItemsTotal, &
                                           0_intType)

         ! On some platforms it is necessary that the offset equals an
         ! integer multiple of sizeP3D_Real. The variable offsetBuffer
         ! is present to accomplish that.

         P3D_IOParts(nn)%offsetBuffer = P3D_IOParts(nn)%offsetIO      &
                                      - mod(P3D_IOParts(nn)%offsetIO, &
                                            sizeP3D_Real)

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

         ! Initialization of the number of local and nonlocal items for
         ! this part.

         P3D_IOParts(nn)%nItemsLocal    = 0
         P3D_IOParts(nn)%nItemsNonLocal = 0
!
!        ****************************************************************
!        *                                                              *
!        * Loop over all the subblocks of global block jj and determine *
!        * the number of overlapping items for each of them. Update     *
!        * nItemsLocal or nItemsNonLocal accordingly.                   *
!        *                                                              *
!        ****************************************************************
!
         allSubBlockLoop: do ll=1,cgnsDoms(jj)%nSubBlocks

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

           ! Adapt the indices iBeg, jBeg, kBeg and lBeg such that they
           ! correspond to the first entity that coincides with the
           ! subblock. This makes the evaluation of the number of
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
             iiBeg = iMin; jjBeg = jMin; kkBeg = kMin; llBeg = llBeg + 1
           endif

           ! Adapt the indices iEnd, jEnd, kEnd and lEnd such that they
           ! correspond to the last entity that coincides with the
           ! subblock. This makes the evaluation of the number of
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
             iiEnd = iMax; jjEnd = jMax; kkEnd = kMax; llEnd = llEnd - 1
           endif

           ! Determine the number of entities of the subblock that
           ! overlap with the subrange of the global block inside
           ! this chunk.

           sizeRel = iMax - iiBeg + 1 + iiEnd - iMin + 1   &
                   + (iMax - iMin + 1)                     &
                   * (jMax - jjBeg + jjEnd - jMin)         &
                   + (iMax - iMin + 1) * (jMax - jMin + 1) &
                   * (kMax - kkBeg + kkEnd - kMin)         &
                   + (iMax - iMin + 1) * (jMax - jMin + 1) &
                   * (kMax - kMin + 1) * (llEnd - llBeg - 1)

           ! Check if the size is positive, because then there is
           ! an overlap and data must be updated.

           if(sizeRel > 0) then

             ! Update either nItemsLocal or nItemsNonLocal, depending
             ! whether this subblock is stored locally or not.

             if(cgnsDoms(jj)%procStored(ll) == myID) then
               P3D_IOParts(nn)%nItemsLocal = &
                  P3D_IOParts(nn)%nItemsLocal + sizeRel
             else

               ! Subblock is not stored on this processor.
               ! Update nItemsNonLocal and nDataComm of the
               ! corresponding processor.

               P3D_IOParts(nn)%nItemsNonLocal = &
                  P3D_IOParts(nn)%nItemsNonLocal + sizeRel

               ii = cgnsDoms(jj)%procStored(ll) + 1
               nDataComm(ii) = nDataComm(ii) + sizeRel

             endif
           endif

         enddo allSubBlockLoop

       enddo globalBlockLoop
!
!      ******************************************************************
!      *                                                                *
!      * Determine the number of processors to/from which I have to     *
!      * send/receive data. Also determine the corresponding processor  *
!      * ID's and the amount of data to be communicated.                *
!      *                                                                *
!      ******************************************************************
!
       testReadWrite: if( dataIsRead ) then

         ! Data is read and must be sent to other processors.
         ! Determine the number of processors to which I have
         ! to sent data.

         P3D_nProcSend = 0
         do nn=1,nProc
           if(nDataComm(nn) > 0) P3D_nProcSend = P3D_nProcSend + 1
         enddo

         ! Allocate the memory for P3D_procSend and P3D_sendSize.

         allocate(P3D_procSend(P3D_nProcSend), &
                  P3D_sendSize(0:P3D_nProcSend), stat=ierr)
         if(ierr /= 0)                            &
           call terminate("determineCommSizesIO", &
                          "Memory allocation failure for P3D_procSend &
                          &and P3D_sendSize")

         ! Determine the values of these arrays.

         P3D_sendSize(0) = 0
         mm = 0
         do nn=1,nProc
           if(nDataComm(nn) > 0) then
             mm = mm + 1
             P3D_procSend(mm) = nn - 1
             P3D_sendSize(mm) = P3D_sendSize(mm-1) + nDataComm(nn)
           endif
         enddo

       else testReadWrite

         ! Data is written and must be received from other processors.
         ! Determine the number of processors from which I have
         ! to receive data.

         P3D_nProcRecv = 0
         do nn=1,nProc
           if(nDataComm(nn) > 0) P3D_nProcRecv = P3D_nProcRecv + 1
         enddo

         ! Allocate the memory for P3D_procRecv and P3D_recvSize.

         allocate(P3D_procRecv(P3D_nProcRecv), &
                  P3D_recvSize(0:P3D_nProcRecv), stat=ierr)
         if(ierr /= 0)                            &
           call terminate("determineCommSizesIO", &
                          "Memory allocation failure for P3D_procRecv &
                          &and P3D_recvSize")

         ! Determine the values of these arrays.

         P3D_recvSize(0) = 0
         mm = 0
         do nn=1,nProc
           if(nDataComm(nn) > 0) then
             mm = mm + 1
             P3D_procRecv(mm) = nn - 1
             P3D_recvSize(mm) = P3D_recvSize(mm-1) + nDataComm(nn)
           endif
         enddo

       endif testReadWrite

       end subroutine determineCommSizesIO

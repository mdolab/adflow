!
!      ******************************************************************
!      *                                                                *
!      * File:          sendBufNonLocalIO.f90                           *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 10-28-2005                                      *
!      * Last modified: 09-13-2007                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine sendBufNonLocalIO(nGlobalBlocks,  nChunks,       &
                                    sizeP3D_Real,   myChunk,       &
                                    globalBlockIDs, startBlock,    &
                                    blockDimCGNS,   cgnsOffset,    &
                                    chunks,         P3D_nProcComm, &
                                    P3D_procComm,   P3D_commSize,  &
                                    sendBuf,        dataIsRead,    &
                                    includeHalos)
!
!      ******************************************************************
!      *                                                                *
!      * sendVarNonLocalIO determines the send buffer of the indices    *
!      * for the non local data needed for the IO. The messages         *
!      * themselves are not send here. It could have been included, but *
!      * then the corresponding calls to mpi_waitany are in a different *
!      * subroutine, which leads to unreadable code.                    *
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
       integer(kind=intType), intent(in) :: P3D_nProcComm

       integer(kind=intType), dimension(nGlobalBlocks), intent(in) :: &
                                                          globalBlockIDs
       integer(kind=intType), dimension(3,nGlobalBlocks), &
                                              intent(in) :: blockDimCGNS

       integer(kind=intType), dimension(nChunks), intent(in) :: &
                                                              startBlock

       integer(kind=intType), dimension(P3D_nProcComm), intent(in) :: &
                                                            P3D_procComm
       integer(kind=intType), dimension(0:P3D_nProcComm), intent(in) :: &
                                                            P3D_commSize

       integer(kind=intType), dimension(5,P3D_commSize(P3D_nProcComm)), &
                                                   intent(out) :: sendBuf

       integer(kind=mpi_offset_kind), dimension(0:nGlobalBlocks), &
                                               intent(in) :: cgnsOffset
       integer(kind=mpi_offset_kind), dimension(0:nChunks), &
                                               intent(in) :: chunks

       logical, intent(in) :: dataIsRead, includeHalos
!
!      Local variables.
!
       integer :: ierr

       integer(kind=intType) :: ii, jj, kk, ll, mm, nn, i, j
       integer(kind=intType) :: offset, iBeg, jBeg, kBeg, lBeg
       integer(kind=intType) :: iGMin, jGMin, kGMin, iGMax, jGMax, kGMax
       integer(kind=intType) :: iMin, jMin, kMin, iMax, jMax, kMax
       integer(kind=intType) :: iOffset, jOffset, kOffset
       integer(kind=intType) :: ic, jc, kc, lc

       integer(kind=intType), dimension(0:P3D_nProcComm) :: counter
       integer(kind=intType), dimension(0:nProc-1)       :: procMapping

       integer(kind=mpi_offset_kind) :: lowerBound
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Initialize counter to P3D_commSize. It is used to determine the
       ! correct position in the send buffer.

       do nn=0,P3D_nProcComm
         counter(nn) = P3D_commSize(nn)
       enddo

       ! Determine the mapping from the processor ID to the position
       ! in P3D_procSend.

       do nn=1,P3D_nProcComm
         procMapping(P3D_procComm(nn)) = nn
       enddo

       ! Loop to determine the indices for copying the data from the
       ! read buffer into the send buffer.

       sendDataLoop: do nn=1,P3D_nIOParts

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

         ! Allocate the memory for posComm and posNonLocal.

         ii = P3D_IOParts(nn)%nItemsNonLocal
         allocate(P3D_IOParts(nn)%posComm(ii), &
                  P3D_IOParts(nn)%posNonLocal(ii), stat=ierr)
         if(ierr /= 0)                         &
           call terminate("sendBufNonLocalIO", &
                          "Memory allocation failure for posComm and &
                          &posNonLocal")

         ! Loop over the subblocks of this global block.
         ! kk is the counter for the nonlocal entries.

         kk = 0

         sendSubBlockLoop: do ll=1,cgnsDoms(jj)%nSubBlocks

           ! Check if this subblock needs to be communicated.

           testComm: if(cgnsDoms(jj)%procStored(ll) /= myID) then

             ! Determine the entry for this processor in P3D_procComm.

             ii = cgnsDoms(mm)%procStored(ll)
             ii = procMapping(ii) - 1

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

             ! Determine the offset between global and local indices.

             iOffset = cgnsDoms(jj)%iBegOr(ll) - 1
             jOffset = cgnsDoms(jj)%jBegOr(ll) - 1
             kOffset = cgnsDoms(jj)%kBegOr(ll) - 1

             ! Loop over the total number of entities in this IO part
             ! and store the indices which coincide.

             ic = iBeg; jc = jBeg; kc = kBeg; lc = lBeg

             itemsNonLocalLoop: do i=1,P3D_IOParts(nn)%nItemsTotal

               ! Check if the current i, j, k indices coincide with the
               ! subblock considered.

               if(ic >= iMin .and. jc >= jMin .and. kc >= kMin .and. &
                  ic <= iMax .and. jc <= jMax .and. kc <= kMax) then

                 ! Update the counter for this processor and store it
                 ! a bit easier in j. Update the counter jj as well.

                 counter(ii) = counter(ii) + 1
                 j  = counter(ii)
                 kk = kk + 1

                 ! Store the position of both the IO and communication
                 ! buffer.

                 P3D_IOParts(nn)%posComm(kk)     = j
                 P3D_IOParts(nn)%posNonLocal(kk) = i

                 ! Store the index information and the local block ID
                 ! on the receiving side in the send buffer.

                 sendBuf(1,j) = ic - iOffset
                 sendBuf(2,j) = jc - jOffset
                 sendBuf(3,j) = kc - kOffset
                 sendBuf(4,j) = lc - 1
                 sendBuf(5,j) = cgnsDoms(jj)%localBlockID(ll)
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
             enddo itemsNonLocalLoop

           endif testComm
         enddo sendSubBlockLoop

         ! Test in debug mode that the desired number of data is
         ! stored in P3D_IOParts(nn).

         if( debug ) then
           if(kk /= P3D_IOParts(nn)%nItemsNonLocal) &
             call terminate("sendBufNonLocalIO",    &
                            "kk differs from nItemsNonLocal")
         endif

       enddo sendDataLoop

       end subroutine sendBufNonLocalIO

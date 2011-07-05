!
!      ******************************************************************
!      *                                                                *
!      * File:          determineMyChunk.f90                            *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 10-27-2005                                      *
!      * Last modified: 10-28-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine determineMyChunk(nBytesCommChunk,  nChunks, &
                                   IParticipateInIO, myChunk)
!
!      ******************************************************************
!      *                                                                *
!      * determineMyChunk determines the chunk the processor is         *
!      * responsible for during the IO. The value of myChunk is         *
!      * irrelevant if the processor does not participate in the IO.    *
!      * This will only happen for the writing if not all blocks are    *
!      * written. In that case myChunk is set to 0.                     *
!      *                                                                *
!      ******************************************************************
!
       use communication
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType), intent(in) :: nChunks
       integer(kind=intType), dimension(nChunks), intent(inout) :: &
                                                          nBytesCommChunk
       logical, intent(in) :: IParticipateInIO

       integer(kind=intType), intent(out) :: myChunk
!
!      Local variables.
!
       integer :: ierr

       integer(kind=intType) :: ii, jj, ll, mm, nn
       integer(kind=intType) :: nDiffChunks, nChunksAssigned

       integer(kind=intType), dimension(2) :: sendBuf

       integer(kind=intType), dimension(0:nChunks) :: nChunkPerBytes
       integer(kind=intType), dimension(nChunks)   :: chunkPerBytes
       integer(kind=intType), dimension(nChunks)   :: tmp
       integer(kind=intType), dimension(nChunks)   :: procChunk
       integer(kind=intType), dimension(2,nProc)   :: recvBuf
!
!      Function definition.
!
       integer(kind=intType) :: bsearchIntegers
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
!      ******************************************************************
!      *                                                                *
!      * Determine a sorted version of nBytesCommChunk and the arrays   *
!      * nChunkPerBytes and chunkPerBytes, which define the mapping to  *
!      * the chunks.                                                    *
!      *                                                                *
!      ******************************************************************
!
       ! Copy nBytesCommChunk to tmp for the bindary search later on and
       ! sort nBytesCommChunk in increasing order.

       tmp = nBytesCommChunk
       call qsortIntegers(nBytesCommChunk, nChunks)

       ! Determine the number of different values as well as the
       ! multiplicity of nBytesCommChunk. This is accomplished by the
       ! array nChunkPerBytes how many chunks are present with that
       ! particular number of bytes to communicate. nChunkPerBytes is
       ! in cumulative storage format.

       nDiffChunks       = 1
       nChunkPerBytes(0) = 0
       nChunkPerBytes(1) = 1

       do nn=2,nChunks
         if(nBytesCommChunk(nn) == nBytesCommChunk(nDiffChunks)) then
           nChunkPerBytes(nDiffChunks) = nChunkPerBytes(nDiffChunks) + 1
         else
           nDiffChunks                  = nDiffChunks + 1
           nBytesCommChunk(nDiffChunks) = nBytesCommChunk(nn)
           nChunkPerBytes(nDiffChunks)  = nChunkPerBytes(nDiffChunks-1) + 1
         endif
       enddo

       ! Determine the chunk ID's per value of nBytesCommChunk.

       do nn=1,nChunks
         ii = bsearchIntegers(tmp(nn), nBytesCommChunk, nDiffChunks) - 1
         nChunkPerBytes(ii)                = nChunkPerBytes(ii) + 1
         chunkPerBytes(nChunkPerBytes(ii)) = nn
       enddo

       ! Restore nChunkPerBytes in its original form again.

       do nn=(nDiffChunks-1),1,-1
         nChunkPerBytes(nn) = nChunkPerBytes(nn-1)
       enddo
       nChunkPerBytes(0) = 0
!
!      ******************************************************************
!      *                                                                *
!      * Determine the chunk this processor will be responsible for.    *
!      *                                                                *
!      ******************************************************************
!
       ! Initialization of nChunksAssigned, the counter nn and
       ! the control variable tmp.

       nChunksAssigned = 0
       nn              = 1
       tmp             = 0

       ! Iterative method to determine the chunk distribution.

       chunkAssignOuter: do ll=1,nDiffChunks
         chunkAssignInner: do mm=(nChunkPerBytes(ll-1)+1), &
                                  nChunkPerBytes(ll)

           ! Criterion to exit the iterative method.

           if(nChunksAssigned == nChunks) exit chunkAssignOuter

           ! Set the chunk I want to read if I have not been assigned
           ! a chunk yet and I participate in the IO.

           if(IParticipateInIO .and. nn == mm) then
             sendBuf(1) = chunkPerBytes(nn)
             sendBuf(2) = nBytesCommChunk(ll)
           else
             sendBuf(1) = 0
             sendBuf(2) = 0
           endif

           ! Gather the values of sendBuf on all processors.

           call mpi_allgather(sendBuf, 2, sumb_integer, recvBuf, 2, &
                              sumb_integer, SUmb_comm_world, ierr)

           ! Loop over the chunks and assign the ones that have
           ! not been assigned yet. In case of multiple candidates for
           ! a chunk take the one with the MAXIMUM number of bytes to
           ! be communicated. This seems strange, but if that one is
           ! not taken in this round it will need to communicate even
           ! more in the next round.

           do ii=1,nProc
             if(recvBuf(1,ii) > 0) then
               if(tmp(recvBuf(1,ii)) == 0) then

                 ! Chunk has not been assigned yet. Do so now.

                 procChunk(recvBuf(1,ii)) = ii - 1
                 tmp(recvBuf(1,ii))       = mm
                 nChunksAssigned          = nChunksAssigned + 1

               else if(tmp(recvBuf(1,ii)) == mm) then

                 ! Chunk has been assigned in this round. Check if the
                 ! current processor is a better candidate. Note that
                 ! nChunksAssigned should NOT be updated.

                 jj = procChunk(recvBuf(1,ii)) + 1
                 if(recvBuf(2,ii) > recvBuf(2,jj)) &
                   procChunk(recvBuf(1,ii)) = ii - 1

               endif
             endif
           enddo

           ! Check if I had been assigned a chunk. If not update nn
           ! for the next round.

           if(IParticipateInIO .and. nn == mm) then
             if(procChunk(chunkPerBytes(nn)) /= myID) nn = nn + 1
           endif

         enddo chunkAssignInner
       enddo chunkAssignOuter

       ! I am responsible for chunkPerBytes(nn) if I participate in the
       ! IO. Store it in a more readable variable. If I do not
       ! participate simply set myChunk to 0.

       myChunk = 0
       if( IParticipateInIO ) myChunk = chunkPerBytes(nn)

       end subroutine determineMyChunk

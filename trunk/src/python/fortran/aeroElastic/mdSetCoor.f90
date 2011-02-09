!
!      ******************************************************************
!      *                                                                *
!      * File:          mdSetCoor.F90                                   *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 03-04-2004                                      *
!      * Last modified: 10-25-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine mdSetCoor(sps, nSubBlocks, cgnsIDs, ranges, &
                            coor, nNode)
!
!      ******************************************************************
!      *                                                                *
!      * MdSetCoor set the coordinates of the current ground level      *
!      * for the given spectral solution to the coordinates given in    *
!      * the header. These are given as subblock(s) of the original     *
!      * cgns blocks. This routine distributes them to the correct      *
!      * processor, after which the coordinates of the blocks are       *
!      * overwritten.                                                   *
!      *                                                                *
!      ******************************************************************
!
       use block
       use cgnsGrid
       use communication
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType), intent(in) :: sps, nSubBlocks, nNode

       integer(kind=intType), dimension(nSubBlocks),     &
                                                   intent(in) :: cgnsIDs
       integer(kind=intType), dimension(3,2,nSubBlocks), &
                                                   intent(in) :: ranges

       real(kind=realType), dimension(3,nNode), intent(in) :: coor
!
!      Local variables.
!
       integer :: ierr, procID, sizeReal

       integer, dimension(mpi_status_size) :: status
       integer, dimension(nProc)           :: sizeMessage

       integer, dimension(:), allocatable :: intRequests, realRequests

       integer(kind=intType) :: nn, mm, ll, kk, ii, jj, j
       integer(kind=intType) :: nMessSend, nMessRecv, ncgnsIDs

       integer(kind=intType), dimension(7) :: intBuf
       integer(kind=intType), dimension(0:nProc-1) :: nMess2Proc

       integer(kind=intType), dimension(0:nDom) :: nBlockPerCGNS
       integer(kind=intType), dimension(nDom)   :: blockPerCGNS
       integer(kind=intType), dimension(nDom)   :: sortedCGNSIDs, tmp

       real(kind=realType), dimension(:,:), allocatable :: realBuf
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
!      * Determine the number and the ID's of the local blocks per cgns *
!      * block present on this processor.                               *
!      *                                                                *
!      ******************************************************************
!
       ! Copy the cgns block ID's of the local blocks in sortedCGNSIDs
       ! and sort them in increasing order.

       do nn=1,nDom
         sortedCGNSIDs(nn) = flowDoms(nn,1,sps)%cgnsBlockID
       enddo

       call qsortIntegers(sortedCGNSIDs, nDom)

       ! Determine the number of blocks in flowDoms per cgns block and
       ! remove the multiple entries in sortedCGNSIDs. Also the number
       ! of different cgns id's is determined.

       nBlockPerCGNS(0) = 0
       nBlockPerCGNS(1) = 1
       ncgnsIDs = 1

       do nn=2,nDom
         if(sortedCGNSIDs(nn) == sortedCGNSIDs(ncgnsIDs)) then
           nBlockPerCGNS(ncgnsIDs) = nBlockPerCGNS(ncgnsIDs) + 1
         else
           ncgnsIDs = ncgnsIDs + 1
           sortedCGNSIDs(ncgnsIDs) = sortedCGNSIDs(nn)
           nBlockPerCGNS(ncgnsIDs) = nBlockPerCGNS(ncgnsIDs-1) + 1
         endif
       enddo

       ! Determine the block ID's per cgns block, which are stored in
       ! blockPerCGNS. First copy the starting places in tmp, which
       ! serves as a counter in the loop below.

       do nn=1,ncgnsIDs
         tmp(nn) = nBlockPerCGNS(nn-1)
       enddo

       ! Loop over the number of local blocks to search for the cgns
       ! block id in sortedCGNSIDs, update the corresponding counter
       ! and store the block ID in blockPerCGNS.

       do nn=1,nDom
         mm = bsearchIntegers(flowDoms(nn,1,sps)%cgnsBlockID, &
                              sortedCGNSIDs, ncgnsIDs)
         tmp(mm) = tmp(mm) + 1
         blockPerCGNS(tmp(mm)) = nn
       enddo
!
!      ******************************************************************
!      *                                                                *
!      * Determine the number of messages to be sent and received.      *
!      *                                                                *
!      ******************************************************************
!
       ! Find out to how many messages i will send to the processors.
       ! Note that it is possible to send multiple messages to one
       ! processor. This is easier, because then it is possible to use
       ! coor directly as a send buffer. It these messages were to be
       ! concatenated into one message, copying the coordinates would
       ! be necessary. As the messages are relatively big anyway, this
       ! is not really a problem.

       nMess2Proc = 0
       nMessSend   = 0

       do nn=1,nSubBlocks
         ii = cgnsIDs(nn)
         do j=1,cgnsDoms(ii)%nSubBlocks
           procID = cgnsDoms(ii)%procStored(j)
           nMess2Proc(procID) = nMess2Proc(procID) + 1
         enddo
       enddo

       ! Make sure that no messages are sent to myself.

       nMess2Proc(myID) = 0

       ! Determine the total messages I have to send.

       nMessSend = 0
       do nn=0,(nProc-1)
         nMessSend = nMessSend + nMess2Proc(nn)
       enddo

       ! Determine the number of messages I will receive.

       sizeMessage = 1
       call mpi_reduce_scatter(nMess2Proc,   nMessRecv, sizeMessage,     &
                               sumb_integer, mpi_sum,   SUmb_comm_world, &
                               ierr)

       ! Allocate the memory for the requests needed for the nonBlocking
       ! sends. As it is theoretically possible that more messages than
       ! processors are sent, sendRequests and revRequests cannot
       ! be used.

       allocate(intRequests(nMessSend), realRequests(nMessSend), &
                stat=ierr)
       if(ierr /= 0)                 &
         call terminate("mdSetCoor", &
                        "Memory allocation failure for requests")
!
!      ******************************************************************
!      *                                                                *
!      * Communicate the coordinates and update them in flowDoms.        *
!      *                                                                *
!      ******************************************************************
!
       ! Loop over the number of subblocks to either send the data or
       ! update my own coordinates.

       mm = 1
       jj = 0

       do nn=1,nSubBlocks
         ii = cgnsIDs(nn)

         ! Determine the number of points for this subblock.

         ll = (ranges(1,2,nn) - ranges(1,1,nn) + 1) &
            * (ranges(2,2,nn) - ranges(2,1,nn) + 1) &
            * (ranges(3,2,nn) - ranges(3,1,nn) + 1)

         ! Loop over the number processors on which this cgns block
         ! is stored.

         do j=1,cgnsDoms(ii)%nSubBlocks
           procID = cgnsDoms(ii)%procStored(j)

           ! Test whether the data is stored locally or must be sent
           ! to an other processor.

           if(procID == myID) then

             ! Local data. Overwrite the coordinates of the
             ! corresponding flowDoms.

             call mdSetMyCoor(sps, ii, ranges(1,1,nn), coor(1,mm), &
                              ncgnsIDs, sortedCGNSIDs,             &
                              nBlockPerCGNS, blockPerCGNS)
           else

             ! Data must be sent. First copy the integer data into
             ! intBuf.

             intBuf(1) = ii
             intBuf(2) = ranges(1,1,nn)
             intBuf(3) = ranges(2,1,nn)
             intBuf(4) = ranges(3,1,nn)
             intBuf(5) = ranges(1,2,nn)
             intBuf(6) = ranges(2,2,nn)
             intBuf(7) = ranges(3,2,nn)

             ! Send the integer and real info. Use nonBlocking sends
             ! to avoid deadlock.

             jj = jj + 1
             call mpi_isend(intBuf, 7, sumb_integer, procID, procID, &
                            SUmb_comm_world, intRequests(jj), ierr)

             kk = 3*ll
             sizeReal = kk
             call mpi_isend(coor(1,mm), sizeReal, sumb_real,     &
                            procID, procID+1, SUmb_comm_world, &
                            realRequests(jj), ierr)
           endif

         enddo

         ! Update mm for the next sub block.

         mm = mm + ll

       enddo

       ! Loop over the number of messages I will receive.

       recvLoop: do nn=1,nMessRecv

         ! Block until a message arrives and determine the
         ! corresponding processor.

         call mpi_probe(mpi_any_source, myID, SUmb_comm_world, &
                       status, ierr)
         procID = status(mpi_source)

         ! Receive the integer message. Use a blocking receive, because
         ! the message is already here.

         call mpi_recv(intBuf, 7, sumb_integer, procID, myID, &
                       SUmb_comm_world, status, ierr)

         ! Determine the size of the coordinate message and allocate
         ! the memory for the receive buffer. To see how intBuf is
         ! constructed, see the send loop. In intBuf(2) till intBuf(7)
         ! the range of the subblock is stored.

         ll = (intBuf(5) - intBuf(2) + 1) &
            * (intBuf(6) - intBuf(3) + 1) &
            * (intBuf(7) - intBuf(4) + 1)

         allocate(realBuf(3,ll), stat=ierr)
         if(ierr /= 0)                   &
           call terminate("mdSetCoor", &
                          "Memory allocation failure for realBuf")

         ! Receive the coordinate message and set the corresponding
         ! coordinates.

         mm = 3*ll
         sizeReal = mm
         call mpi_recv(realBuf, sizeReal, sumb_real, procID, &
                       myID+1, SUmb_comm_world, status, ierr)

         call mdSetMyCoor(sps, intBuf(1), intBuf(2), realBuf, &
                          ncgnsIDs, sortedCGNSIDs,            &
                          nBlockPerCGNS, blockPerCGNS)

         ! Release the memory of realBuf again.

         deallocate(realBuf, stat=ierr)
         if(ierr /= 0)                 &
           call terminate("mdSetCoor", &
                          "Deallocation failure for realBuf")
       enddo recvLoop

       ! Complete the nonBlocking sends. SizeReal is used as working
       ! variable, because it is an integer.

       sizeReal = nMessSend
       do nn=1,nMessSend
         call mpi_waitany(sizeReal, intRequests,  procID, status, ierr)
         call mpi_waitany(sizeReal, realRequests, procID, status, ierr)
       enddo

       ! Release the memory of request.

       deallocate(intRequests, realRequests, stat=ierr)
       if(ierr /= 0)                 &
         call terminate("mdSetCoor", &
                        "Deallocation failure for requests")

       ! Synchronize the processors, because wild cards have been used
       ! in the receives.

       call mpi_barrier(SUmb_comm_world, ierr)

       end subroutine mdSetCoor

       !=================================================================

       subroutine mdSetMyCoor(sps, cgnsID, cgnsRange, coor, &
                              ncgnsIDs, sortedCGNSIDs,      &
                              nBlockPerCGNS, blockPerCGNS)
!
!      ******************************************************************
!      *                                                                *
!      * mdSetMyCoor sets the coordinates of the block(s) for the       *
!      * given spectral solution to coor for the given range of the     *
!      * cgns block. GroundLevel is assumed to be the active grid       *
!      * level in the simulation.                                       *
!      *                                                                *
!      ******************************************************************
!
       use block
       use iteration
       implicit none
!
!      Subroutine arguments
!
       integer(kind=intType), intent(in) :: sps, cgnsID, ncgnsIDs

       integer(kind=intType), dimension(3,2), intent(in) :: cgnsRange
       integer(kind=intType), dimension(0:*), intent(in) :: nBlockPerCGNS
       integer(kind=intType), dimension(*),   intent(in) :: blockPerCGNS
       integer(kind=intType), dimension(*),   intent(in) :: sortedCGNSIDs

       real(kind=realType), dimension(3,*), intent(in) :: coor
!
!      Local variables.
!
       integer(kind=intType) :: i, j, k, ii, jj, kk, nn, mm, ll, count

       logical :: kkOkay, jjOkay, iiOkay
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
       ! Determine the entry of cgnsID in the sorted array of the
       ! locally stored cgns block id's.

       mm = bsearchIntegers(cgnsID, sortedCGNSIDs, ncgnsIDs)

       ! Loop over the number of blocks, which are subblocks of the
       ! given cgns block.

       blockLoop: do nn=(nBlockPerCGNS(mm-1)+1), nBlockPerCGNS(mm)

         ! Store the block ID a bit easier and initialize the counter
         ! count to 1

         ll    = blockPerCGNS(nn)
         count = 1

         ! Loop over the given k-subrange of the cgns block and determine
         ! the possible index in the local block ll.

         kLoop: do k=cgnsRange(3,1),cgnsRange(3,2)
           kk = k - flowDoms(ll,groundLevel,sps)%kBegor + 1

           ! Test whether kk is in the range for the block ll.
           ! Set kkOkay accordingly.

           kkOkay = .false.
           if(kk <= flowDoms(ll,groundLevel,sps)%kl .and. &
              kk >= 1) kkOkay = .true.

           ! Loop over the given j-subrange, etc.

           jLoop: do j=cgnsRange(2,1),cgnsRange(2,2)
             jj = j - flowDoms(ll,groundLevel,sps)%jBegor + 1

             ! Test whether jj is in the range for the block ll.
             ! Set jjOkay accordingly.

             jjOkay = .false.
             if(jj <= flowDoms(ll,groundLevel,sps)%jl .and. &
                jj >= 1) jjOkay = .true.

             ! Loop over the given i-subrange, etc.

             iLoop: do i=cgnsRange(1,1),cgnsRange(1,2)
               ii = i - flowDoms(ll,groundLevel,sps)%iBegor + 1

               ! Test whether ii is in the range for the block ll.
               ! Set iiOkay accordingly.

               iiOkay = .false.
               if(ii <= flowDoms(ll,groundLevel,sps)%il .and. &
                  ii >= 1) iiOkay = .true.

               ! If all three coordinate indices are okay, overwrite
               ! the coordinates of the block.

               if(kkOkay .and. jjOkay .and. iiOkay) then
                 flowDoms(ll,groundLevel,sps)%x(ii,jj,kk,1) = coor(1,count)
                 flowDoms(ll,groundLevel,sps)%x(ii,jj,kk,2) = coor(2,count)
                 flowDoms(ll,groundLevel,sps)%x(ii,jj,kk,3) = coor(3,count)
               endif

               ! Update the counter count.

               count = count + 1

             enddo iLoop
           enddo jLoop
         enddo kLoop

       enddo blockLoop

       end subroutine mdSetMyCoor

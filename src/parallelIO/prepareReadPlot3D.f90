!
!      ******************************************************************
!      *                                                                *
!      * File:          prepareReadPlot3D.f90                           *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 10-20-2005                                      *
!      * Last modified: 09-13-2007                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine prepareReadPlot3D(readHalos)
!
!      ******************************************************************
!      *                                                                *
!      * prepareReadPlot3D sets the Plot3D variables in IOModule such   *
!      * that the reading of the actual data from a Plot3D file can be  *
!      * performed.                                                     *
!      *                                                                *
!      ******************************************************************
!
       use cgnsGrid
       use communication
       use inputIO
       use IOModule
       implicit none
!
!      Subroutine arguments.
!
       logical, intent(in) :: readHalos
!
!      Local variables.
!
       integer :: procID, sizeMessage, ierr

       integer, dimension(mpi_status_size) :: status

       integer(kind=intType) :: ii, nn
       integer(kind=intType) :: nGlobalBlocks, nChunks, myChunk
       integer(kind=intType) :: sizeP3D_Real

       integer(kind=intType), dimension(cgnsNDom)  :: globalBlockIDs
       integer(kind=intType), dimension(nProc)     :: nPartsPerChunk
       integer(kind=intType), dimension(nProc)     :: startBlock
       integer(kind=intType), dimension(nProc)     :: nBytesCommChunk
       integer(kind=intType), dimension(0:nProc-1) :: sendData
       integer(kind=intType), dimension(0:nProc-1) :: recvData

       integer(kind=intType), dimension(3,cgnsNDom) :: blockDimCGNS

       integer(kind=intType), dimension(:,:), allocatable :: sendBuf
       integer(kind=intType), dimension(:,:), allocatable :: recvBuf

       integer(kind=mpi_offset_kind), dimension(0:cgnsNDom) :: cgnsOffset
       integer(kind=mpi_offset_kind), dimension(0:nProc)    :: chunks

       logical :: IParticipateInIO
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Set the values of nGlobalBlocks and globalBlockIDs such that
       ! all blocks are read. This is done such that the routines below
       ! can be used for both reading and writing. For writing it is
       ! possible to write only a part of the blocks.

       nGlobalBlocks = cgnsNDom
       do nn=1,nGlobalBlocks
         globalBlockIDs(nn) = nn
       enddo

       ! Determine the size of the floating point type, depending
       ! on that floating point type.

       select case (P3D_Precision)
         case (precisionSingle)
           sizeP3D_Real = 4

         case (precisionDouble)
           sizeP3D_Real = 8
       end select

       ! Determine for every chunk to be read the number of bytes I
       ! would have to communicate to other processors.

       call determineBytesCommChunk(nGlobalBlocks,   sizeP3D_Real,    &
                                    nChunks,         globalBlockIDs,  &
                                    nPartsPerChunk,  startBlock,      &
                                    nBytesCommChunk, blockDimCGNS,    &
                                    cgnsOffset,      chunks,          &
                                    .true.,          readHalos,       &
                                    IParticipateInIO)

       ! Determine the chunk for which I am responsible.

       call determineMyChunk(nBytesCommChunk,  nChunks, &
                             IParticipateInIO, myChunk)

       ! Determine the number of variables for the local IO as well as
       ! the number to be send to other processors.

       call determineCommSizesIO(nGlobalBlocks,  nChunks,        &
                                 sizeP3D_Real,   myChunk,        &
                                 globalBlockIDs, nPartsPerChunk, &
                                 startBlock,     blockDimCGNS,   &
                                 cgnsOffset,     chunks,         &
                                 .true.,         readHalos)
!
!      ******************************************************************
!      *                                                                *
!      * Determine the number of processors, the corresponding          *
!      * processor ID's and the size of the messages which I have to    *
!      * receive.                                                       *
!      *                                                                *
!      ******************************************************************
!
       ! Set the entries of the processors to which I send data.

       sendData = 0
       do nn=1,P3D_nProcSend
         sendData(P3D_procSend(nn)) = P3D_sendSize(nn) &
                                    - P3D_sendSize(nn-1)
       enddo

       ! Perform an all to all communication, such that the processors
       ! know how many messages will be received as well as their size.

       call mpi_alltoall(sendData, 1, sumb_integer, &
                         recvData, 1, sumb_integer, SUmb_comm_world, ierr)

       ! Determine the number of processors from which I receive data.

       P3D_nProcRecv = 0
       do nn=0,(nProc-1)
         if(recvData(nn) > 0) P3D_nProcRecv = P3D_nProcRecv + 1
       enddo

       ! Allocate the memory for P3D_procRecv and P3D_recvSize.

       allocate(P3D_procRecv(P3D_nProcRecv), &
                P3D_recvSize(0:P3D_nProcRecv), stat=ierr)
       if(ierr /= 0)                         &
         call terminate("prepareReadPlot3D", &
                        "Memory allocation failure for P3D_procRecv &
                        &and P3D_recvSize")

       ! Determine the values of these arrays. Note that P3D_recvSize
       ! is in cumulative storage format.

       ii = 0
       P3D_recvSize(0) = 0

       do nn=0,(nProc-1)
         if(recvData(nn) > 0) then
           ii = ii + 1
           P3D_procRecv(ii) = nn
           P3D_recvSize(ii) = P3D_recvSize(ii-1) + recvData(nn)
         endif
       enddo
!
!      ******************************************************************
!      *                                                                *
!      * The sending of the index data to the processors that will      *
!      * receive data from me during the actual IO.                     *
!      *                                                                *
!      ******************************************************************
!
       ! Allocate the memory for the send buffer.

       ii = P3D_sendSize(P3D_nProcSend)
       allocate(sendBuf(5,ii), stat=ierr)
       if(ierr /= 0)                         &
         call terminate("prepareReadPlot3D", &
                        "Memory allocation failure for sendBuf")

       ! The routine which fill the actual send buffer.

       call sendBufNonLocalIO(nGlobalBlocks,  nChunks,       &
                              sizeP3D_Real,   myChunk,       &
                              globalBlockIDs, startBlock,    &
                              blockDimCGNS,   cgnsOffset,    &
                              chunks,         P3D_nProcSend, &
                              P3D_procSend,   P3D_sendSize,  &
                              sendBuf,        .true.,        &
                              readHalos)

       ! Send the messages to the receiving processors. Use nonblocking
       ! sends to avoid deadlock.

       do nn=1,P3D_nProcSend

         ! Determine the size of the message and store the processor
         ! ID of the target and the starting adress a bit easier.

         procID      = P3D_procSend(nn)
         sizeMessage = 5*(P3D_sendSize(nn) - P3D_sendSize(nn-1))
         ii          = P3D_sendSize(nn-1) + 1

         call mpi_isend(sendBuf(1,ii), sizeMessage, sumb_integer, &
                        procID, procID, SUmb_comm_world,          &
                        sendRequests(nn), ierr)
       enddo
!
!      ******************************************************************
!      *                                                                *
!      * Post the nonblocking receives.                                 *
!      *                                                                *
!      ******************************************************************
!
       ! Allocate the memory for the receive buffer.

       ii = P3D_recvSize(P3D_nProcRecv)
       allocate(recvBuf(5,ii), stat=ierr)
       if(ierr /= 0)                         &
         call terminate("prepareReadPlot3D", &
                        "Memory allocation failure for recvBuf")

       ! Post the receives.

       do nn=1,P3D_nProcRecv

         ! Determine the size of the message and store the processor
         ! ID of the sender and the starting adress a bit easier.

         procID      = P3D_procRecv(nn)
         sizeMessage = 5*(P3D_recvSize(nn) - P3D_recvSize(nn-1))
         ii          = P3D_recvSize(nn-1) + 1

         call mpi_irecv(recvBuf(1,ii), sizeMessage, sumb_integer, &
                        procID, myID, SUmb_comm_world,            &
                        recvRequests(nn), ierr)
       enddo
!
!      ******************************************************************
!      *                                                                *
!      * Determine the variables to perform the local IO.               *
!      * Again a routine can be used which is almost identical for      *
!      * reading and writing.                                           *
!      *                                                                *
!      ******************************************************************
!
       call determineVarLocalIO(nGlobalBlocks,  nChunks,    &
                                sizeP3D_Real,   myChunk,    &
                                globalBlockIDs, startBlock, &
                                blockDimCGNS,   cgnsOffset, &
                                chunks,         .true.,     &
                                readHalos)
!
!      ******************************************************************
!      *                                                                *
!      * The receiving of the messages with the indices and copying of  *
!      * that data in P3D_commPart.                                     *
!      *                                                                *
!      ******************************************************************
!
       ! Complete the nonblocking sends and receives. The variable
       ! sizeMessage is used as a dummy, because it is an integer.

       sizeMessage = P3D_nProcSend
       do nn=1,P3D_nProcSend
         call mpi_waitany(sizeMessage, sendRequests, procID, status, ierr)
       enddo

       sizeMessage = P3D_nProcRecv
       do nn=1,P3D_nProcRecv
         call mpi_waitany(sizeMessage, recvRequests, procID, status, ierr)
       enddo

       ! Release the memory of the send buffer. All the nonblocking sends
       ! have been completed such that it is not needed anymore.

       deallocate(sendBuf, stat=ierr)
       if(ierr /= 0)                         &
         call terminate("prepareReadPlot3D", &
                        "Deallocation failure for sendBuf")

       ! Set the values of nItemsLocal and nItemsNonLocal. Also set
       ! offsetIO and nItemsTotal, although these are not used.

       P3D_commPart%nItemsLocal    = P3D_recvSize(P3D_nProcRecv)
       P3D_commPart%nItemsNonLocal = 0

       P3D_commPart%offsetIO    = 0
       P3D_commPart%nItemsTotal = P3D_commPart%nItemsLocal

       ! Allocate the memory for the block ID's and the indices needed
       ! from the mapping from the receive buffer to IOVar.

       ii = P3D_commPart%nItemsLocal
       allocate(P3D_commPart%blockID(ii),  &
                P3D_commPart%indexW(ii),   &
                P3D_commPart%posLocal(ii), &
                P3D_commPart%indices(ii,3), stat=ierr)
       if(ierr /= 0)                         &
         call terminate("prepareReadPlot3D", &
                        "Memory allocation failure for blockID, etc.")

       ! Copy the data from the receive buffer into P3D_commPart.

       do nn=1,P3D_commPart%nItemsLocal
         P3D_commPart%indices(nn,1) = recvBuf(1,nn)
         P3D_commPart%indices(nn,2) = recvBuf(2,nn)
         P3D_commPart%indices(nn,3) = recvBuf(3,nn)

         P3D_commPart%indexW(nn)  = recvBuf(4,nn)
         P3D_commPart%blockID(nn) = recvBuf(5,nn)

         P3D_commPart%posLocal(nn) = nn
       enddo

       ! Release the memory of the receive buffer.

       deallocate(recvBuf, stat=ierr)
       if(ierr /= 0)                         &
         call terminate("prepareReadPlot3D", &
                        "Deallocation failure for recvBuf")

       ! Correct P3D_sendSize and P3D_recvSize such that it corresponds
       ! to the size in bytes to be communicated later on.

       do nn=0,P3D_nProcSend
         P3D_sendSize(nn) = P3D_sendSize(nn)*sizeP3D_Real
       enddo

       do nn=0,P3D_nProcRecv
         P3D_recvSize(nn) = P3D_recvSize(nn)*sizeP3D_Real
       enddo

       end subroutine prepareReadPlot3D

!
!      ******************************************************************
!      *                                                                *
!      * File:          prepareWritePlot3D.f90                          *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 10-26-2005                                      *
!      * Last modified: 09-13-2007                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine prepareWritePlot3D(nGlobalBlocks, globalBlockIDs)
!
!      ******************************************************************
!      *                                                                *
!      * prepareWritePlot3D sets the Plot3D variables in IOModule such  *
!      * that the writing of the actual data to a Plot3D file can be    *
!      * performed. It is possible to write a selected number of global *
!      * indicated by the two subroutine arguments.                     *
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
       integer(kind=intType), intent(in) :: nGlobalBlocks
       integer(kind=intType), dimension(nGlobalBlocks), intent(in) :: &
                                                         globalBlockIDs
!
!      Local variables.
!
       integer :: sizeMessage, procID, ierr

       integer, dimension(mpi_status_size) :: status

       integer(kind=intType) :: ii, jj, mm, nn
       integer(kind=intType) :: nChunks, myChunk
       integer(kind=intType) :: sizeP3D_Real

       integer(kind=intType), dimension(nProc)     :: nPartsPerChunk
       integer(kind=intType), dimension(nProc)     :: startBlock
       integer(kind=intType), dimension(nProc)     :: nBytesCommChunk
       integer(kind=intType), dimension(0:nProc-1) :: recvData
       integer(kind=intType), dimension(0:nProc-1) :: sendData

       integer(kind=intType), dimension(3,nGlobalBlocks) :: blockDimCGNS

       integer(kind=intType), dimension(:,:), allocatable :: sendBuf
       integer(kind=intType), dimension(:,:), allocatable :: recvBuf

       integer(kind=mpi_offset_kind), dimension(0:nGlobalBlocks) :: &
                                                              cgnsOffset
       integer(kind=mpi_offset_kind), dimension(0:nProc) :: chunks

       logical :: IParticipateInIO
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Determine the size of the floating point type, depending
       ! on that floating point type.

       select case (P3D_Precision)
         case (precisionSingle)
           sizeP3D_Real = 4

         case (precisionDouble)
           sizeP3D_Real = 8
       end select

       ! Determine for every chunk to be written the number of bytes I
       ! would have to communicate to other processors.

       call determineBytesCommChunk(nGlobalBlocks,   sizeP3D_Real,    &
                                    nChunks,         globalBlockIDs,  &
                                    nPartsPerChunk,  startBlock,      &
                                    nBytesCommChunk, blockDimCGNS,    &
                                    cgnsOffset,      chunks,          &
                                    .false.,         .false.,         &
                                    IParticipateInIO)

       ! Determine the chunk for which I am responsible.

       call determineMyChunk(nBytesCommChunk,  nChunks, &
                             IParticipateInIO, myChunk)

       ! Determine the number of variables for the local IO as well as
       ! the number to be received from other processors.

       call determineCommSizesIO(nGlobalBlocks,  nChunks,        &
                                 sizeP3D_Real,   myChunk,        &
                                 globalBlockIDs, nPartsPerChunk, &
                                 startBlock,     blockDimCGNS,   &
                                 cgnsOffset,     chunks,         &
                                 .false.,        .false.)
!
!      ******************************************************************
!      *                                                                *
!      * Determine the number of processors, the corresponding          *
!      * processor ID's and the size of the messages which I have to    *
!      * send.                                                          *
!      *                                                                *
!      ******************************************************************
!
       ! Set the entries of the processors from which I receive data.

       recvData = 0
       do nn=1,P3D_nProcRecv
         recvData(P3D_procRecv(nn)) = P3D_recvSize(nn) &
                                    - P3D_recvSize(nn-1)
       enddo

       ! Perform an all to all communication, such that the processors
       ! know how many messages will be send as well as their size.

       call mpi_alltoall(recvData, 1, sumb_integer, &
                         sendData, 1, sumb_integer, SUmb_comm_world, ierr)

       ! Determine the number of processors to which I send data.

       P3D_nProcSend = 0
       do nn=0,(nProc-1)
         if(sendData(nn) > 0) P3D_nProcSend = P3D_nProcSend + 1
       enddo

       ! Allocate the memory for P3D_procSend and P3D_sendSize.

       allocate(P3D_procSend(P3D_nProcSend), &
                P3D_sendSize(0:P3D_nProcSend), stat=ierr)
       if(ierr /= 0)                          &
         call terminate("prepareWritePlot3D", &
                        "Memory allocation failure for P3D_procSend &
                        &and P3D_sendSize")

       ! Determine the values of these arrays. Note that P3D_sendSize
       ! is in cumulative storage format.

       ii = 0
       P3D_sendSize(0) = 0

       do nn=0,(nProc-1)
         if(sendData(nn) > 0) then
           ii = ii + 1
           P3D_procSend(ii) = nn
           P3D_sendSize(ii) = P3D_sendSize(ii-1) + sendData(nn)
         endif
       enddo
!
!      ******************************************************************
!      *                                                                *
!      * The sending of the index data to the processors that will      *
!      * send data to me during the actual IO.                          *
!      *                                                                *
!      ******************************************************************
!
       ! Allocate the memory for the send buffer.

       ii = P3D_recvSize(P3D_nProcRecv)
       allocate(sendBuf(5,ii), stat=ierr)
       if(ierr /= 0)                          &
         call terminate("prepareWritePlot3D", &
                        "Memory allocation failure for sendBuf")

       ! The routine which fill the actual send buffer.

       call sendBufNonLocalIO(nGlobalBlocks,  nChunks,       &
                              sizeP3D_Real,   myChunk,       &
                              globalBlockIDs, startBlock,    &
                              blockDimCGNS,   cgnsOffset,    &
                              chunks,         P3D_nProcRecv, &
                              P3D_procRecv,   P3D_recvSize,  &
                              sendBuf,        .false.,       &
                              .false.)

       ! Send the messages to the receiving processors. Use nonblocking
       ! sends to avoid deadlock. It may be a bit confusing that 
       ! P3D_nProcRecv etc is used for sending the data, but these names
       ! are chosen for the case when the actual IO takes place. And then
       ! they are used to receive data. Here only the meta-data is built
       ! to facilitate the parallel IO

       do nn=1,P3D_nProcRecv

         ! Determine the size of the message and store the processor
         ! ID of the target and the starting adress a bit easier.

         procID      = P3D_procRecv(nn)
         sizeMessage = 5*(P3D_recvSize(nn) - P3D_recvSize(nn-1))
         ii          = P3D_recvSize(nn-1) + 1

         call mpi_isend(sendBuf(1,ii), sizeMessage, sumb_integer, &
                        procID, procID, SUmb_comm_world,          &
                        sendRequests(nn), ierr)
       enddo
!
!      ******************************************************************
!      *                                                                *
!      * Post the nonblocking receives. Now of course P3D_nProcSend,    *
!      * etc. must be used. See the comments for the sending.           *
!      *                                                                *
!      ******************************************************************
!
       ! Allocate the memory for the receive buffer.

       ii = P3D_sendSize(P3D_nProcSend)
       allocate(recvBuf(5,ii), stat=ierr)
       if(ierr /= 0)                          &
         call terminate("prepareWritePlot3D", &
                        "Memory allocation failure for recvBuf")

       ! Post the receives.

       do nn=1,P3D_nProcSend

         ! Determine the size of the message and store the processor
         ! ID of the sender and the starting adress a bit easier.

         procID      = P3D_procSend(nn)
         sizeMessage = 5*(P3D_sendSize(nn) - P3D_sendSize(nn-1))
         ii          = P3D_sendSize(nn-1) + 1

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
                                chunks,         .false.,    &
                                .false.)
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
       ! Again the P3D_nProcSend correspond to the number of receiving
       ! processors and P3D_nProcRecv to the number of sending ones.

       sizeMessage = P3D_nProcRecv
       do nn=1,P3D_nProcRecv
         call mpi_waitany(sizeMessage, sendRequests, procID, status, ierr)
       enddo

       sizeMessage = P3D_nProcSend
       do nn=1,P3D_nProcSend
         call mpi_waitany(sizeMessage, recvRequests, procID, status, ierr)
       enddo

       ! Release the memory of the send buffer. All the nonblocking sends
       ! have been completed such that it is not needed anymore.

       deallocate(sendBuf, stat=ierr)
       if(ierr /= 0)                          &
         call terminate("prepareWritePlot3D", &
                        "Deallocation failure for sendBuf")

       ! Set the values of nItemsLocal and nItemsNonLocal. Also set
       ! offsetIO and nItemsTotal, although these are not used.

       P3D_commPart%nItemsLocal    = P3D_sendSize(P3D_nProcSend)
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
       if(ierr /= 0)                          &
         call terminate("prepareWritePlot3D", &
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
       if(ierr /= 0)                          &
         call terminate("prepareWritePlot3D", &
                        "Deallocation failure for recvBuf")
!
!      ******************************************************************
!      *                                                                *
!      * The variables needed to write the record integers at the       *
!      * appropriate positions.                                         *
!      *                                                                *
!      ******************************************************************
!
       ! Store the starting and ending index in globalBlockIDs in mm
       ! and nn respectively.

       mm = startBlock(myChunk)
       nn = mm + P3D_nIOParts - 1

       ! Determine the number of record integers I have to write.

       P3D_nRecordIntegersWrite = 2*(P3D_nIOParts - 1)

       if(chunks(myChunk-1) == cgnsOffset(mm-1)) &
         P3D_nRecordIntegersWrite = P3D_nRecordIntegersWrite + 1

       if(chunks(myChunk) == cgnsOffset(nn)) &
         P3D_nRecordIntegersWrite = P3D_nRecordIntegersWrite + 1

       ! Allocate the memory for P3D_recordIntegersWrite and
       ! P3D_recordPosition.

       allocate(P3D_recordIntegersWrite(P3D_nRecordIntegersWrite), &
                P3D_recordPosition(P3D_nRecordIntegersWrite), stat=ierr)
       if(ierr /= 0)                          &
         call terminate("prepareWritePlot3D", &
                        "Memory allocation failure for &
                        &P3D_recordIntegersWrite and P3D_recordPosition")

       ! Determine the values of P3D_recordIntegersWrite and
       ! P3D_recordPosition. The former contains the record size in bytes
       ! and the latter the relative offset in bytes to the beginning
       ! of my part of the write buffer.

       ! First check if I have to write the opening record of the
       ! first block part.

       ii = 0
       if(chunks(myChunk-1) == cgnsOffset(mm-1)) then
         ii = ii + 1
         P3D_recordIntegersWrite(ii) = cgnsOffset(mm) - cgnsOffset(mm-1) &
                                     - 2*nBytesPerRecordIntPLOT3D
         P3D_recordPosition(ii)      = 0
       endif

       ! Loop over the internal record sizes.

       do jj=2,P3D_nIOParts

         ! The ending record of global block mm.

         ii = ii + 1
         P3D_recordIntegersWrite(ii) = cgnsOffset(mm) - cgnsOffset(mm-1) &
                                     - 2*nBytesPerRecordIntPLOT3D
         P3D_recordPosition(ii)      = cgnsOffset(mm) - chunks(myChunk-1) &
                                     - nBytesPerRecordIntPLOT3D

         ! Update mm afterwards for the next global block.

         mm = mm + 1

         ! The starting record of the next global block, which is now
         ! block mm.

         ii = ii + 1
         P3D_recordIntegersWrite(ii) = cgnsOffset(mm)   - cgnsOffset(mm-1) &
                                     - 2*nBytesPerRecordIntPLOT3D
         P3D_recordPosition(ii)      = cgnsOffset(mm-1) - chunks(myChunk-1)

       enddo

       ! Finally check if I have to write the closing record integer
       ! of the last block part.

       if(chunks(myChunk) == cgnsOffset(nn)) then
         ii = ii + 1
         P3D_recordIntegersWrite(ii) = cgnsOffset(nn) - cgnsOffset(nn-1) &
                                     - 2*nBytesPerRecordIntPLOT3D
         P3D_recordPosition(ii)      = cgnsOffset(nn) - chunks(myChunk-1) &
                                     - nBytesPerRecordIntPLOT3D
       endif

       ! Correct P3D_sendSize and P3D_recvSize such that it corresponds
       ! to the size in bytes to be communicated later on.

       do nn=0,P3D_nProcSend
         P3D_sendSize(nn) = P3D_sendSize(nn)*sizeP3D_Real
       enddo

       do nn=0,P3D_nProcRecv
         P3D_recvSize(nn) = P3D_recvSize(nn)*sizeP3D_Real
       enddo

       end subroutine prepareWritePlot3D

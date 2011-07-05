!
!      ******************************************************************
!      *                                                                *
!      * File:          writePlot3DVar.f90                              *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 10-19-2005                                      *
!      * Last modified: 11-02-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine writePlot3DVar(posIO, fact, posW)
!
!      ******************************************************************
!      *                                                                *
!      * writePlot3DVar writes the w variable of IOVar(:,posIO) to a    *
!      * Plot3D file. The starting position in w is posW and P3D_nVar   *
!      * variables are written; P3D_nVar is defined in IOModule. Every  *
!      * processor is responsible for a contiguous chunk of data, which *
!      * is received from the appropriate processor.                    *
!      *                                                                *
!      ******************************************************************
!
       use communication
       use inputIO
       use IOModule
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType), intent(in) :: posIO, posW
       real(kind=realType),   intent(in) :: fact
!
!      Local variables.
!
       integer :: sizeBuf, proc, ierr
       integer, dimension(mpi_status_size) :: status

       integer(kind=intType) :: nn, mm, ii, i
       integer(kind=intType) :: sizeP3D_Real, sizeP3D_Record

       integer(kind=intRecordPLOT3DType) :: sizeRecord

       integer(kind=mpi_offset_kind) :: disp

       character, dimension(:), allocatable :: writeBuf
       character, dimension(:), allocatable :: sendBuf
       character, dimension(:), allocatable :: recvBuf
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

       ! Copy the size of the integer record type into sizeP3D_Record.
       ! The reason is that the variable nBytesPerRecordIntPLOT3D is
       ! of type integer and not integer(kind=intType).

       sizeP3D_Record = nBytesPerRecordIntPLOT3D
!
!      ******************************************************************
!      *                                                                *
!      * Copy the local data in the send buffer and post the            *
!      * non-blocking communication.                                    *
!      *                                                                *
!      ******************************************************************
!
       ! Allocate the memory for the write buffer, the send buffer 
       ! and the receive buffer.

       ii = P3D_mySizeIO
       nn = P3D_sendSize(P3D_nProcSend)
       mm = P3D_recvSize(P3D_nProcRecv)

       allocate(writeBuf(ii), sendBuf(nn), recvBuf(mm), stat=ierr)
       if(ierr /= 0)                      &
         call terminate("writePlot3DVar", &
                        "Memory allocation failure for readBuf, &
                        &sendBuf and recvBuf")

       ! Copy the data from IOVar to the send buffer.
       ! The routine called depends on the floating point type.
       ! The variable mm is a dummy and should not be used.

       select case (P3D_Precision)
         case (precisionSingle)
           call copyIOPartWriteSP(posIO, fact, posW, sendBuf, &
                                  mm, P3D_commPart, .true.)

         case (precisionDouble)
           call copyIOPartWriteDP(posIO, fact, posW, sendBuf, &
                                  mm, P3D_commPart, .true.)
       end select

       ! Send the data I have to send. Use nonblocking sends to avoid
       ! deadlock.

       sendLoop: do nn=1,P3D_nProcSend

         ! Determine the size of the message and its offset in the
         ! send buffer. Send the message afterwards.

         ii      = P3D_sendSize(nn-1) + 1
         sizeBuf = P3D_sendSize(nn) - P3D_sendSize(nn-1)
         proc    = P3D_procSend(nn)

         call mpi_isend(sendBuf(ii), sizeBuf, mpi_character, proc, &
                        proc, SUmb_comm_world, sendRequests(nn), ierr)
       enddo sendLoop

       ! Post the nonblocking receives.

       recvLoop: do nn=1,P3D_nProcRecv

         ! Determine the size of the message and its offset in the
         ! receive buffer. Post the nonblocking receive afterwards.

         ii      = P3D_recvSize(nn-1) + 1
         sizeBuf = P3D_recvSize(nn) - P3D_recvSize(nn-1)
         proc    = P3D_procRecv(nn)

         call mpi_irecv(recvBuf(ii), sizeBuf, mpi_character, proc, &
                        myID, SUmb_comm_world, recvRequests(nn), ierr)
       enddo recvLoop
!
!      ******************************************************************
!      *                                                                *
!      * Copy the local data from IOVar. This is done after the sending *
!      * of the messages to hide the communication as much as possible. *
!      *                                                                *
!      ******************************************************************
!
       ! Loop over the number of IO parts, which correspond to a
       ! global block (part).

       localCopyLoop: do nn=1,P3D_nIOParts

         ! Abbreviate the offset for this part a bit easier.
         ! Use offsetBuffer for alignment reasons. A possible shift is
         ! applied after the nonlocal data has been copied into the
         ! write buffer.

         ii = P3D_IOParts(nn)%offsetBuffer + 1

         ! Copy the local data directly from IOVar in the write buffer.
         ! The routine called depends on the floating point type.
         ! The variable mm is a dummy and should not be used.

         select case (P3D_Precision)
           case (precisionSingle)
             call copyIOPartWriteSP(posIO, fact, posW, writeBuf(ii), &
                                    mm, P3D_IOParts(nn), .true.)

           case (precisionDouble)
             call copyIOPartWriteDP(posIO, fact, posW, writeBuf(ii), &
                                    mm, P3D_IOParts(nn), .true.)
         end select

       enddo localCopyLoop
!
!      ******************************************************************
!      *                                                                *
!      * Copy the data from the receive buffer into the write buffer.   *
!      *                                                                *
!      ******************************************************************
!
       ! Complete the non blocking receives. Use sizeBuf as a help 
       ! variable, because this is an integer variable.

       sizeBuf = P3D_nProcRecv
       do nn=1,P3D_nProcRecv
         call mpi_waitany(sizeBuf, recvRequests, proc, status, ierr)
       enddo

       ! Loop over the number of IO parts, which correspond to a global
       ! block (part), to copy the data from the receive buffer into
       ! the write buffer. In this way it is easy to filter out the
       ! integers of the records.

       recvBufferLoop: do nn=1,P3D_nIOParts

         ! Abbreviate the offset for this part a bit easier.
         ! Use offsetBuffer for alignment reasons.

         ii = P3D_IOParts(nn)%offsetBuffer + 1

         ! Copy the nonlocal data to from the receive the write buffer.
         ! The routine called depends on the floating point type.

         select case (P3D_Precision)
           case (precisionSingle)
             call copyIOPartWriteSP(posIO, fact, posW, writeBuf(ii), &
                                    recvBuf, P3D_IOParts(nn), .false.)

           case (precisionDouble)
             call copyIOPartWriteDP(posIO, fact, posW, writeBuf(ii), &
                                    recvBuf, P3D_IOParts(nn), .false.)
         end select

         ! Apply the byte swapping for the floating point data if needed.

         if( P3D_ByteSwap ) &
           call byteswap(writeBuf(ii), sizeP3D_Real, &
                         P3D_IOParts(nn)%nItemsTotal)

         ! Apply the shift if the character array is not aligned
         ! properly. As a shift to the right is needed here, the
         ! loop starts at the back.

         if(P3D_IOParts(nn)%offsetIO > P3D_IOParts(nn)%offsetBuffer) then

           mm = P3D_IOParts(nn)%offsetIO - P3D_IOParts(nn)%offsetBuffer
           ii = P3D_IOParts(nn)%offsetBuffer &
              + sizeP3D_Real*P3D_IOParts(nn)%nItemsTotal

           do i=ii,(P3D_IOParts(nn)%offsetBuffer+1),-1
             writeBuf(i+mm) = writeBuf(i)
           enddo

         endif

       enddo recvBufferLoop

       ! Write the integers of the record sizes for which I am
       ! responsible. This must be done after the possible shift
       ! is applied to avoid overwriting of data.

       do nn=1,P3D_nRecordIntegersWrite
         sizeRecord = P3D_recordIntegersWrite(nn)
         if( P3D_ByteSwap ) &
           call byteswap(sizeRecord, sizeP3D_Record, 1_intType)

         ii = P3D_recordPosition(nn) + 1
         call writeVarsToBuffer(writeBuf(ii), sizeRecord, &
                                sizeP3D_Record)
       enddo

       ! Complete the nonblocking sends. Use sizeBuf as a help variable,
       ! because this is an integer variable.

       sizeBuf = P3D_nProcSend
       do nn=1,P3D_nProcSend
         call mpi_waitany(sizeBuf, sendRequests, proc, status, ierr)
       enddo

       ! Release the memory of the send and receive buffer.

       deallocate(sendBuf, recvBuf, stat=ierr)
       if(ierr /= 0)                      &
         call terminate("writePlot3DVar", &
                        "Deallocation failure for sendBuf and recvBuf")
!
!      ******************************************************************
!      *                                                                *
!      * Write the data.                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Jump to the correct position in the file and write the data.
       ! This is only done if I have to write something.

       if(P3D_mySizeIO > 0) then
         disp = P3D_Offset + P3D_myOffset
         call mpi_file_seek(fh, disp, mpi_seek_set, ierr)
         call mpi_file_write(fh, writeBuf, P3D_mySizeIO, &
                             mpi_character, status, ierr)
       endif

       ! Deallocate the memory of writeBuf.

       deallocate(writeBuf, stat=ierr)
       if(ierr /= 0)                      &
         call terminate("writePlot3DVar", &
                        "Deallocation failure for writeBuf")

       end subroutine writePlot3DVar

       !=================================================================

       subroutine copyIOPartWriteSP(posIO, fact, posW, writeBuf, &
                                    recvBuf, P3D_IOPart, localCopy)
!
!      ******************************************************************
!      *                                                                *
!      * copyIOPartWriteSP either copies the local IOVar array or       *
!      * recvBuf to writeBuf, depending on the value of the logical     *
!      * localCopy.                                                     *
!      *                                                                *
!      ******************************************************************
!
       use IOModule
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType), intent(in) :: posIO, posW
       real(kind=realType),   intent(in) :: fact

       real(kind=4), dimension(*), intent(in)    :: recvBuf
       real(kind=4), dimension(*), intent(inout) :: writeBuf

       logical, intent(in) :: localCopy

       type(P3D_IOPartType), intent(in) :: P3D_IOPart
!
!      Local variables.
!
       integer(kind=intType) :: nn, mm, i, j, k, l, m
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Check if a local copy or a copy from the receive buffer must
       ! be done.

       testLocalCopy: if( localCopy ) then

         ! Copy the local data in the correct position in IOVar.
         ! Take the scaling factor into account.

         do nn=1,P3D_IOPart%nItemsLocal
           mm = P3D_IOPart%blockID(nn)
           i  = P3D_IOPart%indices(nn,1) + IOVar(mm,posIO)%pointerOffset
           j  = P3D_IOPart%indices(nn,2) + IOVar(mm,posIO)%pointerOffset
           k  = P3D_IOPart%indices(nn,3) + IOVar(mm,posIO)%pointerOffset
           l  = P3D_IOPart%indexW(nn) + posW
           m  = P3D_IOPart%posLocal(nn)

           writeBuf(m) = fact*IOVar(mm,posIO)%w(i,j,k,l)
         enddo

       else testLocalCopy

         ! Copy the nonlocal data in the correct position in writeBuf.
         ! The scaling factor is already taken into account on the
         ! sending side.

         do nn=1,P3D_IOPart%nItemsNonLocal
           l = P3D_IOPart%posComm(nn)
           m = P3D_IOPart%posNonLocal(nn)

           writeBuf(m) = recvBuf(l)
         enddo

       endif testLocalCopy

       end subroutine copyIOPartWriteSP

       !=================================================================

       subroutine copyIOPartWriteDP(posIO, fact, posW, writeBuf, &
                                    recvBuf, P3D_IOPart, localCopy)
!
!      ******************************************************************
!      *                                                                *
!      * copyIOPartWriteDP either copies the local IOVar array or       *
!      * recvBuf to writeBuf, depending on the value of the logical     *
!      * localCopy.                                                     *
!      *                                                                *
!      ******************************************************************
!
       use IOModule
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType), intent(in) :: posIO, posW
       real(kind=realType),   intent(in) :: fact

       real(kind=8), dimension(*), intent(in)    :: recvBuf
       real(kind=8), dimension(*), intent(inout) :: writeBuf

       logical, intent(in) :: localCopy

       type(P3D_IOPartType), intent(in) :: P3D_IOPart
!
!      Local variables.
!
       integer(kind=intType) :: nn, mm, i, j, k, l, m
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Check if a local copy or a copy from the receive buffer must
       ! be done.

       testLocalCopy: if( localCopy ) then

         ! Copy the local data in the correct position in IOVar.
         ! Take the scaling factor into account.

         do nn=1,P3D_IOPart%nItemsLocal
           mm = P3D_IOPart%blockID(nn)
           i  = P3D_IOPart%indices(nn,1) + IOVar(mm,posIO)%pointerOffset
           j  = P3D_IOPart%indices(nn,2) + IOVar(mm,posIO)%pointerOffset
           k  = P3D_IOPart%indices(nn,3) + IOVar(mm,posIO)%pointerOffset
           l  = P3D_IOPart%indexW(nn) + posW
           m  = P3D_IOPart%posLocal(nn)

           writeBuf(m) = fact*IOVar(mm,posIO)%w(i,j,k,l)
         enddo

       else testLocalCopy

         ! Copy the nonlocal data in the correct position in writeBuf.
         ! The scaling factor is already taken into account on the
         ! sending side.

         do nn=1,P3D_IOPart%nItemsNonLocal
           l = P3D_IOPart%posComm(nn)
           m = P3D_IOPart%posNonLocal(nn)

           writeBuf(m) = recvBuf(l)
         enddo

       endif testLocalCopy

       end subroutine copyIOPartWriteDP

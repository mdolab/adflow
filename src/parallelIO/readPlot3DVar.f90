!
!      ******************************************************************
!      *                                                                *
!      * File:          readPlot3DVar.f90                               *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 10-18-2005                                      *
!      * Last modified: 11-02-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine readPlot3DVar(posIO, fact, posW)
!
!      ******************************************************************
!      *                                                                *
!      * readPlot3DVar reads a part of a Plot3D file and puts the data  *
!      * in the w variable of IOVar(:,posIO). The starting position in  *
!      * w is posW and P3D_nVar variables are read; P3D_nVar is defined *
!      * in IOModule. Every processor is responsible for a contiguous   *
!      * chunk of data, which is sent to the appropriate processor.     *
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
       integer(kind=intType) :: sizeP3D_Real

       integer(kind=mpi_offset_kind) :: disp

       character, dimension(:), allocatable :: readBuf
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
!
!      ******************************************************************
!      *                                                                *
!      * Read my buffer, prepare the send buffer and post the           *
!      * non-blocking communication.                                    *
!      *                                                                *
!      ******************************************************************
!
       ! Allocate the memory for the read buffer, the send buffer 
       ! and the receive buffer.

       ii = P3D_mySizeIO
       nn = P3D_sendSize(P3D_nProcSend)
       mm = P3D_recvSize(P3D_nProcRecv)

       allocate(readBuf(ii), sendBuf(nn), recvBuf(mm), stat=ierr)
       if(ierr /= 0)                     &
         call terminate("readPlot3DVar", &
                        "Memory allocation failure for readBuf, &
                        &sendBuf and recvBuf")

       ! Jump to the correct position in the file and read the data.
       ! This is only done if I have to read something.

       if(P3D_mySizeIO > 0) then
         disp = P3D_Offset + P3D_myOffset
         call mpi_file_seek(fh, disp, mpi_seek_set, ierr)
         call mpi_file_read(fh, readBuf, P3D_mySizeIO, mpi_character, &
                            status, ierr)
       endif

       ! Loop over the number of IO parts, which correspond to a global
       ! block (part), to fill the send buffer. In this way it is easy
       ! to filter out the integers of the records and possible iblanks.

       sendBufferLoop: do nn=1,P3D_nIOParts

         ! Apply the shift if the character array is not aligned
         ! properly. On some platforms this causes a problem.

         if(P3D_IOParts(nn)%offsetIO > P3D_IOParts(nn)%offsetBuffer) then

           mm = P3D_IOParts(nn)%offsetIO - P3D_IOParts(nn)%offsetBuffer
           ii = P3D_IOParts(nn)%offsetBuffer &
              + sizeP3D_Real*P3D_IOParts(nn)%nItemsTotal

           do i=(P3D_IOParts(nn)%offsetBuffer+1),ii
             readBuf(i) = readBuf(i+mm)
           enddo

         endif

         ! Abbreviate the offset for this part a bit easier and apply
         ! the byte swapping if needed.

         ii = P3D_IOParts(nn)%offsetBuffer + 1

         if( P3D_ByteSwap )                         &
           call byteswap(readBuf(ii), sizeP3D_Real, &
                         P3D_IOParts(nn)%nItemsTotal)

         ! Copy the nonlocal data to the send buffer.
         ! The routine called depends on the floating point type.

         select case (P3D_Precision)
           case (precisionSingle)
             call copyIOPartReadSP(posIO, fact, posW, readBuf(ii), &
                                   sendBuf, P3D_IOParts(nn), .false.)

           case (precisionDouble)
             call copyIOPartReadDP(posIO, fact, posW, readBuf(ii), &
                                   sendBuf, P3D_IOParts(nn), .false.)
         end select

       enddo sendBufferLoop

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
!      * Copy the local data in IOVar. This is done after the sending   *
!      * of the messages to hide the communication as much as possible. *
!      *                                                                *
!      ******************************************************************
!
       ! Again loop over the number of IO parts, which correspond to a
       ! global block (part).

       localCopyLoop: do nn=1,P3D_nIOParts

         ! Abbreviate the offset for this part a bit easier.
         ! Note that the possible shift has already been applied in
         ! sendBufferLoop.

         ii = P3D_IOParts(nn)%offsetBuffer + 1

         ! Copy the local data directly from the read buffer into IOVar.
         ! The routine called depends on the floating point type.
         ! The variable mm is a dummy and should not be used.

         select case (P3D_Precision)
           case (precisionSingle)
             call copyIOPartReadSP(posIO, fact, posW, readBuf(ii), &
                                   mm, P3D_IOParts(nn), .true.)

           case (precisionDouble)
             call copyIOPartReadDP(posIO, fact, posW, readBuf(ii), &
                                   mm, P3D_IOParts(nn), .true.)
         end select

       enddo localCopyLoop

       ! Release the memory of the read buffer.

       deallocate(readBuf, stat=ierr)
       if(ierr /= 0)                     &
         call terminate("readPlot3DVar", &
                        "Deallocation failure for readBuf")
!
!      ******************************************************************
!      *                                                                *
!      * Copy the received data in IOVar.                               *
!      *                                                                *
!      ******************************************************************
!
       ! Complete the nonblocking receives. Use sizeBuf as a help 
       ! variable, because this is an integer variable.

       sizeBuf = P3D_nProcRecv
       do nn=1,P3D_nProcRecv
         call mpi_waitany(sizeBuf, recvRequests, proc, status, ierr)
       enddo

       ! Copy the data from the receive buffer into IOVar.
       ! The routine called depends on the floating point type.
       ! The variable mm is a dummy and should not be used.

       select case (P3D_Precision)
         case (precisionSingle)
           call copyIOPartReadSP(posIO, fact, posW, recvBuf, &
                                 mm, P3D_commPart, .true.)

         case (precisionDouble)
           call copyIOPartReadDP(posIO, fact, posW, recvBuf, &
                                 mm, P3D_commPart, .true.)
       end select

       ! Complete the nonblocking sends. Use sizeBuf as a help variable,
       ! because this is an integer variable.

       sizeBuf = P3D_nProcSend
       do nn=1,P3D_nProcSend
         call mpi_waitany(sizeBuf, sendRequests, proc, status, ierr)
       enddo

       ! Release the memory of the send and receive buffer.

       deallocate(sendBuf, recvBuf, stat=ierr)
       if(ierr /= 0)                     &
         call terminate("readPlot3DVar", &
                        "Deallocation failure for sendBuf and recvBuf")

       end subroutine readPlot3DVar

       !=================================================================

       subroutine copyIOPartReadSP(posIO, fact, posW, readBuf, &
                                   sendBuf, P3D_IOPart, localCopy)
!
!      ******************************************************************
!      *                                                                *
!      * copyIOPartReadSP copies the read buffer either to the local    *
!      * IOVar array or to sendBuf, depending on the value of the       *
!      * logical localCopy.                                             *
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

       real(kind=4), dimension(*), intent(in)    :: readBuf
       real(kind=4), dimension(*), intent(inout) :: sendBuf

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
       ! Check if a local copy or a copy to the send buffer must
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

           IOVar(mm,posIO)%w(i,j,k,l) = fact*readBuf(m)
         enddo

       else testLocalCopy

         ! Copy the nonlocal data in the correct position in sendBuf.
         ! The scaling factor will be taken into account later on the
         ! receiving side.

         do nn=1,P3D_IOPart%nItemsNonLocal
           l = P3D_IOPart%posComm(nn)
           m = P3D_IOPart%posNonLocal(nn)

           sendBuf(l) = readBuf(m)
         enddo

       endif testLocalCopy

       end subroutine copyIOPartReadSP

       !=================================================================

       subroutine copyIOPartReadDP(posIO, fact, posW, readBuf, &
                                   sendBuf, P3D_IOPart, localCopy)
!
!      ******************************************************************
!      *                                                                *
!      * copyIOPartReadDP copies the read buffer either to the local    *
!      * IOVar array or to sendBuf, depending on the value of the       *
!      * logical localCopy.                                             *
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

       real(kind=8), dimension(*), intent(in)    :: readBuf
       real(kind=8), dimension(*), intent(inout) :: sendBuf

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
       ! Check if a local copy or a copy to the send buffer must
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

           IOVar(mm,posIO)%w(i,j,k,l) = fact*readBuf(m)
         enddo

       else testLocalCopy

         ! Copy the nonlocal data in the correct position in sendBuf.
         ! The scaling factor will be taken into account later on the
         ! receiving side.

         do nn=1,P3D_IOPart%nItemsNonLocal
           l = P3D_IOPart%posComm(nn)
           m = P3D_IOPart%posNonLocal(nn)

           sendBuf(l) = readBuf(m)
         enddo

       endif testLocalCopy

       end subroutine copyIOPartReadDP

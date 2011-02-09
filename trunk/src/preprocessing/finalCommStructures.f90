!
!      ******************************************************************
!      *                                                                *
!      * File:          finalCommStructures.f90                         *
!      * Author:        Edwin van der Weide, Steve Repsher              *
!      * Starting date: 02-20-2003                                      *
!      * Last modified: 06-12-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine finalCommStructures(entityHalo, nHalo, commPattern, &
                                      internalComm, nInterp)
!
!      ******************************************************************
!      *                                                                *
!      * FinalCommStructures determines the communication data          *
!      * structures used in the flow solver, commPattern and            *
!      * internalComm, from the given haloList, entityHalo.             *
!      *                                                                *
!      ******************************************************************
!
       use communication
       use haloList
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType), intent(in)            :: nHalo
       type(haloListType), dimension(*), intent(in) :: entityHalo

       type(commType),         intent(out) :: commPattern
       type(internalCommType), intent(out) :: internalComm

       integer(kind=intType), intent(in) :: nInterp
!
!      Local variables.
!
       integer :: ierr, count, dest, source, sizeRecv

       integer, dimension(mpi_status_size)  :: status

       integer(kind=intType) :: i, j
       integer(kind=intType) :: ii, jj, kk, ll, mm, nn, pp
       integer(kind=intType) :: sizeBuffer

       integer(kind=intType), dimension(:), allocatable :: nHaloPerProc
       integer(kind=intType), dimension(:), allocatable :: sendInfo

       integer(kind=intType), dimension(:,:), allocatable :: buffer
       integer(kind=intType), dimension(:,:), allocatable :: recvBuf

       real(kind=realType), dimension(:,:), allocatable :: bufInt
       real(kind=realType), dimension(:,:), allocatable :: recvBufInt
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Allocate the memory for nHaloPerProc and initialize its
       ! values to 0. Boundary halo's are stored with a processor id -1
       ! and the array will be put in cumulative storage format later,
       ! which explains the boundaries in the allocation.
       ! Also allocate the memory for sendInfo, which is needed in the
       ! all to all communication call.

       allocate(nHaloPerProc(0:nProc), sendInfo(0:nProc-1), stat=ierr)
       if(ierr /= 0)                             &
         call terminate("finalCommStructures", &
                        "Memory allocation error for nHaloPerProc &
                        &and sendInfo")

       nHaloPerProc = 0

       ! Loop over the number of halo's and determine the nHaloPerProc.
       ! Note that the boundary conditions are stored in entityHalo with
       ! a processor id -1.

       do i=1,nHalo
         ii = entityHalo(i)%donorProc + 1
         nHaloPerProc(ii) = nHaloPerProc(ii) + 1
       enddo

       ! Perform an all to all communication, such that the processors
       ! know how many messages will be received as well as their size.

       call mpi_alltoall(nHaloPerProc(1), 1,  sumb_integer, &
                         sendInfo(0),      1, sumb_integer, &
                         SUmb_comm_world, ierr)

       ! Allocate the memory for indexRecvProc, the index in the
       ! receive info for a certain processor. Initialize these value
       ! to 0, indicating that nothing is received from a processor.

       allocate(commPattern%indexRecvProc(0:nProc-1), stat=ierr)
       if(ierr /= 0)                             &
         call terminate("finalCommStructures", &
                        "Memory allocation error for indexRecvProc")

       do i=0,(nProc-1)
         commPattern%indexRecvProc(i) = 0
       enddo

       ! Determine the number of processors from which i receive data.
       ! Receive data here means when an actual exchange takes place.
       ! Furthermore determine the size of the buffer for the
       ! nonblocking sends.

       commPattern%nProcRecv = 0
       sizeBuffer             = 0
       do i=0,(nProc-1)
         ii = i+1
         if(nHaloPerProc(ii) > 0 .and. i /= myId) then
           commPattern%nProcRecv         = commPattern%nProcRecv + 1
           commPattern%indexRecvProc(i) = commPattern%nProcRecv

           sizeBuffer = sizeBuffer + nHaloPerProc(ii)
         endif
       enddo

       ! Allocate the memory for recvProc, nrecv, nrecvCum
       ! and recvList.

       ii = commPattern%nProcRecv
       allocate(commPattern%recvProc(ii),   &
                commPattern%nrecv(ii),       &
                commPattern%nrecvCum(0:ii), &
                commPattern%recvList(ii), stat=ierr)
       if(ierr /= 0)                             &
         call terminate("finalCommStructures", &
                        "Memory allocation error for recvProc, etc")

       ! Allocate memory for buffers, needed for the nonblocking sends.

       allocate(buffer(4,sizeBuffer), bufInt(ninterp,sizeBuffer), &
                stat=ierr)
       if(ierr /= 0)                             &
         call terminate("finalCommStructures", &
                        "Memory allocation error for buffer, bufInt.")

       ! Repeat the loop over the processors, but now store receive info.
       ! Jj is the counter for the current index in the receive processors
       ! and mm the counter in the buffer used to send the halo info.

       jj = 0
       mm = 1
       commPattern%nrecvCum(0) = 0
       do i=0,(nProc-1)
         ii = i+1
         if(nHaloPerProc(ii) > 0 .and. i /= myId) then

           ! Update the counter jj, store the number of halo's a bit
           ! easier in kk and set recvProc and nrecv.

           jj = jj + 1
           kk = nHaloPerProc(ii)

           commPattern%recvProc(jj) = i
           commPattern%nrecv(jj)     = kk
           commPattern%nrecvCum(jj) = commPattern%nrecvCum(jj-1) + kk

           ! Allocate the memory for the receive list.

           allocate(commPattern%recvList(jj)%block(kk),     &
                    commPattern%recvList(jj)%indices(kk,3), &
                    stat=ierr)
           if(ierr /= 0)                             &
             call terminate("finalCommStructures", &
                            "Memory allocation failure for block and &
                            &indices of recvlist")

           ! Copy the halo info in the receive list of commPattern.

           do j=1,kk
             ll = j + nHaloPerProc(i)

             commPattern%recvList(jj)%block(j) = entityHalo(ll)%myBlock

             commPattern%recvList(jj)%indices(j,1) = entityHalo(ll)%myI
             commPattern%recvList(jj)%indices(j,2) = entityHalo(ll)%myJ
             commPattern%recvList(jj)%indices(j,3) = entityHalo(ll)%myK
           enddo

           ! Copy the donor info in the buffer and send it to the
           ! appropriate processor.

           nn = mm - 1
           do j=1,kk
             ll = j + nHaloPerProc(i)
             nn = nn +1

             buffer(1,nn) = entityHalo(ll)%donorBlock
             buffer(2,nn) = entityHalo(ll)%dI
             buffer(3,nn) = entityHalo(ll)%dJ
             buffer(4,nn) = entityHalo(ll)%dK

             do pp = 1,ninterp
               bufInt(pp,nn) = entityHalo(ll)%interp(pp)
             end do
           enddo

           ! Copy some values to be sure integer data is used in the
           ! MPI call.

           count = 4*kk
           dest  = i

           ! And send the stuff.

           call mpi_isend(buffer(1,mm), count, sumb_integer, &
                          dest, dest+2, SUmb_comm_world,     &
                          sendRequests(jj), ierr)

           ! Now send the interpolants, if any are being exchanged.

           if(ninterp > 0) then
             count = nInterp*kk
             call mpi_isend(bufInt(1,mm), count, sumb_real, &
                            dest, dest+3, SUmb_comm_world,  &
                            recvRequests(jj), ierr)
           end if

           ! Update mm to the index in buffer for the next processor.

           mm = mm + kk
         endif

         ! Put nHaloPerProc in cumulative storage format.

         nHaloPerProc(ii) = nHaloPerProc(ii) + nHaloPerProc(i)
       enddo

       ! Determine the number of internal memory to memory copies and
       ! allocate the memory for it.

       ii = nHaloPerProc(myId+1) - nHaloPerProc(myId)
       internalComm%ncopy = ii

       allocate(internalComm%donorBlock(ii),          &
                internalComm%donorIndices(ii,3),      &
                internalComm%donorInterp(ii,ninterp), &
                internalComm%haloBlock(ii),           &
                internalComm%haloIndices(ii,3), stat=ierr)
       if(ierr /= 0)                             &
         call terminate("finalCommStructures", &
                        "Memory allocation failure for internalComm")

       ! Copy the info from the halo list.

       do i=1,internalComm%ncopy
         ii = i + nHaloPerProc(myId)

         internalComm%donorBlock(i) = entityHalo(ii)%donorBlock
         internalComm%haloBlock(i)  = entityHalo(ii)%myBlock

         internalComm%donorIndices(i,1) = entityHalo(ii)%dI
         internalComm%donorIndices(i,2) = entityHalo(ii)%dJ
         internalComm%donorIndices(i,3) = entityHalo(ii)%dK

         do pp = 1,ninterp
           internalComm%donorInterp(i,pp) = entityHalo(ii)%interp(pp)
         end do

         internalComm%haloIndices(i,1) = entityHalo(ii)%myI
         internalComm%haloIndices(i,2) = entityHalo(ii)%myJ
         internalComm%haloIndices(i,3) = entityHalo(ii)%myK
       enddo

       ! Determine from the earlier MPI_alltoall call the processors
       ! to whom i must send data in a normal data exchange. Also
       ! determine the amount i must send.

       ! First determine the number of processors to which i must send
       ! data and the size of the buffer for receiving data.
       ! Allocate and determine indexSendProc as well.

       allocate(commPattern%indexSendProc(0:nProc-1), stat=ierr)
       if(ierr /= 0)                             &
         call terminate("finalCommStructures", &
                        "Memory allocation error for indexSendProc")

       do i=0,(nProc-1)
         commPattern%indexSendProc(i) = 0
       enddo

       commPattern%nProcSend = 0
       sizeBuffer = 0

       do i=0,(nProc-1)
         if(sendInfo(i) > 0 .and. i /= myId) then
           commPattern%nProcSend         = commPattern%nProcSend + 1
           commPattern%indexSendProc(i) = commPattern%nProcSend

           sizeBuffer = max(sizeBuffer,sendInfo(i))
         endif
       enddo

       ! Allocate the memory for sendProc, nsend, nsendCum
       ! and sendList

       ii = commPattern%nProcSend
       allocate(commPattern%sendProc(ii),   &
                commPattern%nsend(ii),       &
                commPattern%nsendCum(0:ii), &
                commPattern%sendList(ii),   stat=ierr)
       if(ierr /= 0)                             &
         call terminate("finalCommStructures", &
                        "Memory allocation error for sendProc, etc")

       ! Repeat the loop over the number of processors, but now store
       ! the sending info. Use ii as a counter.

       ii = 0
       commPattern%nsendCum(0) = 0
       do i=0,(nProc-1)
         if(sendInfo(i) > 0 .and. i /= myId) then
           ii = ii + 1

           commPattern%sendProc(ii) = i
           commPattern%nsend(ii)     = sendInfo(i)
           commPattern%nsendCum(ii) = commPattern%nsendCum(ii-1) &
                                      + sendInfo(i)
         endif
       enddo

       ! Allocate the memory for the send lists.

       do i=1,commPattern%nProcSend
         ii = commPattern%nsend(i)
         allocate(commPattern%sendList(i)%block(ii),          &
                  commPattern%sendList(i)%indices(ii,3),      &
                  commPattern%sendList(i)%interp(ii,ninterp), &
                  stat=ierr)
         if(ierr /= 0)                             &
           call terminate("finalCommStructures", &
                          "Memory allocation failure for block, interp, &
                          &and indices of sendlist")
       enddo

       ! Allocate the memory for the receiving buffers.

       allocate(recvBuf(4,sizeBuffer), &
                recvBufInt(ninterp,sizeBuffer), stat=ierr)
       if(ierr /= 0)                             &
         call terminate("finalCommStructures", &
                        "Memory allocation failure for recvBuffers")

       ! Loop over the number of processors to receive my send info.

       recvSendInfo: do i=1,commPattern%nProcSend

         ! Block until a message arrives.

         call mpi_probe(mpi_any_source, myId+2, SUmb_comm_world, &
                       status, ierr)

         ! Store the source processor a bit easier and determine the
         ! index in the send data structure.

         source = status(mpi_source)
         ii     = commPattern%indexSendProc(source)

         ! Perform some checks in debug mode.

         if( debug ) then
           if(ii <= 0 .or. ii > commPattern%nProcSend) &
             call terminate("finalCommStructures",     &
                            "Send processor not in the list")

           call mpi_get_count(status, sumb_integer, sizeRecv, ierr)
           if(sizeRecv /= 4*commPattern%nsend(ii)) &
             call terminate("finalCommStructures",     &
                            "Unexpected size of message")
         endif

         ! Receive the message. As it has already arrived a blocking
         ! receive can be used.

         sizeRecv = 4*commPattern%nsend(ii)

         call mpi_recv(recvBuf, sizeRecv, sumb_integer, source, &
                       myId+2, SUmb_comm_world, status, ierr)

          ! Now receive the interpolants, if any.

         if(ninterp > 0) then
           sizeRecv = nInterp*commPattern%nsend(ii)
           call mpi_recv(recvBufInt, sizeRecv, sumb_real, source, &
                         myId+3, SUmb_comm_world, status, ierr)
         end if

         ! Store the info I must send to this processor in a normal
         ! data exchange.

         do j=1,commPattern%nsend(ii)
           commPattern%sendList(ii)%block(j) = recvBuf(1,j)

           commPattern%sendList(ii)%indices(j,1) = recvBuf(2,j)
           commPattern%sendList(ii)%indices(j,2) = recvBuf(3,j)
           commPattern%sendList(ii)%indices(j,3) = recvBuf(4,j)

           do pp = 1,ninterp
             commPattern%sendList(ii)%interp(j,pp) = recvBufInt(pp,j)
           end do
         enddo

       enddo recvSendInfo

       ! Complete the nonblocking sends. Dest is used as a temporary
       ! storage such that an integer is passed to the MPI call

       dest = commPattern%nProcRecv
       do i=1,commPattern%nProcRecv
         call mpi_waitany(dest, sendRequests, count, status, ierr)
       enddo

       ! Repeat the call if any interpolants were exchanged.

       if(ninterp > 0) then
         do i=1,commPattern%nProcRecv
           call mpi_waitany(dest, recvRequests, count, status, ierr)
         enddo
       end if

       ! Release the memory of the help variables allocated in this
       ! subroutine.

       deallocate(nHaloPerProc, sendInfo, buffer, recvBuf, &
                  bufInt, recvBufInt, stat=ierr)
       if(ierr /= 0) &
         call terminate("finalCommStructures", &
                        "Deallocation error for nHaloPerProc, &
                        &sendInfo, and temporary buffers")

       end subroutine finalCommStructures

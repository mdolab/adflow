!
!      ******************************************************************
!      *                                                                *
!      * File:          setBlockOversetData.f90                         *
!      * Author:        Steve Repsher                                   *
!      * Starting date: 08-23-2005                                      *
!      * Last modified: 10-14-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine setBlockOversetData(level, sps)
!
!      ******************************************************************
!      *                                                                *
!      * setBlockOversetData uses the overset communication lists for   *
!      * the given level and spectral mode setup by the donor search    *
!      * routines to fill in the overset arrays of the flow domains.    *
!      *                                                                *
!      ******************************************************************
!
       use block
       use communication
       use searchMod
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType), intent(in) :: level, sps
!
!      Local variables.
!
       integer :: size, procId, ierr, index
       integer, dimension(mpi_status_size) :: status

       integer(kind=intType) :: n, i, j, ii, jj
 
       integer(kind=intType), dimension(nDom) :: cc

       integer, dimension(:), allocatable :: recvRequests2, sendRequests2

       integer(kind=intType), dimension(:,:), allocatable :: sendBuf
       integer(kind=intType), dimension(:,:), allocatable :: recvBuf

       real(kind=realType), dimension(:,:), allocatable :: sendBufInt
       real(kind=realType), dimension(:,:), allocatable :: recvBufInt
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Allocate memory for the sending and receiving buffers, and a 
       ! 2nd set of message requests.

       i  = commPatternOverset(level,sps)%nProcSend
       ii = commPatternOverset(level,sps)%nSendCum(i)
       j  = commPatternOverset(level,sps)%nProcRecv
       jj = commPatternOverset(level,sps)%nRecvCum(j)

       allocate(sendBuf(4,ii),          recvBuf(4,jj),          &
                sendBufInt(nInterp,ii), recvBufInt(nInterp,jj), &
                sendRequests2(i),       recvRequests2(j), stat=ierr)
       if(ierr /= 0)                           &
         call terminate("setBlockOversetData", &
                        "Memory allocation failure for buffers")
 
       ! Loop over the processors i need to send to and their lists.

       ii = 1
       jj = 1
       sends: do i=1,commPatternOverset(level,sps)%nProcSend

         do j=1,commPatternOverset(level,sps)%nSend(i)

           ! Copy the block, indices, and interpolants to the buffers.

           sendBuf(1,jj) = &
                  commPatternOverset(level,sps)%sendList(i)%block(j)
           sendBuf(2:4,jj) = &
                  commPatternOverset(level,sps)%sendList(i)%indices(j,:)
           sendBufInt(:,jj) = &
                  commPatternOverset(level,sps)%sendList(i)%interp(j,:)

           jj = jj + 1

         end do

         ! Store the processor id and size then send the messages.

         procID = commPatternOverset(level,sps)%sendProc(i)

         size = 4*commPatternOverset(level,sps)%nSend(i)
         call mpi_isend(sendBuf(1,ii), size, sumb_integer, procId, &
                        1, SUmb_comm_world, sendRequests(i), ierr)

         size = nInterp*commPatternOverset(level,sps)%nSend(i)
         call mpi_isend(sendBufInt(1,ii), size, sumb_real, procId, &
                        2, SUmb_comm_world, sendRequests2(i), ierr)

         ! Update ii for the next processor.

         ii = jj

       end do sends

       ! Post the nonblocking receives.

       ii = 1
       receives: do i=1,commPatternOverset(level,sps)%nProcRecv

         ! Store the processor id and the size of the message, then
         ! post them.

         procID = commPatternOverset(level,sps)%recvProc(i)

         size = 4*commPatternOverset(level,sps)%nRecv(i)
         call mpi_irecv(recvBuf(1,ii), size, sumb_integer, procId, &
                        1, SUmb_comm_world, recvRequests(i), ierr)

         size = nInterp*commPatternOverset(level,sps)%nRecv(i)
         call mpi_irecv(recvBufInt(1,ii), size, sumb_real, procId, &
                        2, SUmb_comm_world, recvRequests2(i), ierr)

         ii = ii + commPatternOverset(level,sps)%nRecv(i)

       end do receives

       ! Initialize the overset cell counter for the domains to 0.

       cc = 0

       ! Allocate the overset data arrays for each block. Note that
       ! ibndry includes space for storing the orphan indices.

       do n = 1,nDom

         i = flowDoms(n,level,sps)%nCellsOverset
         j = i + flowDoms(n,level,sps)%nOrphans

         allocate(flowDoms(n,level,sps)%ibndry(3,j),        &
                  flowDoms(n,level,sps)%idonor(3,i),        &
                  flowDoms(n,level,sps)%neighProcOver(i),   &
                  flowDoms(n,level,sps)%neighBlockOver(i),  &
                  flowDoms(n,level,sps)%overint(nInterp,i), &
                  stat=ierr)
         if(ierr /= 0)                           &
           call terminate("setBlockOversetData", &
                          "Memory allocation failure for block data")
 
       end do

       ! Loop over the local communication.

       localComm: do i=1,internalOverset(level,sps)%nCopy

         ! Copy the data to the appropriate block.

         j     = internalOverset(level,sps)%haloBlock(i)
         cc(j) = cc(j) + 1

         flowDoms(j,level,sps)%ibndry(:,cc(j))       = &
                            internalOverset(level,sps)%haloIndices(i,:)
         flowDoms(j,level,sps)%neighProcOver(cc(j))  = myID
         flowDoms(j,level,sps)%neighBlockOver(cc(j)) = &
                            internalOverset(level,sps)%donorBlock(i)
         flowDoms(j,level,sps)%idonor(:,cc(j))       = &
                            internalOverset(level,sps)%donorIndices(i,:)
         flowDoms(j,level,sps)%overint(:,cc(j))      = &
                            internalOverset(level,sps)%donorInterp(i,:)

       enddo localComm

       ! Complete the nonblocking receives in an arbitrary sequence.

       size = commPatternOverset(level,sps)%nProcRecv
       completeRecvs: do i=1,commPatternOverset(level,sps)%nProcRecv

         ! Complete any of the requests (not the interpolants).

         call mpi_waitany(size, recvRequests, index, status, ierr)

         ! Find the starting point in the buffer for this message and
         ! wait to make sure the interpolants have arrived.

         ii = index
         jj = commPatternOverset(level,sps)%nRecvCum(ii-1)

         call mpi_waitany(1, recvRequests2(ii), index, status, ierr)

         ! Loop over the message.

         do j = 1,commPatternOverset(level,sps)%nRecv(ii)

           ! Copy the data to the appropriate block.

           jj    = jj + 1
           n     = commPatternOverset(level,sps)%recvList(ii)%block(j)
           cc(n) = cc(n) + 1

           flowDoms(n,level,sps)%ibndry(:,cc(n))       = &
                 commPatternOverset(level,sps)%recvList(ii)%indices(j,:)
           flowDoms(n,level,sps)%neighProcOver(cc(n))  = &
                              commPatternOverset(level,sps)%recvProc(ii)

           flowDoms(n,level,sps)%neighBlockOver(cc(n)) = recvBuf(1  ,jj)
           flowDoms(n,level,sps)%idonor(:,cc(n))       = recvBuf(2:4,jj)
           flowDoms(n,level,sps)%overint(:,cc(n))      = recvBufInt(:,jj)

         end do
       enddo completeRecvs

       ! Complete the nonblocking sends.

       size = commPatternOverset(level,sps)%nProcSend
       do i=1,commPatternOverset(level,sps)%nProcSend
         call mpi_waitany(size, sendRequests,  index, status, ierr)
         call mpi_waitany(size, sendRequests2, index, status, ierr)
       enddo

       ! Deallocate the memory for the sending and receiving buffers
       ! and extra set of requests.
 
       deallocate(sendBuf, recvBuf, sendBufInt, recvBufInt, &
                  sendRequests2, recvRequests2, stat=ierr)
       if(ierr /= 0)                        &
         call terminate("setBlockOversetData", &
                        "Deallocation failure for buffers")

       end subroutine setBlockOversetData

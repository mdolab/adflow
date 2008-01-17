!
!      ******************************************************************
!      *                                                                *
!      * File:          coarseDonorInfo.f90                             *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 09-05-2003                                      *
!      * Last modified: 06-12-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine coarseDonorInfo(level)
!
!      ******************************************************************
!      *                                                                *
!      * coarseDonorInfo creates the donor info for the internal        *
!      * block boundaries on the given coarse grid level.               *
!      *                                                                *
!      ******************************************************************
!
       use communication
       use coarse1to1Subface
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType), intent(in) :: level
!
!      Local variables.
!
       integer :: proc, size, ierr

       integer, dimension(mpi_status_size) :: status
       integer, dimension(nProc)           :: sizeMessage

       integer(kind=intType) :: nn, mm, ii, i
       integer(kind=intType) :: nMessageReceive, nMessageSend

       integer(kind=intType), dimension(0:nProc) :: size2proc
       integer(kind=intType), dimension(nProc)   :: tmp

       integer(kind=intType), dimension(:), allocatable :: sendBuf
       integer(kind=intType), dimension(:), allocatable, target :: recvBuf

       type(coarse1to1SubfaceType) :: tmpSubface
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Determine the number of integers sent to every processor.

       size2proc = 0
       do nn=1,nSubface1to1
         mm = subface1to1(nn)%neighProc + 1    ! Proc ID's start at 0.

         size2proc(mm) = size2proc(mm) + 13                          &
                       + subface1to1(nn)%iEnd - subface1to1(nn)%iBeg &
                       + subface1to1(nn)%jEnd - subface1to1(nn)%jBeg &
                       + subface1to1(nn)%kEnd - subface1to1(nn)%kBeg
       enddo

       ! No message is sent to myself, so set the corresponding entry in
       ! size2proc to 0. The processor id's start at 0, which explains
       ! the +1.

       size2proc(myID+1) = 0

       ! Determine the number of messages i have to sent. Store in tmp
       ! a 0 or a 1, depending whether or not a message must be sent
       ! to the corresponding processor. This info is needed to determine
       ! the number of messages i will receive.

       nMessageSend = 0
       do nn=1,nProc
         if(size2proc(nn) > 0) then
           nMessageSend = nMessageSend +1
           tmp(nn) = 1
         else
           tmp(nn) = 0
         endif
       enddo

       ! Determine the number of messages i will receive.

       sizeMessage = 1
       call mpi_reduce_scatter(tmp, nMessageReceive, sizeMessage,      &
                               sumb_integer, mpi_sum, SUmb_comm_world, &
                               ierr)

       ! Put size2proc in cumulative storage format and store the
       ! starting entry in tmp, which is used as a counter.

       do nn=1,nProc
         size2proc(nn) = size2proc(nn) + size2proc(nn-1)
         tmp(nn)       = size2proc(nn-1)
       enddo

       ! Allocate the memory for the send buffer.

       allocate(sendBuf(size2proc(nProc)), stat=ierr)
       if(ierr /= 0)                       &
         call terminate("coarseDonorInfo", &
                        "Memory allocation failure for sendBuf")

       ! Loop over the number of 1 to 1 subfaces to fill the send buffer.

       unOwnedSubfaces: do nn=1,nSubface1to1

         ! Only info that must be sent to other processors must be stored.

         if(subface1to1(nn)%neighProc /= myID) then

           ! Store the entry in tmp in mm.
           ! Note that the proc id's start at 0.

           mm = subface1to1(nn)%neighProc + 1

           ! Store the block id of the donor and the coarse grid range of
           ! the current subface in the send buffer. Update the counter
           ! tmp(mm) accordingly.

           tmp(mm) = tmp(mm) +1
           sendBuf(tmp(mm)) = subface1to1(nn)%neighBlock

           tmp(mm) = tmp(mm) +1
           sendBuf(tmp(mm)) = subface1to1(nn)%iBeg

           tmp(mm) = tmp(mm) +1
           sendBuf(tmp(mm)) = subface1to1(nn)%jBeg

           tmp(mm) = tmp(mm) +1
           sendBuf(tmp(mm)) = subface1to1(nn)%kBeg

           tmp(mm) = tmp(mm) +1
           sendBuf(tmp(mm)) = subface1to1(nn)%iEnd

           tmp(mm) = tmp(mm) +1
           sendBuf(tmp(mm)) = subface1to1(nn)%jEnd

           tmp(mm) = tmp(mm) +1
           sendBuf(tmp(mm)) = subface1to1(nn)%kEnd

           ! Store the number of points in the three coordinate
           ! directions of the coarse grid donor face.

           tmp(mm) = tmp(mm) +1
           sendBuf(tmp(mm)) = subface1to1(nn)%ndi

           tmp(mm) = tmp(mm) +1
           sendBuf(tmp(mm)) = subface1to1(nn)%ndj

           tmp(mm) = tmp(mm) +1
           sendBuf(tmp(mm)) = subface1to1(nn)%ndk

           ! Store the i-donor indices of the fine grid in sendBuf.

           do i=1,subface1to1(nn)%ndi
             tmp(mm) = tmp(mm) +1
             sendBuf(tmp(mm)) = subface1to1(nn)%idfine(i)
           enddo

           ! Store the j-donor indices of the fine grid in sendBuf.

           do i=1,subface1to1(nn)%ndj
             tmp(mm) = tmp(mm) +1
             sendBuf(tmp(mm)) = subface1to1(nn)%jdfine(i)
           enddo

           ! Store the k-donor indices of the fine grid in sendBuf.

           do i=1,subface1to1(nn)%ndk
             tmp(mm) = tmp(mm) +1
             sendBuf(tmp(mm)) = subface1to1(nn)%kdfine(i)
           enddo

         endif

       enddo unOwnedSubfaces

       ! Send the data i have to send.

       mm = 0
       sends: do nn=1,nProc

         if(size2proc(nn) > size2proc(nn-1)) then

           ! Update mm and store the processor id, the size of the
           ! message and the starting index in sendbuf.

           mm   = mm +1
           proc = nn -1
           size = size2proc(nn) - size2proc(nn-1)
           ii   = size2proc(nn-1) +1

           ! Send the message.

           call mpi_isend(sendBuf(ii), size, sumb_integer, proc, proc, &
                          SUmb_comm_world, sendRequests(mm), ierr)

         endif

       enddo sends

       ! Loop over the number of 1 to 1 subfaces stored on this processor
       ! and update the locally stored info.

       do nn=1,nSubface1to1
         if(subface1to1(nn)%neighProc == myID) &
           call update1to1Coarse(level, subface1to1(nn))
       enddo

       ! Release the memory of subface1to1.

       do nn=1,nSubface1to1
         deallocate(subface1to1(nn)%idfine, subface1to1(nn)%jdfine, &
                    subface1to1(nn)%kdfine, stat=ierr)
         if(ierr /= 0)                       &
           call terminate("coarseDonorInfo", &
                          "Deallocation error for idfine, etc")
       enddo

       deallocate(subface1to1, stat=ierr)
       if(ierr /= 0)                       &
         call terminate("coarseDonorInfo", &
                        "Deallocation error for subface1to1")

       ! Loop over the number of messages i must receive to determine
       ! info from externally stored neighboring subfaces.

       receives: do nn=1,nMessageReceive

         ! Block until a message arrives.

         call mpi_probe(mpi_any_source, myID, SUmb_comm_world, &
                       status, ierr)

         ! Find the source and size of the message.

         proc = status(mpi_source)
         call mpi_get_count(status, sumb_integer, size, ierr)

         ! Check in debug mode that the incoming message is of
         ! correct size.

         if( debug ) then
           if(size == mpi_undefined)             &
             call terminate("coarseDonorInfo", &
                            "Unexpected size of message")
         endif

         ! Allocate the memory for the receive buffer.

         allocate(recvBuf(size), stat=ierr)
         if(ierr /= 0)                         &
           call terminate("coarseDonorInfo", &
                          "Memory allocation failure for recvBuf")

         ! Receive the messsage. As it has already arrived a blocking
         ! receive can be used.

         call mpi_recv(recvBuf, size, sumb_integer, proc, myID, &
                       SUmb_comm_world, status, ierr)

         ! Loop to extract the 1 to 1 subface info from the buffer.

         ii = 1
         extractSubface: do
           ! Exit the loop if no more info is present.

           if(ii > size) exit

           ! Store the info of this subface in tmpSubface.

           tmpSubface%neighProc  = proc
           tmpSubface%neighBlock = recvBuf(ii); ii = ii +1

           tmpSubface%iBeg = recvBuf(ii); ii = ii +1
           tmpSubface%jBeg = recvBuf(ii); ii = ii +1
           tmpSubface%kBeg = recvBuf(ii); ii = ii +1

           tmpSubface%iEnd = recvBuf(ii); ii = ii +1
           tmpSubface%jEnd = recvBuf(ii); ii = ii +1
           tmpSubface%kEnd = recvBuf(ii); ii = ii +1

           tmpSubface%ndi = recvBuf(ii); ii = ii +1
           tmpSubface%ndj = recvBuf(ii); ii = ii +1
           tmpSubface%ndk = recvBuf(ii); ii = ii +1

           mm = ii + tmpSubface%ndi
           tmpSubface%idfine => recvBuf(ii:(mm-1))
           ii = mm

           mm = ii + tmpSubface%ndj
           tmpSubface%jdfine => recvBuf(ii:(mm-1))
           ii = mm

           mm = ii + tmpSubface%ndk
           tmpSubface%kdfine => recvBuf(ii:(mm-1))
           ii = mm

           ! Update the 1 to 1 subface info.

           call update1to1Coarse(level, tmpSubface)

         enddo extractSubface

         ! Release the memory of the receive buffer.

         deallocate(recvBuf, stat=ierr)
         if(ierr /= 0)                       &
           call terminate("coarseDonorInfo", &
                          "Deallocation failure for recvBuf")

       enddo receives

       ! Complete the nonblocking sends.

       size = nMessageSend
       do nn=1,nMessageSend
         call mpi_waitany(size, sendRequests, proc, status, ierr)
       enddo

       ! Release the memory of the send buffer.

       deallocate(sendBuf, stat=ierr)
       if(ierr /= 0)                       &
         call terminate("coarseDonorInfo", &
                        "Deallocation failure for sendBuf")

       ! Synchronize the processors, because wild cards have been used.

       call mpi_barrier(SUmb_comm_world, ierr)

       end subroutine coarseDonorInfo

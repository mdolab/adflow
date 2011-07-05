!
!      ******************************************************************
!      *                                                                *
!      * File:          indirectHalosPerLevel.f90                       *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 02-03-2003                                      *
!      * Last modified: 11-30-2007                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine indirectHalosPerLevel(level, iihalo, entityHalo, &
                                        transform, entityIndex)
!
!      ******************************************************************
!      *                                                                *
!      * indirectHalosPerLevel determines the donor cells for the       *
!      * halo's of the given level of indirectness. From the known      *
!      * appropriate direct halo and its donor, the corresponding cell  *
!      * in the donor block is determined.                              *
!      *                                                                *
!      ******************************************************************
!
       use haloList
       use indirectHalo
       use communication
       use BCTypes
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType), intent(in)    :: level
       integer(kind=intType), intent(inout) :: iihalo

       integer(kind=intType), dimension(:,:), intent(in) :: transform

       type(haloListType),  dimension(:), intent(inout) :: entityHalo
       type(indexListType), dimension(:), intent(inout) :: entityIndex
!
!      Local variables.
!
       integer :: ierr
       integer :: count, dest, source, sizeRecv

       integer, dimension(mpi_status_size)  :: status

       integer, allocatable, dimension(:) :: sizeMessage

       integer(kind=intType) :: i, ii, nn, mm, ms
       integer(kind=intType) :: l1, l2, l3, proc
       integer(kind=intType) :: start, end, nProcsSend, nProcsRecv
       integer(kind=intType) :: sizeSendBuf, sizeRecvBuf
       integer(kind=intType) :: nLocalHalos
       integer(kind=intType) :: nItemReturn, nItemSend, nItemAlloc

       integer(kind=intType), dimension(2) :: tmpBuf

       integer(kind=intType), allocatable, dimension(:) :: counter
       integer(kind=intType), allocatable, dimension(:) :: sendBuf
       integer(kind=intType), allocatable, dimension(:) :: recvBuf
       integer(kind=intType), allocatable, dimension(:) :: localHalos
!
!      Interfaces
!
       interface
       subroutine fillSendBuf(sendBuf, proc, entityHalo, transform, &
                              level, mm)
       use haloList
       use indirectHalo
       implicit none
       integer(kind=intType), intent(in)  :: proc, level
       integer(kind=intType), intent(out) :: mm

       integer(kind=intType), dimension(:), intent(out) :: sendBuf

       integer(kind=intType), dimension(:,:), intent(in) :: transform

       type(haloListType), dimension(:), intent(in) :: entityHalo
       end subroutine fillSendBuf

       !=================================================================

       subroutine findDonorsRecvBuffer(recvBuf, nHalos, entityHalo, &
                                       entityIndex, level, nItemReturn)
       use haloList
       use communication
       implicit none
       integer(kind=intType), intent(in) :: nHalos, level, nItemReturn

       integer(kind=intType), dimension(:), intent(inout) :: recvBuf
       type(haloListType),  dimension(:), intent(in) :: entityHalo
       type(indexListType), dimension(:), intent(in) :: entityIndex
       end subroutine findDonorsRecvBuffer

       !=================================================================

       subroutine storeHalosInList(buffer, bufSize, proc, level, &
                                   nItemReturn, entityHalo,      &
                                   entityIndex, iihalo)
       use haloList
       use indirectHalo
       implicit none
       integer(kind=intType), intent(in)    :: bufSize, proc
       integer(kind=intType), intent(in)    :: level, nItemReturn
       integer(kind=intType), intent(inout) :: iihalo

       integer(kind=intType), dimension(:), intent(in) :: buffer

       type(haloListType),  dimension(:), intent(inout) :: entityHalo
       type(indexListType), dimension(:), intent(inout) :: entityIndex
       end subroutine storeHalosInList
       end interface
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Determine the number of items per entity that is returned from
       ! the requested processor. This is 5 + level, because of the
       ! storage of possible periodic transformations.
       ! Also set nItemSend to 7 and determine the number of items that
       ! should be allocated in the communication buffers.

       nItemReturn = 5 + level
       nItemSend   = 7
       nItemAlloc  = max(nItemReturn, nItemSend)

       ! Abbreviate the start and ending index for this level
       ! a bit easier.

       start = nHaloPerLev(level-1) +1
       end   = nHaloPerLev(level)

       ! Allocate the memory for sizeMessage and counter, which will
       ! be used in reduceScatter to determine the amount of messages
       ! and the total size.

       allocate(sizeMessage(nProc), counter(2*nProc), stat=ierr)
       if(ierr /= 0)                             &
         call terminate("indirectHalosPerLevel", &
                        "Memory allocation failure for sizeMessage &
                        &and counter")

       ! Determine the amount of halo's per processors. Here use is made
       ! of the fact that boundary halo's are flagged with a processor
       ! ID of -1.

       nHaloPerProc = 0
       do i=start,end
         ii               = indHalo(i)%donorProc+1
         nHaloPerProc(ii) = nHaloPerProc(ii) +1
       enddo

       ! Determine the amount of messages as well as the total number of
       ! halos I have to send and put nHaloPerProc in cumulative
       ! storage format.

       nProcsSend  = 0
       sizeSendBuf = 0

       do i=1,nProc
         if(nHaloPerProc(i) > 0 .and. i /= (myID+1)) then

           ! Something must be send to this proc. Update nProcsSend
           ! and sizeSendBuf and set the two elements of counter
           ! appropriately.

           nProcsSend     = nProcsSend +1
           sizeSendBuf    = sizeSendBuf + nHaloPerProc(i)
           counter(2*i-1) = 1
           counter(2*i)   = nHaloPerProc(i)

         else

           ! Nothing needs to be sent. Set the two elements of
           ! counter to 0.

           counter(2*i-1) = 0
           counter(2*i)   = 0

         endif

         ! Set the size of the message to 2 and put nHaloPerProc in
         ! cumulative storage format.

         sizeMessage(i)  = 2
         nHaloPerProc(i) = nHaloPerProc(i) + nHaloPerProc(i-1)

       enddo

       ! Call reduceScatter to determine the number of processors
       ! from which I receive data and the total amount of data.
       ! tmpBuf is used as a temporary buffer to receive the data.

       call mpi_reduce_scatter(counter, tmpBuf, sizeMessage,           &
                               sumb_integer, mpi_sum, SUmb_comm_world, &
                               ierr)

       nProcsRecv  = tmpBuf(1)
       sizeRecvBuf = tmpBuf(2)

       ! Allocate the memory for the send buffer.

       allocate(sendBuf(nItemAlloc*sizeSendBuf), stat=ierr)
       if(ierr /= 0)                                           &
         call terminate("indirectHalosPerLevel",               &
                        "Memory allocation error for sendBuf.")

       ! Send the data I must send. Use nonblocking sends to
       ! avoid deadlock. nn is used as a counter for the current
       ! message to be sent and ms as starting index for the
       ! active message in sendBuf.

       nn = 0
       ms = 1
       sendproc: do i=1,nProc
         if(counter(2*i) > 0) then

           ! Something must be send to this processor. Update nn.

           nn = nn +1
 
           ! Fill this part of the send buffer

           proc = i-1
           call fillSendBuf(sendBuf(ms:), proc, entityHalo, &
                            transform, level, mm)

           ! Send this buffer. Make sure that integers are used for
           ! count, destination and tag.

           count = nItemSend*mm
           dest  = i-1
           call mpi_isend(sendBuf(ms), count, sumb_integer, dest, dest, &
                          SUmb_comm_world, sendRequests(nn), ierr)

           ! Update ms to the starting index in buffer for the next
           ! message to be sent.

           ms = ms + nItemAlloc*mm

         endif
       enddo sendproc

       ! Loop over the boundary halos for this level. The value of start
       ! defined earlier is still okay; only end needs to be redefined.

       end = nHaloPerLev(level-1) + nHaloPerProc(0)

       bocos: do i=start,end

         ! Determine the corresponding direct boundary halo.

         ii = indHalo(i)%myDirectHalo

         ! Determine the vector from the direct halo to its donor.
         ! Store this in l1, l2 and l3.

         l1 = entityHalo(ii)%dI - entityHalo(ii)%myI
         l2 = entityHalo(ii)%dJ - entityHalo(ii)%myJ
         l3 = entityHalo(ii)%dK - entityHalo(ii)%myK

         ! Store the info for this boundary halo. Set the boundary
         ! condition, i.e. donorBlock, to the boundary condition of the
         ! closest direct halo. This has been chosen to be the most
         ! important boundary condition.
         ! DonorProc is set to -1 to indicate a boundary halo.

         iihalo = iihalo +1
         entityHalo(iihalo)%myBlock = indHalo(i)%myBlock
         entityHalo(iihalo)%myI     = indHalo(i)%myI
         entityHalo(iihalo)%myJ     = indHalo(i)%myJ
         entityHalo(iihalo)%myK     = indHalo(i)%myK

         entityHalo(iihalo)%donorProc  = -1
         entityHalo(iihalo)%donorBlock = entityHalo(ii)%donorBlock

         ! Determine the vector from my indices to the indices of my
         ! donor. This is the vector from me to the donor indices of
         ! the closest halo, projected on the vector (l1,l2,l3).
         ! Although in general this will not give to indices of my true
         ! donor, this info is stored, because the halo list could be
         ! both a node or a cell halo list. The true donor will be
         ! extracted later.
         ! This implies that my donor could also be a halo, but its
         ! level of indirectness is guaranteed one less than mine.

         ii = l1*(entityHalo(ii)%dI - entityHalo(iihalo)%myI) &
            + l2*(entityHalo(ii)%dJ - entityHalo(iihalo)%myJ) &
            + l3*(entityHalo(ii)%dK - entityHalo(iihalo)%myK)

         entityHalo(iihalo)%dI = ii*l1 + indHalo(i)%myI
         entityHalo(iihalo)%dJ = ii*l2 + indHalo(i)%myJ
         entityHalo(iihalo)%dK = ii*l3 + indHalo(i)%myK

         ! Copy the level of indirectness.

         entityHalo(iihalo)%levOfInd = indHalo(i)%levOfInd

         ! Store the entry of entityHalo in the i,j,k indices
         ! of in entityIndex.

         ii = indHalo(i)%myBlock
         l1 = indHalo(i)%myI
         l2 = indHalo(i)%myJ
         l3 = indHalo(i)%myK

         entityIndex(ii)%entryList(l1,l2,l3) = iihalo

       enddo bocos

       ! Treat the internal block boundary halo's, whose donor is stored
       ! on the same processor. Store this data in localHalos, for
       ! which the memory must be allocated. First determine the number
       ! of local halo's.

       nLocalHalos = nHaloPerProc(myID+1) - nHaloPerProc(myID)
       allocate(localHalos(nItemAlloc*nLocalHalos), stat=ierr)
       if(ierr /= 0)                             &
         call terminate("indirectHalosPerLevel", &
                        "Memory allocation error for localHalos")

       proc = myID
       call fillSendBuf(localHalos, proc, entityHalo, transform, &
                        level, mm)

       ! Receive the messages in arbitrary sequence, determine the
       ! corresponding halo info and send the message back.

       allocate(recvBuf(nItemAlloc*sizeRecvBuf),  stat=ierr)
       if(ierr /= 0)                             &
         call terminate("indirectHalosPerLevel", &
                        "Memory allocation error for recvBuf.")

       ! Initialize ms to 1 and start the loop over the number of
       ! messages to be received.

       ms = 1
       recvproc: do i=1,nProcsRecv

         ! Block until a message arrives.

         call mpi_probe(mpi_any_source, myID, SUmb_comm_world, &
                       status, ierr)

         ! Find the source and size of the message.

         source = status(mpi_source)
         call mpi_get_count(status, sumb_integer, sizeRecv, ierr)

         ! Check in debug mode that the incoming message is of
         ! correct size.

         if( debug ) then
           if(sizeRecv == mpi_undefined .or. &
              mod(sizeRecv,nItemSend) /= 0)  &
             call terminate("indirectHalosPerLevel", &
                            "Unexpected size of message")
         endif

         ! Receive the message. As it has already arrived a blocking
         ! receive can be used.

         call mpi_recv(recvBuf(ms), sizeRecv, sumb_integer, &
                       source, myID, SUmb_comm_world, status, ierr)

         ! Determine the number of halo's in the receive buffer.

         mm = sizeRecv/nItemSend

         ! Find the donors for the halo's in the receive buffer as
         ! well as the periodic info.

         call findDonorsRecvBuffer(recvBuf(ms:), mm, entityHalo, &
                                   entityIndex, level, nItemReturn)

         ! Send the modified receive buffer back to the source processor.

         count = nItemReturn*mm
         call mpi_isend(recvBuf(ms), count, sumb_integer, source, &
                        source+1, SUmb_comm_world, recvRequests(i), ierr)

         ! Update the starting index ms for the next message.

         ms = ms + nItemAlloc*mm

       enddo recvproc

       ! Find the donors for the locally stored halo's and store them
       ! in the list.

       call findDonorsRecvBuffer(localHalos, nLocalHalos, entityHalo, &
                                 entityIndex, level, nItemReturn)

       proc = myID
       call storeHalosInList(localHalos, nLocalHalos, proc, level, &
                             nItemReturn, entityHalo, entityIndex, &
                             iihalo)

       ! Complete the 1st series of nonblocking sends.

       do i=1,nProcsSend
         call mpi_waitany(nProcsSend, sendRequests, count, status, ierr)
       enddo

       ! Loop over the processors to which I sent data to find out the
       ! halo information. Now these messages must be received.

       secondRecv: do i=1,nProcsSend

         ! Block until a message arrives.

         call mpi_probe(mpi_any_source, myID+1, SUmb_comm_world, &
                       status, ierr)

         ! Find the source and size of the message.

         source = status(mpi_source)
         call mpi_get_count(status, sumb_integer, sizeRecv, ierr)

         ! Check in debug mode that the incoming message is of
         ! correct size.

         if( debug ) then
           if(sizeRecv == mpi_undefined .or.         &
              mod(sizeRecv,nItemReturn) /= 0)        &
             call terminate("indirectHalosPerLevel", &
                            "Unexpected size of message")
         endif

         ! Store the number of halo's in mm.

         mm = sizeRecv/nItemReturn

         ! Receive the message. Use a blocking receive, as the message
         ! has already arrived.

         call mpi_recv(sendBuf, sizeRecv, sumb_integer, source, &
                       myID+1, SUmb_comm_world, status, ierr)

         ! Store the donors in the list.

         proc = source
         call storeHalosInList(sendBuf, mm, proc, level, nItemReturn, &
                               entityHalo, entityIndex, iihalo)
       enddo secondRecv

       ! Complete the second series of nonblocking sends.

       do i=1,nProcsRecv
         call mpi_waitany(nProcsRecv, recvRequests, count, status, ierr)
       enddo

       ! Release the memory allocated in this subroutine.

       deallocate(sizeMessage, counter, sendBuf, localHalos, &
                  recvBuf, stat=ierr)
       if(ierr /= 0)                             &
         call terminate("indirectHalosPerLevel", &
                        "Deallocation error for sizeMessage, etc.")

       end subroutine indirectHalosPerLevel

       !=================================================================

       subroutine fillSendBuf(sendBuf, proc, entityHalo, transform, &
                              level, mm)
!
!      ******************************************************************
!      *                                                                *
!      * fillSendBuf fills the buffer, which must be sent to the        *
!      * given processor.                                               *
!      *                                                                *
!      ******************************************************************
!
       use haloList
       use indirectHalo
       implicit none
!
!      Subroutine arguments
!
       integer(kind=intType), intent(in)  :: proc, level
       integer(kind=intType), intent(out) :: mm

       integer(kind=intType), dimension(:), intent(out) :: sendBuf
       integer(kind=intType), dimension(:,:), intent(in) :: transform

       type(haloListType), dimension(:), intent(in) :: entityHalo
!
!      Local variables.
!
       integer(kind=intType) :: i, ii, jj, m
       integer(kind=intType) :: l1, l2, l3

       integer(kind=intType), dimension(3,3) :: tMat
!
!      Function definition.
!
       integer(kind=intType) :: delta
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Initialize m to 0.

       m = 0

       ! Loop over the number of halo's for this processor.

       do i=nHaloPerProc(proc)+1,nHaloPerProc(proc+1)

         ! Abbreviate the current entry in indHalo in ii and the
         ! entry of its corresponding direct halo in jj.

         ii = i + nHaloPerLev(level-1)
         jj = indHalo(ii)%myDirectHalo

         ! Determine the transformation matrix between the block of
         ! the direct halo and its donor block.

         l1 = transform(jj,1)
         l2 = transform(jj,2)
         l3 = transform(jj,3)

         tMat(1,1) = sign(1_intType,l1) * delta(l1,1_intType)
         tMat(2,1) = sign(1_intType,l1) * delta(l1,2_intType)
         tMat(3,1) = sign(1_intType,l1) * delta(l1,3_intType)

         tMat(1,2) = sign(1_intType,l2) * delta(l2,1_intType)
         tMat(2,2) = sign(1_intType,l2) * delta(l2,2_intType)
         tMat(3,2) = sign(1_intType,l2) * delta(l2,3_intType)

         tMat(1,3) = sign(1_intType,l3) * delta(l3,1_intType)
         tMat(2,3) = sign(1_intType,l3) * delta(l3,2_intType)
         tMat(3,3) = sign(1_intType,l3) * delta(l3,3_intType)

         ! Store the direction from the direct to the indirect
         ! halo in l1, l2 and l3

         l1 = indHalo(ii)%myI - entityHalo(jj)%myI
         l2 = indHalo(ii)%myJ - entityHalo(jj)%myJ
         l3 = indHalo(ii)%myK - entityHalo(jj)%myK

         ! Fill the send buffer with block ID and i,j and k
         ! indices of the donor of the direct halo and the
         ! transformed direction to reach the donor of the
         ! indirect halo. For these last three values the
         ! transformation matrix must be applied to l1, l2, l3.

         m = m+1; sendBuf(m) = entityHalo(jj)%donorBlock
         m = m+1; sendBuf(m) = entityHalo(jj)%dI
         m = m+1; sendBuf(m) = entityHalo(jj)%dJ
         m = m+1; sendBuf(m) = entityHalo(jj)%dK

         m = m+1; sendBuf(m) = tMat(1,1)*l1 + tMat(1,2)*l2 + tMat(1,3)*l3
         m = m+1; sendBuf(m) = tMat(2,1)*l1 + tMat(2,2)*l2 + tMat(2,3)*l3
         m = m+1; sendBuf(m) = tMat(3,1)*l1 + tMat(3,2)*l2 + tMat(3,3)*l3

       enddo

       ! Set the return variable mm to the number of halos stored in
       ! the send buffer.

       mm = m/7

       end subroutine fillSendBuf

       !=================================================================

       subroutine findDonorsRecvBuffer(recvBuf, nHalos, entityHalo, &
                                       entityIndex, level, nItemReturn)
!
!      ******************************************************************
!      *                                                                *
!      * findDonorsRecvBuffer finds the donor cells for the halo        *
!      * information stored in recvBuf. On return recvBuf contains      *
!      * for every halo the following information: processor ID,        *
!      * block ID, the i,j,k indices of the donor cell and periodic     *
!      * information. The number of periodic subfaces stored is level,  *
!      * where a 0 indicates that the subface is not periodic.          *
!      *                                                                *
!      ******************************************************************
!
       use haloList
       use communication
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType), intent(in) :: nHalos, level, nItemReturn

       integer(kind=intType), dimension(:), intent(inout) :: recvBuf

       type(haloListType),  dimension(:), intent(in) :: entityHalo
       type(indexListType), dimension(:), intent(in) :: entityIndex
!
!      Local variables.
!
       integer :: ierr

       integer(kind=intType) :: i, j, k, ii, jj, kk, mm, nn
       integer(kind=intType) :: db, l1, L2, l3
       integer(kind=intType) :: nPeriodic

       integer(kind=intType), dimension(:), allocatable :: tmpBuf
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Allocate the memory for tmpBuf to the size needed to store the
       ! return information.

       allocate(tmpBuf(nItemReturn*nHalos), stat=ierr)
       if(ierr /= 0)                            &
         call terminate("findDonorsRecvBuffer", &
                        "Memory allocation failure for tmpBuf");

       ! Initialize nn and mm to 0. nn is the counter for the incoming
       ! halo information (7 per halo) and mm for the outgoing info
       ! (nItemReturn per halo).

       nn = 0
       mm = 0

       ! Loop over the number of halos in the receive buffer.

       do ii=1,nHalos

         ! Store the incoming information a bit easier.

         db = recvBuf(nn+1)
         i  = recvBuf(nn+2)
         j  = recvBuf(nn+3)
         k  = recvBuf(nn+4)
         l1 = recvBuf(nn+5)
         l2 = recvBuf(nn+6)
         l3 = recvBuf(nn+7)

         ! At the moment i,j,k are the indices of the donor of the direct
         ! halo and l1,l2,l3 the path from the direct halo to the
         ! indirect halo. Add l1, L2, l3 to i,j,k such that they store
         ! the indices of the donor of the indirect halo.

         i = i + l1
         j = j + l2
         k = k + l3

         ! Store the entry (if there is one) in entityHalo in jj.

         jj = entityIndex(db)%entryList(i,j,k)

         ! Now determine the situation we are dealing with.

         if(jj == 0) then

           ! Donor is an owned entity of the block. Store the
           ! appropriate info in tmpBuf. There are no periodic
           ! subfaces for this donor.

           tmpBuf(mm+1) = myID
           tmpBuf(mm+2) = db
           tmpBuf(mm+3) = i
           tmpBuf(mm+4) = j
           tmpBuf(mm+5) = k
           nPeriodic    = 0

         else

           ! Donor is also a halo, but its level of indirectness is at
           ! least one less than the one for which information is to
           ! be found. Still there are two possibilities. Either this
           ! is an internal block boundary halo or a physical boundary
           ! halo. In the former case the corresponding halo is stored,
           ! in the latter case the indices of the boundary halo are
           ! returned.

           if(entityHalo(jj)%donorProc == -1) then

             ! Physical boundary halo. Store the appropriate
             ! info in tmpBuf. No periodic subfaces for this donor.

              tmpBuf(mm+1) = myID
              tmpBuf(mm+2) = db
              tmpBuf(mm+3) = i
              tmpBuf(mm+4) = j
              tmpBuf(mm+5) = k
              nPeriodic    = 0

           else

             ! Internal block boundary halo. Store the appropriate
             ! info in tmpBuf, including the possible periodic
             ! subfaces.

              tmpBuf(mm+1) = entityHalo(jj)%donorProc
              tmpBuf(mm+2) = entityHalo(jj)%donorBlock
              tmpBuf(mm+3) = entityHalo(jj)%dI
              tmpBuf(mm+4) = entityHalo(jj)%dJ
              tmpBuf(mm+5) = entityHalo(jj)%dK

              nPeriodic = entityHalo(jj)%nPeriodicSubfaces
              do kk=1,nPeriodic
                tmpBuf(mm+5+kk) = entityHalo(jj)%periodicSubfaces(kk)
              enddo

           endif

         endif

         ! Fill the remaining part reserved for the periodic subfaces
         ! with 0's. A 0 indicates that no periodic subface is crossed.

         do kk=(nPeriodic+1),level
           tmpBuf(mm+5+kk) = 0
         enddo

         ! Update nn and mm for the next halo.

         nn = nn + 7
         mm = mm + nItemReturn

       enddo

       ! Copy the data from tmpBuf into recvBuf and delete tmpBuf
       ! afterwards.

       nn = nItemReturn*nHalos
       do i=1,nn
         recvBuf(i) = tmpBuf(i)
       enddo

       deallocate(tmpBuf, stat=ierr)
       if(ierr /= 0)                            &
         call terminate("findDonorsRecvBuffer", &
                        "Deallocation failure for tmpBuf")

       end subroutine findDonorsRecvBuffer

       !=================================================================

       subroutine storeHalosInList(buffer, bufSize, proc, level, &
                                   nItemReturn, entityHalo,      &
                                   entityIndex, iihalo)
!
!      ******************************************************************
!      *                                                                *
!      * storeHalosInList stores the halo info present in buf, which    *
!      * has been retreived from the given processor, in the correct    *
!      * place in entityHalo and entityIndex.                           *
!      *                                                                *
!      ******************************************************************
!
       use haloList
       use indirectHalo
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType), intent(in)    :: bufSize, proc
       integer(kind=intType), intent(in)    :: level, nItemReturn
       integer(kind=intType), intent(inout) :: iihalo

       integer(kind=intType), dimension(:), intent(in) :: buffer

       type(haloListType),  dimension(:), intent(inout) :: entityHalo
       type(indexListType), dimension(:), intent(inout) :: entityIndex
!
!      Local variables.
!
       integer :: ierr

       integer(kind=intType) :: i, j, k
       integer(kind=intType) :: ii, nn, blockID, iii
       integer(kind=intType) :: nPeriodic
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Store the start index (-1) for this processor in indHalo in nn.

       nn = nHaloPerLev(level-1) + nHaloPerProc(proc)

       ! Loop over the number of halo's stored in the buffer.

       do ii=1,bufSize

         iii = (ii-1)*nItemReturn

         ! Update the counters iihalo and nn.

         iihalo = iihalo +1
         nn     = nn +1

         ! Store the i,j,k indices and the block ID of the current
         ! halo a bit easier.

         blockID = indHalo(nn)%myBlock
         i       = indHalo(nn)%myI
         j       = indHalo(nn)%myJ
         k       = indHalo(nn)%myK

         ! Store the entry of entityHalo in the i,j,k indices
         ! of in entityIndex.

         entityIndex(blockID)%entryList(i,j,k) = iihalo

         ! Store the info of the current halo in entityHalo.

         entityHalo(iihalo)%myBlock = blockID

         entityHalo(iihalo)%myI = i
         entityHalo(iihalo)%myJ = j
         entityHalo(iihalo)%myK = k

         entityHalo(iihalo)%donorProc  = buffer(iii+1)
         entityHalo(iihalo)%donorBlock = buffer(iii+2)

         entityHalo(iihalo)%dI = buffer(iii+3)
         entityHalo(iihalo)%dJ = buffer(iii+4)
         entityHalo(iihalo)%dK = buffer(iii+5)

         ! Determine the number of periodic subfaces in the buffer.

         nPeriodic = 0
         do i=6,nItemReturn
           if(buffer(iii+i) > 0) nPeriodic = nPeriodic + 1
         enddo

         ! Check if the corresponding direct halo borders a periodic
         ! subface and update nPeriodic accordingly.

         j = indHalo(nn)%myDirectHalo
         nPeriodic = nPeriodic + entityHalo(j)%nPeriodicSubfaces

         ! If periodic subfaces are present for this halo, allocate
         ! the memory for periodicSubfaces and copy the data from
         ! both the buffer and the direct halo.

         if(nPeriodic > 0) then
           entityHalo(iihalo)%nPeriodicSubfaces = nPeriodic
           allocate(entityHalo(iihalo)%periodicSubfaces(nPeriodic), &
                    stat=ierr)
           if(ierr /= 0)                        &
             call terminate("storeHalosInList", &
                            "Memory allocation failure for &
                            &periodicSubfaces")
           nPeriodic = 0
           do i=6,nItemReturn
             if(buffer(iii+i) > 0) then
               nPeriodic = nPeriodic + 1
               entityHalo(iihalo)%periodicSubfaces(nPeriodic) = &
                                                           buffer(iii+i)
             endif
           enddo

           do i=1,entityHalo(j)%nPeriodicSubfaces
             nPeriodic = nPeriodic + 1
             entityHalo(iihalo)%periodicSubfaces(nPeriodic) = &
                                   entityHalo(j)%periodicSubfaces(i)
           enddo
         endif
           
       enddo

       ! Check in debug mode if the buffer size was correct.

       if( debug ) then
         if(nn /= nHaloPerLev(level-1) + nHaloPerProc(proc+1)) &
           call terminate("storeHalosInList",                  &
                          "Something wrong with buffer size")
       endif

       end subroutine storeHalosInList

!
!      ******************************************************************
!      *                                                                *
!      * File:          donorSearchCycle.f90                            *
!      * Author:        Steve Repsher                                   *
!      * Starting date: 04-22-2005                                      *
!      * Last modified: 10-16-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine donorSearchCycle(level, sps, nOrphan,   &
                                   commPattern, internal)
!
!      ******************************************************************
!      *                                                                *
!      * donorSearchCycle performs a cycle of a donor search process    *
!      * given the input communication patterns. It will attempt to use *
!      * the estimated donors to find the correct ones if they exist on *
!      * this processor. The failures (meaning those searches that may  *
!      * need to be restarted on another processor, or have crossed a   *
!      * boundary) will be removed from the communication patterns and  *
!      * tacked on to the halo list. The starting position to add to    *
!      * that list is given by the current number of actual orphans,    *
!      * such that all orphans are cumulatively kept throughout the     *
!      * possibly multiple cycles.                                      *
!      *                                                                *
!      ******************************************************************
!
       use block
       use boundaryList
       use communication
       use searchMod
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType), intent(in)    :: level, sps
       integer(kind=intType), intent(inout) :: nOrphan

       type(commType),         intent(inout) :: commPattern
       type(internalCommType), intent(inout) :: internal
!
!      Local variables.
!
       integer :: msize, procId, ierr, index, mcount
       integer, dimension(mpi_status_size) :: mstatus

       integer(kind=intType) :: i, j, k, ii, jj, flag, ihalo, nRecv
       integer(kind=intType) :: blk, nn, fi, fj, fk, n1to1HaloAdd
 
       integer(kind=intType), dimension(3) :: ind

       real(kind=realType), dimension(3)       :: xp
       real(kind=realType), dimension(nInterp) :: interp

       integer(kind=intType), dimension(:,:), allocatable :: sendBuf
       integer(kind=intType), dimension(:,:), allocatable :: recvBuf

       real(kind=realType), dimension(:,:), pointer :: coor

       logical, dimension(:), allocatable :: keep
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Allocate memory for the sending and receiving buffers.

       ii = commPattern%nSendCum(commPattern%nProcSend)
       jj = commPattern%nRecvCum(commPattern%nProcRecv)

       allocate(sendBuf(6,ii), recvBuf(6,jj), stat=ierr)
       if(ierr /= 0)                        &
         call terminate("donorSearchCycle", &
                        "Memory allocation failure for buffers")
 
       ! Loop over the processors i need to send to and try to find the
       ! correct stencil and interpolants for each point in the list.

       ii = 1
       jj = 0
       sends: do i=1,commPattern%nProcSend

         ! Set the pointer to the coordinates (stored in interp) and
         ! reallocate interp to store the donor weights.

         coor => commPattern%sendList(i)%interp

         j = commPattern%nSend(i)
         allocate(commPattern%sendList(i)%interp(j,nInterp), stat=ierr)
         if(ierr /= 0)                        &
           call terminate("donorSearchCycle", &
                          "Memory allocation failure for interp")

         ! Loop over each point in the list.

         k = 0
         do j=1,commPattern%nSend(i)

           ! Store the estimated block id, donor indices, and the
           ! coordinates of the receiver easier.

           blk = commPattern%sendList(i)%block(j)
           ind = commPattern%sendList(i)%indices(j,:)
           xp  = coor(j,:)

           ! Search for the correct donor indices and interpolants. The
           ! loop exists so that if the called routine exits with a flag
           ! = my process id, this means that the search wants to jump to
           ! another block on this processor so repeat the call.

           do
             call stencilSearch(level, sps, blk, ind, xp, interp, flag)
             if (flag /= myId) exit
           end do

           ! Perform some extra tasks based on the output flag of the
           ! stencil search.

           select case (flag)
             case (Success)

               ! Copy correct info back to the list.

               k = k + 1
               commPattern%sendList(i)%block(k)     = blk
               commPattern%sendList(i)%indices(k,:) = ind
               commPattern%sendList(i)%interp(k,:)  = interp

             case (StopAtBoco, HitMaxIter)

               ! The point is an orphan (can't find a donor). Convert
               ! both the block and indices of the terminal info to the
               ! CGNS domains. For a coarse level, the indices are taken
               ! as the fine grid indices for the corner of the coarse
               ! cell.

               do nn = level,2,-1
                 if (ind(1) > 1) &
                           ind(1) = flowDoms(blk,nn,1)%mgIFine(ind(1),1)
                 if (ind(2) > 1) &
                           ind(2) = flowDoms(blk,nn,1)%mgJFine(ind(2),1)
                 if (ind(3) > 1) &
                           ind(3) = flowDoms(blk,nn,1)%mgKFine(ind(3),1)
               end do

               ind(1) = ind(1) + flowDoms(blk,1,1)%iBegor - 2
               ind(2) = ind(2) + flowDoms(blk,1,1)%jBegor - 2
               ind(3) = ind(3) + flowDoms(blk,1,1)%kBegor - 2

               blk = flowDoms(blk,1,1)%cgnsBlockId

               ! Fill the buffer with the results for reporting purposes.

               jj = jj + 1
               sendBuf(1,  jj) = j
               sendBuf(2,  jj) = flag
               sendBuf(3,  jj) = blk
               sendBuf(4:6,jj) = ind

             case default 

               ! The remaining options are a restart or a bad donor
               ! stencil, in which case send back the results as is.

               jj = jj + 1
               sendBuf(1,  jj) = j
               sendBuf(2,  jj) = flag
               sendBuf(3,  jj) = blk
               sendBuf(4:6,jj) = ind

           end select

         end do

         ! Store the processor id and size then send the message.

         procID = commPattern%sendProc(i)
         msize  = 6*(jj - ii + 1)
         call mpi_isend(sendBuf(1,ii), msize, sumb_integer, procId, &
                        1, SUmb_comm_world, sendRequests(i), ierr)

         ! Update the sizes. The number removed from the list is the
         ! difference between the current size and the counter k.

         ii                       = commPattern%nSend(i) - k
         commPattern%nSend(i)     = k
         commPattern%nSendCum(i:) = commPattern%nSendCum(i:) - ii

         ! Update ii for the next processor.

         ii = jj + 1

         ! Deallocate the memory for the coordinates pointer.

         deallocate(coor, stat=ierr)
         if(ierr /= 0)                        &
           call terminate("donorSearchCycle", &
                          "Deallocation failure for coor")

       end do sends

       ! Post the nonblocking receives.

       ii = 1
       receives: do i=1,commPattern%nProcRecv

         ! Store the processor id and the size of the message, then
         ! post them.

         procID = commPattern%recvProc(i)
         msize  = 6*commPattern%nRecv(i)

         call mpi_irecv(recvBuf(1,ii), msize, sumb_integer, procId, &
                        1, SUmb_comm_world, recvRequests(i), ierr)

         ii = ii + commPattern%nRecv(i)

       end do receives

       ! Set the halo index to the current number of failures such
       ! that failures from the previous cycle are not overwritten.

       ihalo = nOrphan

       ! Initialize the number of 1-to-1 cells that will be communicated
       ! to the connecting block for addition to its fringe.

       n1to1HaloAdd = 0

       ! Repeat the process for the local communication.

       ! Set the pointer to the coordinates (stored in interp) and
       ! reallocate interp to store the donor weights.

       coor => internal%donorInterp

       j = internal%nCopy
       allocate(internal%donorInterp(j,nInterp), stat=ierr)
       if(ierr /= 0)                        &
         call terminate("donorSearchCycle", &
                        "Memory allocation failure for donorInterp")

       ! Loop over the cells in the list.

       k = 0
       localSearch: do i=1,internal%nCopy

         ! Store the estimated block id, donor indices, and the
         ! coordinates of the receiver easier.

         blk = internal%donorBlock(i)
         ind = internal%donorIndices(i,:)
         xp  = coor(i,:)

         ! Search for the correct donor indices and interpolants. The
         ! loop exists so that if the called routine exits with a flag
         ! = my process id, this means that the search wants to jump to
         ! another block on this processor so repeat the call.

         do
           call stencilSearch(level, sps, blk, ind, xp, interp, flag)
           if (flag /= myId) exit
         end do

         ! Store the boundary block and indices easier.

         nn = internal%haloBlock(i)
         fi = internal%haloIndices(i,1)
         fj = internal%haloIndices(i,2)
         fk = internal%haloIndices(i,3)

         ! Perform some extra tasks based on the output flag of the
         ! stencil search.

         if (flag == Success) then

           ! Search is complete. Copy correct info back to the list.

           k = k + 1
           internal%donorBlock(k)     = blk
           internal%donorIndices(k,:) = ind
           internal%donorInterp(k,:)  = interp

           ! Also copy the halo info in case k /= i.

           internal%haloBlock(k)     = internal%haloBlock(i)
           internal%haloIndices(k,:) = internal%haloIndices(i,:)

         else if (flag == BadDonor .and. &
                  fi >= 2 .and. fi <= flowDoms(nn,level,1)%il .and. &
                  fj >= 2 .and. fj <= flowDoms(nn,level,1)%jl .and. &
                  fk >= 2 .and. fk <= flowDoms(nn,level,1)%kl) then

           ! An owned boundary cell's donor has poor quality. Switch the
           ! cell to a field cell and add more fringe if necessary.

           call fringeToField(level, sps, nn, fi, fj, fk,        &
                              myID, blk, ind(1), ind(2), ind(3), &
                              ihalo, n1to1HaloAdd)

         else 

           ! Check to make sure the size of oversetHalo is sufficient.

           call checkSizeBoundaryList(ihalo, 1_intType)

           ! Cell was either orphaned or the search needs to be
           ! restarted on another processor. Either way, add its and
           ! the terminal donor info to the halo list.

           ihalo = ihalo + 1

           oversetHalo(ihalo)%myBlock = internal%haloBlock(i)
           oversetHalo(ihalo)%myI     = internal%haloIndices(i,1)
           oversetHalo(ihalo)%myJ     = internal%haloIndices(i,2)
           oversetHalo(ihalo)%myK     = internal%haloIndices(i,3)

           ! Make a distinction here between a restart and an orphan.
           ! In the former, the new donor proc is the flag and the
           ! level of indirectness is not important, while in the
           ! latter the donor processor is set to -1 to indicate an
           ! orphan, and the level of indirectness is set to the flag
           ! such that the failure reason may be reported. An exception
           ! is when the boundary cell is a halo and bad donor quality
           ! is flagged, in which case the receiving processor is kept.

           if (flag >= 0)  then ! Restart

             oversetHalo(ihalo)%donorProc = flag
             oversetHalo(ihalo)%levOfInd  = 0

           else if (flag == BadDonor) then ! Halo with bad donor

             oversetHalo(ihalo)%donorProc = -1
             oversetHalo(ihalo)%levOfInd  = myID
             nOrphan = nOrphan + 1

           else ! Orphan

             oversetHalo(ihalo)%donorProc = -1
             oversetHalo(ihalo)%levOfInd  = flag
             nOrphan = nOrphan + 1

             ! The point is an orphan (can't find a donor). Convert
             ! both the block and indices of the terminal info to the
             ! CGNS domains. For a coarse level, the indices are taken
             ! as the fine grid indices for the corner of the coarse
             ! cell.

             print *, "1"
             do nn = level,2,-1
               if (ind(1) > 1) &
                           ind(1) = flowDoms(blk,nn,1)%mgIFine(ind(1),1)
               if (ind(2) > 1) &
                           ind(2) = flowDoms(blk,nn,1)%mgJFine(ind(2),1)
               if (ind(3) > 1) &
                           ind(3) = flowDoms(blk,nn,1)%mgKFine(ind(3),1)
             end do

             ind(1) = ind(1) + flowDoms(blk,1,1)%iBegor - 2
             ind(2) = ind(2) + flowDoms(blk,1,1)%jBegor - 2
             ind(3) = ind(3) + flowDoms(blk,1,1)%kBegor - 2

             blk = flowDoms(blk,1,1)%cgnsBlockId
             print *, "2"

           end if

           ! Add the donor block, indices, and coordinates to 
           ! the halo list. 

           oversetHalo(ihalo)%donorBlock = blk
           oversetHalo(ihalo)%dI         = ind(1)
           oversetHalo(ihalo)%dJ         = ind(2)
           oversetHalo(ihalo)%dK         = ind(3)
           oversetHalo(ihalo)%interp     = xp

         end if
       enddo localSearch

       ! Update the size to the success counter k.

       internal%nCopy = k

       ! Deallocate the memory for the coordinates pointer.

       deallocate(coor, stat=ierr)
       if(ierr /= 0)                        &
         call terminate("donorSearchCycle", &
                        "Deallocation failure for coor")

       ! Aloocate the mask keep used to determine which parts of each
       ! receive list will be removed.

       ii = maxval(commPattern%nRecv)
       allocate(keep(ii), stat=ierr)
       if(ierr /= 0)                        &
         call terminate("donorSearchCycle", &
                        "Memory allocation failure for keep")

       ! Complete the nonblocking receives in an arbitrary sequence.

       msize = commPattern%nProcRecv
       completeRecvs: do i=1,commPattern%nProcRecv

         ! Wait for any of the requests and get the actual size.

         call mpi_waitany(msize, recvRequests, index, mstatus, ierr)
         call mpi_get_count(mstatus, sumb_integer, mcount, ierr)

         ! Find the starting point in the buffer for this message,
         ! and calculate the number of cells sent back.

         nRecv = mcount/6
         ii    = index
         jj    = commPattern%nRecvCum(ii-1)

         ! Check to make sure the size of oversetHalo is sufficient to
         ! store at list 1 halo for each receive.

         call checkSizeBoundaryList(ihalo, nRecv)

         ! Reset the mask keep to true and loop over the message.

         keep(1:commPattern%nRecv(ii)) = .true.

         do j = 1,nRecv

           ! Update the counter, store the boundary indices and block
           ! easier, and extract the list position and flag from the buffer.

           jj   = jj + 1
           k    = recvBuf(1,jj)
           flag = recvBuf(2,jj)

           keep(k) = .false.

           nn = commPattern%recvList(ii)%block(k)
           fi = commPattern%recvList(ii)%indices(k,1)
           fj = commPattern%recvList(ii)%indices(k,2)
           fk = commPattern%recvList(ii)%indices(k,3)

           ! Perform tasks based on the output flag of the search.

           if (flag == BadDonor .and. &
               fi >= 2 .and. fi <= flowDoms(nn,level,1)%il .and. &
               fj >= 2 .and. fj <= flowDoms(nn,level,1)%jl .and. &
               fk >= 2 .and. fk <= flowDoms(nn,level,1)%kl) then

             ! An owned boundary cell's donor has poor quality. Switch
             ! the cell to a field cell and add more fringe if necessary.

             call fringeToField(level, sps, nn, fi, fj, fk,   &
                                commPattern%recvProc(ii),     &
                                recvBuf(3,jj), recvBuf(4,jj), &
                                recvBuf(5,jj), recvBuf(6,jj), &
                                ihalo, n1to1HaloAdd)

           else

             ! Check to make sure the size of oversetHalo is sufficient.

             call checkSizeBoundaryList(ihalo, 1_intType)

             ! Cell was either orphaned or the search needs to be
             ! restarted on another processor. Either way, add its and
             ! the terminal donor info to the halo list.

             ihalo = ihalo + 1

             oversetHalo(ihalo)%myBlock = nn
             oversetHalo(ihalo)%myI     = fi
             oversetHalo(ihalo)%myJ     = fj
             oversetHalo(ihalo)%myK     = fk

             oversetHalo(ihalo)%donorBlock = recvBuf(3,jj)
             oversetHalo(ihalo)%dI         = recvBuf(4,jj)
             oversetHalo(ihalo)%dJ         = recvBuf(5,jj)
             oversetHalo(ihalo)%dK         = recvBuf(6,jj)

             ! Recompute the cell centroid using the external routine
             ! since the boundary could contain 2nd level halos.

             call cellCentroid(nn, level, sps, fi, fj, fk, &
                               xp(1), xp(2), xp(3))

             oversetHalo(ihalo)%interp(:) = eighth*xp

             ! Make a distinction here between a restart and an orphan.
             ! In the former, the new donor proc is the flag and the
             ! level of indirectness is not important, while in the
             ! latter the donor processor is set to -1 to indicate an
             ! orphan, and the level of indirectness is set to the flag
             ! such that the failure reason may be reported. An exception
             ! is when the boundary cell is a halo and bad donor quality
             ! is flagged, in which case the receiving processor is kept.

             if (flag >= 0)  then ! Restart

               oversetHalo(ihalo)%donorProc = flag
               oversetHalo(ihalo)%levOfInd  = 0

             else if (flag == BadDonor) then ! Halo with bad donor

               oversetHalo(ihalo)%donorProc = -1
               oversetHalo(ihalo)%levOfInd  = commPattern%recvProc(ii)
               nOrphan = nOrphan + 1

             else ! Orphan

               oversetHalo(ihalo)%donorProc = -1
               oversetHalo(ihalo)%levOfInd  = flag
               nOrphan = nOrphan + 1

             end if

           end if
         end do

         ! Remove the cells that received return messages from the list,
         ! and update the size from the counter k. Note nrecvCum is
         ! updated outside the loop so as not to corrupt the buffer
         ! position for the next message.

         k = 0
         do j = 1,commPattern%nRecv(ii)
           if (keep(j)) then
             k = k + 1
             commPattern%recvList(ii)%block(k) = &
             commPattern%recvList(ii)%block(j)
             commPattern%recvList(ii)%indices(k,:) = &
             commPattern%recvList(ii)%indices(j,:)
           end if
         end do

         commPattern%nRecv(ii) = k

       enddo completeRecvs

       ! Reform nRecvCum now that everything has been received.

       do i = 1,commPattern%nProcRecv
         commPattern%nRecvCum(i) = commPattern%nRecvCum(i-1) &
                                 + commPattern%nRecv(i)
       end do

       ! Complete the nonblocking sends.

       msize = commPattern%nProcSend
       do i=1,commPattern%nProcSend
         call mpi_waitany(msize, sendRequests, index, mstatus, ierr)
       enddo

       ! Deallocate the memory for the sending and receiving buffers
       ! and the keep array.
 
       deallocate(sendBuf, recvBuf, keep, stat=ierr)
       if(ierr /= 0)                        &
         call terminate("donorSearchCycle", &
                        "Deallocation failure for buffers and keep")

       ! Some of the fringe added due to bad donors may have been 1-to-1
       ! halos; therefore, scan for these and communicate the donor info
       ! to the connecting blocks.

       call scan1to1HalosForFringe(level, sps, n1to1HaloAdd, ihalo)

       ! Set the new number of halos left to the current halo index.

       nHaloOver = ihalo

       end subroutine donorSearchCycle

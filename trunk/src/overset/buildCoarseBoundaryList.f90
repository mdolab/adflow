!
!      ******************************************************************
!      *                                                                *
!      * File:          buildCoarseBoundaryList.f90                     *
!      * Author:        Steve Repsher                                   *
!      * Starting date: 04-09-2005                                      *
!      * Last modified: 08-26-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine buildCoarseBoundaryList(level, sps, nHaloBndryAdd)
!
!      ******************************************************************
!      *                                                                *
!      * buildCoarseBoundaryList builds the halo (or boundary) list     *
!      * for the coarse grid levels using the output of the creation    *
!      * routine. The block-to-block communication is treated first so  *
!      * that no reallocation of the list will be necessary. The output *
!      * of this routine is a filled boundary list with all estimates   *
!      * for the donor information and the centroid coordinates taking  *
!      * the place of the interpolants to pass to the donor processor.  *
!      * Finally, the flowDom boundary indices are discarded after the  *
!      * copy and will be redetermined later; however the right amount  *
!      * of memory is allocated for the overset info.                   *
!      *                                                                *
!      ******************************************************************
!
       use block
       use boundaryList
       use communication
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType), intent(in) :: level, sps, nHaloBndryAdd
!
!      Local variables.
!
       integer :: size, procId, ierr, index, source
       integer, dimension(mpi_status_size) :: status

       integer(kind=intType) :: i, j, ii, jj, nn, ihalo
       integer(kind=intType) :: d1, i1, j1, k1, d2, i2, j2, k2
 
       integer(kind=intType), dimension(nDom) :: bndryAdds
 
       integer(kind=intType), allocatable :: sendBufInt(:,:)
       integer(kind=intType), allocatable :: recvBufInt(:,:)
       integer(kind=intType), allocatable :: nrecv(:)

       real(kind=realType) :: xc, yc, zc
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Allocate the memory for the sending buffer, whose size is
       ! determined by yhe input nHaloBndryAdd * 6 (i.e. for each cell
       ! send the list position, and estimates for the donor proc,
       ! block, and indices). The +1 is needed since 0 size messages may
       ! be sent after the buffer is maxed out.

       allocate(sendBufInt(6, nHaloBndryAdd + 1), stat=ierr)
       if(ierr /= 0)                               &
         call terminate("buildCoarseBoundaryList", &
                        "Memory allocation failure for send buffer")
 
       ! Loop over the receiving processors and lists for the 1st layer
       ! halos. These are the halo indices and so search for values > 10
       ! as set by the createCoarseBoundary routine. These are the
       ! cells that will need to be part of the connecting block's
       ! boundary.
 
       ii = 1
       jj = 0
       sends: do i=1,commPatternCell_1st(level)%nProcRecv
         do j=1,commPatternCell_1st(level)%nrecv(i)

           ! Store the block id and the indices of the halo
           ! a bit easier.

           d1 = commPatternCell_1st(level)%recvList(i)%block(j)
           i1 = commPatternCell_1st(level)%recvList(i)%indices(j,1)
           j1 = commPatternCell_1st(level)%recvList(i)%indices(j,2)
           k1 = commPatternCell_1st(level)%recvList(i)%indices(j,3)
 
           ! Check the iblank value if this cell is to be communicated.
           ! Note that the createCoarseBoundary routine set the iblank
           ! to 10 + nearest bndry index.
 
           if (flowDoms(d1,level,sps)%iblank(i1,j1,k1) > 10) then

             ! Compute 8 times the centroid coordinates.

             xc = sum(flowDoms(d1,level,sps)%x(i1-1:i1,j1-1:j1,k1-1:k1,1))
             yc = sum(flowDoms(d1,level,sps)%x(i1-1:i1,j1-1:j1,k1-1:k1,2))
             zc = sum(flowDoms(d1,level,sps)%x(i1-1:i1,j1-1:j1,k1-1:k1,3))
 
             ! As these cells were added over holes, determine the index
             ! for the nearest neighbor and redefine the indices.
 
             nn = flowDoms(d1,level,sps)%iblank(i1,j1,k1) - 10
             i1 = flowDoms(d1,level,sps)%ibndry(1,nn)
             j1 = flowDoms(d1,level,sps)%ibndry(2,nn)
             k1 = flowDoms(d1,level,sps)%ibndry(3,nn)

             ! Put the current list position in the buffer so the
             ! receiving processor knows the cell to add/check.
 
             jj = jj + 1
             sendBufInt(1,jj) = j

             ! Estimate the donor info to send and store in the buffer.
 
             call estimateCoarseDonor(d1, level, sps, i1, j1, k1, &
                                      xc, yc, zc,                 &
                                      sendBufInt(2,jj),           &
                                      sendBufInt(3,jj),           &
                                      sendBufInt(4,jj),           &
                                      sendBufInt(5,jj),           &
                                      sendBufInt(6,jj))
 
           end if
         enddo

         ! Compute the size and send the data.
 
         procID = commPatternCell_1st(level)%recvProc(i)
         size = 6*(jj - ii + 1)
         call mpi_isend(sendBufInt(1,ii), size, sumb_integer, procId, &
                       procId, SUmb_comm_world, sendRequests(i),   &
                       ierr)

         ! Set ii to jj+1 for the next processor.

         ii = jj + 1

       end do sends
 
       ! Loop over the number of blocks stored on this processor and
       ! count the total # of overset halos currently stored.

       nHaloOver = 0
       do nn = 1,nDom
         nHaloOver = nHaloOver + flowDoms(nn,level,sps)%nCellsOverset
       end do
 
       ! Add the difference between total sent and nHaloBndryAdd. This
       ! is the possible addition due to internal communication.
 
       nHaloOver = nHaloOver + nHaloBndryAdd - jj

       ! Allocate memory to number of cells being received from each
       ! possible processor.

       ii = commPatternCell_1st(level)%nProcSend
       allocate(nRecv(ii), stat=ierr)
       if(ierr /= 0)                               &
         call terminate("buildCoarseBoundaryList", &
                        "Memory allocation failure for nrecv")

       ! Now block and probe for messages from all the processors I
       ! can receive from. Note these are the sending processors for
       ! this list.
 
       probeSizes: do i=1,commPatternCell_1st(level)%nProcSend

         ! Probe for a message in order. MPI_any_source cannot be used
         ! here because the messages are not received in this loop.

         source = commPatternCell_1st(level)%sendProc(i)
         call mpi_probe(source, myId, SUmb_comm_world, &
                        status, ierr)

         ! Get the message size and store the numer of cells in the
         ! message (i.e. divide by 6).
 
         call mpi_get_count(status, sumb_integer, size, ierr)
 
         nRecv(i) = size/6
 
       end do probeSizes
 
       ! Allocate memory for the receiving buffer. Blocking receives
       ! used so the buffer can be recycled. Don't allow 0 size.
 
       ii = max(maxval(nRecv), 1_intType)
       allocate(recvBufInt(6,ii), stat=ierr)
       if(ierr /= 0)                               &
         call terminate("buildCoarseBoundaryList", &
                        "Memory allocation failure for recv buffer")

       ! Add the total cells being received to the total number of
       ! overset cells on this processor.

       nHaloOver = nHaloOver + sum(nRecv)
 
       ! Allocate memory for the halo list for this processor. Note
       ! that the actual size may be less than this because some
       ! communication could be repetitive. This would happen if a
       ! connecting block has already added the same cell to its
       ! boundary. This is checked by examining the receiving block's
       ! iblank value.
 
       allocate(oversetHalo(nHaloOver), stat=ierr)
       if(ierr /= 0)                              &
        call terminate("buildCoarseBoundaryList", &
                       "Memory allocation failure for oversetHalo")

       ! Initialize the index for building the list to 0.

       ihalo = 0

       ! Initialize the counter for cells added to each local block's
       ! boundary. This will be added to nCellsOverset later.

       bndryAdds = 0
 
       ! Loop over the receiving processors again.
 
       receives: do i=1,commPatternCell_1st(level)%nProcSend
 
         ! Store the source processor a bit easier and set the size
         ! the message.
 
         source = commPatternCell_1st(level)%sendProc(i)
         size   = 6*nRecv(i)
 
         ! Receive the message as it has already arrived a blocking
         ! receive can be used.

         call mpi_recv(recvBufInt, size, sumb_integer, source, &
                       myId, SUmb_comm_world, status, ierr)

         ! Loop over the number of cells being received.

         do j=1,nRecv(i)

           ! Get the list position from the buffer and then store the
           ! block and indices of the cell to be checked.

           jj = recvBufInt(1,j)
           d1 = commPatternCell_1st(level)%sendList(i)%block(jj)
           i1 = commPatternCell_1st(level)%sendList(i)%indices(jj,1)
           j1 = commPatternCell_1st(level)%sendList(i)%indices(jj,2)
           k1 = commPatternCell_1st(level)%sendList(i)%indices(jj,3)

           ! If the iblank value for this cell still indicates it is
           ! a hole, then add it to the list.

           hole: if (flowDoms(d1,level,sps)%iblank(i1,j1,k1) <= 0) then

             ! Reset the iblank to indicate this cell is now part of 
             ! the boundary.

             flowDoms(d1,level,sps)%iblank(i1,j1,k1) = 9

             ! Update the boundary counter and store this cell's info.

             ihalo = ihalo + 1
             oversetHalo(ihalo)%myBlock = d1
             oversetHalo(ihalo)%myI     = i1
             oversetHalo(ihalo)%myJ     = j1
             oversetHalo(ihalo)%myK     = k1

             ! The donor info was sent in the buffer so extract it.

             oversetHalo(ihalo)%donorProc  = recvBufInt(2,j)
             oversetHalo(ihalo)%donorBlock = recvBufInt(3,j)

             oversetHalo(ihalo)%dI = recvBufInt(4,j)
             oversetHalo(ihalo)%dJ = recvBufInt(5,j)
             oversetHalo(ihalo)%dK = recvBufInt(6,j)

             ! Update the # of overset cells added for this block.

             bndryAdds(d1) = bndryAdds(d1) + 1

             ! Now allocate the interpolants for the boundary list and
             ! store the cell centroid coordinates.

             allocate(oversetHalo(ihalo)%interp(3), stat=ierr)
             if(ierr /= 0)                              &
              call terminate("buildCoarseBoundaryList", &
                             "Memory allocation failure for interp")

             xc = sum(flowDoms(d1,level,sps)%x(i1-1:i1,j1-1:j1,k1-1:k1,1))
             yc = sum(flowDoms(d1,level,sps)%x(i1-1:i1,j1-1:j1,k1-1:k1,2))
             zc = sum(flowDoms(d1,level,sps)%x(i1-1:i1,j1-1:j1,k1-1:k1,3))

             oversetHalo(ihalo)%interp(1) = eighth*xc
             oversetHalo(ihalo)%interp(2) = eighth*yc
             oversetHalo(ihalo)%interp(3) = eighth*zc

             ! Set levOfInd to 0 so it doesn't corrupt the sort.

             oversetHalo(ihalo)%levOfInd = 0

           end if hole
         end do

       end do receives

       ! Now treat the internal communication with the same process.

       localComm: do i=1,internalCell_1st(level)%ncopy

         ! Store the block and the indices of the donor a bit easier.

         d1 = internalCell_1st(level)%donorBlock(i)
         i1 = internalCell_1st(level)%donorIndices(i,1)
         j1 = internalCell_1st(level)%donorIndices(i,2)
         k1 = internalCell_1st(level)%donorIndices(i,3)

         ! Idem for the halo's.

         d2 = internalCell_1st(level)%haloBlock(i)
         i2 = internalCell_1st(level)%haloIndices(i,1)
         j2 = internalCell_1st(level)%haloIndices(i,2)
         k2 = internalCell_1st(level)%haloIndices(i,3)

         ! Check that both the halo has the > 10 signal and the receiver
         ! is still a hole.

         checks: if (flowDoms(d1,level,sps)%iblank(i1,j1,k1) <= 0 .and. &
                     flowDoms(d2,level,sps)%iblank(i2,j2,k2) > 10) then

           ! Reset the iblank to indicate this cell is now part of
           ! the boundary.

           flowDoms(d1,level,sps)%iblank(i1,j1,k1) = 9

           ! As these cells were added over holes, determine the index
           ! for the nearest neighbor and redefine the indices.

           nn = flowDoms(d2,level,sps)%iblank(i2,j2,k2) - 10
           i2 = flowDoms(d2,level,sps)%ibndry(1,nn)
           j2 = flowDoms(d2,level,sps)%ibndry(2,nn)
           k2 = flowDoms(d2,level,sps)%ibndry(3,nn)

           ! Update the boundary counter and store this cell's info.

           ihalo = ihalo + 1
           oversetHalo(ihalo)%myBlock = d1
           oversetHalo(ihalo)%myI     = i1
           oversetHalo(ihalo)%myJ     = j1
           oversetHalo(ihalo)%myK     = k1

           ! Update the # of overset cells added for this block.

           bndryAdds(d1) = bndryAdds(d1) + 1

           ! Now allocate the interpolants for the boundary list and
           ! store the cell centroid coordinates.

           allocate(oversetHalo(ihalo)%interp(3), stat=ierr)
           if(ierr /= 0)                              &
            call terminate("buildCoarseBoundaryList", &
                           "Memory allocation failure for interp")

           xc = sum(flowDoms(d1,level,sps)%x(i1-1:i1,j1-1:j1,k1-1:k1,1))
           yc = sum(flowDoms(d1,level,sps)%x(i1-1:i1,j1-1:j1,k1-1:k1,2))
           zc = sum(flowDoms(d1,level,sps)%x(i1-1:i1,j1-1:j1,k1-1:k1,3))

           oversetHalo(ihalo)%interp(1) = eighth*xc
           oversetHalo(ihalo)%interp(2) = eighth*yc
           oversetHalo(ihalo)%interp(3) = eighth*zc

           ! Estimate the donor info to send and store in the list.

           call estimateCoarseDonor(d2, level, sps, i2, j2, k2,    &
                                    xc, yc, zc,                    &
                                    oversetHalo(ihalo)%donorProc,  &
                                    oversetHalo(ihalo)%donorBlock, &
                                    oversetHalo(ihalo)%dI,         &
                                    oversetHalo(ihalo)%dJ,         &
                                    oversetHalo(ihalo)%dK)

           ! Set levOfInd to 0 so it doesn't corrupt the sort.

           oversetHalo(ihalo)%levOfInd = 0

         end if checks
       enddo localComm

       ! Finally, add all of the boundary already created.

       domains: do nn=1,nDom

         boundary: do i=1,flowDoms(nn,level,sps)%ncellsOverset

           ! Store the indices a bit easier.

           i1 = flowDoms(nn,level,sps)%ibndry(1,i)
           j1 = flowDoms(nn,level,sps)%ibndry(2,i)
           k1 = flowDoms(nn,level,sps)%ibndry(3,i)

           ! Update the boundary counter and store this cell's info.

           ihalo = ihalo + 1
           oversetHalo(ihalo)%myBlock = nn
           oversetHalo(ihalo)%myI     = i1
           oversetHalo(ihalo)%myJ     = j1
           oversetHalo(ihalo)%myK     = k1

           ! Now allocate the interpolants for the boundary list to
           ! store the cell centroid coordinates.

           allocate(oversetHalo(ihalo)%interp(3), stat=ierr)
           if(ierr /= 0)                              &
            call terminate("buildCoarseBoundaryList", &
                           "Memory allocation failure for interp")

           ! Compute 8 times the cell centroid using the external
           ! routine since the boundary could contain 2nd level halos.

           call cellCentroid(nn, level, sps, i1, j1, k1, xc, yc, zc)

           oversetHalo(ihalo)%interp(1) = eighth*xc
           oversetHalo(ihalo)%interp(2) = eighth*yc
           oversetHalo(ihalo)%interp(3) = eighth*zc

           ! If the current index is >= istartOverHole, then it won't
           ! contain any of the finer grid boundary so switch to the
           ! nearest neighbor

           if (i >= blockBndry(nn)%istartOverHole) then
             ii = blockBndry(nn)%nearestBndry(i)
             i1 = flowDoms(nn,level,sps)%ibndry(1,ii)
             j1 = flowDoms(nn,level,sps)%ibndry(2,ii)
             k1 = flowDoms(nn,level,sps)%ibndry(3,ii)
           end if

           ! Estimate the donor info to send and store in the list.

           call estimateCoarseDonor(nn, level, sps, i1, j1, k1,    &
                                    xc, yc, zc,                    &
                                    oversetHalo(ihalo)%donorProc,  &
                                    oversetHalo(ihalo)%donorBlock, &
                                    oversetHalo(ihalo)%dI,         &
                                    oversetHalo(ihalo)%dJ,         &
                                    oversetHalo(ihalo)%dK)

           ! Set levOfInd to 0 so it doesn't corrupt the sort.

           oversetHalo(ihalo)%levOfInd = 0

         end do boundary

         ! Release the memory for the block boundary indices. This will
         ! be rebuilt after the donor searches.

         deallocate(flowDoms(nn,level,sps)%ibndry,  &
                    blockBndry(nn)%nearestBndry, stat=ierr)
         if(ierr /= 0)                              &
          call terminate("buildCoarseBoundaryList", &
                      "Deallocation failure for ibndry & nearestBndry")

         ! Calculate the new total number of boundary cells for this
         ! domain. No allocation is made for the flowDom info here.

         flowDoms(nn,level,sps)%nCellsOverset = &
         flowDoms(nn,level,sps)%nCellsOverset + bndryAdds(nn)

       end do domains

       ! Allocate the remaineder of the halo list's interp arrays
       ! since they may be used later, and set the actual value of
       ! nHaloOver to ihalo since the upper bound may be different
       ! than what was used.

       do i=ihalo+1,nHaloOver
         allocate(oversetHalo(i)%interp(3), stat=ierr)
         if(ierr /= 0)                               &
           call terminate("buildCoarseBoundaryList", &
                          "Memory allocation failure for interp")
       end do

       nHaloOver = ihalo

       ! Complete the nonblocking sends.

       size = commPatternCell_1st(level)%nProcRecv
       do i=1,commPatternCell_1st(level)%nProcRecv
         call mpi_waitany(size, sendRequests, index, status, ierr)
       enddo

       ! Deallocate the memory for the sending and receiving buffers.

       deallocate(sendBufInt, recvBufInt, nrecv, stat=ierr)
       if(ierr /= 0)                               &
         call terminate("buildCoarseBoundaryList", &
                        "Deallocation failure for buffers")

       end subroutine buildCoarseBoundaryList

!
!      ******************************************************************
!      *                                                                *
!      * File:          scan1to1HalosForFringe.f90                      *
!      * Author:        Steve Repsher                                   *
!      * Starting date: 08-24-2005                                      *
!      * Last modified: 08-27-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine scan1to1HalosForFringe(level, sps, n1to1HaloAdd, ihalo)
!
!      ******************************************************************
!      *                                                                *
!      * scan1to1HalosForFringe scans over all of the 1-to-1 halo cells *
!      * given in the communication list and looks for cells that have  *
!      * been flagged for communication to the connecting block, i.e.   *
!      * the connecting block will be told to add its corresponding     *
!      * to its fringe. The cells to be communicated are flagged by     *
!      * iblank value from which an index into the donorInfo array may  *
!      * be determined (stoed in the boundaryList module).              *
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
       integer(kind=intType), intent(in)    :: level, sps, n1to1HaloAdd
       integer(kind=intType), intent(inout) :: ihalo
!
!      Local variables.
!
       integer :: msize, procId, ierr, index, mcount
       integer, dimension(mpi_status_size) :: mstatus

       integer(kind=intType) :: i, j, k, ii, jj, nRecv, nn
       integer(kind=intType) :: cb, ci, cj, ck, hb, hi, hj, hk
 
       real(kind=realType) :: xc, yc, zc

       integer(kind=intType), dimension(:,:), allocatable :: sendBuf
       integer(kind=intType), dimension(:,:), allocatable :: recvBuf
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Allocate memory for the sending and receiving buffers. Note
       ! that the maximum possible size is allocated in order to avoid
       ! repeated reallocation and blocking communication.

       ii = commPatternCell_2nd(level)%nProcSend
       ii = commPatternCell_2nd(level)%nsendCum(ii)

       allocate(sendBuf(6,n1to1HaloAdd+1), recvBuf(6,ii), stat=ierr)
       if(ierr /= 0)                              &
         call terminate("scan1to1HalosForFringe", &
                        "Memory allocation failure for buffers")
 
       ! Loop over the receiving processors and lists for the 1to1 halos.
       ! These are the halo indices and so search for values > 10 as set
       ! by the fringeToField routine.

       ii = 1
       jj = 0
       sends: do i=1,commPatternCell_2nd(level)%nProcRecv
         do j=1,commPatternCell_2nd(level)%nRecv(i)

           ! Store the block id and the indices of the halo
           ! a bit easier.

           hb = commPatternCell_2nd(level)%recvList(i)%block(j)
           hi = commPatternCell_2nd(level)%recvList(i)%indices(j,1)
           hj = commPatternCell_2nd(level)%recvList(i)%indices(j,2)
           hk = commPatternCell_2nd(level)%recvList(i)%indices(j,3)

           ! Check the iblank value if this cell is to be communicated,
           ! which is indicated by an iblank > 10.

           if (flowDoms(hb,level,sps)%iblank(hi,hj,hk) > 10) then

             ! Find the index into the donorInfo array and copy it and
             ! the list position to the buffer.

             nn = flowDoms(hb,level,sps)%iblank(hi,hj,hk) - 10

             jj = jj + 1
             sendBuf(1,  jj) = j
             sendBuf(2:6,jj) = donorInfo(:,nn)

           end if
         end do

         ! Store the processor id and size then send the message.

         procID = commPatternCell_2nd(level)%recvProc(i)
         msize  = 6*(jj - ii + 1)

         call mpi_isend(sendBuf(1,ii), msize, sumb_integer, procId, &
                        1, SUmb_comm_world, sendRequests(i), ierr)

         ! Update ii for the next processor.

         ii = jj + 1

       end do sends

       ! Post the nonblocking receives. Note the sending lists are
       ! used since they are the owned cells.

       ii = 1
       receives: do i=1,commPatternCell_2nd(level)%nProcSend

         ! Store the processor id and the size of the message, then
         ! post them.

         procID = commPatternCell_2nd(level)%sendProc(i)
         msize  = 6*commPatternCell_2nd(level)%nSend(i)

         call mpi_irecv(recvBuf(1,ii), msize, sumb_integer, procId, &
                        1, SUmb_comm_world, recvRequests(i), ierr)

         ii = ii + commPatternCell_2nd(level)%nSend(i)

       end do receives

       ! Check to make sure the size of oversetHalo is sufficient.

       call checkSizeBoundaryList(ihalo, n1to1HaloAdd - jj)

       ! Loop over the local communication.

       localComm: do i=1,internalCell_2nd(level)%nCopy

         ! Store the halo and new fringe block and indices easier.

         hb = internalCell_2nd(level)%haloBlock(i)
         hi = internalCell_2nd(level)%haloIndices(i,1)
         hj = internalCell_2nd(level)%haloIndices(i,2)
         hk = internalCell_2nd(level)%haloIndices(i,3)

         cb = internalCell_2nd(level)%donorBlock(i)
         ci = internalCell_2nd(level)%donorIndices(i,1)
         cj = internalCell_2nd(level)%donorIndices(i,2)
         ck = internalCell_2nd(level)%donorIndices(i,3)

         ! Check that both the iblank of the halo indicates it to be
         ! made a fringe and that the corresponding connecting block has
         ! not already made it a fringe or field cell.

         if (flowDoms(hb,level,sps)%iblank(hi,hj,hk) > 10 .and. &
             flowDoms(cb,level,sps)%iblank(ci,cj,ck) <= 0) then

           ! Add the new cell and indices to the boundary list.

           ihalo = ihalo + 1

           oversetHalo(ihalo)%myBlock = cb
           oversetHalo(ihalo)%myI     = ci
           oversetHalo(ihalo)%myJ     = cj
           oversetHalo(ihalo)%myK     = ck

           ! Find the index into the donorInfo array and copy the
           ! donor info to the list.

           nn = flowDoms(hb,level,sps)%iblank(hi,hj,hk) - 10

           oversetHalo(ihalo)%donorProc  = donorInfo(1,nn)
           oversetHalo(ihalo)%donorBlock = donorInfo(2,nn)
           oversetHalo(ihalo)%dI         = donorInfo(3,nn)
           oversetHalo(ihalo)%dJ         = donorInfo(4,nn)
           oversetHalo(ihalo)%dK         = donorInfo(5,nn)

           ! Compute the cell centroid using the external routine
           ! since the boundary cell could be a 2nd level halo.

           call cellCentroid(cb, level, sps, ci, cj, ck, xc, yc, zc)

           oversetHalo(ihalo)%interp(1) = eighth*xc
           oversetHalo(ihalo)%interp(2) = eighth*yc
           oversetHalo(ihalo)%interp(3) = eighth*zc

           ! Set levOfInd to 0 as to not corrupt the sort.

           oversetHalo(ihalo)%levOfInd = 0

           ! Update the number of overset cells in the fringe block and
           ! set the iblank to indicate it's now part of the fringe.

           flowDoms(cb,level,sps)%nCellsOverset = &
           flowDoms(cb,level,sps)%nCellsOverset + 1

           flowDoms(cb,level,sps)%iblank(ci,cj,ck) = 9

         end if
       enddo localComm

       ! Complete the nonblocking receives in an arbitrary sequence.

       msize = commPatternCell_2nd(level)%nProcSend
       completeRecvs: do i=1,commPatternCell_2nd(level)%nProcSend

         ! Wait for any of the requests and get the actual size.

         call mpi_waitany(msize, recvRequests, index, mstatus, ierr)
         call mpi_get_count(mstatus, sumb_integer, mcount, ierr)

         ! Find the starting point in the buffer for this message,
         ! and calculate the number of cells sent back.

         nRecv = mcount/6
         ii    = index
         jj    = commPatternCell_2nd(level)%nSendCum(ii-1)

         ! Check to make sure the size of oversetHalo is sufficient.

         call checkSizeBoundaryList(ihalo, nRecv)

         ! Loop over the message.

         do j = 1,nRecv

           ! Update the counter and get the list position to check.

           jj   = jj + 1
           k    = recvBuf(1,jj)

           ! Store the block and indices a bit easier.

           cb = commPatternCell_2nd(level)%sendList(ii)%block(k)
           ci = commPatternCell_2nd(level)%sendList(ii)%indices(k,1)
           cj = commPatternCell_2nd(level)%sendList(ii)%indices(k,2)
           ck = commPatternCell_2nd(level)%sendList(ii)%indices(k,3)

           ! Check if the iblank to make sure the cell is a hole since
           ! it could have already been made a field or fringe in the
           ! current cycle.

           if (flowDoms(cb,level,sps)%iblank(ci,cj,ck) <= 0) then

             ! Add all info to the boundary list.

             ihalo = ihalo + 1

             oversetHalo(ihalo)%myBlock = cb
             oversetHalo(ihalo)%myI     = ci
             oversetHalo(ihalo)%myJ     = cj
             oversetHalo(ihalo)%myK     = ck

             oversetHalo(ihalo)%donorProc  = recvBuf(2,jj)
             oversetHalo(ihalo)%donorBlock = recvBuf(3,jj)
             oversetHalo(ihalo)%dI         = recvBuf(4,jj)
             oversetHalo(ihalo)%dJ         = recvBuf(5,jj)
             oversetHalo(ihalo)%dK         = recvBuf(6,jj)

             ! Compute the cell centroid using the external routine
             ! since the boundary cell could be a 2nd level halo.

             call cellCentroid(cb, level, sps, ci, cj, ck, xc, yc, zc)

             oversetHalo(ihalo)%interp(1) = eighth*xc
             oversetHalo(ihalo)%interp(2) = eighth*yc
             oversetHalo(ihalo)%interp(3) = eighth*zc

             ! Set levOfInd to 0 as to not corrupt the sort.

             oversetHalo(ihalo)%levOfInd = 0

             ! Update the number of overset cells in the fringe block and
             ! set the iblank to indicate it's now part of the fringe.

             flowDoms(cb,level,sps)%nCellsOverset = &
             flowDoms(cb,level,sps)%nCellsOverset + 1

             flowDoms(cb,level,sps)%iblank(ci,cj,ck) = 9

           end if
         end do

       enddo completeRecvs

       ! Complete the nonblocking sends.

       msize = commPatternCell_2nd(level)%nProcRecv
       do i=1,commPatternCell_2nd(level)%nProcRecv
         call mpi_waitany(msize, sendRequests, index, mstatus, ierr)
       enddo

       ! Deallocate the memory for the sending and receiving buffers.
 
       deallocate(sendBuf, recvBuf, stat=ierr)
       if (ierr /= 0)                             &
         call terminate("scan1to1HalosForFringe", &
                        "Deallocation failure for buffers")

       end subroutine scan1to1HalosForFringe

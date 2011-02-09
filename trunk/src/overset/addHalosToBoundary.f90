!
!      ******************************************************************
!      *                                                                *
!      * File:          addHalosToBoundary.f90                          *
!      * Author:        Steve Repsher                                   *
!      * Starting date: 07-09-2005                                      *
!      * Last modified: 10-15-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine addHalosToBoundary(level, sps)
!
!      ******************************************************************
!      *                                                                *
!      * addHalosToBoundary adds the halo cells that are part of the    *
!      * connecting block's overset boundary to its own boundary. This  *
!      * done for consistency and convergence reasons. In the solver,   *
!      * the halos (1-to-1 and sliding) are communicated first so that  *
!      * donor stencils which straddle the boundary will interpolate    *
!      * using the most recent data. However, this means that halos     *
!      * which are part of the boundary are not the most recent, which  *
!      * will cause different histories with number of processors and   *
!      * has been shown to cause convergence problems. This routine     *
!      * adds to the overset communication pattern such that, where the *
!      * need arises, a donor stencil will interpolate twice: once to   *
!      * the owned cell of a block, and once to the corresponding halo  *
!      * of the connecting block.                                       *
!      *                                                                *
!      ******************************************************************
!
       use block
       use boundaryList
       use communication
       use inputOverset
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType), intent(in) :: level, sps
!
!      Local variables.
!
       integer :: size, procId, ierr, index, source
       integer, dimension(mpi_status_size) :: status

       integer(kind=intType) :: i, j, ii, jj, nn, ihalo
       integer(kind=intType) :: d1, i1, j1, k1, d2, nInterp
 
       integer(kind=intType), dimension(nDom) :: nHaloAdd
 
       integer(kind=intType), allocatable :: sendBuf(:,:), recvBuf(:,:)
       integer(kind=intType), allocatable :: nRecv(:)

       real(kind=realType), allocatable :: sendBufInt(:,:)
       real(kind=realType), allocatable :: recvBufInt(:,:)

       type(commType)         :: commHaloAdds
       type(internalCommType) :: internalHaloAdds
!
!      Interfaces.
!
       interface
         subroutine reallocateInteger2(intArray,            &
                                        newSize1, newSize2, &
                                        oldSize1, oldSize2, &
                                        alwaysFreeMem)
           use precision
           implicit none

           integer(kind=intType), dimension(:,:), pointer :: intArray
           integer(kind=intType), intent(in) :: newSize1, newSize2, &
                                                   oldSize1, oldSize2
           logical, intent(in) :: alwaysFreeMem
         end subroutine reallocateInteger2
       end interface
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Determine the number of interpolants that will be communicated
       ! per boundary cell.

       nInterp = nDonorWeights(oversetInterpTypeCoarse)
       if (level == 1) nInterp = nDonorWeights(oversetInterpType) 

       ! Allocate the memory for the sending buffers. Just assume that
       ! the entire list is to be communicated.

       ii = commPatternCell_2nd(level)%nProcSend
       ii = commPatternCell_2nd(level)%nSendCum(ii)

       allocate(sendBuf(6,ii), sendBufInt(nInterp,ii), stat=ierr)
       if(ierr /= 0)                          &
         call terminate("addHalosToBoundary", &
                        "Memory allocation failure for send buffers")
 
       ! Loop over the sending processors and lists for the 2nd level
       ! cell halo communication pattern.
 
       ii = 1
       jj = 0
       sends: do i=1,commPatternCell_2nd(level)%nProcSend
         do j=1,commPatternCell_2nd(level)%nSend(i)

           ! Store the block id and the indices a bit easier.

           d1 = commPatternCell_2nd(level)%sendList(i)%block(j)
           i1 = commPatternCell_2nd(level)%sendList(i)%indices(j,1)
           j1 = commPatternCell_2nd(level)%sendList(i)%indices(j,2)
           k1 = commPatternCell_2nd(level)%sendList(i)%indices(j,3)
 
           ! Check if the donor is part of the boundary.
 
           if (flowDoms(d1,level,sps)%iblank(i1,j1,k1) >= 10) then

             ! Put the current list position in the buffer so the
             ! receiving processor knows the cell to add.
 
             jj = jj + 1
             sendBuf(1,jj) = j

             ! Compute the index into the boundary information from the
             ! iblank value and store all info in the buffers.
 
             nn = flowDoms(d1,level,sps)%iblank(i1,j1,k1)/10

             sendBuf(2  ,jj) = flowDoms(d1,level,sps)%neighProcOver(nn)
             sendBuf(3  ,jj) = flowDoms(d1,level,sps)%neighBlockOver(nn)
             sendBuf(4:6,jj) = flowDoms(d1,level,sps)%idonor(:,nn)

             sendBufInt(:,jj) = flowDoms(d1,level,sps)%overint(:,nn)

           end if
         enddo

         ! Compute the size and send the data.
 
         procID = commPatternCell_2nd(level)%sendProc(i)
         size = 6*(jj - ii + 1)

         call mpi_isend(sendBuf(1,ii), size, sumb_integer, procId, &
                        1, SUmb_comm_world, sendRequests(i), ierr)

         size = nInterp*(jj - ii + 1)
         call mpi_isend(sendBufInt(1,ii), size, sumb_real, procId, &
                        2, SUmb_comm_world, recvRequests(i), ierr)

         ! Set ii to jj+1 for the next processor.

         ii = jj + 1

       end do sends

       ! Initialize the number of halos to add to the boundary to 0.

       nHaloAdd = 0

       ! Loop over the local communication and count the number of donors
       ! that are also part of the boundary.

       localCount: do i=1,internalCell_2nd(level)%nCopy

         ! Store the block and the indices of the donor a bit easier.

         d1 = internalCell_2nd(level)%donorBlock(i)
         i1 = internalCell_2nd(level)%donorIndices(i,1)
         j1 = internalCell_2nd(level)%donorIndices(i,2)
         k1 = internalCell_2nd(level)%donorIndices(i,3)

         ! Check if the donor is part of the boundary.

         if (flowDoms(d1,level,sps)%iblank(i1,j1,k1) >= 10) then
           d2 = internalCell_2nd(level)%haloBlock(i)
           nHaloAdd(d2) = nHaloAdd(d2) + 1
         end if
       end do localCount
 
       ! Allocate memory to store the number of cells being received from
       ! each possible processor.

       ii = commPatternCell_2nd(level)%nProcRecv
       allocate(nRecv(ii), stat=ierr)
       if(ierr /= 0)                          &
         call terminate("addHalosToBoundary", &
                        "Memory allocation failure for nRecv")

       ! Now block and probe for messages from all the processors I
       ! can receive from.
 
       probeSizes: do i=1,commPatternCell_2nd(level)%nProcRecv

         ! Probe for a message in order. MPI_any_source cannot be used
         ! here because the messages are not received in this loop.

         source = commPatternCell_2nd(level)%recvProc(i)
         call mpi_probe(source, 1, SUmb_comm_world, &
                       status, ierr)

         ! Get the message size and store the numer of cells in the
         ! message (i.e. divide by 6).
 
         call mpi_get_count(status, sumb_integer, size, ierr)
 
         nRecv(i) = size/6
 
       end do probeSizes
 
       ! Allocate memory for the receiving buffer. Blocking receives
       ! are used so the buffer can be recycled. Don't allow 0 size.
       ! Also allocate memory for the halo list whose size will be the
       ! sum of those being received plus local communication.
 
       ii = max(maxval(nRecv), 1_intType)
       nHaloOver = sum(nHaloAdd) + sum(nRecv)

       allocate(recvBuf(6,ii), recvBufInt(nInterp,ii), &
                oversetHalo(nHaloOver), stat=ierr)
       if(ierr /= 0)                          &
         call terminate("addHalosToBoundary", &
                        "Memory allocation failure for recv buffers &
                        &and halo list")

       ! Initialize the index for building the list to 0.

       ihalo = 0

       ! Loop over the local communication again, this time to build
       ! the list.

       localAdd: do i=1,internalCell_2nd(level)%nCopy

         ! Store the block and the indices of the donor a bit easier.

         d1 = internalCell_2nd(level)%donorBlock(i)
         i1 = internalCell_2nd(level)%donorIndices(i,1)
         j1 = internalCell_2nd(level)%donorIndices(i,2)
         k1 = internalCell_2nd(level)%donorIndices(i,3)

         ! Check if the donor is part of the boundary.

         if (flowDoms(d1,level,sps)%iblank(i1,j1,k1) >= 10) then

           ! Update the counter and add the halo information to the list.

           ihalo = ihalo + 1
           oversetHalo(ihalo)%myBlock = &
                                    internalCell_2nd(level)%haloBlock(i)

           oversetHalo(ihalo)%myI = &
                                internalCell_2nd(level)%haloIndices(i,1)
           oversetHalo(ihalo)%myJ = &
                                internalCell_2nd(level)%haloIndices(i,2)
           oversetHalo(ihalo)%myK = &
                                internalCell_2nd(level)%haloIndices(i,3)

           ! Determine the index into the boundary information and copy
           ! it to the list.

           nn = flowDoms(d1,level,sps)%iblank(i1,j1,k1)/10

           oversetHalo(ihalo)%donorProc = &
                             flowDoms(d1,level,sps)%neighProcOver(nn)
           oversetHalo(ihalo)%donorBlock = &
                             flowDoms(d1,level,sps)%neighBlockOver(nn)

           oversetHalo(ihalo)%dI = flowDoms(d1,level,sps)%idonor(1,nn)
           oversetHalo(ihalo)%dJ = flowDoms(d1,level,sps)%idonor(2,nn)
           oversetHalo(ihalo)%dK = flowDoms(d1,level,sps)%idonor(3,nn)

           allocate(oversetHalo(ihalo)%interp(nInterp), stat=ierr)
           if(ierr /= 0)                           &
              call terminate("addHalosToBoundary", &
                             "Memory allocation failure for interp")

           oversetHalo(ihalo)%interp = &
                                    flowDoms(d1,level,sps)%overint(:,nn)
           oversetHalo(ihalo)%levOfInd = 0

         end if
       end do localAdd

       ! Loop over the receiving processors again.
 
       receives: do i=1,commPatternCell_2nd(level)%nProcRecv
 
         ! Receive the messages. Blocking receives are used since they
         ! were already probed above.

         source = commPatternCell_2nd(level)%recvProc(i)
         size   = 6*nRecv(i)
         call mpi_recv(recvBuf, size, sumb_integer, source, &
                       1, SUmb_comm_world, status, ierr)

         size   = nInterp*nRecv(i)
         call mpi_recv(recvBufInt, size, sumb_real, source, &
                       2, SUmb_comm_world, status, ierr)

         ! Loop over the number of cells being received.

         do j=1,nRecv(i)

           ! Get the list position from the buffer and then store the
           ! block and indices of the halo cell in the list.

           jj = recvBuf(1,j)
           ihalo = ihalo + 1

           oversetHalo(ihalo)%myBlock = &
                    commPatternCell_2nd(level)%recvList(i)%block(jj)
           oversetHalo(ihalo)%myI = &
                    commPatternCell_2nd(level)%recvList(i)%indices(jj,1)
           oversetHalo(ihalo)%myJ = &
                    commPatternCell_2nd(level)%recvList(i)%indices(jj,2)
           oversetHalo(ihalo)%myK = &
                    commPatternCell_2nd(level)%recvList(i)%indices(jj,3)

           ! Extract the donor information from the buffers.

           oversetHalo(ihalo)%donorProc  = recvBuf(2,j)
           oversetHalo(ihalo)%donorBlock = recvBuf(3,j)
           oversetHalo(ihalo)%dI         = recvBuf(4,j)
           oversetHalo(ihalo)%dJ         = recvBuf(5,j)
           oversetHalo(ihalo)%dK         = recvBuf(6,j)

           allocate(oversetHalo(ihalo)%interp(nInterp), stat=ierr)
           if(ierr /= 0)                           &
              call terminate("addHalosToBoundary", &
                             "Memory allocation failure for interp")

           oversetHalo(ihalo)%interp = recvBufInt(:,j)
           oversetHalo(ihalo)%levOfInd = 0

           ! Update the counter for the boundary adds.

           d2 = commPatternCell_2nd(level)%recvList(i)%block(jj)
           nHaloAdd(d2) = nHaloAdd(d2) + 1

         end do
       end do receives

       ! Loop over the domains and reallocate the arrays to store the
       ! boundary indices. Note the donor information is not needed.
       ! Set nHaloAdd to the number of owned overset cells so it can
       ! be used as a counter below.

       domains: do nn = 1,nDom

         ii = flowDoms(nn,level,sps)%nCellsOverset &
            + flowDoms(nn,level,sps)%nOrphans
         jj = ii + nHaloAdd(nn)
         nHaloAdd(nn) = ii
         flowDoms(nn,level,sps)%nCellsOversetAll = jj

         call reallocateInteger2(flowDoms(nn,level,sps)%ibndry, &
                                 3_intType, jj, 3_intType, ii, .true.)
       end do domains

       ! Complete the nonblocking sends.

       size = commPatternCell_2nd(level)%nProcSend
       do i=1,commPatternCell_2nd(level)%nProcSend
         call mpi_waitany(size, sendRequests, index, status, ierr)
         call mpi_waitany(size, recvRequests, index, status, ierr)
       enddo

       ! Deallocate the memory for the sending and receiving buffers.

       deallocate(sendBuf, recvBuf, nRecv, &
                  sendBufInt, recvBufInt, stat=ierr)
       if(ierr /= 0)                          &
         call terminate("addHalosToBoundary", &
                        "Deallocation failure for buffers")

       ! Sort the boundary list in increasing order, build the 
       ! communication structures, merge them to the existing overset
       ! structures, and release the memory for the temporary ones.

       call qsortHaloListType(oversetHalo, nHaloOver)
       call finalCommStructures(oversetHalo, nHaloOver, commHaloAdds, &
                                internalHaloAdds, nInterp)
       call mergeComm(commPatternOverset(level,sps), &
                      internalOverset(level,sps),    &
                      commHaloAdds, internalHaloAdds)
       call releaseCommPattern(commHaloAdds)
       call releaseInternalComm(internalHaloAdds)

       ! Loop over the halo list and add the indices to the appropriate
       ! domain's boundary list. Also deallocate the interpolants.

       do i = 1,nHaloOver
         nn = oversetHalo(i)%myBlock
         nHaloAdd(nn) = nHaloAdd(nn) + 1

         flowDoms(nn,level,sps)%ibndry(1,nHaloAdd(nn)) = &
                                                      oversetHalo(i)%myI
         flowDoms(nn,level,sps)%ibndry(2,nHaloAdd(nn)) = &
                                                      oversetHalo(i)%myJ
         flowDoms(nn,level,sps)%ibndry(3,nHaloAdd(nn)) = &
                                                      oversetHalo(i)%myK

         deallocate(oversetHalo(i)%interp, stat=ierr)
         if(ierr /= 0)                          &
           call terminate("addHalosToBoundary", &
                          "Deallocation failure for interp")
       end do

       ! Deallocate the memory for the halo list.

       deallocate(oversetHalo, stat=ierr)
       if(ierr /= 0)                          &
         call terminate("addHalosToBoundary", &
                        "Deallocation failure for oversetHalo")

       end subroutine addHalosToBoundary

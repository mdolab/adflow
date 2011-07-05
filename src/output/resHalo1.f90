!
!      ******************************************************************
!      *                                                                *
!      * File:          resHalo1.f90                                    *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 04-15-2003                                      *
!      * Last modified: 06-12-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine resHalo1(level, start, end)
!
!      ******************************************************************
!      *                                                                *
!      * resHalo1 determines the residuals in the 1st layer of halo     *
!      * cells by applying both the boundary conditions and the         *
!      * exchange. The halo values are needed for post processing       *
!      * reasons.                                                       *
!      *                                                                *
!      ******************************************************************
!
       use BCTypes
       use blockPointers
       use communication
       use inputTimeSpectral
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType) :: level, start, end
!
!      Local variables.
!
       integer :: size, procID, ierr, index
       integer, dimension(mpi_status_size) :: status

       integer(kind=intType) :: nVar, sps
       integer(kind=intType) :: ii, jj, mm, nn, i, j, k, l
       integer(kind=intType) :: dd1, ii1, jj1, kk1, dd2, ii2, jj2, kk2

       real(kind=realType), pointer, dimension(:,:,:) :: ddw1, ddw2
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Determine the number of variables per cell to be sent.

       nVar = max(0_intType,(end - start + 1))
       if(nVar == 0) return

       ! Loop over the spectral solutions and local blocks to apply
       ! the boundary conditions for the residual.

       spectralLoop: do sps=1,nTimeIntervalsSpectral
         domains: do mm=1,nDom

           ! Set the pointers for this block.

           call setPointers(mm, level, sps)

           ! Loop over the boundary condition subfaces of this block
           ! and apply a neumann boundary condition. I know that this is
           ! not entirely correct for some boundary conditions, symmetry,
           ! solid wall, but this is not so important.

           bocos: do nn=1,nBocos

             ! Set the pointer for ddw1 and ddw2, depending on the block
             ! face on which the subface is located.

             select case (BCFaceID(nn))
               case (iMin)
                 ddw1 => dw(1, 1:,1:,:); ddw2 => dw(2, 1:,1:,:)
               case (iMax)
                 ddw1 => dw(ie,1:,1:,:); ddw2 => dw(il,1:,1:,:)
               case (jMin)
                 ddw1 => dw(1:,1, 1:,:); ddw2 => dw(1:,2, 1:,:)
               case (jMax)
                 ddw1 => dw(1:,je,1:,:); ddw2 => dw(1:,jl,1:,:)
               case (kMin)
                 ddw1 => dw(1:,1:,1, :); ddw2 => dw(1:,1:,2, :)
               case (kMax)
                 ddw1 => dw(1:,1:,ke,:); ddw2 => dw(1:,1:,kl,:)
             end select

             ! Loop over the cell range of the subface.

             do j=BCData(nn)%jcBeg, BCData(nn)%jcEnd
               do i=BCData(nn)%icBeg, BCData(nn)%icEnd
                 do l=start,end
                   ddw1(i,j,l) = ddw2(i,j,l)
                 enddo
               enddo
             enddo

           enddo bocos
         enddo domains

         ! Send the variables. The data is first copied into
         ! the send buffer after which the buffer is sent asap.

         ii = 1
         sends: do i=1,commPatternCell_1st(level)%nProcSend

           ! Store the processor id and the size of the message
           ! a bit easier.

           procID = commPatternCell_1st(level)%sendProc(i)
           size    = nVar*commPatternCell_1st(level)%nsend(i)

           ! Copy the data in the correct part of the send buffer.

           jj = ii
           do j=1,commPatternCell_1st(level)%nsend(i)

             ! Store the block id and the indices of the donor
             ! a bit easier.

             dd1 = commPatternCell_1st(level)%sendList(i)%block(j)
             ii1 = commPatternCell_1st(level)%sendList(i)%indices(j,1)
             jj1 = commPatternCell_1st(level)%sendList(i)%indices(j,2)
             kk1 = commPatternCell_1st(level)%sendList(i)%indices(j,3)

             ! Copy the given range of the residuals for this cell
             ! in the buffer. Update the counter jj accordingly.

             do k=start,end
               sendBuffer(jj) = flowDoms(dd1,level,sps)%dw(ii1,jj1,kk1,k)
               jj = jj + 1
             enddo

           enddo

           ! Send the data.

           call mpi_isend(sendBuffer(ii), size, sumb_real, procID,  &
                          procID, SUmb_comm_world, sendRequests(i), &
                          ierr)

           ! Set ii to jj for the next processor.

           ii = jj

         enddo sends

         ! Post the nonblocking receives.

         ii = 1
         receives: do i=1,commPatternCell_1st(level)%nProcRecv

           ! Store the processor id and the size of the message
           ! a bit easier.

           procID = commPatternCell_1st(level)%recvProc(i)
           size    = nVar*commPatternCell_1st(level)%nrecv(i)

           ! Post the receive.

           call mpi_irecv(recvBuffer(ii), size, sumb_real, procID, &
                          myID, SUmb_comm_world, recvRequests(i), ierr)

           ! And update ii.

           ii = ii + size

         enddo receives

         ! Copy the local data.

         localCopy: do i=1,internalCell_1st(level)%ncopy

           ! Store the block and the indices of the donor a bit easier.

           dd1 = internalCell_1st(level)%donorBlock(i)
           ii1 = internalCell_1st(level)%donorIndices(i,1)
           jj1 = internalCell_1st(level)%donorIndices(i,2)
           kk1 = internalCell_1st(level)%donorIndices(i,3)

           ! Idem for the halo's.

           dd2 = internalCell_1st(level)%haloBlock(i)
           ii2 = internalCell_1st(level)%haloIndices(i,1)
           jj2 = internalCell_1st(level)%haloIndices(i,2)
           kk2 = internalCell_1st(level)%haloIndices(i,3)

           ! Copy the given range of residuals.

           do k=start,end
             flowDoms(dd2,level,sps)%dw(ii2,jj2,kk2,k) = &
                  flowDoms(dd1,level,sps)%dw(ii1,jj1,kk1,k)
           enddo

         enddo localCopy

         ! Complete the nonblocking receives in an arbitrary sequence and
         ! copy the variables from the buffer into the halo's.

         size = commPatternCell_1st(level)%nProcRecv
         completeRecvs: do i=1,commPatternCell_1st(level)%nProcRecv

           ! Complete any of the requests.

           call mpi_waitany(size, recvRequests, index, status, ierr)

           ! Copy the data just arrived in the halo's.

           ii = index
           jj = nVar*commPatternCell_1st(level)%nrecvCum(ii-1) +1
           do j=1,commPatternCell_1st(level)%nrecv(ii)

             ! Store the block and the indices of the halo a bit easier.

             dd2 = commPatternCell_1st(level)%recvList(ii)%block(j)
             ii2 = commPatternCell_1st(level)%recvList(ii)%indices(j,1)
             jj2 = commPatternCell_1st(level)%recvList(ii)%indices(j,2)
             kk2 = commPatternCell_1st(level)%recvList(ii)%indices(j,3)

             ! Copy the residuals.

             do k=start,end
               flowDoms(dd2,level,sps)%dw(ii2,jj2,kk2,k) = recvBuffer(jj)
               jj = jj + 1
             enddo

           enddo

         enddo completeRecvs

         ! Complete the nonblocking sends.

         size = commPatternCell_1st(level)%nProcSend
         do i=1,commPatternCell_1st(level)%nProcSend
           call mpi_waitany(size, sendRequests, index, status, ierr)
         enddo

       enddo spectralLoop

       end subroutine resHalo1

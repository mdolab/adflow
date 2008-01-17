!
!      ******************************************************************
!      *                                                                *
!      * File:          f1SST.f90                                       *
!      * Author:        Georgi Kalitzin, Edwin van der Weide,           *
!      *                Steve Repsher                                   *
!      * Starting date: 08-19-2003                                      *
!      * Last modified: 07-03-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine f1SST
!
!      ******************************************************************
!      *                                                                *
!      * f1SST computes the blending function f1 in both the owned      *
!      * cells and the first layer of halo's. The result is stored in   *
!      * dw(:,:,:,if1SST). For the computation of f1 also the cross     *
!      * diffusion term is needed. This is stored in dw(:,:,:,icd) such *
!      * that it can be used in SSTSolve later on.                      *
!      *                                                                *
!      ******************************************************************
!
       use BCTypes
       use blockPointers
       use inputTimeSpectral
       use iteration
       use turbMod
       implicit none
!
!      Local variables.
!
       integer(kind=intType) :: sps, nn, mm, i, j, k

       real(kind=realType) :: t1, t2, arg1
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! First part. Compute the values of the blending function f1
       ! for each block and spectral solution.

       spectralLoop: do sps=1,nTimeIntervalsSpectral
         domains: do mm=1,nDom

           ! Set the pointers to this block.

           call setPointers(mm, currentLevel, sps)

           ! Set the pointers for f1 and kwCD to the correct entries
           ! in dw which are currently not used.

           f1   => dw(1:,1:,1:,if1SST)
           kwCD => dw(1:,1:,1:,icd)

           ! Compute the cross diffusion term.

           call kwCDterm

           ! Compute the blending function f1 for all owned cells.

           do k=2,kl
             do j=2,jl
               do i=2,il

                 t1 = sqrt(w(i,j,k,itu1)) &
                    / (0.09_realType*w(i,j,k,itu2)*d2Wall(i,j,k))
                 t2 = 500.0_realType*rlv(i,j,k) &
                    / (w(i,j,k,irho)*w(i,j,k,itu2)*d2Wall(i,j,k)**2)
                 t1 = max(t1,t2)
                 t2 = two*w(i,j,k,itu1)&
                    / (max(eps,kwCD(i,j,k))*d2Wall(i,j,k)**2)

                 arg1      = min(t1,t2)
                 f1(i,j,k) = tanh(arg1**4)

               enddo
             enddo
           enddo

           ! Loop over the boundary conditions to set f1 in the boundary
           ! halo's. A Neumann boundary condition is used for all BC's.

           bocos: do nn=1,nBocos

             ! Determine the face on which this subface is located, loop
             ! over its range and copy f1. Although the range may include
             ! indirect halo's which are not computed, this is no problem,
             ! because in SSTSolve only direct halo's are used.

             select case (BCFaceID(nn))

               case (iMin)
                 do k=kcBeg(nn),kcEnd(nn)
                   do j=jcBeg(nn),jcEnd(nn)
                     f1(1,j,k) = f1(2,j,k)
                   enddo
                 enddo

!              ==========================================================

               case (iMax)

                 do k=kcBeg(nn),kcEnd(nn)
                   do j=jcBeg(nn),jcEnd(nn)
                     f1(ie,j,k) = f1(il,j,k)
                   enddo
                 enddo

!              ==========================================================

               case (jMin)

                 do k=kcBeg(nn),kcEnd(nn)
                   do i=icBeg(nn),icEnd(nn)
                     f1(i,1,k) = f1(i,2,k)
                   enddo
                 enddo

!              ==========================================================

               case (jMax)

                 do k=kcBeg(nn),kcEnd(nn)
                   do i=icBeg(nn),icEnd(nn)
                     f1(i,je,k) = f1(i,jl,k)
                   enddo
                 enddo

!              ==========================================================

               case (kMin)

                 do j=jcBeg(nn),jcEnd(nn)
                   do i=icBeg(nn),icEnd(nn)
                     f1(i,j,1) = f1(i,j,2)
                   enddo
                 enddo

!              ==========================================================

               case (kMax)

                 do j=jcBeg(nn),jcEnd(nn)
                   do i=icBeg(nn),icEnd(nn)
                     f1(i,j,ke) = f1(i,j,kl)
                   enddo
                 enddo

             end select

           enddo bocos

         enddo domains
       enddo spectralLoop

       ! Exchange the values of f1.

       call exchangeF1SST1to1
       call exchangeF1SSTSliding
       call exchangeF1SSTOverset

       end subroutine f1SST

!      ==================================================================

       subroutine exchangeF1SST1to1
!
!      ******************************************************************
!      *                                                                *
!      * exchangeF1SST1to1 communicates the 1st layer of halo values    *
!      * for the blending function f1 of the SST model for 1 to 1       *
!      * matching halo's. This variable is stored in dw(:,:,:,if1SST).  *
!      *                                                                *
!      ******************************************************************
!
       use block
       use communication
       use inputTimeSpectral
       use iteration
       implicit none
!
!      Local variables.
!
       integer :: size, procID, ierr, index
       integer, dimension(mpi_status_size) :: status

       integer(kind=intType) :: i, j, ii, jj, sps, ll
       integer(kind=intType) :: d1, i1, j1, k1, d2, i2, j2, k2
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Easier storage of the current mg level.

       ll = currentLevel

       ! Loop over the number of spectral solutions.

       spectralModes: do sps=1,nTimeIntervalsSpectral

         ii = 1
         sends: do i=1,commPatternCell_1st(ll)%nProcSend

           ! Store the processor id and the size of the message
           ! a bit easier.

           procID = commPatternCell_1st(ll)%sendProc(i)
           size   = commPatternCell_1st(ll)%nSend(i)

           ! Copy the data in the correct part of the send buffer.

           jj = ii
           do j=1,commPatternCell_1st(ll)%nSend(i)

             ! Store the block id and the indices of the donor a
             !  bit easier.

             d1 = commPatternCell_1st(ll)%sendList(i)%block(j)
             i1 = commPatternCell_1st(ll)%sendList(i)%indices(j,1)
             j1 = commPatternCell_1st(ll)%sendList(i)%indices(j,2)
             k1 = commPatternCell_1st(ll)%sendList(i)%indices(j,3)

             ! Store the value of f1 in the send buffer. Note that the
             ! level is 1 and not ll (= currentLevel).

             sendBuffer(jj) = flowDoms(d1,1,sps)%dw(i1,j1,k1,if1SST)
             jj = jj + 1

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
         receives: do i=1,commPatternCell_1st(ll)%nProcRecv

           ! Store the processor id and the size of the message
           ! a bit easier.

           procID = commPatternCell_1st(ll)%recvProc(i)
           size   = commPatternCell_1st(ll)%nRecv(i)

           ! Post the receive.

           call mpi_irecv(recvBuffer(ii), size, sumb_real, procID, &
                          myID, SUmb_comm_world, recvRequests(i),  &
                          ierr)

           ! And update ii.

           ii = ii + size

         enddo receives

         ! Copy the local data.

         localCopy: do i=1,internalCell_1st(ll)%ncopy

           ! Store the block and the indices of the donor a bit easier.

           d1 = internalCell_1st(ll)%donorBlock(i)
           i1 = internalCell_1st(ll)%donorIndices(i,1)
           j1 = internalCell_1st(ll)%donorIndices(i,2)
           k1 = internalCell_1st(ll)%donorIndices(i,3)

           ! Idem for the halo's.

           d2 = internalCell_1st(ll)%haloBlock(i)
           i2 = internalCell_1st(ll)%haloIndices(i,1)
           j2 = internalCell_1st(ll)%haloIndices(i,2)
           k2 = internalCell_1st(ll)%haloIndices(i,3)

           ! Copy the values. Note that level is 1 and not
           ! ll (= currentLevel).

           flowDoms(d2,1,sps)%dw(i2,j2,k2,if1SST) = &
                   flowDoms(d1,1,sps)%dw(i1,j1,k1,if1SST)

         enddo localCopy

         ! Complete the nonblocking receives in an arbitrary sequence and
         ! copy the variables from the buffer into the halo's.

         size = commPatternCell_1st(ll)%nProcRecv
         completeRecvs: do i=1,commPatternCell_1st(ll)%nProcRecv

           ! Complete any of the requests.

           call mpi_waitany(size, recvRequests, index, status, ierr)

           ! Copy the data just arrived in the halo's.

           ii = index
           jj = commPatternCell_1st(ll)%nRecvCum(ii-1) +1
           do j=1,commPatternCell_1st(ll)%nRecv(ii)

             ! Store the block and the indices of the halo a bit easier.

             d2 = commPatternCell_1st(ll)%recvList(ii)%block(j)
             i2 = commPatternCell_1st(ll)%recvList(ii)%indices(j,1)
             j2 = commPatternCell_1st(ll)%recvList(ii)%indices(j,2)
             k2 = commPatternCell_1st(ll)%recvList(ii)%indices(j,3)

             ! And copy the data in the appropriate place in dw. Note
             ! that level == 1 and not ll (= currentLevel).

             flowDoms(d2,1,sps)%dw(i2,j2,k2,if1SST) = recvBuffer(jj)
             jj = jj + 1

           enddo

         enddo completeRecvs

         ! Complete the nonblocking sends.

         size = commPatternCell_1st(ll)%nProcSend
         do i=1,commPatternCell_1st(ll)%nProcSend
           call mpi_waitany(size, sendRequests, index, status, ierr)
         enddo

       enddo spectralModes

       end subroutine exchangeF1SST1to1

!      ==================================================================

       subroutine exchangeF1SSTSliding
!
!      ******************************************************************
!      *                                                                *
!      * exchangeF1SSTSliding communicates the 1st layer of halo        *
!      * values for the blending function f1 of the SST model for       *
!      * sliding mesh halo's.                                           *
!      *                                                                *
!      ******************************************************************
!
       use block
       use commSliding
       use communication
       use inputTimeSpectral
       use iteration
       use paramTurb
       implicit none
!
!      Local variables.
!
       integer :: size, procID, ierr, index
       integer, dimension(mpi_status_size) :: status

       integer(kind=intType) :: i, j, ii, jj, sps, ll
       integer(kind=intType) :: d1, i1, j1, k1, d2, i2, j2, k2

       real(kind=realType) :: alp
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Easier storage of the current mg level.

       ll = currentLevel

       ! Loop over the number of spectral solutions.

       spectralModes: do sps=1,nTimeIntervalsSpectral

         ii = 1
         sends: do i=1,commSlidingCell_1st(ll,sps)%nProcSend

           ! Store the processor id and the size of the message
           ! a bit easier.

           procID = commSlidingCell_1st(ll,sps)%sendProc(i)
           size   = commSlidingCell_1st(ll,sps)%nSend(i)

           ! Copy the data in the correct part of the send buffer.

           jj = ii
           do j=1,commSlidingCell_1st(ll,sps)%nSend(i)

             ! Store the block id and the indices of the donor
             ! a bit easier.

             d1 = commSlidingCell_1st(ll,sps)%sendList(i)%block(j)
             i1 = commSlidingCell_1st(ll,sps)%sendList(i)%indices(j,1)
             j1 = commSlidingCell_1st(ll,sps)%sendList(i)%indices(j,2)
             k1 = commSlidingCell_1st(ll,sps)%sendList(i)%indices(j,3)

             ! Copy the data in the send buffer. Note that level == 1 and
             ! not currentLevel (== ll) for dw.

             sendBuffer(jj) = flowDoms(d1,1,sps)%dw(i1,j1,k1,if1SST)
             jj = jj + 1

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
         receives: do i=1,commSlidingCell_1st(ll,sps)%nProcRecv

           ! Store the processor id and the size of the message
           ! a bit easier.

           procID = commSlidingCell_1st(ll,sps)%recvProc(i)
           size   = commSlidingCell_1st(ll,sps)%nRecv(i)

           ! Post the receive.

           call mpi_irecv(recvBuffer(ii), size, sumb_real, procID, &
                          myID, SUmb_comm_world, recvRequests(i),  &
                          ierr)

           ! And update ii.

           ii = ii + size

         enddo receives

         ! Initialize the sliding mesh halo's to zero.

         initHalos: do i=1,intSlidingCell_1st(ll,sps)%nslidingHalos

           ! Store the block and the indices of the halo a bit easier.

           d2 = intSlidingCell_1st(ll,sps)%slidingHaloList%block(i)
           i2 = intSlidingCell_1st(ll,sps)%slidingHaloList%indices(i,1)
           j2 = intSlidingCell_1st(ll,sps)%slidingHaloList%indices(i,2)
           k2 = intSlidingCell_1st(ll,sps)%slidingHaloList%indices(i,3)

           ! And set f1 to zero. Again level == 1 for dw.

           flowDoms(d2,1,sps)%dw(i2,j2,k2,if1SST) = zero

         enddo initHalos

         ! Interpolation with the data stored locally.

         localData: do i=1,intSlidingCell_1st(ll,sps)%ncopy

           ! Store the block and the indices of the donor a bit easier.

            d1 = intSlidingCell_1st(ll,sps)%donorList%block(i)
            i1 = intSlidingCell_1st(ll,sps)%donorList%indices(i,1)
            j1 = intSlidingCell_1st(ll,sps)%donorList%indices(i,2)
            k1 = intSlidingCell_1st(ll,sps)%donorList%indices(i,3)

            ! Store the block and the indices of the halo a bit easier.

            d2 = intSlidingCell_1st(ll,sps)%haloList%block(i)
            i2 = intSlidingCell_1st(ll,sps)%haloList%indices(i,1)
            j2 = intSlidingCell_1st(ll,sps)%haloList%indices(i,2)
            k2 = intSlidingCell_1st(ll,sps)%haloList%indices(i,3)

            ! Store the weight a bit easier.

            alp = intSlidingCell_1st(ll,sps)%weight(i)

            ! Update f1 of the halo. Note that level == 1 for dw.

            flowDoms(d2,1,sps)%dw(i2,j2,k2,if1SST) =     &
                flowDoms(d2,1,sps)%dw(i2,j2,k2,if1SST) + &
                     alp*flowDoms(d1,1,sps)%dw(i1,j1,k1,if1SST)

         enddo localData

         ! Complete the nonblocking receives in an arbitrary sequence and
         ! use the variables from the buffer to interpolate the halo's.

         size = commSlidingCell_1st(ll,sps)%nProcRecv
         completeRecvs: do i=1,commSlidingCell_1st(ll,sps)%nProcRecv

           ! Complete any of the requests.

           call mpi_waitany(size, recvRequests, index, status, ierr)

           ! Update the halo's using the data just arrived.

           ii = index

           do j=1,commSlidingCell_1st(ll,sps)%recvList(ii)%ncopy

             ! Store the block and the indices of the halo as well as
             ! the interpolation weight and the starting index in the
             ! receive buffer a bit easier.

             d2  = commSlidingCell_1st(ll,sps)%recvList(ii)%block(j)
             i2  = commSlidingCell_1st(ll,sps)%recvList(ii)%indices(j,1)
             j2  = commSlidingCell_1st(ll,sps)%recvList(ii)%indices(j,2)
             k2  = commSlidingCell_1st(ll,sps)%recvList(ii)%indices(j,3)
             alp = commSlidingCell_1st(ll,sps)%recvList(ii)%weight(j)

             jj = commSlidingCell_1st(ll,sps)%nRecvCum(ii-1) &
                + commSlidingCell_1st(ll,sps)%recvList(ii)%indRecv(j)

             ! Update f1 of the halo. Note that level == 1 and not
             ! currentLevel for dw.

             flowDoms(d2,1,sps)%dw(i2,j2,k2,if1SST) =     &
                 flowDoms(d2,1,sps)%dw(i2,j2,k2,if1SST) + &
                       alp*recvBuffer(jj)

           enddo

         enddo completeRecvs

         ! Complete the nonblocking sends.

         size = commSlidingCell_1st(ll,sps)%nProcSend
         do i=1,commSlidingCell_1st(ll,sps)%nProcSend
           call mpi_waitany(size, sendRequests, index, status, ierr)
         enddo

       enddo spectralModes

       end subroutine exchangeF1SSTSliding

!      ==================================================================

       subroutine exchangeF1SSTOverset
!
!      ******************************************************************
!      *                                                                *
!      * exchangeF1SSTOverset communicates the overset boundary values  *
!      * for the blending function f1 of the SST model. This variable   *
!      * is stored in dw(:,:,:,if1SST).                                 *
!      *                                                                *
!      ******************************************************************
!
       use block
       use communication
       use inputTimeSpectral
       use iteration
       implicit none
!
!      Local variables.
!
       integer :: size, procID, ierr, index
       integer, dimension(mpi_status_size) :: status

       integer(kind=intType) :: i, j, ii, jj, sps, ll
       integer(kind=intType) :: d1, i1, j1, k1, d2, i2, j2, k2

       real(kind=realType), dimension(:), pointer :: weight
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Easier storage of the current mg level.

       ll = currentLevel

       ! Loop over the number of spectral solutions.

       spectralModes: do sps=1,nTimeIntervalsSpectral

         ii = 1
         sends: do i=1,commPatternOverset(ll,sps)%nProcSend

           ! Store the processor id and the size of the message
           ! a bit easier.

           procID = commPatternOverset(ll,sps)%sendProc(i)
           size   = commPatternOverset(ll,sps)%nSend(i)

           ! Copy the data in the correct part of the send buffer.

           jj = ii
           do j=1,commPatternOverset(ll,sps)%nSend(i)

             ! Store the block id and the indices of the donor a
             !  bit easier.

             d1 = commPatternOverset(ll,sps)%sendList(i)%block(j)
             i1 = commPatternOverset(ll,sps)%sendList(i)%indices(j,1)
             j1 = commPatternOverset(ll,sps)%sendList(i)%indices(j,2)
             k1 = commPatternOverset(ll,sps)%sendList(i)%indices(j,3)

             weight => commPatternOverset(ll,sps)%sendList(i)%interp(j,:)

             ! Store the value of f1 in the send buffer. Note that the
             ! level is 1 and not ll (= currentLevel).

             sendBuffer(jj) = &
              weight(1)*flowDoms(d1,1,sps)%dw(i1  ,j1  ,k1  ,if1SST) + &
              weight(2)*flowDoms(d1,1,sps)%dw(i1+1,j1  ,k1  ,if1SST) + &
              weight(3)*flowDoms(d1,1,sps)%dw(i1  ,j1+1,k1  ,if1SST) + &
              weight(4)*flowDoms(d1,1,sps)%dw(i1+1,j1+1,k1  ,if1SST) + &
              weight(5)*flowDoms(d1,1,sps)%dw(i1  ,j1  ,k1+1,if1SST) + &
              weight(6)*flowDoms(d1,1,sps)%dw(i1+1,j1  ,k1+1,if1SST) + &
              weight(7)*flowDoms(d1,1,sps)%dw(i1  ,j1+1,k1+1,if1SST) + &
              weight(8)*flowDoms(d1,1,sps)%dw(i1+1,j1+1,k1+1,if1SST)
             jj = jj + 1

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
         receives: do i=1,commPatternOverset(ll,sps)%nProcRecv

           ! Store the processor id and the size of the message
           ! a bit easier.

           procID = commPatternOverset(ll,sps)%recvProc(i)
           size   = commPatternOverset(ll,sps)%nRecv(i)

           ! Post the receive.

           call mpi_irecv(recvBuffer(ii), size, sumb_real, procID, &
                          myID, SUmb_comm_world, recvRequests(i),  &
                          ierr)

           ! And update ii.

           ii = ii + size

         enddo receives

         ! Copy the local data.

         localCopy: do i=1,internalOverset(ll,sps)%ncopy

           ! Store the block and the indices of the donor a bit easier.

           d1 = internalOverset(ll,sps)%donorBlock(i)
           i1 = internalOverset(ll,sps)%donorIndices(i,1)
           j1 = internalOverset(ll,sps)%donorIndices(i,2)
           k1 = internalOverset(ll,sps)%donorIndices(i,3)

           weight => internalOverset(ll,sps)%donorInterp(i,:)

           ! Idem for the halo's.

           d2 = internalOverset(ll,sps)%haloBlock(i)
           i2 = internalOverset(ll,sps)%haloIndices(i,1)
           j2 = internalOverset(ll,sps)%haloIndices(i,2)
           k2 = internalOverset(ll,sps)%haloIndices(i,3)

           ! Copy the values. Note that level is 1 and not
           ! ll (= currentLevel).

           flowDoms(d2,1,sps)%dw(i2,j2,k2,if1SST) = &
              weight(1)*flowDoms(d1,1,sps)%dw(i1  ,j1  ,k1  ,if1SST) + &
              weight(2)*flowDoms(d1,1,sps)%dw(i1+1,j1  ,k1  ,if1SST) + &
              weight(3)*flowDoms(d1,1,sps)%dw(i1  ,j1+1,k1  ,if1SST) + &
              weight(4)*flowDoms(d1,1,sps)%dw(i1+1,j1+1,k1  ,if1SST) + &
              weight(5)*flowDoms(d1,1,sps)%dw(i1  ,j1  ,k1+1,if1SST) + &
              weight(6)*flowDoms(d1,1,sps)%dw(i1+1,j1  ,k1+1,if1SST) + &
              weight(7)*flowDoms(d1,1,sps)%dw(i1  ,j1+1,k1+1,if1SST) + &
              weight(8)*flowDoms(d1,1,sps)%dw(i1+1,j1+1,k1+1,if1SST)

         enddo localCopy

         ! Complete the nonblocking receives in an arbitrary sequence and
         ! copy the variables from the buffer into the halo's.

         size = commPatternOverset(ll,sps)%nProcRecv
         completeRecvs: do i=1,commPatternOverset(ll,sps)%nProcRecv

           ! Complete any of the requests.

           call mpi_waitany(size, recvRequests, index, status, ierr)

           ! Copy the data just arrived in the halo's.

           ii = index
           jj = commPatternOverset(ll,sps)%nRecvCum(ii-1) +1
           do j=1,commPatternOverset(ll,sps)%nRecv(ii)

             ! Store the block and the indices of the halo a bit easier.

             d2 = commPatternOverset(ll,sps)%recvList(ii)%block(j)
             i2 = commPatternOverset(ll,sps)%recvList(ii)%indices(j,1)
             j2 = commPatternOverset(ll,sps)%recvList(ii)%indices(j,2)
             k2 = commPatternOverset(ll,sps)%recvList(ii)%indices(j,3)

             ! And copy the data in the appropriate place in dw. Note
             ! that level == 1 and not ll (= currentLevel).

             flowDoms(d2,1,sps)%dw(i2,j2,k2,if1SST) = recvBuffer(jj)
             jj = jj + 1

           enddo

         enddo completeRecvs

         ! Complete the nonblocking sends.

         size = commPatternOverset(ll,sps)%nProcSend
         do i=1,commPatternOverset(ll,sps)%nProcSend
           call mpi_waitany(size, sendRequests, index, status, ierr)
         enddo

       enddo spectralModes

       end subroutine exchangeF1SSTOverset

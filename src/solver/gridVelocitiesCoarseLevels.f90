!
!      ******************************************************************
!      *                                                                *
!      * File:          gridVelocitiesCoarseLevels.f90                  *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 02-23-2004                                      *
!      * Last modified: 06-12-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine gridVelocitiesCoarseLevels(sps)
!
!      ******************************************************************
!      *                                                                *
!      * gridVelocitiesCoarseLevels computes the grid velocities for    *
!      * the cell centers and the normal grid velocities for the faces  *
!      * of moving blocks on the coarser grid levels. GroundLevel is    *
!      * considered the fine grid level.                                *
!      *                                                                *
!      ******************************************************************
!
       use blockPointers
       use iteration
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType), intent(in) :: sps
!
!      Local variables.
!
       integer(kind=intType) :: nLevels, level, levm1, nn, mm
       integer(kind=intType) :: i, j, k, iie, jje, kke
       integer(kind=intType) :: ii, ii1, jj, jj1, kk, kk1

       integer(kind=intType), dimension(:,:), pointer :: iFine, jFine
       integer(kind=intType), dimension(:,:), pointer :: kFine

       real(kind=realType) :: jjWeight, kkWeight, weight

       real(kind=realType), dimension(:), pointer :: jWeight, kWeight

       real(kind=realType), dimension(:,:),     pointer :: sFine, sFace
       real(kind=realType), dimension(:,:,:,:), pointer :: sf
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Loop over the number of coarse grid levels, starting at
       ! groundLevel+1,

       nLevels = ubound(flowDoms,2)
       levelLoop: do level=groundLevel+1,nLevels

         ! Loop over the number of local blocks.

         domains: do nn=1,nDom

           ! Set the pointers for this block.

           call setPointers(nn, level, sps)

           ! Check for a moving block. As it is possible that in a
           ! multidisicplinary environment additional grid velocities
           ! are set, the test should be done on addGridVelocities
           ! and not on blockIsMoving.

           testMoving: if( addGridVelocities ) then
!
!            ************************************************************
!            *                                                          *
!            * Grid velocities of the cell centers, including the 1st   *
!            * level halo cells. These are determined by accumulating   *
!            * the fine grid values. At the end the internal halo's are *
!            * communicated to obtain the correct values.               *
!            *                                                          *
!            ************************************************************
!
             levm1 = level - 1

             ! Set the pointer sf to the cell velocities on the fine mesh.

             sf => flowDoms(nn,levm1,sps)%s

             ! Loop over the cells, including the 1st level halo's.
             ! The indices kk, kk1 contain the corresponding fine
             ! grid indices in k-direction. Idem for jj, jj1, ii
             ! and ii1.

             do k=1,ke
               if(k == 1) then
                 kk = 1; kk1 = 1
               else if(k == ke) then
                 kk = flowDoms(nn,levm1,sps)%ke; kk1 = kk
               else
                 kk = mgKFine(k,1); kk1 = mgKFine(k,2)
               endif

               do j=1,je
                 if(j == 1) then
                   jj = 1; jj1 = 1
                 else if(j == je) then
                   jj = flowDoms(nn,levm1,sps)%je; jj1 = jj
                 else
                   jj = mgJFine(j,1); jj1 = mgJFine(j,2)
                 endif

                 do i=1,ie
                   if(i == 1) then
                     ii = 1; ii1 = 1
                   else if(i == ie) then
                     ii = flowDoms(nn,levm1,sps)%ie; ii1 = ii
                   else
                     ii = mgIFine(i,1); ii1 = mgIFine(i,2)
                   endif

                   ! Determine the coarse grid velocity by
                   ! averaging the fine grid values.

                   s(i,j,k,1) = (sf(ii1,jj1,kk1,1) + sf(ii,jj1,kk1,1)  &
                              +  sf(ii1,jj, kk1,1) + sf(ii,jj, kk1,1)  &
                              +  sf(ii1,jj1,kk, 1) + sf(ii,jj1,kk, 1)  &
                              +  sf(ii1,jj, kk, 1) + sf(ii,jj, kk, 1)) &
                              * eighth
                   s(i,j,k,2) = (sf(ii1,jj1,kk1,2) + sf(ii,jj1,kk1,2)  &
                              +  sf(ii1,jj, kk1,2) + sf(ii,jj, kk1,2)  &
                              +  sf(ii1,jj1,kk, 2) + sf(ii,jj1,kk, 2)  &
                              +  sf(ii1,jj, kk, 2) + sf(ii,jj, kk, 2)) &
                              * eighth
                   s(i,j,k,3) = (sf(ii1,jj1,kk1,3) + sf(ii,jj1,kk1,3)  &
                              +  sf(ii1,jj, kk1,3) + sf(ii,jj, kk1,3)  &
                              +  sf(ii1,jj1,kk, 3) + sf(ii,jj1,kk, 3)  &
                              +  sf(ii1,jj, kk, 3) + sf(ii,jj, kk, 3)) &
                              * eighth
                 enddo
               enddo
             enddo
!
!            ************************************************************
!            *                                                          *
!            * Normal grid velocities of the faces.                     *
!            *                                                          *
!            ************************************************************
!
             ! Loop over the three directions.

             loopCoarseDir: do mm=1,3

               ! Set some values depending on the situation.

               select case (mm)

                 case (1_intType)       ! Normals in i-direction
                   iie = ie; jje = je; kke = ke
                   iFine => mgIFine; jFine => mgJFine; kFine => mgKFine
                   jWeight => mgJWeight; kWeight => mgKWeight

                 case (2_intType)       ! Normals in j-direction
                   iie = je; jje = ie; kke = ke
                   iFine => mgJFine; jFine => mgIFine; kFine => mgKFine
                   jWeight => mgIWeight; kWeight => mgKWeight

                 case (3_intType)       ! Normals in k-direction
                   iie = ke; jje = ie; kke = je
                   iFine => mgKFine; jFine => mgIFine; kFine => mgJFine
                   jWeight => mgIWeight; kWeight => mgJWeight

               end select
!
!              **********************************************************
!              *                                                        *
!              * Normal grid velocities in generalized i-direction.     *
!              * mm == 1: i-direction                                   *
!              * mm == 2: j-direction                                   *
!              * mm == 3: k-direction                                   *
!              *                                                        *
!              **********************************************************
!
               do i=0,iie

                 ! Determine the i-index of the corresponding plane on
                 ! the fine grid. Note that halo planes are not entirely
                 ! correct. This is not really a problem.

                 if(i < 2) then
                   ii = i
                 else if(i < iie) then
                   ii = iFine(i,2)
                 else
                   ii = iFine(iie-1,2) + 1
                 endif

                 ! Set the pointers for sFine and sFace, which will
                 ! contain the mesh velocities for this particular
                 ! plane. The index depends on the value of mm.

                 select case (mm)
                   case (1_intType)
                     sFine => flowDoms(nn,levm1,sps)%sFaceI(ii,:,:)
                     sFace => sFaceI(i,:,:)
                   case (2_intType)
                     sFine => flowDoms(nn,levm1,sps)%sFaceJ(:,ii,:)
                     sFace => sFaceJ(:,i,:)
                   case (3_intType)
                     sFine => flowDoms(nn,levm1,sps)%sFaceK(:,:,ii)
                     sFace => sFaceK(:,:,i)
                 end select

                 ! Loop over the k and j faces for this general i-plane.
                 ! Again the halo's are not entirely correct. Kk, kk1,
                 ! jj and jj1 are the children in k and j-direction
                 ! respectively.

                 do k=1,kke
                   if(k == 1) then
                     kk = 1; kk1 = 1
                     kkWeight = kWeight(2)
                   else if(k == kke) then
                     kk = kFine(kke-1,2) + 1; kk1 = kk
                     kkWeight = kWeight(kke-1)
                   else
                     kk = kFine(k,1); kk1 = kFine(k,2)
                     kWeight = kWeight(k)
                   endif

                   do j=1,jje
                     if(j == 1) then
                       jj = 1; jj1 = 1
                       jjWeight = jWeight(2)
                     else if(j == jje) then
                       jj = jFine(jje-1,2) + 1; jj1 = jj
                       jWeight = jWeight(jje-1)
                     else
                       jj = jFine(j,1); jj1 = jFine(j,2)
                       jjWeight = jWeight(j)
                     endif

                     ! Determine the coarse grid normal velocity.
                     ! Take the averaging weight into account; for
                     ! a normal coarsening this weight is 1.0.

                     weight = kkWeight*jjWeight
                     sFace(j,k) = weight*(sFine(jj1,kk1) &
                                +         sFine(jj ,kk1) &
                                +         sFine(jj1,kk)  &
                                +         sFine(jj ,kk))
                   enddo
                 enddo
               enddo

             enddo loopCoarseDir
           endif testMoving
         enddo domains

         ! Exchange the cell centered velocities.

         call exchangeCellGridVelocities(level,sps)

       enddo levelLoop

       end subroutine gridVelocitiesCoarseLevels

!      ==================================================================

       subroutine exchangeCellGridVelocities(level,sps)
!
!      ******************************************************************
!      *                                                                *
!      * exchangeCellGridVelocities exchanges the grid velocities in    *
!      * the cell centers for the given grid level and spectral         *
!      * solution.                                                      *
!      *                                                                *
!      ******************************************************************
!
       use block
       use commSliding
       use communication
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType), intent(in) :: level, sps
!
!      Local variables.
!
       integer :: size, procID, ierr, index
       integer, dimension(mpi_status_size) :: status

       integer(kind=intType) :: i, j, ii, jj
       integer(kind=intType) :: d1, i1, j1, k1, d2, i2, j2, k2

       real(kind=realType) :: alp
       real(kind=realType), dimension(3) :: vv
!
!      Interfaces
!
       interface
         subroutine correctPeriodicGridVel(level, sps, nPeriodic, &
                                           periodicData)
         use block
         use communication
         implicit none

         integer(kind=intType), intent(in) :: level, sps, nPeriodic
         type(periodicDataType), dimension(:), pointer :: periodicData

         end subroutine correctPeriodicGridVel
       end interface
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
!      ******************************************************************
!      *                                                                *
!      * The 1 to 1 communication.                                      *
!      *                                                                *
!      ******************************************************************
!
       ! Send the variables. The data is first copied into
       ! the send buffer after which the buffer is sent asap.

       ii = 1
       sends: do i=1,commPatternCell_1st(level)%nProcSend

         ! Store the processor id and the size of the message
         ! a bit easier.

         procID = commPatternCell_1st(level)%sendProc(i)
         size   = 3*commPatternCell_1st(level)%nSend(i)

         ! Copy the data in the correct part of the send buffer.

         jj = ii
         do j=1,commPatternCell_1st(level)%nSend(i)

           ! Store the block id and the indices of the donor
           ! a bit easier.

           d1 = commPatternCell_1st(level)%sendList(i)%block(j)
           i1 = commPatternCell_1st(level)%sendList(i)%indices(j,1)
           j1 = commPatternCell_1st(level)%sendList(i)%indices(j,2)
           k1 = commPatternCell_1st(level)%sendList(i)%indices(j,3)

           ! Store the grid velocities in sendBuffer, if they
           ! are allocated. Otherwise they are simply zero.

           if( flowDoms(d1,level,sps)%addGridVelocities ) then
             sendBuffer(jj)   = flowDoms(d1,level,sps)%s(i1,j1,k1,1)
             sendBuffer(jj+1) = flowDoms(d1,level,sps)%s(i1,j1,k1,2)
             sendBuffer(jj+2) = flowDoms(d1,level,sps)%s(i1,j1,k1,3)
           else
             sendBuffer(jj)   = zero
             sendBuffer(jj+1) = zero
             sendBuffer(jj+2) = zero
           endif

           jj = jj + 3
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
         size   = 3*commPatternCell_1st(level)%nRecv(i)

         ! Post the receive.

         call mpi_irecv(recvBuffer(ii), size, sumb_real, procID, &
                        myID, SUmb_comm_world, recvRequests(i), ierr)

         ! And update ii.

         ii = ii + size
       enddo receives

       ! Copy the local data.

       localCopy: do i=1,internalCell_1st(level)%nCopy

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

         ! Copy the grid velocities, if they are both allocated.
         ! Otherwise they are either set to zero or nothing is done.

         if( flowDoms(d2,level,sps)%addGridVelocities ) then
           if( flowDoms(d1,level,sps)%addGridVelocities ) then
             flowDoms(d2,level,sps)%s(i2,j2,k2,1) = &
                flowDoms(d1,level,sps)%s(i1,j1,k1,1)
             flowDoms(d2,level,sps)%s(i2,j2,k2,2) = &
                flowDoms(d1,level,sps)%s(i1,j1,k1,2)
             flowDoms(d2,level,sps)%s(i2,j2,k2,3) = &
                flowDoms(d1,level,sps)%s(i1,j1,k1,3)
           else
             flowDoms(d2,level,sps)%s(i2,j2,k2,1) = zero
             flowDoms(d2,level,sps)%s(i2,j2,k2,2) = zero
             flowDoms(d2,level,sps)%s(i2,j2,k2,3) = zero
           endif
         endif

       enddo localCopy

       ! Correct the periodic halo's of the internal communication
       ! pattern, if present.

       if(internalCell_1st(level)%nPeriodic > 0) &
         call correctPeriodicGridVel(level, sps,                        &
                                     internalCell_1st(level)%nPeriodic, &
                                     internalCell_1st(level)%periodicData)

       ! Complete the nonblocking receives in an arbitrary sequence and
       ! copy the variables from the buffer into the halo's.

       size = commPatternCell_1st(level)%nProcRecv
       completeRecvs: do i=1,commPatternCell_1st(level)%nProcRecv

         ! Complete any of the requests.

         call mpi_waitany(size, recvRequests, index, status, ierr)

         ! Copy the data just arrived in the halo's.

         ii = index
         jj = 3*commPatternCell_1st(level)%nRecvCum(ii-1)
         do j=1,commPatternCell_1st(level)%nRecv(ii)

           ! Store the block and the indices of the halo a bit easier.

           d2 = commPatternCell_1st(level)%recvList(ii)%block(j)
           i2 = commPatternCell_1st(level)%recvList(ii)%indices(j,1)
           j2 = commPatternCell_1st(level)%recvList(ii)%indices(j,2)
           k2 = commPatternCell_1st(level)%recvList(ii)%indices(j,3)

           ! Copy the grid velocities from recvBuffer if they are
           ! both allocated.

           if( flowDoms(d2,level,sps)%addGridVelocities ) then
             flowDoms(d2,level,sps)%s(i2,j2,k2,1) = recvBuffer(jj+1)
             flowDoms(d2,level,sps)%s(i2,j2,k2,2) = recvBuffer(jj+2)
             flowDoms(d2,level,sps)%s(i2,j2,k2,3) = recvBuffer(jj+3)
           endif

           jj = jj + 3
         enddo

       enddo completeRecvs

       ! Correct the periodic halo's of the external communication
       ! pattern, if present.

       if(commPatternCell_1st(level)%nPeriodic > 0) &
         call correctPeriodicGridVel(level, sps,                           &
                                     commPatternCell_1st(level)%nPeriodic, &
                                     commPatternCell_1st(level)%periodicData)

       ! Complete the nonblocking sends.

       size = commPatternCell_1st(level)%nProcSend
       do i=1,commPatternCell_1st(level)%nProcSend
         call mpi_waitany(size, sendRequests, index, status, ierr)
       enddo
!
!      ******************************************************************
!      *                                                                *
!      * The sliding mesh communication.                                *
!      *                                                                *
!      ******************************************************************
!
       ! Send the variables. The data is first copied into
       ! the send buffer after which the buffer is sent asap.

       ii = 1
       slidingSends: do i=1,commSlidingCell_1st(level,sps)%nProcSend

         ! Store the processor id and the size of the message
         ! a bit easier.

         procID = commSlidingCell_1st(level,sps)%sendProc(i)
         size   = 3*commSlidingCell_1st(level,sps)%nSend(i)

         ! Copy the data in the correct part of the send buffer.

         jj = ii
         do j=1,commSlidingCell_1st(level,sps)%nSend(i)

           ! Store the block id and the indices of the donor
           ! a bit easier.

           d1 = commSlidingCell_1st(level,sps)%sendList(i)%block(j)
           i1 = commSlidingCell_1st(level,sps)%sendList(i)%indices(j,1)
           j1 = commSlidingCell_1st(level,sps)%sendList(i)%indices(j,2)
           k1 = commSlidingCell_1st(level,sps)%sendList(i)%indices(j,3)

           ! Copy the velocities in the send buffer, if allocated.
           ! Otherwise they are set to zero.

           if( flowDoms(d1,level,sps)%addGridVelocities ) then
             sendBuffer(jj)   = flowDoms(d1,level,sps)%s(i1,j1,k1,1)
             sendBuffer(jj+1) = flowDoms(d1,level,sps)%s(i1,j1,k1,2)
             sendBuffer(jj+2) = flowDoms(d1,level,sps)%s(i1,j1,k1,3)
           else
             sendBuffer(jj)   = zero
             sendBuffer(jj+1) = zero
             sendBuffer(jj+2) = zero
           endif

           jj = jj + 3

         enddo

         ! Send the data.

         call mpi_isend(sendBuffer(ii), size, sumb_real, procID,  &
                        procID, SUmb_comm_world, sendRequests(i), &
                        ierr)

         ! Set ii to jj for the next processor.

         ii = jj

       enddo slidingSends

       ! Post the nonblocking receives.

       ii = 1
       slidingReceives: do i=1,commSlidingCell_1st(level,sps)%nProcRecv

         ! Store the processor id and the size of the message
         ! a bit easier.

         procID = commSlidingCell_1st(level,sps)%recvProc(i)
         size    = 3*commSlidingCell_1st(level,sps)%nRecv(i)

         ! Post the receive.

         call mpi_irecv(recvBuffer(ii), size, sumb_real, procID, &
                        myID, SUmb_comm_world, recvRequests(i), ierr)

         ! And update ii.

         ii = ii + size

       enddo slidingReceives

       ! Initialize the sliding mesh halo's to zero.

       initHalos: do i=1,intSlidingCell_1st(level,sps)%nslidingHalos

         ! Store the block and the indices of the halo a bit easier.

         d2 = intSlidingCell_1st(level,sps)%slidingHaloList%block(i)
         i2 = intSlidingCell_1st(level,sps)%slidingHaloList%indices(i,1)
         j2 = intSlidingCell_1st(level,sps)%slidingHaloList%indices(i,2)
         k2 = intSlidingCell_1st(level,sps)%slidingHaloList%indices(i,3)

         if( flowDoms(d2,level,sps)%addGridVelocities ) then
           flowDoms(d2,level,sps)%s(i2,j2,k2,1) = zero
           flowDoms(d2,level,sps)%s(i2,j2,k2,2) = zero
           flowDoms(d2,level,sps)%s(i2,j2,k2,3) = zero
         endif

       enddo initHalos

       ! Interpolation with the locally stored data.

       localData: do i=1,intSlidingCell_1st(level,sps)%nCopy

         ! Store the block and the indices of the donor a bit easier.

         d1 = intSlidingCell_1st(level,sps)%donorList%block(i)
         i1 = intSlidingCell_1st(level,sps)%donorList%indices(i,1)
         j1 = intSlidingCell_1st(level,sps)%donorList%indices(i,2)
         k1 = intSlidingCell_1st(level,sps)%donorList%indices(i,3)

         ! Store the block and the indices of the halo a bit easier.

         d2 = intSlidingCell_1st(level,sps)%haloList%block(i)
         i2 = intSlidingCell_1st(level,sps)%haloList%indices(i,1)
         j2 = intSlidingCell_1st(level,sps)%haloList%indices(i,2)
         k2 = intSlidingCell_1st(level,sps)%haloList%indices(i,3)

         ! Store the weight a bit easier.

         alp = intSlidingCell_1st(level,sps)%weight(i)

         ! Update the grid velocities. Only if both have been allocated.

         if(flowDoms(d1,level,sps)%addGridVelocities .and. &
            flowDoms(d2,level,sps)%addGridVelocities ) then

           flowDoms(d2,level,sps)%s(i2,j2,k2,1:3) = &
               flowDoms(d2,level,sps)%s(i2,j2,k2,1:3) + &
                               alp*flowDoms(d1,level,sps)%s(i1,j1,k1,1:3)
         endif

       enddo localData

       ! Complete the nonblocking receives in an arbitrary sequence and
       ! use the variables from the buffer to interpolate the halo's.

       size = commSlidingCell_1st(level,sps)%nProcRecv
       slidingComplete: do i=1,commSlidingCell_1st(level,sps)%nProcRecv

         ! Complete any of the requests.

         call mpi_waitany(size, recvRequests, index, status, ierr)

         ! Update the halo's using the data just arrived.

         ii = index
         do j=1,commSlidingCell_1st(level,sps)%recvList(ii)%nCopy

           ! Store the block and the indices of the halo as well as
           ! the interpolation weight and the starting index in the
           ! receive buffer a bit easier.

           d2  = commSlidingCell_1st(level,sps)%recvList(ii)%block(j)
           i2  = commSlidingCell_1st(level,sps)%recvList(ii)%indices(j,1)
           j2  = commSlidingCell_1st(level,sps)%recvList(ii)%indices(j,2)
           k2  = commSlidingCell_1st(level,sps)%recvList(ii)%indices(j,3)
           alp = commSlidingCell_1st(level,sps)%recvList(ii)%weight(j)

           jj = 3*(commSlidingCell_1st(level,sps)%nRecvCum(ii-1)          &
              +    commSlidingCell_1st(level,sps)%recvList(ii)%indRecv(j) &
              -    1)

           ! Update the grid velocities in the halo cells, if allocated.

           if( flowDoms(d2,level,sps)%addGridVelocities ) then
             flowDoms(d2,level,sps)%s(i2,j2,k2,1) = &
                flowDoms(d2,level,sps)%s(i2,j2,k2,1) + alp*recvBuffer(jj+1)
             flowDoms(d2,level,sps)%s(i2,j2,k2,1) = &
                flowDoms(d2,level,sps)%s(i2,j2,k2,2) + alp*recvBuffer(jj+2)
             flowDoms(d2,level,sps)%s(i2,j2,k2,1) = &
                flowDoms(d2,level,sps)%s(i2,j2,k2,3) + alp*recvBuffer(jj+3)
           endif

           jj = jj + 3

         enddo

       enddo slidingComplete

       ! Apply the transformation matrix to the velocities.

       do i=1,intSlidingCell_1st(level,sps)%nslidingHalos

         ! Store the block and the indices of the halo as well as
         ! the index of the rotation matrix a bit easier.

         d2 = intSlidingCell_1st(level,sps)%slidingHaloList%block(i)
         i2 = intSlidingCell_1st(level,sps)%slidingHaloList%indices(i,1)
         j2 = intSlidingCell_1st(level,sps)%slidingHaloList%indices(i,2)
         k2 = intSlidingCell_1st(level,sps)%slidingHaloList%indices(i,3)
         j  = intSlidingCell_1st(level,sps)%rotIndex(i)

         ! Only correct the velocity if j > 0; j == 0 indicates
         ! that no transformation needs to be applied.

         if(j > 0 .and. flowDoms(d2,level,sps)%addGridVelocities ) then

           ! Apply the correct rotation to the velocity.

           vv(1) = flowDoms(d2,level,sps)%s(i2,j2,k2,1)
           vv(2) = flowDoms(d2,level,sps)%s(i2,j2,k2,2)
           vv(3) = flowDoms(d2,level,sps)%s(i2,j2,k2,3)

           flowDoms(d2,level,sps)%s(i2,j2,k2,1) = rotSliding(j,1,1)*vv(1) &
                                                + rotSliding(j,1,2)*vv(2) &
                                                + rotSliding(j,1,3)*vv(3)

           flowDoms(d2,level,sps)%s(i2,j2,k2,2) = rotSliding(j,2,1)*vv(1) &
                                                + rotSliding(j,2,2)*vv(2) &
                                                + rotSliding(j,2,3)*vv(3)

           flowDoms(d2,level,sps)%s(i2,j2,k2,3) = rotSliding(j,3,1)*vv(1) &
                                                + rotSliding(j,3,2)*vv(2) &
                                                + rotSliding(j,3,3)*vv(3)
         endif
       enddo

       ! Complete the nonblocking sends.

       size = commSlidingCell_1st(level,sps)%nProcSend
       do i=1,commSlidingCell_1st(level,sps)%nProcSend
         call mpi_waitany(size, sendRequests, index, status, ierr)
       enddo

       end subroutine exchangeCellGridVelocities

!      ==================================================================

       subroutine correctPeriodicGridVel(level, sps, nPeriodic, &
                                         periodicData)
!
!      ******************************************************************
!      *                                                                *
!      * correctPeriodicGridVel applies the periodic transformation     *
!      * to the grid velocities of the cell halo's in periodicData.     *
!      *                                                                *
!      ******************************************************************
!
       use block
       use communication
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType), intent(in) :: level, sps, nPeriodic
       type(periodicDataType), dimension(:), pointer :: periodicData
!
!      Local variables.
!
       integer(kind=intType) :: nn, mm, ii, i, j, k
       real(kind=realType)   :: vx, vy, vz

       real(kind=realType), dimension(3,3) :: rotMatrix
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Loop over the number of periodic transformations.

       do nn=1,nPeriodic

         ! Store the rotation matrix a bit easier.

         rotMatrix = periodicData(nn)%rotMatrix

         ! Loop over the number of halo cells for this transformation.

         do ii=1,periodicData(nn)%nHalos

           ! Store the block and the indices a bit easier.

           mm = periodicData(nn)%block(ii)
           i  = periodicData(nn)%indices(ii,1)
           j  = periodicData(nn)%indices(ii,2)
           k  = periodicData(nn)%indices(ii,3)

           ! Check if the grid velocities have been allocated.

           if( flowDoms(mm,level,sps)%addGridVelocities ) then

             ! Store the original velocities in vx, vy, vz.

             vx = flowDoms(mm,level,sps)%s(i,j,k,1)
             vy = flowDoms(mm,level,sps)%s(i,j,k,2)
             vz = flowDoms(mm,level,sps)%s(i,j,k,3)

             ! Compute the new velocity vector.

             flowDoms(mm,level,sps)%s(i,j,k,1) = rotMatrix(1,1)*vx &
                                               + rotMatrix(1,2)*vy &
                                               + rotMatrix(1,3)*vz
             flowDoms(mm,level,sps)%s(i,j,k,2) = rotMatrix(2,1)*vx &
                                               + rotMatrix(2,2)*vy &
                                               + rotMatrix(2,3)*vz
             flowDoms(mm,level,sps)%s(i,j,k,3) = rotMatrix(3,1)*vx &
                                               + rotMatrix(3,2)*vy &
                                               + rotMatrix(3,3)*vz

           endif
         enddo

       enddo

       end subroutine correctPeriodicGridVel

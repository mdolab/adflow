!
!      ******************************************************************
!      *                                                                *
!      * File:          whaloSliding.f90                                *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 09-25-2003                                      *
!      * Last modified: 09-16-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine whaloSliding(level, start, end, commPressure,       &
                               commVarGamma, commLamVis, commEddyVis, &
                               commPattern, internal, nlev)
!
!      ******************************************************************
!      *                                                                *
!      * whaloSliding exchanges the sliding mesh halo's for the cell    *
!      * centered variables for the given communication pattern. It     *
!      * is possible to send a range of variables and not the entire    *
!      * set, e.g. only the flow variables or only the turbulent        *
!      * variables. This is controlled by the arguments start, end,     *
!      * commPressure and commViscous. The exchange takes place for     *
!      * the given grid level.                                          *
!      *                                                                *
!      ******************************************************************
!
       use block
       use commSliding
       use communication
       use inputTimeSpectral
       implicit none
!
!      Subroutine arguments
!
       integer(kind=intType), intent(in) :: level, start, end, nlev
       logical, intent(in) :: commPressure, commVarGamma
       logical, intent(in) :: commLamVis, commEddyVis

       type(slidingCommType), dimension(nlev,*), &
                                intent(in) :: commPattern
       type(internalSlidingCommType), dimension(nlev,*), &
                                         intent(in) :: internal
!
!      Local variables.
!
       integer :: size, procID, ierr, index
       integer, dimension(mpi_status_size) :: status

       integer(kind=intType) :: nVar
       integer(kind=intType) :: i, j, k, ii, jj, mm
       integer(kind=intType) :: d1, i1, j1, k1, d2, i2, j2, k2

       real(kind=realType) :: alp
       real(kind=realType), dimension(3) :: vv

       logical :: correctVelocities
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Set the logical correctVelocities. Only if a velocity variable
       ! is communicated it is needed to apply the periodic
       ! transformations.

       correctVelocities = .false.
       if(start <= ivx .and. end >= ivz) correctVelocities = .true.

       ! Determine the number of variables per cell to be sent.

       nVar = max(0_intType,(end - start + 1))
       if( commPressure ) nVar = nVar + 1
       if( commVarGamma ) nVar = nVar + 1
       if( commLamVis )   nVar = nVar + 1
       if( commEddyVis )  nVar = nVar + 1

       if(nVar == 0) return

       ! Loop over the number of spectral solutions.

       spectralModes: do mm=1,nTimeIntervalsSpectral

         ! Send the variables. The data is first copied into
         ! the send buffer after which the buffer is sent asap.

         ii = 1
         sends: do i=1,commPattern(level,mm)%nProcSend

           ! Store the processor id and the size of the message
           ! a bit easier.

           procID = commPattern(level,mm)%sendProc(i)
           size    = nVar*commPattern(level,mm)%nsend(i)

           ! Copy the data in the correct part of the send buffer.

           jj = ii
           do j=1,commPattern(level,mm)%nsend(i)

             ! Store the block id and the indices of the donor
             ! a bit easier.

             d1 = commPattern(level,mm)%sendList(i)%block(j)
             i1 = commPattern(level,mm)%sendList(i)%indices(j,1)
             j1 = commPattern(level,mm)%sendList(i)%indices(j,2)
             k1 = commPattern(level,mm)%sendList(i)%indices(j,3)

             ! Copy the given range of the working variables for this
             ! cell in the buffer. Update the counter jj accordingly.

             do k=start,end
               sendBuffer(jj) = flowDoms(d1,level,mm)%w(i1,j1,k1,k)
               jj = jj + 1
             enddo

             ! The pressure, if needed.

             if( commPressure ) then
               sendBuffer(jj) = flowDoms(d1,level,mm)%p(i1,j1,k1)
               jj = jj + 1
             endif

             ! The specific heat ratio, if needed. Note that level == 1.

             if( commVarGamma ) then
               sendBuffer(jj) = flowDoms(d1,1,mm)%gamma(i1,j1,k1)
               jj = jj + 1
             endif

             ! The laminar viscosity for a viscous computation.
             ! Again level == 1.

             if( commLamVis ) then
               sendBuffer(jj) = flowDoms(d1,1,mm)%rlv(i1,j1,k1)
               jj = jj + 1
             endif

             ! The eddy viscosity ratio for eddy viscosity models.
             ! Level is the true multigrid level, because the eddy
             ! viscosity is allocated on all grid levels.

             if( commEddyVis ) then
               sendBuffer(jj) = flowDoms(d1,level,mm)%rev(i1,j1,k1)
               jj = jj + 1
             endif

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
         receives: do i=1,commPattern(level,mm)%nProcRecv

           ! Store the processor id and the size of the message
           ! a bit easier.

           procID = commPattern(level,mm)%recvProc(i)
           size   = nVar*commPattern(level,mm)%nrecv(i)

           ! Post the receive.

           call mpi_irecv(recvBuffer(ii), size, sumb_real, procID, &
                          myID, SUmb_comm_world, recvRequests(i), ierr)

           ! And update ii.

           ii = ii + size

         enddo receives

         ! Initialize the sliding mesh halo's to zero.

         initHalos: do i=1,internal(level,mm)%nslidingHalos

           ! Store the block and the indices of the halo a bit easier.

           d2 = internal(level,mm)%slidingHaloList%block(i)
           i2 = internal(level,mm)%slidingHaloList%indices(i,1)
           j2 = internal(level,mm)%slidingHaloList%indices(i,2)
           k2 = internal(level,mm)%slidingHaloList%indices(i,3)

           ! Initialize the variables to be communicated. Note that for
           ! gamma and rlv not level but 1 must be used.

           do k=start,end
             flowDoms(d2,level,mm)%w(i2,j2,k2,k) = zero
           enddo

           if( commPressure ) flowDoms(d2,level,mm)%p(i2,j2,k2)   = zero
           if( commVarGamma ) flowDoms(d2,1,mm)%gamma(i2,j2,k2)   = zero
           if( commLamVis   ) flowDoms(d2,1,mm)%rlv(i2,j2,k2)     = zero
           if( commEddyVis  ) flowDoms(d2,level,mm)%rev(i2,j2,k2) = zero

         enddo initHalos

         ! Interpolation with the locally stored data.

         localData: do i=1,internal(level,mm)%ncopy

           ! Store the block and the indices of the donor a bit easier.

           d1 = internal(level,mm)%donorList%block(i)
           i1 = internal(level,mm)%donorList%indices(i,1)
           j1 = internal(level,mm)%donorList%indices(i,2)
           k1 = internal(level,mm)%donorList%indices(i,3)

           ! Store the block and the indices of the halo a bit easier.

           d2 = internal(level,mm)%haloList%block(i)
           i2 = internal(level,mm)%haloList%indices(i,1)
           j2 = internal(level,mm)%haloList%indices(i,2)
           k2 = internal(level,mm)%haloList%indices(i,3)

           ! Store the weight a bit easier.

           alp = internal(level,mm)%weight(i)

           ! Update the halo variables. Note that for gamma and rlv
           ! not level but 1 must be used.

           do k=start,end
             flowDoms(d2,level,mm)%w(i2,j2,k2,k) = &
             flowDoms(d2,level,mm)%w(i2,j2,k2,k) + &
                                 alp*flowDoms(d1,level,mm)%w(i1,j1,k1,k)
           enddo

           if( commPressure )                   &
             flowDoms(d2,level,mm)%p(i2,j2,k2) = &
             flowDoms(d2,level,mm)%p(i2,j2,k2) + &
                                 alp*flowDoms(d1,level,mm)%p(i1,j1,k1)

           if( commVarGamma )                   &
             flowDoms(d2,1,mm)%gamma(i2,j2,k2) = &
             flowDoms(d2,1,mm)%gamma(i2,j2,k2) + &
                             alp*flowDoms(d1,1,mm)%gamma(i1,j1,k1)

           if( commLamVis )                   &
             flowDoms(d2,1,mm)%rlv(i2,j2,k2) = &
             flowDoms(d2,1,mm)%rlv(i2,j2,k2) + &
                             alp*flowDoms(d1,1,mm)%rlv(i1,j1,k1)

           if( commEddyVis )                      &
             flowDoms(d2,level,mm)%rev(i2,j2,k2) = &
             flowDoms(d2,level,mm)%rev(i2,j2,k2) + &
                                 alp*flowDoms(d1,level,mm)%rev(i1,j1,k1)

         enddo localData

         ! Complete the nonblocking receives in an arbitrary sequence and
         ! use the variables from the buffer to interpolate the halo's.

         size = commPattern(level,mm)%nProcRecv
         completeRecvs: do i=1,commPattern(level,mm)%nProcRecv

           ! Complete any of the requests.

           call mpi_waitany(size, recvRequests, index, status, ierr)

           ! Update the halo's using the data just arrived.

           ii = index

           do j=1,commPattern(level,mm)%recvList(ii)%ncopy

             ! Store the block and the indices of the halo as well as
             ! the interpolation weight and the starting index in the
             ! receive buffer a bit easier.

             d2  = commPattern(level,mm)%recvList(ii)%block(j)
             i2  = commPattern(level,mm)%recvList(ii)%indices(j,1)
             j2  = commPattern(level,mm)%recvList(ii)%indices(j,2)
             k2  = commPattern(level,mm)%recvList(ii)%indices(j,3)
             alp = commPattern(level,mm)%recvList(ii)%weight(j)

             jj = nVar*(commPattern(level,mm)%nrecvCum(ii-1)          &
                +       commPattern(level,mm)%recvList(ii)%indRecv(j) &
                -       1)

             ! Update the halo variables. Note that for gamma and rlv
             ! not level but 1 must be used.

             do k=start,end
               jj = jj + 1
               flowDoms(d2,level,mm)%w(i2,j2,k2,k) = &
               flowDoms(d2,level,mm)%w(i2,j2,k2,k) + alp*recvBuffer(jj)
             enddo

             if( commPressure ) then
               jj = jj + 1
               flowDoms(d2,level,mm)%p(i2,j2,k2) = &
               flowDoms(d2,level,mm)%p(i2,j2,k2) + alp*recvBuffer(jj)
             endif

             if( commVarGamma ) then
               jj = jj + 1
               flowDoms(d2,1,mm)%gamma(i2,j2,k2) = &
               flowDoms(d2,1,mm)%gamma(i2,j2,k2) + alp*recvBuffer(jj)
             endif

             if( commLamVis ) then
               jj = jj + 1
               flowDoms(d2,1,mm)%rlv(i2,j2,k2) = &
               flowDoms(d2,1,mm)%rlv(i2,j2,k2) + alp*recvBuffer(jj)
             endif

             if( commEddyVis ) then
               jj = jj + 1
               flowDoms(d2,level,mm)%rev(i2,j2,k2) = &
               flowDoms(d2,level,mm)%rev(i2,j2,k2) + alp*recvBuffer(jj)
             endif

           enddo

         enddo completeRecvs

         ! Apply the transformation matrix to the velocities if the
         ! velocities have been constructed.

         testCorr: if( correctVelocities ) then

           do i=1,internal(level,mm)%nslidingHalos

             ! Store the block and the indices of the halo as well as
             ! the index of the rotation matrix a bit easier.

             d2 = internal(level,mm)%slidingHaloList%block(i)
             i2 = internal(level,mm)%slidingHaloList%indices(i,1)
             j2 = internal(level,mm)%slidingHaloList%indices(i,2)
             k2 = internal(level,mm)%slidingHaloList%indices(i,3)
             j  = internal(level,mm)%rotIndex(i)

             ! Only correct the velocity if j > 0; j == 0 indicates
             ! that no transformation needs to be applied.

             if(j > 0) then

               ! Apply the correct rotation to the velocity.

               vv(1) = flowDoms(d2,level,mm)%w(i2,j2,k2,ivx)
               vv(2) = flowDoms(d2,level,mm)%w(i2,j2,k2,ivy)
               vv(3) = flowDoms(d2,level,mm)%w(i2,j2,k2,ivz)

               flowDoms(d2,level,mm)%w(i2,j2,k2,ivx) = &
                                     rotSliding(j,1,1)*vv(1) + &
                                     rotSliding(j,1,2)*vv(2) + &
                                     rotSliding(j,1,3)*vv(3)

               flowDoms(d2,level,mm)%w(i2,j2,k2,ivy) = &
                                     rotSliding(j,2,1)*vv(1) + &
                                     rotSliding(j,2,2)*vv(2) + &
                                     rotSliding(j,2,3)*vv(3)

               flowDoms(d2,level,mm)%w(i2,j2,k2,ivz) = &
                                     rotSliding(j,3,1)*vv(1) + &
                                     rotSliding(j,3,2)*vv(2) + &
                                     rotSliding(j,3,3)*vv(3)
             endif
           enddo
         endif testCorr

         ! Complete the nonblocking sends.

         size = commPattern(level,mm)%nProcSend
         do i=1,commPattern(level,mm)%nProcSend
           call mpi_waitany(size, sendRequests, index, status, ierr)
         enddo

       enddo spectralModes

       end subroutine whaloSliding

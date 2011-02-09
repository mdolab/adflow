!
!      ******************************************************************
!      *                                                                *
!      * File:          whalo1to1.f90                                   *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 03-07-2003                                      *
!      * Last modified: 09-16-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine whalo1to1(level, start, end, commPressure,       &
                            commVarGamma, commLamVis, commEddyVis, &
                            commPattern, internal)
!
!      ******************************************************************
!      *                                                                *
!      * whalo1to1 exchanges the 1 to 1 internal halo's for the cell    *
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
       use communication
       use inputTimeSpectral
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType), intent(in) :: level, start, end
       logical, intent(in) :: commPressure, commVarGamma
       logical, intent(in) :: commLamVis, commEddyVis

       type(commType), dimension(*), intent(in)         :: commPattern
       type(internalCommType), dimension(*), intent(in) :: internal
!
!      Local variables.
!
       integer :: size, procID, ierr, index
       integer, dimension(mpi_status_size) :: status

       integer(kind=intType) :: nVar, mm
       integer(kind=intType) :: i, j, k, ii, jj
       integer(kind=intType) :: d1, i1, j1, k1, d2, i2, j2, k2

       logical :: correctPeriodic
!
!      Interfaces
!
       interface
         subroutine correctPeriodicVelocity(level, sp, nPeriodic, &
                                            periodicData)
         use block
         use communication
         use constants
         implicit none

         integer(kind=intType), intent(in) :: level, sp, nPeriodic
         type(periodicDataType), dimension(:), pointer :: periodicData

         end subroutine correctPeriodicVelocity
       end interface
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Set the logical correctPeriodic. Only if a momentum variable
       ! is communicated it is needed to apply the periodic
       ! transformations.

       correctPeriodic = .false.
       if(start <= ivx .and. end >= ivz) correctPeriodic = .true.

       ! Determine the number of variables per cell to be sent.

       nVar = max(0_intType,(end - start + 1))
       if( commPressure )  nVar = nVar + 1
       if( commVarGamma ) nVar = nVar + 1
       if( commLamVis )   nVar = nVar + 1
       if( commEddyVis )  nVar = nVar + 1

       if(nVar == 0) return

       ! Loop over the number of spectral solutions.

       spectralModes: do mm=1,nTimeIntervalsSpectral

         ! Send the variables. The data is first copied into
         ! the send buffer after which the buffer is sent asap.

         ii = 1
         sends: do i=1,commPattern(level)%nProcSend

           ! Store the processor id and the size of the message
           ! a bit easier.

           procID = commPattern(level)%sendProc(i)
           size    = nVar*commPattern(level)%nsend(i)

           ! Copy the data in the correct part of the send buffer.

           jj = ii
           do j=1,commPattern(level)%nsend(i)

             ! Store the block id and the indices of the donor
             ! a bit easier.

             d1 = commPattern(level)%sendList(i)%block(j)
             i1 = commPattern(level)%sendList(i)%indices(j,1)
             j1 = commPattern(level)%sendList(i)%indices(j,2)
             k1 = commPattern(level)%sendList(i)%indices(j,3)

             ! Copy the given range of the working variables for
             ! this cell in the buffer. Update the counter jj.

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

             ! The eddy viscosity for eddy viscosity models.
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
         receives: do i=1,commPattern(level)%nProcRecv

           ! Store the processor id and the size of the message
           ! a bit easier.

           procID = commPattern(level)%recvProc(i)
           size    = nVar*commPattern(level)%nrecv(i)

           ! Post the receive.

           call mpi_irecv(recvBuffer(ii), size, sumb_real, procID, &
                          myID, SUmb_comm_world, recvRequests(i), ierr)

           ! And update ii.

           ii = ii + size

         enddo receives

         ! Copy the local data.

         localCopy: do i=1,internal(level)%ncopy

           ! Store the block and the indices of the donor a bit easier.

           d1 = internal(level)%donorBlock(i)
           i1 = internal(level)%donorIndices(i,1)
           j1 = internal(level)%donorIndices(i,2)
           k1 = internal(level)%donorIndices(i,3)

           ! Idem for the halo's.

           d2 = internal(level)%haloBlock(i)
           i2 = internal(level)%haloIndices(i,1)
           j2 = internal(level)%haloIndices(i,2)
           k2 = internal(level)%haloIndices(i,3)

           ! Copy the given range of working variables.

           do k=start,end
             flowDoms(d2,level,mm)%w(i2,j2,k2,k) = &
             flowDoms(d1,level,mm)%w(i1,j1,k1,k)
           enddo

           ! The pressure, if needed.

           if( commPressure )                   &
             flowDoms(d2,level,mm)%p(i2,j2,k2) = &
             flowDoms(d1,level,mm)%p(i1,j1,k1)

           ! The specific heat ratio, if needed. Note that level == 1.

           if( commVarGamma )                   &
             flowDoms(d2,1,mm)%gamma(i2,j2,k2) = &
             flowDoms(d1,1,mm)%gamma(i1,j1,k1)

           ! The laminar viscosity for viscous computations.
           ! Again level == 1.

           if( commLamVis )                   &
             flowDoms(d2,1,mm)%rlv(i2,j2,k2) = &
             flowDoms(d1,1,mm)%rlv(i1,j1,k1)

           ! The eddy viscosity for eddy viscosity models.
           ! Level is the true multigrid level, because the eddy
           ! viscosity is allocated on all grid levels.

           if( commEddyVis )                      &
             flowDoms(d2,level,mm)%rev(i2,j2,k2) = &
             flowDoms(d1,level,mm)%rev(i1,j1,k1)

         enddo localCopy

         ! Correct the periodic halo's of the internal communication
         ! pattern, if needed.

         if(correctPeriodic .and. internal(level)%nPeriodic > 0)   &
           call correctPeriodicVelocity(level, mm,                 &
                                        internal(level)%nPeriodic, &
                                        internal(level)%periodicData)

         ! Complete the nonblocking receives in an arbitrary sequence and
         ! copy the variables from the buffer into the halo's.

         size = commPattern(level)%nProcRecv
         completeRecvs: do i=1,commPattern(level)%nProcRecv

           ! Complete any of the requests.

           call mpi_waitany(size, recvRequests, index, status, ierr)

           ! Copy the data just arrived in the halo's.

           ii = index
           jj = nVar*commPattern(level)%nrecvCum(ii-1)
           do j=1,commPattern(level)%nrecv(ii)

             ! Store the block and the indices of the halo a bit easier.

             d2 = commPattern(level)%recvList(ii)%block(j)
             i2 = commPattern(level)%recvList(ii)%indices(j,1)
             j2 = commPattern(level)%recvList(ii)%indices(j,2)
             k2 = commPattern(level)%recvList(ii)%indices(j,3)

             ! Copy the conservative variables.

             do k=start,end
               jj = jj + 1
               flowDoms(d2,level,mm)%w(i2,j2,k2,k) = recvBuffer(jj)
             enddo

             ! The pressure, if needed.

             if( commPressure ) then
               jj = jj + 1
               flowDoms(d2,level,mm)%p(i2,j2,k2) = recvBuffer(jj)
             endif

             ! The specific heat ratio, if needed. Note that level == 1.

             if( commVarGamma ) then
               jj = jj + 1
               flowDoms(d2,1,mm)%gamma(i2,j2,k2) = recvBuffer(jj)
             endif

             ! The laminar viscosity for viscous computations.
             ! Again level == 1.

             if( commLamVis ) then
               jj = jj + 1
               flowDoms(d2,1,mm)%rlv(i2,j2,k2) = recvBuffer(jj)
             endif

             ! The eddy viscosity ratio for eddy viscosity models.
             ! Level is the true multigrid level, because the eddy
             ! viscosity is allocated on all grid levels.

             if( commEddyVis ) then
               jj = jj + 1
               flowDoms(d2,level,mm)%rev(i2,j2,k2) = recvBuffer(jj)
             endif

            enddo

         enddo completeRecvs

         ! Correct the periodic halo's of the external communication
         ! pattern, if needed.

         if(correctPeriodic .and. commPattern(level)%nPeriodic > 0)   &
           call correctPeriodicVelocity(level, mm,                    &
                                        commPattern(level)%nPeriodic, &
                                        commPattern(level)%periodicData)

         ! Complete the nonblocking sends.

         size = commPattern(level)%nProcSend
         do i=1,commPattern(level)%nProcSend
           call mpi_waitany(size, sendRequests, index, status, ierr)
         enddo

       enddo spectralModes

       end subroutine whalo1to1

!      ==================================================================

       subroutine correctPeriodicVelocity(level, sp, nPeriodic, &
                                          periodicData)
!
!      ******************************************************************
!      *                                                                *
!      * correctPeriodicVelocity applies the periodic transformation    *
!      * to the velocity of the cell halo's in periodicData.            *
!      *                                                                *
!      ******************************************************************
!
       use block
       use communication
       use constants
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType), intent(in) :: level, sp, nPeriodic
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

         do ii=1,periodicData(nn)%nhalos

           ! Store the block and the indices a bit easier.

           mm = periodicData(nn)%block(ii)
           i  = periodicData(nn)%indices(ii,1)
           j  = periodicData(nn)%indices(ii,2)
           k  = periodicData(nn)%indices(ii,3)

           ! Store the original velocities in vx, vy, vz.

           vx = flowDoms(mm,level,sp)%w(i,j,k,ivx)
           vy = flowDoms(mm,level,sp)%w(i,j,k,ivy)
           vz = flowDoms(mm,level,sp)%w(i,j,k,ivz)

           ! Compute the new velocity vector.

           flowDoms(mm,level,sp)%w(i,j,k,ivx) = rotMatrix(1,1)*vx &
                                              + rotMatrix(1,2)*vy &
                                              + rotMatrix(1,3)*vz
           flowDoms(mm,level,sp)%w(i,j,k,ivy) = rotMatrix(2,1)*vx &
                                              + rotMatrix(2,2)*vy &
                                              + rotMatrix(2,3)*vz
           flowDoms(mm,level,sp)%w(i,j,k,ivz) = rotMatrix(3,1)*vx &
                                              + rotMatrix(3,2)*vy &
                                              + rotMatrix(3,3)*vz
         enddo

       enddo

       end subroutine correctPeriodicVelocity

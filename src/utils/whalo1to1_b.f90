!
!      ******************************************************************
!      *                                                                *
!      * File:          whalo1to1b.f90                                   *
!      * Author:        Gaetan K.W. Kenway                              *
!      * Starting date: 01-22-2015                                      *
!      * Last modified: 01-22-2015                                      *
!      *                                                                *
!      ******************************************************************
!
subroutine whalo1to1_b(level, start, end, commPressure,       &
     commVarGamma, commLamVis, commEddyVis, &
     commPattern, internal)
  !
  !      ******************************************************************
  !      *                                                                *
  !      * whalo1to1b performs the *TRANSPOSE* operation of whalo1to1.    *
  !      * It is used for adjoint/reverse mode residual evaluations.      *
  !      * See whalo1to1 for more information. Note that this code does   *
  !      * include the correctPeroidicVelocity computation.               *
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
  !
  !      ******************************************************************
  !      *                                                                *
  !      * Begin execution                                                *
  !      *                                                                *
  !      ******************************************************************
  !
  ! Determine the number of variables per cell to be sent.

  nVar = max(0_intType,(end - start + 1))
  if( commPressure )  nVar = nVar + 1
  if( commVarGamma ) nVar = nVar + 1
  if( commLamVis )   nVar = nVar + 1
  if( commEddyVis )  nVar = nVar + 1

  if(nVar == 0) return

  ! Loop over the number of spectral solutions.
  spectralModes: do mm=1,nTimeIntervalsSpectral

     ! Gather up the seeds into the *recv* buffer. Note we loop
     ! over nProcRECV here! After the buffer is assembled it is
     ! sent off.

     jj = 1
     ii = 1
     recvs: do i=1,commPattern(level)%nProcRecv

        ! Store the processor id and the size of the message
        ! a bit easier.

        procID = commPattern(level)%sendProc(i)
        size    = nVar*commPattern(level)%nrecv(i)

        ! Copy the data into the buffer

        do j=1,commPattern(level)%nrecv(i)

           ! Store the block and the indices of the halo a bit easier.

           d2 = commPattern(level)%recvList(i)%block(j)
           i2 = commPattern(level)%recvList(i)%indices(j,1)
           j2 = commPattern(level)%recvList(i)%indices(j,2)
           k2 = commPattern(level)%recvList(i)%indices(j,3)

           ! Copy the conservative variables.

           do k=start,end
              recvBuffer(jj) = flowDomsd(d2,level,mm)%w(i2,j2,k2,k)
              jj = jj + 1
           enddo

           ! The pressure, if needed.

           if( commPressure ) then
              recvBuffer(jj) = flowDomsd(d2,level,mm)%p(i2,j2,k2)
              jj = jj + 1
           endif

           ! The specific heat ratio, if needed. Note that level == 1.

           if( commVarGamma ) then
              recvBuffer(jj) = flowDomsd(d2,1,mm)%gamma(i2,j2,k2)
              jj = jj + 1
           endif

           ! The laminar viscosity for viscous computations.
           ! Again level == 1.

           if( commLamVis ) then
              recvBuffer(jj) = flowDomsd(d2,1,mm)%rlv(i2,j2,k2)
              jj = jj + 1
           endif

           ! The eddy viscosity ratio for eddy viscosity models.
           ! Level is the true multigrid level, because the eddy
           ! viscosity is allocated on all grid levels.

           if( commEddyVis ) then
              recvBuffer(jj) = flowDomsd(d2,level,mm)%rev(i2,j2,k2)
              jj = jj + 1
           endif

        enddo

        ! Send the data.
        call mpi_isend(recvBuffer(ii), size, sumb_real, procID,  &
             procID, SUmb_comm_world, recvRequests(i), &
             ierr)

        ! Set ii to jj for the next processor.

        ii = jj

     enddo recvs

     ! Post the nonblocking receives.

     ii = 1
     sends: do i=1,commPattern(level)%nProcSend

        ! Store the processor id and the size of the message
        ! a bit easier.

        procID = commPattern(level)%sendProc(i)
        size    = nVar*commPattern(level)%nsend(i)

        ! Post the receive.

        call mpi_irecv(sendBuffer(ii), size, sumb_real, procID, &
             myID, SUmb_comm_world, sendRequests(i), ierr)

        ! And update ii.

        ii = ii + size

     enddo sends
     
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

        ! Sum into the '1' values from the '2' values (halos). 
        do k=start,end
           flowDomsd(d1,level,mm)%w(i1,j1,k1,k) = flowDomsd(d1,level,mm)%w(i1,j1,k1,k) + &
                flowDomsd(d2,level,mm)%w(i2,j2,k2,k)
        enddo

        ! The pressure, if needed.

        if( commPressure )                   &
             flowDomsd(d1,level,mm)%p(i1,j1,k1) = flowDomsd(d1,level,mm)%p(i1,j1,k1) + &
             flowDomsd(d2,level,mm)%p(i2,j2,k2)

        ! The specific heat ratio, if needed. Note that level == 1.

        if( commVarGamma )                   &
             flowDomsd(d1,1,mm)%gamma(i1,j1,k1) = flowDomsd(d1,1,mm)%gamma(i1,j1,k1) + &
             flowDomsd(d2,1,mm)%gamma(i2,j2,k2)


        ! The laminar viscosity for viscous computations.
        ! Again level == 1.

        if( commLamVis )                   &
             flowDomsd(d1,1,mm)%rlv(i1,j1,k1) = flowDomsd(d1,1,mm)%rlv(i1,j1,k1) + &
             flowDomsd(d2,1,mm)%rlv(i2,j2,k2)

        ! The eddy viscosity for eddy viscosity models.
        ! Level is the true multigrid level, because the eddy
        ! viscosity is allocated on all grid levels.

        if( commEddyVis )                      &
             flowDomsd(d1,level,mm)%rev(i1,j1,k1) = flowDomsd(d1,level,mm)%rev(i1,j1,k1) + &
             flowDomsd(d2,level,mm)%rev(i2,j2,k2)


     enddo localCopy

     ! Complete the nonblocking receives in an arbitrary sequence and
     ! copy the variables from the buffer into the halo's.

     size = commPattern(level)%nProcSend
     completeSends: do i=1,commPattern(level)%nProcSend

        ! Complete any of the requests.

        call mpi_waitany(size, sendRequests, index, status, ierr)

        ! ! Copy the data just arrived in the halo's.

        ii = index

        jj = nVar*commPattern(level)%nsendCum(ii-1)
        
        do j=1,commPattern(level)%nsend(ii)

           ! Store the block and the indices of the halo a bit easier.

           d2 = commPattern(level)%sendList(ii)%block(j)
           i2 = commPattern(level)%sendList(ii)%indices(j,1)
           j2 = commPattern(level)%sendList(ii)%indices(j,2)
           k2 = commPattern(level)%sendList(ii)%indices(j,3)

           ! Copy the conservative variables.

           do k=start,end
              jj = jj + 1
              flowDomsd(d2,level,mm)%w(i2,j2,k2,k) = flowDomsd(d2,level,mm)%w(i2,j2,k2,k) + sendBuffer(jj)
           enddo

           ! The pressure, if needed.

           if( commPressure ) then 
              jj = jj + 1
              flowDomsd(d2,level,mm)%p(i2,j2,k2) = flowDomsd(d2,level,mm)%p(i2,j2,k2) + sendBuffer(jj)
           endif

           ! The specific heat ratio, if needed. Note that level == 1.

           if( commVarGamma ) then
              jj = jj + 1
              flowDomsd(d2,1,mm)%gamma(i2,j2,k2) = flowDomsd(d2,1,mm)%gamma(i2,j2,k2) + sendBuffer(jj)
           endif

           ! The laminar viscosity for viscous computations.
           ! Again level == 1.

           if( commLamVis ) then
              jj = jj + 1
              flowDomsd(d2,1,mm)%rlv(i2,j2,k2) = flowDomsd(d2,1,mm)%rlv(i2,j2,k2) + sendBuffer(jj)
           endif

           ! The eddy viscosity ratio for eddy viscosity models.
           ! Level is the true multigrid level, because the eddy
           ! viscosity is allocated on all grid levels.

           if( commEddyVis ) then
              jj = jj + 1
              flowDomsd(d2,level,mm)%rev(i2,j2,k2) =  flowDomsd(d2,level,mm)%rev(i2,j2,k2) + sendBuffer(jj)
           endif

        enddo

     enddo completeSends

     ! Complete the nonblocking sends.

     size = commPattern(level)%nProcRecv
     do i=1,commPattern(level)%nProcRecv
        call mpi_waitany(size, recvRequests, index, status, ierr)
     enddo

  enddo spectralModes

end subroutine whalo1to1_b

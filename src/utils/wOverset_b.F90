subroutine wOverset_b(level, start, end, commPressure,       &
     commVarGamma, commLamVis, commEddyVis, &
     commPattern, internal, nlev)
  !
  !      ******************************************************************
  !      *                                                                *
  !      * wOverset_b performs the *TRANSPOSE* operation of wOveset       *
  !      * It is used for adjoint/reverse mode residual evaluations.      *
  !      * See wOverset  for more information.
  !      ******************************************************************
  !
  use block
  use communication
  use inputTimeSpectral
  implicit none
  !
  !      Subroutine arguments.
  !
  integer(kind=intType), intent(in) :: level, start, end, nlev
  logical, intent(in) :: commPressure, commVarGamma
  logical, intent(in) :: commLamVis, commEddyVis

  type(commType), dimension(nlev,*), intent(in) :: commPattern
  type(internalCommType), dimension(nlev,*), intent(in) :: internal
  !
  !      Local variables.
  !
  integer :: size, procId, ierr, index
  integer, dimension(mpi_status_size) :: status

  integer(kind=intType) :: nVar, mm
  integer(kind=intType) :: i, j, k, ii, jj
  integer(kind=intType) :: d1, i1, j1, k1, d2, i2, j2, k2

  real(kind=realType), dimension(:), pointer :: weight
  !
  !      ******************************************************************
  !      *                                                                *
  !      * Begin execution                                                *
  !      *                                                                *
  !      ******************************************************************
  !
  ! Determine the number of variables per cell to be sent.

  nVar = max(0_intType,(end - start + 1))
  if( commPressure ) nVar = nVar + 1
  if( commVarGamma ) nVar = nVar + 1
  if( commLamVis )   nVar = nVar + 1
  if( commEddyVis )  nVar = nVar + 1

  if(nVar == 0) return

  ! Loop over the number of spectral solutions.

  spectralModes: do mm=1,nTimeIntervalsSpectral

     ! Gather up the seeds into the *recv* buffer. Note we loop over
     ! nProcRECV here! After the buffer is assembled it is send off.

     jj = 1
     ii = 1
     recvs: do i=1,commPattern(level,mm)%nProcRecv

        ! Store the processor id and the size of the message
        ! a bit easier.

        procID = commPattern(level,mm)%recvProc(i)
        size    = nVar*commPattern(level,mm)%nrecv(i)

        ! Copy the data into the buffer

        do j=1,commPattern(level,mm)%nrecv(i)

           ! Store the block and the indices to make code a bit easier to read
           
           d2 = commPattern(level,mm)%recvList(i)%block(j)
           i2 = commPattern(level,mm)%recvList(i)%indices(j,1)
           j2 = commPattern(level,mm)%recvList(i)%indices(j,2)
           k2 = commPattern(level,mm)%recvList(i)%indices(j,3)


           ! Copy the conservative variables.

           do k=start,end
              recvBuffer(jj) = flowDomsd(d2,level,mm)%w(i2,j2,k2,k)
              jj = jj + 1
              flowDomsd(d2,level,mm)%w(i2,j2,k2,k) = zero
           enddo

           ! The pressure, if needed.

           if( commPressure ) then
              recvBuffer(jj) = flowDomsd(d2,level,mm)%p(i2,j2,k2)
              jj = jj + 1
              flowDomsd(d2,level,mm)%p(i2,j2,k2) = zero
           endif

           ! The specific heat ratio, if needed. Note that level == 1.

           if( commVarGamma ) then
              recvBuffer(jj) = flowDomsd(d2,1,mm)%gamma(i2,j2,k2)
              jj = jj + 1
              flowDomsd(d2,1,mm)%gamma(i2,j2,k2) = zero
           endif

           ! The laminar viscosity for viscous computations.
           ! Again level == 1.

           if( commLamVis ) then
              recvBuffer(jj) = flowDomsd(d2,1,mm)%rlv(i2,j2,k2)
              jj = jj + 1
              flowDomsd(d2,1,mm)%rlv(i2,j2,k2) = zero
           endif

           ! The eddy viscosity ratio for eddy viscosity models.
           ! Level is the true multigrid level, because the eddy
           ! viscosity is allocated on all grid levels.

           if( commEddyVis ) then
              recvBuffer(jj) = flowDomsd(d2,level,mm)%rev(i2,j2,k2)
              jj = jj + 1
              flowDomsd(d2,level,mm)%rev(i2,j2,k2) = zero
           endif

        enddo

        ! Send the data.
        call mpi_isend(recvBuffer(ii), size, sumb_real, procID,  &
             procID, SUmb_comm_world, recvRequests(i), &
             ierr)

        ! Set ii to jj for the next processor.

        ii = jj

     end do recvs

     ! Post the nonblocking receives.

     ii = 1
     sends: do i=1,commPattern(level,mm)%nProcSend

        ! Store the processor id and the size of the message
        ! a bit easier.

        procID = commPattern(level,mm)%sendProc(i)
        size    = nVar*commPattern(level,mm)%nsend(i)

        ! Post the receive.

        call mpi_irecv(sendBuffer(ii), size, sumb_real, procId, &
             myId, SUmb_comm_world, sendRequests(i), ierr)

        ! And update ii.

        ii = ii + size

     enddo sends

     ! Do the local interpolation.

     localInterp: do i=1,internal(level,mm)%ncopy

        ! Store the block and the indices of the donor a bit easier.

        d1 = internal(level,mm)%donorBlock(i)
        i1 = internal(level,mm)%donorIndices(i, 1)
        j1 = internal(level,mm)%donorIndices(i, 2)
        k1 = internal(level,mm)%donorIndices(i, 3)

        weight => internal(level,mm)%donorInterp(i, :)

        ! Idem for the halo's.

        d2 = internal(level,mm)%haloBlock(i)
        i2 = internal(level,mm)%haloIndices(i, 1)
        j2 = internal(level,mm)%haloIndices(i, 2)
        k2 = internal(level,mm)%haloIndices(i, 3)

        ! Sum into the '1' values from the '2' values accouting for the weights

        do k=start,end
           flowDomsd(d1,level,mm)%w(i1  , j1  , k1  , k) = flowDomsd(d1, level, mm)%w(i1  , j1  , k1  , k) + &
                weight(1)*flowDomsd(d2, level, mm)%w(i2, j2, k2, k)

           flowDomsd(d1,level,mm)%w(i1+1, j1  , k1  , k) = flowDomsd(d1, level, mm)%w(i1+1, j1  , k1  , k) + &
                weight(2)*flowDomsd(d2, level, mm)%w(i2, j2, k2, k)

           flowDomsd(d1,level,mm)%w(i1  , j1+1, k1  , k) = flowDomsd(d1, level, mm)%w(i1  , j1+1, k1  , k) + &
                weight(3)*flowDomsd(d2, level, mm)%w(i2, j2, k2, k)

           flowDomsd(d1,level,mm)%w(i1+1, j1+1, k1  , k) = flowDomsd(d1, level, mm)%w(i1+1, j1+1, k1  , k) + &
                weight(4)*flowDomsd(d2, level, mm)%w(i2, j2, k2, k)

           flowDomsd(d1,level,mm)%w(i1  , j1  , k1+1, k) = flowDomsd(d1, level, mm)%w(i1  , j1  , k1+1, k) + &
                weight(5)*flowDomsd(d2, level, mm)%w(i2, j2, k2, k)

           flowDomsd(d1,level,mm)%w(i1+1, j1  , k1+1, k) = flowDomsd(d1, level, mm)%w(i1+1, j1  , k1+1, k) + &
                weight(6)*flowDomsd(d2, level, mm)%w(i2, j2, k2, k)

           flowDomsd(d1,level,mm)%w(i1  , j1+1, k1+1, k) = flowDomsd(d1, level, mm)%w(i1  , j1+1, k1+1, k) + &
                weight(7)*flowDomsd(d2, level, mm)%w(i2, j2, k2, k)

           flowDomsd(d1,level,mm)%w(i1+1, j1+1, k1+1, k) = flowDomsd(d1, level, mm)%w(i1+1, j1+1, k1+1, k) + &
                weight(8)*flowDomsd(d2, level, mm)%w(i2, j2, k2, k)

           flowDomsd(d2, level, mm)%w(i2, j2, k2, k) = zero
        enddo

        ! The pressure, if needed.

        if( commPressure ) then
           flowDomsd(d1,level,mm)%p(i1  , j1  , k1  ) = flowDomsd(d1, level, mm)%p(i1  , j1  , k1  ) + &
                weight(1)*flowDomsd(d2, level, mm)%p(i2, j2, k2)

           flowDomsd(d1,level,mm)%p(i1+1, j1  , k1  ) = flowDomsd(d1, level, mm)%p(i1+1, j1  , k1  ) + &
                weight(2)*flowDomsd(d2, level, mm)%p(i2, j2, k2)

           flowDomsd(d1,level,mm)%p(i1  , j1+1, k1  ) = flowDomsd(d1, level, mm)%p(i1  , j1+1, k1  ) + &
                weight(3)*flowDomsd(d2, level, mm)%p(i2, j2, k2)

           flowDomsd(d1,level,mm)%p(i1+1, j1+1, k1  ) = flowDomsd(d1, level, mm)%p(i1+1, j1+1, k1  ) + &
                weight(4)*flowDomsd(d2, level, mm)%p(i2, j2, k2)

           flowDomsd(d1,level,mm)%p(i1  , j1  , k1+1) = flowDomsd(d1, level, mm)%p(i1  , j1  , k1+1) + &
                weight(5)*flowDomsd(d2, level, mm)%p(i2, j2, k2)

           flowDomsd(d1,level,mm)%p(i1+1, j1  , k1+1) = flowDomsd(d1, level, mm)%p(i1+1, j1  , k1+1) + &
                weight(6)*flowDomsd(d2, level, mm)%p(i2, j2, k2)

           flowDomsd(d1,level,mm)%p(i1  , j1+1, k1+1) = flowDomsd(d1, level, mm)%p(i1  , j1+1, k1+1) + &
                weight(7)*flowDomsd(d2, level, mm)%p(i2, j2, k2)

           flowDomsd(d1,level,mm)%p(i1+1, j1+1, k1+1) = flowDomsd(d1, level, mm)%p(i1+1, j1+1, k1+1) + &
                weight(8)*flowDomsd(d2, level, mm)%p(i2, j2, k2)

           flowDomsd(d2, level, mm)%p(i2, j2, k2) = zero

        end if

        ! The specific heat ratio, if needed. Note that level == 1.

        if( commVarGamma ) then

           flowDomsd(d1,1,mm)%gamma(i1  , j1  , k1  ) = flowDomsd(d1, 1, mm)%gamma(i1  , j1  , k1  ) + &
                weight(1)*flowDomsd(d2, 1, mm)%gamma(i2, j2, k2)

           flowDomsd(d1,1,mm)%gamma(i1+1, j1  , k1  ) = flowDomsd(d1, 1, mm)%gamma(i1+1, j1  , k1  ) + &
                weight(2)*flowDomsd(d2, 1, mm)%gamma(i2, j2, k2)

           flowDomsd(d1,1,mm)%gamma(i1  , j1+1, k1  ) = flowDomsd(d1, 1, mm)%gamma(i1  , j1+1, k1  ) + &
                weight(3)*flowDomsd(d2, 1, mm)%gamma(i2, j2, k2)

           flowDomsd(d1,1,mm)%gamma(i1+1, j1+1, k1  ) = flowDomsd(d1, 1, mm)%gamma(i1+1, j1+1, k1  ) + &
                weight(4)*flowDomsd(d2, 1, mm)%gamma(i2, j2, k2)

           flowDomsd(d1,1,mm)%gamma(i1  , j1  , k1+1) = flowDomsd(d1, 1, mm)%gamma(i1  , j1  , k1+1) + &
                weight(5)*flowDomsd(d2, 1, mm)%gamma(i2, j2, k2)

           flowDomsd(d1,1,mm)%gamma(i1+1, j1  , k1+1) = flowDomsd(d1, 1, mm)%gamma(i1+1, j1  , k1+1) + &
                weight(6)*flowDomsd(d2, 1, mm)%gamma(i2, j2, k2)

           flowDomsd(d1,1,mm)%gamma(i1  , j1+1, k1+1) = flowDomsd(d1, 1, mm)%gamma(i1  , j1+1, k1+1) + &
                weight(7)*flowDomsd(d2, 1, mm)%gamma(i2, j2, k2)

           flowDomsd(d1,1,mm)%gamma(i1+1, j1+1, k1+1) = flowDomsd(d1, 1, mm)%gamma(i1+1, j1+1, k1+1) + &
                weight(8)*flowDomsd(d2, 1, mm)%gamma(i2, j2, k2)

           flowDomsd(d2, 1, mm)%gamma(i2, j2, k2) = zero

        end if

        ! The laminar viscosity for viscous computations.
        ! Again level == 1.

        if( commLamVis ) then

           flowDomsd(d1,1,mm)%rlv(i1  , j1  , k1  ) = flowDomsd(d1, 1, mm)%rlv(i1  , j1  , k1  ) + &
                weight(1)*flowDomsd(d2, 1, mm)%rlv(i2, j2, k2)

           flowDomsd(d1,1,mm)%rlv(i1+1, j1  , k1  ) = flowDomsd(d1, 1, mm)%rlv(i1+1, j1  , k1  ) + &
                weight(2)*flowDomsd(d2, 1, mm)%rlv(i2, j2, k2)

           flowDomsd(d1,1,mm)%rlv(i1  , j1+1, k1  ) = flowDomsd(d1, 1, mm)%rlv(i1  , j1+1, k1  ) + &
                weight(3)*flowDomsd(d2, 1, mm)%rlv(i2, j2, k2)

           flowDomsd(d1,1,mm)%rlv(i1+1, j1+1, k1  ) = flowDomsd(d1, 1, mm)%rlv(i1+1, j1+1, k1  ) + &
                weight(4)*flowDomsd(d2, 1, mm)%rlv(i2, j2, k2)

           flowDomsd(d1,1,mm)%rlv(i1  , j1  , k1+1) = flowDomsd(d1, 1, mm)%rlv(i1  , j1  , k1+1) + &
                weight(5)*flowDomsd(d2, 1, mm)%rlv(i2, j2, k2)

           flowDomsd(d1,1,mm)%rlv(i1+1, j1  , k1+1) = flowDomsd(d1, 1, mm)%rlv(i1+1, j1  , k1+1) + &
                weight(6)*flowDomsd(d2, 1, mm)%rlv(i2, j2, k2)

           flowDomsd(d1,1,mm)%rlv(i1  , j1+1, k1+1) = flowDomsd(d1, 1, mm)%rlv(i1  , j1+1, k1+1) + &
                weight(7)*flowDomsd(d2, 1, mm)%rlv(i2, j2, k2)

           flowDomsd(d1,1,mm)%rlv(i1+1, j1+1, k1+1) = flowDomsd(d1, 1, mm)%rlv(i1+1, j1+1, k1+1) + &
                weight(8)*flowDomsd(d2, 1, mm)%rlv(i2, j2, k2)

           flowDomsd(d2, 1, mm)%rlv(i2, j2, k2) = zero

        end if

        ! The eddy viscosity for eddy viscosity models.
        ! Level is the true multigrid level, because the eddy
        ! viscosity is allocated on all grid levels.

        if( commEddyVis ) then
           flowDomsd(d1,level,mm)%rev(i1  , j1  , k1  ) = flowDomsd(d1, level, mm)%rev(i1  , j1  , k1  ) + &
                weight(1)*flowDomsd(d2, level, mm)%rev(i2, j2, k2)

           flowDomsd(d1,level,mm)%rev(i1+1, j1  , k1  ) = flowDomsd(d1, level, mm)%rev(i1+1, j1  , k1  ) + &
                weight(2)*flowDomsd(d2, level, mm)%rev(i2, j2, k2)

           flowDomsd(d1,level,mm)%rev(i1  , j1+1, k1  ) = flowDomsd(d1, level, mm)%rev(i1  , j1+1, k1  ) + &
                weight(3)*flowDomsd(d2, level, mm)%rev(i2, j2, k2)

           flowDomsd(d1,level,mm)%rev(i1+1, j1+1, k1  ) = flowDomsd(d1, level, mm)%rev(i1+1, j1+1, k1  ) + &
                weight(4)*flowDomsd(d2, level, mm)%rev(i2, j2, k2)

           flowDomsd(d1,level,mm)%rev(i1  , j1  , k1+1) = flowDomsd(d1, level, mm)%rev(i1  , j1  , k1+1) + &
                weight(5)*flowDomsd(d2, level, mm)%rev(i2, j2, k2)

           flowDomsd(d1,level,mm)%rev(i1+1, j1  , k1+1) = flowDomsd(d1, level, mm)%rev(i1+1, j1  , k1+1) + &
                weight(6)*flowDomsd(d2, level, mm)%rev(i2, j2, k2)

           flowDomsd(d1,level,mm)%rev(i1  , j1+1, k1+1) = flowDomsd(d1, level, mm)%rev(i1  , j1+1, k1+1) + &
                weight(7)*flowDomsd(d2, level, mm)%rev(i2, j2, k2)

           flowDomsd(d1,level,mm)%rev(i1+1, j1+1, k1+1) = flowDomsd(d1, level, mm)%rev(i1+1, j1+1, k1+1) + &
                weight(8)*flowDomsd(d2, level, mm)%rev(i2, j2, k2)

           flowDomsd(d2, level, mm)%rev(i2, j2, k2) = zero

        end if

     enddo localInterp

     ! Complete the nonblocking receives in an arbitrary sequence and
     ! copy the variables from the buffer into the halo's.

     size = commPattern(level,mm)%nProcSend
     completeSends: do i=1,commPattern(level,mm)%nProcSend

        ! Complete any of the requests.

        call mpi_waitany(size, sendRequests, index, status, ierr)

        ! Copy the data just arrived in the halo's.

        ii = index

        jj = nVar*commPattern(level,mm)%nsendCum(ii-1)
        do j=1,commPattern(level,mm)%nsend(ii)

           ! Store the block and the indices of the halo a bit easier.

           d2 = commPattern(level,mm)%sendList(ii)%block(j)
           i2 = commPattern(level,mm)%sendList(ii)%indices(j,1)
           j2 = commPattern(level,mm)%sendList(ii)%indices(j,2)
           k2 = commPattern(level,mm)%sendList(ii)%indices(j,3)

           weight => commPattern(level, mm)%sendList(ii)%interp(j, :)
           
           ! Copy the conservative variables.

           do k=start,end
              jj = jj + 1

              flowDomsd(d2,level,mm)%w(i2  , j2  , k2  , k) = flowDomsd(d2, level, mm)%w(i2  , j2  , k2  , k) + &
                   weight(1)*sendBuffer(jj)
              
              flowDomsd(d2,level,mm)%w(i2+1, j2  , k2  , k) = flowDomsd(d2, level, mm)%w(i2+1, j2  , k2  , k) + &
                   weight(2)*sendBuffer(jj)
              
              flowDomsd(d2,level,mm)%w(i2  , j2+1, k2  , k) = flowDomsd(d2, level, mm)%w(i2  , j2+1, k2  , k) + &
                   weight(3)*sendBuffer(jj)
              
              flowDomsd(d2,level,mm)%w(i2+1, j2+1, k2  , k) = flowDomsd(d2, level, mm)%w(i2+1, j2+1, k2  , k) + &
                   weight(4)*sendBuffer(jj)
              
              flowDomsd(d2,level,mm)%w(i2  , j2  , k2+1, k) = flowDomsd(d2, level, mm)%w(i2  , j2  , k2+1, k) + &
                   weight(5)*sendBuffer(jj)
              
              flowDomsd(d2,level,mm)%w(i2+1, j2  , k2+1, k) = flowDomsd(d2, level, mm)%w(i2+1, j2  , k2+1, k) + &
                   weight(6)*sendBuffer(jj)
              
              flowDomsd(d2,level,mm)%w(i2  , j2+1, k2+1, k) = flowDomsd(d2, level, mm)%w(i2  , j2+1, k2+1, k) + &
                   weight(7)*sendBuffer(jj)
              
              flowDomsd(d2,level,mm)%w(i2+1, j2+1, k2+1, k) = flowDomsd(d2, level, mm)%w(i2+1, j2+1, k2+1, k) + &
                   weight(8)*sendBuffer(jj)
           enddo

           ! The pressure, if needed.

           if( commPressure ) then
              jj = jj + 1

              flowDomsd(d2,level,mm)%p(i2  , j2  , k2  ) = flowDomsd(d2, level, mm)%p(i2  , j2  , k2  ) + &
                   weight(1)*sendBuffer(jj)
              
              flowDomsd(d2,level,mm)%p(i2+1, j2  , k2  ) = flowDomsd(d2, level, mm)%p(i2+1, j2  , k2  ) + &
                   weight(2)*sendBuffer(jj)
              
              flowDomsd(d2,level,mm)%p(i2  , j2+1, k2  ) = flowDomsd(d2, level, mm)%p(i2  , j2+1, k2  ) + &
                   weight(3)*sendBuffer(jj)
              
              flowDomsd(d2,level,mm)%p(i2+1, j2+1, k2  ) = flowDomsd(d2, level, mm)%p(i2+1, j2+1, k2  ) + &
                   weight(4)*sendBuffer(jj)
              
              flowDomsd(d2,level,mm)%p(i2  , j2  , k2+1) = flowDomsd(d2, level, mm)%p(i2  , j2  , k2+1) + &
                   weight(5)*sendBuffer(jj)
              
              flowDomsd(d2,level,mm)%p(i2+1, j2  , k2+1) = flowDomsd(d2, level, mm)%p(i2+1, j2  , k2+1) + &
                   weight(6)*sendBuffer(jj)
              
              flowDomsd(d2,level,mm)%p(i2  , j2+1, k2+1) = flowDomsd(d2, level, mm)%p(i2  , j2+1, k2+1) + &
                   weight(7)*sendBuffer(jj)
              
              flowDomsd(d2,level,mm)%p(i2+1, j2+1, k2+1) = flowDomsd(d2, level, mm)%p(i2+1, j2+1, k2+1) + &
                   weight(8)*sendBuffer(jj)

           endif

           ! The specific heat ratio, if needed. Note that level == 1.

           if( commVarGamma ) then
              jj = jj + 1

              flowDomsd(d2,1,mm)%gamma(i2  , j2  , k2  ) = flowDomsd(d2, 1, mm)%gamma(i2  , j2  , k2  ) + &
                   weight(1)*sendBuffer(jj)
              
              flowDomsd(d2,1,mm)%gamma(i2+1, j2  , k2  ) = flowDomsd(d2, 1, mm)%gamma(i2+1, j2  , k2  ) + &
                   weight(2)*sendBuffer(jj)
              
              flowDomsd(d2,1,mm)%gamma(i2  , j2+1, k2  ) = flowDomsd(d2, 1, mm)%gamma(i2  , j2+1, k2  ) + &
                   weight(3)*sendBuffer(jj)
              
              flowDomsd(d2,1,mm)%gamma(i2+1, j2+1, k2  ) = flowDomsd(d2, 1, mm)%gamma(i2+1, j2+1, k2  ) + &
                   weight(4)*sendBuffer(jj)
              
              flowDomsd(d2,1,mm)%gamma(i2  , j2  , k2+1) = flowDomsd(d2, 1, mm)%gamma(i2  , j2  , k2+1) + &
                   weight(5)*sendBuffer(jj)
              
              flowDomsd(d2,1,mm)%gamma(i2+1, j2  , k2+1) = flowDomsd(d2, 1, mm)%gamma(i2+1, j2  , k2+1) + &
                   weight(6)*sendBuffer(jj)
              
              flowDomsd(d2,1,mm)%gamma(i2  , j2+1, k2+1) = flowDomsd(d2, 1, mm)%gamma(i2  , j2+1, k2+1) + &
                   weight(7)*sendBuffer(jj)
              
              flowDomsd(d2,1,mm)%gamma(i2+1, j2+1, k2+1) = flowDomsd(d2, 1, mm)%gamma(i2+1, j2+1, k2+1) + &
                   weight(8)*sendBuffer(jj)

           endif

           ! The laminar viscosity for viscous computations.
           ! Again level == 1.

           if( commLamVis ) then
              jj = jj + 1

              flowDomsd(d2,1,mm)%rlv(i2  , j2  , k2  ) = flowDomsd(d2, 1, mm)%rlv(i2  , j2  , k2  ) + &
                   weight(1)*sendBuffer(jj)
              
              flowDomsd(d2,1,mm)%rlv(i2+1, j2  , k2  ) = flowDomsd(d2, 1, mm)%rlv(i2+1, j2  , k2  ) + &
                   weight(2)*sendBuffer(jj)
              
              flowDomsd(d2,1,mm)%rlv(i2  , j2+1, k2  ) = flowDomsd(d2, 1, mm)%rlv(i2  , j2+1, k2  ) + &
                   weight(3)*sendBuffer(jj)
              
              flowDomsd(d2,1,mm)%rlv(i2+1, j2+1, k2  ) = flowDomsd(d2, 1, mm)%rlv(i2+1, j2+1, k2  ) + &
                   weight(4)*sendBuffer(jj)
              
              flowDomsd(d2,1,mm)%rlv(i2  , j2  , k2+1) = flowDomsd(d2, 1, mm)%rlv(i2  , j2  , k2+1) + &
                   weight(5)*sendBuffer(jj)
              
              flowDomsd(d2,1,mm)%rlv(i2+1, j2  , k2+1) = flowDomsd(d2, 1, mm)%rlv(i2+1, j2  , k2+1) + &
                   weight(6)*sendBuffer(jj)
              
              flowDomsd(d2,1,mm)%rlv(i2  , j2+1, k2+1) = flowDomsd(d2, 1, mm)%rlv(i2  , j2+1, k2+1) + &
                   weight(7)*sendBuffer(jj)
              
              flowDomsd(d2,1,mm)%rlv(i2+1, j2+1, k2+1) = flowDomsd(d2, 1, mm)%rlv(i2+1, j2+1, k2+1) + &
                   weight(8)*sendBuffer(jj)

           endif

           ! The eddy viscosity ratio for eddy viscosity models.
           ! Level is the true multigrid level, because the eddy
           ! viscosity is allocated on all grid levels.

           if( commEddyVis ) then
              jj = jj + 1

              flowDomsd(d2,level,mm)%rev(i2  , j2  , k2  ) = flowDomsd(d2, level, mm)%rev(i2  , j2  , k2  ) + &
                   weight(1)*sendBuffer(jj)
              
              flowDomsd(d2,level,mm)%rev(i2+1, j2  , k2  ) = flowDomsd(d2, level, mm)%rev(i2+1, j2  , k2  ) + &
                   weight(2)*sendBuffer(jj)
              
              flowDomsd(d2,level,mm)%rev(i2  , j2+1, k2  ) = flowDomsd(d2, level, mm)%rev(i2  , j2+1, k2  ) + &
                   weight(3)*sendBuffer(jj)
              
              flowDomsd(d2,level,mm)%rev(i2+1, j2+1, k2  ) = flowDomsd(d2, level, mm)%rev(i2+1, j2+1, k2  ) + &
                   weight(4)*sendBuffer(jj)
              
              flowDomsd(d2,level,mm)%rev(i2  , j2  , k2+1) = flowDomsd(d2, level, mm)%rev(i2  , j2  , k2+1) + &
                   weight(5)*sendBuffer(jj)
              
              flowDomsd(d2,level,mm)%rev(i2+1, j2  , k2+1) = flowDomsd(d2, level, mm)%rev(i2+1, j2  , k2+1) + &
                   weight(6)*sendBuffer(jj)
              
              flowDomsd(d2,level,mm)%rev(i2  , j2+1, k2+1) = flowDomsd(d2, level, mm)%rev(i2  , j2+1, k2+1) + &
                   weight(7)*sendBuffer(jj)
              
              flowDomsd(d2,level,mm)%rev(i2+1, j2+1, k2+1) = flowDomsd(d2, level, mm)%rev(i2+1, j2+1, k2+1) + &
                   weight(8)*sendBuffer(jj)

           endif

        enddo
        
     enddo completeSends

     ! Complete the nonblocking sends.

     size = commPattern(level,mm)%nProcRecv
     do i=1,commPattern(level,mm)%nProcRecv
        call mpi_waitany(size, recvRequests, index, status, ierr)
     enddo

  enddo spectralModes

end subroutine wOverset_b

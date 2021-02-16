module haloExchange

contains

  subroutine whalo1(level, start, end, commPressure, commGamma, &
       commViscous)
    !
    !       whalo1 exchanges all the 1st level internal halo's for the
    !       cell centered variables.
    !
    use constants
    use blockPointers
    use communication
    use flowVarRefState, only : viscous, eddyModel
    use inputPhysics
    use inputTimeSpectral, only : ntimeIntervalsSpectral
    use iteration
    use utils, only : setPointers, getCorrectForK
    use flowUtils, only : computeEtotBlock
    implicit none
    !
    !      Subroutine arguments.
    !
    integer(kind=intType), intent(in) :: level, start, end
    logical, intent(in) :: commPressure, commGamma, commViscous
    !
    !      Local variables.
    !
    integer(kind=intType) :: nn, mm, ll

    logical :: correctForK, commLamVis, commEddyVis, commVarGamma

    ! Set the logicals whether or not to communicate the viscosities.

    commLamVis = .false.
    if(viscous .and. commViscous) commLamVis = .true.

    commEddyVis = .false.
    if(eddyModel .and. commViscous) commEddyVis = .true.

    ! Set the logical whether or not to communicate gamma.

    commVarGamma = .false.
    if(commGamma .and. (cpModel == cpTempCurveFits)) &
         commVarGamma = .true.

    ! Exchange the 1 to 1 matching 1st level cell halo's.

    call whalo1to1(level, start, end, commPressure, commVarGamma, &
         commLamVis, commEddyVis, commPatternCell_1st,  &
         internalCell_1st)

    ! Exchange the overset cells
    call wOverset(level, start, end, commPressure, commVarGamma, &
         commLamVis, commEddyVis, commPatternOverset, internalOverset)

    ! Average any overset orphans.

    do ll=1,nTimeIntervalsSpectral
       do nn=1,nDom
          call setPointers(nn,level,ll)
          call orphanAverage(start, end, commPressure, commGamma, &
               commLamVis, commEddyVis)
       end do
    end do

    ! If both the pressure and the total energy has been communicated
    ! compute the energy again. The reason is that both values are
    ! interpolated and consequently the values are not consistent.
    ! The energy depends quadratically on the velocity.

    bothPAndE: if(commPressure .and. start <= irhoE .and. &
         end >= irhoE) then

       ! First determine whether or not the total energy must be
       ! corrected for the presence of the turbulent kinetic energy.
       correctForK = getCorrectForK()

       ! Loop over the blocks to find the sliding mesh subfaces.
       ! Use is made of the fact the boundary conditions are identical
       ! for all spectral solutions. So that loop can be inside the
       ! test for the sliding mesh subface.

       domains: do nn=1,nDom
          do mm=1,flowDoms(nn,level,1)%nBocos
             if(flowDoms(nn,level,1)%BCType(mm) == slidingInterface) then

                ! Loop over the number of spectral solutions.

                do ll=1,nTimeIntervalsSpectral

                   ! Set the pointers for this block and compute the energy
                   ! for the halo cells of this sliding interface subface.

                   call setPointers(nn,level,ll)
                   call computeEtotBlock(icBeg(mm), icEnd(mm), &
                        jcBeg(mm), jcEnd(mm), &
                        kcBeg(mm), kcEnd(mm), correctForK)
                enddo
             endif
          enddo

       enddo domains

    endif bothPAndE

  end subroutine whalo1

  subroutine whalo2(level, start, end, commPressure, commGamma, &
       commViscous)
    !
    !       whalo2 exchanges all the 2nd level internal halo's for the
    !       cell centered variables.
    !
    use constants
    use blockPointers
    use communication
    use flowVarRefState
    use inputPhysics
    use inputTimeSpectral
    use iteration
    use utils, only : setPointers, getCorrectForK
    use flowUtils, only : computeEtotBlock
    implicit none
    !
    !      Subroutine arguments.
    !
    integer(kind=intType), intent(in) :: level, start, end
    logical, intent(in) :: commPressure, commGamma, commViscous
    !
    !      Local variables.
    !
    integer(kind=intType) :: nn, ll
    integer(kind=intType) :: iBeg, iEnd, jBeg, jEnd, kBeg, kEnd

    logical :: correctForK, commLamVis, commEddyVis, commVarGamma

    ! Set the logicals whether or not to communicate the viscosities.

    commLamVis = .false.
    if(viscous .and. commViscous) commLamVis = .true.

    commEddyVis = .false.
    if(eddyModel .and. commViscous) commEddyVis = .true.

    ! Set the logical whether or not to communicate gamma.

    commVarGamma = .false.
    if(commGamma .and. (cpModel == cpTempCurveFits)) &
         commVarGamma = .true.

    ! Exchange the 1 to 1 matching 2nd level cell halo's.

    call whalo1to1(level, start, end, commPressure, commVarGamma, &
         commLamVis, commEddyVis, commPatternCell_2nd,  &
         internalCell_2nd)

    ! Exchange the overset cells
    call wOverset(level, start, end, commPressure, commVarGamma, &
         commLamVis, commEddyVis, commPatternOverset, internalOverset)

    ! Average any overset orphans.

    do ll=1,nTimeIntervalsSpectral
       do nn=1,nDom
          call setPointers(nn,level,ll)
          call orphanAverage(start, end, commPressure, commGamma, &
               commLamVis, commEddyVis)
       end do
    end do

    ! If both the pressure and the total energy has been communicated
    ! compute the energy again. The reason is that both values are
    ! interpolated and consequently the values are not consistent.
    ! The energy depends quadratically on the velocity.

    bothPAndE: if(commPressure .and. start <= irhoE .and. &
         end >= irhoE) then

       ! First determine whether or not the total energy must be
       ! corrected for the presence of the turbulent kinetic energy.

       correctForK = getCorrectForK()

       domains: do nn=1,nDom

          ! Treat the overset blocks. Since we don't have the logic
          ! setup here correctly to only update the overset cells,
          ! just do the whole block, for every block
          do ll=1, nTimeIntervalsSpectral
             call setPointers(nn, level, ll)
             call computeETotBlock(2, il, 2, jl, 2, kl, correctForK)
          end do


       enddo domains

    endif bothPAndE

  end subroutine whalo2

  subroutine orphanAverage(wstart, wend, calcPressure, calcGamma, &
       calcLamVis, calcEddyVis)
    !
    !       orphanAverage uses the neighboring cells of an overset orphan
    !       to set the flow state for the orphan cell by a simple average.
    !       This routine operates on the block given by the block pointers
    !       so it is assumed they are set.
    !
    use constants
    use blockPointers
    use flowVarRefState
    use inputPhysics
    implicit none
    !
    !      Subroutine arguments.
    !
    integer(kind=intType), intent(in) :: wstart, wend

    logical, intent(in) :: calcPressure, calcGamma, calcLamVis
    logical, intent(in) :: calcEddyVis
    !
    !      Local variables.
    !
    integer(kind=intType) :: oi, oj, ok, ni, nj, nk, i, l, m, n, nAvg

    integer(kind=intType), dimension(3) :: del

    real(kind=realType) :: nAvgReal

    ! Return immediately if there are no orphans for this block.

    if (nOrphans == 0) return

    ! Loop over the number of orphans.

    orphanLoop: do n = 1,nOrphans

       ! Store the orphan indices easier.

       oi = orphans(1, n)
       oj = orphans(2, n)
       ok = orphans(3, n)

       ! Initialize the number of neighbors used to 0 and also set
       ! the flow variables to zero such that an average can be
       ! accumulated below.

       nAvg = 0

       do l=wstart, wend
          w(oi,oj,ok,l) = zero
       end do

       if (calcPressure) p(oi,oj,ok)     = zero
       if (calcGamma)    gamma(oi,oj,ok) = zero
       if (calcLamVis)   rlv(oi,oj,ok)   = zero
       if (calcEddyVis)  rev(oi,oj,ok)   = zero

       ! Loop over the 3 coordinate directions, and for both the
       ! positive and negative direction set the delta vector to be a
       ! unit vector in that direction.

       directionLoop: do m = 1,3
          plusMinusLoop: do i = -1,1,2

             del    = 0
             del(m) = i

             ! Compute the neighbor indices and skip if it is outside the
             ! boundaries of the block.

             ni = oi + del(1)
             nj = oj + del(2)
             nk = ok + del(3)

             if (ni < 0 .or. ni > ib .or. &
                  nj < 0 .or. nj > jb .or. &
                  nk < 0 .or. nk > kb) cycle

             ! If the neighboring iblank value indicates a cell that is
             ! part of either the field, fringe, or boundary condition,
             ! then use its flow state in the average.

             if (iblank(ni,nj,nk) == 1) then

                ! Update the number of neighbors used in the average and
                ! compute the flow variables for the given range.

                nAvg = nAvg + 1

                do l=wstart,wend
                   w(oi,oj,ok,l) = w(oi,oj,ok,l) + w(ni,nj,nk,l)
                end do

                ! Check if the pressure, specific heat ratio, laminar
                ! viscosity, and/or eddy viscosity needs to be computed.

                if (calcPressure) &
                     p(oi,oj,ok) = p(oi,oj,ok) + p(ni,nj,nk)
                if (calcGamma) &
                     gamma(oi,oj,ok) = gamma(oi,oj,ok) + gamma(ni,nj,nk)
                if (calcLamVis) &
                     rlv(oi,oj,ok) = rlv(oi,oj,ok) + rlv(ni,nj,nk)
                if (calcEddyVis) &
                     rev(oi,oj,ok) = rev(oi,oj,ok) + rev(ni,nj,nk)

             end if

          end do plusMinusLoop
       end do directionLoop

       ! Check to make sure that at least 1 suitable neighbeor was
       ! found to use in the average.

       checkNoNeighbors: if (nAvg > 0) then

          ! Divide each of the variables being computed by the number
          ! of neighbors used in the average.

          nAvgReal = real(nAvg, realType)

          ! Average the flow variables for the given range.

          do l=wstart,wend
             w(oi,oj,ok,l) = w(oi,oj,ok,l)/nAvgReal
          end do

          ! Check if the pressure, specific heat ratio, laminar
          ! viscosity, and/or eddy viscosity needs to be averaged.

          if (calcPressure) p(oi,oj,ok)     = p(oi,oj,ok)/nAvgReal
          if (calcGamma)    gamma(oi,oj,ok) = gamma(oi,oj,ok)/nAvgReal
          if (calcLamVis)   rlv(oi,oj,ok)   = rlv(oi,oj,ok)/nAvgReal
          if (calcEddyVis)  rev(oi,oj,ok)   = rev(oi,oj,ok)/nAvgReal

       else checkNoNeighbors

          ! No suitable neighbors were found in order to compute an
          ! average. Set the variables back to the the freestream.

          do l=wstart,wend
             w(oi,oj,ok,l) = wInf(l)
          end do

          if (calcPressure) p(oi,oj,ok)     = pInfCorr
          if (calcGamma)    gamma(oi,oj,ok) = gammaInf
          if (calcLamVis)   rlv(oi,oj,ok)   = muInf
          if (calcEddyVis)  rev(oi,oj,ok)   = eddyVisInfRatio*muInf

       end if checkNoNeighbors

    end do orphanLoop

  end subroutine orphanAverage

  subroutine setCommPointers(start, end, commPressure, commVarGamma, commLamVis, &
       commEddyVis, level, sps, derivPointers, nVar, varOffset)

    ! Generic routine for setting pointers to the communication
    ! variables. Can also set pointers to derivatve values if derivPts is True.

    use constants
    use block, only : fLowDoms, blockType, flowDomsd, nDom
    implicit none

    ! Input
    integer(kind=intType), intent(in) :: start, end, level, sps
    logical, intent(in) :: commPressure, commVarGamma, commLamVis, commEddyVis
    logical, intent(in) :: derivPointers
    integer(kind=intType), intent(in) :: varOffset

    ! Output
    integer(kind=intType), intent(out) :: nVar

    ! Working:
    integer(kind=intType) :: nn, k
    type(blockType) , pointer :: blk, blk1, blkLevel

    ! Set the pointers for the required variables
    domainLoop:do nn=1, nDom
       nVar = varOffset
       blk => flowDoms(nn, level, sps)

       if (derivPointers) then
          blkLevel => flowDomsd(nn, level, sps)
          blk1     => flowDomsd(nn, 1    , sps)
       else
          blkLevel => flowDoms(nn, level, sps)
          blk1     => flowDoms(nn, 1    , sps)
       end if

       do k=start, end
          nVar = nVar + 1
          blk%realCommVars(nVar)%var => blkLevel%w(:, :, :, k)
       end do

       if( commPressure )  then
          nVar = nVar + 1
          blk%realCommVars(nVar)%var => blkLevel%P(:, :, :)
       end if

       if( commVarGamma ) then
          nVar = nVar + 1
          blk%realCommVars(nVar)%var => blk1%gamma(:, :, :)
       end if

       if( commLamVis ) then
          nVar = nVar + 1
          blk%realCommvars(nVar)%var => blk1%rlv(:, :, :)
       end if

       if( commEddyVis ) then
          nVar = nVar + 1
          blk%realCommVars(nVar)%var => blkLevel%rev(:, :, :)
       end if

    end do domainLoop
    nVar = nVar - varOffset
  end subroutine setCommPointers

  subroutine whalo1to1(level, start, end, commPressure,       &
       commVarGamma, commLamVis, commEddyVis, &
       commPattern, internal)
    !
    !       whalo1to1 exchanges the 1 to 1 internal halo's for the cell
    !       centered variables for the given communication pattern. It
    !       is possible to send a range of variables and not the entire
    !       set, e.g. only the flow variables or only the turbulent
    !       variables. This is controlled by the arguments start, end,
    !       commPressure and commViscous. The exchange takes place for
    !       the given grid level.
    !
    use constants
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

    integer(kind=intType) :: nVar, nn, k, sps

    logical :: correctPeriodic

    ! Set the logical correctPeriodic. Only if a momentum variable
    ! is communicated it is needed to apply the periodic
    ! transformations.

    correctPeriodic = .false.
    if(start <= ivx .and. end >= ivz) correctPeriodic = .true.

    spectralModes: do sps=1,nTimeIntervalsSpectral

       call setCommPointers(start, end, commPressure, commVarGamma, &
            commLamVis, commEddyVis, level, sps, .False., nVar, 0)

       if (nVar == 0) then
          return
       end if

       ! Run the generic exchange
       call wHalo1to1RealGeneric(nVar, level, sps, commPattern, internal)

       if (correctPeriodic) then
          if ( internal(level)%nPeriodic > 0 )   then
             call correctPeriodicVelocity(level, sps,             &
                  internal(level)%nPeriodic, internal(level)%periodicData)
          end if

          if ( commPattern(level)%nPeriodic > 0 )   then
             call correctPeriodicVelocity(level, sps,             &
                  commPattern(level)%nPeriodic, commPattern(level)%periodicData)
          end if
       end if

    end do spectralModes

  end subroutine whalo1to1

  subroutine correctPeriodicVelocity(level, sp, nPeriodic, &
       periodicData)
    !
    !       correctPeriodicVelocity applies the periodic transformation
    !       to the velocity of the cell halo's in periodicData.
    !
    use constants
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

    ! Loop over the number of periodic transformations.

    do nn=1,nPeriodic

       ! Store the rotation matrix a bit easier.

       rotMatrix = periodicData(nn)%rotMatrix

       ! Loop over the number of halo cells for this transformation.
       !DIR$ NOVECTOR
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

  subroutine whalo1to1RealGeneric(nVar, level, sps, commPattern, internal)
    !
    !       whalo1to1 exchanges the 1 to 1 internal halo's for the cell
    !       centered variables for the given communication pattern.
    !       Pointers must be set for var1, var2...varN
    !
    use constants
    use block
    use communication
    use inputTimeSpectral
    use utils, only : EChk
    implicit none
    !
    !      Subroutine arguments.
    !
    integer(kind=intType), intent(in) :: level, sps

    type(commType), dimension(*), intent(in)         :: commPattern
    type(internalCommType), dimension(*), intent(in) :: internal
    !
    !      Local variables.
    !

    integer :: size, procID, ierr, index
    integer, dimension(mpi_status_size) :: mpiStatus

    integer(kind=intType) :: nVar, mm
    integer(kind=intType) :: i, j, k, ii, jj
    integer(kind=intType) :: d1, i1, j1, k1, d2, i2, j2, k2

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
       !DIR$ NOVECTOR
       do j=1,commPattern(level)%nsend(i)

          ! Store the block id and the indices of the donor
          ! a bit easier.

          d1 = commPattern(level)%sendList(i)%block(j)
          i1 = commPattern(level)%sendList(i)%indices(j,1)+1
          j1 = commPattern(level)%sendList(i)%indices(j,2)+1
          k1 = commPattern(level)%sendList(i)%indices(j,3)+1

          ! Copy the given range of the working variables for
          ! this cell in the buffer. Update the counter jj.
          !DIR$ NOVECTOR
          do k=1, nvar
             sendBuffer(jj) = flowDoms(d1, level, sps)%realCommVars(k)%var(i1, j1, k1)
             jj = jj + 1
          end do
       end do

       ! Send the data.

       call mpi_isend(sendBuffer(ii), size, adflow_real, procID,  &
            procID, ADflow_comm_world, sendRequests(i), &
            ierr)
       call EChk(ierr,__FILE__,__LINE__)

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

       call mpi_irecv(recvBuffer(ii), size, adflow_real, procID, &
            myID, ADflow_comm_world, recvRequests(i), ierr)
       call EChk(ierr,__FILE__,__LINE__)

       ! And update ii.

       ii = ii + size

    enddo receives

    ! Copy the local data.


    !DIR$ NOVECTOR
    localCopy: do i=1,internal(level)%ncopy

       ! Store the block and the indices of the donor a bit easier.

       d1 = internal(level)%donorBlock(i)
       i1 = internal(level)%donorIndices(i,1)+1
       j1 = internal(level)%donorIndices(i,2)+1
       k1 = internal(level)%donorIndices(i,3)+1

       ! Idem for the halo's.

       d2 = internal(level)%haloBlock(i)
       i2 = internal(level)%haloIndices(i,1)+1
       j2 = internal(level)%haloIndices(i,2)+1
       k2 = internal(level)%haloIndices(i,3)+1

       do k=1, nVar
          flowDoms(d2, level, sps)%realCommVars(k)%var(i2, j2, k2) = &
               flowDoms(d1, level, sps)%realCommVars(k)%var(i1, j1, k1)
       end do

    enddo localCopy

    ! Complete the nonblocking receives in an arbitrary sequence and
    ! copy the variables from the buffer into the halo's.

    size = commPattern(level)%nProcRecv
    completeRecvs: do i=1,commPattern(level)%nProcRecv

       ! Complete any of the requests.

       call mpi_waitany(size, recvRequests, index, mpiStatus, ierr)
       call EChk(ierr,__FILE__,__LINE__)

       ! Copy the data just arrived in the halo's.

       ii = index
       jj = nVar*commPattern(level)%nrecvCum(ii-1)
       !DIR$ NOVECTOR
       do j=1,commPattern(level)%nrecv(ii)

          ! Store the block and the indices of the halo a bit easier.

          d2 = commPattern(level)%recvList(ii)%block(j)
          i2 = commPattern(level)%recvList(ii)%indices(j,1)+1
          j2 = commPattern(level)%recvList(ii)%indices(j,2)+1
          k2 = commPattern(level)%recvList(ii)%indices(j,3)+1

          do k=1, nVar
             jj = jj + 1
             flowDoms(d2,level,sps)%realCommVars(k)%var(i2, j2, k2) = recvBuffer(jj)
          end do
       end do
    enddo completeRecvs

    ! Complete the nonblocking sends.

    size = commPattern(level)%nProcSend
    do i=1,commPattern(level)%nProcSend
       call mpi_waitany(size, sendRequests, index, mpiStatus, ierr)
    enddo

  end subroutine whalo1to1RealGeneric

  subroutine whalo1to1RealGeneric_b(nVar, level, sps, commPattern, internal)
    !
    !       whalo1to1RealGeneric_b is a generic implementation
    !       of the reverse mode of whalo1to1RealGeneric.
    !
    use constants
    use block
    use communication
    use inputTimeSpectral
    use utils, only : EChk
    implicit none
    !
    !      Subroutine arguments.
    !
    integer(kind=intType), intent(in) :: level, sps

    type(commType), dimension(*), intent(in)         :: commPattern
    type(internalCommType), dimension(*), intent(in) :: internal
    !
    !      Local variables.
    !

    integer :: size, procID, ierr, index
    integer, dimension(mpi_status_size) :: mpiStatus

    integer(kind=intType) :: nVar, mm
    integer(kind=intType) :: i, j, k, ii, jj
    integer(kind=intType) :: d1, i1, j1, k1, d2, i2, j2, k2

    ! Gather up the seeds into the *recv* buffer. Note we loop
    ! over nProcRECV here! After the buffer is assembled it is
    ! sent off.

    jj = 1
    ii = 1
    recvs: do i=1,commPattern(level)%nProcRecv

       ! Store the processor id and the size of the message
       ! a bit easier.

       procID = commPattern(level)%recvProc(i)
       size    = nVar*commPattern(level)%nrecv(i)

       ! Copy the data into the buffer
       !DIR$ NOVECTOR
       do j=1,commPattern(level)%nrecv(i)

          ! Store the block and the indices of the halo a bit easier.

          d2 = commPattern(level)%recvList(i)%block(j)
          i2 = commPattern(level)%recvList(i)%indices(j,1)+1
          j2 = commPattern(level)%recvList(i)%indices(j,2)+1
          k2 = commPattern(level)%recvList(i)%indices(j,3)+1

          do k=1, nVar
             recvBuffer(jj) = flowDoms(d2,level,sps)%realCommVars(k)%var(i2,j2,k2)
             flowDoms(d2,level,sps)%realCommVars(k)%var(i2, j2, k2) = zero
             jj = jj + 1
          enddo
       end do

       ! Send the data.
       call mpi_isend(recvBuffer(ii), size, adflow_real, procID,  &
            procID, ADflow_comm_world, recvRequests(i), &
            ierr)
       call EChk(ierr,__FILE__,__LINE__)

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

       call mpi_irecv(sendBuffer(ii), size, adflow_real, procID, &
            myID, ADflow_comm_world, sendRequests(i), ierr)
       call EChk(ierr,__FILE__,__LINE__)

       ! And update ii.

       ii = ii + size

    enddo sends

    ! Copy the local data.
    !DIR$ NOVECTOR
    localCopy: do i=1,internal(level)%ncopy

       ! Store the block and the indices of the donor a bit easier.

       d1 = internal(level)%donorBlock(i)
       i1 = internal(level)%donorIndices(i,1)+1
       j1 = internal(level)%donorIndices(i,2)+1
       k1 = internal(level)%donorIndices(i,3)+1

       ! Idem for the halo's.

       d2 = internal(level)%haloBlock(i)
       i2 = internal(level)%haloIndices(i,1)+1
       j2 = internal(level)%haloIndices(i,2)+1
       k2 = internal(level)%haloIndices(i,3)+1

       ! Sum into the '1' values from the '2' values (halos).
       do k=1, nVar
          flowDoms(d1, level, sps)%realCommVars(k)%var(i1, j1, k1) = &
               flowDoms(d1, level, sps)%realCommVars(k)%var(i1, j1, k1) + &
               flowDoms(d2, level, sps)%realCommVars(k)%var(i2, j2, k2)
          flowDoms(d2, level, sps)%realCommVars(k)%var(i2, j2, k2) = zero
       enddo

    enddo localCopy

    ! Complete the nonblocking receives in an arbitrary sequence and
    ! copy the variables from the buffer into the halo's.

    size = commPattern(level)%nProcSend
    completeSends: do i=1,commPattern(level)%nProcSend

       ! Complete any of the requests.

       call mpi_waitany(size, sendRequests, index, mpiStatus, ierr)
       call EChk(ierr,__FILE__,__LINE__)

       ! ! Copy the data just arrived in the halo's.

       ii = index
       jj = nVar*commPattern(level)%nsendCum(ii-1)
       !DIR$ NOVECTOR
       do j=1,commPattern(level)%nsend(ii)

          ! Store the block and the indices of the halo a bit easier.

          d2 = commPattern(level)%sendList(ii)%block(j)
          i2 = commPattern(level)%sendList(ii)%indices(j,1)+1
          j2 = commPattern(level)%sendList(ii)%indices(j,2)+1
          k2 = commPattern(level)%sendList(ii)%indices(j,3)+1

          ! Copy the conservative variables.

          do k=1, nVar
             jj = jj + 1
             flowDoms(d2, level, sps)%realCommVars(k)%var(i2, j2, k2) = &
                  flowDoms(d2, level, sps)%realCommVars(k)%var(i2, j2, k2) + sendBuffer(jj)
          enddo
       enddo

    enddo completeSends

    ! Complete the nonblocking sends.

    size = commPattern(level)%nProcRecv
    do i=1,commPattern(level)%nProcRecv
       call mpi_waitany(size, recvRequests, index, mpiStatus, ierr)
       call EChk(ierr,__FILE__,__LINE__)
    enddo

  end subroutine whalo1to1RealGeneric_b

  subroutine whalo1to1IntGeneric(nVar, level, sps, commPattern, internal)
    !
    !       whalo1to1 exchanges the 1 to 1 internal halo's for the cell
    !       centered variables for the given communication pattern.
    !       Pointers must be set for var1, var2...varN
    !
    use constants
    use block
    use communication
    use inputTimeSpectral
    use utils, only : EChk
    implicit none
    !
    !      Subroutine arguments.
    !
    integer(kind=intType), intent(in) :: level, sps

    type(commType), dimension(*), intent(in)         :: commPattern
    type(internalCommType), dimension(*), intent(in) :: internal
    !
    !      Local variables.
    !

    integer :: size, procID, ierr, index
    integer, dimension(mpi_status_size) :: mpiStatus

    integer(kind=intType) :: nVar, mm
    integer(kind=intType) :: i, j, k, ii, jj
    integer(kind=intType) :: d1, i1, j1, k1, d2, i2, j2, k2

    integer(kind=intType), dimension(:), allocatable :: sendBufInt
    integer(kind=intType), dimension(:), allocatable :: recvBufInt

    ii = commPattern(level)%nProcSend
    ii = commPattern(level)%nsendCum(ii)
    jj = commPattern(level)%nProcRecv
    jj = commPattern(level)%nrecvCum(jj)

    allocate(sendBufInt(ii*nVar), recvBufInt(jj*nVar), stat=ierr)

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
       !DIR$ NOVECTOR
       do j=1,commPattern(level)%nsend(i)

          ! Store the block id and the indices of the donor
          ! a bit easier.

          d1 = commPattern(level)%sendList(i)%block(j)
          i1 = commPattern(level)%sendList(i)%indices(j,1)+1
          j1 = commPattern(level)%sendList(i)%indices(j,2)+1
          k1 = commPattern(level)%sendList(i)%indices(j,3)+1

          ! Copy the given range of the working variables for
          ! this cell in the buffer. Update the counter jj.
          !DIR$ NOVECTOR
          do k=1, nvar
             sendBufInt(jj) = flowDoms(d1, level, sps)%intCommVars(k)%var(i1, j1, k1)
             jj = jj + 1
          end do
       end do

       ! Send the data.

       call mpi_isend(sendBufInt(ii), size, adflow_integer, procID,  &
            procID, ADflow_comm_world, sendRequests(i), &
            ierr)
       call EChk(ierr,__FILE__,__LINE__)

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

       call mpi_irecv(recvBufInt(ii), size, adflow_integer, procID, &
            myID, ADflow_comm_world, recvRequests(i), ierr)
       call EChk(ierr,__FILE__,__LINE__)

       ! And update ii.

       ii = ii + size

    enddo receives

    ! Copy the local data.
    !DIR$ NOVECTOR
    localCopy: do i=1,internal(level)%ncopy

       ! Store the block and the indices of the donor a bit easier.

       d1 = internal(level)%donorBlock(i)
       i1 = internal(level)%donorIndices(i,1)+1
       j1 = internal(level)%donorIndices(i,2)+1
       k1 = internal(level)%donorIndices(i,3)+1

       ! Idem for the halo's.

       d2 = internal(level)%haloBlock(i)
       i2 = internal(level)%haloIndices(i,1)+1
       j2 = internal(level)%haloIndices(i,2)+1
       k2 = internal(level)%haloIndices(i,3)+1

       do k=1, nVar
          flowDoms(d2, level, sps)%intCommVars(k)%var(i2, j2, k2) = &
               flowDoms(d1, level, sps)%intCommVars(k)%var(i1, j1, k1)
       end do

    enddo localCopy

    ! Complete the nonblocking receives in an arbitrary sequence and
    ! copy the variables from the buffer into the halo's.

    size = commPattern(level)%nProcRecv
    completeRecvs: do i=1,commPattern(level)%nProcRecv

       ! Complete any of the requests.

       call mpi_waitany(size, recvRequests, index, mpiStatus, ierr)
       call EChk(ierr,__FILE__,__LINE__)

       ! Copy the data just arrived in the halo's.

       ii = index
       jj = nVar*commPattern(level)%nrecvCum(ii-1)
       !DIR$ NOVECTOR
       do j=1,commPattern(level)%nrecv(ii)

          ! Store the block and the indices of the halo a bit easier.

          d2 = commPattern(level)%recvList(ii)%block(j)
          i2 = commPattern(level)%recvList(ii)%indices(j,1)+1
          j2 = commPattern(level)%recvList(ii)%indices(j,2)+1
          k2 = commPattern(level)%recvList(ii)%indices(j,3)+1

          do k=1, nVar
             jj = jj + 1
             flowDoms(d2,level,sps)%intCommVars(k)%var(i2, j2, k2) = recvBufInt(jj)
          end do
       end do
    enddo completeRecvs

    ! Complete the nonblocking sends.

    size = commPattern(level)%nProcSend
    do i=1,commPattern(level)%nProcSend
       call mpi_waitany(size, sendRequests, index, mpiStatus, ierr)
       call EChk(ierr,__FILE__,__LINE__)
    enddo

    deallocate(recvBufInt, sendBufInt)

  end subroutine whalo1to1IntGeneric

  subroutine whalo1to1IntGeneric_b(nVar, level, sps, commPattern, internal)
    !
    !       whalo1to1IntGeneric_b is a generic implementation of the
    !       reverse mode of whalo1to1IntGeneric. Integers are summed
    !       together in reverse.
    !
    use constants
    use block
    use communication
    use inputTimeSpectral
    use utils, only : EChk
    implicit none
    !
    !      Subroutine arguments.
    !
    integer(kind=intType), intent(in) :: level, sps

    type(commType), dimension(*), intent(in)         :: commPattern
    type(internalCommType), dimension(*), intent(in) :: internal
    !
    !      Local variables.
    !

    integer :: size, procID, ierr, index
    integer, dimension(mpi_status_size) :: mpiStatus

    integer(kind=intType) :: nVar, mm
    integer(kind=intType) :: i, j, k, ii, jj
    integer(kind=intType) :: d1, i1, j1, k1, d2, i2, j2, k2
    integer(kind=intType), dimension(:), allocatable :: sendBufInt
    integer(kind=intType), dimension(:), allocatable :: recvBufInt

    ii = commPattern(level)%nProcSend
    ii = commPattern(level)%nsendCum(ii)
    jj = commPattern(level)%nProcRecv
    jj = commPattern(level)%nrecvCum(jj)

    allocate(sendBufInt(ii*nVar), recvBufInt(jj*nVar), stat=ierr)

    ! Gather up the seeds into the *recv* buffer. Note we loop
    ! over nProcRECV here! After the buffer is assembled it is
    ! sent off.

    jj = 1
    ii = 1
    recvs: do i=1,commPattern(level)%nProcRecv

       ! Store the processor id and the size of the message
       ! a bit easier.

       procID = commPattern(level)%recvProc(i)
       size    = nVar*commPattern(level)%nrecv(i)

       ! Copy the data into the buffer
       !DIR$ NOVECTOR
       do j=1,commPattern(level)%nrecv(i)

          ! Store the block and the indices of the halo a bit easier.

          d2 = commPattern(level)%recvList(i)%block(j)
          i2 = commPattern(level)%recvList(i)%indices(j,1)+1
          j2 = commPattern(level)%recvList(i)%indices(j,2)+1
          k2 = commPattern(level)%recvList(i)%indices(j,3)+1
          !DIR$ NOVECTOR
          do k=1, nVar
             recvBufInt(jj) = flowDoms(d2,level,sps)%intCommVars(k)%var(i2,j2,k2)
             flowDoms(d2,level,sps)%intCommVars(k)%var(i2, j2, k2) = 0
             jj = jj + 1
          enddo
       end do

       ! Send the data.
       call mpi_isend(recvBufInt(ii), size, adflow_integer, procID,  &
            procID, ADflow_comm_world, recvRequests(i), &
            ierr)
       call EChk(ierr,__FILE__,__LINE__)

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

       call mpi_irecv(sendBufInt(ii), size, adflow_integer, procID, &
            myID, ADflow_comm_world, sendRequests(i), ierr)
       call EChk(ierr,__FILE__,__LINE__)

       ! And update ii.

       ii = ii + size

    enddo sends

    ! Copy the local data.
    !DIR$ NOVECTOR
    localCopy: do i=1,internal(level)%ncopy

       ! Store the block and the indices of the donor a bit easier.

       d1 = internal(level)%donorBlock(i)
       i1 = internal(level)%donorIndices(i,1)+1
       j1 = internal(level)%donorIndices(i,2)+1
       k1 = internal(level)%donorIndices(i,3)+1

       ! Idem for the halo's.

       d2 = internal(level)%haloBlock(i)
       i2 = internal(level)%haloIndices(i,1)+1
       j2 = internal(level)%haloIndices(i,2)+1
       k2 = internal(level)%haloIndices(i,3)+1

       ! Sum into the '1' values from the '2' values (halos).
       !DIR$ NOVECTOR
       do k=1, nVar
          flowDoms(d1, level, sps)%intCommVars(k)%var(i1, j1, k1) = &
               flowDoms(d1, level, sps)%intCommVars(k)%var(i1, j1, k1) + &
               flowDoms(d2, level, sps)%intCommVars(k)%var(i2, j2, k2)
          flowDoms(d2, level, sps)%intCommVars(k)%var(i2, j2, k2) = 0
       enddo

    enddo localCopy

    ! Complete the nonblocking receives in an arbitrary sequence and
    ! copy the variables from the buffer into the halo's.

    size = commPattern(level)%nProcSend
    completeSends: do i=1,commPattern(level)%nProcSend

       ! Complete any of the requests.

       call mpi_waitany(size, sendRequests, index, mpiStatus, ierr)
       call EChk(ierr,__FILE__,__LINE__)

       ! ! Copy the data just arrived in the halo's.

       ii = index
       jj = nVar*commPattern(level)%nsendCum(ii-1)
       !DIR$ NOVECTOR
       do j=1,commPattern(level)%nsend(ii)

          ! Store the block and the indices of the halo a bit easier.

          d2 = commPattern(level)%sendList(ii)%block(j)
          i2 = commPattern(level)%sendList(ii)%indices(j,1)+1
          j2 = commPattern(level)%sendList(ii)%indices(j,2)+1
          k2 = commPattern(level)%sendList(ii)%indices(j,3)+1

          ! Copy the conservative variables.
          !DIR$ NOVECTOR
          do k=1, nVar
             jj = jj + 1
             flowDoms(d2, level, sps)%intCommVars(k)%var(i2, j2, k2) = &
                  flowDoms(d2, level, sps)%intCommVars(k)%var(i2, j2, k2) + sendBufInt(jj)
          enddo
       enddo

    enddo completeSends

    ! Complete the nonblocking sends.

    size = commPattern(level)%nProcRecv
    do i=1,commPattern(level)%nProcRecv
       call mpi_waitany(size, recvRequests, index, mpiStatus, ierr)
       call EChk(ierr,__FILE__,__LINE__)
    enddo

    deallocate(recvBufInt, sendBufInt)

  end subroutine whalo1to1IntGeneric_b

  subroutine whalo1to1_d(level, start, end, commPressure,       &
       commVarGamma, commLamVis, commEddyVis, &
       commPattern, internal)
    !
    !       whalo1to1 exchanges the 1 to 1 internal halo's derivatives
    !
    use constants
    use communication, only : commType, internalCommType
    use inputTimeSpectral, only : nTimeIntervalsSpectral
    implicit none
    !
    !      Subroutine arguments.
    !
    integer(kind=intType), intent(in) :: level, start, end
    logical, intent(in) :: commPressure, commVarGamma
    logical, intent(in) :: commLamVis, commEddyVis

    type(commType), dimension(*), intent(in)         :: commPattern
    type(internalCommType), dimension(*), intent(in) :: internal

    integer(kind=intType) :: nVar, nn, k, sps

    spectralModes: do sps=1,nTimeIntervalsSpectral

       call setCommPointers(start, end, commPressure, commVarGamma, &
            commLamVis, commEddyVis, level, sps, .True., nVar, 0)

       if (nVar == 0) then
          return
       end if

       ! Run the generic exchange
       call wHalo1to1RealGeneric(nVar, level, sps, commPattern, internal)

    end do spectralModes

  end subroutine whalo1to1_d

  subroutine whalo1to1_b(level, start, end, commPressure,       &
       commVarGamma, commLamVis, commEddyVis, &
       commPattern, internal)
    !
    !       whalo1to1 exchanges the 1 to 1 internal halo's derivatives
    !
    use constants
    use communication, only : commType, internalCommType
    use inputTimeSpectral, only : nTimeIntervalsSpectral
    implicit none
    !
    !      Subroutine arguments.
    !
    integer(kind=intType), intent(in) :: level, start, end
    logical, intent(in) :: commPressure, commVarGamma
    logical, intent(in) :: commLamVis, commEddyVis

    type(commType), dimension(*), intent(in)         :: commPattern
    type(internalCommType), dimension(*), intent(in) :: internal

    integer(kind=intType) :: nVar, nn, k, sps

    spectralModes: do sps=1,nTimeIntervalsSpectral

       call setCommPointers(start, end, commPressure, commVarGamma, &
            commLamVis, commEddyVis, level, sps, .True., nVar, 0)

       if (nVar == 0) then
          return
       end if

       ! Run the generic exchange
       call wHalo1to1RealGeneric_b(nVar, level, sps, commPattern, internal)

    end do spectralModes

  end subroutine whalo1to1_b

  subroutine wOverset(level, start, end, commPressure,       &
       commVarGamma, commLamVis, commEddyVis, &
       commPattern, internal)
    !
    !       wOverset controls the communication between overset halos
    !       for the cell-centered variables by interpolating the solution
    !       from other blocks consistent with the chimera approach. A tri-
    !       linear interpolation is used as per the input from cgns. It
    !       is possible to send a range of variables and not the entire
    !       set, e.g. only the flow variables or only the turbulent
    !       variables. This is controlled by the arguments start, end,
    !       commPressure and commViscous. The exchange takes place for
    !       the given grid level.
    !
    use constants
    use communication, only : commType, internalCommType
    use inputTimeSpectral, only : nTimeIntervalsSpectral
    implicit none
    !
    !      Subroutine arguments.
    !
    integer(kind=intType), intent(in) :: level, start, end
    logical, intent(in) :: commPressure, commVarGamma
    logical, intent(in) :: commLamVis, commEddyVis

    type(commType), dimension(:, :), intent(in) :: commPattern
    type(internalCommType), dimension(:, :), intent(in) :: internal

    !      Local variables.
    integer(kind=intType) :: nVar, sps

    spectralModes: do sps=1,nTimeIntervalsSpectral

       call setCommPointers(start, end, commPressure, commVarGamma, &
            commLamVis, commEddyVis, level, sps, .False., nVar, 0)

       if (nVar == 0) then
          return
       end if

       ! Run the generic exchange
       call wOversetGeneric(nVar, level, sps, commPattern, internal)
    end do spectralModes

  end subroutine wOverset

  subroutine wOverset_d(level, start, end, commPressure,       &
       commVarGamma, commLamVis, commEddyVis, &
       commPattern, internal)
    !
    !       wOverset controls the communication between overset halos
    !       for the cell-centered variables by interpolating the solution
    !       from other blocks consistent with the chimera approach. A tri-
    !       linear interpolation is used as per the input from cgns. It
    !       is possible to send a range of variables and not the entire
    !       set, e.g. only the flow variables or only the turbulent
    !       variables. This is controlled by the arguments start, end,
    !       commPressure and commViscous. The exchange takes place for
    !       the given grid level.
    !
    use constants
    use communication, only : commType, internalCommType
    use inputTimeSpectral, only : nTimeIntervalsSpectral
    implicit none
    !
    !      Subroutine arguments.
    !
    integer(kind=intType), intent(in) :: level, start, end
    logical, intent(in) :: commPressure, commVarGamma
    logical, intent(in) :: commLamVis, commEddyVis

    type(commType), dimension(:, :), intent(in) :: commPattern
    type(internalCommType), dimension(:, :), intent(in) :: internal
    integer(kind=intType) :: nVar, sps, offset

    spectralModes: do sps=1,nTimeIntervalsSpectral

       ! this one is tricker: We have to set BOTH the real values and
       ! the the derivative values. Set the derivative values first:
       call setCommPointers(start, end, commPressure, commVarGamma, &
            commLamVis, commEddyVis, level, sps, .True., nVar, 0)

       ! And then the original real values
       offset = nVar
       call setCommPointers(start, end, commPressure, commVarGamma, &
            commLamVis, commEddyVis, level, sps, .False., nVar, offset)

       if (nVar == 0) then
          return
       end if

       ! Run the generic exchange
       call wOversetGeneric_d(nVar, level, sps, commPattern, internal)
    end do spectralModes
  end subroutine wOverset_d

  subroutine wOverset_b(level, start, end, commPressure,       &
       commVarGamma, commLamVis, commEddyVis, &
       commPattern, internal)
    !
    !       wOverset_b performs the *TRANSPOSE* operation of wOveset
    !       It is used for adjoint/reverse mode residual evaluations.
    !      * See wOverset  for more information.
    !
    use constants
    use communication, only : commType, internalCommType
    use inputTimeSpectral, only : nTimeIntervalsSpectral
    implicit none
    !
    !      Subroutine arguments.
    !
    integer(kind=intType), intent(in) :: level, start, end
    logical, intent(in) :: commPressure, commVarGamma
    logical, intent(in) :: commLamVis, commEddyVis

    type(commType), dimension(:, :), intent(in) :: commPattern
    type(internalCommType), dimension(:, :), intent(in) :: internal
    integer(kind=intType) :: nVar, sps, offset

    spectralModes: do sps=1,nTimeIntervalsSpectral

       ! this one is tricker: We have to set BOTH the real values and
       ! the the derivative values. Set the derivative values first:
       call setCommPointers(start, end, commPressure, commVarGamma, &
            commLamVis, commEddyVis, level, sps, .True., nVar, 0)

       ! And then the original real values
       offset = nVar
       call setCommPointers(start, end, commPressure, commVarGamma, &
            commLamVis, commEddyVis, level, sps, .False., nVar, offset)

       if (nVar == 0) then
          return
       end if

       ! Run the generic exchange
       call wOversetGeneric_b(nVar, level, sps, commPattern, internal)
    end do spectralModes

  end subroutine wOverset_b

  subroutine wOversetGeneric(nVar, level, sps, commPattern, Internal)
    !
    !       wOverset is the generic halo exhcnage code for the
    !       overset halos.
    !
    use constants
    use block, only : flowDoms
    use communication
    use utils, only : EChk
    implicit none
    !
    !      Subroutine arguments.
    !
    integer(kind=intType), intent(in) :: level, sps

    type(commType), dimension(:, :), intent(in)         :: commPattern
    type(internalCommType), dimension(:, :), intent(in) :: internal
    !
    !      Local variables.
    !
    integer :: size, procId, ierr, index
    integer, dimension(mpi_status_size) :: mpiStatus

    integer(kind=intType) :: nVar
    integer(kind=intType) :: i, j, k, ii, jj
    integer(kind=intType) :: d1, i1, j1, k1, d2, i2, j2, k2
    real(kind=realType), dimension(:), pointer :: weight

    ! Send the variables. The data is first copied into
    ! the send buffer after which the buffer is sent asap.

    ii = 1
    sends: do i=1,commPattern(level, sps)%nProcSend

       ! Store the processor id and the size of the message
       ! a bit easier.

       procID = commPattern(level, sps)%sendProc(i)
       size    = nVar*commPattern(level, sps)%nsend(i)

       ! Copy the data in the correct part of the send buffer.

       jj = ii
       !DIR$ NOVECTOR
       do j=1,commPattern(level, sps)%nsend(i)

          ! Store the block id and the indices of the donor
          ! a bit easier.

          d1 = commPattern(level, sps)%sendList(i)%block(j)
          i1 = commPattern(level, sps)%sendList(i)%indices(j,1)+1
          j1 = commPattern(level, sps)%sendList(i)%indices(j,2)+1
          k1 = commPattern(level, sps)%sendList(i)%indices(j,3)+1
          weight => commPattern(level, sps)%sendList(i)%interp(j, :)

          ! Copy the given range of the working variables for
          ! this cell in the buffer. Update the counter jj.
          !DIR$ NOVECTOR
          do k=1, nvar
             sendBuffer(jj) = &
                  weight(1)*flowDoms(d1,level,sps)%realCommVars(k)%var(i1  , j1,   k1  ) + &
                  weight(2)*flowDoms(d1,level,sps)%realCommVars(k)%var(i1+1, j1,   k1  ) + &
                  weight(3)*flowDoms(d1,level,sps)%realCommVars(k)%var(i1,   j1+1, k1  ) + &
                  weight(4)*flowDoms(d1,level,sps)%realCommVars(k)%var(i1+1, j1+1, k1  ) + &
                  weight(5)*flowDoms(d1,level,sps)%realCommVars(k)%var(i1  , j1,   k1+1) + &
                  weight(6)*flowDoms(d1,level,sps)%realCommVars(k)%var(i1+1, j1,   k1+1) + &
                  weight(7)*flowDoms(d1,level,sps)%realCommVars(k)%var(i1,   j1+1, k1+1) + &
                  weight(8)*flowDoms(d1,level,sps)%realCommVars(k)%var(i1+1, j1+1, k1+1)
             jj = jj + 1
          end do
       enddo

       ! Send the data.

       call mpi_isend(sendBuffer(ii), size, adflow_real, procId,  &
            procId, ADflow_comm_world, sendRequests(i), &
            ierr)
       call EChk(ierr,__FILE__,__LINE__)

       ! Set ii to jj for the next processor.

       ii = jj

    enddo sends

    ! Post the nonblocking receives.

    ii = 1
    receives: do i=1,commPattern(level, sps)%nProcRecv

       ! Store the processor id and the size of the message
       ! a bit easier.

       procID = commPattern(level,sps)%recvProc(i)
       size    = nVar*commPattern(level,sps)%nrecv(i)

       ! Post the receive.

       call mpi_irecv(recvBuffer(ii), size, adflow_real, procId, &
            myId, ADflow_comm_world, recvRequests(i), ierr)
       call EChk(ierr,__FILE__,__LINE__)

       ! And update ii.

       ii = ii + size

    enddo receives

    ! Do the local interpolation.
    !DIR$ NOVECTOR
    localInterp: do i=1,internal(level, sps)%ncopy

       ! Store the block and the indices of the donor a bit easier.

       d1 = internal(level,sps)%donorBlock(i)
       i1 = internal(level,sps)%donorIndices(i, 1)+1
       j1 = internal(level,sps)%donorIndices(i, 2)+1
       k1 = internal(level,sps)%donorIndices(i, 3)+1

       weight => internal(level,sps)%donorInterp(i, :)

       ! Idem for the halo's.

       d2 = internal(level,sps)%haloBlock(i)
       i2 = internal(level,sps)%haloIndices(i, 1)+1
       j2 = internal(level,sps)%haloIndices(i, 2)+1
       k2 = internal(level,sps)%haloIndices(i, 3)+1

       ! Copy the given range of working variables.
       !DIR$ NOVECTOR
       do k=1, nVar
          flowDoms(d2, level, sps)%realCommVars(k)%var(i2, j2, k2) = &
               weight(1)*flowDoms(d1, level, sps)%realCommVars(k)%var(i1,   j1,   k1  ) + &
               weight(2)*flowDoms(d1, level, sps)%realCommVars(k)%var(i1+1, j1,   k1  ) + &
               weight(3)*flowDoms(d1, level, sps)%realCommVars(k)%var(i1,   j1+1, k1  ) + &
               weight(4)*flowDoms(d1, level, sps)%realCommVars(k)%var(i1+1, j1+1, k1  ) + &
               weight(5)*flowDoms(d1, level, sps)%realCommVars(k)%var(i1,   j1,   k1+1) + &
               weight(6)*flowDoms(d1, level, sps)%realCommVars(k)%var(i1+1, j1,   k1+1) + &
               weight(7)*flowDoms(d1, level, sps)%realCommVars(k)%var(i1,   j1+1, k1+1) + &
               weight(8)*flowDoms(d1, level, sps)%realCommVars(k)%var(i1+1, j1+1, k1+1)
       end do
    enddo localInterp

    ! Complete the nonblocking receives in an arbitrary sequence and
    ! copy the variables from the buffer into the halo's.

    size = commPattern(level, sps)%nProcRecv
    completeRecvs: do i=1,commPattern(level, sps)%nProcRecv

       ! Complete any of the requests.

       call mpi_waitany(size, recvRequests, index, mpiStatus, ierr)
       call EChk(ierr,__FILE__,__LINE__)

       ! Copy the data just arrived in the halo's.

       ii = index
       jj = nVar*commPattern(level,sps)%nrecvCum(ii-1)
       !DIR$ NOVECTOR
       do j=1,commPattern(level,sps)%nrecv(ii)

          ! Store the block and the indices of the halo a bit easier.

          d2 = commPattern(level,sps)%recvList(ii)%block(j)
          i2 = commPattern(level,sps)%recvList(ii)%indices(j,1)+1
          j2 = commPattern(level,sps)%recvList(ii)%indices(j,2)+1
          k2 = commPattern(level,sps)%recvList(ii)%indices(j,3)+1
          !DIR$ NOVECTOR
          do k=1, nVar
             jj = jj + 1
             flowDoms(d2,level,sps)%realCommVars(k)%var(i2,j2,k2) = recvBuffer(jj)
          enddo
       enddo
    end do completeRecvs

    ! Complete the nonblocking sends.

    size = commPattern(level,sps)%nProcSend
    do i=1,commPattern(level,sps)%nProcSend
       call mpi_waitany(size, sendRequests, index, mpiStatus, ierr)
       call EChk(ierr,__FILE__,__LINE__)
    enddo

  end subroutine wOversetGeneric

  subroutine wOversetGeneric_d(nVar, level, sps, commPattern, Internal)
    !
    !       wOverset_d is the generic halo forward mode linearized
    !       code for overset halos.
    !
    use constants
    use block, only : flowDoms
    use communication
    use utils, only : EChk
    implicit none
    !
    !      Subroutine arguments.
    !
    integer(kind=intType), intent(in) :: level, sps

    type(commType), dimension(:, :), intent(in)         :: commPattern
    type(internalCommType), dimension(:, :), intent(in) :: internal
    !
    !      Local variables.
    !
    integer :: size, procId, ierr, index
    integer, dimension(mpi_status_size) :: mpiStatus

    integer(kind=intType) :: nVar
    integer(kind=intType) :: i, j, k, ii, jj
    integer(kind=intType) :: d1, i1, j1, k1, d2, i2, j2, k2
    real(kind=realType), dimension(:), pointer :: weight, weightd

    ! Send the variables. The data is first copied into
    ! the send buffer after which the buffer is sent asap.

    ii = 1
    sends: do i=1,commPattern(level, sps)%nProcSend

       ! Store the processor id and the size of the message
       ! a bit easier.

       procID = commPattern(level, sps)%sendProc(i)
       size    = nVar*commPattern(level, sps)%nsend(i)

       ! Copy the data in the correct part of the send buffer.

       jj = ii
       !DIR$ NOVECTOR
       do j=1,commPattern(level, sps)%nsend(i)

          ! Store the block id and the indices of the donor
          ! a bit easier.

          d1 = commPattern(level, sps)%sendList(i)%block(j)
          i1 = commPattern(level, sps)%sendList(i)%indices(j,1)+1
          j1 = commPattern(level, sps)%sendList(i)%indices(j,2)+1
          k1 = commPattern(level, sps)%sendList(i)%indices(j,3)+1
          weight => commPattern(level, sps)%sendList(i)%interp(j, :)
          weightd => commPattern(level, sps)%sendList(i)%interpd(j, :)

          ! Copy the given range of the working variables for
          ! this cell in the buffer. Update the counter jj.
          !DIR$ NOVECTOR
          do k=1, nvar
             sendBuffer(jj) = &
                  weight(1)*flowDoms(d1,level,sps)%realCommVars(k)%var(i1  , j1,   k1  ) + &
                  weight(2)*flowDoms(d1,level,sps)%realCommVars(k)%var(i1+1, j1,   k1  ) + &
                  weight(3)*flowDoms(d1,level,sps)%realCommVars(k)%var(i1,   j1+1, k1  ) + &
                  weight(4)*flowDoms(d1,level,sps)%realCommVars(k)%var(i1+1, j1+1, k1  ) + &
                  weight(5)*flowDoms(d1,level,sps)%realCommVars(k)%var(i1  , j1,   k1+1) + &
                  weight(6)*flowDoms(d1,level,sps)%realCommVars(k)%var(i1+1, j1,   k1+1) + &
                  weight(7)*flowDoms(d1,level,sps)%realCommVars(k)%var(i1,   j1+1, k1+1) + &
                  weight(8)*flowDoms(d1,level,sps)%realCommVars(k)%var(i1+1, j1+1, k1+1) + &
                  weightd(1)*flowDoms(d1,level,sps)%realCommVars(k+nVar)%var(i1  , j1,   k1  ) + &
                  weightd(2)*flowDoms(d1,level,sps)%realCommVars(k+nVar)%var(i1+1, j1,   k1  ) + &
                  weightd(3)*flowDoms(d1,level,sps)%realCommVars(k+nVar)%var(i1,   j1+1, k1  ) + &
                  weightd(4)*flowDoms(d1,level,sps)%realCommVars(k+nVar)%var(i1+1, j1+1, k1  ) + &
                  weightd(5)*flowDoms(d1,level,sps)%realCommVars(k+nVar)%var(i1  , j1,   k1+1) + &
                  weightd(6)*flowDoms(d1,level,sps)%realCommVars(k+nVar)%var(i1+1, j1,   k1+1) + &
                  weightd(7)*flowDoms(d1,level,sps)%realCommVars(k+nVar)%var(i1,   j1+1, k1+1) + &
                  weightd(8)*flowDoms(d1,level,sps)%realCommVars(k+nVar)%var(i1+1, j1+1, k1+1)


             jj = jj + 1
          end do
       enddo

       ! Send the data.

       call mpi_isend(sendBuffer(ii), size, adflow_real, procId,  &
            procId, ADflow_comm_world, sendRequests(i), &
            ierr)
       call EChk(ierr,__FILE__,__LINE__)

       ! Set ii to jj for the next processor.

       ii = jj

    enddo sends

    ! Post the nonblocking receives.

    ii = 1
    receives: do i=1,commPattern(level, sps)%nProcRecv

       ! Store the processor id and the size of the message
       ! a bit easier.

       procID = commPattern(level,sps)%recvProc(i)
       size    = nVar*commPattern(level,sps)%nrecv(i)

       ! Post the receive.

       call mpi_irecv(recvBuffer(ii), size, adflow_real, procId, &
            myId, ADflow_comm_world, recvRequests(i), ierr)
       call EChk(ierr,__FILE__,__LINE__)

       ! And update ii.

       ii = ii + size

    enddo receives

    ! Do the local interpolation.
    !DIR$ NOVECTOR
    localInterp: do i=1,internal(level, sps)%ncopy

       ! Store the block and the indices of the donor a bit easier.

       d1 = internal(level,sps)%donorBlock(i)
       i1 = internal(level,sps)%donorIndices(i, 1)+1
       j1 = internal(level,sps)%donorIndices(i, 2)+1
       k1 = internal(level,sps)%donorIndices(i, 3)+1

       weight => internal(level,sps)%donorInterp(i, :)
       weightd => internal(level, sps)%donorInterpd(i, :)

       ! Idem for the halo's.

       d2 = internal(level,sps)%haloBlock(i)
       i2 = internal(level,sps)%haloIndices(i, 1)+1
       j2 = internal(level,sps)%haloIndices(i, 2)+1
       k2 = internal(level,sps)%haloIndices(i, 3)+1

       ! Copy the given range of working variables.
       !DIR$ NOVECTOR
       do k=1, nVar
          flowDoms(d2, level, sps)%realCommVars(k)%var(i2, j2, k2) = &
               weight(1)*flowDoms(d1,level,sps)%realCommVars(k)%var(i1  , j1,   k1  ) + &
               weight(2)*flowDoms(d1,level,sps)%realCommVars(k)%var(i1+1, j1,   k1  ) + &
               weight(3)*flowDoms(d1,level,sps)%realCommVars(k)%var(i1,   j1+1, k1  ) + &
               weight(4)*flowDoms(d1,level,sps)%realCommVars(k)%var(i1+1, j1+1, k1  ) + &
               weight(5)*flowDoms(d1,level,sps)%realCommVars(k)%var(i1  , j1,   k1+1) + &
               weight(6)*flowDoms(d1,level,sps)%realCommVars(k)%var(i1+1, j1,   k1+1) + &
               weight(7)*flowDoms(d1,level,sps)%realCommVars(k)%var(i1,   j1+1, k1+1) + &
               weight(8)*flowDoms(d1,level,sps)%realCommVars(k)%var(i1+1, j1+1, k1+1) + &
               weightd(1)*flowDoms(d1,level,sps)%realCommVars(k+nVar)%var(i1  , j1,   k1  ) + &
               weightd(2)*flowDoms(d1,level,sps)%realCommVars(k+nVar)%var(i1+1, j1,   k1  ) + &
               weightd(3)*flowDoms(d1,level,sps)%realCommVars(k+nVar)%var(i1,   j1+1, k1  ) + &
               weightd(4)*flowDoms(d1,level,sps)%realCommVars(k+nVar)%var(i1+1, j1+1, k1  ) + &
               weightd(5)*flowDoms(d1,level,sps)%realCommVars(k+nVar)%var(i1  , j1,   k1+1) + &
               weightd(6)*flowDoms(d1,level,sps)%realCommVars(k+nVar)%var(i1+1, j1,   k1+1) + &
               weightd(7)*flowDoms(d1,level,sps)%realCommVars(k+nVar)%var(i1,   j1+1, k1+1) + &
               weightd(8)*flowDoms(d1,level,sps)%realCommVars(k+nVar)%var(i1+1, j1+1, k1+1)

       end do
    enddo localInterp

    ! Complete the nonblocking receives in an arbitrary sequence and
    ! copy the variables from the buffer into the halo's.

    size = commPattern(level, sps)%nProcRecv
    completeRecvs: do i=1,commPattern(level, sps)%nProcRecv

       ! Complete any of the requests.

       call mpi_waitany(size, recvRequests, index, mpiStatus, ierr)
       call EChk(ierr,__FILE__,__LINE__)

       ! Copy the data just arrived in the halo's.

       ii = index
       jj = nVar*commPattern(level,sps)%nrecvCum(ii-1)
       !DIR$ NOVECTOR
       do j=1,commPattern(level,sps)%nrecv(ii)

          ! Store the block and the indices of the halo a bit easier.

          d2 = commPattern(level,sps)%recvList(ii)%block(j)
          i2 = commPattern(level,sps)%recvList(ii)%indices(j,1)+1
          j2 = commPattern(level,sps)%recvList(ii)%indices(j,2)+1
          k2 = commPattern(level,sps)%recvList(ii)%indices(j,3)+1

          do k=1, nVar
             jj = jj + 1
             flowDoms(d2,level,sps)%realCommVars(k)%var(i2,j2,k2) = recvBuffer(jj)
          enddo
       enddo
    end do completeRecvs

    ! Complete the nonblocking sends.

    size = commPattern(level,sps)%nProcSend
    do i=1,commPattern(level,sps)%nProcSend
       call mpi_waitany(size, sendRequests, index, mpiStatus, ierr)
       call EChk(ierr,__FILE__,__LINE__)
   enddo

  end subroutine wOversetGeneric_d

  subroutine wOversetGeneric_b(nVar, level, sps, commPattern, Internal)
    !
    !       wOversetGeneric_b is the generic reverse mode linearized
    !       code for overset halos.
    !
    use constants
    use block, only : flowDoms
    use communication
    use utils, only : EChk
    implicit none
    !
    !      Subroutine arguments.
    !
    integer(kind=intType), intent(in) :: level, sps

    type(commType), dimension(:, :), intent(in)         :: commPattern
    type(internalCommType), dimension(:, :), intent(in) :: internal
    !
    !      Local variables.
    !
    integer :: size, procId, ierr, index
    integer, dimension(mpi_status_size) :: mpiStatus

    integer(kind=intType) :: nVar
    integer(kind=intType) :: i, j, k, ii, jj, kk
    integer(kind=intType) :: d1, i1, j1, k1, d2, i2, j2, k2, iii, jjj, kkk
    real(kind=realType), dimension(:), pointer :: weight, weightd
    real(kind=realType) :: vard
    ! Gather up the seeds into the *recv* buffer. Note we loop over
    ! nProcRECV here! After the buffer is assembled it is send off.

    jj = 1
    ii = 1
    recvs: do i=1, commPattern(level, sps)%nProcRecv

       ! Store the processor id and the size of the message
       ! a bit easier.

       procID = commPattern(level, sps)%recvProc(i)
       size    = nVar*commPattern(level, sps)%nrecv(i)

       ! Copy the data into the buffer
       !DIR$ NOVECTOR
       do j=1,commPattern(level, sps)%nrecv(i)

          ! Store the block and the indices to make code a bit easier to read

          d2 = commPattern(level, sps)%recvList(i)%block(j)
          i2 = commPattern(level, sps)%recvList(i)%indices(j,1)+1
          j2 = commPattern(level, sps)%recvList(i)%indices(j,2)+1
          k2 = commPattern(level, sps)%recvList(i)%indices(j,3)+1
          !DIR$ NOVECTOR
          do k=1, nVar
             recvBuffer(jj) = flowDoms(d2, level, sps)%realCommVars(k)%var(i2, j2, k2)
             jj = jj + 1
             flowDoms(d2, level, sps)%realCommVars(k)%var(i2, j2, k2) = zero
          enddo
       enddo

       ! Send the data.
       call mpi_isend(recvBuffer(ii), size, adflow_real, procID,  &
            procID, ADflow_comm_world, recvRequests(i), &
            ierr)
       call EChk(ierr,__FILE__,__LINE__)

       ! Set ii to jj for the next processor.

       ii = jj

    end do recvs

    ! Post the nonblocking receives.

    ii = 1
    sends: do i=1,commPattern(level, sps)%nProcSend

       ! Store the processor id and the size of the message
       ! a bit easier.

       procID = commPattern(level, sps)%sendProc(i)
       size    = nVar*commPattern(level, sps)%nsend(i)

       ! Post the receive.

       call mpi_irecv(sendBuffer(ii), size, adflow_real, procId, &
            myId, ADflow_comm_world, sendRequests(i), ierr)
       call EChk(ierr,__FILE__,__LINE__)

       ! And update ii.

       ii = ii + size

    enddo sends

    ! Do the local interpolation.
    !DIR$ NOVECTOR
    localInterp: do i=1,internal(level, sps)%ncopy

       ! Store the block and the indices of the donor a bit easier.

       d1 = internal(level, sps)%donorBlock(i)
       i1 = internal(level, sps)%donorIndices(i, 1)+1
       j1 = internal(level, sps)%donorIndices(i, 2)+1
       k1 = internal(level, sps)%donorIndices(i, 3)+1

       weight => internal(level, sps)%donorInterp(i, :)
       weightd => internal(level, sps)%donorInterpd(i, :)

       ! Idem for the halo's.

       d2 = internal(level, sps)%haloBlock(i)
       i2 = internal(level, sps)%haloIndices(i, 1)+1
       j2 = internal(level, sps)%haloIndices(i, 2)+1
       k2 = internal(level, sps)%haloIndices(i, 3)+1

       ! Sum into the '1' values from the '2' values accouting for the weights
       !DIR$ NOVECTOR
       do k=1, nVar
          vard = flowDoms(d2, level, sps)%realCommVars(k)%var(i2, j2, k2)
          kk = 0
          do kkk=k1,k1+1
             do jjj=j1,j1+1
                do iii=i1,i1+1
                   kk = kk + 1
                   flowDoms(d1, level, sps)%realCommVars(k)%var(iii, jjj, kkk) = &
                        flowDoms(d1, level, sps)%realCommVars(k)%var(iii, jjj, kkk) + &
                        weight(kk)*vard

                   weightd(kk) = weightd(kk) + &
                        flowDoms(d1, level, sps)%realCommVars(k+nVar)%var(iii,jjj,kkk)*vard

                end do
             end do
          end do
          flowDoms(d2, level, sps)%realCommVars(k)%var(i2, j2, k2) = zero
       end do
    enddo localInterp

    ! Complete the nonblocking receives in an arbitrary sequence and
    ! copy the variables from the buffer into the halo's.

    size = commPattern(level, sps)%nProcSend
    completeSends: do i=1,commPattern(level, sps)%nProcSend

       ! Complete any of the requests.

       call mpi_waitany(size, sendRequests, index, mpiStatus, ierr)
       call EChk(ierr,__FILE__,__LINE__)

       ! Copy the data just arrived in the halo's.

       ii = index

       jj = nVar*commPattern(level, sps)%nsendCum(ii-1)
       !DIR$ NOVECTOR
       do j=1,commPattern(level, sps)%nsend(ii)

          ! Store the block and the indices of the halo a bit easier.

          d2 = commPattern(level, sps)%sendList(ii)%block(j)
          i2 = commPattern(level, sps)%sendList(ii)%indices(j,1)+1
          j2 = commPattern(level, sps)%sendList(ii)%indices(j,2)+1
          k2 = commPattern(level, sps)%sendList(ii)%indices(j,3)+1

          weight => commPattern(level, sps)%sendList(ii)%interp(j, :)
          weightd => commPattern(level, sps)%sendList(ii)%interpd(j, :)
          !DIR$ NOVECTOR
          do k=1, nVar
             jj =jj + 1
             vard = sendBuffer(jj)
             kk = 0
             do kkk=k2,k2+1
                do jjj=j2,j2+1
                   do iii=i2,i2+1

                      kk = kk + 1
                      flowDoms(d2, level, sps)%realCommVars(k)%var(iii, jjj, kkk) = &
                           flowDoms(d2, level, sps)%realCommVars(k)%var(iii, jjj, kkk) + &
                           weight(kk)*vard

                      weightd(kk) = weightd(kk) + &
                           flowDoms(d2, level, sps)%realCommVars(k+nVar)%var(iii,jjj,kkk)*vard

                   end do
                end do
             end do
          end do
       enddo

    enddo completeSends

    ! Complete the nonblocking sends.

    size = commPattern(level, sps)%nProcRecv
    do i=1,commPattern(level, sps)%nProcRecv
       call mpi_waitany(size, recvRequests, index, mpiStatus, ierr)
       call EChk(ierr,__FILE__,__LINE__)
    enddo

  end subroutine wOversetGeneric_b


#ifndef USE_COMPLEX
  subroutine whalo2_b(level, start, end, commPressure, commGamma, &
       commViscous)
    !
    !       whalo2_b exchanges all the 2nd level internal halo's for the
    !       cell centered variables IN REVERSE MODE
    !
    use constants
    use blockPointers
    use communication
    use flowVarRefState
    use inputPhysics
    use inputTimeSpectral
    use iteration
    use flowUtils_b, only : computeETotBlock_b
    use utils, only : setPointers_b, getCorrectForK
    implicit none
    !
    !      Subroutine arguments.
    !
    integer(kind=intType), intent(in) :: level, start, end
    logical, intent(in) :: commPressure, commGamma, commViscous
    !
    !      Local variables.
    !
    integer(kind=intType) :: nn, mm, ll
    integer(kind=intType) :: iBeg, iEnd, jBeg, jEnd, kBeg, kEnd

    logical :: correctForK, commLamVis, commEddyVis, commVarGamma

    ! Set the logicals whether or not to communicate the viscosities.

    commLamVis = .false.
    if(viscous .and. commViscous) commLamVis = .true.

    commEddyVis = .false.
    if(eddyModel .and. commViscous) commEddyVis = .true.

    ! Set the logical whether or not to communicate gamma.

    commVarGamma = .false.
    if(commGamma .and. (cpModel == cpTempCurveFits)) &
         commVarGamma = .true.

    bothPAndE: if(commPressure .and. start <= irhoE .and. &
         end >= irhoE) then

       ! First determine whether or not the total energy must be
       ! corrected for the presence of the turbulent kinetic energy.

       correctForK = getCorrectForK()

       domains: do nn=1,nDom

          ! Treat the overset blocks. Since we don't have the logic
          ! setup here correctly to only update the overset cells,
          ! just do the whole block, for every block
          do ll=1, nTimeIntervalsSpectral
             call setPointers_b(nn, level, ll)
             call computeETotBlock_b(2, il, 2, jl, 2, kl, correctForK)
          end do
       enddo domains
    endif bothPAndE

    call wOverset_b(level, start, end, commPressure, commVarGamma, &
         commLamVis, commEddyVis, commPatternOverset, internalOverset)

    ! Exchange the 1 to 1 matching 2nd level cell halo's.
    call whalo1to1_b(level, start, end, commPressure, commVarGamma, &
         commLamVis, commEddyVis, commPatternCell_2nd,  &
         internalCell_2nd)

    ! NOTE: Only the 1to1 halo exchange and overset is done. whalosliding,
    ! whalomixing, orphanAverage and PandE corrections
    ! calculation are NOT implementent.

  end subroutine whalo2_b

  subroutine whalo2_d(level, start, end, commPressure, commGamma, &
       commViscous)
    !
    !       whalo2_b exchanges all the 2nd level internal halo's for the
    !       cell centered variables IN FORWARD MODE
    !
    use constants
    use blockPointers
    use communication
    use flowVarRefState
    use inputPhysics
    use inputTimeSpectral
    use iteration
    use flowUtils_d, only : computeETotBlock_d
    use utils, only : setPointers_d, getCorrectForK
    implicit none
    !
    !      Subroutine arguments.
    !
    integer(kind=intType), intent(in) :: level, start, end
    logical, intent(in) :: commPressure, commGamma, commViscous
    !
    !      Local variables.
    !
    integer(kind=intType) :: nn, mm, ll
    integer(kind=intType) :: iBeg, iEnd, jBeg, jEnd, kBeg, kEnd

    logical :: correctForK, commLamVis, commEddyVis, commVarGamma

    ! Set the logicals whether or not to communicate the viscosities.

    commLamVis = .false.
    if(viscous .and. commViscous) commLamVis = .true.

    commEddyVis = .false.
    if(eddyModel .and. commViscous) commEddyVis = .true.

    ! Set the logical whether or not to communicate gamma.

    commVarGamma = .false.
    if(commGamma .and. (cpModel == cpTempCurveFits)) &
         commVarGamma = .true.

    ! Exchange the 1 to 1 matching 2nd level cell halo's.

    call whalo1to1_d(level, start, end, commPressure, commVarGamma, &
         commLamVis, commEddyVis, commPatternCell_2nd,  &
         internalCell_2nd)

    ! Exchange the overset cells
    call wOverset_d(level, start, end, commPressure, commVarGamma, &
         commLamVis, commEddyVis, commPatternOverset, internalOverset)

    ! NOTE: Only the 1to1 halo and wOverset exchange is done. whalosliding,
    ! whalomixing, orphanAverage and PandE corrections
    ! calculation are NOT implementent.

    ! If both the pressure and the total energy has been communicated
    ! compute the energy again. The reason is that both values are
    ! interpolated and consequently the values are not consistent.
    ! The energy depends quadratically on the velocity.

    bothPAndE: if(commPressure .and. start <= irhoE .and. &
         end >= irhoE) then

       ! First determine whether or not the total energy must be
       ! corrected for the presence of the turbulent kinetic energy.

       correctForK = getCorrectForK()

       domains: do nn=1,nDom

          ! Treat the overset blocks. Since we don't have the logic
          ! setup here correctly to only update the overset cells,
          ! just do the whole block, for every block
          do ll=1, nTimeIntervalsSpectral
             call setPointers_d(nn, level, ll)
             call computeETotBlock_d(2, il, 2, jl, 2, kl, correctForK)
          end do
       enddo domains

    endif bothPAndE
  end subroutine whalo2_d
#endif


  subroutine resHalo1(level, start, end)
    !
    !       resHalo1 determines the residuals in the 1st layer of halo
    !       cells by applying both the boundary conditions and the
    !       exchange. The halo values are needed for post processing
    !       reasons.
    !
    use constants
    use blockPointers
    use communication
    use inputTimeSpectral
    use utils, only : setPointers, EChk
    implicit none
    !
    !      Subroutine arguments.
    !
    integer(kind=intType) :: level, start, end
    !
    !      Local variables.
    !
    integer :: size, procID, ierr, index
    integer, dimension(mpi_status_size) :: mpiStatus

    integer(kind=intType) :: nVar, sps
    integer(kind=intType) :: ii, jj, mm, nn, i, j, k, l
    integer(kind=intType) :: dd1, ii1, jj1, kk1, dd2, ii2, jj2, kk2

    real(kind=realType), pointer, dimension(:,:,:) :: ddw1, ddw2

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
             !DIR$ NOVECTOR
             do j=BCData(nn)%jcBeg, BCData(nn)%jcEnd
                !DIR$ NOVECTOR
                do i=BCData(nn)%icBeg, BCData(nn)%icEnd
                   !DIR$ NOVECTOR
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
          !DIR$ NOVECTOR
          do j=1,commPatternCell_1st(level)%nsend(i)

             ! Store the block id and the indices of the donor
             ! a bit easier.

             dd1 = commPatternCell_1st(level)%sendList(i)%block(j)
             ii1 = commPatternCell_1st(level)%sendList(i)%indices(j,1)
             jj1 = commPatternCell_1st(level)%sendList(i)%indices(j,2)
             kk1 = commPatternCell_1st(level)%sendList(i)%indices(j,3)

             ! Copy the given range of the residuals for this cell
             ! in the buffer. Update the counter jj accordingly.
             !DIR$ NOVECTOR
             do k=start,end
                sendBuffer(jj) = flowDoms(dd1,level,sps)%dw(ii1,jj1,kk1,k)
                jj = jj + 1
             enddo

          enddo

          ! Send the data.

          call mpi_isend(sendBuffer(ii), size, adflow_real, procID,  &
               procID, ADflow_comm_world, sendRequests(i), &
               ierr)
          call EChk(ierr,__FILE__,__LINE__)

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

          call mpi_irecv(recvBuffer(ii), size, adflow_real, procID, &
               myID, ADflow_comm_world, recvRequests(i), ierr)
          call EChk(ierr,__FILE__,__LINE__)

          ! And update ii.

          ii = ii + size

       enddo receives

       ! Copy the local data.
       !DIR$ NOVECTOR
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
          !DIR$ NOVECTOR
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

          call mpi_waitany(size, recvRequests, index, mpiStatus, ierr)
          call EChk(ierr,__FILE__,__LINE__)

          ! Copy the data just arrived in the halo's.

          ii = index
          jj = nVar*commPatternCell_1st(level)%nrecvCum(ii-1) +1
          !DIR$ NOVECTOR
          do j=1,commPatternCell_1st(level)%nrecv(ii)

             ! Store the block and the indices of the halo a bit easier.

             dd2 = commPatternCell_1st(level)%recvList(ii)%block(j)
             ii2 = commPatternCell_1st(level)%recvList(ii)%indices(j,1)
             jj2 = commPatternCell_1st(level)%recvList(ii)%indices(j,2)
             kk2 = commPatternCell_1st(level)%recvList(ii)%indices(j,3)

             ! Copy the residuals.
             !DIR$ NOVECTOR
             do k=start,end
                flowDoms(dd2,level,sps)%dw(ii2,jj2,kk2,k) = recvBuffer(jj)
                jj = jj + 1
             enddo

          enddo

       enddo completeRecvs

       ! Complete the nonblocking sends.

       size = commPatternCell_1st(level)%nProcSend
       do i=1,commPatternCell_1st(level)%nProcSend
          call mpi_waitany(size, sendRequests, index, mpiStatus, ierr)
          call EChk(ierr,__FILE__,__LINE__)
    enddo

    enddo spectralLoop

  end subroutine resHalo1

  subroutine exchangeCoor(level)
    !
    !       ExchangeCoor exchanges the coordinates of the given grid
    !       level.
    !
    use block
    use communication
    use inputTimeSpectral
    use utils, only : EChk
    !
    !      Subroutine arguments.
    !
    integer(kind=intType), intent(in) :: level
    !
    !      Local variables.
    !
    integer :: size, procID, ierr, index
    integer, dimension(mpi_status_size) :: mpiStatus

    integer(kind=intType) :: i, j, ii, jj, mm
    integer(kind=intType) :: d1, i1, j1, k1, d2, i2, j2, k2

    ! Loop over the number of spectral solutions.

    spectralLoop: do mm=1,nTimeIntervalsSpectral

       ! Send the coordinates i have to send. The data is first copied
       ! into the send buffer and this buffer is sent.

       ii = 1
       sends: do i=1,commPatternNode_1st(level)%nProcSend

          ! Store the processor id and the size of the message
          ! a bit easier.

          procID = commPatternNode_1st(level)%sendProc(i)
          size   = 3*commPatternNode_1st(level)%nSend(i)

          ! Copy the data in the correct part of the send buffer.

          jj = ii
          !DIR$ NOVECTOR
          do j=1,commPatternNode_1st(level)%nSend(i)

             ! Store the block id and the indices of the donor
             ! a bit easier.

             d1 = commPatternNode_1st(level)%sendList(i)%block(j)
             i1 = commPatternNode_1st(level)%sendList(i)%indices(j,1)
             j1 = commPatternNode_1st(level)%sendList(i)%indices(j,2)
             k1 = commPatternNode_1st(level)%sendList(i)%indices(j,3)

             ! Copy the coordinates of this point in the buffer.
             ! Update the counter jj accordingly.

             sendBuffer(jj)   = flowDoms(d1,level,mm)%x(i1,j1,k1,1)
             sendBuffer(jj+1) = flowDoms(d1,level,mm)%x(i1,j1,k1,2)
             sendBuffer(jj+2) = flowDoms(d1,level,mm)%x(i1,j1,k1,3)
             jj = jj + 3

          enddo

          ! Send the data.

          call mpi_isend(sendBuffer(ii), size, adflow_real, procID,    &
               procID, ADflow_comm_world, sendRequests(i), &
               ierr)
          call EChk(ierr,__FILE__,__LINE__)

          ! Set ii to jj for the next processor.

          ii = jj

       enddo sends

       ! Post the nonblocking receives.

       ii = 1
       receives: do i=1,commPatternNode_1st(level)%nProcRecv

          ! Store the processor id and the size of the message
          ! a bit easier.

          procID = commPatternNode_1st(level)%recvProc(i)
          size   = 3*commPatternNode_1st(level)%nRecv(i)

          ! Post the receive.

          call mpi_irecv(recvBuffer(ii), size, adflow_real, procID, &
               myID, ADflow_comm_world, recvRequests(i), ierr)
          call EChk(ierr,__FILE__,__LINE__)

          ! And update ii.

          ii = ii + size

       enddo receives

       ! Copy the local data.
       !DIR$ NOVECTOR
       localCopy: do i=1,internalNode_1st(level)%nCopy

          ! Store the block and the indices of the donor a bit easier.

          d1 = internalNode_1st(level)%donorBlock(i)
          i1 = internalNode_1st(level)%donorIndices(i,1)
          j1 = internalNode_1st(level)%donorIndices(i,2)
          k1 = internalNode_1st(level)%donorIndices(i,3)
          ! Idem for the halo's.

          d2 = internalNode_1st(level)%haloBlock(i)
          i2 = internalNode_1st(level)%haloIndices(i,1)
          j2 = internalNode_1st(level)%haloIndices(i,2)
          k2 = internalNode_1st(level)%haloIndices(i,3)
          ! Copy the coordinates.
          flowDoms(d2,level,mm)%x(i2,j2,k2,1) = &
               flowDoms(d1,level,mm)%x(i1,j1,k1,1)
          flowDoms(d2,level,mm)%x(i2,j2,k2,2) = &
               flowDoms(d1,level,mm)%x(i1,j1,k1,2)
          flowDoms(d2,level,mm)%x(i2,j2,k2,3) = &
               flowDoms(d1,level,mm)%x(i1,j1,k1,3)

       enddo localCopy

       ! Correct the periodic halos of the internal communication
       ! pattern

       call correctPeriodicCoor(level, mm,                          &
            internalNode_1st(level)%nPeriodic,  &
            internalNode_1st(level)%periodicData)

       ! Complete the nonblocking receives in an arbitrary sequence and
       ! copy the coordinates from the buffer into the halo's.

       size = commPatternNode_1st(level)%nProcRecv
       completeRecvs: do i=1,commPatternNode_1st(level)%nProcRecv

          ! Complete any of the requests.

          call mpi_waitany(size, recvRequests, index, mpiStatus, ierr)
          call EChk(ierr,__FILE__,__LINE__)

          ! Copy the data just arrived in the halo's.

          ii = index
          jj = 3*commPatternNode_1st(level)%nRecvCum(ii-1) +1
          !DIR$ NOVECTOR
          do j=1,commPatternNode_1st(level)%nRecv(ii)

             ! Store the block and the indices of the halo a bit easier.

             d2 = commPatternNode_1st(level)%recvList(ii)%block(j)
             i2 = commPatternNode_1st(level)%recvList(ii)%indices(j,1)
             j2 = commPatternNode_1st(level)%recvList(ii)%indices(j,2)
             k2 = commPatternNode_1st(level)%recvList(ii)%indices(j,3)

             ! Copy the data.

             flowDoms(d2,level,mm)%x(i2,j2,k2,1) = recvBuffer(jj)
             flowDoms(d2,level,mm)%x(i2,j2,k2,2) = recvBuffer(jj+1)
             flowDoms(d2,level,mm)%x(i2,j2,k2,3) = recvBuffer(jj+2)
             jj = jj + 3

          enddo

       enddo completeRecvs

       ! Correct the periodic halos of the external communication
       ! pattern.

       call correctPeriodicCoor(level, mm,                            &
            commPatternNode_1st(level)%nPeriodic, &
            commPatternNode_1st(level)%periodicData)

       ! Complete the nonblocking sends.

       size = commPatternNode_1st(level)%nProcSend
       do i=1,commPatternNode_1st(level)%nProcSend
          call mpi_waitany(size, sendRequests, index, mpiStatus, ierr)
          call EChk(ierr,__FILE__,__LINE__)
    enddo

    enddo spectralLoop

  end subroutine exchangeCoor

  !      ==================================================================

  subroutine correctPeriodicCoor(level, sp, nPeriodic, periodicData)
    !
    !       correctPeriodicCoor applies the periodic transformation to
    !       the coordinates of the nodal halo's in periodicData.
    !
    use block
    use communication
    implicit none
    !
    !      Subroutine arguments
    !
    integer(kind=intType), intent(in) :: level, sp, nPeriodic
    type(periodicDataType), dimension(:), pointer :: periodicData
    !
    !      Local variables.
    !
    integer(kind=intType) :: nn, mm, ii, i, j, k
    real(kind=realType)   :: dx, dy, dz

    real(kind=realType), dimension(3,3) :: rotMatrix
    real(kind=realType), dimension(3)   :: rotCenter, translation

    ! Loop over the number of periodic transformations.

    do nn=1,nPeriodic

       ! Store the rotation matrix, rotation center and translation
       ! vector a bit easier.

       rotMatrix   = periodicData(nn)%rotMatrix
       rotCenter   = periodicData(nn)%rotCenter
       translation = periodicData(nn)%translation + rotCenter

       ! Loop over the number of halo nodes for this transformation.
       !DIR$ NOVECTOR
       do ii=1,periodicData(nn)%nHalos

          ! Store the block and the indices a bit easier.

          mm = periodicData(nn)%block(ii)
          i  = periodicData(nn)%indices(ii,1)
          j  = periodicData(nn)%indices(ii,2)
          k  = periodicData(nn)%indices(ii,3)

          ! Determine the vector from the center of rotation to the
          ! uncorrected halo value.

          dx = flowDoms(mm,level,sp)%x(i,j,k,1) - rotCenter(1)
          dy = flowDoms(mm,level,sp)%x(i,j,k,2) - rotCenter(2)
          dz = flowDoms(mm,level,sp)%x(i,j,k,3) - rotCenter(3)

          ! Compute the corrected coordinates.

          flowDoms(mm,level,sp)%x(i,j,k,1) = rotMatrix(1,1)*dx &
               + rotMatrix(1,2)*dy &
               + rotMatrix(1,3)*dz &
               + translation(1)
          flowDoms(mm,level,sp)%x(i,j,k,2) = rotMatrix(2,1)*dx &
               + rotMatrix(2,2)*dy &
               + rotMatrix(2,3)*dz &
               + translation(2)
          flowDoms(mm,level,sp)%x(i,j,k,3) = rotMatrix(3,1)*dx &
               + rotMatrix(3,2)*dy &
               + rotMatrix(3,3)*dz &
               + translation(3)
       enddo
    enddo

  end subroutine correctPeriodicCoor
  subroutine exchangeCoor_b(level)
    !
    !       ExchangeCoor_b exchanges the *derivatives* of the given grid
    !       level IN REVERSE MODE.
    !
    use constants
    use block
    use communication
    use inputTimeSpectral
    use utils, only : EChk
    implicit none
    !
    !      Subroutine arguments.
    !
    integer(kind=intType), intent(in) :: level
    !
    !      Local variables.
    !
    integer :: size, procID, ierr, index
    integer, dimension(mpi_status_size) :: mpiStatus

    integer(kind=intType) :: i, j, ii, jj, mm, idim
    integer(kind=intType) :: d1, i1, j1, k1, d2, i2, j2, k2


    ! Loop over the number of spectral solutions.

    spectralLoop: do mm=1,nTimeIntervalsSpectral

       ! Send the coordinates i have to send. The data is first copied
       ! into the send buffer and this buffer is sent.

       ii = 1
       jj = 1
       recvs: do i=1,commPatternNode_1st(level)%nProcRecv

          ! Store the processor id and the size of the message
          ! a bit easier.

          procID = commPatternNode_1st(level)%recvProc(i)
          size   = 3*commPatternNode_1st(level)%nRecv(i)

          ! Copy the data in the correct part of the send buffer.
          !DIR$ NOVECTOR
          do j=1,commPatternNode_1st(level)%nRecv(i)

             ! Store the block id and the indices of the donor
             ! a bit easier.

             d1 = commPatternNode_1st(level)%recvList(i)%block(j)
             i1 = commPatternNode_1st(level)%recvList(i)%indices(j,1)
             j1 = commPatternNode_1st(level)%recvList(i)%indices(j,2)
             k1 = commPatternNode_1st(level)%recvList(i)%indices(j,3)

             ! Copy the coordinates of this point in the buffer.
             ! Update the counter jj accordingly.

             recvBuffer(jj)   = flowDomsd(d1,level,mm)%x(i1,j1,k1,1)
             recvBuffer(jj+1) = flowDomsd(d1,level,mm)%x(i1,j1,k1,2)
             recvBuffer(jj+2) = flowDomsd(d1,level,mm)%x(i1,j1,k1,3)
             jj = jj + 3
             flowDomsd(d1, level, mm)%x(i1,j1,k1,:) = zero
          enddo

          ! Send the data.

          call mpi_isend(recvBuffer(ii), size, adflow_real, procID,    &
               procID, ADflow_comm_world, recvRequests(i), &
               ierr)
          call EChk(ierr,__FILE__,__LINE__)

          ! Set ii to jj for the next processor.

          ii = jj

       enddo recvs

       ! Post the nonblocking receives.

       ii = 1
       send: do i=1,commPatternNode_1st(level)%nProcSend

          ! Store the processor id and the size of the message
          ! a bit easier.

          procID = commPatternNode_1st(level)%sendProc(i)
          size   = 3*commPatternNode_1st(level)%nSend(i)

          ! Post the receive.

          call mpi_irecv(sendBuffer(ii), size, adflow_real, procID, &
               myID, ADflow_comm_world, sendRequests(i), ierr)
          call EChk(ierr,__FILE__,__LINE__)

          ! And update ii.

          ii = ii + size

       enddo send

       ! Copy the local data.
       !DIR$ NOVECTOR
       localCopy: do i=1,internalNode_1st(level)%nCopy

          ! Store the block and the indices of the donor a bit easier.

          d1 = internalNode_1st(level)%donorBlock(i)
          i1 = internalNode_1st(level)%donorIndices(i,1)
          j1 = internalNode_1st(level)%donorIndices(i,2)
          k1 = internalNode_1st(level)%donorIndices(i,3)
          ! Idem for the halo's.

          d2 = internalNode_1st(level)%haloBlock(i)
          i2 = internalNode_1st(level)%haloIndices(i,1)
          j2 = internalNode_1st(level)%haloIndices(i,2)
          k2 = internalNode_1st(level)%haloIndices(i,3)

          ! Sum into the '1' values fro the '2' values
          do idim=1,3
             flowDomsd(d1,level,mm)%x(i1,j1,k1,idim) = flowDomsd(d1,level,mm)%x(i1,j1,k1,idim) + &
                  flowDomsd(d2,level,mm)%x(i2,j2,k2,idim)
             flowDomsd(d2, level, mm)%x(i2,j2,k2,idim) = zero
          end do
       enddo localCopy

       ! Correct the periodic halos of the internal communication
       ! pattern

       ! NOT IMPLEMENTED
       ! call correctPeriodicCoor(level, mm,                          &
       !      internalNode_1st(level)%nPeriodic,  &
       !      internalNode_1st(level)%periodicData)

       ! Complete the nonblocking receives in an arbitrary sequence and
       ! copy the coordinates from the buffer into the halo's.

       size = commPatternNode_1st(level)%nProcSend
       completeSends: do i=1,commPatternNode_1st(level)%nProcSend

          ! Complete any of the requests.

          call mpi_waitany(size, sendRequests, index, mpiStatus, ierr)
          call EChk(ierr,__FILE__,__LINE__)

          ! Copy the data just arrived in the halo's.

          ii = index
          jj = 3*commPatternNode_1st(level)%nSendCum(ii-1)
          !DIR$ NOVECTOR
          do j=1,commPatternNode_1st(level)%nSend(ii)

             ! Store the block and the indices of the halo a bit easier.

             d2 = commPatternNode_1st(level)%sendList(ii)%block(j)
             i2 = commPatternNode_1st(level)%sendList(ii)%indices(j,1)
             j2 = commPatternNode_1st(level)%sendList(ii)%indices(j,2)
             k2 = commPatternNode_1st(level)%sendList(ii)%indices(j,3)

             ! Sum into the '2' values from the recv buffer
             do idim=1,3
                flowDomsd(d2,level,mm)%x(i2,j2,k2,idim) = flowDomsd(d2,level,mm)%x(i2,j2,k2,idim) + &
                     sendBuffer(jj + idim )
             end do
             jj = jj + 3

          enddo

       enddo completeSends

       ! Correct the periodic halos of the external communication
       ! pattern.
       ! NOT IMLEMENTED
       ! call correctPeriodicCoor(level, mm,                            &
       !      commPatternNode_1st(level)%nPeriodic, &
       !      commPatternNode_1st(level)%periodicData)

       ! Complete the nonblocking sends.

       size = commPatternNode_1st(level)%nProcRecv
       do i=1,commPatternNode_1st(level)%nProcRecv
          call mpi_waitany(size, recvRequests, index, mpiStatus, ierr)
          call EChk(ierr,__FILE__,__LINE__)
       enddo

    enddo spectralLoop

  end subroutine exchangeCoor_b
  subroutine exchangeCoor_d(level)
    !
    !       ExchangeCoor_d exchanges the *derivatives* of the given grid
    !       level.
    !
    use block
    use communication
    use inputTimeSpectral
    use utils, only : EChk
    implicit none
    !
    !      Subroutine arguments.
    !
    integer(kind=intType), intent(in) :: level
    !
    !      Local variables.
    !
    integer :: size, procID, ierr, index
    integer, dimension(mpi_status_size) :: mpiStatus

    integer(kind=intType) :: i, j, ii, jj, mm
    integer(kind=intType) :: d1, i1, j1, k1, d2, i2, j2, k2


    ! Loop over the number of spectral solutions.

    spectralLoop: do mm=1,nTimeIntervalsSpectral

       ! Send the coordinates i have to send. The data is first copied
       ! into the send buffer and this buffer is sent.

       ii = 1
       sends: do i=1,commPatternNode_1st(level)%nProcSend

          ! Store the processor id and the size of the message
          ! a bit easier.

          procID = commPatternNode_1st(level)%sendProc(i)
          size   = 3*commPatternNode_1st(level)%nSend(i)

          ! Copy the data in the correct part of the send buffer.

          jj = ii
          !DIR$ NOVECTOR
          do j=1,commPatternNode_1st(level)%nSend(i)

             ! Store the block id and the indices of the donor
             ! a bit easier.

             d1 = commPatternNode_1st(level)%sendList(i)%block(j)
             i1 = commPatternNode_1st(level)%sendList(i)%indices(j,1)
             j1 = commPatternNode_1st(level)%sendList(i)%indices(j,2)
             k1 = commPatternNode_1st(level)%sendList(i)%indices(j,3)

             ! Copy the coordinates of this point in the buffer.
             ! Update the counter jj accordingly.

             sendBuffer(jj)   = flowDomsd(d1,level,mm)%x(i1,j1,k1,1)
             sendBuffer(jj+1) = flowDomsd(d1,level,mm)%x(i1,j1,k1,2)
             sendBuffer(jj+2) = flowDomsd(d1,level,mm)%x(i1,j1,k1,3)
             jj = jj + 3

          enddo

          ! Send the data.

          call mpi_isend(sendBuffer(ii), size, adflow_real, procID,    &
               procID, ADflow_comm_world, sendRequests(i), &
               ierr)
          call EChk(ierr,__FILE__,__LINE__)

          ! Set ii to jj for the next processor.

          ii = jj

       enddo sends

       ! Post the nonblocking receives.

       ii = 1
       receives: do i=1,commPatternNode_1st(level)%nProcRecv

          ! Store the processor id and the size of the message
          ! a bit easier.

          procID = commPatternNode_1st(level)%recvProc(i)
          size   = 3*commPatternNode_1st(level)%nRecv(i)

          ! Post the receive.

          call mpi_irecv(recvBuffer(ii), size, adflow_real, procID, &
               myID, ADflow_comm_world, recvRequests(i), ierr)
          call EChk(ierr,__FILE__,__LINE__)

          ! And update ii.

          ii = ii + size

       enddo receives

       ! Copy the local data.
       !DIR$ NOVECTOR
       localCopy: do i=1,internalNode_1st(level)%nCopy

          ! Store the block and the indices of the donor a bit easier.

          d1 = internalNode_1st(level)%donorBlock(i)
          i1 = internalNode_1st(level)%donorIndices(i,1)
          j1 = internalNode_1st(level)%donorIndices(i,2)
          k1 = internalNode_1st(level)%donorIndices(i,3)
          ! Idem for the halo's.

          d2 = internalNode_1st(level)%haloBlock(i)
          i2 = internalNode_1st(level)%haloIndices(i,1)
          j2 = internalNode_1st(level)%haloIndices(i,2)
          k2 = internalNode_1st(level)%haloIndices(i,3)
          ! Copy the coordinates.
          flowDomsd(d2,level,mm)%x(i2,j2,k2,1) = &
               flowDomsd(d1,level,mm)%x(i1,j1,k1,1)
          flowDomsd(d2,level,mm)%x(i2,j2,k2,2) = &
               flowDomsd(d1,level,mm)%x(i1,j1,k1,2)
          flowDomsd(d2,level,mm)%x(i2,j2,k2,3) = &
               flowDomsd(d1,level,mm)%x(i1,j1,k1,3)

       enddo localCopy

       ! Correct the periodic halos of the internal communication
       ! pattern

       ! NOT IMPLEMENTED
       ! call correctPeriodicCoor(level, mm,                          &
       !      internalNode_1st(level)%nPeriodic,  &
       !      internalNode_1st(level)%periodicData)

       ! Complete the nonblocking receives in an arbitrary sequence and
       ! copy the coordinates from the buffer into the halo's.

       size = commPatternNode_1st(level)%nProcRecv
       completeRecvs: do i=1,commPatternNode_1st(level)%nProcRecv

          ! Complete any of the requests.

          call mpi_waitany(size, recvRequests, index, mpiStatus, ierr)
          call EChk(ierr,__FILE__,__LINE__)

          ! Copy the data just arrived in the halo's.

          ii = index
          jj = 3*commPatternNode_1st(level)%nRecvCum(ii-1) +1
          !DIR$ NOVECTOR
          do j=1,commPatternNode_1st(level)%nRecv(ii)

             ! Store the block and the indices of the halo a bit easier.

             d2 = commPatternNode_1st(level)%recvList(ii)%block(j)
             i2 = commPatternNode_1st(level)%recvList(ii)%indices(j,1)
             j2 = commPatternNode_1st(level)%recvList(ii)%indices(j,2)
             k2 = commPatternNode_1st(level)%recvList(ii)%indices(j,3)

             ! Copy the data.

             flowDomsd(d2,level,mm)%x(i2,j2,k2,1) = recvBuffer(jj)
             flowDomsd(d2,level,mm)%x(i2,j2,k2,2) = recvBuffer(jj+1)
             flowDomsd(d2,level,mm)%x(i2,j2,k2,3) = recvBuffer(jj+2)
             jj = jj + 3

          enddo

       enddo completeRecvs

       ! Correct the periodic halos of the external communication
       ! pattern.
       ! NOT IMLEMENTED
       ! call correctPeriodicCoor(level, mm,                            &
       !      commPatternNode_1st(level)%nPeriodic, &
       !      commPatternNode_1st(level)%periodicData)

       ! Complete the nonblocking sends.

       size = commPatternNode_1st(level)%nProcSend
       do i=1,commPatternNode_1st(level)%nProcSend
          call mpi_waitany(size, sendRequests, index, mpiStatus, ierr)
          call EChk(ierr,__FILE__,__LINE__)
    enddo

    enddo spectralLoop

  end subroutine exchangeCoor_d

  ! -----------------------------------------------------------------
  !              Comm routines for zippers
  ! -----------------------------------------------------------------

  subroutine flowIntegrationZipperComm(isInflow, vars, sps)

    ! This routine could technically be inside of the
    ! flowIntegrationZipper subroutine but we split it out becuase it
    ! has petsc comm stuff that will differentiate manually that
    ! tapenade doesn't need to see.

    use constants
    use blockPointers, only : BCFaceID, BCData, addGridVelocities, nDom, nBocos, BCType
    use BCPointers, only : sFace, ww1, ww2, pp1, pp2, gamma1, gamma2, xx
    use oversetData, only : zipperMeshes, zipperMesh
    use surfaceFamilies, only : familyExchange, BCFamExchange
    use utils, only : setPointers, setBCPointers, EChk

#include <petsc/finclude/petsc.h>
    use petsc
    implicit none

    ! Input variables
    logical, intent(in) :: isInflow
    real(kind=realType), dimension(:, :) :: vars
    integer(kind=intType) :: sps

    ! Working variables
    integer(kind=intType) :: ii, iVar, i,j, iBeg, iEnd, jBeg, jEnd, ierr, nn, mm
    real(kind=realType), dimension(:), pointer :: localPtr
    type(zipperMesh), pointer :: zipper
    type(familyExchange), pointer :: exch

    ! Set the zipper pointer to the zipper for inflow/outflow conditions
    if (isInflow) then
        zipper => zipperMeshes(iBCGroupInflow)
        exch => BCFamExchange(iBCGroupInflow, sps)
    else
        zipper => zipperMeshes(iBCGroupOutFlow)
        exch => BCFamExchange(iBCGroupOutflow, sps)
    end if

    ! Note that we can generate all the nodal values we need locally
    ! using simple arithematic averaging since they are all flow
    ! properties. There is no need to use the cellToNodeScatter stuff
    ! here like for the forces.
    varLoop: do iVar=1, nZIppFlowComm
       call vecGetArrayF90(exch%nodeValLocal, localPtr, ierr)
       call EChk(ierr,__FILE__,__LINE__)

       ii = 0
       domainLoop: do nn=1, nDom
          call setPointers(nn, 1, sps)
          bocoLoop: do mm=1, nBocos

             if (((BCType(mm) == SubsonicInflow .or. &
                  BCType(mm) == SupersonicInflow) .and. (isInflow)) &
                  .or. &
                  (BCType(mm) == SubsonicOutflow .or. &
                  BCType(mm) == SupersonicOutflow) .and. (.not. isInflow)) then

                call setBCPointers(mm, .True.)
                iBeg = BCdata(mm)%inBeg; iEnd=BCData(mm)%inEnd
                jBeg = BCdata(mm)%jnBeg; jEnd=BCData(mm)%jnEnd
                do j=jBeg, jEnd
                   do i=iBeg, iEnd
                      ii = ii + 1
                      select case(iVar)
                      case (iRho, iVx, iVy, iVz)
                         localPtr(ii) = eighth*(&
                              ww1(i, j,   iVar) + ww1(i+1, j,   iVar) + &
                              ww1(i, j+1, ivar) + ww1(i+1, j+1, iVar) + &
                              ww2(i, j,   iVar) + ww2(i+1, j,   iVar) + &
                              ww2(i, j+1, ivar) + ww2(i+1, j+1, iVar))
                      case (iZippFlowP)
                         localPtr(ii) = eighth*(&
                              pp1(i, j  ) + pp1(i+1, j  ) + &
                              pp1(i, j+1) + pp1(i+1, j+1) + &
                              pp2(i, j  ) + pp2(i+1, j  ) + &
                              pp2(i, j+1) + pp2(i+1, j+1))
                      case (iZippFlowGamma)
                         localPtr(ii) = eighth*(&
                              gamma1(i, j  ) + gamma1(i+1, j  ) + &
                              gamma1(i, j+1) + gamma1(i+1, j+1) + &
                              gamma2(i, j  ) + gamma2(i+1, j  ) + &
                              gamma2(i, j+1) + gamma2(i+1, j+1))
                      case (iZippFlowSface)
                         if (addGridVelocities) then
                            localPtr(ii) = fourth*(&
                                 sface(i, j  ) + sface(i+1, j  ) + &
                                 sface(i, j+1) + sface(i+1, j+1))
                         else
                            localPtr(ii) = zero
                         end if
                      case (iZippFlowX, iZippFlowY, iZippFlowZ)
                         localPtr(ii) = xx(i+1, j+1, iVar-iZippFlowX+1)
                      end select
                   end do
                end do
             end if
          end do bocoLoop
       end do domainLoop

       ! Return pointer to nodeValLocal
       call vecRestoreArrayF90(exch%nodeValLocal, localPtr, ierr)
       call EChk(ierr,__FILE__,__LINE__)

       call VecPlaceArray(zipper%localVal, vars(:, iVar), ierr)
       call EChk(ierr,__FILE__,__LINE__)

       ! Send these values to the root using the zipper scatter.
       call VecScatterBegin(zipper%scatter, exch%nodeValLocal,&
            zipper%localVal, INSERT_VALUES, SCATTER_FORWARD, ierr)
       call EChk(ierr,__FILE__,__LINE__)

       call VecScatterEnd(zipper%scatter, exch%nodeValLocal,&
            zipper%localVal, INSERT_VALUES, SCATTER_FORWARD, ierr)
       call EChk(ierr,__FILE__,__LINE__)

       ! Reset the original petsc vector.
       call vecResetArray(zipper%localVal, ierr)
       call EChk(ierr,__FILE__,__LINE__)
    end do varLoop

  end subroutine flowIntegrationZipperComm

  subroutine flowIntegrationZipperComm_d(isInflow, vars, varsd, sps)

    ! Forward mode linearization of flowIntegratoinZipperComm

    use constants
    use blockPointers, only : BCFaceID, BCData, addGridVelocities, nDom, nBocos, BCType
    use BCPointers, only : sFaced, ww1d, ww2d, pp1d, pp2d, xxd
    use oversetData, only : zipperMeshes, zipperMesh
    use surfaceFamilies, only : familyExchange, BCFamExchange
    use utils, only : setPointers_d, setBCPointers_d, EChk
#include <petsc/finclude/petsc.h>
    use petsc
    implicit none
    ! Input variables
    logical, intent(in) :: isInflow
    real(kind=realType), dimension(:, :) :: vars, varsd
    integer(kind=intType) :: sps

    ! Working variables
    integer(kind=intType) :: ii, iVar, i,j, iBeg, iEnd, jBeg, jEnd, ierr, nn, mm
    real(kind=realType), dimension(:), pointer :: localPtr
    type(zipperMesh), pointer :: zipper
    type(familyExchange), pointer :: exch



    ! Set the zipper pointer to the zipper for inflow/outflow conditions
    if (isInflow) then
        ! Need to generate the vars themselves.
        call flowIntegrationZipperComm(.true., vars, sps)
        zipper => zipperMeshes(iBCGroupInflow)
        exch => BCFamExchange(iBCGroupInflow, sps)
    else
        ! Need to generate the vars themselves.
        call flowIntegrationZipperComm(.false., vars, sps)
        zipper => zipperMeshes(iBCGroupOutFlow)
        exch => BCFamExchange(iBCGroupOutflow, sps)
    end if

    ! Note that we can generate all the nodal values we need locally
    ! using simple arithematic averaging since they are all flow
    ! properties. There is no need to use the cellToNodeScatter stuff
    ! here like for the forces.
    varLoop: do iVar=1, nZIppFlowComm
       call vecGetArrayF90(exch%nodeValLocal, localPtr, ierr)
       call EChk(ierr,__FILE__,__LINE__)

       ii = 0
       domainLoop: do nn=1, nDom
          call setPointers_d(nn, 1, sps)
          bocoLoop: do mm=1, nBocos
             if (((BCType(mm) == SubsonicInflow .or. &
                  BCType(mm) == SupersonicInflow) .and. isInflow) &
                  .or. &
                  (BCType(mm) == SubsonicOutflow .or. &
                  BCType(mm) == SupersonicOutflow) .and. (.not. isInflow)) then

                call setBCPointers_d(mm, .True.)
                iBeg = BCdata(mm)%inBeg; iEnd=BCData(mm)%inEnd
                jBeg = BCdata(mm)%jnBeg; jEnd=BCData(mm)%jnEnd
                do j=jBeg, jEnd
                   do i=iBeg, iEnd
                      ii = ii + 1
                      select case(iVar)
                      case (iRho, iVx, iVy, iVz)
                         localPtr(ii) = eighth*(&
                              ww1d(i, j,   iVar) + ww1d(i+1, j,   iVar) + &
                              ww1d(i, j+1, ivar) + ww1d(i+1, j+1, iVar) + &
                              ww2d(i, j,   iVar) + ww2d(i+1, j,   iVar) + &
                              ww2d(i, j+1, ivar) + ww2d(i+1, j+1, iVar))
                      case (iZippFlowP)
                         localPtr(ii) = eighth*(&
                              pp1d(i, j  ) + pp1d(i+1, j  ) + &
                              pp1d(i, j+1) + pp1d(i+1, j+1) + &
                              pp2d(i, j  ) + pp2d(i+1, j  ) + &
                              pp2d(i, j+1) + pp2d(i+1, j+1))
                      case (iZippFlowGamma)
                         ! Gamma is not currently an active variable
                         localPtr(ii) = zero
                      case (iZippFlowSface)
                         if (addGridVelocities) then
                            localPtr(ii) = fourth*(&
                                 sfaced(i, j  ) + sfaced(i+1, j  ) + &
                                 sfaced(i, j+1) + sfaced(i+1, j+1))
                         else
                            localPtr(ii) = zero
                         end if
                      case (iZippFlowX, iZippFlowY, iZippFlowZ)
                         localPtr(ii) = xxd(i+1, j+1, iVar-iZippFlowX+1)
                      end select
                   end do
                end do
             end if
          end do bocoLoop
       end do domainLoop

       ! Return pointer to nodeValLocal
       call vecRestoreArrayF90(exch%nodeValLocal, localPtr, ierr)
       call EChk(ierr,__FILE__,__LINE__)

       call VecPlaceArray(zipper%localVal, varsd(:, iVar), ierr)
       call EChk(ierr,__FILE__,__LINE__)

       ! Send these values to the root using the zipper scatter.
       call VecScatterBegin(zipper%scatter, exch%nodeValLocal,&
            zipper%localVal, INSERT_VALUES, SCATTER_FORWARD, ierr)
       call EChk(ierr,__FILE__,__LINE__)

       call VecScatterEnd(zipper%scatter, exch%nodeValLocal,&
            zipper%localVal, INSERT_VALUES, SCATTER_FORWARD, ierr)
       call EChk(ierr,__FILE__,__LINE__)

       ! Reset the original petsc vector.
       call vecResetArray(zipper%localVal, ierr)
       call EChk(ierr,__FILE__,__LINE__)
    end do varLoop

  end subroutine flowIntegrationZipperComm_d

  subroutine flowIntegrationZipperComm_b(isInflow, vars, varsd, sps)

    ! Reverse mode linearization of the flowIntegrationZipperComm routine

    use constants
    use blockPointers, only : BCFaceID, BCData, addGridVelocities, nDom, nBocos, BCType
    use BCPointers, only : sFaced, ww1d, ww2d, pp1d, pp2d, xxd
    use oversetData, only : zipperMeshes, zipperMesh
    use surfaceFamilies, only : familyExchange, BCFamExchange
    use utils, only : setPointers_b, setBCPointers_d, EChk
#include <petsc/finclude/petsc.h>
    use petsc
    implicit none

    ! Input variables
    logical, intent(in) :: isInflow
    real(kind=realType), dimension(:, :) :: vars, varsd
    integer(kind=intType) :: sps

    ! Working variables
    integer(kind=intType) :: ii, iVar, i,j, iBeg, iEnd, jBeg, jEnd, ierr, nn, mm
    real(kind=realType), dimension(:), pointer :: localPtr
    real(kind=realType) ::tmp
    type(zipperMesh), pointer :: zipper
    type(familyExchange), pointer :: exch

    ! Set the zipper pointer to the zipper for inflow/outflow conditions
    if (isInflow) then
        zipper => zipperMeshes(iBCGroupInflow)
        exch => BCFamExchange(iBCGroupInflow, sps)
    else
        zipper => zipperMeshes(iBCGroupOutFlow)
        exch => BCFamExchange(iBCGroupOutflow, sps)
    end if

    ! Run the var exchange loop backwards:
    varLoop: do iVar=1, nZIppFlowComm

       ! Zero the vector we are scatting into:
       call VecSet(exch%nodeValLocal, zero, ierr)
       call EChk(ierr,__FILE__,__LINE__)

       call VecPlaceArray(zipper%localVal, varsd(:, iVar), ierr)
       call EChk(ierr,__FILE__,__LINE__)

       ! Send these values to the root using the zipper scatter.
       call VecScatterBegin(zipper%scatter, zipper%localVal, &
            exch%nodeValLocal, ADD_VALUES, SCATTER_REVERSE, ierr)
       call EChk(ierr,__FILE__,__LINE__)

       call VecScatterEnd(zipper%scatter, zipper%localVal, &
            exch%nodeValLocal, ADD_VALUES, SCATTER_REVERSE, ierr)
       call EChk(ierr,__FILE__,__LINE__)

       ! Reset the original petsc vector.
       call vecResetArray(zipper%localVal, ierr)
       call EChk(ierr,__FILE__,__LINE__)

       ! To be consistent, varsd must be zeroed since it is "used"
       varsd(:, iVar) = zero

       ! Now finish the scatting back to the acutual BCs pointers (and
       ! thus the state variables).

       call vecGetArrayF90(exch%nodeValLocal, localPtr, ierr)
       call EChk(ierr,__FILE__,__LINE__)

       ii = 0
       domainLoop: do nn=1, nDom
          call setPointers_b(nn, 1, sps)
          bocoLoop: do mm=1, nBocos
             if (((BCType(mm) == SubsonicInflow .or. &
                  BCType(mm) == SupersonicInflow) .and. isInflow) &
                  .or. &
                  (BCType(mm) == SubsonicOutflow .or. &
                  BCType(mm) == SupersonicOutflow) .and. (.not. isInflow)) then

                call setBCPointers_d(mm, .True.)
                iBeg = BCdata(mm)%inBeg; iEnd=BCData(mm)%inEnd
                jBeg = BCdata(mm)%jnBeg; jEnd=BCData(mm)%jnEnd
                do j=jBeg, jEnd
                   do i=iBeg, iEnd
                      ii = ii + 1
                      select case(iVar)
                      case (iRho, iVx, iVy, iVz)
                         tmp = eighth * localPtr(ii)
                         ww1d(i  , j  , iVar) = ww1d(i  , j  , iVar) + tmp
                         ww1d(i+1, j  , iVar) = ww1d(i+1, j  , iVar) + tmp
                         ww1d(i  , j+1, iVar) = ww1d(i  , j+1, iVar) + tmp
                         ww1d(i+1, j+1, iVar) = ww1d(i+1, j+1, iVar) + tmp

                         ww2d(i  , j  , iVar) = ww2d(i  , j  , iVar) + tmp
                         ww2d(i+1, j  , iVar) = ww2d(i+1, j  , iVar) + tmp
                         ww2d(i  , j+1, iVar) = ww2d(i  , j+1, iVar) + tmp
                         ww2d(i+1, j+1, iVar) = ww2d(i+1, j+1, iVar) + tmp
                      case (iZippFlowP)
                         tmp = eighth * localPtr(ii)
                         pp1d(i  , j  ) = pp1d(i  , j  ) + tmp
                         pp1d(i+1, j  ) = pp1d(i+1, j  ) + tmp
                         pp1d(i  , j+1) = pp1d(i  , j+1) + tmp
                         pp1d(i+1, j+1) = pp1d(i+1, j+1) + tmp

                         pp2d(i  , j  ) = pp2d(i  , j  ) + tmp
                         pp2d(i+1, j  ) = pp2d(i+1, j  ) + tmp
                         pp2d(i  , j+1) = pp2d(i  , j+1) + tmp
                         pp2d(i+1, j+1) = pp2d(i+1, j+1) + tmp

                      case (iZippFlowGamma)
                         ! gamma is not currently active


                      case (iZippFlowSFace)
                         if (addGridVelocities) then
                            tmp = fourth*localPtr(ii)
                            sfaced(i  , j  ) = sfaced(i  , j  ) + tmp
                            sfaced(i+1, j  ) = sfaced(i+1, j  ) + tmp
                            sfaced(i  , j+1) = sfaced(i  , j+1) + tmp
                            sfaced(i+1, j+1) = sfaced(i+1, j+1) + tmp
                         else
                            localPtr(ii) = zero
                         end if
                      case (iZippFlowX, iZippFlowY, iZippFlowZ)
                         xxd(i+1, j+1, iVar-iZippFlowX+1) = xxd(i+1, j+1, iVar-iZippFlowX+1) + localPtr(ii)
                      end select
                   end do
                end do
             end if
          end do bocoLoop
       end do domainLoop

       ! Return pointer to nodeValLocal
       call vecRestoreArrayF90(exch%nodeValLocal, localPtr, ierr)
       call EChk(ierr,__FILE__,__LINE__)

    end do varLoop

  end subroutine flowIntegrationZipperComm_b

  subroutine wallIntegrationZipperComm(vars, sps)

    ! This routine could technically be inside of the
    ! flowIntegrationZipper subroutine but we split it out becuase it
    ! has petsc comm stuff that will differentiate manually that
    ! tapenade doesn't need to see.

    use constants
    use blockPointers, only : BCData, nDom, BCType, nBocos
    use BCPointers, only : xx
    use oversetData, only : zipperMeshes, zipperMesh
    use surfaceFamilies, only : familyExchange, BCFamExchange
    use utils, only : setPointers, setBCPointers, EChk, isWallType
#include <petsc/finclude/petsc.h>
    use petsc
    implicit none
    ! Input variables
    real(kind=realType), dimension(:, :) :: vars
    integer(kind=intType) :: sps

    ! Working variables
    integer(kind=intType) :: ii, iVar, i,j, iBeg, iEnd, jBeg, jEnd, ierr, nn, mm
    real(kind=realType), dimension(:), pointer :: localPtr
    type(zipperMesh), pointer :: zipper
    type(familyExchange), pointer :: exch

    ! Set the zipper pointer to the zipper for inflow/outflow conditions
    zipper => zipperMeshes(iBCGroupWalls)
    exch => BCFamExchange(iBCGroupWalls, sps)

    ! Make sure the nodal tractions are computed
    call computeNodalTractions(sps)

    varLoop: do iVar=1, nZippWallComm
       call vecGetArrayF90(exch%nodeValLocal, localPtr, ierr)
       call EChk(ierr,__FILE__,__LINE__)

       ii = 0
       domainLoop: do nn=1, nDom
          call setPointers(nn, 1, sps)
          bocoLoop: do mm=1, nBocos
             if (isWallType(BCType(mm))) then
                call setBCPointers(mm, .True.)
                iBeg = BCdata(mm)%inBeg; iEnd=BCData(mm)%inEnd
                jBeg = BCdata(mm)%jnBeg; jEnd=BCData(mm)%jnEnd
                do j=jBeg, jEnd
                   do i=iBeg, iEnd
                      ii = ii + 1
                      select case(iVar)

                      case (iZippWallTpx, iZippWallTpy, iZippWallTpz)

                         localPtr(ii) = BCData(mm)%Tp(i, j, iVar)

                      case (iZippWallTvx, iZippWallTvy, iZippWallTvz)

                         localPtr(ii) = BCData(mm)%Tv(i, j, iVar-iZippWallTvx+1)

                      case (iZippWallX, iZippWallY, iZippWallZ)

                         ! The +1 is due to pointer offset
                         localPtr(ii) = xx(i+1, j+1, iVar-iZippWallX+1)

                      end select
                   end do
                end do
             end if
          end do bocoLoop
       end do domainLoop

       ! Return pointer to nodeValLocal
       call vecRestoreArrayF90(exch%nodeValLocal, localPtr, ierr)
       call EChk(ierr,__FILE__,__LINE__)

       call VecPlaceArray(zipper%localVal, vars(:, iVar), ierr)
       call EChk(ierr,__FILE__,__LINE__)

       ! Send these values to the root using the zipper scatter.
       call VecScatterBegin(zipper%scatter, exch%nodeValLocal,&
            zipper%localVal, INSERT_VALUES, SCATTER_FORWARD, ierr)
       call EChk(ierr,__FILE__,__LINE__)

       call VecScatterEnd(zipper%scatter, exch%nodeValLocal,&
            zipper%localVal, INSERT_VALUES, SCATTER_FORWARD, ierr)
       call EChk(ierr,__FILE__,__LINE__)

       ! Reset the original petsc vector.
       call vecResetArray(zipper%localVal, ierr)
       call EChk(ierr,__FILE__,__LINE__)
    end do varLoop

  end subroutine wallIntegrationZipperComm

  subroutine wallIntegrationZipperComm_d(vars, varsd, sps)

    ! Forward mode linearization of the wallIntegrationZipperComm

    use constants
    use blockPointers, only : BCDatad, BCData, nBocos, nDom, BCType
    use BCPointers, only : xxd
    use oversetData, only : zipperMeshes, zipperMesh
    use surfaceFamilies, only : familyExchange, BCFamExchange
    use utils, only : setPointers_d, setBCPointers_d, EChk, isWallType
#include <petsc/finclude/petsc.h>
    use petsc
    implicit none

    ! Input variables
    real(kind=realType), dimension(:, :):: vars, varsd
    integer(kind=intType) :: sps

    ! Working variables
    integer(kind=intType) :: ii, iVar, i,j, iBeg, iEnd, jBeg, jEnd, ierr, nn, mm
    real(kind=realType), dimension(:), pointer :: localPtr
    type(zipperMesh), pointer :: zipper
    type(familyExchange), pointer :: exch

    ! Need to set the actual variables first.
    call wallIntegrationZIpperCOmm(vars, sps)

    ! Compute the derivative of the nodal tractions
    call computeNodalTractions_d(sps)

    ! Set the zipper pointer to the zipper for inflow/outflow conditions
    zipper => zipperMeshes(iBCGroupWalls)
    exch => BCFamExchange(iBCGroupWalls, sps)

    varLoop: do iVar=1, nZippWallComm
       call vecGetArrayF90(exch%nodeValLocal, localPtr, ierr)
       call EChk(ierr,__FILE__,__LINE__)

       ii = 0
       domainLoop: do nn=1, nDom
          call setPointers_d(nn, 1, sps)
          bocoLoop: do mm=1, nBocos
             if (isWallType(BCType(mm))) then
                call setBCPointers_d(mm, .True.)
                iBeg = BCdata(mm)%inBeg; iEnd=BCData(mm)%inEnd
                jBeg = BCdata(mm)%jnBeg; jEnd=BCData(mm)%jnEnd
                do j=jBeg, jEnd
                   do i=iBeg, iEnd
                      ii = ii + 1
                      select case(iVar)

                      case (iZippWallTpx, iZippWallTpy, iZippWallTpz)

                         localPtr(ii) = BCDatad(mm)%Tp(i, j, iVar)

                      case (iZippWallTvx, iZippWallTvy, iZippWallTvz)

                         localPtr(ii) = BCDatad(mm)%Tv(i, j, iVar-iZippWallTvx+1)

                      case (iZippWallX, iZippWallY, iZippWallZ)

                         ! The +1 is due to pointer offset
                         localPtr(ii) = xxd(i+1, j+1, iVar-iZippWallX+1)

                      end select
                   end do
                end do
             end if
          end do bocoLoop
       end do domainLoop

       ! Return pointer to nodeValLocal
       call vecRestoreArrayF90(exch%nodeValLocal, localPtr, ierr)
       call EChk(ierr,__FILE__,__LINE__)

       call VecPlaceArray(zipper%localVal, varsd(:, iVar), ierr)
       call EChk(ierr,__FILE__,__LINE__)

       ! Send these values to the root using the zipper scatter.
       call VecScatterBegin(zipper%scatter, exch%nodeValLocal,&
            zipper%localVal, INSERT_VALUES, SCATTER_FORWARD, ierr)
       call EChk(ierr,__FILE__,__LINE__)

       call VecScatterEnd(zipper%scatter, exch%nodeValLocal,&
            zipper%localVal, INSERT_VALUES, SCATTER_FORWARD, ierr)
       call EChk(ierr,__FILE__,__LINE__)

       ! Reset the original petsc vector.
       call vecResetArray(zipper%localVal, ierr)
       call EChk(ierr,__FILE__,__LINE__)
    end do varLoop

  end subroutine wallIntegrationZipperComm_d

  subroutine wallIntegrationZipperComm_b(vars, varsd, sps)

    ! Reverse mode linearization of the wallIntegrationZipperComm

    use constants
    use blockPointers, only : BCDatad, BCData, nBocos, nDom, BCType
    use BCPointers, only : xxd
    use oversetData, only : zipperMeshes, zipperMesh
    use surfaceFamilies, only : familyExchange, BCFamExchange
    use utils, only : setPointers_b, setBCPointers_d, EChk, isWallType
#include <petsc/finclude/petsc.h>
    use petsc
    implicit none

    ! Input variables
    real(kind=realType), dimension(:, :) :: vars, varsd
    integer(kind=intType) :: sps

    ! Working variables
    integer(kind=intType) :: ii, iVar, i,j, iBeg, iEnd, jBeg, jEnd, ierr, nn, mm
    real(kind=realType), dimension(:), pointer :: localPtr
    type(zipperMesh), pointer :: zipper
    type(familyExchange), pointer :: exch

    ! Set the zipper pointer to the zipper for inflow/outflow conditions
    zipper => zipperMeshes(iBCGroupWalls)
    exch => BCFamExchange(iBCGroupWalls, sps)

    ! Run the var exchange loop backwards:
    varLoop: do iVar=1, nZippWallComm

       ! Zero the vector we are scatting into:
       call VecSet(exch%nodeValLocal, zero, ierr)
       call EChk(ierr,__FILE__,__LINE__)

       call VecPlaceArray(zipper%localVal, varsd(:, iVar), ierr)
       call EChk(ierr,__FILE__,__LINE__)

       ! Send these values to the root using the zipper scatter.
       call VecScatterBegin(zipper%scatter, zipper%localVal, &
            exch%nodeValLocal, ADD_VALUES, SCATTER_REVERSE, ierr)
       call EChk(ierr,__FILE__,__LINE__)

       call VecScatterEnd(zipper%scatter, zipper%localVal, &
            exch%nodeValLocal, ADD_VALUES, SCATTER_REVERSE, ierr)
       call EChk(ierr,__FILE__,__LINE__)

       ! Reset the original petsc vector.
       call vecResetArray(zipper%localVal, ierr)
       call EChk(ierr,__FILE__,__LINE__)

       ! To be consistent, varsd must be zeroed since it is "used"
       varsd(:, iVar) = zero

       ! Now finish the scatting back to the acutual BCs pointers (and
       ! thus the state variables).

       call vecGetArrayF90(exch%nodeValLocal, localPtr, ierr)
       call EChk(ierr,__FILE__,__LINE__)

       ii = 0
       domainLoop: do nn=1, nDom
          call setPointers_b(nn, 1, sps)
          bocoLoop: do mm=1, nBocos
             if (isWallType(BCType(mm))) then
                call setBCPointers_d(mm, .True.)
                iBeg = BCdata(mm)%inBeg; iEnd=BCData(mm)%inEnd
                jBeg = BCdata(mm)%jnBeg; jEnd=BCData(mm)%jnEnd
                do j=jBeg, jEnd
                   do i=iBeg, iEnd
                      ii = ii + 1
                      select case(iVar)
                      case (iZippWallTpx, iZippWallTpy, iZippWallTpz)

                         BCDatad(mm)%Tp(i, j, iVar) = localPtr(ii)

                      case (iZippWallTvx, iZippWallTvy, iZippWallTvz)

                         BCDatad(mm)%Tv(i, j, iVar-iZippWallTvx+1) = localPtr(ii)

                      case (iZippWallX, iZippWallY, iZippWallZ)

                         ! The +1 is due to pointer offset
                         xxd(i+1, j+1, iVar-iZippWallX+1) = xxd(i+1, j+1, iVar-iZippWallX+1) + localPtr(ii)

                      end select
                   end do
                end do
             end if
          end do bocoLoop
       end do domainLoop

       ! Return pointer to nodeValLocal
       call vecRestoreArrayF90(exch%nodeValLocal, localPtr, ierr)
       call EChk(ierr,__FILE__,__LINE__)

    end do varLoop

    ! Compute the derivatives of the nodal tractions. The will
    !accumulate the seeds onto bcDatad%Fv, bcDatad%Fv and bcDatad%area

    call computeNodalTractions_b(sps)

  end subroutine wallIntegrationZipperComm_b
end module haloExchange

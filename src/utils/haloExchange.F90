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

    call whalo1to1(level, start, end, commPressure, commVarGamma, &
         commLamVis, commEddyVis, commPatternCell_2nd,  &
         internalCell_2nd)

    ! Exchange the overset cells

    mm = ubound(commPatternOverset, 1)
    call wOverset(level, start, end, commPressure, commVarGamma, &
         commLamVis, commEddyVis, commPatternOverset, internalOverset, mm)

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

       ! Set the pointers for the required variables
       domainLoop:do nn=1, nDom
          nVar = 0

          do k=start, end
             nVar = nVar + 1 
             flowDoms(nn, level, sps)%realCommVars(nVar)%var => &
                  flowDoms(nn, level, sps)%w(:, :, :, k)
          end do

          if( commPressure )  then 
             nVar = nVar + 1
             flowDoms(nn, level, sps)%realCommVars(nVar)%var => &
                  flowDoms(nn, level, sps)%P(:, :, :)
          end if

          if( commVarGamma ) then 
             nVar = nVar + 1
             flowDoms(nn, level, sps)%realCommVars(nVar)%var => &
                  flowDoms(nn, 1, sps)%gamma(:, :, :)
          end if

          if( commLamVis ) then 
             nVar = nVar + 1
             flowDoms(nn, level, sps)%realCommVars(nVar)%var => &
                  flowDoms(nn, 1, sps)%rlv(:, :, :)
          end if

          if( commEddyVis ) then 
             nVar = nVar + 1
             flowDoms(nn, level, sps)%realCommVars(nVar)%var => &
                  flowDoms(nn, level, sps)%rev(:, :, :)
          end if

          if(nVar == 0) return
       end do domainLoop

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
    integer, dimension(mpi_status_size) :: status

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
       do j=1,commPattern(level)%nsend(i)

          ! Store the block id and the indices of the donor
          ! a bit easier.

          d1 = commPattern(level)%sendList(i)%block(j)
          i1 = commPattern(level)%sendList(i)%indices(j,1)+1
          j1 = commPattern(level)%sendList(i)%indices(j,2)+1
          k1 = commPattern(level)%sendList(i)%indices(j,3)+1

          ! Copy the given range of the working variables for
          ! this cell in the buffer. Update the counter jj.
          do k=1, nvar
             sendBuffer(jj) = flowDoms(d1, level, sps)%realCommVars(k)%var(i1, j1, k1)
             jj = jj + 1
          end do
       end do

       ! Send the data.

       call mpi_isend(sendBuffer(ii), size, adflow_real, procID,  &
            procID, ADflow_comm_world, sendRequests(i), &
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

       call mpi_irecv(recvBuffer(ii), size, adflow_real, procID, &
            myID, ADflow_comm_world, recvRequests(i), ierr)

       ! And update ii.

       ii = ii + size

    enddo receives

    ! Copy the local data.

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

       call mpi_waitany(size, recvRequests, index, status, ierr)

       ! Copy the data just arrived in the halo's.

       ii = index
       jj = nVar*commPattern(level)%nrecvCum(ii-1)
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
       call mpi_waitany(size, sendRequests, index, status, ierr)
    enddo

  end subroutine whalo1to1RealGeneric


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
    integer, dimension(mpi_status_size) :: status

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
       do j=1,commPattern(level)%nsend(i)

          ! Store the block id and the indices of the donor
          ! a bit easier.

          d1 = commPattern(level)%sendList(i)%block(j)
          i1 = commPattern(level)%sendList(i)%indices(j,1)+1
          j1 = commPattern(level)%sendList(i)%indices(j,2)+1
          k1 = commPattern(level)%sendList(i)%indices(j,3)+1

          ! Copy the given range of the working variables for
          ! this cell in the buffer. Update the counter jj.
          do k=1, nvar
             sendBufInt(jj) = flowDoms(d1, level, sps)%intCommVars(k)%var(i1, j1, k1)
             jj = jj + 1
          end do
       end do

       ! Send the data.

       call mpi_isend(sendBufInt(ii), size, adflow_integer, procID,  &
            procID, ADflow_comm_world, sendRequests(i), &
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

       call mpi_irecv(recvBufInt(ii), size, adflow_integer, procID, &
            myID, ADflow_comm_world, recvRequests(i), ierr)

       ! And update ii.

       ii = ii + size

    enddo receives

    ! Copy the local data.

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

       call mpi_waitany(size, recvRequests, index, status, ierr)

       ! Copy the data just arrived in the halo's.

       ii = index
       jj = nVar*commPattern(level)%nrecvCum(ii-1)
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
       call mpi_waitany(size, sendRequests, index, status, ierr)
    enddo

    deallocate(recvBufInt, sendBufInt)

  end subroutine whalo1to1IntGeneric

  subroutine whalo1to1_b(level, start, end, commPressure,       &
       commVarGamma, commLamVis, commEddyVis, &
       commPattern, internal)
    !
    !       whalo1to1b performs the *TRANSPOSE* operation of whalo1to1.    
    !       It is used for adjoint/reverse mode residual evaluations.      
    !       See whalo1to1 for more information. Note that this code does   
    !       include the correctPeroidicVelocity computation.               
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
    !
    !      Local variables.
    !
    integer :: size, procID, ierr, index
    integer, dimension(mpi_status_size) :: status

    integer(kind=intType) :: nVar, mm
    integer(kind=intType) :: i, j, k, ii, jj
    integer(kind=intType) :: d1, i1, j1, k1, d2, i2, j2, k2

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

          procID = commPattern(level)%recvProc(i)
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
          call mpi_isend(recvBuffer(ii), size, adflow_real, procID,  &
               procID, ADflow_comm_world, recvRequests(i), &
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

          call mpi_irecv(sendBuffer(ii), size, adflow_real, procID, &
               myID, ADflow_comm_world, sendRequests(i), ierr)

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
             flowDomsd(d2,level,mm)%w(i2,j2,k2,k) = zero
          enddo

          ! The pressure, if needed.

          if( commPressure ) then
             flowDomsd(d1,level,mm)%p(i1,j1,k1) = flowDomsd(d1,level,mm)%p(i1,j1,k1) + &
                  flowDomsd(d2,level,mm)%p(i2,j2,k2)
             flowDomsd(d2,level,mm)%p(i2,j2,k2) = zero
          end if

          ! The specific heat ratio, if needed. Note that level == 1.

          if( commVarGamma ) then
             flowDomsd(d1,1,mm)%gamma(i1,j1,k1) = flowDomsd(d1,1,mm)%gamma(i1,j1,k1) + &
                  flowDomsd(d2,1,mm)%gamma(i2,j2,k2)
             flowDomsd(d2,1,mm)%gamma(i2,j2,k2) = zero
          end if

          ! The laminar viscosity for viscous computations.
          ! Again level == 1.

          if( commLamVis ) then
             flowDomsd(d1,1,mm)%rlv(i1,j1,k1) = flowDomsd(d1,1,mm)%rlv(i1,j1,k1) + &
                  flowDomsd(d2,1,mm)%rlv(i2,j2,k2)
             flowDomsd(d2,1,mm)%rlv(i2,j2,k2) = zero
          end if

          ! The eddy viscosity for eddy viscosity models.
          ! Level is the true multigrid level, because the eddy
          ! viscosity is allocated on all grid levels.

          if( commEddyVis ) then
             flowDomsd(d1,level,mm)%rev(i1,j1,k1) = flowDomsd(d1,level,mm)%rev(i1,j1,k1) + &
                  flowDomsd(d2,level,mm)%rev(i2,j2,k2)
             flowDomsd(d2,level,mm)%rev(i2,j2,k2) = zero
          end if

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

  subroutine whalo1to1_d(level, start, end, commPressure,       &
       commVarGamma, commLamVis, commEddyVis, &
       commPattern, internal)
    !
    !       whalo1to1 exchanges the 1 to 1 internal halo's derivatives     
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

    spectralModes: do sps=1,nTimeIntervalsSpectral

       ! Set the pointers for the required variables
       domainLoop:do nn=1, nDom
          nVar = 0

          do k=start, end
             nVar = nVar + 1 
             flowDoms(nn, level, sps)%realCommVars(nVar)%var => &
                  flowDomsd(nn, level, sps)%w(:, :, :, k)
          end do

          if( commPressure )  then 
             nVar = nVar + 1
             flowDoms(nn, level, sps)%realCommVars(nVar)%var => &
                  flowDomsd(nn, level, sps)%P(:, :, :)
          end if

          if( commVarGamma ) then 
             nVar = nVar + 1
             flowDoms(nn, level, sps)%realCommVars(nVar)%var => &
                  flowDomsd(nn, 1, sps)%gamma(:, :, :)
          end if

          if( commLamVis ) then 
             nVar = nVar + 1
             flowDoms(nn, level, sps)%realCommVars(nVar)%var => &
                  flowDomsd(nn, 1, sps)%rlv(:, :, :)
          end if

          if( commEddyVis ) then 
             nVar = nVar + 1
             flowDoms(nn, level, sps)%realCommVars(nVar)%var => &
                  flowDomsd(nn, level, sps)%rev(:, :, :)
          end if

          if(nVar == 0) return
       end do domainLoop

       ! Run the generic exchange
       call wHalo1to1RealGeneric(nVar, level, sps, commPattern, internal)

    end do spectralModes

  end subroutine whalo1to1_d

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
             call computeETotBlock_b(2, il, 2, jl, 2, kl, correctForK)
          end do
       enddo domains
    endif bothPAndE

    mm = ubound(commPatternOverset, 1)
    call wOverset_b(level, start, end, commPressure, commVarGamma, &
         commLamVis, commEddyVis, commPatternOverset, internalOverset, mm)

    ! Exchange the 1 to 1 matching 2nd level cell halo's.

    call whalo1to1_b(level, start, end, commPressure, commVarGamma, &
         commLamVis, commEddyVis, commPatternCell_2nd,  &
         internalCell_2nd)

    ! NOTE: Only the 1to1 halo exchange and overset is done. whalosliding,
    ! whalomixing, orphanAverage and PandE corrections
    ! calculation are NOT implementent. 

  end subroutine whalo2_b
#endif
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
    mm = ubound(commPatternOverset, 1)
    call wOverset_d(level, start, end, commPressure, commVarGamma, &
         commLamVis, commEddyVis, commPatternOverset, internalOverset, mm)

    ! NOTE: Only the 1to1 halo and wOverset exchange is done. whalosliding,
    ! whalomixing, orphanAverage and PandE corrections
    ! calculation are NOT implementent. 

  end subroutine whalo2_d

  subroutine wOverset(level, start, end, commPressure,       &
       commVarGamma, commLamVis, commEddyVis, &
       commPattern, internal, nlev)
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

             weight => commPattern(level,mm)%sendList(i)%interp(j,:)

             ! Copy the given range of the working variables for
             ! this cell in the buffer. Update the counter jj.

             do k=start,end
                sendBuffer(jj) = &
                     weight(1)*flowDoms(d1,level,mm)%w(i1  ,j1  ,k1  ,k) + &
                     weight(2)*flowDoms(d1,level,mm)%w(i1+1,j1  ,k1  ,k) + &
                     weight(3)*flowDoms(d1,level,mm)%w(i1  ,j1+1,k1  ,k) + &
                     weight(4)*flowDoms(d1,level,mm)%w(i1+1,j1+1,k1  ,k) + &
                     weight(5)*flowDoms(d1,level,mm)%w(i1  ,j1  ,k1+1,k) + &
                     weight(6)*flowDoms(d1,level,mm)%w(i1+1,j1  ,k1+1,k) + &
                     weight(7)*flowDoms(d1,level,mm)%w(i1  ,j1+1,k1+1,k) + &
                     weight(8)*flowDoms(d1,level,mm)%w(i1+1,j1+1,k1+1,k)
                jj = jj + 1
             enddo

             ! The pressure, if needed.

             if( commPressure ) then
                sendBuffer(jj) = &
                     weight(1)*flowDoms(d1,level,mm)%p(i1  ,j1  ,k1  ) + &
                     weight(2)*flowDoms(d1,level,mm)%p(i1+1,j1  ,k1  ) + &
                     weight(3)*flowDoms(d1,level,mm)%p(i1  ,j1+1,k1  ) + &
                     weight(4)*flowDoms(d1,level,mm)%p(i1+1,j1+1,k1  ) + &
                     weight(5)*flowDoms(d1,level,mm)%p(i1  ,j1  ,k1+1) + &
                     weight(6)*flowDoms(d1,level,mm)%p(i1+1,j1  ,k1+1) + &
                     weight(7)*flowDoms(d1,level,mm)%p(i1  ,j1+1,k1+1) + &
                     weight(8)*flowDoms(d1,level,mm)%p(i1+1,j1+1,k1+1)
                jj = jj + 1
             endif

             ! The specific heat ratio, if needed. Note that level == 1.

             if( commVarGamma ) then
                sendBuffer(jj) = &
                     weight(1)*flowDoms(d1,1,mm)%gamma(i1  ,j1  ,k1  ) + &
                     weight(2)*flowDoms(d1,1,mm)%gamma(i1+1,j1  ,k1  ) + &
                     weight(3)*flowDoms(d1,1,mm)%gamma(i1  ,j1+1,k1  ) + &
                     weight(4)*flowDoms(d1,1,mm)%gamma(i1+1,j1+1,k1  ) + &
                     weight(5)*flowDoms(d1,1,mm)%gamma(i1  ,j1  ,k1+1) + &
                     weight(6)*flowDoms(d1,1,mm)%gamma(i1+1,j1  ,k1+1) + &
                     weight(7)*flowDoms(d1,1,mm)%gamma(i1  ,j1+1,k1+1) + &
                     weight(8)*flowDoms(d1,1,mm)%gamma(i1+1,j1+1,k1+1)
                jj = jj + 1
             endif

             ! The laminar viscosity for a viscous computation.
             ! Again level == 1.

             if( commLamVis ) then
                sendBuffer(jj) = &
                     weight(1)*flowDoms(d1,1,mm)%rlv(i1  ,j1  ,k1  ) + &
                     weight(2)*flowDoms(d1,1,mm)%rlv(i1+1,j1  ,k1  ) + &
                     weight(3)*flowDoms(d1,1,mm)%rlv(i1  ,j1+1,k1  ) + &
                     weight(4)*flowDoms(d1,1,mm)%rlv(i1+1,j1+1,k1  ) + &
                     weight(5)*flowDoms(d1,1,mm)%rlv(i1  ,j1  ,k1+1) + &
                     weight(6)*flowDoms(d1,1,mm)%rlv(i1+1,j1  ,k1+1) + &
                     weight(7)*flowDoms(d1,1,mm)%rlv(i1  ,j1+1,k1+1) + &
                     weight(8)*flowDoms(d1,1,mm)%rlv(i1+1,j1+1,k1+1)
                jj = jj + 1
             endif

             ! The eddy viscosity for eddy viscosity models.
             ! Level is the true multigrid level, because the eddy
             ! viscosity is allocated on all grid levels.

             if( commEddyVis ) then
                sendBuffer(jj) = &
                     weight(1)*flowDoms(d1,level,mm)%rev(i1  ,j1  ,k1  ) + &
                     weight(2)*flowDoms(d1,level,mm)%rev(i1+1,j1  ,k1  ) + &
                     weight(3)*flowDoms(d1,level,mm)%rev(i1  ,j1+1,k1  ) + &
                     weight(4)*flowDoms(d1,level,mm)%rev(i1+1,j1+1,k1  ) + &
                     weight(5)*flowDoms(d1,level,mm)%rev(i1  ,j1  ,k1+1) + &
                     weight(6)*flowDoms(d1,level,mm)%rev(i1+1,j1  ,k1+1) + &
                     weight(7)*flowDoms(d1,level,mm)%rev(i1  ,j1+1,k1+1) + &
                     weight(8)*flowDoms(d1,level,mm)%rev(i1+1,j1+1,k1+1)
                jj = jj + 1
             endif

          enddo

          ! Send the data.

          call mpi_isend(sendBuffer(ii), size, adflow_real, procId,  &
               procId, ADflow_comm_world, sendRequests(i), &
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
          size    = nVar*commPattern(level,mm)%nrecv(i)

          ! Post the receive.

          call mpi_irecv(recvBuffer(ii), size, adflow_real, procId, &
               myId, ADflow_comm_world, recvRequests(i), ierr)

          ! And update ii.

          ii = ii + size

       enddo receives

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

          ! Copy the given range of working variables.

          do k=start,end
             flowDoms(d2,level,mm)%w(i2,j2,k2,k) = &
                  weight(1)*flowDoms(d1,level,mm)%w(i1  ,j1  ,k1  ,k) + &
                  weight(2)*flowDoms(d1,level,mm)%w(i1+1,j1  ,k1  ,k) + &
                  weight(3)*flowDoms(d1,level,mm)%w(i1  ,j1+1,k1  ,k) + &
                  weight(4)*flowDoms(d1,level,mm)%w(i1+1,j1+1,k1  ,k) + &
                  weight(5)*flowDoms(d1,level,mm)%w(i1  ,j1  ,k1+1,k) + &
                  weight(6)*flowDoms(d1,level,mm)%w(i1+1,j1  ,k1+1,k) + &
                  weight(7)*flowDoms(d1,level,mm)%w(i1  ,j1+1,k1+1,k) + &
                  weight(8)*flowDoms(d1,level,mm)%w(i1+1,j1+1,k1+1,k)
          enddo

          ! The pressure, if needed.

          if( commPressure ) then
             flowDoms(d2,level,mm)%p(i2,j2,k2) = &
                  weight(1)*flowDoms(d1,level,mm)%p(i1  ,j1  ,k1  ) + &
                  weight(2)*flowDoms(d1,level,mm)%p(i1+1,j1  ,k1  ) + &
                  weight(3)*flowDoms(d1,level,mm)%p(i1  ,j1+1,k1  ) + &
                  weight(4)*flowDoms(d1,level,mm)%p(i1+1,j1+1,k1  ) + &
                  weight(5)*flowDoms(d1,level,mm)%p(i1  ,j1  ,k1+1) + &
                  weight(6)*flowDoms(d1,level,mm)%p(i1+1,j1  ,k1+1) + &
                  weight(7)*flowDoms(d1,level,mm)%p(i1  ,j1+1,k1+1) + &
                  weight(8)*flowDoms(d1,level,mm)%p(i1+1,j1+1,k1+1)
          end if

          ! The specific heat ratio, if needed. Note that level == 1.

          if( commVarGamma ) then
             flowDoms(d2,1,mm)%gamma(i2,j2,k2) = &
                  weight(1)*flowDoms(d1,1,mm)%gamma(i1  ,j1  ,k1  ) + &
                  weight(2)*flowDoms(d1,1,mm)%gamma(i1+1,j1  ,k1  ) + &
                  weight(3)*flowDoms(d1,1,mm)%gamma(i1  ,j1+1,k1  ) + &
                  weight(4)*flowDoms(d1,1,mm)%gamma(i1+1,j1+1,k1  ) + &
                  weight(5)*flowDoms(d1,1,mm)%gamma(i1  ,j1  ,k1+1) + &
                  weight(6)*flowDoms(d1,1,mm)%gamma(i1+1,j1  ,k1+1) + &
                  weight(7)*flowDoms(d1,1,mm)%gamma(i1  ,j1+1,k1+1) + &
                  weight(8)*flowDoms(d1,1,mm)%gamma(i1+1,j1+1,k1+1)
          end if

          ! The laminar viscosity for viscous computations.
          ! Again level == 1.

          if( commLamVis ) then
             flowDoms(d2,1,mm)%rlv(i2,j2,k2) = &
                  weight(1)*flowDoms(d1,1,mm)%rlv(i1  ,j1  ,k1  ) + &
                  weight(2)*flowDoms(d1,1,mm)%rlv(i1+1,j1  ,k1  ) + &
                  weight(3)*flowDoms(d1,1,mm)%rlv(i1  ,j1+1,k1  ) + &
                  weight(4)*flowDoms(d1,1,mm)%rlv(i1+1,j1+1,k1  ) + &
                  weight(5)*flowDoms(d1,1,mm)%rlv(i1  ,j1  ,k1+1) + &
                  weight(6)*flowDoms(d1,1,mm)%rlv(i1+1,j1  ,k1+1) + &
                  weight(7)*flowDoms(d1,1,mm)%rlv(i1  ,j1+1,k1+1) + &
                  weight(8)*flowDoms(d1,1,mm)%rlv(i1+1,j1+1,k1+1)
          end if

          ! The eddy viscosity for eddy viscosity models.
          ! Level is the true multigrid level, because the eddy
          ! viscosity is allocated on all grid levels.

          if( commEddyVis ) then
             flowDoms(d2,level,mm)%rev(i2,j2,k2) = &
                  weight(1)*flowDoms(d1,level,mm)%rev(i1  ,j1  ,k1  ) + &
                  weight(2)*flowDoms(d1,level,mm)%rev(i1+1,j1  ,k1  ) + &
                  weight(3)*flowDoms(d1,level,mm)%rev(i1  ,j1+1,k1  ) + &
                  weight(4)*flowDoms(d1,level,mm)%rev(i1+1,j1+1,k1  ) + &
                  weight(5)*flowDoms(d1,level,mm)%rev(i1  ,j1  ,k1+1) + &
                  weight(6)*flowDoms(d1,level,mm)%rev(i1+1,j1  ,k1+1) + &
                  weight(7)*flowDoms(d1,level,mm)%rev(i1  ,j1+1,k1+1) + &
                  weight(8)*flowDoms(d1,level,mm)%rev(i1+1,j1+1,k1+1)
          end if

       enddo localInterp

       ! Complete the nonblocking receives in an arbitrary sequence and
       ! copy the variables from the buffer into the halo's.

       size = commPattern(level,mm)%nProcRecv
       completeRecvs: do i=1,commPattern(level,mm)%nProcRecv

          ! Complete any of the requests.

          call mpi_waitany(size, recvRequests, index, status, ierr)

          ! Copy the data just arrived in the halo's.

          ii = index
          jj = nVar*commPattern(level,mm)%nrecvCum(ii-1)
          do j=1,commPattern(level,mm)%nrecv(ii)

             ! Store the block and the indices of the halo a bit easier.

             d2 = commPattern(level,mm)%recvList(ii)%block(j)
             i2 = commPattern(level,mm)%recvList(ii)%indices(j,1)
             j2 = commPattern(level,mm)%recvList(ii)%indices(j,2)
             k2 = commPattern(level,mm)%recvList(ii)%indices(j,3)

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

       ! Complete the nonblocking sends.

       size = commPattern(level,mm)%nProcSend
       do i=1,commPattern(level,mm)%nProcSend
          call mpi_waitany(size, sendRequests, index, status, ierr)
       enddo

    enddo spectralModes

  end subroutine wOverset
  subroutine wOverset_b(level, start, end, commPressure,       &
       commVarGamma, commLamVis, commEddyVis, &
       commPattern, internal, nlev)
    !
    !       wOverset_b performs the *TRANSPOSE* operation of wOveset       
    !       It is used for adjoint/reverse mode residual evaluations.      
    !      * See wOverset  for more information.
    !
    use constants
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
          call mpi_isend(recvBuffer(ii), size, adflow_real, procID,  &
               procID, ADflow_comm_world, recvRequests(i), &
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

          call mpi_irecv(sendBuffer(ii), size, adflow_real, procId, &
               myId, ADflow_comm_world, sendRequests(i), ierr)

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

  subroutine wOverset_d(level, start, end, commPressure,       &
       commVarGamma, commLamVis, commEddyVis, &
       commPattern, internal, nlev)
    !
    !       Modification to wOverset that communicates derivative values   
    !       in forward mode.                                               
    !
    use constants
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

             weight => commPattern(level,mm)%sendList(i)%interp(j,:)

             ! Copy the given range of the working variables for
             ! this cell in the buffer. Update the counter jj.

             do k=start,end
                sendBuffer(jj) = &
                     weight(1)*flowDomsd(d1,level,mm)%w(i1  ,j1  ,k1  ,k) + &
                     weight(2)*flowDomsd(d1,level,mm)%w(i1+1,j1  ,k1  ,k) + &
                     weight(3)*flowDomsd(d1,level,mm)%w(i1  ,j1+1,k1  ,k) + &
                     weight(4)*flowDomsd(d1,level,mm)%w(i1+1,j1+1,k1  ,k) + &
                     weight(5)*flowDomsd(d1,level,mm)%w(i1  ,j1  ,k1+1,k) + &
                     weight(6)*flowDomsd(d1,level,mm)%w(i1+1,j1  ,k1+1,k) + &
                     weight(7)*flowDomsd(d1,level,mm)%w(i1  ,j1+1,k1+1,k) + &
                     weight(8)*flowDomsd(d1,level,mm)%w(i1+1,j1+1,k1+1,k)
                jj = jj + 1
             enddo

             ! The pressure, if needed.

             if( commPressure ) then
                sendBuffer(jj) = &
                     weight(1)*flowDomsd(d1,level,mm)%p(i1  ,j1  ,k1  ) + &
                     weight(2)*flowDomsd(d1,level,mm)%p(i1+1,j1  ,k1  ) + &
                     weight(3)*flowDomsd(d1,level,mm)%p(i1  ,j1+1,k1  ) + &
                     weight(4)*flowDomsd(d1,level,mm)%p(i1+1,j1+1,k1  ) + &
                     weight(5)*flowDomsd(d1,level,mm)%p(i1  ,j1  ,k1+1) + &
                     weight(6)*flowDomsd(d1,level,mm)%p(i1+1,j1  ,k1+1) + &
                     weight(7)*flowDomsd(d1,level,mm)%p(i1  ,j1+1,k1+1) + &
                     weight(8)*flowDomsd(d1,level,mm)%p(i1+1,j1+1,k1+1)
                jj = jj + 1
             endif

             ! The specific heat ratio, if needed. Note that level == 1.

             if( commVarGamma ) then
                sendBuffer(jj) = &
                     weight(1)*flowDomsd(d1,1,mm)%gamma(i1  ,j1  ,k1  ) + &
                     weight(2)*flowDomsd(d1,1,mm)%gamma(i1+1,j1  ,k1  ) + &
                     weight(3)*flowDomsd(d1,1,mm)%gamma(i1  ,j1+1,k1  ) + &
                     weight(4)*flowDomsd(d1,1,mm)%gamma(i1+1,j1+1,k1  ) + &
                     weight(5)*flowDomsd(d1,1,mm)%gamma(i1  ,j1  ,k1+1) + &
                     weight(6)*flowDomsd(d1,1,mm)%gamma(i1+1,j1  ,k1+1) + &
                     weight(7)*flowDomsd(d1,1,mm)%gamma(i1  ,j1+1,k1+1) + &
                     weight(8)*flowDomsd(d1,1,mm)%gamma(i1+1,j1+1,k1+1)
                jj = jj + 1
             endif

             ! The laminar viscosity for a viscous computation.
             ! Again level == 1.

             if( commLamVis ) then
                sendBuffer(jj) = &
                     weight(1)*flowDomsd(d1,1,mm)%rlv(i1  ,j1  ,k1  ) + &
                     weight(2)*flowDomsd(d1,1,mm)%rlv(i1+1,j1  ,k1  ) + &
                     weight(3)*flowDomsd(d1,1,mm)%rlv(i1  ,j1+1,k1  ) + &
                     weight(4)*flowDomsd(d1,1,mm)%rlv(i1+1,j1+1,k1  ) + &
                     weight(5)*flowDomsd(d1,1,mm)%rlv(i1  ,j1  ,k1+1) + &
                     weight(6)*flowDomsd(d1,1,mm)%rlv(i1+1,j1  ,k1+1) + &
                     weight(7)*flowDomsd(d1,1,mm)%rlv(i1  ,j1+1,k1+1) + &
                     weight(8)*flowDomsd(d1,1,mm)%rlv(i1+1,j1+1,k1+1)
                jj = jj + 1
             endif

             ! The eddy viscosity for eddy viscosity models.
             ! Level is the true multigrid level, because the eddy
             ! viscosity is allocated on all grid levels.

             if( commEddyVis ) then
                sendBuffer(jj) = &
                     weight(1)*flowDomsd(d1,level,mm)%rev(i1  ,j1  ,k1  ) + &
                     weight(2)*flowDomsd(d1,level,mm)%rev(i1+1,j1  ,k1  ) + &
                     weight(3)*flowDomsd(d1,level,mm)%rev(i1  ,j1+1,k1  ) + &
                     weight(4)*flowDomsd(d1,level,mm)%rev(i1+1,j1+1,k1  ) + &
                     weight(5)*flowDomsd(d1,level,mm)%rev(i1  ,j1  ,k1+1) + &
                     weight(6)*flowDomsd(d1,level,mm)%rev(i1+1,j1  ,k1+1) + &
                     weight(7)*flowDomsd(d1,level,mm)%rev(i1  ,j1+1,k1+1) + &
                     weight(8)*flowDomsd(d1,level,mm)%rev(i1+1,j1+1,k1+1)
                jj = jj + 1
             endif

          enddo

          ! Send the data.

          call mpi_isend(sendBuffer(ii), size, adflow_real, procId,  &
               procId, ADflow_comm_world, sendRequests(i), &
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
          size    = nVar*commPattern(level,mm)%nrecv(i)

          ! Post the receive.

          call mpi_irecv(recvBuffer(ii), size, adflow_real, procId, &
               myId, ADflow_comm_world, recvRequests(i), ierr)

          ! And update ii.

          ii = ii + size

       enddo receives

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

          ! Copy the given range of working variables.

          do k=start,end
             flowDomsd(d2,level,mm)%w(i2,j2,k2,k) = &
                  weight(1)*flowDomsd(d1,level,mm)%w(i1  ,j1  ,k1  ,k) + &
                  weight(2)*flowDomsd(d1,level,mm)%w(i1+1,j1  ,k1  ,k) + &
                  weight(3)*flowDomsd(d1,level,mm)%w(i1  ,j1+1,k1  ,k) + &
                  weight(4)*flowDomsd(d1,level,mm)%w(i1+1,j1+1,k1  ,k) + &
                  weight(5)*flowDomsd(d1,level,mm)%w(i1  ,j1  ,k1+1,k) + &
                  weight(6)*flowDomsd(d1,level,mm)%w(i1+1,j1  ,k1+1,k) + &
                  weight(7)*flowDomsd(d1,level,mm)%w(i1  ,j1+1,k1+1,k) + &
                  weight(8)*flowDomsd(d1,level,mm)%w(i1+1,j1+1,k1+1,k)
          enddo

          ! The pressure, if needed.

          if( commPressure ) then
             flowDomsd(d2,level,mm)%p(i2,j2,k2) = &
                  weight(1)*flowDomsd(d1,level,mm)%p(i1  ,j1  ,k1  ) + &
                  weight(2)*flowDomsd(d1,level,mm)%p(i1+1,j1  ,k1  ) + &
                  weight(3)*flowDomsd(d1,level,mm)%p(i1  ,j1+1,k1  ) + &
                  weight(4)*flowDomsd(d1,level,mm)%p(i1+1,j1+1,k1  ) + &
                  weight(5)*flowDomsd(d1,level,mm)%p(i1  ,j1  ,k1+1) + &
                  weight(6)*flowDomsd(d1,level,mm)%p(i1+1,j1  ,k1+1) + &
                  weight(7)*flowDomsd(d1,level,mm)%p(i1  ,j1+1,k1+1) + &
                  weight(8)*flowDomsd(d1,level,mm)%p(i1+1,j1+1,k1+1)
          end if

          ! The specific heat ratio, if needed. Note that level == 1.

          if( commVarGamma ) then
             flowDomsd(d2,1,mm)%gamma(i2,j2,k2) = &
                  weight(1)*flowDomsd(d1,1,mm)%gamma(i1  ,j1  ,k1  ) + &
                  weight(2)*flowDomsd(d1,1,mm)%gamma(i1+1,j1  ,k1  ) + &
                  weight(3)*flowDomsd(d1,1,mm)%gamma(i1  ,j1+1,k1  ) + &
                  weight(4)*flowDomsd(d1,1,mm)%gamma(i1+1,j1+1,k1  ) + &
                  weight(5)*flowDomsd(d1,1,mm)%gamma(i1  ,j1  ,k1+1) + &
                  weight(6)*flowDomsd(d1,1,mm)%gamma(i1+1,j1  ,k1+1) + &
                  weight(7)*flowDomsd(d1,1,mm)%gamma(i1  ,j1+1,k1+1) + &
                  weight(8)*flowDomsd(d1,1,mm)%gamma(i1+1,j1+1,k1+1)
          end if

          ! The laminar viscosity for viscous computations.
          ! Again level == 1.

          if( commLamVis ) then
             flowDomsd(d2,1,mm)%rlv(i2,j2,k2) = &
                  weight(1)*flowDomsd(d1,1,mm)%rlv(i1  ,j1  ,k1  ) + &
                  weight(2)*flowDomsd(d1,1,mm)%rlv(i1+1,j1  ,k1  ) + &
                  weight(3)*flowDomsd(d1,1,mm)%rlv(i1  ,j1+1,k1  ) + &
                  weight(4)*flowDomsd(d1,1,mm)%rlv(i1+1,j1+1,k1  ) + &
                  weight(5)*flowDomsd(d1,1,mm)%rlv(i1  ,j1  ,k1+1) + &
                  weight(6)*flowDomsd(d1,1,mm)%rlv(i1+1,j1  ,k1+1) + &
                  weight(7)*flowDomsd(d1,1,mm)%rlv(i1  ,j1+1,k1+1) + &
                  weight(8)*flowDomsd(d1,1,mm)%rlv(i1+1,j1+1,k1+1)
          end if

          ! The eddy viscosity for eddy viscosity models.
          ! Level is the true multigrid level, because the eddy
          ! viscosity is allocated on all grid levels.

          if( commEddyVis ) then
             flowDomsd(d2,level,mm)%rev(i2,j2,k2) = &
                  weight(1)*flowDomsd(d1,level,mm)%rev(i1  ,j1  ,k1  ) + &
                  weight(2)*flowDomsd(d1,level,mm)%rev(i1+1,j1  ,k1  ) + &
                  weight(3)*flowDomsd(d1,level,mm)%rev(i1  ,j1+1,k1  ) + &
                  weight(4)*flowDomsd(d1,level,mm)%rev(i1+1,j1+1,k1  ) + &
                  weight(5)*flowDomsd(d1,level,mm)%rev(i1  ,j1  ,k1+1) + &
                  weight(6)*flowDomsd(d1,level,mm)%rev(i1+1,j1  ,k1+1) + &
                  weight(7)*flowDomsd(d1,level,mm)%rev(i1  ,j1+1,k1+1) + &
                  weight(8)*flowDomsd(d1,level,mm)%rev(i1+1,j1+1,k1+1)
          end if

       enddo localInterp

       ! Complete the nonblocking receives in an arbitrary sequence and
       ! copy the variables from the buffer into the halo's.

       size = commPattern(level,mm)%nProcRecv
       completeRecvs: do i=1,commPattern(level,mm)%nProcRecv

          ! Complete any of the requests.

          call mpi_waitany(size, recvRequests, index, status, ierr)

          ! Copy the data just arrived in the halo's.

          ii = index
          jj = nVar*commPattern(level,mm)%nrecvCum(ii-1)
          do j=1,commPattern(level,mm)%nrecv(ii)

             ! Store the block and the indices of the halo a bit easier.

             d2 = commPattern(level,mm)%recvList(ii)%block(j)
             i2 = commPattern(level,mm)%recvList(ii)%indices(j,1)
             j2 = commPattern(level,mm)%recvList(ii)%indices(j,2)
             k2 = commPattern(level,mm)%recvList(ii)%indices(j,3)

             ! Copy the conservative variables.

             do k=start,end
                jj = jj + 1
                flowDomsd(d2,level,mm)%w(i2,j2,k2,k) = recvBuffer(jj)
             enddo

             ! The pressure, if needed.

             if( commPressure ) then
                jj = jj + 1
                flowDomsd(d2,level,mm)%p(i2,j2,k2) = recvBuffer(jj)
             endif

             ! The specific heat ratio, if needed. Note that level == 1.

             if( commVarGamma ) then
                jj = jj + 1
                flowDomsd(d2,1,mm)%gamma(i2,j2,k2) = recvBuffer(jj)
             endif

             ! The laminar viscosity for viscous computations.
             ! Again level == 1.

             if( commLamVis ) then
                jj = jj + 1
                flowDomsd(d2,1,mm)%rlv(i2,j2,k2) = recvBuffer(jj)
             endif

             ! The eddy viscosity ratio for eddy viscosity models.
             ! Level is the true multigrid level, because the eddy
             ! viscosity is allocated on all grid levels.

             if( commEddyVis ) then
                jj = jj + 1
                flowDomsd(d2,level,mm)%rev(i2,j2,k2) = recvBuffer(jj)
             endif

          enddo

       enddo completeRecvs

       ! Complete the nonblocking sends.

       size = commPattern(level,mm)%nProcSend
       do i=1,commPattern(level,mm)%nProcSend
          call mpi_waitany(size, sendRequests, index, status, ierr)
       enddo

    enddo spectralModes

  end subroutine wOverset_d

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
    use utils, only : setPointers
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
                !print *,'sendbuf',jj,dd1,level,s,ii1,jj1,kk1,k, shape(sendbuffer), shape(flowDoms),shape(flowDoms(dd1,level,sps)%dw)
                sendBuffer(jj) = flowDoms(dd1,level,sps)%dw(ii1,jj1,kk1,k)
                jj = jj + 1
             enddo

          enddo

          ! Send the data.

          call mpi_isend(sendBuffer(ii), size, adflow_real, procID,  &
               procID, ADflow_comm_world, sendRequests(i), &
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

          call mpi_irecv(recvBuffer(ii), size, adflow_real, procID, &
               myID, ADflow_comm_world, recvRequests(i), ierr)

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

  subroutine exchangeCoor(level)
    !
    !       ExchangeCoor exchanges the coordinates of the given grid       
    !       level.                                                         
    !
    use block
    use communication
    use inputTimeSpectral
    !
    !      Subroutine arguments.
    !
    integer(kind=intType), intent(in) :: level
    !
    !      Local variables.
    !
    integer :: size, procID, ierr, index
    integer, dimension(mpi_status_size) :: status

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

          ! And update ii.

          ii = ii + size

       enddo receives

       ! Copy the local data.

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

          call mpi_waitany(size, recvRequests, index, status, ierr)

          ! Copy the data just arrived in the halo's.

          ii = index
          jj = 3*commPatternNode_1st(level)%nRecvCum(ii-1) +1
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
          call mpi_waitany(size, sendRequests, index, status, ierr)
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
  implicit none
  !
  !      Subroutine arguments.
  !
  integer(kind=intType), intent(in) :: level
  !
  !      Local variables.
  !
  integer :: size, procID, ierr, index
  integer, dimension(mpi_status_size) :: status

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

        ! And update ii.

        ii = ii + size

     enddo send

     ! Copy the local data.

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

        call mpi_waitany(size, sendRequests, index, status, ierr)

        ! Copy the data just arrived in the halo's.

        ii = index
        jj = 3*commPatternNode_1st(level)%nSendCum(ii-1)
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
        call mpi_waitany(size, recvRequests, index, status, ierr)
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
  implicit none
  !
  !      Subroutine arguments.
  !
  integer(kind=intType), intent(in) :: level
  !
  !      Local variables.
  !
  integer :: size, procID, ierr, index
  integer, dimension(mpi_status_size) :: status

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

        ! And update ii.

        ii = ii + size

     enddo receives

     ! Copy the local data.

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

        call mpi_waitany(size, recvRequests, index, status, ierr)

        ! Copy the data just arrived in the halo's.

        ii = index
        jj = 3*commPatternNode_1st(level)%nRecvCum(ii-1) +1
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
        call mpi_waitany(size, sendRequests, index, status, ierr)
     enddo

  enddo spectralLoop

end subroutine exchangeCoor_d

end module haloExchange

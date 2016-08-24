module haloExchange 

contains

  subroutine whalo1(level, start, end, commPressure, commGamma, &
       commViscous)
    !
    !      ******************************************************************
    !      *                                                                *
    !      * whalo1 exchanges all the 1st level internal halo's for the     *
    !      * cell centered variables.                                       *
    !      *                                                                *
    !      ******************************************************************
    !
    use BCTypes
    use blockPointers
    use commSliding
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

    ! Exchange the sliding mesh 1st level cell halo's.

    mm = ubound(commSlidingCell_1st,1)
    call whaloSliding(level, start, end, commPressure,       &
         commVarGamma, commLamVis, commEddyVis, &
         commSlidingCell_1st, intSlidingCell_1st, mm)

    ! Exchange the mixing plane 1st level cell halo's.

    call whaloMixing(level, start, end, commPressure, commVarGamma, &
         commLamVis, commEddyVis, 1_intType)


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
    !      ******************************************************************
    !      *                                                                *
    !      * whalo2 exchanges all the 2nd level internal halo's for the     *
    !      * cell centered variables.                                       *
    !      *                                                                *
    !      ******************************************************************
    !
    use BCTypes
    use blockPointers
    use commSliding
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

    ! Exchange the sliding mesh 2nd level cell halo's.

    mm = ubound(commSlidingCell_2nd,1)
    call whaloSliding(level, start, end, commPressure,       &
         commVarGamma, commLamVis, commEddyVis, &
         commSlidingCell_2nd, intSlidingCell_2nd, mm)

    ! Exchange the mixing plane 2nd level cell halo's.

    call whaloMixing(level, start, end, commPressure, commVarGamma, &
         commLamVis, commEddyVis, 2_intType)

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

       ! Loop over the blocks to find the sliding mesh subfaces.
       ! Use is made of the fact the boundary conditions are identical
       ! for all spectral solutions. So that loop can be inside the
       ! test for the sliding mesh subface.

       domains: do nn=1,nDom
          do mm=1,flowDoms(nn,level,1)%nBocos
             if(flowDoms(nn,level,1)%BCType(mm) == slidingInterface) then

                ! Loop over the number of spectral solutions.

                do ll=1,nTimeIntervalsSpectral

                   ! Set the pointers for this block.

                   call setPointers(nn,level,ll)

                   ! Set the range, depending on the block face on which
                   ! the subface is located.

                   select case (bcFaceID(mm))
                   case (iMin)
                      iBeg = 0;         jBeg = jcBeg(mm); kBeg = kcBeg(mm)
                      iEnd = 1;         jEnd = jcEnd(mm); kEnd = kcEnd(mm)
                   case (iMax)
                      iBeg = ie;        jBeg = jcBeg(mm); kBeg = kcBeg(mm)
                      iEnd = ib;        jEnd = jcEnd(mm); kEnd = kcEnd(mm)
                   case (jMin)
                      iBeg = icBeg(mm); jBeg = 0;         kBeg = kcBeg(mm)
                      iEnd = icEnd(mm); jEnd = 1;         kEnd = kcEnd(mm)
                   case (jMax)
                      iBeg = icBeg(mm); jBeg = je;        kBeg = kcBeg(mm)
                      iEnd = icEnd(mm); jEnd = jb;        kEnd = kcEnd(mm)
                   case (kMin)
                      iBeg = icBeg(mm); jBeg = jcBeg(mm); kBeg = 0
                      iEnd = icEnd(mm); jEnd = jcEnd(mm); kEnd = 1
                   case (kMax)
                      iBeg = icBeg(mm); jBeg = jcBeg(mm); kBeg = ke
                      iEnd = icEnd(mm); jEnd = jcEnd(mm); kEnd = kb
                   end select

                   ! Compute the total energy for the sliding mesh halo's

                   call computeEtotBlock(iBeg, iEnd, jBeg, jEnd, kBeg, kEnd, &
                        correctForK)
                enddo
             endif

             ! Treat the overset blocks. Since we don't have the logic
             ! setup here correctly to only update the overset cells,
             ! just do the whole block, for every block
             do ll=1, nTimeIntervalsSpectral
                call setPointers(nn, level, ll)
                call computeETotBlock(2, il, 2, jl, 2, kl, correctForK)
             end do

          enddo

       enddo domains

    endif bothPAndE

  end subroutine whalo2

  subroutine orphanAverage(wstart, wend, calcPressure, calcGamma, &
       calcLamVis, calcEddyVis)
    !
    !      ******************************************************************
    !      *                                                                *
    !      * orphanAverage uses the neighboring cells of an overset orphan  *
    !      * to set the flow state for the orphan cell by a simple average. *
    !      * This routine operates on the block given by the block pointers *
    !      * so it is assumed they are set.                                 *
    !      *                                                                *
    !      ******************************************************************
    !
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
    !      ******************************************************************
    !      *                                                                *
    !      * whalo1to1 exchanges the 1 to 1 internal halo's for the cell    *
    !      * centered variables for the given communication pattern.        *
    !      * Pointers must be set for var1, var2...varN                     *
    !      *                                                                *
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
    !      ******************************************************************
    !      *                                                                *
    !      * whalo1to1 exchanges the 1 to 1 internal halo's for the cell    *
    !      * centered variables for the given communication pattern.        *
    !      * Pointers must be set for var1, var2...varN                     *
    !      *                                                                *
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

       call mpi_isend(sendBufInt(ii), size, sumb_integer, procID,  &
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

       call mpi_irecv(recvBufInt(ii), size, sumb_integer, procID, &
            myID, SUmb_comm_world, recvRequests(i), ierr)

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
    !      ******************************************************************
    !      *                                                                *
    !      * whalo1to1 exchanges the 1 to 1 internal halo's derivatives     *
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

  !
  !      ******************************************************************
  !      *                                                                *
  !      * File:          whaloMixing.f90                                 *
  !      * Author:        Edwin van der Weide                             *
  !      * Starting date: 01-26-2005                                      *
  !      * Last modified: 09-16-2005                                      *
  !      *                                                                *
  !      ******************************************************************
  !
  subroutine whaloMixing(level, startIn, endIn, commPressure,   &
       commVarGamma, commLamVis, commEddyVis, &
       nLayers)
    !
    !      ******************************************************************
    !      *                                                                *
    !      * whaloMixing determines the halo values for the cells adjacent  *
    !      * to sliding interfaces for the given communication pattern.     *
    !      * The mixing plane approximation is just a azimuthal averaging   *
    !      * of the variables just that a steady flow can be computed for   *
    !      * the inherently unsteady problem of rotor/stator interaction.   *
    !      * It is possible to exchange a range of variables and not the    *
    !      * entire set, e.g. only the flow variables or only the turbulent *
    !      * variables. This is controlled by the arguments start, end,     *
    !      * commPressure and commViscous. The exchange takes place for     *
    !      * the given grid level.                                          *
    !      *                                                                *
    !      ******************************************************************
    !
    use communication
    use commMixing
    use constants
    use inputPhysics
    use inputTimeSpectral
    use interfaceGroups
    implicit none
    !
    !      Subroutine arguments.
    !
    integer(kind=intType), intent(in) :: level, startIn, endIn
    integer(kind=intType), intent(in) :: nLayers

    logical, intent(in) :: commPressure, commVarGamma
    logical, intent(in) :: commLamVis, commEddyVis
    !
    !      Local variables.
    !
    integer :: comm, size, ierr

    integer(kind=intType) :: nVar, start, end
    integer(kind=intType) :: sps, nn, mm, ll
    integer(kind=intType) :: nInter, nDonor, nHalo

    integer(kind=intType), dimension(:),   pointer :: bd, bh
    integer(kind=intType), dimension(:),   pointer :: indListD
    integer(kind=intType), dimension(:),   pointer :: nintD
    integer(kind=intType), dimension(:,:), pointer :: indListH
    integer(kind=intType), dimension(:,:), pointer :: indD, indH

    real(kind=realType), dimension(:),     pointer :: wd
    real(kind=realType), dimension(:,:),   pointer :: wh
    real(kind=realType), dimension(:,:,:), pointer :: matD, matH

    logical :: commVel

    ! Return immediately if either not a steady problem is solved or
    ! if no sliding mesh interfaces are present.

    if(equationMode /= steady .or. nInterfaceGroups == 0) return

    ! Set the logical whether or not the velocity components are
    ! communicated. If they have to be communicated, possibly correct
    ! the start and end indices.

    commVel = .false.
    if(startIn <= ivz .and. endIn >= ivx) commVel = .true.

    start = startIn;  end = endIn

    if( commVel ) then
       if(start > ivx) start = ivx
       if(end   < ivz) end   = ivz
    endif

    ! Determine the number of variables per cell to be set.

    nVar = max(0_intType,(end - start + 1))
    if( commPressure ) nVar = nVar + 1
    if( commVarGamma ) nVar = nVar + 1
    if( commLamVis )   nVar = nVar + 1
    if( commEddyVis )  nVar = nVar + 1

    if(nVar == 0) return

    ! Loop over the number of spectral solutions. For a steady state
    ! computation this will be one, but it is added for consistency.

    spectralLoop: do sps=1,nTimeIntervalsSpectral

       ! Loop over the number of interface groups.

       interfaceLoop: do nn=1,nInterfaceGroups

          ! If i do not contribute to the current "color" continue
          ! to the next color.

          if(.not. myInterfaces(nn)%procContributes) cycle

          ! Store the communicator a bit easier.

          comm = myInterfaces(nn)%commSlide

          ! Loop over the two sides of the interface.

          sideLoop: do ll=1,2

             ! Some abbreviations to make the code more readable.

             nInter = commPatternMixing(level,nn,ll)%nInter
             nDonor = commPatternMixing(level,nn,ll)%nDonor
             nHalo  = commPatternMixing(level,nn,ll)%nHalo

             wd       => commPatternMixing(level,nn,ll)%weightDonor
             matD     => commPatternMixing(level,nn,ll)%rotMatDonor
             bd       => commPatternMixing(level,nn,ll)%blockDonor
             nintD    => commPatternMixing(level,nn,ll)%nIntervalsDonor
             indListD => commPatternMixing(level,nn,ll)%indListDonor

             wh       => commPatternMixing(level,nn,ll)%weightHalo
             matH     => commPatternMixing(level,nn,ll)%rotMatHalo
             bh       => commPatternMixing(level,nn,ll)%blockHalo
             indListH => commPatternMixing(level,nn,ll)%indListHalo

             ! Loop over the number of halo layers to be set. This is
             ! either 1 or 2.

             layerLoop: do mm=1,nLayers

                ! Set the pointers for the donor and halo indices to make
                ! the code more readable.

                indD => commPatternMixing(level,nn,ll)%indD(:,:,mm)
                indH => commPatternMixing(level,nn,ll)%indH(:,:,mm)

                ! Determine the local contribution to interpolation
                ! list. This data is stored in sendBuffer.

                call localPartMixingPlane(sendBuffer)

                ! Determine the sum of the local contribution for the
                ! processor group of this sliding interface. Store the
                ! global sum in recvBuffer.

                size = nVar*nInter
                call mpi_allreduce(sendBuffer, recvBuffer, size, &
                     sumb_real, mpi_sum, comm, ierr)

                ! Determine the local mixing plane halo cells.

                call setMixingPlaneHalos(recvBuffer)

             enddo layerLoop
          enddo sideLoop
       enddo interfaceLoop
    enddo spectralLoop

  contains

    subroutine localPartMixingPlane(buffer)
      !
      !        ****************************************************************
      !        *                                                              *
      !        * LocalPartMixingPlane stores the local contribution to the    *
      !        * currently active mixing plane into buffer.                   *
      !        *                                                              *
      !        ****************************************************************
      !
      use block
      implicit none
      !
      !        Subroutine arguments.
      !
      real(kind=realType), dimension(nVar,*), intent(out) :: buffer
      !
      !        Local variables.
      !
      integer(kind=intType) :: i, j, k, b, ii, jj, nn, mm
      integer(kind=intType) :: iax, irad, itheta, iof, iof2, start2

      real(kind=realType) :: vx, vy, vz, vax, vrad, vtheta

      ! Initialize the buffer to zero; it will contain the sum of the
      ! local values.

      do i=1,nInter
         do j=1,nVar
            buffer(j,i) = zero
         enddo
      enddo

      ! First treat the velocities, because these require a
      ! transformation from the global cartesian to the local
      ! cylindrical frame.

      testCommVel: if( commVel ) then

         ! Determine the indices in the buffer to store the three
         ! cylindrical velocity components.

         iax    = ivx - start + 1
         irad   = iax + 1
         itheta = irad + 1

         ! Loop over the number of donors.

         velDonorLoop: do nn=1,nDonor

            ! Store the indices of the donor a bit easier.

            b = bd(nn)
            i = indD(nn,1)
            j = indD(nn,2)
            k = indD(nn,3)

            ! Store the 3 cartesian velocity components a bit easier.

            vx = flowDoms(b,level,sps)%w(i,j,k,ivx)
            vy = flowDoms(b,level,sps)%w(i,j,k,ivy)
            vz = flowDoms(b,level,sps)%w(i,j,k,ivz)

            ! Transform the velocity to the components of the local
            ! cylindrical coordinate system

            vax    = matD(nn,1,1)*vx + matD(nn,1,2)*vy &
                 + matD(nn,1,3)*vz
            vrad   = matD(nn,2,1)*vx + matD(nn,2,2)*vy &
                 + matD(nn,2,3)*vz
            vtheta = matD(nn,3,1)*vx + matD(nn,3,2)*vy &
                 + matD(nn,3,3)*vz

            ! Loop over the number of intervals to which this
            ! donor contributes.

            do jj=(nintD(nn-1)+1),nintD(nn)

               ! Store the index in the buffer a bit easier and
               ! update the velocity entries in the buffer.

               ii = indListD(jj)

               buffer(iax,ii)    = buffer(iax,ii)    + wd(jj)*vax
               buffer(irad,ii)   = buffer(irad,ii)   + wd(jj)*vrad
               buffer(itheta,ii) = buffer(itheta,ii) + wd(jj)*vtheta
            enddo

         enddo velDonorLoop
      endif testCommVel

      ! Determine the offset such that the index start will be stored
      ! at the first position in the buffer. Also set start2, needed
      ! to store the variables with an index larger than the velocity.
      ! Iof2 is the offset to store the other variables than the
      ! working ones at the correct location in buffer

      iof    = 1 - start
      iof2   = max(0_intType,(end - start + 1))
      start2 = max(start,ivz+1)

      ! Store the other variables in the buffer.

      donorLoop: do nn=1,nDonor

         ! Store the indices of the donor a bit easier.

         b = bd(nn)
         i = indD(nn,1)
         j = indD(nn,2)
         k = indD(nn,3)

         ! Loop over the number of intervals to which this
         ! donor contributes.

         contribLoop: do jj=(nintD(nn-1)+1),nintD(nn)

            ! Store the index in the list a bit easier.

            ii = indListD(jj)

            ! The working variables. Be careful not to overwrite the
            ! already computed velocity components.

            do mm=start,(ivx-1)
               buffer(iof+mm,ii) = buffer(iof+mm,ii) &
                    + wd(jj)*flowDoms(b,level,sps)%w(i,j,k,mm)
            enddo

            do mm=start2,end
               buffer(iof+mm,ii) = buffer(iof+mm,ii) &
                    + wd(jj)*flowDoms(b,level,sps)%w(i,j,k,mm)
            enddo

            ! The other variables. Note that for gamma and rlv the level
            ! is 1, because these variables are only allocated on the
            ! finest grid; they are not considered true mg variables in
            ! the sense that they depend on other variables.

            mm = iof2

            if( commPressure ) then
               mm = mm + 1
               buffer(mm,ii) = buffer(mm,ii) &
                    + wd(jj)*flowDoms(b,level,sps)%p(i,j,k)
            endif

            if( commVarGamma ) then
               mm = mm + 1
               buffer(mm,ii) = buffer(mm,ii) &
                    + wd(jj)*flowDoms(b,1,sps)%gamma(i,j,k)
            endif

            if( commLamVis ) then
               mm = mm + 1
               buffer(mm,ii) = buffer(mm,ii) &
                    + wd(jj)*flowDoms(b,1,sps)%rlv(i,j,k)
            endif

            if( commEddyVis ) then
               mm = mm + 1
               buffer(mm,ii) = buffer(mm,ii) &
                    + wd(jj)*flowDoms(b,level,sps)%rev(i,j,k)
            endif

         enddo contribLoop
      enddo donorLoop

    end subroutine localPartMixingPlane

    !===============================================================

    subroutine setMixingPlaneHalos(buffer)
      !
      !        ****************************************************************
      !        *                                                              *
      !        * SetMixingPlaneHalos set the values in the halo cells         *
      !        * adjacent to the active mixing plane. The variables are       *
      !        * interpolated from buffer.                                    *
      !        *                                                              *
      !        ****************************************************************
      !
      use block
      implicit none
      !
      !        Subroutine arguments.
      !
      real(kind=realType), dimension(nVar,*), intent(in) :: buffer
      !
      !        Local variables.
      !
      integer(kind=intType) :: i, j, k, b, nn, mm, i1, i2
      integer(kind=intType) :: iof, iof2

      real(kind=realType) :: w1, w2, vax, vrad, vtheta, vx, vy, vz

      ! Determine the offset iof such that the index start will get
      ! its data from index 1 in the buffer. Also set the offset iof2
      ! for the "other" variables.

      iof  = 1 - start
      iof2 = max(0_intType,(end - start + 1))

      ! Loop over the number of halos to be set.

      haloLoop: do nn=1,nHalo

         ! Store the indices of the halo as well as the indices and
         ! the weights in the list a bit easier.

         b = bh(nn)
         i = indH(nn,1)
         j = indH(nn,2)
         k = indH(nn,3)

         i1 = indListH(nn,1)
         i2 = indListH(nn,2)

         w1 = wh(nn,1)
         w2 = wh(nn,2)

         ! Loop over the working variables to be set. Note that the
         ! velocity components of the local cylindrical coordinate
         ! system are interpolated. These are transformed to the
         ! cartesian frame later on.

         do mm=start,end
            flowDoms(b,level,sps)%w(i,j,k,mm) = w1*buffer(iof+mm,i1) &
                 + w2*buffer(iof+mm,i2)
         enddo

         ! The other variables. Again level = 1 for gamma and rlv,
         ! see the explanation earlier.

         mm = iof2

         if( commPressure ) then
            mm = mm + 1
            flowDoms(b,level,sps)%p(i,j,k) = w1*buffer(mm,i1) &
                 + w2*buffer(mm,i2)
         endif

         if( commVarGamma ) then
            mm = mm + 1
            flowDoms(b,1,sps)%gamma(i,j,k) = w1*buffer(mm,i1) &
                 + w2*buffer(mm,i2)
         endif

         if( commLamVis ) then
            mm = mm + 1
            flowDoms(b,1,sps)%rlv(i,j,k) = w1*buffer(mm,i1) &
                 + w2*buffer(mm,i2)
         endif

         if( commEddyVis ) then
            mm = mm + 1
            flowDoms(b,level,sps)%rev(i,j,k) = w1*buffer(mm,i1) &
                 + w2*buffer(mm,i2)
         endif

      enddo haloLoop

      ! Transform the velocities to the cartesian components if
      ! velocities must be communicated.

      testCommVel: if( commVel ) then

         ! Loop over the number of halo cells to be set.

         velHaloLoop: do nn=1,nHalo

            ! Store the indices of the halo a bit easier.

            b = bh(nn)
            i = indH(nn,1)
            j = indH(nn,2)
            k = indH(nn,3)

            ! Compute the cartesian components of the velocity.

            vax    = flowDoms(b,level,sps)%w(i,j,k,ivx)
            vrad   = flowDoms(b,level,sps)%w(i,j,k,ivy)
            vtheta = flowDoms(b,level,sps)%w(i,j,k,ivz)

            vx = matH(nn,1,1)*vax    + matH(nn,1,2)*vrad &
                 + matH(nn,1,3)*vtheta
            vy = matH(nn,2,1)*vax    + matH(nn,2,2)*vrad &
                 + matH(nn,2,3)*vtheta
            vz = matH(nn,3,1)*vax    + matH(nn,3,2)*vrad &
                 + matH(nn,3,3)*vtheta

            flowDoms(b,level,sps)%w(i,j,k,ivx) = vx
            flowDoms(b,level,sps)%w(i,j,k,ivy) = vy
            flowDoms(b,level,sps)%w(i,j,k,ivz) = vz

         enddo velHaloLoop
      endif testCommVel

    end subroutine setMixingPlaneHalos

  end subroutine whaloMixing

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

  subroutine whalo2_b(level, start, end, commPressure, commGamma, &
       commViscous)
    !
    !      ******************************************************************
    !      *                                                                *
    !      * whalo2_b exchanges all the 2nd level internal halo's for the   *
    !      * cell centered variables IN REVERSE MODE                        *
    !      *                                                                *
    !      ******************************************************************
    !
    use BCTypes
    use blockPointers
    use commSliding
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


  subroutine whalo2_d(level, start, end, commPressure, commGamma, &
       commViscous)
    !
    !      ******************************************************************
    !      *                                                                *
    !      * whalo2_b exchanges all the 2nd level internal halo's for the   *
    !      * cell centered variables IN FORWARD MODE                        *
    !      *                                                                *
    !      ******************************************************************
    !
    use BCTypes
    use blockPointers
    use commSliding
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

  !
  !      ******************************************************************
  !      *                                                                *
  !      * File:          wOverset.f90                                    *
  !      * Author:        Steve Repsher, Gaetan Kenway                    *
  !      * Starting date: 02-04-2005                                      *
  !      * Last modified: 10-01-2015                                      *
  !      *                                                                *
  !      ******************************************************************
  !
  subroutine wOverset(level, start, end, commPressure,       &
       commVarGamma, commLamVis, commEddyVis, &
       commPattern, internal, nlev)
    !
    !      ******************************************************************
    !      *                                                                *
    !      * wOverset controls the communication between overset halos      *
    !      * for the cell-centered variables by interpolating the solution  *
    !      * from other blocks consistent with the chimera approach. A tri- *
    !      * linear interpolation is used as per the input from cgns. It    *
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

          call mpi_isend(sendBuffer(ii), size, sumb_real, procId,  &
               procId, SUmb_comm_world, sendRequests(i), &
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

          call mpi_irecv(recvBuffer(ii), size, sumb_real, procId, &
               myId, SUmb_comm_world, recvRequests(i), ierr)

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

  subroutine wOverset_d(level, start, end, commPressure,       &
       commVarGamma, commLamVis, commEddyVis, &
       commPattern, internal, nlev)
    !
    !      ******************************************************************
    !      *                                                                *
    !      * Modification to wOverset that communicates derivative values   *
    !      * in forward mode.                                               *
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

          call mpi_isend(sendBuffer(ii), size, sumb_real, procId,  &
               procId, SUmb_comm_world, sendRequests(i), &
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

          call mpi_irecv(recvBuffer(ii), size, sumb_real, procId, &
               myId, SUmb_comm_world, recvRequests(i), ierr)

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
                !print *,'sendbuf',jj,dd1,level,s,ii1,jj1,kk1,k, shape(sendbuffer), shape(flowDoms),shape(flowDoms(dd1,level,sps)%dw)
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


end module haloExchange

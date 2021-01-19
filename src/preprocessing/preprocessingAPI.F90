module preprocessingAPI

contains


  subroutine preprocessing
    !
    !       preprocessing determines the communication patterns between
    !       the processors for all the mg levels, computes the wall
    !       distances and the metrics.
    !
    use constants
    use block
    use blockPointers, only : BCType, nBocos
    use cgnsGrid
    use bcdata, only : initBCData, allocMemBCData
    use communication, only : adflow_comm_world, commPatternCell_1st, &
         commPatternCell_2nd, commPatternNode_1st, internalCell_1st, &
         internalCell_2nd, internalNode_1st, myid, nProc, &
         recvBufferSize_1to1, sendBufferSize_1to1, sendBufferSIzeOver,&
         recvBufferSizeOver, commPatternOverset, internalOverset, sendBuffer, &
         recvBuffer, sendBufferSize, recvBufferSize
    use inputPhysics
    use inputTimeSpectral
    use section
    use wallDistance, only : xVolumeVec, xSurfVec, wallScatter, &
         wallDistanceDataAllocated, updateWallAssociation, &
         computeWallDistance
    use oversetData, only : cumDomProc, nDomProc, wallFringes, nDomTotal, &
         overlapMatrix, oversetPresent, localWallFringes
    use utils, only : setPointers, EChk, setBufferSizes, terminate
    use coarseUtils, only : createCoarseBlocks
    use pointMatchedCommPattern, only : determineCommPattern
    use oversetAPI, only : oversetComm, determineClusters, determineViscousDirs
    implicit none
    !
    !      Local variables.
    !
    integer :: ierr

    integer(kind=intType) :: nLevels, level, nn, mm, nsMin, nsMax, i, iProc
    logical :: local
    !
    ! Check that for the unsteady modes the number of periodic slices
    ! is identical for all sections.

    nsMin = sections(1)%nSlices
    nsMax = sections(1)%nSlices

    do nn=2,nSections
       nsMin = min(nsMin,sections(nn)%nSlices)
       nsMax = max(nsMax,sections(nn)%nSlices)
    enddo

    if((equationMode == unsteady .or. &
         equationMode == timeSpectral) .and. nsMin < nsMax) then

       if(myID == 0)                     &
            call terminate("preprocessing", &
            "Different rotational periodicity encountered &
            &for time accurate computation")
       call mpi_barrier(ADflow_comm_world, ierr)

    endif

    ! Determine the number of multigrid levels needed in the
    ! computation, and allocate the memory for the node and cell
    ! communication patterns, including the memory copies for internal
    ! communication, and the global number of cells on each level.
    ! Note that this communication pattern does not change in time.

    nLevels = ubound(flowDoms,2)
    nn      = nLevels
    allocate(commPatternCell_1st(nn), commPatternCell_2nd(nn), &
         commPatternNode_1st(nn), internalCell_1st(nn),    &
         internalCell_2nd(nn),    internalNode_1st(nn),    &
         nCellGlobal(nn),         stat=ierr)
    if(ierr /= 0)                     &
         call terminate("preprocessing", &
         "Memory allocation failure for commPatterns")

    ! Set the sizes here so that we know how to dealloc the stuff
    ! later on
    do i=1,nLevels
       commPatternCell_1st(i)%nPeriodic = 0
       commPatternCell_1st(i)%nProcSend = 0
       commPatternCell_1st(i)%nProcRecv = 0

       commPatternCell_2nd(i)%nPeriodic = 0
       commPatternCell_2nd(i)%nProcSend = 0
       commPatternCell_2nd(i)%nProcRecv = 0

       commPatternNode_1st(i)%nPeriodic = 0
       commPatternNode_1st(i)%nProcSend = 0
       commPatternNode_1st(i)%nProcRecv = 0

       internalCell_1st(i)%nPeriodic = 0
       internalCell_2nd(i)%nPeriodic = 0
       internalNode_1st(i)%nPeriodic = 0
    end do

    ! Allocate the memory for the overset mesh communication pattern.
    ! This pattern changes in time and therefore each spectral time
    ! value has its own sliding mesh communication pattern.

    mm = nTimeIntervalsSpectral
    allocate(commPatternOverset(nn,mm), internalOverset(nn,mm), &
         overlapMatrix(nn, mm), stat=ierr)
    if(ierr /= 0)                     &
         call terminate("preprocessing", &
         "Memory allocation failure for commOverset")

    ! Determine the fine grid 1 to 1 matching communication pattern.

    call determineCommPattern(1_intType)

    ! Initialize the send and receive buffer sizes to 0 and determine
    ! its size for the finest grid for the 1 to 1 communication.

    sendBufferSize_1to1  = 0
    recvBufferSize_1to1  = 0
    sendBufferSizeOver   = 0
    recvBufferSizeOver   = 0

    call setBufferSizes(1_intType, 1_intType, .true., .false.)

    ! Loop to create the coarse grid levels.

    do level=2,nLevels

       ! Create the coarse grid blocks, its communication pattern, the
       ! coarse grid level 0 cooling parameters and check the
       ! communication buffer sizes.

       call createCoarseBlocks(level)
       call determineCommPattern(level)
       call setBufferSizes(level, 1_intType, .true., .false.)

    enddo

    ! Synchronize the processors, just to be sure.

    call mpi_barrier(ADflow_comm_world, ierr)

    ! Allocate memory for the nonblocking point to point communication.

    allocate(sendBuffer(sendBufferSize), &
         recvBuffer(recvBufferSize), stat=ierr)
    if(ierr /= 0)                     &
         call terminate("preprocessing", &
         "Memory allocation failure for sendBuffer &
         &and recvBuffer")

    ! Determine the cell range for the subfaces and initialize the
    ! arrays for the boundary condition data.
    ! Done for all grid levels.

    call cellRangeSubface
    call initBcdata

    ! Allocate some data of size nLevels for the fast wall distance calc
    allocate(xVolumeVec(nLevels), xSurfVec(nLevels, mm), wallScatter(nLevels, mm), &
         wallDistanceDataAllocated(nLevels), updateWallAssociation(nLevels))
    wallDistanceDataAllocated = .False.
    updateWallAssociation = .True.

    ! Nullify the wallFringe poiter as initialization
    nullify(wallFringes, localWallFringes)

    ! Allocate nDomProc: the number of domains on each processor
    ! and a cumulative form.
    allocate(nDomProc(0:nProc-1), cumDomProc(0:nProc))

    ! Gather the dimensions of all blocks to everyone
    call mpi_allreduce(nDom, nDomTotal, 1, adflow_integer, MPI_SUM, &
         adflow_comm_world, ierr)

    ! Receive the number of domains from each proc using an allgather.
    call mpi_allgather(nDom, 1, adflow_integer, nDomProc, 1, adflow_integer, &
         adflow_comm_world, ierr)

    ! Compute the cumulative format:
    cumDomProc(0) = 0
    do iProc=1, nProc
       cumDomProc(iProc) = cumDomProc(iProc-1) + nDomProc(iProc-1)
    end do

    ! Determine the number of grid clusters
    call determineClusters()

    ! Detertmine the viscous directions in the CGNSBlocks
    call determineViscousDirs()

    ! Determine if we have overset mesh present:
    local = .False.
    do nn=1,nDom
       call setPointers(nn, 1_intType, 1_intType)

       do mm=1, nBocos
          if (BCType(mm) == OversetOuterBound) then
             local = .True.
          end if
       end do
    end do

    call mpi_allreduce(local, oversetPresent, 1, MPI_LOGICAL, MPI_LOR, ADflow_comm_world, ierr)

    ! Loop over the number of levels and perform a lot of tasks.
    ! See the corresponding subroutine header, although the
    ! names are pretty self-explaining


    do level=1,nLevels
       call xhalo(level)
       call allocateMetric(level)
       call metric(level)
       call setPorosities(level)
       call setFamilyInfoFaces(level)
       call faceRotationMatrices(level, .true.)
       call checkSymmetry(level)
       call viscSubfaceInfo(level)
       call determineNcellGlobal(level)
       call setGlobalCellsAndNodes(level)
       call setReferenceVolume(level)
    end do

    ! BC Data must be alloaced (for surface iblank) before we can do
    ! the overset computation.
    call allocMemBCData

    call setSurfaceFamilyInfo

    ! Surface info needs to be computed before the wall distance can
    ! be done and overset connectivity computed
    do level=1,nLevels
       call computeWallDistance(level, .True.)
    end do
    call preprocessingADjoint

  end subroutine preprocessing

  subroutine preprocessingoverset(flag, n, closedFamList, nFam)

    use constants
    use block, only : flowDoms
    use oversetAPI, only : oversetComm, setExplicitHoleCut

    implicit none

    ! Input/Output
    integer(kind=intType), dimension(n) :: flag
    integer(kind=intType), dimension(nFam) :: closedFamList
    integer(kind=intType), intent(in) :: n, nFam

    ! Working
    integer(kind=intType) :: level, nLevels

    nLevels = ubound(flowDoms,2)

    do level=1,nLevels
       if (level == 1) then
          call setExplicitHoleCut(flag)
          call oversetComm(level, .true., .false., closedFamList)
       else
          call oversetComm(level, .True., .True., closedFamList)
       end if
    end do
  end subroutine preprocessingoverset

  subroutine cellRangeSubface
    !
    !       cellRangeSubface determines the cell range for every subface
    !       of every block all grid levels. This subrange can include one
    !       cell of overlap if the boundary coincides with the block
    !       boundary.
    !
    use constants
    use block
    use utils, only : terminate
    implicit none
    !
    !      Local variables.
    !
    integer :: ierr

    integer(kind=intType) :: nLevels, level
    integer(kind=intType) :: nn, mm, il, jl, kl, ie, je, ke
    integer(kind=intType) :: iBeg, jBeg, kBeg, iEnd, jEnd, kEnd

    ! Determine the number of grid levels.

    nLevels = ubound(flowDoms,2)

    ! Loop over the number of grid levels.

    levelLoop: do level=1,nLevels

       ! Loop over the blocks.

       domains: do nn=1,nDom

          ! Allocate the memory for the variables defining the cell
          ! range of the subfaces. Only allocated for the 1st spectral
          ! solution, because this info is identical for all of them.

          mm = flowDoms(nn,level,1)%nSubface
          allocate(flowDoms(nn,level,1)%icBeg(mm), &
               flowDoms(nn,level,1)%jcBeg(mm), &
               flowDoms(nn,level,1)%kcBeg(mm), &
               flowDoms(nn,level,1)%icEnd(mm), &
               flowDoms(nn,level,1)%jcEnd(mm), &
               flowDoms(nn,level,1)%kcEnd(mm), stat=ierr)
          if(ierr /= 0)                        &
               call terminate("cellRangeSubface", &
               "Memory allocation failure for &
               &cell subranges")

          ! Store the nodal dimensions of the block a bit easier.

          il = flowDoms(nn,level,1)%il
          jl = flowDoms(nn,level,1)%jl
          kl = flowDoms(nn,level,1)%kl

          ie = flowDoms(nn,level,1)%ie
          je = flowDoms(nn,level,1)%je
          ke = flowDoms(nn,level,1)%ke

          ! Loop over the number of subfaces for this block.

          subfaces: do mm=1,flowDoms(nn,level,1)%nSubface

             ! Store the nodal range of the subface a bit easier.
             ! Make sure that iBeg, jBeg and kBeg contain the lowest and
             ! iEnd, jEnd and kEnd the highest node numbers.

             iBeg = min(flowDoms(nn,level,1)%inBeg(mm), &
                  flowDoms(nn,level,1)%inEnd(mm))
             iEnd = max(flowDoms(nn,level,1)%inBeg(mm), &
                  flowDoms(nn,level,1)%inEnd(mm))

             jBeg = min(flowDoms(nn,level,1)%jnBeg(mm), &
                  flowDoms(nn,level,1)%jnEnd(mm))
             jEnd = max(flowDoms(nn,level,1)%jnBeg(mm), &
                  flowDoms(nn,level,1)%jnEnd(mm))

             kBeg = min(flowDoms(nn,level,1)%knBeg(mm), &
                  flowDoms(nn,level,1)%knEnd(mm))
             kEnd = max(flowDoms(nn,level,1)%knBeg(mm), &
                  flowDoms(nn,level,1)%knEnd(mm))

             ! Determine the block face on which the subface is located
             ! and set the range accordingly.

             select case (flowDoms(nn,level,1)%BCFaceID(mm))

             case (iMin)
                flowDoms(nn,level,1)%icBeg(mm) = 1
                flowDoms(nn,level,1)%icEnd(mm) = 1

                flowDoms(nn,level,1)%jcBeg(mm) = jBeg +1
                if(jBeg == 1) flowDoms(nn,level,1)%jcBeg(mm) = 1

                flowDoms(nn,level,1)%jcEnd(mm) = jEnd
                if(jEnd == jl) flowDoms(nn,level,1)%jcEnd(mm) = je

                flowDoms(nn,level,1)%kcBeg(mm) = kBeg +1
                if(kBeg == 1) flowDoms(nn,level,1)%kcBeg(mm) = 1

                flowDoms(nn,level,1)%kcEnd(mm) = kEnd
                if(kEnd == kl) flowDoms(nn,level,1)%kcEnd(mm) = ke

                !=========================================================

             case (iMax)
                flowDoms(nn,level,1)%icBeg(mm) = ie
                flowDoms(nn,level,1)%icEnd(mm) = ie

                flowDoms(nn,level,1)%jcBeg(mm) = jBeg +1
                if(jBeg == 1) flowDoms(nn,level,1)%jcBeg(mm) = 1

                flowDoms(nn,level,1)%jcEnd(mm) = jEnd
                if(jEnd == jl) flowDoms(nn,level,1)%jcEnd(mm) = je

                flowDoms(nn,level,1)%kcBeg(mm) = kBeg +1
                if(kBeg == 1) flowDoms(nn,level,1)%kcBeg(mm) = 1

                flowDoms(nn,level,1)%kcEnd(mm) = kEnd
                if(kEnd == kl) flowDoms(nn,level,1)%kcEnd(mm) = ke

                !=========================================================

             case (jMin)
                flowDoms(nn,level,1)%icBeg(mm) = iBeg +1
                if(iBeg == 1) flowDoms(nn,level,1)%icBeg(mm) = 1

                flowDoms(nn,level,1)%icEnd(mm) = iEnd
                if(iEnd == il) flowDoms(nn,level,1)%icEnd(mm) = ie

                flowDoms(nn,level,1)%jcBeg(mm) = 1
                flowDoms(nn,level,1)%jcEnd(mm) = 1

                flowDoms(nn,level,1)%kcBeg(mm) = kBeg +1
                if(kBeg == 1) flowDoms(nn,level,1)%kcBeg(mm) = 1

                flowDoms(nn,level,1)%kcEnd(mm) = kEnd
                if(kEnd == kl) flowDoms(nn,level,1)%kcEnd(mm) = ke

                !=========================================================

             case (jMax)
                flowDoms(nn,level,1)%icBeg(mm) = iBeg +1
                if(iBeg == 1) flowDoms(nn,level,1)%icBeg(mm) = 1

                flowDoms(nn,level,1)%icEnd(mm) = iEnd
                if(iEnd == il) flowDoms(nn,level,1)%icEnd(mm) = ie

                flowDoms(nn,level,1)%jcBeg(mm) = je
                flowDoms(nn,level,1)%jcEnd(mm) = je

                flowDoms(nn,level,1)%kcBeg(mm) = kBeg +1
                if(kBeg == 1) flowDoms(nn,level,1)%kcBeg(mm) = 1

                flowDoms(nn,level,1)%kcEnd(mm) = kEnd
                if(kEnd == kl) flowDoms(nn,level,1)%kcEnd(mm) = ke

                !=========================================================

             case (kMin)
                flowDoms(nn,level,1)%icBeg(mm) = iBeg +1
                if(iBeg == 1) flowDoms(nn,level,1)%icBeg(mm) = 1

                flowDoms(nn,level,1)%icEnd(mm) = iEnd
                if(iEnd == il) flowDoms(nn,level,1)%icEnd(mm) = ie

                flowDoms(nn,level,1)%jcBeg(mm) = jBeg +1
                if(jBeg == 1) flowDoms(nn,level,1)%jcBeg(mm) = 1

                flowDoms(nn,level,1)%jcEnd(mm) = jEnd
                if(jEnd == jl) flowDoms(nn,level,1)%jcEnd(mm) = je

                flowDoms(nn,level,1)%kcBeg(mm) = 1
                flowDoms(nn,level,1)%kcEnd(mm) = 1

                !=========================================================

             case (kMax)
                flowDoms(nn,level,1)%icBeg(mm) = iBeg +1
                if(iBeg == 1) flowDoms(nn,level,1)%icBeg(mm) = 1

                flowDoms(nn,level,1)%icEnd(mm) = iEnd
                if(iEnd == il) flowDoms(nn,level,1)%icEnd(mm) = ie

                flowDoms(nn,level,1)%jcBeg(mm) = jBeg +1
                if(jBeg == 1) flowDoms(nn,level,1)%jcBeg(mm) = 1

                flowDoms(nn,level,1)%jcEnd(mm) = jEnd
                if(jEnd == jl) flowDoms(nn,level,1)%jcEnd(mm) = je

                flowDoms(nn,level,1)%kcBeg(mm) = ke
                flowDoms(nn,level,1)%kcEnd(mm) = ke

             end select

          enddo subfaces
       enddo domains
    enddo levelLoop

  end subroutine cellRangeSubface


  subroutine determineNcellGlobal(level)
    !       determineNcellGlobal determines the global number of cells
    !       the given grid level. This info is needed to compute the L2
    !       norm of the residuals in the flow solver.
    !       Only the 1st spectral solution needs to be considered, because
    !       this info is identical for all of them.
    !
    use constants
    use block
    use communication
    implicit none
    !
    !      Subroutine arguments.
    !
    integer(kind=intType), intent(in) :: level
    !
    !      Local variables.
    !
    integer :: ierr
    integer(kind=intType) :: nn, nCellLocal

    character(len=12) :: int1String, int2String
    !
    ! Determine the local number of cells by looping over the blocks.

    nCellLocal = 0
    do nn=1,nDom
       nCellLocal = nCellLocal + flowDoms(nn,level,1)%nx &
            *              flowDoms(nn,level,1)%ny &
            *              flowDoms(nn,level,1)%nz
    enddo

    ! And determine the global sum.

    call mpi_allreduce(nCellLocal, nCellGlobal(level), 1, &
         adflow_integer, mpi_sum, ADflow_comm_world, ierr)

    ! Write the total number of cells to stdout; only done by
    ! processor 0 to avoid a messy output.

    if(myID == 0) then

       write(int1String,"(i12)") level
       write(int2String,"(i12)") nCellGlobal(level)
       int1String = adjustl(int1String)
       int2String = adjustl(int2String)

       print "(a)", "#"
       print 101, trim(int1String), trim(int2String)
       print "(a)", "#"
101    format("# Grid level: ", a,", Total number of cells: ", a)

    endif

  end subroutine determineNcellGlobal

  subroutine setPorosities(level)
    !
    !       setPorosities sets the porosities for the faces to a certain
    !       flag. Default is normalFlux. The two other possibilities are
    !       boundFlux, used for solid wall boundaries, and noFlux for a
    !       conservative treatment of non matching block boundaries. In
    !       the latter case the flux is constructed differently and the
    !       flux computation in the block must be neglected.
    !       Note that only the 1st spectral solution is treated, because
    !       this informations is the same for all of them.
    !
    use blockPointers
    use constants
    use inputDiscretization
    use utils, only : terminate, setPointers
    implicit none
    !
    !      Subroutine arguments.
    !
    integer(kind=intType), intent(in) :: level
    !
    !      Local variables.
    !
    integer :: ierr

    integer(kind=intType) :: nn, mm, i, j, k

    integer(kind=intType), dimension(2) :: ri, rj, rk

    integer(kind=porType) :: por

    ! Loop over the number of domains.

    domains: do nn=1,nDom

       ! Store the number of nodes in this block a bit easier.

       il = flowDoms(nn,level,1)%il
       jl = flowDoms(nn,level,1)%jl
       kl = flowDoms(nn,level,1)%kl

       ! Allocate the memory for the porosities.

       allocate(flowDoms(nn,level,1)%porI(1:il,2:jl,2:kl), &
            flowDoms(nn,level,1)%porJ(2:il,1:jl,2:kl), &
            flowDoms(nn,level,1)%porK(2:il,2:jl,1:kl), stat=ierr)
       if(ierr /= 0)                     &
            call terminate("setPorosities", &
            "Memory allocation failure for porosities")

       ! Set the pointers for this block to make the source
       ! more readable.

       call setPointers(nn, level, 1_intType)

       ! Initialize the porosities to normalFlux.

       porI = normalFlux
       porJ = normalFlux
       porK = normalFlux

       ! Loop over the subfaces to alter the porosities.

       subface: do mm=1,nsubface

          ! Set the porosity for this subface or continue with the
          ! next if the porosity should not be changed.

          if(BCType(mm) == NSWallAdiabatic  .or. &
               BCType(mm) == NSWallIsothermal .or. &
               BCType(mm) == EulerWall        .or. &
               BCType(mm) == Extrap)    then
             por = boundFlux
          else if(BCType(mm)        == B2BMismatch .and. &
               nonMatchTreatment == Conservative) then
             por = noFlux
          else
             cycle
          endif

          ! Set the range for the faces on this subface.

          ri(1) = min(inBeg(mm), inEnd(mm)) +1
          ri(2) = max(inBeg(mm), inEnd(mm))

          rj(1) = min(jnBeg(mm), jnEnd(mm)) +1
          rj(2) = max(jnBeg(mm), jnEnd(mm))

          rk(1) = min(knBeg(mm), knEnd(mm)) +1
          rk(2) = max(knBeg(mm), knEnd(mm))

          ! Determine the block face this subface is located on and
          ! set the corresponding porosities correctly.

          select case( BCFaceID(mm) )

          case (iMin)
             do k=rk(1),rk(2)
                do j=rj(1),rj(2)
                   porI(1,j,k) = por
                enddo
             enddo

             !===========================================================

          case (iMax)
             do k=rk(1),rk(2)
                do j=rj(1),rj(2)
                   porI(il,j,k) = por
                enddo
             enddo

             !===========================================================

          case (jMin)
             do k=rk(1),rk(2)
                do i=ri(1),ri(2)
                   porJ(i,1,k) = por
                enddo
             enddo

             !===========================================================

          case (jMax)
             do k=rk(1),rk(2)
                do i=ri(1),ri(2)
                   porJ(i,jl,k) = por
                enddo
             enddo

             !===========================================================

          case (kMin)
             do j=rj(1),rj(2)
                do i=ri(1),ri(2)
                   porK(i,j,1) = por
                enddo
             enddo

             !===========================================================

          case (kMax)
             do j=rj(1),rj(2)
                do i=ri(1),ri(2)
                   porK(i,j,kl) = por
                enddo
             enddo

          end select

       enddo subface

    enddo domains

  end subroutine setPorosities

  subroutine exchangeGlobalCells(level, sps, commPattern, internal)
    !
    !       ExchangeIblank exchanges the 1 to 1 internal halo's for the
    !       given level and sps instance.
    !
    use constants
    use block
    use communication
    use utils, only : terminate
    implicit none
    !
    !      Subroutine arguments.
    !
    integer(kind=intType), intent(in) :: level, sps

    type(commType),          dimension(*), intent(in) :: commPattern
    type(internalCommType), dimension(*), intent(in) :: internal
    !
    !      Local variables.
    !
    integer :: size, procId, ierr, index
    integer, dimension(mpi_status_size) :: mpiStatus

    integer(kind=intType) :: i, j, ii, jj
    integer(kind=intType) :: d1, i1, j1, k1, d2, i2, j2, k2

    integer(kind=intType), dimension(:), allocatable :: sendBufInt
    integer(kind=intType), dimension(:), allocatable :: recvBufInt

    ! Allocate the memory for the sending and receiving buffers.

    ii = commPattern(level)%nProcSend
    ii = commPattern(level)%nsendCum(ii)
    jj = commPattern(level)%nProcRecv
    jj = commPattern(level)%nrecvCum(jj)

    allocate(sendBufInt(ii), recvBufInt(jj), stat=ierr)
    if(ierr /= 0)                       &
         call terminate("exchangeIblank", &
         "Memory allocation failure for buffers")

    ! Send the variables. The data is first copied into
    ! the send buffer after which the buffer is sent asap.

    ii = 1
    sends: do i=1,commPattern(level)%nProcSend

       ! Store the processor id and the size of the message
       ! a bit easier.

       procID = commPattern(level)%sendProc(i)
       size    = commPattern(level)%nsend(i)

       ! Copy the data in the correct part of the send buffer.

       jj = ii
       do j=1,commPattern(level)%nsend(i)

          ! Store the block id and the indices of the donor
          ! a bit easier.

          d1 = commPattern(level)%sendList(i)%block(j)
          i1 = commPattern(level)%sendList(i)%indices(j,1)
          j1 = commPattern(level)%sendList(i)%indices(j,2)
          k1 = commPattern(level)%sendList(i)%indices(j,3)

          ! Copy globalCell values to buffer.

          sendBufInt(jj) = flowDoms(d1,level,sps)%globalCell(i1,j1,k1)
          jj = jj + 1

       enddo

       ! Send the data.

       call mpi_isend(sendBufInt(ii), size, adflow_integer, procId, &
            procId, ADflow_comm_world, sendRequests(i),   &
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
       size    = commPattern(level)%nrecv(i)

       ! Post the receive.

       call mpi_irecv(recvBufInt(ii), size, adflow_integer, procId, &
            myId, ADflow_comm_world, recvRequests(i), ierr)

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

       ! Copy the globalCell value

       flowDoms(d2,level,sps)%globalCell(i2,j2,k2) = flowDoms(d1,level,sps)%globalCell(i1,j1,k1)

    enddo localCopy

    ! Complete the nonblocking receives in an arbitrary sequence and
    ! copy the variables from the buffer into the halo's.

    size = commPattern(level)%nProcRecv
    completeRecvs: do i=1,commPattern(level)%nProcRecv

       ! Complete any of the requests.

       call mpi_waitany(size, recvRequests, index, mpiStatus, ierr)

       ! Copy the data just arrived in the halo's.

       ii = index
       jj = commPattern(level)%nrecvCum(ii-1)
       do j=1,commPattern(level)%nrecv(ii)

          ! Store the block and the indices of the halo a bit easier.

          d2 = commPattern(level)%recvList(ii)%block(j)
          i2 = commPattern(level)%recvList(ii)%indices(j,1)
          j2 = commPattern(level)%recvList(ii)%indices(j,2)
          k2 = commPattern(level)%recvList(ii)%indices(j,3)

          ! Copy the globalCell value

          jj = jj + 1
          flowDoms(d2,level,sps)%globalCell(i2,j2,k2) = recvBufInt(jj)

       enddo

    enddo completeRecvs

    ! Complete the nonblocking sends.

    size = commPattern(level)%nProcSend
    do i=1,commPattern(level)%nProcSend
       call mpi_waitany(size, sendRequests, index, mpiStatus, ierr)
    enddo

    ! Deallocate the memory for the sending and receiving buffers.

    deallocate(sendBufInt, recvBufInt, stat=ierr)
    if(ierr /= 0)                       &
         call terminate("exchangeGlobalCell", &
         "Deallocation failure for buffers")

  end subroutine exchangeGlobalCells

  subroutine checkSymmetry(level)
    !
    !       checkSymmetry checks whether or not the symmetry planes are
    !       really planar (within a certain tolerance). If this is not the
    !       case for the finest level, a warning is printed. In all cases
    !       the unit normals are replaced by the face averaged unit
    !       normal.
    !
    use constants
    use blockPointers
    use cgnsGrid
    use inputTimeSpectral
    use utils, only : setPointers
    implicit none
    !
    !      Subroutine arguments.
    !
    integer(kind=intType), intent(in) :: level
    !
    !      Local parameter, tolerance for planar, 0.1 degrees.
    !
    real(kind=realType), parameter :: tolDotmin = 0.9999985_realType
    !
    !      Local variables.
    !
    integer(kind=intType) :: nn, mm, sps, i, j

    real(kind=realType) :: fact, dotMin, dot, mult

    real(kind=realType), dimension(3) :: faceNorm
    real(kind=realType), dimension(:,:,:), pointer :: ss

    ! Loop over the number of spectral solutions and local domains.

    spectral: do sps=1,nTimeIntervalsSpectral
       domains: do nn=1,nDom

          ! Set the pointers for this block.

          call setPointers(nn, level, sps)

          ! Loop over the number of boundary subfaces for this block.

          bocos: do mm=1,nBocos

             ! Check for symmetry boundary condition.

             symmetry: if(BCType(mm) == symm) then

                ! Determine the block face on which this subface is
                ! located and set some variables accordingly.

                select case (BCFaceID(mm))

                case (iMin)
                   mult = -one; ss => si(1,:,:,:)
                case (iMax)
                   mult = one; ss => si(il,:,:,:)
                case (jMin)
                   mult = -one; ss => sj(:,1,:,:)
                case (jMax)
                   mult = one; ss => sj(:,jl,:,:)
                case (kMin)
                   mult = -one; ss => sk(:,:,1,:)
                case (kMax)
                   mult = one; ss => sk(:,:,kl,:)

                end select

                ! Loop over the range of the subface compute the face
                ! normal. The halo cells should not be taken into account,
                ! which explains why the nodal range of BCData is used.
                ! As the starting index of the cell range is shifted 1,
                ! (inBeg+1) and (jnBeg+1) are the starting indices for the
                ! owned cell range.

                faceNorm = zero

                do j=(BCData(mm)%jnBeg+1), BCData(mm)%jnEnd
                   do i=(BCData(mm)%inBeg+1), BCData(mm)%inEnd
                      faceNorm(1) = faceNorm(1) + ss(i,j,1)
                      faceNorm(2) = faceNorm(2) + ss(i,j,2)
                      faceNorm(3) = faceNorm(3) + ss(i,j,3)
                   enddo
                enddo

                ! Create the unit normal for faceNorm. Make sure it
                ! is outward pointing by multiplying it by mult;
                ! mult is either 1.0 or -1.0.

                fact = sqrt(faceNorm(1)*faceNorm(1) &
                     +      faceNorm(2)*faceNorm(2) &
                     +      faceNorm(3)*faceNorm(3))
                if(fact > zero) fact = mult/fact

                faceNorm(1) = faceNorm(1)*fact
                faceNorm(2) = faceNorm(2)*fact
                faceNorm(3) = faceNorm(3)*fact

                ! Check if the symmetry plane is really planar. This is
                ! only done on the finest mesh and for the 1st spectral
                ! solution, because it is only to inform the user.
                ! Afterwards the normals will be reset to the unit
                ! normal of the face anyway.

                fineLevelTest: if(level == 1 .and. sps == 1) then

                   ! Initialize dotMin such that it will always
                   ! be overwritten.

                   dotMin = one

                   ! Loop over the physical faces of the symmetry plane,
                   ! i.e. no halo's.

                   do j=(BCData(mm)%jnBeg+1), BCData(mm)%jnEnd
                      do i=(BCData(mm)%inBeg+1), BCData(mm)%inEnd

                         ! Compute the dot product between the normal of
                         ! this face and the averaged normal of the plane.

                         dot = BCData(mm)%norm(i,j,1)*faceNorm(1) &
                              + BCData(mm)%norm(i,j,2)*faceNorm(2) &
                              + BCData(mm)%norm(i,j,3)*faceNorm(3)

                         ! And determine the minimum of dot and dotMin

                         dotMin = min(dot,dotMin)
                      enddo
                   enddo

                   ! Test if the minimum dot product is smaller than the
                   ! tolerance. If so, the plane is considered as not
                   ! planar.

                   if(dotMin < tolDotmin) then

                      ! Determine the corresponding angle in degrees of
                      ! dotmin.

                      fact = acos(dotMin)*180.0_realType/pi

                      ! Store the corresponding cgns block id and the
                      ! subface in this block a bit easier.

                      i = nbkGlobal
                      j = cgnsSubface(mm)

                      ! Print a warning.

                      print "(a)", "#"
                      print "(a)", "#                      Warning"
                      print 100,                              &
                           trim(cgnsDoms(i)%bocoInfo(j)%bocoName), &
                           trim(cgnsDoms(i)%zonename)
                      print 110, fact
                      print "(a)", "#"
100                   format("# Symmetry boundary face",1X,A,1X,"of zone", &
                           1x,a,1x, "is not planar.")
110                   format("# Maximum deviation from the mean normal: ", &
                           e12.5, " degrees")

                   endif

                endif fineLevelTest

                !removed as this would cause issues with the ADjoint
!!$               ! Set the unit normals to the unit normal of the entire
!!$               ! plane. All the cells, also possible halo's, are treated.
!!$
!!$               do j=BCData(mm)%jcBeg, BCData(mm)%jcEnd
!!$                 do i=BCData(mm)%icBeg, BCData(mm)%icEnd
!!$
!!$                   BCData(mm)%norm(i,j,1) = faceNorm(1)
!!$                   BCData(mm)%norm(i,j,2) = faceNorm(2)
!!$                   BCData(mm)%norm(i,j,3) = faceNorm(3)
!!$
!!$                 enddo
!!$               enddo

             endif symmetry
          enddo bocos
       enddo domains
    enddo spectral

  end subroutine checkSymmetry
  subroutine xhalo(level)
    !
    !       xhalo determines the coordinates of the nodal halo's.
    !       First it sets all halo coordinates by simple extrapolation,
    !       then the symmetry planes are treated (also the unit normal of
    !       symmetry planes are determined) and finally an exchange is
    !       made for the internal halo's.
    !
    use constants
    use blockPointers
    use communication
    use inputTimeSpectral
    use utils, only : setPointers
    use haloExchange, only : exchangeCoor
    implicit none
    !
    !      Subroutine arguments.
    !
    integer(kind=intType) :: level
    !
    !      Local variables.
    !
    integer(kind=intType) :: nn, mm, sps, i, j, k
    integer(kind=intType) :: iBeg, iEnd, jBeg, jEnd, iiMax, jjMax

    real(kind=realType), dimension(:,:,:), pointer :: x0, x1, x2

    real(kind=realType) :: length, dot

    real(kind=realType), dimension(3) :: v1, v2, norm, tmp, tmp2
    real(kind=realType), parameter :: tolDotmin = 0.99_realType


    ! Loop over the number of spectral solutions and the local
    ! number of blocks.

    spectralLoop: do sps=1,nTimeIntervalsSpectral
       domains: do nn=1,nDom

          ! Set the pointers to this block.

          call setPointers(nn, level, sps)

          !
          !           Extrapolation of the coordinates. First extrapolation in
          !           i-direction, without halo's, followed by extrapolation in
          !           j-direction, with i-halo's and finally extrapolation in
          !           k-direction, with both i- and j-halo's. In this way also
          !           the indirect halo's get a value, albeit a bit arbitrary.
          !
          ! Extrapolation in i-direction.

          do k=1,kl
             do j=1,jl
                x(0,j,k,1) = two*x(1,j,k,1) - x(2,j,k,1)
                x(0,j,k,2) = two*x(1,j,k,2) - x(2,j,k,2)
                x(0,j,k,3) = two*x(1,j,k,3) - x(2,j,k,3)

                x(ie,j,k,1) = two*x(il,j,k,1) - x(nx,j,k,1)
                x(ie,j,k,2) = two*x(il,j,k,2) - x(nx,j,k,2)
                x(ie,j,k,3) = two*x(il,j,k,3) - x(nx,j,k,3)
             enddo
          enddo

          ! Extrapolation in j-direction.

          do k=1,kl
             do i=0,ie
                x(i,0,k,1) = two*x(i,1,k,1) - x(i,2,k,1)
                x(i,0,k,2) = two*x(i,1,k,2) - x(i,2,k,2)
                x(i,0,k,3) = two*x(i,1,k,3) - x(i,2,k,3)

                x(i,je,k,1) = two*x(i,jl,k,1) - x(i,ny,k,1)
                x(i,je,k,2) = two*x(i,jl,k,2) - x(i,ny,k,2)
                x(i,je,k,3) = two*x(i,jl,k,3) - x(i,ny,k,3)
             enddo
          enddo

          ! Extrapolation in k-direction.

          do j=0,je
             do i=0,ie
                x(i,j,0,1) = two*x(i,j,1,1) - x(i,j,2,1)
                x(i,j,0,2) = two*x(i,j,1,2) - x(i,j,2,2)
                x(i,j,0,3) = two*x(i,j,1,3) - x(i,j,2,3)

                x(i,j,ke,1) = two*x(i,j,kl,1) - x(i,j,nz,1)
                x(i,j,ke,2) = two*x(i,j,kl,2) - x(i,j,nz,2)
                x(i,j,ke,3) = two*x(i,j,kl,3) - x(i,j,nz,3)
             enddo
          enddo
          !
          !           Mirror the halo coordinates adjacent to the symmetry
          !           planes
          !
          ! Loop over boundary subfaces.

          loopBocos: do mm=1,nBocos
             ! The actual correction of the coordinates only takes
             ! place for symmetry planes.

             testSymmetry: if(BCType(mm) == Symm) then
                ! Set some variables, depending on the block face on
                ! which the subface is located.

                select case (BCFaceID(mm))
                case (iMin)
                   iBeg = jnBeg(mm); iEnd = jnEnd(mm); iiMax = jl
                   jBeg = knBeg(mm); jEnd = knEnd(mm); jjMax = kl
                   x0 => x(0,:,:,:); x1 => x(1,:,:,:); x2 => x(2,:,:,:)

                case (iMax)
                   iBeg = jnBeg(mm); iEnd = jnEnd(mm); iiMax = jl
                   jBeg = knBeg(mm); jEnd = knEnd(mm); jjMax = kl
                   x0 => x(ie,:,:,:); x1 => x(il,:,:,:); x2 => x(nx,:,:,:)

                case (jMin)
                   iBeg = inBeg(mm); iEnd = inEnd(mm); iiMax = il
                   jBeg = knBeg(mm); jEnd = knEnd(mm); jjMax = kl
                   x0 => x(:,0,:,:); x1 => x(:,1,:,:); x2 => x(:,2,:,:)

                case (jMax)
                   iBeg = inBeg(mm); iEnd = inEnd(mm); iiMax = il
                   jBeg = knBeg(mm); jEnd = knEnd(mm); jjMax = kl
                   x0 => x(:,je,:,:); x1 => x(:,jl,:,:); x2 => x(:,ny,:,:)

                case (kMin)
                   iBeg = inBeg(mm); iEnd = inEnd(mm); iiMax = il
                   jBeg = jnBeg(mm); jEnd = jnEnd(mm); jjMax = jl
                   x0 => x(:,:,0,:); x1 => x(:,:,1,:); x2 => x(:,:,2,:)

                case (kMax)
                   iBeg = inBeg(mm); iEnd = inEnd(mm); iiMax = il
                   jBeg = jnBeg(mm); jEnd = jnEnd(mm); jjMax = jl
                   x0 => x(:,:,ke,:); x1 => x(:,:,kl,:); x2 => x(:,:,nz,:)
                end select


                ! Determine the vector from the lower left corner to
                ! the upper right corner. Due to the usage of pointers
                ! an offset of +1 must be used, because the original
                ! array x start at 0.

                v1(1) = x1(iimax+1,jjmax+1,1) - x1(1+1,1+1,1)
                v1(2) = x1(iimax+1,jjmax+1,2) - x1(1+1,1+1,2)
                v1(3) = x1(iimax+1,jjmax+1,3) - x1(1+1,1+1,3)

                ! And the vector from the upper left corner to the
                ! lower right corner.

                v2(1) = x1(iimax+1,1+1,1) - x1(1+1,jjmax+1,1)
                v2(2) = x1(iimax+1,1+1,2) - x1(1+1,jjmax+1,2)
                v2(3) = x1(iimax+1,1+1,3) - x1(1+1,jjmax+1,3)

                ! Determine the normal of the face by taking the cross
                ! product of v1 and v2 and add it to norm.

                norm(1) = v1(2)*v2(3) - v1(3)*v2(2)
                norm(2) = v1(3)*v2(1) - v1(1)*v2(3)
                norm(3) = v1(1)*v2(2) - v1(2)*v2(1)

                ! Check if BCData is allocated yet:
                if (.not. bcData(mm)%symNormSet) then
                   length = sqrt(norm(1)**2 + norm(2)**2 + norm(3)**2)
                   if (length == 0) then
                      length = eps
                   end if
                   bcData(mm)%symNorm(1) = norm(1)/length
                   bcData(mm)%symNorm(2) = norm(2)/length
                   bcData(mm)%symNorm(3) = norm(3)/length
                   bcData(mm)%symNormSet = .True.
                else

                   ! Check that the orientation of norm() is not
                   ! different from the stored one:
                   length = sqrt(norm(1)**2 + norm(2)**2 + norm(3)**2)
                   if (length > eps) then
                      tmp = norm / length
                      tmp2 = bcData(mm)%symNorm
                      dot = dot_product(tmp, tmp2)
                      if (abs(dot) < tolDotmin) then
                         print *, 'Symmetry Plane normal has changed from initial configuration. Resetting.'
                         print *, 'This may cause a slightly inaccurate gradient!'
                         bcData(mm)%symNorm(1) = norm(1)
                         bcData(mm)%symNorm(2) = norm(2)
                         bcData(mm)%symNorm(3) = norm(3)
                      end if
                   end if

                   ! Copy out the saved symNorm
                   norm(1) = bcData(mm)%symNorm(1)
                   norm(2) = bcData(mm)%symNorm(2)
                   norm(3) = bcData(mm)%symNorm(3)
                end if

                ! Compute the length of the normal and test if this is
                ! larger than eps. If this is the case this means that
                ! it is a nonsingular subface and the coordinates are
                ! corrected.

                length = sqrt(norm(1)**2 + norm(2)**2 + norm(3)**2)

                testSingular: if(length > eps) then

                   ! Compute the unit normal of the subface.

                   norm(1) = norm(1)/length
                   norm(2) = norm(2)/length
                   norm(3) = norm(3)/length

                   ! Add an overlap to the symmetry subface if the
                   ! boundaries coincide with the block boundaries.
                   ! This way the indirect halo's are treated properly.

                   if(iBeg == 1)     iBeg = 0
                   if(iEnd == iiMax) iEnd = iiMax + 1

                   if(jBeg == 1)     jBeg = 0
                   if(jEnd == jjMax) jEnd = jjMax + 1

                   ! Loop over the nodes of the subface and set the
                   ! corresponding halo coordinates.

                   do j=jBeg,jEnd
                      do i=iBeg,iEnd

                         ! Determine the vector from the internal node to the
                         ! node on the face. Again an offset of +1 must be
                         ! used, due to the usage of pointers.

                         v1(1) = x1(i+1,j+1,1) - x2(i+1,j+1,1)
                         v1(2) = x1(i+1,j+1,2) - x2(i+1,j+1,2)
                         v1(3) = x1(i+1,j+1,3) - x2(i+1,j+1,3)

                         ! Determine two times the normal component of this
                         ! vector; this vector must be added to the
                         ! coordinates of the internal node to obtain the
                         ! halo coordinates. Again the offset of +1.

                         dot = two*(v1(1)*norm(1) + v1(2)*norm(2) &
                              +      v1(3)*norm(3))

                         x0(i+1,j+1,1) = x2(i+1,j+1,1) + dot*norm(1)
                         x0(i+1,j+1,2) = x2(i+1,j+1,2) + dot*norm(2)
                         x0(i+1,j+1,3) = x2(i+1,j+1,3) + dot*norm(3)

                      enddo
                   enddo
                endif testSingular
             endif testSymmetry
          enddo loopBocos
       enddo domains
    enddo spectralLoop

    !
    !       Exchange the coordinates for the internal halo's.
    !
    call exchangeCoor(level)

  end subroutine xhalo

  subroutine setSurfaceFamilyInfo

    use constants
    use su_cgns
    use blockPointers, onlY : nDom, flowDoms, nBocos, cgnsSubFace, BCType, BCData
    use cgnsGrid, onlY : cgnsDoms
    use communication, only : myid, adflow_comm_world, nProc
    use inputTimeSpectral, only : nTimeIntervalsSpectral
    use surfaceFamilies, only : BCFamExchange, famNames, fullFamList, &
         zeroCellVal, zeroNodeVal, oneCellVal, BCFamgroups
    use utils, only : setPointers, EChk, pointReduce, terminate, convertToLowerCase
    use sorting, only : qsortStrings, bsearchStrings, famInList
    use surfaceUtils, only : getSurfaceSize
    implicit none

    integer :: ierr
    integer(kind=intType) :: nLevels, level, nn, mm, nsMin, nsMax, i, j, k, nFam, famID, cgb, iFam
    integer(kind=intType) :: sps, isizemax, jsizemax, totalFamilies, totalWallFamilies
    integer(kind=intType) :: iBeg, iEnd, jBeg, jEnd, ii, iBCGroup, totalBCFamilies
    character(maxCGNSNameLen), dimension(25) :: defaultFamName
    character(maxCGNSNameLen) :: curStr, family
    character(maxCGNSNameLen), dimension(:), allocatable :: uniqueFamListNames
    integer(kind=intType), dimension(:), allocatable :: localFlag, famIsPartOfBCGroup
    integer(kind=intType), dimension(:), allocatable :: localIndices, nodeSizes, nodeDisps
    integer(kind=intType) :: iProc, nodeSize, cellSize

    ! Process out the family information. The goal here is to
    ! assign a unique integer to each family in each boundary
    ! condition. The CGNS grid has all the information we need.

    ! Firstly make sure that there is actual family specified for
    ! each BC. If there isn't, we will provide one for you.
    defaultFamName(BCAxisymmetricWedge) = 'axi'
    defaultFamName(BCDegenerateLine) = 'degenerate'
    defaultFamName(BCDegeneratePoint) ='degenerate'
    defaultFamName(BCDirichlet) = 'dirichlet'
    defaultFamName(BCExtrapolate) = 'extrap'
    defaultFamName(BCFarfield) = 'far'
    defaultFamName(BCGeneral) = 'general'
    defaultFamName(BCInflow) = 'inflow'
    defaultFamName(BCInflowSubsonic) = 'inflow'
    defaultFamName(BCInflowSupersonic) = 'inflow'
    defaultFamName(BCNeumann) = 'neumann'
    defaultFamName(BCOutflow) = 'outflow'
    defaultFamName(BCOutflowSubsonic) = 'outflow'
    defaultFamName(BCOutflowSupersonic)  ='outflow'
    defaultFamName(BCSymmetryPlane) = 'sym'
    defaultFamName(BCSymmetryPolar) = 'sympolar'
    defaultFamName(BCTunnelInflow) = 'inflow'
    defaultFamName(BCTunnelOutflow) = 'outflow'
    defaultFamName(BCWall) = 'wall'
    defaultFamName(BCWallInviscid) = 'wall'
    defaultFamName(BCWallViscous) = 'wall'
    defaultFamName(BCWallViscousHeatFlux) = 'wall'
    defaultFamName(BCWallViscousIsothermal) = 'wall'
    defaultFamName(UserDefined) = 'userDefined'

101 format("CGNS Block ",I4,", boundary condition ",I4, ", of type ",a, &
         " does not have a family. Based on the boundary condition type," &
         " a name of: '", a, "' will be used.")

    nFam = 0
    do i=1, size(cgnsDoms)
       do j=1, size(cgnsDoms(i)%bocoInfo)
          if (cgnsDoms(i)%bocoInfo(j)%actualFace) then
             if (trim(cgnsDoms(i)%bocoInfo(j)%wallBCName) == "") then
                if (myid == 0) then
                   ! Tell the user we are adding an automatic family name
                   write(*, 101) i, j, trim(BCTypeName(cgnsDoms(i)%bocoInfo(j)%BCTypeCGNS)), &
                        trim(defaultFamName(cgnsDoms(i)%bocoInfo(j)%BCTypeCGNS))
                end if
                cgnsDoms(i)%bocoInfo(j)%wallBCName = trim(defaultFamName(cgnsDoms(i)%bocoInfo(j)%BCTypeCGNS))
             end if
             nFam = nFam + 1
          end if
       end do
    end do

    ! Allocate space for the full family list
    allocate(famNames(nFam))
    nFam = 0
    do i=1, size(cgnsDoms)
       do j=1, size(cgnsDoms(i)%bocoInfo)
          if (cgnsDoms(i)%bocoInfo(j)%actualFace) then
             nFam = nFam + 1
             famNames(nfam) = cgnsDoms(i)%bocoInfo(j)%wallBCName
             call convertToLowerCase(famNames(nFam))
          end if
       end do
    end do

    ! Now sort the family names:
    call qsortStrings(famNames, nFam)

    ! Next we need to generate a unique set of names.
    allocate(uniqueFamListNames(nFam))

    curStr = famNames(1)
    uniqueFamListNames(1) = curStr
    j = 1
    i = 1
    do while(i < nFam)

       i = i + 1
       if (famNames(i) == curStr) then
          ! Same str, do nothing.
       else
          j = j + 1
          curStr = famNames(i)
          uniqueFamListNames(j) = curStr
       end if
    end do


    totalFamilies = j
    ! Now copy the uniqueFamListNames back to "famNames" and allocate
    ! exactly the right size.
    deallocate(famNames)
    allocate(famNames(totalFamilies))
    famNames(1:totalFamilies) = uniqueFamListNames(1:totalFamilies)
    deallocate(uniqueFamListNames)

    ! Now each block boundary condition can uniquely determine it's
    ! famID. We do all BC on all blocks and levels.
    nLevels = ubound(flowDoms,2)
    do nn=1, nDom
       call setPointers(nn, 1_intType, 1_intType)
       do mm=1, nBocos

          cgb = flowDoms(nn, 1, 1)%cgnsBlockID
          family = cgnsDoms(cgb)%bocoInfo(cgnsSubface(mm))%wallBCName
          call convertToLowerCase(family)

          famID = bsearchStrings(family, famNames)
          if (famID == 0) then
             ! Somehow we never found the family...
             call terminate("setSurfaceFamilyInfo", &
                  "An error occuring in assigning families")
          end if

          ! Now set the data on each of the level/sps instances
          do sps=1, nTimeIntervalsSpectral
             do level=1,nlevels

                flowDoms(nn, level, sps)%bcData(mm)%famID = famID
                flowDoms(nn, level, sps)%bcData(mm)%family = family

             end do
          end do
       end do
    end do

    ! Next we need to group the families based on their boundary
    ! condition. The reason for this is that we generate the reduction
    ! scatterd based on groups of BC types. Specifically the following groups:

    ! 1. Walls : EulerWall, NSWallAdiabatic, NSWallIsothermal
    ! 2. Symm : Symm, SymmPolar
    ! 3. Inflow/Outflow : subSonicInflow, subSonicOutflow, supersonicInflow, superSonicOutflow
    ! 4. Farfield : Farfield
    ! 5. Overset : OversetouterBound
    ! 6. Others : All remaining BCs

    ! The final familyExchange structure.
    allocate(BCFamExchange(nFamExchange, nTimeIntervalsSpectral), localFlag(totalFamilies))

    BCGroupLoop: do iBCGroup=1, nfamExchange
       localFlag = 0
       ! Determine which of the unique families match the specific
       ! BCGroup.  This is slightly inefficient but not it isn't
       ! performance critical.
       famLoop: do iFam=1, totalFamilies
          domainLoop: do nn=1,nDom
             call setPointers(nn, 1_intType, 1_intType)
             bocoLoop: do mm=1, nBocos
                matchiFam: if (flowDoms(nn, 1, 1)%bcData(mm)%famID == iFam) then
                   select case(iBCGroup)

                   case (iBCGroupWalls)
                      if (BCType(mm) == EulerWall .or. &
                           BCType(mm) == NSWallAdiabatic .or. &
                           BCType(mm) == NSwallIsoThermal) then
                         localFlag(iFam) = 1
                      end if

                   case (iBCGroupInflow)
                      if (BCType(mm) == SubsonicInflow .or. &
                           BCType(mm) == SupersonicInflow) then
                         localFlag(iFam) = 1
                      end if

                   case (iBCGroupOutflow)
                      if (BCType(mm) == SubsonicOutflow .or. &
                           BCType(mm) == SupersonicOutflow) then
                         localFlag(iFam) = 1
                      end if

                   case (iBCGroupSymm)
                      if (BCType(mm) == Symm .or. BCType(mm) == SymmPolar) then
                         localFlag(iFam) = 1
                      end if

                   case (iBCGroupFarfield)
                      if (BCType(mm) == Farfield) then
                         localFlag(iFam) = 1
                      end if

                   case (iBCGroupOverset)
                      if (BCType(mm) == OversetOuterBound) then
                         localFlag(iFam) = 1
                      end if

                   case (iBCGroupOther)
                      ! All other boundary conditions. Note that some
                      ! of these are not actually implemented
                      if (BCType(mm) == BCNull .or. &
                           BCType(mm) == MassBleedInflow .or. &
                           BCType(mm) == MassbleedOutflow .or. &
                           BCType(mm) == mDot .or. &
                           BCType(mm) == BCThrust .or. &
                           BCType(mm) == Extrap .or. &
                           BCType(mm) == B2BMatch .or. &
                           BCType(mm) == B2BMisMatch .or. &
                           BCType(mm) == SlidingInterface .or. &
                           BCType(mm) == DomainInterfaceAll .or. &
                           BCType(mm) == DomainInterfaceRhoUVW .or. &
                           BCType(mm) == DomainInterfaceP .or. &
                           BCType(mm) == DomainInterfaceRho .or. &
                           BCType(mm) == DomainInterfaceTotal) then
                         localFlag(iFam) = 1
                      end if
                   end select
                end if matchiFam
             end do bocoLoop
          end do domainLoop
       end do famLoop

       ! All Reduce so all procs know the same information.
       allocate(famIsPartOfBCGroup(totalFamilies))
       call mpi_allreduce(localFlag, famIsPartOfBCGroup, totalFamilies, &
            adflow_integer, MPI_SUM, adflow_comm_world, ierr)
       call EChk(ierr, __FILE__, __LINE__)

       ! Count up the number of families this BC has.
       totalBCFamilies = 0
       do i=1, totalFamilies
          if (famIsPartOfBCGroup(i) > 0) then
             totalBCFamilies = totalBCFamilies + 1
          end if
       end do

       ! Allocate the space for the list of fam for each BC and set.
       allocate(BCFamGroups(iBCGroup)%famList(totalBCFamilies))
       k = 0
       do i=1, totalFamilies
          if (famIsPartOfBCGroup(i) > 0) then
             k = k + 1
             BCFamGroups(iBCGroup)%famList(k) = i
          end if
       end do
       deallocate(famIsPartOfBCGroup)
    end do BCGroupLoop

    ! Dump a little information out to the user giving the family and
    ! the BC types. This will probably be useful in general.

    if (myid == 0) then
       write(*, "(a)") '+--------------------------------------------------+'
       write(*, "(a)") '  CGNS Surface Families by Boundary Condition Type'
       write(*, "(a)") '+--------------------------------------------------+'

       do iBCGroup=1,6
          select case(iBCGroup)
          case (iBCGroupWalls)
             write(*,"(a)",advance="no") '| Wall Types           : '
          case (iBCGroupInflow)
             write(*,"(a)",advance="no") '| Inflow Types : '
          case (iBCGroupOutflow)
             write(*,"(a)",advance="no") '| Outflow Types : '
          case (iBCGroupSymm)
             write(*,"(a)",advance="no") '| Symmetry Types       : '
          case (iBCGroupFarfield)
             write(*,"(a)",advance="no") '| Farfield Types       : '
          case (iBCGroupOverset)
             write(*,"(a)",advance="no") '| Oveset Types         : '
          case (iBCGroupOther)
             write(*,"(a)",advance="no") '| Other Types          : '
          end select

          do i=1,size(BCFamGroups(iBCGroup)%famList)
             write(*,"(a,1x)",advance="no") trim(famNames(BCFamGroups(iBCGroup)%famList(i)))
          end do
          print "(1x)"
       end do
       write(*, "(a)") '+--------------------------------------------------+'
    end if

    ! Generate the node scatters for each family. This will also tell
    ! us the surfaceIndex for each BC. This is the index into the
    ! *gloablly reduced vector*. This is what we will need for tecplot
    ! output as well as the zipper mesh computations.

    do iBCGroup=1, nFamExchange
       do sps=1,nTimeIntervalsSpectral
          call createNodeScatterForFamilies(&
               BCFamGroups(iBCGroup)%famList, BCFamExchange(iBCGroup, sps), sps, localIndices)

          ! this won't include the zipper nodes since that isn't done yet.
          call getSurfaceSize(nodeSize, cellSize, BCFamGroups(iBCGroup)%famList, &
               size(BCFamGroups(iBCGroup)%famlist), .False.)
          allocate(nodeSizes(nProc), nodeDisps(0:nProc))
          nodeSizes = 0
          nodeDisps = 0

          call mpi_allgather(nodeSize, 1, adflow_integer, nodeSizes, 1, adflow_integer, &
               adflow_comm_world, ierr)
          call EChk(ierr,__FILE__,__LINE__)
          nodeDisps(0) = 0
          do iProc=1, nProc
             nodeDisps(iProc) = nodeDisps(iProc-1) + nodeSizes(iProc)
          end do


          ii = 0
          do nn=1, nDom
             call setPointers(nn, 1, sps)
             do mm=1, nBocos
                famInclude: if (famInList(BCData(mm)%famId, BCFamGroups(iBCGroup)%famList)) then
                   iBeg = BCData(mm)%inbeg; iEnd = BCData(mm)%inend
                   jBeg = BCData(mm)%jnbeg; jEnd = BCData(mm)%jnend
                   do j=jBeg, jEnd
                      do i=iBeg, iEnd
                         ii = ii + 1
                         BCData(mm)%surfIndex(i,j) = ii + nodeDisps(myid)
                      end do
                   end do
                end if famInclude
             end do
          end do
          deallocate(localIndices, nodeSizes, nodeDisps)
       end do
    end do
    ! Allocate  arrays that have the maximum face size. These may
    ! be slightly larger than necessary, but that's ok. We just need
    ! somethwere to point the pointers.
    isizemax = 0
    jsizemax = 0
    do nn=1,nDom
       isizemax = max(isizemax, flowDoms(nn, 1, 1)%ie)
       isizemax = max(isizemax, flowDoms(nn, 1, 1)%je)

       jsizemax = max(jsizemax, flowDoms(nn, 1, 1)%je)
       jsizemax = max(jsizemax, flowDoms(nn, 1, 1)%ke)
    end do


    ! Allocate generic arrays for the cell and nodes. These will be
    ! used when a BC is not included in a computed but needs to be
    ! point somehwere.
    allocate(zeroCellVal(isizemax, jsizemax), zeroNodeVal(isizemax, jsizemax), oneCellVal(isizemax, jsizemax))
    oneCellVal = one
    zeroCellVal = zero
    zeroNodeVal = zero

    ! Finally, create the shortcut array for all families. This is just
    ! 1,2,3..totalFamilies.
    allocate(fullFamList(totalFamilies))
    do i=1, totalFamilies
       fullFamList(i) = i
    end do

  end subroutine setSurfaceFamilyInfo

  subroutine createNodeScatterForFamilies(famList, exch, sps, localIndices)

    ! The purpose of this routine is to create the appropriate data
    ! structures that allow for the averaging of cell based surface
    ! quantities to node-based quantities. The primary reason for this
    ! is that the viscous stress tensor is not available at halo cells
    ! and therefore it is not possible to create consistent node-based
    ! values locallly. What the scatter does is allows us to sum the
    ! nodal values across processors, average them and finally update
    ! the node based values to be consistent. This operation is
    ! necessary for several operations:

    ! 1. Integration of forces over zipper triangles requires force/area
    ! at nodes.
    ! 2. Lift distributions/slices also requires node-based tractions
    ! 3. Node-based output for tecplot files.
    use constants
    use communication, only : adflow_comm_world, myid, nProc
    use surfaceFamilies, only : familyExchange, IS1, IS2!, PETSC_COPY_VALUES, PETSC_DETERMINE
    use utils, only : pointReduce, eChk
    use surfaceUtils
    implicit none

    ! Input Parameters
    integer(kind=intType) , dimension(:), intent(in) :: famList
    integer(kind=intType) , intent(in) :: sps
    type(familyExchange), intent(inout) :: exch
    integer(kind=intType), dimension(:), intent(out), allocatable :: localIndices

    ! Working param
    integer(kind=intType) :: i,j, ierr, nNodesLocal, nNodesTotal, nCellsLocal, nFam
    integer(kind=intType) :: nUnique, iSize, iStart, iEnd, iProc
    real(kind=realType), dimension(:, :), allocatable :: localNodes, allNodes
    real(kind=realType), dimension(:, :), allocatable :: uniqueNodes
    integer(kind=intType), dimension(:), allocatable :: link, startIndices, endIndices
    integer(kind=intType), dimension(:), allocatable :: nNodesProc, cumNodesProc
    real(kind=realType) :: tol
    integer(kind=intType) :: mpiStatus(MPI_STATUS_SIZE)

    ! Save the family list.
    nFam = size(famList)
    allocate(exch%famList(nFam))
    exch%famList = famList
    exch%sps = sps

    ! Determine the total number of nodes and cells on this
    ! processor. This will include the zipper mesh if there is one.
    call getSurfaceSize(nNodesLocal, nCellsLocal, famList, nFam, .True.)

    ! Allocate the space to store nodal values, connectivity and the
    ! family of each element
    exch%nNodes = nNodesLocal

    ! Allocate space for the some arrays
    allocate(localNodes(3, nNodesLocal), nNodesProc(nProc), cumNodesProc(0:nProc))
    call getSurfacePoints(localNodes, nNodesLocal, sps, famList, nFam, .True.)

    ! Determine the total number of nodes on each proc
    call mpi_allgather(nNodesLocal, 1, adflow_integer, nNodesProc, 1, adflow_integer, &
         adflow_comm_world, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    ! Determine cumulative version
    cumNodesProc(0) = 0_intType
    nNodesTotal = 0
    do i=1, nProc
       nNodesTotal = nNodesTotal + nNodesProc(i)
       cumNodesProc(i) = cumNodesProc(i-1) + nNodesProc(i)
    end do

    ! Send all the nodes to everyone
    allocate(allNodes(3, nNodesTotal))
    call mpi_allgatherv(localNodes, nNodesLocal*3, adflow_real, allNodes, &
         nNodesProc*3, cumNodesProc*3, adflow_real, adflow_comm_world, ierr)
    call EChk(ierr, __FILE__, __LINE__)

    ! Local nodes is no longer necessary
    deallocate(localNodes)

    ! Now point reduce
    allocate(uniqueNodes(3, nNodestotal), link(nNodestotal))
    tol = 1e-12

    call pointReduce(allNodes, nNodesTotal, tol, uniqueNodes, link, nUnique)

    ! We can immediately discard everything but link since we are only
    ! doing logical operations here:
    deallocate(uniqueNodes, allNodes)

    ! Now back out the global indices for our local points
    if (allocated(localIndices)) then
       deallocate(localIndices)
    end if
    allocate(localIndices(nNodesLocal))
    do i=1, nNodesLocal
       ! The -1 is to convert to 0-based ordering for petsc
       localIndices(i) = link(cumNodesProc(myid) + i)-1
    end do

    ! Create the basic (scalar) local vector
    call VecCreateMPI(ADFLOW_COMM_WORLD, nNodesLocal, PETSC_DETERMINE, &
         exch%nodeValLocal, ierr)
    call EChk(ierr,__FILE__,__LINE__)

    ! Create the basic global vector. This is slightly tricker than it
    ! sounds. We could just make it uniform, but then there would be
    ! more communicaiton than necessary. Instead what we do is determine
    ! the min and max range of local indices on the proc and the one
    ! before it. A little diagram will help
    !
    ! Proc 0 +---------------------+
    ! Proc 1                    +-------------+
    ! Proc 2                                +----------------+
    !
    ! Proc zero has a many global nodes as local since they are by
    ! definition all unqiue. Proc 1 then will start at 1 more than the
    ! proc 0 and continue to it's maximum value. Proc 2 starts at the
    ! end of proc 1 etc. This way the vast majority of the global nodes
    ! are owned locally.

    ! In order to determine the owning range for each processor, it is
    ! much trickier than it sounds. We do a linear cascasde through the
    ! procs sending the upper range from proc 0 to proc 1, then proc1 to
    ! proc 2 and so on.

    ! Proc zero owns all of it's nodes.
    if (myid == 0) then
       iStart = 0
       if (nNodesLocal == 0) then
          iEnd = 0
       else
          iEnd = maxval(localIndices) + 1
       end if
    end if

    do iProc=0, nProc-2
       if (myid == iProc) then
          ! I need to send my iEnd to proc+1
          call mpi_send(iEnd, 1, adflow_integer, iProc+1, iProc, adflow_comm_world, ierr)
          call EChk(ierr,__FILE__,__LINE__)
       else if(myid == iProc+1) then

          ! Receive the value from the proc below me:
          call mpi_recv(iEnd, 1, adflow_integer, iProc, iProc, adflow_comm_world, mpiStatus, ierr)
          call EChk(ierr,__FILE__,__LINE__)

          ! On this proc, the start index is the
          iStart = iEnd
          if (nNodesLOCAl == 0) then
             iEnd = iStart
          else
             iEnd = max(iStart, maxval(localIndices)+1)
          end if
       end if
    end do

    iSize = iEnd-iStart
    ! Create the actual global vec. Note we also include nUnique to make
    ! sure we have all the local sizes correct.
    call VecCreateMPI(ADFLOW_COMM_WORLD, iSize, nUnique, &
         exch%nodeValGlobal, ierr)
    call EChk(ierr,__FILE__,__LINE__)

    call VecDuplicate(exch%nodeValGlobal, exch%sumGlobal, ierr)
    call EChk(ierr,__FILE__,__LINE__)

    ! Now create the scatter that goes from the local vec to the global
    ! vec.

    ! Indices for the local vector is just a stride, starting at the
    ! offset
    call ISCreateStride(ADFLOW_COMM_WORLD, nNodesLocal, cumNodesProc(myid), &
         1, IS1, ierr)
    call EChk(ierr,__FILE__,__LINE__)

    ! Indices for the global vector are the "localIndices" we previously
    ! computed.
    call ISCreateGeneral(adflow_comm_world, nNodesLocal, localIndices, &
         PETSC_COPY_VALUES, IS2, ierr)
    call EChk(ierr,__FILE__,__LINE__)

    call VecScatterCreate(exch%nodeValLocal, IS1, exch%nodeValGlobal, IS2, &
         exch%scatter, ierr)
    call EChk(ierr,__FILE__,__LINE__)

    ! And dont' forget to destroy the index sets
    call ISDestroy(IS1, ierr)
    call EChk(ierr,__FILE__,__LINE__)

    call ISDestroy(IS2, ierr)
    call EChk(ierr,__FILE__,__LINE__)

    exch%allocated = .True.
  end subroutine createNodeScatterForFamilies


  subroutine setReferenceVolume(level)

    use constants
    use blockPointers, only : nDom, flowDoms, ib, jb, kb
    use inputTimeSpectral, only : nTimeIntervalsSpectral
    use utils, only : setPointers
    implicit none
    integer :: ierr

    integer(kind=intType), intent(in) :: level
    integer(kind=intType) :: nn, sps
    integer(kind=intType) :: i,j, k

    spectral: do sps=1,nTimeIntervalsSpectral
       domains: do nn=1,nDom
          call setPointers(nn, level, sps)
          allocate(flowDoms(nn, level, sps)%volRef(0:ib, 0:jb, 0:kb))

          do k=0, kb
             do j=0, jb
                do i=0, ib
                   flowDoms(nn, level, sps)%volRef(i, j, k) = &
                        flowDoms(nn, level, sps)%vol(i, j, k)
                end do
             end do
          end do
       end do domains
    end do spectral
  end subroutine setReferenceVolume

  subroutine setGlobalCellsAndNodes(level)
    !
    !      Determine the global node numbering that is used to assemble
    !      the adjoint system of equations. It take cares of all the halo
    !      nodes between the blocks.
    !      The nodes are numbered according to the following sequence:
    !      loop processor = 1, nProc
    !        loop domain = 1, nDom
    !          loop k = 2, kl
    !            loop j = 2, jl
    !              loop i = 2, il
    !      Only the onwned nodes are numbered, meaning i/j/k span from 2
    !      to il/jl/kl. The halo nodes receive the numbering from the
    !      neighboring block that owns them.
    !      These variables are the same for all spectral modes, therefore
    !      only the 1st mode needs to be communicated.
    !      This function will also set FMPointer which is only defined
    !      on wall boundary conditions and points to the correct index
    !      for the vectors that are of shape nsurface nodes
    !
    use ADjointVars
    use blockpointers
    use communication
    use inputTimeSpectral
    use utils, only: setPointers, terminate
    use haloExchange, only : whalo1to1intgeneric
    implicit none

    ! Input variables
    integer(kind=intType), intent(in) :: level

    ! Local variables
    integer(kind=intType) :: nn, i, j, k, sps, iDim
    integer(kind=intType) :: ierr, istart
    logical :: commPressure, commLamVis, commEddyVis, commGamma
    integer(kind=intType), dimension(nProc) :: nNodes, nCells, nCellOffset, nNodeOffset
    integer(kind=intType), dimension(nDom) :: nCellBLockOffset,nNodeBLockOffset
    integer(kind=intType) :: npts, nCell, nNode
    integer(kind=intType), dimension(:), allocatable :: nNodesProc, cumNodesProc
    integer(kind=intTYpe), dimension(:), allocatable :: nCellsProc, cumCellsProc
    integer(kind=intType) :: iBeg, iEnd, jBeg, jEnd, ii, jj,mm

    do sps=1, nTimeIntervalsSpectral
       do nn=1, nDom
          call setPointers(nn, level, sps)
          ! Allocate memory for the cell and node indexing...only on sps=1
          allocate(flowDoms(nn,level,sps)%globalCell(0:ib,0:jb,0:kb), &
               flowDoms(nn,level,sps)%globalNode(0:ie,0:je,0:ke), stat=ierr)
          if (ierr /=0) then
             call terminate("setGlobalCellsAndNodes", "Allocation failure for globalCell/Node")
          end if
          ! Assign a 'magic number' of -5 to globalCell and global Node:
          flowDoms(nn,level,sps)%globalCell = -5
          flowDoms(nn,level,sps)%globalNode = -5
       end do
    end do

    ! Determine the number of nodes and cells owned by each processor
    ! by looping over the local block domains.
    nCellsLocal(level) = 0
    nNodesLocal(level) = 0
    do nn=1,nDom
       ! Set to first spectral instance since we only need sizes
       call setPointers(nn, level, 1_intType)
       nCellsLocal(level) = nCellsLocal(level) + nx*ny*nz
       nNodesLocal(level) = nNodesLocal(level) + il*jl*kl
    enddo

    ! Reduce the number of cells in all processors: add up nCellsLocal
    ! into nCellsGlobal and sends the result to all processors.
    ! (use mpi sum operation)

    call mpi_allreduce(nCellsLocal(level), nCellsGlobal(level), 1, adflow_integer, &
         mpi_sum, ADflow_comm_world, ierr)

    ! Gather the number of Cells per processor in the root processor.
    call mpi_gather(nCellsLocal(level), 1, adflow_integer, nCells, 1, &
         adflow_integer, 0, ADflow_comm_world, ierr)

    ! Repeat for the number of nodes.
    ! (use mpi sum operation)
    call mpi_allreduce(nNodesLocal(level), nNodesGlobal(level), 1, adflow_integer, &
         mpi_sum, ADflow_comm_world, ierr)

    ! Gather the number of nodes per processor in the root processor.
    call mpi_gather(nNodesLocal(level), 1, adflow_integer, nNodes, 1, &
         adflow_integer, 0, ADflow_comm_world, ierr)

    ! Determine the global cell number offset for each processor.
    rootProc: if( myID==0) then
       nCellOffset(1) = 0
       nNodeOffset(1) = 0
       do nn=2,nProc
          nCellOffset(nn) = nCellOffset(nn-1) + nCells(nn-1)
          nNodeOffset(nn) = nNodeOffset(nn-1) + nNodes(nn-1)
       enddo
    endif rootProc

    ! Scatter the global cell number offset per processor.
    call mpi_scatter(nCellOffset, 1, adflow_integer, nCellOffsetLocal(level), 1, &
         adflow_integer, 0, ADflow_comm_world, ierr)

    ! Determine the global cell number offset for each local block.
    nCellBlockOffset(1) = nCellOffsetLocal(level)
    do nn=2,nDom
       call setPointers(nn-1, level, 1)
       nCellBlockOffset(nn) = nCellBlockOffset(nn-1)          &
            + nx*ny*nz
    enddo

    ! Repeat for nodes.
    call mpi_scatter(nNodeOffset, 1, adflow_integer, nNodeOffsetLocal(level), 1, &
         adflow_integer, 0, ADflow_comm_world, ierr)

    ! Determine the global node number offset for each local block.
    nNodeBlockOffset(1) = nNodeOffsetLocal(level)
    do nn=2,nDom
       call setPointers(nn-1, level, 1)
       nNodeBlockOffset(nn) = nNodeBLockOffset(nn-1) + il*jl*kl
    enddo

    ! Determine the global block row index for each (i,j,k) cell in
    ! each local block.

    do nn=1, nDom
       do sps=1, nTimeIntervalsSpectral
          call setPointers(nn, level, sps)
          do k=2, kl
             do j=2, jl
                do i=2, il
                   ! modified Timespectral indexing. Put all time
                   ! instances of a give block adjacent to each other in
                   ! the matrix
                   globalCell(i, j, k) = &
                        nCellBLockOffset(nn)*nTimeIntervalsSpectral+nx*ny*nz*(sps-1)+&
                        (i-2) +(j-2)*nx +(k-2)*nx*ny
                enddo
             enddo
          enddo
       enddo
    end do

    ! Determine the global block row index for each (i,j,k) node in
    ! each local block.
    do sps=1, nTimeIntervalsSpectral
       do nn=1, nDom
          call setPointers(nn, level, sps)
          do k=1, kl
             do j=1, jl
                do i=1, il
                   !modified Timespectral indexing. Put all time
                   !instances of a give block adjacent to each other in
                   !the matrix
                   globalNode(i, j, k) = &
                        nNodeBLockOffset(nn)*nTimeIntervalsSpectral + &
                        il*jl*kl*(sps-1) + (i-1)+(j-1)*il + (k-1)*il*jl

                end do
             end do
          end do
       end do
    end do

    ! The above procedure has uniquely numbered all cells and nodes
    ! owned on each processor. However we must also determine the
    ! indices of the halo cells/nodes from other processors. To do this
    ! we just run the specific halo exchanges for the cells and one for
    ! the nodes

    spectralModes: do sps=1,nTimeIntervalsSpectral
       domainLoop:do nn=1, nDom
          flowDoms(nn, level, sps)%intCommVars(1)%var => &
               flowDoms(nn, level, sps)%globalNode(:, :, :)
       end do domainLoop

       ! Run the generic integer exchange
       call wHalo1to1IntGeneric(1, level, sps, commPatternNode_1st, internalNode_1st)
    end do spectralModes

    spectralModes2: do sps=1,nTimeIntervalsSpectral
       domainLoop2:do nn=1, nDom
          flowDoms(nn, level, sps)%intCommVars(1)%var => &
               flowDoms(nn, level, sps)%globalCell(:, :, :)
       end do domainLoop2

       ! Run the generic integer exchange
       call wHalo1to1IntGeneric(1, level, sps, commPatternCell_2nd, internalCell_2nd)
    end do spectralModes2

  end subroutine setGlobalCellsAndNodes
  subroutine setFamilyInfoFaces(level)
    !
    !       setFamilyInfoFaces sets the values of the family parameters
    !       for faces on the given multigrid level. The default values for
    !       indFamily is 0, which means that the mass flow through that
    !       face does not contribute to the mass flow that must be
    !       monitored. For sliding mesh interfaces both sides of the
    !       interface are monitored and the value of indFamily corresponds
    !       to one of the two entries in the local monitoring arrays. The
    !       values of factFamily are such that the mass flow entering the
    !       block is defined positive.
    !       Note that only the 1st spectral solution is treated, because
    !       this informations is the same for all of them.
    !
    use constants
    use blockPointers
    use cgnsGrid
    use inputTimeSpectral
    use monitor
    use section
    use utils, only : setPointers, terminate
    implicit none
    !
    !      Subroutine arguments.
    !
    integer(kind=intType), intent(in) :: level
    !
    !      Local variables.
    !
    integer :: ierr

    integer(kind=intType) :: nn, mm, i, j, k, ii, nSlices
    integer(kind=intType) :: iBeg, iEnd, jBeg, jEnd, kBeg, kEnd

    integer(kind=intType), dimension(cgnsNFamilies) :: orToMassFam

    ! Set the offset ii for the values of orToMassFam. If the mass
    ! flow through sliding mesh interfaces must be monitored this
    ! offset if 2*cgnsNSliding. This means that in the arrays to store
    ! the mass flow the sliding mesh interfaces are stored first,
    ! followed by the families.

    if( monMassSliding ) then
       ii = 2*cgnsNSliding
    else
       ii = 0
    endif

    ! Determine the number of families for which the mass flow must
    ! be monitored and set the entries of orToMassFam accordingly,
    ! i.e. the offset ii is included.

    mm = ii
    do nn=1,cgnsNFamilies
       if(cgnsFamilies(nn)%monitorMassflow            .and. &
            cgnsFamilies(nn)%BCType /= MassBleedInflow  .and. &
            cgnsFamilies(nn)%BCType /= MassBleedOutflow .and. &
            cgnsFamilies(nn)%BCType /= SlidingInterface) then
          mm = mm + 1
          orToMassFam(nn) = mm
       else
          orToMassFam(nn) = 0
       endif
    enddo

    ! Set monMassFamilies to .true. if the mass flow of at least one
    ! family must be monitored. Otherwise set it to .false.

    if(mm > ii) then
       monMassFamilies = .true.
    else
       monMassFamilies = .false.
    endif

    ! If this is the first level, allocate the memory for
    ! massFlowFamilyInv and massFlowFamilyDiss.

    if(level == 1) then
       nn = nTimeIntervalsSpectral

       allocate(massFlowFamilyInv(0:mm,nn), &
            massFlowFamilyDiss(0:mm,nn), stat=ierr)
       if(ierr /= 0) &
            call terminate("setFamilyInfoFaces", &
            "Memory allocation failure for &
            &massFlowFamilyInv and massFlowFamilyDiss")
    endif

    ! Loop over the number of domains.

    domains: do nn=1,nDom

       ! Allocate the memory for indFamily and factFamily.

       il = flowDoms(nn,level,1)%il
       jl = flowDoms(nn,level,1)%jl
       kl = flowDoms(nn,level,1)%kl

       allocate(flowDoms(nn,level,1)%indFamilyI (1:il,2:jl,2:kl), &
            flowDoms(nn,level,1)%indFamilyJ (2:il,1:jl,2:kl), &
            flowDoms(nn,level,1)%indFamilyK (2:il,2:jl,1:kl), &
            flowDoms(nn,level,1)%factFamilyI(1:il,2:jl,2:kl), &
            flowDoms(nn,level,1)%factFamilyJ(2:il,1:jl,2:kl), &
            flowDoms(nn,level,1)%factFamilyK(2:il,2:jl,1:kl), &
            stat=ierr)
       if(ierr /= 0)                          &
            call terminate("setFamilyInfoFaces", &
            "Memory allocation failure for indFamily &
            &and factFamily")

       ! Set the pointers for this domain.

       call setPointers(nn, level, 1_intType)

       ! Determine the number of slices for this block to make
       ! the full wheel.

       nSlices = sections(sectionID)%nSlices

       ! Initialize the values of indFamily and factFamily.

       indFamilyI = 0_intType
       indFamilyJ = 0_intType
       indFamilyK = 0_intType

       factFamilyI = 0_intType
       factFamilyJ = 0_intType
       factFamilyK = 0_intType

       ! Loop over the boundary conditions.

       boco: do mm=1,nBocos

          ! Test for the boundary condition.

          select case (BCType(mm))
          case (SlidingInterface)

             ! Sliding mesh boundary.
             ! If the mass flow through sliding interfaces must be monitored,
             ! determine the index in the arrays massFlowFamilyInv and
             ! massFlowFamilyDiss where to store the contribution of this
             ! subface. If the sliding mesh mass flows are not monitored,
             ! set the index to 0.

             if( monMassSliding ) then
                ii = 2*abs(groupNum(mm))
                if(groupNum(mm) < 0) ii = ii - 1
             else
                ii = 0
             endif

          case (MassBleedInflow, MassBleedOutflow)

             ! Inflow or outflow bleed. These boundary conditions are
             ! handled separetely and need not be monitored.

             ii = 0

          case default

             ! Subface is an ordinary boundary condition. Determine the
             ! family ID and set the index ii in the arrays massFlowFamilyInv
             ! and massFlowFamilyDiss accordingly.

             if(groupNum(mm) > 0) then
                ii = orToMassFam(groupNum(mm))
             else
                ii = 0
             endif

          end select

          ! Set the owned cell range for the faces on this subface.
          ! As icBeg, etc. may contain halo cells, inBeg, etc. is
          ! used.

          iBeg = min(inBeg(mm), inEnd(mm)) +1
          iEnd = max(inBeg(mm), inEnd(mm))

          jBeg = min(jnBeg(mm), jnEnd(mm)) +1
          jEnd = max(jnBeg(mm), jnEnd(mm))

          kBeg = min(knBeg(mm), knEnd(mm)) +1
          kEnd = max(knBeg(mm), knEnd(mm))

          ! Determine the block this subface is located on and set
          ! the corresponding values of indFamily and factFamily.
          ! Note that factFamily is set to nSlices on min faces and to
          ! -nSlices on max faces, such that the mass flow entering the
          ! domain is defined positive and the mass flow of the entire
          ! wheel is monitored.

          select case( BCFaceID(mm) )
          case (iMin)
             do k=kBeg,kEnd
                do j=jBeg,jEnd
                   indFamilyI (1,j,k) = ii
                   factFamilyI(1,j,k) = nSlices
                enddo
             enddo

             !===========================================================

          case (iMax)
             do k=kBeg,kEnd
                do j=jBeg,jEnd
                   indFamilyI (il,j,k) = ii
                   factFamilyI(il,j,k) = -nSlices
                enddo
             enddo

             !===========================================================

          case (jMin)
             do k=kBeg,kEnd
                do i=iBeg,iEnd
                   indFamilyJ (i,1,k) = ii
                   factFamilyJ(i,1,k) = nSlices
                enddo
             enddo

             !===========================================================

          case (jMax)
             do k=kBeg,kEnd
                do i=iBeg,iEnd
                   indFamilyJ (i,jl,k) = ii
                   factFamilyJ(i,jl,k) = -nSlices
                enddo
             enddo

             !===========================================================

          case (kMin)
             do j=jBeg,jEnd
                do i=iBeg,iEnd
                   indFamilyK (i,j,1) = ii
                   factFamilyK(i,j,1) = nSlices
                enddo
             enddo

             !===========================================================

          case (kMax)
             do j=jBeg,jEnd
                do i=iBeg,iEnd
                   indFamilyK (i,j,kl) = ii
                   factFamilyK(i,j,kl) = -nSlices
                enddo
             enddo

          end select

       enddo boco
    enddo domains

  end subroutine setFamilyInfoFaces
  subroutine shiftCoorAndVolumes
    !
    !       shiftCoorAndVolumes shifts the owned coordinates and
    !       volumes in case of a deforming mesh for an unsteady
    !       computation. In this case the old coordinates are needed to
    !       determine the mesh velocities. The loop over the number of
    !       spectral solutions is present for consistency, but this number
    !       will be 1 when this routine is called.
    !
    use blockPointers
    use inputTimeSpectral
    use iteration
    use utils, only : setPointers
    implicit none
    !
    !      Local variables.
    !
    integer(kind=intType) :: i, j, k, nn, mm, ll, kk

    ! Loop over the number of spectral solutions and local blocks.

    spectralLoop: do kk=1,nTimeIntervalsSpectral
       domains: do nn=1,nDom

          ! Set the pointers for this block on the ground level.

          call setPointers(nn, groundLevel,kk)

          ! Shift the coordinates already stored in xOld and the
          ! volumes stored in volOld.

          loopOldLevels: do mm=nOldLevels,2,-1

             ! Shift the coordinates from level mm-1 to mm, including
             ! the halo's.

             ll = mm - 1

             do k=0,ke
                do j=0,je
                   do i=0,ie
                      xOld(mm,i,j,k,1) = xOld(ll,i,j,k,1)
                      xOld(mm,i,j,k,2) = xOld(ll,i,j,k,2)
                      xOld(mm,i,j,k,3) = xOld(ll,i,j,k,3)
                   enddo
                enddo
             enddo

             ! Shift the old volumes from level mm-1 to mm.
             ! Only the owned ones need to be considered.

             do k=2,kl
                do j=2,jl
                   do i=2,il
                      volOld(mm,i,j,k) = volOld(ll,i,j,k)
                   enddo
                enddo
             enddo

          enddo loopOldLevels

          ! Shift the current coordinates into the 1st level of xOld.

          do k=0,ke
             do j=0,je
                do i=0,ie
                   xOld(1,i,j,k,1) = x(i,j,k,1)
                   xOld(1,i,j,k,2) = x(i,j,k,2)
                   xOld(1,i,j,k,3) = x(i,j,k,3)
                enddo
             enddo
          enddo

          ! Shift the current volumes into the 1st level of volOld.

          do k=2,kl
             do j=2,jl
                do i=2,il
                   volOld(1,i,j,k) = vol(i,j,k)
                enddo
             enddo
          enddo

       enddo domains
    enddo spectralLoop

  end subroutine shiftCoorAndVolumes
  subroutine viscSubfaceInfo(level)
    !
    !       viscSubfaceInfo allocates the memory for the storage of the
    !       stress tensor and heat flux vector of viscous subfaces for the
    !       given multigrid level and all spectral solutions. Furthermore
    !       the pointers viscIminPointer, etc. Are allocated and set.
    !       These pointers contain info to which viscous subface the faces
    !       of the block faces possibly belong. If not part of a viscous
    !       subface these values are set to 0. Note that these pointers
    !       are only allocated and determined for the 1st spectral
    !       solution, because the info is the same for all of them.
    !
    use constants
    use blockPointers
    use inputTimeSpectral
    use utils, only : setPointers, terminate
    implicit none
    !
    !      Subroutine argument.
    !
    integer(kind=intType), intent(in) :: level
    !
    !      Local variables.
    !
    integer :: ierr

    integer(kind=intType) :: nn, mm, sps, i, j
    integer(kind=intType) :: iBeg, iEnd, jBeg, jEnd

    integer(kind=intType), dimension(:,:), pointer :: viscPointer

    ! Loop over the number blocks stored on this processor.

    domains: do nn=1,nDom

       ! Set the pointers to the block of the 1st spectral solution.

       call setPointers(nn, level, 1_intType)

       ! Allocate the memory for viscSubface and the pointers
       ! viscIminPointer, etc. ViscSubface must be allocated for
       ! all spectral solutions, the pointers only for the 1st.

       do sps=1,nTimeIntervalsSpectral
          allocate(flowDoms(nn,level,sps)%viscSubface(nViscBocos), &
               stat=ierr)
          if(ierr /= 0)                         &
               call terminate("viscSubfaceInfo", &
               "Memory allocation failure for viscSubface")
       enddo

       allocate(flowDoms(nn,level,1)%viscIminPointer(2:jl,2:kl), &
            flowDoms(nn,level,1)%viscImaxPointer(2:jl,2:kl), &
            flowDoms(nn,level,1)%viscJminPointer(2:il,2:kl), &
            flowDoms(nn,level,1)%viscJmaxPointer(2:il,2:kl), &
            flowDoms(nn,level,1)%viscKminPointer(2:il,2:jl), &
            flowDoms(nn,level,1)%viscKmaxPointer(2:il,2:jl), &
            stat=ierr)
       if(ierr /= 0)                         &
            call terminate("viscSubfaceInfo", &
            "Memory allocation failure for subface info")

       ! Reset the pointers viscIminPointer, etc. to make it more
       ! readable and initialize them to 0. This indicates that
       ! the faces are not part of a viscous wall subfaces.

       viscIminPointer => flowDoms(nn,level,1)%viscIminPointer
       viscImaxPointer => flowDoms(nn,level,1)%viscImaxPointer
       viscJminPointer => flowDoms(nn,level,1)%viscJminPointer
       viscJmaxPointer => flowDoms(nn,level,1)%viscJmaxPointer
       viscKminPointer => flowDoms(nn,level,1)%viscKminPointer
       viscKmaxPointer => flowDoms(nn,level,1)%viscKmaxPointer

       viscIminPointer = 0
       viscImaxPointer = 0
       viscJminPointer = 0
       viscJmaxPointer = 0
       viscKminPointer = 0
       viscKmaxPointer = 0

       ! Loop over the viscous subfaces to allocate the memory for the
       ! stress tensor and the heat flux vector and to set the range
       ! in viscIminPointer, etc.

       viscSubfaces: do mm=1,nViscBocos

          ! Store the cell range in iBeg, iEnd, etc. As the viscous data
          ! do not allow for an overlap, the nodal range of the
          ! subface must be used.

          iBeg = BCData(mm)%inBeg + 1
          iEnd = BCData(mm)%inEnd

          jBeg = BCData(mm)%jnBeg + 1
          jEnd = BCData(mm)%jnEnd

          ! Loop over the spectral solutions and allocate the memory
          ! for the stress tensor, heat flux and friction velocity.

          do sps=1,nTimeIntervalsSpectral

             ! Set the pointer for viscSubface to make the code
             ! more readable and allocate the memory.

             viscSubface => flowDoms(nn,level,sps)%viscSubface

             allocate(viscSubface(mm)%tau( iBeg:iEnd,jBeg:jEnd,6), &
                  viscSubface(mm)%q(   iBeg:iEnd,jBeg:jEnd,3), &
                  viscSubface(mm)%utau(iBeg:iEnd,jBeg:jEnd),   &
                  stat=ierr)
             if(ierr /= 0)                       &
                  call terminate("viscSubfaceInfo", &
                  "Memory allocation failure for tau, q &
                  &and utau.")
          enddo

          ! Set the pointer viscPointer, depending on the block face
          ! on which the subface is located.

          select case (BCFaceID(mm))
          case (iMin)
             viscPointer => viscIminPointer

          case (iMax)
             viscPointer => viscImaxPointer

          case (jMin)
             viscPointer => viscJminPointer

          case (jMax)
             viscPointer => viscJmaxPointer

          case (kMin)
             viscPointer => viscKminPointer

          case (kMax)
             viscPointer => viscKmaxPointer
          end select

          ! Set this range in viscPointer to viscous subface mm.

          do j=jBeg,jEnd
             do i=iBeg,iEnd
                viscPointer(i,j) = mm
             enddo
          enddo

       enddo viscSubfaces

    enddo domains

  end subroutine viscSubfaceInfo

  !      ==================================================================

  subroutine allocateMetric(level)
    !
    !       allocateMetric allocates the memory for the metric variables
    !       on the given grid level for all spectral solutions.
    !
    use constants
    use block
    use inputPhysics
    use inputTimeSpectral
    use iteration
    use inputUnsteady
    use utils, only : terminate, setPointers
    implicit none
    !
    !      Subroutine arguments.
    !
    integer(kind=intType), intent(in) :: level
    !
    !      Local variables.
    !
    integer :: ierr

    integer(kind=intType) :: nn, mm, sps
    integer(kind=intType) :: il, jl, kl, ie, je, ke, ib, jb, kb

    type(BCDataType), dimension(:), pointer :: BCData

    ! Loop over the number of spectral solutions and local blocks.

    spectral: do sps=1,nTimeIntervalsSpectral
       domains: do nn=1,nDom

          ! Store the the upper boundaries of the block a bit easier.

          il = flowDoms(nn,level,sps)%il
          jl = flowDoms(nn,level,sps)%jl
          kl = flowDoms(nn,level,sps)%kl

          ie = flowDoms(nn,level,sps)%ie
          je = flowDoms(nn,level,sps)%je
          ke = flowDoms(nn,level,sps)%ke

          ib = flowDoms(nn,level,sps)%ib
          jb = flowDoms(nn,level,sps)%jb
          kb = flowDoms(nn,level,sps)%kb

          ! Allocate the memory for the volumes and the face normals.

          allocate(flowDoms(nn,level,sps)%si(0:ie,1:je,1:ke,3), &
               flowDoms(nn,level,sps)%sj(1:ie,0:je,1:ke,3), &
               flowDoms(nn,level,sps)%sk(1:ie,1:je,0:ke,3), &
               flowDoms(nn,level,sps)%vol(0:ib,0:jb,0:kb),  &
               stat=ierr)
          if(ierr /= 0)                      &
               call terminate("allocateMetric", &
               "Memory allocation failure for &
               &normals and volumes")

          ! Added by HDN
          ! Added s[I,J,K]ALE
          if (equationMode == unSteady .and. useALE) then
             allocate(flowDoms(nn,level,sps)%sIALE(0:nALEsteps,0:ie,1:je,1:ke,3), &
                  flowDoms(nn,level,sps)%sJALE(0:nALEsteps,1:ie,0:je,1:ke,3), &
                  flowDoms(nn,level,sps)%sKALE(0:nALEsteps,1:ie,1:je,0:ke,3), &
                  stat=ierr)
             if(ierr /= 0)                      &
                  call terminate("allocateMetric", &
                  "Memory allocation failure for &
                  &sIALE, sJALE, and sKALE")
          end if

          ! Allocate the memory for the unit normals of the boundary
          ! faces. First set the pointer to make it more readable.

          BCData => flowDoms(nn,level,sps)%BCData

          do mm=1,flowDoms(nn,level,sps)%nBocos

             ! Store the size of the subface in ie:ib and je:jb, because
             ! these variables are not needed anymore.

             ie = BCData(mm)%icBeg
             ib = BCData(mm)%icEnd
             je = BCData(mm)%jcBeg
             jb = BCData(mm)%jcEnd

             ! Allocate the memory for the unit normals.

             ! Added by HDN
             ! Added normALE
             allocate(                                           &
                  BCData(mm)%norm(ie:ib,je:jb,3),                &
                  BCData(mm)%normALE(0:nALEsteps,ie:ib,je:jb,3), &
                  stat=ierr)
             if(ierr /= 0)                      &
                  call terminate("allocateMetric", &
                  "Memory allocation failure for norm")
          enddo

          ! Allocate the memory for the old volumes; only for unsteady
          ! problems on deforming meshes on the finest grid level.

          if(level == 1 .and. deforming_Grid .and. &
               equationMode == unsteady) then

             allocate(                                                   &
                  flowDoms(nn,level,sps)%volOld(nOldLevels,2:il,2:jl,2:kl), &
                  stat=ierr)
             if(ierr /= 0)                      &
                  call terminate("allocateMetric", &
                  "Memory allocation failure for volOld")
          endif

       enddo domains
    enddo spectral

  end subroutine allocateMetric

  !      ==================================================================

  subroutine metric(level)
    !
    !       metric computes the face normals and the volume for the given
    !       grid level for all spectral solutions. First the volumes are
    !       computed assuming that the block is right handed. Then the
    !       number of positive and negative volumes are determined. If all
    !       volumes are positive the block is indeed right handed; if all
    !       volumes are negative the block is left handed and both the
    !       volumes and the normals must be negated (for the normals this
    !       is done by the introduction of fact, which is either -0.5 or
    !       0.5); if there are both positive and negative volumes the mesh
    !       is not valid.
    !
    use constants
    use blockPointers
    use cgnsGrid
    use communication
    use inputTimeSpectral
    use checkVolBlock
    use inputIteration
    use utils, only : setPointers, terminate, returnFail
    implicit none
    !
    !      Subroutine arguments.
    !
    integer(kind=intType), intent(in) :: level
    !
    !      Local parameter.
    !
    real(kind=realType), parameter :: thresVolume = 1.e-2_realType
    real(kind=realType), parameter :: haloCellRatio = 1e-10_realType
    !
    !      Local variables.
    !
    integer :: ierr

    integer(kind=intType) :: i, j, k, n, m, l
    integer(kind=intType) :: nn, mm, sps
    integer(kind=intType) :: nVolNeg,   nVolPos
    integer(kind=intType) :: nVolBad,   nVolBadGlobal
    integer(kind=intType) :: nBlockBad, nBlockBadGlobal

    real(kind=realType) :: fact, mult
    real(kind=realType) :: xp, yp, zp, vp1, vp2, vp3, vp4, vp5, vp6

    real(kind=realType), dimension(3) :: v1, v2

    real(kind=realType), dimension(:,:,:), pointer :: ss

    character(len=10) :: integerString

    logical :: checkK, checkJ, checkI, checkAll, checkBlank
    logical :: badVolume, iBlankAllocated

    logical, dimension(:,:,:), pointer :: volumeIsNeg

    type(checkVolBlockType), &
         dimension(nDom,nTimeIntervalsSpectral) :: checkVolDoms

    ! Initialize the number of bad volumes and bad blocks to 0.

    nVolBad   = 0
    nBlockBad = 0

    ! Loop over the number of spectral solutions and local blocks.

    spectral: do sps=1,nTimeIntervalsSpectral
       domains: do nn=1,nDom

          ! Set the pointers to this block and allocate the memory for
          ! volumeIsNeg. Set a pointer to this entry afterwards to make
          ! the code more readable.

          call setPointers(nn, level, sps)
          if (associated(flowDoms(nn, level, sps)%iblank)) then
             iBlankAllocated = .True.
          else
             iBlankAllocated = .False.
          end if

          allocate(checkVolDoms(nn,sps)%volumeIsNeg(2:il,2:jl,2:kl), &
               stat=ierr)
          if(ierr /= 0)              &
               call terminate("metric", &
               "Memory allocation failure for volumeIsNeg")
          volumeIsNeg => checkVolDoms(nn,sps)%volumeIsNeg
          !
          !           Volume and block orientation computation.
          !
          ! Initialize the number of positive and negative volumes for
          ! this block to 0.

          nVolNeg = 0
          nVolPos = 0

          ! Compute the volumes. The hexahedron is split into 6 pyramids
          ! whose volumes are computed. The volume is positive for a
          ! right handed block.
          ! Initialize the volumes to zero. The reasons is that the second
          ! level halo's must be initialized to zero and for convenience
          ! all the volumes are set to zero.

          vol = zero

          do k=1,ke
             n = k -1

             checkK = .true.
             if(k == 1 .or. k == ke) checkK = .false.

             do j=1,je
                m = j -1

                checkJ = .true.
                if(j == 1 .or. j == je) checkJ = .false.

                do i=1,ie
                   l = i -1

                   checkI = .true.
                   if(i == 1 .or. i == ie) checkI = .false.

                   ! Determine whether or not the voluem must be checked for
                   ! quality. Only owned volumes are checked, not halo's.

                   checkAll = .false.

                   ! Only care about the quality of compute cells (1)
                   ! and fringe cells (-1)
                   checkBlank = .False.
                   if (iblankAllocated) then
                      if (abs(iblank(i, j, k)) == 1) then
                         checkBlank = .True.
                      end if
                   end if

                   if (checkK .and. checkJ .and. checkI .and. checkBlank) then
                      checkAll = .true.
                   end if

                   ! Compute the coordinates of the center of gravity.

                   xp = eighth*(x(i,j,k,1) + x(i,m,k,1) &
                        +         x(i,m,n,1) + x(i,j,n,1) &
                        +         x(l,j,k,1) + x(l,m,k,1) &
                        +         x(l,m,n,1) + x(l,j,n,1))
                   yp = eighth*(x(i,j,k,2) + x(i,m,k,2) &
                        +         x(i,m,n,2) + x(i,j,n,2) &
                        +         x(l,j,k,2) + x(l,m,k,2) &
                        +         x(l,m,n,2) + x(l,j,n,2))
                   zp = eighth*(x(i,j,k,3) + x(i,m,k,3) &
                        +         x(i,m,n,3) + x(i,j,n,3) &
                        +         x(l,j,k,3) + x(l,m,k,3) &
                        +         x(l,m,n,3) + x(l,j,n,3))

                   ! Compute the volumes of the 6 sub pyramids. The
                   ! arguments of volpym must be such that for a (regular)
                   ! right handed hexahedron all volumes are positive.

                   vp1 = volpym(x(i,j,k,1), x(i,j,k,2), x(i,j,k,3), &
                        x(i,j,n,1), x(i,j,n,2), x(i,j,n,3), &
                        x(i,m,n,1), x(i,m,n,2), x(i,m,n,3), &
                        x(i,m,k,1), x(i,m,k,2), x(i,m,k,3))

                   vp2 = volpym(x(l,j,k,1), x(l,j,k,2), x(l,j,k,3), &
                        x(l,m,k,1), x(l,m,k,2), x(l,m,k,3), &
                        x(l,m,n,1), x(l,m,n,2), x(l,m,n,3), &
                        x(l,j,n,1), x(l,j,n,2), x(l,j,n,3))

                   vp3 = volpym(x(i,j,k,1), x(i,j,k,2), x(i,j,k,3), &
                        x(l,j,k,1), x(l,j,k,2), x(l,j,k,3), &
                        x(l,j,n,1), x(l,j,n,2), x(l,j,n,3), &
                        x(i,j,n,1), x(i,j,n,2), x(i,j,n,3))

                   vp4 = volpym(x(i,m,k,1), x(i,m,k,2), x(i,m,k,3), &
                        x(i,m,n,1), x(i,m,n,2), x(i,m,n,3), &
                        x(l,m,n,1), x(l,m,n,2), x(l,m,n,3), &
                        x(l,m,k,1), x(l,m,k,2), x(l,m,k,3))

                   vp5 = volpym(x(i,j,k,1), x(i,j,k,2), x(i,j,k,3), &
                        x(i,m,k,1), x(i,m,k,2), x(i,m,k,3), &
                        x(l,m,k,1), x(l,m,k,2), x(l,m,k,3), &
                        x(l,j,k,1), x(l,j,k,2), x(l,j,k,3))

                   vp6 = volpym(x(i,j,n,1), x(i,j,n,2), x(i,j,n,3), &
                        x(l,j,n,1), x(l,j,n,2), x(l,j,n,3), &
                        x(l,m,n,1), x(l,m,n,2), x(l,m,n,3), &
                        x(i,m,n,1), x(i,m,n,2), x(i,m,n,3))

                   ! Set the volume to 1/6 of the sum of the volumes of the
                   ! pyramid. Remember that volpym computes 6 times the
                   ! volume.

                   vol(i,j,k) = sixth*(vp1 + vp2 + vp3 + vp4 + vp5 + vp6)

                   ! Check the volume and update the number of positive
                   ! and negative volumes if needed.

                   if( checkAll ) then

                      ! Update either the number of negative or positive
                      ! volumes. Negative volumes should only occur for left
                      ! handed blocks. This is checked later.
                      ! Set the logical volumeIsNeg accordingly.

                      if(vol(i,j,k) < zero) then
                         nVolNeg            = nVolNeg + 1
                         volumeIsNeg(i,j,k) = .true.
                      else
                         nVolPos            = nVolPos + 1
                         volumeIsNeg(i,j,k) = .false.
                      endif

                      ! Set the threshold for the volume quality.

                      fact = thresVolume*abs(vol(i,j,k))

                      ! Check the quality of the volume.

                      badVolume = .false.
                      if(vp1*vol(i,j,k) < zero .and. &
                           abs(vp1)       > fact) badVolume = .true.
                      if(vp2*vol(i,j,k) < zero .and. &
                           abs(vp2)       > fact) badVolume = .true.
                      if(vp3*vol(i,j,k) < zero .and. &
                           abs(vp3)       > fact) badVolume = .true.
                      if(vp4*vol(i,j,k) < zero .and. &
                           abs(vp4)       > fact) badVolume = .true.
                      if(vp5*vol(i,j,k) < zero .and. &
                           abs(vp5)       > fact) badVolume = .true.
                      if(vp6*vol(i,j,k) < zero .and. &
                           abs(vp6)       > fact) badVolume = .true.

                      ! Update nVolBad if this is a bad volume.

                      if( badVolume ) nVolBad = nVolBad + 1

                   endif

                   ! Set the volume to the absolute value.

                   vol(i,j,k) = abs(vol(i,j,k))

                enddo
             enddo
          enddo

          ! Some additional safety stuff for halo volumes.

          do k=2,kl
             do j=2,jl
                if(vol(1, j,k)/vol(2, j, k) < haloCellRatio) then
                   vol(1, j,k) = vol(2, j,k)
                end if
                if(vol(ie,j,k)/vol(il,j,k)  < haloCellRatio) then
                   vol(ie,j,k) = vol(il,j,k)
                end if
             enddo
          enddo

          do k=2,kl
             do i=1,ie
                if(vol(i,1, k)/vol(i,2,k) < haloCellRatio) then
                   vol(i,1, k) = vol(i,2, k)
                end if
                if(vol(i,je,k)/voL(i,jl,k) < haloCellRatio) then
                   vol(i,je,k) = vol(i,jl,k)
                end if
             enddo
          enddo

          do j=1,je
             do i=1,ie
                if(vol(i,j,1)/vol(i,j,2)  < haloCellRatio) then
                   vol(i,j,1)  = vol(i,j,2)
                end if
                if(vol(i,j,ke)/vol(i,j,kl) < haloCellRatio) then
                   vol(i,j,ke) = vol(i,j,kl)
                end if
             enddo
          enddo

          ! Determine the orientation of the block. For the fine level
          ! this is based on the number of positive and negative
          ! volumes; on the coarse levels the corresponding fine level
          ! value is taken. If both positive and negative volumes are
          ! present it is assumed that the block was intended to be
          ! right handed. The code will terminate later on anyway.

          if(level == 1) then
             if(nVolPos == 0) then       ! Left handed block.
                flowDoms(nn,level,sps)%rightHanded = .false.
             else                        ! Right handed (or bad) block.
                flowDoms(nn,level,sps)%rightHanded = .true.
             endif
          else
             flowDoms(nn,level,sps)%rightHanded = &
                  flowDoms(nn,1,sps)%rightHanded
          endif

          ! Set the factor in the surface normals computation. For a
          ! left handed block this factor is negative, such that the
          ! normals still point in the direction of increasing index.
          ! The formulae used later on assume a right handed block
          ! and fact is used to correct this for a left handed block,
          ! as well as the scaling factor of 0.5

          if( flowDoms(nn,level,sps)%rightHanded ) then
             fact =  half
          else
             fact = -half
          endif

          ! Check if both positive and negative volumes occur. If so,
          ! the block is bad and the counter nBlockBad is updated.

          if(nVolNeg > 0 .and. nVolPos > 0) then
             checkVolDoms(nn,sps)%blockHasNegVol = .true.
             nBlockBad = nBlockBad + 1
          else
             checkVolDoms(nn,sps)%blockHasNegVol = .false.
          endif
          !
          !           Computation of the face normals in i-, j- and k-direction.
          !           Formula's are valid for a right handed block; for a left
          !           handed block the correct orientation is obtained via fact.
          !           The normals point in the direction of increasing index.
          !           The absolute value of fact is 0.5, because the cross
          !           product of the two diagonals is twice the normal vector.
          !           Note that also the normals of the first level halo cells
          !           are computed. These are needed for the viscous fluxes.
          !
          ! Projected areas of cell faces in the i direction.

          do k=1,ke
             n = k -1
             do j=1,je
                m = j -1
                do i=0,ie

                   ! Determine the two diagonal vectors of the face.

                   v1(1) = x(i,j,n,1) - x(i,m,k,1)
                   v1(2) = x(i,j,n,2) - x(i,m,k,2)
                   v1(3) = x(i,j,n,3) - x(i,m,k,3)

                   v2(1) = x(i,j,k,1) - x(i,m,n,1)
                   v2(2) = x(i,j,k,2) - x(i,m,n,2)
                   v2(3) = x(i,j,k,3) - x(i,m,n,3)

                   ! The face normal, which is the cross product of the two
                   ! diagonal vectors times fact; remember that fact is
                   ! either -0.5 or 0.5.

                   si(i,j,k,1) = fact*(v1(2)*v2(3) - v1(3)*v2(2))
                   si(i,j,k,2) = fact*(v1(3)*v2(1) - v1(1)*v2(3))
                   si(i,j,k,3) = fact*(v1(1)*v2(2) - v1(2)*v2(1))

                enddo
             enddo
          enddo

          ! Projected areas of cell faces in the j direction.

          do k=1,ke
             n = k -1
             do j=0,je
                do i=1,ie
                   l = i -1

                   ! Determine the two diagonal vectors of the face.

                   v1(1) = x(i,j,n,1) - x(l,j,k,1)
                   v1(2) = x(i,j,n,2) - x(l,j,k,2)
                   v1(3) = x(i,j,n,3) - x(l,j,k,3)

                   v2(1) = x(l,j,n,1) - x(i,j,k,1)
                   v2(2) = x(l,j,n,2) - x(i,j,k,2)
                   v2(3) = x(l,j,n,3) - x(i,j,k,3)

                   ! The face normal, which is the cross product of the two
                   ! diagonal vectors times fact; remember that fact is
                   ! either -0.5 or 0.5.

                   sj(i,j,k,1) = fact*(v1(2)*v2(3) - v1(3)*v2(2))
                   sj(i,j,k,2) = fact*(v1(3)*v2(1) - v1(1)*v2(3))
                   sj(i,j,k,3) = fact*(v1(1)*v2(2) - v1(2)*v2(1))

                enddo
             enddo
          enddo

          ! Projected areas of cell faces in the k direction.

          do k=0,ke
             do j=1,je
                m = j -1
                do i=1,ie
                   l = i -1

                   ! Determine the two diagonal vectors of the face.

                   v1(1) = x(i,j,k,1) - x(l,m,k,1)
                   v1(2) = x(i,j,k,2) - x(l,m,k,2)
                   v1(3) = x(i,j,k,3) - x(l,m,k,3)

                   v2(1) = x(l,j,k,1) - x(i,m,k,1)
                   v2(2) = x(l,j,k,2) - x(i,m,k,2)
                   v2(3) = x(l,j,k,3) - x(i,m,k,3)

                   ! The face normal, which is the cross product of the two
                   ! diagonal vectors times fact; remember that fact is
                   ! either -0.5 or 0.5.

                   sk(i,j,k,1) = fact*(v1(2)*v2(3) - v1(3)*v2(2))
                   sk(i,j,k,2) = fact*(v1(3)*v2(1) - v1(1)*v2(3))
                   sk(i,j,k,3) = fact*(v1(1)*v2(2) - v1(2)*v2(1))

                enddo
             enddo
          enddo
          !
          !           The unit normals on the boundary faces. These always point
          !           out of the domain, so a multiplication by -1 is needed for
          !           the iMin, jMin and kMin boundaries.
          !
          ! Loop over the boundary subfaces of this block.

          bocoLoop: do mm=1,nBocos

             ! Determine the block face on which this subface is located
             ! and set ss and mult accordingly.

             select case (BCFaceID(mm))

             case (iMin)
                mult = -one; ss => si(1,:,:,:)

             case (iMax)
                mult = one;  ss => si(il,:,:,:)

             case (jMin)
                mult = -one; ss => sj(:,1,:,:)

             case (jMax)
                mult = one;  ss => sj(:,jl,:,:)

             case (kMin)
                mult = -one; ss => sk(:,:,1,:)

             case (kMax)
                mult = one;  ss => sk(:,:,kl,:)

             end select

             ! Loop over the boundary faces of the subface.

             do j=BCData(mm)%jcBeg, BCData(mm)%jcEnd
                do i=BCData(mm)%icBeg, BCData(mm)%icEnd

                   ! Compute the inverse of the length of the normal vector
                   ! and possibly correct for inward pointing.

                   xp = ss(i,j,1);  yp = ss(i,j,2);  zp = ss(i,j,3)
                   fact = sqrt(xp*xp + yp*yp + zp*zp)
                   if(fact > zero) fact = mult/fact

                   ! Compute the unit normal.

                   BCData(mm)%norm(i,j,1) = fact*xp
                   BCData(mm)%norm(i,j,2) = fact*yp
                   BCData(mm)%norm(i,j,3) = fact*zp

                enddo
             enddo

          enddo bocoLoop
          !
          !           Check in debug mode the sum of the normals of the cells.
          !           If everything is correct this should sum up to zero.
          !
          debugging: if( debug ) then

             ! Loop over the cells including the 1st level halo's.

             do k=2,kl
                n = k -1
                do j=2,jl
                   m = j -1
                   do i=2,il
                      l = i -1

                      ! Store the sum of the outward pointing surrounding
                      ! normals in v1. Due to the outward convention the
                      ! normals with the lowest index get a negative sign;
                      ! normals point in the direction of the higher index.

                      v1(1) = si(i,j,k,1) + sj(i,j,k,1) + sk(i,j,k,1) &
                           - si(l,j,k,1) - sj(i,m,k,1) - sk(i,j,n,1)
                      v1(2) = si(i,j,k,2) + sj(i,j,k,2) + sk(i,j,k,2) &
                           - si(l,j,k,2) - sj(i,m,k,2) - sk(i,j,n,2)
                      v1(3) = si(i,j,k,3) + sj(i,j,k,3) + sk(i,j,k,3) &
                           - si(l,j,k,3) - sj(i,m,k,3) - sk(i,j,n,3)

                      ! Store the inverse of the sum of the areas of the
                      ! six faces in fact.

                      fact = one/(sqrt(si(i,j,k,1)*si(i,j,k,1)  &
                           +           si(i,j,k,2)*si(i,j,k,2)  &
                           +           si(i,j,k,3)*si(i,j,k,3)) &
                           +      sqrt(si(l,j,k,1)*si(l,j,k,1)  &
                           +           si(l,j,k,2)*si(l,j,k,2)  &
                           +           si(l,j,k,3)*si(l,j,k,3)) &
                           +      sqrt(sj(i,j,k,1)*sj(i,j,k,1)  &
                           +           sj(i,j,k,2)*sj(i,j,k,2)  &
                           +           sj(i,j,k,3)*sj(i,j,k,3)) &
                           +      sqrt(sj(i,m,k,1)*sj(i,m,k,1)  &
                           +           sj(i,m,k,2)*sj(i,m,k,2)  &
                           +           sj(i,m,k,3)*sj(i,m,k,3)) &
                           +      sqrt(sk(i,j,k,1)*sk(i,j,k,1)  &
                           +           sk(i,j,k,2)*sk(i,j,k,2)  &
                           +           sk(i,j,k,3)*sk(i,j,k,3)) &
                           +      sqrt(sk(i,j,n,1)*sk(i,j,n,1)  &
                           +           sk(i,j,n,2)*sk(i,j,n,2)  &
                           +           sk(i,j,n,3)*sk(i,j,n,3)))

                      ! Multiply v1 by fact to obtain a nonDimensional
                      ! quantity and take tha absolute value of it.

                      v1(1) = abs(v1(1)*fact)
                      v1(2) = abs(v1(2)*fact)
                      v1(3) = abs(v1(3)*fact)

                      ! Check if the control volume is closed.

                      if(v1(1) > thresholdReal .or. &
                           v1(2) > thresholdReal .or. &
                           v1(3) > thresholdReal)     &
                           call terminate("metric", &
                           "Normals do not sum up to 0")

                   enddo
                enddo
             enddo

          endif debugging

       enddo domains
    enddo spectral

    ! Determine the global number of bad blocks. The result must be
    ! known on all processors and thus an allreduce is needed.

    call mpi_allreduce(nBlockBad, nBlockBadGlobal, 1, adflow_integer, &
         mpi_sum, ADflow_comm_world, ierr)

    ! Test if bad blocks are present in the grid. If so, the action
    ! taken depends on the grid level.

    if(nBlockBadGlobal > 0) then
       if(level == 1) then

          ! Negative volumes present on the fine grid level. Print a
          ! list of the bad volumes and terminate executation.

          call writeNegVolumes(checkVolDoms)

          call returnFail("metric", &
            "Negative volumes present in grid.")
          call mpi_barrier(ADflow_comm_world, ierr)

       else

          ! Coarser grid level. The fine grid is okay, but due to the
          ! coarsening negative volumes are introduced. Print a warning.

          if(myID == 0) then
             print "(a)", "#"
             print "(a)", "#                      Warning"
             print 100, level
100          format("#* Negative volumes present on coarse grid &
                  &level",1x,i1,".")
             print "(a)", "#* Computation continues, &
                  &but be aware of this"
             print "(a)", "#"
          endif

       endif
    endif

    ! Determine the global number of bad volumes. The result will
    ! only be known on processor 0. The quality volume check will
    ! only be done for the finest grid level.

    if(level == 1) then
       call mpi_reduce(nVolBad, nVolBadGlobal, 1, adflow_integer, &
            mpi_sum, 0, ADflow_comm_world, ierr)

       ! Print a warning in case bad volumes were found. Only processor
       ! 0 prints this warning.

       if(myID == 0 .and. nVolBadGlobal > 0 .and. printWarnings) then
          write(integerString,"(i10)") nVolBadGlobal
          integerString = adjustl(integerString)
          integerString = trim(integerString)
          print "(a)", "#"
          print "(a)", "#                      Warning"
          print 101, trim(integerString)
          print 102
          print "(a)", "#"
101       format("# ",a," bad quality volumes found.")
102       format("# Computation will continue, but be aware of this")
       endif
    endif

    ! Release the memory of volumeIsNeg of all local blocks again.

    do sps=1,nTimeIntervalsSpectral
       do nn=1,nDom
          deallocate(checkVolDoms(nn,sps)%volumeIsNeg, stat=ierr)
          if(ierr /= 0)              &
               call terminate("metric", &
               "Deallocation failure for volumeIsNeg")
       enddo
    enddo

  contains

    !        ================================================================

    function volpym(xa,ya,za,xb,yb,zb,xc,yc,zc,xd,yd,zd)
      !
      !         volpym computes 6 times the volume of a pyramid. Node p,
      !         whose coordinates are set in the subroutine metric itself,
      !         is the top node and a-b-c-d is the quadrilateral surface.
      !         It is assumed that the cross product vCa * vDb points in
      !         the direction of the top node. Here vCa is the diagonal
      !         running from node c to node a and vDb the diagonal from
      !         node d to node b.
      !
      use precision
      implicit none
      !
      !        Function type.
      !
      real(kind=realType) :: volpym
      !
      !        Function arguments.
      !
      real(kind=realType), intent(in) :: xa, ya, za, xb, yb, zb
      real(kind=realType), intent(in) :: xc, yc, zc, xd, yd, zd

      volpym = (xp - fourth*(xa + xb  + xc + xd))              &
           * ((ya - yc)*(zb - zd) - (za - zc)*(yb - yd))   + &
           (yp - fourth*(ya + yb  + yc + yd))              &
           * ((za - zc)*(xb - xd) - (xa - xc)*(zb - zd))   + &
           (zp - fourth*(za + zb  + zc + zd))              &
           * ((xa - xc)*(yb - yd) - (ya - yc)*(xb - xd))

    end function volpym

  end subroutine metric

  !      ==================================================================

  subroutine writeNegVolumes(checkVolDoms)
    !
    !       writeNegVolumes writes the negative volumes of a block to
    !       stdout. If a block is flagged to have negative volumes it is
    !       assumed that the block is intended to be a right handed block.
    !
    use constants
    use blockPointers
    use cgnsGrid
    use communication
    use inputPhysics
    use inputTimeSpectral
    use checkVolBlock
    use utils, only : setPointers, terminate
    implicit none
    !
    !      Subroutine arguments.
    !
    type(checkVolBlockType), &
         dimension(nDom,nTimeIntervalsSpectral), intent(in) :: checkVolDoms
    !
    !      Local variables.
    !
    integer :: proc, ierr

    integer(kind=intType) :: nn, sps, i, j, k

    real(kind=realType), dimension(3) :: xc

    character(len=10) :: intString1, intString2, intString3

    ! Processor 0 prints a message that negative volumes are present
    ! in the grid.

    if(myID == 0) then
       print "(a)", "#"
       print "(a)", "#                      Error"
       print "(a)", "# Negative volumes found in the grid."
       print "(a)", "# A list of the negative volumes is &
            &printed below"
       print "(a)", "#"
    endif

    ! Loop over the processors such that a clean output is obtained.
    ! This may not be the most efficient solution, but that is not
    ! an issue here.

    procLoop: do proc=0,(nProc-1)

       ! Test if I'm to one that must write my bad volumes.

       testIWrite: if(proc == myID) then

          ! Loop over the number of spectral solutions and local blocks.

          spectral: do sps=1,nTimeIntervalsSpectral
             domains: do nn=1,nDom

                ! Test for a bad block.

                testBad: if( checkVolDoms(nn,sps)%blockHasNegVol ) then

                   ! Set the pointers for this block.

                   call setPointers(nn, 1_intType, sps)

                   ! Write the name of the block. The error message
                   ! depends a bit on the case computed. For a time
                   ! spectral solution also the spectral solution is
                   ! mentioned, for steady and unsteady this is not
                   ! the case, because there is only one.

                   select case (equationMode)
                   case (steady, unsteady)

                      print "(a)", "#"
                      print 100, trim(cgnsDoms(nbkGlobal)%zoneName)
100                   format("# Block",1x,a,1x,"contains the following &
                           &negative volumes")
                      print "(a)", "#=================================&
                           &==================================="
                      print "(a)", "#"

                      !====================================================

                   case (timeSpectral)

                      write(intString1,"(i10)") sps
                      intString1 = adjustl(intString1)

                      print "(a)", "#"
                      print 101, trim(intString1), &
                           trim(cgnsDoms(nbkGlobal)%zoneName)
101                   format("# Spectral solution",1x,a, ":block", &
                           1x,a,1x,"contains the following negative &
                           &volumes")
                      print "(a)", "#=================================&
                           &==================================="
                      print "(a)", "#"

                   end select

                   ! Loop over the owned volumes and write the
                   ! negative ones.

                   do k=2,kl
                      do j=2,jl
                         do i=2,il
                            if(checkVolDoms(nn,sps)%volumeIsNeg(i,j,k)) then

                               xc(1:3) = eighth*(x(i-1,j-1,k-1,1:3) &
                                    +         x(i,  j-1,k-1,1:3) &
                                    +         x(i-1,j  ,k-1,1:3) &
                                    +         x(i,  j  ,k-1,1:3) &
                                    +         x(i-1,j-1,k  ,1:3) &
                                    +         x(i,  j-1,k  ,1:3) &
                                    +         x(i-1,j  ,k  ,1:3) &
                                    +         x(i,  j  ,k  ,1:3))

                               write(intString1,"(i10)") i
                               write(intString2,"(i10)") j
                               write(intString3,"(i10)") k

                               intString1 = adjustl(intString1)
                               intString2 = adjustl(intString2)
                               intString3 = adjustl(intString3)

                               print 102, trim(intString1), &
                                    trim(intString2), &
                                    trim(intString3), &
                                    xc(1), xc(2), xc(3), -vol(i,j,k)
102                            format("# Indices (",a,",",a,",",a,"), &
                                    &coordinates (",e10.3,",",e10.3,",", &
                                    e10.3,"), Volume: ",e10.3)

                            endif
                         enddo
                      enddo
                   enddo


                endif testBad
             enddo domains
          enddo spectral

       endif testIWrite

       ! Synchronize the processors to avoid a messy output.

       call mpi_barrier(ADflow_comm_world, ierr)

    enddo procLoop

  end subroutine writeNegVolumes
  subroutine faceRotationMatrices(level, allocMem)
    !
    !       faceRotationMatrices computes the rotation matrices on the
    !       faces, such that for a rotationally periodic the nonlinear
    !       reconstruction in the upwind schemes is consistent with its
    !       periodic neighbor.
    !
    use constants
    use blockPointers
    use inputDiscretization
    use inputTimeSpectral
    use section
    use utils, only : setPointers, terminate
    implicit none
    !
    !      Subroutine arguments.
    !
    integer(kind=intType), intent(in) :: level
    logical,               intent(in) :: allocMem
    !
    !      Local variables.
    !
    integer :: ierr

    integer(kind=intType) :: nn, sps, mm

    real(kind=realType), dimension(:,:,:),   pointer :: xFace
    real(kind=realType), dimension(:,:,:,:), pointer :: rotFace

    real(kind=realType), dimension(3) :: axis, vecR1, vecR2, rotCenter

    ! If this is not the finest level, return. The matrices are
    ! only needed on the finest grid.

    if(level /= 1) return

    ! Check if an upwind scheme is used. If not, return.

    if(spaceDiscr /= upwind) return

    ! Check if a limiter is used at all. The rotation matrices are
    ! only needed when a nonlinear construction is used.

    if(limiter == firstOrder .or. limiter == noLimiter) return

    ! Check if rotational periodicity occurs at all. If not, there
    ! is no need for the rotation matrices either.

    do nn=1,nSections
       if(sections(nn)%nSlices > 1) exit
    enddo
    if(nn > nSections) return

    ! Loop over the number of blocks..

    domains: do nn=1,nDom

       ! Hard coded section ID to 1
       sectionID = 1

       ! Determine the two unit vectors in the plane normal to
       ! the rotation axis of this section.

       axis      = sections(sectionID)%rotAxis
       rotCenter = sections(sectionID)%rotCenter
       call unitVectorsInAxialPlane(axis, vecR1, vecR2)

       ! Loop over the number of time instances.

       spectral: do sps=1,nTimeIntervalsSpectral

          ! Check if the memory for the rotation matrices must be
          ! allocated and do so if needed.

          if( allocMem ) then

             il = flowDoms(nn,1,sps)%il
             jl = flowDoms(nn,1,sps)%jl
             kl = flowDoms(nn,1,sps)%kl

             allocate(flowDoms(nn,1,sps)%rotMatrixI(il,2:jl,2:kl,3,3), &
                  flowDoms(nn,1,sps)%rotMatrixJ(2:il,jl,2:kl,3,3), &
                  flowDoms(nn,1,sps)%rotMatrixK(2:il,2:jl,kl,3,3), &
                  stat=ierr)
             if(ierr /= 0)                            &
                  call terminate("faceRotationMatrices", &
                  "Memory allocation failure for the &
                  &rotation matrices.")
          endif

          ! Set the pointers to this block.

          call setPointers(nn,level,sps)

          ! The rotation matrices for the i-faces.

          do mm=1,il
             xFace   => x(mm,1:,1:,:);
             rotFace => rotMatrixI(mm,:,:,:,:)

             call computeRotMatrixFace(xFace, rotFace, jl, kl)
          enddo

          ! The rotation matrices for the j-faces.

          do mm=1,jl
             xFace   => x(1:,mm,1:,:);
             rotFace => rotMatrixJ(:,mm,:,:,:)

             call computeRotMatrixFace(xFace, rotFace, il, kl)
          enddo

          ! The rotation matrices for the k-faces.

          do mm=1,kl
             xFace   => x(1:,1:,mm,:);
             rotFace => rotMatrixK(:,:,mm,:,:)

             call computeRotMatrixFace(xFace, rotFace, il, jl)
          enddo

       enddo spectral
    enddo domains

    !=================================================================

  contains

    !===============================================================

    subroutine computeRotMatrixFace(xx, rotMat, iil, jjl)
      !
      !         computeRotMatrixFace is an internal subroutine, which
      !         computes the rotation matrix from Cartesian to local
      !         cylindrical velocity components for the face centers.
      !
      implicit none
      !
      !        Subroutine arguments.
      !
      integer(kind=intType), intent(in) :: iil, jjl

      real(kind=realType), dimension(:,:,:),     intent(in)  :: xx
      real(kind=realType), dimension(2:,2:,:,:), intent(out) :: rotMat
      !
      !        Local variables.
      !
      integer(kind=intType) :: i, j

      real(kind=realType) :: r1, r2, rInv, cosTheta, sinTheta

      real(kind=realType), dimension(3) :: xF

      ! Loop over the face centers.

      do j=2,jjl
         do i=2,iil

            ! Compute the coordinates of the face center relative to
            ! the center of rotation.

            xF(1) = fourth*(xx(i-1,j-1,1) + xx(i-1,j,1) &
                 +         xx(i,  j-1,1) + xx(i,  j,1)) - rotCenter(1)
            xF(2) = fourth*(xx(i-1,j-1,2) + xx(i-1,j,2) &
                 +         xx(i,  j-1,2) + xx(i,  j,2)) - rotCenter(2)
            xF(3) = fourth*(xx(i-1,j-1,3) + xx(i-1,j,3) &
                 +         xx(i,  j-1,3) + xx(i,  j,3)) - rotCenter(3)

            ! Determine the two radial components for this point.

            r1 = xF(1)*vecR1(1) + xF(2)*vecR1(2) + xF(3)*vecR1(3)
            r2 = xF(1)*vecR2(1) + xF(2)*vecR2(2) + xF(3)*vecR2(3)

            ! Determine the sine and cosine of the polar angle.

            rInv     = one/sqrt(r1*r1 + r2*r2)
            cosTheta = r1*rInv
            sinTheta = r2*rInv

            ! Compute the transformation matrix.

            rotMat(i,j,1,1) = axis(1)
            rotMat(i,j,1,2) = axis(2)
            rotMat(i,j,1,3) = axis(3)

            rotMat(i,j,2,1) = cosTheta*vecR1(1) + sinTheta*vecR2(1)
            rotMat(i,j,2,2) = cosTheta*vecR1(2) + sinTheta*vecR2(2)
            rotMat(i,j,2,3) = cosTheta*vecR1(3) + sinTheta*vecR2(3)

            rotMat(i,j,3,1) = cosTheta*vecR2(1) - sinTheta*vecR1(1)
            rotMat(i,j,3,2) = cosTheta*vecR2(2) - sinTheta*vecR1(2)
            rotMat(i,j,3,3) = cosTheta*vecR2(3) - sinTheta*vecR1(3)

         enddo
      enddo

    end subroutine computeRotMatrixFace

  end subroutine faceRotationMatrices
  subroutine updateCoordinatesAllLevels
    !
    !       updateCoordinatesAllLevels updates the coordinates of all
    !       grid levels, assuming that the owned coordinates of the fine
    !       grid are known.
    !
    use constants
    use block
    use iteration
    use inputTimeSpectral
    use blockPointers
    use coarseUtils, only :coarseOwnedCoordinates
    implicit none
    !
    !      Local variables.
    !
    integer(kind=intType) :: nLevels, nn
    real(kind=realType)   :: origGroundLevel

    ! Determine the halo coordinates of the fine level.
    origGroundLevel = groundLevel
    groundLevel     = 1
    call xhalo(groundLevel)

    ! Loop over the coarse grid levels; first the owned coordinates
    ! are determined, followed by the halo's.

    nLevels = ubound(flowDoms,2)
    do nn=(groundLevel+1),nLevels
       call coarseOwnedCoordinates(nn)
       call xhalo(nn)
    enddo

    groundLevel = origGroundLevel

  end subroutine updateCoordinatesAllLevels

  !      ==================================================================

  subroutine updateMetricsAllLevels
    !
    !       updateMetricsAllLevels recomputes the metrics on all grid
    !       levels. This routine is typically called when the coordinates
    !       have changed, but the connectivity remains the same, i.e. for
    !       moving or deforming mesh problems.
    !
    use constants
    use block
    use iteration
    use inputphysics
    use inputIteration
    implicit none
    !
    !      Local variables.
    !
    integer(kind=intType) :: nLevels, nn

    ! Loop over the grid levels and call metric and checkSymmetry.

    nLevels = ubound(flowDoms,2)
    do nn=groundLevel,nLevels
       if(equationMode == unsteady) then
          call metric(nn)
       else
          call metric(nn)
       end if

       if (printWarnings) then
          call checkSymmetry(nn)
       end if
    enddo

  end subroutine updateMetricsAllLevels

  subroutine updateGridVelocitiesAllLevels

    !
    !       updateGridVelocitesAllLevels recomputes the rotational
    !       parameters on all grid
    !       levels. This routine is typically called when the coordinates
    !       have changed, but the connectivity remains the same, i.e. for
    !       moving or deforming mesh problems.
    !
    use constants
    use block
    use iteration
    use section
    use monitor
    use inputTimeSpectral
    use inputPhysics
    use solverUtils
    implicit none

    !subroutine variables

    !Local Variables

    integer(kind=inttype):: mm,nnn

    real(kind=realType), dimension(nSections) :: t

    do mm=1,nTimeIntervalsSpectral

       ! Compute the time, which corresponds to this spectral solution.
       ! For steady and unsteady mode this is simply the restart time;
       ! for the spectral mode the periodic time must be taken into
       ! account, which can be different for every section.

       t = timeUnsteadyRestart

       if(equationMode == timeSpectral) then
          do nnn=1,nSections
             t(nnn) = t(nnn) + (mm-1)*sections(nnn)%timePeriod &
                  /         real(nTimeIntervalsSpectral,realType)
          enddo
       endif

       call gridVelocitiesFineLevel(.false., t, mm)
       call gridVelocitiesCoarseLevels(mm)
       call normalVelocitiesAllLevels(mm)

       call slipVelocitiesFineLevel(.false., t, mm)
       call slipVelocitiesCoarseLevels(mm)

    enddo

  end subroutine updateGridVelocitiesAllLevels

  subroutine updatePeriodicInfoAllLevels


    !
    !       updatePeriodicInfoAllLevels recomputes the spectral parameters
    !       on all grid levels. This routine is typically called when the
    !       frequnecy or amplitude of the oscillation in the time spectral
    !       computation has changed
    !
    use block
    use iteration
    use section
    use monitor
    use inputTimeSpectral
    use inputPhysics
    use communication
    use initializeFlow, onlY : timeSpectralMatrices
    use partitioning, only : fineGridSpectralCoor, timeRotMatricesSpectral, &
         timePeriodSpectral
    !
    implicit none

    ! Determine for the time spectral mode the time of one period,
    ! the rotation matrices for the velocity components and
    ! create the fine grid coordinates of all time spectral locations.

    call timePeriodSpectral
    call timeRotMatricesSpectral
    call fineGridSpectralCoor
    call timeSpectralMatrices


  end subroutine updatePeriodicInfoAllLevels
  subroutine unitVectorsInAxialPlane(axis, vecR1, vecR2)
    !
    !       unitVectorsInAxialPlane computes from the given unit vector
    !       axis the two unit vectors which describe the plane normal to
    !       axis. There is of course an ambiguity in this choice, but this
    !       is not a problem as long as the choice is consistent
    !       throughout the code.
    !
    use constants
    implicit none
    !
    !      Subroutine arguments.
    !
    real(kind=realType), dimension(3), intent(in)  :: axis
    real(kind=realType), dimension(3), intent(out) :: vecR1, vecR2
    !
    !      Local variables.
    !
    real(kind=realType) :: dot

    ! The vectors which span the axial plane must be normal to axis.
    ! For the first vector try first the y-axis. If not good enough
    ! use the z-axis.

    if(abs(axis(2)) < 0.707107_realType) then
       vecR1(1) = zero
       vecR1(2) = one
       vecR1(3) = zero
    else
       vecR1(1) = zero
       vecR1(2) = zero
       vecR1(3) = one
    endif

    ! Make sure that vecR1 is normal to axis. Create a unit
    ! vector again.

    dot = vecR1(1)*axis(1) + vecR1(2)*axis(2) + vecR1(3)*axis(3)
    vecR1(1) = vecR1(1) - dot*axis(1)
    vecR1(2) = vecR1(2) - dot*axis(2)
    vecR1(3) = vecR1(3) - dot*axis(3)

    dot = one/sqrt(vecR1(1)**2 + vecR1(2)**2 + vecR1(3)**2)
    vecR1(1) = vecR1(1)*dot
    vecR1(2) = vecR1(2)*dot
    vecR1(3) = vecR1(3)*dot

    ! Create the second vector which spans the axial plane. This must
    ! be normal to both axis and vecR1, i.e. the cross-product.

    vecR2(1) = axis(2)*vecR1(3) - axis(3)*vecR1(2)
    vecR2(2) = axis(3)*vecR1(1) - axis(1)*vecR1(3)
    vecR2(3) = axis(1)*vecR1(2) - axis(2)*vecR1(1)

  end subroutine unitVectorsInAxialPlane

  subroutine preprocessingADjoint
    !
    !      Perform the preprocessing tasks for the adjoint solver. This
    !      routine is called only once. The memory allcoated here is
    !      deallocated in src/utils/releaseMemory.f90
    !
    use constants
    use communication, only : adflow_comm_world
    use adjointVars, only :nCellsLocal, nNOdesLocal
    use flowVarRefState, only : nw, nwf
    use inputTimeSpectral, only : nTimeIntervalsSpectral
    use inputAdjoint, only : frozenTurbulence
    use ADjointPETSc, only: w_like1, w_like2, PETScIerr, &
         psi_like1, psi_like2, x_like, psi_like3
    use utils, only : setPointers, EChk
#include <petsc/finclude/petsc.h>
    use petsc
    implicit none

    !     Local variables.
    !
    integer(kind=intType) :: ndimW, ncell, nState, nDimPsi, nDimX
    !
    !
    ! Create PETSc Vectors that are actually empty. These do NOT take
    ! any (substantial) memory. We want to keep these around inbetween
    ! creations/deletions of adjoint/NKsolver memory

    ! Setup number of state variable based on turbulence assumption
    if ( frozenTurbulence ) then
       nState = nwf
    else
       nState = nw
    endif

    nDimW = nw * nCellsLocal(1_intType)*nTimeIntervalsSpectral
    nDimPsi = nState*  nCellsLocal(1_intType)*nTimeIntervalsSpectral
    nDimX = 3 * nNodesLocal(1_intType)*nTimeIntervalsSpectral

    ! Two w-like vectors.
    call VecCreateMPIWithArray(ADFLOW_COMM_WORLD,nw,ndimW,PETSC_DECIDE, &
         PETSC_NULL_SCALAR,w_like1,PETScIerr)
    call EChk(PETScIerr,__FILE__,__LINE__)

    call VecCreateMPIWithArray(ADFLOW_COMM_WORLD,nw,ndimW,PETSC_DECIDE, &
         PETSC_NULL_SCALAR,w_like2,PETScIerr)
    call EChk(PETScIerr,__FILE__,__LINE__)

    ! Two psi-like vectors.
    call VecCreateMPIWithArray(ADFLOW_COMM_WORLD,nState,ndimPsi,PETSC_DECIDE, &
         PETSC_NULL_SCALAR,psi_like1,PETScIerr)
    call EChk(PETScIerr,__FILE__,__LINE__)

    call VecCreateMPIWithArray(ADFLOW_COMM_WORLD,nstate,ndimPsi,PETSC_DECIDE, &
         PETSC_NULL_SCALAR,psi_like2,PETScIerr)
    call EChk(PETScIerr,__FILE__,__LINE__)

    call VecCreateMPIWithArray(ADFLOW_COMM_WORLD,nstate,ndimPsi,PETSC_DECIDE, &
         PETSC_NULL_SCALAR,psi_like3,PETScIerr)
    call EChk(PETScIerr,__FILE__,__LINE__)

    call VecCreateMPIWithArray(ADFLOW_COMM_WORLD,3,ndimX,PETSC_DECIDE, &
         PETSC_NULL_SCALAR,x_like,PETScIerr)
    call EChk(PETScIerr,__FILE__,__LINE__)

    ! Need to initialize the stencils as well, only once:
    call initialize_stencils

  end subroutine preprocessingADjoint

  subroutine updateReferencePoint
    !
    !       reruns the initialization routines to update AOA and other
    !       flow variables after a design change
    !
    use constants
    use inputTimeSpectral, only : nTimeINTervalsSpectral
    use section, only : sections, nSections
    use inputPhysics, only : equationMode
    use inputMotion, only : rotPoint
    use cgnsGrid, only : cgnsDoms
    use monitor, only : timeUnsteadyRestart
    use iteration, only : groundLevel
    use blockpointers, only : nDom, nbkGlobal
    use utils, onlY : setPointers
    use solverUtils
    implicit none

    ! Working variables
    integer(kind=intType) ::mm,nnn,nn
    real(kind=realType), dimension(nSections) :: t

    groundlevel = 1
    do mm=1,nTimeIntervalsSpectral
       do nn=1,nDom
          ! Set the pointers for this block.
          call setPointers(nn, groundLevel, mm)
          !lref is outside
          cgnsDoms(nbkglobal)%rotCenter = rotPoint
       enddo
    enddo

    groundlevel = 1
    do mm=1,nTimeIntervalsSpectral

       ! Compute the time, which corresponds to this spectral solution.
       ! For steady and unsteady mode this is simply the restart time;
       ! for the spectral mode the periodic time must be taken into
       ! account, which can be different for every section.
       t = timeUnsteadyRestart

       if(equationMode == timeSpectral) then
          do nnn=1,nSections
             t(nnn) = t(nnn) + (mm-1)*sections(nnn)%timePeriod &
                  /         real(nTimeIntervalsSpectral,realType)
          enddo
       endif

       call gridVelocitiesFineLevel(.false., t, mm)
       call gridVelocitiesCoarseLevels(mm)
       call normalVelocitiesAllLevels(mm)
       call slipVelocitiesFineLevel(.false., t, mm)
       call slipVelocitiesCoarseLevels(mm)
    enddo
  end subroutine updateReferencePoint

  subroutine updateRotationRate(rotCenter, rotRate, blocks, nblocks)

    use constants
    use inputTimeSpectral, only : nTimeInTervalsSpectral
    use section, only : sections, nSections
    use inputPhysics, only : equationMode
    use inputMotion, only : rotPoint
    use cgnsGrid, only : cgnsDoms
    use monitor, only : timeUnsteadyRestart
    use iteration, only : groundLevel
    use blockpointers, only : nDom, nbkGlobal, flowDoms
    use solverUtils
    implicit none

    real(kind=realType),intent(in)::rotCenter(3), rotRate(3)
    integer(kind=intType), intent(in) :: nblocks
    integer(kind=intType), intent(in) :: blocks(nblocks)

    integer(kind=intType) ::mm,nnn,nn, level, sps, i
    real(kind=realType), dimension(nSections) :: t

    groundlevel = 1

    do nn=1,nblocks
       cgnsDoms(nn)%rotRate = rotRate
       cgnsDoms(nn)%rotCenter = rotCenter

       do i=1,cgnsDoms(nn)%nBocos
          cgnsDoms(nn)%bocoInfo(i)%rotRate = rotRate
       end do
    enddo

    do sps=1,nTimeIntervalsSpectral
       do level=1,ubound(flowDoms,2)
          do nn=1,nDom
             flowDoms(nn,level,sps)%blockIsMoving = .True.
             flowDoms(nn,level,sps)%addGridVelocities = .True.
          end do
       end do
    end do
    groundlevel = 1
    do mm=1,nTimeIntervalsSpectral

       ! Compute the time, which corresponds to this spectral solution.
       ! For steady and unsteady mode this is simply the restart time;
       ! for the spectral mode the periodic time must be taken into
       ! account, which can be different for every section.

       t = timeUnsteadyRestart

       if(equationMode == timeSpectral) then
          do nnn=1,nSections
             t(nnn) = t(nnn) + (mm-1)*sections(nnn)%timePeriod &
                  /         real(nTimeIntervalsSpectral,realType)
          enddo
       endif

       call gridVelocitiesFineLevel(.false., t, mm)
       call gridVelocitiesCoarseLevels(mm)
       call normalVelocitiesAllLevels(mm)
       call slipVelocitiesFineLevel(.false., t, mm)
       call slipVelocitiesCoarseLevels(mm)

    enddo

  end subroutine updateRotationRate

end module preprocessingAPI

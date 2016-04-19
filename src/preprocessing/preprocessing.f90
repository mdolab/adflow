!
!      ******************************************************************
!      *                                                                *
!      * File:          preprocessing.f90                               *
!      * Author:        Edwin van der Weide, Steve Repsher              *
!      * Starting date: 01-20-2003                                      *
!      * Last modified: 11-30-2007                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine preprocessing
!
!      ******************************************************************
!      *                                                                *
!      * preprocessing determines the communication patterns between    *
!      * the processors for all the mg levels, computes the wall        *
!      * distances and the metrics.                                     *
!      *                                                                *
!      ******************************************************************
!
       use block
       use blockPointers
       use cgnsGrid
       use commMixing
       use commSliding
       use communication
       use inputPhysics
       use inputTimeSpectral
       use interfaceGroups
       use section
       use wallDistanceData
       use overset
       implicit none
!
!      Local variables.
!
       integer :: ierr

       integer(kind=intType) :: nLevels, level, nn, mm, nsMin, nsMax, i, iProc
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
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
           call returnFail("preprocessing", &
                          "Different rotational periodicity encountered &
                          &for time accurate computation")
         call mpi_barrier(SUmb_comm_world, ierr)

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
         call returnFail("preprocessing", &
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

       ! Allocate the memory for the sliding mesh communication pattern.
       ! This pattern changes in time and therefore each spectral time
       ! value has its own sliding mesh communication pattern.

       mm = nTimeIntervalsSpectral

       allocate(commSlidingCell_1st(nn,mm), &
                commSlidingCell_2nd(nn,mm), &
                intSlidingCell_1st(nn,mm),  &
                intSlidingCell_2nd(nn,mm),  stat=ierr)
       if(ierr /= 0)                     &
         call returnFail("preprocessing", &
                        "Memory allocation failure for &
                        &slidingCommPatterns")

       ! Allocate the memory for the communication of the mixing plane
       ! halo cells. As this type of "boundary condition" can only
       ! be applied for steady state computations, there is no
       ! dependency on the number of time spectral solutions.

       mm = nInterfaceGroups
       allocate(commPatternMixing(nn,mm,2), stat=ierr)
       if(ierr /= 0)                     &
         call returnFail("preprocessing", &
                        "Memory allocation failure for &
                        &commPatternMixing")

       ! Allocate the memory for the overset mesh communication pattern.
       ! This pattern changes in time and therefore each spectral time
       ! value has its own sliding mesh communication pattern.

       mm = nTimeIntervalsSpectral
       allocate(commPatternOverset(nn,mm), internalOverset(nn,mm), &
            overlapMatrix(nn, mm), stat=ierr)
       if(ierr /= 0)                     &
         call returnFail("preprocessing", &
                        "Memory allocation failure for commOverset")

       ! Determine the fine grid 1 to 1 matching communication pattern.

       call determineCommPattern(1_intType)

       ! Initialize the send and receive buffer sizes to 0 and determine
       ! its size for the finest grid for the 1 to 1 communication.

       sendBufferSize_1to1  = 0
       recvBufferSize_1to1  = 0
       sendBufferSizeSlide  = 0
       recvBufferSizeSlide  = 0
       sendBufferSizeOver   = 0
       recvBufferSizeOver   = 0

       call setBufferSizes(1_intType, 1_intType, .true., .false., .false.)

       ! Loop to create the coarse grid levels.

       do level=2,nLevels

         ! Create the coarse grid blocks, its communication pattern, the
         ! coarse grid level 0 cooling parameters and check the
         ! communication buffer sizes.

         call createCoarseBlocks(level)
         call determineCommPattern(level)
         call coarseLevel0CoolingParameters(level)
         call setBufferSizes(level, 1_intType, .true., .false., .false.)

       enddo

       ! Synchronize the processors, just to be sure.

       call mpi_barrier(SUmb_comm_world, ierr)

       ! Allocate memory for the nonblocking point to point communication.

       allocate(sendBuffer(sendBufferSize), &
                recvBuffer(recvBufferSize), stat=ierr)
       if(ierr /= 0)                     &
         call returnFail("preprocessing", &
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
       call mpi_allreduce(nDom, nDomTotal, 1, sumb_integer, MPI_SUM, &
            sumb_comm_world, ierr)
       
       ! Receive the number of domains from each proc using an allgather.
       call mpi_allgather(nDom, 1, sumb_integer, nDomProc, 1, sumb_integer, &
            sumb_comm_world, ierr)
       
       ! Compute the cumulative format:
       cumDomProc(0) = 0
       do iProc=1, nProc
          cumDomProc(iProc) = cumDomProc(iProc-1) + nDomProc(iProc-1)
       end do

       ! Determine the number of grid clusters
       call determineClusters()


       ! Loop over the number of levels and perform a lot of tasks.
       ! See the corresponding subroutine header, although the
       ! names are pretty self-explaining


       do level=1,nLevels
         call xhalo(level)
         call slidingComm(level, .true.)
         call allocateMetric(level)
         call metric(level)
         call setPorosities(level)
         call setFamilyInfoFaces(level)
         call faceRotationMatrices(level, .true.)
         call checkSymmetry(level)
         call viscSubfaceInfo(level)
         call determineAreaLevel0Cooling(level)
         call determineNcellGlobal(level)
         call setGlobalCellsAndNodes(level)
         call wallDistance(level, .True.)
      end do

      ! BC Data must be alloaced (for surface iblank) before we can do
      ! the overset computation.
      call allocMemBCData

      do level=1,nLevels
         if (level == 1) then
            call oversetComm(level, .true., .false.)

         else
            call oversetComm(level, .true., .true.)
         end if
      end do

       
     end subroutine preprocessing


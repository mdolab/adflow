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
       implicit none
!
!      Local variables.
!
       integer :: ierr

       integer(kind=intType) :: nLevels, level, nn, mm, nsMin, nsMax
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
           call terminate("preprocessing", &
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
         call terminate("preprocessing", &
                        "Memory allocation failure for commPatterns")

       ! Allocate the memory for the sliding mesh communication pattern.
       ! This pattern changes in time and therefore each spectral time
       ! value has its own sliding mesh communication pattern.

       mm = nTimeIntervalsSpectral

       allocate(commSlidingCell_1st(nn,mm), &
                commSlidingCell_2nd(nn,mm), &
                intSlidingCell_1st(nn,mm),  &
                intSlidingCell_2nd(nn,mm),  stat=ierr)
       if(ierr /= 0)                     &
         call terminate("preprocessing", &
                        "Memory allocation failure for &
                        &slidingCommPatterns")

       ! Allocate the memory for the communication of the mixing plane
       ! halo cells. As this type of "boundary condition" can only
       ! be applied for steady state computations, there is no
       ! dependency on the number of time spectral solutions.

       mm = nInterfaceGroups
       allocate(commPatternMixing(nn,mm,2), stat=ierr)
       if(ierr /= 0)                     &
         call terminate("preprocessing", &
                        "Memory allocation failure for &
                        &commPatternMixing")

       ! Allocate the memory for the overset mesh communication pattern.
       ! This pattern changes in time and therefore each spectral time
       ! value has its own sliding mesh communication pattern.

       mm = nTimeIntervalsSpectral
       allocate(commPatternOverset(nn,mm), internalOverset(nn,mm), &
                stat=ierr)
       if(ierr /= 0)                     &
         call terminate("preprocessing", &
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
         call terminate("preprocessing", &
                        "Memory allocation failure for sendBuffer &
                        &and recvBuffer")

       ! Determine the cell range for the subfaces and initialize the
       ! arrays for the boundary condition data.
       ! Done for all grid levels.

       call cellRangeSubface
       call initBcdata

       ! Loop over the number of levels and perform a lot of tasks.
       ! See the corresponding subroutine header, although the
       ! names are pretty self-explaining

       do level=1,nLevels
         call xhalo(level)
         if (level == 1) then
           call oversetComm(level, .true., .false.)
         else
           call oversetComm(level, .true., .true.)
         end if
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
         call wallDistance(level, .true.)
       enddo

       ! Before heading to the solver, set all the boundary iblanks
       ! for the levels just updated to 0.

       do mm=1,nTimeIntervalsSpectral
         do level=1,nLevels
           do nn=1,nDom
             call setPointers(nn, level, mm)
             call changeIblanks(.false., 0_intType)
           end do
         end do
       end do

       end subroutine preprocessing

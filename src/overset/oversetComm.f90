!
!      ******************************************************************
!      *                                                                *
!      * File:          oversetComm.f90                                 *
!      * Author:        Steve Repsher                                   *
!      * Starting date: 04-09-2005                                      *
!      * Last modified: 10-18-2005                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine oversetComm(level, firstTime, coarseLevel)
!
!      ******************************************************************
!      *                                                                *
!      * oversetComm determines the communication pattern for the       *
!      * overset mesh boundaries on the given grid level.               *
!      *                                                                *
!      ******************************************************************
!
       use block
       use boundaryList
       use communication
       use inputIteration
       use inputOverset
       use inputTimeSpectral
       use iteration
       use searchMod
       implicit none
!
!      Subroutine arguments.
!
       integer(kind=intType), intent(in) :: level
       logical,               intent(in) :: firstTime, coarseLevel
!
!      Local variables.
!
       integer :: ierr

       integer(kind=intType) :: nn, sps, nHaloBndryAdd
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Release temporarily some memory such that the overall memory
       ! requirement is not dictated by this routine. What memory is
       ! released depends on firstTime. If firstTime is .True. This
       ! means that the routine is called in the preprocessing phase
       ! and thus only the send and receive buffer can be released.
       ! Otherwise this routine is called in a moving mesh computation
       ! and some more memory can be released.

       if( firstTime ) then
         deallocate(sendBuffer, recvBuffer, stat=ierr)
         if(ierr /= 0)                   &
           call terminate("oversetComm", &
                         "Deallocation error for communication buffers")
       else
         call deallocateTempMemory(.false.)
       endif

       ! Allocate memory for the derived type to help in coarse boundary
       ! creation only if this is a coarse level.
 
       if (coarseLevel) then
         allocate(blockBndry(nDom), stat=ierr)
         if(ierr /= 0)                   &
           call terminate("oversetComm", &
                         "Allocation failure for blockBndry")
       end if

       ! Set the fringe size, which is the number of cells in any
       ! direction that are used to compute a residual (e.g. just
       ! central differencing is a fringe size of 1).

       nn = ubound(flowDOms,2)
       fringeSize = 2
       if (level > mgStartLevel .and. level == nn) fringeSize = 1

       ! Set the interpolation type and number of weights based on the 
       ! current level. These are stored in the search module.

       if (level == 1) then
         interpolationType = oversetInterpType
         nInterp = nDonorWeights(oversetInterpType)
       else
         interpolationType = oversetInterpTypeCoarse
         nInterp = nDonorWeights(oversetInterpTypeCoarse)
       end if

       ! Loop over the number of spectral modes.

       spectralLoop: do sps=1,nTimeIntervalsSpectral

         ! If the overset connectivities are not changing in time, then
         ! all communication and boundary info is identical to the first
         ! mode, so just point everything to the first mode in this case.

         skipMode: if (.not. changingOverset .and. sps > 1) then

           ! Point the iblank and boundary information for the flowDoms.

           do nn = 1,nDom
             flowDoms(nn,level,sps)%nCellsOverset = &
             flowDoms(nn,level,  1)%nCellsOverset
             flowDoms(nn,level,sps)%nCellsOversetAll = &
             flowDoms(nn,level,  1)%nCellsOversetAll
             flowDoms(nn,level,sps)%nHoles = &
             flowDoms(nn,level,  1)%nHoles
             flowDoms(nn,level,sps)%nOrphans = &
             flowDoms(nn,level,  1)%nOrphans

             flowDoms(nn,level,sps)%iblank  => &
             flowDoms(nn,level,  1)%iblank
             flowDoms(nn,level,sps)%ibndry  => &
             flowDoms(nn,level,  1)%ibndry
             flowDoms(nn,level,sps)%idonor  => &
             flowDoms(nn,level,  1)%idonor
             flowDoms(nn,level,sps)%overint => &
             flowDoms(nn,level,  1)%overint

             flowDoms(nn,level,sps)%neighProcOver => &
             flowDoms(nn,level,  1)%neighProcOver
             flowDoms(nn,level,sps)%neighBlockOver => &
             flowDoms(nn,level,  1)%neighBlockOver
           end do

           ! Point the communication info.

           commPatternOverset(level,sps)%nProcSend = &
           commPatternOverset(level,  1)%nProcSend
           commPatternOverset(level,sps)%nProcRecv = &
           commPatternOverset(level,  1)%nProcRecv

           commPatternOverset(level,sps)%sendProc => &
           commPatternOverset(level,  1)%sendProc
           commPatternOverset(level,sps)%recvProc => &
           commPatternOverset(level,  1)%recvProc

           commPatternOverset(level,sps)%indexSendProc => &
           commPatternOverset(level,  1)%indexSendProc
           commPatternOverset(level,sps)%indexRecvProc => &
           commPatternOverset(level,  1)%indexRecvProc

           commPatternOverset(level,sps)%nSend => &
           commPatternOverset(level,  1)%nSend
           commPatternOverset(level,sps)%nRecv => &
           commPatternOverset(level,  1)%nRecv

           commPatternOverset(level,sps)%nSendCum => &
           commPatternOverset(level,  1)%nSendCum
           commPatternOverset(level,sps)%nRecvCum => &
           commPatternOverset(level,  1)%nRecvCum

           commPatternOverset(level,sps)%sendList => &
           commPatternOverset(level,  1)%sendList
           commPatternOverset(level,sps)%recvList => &
           commPatternOverset(level,  1)%recvList

           commPatternOverset(level,sps)%nPeriodic = &
           commPatternOverset(level,  1)%nPeriodic
           commPatternOverset(level,sps)%periodicData => &
           commPatternOverset(level,  1)%periodicData

           ! Point the internal communication.

           internalOverset(level,sps)%nCopy = &
           internalOverset(level,  1)%nCopy

           internalOverset(level,sps)%haloBlock => &
           internalOverset(level,  1)%haloBlock
           internalOverset(level,sps)%haloIndices => &
           internalOverset(level,  1)%haloIndices

           internalOverset(level,sps)%donorBlock => &
           internalOverset(level,  1)%donorBlock
           internalOverset(level,sps)%donorIndices => &
           internalOverset(level,  1)%donorIndices
           internalOverset(level,sps)%donorInterp => &
           internalOverset(level,  1)%donorInterp

           internalOverset(level,sps)%nPeriodic = &
           internalOverset(level,  1)%nPeriodic
           internalOverset(level,sps)%periodicData => &
           internalOverset(level,  1)%periodicData

           ! Cycle to the next mode.

           cycle

         end if skipMode

         ! Release the memory of the communication pattern and possibly
         ! the flow domain data if this is not the first time.

         if (.not. firstTime) call releaseMemOverset(level, sps, &
                                                     coarseLevel)

         ! How to build the overset boundary list depends if it is a
         ! coarse level or not.
 
         levelCheck: if (coarseLevel) then

           ! Since this is a coarse level, create the boundary based on
           ! the next finer level and build the fringe list.

           call createCoarseBoundary(level, sps, nHaloBndryAdd)
           call buildCoarseBoundaryList(level, sps, nHaloBndryAdd)

         else levelCheck

           ! This is a fine level so it is expected that the flow domain
           ! overset data is available. Build the list based on this data.

           call buildFineBoundaryList(level, sps)
 
         end if levelCheck

         ! Exchange the iblanks on the 1-to-1 halos. If a search is
         ! performed, then this exchange will be repeated as needed.

         call exchangeIblanks(level, sps, commPatternCell_2nd, &
                              internalCell_2nd)

         ! Determine the communication pattern (not quite final).

         call determineOversetComm(level, sps, coarseLevel)

         ! Check the boundary and overlap.

         call checkOverset(level, sps, level==1, .true., fringeSize)

         ! Add the halo fringe from 1-to-1 connecting blocks to the
         ! overset communication pattern. See the routine for details.

         call addHalosToBoundary(level, sps)

         ! Redetermine the size of the communication buffers.

         call setBufferSizes(level, sps, .false., .false., .true.)

       enddo spectralLoop

       ! Deallocate memory for the coarse boundary creation helper.
 
       if (coarseLevel) then
         deallocate(blockBndry, stat=ierr)
         if(ierr /= 0)                   &
           call terminate("oversetComm", &
                         "Deallocation error for blockBndry")
       end if
 
       ! Allocate the temporarily released memory again. For more info
       ! see the comments at the beginning of this routine.

       if( firstTime ) then
         allocate(sendBuffer(sendBufferSize), &
                  recvBuffer(recvBufferSize), stat=ierr)
         if(ierr /= 0)                    &
           call terminate("oversetComm", &
                          "Memory allocation failure for comm buffers")
       else
         call allocateTempMemory(.false.)
       endif

       end subroutine oversetComm

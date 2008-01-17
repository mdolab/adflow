!
!      ******************************************************************
!      *                                                                *
!      * File:          initFlow.f90                                    *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 03-06-2003                                      *
!      * Last modified: 09-13-2007                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine initFlow
!
!      ******************************************************************
!      *                                                                *
!      * initFlow sets the prescribed boundary data, allocates the      *
!      * memory for and initializes the flow variables. In case a       *
!      * restart is performed the owned variables are read from the     *
!      * previous solution file(s).                                     *
!      *                                                                *
!      ******************************************************************
!
       use bleedFlows
       use block
       use communication
       use inputIO
       use inputTimeSpectral
       use IOModule
       use iteration
       use restartMod
       implicit none
!
!      Local variables.
!
       integer :: ierr

       integer(kind=intType) :: sps, level, nLevels

       logical :: halosRead
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Determine the number of multigrid levels.

       nLevels = ubound(flowDoms,2)

       ! Allocate the memory for the prescribed boundary data at the
       ! boundary faces and determine the data for both the fine grid.

       groundLevel = 1
       call allocMemBcdata
       call setBCDataFineGrid(.true.)

       ! As some boundary conditions can be treated in multiple ways,
       ! some memory allocated must be released again.

       call releaseExtraMemBcs

       ! Determine the reference state.

       call referenceState

       ! Determine the prescribed data on the coarse grid levels
       ! by interpolation.

       call setBCDataCoarseGrid

       ! Set the infinity state.

       call setFlowInfinityState

       ! Non-dimensionalize the boundary data.

       call nonDimBoundData

       ! Set the primitive variables to the free stream values for
       ! the supersonic inflow faces for which this data has not
       ! been prescribed.

       call setSupersonicInletFreeStream

       ! Set the turbulent quantities to the free stream values for
       ! the inflow faces for which this data has not been prescribed.
       ! These can be both subsonic and supersonic inflow faces.

       call setInletFreestreamTurb

       ! Determine for the time spectral mode the matrices for the
       ! time derivatives.

       call timeSpectralMatrices

       ! Loop over the number of spectral solutions to allocate
       ! the memory for the w-variables and p on the fine grid.

       do sps=1,nTimeIntervalsSpectral
         call allocMemFlovarPart1(sps, 1_intType)
       enddo

       ! Initialize halosRead to .false. This will be overwritten
       ! if halo values are read during the restart.

       halosRead = .false.

       ! If a restart is performed, read the solution from file.

       testRestart: if( restart ) then

         ! Determine the number and names of the solution files and
         ! determine the contents of IOVar, the data structure to
         ! generalize the IO.

         call determineSolFileNames
         call setIOVar

         ! Determine the format of the files and read them.

         select case (fileFormatRead)
           case (cgnsFormat)
             call readRestartFile(halosRead)

           case (plot3DFormat)
             call readRestartFilePlot3D(halosRead)
         end select

         ! Copy or interpolate the spectral solution, if needed.

         if( copySpectral )     call copySpectralSolution
         if( interpolSpectral ) call interpolateSpectralSolution

         ! Release the memory of the file names and IOVar.

         deallocate(solFiles, IOVar, stat=ierr)
         if(ierr /= 0)                &
           call terminate("initFlow", &
                          "Deallocation failure for solFiles and IOVar")

         ! At the moment the pressure is stored at the location of the
         ! total energy. Copy the pressure to its own arrays and
         ! compute the total energy.

         call setPressureAndComputeEnergy(halosRead)

         ! If no halo values were read and inflow or outflow bleed
         ! regions are present, print a warning that the restart will
         ! not be perfect.

         if((.not. halosRead) .and. myID == 0) then
           if(nInflowBleeds > 0 .or. nOutflowBleeds > 0) then
             print "(a)", "#"
             print "(a)", "#                 Warning"
             print "(a)", "# Halo values could not read during the &
                          &restart."
             print "(a)", "# The restart is not consistent, because &
                          &the halo info is needed for a"
             print "(a)", "# consistent initialization of the bleed &
                          &regions."
             print "(a)", "# A constant extrapolation is used instead."
             print "(a)", "#"
           endif
         endif

       endif testRestart

       ! Initialize the bleed regions from the halos if the halos
       ! were read.

       if( halosRead ) call initBleedsFromHalos

       ! Allocate the memory for the solution variables on the coarse
       ! grid levels and the memory for the dependent flow variables,
       ! residuals, etc, on all multigrid levels.

       do sps=1,nTimeIntervalsSpectral
         call allocMemFlovarPart2(sps, 1_intType)

         do level=2,nLevels
           call allocMemFlovarPart1(sps, level)
           call allocMemFlovarPart2(sps, level)
         enddo
       enddo

       ! Initialize the halo cells if a restart is performed and the
       ! entire flow field if this is not the case.

       if( restart ) then
         call initializeHalos(halosRead)
       else
         call initFlowfield
       endif

       ! Initialize the prescribed boundary data for domain interfaces.
       ! This is just to avoid weird behavior in initDepvarAndHalos due
       ! to a wrong initialization. The actual values should be
       ! overwritten by the coupler.

       call initBCDataDomainInterfaces

       ! Initialize the dependent flow variables and the halo values.

       call initDepvarAndHalos(halosRead)

       end subroutine initFlow

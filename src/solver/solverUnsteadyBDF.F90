!
!      ******************************************************************
!      *                                                                *
!      * File:          solverUnsteadyBDF.F90                           *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 03-13-2003                                      *
!      * Last modified: 11-21-2007                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine solverUnsteadyBDF
!
!      ******************************************************************
!      *                                                                *
!      * solverUnsteadyBDF solves the unsteady equations using the BDF  *
!      * schemes for the multigrid level groundLevel.                   *
!      *                                                                *
!      ******************************************************************
!
       use inputIteration
       use inputUnsteady
       use iteration
       use killSignals
       use monitor
       implicit none
!
!      Local variables.
!
       integer(kind=intType) :: iter, nTimeSteps
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Initializations of the write parameters.

       writeVolume  = .false.
       writeSurface = .false.
       writeGrid    = .false.

       ! Initialize timeStepUnsteady to 0 and set the number of
       ! time steps depending on the grid level.

       timeStepUnsteady = 0

       nTimeSteps = nTimeStepsCoarse
       if(groundLevel == 1) nTimeSteps = nTimeStepsFine

       ! Loop over the number of time steps to be computed.

       timeLoop: do iter=1,nTimeSteps

         ! Perform the two initialization tasks for this time step.
         ! These are split into two, such that in python mode something
         ! could be added, i.e. an additional grid velocity.

         call initTimeStepPart1
         call initTimeStepPart2

         ! Solve the state for the current time step and 
         ! update nOldSolAvail.

         call solveState
         nOldSolAvail = nOldSolAvail + 1

         ! Update PV3 solution if PV3 support is required.

#ifdef USE_PV3
         call pv_update(real(timeStepUnsteady,realPV3Type))
#endif

         ! Determine whether or not solution files must be written.

         call checkWriteUnsteadyInLoop

         ! Exit the loop if the corresponding kill signal
         ! has been received.

         if(globalSignal == signalWriteQuit) exit

       enddo timeLoop

       ! Determine whether or not the final solution must be written.

       call checkWriteUnsteadyEndLoop

       end subroutine solverUnsteadyBDF

!      ==================================================================

       subroutine initTimeStepPart1
!
!      ******************************************************************
!      *                                                                *
!      * initTimeStepPart1 performs the first part of the               *
!      * initialization tasks before the actual computation of an       *
!      * unsteady time step is performed. It is split into two parts,   *
!      * such that some additional data can be set in multidisciplinary *
!      * mode, i.e. this routine may be used in python mode, although   *
!      * it is certainly possible to do this task via a python script.  *
!      *                                                                *
!      ******************************************************************
!
       use communication
       use inputMotion
       use inputUnsteady
       use iteration
       use monitor
       use section
       implicit none
!
!      Local variables.
!
       integer(kind=intType) :: nn

       real(kind=realType), dimension(nSections) :: tNewSec, deltaTSec
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Increment timeStepUnsteady and update
       ! timeUnsteady with the current time step.

       timeStepUnsteady = timeStepUnsteady + 1
       timeUnsteady     = timeUnsteady     + deltaT

       ! Write the unsteady header. Only done by processor 0
       ! to avoid a messy output.

       if(myID == 0) call unsteadyHeader

       ! If the grid is changing a whole lot of geometric
       ! info must be adapted.

       testChanging: if(changing_Grid .or. gridMotionSpecified) then

         ! Set the new time for all sections; also store their
         ! time step. They are the same for all sections, but all
         ! of them should be updated because of consistency.

         do nn=1,nSections
           tNewSec(nn)   = timeUnsteady + timeUnsteadyRestart
           deltaTSec(nn) = deltaT
         enddo

         ! Shift the coordinates and volumes and advance the
         ! coordinates 1 time step for deforming meshes.
         ! The shift only takes place for deforming meshes.

         if( deforming_Grid ) call shiftCoorAndVolumes
         call updateCoorFineMesh(deltaTSec, 1_intType)

         ! Adapt the geometric info on all grid levels needed for the
         ! current ground level and multigrid cycle.
         ! The wall distance only needs to be recomputed when the
         ! grid is changing; not when a rigid body motion is
         ! specified. Furthermore, the user can choose not to update
         ! the wall distance, because he may know a priori that the
         ! changes in geometry happen quite far away from the boundary
         ! layer. This is accomplished via updateWallDistanceUnsteady.
         ! Also note changingOverset can be true only when the grid
         ! is changing.

         call updateCoordinatesAllLevels
         if(changing_Grid .and. updateWallDistanceUnsteady) &
           call updateWallDistanceAllLevels

         if(changingOverset) call updateOversetAllLevels
         call updateSlidingAllLevels
         call updateMetricsAllLevels

         ! Update the rotation matrices of the faces. Only needed
         ! on the finest grid level.

         call faceRotationMatrices(currentLevel, .false.)

         ! Determine the velocities of the cell centers and faces
         ! for the current ground level. Note that the spectral mode
         ! is always 1 for unsteady mode.

         call gridVelocitiesFineLevel(deforming_Grid, tNewSec, 1_intType)

         ! Determine the new slip velocities on the viscous walls.

         call slipVelocitiesFineLevel(deforming_Grid, tNewSec, 1_intType)

       endif testChanging

       end subroutine initTimeStepPart1

!      ==================================================================

       subroutine initTimeStepPart2
!
!      ******************************************************************
!      *                                                                *
!      * initTimeStepPart2 performs the second part of the              *
!      * initialization tasks before the actual computation of an       *
!      * unsteady time step is performed. It is split into two parts,   *
!      * such that some additional data can be set in multidisciplinary *
!      * mode, i.e. this routine may be used in python mode, although   *
!      * it is certainly possible to do this task via a python script.  *
!      *                                                                *
!      ******************************************************************
!
       use inputMotion
       use iteration
       use monitor
       implicit none
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! If the grid is changing a whole lot of geometric
       ! info must be adapted.

       testChanging: if(changing_Grid .or. gridMotionSpecified) then

         ! Determine the velocities of the cell centers and faces and
         ! the slip velocities on the coarse grid levels . Note that the
         ! spectral mode is always 1 for unsteady mode.

         call gridVelocitiesCoarseLevels(1_intType)
         call slipVelocitiesCoarseLevels(1_intType)

         ! Compute the normal velocities of the boundaries, if
         ! needed for the corresponding boundary condition.

         call normalVelocitiesAllLevels(1_intType)

         ! Determine the prescribed boundary condition data for the
         ! boundary subfaces. First for the (currently) finest grid
         ! and then iteratively for the coarser grid levels.

         call setBCDataFineGrid(.false.)
         call setBCDataCoarseGrid

         ! Non-dimensionalize the boundary data.

         call nonDimBoundData

         ! Set the turbulent quantities to the free stream values for
         ! the inflow faces for which this data has not been
         ! prescribed. As the free stream values are nonDimensional,
         ! this call must be after the nonDimensionalization.

         call setInletFreestreamTurb

       endif testChanging

       ! Shift the old solution for the new time step.

       call shiftSolution

       ! Set the coefficients for the time integrator.

       call setCoefTimeIntegrator

       ! Determine the data for the outflow bleeds.

       call BCDataMassBleedOutflow(.false., .false.)

       end subroutine initTimeStepPart2

!      ==================================================================

       subroutine checkWriteUnsteadyInLoop
!
!      ******************************************************************
!      *                                                                *
!      * checkWriteUnsteadyInLoop checks if a solution must be          *
!      * written inside the time loop and if so write it.               *
!      *                                                                *
!      ******************************************************************
!
       use communication
       use inputIteration
       use inputMotion
       use iteration
       use killSignals
       use monitor
       implicit none
!
!      Local variables.
!
       integer :: ierr

       integer(kind=intType) :: nn
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Determine the global kill parameter if signals are supported.

#ifndef USE_NO_SIGNALS
       call mpi_allreduce(localSignal, globalSignal, 1, sumb_integer, &
                          mpi_max, SUmb_comm_world, ierr)
#endif

       ! Initialize the logicals for the writing to .false.

       writeVolume  = .false.
       writeSurface = .false.
       writeGrid    = .false.

       ! Check whether a solution file, either volume or surface, must
       ! be written. Only on the finest grid level in stand alone mode.

       if(standAloneMode .and. groundLevel == 1) then

         if(mod(timeStepUnsteady, nSaveVolume) == 0)  &
           writeVolume  = .true.
         if(mod(timeStepUnsteady, nSaveSurface) == 0) &
           writeSurface = .true.

         if(globalSignal == signalWrite .or. &
            globalSignal == signalWriteQuit) then
           writeVolume  = .true.
           writeSurface = .true.
         endif

         ! Determine whether or not a grid file must be written.

         if(changing_Grid .or. gridMotionSpecified) &
           writeGrid = writeVolume

         ! Write the solution.

         if(writeGrid .or. writeVolume .or. writeSurface) &
           call writeSol

       endif

       ! Update the variable oldSolWritten.

       do nn=(nOldLevels-1),2,-1
         oldSolWritten(nn) = oldSolWritten(nn-1)
       enddo

       oldSolWritten(1) = writeVolume

       end subroutine checkWriteUnsteadyInLoop

!      ==================================================================

       subroutine checkWriteUnsteadyEndLoop
!
!      ******************************************************************
!      *                                                                *
!      * checkWriteUnsteadyEndLoop checks if a solution must be         *
!      * written at the end of the time loop and if so write it.        *
!      *                                                                *
!      ******************************************************************
!
       use inputMotion
       use iteration
       use monitor
       implicit none
!
!      Local variables.
!
       integer(kind=intType) :: nn
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Write a solution. Only on the finest grid in stand alone mode
       ! and if some time steps have been performed since the last write.

       if(standAloneMode .and. groundLevel == 1) then

         writeVolume  = .not. writeVolume
         writeSurface = .not. writeSurface

         if(writeVolume .or. writeSurface) then

           ! Determine whether or not a grid file must be written.

           writeGrid = .false.
           if(changing_Grid .or. gridMotionSpecified) &
             writeGrid = writeVolume

           ! Write the solution.

           call writeSol

           ! Update the variable oldSolWritten.

           do nn=(nOldLevels-1),2,-1
             oldSolWritten(nn) = oldSolWritten(nn-1)
           enddo

           oldSolWritten(1) = writeVolume

         endif

       endif

       end subroutine checkWriteUnsteadyEndLoop

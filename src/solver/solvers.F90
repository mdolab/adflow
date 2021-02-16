module solvers
contains

  subroutine solver
    !
    !       solver is the main subroutine of the solver library.
    !       It controls the full multigrid and the kill signals.
    !
    use constants
    use communication, only : myid
    use inputDiscretization, only : eulerWallBCTreatment
    use inputIteration, only: mgStartLevel, printIterations
    use inputPhysics, only : equationMode
    use inputUnsteady, only : timeIntegrationScheme
    use killSignals, only : localSignal, noSignal
    use iteration, only : changing_grid, currentLevel, exchangePressureEarly, &
         groundLevel, nOldSolAvail, t0Solver
    use monitor, only : timeUnsteady
    use section, only : nSections
    use utils, only : eulerWallsPresent
    use multiGrid, only : transferToFineGrid
    use partitioning, only : updateCoorFineMesh
    implicit none
    !
    !      Local variables.
    !
    real(kind=realType), dimension(nSections) :: dtAdvance

    ! If the normal momentum equation should be used to determine
    ! the pressure in the halo for inviscid walls, find out if there
    ! actually are inviscid walls. If so, set the logical
    ! exchangePressureEarly to .true.; otherwise set it to .false.

    if(eulerWallBcTreatment == normalMomentum) then
       exchangePressureEarly = EulerWallsPresent()
    else
       exchangePressureEarly = .false.
    endif

    ! Connect the kill signals with the appropriate functions.
    ! Initialize localSignal for safety.
    ! Only if signalling is supported.

#ifndef USE_NO_SIGNALS
    localSignal = noSignal
    call connect_signals
#endif

    ! Determine the reference time for the solver.

    t0Solver = mpi_wtime()

    ! Set timeUnsteady to zero; this is amount of time simulated
    ! in unsteady mode.

    timeUnsteady = zero

    ! Loop over the number of grid levels in the current computation.
    ! Note that the counter variable, groundLevel, is defined in
    ! the module iteration.

    do groundLevel=mgStartlevel,1,-1

       ! Solve either the steady or the unsteady equations for this
       ! grid level. The time spectral method can be considered as
       ! a special kind of steady mode.
       select case (equationMode)
       case (steady, timeSpectral)
          call solveState
       case (unsteady)
          select case (timeIntegrationScheme)
          case (explicitRK)
             call solverUnsteadyExplicitRK
          end select
       end select

       ! If this is not the finest grid level, interpolate the
       ! solution to the next finer mesh and write a message that
       ! this is happening. Only processor 0 performs this latter
       ! task.

       if(groundLevel > 1) then

          currentLevel = groundLevel - 1

          if(myID == 0) then
             if (printIterations) then
                print "(a)", "#"
                print 100, currentLevel
                print "(a)", "#"
100             format("# Going down to grid level",1X,I1)
             end if
          endif

          call transferToFineGrid(.false.)

          ! Move the coordinates of the new fine grid level into the
          ! correct position. Only for unsteady problems with changing
          ! meshes. Note that the first argument of updateCoorFineMesh
          ! is an array with the time step for section.

          if(equationMode == unsteady .and. changing_Grid) then
             dtAdvance = timeUnsteady
             call updateCoorFineMesh(dtAdvance, 1_intType)
          endif

          ! Reset nOldSolAvail to 1, such that the unsteady
          ! computation on the finer mesh starts with a lower
          ! order scheme.

          nOldSolAvail = 1
       endif

    enddo

    ! Explictly set groundlevel to 1
    groundLevel = 1

  end subroutine solver

  ! ================================================
  ! Utilities for unsteady simulation
  ! ================================================
  subroutine solverUnsteadyInit
    !
    !       Initialize variables related to unsteady simulation.
    !       Some are the same as those in the *solver* subroutine,
    !       while others are specific for unsteady problems.
    !
    use ALEUtils, only : fillCoor, setLevelALE
    use constants
    use inputDiscretization, only : eulerWallBCTreatment
    use iteration, only : exchangePressureEarly, t0Solver
    use killSignals, only : localSignal, noSignal
    use monitor, only : timeUnsteady, timeStepUnsteady, writeVolume, writeSurface, writeGrid
    use utils, only : eulerWallsPresent
    implicit none

    ! BC treatment for normal momentum equation
    if(eulerWallBcTreatment == normalMomentum) then
       exchangePressureEarly = EulerWallsPresent()
    else
       exchangePressureEarly = .false.
    endif

    ! Connect the kill signals
#ifndef USE_NO_SIGNALS
    localSignal = noSignal
    call connect_signals
#endif

    ! Determine the reference time for the solver.
    t0Solver = mpi_wtime()

    ! Set time to zero
    timeUnsteady = zero
    timeStepUnsteady = 0

    ! Fill up old, xold and volold
    call fillCoor

    ! Set all ALE levels by initial configuration
    call setLevelALE(-1_intType)
  end subroutine solverUnsteadyInit

  subroutine updateUnsteadyGeometry
    !
    !       Update quantities related to geometry due to modification of mesh
    !       That could happen when
    !       - Steady mode with mesh modification
    !       - Unsteady mode with non-moving mesh but with prescribed mesh motion
    !       - Unsteady mode with warping and/or rigidly moving mesh
    !       - Unsteady mode coupled with an external solver
    !
    use ALEUtils, only : storeCoor, interpCoor, recoverCoor, setLevelALE, &
         slipVelocitiesFineLevel_ALE
    use bcdata, only : setbcdataFineGrid, setBCDataCoarseGrid
    use constants
    use inputMotion, only : gridMotionSpecified
    use inputUnsteady, only : deltaT, updateWallDistanceUnsteady, useALE
    use iteration, only : changing_grid, deforming_grid, currentLevel, groundLevel, nALEMeshes
    use monitor, only : timeUnsteady, timeUnsteadyRestart
    use partitioning, only : updateCoorFineMesh
    use preprocessingAPI, only : shiftCoorAndVolumes, metric, &
         updateCoordinatesAllLevels, updateMetricsAllLevels, faceRotationMatrices
    use section, only : nSections
    use solverUtils, only : gridVelocitiesFineLevel, gridVelocitiesCoarseLevels, &
         gridvelocitiesfinelevelpart1, gridvelocitiesfinelevelpart2, &
         normalVelocitiesAllLevels, slipVelocitiesFineLevel, slipVelocitiesCoarseLevels
    use wallDistance, only : updateWallDistanceAllLevels
    implicit none
    !
    !      Local variables.
    !
    integer(kind=intType) :: nn
    real(kind=realType), dimension(nSections) :: tNewSec, deltaTSec
    integer(kind=intType) :: lale

    testChanging: if(changing_Grid .or. gridMotionSpecified) then

       ! Set the new time for all sections
       do nn=1,nSections
          tNewSec(nn)   = timeUnsteady + timeUnsteadyRestart
          deltaTSec(nn) = deltaT
       enddo

       ! For prescribed motion only
       if (gridMotionSpecified) &
            call updateCoorFineMesh(deltaTSec, 1_intType)

       ! Adapt the geometric info on all grid levels needed for the
       ! current ground level and multigrid cycle.
       ! The wall distance only needs to be recomputed when the
       ! grid is changing; not when a rigid body motion is
       ! specified. Furthermore, the user can choose not to update
       ! the wall distance, because he may know a priori that the
       ! changes in geometry happen quite far away from the boundary
       ! layer. This is accomplished via updateWallDistanceUnsteady.

       call updateCoordinatesAllLevels
       if (changing_Grid .and. updateWallDistanceUnsteady) &
            call updateWallDistanceAllLevels
       call updateMetricsAllLevels

       ! Update the rotation matrices of the faces. Only needed
       ! on the finest grid level.
       ! For prescribed motion only

       if (gridMotionSpecified) &
            call faceRotationMatrices(currentLevel, .false.)

       if (useALE) then
          ! Update the velocities using ALE scheme if moving mesh is present

          ! First update cell and surface velocity, both are vectors
          ! Only quantities in blocks are updated, and they will not
          ! be interpolated

          call gridVelocitiesFineLevelPart1(deforming_Grid, tNewSec, 1_intType)

          ! Secondly store x to a temporary variable xALE

          call storeCoor

          ! Thirdly update surface normal and normal velocity

          ALEloop : do lale = 1, nALEMeshes
             ! Interpolate mesh over latest time step for all ALE Meshes
             call interpCoor(lale)

             ! Update s[I,J,K], norm
             call metric(groundLevel)

             ! Update sFace[I,J,K]
             call gridVelocitiesFineLevelPart2(deforming_Grid, tNewSec, 1_intType)

             ! Update uSlip
             call slipVelocitiesFineLevel_ALE(deforming_Grid, tNewSec, 1_intType)

             ! Update coarse level quantities to make sure multigrid is working
             call gridVelocitiesCoarseLevels(1_intType)
             call slipVelocitiesCoarseLevels(1_intType)

             ! Update rFace
             call normalVelocitiesAllLevels(1_intType)

             ! Store data to *lale* ALE level
             call setLevelALE(lale)
          enddo ALEloop

          ! Lastly recover x from temporary variable
          ! Then compute data for current level

          call recoverCoor

          ! Finish the rest of the update
          call metric(groundLevel)
          call gridVelocitiesFineLevelPart2(deforming_Grid, tNewSec, 1_intType)
          call slipVelocitiesFineLevel_ALE(deforming_Grid, tNewSec, 1_intType)

       else
          ! Otherwise update the velocities naively

          ! Determine the velocities of the cell centers and faces
          ! for the current ground level. Note that the spectral mode
          ! is always 1 for unsteady mode.

          call gridVelocitiesFineLevel(deforming_Grid, tNewSec, 1_intType)

          ! Determine the new slip velocities on the viscous walls.

          call slipVelocitiesFineLevel(deforming_Grid, tNewSec, 1_intType)
       endif

       ! After velocity computations on finest level are done,
       ! Update those on coarser levels

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

    endif testChanging

  end subroutine updateUnsteadyGeometry

  subroutine solverUnsteadyStep
    !
    !       Solve for next time step in unsteady simulation
    !
    use iteration, only : nOldSolAvail
    use solverUtils, only : shiftsolution
    use utils, only : setCoefTimeIntegrator
    use wallDistance, only : updateWallDistanceAllLevels
    implicit none

    ! Shift the old solution for the new time step.

    call shiftSolution

    ! Set the coefficients for the time integrator.

    call setCoefTimeIntegrator

    ! Solve the state for the current time step and
    ! update nOldSolAvail.
    call solveState
    nOldSolAvail = nOldSolAvail + 1

  end subroutine solverUnsteadyStep


  ! ================================================
  ! The following are not interfaced with Python
  ! ================================================

  subroutine checkWriteUnsteadyInLoop
    !
    !       checkWriteUnsteadyInLoop checks if a solution must be
    !       written inside the time loop and if so write it.
    !
    use communication, only : ADflow_comm_world
    use constants
    use inputIO, only : liftDistributionFile, sliceSolFile
    use inputIteration, only : nSaveSurface, nSaveVolume
    use inputMotion, only : gridMotionSpecified
    use iteration, only : changing_grid, groundLevel, nOldLevels, oldSolWritten
    use killSignals, only : localSignal, globalSignal, signalWrite, signalWriteQuit
    use monitor, only : nTimeStepsRestart, timeStepUnsteady, timeUnsteadyRestart, &
         writeVolume, writeSurface, writeGrid
    use tecplotIO, only : writeLiftDistributionFile, writeSlicesFile
    use surfaceFamilies, only : fullFamList
    implicit none
    !
    !      Local variables.
    !
    integer :: ierr

    integer(kind=intType) :: nn
    character(len=7) :: intString

    ! Determine the global kill parameter if signals are supported.

#ifndef USE_NO_SIGNALS
    call mpi_allreduce(localSignal, globalSignal, 1, adflow_integer, &
         mpi_max, ADflow_comm_world, ierr)
#endif

    ! Initialize the logicals for the writing to .false.

    writeVolume  = .false.
    writeSurface = .false.
    writeGrid    = .false.

    ! Check whether a solution file, either volume or surface, must
    ! be written. Only on the finest grid level in stand alone mode.

    !if(standAloneMode .and. groundLevel == 1) then
    if (groundLevel == 1) then
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
            call writeSol(fullFamList, size(fullFamList))

       ! Write the slice files for this timestep if they have been specified. TEMPORARY
       ! Write lift distribution TEMPORARY
       write(intString,"(i4.4)") timeStepUnsteady + nTimeStepsRestart
       intString = adjustl(intString)
       !call writeSlicesFile(trim(slicesolfile)//"_Timestep"//trim(intString)//".dat", .True.)
       !call writeLiftDistributionFile(trim(liftDistributionFile)//"_Timestep"//trim(intString)//".dat", .True.)

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
    !       checkWriteUnsteadyEndLoop checks if a solution must be
    !       written at the end of the time loop and if so write it.
    !
    use constants
    use inputMotion, only : gridMotionSpecified
    use iteration, only : changing_grid, groundLevel, nOldLevels, &
         oldSolWritten, standAloneMode
    use monitor, only : writeVolume, writeSurface, writeGrid
    use surfaceFamilies, only : fullFamList
    implicit none
    !
    !      Local variables.
    !
    integer(kind=intType) :: nn

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

          call writeSol(fullFamList, size(fullFamList))

          ! Update the variable oldSolWritten.

          do nn=(nOldLevels-1),2,-1
             oldSolWritten(nn) = oldSolWritten(nn-1)
          enddo

          oldSolWritten(1) = writeVolume

       endif

    endif

  end subroutine checkWriteUnsteadyEndLoop

  ! ================================================
  ! Utilities for explicit RK solver
  ! ================================================

  subroutine solverUnsteadyExplicitRK
    !
    !       solverUnsteadyExplicitRK solves the unsteady equations using
    !       the explicit Runge-Kutta schemes for the multigrid level
    !       groundLevel.
    !
    use constants
    use blockPointers, only: nDom, il, jl, kl, vol, dw, w, dwOldRK, p
    use communication
    use flowVarRefState
    use inputPhysics
    use inputUnsteady
    use iteration
    use killSignals
    use monitor
    use utils, only : setPointers
    use flowUtils, only : computePressure
    use haloExchange, only : whalo1, whalo2
    use turbutils, only : computeeddyviscosity
    use turbAPI, only : turbResidual
    use turbBCRoutines, only : applyAllTurbBC
    use utils, only : convergenceHeader
    use residuals, only :initRes, residual, sourceTerms
    use flowUtils, only : computeLamViscosity
    use BCRoutines, only : applyAllBC
    use preprocessingAPI, only : updateCoordinatesAllLevels
    implicit none
    !
    !      Local parameter.
    !
    integer(kind=intType), parameter :: nWriteConvHeader = 50
    !
    !      Local variables.
    !
    integer(kind=intType) :: iter, nTimeSteps, stage
    integer(kind=intType) :: i, j, k, l, m, nn

    real(kind=realType) :: tmp, timeUnsteadyOld

    character(len=7) :: numberString

    ! Set rkStage to 0. This variable is only relevant if Runge-Kutta
    ! smoothers are used as an iterative algorithm, but rkStage needs
    ! to be initialized.

    rkStage = 0

    ! Initializations of the write and signal parameters.

    writeVolume  = .false.
    writeSurface = .false.
    writeGrid    = .false.
    globalSignal = noSignal
    localSignal  = noSignal

    ! Initialize timeStepUnsteady to 0 and set the number of
    ! time steps depending on the grid level.

    timeStepUnsteady = 0

    nTimeSteps = nTimeStepsCoarse
    if(groundLevel == 1) nTimeSteps = nTimeStepsFine

    ! Write a message to stdout about how many time steps will
    ! be taken on the current grid. Only processor 0 does this.

    if(myID == 0) then
       write(numberString, "(i6)") nTimeSteps
       numberString = adjustl(numberString)

       print "(a)", "#"
       print 100, groundLevel, trim(numberString)
       print "(a)", "#"
100    format("# Grid",1x,i1,": Performing",1x,a,1x, &
            "explicit Runge Kutta time steps.")

       ! Also write the convergence header. Technically this is
       ! not really a convergence history for explicit RK, but
       ! the routine is called this way.

       call convergenceHeader
    endif

    ! Determine and write the initial convergence info. Again not
    ! really convergence for explicit RK.

    call convergenceInfo
    !
    !       Loop over the number of time steps to be computed.
    !
    ! Initialize the unsteady time.

    timeUnsteadyOld = timeUnsteady
    timeUnsteady    = timeUnsteadyOld + gammaRKUnsteady(1)*deltaT

    timeLoop: do iter=1,nTimeSteps

       ! Rewrite the convergence header after a certain number of
       ! iterations. Only processor 0 does this.

       if(myID == 0 .and. mod(iter,nWriteConvHeader) == 0) &
            call convergenceHeader

       ! Loop over the number of RK stages.

       stageLoop: do stage=1,nRKStagesUnsteady

          ! Compute the residual. Note that the turbulent residual must
          ! be computed first, because it uses the array of the flow
          ! field residuals as temporary buffers.

          if(equations == RANSEquations) then
             call initres(nt1, nt2)
             call turbResidual
          endif

          call initres(1_intType, nwf)
          call sourceTerms()
          call residual

          ! Loop over the number of domains.

          domainLoop: do nn=1,nDom

             ! Set the pointers to this domain. Note that there is only
             ! one time instance for a time accurate computation.

             call setPointers(nn,currentLevel,1)

             ! Step 1: Store in dw the conservative update for the
             !         flow field variables. This means a multiplication
             !         with deltaT/vol. Also store in w the conservative
             !         flow field variables.

             do k=2,kl
                do j=2,jl
                   do i=2,il
                      tmp = deltaT/vol(i,j,k)
                      do l=1,nw
                         dw(i,j,k,l) = tmp*dw(i,j,k,l)
                      enddo

                      w(i,j,k,ivx) = w(i,j,k,ivx)*w(i,j,k,irho)
                      w(i,j,k,ivy) = w(i,j,k,ivy)*w(i,j,k,irho)
                      w(i,j,k,ivz) = w(i,j,k,ivz)*w(i,j,k,irho)

                   enddo
                enddo
             enddo

             ! Step 2: Compute the new conservative flow field variables
             !         and primitive turbulence variables.

             do m=1,(stage-1)
                tmp = betaRKUnsteady(stage,m)
                if(tmp /= zero) then
                   do l=1,nw
                      do k=2,kl
                         do j=2,jl
                            do i=2,il
                               w(i,j,k,l) = w(i,j,k,l) - tmp*dwOldRK(m,i,j,k,l)
                            enddo
                         enddo
                      enddo
                   enddo
                endif
             enddo

             tmp = betaRKUnsteady(stage,stage)
             do l=1,nw
                do k=2,kl
                   do j=2,jl
                      do i=2,il
                         w(i,j,k,l) = w(i,j,k,l) - tmp*dw(i,j,k,l)
                      enddo
                   enddo
                enddo
             enddo

             ! Step 3. Convert the conservative variables back to
             !         primitive variables for this block. Compute the
             !         laminar and eddy viscosities as well.

             ! Convert the momentum into velocities.
             ! Also store the total energy in the pressure, such that
             ! the routine computePressure can be used.

             do k=2,kl
                do j=2,jl
                   do i=2,il
                      tmp          = one/w(i,j,k,irho)
                      w(i,j,k,ivx) = w(i,j,k,ivx)*tmp
                      w(i,j,k,ivy) = w(i,j,k,ivy)*tmp
                      w(i,j,k,ivz) = w(i,j,k,ivz)*tmp
                      p(i,j,k)     = w(i,j,k,irhoE)
                   enddo
                enddo
             enddo

             ! Compute the pressure.

             call computePressure(2_intType, il, 2_intType, jl, &
                  2_intType, kl, 0_intType)

             ! Swap the pressure and total energy, because
             ! computePressure stores the pressure in the position of
             ! the total energy.

             do k=2,kl
                do j=2,jl
                   do i=2,il
                      tmp            = p(i,j,k)
                      p(i,j,k)       = w(i,j,k,irhoE)
                      w(i,j,k,irhoE) = tmp
                   enddo
                enddo
             enddo

             ! Compute the laminar and eddy viscosities.

             call computeLamViscosity(.False.)
             call computeEddyViscosity(.False.)

             ! Step 4. Store dw in dwOldRK if this is not the last
             !         RK stage.

             if(stage < nRKStagesUnsteady) then
                do l=1,nw
                   do k=2,kl
                      do j=2,jl
                         do i=2,il
                            dwOldRK(stage,i,j,k,l) = dw(i,j,k,l)
                         enddo
                      enddo
                   enddo
                enddo
             endif

          enddo domainLoop

          ! Exchange the pressure if the pressure must be exchanged
          ! early. Only the first halo's are needed, thus whalo1 is
          ! called.

          call whalo1(currentLevel, 1_intType, 0_intType, .true.,&
               .false., .false.)

          ! Apply the boundary conditions to all blocks.

          call applyAllBC(.true.)
          if(equations == RANSEquations) call applyAllTurbBC(.true.)

          ! Exchange the halo data. As we are on the fine grid
          ! the second halo is needed.

          call whalo2(currentLevel, 1_intType, nw, .true., &
               .true., .true.)

          ! Determine the time and initialize the geometrical and
          ! boundary info for next stage, if needed.

          if(stage < nRKStagesUnsteady) then
             timeUnsteady = timeUnsteadyOld &
                  + gammaRKUnsteady(stage+1)*deltaT

             call initStageRK(stage+1)
          endif

       enddo stageLoop

       ! Increment timeStepUnsteady and update
       ! timeUnsteady with the current time step.

       timeStepUnsteady = timeStepUnsteady + 1
       timeUnsteady     = timeUnsteadyOld  + deltaT

       ! Write the convergence info.

       call convergenceInfo

       ! Determine the time and initialize the geometrical and
       ! boundary info for the next time step, such that the writing
       ! of grid files with moving geometries is done correctly.

       timeUnsteadyOld = timeUnsteady
       timeUnsteady    = timeUnsteadyOld + gammaRKUnsteady(1)*deltaT

       call initStageRK(1_intType)

       ! Determine whether or not solution files must be written.

       call checkWriteUnsteadyInLoop

       ! Exit the loop if the corresponding kill signal
       ! has been received.

       if(globalSignal == signalWriteQuit) exit

    enddo timeLoop

    ! Determine whether or not the final solution must be written.

    call checkWriteUnsteadyEndLoop

  end subroutine solverUnsteadyExplicitRK

  !=================================================================

  subroutine initStageRK(stage)
    !
    !       initStageRK performs the initialization tasks for the
    !       Runge-Kutta schemes in unsteady mode.
    !
    use inputMotion
    use inputUnsteady
    use iteration
    use monitor
    use section
    use BCData, only : setBCDataFineGrid
    use wallDistance, only : updateWallDistanceAllLevels
    use solverUtils
    use ALEUtils
    use preprocessingAPI, only :  updateMetricsAllLevels, &
         updateCoordinatesAllLevels, faceRotationMatrices
    use partitioning, only : updateCoorFineMesh
    implicit none
    !
    !      Subroutine arguments.
    !
    integer(kind=intType), intent(in) :: stage
    !
    !      Local variables.
    !
    integer(kind=intType) :: nn
    real(kind=realType)   :: fact

    real(kind=realType), dimension(nSections) :: tNewSec, deltaTSec

    ! If the grid is changing a whole lot of geometric
    ! info must be adapted.

    testChanging: if(changing_Grid .or. gridMotionSpecified) then

       ! Determine the time step relative to the previous stage.

       if(stage == 1) then
          fact = deltaT*(one - gammaRKUnsteady(nRKStagesUnsteady))
       else
          fact = deltaT*(gammaRKUnsteady(stage) &
               -         gammaRKUnsteady(stage-1))
       endif

       ! Set the time for all sections; also store their
       ! time step. They are the same for all sections, but all
       ! of them should be updated because of consistency.

       do nn=1,nSections
          tNewSec(nn)   = timeUnsteady + timeUnsteadyRestart
          deltaTSec(nn) = fact
       enddo

       ! Advance the coordinates 1 time step.

       call updateCoorFineMesh(deltaTSec, 1_intType)

       ! Adapt the geometric info on all grid levels needed for the
       ! current ground level and multigrid cycle.
       ! The wall distance only needs to be recomputed when the
       ! grid is changing; not when a rigid body motion is
       ! specified. Furthermore, the user can choose not to update
       ! the wall distance, because he may know a priori that the
       ! changes in geometry happen quite far away from the boundary
       ! layer. This is accomplished via updateWallDistanceUnsteady.

       call updateCoordinatesAllLevels
       if(changing_Grid .and. updateWallDistanceUnsteady) &
            call updateWallDistanceAllLevels

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

    ! Determine the prescribed boundary condition data for the
    ! boundary subfaces for the (currently) finest grid.

    call setBCDataFineGrid(.false.)

  end subroutine initStageRK

  ! ================================================
  ! Internal utilities
  ! ================================================
  subroutine solveState
    !
    !       solveState computes either the steady or unsteady state
    !       vector for the multigrid level groundLevel, which can be
    !       found in the module iteration.
    !
    use constants
    use communication, only : myID, adflow_comm_world
    use NKSolver, only : NKLSFuncEvals, freestreamResset, NK_LS, &
         NK_switchTol, useNKSolver, NK_CFL, getFreeStreamResidual, &
         getCurrentResidual, NKStep, computeResidualNK
    use anksolver, only : ANK_switchTol, useANKSolver, ANK_CFL, ANKStep, destroyANKSolver
    use inputio, only : forcedLiftFile, forcedSliceFile, forcedVolumeFile, &
         forcedSurfaceFile, solFile, newGridFile, surfaceSolFile
    use inputIteration, only: CFL, CFLCoarse, minIterNum, nCycles, &
         nCyclesCoarse, nMGSteps, nUpdateBleeds, printIterations, rkReset, timeLimit
    use iteration, only : cycling, approxTotalIts, converged, CFLMonitor, &
         groundLevel, iterTot, iterType, currentLevel, rhoRes0, totalR, t0Solver,&
         rhoResStart, totalR0, totalRFinal, totalRStart, stepMonitor, linResMonitor, ordersConverged
    use killSignals, only : globalSignal, localSignal, noSignal, routineFailed, signalWrite, &
         signalWriteQuit
    use monitor, only : writeGrid, writeSurface, writeVolume
    use utils, only: allocConvArrays, convergenceHeader
    use tecplotIO, only : writeTecplot
    use multiGrid, only : setCycleStrategy, executeMGCycle
    use surfaceFamilies, only : fullFamList
    use flowVarRefState, only : nwf
    use residuals, only : residual, initres, sourceTerms
    use solverUtils, only : timeStep
    implicit none
    !
    !      Local parameter
    !
    integer(kind=intType), parameter :: nWriteConvHeader = 50
    !
    !      Local variables.
    !
    integer :: ierr
    integer(kind=intType) ::  nMGCycles
    character (len=7) :: numberString
    logical :: absConv, relConv, firstNK, firstANK
    real(kind=realType) :: nk_switchtol_save, curTime, ordersConvergedOld

    ! Allocate the memory for cycling.
    if (allocated(cycling)) then
       deallocate(cycling)
    end if
    allocate(cycling(nMGSteps), stat=ierr)

    ! Some initializations.

    writeVolume  = .false.
    writeSurface = .false.
    converged    = .false.
    globalSignal = noSignal
    localSignal  = noSignal
    iterTot      = 0

    ! Determine the cycling strategy for this level.

    call setCycleStrategy

    ! Determine the number of multigrid cycles, which depends
    ! on the current multigrid level.

    nMGCycles = nCycles
    if(groundLevel > 1) nMGCycles = nCyclesCoarse

    ! Allocate (or reallocate) the convergence arry for this solveState.
    call allocConvArrays(nMGCycles)

    ! Allocate space for storing hisotry of function evaluations for NK
    ! solver with non-monotone line search if necessary
    if (currentLevel == 1 .and. useNKSolver .and. NK_LS==nonMonotoneLineSearch) then
       allocate(NKLSFuncEvals(nMGCycles))
    end if

    ! Only compute the free stream resisudal once for efficiency on the
    ! fine grid only.
    if (groundLevel == 1)  then
       if (.not. freeStreamResSet)  then
          call getFreeStreamResidual(rhoRes0, totalR0)
          freeStreamResSet = .True.
       end if
    end if
    ! Write a message. Only done by processor 0.

    if(myID == 0) then

       ! Write a message about the number of multigrid iterations
       ! to be performed.

       write(numberString,"(i6)") nMGCycles
       numberString = adjustl(numberString)
       numberString = trim(numberString)
       if (printIterations) then
          print "(a)", "#"
          print 102, groundLevel, trim(numberString),minIterNum,(NK_switchTol * totalR0)
          print "(a)", "#"
       end if
#ifndef USE_COMPLEX
102    format("# Grid",1X,I1,": Performing",1X,A,1X, &
            "iterations, unless converged earlier.",&
            " Minimum required iteration before NK switch: ",&
            I6,". Switch to NK at totalR of:",1X,e10.2)
#else
102    format("# Grid",1X,I1,": Performing",1X,A,1X, &
            "iterations, unless converged earlier.",&
            " Minimum required iteration before NK switch: ",&
            I6,". Switch to NK at totalR of:",1X,2e10.2)
#endif

       if (printIterations)  then
          call convergenceHeader
       end if
    end if


    ! Initialize the approxiate iteration count. We won't count this
    ! first residual evaluation. This way on the first convergence info
    ! call, "Iter" and "Iter Total" will both display 0.
    approxTotalIts = 0

    ! we need to re-set the orders converged to 16 as it might have been modified in the previous iteration
    ordersConverged = 16.0_realType

    ! Evaluate the initial residual
    call computeResidualNK

    ! Need to run the time step here since the RK/DADI is expecting
    ! the rad{i,j,k} to be computed.
    call timeStep(.False.)

    ! Extract the rhoResStart and totalRStart
    call getCurrentResidual(rhoResStart, totalRStart)

    ! No iteration type for first residual evaluation
    iterType = "    None"

    ! Determine and write the initial convergence info.
    call convergenceInfo

    ! Loop over the maximum number of nonlinear iterations
    firstANK = .True.
    firstNK = .True.

    ! Save the NKSwitch tol since it may be modified in the loop
    NK_SwitchTol_save = NK_switchtol

    ! Set the converged order factor now:
    ordersConverged = log10(totalR0/totalRStart)
    ordersConvergedOld = ordersConverged

    nonlinearIteration: do while (approxTotalIts < nMGCycles)

       ! Update iterTot
       iterTot = iterTot + 1

       if(mod(iterTot, nWriteConvHeader) == 0 .and. &
            myID == 0 .and. printIterations) then
          call convergenceHeader
       endif

       ! Defaults for the monitor
       CFLMonitor = CFL
       stepMonitor = 1.0
       linResMonitor = -1

       ! Determine what type of update to do:
       if (currentLevel > 1) then

          ! Coarse grids do RK/DADI always
          call executeMGCycle
          CFLMonitor = CFLCoarse
       else
          if (.not. useANKSolver .and. .not. useNKSolver .or. (iterTot <= minIterNum .and. rkreset)) then

             ! Always RK/DADI or a RK startup. Run the MG Cycle

             call executeMGCycle

          else if (useANKSolver .and. .not. useNKSolver) then

             ! Approx solver, no NKSolver

             if (totalR > ANK_switchTol * totalR0) then

                call executeMGCycle

             else
                call ANKStep(firstANK)
                firstANK = .False.
                CFLMonitor = ANK_CFL

             end if

          else if (.not. useANKSolver .and. useNKSolver) then

             ! NK Solver no approx solver

             if (totalR > NK_switchTol * totalR0) then

                call executeMGCycle

             else

                call NKStep(firstNK)
                firstNK = .False.
                CFLMonitor = -1

             end if

          else if (useANKSolver .and. useNKSolver) then

             ! Both approximate and NK solver.

             if (totalR > ANK_switchTol*totalR0) then

                call executeMGCycle

             else if (totalR <= ANK_switchTol*totalR0 .and. &
                  totalR > NK_switchTol*totalR0) then

                call ANKStep(firstANK)
                firstANK = .False.
                firstNK = .True.
                CFLMonitor = ANK_CFL

             else

                ! We free the memory for the ANK solver here because ANK solver
                ! module depends on the NK solver module and we cannot call
                ! destroyANKSolver within NKSolver itself.
                if (firstNK) then
                    call destroyANKSolver()
                end if

                call NKStep(firstNK)
                firstNK = .False.
                firstANK = .True.
                CFLMonitor = -1

             end if
          end if
       end if

       if (timeLimit > zero) then
          ! Check if we ran out of time but only if we are required to use the timeLimit
          if (myid == 0) then
             curTime = mpi_wtime() - t0solver
          end if

          call mpi_bcast(curTime, 1, adflow_real, 0, adflow_comm_world, ierr)

          if (curTime > timeLimit) then
             ! Set the iterTot to the limit directly so that convergence
             ! info thinks we are just out of cycles
             approxTotalIts = nMGCycles
          end if
       end if

       ! Determine and write the convergence info.
       call convergenceInfo

       totalRFinal = totalR

       ! Update how far we are converged:
       ordersConverged = max(log10(totalR0/totalR), ordersConvergedOld)
       ordersConvergedOld = ordersConverged

       ! Check for divergence or nan here
       if(routineFailed) then
          exit NonLinearIteration
       endif

       ! Exit the loop if we are converged
       if (converged) then
          exit nonLinearIteration
       end if

       ! Check if the bleed boundary conditions must be updated and
       ! do so if needed.


       ! Check if we've received a signal:
#ifndef USE_NO_SIGNALS
       call mpi_allreduce(localSignal, globalSignal, 1,         &
            adflow_integer, mpi_max, ADflow_comm_world, &
            ierr)
#endif

       if (globalSignal == signalWrite) then

          ! We have been told to write the solution

          writeGrid = .True.
          writeVolume = .True.
          writeSurface = .True.

          surfaceSolFile = forcedSurfaceFile
          newGridFile = forcedVolumeFile
          solFile = forcedVolumeFile

          call writeSol(fullFamList, size(fullFamList))

          ! Also write potential tecplot files. Note that we are not
          ! writing the tecplot surface file so pass in an empty file
          ! name and a nonsense family list.
          call writeTecplot(forcedSliceFile, .True., forcedLiftFile, .True., &
               "", .False., [0], 1)

          ! Reset the signal
          localSignal = noSignal
       end if

       if (globalSignal == signalWriteQuit) then
          exit nonLinearIteration
       end if


    enddo nonLinearIteration

    ! Restore the switch tol in case it was changed
    NK_switchtol = NK_SwitchTol_save

    ! deallocate space for storing hisotry of function evaluations for NK
    ! solver with non-monotone line search if necessary
    if (currentLevel == 1 .and. useNKSolver .and. NK_LS==nonMonotoneLineSearch) then
       deallocate(NKLSFuncEvals)
    end if

  end subroutine solveState


  subroutine convergenceInfo
    !
    !       convergenceInfo computes and writes the convergence info to
    !       standard output. In spectral mode a convergence history for
    !       every mode is written.
    !
    use constants
    use cgnsNames
    use block, only : nCellGlobal
    use blockPointers, only : nDom
    use communication, only : adflow_comm_world, myid
    use inputIteration, only : printIterations, l2convcoarse, l2conv, l2convrel, &
         minIterNum, maxL2DeviationFactor, ncycles, RKReset
    use inputPhysics, only : liftDirection, dragDirection, equationMode, &
         lengthRef, machCoef, surfaceRef
    use flowVarRefState, only : pRef, lRef, gammaInf
    use inputIO, only : storeConvInnerIter
    use inputTimeSpectral, only : nTimeIntervalsSpectral
    use inputUnsteady, only : timeIntegrationScheme
    use monitor, only : monLoc, monGlob, nMon, nMonMax, nMonSum, monNames, timeDataArray, &
         showCPU, monRef, convArray, timeUnsteadyRestart, timeArray, timeStepUnsteady, &
         timeUnsteady, nTimeStepsRestart
    use iteration, only : groundLevel, currentLevel, iterTot, iterType, approxTotalIts, &
         CFLMonitor, stepMonitor, t0solver, converged, linResMonitor
    use killSignals, only : routineFailed, fromPython
    use iteration, only : rhoRes, rhoResStart, totalR, totalRStart, totalR0
    use oversetData, only: oversetPresent
    use utils, only : setPointers, myisnan, returnFail, maxHDiffMach, maxEddyv, &
         sumResiduals, sumAllResiduals
    use surfaceIntegrations, only : integrateSurfaces
    use zipperIntegrations, only : integrateZippers
    use surfaceFamilies, only : fullFamLIst
    implicit none
    !
    !      Local variables.
    !
    integer :: ierr, iConvStdout

    integer(kind=intType) :: sps, nn, mm, iConv

    real(kind=realType) :: hdiffMax, MachMax
    real(kind=realType) :: eddyvisMax, yplusMax, sepSensor, Cavitation, axisMoment
    real(kind=realType) :: sepSensorAvg(3)

    real(kind=realType) :: L2ConvThisLevel, fact
    real(kind=realType), dimension(3) :: cfp, cfv, cmp, cmv
    real(kind=realType) :: cmpaxis, cmvaxis
    logical :: nanOccurred, writeIterations
    logical :: absConv, relConv
    real(kind=realType) :: localValues(nLocalValues)
    real(kind=realType) :: funcValues(nCostFunction)
    ! Determine whether or not the iterations must be written.

    writeIterations = .true.
    if(equationMode          == unsteady .and. &
         timeIntegrationScheme == explicitRK) writeIterations = .false.

    ! Initializations
    converged   = .False.
    nanOccurred = .false.

    ! Set the L2 norm for convergence for this level.

    L2ConvThisLevel = L2ConvCoarse
    if(groundLevel == 1) L2ConvThisLevel = L2Conv

    ! Loop over the number of spectral solutions.

    spectralLoop: do sps=1,nTimeIntervalsSpectral

       ! Initialize the local monitoring variables to zero.

       monLoc = zero

       ! Loop over the blocks.

       domains: do nn=1,nDom

          ! Set the pointers for this block.

          call setPointers(nn, groundLevel, sps)

          ! Compute the forces and moments for this block.  Note that
          ! we zero localValues before each call becuase we are
          ! summing into momLocal.
          localvalues = zero
          call integrateSurfaces(localValues, fullFamList)

          ! Convert to coefficients for monitoring:
          fact = two/(gammaInf*MachCoef*MachCoef &
               *surfaceRef*LRef*LRef*pRef)
          cfp = fact*localValues(iFp:iFp+2)
          cfv = fact*localValues(iFv:iFv+2)
          fact = fact/(lengthRef*Lref)
          cmp = fact*localValues(iMp:iMp+2)
          cmv = fact*localValues(iMv:iMv+2)
          yplusmax = localValues(iYplus)
          ! Determine the maximum values of the monitoring variables
          ! of this block.

          call maxHdiffMach(hdiffMax, MachMax)
          call maxEddyv(eddyvisMax)

          ! Loop over the number of monitoring variables.
          nMonitoringVar: do mm=1,nMon

             ! Determine the monitoring variable and act accordingly.

             select case (monNames(mm))

             case ('totalR')
                call sumAllResiduals(mm)

             case (cgnsL2resRho)
                call sumResiduals(irho, mm)

             case (cgnsL2resMomx)
                call sumResiduals(imx, mm)

             case (cgnsL2resMomy)
                call sumResiduals(imy, mm)

             case (cgnsL2resMomz)
                call sumResiduals(imz, mm)

             case (cgnsL2resRhoe)
                call sumResiduals(irhoE, mm)

             case (cgnsL2resNu, cgnsL2resK)
                call sumResiduals(itu1, mm)

             case (cgnsL2resOmega, cgnsL2resTau, cgnsL2resEpsilon)
                call sumResiduals(itu2, mm)

             case (cgnsL2resV2)
                call sumResiduals(itu3, mm)

             case (cgnsL2resF)
                call sumResiduals(itu4, mm)

             case (cgnsCl)
                monLoc(mm) = monLoc(mm)                         &
                     + (cfp(1) + cfv(1))*liftDirection(1) &
                     + (cfp(2) + cfv(2))*liftDirection(2) &
                     + (cfp(3) + cfv(3))*liftDirection(3)

             case (cgnsClp)
                monLoc(mm) = monLoc(mm) + cfp(1)*liftDirection(1) &
                     +              cfp(2)*liftDirection(2) &
                     +              cfp(3)*liftDirection(3)

             case (cgnsClv)
                monLoc(mm) = monLoc(mm) + cfv(1)*liftDirection(1) &
                     +              cfv(2)*liftDirection(2) &
                     +              cfv(3)*liftDirection(3)

             case (cgnsCd)
                monLoc(mm) = monLoc(mm)                         &
                     + (cfp(1) + cfv(1))*dragDirection(1) &
                     + (cfp(2) + cfv(2))*dragDirection(2) &
                     + (cfp(3) + cfv(3))*dragDirection(3)

             case (cgnsCdp)
                monLoc(mm) = monLoc(mm) + cfp(1)*dragDirection(1) &
                     +              cfp(2)*dragDirection(2) &
                     +              cfp(3)*dragDirection(3)

             case (cgnsCdv)
                monLoc(mm) = monLoc(mm) + cfv(1)*dragDirection(1) &
                     +              cfv(2)*dragDirection(2) &
                     +              cfv(3)*dragDirection(3)

             case (cgnsCfx)
                monLoc(mm) = monLoc(mm) + cfp(1) + cfv(1)

             case (cgnsCfy)
                monLoc(mm) = monLoc(mm) + cfp(2) + cfv(2)

             case (cgnsCfz)
                monLoc(mm) = monLoc(mm) + cfp(3) + cfv(3)

             case (cgnsCmx)
                monLoc(mm) = monLoc(mm) + cmp(1) + cmv(1)

             case (cgnsCmy)
                monLoc(mm) = monLoc(mm) + cmp(2) + cmv(2)

             case (cgnsCmz)
                monLoc(mm) = monLoc(mm) + cmp(3) + cmv(3)

             case (cgnsHdiffMax)
                monLoc(mm) = max(monLoc(mm), hdiffMax)

             case (cgnsMachMax)
                monLoc(mm) = max(monLoc(mm), MachMax)

             case (cgnsYplusMax)
                monLoc(mm) = max(monLoc(mm), yplusMax)

             case (cgnsEddyMax)
                monLoc(mm) = max(monLoc(mm), eddyvisMax)

             case (cgnsSepSensor)
                monLoc(mm) = monLoc(mm) + localValues(isepSensor)

             case (cgnsCavitation)
                monLoc(mm) = monLoc(mm) + localValues(iCavitation)

             case (cgnsAxisMoment)
                monLoc(mm) = monLoc(mm) + localValues(iaxisMoment)

             end select

          end do nMonitoringVar
       end do domains

       ! Add the corrections from zipper meshes from proc 0
       if (oversetPresent) then
          localValues = zero
          call integrateZippers(localValues, fullFamList, sps)

          fact = two/(gammaInf*MachCoef*MachCoef &
               *surfaceRef*LRef*LRef*pRef)
          cfp = localValues(iFp:iFp+2)*fact
          cfv = localValues(iFv:iFv+2)*fact
          fact = fact/(lengthRef*Lref)
          cmp = localValues(iMp:iMp+2)*fact
          cmv = localValues(iMv:iMv+2)*fact

          !Loop over the number of monitoring variables and just modify
          !the ones that need to be updated with the zipper forces we just
          !computed.
          nMonitoringVarZip: do mm=1,nMon

             ! Determine the monitoring variable and act accordingly.

             select case (monNames(mm))

             case (cgnsCl)
                monLoc(mm) = monLoc(mm)                         &
                     + (cfp(1) + cfv(1))*liftDirection(1) &
                     + (cfp(2) + cfv(2))*liftDirection(2) &
                     + (cfp(3) + cfv(3))*liftDirection(3)

             case (cgnsClp)
                monLoc(mm) = monLoc(mm) + cfp(1)*liftDirection(1) &
                     +              cfp(2)*liftDirection(2) &
                     +              cfp(3)*liftDirection(3)

             case (cgnsClv)
                monLoc(mm) = monLoc(mm) + cfv(1)*liftDirection(1) &
                     +              cfv(2)*liftDirection(2) &
                     +              cfv(3)*liftDirection(3)

             case (cgnsCd)
                monLoc(mm) = monLoc(mm)                         &
                     + (cfp(1) + cfv(1))*dragDirection(1) &
                     + (cfp(2) + cfv(2))*dragDirection(2) &
                     + (cfp(3) + cfv(3))*dragDirection(3)

             case (cgnsCdp)
                monLoc(mm) = monLoc(mm) + cfp(1)*dragDirection(1) &
                     +              cfp(2)*dragDirection(2) &
                     +              cfp(3)*dragDirection(3)

             case (cgnsCdv)
                monLoc(mm) = monLoc(mm) + cfv(1)*dragDirection(1) &
                     +              cfv(2)*dragDirection(2) &
                     +              cfv(3)*dragDirection(3)

             case (cgnsCfx)
                monLoc(mm) = monLoc(mm) + cfp(1) + cfv(1)

             case (cgnsCfy)
                monLoc(mm) = monLoc(mm) + cfp(2) + cfv(2)

             case (cgnsCfz)
                monLoc(mm) = monLoc(mm) + cfp(3) + cfv(3)

             case (cgnsCmx)
                monLoc(mm) = monLoc(mm) + cmp(1) + cmv(1)

             case (cgnsCmy)
                monLoc(mm) = monLoc(mm) + cmp(2) + cmv(2)

             case (cgnsCmz)
                monLoc(mm) = monLoc(mm) + cmp(3) + cmv(3)

             end select

          end do nMonitoringVarZip
       end if
       ! Determine the global sum of the summation monitoring
       ! variables. This is an all reduce since every processor needs to
       ! know the residual to make the same descisions.

       if(nMonSum > 0) &
            call mpi_allreduce(monLoc, monGlob, nMonSum, adflow_real, &
            mpi_sum, ADflow_comm_world, ierr)

       ! Idem for the maximum monitoring variables.
#ifndef USE_COMPLEX
       if(nMonMax > 0) &
            call mpi_allreduce(monLoc(nMonSum+1), monGlob(nMonSum+1), &
            nMonMax, adflow_real, mpi_max, ADflow_comm_world, ierr)
#else
       if (nMonMax < 0) &
            monGlob(nMonSum+1) = zero
#endif

       ! Write the convergence info; only processor 0 does this.

       testRootProc: if(myID == 0) then

          ! The variables which must always be written.

          if(printIterations) then
             write(*,"(1x,i6,2x)",advance="no") groundLevel

             if(equationMode == unsteady) then

                write(*,"(i6,1x)",advance="no") timeStepUnsteady + &
                     nTimeStepsRestart
                write(*,"(e12.5,1x)",advance="no") timeUnsteady + &
                     timeUnsteadyRestart

             else if(equationMode == timeSpectral) then

                write(*,"(i8,3x)",advance="no") sps

             endif

             if( writeIterations ) then
                write(*,"(i6,1x)",advance="no") iterTot
                write(*,"(i6,1x)",advance="no") approxTotalIts
                write(*,"(a,1x)", advance="no") iterType

                if (CFLMonitor < zero) then
                   ! Print dashes if no cfl term is used, i.e. NK solver
                   write(*,"(a,1x)", advance="no") "    ----  "
                else
#ifndef USE_COMPLEX
                   write(*,"(e10.2,1x)",advance="no") CFLMonitor
#else
                   write(*,"(e10.2,1x)",advance="no") real(CFLMonitor)
#endif
                end if
#ifndef USE_COMPLEX
                write(*,"(f5.2,2x)",advance="no") stepMonitor
#else
                write(*,"(f5.2,2x)",advance="no") real(stepMonitor)
#endif
                if (linResMonitor < zero) then
                   ! For RK/DADI just print dashes
                   write(*,"(a,1x)", advance="no") " ----"
                else

#ifndef USE_COMPLEX
                   write(*,"(f5.3,1x)",advance="no") linResMonitor
#else
                   write(*,"(f5.3,1x)",advance="no") real(linResMonitor)
#endif
                end if

                if( showCPU ) then
#ifndef USE_COMPLEX
                   write(*,"(e12.5,1x)",advance="no") mpi_wtime() - t0Solver
#else
                   write(*,"(e12.5,1x)",advance="no") real(mpi_wtime() - t0Solver)
#endif
                end if
             end if
          end if
       end if testRootProc

       ! Loop over the number of monitoring values.
       variableLoop: do mm=1, nMon

          ! The residual variables must be corrected.

          select case (monNames(mm))

          case (cgnsL2resRho,  cgnsL2resMomx,    &
               cgnsL2resMomy, cgnsL2resMomz,    &
               cgnsL2resRhoe, cgnsL2resNu,      &
               cgnsL2resK,    cgnsL2resOmega,   &
               cgnsL2resTau,  cgnsL2resEpsilon, &
               cgnsL2resV2,   cgnsL2resF        )
             monGlob(mm) = sqrt(monGlob(mm)/nCellGlobal(groundLevel))
             if (monNames(mm) == cgnsL2resRho) then
                rhoRes = monGlob(mm)
             end if
          case ('totalR')
             monGlob(mm) = sqrt(monGlob(mm))
             totalR = monGlob(mm)
          end select

          if( myIsNAN(monGlob(mm)) ) nanOccurred = .true.

          if (myid == 0 .and. printIterations) then
             ! Write the convergence info to stdout.
#ifndef USE_COMPLEX
             write(*,"(e24.16,1x)",advance="no") monGlob(mm)
#else
             write(*,"(2e24.16,1x)",advance="no") monGlob(mm)
#endif
          end if

          ! Store the convergence info in convArray, if desired.
          if( storeConvInnerIter ) then
             convArray(iterTot, sps, mm) = monGlob(mm)
          endif
       end do variableLoop

       if (myid == 0 .and. printIterations) then
          ! Write the carriage return.
          print "(1x)"
       end if

       ! Determine whether or not the solution is converged.
       ! A distinction must be made between unsteady mode and the
       ! other modes, because in unsteady mode the convergence of
       ! the inner iterations is not necessarily stored.

       select case (equationMode)
       case (steady, timeSpectral)

          ! Steady or time spectral mode. The convergence histories
          ! are stored and this info can be used. The logical
          ! converged is set to .false. if the density residual
          ! has not converged yet.

          absConv = .False.
          relConv = .False.

          ! We make a split here based on if we are operating on a
          ! coarse grid or the finest. On the coarse grid levels, we
          ! l2convThisLevel refers to the relative convergence of
          ! rhoRes.

          if (currentLevel /= 1) then
             if (rhoRes < L2ConvThisLevel * rhoResStart) then
                relConv = .True.
             end if
          else
             ! We are on the fine level. All checking is done using
             ! the total residual.

             ! Absolute convergence check
             if (totalR < L2ConvThisLevel * totalR0) then
                absConv = .True.
             end if

             ! Relative check only done on finest level
             if (totalR < L2ConvRel*totalRStart) then
                relConv = .True.
             end if

             ! If the totla number of iterations is less than the
             ! RKReset, don't check the residual
             if (iterTot < minIterNum .and. rkreset) then
                relConv = .False.
                absConv = .False.
             end if
          end if

          ! Combine the two flags.
          if (absConv .or. relConv)  then
             converged = .True.
          end if

          !===========================================================

       case (unsteady)

          ! Unsteady mode. The array convArray may not be present
          ! and therefore something else must be done.
          ! First determine the position in the array timeDataArray
          ! that can be used. For the coarser grids this is 1,
          ! because the time evolution is overwritten. For the fine
          ! mesh the actual position is determined.

          nn = 1
          if(groundLevel == 1) &
               nn = timeStepUnsteady + nTimeStepsRestart
          nn = max(nn,1_intType)

          ! Make a distinction between the time integration
          ! schemes.

          select case(timeIntegrationScheme)

          case (explicitRK)

             ! Explicit scheme. Simply store the data in the
             ! convergence arrays.

             timeArray(nn) = timeUnsteady + timeUnsteadyRestart

             ! For explicit schemes the residuals are not
             ! monitored and therefore the monitoring variables
             ! can simply be copied.

             do mm=1,nMon
                timeDataArray(nn,mm) = monGlob(mm)
             enddo

             !=======================================================

          case (BDF,implicitRK)

             ! An implicit scheme is used and therefore an
             ! iterative algorithm within every time step.
             ! The array convArray may not be present and
             ! therefore
             ! something else must be done. First determine the
             ! position in the array timeDataArray that can be
             ! used.
             ! For the coarser grids this is 1, because the time
             ! evolution is overwritten. For the fine mesh the
             ! actual position is determined.

             ! Determine the situation we have here.

             testInitUnsteady: if(iterTot == 0) then

                ! This is the initialization phase for this time step.
                ! Simply copy monGlob into monRef, store the value
                ! of the physical time and set converged to .false..


                do mm=1,nMon
                   monRef(mm) = monGlob(mm)
                enddo

                timeArray(nn) = timeUnsteady + timeUnsteadyRestart
                converged     = .false.

             else testInitUnsteady

                ! Iteration for this time step. Store the relative
                ! convergence for the residual compared to the start
                ! of this time step; an absolute norm does not give
                ! any information here. For all other monitoring
                ! variables the current value is stored.

                do mm=1,nMon
                   select case (monNames(mm))

                   case (cgnsL2resRho,  cgnsL2resMomx,    &
                        cgnsL2resMomy, cgnsL2resMomz,    &
                        cgnsL2resRhoE, cgnsL2resNu,      &
                        cgnsL2resK,    cgnsL2resOmega,   &
                        cgnsL2resTau,  cgnsL2resEpsilon, &
                        cgnsL2resV2,   cgnsL2resF)

                      timeDataArray(nn,mm) = monGlob(mm) &
                           / max(monRef(mm),eps)

                      !===============================================

                   case default


                      timeDataArray(nn,mm) = monGlob(mm)

                   end select
                enddo

                ! Set the logical converged to .false. if the density
                ! residual has not converged yet.

                if(timeDataArray(nn,1) < L2ConvThisLevel) &
                     converged = .True.

             endif testInitUnsteady
          end select ! temporal integration scheme

       end select ! unsteady

    end do spectralLoop


    if( nanOccurred )then

       if (myid == 0) then
          print *,'Nan occured in Convergence Info on proc:',myid
       end if

       routineFailed = .True.

       call returnFail("convergenceInfo", &
            "A NaN occurred during the computation.")

       ! in a normal computation, code will simply exit.
       ! in a python based computation, code will set
       ! routinedFailed to .True. and return to the
       ! python level...
       return
    endif

    ! ! If we are at the max iteration limit but the residual is
    ! ! *close*, ie within maxL2DeviationFactor we say that's fine

    if(fromPython .and. groundLevel ==1 .and. approxTotalIts >= nCycles) then

       !Check to see if residuals are diverging or stalled for python
       select case (equationMode)

       case (steady, timeSpectral)

          ! Steady or time spectral mode. The convergence histories
          ! are stored and this info can be used. If the residuals
          ! are diverging the, logical routineFailed in killSignals
          ! is set to true and the progress is halted.
          !only check on root porcessor
          if (myID == 0) then

             ! If we made it to ncycles, check to see if we're
             ! "close" to being converged.
             do sps = 1,nTimeIntervalsSpectral
                if (totalR > maxL2DeviationFactor * totalR0 * L2ConvThisLevel) then
                   routineFailed = .True.
                end if
             enddo
          endif
          call mpi_bcast(routineFailed, 1, MPI_LOGICAL, 0, ADflow_comm_world, ierr)

       end select
    end if

    ! Determine whether or not the solution is considered
    ! converged.  This info is only known at processor 0 and must
    ! therefore be broadcast to the other processors.

    call mpi_bcast(converged, 1, MPI_LOGICAL, 0, ADflow_comm_world, ierr)

  end subroutine convergenceInfo

end module solvers

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
  !      ******************************************************************
  !      *                                                                *
  !      * solverUnsteadyWrap is a wrapper of solverUnsteady_ALE for      *
  !      * MD coupling at python level.                                   *
  !      * The original solver is dismantled to seperate out              *
  !      * time-loop-related procedures.                                  *
  !      * Expected structure in Python coupler:                          *
  !      * - solverUnsteadyWrapBegin                                      *
  !      * - Update physical time step                                    *
  !      * - Time loop begins                                             *
  !      *   - Time increment                                             *
  !      *   - Mesh deformation                                           *
  !      *   - solverUnsteadyWrapInLoop                                   *
  !      * - Time loop ends                                               *
  !      * - solverUnsteadyWrapEnd: checkWriteUnsteadyEndLoop             *
  !      *                                                                *
  !      ******************************************************************
  !

! ******************************************************************
! Wrappers
! ******************************************************************
subroutine solverUnsteadyWrapBegin
  use inputIteration
  use inputUnsteady
  use blockPointers
  use inputDiscretization
  use iteration
  use killSignals
  use monitor
  use communication 
  use section
  use inputMotion
  use inputTimeSpectral
  implicit none
  !
  !      Local variables.
  
  logical :: EulerWallsPresent
  
  ! If the normal momentum equation should be used to determine
  ! the pressure in the halo for inviscid walls, find out if there
  ! actually are inviscid walls. If so, set the logical
  ! exchangePressureEarly to .true.; otherwise set it to .false.
  
  if(wallBcTreatment == normalMomentum) then
     exchangePressureEarly = EulerWallsPresent()
  else
     exchangePressureEarly = .false.
  endif
  
  ! Determine the reference time for the solver.
  
  t0Solver = mpi_wtime()
  
  ! Set timeUnsteady to zero; this is amount of time simulated
  ! in unsteady mode.
  
  ! timeUnsteady = zero

  ! Initializations of the write parameters.

  writeVolume  = .false.
  writeSurface = .false.
  writeGrid    = .false.

  timeUnsteady     = 0.0
  timeStepUnsteady = 0

  ! Fill up old xold and volold
  call fillCoor

  ! Set all ALE levels by initial configuration
  call setLevelALE(-1_intType)

end subroutine solverUnsteadyWrapBegin

! ******************************************************************
subroutine solverUnsteadyWrapInLoop
  use iteration
  implicit none

  call initTimeStepWrap

  ! Solve the state for the current time step and 
  ! update nOldSolAvail.

  call solveState

  nOldSolAvail = nOldSolAvail + 1

  ! Determine whether or not solution files must be written.
  ! This will be done in Python
  ! call checkWriteUnsteadyInLoop

  ! --------------------------------
  ! Unnecessary since there is no loop
  ! --------------------------------
  !if(globalSignal == signalWriteQuit) exit     

end subroutine solverUnsteadyWrapInLoop

! ******************************************************************
subroutine solverUnsteadyWrapEnd
  ! Determine whether or not the final solution must be written.
  call checkWriteUnsteadyEndLoop
end subroutine solverUnsteadyWrapEnd


! ******************************************************************
! Inner utilities
! ******************************************************************
subroutine initTimeStepWrap
!
!      ******************************************************************
!      *                                                                *
!      * Second part of initTimeStepPart1_ALE and                       *
!      * initTimeStepPart2_ALE                                          *
!      *                                                                *
!      ******************************************************************
!
  use blockPointers
  use communication
  use inputMotion
  use inputUnsteady
  use iteration
  use monitor
  use section
  use flowVarRefState
  implicit none
  !
  !      Local variables.
  !
  integer(kind=intType) :: nn
  real(kind=realType), dimension(nSections) :: tNewSec, deltaTSec
  integer(kind=intType) :: lale ! for ALE

  testChanging: if(changing_Grid .or. gridMotionSpecified) then

     do nn=1,nSections
        tNewSec(nn)   = timeUnsteady + timeUnsteadyRestart
        deltaTSec(nn) = deltaT
     enddo


     ! --------------------------------
     ! These steps are done in Python:
         ! if( deforming_Grid )  then
         !    call shiftCoorAndVolumes
         !    call shiftLevelALE
         ! endif
         ! call callback(timeUnsteady)
     ! --------------------------------


     ! Added 100915
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

     ! Added 100915
     call updateSlidingAllLevels
     call updateMetricsAllLevels

     ! --------------------------------
     ! Next do a loop to interpolated mesh for intermediate configurations
     ! So that surface normals, normal velocities and slip velocities can
     ! be generated
     ! --------------------------------

     ! Determine the velocities of the cell centers and faces
     ! for the current ground level. Note that the spectral mode
     ! is always 1 for unsteady mode.

     ! --------------------------------
     ! First update cell and surface velocity, both are vectors
     ! Only quantities in blocks are updated, and they will not
     ! be interpolated
     ! --------------------------------
     call gridVelocitiesFineLevelPart1(deforming_Grid, tNewSec, 1_intType)

     ! --------------------------------
     ! Secondly store x to a temporary variable
     ! --------------------------------
     call storeCoor

     ! --------------------------------
     ! Thirdly update surface normal and normal velocity
     ! --------------------------------
     ALEloop : do lale = 1,nALEMeshes
        ! Interpolate mesh over latest time step
        call interpCoor(lale)

        ! Update s[I,J,K], norm
        ! call metric(currentLevel)
        call metric(groundLevel)

        ! Update sFace[I,J,K]
        call gridVelocitiesFineLevelPart2(deforming_Grid, tNewSec, 1_intType)
        ! Update uSlip
        call slipVelocitiesFineLevel_ALE(deforming_Grid, tNewSec, 1_intType) ! Maybe Unnecessary - to be verified

        ! Added 101115
        call gridVelocitiesCoarseLevels(1_intType)
        call slipVelocitiesCoarseLevels(1_intType)



        ! Update rFace
        call normalVelocitiesAllLevels(1_intType)

        ! Store data to *lale* ALE level
        call setLevelALE(lale)
     enddo ALEloop

     ! --------------------------------
     ! Lastly recover x from temporary variable
     ! Then compute data for current level
     ! --------------------------------
     call recoverCoor

     ! call metric(currentLevel)
     call metric(groundLevel)
     call gridVelocitiesFineLevelPart2(deforming_Grid, tNewSec, 1_intType)
     call slipVelocitiesFineLevel_ALE(deforming_Grid, tNewSec, 1_intType) ! Maybe Unnecessary - to be verified
     ! call slipVelocitiesFineLevel(deforming_Grid, tNewSec, 1_intType)

     ! Moved down
     ! call normalVelocitiesAllLevels(1_intType)

     ! Determine the velocities of the cell centers and faces and
     ! the slip velocities on the coarse grid levels . Note that the
     ! spectral mode is always 1 for unsteady mode.

     call gridVelocitiesCoarseLevels(1_intType)
     call slipVelocitiesCoarseLevels(1_intType)
     call normalVelocitiesAllLevels(1_intType)


     ! Compute the normal velocities of the boundaries, if
     ! needed for the corresponding boundary condition.

     ! --------------------------------
     ! Moved to Step1
     ! --------------------------------
     ! call normalVelocitiesAllLevels(1_intType)

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

end subroutine initTimeStepWrap

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
subroutine solverUnsteady_ALE(mdcallback_python)
  !
  !      ******************************************************************
  !      *                                                                *
  !      * solverUnsteadyMD  solves the unsteady equations using the BDF  *
  !      * schemes for the multigrid level groundLevel. However, it       *
  !      * this version includes a call-back to python for aero-elastic   *
  !      * simulations                                                    *
  !      *                                                                *
  !      * Modified for ALE scheme, by HDN                                *
  !      *                                                                *
  !      ******************************************************************
  !
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
  integer(kind=intType) :: iter, nTimeSteps
  integer(kind=intType) :: i,j,k,nn,kk
  external mdcallback_python

  !      ******************************************************************
  !      *                                                                *
  !      * Begin execution                                                *
  !      *                                                                *
  !      ******************************************************************
  
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

  ! Initialize timeStepUnsteady to 0 and set the number of
  ! time steps depending on the grid level.

  timeStepUnsteady = 0

  nTimeSteps = nTimeStepsCoarse
  if(groundLevel == 1) nTimeSteps = nTimeStepsFine


  ! Fill up old xold and volold
  call fillCoor

  ! Set all ALE levels by initial configuration
  call setLevelALE(-1_intType)

  
  ! Loop over the number of time steps to be computed.
 
  timeLoop: do iter=1,nTimeSteps

     ! Perform the two initialization tasks for this time step.
     ! These are split into two, such that in python mode something
     ! could be added, i.e. an additional grid velocity.

     call initTimeStepPart1_ALE(mdcallback_python)

     call initTimeStepPart2_ALE

     ! Solve the state for the current time step and 
     ! update nOldSolAvail.



!------------


     call solveState
     nOldSolAvail = nOldSolAvail + 1

     ! Determine whether or not solution files must be written.
     call checkWriteUnsteadyInLoop

     ! Exit the loop if the corresponding kill signal
     ! has been received.
     if(globalSignal == signalWriteQuit) exit     
  enddo timeLoop

  ! Determine whether or not the final solution must be written.
  
  call checkWriteUnsteadyEndLoop

end subroutine solverUnsteady_ALE

subroutine initTimeStepPart1_ALE(callback)
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
!      * Modified for ALE scheme, by HDN                                *
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
       external callback
       real(kind=realType), dimension(nSections) :: tNewSec, deltaTSec
       integer(kind=intType) :: bale, iale, jale, kale, lale ! for ALE
!
       integer(kind=intType) :: itest, jtest, ktest, ltest ! for ALE testing
       real(kind=realType)   :: coef
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

         if( deforming_Grid )  then
            call shiftCoorAndVolumes
            call shiftLevelALE
         endif

         !-------------------------------------------------------
         call callback(timeUnsteady)
         !-------------------------------------------------------
         ! Hardcode the deformation here
         ! do ktest = 0,3
         !    do jtest = 64,83
         !       do itest = 183,205
         !          x(itest, jtest, ktest, 1) = x(itest, jtest, ktest, 1) &
         !               + 0.004*(1.0-abs(x(itest, jtest, ktest, 1)-0.34)/0.14)*&
         !               (2.0*pi*10.0)*sin(2.0*pi*10.0*timeUnsteady)*deltaT
         !          x(itest, jtest, ktest, 2) = x(itest, jtest, ktest, 2) &
         !               + 0.005*(1.0-abs(x(itest, jtest, ktest, 2)-0.25)/0.15)*&
         !               (2.0*pi*10.0)*sin(2.0*pi*10.0*timeUnsteady)*deltaT
         !       enddo
         !    enddo
         ! enddo
         !
         ! !----blk16
         ! do ktest = 2,16
         !    do jtest = 2,16
         !       do itest = 2,16
         !          coef = 0.4*(1.0-abs(0.625*(itest-1.0)-5.0)/5.0)*(1.0-abs(0.625*(jtest-1.0)-5.0)/5.0)*&
         !               (1.0-abs(0.625*(ktest-1.0)-5.0)/5.0)
         !          x(itest, jtest, ktest, 1) = 0.625*(itest-1.0)+coef*sin(2*pi/0.5*timeUnsteady)
         !          x(itest, jtest, ktest, 2) = 0.625*(jtest-1.0)+coef*sin(2*pi/0.5*timeUnsteady)
         !          x(itest, jtest, ktest, 3) = 0.625*(ktest-1.0)+coef*sin(2*pi/0.5*timeUnsteady)
         !       enddo
         !    enddo
         ! enddo
         ! print *, 'At time ', timeUnsteady
         !
         !----pistf
         ! do itest = 2,100
         !    coef = 0.8*(1.0-abs((itest-1.0)-50.0)/50.0)
         !    do ktest = 1,4
         !       do jtest = 1,4
         !          x(itest, jtest, ktest, 1) = (itest-1.0)+coef*sin(2*pi/0.08*timeUnsteady)
         !       enddo
         !    enddo
         ! enddo
         ! print *, 'At time ', timeUnsteady
         !
         !----pistw
         ! do itest = 2,101
         !    coef = 50*(1.0-abs((itest-1.0)-100.0)/100.0)*timeUnsteady
         !    do ktest = 1,4
         !       do jtest = 1,4
         !          x(itest, jtest, ktest, 1) = (itest-1.0)-coef
         !       enddo
         !    enddo
         ! enddo
         ! print *, 'At time ', timeUnsteady
         !-------------------------------------------------------

         ! --------------------------------
         ! Assuming this subroutine is irrelevant to ALE
         ! --------------------------------
         ! call updateCoorFineMesh(deltaTSec, 1_intType)

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

! ===========================================================
!
! Assuming ALE has nothing to do with mesh updating
!
! Next do a loop to interpolated mesh for intermediate configurations
! So that surface normals, normal velocities and slip velocities can
! be generated
!
! ===========================================================

         ! --------------------------------
         ! Assuming this subroutine is irrelevant to ALE
         ! --------------------------------
         ! call faceRotationMatrices(currentLevel, .false.)

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
            call metric(currentLevel)
            ! Update sFace[I,J,K]
            call gridVelocitiesFineLevelPart2(deforming_Grid, tNewSec, 1_intType)
            ! Update uSlip
            call slipVelocitiesFineLevel_ALE(deforming_Grid, tNewSec, 1_intType)
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

         call metric(currentLevel)
         call gridVelocitiesFineLevelPart2(deforming_Grid, tNewSec, 1_intType)
         call slipVelocitiesFineLevel_ALE(deforming_Grid, tNewSec, 1_intType)
         call normalVelocitiesAllLevels(1_intType)

! ===========================================================
!
! Modification for ALE ends here
!
! ===========================================================

      endif testChanging
     
end subroutine initTimeStepPart1_ALE


subroutine initTimeStepPart2_ALE
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

  testChanging2: if(changing_Grid .or. gridMotionSpecified) then

     ! Determine the velocities of the cell centers and faces and
     ! the slip velocities on the coarse grid levels . Note that the
     ! spectral mode is always 1 for unsteady mode.

     call gridVelocitiesCoarseLevels(1_intType)
     call slipVelocitiesCoarseLevels(1_intType)

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

  endif testChanging2

  ! Shift the old solution for the new time step.

  call shiftSolution

  ! Set the coefficients for the time integrator.

  call setCoefTimeIntegrator_ALE

  ! Determine the data for the outflow bleeds.

  call BCDataMassBleedOutflow(.false., .false.)

end subroutine initTimeStepPart2_ALE

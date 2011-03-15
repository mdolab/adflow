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
subroutine solverUnsteadyMD(mdcallback_python)
  !
  !      ******************************************************************
  !      *                                                                *
  !      * solverUnsteadyMD  solves the unsteady equations using the BDF  *
  !      * schemes for the multigrid level groundLevel. However, it       *
  !      * this version includes a call-back to python for aero-elastic   *
  !      * simulations         
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
  use inputTimeSpectral
  implicit none
  !
  !      Local variables.
  
  logical :: EulerWallsPresent
  integer(kind=intType) :: iter, nTimeSteps
  real(kind=realType), dimension(nSections) :: tNewSec, deltaTSec
  real(kind=realType) :: trueTime
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
  
  timeUnsteady = zero

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
       
  spectralLoop: do kk=1,nTimeIntervalsSpectral
     domains: do nn=1,nDom
             
        ! Set the pointers for this block on the ground level.
        
        call setPointers(nn, groundLevel,kk)
        
        do k=0,ke
           do j=0,je
              do i=0,ie
                 xOld(:,i,j,k,1) = x(i,j,k,1)
                 xOld(:,i,j,k,2) = x(i,j,k,2)
                 xOld(:,i,j,k,3) = x(i,j,k,3)
              enddo
           enddo
        enddo
             
        do k=2,kl
           do j=2,jl
              do i=2,il
                 volOld(:,i,j,k) = vol(i,j,k)
              enddo
           enddo
        enddo
     end do domains
  end do spectralLoop
  
  ! Loop over the number of time steps to be computed.
 
  timeLoop: do iter=1,nTimeSteps

     ! Increment timeStepUnsteady and update
     ! timeUnsteady with the current time step.

     timeStepUnsteady = timeStepUnsteady + 1
     timeUnsteady     = timeUnsteady     + deltaT

     trueTime = (timeStepUnsteady)*deltaT

     ! Write the unsteady header. Only done by processor 0
     ! to avoid a messy output.

     if(myID == 0) call unsteadyHeader

     ! Set the new time for all sections; also store their
     ! time step. They are the same for all sections, but all
     ! of them should be updated because of consistency.

     if (deforming_grid) then
        call shiftCoorAndVolumes
     end if

     ! Shift the old solution for the new time step.

     call shiftSolution

     ! Set the coefficients for the time integrator.

     call setCoefTimeIntegrator

     ! Determine the data for the outflow bleeds.

     call BCDataMassBleedOutflow(.false., .false.)

     ! Solve the state for the current time step and 
     ! update nOldSolAvail.

     ! The "Multi-displinary" loop is actually here. We will call
     ! solve state with a (small) number of iterations, then solve the
     ! structures and solve the state again

     nCycles = 4

     do i=1,20
        
        ! Call the solve state func
        call solveState

        if (i == 1) then
           call MDUpdateAllUnsteady(mdcallback_python,trueTime,.True.)
        else
           call MDUpdateAllUnsteady(mdcallback_python,trueTime,.False.)
        end if

     end do

     nOldSolAvail = nOldSolAvail + 1

     ! Determine whether or not solution files must be written.

     call checkWriteUnsteadyInLoop

     ! Exit the loop if the corresponding kill signal
     ! has been received.

     if(globalSignal == signalWriteQuit) exit

  enddo timeLoop

  ! Determine whether or not the final solution must be written.

  call checkWriteUnsteadyEndLoop

end subroutine solverUnsteadyMD


subroutine MDUpdateAllUnsteady(mdcallback_python,t,firstIter)

  use inputIteration
  use inputUnsteady
  use blockPointers
  use inputDiscretization
  use iteration
  use killSignals
  use monitor
  use communication 
  use section
  use inputTimeSpectral

  implicit none
  real(kind=realType) :: t
  logical :: FirstIter
  integer(kind=intType) :: nn
  real(kind=realType), dimension(nSections) :: tNewSec, deltaTSec
  external mdcallback_python


  do nn=1,nSections
     tNewSec(nn)   = timeUnsteady + timeUnsteadyRestart
     deltaTSec(nn) = deltaT
  enddo

  call mdcallback_python(t,firstIter)

  call updateCoordinatesAllLevels
  if(updateWallDistanceUnsteady) then
     call updateWallDistanceAllLevels
  end if

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

end subroutine MDUpdateAllUnsteady

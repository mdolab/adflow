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
  use inputMotion
  use inputTimeSpectral
  implicit none
  !
  !      Local variables.
  
  logical :: EulerWallsPresent
  integer(kind=intType) :: iter, nTimeSteps
  real(kind=realType), dimension(nSections) :: tNewSec, deltaTSec
  real(kind=realType) :: temp
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

     ! Perform the two initialization tasks for this time step.
     ! These are split into two, such that in python mode something
     ! could be added, i.e. an additional grid velocity.

     call initTimeStepPart1_md(mdcallback_python)
     !call mdcallback_python(timeUnsteady)
     call initTimeStepPart2

     ! Solve the state for the current time step and 
     ! update nOldSolAvail.
     call solveState

     nOldSolAvail = nOldSolAvail + 1
     
  enddo timeLoop

  ! Determine whether or not the final solution must be written.
  
  call checkWriteUnsteadyEndLoop

end subroutine solverUnsteadyMD

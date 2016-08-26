!
!       File:          solveState.F90                                  
!       Author:        Edwin van der Weide                             
!       Starting date: 01-29-2004                                      
!       Last modified: 10-31-2007                                      
!
subroutine solveState
  !
  !       solveState computes either the steady or unsteady state        
  !       vector for the multigrid level groundLevel, which can be       
  !       found in the module iteration.                                 
  !
  use constants
  use communication, only : myID, sumb_comm_world
  use NKSolverVars, only : NKLSFuncEvals, freestreamResset, NK_LS, &
       rhoRes0, totalR, rhoResStart, totalR0, totalRFinal, totalRStart
  use inputio, only : forcedLiftFile, forcedSliceFile, forcedVolumeFile, &
       forcedSurfaceFile, solFile, newGridFile, surfaceSolFile
  use inputIteration, only: CFL, CFLCoarse, minIterNum, nCycles, &
       nCyclesCoarse, nMGSteps, nUpdateBleeds, printIterations
  use iteration, only : cycling, approxTotalIts, converged, CFLMonitor, &
       groundLevel, iterTot, iterType, currentLevel
  use killSignals, only : globalSignal, localSignal, noSignal, routineFailed, signalWrite, &
       signalWriteQuit
  use monitor, only : writeGrid, writeSurface, writeVolume
  use nksolvervars, only : NK_switchTol, useNKSolver, NK_CFL, rkREset
  use anksolvervars, only : ANK_switchTol, useANKSolver, ANK_CFL

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

  ! Write a message. Only done by processor 0.

  if(myID == 0) then

     ! Write a message about the number of multigrid iterations
     ! to be performed.
     
     write(numberString,"(i6)") nMGCycles		
     numberString = adjustl(numberString)
     numberString = trim(numberString)
     if (printIterations) then
        print "(a)", "#"
        print 102, groundLevel, trim(numberString)
        print "(a)", "#"
     end if
102  format("# Grid",1X,I1,": Performing",1X,A,1X, &
          "iterations, unless converged earlier")

     if (printIterations)  then
        call convergenceHeader
     end if
  end if

  ! Only compute the free stream resisudal once for efficiency on the
  ! fine grid only. 
  if (groundLevel == 1)  then
     if (.not. freeStreamResSet)  then
        call getFreeStreamResidual(rhoRes0, totalR0)
        freeStreamResSet = .True.
     end if
  end if

  ! Initialize the approxiate iteration count. We won't count this
  ! first residual evaluation. This way on the first convergence info
  ! call, "Iter" and "Iter Total" will both display 0. 
  approxTotalIts = 0
 
  ! Evaluate the initial residual
  call computeResidualNK

  ! Extract the rhoResStart and totalRStart
  call getCurrentResidual(rhoResStart, totalRStart)
     
  ! No iteration type for first residual evaluation
  iterType = "  None"

  ! Determine and write the initial convergence info.
  call convergenceInfo

  ! Loop over the maximum number of nonlinear iterations
  firstANK = .True.
  firstNK = .True. 

  nonlinearIteration: do while (approxTotalIts < nMGCycles)

     ! Update iterTot 
     iterTot = iterTot + 1

     if(mod(iterTot, nWriteConvHeader) == 0 .and. &
          myID == 0 .and. printIterations) then
        call convergenceHeader
     endif

     ! Determine what type of update to do:
      if (currentLevel > 1) then 

         ! Coarse grids do RK/DADI always
         call executeMGCycle
         CFLMonitor = CFLCoarse
     else
        if (.not. useANKSolver .and. .not. useNKSolver .or. (iterTot < minIterNum .and. rkreset)) then 

           ! Always RK/DADI or a RK startup. Run the MG Cycle
           
           call executeMGCycle
           CFLMonitor = CFL
           
        else if (useANKSolver .and. .not. useNKSolver) then 

           ! Approx solver, no NKSolver
           
           if (totalR > ANK_switchTol * totalR0) then 
                
              call executeMGCycle
              CFLMonitor = CFL

           else
              call ANKStep(firstANK)
              firstANK = .False.
              iterType = "   ANK"
              CFLMonitor = ANK_CFL
           end if

        else if (.not. useANKSolver .and. useNKSolver) then 

           ! NK Solver no approx solver

           if (totalR > NK_switchTol * totalR0) then 
                
              call executeMGCycle
              CFLMonitor = CFL
           else
              
              call NKStep(firstNK)
              firstNK = .False.
              iterType = "    NK"
              CFLMonitor = NK_CFL
           end if

        else if (useANKSolver .and. useNKSolver) then 

           ! Both approximate and NK solver. 

           if (totalR > ANK_switchTol*totalR0) then 

              call executeMGCycle
              CFLMonitor = CFL

           else if (totalR <= ANK_switchTol*totalR0 .and. &
                totalR > NK_switchTol*totalR0) then 

              iterType = "   ANK"
              call ANKStep(firstANK)
              firstANK = .False.
              firstNK = .True.
              CFLMonitor = ANK_CFL

           else

              iterType = "    NK"
              call NKStep(firstNK)
              firstNK = .False.
              firstANK = .True.
              CFLMonitor = NK_CFL

           end if
        end if
     end if

     ! Determine and write the convergence info.
     call convergenceInfo

     totalRFinal = totalR

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
          sumb_integer, mpi_max, SUmb_comm_world, &
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
        
        call writeSol()

        ! Also write potential tecplot files
        call writeTecplot(forcedSliceFile, .True., forcedLiftFile, .True., "", .False.)
        
        ! Reset the signal 
        localSignal = noSignal
     end if

     if (globalSignal == signalWriteQuit) then 
        exit nonLinearIteration
     end if


  enddo nonLinearIteration

  ! deallocate space for storing hisotry of function evaluations for NK
  ! solver with non-monotone line search if necessary
  if (currentLevel == 1 .and. useNKSolver .and. NK_LS==nonMonotoneLineSearch) then 
     deallocate(NKLSFuncEvals)
  end if

end subroutine solveState


!
!      ******************************************************************
!      *                                                                *
!      * File:          solveState.F90                                  *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 01-29-2004                                      *
!      * Last modified: 10-31-2007                                      *
!      *                                                                *
!      ******************************************************************
!
subroutine solveState
  !
  !      ******************************************************************
  !      *                                                                *
  !      * solveState computes either the steady or unsteady state        *
  !      * vector for the multigrid level groundLevel, which can be       *
  !      * found in the module iteration.                                 *
  !      *                                                                *
  !      ******************************************************************
  !
  use communication
  use constants
  use NKSolverVars
  use flowVarRefState
  use inputIteration
  use inputPhysics
  use iteration
  use killSignals
  use monitor
  use blockPointers
  implicit none
  !
  !      Local parameter
  !
  integer(kind=intType), parameter :: nWriteConvHeader = 50
  !
  !      Local variables.
  !
  integer :: ierr


  integer(kind=intType) :: iter, nMGCycles
  character (len=7) :: numberString
  logical :: solve_NK, solve_RK
  real(kind=realType) :: rhoRes1,totalRRes1
  real(kind=realType) :: l2convsave
  !
  !      ******************************************************************
  !      *                                                                *
  !      * Begin execution                                                *
  !      *                                                                *
  !      ******************************************************************
  !
  ! Allocate the memory for cycling.


  allocate(cycling(nMGSteps), stat=ierr)
  if(ierr /= 0)                  &
       call terminate("solveState", &
       "Memory allocation failure for cycling")

  ! Some initializations.
  
  writeVolume  = .false.
  writeSurface = .false.
  converged    = .false.
  globalSignal = noSignal
  localSignal  = noSignal
  iterTot      = 0

  ! Determine the initial residual.
  rkStage = 0
  currentLevel = groundLevel
  call timeStep(.false.)

  if(turbSegregated .or. turbCoupled) then
     call computeUtau
     call initres(nt1MG, nMGVar)
     call turbResidual
  endif

  call initres(1_intType, nwf)
  call residual

  ! Set the cycle strategy to a single grid iteration.

  nStepsCycling = 1
  cycling(1)    = 0

  ! If single grid iterations must be performed, write a message.

  if(nsgStartup > 0) then

     ! Write a message to stdout that some single grid
     ! iterations will be done. Only processor 0 does this.

     if(myID == 0) then
        write(numberString,"(i6)") nsgStartup
        numberString = adjustl(numberString)
        numberString = trim(numberString)

        if (printIterations) then
           print "(a)", "#"
           print 100, groundLevel, trim(numberString)
           print "(a)", "#"

100        format("# Grid",1X,I1,": Performing",1X,A,1X, "single grid &
                &startup iterations, unless converged earlier")
        end if
     endif

     ! Write the sliding mesh mass flow parameters (not for
     ! unsteady) and the convergence header.


     if(equationMode == steady .or. &
          equationMode == timeSpectral) call writeFamilyMassflow
     if(myID == 0) call convergenceHeader

     ! Determine and write the initial convergence info.

     call convergenceInfo

  endif

  ! Loop over the number of single grid start up iterations.


  singleGrid: do iter=1,nsgStartup

     ! Rewrite the sliding mesh mass flow parameters (not for
     ! unsteady) and the convergence header after a certain number
     ! of iterations. The latter is only done by processor 0.

     if(mod(iter,nWriteConvHeader) == 0) then
        if(equationMode == steady .or. &
             equationMode == timeSpectral) call writeFamilyMassflow
        if(myID == 0) call convergenceHeader
     endif

     ! Update iterTot and call executeMGCycle for executing
     ! a single grid cycle.

     iterTot = iterTot + 1
     call executeMGCycle

     ! Update PV3 solution for a steady computation if PV3 support
     ! is required. For unsteady mode the PV3 stuff is updated
     ! after every physical time step.


#ifdef USE_PV3
     !eran-viz         if(equationMode == steady .or.    &
     !eran-viz            equationMode == timeSpectral) &
     call pv_update(real(iterTot,realPV3Type))
#endif

     ! Determine and write the convergence info.

     call convergenceInfo

     ! The signal stuff must be done after every iteration only in
     ! steady spectral mode. In unsteady mode the signals are
     ! handled in the routine solverUnsteady.

     testSteady1: if(equationMode == steady .or. &
          equationMode == timeSpectral) then

        ! Determine the global kill parameter
        ! if signals are supported.

#ifndef USE_NO_SIGNALS

        call mpi_allreduce(localSignal, globalSignal, 1,         &
             sumb_integer, mpi_max, SUmb_comm_world, &
             ierr)
#endif

        ! Check whether a solution file, either volume or surface,
        ! must be written. Only on the finest grid level in stand
        ! alone mode.

        if(standAloneMode .and. groundLevel == 1) then
           writeVolume  = .false.
           writeSurface = .false.

           if(mod(iterTot, nSaveVolume)  == 0) writeVolume  = .true.
           if(mod(iterTot, nSaveSurface) == 0) writeSurface = .true.

           if(globalSignal == signalWrite .or. &
                globalSignal == signalWriteQuit) then

              writeVolume  = .true.
              writeSurface = .true.
           endif

           ! Make a distinction between steady and spectral
           ! mode. In the former case the grid will never
           ! be written; in the latter case when only when
           ! the volume is written.

           ! Check this: NOT true for optimization?

           writeGrid = .false.
           if(equationMode == timeSpectral .and. writeVolume) &
                writeGrid = .true.

           if(writeGrid .or. writeVolume .or. writeSurface) &
                call writeSol

        endif

        ! Reset the value of localSignal.

        localSignal = noSignal

     endif testSteady1

     ! Exit the loop if the corresponding kill signal
     ! has been received or if the solution is converged.

     if(globalSignal == signalWriteQuit .or. converged) exit

     ! Check if the bleed boundary conditions must be updated and
     ! do so if needed.

     if(mod(iter, nUpdateBleeds) == 0) &
          call BCDataMassBleedOutflow(.false., .false.)

  enddo singleGrid

  ! If the previous loop was ended due to a writeQuit signal
  ! or if the solution is converged go to the end of this routine.

  if(globalSignal == signalWriteQuit .or. converged) goto 99

  ! Determine the cycling strategy for this level.

  call setCycleStrategy

  ! Determine the number of multigrid cycles, which depends
  ! on the current multigrid level.

  nMGCycles = nCycles
  if(groundLevel > 1) nMGCycles = nCyclesCoarse

  ! Write a message. Only done by processor 0.

  if(myID == 0) then

     ! If single grid iterations were performed, write a message
     ! that from now on multigrid is turned on.

     if(nsgStartup > 0) then
        print "(a)", "#"
        print 101, groundLevel
        print "(a)", "#"
101     format("# Grid",1X,I1,": Switching to multigrid iterations")
     endif

     notCDB1 : if(.not.componentsBreakDown)then  !  eran-CBD
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
102     format("# Grid",1X,I1,": Performing",1X,A,1X, &
             "multigrid iterations, unless converged earlier")
     end if notCDB1  !  eran-CBD
  endif

  notCDB2 : if(.not.componentsBreakDown)then  !  eran-CBD

     ! Write the sliding mesh mass flow parameters (not for unsteady)
     ! and the convergence header.

     if(equationMode == steady .or. &
          equationMode == timeSpectral) call writeFamilyMassflow
     if(myID == 0) then 
        if (printIterations)  then
           call convergenceHeader
        end if
     end if
  end if	notCDB2 !  eran-CBD

  ! Determine and write the initial convergence info.
  ! Note that if single grid startup iterations were made, this
  ! value is printed twice.
  call convergenceInfo
  if(routineFailed)then
     ! Release the memory of cycling.
     deallocate(cycling, stat=ierr)
     if(ierr /= 0)                 &
          call terminate("solveState", &
          "Deallocation failure for cycling")
     ! Reset the L2Convergence in case it had been changed above
     L2Conv = L2ConvSave
     return
  endif

  !
  !         In a components-break-down run, do not execute MG, but exit
  !
  if(componentsBreakDown)go to 99 !  eran-CBD

  ! Determine if we need to run the RK solver, the NK solver or Both.

  if ( groundLevel .ne. 1_intType .or. .not. useNKSolver) then
     ! If we are not on the finest level, or if we don't want to use
     ! the NKsolver we will ALWAYS solve RK and NEVER NK.
     solve_RK = .True.
     solve_NK = .False.
     L2ConvSave = L2Conv

     if (groundLevel == 1 .and. fromPython) then
        if (.not. freeStreamResSet)  then
           call getFreeStreamResidual(rhoRes0,totalR0)
           freeStreamResSet = .True.
        end if
        call getCurrentResidual(rhoResStart,totalRStart)
     end if

  else ! We want to use the NKsolver AND we are on the fine grid. 

     ! Now we must determine if the solution is converged sufficently
     ! to start directly into the NKsolver OR if we have to run the
     ! RKsolver first to improve the starting point and then the NKsolver

     ! Get Frestream and starting residuals
     if (.not. freeStreamResSet) then
        call getFreeStreamResidual(rhoRes0, totalR0)
        freeStreamResSet = .True.
     end if
     call getCurrentResidual(rhoResStart,totalRStart)

     ! Determine if we need to run the RK solver, before we can
     ! run the NK solver 
     L2ConvSave = L2Conv
     if (totalRStart/totalR0 > NK_Switch_Tol) then

        ! If we haven't run anything yet, set rhoResL1Start...this is
        ! the intial rho residual on the fine grid. This is what is
        ! used in convergence info
        if (rhoResL1Start < zero) then 
          rhoResL1Start = rhoResStart
        end if
        
        L2Conv = NK_Switch_Tol*rhoRes0/rhoResL1Start

        solve_RK = .True.
        solve_NK = .True.
     else 
        ! Nominally we don't run any RK iterations....
        solve_NK = .True.
        solve_RK = .False.
        if (RKreset) then
           ! EXCEPT if we explictly want the solver to run a (small)
           ! fixed number of iterations to re-globlaize the solution
           solve_RK = .True.
           nMGCycles = nRKreset
        end if
     end if
  end if

  if (solve_RK) then

     ! Loop over the number of multigrid cycles.

     multiGrid: do iter=1,nMGCycles

        ! Rewrite the sliding mesh mass flow parameters (not for
        ! unsteady) and the convergence header after a certain number
        ! of iterations. The latter is only done by processor 0.

        if(mod(iter,nWriteConvHeader) == 0) then
           if(equationMode == steady .or. &
                equationMode == timeSpectral) call writeFamilyMassflow
           if(myID == 0) then
              if (printIterations) then
                 call convergenceHeader
              end if
           end if
        endif

        ! Update iterTot and call executeMGCycle.

        iterTot = iterTot + 1
        call executeMGCycle

        ! Update PV3 solution for a steady computation if PV3 support
        ! is required. For unsteady mode the PV3 stuff is updated
        ! after every physical time step.

#ifdef USE_PV3
        if(equationMode == steady .or.    &
             equationMode == timeSpectral) &
             call pv_update(real(iterTot,realPV3Type))
#endif

        ! Determine and write the convergence info.
        call convergenceInfo

        !check for divergence or nan here
        if(routineFailed)then
           ! Release the memory of cycling.
           deallocate(cycling, stat=ierr)
           if(ierr /= 0)                 &
                call terminate("solveState", &
                "Deallocation failure for cycling")
           ! Reset the L2Convergence in case it had been changed above
           L2Conv = L2ConvSave
           return
        endif

        ! The signal stuff must be done after every iteration only in
        ! steady or spectral mode. In unsteady mode the signals are
        ! handled in the routine solverUnsteady.

        testSteady2: if(equationMode == steady .or. &
             equationMode == timeSpectral) then

           ! Determine the global kill parameter.

           call mpi_allreduce(localSignal, globalSignal, 1,         &
                sumb_integer, mpi_max, SUmb_comm_world, &
                ierr)

           ! Check whether a solution file, either volume or surface,
           ! must be written. Only on the finest grid level in stand
           ! alone mode.

           if(standAloneMode .and. groundLevel == 1) then

              writeVolume  = .false.
              writeSurface = .false.

              if(mod(iterTot, nSaveVolume)  == 0) writeVolume  = .true.
              if(mod(iterTot, nSaveSurface) == 0) writeSurface = .true.

              if(globalSignal == signalWrite .or. &
                   globalSignal == signalWriteQuit) then
                 writeVolume  = .true.
                 writeSurface = .true.
              endif

              ! Make a distinction between steady and spectral
              ! mode. In the former case the grid will never
              ! be written; in the latter case when only when
              ! the volume is written.

              writeGrid = .false.
              if(equationMode == timeSpectral .and. writeVolume) &
                   writeGrid = .true.

              if(writeGrid .or. writeVolume .or. writeSurface) &
                   call writeSol

           endif

           ! Reset the value of localSignal.

           localSignal = noSignal

        endif testSteady2

        ! Exit the loop if the corresponding kill signal
        ! has been received or if the solution is converged.

        if(globalSignal == signalWriteQuit .or. converged) exit

        ! Check if the bleed boundary conditions must be updated and
        ! do so if needed.

        if(mod(iter, nUpdateBleeds) == 0) &
             call BCDataMassBleedOutflow(.false., .false.)
        
     enddo multiGrid
  end if

  ! Reset the L2Convergence in case it had been changed above
  L2Conv = L2ConvSave

  ! Now we have run the RK solver. We will run the NK solver only if
  ! the solve_NK flag is set:
  if (solve_NK) then
     
     ! We have to check to see if NKSwitchtol was LOWER than
     ! l2convrel. This means that the RK solution is already good
     ! enough for what we want so we're done
    
     call getcurrentResidual(rhoRes1, totalRRes1)
     if (rhoRes1 < rhoResStart * l2convrel) then
        ! We no not have to run the NK solver so do nothing
     else
        ! Run the NK solver down as far we need to go
        call setupNKsolver() 
        call NKsolver()
     end if
  end if

  ! Finally...if we're on on the fine grid get the final residual
  if (groundLevel == 1) then
     call getCurrentResidual(rhoResFinal,totalRFinal)
  end if

  ! Release the memory of cycling.

99 deallocate(cycling, stat=ierr)
  if(ierr /= 0)                 &
       call terminate("solveState", &
       "Deallocation failure for cycling")

  ! Write the final values of the sliding mesh mass flow.

  call writeFamilyMassflow

  ! Write a solution. Only on the finest grid and if some iterations
  ! have been done since the last write. The solution is only
  ! written in stand alone mode for a steady or time spectral
  ! computation.

  testSteady3: if(equationMode == steady .or. &
       equationMode == timeSpectral) then

     notCBD : if(.not.componentsBreakDown)then   !  eran-CDB
       if(standAloneMode .and. groundLevel == 1) then

          writeVolume  = .not. writeVolume
          writeSurface = .not. writeSurface

        ! Make a distinction between steady and spectral
        ! mode. In the former case the grid will never
        ! be written; in the latter case when only when
        ! the volume is written.

          writeGrid = .false.
          if(equationMode == timeSpectral .and. writeVolume) &
             writeGrid = .true.

          if(writeGrid .or. writeVolume .or. writeSurface) &
               call writeSol


	   if(genCBDOUT)call componentsBreakDownPrintout(1) ! eran-CBD
       endif  ! standAloneMode .and. groundLevel 
     end if notCBD    !   eran-CDB
  endif testSteady3

end subroutine solveState

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
       use flowVarRefState
       use inputIteration
       use inputPhysics
       use iteration
       use killSignals
       use monitor
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

           print "(a)", "#"
           print 100, groundLevel, trim(numberString)
           print "(a)", "#"
 100       format("# Grid",1X,I1,": Performing",1X,A,1X, "single grid &
                  &startup iterations, unless converged earlier")
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
         if(equationMode == steady .or.    &
            equationMode == timeSpectral) &
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
 101       format("# Grid",1X,I1,": Switching to multigrid iterations")
         endif

         ! Write a message about the number of multigrid iterations
         ! to be performed.

         write(numberString,"(i6)") nMGCycles
         numberString = adjustl(numberString)
         numberString = trim(numberString)

         print "(a)", "#"
         print 102, groundLevel, trim(numberString)
         print "(a)", "#"
 102     format("# Grid",1X,I1,": Performing",1X,A,1X, &
                "multigrid iterations, unless converged earlier")
       endif

       ! Write the sliding mesh mass flow parameters (not for unsteady)
       ! and the convergence header.

       if(equationMode == steady .or. &
          equationMode == timeSpectral) call writeFamilyMassflow
       if(myID == 0) call convergenceHeader

       ! Determine and write the initial convergence info.
       ! Note that if single grid startup iterations were made, this
       ! value is printed twice.

       call convergenceInfo

       ! Loop over the number of multigrid cycles.

       multiGrid: do iter=1,nMGCycles

         ! Rewrite the sliding mesh mass flow parameters (not for
         ! unsteady) and the convergence header after a certain number
         ! of iterations. The latter is only done by processor 0.

         if(mod(iter,nWriteConvHeader) == 0) then
           if(equationMode == steady .or. &
              equationMode == timeSpectral) call writeFamilyMassflow
           if(myID == 0) call convergenceHeader
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

       ! Release the memory of cycling.

 99    deallocate(cycling, stat=ierr)
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

         endif
       endif testSteady3

       end subroutine solveState

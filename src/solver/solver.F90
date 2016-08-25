       subroutine solver
!
!      ******************************************************************
!      *                                                                *
!      * solver is the main subroutine of the solver library.           *
!      * It controls the full multigrid and the kill signals.           *
!      *                                                                *
!      ******************************************************************
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
             call solverSteady
          case (unsteady)
             select case (timeIntegrationScheme)
               case (BDF)
                 call solverUnsteadyBDF

               case (explicitRK)
                 call solverUnsteadyExplicitRK

              case (MD)
                 call solverUnsteady_ALE
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
100              format("# Going down to grid level",1X,I1)
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

!
!      ******************************************************************
!      *                                                                *
!      * File:          solver.F90                                      *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 03-12-2003                                      *
!      * Last modified: 11-21-2007                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine solver
!
!      ******************************************************************
!      *                                                                *
!      * solver is the main subroutine of the solver library.           *
!      * It controls the full multigrid and the kill signals.           *
!      *                                                                *
!      ******************************************************************
!
       use communication
       use constants
       use inputDiscretization
       use inputIteration
       use inputPhysics
       use inputTimeSpectral
       use inputUnsteady
       use killSignals
       use iteration
       use monitor
       use section
       implicit none
!
!      Local variables.
!
       real(kind=realType), dimension(nSections) :: dtAdvance
!
!      Function definitions.
!
       logical :: EulerWallsPresent
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! If the normal momentum equation should be used to determine
       ! the pressure in the halo for inviscid walls, find out if there
       ! actually are inviscid walls. If so, set the logical
       ! exchangePressureEarly to .true.; otherwise set it to .false.

       if(wallBcTreatment == normalMomentum) then
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

       ! Initialize PV3 routines if PV3 support is required and has not
       ! been called before. In order to make sure that the overset 
       ! iblanking works properly with pV3, set groundLevel to 1 so that
       ! the iblank arrays are allocated with the maximum size.

#ifdef USE_PV3
       if (.not. PV3Initialized) then
         groundLevel = 1
         call initializePV3
         PV3Initialized = .true.
       end if
#endif

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

               case (implicitRK)
                 call solverUnsteadyImplicitRK
             end select
         end select

         ! If this is not the finest grid level, interpolate the
         ! solution to the next finer mesh and write a message that
         ! this is happening. Only processor 0 performs this latter
         ! task.

         if(groundLevel > 1) then

           currentLevel = groundLevel - 1

           if(myID == 0) then
             print "(a)", "#"
             print 100, currentLevel
             print "(a)", "#"
 100         format("# Going down to grid level",1X,I1)
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

       end subroutine solver

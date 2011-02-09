!
!      ******************************************************************
!      *                                                                *
!      * File:          checkOutput.f90                                 *
!      * Author:        Edwin van der Weide                             *
!      * Starting date: 03-26-2003                                      *
!      * Last modified: 11-27-2007                                      *
!      *                                                                *
!      ******************************************************************
!
       subroutine checkOutput
!
!      ******************************************************************
!      *                                                                *
!      * checkOutput checks and possibly corrects the and output        *
!      * variables. This depends on the set of governing equations to   *
!      * be solved.                                                     *
!      *                                                                *
!      ******************************************************************
!
       use inputPhysics
       use inputUnsteady
       use flowVarRefState
       use extraOutput
       implicit none
!
!      ******************************************************************
!      *                                                                *
!      * Begin execution                                                *
!      *                                                                *
!      ******************************************************************
!
       ! Determine the governing equations to be solved and set the
       ! variables accordingly.

       select case (equations)
         case (EulerEquations)
           surfWriteCf    = .false.
           surfWriteCh    = .false.
           surfWriteYplus = .false.
           surfWriteCfx   = .false.
           surfWriteCfy   = .false.
           surfWriteCfz   = .false.

           volWriteMachTurb     = .false.
           volWriteEddyVis      = .false.
           volWriteRatioEddyVis = .false.
           volWriteDist         = .false.
           volWriteResTurb      = .false.

         case (NSEquations)
           surfWriteYplus = .false.

           volWriteMachTurb     = .false.
           volWriteEddyVis      = .false.
           volWriteRatioEddyVis = .false.
           volWriteDist         = .false.
           volWriteResTurb      = .false.

         case (RANSEquations)

           ! Check if it is possible to write a turbulent Mach
           ! number and eddy viscosity; this depends on the turbulence
           ! model used.

           if(.not. kPresent) volWriteMachTurb = .false.
           if(.not. eddyModel) then
             volWriteEddyVis      = .false.
             volWriteRatioEddyVis = .false.
           endif

           ! If a wall distance free turbulence model is used,
           ! set volWriteDist to .false.

           if(.not. wallDistanceNeeded) volWriteDist = .false.

       end select

       if(equationMode          == unsteady .and. &
          timeIntegrationScheme == explicitRK) then
         volWriteResRho  = .false.
         volWriteResMom  = .false.
         volWriteResRhoE = .false.
         volWriteResTurb = .false.
       endif

       end subroutine checkOutput
